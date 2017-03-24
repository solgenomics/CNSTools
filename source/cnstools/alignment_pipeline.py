from _utils import Progress_tracker, call_commands_async, safe_print, header_print
import file_handlers as fhs

import os
import shlex
import json

config_defaults = {
    #"reference":"A",
    #"query_genomes":{"A":/home/dal333/haudrey_test_data/crubella_183_v1.fa","B":"/home/dal333/haudrey_test_data/Alyrata_107_v1.fa"},
    #optional:
    "out_folder":"./",
    "chaining_script_bin_path":"",
    "num_processes":1,
    "lastz_options":"C=0 E=30 H=2000 K=2200 L=6000 M=50 O=400 T=2 Y=3400",#Q=/home/dal333/alignment_step/test_input/HoxD55.q",
    "chainNet_options":"-minSpace=25",
    "lastz_bin_path":"",
}

def parser(parser_add_func,name):
    p = parser_add_func(name,description="Aligns a genome to a reference")
    p.add_argument("config_path", help="Absolute(!) path to config file")
    return p

def run(config_path):

    with open(config_path) as config_file:
      config = json.loads(config_file.read())
    for key in config_defaults:
        config.setdefault(key, config_defaults[key])

    original_wd = os.getcwd()

    #grab vars from config and make sure the paths are absolute and interperted as such relative to the config file.
    os.chdir(os.path.dirname(os.path.abspath(config_path)))

    num_processes = config["num_processes"]
    lastz_bin_path = config["lastz_bin_path"]
    lastz_options = shlex.split(config["lastz_options"])
    chaining_script_bin_path = config["chaining_script_bin_path"]
    chainNet_options = shlex.split(config["chainNet_options"])
    out_folder = os.path.abspath(config["out_folder"])

    reference = config["reference"]
    ref_genome = os.path.abspath(config["genomes"][reference])
    query_genomes = {key:os.path.abspath(config["genomes"][key]) for key in config["genomes"] if key!=reference}

    #set the working directory to the out folder
    try:
        os.makedirs(out_folder)
    except OSError:
        if not os.path.isdir(out_folder):
            raise
    os.chdir(out_folder)
    #copy env so we can modify PATH for subprocesses
    cmd_env = os.environ.copy()

    if chaining_script_bin_path!="":
        cmd_env["PATH"] = chaining_script_bin_path+":" + cmd_env["PATH"]
    if lastz_bin_path!="":
        cmd_env["PATH"] = lastz_bin_path+":" + cmd_env["PATH"]

    header_print("Preparing reference (%s)..." % reference,h_type=1)

    #split fasta files to on file per sequence
    orig_fasta = fhs.fasta.Handler(ref_genome)
    split_fastas = [file.path for file in orig_fasta.split(1,file_prefix="%s."%reference,file_suffix="",out_folder=out_folder,tracker=Progress_tracker("Splitting %s"%os.path.basename(ref_genome),orig_fasta.size).auto_display(1))]
    
    #distribute the hard jobs (mostly so that the status bar moves more uniformly). The list is still weighted towards the left though. :)
    split_fastas.sort(key = lambda name: os.stat(name).st_size, reverse = True)
    split_fastas = [v for m in [split_fastas[i::num_processes] for i in range(num_processes)] for v in m] 

    #faSize, gets sizes of ref genomes, makes .sizes files
    #a similar sections to this section will repeat throughout the file, as much of the script is running other scripts!
    sizes_files = [os.path.splitext(file)[0]+".sizes" for file in split_fastas] #generate new names for the output files
    faSize_commandlists = [["faSize", in_file, "-detailed", ">", out_file] for in_file,out_file in zip(split_fastas,sizes_files)] #create a command/arguemnt list for each command to be run in a batch
    call_commands_async(faSize_commandlists,num_processes,shell=True,tracker_name="Sizing (faSize)",env=cmd_env) #runs commands asynchronously with maximum simultanious process count, waits for them to finish [shell=True indicates the command list should be joined with spaces and interpreted by a shell, the shell can be changed in the config file]

    #faToTwoBit makes .2bit files
    twobit_files = [os.path.splitext(file)[0]+".2bit" for file in split_fastas]
    faToTwoBit_commandlists = [["faToTwoBit", in_file, out_file] for in_file,out_file in zip(split_fastas,twobit_files)]
    call_commands_async(faToTwoBit_commandlists,num_processes,tracker_name="Converting to 2bit (faToTwoBit)",env=cmd_env)
    
    config["aligned_query_genomes"] = {}

    for query_genome_name in query_genomes:
        header_print("Aligning %s..."%query_genome_name,h_type=1)
        query_genome = query_genomes[query_genome_name]

        #faSize, gets size of genome fasta file and makes .sizes file
        query_genome_sizes = os.path.splitext(query_genome)[0]+".sizes"
        faSize_commandlists = [["faSize", query_genome, "-detailed", ">" , query_genome_sizes]]
        call_commands_async(faSize_commandlists,1,shell=True,tracker_name="Sizing (faSize)",env=cmd_env)

        #faToTwoBit makes genome .2bit file
        query_genome_twobit = os.path.splitext(query_genome)[0]+".2bit"
        faToTwoBit_commandlists = [["faToTwoBit", query_genome, query_genome_twobit]]
        call_commands_async(faToTwoBit_commandlists,1,tracker_name="Converting to 2bit (faToTwoBit)",env=cmd_env)

        #lastz, aligns and makes .lav files
        prefix = query_genome_name+".to."
        lav_files = [os.path.join(os.path.dirname(file),prefix+os.path.splitext(os.path.basename(file))[0]+".lav") for file in split_fastas]
        lastz_commandlists = [["lastz", in_file, query_genome, "--format=lav"] + lastz_options + [">", out_file] for in_file,out_file in zip(split_fastas,lav_files)]
        call_commands_async(lastz_commandlists,num_processes,shell=True,tracker_name="Aligning (lastz)",env=cmd_env)
        
        #lavToAxt, converts .lav files and adds sequences, makes .axt files
        axt_files = [os.path.splitext(file)[0]+".axt" for file in lav_files]
        lavToAxt_commandlists = [["lavToAxt","-fa","-tfa", in_file, ref_fasta, query_genome, out_file] for in_file,ref_fasta,out_file in zip(lav_files,split_fastas,axt_files)]
        call_commands_async(lavToAxt_commandlists,num_processes,shell=True,tracker_name="lavToAxt",env=cmd_env)

        #axtChain, Two matching alignments next to each other are joined into one fragment if they are close enough, and makes .chain files
        chain_files = [os.path.splitext(file)[0]+".chain" for file in axt_files]
        axtChain_commandlists = [["axtChain","-linearGap=loose","-faQ","-faT",in_file, ref_fasta, query_genome, out_file] for in_file,ref_fasta,out_file in zip(axt_files,split_fastas,chain_files)]
        call_commands_async(axtChain_commandlists,num_processes,tracker_name="axtChain",stderr=False,env=cmd_env)

        #chainPreNet makes .prenet files
        prenet_files = [os.path.splitext(file)[0]+".prenet" for file in chain_files]
        chainPreNet_commandlists = [["chainPreNet", in_file, size_file, query_genome_sizes, out_file] for in_file,size_file,out_file in zip(chain_files,sizes_files,prenet_files)]
        call_commands_async(chainPreNet_commandlists,num_processes,tracker_name="chainPreNet",env=cmd_env)

        #chainNet makes .no_synt.chainnet files
        no_synt_chainnet_files = [os.path.splitext(file)[0]+".chainnet" for file in prenet_files]
        chainNet_commandlists = [["chainNet"] + chainNet_options + [in_file, size_file, query_genome_sizes, out_file, "/dev/null"] for in_file,size_file,out_file in zip(prenet_files,sizes_files,no_synt_chainnet_files)]
        call_commands_async(chainNet_commandlists,num_processes,tracker_name="chainNet",stderr=False,env=cmd_env)

        #netSynteny makes .chainnet files
        chainnet_files = [os.path.splitext(file)[0]+".syntenic.chainnet" for file in prenet_files]
        netSynteny_commandlists = [["netSyntenic", in_file, out_file] for in_file,out_file in zip(no_synt_chainnet_files,chainnet_files)]
        call_commands_async(netSynteny_commandlists,num_processes,tracker_name="netSynteny",env=cmd_env)

        #netFilter makes .filter files
        filter_files = [os.path.splitext(file)[0]+".filter" for file in chainnet_files]
        netFilter_commands = [["netFilter ", in_file, ">", out_file] for in_file,out_file in zip(chainnet_files,filter_files)]
        call_commands_async(netFilter_commands,num_processes,shell=True,tracker_name="netFilter",env=cmd_env)

        ##
        ## Aligned and chained! Now convert to .maf file!
        ##

        #netToAxt makes .filter.axt files
        filtered_axt_files = [os.path.splitext(file)[0]+".filter.axt" for file in axt_files]
        netFilter_commands = [["netToAxt", in_file, prenet_file, twobit_file, query_genome_twobit, out_file] for in_file,prenet_file,twobit_file,out_file in zip(filter_files,prenet_files,twobit_files,filtered_axt_files)]
        call_commands_async(netFilter_commands,num_processes,tracker_name="netToAxt",env=cmd_env)

        #axtSort makes .sort.axt files
        sorted_axt_files = [os.path.splitext(file)[0]+".sort.axt" for file in filtered_axt_files]
        axtSort_commands = [["axtSort", in_file, out_file] for in_file,out_file in zip(filtered_axt_files,sorted_axt_files)]
        call_commands_async(axtSort_commands,num_processes,tracker_name="axtSort",env=cmd_env)

        #axtToMaf makes .maf files !!!!!!!
        maf_files = [os.path.splitext(file)[0]+".maf" for file in sorted_axt_files]
        axtToMaf_commands = [["axtToMaf", in_file, ref_size, query_genome_sizes, out_file] for in_file,ref_size,out_file in zip(sorted_axt_files,sizes_files,maf_files)]
        call_commands_async(axtToMaf_commands,num_processes,tracker_name="axtToMaf",env=cmd_env)
        
        config["aligned_query_genomes"][query_genome_name] = maf_files
        safe_print("%s Aligned."%query_genome_name)

    header_print("Alignments Complete.",h_type=0)
    with open(os.path.join("results.config.json"),"w") as out:
        json.dump(config,out,sort_keys=True,indent=4)

    os.chdir(original_wd)
