import subprocess
import os
import shlex
import json
from file_handler import Fasta_handler
from _utils import Progress_tracker

config_defaults = {
    #"ref_genome":"/home/dal333/haudrey_test_data/TAIR10_chr_all.fas",
    #"query_genomes":["/home/dal333/haudrey_test_data/crubella_183_v1.fa","/home/dal333/haudrey_test_data/Alyrata_107_v1.fa"],
    #"out_folder":"./",
    #optional:
    "chaining_script_directory":"",
    "num_processes":1,
    "lastz_options":"C=0 E=30 H=2000 K=2200 L=6000 M=50 O=400 T=2 Y=3400",#Q=/home/dal333/alignment_step/test_input/HoxD55.q",
    "chainNet_options":"-minSpace=25",
    "lastz_path":"lastz",
    "alt_shell_path":None,
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

    lastz_path = config["lastz_path"]
    chaining_script_directory = config["chaining_script_directory"]
    lastz_options = shlex.split(config["lastz_options"])
    chainNet_options = shlex.split(config["chainNet_options"])
    query_genomes = config["query_genomes"]
    ref_genome = config["ref_genome"]
    out_folder = config["out_folder"]
    if config["alt_shell_path"]!=None:
        call_commands_async.alt_shell_path = config["alt_shell_path"]

    num_processes = config["num_processes"] #make param

    split_fasta = [file.path for file in Fasta_handler(ref_genome).split(1,file_prefix="ref.",file_suffix=".split",out_folder=out_folder)]

    print "Reference file has %s sequences" % len(split_fasta)
    split_fasta.sort(key=lambda name:os.stat(name).st_size)
    split_fasta = [v for m in [split_fasta[i::num_processes] for i in range(num_processes)] for v in m] #distribute the hard jobs (mostly so that the status bar moves more uniformly). The list is still weighted towards the right. :)
    split_fasta_sizes = [os.stat(name).st_size for name in split_fasta]

    #faSize, gets sizes of ref genomes, makes .sizes files
    #a similar sections to this section will repeat throughout the file, as much of the script is running other scripts!
    sizes_files = [os.path.splitext(file)[0]+".sizes" for file in split_fasta] #generate new names for the output files
    faSize_commandlists = [[chaining_script_directory+"faSize", in_file, "-detailed", ">", out_file] for in_file,out_file in zip(split_fasta,sizes_files)] #create a command/arguemnt list for each command to be run in a batch
    call_commands_async(faSize_commandlists,num_processes,shell=True,tracker_name="faSize") #runs commands asynchronously with maximum simultanious process count, waits for them to finish [shell=True indicates the command list should be joined with spaces and interpreted by a shell, the shell can be changed in the config file]

    #faToTwoBit makes .2bit files
    twobit_files = [os.path.splitext(file)[0]+".2bit" for file in split_fasta]
    faToTwoBit_commandlists = [[chaining_script_directory+"faToTwoBit", in_file, out_file] for in_file,out_file in zip(split_fasta,twobit_files)]
    call_commands_async(faToTwoBit_commandlists,num_processes,tracker_name="faToTwoBit")

    for query_genome in query_genomes:

        #faSize, gets size of genome fasta file and makes .sizes file
        query_genome_sizes = os.path.splitext(query_genome)[0]+".sizes"
        faSize_commandlists = [[chaining_script_directory+"faSize", query_genome, "-detailed", ">" , query_genome_sizes]]
        call_commands_async(faSize_commandlists,1,shell=True,tracker_name="faSize")

        #faToTwoBit makes genome .2bit file
        query_genome_twobit = os.path.splitext(query_genome)[0]+".2bit"
        faToTwoBit_commandlists = [[chaining_script_directory+"faToTwoBit", query_genome, query_genome_twobit]]
        call_commands_async(faToTwoBit_commandlists,1,tracker_name="faToTwoBit")

        #lastz, aligns and makes .lav files
        prefix = os.path.splitext(os.path.basename(query_genome))[0]+"_to_"
        lav_files = [os.path.join(os.path.dirname(file),prefix+os.path.splitext(os.path.basename(file))[0]+".lav") for file in split_fasta]
        lastz_commandlists = [[lastz_path, in_file, query_genome, "--format=lav"] + lastz_options + [">", out_file] for in_file,out_file in zip(split_fasta,lav_files)]
        k = call_commands_async(lastz_commandlists,num_processes,shell=True,tracker_name="lastz")

        #print len(k),"/",len(split_fasta)
        
        #lavToAxt, converts .lav files and adds sequences, makes .axt files
        axt_files = [os.path.splitext(file)[0]+".axt" for file in lav_files]
        lavToAxt_commandlists = [[chaining_script_directory+"lavToAxt","-fa","-tfa", in_file, ref_fasta, query_genome, out_file] for in_file,ref_fasta,out_file in zip(lav_files,split_fasta,axt_files)]
        call_commands_async(lavToAxt_commandlists,num_processes,shell=True,tracker_name="lavToAxt")

        #axtChain, Two matching alignments next to each other are joined into one fragment if they are close enough, and makes .chain files
        chain_files = [os.path.splitext(file)[0]+".chain" for file in axt_files]
        axtChain_commandlists = [[chaining_script_directory+"axtChain","-linearGap=loose","-faQ","-faT",in_file, ref_fasta, query_genome, out_file] for in_file,ref_fasta,out_file in zip(axt_files,split_fasta,chain_files)]
        call_commands_async(axtChain_commandlists,num_processes,tracker_name="axtChain")

        #chainPreNet makes .prenet files
        prenet_files = [os.path.splitext(file)[0]+".prenet" for file in chain_files]
        chainPreNet_commandlists = [[chaining_script_directory+"chainPreNet", in_file, size_file, query_genome_sizes, out_file] for in_file,size_file,out_file in zip(chain_files,sizes_files,prenet_files)]
        call_commands_async(chainPreNet_commandlists,num_processes,tracker_name="chainPreNet")

        #chainNet makes .no_synt.chainnet files
        no_synt_chainnet_files = [os.path.splitext(file)[0]+".no_synt.chainnet" for file in prenet_files]
        chainNet_commandlists = [[chaining_script_directory+"chainNet"] + chainNet_options + [in_file, size_file, query_genome_sizes, out_file, "/dev/null"] for in_file,size_file,out_file in zip(prenet_files,sizes_files,no_synt_chainnet_files)]
        call_commands_async(chainNet_commandlists,num_processes,tracker_name="chainNet")

        #netSynteny makes .chainnet files
        chainnet_files = [os.path.splitext(file)[0]+".chainnet" for file in prenet_files]
        netSynteny_commandlists = [[chaining_script_directory+"netSyntenic", in_file, out_file] for in_file,out_file in zip(no_synt_chainnet_files,chainnet_files)]
        call_commands_async(netSynteny_commandlists,num_processes,tracker_name="netSynteny")

        #netFilter makes .filter files
        filter_files = [os.path.splitext(file)[0]+".filter" for file in chainnet_files]
        netFilter_commands = [[chaining_script_directory+"netFilter ", in_file, ">", out_file] for in_file,out_file in zip(chainnet_files,filter_files)]
        call_commands_async(netFilter_commands,num_processes,shell=True,tracker_name="netFilter")

        ##
        ## Aligned and chained! Now convert to .maf file!
        ##

        #netToAxt makes .filter.axt files
        filtered_axt_files = [os.path.splitext(file)[0]+".filter.axt" for file in axt_files]
        netFilter_commands = [[chaining_script_directory+"netToAxt", in_file, prenet_file, twobit_file, query_genome_twobit, out_file] for in_file,prenet_file,twobit_file,out_file in zip(filter_files,prenet_files,twobit_files,filtered_axt_files)]
        call_commands_async(netFilter_commands,num_processes,tracker_name="netToAxt")

        #axtSort makes .sort.axt files
        sorted_axt_files = [os.path.splitext(file)[0]+".sort.axt" for file in filtered_axt_files]
        axtSort_commands = [[chaining_script_directory+"axtSort", in_file, out_file] for in_file,out_file in zip(filtered_axt_files,sorted_axt_files)]
        call_commands_async(axtSort_commands,num_processes,tracker_name="axtSort")

        #axtToMaf makes .maf files !!!!!!!
        maf_files = [os.path.splitext(file)[0]+".maf" for file in sorted_axt_files]
        axtToMaf_commands = [[chaining_script_directory+"axtToMaf", in_file, ref_size, query_genome_sizes, out_file] for in_file,ref_size,out_file in zip(sorted_axt_files,sizes_files,maf_files)]
        call_commands_async(axtToMaf_commands,num_processes,tracker_name="axtToMaf")


def call_commands_async(command_iterable,num,shell=False,tracker_name="Running command"):
    process_list = []
    finished = []
    try:
        tracker = Progress_tracker(tracker_name,len(command_iterable)).estimate(False)
    except:
        tracker = None
    for command in command_iterable:
        if shell==True:
            if call_commands_async.alt_shell_path!=None:
                process_list.append(subprocess.Popen(" ".join(command),shell=True,executable=call_commands_async.alt_shell_path))
            else:
                process_list.append(subprocess.Popen(" ".join(command),shell=True))
        elif shell==False:
            process_list.append(subprocess.Popen(command))
        if tracker and len(process_list) >= num: tracker.status("%s/%s processes active"%(len(process_list),num))
        while len(process_list) >= num:
            pid,exitstat = os.waitpid(-1,0)
            for i in range(len(process_list)-1,-1,-1):
                if process_list[i].pid == pid or process_list[i].poll() != None:
                    finished.append(process_list.pop(i))
                    if tracker: tracker.step()
    for i in range(len(process_list)):
        proc = process_list.pop(0)
        proc.wait()
        finished.append(proc)
        if tracker: tracker.step().display().status("%s/%s processes active"%(len(process_list),num))
    if tracker: tracker.status().done()
    return finished
call_commands_async.alt_shell_path=None
