import os, json, random, copy
from ._utils import call_commands_async, call_command, MultiTracker
from .file_handlers import ss

def score(reference,
          aligned_query_genomes,
          tree,
          annotations,
          msa_view_path = "",
          multiz_path   = "",
          num_processes = 1,
          out_folder    = "./",
          phast_path    = "", 
          cds_seqtypes = ["CDS"],
          **kwargs):
    '''ROASTs and scores files for conservation'''

    tracker = MultiTracker("Scoring alignment",1,estimate=False,style="noProg")
    tracker.auto_display(1)

    ref_genome_gff =   annotations[reference]

    try:
        os.makedirs(out_folder)
    except OSError:
        if not os.path.isdir(out_folder):
            raise
    cmd_env = os.environ.copy()

    if multiz_path!="":
            cmd_env["PATH"] = multiz_path+":" + cmd_env["PATH"]
    if phast_path!="":
            cmd_env["PATH"] = phast_path+":" + cmd_env["PATH"]

    per_chrom_labeled_mafs = {}
    per_chrom_existing_genomes = {}
    out_name_template = os.path.join(out_folder,"{chrom}.{query_genome}.sing.maf")
    pclm_count = sum(len(aligned_query_genomes[query_genome]) for query_genome in aligned_query_genomes)
    roast_format_tracker = tracker.subTracker("Format MAFs for ROAST",pclm_count,estimate=True,style="fraction")
    for query_genome in aligned_query_genomes:
        for maf_name in aligned_query_genomes[query_genome]:
            roast_format_tracker.step()
            out_name = os.path.join(out_folder,"temp.sing.maf")
            chrom, num_entries = prefix_and_get_chrom_and_count(maf_name,out_name,[reference,query_genome])
            new_name = out_name_template.format(chrom=chrom,query_genome=query_genome)
            os.rename(out_name,new_name)
            if num_entries>0: #We dont need to do anything with the empty files!
                if chrom not in per_chrom_labeled_mafs: 
                    per_chrom_labeled_mafs[chrom] = [] #a set here because we need dont want duplicates from the next step
                    per_chrom_existing_genomes[chrom] = []
                per_chrom_labeled_mafs[chrom].append(new_name)
                per_chrom_existing_genomes[chrom].append(query_genome)
    roast_format_tracker.done()

    roast_commandlists = []
    roast_files = {}
    currentwd = os.getcwd()
    os.chdir(out_folder)
    for chrom in per_chrom_labeled_mafs:
        outfile = chrom+".roast.maf"
        if len(per_chrom_existing_genomes[chrom])==len(aligned_query_genomes):
            chr_tree = '"%s"'%(tree.replace("*",chrom).replace(","," "))
            file_glob = chrom+".*.sing.maf"
            roast_commandlists.append(["roast",'E="%s"'%chrom,"X=0", chr_tree, file_glob, outfile])
            roast_files[chrom] = os.path.join(out_folder,outfile)
    random.shuffle(roast_commandlists)
    call_commands_async(roast_commandlists,num_processes,shell=True,parent=tracker,tracker_name="ROAST",env=cmd_env) #runs commands asynchronously with a maximum simultanious process count
    os.chdir(currentwd)

    msa_format_tracker = tracker.subTracker("Format MAFs for msa_view",len(roast_files),estimate=True,style="fraction")
    prepared_for_msa = {}
    for chrom in roast_files:
        msa_format_tracker.step()
        maf_name = roast_files[chrom]
        out_name = os.path.splitext(maf_name)[0]+".nochrom.maf"
        chrom,num_entries = remove_target_chrom_get_target_and_count(maf_name,out_name,reference)
        if num_entries > 0:
            prepared_for_msa[chrom] = out_name
    msa_format_tracker.done()
    
    msa_codon_commandlists = []
    codon_4d_names = {}
    gff_split_tracker = tracker.subTracker("Splitting Reference GFF",1,estimate=False,style="noProg")
    chrom_gffs = split_gff(ref_genome_gff,out_folder,cds_seqtypes)
    gff_split_tracker.done()
    for chrom in prepared_for_msa:
        if chrom in chrom_gffs:
            maf_name = prepared_for_msa[chrom]
            out_name = os.path.join(out_folder,chrom+".4d-codons.ss")
            msa_codon_commandlists.append(["msa_view","--in-format","MAF","--4d",maf_name,"--features",chrom_gffs[chrom],">",out_name])
            codon_4d_names[chrom] = out_name
    call_commands_async(msa_codon_commandlists,num_processes,shell=True,parent=tracker,tracker_name="msa_view step 1",env=cmd_env)

    msa_site_commandlists = []
    site_4d_names = {}
    for chrom in codon_4d_names:
        codon_ss_name = codon_4d_names[chrom]
        out_name = os.path.join(out_folder,chrom+".4d-sites.ss")
        msa_site_commandlists.append(["msa_view",codon_ss_name,"--in-format","SS","--out-format","SS","--tuple-size","1",">",out_name])
        site_4d_names[chrom] = out_name
    call_commands_async(msa_site_commandlists,num_processes,shell=True,parent=tracker,tracker_name="msa_view step 2",env=cmd_env)
    
    #fix this!
    site_4d_names = {key:site_4d_names[key] for key in site_4d_names if os.path.isfile(site_4d_names[key]) and os.stat(site_4d_names[key]).st_size>20}
    ss_folder = os.path.join(out_folder,"anonymized_ss")
    call_command(["mkdir","-p",ss_folder],shell=True,env=cmd_env)
    to_combine = []
    ss_mod_func.names = [name for name in aligned_query_genomes]
    for chrom in site_4d_names:
        perchrom_ss = ss.Handler(site_4d_names[chrom])
        out_name = os.path.join(ss_folder,chrom+".anon.ss")
        ss_mod_func.first = True
        perchrom_ss.modify(ss_mod_func,out_name)
        to_combine.append(out_name)

    combined_ss_file = os.path.join(out_folder,"combined.ss")
    aggregate_names = ",".join(["REFERENCE"]+ss_mod_func.names)
    call_command(["msa_view","--unordered-ss","--aggregate",aggregate_names,"--out-format","SS"]+to_combine+[">",combined_ss_file],shell=True,env=cmd_env,parent=tracker,tracker_name="Combining Models")
    out_root = os.path.join(out_folder,"conservation_model")
    chr_tree = '"%s"'%(tree.replace("*","REFERENCE"))
    pf_command = ["phyloFit","--min-informative","10","--tree",chr_tree,"--msa-format","SS","--out-root",out_root,combined_ss_file]
    call_command(pf_command,shell=True,env=cmd_env,parent=tracker,tracker_name="Running PhyloFit")
    consv_model = out_root+".mod"

    phastcons_commandlists = []
    chrom_beds = {}
    chrom_wigs = {}

    temp_model_folder = os.path.join(out_folder,"temp_models")
    call_command(["mkdir","-p",temp_model_folder],shell=True,env=cmd_env)
    e_rho_tree_folder = os.path.join(out_folder,"rho.estimated.trees")
    call_command(["mkdir","-p",e_rho_tree_folder],shell=True,env=cmd_env)

    for chrom in prepared_for_msa:
        maf_name = prepared_for_msa[chrom]

        bed_out = os.path.join(out_folder,chrom+".most-cons.bed")
        chrom_beds[chrom] = bed_out
        wig_out = os.path.join(out_folder,chrom+".scores.wig")
        chrom_wigs[chrom] = wig_out
        
        temp_model = os.path.join(temp_model_folder,chrom+".mod")
        create_model = ['sed','"s/REFERENCE/'+chrom+'/"',consv_model,'>',temp_model]
        delete_model = ['rm',temp_model]
        phast_opts = ["--target-coverage","0.25","--expected-length","12","--rho","0.4"]
        phast_command = ["phastCons"]+phast_opts+["--estimate-rho",e_rho_tree_folder, "--msa-format","MAF",maf_name,temp_model,"--most-conserved", bed_out,"--score","--seqname",chrom,">",wig_out]
        phastcons_commandlists.append(create_model+[";"]+phast_command+[";"]+delete_model)
    call_commands_async(phastcons_commandlists,num_processes,shell=True,parent=tracker,tracker_name="Running PhastCons",env=cmd_env)
    
    results = {}
    results["chrom_data"] = {chrom:{"chrom_seq_maf":roast_files[chrom],
                                    "chrom_conservation_wig":chrom_wigs[chrom],
                                    "chrom_conserved_bed":chrom_beds[chrom]} for chrom in chrom_wigs}
    tracker.freeze()
    return results

def ss_mod_func(entry):
    if ss_mod_func.first:
        ss_mod_func.first = not ss_mod_func.first
        entry.NAMES = [name if name in ss_mod_func.names else "REFERENCE" for name in entry.NAMES]

def prefix_and_get_chrom_and_count(maf_name,out_maf,names):
    with open(maf_name) as maf, open(out_maf,"w") as out:
        a_count = 0
        s_count = -1
        chrom = None
        for line in maf:
            if line.startswith("a"):
                s_count = 0
                a_count+= 1
            if line.startswith("s"):
                line_arr = line.split()
                if s_count==0: chrom = line_arr[1]
                else: line_arr[1] = "%s:%s" % (names[1],line_arr[1])
                line = " ".join(line_arr)+"\n"
                s_count+=1
            out.write(line)
    return chrom,a_count

def remove_target_chrom_get_target_and_count(maf_name,out_maf,ref_name):
    with open(maf_name) as maf, open(out_maf,"w") as out:
        a_count = 0
        s_count = 0
        chrom = None
        for line in maf:
            if line.startswith("a"):
                s_count = 0
                a_count+= 1
            if line.startswith("s"):
                line_arr = line.split()
                if line_arr[1].startswith(ref_name): #if it is the reference, just use the chromosome name
                    line_arr[1] = line_arr[1].split(":",1)[1]
                else:                                #otherwise, use only the species identifier (ignore chromosome).
                    line_arr[1] = line_arr[1].split(":",1)[0]
                if s_count==0: chrom = line_arr[1]
                s_count+= 1
                line = " ".join(line_arr)+"\n"
            out.write(line)
    return chrom,a_count

def split_gff(ref_genome_gff,out_foler,cds_seqtypes):
    out_form = os.path.join(out_foler,"{chrom}.gff")
    file_paths = {}
    with open(ref_genome_gff) as gff:
        for line in gff:
            if line.startswith("##FASTA"): break
            elif line.startswith("#") or line.strip()=="": continue
            list = line.split()
            chrom = list[0].strip()
            seqtype = list[2].strip()
            list[1] = "splitGFF" #insure identical source
            if seqtype in cds_seqtypes:
                list[2] = "CDS"
            if not chrom in file_paths:
                file_paths[chrom] = out_form.format(chrom=chrom)
                open(file_paths[chrom],"w").close()
            with open(file_paths[chrom],"a") as file:
                out = "\t".join(list)
                file.write(out if out.endswith("\n") else out+"\n")
    return file_paths

def config_score(config_path): 
    '''Runs the `score` function but loads arguements from a config file and uses the directory of that file as the working directory. It also saves a file called "align.results.json" in that directory.'''
    with open(config_path) as config_file:
        config = json.loads(config_file.read())
    original_wd = os.getcwd()
    config_directory = os.path.dirname(os.path.abspath(config_path))
    os.chdir(config_directory)
    score_results = score(**config)
    # combine results dict with config and output as JSON
    score_results = copy.deepcopy(config).update(score_results)
    results_path = os.path.join(config_directory,"score.results.json")
    with open(results_path,"w") as results_file:
        json.dump(score_results,results_file,sort_keys=True,indent=4)
    os.chdir(original_wd)
    
_cl_entry = config_score #function that should be run on command line entry to this subcommand
def _parser(parser_add_func,name):
    p = parser_add_func(name,description="Aligns a genome to a reference")
    p.add_argument("config_path", help="Absolute(!) path to config file")
    return p
