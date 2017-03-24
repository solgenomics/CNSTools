import os, json, random
from _utils import call_commands_async

def score(reference,
          aligned_query_genomes,
          tree,
          annotations,
          msa_view_path = "",
          multiz_path   = "",
          num_processes = 1,
          out_folder    = "./",
          phast_path    = ""):
    '''ROASTs and scores files for conservation'''
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
    for query_genome in aligned_query_genomes:
        for maf_name in aligned_query_genomes[query_genome]:
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

    roast_commandlists = []
    roast_files = {}
    for chrom in per_chrom_labeled_mafs:
        outfile = os.path.join(out_folder,chrom+".roast.maf")
        if len(per_chrom_existing_genomes[chrom])==len(aligned_query_genomes):
            chr_tree = '"%s"'%(tree.replace("*",chrom).replace(","," "))
            roast_commandlists.append(["roast",'E="%s"'%chrom,"X=0", chr_tree, chrom+".*.sing.maf", outfile])
            roast_files[chrom] = outfile
    random.shuffle(roast_commandlists)
    for line in roast_commandlists: print line
    call_commands_async(roast_commandlists,num_processes,shell=True,tracker_name="roast",env=cmd_env) #runs commands asynchronously with a maximum simultanious process count

    prepared_for_msa = {}
    for chrom in roast_files:
        maf_name = roast_files[chrom]
        out_name = os.path.splitext(maf_name)[0]+".nochrom.maf"
        chrom,num_entries = remove_target_chrom_get_target_and_count(maf_name,out_name,reference)
        if num_entries > 0:
            prepared_for_msa[chrom] = out_name
            
    print("MSA_CODON")
    msa_codon_commandlists = []
    codon_4d_names = {}
    chrom_gffs = split_gff(ref_genome_gff,out_folder)
    for chrom in prepared_for_msa:
        if chrom in chrom_gffs:
            maf_name = prepared_for_msa[chrom]
            out_name = os.path.join(out_folder,chrom+".4d-codons.ss")
            msa_codon_commandlists.append(["msa_view","--in-format","MAF","--4d",maf_name,"--features",chrom_gffs[chrom],">",out_name])
            codon_4d_names[chrom] = out_name
    call_commands_async(msa_codon_commandlists,num_processes,shell=True,tracker_name=None,env=cmd_env)

    print("MSA_SITE")
    msa_site_commandlists = []
    site_4d_names = {}
    for chrom in codon_4d_names:
        codon_ss_name = codon_4d_names[chrom]
        out_name = os.path.join(out_folder,chrom+".4d-sites.ss")
        msa_site_commandlists.append(["msa_view",codon_ss_name,"--in-format","SS","--out-format","SS","--tuple-size","1",">",out_name])
        site_4d_names[chrom] = out_name
    call_commands_async(msa_site_commandlists,num_processes,shell=True,tracker_name=None,env=cmd_env)
    
    #fix this!
    site_4d_names = {key:site_4d_names[key] for key in site_4d_names if os.path.isfile(site_4d_names[key])}

    print("PHYLOFIT")
    phylofit_commandlists = []
    mod_names = {}
    for chrom in site_4d_names:
        site_4d_name = site_4d_names[chrom]
        chr_tree = '"%s"'%(tree.replace("*",chrom))
        out_root = os.path.join(out_folder,chrom+".nonconserved-4d2")
        phylofit_commandlists.append(["phyloFit","--min-informative","10","--tree",chr_tree,"--msa-format","SS","--out-root",out_root,site_4d_name])
        mod_names[chrom] = out_root+".mod"
    call_commands_async(phylofit_commandlists,num_processes,shell=True,tracker_name=None,env=cmd_env)

    phastcons_commandlists = []
    chrom_beds = {}
    chrom_wigs = {}

    for chrom in mod_names:
        mod_name = mod_names[chrom]
        maf_name = prepared_for_msa[chrom]

        bed_out = os.path.join(out_folder,chrom+".most-cons.bed")
        chrom_beds[chrom] = bed_out
        wig_out = os.path.join(out_folder,chrom+".scores.wig")
        chrom_wigs[chrom] = wig_out
        
        e_rho_tree_folder = os.path.join(out_folder,"rho.estimated.trees")
        call_commands_async([["mkdir","-p",e_rho_tree_folder]],1,shell=True,tracker_name=None,env=cmd_env)
        phast_opts = ["--target-coverage","0.25","--expected-length","12","--rho","0.4"]
        phastcons_commandlists.append(["phastCons"]+phast_opts+["--estimate-rho",e_rho_tree_folder,
                                       "--msa-format","MAF",maf_name,mod_name,"--most-conserved",
                                       bed_out,"--score","--seqname",chrom,">",wig_out])
    call_commands_async(phastcons_commandlists,num_processes,shell=True,tracker_name=None,env=cmd_env)
    
    results = {}
    results["chrom_data"] = {chrom:{"chrom_seq_maf":roast_files[chrom],
                                    "chrom_conservation_wig":chrom_wigs[chrom],
                                    "chrom_conserved_bed":chrom_beds[chrom]} for chrom in mod_names}
    return results

def config_score(config_path): 
    '''Runs the `score` function but loads arguements from a config file and uses the directory of that file as the working directory. It also saves a file called "align.results.json" in that directory.'''
    with open(config_path) as config_file:
        config = json.loads(config_file.read())
    original_wd = os.getcwd()
    config_directory = os.path.dirname(config_path)
    os.chdir(config_directory)
    score_results = score(**config)
    # combine results dict with config and output as JSON
    results = score_results.update(copy.deepcopy(config))
    results_path = os.path.join(config_directory,"score.results.json")
    with open(results_path) as results_file:
        json.dump(results,results_file,sort_keys=True,indent=4)
    os.chdir(original_wd)
    
_cl_entry = config_score #function that should be run on command line entry to this subcommand
def _parser(parser_add_func,name):
    p = parser_add_func(name,description="Aligns a genome to a reference")
    p.add_argument("config_path", help="Absolute(!) path to config file")
    return p