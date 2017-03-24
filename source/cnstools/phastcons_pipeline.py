import os
import json
from _utils import call_commands_async
import random

config_defaults = {
    #"aligned_query_genomes": {"name":["chr_path","chr_path]},
    #"ref_genome_gff": "path",
    #"reference": "Metru",
    #"tree": "",
    "out_folder":".",
    "multiz_bin_path": "",
    "phast_bin_path": "",
    "msa_view_bin_path": "",
    "num_processes": 1,
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


    reference =     config["reference"]
    tree =          config["tree"]
    num_processes = config["num_processes"]

    aligned_query_genomes = config["aligned_query_genomes"]
    for genome in aligned_query_genomes:
        aligned_query_genomes[genome][:] = (os.path.abspath(p) for p in aligned_query_genomes[genome])

    ref_genome_gff =   os.path.abspath(config["annotations"][reference])
    out_folder =        os.path.abspath(config["out_folder"])
    multiz_bin_path =   os.path.abspath(config["multiz_bin_path"])
    msa_view_bin_path = os.path.abspath(config["msa_view_bin_path"])
    phast_bin_path =    os.path.abspath(config["phast_bin_path"])

    try:
        os.makedirs(out_folder)
    except OSError:
        if not os.path.isdir(out_folder):
            raise
    os.chdir(out_folder)
    cmd_env = os.environ.copy()

    if multiz_bin_path!="":
            cmd_env["PATH"] = multiz_bin_path+":" + cmd_env["PATH"]
    if phast_bin_path!="":
            cmd_env["PATH"] = phast_bin_path+":" + cmd_env["PATH"]

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
    #call_commands_async(roast_commandlists,num_processes,shell=True,tracker_name="roast",env=cmd_env) #runs commands asynchronously with a maximum simultanious process count

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
    #call_commands_async(msa_codon_commandlists,num_processes,shell=True,tracker_name=None,env=cmd_env)

    print("MSA_SITE")
    msa_site_commandlists = []
    site_4d_names = {}
    for chrom in codon_4d_names:
        codon_ss_name = codon_4d_names[chrom]
        out_name = os.path.join(out_folder,chrom+".4d-sites.ss")
        site_4d_names[chrom] = out_name
        msa_site_commandlists.append(["msa_view",codon_ss_name,"--in-format","SS","--out-format","SS","--tuple-size","1",">",out_name])
    #call_commands_async(msa_site_commandlists,num_processes,shell=True,tracker_name=None,env=cmd_env)
    
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

    config["chrom_data"] = {chrom:{"chrom_seq_maf":roast_files[chrom],
                                    "chrom_conservation_wig":chrom_wigs[chrom],
                                    "chrom_conserved_bed":chrom_beds[chrom]} for chrom in site_4d_names}
    with open("results.config.json","w") as out:
        json.dump(config,out,sort_keys=True,indent=4)



    os.chdir(original_wd)


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

def split_gff(ref_genome_gff,out_foler):
    out_form = os.path.join(out_foler,"{chrom}.gff")
    file_paths = {}
    with open(ref_genome_gff) as gff:
        for line in gff:
            if line.startswith("##FASTA"): break
            elif line.startswith("#") or line.strip()=="": continue
            chrom = line.split(None,1)[0].strip()
            if not chrom in file_paths:
                file_paths[chrom] = out_form.format(chrom=chrom)
            with open(file_paths[chrom],"a") as file:
                file.write(line)
    return file_paths
