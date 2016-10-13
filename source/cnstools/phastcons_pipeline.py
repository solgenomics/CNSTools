import os
import json
from alignment_pipeline import call_commands_async

def run():
    per_genome_input_mafs = {"name":["filename"]}
    per_genome_input_gffs = {"name":["filename"]}
    roast_path = "/home/dal333/alignment_step/multiz-tba.012109/x86_64/bin/roast"
    reference = "Metru"
    msa_view_path = "msa_view"
    tree = ""
    num_processes = 9

    per_chrom_labeled_mafs = {}
    for genome in per_genome_input_mafs:
        for maf_name in per_genome_input_mafs[genome]:
            out_name = os.path.join(out,os.path.splitext(os.path.basename(maf_name))[0]+".prefixed.maf")
            chrom, num_entries = prefix_and_get_chrom_and_count(maf_name,out_name,[reference,genome])
            if num_entries>0: #We dont need to do anything with the empty files!
                if chrom not in per_chrom_labeled_mafs: per_chrom_labeled_mafs[chrom] = []
                per_chrom_labeled_mafs[chrom].append(out_name)

    roast_commandlists = []
    for chrom in per_chrom_labeled_mafs:
        folder = os.path.dirname(per_chrom_labeled_mafs[chrom][0])
        new_names = [os.path.join(out,".".join(chrom,os.path.basename(maf_name))) for maf_name in per_chrom_labeled_mafs[chrom]]
        for old_name,new_name in zip(per_chrom_labeled_mafs[chrom],new_names):
            os.rename(old_name,new_name)
        per_chrom_labeled_mafs[chrom] = new_names
        roast_commandlists.append([roast_path,"E="+chrom,"X=0", tree, chrom+".*.prefixed.maf", os.path.join(out,chrom+".roast.maf")])
    roast_files = [l[-1] for l in roast_commandlists]
    call_commands_async(roast_commandlists,num_processes,shell=True,tracker_name="roast") #runs commands asynchronously with maximum simultanious process count

    prepared_for_msa = []
    for maf_name in roast_files:
        out_name = os.path.splitext(maf_name)+".qchrom.maf"
        num_entries = remove_target_chrom_and_count(maf_name,out_name,reference)
        if num_entries > 0:
            prepared_for_msa.append(out_name)

    mas_commandlists = []
    for file in prepared_for_msa:
        genome = "???"

    "msa_view","--in-format","MAF","--4d","$file3","--features","$outdir/$ref.gff",">","$outdir/$ref.4d-codons.ss"
    call_commands_async(mas_commandlists,num_processes,shell=True,tracker_name="msa")


def prefix_and_get_chrom_and_count(maf_name,out_maf,names):
    with open(maf_name) as maf, open(out_maf,"w") as out:
        a_count = 0
        s_count = -1
        for line in maf:
            if line.startswith("a"):
                s_count = 0
                a_count+= 1
            if line.startswith("s") and s_count < len(names):
                line_arr = line.split()
                if s_count==0: chrom = line_arr[1]
                line_arr[1] = "%s:%s" % (names[s_count],line_arr[1])
                line = " ".join(line_arr)+"\n"
                s_count+=1
            out.write(line)
    return chrom,a_count

def remove_target_chrom_and_count(maf_name,out_maf,ref_name):
    with open(maf_name) as maf, open(out_maf,"w") as out:
        a_count = 0
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
                line = " ".join(line_arr)+"\n"
            out.write(line)
    return a_count

run()