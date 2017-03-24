from _utils import create_path, JSON_saver, safe_print, header_print, Progress_tracker
from copy import deepcopy
from filetype_classes import Maf, Gff3, Bed6, Bed13, Cns
from gff3_to_bed import run as gff3_to_bed
from maf_to_bed import run as maf_to_bed
from slice_maf_by_bed import run as slice_maf_by_bed
from wiggle_to_bed import run as wiggle_to_bed
import argparse
import json
import os
import subprocess

config_defaults = {

}

# data = {
#     "ref_aligned_chroms":{
#         "CHROM":{
#             "chrom_seq_maf":"PATH",
#             "chrom_conservation_wig":"PATH"
#         }
#     },
#     "ref_genome":"PREFIX",
#     "genomes":{
#         "PREFIX":{
#             "annot_gff3":"PATH"
#         }
#     }
# }
# 

def parser(parser_add_func,name):
    p = parser_add_func(name)
    p.add_argument("config_path", help="Path to the configuration file.")
    return p

def run(config_path):

    with open(config_path) as config_file:
      config = json.loads(config_file.read())
    for key in config_defaults:
        config.setdefault(key, config_defaults[key])

    os.chdir(os.path.dirname(os.path.abspath(config_path)))

    out_folder = os.path.abspath(config["out_folder"])

    try:
        os.makedirs(out_folder)
    except OSError:
        if not os.path.isdir(out_folder):
            raise
    os.chdir(out_folder)

    header_print("Running full CNS identification pipeline on %s alignment files" % len(config["chrom_data"]),h_type=1)
    create_genome_beds(config)

    for chromosome in config['chrom_data']:
        header_print("Identifying CNS on %s" % chromosome,h_type=2)
        config['chrom_data'][chromosome]['chrom_out'] = os.path.join(out_folder,chromosome+"_cns")
        try:
            os.makedirs(config['chrom_data'][chromosome]['chrom_out'])
        except OSError:
            if not os.path.isdir(config['chrom_data'][chromosome]['chrom_out']):
                raise
        chrom_cns_identify(config,chrom_name=chromosome)


def create_genome_beds(config):
    gb_out = os.path.abspath("./genome_beds")
    try:
        os.makedirs(gb_out)
    except OSError:
        if not os.path.isdir(gb_out):
            raise

    #gff3_to_bed
    info = "Convert coding regions to .bed:"
    header_print(info)
    config["reference_coding_bed"] = os.path.join(gb_out,"ref_coding.bed")
    gff3_to_bed(gff3_file   = config['annotations'][config['reference']],
                bed_out     = config["reference_coding_bed"],
                type_list   = config["coding_features"],
                sequence_prefix = config['reference']+":")

    #gff3_to_bed
    info = "Convert per-genome gene regions to .bed:"
    header_print(info)
    config["genome_gene_beds"] = {}
    for genome in config['annotations']:
        config["genome_gene_beds"][genome] = os.path.join(gb_out,genome+"_genes.bed")
        Gff3(file_name=config['annotations'][genome]) \
            .to_bed(type_list=['gene'],genome=genome) \
            .save_file(config["genome_gene_beds"][genome])



def chrom_cns_identify(config,chrom_name):

    chrom = config['chrom_data'][chrom_name]

    #maf_to_needed_format
    chrom['formatted_chrom_seq_maf'] = os.path.join(chrom["chrom_out"],"chrom_seqs_formatted.maf")
    maf = Maf(file_name=chrom['chrom_seq_maf'])
    for entry in maf.entries:
        for sequence in entry.sequences:
            if sequence.src.count(":")<1:
                sequence.src = ":".join((config['reference'],sequence.src))
    maf.save_file(chrom['formatted_chrom_seq_maf'])

    #maf_to_bed
    info = "Convert aligned sequences to .bed:"
    header_print(info)
    chrom['ref_seq_bed'] = os.path.join(chrom["chrom_out"],"ref_seq.bed")
    maf_to_bed(maf_file    = chrom['formatted_chrom_seq_maf'],
               bed_out     = chrom['ref_seq_bed'],
               ref_genome  = config['reference'],
               index_tag   = "chrom_maf_index")

    #$bedtools subtract
    info = "Subtract coding regions from aligned regions:"
    header_print(info)
    chrom['aligned_noncoding_bed'] = os.path.join(chrom["chrom_out"],"aligned_noncoding_bed.bed")
    cmd = "bedtools subtract -a %s -b %s > %s" % (chrom['ref_seq_bed'],config['reference_coding_bed'],chrom['aligned_noncoding_bed'])
    tracker = Progress_tracker("Running bedtools subtract",1).estimate(False).display()
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    tracker.done()

    #wiggle_to_bed
    info = "Converting especially conserved regions in wiggle file to bed"
    header_print(info)
    chrom['best_conserved_bed'] = os.path.join(chrom["chrom_out"],"best_conserved.bed")
    wiggle_to_bed(wig_file=chrom['chrom_conservation_wig'],
                  out_file=chrom['best_conserved_bed'],
                  genome_name=config['reference'])

    #filter_bed_with_wiggle
    info = "Intersecting wiggle bed with the noncoding bed"
    header_print(info)
    chrom['cns_bed'] = os.path.join(chrom["chrom_out"],"cns.bed")
    cmd = "bedtools intersect -a %s -b %s > %s" % (chrom['aligned_noncoding_bed'],chrom['best_conserved_bed'],chrom['cns_bed'])
    tracker = Progress_tracker("Running bedtools intersect",1).estimate(False).display()
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    tracker.done()

    #slice_maf_by_bed
    info = "Slice multi-alignment file based on identified conserved non-coding regions:"
    header_print(info)
    chrom['cns_maf'] = os.path.join(chrom["chrom_out"],"cns.maf")
    slice_maf_by_bed(maf_file       = chrom['formatted_chrom_seq_maf'],
                     bed_file       = chrom['cns_bed'],
                     index_tag      = "chrom_maf_index",
                     ref_genome     = config['reference'],
                     out_file       = chrom['cns_maf'],
                     max_N_ratio    = config['max_cns_N_ratio'],#0.5,
                     max_gap_ratio  = config['max_cns_gap_ratio'],#0.5,
                     min_len        = config['min_cns_length'])#15)

    #for iteration over genomes regardless of reference
    all_genomes = config["genomes"]

    #maf_to_bed
    info = "Convert per-genome CNS regions to .bed:"
    header_print(info)
    chrom['genome_cns_beds_folder'] = os.path.join(chrom["chrom_out"],"genome_cns_beds")
    try:
        os.makedirs(chrom['genome_cns_beds_folder'])
    except OSError:
        if not os.path.isdir(chrom['genome_cns_beds_folder']):
            raise
    chrom['genome_cns_beds'] = {}
    cns_maf = Maf(file_name=chrom['cns_maf'])
    for genome in all_genomes:
        chrom['genome_cns_beds'][genome] = os.path.join(chrom['genome_cns_beds_folder'],genome+"_cns_"+chrom_name+".bed")
        bed = cns_maf.to_bed(genome_name=genome,index_tag="cns_maf_index")
        bed.save_file(chrom['genome_cns_beds'][genome])
    del cns_maf


    #$bedtools closest
    info = "Find closest gene for each CNS region:"
    header_print(info)
    chrom['gene_proximity_beds_folder'] = os.path.join(chrom["chrom_out"],"gene_proximity_beds")
    try:
        os.makedirs(chrom['gene_proximity_beds_folder'])
    except OSError:
        if not os.path.isdir(chrom['gene_proximity_beds_folder']):
            raise
    chrom['gene_proximity_beds'] = {}
    for genome in all_genomes:
        chrom['gene_proximity_beds'][genome] = os.path.join(chrom['gene_proximity_beds_folder'],genome+"_proxim.bed")
        cmd = "bedtools closest -D a -a %s -b %s > %s" % \
            (chrom['genome_cns_beds'][genome],
             config['genome_gene_beds'][genome],
             chrom['gene_proximity_beds'][genome])
        process = subprocess.Popen(cmd, shell=True)
        process.wait()

    #maf_and_proxim_bed_to_cns
    info = "Process proximity and maf files into .cns file:"
    header_print(info)
    chrom['results'] = os.path.join(chrom["chrom_out"],"identified_CNSs.cns")
    cns_proxim_beds = {genome:Bed13(chrom['gene_proximity_beds'][genome]) for genome in all_genomes}
    Maf(file_name=chrom['cns_maf'])\
        .cns_from_proxim_beds(cns_proxim_beds,"cns_maf_index",prefix=chrom_name+".")\
        .save_file(chrom['results'])    
        
def combine_cns():
    return None