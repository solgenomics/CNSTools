import json
import os
from copy import deepcopy
from _filetypes import Maf, Gff3, Bed6, Bed13
from _utils import create_path, JSON_saver, safe_print, header_print, Progress_tracker
import argparse
import subprocess

from maf_to_bed import run as maf_to_bed
from gff3_to_bed import run as gff3_to_bed
from slice_maf_by_bed import run as slice_maf_by_bed

# json_file_format = {
#     "all_seq_maf":"PATH",
#     "all_conserved_bed":"PATH",
#     "ref_genome":"PREFIX",
#     "genomes":{
#         "PREFIX":{
#             "annot_gff3":"PATH"
#         }
#     }
# }

def main(parameter_dict,output_folder,num_threads,overwrite=False):
    datasaver = JSON_saver(create_path(output_folder,"record","json",overwrite=overwrite))
    data = deepcopy(parameter_dict)
    datasaver.save(data)

    #maf_to_bed
    info = "Convert aligned sequences to .bed:"
    header_print(info)
    data['all_seq_bed'] = create_path(output_folder,"all_seq","bed",overwrite=overwrite)
    maf_to_bed(maf_file    = data['all_seq_maf'],
               bed_out     = data['all_seq_bed'],
               ref_genome  = data['ref_genome'],
               index_tag   = "all_maf_index")
    datasaver.save(data)

    #$bedtools intersect
    info = "Intersect aligned regions with conserved regions:"
    header_print(info)
    data['conserved_bed'] = create_path(output_folder,"conserved","bed",overwrite=overwrite)
    cmd = "bedtools intersect -a %s -b %s > %s" % (data['all_seq_bed'],data['all_conserved_bed'],data['conserved_bed'])
    tracker = Progress_tracker("Running bedtools intersect",1).estimate(False).display()
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    tracker.done()
    datasaver.save(data)

    #gff3_to_bed
    info = "Convert coding regions to .bed:"
    header_print(info)
    data['coding_bed'] = create_path(output_folder,"coding","bed",overwrite=overwrite)
    gff3_to_bed(gff3_file   = data['genomes'][data['ref_genome']]['annot_gff3'],
                bed_out     = data['coding_bed'],
                type_list   = ['CDS'],
                sequence_prefix = data['ref_genome']+":")
    datasaver.save(data)

    #$bedtools subtract
    info = "Subtract coding regions from conserved aligned regions:"
    header_print(info)
    data['cns_bed'] = create_path(output_folder,"cns","bed",overwrite=overwrite)
    cmd = "bedtools subtract -a %s -b %s > %s" % (data['conserved_bed'],data['coding_bed'],data['cns_bed'])
    tracker = Progress_tracker("Running bedtools subtract",1).estimate(False).display()
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    tracker.done()
    datasaver.save(data)

    #slice_maf_by_bed
    info = "Slice multi-alignment file based on identified conserved non-coding regions:"
    header_print(info)
    data['cns_maf'] = create_path(output_folder,"cns","maf",overwrite=overwrite)
    slice_maf_by_bed(maf_file       = data['all_seq_maf'],
                     bed_file       = data['cns_bed'],
                     index_tag      = "all_maf_index",
                     ref_genome     = data['ref_genome'],
                     out_file       = data['cns_maf'],
                     max_N_ratio    = 0.5,
                     max_gap_ratio  = 0.5,
                     min_len        = 15)
    datasaver.save(data)

    #maf_to_bed
    info = "Convert per-genome CNS regions to .bed:"
    header_print(info)
    data['genome_cns_beds_folder'] = create_path(output_folder+"genome_cns_beds",overwrite=overwrite)
    cns_maf = Maf(file_name=data['cns_maf'])
    for genome in data['genomes']:
        data['genomes'][genome]['cns_bed'] = create_path(data['genome_cns_beds_folder'],genome+"_cns","bed",overwrite=overwrite)
        bed = cns_maf.to_bed(genome_name=genome,index_tag="cns_maf_index")
        bed.save_file(data['genomes'][genome]['cns_bed'])
    del cns_maf
    datasaver.save(data)

    #gff3_to_bed
    info = "Convert per-genome gene regions to .bed:"
    header_print(info)
    data['genome_annot_beds_folder'] = create_path(output_folder+"genome_annot_beds",overwrite=overwrite)
    for genome in data['genomes']:
        data['genomes'][genome]['annot_bed'] = create_path(data['genome_annot_beds_folder'],genome+"_annot","bed",overwrite=overwrite)
        Gff3(file_name=data['genomes'][genome]['annot_gff3']) \
            .to_bed(type_list=['gene'],genome=genome) \
            .save_file(data['genomes'][genome]['annot_bed'])
    datasaver.save(data)

    #$bedtools closest
    info = "Find closest gene for each CNS region:"
    header_print(info)
    data['gene_proximity_beds_folder'] = create_path(output_folder+"gene_proximity_beds",overwrite=overwrite)
    for genome in data['genomes']:
        data['genomes'][genome]['gene_proximity_bed'] = \
            create_path(data['gene_proximity_beds_folder'],genome+"_proxim","bed",overwrite=overwrite)
        cmd = "bedtools closest -D a -a %s -b %s > %s" % \
            (data['genomes'][genome]['cns_bed'],
             data['genomes'][genome]['annot_bed'],
             data['genomes'][genome]['gene_proximity_bed'])
        process = subprocess.Popen(cmd, shell=True)
        process.wait()
    datasaver.save(data)

    #maf_and_proxim_bed_to_cns
    info = "Process proximity and maf files into .cns file:"
    header_print(info)
    data['results'] = create_path(output_folder,"identified_CNSs","cns",overwrite=overwrite)
    cns_proxim_beds = {genome:Bed13(data['genomes'][genome]['gene_proximity_bed']) for genome in data['genomes']}
    Maf(file_name=data['cns_maf'])\
        .cns_from_proxim_beds(cns_proxim_beds,"cns_maf_index")\
        .save_file(data['results'])
    datasaver.save(data)
    
parser = argparse.ArgumentParser("identify")
parser.add_argument("parameter_file", help="Path to the parameter file.")
parser.add_argument("-o", "--output_folder", default="./cnstools_out/", help="Path to output folder.")
parser.add_argument("-t", "--num_threads", type=int, default=1, help="Reference sequence in the maf file to take location information from.")
parser.add_argument("-f", "--overwrite", action='store_true', help="If present, the program overwrites data in the output folder.")

def run(parameter_file,output_folder,num_threads,overwrite=False):
    config = None
    with open(parameter_file) as intructionJSON:
        config = json.load(intructionJSON)
    output_folder = create_path(output_folder,overwrite=overwrite)
    main(config,output_folder,num_threads,overwrite)
    
        
