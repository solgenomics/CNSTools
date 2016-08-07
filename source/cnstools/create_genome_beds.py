import json
from filetype_classes import Gff3, Bed6
from _utils import create_path, JSON_saver, safe_print, header_print, Progress_tracker
import argparse

from gff3_to_bed import run as gff3_to_bed

# data = {
#     "ref_aligned_chroms":{
#         "CHROM":{
#             "chrom_seq_maf":"PATH",
#             "chrom_conserved_bed":"PATH"
#         }
#     },
#     "ref_genome":"PREFIX",
#     "genomes":{
#         "PREFIX":{
#             "annot_gff3":"PATH"
#         }
#     }
# }

def _main(data,output_folder,overwrite=False):
    datasaver = JSON_saver(create_path(output_folder,"record","json",overwrite=overwrite))
    datasaver.save(data)

    #gff3_to_bed
    info = "Convert coding regions to .bed:"
    header_print(info)
    data['ref_coding_bed'] = create_path(output_folder,"ref_coding","bed",overwrite=overwrite)
    gff3_to_bed(gff3_file   = data['genomes'][data['ref_genome']]['annot_gff3'],
                bed_out     = data['ref_coding_bed'],
                type_list   = ['CDS'],
                sequence_prefix = data['ref_genome']+":")
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

    return data

def parser(parser_add_func,name):
    p = parser_add_func(name)
    p.add_argument("parameter_file", help="Path to the parameter file.")
    p.add_argument("-o", "--output_folder", default="./cnstools_out/", help="Path to output folder.")
    p.add_argument("-f", "--overwrite", action='store_true', help="If present, the program overwrites data in the output folder.")
    return p

def run(parameter_file,output_folder,overwrite=False):
    config = None
    with open(parameter_file) as intructionJSON:
        config = json.load(intructionJSON)
    output_folder = create_path(output_folder,overwrite=overwrite)
    _main(config,output_folder,num_threads,overwrite)
