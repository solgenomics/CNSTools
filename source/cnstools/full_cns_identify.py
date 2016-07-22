import json
from _filetypes import Gff3, Bed6
from _utils import create_path, JSON_saver, safe_print, header_print, Progress_tracker
import argparse
from copy import deepcopy

from create_genome_beds import _main as create_genome_beds
from chrom_cns_identify import _main as chrom_cns_identify

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

def _main(data,output_folder,num_threads,overwrite=False):
    data['out'] = output_folder
    datasaver = JSON_saver(create_path(data['out'],"record","json",overwrite=overwrite))
    datasaver.save(data)

    header_print("Running full CNS identification pipeline on %s chrom" % len(data["ref_aligned_chroms"]))
    data['genome_beds'] = create_path(data['out']+"genome_beds",overwrite=overwrite)
    data = create_genome_beds(data,data['genome_beds'],overwrite=overwrite)

    for chromosome in data['ref_aligned_chroms']:
        header_print("Identify CNS on %s" % chromosome)
        chromDat = {key:data[key] for key in data if key!="ref_aligned_chroms"}
        chromDat['chrom_seq_maf'] = data['ref_aligned_chroms'][chromosome]['chrom_seq_maf']
        chromDat['chrom_conserved_bed'] = data['ref_aligned_chroms'][chromosome]['chrom_conserved_bed']
        chromDat['out'] = create_path(data['out']+"chrom/"+chromosome,overwrite=overwrite)
        chromDat = chrom_cns_identify(chromDat,chromDat['out'],num_threads,overwrite=overwrite,chrom_name=chromosome)

        data['ref_aligned_chroms'][chromosome] = {key:chromDat[key] for key in chromDat if not key.startswith("ref_")}
        datasaver.save(data)

    return data

def parser(parser_add_func,name):
    p = parser_add_func(name)
    p.add_argument("parameter_file", help="Path to the parameter file.")
    p.add_argument("-o", "--output_folder", default="./cnstools_out/", help="Path to output folder.")
    p.add_argument("-t", "--num_threads", type=int, default=1, help="Reference sequence in the maf file to take location information from.")
    p.add_argument("-f", "--overwrite", action='store_true', help="If present, the program overwrites data in the output folder.")
    return p

def run(parameter_file,output_folder,num_threads,overwrite=False):
    config = None
    with open(parameter_file) as intructionJSON:
        config = json.load(intructionJSON)
    output_folder = create_path(output_folder,overwrite=overwrite)
    _main(config,output_folder,num_threads,overwrite=overwrite)
