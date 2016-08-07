import json
from filetype_classes import Cns
from _utils import create_path, JSON_saver, safe_print, header_print, Progress_tracker
import argparse

def _main(data,output_folder,overwrite=False):
    datasaver = JSON_saver(create_path(output_folder,"record","json",overwrite=overwrite))
    datasaver.save(data)

    header_print("Combining %s CNS files"%len(data['ref_aligned_chroms']))

    data["combined_cns"] = create_path(output_folder,"combined_identified","cns",overwrite=overwrite)
    print len(data['ref_aligned_chroms'])
    cns = Cns()
    i = 0
    for chrom in data['ref_aligned_chroms']:
        i+=1
        print i
        chrom_cns = Cns(file_name=data['ref_aligned_chroms'][chrom]['results'])
        for entry in chrom_cns.entries:
            entry.cns_ID = "%s:%s" %(chrom,str(entry.cns_ID))
            for genome in entry.sequences:
                for seq in entry.sequences[genome]:
                    seq.cns_ID = entry.cns_ID
            cns.entries.append(entry)
    cns.save_file(data["combined_cns"])
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
    _main(config,output_folder,overwrite)
