from filetype_classes import Cns, Bed6
import sys

import argparse

def parser(parser_add_func,name):
    p = parser_add_func(name,description="Converts the set of sequences in a .cns file into a bed file for a genome. If no output path is provided, outputs to stdout.")
    p.add_argument("cns_file", help="Path to .cns file.")
    p.add_argument("genome", help="Reference genome in the maf file to convert locations for.")
    p.add_argument("-o","--out_file", help="Path to maf file.")
    return p

def run(cns_file,genome,out_file=None):
    bed = Cns(file_name=cns_file).to_bed(genome)
    if out_file:
        bed.save_file(out_file)
    else:
        sys.stdout.write("\n".join(bed.get_lines())+"\n")
        sys.stdout.flush()