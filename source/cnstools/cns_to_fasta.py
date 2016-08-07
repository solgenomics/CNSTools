from filetype_classes import Cns, Fasta
import sys

import argparse

def parser(parser_add_func,name):
    p = parser_add_func(name,description="Converts the set of sequences in a .cns file into a fasta file for a genome. If no output path is provided, outputs to stdout.")
    p.add_argument("cns_file", help="Path to .cns file.")
    p.add_argument("genome", help="Reference genome in the maf file to convert locations for.")
    p.add_argument("-o","--out_file", help="Path to maf file.")
    return p

def run(cns_file,genome,out_file=None):
    fasta = Cns(file_name=cns_file).to_fasta(genome)
    if out_file:
        fasta.save_file(out_file)
    else:
        sys.stdout.write("\n".join(fasta.get_lines())+"\n")
        sys.stdout.flush()