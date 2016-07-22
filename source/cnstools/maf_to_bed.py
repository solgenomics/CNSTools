from _filetypes import Maf, Bed6
import sys

import argparse

def parser(parser_add_func,name):
    p = parser_add_func(name,description="Converts a maf file into a bed file. If no output path is provided, outputs to stdout.")
    p.add_argument("maf_file", help="Path to maf file.")
    p.add_argument("-o", "--bed_out", default=None, help="Path to output (bed file). If not provided, program outputs to stdout.")
    p.add_argument("-r", "--ref_genome", default=None, help="Reference genome in the maf file to convert locations for.")
    p.add_argument("-t", "--index_tag", default=None, help="Adds TAG=INDEX to the bed entry name where INDEX is the index of maf entry in the input file and TAG is the provided string.")
    return p

def run(maf_file,bed_out=None,ref_genome=None,index_tag=None):
    """Converts a maf file into a bed file. If no output path is provided, outputs to stdout.

    :param string maf_file: Path to maf file.
    :param  string bed_out: Path to output (bed file). If not provided, program outputs to stdout.
    :param  string ref_genome: Reference genome in the maf file to convert locations for.
    :param  string index_tag: Adds TAG=INDEX to the bed entry name where INDEX is the index of maf entry in the input file and TAG is the provided string.
    :returns:  `None`
    """ 
    kwargs = {}
    if ref_genome: kwargs["genome_name"] = ref_genome
    if index_tag: kwargs["index_tag"] = index_tag
    bed = Maf(file_name=maf_file).to_bed(**kwargs)
    if bed_out:
        bed.save_file(bed_out)
    else:
        sys.stdout.write("\n".join(bed.get_lines())+"\n")
        sys.stdout.flush()