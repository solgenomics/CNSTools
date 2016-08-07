from filetype_classes import Bed6,Wiggle
from _utils import get_docstring_info
import sys

import argparse

def parser(parser_add_func,name):
    p = parser_add_func(name,description=get_docstring_info(run))
    p.add_argument("wig_file", help=get_docstring_info(run,"wig_file"))
    p.add_argument("-o", "--out_file", default=None, help=get_docstring_info(run,"out_file"))
    p.add_argument("-l", "--min_seg_length", default=None, help=get_docstring_info(run,"min_seg_length"))
    p.add_argument("-s", "--min_seg_score", default=None, help=get_docstring_info(run,"min_seg_score"))
    p.add_argument("-L", "--rejection_seg_len", default=None, help=get_docstring_info(run,"rejection_seg_len"))
    p.add_argument("-S", "--rejection_seg_score", default=None, help=get_docstring_info(run,"rejection_seg_score"))
    return p

def run(wig_file,out_file=None,min_seg_length=7,min_seg_score=0.82,rejection_seg_len=12,rejection_seg_score=0.55,genome_name=""):
    """Converts a wig file into a bed file, filtered by conservation. If no output path is provided, outputs to stdout.

    :param string wig_file: Path to .wig file.
    :param string out_file: Path to the output .bed file.
    :param int min_seg_length: Minimum length for a returned region.
    :param float min_seg_score: Minimum average locus score for a returned region.
    :param int rejection_seg_len: Minimum size of low score region to be sliced out of a region.
    :param float rejection_seg_score: Maximum score that is considered a candidate for removal in a contigous segment of at least length 'rejection_seg_len'.
    """
    wig = Wiggle(file_name=wig_file)
    bed = wig.to_bed(min_seg_length,min_seg_score,rejection_seg_len,rejection_seg_score)
    if genome_name:
        for entry in bed.entries:
            entry.chrom = genome_name+":"+str(entry.chrom)
    if out_file:
        bed.save_file(out_file)
    else:
        sys.stdout.write("\n".join(bed.get_lines())+"\n")
        sys.stdout.flush()