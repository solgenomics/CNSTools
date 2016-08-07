from filetype_classes import Gff3, Bed6 
import sys
from _utils import Progress_tracker

import argparse

def parser(parser_add_func,name):
    p = parser_add_func(name,description="Converts a gff3 file into a bed file. If no output path is provided, outputs to stdout.")
    p.add_argument("gff3_file", help="Path to gff3 file.")
    p.add_argument("-o", "--bed_out", default=None, help="Path to output. (bed file)")
    p.add_argument("-t", "--type_list", nargs='*', default=[], help="Specifies that only sequences of these types should be added to the bed.")
    p.add_argument("-p", "--sequence_prefix", nargs='*', default=[], help="Appends provided string to the front of the chromosome names in the bed file.")
    return p

def run(gff3_file,bed_out=None,type_list=[],sequence_prefix=None):
    """Converts a gff3 file into a bed file. If no output path is provided, outputs to stdout.

    :param string gff3_file: Path to gff3 file.
    :param  string bed_out: Path to output. (bed file)
    :param list[string] type_list: Specifies that only sequences of these types should be added to the bed.
    :param string sequence_prefix: Appends provided string to the front of the chromosome names in the bed file.
    :returns:  `None`
    """ 
    bed = Gff3(file_name=gff3_file).to_bed(type_list)
    if sequence_prefix:
        tracker = Progress_tracker("Prefixing bed chrom IDs",len(bed.entries)).auto_display().start()
        for entry in bed.entries:
            entry.chrom = sequence_prefix+entry.chrom
            tracker.step()
        tracker.done()
    if bed_out:
        bed.save_file(bed_out)
    else:
        sys.stdout.write("\n".join(bed.get_lines())+"\n")
        sys.stdout.flush()