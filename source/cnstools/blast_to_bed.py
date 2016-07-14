"""Script to convert a blast output file (outfmt 6) to a .bed file.
"""
from _filetypes import BlastF6, Bed6 

def main(blast_file,outFile):
    """Runs the main workflow for the script

    :param str blast_file: Path to the blast results to be converted.
    :param str out_file: Path where the bed file should be saved.
    :returns: `None`
    """
    with open(blast_file) as f, open(outFile,"w") as out:
        blast = BlastF6(lines=f.readlines())
        bed = blast.to_bed()
        lines = bed.get_lines()
        out.write("\n".join(lines))

def file_run(blast_file,outFile):
    """ Recieves parameters and forwards them to :func:`.main`.
    
    :param str blast_file: Path to the blast results to be converted.
    :param str out_file: Path where the bed file should be saved.
    :returns: `None`
    """
    main(blast_file,outFile)