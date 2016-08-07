"""The cnstools package can be run in three ways. It can be run as an executable directly from the commandline, 
run as a python program with ``$python cnstools`` or it can be imported as a python module :class:`cnstools` and used with other 
python scripts. Currently this documentation only covers the imported module, more to come later."""

#A list of all task modules to be imported and added to __all__.
names = [
    "analyze_cns",
    "chrom_cns_identify",
    "create_genome_beds",
    "full_cns_identify",
    "gff3_to_bed",
    "maf_to_bed",
    "slice_maf_by_bed",
    "combine_cns",
    "wiggle_to_bed",
    "doc_testing",
    "cns_to_fasta",
    "cns_to_bed"
]
__all__ = names
for name in names:
    exec "import %s"%(name)