'''Responsible for importing and adding to package __all__ list all task modules, the modules listed in the tasks list will be the only ones that a user of the commandline tool can specify'''

#A list of all task modules to be imported and added to __all__.
names = [
    "analyze_cns",
    "bed_maf_parse",
    "blast_to_bed",
    "cns_to_fasta",
    "gff3_to_bed",
    "identify",
    "maf_to_bed",
    "maf_to_fasta",
    "parse_cns_data",
]
__all__ = names
for name in names:
    exec "import %s"%(name)