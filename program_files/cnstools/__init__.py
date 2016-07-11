'''Responsible for importing and adding to package __all__ list all task modules, the modules listed in the tasks list will be the only ones that a user of the commandline tool can specify'''

#A list of all task modules to be imported and added to __all__. If present, the first instance of "task_" will be removed from the name.
tasks = [
    "task_analyze_cns",
    "task_bed_maf_parse",
    "task_blast_to_bed",
    "task_cns_to_fasta",
    "task_gff3_to_bed",
    "task_identify",
    "task_maf_to_bed",
    "task_maf_to_fasta",
    "task_parse_cns_data",
]
names = [file.replace("task_","",1) for file in tasks]
__all__ = names
task_name_pairs = zip(tasks,__all__)
for task,name in task_name_pairs:
    exec "import %s as %s"%(task,name)