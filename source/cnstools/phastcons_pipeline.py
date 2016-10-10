from file_handler import Maf_handler
import os

def run():
    per_genome_input_mafs = {"name":["filename"]}
    per_genome_labeled_mafs = {}
    for genome in per_genome_input_mafs:
        per_genome_labeled_mafs[genome] = []
        for maf_file in per_genome_input_mafs[genome]:
            name = ______
            Maf_handler(maf_file).modify_entries(add_genome_prefix(genome),path=name)
            per_genome_labeled_mafs[genome].append(name)
    for genome in per_genome_input_mafs:
        per_genome_labeled_mafs[genome] = []
        for maf_file in per_genome_input_mafs[genome]:

def add_genome_prefix(genome):
    first_sequence = True
    def modify_func(entry):
        if first_sequence: #the first sequence of each maf file is the 
            first_sequence = False
            return entry
        else:
            entry.name = "%s:%s"%(genome,entry.name.split(":")[-1])
            return entry
    return modify_func

run()
