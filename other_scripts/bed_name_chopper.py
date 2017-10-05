import cnstools.file_handlers.bed6 as bed
import sys

def mod(entry):
    entry.chrom = entry.chrom.split(":")[-1]
    return entry

bed.Handler(sys.argv[1]).modify(mod, 
        path=sys.argv[2], 
        add_end_break=False)