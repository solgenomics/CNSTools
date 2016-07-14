
import _progress_tracker as pt
from _filetypes import Maf, Bed6

def main(mafFile,outFile):
    with open(mafFile) as file, open(outFile,"w") as out:
        maf = Maf(lines=file.readlines())
        bed = maf.to_bed()
        for i in range(len(bed.entries)):
            id_string = "ID=%s" % i
            if bed.entries[i].name:
                bed.entries[i].name += ";"+id_string
            else:
                bed.entries[i].name = id_string

        lines = bed.get_lines()
        out.write("\n".join(lines))

def file_run(mafFile,outFile): 
    main(mafFile,outFile)