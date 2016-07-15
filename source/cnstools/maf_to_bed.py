
from _filetypes import Maf, Bed6

def main(mafFile,outFile,seq_name=None,index_tag="ID"):
    with open(mafFile) as file, open(outFile,"w") as out:
        maf = Maf(lines=file.readlines())
        bed = maf.to_bed(seq_name)
        for i in range(len(bed.entries)):
            id_string = "%s=%s" % (index_tag,i) if index_tag!=None else ""
            if bed.entries[i].name:
                bed.entries[i].name += ";"+id_string
            else:
                bed.entries[i].name = id_string

        lines = bed.get_lines()
        out.write("\n".join(lines))

def file_run(mafFile,outFile,seq_name=None,index_tag=None): 
    main(mafFile,outFile,seq_name,index_tag)

def new_main(file,seq_name=None,index_tag="ID"):
    maf = Maf(lines=file.readlines())
    bed = maf.to_bed(seq_name)
    for i in range(len(bed.entries)):
        id_string = "%s=%s" % (index_tag,i) if index_tag!=None else ""
        if bed.entries[i].name:
            bed.entries[i].name += ";"+id_string
        else:
            bed.entries[i].name = id_string
    return bed