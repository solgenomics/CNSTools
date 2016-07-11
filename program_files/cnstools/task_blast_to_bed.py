from filetypes import BlastF6, Bed6 

def main(blast_file,outFile):
    with open(blast_file) as f, open(outFile,"w") as out:
        blast = BlastF6(lines=f.readlines())
        bed = blast.to_bed()
        lines = bed.get_lines()
        out.write("\n".join(lines))

def file_run(blast_file,outFile):
    main(blast_file,outFile)