from _filetypes import Gff3, Bed6 

def main(gff3_file,sequence_types,out_file):
    with open(gff3_file) as file, open(out_file,"w") as out:
        gff3 = Gff3(lines=file.readlines())
        bed = gff3.to_bed(type_list=sequence_types)
        out.write("\n".join(bed.get_lines()))

def file_run(gff3_file,out_file,*args): 
    sequence_types = list(args) if len(args)>0 else None
    main(gff3_file,sequence_types,out_file)