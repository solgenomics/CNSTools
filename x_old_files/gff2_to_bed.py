from _filetypes import Gff2, Bed6 

def main(gff2_file,sequence_types,out_file):
    with open(gff2_file) as file, open(out_file,"w") as out:
        gff2 = Gff2(lines=file.readlines())
        bed = gff2.to_bed(type_list=sequence_types)
        out.write("\n".join(bed.get_lines()))

def file_run(gff2_file,out_file,*args): 
    sequence_types = list(args) if len(args)>0 else None
    main(gff2_file,sequence_types,out_file)