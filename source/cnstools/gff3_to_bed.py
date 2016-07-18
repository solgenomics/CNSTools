from _filetypes import Gff3, Bed6 

def file_run(gff3_file,out_file,*args): 
    type_list = list(args) if len(args)>0 else None

    Gff3(file_name=gff3_file) \
        .to_bed(type_list) \
        .save_file(out_file)