from _filetypes import Cns, Fasta 
import os

def main(cns_file,sequence_types,out_folder):
    with open(cns_file) as file:
        cns = Cns(lines=file.readlines())
        fastas = cns.to_fasta(sequence_types)
        for key in fastas:
            with open("%s%s.fasta"%(out_folder,key),"w") as out:
                out.write("\n".join(fastas[key].get_lines()))

def file_run(cns_file,out_folder,*args): 
    sequence_types = list(*args)
    if not os.path.exists(out_folder): os.makedirs(out_folder)
    main(cns_file,sequence_types,out_folder)