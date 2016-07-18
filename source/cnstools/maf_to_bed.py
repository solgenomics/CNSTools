from _filetypes import Maf, Bed6

def file_run(mafFile,outFile,seq_name=None,index_tag=None):
    return Maf(file_name=mafFile).to_bed(seq_name,index_tag).save_file(outFile)