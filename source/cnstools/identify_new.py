    #maf_to_bed(roasted_maf_file) -> seq_bed w/maf_seq_ids

    #gff3_to_bed(main_gff3,"CDS") -> cds_bed

    #$bedtools.subtract(seq_bed - cds_bed) -> cns_bed w/maf_seq_ids

    # bed_maf_parse(cns_bed,roasted_maf_file) -> cns_maf w/cns_seq_ids

    # maf_to_fasta(cns_maf) -> cns_fastas w/cns_seq_ids

    #$makeblastdb(gff3_files) -> blasts_dbs

    #$blast(cns_fastas@blast_dbs) -> cns_blasts w/cns_seq_ids

    #blast_to_bed(cns_blasts) -> cns_locs w/cns_seq_ids

    #gff3_to_bed(gff3_files,"gene") -> gene_beds

    #$bedtools.closest(cns_locs@gene_beds) -> cns_assoc_data w/distance & cns_seq_ids

    #!parse_cns_data(data,out) -> cns_assoc_info w/cns_seq_ids
from _filetypes import Maf,Gff3,Fasta
from copy import deepcopy
def main(input_data,num_threads):
    data = 

def file_run(json_file,out_folder,num_threads):

    config = None
    with open(json_file) as intructionJSON:
        config = json.load(intructionJSON)
    if not out_folder.endswith("/"):
        out_folder+="/"
    num_threads = int(num_threads)
    main(input_data,num_threads)


