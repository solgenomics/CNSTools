import progress_tracker as pt
import subprocess

from __consts import listOfModules
for module in listOfModules:
    if module!="identify":
        exec "from "+module+" import main as " + module
        

def main(roasted_maf_files,gff3_files,out_foler):
    #   maf_to_bed(roasted_maf_files) -> seq_beds w/maf_seq_ids
    seq_beds = []
    for file in roasted_maf_files:
        seq_beds.append(out_foler+file.split("/")[-1].split(".maf")[0]+".bed")
        maf_to_bed(file,seq_beds[-1])
    #   gff3_to_bed(gff3_files,"CDS") -> cds_beds w/ alinged_to_cds_bed file
    cds_beds = []
    for file in gff3_files:
        cds_beds.append(out_foler+file.split("/")[-1].split(".gff3")[0]+"_CDSs.bed")
        gff3_to_bed(file,['CDS'],cds_beds[-1])
    #   $bedtools.subtract(seq_beds - alinged_to_cds_bed) -> ncs_beds w/maf_seq_ids
    ncs_beds = []
    for i in range(len(seq_beds)):
        ncs_beds.append(out_foler+file.split("/")[-1].split(".bed")[0]+"_masked.bed")
        subprocess.Popen("bedtools subtract -a "+seq_beds[i]+" -b "+cds_beds[i]+" > "+ncs_beds[-1], shell=True)
    #   !bed_maf_parse(ncs_beds,roasted_maf_files) -> ncs_maf w/ncs_seq_ids
    #   maf_to_fasta(ncs_maf) -> ncs_fastas w/ncs_seq_ids
    #   $makeblastdb(gff3_files) -> blasts_dbs
    #   $blast(ncs_fastas@blast_dbs) -> ncs_blasts w/ncs_seq_ids
    #   blast_to_bed(ncs_blasts) -> ncs_locs w/ncs_seq_ids
    #   gff3_to_bed(gff3_files,"gene") -> gene_beds
    #   $bedtools.closest(ncs_locs@gene_beds) -> ncs_assoc_data w/distance & ncs_seq_ids
    #   !parse_assoc_info(ncs_assoc_data) -> ncs_assoc_info w/ncs_seq_ids
    #   OUTPUT:
    #       NCS alignment from ncs_maf
    #       location in each species from ncs_locs
    #       closest genes and distance for each species from ncs_assoc_info
    #       catagory from ncs_assoc_info
    pass

def run(argv):
    out_foler = argv[3] if argv[3].endswith("/") else argv[3]+"/"
    main([argv[1]],[argv[2]],out_foler)

if __name__ == '__main__':
    import sys
    run(sys.argv)