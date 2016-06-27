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
    #   $bedtools.subtract(seq_beds - alinged_to_cds_bed) -> cns_beds w/maf_seq_ids
    cns_beds = []
    for i in range(len(seq_beds)):
        cns_beds.append(out_foler+seq_beds[i].split("/")[-1].split(".bed")[0]+"_masked.bed")
        process = subprocess.Popen("bedtools subtract -a "+seq_beds[i]+" -b "+cds_beds[i]+" > "+cns_beds[-1], shell=True)
        process.wait()
    #   !bed_maf_parse(cns_beds,roasted_maf_files) -> cns_maf w/cns_seq_ids
    cns_mafs = []
    for i in range(len(cns_beds)):
        cns_mafs.append(out_foler+cns_beds[i].split("/")[-1].split(".bed")[0]+"_cns.maf")
        bed_maf_parse(cns_beds[i],roasted_maf_files[i],cns_mafs[-1])
    #   maf_to_fasta(cns_maf) -> cns_fastas w/cns_seq_ids
    cns_fastas = []
    for i in range(len(cns_mafs)):
        cns_fastas.append(out_foler+cns_mafs[i].split("/")[-1].split(".maf")[0]+"_fastas/")
        process = subprocess.Popen("mkdir "+cns_fastas[-1], shell=True)
        process.wait()
        maf_to_fasta(cns_mafs[i],cns_fastas[-1])
    #   $makeblastdb(gff3_files) -> blasts_dbs
        #"makeblastdb -in file -out name -dbtype nucl -hash_index"
    #   $blast(cns_fastas@blast_dbs) -> cns_blasts w/cns_seq_ids
        for i in range(len(cns_fastas)):
            process = subprocess.Popen("find "+cns_fastas[i]+"*", shell=True,stdout=subprocess.PIPE)
            process.wait()
            files = [file.strip() for file in process.stdout.readlines()]
    #   blast_to_bed(cns_blasts) -> cns_locs w/cns_seq_ids
    #   gff3_to_bed(gff3_files,"gene") -> gene_beds
    #   $bedtools.closest(cns_locs@gene_beds) -> cns_assoc_data w/distance & cns_seq_ids
    #   !parse_assoc_info(cns_assoc_data) -> cns_assoc_info w/cns_seq_ids
    #   OUTPUT:
    #       cns alignment from cns_maf
    #       location in each species from cns_locs
    #       closest genes and distance for each species from cns_assoc_info
    #       catagory from cns_assoc_info
    pass

def run(argv):
    out_foler = argv[3] if argv[3].endswith("/") else argv[3]+"/"
    main([argv[1]],[argv[2]],out_foler)

if __name__ == '__main__':
    import sys
    run(sys.argv)