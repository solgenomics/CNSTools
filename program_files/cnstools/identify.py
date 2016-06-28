import progress_tracker as pt
import subprocess
import json

from __consts import listOfModules
for module in listOfModules:
    if module!="identify":
        exec "from "+module+" import main as " + module
        

def main(data,out_folder):

    process = subprocess.Popen("mkdir -p "+out_folder, shell=True)
    process.wait()

    #maf_to_bed(roasted_maf_file) -> seq_bed w/maf_seq_ids
    header_print("Converting .maf file main sequence to .bed file locations.")
    data['seq_bed'] = out_folder+data['maf_file'].split("/")[-1].split(".maf")[0]+".bed"
    maf_to_bed(data['maf_file'],data['seq_bed'])

    #gff3_to_bed(main_gff3,"CDS") -> cds_bed
    header_print("Converting main .gff3 CDS ranges to .bed file locations.")
    data['main_gff3'] = data['seqs'][0]['gff3_file_name']
    data['cds_bed'] = out_folder+data['main_gff3'].split("/")[-1].split(".gff3")[0]+"_CDS.bed"
    gff3_to_bed(data['main_gff3'],['CDS'],data['cds_bed'])

    #$bedtools.subtract(seq_bed - cds_bed) -> cns_bed w/maf_seq_ids
    header_print("Subtracting CDS seqs from alignment.") 
    data['cns_bed'] = out_folder+data['seq_bed'].split("/")[-1].split(".bed")[0]+"_CNS.bed"
    cmd = "bedtools subtract -a %s -b %s > %s" % (data['seq_bed'],data['cds_bed'],data['cns_bed'])
    process = subprocess.Popen(cmd, shell=True)
    process.wait()

    #bed_maf_parse(cns_bed,roasted_maf_file) -> cns_maf w/cns_seq_ids
    header_print("Making a new .maf of only CNSs")
    data['cns_maf'] = out_folder+data['cns_bed'].split("/")[-1].split(".bed")[0]+".maf"
    bed_maf_parse(data['cns_bed'],data['maf_file'],data['cns_maf'])

    #maf_to_fasta(cns_maf) -> cns_fastas w/cns_seq_ids
    header_print("Making fasta files from .maf")
    data['cns_fasta_folder'] = out_folder+data['cns_maf'].split("/")[-1].split(".maf")[0]+"_fastas/"
    process = subprocess.Popen("mkdir "+data['cns_fasta_folder'], shell=True)
    process.wait()
    for seq in data['seqs']:
        seq['cns_fasta'] = data['cns_fasta_folder']+seq['maf_name']+".fasta"
    maf_to_fasta(data['cns_maf'],data['cns_fasta_folder'])

    #$makeblastdb(gff3_files) -> blasts_dbs
    header_print("Building blast databases for each genome.")
    data['blast_db_folder'] = out_folder+data['cns_maf'].split("/")[-1].split(".maf")[0]+"_blast_dbs/"
    cmd = "mkdir %s" % (data['blast_db_folder'])
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    for seq in data['seqs']:
        if seq['genome_fasta_name']=="": 
            break
        genome_file = seq['genome_fasta_name']
        seq['db_name'] = data['blast_db_folder']+seq['maf_name']+"_db"
        cmd = "makeblastdb -in %s -out %s -dbtype nucl -hash_index" % (genome_file,seq['db_name'])
        print 'builidng blast database for %s...' % (seq['maf_name'])
        process = subprocess.Popen(cmd, shell=True,stdout=subprocess.PIPE)
        process.wait()
        if process.stderr!=None: 
            for line in process.stderr.readlines(): print line
        for line in process.stdout.readlines(): print line

    #$blast(cns_fastas@blast_dbs) -> cns_blasts w/cns_seq_ids
    header_print("Blasting CNSs for each genome.")
    data['blast_results_folder'] = out_folder+data['cns_maf'].split("/")[-1].split(".maf")[0]+"_blast_results/"
    process = subprocess.Popen("mkdir "+data['blast_results_folder'], shell=True)
    process.wait()
    for seq in data['seqs']:
        if seq['genome_fasta_name']=="": 
            break
        print 'blasting CNSs in %s against %s...' % (seq['cns_fasta'],seq['db_name'])
        depth = len(data['blast_results_folder'].split("/"))
        up = "../"*depth
        out = data['blast_results_folder']+seq['maf_name']+".txt"
        cmd = "cd %s; blastn -db %s -query %s -out %s -outfmt 6; cd %s" % (data['blast_results_folder'],seq['db_name'].split("/")[-1],up+seq['cns_fasta'],up+out,up)
        process = subprocess.Popen(cmd, shell=True,stdout=subprocess.PIPE)
        process.wait()
        if process.stderr!=None: 
            for line in process.stderr.readlines(): print line
        for line in process.stdout.readlines(): print line

    #blast_to_bed(cns_blasts) -> cns_locs w/cns_seq_ids
    #gff3_to_bed(gff3_files,"gene") -> gene_beds
    #$bedtools.closest(cns_locs@gene_beds) -> cns_assoc_data w/distance & cns_seq_ids
    #!parse_assoc_info(cns_assoc_data) -> cns_assoc_info w/cns_seq_ids
    #OUTPUT:
    #    cns alignment from cns_maf
    #    location in each species from cns_locs
    #    closest genes and distance for each species from cns_assoc_info
    #    catagory from cns_assoc_info
    pass

def header_print(header):
    wall = "-"*(len(header)+2)
    print "%s\n %s \n%s" % (wall,header,wall)
    return

def run(argv):
    config = None
    with open(argv[1]) as intructionJSON:
        config = json.load(intructionJSON)
    out_folder = argv[2] if argv[2].endswith("/") else argv[2]+"/"
    main(config,out_folder)

if __name__ == '__main__':
    import sys
    run(sys.argv)