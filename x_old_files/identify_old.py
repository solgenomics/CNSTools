from _utils import Progress_tracker,header_print
import subprocess
import json
import os

from __init__ import names as listOfModules
for name in listOfModules:
    if name!="identify":
        exec "from "+name+" import main as " + name
        

def main(data,out_folder,num_threads):
    num_availible_threads = num_threads-1

    process = subprocess.Popen("mkdir -p "+out_folder, shell=True)
    process.wait()

    #print(json.dumps(data, sort_keys=True, indent=4, separators=(',', ': ')))

    #maf_to_bed(roasted_maf_file) -> seq_bed w/maf_seq_ids
    if not 'seq_bed' in data:
        header_print("Converting .maf file main sequence to .bed file ranges.")
        data['seq_bed'] = out_folder+data['maf_file'].split("/")[-1].split(".maf")[0]+".bed"
        maf_to_bed(data['maf_file'],data['seq_bed'])

    #gff3_to_bed(main_gff3,"CDS") -> cds_bed
    header_print("Converting main .gff3 CDS ranges to .bed file ranges.")
    if not 'main_gff3' in data:
        data['main_gff3'] = data['seqs'][0]['gff3_file_name']
        data['cds_bed'] = out_folder+data['main_gff3'].split("/")[-1].split(".gff3")[0]+"_CDS.bed"
        gff3_to_bed(data['main_gff3'],['CDS'],data['cds_bed'])

    #$bedtools.subtract(seq_bed - cds_bed) -> cns_bed w/maf_seq_ids
    header_print("Subtracting CDS seqs from alignment.")
    if not 'cns_bed' in data: 
        data['cns_bed'] = out_folder+data['seq_bed'].split("/")[-1].split(".bed")[0]+"_CNS.bed"
        cmd = "bedtools subtract -a %s -b %s > %s" % (data['seq_bed'],data['cds_bed'],data['cns_bed'])
        tracker = Progress_tracker("Running bedtools",1).estimate(False).display()
        process = subprocess.Popen(cmd, shell=True)
        process.wait()
        tracker.done()

    # bed_maf_parse(cns_bed,roasted_maf_file) -> cns_maf w/cns_seq_ids
    header_print("Making a new .maf of only CNSs")
    if not 'cns_maf' in data:
        data['cns_maf'] = out_folder+data['cns_bed'].split("/")[-1].split(".bed")[0]+".maf"
        bed_maf_parse(data['cns_bed'],data['maf_file'],data['cns_maf'])

    # maf_to_fasta(cns_maf) -> cns_fastas w/cns_seq_ids
    header_print("Making fasta files from .maf")
    if not 'cns_fasta_folder' in data:
        data['cns_fasta_folder'] = out_folder+data['cns_maf'].split("/")[-1].split(".maf")[0]+"_fastas/"
        process = subprocess.Popen("mkdir "+data['cns_fasta_folder'], shell=True)
        process.wait()
        for seq in data['seqs']:
            seq['cns_fasta'] = data['cns_fasta_folder']+seq['maf_name']+".fasta"
        maf_to_fasta(data['cns_maf'],data['cns_fasta_folder'])

    #$makeblastdb(gff3_files) -> blasts_dbs
    header_print("Building BLAST databases for each genome.")
    if not 'blast_db_folder' in data:
        data['blast_db_folder'] = out_folder+data['cns_maf'].split("/")[-1].split(".maf")[0]+"_blast_dbs/"
        cmd = "mkdir %s" % (data['blast_db_folder'])
        process = subprocess.Popen(cmd, shell=True)
        process.wait()
        tracker = Progress_tracker("Building BLAST databases",len(data['seqs'])).estimate(False).display()
        for seq in data['seqs']:
            genome_file = seq['genome_fasta_name']
            seq['db_name'] = data['blast_db_folder']+seq['maf_name']+"_db"
            cmd = "makeblastdb -in %s -out %s -dbtype nucl -hash_index" % (genome_file,seq['db_name'])
            tracker.status('building %s database' % (seq['maf_name']))
            process = subprocess.Popen(cmd, shell=True,stdout=subprocess.PIPE)
            process.wait()
            if process.stderr!=None: 
                for line in process.stderr.readlines(): print line+"\n\n"
            tracker.step().display()
        tracker.status(None).done()

    #$blast(cns_fastas@blast_dbs) -> cns_blasts w/cns_seq_ids
    header_print("BLASTing CNSs for each genome.")
    if not 'blast_results_folder' in data:
        data['blast_results_folder'] = out_folder+data['cns_maf'].split("/")[-1].split(".maf")[0]+"_blast_results/"
        process = subprocess.Popen("mkdir "+data['blast_results_folder'], shell=True)
        process.wait()
        tracker = Progress_tracker("Running BLAST",len(data['seqs'])).estimate(False).display()
        for seq in data['seqs']:
            tracker.status('running BLAST for %s' % (seq['maf_name']))
            depth = len(data['blast_db_folder'].split("/"))-1
            up = "../"*depth
            seq['cns_blast'] = data['blast_results_folder']+seq['maf_name']+".txt"
            cmd = "cd %s; blastn -db %s -query %s -out %s -outfmt 6 -word_size 14 -dust no -ungapped -perc_identity 100 -qcov_hsp_perc 100 -num_threads %s -culling_limit 1 -penalty -100; cd %s" \
            % (data['blast_db_folder'],seq['db_name'].split("/")[-1],up+seq['cns_fasta'],up+seq['cns_blast'],num_availible_threads,up)
            process = subprocess.Popen(cmd, shell=True,stdout=subprocess.PIPE)
            process.wait()
            if process.stderr!=None: 
                for line in process.stderr.readlines(): print line
            tracker.step().display()
        tracker.status(None).done()

    #blast_to_bed(cns_blasts) -> cns_locs w/cns_seq_ids
    header_print("Filtering and converting Blast results to .bed files for each genome.")
    if not 'cns_locs_folder' in data:
        data['cns_locs_folder'] = out_folder+data['cns_maf'].split("/")[-1].split(".maf")[0]+"_cns_locs/"
        process = subprocess.Popen("mkdir "+data['cns_locs_folder'], shell=True)
        process.wait()
        for seq in data['seqs']:
            print "%s:"%seq['maf_name']
            seq['cns_locs'] = data['cns_locs_folder']+seq['maf_name']+".bed"
            blast_to_bed(seq['cns_blast'],seq['cns_locs'])

    #gff3_to_bed(gff3_files,"gene") -> gene_beds
    header_print("Converting genomes' .gff gene ranges to .bed file locations.")
    if not 'gene_bed_folder' in data:
        data['gene_bed_folder'] = out_folder+data['cns_maf'].split("/")[-1].split(".maf")[0]+"_gene_beds/"
        process = subprocess.Popen("mkdir "+data['gene_bed_folder'], shell=True)
        process.wait()
        for seq in data['seqs']:
            if seq['gff3_file_name']!="":
                print "%s:"%seq['maf_name']
                seq['gene_bed'] = data['gene_bed_folder']+seq['maf_name'].split("/")[-1].split(".gff")[0]+"_genes.bed"
                gff3_to_bed(seq['gff3_file_name'],['gene'],seq['gene_bed'])

    #$bedtools.closest(cns_locs@gene_beds) -> cns_assoc_data w/distance & cns_seq_ids
    header_print("Checking CNS locations against gene locations.")
    if not 'cns_assoc_folder' in data:
        data['cns_assoc_folder'] = out_folder+data['cns_maf'].split("/")[-1].split(".maf")[0]+"_cns_assoc/"
        process = subprocess.Popen("mkdir "+data['cns_assoc_folder'], shell=True)
        process.wait()
        tracker = Progress_tracker("Finding associations",len(data['seqs'])).estimate(False).display()
        for seq in data['seqs']:
            tracker.status(seq['maf_name'])
            if seq['gff3_file_name']!="":
                seq['cns_assoc'] = data['cns_assoc_folder']+seq['maf_name'].split("/")[-1].split(".gff")[0]+"_cns_assoc.bed"
                cmd = "bedtools closest -D a -a %s -b %s > %s" % (seq['cns_locs'],seq['gene_bed'],seq['cns_assoc'])
                process = subprocess.Popen(cmd, shell=True)
                process.wait()
            tracker.step().display()
        tracker.status(None).done()

    #!parse_cns_data(data,out) -> cns_assoc_info w/cns_seq_ids
    header_print("Parsing association data.")
    if not 'final_results_folder' in data:
        data['final_results_folder'] = out_folder+data['cns_maf'].split("/")[-1].split(".maf")[0]+"_results/"
        process = subprocess.Popen("mkdir "+data['final_results_folder'], shell=True)
        process.wait()
        parse_cns_data(data,data['final_results_folder'])

    with open(out_folder+"data.json","w") as out:
        out.write(json.dumps(data, sort_keys=True, indent=4, separators=(',', ': ')))

    #make_figures(data['final_results_folder'])
    #OUTPUT:
    #    cns alignment from cns_maf
    #    location in each species from cns_locs
    #    closest genes and distance for each species from cns_assoc_info
    #    catagory from cns_assoc_info
    header_print("Finished! Results in:%s"%data['final_results_folder'])

def file_run(json_file,out_folder,num_threads_in):
    config = None
    with open(json_file) as intructionJSON:
        config = json.load(intructionJSON)
    if not out_folder.endswith("/"):
        out_folder+="/"
    num_threads = int(num_threads_in)
    main(config,out_folder,num_threads)