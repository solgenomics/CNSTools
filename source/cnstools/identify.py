import json
import os
from copy import deepcopy
from _filetypes import Maf, Gff3, Bed6
from build_cns_file import main as build_cns_file
from _utils import create_path, JSON_saver, safe_print, header_print, Progress_tracker

json_file_format = {
    "all_seq_maf":"PATH",
    "all_conserved_bed":"PATH",
    "ref_genome":"PREFIX",
    "genomes":{
        "PREFIX":{
            "annot_gff3":"PATH"
        }
    }
}

def main(json_data,work_folder,num_threads,overwrite=False):
    datasaver = JSON_saver(create_path(work_folder,"record","json",overwrite))
    data = deepcopy(json_data)
    datasaver.save(data)

    info = "Convert aligned sequences to .bed:"
    header_print(info)
    data['all_seq_bed'] = create_path(work_folder,"all_seq","bed",overwrite)
    Maf(file_name=data['all_seq_maf']) \
        .to_bed(seq_name=None,index_tag="all_maf_index") \
        .save_file(data['all_seq_bed'])
    datasaver.save(data)

    info = "Intersect aligned regions with conserved regions:"
    header_print(info)
    data['conserved_bed'] = create_path(work_folder,"conserved","bed",overwrite)
    cmd = "bedtools intersect -a %s -b %s > %s" % (data['all_seq_bed'],data['all_conserved_bed'],data['conserved_bed'])
    tracker = Progress_tracker("Running bedtools intersect",1).estimate(False).display()
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    tracker.done()
    datasaver.save(data)

    info = "Convert coding regions to .bed:"
    header_print(info)
    data['coding_bed'] = create_path(work_folder,"coding","bed",overwrite)
    Gff3(file_name=data['genomes'][data['ref_genome']]['annot_gff3']) \
        .to_bed(['CDS']) \
        .save_file(data['coding_bed'])
    datasaver.save(data)

    info = "Subtract coding regions from conserved aligned regions:"
    header_print(info)
    data['cns_bed'] = create_path(work_folder,"cns","bed",overwrite)
    cmd = "bedtools subtract -a %s -b %s > %s" % (data['conserved_bed'],data['coding_bed'],data['cns_bed'])
    tracker = Progress_tracker("Running bedtools subtract",1).estimate(False).display()
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    tracker.done()
    datasaver.save(data)

    info = "Slice multi-alignment file based on identified conserved non-coding regions:"
    header_print(info)
    data['cns_maf'] = create_path(work_folder,"cns","maf",overwrite)
    cns_bed = Bed6(file_name=data['cns_bed'])
    Maf(file_name=data['all_seq_maf']) \
        .slice_with_bed(cns_bed,data['ref_genome'],"all_maf_index",max_N_ratio=0.5,max_gap_ratio=0.5,min_len=15) \
        .save_file(data['cns_maf'])
    del cns_bed
    datasaver.save(data)

    info = "Convert per-genome CNS regions to .bed:"
    header_print(info)
    data['genome_cns_beds_folder'] = create_path(work_folder+"genome_cns_beds",overwrite)
    cns_maf = Maf(file_name=data['cns_maf'])
    for genome in data['genomes']:
        data['genomes'][genome]['cns_bed'] = create_path(data['genome_cns_beds_folder'],genome+"_cns","bed",overwrite)
        cns_maf.to_bed(seq_name=genome,index_tag="cns_maf_index") \
            .save_file(data['genomes'][genome]['cns_bed'])
    del cns_maf
    datasaver.save(data)

    info = "Convert per-genome gene regions to .bed:"
    header_print(info)
    data['genome_annot_beds_folder'] = create_path(work_folder+"genome_annot_beds",overwrite)
    for genome in data['genomes']:
        data['genomes'][genome]['annot_bed'] = create_path(data['genome_annot_beds_folder'],genome+"_annot","bed",overwrite)
        Gff3(file_name=data['genomes'][genome]['annot_gff3']) \
            .to_bed(type_list=['gene'],genome=genome) \
            .save_file(data['genomes'][genome]['annot_bed'])
    datasaver.save(data)

    info = "Find closest gene for each CNS region:"
    header_print(info)
    data['gene_proximity_beds_folder'] = create_path(work_folder+"gene_proximity_beds",overwrite)
    for genome in data['genomes']:
        data['genomes'][genome]['gene_proximity_bed'] = \
            create_path(data['gene_proximity_beds_folder'],genome+"_proxim","bed",overwrite)
        cmd = "bedtools closest -D a -a %s -b %s > %s" % \
            (data['genomes'][genome]['cns_bed'],
             data['genomes'][genome]['annot_bed'],
             data['genomes'][genome]['gene_proximity_bed'])
        process = subprocess.Popen(cmd, shell=True)
        process.wait()
    datasaver.save(data)

    info = "Process proximity and maf files into .cns file:"
    header_print(info)
    data['results'] = create_path(work_folder,"identified_CNSs","cns",overwrite)
    cns_proxim_beds = {genome:Bed6(data['genomes'][genome]['gene_proximity_bed']) for genome in data['genomes']}
    Maf(file_name=data['cns_maf'])\
        .cns_from_proxim_beds(cns_proxim_beds,"cns_maf_index")
        .save_file(data['results'])
    datasaver.save(data)
    
def file_run(json_file,work_folder,num_threads,overwrite=False):
    config = None
    with open(json_file) as intructionJSON:
        config = json.load(intructionJSON)
    work_folder = create_path(work_folder,overwrite=overwrite)
    num_threads = int(num_threads)
    main(input_data,work_folder,num_threads,overwrite)
    
        
