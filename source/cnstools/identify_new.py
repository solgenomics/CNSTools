# *all_seq_maf

# *all_conserved_bed

# *ref_genome

# *genome_annot_gff3s

json_file_format = {
    "all_seq_maf":"PATH"
    "all_conserved_bed":"PATH",
    "ref_genome":"PREFIX",
    "genomes":{
        "PREFIX":{
            "annot_gff3":"PATH"
        }
    }
}

import json
import os
from copy import deepcopy
from _filetypes import Maf,Gff3,Fasta
import json

def main(json_data,work_folder,num_threads,overwrite=False):
    datasaver = JSON_saver(create_path(work_folder,"record","json",overwrite))
    data = deepcopy(json_data)
    datasaver.set(data).save()

# all_seq_maf -> all_seq_bed.all_seq_maf_indexes(all_maf)
    data['all_seq_bed'] = create_path(work_folder,"all_seq_bed","bed",overwrite)
    maf_to_bed(Maf(file_name=data['all_seq_maf']),index_tag="all_maf").save_file(data['all_seq_bed'])
    datasaver.set(data).save()

# intersect(all_conserved_bed,all_seq_bed.all_seq_maf_indexes) -> index_conserved_bed.all_seq_maf_indexes
    data['index_conserved_bed'] = create_path(work_folder,"index_conserved_bed","bed",overwrite)
    cmd = "bedtools intersect -a %s -b %s > %s" % (data['all_seq_bed'],data['all_conserved_bed'],data['index_conserved_bed'])
    tracker = Progress_tracker("Running bedtools",1).estimate(False).display()
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    tracker.done()
    datasaver.set(data).save()


# genome_annot_gff3s[ref_genome]('CDS') -> coding_bed

# subtract(index_conserved_bed.all_seq_maf_indexes - coding_bed) -> cns_bed.all_seq_maf_indexes

# subset_parse_to_maf(cns_bed.all_seq_maf_indexes from all_seq_maf) -> cns_maf

# split_to_bed(cns_maf) -> per_genome_cns_beds.cns_maf_indexes

# genome_annot_gff3s('gene') -> per_genome_gene_beds

# closest(per_genome_gene_beds to per_genome_cns_beds.cns_maf_indexes) -> cns_assoc_raw_bed

# make_cns_file(cns_assoc_raw_bed,cns_maf,genome_annot_gff3s)




def file_run(json_file,work_folder,num_threads,overwrite=False):
    config = None
    with open(json_file) as intructionJSON:
        config = json.load(intructionJSON)
    if not work_folder.endswith("/"):
        work_folder+="/"
    num_threads = int(num_threads)
    main(input_data,num_threads)

def create_path(path,name=None,extension=None,overwrite=True):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

    if not path.endswith("/"):
        path+="/"
    if name:
        path+=name
    if extension:
        if not extension.startswith("."):
            path+="."
        path+=extension

    if (not overwrite) and os.path.exists(path):
        raise ValueError("You provided an existing path but asked not to overwrite!")

    return path

class JSON_saver(object):
    """docstring for JSON_saver"""
    def __init__(self, path):
        self.data = None
        self.path = path
    def set(self,data):
        self.data = deepcopy(data)
        return self
    def access(self):
        return deepcopy(self.data)
    def save(self):
        with open(self.path,'w') as out:
            json.dump(self.data,out,sort_keys=True,indent=4)
        return self
    

        

        
