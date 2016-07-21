from _filetypes import Cns 
def main(cns_maf,cns_proxim_beds,index_tag="cns_maf_index"):
    cns = Cns()
    for i in range(len(cns_maf.entries)):
        cns.add_entry(i)
    for genome in cns_proxim_beds:
        for entry in cns_proxim_beds[genome].entries:
            try:
                kv = (nam for nam in entry.first.name.split(";") if nam.find(index_tag)!=-1).next()
            except StopIteration:
                raise ValueError("The provided index_tag was not found in a bed entry.")

            cns_index = int(kv.split('=')[1].strip())
            cns_entry = cns.entries[cns_index]

            loc_chrom = entry.first.chrom
            closest_gene = entry.second.name
            start = entry.first.chromStart
            stop = entry.first.chromEnd
            gene_start = entry.second.chromStart
            gene_stop = entry.second.chromEnd

            try:
                sequence_text = (sequence.text for sequence in cns_maf.entries[cns_index].sequences if sequence.src == loc_chrom).next()
            except StopIteration:
                raise ValueError("Missing sequence information in .maf file.")

            dist = entry.score
            cns_type = "???"
            if(abs(dist)<=1000):
                if   dist>0 or start==gene_stop: cns_type = "downstream"
                elif dist<0 or stop==gene_start: cns_type = "upstream"
                elif dist==0: cns_type = "intronic"
            else:
                cns_type = "intergenic"

            cns_entry.add_seq(genome,cns_type,dist,loc_chrom,closest_gene,start,stop,gene_start,gene_stop,sequence_text)
    return cns