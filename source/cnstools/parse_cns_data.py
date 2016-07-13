import progress_tracker as pt
from filetypes import Maf, Cns, Bed13

def main(data,out_folder):

    maf = None
    cns_out = Cns()
    with open(data['cns_maf']) as maf_file:
        maf = Maf(lines = maf_file.readlines())


    cns_entry_by_ID = {}
    tracker = pt.Progress_tracker("Converting to .cns",len(maf.entries)).display(estimate=False,rate=1)
    for entry in maf.entries:
        ID = entry.a_line.split("cns_ID=")[1].split()[0]
        if not ID in cns_entry_by_ID: cns_entry_by_ID[ID] = cns_out.add_entry(ID)
        for sequence in entry.sequences:
            cns_entry_by_ID[ID].add_seq(sequence.src,None,None,None,None,None,None,None,None,sequence.text)
        tracker.step()
    tracker.display()
    del tracker


    assoc_files = [(seq['maf_name'],open(seq['cns_assoc'])) for seq in data['seqs'] if 'cns_assoc' in seq]
    tracker = pt.Progress_tracker("Parsing association files",len(assoc_files)).display(estimate=False,rate=1)
    for name,file in assoc_files:
        a_dat = Bed13(lines=file.readlines())
        for entry in a_dat.entries:
            info_list = entry.first.name.split('|')
            ID = info_list[2].split('cns_ID=')[1]
            genome = info_list[0]
            sequences = cns_entry_by_ID[ID].get_seqs(genome)
            if sequences[-1].closest_gene!=None:
                sequences.append(sequences[-1].duplicate())
            modSeq = sequences[-1]

            modSeq.loc_chrom = entry.first.chrom
            modSeq.closest_gene = entry.second.name
            modSeq.start = entry.first.chromStart
            modSeq.stop = entry.first.chromEnd
            modSeq.gene_start = entry.second.chromStart
            modSeq.gene_stop = entry.second.chromEnd

            modSeq.dist = entry.score
            if(abs(modSeq.dist)<=1000):
                if   modSeq.dist>0 or modSeq.start==modSeq.gene_stop: modSeq.type = "downstream"
                elif modSeq.dist<0 or modSeq.stop==modSeq.gene_start: modSeq.type = "upstream"
                elif modSeq.dist==0: modSeq.type = "intronic"
                else: modSeq.type = "???"
            else:
                modSeq.type = "intergenic"

        tracker.step()
    tracker.display()
    del tracker

    with open(out_folder+"results.cns","w") as out:
        out.write("\n".join(cns_out.get_lines()))

def file_run(data,out_folder):
    if not outFolder.endswith("/"):
        outFolder+="/"
    main(data,out_folder)

