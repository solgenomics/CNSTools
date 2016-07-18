from .serial_filetypes import Serial_Filetype
from .._utils import Progress_tracker

class Cns_sequence(object):
    def __init__(self,genome,type,dist,loc_chrom,closest_gene,start,stop,gene_start,gene_stop,sequence,cns_ID=None):
        self.cns_ID = cns_ID
        self.genome = genome
        self.type = type
        self.dist = int(dist) if dist else None
        self.loc_chrom = loc_chrom
        self.closest_gene = closest_gene
        self.start = int(start) if start else None
        self.stop = int(stop) if stop else None
        self.gene_start = int(gene_start) if gene_start else None
        self.gene_stop = int(gene_stop) if gene_stop else None
        self.sequence = sequence
    def get_line(self):
        strs = (str(i) if i!=None else '.' for i in (self.cns_ID,self.genome,self.type,self.dist,self.loc_chrom,self.closest_gene,self.start,self.stop,self.gene_start,self.gene_stop,self.sequence))
        return "\t".join(strs)
    def duplicate(self):
        return Cns_sequence(self.genome,self.type,self.dist,self.loc_chrom,self.closest_gene,self.start,self.stop,self.gene_start,self.gene_stop,self.sequence,cns_ID=self.cns_ID)

class Cns_entry(object):
    def __init__(self,cns_ID):
        self.cns_ID = cns_ID
        self.sequences = {}
    def add_seq(self,genome,*args):
        if not genome in self.sequences: self.sequences[genome] = []
        self.sequences[genome].append(Cns_sequence(genome,*args,cns_ID=self.cns_ID))
        return self.sequences[genome][-1]
    def get_seqs(self,genome=None):
        if genome==None: return [seq for key in self.sequences for seq in self.sequences[key]]
        elif not genome in self.sequences: return None
        else: return self.sequences[genome][:]
    def get_lines(self):
        seqs = []
        for key in self.sequences:
            seqs+=self.sequences[key]
        return [seq.get_line() for seq in seqs]


class Cns(Serial_Filetype):
    Entry_class = Cns_entry
    def add_lines(self,lines):
        ID=None
        tracker = Progress_tracker("Parsing .cns",len(lines)).auto_display().start()
        for line in lines:
            list = [item if item!='.' else None for item in line.split('\t')]
            if list[0]!=ID:
                ID = list[0]
                self.entries.append(Cns_entry(ID))
            self.entries[-1].add_seq(*(list[1:]))
            tracker.step()
        tracker.done()
    def get_lines(self):
        return [line for entry in self.entries for line in entry.get_lines()]
    def to_fasta(self,sequences=False):
        """converts the .cns data to multiple .fasta data classes and returns a dict with sequence origin names as keys"""
        fastas = {}
        tracker = Progress_tracker("Converting to .fasta files",len(self.entries)).auto_display().start()
        for entry in self.entries:
            for seq in [seq for key in entry.sequences for seq in entry.sequences[key]]:
                if seq.genome not in fastas: fastas[seq.genome] = Fasta()
                description = "|".join((str(a) for a in (seq.cns_ID,seq.type,seq.loc_chrom,seq.start,seq.stop)))
                fastas[seq.genome].add_entry(description,seq.sequence.replace("-", ""))
            tracker.step()
        tracker.done()
        return fastas