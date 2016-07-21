from .serial_filetypes import Serial_Filetype
from .bed_filetypes import Bed6
from .cns_filetypes import Cns
from .._utils import Progress_tracker

class Maf_sequence(object):
    def __init__(self,src,start,size,strand,srcSize,text,metadata=None):
        self.src = src
        self.start = int(start)
        self.size = int(size)
        self.strand = strand
        self.srcSize = int(srcSize)
        self.text = text
        self.metadata = metadata
    def get_lines(self):
        lines = ['##'+metadata] if metadata else []
        lines.append('\t'.join([str(item) for item in [self.src,self.start,self.size,self.strand,self.srcSize,self.text]]))
        return lines
        
class Maf_entry(object):
    def __init__(self,paragraph=None):
        self.a_meta = None
        self.a_line = None
        self.sequences = []
        if paragraph:
            rec_metadata = None
            for line in paragraph:
                if line.startswith('##'):
                    rec_metadata = line.split("##")[1]
                elif line.startswith('a'):
                    self.a_line = line
                    if rec_metadata:
                        self.a_meta = rec_metadata
                        rec_metadata = None
                elif line.startswith('s'):
                    vals = line.split()[1:]
                    if rec_metadata:
                        vals.append(rec_metadata)
                        rec_metadata = None
                    self.sequences.append(Maf_sequence(*vals))
    def get_lines(self):
        lines = []
        if self.a_meta: lines.append('##'+self.a_meta)
        if self.a_line: lines.append(self.a_line)
        for sequence in self.sequences:
            lines+=sequence.get_lines()
        return lines

class Maf(Serial_Filetype):
    """0-based"""
    Entry_class = Maf_entry

    def add_lines(self,lines):
        if not hasattr(self, 'headerLines'): self.headerLines = []
        paragraph = []
        tracker = Progress_tracker("Parsing .maf",len(lines)).auto_display().start()
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("#") and not stripped.startswith("##--"):
                self.headerLines.append(stripped)
            elif stripped=="":
                if len(paragraph)>1:
                    self.entries.append(Maf_entry(paragraph))
                paragraph = []
            else:
                paragraph.append(stripped)
            tracker.step()
        tracker.done()

    def get_lines(self):
        lines = []
        for entry in self.entries:
            lines+=entry.get_lines()
            lines.append("")
        return lines

    def to_bed(self,seq_name=None,index_tag="maf_index"):
        new_bed = Bed6()
        tracker = Progress_tracker("Converting to .bed",len(self.entries)).auto_display().start()
        if not seq_name: 
            seq_name=self.entries[0].sequences[0].src.split(":")[0]
        index = 0
        for entry in self.entries:
            seq_to_convert = (seq for seq in entry.sequences if seq.src.split(":")[0]==seq_name)
            for sequence in seqs_to_convert:
                id_string = "%s=%s" % (index_tag,index) if index_tag!=None else None
                if sequence.metadata:
                    id_string = (sequence.metadata+";"+id_string) if id_string else (sequence.metadata)
                new_bed.add_entry(sequence.src, sequence.start, sequence.start+sequence.size, name=id_string, strand=sequence.strand)
            index+=1
            tracker.step()
        tracker.done()
        return new_bed

    def slice_with_bed(self,bed,ref_genome,index_tag,max_N_ratio=0.5,max_gap_ratio=0.5,min_len=15):
        new_maf = Maf()
        new_maf.headerlines = self.headerLines[:]
        for bed_entry in bed.entries:
            try:
                keyVal = (kv for kv in bed_entry.name.split(";") if kv.find(index_tag)!=-1).next()
            except StopIteration:
                raise ValueError("The provided index_tag was not found in a bed entry.")

            parent_index = int(keyVal.split(index_tag+"=")[1])
            original_maf_entry = self.entries[parent_index]
            new_maf_entry = Maf_entry()
            new_maf_entry.a_line = original_maf_entry.a_line
            new_maf_entry.a_meta = original_maf_entry.a_meta

            try:
                mainseq = (seq for seq in original_maf_entry.sequences if seq.chrom.split(":")[0].strip()==ref_genome).next()
            except StopIteration:
                raise ValueError("The provided ref_genome was not found in an alignment.")

            front_offset = bed.chromStart - ref_seq.start
            back_offset = ref_seq.stop - bed.chromEnd
            front_cut_num = Maf._gap_cut_loc(ref_seq.text,front_offset)
            back_cut_num = Maf._gap_cut_loc(reversed(ref_seq.text),back_offset)

            for seq in original_maf_entry.sequences:
                front_removed_seq = original_maf_entry.text[:front_cut_num]
                new_seq = (original_maf_entry.text[front_cut_num:-back_cut_num]) if back_cut_num!=0 else (original_maf_entry.text[front_cut_num:])
                new_len = Maf._no_gap_len(new_seq)
                new_start = original_maf_entry.start + Maf._no_gap_len(front_removed_seq)
                if 1-(new_len/len(new_seq)) <= max_gap_ratio:
                    if 1-Maf._no_gap_len(new_seq.replace('N',''))/new_len <= max_N_ratio:    
                        new_seq = Maf_sequence(
                            original_maf_entry.src,
                            new_start,
                            new_len,
                            original_maf_entry.strand,
                            original_maf_entry.srcSize,
                            new_seq,
                            original_maf_entry.metadata)
                        new_maf_entry.sequences.append(new_seq)
            if len(new_maf_entry.sequences) > 1 and any((seq.size>=min_len for seq in new_maf_entry.sequences)):
                new_maf.entries.append(new_maf_entry)
        return new_maf

    def cns_from_proxim_beds(self,cns_proxim_beds,index_tag="maf_index"):
        cns = Cns()
        for i in range(len(self.entries)):
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
                    sequence_text = (sequence.text for sequence in self.entries[cns_index].sequences if sequence.src == loc_chrom).next()
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

    @staticmethod
    def _gap_cut_loc(seq,cutNum):
        anonymized = "".join(["x" if i!="-" else i for i in seq]) #use x to represent valid chars
        eliminated = anonymized.replace("x","o",cutNum) #change the first n x's to o's
        return eliminated.find("x") # finds the index of the next x, which is how many chars to cut off
    @staticmethod
    def _no_gap_len(seq): return len(seq[:].replace("-",""))
