from .serial_filetypes import Serial_Filetype
from .bed_filetypes import Bed6
from .._utils import Progress_tracker

class Gff3_entry(object):
    """1-based"""
    def __init__(self,seqid,source,type,start,end,score,strand,phase,attributes):
        #super(Gff3_entry, self).__init__()
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = int(start)
        self.end = int(end)
        try:
            self.score = float(score)
        except:
            self.score = None
        self.strand = strand
        self.phase = phase
        self.attributes = attributes
    def get_line(self):
        return "\t".join([str(item) for item in (self.seqid,self.source,self.type,self.start,self.end,self.score,self.strand,self.phase,self.attributes)])

class Gff3(Serial_Filetype):
    """1-based"""
    Entry_class = Gff3_entry
    def add_lines(self,lines):
        tracker = Progress_tracker("Parsing .gff3",len(lines)).auto_display().start()
        for line in lines:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if (len(fields)==9):
                    self.entries.append(Gff3_entry(*fields))
            tracker.step()
        tracker.done()
    def get_lines(self):
        lines = []
        for entry in self.entries:
            lines.append(entry.get_line())
        return lines
    def to_bed(self,type_list=None,genome=None):
        new_bed = Bed6()
        entry_selection = None
        if(type_list):
            entry_selection = [entry for entry in self.entries if entry.type in type_list]
        else:
            entry_selection = self.entries
        tracker = Progress_tracker("Converting to .bed",len(entry_selection)).auto_display().start()
        for entry in entry_selection:
            if(entry.start<entry.end):
                chromStart,chromEnd = entry.start,entry.end
            else:
                chromStart,chromEnd = entry.end,entry.start
            id_with_type = entry.attributes+";seqType="+entry.type
            chrom = entry.seqid if not genome else genome+":"+entry.seqid
            new_bed.add_entry(chrom, chromStart-1, chromEnd, name=id_with_type, score=entry.score, strand=entry.strand)
            tracker.step()
        tracker.done()
        return new_bed

class Gff2(Gff3):
    def add_lines(self,lines):
        tracker = Progress_tracker("Parsing .gff3",len(lines)).auto_display().start()
        for line in lines:
            if not line.startswith('#'):
                fields = [i for i in line.strip().split('\t') if i!='']
                if (len(fields)==8):
                    fields = fields[:2]+["?"]+fields[2:] #unknown type!
                    self.entries.append(Gff3_entry(*fields))
            tracker.step()
        tracker.done()