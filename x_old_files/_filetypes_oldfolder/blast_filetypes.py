from .serial_filetypes import Serial_Filetype
from .bed_filetypes import Bed6
from .._utils import Progress_tracker

class BlastF6_entry(object):
    """1-based"""
    def __init__(self,query,target,identity,length,mismatches,gapOpens,queryStart,queryEnd,targetStart,targetEnd,eVal,bitScore):
        #super(BlastF6_entry, self).__init__()
        self.query = query
        self.target = target
        self.identity = float(identity)
        self.length = int(length)
        self.mismatches = int(mismatches)
        self.gapOpens = int(gapOpens)
        self.queryStart = int(queryStart)
        self.queryEnd = int(queryEnd)
        self.targetStart = int(targetStart)
        self.targetEnd = int(targetEnd)
        self.eVal = float(eVal)
        self.bitScore = float(bitScore)
    def get_line(self):
        return "\t".join([str(item) for item in (self.query,self.target,self.identity,self.length,self.mismatches,self.gapOpens,self.queryStart,self.queryEnd,self.targetStart,self.targetEnd,self.eVal,self.bitScore)])

class BlastF6(Serial_Filetype):
    """1-based"""
    Entry_class = BlastF6_entry
    def add_lines(self,lines):
        tracker = Progress_tracker("Parsing blast output",len(lines)).auto_display().start()
        for line in lines:
            fields = line.strip().split('\t')
            if (len(fields)==12):
                self.entries.append(BlastF6_entry(*fields))
            tracker.step()
        tracker.done()
    def get_lines(self):
        lines = []
        for entry in self.entries:
            lines.append(entry.get_line())
        return lines
    def to_bed(self):
        new_bed = Bed6()
        tracker = Progress_tracker("Converting to .bed",len(self.entries)).auto_display().start()
        for entry in self.entries:
            strand = None
            if entry.targetStart < entry.targetEnd:
                start, end = entry.targetStart, entry.targetEnd
                strand = '+'
            else:
                start, end = entry.targetEnd, entry.targetStart
                strand = '-'
            #convert to 0-based!
            new_bed.add_entry(entry.target, start-1, end, entry.query, entry.bitScore, strand)
            tracker.step()
        tracker.done()
        return new_bed