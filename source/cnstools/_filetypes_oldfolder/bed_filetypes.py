from .serial_filetypes import Serial_Filetype
from __init__.cnstools._utils import Progress_tracker
class Bed6_entry(object):
    """0-based"""
    def __init__(self, chrom, chromStart, chromEnd, name=None, score=None, strand=None):
        #super(Bed6_entry, self).__init__()
        self.chrom = chrom
        self.chromStart = int(chromStart)
        self.chromEnd = int(chromEnd)
        self.name = name
        try:
            self.score = float(score)
        except:
            self.score = None
        self.strand = strand
    def get_line(self):
        strs = (str(i) if i!=None else '.' for i in (self.chrom, self.chromStart, self.chromEnd, self.name, self.score, self.strand))
        return "\t".join(strs)

class Bed6(Serial_Filetype):
    """0-based"""
    Entry_class = Bed6_entry
    def add_lines(self,lines):
        tracker = Progress_tracker("Parsing 6 column .bed",len(lines)).auto_display().start()
        for line in lines:
            fields = line.strip().split('\t')
            if len(fields)>1:
                if len(fields)<6:
                    fields.append([None]*(6-len(fields)))
                fields[:] = [item if item!='.' else None for item in fields]
                self.entries.append(Bed6_entry(*fields))
            tracker.step()
        tracker.done()
    def get_lines(self):
        lines = []
        for entry in self.entries:
            lines.append(entry.get_line())
        return lines

class Bed13_entry(object):
    def __init__(self,*fields):
        if len(fields)!=13:
            raise Exception("not 13 fields")
        self.first = Bed6_entry(*(fields[0:6]))
        self.second = Bed6_entry(*(fields[6:12]))
        self.score = int(fields[12])
    def get_line(self):
        return "\t".join([self.first.get_line(),self.second.get_line(),str(self.score)])

class Bed13(Serial_Filetype):
    """0-based"""
    Entry_class = Bed13_entry
    def add_lines(self,lines):
        for line in lines:
            fields = line.strip().split('\t')
            if len(fields)>1:
                if len(fields)<13:
                    fields.append([None]*(13-len(fields)))
                fields[:] = [item if item!='.' else None for item in fields]
                self.entries.append(Bed13_entry(*fields))
    def get_lines(self):
        lines = []
        for entry in self.entries:
            lines.append(entry.get_line())
        return lines