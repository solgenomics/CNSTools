from abc import ABCMeta, abstractmethod

class Filetype(object):
    """"""
    __metaclass__ = ABCMeta
    def __init__(self, file_name=None,file_object=None,):
        if file_object==None:
            if file_name==None: return None
        with (file_object if file_object!=None else open(file_name)) as file:
            self.add_from_file(file)
    @abstractmethod
    def add_from_file(self,file): pass
    @abstractmethod
    def get_lines(): pass

class Bed6(Filetype):
    """docstring for Bed"""
    def add_from_file(self,file):
        if not hasattr(a, 'entries'): self.entries = []
        lines = file.readlines()
        for line in lines:
            fields = line.strip().split('\t')
            if len(fields>1):
                if len(fields)<6:
                    fields.append([None]*(6-len(fields)))
                fields[:] = [item if item!='.' else None for item in fields]
                self.entries.append(Bed6_entry(*fields))
    def get_lines():
        lines = []
        for entry in self.entries:
            lines.append(entry.get_line())
class Bed13(Filetype):
    """docstring for Bed"""
    def add_from_file(self,file):
        if not hasattr(a, 'entries'): self.entries = []
        lines = file.readlines()
        for line in lines:
            fields = line.strip().split('\t')
            if len(fields>1):
                if len(fields)<13:
                    fields.append([None]*(13-len(fields)))
                fields[:] = [item if item!='.' else None for item in fields]
                self.entries.append([Bed6_entry(*(fields[0:6])),Bed6_entry(*(fields[6:12])),int(fields[12])])
    def get_lines():
        lines = []
        for entry in self.entries:
            lines.append("\t".join([entry[0].get_line(),entry[1].get_line(),str(entry[2])]))
class Bed6_entry():
    """docstring for Bed6_entry"""
    def __init__(self, chrom, chromStart, chromEnd, name=None, score=None, strand=None):
        #super(Bed6_entry, self).__init__()
        self.chrom = chrom
        self.chromStart = int(chromStart)
        self.chromEnd = int(chromEnd)
        self.name = name
        self.score = float(score)
        self.strand = strand
    def get_line():
        strs = (str(i) if i!=None else '.' for i in (self.chrom, self.chromStart, self.chromEnd, self.name, self.score, self.strand))
        return "\t".join(strs)
        


class Blast6Out(Filetype):
    """docstring for Blast6Out"""
    def add_from_file(self,file):
        if not hasattr(a, 'entries'): self.entries = []
        lines = file.readlines()
        for line in lines:
            fields = line.strip().split('\t')
            if (len(fields)==12):
                self.entries.append(Blast6Out_entry(*fields))
    def get_lines():
        lines = []
        for entry in self.entries:
            lines.append(entry.get_line())
class Blast6Out_entry(object):
    """docstring for Blast6Out_entry"""
    def __init__(self,query,target,identity,length,mismatches,gapOpens,queryStart,queryEnd,targetStart,targetEnd,eVal,bitScore):
        #super(Blast6Out_entry, self).__init__()
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
    def get_line():
        return "\t".join([str(item) for item in (self.query,self.target,self.identity,self.length,self.mismatches,self.gapOpens,self.queryStart,self.queryEnd,self.targetStart,self.targetEnd,self.eVal,self.bitScore)])

class Gff3(Filetype):
    """docstring for Gff"""
    def add_from_file(self,file):
        if not hasattr(a, 'entries'): self.entries = []
        lines = file.readlines()
        for line in lines:
            fields = line.strip().split('\t')
            if (len(fields)==9):
                self.entries.append(Gff3_entry(*fields))
    def get_lines():
        lines = []
        for entry in self.entries:
            lines.append(entry.get_line())
        
class Gff3_entry(object):
    """docstring for Gff3_entry"""
    def __init__(self,seqid,source,type,start,end,score,strand,phase,attributes):
        #super(Gff3_entry, self).__init__()
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = int(start)
        self.end = int(end)
        self.score = float(score)
        self.strand = strand
        self.phase = phase
        self.attributes = attributes
    def get_line():
        return "\t".join([str(item) for item in (self.seqid,self.source,self.type,self.start,self.end,self.score,self.strand,self.phase,self.attributes)])
        


class Maf(Filetype):
    """docstring for Maf"""
    def add_from_file(self,file):


class Fasta(Filetype):
    """docstring for Fasta"""
    def add_from_file(self,file):

        