import abstract_handler as ah

class Entry(ah.Entry):
    """doc"""
    def __init__(self, chrom, chromStart, chromEnd, name=None, score=None, strand=None):
        #super(_Bed6_entry, self).__init__()
        self.chrom = chrom
        self.chromStart = int(chromStart)
        self.chromEnd = int(chromEnd)
        self.name = name
        try:
            self.score = float(score)
        except:
            self.score = None
        self.strand = strand
    def get_lines(self):
        strs = (str(i) if i!=None else '.' for i in (self.chrom, self.chromStart, self.chromEnd, self.name, self.score, self.strand))
        return ["\t".join(strs)]
class Handler(ah.Handler):
    def _entry_generator(self):
        with open(self.path,"r") as file_object:
            for line in file_object:
                fields = line.strip().split('\t')
                if len(fields)>1:
                    if len(fields)<6:
                        fields.append([None]*(6-len(fields)))
                    fields[:] = [item if item!='.' else None for item in fields]
                    yield Bed6_entry(*fields)