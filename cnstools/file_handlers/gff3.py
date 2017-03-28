import os
import abstract_handler as ah
from .._utils import MultiTracker
import bed6

class Entry(ah.Entry):
    """1-based"""
    def __init__(self,seqid,source,seq_type,start,end,score,strand,phase,attributes):
        #super(_Gff3_entry, self).__init__()
        self.seqid = seqid
        self.source = source
        self.seq_type = seq_type
        self.start = int(start)
        self.end = int(end)
        try:
            self.score = float(score)
        except:
            self.score = None
        self.strand = strand
        self.phase = phase
        self.attributes = attributes
    def get_lines(self):
        return ["\t".join([str(item) for item in (self.seqid,self.source,self.seq_type,self.start,self.end,self.score,self.strand,self.phase,self.attributes)])]

class Comment(ah.Comment):
    def __init__(self,content):
        self.content = content.strip()
    def get_lines(self):
        return self.content.split("\n")

class Handler(ah.Handler):

    def __init__(self, path):
        super(Handler, self).__init__(path)
        self._modify = self.modify
        self.modify = self._modify_wrap

    def _modify_wrap(self,*args,**kwargs):
        kwargs["include_comments"]=True
        return self._modify(*args,**kwargs)

    def _entry_generator(self,include_comments=False,parent=None,tracker_name=None):
        if tracker_name:
            size = os.stat(self.path).st_size
            if parent!=None:
                tracker = parent.subTracker(tracker_name,size,estimate=False,style="percent")
            else:
                tracker = MultiTracker(tracker_name,size,estimate=False,style="percent").auto_display(1)
        with open(self.path,"r") as file_object:
            for line in file_object:
                if tracker_name:
                    tracker.step(len(line))
                if not line.startswith('#'):
                    fields = line.strip().split('\t')
                    if (len(fields)==9):
                        yield Entry(*fields)
                        continue
                if include_comments: yield Comment(line)
        if tracker_name:
            tracker.done()

    def to_bed(self,path,type_list=None,genome=None):
        """"""
        bed = self._modify(_to_bed_entry(type_list,genome), 
                           path, target_type=bed6.Handler, 
                           add_end_break=False,
                           include_comments=False)
        return bed

def _to_bed_entry(type_list,genome):
    def func(entry):
        if type_list!=None:
            if entry.seq_type not in type_list:
                return False
            if(entry.start<entry.end):
                chromStart,chromEnd = entry.start,entry.end
            else:
                chromStart,chromEnd = entry.end,entry.start
            id_with_type = entry.attributes+";seqType="+entry.seq_type
            chrom = entry.seqid if genome==None else genome+":"+entry.seqid
            return bed6.Entry(chrom, chromStart-1, chromEnd, name=id_with_type, score=entry.score, strand=entry.strand)
    return func