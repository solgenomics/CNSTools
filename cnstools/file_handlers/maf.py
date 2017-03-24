import abstract_handler as ah
from .._utils import MultiTracker
import bed6
#

class Sequence(ah.Entry):
    def __init__(self,src,start,size,strand,srcSize,text):
        self.src = src
        self.start = int(start)
        self.size = int(size)
        self.strand = strand
        self.srcSize = int(srcSize)
        self.text = text
    def get_line(self):
        """"""
        return lines ('s '+'\t'.join([str(item) for item in [self.src,self.start,self.size,self.strand,self.srcSize,self.text]]))
        

class Entry(ah.Entry):
    """0-based"""
    def __init__(self,paragraph=None):
        self.a_line = None
        self.sequences = []
        if paragraph:
            for line in paragraph:
                line=line.strip()
                if line.startswith('a'):
                    self.a_line = line
                elif line.startswith('s'):
                    vals = line.split()[1:]
                    self.sequences.append(_Maf_sequence(*vals))
    def get_lines(self):
        """ See :py:func:`Filetype.get_lines`.

        :returns: `list[str]`"""
        lines = []
        if self.a_line: lines.append(self.a_line)
        for sequence in self.sequences:
            lines.appen(sequence.get_line())
        return lines


class Comment(ah.Comment):
    def __init__(self,content):
        self.content = content.strip()
    def get_lines(self):
        return self.content.split("\n")

class Handler(ah.Handler):
    """0-based"""
    def __init__(self, path):
        super(Handler, self).__init__(path)
        self._modify = self.modify
        self.modify = self._modify_wrap

    def _modify_wrap(self,*args,**kwargs):
        kwargs["include_comments"]=True
        kwargs["add_end_break"]=True
        return self._modify(*args,**kwargs)

    def _entry_generator(self,include_comments=False,parent=None,tracker_name=None):
        if tracker_name:
            size = os.stat(self.path).st_size
            if parent!=None:
                tracker = parent.subTracker(tracker_name,size,estimate=False,style="percent")
            else:
                tracker = MultiTracker(tracker_name,size,estimate=False,style="percent").auto_display(1)
        with open(self.path,"r") as file_object:
            paragraph = []
            for line in file_object:
                if tracker_name:
                    tracker.step(len(line))
                stripped = line.strip()
                if stripped.startswith("#") and include_comments:
                    yield Comment(stripped)
                elif stripped=="":
                    if len(paragraph)>1:
                        yield Entry(paragraph)
                    paragraph = []
                else:
                    paragraph.append(stripped)
            if len(paragraph)>1:
                yield Entry(paragraph)
        if tracker_name:
            tracker.done()

    def to_bed(self,path,genome_name=None,index_tag=None,parent=None,tracker_name=None):
        """"""
        bed = self._modify(_to_bed_entry(genome_name,index_tag), 
                           path, target_type=bed6.Handler, 
                           add_end_break=False,
                           include_comments=False,
                           parent=parent,
                           tracker_name=tracker_name)
        return bed

    def bed_intersect(self,bed,path,genome_name,index_tag,max_N_ratio=0.5,max_gap_ratio=0.5,min_len=15,parent=None,tracker_name=None):
        """Intersects the regions from a :py:class:`Bed6` object with the :py:class:`Maf` entries. 
        Also filters for unwanted sequences.
        
        :param bed: BED object to use regions from.
        :type bed: :py:class:`Bed6`
        :param string genome_name: Genome to use MAF regions from.
        :param string index_tag: Tag name in the BED entry which conatins the relevant MAF entry's index.
        :param float max_N_ratio: Maximum ratio of Ns to ATCG nucleotides.
        :param float max_gap_ratio: Maximum ratio of gaps to nucleotides of a sequence in an alignment.
        :param int min_len: Minimum length of a sequence in an alignment.
        """
        newMaf = self._modify(_bed_intersect(bed,genome_name,index_tag,max_N_ratio,max_gap_ratio,min_len), 
                           path, target_type=Handler, 
                           add_end_break=True,
                           include_comments=True,
                           parent=parent,
                           tracker_name=tracker_name)
        return newMaf

def _bed_intersect(bed,genome_name,index_tag,max_N_ratio,max_gap_ratio,min_len):
    bed_entry_iterator = bed._entry_generator()
    bed_empty = False
    try:
        bed_entry = bed_entry_iterator.next()
        bed_id_string = (kv for kv in bed_entry.name.split(";") if kv.startswith(index_tag))[0]
        bed_index = int(bed_id_string.split("=")[1])
    except StopIteration:
        bed_empty=True
    maf_index = 0
    def maf_entry_filter_func(maf_entry):
        if bed_empty:
            return False
        if isinstance(maf_entry, Comment):
            return maf_entry
        else:
            if maf_index>bed_index:
                while maf_index>bed_index:
                    try:
                        bed_entry = bed_entry_iterator.next()
                        bed_id_string = (kv for kv in bed_entry.name.split(";") if kv.startswith(index_tag))[0]
                        bed_index = int(bed_id_string.split("=")[1])
                    except StopIteration:
                        bed_empty=True
                        return False
            if maf_index!=bed_index: raise Error("This should not have happened, bed has outrun maf.")
        new_maf_entries = []
        while maf_index==bed_index:
            new_maf_entry = Entry()
            new_maf_entry.a_line = maf_entry.a_line
            new_maf_entry.a_meta = maf_entry.a_meta

            try:
                mainseq = (seq for seq in maf_entry.sequences if seq.src.split(":")[0].strip()==ref_genome).next()
            except StopIteration:
                raise ValueError("The provided ref_genome was not found in an alignment.")

            front_offset = bed_entry.chromStart - mainseq.start
            back_offset = mainseq.start+mainseq.size - bed_entry.chromEnd
            front_cut_num = _gap_cut_loc(mainseq.text,front_offset)
            back_cut_num = _gap_cut_loc(reversed(mainseq.text),back_offset)

            for seq in maf_entry.sequences:
                front_removed_seq = seq.text[:front_cut_num]
                new_seq = (seq.text[front_cut_num:-back_cut_num]) if back_cut_num!=0 else (seq.text[front_cut_num:])
                new_len = _no_gap_len(new_seq)
                new_start = seq.start + _no_gap_len(front_removed_seq)
                is_ref_seq = seq.src.split(":")[0].strip()==ref_genome
                is_within_params = (_no_gap_len(new_seq) >= min_len) and \
                                   (1-(new_len/float(len(new_seq))) <= max_gap_ratio) and \
                                   (1-_no_gap_len(new_seq.replace('N',''))/float(new_len) <= max_N_ratio)
                if is_ref_seq or is_within_params:    
                    new_seq = Sequence(
                        seq.src,
                        new_start,
                        new_len,
                        seq.strand,
                        seq.srcSize,
                        new_seq)
                    new_maf_entry.sequences.append(new_seq)
            if sum([int(seq.size >= min_len) for seq in new_maf_entry.sequences])>=2:
                new_seq_text = _reduce_gaps([seq.text for seq in new_maf_entry.sequences])
                for i in range(len(new_maf_entry.sequences)):
                    new_maf_entry.sequences[i].text = new_seq_text[i]
                new_maf_entries.append(new_maf_entry)
            try:
                bed_entry = bed_entry_iterator.next()
                bed_id_string = (kv for kv in bed_entry.name.split(";") if kv.startswith(index_tag))[0]
                bed_index = int(bed_id_string.split("=")[1])
            except StopIteration:
                bed_empty=True
                break;
        return new_maf_entries

def _to_bed_entry(genome_name=None,index_tag=None):
    index = 0
    def func(entry):
        if genome_name==None: 
            genome_name = entry.sequences[0].src.split(":")[0].strip()
        seq_to_convert = (seq for seq in entry.sequences if seq.src.split(":")[0].strip()==genome_name)
        for sequence in seq_to_convert:
            id_string = "%s=%s" % (index_tag,index) if index_tag!=None else None
            index+=1
            new_entry = Entry(sequence.src, sequence.start, sequence.start+sequence.size, name=id_string, strand=sequence.strand)
        return new_entry
    return func

def _gap_cut_loc(seq,cutNum):
    anonymized = "".join(["x" if i!="-" else i for i in seq]) #use x to represent valid chars
    eliminated = anonymized.replace("x","o",cutNum) #change the first n x's to o's
    return eliminated.find("x") # finds the index of the next x, which is how many chars to cut off
def _no_gap_len(seq): 
    return len(seq[:].replace("-",""))
def _reduce_gaps(seqlist):
    lists = [list(seq) for seq in seqlist]
    for i in range(len(lists[0])-1,-1,-1):
        if all(seq[i]=="-" for seq in lists):
            for seq in lists: 
                del seq[i]
    return ["".join(seq) for seq in lists]



class _old(object):
    def cns_from_proxim_beds(self,cns_proxim_beds,index_tag="maf_index",prefix=""):
        """Creates a :py:class:`Cns` object from a :py:class:`Maf` object and multiple :py:class:`Bed13` object conatining gene proximity information for each genome. 
        Also filters for unwanted sequences.
        
        :param cns_proxim_beds: BED objects to use proximity info from.
        :type cns_proxim_beds: dict{genome_name: :py:class:`Bed6`}
        :param string index_tag: Tag name in the BED entry which conatins the relevant MAF entry's index.
        """
        cns = Cns()
        for i in range(len(self.entries)):
            cns.add_entry(prefix+str(i))
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
        cns.entries = [entry for entry in cns.entries if len(entry.sequences)>1]
        for entry in cns.entries:
            cns_seq_keys = [(key,len(entry.sequences[key])) for key in entry.sequences]
            texts = reduce_gaps([seq.sequence for key,count in cns_seq_keys for seq in entry.sequences[key]])[::-1]
            for key,count in cns_seq_keys:
                for i in range(count):
                    entry.sequences[key][i].sequence = texts.pop()
        return cns
