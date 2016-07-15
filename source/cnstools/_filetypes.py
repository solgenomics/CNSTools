"""Module which contains classes for storing the data for each filetype used in cnstools. The filetype classes also have some conversion functions for changing file formats."""
from abc import ABCMeta, abstractmethod
from _utils import Progress_tracker

class Serial_Filetype(object):
    """Abstract parent class of all of the filetype functions. The __init__ function here should be called by each subclass."""
    __metaclass__ = ABCMeta

    Entry_class = None

    def __init__(self,file_name=None,file_object=None,lines=None):
        """Accepts a file object, file name, or list of line strings. Sets up an entry list (self.entries) and then calls the abstract add_lines() function with a list of lines to be parsed into entrys."""
        self.entries = []
        if lines:
            self.add_lines(lines)
        elif file_object==None and file_name==None:
            return None
        else:
            with (file_object if file_object!=None else open(file_name)) as file:
                self.add_lines(file.readlines())

    def add_entry(self,*args,**kwargs): 
        """Creates a new entry by passing along the arguements to instantiate the Entry_class then appends that new instance to self.entries. Entry_class must be defined in the subclass for the add_entry Serial_Filetype function to work, otherwise it will return None.""" 
        if self.Entry_class:
            new_entry = self.Entry_class(*args,**kwargs) 
            self.entries.append(new_entry)
            return new_entry
        else: 
            return None

    def save_file(self,save_name):
        with open(save_name,"w") as out:
            out.write("\n".join(self.get_lines()))

    @abstractmethod
    def add_lines(self,file): 
        """Abstract. Should take a list of line strings and parse them into instances of the Entry_class then append them to self.entries"""
        pass
    @abstractmethod
    def get_lines(self): 
        """Abstract. Should return the entry data formatted as a list of line strings"""
        pass

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
    def to_bed(self,type_list=None):
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
            new_bed.add_entry(entry.seqid, chromStart-1, chromEnd, name=id_with_type, score=entry.score, strand=entry.strand)
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
    def to_bed(self,seq_name=None):
        new_bed = Bed6()
        tracker = Progress_tracker("Converting to .bed",len(self.entries)).auto_display().start()
        if not seq_name: 
            seq_name=self.entries[0].sequences[0].src
        for entry in self.entries:
            for sequence in (seq for seq in entry.sequences if seq.src==seq_name):
                new_bed.add_entry(sequence.src, sequence.start, sequence.start+sequence.size, name=sequence.metadata, strand=sequence.strand)
            tracker.step()
        tracker.done()
        return new_bed

class Fasta_entry(object):
    def __init__(self,description,sequence):
        self.description = description
        self.sequence = sequence
    def get_lines(self):
        return [">"+self.description]+[self.sequence[i:i+70] for i in range(0,len(self.sequence),70)]

class Fasta(Serial_Filetype):
    """docstring for Fasta"""
    Entry_class = Fasta_entry
    def add_lines(self,lines):
        paragraphs = [[]]
        first_found = False
        tracker = Progress_tracker("Parsing .fasta data",len(lines*2)).auto_display().start()
        for line in lines:
            stripped = line.strip()
            if stripped.startswith('>'):
                first_found = True
                paragraphs[-1].append(stripped[1:])
            elif first_found:
                paragraphs[-1].append(stripped)
            tracker.step()
        for paragraph in paragraphs:
            for item in paragraph[1:]:
                if not type(item) is str:
                    raise TypeError (str(item))
            self.entries.append(Fasta_entry(paragraph[0],"".join(paragraph[1:])))
            tracker.step(len(paragraph))
        tracker.done()
    def get_lines(self):
        lines = []
        for entry in self.entries:
            lines+= entry.get_lines()
        return lines

class Tomtom_match(object):
    def __init__(self,query_id,target_id,optimal_offset,p_value,e_value,q_value,overlap,query_consensus,target_consensus,orientation):
        self.query_id = query_id
        self.target_id = target_id
        self.optimal_offset = optimal_offset
        self.p_value = p_value
        self.e_value = e_value
        self.q_value = q_value
        self.overlap = overlap
        self.query_consensus = query_consensus
        self.target_consensus = target_consensus
        self.orientation = orientation
    def get_line(self):
        return "\t".join((self.query_id,self.target_id,self.optimal_offset,self.p_value,self.e_value,self.q_value,self.overlap,self.query_consensus,self.target_consensus,self.orientation))

class Tomtom_entry(object):
    def __init__(self,lines):
        lists = [line.split("\t") for line in lines]
        self.query_id = lists[0][0]
        self.query_consensus = lists[0][7]
        self.matches = [Tomtom_match(*list) for list in lists]
    def get_lines(self):
        return [match.get_line() for match in self.matches]

class Tomtom(Serial_Filetype):
    """docstring for Tomtom"""
    Entry_class = Tomtom_entry
    def add_lines(self,lines):
        lines = [line.strip() for line in lines]
        paragraphs = []
        self.header = []
        last_ID = None
        for line in lines:
            if not line.startswith("#"):
                this_ID = line[:line.find('\t')]
                if this_ID!=last_ID:
                    last_ID = this_ID
                    paragraphs.append([])
                if not this_ID==None:
                    paragraphs[-1].append(line)
            else:
                self.header.append(line)
        for paragraph in paragraphs:
            if(len(paragraph)>0):
                self.entries.append(Tomtom_entry(paragraph))
    def get_lines(self):
        lines = self.header
        for entry in self.entries:
            lines+= entry.get_lines()
        return lines

class Meme_v_4_entry(object):
    def __init__(self,lines):
        start_index = 0
        while not lines[start_index].startswith('MOTIF'):
            start_index+=1
        idList = lines[start_index].split()
        self.identifier = idList[1]
        self.alt_name = idList[2]
        self.lines = lines
    def get_lines(self):
        return self.lines

class Meme_v_4(Serial_Filetype):
    Entry_class = Meme_v_4_entry
    def __init__(self,*args,**kwargs):
        if not hasattr(self,"entry_dict"):
            self.entry_dict = {}
        if not hasattr(self,"header"):
            self.header = []
        super(Meme_v_4, self).__init__(*args,**kwargs)
    def add_lines(self,lines):
        lines = [line.strip() for line in lines]
        header_end = 0
        while not lines[header_end].startswith('MOTIF'):
            header_end+=1
        self.header+= lines[:header_end]
        paragraphs = []
        first_found = False
        for line in lines[header_end:]:
            if line.startswith('MOTIF'):
                paragraphs.append([])
                first_found = True
            if first_found:
                paragraphs[-1].append(line)
        for paragraph in paragraphs:
            if(len(paragraph)>1):
                self.entries.append(self.Entry_class(paragraph))
                self.entry_dict[self.entries[-1].identifier] = self.entries[-1]
    def get_lines(self):
        lines = self.header
        for entry in self.entries:
            lines+= entry.get_lines()
        return lines
    def add_entry(self,lines=None,entry=None): 
        if entry!=None:
            self.entries.append(entry)
        elif lines!=None:
            self.entries.append(self.Entry_class(lines))
        else:
            return None
        self.entry_dict[self.entries[-1].identifier] = self.entries[-1]
        return self.entries[-1]
    def lookup(self,identifier):
        if identifier in self.entry_dict:
            return self.entry_dict[identifier]
        else:
            return None
        

        