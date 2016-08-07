from _utils import Progress_tracker,safe_print,reduce_gaps
from abc import ABCMeta, abstractmethod

class Filetype(object):
    """The base filetype class, this is an abstract class whose constructer takes either nothing,
     a filename, a file object, or lines then, if necessary, loads them into memory 
     and passes them to the add_lines abstract method.

     :param string file_name: Path to a file.
     :param object file_object: Any object with a :py:func:`readlines()` method.
     :param list[string] lines: A list of line strings."""
    __metaclass__ = ABCMeta
    def __init__(self,file_name=None,file_object=None,lines=None):
        if lines:
            self.add_lines(lines)
        elif file_object==None and file_name==None:
            return None
        else:
            with (file_object if file_object!=None else open(file_name)) as file:
                self.add_lines(file.readlines())
    def save_file(self,save_name):
        """Calls the :py:func:`Filetype.get_lines` method then saves the result to a file."""
        with open(save_name,"w") as out:
            out.write("\n".join(self.get_lines()))
    @abstractmethod
    def add_lines(self,lines):
        """(`abstract method`) Should take a list of line strings and then parse and store the data.

        :param list[string] lines: Lines in the filetype's format to add to stored data."""
        pass
    @abstractmethod
    def get_lines(self): 
        """(`abstract method`) Should return the entry data as a list of formatted line strings without line breaks.

        :returns: `list[str]`"""
        pass


class Serial_Filetype(Filetype):
    """This is an abstract subclass of :py:class:`Filetype` that sets up the common mechanisms
     for filetypes with contain multiple entries. Subclasses should set the class variable
     _Entry_class equal to the class of the entry obects which the subclass will use.

     :param str file_name: Path to a file.
     :param object file_object: Any object with a :py:func:`readlines()` method.
     :param list[string] lines: A list of line strings."""
    __metaclass__ = ABCMeta

    _Entry_class = None

    def __init__(self,file_name=None,file_object=None,lines=None):
        """Accepts a file object, file name, or list of line strings. Sets up an entry list (self.entries) and then calls the :func:'Filetype.__init__'"""
        self.entries = []
        super(Serial_Filetype, self).__init__(file_name,file_object,lines)

    def add_entry(self,*args,**kwargs): 
        """Creates a new entry by passing along the arguements to instantiate the _Entry_class then appends that new instance to self.entries. _Entry_class must be defined in the subclass for the add_entry Serial_Filetype function to work, otherwise it will return None.""" 
        if self._Entry_class:
            new_entry = self._Entry_class(*args,**kwargs) 
            self.entries.append(new_entry)
            return new_entry
        else: 
            return None

class _Wiggle_entry(object):
    """docstring for _Wiggle_entry"""
    def __init__(self, step_type, chrom, start, step, val_list):
        self.step_type = step_type
        self.chrom = chrom
        self.start = start
        self.step = step
        self.val_list = val_list[:]
    def get_lines(self):
        """ See :py:func:`Filetype.get_lines`.

        :returns: `list[str]`"""
        info_line = "%s chrom=%s start=%s step=%s" % (self.step_type,self.chrom,self.start,self.step)
        lines = [info_line]+[str(val) for val in self.val_list]
        return lines

class Wiggle(Serial_Filetype):
    """Filetype for loading a wiggle (.wig) file containing scores for loci. Currently only supports `fixedStep` with `step=1`.

     :param str file_name: Path to a file.
     :param object file_object: Any object with a :py:func:`readlines()` method.
     :param list[string] lines: A list of line strings."""
    _Entry_class = _Wiggle_entry
    _score_starting_chars = '.-1234567890'
    def add_lines(self,lines):
        """See :py:func:`Filetype.add_lines()`.

        :param list[string] lines: Lines in the filetype's format to add to stored data."""
        paragraphs = []
        found_first=False
        tracker = Progress_tracker("Parsing wiggle file",len(lines)*2).auto_display().start()
        for line in lines+['']:
            stripped = line.strip()
            if  (not stripped) or stripped[0] not in self._score_starting_chars:
                if not found_first: 
                    found_first=True
                else:
                    info = [item.split("=") for item in paragraphs[-1][0].split()]
                    step_type = info[0][0]
                    info_dict = {pair[0]:pair[1] for pair in info[1:]}
                    val_list = [float(item) for item in paragraphs[-1][1:]]
                    self.entries.append(_Wiggle_entry(step_type,info_dict['chrom'],int(info_dict['start']),int(info_dict['step']),val_list))
                    tracker.step(len(paragraphs[-1]))
                paragraphs.append([stripped])
            elif found_first:
                paragraphs[-1].append(stripped)
            tracker.step()
        tracker.done()      
    def to_bed(self,min_seg_length=7,min_seg_score=0.82,rejection_seg_len=12,rejection_seg_score=0.55):
        """Converts :py:class:`Wiggle` data to :py:class:`Bed6` data after trimming low score regions.

        :param int min_seg_length: Minimum length for a returned region.
        :param float min_seg_score: Minimum average locus score for a returned region.
        :param int rejection_seg_len: Minimum size of low score region to be sliced out of a region.
        :param float rejection_seg_score: Maximum score that is considered a candidate for removal in a contigous segment of at least length 'rejection_seg_len'."""

        new_bed = Bed6()
        for entry in self.entries:
            regions = [[0,len(entry.val_list)]]
            if len(entry.val_list) >= rejection_seg_len:
                start = 0
                end = rejection_seg_len-1
                below_rejec_score = [val<0.55 for val in entry.val_list[start:end]]
                #scans for regions of at least the rejection length in which all scores are lower that the rejection score
                #and cuts them out of the region
                while end<regions[-1][1]:
                    below_rejec_score.append(entry.val_list[end]<0.55)
                    end+=1
                    if all(below_rejec_score):
                        while end<regions[-1][1] and entry.val_list[end] < rejection_seg_score:
                            end+=1

                        cut_region = regions.pop()
                        regions.append([cut_region[0],start])
                        regions.append([end,cut_region[1]])

                        start = regions[-1][0]
                        end = regions[-1][0]+rejection_seg_len-1
                        below_rejec_score = [val<0.55 for val in entry.val_list[start:end]]
                        if regions[-1][1]-regions[-1][0] < min_seg_length:
                            break
                    below_rejec_score.pop(0)
                    start+=1
            for region in regions:
                region_sum = sum(val for val in entry.val_list[region[0]:region[1]])
                while region_sum < min_seg_score*(region[1]-region[0]) and region[1]-region[0] >= min_seg_length:
                    if entry.val_list[region[1]-1]<entry.val_list[region[0]]:
                        region_sum -= entry.val_list[region[1]-1]
                        region[1] = region[1]-1
                    else:
                        region_sum -= entry.val_list[region[0]]
                        region[0] = region[0]+1
                if region_sum >= min_seg_score*(region[1]-region[0]) and region[1]-region[0] >= min_seg_length:
                    new_bed.add_entry(entry.chrom, entry.start+region[0], entry.start+region[1], name=str(region_sum/float(region[1]-region[0])), score=(region[1]-region[0]), strand="+")
        return new_bed


    def get_lines(self):
        """ See :py:func:`Filetype.get_lines`.

        :returns: `list[str]`"""
        lines = []
        for entry in self.entries:
            lines+=entry.get_lines()
        return lines

class _Bed6_entry(object):
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
    def get_line(self):
        strs = (str(i) if i!=None else '.' for i in (self.chrom, self.chromStart, self.chromEnd, self.name, self.score, self.strand))
        return "\t".join(strs)

class Bed6(Serial_Filetype):
    """Filetype for loading a six field BED (.bed) file containing sequence regions.

    :param str file_name: Path to a file.
    :param object file_object: Any object with a :py:func:`readlines()` method.
    :param list[string] lines: A list of line strings."""
    _Entry_class = _Bed6_entry
    def add_lines(self,lines):
        """See :py:func:`Filetype.add_lines()`.

        :param list[string] lines: Lines in the filetype's format to add to stored data."""
        tracker = Progress_tracker("Parsing 6 column .bed",len(lines)).auto_display().start()
        for line in lines:
            fields = line.strip().split('\t')
            if len(fields)>1:
                if len(fields)<6:
                    fields.append([None]*(6-len(fields)))
                fields[:] = [item if item!='.' else None for item in fields]
                self.entries.append(_Bed6_entry(*fields))
            tracker.step()
        tracker.done()
    def get_lines(self):
        """ See :py:func:`Filetype.get_lines`.

        :returns: `list[str]`"""
        lines = []
        for entry in self.entries:
            lines.append(entry.get_line())
        return lines

class _Bed13_entry(object):
    def __init__(self,*fields):
        if len(fields)!=13:
            raise Exception("not 13 fields")
        self.first = _Bed6_entry(*(fields[0:6]))
        self.second = _Bed6_entry(*(fields[6:12]))
        self.score = int(fields[12])
    def get_line(self):
        return "\t".join([self.first.get_line(),self.second.get_line(),str(self.score)])

class Bed13(Serial_Filetype):
    """Filetype for loading a thirteen field BED (.bed) file containing sequence regions.

    :param str file_name: Path to a file.
    :param object file_object: Any object with a :py:func:`readlines()` method.
    :param list[string] lines: A list of line strings."""
    _Entry_class = _Bed13_entry
    def add_lines(self,lines):
        """See :py:func:`Filetype.add_lines()`.

        :param list[string] lines: Lines in the filetype's format to add to stored data."""
        for line in lines:
            fields = line.strip().split('\t')
            if len(fields)>1:
                if len(fields)<13:
                    fields.append([None]*(13-len(fields)))
                fields[:] = [item if item!='.' else None for item in fields]
                self.entries.append(_Bed13_entry(*fields))
    def get_lines(self):
        """ See :py:func:`Filetype.get_lines`.

        :returns: `list[str]`"""
        lines = []
        for entry in self.entries:
            lines.append(entry.get_line())
        return lines

class _Cns_sequence(object):
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
        return _Cns_sequence(self.genome,self.type,self.dist,self.loc_chrom,self.closest_gene,self.start,self.stop,self.gene_start,self.gene_stop,self.sequence,cns_ID=self.cns_ID)

class _Cns_entry(object):
    def __init__(self,cns_ID):
        self.cns_ID = cns_ID
        self.sequences = {}
    def add_seq(self,genome,*args):
        if not genome in self.sequences: self.sequences[genome] = []
        self.sequences[genome].append(_Cns_sequence(genome,*args,cns_ID=self.cns_ID))
        return self.sequences[genome][-1]
    def get_seqs(self,genome=None):
        if genome==None: return [seq for key in self.sequences for seq in self.sequences[key]]
        elif not genome in self.sequences: return None
        else: return self.sequences[genome][:]
    def get_lines(self):
        """ See :py:func:`Filetype.get_lines`.

        :returns: `list[str]`"""
        seqs = []
        for key in self.sequences:
            seqs+=self.sequences[key]
        return [seq.get_line() for seq in seqs]

class Cns(Serial_Filetype):
    """Filetype for loading a Conserved Noncoding Sequence (.cns) file, like the ones cnstools generates.

    :param str file_name: Path to a file.
    :param object file_object: Any object with a :py:func:`readlines()` method.
    :param list[string] lines: A list of line strings."""
    _Entry_class = _Cns_entry
    def add_lines(self,lines):
        """See :py:func:`Filetype.add_lines()`.

        :param list[string] lines: Lines in the filetype's format to add to stored data."""
        ID=None
        tracker = Progress_tracker("Parsing .cns",len(lines)).auto_display().start()
        for line in lines:
            list = [item if item!='.' else None for item in line.strip().split('\t')]
            if list[0]!=ID:
                ID = list[0]
                self.entries.append(_Cns_entry(ID))
            self.entries[-1].add_seq(*(list[1:]))
            tracker.step()
        tracker.done()
    def get_lines(self):
        """ See :py:func:`Filetype.get_lines`.

        :returns: `list[str]`"""
        return [line for entry in self.entries for line in entry.get_lines()]
    def to_fasta(self,genome):
        new_fasta = Fasta()
        tracker = Progress_tracker("Converting to .fasta",len(self.entries)).auto_display().start()
        for entry in self.entries:
            for seq in entry.get_seqs(genome):
                new_fasta.add_entry(entry.cns_ID,seq.sequence.replace("-",""))
            tracker.step()
        tracker.done()
        return new_fasta
    def to_bed(self,genome):
        new_bed = Bed6()
        tracker = Progress_tracker("Converting to .bed",len(self.entries)).auto_display().start()
        for entry in self.entries:
            for seq in entry.get_seqs(genome):
                new_bed.add_entry(seq.loc_chrom, seq.start, seq.stop, seq.cns_ID, seq.dist)
            tracker.step()
        tracker.done()
        return new_bed


class _Gff3_entry(object):
    """1-based"""
    def __init__(self,seqid,source,type,start,end,score,strand,phase,attributes):
        #super(_Gff3_entry, self).__init__()
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
    """Filetype for loading a GFF3 (.gff/.gff3) annotation file.

    :param str file_name: Path to a file.
    :param object file_object: Any object with a :py:func:`readlines()` method.
    :param list[string] lines: A list of line strings."""
    _Entry_class = _Gff3_entry
    def add_lines(self,lines):
        """See :py:func:`Filetype.add_lines()`.

        :param list[string] lines: Lines in the filetype's format to add to stored data."""
        tracker = Progress_tracker("Parsing .gff3",len(lines)).auto_display().start()
        for line in lines:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if (len(fields)==9):
                    self.entries.append(_Gff3_entry(*fields))
            tracker.step()
        tracker.done()
    def get_lines(self):
        """ See :py:func:`Filetype.get_lines`.

        :returns: `list[str]`"""
        lines = []
        for entry in self.entries:
            lines.append(entry.get_line())
        return lines
    def to_bed(self,type_list=None,genome=None):
        """Converts :py:class:`Gff3` data to :py:class:`Bed6` data after selecting for sequence types.

        :params list[string] type_list: List of sequence types which should be converted to BED regions.
        :params string genome: Genome name to prepend to the chromosome name when it is converted to BED."""
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

class _Maf_sequence(object):
    def __init__(self,src,start,size,strand,srcSize,text,metadata=None):
        self.src = src
        self.start = int(start)
        self.size = int(size)
        self.strand = strand
        self.srcSize = int(srcSize)
        self.text = text
        self.metadata = metadata
    def get_lines(self):
        """ See :py:func:`Filetype.get_lines`.

        :returns: `list[str]`"""
        lines = ['##'+self.metadata] if self.metadata else []
        lines.append('s '+'\t'.join([str(item) for item in [self.src,self.start,self.size,self.strand,self.srcSize,self.text]]))
        return lines
        
class _Maf_entry(object):
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
                    self.sequences.append(_Maf_sequence(*vals))
    def get_lines(self):
        """ See :py:func:`Filetype.get_lines`.

        :returns: `list[str]`"""
        lines = []
        if self.a_meta: lines.append('##'+self.a_meta)
        if self.a_line: lines.append(self.a_line)
        for sequence in self.sequences:
            lines+=sequence.get_lines()
        return lines

class Maf(Serial_Filetype):
    """0-based

    :param str file_name: Path to a file.
    :param object file_object: Any object with a :py:func:`readlines()` method.
    :param list[string] lines: A list of line strings."""
    _Entry_class = _Maf_entry

    def add_lines(self,lines):
        """See :py:func:`Filetype.add_lines()`.

        :param list[string] lines: Lines in the filetype's format to add to stored data."""
        if not hasattr(self, 'headerLines'): self.headerLines = []
        paragraph = []
        tracker = Progress_tracker("Parsing .maf",len(lines)).auto_display().start()
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("#") and not stripped.startswith("##--"):
                self.headerLines.append(stripped)
            elif stripped=="":
                if len(paragraph)>1:
                    self.entries.append(_Maf_entry(paragraph))
                paragraph = []
            else:
                paragraph.append(stripped)
            tracker.step()
        tracker.done()

    def get_lines(self):
        """ See :py:func:`Filetype.get_lines`.

        :returns: `list[str]`"""
        lines = []
        for entry in self.entries:
            lines+=entry.get_lines()
            lines.append("")
        return lines

    def to_bed(self,genome_name=None,index_tag="maf_index"):
        """Converts :py:class:`Maf` data to :py:class:`Bed6` using the specified genome for location information.
        
        :params string genome_name: Name of the genome to use to make BED regions. If `genome_name==None` then it will use the first genome it encounters in the file.
        :params string index_tag: Tag name which will be use to mark the BED entry with the MAF entry's index."""
        new_bed = Bed6()
        tracker = Progress_tracker("Converting to .bed",len(self.entries)).auto_display().start()
        if not genome_name: 
            genome_name=self.entries[0].sequences[0].src.split(":")[0].strip()
        index = 0
        for entry in self.entries:
            seq_to_convert = (seq for seq in entry.sequences if seq.src.split(":")[0].strip()==genome_name)
            for sequence in seq_to_convert:
                id_string = "%s=%s" % (index_tag,index) if index_tag!=None else None
                if sequence.metadata:
                    id_string = (sequence.metadata+";"+id_string) if id_string else (sequence.metadata)
                new_bed.add_entry(sequence.src, sequence.start, sequence.start+sequence.size, name=id_string, strand=sequence.strand)
            index+=1
            tracker.step()
        tracker.done()
        return new_bed

    def slice_with_bed(self,bed,ref_genome,index_tag,max_N_ratio=0.5,max_gap_ratio=0.5,min_len=15):
        """Intersects the regions from a :py:class:`Bed6` object with the :py:class:`Maf` entries. 
        Also filters for unwanted sequences.
        
        :param bed: BED object to use regions from.
        :type bed: :py:class:`Bed6`
        :param string ref_genome: Genome to use MAF regions from.
        :param string index_tag: Tag name in the BED entry which conatins the relevant MAF entry's index.
        :param float max_N_ratio: Maximum ratio of Ns to ATCG nucleotides.
        :param float max_gap_ratio: Maximum ratio of gaps to nucleotides of a sequence in an alignment.
        :param int min_len: Minimum length of a sequence in an alignment.
        """
        new_maf = Maf()
        new_maf.headerlines = self.headerLines[:]
        for bed_entry in bed.entries:
            try:
                keyVal = (kv for kv in bed_entry.name.split(";") if kv.find(index_tag)!=-1).next()
            except StopIteration:
                raise ValueError("The provided index_tag was not found in a bed entry.")

            parent_index = int(keyVal.split(index_tag+"=")[1])
            original_maf_entry = self.entries[parent_index]
            new_maf_entry = _Maf_entry()
            new_maf_entry.a_line = original_maf_entry.a_line
            new_maf_entry.a_meta = original_maf_entry.a_meta

            try:
                mainseq = (seq for seq in original_maf_entry.sequences if seq.src.split(":")[0].strip()==ref_genome).next()
            except StopIteration:
                raise ValueError("The provided ref_genome was not found in an alignment.")

            front_offset = bed_entry.chromStart - mainseq.start
            back_offset = mainseq.start+mainseq.size - bed_entry.chromEnd
            front_cut_num = Maf._gap_cut_loc(mainseq.text,front_offset)
            back_cut_num = Maf._gap_cut_loc(reversed(mainseq.text),back_offset)

            for seq in original_maf_entry.sequences:
                front_removed_seq = seq.text[:front_cut_num]
                new_seq = (seq.text[front_cut_num:-back_cut_num]) if back_cut_num!=0 else (seq.text[front_cut_num:])
                new_len = Maf._no_gap_len(new_seq)
                new_start = seq.start + Maf._no_gap_len(front_removed_seq)
                is_ref_seq = seq.src.split(":")[0].strip()==ref_genome
                is_within_params = Maf._no_gap_len(new_seq) >= min_len and \
                                   1-(new_len/float(len(new_seq))) <= max_gap_ratio and \
                                   1-Maf._no_gap_len(new_seq.replace('N',''))/float(new_len) <= max_N_ratio
                if is_ref_seq or is_within_params:    
                    new_seq = _Maf_sequence(
                        seq.src,
                        new_start,
                        new_len,
                        seq.strand,
                        seq.srcSize,
                        new_seq,
                        seq.metadata)
                    new_maf_entry.sequences.append(new_seq)
            if sum([int(seq.size >= min_len) for seq in new_maf_entry.sequences])>=2:
                new_seq_text = reduce_gaps([seq.text for seq in new_maf_entry.sequences])
                for i in range(len(new_maf_entry.sequences)):
                    new_maf_entry.sequences[i].text = new_seq_text[i]
                new_maf.entries.append(new_maf_entry)
        return new_maf

    def cns_from_proxim_beds(self,cns_proxim_beds,index_tag="maf_index"):
        """Creates a :py:class:`Cns` object from a :py:class:`Maf` object and multiple :py:class:`Bed13` object conatining gene proximity information for each genome. 
        Also filters for unwanted sequences.
        
        :param cns_proxim_beds: BED objects to use proximity info from.
        :type cns_proxim_beds: dict{genome_name: :py:class:`Bed6`}
        :param string index_tag: Tag name in the BED entry which conatins the relevant MAF entry's index.
        """
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
        cns.entries = [entry for entry in cns.entries if len(entry.sequences)>1]
        for entry in cns.entries:
            cns_seq_keys = [(key,len(entry.sequences[key])) for key in entry.sequences]
            texts = reduce_gaps([seq.sequence for key,count in cns_seq_keys for seq in entry.sequences[key]])[::-1]
            for key,count in cns_seq_keys:
                for i in range(count):
                    entry.sequences[key][i].sequence = texts.pop()
        return cns

    @staticmethod
    def _gap_cut_loc(seq,cutNum):
        anonymized = "".join(["x" if i!="-" else i for i in seq]) #use x to represent valid chars
        eliminated = anonymized.replace("x","o",cutNum) #change the first n x's to o's
        return eliminated.find("x") # finds the index of the next x, which is how many chars to cut off
    @staticmethod
    def _no_gap_len(seq): return len(seq[:].replace("-",""))

# class Tomtom_match(object):
#     def __init__(self,query_id,target_id,optimal_offset,p_value,e_value,q_value,overlap,query_consensus,target_consensus,orientation):
#         self.query_id = query_id
#         self.target_id = target_id
#         self.optimal_offset = optimal_offset
#         self.p_value = p_value
#         self.e_value = e_value
#         self.q_value = q_value
#         self.overlap = overlap
#         self.query_consensus = query_consensus
#         self.target_consensus = target_consensus
#         self.orientation = orientation
#     def get_line(self):
#         return "\t".join((self.query_id,self.target_id,self.optimal_offset,self.p_value,self.e_value,self.q_value,self.overlap,self.query_consensus,self.target_consensus,self.orientation))

# class Tomtom_entry(object):
#     def __init__(self,lines):
#         lists = [line.split("\t") for line in lines]
#         self.query_id = lists[0][0]
#         self.query_consensus = lists[0][7]
#         self.matches = [Tomtom_match(*list) for list in lists]
#     def get_lines(self):
#         return [match.get_line() for match in self.matches]

# class Tomtom(Serial_Filetype):
#     """docstring for Tomtom"""
#     _Entry_class = Tomtom_entry
#     def add_lines(self,lines):
#         lines = [line.strip() for line in lines]
#         paragraphs = []
#         self.header = []
#         last_ID = None
#         for line in lines:
#             if not line.startswith("#"):
#                 this_ID = line[:line.find('\t')]
#                 if this_ID!=last_ID:
#                     last_ID = this_ID
#                     paragraphs.append([])
#                 if not this_ID==None:
#                     paragraphs[-1].append(line)
#             else:
#                 self.header.append(line)
#         for paragraph in paragraphs:
#             if(len(paragraph)>0):
#                 self.entries.append(Tomtom_entry(paragraph))
#     def get_lines(self):
#         lines = self.header
#         for entry in self.entries:
#             lines+= entry.get_lines()
#         return lines

# class Meme_v_4_entry(object):
#     def __init__(self,lines):
#         start_index = 0
#         while not lines[start_index].startswith('MOTIF'):
#             start_index+=1
#         idList = lines[start_index].split()
#         self.identifier = idList[1]
#         self.alt_name = idList[2]
#         self.lines = lines
#     def get_lines(self):
#         return self.lines

# class Meme_v_4(Serial_Filetype):
#     _Entry_class = Meme_v_4_entry
#     def __init__(self,*args,**kwargs):
#         if not hasattr(self,"entry_dict"):
#             self.entry_dict = {}
#         if not hasattr(self,"header"):
#             self.header = []
#         super(Meme_v_4, self).__init__(*args,**kwargs)
#     def add_lines(self,lines):
#         lines = [line.strip() for line in lines]
#         header_end = 0
#         while not lines[header_end].startswith('MOTIF'):
#             header_end+=1
#         self.header+= lines[:header_end]
#         paragraphs = []
#         first_found = False
#         for line in lines[header_end:]:
#             if line.startswith('MOTIF'):
#                 paragraphs.append([])
#                 first_found = True
#             if first_found:
#                 paragraphs[-1].append(line)
#         for paragraph in paragraphs:
#             if(len(paragraph)>1):
#                 self.entries.append(self._Entry_class(paragraph))
#                 self.entry_dict[self.entries[-1].identifier] = self.entries[-1]
#     def get_lines(self):
#         lines = self.header
#         for entry in self.entries:
#             lines+= entry.get_lines()
#         return lines
#     def add_entry(self,lines=None,entry=None): 
#         if entry!=None:
#             self.entries.append(entry)
#         elif lines!=None:
#             self.entries.append(self._Entry_class(lines))
#         else:
#             return None
#         self.entry_dict[self.entries[-1].identifier] = self.entries[-1]
#         return self.entries[-1]
#     def lookup(self,identifier):
#         if identifier in self.entry_dict:
#             return self.entry_dict[identifier]
#         else:
#             return None

# class BlastF6_entry(object):
#     """1-based"""
#     def __init__(self,query,target,identity,length,mismatches,gapOpens,queryStart,queryEnd,targetStart,targetEnd,eVal,bitScore):
#         #super(BlastF6_entry, self).__init__()
#         self.query = query
#         self.target = target
#         self.identity = float(identity)
#         self.length = int(length)
#         self.mismatches = int(mismatches)
#         self.gapOpens = int(gapOpens)
#         self.queryStart = int(queryStart)
#         self.queryEnd = int(queryEnd)
#         self.targetStart = int(targetStart)
#         self.targetEnd = int(targetEnd)
#         self.eVal = float(eVal)
#         self.bitScore = float(bitScore)
#     def get_line(self):
#         return "\t".join([str(item) for item in (self.query,self.target,self.identity,self.length,self.mismatches,self.gapOpens,self.queryStart,self.queryEnd,self.targetStart,self.targetEnd,self.eVal,self.bitScore)])

# class BlastF6(Serial_Filetype):
#     """1-based"""
#     _Entry_class = BlastF6_entry
#     def add_lines(self,lines):
#         tracker = Progress_tracker("Parsing blast output",len(lines)).auto_display().start()
#         for line in lines:
#             fields = line.strip().split('\t')
#             if (len(fields)==12):
#                 self.entries.append(BlastF6_entry(*fields))
#             tracker.step()
#         tracker.done()
#     def get_lines(self):
#         lines = []
#         for entry in self.entries:
#             lines.append(entry.get_line())
#         return lines
#     def to_bed(self):
#         new_bed = Bed6()
#         tracker = Progress_tracker("Converting to .bed",len(self.entries)).auto_display().start()
#         for entry in self.entries:
#             strand = None
#             if entry.targetStart < entry.targetEnd:
#                 start, end = entry.targetStart, entry.targetEnd
#                 strand = '+'
#             else:
#                 start, end = entry.targetEnd, entry.targetStart
#                 strand = '-'
#             #convert to 0-based!
#             new_bed.add_entry(entry.target, start-1, end, entry.query, entry.bitScore, strand)
#             tracker.step()
#         tracker.done()
#         return new_bed

class Fasta_entry(object):
    def __init__(self,description,sequence):
        self.description = description
        self.sequence = sequence
    def get_lines(self):
        return [">"+self.description]+[self.sequence[i:i+70] for i in range(0,len(self.sequence),70)]

class Fasta(Serial_Filetype):
    """docstring for Fasta"""
    _Entry_class = Fasta_entry
    def add_lines(self,lines):
        """See :py:func:`Filetype.add_lines()`.

        :param list[string] lines: Lines in the filetype's format to add to stored data."""
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
        """ See :py:func:`Filetype.get_lines`.

        :returns: `list[str]`"""
        lines = []
        for entry in self.entries:
            lines+= entry.get_lines()
        return lines

# class Gff2(Gff3):
#     def add_lines(self,lines):
#         tracker = Progress_tracker("Parsing .gff3",len(lines)).auto_display().start()
#         for line in lines:
#             if not line.startswith('#'):
#                 fields = [i for i in line.strip().split('\t') if i!='']
#                 if (len(fields)==8):
#                     fields = fields[:2]+["?"]+fields[2:] #unknown type!
#                     self.entries.append(_Gff3_entry(*fields))
#             tracker.step()
#         tracker.done()
