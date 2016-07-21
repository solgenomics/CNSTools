from .serial_filetypes import Serial_Filetype
from .._utils import Progress_tracker

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