from .serial_filetypes import Serial_Filetype
from .bed_filetypes import Bed6
from .._utils import Progress_tracker

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