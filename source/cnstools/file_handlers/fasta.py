import abstract_handler as ah

import os

class Entry(ah.Entry):
    def __init__(self, description, data, line_length):
        self.description,self.data,self.line_length = description,data,line_length

    def get_lines(self):
        lines = [">"+self.description]
        lines+= [self.data[i:i+self.line_length] for i in range(0,len(self.data),self.line_length)]
        return lines

class Handler(ah.Handler):
    in_place = object()
    accepted_symbols = set(list("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ*-"))

    def __init__(self, path,line_length=80):
        self.line_length = line_length
        super(Handler, self).__init__(path)

    def _entry_generator(self,tracker=None):
        with open(self.path,"r") as file_object:
            entry_description = None
            entry_data = ""
            for line in file_object:
                if tracker: tracker.step(len(line))
                line = line.strip()
                if line.startswith(">"):
                    if entry_description!=None:
                        yield Fasta_entry(entry_description,entry_data,self.line_length)
                    entry_data = ""
                    entry_description = line[1:]
                    continue
                elif entry_description!=None:
                    entry_data+="".join([symbol for symbol in line if symbol in self.accepted_symbols])
            yield Fasta_entry(entry_description,entry_data)
            if tracker: tracker.done()

    def split(self,num_per_file=1,out_folder=in_place,file_prefix=in_place,file_suffix="",file_extension=".fa",tracker=None):
        file_prefix = (file_prefix) if file_prefix != self.in_place else os.path.splitext(os.path.basename(self.path))[0]+"."
        out_folder = (out_folder) if out_folder != self.in_place else os.path.dirname(self.path)
        split_list = []
        with open(self.path,"r") as file_obj:
            start_found = False
            count = 0
            out_obj = None
            for line in file_obj:
                if tracker: tracker.step(len(line))
                if (not start_found) and (not line.startswith(">")):
                    continue
                elif line.startswith(">"):
                    count+=1
                    if (not start_found) or ((count-1)%num_per_file)==0:
                        if not start_found: start_found = True
                        if out_obj: out_obj.close()
                        out_path = os.path.join(out_folder,file_prefix+str((count-1)//num_per_file)+file_suffix+file_extension)
                        out_obj = open(out_path,"w")
                        split_list.append(Handler(out_path,self.line_length))
                    out_obj.write(line)
                else:
                    out_obj.write(line)
        if out_obj: out_obj.close()
        if tracker: tracker.done()
        return split_list