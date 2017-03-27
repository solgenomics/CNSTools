import abstract_handler as ah
from .._utils import MultiTracker
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

    def _entry_generator(self,parent=None,tracker_name=None):
        if tracker_name:
            size = os.stat(self.path).st_size
            if parent!=None:
                tracker = parent.subTracker(tracker_name,size,estimate=False,style="percent")
            else:
                tracker = MultiTracker(tracker_name,size,estimate=False,style="percent").auto_display(1)
        with open(self.path,"r") as file_object:
            entry_description = None
            entry_data = ""
            for line in file_object:
                if tracker_name: tracker.step(len(line))
                line = line.strip()
                if line.startswith(">"):
                    if entry_description!=None:
                        yield Entry(entry_description,entry_data,self.line_length)
                    entry_data = ""
                    entry_description = line[1:]
                    continue
                elif entry_description!=None:
                    entry_data+="".join([symbol for symbol in line if symbol in self.accepted_symbols])
            yield Entry(entry_description,entry_data)
            if tracker_name: tracker.done()

    def split(self,num_per_file=1,out_folder=in_place,file_prefix=in_place,file_suffix="",file_extension=".fa",parent=None,tracker_name=None):
        file_prefix = (file_prefix) if file_prefix != self.in_place else os.path.splitext(os.path.basename(self.path))[0]+"."
        out_folder = (out_folder) if out_folder != self.in_place else os.path.dirname(self.path)
        split_list = []
        size = os.stat(self.path).st_size
        if tracker_name:
            size = os.stat(self.path).st_size
            if parent!=None:
                tracker = parent.subTracker(tracker_name,size,estimate=False,style="percent")
            else:
                tracker = MultiTracker(tracker_name,size,estimate=False,style="percent").auto_display(1)
        with open(self.path,"r") as file_obj:
            start_found = False
            count = 0
            out_obj = None
            for line in file_obj:
                if tracker_name: tracker.step(len(line))
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
        if tracker_name: tracker.done()
        return split_list