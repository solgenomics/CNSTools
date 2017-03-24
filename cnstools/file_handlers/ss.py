import abstract_handler as ah
from collections import OrderedDict
from .._utils import MultiTracker
import os

class Entry(ah.Entry):
    def __init__(self, index, tuples_string, count):
        self.index, self.tuples_string, self.count = index, tuples_string, count

    def get_lines(self):
        return [str(self.index)+"\t"+self.tuples_string+"\t"+str(self.count)]

class KVinfo(ah.Entry):
    def __init__(self, kvdict):
        self._keys = []
        for key in kvdict:
            self._keys.append(key)
            setattr(self, key, kvdict[key])

    def get_lines(self):
        lines = []
        for key in self._keys:
            if key!="NAMES":
                lines.append(key+" = "+str(getattr(self,key)))
            else:
                lines.append(key+" = "+",".join(getattr(self,key)))
        lines.append("")
        return lines

class Handler(ah.Handler):
    in_place = object()

    def __init__(self, path):
        super(Handler, self).__init__(path)
        self._modify = self.modify
        self.modify = self._modify_wrap

    def _modify_wrap(self,*args,**kwargs):
        kwargs["kvinfo_first"]=True
        kwargs["add_end_break"]=True
        return self._modify(*args,**kwargs)

    def _entry_generator(self,kvinfo_first=False,parent=None,tracker_name=None):
        if tracker_name:
            size = os.stat(self.path).st_size
            if parent!=None:
                tracker = parent.subTracker(tracker_name,size,estimate=False,style="percent")
            else:
                tracker = MultiTracker(tracker_name,size,estimate=False,style="percent").auto_display(1)
        with open(self.path,"r") as file_object:

            kvinfo = OrderedDict()
            for line in file_object:
                if tracker_name:
                    tracker.step(len(line))
                line = line.strip()
                if line=="":
                    break
                elif kvinfo_first:
                    key,value = line.split("=")
                    key,value = key.strip(),value.strip()
                    if key=="NAMES":
                        value = value.split(",")
                    else:
                        try:
                            value = int(value)
                        except ValueError:
                            pass
                    kvinfo[key] = value
            if kvinfo_first:
                yield KVinfo(kvinfo)

            for line in file_object:
                if tracker_name:
                    tracker.step(len(line))
                line = line.strip()
                if line!="":
                    index, tuples_string, count = line.split("\t")
                    index, count = int(index), int(count)
                    yield Entry(index, tuples_string, count)

