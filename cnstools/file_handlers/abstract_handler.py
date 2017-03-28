import os
import itertools
from abc import ABCMeta, abstractmethod

class Entry(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_lines(self):
        pass


class Comment(Entry):
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_lines(self):
        pass


class Handler(object):
    __metaclass__ = ABCMeta

    def __init__(self, path):
        self.path = path

    def entries(self, **kwargs):
        return self._entry_generator(**kwargs)

    def modify(self, modify_func, path=None, target_type=None, add_end_break=False, **kwargs):
        if target_type == None:
            target_type = type(self)
        save_path = (path) if path != None else (self.path)
        temp_file_path = save_path+".fhmtemp"
        with open(temp_file_path,"w") as out_obj:
            first = True
            for entry in self._entry_generator(**kwargs):
                modified = modify_func(entry) if not isinstance(entry, Comment) else None
                if not isinstance(modified, (list,tuple)):
                    if modified == None:
                        modified = [entry]
                    elif modified==False:
                        modified = []
                    else:
                        modified = [modified]
                for result in modified:
                    for line in result.get_lines():
                        if first:
                            out_obj.write(line)
                            first = False
                        else:
                            out_obj.write("\n"+line)
            if add_end_break:
                out_obj.write("\n")
        os.rename(temp_file_path, save_path)
        return target_type(save_path)

    @abstractmethod
    def _entry_generator(self,include_comments=False):
        with open(self.path) as f:
            for line in f: 
                entry = None
                yield entry