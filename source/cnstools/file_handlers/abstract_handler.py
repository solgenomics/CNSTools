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
    in_place = object()

    def __init__(self, path):
        self.path = path
        self.size = os.stat(path).st_size

    def entries(self, **kwargs):
        return self._entry_generator(**kwargs)

    def modify(self, modfiy_func, path=in_place, target_type=in_place, **kwargs):
        if target_type == self.in_place:
            target_type = type(self)
        save_path = (path) if path != self.in_place else (self.path)
        temp_file_path = save_path+".fhmtemp"
        with open(temp_file_path,"w") as out_obj:
            first = True
            for entry in self._entry_generator(**kwargs):
                modified = modfiy_func(entry) if not isinstance(entry, Comment) else None
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
        os.rename(temp_file_path, save_path)
        return target_type(save_path)

    @abstractmethod
    def _entry_generator(self,include_comments=False,tracker=None):
        self.progress = 0
        with open(self.path) as f:
            for line in f: 
                self.progress += len(line)
                entry = None
                yield entry
        self.progress = self.size