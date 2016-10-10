import os
import itertools
from abc import ABCMeta, abstractmethod

class File_entry(object):
    """docstring for File_entry"""
    __metaclass__ = ABCMeta
    @abstractmethod
    def get_lines(self):
        pass
class File_handler(object):
    """docstring for File_handler"""
    __metaclass__ = ABCMeta
    in_place = object()
    def __init__(self, path):
        self.path = path
    def entries(self):
        return self._entry_generator()
    def modify_entries(self, modfiy_func, path=in_place, target_type=in_place):
        if target_type == self.in_place:
            target_type = type(self)
        save_path = (path) if path != self.in_place else (self.path)
        temp_file_path = save_path+".fhmtemp"
        with open(temp_file_path,"w") as out_obj:
            first = True
            for entry in self.entries():
                results = modfiy_func(entry)
                if not isinstance(results, (list,tuple)):
                    if results == None:
                        results = [entry]
                    elif results==False:
                        results = []
                    else:
                        results = [results]
                for result in results:
                    for line in result.get_lines():
                        if first:
                            out_obj.write(line)
                            first = False
                        else:
                            out_obj.write("\n"+line)
        os.rename(temp_file_path, save_path)
        return target_type(save_path)

    @abstractmethod
    def _entry_generator(self,file_path):
        pass

        
class Wiggle_entry(File_entry):
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
class Wiggle_handler(File_handler):
    """docstring for Test_handler"""
    def _entry_generator(self):
        with open(self.path,"r") as file_object:
            paragraph = []
            for line in file_object:
                line = line.strip()
                if line.startswith("fixedStep"):
                    if len(paragraph)>0:
                        info_list = [item.split("=") for item in paragraph[0].split()]
                        step_type = info_list[0][0]
                        info_dict = {pair[0]:pair[1] for pair in info_list[1:]}
                        val_list = [float(item) for item in paragraph[1:]]
                        yield Wiggle_entry(step_type,info_dict['chrom'],int(info_dict['start']),int(info_dict['step']),val_list)
                    paragraph = [line]
                else:
                    paragraph.append(line)
            if len(paragraph)>0:
                info_list = [item.split("=") for item in paragraph[0].split()]
                step_type = info_list[0][0]
                info_dict = {pair[0]:pair[1] for pair in info_list[1:]}
                val_list = [float(item) for item in paragraph[1:]]
                yield Wiggle_entry(step_type,info_dict['chrom'],int(info_dict['start']),int(info_dict['step']),val_list)
    def to_bed(self,path,min_seg_length=7,min_seg_score=0.82,rejection_seg_len=12,rejection_seg_score=0.55):
        def mod_func(entry):
            return self._wig_to_bed(entry,min_seg_length,min_seg_score,rejection_seg_len,rejection_seg_score)
        return self.modify_entries(mod_func,path,Bed6_handler)
    @staticmethod
    def _wig_to_bed(entry,min_seg_length,min_seg_score,rejection_seg_len,rejection_seg_score):
        regions = [[0,len(entry.val_list)]]
        bed_entries = []
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
        for region in sorted(regions):
            region_sum = sum(val for val in entry.val_list[region[0]:region[1]])
            while region_sum < min_seg_score*(region[1]-region[0]) and region[1]-region[0] >= min_seg_length:
                if entry.val_list[region[1]-1]<entry.val_list[region[0]]:
                    region_sum -= entry.val_list[region[1]-1]
                    region[1] = region[1]-1
                else:
                    region_sum -= entry.val_list[region[0]]
                    region[0] = region[0]+1
            if region_sum >= min_seg_score*(region[1]-region[0]) and region[1]-region[0] >= min_seg_length:
                bed_entries.append(Bed6_entry(entry.chrom, entry.start+region[0], entry.start+region[1], name=str(region_sum/float(region[1]-region[0])), score=(region[1]-region[0]), strand="+"))
        return bed_entries

class Bed6_entry(File_entry):
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
    def get_lines(self):
        strs = (str(i) if i!=None else '.' for i in (self.chrom, self.chromStart, self.chromEnd, self.name, self.score, self.strand))
        return ["\t".join(strs)]
class Bed6_handler(File_handler):
    def _entry_generator(self):
        with open(self.path,"r") as file_object:
            for line in file_object:
                fields = line.strip().split('\t')
                if len(fields)>1:
                    if len(fields)<6:
                        fields.append([None]*(6-len(fields)))
                    fields[:] = [item if item!='.' else None for item in fields]
                    yield Bed6_entry(*fields)

class Bed13_entry(File_entry):
    def __init__(self,*fields):
        if len(fields)!=13:
            raise Exception("not 13 fields")
        self.first = Bed6_entry(*(fields[0:6]))
        self.second = Bed6_entry(*(fields[6:12]))
        self.score = int(fields[12])
    def get_lines(self):
        return ["\t".join([self.first.get_lines()[0],self.second.get_lines()[0],str(self.score)])]
class Bed13_handler(File_handler):
    """docstring for Bed13_handler"""
    def _entry_generator(self):
        with open(self.path,"r") as file_object:
            for line in file_object:
                fields = line.strip().split('\t')
                if len(fields)>1:
                    if len(fields)<13:
                        fields.append([None]*(13-len(fields)))
                    fields[:] = [item if item!='.' else None for item in fields]
                    yield Bed13_entry(*fields)

class Fasta_entry(File_entry):
    """docstring for Test_entry"""
    def __init__(self, description, data, line_length):
        self.description,self.data,self.line_length = description,data,line_length
    def get_lines(self):
        lines = [">"+self.description]
        lines+= [self.data[i:i+self.line_length] for i in range(0,len(self.data),self.line_length)]
        return lines

class Fasta_handler(File_handler):
    in_place = object()
    accepted_symbols = set(list("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ*-"))
    """docstring for Test_handler"""
    def __init__(self, path,line_length=80):
        self.line_length = line_length
        super(Fasta_handler, self).__init__(path)
    def _entry_generator(self):
        with open(self.path,"r") as file_object:
            entry_description = None
            entry_data = ""
            for line in file_object:
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
    def split(self,num_per_file=1,out_folder=in_place,file_prefix=in_place,file_suffix="",file_extension=".fa"):
        file_prefix = (file_prefix) if file_prefix != self.in_place else os.path.splitext(os.path.basename(self.path))[0]+"."
        out_folder = (out_folder) if out_folder != self.in_place else os.path.dirname(self.path)
        split_list = []
        with open(self.path,"r") as file_obj:
            start_found = False
            count = 0
            out_obj = None
            for line in file_obj:
                if (not start_found) and (not line.startswith(">")):
                    continue
                elif line.startswith(">"):
                    count+=1
                    if (not start_found) or ((count-1)%num_per_file)==0:
                        if not start_found: start_found = True
                        if out_obj: out_obj.close()
                        out_path = os.path.join(out_folder,file_prefix+str((count-1)//num_per_file)+file_suffix+file_extension)
                        out_obj = open(out_path,"w")
                        split_list.append(Fasta_handler(out_path,self.line_length))
                    out_obj.write(line)
                else:
                    out_obj.write(line)
        if out_obj: out_obj.close()
                
        return split_list







class Test_entry(File_entry):
    """docstring for Test_entry"""
    def __init__(self, a, b):
        self.a,self.b = a,b
    def get_lines(self):
        return [("%s\t%s"%(self.a,self.b))]
class Test_handler(File_handler):
    """docstring for Test_handler"""
    def _entry_generator(self):
        with open(self.path,"r") as file_object:
            for line in file_object:
                a,b = line.strip().split("\t")
                yield Test_entry(a,b)
        