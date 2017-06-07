import os
import abstract_handler as ah
from .._utils import MultiTracker
import bed6

class Entry(ah.Entry):
    def __init__(self, step_type, chrom, start, step, val_list):
        self.step_type = step_type
        self.chrom = chrom
        self.start = start
        self.step = step
        self.val_list = val_list[:]
    def get_lines(self):
        info_line = "%s chrom=%s start=%s step=%s" % (self.step_type,self.chrom,self.start,self.step)
        lines = [info_line]+[str(val) for val in self.val_list]
        return lines
    def to_bed6_entries(self,min_seg_length=7,min_seg_score=0.82,max_conservation_gap=12,min_site_score=0.55):
        bed_entries = []
        if len(self.val_list) >= min_seg_length:
            #split into lists which are split by contiguous sections of size max_conservation_gap where every value is less than min_site_score
            nogap_lists = [[]]
            nogap_offsets = [0]
            consv_gap_count = 0
            for i in range(len(self.val_list)):
                nogap_lists[-1].append(self.val_list[i])
                if self.val_list[i]>=min_site_score:
                    consv_gap_count = 0
                else:
                    consv_gap_count += 1
                if consv_gap_count>=max_conservation_gap:
                    nogap_offsets.append(nogap_offsets[-1]+len(nogap_lists[-1]))
                    nogap_lists.append([])
                    consv_gap_count = 0

            #identify contiguous regions where every value is greater than or equal to min_seg_score
            regions = []
            for nogap_list,nogap_offset in zip(nogap_lists,nogap_offsets):
                i=0
                while i < len(nogap_list):
                    while i < len(nogap_list) and nogap_list[i]<min_seg_score:
                        i+=1
                    if not i < len(nogap_list):
                        continue
                    seg_start = i
                    while i < len(nogap_list) and nogap_list[i]>=min_seg_score:
                        i+=1
                    seg_end = i
                    regions.append([seg_start+nogap_offset,seg_end+nogap_offset])

            #combine sequencial regions such that the resulting regions have >=min_seg_score average value.
            regionsID = None
            currID = tuple(i for pair in regions for i in pair)
            direction = 1
            while regionsID!=currID:
                regions[:] = regions[::-1]
                direction*=-1
                for i in range(len(regions)-1,-1,-1):
                    if i>0:
                        tot_score = sum(self.val_list[regions[i-1][0]:regions[i][1]])/float(regions[i][1]-regions[i-1][0])
                        if tot_score >= min_seg_score:
                            regions[i-1:i+1] = ([regions[i-1][0],regions[i][1]],)
                            continue
                if direction==1:
                    regionsID = currID
                    currID = tuple(i for pair in regions for i in pair)

            #expand regions so that they cover the most contiguous values without dropping below an average value >=min_seg_score
            for reg in regions:
                done = False
                while not done:
                    expand_left =  (sum(self.val_list[reg[0]-1:reg[1]])/float(reg[1]-(reg[0]-1)), (reg[0]-1,reg[1]))
                    expand_right = (sum(self.val_list[reg[0]:reg[1]+1])/float(reg[1]+1-reg[0]), (reg[0],reg[1]+1))
                    best_expansion = max(expand_left,expand_right)
                    if (best_expansion[0]>=min_seg_score and
                        best_expansion[1][0] >= 0 and self.val_list[best_expansion[1][0]]>=min_site_score and
                        best_expansion[1][1] < len(self.val_list) and self.val_list[best_expansion[1][1]]>=min_site_score):
                        reg[:] = best_expansion[1]
                    else:
                        done = True
            # remove any regions too short to qualify
            regions[:] = (reg for reg in regions if reg[1]-reg[0]>=min_seg_length)
            #convert to bed
            for region in regions:
                region_sum = sum(self.val_list[region[0]:region[1]])
                bed_entries.append(bed6.Entry(self.chrom, self.start+region[0], self.start+region[1], name=str(region_sum/float(region[1]-region[0])), score=(region[1]-region[0]), strand="+"))
        return bed_entries


class Handler(ah.Handler):
    def _entry_generator(self,include_comments=False,genome=None,parent=None,tracker_name=None):
        if tracker_name:
            size = os.stat(self.path).st_size
            if parent!=None:
                tracker = parent.subTracker(tracker_name,size,estimate=False,style="percent")
            else:
                tracker = MultiTracker(tracker_name,size,estimate=False,style="percent").auto_display(1)
        with open(self.path,"r") as file_object:
            paragraph = []
            for line in file_object:
                if tracker_name:
                    tracker.step(len(line))
                line = line.strip()
                if line.startswith("fixedStep"):
                    if len(paragraph)>0:
                        info_list = [item.split("=") for item in paragraph[0].split()]
                        step_type = info_list[0][0]
                        info_dict = {pair[0]:pair[1] for pair in info_list[1:]}
                        val_list = [float(item) for item in paragraph[1:]]
                        if genome:
                            info_dict['chrom'] = genome+":"+info_dict['chrom']
                        yield Entry(step_type,info_dict['chrom'],int(info_dict['start']),int(info_dict['step']),val_list)
                    paragraph = [line]
                else:
                    paragraph.append(line)
            if len(paragraph)>0:
                info_list = [item.split("=") for item in paragraph[0].split()]
                step_type = info_list[0][0]
                info_dict = {pair[0]:pair[1] for pair in info_list[1:]}
                val_list = [float(item) for item in paragraph[1:]]
                yield Entry(step_type,info_dict['chrom'],int(info_dict['start']),int(info_dict['step']),val_list)
            if tracker_name:
                tracker.done()
    def to_bed(self,path,include_comments=False,genome=None,parent=None,tracker_name=None,**kwargs):
        def mod_func(entry):
            return entry.to_bed6_entries(**kwargs)
        return self.modify(mod_func,path=path,target_type=bed6.Handler,include_comments=include_comments,genome=genome,parent=parent,tracker_name=tracker_name)