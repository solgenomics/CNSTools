import abstract_handler as ah
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
    def to_bed6_entries(self,min_seg_length=7,min_seg_score=0.82,rejection_seg_len=12,rejection_seg_score=0.55):
        regions = [[0,len(self.val_list)]]
        bed_entries = []
        if len(self.val_list) >= rejection_seg_len:
            start = 0
            end = rejection_seg_len-1
            below_rejec_score = [val<0.55 for val in self.val_list[start:end]]
            #scans for regions of at least the rejection length in which all scores are lower that the rejection score
            #and cuts them out of the region
            while end<regions[-1][1]:
                below_rejec_score.append(self.val_list[end]<0.55)
                end+=1
                if all(below_rejec_score):
                    while end<regions[-1][1] and self.val_list[end] < rejection_seg_score:
                        end+=1
                    cut_region = regions.pop()
                    regions.append([cut_region[0],start])
                    regions.append([end,cut_region[1]])
                    start = regions[-1][0]
                    end = regions[-1][0]+rejection_seg_len-1
                    below_rejec_score = [val<0.55 for val in self.val_list[start:end]]
                    if regions[-1][1]-regions[-1][0] < min_seg_length:
                        break
                below_rejec_score.pop(0)
                start+=1
        for region in sorted(regions):
            region_sum = sum(val for val in self.val_list[region[0]:region[1]])
            while region_sum < min_seg_score*(region[1]-region[0]) and region[1]-region[0] >= min_seg_length:
                if self.val_list[region[1]-1]<self.val_list[region[0]]:
                    region_sum -= self.val_list[region[1]-1]
                    region[1] = region[1]-1
                else:
                    region_sum -= self.val_list[region[0]]
                    region[0] = region[0]+1
            if region_sum >= min_seg_score*(region[1]-region[0]) and region[1]-region[0] >= min_seg_length:
                bed_entries.append(bed6.Entry(self.chrom, self.start+region[0], self.start+region[1], name=str(region_sum/float(region[1]-region[0])), score=(region[1]-region[0]), strand="+"))
        return bed_entries


class Handler(ah.Handler):
    def _entry_generator(self,include_comments=False,tracker=None):
        with open(self.path,"r") as file_object:
            paragraph = []
            for line in file_object:
                if tracker: tracker.step(len(line))
                line = line.strip()
                if line.startswith("fixedStep"):
                    if len(paragraph)>0:
                        info_list = [item.split("=") for item in paragraph[0].split()]
                        step_type = info_list[0][0]
                        info_dict = {pair[0]:pair[1] for pair in info_list[1:]}
                        val_list = [float(item) for item in paragraph[1:]]
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
            if tracker: tracker.done()
    def to_bed(self,path,include_comments=False,tracker=None,**kwargs):
        def mod_func(entry):
            return entry.to_bed6_entries(**kwargs)
        return self.modify(mod_func,path=path,target_type=bed6.Handler,include_comments=include_comments,tracker=tracker)