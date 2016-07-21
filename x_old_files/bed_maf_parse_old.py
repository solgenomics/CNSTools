"""A script for creating an intersection between a .maf multiple alignment file and a .bed file containing ranges within those alignments.
"""

from _utils import Progress_tracker
def main(bed_file,maf_file,out_file,min_size=15,max_gap_ratio=0.5,max_N_ratio=0.5):
    """This function runs the main workflow for the script.

    :param str bed_file: Path of the .bed file to be used.
    :param str maf_file: Path of the .maf file to be used.
    :param str out_file: Path to save the .maf format output to.
    :param int min_size: Minimum size (including gaps) for CNS returned in the resulting .maf file.
    :param float max_gap_ratio: Maximum ratio of gaps to nucleotides that may be present in a returned CNS.
    :param float max_gap_ratio: Maximum ratio of N positions to known nucleotides that may be present in a returned CNS.
    :returns:  `None`
    """

    bed_entries = []
    with open(bed_file) as f:
        filelines = f.readlines()
        tracker = Progress_tracker("Parsing .bed",len(filelines)).estimate(False).auto_display()
        for line in filelines:
            list = line.strip().split("\t")
            list[3] = list[3].split("ID=")[1]
            dict = {
                    "seq":list[0],
                    "start":int(list[1]),
                    "stop":int(list[2]),
                    "parent":int(list[3]),
                    "list":list,
                    }
            if(dict['stop']-dict['start']>=min_size):
                bed_entries.append(dict)
            tracker.step()
        tracker.done()


    maf_entries = []
    header = []
    with open(maf_file) as f:
        body = [[]]
        filelines = f.readlines()
        tracker = Progress_tracker("Parsing .maf",len(filelines)).estimate(False).auto_display()
        if filelines[-1].strip()!="":
            filelines.append("")
        for line in filelines:
            stripped = line.strip()
            if stripped.startswith("#"):
                header.append(stripped)
            elif stripped=="":
                chunk = body[-1]
                if len(chunk)>1:
                    a_line = _load_a_line(chunk[0]) if chunk[0].startswith("a") else None
                    s_lines = [_load_s_line(line) for line in chunk if line.startswith("s")]
                    maf_entries.append({
                        "a_line":a_line,
                        "s_lines":s_lines,
                        "main_start":s_lines[0][2],
                        "main_stop":s_lines[0][2]+s_lines[0][3],
                        "main_seq":s_lines[0][1]
                        })
                body.append([])
            else:
                body[-1].append(stripped)
            tracker.step()
        tracker.done()
            
    tracker = Progress_tracker("Trimming .maf",len(bed_entries)).estimate(False).auto_display()
    index = 0
    new_maf_entries = []
    for bed in bed_entries:
        parent = maf_entries[bed['parent']]
        new_ma = {
            "a_line":parent['a_line'][:],
            "s_lines":[s_line[:6]+[s_line[6][:]] for s_line in parent['s_lines']]
        }
        new_ma['a_line'].append("cns_ID="+str(index))
        new_ma['a_line'].append("original_maf_ID="+str(bed['parent']))
        index+=1
        front_offset = bed['start']-parent['main_start']
        back_offset = parent['main_stop']-bed['stop']
        front_cut_num = _gap_cut_loc(new_ma['s_lines'][0][6],front_offset)
        back_cut_num = _gap_cut_loc(reversed(new_ma['s_lines'][0][6]),back_offset)
        for s_line in new_ma['s_lines']:
            s_line[2] = s_line[2]+_no_gap_len(s_line[6][:front_cut_num])
            s_line[6] = s_line[6][front_cut_num:-back_cut_num] if back_cut_num!=0 else s_line[6][front_cut_num:]
            s_line[3] = _no_gap_len(s_line[6])
        new_ma['s_lines'][:] = [line for line in new_ma['s_lines'] if _s_line_valid(line,max_gap_ratio,max_N_ratio)]
        if len(new_ma['s_lines'])>1 and max((line[3] for line in new_ma['s_lines'])) >= min_size:
            new_maf_entries.append(new_ma)
        tracker.step()
    tracker.done()

    with open(out_file,"w") as out:
        tracker = Progress_tracker("Saving new .maf",len(new_maf_entries)).estimate(False).auto_display()
        out_lines = []
        for entry in new_maf_entries:
            out.write(" ".join(entry['a_line'])+"\n") 
            for s_line in entry['s_lines']:
                out.write(" ".join([str(item) for item in s_line])+"\n")
            out.write("\n")
            tracker.step()
        tracker.done()

def _s_line_valid(x,max_gap_ratio,max_N_ratio):
    if x[3]/float(len(x[6])) > 1-max_gap_ratio and _no_gap_len(x[6][:].replace('N',''))/float(x[3]) > 1-max_N_ratio:
        return True
    else : return False


def _no_gap_len(seq): return len(seq[:].replace("-",""))

def _gap_cut_loc(seq,cutNum):
    anonymized = "".join(["x" if i!="-" else i for i in seq]) #use x to represent valid chars
    eliminated = anonymized.replace("x","o",cutNum) #change the first n x's to o's
    return eliminated.find("x") # finds the index of the next x, which is how many chars to cut off

def _load_s_line(line):
    arr = [item for item in line.strip().split(" ") if item!=""]
    arr[2] = int(arr[2])
    arr[3] = int(arr[3])
    arr[5] = int(arr[5])
    return arr

def _load_a_line(line):
    arr = [item for item in line.strip().split() if item!=""]
    return arr
        
def file_run(bed_file,maf_file,out_file,min_size=15,max_gap_ratio=0.5,max_N_ratio=0.5):
    """This function parses the arguements it recieves then calls :func:`.main` using them.
    
    :param str bed_file: Path of the .bed file to be used.
    :param str maf_file: Path of the .maf file to be used.
    :param str out_file: Path to save the .maf format output to.
    :param int min_size: Minimum size (including gaps) for CNS returned in the resulting .maf file.
    :param float max_gap_ratio: Maximum ratio of gaps to nucleotides that may be present in a returned CNS.
    :param float max_gap_ratio: Maximum ratio of N positions to known nucleotides that may be present in a returned CNS.
    :returns:  `None`
    """
    min_size = int(min_size)
    max_gap_ratio = float(max_gap_ratio)
    max_N_ratio = float(max_N_ratio)
    main(bed_file,maf_file,out_file,min_size,max_gap_ratio,max_N_ratio)