import progress_tracker as pt

def main(bed_file,maf_file,out_file,min_size=1,max_gap_ratio=0.75,max_N_ratio=0.75):

    bed_entries = []
    with open(bed_file) as f:
        filelines = f.readlines()
        tracker = pt.Progress_tracker("Parsing .bed",len(filelines))
        tracker.display(estimate=False,rate=2)
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
        tracker.display()
        del tracker


    maf_entries = []
    header = []
    with open(maf_file) as f:
        body = [[]]
        filelines = f.readlines()
        tracker = pt.Progress_tracker("Parsing .maf",len(filelines))
        tracker.display(estimate=False,rate=2)
        if filelines[-1].strip()!="":
            filelines.append("")
        for line in filelines:
            stripped = line.strip()
            if stripped.startswith("#"):
                header.append(stripped)
            elif stripped=="":
                chunk = body[-1]
                if len(chunk)>1:
                    a_line = load_a_line(chunk[0]) if chunk[0].startswith("a") else None
                    s_lines = [load_s_line(line) for line in chunk if line.startswith("s")]
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
        tracker.display()
        del tracker
            
    s_line_valid = lambda x: x[3]/len(x[6]) > 1-max_gap_ratio \
                         and no_gap_len(x[6].replace('N',''))/x[3] > 1-max_N_ratio
    tracker = pt.Progress_tracker("Trimming .maf",len(bed_entries))
    tracker.display(estimate=False,rate=1)
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
        front_cut_num = gap_cut_loc(new_ma['s_lines'][0][6],front_offset)
        back_cut_num = gap_cut_loc(reversed(new_ma['s_lines'][0][6]),back_offset)
        for s_line in new_ma['s_lines']:
            s_line[2] = s_line[2]+no_gap_len(s_line[6][:front_cut_num])
            s_line[6] = s_line[6][front_cut_num:-back_cut_num] if back_cut_num!=0 else s_line[6][front_cut_num:]
            s_line[3] = no_gap_len(s_line[6])
        new_ma['s_lines'][:] = [line for line in new_ma['s_lines'] if s_line_valid(line)]
        if len(new_ma['s_lines'])>1 and max((line[3] for line in new_ma['s_lines'])) >= min_size:
            new_maf_entries.append(new_ma)
        tracker.step()
    tracker.display()
    del tracker

    with open(out_file,"w") as out:
        tracker = pt.Progress_tracker("Saving new .maf",len(new_maf_entries))
        tracker.display(estimate=False,rate=1)
        out_lines = []
        for entry in new_maf_entries:
            out.write(" ".join(entry['a_line'])+"\n") 
            for s_line in entry['s_lines']:
                out.write(" ".join([str(item) for item in s_line])+"\n")
            out.write("\n")
            tracker.step()
        tracker.display()
        del tracker

def no_gap_len(seq): return len(seq[:].replace("-",""))

def gap_cut_loc(seq,cutNum):
    anonymized = "".join(["x" if i!="-" else i for i in seq]) #use x to represent valid chars
    eliminated = anonymized.replace("x","o",cutNum) #change the first n x's to o's
    return eliminated.find("x") # finds the index of the next x, which is how many chars to cut off

def load_s_line(line):
    arr = [item for item in line.strip().split(" ") if item!=""]
    arr[2] = int(arr[2])
    arr[3] = int(arr[3])
    arr[5] = int(arr[5])
    return arr

def load_a_line(line):
    arr = [item for item in line.strip().split() if item!=""]
    return arr
        
def run(argv): 
    bed_file = argv[1]
    maf_file = argv[2]
    out_file = argv[3]
    main(bed_file,maf_file,out_file)