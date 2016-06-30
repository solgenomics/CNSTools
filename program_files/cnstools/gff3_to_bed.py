
import progress_tracker as pt

def main(gff3_file,sequence_types,out_file):
    entries = []
    with open(gff3_file) as f:
        lines = f.readlines()
        tracker = pt.Progress_tracker("Loading .gff",len(lines),True)
        tracker.display(estimate=False,rate=1)
        for line in lines:
            if not line.startswith("#"):
                list = line.strip().split('\t')
                if list[2] in sequence_types:
                    entries.append(list)
            tracker.step()
        tracker.display()
        del tracker

    with open(out_file,"w") as out:
        tracker = pt.Progress_tracker("Saving .bed",len(entries),True)
        tracker.display(estimate=False,rate=1)
        for line in entries:
            bed_line = [0]*6
            bed_line[0] = line[0]
            bed_line[1] = str(int(line[3])-1)
            bed_line[2] = line[4]
            bed_line[3] = "seqType="+line[2]+";"+line[8]
            bed_line[4] = get_score(line[5])
            bed_line[5] = line[6]
            out.write("\t".join(bed_line)+"\n")
            tracker.step()
        tracker.display()
        del tracker

def get_score(score):
    try:
        int(score)
        return score
    except ValueError:
        return "0"

def run(argv): 
    gff3_file = argv[1]
    sequence_types = argv[2].split(",")
    out_file = argv[3]
    main(gff3_file,sequence_types,out_file)