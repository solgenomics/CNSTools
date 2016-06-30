
import progress_tracker as pt

def main(blast_file,outFile):
    #import blast
    best_blast_info = [[None]]
    with open(blast_file) as f:
        lines = f.readlines()
        tracker = pt.Progress_tracker("Loading BLAST results",len(lines),True).display(estimate=False,rate=1)
        for line in lines:
            list = line.strip().split("\t")
            if list[0] != best_blast_info[-1][0]:
                best_blast_info.append(list)
            tracker.step()
        best_blast_info.pop(0)
        tracker.display()
        del tracker

    with open(outFile,"w") as out:
        bed_lines = []
        tracker = pt.Progress_tracker("Saving .bed",len(best_blast_info),True)
        tracker.display(estimate=False,rate=1)
        for line in best_blast_info:
            tracker.step()
            bed_line = [0]*6
            bed_line[0] = line[1]
            if int(line[8])<=int(line[9]):
                bed_line[1],bed_line[2],bed_line[5] = line[8],line[9],"+"
            else:
                bed_line[1],bed_line[2],bed_line[5] = line[9],line[8],"-"
            bed_line[3] = line[0]
            bed_line[4] = line[11]
            out.write("\t".join(bed_line)+"\n")
        tracker.display()
        del tracker

def run(argv):
    blast_file = argv[1]
    outFile = argv[2]
    main(blast_file,outFile)