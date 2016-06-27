
import progress_tracker as pt

def main(blast_file,outFile):
    #import blast
    best_blast_info = [[None]]
    with open(blast_file) as f:
        lines = f.readlines()
        tracker = pt.Progress_tracker("importing blast",len(lines),True)
        tracker.display(False,2)
        for line in lines:
            list = line.strip().split("\t")
            if list[0] != best_blast_info[-1][0]:
                best_blast_info.append(list)
            tracker.step()
        best_blast_info.pop(0)
        tracker.display()
        del tracker

    print best_blast_info[:10]

    bed_lines = []
    tracker = pt.Progress_tracker("formatting as bed",len(best_blast_info),True)
    tracker.display(estimate=False,rate=1)
    for line in best_blast_info:
        tracker.step()
        bed_line = [0]*6
        bed_line[0] = line[1]
        if line[8]<=line[9]:
            bed_line[1],bed_line[2],bed_line[5] = line[8],line[9],"+"
        else:
            bed_line[1],bed_line[2],bed_line[5] = line[9],line[8],"-"
        bed_line[3] = line[0]
        bed_line[4] = line[11]
        bed_lines.append(bed_line)
    tracker.display()
    del tracker
    print bed_lines[:10]
    with open(outFile,"w") as out:
        out.write("\n".join(("\t".join(line) for line in bed_lines)))
    print "Done."


def run(argv):
    blast_file = argv[1]
    outFile = argv[2]
    main(blast_file,outFile)