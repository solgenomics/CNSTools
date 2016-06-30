
import progress_tracker as pt

def main(mafFile,outFile):
    chunks = [[]]
    with open(mafFile) as f:
        lines = f.readlines()
        tracker = pt.Progress_tracker("Loading .maf",len(lines),True)
        tracker.display(estimate=False,rate=1)
        for line in lines:
            tracker.step()
            if not line.startswith("#"):
                stripped = line.strip()
                if stripped!="":
                    chunks[-1].append(stripped.split())
                else:
                    chunks.append([])
        tracker.display()
        del tracker

    master_alignment_lines = []
    tracker = pt.Progress_tracker("Parsing .maf",len(chunks),True)
    tracker.display(estimate=False,rate=1)
    for chunk in chunks:
        tracker.step()
        for line in chunk:
            if line[0]=="s":
                master_alignment_lines.append(line)
                break
    tracker.display()
    del tracker

    bed_lines = []
    with open(outFile,"w") as out:
        tracker = pt.Progress_tracker("Exporting .bed",len(master_alignment_lines),True)
        tracker.display(estimate=False,rate=1)
        for i in xrange(len(master_alignment_lines)):
            tracker.step()
            line = master_alignment_lines[i]
            bed_line = [0]*6
            bed_line[0] = line[1]
            bed_line[1] = line[2]
            bed_line[2] = str(int(line[2])+int(line[3]))
            bed_line[3] = "ID="+str(i)
            bed_line[4] = "0"
            bed_line[5] = line[4]
            out.write("\t".join(bed_line)+"\n")
        tracker.display()
        del tracker
    print "Done."

def run(argv): 
    mafFile = argv[1]
    outFile = argv[2]
    main(mafFile,outFile)