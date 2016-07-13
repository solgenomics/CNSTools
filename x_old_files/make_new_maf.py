
import progress_tracker as pt

def main(rangeListFile,originalMaf,outFile):

    newRangeList = None

    print "Loading and parsing "+rangeListFile+"..."
    with open(rangeListFile) as f:
        newRangeList = [[int(item) for item in line.strip().split('\t')] for line in f.readlines()]
        newRangeList.sort()
    print "Done"

    alignments = []
    header = None
    print "Loading and parsing "+originalMaf+"..."
    with open(originalMaf) as f:
        header = []
        body = [[]]

        for line in f.readlines():
            stripped = line.strip()
            if stripped.startswith("#"):
                header.append(stripped)
            elif stripped=="":
                body.append([])
            else:
                body[-1].append(stripped)

        for chunk in body:
            if len(chunk)>1:
                score = float(chunk[0].split("=")[1].strip()) if chunk[0].startswith("a") else 0
                s_lines = [load_s_line(line) for line in chunk[1:] if line.startswith("s")]
                alignments.append({"score":score,"s_lines":s_lines})
    print "Done"

    #make and print newAlignments
    newAlignments = []
    with open(outFile,"w") as out:
        out.write("\n".join(header)+"\n\n")
        displayLoop = 0
        tracker = pt.Progress_tracker("Slicing Alignments",len(newRangeList),True)
        tracker.display(estimate=False)
        for newRange in newRangeList:
            parentSeq = alignments[newRange[2]]
            newAlignments.append({"score":0,"s_lines":[]})
            frontOffsetUngapped = newRange[0]-parentSeq["s_lines"][0][2]
            backOffsetUngapped = (parentSeq["s_lines"][0][2]+parentSeq["s_lines"][0][3])-newRange[1]
            frontOffset = gapCutLoc(parentSeq["s_lines"][0][6],frontOffsetUngapped)
            backOffset = gapCutLoc(reversed(parentSeq["s_lines"][0][6]),backOffsetUngapped)
            for s_line in parentSeq["s_lines"]:
                newAlignments[-1]["s_lines"].append(s_line[:2]+[0]*5)

                newAlignments[-1]["s_lines"][-1][6] = (s_line[6][frontOffset:len(s_line[6])-backOffset+1])
                newAlignments[-1]["s_lines"][-1][3] = noGapLen(newAlignments[-1]["s_lines"][-1][6])

                newAlignments[-1]["s_lines"][-1][2] = (s_line[2]+noGapLen(s_line[6][:frontOffset]))
                newAlignments[-1]["s_lines"][-1][4] = (s_line[4])
                newAlignments[-1]["s_lines"][-1][5] = (s_line[5])
            if isGoodAlignment(newAlignments[-1]):
                out.write("a score="+str(newAlignments[-1]["score"])+"\n")
                for s_line in newAlignments[-1]["s_lines"]:
                    out.write(" ".join([str(item) for item in s_line])+"\n")
                out.write("\n")
            tracker.step()
            displayLoop+=1
            if(displayLoop>=10000):
                displayLoop = 0
                tracker.display()
        tracker.display()

def load_s_line(line):
    arr = [item for item in line.strip().split(" ") if item!=""]
    arr[2] = int(arr[2])
    arr[3] = int(arr[3])
    arr[5] = int(arr[5])
    return arr

def isGoodAlignment(alignment):
    return noGapLen(alignment["s_lines"][-1][-1]) > 10

def gapCutLoc(seq,cutNum):
    anonymized = "".join(["x" if i!="-" else i for i in seq]) #use x to represent valid chars
    eliminated = anonymized.replace("x","o",cutNum) #change the first n x's to o's
    return eliminated.find("x") # finds the index of the next x, which is how many chars to cut off

def noGapLen(seq): return len(seq[::1].replace("-",""))

def run(argv):
    rangeListFile = argv[1]#"outtest.txt"
    originalMaf = argv[2]#"roast2.maf"
    outFile = argv[3]#"outFileTest.maf"
    main(rangeListFile,originalMaf,outFile)
