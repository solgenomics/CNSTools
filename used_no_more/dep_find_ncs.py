
import progress_tracker as pt

def find_ncs(chromosomeName,alignmentFile,codingFile,outFile):

    # Scrapes alignmentFile for alignment starts and stops
    unchecked = []
    print "Loading and parsing "+alignmentFile+"..."
    with open(alignmentFile) as f:
        done = False
        while not done:
            curr = f.readline()
            while not curr.startswith('a'):
                curr = f.readline()
                if not curr:
                    done = True
                    break
            if done: break
            else: #found an Alignment!
                infoLine = f.readline()
                info = [item for item in infoLine.split(" ") if item != ""]
                start = int(info[2])
                length = int(info[3])
                unchecked.append([start,start+length-1,len(unchecked)])# list of start (int), end (int), index of the original alignment in the file (int)
    print "Done"

    print "Loading and parsing "+codingFile+"..."
    #Scrapes Coding seq regions from codingFile
    codingRegions = []
    with open(codingFile) as f:
        for line in f:
            if line.startswith(chromosomeName):
                info = [item for item in line.split("\t") if item != ""]
                if len(info)<3: print info
                if(info[2]=="CDS"):
                    codingRegions.append((int(info[3]),int(info[4])))
        codingRegions.sort()
    print "Done"

    with open(outFile,"w") as out:
        nonCoding = []
        tracker = pt.Progress_tracker("Eliminating coding regions",len(unchecked),True)
        displayLoop = 0
        loopSize = len(unchecked)/100-1
        tracker.display(estimate=False)
        for alignment in unchecked:
            # check if each alignment is in a coding region
            for r in getNoCollisionRanges(alignment,codingRegions):
                out.write("\t".join([str(i) for i in r])+"\n")
                nonCoding.append(r)
            tracker.step()
            displayLoop+=1
            if(displayLoop>=loopSize):
                displayLoop = 0
                tracker.display()
        tracker.display()

def getNoCollisionRanges(toCheck,regions,start=0):
    newRange = [val for val in toCheck] if start==0 else toCheck # duplicates the checked region but does not duplicate it during recursion
    for i in range(start,len(regions)):
        if(newRange[0]<regions[i][0] and newRange[1]<regions[i][0]):
            continue
        elif(newRange[0]>regions[i][1] and newRange[1]>regions[i][1]):
            continue
        elif(newRange[0]>=regions[i][0] and newRange[1]<=regions[i][1]):
            return []
        elif(newRange[0]>=regions[i][0] and newRange[1]>regions[i][1]):
            newRange[0] = regions[i][1]+1
        elif(newRange[0]<regions[i][0] and newRange[1]<=regions[i][1]):
            newRange[1] = regions[i][0]-1
        elif(newRange[0]<regions[i][0] and newRange[1]>regions[i][1]):
            return getNoCollisionRanges([newRange[0],regions[i][0]-1,newRange[2]],regions,start=i+1)+getNoCollisionRanges([regions[i][1]+1,newRange[1],newRange[2]],regions,start=i+1)
        else: print "Wut. (Check for non-numbers in range data)"
    return [newRange]

def run(argv):
    chromosomeName = argv[1] # "chr1"
    alignmentFile = argv[2] # "roast2.maf"
    codingFile = argv[3] # "Mt4.0v1_genome.gff3"
    outFile = argv[4] #"out.txt"
    find_ncs(chromosomeName,alignmentFile,codingFile,outFile)
