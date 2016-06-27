
import progress_tracker as pt

def main(mafFile,outFolder):

    fastaChunkBySpecies = {}
    #{'index': 8, 'score': 0.0, 's_lines': [['s', 'chr1', '1351', '15', '+', '52991155', 'TGAAAGGAAAATACC'], ['s', 'Cicer', '7292', '14', '+', '13526', 'TGAGATGAGACTAC-']]}
    with open(mafFile) as f:
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

        index = 0
        track = pt.Progress_tracker("Converting...",len(body),True)
        track.display(estimate=False, rate=5)
        for chunk in body:
            if len(chunk)>1:
                score = 0
                s_lines = [[item for item in line.split(" ") if item!=""] for line in chunk[1:] if line.startswith("s")]
                index+=1
                for line in s_lines:
                    newChunk = []
                    newChunk += [([">",line[1],"alignment#"+str(index),"score="+str(score)]+line[2:6])]
                    newChunk += formatSeq(line[6])
                    if newChunk[0][1] not in fastaChunkBySpecies: fastaChunkBySpecies[newChunk[0][1]] = []
                    fastaChunkBySpecies[newChunk[0][1]].append(newChunk)
            track.step()
        track.display()
    for i in fastaChunkBySpecies: 
        sp = fastaChunkBySpecies[i]
        with open(outFolder+i+".fasta","w") as out:
            for chunk in sp:
                out.write(chunk[0][0]+"|".join(chunk[0][1:])+"\n")
                out.write("\n".join(chunk[1:])+"\n")


def formatSeq(seq,lLen=70):
    nogaps = seq.replace("-","")
    lines = [nogaps[i:i+lLen] for i in xrange(0,len(nogaps),lLen)]
    return lines

def run(argv): 
    mafFile = argv[1]
    outFolder = argv[2] if outFolder.endswith("/") else argv[2]+"/"
    main(mafFile,outFolder)
