
from _utils import Progress_tracker
def main(mafFile,outFolder):

    fastaChunkBySpecies = {}
    #{'index': 8, 'score': 0.0, 's_lines': [['s', 'chr1', '1351', '15', '+', '52991155', 'TGAAAGGAAAATACC'], ['s', 'Cicer', '7292', '14', '+', '13526', 'TGAGATGAGACTAC-']]}
    with open(mafFile) as f:
        header = []
        body = [[]]
        lines = f.readlines()
        track = Progress_tracker("Loading .maf",len(lines),True).display(estimate=False, rate=1)
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("#"):
                header.append(stripped)
            elif stripped=="":
                body.append([])
            else:
                body[-1].append(stripped)
            track.step()
        track.display()
        del track

        index = 0
        track = Progress_tracker("Converting to .fasta",len(body),True).display(estimate=False, rate=1)
        for chunk in body:
            chunkInfo = []
            if len(chunk)>1:
                for line in chunk:
                    if line.startswith('a'):
                        chunkInfo = line.split()[1:]
                    if line.startswith('s'):
                        line = line.split()
                        newChunk = []
                        newChunk += [([">",line[1]]+chunkInfo)]
                        newChunk += formatSeq(line[6])
                        if newChunk[0][1] not in fastaChunkBySpecies: fastaChunkBySpecies[newChunk[0][1]] = []
                        fastaChunkBySpecies[newChunk[0][1]].append(newChunk)   
            track.step()
        track.display()
        del track

    for i in fastaChunkBySpecies: 
        sp = fastaChunkBySpecies[i]
        track = Progress_tracker("Saving %s.fasta"%(i),len(fastaChunkBySpecies[i]),True).display(estimate=False, rate=1)
        with open(outFolder+i+".fasta","w") as out:
            for chunk in sp:
                out.write(chunk[0][0]+"|".join(chunk[0][1:])+"\n")
                out.write("\n".join(chunk[1:])+"\n")
                track.step()
        track.display()
        del track


def formatSeq(seq,lLen=70):
    nogaps = seq.replace("-","")
    lines = [nogaps[i:i+lLen] for i in xrange(0,len(nogaps),lLen)]
    return lines

def file_run(mafFile,outFolder):
    if not outFolder.endswith("/"):
        outFolder+="/"
    main(mafFile,outFolder)
