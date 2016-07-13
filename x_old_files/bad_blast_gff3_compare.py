
from progress_tracker import Progress_tracker
from interval_tree import Interval_tree

def find_ncs(blast_file,gff3_file,out_file):
    #import gff3
    genome_info = None
    print "loading "+gff3_file+"..."
    with open(gff3_file) as f:
        genome_info = f.readlines()
    chromosomes = {}
    tracker = Progress_tracker("making genome intervals",len(genome_info),True)
    tracker.display(False,1)
    for i in range(len(genome_info)-1,-1,-1):
        if genome_info[i].startswith("#"):
            del genome_info[i]
        else:
            genome_info[i] = genome_info[i].strip().split("\t")
            if not (genome_info[i][0] in chromosomes): 
                chromosomes[genome_info[i][0]] = [[int(genome_info[i][3]),int(genome_info[i][4]),genome_info[i]]]
            else: 
                chromosomes[genome_info[i][0]].append([int(genome_info[i][3]),int(genome_info[i][4]),genome_info[i]])
        tracker.step()
    tracker.display()
    del tracker
    trees = {}
    #for key in chromosomes: print key +" "+str(len(chromosomes[key]))
    #return
    tracker = Progress_tracker("building trees",len(genome_info),True)
    tracker.display(False,1)
    for key in chromosomes:
        trees[key] = Interval_tree()
        lower = float("inf")
        upper = 0
        for info_line in chromosomes[key]:
            for num in info_line[:2]:
                if num < lower: lower = num
                if num > upper: upper = num
        trees[key].build_from_list(chromosomes[key],lower,upper)
        tracker.step(len(chromosomes[key]))
    tracker.display()
    del tracker

    #import blast
    best_blast_info = [[0]]
    with open(blast_file) as f:
        lines = f.readlines()
        tracker = Progress_tracker("importing blast",len(lines),True)
        tracker.display(False,2)
        for line in lines:
            list = line.strip().split("\t")
            if list[0] != best_blast_info[-1][0]:
                best_blast_info.append(list)
            line = f.readline()
            tracker.step()
        best_blast_info.pop(0)
        tracker.display()
        del tracker

    out_lines = []
    for alignment in best_blast_info[:100]:
        out_lines.append("blast_result\t"+"\t".join(alignment))
        for line in trees[alignment[1]].range_query(int(alignment[8]),int(alignment[9])):
            out_lines.append("loc_match\t"+str("\t".join(line[2])))
    print "\n".join(out_lines)

def run(argv): 
    blast_file = argv[1] # "in.txt"
    gff3_file = argv[2] # "in.gff3"
    out_file = argv[3] # "out.txt"
    find_ncs(blast_file,gff3_file,out_file)