
import progress_tracker as pt

def main(gff3_file,sequence_types,out_file):
    print sequence_types
    entries = []
    with open(gff3_file) as f:
        for line in f.readlines():
            if not line.startswith("#"):
                list = line.strip().split('\t')
                if list[2] in sequence_types:
                    entries.append(list)
    print entries[:10]

    bed_lines = []
    for line in entries:
        bed_line = [0]*6
        bed_line[0] = line[0]
        bed_line[1] = str(int(line[3])-1)
        bed_line[2] = line[4]
        bed_line[3] = "seqType="+line[2]+","+line[8]
        bed_line[4] = get_score(line[5])
        bed_line[5] = line[6]
        bed_lines.append(bed_line)

    print bed_lines[:10]

    with open(out_file,"w") as out:
        out.write("\n".join(["\t".join(line) for line in bed_lines]))

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