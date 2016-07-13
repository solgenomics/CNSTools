import progress_tracker as pt

def main(data,out_folder):
    maf_entries_by_cns_id = {}
    with open(data['cns_maf']) as maf_file:
        lines = [line.strip() for line in maf_file.readlines() if not line.startswith("#")]
        new = True
        line_list = []
        tracker = pt.Progress_tracker("Parsing .maf",len(lines)).display(estimate=False,rate=1)
        for line in lines:
            if line =="" and len(line_list)>0:
                dict = {
                "a_line":None,
                "seqs":[]
                }
                for line_2 in line_list:
                    if line_2.startswith('a'):
                        dict['a_line'] = line_2.split()
                    if line_2.startswith('s'):
                        dict['seqs'].append({'s_line':line_2.split()})
                key = int(dict['a_line'][2].split('=')[1])
                maf_entries_by_cns_id[key] = dict
                line_list = []
            else:
                line_list.append(line)
            tracker.step()
        tracker.display()
        del tracker
        
    assoc_files = [(seq['maf_name'],open(seq['cns_assoc'])) for seq in data['seqs'] if 'cns_assoc' in seq]
    tracker = pt.Progress_tracker("Parsing association files",len(assoc_files)).display(estimate=False,rate=1)
    for name,file in assoc_files:
        lines = file.readlines()
        for line in lines:
            list = line.strip().split('\t')
            info = list[3].split("|")
            seq_key = int(info[2].split("=")[1])
            gen_key = info[0]
            closest_info = list[9]
            closest_dist = int(list[-1])
            CNS_type = "intergenic"
            if(abs(closest_dist)<=1000):
                if closest_dist==0: CNS_type = "intronic"
                elif closest_dist>0: CNS_type = "downstream"
                else: CNS_type = "upstream"
            if 'assoc_info' not in maf_entries_by_cns_id[seq_key]:
                maf_entries_by_cns_id[seq_key]['assoc_info'] = {}
            maf_entries_by_cns_id[seq_key]['assoc_info'][gen_key] = {
                "closest_info":closest_info,
                "closest_dist":closest_dist,
                "CNS_type": CNS_type
            }
        tracker.step()
    tracker.display()
    del tracker

    with open(out_folder+"results.maf","w") as out:
        tracker = pt.Progress_tracker("Exporting associations .maf",len(maf_entries_by_cns_id)).display(estimate=False,rate=1)
        for key in maf_entries_by_cns_id:
            entry = maf_entries_by_cns_id[key]
            out.write(" ".join(entry['a_line'])+"\n")
            for seq in entry['seqs']:
                istring = "i"
                if 'assoc_info' in entry and seq['s_line'][1] in entry['assoc_info']:
                    for key, value in sorted(entry['assoc_info'][seq['s_line'][1]].items()):
                        istring +=" %s='%s'" % (key,value)
                else:
                    istring += " noAssociationFound=True"
                out.write(istring+"\n")
                out.write(" ".join(seq['s_line'])+"\n")
            out.write("\n")
            tracker.step()
        tracker.display()
        del tracker

def run(argv):
    data = argv[0]
    out_folder = argv[1] if argv[1].endswith("/") else argv[1]+"/"
    main(data,out_folder)


