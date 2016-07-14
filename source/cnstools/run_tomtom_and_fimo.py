import subprocess
from _filetypes import Tomtom, Meme_v_4, Bed6
import sys
import os
from gff2_to_bed import file_run as gff2_to_bed

def main(dreme_out,motif_db,bfile,out_folder,cns_fasta,threads=1):
    if not out_folder.endswith("/"):
        out_folder+="/"
    cmd = "mkdir -p %s" % (out_folder)
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    tomtom_out = out_folder+"tom_tom_out/"
    cmd = "tomtom -oc %s -bfile %s -internal -incomplete-scores -no-ssc %s %s" % (tomtom_out,bfile,dreme_out,motif_db)
    "tomtom -oc motif_scanning/tomtom_out/ -internal -incomplete-scores -bfile fastas_for_cns/chr1.markov motif_scanning/dreme_out/dreme.txt data_files/JASPAR_CORE/jaspar.meme"
    print cmd
    # process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1)
    # while process.poll() != None:
    #     for line in iter(proc.stdout.readline, b''):
    #         sys.stdout.write(line)
    #         sys.stdout.flush()
    # process.wait()
    tom_tom_txt = tomtom_out+"tomtom.txt"
    tom_tom_data = Tomtom(file_name=tom_tom_txt)
    meme_target = Meme_v_4(file_name=motif_db)
    dreme_query = Meme_v_4(file_name=dreme_out)

    files_for_fimo_folder = out_folder+"files_for_fimo_folder/"
    try:
        os.makedirs(files_for_fimo_folder)
    except OSError:
        pass

    files_for_fimo = []
    for entry in tom_tom_data.entries:
        for match in entry.matches[:3]:
            paired_file = Meme_v_4()
            paired_file.add_entry(entry=dreme_query.lookup(match.query_id))
            paired_file.add_entry(entry=meme_target.lookup(match.target_id))
            paired_file.header = [line for line in dreme_query.header if not line.startswith('#')]
            paired_file.header.append("#tomtom match info:\n#%s\n"%match.get_line())
            f_name = files_for_fimo_folder+"%s_%s.meme.txt"%(match.query_id,match.target_id)
            paired_file.save_file(f_name)
            files_for_fimo.append(f_name)

    fimo_out = out_folder+"fimo_out/"
    try:
        os.makedirs(fimo_out)
    except OSError:
        pass

    fimo_out_file_list = []
    for filename in files_for_fimo:
        fimo_out_file = fimo_out+filename.split("/")[-1]+".gff2"
        fimo_out_file_list.append(fimo_out_file)
        cmd = "fimo --text --thresh 0.00001  --bgfile %s %s %s > %s" % (bfile,filename,cns_fasta,fimo_out_file)
        print cmd
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1)
        while process.poll() != None:
            for line in iter(proc.stdout.readline, b''):
                sys.stdout.write(line)
                sys.stdout.flush()
        process.wait()

    paired_beds = []
    for file in fimo_out_file_list:
        paired_bed = file.split(".gff")[0]+".bed"
        paired_beds.append(paired_bed)
        gff2_to_bed(file,paired_bed)
        
    temp = out_folder+"temp/"
    try:
        os.makedirs(temp)
    except OSError:
        pass

    for bed in paired_beds:
        line_dict = {}
        with open(bed) as bed_f:
            for line in bed_f.readlines():
                chro = line.split()[0]
                if chro not in line_dict:
                    line_dict[chro] = []
                line_dict[chro].append(line)
        bed_objs = [Bed6(lines=line_dict[chro]) for chro in line_dict]
        print bed_objs
        ftup = [temp+"1.bed",temp+"2.bed",temp+"1i2.bed"]
        bed_objs[0].save_file(ftup[0])
        bed_objs[0].save_file(ftup[1])
        cmd = "bedtools intersect -s -v -a %s -b %s > %s" % (ftup[0],ftup[1],ftup[2])
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1)
        while process.poll() != None:
            for line in iter(proc.stdout.readline, b''):
                sys.stdout.write(line)
                sys.stdout.flush()
        process.wait()
        ftupBedObjs = [Bed6(file_name).entries for file_name in ftup]
        counts = [len(b.entries) for b in ftupBedObjs]
        print zip(ftup,counts)







def file_run(dreme_out,motif_db,bfile,out_folder,cns_fasta,threads=1):
    main(dreme_out,motif_db,bfile,out_folder,cns_fasta,threads)