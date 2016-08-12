import matplotlib as mpl
mpl.use('SVG')
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
import numpy as np
from filetype_classes import Cns
import sys
import argparse

def parser(parser_add_func,name):
    p = parser_add_func(name,description="")
    p.add_argument("cns_file", help="Path to .cns file.")
    p.add_argument("genome", help="Reference genome.")
    p.add_argument("fimo_file", help="Path to fimo out file.")
    p.add_argument('-o','--outfile',default="fimo_compare_out.svg")
    return p

def run(cns_file,genome,fimo_file,outfile="fimo_compare_out.svg"):
    with open(cns_file) as line_file:
        cns = Cns(lines=line_file.readlines())
    fimo_data_dict = {}
    motif_best_sequece = {}
    with open(fimo_file) as fimo:
        for line in fimo.readlines():
            if not line.startswith("#"):
                full = line.strip().split('\t')
                motif_id,cns_id,motif_seq = (full[i] for i in [0,1,8])
                if not cns_id in fimo_data_dict:
                    fimo_data_dict[cns_id] = []
                if not motif_id in motif_best_sequece:
                    motif_best_sequece[motif_id] = motif_seq
                fimo_data_dict[cns_id].append(motif_id)
    type_counts = {}
    motif_type_dist = {}
    for entry in cns.entries:
        for seq in entry.get_seqs(genome):
            if seq.cns_ID in fimo_data_dict:
                for motif in fimo_data_dict[seq.cns_ID]:
                    if motif not in motif_type_dist:
                        motif_type_dist[motif] = {}
                    if seq.type not in motif_type_dist[motif]:
                        motif_type_dist[motif][seq.type] = 0
                    if seq.type not in type_counts:
                        type_counts[seq.type] = 0
                    motif_type_dist[motif][seq.type] += 1
                    type_counts[seq.type] += 1
    motif_names = motif_type_dist.keys()
    type_names = list(set([name for key in motif_type_dist for name in motif_type_dist[key]]))
    for dict in motif_type_dist:
        motif_type_dist[dict].update({t:0 for t in type_names if t not in motif_type_dist[dict]})
    type_names.sort(key=lambda x:-sum(motif_type_dist[g][x] for g in motif_names))
    data_array = [[motif_type_dist[x][y] for y in type_names] for x in motif_names]
    for col in range(len(data_array[0])):
        for line in data_array:
            line[col]/=float(type_counts[type_names[col]])
    total_occurences = sum(type_counts[t] for t in type_counts)
    for row in data_array:
        s = sum([val*type_counts[type_names[i]]/total_occurences for i,val in enumerate(row)]) 
        row[:] = [n/float(s) for n in row]
    motif_names = sorted(motif_names,key=lambda m:max((data_array[motif_names.index(m)][t],-t) for t in range(len(type_names)))[::-1])[::-1]
    data_array.sort(key = lambda l: max((l[t],-t) for t in range(len(type_names)))[::-1])
    data_array.reverse()
    for line in data_array: print line
    _heatmap(data_array,"fimo_comap.svg","type","motif",type_names,[motif_best_sequece[motif_id] for motif_id in motif_names])

def _heatmap(data,out_file_path,xName="",yName="",xlabels=None,ylabels=None,y2labels=None,plt_title=""):
    plt.clf()
    plt.figure(figsize=(len(data[0]),len(data)/6.0))
    highest = max(max(row) for row in data)
    ax = plt.gca()
    im = ax.imshow(data, extent=[0,len(data[0]),0,len(data)], cmap='cool',aspect='auto', interpolation='nearest',vmin=0)
    ax.set_title(plt_title,fontsize=14, y=1.2)
    ax.set_xlabel(xName)
    ax.set_ylabel(yName)
    plt.colorbar(im,ticks=[0, highest/2.0,highest])
    if ylabels:
        ax.set_yticks(np.arange(len(ylabels)) + 0.5)
        ax.set_yticklabels(ylabels[::-1], fontsize=10)
    if y2labels:
        ax2 = ax.twinx()
        ax2.set_yticks(np.arange(len(y2labels)) + 0.5)
        ax2.set_yticklabels(y2labels[::-1], fontsize=10)
    if xlabels:
        ax.set_xticks(np.arange(len(xlabels)) + 0.5)
        ax.xaxis.tick_top()
        ax.set_xticklabels(xlabels,rotation=90)
        ax.xaxis.set_label_position('top') 
    plt.savefig(out_file_path, bbox_inches='tight')
