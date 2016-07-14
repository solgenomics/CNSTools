from filetypes import Cns
import numpy as np
import scipy as sp
import matplotlib as mpl
mpl.use('SVG')
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors

def main(do_list,cns):

    if "count" in do_list: print len("%s CNSs" % cns.entries)

    if "types" in do_list:
        list_of_cns_types = ["intergenic","intronic","downstream","upstream"]
        numGenoms = 2

        type_list=list_of_cns_types+["unknown","DNE"]
        # list_of_cns_types = type_list
        type_count_dict = {key:0 for key in type_list}
        type_pair_count_dict = {(key1,key2):0 for key1 in type_list for key2 in type_list}

        for entry in cns.entries:
            entry_count_dict = {}
            entry_pair_counts = {}
            for seq in entry.get_seqs():
                if seq.genome=="chr1" or seq.genome=="Gmax":
                    t = seq.type if seq.type!=None else "unknown"
                    if t not in entry_count_dict: 
                        entry_count_dict[t] = 0
                    entry_count_dict[t]+=1
            entry_count_dict["DNE"] = numGenoms - sum((entry_count_dict[key] for key in entry_count_dict))
            entry_types = entry_count_dict.keys()
            for i in range(len(entry_types)):
                for j in range(len(entry_types)):
                    key = (entry_types[i],entry_types[j])
                    if key not in entry_pair_counts: 
                        entry_pair_counts[key] = entry_count_dict[key[0]] * (entry_count_dict[key[1]] - (1 if key[0]==key[1] else 0))
            for key in entry_count_dict:
                type_count_dict[key]+=entry_count_dict[key]
            for key in entry_pair_counts:
                type_pair_count_dict[key]+=entry_pair_counts[key]


        type_prob_dict = {key:type_count_dict[key]/float(numGenoms*len(cns.entries)) for key in type_count_dict}

        type_pair_prob_dict = {key:type_pair_count_dict[key]/float(numGenoms*(numGenoms-1)*len(cns.entries)) for key in type_pair_count_dict}

        type_conditional_prob_dict = {key:type_pair_prob_dict[key]/type_prob_dict[key[0]] for key in type_pair_prob_dict}

        conditional_prob_array = [[0]*len(list_of_cns_types) for i in range(len(list_of_cns_types))]
        for i in range(len(list_of_cns_types)):
                for j in range(len(list_of_cns_types)):
                    key = (list_of_cns_types[i],list_of_cns_types[j])
                    probability = type_conditional_prob_dict[key]
                    conditional_prob_array[i][j] = probability

        print
        print ("px given y" + " ".join(['%14s'%(int(val)) for val in range(len(conditional_prob_array))]))
        for i in range(len(conditional_prob_array)):
            row = conditional_prob_array[i]
            print ("%s %-10s " % (i,list_of_cns_types[i])) + " ".join(['%14s'%(val) for val in row])


        heatmap(conditional_prob_array,"seq_conditional_heatmap.svg","type X","type Y",list_of_cns_types,list_of_cns_types,"Probability that a CNS is known to be of type X on any genome, given that it is known to be type Y on one.")

    if "len_dist" in do_list:
        dists = []
        lens = []
        typeLens = {}
        for entry in cns.entries:
            for seq in entry.get_seqs():
                if seq.dist!=None:
                    if seq.dist!=0:
                        lens.append(seq.stop-seq.start)
                        dists.append(seq.dist)
                    if seq.type not in typeLens: typeLens[seq.type] = []
                    typeLens[seq.type].append(seq.stop-seq.start)

        heatmap_histogram((dists,lens),"dist_len_heatmap.svg")
        histogram(typeLens,"type_len_histogram.svg")

    if "nucleotide_freq" in do_list:
        nucs = ('A','T','C','G')
        typefreqs = {}
        for entry in cns.entries:
            for seq in entry.get_seqs():
                if seq.type!=None:
                    if seq.type not in typefreqs: 
                        typefreqs[seq.type] = {'A':0,'T':0,'C':0,'G':0,'__denom':0}
                    typefreqs[seq.type]['__denom']+= 1
                    for nuc in nucs:
                        typefreqs[seq.type][nuc]+= seq.sequence.count(nuc)/float(len(seq.sequence))
        for t in typefreqs:
            for nuc in nucs:
                typefreqs[t][nuc] = typefreqs[t][nuc]/float(typefreqs[t]['__denom'])
            del typefreqs[t]['__denom']
        multi_bar_graph(typefreqs,"nucleotide_freq_graph.svg")

def scatter(data,out_file_path):
    plt.clf()
    twoDArr = data
    unzipped = zip(*twoDArr)
    plt.scatter(unzipped[1],unzipped[0])
    plt.xlabel('')
    plt.ylabel('')
    plt.grid(True)
    plt.savefig(out_file_path, bbox_inches='tight')

def heatmap_histogram(data,out_file_path):
    plt.clf()
    plt.xlabel('')
    plt.ylabel('')
    plt.hexbin(data[0], data[1], cmap=plt.cm.jet, norm=mpl.colors.LogNorm())
    plt.savefig(out_file_path, bbox_inches='tight')

def heatmap(data,out_file_path,xName="",yName="",xlabels=None,ylabels=None,plt_title=""):
    plt.clf()
    biggest = max((val for row in data for val in row))
    ax = plt.gca()
    im = ax.imshow(data, extent=[0,len(data),0,len(data[0])], aspect=1.0, interpolation='nearest', vmin=0,vmax=biggest)
    ax.set_title(plt_title,fontsize=14, y=1.1)
    ax.set_xlabel(xName)
    ax.set_ylabel(yName)
    plt.colorbar(im)
    if xlabels:
        ax.set_xticks(np.arange(len(xlabels)) + 0.5)
        ax.xaxis.tick_top()
        ax.set_xticklabels(xlabels)
        ax.xaxis.set_label_position('top') 
    if ylabels:
        ax.set_yticks(np.arange(len(ylabels)) + 0.5)
        ax.set_yticklabels(ylabels[::-1])
    for i in range(len(data)):
        for j in range(len(data[i])):
            ax.annotate(('%s'%(data[i][j]))[:5], xy = (j+0.5, (len(data[i])-i)-0.5),
                ha = 'center', va = 'center')
                #bbox = dict(boxstyle = 'square,pad=0.5', fc = 'white', alpha = 0.3))
    plt.savefig(out_file_path, bbox_inches='tight')

def histogram(data_dict,out_file_path):
    plt.clf()
    plt.hold(True)
    all_min = min((min(data_dict[key]) for key in data_dict))
    all_max = max((max(data_dict[key]) for key in data_dict))
    bins = np.linspace(all_min, all_max, 50)
    for key in data_dict:
        plt.hist(data_dict[key], bins, alpha=0.3, label=key,log=True)
    plt.legend(loc='upper right')
    plt.savefig(out_file_path, bbox_inches='tight')

def bar_graph(dict,out_file_path):
    plt.clf()
    plt.bar(range(len(dict)), dict.values(), align='center')
    plt.xticks(range(len(dict)), dict.keys())
    plt.savefig(out_file_path, bbox_inches='tight')

def multi_bar_graph(dict_dict,out_file_path):
    plt.clf()
    plt.hold(True)
    keys = dict_dict.keys()
    group_pos = np.arange(len(dict_dict[keys[0]]))
    bar_width = 0.20
    for i in range(len(keys)):
        dict = dict_dict[keys[i]]
        vals = np.array(dict.values())
        stdev = sp.stats.sem(vals)
        plt.bar(group_pos+(i*bar_width), vals, bar_width, label=keys[i])
        plt.errorbar(group_pos+((i+0.5)*bar_width), vals, yerr=stdev, fmt='o',capsize=5,capthick=1,elinewidth=1,markersize=0,color='black')
    plt.xticks(group_pos+(len(keys)/2)*bar_width, dict_dict[keys[0]].keys())
    plt.legend(loc='upper right')
    plt.savefig(out_file_path, bbox_inches='tight')



def file_run(cns_file,*args): 

    do_list = list(args)
    print do_list
    cns = Cns()
    with open(cns_file) as file:
        cns.add_lines(file.readlines())

    main(do_list,cns)

