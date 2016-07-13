from filetypes import Cns
import numpy as np
import scipy as sp
import matplotlib as mpl
mpl.use('SVG')
import matplotlib.pyplot as plt
from matplotlib import cm

def main(do_list,cns):

    if "count" in do_list: print len("%s CNSs" % cns.entries)

    if "types" in do_list:
        type_list = ["intergenic","intronic","downstream","upstream"]
        type_count_dict = {key:0 for key in type_list}
        type_exist_dict = {key:0 for key in type_list}
        type_prob_dict = {key:0 for key in type_list}
        type_conditional_count_dict = {(key1,key2):0 for key1 in type_list for key2 in type_list}
        type_conditional_prob_dict = {(key1,key2):0 for key1 in type_list for key2 in type_list}
        numGenoms = 10
        for entry in cns.entries:
            entry_count_dict = {}
            entry_prob_dict = {}
            entry_conditional_counts = {}
            for seq in entry.get_seqs():
                if seq.type!=None:
                    if seq.type not in entry_count_dict: 
                        entry_count_dict[seq.type] = 0
                    entry_count_dict[seq.type]+=1
            entry_types = entry_count_dict.keys()
            for t in set(entry_types):
                type_exist_dict[t]+=1
            for i in range(len(entry_types)):
                for j in range(len(entry_types)):
                    key = (entry_types[i],entry_types[j])
                    if key not in entry_conditional_counts: 
                        entry_conditional_counts[key] = 0
                    entry_conditional_counts[key]+= entry_count_dict[key[1]] - (1 if key[0]==key[1] else 0)
            for key in entry_count_dict:
                type_count_dict[key]+=entry_count_dict[key]
                entry_prob_dict[key] = entry_count_dict[key]/float(numGenoms)
                type_prob_dict[key]+=entry_prob_dict[key]
            for key in entry_conditional_counts:
                type_conditional_count_dict[key]+=entry_conditional_counts[key]
                type_conditional_prob_dict[key]+=entry_conditional_counts[key]/float(numGenoms-1)

        type_prob_dict = {key:type_prob_dict[key]/float(len(cns.entries)) for key in type_prob_dict}
        type_conditional_prob_dict = {key:type_conditional_prob_dict[key]/type_exist_dict[key[0]] for key in type_conditional_prob_dict}

        conditional_prob_array = [[0]*len(type_list) for i in range(len(type_list))]
        for i in range(len(type_list)):
                for j in range(len(type_list)):
                    key = (type_list[i],type_list[j])
                    probability = type_conditional_prob_dict[key]
                    conditional_prob_array[i][j] = probability

        print "Conditional probs"
        print ("%s %-10s " % (i,"")) + " ".join(['%s'%(int(val)) for val in range(len(conditional_prob_array))])
        for i in range(len(conditional_prob_array)):
            row = conditional_prob_array[i]
            print ("%s %-10s " % (i,type_list[i])) + " ".join(['%s'%(val) for val in row])

        heatmap(conditional_prob_array,"seq_conditional_heatmap.svg",type_list,type_list)

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

def heatmap(data,out_file_path,xlabels=None,ylabels=None):
    plt.clf()
    im = plt.imshow(data, extent=[0,len(data),0,len(data[0])], aspect=1.0, cmap=cm.hot, interpolation='nearest', vmin=0, vmax=0.3)
    plt.xlabel('Given:')
    plt.ylabel('Probability of:')
    plt.colorbar(im)
    ax = plt.gca()
    if xlabels:
        plt.xticks([i+0.5 for i in range(len(xlabels))], xlabels)
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top') 
    if ylabels:
        plt.yticks([i+0.5 for i in range(len(ylabels))], ylabels[::-1])
    for i in range(len(data)):
        for j in range(len(data[i])):
            plt.annotate(('%s'%(data[i][j]))[:4], xy = (j+0.5, (len(data[i])-i)-0.5),
                ha = 'center', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.5', fc = 'white', alpha = 0.3))
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
        cns.add_lines(file.readlines()[:])

    main(do_list,cns)

