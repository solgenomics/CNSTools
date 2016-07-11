from filetypes import Cns
import numpy as np
import scipy as sp
import matplotlib as mpl
mpl.use('SVG')
import matplotlib.pyplot as plt

def main(do_list,cns):

    if "count" in do_list: print len("%s CNSs" % cns.entries)

    if "types" in do_list:
        seq_type_count = {}
        paired_seq_type = {}
        entry_counted = 0
        for entry in cns.entries:
            type_list = [seq.type for seq in entry.get_seqs() if seq.type!=None]
            types = tuple(set(type_list))
            if len(types)>0:
                entry_counted+=1
                for i in range(len(types)):
                    if not types[i] in seq_type_count: 
                        seq_type_count[types[i]] = 1
                    else: 
                        seq_type_count[types[i]]+= 1
                    exclusion_list = type_list[:]
                    exclusion_list.remove(types[i])
                    exclusion_set = tuple(set(exclusion_list))
                    for j in range(len(exclusion_set)):
                        key = tuple(sorted([types[i],exclusion_set[j]]))
                        if not key in paired_seq_type: 
                            paired_seq_type[key] = 1
                        else: 
                            paired_seq_type[key]+= 1
                    type_list[:] = [t for t in type_list if t!=types[i]]

        seq_type_ratio = {key:seq_type_count[key]/float(entry_counted) for key in seq_type_count}

        paired_seq_type_ratio = {key:paired_seq_type[key]/float(entry_counted) for key in paired_seq_type}

        type_to_index = {}
        index_to_type = {}
        i = 0
        for key in seq_type_ratio:
            type_to_index[key] = i
            index_to_type[i] = key
            i+=1

        conditional_prob_array = [[0]*len(type_to_index) for i in range(len(type_to_index))]
        for two in paired_seq_type_ratio:
            for pair in (two,two[::-1]):
                probability = paired_seq_type_ratio[two]/seq_type_ratio[pair[1]]
                p0_index = type_to_index[pair[0]]
                p1_index = type_to_index[pair[1]]
                conditional_prob_array[p0_index][p1_index] = probability

        print "probability that a CNS is of type on at least one seq", seq_type_ratio
        print "probability of type given type", paired_seq_type_ratio

        print "conditional probability of type given type, normalized by type freq"
        for i in range(len(conditional_prob_array)):
            row = conditional_prob_array[i]
            print ("%s %-10s " % (i,index_to_type[i])) + " ".join(['%3s'%(int(val*100)) for val in row])

        labels = [index_to_type[i] for i in range(len(index_to_type))]
        heatmap(conditional_prob_array,"seq_conditional_heatmap.svg",labels,labels)

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
    plt.imshow(data, extent=[0,len(data),0,len(data[0])], aspect=1.0, cmap=plt.cm.inferno, interpolation='nearest')
    plt.xlabel('Given:')
    plt.ylabel('Probability of:')
    plt.colorbar()
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

    do_list = list(*args)

    cns = Cns()
    with open(cns_file) as file:
        cns.add_lines(file.readlines())

    main(do_list,cns)

