"""A script for analyzing the contents of a .cns file or :class:`.Cns` object.
"""

from _filetypes import Cns
import _utils
import numpy as np
import scipy as sp
import matplotlib as mpl
mpl.use('SVG')
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors


def _main(cns,task_list):
    """This function runs the main workflow for the script.

    :param cns: The :class:`.Cns` object to be used in the program.
    :param  list[str] task_list: A list of tasks to perform, options are ``["count","types","len_dist","nucleotide_freq"]``
    :type cns: :class:`.Cns`
    :returns:  `None`
    """
    if "count" in task_list: 
        print len("%s CNSs" % cns.entries)

    if "type_conditional_probs" in task_list:
        list_of_cns_types = ["intergenic","intronic","downstream","upstream"]
        numGenoms = 7

        type_list=list_of_cns_types+["unknown","DNE"]
        # list_of_cns_types = type_list
        type_count_dict = {key:0 for key in type_list}
        type_pair_count_dict = {(key1,key2):0 for key1 in type_list for key2 in type_list}

        for entry in cns.entries:
            entry_count_dict = {}
            entry_pair_counts = {}
            for seq in entry.get_seqs():
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

        print type_count_dict
        type_prob_dict = {key:type_count_dict[key]/float(numGenoms*len(cns.entries)) for key in type_count_dict}

        type_pair_prob_dict = {key:type_pair_count_dict[key]/float(numGenoms*(numGenoms-1)*len(cns.entries)) for key in type_pair_count_dict}

        type_conditional_prob_dict = {key:type_pair_prob_dict[key]/(type_prob_dict[key[0]] if type_prob_dict[key[0]]!=0 else float('inf')) for key in type_pair_prob_dict}

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

        _heatmap(conditional_prob_array,"seq_conditional_heatmap.svg","type X","type Y",list_of_cns_types,list_of_cns_types,"Probability that a CNS is known to be of type X on any genome, given that it is known to be type Y on one.")

    if "len_dist" in task_list:
        dists = []
        lens = []
        typeLens = {}
        for entry in cns.entries:
            for seq in entry.get_seqs():
                if seq.stop and seq.start:
                    if seq.dist and abs(seq.dist)>0:
                        lens.append(seq.stop-seq.start)
                        dists.append(seq.dist)
                    if seq.type not in typeLens: typeLens[seq.type] = []
                    typeLens[seq.type].append(seq.stop-seq.start)

        _heatmap_histogram((dists,lens),"dist_len_heatmap.svg",xName="Downstream distance of CNS from nearest gene (bp)",yName="Length of CNS (bp)",plt_title="Distribution of CNS by Length and Distance")
        _histogram(typeLens,"type_len_histogram.svg",title="CNS Length Histogram",xName="CNS length (bp)")

    if "nucleotide_freq" in task_list:
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
        _multi_bar_graph(typefreqs,"nucleotide_freq_graph.svg")
    if "types" in task_list:
        type_count_dict = {}
        for entry in cns.entries:
            for seq in entry.get_seqs():
                if seq.type not in type_count_dict:
                    type_count_dict[seq.type] = 0
                type_count_dict[seq.type]+=1
        _bar_graph(type_count_dict,"type_bar_graph.svg",title="CNS Type Count",xName="CNS type")


def _heatmap_histogram(data,out_file_path,xName="",yName="",plt_title=""):
    plt.clf()
    ax = plt.gca()
    ax.set_title(plt_title,fontsize=14, y=1.1)
    ax.set_xlabel(xName)
    ax.set_ylabel(yName)
    im = ax.hexbin(data[0], data[1], cmap=plt.cm.jet, norm=mpl.colors.LogNorm())
    plt.colorbar(im)
    plt.savefig(out_file_path, bbox_inches='tight')

def _heatmap(data,out_file_path,xName="",yName="",xlabels=None,ylabels=None,plt_title=""):
    plt.clf()
    biggest = max((val for row in data for val in row))
    ax = plt.gca()
    im = ax.imshow(data, extent=[0,len(data),0,len(data[0])], cmap=plt.cm.jet,aspect=1.0, interpolation='nearest', vmin=0,vmax=biggest)
    ax.set_title(plt_title,fontsize=14, y=1.2)
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
    plt.savefig(out_file_path, bbox_inches='tight')

def _histogram(data_dict,out_file_path,title="",xName="",yName="Count"):
    plt.clf()
    plt.hold(True)
    n_bins = 20 
    ax = plt.gca()
    data = zip(*[[key,data_dict[key]] for key in data_dict])
    labels = data[0]
    list_list = data[1]
    ax.hist(list_list, n_bins, histtype='bar' ,label=labels,rwidth=0.9,log=True)
    ax.legend(prop={'size': 10})
    ax.set_title(title)
    ax.set_xlabel(xName)
    ax.set_ylabel(yName)
    plt.savefig(out_file_path, bbox_inches='tight')

def _bar_graph(dict,out_file_path,title="",xName="",yName="Count"):
    plt.clf()
    ax = plt.gca()
    ax.bar(range(len(dict)), dict.values(), align='center')
    ax.set_xticks(range(len(dict)))
    ax.set_xticklabels(dict.keys())
    ax.set_title(title)
    ax.set_xlabel(xName)
    ax.set_ylabel(yName)
    plt.savefig(out_file_path, bbox_inches='tight')

def _multi_bar_graph(dict_dict,out_file_path):
    plt.clf()
    plt.hold(True)
    keys = dict_dict.keys()
    group_pos = np.arange(len(dict_dict[keys[0]]))
    bar_width = 0.20
    for i in range(len(keys)):
        dict = dict_dict[keys[i]]
        vals = np.array(dict.values())
        stdev = np.std(vals)
        plt.bar(group_pos+(i*bar_width), vals, bar_width, label=keys[i])
        plt.errorbar(group_pos+((i+0.5)*bar_width), vals, yerr=stdev, fmt='o',capsize=5,capthick=1,elinewidth=1,markersize=0,color='black')
    plt.xticks(group_pos+(len(keys)/2)*bar_width, dict_dict[keys[0]].keys())
    plt.legend(loc='upper right')
    plt.savefig(out_file_path, bbox_inches='tight')

def run(cns_file,*tasks):
    """This function is for running the workflow on a file and not a :class:`.Cns` object. 
    It creates a :class:`.Cns` object from the file `cns_file` and compiles the `tasks` arguements into a list of strings then calls :func:`.main` using them.
    
    :param str cns_file: A path to the .cns file to analyze.
    :param *args tasks: Any number of arguements, each of which is the name of a task to perform on `cns_file`.
    :returns:  `None`
    """
    task_list = list(tasks)
    _utils.safe_print(task_list)
    cns = Cns()
    with open(cns_file) as file:
        cns.add_lines(file.readlines())

    _main(cns,task_list)

