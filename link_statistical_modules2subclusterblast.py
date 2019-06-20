#!/usr/bin/env python3
"""
Author: Joris Louwen
Script to find modules with LDA algorithm.
"""

import argparse
from collections import Counter, defaultdict
from functools import partial
import matplotlib
matplotlib.use('Agg') #to not rely on X-forwarding (not available in screen)
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
from numpy import sqrt
import numpy as np
from operator import itemgetter
import os
import pandas as pd
# import scipy.cluster.hierarchy as sch
# import seaborn as sns
from statistics import mean,median
import subprocess
from sys import argv
import time


def get_commands():
    parser = argparse.ArgumentParser(description="A script to link modules\
        from the statistical method with known subclusters.")
    parser.add_argument("-i", "--bgcfile", dest="bgcfile", help="Input \
        csv file of BGCs with genes as domain combinations", required=True)
    parser.add_argument("-m", "--modfile", dest="modfile", help="Input \
        txt file of putative modules. Last column should contain\
        modules", default=False)
    parser.add_argument("-o", "--out_folder", dest="out_folder", help="Output\
        folder", required=True)
    parser.add_argument("-c", "--cores", dest="cores", help="Amount \
        of cores to use, default = all available cores",\
        default=cpu_count(), type=int)
    parser.add_argument("--classes", help="A file containing classes of the \
        BGCs used in the analysis. First column should contain matching BGC\
        names. Consecutive columns should contain classes.", default=False)
    parser.add_argument("--plot", help="If provided: make plots about \
        several aspects of the output. Default is off.", default=False, \
        action="store_true")
    parser.add_argument("--known_subclusters", help="A tab delimited file \
        with known subclusters. Should contain subclusters in the last column\
        and BGC identifiers in the first column. Subclusters are comma \
        separated genes represented as domains. Multiple domains in a gene \
        are separated by semi-colon.")
    parser.add_argument("--bgc_with_mods",help="A tab delimited file linking\
        bgcs to modules. First column should be the name of a bgc and the \
        last column should be the module numbers separated by '; '. These \
        numbers should be the first column in the modfile")
    return parser.parse_args()

def line_plot_known_matches(known_subcl_matches, outname, cutoff,steps=0.1):
    '''Plot a line of the amount of known_subcl matches with different cutoffs


    Matches are only reported if at least two genes match, these can be two
    of the same genes if the prob is 1.5 or higher (close enough to two)
    '''
    ys=[cutoff+i*steps for i in range(round((1.0-cutoff)/steps)+1)]
    xs=[0]*len(ys)
    for info in known_subcl_matches.values():
        if len(info) > 0:
            for i,thresh in enumerate(ys):
                for overlap in info:
                    if overlap[0] >= thresh and overlap[1] > 1:
                        xs[i]+=1
                        break
    print('There are {} known sub_clusters with at least one'.format(xs[0])+
        ' overlap of two genes with a minimum overlap % of {}'.format(ys[0]))
    plt.plot(ys, xs)
    plt.xlabel('Overlap threshold')
    plt.ylabel('Characterised subclusters with an overlap')
    plt.title(\
    'Number of characterised subclusters overlapping with the statistical \n\
        method according to different overlap thresholds')
    plt.savefig(outname)
    plt.close()

def read2dict(filepath, sep=',',header=False):
    '''Read file into a dict {first_column:[other_columns]}

    filepath: str
    sep: str, delimiter in the file
    header: bool, ignore first line
    '''
    output = {}
    with open(filepath,'r') as inf:
        if header:
            inf.readline()
        for line in inf:
            line = line.strip().split(sep)
            output[line[0]] = line[1:]
    return output

def find_stat_method_overlap(mods, known_subcl, bgc_classes, cutoff, \
    mods2bgc, outfile):
    '''

    mods: dict of {module_tuple:[number,other_info]}
    known_subcl: {bgc: [[info,domains]]}
    bgc_classes: {bgc: [class1,class2]}
    cutoff: float
    mods2bgc: {mod_number:[bgcs]}
    outfolder: str,filepath
    match_dict: {known_subcl_info:[%overlap,len_overlap, mod_number,\
        overlapping_genes,non_overlapping_genes,[bgc_and_class_str]]} 
        known_subcl_info in this case is the first info element from info in
        known_subcl
    '''
    print('\nLinking modules to known subclusters')
    match_dict = defaultdict(list)
    #loop known_subcl and then each mod
    for infos in known_subcl.values():
        for info in infos:
            name_k = info[0]
            doms_k = info[-1].split(',')
            doms_k_set = set(doms_k)
            if '-' in doms_k_set:
                doms_k_set.remove('-')
                doms_k = [k for k in doms_k if not k =='-']
            if doms_k:
                for mod in mods:
                    mod_set = set(mod)
                    dom_overlap_set = doms_k_set&mod_set
                    l_overlap = len(dom_overlap_set)
                    overlap = round(l_overlap / len(doms_k_set),3)
                    if overlap >= cutoff:
                        mod_num = mods[mod][0]
                        bgcs = sorted(mods2bgc[mod_num])
                        overl_genes = ','.join(sorted(dom_overlap_set))
                        non_overl_genes = ','.join(sorted(\
                            mod_set.difference(doms_k_set)))
                        bgc_list=[]
                        for bgc in bgcs:
                            bgc_class = bgc_classes.get(bgc,['None'])[0]
                            bgc_list.append('{}\t{}'.format(bgc,bgc_class))
                        match = [overlap,l_overlap,mod_num, overl_genes,\
                            non_overl_genes,bgc_list]
                        match_dict[name_k].append(match)
    write_stat_method_overlap(match_dict,known_subcl,cutoff,outfile)
    plotname = outfile.split('.txt')[0] + '_vs_cutoff.pdf'
    line_plot_known_matches(match_dict,plotname,cutoff)

def write_stat_method_overlap(match_dict, known_subcl, cutoff, outfile):
    '''Write overlap to file
    '''
    print('  writing overlap to {}'.format(outfile))
    with open(outfile,'w') as outf:
        outf.write('##Values below each subcluster: %overlap len_overlap'+
            ' mod_num overlap_genes non_overlap_genes\n'+
            '##Below that are the BGCs that have the module with class\n')
        for bgc,k_subcl_info in sorted(known_subcl.items(),\
            key=lambda x: x[1][0][0]):
            for k_subcl in k_subcl_info:
                outf.write('#{}\t{}\n'.format(bgc,'\t'.join(map(str,\
                    k_subcl))))
                overlap_list = match_dict[k_subcl[0]]
                for m_overlap in sorted(overlap_list, key=lambda x:\
                    (-x[0],len(x[4]),int(x[2]))):
                    bgc_list = m_overlap.pop(-1)
                    outf.write('{}\n'.format('\t'.join(map(str,m_overlap))))
                    for bgc_class in bgc_list:
                        outf.write('\t{}\n'.format(bgc_class))


# with open(subcl_out,'w') as outf:
            # #sort the subclusters alphabetically on first info element
            # outf.write('##Values below each subcluster: %overlap len_overlap'+
                # ' bgc class topic topic_probability overlap_genes'+
                # ' non_overlap_genes\n')
            # for bgc, info in sorted(known_subcl.items(),\
                # key=lambda x: x[1][0][0]):
                # for k_subclust in info:
                    # outf.write('#{}\t{}\n'.format(bgc,'\t'.join(map(str,\
                        # k_subclust))))
                    # overlap_list = known_subcl_matches[k_subclust[0]]
                    # #give summary per topic?
                    # #e.g. #topic x: 12 avg_overlap: 0.403 
                    # #sort from high to low overlap,topic,bgc
                    # for m_overlap in sorted(overlap_list, key=lambda x: \
                        # (-x[0],x[4],x[2])):
                        # #overlap bgc class topic prob genes:prob
                        # outf.write('{}\n'.format('\t'.join(\
                            # map(str,m_overlap))))

if __name__ == '__main__':
    start = time.time()
    #files provided should be filtered bgc csv file and filtered module file

    print('\nStart')
    cmd = get_commands()
    bgcs = read2dict(cmd.bgcfile)
    with open(cmd.modfile, 'r') as inf:
        modules = {}
        #{modules:[info]}
        for line in inf:
            line = line.strip().split('\t')
            mod = tuple(line[-1].split(',')) #now a tuple of str
            modules[mod] = line[:-1]
    if cmd.classes:
        bgc_classes_dict = read2dict(cmd.classes, sep='\t',header=True)
    else:
        bgc_classes_dict = {bgc:'None' for bgc in bgcs}
    if not os.path.isdir(cmd.out_folder):
        subprocess.check_call('mkdir {}'.format(cmd.out_folder), shell=True)

    if cmd.known_subclusters:
        known_subclusters = defaultdict(list)
        with open(cmd.known_subclusters,'r') as inf:
            for line in inf:
                line = line.strip().split('\t')
                known_subclusters[line[0]].append(line[1:])
    else:
        known_subclusters = False

    mod_nums2bgc=defaultdict(list)
    with open(cmd.bgc_with_mods,'r') as inf:
        for line in inf:
            line = line.strip('\n').split('\t')
            mod_nums = line[-1].split('; ')
            if mod_nums:
                for m_num in mod_nums:
                    mod_nums2bgc[m_num].append(line[0])

    prefix = cmd.bgcfile.split('.csv')[0].split('filtered_clusterfile')[0]
    out_file = os.path.join(cmd.out_folder,prefix+\
        'known_subcluster_overlap_statistical_method.txt')
    find_stat_method_overlap(modules, known_subclusters, bgc_classes_dict,\
        0.4, mod_nums2bgc, out_file)

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
