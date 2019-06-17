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
import scipy.cluster.hierarchy as sch
import seaborn as sns
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
        ' match of two genes with a minimum overlap of {}'.format(ys[0]))
    plt.plot(ys, xs)
    plt.xlabel('Overlap threshold')
    plt.ylabel('Characterised subclusters with a match')
    plt.title(\
    'Number of characterised subclusters with a match according\n\
        to different overlap thresholds')
    plt.savefig(outname)
    plt.close()

def compare_known_subclusters(known_subcl, bgc, bgc_class, matches,cutoff):
    '''Find % overlap with known subclusters and returns it as a list

    known_subcl: {bgc: [[info,domains]]}
    bgc: str, bgcname
    bgc_class: str, class of bgc
    matches: [[domain_combinations]]
    cutoff: float, overlap cutoff used for reporting
    matches_overlap: [[first_info_element,%overlap,len_overlap,bgc,bgc_class,
        topic_num,prob,overlapping_genes,non_overlapping_genes]]
    '''
    matches_overlap = []
    for match in matches:
        g_list = match[2]
        doms = set(list(zip(*g_list))[0])
        for k_subs in known_subcl.values():
            for k_sub in k_subs:
                k_list = k_sub[-1].split(',')
                k_sub_doms = set(k_sub[-1].split(','))
                if '-' in k_sub_doms:
                    k_sub_doms.remove('-')
                    k_list = [k for k in k_list if not k =='-']
                overl_d_set = doms&k_sub_doms
                l_overlap = len(overl_d_set)
                if not len(k_sub_doms) - len(k_list) == 0:
                    #there are doms in the k-subcl that are duplicated
                    dupls = [kc for kc in Counter(k_list).items() if kc[1]>1]
                    add_overl = 0
                    for dom,count in dupls:
                        if dom in doms:
                            overl_domtups = [domt for domt in g_list \
                                if domt[0]==dom]
                            for overl_domtup in overl_domtups:
                                if round(overl_domtup[1]) >= count:
                                    l_overlap += count-1
                overlap = l_overlap / len(k_list)
                if overlap > cutoff and len(k_list) > 1:
                    match_overl_genes = [(g,p,) for g,p in\
                        g_list if g in overl_d_set]
                    overl_d = ','.join(sorted([g+':'+str(p) for g,p in\
                        match_overl_genes]))
                    non_overl_d = ','.join(sorted([g+':'+str(p) for g,p in\
                        g_list if not g in overl_d_set]))

                    matches_overlap.append([k_sub[0],round(overlap,3),\
                        l_overlap,bgc,bgc_class,match[0],match[1],overl_d,\
                        non_overl_d])
    return matches_overlap

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
    mods2bgc, outfolder):
    '''

    mods: dict of {module_tuple:[number,other_info]}
    known_subcl: {bgc: [[info,domains]]}
    bgc_classes: {bgc: [class1,class2]}
    cutoff: float
    mods2bgc: {mod_number:[bgcs]}
    outfolder: str,filepath
    match_dict: {known_subcl:[%overlap,len_overlap,bgc,bgc_class,
        topic_num,prob,overlapping_genes,non_overlapping_genes]}
    '''
    match_dict = defaultdict(list)
    #loop known_subcl and then each mod
    for info in known_subcl.values():
        matches_k = []
        name_k = info[0]
        doms_k = info[-1].split(',')
        doms_k_set = set(k_sub[-1].split(','))
        if '-' in doms_k_set:
            doms_k_set.remove('-')
            doms_k = [k for k in doms_k if not k =='-']
        if doms_k:
            for mod in mods:
                mod_set = set(mod)
                dom_overlap_set = doms_k_set&mod_set
                l_overlap = len(dom_overlap_set)
                overlap = l_overlap / len(doms_k_set)
                if overlap >= cutoff:
                    mod_num = mods[mod][0]
                    bgcs = sorted(mods2bgc[mod_num])
                    overl_genes = sorted(dom_overlap_set)
                    non_overl_genes = sorted(mod_set|doms_k_set)
                    for bgc in bgcs:
                        bgc_class = bgc_classes.get(bgc,['None'])[0]
                        match = [overlap,l_overlap,bgc,bgc_class,mod_num,\
                            overl_genes,non_overl_genes]
                        match_dict[name_k].append(match)
            print(match_dict)
            break

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
    # process_lda(lda, lda_dict, bow_corpus, modules, cmd.feat_num, bgcs,
        # cmd.min_feat_score, bgclist, cmd.out_folder, bgc_classes_dict, \
        # amplif=cmd.amplify, plot=cmd.plot, known_subcl=known_subclusters)
    
    find_stat_method_overlap(modules, known_subclusters, bgc_classes_dict,\
        0.4, mod_nums2bgc, cmd.out_folder)

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
