#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to detect sub-clusters in bgcs represented as strings of Pfams.

Usage:
python3 module_detection.py -h

Example usage:
python3 module_detection.py -i ../testdata -o ../testdata_domains -c 16
Notes:

Layout:
get_commands

Required:
python 3.6
Biopython
networkx

To do:
preprocess clusters (remove duplicates pfams, remove domains that only occur
    once or twice (less than thrice), remove bgcs that have less than 2 
    domains as a result)
    in the md script they remove duplicate pfams, replace them with - and
    leave one copy at the end
generate dom interactions in a matrix (adj and coloc) (using dom_list as index
    and np arrays?)
computing p-values by randomly distributing data
multiple testing corrections on (BY)
graph generation with sign interactions and iteration
'''
import argparse
from Bio import SeqIO
from Bio import SearchIO
from collections import OrderedDict, Counter, defaultdict
from functools import partial
from glob import glob, iglob
from itertools import combinations
from multiprocessing import Pool, cpu_count
#import networkx as nx
import os
import random
import subprocess

def get_commands():
    parser = argparse.ArgumentParser(description="Detects sub-clusters in \
        bgcs represented as strings of Pfams according to Del Carratore et. \
        al (2019)")
    parser.add_argument("-i", "--in_folder", dest="in_folder", help="Input \
        directory of gbk files", required=True)
    parser.add_argument("-o", "--out_folder", dest="out_folder", 
        required=True, help="Output directory, this will contain all output \
        data files.")
    parser.add_argument("-c", "--cores", dest="cores", default=cpu_count(), 
        help="Set the number of cores the script may use (default: use all \
        available cores)", type=int)
    parser.add_argument("-v", "--verbose", dest="verbose", required=False,
        action="store_true", default=False, help="Prints more detailed \
        information.")
    parser.add_argument("--min_doms", dest="min_doms", default=2,
        help="The minimum amount of domains in a BGC to be included in the \
        analysis. Default is 2 domains", type=int)
    parser.add_argument("--sim_cutoff", dest="sim_cutoff", default=0.95,
        help="Cutoff for cluster similarity in redundancy filtering (default:\
        0.95)", type=float)
    return parser.parse_args()

def read_clusterfile(infile, m_doms, verbose):
    """Reads a clusterfile into a dictionary

    infile: str, filepath
    m_doms: int, minimum of domains a cluster should have
    verbose: bool, if True print additional info
    clusters with less than m_doms domains are not returned
    """
    print("\nReading {}".format(infile))
    filtered = 0
    with open(infile, 'r') as inf:
        clus_dict = OrderedDict()
        len_dict = OrderedDict()
        for line in inf:
            line = line.strip().split(',')
            clus = line[0]
            doms = line[1:]
            ldoms = len([dom for dom in doms if not dom == '-'])
            if ldoms < m_doms:
                filtered +=1
                if verbose:
                    print("  excluding {} less than min domains".format(clus))
                continue
            if not clus in clus_dict.keys():
                clus_dict[clus] = doms
                len_dict[clus] = ldoms
            else:
                print("Clusternames not unique, {} read twice".format(clus))
    print("Done. Read {} clusters".format(len(clus_dict)))
    print(" {} clusters have less than {} domains and are excluded".format(\
        filtered,m_doms))
    return clus_dict, len_dict

def remove_infr_doms(clusdict, verbose):
    '''Returns clusdict with domains replaced  with - if they occur < 3

    clusdict: dict of {cluster:[domains]}
    verbose: bool, if True print additional info

    Deletes clusters with 1 unique dom
    '''
    print('\nRemoving domains that occur less than 3 times')
    domcounter = Counter()
    domcounter.update([v for vals in clusdict.values() for v in vals \
        if not v == '-'])
    deldoms = [key for key in domcounter if domcounter[key] <= 2]
    clus_no_deldoms = {}
    for k,v in clusdict.items():
        newv = ['-' if dom in deldoms else dom for dom in v]
        doml = len({v for v in newv if not v == '-'})
        if doml > 1:
            clus_no_deldoms[k] = newv
        else:
            if verbose:
                print('  {} removed as it has less than 2 domains'.format(k))
    print('Continuing with {} bgcs that have more than 1 domain left'.format(\
        len(clus_no_deldoms)))
    return clus_no_deldoms

def remove_dupl_doms(clusdict):
    '''
    Replaces duplicate domains in a cluster with '-', writes domain at the end

    clusdict: dict of {cluster:[domains]}
    '''
    clus_no_dupl = {}
    for k,v in clusdict.items():
        domc = Counter(v)
        dupl = [dom for dom in domc if domc[dom] > 1 if not dom == '-']
        if dupl:
            newv = ['-' if dom in dupl else dom for dom in v]
            for dom in dupl:
                newv += ['-',dom]
        else:
            newv = v
        clus_no_dupl[k] = newv
    return clus_no_dupl

def count_adj(counts, cluster, verbose):
    '''Counts all adjacency interactions between domains in a cluster

    cluster: list of strings, domains
    verbose: bool, if True print additional info
    '''
    #for multiprocessing, from concurrent.futures ProcessPoolExecutor?
    for i, dom in enumerate(cluster):
        if not dom == '-':
            if i == 0:
                edge = 1
                adj = [cluster[1]]
            elif i == len(cluster)-1:
                edge = 1
                adj = [cluster[i-1]]
            else:
                edge = 2
                adj = [cluster[i-1],cluster[i+1]]
            print('Domain: {} , Edges: {}, Adj:{}'.format(dom,edge,', '.join(adj)))
            counts[dom]['count'] += 1
            counts[dom]['edge'] += edge
            for ad in adj:
                if not ad == '-':
                    try:
                        counts[dom][ad] += 1
                    except TypeError:
                        counts[dom][ad] = 1
    print(cluster)

def makehash():
    '''Function to initialise nested dict
    '''
    return defaultdict(makehash)

def count_interactions(clusdict, verbose):
    '''
    '''
    all_doms = {v for val in clusdict.values() for v in val}
    all_doms.remove('-')
    dcounts = makehash()
    for d in all_doms:
        for v in ['count','edge']:
            dcounts[d][v] = 0
    print(dcounts)
    for clus in clusdict.values():
        count_adj(dcounts, clus, verbose)
    print(dcounts)

if __name__ == "__main__":
    cmd = get_commands()
    #adjust later
    infile = os.path.join(cmd.out_folder, 'testdata_filtered_clusterfile.csv')

    f_clus_dict, f_doml_dict = read_clusterfile(infile, cmd.min_doms, \
        cmd.verbose)
    f_clus_dict_rem1 = remove_infr_doms(f_clus_dict, cmd.verbose)
    clusters = list(f_clus_dict_rem1.keys())
    # for i in range(10):
        # subdict = {clusters[i] : f_clus_dict_rem1[clusters[i]]}
        # print(remove_dupl_doms(subdict))
    x=5
    testdict = {clus : f_clus_dict_rem1[clus] for clus in clusters[:x]}
    # testclus = f_clus_dict_rem1[clusters[0]]
    # def_val = ['count', 'edge']
    # t=set(testclus)
    # t.remove('-')
    # testclus_c = makehash()
    # for c in t:
        # for v in def_val:
            # testclus_c[c][v] = 0
    # print(testclus_c)
    # count_adj(testclus_c, testclus,cmd.verbose)
    # print(testclus_c)
    count_interactions(testdict, cmd.verbose)
    
