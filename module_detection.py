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
read_clusterfile
remove_infr_doms

Required:
python 3.6
Biopython
networkx
scipy

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
from itertools import combinations, product
from multiprocessing import Pool, cpu_count
import networkx as nx
import os
import random
from statsmodels.stats.multitest import multipletests
import subprocess
from sympy import binomial as ncr

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

def remove_dupl_doms(cluster):
    '''
    Replaces duplicate domains in a cluster with '-', writes domain at the end

    cluster: list of strings, domain list
    '''
    domc = Counter(cluster)
    dupl = [dom for dom in domc if domc[dom] > 1 if not dom == '-']
    if dupl:
        newclus = ['-' if dom in dupl else dom for dom in cluster]
        for dom in dupl:
            newclus += ['-',dom]
    else:
        newclus = cluster
    return newclus

def count_adj(counts, cluster):
    '''Counts all adjacency interactions between domains in a cluster

    counts: nested dict { dom1:{ count:x,N1:y,N2:z,B1:{dom2:v},B2:{dom2:w} } }
    cluster: list of strings, domains
    '''
    for i, dom in enumerate(cluster):
        if i == 0:
            edge = 1
            adj = [cluster[1]]
        elif i == len(cluster)-1:
            edge = 1
            adj = [cluster[i-1]]
        else:
            edge = 2
            adj = [cluster[i-1],cluster[i+1]]
            if adj[0] == adj[1] and adj[0] != '-':
                #B2 and N2 counts
                prevdom = cluster[i-1]
                counts[prevdom]['N1'] -= 2
                counts[prevdom]['N2'] += 1
                if dom != '-' and dom != prevdom:
                    counts[prevdom]['B1'][dom] -= 2
                    try:
                        counts[prevdom]['B2'][dom] += 1
                    except TypeError:
                        counts[prevdom]['B2'][dom] = 1
        if not dom == '-':
            counts[dom]['count'] += 1
            counts[dom]['N1'] += edge
            for ad in adj:
                if ad != '-' and ad != dom:
                    try:
                        counts[dom]['B1'][ad] += 1
                    except TypeError:
                        counts[dom]['B1'][ad] = 1

def count_coloc(counts, cluster):
    '''Counts all colocalisation interactions between domains in a cluster

    counts: nested dict { dom1:{ count:x,N1:y,B1:{dom2:v,dom3:w } } }
    cluster: list of strings, domains
    verbose: bool, if True print additional info
    '''
    N1 = len(cluster)-1
    for dom in cluster:
        if not dom == '-':
            counts[dom]['count'] += 1
            counts[dom]['N1'] += N1
            coloc = set(cluster)
            try:
                coloc.remove('-')
            except KeyError:
                pass
            coloc.remove(dom)
            for colo in coloc:
                try:
                    counts[dom]['B1'][colo] += 1
                except TypeError:
                    counts[dom]['B1'][colo] = 1

def makehash():
    '''Function to initialise nested dict
    '''
    return defaultdict(makehash)

def count_interactions(clusdict, verbose):
    '''Count all adj and coloc interactions between all domains in clusdict

    clusdict: dict of {cluster:[domains]}
    verbose: bool, if True print additional info
    Returns two dicts, one dict with adj counts and one with coloc counts
    adj counts:
        { dom1:{ count:x,N1:y,N2:z,B1:{dom2:v},B2:{dom2:w} } }
    coloc counts:
        { dom1:{ count:x,N1:y,B1:{dom2:v,dom3:w } } }
    '''
    print('\nCounting colocalisation and adjacency interactions')
    all_doms = {v for val in clusdict.values() for v in val}
    all_doms.remove('-')
    #initialising count dicts
    adj_counts = makehash()
    for d in all_doms:
        for v in ['count','N1','N2']:
            adj_counts[d][v] = 0
        for w in ['B1','B2']:
            adj_counts[d][w] = makehash()
        #N1: positions adj to one domA, N2: positions adj to two domA
        #B1: amount of domB adj to one domA, B2: positions adj to two domA

    coloc_counts = makehash()
    for d in all_doms:
        for v in ['count','N1']:
            coloc_counts[d][v] = 0
        coloc_counts[d]['B1'] = makehash()
        #N1: all possible coloc positions in a cluster, cluster lenght - 1
        #B1: amount of domB coloc with domA

    for clus in clusdict.values():
        count_adj(adj_counts, clus)
        filt_clus = remove_dupl_doms(clus)
        count_coloc(coloc_counts, filt_clus)
    return(adj_counts, coloc_counts)

def calc_adj_pval_wrapper(count_dict, clusdict, cores):
    '''Returns list of tuples of corrected pvals for each domain pair

    counts: nested dict { dom1:{ count:x,N1:y,N2:z,B1:{dom2:v},B2:{dom2:w} } }
    clusdict: dict of {cluster:[domains]}
    '''
    #NB. Nall for coloc_pval should be len(remove_dupl_doms(values))
    print('\nCalculating adjacency pvalues')
    N = sum([len(values) for values in clusdict.values()])
    pool = Pool(cores, maxtasksperchild=20)
    pvals = pool.map(partial(calc_adj_pval, counts=count_dict, Nall=N), \
        count_dict.items())
    #remove Nones and unlist
    pvals = [lst for lst in pvals if lst]
    pvals = sorted([tup for lst in pvals for tup in lst])
    #Benjamini-Yekutieli correction
    pvals_adj = multipletests(list(zip(*pvals))[2], method='fdr_by')
    #adding adjusted pvals and choosing max
    print(pvals)
    return pvals

def calc_adj_pval(domval_pair, counts, Nall):
    '''Returns a list of sorted tuples (domA,domB,pval)
    
    domA: string of domain name
    vals: dict of domA interaction info
        {count:x,N1:y,N2:z,B1:{dom2:v},B2:{dom2:w}}
    counts: nested dict { domA:{ count:x,N1:y,N2:z,B1:{dom2:v},B2:{dom2:w} } }
    '''
    domA, vals = domval_pair
    #domains without interactions do not end up in pvals
    if not vals['B1'] and not vals['B2']:
        return
    pvals = []
    count = vals['count']
    Ntot = Nall - count
    N1 = vals['N1']
    N2 = vals['N2']
    N0 = Ntot - N1 - N2
    interactions = vals['B1'].keys() | vals['B2'].keys()
    for domB in interactions:
        if domB not in vals['B2']:
            B1 = vals['B1'][domB]
            Btot = counts[domB]['count']
            pval = float(1 - sum([ncr(N0,(Btot-d)) * ncr(N1, d) \
                for d in range(B1)]) / ncr(Ntot,Btot))
        elif vals['B1'][domB] == 0:
            B2 = vals['B2'][domB]
            Btot = counts[domB]['count']
            pval = float(1 - sum([ncr(N0,(Btot-d)) * ncr(N2, d) \
                for d in range(B2)]) / ncr(Ntot,Btot))
        else:
            B1 = vals['B1'][domB]
            B2 = vals['B2'][domB]
            Btot = counts[domB]['count']
            pval = float(\
                1 - sum([ncr(N0,Btot-d1-d2) * ncr(N1,d1) * ncr(N2,d2) \
                for d1,d2 in product(range(B1+1),range(B2+1)) \
                if d1+d2 != B1+B2]) / ncr(Ntot,Btot))
        ab_int = sorted((domA,domB))
        pvals.append((ab_int[0],ab_int[1],pval))
    return pvals

if __name__ == "__main__":
    cmd = get_commands()
    #adjust later
    infile = os.path.join(cmd.out_folder, 'testdata_filtered_clusterfile.csv')

    f_clus_dict = read_clusterfile(infile, cmd.min_doms, \
        cmd.verbose)[0]
    f_clus_dict_rem = remove_infr_doms(f_clus_dict, cmd.verbose)
    clusters = list(f_clus_dict_rem.keys())
    # for i in range(10):
        # subdict = {clusters[i] : f_clus_dict_rem1[clusters[i]]}
        # print(remove_dupl_doms(subdict))
    x=10
    testdict = {clus : f_clus_dict_rem[clus] for clus in clusters[:x]}

    adj_counts, c_counts = count_interactions(f_clus_dict_rem, cmd.verbose)
    adj_pvals = calc_adj_pval_wrapper(adj_counts, f_clus_dict_rem, cmd.cores)
    # print(adj_counts)
    #print(adj_counts['RCC1'])
    #print(adj_counts['Big_3_5'])
    # adj_counts, c_counts = count_interactions(testdict, cmd.verbose)
    # print('\n',sorted([(key,c_counts[key]['count'],c_counts[key]['N1']) for key in c_counts],key=lambda x: x[1]))
    # print([(key,c_counts[key]['count'],c_counts[key]['N1']) for key in c_counts \
        # if c_counts[key]['count'] > len(f_clus_dict_rem.keys())])
