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
from copy import deepcopy
from functools import partial
from glob import glob, iglob
from itertools import combinations, product
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count, Manager
import networkx as nx
from operator import itemgetter
import os
import random
from statsmodels.stats.multitest import multipletests
import subprocess
from sympy import binomial as ncr
import time

def get_commands():
    parser = argparse.ArgumentParser(description="Detects sub-clusters in \
        BGCs represented as strings of Pfams according to Del Carratore et. \
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
        analysis (default: 2)", type=int)
    parser.add_argument("--sim_cutoff", dest="sim_cutoff", default=0.95,
        help="Cutoff for cluster similarity in redundancy filtering (default:\
        0.95)", type=float)
    parser.add_argument("-p", "--pval_cutoff", dest="pval_cutoff", \
        default = 0.1, type=float, help='P-value cutoff for determining a \
            significant interaction in module detection (default: 0.1)')
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

def remove_infr_doms(clusdict, m_doms, verbose):
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
    print('  {} domains are removed, {} domains are left'.format(\
        len(deldoms),len(domcounter.keys())-len(deldoms)))
    clus_no_deldoms = {}
    for k,v in clusdict.items():
        newv = ['-' if dom in deldoms else dom for dom in v]
        doml = len({v for v in newv if not v == '-'})
        if doml > 1:
            clus_no_deldoms[k] = newv
        else:
            if verbose:
                print('  {} removed as it has less than 2 domains'.format(k))
    print(' {} clusters have less than {} domains and are excluded'.format(\
        len(clusdict.keys()) - len(clus_no_deldoms), m_doms))
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

def calc_adj_pval_wrapper(count_dict, clusdict, cores, verbose):
    '''Returns list of tuples of corrected pvals for each domain pair

    counts: nested dict { dom1:{ count:x,N1:y,N2:z,B1:{dom2:v},B2:{dom2:w} } }
    clusdict: dict of {cluster:[domains]}
    cores: int, amount of cores to use
    verbose: bool, if True print additional information
    '''
    #NB. Nall for coloc_pval should be len(remove_dupl_doms(values))
    print('\nCalculating adjacency pvalues')
    N = sum([len(values) for values in clusdict.values()])
    pool = Pool(cores, maxtasksperchild=20)
    pvals_ori = pool.map(partial(calc_adj_pval, counts=count_dict, Nall=N), \
        count_dict.items())
    #remove Nones, unlist and sort
    pvals_ori = [lst for lst in pvals_ori if lst]
    pvals_ori = sorted([tup for lst in pvals_ori for tup in lst])
    #to check if there are indeed 2 pvalues for each combination
    check_ps = [(tup[0],tup[1]) for tup in pvals_ori]
    check_c = Counter(check_ps)
    pvals = [p for p in pvals_ori if check_c[(p[0],p[1])] == 2]
    if not len(pvals) == len(pvals_ori):
        if verbose:
            p_excl = [p for p in pvals if check_c[(p[0],p[1])] != 2]
            print('  error with domain pairs {}'.format(', '.join(p_excl)))
            print('  these are excluded')
    #Benjamini-Yekutieli multiple testing correction
    pvals_adj = multipletests(list(zip(*pvals))[2], method='fdr_by')[1]
    #adding adjusted pvals and choosing max
    ptups = []
    for ab1, ab2, p1, p2 in \
        zip(pvals[::2], pvals[1::2], pvals_adj[::2], pvals_adj[1::2]):
        assert(ab1[0]==ab2[0] and ab1[1]==ab2[1])
        pmax = max([p1,p2])
        ptups.append(((ab1[0],ab1[1]),pmax))
    return ptups

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

def calc_coloc_pval(domval_pair, counts, Nall):
    '''Returns a list of sorted tuples (domA,domB,pval)
    
    domval_pair: tuple of (domA, { count:x,N1:y,B1:{dom2:v,dom3:w } })
    counts: nested dict { domA:{ count:x,N1:y,B1:{dom2:v,dom3:w } } }
    Nall: int, all possible positions in all clusters
    '''
    domA, vals = domval_pair
    #domains without interactions do not end up in pvals
    if not vals['B1']:
        return
    pvals = []
    count = vals['count']
    Ntot = Nall - count
    N1 = vals['N1']
    N0 = Ntot - N1
    interactions = vals['B1'].keys()
    for domB in interactions:
        B1 = vals['B1'][domB]
        Btot = counts[domB]['count']
        pval = float(1 - sum([ncr(N0,(Btot-d)) * ncr(N1, d) \
            for d in range(B1)]) / ncr(Ntot,Btot))
        ab_int = sorted((domA,domB))
        pvals.append((ab_int[0],ab_int[1],pval))
    return pvals

def calc_coloc_pval_wrapper(count_dict, clusdict, cores, verbose):
    '''Returns list of tuples of corrected pvals for each domain pair

    counts: nested dict { domA:{ count:x,N1:y,B1:{dom2:v,dom3:w } } }
    clusdict: dict of {cluster:[domains]}
    cores: int, amount of cores to use
    verbose: bool, if True print additional information
    '''
    #NB. Nall for coloc_pval should be len(remove_dupl_doms(values))
    print('\nCalculating colocalisation pvalues')
    N = sum([len(remove_dupl_doms(values)) for values in clusdict.values()])
    pool = Pool(cores, maxtasksperchild=50)
    pvals_ori = pool.map(partial(calc_coloc_pval, counts=count_dict, Nall=N), \
        count_dict.items())
    #remove Nones, unlist and sort
    pvals_ori = [lst for lst in pvals_ori if lst]
    pvals_ori = sorted([tup for lst in pvals_ori for tup in lst])
    #to check if there are indeed 2 pvalues for each combination
    check_ps = [(tup[0],tup[1]) for tup in pvals_ori]
    check_c = Counter(check_ps)
    pvals = [p for p in pvals_ori if check_c[(p[0],p[1])] == 2]
    if not len(pvals) == len(pvals_ori):
        if verbose:
            p_excl = [p for p in pvals if check_c[(p[0],p[1])] != 2]
            print('  error with domain pairs {}'.format(', '.join(p_excl)))
            print('  these are excluded')
    #Benjamini-Yekutieli multiple testing correction
    pvals_adj = multipletests(list(zip(*pvals))[2], method='fdr_by')[1]
    #adding adjusted pvals and choosing max
    ptups = []
    for ab1, ab2, p1, p2 in \
        zip(pvals[::2], pvals[1::2], pvals_adj[::2], pvals_adj[1::2]):
        assert(ab1[0]==ab2[0] and ab1[1]==ab2[1])
        pmax = max([p1,p2])
        ptups.append(((ab1[0],ab1[1]),pmax))
    return ptups

def keep_lowest_pval(colocs, adjs):
    '''
    Returns all domain pairs with their lowest pvalue as an edge for n

    colocs, adjs: list of tuples [((dom1,dom2),pval)]
    Tuples look like (dom1,dom2,{pval:x})
    '''
    pvals = colocs+adjs
    counter = Counter(list(zip(*pvals))[0])
    dupl = sorted([tup for tup in pvals if counter[tup[0]] == 2])
    uniques = [tup for tup in pvals if counter[tup[0]] == 1]
    lowest = []
    for p1,p2 in zip(dupl[::2],dupl[1::2]):
        pmin = min([p1[1],p2[1]])
        lowest.append((p1[0][0],p1[0][1],{'pval':pmin}))
    uniques = [(tup[0][0],tup[0][1],{'pval':tup[1]}) for tup in uniques]
    return lowest+uniques

def generate_graph(edges, verbose):
    '''Returns a networkx graph

    edges: list of tuples, (pair1,pair2,{attributes})
    '''
    g = nx.Graph()
    g.add_edges_from(edges)
    if verbose:
        print('\nGenerated graph with:')
        print(' {} nodes'.format(g.number_of_nodes()))
        print(' {} edges'.format(g.number_of_edges()))
    return g

def visualise_graph(graph, subgraph_list = None, groups = True):
    '''Plots a graph with possible subgraphs in different colours

    graph: networkx graph
    subgraph_list: list of lists of node names that should be coloured
        differently, default = None
    groups: bool, are there groups in the subgraph_list that you want to
        colour differently (True)? or are nodes in subgraph list one seperate
        group (False)
    '''
    cols = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    options = {'node_size': 2,'width': 0.2}
    pos = nx.spring_layout(graph)
    plt.figure()
    nx.draw_networkx(graph, pos=pos, with_labels = False, node_color='black',\
        **options)
    if subgraph_list:
        if groups:
            for sub in subgraph_list:
                nx.draw_networkx_nodes(graph, pos=pos, nodelist=sub, \
                    node_color=random.choice(cols), **options)
        else:
            nx.draw_networkx_nodes(graph, pos=pos, nodelist=subgraph_list, \
                    node_color='#91bfdb', marker='s', **options)
    plt.show()

def generate_modules_wrapper(pval_edges, sign_cutoff, outfile, cores, \
    verbose):
    '''
    Returns a dict with all modules {(module):strictest_pval_cutoff}

    pval_edges: list of tuples, [(dom1,dom2,pval)]
    sign_cutoff: float, pvalue cutoff
    cores: int, number of cores to use
    verbose: bool, if True print additional information
    '''
    print('\nFinding all modules with a pvalue lower than {}'.format(\
        sign_cutoff))
    sign_pvs = [ptup for ptup in pval_edges if ptup[2]['pval'] <= sign_cutoff]
    print('  {} significant domain pair interactions'.format(len(sign_pvs)))
    pv_values = {pv['pval'] for pv in list(zip(*sign_pvs))[2]}
    #watch out if pv_values gets really big, maybe get 100,000 fixed numbers
    #to loop over
    print('  looping through {} pvalue cutoffs'.format(len(pv_values)))
    pool = Pool(cores, maxtasksperchild = 50)
    modules = pool.map(partial(generate_modules, dom_pairs=sign_pvs), \
        pv_values)
    modules_dict = {}
    for p_mods_pair in modules:
        p = list(p_mods_pair)[0]
        mods = p_mods_pair[p]
        for mod in mods:
            try:
                prev_val = modules_dict[mod]
            except KeyError:
                modules_dict[mod] = p
            else:
                if p < prev_val:
                    modules_dict[mod] = p
    print('{} modules detected'.format(len(modules_dict)))
    with open(outfile, 'w') as out:
        #sort them first, and maybe print more info length for example
        for mod,p in modules_dict.items():
            out.write('{}\t{}\n'.format(','.join(mod),p))
    return modules_dict

def generate_modules(sign_cutoff, dom_pairs):
    '''Updates main_dict with all modules and the lowest cutoff to detect them
    '''
    edges = (edge for edge in dom_pairs if edge[2]['pval'] <= sign_cutoff)
    mod_graph = generate_graph(edges, False)
    cliqs = nx.algorithms.clique.find_cliques(mod_graph)
    cliqs = {tuple(sorted(clq)) for clq in cliqs if len(clq) > 2}
    return {sign_cutoff:cliqs}

if __name__ == "__main__":
    cmd = get_commands()
    #adjust later
    infile = os.path.join(cmd.out_folder, 'testdata_filtered_clusterfile.csv')

    f_clus_dict = read_clusterfile(infile, cmd.min_doms, \
        cmd.verbose)[0]
    f_clus_dict_rem = remove_infr_doms(f_clus_dict, cmd.min_doms, cmd.verbose)
    clusters = list(f_clus_dict_rem.keys())
    adj_counts, c_counts = count_interactions(f_clus_dict_rem, cmd.verbose)
    adj_pvals = calc_adj_pval_wrapper(adj_counts, f_clus_dict_rem, cmd.cores,\
        cmd.verbose)
    col_pvals = calc_coloc_pval_wrapper(c_counts, f_clus_dict_rem, cmd.cores,\
        cmd.verbose)

    pvals = keep_lowest_pval(col_pvals,adj_pvals)

    start = time.time()
    out_file = os.path.join(cmd.out_folder, 'testdata_modules.txt')
    mods = generate_modules_wrapper(pvals, cmd.pval_cutoff, out_file, cmd.cores,\
        cmd.verbose)
    end = time.time()
    print('\nFound modules in {} seconds'.format(end-start))
    
