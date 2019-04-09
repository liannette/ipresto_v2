#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to convert domainlist of bgcs into graphs based on a Adjacency index
and filter out redundant bgcs.

Usage:
python3 bgc_to_pfam.py -h

Example usage:
python3 filter_bgcs.py

Notes:

Layout:
get_commands
read_clusterfile
calc_adj_index

TODO:
make functions that:
-read clusterfile and write a new filtered file without clusters with less
    domains than min_doms
-calculate Adjacency dist between two domains/see if one is contained within
    the other
-loop through every combination of domains and return dist in a matrix
-create graphs from the matrix using some cutoff (using igraph, networkx?)
-choose representatives
-write file with just representatives and one with representatives linked
    to the bgcs that are filtered out

Required:
networkx (https://github.com/networkx/networkx)
'''
import os
from glob import glob, iglob
import subprocess
from collections import OrderedDict, Counter
import argparse
from multiprocessing import Pool, cpu_count
from functools import partial
from itertools import combinations
import networkx as nx
import matplotlib.pyplot as plt
import random

def get_commands():
    parser = argparse.ArgumentParser(description="A script to turn a \
        clusterfile with domains into graphs and filter out redundant bgcs")
    parser.add_argument("-i", "--in_file", dest="in_file", help="Input \
        file", required=True)
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
    print("\nFiltering clusterfile")
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
    print("Done. Keeping {} clusters".format(len(clus_dict)))
    print(" {} clusters have less than {} domains".format(filtered,m_doms))
    return clus_dict, len_dict

def calc_adj_index(clus1, clus2):
    '''Returns the adjacency index between two clusters

    clus1, clus2: list of strings, domainlist of a cluster

    If there is an empty gene between two domains these two domains are not
        adjacent
    '''
    #generate all unique domain pairs
    dom_p1 = {tuple(sorted(dp)) for dp in zip(*(clus1[:-1],clus1[1:])) \
        if not '-' in dp}
    dom_p2 = {tuple(sorted(dp)) for dp in zip(*(clus2[:-1],clus2[1:])) \
        if not '-' in dp}
    #if doms are separated by '-' then there are no dom pairs. if happens ai=0
    if not dom_p1 or not dom_p2:
        return 0.0        
    ai = len(dom_p1 & dom_p2)/len(dom_p1 | dom_p2)
    return ai

def is_contained(clus1, clus2):
    '''
    Returns a bool if all domains from one of the clusters are in the other

    clus1, clus2: list of strings, domainlist of a cluster
    '''
    one_in_two = all([dom in clus2 for dom in clus1 if not dom == '-'])
    two_in_one = all([dom in clus1 for dom in clus2 if not dom == '-'])
    if one_in_two or two_in_one:
        return True
    return False

def generate_edges(nodes, dom_dict, cutoff, cores):
    '''Returns a pair of clusters in a tuple if ai/contained above cutoff

    nodes: list of strings, clusternames
    dom_dict: dict {clus1:[domains]}, clusters linked to domains
    cutoff: float, between 0-1, when clusters are similar
    cores: int, amount of cores used for calculation
    '''
    print("\nGenerating similarity scores")
    pairs = combinations(clus_names, 2)
    pool = Pool(cores, maxtasksperchild = 100)
    #I could add imap if this is still too slow for antismashdb
    #with imap specify chunksize as number of pairs?
    edges = pool.map(partial(generate_edge, d_dict = dom_dict, \
        cutoff = cutoff), pairs)
    edges = [edge for edge in edges if not edge == None]
    print("Done. {} pairs above threshold".format(len(edges)))
    return edges

def generate_edge(pair, d_dict, cutoff):
    '''
    Calculate similarity scores between two bgcs and return if above cutoff

    pair: tuple of 2 strings, 2 clusternames
    d_dict: dict of {clustername:domains}
    cutoff: float
    A tuple is returned that can be read as an edge by nx.Graph.add_edges_from
    '''
    p1,p2 = pair
    contained = is_contained(d_dict[p1], d_dict[p2])
    ai = calc_adj_index(d_dict[p1],d_dict[p2])
    if ai == None:
        print(ai,pair)
    if contained or ai > cutoff:
        # print(pair,ai,contained)
        return(p1,p2,{'ai':ai,'contained':contained})

def generate_graph(edges):
    '''Returns a networkx graph

    edges: list of tuples, (pair1,pair2,{attributes})
    '''
    g = nx.Graph()
    g.add_edges_from(edges)
    print('Generated graph with:')
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
            for sub in subgraph_list:
                nx.draw_networkx_nodes(graph, pos=pos, nodelist=sub, \
                    node_color='#91bfdb', marker='s', **options)
    plt.show()

def find_representatives(clqs, d_l_dict, graph):
    '''
    Returns {representative:[clique]} based on bgc with most domains in clique

    clqs: list of lists of strings, cliques of clusters
    d_l_dict: dict of {clus_name:amount_of_domains(int)}
    graph: networkx graph structure of the cliques
    Returns also the representative in the clique list (dict values)
    The longest cluster is chosen (most domains). If there are multiple
        longest clusters then the cluster with the least connections is
        chosen (to preserve most information).
    '''
    random.seed(0)
    reps_dict = OrderedDict()
    dels = set() #set of nodes for which a representative has been found
    clqs = sorted(clqs, key=len, reverse=True)
    for cliq in clqs:
        cliq = [clus for clus in cliq if not clus in dels]
        if cliq:
            domlist = [(clus,d_l_dict[clus]) for clus in cliq]
            maxdoml = max([doms[1] for doms in domlist])
            clus_maxlen = [clus for clus, doml in domlist \
                if doml == maxdoml]
            if len(clus_maxlen) > 1:
                min_degr = min([deg for clus, deg in graph.degree(clus_maxlen)])
                rep = random.choice([clus for clus in clus_maxlen \
                    if graph.degree(clus) == min_degr])
            else:
                rep = clus_maxlen[0]
            try:
                reps_dict[rep].update(cliq)
            except KeyError:
                reps_dict[rep] = set(cliq)
            cliq.remove(rep)
            dels.update(cliq)
    return reps_dict

def find_all_representatives(d_l_dict, g):
    '''Iterates find_representatives until there are no similar bgcs

    d_l_dict: dict of {clus_name:amount_of_domains(int)}
    graph: networkx graph structure containing the cliques
    all_reps_dict: dict of {representative:[represented]}
    '''
    all_reps_dict = {}
    subg = g.subgraph(g.nodes)
    #make a while loop after testing while subg.number_of_edges != 0
    
    for i in range(5):
        print('iteration {}, edges left: {}'.format(i,subg.number_of_edges()))
        cliqs = nx.algorithms.clique.find_cliques(subg)
        cliqs = sorted(cliqs, key=len)
        print(len(cliqs))
        cliqs = [cl for cl in cliqs if len(cl) > 1]
        print(len(cliqs))
        reps_dict = find_representatives(cliqs, d_l_dict, subg)
        subg = subg.subgraph(reps_dict.keys())
        #merge reps_dict with all_reps_dict
        for key, vals in reps_dict.items():
            if not key in all_reps_dict:
                all_reps_dict[key] = vals
            else:
                #merge represented clusters in a new representative
                newvals = [v for val in vals for v in all_reps_dict[val] \
                    if val in all_reps_dict.keys()]
                #if statement for bgcs already represented by this 
                #representative and thus no longer in all_reps_dict
                for remove_key in vals:
                    del all_reps_dict[remove_key]
                all_reps_dict[key] = newvals
    return all_reps_dict

if __name__ == "__main__":
    cmd = get_commands()
    random.seed(0)
    #read_clusterfile filters out bgcs with less than min_dom domains
    dom_dict, doml_dict = read_clusterfile(cmd.in_file, cmd.min_doms, \
        cmd.verbose)
    clus_names = list(dom_dict.keys())#[0:300] #make list so bgcs have an index
    edges = generate_edges(clus_names, dom_dict, cmd.sim_cutoff, cmd.cores)
    graph = generate_graph(edges)
    cliqs = nx.algorithms.clique.find_cliques(graph)
    cliqs = sorted(cliqs, key=len)
    reps_dict = find_representatives(cliqs, doml_dict, graph)
    print(len(reps_dict.keys()))
    reps = [[rep] for rep in reps_dict.keys()]
    # visualise_graph(graph,reps)
    rep_g = graph.subgraph(reps_dict.keys())
    r_cl = nx.algorithms.clique.find_cliques(rep_g)
    r_cl = [cl for cl in r_cl if len(cl) > 1]
    r_cl = sorted(r_cl,key=len)
    print([len(cl) for cl in r_cl],len(r_cl))
    # visualise_graph(rep_g,r_cl)
    reps_reps_dict = find_representatives(r_cl, doml_dict, graph)
    print(reps_reps_dict)
    # print('ARGH01000000_KB894962.1.cluster038:',graph.adj['ARGH01000000_KB894962.1.cluster038'])
    all_reps = find_all_representatives(doml_dict, graph)
    maxrepr = max([len(all_reps[key]) for key in all_reps.keys()])
    rep_max = [rep for rep in all_reps.keys() if len(all_reps[rep]) == maxrepr][0]
    print(rep_max,doml_dict[rep_max],all_reps[rep_max],maxrepr)
    all_reps1 = [[rep] for rep in all_reps.keys()]
    visualise_graph(graph, all_reps1, groups = False)
    #choose from each cliq the cluster with the most domains as representative
    # reps = []
    # for cliq in cliqs:
        # domlist = [(clus,doml_dict[clus]) for clus in cliq]
        # # domlens = [doml_dict[clus] for clus in cliq]
        # maxdoml = max([doms[1] for doms in domlist])
        # clus_maxlen = random.choice([clus for clus, doml in domlist \
            # if doml == maxdoml])
        # reps.append(clus_maxlen)
        # # print(maxdoml, clus_maxlen, domlist)
    # dupl_reps = set(rep for rep in reps if reps.count(rep) > 1)
    # print(len(set(reps)), len(dupl_reps), len(reps))
    # reps = list(set(reps))
    # # reps_edges = generate_edges(reps,dom_dict,cmd.sim_cutoff,cmd.cores)
    # # reps_graph = generate_graph(reps_edges)
    # reps_graph = graph.subgraph(reps)
    # print("nodes, edges of reps_graph: {}, {}".format(reps_graph.number_of_nodes(),reps_graph.number_of_edges()))
    # reps_cliqs = nx.algorithms.clique.find_cliques(reps_graph)
    # visualising, maybe with subplot plot all different cliques/networks
    
    # visualise_graph(reps_graph, reps_cliqs)
    

