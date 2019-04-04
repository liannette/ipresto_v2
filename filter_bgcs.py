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
from random import choice

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
    parser.add_argument("--min_doms", dest="min_doms", default=5,
        help="The minimum amount of domains in a BGC to be included in the \
        analysis. Default is 0 domains", type=int)
    parser.add_argument("--sim_cutoff", dest="sim_cutoff", default=0.95,
        help="Cutoff for cluster similarity in redundancy filtering (default:\
        0.95)", type=float)
    return parser.parse_args()

def read_clusterfile(infile, m_doms, verbose):
    """Reads a clusterfile into a dictionary

    infile: str, filepath
    clusters with less than m_doms domains are not returned
    """
    print("Filtering clusterfile")
    filtered = 0
    with open(infile, 'r') as inf:
        clus_dict = OrderedDict()
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
                clus_dict[line[0]] = line[1:]
            else:
                print("Clusternames not unique, {} read twice".format(clus))
    print("Done. Keeping {} clusters".format(len(clus_dict)))
    print(" {} clusters have less than {} domains".format(filtered,m_doms))
    return clus_dict

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
    if not dom_p1 or not dom_p2:
        return
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
    #with imap specify chunksize as number of pairs-1
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

if __name__ == "__main__":
    cmd = get_commands()
    #read_clusterfile filters out bgcs with less than min_dom domains
    dom_dict = read_clusterfile(cmd.in_file, cmd.min_doms, cmd.verbose)
    clus_names = list(dom_dict.keys())#[0:300] #make list so bgcs have an index
    edges = generate_edges(clus_names, dom_dict, cmd.sim_cutoff, cmd.cores)
    graph = generate_graph(edges)
    print('ARGH01000000_KB894962.1.cluster038:',graph.adj['ARGH01000000_KB894962.1.cluster038'])
    
    #find community structure using some some community/clique algorithm
    cliqs = nx.algorithms.clique.find_cliques(graph)
    # cliqs = nx.algorithms.community.greedy_modularity_communities(graph)
    cliqs = sorted(cliqs, key=len, reverse = True)
    # print(cliqs[0])
    for c in cliqs:
        #incorporate '-'s
        domlist = [dom_dict[clus] for clus in c]
        domlist_del_empty = []
        for doms in domlist:
            doms_del_empty = [d for d in doms if not d == '-']
        domlens = [len(doms) for doms in domlist_del_empty]
        #if there are multiple clus with biggest size, choose one?
        clus_maxlen = [clus for clus, dlen in zip(c,domlens) \
            if dlen == max(domlens)]
        print(max(domlens), clus_maxlen, domlens)
    
    # visualising, maybe with subplot plot all different cliques/networks
    # cols = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    # options = {'node_size': 2,'width': 0.2}
    # pos = nx.spring_layout(graph)
    # plt.figure()
    # nx.draw_networkx(graph, pos=pos, with_labels = False, node_color='black',\
        # **options)
    # for c in cliqs:
        # nx.draw_networkx_nodes(graph, pos=pos, nodelist=c, \
            # node_color=choice(cols), **options)
    # plt.show()

