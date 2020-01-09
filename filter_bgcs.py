#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to convert tokenised bgcs into graphs based on a Adjacency index
and filter out redundant bgcs.

Usage:
python3 bgc_to_pfam.py -h

Example usage:
python3 filter_bgcs.py -i ../testdata_domains/testdata_clusterfile.csv
    -o ../testdata_domains/ -c 20

Required:
networkx (https://github.com/networkx/networkx)
'''
from presto_stat.presto_stat import *

import sys
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
import time

def get_commands_local():
    parser = argparse.ArgumentParser(description="Script to turn a \
        clusterfile with tokenised BGCs into similarity networks and \
        filter out redundant bgcs.")
    parser.add_argument("-i", "--in_file", dest="in_file", help="Input \
        clusterfile: GbkName,dom1_gene1,dom1_gene2;dom2:gene2", required=True)
    parser.add_argument("-o", "--out_folder", dest="out_folder", 
        required=True, help="Output directory, this will contain all output \
        data files.")
    parser.add_argument("-c", "--cores", dest="cores", default=cpu_count(), 
        help="Set the number of cores the script may use (default: use all \
        available cores)", type=int)
    parser.add_argument("-v", "--verbose", dest="verbose", required=False,
        action="store_true", default=False, help="Prints more detailed \
        information.")
    parser.add_argument("-m", "--min_genes", dest="min_genes", default=0,
        help="Provide the minimum size of a BGC to be included in the \
        analysis. Default is 0 genes", type=int)
    parser.add_argument("--sim_cutoff", dest="sim_cutoff", default=0.95,
        help="Cutoff for cluster similarity in redundancy filtering (default:\
        0.95)", type=float)
    parser.add_argument("--include_list", dest="include_list", default=None, \
        help="If provided only the domains in this file will be taken into \
        account in the analysis. One line should contain one Pfam ID \
        (default: None - meaning all Pfams from database)")
    return parser.parse_args()

if __name__ == "__main__":
    start = time.time()
    cmd = get_commands_local()

    #filtering clusters based on similarity
    random.seed(595)
    clus_file = cmd.in_file
    dom_dict = read_clusterfile(clus_file, cmd.min_genes, \
        cmd.verbose)
    doml_dict = {bgc: sum(len(g) for g in genes if not g == ('-',)) \
        for bgc,genes in dom_dict.items()}
    filt_file = '{}_filtered_clusterfile.csv'.format(\
        os.path.join(cmd.out_folder,\
        os.path.split(clus_file)[1].split('_clusterfile.csv')[0]))
    if not os.path.isfile(filt_file):
        #do not perform redundancy filtering if it already exist
        edges_file = generate_edges(dom_dict, cmd.sim_cutoff,\
            cmd.cores, cmd.out_folder)
        similar_bgcs = read_edges_from_temp(edges_file)
        graph = generate_graph(similar_bgcs, True)
        uniq_bgcs = [clus for clus in dom_dict.keys() if not clus in \
            graph.nodes()]
        all_reps = find_all_representatives(doml_dict, graph)
        if cmd.include_list:
            include_list = read_txt(cmd.include_list)
            dom_dict = filter_out_domains(dom_dict, include_list)
        all_reps_file = write_filtered_bgcs(uniq_bgcs, all_reps, \
            dom_dict, filt_file)
    else:
        print('\nFiltered clusterfile existed, redundancy filtering not'+
            ' performed again')

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
