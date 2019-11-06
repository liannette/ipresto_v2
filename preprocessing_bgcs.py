#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to convert BGCs into strings of domains and filter the domains based
on similarity.

Usage:
python3 preprocessing_bgcs.py -h

Example usage:
python3 bgc_to_pfam.py -i ../testdata -o ../testdata_domains --hmm_path 
    ../domains/Pfam_100subs_tc.hmm --exclude final -c 12 -e True

Notes:
Only handles gbk files with one cluster

Required:
python 3.6
hmmscan
Biopython
networkx
'''
import sys
sys.path.insert(1, sys.path[0]+'/PRESTO-STAT/')
from presto_stat import *

import argparse
from Bio import SeqIO
from Bio import SearchIO
from collections import OrderedDict, Counter
from functools import partial
from glob import glob, iglob
from itertools import combinations
from multiprocessing import Pool, cpu_count
import networkx as nx
import os
import random
import subprocess
import time


def get_commands_local():
    parser = argparse.ArgumentParser(description="A script to turn bgcs from \
        gbk files into strings of domains using a domain hmm database and to \
        reduce redundancy by filtering out similar bgcs.")
    parser.add_argument("-i", "--in_folder", dest="in_folder", help="Input \
        directory of gbk files", required=True)
    parser.add_argument("--exclude", dest="exclude", default=["final"],
        nargs="+", help="If any string in this list occurs in the gbk \
        filename, this file will not be used for the analysis. \
        (default: final)")
    parser.add_argument("-o", "--out_folder", dest="out_folder", 
        required=True, help="Output directory, this will contain all output \
        data files.")
    parser.add_argument("--hmm_path", dest="hmm_path", required=True,
        help="File containing domain hmms that is hmmpress-processed.")
    parser.add_argument("-c", "--cores", dest="cores", default=cpu_count(), 
        help="Set the number of cores the script may use (default: use all \
        available cores)", type=int)
    parser.add_argument("-v", "--verbose", dest="verbose", required=False,
        action="store_true", default=False, help="Prints more detailed \
        information.")
    parser.add_argument("-d", "--domain_overlap_cutoff", 
        dest="domain_overlap_cutoff", default=0.1, help="Specify at which \
        overlap percentage domains are considered to overlap. Domain with \
        the best score is kept (default=0.1).")
    parser.add_argument("-e", "--exclude_contig_edge",
        dest="exclude_contig_edge", default=True, type=bool, help="\
        Exclude clusters that lie on a contig edge")
    parser.add_argument("-m", "--min_genes", dest="min_genes", default=0,
        help="Provide the minimum size of a BGC to be included in the \
        analysis. Default is 0 genes", type=int)
    parser.add_argument("--min_doms", dest="min_doms", default=0,
        help="The minimum amount of domains in a BGC to be included in the \
        analysis. Default is 0 domains", type=int)
    parser.add_argument("--sim_cutoff", dest="sim_cutoff", default=0.95,
        help="Cutoff for cluster similarity in redundancy filtering (default:\
        0.95)", type=float)
    parser.add_argument("--use_fastas", dest="use_fastas", default=None, \
        help="Use already created fasta files from some folder")
    parser.add_argument("--use_domtabs", dest="use_domtabs", default=None, \
        help="Use already created domtables from some folder")
    parser.add_argument("--include_list", dest="include_list", default=None, \
        help="If provided only the domains in this file will be taken into \
        account in the analysis. One line should contain one Pfam ID \
        (default: None - meaning all Pfams from database)")
    parser.add_argument("--start_from_clusterfile", default=None, help="A file\
        with BGCs and domain-combinations to start with (csv and domains in a\
        gene separated by ';'). This overwrites in_folder (which still has to\
        be supplied symbolically) and use_domtabs/use_fastas.")
    parser.add_argument("--no_redundancy_filtering",default=False,help="If \
        provided, redundancy filtering will not be performed",\
        action="store_true")
    return parser.parse_args()

if __name__ == "__main__":
    start = time.time()
    cmd = get_commands_local()

    #converting genes in each bgc to a combination of domains
    if cmd.start_from_clusterfile:
        if not os.path.isdir(cmd.out_folder):
            f_command = 'mkdir {}'.format(cmd.out_folder)
            subprocess.check_call(f_command,shell=True)
        filepre = os.path.split(cmd.start_from_clusterfile)[-1].split(\
            '.csv')[0]
        clus_file = os.path.join(cmd.out_folder, filepre+'_clusterfile.csv')
        c_command = 'cp {} {}'.format(cmd.start_from_clusterfile,clus_file)
        subprocess.check_call(c_command, shell=True)
    else:
        fasta_folder, exist_fastas = process_gbks(cmd.in_folder, \
            cmd.out_folder, cmd.exclude, cmd.exclude_contig_edge, \
            cmd.min_genes, cmd.cores, cmd.verbose, cmd.use_fastas)
        dom_folder, exist_doms = hmmscan_wrapper(fasta_folder, cmd.hmm_path,\
            cmd.verbose, cmd.cores, exist_fastas, cmd.use_domtabs)
        clus_file = parse_dom_wrapper(dom_folder, cmd.out_folder, \
            cmd.domain_overlap_cutoff, cmd.verbose, exist_doms)

    #filtering clusters based on similarity
    random.seed(595)
    dom_dict = read_clusterfile(clus_file, cmd.min_genes, \
        cmd.verbose)
    doml_dict = {bgc: sum(len(g) for g in genes if not g == ('-',)) \
        for bgc,genes in dom_dict.items()}
    filt_file = '{}_filtered_clusterfile.csv'.format(\
        clus_file.split('_clusterfile.csv')[0])
    if not os.path.isfile(filt_file):
        #do not perform redundancy filtering if it already exist
        if not cmd.no_redundancy_filtering:
            edges_file = generate_edges(dom_dict, cmd.sim_cutoff,\
                cmd.cores, cmd.out_folder)
            similar_bgcs = read_edges_from_temp(edges_file)
            graph = generate_graph(similar_bgcs, True)
            uniq_bgcs = [clus for clus in dom_dict.keys() if not clus in \
                graph.nodes()]
            all_reps = find_all_representatives(doml_dict, graph)
        else:
            #dont perform redundancy filtering and duplicate clus_file to
            #filt file, representative file is created but this is symbolic
            uniq_bgcs = list(dom_dict.keys())
            all_reps = {}
            print('\nRedundancy filtering is turned off.')
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
