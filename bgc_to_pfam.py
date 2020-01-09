#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to convert BGCs into strings of domains.

Usage:
python3 bgc_to_pfam.py -h

Example usage:
python3 bgc_to_pfam.py -i ../testdata -o ../testdata_domains --hmm_path 
    ../domains/Pfam_100subs_tc.hmm --exclude final -c 12 -e True

Notes:
Only handles gbk files with one cluster

Required:
python 3
hmmscan
Biopython
'''
from presto_stat.presto_stat import *

import sys
import os
from glob import glob, iglob
import subprocess
from Bio import SeqIO
from Bio import SearchIO
from collections import OrderedDict, Counter
from multiprocessing import Pool, cpu_count
import argparse
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
    parser.add_argument("--use_fastas", dest="use_fastas", default=None, \
        help="Use already created fasta files from some folder")
    parser.add_argument("--use_domtabs", dest="use_domtabs", default=None, \
        help="Use already created domtables from some folder")
    return parser.parse_args()

if __name__ == "__main__":
    start = time.time()
    cmd = get_commands_local()

    fasta_folder, exist_fastas = process_gbks(cmd.in_folder, \
        cmd.out_folder, cmd.exclude, cmd.exclude_contig_edge, \
        cmd.min_genes, cmd.cores, cmd.verbose, cmd.use_fastas)
    dom_folder, exist_doms = hmmscan_wrapper(fasta_folder, cmd.hmm_path,\
        cmd.verbose, cmd.cores, exist_fastas, cmd.use_domtabs)
    parse_dom_wrapper(dom_folder, cmd.out_folder, cmd.domain_overlap_cutoff,\
        cmd.verbose, exist_doms)

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
