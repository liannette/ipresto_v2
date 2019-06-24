#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to convert fastas and dom_hits file into json files that can be
used for visualisations of bgcs. A nice extension would be to include
which genes belong to putative modules and somehow visualise that too.

Usage:
python3 alter_pfamsjs.py -h

Example usage:
python3 alter_pfamsjs.py -i pfams.js -o outfile.js --hmm hmm_file_with_subPfam

'''

import argparse
from collections import OrderedDict, Counter, defaultdict
# from copy import deepcopy
# from functools import partial
from glob import glob, iglob
from itertools import combinations, product, islice, chain,groupby
import json
# from math import floor, log10
from multiprocessing import Pool, cpu_count
# from operator import itemgetter
import os
# import random
# from statsmodels.stats.multitest import multipletests
import subprocess
import time

def get_commands():
    parser = argparse.ArgumentParser(description="A script to alter pfams.js\
        so that it contains Pfam IDs instead of accessions including subPfam")
    parser.add_argument("-i", "--in_file", help="Input file pfam.js",
        required=True)
    parser.add_argument("--hmm", help="Uncompressed Pfam database file .hmm")
    parser.add_argument("-o", "--out_file", dest="out_file", 
        required=True, help="Output file")
    parser.add_argument("--sub_count",help='Counts of subPfam clades so that\
        the right amount of subPfams can be added')
    return parser.parse_args()

def parse_hmm(hmm_file):
    '''
    '''
    acc2name = {}
    with open(hmm_file, 'r') as inf:
    #returns a bool in boolval and an iterable in record
        for boolval, record in groupby(inf, lambda line:\
            line.startswith('HMMER3/f')):
            if boolval:
                header = next(record) #will only contain 1 line
            else:
                name = next(record).rstrip('\n').split()[1]
                acc = next(record).rstrip('\n').split()[1].split('.')[0]
                acc2name[acc]=name
    return acc2name

def process_pfamjs(pfam_file, out_file, acc2name, sub_count):
    '''
    '''
    #copy whole file to out_file but wihout the var pfams= so it is a json

if __name__ == "__main__":
    cmd = get_commands()

    # acc_to_name = parse_hmm(cmd.hmm)
    # print(acc_to_name,len(acc_to_name))
    read_pfamjs(cmd.in_file,cmd.out_file,{},{})
