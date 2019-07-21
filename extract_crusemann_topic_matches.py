#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to get a matrix of the absence/presence of each topic(subcluster-motif)
in each strain of the Crusemann dataset.

Usage:
python3 extract_crusemann_topic_matches.py -h

Example usage:
python3 extract_crusemann_topic_matches.py -i bgc_topic_filtered.txt
'''

import argparse
from collections import OrderedDict, Counter, defaultdict
from itertools import groupby
from math import floor, log10
from multiprocessing import Pool, cpu_count
from operator import itemgetter
import os
import random
import re
import subprocess
import time

def get_commands():
    parser = argparse.ArgumentParser(description="A script to extract topic\
        matches for the crusemann bgcs and group them per strain by topic.")
    parser.add_argument("-i", "--in_file", help="Input file containing topic\
        matches per bgc in fasta like format", required=True)
    parser.add_argument('-s','--strain_ids_file',help='Input file linking\
        strain ids to accessions, "id, accession" on every line',\
        required=True)
    parser.add_argument('-o','--out_file',help='Output file',required=True)
    return parser.parse_args()

def read_matches(infile):
    '''Read bgc_topics_filtered file to dict {bgc:[[matches]]}

    infile: str, filepath of bgc_topics_filtered file
    '''
    matches = defaultdict(list)
    with open(infile, 'r') as inf:
        #bval is true if lines are header, false otherwise
        for bval, lines in groupby(inf,key=lambda line: line.startswith('>')):
            if bval:
                bgc = next(lines).strip()[1:]
            else:
                for line in lines:
                    if not line[0] == 'c' and not line[0] == 'k':
                        #ignore class and known_subcluster lines
                        line = line.strip().split('\t')
                        ####apply some filtering here
                        matches[bgc].append(line)
    return matches

def link_matches_to_strain(bgc_matches, strains):
    '''Write matrix of topic presence/absence for each strain

    bgc_matches: dict of {bgc: [[matches]]}
    strains: dict of {id: accession}
    '''
    strain_set = set()
    for bgc, matches in bgc_matches.items():
        strain_id = bgc.split('.')[0].split('_')[0]
        
        strain_set.add(strain_id)
    print(sorted(strain_set),len(strain_set))
    print(sorted(strains),len(strains))
    non_over = [st for st in strain_set if st not in strains]
    print(sorted(non_over),len(non_over))
    

if __name__ == '__main__':
    cmd = get_commands()

    strain2acc = {}
    with open(cmd.strain_ids_file, 'r') as inf:
        for line in inf:
            line = line.strip().split(',')
            if len(line)==2:
                strain2acc[line[0]] = line[1]

    bgc2matches = read_matches(cmd.in_file)
    link_matches_to_strain(bgc2matches, strain2acc)
