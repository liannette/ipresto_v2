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
from glob import glob, iglob
from itertools import combinations, product, islice, chain
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

if __name__ == '__main__':
    cmd = get_commands()

    strain2acc = {}
    with open(cmd.strain_ids_file, 'r') as inf:
        for line in inf:
            line = line.strip().split(',')
            if len(line)==2:
                strain2acc[line[0]] = line[1]
        print(strain2acc)

    
