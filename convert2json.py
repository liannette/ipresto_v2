#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to convert fastas and dom_hits file into json files that can be
used for visualisations of bgcs. A nice extension would be to include
which genes belong to putative modules and somehow visualise that too.

Usage:
python3 convert2json.py -h

Example usage:
python3 convert2json.py -d dom_hits.txt -f fasta_folder/ -o output_folder

'''

import argparse
from collections import OrderedDict, Counter, defaultdict
from copy import deepcopy
from functools import partial
from glob import glob, iglob
from itertools import combinations, product, islice, chain
from math import floor, log10
from multiprocessing import Pool, cpu_count
from operator import itemgetter
import os
import random
from statsmodels.stats.multitest import multipletests
import subprocess
import time

def get_commands():
    parser = argparse.ArgumentParser(description="A script to turn bgcs in \
        fasta files and dom_hits file into json files used for visualisation")
    parser.add_argument("-d", "--dom_hits", help="Input file of domain hits",
        required=True)
    parser.add_argument("-f","--fasta_folder", help="Folder of fasta files")
    parser.add_argument("-o", "--out_folder", dest="out_folder", 
        required=True, help="Output directory, this will contain all output \
        data files.")
    return parser.parse_args()

if __name__ == "__main__":
    #parse fasta files to dict format for visualisation with empty domains
    #add domains from parsed dom_hits file
    #for ending point of a gbk just stop at last gene? or maybe add 100
    #neater to directly parse gbks maybe
    cmd = get_commands()
