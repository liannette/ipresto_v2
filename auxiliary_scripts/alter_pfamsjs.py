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
from itertools import combinations, product, islice, chain
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
    return parser.parse_args()

if __name__ == "__main__":
    
