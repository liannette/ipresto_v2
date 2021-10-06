#!/usr/bin/env python3
"""
Author: Joris Louwen
Script to find modules with LDA algorithm.
"""

import argparse
from collections import Counter, defaultdict
from functools import partial
import matplotlib
matplotlib.use('Agg') #to not rely on X-forwarding (not available in screen)
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
from numpy import sqrt
import numpy as np
from operator import itemgetter
import os
import pandas as pd
# import scipy.cluster.hierarchy as sch
# import seaborn as sns
from statistics import mean,median
import subprocess
from sys import argv
import time


def get_commands():
    parser = argparse.ArgumentParser(description="A script to link modules\
        from the statistical method with known subclusters.")
    parser.add_argument("-m", "--modfile", dest="modfile", help="Input \
        txt file of putative modules. Last column should contain\
        modules", default=False)
    parser.add_argument("-o", "--out_folder", dest="out_folder", help="Output\
        folder", required=True)
    parser.add_argument("-c", "--cores", dest="cores", help="Amount \
        of cores to use, default = all available cores",\
        default=cpu_count(), type=int)
    parser.add_argument("--classes", help="A file containing classes of the \
        BGCs used in the analysis. First column should contain matching BGC\
        names. Consecutive columns should contain classes.", default=False)
    parser.add_argument("--plot", help="If provided: make plots about \
        several aspects of the output. Default is off.", default=False, \
        action="store_true")
    parser.add_argument("--known_subclusters", help="A tab delimited file \
        with known subclusters. Should contain subclusters in the last column\
        and BGC identifiers in the first column. Subclusters are comma \
        separated genes represented as domains. Multiple domains in a gene \
        are separated by semi-colon.")
    parser.add_argument("--bgc_with_mods",help="A tab delimited file linking\
        bgcs to modules. First column should be the name of a bgc and the \
        last column should be the module numbers separated by '; '. These \
        numbers should be the first column in the modfile")
    return parser.parse_args()


def read2dict(filepath, sep=',',header=False):
    '''Read file into a dict {first_column:[other_columns]}

    filepath: str
    sep: str, delimiter in the file
    header: bool, ignore first line
    '''
    output = {}
    with open(filepath,'r') as inf:
        if header:
            inf.readline()
        for line in inf:
            line = line.strip().split(sep)
            output[line[0]] = line[1:]
    return output



if __name__ == '__main__':
    start = time.time()
    #files provided should be filtered bgc csv file and filtered module file

    print('\nStart')
    cmd = get_commands()
    bgcs = read2dict(cmd.bgcfile)
    with open(cmd.modfile, 'r') as inf:
        modules = {}
        #{modules:[info]}
        for line in inf:
            line = line.strip().split('\t')
            mod = tuple(line[-1].split(',')) #now a tuple of str
            modules[mod] = line[:-1]
    if cmd.classes:
        bgc_classes_dict = read2dict(cmd.classes, sep='\t',header=True)
    else:
        bgc_classes_dict = {bgc:'None' for bgc in bgcs}
    if not os.path.isdir(cmd.out_folder):
        subprocess.check_call('mkdir {}'.format(cmd.out_folder), shell=True)

    if cmd.known_subclusters:
        known_subclusters = defaultdict(list)
        with open(cmd.known_subclusters,'r') as inf:
            for line in inf:
                line = line.strip().split('\t')
                known_subclusters[line[0]].append(line[1:])
    else:
        known_subclusters = False

    mod_nums2bgc=defaultdict(list)
    with open(cmd.bgc_with_mods,'r') as inf:
        for line in inf:
            line = line.strip('\n').split('\t')
            mod_nums = line[-1].split('; ')
            if mod_nums:
                for m_num in mod_nums:
                    mod_nums2bgc[m_num].append(line[0])

    outfile = os.path.join(cmd.out_folder, \
        os.path.split(cmd.modfile)[1].split('.txt')[0]+'_reduced.txt')
    

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
