#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to convert the genes in BGCs into strings of domains, filter BGCs
based on similarity, and detect sub-clusters based on a statistical method
and an LDA alorithm.

Usage:
python3 get_antismash_classes.py -h

Example usage:
python3 get_antismash_classes.py
'''

import argparse
from Bio import SeqIO
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
import re
import time

def get_commands():
    parser = argparse.ArgumentParser(description="A script to retrieve\
        classes from antiSMASH 4 output gbk files.")
    parser.add_argument("-i", "--in_folder", dest="in_folder", help="Input \
        directory of gbk files", required=True)
    parser.add_argument("-o", "--out_file", dest="out_file", 
        required=True, help="Output file, this will contain all BGC names \
        linked to a class as tsv")
    parser.add_argument("-c", "--cores", dest="cores", default=cpu_count(), 
        help="Set the number of cores the script may use (default: use all \
        available cores)", type=int)
    return parser.parse_args()

def process_gbks(input_folder, output_file, exclude_contig_edge,\
    cores, verbose):
    '''Convert gbk files from input folder to fasta files for each gbk file

    input_folder, output_file: str
    exclude_contig_edge: bool
    verbose: bool, print additional info to stdout
    '''
    print("\nRetreiving classes from gbk files.")
    files = iglob(os.path.join(input_folder, "*.gbk"))
    done = []
    pool = Pool(cores, maxtasksperchild=250)
    for file_path in files:
        pool.apply_async(convert_gbk2fasta, args=(file_path, True, verbose),\
            callback=lambda x: done.append(x))
    pool.close()
    pool.join()
    class_tups = [val for val in done if val]
    filtered = len([val for val in done if val == None])
    print("Processed {} gbk files into {} classes.".format(\
        len(class_tups)+filtered, len(class_tups)))
    print(" filtered {} files that are on contig edge".format(\
        filtered))
    with open(output_file,'w') as outf:
        for class_tup in class_tups:
            outf.write('{}\n'.format('\t'.join(class_tup)))

def convert_gbk2fasta(file_path, exclude_contig_edge, verbose):
    '''Retrieve a tuple of (name,class) from one gbk file

    file_path, out_folder: strings
    exclude_contig_edge: bool
    verbose: bool, print additional info to stdout

    Returns None if there is a contig edge.
    '''
    file_name = os.path.split(file_path)[1]
    name = file_name.split('.gbk')[0]
    try:
        record = next(SeqIO.parse(file_path, 'genbank'))
    except ValueError as e:
        print(" Excluding {}: {}".format(file_path, e))
        return
    for feature in record.features:
        if feature.type == 'cluster':
            if "contig_edge" in feature.qualifiers:
                    if feature.qualifiers["contig_edge"][0] == "True":
                        if exclude_contig_edge:
                            if verbose:
                                print("  excluding {}: {}".format(file_name,\
                                    "contig edge"))
                            return
            if "product" in feature.qualifiers:
                product = feature.qualifiers.get('product',['None'])[0]
                return name,product
    return

if __name__ == '__main__':
    start = time.time()
    cmd = get_commands()

    process_gbks(cmd.in_folder, cmd.out_file, True, cmd.cores, False)

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
