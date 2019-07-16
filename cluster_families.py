#!/usr/bin/env python3
'''
Find clans of families of statistical method modules.
Author: Joris Louwen
'''
import os
import argparse
from collections import defaultdict, Counter
from itertools import groupby, combinations
from multiprocessing import Pool
import numpy as np
import time


def get_commands():
    parser = argparse.ArgumentParser(description="A script to group the\
        families of the statistical modules in clans (related families) using\
        a soergel distance metric (weighted jaccard distance)")
    parser.add_argument("-i", "--infile", help="Input file with modules\
        grouped by family, headerlines are # and each module is on a separate\
        line where the domains are the last element on the tsv line", \
        required=True)
    parser.add_argument('-l','--list_file',help='File where modules are\
        listed which will be written again but with clans added',
        required=True)
    parser.add_argument('-c', '--cores', help='Cores to use, default = 1',
        default=1, type=int)
    parser.add_argument('--cutoff', help='Cutoff for when two families are\
        similar, default = 0.5', default=0.5)
    return parser.parse_args()

def read_families(infile):
    '''
    Read family file to three dicts: header_descr, fams with modules, mod_info

    infile: str, filepath
    family_dict: dict of {family_number:[header_info]}
    family_module: dict of {family_number:[module_tuples]}
    modules_info: dict of {module_tup:[module_info]}
    '''
    family_dict = {}
    feat_dict = defaultdict(dict)
    family_modules = defaultdict(list)
    modules_info = {}
    with open(infile, 'r') as inf:
        #bval is true if lines are header, false otherwise
        for bval, lines in groupby(inf,key=lambda line: line.startswith('#')):
            if bval:
                #three lines: description, occurences, relative-abundance
                desc = next(lines).strip().split(' ')
                num = int(desc[1].strip(','))
                len_fam = desc[2:3]
                family_dict[num] = len_fam + [l.strip() for l in lines]
            else:
                for line in lines:
                    line = line.strip().split('\t')
                    mod = tuple(line[-1].split(','))
                    modules_info[mod] = line[:-1]
                    family_modules[num].append(mod)
    return family_dict,family_modules,modules_info

if __name__ == '__main__':
    cmd = get_commands()

    #make three outfiles: modified list_file with an extra column with clans,
    #fasta_file with >clan to families, file with #clan to modules
