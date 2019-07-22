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
import numpy as np
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
    print('\nReading matches from {}'.format(infile))
    matches = defaultdict(list)
    match_per_bgc = []
    total_matches = 0
    zeros = 0
    with open(infile, 'r') as inf:
        #bval is true if lines are header, false otherwise
        for bval, lines in groupby(inf,key=lambda line: line.startswith('>')):
            if bval:
                bgc = next(lines).strip()[1:]
            else:
                i=0
                for line in lines:
                    if not line[0] == 'c' and not line[0] == 'k':
                        #ignore class and known_subcluster lines
                        line = line.strip().split('\t')
                        ####only consider matches of 2 genes or more
                        feats = [tuple(f.split(':')) for f in \
                            line[-1].split(',')]
                        if len(feats) > 1 or round(float(feats[0][1])) > 1:
                            matches[bgc].append(line)
                            i+=1
                if i != 0:
                    match_per_bgc.append(i)
                    total_matches += i
                else:
                    zeros += 1

    print('{} matches, on average {:.2f} matches per BGC, stdev {}'.format(\
        total_matches,np.mean(match_per_bgc),np.std(match_per_bgc)))
    print('  {} BGCs without any match'.format(zeros))
    return matches

def link_matches_to_strain(bgc_matches, strains, outfile):
    '''Write matrix of topic presence/absence for each strain

    bgc_matches: dict of {bgc: [[matches]]}
    strains: dict of {id: accession}
    '''
    print('\nLinking matches to strains')
    strain_set = set()
    topic_dict = defaultdict(set)
    strain_dict = defaultdict(set)
    for bgc, matches in bgc_matches.items():
        strain_id = bgc.split('.')[0].split('_')[0]
        strain_acc = strains[strain_id]
        for match in matches:
            topic = match[0] #topic id
            topic_dict[topic].add(strain_acc)
            strain_dict[strain_acc].add(topic)
    topics_per_strain = [len(vals) for vals in strain_dict.values()]
    print('{} strains with topics linked to them, on average {:.2f}'.format(\
        len(strain_dict.keys()),np.mean(topics_per_strain)) +\
        ' topics per strain, stdev {:.2f}'.format(np.std(topics_per_strain)))
    print('  {} strains with 2 or less topics linked'.format(len([1 for \
        num_topics in topics_per_strain if num_topics <= 2])))

    #same out format as motiftable.tsv from justin. cols=strains, rows=topics
    print('\nWriting motiftable to {}'.format(outfile))
    sorted_cols = sorted(strain_dict)
    sorted_rows = sorted(topic_dict.keys(),key=lambda x: int(x))
    with open(outfile,'w') as outf:
        outf.write('id\t{}\n'.format('\t'.join(sorted_cols)))
        for row in sorted_rows:
            presence = topic_dict[row]
            presence_absence = ''
            for col in sorted_cols:
                pr_ab = '\t1' if col in presence else '\t0'
                presence_absence += pr_ab
            outf.write('{}{}\n'.format(row,presence_absence))

if __name__ == '__main__':
    cmd = get_commands()

    strain2acc = {}
    with open(cmd.strain_ids_file, 'r') as inf:
        for line in inf:
            line = line.strip().split(',')
            if len(line)==2:
                strain2acc[line[0]] = line[1]
    #looked for manually
    strain2acc['AZWM01000001'] = 'CNT796'
    strain2acc['AZWP01000001'] = 'CNS103'

    bgc2matches = read_matches(cmd.in_file)
    link_matches_to_strain(bgc2matches, strain2acc, cmd.out_file)
