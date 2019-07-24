#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to correlate subclusters to molecular motifs in a group of strains.
FDR is estimated with target-decoy approach.
Matrix tsv files look like:
id\tStrain_1\tStrain2\n
motif1\t0\t1\n
motif2\t1\t0\n

Usage:
python3 correlate_molecular_motifs_to_subclusters.py -h

Example usage:
python3 correlate_molecular_motifs_to_subclusters.py -m molecular_motifs.tsv
    -s subcluster_motifs.tsv -o out_file
'''

import argparse
from collections import OrderedDict, Counter, defaultdict
from itertools import groupby, product
from math import floor, log10
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
import numpy as np
from operator import itemgetter
import os
import random
import re
import seaborn as sns
import subprocess
import time

def get_commands():
    parser = argparse.ArgumentParser(description="A script to correlate\
        subcluster motifs to molecular motifs that are present/absent in\
        a group of organisms. FDR is estimated with target-decoy approach")
    parser.add_argument("-m", "--molecular_motifs", help="Input file which is\
        a presence/absence matrix for molecular motifs in strains in tsv\
        format. The first line contains id followed by strain names\
        (columns). The rows are the motifs, each value is either 0 or 1",\
        required=True)
    parser.add_argument('-s','--subcluster_motifs',help="Input file which is\
        a presence/absence matrix for subcluster motifs in strains in tsv\
        format. The first line contains id followed by strain names\
        (columns). The rows are the motifs, each value is either 0 or 1",\
        required=True)
    parser.add_argument('-o','--out_file',help='Output file',required=True)
    parser.add_argument('-f','--fdr_cutoff',help='FDR cutoff in percentage,\
        default = 1',type=float, default=1)
    return parser.parse_args()

def read_matrix(infile):
    '''
    Returns motif matrix {strain:set(present_motifs)} and set(motifs)

    infile: str, filepath
    '''
    print('\nReading motif file from {}'.format(infile))
    strain2motif = defaultdict(set)
    rownames = set()
    with open(infile,'r') as inf:
        colnames = inf.readline().strip().split('\t')[1:] #ignore 'id'
        for row in inf:
            row = row.strip().split('\t')
            motif = row.pop(0)
            rownames.add(motif)
            for presence,col in zip(row,colnames):
                if presence == '1':
                    strain2motif[col].add(motif)
    filtered_strain2motif = filter_motif_dict(strain2motif)
    print('  filtered out {} strains containing only one motif'.format(\
        len(strain2motif) - len(filtered_strain2motif)))
    print("Motif matrix contains {} strains, {} motifs and {} 1's".format(\
        len(filtered_strain2motif), len(rownames), sum(\
        len(vals) for vals in filtered_strain2motif.values())))
    return filtered_strain2motif, rownames

def filter_motif_dict(motif_dict):
    '''Returns same dict but strains are removed if they contain < 2 motifs

    motif_dict: {strain:[present_motifs]}
    '''
    filtered_motifs = {strain:motifs for strain,motifs in motif_dict.items()\
        if len(motifs) > 1}
    return filtered_motifs

def make_scoring_matrix(m_motifs, s_motifs, m_m_names, s_m_names, strains_used):
    '''Returns list of tuples [(m_motif, s_motif, score)]

    m_motifs, s_motifs: {strain:[present_motifs]} m for molecular-and s
        for subcluster motifs
    m_m_names, s_m_names: list of str, all used motif names for molecular-and
        subcluster motifs respectively

    Score is +10 is both present in a strain, +1 if both absent, 0 if the
    s_motif is there but not m_motif, -10 if m_motif is there but not s_motif
    '''
    #keep track of score as dict of dict {m_motif:{s_motif:score}}
    scoring_matrix = {m_m:defaultdict(int) for m_m in m_m_names}
    for strain in strains_used:
        molec_strain = m_motifs[strain]
        subcl_strain = s_motifs[strain]
        not_molec = m_m_names - molec_strain
        not_subcl = s_m_names - subcl_strain
        for both in product(molec_strain, subcl_strain):
            scoring_matrix[both[0]][both[1]] += 10
        for only_m in product(molec_strain, not_subcl):
            scoring_matrix[only_m[0]][only_m[1]] -= 10
        #only_s not necessary as it is 0
        for neither in product(not_molec,not_subcl):
            scoring_matrix[neither[0]][neither[1]] += 1
    return scoring_matrix

def create_decoy_matrix(motif_matrix, motif_names, strains_used):
    '''Returns a scrambled version of motif_matrix

    motif_matrix: {strain:set(present_motifs)}
    motif_names: list of str, all used motif names
    strains_used: set of str, all available strains to choose from

    Scrambled matrix will contain same amount of strains per motif, but
    just randomly chosen from the strains, so scrambling the presence absence
    in across strains but keeping same presence/absence distribution.
    '''
    print('\nScrambling motif matrix')
    #scramble all the ones across strains and collect in decoy
    decoy = defaultdict(set)
    #first loop through the matrix to get all strains per motif
    motif_dict = defaultdict(set) #set so they can be compared later
    strain_lens = {}
    for strain in strains_used:
        motifs = motif_matrix[strain]
        strain_lens[strain] = len(motifs)
        for motif in motifs:
            motif_dict[motif].add(strain)
    motif_list = list(motif_dict.keys())
    random.shuffle(motif_list)
    for motif in motif_list:
        target_strains = motif_dict[motif]
        length = len(target_strains)
        #choose the same amount of random strains
        scramble = set(random.sample(strains_used,length))
        #make sure target vector does not get in decoy matrix
        while target_strains == scramble:
            print(motif, 'while_loop')
            scramble = set(random.sample(strains_used,length))
        print(motif, length, len(scramble),len(target_strains & scramble))
        for decoy_strain in scramble:
            decoy[decoy_strain].add(motif)
    return decoy

def create_decoy_matrix2(motif_matrix, motif_names, strains_used):
    '''Returns a scrambled version of motif_matrix

    motif_matrix: {strain:set(present_motifs)}
    motif_names: list of str, all used motif names
    strains_used: list of str, strains that are used

    Sort motifs from high to low presence and then for each motif it fills
    the matrix again by sampling the strains that do not have their original
    length yet. In this way it is semi random: choosing randomly the strains
    while keeping same amount of motifs per strain and the same
    present motifs, but the motifs have to be sorted for it to converge.

    Two problems: it only works when having weights for the bigger strains,
    and scrambling rows/colomns with almost all 1-s has little effect on
    end distribution (they are really similar)
    '''
    print('\nScrambling motif matrix')
    #scramble all the ones across strains and collect in decoy
    decoy = defaultdict(set)
    #first loop through the matrix to get all strains per motif
    motif_dict = defaultdict(list)
    strain_lens = {}
    for strain in strains_used:
        motifs = motif_matrix[strain]
        strain_lens[strain] = len(motifs)
        for motif in motifs:
            motif_dict[motif].append(strain)

    used_motifs = list(motif_dict)
    random.shuffle(used_motifs) #shuffle the motifs

    new_strain_lens = defaultdict(int)
    for i,(motif,strns) in enumerate(sorted(motif_dict.items(),\
        key=lambda x: -len(x[1]))):
    # for i, motif in enumerate(used_motifs): #this never works
        length = len(motif_dict[motif])
        possible_strains = []
        weights = []
        for strain,max_len in strain_lens.items():
            current_len = new_strain_lens[strain]
            weight = max_len - current_len
            if weight > 0:
                possible_strains.append(strain)
                weights.append(weight)
        sum_w = sum(weights)
        new_weights = [wght/sum_w for wght in weights]
        # print(i,motif,length,len(possible_strains))
        try:
            # scrambled_strains = random.sample(possible_strains, k=length)
            scrambled_strains = list(np.random.choice(possible_strains,\
                size=length, replace=False, p=new_weights)) #, p=new_weights
        except ValueError:
            create_decoy_matrix2(motif_matrix, motif_names, strains_used)
        # print(len(scrambled_strains))
        j=0
        while scrambled_strains == motif_dict[motif]:
            print(motif,'while_loop')
            if j>=10:
                #probably this version is only option
                create_decoy_matrix2(motif_matrix, motif_names, strains_used)
            try:
                # scrambled_strains = random.sample(possible_strains, k=length)
                scrambled_strains = list(np.random.choice(possible_strains,\
                    size=length, replace=False, p=new_weights)) #, p=new_weights
            except ValueError:
                create_decoy_matrix2(motif_matrix, motif_names, strains_used)
            j+=1
        for scr_strain in scrambled_strains:
            # print(scr_strain,motif)
            new_strain_lens[scr_strain] += 1
            decoy[scr_strain].add(motif)
    return decoy


def plot_scoring_matrix(target, max_len, decoy = False, plot_name = False):
    '''
    '''
    #get all values and convert to numpy array
    scores = np.array(list(zip(*target))[2])
    #could be that some zero values are not in matrix as i didn't initialise
    #does not matter much but for correctness
    if len(scores) != max_len:
        difference = max_len - len(scores)
        scores = np.append(scores, [0]*difference)
    sns.set_style('darkgrid')
    sns.distplot(scores,kde_kws={'color':'#004B96','label':'Target'},\
        hist_kws= {'color':'#004B96','alpha':0.4})
    if decoy:
        s_decoy = np.array(list(zip(*decoy))[2])
        if len(s_decoy) != max_len:
            print(len(s_decoy),len(scores))
            difference = max_len - len(s_decoy)
            s_decoy = np.append(s_decoy, [0]*difference)
        print(len(s_decoy),len(scores))
        sns.distplot(s_decoy,kde_kws={'color':'#DC3220','label':'Decoy'},\
            hist_kws= {'color':'#DC3220','alpha':0.4})
    plt.title('Target and decoy distribution of correlation scores')
    plt.xlabel('Correlation scores')
    plt.ylabel('Density')
    if plot_name:
        plt.savefig(plot_name)
    else:
        plt.show()


if __name__ == '__main__':
    cmd = get_commands()

    molecular_motifs, m_motif_names = read_matrix(cmd.molecular_motifs)
    subcluster_motifs, s_motif_names = read_matrix(cmd.subcluster_motifs)
    #only compare strains present in both matrices
    used_strains = set(molecular_motifs) & set(subcluster_motifs)
    target_matrix = make_scoring_matrix(molecular_motifs, subcluster_motifs,\
        m_motif_names, s_motif_names, used_strains)
    # plot_scoring_matrix(target_matrix)
    molecular_decoy = create_decoy_matrix(molecular_motifs,\
        list(m_motif_names),used_strains)
    subcluster_decoy = create_decoy_matrix(subcluster_motifs,\
        list(s_motif_names),used_strains)
    decoy_matrix = make_scoring_matrix(molecular_decoy, subcluster_decoy,\
        m_motif_names, s_motif_names, used_strains)

    #get all values and convert tuples
    target_tuples = sorted([(mm,sm,score) for mm, val_dict in \
        target_matrix.items() for sm,score in val_dict.items()],\
        key=lambda x: -x[2])
    decoy_tuples = sorted([(mm,sm,score) for mm, val_dict in \
        decoy_matrix.items() for sm,score in val_dict.items()],\
        key=lambda x: -x[2])
    print(target_tuples[:10])
    print(decoy_tuples[:10])

    with open(cmd.out_file,'w') as outf:
        for tup in target_tuples:
            outf.write('{}\n'.format('\t'.join(map(str,tup))))

    # print(decoy_matrix)
    max_length = len(m_motif_names) * len(s_motif_names)
    plot_file = cmd.out_file.split('.')[0] + '.pdf'
    plot_scoring_matrix(target_tuples, max_length, decoy_tuples, plot_file)
    
