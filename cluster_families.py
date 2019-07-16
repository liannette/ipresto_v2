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

def read_families(infile, occ=False):
    '''
    Read family file to three dicts: header_descr, fams with modules, mod_info

    infile: str, filepath
    family_dict: dict of {family_number:len_family}
    feat_dict: dict of {family_number:{feat1:score_feature,
        feat2:score_feature} } weights for each feature (relative abundance of
        a domain-combinations in a family)
    family_module: dict of {family_number:[module_tuples]}
    modules_info: dict of {module_tup:[module_info]}
    '''
    family_dict = {}
    feat_dict = {}
    family_modules = defaultdict(list)
    modules_info = {}
    with open(infile, 'r') as inf:
        #bval is true if lines are header, false otherwise
        for bval, lines in groupby(inf,key=lambda line: line.startswith('#')):
            if bval:
                #three lines: description, occurences, relative-abundance
                desc = next(lines).strip().split(' ')
                num = int(desc[1].strip(','))
                len_fam = desc[2]
                family_dict[num] = int(len_fam)
                if occ:
                    feats = [dom.split(':') for dom in \
                        next(lines).strip().split('#Occurrences: ')[1].split(', ')]
                    tot_occ = sum(map(int,list(zip(*feats))[1]))
                    feat_dict[num] = {f:int(score)/tot_occ for f,score in feats}
                else:
                    next(lines) #discard occurences
                    feats = (dom.split(':') for dom in \
                        next(lines).strip().split('#Features: ')[1].split(', '))
                    feat_dict[num] = {f:float(score) for f,score in feats}
            else:
                for line in lines:
                    line = line.strip().split('\t')
                    mod = tuple(line[-1].split(','))
                    modules_info[mod] = line[:-1]
                    family_modules[num].append(mod)
    return family_dict,feat_dict,family_modules,modules_info

def calc_soergel_dist(fam1,feat_dict1,fam2,feat_dict2):
    '''
    Returns tuple of (fam1,fam2,dist), soergel distance between two families

    fam1,fam2: int, family numbers
    feat_dict1, feat_dict2: {feat1:score_feature, feat2:score_feature, ..}
        weights for each feature (relative abundance of a domain-combinations
        in a family)
    '''
    s1 = set(feat_dict1)
    s2 = set(feat_dict2)
    overl = s1 & s2
    mins = []
    maxs = []
    #calc min and max for each overlapping feature
    for dom in overl:
        both_scores = (feat_dict1[dom],feat_dict2[dom])
        mins.append(min(both_scores))
        maxs.append(max(both_scores))
        print(dom,both_scores)
    # print(s1)
    # print(s2)
    non_over1_score = [feat_dict1[feat] for feat in s1 if feat not in overl]
    non_over2_score = [feat_dict2[feat] for feat in s2 if feat not in overl]
    print(non_over1_score)
    print(non_over2_score)
    print(len(overl),len(non_over1_score),len(non_over2_score))
    #divide sum of mins by sum of (maxs + scores of non_overlapping features)
    numerator = sum(mins)
    denominator = sum(maxs + non_over1_score + non_over2_score)
    print(numerator,denominator)
    soerg = 1 - (numerator / denominator)
    print(soerg)

if __name__ == '__main__':
    cmd = get_commands()

    #make three outfiles: modified list_file with an extra column with clans,
    #fasta_file with >clan to families, file with #clan to modules
    fam_dict,feat_dict,fam_modules,mod_info = read_families(cmd.infile)
    # print(feat_dict)
    f1=3391
    f2=6113
    f3=231
    calc_soergel_dist(f1,feat_dict[f1],f2,feat_dict[f2])
    calc_soergel_dist(f1,feat_dict[f1],f3,feat_dict[f3])
