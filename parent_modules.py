#!/usr/bin/env python3
'''
Find parent modules in statistical method modules.
Author: Joris Louwen
'''
import os
import argparse
from collections import defaultdict, Counter
from itertools import groupby, combinations
import numpy as np
import time


def get_commands():
    parser = argparse.ArgumentParser(description="A script to group the\
        statistical modules in each family into parent modules")
    parser.add_argument("-i", "--infile", help="Input file with modules\
        grouped by family, headerlines are # and each module is on a separate\
        line where the domains are the last element on the tsv line", \
        required=True)
    parser.add_argument('-c', '--cores', help='Cores to use, default = 1',
        default=1, type=int)
    return parser.parse_args()

def read_families(infile):
    '''
    Read family file to three dicts: header_descr, fams with modules, mod_info

    infile: str, filepath
    '''
    family_dict = {}
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

def find_redundant_modules(fam, mod

if __name__ == '__main__':
    start = time.time()
    print('Start')
    cmd = get_commands()

    fam_dict,fam_modules,mod_info = read_families(cmd.infile)
    # print(fam_modules)
    #check intersection between all, is intersection length of smaller one?
    #then the bigger is parents
    dels_all = []
    for fam, mod in fam_modules.items():
        dels = []
        parent_dict = defaultdict(list) #record parents for each module
        pairs = combinations(mod,2)
        for pair in pairs:
            p1,p2 = pair
            # print('pairs',p1,p2)
            inter = set(p1).intersection(set(p2))
            if len(inter) == len(p1):
                #p2 is parent of p1
                parent_dict[p1].append(p2)
            if len(inter) == len(p2):
                parent_dict[p2].append(p1)
        for child,parents in sorted(parent_dict.items(),key=len):
            print(child,parents)
            occ_c = int(mod_info[child][1])
            occ_p = max(map(int,[mod_info[p][1] for p in parents]))
            print(occ_c,occ_p)
            if occ_c <= occ_p:
                dels.append(child)
        print(dels,len(dels))
            # print('overlap',inter1,inter2)
            # if inter1 == len(p1):
                # print(1)
    # for k,v in parent_dict.items():
        # if len(v) != 1:
            # print('\nchild')
            # print(k)
            # print('parents')
            # for par in sorted(v,key=len):
                # print(par)
    # print(np.mean([len(val) for val in parent_dict.values()]))

    # print(fam_dict)
    #calc clans?
    # fam_feats = {}
    # for fam,val in fam_dict.items():
        # feat = val[2]
        # feats = [tuple(dom.split(':')) for dom in \
            # feat.split('#Features: ')[1].split(', ')]
        # fam_feats[fam] = feats
    # fam_pairs = combinations(fam_feats.keys())
    

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
