#!/usr/bin/env python3
'''
Find parent modules in statistical method modules.
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
        statistical modules in each family into parent modules")
    parser.add_argument("-i", "--infile", help="Input file with modules\
        grouped by family, headerlines are # and each module is on a separate\
        line where the domains are the last element on the tsv line", \
        required=True)
    parser.add_argument('-l','--list_file',help='File where modules are\
        listed which will be written again but without redundant modules',
        required=True)
    parser.add_argument('-c', '--cores', help='Cores to use, default = 1',
        default=1, type=int)
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

def find_redundant_modules(fam, mods, mod_info):
    '''Returns tuple of fam, list of redundant modules

    fam: int, family number
    mods: list of modules in the family (tuples of domains)
    mod_info: dict of {module:[info]}, second elem should be occurrence
    '''
    dels = []
    parent_dict = defaultdict(list) #record parents for each module
    pairs = combinations(mods,2)
    for pair in pairs:
        p1,p2 = pair
        inter = set(p1).intersection(set(p2))
        if len(inter) == len(p1):
            #p2 is parent of p1
            parent_dict[p1].append(p2)
        if len(inter) == len(p2):
            parent_dict[p2].append(p1)
    for child,parents in sorted(parent_dict.items(),key=len):
        # print(child,parents)
        occ_c = int(mod_info[child][1])
        occ_p = max(map(int,[mod_info[p][1] for p in parents]))
        # print(occ_c,occ_p)
        if occ_c <= occ_p:
            dels.append(child)
    return (fam,dels)

def find_redundant_modules_yield(in_tuple):
    '''Returns tuple of fam, list of redundant modules

    fam: int, family number
    mods: list of modules in the family (tuples of domains)
    mod_info: dict of {module:[info]}, second elem should be occurrence
    '''
    fam, mods, mod_info = in_tuple
    dels = []
    parent_dict = defaultdict(list) #record parents for each module
    pairs = combinations(mods,2)
    for pair in pairs:
        p1,p2 = pair
        inter = set(p1).intersection(set(p2))
        if len(inter) == len(p1):
            #p2 is parent of p1
            parent_dict[p1].append(p2)
        if len(inter) == len(p2):
            parent_dict[p2].append(p1)
    for child,parents in sorted(parent_dict.items(),key=len):
        # print(child,parents)
        occ_c = int(mod_info[child][1])
        occ_p = max(map(int,[mod_info[p][1] for p in parents]))
        # print(occ_c,occ_p)
        if occ_c <= occ_p:
            dels.append(child)
    return (fam,dels)

def yield_loop_info(fam_modules, mod_info):
    '''Returns generator of (family, modules, specific_module_info)
    '''
    for fam, mods in fam_modules.items():
        specific_info = {mod:mod_info[mod] for mod in mods}
        yield (fam,mods,specific_info)

if __name__ == '__main__':
    start = time.time()
    print('Start')
    cmd = get_commands()

    out_by_fam = cmd.infile.split('.txt')[0] + '_reduced.txt'
    out_list = cmd.infile.split('_by_family.txt')[0] + '_reduced.txt'
    fam_dict,fam_modules,mod_info = read_families(cmd.infile)
    # print(fam_modules)
    #check intersection between all, is intersection length of smaller one?
    #then the bigger is parents
    dels_all = []
    # for fam, mods in fam_modules.items():
        # print(find_redundant_modules(fam, mods, mod_info))
    info = yield_loop_info(fam_modules,mod_info)
    results = []
    # for inf in info:
        # result = find_redundant_modules_yield(inf)
        # results.append(result)
    pool = Pool(cmd.cores)
    results = pool.map(find_redundant_modules_yield, info)
    dels_dict = {fam:set(dels) for fam,dels in results}
    with open(out_list,'w') as outf, open(cmd.list_file,'r') as inf:
        header=inf.readline()
        outf.write(header)
        for line in inf:
            splitline = line.strip().split('\t')
            mod = tuple(splitline[-2].split(','))
            family = int(splitline[-1])
            try:
                mods = dels_dict[family]
            except KeyError:
                outf.write(line)
            else:
                if not mod in mods:
                    outf.write(line)

    with open(out_by_fam,'w') as outf:
        for fam, mods in sorted(fam_modules.items()):
            dels = dels_dict[fam]
            mods = [mod for mod in mods if not mod in dels]
            counts = Counter([dom for mod in mods for dom in mod])
            outf.write('#Subcluster-family {}, {} subclusters\n'.format(fam,\
                len(mods)))
            outf.write('#Occurrences: {}\n'.format(', '.join(\
                [dom+':'+str(c) for dom,c in counts.most_common()])))
            outf.write('#Features: {}\n'.format(', '.join(\
                ['{}:{:.2f}'.format(dom,c/len(mods)) for dom,c in \
                counts.most_common()])))
            for mod in mods:
                outf.write('{}\t{}\n'.format('\t'.join(mod_info[mod]),\
                    ','.join(sorted(mod))))

    print('Removed {} redundant modules'.format(sum(\
        [len(vals) for vals in dels_dict.values()])))


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
