#!/usr/bin/env python3
'''
Find parent modules in clans of statistical method modules.
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
        statistical modules in each clan into parent modules")
    parser.add_argument("-i", "--infile", help="Input file with modules\
        grouped by clan, headerlines are # and each module is on a separate\
        line where the domains are the -2 element on the tsv line", \
        required=True)
    parser.add_argument('-l','--list_file',help='File where modules are\
        listed which will be written again but without redundant modules',
        required=True)
    parser.add_argument('-c', '--cores', help='Cores to use, default = 1',
        default=1, type=int)
    return parser.parse_args()

def read_clans(infile):
    '''
    Read clan file to three dicts: fams per clan, clans with modules, mod_info

    infile: str, filepath
    clan_dict: dict of {clan_number:len_families}
    clan_module: dict of {clan_number:[module_tuples]}
    modules_info: dict of {module_tup:[module_info]}
    '''
    clan_dict = {}
    clan_modules = defaultdict(list)
    modules_info = {}
    with open(infile, 'r') as inf:
        #bval is true if lines are header, false otherwise
        for bval, lines in groupby(inf,key=lambda line: line.startswith('#')):
            if bval:
                #three lines: description, occurences, relative-abundance
                desc = next(lines).strip().split(' ')
                num = int(desc[1].strip(','))
                len_clan = int(desc[2])
                clan_dict[num] = len_clan
            else:
                for line in lines:
                    line = line.strip().split('\t')
                    mod = tuple(line[-2].split(','))
                    modules_info[mod] = line
                    clan_modules[num].append(mod)
    return clan_dict,clan_modules,modules_info

def find_redundant_modules(clan, mods, mod_info):
    '''Returns tuple of clan, list of redundant modules

    clan: int, clan number
    mods: list of modules in the clan (tuples of domains)
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
    return (clan,dels)

def find_redundant_modules_tup(in_tuple):
    '''Returns tuple of clan, list of redundant modules

    clan: int, clan number
    mods: list of modules in the clan (tuples of domains)
    mod_info: dict of {module:[info]}, second elem should be occurrence
    '''
    clan, mods, mod_info = in_tuple
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
    return (clan,dels)

def return_loop_info(clan_modules, mod_info):
    '''Returns list of (clan, modules, specific_module_info)
    '''
    all_tups = []
    for clan, mods in clan_modules.items():
        specific_info = {mod:mod_info[mod] for mod in mods}
        all_tups.append((clan,mods,specific_info))
    return all_tups

if __name__ == '__main__':
    start = time.time()
    print('Start')
    cmd = get_commands()

    out_by_clan = cmd.infile.split('.txt')[0] + '_reduced.txt'
    out_list = cmd.infile.split('_by_clan.txt')[0] + '_reduced.txt'
    clan_dict,clan_modules,mod_info = read_clans(cmd.infile)
    #check intersection between all, is intersection length of smaller one?
    #then the bigger is parents
    dels_all = []
    # for clan, mods in clan_modules.items():
        # print(find_redundant_modules(clan, mods, mod_info))
    info = return_loop_info(clan_modules,mod_info)
    results = []
    # for inf in info:
        # result = find_redundant_modules_tup(inf)
        # results.append(result)
    print('\nFinding redundant modules')
    pool = Pool(cmd.cores,maxtasksperchild=100)
    results = pool.map(find_redundant_modules_tup, info, chunksize=10)

    dels_dict = {clan:set(dels) for clan,dels in results}
    with open(out_list,'w') as outf, open(cmd.list_file,'r') as inf:
        header=inf.readline()
        outf.write(header)
        for line in inf:
            splitline = line.strip().split('\t')
            mod = tuple(splitline[-2].split(','))
            clan = int(splitline[-1])
            try:
                mods = dels_dict[clan]
            except KeyError:
                outf.write(line)
            else:
                if not mod in mods:
                    outf.write(line)

    with open(out_by_clan,'w') as outf:
        for clan, mods in sorted(clan_modules.items()):
            dels = dels_dict[clan]
            mods = [mod for mod in mods if not mod in dels]
            counts = Counter([dom for mod in mods for dom in mod])
            outf.write(\
                '#Subcluster-clan {}, {} families, {} subclusters\n'.format(\
                clan, clan_dict[clan],len(mods)))
            outf.write('#Occurrences: {}\n'.format(', '.join(\
                [dom+':'+str(c) for dom,c in counts.most_common()])))
            outf.write('#Features: {}\n'.format(', '.join(\
                ['{}:{:.2f}'.format(dom,c/len(mods)) for dom,c in \
                counts.most_common()])))
            for mod in mods:
                outf.write('{}\n'.format('\t'.join(mod_info[mod])))

    print('Removed {} redundant modules'.format(sum(\
        [len(vals) for vals in dels_dict.values()])))

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
