#!/usr/bin/env python3
"""
Author: Joris Louwen
Script to link modules back to BGCs.
"""

from sys import argv
from multiprocessing import Pool
from functools import partial
from collections import Counter
import time

from copy import deepcopy

def read_mods_bgcs(modsfile,bgcsfile):
    '''Returns two dicts, one with {mod:[info]} and one with {bgc:[domains]}

    modsfile, bgcsfile: strings, filenames
    '''
    with open(mods_file, 'r') as mods_in, open(bgcs_file, 'r')  as bgcs_in:
        mods = {}
        mods_in.readline()
        for line in mods_in:
            line = line.strip().split('\t')
            mod = tuple(line[-1].split(','))
            mods[mod] = line[:-1]
        bgcs = {}
        for line in bgcs_in:
            line = line.strip().split(',')
            bgcs[line[0]] = line[1:]
    return mods, bgcs

def link_mods2bgc(bgc, doms, modules):
    '''Returns a tuple of (bgc, [(modules)])

    bgc: string, bgc name
    doms: list of strings, all domain names in bgc
    modules: list of tuples of strings, each tuple contains the domains of a
        module
    '''
    modlist = []
    for mod in modules:
        if set(doms).intersection(mod) == set(mod):
            modlist.append(mod)
    return (bgc,sorted(modlist,key=len))

def link_all_mods2bgcs(bgcs, modules, cores):
    '''Returns a dict of {bgc: [(modules)]}

    bgcs: dict of {bgc: [domains])
    modules: list of module tuples
    cores: int, amount of cores to use
    '''
    pool = Pool(cores, maxtasksperchild=100)
    bgcs_mod = pool.starmap(partial(link_mods2bgc, modules=modules), \
        bgcs.items())
    bgc_mod_dict = {pair[0]:pair[1] for pair in bgcs_mod}
    return bgc_mod_dict

def remove_infr_mods(bgc_mod_dict, modules_dict):
    '''Returns updated input where modules that occur < 2 times are removed

    bgc_mod_dict: dict of {bgc: [(modules)]}
    modules: dict  of {mod:[info]}
    '''
    new_bgc_dict = deepcopy(bgc_mod_dict)
    new_mods_dict = deepcopy(modules_dict)
    mod_counts = Counter(modules_dict.keys())
    mod_counts.update([mod for modlist in bgc_mod_dict.values() \
        for mod in modlist])
    #I initialised all mods with 1 so the if statements says < 3 instead of 2
    infr_mods = [mod for mod, count in mod_counts.items() if count < 3]
    print('\nRemoving {:.1f}% of modules that occur less than twice'.format(\
        len(infr_mods)/len(modules_dict)*100))
    for infr_mod in infr_mods:
        del new_mods_dict[infr_mod]
    for bgc,mods in new_bgc_dict.items():
        new_bgc_dict[bgc] = [mod for mod in mods if not mod in infr_mods]
    return new_bgc_dict,new_mods_dict

def write_bgcs_and_modules(outfile, clusters, bgc_mod_dict):
    '''Writes a files containing Bgcname\tDomains\tModules

    clusters: dict linking bgc name to its domains {bgc:[domains]}
    bgc_mod_dict: dict linking bgc to its modules {bgc: [(modules)]}
    Before each domain its index is written seperated with a .
    If there are multiple indeces they are separated with a -
    Modules are seperated by ';'
    outfile:
    clustername    dom1,dom2,dom3,dom2    0.dom1,1-3.dom2,2.dom3;0.dom1,2.dom3
    '''
    print('\nWriting file that links bgcs to modules at {}'.format(outfile))
    with open(outfile, 'w') as outf:
        for bgc, doms in clusters.items():
            dom_i = {}
            for i,dm in enumerate(doms):
                try:
                    dom_i[dm].append(str(i))
                except KeyError:
                    dom_i[dm] = [str(i)]
            dom_str = {dm:'-'.join(inds)+'.'+dm for dm,inds in dom_i.items()}
            mods_str = '; '.join([','.join([dom_str[d] for d in mod]) \
                for mod in bgc_mod_dict[bgc]])
            #I could add also as a last column the ids of the modules
            outf.write('{}\t{}\t{}\n'.format(bgc,','.join(doms),mods_str))

if __name__ == '__main__':
    mods_file = argv[1]
    bgcs_file = argv[2]
    cors = int(argv[3])
    print('\nStart')
    mods, bgcs = read_mods_bgcs(mods_file, bgcs_file)
    bgcs_with_mods = link_all_mods2bgcs(bgcs, mods, cors)
    bgcs_with_mods, mods = remove_infr_mods(bgcs_with_mods, mods)
    bgcmodfile = '{}_bgcs_with_domains.txt'.format(\
        mods_file.split('_modules.txt')[0])
    write_bgcs_and_modules(bgcmodfile, bgcs, bgcs_with_mods)
    

    #after finding which occur, I can also figure out where the domains are
