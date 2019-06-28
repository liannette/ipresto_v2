#!/usr/bin/env python3
"""
Author: Joris Louwen
Script to find the differences between two domstats files
"""

from sys import argv
import re
from collections import defaultdict

def read_f(filepath):
    '''reads file in a dict {line[0]:line[1]}
    '''
    with open(filepath,'r') as inf:
        header = inf.readline()
        dct = {}
        for line in inf:
            line = line.strip('\n').split('\t')
            dct[line[0]] = int(line[1])
    return dct

if __name__ == '__main__':
    f_ori = argv[1]
    f_comp = argv[2]
    f_ecdm = argv[3]
    #find doms in original file not in compare file
    #output also a file with dif domains not in ecdm list

    ori = read_f(f_ori)
    comp = read_f(f_comp)
    with open(f_ecdm,'r') as inf:
        ecdm = [line.strip() for line in inf]
    dif = {dom:ori[dom] for dom in ori if not dom in comp}
    # print(sorted(dif.items(), key=lambda x: x[1]))
    with open('dom_difference.txt','w') as outf:
        outf.write('#Total different domain combinations: {}\n'.format(\
            len(dif)))
        for key,val in sorted(dif.items(), key=lambda x: x[1],reverse=True):
            outf.write('{}\t{}\n'.format(key,val))
    unlist_comp = {d for dom in comp for d in dom.split(';')}
    dif_ecdm = {}
    dif_ecdm_list = defaultdict(int)
    for dom in dif:
        no_ecdm = False
        doml = dom.split(';')
        for d in doml:
            m = re.search(r'_c\d+$',d)
            if m:
                if d[:m.start()] not in ecdm:
                    no_ecdm = True
                    if d not in unlist_comp:
                        #add all counts
                        dif_ecdm_list[d]+=dif[dom]
            else:
                if d not in ecdm:
                    no_ecdm = True
                    if d not in unlist_comp:
                        dif_ecdm_list[d]+=int(dif[dom])
        if no_ecdm:
            dif_ecdm[dom] = dif[dom]
    # print(sorted(dif_ecdm.items(),key=lambda x: x[1]),len(dif),len(dif_ecdm))
    with open('dom_difference_ecdm.txt','w') as outf:
        outf.write('#Total different domain combinations: {}\n'.format(\
            len(dif_ecdm)))
        for key,val in sorted(dif_ecdm.items(), key=lambda x: x[1],\
            reverse=True):
            outf.write('{}\t{}\n'.format(key,val))
    for dm,c in sorted(dif_ecdm_list.items(),key=lambda x: x[1],reverse=True):
        print(dm+'\t'+str(c))
    print(len(dif_ecdm_list))
