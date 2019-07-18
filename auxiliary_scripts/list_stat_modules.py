#!/usr/bin/env python3
'''
Make a list of statistical method modules like for lda.
Author: Joris Louwen
'''

from collections import defaultdict
import os
from sys import argv

if __name__ == '__main__':
    bgc_with_mods = argv[1]
    mod_info = argv[2]

    out_bgc = bgc_with_mods.split('.txt')[0]+'_list.txt'
    out_mods = mod_info.split('.txt')[0]+'_with_bgcs.txt'
    bgc_dct = {}
    mod_dct = defaultdict(list)
    with open(bgc_with_mods,'r') as inf:
        for line in inf:
            line = line.strip('\n').split('\t')
            bgc = line[0]
            mod_nums = line[-1].split('; ')
            empty = ['']
            if mod_nums != empty:
                mod_nums = list(map(int,mod_nums))
                bgc_dct[bgc] = mod_nums
                for mod_num in mod_nums:
                    mod_dct[mod_num].append(bgc)

    with open(mod_info,'r') as inf:
        modules = {}
        #{mod_num:[info]}
        inf.readline() #header
        for line in inf:
            line = line.strip().split('\t')
            # mod = tuple(line[-1].split(',')) #now a tuple of str
            modules[int(line[0])] = line[1:]

    with open(out_bgc,'w') as outf:
        for bgc,mod_nums in sorted(bgc_dct.items()):
            outf.write('>{}\n'.format(bgc))
            for mod_num in sorted(mod_nums):
                try:
                    mod_info = modules[mod_num]
                except KeyError:
                    #this module was filtered out somewhere else
                    pass
                else:
                    outf.write('{}\t{}\n'.format(mod_num,'\t'.join(mod_info)))

    with open(out_mods,'w') as outf:
        for mod_num,bgcs in sorted(mod_dct.items()):
            try:
                mod_info = modules[mod_num]
            except KeyError:
                pass
            else:
                outf.write('#Statistical-module {}, N_Occurences {},'.format(\
                    mod_num, mod_info[0])+
                    ' N_Genes {}, N_Domains {}, Strictest_detection_'.format(\
                    mod_info[1],mod_info[2])+'cutoff {}'.format(mod_info[3]))
                if len(mod_info) == 6:
                    outf.write(', subcluster-family {}\n'.format(mod_info[5]))
                elif len(mod_info) == 7:
                    outf.write(', subcluster-clan {}\n'.format(mod_info[6]))
                else:
                    outf.write('\n')
                outf.write('#{}\n'.format(mod_info[4]))
                for bgc in sorted(bgcs):
                    outf.write('{}\n'.format(bgc))
