#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to get the TC scores for certain domains from the Pfam hmm database
and then write the TC scores in the subdomain files for those domains

Usage:
    python3 replace_domains.py <main_hmm> <get_domains.txt> <sub_dom_hmm>

Notes:
right all three scores are added
'''

from sys import argv
from itertools import groupby
from glob import iglob
import os
from subprocess import check_call

def read_doms(in_file):
    '''Reads text file into list

    in_file: str, file path
    '''
    with open(in_file, 'r') as inf:
        lines = [line.strip() for line in inf]
    return lines

def get_scores(in_file, doms):
    '''
    '''
    score_dict = {}
    with open(in_file, 'r') as inf:
        #returns a bool in boolval and an iterable in lines
        for boolval, lines in groupby(inf, lambda line:\
            line.startswith('HMMER3/f')):
            if not boolval:
                record = [line for line in lines]
                name = record[0].rstrip('\n').split()[1]
                if name in doms:
                    print("Getting score from {}".format(name))
                    cksum = ''
                    for i in range(25): #should be enough
                        #find CKSUM row
                        l_id = record[i].strip().split()[0]
                        if l_id == 'CKSUM':
                            cksum = i+1
                            break
                    if not cksum:
                        print('Warning: no CKSUM row in {}'.format(name))
                    score_dict[name] = record[cksum:cksum+3]
    return score_dict

def update_entries(in_file, scores, out_folder):
    '''Writes a new file with elements from scores added in in_file

    in_file: string, file path
    scores: dict of list of str {dom: ['line with score','..'],..}
    '''
    out_name = os.path.split(in_file)[1][:-4]+'_addscore.hmm'
    out_file = os.path.join(out_folder, out_name)
    dom = os.path.split(in_file)[1].split('-subdomains.hmm')[0]
    with open(in_file, 'r') as inf, open(out_file, 'w') as outf:
        #returns a bool in boolval and an iterable in lines
        for boolval, lines in groupby(inf, lambda line:\
            line.startswith('HMMER3/f')):
            if boolval:
                header = next(lines) #will only contain 1 line
            else:
                record = [line for line in lines]
                cksum = ''
                for i in range(25): #should be enough
                    #find CKSUM row
                    l_id = record[i].strip().split()[0]
                    if l_id == 'CKSUM':
                        cksum = i+1
                        break
                if not cksum:
                    cksum = 9 #just in case CKSUM is not there
                newrecord = record[:cksum]+scores[dom]+record[cksum:]
                outf.write(header)
                outf.write('{}'.format(''.join(newrecord)))

if __name__ == "__main__":
    hmm_file = argv[1]
    dom_file = argv[2]
    sub_folder = argv[3]
    out_folder = argv[4]
    if not os.path.isdir(out_folder):
        check_call('mkdir {}'.format(out_folder), shell = True)

    dom_list = read_doms(dom_file)
    scores = get_scores(hmm_file, dom_list)
    sub_files = iglob(os.path.join(sub_folder, '*subdomains.hmm'))
    for sub_file in sub_files:
        update_entries(sub_file, scores, out_folder)
