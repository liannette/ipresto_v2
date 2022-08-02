#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to convert fastas and dom_hits file into json files that can be
used for visualisations of bgcs. A nice extension would be to include
which genes belong to putative modules and somehow visualise that too.

Usage:
python3 convert2json.py -h

Example usage:
python3 convert2json.py -d dom_hits.txt -f fasta_folder/ -o output_folder

'''

import argparse
from collections import OrderedDict, Counter, defaultdict
# from copy import deepcopy
# from functools import partial
from glob import glob, iglob
from itertools import combinations, product, islice, chain
import json
# from math import floor, log10
from multiprocessing import Pool, cpu_count
# from operator import itemgetter
import os
# import random
# from statsmodels.stats.multitest import multipletests
import subprocess
import time

def get_commands():
    parser = argparse.ArgumentParser(description="A script to turn bgcs in \
        fasta files and dom_hits file into json files used for visualisation")
    parser.add_argument("-d", "--dom_hits", help="Input file of domain hits",
        required=True)
    parser.add_argument("-f","--fasta_folder", help="Folder of fasta files")
    parser.add_argument("-o", "--out_folder", dest="out_folder", 
        required=True, help="Output directory, this will contain all output \
        data files.")
    parser.add_argument("-c", "--cores", dest="cores", default=cpu_count(), 
        help="Set the number of cores the script may use (default: use all \
        available cores)", type=int)
    return parser.parse_args()

def parse_fasta(fasta_file):
    '''Parse the headers from one fasta file into a dict

    fasta_file: str, file path
    bgc_dict = {'desc':'possible description','end':int,'id':'bgc_name',
        'orfs':[{'domains':[{'bitscore':float,'code':'pfam','end':int,
        'start':int}],"end":int,"id":"protein_id","start":int,"strand":1/-1}]
        "start":1}
    '''
    if not os.path.isfile(fasta_file):
        raise SystemExit('{} does not exist'.format(fasta_file))
    bgc_name = os.path.split(fasta_file)[1].split('.fasta')[0]
    bgc_dict = {'desc':'','end':1,'id':bgc_name,'orfs':[],'start':1}
    with open(fasta_file, 'r') as inf:
        for line in inf:
            if line.startswith('>'):
                bgc,line = line.strip('\n').split('_gid')
                line = line.split('_')
                gid,pid,loc,num = [elem.split(':')[-1] for elem in line]
                if not pid:
                    pid = num.split('/')[0]
                #< could still be in there, this is corrected in newer version
                start,end,strand = loc.replace('<','').split(';')
                if strand == '+':
                    strand = 1
                else:
                    strand = -1
                orf_dict = {'domains':[],"end":int(end),"id":pid,"start":\
                    int(start), "strand":strand}
                bgc_dict['orfs'].append(orf_dict)
    #arbitrary value for the ending of the cluster
    bgc_dict['end'] = int(end)+100
    # print(bgc_dict)
    return(bgc_dict)

def parse_dom_hits(dom_hits_file):
    '''Parses dom hits file to dict

    dom_hits_file: str, file path
    dom_hits: {BGC:{location_start:[{'bitscore':float,'code':'pfam','end':int,
        'start':int}]}}
    '''
    if not os.path.isfile(dom_hits_file):
        raise SystemExit('{} does not exist'.format(dom_hits_file))
    dom_hits = defaultdict(dict)
    with open(dom_hits_file, 'r') as inf:
        inf.readline()
        for line in inf:
            line = line.strip('\n').split('\t')
            #i assume that the start of the gene is unique
            start_gene = int(line[3].replace('<','').split(';')[0])
            st,end = map(int,line[-2].split(';'))
            dom_dict = {'bitscore':float(line[-1]),'code':line[-3],'start':st,
                'end':end}
            try:
                dom_hits[line[0]][start_gene].append(dom_dict)
            except KeyError:
                dom_hits[line[0]][start_gene] = [dom_dict]
    return dom_hits

def make_jsons(fasta_folder, out_folder, dom_hits_file):
    '''
    Get bgc_dicts from fastas, fill domains in from dom_hits and write jsons

    '''
    if not os.path.isdir(fasta_folder):
        raise SystemExit('{} does not exist'.format(fasta_folder))
    fastas = iglob(os.path.join(fasta_folder,'*.fasta'))
    bgc_dicts = (parse_fasta(fasta) for fasta in fastas)
    dom_hits = parse_dom_hits(dom_hits_file)
    for bgc_dict in bgc_dicts:
        bgc_name = bgc_dict['id']
        for i in range(len(bgc_dict['orfs'])):
            start = bgc_dict['orfs'][i]['start']
            try:
                dom_list = dom_hits[bgc_name][start]
            except KeyError:
                continue
            else:
                bgc_dict['orfs'][i]['domains'] = dom_list
        out = os.path.join(out_folder,bgc_name+'.json')
        with open(out,'w') as outf:
            json.dump(bgc_dict, outf)

if __name__ == "__main__":
    print('Start')
    start=time.time()
    #parse fasta files to dict format for visualisation with empty domains
    #add domains from parsed dom_hits file
    #for ending point of a gbk just stop at last gene? or maybe add 100
    #neater to directly parse gbks maybe
    cmd = get_commands()
    if not os.path.isdir(cmd.out_folder):
        subprocess.check_call('mkdir {}'.format(cmd.out_folder),shell=True)
    make_jsons(cmd.fasta_folder,cmd.out_folder,cmd.dom_hits)

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
