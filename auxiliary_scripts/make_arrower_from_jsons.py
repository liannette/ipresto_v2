#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to convert fastas and dom_hits file into json files that can be
used for visualisations of bgcs. A nice extension would be to include
which genes belong to putative modules and somehow visualise that too.

Usage:
python3 make_arrower_from_jsons.py -h

Example usage:
python3 make_arrower_from_jsons.py -i bgc_name_file -p pfam_list_file
    -o out_name -j json_folder -t template_html --pfam_json pfam_json.json

'''

import argparse
from collections import OrderedDict, Counter, defaultdict
# from copy import deepcopy
# from functools import partial
from glob import glob, iglob
from itertools import combinations, product, islice, chain,groupby
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
    parser = argparse.ArgumentParser(description="A script to draw an html\
        representations of BGCs colouring only certain domains")
    parser.add_argument("-i", "--in_file", help="input file containing list\
        of bgcs to draw", required=True)
    parser.add_argument("-p","--pfam_list", help="file containing (sub)Pfams\
        to use")
    parser.add_argument("-o", "--out_file", dest="out_file", 
        required=True, help="extensionless utput file, also html title")
    parser.add_argument("-j","--json_folder",help='folder containing json \
        files')
    parser.add_argument("-t","--template_html",help='template html to build\
        on')
    parser.add_argument("--pfam_json",help='json file containing domains\
        with colouring info')
    return parser.parse_args()

def read_txt(file_name):
    '''
    '''
    out_list = []
    with open(file_name,'r') as inf:
        for line in inf:
            out_list.append(line.strip())
    out_list.reverse()
    return out_list

def write_html(template, outname, pfam_file):
    '''
    '''
    with open(template,'r') as inf:
        html = inf.readlines()
    # title = os.path.split(outname)[1]
    string_list = [pfam_file,outname+'_data.js',outname]
    for i,string in zip([4,6,19],string_list):
        html[i] = html[i].format(os.path.split(string)[1])
    htmlfile = outname+'.html'
    with open(htmlfile,'w') as outf:
        for line in html:
            outf.write(line)

if __name__ == "__main__":
    cmd = get_commands()
    bgcs = read_txt(cmd.in_file)
    doms = read_txt(cmd.pfam_list)
    data_js = cmd.out_file + '_data.js'
    with open(data_js, 'w') as outf:
        outf.write('var bs_data=')
        json_list = []
        for bgc in bgcs:
            filename = os.path.join(cmd.json_folder,bgc+'.json')
            with open(filename, 'r') as inf:
                json_file = json.load(inf)
                json_list.append(json_file)
        json.dump(json_list,outf)
        outf.write(';')
    pfam_js_file = cmd.out_file + '_pfams.js'
    with open(cmd.pfam_json , 'r') as inf, open(pfam_js_file,'w') as outf:
        pfams = json.load(inf)
        newpfams = {dom:pfams[dom] for dom in doms}
        outf.write('var pfams=')
        json.dump(newpfams,outf)
        outf.write(';')
    write_html(cmd.template_html,cmd.out_file,pfam_js_file)
