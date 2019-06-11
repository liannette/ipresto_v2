#!/usr/bin/env python3
"""
Author: Joris Louwen
Script to parse JSON files into a clusterfile: a file where BGCs are
represented by domains. Domains are kept in their gene structure: each ,
separates a gene while the domains within a gene are separated by ;.
"""

import argparse
from glob import iglob
from multiprocessing import Pool,cpu_count
import os
import subprocess
import json

def get_commands():
    parser = argparse.ArgumentParser(description="A script to parse JSON\
        files into a csv file containing each BGC and its domains.")
    parser.add_argument("-i", "--in_folder", dest="in_folder", help="Input \
        folder with JSON files", required=True)
    parser.add_argument("-o", "--out_folder", dest="out_folder", help="Output\
        file", required=True)
    parser.add_argument("-p","--prefix", dest="prefix",help="If provided,\
        only files starting with prefix in the in_folder will be taken into \
        account",default="")
    parser.add_argument("-d", "--domain_overlap_cutoff", 
        dest="domain_overlap_cutoff", default=0.1, help="Specify at which \
        overlap percentage domains are considered to overlap. Domain with \
        the best score is kept (default=0.1).")
    return parser.parse_args()

def sign_overlap(tup1, tup2, cutoff):
    '''
    Returns true if there is an overlap between two ranges higher than cutoff

    tup1, tup2: tuples of two ints, start and end of alignment
    cutoff: float, fraction that two alignments are allowed to overlap

    Overlap is be calculated with the smallest domain alignment to be strict
    '''
    overlap = len(range(max(tup1[0], tup2[0]), min(tup1[1], tup2[1])))
    if overlap > 0:
        if overlap > min(abs(tup1[0]-tup1[1]), abs(tup2[0]-tup2[1]))*cutoff:
            return True
    return False

def parse_json_wrapper(in_folder,out_folder,min_overlap):
    '''Parses all json files to a csv clusterfile and a dom_hits file

    file_path: str, path to file
    '''
    

def parse_json(file_path,min_overlap):
    '''Parses a json file to a list of lists

    file_path: str, path to file
    min_overlap: float, the amount of overlap two domains are allowed to have
    bgc_info_sorted: list of lists: [bgc, g_id, p_id, loc, orf_num, tot_orf,
        dom, range, bitscore]
    '''
    if not os.path.isfile(file_path):
        print('Ignored {} file does not exist'.format(file_path))
        return
    bgc_info = []
    with open(file_path, 'r') as inf:
        doc = json.load(inf)
    bgc = doc['name']
    genes_doc = doc['clusters'][0]['genes']
    tot_genes = len(genes_doc)
    for gene_doc in genes_doc:
        gene = [bgc]
        for x in ['locus_tag','protein_id']:
            gene.append(gene_doc['info'][x])
        for x in ['start','end']:
            gene.append(gene_doc[x])
        pfams_doc = gene_doc['pfams']
        if not pfams_doc:
            gene += ['-',0,0,'']
            bgc_info.append(gene)
        else:
            pfams = []
            for pfam in pfams_doc:
                try:
                    #take only the best subPfam hit
                    sub = pfam['subdomains'][0]
                except (KeyError,IndexError):
                    hit_info = [pfam[x] for x in \
                        ['name','hit_start','hit_end','bitscore']]
                else:
                    hit_info = [sub[x] for x in \
                        ['name','start','end','bitscore']]
                pfams.append(hit_info)
            #correct for overlapping domains
            dels = []
            if len(pfams) > 1:
                for i in range(len(pfams)-1):
                    for j in range(i+1, len(pfams)):
                        #if there is a significant overlap delete the one with
                        #the lower bitscore
                        if sign_overlap(pfams[i][1:3],pfams[j][1:3],
                            min_overlap):
                            if pfams[i][2] >= pfams[j][2]:
                                dels.append(j)
                            else:
                                dels.append(i)
            bgc_info += [gene+pfams[i] for i in range(len(pfams)) if not i in\
                dels]
    #sort for original gene/domain order
    #end result bgc_info: bgc g_id p_id loc orf_num tot_orf dom range bitscore
    #in this case g_id is locus_tag
    gene_num = 0
    prev= ''
    bgc_info_sorted = []
    for bgc_i in sorted(bgc_info,key=lambda x: (x[3],x[6])):
        gene_curr = bgc_i[1]
        if gene_curr != prev:
            gene_num+=1
            prev = gene_curr
        new_bgc = bgc_i[:3] + [';'.join(map(str,bgc_i[3:5])),gene_num,\
            tot_genes,bgc_i[5],';'.join(map(str,bgc_i[6:8])),bgc_i[8]]
        bgc_info_sorted.append(new_bgc)
    return bgc_info_sorted

if __name__ == "__main__":
    cmd = get_commands()
    file_path = os.path.join(cmd.in_folder,\
        'antismashdb_v2_all.NC_014500.1.cluster006.json')
    if not os.path.isdir(cmd.out_folder):
        subprocess.check_call('mkdir {}'.format(cmd.out_folder))
    cluster = parse_json(file_path, cmd.domain_overlap_cutoff)
    parse_json_wrapper(cmd.in_folder,cmd.out_folder,cmd.domain_overlap_cutoff)
