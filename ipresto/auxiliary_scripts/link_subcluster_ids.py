#!/usr/bin/env python3
'''
Author: Joris Louwen

Fast script to turn the accessions from the SubClusterBlast into domains
with a dom_hits.txt file. Returns the SubClusterBlast file with an added
column with the domain combinations.
'''
from sys import argv
import os
from collections import defaultdict
import re

def filter_out_domains(genes,include_list):
    '''Returns the same list only keeping domains from the include_list

    genes: list of str, [domains_in_genes]
    include_list: list of str, domains to keep
    '''
    newgenes = []
    for gene in genes:
        if gene == '-':
            newgenes.append(gene)
        else:
            ngene = []
            for dom in gene.split(';'):
                #check if domain is a subPfam
                m = re.search(r'_c\d+$',dom)
                if m:
                    if dom[:m.start()] in include_list:
                        ngene.append(dom)
                else:
                    if dom in include_list:
                        ngene.append(dom)
            #if gene becomes empty print -
            if not ngene:
                newgenes.append('-')
            else:
                newgenes.append(';'.join(ngene))
    return newgenes

def read_txt(in_file):
    '''Reads text file into list

    in_file: str, file path
    '''
    with open(in_file, 'r') as inf:
        lines = [line.strip() for line in inf]
    return lines

if __name__ == '__main__':
    subclustblast = argv[1]
    dom_hits_file = argv[2]
    if len(argv) > 3:
        include_file = argv[3]
        suffix = '_domains_synt_subset.txt'
    else:
        include_file = False
        suffix = '_domains.txt'
    out = subclustblast.split('.txt')[0] + suffix

    subclusts=[]
    ids_doms = defaultdict(list)
    ids_bgcs = defaultdict(str)
    with open(subclustblast,'r') as inf:
        for line in inf:
            subclusts.append(line.strip().split('\t'))
    with open(dom_hits_file, 'r') as inf:
        inf.readline()
        for line in inf:
            line = line.strip().split('\t')
            prot_id = line[-7]
            if prot_id:
                prot_id_splt = prot_id.split('.')[0]
                ids_doms[prot_id].append(line[-3])
                ids_doms[prot_id_splt].append(line[-3])
                ids_bgcs[prot_id] = line[0]
                ids_bgcs[prot_id_splt] = line[0]
    if include_file:
        include_list = read_txt(include_file)

    subclusts_doms = []
    for subcl in subclusts:
        genes = []
        bgc=''
        for prot_id in subcl[-1].split(';')[:-1]:
            if not bgc:
                #try to get the bgc for all protein ids
                bgc = ids_bgcs[prot_id]
            doms_list = ids_doms.get(prot_id, ['-'])
            #get doms from ids_doms
            doms = ';'.join(doms_list)
            genes.append(doms)
        if include_file:
            genes = filter_out_domains(genes,include_list)
        gene_str = ','.join(genes)
        not_empty = [d for d in genes if not d == '-']
        if not not_empty:
            print('Subcluster {} has no prot_id match'.format(' '.join(\
                subcl[:2])))
            gene_str = ''
        subcl = [bgc]+subcl
        subcl.append(gene_str)
        subclusts_doms.append(subcl)

    with open(out, 'w') as outf:
        for line in subclusts_doms:
            try:
                outline = [line[i] for i in [0,2,4,-1]]
            except IndexError:
                outline = line
            outf.write('{}\n'.format('\t'.join(outline)))
