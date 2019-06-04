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

if __name__ == '__main__':
    subclustblast = argv[1]
    dom_hits_file = argv[2]
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
        gene_str = ','.join(genes)
        not_empty = [d for d in genes if not d == '-']
        if not not_empty:
            print('Subcluster {} has no prot_id match'.format(' '.join(\
                subcl[:2])))
            gene_str = ''
        subcl = [bgc]+subcl
        subcl.append(gene_str)
        subclusts_doms.append(subcl)

    out = subclustblast.split('.txt')[0] + '_domains.txt'
    with open(out, 'w') as outf:
        for line in subclusts_doms:
            outf.write('{}\n'.format('\t'.join(line)))
