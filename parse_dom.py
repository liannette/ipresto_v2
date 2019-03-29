#!/usr/bin/env python3
'''
Author: Joris Louwen

Contains functions to parse a sigle domtable into a clus csv table
'''

from Bio import SearchIO
from sys import argv
import os

def parse_domtab(domfile, clus_file, min_overlap):
    '''Parses domtab into a cluster domain file (csv)

    domfile: string, file path
    clus_file: open file for writing
    Clus1,dom1,dom2,-(gene without domain)\\nClus2,dom1..
    '''
    queries = SearchIO.parse(domfile, 'hmmscan3-domtab')
    cds_before = 0
    cluster_doms = []
    for query in queries:
        dom_matches = []
        cds_num = query.id.split('_')[-1]
        for hit in query:
            match = hit[0]
            #print(hit) #get Query range (corresponds to ali from to)
            domain = match.hit_id
            range_q = match.query_range
            bitsc = match.bitscore
            dom_matches.append((domain, range_q, bitsc))
        dels = []
        if len(query) > 1:
            for i in range(len(query)-1):
                for j in range(i+1, len(query)):
                    if sign_overlap(dom_matches[i][1],dom_matches[j][1],
                        min_overlap):
                        if dom_matches[i][2] >=dom_matches[j][2]:
                            dels.append(j)
                        else:
                            dels.append(i)
        cds_doms = [dom_matches[i][0] for i in range(len(query)) if i not in dels]
        #find overlap between all hits of a domain
        #if a cds has no domains print - in output
        gene_gap = cds_num - cds_before -1
        if gene_gap > 0:
            cds_doms = ['-']*gene_gap + cds_doms
        cluster_doms += cds_doms
        cds_before = cds_num
    clus_file.write('{},{}\n'.format(\
        os.path.split(domfile)[-1].split('.domtable')[0],
        ','.join(cluster_doms)))


def sign_overlap(tup1, tup2, cutoff):
    '''Returns true if there is an overlap between two ranges higher than cutoff

    tup1, tup2: tuples of two ints, start and end of alignment
    cutoff: float, fraction that two alignments are allowed to overlap

    Overlap is be calculated with the smallest domain alignment to be strict
    '''
    overlap = len(range(max(tup1[0], tup2[0]), min(tup1[1], tup2[1])))
    print(overlap)
    if overlap > 0:
        if overlap > min(abs(tup1[0]-tup1[1]), abs(tup2[0]-tup2[1]))*cutoff:
            return True
    return False

if __name__ == '__main__':
    domfile = argv[1]
    outfile = argv[2]
    min_overlap = 0.1
    with open(outfile) as out:
        parse_domtab(domfile, out, min_overlap)

    

