#!/usr/bin/env python3
'''
'''

from Bio import SearchIO
from sys import argv

def sign_overlap(tup1, tup2, cutoff):
    overlap = len(range(max(tup1[0], tup2[0]), min(tup1[1], tup2[1])))
    print(overlap)
    if overlap > 0:
        if overlap > min(abs(tup1[0]-tup1[1]), abs(tup2[0]-tup2[1]))*cutoff:
            return True
    return False

if __name__ == '__main__':
    domfile = argv[1]
    min_overlap = 0.1
    queries = list(SearchIO.parse(domfile, 'hmmscan3-domtab'))
    print(queries)

    for query in queries:
        dom_matches = []
        gene_name = query.id
        for hit in query:
            match = hit[0]
            #print(hit) #get Query range (corresponds to ali from to)
            domain = match.hit_id
            range_q = match.query_range
            bitsc = match.bitscore
            dom_matches.append((gene_name.split('_')[-1], domain, range_q,
                bitsc))
        print(dom_matches)
        dels = []
        if len(query) > 1:
            for i in range(len(query)-1):
                for j in range(i+1, len(query)):
                    if sign_overlap(dom_matches[i][2],dom_matches[j][2],
                        min_overlap):
                        if dom_matches[i][3] >=dom_matches[j][3]:
                            dels.append(j)
                        else:
                            dels.append(i)
        print(dels)
        doms = [dom_matches[i][1] for i in range(len(query)) if i not in dels]
        print(doms)
        #find overlap between all hits of a domain
    print(sign_overlap((1,10), (8,20), min_overlap))

