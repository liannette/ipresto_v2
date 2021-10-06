#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to get certain domain hits with their bitscores from a folder
containing domtables.
The domains to match are read from an input txt file, which has one (partial)
domain name per line. Output is a file with all the matches and their
bitscores.

Usage:
    python3 replace_domains.py <dom_names.txt> <domtab_folder> <output_folder>
'''

from glob import iglob
from itertools import groupby
from sys import argv
import os
from numpy import std, mean

def read_doms(in_file):
    '''Reads text file into list

    in_file: str, file path
    '''
    with open(in_file, 'r') as inf:
        lines = [line.strip() for line in inf]
    return lines

if __name__ == "__main__":
    dom_names_file = argv[1]
    domtab_folder = argv[2]
    out_folder = argv[3]

    domtables = iglob(os.path.join(domtab_folder,"*.domtable"))
    # ~ domtables = ["../testdata_domains/testdata_domtables/35131.assembled_unknown.cluster033.domtable"]
    doms = read_doms(dom_names_file)
    for dom in doms:
        dom_dict = {}
        outfile = os.path.join(out_folder, dom+"_bitscores.txt")
        print("Writing {}".format(outfile))
        dom += '_c'
        domtables = iglob(os.path.join(domtab_folder,"*.domtable"))
        for domtab in domtables:
            with open(domtab, 'r') as domt:
                for b, lines in groupby(domt,key=lambda l: l.startswith(dom)):
                    if b:
                        for line in lines:
                            line = line.split()
                            name = line[0]
                            score = float(line[7])
                            try:
                                dom_dict[name].append(score)
                            except KeyError:
                                dom_dict[name] = [score]
        with open(outfile, 'w') as outf:
            for k,v in sorted(dom_dict.items()):
                outf.write("{},mean={},sd={},number={}\n{}\n".format(k,
                    mean(v),std(v),len(v),','.join(map(str,v))))

