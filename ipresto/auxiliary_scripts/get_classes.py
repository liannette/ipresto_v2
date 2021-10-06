#!/usr/bin/env python3
'''
Author: Joris Louwen

Fast script to get classes from network_annotations.tsv from big-scape output
'''
from sys import argv
import os

if __name__ == '__main__':
    net_ann_file = argv[1]
    out = os.path.join(os.path.split(net_ann_file)[0],'big-scape_classes.txt')
    with open(net_ann_file,'r') as inf, open(out,'w') as outf:
        for line in inf:
            line = line.strip().split('\t')
            outf.write('{}\t{}\n'.format(line[0],line[4]))
