#!/usr/bin/env python3
"""
Author: Joris Louwen
Script to find modules with LDA algorithm.
"""

from sys import argv
from multiprocessing import Pool
from functools import partial
from collections import Counter
import time

from gensim.models import LdaModel
from gensim.corpora.dictionary import Dictionary

if __name__ == '__main__':
    bgcfile = argv[1]
    with open(bgcfile, 'r') as inf:
        bgcs = {}
        for line in inf:
            line = line.strip().split(',')
            bgcs[line[0]] = line[1:]
    domlist = list(bgcs.values())
    lda_dict = Dictionary(domlist)
    print(lda_dict)
    corpus = [lda_dict.doc2bow(doms) for doms in domlist]
    print(corpus)
