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
from operator import itemgetter

from gensim.models.ldamulticore import LdaMulticore
from gensim.corpora.dictionary import Dictionary

if __name__ == '__main__':
    bgcfile = argv[1]
    cors = int(argv[2])
    modfile = argv[3]
    print('\nStart')
    with open(bgcfile, 'r') as inf:
        bgcs = {}
        for line in inf:
            line = line.strip().split(',')
            bgcs[line[0]] = line[1:]
    with open(modfile, 'r') as inf:
        modules = {}
        #{modules:[info]}
        for line in inf:
            line = line.strip().split('\t')
            mod = tuple(line[-1].split(','))
            modules[mod] = line[:-1]
    domlist = list(bgcs.values())
    lda_dict = Dictionary(domlist)
    print(lda_dict)
    #Could also use rem_infr_doms, which filters on absolute occurence
    #instead of the amount of documents containing the dom
    lda_dict.filter_extremes(no_below=3, no_above=0.5)
    print(lda_dict)
    corpus_bow = [lda_dict.doc2bow(doms) for doms in domlist]
    # print(corpus_bow)
    lda = LdaMulticore(corpus=corpus_bow, num_topics=100, id2word=lda_dict,\
        workers=cors)
    print(lda)
    cs = []
    for top, mod in lda.print_topics(-1, 15):
        modstr = {m.split('*')[1].strip('"') for m in mod.split(' + ')}
        coincide = []
        for m in modules:
            if modstr.intersection(m) == set(m):
                coincide.append((modules[m][0],m))
        if coincide:
            print('{}:\n{}'.format(top,mod))
            print(coincide)
            cs.append(coincide)
    print(len(cs))
    
    
