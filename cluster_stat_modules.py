#!/usr/bin/env python3
'''
Perform clustering of the statistical method modules.
Author: Joris Louwen
'''

from collections import defaultdict
import numpy as np
import os
from scipy.sparse import csr_matrix
from sklearn.feature_extraction.text import CountVectorizer
from sys import argv


if __name__ == '__main__':
    mod_info = argv[1]

    out_mods = mod_info.split('.txt')[0]+'clustering.txt'
    modules = {} #keep track of info
    corpus = []
    vectorizer = CountVectorizer()
    with open(mod_info,'r') as inf:
        #{mod_num:[info]}
        inf.readline() #header
        for line in inf:
            line = line.strip().split('\t')
            mod_num = int(line[0])
            #NB: mod_nums are 1-based while index is 0-based
            modules[mod_num] = line[1:]
            corpus.append(line[-1].replace(',',' '))
    sparse_feat_matrix = vectorizer.fit_transform(corpus)
    print(sparse_feat_matrix)

