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
import os
from statistics import mean,median

from gensim.models.ldamulticore import LdaMulticore
from gensim.models.coherencemodel import CoherenceModel
from gensim.corpora.dictionary import Dictionary

import pyLDAvis
import pyLDAvis.gensim

def remove_infr_doms_str(clusdict, m_gens, verbose):
    '''Returns clusdict with genes replaced  with - if they occur < 3

    clusdict: dict of {cluster:[domains_of_a_gene]}
    m_gens: int, minimal distinct genes a cluster must have to be included
    verbose: bool, if True print additional info

    Deletes clusters with 1 unique gene
    '''
    print('\nRemoving domain combinations that occur less than 3 times')
    domcounter = Counter()
    domcounter.update([v for vals in clusdict.values() for v in vals \
        if not v == '-'])
    deldoms = [key for key in domcounter if domcounter[key] <= 2]
    print('  {} domain combinations are left, {} are removed'.format(\
        len(domcounter.keys())-len(deldoms),len(deldoms)))
    clus_no_deldoms = {}
    for k,v in clusdict.items():
        newv = ['-' if dom in deldoms else dom for dom in v]
        doml = len({v for v in newv if not v == '-'})
        if doml >= m_gens:
            clus_no_deldoms[k] = newv
        else:
            if verbose:
                print('  {} removed as it has less than min_genes'.format(k))
    print(' {} clusters have less than {} domains and are excluded'.format(\
        len(clusdict.keys()) - len(clus_no_deldoms), m_gens))
    return clus_no_deldoms

def run_lda(dict_lda, domlist, no_below, no_above, num_topics, modules,cores,\
    min_f_score, outfolder, ldavis=True):
    '''
    '''
    dict_lda.filter_extremes(no_below=no_below, no_above=no_above)
    print(dict_lda)
    corpus_bow = [dict_lda.doc2bow(doms) for doms in domlist]
    model = os.path.join(outfolder,'lda_model')
    if not os.path.exists(model):
        lda = LdaMulticore(corpus=corpus_bow, num_topics=num_topics, \
            id2word=dict_lda, workers=cores, per_word_topics=True)
        lda.save(model)
    else:
        lda = LdaModel.load(model)
    # cs = []
    ldamods = []
    num_sums = []
    select_features = [] #list with how many features to select to get to
    #a certain score like 0.3 for every topic
    for top, mod in lda.print_topics(-1, 50):
        modstr = {m.split('*')[1].strip('"') for m in mod.split(' + ')}
        nums = []
        doms = []
        for m in mod.split(' + '):
            num, dom = m.split('*')
            dom = dom.strip('"')
            nums.append(float(num))
            doms.append(dom)
        # print('{}:\n{}'.format(sum(nums),mod))
        # print(sum(nums))
        s=[]
        m_len = len([s.append(num) for num in nums if sum(s) < min_f_score])
        select_features.append(m_len)
        num_sums.append(sum(nums))
        ldamods.append(tuple(doms[:m_len]))
        # coincide = []
        # for m in modules:
            # if modstr.intersection(m) == set(m):
                # coincide.append((modules[m][0],m))
        # if coincide:
            # print(coincide)
            # cs.append(coincide)
    #amount of modules overlapping with statistical approach
    # mods_overlap = len({m for mods in cs for m in mods})
    cm = CoherenceModel(model=lda, corpus=corpus_bow, dictionary=dict_lda,\
        coherence='c_v', texts=domlist)
    coherence = cm.get_coherence()
    print('Coherence: {}, num_topics: {}'.format(coherence, num_topics))
    # print(num_sums)
    # print(mean(num_sums),median(num_sums),max(num_sums))
    l_nonsense = len([zero for zero in num_sums if zero == 0])
    print('Length of 0 topics:',l_nonsense)
    sums = [nzero for nzero in num_sums if nzero != 0]
    print(select_features)
    print('Amount of modules below 10:',len([f for f in select_features if \
        f < 10]))
    for l,mod in zip(select_features,ldamods):
        if l < 10:
            print(','.join(sorted(mod)),l)
    # print(mean(sums),median(sums),max(sums))
    # print(len(cs),len({m for mods in cs for m in mods}))
    # print(\
    # 'Coherence: {}, num_topics: {}, topic_len: {}, mods_overlap: {}'.format(\
        # coherence, num_topics, topic_len, mods_overlap))
    # print(lda)
    #the per document topic matrix, a list of lists with (topic_n, score)
    doc_lda = lda[corpus_bow]
    print(doc_lda[0])
    if ldavis:
        visname = os.path.join(outfolder,'lda.html')
        vis = pyLDAvis.gensim.prepare(lda, corpus_bow, dict_lda)
        pyLDAvis.save_html(vis, visname)

if __name__ == '__main__':
    #filtered bgc csv file
    bgcfile = argv[1]
    cors = int(argv[2])
    #filtered module file
    modfile = argv[3]
    topics = int(argv[4])
    min_feat_score = float(argv[5])
    out_folder = argv[6]
    if len(argv) == 8:
        vis = True
    else:
        ldavisname = None
        vis = False

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
            mod = tuple(line[-1].split(',')) #now a tuple of str
            modules[mod] = line[:-1]
    bgcs = remove_infr_doms_str(bgcs, 0, False)
    domlist = list(bgcs.values())
    lda_dict = Dictionary(domlist)
    # print(lda_dict)
    run_lda(lda_dict, domlist, no_below=1, no_above=0.5, num_topics=topics,\
            modules=modules, cores=cors, min_f_score=min_feat_score, \
            outfolder=out_folder, ldavis=vis)
    
