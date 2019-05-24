#!/usr/bin/env python3
"""
Author: Joris Louwen
Script to find modules with LDA algorithm.
"""

from sys import argv
from multiprocessing import Pool
from functools import partial
from collections import Counter, defaultdict
import matplotlib
matplotlib.use('Agg') #to not rely on X-forwarding (not available in screen)
import matplotlib.pyplot as plt
import time
from operator import itemgetter
import os
from statistics import mean,median
import subprocess

from numpy import sqrt

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

def run_lda(domlist, no_below, no_above, num_topics, cores, outfolder, \
    ldavis=True):
    '''
    Returns LDA model with the Dictionary and the corpus, LDAvis is optional

    domlist: list of list of str, list of the bgc domain-combinations
    no_below: int, domain-combinations that occur in less than no_below
        bgcs will be removed
    no_above: float, remove domain-combinations that occur in more than
        no_above fraction of the dataset
    num_topics: int, number of topics
    cores: int, number of cores to use
    outfolder: str, filepath
    ldavis: bool, if true save LDAvis visualisation of model
    '''
    dict_lda = Dictionary(domlist)
    dict_lda.filter_extremes(no_below=no_below, no_above=no_above)
    print(dict_lda)
    corpus_bow = [dict_lda.doc2bow(doms) for doms in domlist]
    model = os.path.join(outfolder,'lda_model')
    if not os.path.exists(model):
        lda = LdaMulticore(corpus=corpus_bow, num_topics=num_topics, \
            id2word=dict_lda, workers=cores, per_word_topics=True)
        lda.save(model)
    else:
        print('Loaded existing LDA model')
        lda = LdaMulticore.load(model)
    cm = CoherenceModel(model=lda, corpus=corpus_bow, dictionary=dict_lda,\
        coherence='c_v', texts=domlist)
    coherence = cm.get_coherence()
    print('Coherence: {}, num_topics: {}'.format(coherence, num_topics))
    if ldavis:
        visname = os.path.join(outfolder,'lda.html')
        print('Running pyLDAvis for visualisation')
        vis = pyLDAvis.gensim.prepare(lda, corpus_bow, dict_lda)
        print('  saving visualisation to html')
        pyLDAvis.save_html(vis, visname)
    return lda, dict_lda, corpus_bow

def process_lda(lda, dict_lda, corpus_bow, modules, min_f_score, bgcs, \
    outfolder, amplif=False):
    '''
    '''
    # cs = []
    # ldamods = []
    num_sums = []
    select_features = [] #list with how many features to select to get to
    #a certain score like 0.3 for every topic
    out_topics = os.path.join(outfolder, 'topics.txt')
    with open(out_topics,'w') as outf:
        outf.write('Topic\tDomain_combinations\tScores\n')
        for top, mod in lda.print_topics(-1, 50):
            # modstr = {m.split('*')[1].strip('"') for m in mod.split(' + ')}
            nums = []
            doms = []
            for m in mod.split(' + '):
                num, dom = m.split('*')
                dom = dom.strip('"')
                nums.append(float(num))
                doms.append(dom)
            s=[]
            m_len = len([s.append(num) for num in nums \
                if sum(s) < min_f_score])
            select_features.append(m_len)
            num_sums.append(sum(nums))
            # ldamods.append(tuple(doms[:m_len]))
            #now all features, change to doms[:m_len] for features until min_f
            outf.write('{}\t{}\t{}\n'.format(top,','.join(doms),\
                ','.join(map(str,nums))))
            # coincide = []
            # for m in modules:
                # if modstr.intersection(m) == set(m):
                    # coincide.append((modules[m][0],m))
            # if coincide:
                # print(coincide)
                # cs.append(coincide)
    #amount of modules overlapping with statistical approach
    # mods_overlap = len({m for mods in cs for m in mods})
    # print(num_sums)
    # print(mean(num_sums),median(num_sums),max(num_sums))
    l_nonsense = len([zero for zero in num_sums if zero == 0])
    print('Length of 0 topics:',l_nonsense)
    sums = [nzero for nzero in num_sums if nzero != 0]
    # print(select_features)
    # print('Amount of modules below 10:',len([f for f in select_features if \
        # f < 10]))
    # for l,mod in zip(select_features,ldamods):
        # if l < 10:
            # print(','.join(sorted(mod)),l)
    # print(mean(sums),median(sums),max(sums))
    # print(len(cs),len({m for mods in cs for m in mods}))
    # print(\
    # 'Coherence: {}, num_topics: {}, topic_len: {}, mods_overlap: {}'.format(\
        # coherence, num_topics, topic_len, mods_overlap))
    # print(lda)
    #the per document topic matrix, a list of lists with (topic_n, score)
    # print(doc_lda[0])
    bgc2topic = link_bgc_topics(lda, dict_lda, corpus_bow, bgcs, outfolder,\
        True, amplif)
    link_genes2topic(lda, dict_lda, corpus_bow, bgcs, outfolder)
    t_matches = retrieve_topic_matches(bgc2topic)


def link_bgc_topics(lda, dict_lda, corpus_bow, bgcs, outfolder, plot=True, \
    amplif=False):
    '''Returns dict of {bgc:{'len':bgc_len,topic_num:[prob,[(gene,prob)]]}}
    
    
    Writes file to outfolder/bgc_topics.txt and supplies plots if plot=True
    '''
    print('\nLinking topics to BGCs')
    doc_lda = lda[corpus_bow]
    doc_topics = os.path.join(outfolder, 'bgc_topics.txt')
    bgc2topic = {}
    #filter the zip(bgcs,doc_lda) here
    if amplif:
        get_index = set(range(0,len(bgcs),amplif))
        bgc_docs = [(bgcs[i],doc_lda[i]) for i in get_index]
    else:
        bgc_docs = zip(bgcs,doc_lda)
    with open(doc_topics,'w') as outf:
        for bgc, doc_doms in bgc_docs:
            #doc_doms consists of three lists:
            #1 topics 2 word2topic 3 word2topic with probability
            topd = {tpc:[prob,[]] for tpc,prob in doc_doms[0]}
            for domcomb in doc_doms[2]:
                #find matching words with probabilities
                name = dict_lda[domcomb[0]]
                toptup = domcomb[1] #all topic assignments for a word
                for t in toptup:
                    #each t is a tuple of (topic, probability)
                    topd[t[0]][1].append((name,t[1]))
            outf.write('>{}\n'.format(bgc))
            outf.write('len={}\n'.format(len(doc_doms[1])))
            for top, info in sorted(topd.items(),key=lambda x: x[1][0],\
                reverse=True):
                s_genes = sorted(info[1],key=lambda x: x[1],reverse=True)
                topd[top][1] = s_genes #sort matching genes - high to low p
                genes = ','.join(['{}:{:.2f}'.format(g,p) for g,p in s_genes])
                string='topic={}\n\tp={}\n\tlen={}\n\tgenes={}\n'.format(\
                    top,info[0], len(info[1]), genes)
                outf.write(string)
            topd['len'] = len(doc_doms[1])
            bgc2topic[bgc] = topd
    if plot:
        print('Plotting')
        #extract length of each bgc vs len of topic in each bgc
        lengths = ((val['len'],len(val[t][1])) for val in\
            bgc2topic.values() for t in val if not t == 'len')
        len_counts = Counter(lengths)
        x_y, counts = zip(*len_counts.items())
        bgc_len, topic_len = zip(*x_y)
        m_counts = max(len_counts.values())
        print(m_counts)
        fig, ax = plt.subplots()
        scatter = ax.scatter(bgc_len, topic_len, c=sqrt(counts), s=2.5,vmin=1,\
            vmax=sqrt(m_counts), cmap='hot')
        leg_range = [1]+[round(x,-1) for x in \
            range(20,m_counts,int(m_counts/4))]
        if len(leg_range) == 4:
            leg_range.append(round(m_counts,-1))
        kw = dict(num=leg_range,func=lambda c: c**2)
        legend = ax.legend(*scatter.legend_elements(**kw), loc='upper left',\
            title='Occurrence')
        ax.add_artist(legend)
        plt.xlabel('Length BGC')
        plt.ylabel('Length topic match')
        plt.title('Length of a BGC vs length of matching topic')
        plt.savefig(os.path.join(outfolder,'len_bgcs_vs_len_topic_match.pdf'))

        #count amount of topics per bgc
        topics_per_bgc = Counter([len(vals)-1 for vals in bgc2topic.values()])
        xs = range(max(topics_per_bgc)+1)
        h = [topics_per_bgc[x] if x in topics_per_bgc else 0 for x in xs]
        plt.close()
        plt.bar(xs, h)
        plt.xlabel('Number of topics per BGC')
        plt.ylabel('Occurence')
        plt.title('Topics per BGC')
        plt.savefig(os.path.join(outfolder,'topics_per_bgc.pdf'))
        plt.close()
    return bgc2topic

def link_genes2topic(lda, dict_lda, corpus_bow, bgcs, outfolder):
    '''
    '''
    outfile = os.path.join(outfolder, 'terms_to_topic.txt')
    with open(outfile, 'w') as outf:
        for d_id in dict_lda:
            d_name = dict_lda[d_id]
            domc_topics = sorted(lda.get_term_topics(d_name,0.001), key=\
                lambda x: x[1], reverse=True)
            dom_top_str = '\t'.join(';'.join(map(str,d)) for d in domc_topics)
            outf.write('{}\t{}\n'.format(d_name, dom_top_str))
        #visualise amount of topics per term

def retrieve_topic_matches(bgc2topic):
    '''

    bgc2topic: dict of {bgc:{'len':bgc_len,topic_num:[prob,[(gene,prob)]]}}
    '''
    #get all topic matches per topic
    topic_matches = defaultdict(list)
    for dc in bgc2topic.values():
        for k,v in dc.items():
            topic_matches[k].append(v)
    del topic_matches['len']
    prevalence = {t:len(vals) for t,vals in topic_matches.items()}
    print(sorted(prevalence.items(),key=lambda x: x[1],reverse=True))
    print(len(prevalence))
    for match in sorted(topic_matches[6],key=lambda x: len(x[1])):
        print(match)
    # print(topic_matches[6])


if __name__ == '__main__':
    start = time.time()
    #filtered bgc csv file
    bgcfile = argv[1]
    cors = int(argv[2])
    #filtered module file
    modfile = argv[3]
    topics = int(argv[4])
    min_feat_score = float(argv[5])
    out_folder = argv[6]
    if len(argv) > 7:
        amplify = int(argv[7])
    else:
        amplify = False
    if len(argv) == 9:
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
    if not os.path.isdir(out_folder):
        subprocess.check_call('mkdir {}'.format(out_folder), shell=True)
    bgcs = remove_infr_doms_str(bgcs, 0, False)
    if amplify:
        bgc_items = []
        for bgc in bgcs.items():
            bgc_items += [bgc]*amplify
        bgclist, domlist = zip(*bgc_items)
    else:
        bgclist, domlist = zip(*bgcs.items())

    lda, lda_dict, bow_corpus = run_lda(domlist, no_below=1, no_above=0.5, \
    num_topics=topics, cores=cors, outfolder=out_folder, ldavis=vis)
    process_lda(lda, lda_dict, bow_corpus, modules, min_feat_score, \
        bgclist, out_folder, amplif=amplify)

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
