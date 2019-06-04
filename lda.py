#!/usr/bin/env python3
"""
Author: Joris Louwen
Script to find modules with LDA algorithm.
"""
import argparse
from collections import Counter, defaultdict
from functools import partial
import matplotlib
matplotlib.use('Agg') #to not rely on X-forwarding (not available in screen)
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
from numpy import sqrt
import numpy as np
from operator import itemgetter
import os
import pandas as pd
import scipy.cluster.hierarchy as sch
import seaborn as sns
from statistics import mean,median
import subprocess
from sys import argv
import time

from gensim.models.ldamulticore import LdaMulticore
from gensim.models.coherencemodel import CoherenceModel
from gensim.corpora.dictionary import Dictionary

import pyLDAvis
import pyLDAvis.gensim

def get_commands():
    parser = argparse.ArgumentParser(description="A script to cluster genes\
        from BGCs represented as strings of domains with the LDA algorithm\
        to discover sub-clusters of genes which putatively synthesise a\
        chemical moiety in the natural product.")
    parser.add_argument("-i", "--bgcfile", dest="bgcfile", help="Input \
        csv file of BGCs with genes as domain combinations", required=True)
    parser.add_argument("-m", "--modfile", dest="modfile", help="Input \
        txt file of putative modules to compare. Last column should contain\
        modules", default=False)
    parser.add_argument("-o", "--out_folder", dest="out_folder", help="Output\
        folder", required=True)
    parser.add_argument("-c", "--cores", dest="cores", help="Amount \
        of cores to use for the LDA model, default = all available cores",\
        default=cpu_count(), type=int)
    parser.add_argument("-t", "--topics", dest="topics", help="Amount \
        of topics to use for the LDA model", required=True, type=int)
    parser.add_argument("-f", "--min_feat_score", dest="min_feat_score",
        help="Only include features until their scores add up to this number.\
        Default = 0.9. Can be combined with feat_num, where feat_num features\
        are selected or features that add up to min_feat_score",type=float, \
        default=0.9)
    parser.add_argument("-n", "--feat_num", dest="feat_num",
        help="Include the first feat_num features for each topic, \
        default = 15.",type=int, default=15)
    parser.add_argument("-a", "--amplify", dest="amplify", help="Amplify \
        the dataset in order to achieve a better LDA model. Each BGC will be\
        present amplify times in the dataset. After calculating the LDA model \
        the dataset will be scaled back to normal.",type=int, default=0)
    parser.add_argument("-v", "--visualise", help="Make a visualation of the\
        LDA model with pyLDAvis (html file). If number of topics is too big\
        this might fail. No visualisation will then be made", default=False,
        action="store_true")
    parser.add_argument("--classes", help="A file containing classes of the \
        BGCs used in the analysis. First column should contain matching BGC\
        names. Consecutive columns should contain classes.", default=False)
    parser.add_argument("--plot", help="If provided: make plots about \
        several aspects of the output. Default is off.", default=False, \
        action="store_true")
    parser.add_argument("--known_subclusters", help="A tab delimited file \
    with known subclusters. Should contain subclusters in the last column and\
    BGC identifiers in the first column. Subclusters are comma separated\
    genes represented as domains. Multiple domains in a gene are separated by\
    semi-colon.")
    return parser.parse_args()

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

def process_lda(lda, dict_lda, corpus_bow, modules, feat_num, bgc_dict,
    min_f_score, bgcs, outfolder, bgc_classes, amplif=False, min_t_match=0.1,\
    min_feat_match=0.3, plot=True, known_subcl=False):
    '''Analyses the topics in the bgcs
    '''
    #this is a list of tuple (topic_num, 'features_with_scores')
    lda_topics = lda.print_topics(-1, 50)
    topic_num = len(lda_topics)
    #get the topic names from the lda html visualisation
    ldahtml = os.path.join(outfolder, 'lda.html')
    if os.path.isfile(ldahtml):
        with open(ldahtml, 'r') as inf:
            for line in inf:
                if line.startswith('var lda'):
                    lst_str = line.strip().split('"topic.order": ')[-1]
                    nums = map(int, lst_str.strip('[]};').split(', '))
                    trans = {i_lda-1:i_vis+1 for i_vis,i_lda in \
                        zip(range(topic_num), nums)}
    else:
        trans = {x:'-' for x in range(topic_num)}
    select_features = [] #list with how many features to select to get to
    #a certain score like 0.3 for every topic
    out_topics = os.path.join(outfolder, 'topics.txt')
    #to record the features as {topic:[(gene,prob)]}, features are selected
    #until the min_f_score. coul also just select a certain number like 15
    filt_features = {}
    with open(out_topics,'w') as outf:
        outf.write('Topic\tNumber_LDAvis\tDomain_combinations\tScores\n')
        for top, mod in lda_topics:
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
            if m_len > feat_num:
                filt_features[top] = set(doms[:feat_num])
            else:
                filt_features[top] = set(doms[:m_len])
            #now all features, change to doms[:m_len] for features until min_f
            outf.write('{}\t{}\t{}\t{}\n'.format(top,trans[top],\
                ','.join(doms),','.join(map(str,nums))))
    # print(select_features, mean(select_features))
    bgcl_dict = {bgc: sum(1 for g in genes if not g == '-') \
        for bgc,genes in bgc_dict.items()}
    bgc2topic = link_bgc_topics(lda, dict_lda, corpus_bow, bgcs, outfolder,\
        bgcl_dict, plot=plot, amplif=amplif)
    link_genes2topic(lda, dict_lda, corpus_bow, bgcs, outfolder)
    t_matches = retrieve_topic_matches(bgc2topic)
    top_match_f = os.path.join(outfolder,'matches_per_topic.txt')
    t_matches = write_topic_matches(t_matches, bgc_classes, top_match_f)
    t_matches = filter_matches(t_matches, filt_features, min_t_match, \
        min_feat_match)
    top_match_f_filt = top_match_f.split('.txt')[0]+'_filtered.txt'
    write_topic_matches(t_matches, bgc_classes, top_match_f_filt)
    bgc_with_topics = retrieve_match_per_bgc(t_matches, bgc_classes, \
        known_subcl,outfolder)
    if plot:
        bgc_topic_heatmap(bgc_with_topics, bgc_classes, topic_num, outfolder,\
            metric='euclidean')
        bgc_topic_heatmap(bgc_with_topics, bgc_classes, topic_num, outfolder,\
            metric='correlation')
        bgc_class_heatmap(bgc_with_topics, bgc_classes, topic_num, outfolder,\
            metric='correlation')


def link_bgc_topics(lda, dict_lda, corpus_bow, bgcs, outfolder, bgcl_dict,
    plot=True, amplif=False):
    '''Returns dict of {bgc:{topic_num:[prob,[(gene,prob)]]}}
    
    
    Writes file to outfolder/bgc_topics.txt and supplies plots if plot=True
    '''
    print('\nLinking topics to BGCs')
    doc_lda = lda[corpus_bow]
    doc_topics = os.path.join(outfolder, 'bgc_topics.txt')
    bgc2topic = {}
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
            outf.write('len={}\n'.format(bgcl_dict[bgc]))
            for top, info in sorted(topd.items(),key=lambda x: x[1][0],\
                reverse=True):
                s_genes = sorted(info[1],key=lambda x: x[1],reverse=True)
                topd[top][1] = s_genes #sort matching genes - high to low p
                genes = ','.join(['{}:{:.2f}'.format(g,p) for g,p in s_genes])
                string='topic={}\n\tp={}\n\tlen={}\n\tgenes={}\n'.format(\
                    top,info[0], len(info[1]), genes)
                outf.write(string)
            bgc2topic[bgc] = topd
    if plot:
        print('Plotting')
        #extract length of each bgc vs len of topic in each bgc
        lengths = ((bgcl_dict[bgc],len(val[t][1])) for bgc,val in\
            bgc2topic.items() for t in val)
        len_name = os.path.join(outfolder,'len_bgcs_vs_len_topic_match.pdf')
        plot_topic_matches_lengths(lengths, len_name)

        #count amount of topics per bgc
        tpb_name = os.path.join(outfolder,'topics_per_bgc.pdf')
        topics_per_bgc = Counter([len(vals) for vals in bgc2topic.values()])
        plot_topics_per_bgc(topics_per_bgc,tpb_name)
    return bgc2topic

def plot_topic_matches_lengths(lengths, outname):
    '''
    Make a scatterplot of the lengths of the topic matches vs the bgc lengths

    lengths: list of tuples, [(bgc_len,match_len)]
    outname: str, filepath
    '''
    len_counts = Counter(lengths)
    x_y, counts = zip(*len_counts.items())
    bgc_len, topic_len = zip(*x_y)
    m_counts = max(len_counts.values())
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
    plt.savefig(outname)
    plt.close()

def plot_topics_per_bgc(topics_per_bgc, outname):
    '''Make a barplot of the amount of topics per bgc

    topics_per_bgc: dict/counter object, {n:bgcs_with_n_topics}
    outname: str
    '''
    xs = range(max(topics_per_bgc)+1)
    h = [topics_per_bgc[x] if x in topics_per_bgc else 0 for x in xs]
    plt.close()
    plt.bar(xs, h)
    plt.xlabel('Number of topics per BGC')
    plt.ylabel('Occurence')
    plt.title('Topics per BGC')
    plt.savefig(outname)
    plt.close()

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
    '''Turns bgcs with matching topics to topics with matches from bgc

    bgc2topic: dict of {bgc:{'len':bgc_len,topic_num:[prob,[(gene,prob)]]}}
    topic_matches: {topic:[[prob,[(gene,prob)],bgc]]}
    '''
    #get all topic matches per topic
    topic_matches = defaultdict(list)
    for bgc,dc in bgc2topic.items():
        for k,v in dc.items():
            if not k == 'len':
                newv = v+[bgc]
                topic_matches[k].append(newv)
    return topic_matches

def retrieve_match_per_bgc(topic_matches,bgc_classes,known_subcl,outfolder):
    '''Turns topics with matches back into bgc with matches

    topic_matches: {topic:[[prob,(gene,prob),bgc]]}
    bgc_classes: {bgc:[class1,class2]}
    known_subcl: {bgc: [[info,about,subcl]]}
    bgc2topic: dict of {bgc:[[topic_num,prob,[(gene,prob)]]]}
    '''
    bgc2topic = defaultdict(list)
    for topic,info in topic_matches.items():
        for match in info:
            bgc2topic[match[2]].append([topic]+match[:2])
    with open(os.path.join(outfolder, 'bgc_topics_filtered.txt'),'w') as outf:
        for bgc,info in sorted(bgc2topic.items()):
            outf.write('>{}\nclass={}\n'.format(bgc,\
                bgc_classes.get(bgc,['None'])[0]))
            if known_subcl:
                if bgc in known_subcl:
                    for i,subcl in enumerate(known_subcl[bgc]):
                        outf.write('known_sub-cluster={}\n'.format(', '.join(\
                            subcl)))
            for match in sorted(info, key=lambda x: x[1],reverse=True):
                outf.write('{}\t{}\t{}\n'.format(match[0],match[1],\
                ','.join(['{}:{:.2f}'.format(m[0],m[1]) for m in match[2]])))
    return bgc2topic

def write_topic_matches(topic_matches, bgc_classes, outname):
    '''Writes topic matches to a file sorted on length and alphabet

    topic_matches: {topic:[[prob,(gene,prob),bgc]]}
    bgc_classes: {bgc: [class1,class2]}
    outname: str, filepath
    '''
    plotlines = []
    s_b_c = set(bgc_classes)
    #occurence of each topic
    prevl = {t:len(vals) for t,vals in topic_matches.items()}
    sumfile = outname.split('.txt')[0]+'_summary.txt'
    with open(outname,'w') as outf, open(sumfile,'w') as sumf:
        sumf.write('Topic\tmatches\tmatches_len>1\tclasses\n')
        for topic, matches in sorted(topic_matches.items()):
            classes = Counter([bgc_classes.get(bgc,['None'])[0] for bgc in \
                list(zip(*matches))[2]])
            class_str = ','.join([':'.join(map(str,cls)) for cls in \
                sorted(classes.items())])
            prevl = len(matches)
            prevl_bigger_1 = len([m for m in matches if len(m[1]) > 1])
            #topicnr #matches #matches>1 classes
            outf.write('#Topic {}, matches:{}, matches_len>1:{}, '+
                'classes:{}\n'.format(topic,prevl,prevl_bigger_1, class_str))
            sum_line = [topic, prevl, prevl_bigger_1, class_str]
            plotlines.append(sum_line)
            sumf.write('{}\n'.format('\t'.join(map(str,sum_line))))
            #sort the matches by length and then by alphabet
            if matches:
                sorted_matches = sorted(matches,key=lambda x: \
                    (len(x[1]),list(zip(*x[1]))[0]))
                topic_matches[topic] = sorted_matches
                for match in sorted_matches:
                    outf.write('{:.3f}\t{}\t{}\t{}\n'.format(match[0], ','.join(\
                    ['{}:{:.2f}'.format(m[0],m[1]) for m in match[1]]),\
                    match[2],bgc_classes.get(match[2],['None'])[0]))
    return topic_matches

def filter_matches(topic_matches, topic_features, min_t_match,min_feat_match):
    '''Filters topic_matches based on cutoffs

    topic_matches: {topic:[[prob,(gene,prob)],bgc]}, topic linked to matches
    topic_features: {topic:set(genes)}, dict of features to use
    min_t_match: float, minimal score of a topic matching a bgc
    min_feat_match: float, minimal score of a feature matching in a topic in
        a bgc
    '''
    filt_topic_matches = defaultdict(list)
    for topic, matches in topic_matches.items():
        # filt_topic_matches[topic] = []
        use_feats = topic_features[topic]
        for match in matches:
            match_p = match[0]
            bgc = match[2]
            if match_p > min_t_match:
                newfeats = []
                for feat in match[1]:
                    if feat[0] in use_feats and feat[1] >= min_feat_match:
                        newfeats.append(feat)
                if newfeats:
                    filt_topic_matches[topic].append([match_p,newfeats,bgc])
    return filt_topic_matches

def bgc_topic_heatmap(bgc_with_topic, bgc_classes, topic_num, outfolder,
    metric='euclidean'):
    '''Make a clustered heatmap of bgcs and topics, and optional bgc_classes

    bgc_with_topic: dict of {bgc:[[topic_num,prob,[(gene,prob)]]]}
    bgc_classes: dict of {bgc:[[class1,class2]]}
    topic_num: int, number of topics in the model
    
    '''
    print('\nMaking clustered heatmap, metric: {}'.format(metric))
    #make pd dataframe from bgc with topic with prob as value for present tpic
    bgcs, topics = zip(*bgc_with_topic.items())
    data = [{v[0]:v[1] for v in val} for val in topics]
    df = pd.DataFrame(data,index=bgcs,columns=list(range(topic_num)))
    df = df.fillna(0)
    #colour rows by bgc class
    class_set = set(bgc_classes.keys())
    labels = [bgc_classes[bgc][0] if bgc in class_set else 'None' for bgc \
        in bgcs]
    s_labels = sorted(set(labels))
    #get colours
    if 'None' in s_labels:
        s_labels.remove("None")
    if len(s_labels) > 10:
        lut = dict(zip(s_labels, sns.cubehelix_palette(len(\
            s_labels),start=1.2,rot=2,dark=0.11,light=0.85)))
    else:
        lut = dict(zip(s_labels, sns.color_palette()))
    lut['None'] = 'w' #make None always white
    s_labels = ['None']+s_labels
    row_labs = pd.DataFrame(labels,index=bgcs,columns=['BGC classes'])
    row_colours = row_labs['BGC classes'].map(lut) #map colour to a label

    g = sns.clustermap(df, cmap = 'nipy_spectral', row_colors = row_colours, \
        linewidths = 0, metric=metric, yticklabels=False, xticklabels=True, \
        cbar_kws = {'orientation':'horizontal'},vmin=0,vmax=1)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(),\
        fontsize = 5)
    #don't show dendrograms
    g.ax_col_dendrogram.set_visible(False)
    g.ax_row_dendrogram.set_ylim([0,0.00001])
    g.ax_row_dendrogram.set_xlim([0,0.00001])
    #make legend for classes
    for label in s_labels:
        g.ax_row_dendrogram.bar(0,0,color=lut[label], label=label,linewidth=0)
    g.ax_row_dendrogram.legend(loc="center left",fontsize='small',\
        title='BGC classes')
    #move colourbar
    g.cax.set_position([.35, .78, .45, .0225])
    plt.savefig(\
        os.path.join(outfolder, 'topic_heatmap_{}.pdf'.format(metric)))
    plt.close()

def bgc_class_heatmap(bgc_with_topic, bgc_classes, topic_num, outfolder,
    metric='correlation'):
    '''Make a clustered heatmap of bgcs and topics, and optional bgc_classes

    bgc_with_topic: dict of {bgc:[[topic_num,prob,[(gene,prob)]]]}
    bgc_classes: dict of {bgc:[[class1,class2]]}
    topic_num: int, number of topics in the model
    
    '''
    print('\nMaking clustered heatmap of classes, metric: {}'.format(metric))
    #make pd dataframe from bgc with topic with prob as value for present tpic
    bgcs, topics = zip(*bgc_with_topic.items())
    data = [{v[0]:v[1] for v in val} for val in topics]
    df = pd.DataFrame(data,index=bgcs,columns=list(range(topic_num)))
    df = df.fillna(0)
    #colour rows by bgc class
    class_set = set(bgc_classes.keys())
    labels = [bgc_classes[bgc][0] if bgc in class_set else 'None' for bgc \
        in bgcs]
    s_labels = sorted(set(labels))
    #cluster each class (hierarchical, correlation)
    class_i = clust_class_bgcs(df, labels, s_labels)
    #get colours
    if 'None' in s_labels:
        s_labels.remove("None")
    if len(s_labels) > 10:
        lut = dict(zip(s_labels, sns.cubehelix_palette(len(\
            s_labels),start=1.2,rot=2,dark=0.11,light=0.85)))
    else:
        lut = dict(zip(s_labels, sns.color_palette()))
    lut['None'] = 'w' #make None always white
    s_labels = ['None']+s_labels
    row_labs = pd.DataFrame(labels,index=bgcs,columns=['BGC classes'])
    row_colours = row_labs.iloc[class_i,0].map(lut) #map colour to a label

    g = sns.clustermap(df.iloc[class_i,:], cmap = 'nipy_spectral', \
        row_colors = row_colours, linewidths = 0, metric=metric, \
        yticklabels=False, xticklabels=True, cbar_kws = \
        {'orientation':'horizontal'},vmin=0,vmax=1, row_cluster=False)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(),\
        fontsize = 5)
    #don't show dendrograms
    g.ax_col_dendrogram.set_visible(False)
    g.ax_row_dendrogram.set_ylim([0,0.00001])
    g.ax_row_dendrogram.set_xlim([0,0.00001])
    #make legend for classes
    for label in s_labels:
        g.ax_row_dendrogram.bar(0,0,color=lut[label], label=label,linewidth=0)
    g.ax_row_dendrogram.legend(loc="center left",fontsize='small',\
        title='BGC classes')
    #move colourbar
    g.cax.set_position([.35, .78, .45, .0225])
    plt.savefig(\
        os.path.join(outfolder, 'class-topic_heatmap_{}.pdf'.format(metric)))
    plt.close()

def clust_class_bgcs(df, labels, s_labels):
    '''Returns a list of indeces ordered on clustered classes
    '''
    #get a list of clustered indexes for all and then add them
    inds = np.array([],dtype='int32')
    for bgc_class in s_labels:
        c_i = [i for i,cls in enumerate(labels) if cls == bgc_class]
        dist = sch.distance.pdist(df.iloc[c_i,:], metric = 'correlation')
        clust = sch.linkage(dist, metric='correlation')
        ind = sch.leaves_list(clust)
        # print(ind)
        ind_reorder = [c_i[i] for i in ind]
        inds = np.append(inds,ind_reorder)
    return inds

def read2dict(filepath, sep=','):
    '''Read file into a dict {first_column:[other_columns]}

    filepath: str
    sep: str, delimiter in the file
    '''
    output = {}
    with open(filepath,'r') as inf:
        for line in inf:
            line = line.strip().split(sep)
            output[line[0]] = line[1:]
    return output

if __name__ == '__main__':
    start = time.time()
    #files provided should be filtered bgc csv file and filtered module file

    print('\nStart')
    cmd = get_commands()
    bgcs = read2dict(cmd.bgcfile)
    with open(cmd.modfile, 'r') as inf:
        modules = {}
        #{modules:[info]}
        for line in inf:
            line = line.strip().split('\t')
            mod = tuple(line[-1].split(',')) #now a tuple of str
            modules[mod] = line[:-1]
    if cmd.classes:
        bgc_classes_dict = read2dict(cmd.classes, sep='\t')
    else:
        bgc_classes_dict = {bgc:'None' for bgc in bgcs}
    if not os.path.isdir(cmd.out_folder):
        subprocess.check_call('mkdir {}'.format(cmd.out_folder), shell=True)

    bgcs = remove_infr_doms_str(bgcs, 0, False)
    if cmd.amplify:
        bgc_items = []
        for bgc in bgcs.items():
            bgc_items += [bgc]*cmd.amplify
        bgclist, domlist = zip(*bgc_items)
    else:
        bgclist, domlist = zip(*bgcs.items())

    if cmd.known_subclusters:
        known_subclusters = read2dict(cmd.known_subclusters,sep='\t')
    else:
        known_subclusters = False

    lda, lda_dict, bow_corpus = run_lda(domlist, no_below=1, no_above=0.5, \
        num_topics=cmd.topics, cores=cmd.cores, outfolder=cmd.out_folder, \
        ldavis=cmd.visualise)
    process_lda(lda, lda_dict, bow_corpus, modules, cmd.feat_num, bgcs,
        cmd.min_feat_score, bgclist, cmd.out_folder, bgc_classes_dict, \
        amplif=cmd.amplify, plot=cmd.plot, known_subcl=known_subclusters)

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
