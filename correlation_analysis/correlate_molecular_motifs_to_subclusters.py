#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to correlate subclusters to molecular motifs in a group of strains.
FDR is estimated with target-decoy approach.
Matrix tsv files look like:
id\tStrain_1\tStrain2\n
motif1\t0\t1\n
motif2\t1\t0\n

Usage:
python3 correlate_molecular_motifs_to_subclusters.py -h

Example usage:
python3 correlate_molecular_motifs_to_subclusters.py -m molecular_motifs.tsv
    -s subcluster_motifs.tsv -o out_file
'''

import argparse
from collections import OrderedDict, Counter, defaultdict
from itertools import groupby, product
from math import floor, log10
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
import numpy as np
from operator import itemgetter
import os
import random
import re
import scipy.stats as scs
import seaborn as sns
import subprocess
import time

def get_commands():
    parser = argparse.ArgumentParser(description="A script to correlate\
        subcluster motifs to molecular motifs that are present/absent in\
        a group of organisms. FDR can be estimated with a target-decoy\
        approach, a maximum possible score is calculated and a permutation\
        test is performed.")
    parser.add_argument("-m", "--molecular_motifs", help="Input file which is\
        a presence/absence matrix for molecular motifs in strains in tsv\
        format. The first line contains id followed by strain names\
        (columns). The rows are the motifs, each value is either 0 or 1",\
        required=True)
    parser.add_argument('-s','--subcluster_motifs',help="Input file which is\
        a presence/absence matrix for subcluster motifs in strains in tsv\
        format. The first line contains id followed by strain names\
        (columns). The rows are the motifs, each value is either 0 or 1",\
        required=True)
    parser.add_argument('-o','--out_file',help='Output file',required=True)
    parser.add_argument('-f','--fdr_cutoff',help='FDR cutoff in percentage,\
        default = 1',type=float, default=1)
    parser.add_argument('-v','--verbose',action='store_true',default=False,\
        help='Print additional info')
    parser.add_argument('--filter',default=False, action='store_true',\
        help='Filter out motifs that occur in more than -r of the dataset')
    parser.add_argument('-r','--remove',default=0.5, help='Percentage of\
        dataset that will be taken for removing motifs, default = 0.5',\
        type=float)
    parser.add_argument('-c','--cores', help='Number of cores to use',
        default=5,type=int)
    parser.add_argument('-n','--n_resamples',help='The number of permutations\
        taken, default=500',default=500,type=int)
    parser.add_argument('--annotation_m',help='Annotation csv file for\
        molecular motifs, optional. Contains a header, first column motif\
        names, last column annotation',default=False)
    parser.add_argument('--annotation_s',help='Annotation csv file for\
        subcluster motifs, optional. Contains a header, first column motif\
        names, last column annotation',default=False)
    return parser.parse_args()

def read_matrix(infile, remove=0.5, filtering=False):
    '''
    Returns motif matrix {strain:set(present_motifs)} and set(motifs)

    infile: str, filepath
    filtering: bool, filter out motifs occurring in more than half the strains
    '''
    print('\nReading motif file from {}'.format(infile))
    strain2motif = defaultdict(set)
    rownames = set()
    with open(infile,'r') as inf:
        colnames = inf.readline().strip().split('\t')[1:] #ignore 'id'
        for row in inf:
            row = row.strip().split('\t')
            motif = row.pop(0)
            rownames.add(motif)
            for presence,col in zip(row,colnames):
                if presence == '1':
                    strain2motif[col].add(motif)
    if filtering:
        strain2motif,rownames = filter_motifs(strain2motif,\
            remove)
    filtered_strain2motif,rownames = filter_strains(strain2motif)
    print('  filtered out {} strains containing only one motif'.format(\
        len(strain2motif) - len(filtered_strain2motif)))
    print("Motif matrix contains {} strains, {} motifs and {} 1's".format(\
        len(filtered_strain2motif), len(rownames), sum(\
        len(vals) for vals in filtered_strain2motif.values())))
    return filtered_strain2motif, rownames

def filter_strains(motif_dict):
    '''Returns same dict but strains are removed if they contain < 2 motifs

    motif_dict: {strain:set(present_motifs)}
    '''
    filtered_motifs = {strain:motifs for strain,motifs in motif_dict.items()\
        if len(motifs) > 1}
    rownames = {val for values in filtered_motifs.values() for val in values}
    return filtered_motifs,rownames

def filter_motifs(motif_dict,remove):
    '''
    Returns same dict, motifs removed if they occurr in > remove% of strains

    motif_dict: {strain:set(present_motifs)}
    remove: remove motifs that occur in more than remove*len(strains)
    keep_motifs: set of str, motifs to keep
    '''
    motif_counts = Counter([m for motifs in motif_dict.values() for m in\
        motifs])
    remove_len = len(motif_dict)*remove
    keep_motifs = {m for m,count in motif_counts.items() if count<=remove_len}
    print('  filtered out {} motifs that occur in > {} of strains'.format(\
        len(motif_counts.keys())-len(keep_motifs),remove))
    filtered_motifs = {strain : set(\
        [m for m in motifs if motif_counts[m] <= remove_len])\
        for strain,motifs in motif_dict.items()
        }
    return filtered_motifs, keep_motifs

def make_scoring_matrix(m_motifs, s_motifs, m_used, s_used,\
    strains_used):
    '''Returns dict of {m_motif: {s_motif:score} }

    m_motifs, s_motifs: {strain:[present_motifs]} m for molecular-and s
        for subcluster motifs
    m_used, s_used: set of str, all used motif names for molecular-and
        subcluster motifs respectively
    strains_used: list of str, all available strains to choose from

    Score is +10 is both present in a strain, +1 if both absent, 0 if the
    s_motif is there but not m_motif, -10 if m_motif is there but not s_motif
    '''
    #keep track of score as dict of dict {m_motif:{s_motif:score}}
    scoring_matrix = {m_m:defaultdict(int) for m_m in m_used}
    for strain in strains_used:
        molec_strain = m_motifs[strain]
        subcl_strain = s_motifs[strain]
        not_molec = m_used - molec_strain
        not_subcl = s_used - subcl_strain
        for both in product(molec_strain, subcl_strain):
            scoring_matrix[both[0]][both[1]] += 10
        for only_m in product(molec_strain, not_subcl):
            scoring_matrix[only_m[0]][only_m[1]] -= 10
        #only_s not necessary as it is 0
        for neither in product(not_molec,not_subcl):
            scoring_matrix[neither[0]][neither[1]] += 1
    return scoring_matrix

def get_max_scores(m_motifs, s_motifs, m_m_names, s_m_names,\
    strains_used):
    '''Returns dict of {m_motif: {s_motif:max_score} }

    m_motifs, s_motifs: {strain:[present_motifs]} m for molecular-and s
        for subcluster motifs
    m_m_names, s_m_names: set of str, all used motif names for molecular-and
        subcluster motifs respectively
    strains_used: list of str, all available strains to choose from
    '''
    #first loop through the matrix to get all strains per motif
    m_motif_dict = defaultdict(int)
    s_motif_dict = defaultdict(int)
    strains_len = len(strains_used)
    for strain in strains_used:
        m_mots = m_motifs[strain]
        s_mots = s_motifs[strain]
        for m_mot in m_mots:
            m_motif_dict[m_mot] += 1
        for s_mot in s_mots:
            s_motif_dict[s_mot] += 1
    max_scores_dict = defaultdict(dict)
    i=0
    for combi in product(m_motif_dict.items(), s_motif_dict.items()):
        (m,num_m),(s,num_s) = combi
        ind, max_num = [(i,num) for i,num in enumerate([num_m,num_s]) if\
            num == max(num_m,num_s)][0]
        plus_one = strains_len - max_num
        if ind == 0:
            minus_ten = num_m - num_s
            plus_ten = num_s
        else:
            minus_ten = 0
            plus_ten = num_m
        max_sc = 1 * plus_one + 10 * plus_ten - 10 * minus_ten
        max_scores_dict[m][s] = max_sc
        i+=1
    return max_scores_dict


def create_decoy_matrix(motif_matrix, motif_names, strains_used, verbose):
    '''Returns a scrambled version of motif_matrix

    motif_matrix: {strain:set(present_motifs)}
    motif_names: set of str, all used motif names
    strains_used: list of str, all available strains to choose from

    Scrambled matrix will contain same amount of strains per motif, but
    just randomly chosen from the strains, so scrambling the presence absence
    in across strains but keeping same presence/absence distribution.
    '''
    #scramble all the ones across strains and collect in decoy
    decoy = defaultdict(set)
    #first loop through the matrix to get all strains per motif
    motif_dict = defaultdict(set) #set so they can be compared later
    motif_list = []
    strain_lens = {}
    for strain in strains_used:
        motifs = motif_matrix[strain]
        strain_lens[strain] = len(motifs)
        for motif in sorted(motifs): #sort to make deterministic
            motif_dict[motif].add(strain)
            if not motif in motif_list:
                motif_list.append(motif)
    random.shuffle(motif_list)
    for motif in motif_list:
        target_strains = motif_dict[motif]
        length = len(target_strains)
        #choose the same amount of random strains
        random.shuffle(strains_used)
        scramble = set(random.sample(strains_used,length))
        #make sure target vector does not get in decoy matrix
        while target_strains == scramble:
            scramble = set(random.sample(strains_used,length))
        if verbose:
            print(motif, length, len(scramble),len(target_strains & scramble))
        for decoy_strain in scramble:
            decoy[decoy_strain].add(motif)
    return decoy


def plot_scoring_matrix(target, max_len, decoy = False, plot_name = False,
    cutoff_score = False, fdr= False, above_fdr=False):
    '''Save/show a plot of target distribution overlayed with decoy with FDR
    '''
    #get all values and convert to numpy array
    scores = np.array(list(zip(*target))[2])
    #could be that some zero values are not in matrix as i didn't initialise
    #does not matter much but for correctness
    if len(scores) != max_len:
        difference = max_len - len(scores)
        scores = np.append(scores, [0]*difference)
    sns.set_style('darkgrid')
    ax = sns.distplot(scores,kde_kws={'color':'#004B96','label':'Target'},\
        hist_kws= {'color':'#004B96','alpha':0.4})
    if decoy:
        s_decoy = np.array(list(zip(*decoy))[2])
        if len(s_decoy) != max_len:
            print(len(s_decoy),len(scores))
            difference = max_len - len(s_decoy)
            s_decoy = np.append(s_decoy, [0]*difference)
        assert len(s_decoy)==len(scores)
        sns.distplot(s_decoy,kde_kws={'color':'#DC3220','label':'Decoy',\
            'bw':.08}, hist_kws= {'color':'#DC3220','alpha':0.4})
    if cutoff_score:
        ax.axvline(cutoff_score, color = 'black', ls = '-.', lw=0.5)
        bot, top = ax.get_ylim()
        ax.text(x = cutoff_score+5,y=top/2, s='{} above\n{}% FDR'.format(\
            above_fdr,fdr), multialignment='left',size=8)
    plt.title('Target and decoy distribution of correlation scores')
    plt.xlabel('Correlation scores')
    plt.ylabel('Density')
    if plot_name:
        plt.savefig(plot_name)
    else:
        plt.show()
    plt.close()

def get_random_scores(m_motifs, m_used, s_motifs, s_used, strains_used,\
    verbose, cores, loops = 500):
    '''Returns {m_motif: {s_motif:[list_of_random_scores]} }

    m_motifs, s_motifs: {strain:set(present_motifs)}
    m_used, s_used: set of motifs used
    strains_used: set of strains used
    verbose: bool, print extra information
    loops: int, number of loops
    '''
    print('\nResampling matrix',loops,'times to get random distribution')
    all_matrix = {m_m:defaultdict(list) for m_m in m_used}
    pool = Pool(cores,maxtasksperchild=10)
    done_list = []
    for new_seed in random.sample(range(100000),k=loops):
        pool.apply_async(get_random_score,args=(m_motifs, m_used, s_motifs,\
            s_used, strains_used, verbose,new_seed),\
            callback=lambda x: done_list.append(x))
    pool.close()
    pool.join()
    for decoy_matrix in done_list:
        for m_m,s_ms in decoy_matrix.items():
            for s_m,score in s_ms.items():
                all_matrix[m_m][s_m].append(score)
    return all_matrix

def get_random_score(m_motifs, m_used, s_motifs, s_used, strains_used,\
    verbose,newseed):
    '''Returns {m_motif: {s_motif:random_score} }

    m_motifs, s_motifs: {strain:set(present_motifs)}
    m_used, s_used: set of motifs used
    strains_used: set of strains used
    verbose: bool, print extra information
    newseed: a seed number to randomise the decoy
    '''
    random.seed(newseed) #make random (there is a seed for making one
    #decoy so choose a new seed each time here
    mol_decoy = create_decoy_matrix(m_motifs,\
        list(m_used),strains_used,verbose)
    sub_decoy = create_decoy_matrix(s_motifs,\
        list(s_used),strains_used,verbose)
    decoy_matrix = make_scoring_matrix(mol_decoy, sub_decoy,\
        m_used, s_used, strains_used)
    return decoy_matrix

def compare_target_to_random(target, list_decoy_matrix, m_used,s_used,\
    n_resamples=500):
    '''Returns {m_motif: {s_motif:pvalue} }

    target: {m_motif: {s_motif:score} }
    list_decoy_matrix: {m_motif: {s_motif:[list_of_random_scores]} }
    m_used, s_used: set of motifs used
    p-value comes from one sided resampling test (amount of values from
        resampled distribution larger or equal to target value divided by
        amount of resamplings), so basically checking if target value is on
        edge of distribution
    '''
    all_matrix = defaultdict(dict)
    for m_m,s_ms in target.items():
        for s_m,score in s_ms.items():
            distr = list_decoy_matrix[m_m][s_m]
            above_or_equal = len([0 for i in distr if i >= score])
            pval = above_or_equal / n_resamples if above_or_equal > 0 else \
                1/n_resamples
            all_matrix[m_m][s_m] = pval
    return all_matrix

def read_annotation(infile):
    '''Returns dict with {motif_name:'annotation'}, annotation is last column

    infile: str, file path to csv
    '''
    annot_dict = {}
    with open(infile,'r') as inf:
        inf.readline() #header
        for line in inf:
            line = line.strip('\r\n').split('","')
            if len(line) == 1:
                line = line[0].split(',') #account for quotes
            mot = line[0].strip('"')
            ann = line[-1].strip('"')
            if ann == 'None':
                ann = ''
            annot_dict[mot] = ann
    return annot_dict

def correlation_analysis(molecular_infile, sub_cluster_infile, outfile,\
    filtering, remove, fdr, verbose, cores, resamples_n, annotation_m,\
    annotation_s):
    '''Combines all functions to make plot of target-decoy distr and outfile

    molecular_infile, sub_cluster_infile, outfile: str, filepaths to matrix
        files and to output file
    filtering: bool, to filter motifs occurring in 'remove' or more strains
    fdr: float, FDR cutoff in percentage
    verbose: bool, print more info
    cores: int, cores to use
    resamples_n: int, rounds of resampling
    annotation_m,annotation_s: str, filepaths to annotation files
    '''
    molecular_motifs, m_motif_names = read_matrix(molecular_infile,remove,\
        filtering)
    subcluster_motifs, s_motif_names = read_matrix(sub_cluster_infile,\
        remove, filtering)

    #only compare strains present in both matrices and motifs in these strains
    used_strains = sorted(set(molecular_motifs) & set(subcluster_motifs))
    m_motif_names = {m for strain, motifs in molecular_motifs.items() for m\
        in motifs if strain in used_strains}
    s_motif_names = {m for strain, motifs in subcluster_motifs.items() for m\
        in motifs if strain in used_strains}

    target_matrix = make_scoring_matrix(molecular_motifs, subcluster_motifs,\
        m_motif_names, s_motif_names, used_strains)
    max_scores = get_max_scores(molecular_motifs, subcluster_motifs,\
        m_motif_names, s_motif_names, used_strains)

    print('\nScrambling motif matrix')
    molecular_decoy = create_decoy_matrix(molecular_motifs,\
        list(m_motif_names),used_strains,verbose)
    subcluster_decoy = create_decoy_matrix(subcluster_motifs,\
        list(s_motif_names),used_strains,verbose)
    decoy_matrix = make_scoring_matrix(molecular_decoy, subcluster_decoy,\
        m_motif_names, s_motif_names, used_strains)

    #Permutation test by simulating each target distribution by constructing
    #many decoys to compare target value to
    multiple_decoy_matrix = get_random_scores(molecular_motifs,\
        m_motif_names, subcluster_motifs, s_motif_names, used_strains,verbose,\
        cores,loops=resamples_n)
    target_pvals = compare_target_to_random(target_matrix,\
        multiple_decoy_matrix,m_motif_names,s_motif_names,resamples_n)

    #get all values and convert tuples
    target_tuples = sorted([(mm,sm,score) for mm, val_dict in \
        target_matrix.items() for sm,score in val_dict.items()],\
        key=lambda x: -x[2])
    decoy_tuples = sorted([(mm,sm,score) for mm, val_dict in \
        decoy_matrix.items() for sm,score in val_dict.items()],\
        key=lambda x: -x[2])
    if verbose:
        print(target_tuples[:10])
        print(decoy_tuples[:10])

    cutoff_score, len_above_cutoff = calc_fdr_score(target_tuples,\
        decoy_tuples, fdr)

    #get annotations if there are any
    if annotation_m:
        annot_m_m = read_annotation(annotation_m)
    else:
        annot_m_m = {}
    if annotation_s:
        annot_s_m = read_annotation(annotation_s)
    else:
        annot_s_m = {}

    with open(outfile,'w') as outf:
        outf.write('#{}% FDR cutoff: {}. {} values above cutoff\n'.format(\
            fdr, cutoff_score, len_above_cutoff))
        outf.write('#{} molecular-and {} subcluster motifs used'.format(\
            len(m_motif_names),len(s_motif_names)) +\
            ' across {} strains\n'.format(len(used_strains)))
        header = ['Molecular_motif', 'Subcluster_motif', 'Score',\
            'Max_score', 'Ratio_max_score','P-value_permutation_test']
        outf.write('{}\n'.format('\t'.join(header)))
        ratios = []
        for tup in target_tuples:
            m_m,s_m,sc = tup
            max_s = max_scores[m_m][s_m]
            ratio = sc/max_s if max_s > 0 else 0
            ratio = ratio if ratio > 0 else 0
            ratios.append(ratio)
            pval = target_pvals[m_m][s_m]
            m_m_an = annot_m_m.get(m_m,'')
            s_m_an = annot_s_m.get(s_m,'')
            outf.write('{}\t{}\t{:.2f}\t{:.3f}\t{}\t{}\n'.format(\
                '\t'.join(map(str,tup)), max_s, ratio, pval, m_m_an, s_m_an))
    decoy_out = outfile.split('.txt')[0] + '_decoy_scores.txt'
    with open(decoy_out,'w') as outf:
        for tup in decoy_tuples:
            outf.write('{}\n'.format('\t'.join(map(str,tup))))

    # sns.set_style('darkgrid')
    # ax = sns.distplot(ratios,kde_kws={'color':'#004B96','label':\
        # 'Ratio_of_max_score'}, hist_kws= {'color':'#004B96','alpha':0.4})
    # plt.show()

    max_length = len(m_motif_names) * len(s_motif_names)
    plot_file = outfile.split('.txt')[0] + '.pdf'
    plot_scoring_matrix(target_tuples, max_length, decoy_tuples, plot_file,
        cutoff_score, fdr, len_above_cutoff)
    n_plot_file = outfile.split('.txt')[0] + '_only_target.pdf'
    plot_scoring_matrix(target_tuples, max_length, plot_name=n_plot_file)

def calc_fdr_score(targets, decoys, fdr = 1):
    '''Calculates the correlation score cutoff where there is FDR of fdr%

    targets, decoys: list of tuples [(m_motif,s_motif,score)] of target and
        decoy distribution respectively
    '''
    #print also the amount of targets above cutoff
    scores = list(zip(*targets))[2]
    s_decoy = list(zip(*decoys))[2]
    if max(scores) <= max(s_decoy):
        print('\nFDR is 100% at maximum value of target distribution')
        return max(s_decoy),0
    fdr_score = {} # keeps track of fdr for each score {fdr:which_score}
    decoy_select = [value for value in s_decoy if value >= 0]
    scores_select = [value for value in scores if value >= 0]
    for scr in range(max(scores)+1):
        decoy_select = [value for value in decoy_select if value >= scr]
        scores_select = [value for value in scores_select if value >= scr]
        false_pos_estimate = len(decoy_select)
        len_values = len(scores_select)
        fdr_estimate = false_pos_estimate / len_values * 100
        if not fdr_estimate in fdr_score:
            #keep lowest score associated to a FDR
            fdr_score[fdr_estimate] = (scr,len_values)
    # print(sorted(fdr_score.items(),key=lambda x: x[1][0]))
    #find cutoff
    prev = 0
    for fdr_sc in sorted(fdr_score):
        print(fdr_sc, fdr_score[fdr_sc])
        if fdr_sc > fdr:
            chosen_cutoff = prev
            break
        prev = fdr_sc
    print(chosen_cutoff, fdr_score[chosen_cutoff])
    return fdr_score[chosen_cutoff]

if __name__ == '__main__':
    start = time.time()
    random.seed(595) #get same output every time
    cmd = get_commands()

    correlation_analysis(cmd.molecular_motifs,cmd.subcluster_motifs,\
        cmd.out_file, cmd.filter, cmd.remove, cmd.fdr_cutoff, cmd.verbose,\
        cmd.cores, cmd.n_resamples, cmd.annotation_m, cmd.annotation_s)

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
