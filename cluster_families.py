#!/usr/bin/env python3
'''
Find clans of families of statistical method modules.
Author: Joris Louwen
'''
import os
import argparse
from collections import defaultdict, Counter
from functools import partial
from itertools import groupby, combinations, repeat
from multiprocessing import Pool
import numpy as np
import scipy.sparse as sp
from sklearn.cluster import KMeans, DBSCAN
import time


def get_commands():
    parser = argparse.ArgumentParser(description="A script to group the\
        families of the statistical modules in clans (related families) using\
        a soergel distance metric (weighted jaccard distance)")
    parser.add_argument("-i", "--infile", help="Input file with modules\
        grouped by family, headerlines are # and each module is on a separate\
        line where the domains are the last element on the tsv line", \
        required=True)
    parser.add_argument('-l','--list_file',help='File where modules are\
        listed which will be written again but with clans added',
        required=True)
    parser.add_argument('-c', '--cores', help='Cores to use, default = 1',
        default=1, type=int)
    parser.add_argument('--cutoff', help='Cutoff for when to stop saving\
        distance between two families in distance matrix, default = 0.9999',\
        default=0.9999, type=float)
    parser.add_argument('--dbscan_cutoff',help='Distance cutoff for dbscan,\
        default = 0.95',
        type=float, default = 0.95)
    return parser.parse_args()

def read_families(infile, occ=False):
    '''
    Read family file to three dicts: header_descr, fams with modules, mod_info

    infile: str, filepath
    family_dict: dict of {family_number:len_family}
    feat_dict: dict of {family_number:{feat1:score_feature,
        feat2:score_feature} } weights for each feature (relative abundance of
        a domain-combinations in a family)
    family_modules: dict of {family_number:[module_tuples]}
    modules_info: dict of {module_tup:[module_info]}
    '''
    family_dict = {}
    feat_dict = {}
    family_modules = defaultdict(list)
    modules_info = {}
    with open(infile, 'r') as inf:
        #bval is true if lines are header, false otherwise
        for bval, lines in groupby(inf,key=lambda line: line.startswith('#')):
            if bval:
                #three lines: description, occurences, relative-abundance
                desc = next(lines).strip().split(' ')
                num = int(desc[1].strip(','))
                len_fam = desc[2]
                family_dict[num] = int(len_fam)
                if occ:
                    feats = [dom.split(':') for dom in \
                        next(lines).strip().split('#Occurrences: ')[1].split(', ')]
                    tot_occ = sum(map(int,list(zip(*feats))[1]))
                    feat_dict[num] = {f:int(score)/tot_occ for f,score in feats}
                else:
                    next(lines) #discard occurences
                    feats = (dom.split(':') for dom in \
                        next(lines).strip().split('#Features: ')[1].split(', '))
                    feat_dict[num] = {f:float(score) for f,score in feats}
            else:
                for line in lines:
                    line = line.strip().split('\t')
                    mod = tuple(line[-1].split(','))
                    modules_info[mod] = line[:-1]
                    family_modules[num].append(mod)
    return family_dict,feat_dict,family_modules,modules_info

def calc_soergel_dist(in_tuple, cutoff, jacc=False):
    '''
    Returns tuple of (fam1,fam2,dist), soergel distance between two families

    in_tuple contains (fam1,feat_dict1,fam2,feat_dict2):
    -fam1,fam2: int, family numbers
    -feat_dict1, feat_dict2: {feat1:score_feature, feat2:score_feature, ..}
        weights for each feature (relative abundance of a domain-combinations
        in a family)
    '''
    fam1,feat_dict1,fam2,feat_dict2 = in_tuple
    s1 = set(feat_dict1)
    s2 = set(feat_dict2)
    overl = s1 & s2
    if jacc:
        jacc = 1 - (len(overl) / len(s1 | s2))
        if jacc <= cutoff:
            return (fam1,fam2,jacc)
        else:
            return
    mins = []
    maxs = []
    #calc min and max for each overlapping feature
    for dom in overl:
        both_scores = (feat_dict1[dom],feat_dict2[dom])
        mins.append(min(both_scores))
        maxs.append(max(both_scores))
    non_over1_score = [feat_dict1[feat] for feat in s1 if feat not in overl]
    non_over2_score = [feat_dict2[feat] for feat in s2 if feat not in overl]
    #divide sum of mins by sum of (maxs + scores of non_overlapping features)
    numerator = sum(mins)
    denominator = sum(maxs + non_over1_score + non_over2_score)
    soerg = 1 - (numerator / denominator)
    if soerg <= cutoff:
        return (fam1,fam2,soerg)

def calc_all_dists(feat_dict,cutoff,cores):
    '''
    Calculate all distances and returns tuples [(fam_num1,fam_num2,distance)]

    feat_dict: dict of {family_number:{feat1:score_feature,
        feat2:score_feature} }
    cutoff: float, cutoff for distance
    cores: int, amount of cores
    '''
    pairs = combinations(feat_dict.keys(),2)
    pair_tups = ((p1,feat_dict[p1],p2,feat_dict[p2]) for p1,p2 in pairs)
    distances = []
    #use imap as resulting list will be big
    pool = Pool(cores,maxtasksperchild=100)
    distances = pool.imap(partial(calc_soergel_dist,cutoff=cutoff),pair_tups,\
        chunksize = 10000)
    # for pair_tup in pair_tups:
        # distances.append(calc_soergel_dist(pair_tup, cutoff))
    distances = [dist for dist in distances if dist]
    return distances

def make_dist_matrix(dist_tuples, len_rows, square=False):
    '''Returns a distance matrix

    dist_tuples: list of tuples [(fam_num1,fam_num2,distance)]
    len_rows: amount of rows, so the matrix has correct shape
    '''
    rows = []
    cols = []
    dists = []
    for tup in dist_tuples:
        rows.append(tup[0])
        cols.append(tup[1])
        dists.append(tup[2])
        if square:
            rows.append(tup[1])
            cols.append(tup[0])
            dists.append(tup[2])
    sp_mat = sp.csr_matrix((dists,(rows,cols)),shape=(len_rows,len_rows))
    return sp_mat

def run_dbscan(sparse_m, dist_cutoff, rownames, list_file, prefix,\
    cores, family_modules, module_info):
    '''
    Nearest neighbour algorithm with soergel distance

    sparse_m: csr_matrix, distance matrix of shape(n_samples, n_samples)
    dist_cutoff: float, between 0 and 1 that denotes when fams are
        in a cluster based on distance
    rownames: list of ints, [mod_nums], sequential fam_nums, keeping track of
        rows of sparse_m
    prefix: str, prefix of outfile
    cores: int, amount of cores to use
    '''
    print('\nRunning DBSCAN')
    clustering = DBSCAN(eps=dist_cutoff, metric='precomputed',\
        min_samples=2, n_jobs=cores).fit(sparse_m)
    print(clustering)
    labels = clustering.labels_
    n_labels = len(set(labels))
    n_noise = list(labels).count(-1)
    if n_noise:
        n_labels -= 1
    n_clans = n_labels + n_noise-1 # #_clusters + #_singletons = #_clans
    print('Found {} clusters and {} noise points'.format(n_labels,n_noise))
    print('  writing {} clans'.format(n_clans))

    # #measure for average dist in cluster vs average dist to nearest clust
    # #1 is best -1 is worst
    # sil = silhouette_score(dist_m, labels, metric='precomputed')
    # print("  silhouette Coefficient: {0.3f}".format(sil))

    #outfiles
    dbscan_pre = '_dbscan{}_{}_clans'.format(dist_cutoff,n_clans)
    out_mods = prefix + dbscan_pre + '.txt'
    out_clusts = prefix + dbscan_pre + '_by_clan.txt'
    out_fams = prefix + dbscan_pre + '_to_family.txt'

    #link each fam to a clan/cluster
    cluster_dict = defaultdict(list)
    fam_dict = {}
    added_label = 0
    max_l = max(labels)
    for fam,cl in zip(rownames,labels):
        if cl == -1:
            #make unique clan for each noise point
            cl = max_l + added_label + 1
            added_label += 1
        cluster_dict[cl].append(fam)
        fam_dict[fam] = cl

    with open(out_fams,'w') as outf:
        for cl,fams in sorted(cluster_dict.items()):
            outf.write('>{}\n{}\n'.format(cl,','.join(map(str,fams))))

    with open(out_mods,'w') as outf, open(list_file,'r') as inf:
        header=inf.readline()
        outf.write(header+'\tClan\n')
        for line in inf:
            line = line.strip()
            splitline = line.split('\t')
            # mod = tuple(splitline[-2].split(','))
            family = int(splitline[-1])
            clan = fam_dict[family]
            outf.write('{}\t{}\n'.format(line,clan))

    #write file listing all clan/clusters by clan with their modules
    avg_cln_size = []
    avg_mods_in_cln = []
    with open(out_clusts,'w') as outf_c:
        for i in sorted(cluster_dict.keys()):
            matching_fams = cluster_dict[i]
            l_fam_matches = len(matching_fams)
            avg_cln_size.append(l_fam_matches)
            #list of [[mod_info,mod_tup,fam]]
            matches = [module_info[mod]+[mod,fam] for fam in \
                matching_fams for mod in family_modules[fam]]
            l_matches = len(matches)
            avg_mods_in_cln.append(l_matches)
            counts = Counter([dom for m in matches for dom in \
                m[-2]])
            outf_c.write(\
                '#Subcluster-clan {}, {} families, {} subclusters\n'.format(\
                i, l_fam_matches,l_matches))
            outf_c.write('#Occurrences: {}\n'.format(', '.join(\
                [dom+':'+str(c) for dom,c in counts.most_common()])))
            outf_c.write('#Features: {}\n'.format(', '.join(\
                ['{}:{:.2f}'.format(dom,c/l_matches) for dom,c in \
                counts.most_common()])))
            #maybe as a score the avg distance?
            for match in matches:
                outf_c.write('{}\t{}\t{}\n'.format('\t'.join(match[:-2]),\
                    ','.join(match[-2]),match[-1]))
    print('\nAverage clansize:', np.mean(avg_cln_size))
    print('Average amount of modules per clan:',np.mean(avg_mods_in_cln))


if __name__ == '__main__':
    start = time.time()
    cmd = get_commands()

    #make three outfiles: modified list_file with an extra column with clans,
    #fasta_file with >clan to families, file with #clan to modules
    fam_dict,feat_dict,fam_modules,mod_info = read_families(cmd.infile)

    out_sparse_matrix = cmd.infile.split('.txt')[0]+str(cmd.cutoff)+'.npz'
    if not os.path.isfile(out_sparse_matrix):
        print('\nCalculating distance matrix')
        dists = calc_all_dists(feat_dict, cmd.cutoff, cmd.cores)
        dist_matrix = make_dist_matrix(dists,len(feat_dict),square=True)
        print('  saving distance_matrix to',out_sparse_matrix)
        sp.save_npz(out_sparse_matrix, dist_matrix)
    else:
        print('\nLoaded distance matrix from',out_sparse_matrix)
        dist_matrix = sp.load_npz(out_sparse_matrix)
    print('  {} pairs have a distance below cutoff'.format(\
        len(dist_matrix.data)))

    rownms = sorted(fam_dict.keys())
    prefx = cmd.list_file.split('.txt')[0]
    run_dbscan(dist_matrix, cmd.dbscan_cutoff, rownms, cmd.list_file, prefx,\
        cmd.cores, fam_modules, mod_info)

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
