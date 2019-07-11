#!/usr/bin/env python3
'''
Perform clustering of the statistical method modules.
Author: Joris Louwen
'''
import os
#make sure that numpy only uses one thread
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import numpy as np
# import os
import scipy.cluster.hierarchy as sch
import scipy.sparse as sp
from scipy.spatial.distance import cdist
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import TruncatedSVD
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics import silhouette_score, pairwise_distances_chunked
from sys import argv
import time

def get_commands():
    parser = argparse.ArgumentParser(description="A script to cluster the\
        statistical modules into families with different algorithms.")
    parser.add_argument("-i", "--infile", help="Input tsv\
        file with module information (filtered_modules.txt). A header is\
        present, one module per line where a module is the last element in\
        the line. Genes separated by , and domains by ;", required=True)
    parser.add_argument('-c', '--cores', help='Cores to use, default = 1',
        default=1, type=int)
    parser.add_argument('-k', '--k_clusters',help='Amount of clusters to use\
        with k-means clustering. Only used if k-means is specified, default =\
        1000',
        default=1000, type=int)
    parser.add_argument('-m', '--method', help='Method for clustering. Should\
        be a number:\n\t0 = all methods\n\t1 = k-means\n\t2 = DBSCAN', \
        default = 0, type=int)
    return parser.parse_args()


def plot_svd_components(sparse_m):
    '''Plots first two components of truncatedSVD analysis (PCA)

    sparse_m: scipy_sparse_matrix
    '''
    svd = TruncatedSVD(n_components=2, n_iter=7, random_state=595)
    components = svd.fit_transform(sparse_m)
    # print(components)
    print(svd.explained_variance_ratio_)
    x,y=zip(*components)
    plt.scatter(x,y)
    plt.show()

def calc_jacc_distance_matrix(sparse_m):
    '''
    Returns distance matrix (jaccard distance) of subclusters in sparse_m

    sparse_m: scipy_sparse_matrix
    distance_m: np.array of condensed dist matrix, each value between
        0 (min distance) and 1 (max distance) (1d array that represents upper
        triangle)
    '''
    #this will construct the upper triangle of the matrix
    len_rows = sparse_m.shape[0]
    distance_m = np.empty(shape=(0,1),dtype=np.float32)
    for i in range(len_rows-1):
        for j in range(i+1, len_rows):
            doms_a = set(sparse_m[i].nonzero()[1])
            doms_b = set(sparse_m[j].nonzero()[1])
            jd = 1 - (len(doms_a & doms_b) / len(doms_a | doms_b))
            distance_m = np.append(distance_m, [jd])
    return distance_m

def calc_jacc_index_matrix(sparse_m):
    '''
    Returns sparse distance matrix (jaccard index) of subclusters in sparse_m

    sparse_m: csr_matrix shape(n_samples,n_features)
    distance_m: csr_matrix, shape(n_samples,n_samples)
    '''
    #this will construct a square matrix
    len_rows = sparse_m.shape[0]
    row = []
    col = []
    data = []
    for i in range(len_rows-1):
        for j in range(i+1, len_rows):
            doms_a = set(sparse_m[i].nonzero()[1])
            doms_b = set(sparse_m[j].nonzero()[1])
            ji = len(doms_a & doms_b) / len(doms_a | doms_b)
            if ji > 0:
                jd = 1-ji
                row.append(i)
                row.append(j) #make it a square matrix
                col.append(j)
                col.append(i)
                data.append(jd)
                data.append(jd)
                #assume there are no identical subclusters 1-ji would become 0
    distance_m = sp.csr_matrix((data, (row,col)), shape=(len_rows,len_rows))
    return distance_m

def new_euclidean_distances(X, Y=None, Y_norm_squared=None, squared=False):
    '''

    X, Y: sparse matrice or arrays with shape (1,n)
    '''
    if sp.issparse(X):
        x = X.A
    else:
        x = X
    if sp.issparse(Y):
        y = Y.A
    else:
        if Y != None:
            y = Y
        else:
            y = x
    print(x)
    print(y)
    #calc Soergel distance, variant on Jaccard distance, but with weights
    #so that centroids can be floats
    intersection = np.array(0,dtype=np.float64)
    union = np.array(0,dtype=np.float64)
    for tup in zip(x[0],y[0]):
        intersection += min(tup)
        union += max(tup)
    jd = np.array([[1-intersection/union]],dtype=np.float64)
    return jd

def cluster_hierarchical(data):
    '''
    data: sparse matrix
    '''
    dist = calc_jacc_distance_matrix(data)
    print(dist)
    clust = sch.linkage(dist)
    print(clust)

def cluster_kmeans(sparse_m, modules, num_clusters, rownames, colnames, \
    prefix, cores=1):
    '''Kmeans clustering on sparse_m with num_clusters and writes to file

    sparse_m: csr_matrix, shape(n_samples, n_features)
    modules: dict {mod_num:[info,modules]}
    num_clusters: int, number of clusters
    rownames: list of ints, [mod_nums], sequential mod_nums, keeping track of
        rows of sparse_m
    colnames: list of str, all domains sequential order to keep track of
        columns of sparse_m
    prefix: str, prefix of outfile
    cores: int, amount of cores to use
    '''
    print('\nRunning k-means')
    #outfiles
    kmeans_pre = '_kmeans_{}_families'.format(num_clusters)
    out_mods = prefix + kmeans_pre + '.txt'
    # out_clust_centers = prefix + '_cluster_centers.txt'
    out_clusts = prefix + kmeans_pre + '_by_family.txt'

    #running algorithm
    kmeans = KMeans(n_clusters=num_clusters, n_init=20, max_iter=1000, \
        random_state=595, verbose=0, tol=0.000001,n_jobs=cores).fit(sparse_m)
    print(kmeans)
    clust_centers = sp.csr_matrix(kmeans.cluster_centers_)
    labels = kmeans.labels_
    cluster_dict = defaultdict(list)
    np.set_printoptions(precision=2)
    print('Within-cluster sum-of-squares (inertia):', kmeans.inertia_)
    #link each module to a family/k-means cluster
    with open(out_mods,'w') as outf:
        outf.write(header+'\tFamily\n')
        for subcl,cl in zip(rownames,labels):
            cluster_dict[cl].append(subcl)
            outf.write('{}\t{}\t{}\n'.format(subcl,'\t'.join(modules[subcl]),\
                cl))
    #write file of listing all families/clusters by family with their modules
    with open(out_clusts,'w') as outf_c:
        for i in range(clust_centers.shape[0]):
            matches = cluster_dict[i]
            counts = Counter([dom for m in matches for dom in \
                modules[m][-1].split(',')])
            spars = clust_centers[i]
            feat_inds = spars.nonzero()[1]
            feat_tups = [(spars[0,ind],colnames[ind]) for ind in \
                feat_inds]
            feat_format = ['{}:{:.2f}'.format(dom,float(score)) for score,dom\
                in sorted(feat_tups,reverse=True)]
            outf_c.write('#Subcluster-family {}, {} subclusters\n'.format(i,\
                len(matches)))
            outf_c.write('#Occurrences: {}\n'.format(', '.join(\
                [dom+':'+str(c) for dom,c in counts.most_common()])))
            outf_c.write('#Features: {}\n'.format(', '.join(feat_format)))
            #maybe as a score the distance to the cluster center?
            for match in matches:
                outf_c.write('{}\t{}\n'.format(match,\
                    '\t'.join(modules[match])))

def jaccard_on_sparse(row1, row2):
    '''Use jaccard distance on sparse matrix rows, returns float in array
    '''
    dist = cdist(row1.A,row2.A,metric='jaccard')
    return dist


def run_dbscan(sparse_m, dist_cutoff, modules, rownames, prefix, cores):
    '''
    Nearest neighbour algorithm with jaccard distance

    sparse_m: csr_matrix, shape(n_samples, n_features)
    dist_cutoff: float, between 0 and 1 that denotes when modules are
        in a cluster based on jaccard distance
    modules: dict {mod_num:[info,modules]}
    rownames: list of ints, [mod_nums], sequential mod_nums, keeping track of
        rows of sparse_m
    prefix: str, prefix of outfile
    cores: int, amount of cores to use
    '''
    # print('\nCalculating distance matrix')
    # dist_m = calc_jacc_index_matrix(sparse_m)
    # print(dist_m)
    print('\nRunning DBSCAN')
    clustering = DBSCAN(eps=dist_cutoff, metric=jaccard_on_sparse,\
        min_samples=5, n_jobs=cores).fit(sparse_m)
    print(clustering)
    labels = clustering.labels_
    n_labels = len(set(labels))
    n_noise = list(labels).count(-1)
    if n_noise:
        n_labels -= 1
    print('Found {} clusters and {} noise points'.format(n_labels,n_noise))

    # #measure for average dist in cluster vs average dist to nearest clust
    # #1 is best -1 is worst
    # sil = silhouette_score(dist_m, labels, metric='precomputed')
    # print("  silhouette Coefficient: {0.3f}".format(sil))

    #outfiles
    dbscan_pre = '_dbscan_{}_families'.format(n_labels)
    out_mods = prefix + dbscan_pre + '.txt'
    out_clusts = prefix + dbscan_pre + '_by_family.txt'

    #link each module to a family/cluster
    cluster_dict = defaultdict(list)
    with open(out_mods,'w') as outf:
        outf.write(header+'\tFamily\n')
        for subcl,cl in zip(rownames,labels):
            cluster_dict[cl].append(subcl)
            outf.write('{}\t{}\t{}\n'.format(subcl,'\t'.join(modules[subcl]),\
                cl))
    #write file of listing all families/clusters by family with their modules
    with open(out_clusts,'w') as outf_c:
        for i in range(len(labels)):
            matches = cluster_dict[i]
            counts = Counter([dom for m in matches for dom in \
                modules[m][-1].split(',')])
            outf_c.write('#Subcluster-family {}, {} subclusters\n'.format(i,\
                len(matches)))
            outf_c.write('#Occurrences: {}\n'.format(', '.join(\
                [dom+':'+str(c) for dom,c in counts.most_common()])))
            #maybe as a score the avg distance?
            for match in matches:
                outf_c.write('{}\t{}\n'.format(match,\
                    '\t'.join(modules[match])))

if __name__ == '__main__':
    print('Start')
    start = time.time()
    # mod_info = argv[1]
    # number_clusters = int(argv[2])
    # num_cores = int(argv[3])

    cmd = get_commands()

    #outfiles
    out_prefix = cmd.infile.split('.txt')[0]

    #construct feature matrix
    modules = {} #keep track of info
    rownames = [] #in case mod_nums are not sequential
    corpus = [] #list of strings
    vectorizer = CountVectorizer(lowercase=False,binary=True,dtype=np.int32,\
        token_pattern=r"(?u)[^,]+") #finds everything separated by ','
    with open(cmd.infile,'r') as inf:
        print('\nReading module file')
        #{mod_num:[info]}
        header = inf.readline().strip('\n') #header
        for line in inf:
            line = line.strip().split('\t')
            mod_num = int(line[0])
            modules[mod_num] = line[1:]
            rownames.append(mod_num)
            corpus.append(line[-1])
    print('\nBuilding sparse matrix representation of module features')
    sparse_feat_matrix = vectorizer.fit_transform(corpus)
    colnames = vectorizer.get_feature_names()
    print('  {} features'.format(len(colnames)))

    if cmd.method in [0,1]:
        cluster_kmeans(sparse_feat_matrix, modules, cmd.k_clusters, rownames,\
            colnames, out_prefix, cores=cmd.cores)
    if cmd.method in [0,2]:
        neighbour_cutoff = 0.5
        run_dbscan(sparse_feat_matrix, neighbour_cutoff, modules, rownames,\
            out_prefix, cmd.cores)

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
