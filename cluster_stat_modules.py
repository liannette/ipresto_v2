#!/usr/bin/env python3
'''
Perform clustering of the statistical method modules.
Author: Joris Louwen
'''
import os
#make sure that numpy only uses one thread
os.environ['OMP_NUM_THREADS'] = '1'
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import numpy as np
# import os
import scipy.cluster.hierarchy as sch
import scipy.sparse as sp
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import TruncatedSVD
from sklearn.feature_extraction.text import CountVectorizer
import sklearn
from sys import argv
import time


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

    sparse_m: scipy_sparse_matrix
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
    out_mods, out_clusts, cores=1):
    '''Kmeans clustering on sparse_m with num_clusters and writes to file

    sparse_m: csr_matrix, shape(n_samples, n_features)
    num_clusters: int, number of clusters
    out_mods, out_clusts: str, filepaths
    cores: int, amount of cores to use
    '''
    print('\nRunning k-means')
    kmeans = KMeans(n_clusters=num_clusters, n_init=20, max_iter=1000, \
        random_state=595, verbose=0, tol=0.000001,n_jobs=cores).fit(sparse_m)
    print(kmeans)
    clust_centers = sp.csr_matrix(kmeans.cluster_centers_)
    labels = kmeans.labels_
    cluster_dict = defaultdict(list)
    np.set_printoptions(precision=2)
    print('Within-cluster sum-of-squares (inertia):', kmeans.inertia_)
    with open(out_mods,'w') as outf:
        outf.write(header+'\tFamily\n')
        for subcl,cl in zip(rownames,labels):
            cluster_dict[cl].append(subcl)
            outf.write('{}\t{}\t{}\n'.format(subcl,'\t'.join(modules[subcl]),\
                cl))
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

def run_dbscan(sparse_m, dist_cutoff, cores):
    '''
    '''
    print('\nCalculating distance matrix')
    dist_m = calc_jacc_index_matrix(sparse_m)
    print('\nRunning DBSCAN')
    clustering = DBSCAN(eps=dist_cutoff, metric='precomputed',min_samples=5,\
        n_jobs=cores)
    labels = clustering.labels_

if __name__ == '__main__':
    print('Start')
    start = time.time()
    mod_info = argv[1]
    number_clusters = int(argv[2])
    num_cores = int(argv[3])

    #outfiles
    prefix = mod_info.split('.txt')[0]+'_kmeans_'.format(number_clusters)
    out_mods = prefix +'_clustering.txt'
    # out_clust_centers = prefix + '_cluster_centers.txt'
    out_clusts = prefix + '_clusters_with modules.txt'
    #construct feature matrix
    modules = {} #keep track of info
    rownames = [] #in case mod_nums are not sequential
    corpus = [] #list of strings
    vectorizer = CountVectorizer(lowercase=False,binary=True,dtype=np.int32,\
        token_pattern=r"(?u)[^,]+") #finds everything separated by ','
    with open(mod_info,'r') as inf:
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

    cluster_kmeans(sparse_feat_matrix, modules, number_clusters, rownames,\
        colnames, out_mods, out_clusts, cores=num_cores)

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
