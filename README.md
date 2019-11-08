# iPRESTO

iPRESTO (integrated Prediction and Rigorous Exploration of biosynthetic
Sub-clusters Tool)
is a collection of python scripts for the detection of gene sub-clusters in
a set of Biosynthetic Gene Clusters (BGCs) in GenBank format. BGCs are tokenised
by representing each gene as a combination of its Pfam domains, where subPfams
are used to increase resolution. Tokenised BGCs are filtered for redundancy
using similarity network with an Adjacency Index of domains as a distance metric.
For the detection of sub-clusters two methods are used: PRESTO-STAT, which is
based on the statistical algorithm from Del Carratore et al. (2019), and the
novel method PRESTO-TOP, which uses topic modelling with Latent Dirichlet
Allocation. The sub-clusters found with iPRESTO can then be linked to Natural
Product substructures.

Developed by Joris Louwen.
Supervisors: Marnix Medema (PI), Justin van der Hooft and Satria Kautsar.
All from the Bioinformatics group at Wageningen University. 

![Workflow](final_workflow_black_900ppi.png)

## Usage

To use iPRESTO, there are some main scripts to use, which are explained with
example commands below. All main scripts have a -h or --help option for
additional command line arguments and default values. Generally, the input for
iPRESTO analysis is a directory with BGCs in GenBank format, and a hmmpressed
pHMM database.

preprocessing.py turns the input directory into a csv file with
tokenised BGCs (called clusterfile.csv) and filters out redundant BGCs.
```
python3 preprocessing.py -i my_gbk_dir -o output_dir --hmm_path Pfam_A.hmm
        --exclude final -c 12 -e True
```

presto_stat.py performs the PRESTO-STAT method. It can start with the same
input as preprocessing.py, but it is also possible to start from a
clusterfile.csv with the flag --start_from_clusterfile. Redundancy filtering is
on by default but can be turned of by toggling --no_redundancy_filtering.
```
#presto-stat with GBK folder input
python3 presto_stat.py -i my_gbk_dir -o output_dir --hmm_path Pfam_A.hmm
        --exclude final -c 12 -e True -p 0.1

#presto_stat with clusterfile input
# -i -o and --hmm_path have to be supplied symbolically
python3 presto_stat.py--start_from_clusterfile my_clusterfile.csv -c 12
        --no_redundancy_filtering -i symbolic -o symbolic --hmm_path symbolic
```

query_statistical_modules.py allows for querying a list of statistical
sub-clusters as produced by presto_stat.py. Input should be a clusterfile.
```
python3 query_statistical_modules.py -i my_clusterfile.csv -m my_modules.txt
        -c 12 -o my_clusterfile
```

presto_top.py performs the PRESTO-TOP method. It takes a clusterfile as input
and has many commandline options to modify its behaviour in for example the
construction of the LDA model. With the -r one can query an existing LDA model.
```
#Creating an LDA model and querying it at the same time.
python3 presto_top.py -i my_clusterfile.csv -o my_output_folder -c 10 -t 1000 -C 3000
        -I 2000 --min_genes 2 -f 0.95 -n 75 --classes my_bgc_classes.txt
        --known_subclusters known_subcl.txt
#Querying an existing model with -r
python3 presto_top.py -i my_clusterfile.csv -o my_output_folder -c 10 -t 1000
        --min_genes 2 -f 0.95 -n 75 --classes my_bgc_classes.txt
        --known_subclusters known_subcl.txt -r my_lda_model_location
```

subcluster_arrower.py creates powerful visualisations of the sub-cluster output.
One can provide one or more BGCs in GenBank format.
```
#one BGC
python3 subcluster_arrower.py
#multiple BGCs
python3 subcluster_arrower.py
```

An example clusterfile:
```
BGC_name1,Lactamase_B,adh_short,ketoacyl-synt;Ketoacyl-synt_C,-\n
BGC_name2,-,Lant_dehydr_N;Lant_dehydr_C,LANC_like\n
```

Other scripts fullfill additional roles for more functionality. subPfams can be
created with https://github.com/satriaphd/build_subpfam.

## Dependencies

iPRESTO is build in python3.6. It requires hmmscan and some python packages:
* Biopython
* gensim
* networkx
* numpy
* pandas
* statsmodels
