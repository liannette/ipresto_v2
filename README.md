# iPRESTO

iPRESTO (integrated Prediction and Rigorous Exploration of biosynthetic
Sub-clusters Tool)
is a command line tool for the detection of gene sub-clusters in
a set of Biosynthetic Gene Clusters (BGCs) in GenBank format. BGCs are tokenised
by representing each gene as a combination of its Pfam domains, where subPfams
are used to increase resolution. Tokenised BGCs are filtered for redundancy
using similarity network with an Adjacency Index of domains as a distance metric.
For the detection of sub-clusters two methods are used: PRESTO-STAT, which is
based on the statistical algorithm from Del Carratore et al. (2019), and the
novel method PRESTO-TOP, which uses topic modelling with Latent Dirichlet
Allocation. The sub-clusters found with iPRESTO can then be linked to Natural
Product substructures.

Developed by:

Joris J.R. Louwen1, Satria A. Kautsar1, Sven van der Burg2, Marnix H. Medema1*, Justin J.J. van der Hooft1,3*
1. Bioinformatics Group, Wageningen University, Wageningen, the Netherlands
2. Netherlands eScience Center, Amsterdam, the Netherlands
3. Department of Biochemistry, University of Johannesburg, Johannesburg, South Africa

*Corresponding authors

E-mail: marnix.medema@wur.nl, justin.vanderhooft@wur.nl 

![Workflow](final_workflow_black_900ppi.png)

# Citing iPRESTO

Did you use iPRESTO for your work? Please cite our paper: https://doi.org/10.1371/journal.pcbi.1010462
In there, you can also find more background information on the tool.

## Dependencies

iPRESTO is build and tested in python3.6. The required python packages are
automatically installed when using pip or setup.py. We recommend installing
iPRESTO in a conda environment like so:
```
# create new environment
conda create -n ipresto python=3.6

# activate new environment
conda activate ipresto

# install ipresto and dependencies
python -m pip install iPRESTO
```

iPRESTO also requires the HMMER suite. If HMMER is not installed on your
system, it can be installed with conda in your ipresto environment:
```
conda install -c bioconda hmmer
```

Then you can download this repo to be able to run iPRESTO.
```
# navigate to where you want to download the ipresto repo
git clone https://git.wur.nl/bioinformatics/iPRESTO.git
```

## Querying existing sub-cluster models

The sub-cluster models that were created in the publication can be downloaded from
Zenodo at https://doi.org/10.5281/zenodo.6953656. They can then be used to query your own BGCs
for sub-clusters. At that link you can also find the HMMs used in the publication (with the subPfam HMMs)
and an example clusterfile with tokenised BGCs from the antiSMASH-DB dataset.

## Usage

ipresto.py executes the main functionalities of iPRESTO. See below for
example commands. Toggle -h or --help for additional command line arguments
and default values (also available for many other additional scripts). Generally,
the input for iPRESTO analysis is a directory with BGCs in GenBank format, and 
a hmmpressed pHMM database.

ipresto.py performs pre-processing, the PRESTO-STAT and PRESTO-TOP methods, and 
visualisation. It takes a directory of BGCs in GenBank format and a hmmpressed
pHMM database as input. Redundancy filtering is
on by default but can be turned of by toggling --no_redundancy_filtering.
To query existing sub-cluster models based on the publication:
```
#ipresto with GBK folder input querying existing models
python ipresto.py -i my_gbk_dir -o output_dir --hmm_path Pfam_100subs_tc.hmm -c 12
        --include_list files/biosynthetic_domains.txt  # only include biosyntetic domains (see publication)
        --stat_subclusters PRESTO-STAT_subclusters.txt
        --top_motifs_model lda_model
        --no_redundancy_filtering  # probably you want to query all BGCs
```

It is also possible to start from a
clusterfile.csv with the flag --start_from_clusterfile to not have to read the domtables everytime.
```
#ipresto with clusterfile input querying existing models
# -i and --hmm_path have to be supplied symbolically
python ipresto.py --start_from_clusterfile my_clusterfile.csv -o output_dir -c 12
        --include_list files/biosynthetic_domains.txt  # only include biosyntetic domains (see publication)
        --stat_subclusters PRESTO-STAT_subclusters.txt
        --top_motifs_model lda_model
        --no_redundancy_filtering  # probably you want to query all BGCs
        -i bla --hmm_path bla  # to prevent ipresto from crashing (will get updated in the future)
```

Omitting --stat_subclusters or --top_motifs_model will result in PRESTO-STAT and PRESTO-TOP
detecting new sub-cluster models in your set of input BGCs. When building new models you
probably want to filter out redundant BGCs (so omit the --no_redundancy_filtering flag).
```
#ipresto with GBK folder input creating new sub-cluster models
python ipresto.py -i my_gbk_dir -o output_dir --hmm_path Pfam_100subs_tc.hmm -c 12
        --include_list files/biosynthetic_domains.txt  # only include biosyntetic domains (see publication)
        -p 0.2  # you can override the defualt p-value for PRESTO-STAT
        -t 500  # you can override the defualt number of topics for PRESTO-TOP
```

Visualisations of the sub-cluster output are made by ipresto.py, but visualisations can
also be made separately using subcluster_arrower.py.
One can provide one or more BGCs in GenBank format.
```
#one BGC
python3 ipresto/subcluster_arrower.py --one -f BGC0000052.gbk
        -c files/domains_colour_file.tsv -d preprocessing_domhits_file.txt
        -o BGC0000052.html -s input_presto_stat_subclusters.txt
        -l bgc_topics_filtered.txt --include_list files/biosynthetic_domains.txt
#multiple BGCs
python3 ipresto/subcluster_arrower.py -f file_with_gbk_locations.txt
        -c files/domains_colour_file.tsv -d preprocessing_domhits_file.txt
        -o BGC0000052.html -s input_presto_stat_subclusters.txt
        -l bgc_topics_filtered.txt --include_list files/biosynthetic_domains.txt
```

See below for an example clusterfile. Genes (and BGC names) are separated by
commas, domains in the same gene by semi-colons and genes without domains are
represented by a dash.
```
BGC_name1,Lactamase_B,adh_short,ketoacyl-synt;Ketoacyl-synt_C,-\n
BGC_name2,-,Lant_dehydr_N;Lant_dehydr_C,LANC_like\n
```

Other scripts fullfill additional roles for more functionality. subPfams can be
created with https://github.com/satriaphd/build_subpfam. subPfams used in the
publication can be downloaded from the Zenodo https://doi.org/10.5281/zenodo.6953656.

## Output

Most important outputs are *_presto_stat_subclusters.txt, presto_top/bgc_topics_filtered.txt and
*_ipresto_output_visualisation.html. They contain the input BGCs queried to PRESTO-STAT and the
PRESTO-TOP sub-clusers, and visualisation of the sub-cluster results, respectively. The
annotations for the sub-cluster motifs from PRESTO-TOP can also be found on Zenodo.
