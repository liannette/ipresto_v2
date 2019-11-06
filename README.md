# iPRESTO

iPRESTO (integrated Predictions for the Rapid Elucidation of Sub-clusters Tool)
is a collection of python scripts for the detection of gene sub-clusters in
a set of Biosynthetic Gene Clusters (BGCs) in GenBank format. BGCs are tokenised
by representing each gene as a combination of its Pfam domains, where subPfams
are used to increase resolution. Tokenised BGCs are filtered for redundancy
using an Adjacency Index of domains. For the detection of sub-clusters the
PRESTO-STAT method is used, which is based on the statistical algorithm from
Del Carratore et al., (2019). Additionally we use the novel method PRESTO-TOP
for sub-cluster detection, which uses topic modelling with Latent Dirichlet
Allocation. The sub-clusters found with iPRESTO can then be linked to Natural
Product substructures.

Written by Joris Louwen.
Supervisors: Marnix Medema (PI), Justin van der Hooft and Satria Kautsar.
All from the Bioinformatics group at Wageningen University. 
