#!/usr/bin/env python3

import argparse
from multiprocessing import cpu_count

from ipresto.pipeline import IprestoPipeline
import sys


def get_commands():
    parser = argparse.ArgumentParser(
        description="iPRESTO uses topic modelling and statistical analyses "
        "to detect sub-clusters of co-evolving genes in Gene Clusters, which "
        "can be linked to substructures of Natural Products. This script is "
        "the main functionality of iPRESTO. It can build new sub-cluster "
        "models from gbks or use previously constructed models to detect "
        "sub-clusters in unseen gbks.",
    )
    # Output directory
    parser.add_argument(
        "--out",
        dest="out_dir_path",
        required=True,
        metavar="<dir>",
        help="Output directory, this will contain all output data files.",
    )

    # Input options
    parser.add_argument(
        "--gbks",
        dest="gbk_dir_path",
        metavar="<dir>",
        help="Input directory containing gbk files of the gene clusters.",
    )
    parser.add_argument(
        "--clusters",
        dest="existing_clusterfile",
        default=None,
        metavar="<file>",
        help="Path to a previously created CSV file of tokenized clusters. "
        "This option overrides the --gbks argument, although --gbks must "
        "still be provided.",
    )
    # Preprocessing settings
    parser.add_argument(
        "--hmm",
        dest="hmm_file_path",
        metavar="<file>",
        help="Path to the HMM file containing protein domain HMMs that has been "
        "processed with hmmpress.",
    )
    parser.add_argument(
        "--exclude_name",
        dest="exclude_name",
        default=["final"],
        nargs="+",
        metavar="<str>",
        help="If any string in this list occurs in the gbk filename, this "
        "file will not be used for the analysis (default: ['final']).",
    )
    parser.add_argument(
        "--incl_contig_edge",
        dest="include_contig_edge_clusters",
        action="store_true",
        default=False,
        help="Include clusters that lie on a contig edge. (default = false)",
    )
    parser.add_argument(
        "--max_domain_overlap",
        dest="max_domain_overlap",
        default=0.1,
        metavar="<float>",
        help="Specify at which overlap percentage domains are considered to overlap. "
        "Domain with the best score is kept (default=0.1).",
    )
    parser.add_argument(
        "--min_genes_per_bgc",
        dest="min_genes_per_bgc",
        default=2,
        type=int,
        metavar="<int>",
        help="Minimum number of non-empty genes required in a cluster to "
        "be included in the analysis. A gene is considered empty if it "
        "lacks any protein domains (default: 2).",
    )
    parser.add_argument(
        "--domain_filtering",
        action="store_true",
        default=False,
        help="If provided, domain filtering will be performed, to remove "
        "any non-biosynthetic protein domains from the tokenized clusters.",
    )
    parser.add_argument(
        "--similarity_filtering",
        action="store_true",
        default=False,
        help="If provided, similarity filtering of the tokenized clusters "
        "will be performed",
    )
    parser.add_argument(
        "--similarity_cutoff",
        dest="similarity_cutoff",
        default=0.95,
        type=float,
        metavar="<float>",
        help="Cutoff for cluster similarity in redundancy filtering (default:0.95).",
    )
    # PRESTO-STAT
    parser.add_argument(
        "--stat_modules",
        dest="stat_modules_file_path",
        default=None,
        metavar="<file>",
        help="Text file containing previously inferred subclusters to detect "
        "in the input. If not provided, PRESTO-STAT will run to detect new "
        "subclusters in the input (default: None)",
    )
    parser.add_argument(
        "--stat_min_gene_occurrence",
        default=3,
        type=int,
        metavar="<int>",
        help="Remove tokenized genes that occur fewer than the specified "
        "number of times in the data (default: 3).",
    )
    parser.add_argument(
        "-p",
        "--stat_pval_cutoff",
        dest="stat_pval_cutoff",
        default=0.01,
        type=float,
        help="P-value cutoff for determining a significant interaction in module detection "
        "(default: 0.01)",
        metavar="<float>",
    )
    # # PRESTO-TOP
    # parser.add_argument(
    #     '--lda_model',
    #     dest="top_lda_model",
    #     default=None,
    #     metavar="<file>",
    #     help='Use PRESTO-TOP with existing sub-cluster motifs in an LDA model. '
    #          'Supply here the path to the model. In that location there should be also '
    #          'model.dict, model.expElogbeta.npy, model.id2word, model.state, '
    #          'model.state.sstats.npy',
    # )
    # parser.add_argument(
    #     "-t", "--topics",
    #     dest="n_topics",
    #     help="Amount of topics to use for the LDA model in PRESTO-TOP (default: 1000)",
    #     default=1000,
    #     type=int,
    #     metavar="<int>"
    # )

    # parser.add_argument(
    #     "-f", "--min_feat_score",
    #     dest="min_feat_score",
    #     help="Only include features until their scores add up to this number (default: 0.95) Can "
    #          "be combined with feat_num, where feat_num features are selected or features that "
    #          "add up to min_feat_score",
    #     type=float,
    #     default=0.95,
    #     metavar="<float>"
    # )
    # parser.add_argument(
    #     "-n", "--feat_num",
    #     dest="feat_num",
    #     help="Include the first feat_num features for each topic (default: 75)",
    #     type=int,
    #     default=75,
    #     metavar="<int>"
    # )
    # parser.add_argument(
    #     "-a", "--amplify",
    #     dest="amplify",
    #     help="Amplify the dataset in order to achieve a better LDA model. Each BGC will be present "
    #          "amplify times in the dataset. After calculating the LDA model the dataset will be "
    #          "scaled back to normal.",
    #     type=int,
    #     default=None,
    #     metavar="<int>"
    # )
    # parser.add_argument(
    #     "--visualise",
    #     help="Make a visualation of the LDA model with pyLDAvis (html file). If number of topics "
    #          "is too big this might fail. No visualisation will then be made",
    #     default=False,
    #     action="store_true"
    # )
    # parser.add_argument(
    #     "--classes",
    #     help="A file containing classes of the BGCs used in the analysis. First column should "
    #          "contain matching BGC names. Consecutive columns should contain classes.",
    #     default=False,
    #     metavar="<file>"
    # )
    # parser.add_argument(
    #     "--plot",
    #     help="If provided: make plots about several aspects of the presto-top output",
    #     default=False,
    #     action="store_true"
    # )
    # parser.add_argument(
    #     "--known_subclusters",
    #     help="A tab delimited file with known subclusters. Should contain subclusters in the last "
    #          "column and BGC identifiers in the first column. Subclusters are comma separated "
    #          "genes represented as domains. Multiple domains in a gene are separated by "
    #          "semi-colon.",
    #     metavar="<file>"
    # )
    # parser.add_argument(
    #     "-I", "--iterations",
    #     help="Amount of iterations for training the LDA model (default: 1000)",
    #     default=1000,
    #     type=int,
    #     metavar="<int>"
    # )
    # parser.add_argument(
    #     "-C", "--chunksize",
    #     default=2000,
    #     type=int,
    #     help='The chunksize used to train the model (default: 2000)',
    #     metavar="<int>"
    # )
    # parser.add_argument(
    #     "-u", "--update",
    #     help="If provided and a model already exists, the existing model will be updated with "
    #          "original parameters, new parameters cannot be passed in the LdaMulticore version.",
    #     default=False,
    #     action="store_true"
    # )
    # parser.add_argument(
    #     "--alpha",
    #     default="symmetric",
    #     help="alpha parameter for the LDA model, see gensim. Options: (a)symmetric, auto, or <int>"
    # )
    # parser.add_argument(
    #     "--beta",
    #     default="symmetric",
    #     help="beta parameter for the LDA model, see gensim. Options: (a)symmetric, auto, or <int>"
    # )
    # # Visualisation
    # parser.add_argument(
    #     "--visualise_subclusters",
    #     default=False,
    #     help="If provided, subclusters will be visualised for all gbk inputs, otherwise just the "
    #          "1000 first bgcs of the data will be visualised to consider time/space",
    #     action="store_true"
    # )
    # Other arguments
    parser.add_argument(
        "-c",
        "--cores",
        dest="cores",
        default=cpu_count(),
        help="Set the number of cores the script may use (default: use all "
        "available cores)",
        type=int,
        metavar="<int>",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Prints more detailed information.",
    )

    args = parser.parse_args()

    # Validation logic
    if not args.existing_clusterfile:
        if not args.gbk_dir_path or not args.hmm_file_path:
            parser.error("Either --gbks and --hmm must be provided, or --clusters must be specified.")

    return args


def main():
    """
    Main function to execute the iPRESTO pipeline.

    This function retrieves command line arguments, prints them if verbose mode is enabled,
    and then runs the main pipeline with the provided arguments.

    Parameters:
    None

    Returns:
    None
    """
    # Get the command line arguments
    cmd = get_commands()

    # Print the command line arguments if verbose
    if cmd.verbose:
        print("Command:", " ".join(sys.argv))
        print(cmd)

    # Execute the main pipeline with provided arguments
    IprestoPipeline().run(
        cmd.out_dir_path,
        cmd.gbk_dir_path,
        cmd.existing_clusterfile,
        cmd.exclude_name,
        cmd.include_contig_edge_clusters,
        cmd.hmm_file_path,
        cmd.max_domain_overlap,
        cmd.min_genes_per_bgc,
        cmd.domain_filtering,
        cmd.similarity_filtering,
        cmd.similarity_cutoff,
        cmd.stat_modules_file_path,
        cmd.stat_min_gene_occurrence,
        cmd.stat_pval_cutoff,
        cmd.cores,
        cmd.verbose,
    )


if __name__ == "__main__":
    main()
