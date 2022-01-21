#!/usr/bin/env python3
"""
Author: Joris Louwen (joris.louwen@wur.nl)

Part of iPRESTO, Bioinformatics group Wageningen University.
PIs: Marnix Medema, Justin van der Hooft
Collaborators: Satria Kautsar

usage:
python ipresto.py -h
"""
import argparse
from ipresto.presto_stat.presto_stat import *
from ipresto.presto_stat import query_statistical_modules as q_stat
from ipresto.presto_top.presto_top import *
from multiprocessing import cpu_count
from sys import argv
import logging
from typing import Union, List


def get_commands():
    parser = argparse.ArgumentParser(
        description="iPRESTO uses topic modelling and statistical analyses \
        to detect sub-clusters of co-evolving genes in Gene Clusters, which \
        can be linked to substructures of Natural Products. This script is \
        the main functionality of iPRESTO. It can build new sub-cluster \
        models from gbks or use previously constructed models to detect \
        sub-clusters in unseen gbks.")
    parser.add_argument(
        "-i", "--in_folder", dest="in_folder", help="Input directory of gbk \
        files", required=True, metavar="<dir>")
    parser.add_argument(
        "-o", "--out_folder", dest="out_folder", required=True,
        help="Output directory, this will contain all output data files.",
        metavar="<dir>")
    parser.add_argument(
        "--hmm_path", dest="hmm_path", required=True, metavar="<file>",
        help="File containing domain hmms that is hmmpress-processed.")
    parser.add_argument(
        "--stat_subclusters", default=None, metavar="<file>", help="Txt file \
        containing previously inferred subclusters to detect in the input - \
        if not provided, PRESTO-STAT will run to detect new subclusters in \
        the input (default: None)")
    parser.add_argument(
        '--top_motifs_model', help='Use PRESTO-TOP with existing \
        sub-cluster motifs in an LDA model. Supply here the path to the \
        model. In that location there should be also model.dict, \
        model.expElogbeta.npy, model.id2word, model.state, \
        model.state.sstats.npy', required=False, default=False,
        metavar="<file>")
    parser.add_argument(
        "--include_list", dest="include_list", default=None, help="If \
        provided only the domains in this file will be taken into account in \
        the analysis. One line should contain one Pfam ID (default: None - \
        meaning all Pfams from database)", metavar="<file>")
    parser.add_argument(
        "--start_from_clusterfile", default=None, help="A file with BGCs and \
        domain-combinations to start with (csv and domains in a gene \
        separated by ';'). This overwrites in_folder (which still has to be \
        supplied symbolically) and use_domtabs/use_fastas.",
        metavar="<file>")
    parser.add_argument(
        "-c", "--cores", dest="cores", default=cpu_count(),
        help="Set the number of cores the script may use (default: use all \
        available cores)", type=int, metavar="<int>")
    parser.add_argument(  # todo: make invalid if only querying models
        "--no_redundancy_filtering", default=False, help="If provided, \
            redundancy filtering will not be performed", action="store_true")
    parser.add_argument(
        "--visualise_subclusters", default=False, help="If provided, \
        subclusters will be visualised for all gbk inputs, otherwise just the \
        1000 first bgcs of the data will be visualised to consider time/space",
        action="store_true")
    parser.add_argument(
        "--exclude", dest="exclude", default=["final"], nargs="+",
        help="If any string in this list occurs in the gbk filename, this \
        file will not be used for the analysis. (default: [final])",
        metavar="<str>")
    parser.add_argument(
        "-v", "--verbose", dest="verbose", required=False, action="store_true",
        default=False, help="Prints more detailed information.")
    parser.add_argument(
        "-d", "--domain_overlap_cutoff", dest="domain_overlap_cutoff",
        default=0.1, help="Specify at which overlap percentage domains are \
        considered to overlap. Domain with the best score is kept \
        (default=0.1).", metavar="<float>")
    parser.add_argument(  # todo: again include query edge bgcs when querying
        "-e", "--exclude_contig_edge", dest="exclude_contig_edge",
        default=False, help="Exclude clusters that lie on a contig edge \
        (default = false)", action="store_true")
    parser.add_argument(
        "-m", "--min_genes", dest="min_genes", default=0, help="Provide the \
        minimum size of a BGC to be included in the analysis. Default is 0 \
        genes", type=int, metavar="<int>")
    parser.add_argument(
        "--min_doms", dest="min_doms", default=0, help="The minimum amount of \
        domains in a BGC to be included in the analysis. Default is 0 domains",
        type=int, metavar="<int>")
    parser.add_argument(
        "--sim_cutoff", dest="sim_cutoff", default=0.95, help="Cutoff for \
        cluster similarity in redundancy filtering (default:0.95)", type=float,
        metavar="<float>")
    parser.add_argument(
        "--remove_genes_below_count", default=3, type=int, help="Remove genes \
        (domain combinations) when they occur less than <int> times in the \
        data (default: 3)", metavar="<int>")
    parser.add_argument(
        "-p", "--pval_cutoff", dest="pval_cutoff", default=0.1, type=float,
        help="P-value cutoff for determining a significant interaction in \
        module detection (default: 0.1)", metavar="<float>")
    parser.add_argument(
        "--use_fastas", dest="use_fastas", default=None, help="Use already \
        created fasta files from some folder", metavar="<dir>")
    parser.add_argument(
        "--use_domtabs", dest="use_domtabs", default=None, help="Use already \
        created domtables from some folder", metavar="<dir>")
    parser.add_argument(
        "-t", "--topics", dest="topics", help="Amount of topics to use for \
        the LDA model in PRESTO-TOP (default: 1000)", default=1000, type=int,
        metavar="<int>")
    parser.add_argument(
        "-f", "--min_feat_score", dest="min_feat_score", help="Only include \
        features until their scores add up to this number (default: 0.95) Can \
        be combined with feat_num, where feat_num features are selected or \
        features that add up to min_feat_score", type=float, default=0.95,
        metavar="<float>")
    parser.add_argument(
        "-n", "--feat_num", dest="feat_num", help="Include the first feat_num \
        features for each topic (default: 75)", type=int, default=75,
        metavar="<int>")
    parser.add_argument(
        "-a", "--amplify", dest="amplify", help="Amplify the dataset in order \
        to achieve a better LDA model. Each BGC will be present amplify times \
        in the dataset. After calculating the LDA model the dataset will be \
        scaled back to normal.", type=int, default=None, metavar="<int>")
    parser.add_argument(
        "--visualise", help="Make a visualation of the LDA model with \
        pyLDAvis (html file). If number of topics is too big this might fail. \
        No visualisation will then be made", default=False,
        action="store_true")
    parser.add_argument(
        "--classes", help="A file containing classes of the BGCs used in the \
        analysis. First column should contain matching BGC names. Consecutive \
        columns should contain classes.", default=False, metavar="<file>")
    parser.add_argument(
        "--plot", help="If provided: make plots about several aspects of the \
        presto-top output", default=False, action="store_true")
    parser.add_argument(
        "--known_subclusters", help="A tab delimited file with known \
        subclusters. Should contain subclusters in the last column and BGC \
        identifiers in the first column. Subclusters are comma separated \
        genes represented as domains. Multiple domains in a gene are \
        separated by semi-colon.", metavar="<file>")
    parser.add_argument(
        "-I", "--iterations", help="Amount of iterations for training the \
        LDA model (default: 1000)", default=1000, type=int, metavar="<int>")
    parser.add_argument(
        "-C", "--chunksize", default=2000, type=int, help='The chunksize \
        used to train the model (default: 2000)', metavar="<int>")
    parser.add_argument(
        "-u", "--update", help="If provided and a model already exists, the \
        existing model will be updated with original parameters, new \
        parameters cannot be passed in the LdaMulticore version.",
        default=False, action="store_true")

    return parser.parse_args()


def preprocessing_bgcs_to_dom_combinations(
        out_folder: str,
        in_folder: str,
        hmm_path: str,
        start_from_clusterfile: Union[str, None],
        exclude: List[str],
        exclude_contig_edge: bool,
        min_genes: int,
        cores: int,
        verbose: bool,
        use_fastas: Union[str, None],
        use_domtabs: Union[str, None],
        domain_overlap_cutoff: float) -> str:
    """Processes BGCs (gbks) into list of domain (Pfams) combinations

    :param out_folder:
    :param in_folder:
    :param hmm_path:
    :param start_from_clusterfile:
    :param exclude:
    :param exclude_contig_edge:
    :param min_genes:
    :param cores:
    :param verbose:
    :param use_fastas:
    :param use_domtabs:
    :param domain_overlap_cutoff:
    :return: path to clusterfile - csv of domain combinations:
        bgc_name,dom1;dom2,dom1,dom3;dom4;dom1
    """
    if start_from_clusterfile:
        if not os.path.isdir(out_folder):
            f_command = 'mkdir {}'.format(out_folder)
            subprocess.check_call(f_command, shell=True)
        filepre = os.path.split(start_from_clusterfile)[-1].split(
            '.csv')[0]
        clus_file = os.path.join(out_folder, filepre + '_clusterfile.csv')
        c_command = 'cp {} {}'.format(start_from_clusterfile, clus_file)
        subprocess.check_call(c_command, shell=True)
    else:
        fasta_folder, exist_fastas = process_gbks(
            in_folder, out_folder, exclude,
            exclude_contig_edge, min_genes, cores, verbose,
            use_fastas)
        dom_folder, exist_doms = hmmscan_wrapper(
            fasta_folder, hmm_path, verbose, cores, exist_fastas,
            use_domtabs)
        clus_file = parse_dom_wrapper(dom_folder, out_folder,
                                      domain_overlap_cutoff, verbose,
                                      exist_doms)
    return clus_file


def filtering_cluster_representations(
        clus_file: str,
        out_folder: str,
        no_redundancy_filtering: bool,
        min_genes: int,
        cores: int,
        verbose: bool,
        sim_cutoff: float,
        include_list: Union[str, None]) -> str:
    """Wrapper for doing redundancy filtering and domain filtering of clusters

    :param clus_file:
    :param out_folder:
    :param no_redundancy_filtering:
    :param min_genes:
    :param cores:
    :param verbose:
    :param sim_cutoff:
    :param include_list:
    :return: path to filtered clusterfile, containing the domain combinations
        of the filtered bgcs

    Redundancy filtering is based on jaccard overlap of adjacent domain pairs,
    and graph based filtering techniques
    Domain filtering is based on --include_list, e.a. only the biosynthetic
    domains that are used in the paper (biosynthetic_domains.txt)
    """
    random.seed(595)
    dom_dict = read_clusterfile(clus_file, min_genes,
                                verbose)
    doml_dict = {bgc: sum(len(g) for g in genes if not g == ('-',))
                 for bgc, genes in dom_dict.items()}
    filt_file = '{}_filtered_clusterfile.csv'.format(
        clus_file.split('_clusterfile.csv')[0])
    if not os.path.isfile(filt_file):
        # do not perform redundancy filtering if it already exist
        if not no_redundancy_filtering:
            edges_file = generate_edges(dom_dict, sim_cutoff,
                                        cores, out_folder)
            similar_bgcs = read_edges_from_temp(edges_file)
            graph = generate_graph(similar_bgcs, True)
            uniq_bgcs = [clus for clus in dom_dict.keys() if clus not in
                         graph.nodes()]
            all_reps = find_all_representatives(doml_dict, graph)
        else:
            # dont perform redundancy filtering and duplicate clus_file to
            # filt file, representative file is created but this is symbolic
            # todo: remove symbolic (text)
            uniq_bgcs = list(dom_dict.keys())
            all_reps = {}
            print('\nRedundancy filtering is turned off.')
        if include_list:
            print(f"\nOnly domains from {include_list} are included, other "
                  "domains filtered out.")
            include_list = read_txt(include_list)
            dom_dict = filter_out_domains(dom_dict, include_list)
        write_filtered_bgcs(uniq_bgcs, all_reps,
                            dom_dict, filt_file)
    else:
        print('\nFiltered clusterfile existed, (redundancy) filtering not' +
              ' performed again')
    return filt_file


def presto_stat_build_subclusters(
        filt_file: str,
        stat_subclusters_file: str,
        remove_genes_below_count: int,
        min_genes: int,
        cores: int,
        verbose: bool,
        pval_cutoff: float) -> str:
    """Build presto-stat subclusters, and query them to (filtered) train set

    :param filt_file:
    :param stat_subclusters_file:
    :param remove_genes_below_count:
    :param min_genes:
    :param cores:
    :param verbose:
    :param pval_cutoff:
    :return: file containing the filtered final detected modules
    """
    f_clus_dict = read_clusterfile(filt_file, min_genes, verbose)
    if not stat_subclusters_file:
        # run presto-stat to infer sub-clusters from input clusters
        print("\nBuilding PRESTO-STAT sub-clusters from input")
        f_clus_dict_rem = remove_infr_doms(f_clus_dict, min_genes, verbose,
                                           remove_genes_below_count)
        adj_counts, c_counts = count_interactions(f_clus_dict_rem, verbose)
        adj_pvals = calc_adj_pval_wrapper(adj_counts, f_clus_dict_rem, cores,
                                          verbose)
        col_pvals = calc_coloc_pval_wrapper(c_counts, f_clus_dict_rem, cores,
                                            verbose)
        pvals = keep_lowest_pval(col_pvals, adj_pvals)
        # todo: keep from crashing when there are no significant modules
        mods = generate_modules_wrapper(pvals, pval_cutoff, cores,
                                        verbose)
        mod_file = '{}_modules.txt'.format(
            filt_file.split('_filtered_clusterfile.csv')[0])
        write_module_file(mod_file, mods)
        # linking modules to bgcs and filtering mods that occur less than twice
        bgcs_with_mods_ori = q_stat.link_all_mods2bgcs(f_clus_dict_rem, mods,
                                                       cores)
        bgcs_with_mods, modules = remove_infr_mods(bgcs_with_mods_ori, mods)
        mod_file_f = '{}_filtered_modules.txt'.format(
            filt_file.split('_filtered_clusterfile.csv')[0])
        write_module_file(mod_file_f, modules, bgcs_with_mods)

        # todo: assess if this output is necessary
        bgcmodfile = '{}_bgcs_with_mods.txt'.format(
            filt_file.split('_filtered_clusterfile.csv')[0])
        rank_mods = {pair[0]: i + 1 for i, pair in
                     enumerate(sorted(modules.items(),
                                      key=itemgetter(1)))}
        write_bgcs_and_modules(bgcmodfile, f_clus_dict_rem, bgcs_with_mods,
                               rank_mods)
    else:
        # read previously inferred subclusters from file
        print("\nReading PRESTO-STAT subclusters from file:",
              stat_subclusters_file)
        modules = q_stat.read_mods(stat_subclusters_file)
        bgcs_with_mods = q_stat.link_all_mods2bgcs(f_clus_dict, list(modules),
                                                   cores)

    out_file = '{}_presto_stat_subclusters.txt'.format(
        filt_file.split('_filtered_clusterfile.csv')[0])
    print("\nWriting clusters with detected subclusters to", out_file)
    q_stat.write_bgc_mod_fasta(bgcs_with_mods, modules, out_file)
    return out_file


if __name__ == "__main__":
    start = time.time()
    cmd = get_commands()

    # init messages
    if not cmd.include_list:
        print("\n#Warning#: for using models from/replicating the paper, "
              "biosynthetic_domains.txt should be supplied with"
              "--include_list")

    # converting genes in each bgc to a combination of domains
    print("\n1. Preprocessing BGCs into domain combinations")
    cluster_file = preprocessing_bgcs_to_dom_combinations(
        cmd.out_folder,
        cmd.in_folder,
        cmd.hmm_path,
        cmd.start_from_clusterfile,
        cmd.exclude,
        cmd.exclude_contig_edge,
        cmd.min_genes,
        cmd.cores,
        cmd.verbose,
        cmd.use_fastas,
        cmd.use_domtabs,
        cmd.domain_overlap_cutoff)

    # filtering clusters based on similarity
    print("\n2. Filtering the clusters (represented as domain combinations)")
    filtered_cluster_file = filtering_cluster_representations(
        cluster_file,
        cmd.out_folder,
        cmd.no_redundancy_filtering,
        cmd.min_genes,
        cmd.cores,
        cmd.verbose,
        cmd.sim_cutoff,
        cmd.include_list)

    # detecting modules with statistical approach
    print("\n3. PRESTO-STAT - statistical subcluster detection")
    # todo: keep from crashing when no gbks/sufficient doms are present
    bgcs_w_stat_subclusters_file = presto_stat_build_subclusters(
        filtered_cluster_file,
        cmd.stat_subclusters,
        cmd.remove_genes_below_count,
        cmd.min_genes,
        cmd.cores,
        cmd.verbose,
        cmd.pval_cutoff)

    # detecting sub-cluster motifs with topic modelling
    print("\n4. PRESTO-TOP - sub-cluster motif detection with topic modelling")
    presto_top_dir = os.path.join(cmd.out_folder, "presto_top")
    if not os.path.isdir(presto_top_dir):
        os.mkdir(presto_top_dir)

    if not cmd.top_motifs_model:
        print(
            'Parameters: {} topics, {} amplification, '.format(cmd.topics,
                                                               cmd.amplify) +
            '{} iterations of chunksize {}'.format(cmd.iterations,
                                                   cmd.chunksize))
    else:
        print('Parameters: running on existing model at {}'.format(
            cmd.top_motifs_model))

    # writing log information to log.txt
    log_out = os.path.join(presto_top_dir, 'log.txt')
    with open(log_out, 'a') as outf:
        for arg in argv:
            outf.write(arg + '\n')
    logging.basicConfig(filename=log_out,
                        format="%(asctime)s:%(levelname)s:%(message)s",
                        level=logging.INFO)

    bgcs = read2dict(filtered_cluster_file)

    if cmd.classes:
        bgc_classes_dict = read2dict(cmd.classes, sep='\t', header=True)
    else:
        bgc_classes_dict = {bgc: 'None' for bgc in bgcs}

    if not cmd.top_motifs_model:
        bgcs = remove_infr_doms_str(bgcs, cmd.min_genes, False,
                                    cmd.remove_genes_below_count)

    if cmd.amplify:
        bgc_items = []
        for bgc in bgcs.items():
            bgc_items += [bgc] * cmd.amplify
        bgclist, dom_list = zip(*bgc_items)
    else:
        bgclist, dom_list = zip(*bgcs.items())

    if cmd.known_subclusters:
        known_subclusters = defaultdict(list)
        with open(cmd.known_subclusters, 'r') as inf:
            for line in inf:
                line = line.strip().split('\t')
                known_subclusters[line[0]].append(line[1:])
    else:
        known_subclusters = False

    if not cmd.top_motifs_model:
        lda, lda_dict, bow_corpus = run_lda(
            dom_list, no_below=cmd.remove_genes_below_count, no_above=0.5,
            num_topics=cmd.topics, cores=cmd.cores, outfolder=presto_top_dir,
            iters=cmd.iterations, chnksize=cmd.chunksize,
            update_model=cmd.update, ldavis=cmd.visualise)
    else:
        with open(log_out, 'w') as outf:
            outf.write('\nUsing model from {}'.format(cmd.top_motifs_model))
        lda, lda_dict, bow_corpus = run_lda_from_existing(
            cmd.top_motifs_model, dom_list, no_below=1, no_above=0.5)

    process_lda(lda, lda_dict, bow_corpus, cmd.feat_num, bgcs,
                cmd.min_feat_score, bgclist, presto_top_dir, bgc_classes_dict,
                num_topics=cmd.topics, amplif=cmd.amplify, plot=cmd.plot,
                known_subcl=known_subclusters)

    if not cmd.top_motifs_model:
        plot_convergence(log_out, cmd.iterations)

    # visualise subclusters in the output for first 1000 bgcs, otherwise this
    # be very lengthy, surpress with --visualise_subclusters
    print("\n5. Visualising sub-clusters")
    if not cmd.visualise_subclusters:
        print(
            "  use --visualise_subclusters to vis more than first 1000 gbks")
    ipresto_dir = os.path.dirname(os.path.realpath(__file__))
    subcl_arrower = os.path.join(ipresto_dir, "ipresto",
                                 "subcluster_arrower.py")
    # write names of all input gbk to file
    # todo: return from step 1.
    in_gbks = glob(os.path.join(cmd.in_folder, "*.gbk"))
    gbks_path_file = os.path.join(cmd.out_folder, "input_gbks.txt")
    with open(gbks_path_file, "w") as gbk_out:
        if cmd.visualise_subclusters:
            for gbk in in_gbks:
                gbk_out.write(f"{gbk}\n")
        else:
            for gbk in in_gbks[:1000]:
                gbk_out.write(f"{gbk}\n")
    dom_col_file = os.path.join(ipresto_dir, "files",
                                "domains_colour_file.tsv")
    dom_hits_file = '{}_dom_hits.txt'.format(
            filtered_cluster_file.split('_filtered_clusterfile.csv')[0])
    vis_out_file = '{}_ipresto_output_visualisation.html'.format(
            filtered_cluster_file.split('_filtered_clusterfile.csv')[0])
    bgc_topics_filtered = os.path.join(presto_top_dir,
                                       "bgc_topics_filtered.txt")

    vis_cmd = f"python {subcl_arrower} -f {gbks_path_file} " \
              f"-c {dom_col_file} -d {dom_hits_file} " \
              f"-o {vis_out_file} -s {bgcs_w_stat_subclusters_file} " \
              f"-l {bgc_topics_filtered}"
    if cmd.include_list:
        vis_cmd += f" --include_list {cmd.include_list}"

    try:
        subprocess.check_call(vis_cmd, shell=True)
    except subprocess.CalledProcessError as error:
        print("\nVisualising subclusters failed, make sure dom_hits.txt is "
              "present in output as well as the gbks in the input folder. "
              "Alternatively: run ipresto/subcluster_arrower.py yourself")

    end = time.time()
    t = end - start
    t_str = '{}h{}m{}s'.format(int(t / 3600), int(t % 3600 / 60),
                               int(t % 3600 % 60))
    print('\nScript completed in {}'.format(t_str))
