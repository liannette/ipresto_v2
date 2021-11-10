#!/usr/bin/env python3
"""
Author: Joris Louwen (joris.louwen@wur.nl)

Part of iPRESTO, Bioinformatics group Wageningen University.
PI: Marnix Medema, Justin van der Hooft

usage:
python ipresto.py -h
"""

from ipresto.presto_stat.presto_stat import *
from ipresto.presto_top.presto_top import *

if __name__ == "__main__":
    start = time.time()
    cmd = get_commands()

    # converting genes in each bgc to a combination of domains
    if cmd.start_from_clusterfile:
        if not os.path.isdir(cmd.out_folder):
            f_command = 'mkdir {}'.format(cmd.out_folder)
            subprocess.check_call(f_command, shell=True)
        filepre = os.path.split(cmd.start_from_clusterfile)[-1].split(
            '.csv')[0]
        clus_file = os.path.join(cmd.out_folder, filepre + '_clusterfile.csv')
        c_command = 'cp {} {}'.format(cmd.start_from_clusterfile, clus_file)
        subprocess.check_call(c_command, shell=True)
    else:
        fasta_folder, exist_fastas = process_gbks(
            cmd.in_folder, cmd.out_folder, cmd.exclude,
            cmd.exclude_contig_edge, cmd.min_genes, cmd.cores, cmd.verbose,
            cmd.use_fastas)
        dom_folder, exist_doms = hmmscan_wrapper(
            fasta_folder, cmd.hmm_path, cmd.verbose, cmd.cores, exist_fastas,
            cmd.use_domtabs)
        clus_file = parse_dom_wrapper(dom_folder, cmd.out_folder,
                                      cmd.domain_overlap_cutoff, cmd.verbose,
                                      exist_doms)

    # filtering clusters based on similarity
    random.seed(595)
    dom_dict = read_clusterfile(clus_file, cmd.min_genes,
                                cmd.verbose)
    doml_dict = {bgc: sum(len(g) for g in genes if not g == ('-',))
                 for bgc, genes in dom_dict.items()}
    filt_file = '{}_filtered_clusterfile.csv'.format(
        clus_file.split('_clusterfile.csv')[0])
    if not os.path.isfile(filt_file):
        # do not perform redundancy filtering if it already exist
        if not cmd.no_redundancy_filtering:
            edges_file = generate_edges(dom_dict, cmd.sim_cutoff,
                                        cmd.cores, cmd.out_folder)
            similar_bgcs = read_edges_from_temp(edges_file)
            graph = generate_graph(similar_bgcs, True)
            uniq_bgcs = [clus for clus in dom_dict.keys() if not clus in
                                                                 graph.nodes()]
            all_reps = find_all_representatives(doml_dict, graph)
        else:
            # dont perform redundancy filtering and duplicate clus_file to
            # filt file, representative file is created but this is symbolic
            uniq_bgcs = list(dom_dict.keys())
            all_reps = {}
            print('\nRedundancy filtering is turned off.')
        if cmd.include_list:
            include_list = read_txt(cmd.include_list)
            dom_dict = filter_out_domains(dom_dict, include_list)
        all_reps_file = write_filtered_bgcs(uniq_bgcs, all_reps,
                                            dom_dict, filt_file)
    else:
        print('\nFiltered clusterfile existed, redundancy filtering not' +
              ' performed again')

    # detecting modules with statistical approach
    f_clus_dict = read_clusterfile(filt_file, cmd.min_genes, cmd.verbose)
    f_clus_dict_rem = remove_infr_doms(f_clus_dict, cmd.min_genes, cmd.verbose)
    adj_counts, c_counts = count_interactions(f_clus_dict_rem, cmd.verbose)
    adj_pvals = calc_adj_pval_wrapper(adj_counts, f_clus_dict_rem, cmd.cores, \
                                      cmd.verbose)
    col_pvals = calc_coloc_pval_wrapper(c_counts, f_clus_dict_rem, cmd.cores, \
                                        cmd.verbose)
    pvals = keep_lowest_pval(col_pvals, adj_pvals)
    mods = generate_modules_wrapper(pvals, cmd.pval_cutoff, cmd.cores, \
                                    cmd.verbose)
    mod_file = '{}_modules.txt'.format( \
        filt_file.split('_filtered_clusterfile.csv')[0])
    write_module_file(mod_file, mods)
    # linking modules to bgcs and filtering mods that occur less than twice
    bgcs_with_mods_ori = link_all_mods2bgcs(f_clus_dict_rem, mods, cmd.cores)
    bgcs_with_mods, modules = remove_infr_mods(bgcs_with_mods_ori, mods)
    mod_file_f = '{}_filtered_modules.txt'.format( \
        filt_file.split('_filtered_clusterfile.csv')[0])
    write_module_file(mod_file_f, modules, bgcs_with_mods)
    bgcmodfile = '{}_bgcs_with_mods.txt'.format( \
        mod_file.split('_modules.txt')[0])
    rank_mods = {pair[0]: i + 1 for i, pair in
                 enumerate(sorted(modules.items(), \
                                  key=itemgetter(1)))}
    write_bgcs_and_modules(bgcmodfile, f_clus_dict_rem, bgcs_with_mods, \
                           rank_mods)

    end = time.time()
    t = end - start
    t_str = '{}h{}m{}s'.format(int(t / 3600), int(t % 3600 / 60),
                               int(t % 3600 % 60))
    print('\nScript completed in {}'.format(t_str))
