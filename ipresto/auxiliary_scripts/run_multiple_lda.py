#!/usr/bin/env python3
"""
Author: Joris Louwen
Script to run LDA multiple times to find good parameters for applying LDA
to identify sub-clusters.
"""

import subprocess
from sys import argv
import os
import time

if __name__ == '__main__':
    print("Start - output will be written to", os.getcwd())
    start = time.time()
    if len(argv) < 2:
        raise ValueError("Usage: python run_multiple_lda.py <path_to_ipresto>")
    ipresto_path = argv[1]

    # see if ipresto is found correctly
    test_command = f"python {ipresto_path} -h"
    try:
        subprocess.check_call(test_command, shell=True)
    except subprocess.CalledProcessError as e:
        raise FileNotFoundError("Failed to find <ipresto_path>:", e)

    # hardcoding antismashdb files
    base = "/mnt/scratch/louwe015/subcluster_data"
    clust_file = os.path.join(
        base, "antismashdb_crusemann_mibig_filtered_clusterfile.csv")
    hmm = os.path.join(base, 'domains', 'Pfam_100subs_tc.hmm')
    chunksize = '3000'
    cores = '10'
    iters = '500'
    classes = os.path.join(base, 'all_classes.txt')
    known_subcl = '/mnt/scratch/louwe015/subcluster_data/subclusterblast' \
                  '_data/subclusters_subclusterblast_domains_synt_subset.txt'
    biosynt_domains = os.path.join(os.path.split(ipresto_path)[0],
                                   'files', 'biosynthetic_domains.txt')
    stat_subcl = '/mnt/scratch/louwe015/iPRESTO/test_ipresto_params/' \
                 'dummy_presto_stat_subclusters.txt'
    logfile = 'log_{}.txt'
    command = (
        'python {} -i bla -o {} --start_from_clusterfile {} -t {} -C {} -I {} '
        '-c {} --classes {} --known_subclusters {} --include_list {} '
        '--min_genes 2 --no_redundancy_filtering --visualise '
        '--stat_subclusters {} --hmm_path {} --alpha {} --beta {} 2>&1 '
        '| tee -a {}')
    topic_range = []  # [1, 10, 100, 250, 500, 750, 1000, 1500]

    for topics in topic_range:
        alpha = 'symmetric'
        beta = None
        output_dir = (
            f'ipresto_out_topics-{topics}_alpha-{alpha}_beta-{beta}')
        if os.path.isdir(output_dir):
            print(f'\n###{output_dir} already exists, '
                  f'model not trained again\n')
        form_logfile = logfile.format(output_dir)
        formatted_command = command.format(
            ipresto_path, output_dir, clust_file, topics,
            chunksize, iters, cores, classes, known_subcl,
            biosynt_domains, stat_subcl, hmm, alpha, beta, form_logfile)

        # run the command
        print(f"\n###Running, {formatted_command}\n")
        try:
            subprocess.check_call(formatted_command, shell=True)
        except subprocess.CalledProcessError as com_err:
            print(
                f"Command '{formatted_command}' failed: {com_err}")

    # alpha beta opt
    alpha_range = ['symmetric', 'asymmetric', '1']
    beta_range = [None, 'auto', '1']  # also called eta
    topics = 1000
    for alpha in alpha_range:
        # beta = None
        for beta in beta_range:
            output_dir = (
                f'ipresto_out_topics-{topics}_alpha-{alpha}_beta-{beta}')
            if os.path.isdir(output_dir):
                print(f'\n###{output_dir} already exists, '
                      f'model not trained again\n')
            form_logfile = logfile.format(output_dir)
            formatted_command = command.format(
                ipresto_path, output_dir, clust_file, topics,
                chunksize, iters, cores, classes, known_subcl,
                biosynt_domains, stat_subcl, hmm, alpha, beta, form_logfile)

            # run the command
            print(f"\n###Running, {formatted_command}\n")
            try:
                subprocess.check_call(formatted_command, shell=True)
            except subprocess.CalledProcessError as com_err:
                print(
                    f"Command '{formatted_command}' failed: {com_err}")

    # for beta in beta_range:
    #     alpha = 'symmetric'
    #     output_dir = (
    #         f'ipresto_out_topics-{topics}_alpha-{alpha}_beta-{beta}')
    #     if os.path.isdir(output_dir):
    #         print(f'\n###{output_dir} already exists, '
    #               f'model not trained again\n')
    #     form_logfile = logfile.format(output_dir)
    #     formatted_command = command.format(
    #         ipresto_path, output_dir, clust_file, topics,
    #         chunksize, iters, cores, classes, known_subcl,
    #         biosynt_domains, stat_subcl, hmm, alpha, beta, form_logfile)
    #
    #     # run the command
    #     print(f"\n###Running, {formatted_command}\n")
    #     try:
    #         subprocess.check_call(formatted_command, shell=True)
    #     except subprocess.CalledProcessError as com_err:
    #         print(
    #             f"Command '{formatted_command}' failed: {com_err}")

    end = time.time()
    print(f"\nCompleted in {end-start:.1f}s")
