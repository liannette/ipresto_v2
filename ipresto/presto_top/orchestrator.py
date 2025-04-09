from pathlib import Path
from typing import Optional
from sys import argv
import logging
from ipresto.presto_top.presto_top import (
    read2dict,
    run_lda,
    run_lda_from_existing,
    process_lda,
    plot_convergence,
    remove_infr_doms_str,
)


class TopOrchestrator:
    def run(
        self,
        out_dir_path: str,
        cluster_file: str,
        top_model_file_path: Optional[str],
        min_gene_occurence: int,
        pval_cutoff: float,
        min_genes_per_bgc: int,
        cores: int,
        verbose: bool,
    ) -> str:
        """ """
        # Create output directory
        out_dir = Path(out_dir_path)
        out_dir.mkdir()

        if top_model_file_path is not None:
            if verbose:
                print(
                    f"\nUsing existing model from: {top_model_file_path}"
                )
        else:
            if verbose:
                print(
                    "\nComputing new PRESTO-TOP model from the input clusters: "
                    f"{cluster_file}"
                )
                print(
                    f"Parameters: {topics} topics, {amplify} amplification, "
                    f"{iterations} iterations of chunksize {chunksize}"
                )

        # writing log information to log.txt
        log_out = out_dir / 'log.txt'
        with open(log_out, 'a') as outfile:
            for arg in argv:
                outfile.write(arg + '\n')
        logging.basicConfig(filename=log_out,
                            format="%(asctime)s:%(levelname)s:%(message)s",
                            level=logging.INFO)

        bgcs = read2dict(filtered_cluster_file)

        if cmd.classes:
            bgc_classes_dict = read2dict(cmd.classes, sep='\t', header=True)
        else:
            bgc_classes_dict = {bgc: 'None' for bgc in bgcs}

        if not cmd.top_motifs_model:
            bgcs = remove_infr_doms_str(bgcs, cmd.min_genes, cmd.verbose,
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
                update_model=cmd.update, ldavis=cmd.visualise, alpha=cmd.alpha,
                beta=cmd.beta)
        else:
            with open(log_out, 'w') as outf:
                outf.write('\nUsing model from {}'.format(cmd.top_motifs_model))
            lda, lda_dict, bow_corpus = run_lda_from_existing(
                cmd.top_motifs_model, dom_list, presto_top_dir,
                no_below=1, no_above=0.5)

        process_lda(lda, lda_dict, bow_corpus, cmd.feat_num, bgcs,
                    cmd.min_feat_score, bgclist, presto_top_dir, bgc_classes_dict,
                    num_topics=cmd.topics, amplif=cmd.amplify, plot=cmd.plot,
                    known_subcl=known_subclusters)

        if not cmd.top_motifs_model:
            plot_convergence(log_out, cmd.iterations)