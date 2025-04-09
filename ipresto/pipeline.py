from pathlib import Path

from ipresto.preprocess.orchestrator import PreprocessOrchestrator
from ipresto.presto_stat.orchestrator import StatOrchestrator
#from ipresto.presto_top.orchestrator import TopOrchestrator

class IprestoPipeline:
    def run(
        self,
        out_dir_path,
        gbks_dir_path,
        existing_clusterfile,
        exclude_name,
        include_contig_edge_clusters,
        hmm_file_path,
        max_domain_overlap,
        min_genes_per_bgc,
        domain_filtering,
        similarity_filtering,
        similarity_cutoff,
        stat_modules_file_path,
        stat_min_gene_occurrence,
        stat_pval_cutoff,
        cores,
        verbose,
    ):
        """
        Runs the entire pipeline.
        """
        out_dir_path = Path(out_dir_path)
        out_dir_path.mkdir(parents=True, exist_ok=True)

        # Step 1: Preprocessing clusters, skip if a cluster file is provided
        if verbose:
            print("\n=== Preprocessing BGCs ===")

        if existing_clusterfile:
            clusters_file_path = Path(existing_clusterfile)
            if verbose:
                print(
                    f"\nUsing provided file of preprocessed BGCs: {clusters_file_path}."
                )
        else:
            preprocess_dir_path = out_dir_path / "preprocess"
            clusters_file_path = PreprocessOrchestrator().run(
                preprocess_dir_path,
                gbks_dir_path,
                hmm_file_path,
                exclude_name,
                include_contig_edge_clusters,
                min_genes_per_bgc,
                max_domain_overlap,
                domain_filtering,
                similarity_filtering,
                similarity_cutoff,
                cores,
                verbose,
            )

        # Step 2: Statistical subcluster detection (PRESTO-STAT)
        if verbose:
            print("\n=== PRESTO-STAT: statistical subcluster detection ===")

        stat_dir_path = out_dir_path / "stat_subclusters"
        StatOrchestrator().run(
            stat_dir_path,
            clusters_file_path,
            stat_modules_file_path,
            stat_min_gene_occurrence,
            stat_pval_cutoff,
            min_genes_per_bgc,
            cores,
            verbose,
        )

        # if verbose:
        #     print("\n=== PRESTO-TOP: sub-cluster motif detection with topic modelling ===")

        # top_dir_path = out_dir_path / "top_subclusters"
        # TopOrchestrator().run(
        #     top_dir_path,
        #     clusters_file_path,
        #     stat_dir_path / "stat_modules",
        #     stat_modules_file_path,
        #     stat_min_gene_occurrence,
        #     cores,
        #     verbose,
        # )

        # save_results(output_path, results)

        # if visualize:
        #     if verbose:
        #         print("=== Visualizations ===")
        #     plot_histogram(load_data(preprocessed_data_path), output_path="histogram.png")
        #     plot_comparison(results, output_path="comparison.png")
