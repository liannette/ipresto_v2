import os
from pathlib import Path

from ipresto.preprocess.tokenize.process_domtable import process_domtables
from ipresto.preprocess.tokenize.process_fasta import process_fastas
from ipresto.preprocess.tokenize.process_gbk import process_gbks


class TokenizeOrchestrator:
    def run(
        self,
        outdir_path,
        gbks_dir_path,
        hmm_file_path,
        exclude_name,
        include_contig_edge_clusters,
        min_genes,
        max_domain_overlap,
        cores,
        verbose,
    ):
        """Wrapper for tokenization of clusters.

        This function tokenizes clusters by processing gbk files into fasta
        files, running hmmscan on the fastas to generate domtables, parsing
        the domtables into tokenized clusters, and writing the tokenized
        clusters to a file.

        Args:
            outdir_path (str): Path to the output directory where intermediate files will be saved.
            gbks_dir_path (str): Path to the folder containing gbk files.
            hmm_file_path (str): Path to the HMM file to be used as the database.
            exclude_name (list of str): List of substrings; files will be excluded if part of the
                file name is present in this list.
            include_contig_edge_clusters (bool): Whether to include clusters on contig edges.
            min_genes (int): Minimum number of genes required for processing.
            max_domain_overlap (float): If two domains overlap more than this
                value, only the domain with the highest score is kept.
            cores (int): Number of CPU cores to use for parallel processing.
            verbose (bool): If True, print additional info to stdout.

        Returns:
            str: Path to the tokenized clusters file.
        """
        # Create output directory
        outdir_path = Path(outdir_path)
        os.makedirs(outdir_path, exist_ok=True)

        # Skip tokenization step if clusters file already exists
        clusters_file_path = outdir_path / "clusters_unfiltered.csv"
        if clusters_file_path.exists():
            if verbose:
                print(
                    "Skipping tokenization step, because clusters file "
                    f"already exists at {clusters_file_path}."
                )
            return clusters_file_path

        # Step 1: Processing gbk files into fasta files
        fastas_dir_path = outdir_path / "fastas"
        os.makedirs(fastas_dir_path, exist_ok=True)
        process_gbks(
            gbks_dir_path,
            fastas_dir_path,
            exclude_name,
            include_contig_edge_clusters,
            cores,
            verbose,
        )

        # Step 2: Processing fastas with hmmscan to generate domtables
        domtables_dir_path = outdir_path / "domtables"
        os.makedirs(domtables_dir_path, exist_ok=True)
        process_fastas(
            fastas_dir_path,
            domtables_dir_path,
            hmm_file_path,
            cores,
            verbose,
        )

        # Step 3: Processing domtables into tokenized clusters
        domain_hits_file_path = outdir_path / "all_domain_hits.txt"
        gene_counts_file_path = outdir_path / "clusters_gene_counts.csv"
        process_domtables(
            domtables_dir_path,
            clusters_file_path,
            gene_counts_file_path,
            domain_hits_file_path,
            min_genes,
            max_domain_overlap,
            cores,
            verbose,
        )

        # Print paths
        if verbose:
            print("\nTokenization complete.")
            print(f"Tokenized clusters have been saved to {clusters_file_path}")
            print(f"Gene counts have been saved to {gene_counts_file_path}")
            print(f"Summary of domain hits has been saved to {domain_hits_file_path}")

        return clusters_file_path
