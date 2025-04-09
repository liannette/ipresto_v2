import random
from pathlib import Path
from typing import List

from ipresto.preprocess.tokenize.orchestrator import TokenizeOrchestrator
from ipresto.preprocess.domain_filtering import perform_domain_filtering
from ipresto.preprocess.similarity_filtering import filter_similar_clusters


class PreprocessOrchestrator:
    def run(
        self,
        out_dir_path: str,
        gbks_dir_path: str,
        hmm_file_path: str,
        exclude_name: List[str],
        include_contig_edge_clusters: bool,
        min_genes: int,
        max_domain_overlap: float,
        domain_filtering: bool,
        similarity_filtering: bool,
        sim_cutoff: float,
        cores: int,
        verbose: bool,
    ) -> str:
        """Orchestrate all preprocessing steps.

        Preprocess BGCs by tokenizing, filtering non-biosynthetic domains, and removing similar clusters.

        Args:
            out_dir_path (str): Output directory path.
            gbks_dir_path (str): Directory path containing GenBank files.
            hmm_file_path (str): Path to HMM file.
            exclude_name (List[str]): If any string in this list occurs in the gbk filename, this
                file will not be used for the analysis.
            include_contig_edge_clusters (bool): Whether to include clusters that lie on a contig edge.
            min_genes (int): Minimum number of genes.
            cores (int): Number of cores to use.
            verbose (bool): Whether to print verbose output.
            max_domain_overlap (float): If two domains overlap more than this value, only the domain with the highest score is kept.
            domain_filtering (bool): Whether to filter non-biosynthetic protein domains.
            similarity_filtering (bool): Whether to filter similar clusters.
            sim_cutoff (float): Similarity cutoff value for filtering clusters.

        Returns:
            str: Path to the final preprocessed clusters file.
        """
        # Set random seed for reproducibility
        random.seed(595)

        # Create output directory
        out_dir = Path(out_dir_path)
        out_dir.mkdir(parents=True, exist_ok=True)

        # Step 1: Tokenize the genes of the clusters into protein domain combinations
        if verbose:
            print("\nTokenizing the BGC genes into protein domain combinations")
        clusters_file_path = TokenizeOrchestrator().run(
            out_dir_path,
            gbks_dir_path,
            hmm_file_path,
            exclude_name,
            include_contig_edge_clusters,
            min_genes,
            max_domain_overlap,
            cores,
            verbose,
        )

        # Step 2: Filter non-biosynthetic protein domains
        if domain_filtering is False:
            if verbose:
                print("\nSkipping domain filtering, because it has been turned off.")
        else:
            out_file_path = out_dir / "clusters_domfiltered.csv"
            if out_file_path.is_file():
                if verbose:
                    print(
                        f"\nSkipping domain filtering, because the file already exists: {out_file_path}"
                    )
            else:
                domain_filtering_file_path = (
                    Path(__file__).parent.parent.parent
                    / "data"
                    / "biosynthetic_domains.txt"
                )
                counts_file_path = out_dir / "clusters_domfiltered_gene_counts.txt"
                perform_domain_filtering(
                    clusters_file_path,
                    domain_filtering_file_path,
                    out_file_path,
                    counts_file_path,
                    min_genes,
                    cores,
                    verbose,
                )
            clusters_file_path = out_file_path

        # Step 3: Filter similar clusters
        if similarity_filtering is False:
            if verbose:
                print(
                    "\nSkipping similarity filtering, because it has been turned off."
                )
        else:
            out_file_path = out_dir / "clusters_representatives.csv"
            if out_file_path.is_file():
                if verbose:
                    print(
                        f"\nSkipping similarity filtering, because the file already exists: {out_file_path}"
                    )
            else:
                counts_file_path = out_dir / "clusters_representatives_gene_counts.txt"
                representatives_file_path = out_dir / "representative_clusters.txt"
                edge_file_path = out_dir / "edges_clusters.txt"
                clusters_file_path = filter_similar_clusters(
                    clusters_file_path,
                    out_file_path,
                    counts_file_path,
                    representatives_file_path,
                    edge_file_path,
                    sim_cutoff,
                    cores,
                    verbose,
                )
            clusters_file_path = out_file_path

        return clusters_file_path
