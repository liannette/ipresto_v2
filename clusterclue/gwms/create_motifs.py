import logging
from collections import Counter
from itertools import product
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
import pickle
import json

from clusterclue.classes.subcluster_motif import SubclusterMotif
from clusterclue.clusters.utils import read_clusters
from clusterclue.gwms.create.build_gwms import build_motif_gwms, write_motif_gwms
from clusterclue.gwms.create.clustering import ClusteringComparison
from clusterclue.gwms.create.combine_matches import combine_presto_matches
from clusterclue.gwms.create.merge_motifs import merge_similar_motifs
from clusterclue.gwms.create.reports import generate_evaluation_report
from clusterclue.gwms.create.plots import generate_clustering_plots, plot_evaluation_summary
from clusterclue.evaluate.evaluate_hits import (
    assign_best_hit,
    read_reference_subclusters_and_tokenize_genes
)
from clusterclue.gwms.detect_motifs import (
    detect_motifs, write_motif_hits, parse_clusters_file, parse_motifs_file
)


logger = logging.getLogger(__name__)


def get_gene_background_count(clusters: dict) -> Counter:
    """Counts how many BGCs each tokenised gene occurs in."""
    gene_counts = Counter()
    for genes in clusters.values():
        tokenized_genes = set([';'.join(gene) for gene in genes])
        gene_counts.update(tokenized_genes)
    gene_counts.pop("-", None) # remove genes without biosynthetic domains
    return gene_counts


def convert_numpy_types(obj):
    """Recursively convert numpy types to native Python types for JSON serialization."""
    if isinstance(obj, dict):
        return {key: convert_numpy_types(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(item) for item in obj]
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj


def save_module_label_mapping(modules: List[str], labels: np.ndarray, filepath: Path):
    """Save module to cluster label mapping."""
    mapping = pd.DataFrame({
        'module': [','.join(mod) for mod in modules],
        'cluster_label': labels
    })
    mapping.to_csv(filepath, sep='\t', index=False)
    logger.info(f"Saved module-label mapping to {filepath}")
    return mapping


def analyze_motifs(motifs: Dict[str, SubclusterMotif]) -> Dict:
    """Analyze motif statistics."""
    if not motifs:
        return {}
    
    match_counts = [m.n_matches for m in motifs.values()]
    all_genes = set()
    for m in motifs.values():
        all_genes.update(m.tokenized_genes)
    
    stats = {
        'total_motifs': len(motifs),
        'total_unique_genes': len(all_genes),
        'match_count_mean': np.mean(match_counts),
        'match_count_median': np.median(match_counts),
        'match_count_std': np.std(match_counts),
        'match_count_min': np.min(match_counts),
        'match_count_max': np.max(match_counts),
    }
    
    # Only include core gene stats if core_genes attribute exists (after GWM building)
    if hasattr(next(iter(motifs.values())), 'core_genes'):
        core_gene_counts = [len(m.core_genes) for m in motifs.values()]
        stats.update({
            'core_genes_mean': np.mean(core_gene_counts),
            'core_genes_median': np.median(core_gene_counts),
            'core_genes_std': np.std(core_gene_counts),
            'core_genes_min': np.min(core_gene_counts),
            'core_genes_max': np.max(core_gene_counts),
        })
    
    return stats



# ============================================================================
# STEP 1: CLUSTERING
# ============================================================================

def run_clustering(
    configs: List[dict],
    stat_matches_filepath: Path,
    top_matches_filepath: Path,
    out_dirpath: Path,
    n_jobs: int
) -> dict:
    """
    Run clustering with multiple configurations.
    
    Args:
        configs: List of clustering configurations
        stat_matches_filepath: Path to STAT matches
        top_matches_filepath: Path to TOP matches
        out_dirpath: Output directory
        n_jobs: Number of parallel jobs
    
    Returns:
        dict: Clustering results for each configuration
    """
    logger.info("=== Clustering subcluster modules ===")

    # Combine PRESTO matches
    combined_matches = combine_presto_matches(
        stat_matches_filepath, 
        top_matches_filepath,
    )
    
    # Create module to BGC mapping
    module2bgcs = defaultdict(list)
    for bgc_id, module in combined_matches:
        module2bgcs[module].append(bgc_id)
    modules = sorted(module2bgcs.keys())

    logger.info(f"Loaded {len(modules)} unique modules")

    # Initialize comparison
    comparison = ClusteringComparison(modules, module2bgcs)

    results = {}

    for config in configs:
        name = config['name']
        logger.info(f"Running clustering: {name}")

        # Create config-specific output directory
        config_dir = out_dirpath / name
        config_dir.mkdir(parents=True, exist_ok=True)

        try:
            # Run clustering
            if config['method'] == 'kmeans':
                labels, metadata = comparison.run_kmeans(
                    k=config['k'],
                    use_svd=config.get('use_svd', False),
                    target_variance=config.get('target_variance', 0.30)
                )
            elif config['method'] == 'hdbscan':
                labels, metadata = comparison.run_hdbscan(
                    min_cluster_size=config['min_cluster_size'],
                    use_svd=config.get('use_svd', True),
                    target_variance=config.get('target_variance', 0.30),
                    min_samples=config.get('min_samples', None),
                    cluster_selection_epsilon=config.get('cluster_selection_epsilon', 0.0),
                    n_jobs=n_jobs
                )
            else:
                logger.error(f"Unknown method: {config['method']}")
                continue

            # Save module-label mapping
            save_module_label_mapping(modules, labels, config_dir / "module_labels.tsv")

            # Calculate cluster quality metrics
            cluster_metrics = comparison.calculate_cluster_metrics(labels)
            metadata.update(cluster_metrics)

            # Convert to motifs
            motifs = comparison.labels_to_motifs(labels)
            logger.info(f"Created {len(motifs)} initial motifs")
            
            # Save motifs
            motifs_filepath = config_dir / f"motifs_{name}.pkl"
            with open(motifs_filepath, "wb") as f:
                pickle.dump(motifs, f)
            logger.info(f"Saved initial motifs to {motifs_filepath}")

            # Analyze motifs
            motifs_stats = analyze_motifs(motifs)

            # Store results
            results[name] = {
                'config': config,
                'metadata': metadata,
                'n_motifs': len(motifs),
                'motifs_stats': motifs_stats,
                'motifs_filepath': str(motifs_filepath),
                'config_dir': str(config_dir),
            }

            # Generate plots
            try:
                if metadata.get('use_svd', False):
                    X_plot, _, _ = comparison.apply_svd(metadata.get('target_variance', 0.30))
                else:
                    X_plot = comparison.X.toarray()
                
                generate_clustering_plots(name, X_plot, labels, metadata, config_dir)
            except Exception as e:
                logger.warning(f"Error generating plots for {name}: {e}")
        
        except Exception as e:
            logger.error(f"Error in config {name}: {e}", exc_info=True)
            results[name] = {'error': str(e)}

    # Save clustering results
    results_file = out_dirpath / "results_clustering.json"
    with open(results_file, 'w') as f:
        json.dump(convert_numpy_types(results), f, indent=2)
    logger.info(f"Saved clustering results to {results_file}")

    return results


# ============================================================================
# STEP 2: MERGING
# ============================================================================

def run_merging(
    clustering_results: dict,
    merge_params: List[Tuple[float, float, str]],
    out_dirpath: Path
) -> dict:
    """
    Merge motifs with multiple parameter configurations.
    
    Args:
        clustering_results: Results from run_clustering
        merge_params: List of (similarity_threshold, gene_threshold, metric) tuples
        out_dirpath: Output directory
    
    Returns:
        dict: Merge results for each clustering + merge combination
    """
    logger.info("=== Running motif merging ===")
    
    results = {}
    
    for cluster_name, cluster_result in clustering_results.items():
        if 'error' in cluster_result:
            logger.warning(f"Skipping {cluster_name} due to clustering error")
            continue
        
        # Load initial motifs
        motifs_filepath = Path(cluster_result['motifs_filepath'])
        with open(motifs_filepath, "rb") as f:
            motifs = pickle.load(f)
        logger.info(f"Loaded {len(motifs)} motifs for {cluster_name}")
        
        config_dir = Path(cluster_result['config_dir'])
        
        # Test each merge parameter combination
        for sim_threshold, gene_threshold, metric in merge_params:
            merge_name = f"sim{int(sim_threshold*100)}_gene{int(gene_threshold*100)}_{metric}"
            full_name = f"{cluster_name}_{merge_name}"
            
            logger.info(f"Merging with {merge_name}")
            
            # Merge similar motifs
            merged_motifs = merge_similar_motifs(
                motifs, 
                similarity_threshold=sim_threshold,
                gene_prob_threshold=gene_threshold,
                similarity_metric=metric
            )
            logger.info(f"Reduced {len(motifs)} → {len(merged_motifs)} motifs")

            # Save merged motifs
            merged_motifs_filepath = config_dir / f"motifs_merged_{merge_name}.pkl"
            with open(merged_motifs_filepath, "wb") as f:
                pickle.dump(merged_motifs, f)

            # Analyze merged motifs
            merged_stats = analyze_motifs(merged_motifs)

            # Store results
            results[full_name] = {
                'config': cluster_result['config'],
                'metadata': cluster_result['metadata'],
                'merge_params': {
                    'similarity_threshold': sim_threshold,
                    'gene_threshold': gene_threshold,
                    'metric': metric,
                },
                'n_initial_motifs': cluster_result['n_motifs'],
                'n_merged_motifs': len(merged_motifs),
                'merged_motifs_stats': merged_stats,
                'motifs_filepath': str(motifs_filepath),
                'merged_motifs_filepath': str(merged_motifs_filepath),
                'config_dir': str(config_dir),
            }
    
    # Save merge results
    results_file = out_dirpath / "results_merging.json"
    with open(results_file, 'w') as f:
        json.dump(convert_numpy_types(results), f, indent=2)
    logger.info(f"Saved merge results to {results_file}")
    
    return results


# ============================================================================
# STEP 3: GWM BUILDING
# ============================================================================

def run_gwm_building(
    merge_results: dict,
    gwm_base_params: List[Tuple[int, int, float]],
    clusters_filepath: Path,
    out_dirpath: Path
) -> dict:
    """
    Build GWMs for all merge configurations with multiple GWM parameters.
    
    Args:
        merge_results: Results from run_merging
        gwm_base_params: List of (min_matches, min_core_genes, core_threshold) tuples
        clusters_filepath: Path to tokenized clusters
        out_dirpath: Output directory
    
    Returns:
        dict: GWM results for each configuration
    """
    logger.info("=== Building GWMs ===")
    
    gwms_dirpath = out_dirpath / "gwms"
    gwms_dirpath.mkdir(parents=True, exist_ok=True)

    # Load training data
    training_clusters = read_clusters(clusters_filepath)
    n_clusters_total = len(training_clusters)
    background_counts = get_gene_background_count(training_clusters)

    results = {}
    
    for merge_name, merge_result in merge_results.items():
        if 'error' in merge_result:
            logger.warning(f"Skipping {merge_name} due to previous error")
            continue

        # Extract gene_prob used in merging
        merge_gene_prob = merge_result['merge_params']['gene_threshold']

        # Load merged motifs
        merged_motifs_filepath = Path(merge_result['merged_motifs_filepath'])
        with open(merged_motifs_filepath, "rb") as f:
            merged_motifs = pickle.load(f)
        logger.info(f"Loaded {len(merged_motifs)} merged motifs for {merge_name}")

        # Test GWM params, using the same gene_prob as merging
        for mm, mgc, ct in gwm_base_params:
            params_str = f"mm{mm}_mgc{mgc}_ct{int(ct * 100)}_mgp{int(merge_gene_prob * 100)}"
            full_name = f"{merge_name}_{params_str}"

            logger.info(f"Building GWMs: {params_str}")
                
            gwms = build_motif_gwms(
                merged_motifs,
                background_counts,
                n_clusters_total,
                min_matches=mm,
                min_core_genes=mgc,
                core_threshold=ct,
                min_gene_prob=merge_gene_prob,
            )
            logger.info(f"Built {len(gwms)} GWMs")

            # Save GWMs
            gwm_filepath = gwms_dirpath / f"GWMs_{full_name}.txt"
            write_motif_gwms(gwms, gwm_filepath)

            # Analyze GWMs
            gwms_stats = analyze_motifs(gwms)
        
            # Store results
            results[full_name] = {
                'config': merge_result['config'],
                'metadata': merge_result['metadata'],
                'merge_params': merge_result['merge_params'],
                'n_initial_motifs': merge_result['n_initial_motifs'],
                'n_merged_motifs': merge_result['n_merged_motifs'],
                'motifs_filepath': merge_result['motifs_filepath'],
                'merged_motifs_filepath': merge_result['merged_motifs_filepath'],
                'gwm_hyperparameter': {
                    'min_matches': mm,
                    'min_core_genes': mgc,
                    'core_threshold': ct,
                    'min_gene_prob': merge_gene_prob,
                },
                'n_gwms': len(gwms),
                'gwms_stats': gwms_stats,
                'gwms_filepath': str(gwm_filepath),
            }
    
    # Save GWM results
    results_file = out_dirpath / "results_gwms.json"
    with open(results_file, 'w') as f:
        json.dump(convert_numpy_types(results), f, indent=2)
    logger.info(f"Saved GWM results to {results_file}")
    
    return results


# ============================================================================
# STEP 4: EVALUATION
# ============================================================================

def run_evaluation(
    gwm_results: dict,
    annotated_subclusters_filepath: Path,
    reference_clusters_filepath: Path,
    overlap_penalty_alpha: float,
    overlap_penalty_beta: float,
    out_dirpath: Path
) -> dict:
    """
    Evaluate GWMs against reference subclusters.
    
    Args:
        gwm_results: Results from run_gwm_building
        annotated_subclusters_filepath: Path to annotated subclusters
        reference_clusters_filepath: Path to reference clusters
        overlap_penalty_alpha: Alpha parameter for overlap penalty
        overlap_penalty_beta: Beta parameter for overlap penalty
        out_dirpath: Output directory
    
    Returns:
        dict: Best result based on MRPOS
    """
    logger.info("=== Evaluating GWMs ===")
    
    # Load reference data
    annotated_subclusters = read_reference_subclusters_and_tokenize_genes(
        annotated_subclusters_filepath, 
        reference_clusters_filepath.parent / "all_domain_hits.txt"
    )
    logger.info(f"Loaded {len(annotated_subclusters)} annotated subclusters")

    ref_clusters = parse_clusters_file(reference_clusters_filepath)
    logger.info(f"Loaded {len(ref_clusters)} reference clusters")

    # Create output directories
    hits_dirpath = out_dirpath / "evaluation_hits"
    hits_dirpath.mkdir(parents=True, exist_ok=True)

    evaluation_dirpath = out_dirpath / "evaluation_best_hits"
    evaluation_dirpath.mkdir(parents=True, exist_ok=True)

    results = gwm_results.copy()

    for gwm_name, gwm_result in results.items():
        logger.info(f"Evaluating {gwm_name}")

        # Load GWMs
        gwm_filepath = Path(gwm_result['gwms_filepath'])
        gwms = parse_motifs_file(gwm_filepath)
        
        # Detect motifs
        hits = detect_motifs(ref_clusters, gwms)
        logger.info(f"Found {len(hits)} motif hits")

        # Save hits
        eval_hits_filepath = hits_dirpath / f"hits_{gwm_name}.txt"
        write_motif_hits(hits, gwms, eval_hits_filepath)

        # Evaluate against annotated subclusters
        eval_df = pd.DataFrame(
            annotated_subclusters.apply(
                lambda row: assign_best_hit(
                    row, hits, 
                    alpha=overlap_penalty_alpha, 
                    beta=overlap_penalty_beta
                ), 
                axis=1
            ).tolist()
        )
        
        # Save evaluation
        eval_best_hits_filepath = evaluation_dirpath / f"best_hits_{gwm_name}.txt"
        eval_df.to_csv(eval_best_hits_filepath, sep="\t", index=False)

        # Calculate scores
        mean_overlap = eval_df["overlap_score"].mean()
        mean_penalized = eval_df["penalized_score"].mean()
        
        logger.info(f"Overlap: {mean_overlap:.4f}, MRPOS: {mean_penalized:.4f}")
        
        # Update results
        results[gwm_name].update({
            'mean_overlap_score': float(mean_overlap),
            'mean_penalized_score': float(mean_penalized),
            'alpha': overlap_penalty_alpha,
            'beta': overlap_penalty_beta,
            'eval_hits_filepath': str(eval_hits_filepath),
            'eval_best_hits_filepath': str(eval_best_hits_filepath),
        })

    # Save final results
    results_file = out_dirpath / "results_final.json"
    with open(results_file, 'w') as f:
        json.dump(convert_numpy_types(results), f, indent=2)
    logger.info(f"Saved final results to {results_file}")

    # Find best configuration
    best_name, best_result = max(
        results.items(), 
        key=lambda x: x[1].get('mean_penalized_score', 0)
    )
    logger.info(f"Best configuration: {best_name}")
    logger.info(f"  Overlap: {best_result['mean_overlap_score']:.4f}")
    logger.info(f"  MRPOS: {best_result['mean_penalized_score']:.4f}")
    logger.info(f"  GWMs: {best_result['n_gwms']:,}")

    # Generate report
    report_filepath = out_dirpath / "evaluation_report.txt"
    generate_evaluation_report(results, report_filepath)
    logger.info(f"Saved evaluation report to {report_filepath}")

    best_result["name"] = best_name
    return best_result


# ============================================================================
# MAIN PIPELINE FUNCTIONS
# ============================================================================

def create_and_evaluate_motif_gwms(
    stat_matches_filepath: Path,
    top_matches_filepath: Path,
    clusters_filepath: Path,
    annotated_subclusters_filepath: Path,
    reference_clusters_filepath: Path,
    overlap_penalty_alpha: float,
    overlap_penalty_beta: float,
    out_dirpath: Path,
    n_jobs: int,
    clustering_configs: List[dict],
    merge_similarity_thresholds: List[float],
    merge_gene_thresholds: List[float],
    merge_metrics: List[str],
    gwm_min_matches: List[int],
    gwm_min_core_genes: List[int],
    gwm_core_thresholds: List[float],
    gwm_min_gene_probs: List[float],
):
    """
    Full pipeline: clustering → merging → GWM building → evaluation.
    
    All parameters MUST be provided - no defaults.
    Parameters come from CLI arguments in pipeline.py.
    """
    
    # Create parameter combinations
    merge_params = list(product(
        merge_similarity_thresholds,
        merge_gene_thresholds,
        merge_metrics
    ))
    
    gwm_base_params = list(product(
        gwm_min_matches,
        gwm_min_core_genes,
        gwm_core_thresholds,
    ))
    
    logger.info(f"Testing {len(clustering_configs)} clustering × "
               f"{len(merge_params)} merge × {len(gwm_base_params)} GWM = "
               f"{len(clustering_configs) * len(merge_params) * len(gwm_base_params)} total configurations")

    # Step 1: Clustering
    clustering_results = run_clustering(
        clustering_configs,
        stat_matches_filepath,
        top_matches_filepath,
        out_dirpath,
        n_jobs
    )

    # Step 2: Merging
    merge_results = run_merging(
        clustering_results,
        merge_params,
        out_dirpath
    )

    # Step 3: GWM building
    gwm_results = run_gwm_building(
        merge_results,
        gwm_base_params,
        clusters_filepath,
        out_dirpath,
    )

    # Step 4: Evaluation
    best_result = run_evaluation(
        gwm_results,
        annotated_subclusters_filepath,
        reference_clusters_filepath,
        overlap_penalty_alpha,
        overlap_penalty_beta,
        out_dirpath
    )

    return best_result


def create_and_evaluate_gwms_from_motifs(
    motifs_dirpath: Path,
    clusters_filepath: Path,
    annotated_subclusters_filepath: Path,
    reference_clusters_filepath: Path,
    overlap_penalty_alpha: float,
    overlap_penalty_beta: float,
    out_dirpath: Path,
    merge_similarity_thresholds: List[float],
    merge_gene_thresholds: List[float],
    merge_metrics: List[str],
    gwm_min_matches: List[int],
    gwm_min_core_genes: List[int],
    gwm_core_thresholds: List[float],
    gwm_min_gene_probs: List[float],
):
    """
    Partial pipeline: start from existing motifs → merging → GWM building → evaluation.
    
    Skips clustering step and loads from pickled motif files.
    
    Args:
        motifs_dirpath: Directory containing config subdirs with motifs_*.pkl files
        ... (other args same as full pipeline)
    """
    
    # Create parameter combinations
    merge_params = list(product(
        merge_similarity_thresholds,
        merge_gene_thresholds,
        merge_metrics
    ))
    
    gwm_base_params = list(product(
        gwm_min_matches,
        gwm_min_core_genes,
        gwm_core_thresholds,
    ))
    
    logger.info(f"Testing {len(merge_params)} merge × {len(gwm_base_params)} GWM = "
               f"{len(merge_params) * len(gwm_base_params)} configurations from existing motifs")
    
    # Load existing clustering results from motifs
    clustering_results = {}
    
    for config_dir in sorted(motifs_dirpath.iterdir()):
        if not config_dir.is_dir():
            continue
        
        config_name = config_dir.name
        motifs_pkl = config_dir / f"motifs_{config_name}.pkl"
        
        if not motifs_pkl.exists():
            logger.warning(f"No motifs file at {motifs_pkl}, skipping {config_name}")
            continue
        
        # Load motifs
        with open(motifs_pkl, "rb") as f:
            motifs = pickle.load(f)
        
        logger.info(f"Loaded {len(motifs)} motifs for {config_name}")
        
        # Load metadata if available
        metadata = {}
        metadata_file = config_dir / f"metadata_{config_name}.json"
        if metadata_file.exists():
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
        
        # Analyze motifs
        motifs_stats = analyze_motifs(motifs)
        
        # Store results
        clustering_results[config_name] = {
            'config': {'name': config_name},
            'metadata': metadata,
            'n_motifs': len(motifs),
            'motifs_stats': motifs_stats,
            'motifs_filepath': str(motifs_pkl),
            'config_dir': str(config_dir),
        }
    
    if not clustering_results:
        raise ValueError(f"No valid motifs found in {motifs_dirpath}. "
                        f"Expected: {motifs_dirpath}/<config_name>/motifs_<config_name>.pkl")
    
    # Save loaded clustering results
    results_file = out_dirpath / "results_clustering_loaded.json"
    with open(results_file, 'w') as f:
        json.dump(convert_numpy_types(clustering_results), f, indent=2)
    logger.info(f"Saved loaded clustering results to {results_file}")
    
    # Continue with standard pipeline steps
    merge_results = run_merging(
        clustering_results,
        merge_params,
        out_dirpath
    )
    
    gwm_results = run_gwm_building(
        merge_results,
        gwm_base_params,
        clusters_filepath,
        out_dirpath,
    )
    
    best_result = run_evaluation(
        gwm_results,
        annotated_subclusters_filepath,
        reference_clusters_filepath,
        overlap_penalty_alpha,
        overlap_penalty_beta,
        out_dirpath
    )
    
    return best_result