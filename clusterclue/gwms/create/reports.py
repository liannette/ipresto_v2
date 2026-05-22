from pathlib import Path


def generate_evaluation_report(results: dict, output_file: Path):
    """Generate human-readable evaluation report."""
    with open(output_file, 'w') as f:
        f.write("="*140 + "\n")
        f.write("GWM EVALUATION REPORT\n")
        f.write("="*140 + "\n\n")
        
        f.write("This report compares different configurations based on their ability to\n")
        f.write("produce biologically valid Gene Weight Modules (GWMs) that match known subclusters.\n\n")
        
        # Summary table
        f.write("SUMMARY TABLE\n")
        f.write("-"*140 + "\n")
        f.write(f"{'Config':<40} {'Method':<10} {'Initial':<8} {'Merged':<8} {'GWMs':<8} "
                f"{'Overlap':<10} {'MRPOS':<10}\n")
        f.write(f"{'':40} {'':10} {'Motifs':8} {'Motifs':8} {'':8} {'Score':10} {'Score':10}\n")
        f.write("-"*140 + "\n")
        
        # Sort by MRPOS
        sorted_results = sorted(
            [(name, r) for name, r in results.items() if 'error' not in r],
            key=lambda x: x[1]['mean_penalized_score'],
            reverse=True
        )
        
        for name, result in sorted_results:
            config = result.get('config', {})
            if 'method' in config:
                method = config['method'].upper()
            else:
                config_name = config.get('name', name)
                if config_name.startswith('hdb'):
                    method = 'HDBSCAN'
                elif config_name.startswith('kmeans'):
                    method = 'KMEANS'
                else:
                    method = 'UNKNOWN'
            
            n_initial = result['n_initial_motifs']
            n_merged = result['n_merged_motifs']
            n_gwms = result['n_gwms']
            overlap = result['mean_overlap_score']
            penalized = result['mean_penalized_score']
            
            f.write(f"{name:<40} {method:<10} {n_initial:<8} {n_merged:<8} {n_gwms:<8} "
                    f"{overlap:<10.4f} {penalized:<10.4f}\n")
        
        f.write("-"*140 + "\n\n")
    
        # Winner section
        winner_name, winner = sorted_results[0]
        f.write("="*140 + "\n")
        f.write("BEST CONFIGURATION\n")
        f.write("="*140 + "\n\n")
        f.write(f"Winner: {winner_name}\n\n")
        
        config = winner.get('config', {})
        config_name = config.get('name', winner_name)
        
        f.write("Clustering Configuration:\n")
        
        # Determine method
        if 'method' in config:
            method = config['method']
        elif config_name.startswith('hdb'):
            method = 'hdbscan'
        elif config_name.startswith('kmeans'):
            method = 'kmeans'
        else:
            method = 'unknown'
        
        f.write(f"  Method: {method.upper()}\n")
        
        # Parse config for parameters
        import re
        if method == 'kmeans':
            if 'k' in config:
                f.write(f"  K: {config['k']}\n")
            else:
                k_match = re.search(r'_k(\d+)', config_name)
                if k_match:
                    f.write(f"  K: {k_match.group(1)}\n")
        else:  # hdbscan
            if 'min_cluster_size' in config:
                f.write(f"  Min cluster size: {config['min_cluster_size']}\n")
            else:
                mcs_match = re.search(r'_mcs(\d+)', config_name)
                if mcs_match:
                    f.write(f"  Min cluster size: {mcs_match.group(1)}\n")
            
            if 'cluster_selection_epsilon' in config:
                f.write(f"  Cluster selection epsilon: {config['cluster_selection_epsilon']}\n")
            else:
                eps_match = re.search(r'_eps(\d+)', config_name)
                if eps_match:
                    eps_val = int(eps_match.group(1)) / 100.0
                    f.write(f"  Cluster selection epsilon: {eps_val}\n")
        
        # SVD info
        use_svd = config.get('use_svd', 'raw' not in config_name)
        f.write(f"  SVD: {use_svd}\n")
        
        if use_svd:
            if 'target_variance' in config:
                f.write(f"  Target variance: {config['target_variance']}\n")
            else:
                svd_match = re.search(r'_svd(\d+)', config_name)
                if svd_match:
                    tv_val = int(svd_match.group(1)) / 100.0
                    f.write(f"  Target variance: {tv_val}\n")
        
        f.write("\nMerge Parameters:\n")
        merge_params = winner.get('merge_params', {})
        f.write(f"  Similarity threshold: {merge_params.get('similarity_threshold', 'N/A')}\n")
        f.write(f"  Gene prob threshold: {merge_params.get('gene_threshold', 'N/A')}\n")
        f.write(f"  Similarity metric: {merge_params.get('metric', 'N/A')}\n")
        
        f.write("\nGWM Hyperparameters:\n")
        params = winner['gwm_hyperparameter']
        f.write(f"  Min matches: {params['min_matches']}\n")
        f.write(f"  Min core genes: {params['min_core_genes']}\n")
        f.write(f"  Core threshold: {params['core_threshold']}\n")
        f.write(f"  Min gene probability: {params['min_gene_prob']}\n")
        
        f.write("\nPipeline Results:\n")
        f.write(f"  Initial motifs: {winner['n_initial_motifs']}\n")
        f.write(f"  Merged motifs: {winner['n_merged_motifs']} "
                f"({(winner['n_initial_motifs']-winner['n_merged_motifs'])/winner['n_initial_motifs']*100:.1f}% reduction)\n")
        f.write(f"  Final GWMs: {winner['n_gwms']} "
                f"({winner['n_gwms']/winner['n_merged_motifs']*100:.1f}% build success)\n")
        
        f.write("\nBiological Validation:\n")
        f.write(f"  Mean overlap score: {winner['mean_overlap_score']:.4f}\n")
        f.write(f"  Mean MRPOS (penalized): {winner['mean_penalized_score']:.4f}\n")
        
        # Cluster quality metrics
        metadata = winner.get('metadata', {})
        if metadata:
            f.write("\nCluster Quality Metrics:\n")
            if 'silhouette' in metadata:
                f.write(f"  Silhouette score: {metadata['silhouette']:.4f}\n")
            if 'davies_bouldin' in metadata:
                f.write(f"  Davies-Bouldin score: {metadata['davies_bouldin']:.4f}\n")
            if 'calinski_harabasz' in metadata:
                f.write(f"  Calinski-Harabasz score: {metadata['calinski_harabasz']:.2f}\n")
            if 'noise_fraction' in metadata:
                f.write(f"  Noise fraction: {metadata['noise_fraction']:.2%}\n")
        
        f.write("\n")
        
        f.write("="*140 + "\n")
        f.write("\nKEY METRICS EXPLAINED\n")
        f.write("-"*140 + "\n")
        f.write("Overlap Score:      How well GWMs match reference subclusters (F1 score)\n")
        f.write("MRPOS Score:        Mean Redundancy Penalized Overlap Score - overlap with penalty for cluster size\n")
        f.write("Silhouette Score:   Cluster separation quality (-1 to 1, higher is better)\n")
        f.write("Davies-Bouldin:     Cluster compactness (lower is better)\n")
        f.write("Calinski-Harabasz:  Cluster definition quality (higher is better)\n")
        f.write("Noise Fraction:     Proportion of modules classified as noise (HDBSCAN only)\n")
        f.write("="*140 + "\n")