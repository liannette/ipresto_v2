from pathlib import Path
from typing import Optional

from ipresto.presto_stat.utils import (
    read_clusterfile,
    read_stat_modules,
    write_detected_stat_modules,
    write_detected_stat_module_ids,
)
from ipresto.presto_stat.build import generate_stat_modules, filter_infrequent_modules
from ipresto.presto_stat.detect import detect_modules_in_bgcs
from ipresto.presto_stat.utils import write_stat_modules


class StatOrchestrator:
    def run(
        self,
        out_dir_path: str,
        cluster_file: str,
        stat_modules_file_path: Optional[str],
        min_gene_occurence: int,
        pval_cutoff: float,
        min_genes_per_bgc: int,
        cores: int,
        verbose: bool,
    ) -> str:
        """ """
        # Create output directory
        out_dir = Path(out_dir_path)
        out_dir_path.mkdir(parents=True, exist_ok=True)

        # Read clusters
        bgcs = read_clusterfile(cluster_file, min_genes_per_bgc, verbose)

        # Generate subcluster modules, if not provided
        if stat_modules_file_path is None:
            if verbose:
                print(
                    "\nComputing new PRESTO-STAT subcluster modules using the  "
                    f"clusters from {cluster_file}"
                )
            modules = generate_stat_modules(
                out_dir,
                bgcs,
                min_genes_per_bgc,
                pval_cutoff,
                min_gene_occurence,
                cores,
                verbose,
            )
            module_file_path = out_dir / "stat_modules_unfiltered.txt"
            write_stat_modules(modules, module_file_path)
            if verbose:
                print(
                    f"Unfiltered PRESTO-STAT subcluster modules have been "
                    f"saved to: {module_file_path}"
                )

            # remove modules that occur too infrequently in the bgc dataset
            min_occurence = 2
            if verbose:
                print(
                    f"\nRemoving modules that occur less than {min_occurence} "
                    "times in the BGC dataset."
                )
            filtered_modules = filter_infrequent_modules(modules, min_occurence)
            if verbose:
                percent_removed = (
                    (len(modules) - len(filtered_modules)) / len(modules) * 100
                )
                print("  removed {:.1f}% of modules".format(percent_removed))

            # write filtered modules to file
            module_file_path = out_dir / "stat_modules_filtered.txt"
            write_stat_modules(filtered_modules, module_file_path)
            if verbose:
                print(
                    f"Filtered PRESTO-STAT subcluster modules have been saved to: "
                    f"{module_file_path}"
                )

        # # Use subcluster modules, if provided
        # else:
        #     if verbose:
        #         print(
        #             "\nUsing PRESTO-STAT subcluster modules from: "
        #             f"{stat_modules_file_path}"
        #         )
        #     modules = read_stat_modules(stat_modules_file_path)

        # Detect modules in bgcs
        if verbose:
            print("\nDetecting PRESTO-STAT subcluster modules in the input clusters")
        detected_modules = detect_modules_in_bgcs(bgcs, modules, cores)

        if verbose:
            print("\nPRESTO-STAT subcluster detection completed.")

        # Write the detected modules IDs to a text file
        stat_subclusters_path = out_dir / "detected_stat_module_ids.txt"
        write_detected_stat_module_ids(detected_modules, stat_subclusters_path)
        if verbose:
            print(
                f"Detected STAT module IDs have been saved to: {stat_subclusters_path}"
            )

        # Write the detected modules to a fasta file
        stat_subclusters_path = out_dir / "detected_stat_modules.txt"
        write_detected_stat_modules(detected_modules, stat_subclusters_path)
        if verbose:
            print(f"Detected STAT modules have been saved to: {stat_subclusters_path}")

        return stat_subclusters_path
