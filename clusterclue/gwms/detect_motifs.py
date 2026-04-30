import logging
from collections import defaultdict
from clusterclue.classes.hits import MotifHit
from clusterclue.classes.subcluster_motif import SubclusterMotif
from pathlib import Path
from importlib.resources import files
from multiprocessing import Pool
from functools import partial


logger = logging.getLogger(__name__)


def parse_clusters_file(clusters_file):
    with open(clusters_file, "r") as infile:
        clusters = [read_cluster(line) for line in infile]
        clusters = {bgc_id: bgc_genes for bgc_id, bgc_genes in clusters}
    return clusters


def read_cluster(line):
    bgc_id, tokenized_genes = line.rstrip().split(",", 1)
    tokenized_genes = set(tokenized_genes.split(","))
    tokenized_genes.discard("-")  # genes without biosynthetic domains
    return bgc_id, tokenized_genes


def parse_motifs_file(motifs_file):
    subcluster_motifs = dict()
    infile = open(motifs_file, "r")
    while True:
        # read 4 lines at a time
        lines = [infile.readline().rstrip() for _ in range(4)]
        # stop it end of file
        if not lines[0]:
            break
        # add subcluster motif
        motif = SubclusterMotif.from_lines(lines)
        subcluster_motifs[motif.motif_id] = motif
    infile.close()
    return subcluster_motifs


def write_motif_hits(motif_hits, motifs, output_filepath):
    with open(output_filepath, "w") as outfile:
        # print header
        header_fields = [
            "bgc_id",
            "motif_id",
            "n_training_matches",
            "score_threshold",
            "score",
            "hit_genes",
        ]
        print("\t".join(header_fields), file=outfile)
        # print lines
        for bgc_id, hits in motif_hits.items():
            for motif_hit in hits:
                n_training_matches = motifs[motif_hit.motif_id].n_matches
                threshold = motifs[motif_hit.motif_id].threshold
                line_fields = [
                    motif_hit.bgc_id,
                    motif_hit.motif_id,
                    str(n_training_matches),
                    str(threshold),
                    str(motif_hit.score),
                    ",".join(sorted(motif_hit.tokenized_genes)),
                ]
                print("\t".join(line_fields), file=outfile)


def _process_bgc_motifs(bgc_data, motifs):
    """Worker function to detect motifs in a single BGC.
    
    Args:
        bgc_data: Tuple of (bgc_id, bgc_genes)
        motifs: Dictionary of motif objects
    
    Returns:
        Tuple of (bgc_id, list of MotifHit objects)
    """
    bgc_id, bgc_genes = bgc_data
    hits = []
    
    for motif in motifs.values():
        score = motif.calculate_score(bgc_genes)
        if score < motif.threshold:
            continue
        
        hit_genes = set(motif.tokenized_genes) & set(bgc_genes)
        if len(hit_genes) < 2:
            continue
        
        hits.append(MotifHit(bgc_id, motif.motif_id, score, hit_genes))
    
    return bgc_id, hits


def detect_motifs(clusters, motifs, n_jobs=1):
    """Detect motifs in clusters with optional parallelization.
    
    Args:
        clusters: Dictionary mapping bgc_id to sets of genes
        motifs: Dictionary of SubclusterMotif objects
        n_jobs: Number of parallel processes (1 = sequential)
    
    Returns:
        Dictionary mapping bgc_id to list of MotifHit objects
    """
    if n_jobs == 1:
        # Sequential processing (original behavior)
        motif_hits = defaultdict(list)
        for bgc_id, bgc_genes in clusters.items():
            for motif in motifs.values():
                score = motif.calculate_score(bgc_genes)
                if score < motif.threshold:
                    continue

                hit_genes = set(motif.tokenized_genes) & set(bgc_genes)
                if len(hit_genes) < 2:
                    continue

                motif_hits[bgc_id].append(
                    MotifHit(bgc_id, motif.motif_id, score, hit_genes)
                )
    else:
        # Parallel processing
        logger.info(f"Detecting motifs using {n_jobs} processes")
        worker_func = partial(_process_bgc_motifs, motifs=motifs)
        
        with Pool(processes=n_jobs) as pool:
            results = pool.map(worker_func, clusters.items())
        
        # Collect results
        motif_hits = defaultdict(list)
        for bgc_id, hits in results:
            if hits:
                motif_hits[bgc_id] = hits

    return motif_hits


def detect_gwms_in_clusters(
    clusters_filepath, 
    motifs_filepath, 
    output_filepath=None,
    n_jobs=1,
    ):
    clusters = parse_clusters_file(clusters_filepath)
    logger.info(f"Parsed {len(clusters)} clusters from {clusters_filepath}")
    motifs = parse_motifs_file(motifs_filepath)
    logger.info(f"Parsed {len(motifs)} motifs from {motifs_filepath}")
    motif_hits = detect_motifs(clusters, motifs, n_jobs=n_jobs)
    logger.info(f"Detected {sum(len(hits) for hits in motif_hits.values())} motif hits across {len(motif_hits)} clusters")

    if output_filepath is not None:
        write_motif_hits(motif_hits, motifs, output_filepath)
        logger.info(f"Wrote motif hits to {output_filepath}")

    return motif_hits
    

def visualise_gwm_hits(
    motif_gwms_filepath: str | Path,
    motif_hits_filepath: str | Path,
    genbank_dirpath: str | Path,
    domain_hits_filepath: str | Path,
    compound_structures_filepath: str | Path | None,
    output_dirpath: str | Path,
    n_jobs: int = 1,
):
    """Generate comprehensive HTML reports for BGCs and motifs with a master index.
    
    Args:
        motif_gwms_filepath: Path to motif GWMs file
        motif_hits_filepath: Path to motif hits file
        genbank_dirpath: Path to directory containing GenBank files
        domain_hits_filepath: Path to domain hits file
        compound_structures_filepath: Path to compound structures file (optional)
        output_dirpath: Path to output directory for HTML reports
        n_jobs: Number of parallel processes (1 = sequential, default: 1)
    """
    import subsketch as subsk
    
    session = subsk.SubSketchSession(
        motifs_file=motif_gwms_filepath,
        genbank_dir=genbank_dirpath,
        domain_hits_file=domain_hits_filepath,
        motif_hits_file=motif_hits_filepath,
        compounds_file=compound_structures_filepath,
        domain_colors_file=Path(files("clusterclue").joinpath("data").joinpath("domain_colors.txt"))
    )

    logger.info("Loading SubSketch session...")
    session.load(n_jobs=n_jobs, load_bgcs_upfront=(n_jobs > 1))
    
    # Generate all reports with master index in one call
    num_bgcs = len(session.list_genbanks())
    num_motifs = len(session.data.motifs)
    logger.info(f"Generating comprehensive report with {num_bgcs} BGCs and {num_motifs} motifs")
    
    session.generate_all_reports_with_master_index(
        output_dir=output_dirpath,
        gene_arrow_scaling=60,
        include_compound_plots=True,
        include_motif_plots=True,
        n_jobs=n_jobs,
    )
    
    logger.info(f"Reports generated successfully. Open {Path(output_dirpath) / 'index.html'} to view.")