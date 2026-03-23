"""Project paths.

Small reference files and processed outputs live inside the project under data/.
Bulk downloaded data (eggnog annotations, gff files, KEGG completeness DB, KEGG cache)
lives outside the project at BULK_DATA_DIR to avoid IDE indexing issues.

Set PHYLO_BULK_DATA env variable to override the default bulk data location.
"""
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"
RUNS_DIR = RESULTS_DIR / "runs"
CHECKPOINTS_DIR = RESULTS_DIR / "checkpoints"
MODEL_CONFIGS_DIR = PROJECT_ROOT / "model_configs"

# Bulk data lives outside the project (large files, not committed)
BULK_DATA_DIR = Path(os.environ.get("PHYLO_BULK_DATA", str(Path.home() / "Desktop" / "data" / "bio_research_data")))
UHGG_BULK_DIR = BULK_DATA_DIR / "uhgg"
KEGG_BULK_DIR = BULK_DATA_DIR / "kegg"

# Small reference files (committed, inside project)
UHGG_FOLDER = RAW_DIR / "uhgg"
KEGG_FOLDER = RAW_DIR / "kegg"
FIGURES_DIR = RAW_DIR / "figures"

# KEGG API cache (large, outside project)
MODULES_DIR = KEGG_BULK_DIR / "modules_by_org"
MODULE_ENTRY_DIR = KEGG_FOLDER / "module_entries"
GENOMES_DIR = KEGG_BULK_DIR / "genome_entries"

# Key data files
RAW_AGORA_METADATA = RAW_DIR / "agora2" / "metadata.xlsx"
CLEAN_AGORA_METADATA = PROCESSED_DIR / "agora2" / "metadata.csv"
UHGG_METADATA = DATA_DIR / "genomes-all_metadata.tsv"
KEGG_KOs = PROCESSED_DIR / "kegg" / "kegg_kos.json"
AGORA_LABELS = RAW_DIR / "agora2" / "metadata.csv"
UHGG_FILTERED_AGORA_MATCHED = PROCESSED_DIR / "uhgg_filtered_agora_matched.csv"
COMPLETENESS_DF = PROCESSED_DIR / "completeness_df.csv"
FINAL_LABELS = PROCESSED_DIR / "metabolite_labels.csv"
FINAL_LABELS_FILTERED = PROCESSED_DIR / "metabolite_labels_filtered.csv"

# Large data outside project
UHGG_KEGG_COMPLETENESS = UHGG_BULK_DIR / "kegg_completeness_DB"
UHGG_GFF_DIR = UHGG_BULK_DIR / "gff_files"
UHGG_EGGNOG_FOLDER = UHGG_BULK_DIR / "eggnog_annotations"

PHYLO_TREE = UHGG_FOLDER / "bac120_iqtree.nwk"
PHYLO_DIST_MATRIX = PROCESSED_DIR / "phylo_distance_matrix.npy"

ALL_RM_LABELS = PROCESSED_DIR / "all_rm_labels"
LABELS_RM_METABOLITES = PROCESSED_DIR / "labels.csv"

SIMPLE_GNN_XP_RESULTS = RUNS_DIR / "SimpleGNN"
