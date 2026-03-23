#!/bin/bash
# setup_data.sh - Download large data files not tracked in git.
# Run once after cloning on a new machine.
#
# Usage: bash scripts/setup_data.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
UHGG_DIR="$PROJECT_ROOT/data/raw/uhgg"
UHGG_FTP="ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0"

mkdir -p "$UHGG_DIR"

# ── UHGG reference files ────────────────────────────────────────────────────

echo "Downloading genomes-all_metadata.tsv..."
if [ ! -f "$UHGG_DIR/genomes-all_metadata.tsv" ]; then
    wget -q -O "$UHGG_DIR/genomes-all_metadata.tsv" \
        "$UHGG_FTP/genomes-all_metadata.tsv"
    echo "  done."
else
    echo "  already exists, skipping."
fi

echo "Downloading KEGG completeness DB (53MB)..."
if [ ! -d "$UHGG_DIR/kegg_completeness_DB" ] || [ -z "$(ls -A "$UHGG_DIR/kegg_completeness_DB")" ]; then
    mkdir -p "$UHGG_DIR/kegg_completeness_DB"
    # Each species has a folder with a KEGG completeness file
    # Re-run the pipeline notebook 01_uhgg_labels.ipynb which reads these,
    # or use the UHGG FTP directly:
    echo "  NOTE: kegg_completeness_DB must be downloaded manually from UHGG FTP."
    echo "  See: $UHGG_FTP/species_catalogue/<prefix>/<species>/genome/<species>_kegg.tsv"
    echo "  Or run notebooks/pipeline/01_uhgg_labels.ipynb after downloading metadata."
else
    echo "  already exists, skipping."
fi

# ── eggNOG annotations (~4.2GB) ────────────────────────────────────────────

echo ""
echo "Downloading eggNOG annotations (~4.2GB, 4744 species)..."
echo "  This uses the existing download script."
if [ ! -f "$PROJECT_ROOT/data/processed/species_ids.txt" ]; then
    echo "  ERROR: data/processed/species_ids.txt not found. Make sure to pull from git first."
    exit 1
fi

EGGNOG_DIR="$UHGG_DIR/eggnog_annotations"
mkdir -p "$EGGNOG_DIR"

FTP_BASE="$UHGG_FTP/species_catalogue"
TOTAL=$(wc -l < "$PROJECT_ROOT/data/processed/species_ids.txt")
COUNT=0

while read species_id; do
    COUNT=$((COUNT + 1))
    OUT_FILE="$EGGNOG_DIR/${species_id}/${species_id}_eggNOG.tsv"
    if [ -f "$OUT_FILE" ]; then
        continue
    fi
    prefix="${species_id:0:11}"
    mkdir -p "$EGGNOG_DIR/$species_id"
    wget -q -O "$OUT_FILE" \
        "${FTP_BASE}/${prefix}/${species_id}/genome/${species_id}_eggNOG.tsv" || \
        echo "  WARNING: failed to download $species_id"
    sleep 0.1
    if (( COUNT % 100 == 0 )); then
        echo "  $COUNT / $TOTAL species downloaded..."
    fi
done < "$PROJECT_ROOT/data/processed/species_ids.txt"

echo "  eggNOG annotations complete."

echo ""
echo "Setup complete. Missing files summary:"
echo "  genomes-all_metadata.tsv : $([ -f "$UHGG_DIR/genomes-all_metadata.tsv" ] && echo 'OK' || echo 'MISSING')"
echo "  bac120_iqtree.nwk        : $([ -f "$UHGG_DIR/bac120_iqtree.nwk" ] && echo 'OK' || echo 'MISSING')"
echo "  eggnog_annotations/      : $([ -d "$UHGG_DIR/eggnog_annotations" ] && echo 'exists (may be partial)' || echo 'MISSING')"
echo "  kegg_completeness_DB/    : $([ -d "$UHGG_DIR/kegg_completeness_DB" ] && echo 'exists (may be partial)' || echo 'MISSING')"
