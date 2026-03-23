#!/bin/bash
# download_missing_eggnog.sh

# Use HTTPS instead of FTP
BASE_URL="https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/species_catalogue"

while read species_id; do
    prefix="${species_id:0:11}"

    output_dir="/Users/encordsf/Desktop/data/bio_research_data/uhgg/eggnog_annotations/${species_id}"
    mkdir -p "$output_dir"

    wget -q "${BASE_URL}/${prefix}/${species_id}/genome/${species_id}_eggNOG.tsv" \
        -O "${output_dir}/${species_id}_eggNOG.tsv"

    if [ $? -eq 0 ]; then
        echo "✓ Downloaded ${species_id}"
    else
        echo "✗ Failed ${species_id}"
    fi

    sleep 0.1
done < missing_species.txt