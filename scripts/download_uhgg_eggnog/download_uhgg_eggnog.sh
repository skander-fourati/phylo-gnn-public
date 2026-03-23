#!/bin/bash

FTP_BASE="ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/species_catalogue"

while read species_id; do
    prefix="${species_id:0:11}"
    wget -q -P "/Users/encordsf/Desktop/data/bio_research_data/uhgg/eggnog_annotations/${species_id}/" \
        "${FTP_BASE}/${prefix}/${species_id}/genome/${species_id}_eggNOG.tsv"
    sleep 0.1
done < /Users/encordsf/Desktop/phylo-gnn/data/processed/species_ids.txt