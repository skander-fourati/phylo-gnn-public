import os
from pathlib import Path
from src.utils.config import UHGG_EGGNOG_FOLDER, FINAL_LABELS_FILTERED
import pandas as pd

# 1. How many species do we expect?
final_labels = pd.read_csv(FINAL_LABELS_FILTERED)
expected_species = set(final_labels['genome_id'].tolist())
print(f"Expected species: {len(expected_species)}")

# 2. How many directories were created?
downloaded_dirs = [d.name for d in Path(UHGG_EGGNOG_FOLDER).iterdir() if d.is_dir()]
print(f"Downloaded directories: {len(downloaded_dirs)}")

# 3. How many actually have eggNOG files?
has_eggnog = []
for species_id in downloaded_dirs:
    eggnog_file = Path(UHGG_EGGNOG_FOLDER) / species_id / f"{species_id}_eggNOG.tsv"
    if eggnog_file.exists():
        has_eggnog.append(species_id)

print(f"With eggNOG files: {len(has_eggnog)}")

# 4. Which species are missing?
missing = expected_species - set(has_eggnog)
print(f"\nMissing: {len(missing)} species")
print(f"First 10 missing: {list(missing)[:10]}")

# 5. Check if any files are size 0 (failed downloads)
zero_size = []
for species_id in has_eggnog:
    eggnog_file = Path(UHGG_EGGNOG_FOLDER) / species_id / f"{species_id}_eggNOG.tsv"
    if eggnog_file.stat().st_size == 0:
        zero_size.append(species_id)

print(f"\nZero-size files: {len(zero_size)}")
print(f"First 10 zero-size: {list(zero_size)[:10]}")