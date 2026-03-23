import pandas as pd
from pathlib import Path
from src.utils.config import UHGG_EGGNOG_FOLDER, FINAL_LABELS_FILTERED

final_labels = pd.read_csv(FINAL_LABELS_FILTERED)
expected_species = set(final_labels['genome_id'].tolist())

downloaded = [d.name for d in Path(UHGG_EGGNOG_FOLDER).iterdir() if d.is_dir()]
has_eggnog = []
for species_id in downloaded:
    eggnog_file = Path(UHGG_EGGNOG_FOLDER) / species_id / f"{species_id}_eggNOG.tsv"
    if eggnog_file.exists() and eggnog_file.stat().st_size > 0:
        has_eggnog.append(species_id)

missing = expected_species - set(has_eggnog)

with open('missing_species.txt', 'w') as f:
    f.write('\n'.join(sorted(missing)))

print(f"Wrote {len(missing)} missing species to missing_species.txt")