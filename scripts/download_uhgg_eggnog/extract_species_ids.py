import pandas as pd
from src.utils.config import FINAL_LABELS_FILTERED

df = pd.read_csv(FINAL_LABELS_FILTERED)
species_ids = df['genome_id'].tolist()

# Save to file for bash script
with open('../../data/processed/species_ids.txt', 'w') as f:
    f.write('\n'.join(species_ids))