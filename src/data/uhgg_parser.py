from typing import Optional, Literal
import re, time
from src.utils.config import UHGG_METADATA, UHGG_GFF_DIR,UHGG_KEGG_COMPLETENESS
from src.utils.https_utils import sleep_if_needed
import pandas as pd
from urllib.request import urlretrieve
from urllib.error import URLError, ContentTooShortError
from pathlib import Path

def _clean_name_token(token: str) -> Optional[str]:
    if not isinstance(token, str) or not token.strip():
        return None
    # remove GTDB-style rank prefixes if present (e.g., "s__")
    token = re.sub(r"^[dpcofgs]__", "", token)

    # drop parenthetical notes and anything after "="
    token = re.sub(r"\(.*?\)", "", token)
    token = re.sub(r"=.*", "", token)

    # normalize underscores/whitespace
    token = token.replace("_", " ")
    token = re.sub(r"\s+", " ", token).strip()

    parts = token.split()
    if not parts:
        return None

    # drop single-letter GTDB clade tag after genus (e.g., "Blautia A" -> "Blautia") ---
    if (
        len(parts) >= 2
        and len(parts[1]) == 1
        and parts[1].isalpha()
        and parts[1].isupper()
    ):
        parts.pop(1)

    # "Genus sp.*" → standardize to "Genus sp."
    if len(parts) >= 2 and parts[1].lower().startswith("sp"):
        return f"{parts[0]} sp."

    # default to Genus + species epithet
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"

    return parts[0]


def clean_species_from_uhgg_lineage(lineage: str) -> Optional[str]:
    """
    UHGG lineage like '...;s__Blautia_A faecis' → 'Blautia faecis'
    We extract the LAST semicolon-separated segment (the s__ entry), then clean it.
    """
    if not isinstance(lineage, str) or not lineage.strip():
        return None
    # split on semicolons, take the last non-empty segment
    parts = [p.strip() for p in lineage.split(";") if p.strip()]
    if not parts:
        return None
    last = parts[-1]  # s__...
    return _clean_name_token(last)


def download_uhgg_genomes(
        genome_ids: Optional[set[str]] = None,
        length_range: Optional[tuple[int, int]] = None,
        n_contigs_range: Optional[tuple[int, int]] = None,
        n50_range: Optional[tuple[int, int]] = None,
        gc_range: Optional[tuple[float, float]] = None,
        completeness_range: Optional[tuple[float, float]] = None,
        contamination_range: Optional[tuple[float, float]] = None,
        rrna_5s_range: Optional[tuple[int, int]] = None,
        rrna_16s_range: Optional[tuple[int, int]] = None,
        rrna_23s_range: Optional[tuple[int, int]] = None,
        trnas_range: Optional[tuple[int, int]] = None,
        ncrnas_range: Optional[tuple[int, int]] = None,
        genome_type: Optional[Literal["Isolate", "MAG"]] = None,
        taxonomic_lineage: Optional[str] = None
) -> list[str]:
    """Download UHGG genomes matching criteria."""

    if genome_type is not None and genome_type not in ("Isolate", "MAG"):
        raise ValueError(f"genome_type must be 'Isolate' or 'MAG', got: {genome_type}")

    # Load metadata
    df = pd.read_csv(UHGG_METADATA, sep='\t')

    # Apply filters
    mask = pd.Series([True] * len(df), index=df.index)

    if genome_ids is not None:
        mask &= df['Genome'].isin(genome_ids)

    if length_range is not None:
        mask &= (df['Length'] >= length_range[0]) & (df['Length'] <= length_range[1])

    if n_contigs_range is not None:
        mask &= (df['N_contigs'] >= n_contigs_range[0]) & (df['N_contigs'] <= n_contigs_range[1])

    if n50_range is not None:
        mask &= (df['N50'] >= n50_range[0]) & (df['N50'] <= n50_range[1])

    if gc_range is not None:
        mask &= (df['GC_content'] >= gc_range[0]) & (df['GC_content'] <= gc_range[1])

    if completeness_range is not None:
        mask &= (df['Completeness'] >= completeness_range[0]) & (df['Completeness'] <= completeness_range[1])

    if contamination_range is not None:
        mask &= (df['Contamination'] >= contamination_range[0]) & (df['Contamination'] <= contamination_range[1])

    if rrna_5s_range is not None:
        mask &= (df['rRNA_5S'] >= rrna_5s_range[0]) & (df['rRNA_5S'] <= rrna_5s_range[1])

    if rrna_16s_range is not None:
        mask &= (df['rRNA_16S'] >= rrna_16s_range[0]) & (df['rRNA_16S'] <= rrna_16s_range[1])

    if rrna_23s_range is not None:
        mask &= (df['rRNA_23S'] >= rrna_23s_range[0]) & (df['rRNA_23S'] <= rrna_23s_range[1])

    if trnas_range is not None:
        mask &= (df['tRNAs'] >= trnas_range[0]) & (df['tRNAs'] <= trnas_range[1])

    if ncrnas_range is not None:
        mask &= (df['ncRNAs'] >= ncrnas_range[0]) & (df['ncRNAs'] <= ncrnas_range[1])

    if genome_type is not None:
        mask &= df['Genome_type'] == genome_type

    if taxonomic_lineage is not None:
        mask &= df['Taxonomic_lineage'].str.contains(taxonomic_lineage, case=False, na=False)

    # Get filtered genomes
    filtered = df[mask]

    # Download each genome
    downloaded = []
    for idx, row in filtered.iterrows():
        genome_id = row['Genome']
        ftp_url = row['FTP_download']
        url = ftp_url.replace("ftp://", "https://")
        output_path = UHGG_GFF_DIR / f"{genome_id}.gff.gz"

        if output_path.exists():
            downloaded.append(genome_id)
            continue

        for attempt in range(5):
            try:
                sleep_if_needed(min_interval=0.25)
                urlretrieve(url, str(output_path))
                downloaded.append(genome_id)
                break
            except (URLError, ContentTooShortError, ConnectionError, TimeoutError):
                if attempt == 4:
                    raise
                time.sleep(2 * (attempt + 1))  # 2s,4s,6s,8s

    return downloaded


def list_species_reps(metadata_tsv_path: Path = UHGG_METADATA)-> set[str]:
    uhgg_metadata = pd.read_csv(metadata_tsv_path, sep="\t")
    return sorted(uhgg_metadata["Species_rep"].unique())


def parse_all_kegg_completeness(
        kegg_completeness_path: Path = UHGG_KEGG_COMPLETENESS,
        use_core: bool = False
) -> pd.DataFrame:
    """
    Parse all KEGG completeness files at once.

    Returns:
        DataFrame with shape (n_species, n_modules)
        Rows = genome IDs, Columns = module IDs, Values = completeness
    """
    completeness_col = 'core' if use_core else 'pangenome'

    all_data = []
    missing_files = []

    # Get all completeness files
    tsv_files = list(kegg_completeness_path.glob("*_clstr_kegg_comp.tsv"))

    for file_path in tsv_files:
        # Extract genome ID from filename
        genome_id = file_path.stem.replace('_clstr_kegg_comp', '')

        try:
            df = pd.read_csv(file_path, sep="\t")

            # Split module column
            df[['module_id', 'desc']] = df["#module"].str.split('|', expand=True)

            # Create dict for this genome
            comp_dict = dict(zip(df['module_id'], df[completeness_col]))
            comp_dict['genome_id'] = genome_id

            all_data.append(comp_dict)

        except Exception as e:
            missing_files.append((genome_id, str(e)))

    # Convert to wide DataFrame
    completeness_df = pd.DataFrame(all_data).set_index('genome_id')

    # Fill NaN with 0 (module not present)
    completeness_df = completeness_df.fillna(0.0)

    return completeness_df
