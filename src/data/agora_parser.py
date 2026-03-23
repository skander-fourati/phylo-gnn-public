from typing import Optional
import re


def clean_species_from_agora(full_name: str) -> Optional[str]:
    """
    Clean AGORA2 species names.

    Examples:
        "Escherichia coli K-12" → "Escherichia coli"
        "Bacteroides fragilis NCTC 9343" → "Bacteroides fragilis"
        "Alistipes sp. CAG:831" → "Alistipes sp."
    """
    if not isinstance(full_name, str) or not full_name.strip():
        return None

    # Normalize whitespace
    full_name = re.sub(r"\s+", " ", full_name).strip()

    parts = full_name.split()
    if not parts:
        return None

    # Handle "Genus sp.*" → "Genus sp."
    if len(parts) >= 2 and parts[1].lower().startswith("sp"):
        return f"{parts[0]} sp."

    # Take genus + species only (first 2 words)
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"

    return parts[0]