from rapidfuzz import fuzz, process

def best_match(lookup: str, reference_list: list):
    if not isinstance(lookup, str) or not lookup.strip():
        return None, None, 0.0
    hit, score, _ = process.extractOne(lookup, reference_list, scorer=fuzz.token_sort_ratio)
    return lookup, hit, float(score)