import time, random, requests
from pathlib import Path
from src.utils.https_utils import get_text
from src.utils.config import MODULES_DIR, MODULE_ENTRY_DIR, GENOMES_DIR

KEGG_API_URL = "https://rest.kegg.jp"


def _fetch_with_cache_and_retry(url: str, cache_path: Path, retries: int = 3, min_interval: float = 0.25) -> tuple[
    str, str | None]:
    """Generic fetch with caching and retry logic."""

    # Create directory only when actually needed ✅
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    if cache_path.exists():
        return cache_path.read_text(encoding="utf-8"), None

    backoff = 0.5
    for attempt in range(1, retries + 1):
        try:
            text = get_text(url, min_interval=min_interval)
            if not text.strip():
                return "", "EMPTY"

            cache_path.write_text(text, encoding="utf-8")
            return text, None

        except requests.HTTPError as e:
            code = e.response.status_code if e.response else None
            if code == 403:
                return "", "HTTP_403"
            if not (500 <= (code or 0) < 600) and attempt == retries:
                return "", f"HTTP_{code or 'ERROR'}"

        except requests.Timeout:
            if attempt == retries:
                return "", "TIMEOUT"

        except requests.RequestException:
            if attempt == retries:
                return "", "NETWORK"

        time.sleep(backoff + random.uniform(0, 0.25))
        backoff *= 2

    return "", "MAX_RETRIES"


def fetch_modules_for_org(org_code: str, retries: int = 3, min_interval: float = 0.25) -> tuple[str, str | None]:
    cache_path = MODULES_DIR / f"{org_code}.txt"
    url = f"{KEGG_API_URL}/link/module/{org_code}"
    return _fetch_with_cache_and_retry(url, cache_path, retries, min_interval)


def fetch_genome_entry(t_id: str, retries: int = 3, min_interval: float = 0.25) -> tuple[str, str | None]:
    cache_path = GENOMES_DIR / f"{t_id}.txt"
    url = f"{KEGG_API_URL}/get/gn:{t_id}"
    return _fetch_with_cache_and_retry(url, cache_path, retries, min_interval)


def fetch_module_entry(module_id: str, retries: int = 3, min_interval: float = 0.25) -> tuple[str, str | None]:
    cache_path = MODULE_ENTRY_DIR / f"{module_id}.txt"
    url = f"{KEGG_API_URL}/get/{module_id}"
    return _fetch_with_cache_and_retry(url, cache_path, retries, min_interval)


def parse_module_definition(module_txt: str) -> list[set[str]]:
    """
    Parse DEFINITION field into steps.
    DEFINITION can span multiple lines (continuation lines are indented).
    """
    # Find and concatenate all DEFINITION lines
    definition_lines = []
    in_definition = False

    for line in module_txt.strip().splitlines():
        if line.startswith("DEFINITION"):
            in_definition = True
            definition_lines.append(line.replace("DEFINITION", "").strip())
        elif in_definition:
            # Continuation line (starts with whitespace)
            if line and line[0].isspace():
                definition_lines.append(line.strip())
            else:
                # Hit next section, stop
                break

    # Join all definition lines into one string
    definition = " ".join(definition_lines)

    if not definition:
        return []

    # Now parse the combined definition
    steps = []
    current_step = set()
    in_parens = 0
    current_ko = ""

    for char in definition:
        if char == '(':
            in_parens += 1
        elif char == ')':
            in_parens -= 1
            if current_ko.startswith('K'):
                current_step.add(current_ko)
            current_ko = ""
            if in_parens == 0 and current_step:
                steps.append(current_step)
                current_step = set()
        elif char == ',' and in_parens > 0:
            if current_ko.startswith('K'):
                current_step.add(current_ko)
            current_ko = ""
        elif char == ' ':
            if in_parens > 0:
                # Inside parentheses, space just separates (shouldn't happen)
                if current_ko.startswith('K'):
                    current_step.add(current_ko)
                current_ko = ""
            else:
                # Outside parentheses, space = end of step
                if current_ko.startswith('K'):
                    steps.append({current_ko})
                current_ko = ""
        elif char.isalnum():
            current_ko += char

    # Don't forget last KO if any
    if current_ko.startswith('K'):
        if in_parens > 0:
            current_step.add(current_ko)
        else:
            steps.append({current_ko})

    # Add final step if we were building one
    if current_step:
        steps.append(current_step)

    return steps


def calculate_pathway_completeness(module_steps: list[set[str]], ko_hits: set[str]) -> float:
    """
    Calculate what fraction of pathway steps are satisfied.

    Args:
        module_steps: List of sets, each set = alternative KOs for that step. Ex: [{'K01','K02'}, {'K03','K04'}]
        ko_hits: Set of KOs present in the genome

    Returns:
        Fraction of steps satisfied (0.0 to 1.0)
    """
    if not module_steps:
        return 0.0

    satisfied_steps = 0
    for step_alternatives in module_steps:
        # Step is satisfied if ANY alternative KO is present
        if any(ko in ko_hits for ko in step_alternatives):
            satisfied_steps += 1

    return satisfied_steps / len(module_steps)