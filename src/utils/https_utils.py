import requests
import time

session = requests.Session()
session.headers.update(
    {"User-Agent": "phylo-gnn/1.0"}
)

_last_call = 0


def sleep_if_needed(min_interval: float = 0.25):
    global _last_call
    now = time.monotonic()
    elapsed = now - _last_call
    if elapsed < min_interval:
        time.sleep(min_interval - elapsed)
    _last_call = time.monotonic()
    return


def get_text(url: str, min_interval: float = 0.25) -> str:
    sleep_if_needed(min_interval=min_interval)
    response = session.get(url, timeout=20)
    response.raise_for_status()
    return response.text
