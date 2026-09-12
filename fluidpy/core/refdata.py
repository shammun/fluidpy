"""Locating reference data — benchmark tables, cached downloads, manual files (seed file).

Resolution order for ``load(chapter, name)``:

1. ``reference/<chapter>/<name>``          — committed, cited benchmark data (Tier 2 of the data-and-benchmarks skill)
2. ``data/online/<chapter>/<name>``        — previously downloaded public dataset (git-ignored cache)
3. the registered URL for that (chapter, name), downloaded into ``data/online/<chapter>/`` now
4. ``data/manual/<chapter>/<name>``        — a file the reader downloaded by hand
5. ``/content/drive/MyDrive/fluidpy/data/manual/<chapter>/<name>`` on Colab

If nothing is found it raises with the exact instructions for getting the file — never a silent fallback, never a
fabricated array.
"""
from __future__ import annotations

from pathlib import Path

REGISTRY: dict[tuple[str, str], dict[str, str]] = {
    # ("ch07", "ghia_re100.csv"): {
    #     "url": "https://…",
    #     "credit": "Ghia, Ghia & Shin, J. Comput. Phys. 48, 387 (1982), table I",
    #     "licence": "reproduced values, cited",
    #     "verified": "YYYY-MM-DD",
    # },
}

DRIVE = Path("/content/drive/MyDrive/fluidpy/data/manual")


def _root() -> Path:
    here = Path.cwd().resolve()
    for cand in (here, *here.parents):
        if (cand / "fluidpy").is_dir():
            return cand
    return here


def fetch(url: str, dest: Path, timeout: int = 60) -> Path:
    """Download ``url`` to ``dest`` once; return the cached path on later calls."""
    import requests

    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    dest.write_bytes(r.content)
    return dest


def load(chapter: str, name: str) -> Path:
    """Return a local path to the reference file, downloading it if it is registered. Raises with instructions."""
    root = _root()
    candidates = [
        root / "reference" / chapter / name,
        root / "data" / "online" / chapter / name,
        root / "data" / "manual" / chapter / name,
        DRIVE / chapter / name,
    ]
    for c in candidates:
        if c.is_file():
            return c
    entry = REGISTRY.get((chapter, name.lower()))
    if entry:
        return fetch(entry["url"], root / "data" / "online" / chapter / name)
    raise FileNotFoundError(
        f"reference data '{name}' for {chapter} not found.\n"
        f"Looked in: {[str(c) for c in candidates]}\n"
        f"Either add it to reference/{chapter}/ via reference/{chapter}/make_refs.py (with its citation in "
        f"reference/{chapter}/SOURCES.md), or register a public URL in fluidpy/core/refdata.py REGISTRY."
    )
