"""Fail if anything book-derived is about to be (or has been) committed to the PUBLIC repository (CLAUDE.md rule 9).

Checks every tracked file (``git ls-files``; after ``git add`` that includes staged files), or the paths given:

* forbidden paths: ``*.pdf``, ``chapters/`` (except ``chapters/_page_map.json``: page numbers only), ``data/book/``,
  ``data/manual/``, ``outputs/``, ``reports/viz/``, ``reports/**/figures/``, ``notebooks/executed_*``,
  ``tests/book_values_*``;
* ``*_colab.ipynb`` must have **no outputs** (Colab's "Save a copy in GitHub" writes them);
* no notebook or HTML page may carry the private-run banner ``book values active``;
* explainers and pages must not reference rendered book pages (``chapters/pages``);
* any tracked file above 50 MB (GitHub refuses 100 MB; pages that large are also unusable on phones).

Usage::

    .venv/Scripts/python.exe tools/check_public.py            # all tracked files
    .venv/Scripts/python.exe tools/check_public.py a.ipynb b.html
    git config core.hooksPath .githooks                       # installs the pre-push hook that runs this

Exit code 0 = clean, 1 = violations (listed).
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN = re.compile(
    r"(\.pdf$|^chapters/(?!_page_map\.json$)|^data/book/|^data/manual/|^outputs/|^reports/viz/|^reports/.*/figures/|"
    r"^notebooks/executed_|^tests/book_values_)", re.I)
BANNER = re.compile(r"book values active", re.I)
PAGE_REF = re.compile(r"chapters/pages/")
MAX_MB = 50


def tracked_files() -> list[str]:
    out = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True).stdout
    return [line for line in out.splitlines() if line]


def check(paths: list[str]) -> list[str]:
    problems: list[str] = []
    for rel in paths:
        rel = rel.replace("\\", "/")
        if FORBIDDEN.search(rel):
            problems.append(f"{rel}: book-derived or local-only path must not be tracked")
            continue
        p = ROOT / rel
        if not p.is_file():
            continue
        size_mb = p.stat().st_size / 1e6
        if size_mb > MAX_MB:
            problems.append(f"{rel}: {size_mb:.0f} MB (> {MAX_MB} MB)")
        if p.suffix == ".ipynb":
            text = p.read_text(encoding="utf-8", errors="replace")
            try:
                nb = json.loads(text)
            except Exception as exc:  # noqa: BLE001
                problems.append(f"{rel}: not valid JSON ({exc})")
                continue
            n_out = sum(len(c.get("outputs", [])) for c in nb.get("cells", []) if c.get("cell_type") == "code")
            if rel.endswith("_colab.ipynb") and n_out:
                problems.append(f"{rel}: {n_out} cell output(s) in the Colab twin — re-run tools/publish_notebook.py "
                                "(was it saved from Colab with 'Save a copy in GitHub'?)")
            if BANNER.search(text):
                problems.append(f"{rel}: contains the private-run banner ('book values active')")
            if PAGE_REF.search(text):
                problems.append(f"{rel}: references rendered book pages (chapters/pages/)")
        elif p.suffix == ".html":
            text = p.read_text(encoding="utf-8", errors="replace")
            if BANNER.search(text):
                problems.append(f"{rel}: built from a private run — rebuild without the book values")
            if PAGE_REF.search(text):
                problems.append(f"{rel}: references rendered book pages (chapters/pages/)")
    return problems


def main(argv: list[str]) -> int:
    paths = argv or tracked_files()
    problems = check(paths)
    if problems:
        print("PUBLIC-REPO CHECK FAILED (CLAUDE.md rule 9):")
        for msg in problems:
            print("  -", msg)
        return 1
    print(f"public-repo check OK ({len(paths)} files)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
