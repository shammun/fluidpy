"""Fail if anything book-derived is about to be (or has been) committed to the PUBLIC repository (CLAUDE.md rule 9).

Checks every tracked file (``git ls-files``; after ``git add`` that includes staged files), or the paths given:

* forbidden paths: ``*.pdf``, ``chapters/`` (except ``chapters/_page_map.json``: page numbers only), ``data/book/``,
  ``data/manual/``, ``outputs/``, ``reports/viz/``, ``reports/**/figures/``, ``notebooks/executed_*``,
  ``tests/book_values_*``;
* ``*_colab.ipynb`` must have **no outputs** (Colab's "Save a copy in GitHub" writes them);
* no notebook or HTML page may carry the private-run banner ``book values active``;
* explainers and pages must not reference rendered book pages (``chapters/pages``);
* any tracked file above 50 MB (GitHub refuses 100 MB; pages that large are also unusable on phones);
* no tracked text file may contain a string listed under ``"_forbidden_public"`` in a local (git-ignored)
  ``tests/book_values_*.json`` — book-quoted run parameters and printed values that must stay private;
* no tracked text file may match a regular expression listed under ``"_forbidden_public_regex"`` in one of the local
  (git-ignored) private lists, ``tests/book_values_*.json`` or ``data/book/forbidden_public*.json`` — for book-quoted
  *constants*, which are too short to be literal strings (a bare one- or two-digit number) and are therefore keyed by
  the name they are assigned to (``<name> = <value>``, ``"<name>": <value>``, a pair of values on one line).  An entry
  is either a pattern string or
  ``{"regex": "...", "why": "...", "paths": "<regex on the repo-relative path, optional>"}``.  The patterns themselves
  are book-derived and therefore never written in this (public) file.

Usage::

    .venv/Scripts/python.exe tools/check_public.py            # all tracked files
    .venv/Scripts/python.exe tools/check_public.py a.ipynb b.html
    git config core.hooksPath .githooks                       # installs the pre-push hook that runs this

Exit code 0 = clean, 1 = violations (listed), 2 = a local private list could not be parsed (nothing was checked).
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
TEXT_SUFFIXES = {".md", ".py", ".ipynb", ".html", ".csv", ".json", ".yaml", ".yml", ".txt", ".js", ".css"}


PRIVATE_LISTS = ("tests/book_values_*.json", "data/book/forbidden_public*.json")  # git-ignored, local only
MAX_HITS = 3  # matches reported per (file, pattern)


def _private_lists() -> list[dict]:
    """The parsed local private lists (absent on a fresh clone → no content check).

    A private list that exists but cannot be read or parsed stops the check (exit code 2, file and error printed):
    skipping it would drop every forbidden string and pattern it holds without a word.
    """
    out: list[dict] = []
    for pattern in PRIVATE_LISTS:
        for f in sorted(ROOT.glob(pattern)):
            try:
                data = json.loads(f.read_text(encoding="utf-8"))
            except Exception as exc:  # noqa: BLE001
                rel = f.relative_to(ROOT).as_posix()
                print(f"PUBLIC-REPO CHECK ABORTED: private list {rel} cannot be parsed "
                      f"({type(exc).__name__}: {exc}).\n  Its forbidden strings and patterns would be skipped — "
                      "fix the JSON and run the check again.", file=sys.stderr)
                raise SystemExit(2) from exc
            if isinstance(data, dict):
                out.append(data)
    return out


def forbidden_strings() -> list[str]:
    """Private literal strings from the local book-value files (absent on a fresh clone → no content check)."""
    out: list[str] = []
    for data in _private_lists():
        out += [s for s in data.get("_forbidden_public", []) if s and isinstance(s, str)]
    return out


def forbidden_patterns() -> list[tuple[re.Pattern, re.Pattern | None, str]]:
    """Private regular expressions for book-quoted constants: (pattern, path filter or None, why).

    Read from ``"_forbidden_public_regex"`` of the local private lists; a malformed entry raises (a checker that
    silently drops a rule is worse than one that stops).
    """
    out: list[tuple[re.Pattern, re.Pattern | None, str]] = []
    for data in _private_lists():
        for entry in data.get("_forbidden_public_regex", []):
            if isinstance(entry, str):
                entry = {"regex": entry}
            paths = entry.get("paths")
            out.append((re.compile(entry["regex"]), re.compile(paths) if paths else None, entry.get("why", "")))
    return out


def pattern_hits(rel: str, body: str, patterns: list[tuple[re.Pattern, re.Pattern | None, str]]) -> list[str]:
    """Problem lines for every private pattern that matches ``body`` (at most MAX_HITS per pattern, with line numbers)."""
    problems: list[str] = []
    for pat, paths, why in patterns:
        if paths is not None and not paths.search(rel):
            continue
        for n, m in enumerate(pat.finditer(body)):
            if n >= MAX_HITS:
                problems.append(f"{rel}: … further matches of the same private pattern not listed")
                break
            line = body.count("\n", 0, m.start()) + 1
            problems.append(f"{rel}:{line}: matches a private book-constant pattern ({why or 'see the private list'}): "
                            f"{m.group(0)!r}")
    return problems


def tracked_files() -> list[str]:
    out = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True).stdout
    return [line for line in out.splitlines() if line]


def check(paths: list[str]) -> list[str]:
    problems: list[str] = []
    private = forbidden_strings()
    patterns = forbidden_patterns()
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
        if (private or patterns) and p.suffix.lower() in TEXT_SUFFIXES:
            body = p.read_text(encoding="utf-8", errors="replace")
            for s in private:
                if s in body:
                    problems.append(f"{rel}: contains a private book value listed in tests/book_values_*.json: {s!r}")
            problems += pattern_hits(rel, body, patterns)
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
