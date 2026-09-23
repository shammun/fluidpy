"""List explainer text that cites a book equation by number without showing the equation.

House rule (CLAUDE.md rule 3): wherever a book equation is referred to, the equation itself is written (TeX) next to its
number. This scans the JS string literals of an explainer (after the inlined library) and prints every string that
contains an equation number such as ``(3.5)`` / ``Eq. 3.5`` but no TeX (no ``$`` and no backslash command). It is a
heuristic — the viz-reviewer judges each hit: metadata lists (``<meta name="viz:equations">``), a bare ``ref:`` label
on an equation card that already shows its TeX, and code strings are fine; tour, Explain, Derivation, quiz, notes and
status text are not.

Usage::

    .venv/Scripts/python.exe tools/eq_refs.py viz/ch02/rotation_of_axes.html [...]
    .venv/Scripts/python.exe tools/eq_refs.py --chapter ch02
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EQ = re.compile(r"(?:\bEqs?\.\s*\(?|(?<![\w.])\()((?:\d+|[A-D])\.\d+)(?=[\s,)–-])")
STR = re.compile(r"'((?:[^'\\\n]|\\.)*)'|\"((?:[^\"\\\n]|\\.)*)\"|`((?:[^`\\]|\\.)*)`")
ONLY_REF = re.compile(r"\s*(?:Eqs?\.\s*)?\(?(?:\d+|[A-D])\.\d+\)?(?:\s*[,;]\s*\(?(?:\d+|[A-D])\.\d+\)?)*\s*")
CODE_LIKE = re.compile(r"\b(?:np|sp|ch\d\d|fluidpy|Math)\.|\w\(\d")


def bare_refs(path: Path) -> list[tuple[int, str]]:
    """Return (line number, string) for every JS string that cites an equation number but shows no TeX."""
    s = path.read_text(encoding="utf-8")
    start = s.find("fluidpy-viz-lib:END")
    offset = s[:start].count("\n") if start >= 0 else 0
    body = s[start:] if start >= 0 else s
    hits = []
    for m in STR.finditer(body):
        x = next(g for g in m.groups() if g is not None)
        if not EQ.search(x) or "$" in x or "\\" in x or ONLY_REF.fullmatch(x) or CODE_LIKE.search(x):
            continue
        hits.append((offset + body[: m.start()].count("\n") + 1, x))
    return hits


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("files", nargs="*")
    ap.add_argument("--chapter")
    a = ap.parse_args(argv)
    files = [Path(f) for f in a.files] + (sorted((ROOT / "viz" / a.chapter).glob("*.html")) if a.chapter else [])
    total = 0
    for f in files:
        hits = bare_refs(f)
        total += len(hits)
        print(f"{len(hits):4d}  {f}")
        for line, x in hits:
            print(f"      L{line}: {x[:160]}")
    print(f"{total} string(s) cite an equation number without showing the equation")
    return 1 if total else 0


if __name__ == "__main__":
    raise SystemExit(main())
