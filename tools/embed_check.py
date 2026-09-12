"""Integration gate: are a chapter's explainers built, current, lint-clean and embedded exactly once?

Checks for chapter ``chNN``:
  1. ``viz/chNN/*.html`` holds 1–5 explainers (``book.yaml → project.max_explainers_per_chapter``);
  2. each is inlined with the current ``assets/viz_lib.js`` + ``viz_base.css`` (``tools/viz_inline.py --check``) and
     passes ``tools/viz_lint.py``;
  3. the notebook ``notebooks/chNN_<slug>.ipynb`` calls ``show_viz("chNN", "<slug>")`` exactly once per explainer,
     and never for a slug that does not exist;
  4. the design document's explainer shortlist (``analysis/chNN_design.md``, lines like ``### E1 · <slug>``) matches
     the files on disk (warning only — the design may legitimately drop one).

Usage::  .venv/Scripts/python.exe tools/embed_check.py ch07        (exit 0 = OK)
"""
from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))


def check(chapter: str) -> list[str]:
    import yaml

    import viz_inline
    import viz_lint
    from nbkit import explainer_calls

    book = yaml.safe_load((ROOT / "book.yaml").read_text(encoding="utf-8"))
    limit = int(book.get("project", {}).get("max_explainers_per_chapter", 5))
    row = next((c for c in book["chapters"] if c["id"] == chapter), None)
    if row is None:
        return [f"{chapter} not in book.yaml"]
    problems: list[str] = []
    files = sorted(p for p in (ROOT / "viz" / chapter).glob("*.html") if p.name != "index.html")
    slugs = [p.stem for p in files]
    if not files:
        problems.append(f"no explainers in viz/{chapter}/")
    if len(files) > limit:
        problems.append(f"{len(files)} explainers in viz/{chapter}/ (limit {limit})")

    if files:
        bodies = {k: viz_inline.asset_body(k) for k in viz_inline.BLOCKS}
        for f in files:
            new, missing = viz_inline.inline_text(f.read_text(encoding="utf-8"), bodies)
            if missing:
                problems.append(f"{f.name}: missing VIZ marker blocks {missing}")
            elif new != f.read_text(encoding="utf-8"):
                problems.append(f"{f.name}: library copy is stale — run tools/viz_inline.py --all")
        for f, msgs in viz_lint.lint_files(files).items():
            problems += [f"{f.name}: {m}" for m in msgs if not m.startswith("WARN")]

    nb = ROOT / "notebooks" / f"{chapter}_{row['slug']}.ipynb"
    if not nb.exists():
        problems.append(f"{nb.relative_to(ROOT).as_posix()} does not exist")
    else:
        calls = explainer_calls(nb)
        counts = Counter(calls)
        for s in slugs:
            if counts[s] == 0:
                problems.append(f"explainer {s} is not embedded in {nb.name}")
            elif counts[s] > 1:
                problems.append(f"explainer {s} is embedded {counts[s]} times in {nb.name}")
        for s in counts:
            if s not in slugs:
                problems.append(f"{nb.name} calls show_viz for missing explainer '{s}'")

    design = ROOT / "analysis" / f"{chapter}_design.md"
    if design.exists():
        planned = re.findall(r"^###\s+E\d+\s*[·\-:]\s*`?([\w\-]+)`?", design.read_text(encoding="utf-8"), flags=re.M)
        extra = sorted(set(slugs) - set(planned))
        dropped = sorted(set(planned) - set(slugs))
        if extra:
            print(f"WARN explainers not in the design shortlist: {extra}")
        if dropped:
            print(f"WARN design shortlist items with no file: {dropped}")
    return problems


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    problems = check(argv[0])
    if problems:
        print(f"EMBED CHECK FAILED for {argv[0]}:")
        for p in problems:
            print("  -", p)
        return 1
    print(f"embed check OK for {argv[0]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
