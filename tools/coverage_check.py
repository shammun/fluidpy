"""Coverage gate for a chapter notebook: every new idea taught with code and a visual, nothing used unexplained.

Checks (exit 1 on any error):
  1. every book section from ``book.yaml`` has a notebook section (``## N.M …`` heading);
  2. every CORE row of ``analysis/chNN_curation.md`` (IDs ``C01``…) has cells tagged with ``metadata.fluidpy.core``, with
     at least one code cell and at least one visual — in an EXECUTED notebook the visual must be a real output
     (PNG image, video, plotly figure or explainer block), otherwise a figure/animation/plotly/explainer cell;
  3. every RECAP row (``R01``…) appears as a recap cell;
  4. 4–5 explainers are embedded (``book.yaml → project.min/max_explainers_per_chapter``);
  5. the prerequisite ledger in ``analysis/chNN_design.md`` (``## Part E`` table: ``| Concept | First used in | Explained by |``)
     — every row whose "Explained by" says ``primer`` has a primer cell naming that concept (warning if not found);
     rows that say ``C07``/``R02`` must point at an existing CORE/RECAP id (error);
  6. no code cell is left uncommented: code cells with ≥ 4 non-blank lines need a comment on at least half of them
     (warning — the lesson-reviewer judges quality; this catches the obvious case).

Usage::

    .venv/Scripts/python.exe tools/coverage_check.py ch07                 # source notebook (tags only)
    .venv/Scripts/python.exe tools/coverage_check.py ch07 --executed      # executes first (fluidpy-venv), checks outputs
    .venv/Scripts/python.exe tools/coverage_check.py ch07 --nb outputs/ch07/executed.ipynb
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))


def _visual_output(cell) -> bool:
    for o in cell.get("outputs", []):
        data = o.get("data", {})
        if "image/png" in data or "image/svg+xml" in data or "application/vnd.plotly.v1+json" in data:
            return True
        html = data.get("text/html", "")
        html = "".join(html) if isinstance(html, list) else html
        if "<video" in html or "fluidpy-viz:BEGIN" in html or "animation" in html and "<img" in html or "plotly" in html.lower():
            return True
    return False


def check(chapter: str, nb_path: Path | None = None, executed: bool = False) -> tuple[list[str], list[str]]:
    import nbformat
    import yaml

    import nbkit

    book = yaml.safe_load((ROOT / "book.yaml").read_text(encoding="utf-8"))
    row = next(c for c in book["chapters"] if c["id"] == chapter)
    errors: list[str] = []
    warns: list[str] = []
    path = nb_path or (ROOT / "notebooks" / f"{chapter}_{row['slug']}.ipynb")
    if executed and nb_path is None:
        import publish_notebook as pub
        nb, _ = pub.execute(path)
    else:
        nb = nbformat.read(path, as_version=4)
    has_outputs = any(c.get("outputs") for c in nb.cells if c.cell_type == "code")

    # 1 sections
    headings = " ".join(c.source for c in nb.cells if c.cell_type == "markdown")
    for s in row.get("sections", []):
        m = re.match(r"\s*(\d+\.\d+)", s)
        if m and not re.search(rf"^##\s+[^\n]*\b{re.escape(m.group(1))}\b", headings, flags=re.M):
            errors.append(f"book section {m.group(1)} has no notebook section")

    # 2/3 CORE and RECAP
    planned = nbkit.curation_items(chapter)
    if not planned:
        warns.append(f"analysis/{chapter}_curation.md has no tiers table with IDs (C01…) — CORE coverage not checked")
    blocks: dict[str, dict] = {}
    for c in nb.cells:
        cid = c.metadata.get("fluidpy", {}).get("core")
        tags = set(c.metadata.get("tags", []))
        if not cid:
            continue
        b = blocks.setdefault(cid, {"code": 0, "visual": 0})
        if c.cell_type == "code":
            b["code"] += 1
            if has_outputs:
                b["visual"] += _visual_output(c)
            elif tags & nbkit.VISUAL_TAGS or nbkit.VISUAL_CODE.search(c.source):
                b["visual"] += 1
    meta = nb.metadata.get("fluidpy", {})
    for cid, r in planned.items():
        if r["tier"] == "CORE":
            b = blocks.get(cid)
            if not b:
                errors.append(f"CORE {cid} ({r['item']}) is not taught in the notebook")
            else:
                if not b["code"]:
                    errors.append(f"CORE {cid} ({r['item']}) has no code")
                if not b["visual"]:
                    errors.append(f"CORE {cid} ({r['item']}) has no visual{' output' if has_outputs else ''}")
        if r["tier"] == "RECAP" and cid not in meta.get("recaps", []):
            errors.append(f"RECAP {cid} ({r['item']}) is not recapped in the notebook")

    # 4 explainers
    lo, hi = nbkit.explainer_limits(book)
    n_viz = len(set(nbkit.explainer_calls(path))) if nb_path is None else sum(
        1 for c in nb.cells if c.cell_type == "code" and "show_viz(" in c.source)
    if not (lo <= n_viz <= hi):
        errors.append(f"{n_viz} explainers embedded (need {lo}–{hi})")

    # 5 ledger
    design = ROOT / "analysis" / f"{chapter}_design.md"
    primers = [p.lower() for p in meta.get("primers", [])]
    if design.exists():
        text = design.read_text(encoding="utf-8")
        part = text.split("## Part E", 1)[1] if "## Part E" in text else ""
        rows = [[x.strip() for x in l.strip().strip("|").split("|")] for l in part.splitlines() if l.strip().startswith("|")]
        rows = [r for r in rows if len(r) >= 3 and not set("".join(r)) <= set("-: ") and r[0].lower() != "concept"]
        if not rows:
            warns.append("design Part E (prerequisite ledger) is missing or empty")
        for concept, _first, by, *rest in rows:
            ids = re.findall(r"\b([CR]\d{2,3})\b", by)
            for i in ids:
                if i not in planned:
                    errors.append(f"ledger: '{concept}' is explained by {i}, which is not in the curation")
            if "primer" in by.lower():
                key = re.sub(r"[`*$\\]", "", concept).lower().split("(")[0].strip()
                if key and not any(key[:18] in p or p[:18] in key for p in primers):
                    warns.append(f"ledger: primer for '{concept}' not found among the notebook's primers")
    else:
        warns.append(f"analysis/{chapter}_design.md not found — ledger not checked")

    # 6 comments
    for i, c in enumerate(nb.cells):
        if c.cell_type != "code" or "setup" in c.metadata.get("tags", []):
            continue
        lines = [l for l in c.source.splitlines() if l.strip()]
        if len(lines) >= 4 and sum("#" in l for l in lines) < len(lines) / 2:
            warns.append(f"cell {i}: {len(lines)} lines but only {sum('#' in l for l in lines)} commented")
    return errors, warns


def main(argv: list[str] | None = None) -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("chapter")
    ap.add_argument("--executed", action="store_true", help="execute the notebook first and check real outputs")
    ap.add_argument("--nb", help="check this notebook file instead (e.g. an executed copy)")
    a = ap.parse_args(argv)
    errors, warns = check(a.chapter, Path(a.nb) if a.nb else None, a.executed)
    for w in warns:
        print("WARN ", w)
    for e in errors:
        print("ERROR", e)
    print(("COVERAGE OK" if not errors else "COVERAGE FAILED") + f" for {a.chapter} ({len(errors)} errors, {len(warns)} warnings)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
