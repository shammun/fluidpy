"""Coverage gate for a chapter notebook: every new idea taught with code and a visual, nothing used unexplained.

Checks (exit 1 on any error):
  1. every book section from ``book.yaml`` has a notebook section (``## N.M …`` heading);
  2. every CORE row of ``analysis/chNN_curation.md`` (IDs ``C01``…) has cells tagged with ``metadata.fluidpy.core``, with
     at least one code cell and at least one visual — in an EXECUTED notebook the visual must be a real output
     (PNG image, video, plotly figure or explainer block), otherwise a figure/animation/plotly/explainer cell;
  3. every RECAP row (``R01``…) appears as a recap cell;
  4. 5–10 explainers are embedded (``book.yaml → project.min/max_explainers_per_chapter``);
  5. the prerequisite ledger in ``analysis/chNN_design.md`` (``## Part E`` table: ``| Concept | First used in | Explained by |``)
     — every row whose "Explained by" says ``primer`` has a primer cell naming that concept (warning if not found);
     rows that say ``C07``/``R02`` must point at an existing CORE/RECAP id (error);
  6. every DERIVATION row (``D01``… in the curation's derivations table) is derived step by step (``nb.derivation``)
     — ★★★ rows also need a sympy check cell (no error when executed), and every explainer named in the row's last
     column lists the id in its ``viz:derivations`` meta;
  7. no code cell is left uncommented: code cells with ≥ 4 non-blank lines need a comment on at least half of them
     (warning — the lesson-reviewer judges quality; this catches the obvious case);
  8. equations are shown, not just cited: a markdown cell that names a book equation number — ``(3.5)``, ``Eq. 3.5`` —
     must also contain the equation in LaTeX (``$…$``/``$$…$$``) (error). The reviewer checks it is the RIGHT equation;
  9. no cell contains a control character (\\f, \\b, \\v, \\a) — the sign of a LaTeX backslash eaten by a non-raw string (error).

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

    # 6 derivations: every D row is written out step by step inside its CORE block; ★★★ ones carry a sympy check;
    #   an explainer named in the row lists the id in its viz:derivations meta (and so has the Derivation tab)
    der_cells: dict[str, set[str]] = {}
    for c in nb.cells:
        did = c.metadata.get("fluidpy", {}).get("derivation")
        if did:
            der_cells.setdefault(did, set()).update(c.metadata.get("tags", []))
            if c.cell_type == "code" and has_outputs and any(o.get("output_type") == "error" for o in c.get("outputs", [])):
                errors.append(f"DERIVATION {did}: the sympy check cell raised an error")
    for did, r in planned.items():
        if r["tier"] != "DERIVATION":
            continue
        tags = der_cells.get(did)
        if not tags or "derivation" not in tags:
            errors.append(f"DERIVATION {did} ({r['item']}) is not derived step by step in the notebook")
            continue
        if r.get("hard") and "derivation-check" not in tags:
            errors.append(f"DERIVATION {did} ({r['item']}) is ★★★ but has no sympy check cell")
        for slug in r.get("explainers", []):
            html = ROOT / "viz" / chapter / f"{slug}.html"
            if not html.exists():
                warns.append(f"DERIVATION {did}: explainer {slug} not built yet (its Derivation tab is not checked)")
                continue
            m = re.search(r'<meta\s+name="viz:derivations"\s+content="([^"]*)"', html.read_text(encoding="utf-8"))
            if not m or did not in re.split(r"[\s,;]+", m.group(1)):
                errors.append(f"DERIVATION {did}: explainer {slug} does not list it in viz:derivations (no Derivation tab for it)")

    # 7 comments
    for i, c in enumerate(nb.cells):
        if c.cell_type != "code" or "setup" in c.metadata.get("tags", []):
            continue
        lines = [l for l in c.source.splitlines() if l.strip()]
        if len(lines) >= 4 and sum("#" in l for l in lines) < len(lines) / 2:
            warns.append(f"cell {i}: {len(lines)} lines but only {sum('#' in l for l in lines)} commented")

    # 8 equations shown, not just cited: a markdown cell naming a book equation number must also show maths
    for i, c in enumerate(nb.cells):
        if c.cell_type != "markdown":
            continue
        nums = sorted(set(EQ_REF.findall(c.source)))
        if nums and "$" not in c.source and "\\begin{" not in c.source:
            errors.append(f"cell {i}: cites Eq. {', '.join(nums)} by number only — show the equation next to its number")

    # 9 no control characters: a non-raw builder string turns \frac, \beta, \vec, \nabla… into \f, \b, \v, \n…
    for i, c in enumerate(nb.cells):
        bad = sorted({repr(ch) for ch in c.source if ch in "\x07\x08\x0b\x0c"})
        if bad:
            errors.append(f"cell {i}: control character(s) {', '.join(bad)} — a LaTeX backslash was eaten; use a raw string r\"…\"")
    return errors, warns


# "(3.5)", "Eq. 3.5", "Eqs. (3.5)", "(3.11, 3.12)" — a book equation number (chapter digits or an appendix letter)
EQ_REF = re.compile(r"(?:\bEqs?\.\s*\(?|\()((?:\d+|[A-D])\.\d+)(?=[\s,)–-])")


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
