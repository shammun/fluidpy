"""Inline the shared explainer library into explainer HTML files, so every explainer stays ONE self-contained file.

Each explainer (``viz/chNN/<slug>.html``, ``templates/viz_*.html``) carries two marker blocks::

    <style id="viz-base">/* VIZ_BASE_CSS:BEGIN */ … /* VIZ_BASE_CSS:END */</style>
    <script id="viz-lib">/* VIZ_LIB_JS:BEGIN */ … /* VIZ_LIB_JS:END */</script>

This tool replaces whatever sits between the markers with the current body of ``assets/viz_base.css`` and
``assets/viz_lib.js`` (the part between the same markers in those files). Edit the assets, never the inlined copies.

Why inline instead of ``<script src="../../assets/viz_lib.js">``: the explainer is shown through ``<iframe srcdoc>`` in
Jupyter/VS Code/Colab, where relative URLs do not resolve; a self-contained file works in every host, offline too.

Usage (repo root)::

    .venv/Scripts/python.exe tools/viz_inline.py viz/ch07/dispersion_relation.html
    .venv/Scripts/python.exe tools/viz_inline.py --all          # every viz/**/*.html + templates/viz_*.html
    .venv/Scripts/python.exe tools/viz_inline.py --all --check  # exit 1 if any file is out of date (used by viz_lint)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "assets"

BLOCKS = {
    "css": ("/* VIZ_BASE_CSS:BEGIN", "/* VIZ_BASE_CSS:END */", ASSETS / "viz_base.css"),
    "js": ("/* VIZ_LIB_JS:BEGIN */", "/* VIZ_LIB_JS:END */", ASSETS / "viz_lib.js"),
}


def asset_body(kind: str) -> str:
    """The text between the markers in the asset file (markers excluded, surrounding newlines trimmed)."""
    begin, end, path = BLOCKS[kind]
    text = path.read_text(encoding="utf-8")
    i = text.find(begin)
    j = text.find(end)
    if i < 0 or j < 0:
        raise SystemExit(f"{path}: markers {begin!r} / {end!r} not found")
    i = text.find("*/", i) + 2  # skip the rest of the BEGIN comment line
    return text[i:j].strip("\n")


def inline_text(html: str, bodies: dict[str, str]) -> tuple[str, list[str]]:
    """Return (new_html, missing_marker_kinds)."""
    missing = []
    for kind, (begin, end, _path) in BLOCKS.items():
        start = html.find(begin)
        stop = html.find(end, start + 1) if start >= 0 else -1
        if start < 0 or stop < 0:
            missing.append(kind)
            continue
        if kind == "css":
            head = "/* VIZ_BASE_CSS:BEGIN (inlined from assets/viz_base.css by tools/viz_inline.py) */"
        else:
            head = "/* VIZ_LIB_JS:BEGIN */"
        html = html[:start] + head + "\n" + bodies[kind] + "\n" + html[stop:]
    return html, missing


def targets(paths: list[str], all_: bool) -> list[Path]:
    out: list[Path] = []
    if all_:
        out += sorted(p for p in (ROOT / "viz").glob("**/*.html") if p.name != "index.html")
        out += sorted((ROOT / "templates").glob("viz_*.html"))
    for p in paths:
        q = Path(p)
        out.append(q if q.is_absolute() else ROOT / q)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("paths", nargs="*")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--check", action="store_true", help="do not write; exit 1 if a file is out of date")
    a = ap.parse_args(argv)
    files = targets(a.paths, a.all)
    if not files:
        ap.error("give explainer paths or --all")
    bodies = {k: asset_body(k) for k in BLOCKS}
    stale, bad = [], []
    for f in files:
        html = f.read_text(encoding="utf-8")
        new, missing = inline_text(html, bodies)
        rel = f.relative_to(ROOT).as_posix() if f.is_relative_to(ROOT) else str(f)
        if missing:
            bad.append(f"{rel}: missing marker block(s) {missing}")
            continue
        if new != html:
            stale.append(rel)
            if not a.check:
                f.write_text(new, encoding="utf-8")
    for b in bad:
        print("ERROR", b)
    if a.check:
        for s in stale:
            print("STALE", s, "(run tools/viz_inline.py)")
        return 1 if (stale or bad) else 0
    print(f"inlined viz_base.css + viz_lib.js into {len(stale)} file(s); {len(files) - len(stale) - len(bad)} already current")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
