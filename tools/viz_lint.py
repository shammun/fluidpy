"""Static lint for interactive explainers — fast checks that run before the (slower) browser audit in tools/shot.py.

Errors (fail):
  * missing or placeholder ``<meta name="viz:*">`` tags (chapter, slug, order, title, summary, concept, sections,
    equations, fluidpy); ``viz:chapter``/``viz:slug`` must match the file's folder/name;
  * missing VIZ_BASE_CSS / VIZ_LIB_JS marker blocks, or no ``Viz.app(`` call;
  * external resources other than KaTeX (which viz_lib loads itself): ``<script src=``, ``<link href=``, ``@import``,
    ``fetch(``, ``XMLHttpRequest``, ``<img src="http``;
  * no ``tour:``, ``equations:`` or ``selftest:`` section; no ``ref:`` (book equation number) in the equations;
    no ``py:`` parity row in selftest (explainers in a chapter folder must mirror a fluidpy function);
  * chapter CSS that re-introduces scrolling (``overflow: auto|scroll``, ``overflow-y``) or tiny text
    (``font-size`` below 12px) or fixed pixel widths above 360px;
  * file larger than 600 kB;
  * more than ``max_explainers_per_chapter`` explainers in the folder.
Warnings: ``console.log`` left in; very long walkthrough texts (> 90 words in one step, heuristic).

Usage::  .venv/Scripts/python.exe tools/viz_lint.py viz/ch07/*.html      or   --chapter ch07   or   --all
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REQUIRED_META = ("chapter", "slug", "order", "title", "summary", "concept", "sections", "equations", "fluidpy")
MAX_BYTES = 600_000


def _chapter_part(html: str) -> tuple[str, str]:
    """(chapter CSS, chapter JS) — i.e. everything that is not the inlined library."""
    css = re.sub(r"/\* VIZ_BASE_CSS:BEGIN.*?/\* VIZ_BASE_CSS:END \*/", "", html, flags=re.S)
    js = re.sub(r"/\* VIZ_LIB_JS:BEGIN \*/.*?/\* VIZ_LIB_JS:END \*/", "", css, flags=re.S)
    styles = "\n".join(re.findall(r"<style[^>]*>(.*?)</style>", js, flags=re.S))
    scripts = "\n".join(re.findall(r"<script[^>]*>(.*?)</script>", js, flags=re.S))
    return styles, scripts


def lint_text(path: Path, html: str) -> list[str]:
    msgs: list[str] = []
    meta = dict(re.findall(r'<meta\s+name="viz:([\w-]+)"\s+content="([^"]*)"', html))
    for k in REQUIRED_META:
        if not meta.get(k) or meta[k].startswith("__"):
            msgs.append(f"meta viz:{k} missing or placeholder")
    in_chapter = re.fullmatch(r"ch\d\d", path.parent.name) is not None
    if in_chapter:
        if meta.get("chapter") and meta["chapter"] != path.parent.name:
            msgs.append(f"viz:chapter={meta['chapter']!r} but the file is in {path.parent.name}/")
        if meta.get("slug") and meta["slug"] != path.stem:
            msgs.append(f"viz:slug={meta['slug']!r} but the file is {path.name}")
    if "/* VIZ_BASE_CSS:BEGIN" not in html or "/* VIZ_LIB_JS:BEGIN */" not in html:
        msgs.append("missing VIZ marker blocks (start from templates/viz_template.html via tools/new_viz.py)")
    css, js = _chapter_part(html)
    if "Viz.app(" not in js:
        msgs.append("no Viz.app({...}) call")
    rest = re.sub(r"/\* VIZ_LIB_JS:BEGIN \*/.*?/\* VIZ_LIB_JS:END \*/", "", html, flags=re.S)
    for pat, what in [(r"<script[^>]+src=", "<script src=…> (inline everything; KaTeX is loaded by viz_lib)"),
                      (r"<link[^>]+href=", "<link href=…>"), (r"@import", "@import"), (r"\bfetch\(", "fetch("),
                      (r"XMLHttpRequest", "XMLHttpRequest"), (r"<img[^>]+src=\"https?:", "remote <img>")]:
        if re.search(pat, rest):
            msgs.append(f"external resource: {what}")
    for section in ("tour:", "equations:", "selftest:"):
        if section not in js:
            msgs.append(f"no {section[:-1]} section")
    if "equations:" in js and "ref:" not in js:
        msgs.append("equations carry no ref: 'Eq. (N.M)' (cite the book's equation numbers)")
    if in_chapter and "selftest:" in js and not re.search(r"\bpy\s*:", js):
        msgs.append("selftest has no py: parity row against a fluidpy function")
    if re.search(r"overflow(-[xy])?\s*:\s*(auto|scroll)", css):
        msgs.append("chapter CSS sets overflow auto/scroll (the no-scroll rule: use a tab or the pager)")
    for m in re.finditer(r"font-size\s*:\s*([\d.]+)px", css):
        if float(m.group(1)) < 12:
            msgs.append(f"chapter CSS font-size {m.group(1)}px < 12px")
    for m in re.finditer(r"(?<![-\w])(min-)?width\s*:\s*(\d+)px", css):
        if int(m.group(2)) > 360:
            msgs.append(f"chapter CSS fixed width {m.group(2)}px > 360px (breaks phones)")
    size = len(html.encode("utf-8"))
    if size > MAX_BYTES:
        msgs.append(f"file is {size / 1e3:.0f} kB (> {MAX_BYTES / 1e3:.0f} kB)")
    if "console.log(" in js:
        msgs.append("WARN console.log left in the chapter script")
    for m in re.finditer(r"text\s*:\s*'((?:[^'\\]|\\.)*)'", js):
        if len(m.group(1).split()) > 90:
            msgs.append(f"WARN a walkthrough/explore text has {len(m.group(1).split())} words (keep steps short; phones)")
    return msgs


def lint_files(files: list[Path]) -> dict[Path, list[str]]:
    out = {f: lint_text(f, f.read_text(encoding="utf-8")) for f in files}
    try:
        import yaml
        limit = int(yaml.safe_load((ROOT / "book.yaml").read_text(encoding="utf-8"))["project"]["max_explainers_per_chapter"])
    except Exception:  # noqa: BLE001
        limit = 5
    by_folder: dict[Path, int] = {}
    for f in files:
        if re.fullmatch(r"ch\d\d", f.parent.name):
            by_folder[f.parent] = len([p for p in f.parent.glob("*.html") if p.name != "index.html"])
    for folder, n in by_folder.items():
        if n > limit:
            first = next(f for f in files if f.parent == folder)
            out[first].append(f"{folder.name} has {n} explainers (limit {limit})")
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("paths", nargs="*")
    ap.add_argument("--chapter")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args(argv)
    files = [Path(p) if Path(p).is_absolute() else ROOT / p for p in a.paths]
    if a.chapter:
        files += sorted(p for p in (ROOT / "viz" / a.chapter).glob("*.html") if p.name != "index.html")
    if a.all:
        files += sorted(p for p in (ROOT / "viz").glob("ch??/*.html") if p.name != "index.html")
    if not files:
        ap.error("give files, --chapter or --all")
    bad = 0
    for f, msgs in lint_files(files).items():
        errors = [m for m in msgs if not m.startswith("WARN")]
        rel = f.relative_to(ROOT).as_posix() if f.is_relative_to(ROOT) else str(f)
        print(("FAIL " if errors else "ok   ") + rel)
        for m in msgs:
            print("   -", m)
        bad += bool(errors)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
