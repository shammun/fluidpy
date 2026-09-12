"""Read the PDF's bookmark outline and compare it with the chapter map in book.yaml (read-only).

book.yaml was generated from this outline once and then annotated (slugs, blurbs, viz_seeds), so this tool never
rewrites it. It prints the outline's chapters and sections with exact PDF pages and flags every difference from
book.yaml — use it if you replace the PDF (another printing) or suspect a page range.

Usage::  .venv/Scripts/python.exe tools/detect_chapters.py            (exit 0 = book.yaml matches the outline)
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    import fitz

    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    pdf = max(ROOT.glob("*.pdf"), key=lambda p: p.stat().st_size)
    doc = fitz.open(pdf)
    toc = doc.get_toc()  # [level, title, page]
    book = yaml.safe_load((ROOT / "book.yaml").read_text(encoding="utf-8"))
    starts = []
    for lvl, title, page in toc:
        m = re.match(r"^\s*(?:Chapter\s*)?(\d{1,2})\s*[.:\-]?\s+(.*)$", title)
        if m and lvl <= 2 and not re.match(r"^\d+\.\d+", title.strip()):
            starts.append((int(m.group(1)), m.group(2).strip(), page))
    print(f"{pdf.name}: {doc.page_count} pages, {len(toc)} outline entries, {len(starts)} chapter-like entries")
    problems = 0
    by_num = {n: (t, p) for n, t, p in starts}
    for c in book["chapters"]:
        n = c["number"]
        if n not in by_num:
            print(f"  ? ch{n:02d} not found in the outline (titles may be formatted differently) — check by eye")
            continue
        title, page = by_num[n]
        ok = page == c["pdf_start"]
        problems += not ok
        print(f"  {'ok ' if ok else 'DIFF'} ch{n:02d} outline p{page} vs book.yaml p{c['pdf_start']} | {title[:60]}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
