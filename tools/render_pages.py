"""Render book pages as PNG images so agents can READ equations, figures and tables (the extracted text garbles maths).

Pages are addressed either by chapter-relative page number or by printed page number, and written to
chapters/pages/<id>/p<pdfpage>.png (git-ignored: page images are the book's content and are never committed).

Usage (from the repo root):
    python tools/render_pages.py ch07 --printed 283-285        # printed page numbers, as the book cites them
    python tools/render_pages.py ch07 --pdf 312                # 1-based PDF page numbers
    python tools/render_pages.py ch07 --eq 7.27                # only the page where (7.27) is defined
    python tools/render_pages.py ch07 --find "(7.27)"          # every page whose text contains the label (also citations)
    python tools/render_pages.py ch07 --printed 290 --dpi 200 --clip 0.0-0.5   # top half only, sharper

Then Read the PNG. Transcribe equations into LaTeX from the image, never from chapters/<id>.txt.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def parse_range(spec: str) -> list[int]:
    pages: list[int] = []
    for part in spec.split(","):
        if "-" in part:
            a, b = part.split("-")
            pages += list(range(int(a), int(b) + 1))
        else:
            pages.append(int(part))
    return pages


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("id", help="chapter id, e.g. ch07")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--printed", help="printed page numbers, e.g. 283-285,290")
    g.add_argument("--pdf", help="1-based PDF page numbers")
    g.add_argument("--find", help="render pages whose text contains this string, e.g. '(7.27)' (also pages citing it)")
    g.add_argument("--eq", help="render only the page(s) where equation label N.M is DEFINED (a line that is exactly '(N.M)')")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--clip", default=None, help="vertical fraction range, e.g. 0.0-0.5 for the top half")
    a = ap.parse_args()

    try:
        import fitz  # pymupdf
    except ImportError:
        sys.exit("pymupdf is required: .venv/Scripts/python.exe -m pip install pymupdf")

    book = yaml.safe_load((ROOT / "book.yaml").read_text(encoding="utf-8"))
    off = int(book["book"].get("printed_offset", 0))
    entry = next((c for c in book["chapters"] + book.get("appendices", []) if c["id"] == a.id), None)
    if entry is None:
        sys.exit(f"{a.id} not in book.yaml")
    pdf = max((p for p in ROOT.glob("*.pdf")), key=lambda p: p.stat().st_size)
    doc = fitz.open(pdf)

    if a.printed:
        pages = [p + off for p in parse_range(a.printed)]
    elif a.pdf:
        pages = parse_range(a.pdf)
    elif a.eq:
        label = f"({a.eq.strip('()')})"
        pages = [p for p in range(entry["pdf_start"], entry["pdf_end"] + 1)
                 if any(line.strip() == label for line in doc[p - 1].get_text("text").splitlines())]
        if not pages:
            print(f"no page in {a.id} defines {label} on its own line; falling back to --find")
            pages = [p for p in range(entry["pdf_start"], entry["pdf_end"] + 1) if label in doc[p - 1].get_text("text")]
    else:
        pages = [p for p in range(entry["pdf_start"], entry["pdf_end"] + 1)
                 if a.find in doc[p - 1].get_text("text")]
        if not pages:
            print(f"no page in {a.id} contains {a.find!r} (extraction may have mangled it; try a shorter string)")
            return 1

    out = ROOT / "chapters" / "pages" / a.id
    out.mkdir(parents=True, exist_ok=True)
    for p in pages:
        if not (entry["pdf_start"] <= p <= entry["pdf_end"]):
            print(f"warning: pdf page {p} is outside {a.id} ({entry['pdf_start']}-{entry['pdf_end']})")
        page = doc[p - 1]
        clip = None
        if a.clip:
            y0, y1 = (float(v) for v in a.clip.split("-"))
            r = page.rect
            clip = fitz.Rect(r.x0, r.y0 + y0 * r.height, r.x1, r.y0 + y1 * r.height)
        dest = out / f"p{p:03d}.png"
        page.get_pixmap(dpi=a.dpi, clip=clip).save(dest)
        print(f"{dest.relative_to(ROOT).as_posix()}  (printed page {p - off})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
