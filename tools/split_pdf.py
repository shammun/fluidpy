"""Split the book PDF into one PDF + one text file per chapter, using the exact PDF pages recorded in book.yaml.

book.yaml was generated from the PDF's own bookmark outline, so every chapter carries ``pdf_start``/``pdf_end``
(1-based, inclusive). There is no offset guessing (the failure that cost the sea-ice project a day).

Writes (all git-ignored, the book never leaves this machine):
    chapters/chNN.pdf            the chapter's pages
    chapters/chNN.txt            extracted text, with a header per page:  ===== PDF page 123 (printed page 96) =====
    chapters/_page_map.json      id -> pdf/printed ranges (safe to commit: numbers only)

Usage (from the repo root):
    python tools/split_pdf.py                 # all chapters (+ appendices with --appendices)
    python tools/split_pdf.py ch04 ch07       # only these
    python tools/split_pdf.py --dry-run       # print the mapping and each chapter's first text lines
    python tools/split_pdf.py --check         # verify every chapter's first page really is its opener

The text extraction garbles mathematics (minus signs become "/C0", "=" becomes "¼", Greek letters shift). Agents must
render equation pages as images with tools/render_pages.py and read the image, never trust the .txt for an equation.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def load_book() -> dict:
    return yaml.safe_load((ROOT / "book.yaml").read_text(encoding="utf-8"))


def locate_book_pdf(book: dict) -> Path:
    spec = str(book.get("book", {}).get("pdf", "auto") or "auto")
    if spec != "auto":
        p = ROOT / spec
        if not p.exists():
            sys.exit(f"ERROR: book.yaml -> book.pdf = {spec!r} does not exist in {ROOT}")
        return p
    pdfs = [p for p in ROOT.glob("*.pdf") if p.is_file()]
    if not pdfs:
        sys.exit(f"ERROR: no *.pdf in the repo root {ROOT}. Put the book PDF there or set book.pdf in book.yaml.")
    return max(pdfs, key=lambda p: p.stat().st_size)


LIGATURES = {"ﬀ": "ff", "ﬁ": "fi", "ﬂ": "fl", "ﬃ": "ffi", "ﬄ": "ffl", "ﬅ": "st", "ﬆ": "st"}


def page_text(doc, i: int) -> str:
    """Text of 0-based page i (pymupdf preferred, pypdf fallback). NUL bytes removed so grep treats files as text;
    typographic ligatures (ﬁ ﬂ ﬀ ﬃ ﬄ) expanded so that grep "fluid" finds "ﬂuid"."""
    try:
        text = doc[i].get_text("text")
    except TypeError:  # pypdf reader
        text = doc.pages[i].extract_text() or ""
    text = text.replace("\x00", "")
    for lig, rep in LIGATURES.items():
        text = text.replace(lig, rep)
    return text


def open_doc(pdf: Path):
    try:
        import fitz  # pymupdf
        return fitz.open(pdf), "pymupdf"
    except ImportError:
        from pypdf import PdfReader
        return PdfReader(str(pdf)), "pypdf"


def n_pages(doc) -> int:
    return doc.page_count if hasattr(doc, "page_count") else len(doc.pages)


def write_pdf(pdf: Path, start: int, end: int, dest: Path) -> None:
    from pypdf import PdfReader, PdfWriter
    reader, writer = PdfReader(str(pdf)), PdfWriter()
    for i in range(start - 1, end):
        writer.add_page(reader.pages[i])
    with open(dest, "wb") as f:
        writer.write(f)


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")   # Windows console: ligatures like 'ﬂ' in titles
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("ids", nargs="*", help="chapter ids (ch04); default all with do_chapter: true")
    ap.add_argument("--appendices", action="store_true", help="also split the appendices (reference material)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--check", action="store_true", help="assert each chapter opener page says CHAPTER <n>")
    a = ap.parse_args()

    book = load_book()
    pdf = locate_book_pdf(book)
    doc, engine = open_doc(pdf)
    total = n_pages(doc)
    off = int(book.get("book", {}).get("printed_offset", 0))
    print(f"Book PDF: {pdf.name}  ({total} pages, text engine: {engine})")

    entries = [c for c in book["chapters"] if c.get("do_chapter", True)]
    if a.appendices:
        entries += book.get("appendices", [])
    if a.ids:
        entries = [c for c in entries if c["id"] in a.ids]
        if not entries:
            sys.exit(f"ERROR: none of {a.ids} found in book.yaml")

    out = ROOT / "chapters"
    out.mkdir(exist_ok=True)
    page_map, problems = {}, []
    for c in entries:
        s, e = int(c["pdf_start"]), int(c["pdf_end"])
        if not (1 <= s <= e <= total):
            problems.append(f"{c['id']}: bad page range {s}-{e} (PDF has {total})")
            continue
        first = [ln.strip() for ln in page_text(doc, s - 1).splitlines() if ln.strip()][:3]
        print(f"  {c['id']}: pdf {s}-{e} ({e - s + 1} pp) | {c['title'][:40]:40s} | first lines: {first}")
        if a.check and "number" in c:
            opener = re.sub(r"\s+", "", "".join(first)).upper()     # the book sets it as "C H A P T E R / 7"
            if not re.search(rf"CHAPTER{c['number']}(?!\d)", opener):
                problems.append(f"{c['id']}: pdf page {s} does not look like the CHAPTER {c['number']} opener: {first}")
        page_map[c["id"]] = {"title": c["title"], "pdf_start": s, "pdf_end": e,
                             "printed_start": c.get("printed_start"), "printed_end": c.get("printed_end")}
        if a.dry_run or a.check:
            continue
        write_pdf(pdf, s, e, out / f"{c['id']}.pdf")
        with open(out / f"{c['id']}.txt", "w", encoding="utf-8") as f:
            for p in range(s, e + 1):
                f.write(f"\n\n===== PDF page {p} (printed page {p - off}) =====\n")
                f.write(page_text(doc, p - 1))

    if not (a.dry_run or a.check):
        (out / "_page_map.json").write_text(json.dumps({"pdf_pages": total, "printed_offset": off,
                                                        "chapters": page_map}, indent=2), encoding="utf-8")
        print(f"Wrote chapters/<id>.pdf, chapters/<id>.txt and chapters/_page_map.json for {len(page_map)} entries.")
    if problems:
        print("PROBLEMS:")
        for p in problems:
            print("  -", p)
        return 1
    if a.check:
        print("check OK: every chapter starts on its CHAPTER opener page")
    return 0


if __name__ == "__main__":
    sys.exit(main())
