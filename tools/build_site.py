"""Generate the site's navigation pages from book.yaml + what exists on disk.

* ``index.html``      — home page: what the site is, how to use it, one card per chapter (Read · Colab · .ipynb ·
                         its explainers), mobile-first grid. Chapters that are not published yet show as "coming".
* ``viz/index.html``  — gallery of every interactive explainer, grouped by chapter, each opening full-window.
* ``README.md``       — the chapter table between ``<!-- INDEX_TABLE_START -->`` and ``<!-- INDEX_TABLE_END -->``.

Called by ``tools/publish_notebook.py``; run alone with ``.venv/Scripts/python.exe tools/build_site.py``.
"""
from __future__ import annotations

import html
import sys
import time
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fluidpy.core.embed import list_viz  # noqa: E402
from fluidpy.core.project import colab_url, config, pages_url, repo_url  # noqa: E402

E = html.escape


def chapter_rows() -> list[dict]:
    cfg = config(ROOT)
    rows = []
    for ch in cfg["chapters"]:
        if not ch.get("do_chapter", True):
            continue
        cid, slug = ch["id"], ch["slug"]
        name = f"{cid}_{slug}"
        nb = ROOT / "notebooks" / f"{name}.ipynb"
        rows.append({
            "id": cid, "number": ch["number"], "title": ch["title"], "slug": slug, "name": name,
            "blurb": ch.get("blurb", ""),
            "html": (ROOT / "notebooks" / f"{name}.html").exists(),
            "ipynb": nb.exists(),
            "colab": (ROOT / "notebooks" / f"{name}_colab.ipynb").exists(),
            "viz": list_viz(cid, ROOT),
        })
    return rows


def _head(title: str, prefix: str, desc: str) -> str:
    v = int(time.time())
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{E(title)}</title>
<meta name="description" content="{E(desc)}">
<link rel="stylesheet" href="{prefix}assets/clean-educational.css?v={v}">
<link rel="stylesheet" href="{prefix}assets/fluidpy-site.css?v={v}">
</head>
<body>
"""


def write_index(rows: list[dict]) -> Path:
    cfg = config(ROOT)
    proj, book = cfg["project"], cfg["book"]
    cards = []
    for r in rows:
        viz = "".join(f'<li><a href="viz/{r["id"]}/{E(m["slug"])}.html">{E(m["title"])}</a></li>' for m in r["viz"])
        if r["html"]:
            buttons = (f'<a class="fp-btn primary" href="notebooks/{r["name"]}.html">Read the chapter</a>'
                       + (f'<a class="fp-btn colab" href="{E(colab_url(r["name"] + "_colab", ROOT))}" target="_blank" rel="noopener">Open in Colab</a>' if r["colab"] else "")
                       + (f'<a class="fp-btn" href="notebooks/{r["name"]}.ipynb">.ipynb</a>' if r["ipynb"] else ""))
            status = f'<span class="fp-pill">{len(r["viz"])} interactive explainer{"s" if len(r["viz"]) != 1 else ""}</span>'
        else:
            buttons = ""
            status = '<span class="fp-pill todo">coming</span>'
        cards.append(f"""<article class="fp-ch{'' if r['html'] else ' pending'}" id="{r['id']}">
<div class="fp-ch-num">Chapter {r['number']} {status}</div>
<h3>{E(r['title'])}</h3>
<p>{E(r['blurb'])}</p>
{f'<ul>{viz}</ul>' if viz else ''}
<div class="fp-row">{buttons}</div>
</article>""")
    n_pub = sum(1 for r in rows if r["html"])
    page = _head(proj.get("site_title", "fluidpy"), "", "Kundu, Cohen & Dowling, Fluid Mechanics, learned through Python and interactive explainers") + f"""<div class="ce-container">
<header class="fp-hero">
<h1>{E(book.get('title', 'Fluid Mechanics'))} — learned through Python</h1>
<p class="lead">Every chapter of <em>{E(book.get('title', ''))}</em> ({E(book.get('edition', ''))} ed., {E(', '.join(book.get('authors', [])))},
{E(str(book.get('year', '')))}) turned into one learning package: a teaching notebook with plain-words explanations,
step-by-step mathematics, tested Python, figures and animations — plus up to five full-window interactive explainers
for the ideas that are easier to <em>feel</em> than to read.</p>
<div class="fp-actions"><a class="fp-btn primary" href="viz/index.html">&#127918; All interactive explainers</a>
<a class="fp-btn" href="{E(repo_url(ROOT))}">Source code on GitHub</a></div>
</header>
<section class="fp-how">
<div><b>&#128214; Read</b>Open a chapter page: explanations, code and figures, already run. The explainers fill your window right where their idea is taught.</div>
<div><b>&#127918; Play</b>Each explainer has a guided walkthrough, sliders, the equations with your numbers, and questions to check yourself. Works on a phone.</div>
<div><b>&#9654;&#65039; Run</b>Press <em>Open in Colab</em> to run and change every cell in your browser — nothing to install.</div>
<div><b>&#9989; Trust</b>Every function is tested against exact solutions, symbolic derivations, convergence studies or published benchmarks.</div>
</section>
<h2>Chapters <small style="font-size:.55em;color:var(--ce-text-faint)">{n_pub} of {len(rows)} published</small></h2>
<div class="fp-chapters">
{chr(10).join(cards)}
</div>
<p class="ce-source" style="margin-top:28px">The book is not part of this site or its repository: explanations are written in our own
words, figures are generated by our code, and equations are cited by number so you can follow along in your copy.
Code: MIT licence.</p>
</div>
</body>
</html>
"""
    out = ROOT / "index.html"
    out.write_text(page, encoding="utf-8")
    (ROOT / ".nojekyll").touch()
    return out


def write_gallery(rows: list[dict]) -> Path:
    sections = []
    for r in rows:
        if not r["viz"]:
            continue
        cards = []
        for i, m in enumerate(r["viz"]):
            in_ch = f'<a class="fp-btn" href="../notebooks/{r["name"]}.html#viz-{E(m["slug"])}">In the chapter</a>' if r["html"] else ""
            cards.append(f"""<div class="fp-gal">
<small>Ch. {r['number']} · explainer {i + 1} · §{E(m.get('sections', ''))}</small>
<b>{E(m['title'])}</b>
<p>{E(m.get('summary', ''))}</p>
<div class="fp-row" style="display:flex;gap:6px;flex-wrap:wrap"><a class="fp-btn primary" href="{r['id']}/{E(m['slug'])}.html">Open full window</a>{in_ch}</div>
</div>""")
        sections.append(f'<h2 id="{r["id"]}">Chapter {r["number"]} — {E(r["title"])}</h2>\n<div class="fp-gallery">{"".join(cards)}</div>')
    body = "\n".join(sections) or "<p>No explainers published yet.</p>"
    page = _head("Interactive explainers — fluidpy", "../", "Every interactive explainer, by chapter") + f"""<div class="ce-container">
<div class="fp-topbar"><div class="fp-crumbs"><a href="../index.html">&larr; All chapters</a></div></div>
<h1>Interactive explainers</h1>
<p>Each explainer fits your window with no scrolling, on a phone or a desktop. Start with the <b>Walkthrough</b> tab,
then <b>Explore</b> with the sliders, read the <b>Equations</b> with your own numbers substituted, and try the
<b>Check yourself</b> questions.</p>
{body}
</div>
</body>
</html>
"""
    out = ROOT / "viz" / "index.html"
    out.parent.mkdir(exist_ok=True)
    out.write_text(page, encoding="utf-8")
    return out


def write_readme(rows: list[dict]) -> Path:
    lines = ["| Chapter | Read (HTML) | Notebook | Colab | Explainers |", "|---|---|---|---|---|"]
    for r in rows:
        label = f"Ch. {r['number']}: {r['title']}"
        view = f"[View]({pages_url('notebooks/' + r['name'] + '.html', ROOT)})" if r["html"] else "—"
        ipynb = f"[.ipynb](notebooks/{r['name']}.ipynb)" if r["ipynb"] else "—"
        colab = f"[Open]({colab_url(r['name'] + '_colab', ROOT)})" if r["colab"] else "—"
        viz = ", ".join(f"[{m['title']}]({pages_url('viz/' + r['id'] + '/' + m['slug'] + '.html', ROOT)})" for m in r["viz"]) or "—"
        lines.append(f"| {label} | {view} | {ipynb} | {colab} | {viz} |")
    table = "\n".join(lines)
    readme = ROOT / "README.md"
    start, end = "<!-- INDEX_TABLE_START -->", "<!-- INDEX_TABLE_END -->"
    text = readme.read_text(encoding="utf-8") if readme.exists() else f"# fluidpy\n\n{start}\n{end}\n"
    if start not in text or end not in text:
        text += f"\n\n## Chapters\n\n{start}\n{end}\n"
    pre, rest = text.split(start, 1)
    _, post = rest.split(end, 1)
    readme.write_text(f"{pre}{start}\n{table}\n{end}{post}", encoding="utf-8")
    return readme


def build_all() -> list[Path]:
    rows = chapter_rows()
    return [write_index(rows), write_gallery(rows), write_readme(rows)]


if __name__ == "__main__":
    for p in build_all():
        print("wrote", p.relative_to(ROOT).as_posix())
