"""Publish one chapter: the executed notebook, its Colab twin, the web page, and the site navigation.

For chapter ``chNN`` (slug from ``book.yaml``) and its source notebook ``notebooks/chNN_<slug>.ipynb`` (written by
``notebooks/build_chNN.py``, no outputs) this produces:

``notebooks/chNN_<slug>.ipynb``        the SAME notebook, executed (public-safe run: no private book values), so it
                                      shows its figures on GitHub and in Jupyter. Explainer outputs become a static
                                      link card (GitHub's viewer strips iframes and scripts).
``notebooks/chNN_<slug>_colab.ipynb``  outputs stripped, Colab metadata; its setup cell clones the repository, and
                                      ``show_viz`` loads the explainers from GitHub Pages inside Colab.
``notebooks/chNN_<slug>.html``         the executed notebook rendered with nbconvert's lab template inside the site
                                      chrome: top bar (All chapters · Explainers · previous/next chapter), buttons
                                      (Open in Colab · Download .ipynb · GitHub), a strip of explainer cards, and every
                                      explainer as a FULL-WINDOW block (100 % width × 100 % height of the browser
                                      window, ``src="../viz/chNN/<slug>.html"``); live-only (ipywidgets) outputs become a
                                      "run it in Colab" note; a floating Contents button for phones.
``index.html``, ``viz/index.html``, ``README.md`` table   via ``tools/build_site.py``.

Gates (the tool refuses to publish when any fails): ``tools/embed_check.py`` (every explainer of the chapter exists,
is inlined with the current library, and is embedded exactly once), 0 execution errors, no private-run banner.

Usage (repo root)::

    .venv/Scripts/python.exe tools/publish_notebook.py ch07                # execute + publish one chapter
    .venv/Scripts/python.exe tools/publish_notebook.py ch07 --fast         # FLUIDPY_FAST=1 for the run
    .venv/Scripts/python.exe tools/publish_notebook.py --site-only         # only index.html, gallery, README
"""
from __future__ import annotations

import argparse
import copy
import html
import json
import os
import re
import sys
import time
from pathlib import Path

import nbformat

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from fluidpy.core.embed import MARKER_RE, list_viz, parse_meta  # noqa: E402
from fluidpy.core.project import colab_url, config, github_blob_url, pages_url, repo_url  # noqa: E402

E = html.escape
BANNER = re.compile(r"book values active", re.I)


# ---------------------------------------------------------------------------------------------------------------------
# execution
# ---------------------------------------------------------------------------------------------------------------------
KERNEL = "fluidpy-venv"


def ensure_kernel() -> str:
    """Register (once) a kernelspec that runs THIS interpreter, so execution never lands in another Python
    (Anaconda's ``python3`` kernel can shadow the venv's — lesson A5)."""
    from jupyter_client.kernelspec import KernelSpecManager, NoSuchKernel

    ksm = KernelSpecManager()
    try:
        spec = ksm.get_kernel_spec(KERNEL)
        if Path(spec.argv[0]).resolve() == Path(sys.executable).resolve():
            return KERNEL
    except NoSuchKernel:
        pass
    import subprocess

    subprocess.run([sys.executable, "-m", "ipykernel", "install", "--sys-prefix", "--name", KERNEL,
                    "--display-name", "Python (fluidpy .venv)"], check=True, capture_output=True)
    return KERNEL


def execute(nb_path: Path, *, fast: bool = False, timeout: int = 1800) -> tuple[nbformat.NotebookNode, float]:
    from nbconvert.preprocessors import ExecutePreprocessor

    nb = nbformat.read(nb_path, as_version=4)
    for c in nb.cells:                     # always start from a clean slate
        if c.cell_type == "code":
            c.outputs, c.execution_count = [], None
    old = os.environ.get("FLUIDPY_FAST")
    os.environ["FLUIDPY_FAST"] = "1" if fast else "0"
    t0 = time.perf_counter()
    try:
        ExecutePreprocessor(timeout=timeout, kernel_name=ensure_kernel()).preprocess(nb, {"metadata": {"path": str(ROOT)}})
    finally:
        if old is None:
            os.environ.pop("FLUIDPY_FAST", None)
        else:
            os.environ["FLUIDPY_FAST"] = old
    return nb, time.perf_counter() - t0


def execution_problems(nb) -> list[str]:
    out = []
    for i, c in enumerate(nb.cells):
        if c.cell_type != "code":
            continue
        for o in c.get("outputs", []):
            if o.get("output_type") == "error":
                out.append(f"cell {i}: {o.get('ename')}: {o.get('evalue')}")
            text = json.dumps(o)
            if BANNER.search(text):
                out.append(f"cell {i}: private-run banner 'book values active' in the output")
    return out


# ---------------------------------------------------------------------------------------------------------------------
# output transforms
# ---------------------------------------------------------------------------------------------------------------------
def _html_of(output) -> str | None:
    data = output.get("data", {}) if output.get("output_type") in ("display_data", "execute_result") else {}
    v = data.get("text/html")
    if v is None:
        return None
    return "".join(v) if isinstance(v, list) else v


def _viz_meta(key: str) -> dict:
    chapter, slug = key.split("/")
    p = ROOT / "viz" / chapter / f"{slug}.html"
    meta = parse_meta(p.read_text(encoding="utf-8")) if p.exists() else {}
    meta.setdefault("title", slug.replace("_", " "))
    return meta


def page_block(key: str, index: int, total: int) -> str:
    chapter, slug = key.split("/")
    m = _viz_meta(key)
    src = f"../viz/{chapter}/{slug}.html"
    return f"""<section class="fluidpy-viz fp-viz" id="viz-{E(slug)}" data-viz="{E(key)}">
<div class="fluidpy-viz-bar">
<span class="fp-viz-label">&#127918; <span class="fp-long">Interactive </span>{index}/{total}</span>
<span class="fp-viz-title">{E(m['title'])}</span>
<a class="fp-viz-btn fp-viz-open" href="{src}" target="_blank" rel="noopener">Open alone &#8599;</a>
<button type="button" class="fp-viz-btn fp-viz-fs">&#10530; Full screen</button>
</div>
<iframe src="{src}" title="{E(m['title'])}" loading="lazy" allow="fullscreen" allowfullscreen></iframe>
</section>"""


def link_card(key: str) -> str:
    chapter, slug = key.split("/")
    m = _viz_meta(key)
    url = pages_url(f"viz/{chapter}/{slug}.html", ROOT)
    return (f'<div style="border:1px solid #e3e6ee;border-radius:12px;padding:12px 14px;background:#f8f7ff">'
            f'<b>&#127918; Interactive explainer: {E(m["title"])}</b><br>{E(m.get("summary", ""))}<br>'
            f'<a href="{E(url)}" target="_blank" rel="noopener">Open it full-window &#8599;</a> '
            f'(re-run this cell to show it here)</div>')


LIVE_HTML = ('<div class="fp-live-note">&#9654;&#65039; <b>Live cell.</b> These sliders need a running Python kernel, so '
             'they are frozen on this page. Press <b>Open in Colab</b> at the top (or run the notebook locally) to use them.</div>')


def transform(nb, mode: str) -> tuple[nbformat.NotebookNode, list[str]]:
    """mode="page": explainer outputs → full-window blocks, live-only outputs → note.
    mode="ipynb": explainer outputs → static link cards (GitHub view); everything else unchanged."""
    nb = copy.deepcopy(nb)
    keys: list[str] = []
    for c in nb.cells:
        if c.cell_type != "code":
            continue
        tags = c.metadata.get("tags", [])
        new_outputs = []
        for o in c.get("outputs", []):
            h = _html_of(o)
            if h and MARKER_RE.search(h):
                for m in MARKER_RE.finditer(h):
                    keys.append(m.group("key"))
                    new_outputs.append(("viz", m.group("key")))
                continue
            new_outputs.append(("keep", o))
        outs = []
        for kind, val in new_outputs:
            if kind == "keep":
                outs.append(val)
            else:
                outs.append(nbformat.v4.new_output("display_data", data={"text/html": f"<!--VIZ:{val}-->", "text/plain": f"[explainer {val}]"}))
        if mode == "page" and "live-only" in tags:
            outs = [nbformat.v4.new_output("display_data", data={"text/html": LIVE_HTML, "text/plain": "[live cell]"})]
        c.outputs = outs
    total = len(keys)
    for c in nb.cells:
        if c.cell_type != "code":
            continue
        for o in c.outputs:
            h = _html_of(o)
            if h and h.startswith("<!--VIZ:"):
                key = h[len("<!--VIZ:"):-3]
                o["data"]["text/html"] = page_block(key, keys.index(key) + 1, total) if mode == "page" else link_card(key)
    return nb, keys


def colab_variant(src_nb, name: str) -> nbformat.NotebookNode:
    nb = copy.deepcopy(src_nb)
    for c in nb.cells:
        if c.cell_type == "code":
            c.outputs, c.execution_count = [], None
        c.metadata.pop("execution", None)
    nb.metadata["colab"] = {"provenance": [], "name": f"{name}_colab.ipynb", "toc_visible": True}
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3"}
    nb.metadata["language_info"] = {"name": "python"}
    nbformat.validate(nb)
    return nb


# ---------------------------------------------------------------------------------------------------------------------
# page
# ---------------------------------------------------------------------------------------------------------------------
PAGE_JS = r"""
<script id="fp-page-js">
(function () {
  // 1. move every explainer block out of its notebook cell so it can span the whole window
  document.querySelectorAll('section.fp-viz').forEach(function (el) {
    var cell = el.closest('.jp-Cell'); if (cell && cell.parentNode) cell.parentNode.insertBefore(el, cell.nextSibling);
  });
  // 2. full-bleed: exactly the viewport width, whatever the page layout around it
  function bleed() {
    var W = document.documentElement.clientWidth;
    document.querySelectorAll('section.fp-viz').forEach(function (el) {
      el.style.marginLeft = '0px'; el.style.width = W + 'px';
      el.style.marginLeft = (-el.getBoundingClientRect().left) + 'px';
    });
  }
  bleed(); window.addEventListener('resize', bleed); window.addEventListener('load', bleed);
  // 3. full-screen buttons (hidden where the browser does not allow it)
  var canFs = document.fullscreenEnabled || document.webkitFullscreenEnabled;
  document.querySelectorAll('.fp-viz-fs').forEach(function (b) {
    var f = b.closest('section').querySelector('iframe');
    if (!canFs) { b.style.display = 'none'; return; }
    b.onclick = function () { (f.requestFullscreen || f.webkitRequestFullscreen).call(f); };
  });
  // 4. floating contents (sections + explainers)
  var items = [];
  document.querySelectorAll('.jp-RenderedMarkdown h2, .jp-RenderedMarkdown h3, section.fp-viz').forEach(function (h) {
    if (h.tagName === 'SECTION') { items.push({ id: h.id, text: '🎮 ' + h.querySelector('.fp-viz-title').textContent, cls: 'viz' }); return; }
    if (!h.id) return;
    items.push({ id: h.id, text: h.textContent.replace(/¶$/, '').trim(), cls: h.tagName.toLowerCase() });
  });
  if (items.length > 2) {
    var btn = document.createElement('button'); btn.className = 'fp-toc-btn'; btn.type = 'button'; btn.textContent = '☰ Contents';
    var box = document.createElement('nav'); box.className = 'fp-toc'; box.setAttribute('aria-label', 'Contents');
    items.forEach(function (it) { var a = document.createElement('a'); a.href = '#' + it.id; a.textContent = it.text; a.className = it.cls; a.onclick = function () { box.classList.remove('open'); }; box.appendChild(a); });
    btn.onclick = function () { box.classList.toggle('open'); };
    document.body.appendChild(box); document.body.appendChild(btn);
    // never cover an explainer's own buttons: hide the floating button while an explainer fills the screen
    if ('IntersectionObserver' in window) {
      var inView = new Set();
      var io = new IntersectionObserver(function (entries) {
        entries.forEach(function (e) { if (e.intersectionRatio > 0.35) inView.add(e.target); else inView.delete(e.target); });
        var hide = inView.size > 0; btn.style.display = hide ? 'none' : ''; if (hide) box.classList.remove('open');
      }, { threshold: [0, 0.35, 0.6, 1] });
      document.querySelectorAll('section.fp-viz').forEach(function (s) { io.observe(s); });
    }
  }
})();
</script>
"""


def render_page(nb_page, *, row: dict, name: str, prev_row: dict | None, next_row: dict | None, title: str) -> str:
    from nbconvert import HTMLExporter

    exporter = HTMLExporter(template_name="lab")
    body, _ = exporter.from_notebook_node(nb_page)
    v = int(time.time())
    head_extra = (f'\n<link rel="stylesheet" href="../assets/clean-educational.css?v={v}">'
                  f'\n<link rel="stylesheet" href="../assets/fluidpy-site.css?v={v}">'
                  f'\n<meta name="description" content="{E(row.get("blurb", ""))}">'
                  "\n<style>.ce-container table code { white-space: pre-wrap; overflow-wrap: anywhere; }</style>")
    body = re.sub(r'(<meta charset="utf-8"\s*/?>)', lambda m: m.group(1) + head_extra, body, count=1)
    body = re.sub(r"<title>.*?</title>", f"<title>{E(title)}</title>", body, count=1, flags=re.S)

    def nav(r, label):
        if not r:
            return f"<span>{label}</span>"
        return f'<a href="{r["id"]}_{r["slug"]}.html" title="{E(r["title"])}">{label}</a>'

    viz = list_viz(row["id"], ROOT)
    cards = "".join(f'<a class="fp-viz-card" href="#viz-{E(m["slug"])}"><span class="fp-num">&#127918; {i + 1}</span>'
                    f'<b>{E(m["title"])}</b><span>{E(m.get("summary", ""))}</span></a>' for i, m in enumerate(viz))
    strip = (f'<div class="fp-viz-strip"><h2 id="fp-explainers">Interactive explainers in this chapter</h2>'
             f'<div class="fp-viz-cards">{cards}</div></div>') if viz else ""
    header = f"""
<div class="ce-container">
<div class="fp-topbar">
  <div class="fp-crumbs"><a href="../index.html">&larr; All chapters</a><a href="../viz/index.html">&#127918; All explainers</a></div>
  <div class="fp-chapnav">{nav(prev_row, "&larr; Previous")}{nav(next_row, "Next &rarr;")}</div>
</div>
<div class="fp-actions">
  <a class="fp-btn colab" href="{E(colab_url(name + '_colab', ROOT))}" target="_blank" rel="noopener">&#9654; Open in Colab</a>
  <a class="fp-btn" href="{name}.ipynb" download>&darr; Download .ipynb</a>
  <a class="fp-btn" href="{E(github_blob_url('notebooks/' + name + '.ipynb', ROOT))}" target="_blank" rel="noopener">View on GitHub</a>
</div>
{strip}
"""
    footer = f"""
<div class="fp-topbar" style="margin-top:40px">
  <div class="fp-crumbs"><a href="../index.html">&larr; All chapters</a></div>
  <div class="fp-chapnav">{nav(prev_row, "&larr; Previous")}{nav(next_row, "Next &rarr;")}</div>
</div>
<p class="ce-source">Code and explanations: <a href="{E(repo_url(ROOT))}">{E(config(ROOT)['project']['repo'])}</a> (MIT).
The book is not reproduced here; equations are cited by number.</p>
</div>
"""
    body = re.sub(r"(<body[^>]*>)", lambda m: m.group(1) + header, body, count=1)
    body = body.replace("</body>", footer + PAGE_JS + "\n</body>", 1)
    return body


def notebook_title(nb) -> str:
    for c in nb.cells:
        if c.cell_type == "markdown":
            for line in c.source.splitlines():
                if line.startswith("# "):
                    return line[2:].strip()
    return "fluidpy"


# ---------------------------------------------------------------------------------------------------------------------
def publish(chapter: str, *, fast: bool = False, skip_checks: bool = False) -> dict:
    cfg = config(ROOT)
    rows = [c for c in cfg["chapters"] if c.get("do_chapter", True)]
    idx = next(i for i, c in enumerate(rows) if c["id"] == chapter)
    row = rows[idx]
    name = f"{chapter}_{row['slug']}"
    src_path = ROOT / "notebooks" / f"{name}.ipynb"
    if not src_path.exists():
        raise SystemExit(f"{src_path.relative_to(ROOT)} not found — run notebooks/build_{chapter}.py first")

    if not skip_checks:
        import embed_check
        problems = embed_check.check(chapter)
        if problems:
            raise SystemExit("embed_check failed:\n  - " + "\n  - ".join(problems))

    print(f"== {name}: executing (FAST={fast}) …")
    executed, secs = execute(src_path, fast=fast)
    probs = execution_problems(executed)
    if probs:
        raise SystemExit("execution problems:\n  - " + "\n  - ".join(probs))
    print(f"   executed in {secs:.0f} s")

    source = copy.deepcopy(executed)
    for c in source.cells:
        if c.cell_type == "code":
            c.outputs, c.execution_count = [], None
    colab = colab_variant(source, name)
    nbformat.write(colab, ROOT / "notebooks" / f"{name}_colab.ipynb")

    page_nb, keys = transform(executed, "page")
    ipynb_nb, _ = transform(executed, "ipynb")
    nbformat.write(ipynb_nb, src_path)

    prev_row = rows[idx - 1] if idx > 0 else None
    next_row = rows[idx + 1] if idx + 1 < len(rows) else None
    page = render_page(page_nb, row=row, name=name, prev_row=prev_row, next_row=next_row,
                       title=f"Ch. {row['number']} {row['title']} — fluidpy")
    out = ROOT / "notebooks" / f"{name}.html"
    out.write_text(page, encoding="utf-8")

    import build_site as _site
    _site.build_all()
    rep = {"chapter": chapter, "seconds": round(secs, 1), "explainers": keys,
           "html_mb": round(out.stat().st_size / 1e6, 2), "ipynb_mb": round(src_path.stat().st_size / 1e6, 2),
           "page_url": pages_url(f"notebooks/{name}.html", ROOT), "colab_url": colab_url(f"{name}_colab", ROOT)}
    print(json.dumps(rep, indent=2))
    return rep


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("chapter", nargs="?")
    ap.add_argument("--fast", action="store_true", help="run with FLUIDPY_FAST=1")
    ap.add_argument("--site-only", action="store_true")
    ap.add_argument("--skip-checks", action="store_true", help="(debug only) skip embed_check")
    a = ap.parse_args(argv)
    if a.site_only:
        import build_site as _site
        for p in _site.build_all():
            print("wrote", p.relative_to(ROOT).as_posix())
        return 0
    if not a.chapter:
        ap.error("give a chapter id (ch07) or --site-only")
    publish(a.chapter, fast=a.fast, skip_checks=a.skip_checks)
    return 0


if __name__ == "__main__":
    sys.exit(main())
