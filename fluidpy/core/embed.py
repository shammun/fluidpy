"""Show an interactive explainer (``viz/chNN/<slug>.html``) inside a notebook, a Colab cell and the published page.

One call in a notebook code cell::

    from fluidpy.core.embed import show_viz
    show_viz("ch07", "dispersion_relation")

Where it is shown, and how
--------------------------
=====================  ====================================================================================
host                   what ``mode="auto"`` emits
=====================  ====================================================================================
local Jupyter/VS Code  ``<iframe srcdoc="…">`` with the whole self-contained explainer from the repo on disk
                       (works offline and shows your local edits immediately)
Google Colab           ``<iframe src="https://<user>.github.io/fluidpy/viz/chNN/<slug>.html">`` — the proven
                       pattern from the fast.ai notes — when that URL answers 200; otherwise ``srcdoc`` from the
                       repository clone in ``/content/fluidpy`` (so an unpublished explainer still works)
published HTML page    ``tools/publish_notebook.py`` finds the ``fluidpy-viz`` markers and replaces the whole
                       output with a full-bleed block: 100 % of the window width and height, loading
                       ``../viz/chNN/<slug>.html`` (same origin, lazy)
=====================  ====================================================================================

Size: the iframe is as wide as the output area and as tall as the visible window — the browser window when the
output lives in the top-level page (JupyterLab, Notebook 7), or ``screen.availHeight − 170`` (clamped to 520…1000 px)
inside sandboxed output frames (Colab, VS Code), whose own height grows with their content. The explainers themselves
are built to fit any box without scrolling (``assets/viz_lib.js``), so whatever height we give them is used fully.
A slim toolbar offers **Full screen** (only when the host allows it) and **Open in new tab**.

If the explainer file does not exist yet (explainers and the notebook are built in parallel), a clearly marked
placeholder is shown instead of raising, and ``tools/embed_check.py`` refuses to publish until the file exists.
"""
from __future__ import annotations

import html
import re
import sys
import urllib.request
import uuid
from functools import lru_cache
from html.parser import HTMLParser
from pathlib import Path

from .project import pages_url, repo_root

MIN_HEIGHT = 520
MAX_FRAME_HEIGHT = 1000
BEGIN = "<!-- fluidpy-viz:BEGIN {key} -->"
END = "<!-- fluidpy-viz:END -->"
MARKER_RE = re.compile(r"<!-- fluidpy-viz:BEGIN (?P<key>[\w\-]+/[\w\-]+) -->.*?<!-- fluidpy-viz:END -->", re.S)


# ---------------------------------------------------------------------------------------------------------------------
# metadata
# ---------------------------------------------------------------------------------------------------------------------
class _MetaParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.meta: dict[str, str] = {}
        self._in_title = False
        self.title = ""

    def handle_starttag(self, tag, attrs):
        if tag == "meta":
            a = dict(attrs)
            name = a.get("name") or ""
            if name.startswith("viz:"):
                self.meta[name[4:]] = a.get("content") or ""
        elif tag == "title":
            self._in_title = True

    def handle_endtag(self, tag):
        if tag == "title":
            self._in_title = False
        if tag == "head":
            raise _StopParsing

    def handle_data(self, data):
        if self._in_title:
            self.title += data


class _StopParsing(Exception):
    pass


def parse_meta(text: str) -> dict[str, str]:
    """The ``<meta name="viz:*">`` tags (and ``<title>``) of an explainer's HTML text."""
    p = _MetaParser()
    try:
        p.feed(text)
    except _StopParsing:
        pass
    meta = dict(p.meta)
    if "title" not in meta and p.title.strip():
        meta["title"] = p.title.strip()
    return meta


# ---------------------------------------------------------------------------------------------------------------------
# titles with a little maths: "Where does $6\pi\mu a U$ come from?"
# ---------------------------------------------------------------------------------------------------------------------
# Explainer titles and summaries may carry inline TeX between $…$ (the explainer itself renders it with KaTeX). The
# site index, the gallery, the README, the notebook toolbar and the page blocks have no maths renderer, so they use
# these two helpers: ``title_html`` turns the small TeX subset used in titles into styled HTML (italic letters,
# <sup>/<sub>, Greek letters, a real minus sign) and ``title_text`` gives a clean plain-text form (browser tabs,
# iframe titles, Markdown). Anything outside the subset is shown literally, never dropped.
_TEX_SYMBOLS = {
    "alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "epsilon": "ε", "varepsilon": "ε", "zeta": "ζ",
    "eta": "η", "theta": "θ", "kappa": "κ", "lambda": "λ", "mu": "μ", "nu": "ν", "xi": "ξ", "pi": "π", "rho": "ρ",
    "sigma": "σ", "tau": "τ", "phi": "φ", "varphi": "φ", "chi": "χ", "psi": "ψ", "omega": "ω",
    "Gamma": "Γ", "Delta": "Δ", "Theta": "Θ", "Lambda": "Λ", "Pi": "Π", "Sigma": "Σ", "Phi": "Φ", "Psi": "Ψ",
    "Omega": "Ω", "nabla": "∇", "partial": "∂", "infty": "∞", "cdot": "·", "times": "×", "approx": "≈",
    "sim": "~", "le": "≤", "leq": "≤", "ge": "≥", "geq": "≥", "ne": "≠", "neq": "≠", "to": "→", "sqrt": "√",
    "pm": "±", ",": " ", ";": " ", "!": "", " ": " ",
}
_MATH_FONT = "font-family:'KaTeX_Main','Cambria Math','Times New Roman',serif;font-size:1.08em;white-space:nowrap"
_SCRIPT = re.compile(r"([\^_])(\{([^{}]*)\}|\\[A-Za-z]+|.)")


def _group(x: str) -> str:
    return x if len(x) <= 1 else f"({x})"


def _tex_structures(tex: str) -> str:
    """\\frac{a}{b} → a/b and \\sqrt{a} → √(a) (innermost first), then drop grouping braces not owned by ^ or _."""
    for _ in range(4):
        new = re.sub(r"\\frac\{([^{}]*)\}\{([^{}]*)\}", lambda m: f"{_group(m.group(1))}/{_group(m.group(2))}", tex)
        new = re.sub(r"\\sqrt\{([^{}]*)\}", lambda m: f"√{_group(m.group(1))}", new)
        if new == tex:
            break
        tex = new
    return re.sub(r"(?<![\^_])\{([^{}]*)\}", r"\1", tex)


_RELATIONS = "=≈<>≤≥≠~→×"
_SIGNS = "−+±"


def _tex_symbols(s: str, after_operand: bool = False) -> str:
    """Replace \\command by its symbol; "-" by a real minus; drop TeX's ignored spaces; space binary operators.

    Relations (=, ≈, …) are always spaced. A sign (−, +, ±) is spaced when it follows an operand (a − b) and left
    tight when unary (−b); ``after_operand`` says whether the text just before ``s`` (e.g. a superscript) was one."""
    s = re.sub(r"\\([A-Za-z]+|[,;! ])", lambda m: _TEX_SYMBOLS.get(m.group(1), m.group(0)), s)
    s = s.replace("-", "−").replace(" ", "")
    out, prev_operand = [], after_operand
    for ch in s:
        if ch in _RELATIONS or (ch in _SIGNS and prev_operand):
            out.append(f" {ch} ")
            prev_operand = False
        else:
            out.append(ch)
            prev_operand = ch not in _SIGNS and ch not in "(["
    return "".join(out)


def _math_html(tex: str) -> str:
    tex = _tex_structures(tex)
    out: list[str] = []
    pos = 0
    for m in _SCRIPT.finditer(tex):
        out.append(_letters_html(tex[pos:m.start()], after_operand=pos > 0))
        body = m.group(3) if m.group(3) is not None else m.group(2)
        tag = "sup" if m.group(1) == "^" else "sub"
        gap = ' style="margin-left:.1em"' if tag == "sup" else ""  # italic correction: an italic base leans into it
        out.append(f"<{tag}{gap}>{_letters_html(body)}</{tag}>")
        pos = m.end()
    out.append(_letters_html(tex[pos:], after_operand=pos > 0))
    return f'<span class="fp-math" style="{_MATH_FONT}">{"".join(out)}</span>'


def _letters_html(tex: str, after_operand: bool = False) -> str:
    s = _tex_symbols(tex, after_operand)
    parts = []
    for ch in s:
        e = html.escape(ch)
        # italic for Latin and lower-case Greek letters (as KaTeX sets them), upright for digits, capitals Greek, signs
        parts.append(f"<i>{e}</i>" if (ch.isascii() and ch.isalpha()) or ("α" <= ch <= "ω") else e)
    return "".join(parts)


def title_html(title: str) -> str:
    """HTML for a title that may contain inline $…$ TeX (escaped text outside the maths)."""
    pieces = re.split(r"(?<!\\)\$", title or "")
    if len(pieces) % 2 == 0:  # unbalanced $: show it literally
        return html.escape(title)
    return "".join(_math_html(p) if i % 2 else html.escape(p) for i, p in enumerate(pieces))


def title_text(title: str) -> str:
    """Plain text for a title that may contain inline $…$ TeX: t^{1/5} → t^(1/5), 6\\pi\\mu a U → 6πμaU."""
    pieces = re.split(r"(?<!\\)\$", title or "")
    if len(pieces) % 2 == 0:
        return title

    def plain(tex: str) -> str:
        tex = _tex_structures(tex)
        tex = _SCRIPT.sub(lambda m: m.group(1) + (f"({m.group(3)})" if m.group(3) is not None and len(m.group(3)) > 1
                                                   else (m.group(3) if m.group(3) is not None else m.group(2))), tex)
        return _tex_symbols(tex)

    return "".join(plain(p) if i % 2 else p for i, p in enumerate(pieces))


def viz_meta(path: str | Path) -> dict[str, str]:
    """Metadata of one explainer file: chapter, slug, order, title, summary, concept, sections, equations, fluidpy."""
    path = Path(path)
    meta = parse_meta(path.read_text(encoding="utf-8", errors="replace"))
    meta.setdefault("chapter", path.parent.name)
    meta.setdefault("slug", path.stem)
    meta.setdefault("title", path.stem.replace("_", " ").capitalize())
    meta["path"] = path.as_posix()
    return meta


def list_viz(chapter: str, root: str | Path | None = None) -> list[dict[str, str]]:
    """All explainers of a chapter (``viz/<chapter>/*.html`` except ``index.html``), sorted by ``viz:order`` then slug."""
    base = Path(root) if root else repo_root()
    folder = base / "viz" / chapter
    if not folder.is_dir():
        return []
    items = [viz_meta(p) for p in sorted(folder.glob("*.html")) if p.name != "index.html"]

    def order(m: dict) -> tuple:
        try:
            return (float(m.get("order") or "inf"), m["slug"])
        except ValueError:
            return (float("inf"), m["slug"])

    return sorted(items, key=order)


# ---------------------------------------------------------------------------------------------------------------------
# environment
# ---------------------------------------------------------------------------------------------------------------------
def in_colab() -> bool:
    return "google.colab" in sys.modules


@lru_cache(maxsize=64)
def _url_ok(url: str, timeout: float = 6.0) -> bool:
    try:
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "fluidpy"})
        with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310 - fixed https URL
            return 200 <= r.status < 300
    except Exception:  # noqa: BLE001
        return False


def _local_file(chapter: str, slug: str, root: str | Path | None) -> Path | None:
    rel = Path("viz") / chapter / f"{slug}.html"
    base = Path(root) if root else repo_root()
    for c in (base / rel, Path.cwd() / rel, Path("/content/fluidpy") / rel):
        if c.is_file():
            return c
    return None


# ---------------------------------------------------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------------------------------------------------
_BAR = ("display:flex;align-items:center;gap:8px;flex-wrap:wrap;padding:6px 10px;"
        "font:500 13px/1.3 Inter,-apple-system,'Segoe UI',Roboto,sans-serif;"
        "background:linear-gradient(90deg,#f3f1fd,#eefaf9);color:#1a1f36;border:1px solid #e3e6ee;border-bottom:0;"
        "border-radius:12px 12px 0 0")
_BTN = ("font:inherit;font-size:12px;padding:3px 10px;border-radius:999px;border:1px solid #6c5ce7;background:#fff;"
        "color:#4f3dd0;cursor:pointer;text-decoration:none;white-space:nowrap")


def viz_html(chapter: str, slug: str, *, mode: str = "auto", height: int | None = None, title: str | None = None,
             root: str | Path | None = None) -> str:
    """The HTML that :func:`show_viz` displays (also used by tests and by the publish tool).

    Parameters
    ----------
    chapter : "ch07"-style chapter id.
    slug : explainer file name without ``.html``.
    mode : "auto" (see the module docstring), "srcdoc" (inline the file) or "src" (load the GitHub Pages copy).
    height : fixed iframe height in px; by default the iframe fills the visible window height.
    title : toolbar title; default is the explainer's ``viz:title`` meta tag.
    """
    if not re.fullmatch(r"[A-Za-z0-9_\-]+", chapter) or not re.fullmatch(r"[A-Za-z0-9_\-]+", slug):
        raise ValueError(f"bad chapter/slug: {chapter!r}/{slug!r}")
    key = f"{chapter}/{slug}"
    rel = f"viz/{chapter}/{slug}.html"
    page = pages_url(rel, root=root)
    local = _local_file(chapter, slug, root)

    if mode == "auto":
        mode = "src" if (in_colab() and _url_ok(page)) else "srcdoc"
    if mode == "srcdoc" and local is None:
        if in_colab() and _url_ok(page):
            mode = "src"
        else:
            return _placeholder(key, rel, page)

    text = local.read_text(encoding="utf-8") if local is not None else ""
    meta = parse_meta(text) if text else {}
    title = title or meta.get("title") or slug.replace("_", " ")
    uid = f"fpviz-{uuid.uuid4().hex[:10]}"
    fixed = "true" if height else "false"
    h0 = int(height) if height else 760
    if mode == "src":
        source = f'src="{html.escape(page)}"'
    else:
        source = f'srcdoc="{html.escape(text, quote=True)}"'
    return f"""{BEGIN.format(key=key)}
<div class="fluidpy-viz" id="{uid}" data-viz="{key}" data-mode="{mode}" style="width:100%;margin:10px 0 18px">
<div class="fluidpy-viz-bar" style="{_BAR}">
<span style="font-weight:650">&#127918; Interactive explainer</span>
<span class="fluidpy-viz-title" style="flex:1 1 auto;min-width:8em;color:#3c4257">{title_html(title)}</span>
<button type="button" class="fluidpy-viz-fs" style="{_BTN}">&#10530; Full screen</button>
<a class="fluidpy-viz-open" href="{html.escape(page)}" target="_blank" rel="noopener" style="{_BTN}">Open in new tab &#8599;</a>
</div>
<iframe class="fluidpy-viz-frame" title="{html.escape(title_text(title))}" {source} allow="fullscreen" allowfullscreen
 style="width:100%;height:{h0}px;border:1px solid #e3e6ee;border-radius:0 0 12px 12px;display:block;background:#f5f6fb"></iframe>
<script>
(function () {{
  var w = document.getElementById("{uid}"); if (!w) return;
  var f = w.querySelector("iframe"), bar = w.querySelector(".fluidpy-viz-bar"), fs = w.querySelector(".fluidpy-viz-fs");
  var FIXED = {fixed}, isTop = false;
  try {{ isTop = (window.top === window); }} catch (e) {{ isTop = false; }}
  function size() {{
    if (FIXED) return;
    var H = isTop ? window.innerHeight : Math.min(Math.max(screen.availHeight - 170, {MIN_HEIGHT}), {MAX_FRAME_HEIGHT});
    f.style.height = Math.max({MIN_HEIGHT}, Math.round(H - (bar ? bar.offsetHeight : 0) - 4)) + "px";
  }}
  size(); if (isTop) window.addEventListener("resize", size);
  var can = (document.fullscreenEnabled || document.webkitFullscreenEnabled) && (f.requestFullscreen || f.webkitRequestFullscreen);
  if (!can) {{ fs.style.display = "none"; }}
  else {{ fs.onclick = function () {{ (f.requestFullscreen || f.webkitRequestFullscreen).call(f); }}; }}
}})();
</script>
</div>
{END}"""


def _placeholder(key: str, rel: str, page: str) -> str:
    return f"""{BEGIN.format(key=key)}
<div class="fluidpy-viz fluidpy-viz-missing" data-viz="{key}" style="margin:10px 0 18px;padding:14px 16px;border:1px dashed #d97706;
 border-radius:12px;background:#fffbeb;font:14px/1.5 Inter,-apple-system,'Segoe UI',sans-serif;color:#422006">
&#9888;&#65039; Interactive explainer <code>{html.escape(rel)}</code> is not in this copy of the repository yet.
Once it is published it lives at <a href="{html.escape(page)}" target="_blank" rel="noopener">{html.escape(page)}</a>.
</div>
{END}"""


def show_viz(chapter: str, slug: str, *, mode: str = "auto", height: int | None = None, title: str | None = None,
             root: str | Path | None = None) -> None:
    """Display the explainer ``viz/<chapter>/<slug>.html`` in the current notebook output (see module docstring).

    Displays the output itself and returns ``None`` (so it is never shown twice).
    """
    from IPython.display import HTML, display

    display(HTML(viz_html(chapter, slug, mode=mode, height=height, title=title, root=root)))


def embedded_keys(text: str) -> list[str]:
    """``chNN/slug`` keys of every explainer block inside an HTML page or notebook JSON text."""
    return [m.group("key") for m in MARKER_RE.finditer(text)]


__all__ = ["show_viz", "viz_html", "viz_meta", "list_viz", "parse_meta", "title_html", "title_text", "embedded_keys", "in_colab",
           "BEGIN", "END", "MARKER_RE"]
