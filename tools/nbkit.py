"""Building blocks for chapter notebooks (used by ``notebooks/build_chNN.py``), so every chapter has the same shape.

    import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "tools"))
    from nbkit import ChapterNotebook

    nb = ChapterNotebook("ch07")
    nb.title(big_idea="Waves on water are a balance between gravity pulling the surface flat and inertia overshooting.",
             roadmap=["7.2 linear waves", "7.5 group velocity", ...])
    nb.explainer_index([("dispersion_relation", "Why long waves outrun short ones"), ...])
    nb.setup()
    nb.section("7.2", "Linear liquid-surface gravity waves")   # one per book section; save() checks all are present
    nb.note("…", equation=r"c = \sqrt{gH}", ref="7.xx")  # NOTE tier;  nb.pointer("…") for a SKIP item
    nb.md("### The problem in plain words\n...")
    nb.code("...", explain="**What does the code above do?** ...")
    nb.figure_notes(see="...", read="...", change="...")
    nb.explainer("dispersion_relation", heading="Why long waves outrun short ones", why="...", tries=["...", "..."])
    nb.animation("...")                    # a FuncAnimation shown with fluidpy.core.anim.show_animation
    nb.live("...")                          # ipywidgets cell: tagged live-only (the published page shows a note)
    nb.summary(clicked=[...], feeds_forward=[...])
    nb.save()                               # notebooks/ch07_gravity_waves.ipynb

Cell tags written here and honoured by ``tools/publish_notebook.py``:
  ``setup``      the first code cell (clone on Colab / find the repo root locally)
  ``explainer``  a ``show_viz`` cell (its output becomes a full-window block on the published page)
  ``live-only``  output needs a running kernel (ipywidgets): the page replaces it with a "run it in Colab" note
  ``slow``       skipped in FAST runs by the cell's own ``if not FAST:`` guard (documentation only)
"""
from __future__ import annotations

import re
import textwrap
from pathlib import Path

import nbformat as nbf
import yaml

ROOT = Path(__file__).resolve().parents[1]
MAX_EXPLAINERS = 5


def _book() -> dict:
    return yaml.safe_load((ROOT / "book.yaml").read_text(encoding="utf-8"))


SETUP_MD = """## ⚙️ Setup — run this cell first

It works in three places without changes:
* **Google Colab** — clones the public repository into `/content/fluidpy` (the `fluidpy` package, its tests and the
  interactive explainers) and installs the one package Colab lacks (`pint`). No GPU is needed.
* **Jupyter / VS Code on your computer** — finds the repository root (the folder with `book.yaml`) and imports from it.
* **The published web page** — already executed; nothing to run.
"""

SETUP_CODE = '''# ── Setup: find (or fetch) the fluidpy repository, then load the house style ─────────────────────────────
import os, sys, subprocess, pathlib                     # standard library only, so this cell never fails on imports

IN_COLAB = "google.colab" in sys.modules                 # True when this notebook runs inside Google Colab
REPO_URL = "https://github.com/{repo}.git"               # the public repository with the fluidpy package

if IN_COLAB:
    ROOT = pathlib.Path("/content/fluidpy")              # where the repository is cloned on Colab
    if not ROOT.exists():                                 # first run in this Colab session: clone it
        subprocess.run(["git", "clone", "--depth", "1", REPO_URL, str(ROOT)], check=True)
    else:                                                 # later runs: fetch the newest version
        subprocess.run(["git", "-C", str(ROOT), "pull", "--ff-only", "-q"], check=False)
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", "-r", str(ROOT / "requirements-colab.txt")], check=True)
else:
    ROOT = pathlib.Path.cwd().resolve()                   # start where Jupyter was launched …
    while not (ROOT / "book.yaml").exists() and ROOT != ROOT.parent:
        ROOT = ROOT.parent                                # … and walk up until we find the repository root

os.chdir(ROOT)                                            # relative paths (outputs/, viz/) now resolve from the root
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))                         # make `import fluidpy` work

from fluidpy.core.style import setup_notebook            # matplotlib + plotly style shared with the explainers
from fluidpy.core.embed import show_viz                  # shows an interactive explainer in the notebook
from fluidpy.core.anim import show_animation             # plays a matplotlib animation inline
FAST = setup_notebook()                                   # FAST=True (env FLUIDPY_FAST=1) shrinks grids for quick runs
print(f"IN_COLAB = {{IN_COLAB}} | repository root: {{ROOT}} | FAST = {{FAST}}")
'''

LIVE_NOTE = ("> ▶️ **Live cell.** The widget below needs a running Python kernel. On the web page it is frozen — open this "
             "notebook in Colab (button at the top) or run it locally to drag the sliders.")


class ChapterNotebook:
    """Accumulates cells for one chapter notebook and writes ``notebooks/chNN_<slug>.ipynb``."""

    def __init__(self, chapter: str):
        book = _book()
        row = next((c for c in book["chapters"] if c["id"] == chapter), None)
        if row is None:
            raise KeyError(f"{chapter} not in book.yaml")
        self.chapter, self.row, self.book = chapter, row, book
        self.slug = row["slug"]
        self.nb = nbf.v4.new_notebook()
        self.cells = self.nb.cells
        self.explainers: list[str] = []
        self._has_setup = False
        self.covered: set[str] = set()

    # ---- primitives -------------------------------------------------------------------------------------------------
    def md(self, text: str, tags: list[str] | None = None) -> "ChapterNotebook":
        c = nbf.v4.new_markdown_cell(textwrap.dedent(text).strip("\n"))
        if tags:
            c.metadata["tags"] = tags
        self.cells.append(c)
        return self

    def code(self, src: str, explain: str | None = None, tags: list[str] | None = None) -> "ChapterNotebook":
        """A code cell (every line commented for a novice) optionally followed by its plain-words explanation."""
        c = nbf.v4.new_code_cell(textwrap.dedent(src).strip("\n"))
        if tags:
            c.metadata["tags"] = tags
        self.cells.append(c)
        if explain:
            body = explain.strip()
            if not body.startswith("**"):
                body = "**What does the code above do?**\n\n" + body
            self.md(body)
        return self

    # ---- structure ----------------------------------------------------------------------------------------------------
    def title(self, big_idea: str, roadmap: list[str], prerequisites: list[str] | None = None) -> "ChapterNotebook":
        r, b = self.row, self.book["book"]
        lines = [
            f"# Chapter {r['number']} — {r['title']}",
            "",
            f"*{b['title']}, {b['edition']} edition — {', '.join(b['authors'])} ({b['year']}).*",
            "This notebook teaches the chapter's key ideas with plain words, step-by-step mathematics, small worked "
            "examples, tested Python from the `fluidpy` package, figures, animations and interactive explainers. It does "
            "not reproduce the book's text or figures; equations are cited by their numbers so you can follow along in "
            "your copy.",
            "",
            "---",
            "",
            "**The big idea in one paragraph**",
            "",
            big_idea.strip(),
            "",
            "**Road map** — the ideas that, once they click, make the rest of the chapter easy:",
            "",
        ]
        lines += [f"{i + 1}. {s}" for i, s in enumerate(roadmap)]
        if prerequisites:
            lines += ["", "**What you should already know** (from earlier chapters):", ""] + [f"- {x}" for x in prerequisites]
        self.cells.append(nbf.v4.new_markdown_cell("\n".join(lines)))
        return self

    def explainer_index(self, items: list[tuple[str, str, str]]) -> "ChapterNotebook":
        """Table of the chapter's interactive explainers: (slug, title, what it makes clear)."""
        if len(items) > MAX_EXPLAINERS:
            raise ValueError(f"at most {MAX_EXPLAINERS} explainers per chapter (got {len(items)})")
        lines = ["## 🎮 Interactive explainers in this chapter", "", "| # | Explainer | What it makes clear |", "|---|---|---|"]
        lines += [f"| {i + 1} | **{t}** | {why} |" for i, (_s, t, why) in enumerate(items)]
        lines += ["", "Each one opens right where its idea is taught, fills the window (no scrolling), starts with a guided "
                  "**Walkthrough**, lets you **Explore** with sliders, shows the **Equations** with your numbers "
                  "substituted, and ends with **Check yourself** questions. On the web page they are full-window; in Colab "
                  "and Jupyter they appear in the cell output (use ⤢ *Full screen* or *Open in new tab* for more room)."]
        self.cells.append(nbf.v4.new_markdown_cell("\n".join(lines)))
        return self

    def setup(self) -> "ChapterNotebook":
        repo = self.book.get("project", {}).get("repo", "shammun/fluidpy")
        self.md(SETUP_MD)
        self.code(SETUP_CODE.format(repo=repo), tags=["setup"])
        self._has_setup = True
        return self

    def section(self, number: str, title: str, intro: str = "") -> "ChapterNotebook":
        """A notebook section for one book section ("7.2") or two merged short ones ("7.3–7.4" / "7.3, 7.4")."""
        self.covered.update(re.findall(r"\d+\.\d+", number))
        self.md(f"---\n\n## {number} {title}\n\n{textwrap.dedent(intro).strip()}")
        return self

    def note(self, text: str, equation: str | None = None, ref: str | None = None) -> "ChapterNotebook":
        """A NOTE-tier item: a short paragraph in our words, optionally the equation displayed with its book number."""
        body = textwrap.dedent(text).strip()
        if equation:
            tag = f"\\qquad \\text{{({ref})}}" if ref else ""
            body += f"\n\n$$ {equation.strip()} {tag} $$"
        self.cells.append(nbf.v4.new_markdown_cell(f"> 📝 **Note.** {body}" if "\n" not in body else f"📝 **Note.** {body}"))
        return self

    def pointer(self, text: str) -> "ChapterNotebook":
        """A SKIP-tier item: one line saying where the idea is covered instead."""
        self.cells.append(nbf.v4.new_markdown_cell(f"*↪ {textwrap.dedent(text).strip()}*"))
        return self

    def worked_example(self, title: str, steps_md: str) -> "ChapterNotebook":
        """A tiny example with easy numbers, traced step by step (before the general code)."""
        self.md(f"### ✏️ Tiny example: {title}\n\n{textwrap.dedent(steps_md).strip()}")
        return self

    def figure_notes(self, see: str, read: str, change: str) -> "ChapterNotebook":
        text = "\n\n".join([f"**What you see.** {textwrap.dedent(see).strip()}",
                            f"**How to read it.** {textwrap.dedent(read).strip()}",
                            f"**What would change if…** {textwrap.dedent(change).strip()}"])
        self.cells.append(nbf.v4.new_markdown_cell(text))
        return self

    def explainer(self, slug: str, heading: str, why: str, tries: list[str]) -> "ChapterNotebook":
        if slug in self.explainers:
            raise ValueError(f"explainer {slug} embedded twice")
        if len(self.explainers) >= MAX_EXPLAINERS:
            raise ValueError(f"at most {MAX_EXPLAINERS} explainers per chapter")
        self.explainers.append(slug)
        try_md = "\n".join(f"- {t}" for t in tries)
        self.md(f"### 🎮 Interactive: {heading}\n\n{textwrap.dedent(why).strip()}\n\n**What to try:**\n{try_md}")
        self.code(f'show_viz("{self.chapter}", "{slug}")   # full-width explainer; ⤢ Full screen for more room',
                  tags=["explainer"])
        return self

    def animation(self, src: str, explain: str | None = None) -> "ChapterNotebook":
        return self.code(src, explain=explain, tags=["animation"])

    def plotly(self, src: str, explain: str | None = None) -> "ChapterNotebook":
        return self.code(src, explain=explain, tags=["plotly"])

    def live(self, src: str, explain: str | None = None) -> "ChapterNotebook":
        self.md(LIVE_NOTE)
        return self.code(src, explain=explain, tags=["live-only"])

    def check_agree(self, src: str) -> "ChapterNotebook":
        """The from-scratch version next to the tested fluidpy function, with an assertion that both agree."""
        return self.code(src, tags=["from-scratch"])

    def summary(self, clicked: list[str], feeds_forward: list[str], left_out: list[str] | None = None) -> "ChapterNotebook":
        c = "\n".join(f"- {x}" for x in clicked)
        f = "\n".join(f"- {x}" for x in feeds_forward)
        lo = ("\n\n**Deliberately left out (and why that is safe)**\n" + "\n".join(f"- {x}" for x in left_out)) if left_out else ""
        self.md(f"---\n\n## ✅ What should have clicked\n\n{c}\n\n**What later chapters build on**\n\n{f}{lo}")
        return self

    # ---- output -------------------------------------------------------------------------------------------------------
    def missing_sections(self) -> list[str]:
        """Book sections (from book.yaml) that have no ``section(...)`` in this notebook yet."""
        numbers = [re.match(r"\s*(\d+\.\d+)", s).group(1) for s in self.row.get("sections", []) if re.match(r"\s*\d+\.\d+", s)]
        return [n for n in numbers if n not in self.covered]

    def save(self, allow_missing: list[str] | None = None) -> Path:
        """Write the notebook. Refuses if the setup cell is missing or a book section is not covered (every section of
        the chapter must appear — CORE/SUPPORT in depth, NOTE briefly, SKIP as a pointer line)."""
        if not self._has_setup:
            raise ValueError("call .setup() — every notebook needs the setup cell")
        missing = [n for n in self.missing_sections() if n not in (allow_missing or [])]
        if missing:
            raise ValueError(f"book sections with no nb.section(...): {missing} — every section must appear in the notebook")
        self.nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
        self.nb.metadata["language_info"] = {"name": "python"}
        self.nb.metadata["fluidpy"] = {"chapter": self.chapter, "explainers": self.explainers}
        for i, c in enumerate(self.cells):           # stable, readable cell ids
            c["id"] = f"{self.chapter}-{i:03d}"
        nbf.validate(self.nb)
        dest = ROOT / "notebooks" / f"{self.chapter}_{self.slug}.ipynb"
        dest.parent.mkdir(exist_ok=True)
        nbf.write(self.nb, dest)
        return dest


def notebook_path(chapter: str) -> Path:
    row = next(c for c in _book()["chapters"] if c["id"] == chapter)
    return ROOT / "notebooks" / f"{chapter}_{row['slug']}.ipynb"


def explainer_calls(nb_path: Path) -> list[str]:
    """``slug`` of every ``show_viz("chNN", "slug")`` call in a notebook, in order."""
    nb = nbf.read(nb_path, as_version=4)
    out = []
    for c in nb.cells:
        if c.cell_type == "code":
            out += [m.group(2) for m in re.finditer(r"show_viz\(\s*[\"'](\w+)[\"']\s*,\s*[\"']([\w\-]+)[\"']", c.source)]
    return out
