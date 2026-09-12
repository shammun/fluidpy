# fluidpy — Kundu, Cohen & Dowling, *Fluid Mechanics* (5th ed., 2012), learned through Python and interactive explainers

## Mission
Turn each chapter of the book into **one learning package** that makes the chapter *click* for a first-time reader
(and for Shammunul, who is working through the book for a PhD in climate dynamics):

| Per chapter | What it is |
|---|---|
| `notebooks/chNN_<slug>.ipynb` | The teaching notebook, **executed**: plain words → step-by-step maths → tiny worked example → commented Python → figure → how to read it; Python animations, plotly slider figures, live widgets; the chapter's explainers embedded |
| `notebooks/chNN_<slug>_colab.ipynb` | The same notebook for Google Colab (outputs stripped; its setup cell clones this repo; explainers load from GitHub Pages) |
| `notebooks/chNN_<slug>.html` | The published page: executed notebook, every explainer a **full-window block** (100 % width × 100 % height), Open-in-Colab button, prev/next, contents |
| `viz/chNN/<slug>.html` (**≤ 5**) | Self-contained vanilla-JS explainers: phenomenon + controls + guided walkthrough + live equations + quiz, fitting any window with **no scrolling** |

Plus the machinery that makes it trustworthy and cumulative: `fluidpy/` (the physics as tested functions), `tests/`,
`reports/`, `knowledge/`. Site: `index.html` (chapters) and `viz/index.html` (explainer gallery), served at
`https://shammun.github.io/fluidpy/` (all URLs derive from `book.yaml → project`).

**Clicking beats coverage.** A chapter has a handful of ideas that, once understood, make dozens of smaller results
obvious. The `concept-curator` picks them and writes down what it leaves out and why that is safe.

## One command
`/do-chapter N` runs the whole pipeline for chapter N (10 phases, parallel where safe, resumable with `--from PHASE`).
`/status` says where things are and what to run next. `GUIDE.md` explains everything else.

| # | Phase | Agent(s) | Writes |
|---|---|---|---|
| 1 | analyze | `concept-analyst` | `analysis/chNN.md` |
| 2 | curate | `concept-curator` | `analysis/chNN_curation.md` |
| 3 | design ∥ 4 implement | `lesson-designer` ∥ `concept-implementer` | `analysis/chNN_design.md` ∥ `fluidpy/`, `scripts/` |
| 5 | verify | `math-verifier` (⟲ implementer, ≤ 3 loops) | `tests/`, `reference/`, `reports/chNN_verification.md` |
| 6 | review ∥ 7 viz ∥ 8 notebook | `derivation-reviewer` ∥ `viz-builder` ×≤5 → `viz-reviewer` ∥ `notebook-builder` | `reports/chNN_review.md` ∥ `viz/chNN/`, `reports/chNN_viz.md` ∥ `notebooks/` |
| 9 | knowledge ∥ 10 publish | `knowledge-keeper` ∥ `site-publisher` | `knowledge/` ∥ pages, Colab twin, index, gallery |

## Where things live
| Path | What | Committed? |
|---|---|---|
| `*.pdf` (repo root) | the book — read-only, never modified | **never** |
| `book.yaml` | project config + chapter map (exact PDF pages from the outline, sections, explainer seeds) | yes |
| `chapters/chNN.{pdf,txt}`, `chapters/pages/` | split chapter, extracted text, rendered page images | **never** |
| `analysis/chNN.md`, `chNN_curation.md`, `chNN_design.md` | inventory · teaching spine + ≤5 explainer picks · storyboards | yes |
| `fluidpy/core/` | primitives: `project`, `embed` (show_viz), `anim`, `interact` (plotly sliders), `style`, `units`, `refdata`, + physics reused by ≥2 chapters | yes |
| `fluidpy/chNN_<slug>.py`, `scripts/chNN_*.py` | chapter physics + runnable demos | yes |
| `viz/chNN/<slug>.html` | explainers (library inlined from `assets/viz_lib.js` + `assets/viz_base.css` by `tools/viz_inline.py`) | yes |
| `tests/test_chNN.py`, `reference/chNN/` | evidence (+ cited public benchmark data) | yes |
| `reports/chNN_{verification,review,viz}.md` | verdicts | yes |
| `reports/viz/**` | screenshots + audit JSON from `tools/shot.py` | no |
| `notebooks/build_chNN.py` → `chNN_<slug>.ipynb` (+ `_colab.ipynb`, `.html`) | builder (uses `tools/nbkit.py`) and its products | yes |
| `knowledge/` | `CUMULATIVE.md`, `concept_map.md`, `notation.md`, `viz_patterns.md`, `chNN.md` | yes |
| `progress.json` | phase status per chapter | yes |
| `outputs/` | figures, GIFs, self-test pages | no |

## Non-negotiable rules
1. **Sequential chapters.** Start chapter N only when chapter N−1 has `publish: pass` and `knowledge: done` (or the
   user explicitly overrides). Every agent reads `knowledge/CUMULATIVE.md`, `knowledge/notation.md` and
   `knowledge/viz_patterns.md` before chapter work.
2. **Equations come from page images.** `chapters/chNN.txt` garbles maths (`¼` is `=`, `ð…Þ` are parentheses, a minus
   vanishes, ω prints as `u`). Render the defining page (`.venv/Scripts/python.exe tools/render_pages.py chNN --eq 7.27`)
   and read the PNG before transcribing, coding or displaying an equation. Never code an equation from memory.
3. **Teach in Shammunul's style** (skill `teaching-style`): plain words first; *problem → idea → maths → code → what the
   output shows*; a tiny example with easy numbers traced step by step; every code line commented for a novice; every
   figure followed by *what you see / how to read it / what would change if…*; "What does the code above do?" blocks.
4. **At most 5 interactive explainers per chapter**, each attached to a CORE concept and justified by "why interaction
   beats a static figure here". Each fits the window with **no scrolling** at 360×640 … 1920×1080 and in a notebook
   frame, shows and explains its equations, has a 4–8-step walkthrough. Content that does not fit becomes a tab or a
   page — never a scrollbar, never text below 12 px. Enforced by `tools/viz_lint.py` + `tools/shot.py` and judged by eye
   by the `viz-reviewer`.
5. **Traceable physics.** Every `fluidpy` function's docstring cites book section + equation number(s), lists symbols
   with SI units and assumptions, and carries its validation label. Explainer JS physics mirrors a `fluidpy` function
   and proves it with `selftest()` parity rows (`py:` expressions evaluated by `tools/shot.py`).
6. **Evidence, proportional to importance** (skill `verify-implementation`): CORE items ≥ 2 independent evidence levels
   (one of V1 analytic / V2 symbolic / V3 convergence / V5 benchmark); SUPPORT ≥ 1. Never loosen a tolerance to pass;
   fix the physics or report the discrepancy. Max 3 fix loops, then stop and tell the user.
7. **Subagents do the work; the main session orchestrates** and keeps its context small. Agents never commit — the
   orchestrator commits after each phase (`chNN: <phase> — <one line>`). Never run two agents that write the same
   folder at the same time (safe parallel groups are the ∥ rows above; viz-builders each write one file).
8. **Reuse before re-implementing.** Check `knowledge/concept_map.md`, `fluidpy/core/`, `assets/viz_lib.js`,
   `knowledge/viz_patterns.md` first. Anything used by two chapters moves to `core/` (Python) or `viz_lib.js` (JS; only
   the orchestrator or knowledge-keeper edits the library, then `tools/viz_inline.py --all`).
9. **The repo is public; the book is not.** Never commit or publish book text, page images, crops or scans of book
   figures, transcribed tables, or exercise text. Equations may be shown (cited by number); explanations are in our own
   words; figures are generated by our code. Book-quoted numbers used as test evidence stay in git-ignored
   `tests/book_values_chNN.json`. `tools/check_public.py` runs before every push (`.githooks/pre-push`). Never use
   Colab's *Save a copy in GitHub*.
10. **Windows host, path with spaces.** Quote paths; `pathlib` in Python; run Python as `.venv/Scripts/python.exe`; use
    `python -m nbconvert` (never `jupyter nbconvert`, which can dispatch to Anaconda); notebooks execute on the
    `fluidpy-venv` kernel (`tools/publish_notebook.py` registers it). Long commands go in `tools/*.py`, not one-liners.
    Write code files with the Write/Edit tools — shell heredocs mangle backslashes on this machine.
11. **Outward-facing actions** (`gh repo create`, first push, enabling Pages) — ask the user, or have them run the exact
    command with a leading `!` if the permission classifier blocks it. Never force-push.
12. **Never hand back placeholders.** No agent returns while a run is still going; report real numbers.

## Python stack
numpy · scipy · sympy · matplotlib · plotly · ipywidgets · pandas · pint · pytest · nbformat/nbconvert · pymupdf/pypdf
· pyyaml · playwright (browser audits via the installed Edge/Chrome) · imageio-ffmpeg (MP4 animations).
Explainers: vanilla JS + `assets/viz_lib.js` + KaTeX (CDN, with text fallback). No build step anywhere.

## Commands (skills with `/`)
`/setup-project` (once) · `/do-chapter N [--from PHASE]` · `/status` · `/build-viz N <slug>` · `/verify-chapter N` ·
`/notebook-chapter N` · `/publish-chapter N` · `/explain-concept N "<concept>"`.
Background skills loaded by agents: `fluids-book`, `teaching-style`, `interactive-viz`, `python-viz`, `math-to-python`,
`verify-implementation`, `colab-notebook`, `chapter-knowledge`, `data-and-benchmarks`.

## Machinery self-tests (run after changing any tool, the viz library or the embed code)
`.venv/Scripts/python.exe -m pytest -q tests/test_machinery.py` · `.venv/Scripts/python.exe tools/shot.py templates/viz_example.html` ·
`.venv/Scripts/python.exe tools/pipeline_selftest.py` (builds, executes, publishes and browser-audits a throw-away notebook).

## Style of code
Python ≥ 3.11, type hints, numpy docstrings, SI units, `pathlib`, functions over scripts, `if __name__ == "__main__":`
guards, figures to `outputs/chNN/`. Notebooks import from `fluidpy/`; they show short, heavily commented code that
*calls* the tested functions — and, for the key idea of a section, a transparent from-scratch version next to it so the
reader sees the mechanics, followed by an assertion that both agree.
