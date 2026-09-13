# fluidpy — Kundu, Cohen & Dowling, *Fluid Mechanics* (5th ed., 2012), learned through Python and interactive explainers

## Mission
Turn each chapter of the book into **one learning package** that makes the chapter *click* for a first-time reader
(and for Shammunul, who is working through the book for a PhD in climate dynamics):

| Per chapter | What it is |
|---|---|
| `notebooks/chNN_<slug>.ipynb` | The teaching notebook, **executed**: plain words → step-by-step maths (hard results as full **derivations**, one move per step, each explained) → tiny worked example → commented Python → figure → how to read it; Python animations, plotly slider figures, live widgets; the chapter's explainers embedded |
| `notebooks/chNN_<slug>_colab.ipynb` | The same notebook for Google Colab (outputs stripped; its setup cell clones this repo; explainers load from GitHub Pages) |
| `notebooks/chNN_<slug>.html` | The published page: executed notebook, every explainer a **full-window block** (100 % width × 100 % height), Open-in-Colab button, prev/next, contents |
| `viz/chNN/<slug>.html` (**4–5**) | Self-contained vanilla-JS explainers at the depth of Shammunul's reference explainers (preferred: the Unit 4 mathlets — forced damped vibrations, amplitude & phase, angular frequency explorer): linked views of the phenomenon, controls and presets, guided walkthrough, a live **Explain** tab (every number worked out with your settings + interpretation), a step-by-step **Derivation** tab for derivation-heavy ideas, synced code, live equations, quiz — fitting any window with **no scrolling** |

Plus the machinery that makes it trustworthy and cumulative: `fluidpy/` (the physics as tested functions), `tests/`,
`reports/`, `knowledge/`. Site: `index.html` (chapters) and `viz/index.html` (explainer gallery), served at
`https://shammun.github.io/fluidpy/` (all URLs derive from `book.yaml → project`).

**Coverage exhaustive, depth tiered; nothing used unexplained.** Every book section appears in the notebook. The
`concept-curator` gives every inventory item (definition, theorem, numbered equation, example) an ID and a depth in an
auditable chapter map (`book.yaml → policy`): **A · full treatment** = **CORE** (`C01…`), the 12–18 load-bearing ideas
(hard max 18) — picture → question → step-by-step derivation → worked number → commented Python → **at least one
visualization**; **B · stated and explained** = **NOTE** (`N01…`) inside the nearest A block — a paragraph, the equation,
a number; **C · named** = **NOTE** with a sentence and a pointer to where it is used later; **RECAP** (`R01…`) = taught in
an earlier chapter, reminded where needed; **SKIP** (`S01…`) = exercises, bibliography, deferred material — one pointer
line. A result obtained by manipulating equations is written out as a **DERIVATION** (`D01…`, ★–★★★) only if it belongs
to an A item or the book never writes it out (other results are stated, not derived): one small move per step (what we
did · the line · why it is allowed · in words), the book's skipped moves filled in, then checked (units, limits, sympy
for ★★★) and interpreted — in the notebook (hence Colab and the page) and in the Derivation tab of the explainer for that
idea. Every concept, symbol, maths tool and
Python function is explained where it is first used (CORE block, recap, or 📎 primer), tracked in a prerequisite ledger. 4–5 interactive explainers
are chosen from the CORE ideas where interaction teaches most. `tools/nbkit.py`, `tools/coverage_check.py`, the
`lesson-reviewer` and `tools/shot.py` enforce all of this.

## One command
`/do-chapter N` runs the whole pipeline for chapter N (10 phases, parallel where safe, resumable with `--from PHASE`).
`/status` says where things are and what to run next. `GUIDE.md` explains everything else.

| # | Phase | Agent(s) | Writes |
|---|---|---|---|
| 1 | analyze | `concept-analyst` | `analysis/chNN.md` |
| 2 | curate | `concept-curator` | `analysis/chNN_curation.md` |
| 3 | design ∥ 4 implement | `lesson-designer` ∥ `concept-implementer` | `analysis/chNN_design.md` ∥ `fluidpy/`, `scripts/` |
| 5 | verify | `math-verifier` (⟲ implementer, ≤ 3 loops) | `tests/`, `reference/`, `reports/chNN_verification.md` |
| 6 | review ∥ 7 viz ∥ 8 notebook | `derivation-reviewer` ∥ `viz-builder` ×4–5 → `viz-reviewer` ∥ `notebook-builder` → `lesson-reviewer` | `reports/chNN_review.md` ∥ `viz/chNN/`, `reports/chNN_viz.md` ∥ `notebooks/`, `reports/chNN_lesson.md` |
| 9 | knowledge ∥ 10 publish | `knowledge-keeper` ∥ `site-publisher` | `knowledge/` ∥ pages, Colab twin, index, gallery |

## Where things live
| Path | What | Committed? |
|---|---|---|
| `*.pdf` (repo root) | the book — read-only, never modified | **never** |
| `book.yaml` | project config + chapter map (exact PDF pages from the outline, sections, explainer seeds) | yes |
| `chapters/chNN.{pdf,txt}`, `chapters/pages/` | split chapter, extracted text, rendered page images | **never** |
| `analysis/chNN.md`, `chNN_curation.md`, `chNN_design.md` | inventory (NEW/SEEN, prerequisites, derivations and the moves the book skips) · IDs + tiers + derivations (D rows) + section coverage + 4–5 explainers · storyboards + prerequisite ledger + derivations written out step by step | yes |
| `fluidpy/core/` | primitives: `project`, `embed` (show_viz), `anim`, `interact` (plotly sliders), `style`, `units`, `refdata`, + physics reused by ≥2 chapters | yes |
| `fluidpy/chNN_<slug>.py`, `scripts/chNN_*.py` | chapter physics + runnable demos | yes |
| `viz/chNN/<slug>.html` | explainers (library inlined from `assets/viz_lib.js` + `assets/viz_base.css` by `tools/viz_inline.py`) | yes |
| `tests/test_chNN.py`, `reference/chNN/` | evidence (+ cited public benchmark data) | yes |
| `reports/chNN_{verification,review,viz,lesson}.md` | verdicts | yes |
| `reports/viz/**` | screenshots + audit JSON from `tools/shot.py` | no |
| `notebooks/build_chNN.py` → `chNN_<slug>.ipynb` (+ `_colab.ipynb`, `.html`) | builder (uses `tools/nbkit.py`) and its products | yes |
| `knowledge/` | `CUMULATIVE.md`, `concept_map.md`, `notation.md`, `primers.md`, `viz_patterns.md`, `chNN.md` | yes |
| `progress.json` | phase status per chapter | yes |
| `outputs/` | figures, GIFs, self-test pages | no |

## Non-negotiable rules
1. **Sequential chapters.** Start chapter N only when chapter N−1 has `publish: pass` and `knowledge: done` (or the
   user explicitly overrides). Every agent reads `knowledge/CUMULATIVE.md`, `knowledge/notation.md` and
   `knowledge/viz_patterns.md` before chapter work.
2. **Equations come from page images.** `chapters/chNN.txt` garbles maths (`¼` is `=`, `ð…Þ` are parentheses, a minus
   vanishes, ω prints as `u`). Render the defining page (`.venv/Scripts/python.exe tools/render_pages.py chNN --eq 7.27`)
   and read the PNG before transcribing, coding or displaying an equation. Never code an equation from memory.
3. **Teach in Shammunul's style** (skill `teaching-style`): every new idea is a CORE block with plain words first,
   *problem → idea → maths → code → what the output shows*, a tiny example with easy numbers traced step by step, every
   code line commented for a novice, **at least one visualization**, and *what you see / how to read it / what would
   change if…*; nothing (concept, symbol, maths tool, Python function) is used before it is explained — CORE block,
   recap or 📎 primer. **Derivations are never skipped or compressed**: goal and plan in plain words, one small move
   per step with *what we did*, *why we can do this* (the rule, the assumption) and *in words*, then a check and what
   the result means (`teaching-style` §1c, `nb.derivation`).
4. **4–5 interactive explainers per chapter** (plus a backup idea), each attached to CORE ideas and justified by "why
   interaction beats a static figure here", built to the depth of the reference explainers (skill `interactive-viz`
   §4–§5): a required live **Explain** tab ("Explanation & interpretation" in numbered sections, as in
   `forced_damped_vibrations.html`), a **Derivation** tab for every derivation the curation assigns to it, a synced
   Code tab, plus at least two depth features (linked views, transport, presets, status, term bars, inspector, notes,
   modes, 3-D). Each fits the window with **no scrolling** at 360×640 …
   1920×1080 and in a notebook frame, shows and explains its equations, has a 4–8-step walkthrough and ≥ 3 check
   questions. Content that does not fit becomes a tab or a page — never a scrollbar, never text below 12 px. Enforced by
   `tools/viz_lint.py` + `tools/shot.py` and judged by eye by the `viz-reviewer`.
5. **Traceable physics.** Every `fluidpy` function's docstring cites book section + equation number(s), lists symbols
   with SI units and assumptions, and carries its validation label. Explainer JS physics mirrors a `fluidpy` function
   and proves it with `selftest()` parity rows (`py:` expressions evaluated by `tools/shot.py`).
6. **Evidence, proportional to importance** (skill `verify-implementation`): computable CORE items ≥ 2 independent
   evidence levels (one of V1 analytic / V2 symbolic / V3 convergence / V5 benchmark); coded NOTE items and every
   function the notebook or explainers call ≥ 1. Never loosen a tolerance to pass;
   fix the physics or report the discrepancy. Max 3 fix loops, then stop and tell the user.
   **Depth is tiered, coverage is exhaustive**: every inventory item gets an ID and a depth (A full treatment, 12–18 per
   chapter, hard max 18 · B stated inside an A block · C named with a pointer); derivations are written out only for A
   items or results the book never writes out (`book.yaml → policy`).
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
