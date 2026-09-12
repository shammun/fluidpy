---
name: colab-notebook
description: How the per-chapter notebook, its Colab twin and its web page are built, executed and published in fluidpy - notebooks/build_chNN.py with tools/nbkit.py, the setup cell that works locally and on Colab, cell tags, FAST runs, headless execution on the fluidpy-venv kernel (tools/run_notebook.py), what tools/publish_notebook.py produces, public-repo safety and the Colab save hazard. Load when creating, fixing or publishing any .ipynb.
---

# colab-notebook — one notebook, three deliverables

## 1. Build with code, never by hand
`notebooks/build_chNN.py` (committed) uses `tools/nbkit.py`:
```python
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "tools"))
from nbkit import ChapterNotebook

nb = ChapterNotebook("ch07")
nb.title(big_idea="…", roadmap=["…", "…"], prerequisites=["Bernoulli (Ch. 4)", "potential flow (Ch. 6)"])
nb.explainer_index([("dispersion_relation", "Why long waves outrun short ones", "phase vs group speed"), …])   # 4–5 rows
nb.setup()                                          # ⚙️ setup markdown + the setup code cell (tag: setup)
nb.section("7.2", "Linear surface waves", intro="**What is this section about?** …")   # one per book section
nb.recap("R02", "Bernoulli's equation", "…one paragraph…", where="Ch. 4 §4.9")        # RECAP tier
nb.core("C03", "The dispersion relation", question="Why do long waves travel faster than short ones?")   # CORE block
nb.md("#### The problem in plain words\n…")
nb.primer("tanh(x)", "A smooth step: ≈ x for small x, → 1 for large x. Here it switches between shallow and deep water.",
          code="import numpy as np                 # numbers\nprint(np.tanh([0.1, 1, 5]))          # [0.0997 0.7616 0.9999]")
nb.worked_example("a 10 m wave in 2 m of water", "1. k = 2π/10 ≈ 0.63 m⁻¹ …")
nb.code("""
from fluidpy import ch07_gravity_waves as ch07      # the tested chapter module
omega = ch07.omega(k, H)                            # Eq. (7.36): frequency for wavenumber k in depth H
""", explain="1. … 2. …")
nb.figure("…plotting code, every line commented…", see="…", read="…", change="…")   # at least one visual per CORE block
nb.explainer("dispersion_relation", heading="Why long waves outrun short ones", why="…", tries=["Drag H…", "…"])
nb.note("The deep-water limit is the same relation with tanh → 1.", equation=r"\omega^2 = gk", ref="7.xx")   # NOTE
nb.animation("…"); nb.plotly("…"); nb.live("…"); nb.check_agree("…")
nb.pointer("Nonlinear shallow-water waves are only mentioned here; Ch. 13 returns to shallow water.")  # SKIP
nb.summary(clicked=["one line per CORE idea"], feeds_forward=["…"], left_out=["… (where to find it)"])
nb.save()                                           # → notebooks/ch07_gravity_waves.ipynb (no outputs)
```
`save()` refuses to write the notebook when a book section has no `nb.section(...)`, a CORE/RECAP id from the curation
has no block, a CORE block has no code or no visual, or the explainer count is outside 4–5. After executing,
`tools/coverage_check.py chNN --nb outputs/chNN/executed.ipynb` repeats the checks on real outputs and checks the
prerequisite ledger (design Part E).
Markdown is our own words; equations in LaTeX with their book numbers. Physics lives in `fluidpy/`, never only in a cell.

## 2. The setup cell (written by `nb.setup()`; identical in every chapter)
Colab: clone `https://github.com/<repo>.git` into `/content/fluidpy` (or `git pull`), `pip install -r
requirements-colab.txt` (only what Colab lacks — currently `pint`; never re-pin numpy/scipy, which would force a runtime
restart). Local: walk up to the folder with `book.yaml`. Then `os.chdir(ROOT)`, `sys.path.insert(0, ROOT)`, import
`setup_notebook`, `show_viz`, `show_animation`, and set `FAST = setup_notebook()` (env `FLUIDPY_FAST=1`).

## 3. Cell tags (honoured by the publish tool)
`setup` · `explainer` (show_viz output → full-window block on the page, link card in the GitHub .ipynb) ·
`live-only` (ipywidgets → "run it in Colab" note on the page) · `animation` · `plotly` · `from-scratch` · `slow`.

## 4. Execute headlessly = verification
`.venv/Scripts/python.exe tools/run_notebook.py chNN [--fast] [--save]` — runs on the `fluidpy-venv` kernel (registered
automatically; it guarantees the venv interpreter, not Anaconda's `python3` kernel), prints runtime, slowest cells,
errors with the last traceback line, and outputs larger than 3 MB. Never `jupyter nbconvert` (PATH may dispatch to
Anaconda). Budget: < 5 min full run; `FAST` for anything heavier.

## 5. Publish (`tools/publish_notebook.py chNN`, run by the site-publisher)
1. Gate: `tools/embed_check.py chNN` (explainers exist, are current, lint-clean, embedded exactly once).
2. Execute (0 errors; no `book values active` banner in any output).
3. Write `notebooks/chNN_<slug>_colab.ipynb` — outputs stripped, Colab metadata.
4. Write `notebooks/chNN_<slug>.ipynb` — **executed** (figures visible on GitHub), explainer outputs as link cards.
5. Write `notebooks/chNN_<slug>.html` — nbconvert lab template + site chrome (All chapters · All explainers ·
   prev/next · Open in Colab · Download .ipynb · View on GitHub · explainer cards), explainers full-window
   (`section.fp-viz`, 100 % × 100 % of the window), live-only notes, floating Contents (hidden while an explainer fills
   the screen), MathJax for markdown maths, plotly via require.js + CDN.
6. Regenerate `index.html`, `viz/index.html`, the README table (`tools/build_site.py`).
7. Audit: `tools/shot.py --page notebooks/chNN_<slug>.html`; then `tools/check_public.py`.

## 6. Public-repo safety and the Colab save hazard
- The executed `.ipynb` and the page come from the public-safe run (no private book values). `tools/check_public.py`
  fails on: any `_colab.ipynb` with outputs, the private-run banner anywhere, tracked `chapters/`, `outputs/`,
  `reports/viz/`, `tests/book_values_*`, page images.
- **Never use Colab's File → Save a copy in GitHub**: it commits outputs into `_colab.ipynb` on `main`. Save to Drive.
  If it happens: fetch, rebase, re-run `tools/publish_notebook.py chNN`, commit
  `chNN: restore outputs-stripped _colab.ipynb`, tell the user the old commit stays reachable by SHA; never force-push.

## 7. Lessons (the knowledge-keeper appends; one line each, dated by chapter)
- (machinery) `import site` inside tools resolves to the standard library module — the site builder is `tools/build_site.py`.
- (machinery) Anaconda's `python3` kernelspec can shadow the venv's; execution uses the dedicated `fluidpy-venv` kernel.
