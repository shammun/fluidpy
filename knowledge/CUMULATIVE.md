# CUMULATIVE knowledge — fluidpy
(Rewritten by the knowledge-keeper after every chapter. Every agent reads this first.)

## Physics pipeline so far
(no chapter completed yet) — the book's structure: Ch1 vocabulary → Ch2 tensors → Ch3 kinematics → **Ch4 conservation
laws / Navier–Stokes (the trunk)** → Ch5 vorticity · Ch6 ideal flow · Ch7 gravity waves · Ch8 laminar flow → Ch9 boundary
layers → Ch10 CFD → Ch11 instability → Ch12 turbulence → **Ch13 GFD (uses 4, 5, 7, 8, 11, 12)** → Ch14 aerodynamics ·
Ch15 compressible · Ch16 biofluids.

## Available primitives
### Machinery (`fluidpy/core/`, ready before chapter 1)
| module | what | used for |
|---|---|---|
| `project` | repo root, `book.yaml` config, Pages/raw/Colab URLs | every tool and notebook |
| `style` | `setup_notebook()` (matplotlib style, plotly renderer, FAST), `COLORS`, `savefig` | every notebook figure |
| `embed` | `show_viz(chapter, slug)` — explainer in Jupyter (srcdoc) / Colab (Pages src) / page | explainer cells |
| `anim` | `animate`, `show_animation(player="video"/"frames")` | animations |
| `interact` | `slider_figure`, `animate_figure` (plotly, work on the page), `live` (ipywidgets) | Python interactives |
| `units` | pint registry, `dimensional_check`, `DIM` | V2 tests |
| `refdata` | reference-data resolution + registry | benchmarks |
### Physics primitives
| function | module | book § / Eq. | used by chapters | validation label |
|---|---|---|---|---|

### Explainer engine (`assets/viz_lib.js`)
`Viz.app` (tabs: Walkthrough / Explore / Equations / Check; fit-to-window with density levels and pagers), `Plot`
(world-coordinate plotting), `Viz.field` (grid, contours, heatmap, streamlines, quiver, particles), `Viz.num` (RK4,
odeint, brentq, erf, complex), KaTeX with fallback, selftest parity. Reference: `templates/viz_example.html`.

## Explainer inventory
| chapter | slug | CORE idea | reusable stage pattern |
|---|---|---|---|

## Notation
See `knowledge/notation.md` (symbol register) — conventions that change between chapters are flagged there.

## Teaching lessons
- **Depth is tiered; coverage is exhaustive** (adopted 2026-09-12, during ch01). The first ch01 curation made every new
  idea CORE at full depth: 76 items, 37 derivations. The design and implement agents then each ran ~50 min and ~500k
  tokens without finishing, and everything would have arrived at the same volume, so the ideas that carry the rest of
  the book did not stand out. Now every inventory row still gets an ID and a depth in an auditable chapter map, but
  only 12–18 load-bearing ideas get the full treatment (A: picture → question → derivation → worked number → code →
  figure). The rest are stated and explained inside the nearest A block (B), or named with a pointer to where they
  are used later (C). A derivation is written out only if it belongs to an A item or the book never writes it out.
  Policy: `book.yaml → policy`; rule: `.claude/agents/concept-curator.md`, `/do-chapter` §2, CLAUDE.md rule 6.

## Global pitfalls confirmed in this project
- Extracted text garbles maths (`¼` = `=`, `ð…Þ` = parentheses, missing minus, ω → `u`): read equations from rendered pages.
- Extracted text contains ligatures in the PDF; `tools/split_pdf.py` expands them (grep "fluid" works).
- Windows: shell heredocs mangle backslashes — write code files with Write/Edit; console needs UTF-8 (settings.json env).
- Anaconda's `python3` kernelspec can shadow the venv's — notebooks execute on `fluidpy-venv`.

## Benchmark inventory
| value / table | source (citation) | stored in | used by |
|---|---|---|---|

## Validation summary per chapter
| chapter | analytic | symbolic | converged | conserved | benchmark | book-value | qualitative | unverified | verdict |
|---|---|---|---|---|---|---|---|---|---|
