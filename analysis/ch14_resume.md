# Chapter 14 — resume note (written 2026-10-12, session stopped by the user mid-phase)

When the user writes **resume**: continue `/do-chapter 14` from phases **7 (viz) ∥ 8 (notebook)**, following
`.claude/skills/do-chapter/SKILL.md` §6∥7∥8 onward. Delete this file once ch14 has `publish: pass`.

## Where things stand
| Phase | Status | Evidence |
|---|---|---|
| analyze, curate, design, implement, verify | pass | commits up to `42d930d` |
| 6 review | **pass**, committed `42472dc` | `reports/ch14_review.md` (1 must-fix + 15 should-fix, all resolved); `tests/test_ch14.py` 95 passed |
| 7 viz | **running — 3 of 9 built** | see table below |
| 8 notebook | **running — builder stopped mid-work** | `notebooks/build_ch14.py` (3470 lines), unfinished |
| 9 knowledge, 10 publish | todo | — |

The full test suite has NOT been re-run since the review fixes (only `tests/test_ch14.py`); it runs at the merge gate.

## Explainers (storyboards: `analysis/ch14_design.md` `### E<k> · <slug>`; summaries: `analysis/ch14_curation.md` §5)
| # | slug | State on disk | Next action |
|---|---|---|---|
| E1 | `zhukhovsky_airfoil_builder` | **built**: lint ok, shot PASS 8 sizes, 35/35 selftest rows | viz-reviewer |
| E2 | `kutta_condition` | **built**: lint ok, shot PASS 8 sizes, 38/38 rows | viz-reviewer (uses internal `app._fromLoop` for the 0→1 ramp — fragile spot) |
| E3 | `starting_vortex` | **built**: lint ok, shot PASS 8 sizes, 20/20 rows | viz-reviewer; its text says "90 % after 6.2 chords" — the verified value is 6.1 at dt = 0.05 c/U (6.2 is dt = 0.1): check whether the page computes it or has it typed, fix if typed |
| E4 | `airfoil_pressure_lift` | **partial** (3104 lines, builder stopped mid-build; it was about to apply small fixes — ghost NaN, quiz Q1, a `float` builtin in a py: row — from a scratch script and splice them in) | re-brief a `viz-builder`: "file is an unfinished build, continue it; run lint + shot until PASS" |
| E5 | `lift_curve_stall` | **scaffold only** (new_viz.py template, TODO marker) | brief a `viz-builder` from the start |
| E6 | `downwash_induced_drag` | **scaffold only** | brief a `viz-builder` from the start |
| E7 | `lifting_line_solver` | not started | brief a `viz-builder` (new_viz.py scaffold first) |
| E8 | `elliptic_loading_aspect_ratio` | not started | same |
| E9 | `thrust_from_lift` | not started | same |
| B1 | `wake_rollup` (backup) | not built | only if an explainer fails review |

At most 3 viz-builders at a time (`book.yaml → project.max_parallel_viz`). After all nine: `viz-reviewer` →
`reports/ch14_viz.md`.

**Things every remaining viz brief must say** (learned this session):
- call fluidpy functions **by keyword** in Code tabs and `py:` parity rows (sibling functions have different positional
  orders — review item 11);
- `tools/shot.py` does not detect KaTeX errors: the builder checks for `.katex-error` on every tab itself (an
  under-escaped `\beta` in a JS string becomes a backspace);
- `shot.py`'s evaluator does not expose `np` inside a lambda body in `py:` expressions (bind it as a default argument);
- finished siblings to copy local helpers from (`cx` complex arithmetic, `measure(contour)`, `arrowPx`, `spinArrow`,
  phone layouts): E1, E2, E3; do not edit them or `assets/viz_lib.js`;
- E6: staircase of horseshoes converges at order 1.5 (uniform panels) / 2.0 (cosine) — never "first order"; 3-D view
  must fall back to the oblique 2-D drawing without WebGL;
- E5: the measured-looking curve and stall are a drawn MODEL, labelled so; no book numbers.

## Notebook (phase 8)
`notebooks/build_ch14.py` and `notebooks/ch14_aerodynamics.ipynb` exist but the builder was stopped before it reported:
treat them as **inputs, not results**. Re-brief the `notebook-builder`: "an unfinished builder exists — read it, finish
it against `analysis/ch14_design.md` Parts A/C/D/E/F, then loop build → `tools/run_notebook.py ch14 --save` →
`tools/coverage_check.py ch14 --nb outputs/ch14/executed.ipynb` until clean". The brief must repeat the two corrections
the builder was sent (it may or may not have applied them — have it check):
1. **Wagner's function is a lift ratio** (circulatory lift / steady lift, starts at ½); the exact bound circulation
   starts at 0. In C03 it is plotted only against `lift_circulatory/(ρUΓ_steady)` on a lift panel, never on the
   circulation axis. `scripts/ch14_starting_vortex.py`'s figure now has 4 panels (was 3). Say "the model reaches 90 % of
   the steady value after about 6.1 chords" (dt = 0.05 c/U).
2. **Keyword calls** for all fluidpy functions.
Explainer slugs to embed once each: the nine above (CORE attachment in curation §5).
Then `lesson-reviewer` → `reports/ch14_lesson.md` (max 2 fix rounds).

## Then
Merge gate (SKILL.md): full `pytest -q` · `tools/shot.py --chapter ch14 --quick` · `tools/run_notebook.py ch14 --save` ·
`tools/coverage_check.py` · `tools/embed_check.py ch14` · both reports PASS → `viz: pass`, `notebook: pass`, commit →
phases 9 ∥ 10 (knowledge-keeper ∥ site-publisher) → publish gate (push needs the user's go-ahead as usual).

## Carry into the knowledge phase (promotion candidates and tool gaps reported by builders)
- `assets/viz_lib.js` candidates: `cx` (complex arithmetic incl. principal sqrt; E1, E2, E4 each carry a copy),
  `measure(contour)` / `parabolicPeak` / `interpSorted`, `buildNet` (multi-level marching squares clipped at a circle),
  `spinArrow`, `arrowPx`, `arcPx`, `gutterPlot`, `planePlot`, `layoutOf()`, signed formatters `sg` / `z3`, a log-axis
  option for `Plot`, a `transport.max` / `holdAt` option (would remove E2's use of `app._fromLoop`).
- `tools/shot.py` gaps: no `.katex-error` probe; no click/drag pass; `np` not visible inside lambdas in `py:` rows.
- Lesson for the skills: a cited "growth" curve must be checked for *which quantity* it is a ratio of before it is
  overlaid (the Wagner mislabel passed verify because lift ratio = circulation ratio numerically in the model).
- The exact-circulation numbers quoted in docstrings (0 at start, 0.55 after one chord, 90 % at 7.1 chords) are the
  reviewer's own computation — attributed, not independently tested.

## Scratch files
The stopped agents' scratch work (E4 splice scripts, notebook part files) was copied to the git-ignored
`outputs/ch14/_resume_scratch/` — optional input for the re-briefed builders, safe to delete.
