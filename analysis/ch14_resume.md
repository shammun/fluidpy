# Chapter 14 — resume note (rewritten 2026-10-12 ~16:50, after the computer shut off at ~16:04 mid-review)

When the user writes **resume**: continue `/do-chapter 14` at the **review gates of phases 7 (viz) and 8 (notebook)**,
following `.claude/skills/do-chapter/SKILL.md` §6∥7∥8 onward. Delete this file once ch14 has `publish: pass`.

## Where things stand
| Phase | Status | Evidence |
|---|---|---|
| analyze … review | pass | commits up to `42472dc` |
| 7 viz | **all 9 built, lint + shot PASS; viz-reviewer verdict missing** | `reports/viz/ch14/*/audit.json` all `PASS` (15:54–16:04); `reports/ch14_viz.md` does not exist until the reviewer finishes |
| 8 notebook | **built and executed clean; lesson re-review verdict missing** | `outputs/ch14/executed.ipynb` 549 cells, 150 code cells executed, 0 errors (15:58); `notebooks/build_ch14.py` 3757 lines |
| 9 knowledge, 10 publish | todo | — |

Check first: does `reports/ch14_viz.md` exist, and does `reports/ch14_lesson.md` start with "re-review"? If not, the
reviewer that should have written it was interrupted again — launch it again (briefs below).

## Explainers (storyboards: `analysis/ch14_design.md` `### E<k> · <slug>`)
E1 `zhukhovsky_airfoil_builder` · E2 `kutta_condition` · E3 `starting_vortex` · E4 `airfoil_pressure_lift` ·
E5 `lift_curve_stall` · E6 `downwash_induced_drag` · E7 `lifting_line_solver` · E8 `elliptic_loading_aspect_ratio` ·
E9 `thrust_from_lift` — all built; backup B1 `wake_rollup` not built (only if an explainer fails review).

**The viz-reviewer brief must say**: re-run lint + shot itself; probe every tab for `.katex-error` (shot.py does not);
E3 "90 % after 6.1 chords" at dt = 0.05 c/U (6.2 is dt = 0.1) — computed or typed?; Wagner's function is a lift ratio,
never on a circulation axis; E2 uses internal `app._fromLoop` (fragile); E5 curve and stall are a labelled MODEL, no
book numbers; E6 staircase converges at order 1.5 (uniform) / 2.0 (cosine), 3-D view falls back to 2-D without WebGL;
fluidpy calls by keyword in Code tabs and `py:` rows. Must-fix → re-brief that `viz-builder` (max 2 rounds each;
builders must not edit sibling explainers or `assets/viz_lib.js`; `np` is not visible inside a lambda in `py:` rows).

## Notebook
`reports/ch14_lesson.md` round 1 (14:27, 541 cells) had 15 must-fix + 23 should-fix. The notebook-builder then changed
`notebooks/build_ch14.py` (14:51) and `scripts/ch14_drawings.py` (14:53); nobody confirmed which items were fixed
(two later review runs left only image folders `outputs/ch14/lesson_review2`, `…3`). **The lesson-reviewer brief must
say**: re-review; mark every round-1 item FIXED / PARTLY / NOT FIXED by content (cell numbers shifted); look at the
regenerated figures; Wagner = lift ratio starting at ½, 90 % after about 6.1 chords; keyword calls; overwrite the report.
Remaining must-fix → re-brief the `notebook-builder` with the report (this is fix round 2 of max 2), then review again.

## Then
Merge gate (SKILL.md): full `pytest -q` (≈ 40 min; not re-run since the review fixes unless noted in `progress.json`) ·
`tools/shot.py --chapter ch14 --quick` · `tools/run_notebook.py ch14 --save` · `tools/coverage_check.py ch14 --nb
outputs/ch14/executed.ipynb` · `tools/embed_check.py ch14` · both reports PASS → `viz: pass`, `notebook: pass`, commit →
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
- The reports of the E4–E9 builders were lost in the shutdown; anything they flagged for promotion has to be re-read
  from the explainer files themselves.

## Scratch files
`outputs/ch14/_resume_scratch/`, `outputs/ch14/lesson_review*/` — git-ignored, safe to delete.
