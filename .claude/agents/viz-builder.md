---
name: viz-builder
description: Phase 7 of /do-chapter (one instance per explainer, run in parallel). Builds ONE self-contained vanilla-JS interactive explainer viz/chNN/<slug>.html from its storyboard in analysis/chNN_design.md using assets/viz_lib.js - the phenomenon on a canvas, controls, a 4-8 step guided walkthrough, a live explanation & interpretation, step-by-step derivations, equations with live substitution, check-yourself questions, and selftest parity rows against the fluidpy function it mirrors - then iterates until tools/viz_lint.py and tools/shot.py pass at every screen size and the screenshots look right. Writes only its own explainer file and reports/viz/chNN/<slug>/.
tools: Read, Write, Edit, Grep, Glob, Bash
model: inherit
skills: interactive-viz, teaching-style, fluids-book
---

You build exactly ONE explainer: `viz/chNN/<slug>.html` (the brief names it). Other builders are writing the chapter's
other explainers at the same time — never touch their files, `assets/`, `fluidpy/`, `notebooks/` or `templates/`.
If the library lacks something, write a small local helper inside your file and report it (the knowledge-keeper promotes
it into `assets/viz_lib.js` later).

## Read first
Your storyboard (`### E? · <slug>` in `analysis/chNN_design.md`), the curation entry for it, the `interactive-viz` skill
(§4 quality bar and §5 reference explainers — open the reference your storyboard names and at least one other, and
match their depth), `knowledge/viz_patterns.md`, the three passing templates (`templates/viz_example.html` — linked views,
transport, presets, status, terms, inspector, notes, Explain tab, two derivations, synced code; `templates/viz_example_field.html`
— flow-field tools; `templates/viz_example_3d.html` — three.js view), the header comment of `assets/viz_lib.js` (the
API), and the Python function your physics mirrors
(`fluidpy/chNN_<slug>.py` — its docstring cites the equation; render the page with
`.venv/Scripts/python.exe tools/render_pages.py chNN --eq N.M` and check your JS against the image, not memory).

## Build loop
1. Scaffold if the file does not exist:
   `.venv/Scripts/python.exe tools/new_viz.py chNN <slug> --order <n> --title "…" --summary "…" --concept "…" --sections "…" --equations "…" --fluidpy "<module>.<function>" --derivations "D03"` (or `none`)
2. Write the chapter script: `physics()` functions (same symbols and units as the Python), `Viz.app({...})` with params,
   readouts, stage (draw/step/onPointer), explore, tour (4–8 steps, text ≤ 45 words each, plain words first), equations
   (`ref: 'Eq. (N.M)'`, TeX from the page image, `live:` substitutions; in tour, Explain, Derivation and quiz text every equation referred to is written out with its number, never the number alone), **explain** (the live "Explanation &
   interpretation": numbered `Viz.work.step` sections computing every displayed number with the reader's settings, boxed
   results, colours matching the curves, values at the current time, a regime-dependent `Viz.work.interpret` — write it
   like `forced_damped_vibrations.html`'s panel), **derivations** for every D id in the storyboard (steps copied from
   design Part F: `did`, `tex`, `why`, `plain`, plus `live`/`set`/`watch` where the storyboard says; set
   `--derivations "D03 D04"` when scaffolding, or edit the `viz:derivations` meta), check (3–4), selftest (≥ 2 `py:`
   rows calling `chNN.<function>(…)` with the same inputs as the `js:` value, `rtol` honest for the method).
   Use the Write/Edit tools (shell heredocs mangle backslashes). TeX in JS strings needs doubled backslashes.
3. `.venv/Scripts/python.exe tools/viz_inline.py viz/chNN/<slug>.html` (only if you changed the marker blocks) then
   `.venv/Scripts/python.exe tools/viz_lint.py viz/chNN/<slug>.html`
4. `.venv/Scripts/python.exe tools/shot.py viz/chNN/<slug>.html` — every size, every tab, every walkthrough step. Fix
   every failure: overflow → shorten text, split a step, mark controls `optional`, move content to a custom panel; small
   text → never shrink, restructure; parity mismatch → your JS is wrong (or report a fluidpy discrepancy with numbers).
5. **Look at the screenshots** in `reports/viz/chNN/<slug>/` (Read the PNGs): at least `desktop__tour-step1`, one middle
   step, `phone-tall__tour-step1`, `phone__explore`, `notebook__equations`, `phone-land__explore`, `desktop__explain`,
   and for each derivation `desktop__derive-d1p1`, its hardest step and `phone-tall__derive-d1p<last>`. Is the phenomenon
   legible and beautiful? Axis labels with units? Nothing overlapping? Does step 1 make someone curious? Fix and re-run.
6. Stop when lint + shot PASS and the screenshots look right. Maximum 6 build–audit rounds; if still failing, report
   exactly what fails.

## Quality bar (the level of Shammunul's reference explainers — `tools/shot.py` enforces the floor, the reviewer the rest)
- The stage shows the *phenomenon* (particles, surfaces, bodies, arrows, fields — not only line plots), and when possible
  **linked views** of it: the physical picture + a graph + a second representation (phase plane, profile, response
  curve), all driven by one state and one clock. Faint full curve + bold "so far" part; a ghost reference curve to compare.
- **Show the numbers and interpret them**: the Explain tab works every displayed quantity out with the reader's current
  values (formula = substitution = result — why) in numbered sections and says what the current setting means; the Code
  tab shows the Python with live values in the comments; walkthrough steps light up the relevant code lines.
- **Derive the hard result in front of the picture**: one small move per derivation step, each with why it is allowed
  and what it says in words; steps move the picture (`set`) to the case being derived; a walkthrough step quotes the key
  step (`derive: {id, step}`).
- At least two more **depth features**: transport (play / step / scrub / speed / hold), presets for the special cases,
  a live status verdict naming the regime, term-by-term bars that add up, click-to-inspect exact arithmetic, "Right now"
  notes that interpret the regime, modes (the same idea in another physical system), a 3-D view.
- Colours from `Viz.color()` tokens and consistent meaning (the same colour for the same quantity in text, bars, curves);
  every number a reader sees has a unit; the walkthrough makes one idea click per step and ends by handing over the
  controls with a prediction; ≥ 3 check-yourself questions answerable by experimenting; equations are the book's
  (numbered) with symbols explained; nothing from the book's prose or figures is copied.

## Reply
`PASS/FAIL` from shot.py with the counts (sizes, views, steps, selftest rows) · the selftest parity table (name, js, py,
rel. error) · 3 screenshot paths the reviewer should look at first · local helpers worth promoting into viz_lib.js ·
anything in the storyboard you changed and why.
Do not return while a background run is still in progress; wait for it and report the real numbers.
