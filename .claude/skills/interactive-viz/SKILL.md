---
name: interactive-viz
description: How fluidpy interactive explainers are designed, built, fitted to the window, verified and embedded - the explainer contract, the assets/viz_lib.js API (Viz.app, Plot, field tools, particles, pager, walkthrough, equations, selftest), no-scroll layout rules, design patterns for fluid phenomena, the lint/shot audit loop, and how show_viz and the publish tool put explainers full-window into notebooks, Colab and the web page. Load when choosing, storyboarding, building, reviewing or embedding an explainer.
---

# interactive-viz — explainers that fit the window and make an idea click

## 1. The contract (checked by `tools/viz_lint.py`, `tools/shot.py`, `tools/embed_check.py`)
- One self-contained file `viz/chNN/<slug>.html`, created with `tools/new_viz.py`; `assets/viz_base.css` and
  `assets/viz_lib.js` are inlined between marker comments by `tools/viz_inline.py` (edit the assets, never the copies).
  Only external resource: KaTeX, loaded by the library from jsDelivr/cdnjs with a text fallback.
- `<meta name="viz:*">`: chapter, slug, order, title, summary, concept, sections, equations, fluidpy.
- ≤ 5 per chapter, each tied to a CORE idea from `analysis/chNN_curation.md`.
- **Fits the window with no scrolling** at 360×640, 390×844, 844×390, 768×1024, 1280×720, 1000×700 (notebook),
  1366×768, 1920×1080 — the app fills the window exactly, nothing overflows, no text below 12 px, tap targets ≥ 24 px.
- Tabs: **Walkthrough** (4–8 steps) · **Explore** (controls + live numbers + callouts) · **Equations** (book-numbered,
  live substitutions, symbols) · **Check yourself** (3–4 questions) · optional custom panels.
- `selftest()` rows: `{name, js, py: "chNN.fn(…)", rtol}` (JS ↔ fluidpy parity, evaluated by shot.py) and
  `{name, js, expect, rtol}` (JS-only invariants). At least 2 `py` rows.

## 2. `Viz.app(config)` — the API in one screen (full doc: header of `assets/viz_lib.js`; example: `templates/viz_example.html`)
```js
const app = Viz.app({
  title, subtitle,
  params: { k: { label: 'Wavenumber $k$', min: 0.05, max: 3, step: 0.01, value: 1, unit: '1/m', help: '…', log: false, optional: false },
            deep: { type: 'toggle', label: 'Deep water', value: false },
            mode: { type: 'select', options: [['a', 'A'], ['b', 'B']], value: 'a' },
            kick: { type: 'button', label: 'Kick', action: app => … },
            speed: { label: 'Speed', min: 0.1, max: 3, value: 1, optional: true } },   // 'speed' scales animation time
  readouts: [{ id: 'c', label: 'Phase speed $c$', unit: 'm/s', value: s => …, tone: v => v > 0 ? 'pos' : 'neg' }],
  stage: { animate: true, setup(g, s) {}, step(g, s, dt) {}, draw(g, s) {}, reset(g, s) {}, onPointer(g, ev, s) {} },
  explore: { intro: 'html with $tex$', controls: ['k', 'deep'], readouts: ['c'], callouts: [{ kind: 'try', html: '…' }] },
  tour: [{ title, text, set: { k: 1 }, controls: ['k'], readouts: ['c'], eq: 'disp', highlight: ['readout:c'],
           callout: { kind: 'key', html }, play: true, enter(app) {} }],
  equations: [{ id: 'disp', title: 'Dispersion relation', ref: 'Eq. (7.36)', tex: '\\omega^2 = g k \\tanh(kH)',
                live: s => `\\omega = ${Viz.tnum(omega(s))}\\ \\text{s}^{-1}`, note: '…', symbols: [['\\omega', 'angular frequency', 'rad/s']] }],
  check: [{ q: '…', a: '…', set: { … } }],
  panels: [{ id: 'orbits', label: 'Orbits', short: 'Orbits', render(el, app) {} }],
  selftest: () => [{ name: 'omega(k=1,H=10)', js: omega({ k: 1, H: 10 }), py: 'ch07.omega(1.0, 10.0)', rtol: 1e-10 }],
  onChange(s, key, app) {}
});
```
- Inside callbacks use `g.app` (not the outer `app` constant — it does not exist yet during the first draw).
- Graphics: `const P = g.plot({ xlim, ylim, equal, xlabel: 'x [m]', ylabel, pad, rect })` then `P.axes()`,
  `P.line(xs, ys, {color, width, dash})`, `P.fn(f)`, `P.fillBetween`, `P.dot`, `P.arrow(x, y, dx, dy)`, `P.text`,
  `P.clip(fn)`, `P.X/P.Y/P.ix/P.iy` (world ↔ pixel). Several plots in one stage: pass `rect: {x, y, w, h}`.
- Fields: `Viz.field.grid(f, xlim, ylim, nx, ny)`, `.contour(G, level)`, `.drawContours(P, G, levels)`,
  `.heatmap(P, f, {cmap:'viridis'|'coolwarm'|'blues'|'ocean', vmin, vmax})`, `.streamline(vel, x0, y0, {ds, n, bounds})`,
  `.quiver(P, vel, {nx})`, `.particles(vel, {n, bounds, life, trail, seed, spawn, kill})` → `.step(dt, t)`, `.draw(P)`.
- Numerics: `Viz.num.linspace, trapz, rk4Step, odeint, brentq, erf, erfc, niceTicks, C` (complex ops). `Viz.rng(seed)`.
- Formatting: `Viz.fmt(v, {sig, unit})` for labels, `Viz.tnum(v)` inside TeX. Colours: `Viz.color('accent'|'teal'|
  'orange'|'rose'|'blue'|'amber'|'muted'|'text')`, `Viz.alpha(c, a)` — never hard-code colours (dark theme).
- Inline maths in any text: `$…$`; display: `$$…$$`. In JS strings double the backslashes: `'$\\omega$'`.

## 3. Fit rules (how to pass the audit without shrinking anything)
- The library already: picks layout (`wide` stage|side; `portrait` stage over side; `landscape` compact), raises
  density 0–3 (tighter spacing; level 2 hides help lines; level 3 hides `optional` items and the explore intro), pages
  long lists (controls, equations, questions) with ‹ 1/2 ›, shortens tab labels.
- You: keep step text ≤ 45 words; ≤ 2 extras per step (controls / readouts / eq / callout); ≤ 5 controls in Explore
  (mark the rest `optional: true`); readout labels ≤ 22 chars; equations that are wide use `\begin{aligned}` over two
  lines; put a second visual in a custom panel (tab) instead of squeezing two plots on a phone.
- Never add `overflow: auto`, fixed widths > 360 px, or font sizes < 12 px in chapter CSS.
- The stage must look right in any aspect ratio: compute layout from `g.w`/`g.h` in `draw`; use `equal: true` for
  geometry; for phones (`g.w < 420`) reduce annotation density rather than font size.

## 4. Patterns that suit fluid mechanics
| Idea type | Stage pattern |
|---|---|
| a field (velocity, pressure, vorticity) | heatmap + contours/streamlines + tracer particles; a draggable probe showing the local vector |
| superposition / linear theory | toggles for each component + the sum drawn thicker; readout of the combined quantity |
| waves / dispersion | animated surface + particles in orbits; a second rect with the ω(k) curve and a moving dot |
| boundary / similarity layers | profiles at several times or stations that collapse onto one curve when rescaled (toggle "rescale") |
| stability | growth-rate curve vs wavenumber + an animated perturbation that grows/decays for the selected k |
| balances (geostrophy, Ekman, drag) | force arrows at a parcel that update with sliders; "turn off term X" toggles |
| control volumes / budgets | a box with inflow/outflow arrows proportional to fluxes and a live budget table |
| numerics (Ch. 10) | the grid and stencil drawn; step-by-step time marching with play/step; error vs h on log axes |
Each explainer should have at least one *motion* or *direct manipulation* that a static figure cannot give.

## 5. Build → audit loop
`tools/new_viz.py` → write → `tools/viz_lint.py <file>` → `tools/shot.py <file>` (≈ 1 min; screenshots in
`reports/viz/chNN/<slug>/`) → **Read the PNGs** (desktop step 1, phone step 1, phone explore, notebook equations,
phone-land explore) → fix → repeat (≤ 6 rounds). `--quick` audits only 3 sizes while iterating; the final run is full.

## 6. Embedding (done by the notebook-builder and the publish tool — explainer authors need not do anything)
- Notebook cell: `show_viz("chNN", "<slug>")` (via `nb.explainer(...)`). Local Jupyter/VS Code → `<iframe srcdoc>` of the
  file on disk; Colab → `<iframe src>` of the GitHub Pages copy (srcdoc from the clone if Pages is not live yet);
  the iframe is output-width × visible-window-height with Full screen / Open in new tab.
- Published page: `tools/publish_notebook.py` swaps each output for a `section.fp-viz` block that its page script moves
  out of the cell and sizes to exactly 100 % of the viewport width and height (`../viz/chNN/<slug>.html`, lazy).
  `tools/shot.py --page` verifies it.
- The GitHub `.ipynb` view shows a link card instead (GitHub strips iframes).

## 7. Lessons (the knowledge-keeper appends; one line each, dated by chapter)
- (machinery) The `#app` mount must be 100 % high or the no-overflow audit cannot see clipping — fixed in viz_base.css;
  the audit also fails `app-does-not-fill-window`.
- (machinery) Explore callouts are optional by default; on phones they are the first thing to go.
