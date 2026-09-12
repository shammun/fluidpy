---
name: interactive-viz
description: How fluidpy interactive explainers are designed, built, fitted to the window, verified and embedded - the explainer contract (4–5 per chapter), the assets/viz_lib.js API (Viz.app with linked views, transport, presets, status, terms, inspector, notes, step-by-step working, synced code, modes, 3-D via three.js; Plot, field tools, particles, pager), the quality bar distilled from Shammunul's reference explainers (with their paths), no-scroll layout rules, patterns for fluid phenomena, the lint/shot audit loop, and how show_viz and the publish tool put explainers full-window into notebooks, Colab and the web page. Load when choosing, storyboarding, building, reviewing or embedding an explainer.
---

# interactive-viz — explainers that fit the window and make an idea click

## 1. The contract (checked by `tools/viz_lint.py`, `tools/shot.py`, `tools/embed_check.py`)
- **4–5 per chapter** (`book.yaml → project.min/max_explainers_per_chapter`), each tied to CORE ideas of
  `analysis/chNN_curation.md`, plus one backup idea in the curation in case one fails review.
- One self-contained file `viz/chNN/<slug>.html`, created with `tools/new_viz.py`; `assets/viz_base.css` and
  `assets/viz_lib.js` are inlined between marker comments by `tools/viz_inline.py` (edit the assets, never the copies).
  External resources: only KaTeX and (for 3-D) three.js, both loaded by the library from CDNs with fallbacks.
- `<meta name="viz:*">`: chapter, slug, order, title, summary, concept, sections, equations, fluidpy.
- **Fits the window with no scrolling** at 360×640, 390×844, 844×390, 768×1024, 1280×720, 1000×700 (notebook),
  1366×768, 1920×1080 — fills the window exactly, nothing overflows, no text below 12 px, tap targets ≥ 24 px.
- Tabs: **Walkthrough** (4–8 steps) · **Explore** · **Step by step** (required) · **Equations** (book-numbered, live) ·
  **Code** (required) · **Check yourself** (≥ 3) · optional custom panels.
- **Depth features: at least 2 of** linked views (≥ 2), transport, presets, status, terms, inspector, notes, modes.
- `selftest()` rows: ≥ 2 `{name, js, py: "chNN.fn(…)", rtol}` (JS ↔ fluidpy parity, evaluated by shot.py) plus optional
  `{name, js, expect, rtol}` invariants.

## 2. `Viz.app(config)` — the API (full reference: header of `assets/viz_lib.js`; complete example: `templates/viz_example.html`)
```js
const app = Viz.app({
  title, subtitle,
  params: { k: { label: 'Wavenumber $k$', min: 0.05, max: 3, step: 0.01, value: 1, unit: '1/m', help: '…', optional: false },
            deep: { type: 'toggle', label: 'Deep water', value: false },
            sys: { type: 'select', options: [['a', 'A'], ['b', 'B']], value: 'a', display: 'chips' } },
  stage: {
    rows: [1.1, 1],                                        // relative row heights
    views: [ { id: 'flow', row: 0, flex: 1.4, title: s => `<b>The flow</b> · c = ${Viz.fmt(c(s))} m/s` },
             { id: 'curve', row: 0, flex: 1, hidePortrait: true, title: '<b>ω(k)</b>' },
             { id: 'graph', row: 1, title: '<b>η(x, t)</b>' } ],
    setup(g, s) {}, step(g, s, dt) {}, reset(g, s) {},
    draw(g, s) { const v = g.view('graph'); v.clear(); const P = v.plot({ xlim, ylim, xlabel: 'x [m]', ylabel: 'η [m]' }).axes(); … },
    onPointer(g, ev, s) { if (ev.viewId !== 'graph' || ev.type !== 'pointerdown') return false; … g.app.set({probe: …}); return true; }
  },
  transport: { param: 't', min: 0, max: 20, rate: 2, end: 'hold' /* | 'loop' | 'stop' */, hold: 1.5, unit: 's' },
  modes: { param: 'system', options: [['wheel', 'Rotating wheel'], ['spring', 'Spring–mass']] },
  presets: [{ label: 'deep', set: { H: 200 }, title: 'kH ≫ 1' }, { label: 'shallow', set: { H: 1 } }],
  status: s => ({ text: s.k * s.H > 3 ? '🌊 deep water' : '🏖️ shallow water', tone: 'info' }),
  readouts: [{ id: 'c', label: 'Phase speed $c$', unit: 'm/s', value: s => c(s) }],
  terms: { title: 'Energy budget', unit: 'J/m²', items: [{ id: 'KE', label: 'kinetic', color: '#f97316', value: s => … }], total: { label: 'total' } },
  inspect: s => s.probe ? `$\\eta = a\\cos(kx-\\omega t) = ${Viz.tnum(a)} \\times ${Viz.tnum(cosv)} = \\mathbf{${Viz.tnum(eta)}}$ m` : null,
  notes: s => `<p>kH = ${Viz.fmt(s.k * s.H)}: …regime-dependent interpretation…</p>`,
  calc: { html: s => Viz.work.head('…') + Viz.work.line('$\\omega^2 = gk\\tanh kH$', Viz.fmt(w) + ' rad/s', 'why') + Viz.work.line('$t$', Viz.live('t')) + Viz.work.result('…'),
          live: s => ({ t: s.t.toFixed(2) }) },
  code: [{ id: 'disp', title: 'Dispersion in Python', ref: 'Eq. (7.36)', src: `omega = np.sqrt(g * k * np.tanh(k * H))   # = {{w}} rad/s`, live: s => ({ w: w(s).toFixed(3) }) }],
  explore: { intro: 'html with $tex$', controls: ['k', 'deep'], readouts: ['c'], callouts: [{ kind: 'try', html: '…' }] },
  tour: [{ title, text, set: { k: 1 }, play: true, controls: ['k'], readouts: ['c'], eq: 'disp', code: { id: 'disp', lines: [1, 1] },
           terms: true, inspect: true, notes: true, callout: { kind: 'key', html }, highlight: ['readout:c'], enter(app) {} }],
  equations: [{ id: 'disp', title: 'Dispersion relation', ref: 'Eq. (7.36)', tex: '\\omega^2 = gk\\tanh(kH)', live: s => …, note, symbols: [['\\omega', 'angular frequency', 'rad/s']] }],
  check: [{ q: '…', a: '…', set: { … } }],
  selftest: () => [{ name: 'omega(k=1,H=10)', js: omega({ k: 1, H: 10 }), py: 'ch07.omega(1.0, 10.0)', rtol: 1e-10 }]
});
```
- Inside callbacks use `g.app` (the outer `app` constant does not exist yet during the first draw). `g` is the first view;
  `g.view(id)` returns any view (`v.ctx, v.w, v.h, v.plot(), v.clear(), v.pointer`). `g.app.holding` is true while the
  transport holds at the end of a run (draw an end-of-run summary card then).
- Plot: `P.axes()`, `P.line`, `P.fn(f, {x0, x1, n, color, dash})`, `P.fillBetween`, `P.dot`, `P.arrow`, `P.text`, `P.clip`,
  `P.X/P.Y/P.ix/P.iy`. Fields: `Viz.field.grid/contour/drawContours/heatmap/streamline/quiver/particles`.
  Numerics: `Viz.num.linspace, trapz, rk4Step, odeint, brentq, erf, erfc, niceTicks, C`. `Viz.rng(seed)`.
- 3-D: `Viz.three(g.view('scene'), {distance, theta, phi, onPick(hit, ev)}).then(T => { if (!T) return; T.scene.add(mesh); T.render(); })`
  — drag orbits, wheel zooms, click picks; `T.project(vec3)` gives 2-D positions for labels on the view's 2-D canvas.
- Formatting: `Viz.fmt(v, {sig, unit})`, `Viz.tnum(v)` in TeX (round-off below 1e-12 prints as 0). Colours:
  `Viz.color('accent'|'teal'|'orange'|'rose'|'blue'|'amber'|'muted'|'text')`, `Viz.alpha(c, a)`; coloured words in HTML
  with `<b class="c-teal">`. In JS strings double every TeX backslash.

## 3. Fit rules (how to pass the audit without shrinking anything)
- The engine: picks the layout (`wide` | `portrait` | `landscape`), raises density 0–3 (hides help lines, `optional`
  items, the Explore intro, view titles on phones), pages every long list (Explore, Step by step, Equations, Code in
  5-line chunks, Check, and the walkthrough card itself: a step's extras continue on "›"), shortens tab labels, hides the
  transport's step buttons and speed on phones, and flags views smaller than 60 px.
- You: step text ≤ 45 words; ≤ 2 extras per step; ≤ 5 controls in Explore (mark the rest `optional`); readout labels ≤ 22
  chars; at most 3 views, `hidePortrait` on the least important; wide equations use `\begin{aligned}`; code lines ≤ 70 chars.
- Never add `overflow: auto`, fixed widths > 360 px or font sizes < 12 px in chapter CSS.
- Draw from `v.w`/`v.h`; `equal: true` for geometry; on small views reduce annotation density, not font size.

## 4. What a great explainer has (the quality bar — distilled from the references in §5)
1. **A question as the title** and a subtitle that says what is on screen.
2. **The phenomenon itself**, not just a chart — and ideally **several linked views of it** driven by one clock: the
   physical system, a graph over time/space, and a second representation (phase plane, response curve, profile, spectrum).
3. **A faint whole curve + a bold "so far" part**, a **ghost reference** (e.g. ω = 1, the undamped case, the inviscid
   limit) and **markers of the special points** (nodes, stagnation points, separation, critical values).
4. **Stages the reader can replay** — the walkthrough steps set the picture up; transport with play / step / scrub /
   speed; an **end-of-run summary** card.
5. **Show every number**: the Step-by-step tab works the formulas out with the current values ("formula = substituted =
   result — why"); an **inspector** traces the exact arithmetic for a clicked point.
6. **Show the code**: the Python behind the picture with live values in the comments, lighting up line by line in the
   walkthrough.
7. **Decompose**: term-by-term bars that add up to the total (click a term to isolate it).
8. **Presets at the special cases** and a **status verdict** naming the regime ("🎯 resonance", "⚠️ overturning").
9. **"Right now" interpretation** that changes with the regime, a small table of real-world values with the current row
   highlighted, and "Things to try".
10. **The same idea in another system** (modes/tabs: a wheel, a spring, a house's temperature — or for fluids: a pipe, a
    river, the atmosphere) when that makes the idea general.
11. Consistent colour meaning across text, bars and curves; a legend line; every number with a unit.

## 5. Reference explainers (open the one your storyboard names before building; match its depth)
| File | Study it for |
|---|---|
| `G:\Differential Equations_OCW\Course I am Completing\1_Introduction to Differential Equations\Unit 4_Exponential Response and Resonance\angular_frequency_explorer_1.html` | the gold standard: linked views (system animation · rotating pointer · 2π-window graph · response curves) on one clock, modes (wheel / spring / heat), presets, live "calc" working and regime-dependent notes with a highlighted table, ghost reference, end-of-run summary, fits 100dvh |
| `G:\Differential Equations_OCW\…\Unit 4_Exponential Response and Resonance\amplitude_phase_second_order_II_3.html` | system + graph + Bode + Nyquist windows linked by one state; crosshair readouts; a long live explanation that walks the complex-gain arithmetic step by step |
| `C:\Users\sislam27\Work\Climate Dynamics PHD\fast.ai\shammunul-fastai-notes\notebooks\interactive_viz\fid_formula_lab.html` | a formula whose terms are clickable, term-by-term bars and a total, stage chips, story text + code with the active line highlighted, a "make it perfect" preset |
| `…\interactive_viz\forward_noising_lab.html` | a closed-form jump shown as image + image = image, live coefficients, click a pixel to see its exact arithmetic, "variance is conserved" check line |
| `…\interactive_viz\pixels_as_parameters.html` | an algorithm in slow motion: forward → loss → gradient → step stages, step back/forward, a clickable loss curve to jump in time |
| `…\interactive_viz\random_copy_lab.html` | click-to-place interaction, before/after panels, code whose comments show the live random draws |
| `…\interactive_viz\overfitting_curves.html` | a minimal, very clear two-slider figure: dashed "future" remainder, a marker for the optimum, a verdict line that changes with the regime |
| `…\interactive_viz\stride_padding_playground.html` | presets that are the classic settings, the size formula with numbers plugged in, badges explaining the result, click an output cell to move the window |
| `…\interactive_viz\ddpm3_unet_3d.html`, `…\interactive_viz\np_resnet_3d.html` | three.js scenes: blocks you can orbit and click, a data token travelling through, pills to jump, an inspector card and code synced to the selected block, offline fallback note |
(`…` = `C:\Users\sislam27\Work\Climate Dynamics PHD\fast.ai\shammunul-fastai-notes\notebooks`. The published copies are at
`https://shammun.github.io/shammunul-fastai-notes/notebooks/interactive_viz/<file>`.) These references scroll or post their
height; ours must fit the window — take their *depth*, keep our layout.

## 6. Patterns that suit fluid mechanics
| Idea type | Views and features |
|---|---|
| a field (velocity, pressure, vorticity) | heatmap + streamlines + tracer particles; probe inspector with the velocity arithmetic; presets for classic flows |
| superposition / potential flow | toggles per element + the sum; stagnation points marked; terms for the complex potential parts |
| waves / dispersion | animated surface with particle orbits · ω(k) curve with a moving dot · a packet showing phase vs group speed; transport |
| boundary / similarity layers | profiles at several stations that collapse when rescaled (mode: raw / rescaled); status for separation |
| stability | growth-rate curve vs wavenumber · animated perturbation for the chosen k; status stable/unstable; presets at critical values |
| balances (geostrophy, Ekman, drag) | force arrows on a parcel with term bars that sum to zero; "turn off term X" presets |
| control volumes / budgets | a box with flux arrows, terms for inflow/outflow/storage, calc tab with the budget arithmetic |
| rotating / stratified (Ch. 13) | 3-D view of the spiral/wave · a hodograph · a profile vs depth; modes: ocean / atmosphere |
| numerics (Ch. 10) | grid + stencil view · solution vs exact · error vs h on log axes; transport steps the scheme |

## 7. Build → audit loop
`tools/new_viz.py` → write → `tools/viz_lint.py <file>` → `tools/shot.py <file>` (full: 8 sizes, every tab and step;
`--quick` while iterating) → **Read the PNGs** in `reports/viz/chNN/<slug>/` → compare with the reference → fix → repeat
(≤ 6 rounds).

## 8. Embedding (done by the notebook-builder and the publish tool)
- Notebook cell: `show_viz("chNN", "<slug>")` (via `nb.explainer(...)`). Local Jupyter/VS Code → `<iframe srcdoc>` of the
  file on disk; Colab → `<iframe src>` of the GitHub Pages copy (srcdoc from the clone if Pages is not live yet); the
  iframe is output-width × visible-window-height with Full screen / Open in new tab.
- Published page: `tools/publish_notebook.py` swaps each output for a `section.fp-viz` block sized to exactly 100 % of
  the window width and height; `tools/shot.py --page` verifies it. The GitHub `.ipynb` view shows a link card.

## 9. Lessons (the knowledge-keeper appends; one line each, dated by chapter)
- (machinery) The `#app` mount must be 100 % high or the no-overflow audit cannot see clipping; the audit also fails
  `app-does-not-fill-window`.
- (machinery) On 360×640 phones a multi-view stage needs the transport's step/speed controls hidden and the walkthrough
  card paged; the engine does both.
