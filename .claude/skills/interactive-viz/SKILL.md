---
name: interactive-viz
description: How fluidpy interactive explainers are designed, built, fitted to the window, verified and embedded - the explainer contract (5–10 per chapter), the assets/viz_lib.js API (Viz.app with linked views, transport, presets, status, terms, inspector, notes, the live Explain tab (explanation & interpretation), the step-by-step Derivation tab, synced code, modes, 3-D via three.js; Plot, field tools, particles, pager), the quality bar distilled from Shammunul's reference explainers (with their paths; the three Unit 4 mathlets are preferred), no-scroll layout rules, patterns for fluid phenomena, the lint/shot audit loop, and how show_viz and the publish tool put explainers full-window into notebooks, Colab and the web page. Load when choosing, storyboarding, building, reviewing or embedding an explainer.
---

# interactive-viz — explainers that fit the window and make an idea click

## 1. The contract (checked by `tools/viz_lint.py`, `tools/shot.py`, `tools/embed_check.py`)
- **5–10 per chapter** (`book.yaml → project.min/max_explainers_per_chapter`), each tied to CORE ideas of
  `analysis/chNN_curation.md`, plus one backup idea in the curation in case one fails review.
- Equations are **shown, not just numbered**: any tour, Explain, Derivation or quiz text that refers to a book equation writes
  it out (TeX) next to its number.
- One self-contained file `viz/chNN/<slug>.html`, created with `tools/new_viz.py`; `assets/viz_base.css` and
  `assets/viz_lib.js` are inlined between marker comments by `tools/viz_inline.py` (edit the assets, never the copies).
  External resources: only KaTeX and (for 3-D) three.js, both loaded by the library from CDNs with fallbacks.
- `<meta name="viz:*">`: chapter, slug, order, title, summary, concept, sections, equations, fluidpy, **derivations**
  (the curation's `D` ids this explainer steps through, e.g. `D03 D04`, or `none`).
- **Fits the window with no scrolling** at 360×640, 390×844, 844×390, 768×1024, 1280×720, 1000×700 (notebook),
  1366×768, 1920×1080 — fills the window exactly, nothing overflows, no text below 12 px, tap targets ≥ 24 px.
- Tabs: **Walkthrough** (4–8 steps) · **Explore** · **Explain** (required: "Explanation & interpretation" — ≥ 3
  numbered sections + an interpretation) · **Derivation** (required when `viz:derivations` lists ids: every page of every
  derivation is audited) · **Equations** (book-numbered, live) · **Code** (required) · **Check yourself** (≥ 3) ·
  optional custom panels.
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
  explain: { html: s => Viz.work.step(1, 'The dispersion relation with your numbers') +
                   Viz.work.line('$\\omega^2 = gk\\tanh kH = 9.81 \\times ' + Viz.tnum(s.k) + ' \\times \\tanh(' + Viz.tnum(s.k * s.H) + ')$', Viz.fmt(w) + ' rad/s', 'why') +
                   Viz.work.box('$c = \\omega/k = \\mathbf{' + Viz.tnum(c) + '}$ m/s') + Viz.work.step(2, '…') + Viz.work.step(3, 'At the current time') +
                   Viz.work.line('$t$', Viz.live('t')) + Viz.work.hint('On phones the ω(k) view is hidden; turn the phone sideways.') +
                   Viz.work.interpret('regime-dependent meaning of the current setting', 'Reading the current setting'),
             live: s => ({ t: s.t.toFixed(2) }) },       // ↗ opens it in its own live browser tab on big screens
  derivations: [{ id: 'D03', title: 'Where the dispersion relation comes from', short: 'ω(k)', ref: 'Eq. (7.36)',
                  goal: 'plain words', start: { tex: '…', plain: '…' }, plan: ['…'], uses: ['…'], view: 'curve', set: { … },
                  steps: [{ did: 'Substitute the wave form', tex: '…', why: 'allowed because … ; we do it to …', plain: '…',
                            live: s => 'the line with numbers', set: { … }, highlight: ['readout:c'], watch: '…' }],
                  result: { tex: '…', plain: '…' }, interpret: s => 'with your numbers …', check: 'units · limits · selftest' }],
  code: [{ id: 'disp', title: 'Dispersion in Python', ref: 'Eq. (7.36)', src: `omega = np.sqrt(g * k * np.tanh(k * H))   # = {{w}} rad/s`, live: s => ({ w: w(s).toFixed(3) }) }],
  explore: { intro: 'html with $tex$', controls: ['k', 'deep'], readouts: ['c'], callouts: [{ kind: 'try', html: '…' }] },
  tour: [{ title, text, set: { k: 1 }, play: true, controls: ['k'], readouts: ['c'], eq: 'disp', code: { id: 'disp', lines: [1, 1] },
           terms: true, inspect: true, notes: true, derive: { id: 'D03', step: 2 }, callout: { kind: 'key', html }, highlight: ['readout:c'], enter(app) {} }],
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
- Explain helpers: `Viz.work.step(n, title)` (numbered section) · `line(formula, value, why)` · `box(html)` /
  `result(html)` (boxed result) · `say(html)` · `hint(html)` · `interpret(html, label)` · `table(headers, rows, current)` ·
  `Viz.live(name)` (a value refreshed every frame from `explain.live`). Write it like the reference panels: what the views
  are (colours!) → each quantity computed from the controls → the values at the current time → what it means now.
- Derivation tab: pages are *the goal* (goal, start, plan, tools) · *step k* (we had → the move → now → Why → In words →
  your numbers → Watch) · *the result* (whole chain, result, what it means right now, check). A step's `set` moves the
  picture to the case being derived; on phones only `view` (else the first view) stays visible. A walkthrough step quotes
  one step with `derive: {id, step}` and links to the full derivation. Deep link: `#tab=derive&d=1&ds=3`.
- Pixel-space canvas helpers (after ch01): `Viz.font(px, weight, mono)`, `Viz.roundRect(ctx, x, y, w, h, r)`,
  `Viz.text(ctx, str, x, y, {size, weight, color, align, baseline, halo, bg, maxWidth})` (label on a solid box with
  `bg: true`), `Viz.card(ctx, cx, cy, lines, {color, colors, maxWidth, align})` (end-of-run card, wraps and stays inside
  the canvas).
- Formatting: `Viz.fmt(v, {sig, unit, keepTiny})`, `Viz.tnum(v, sig, {keepTiny})` in TeX (round-off below 1e-12 prints
  as 0 unless `keepTiny`), `Viz.fmtTime(seconds)` (µs … years). Colours:
  `Viz.color('accent'|'teal'|'orange'|'rose'|'blue'|'amber'|'muted'|'text')`, `Viz.alpha(c, a)`; coloured words in HTML
  with `<b class="c-teal">`. In JS strings double every TeX backslash.

## 3. Fit rules (how to pass the audit without shrinking anything)
- The engine: picks the layout (`wide` | `portrait` | `landscape`), raises density 0–3 (hides help lines, `optional`
  items, the Explore intro, view titles on phones), pages every long list (Explore, Explain, Derivation steps, Equations, Code in
  5-line chunks, Check, and the walkthrough card itself: a step's extras continue on "›"), shortens tab labels, hides the
  transport's step buttons and speed on phones, and flags views smaller than 60 px.
- You: step text ≤ 45 words; ≤ 2 extras per step; ≤ 5 controls in Explore (mark the rest `optional`); readout labels ≤ 22
  chars; at most 3 views, `hidePortrait` on the least important; wide equations use `\begin{aligned}`; code lines ≤ 70 chars;
  derivation lines short enough for a 1000×700 notebook iframe (the audit fails `equation-too-wide`) — split a long line
  into two steps or move the definition into *why*.
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
5. **Show every number, then say what it means**: the **Explain** tab ("Explanation & interpretation", the
   `forced_damped_vibrations.html` panel) works every displayed quantity out with the current values in numbered steps
   ("formula = substituted = result — why", results boxed, colours matching the curves), gives the values at the current
   time, and ends with a regime-dependent interpretation; an **inspector** traces the exact arithmetic for a clicked point.
5b. **Derive the hard formula in front of the picture**: for derivation-heavy ideas, the **Derivation** tab builds the
   result one move at a time (we had → move → now → why → in words → your numbers), each step moving the picture to what
   it is about, ending with the whole chain, a check and the interpretation.
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
**Shammunul's preferred models are the three Unit 4 mathlets** (★ below). Every explainer should feel like them: the
phenomenon and its graphs side by side on one clock, sliders and toggles that change everything at once, and a live
"Explanation & interpretation" that computes every number step by step and says what it means.
| File | Study it for |
|---|---|
| ★ `G:\Differential Equations_OCW\Course I am Completing\1_Introduction to Differential Equations\Unit 4_Exponential Response and Resonance\forced_damped_vibrations.html` | **the explanation panel to copy**: phase plane (drag the initial state) + x(t) graph with steady state / transient / solution toggles (green / blue / orange) on one time slider; the "+ explanation" panel computes ω₀ and ζ with numbers, the regime and characteristic roots, p(iω) → amplitude → phase lag in boxes, the transient from the initial condition (half-life, 1 % time), the values "at the current time" (x = steady + transient), what each window shows, and a regime-dependent interpretation (below / near / above resonance); floats, minimises and pops out into a live tab |
| ★ `G:\Differential Equations_OCW\…\Unit 4_Exponential Response and Resonance\angular_frequency_explorer_1.html` | linked views (system animation · rotating pointer · 2π-window graph · response curves) on one clock, modes (wheel / spring / heat), presets, live working and regime-dependent notes with a highlighted table ("From ω to period", "Right now", "The 2π puzzle", "Things to try"), ghost reference, end-of-run summary, fits 100dvh |
| ★ `G:\Differential Equations_OCW\…\Unit 4_Exponential Response and Resonance\amplitude_phase_second_order_II_3.html` | system + graph + Bode + Nyquist windows linked by one state; crosshair readouts; the explanation is a **numbered derivation with live numbers** (1. evaluate p(iω) · 2. complex gain via the conjugate · 3. amplitude and phase · 4. period and time lag · 5. reading the current setting · 6. at time t) and each optional window gets its own section, or a hint to open it |
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
| control volumes / budgets | a box with flux arrows, terms for inflow/outflow/storage, Explain tab with the budget arithmetic, Derivation tab from the integral to the differential form |
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
- (machinery) Seven tabs squeeze the title on a 768 px tablet: the header now falls back to short labels, then to a tab
  row of its own. On phones the Derivation tab keeps one view and hides the preset strip and transport.
- (machinery) Derivation lines wider than the side panel of a 1000×700 notebook iframe fail `equation-too-wide`: keep a
  line to one relation; put definitions (`ω_d = …`) in their own step or in *why*.
- (ch01) Never assign `ev.view` in a pointer handler (read-only `UIEvent` getter; it threw in every `onPointer`); use
  `ev.viewId` / `ev.vizView`. `shot.py` does not click yet: run a click/drag pass over every view and count page errors.
- (ch01) Two sign conventions: the status badge reads verdict word → the fluidpy criterion text in the chosen
  convention → the other convention's bare relation, at every size; pin it with exact-text selftest rows (E4).
- (ch01) A view hidden on portrait phones must not carry the step's key number: repeat it in a visible view's title or
  a readout ("Kn(body) 0.064 → slip flow"), and never tell phone readers to use controls that are hidden there.
- (ch01) A step that traces one sample's arithmetic sets the global parameter to that sample's value, and the badge,
  equation card and inspector compute the same quantity on the same basis (mean vs crest count disagreed in E1).
- (ch01) The Derivation tab has the notebook's step count and move titles; when either side splits a step, renumber
  every walkthrough `derive: {id, step}` link in the same edit (reviewer compares `app.cfg.derivations` with the builder).
- (ch01) Plain-text *why*/*watch* strings never contain raw TeX (write e^(±λt)); *why* ≤ 35 words.
- (ch01) Narrow views (< 420 px): draw the rotated y title in its own strip (the library clamps the left pad to 40 px
  and the title covers minus signs of ticks); move legends into the view title; below 60 px draw only the 0 tick.
- (ch01) Use `Viz.card` for end-of-run cards (wraps inside the canvas; a one-line card clipped on phones), `Viz.text`
  with `bg: true` for labels crossing curves, `Viz.fmtTime` for durations, `Viz.tnum(v, sig, {keepTiny: true})` for
  molecular-scale numbers. Units inside TeX use `m^{3}`, not "m³" (KaTeX warns).
