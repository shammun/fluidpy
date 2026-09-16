# Explainer patterns — what worked, what failed, what to promote

Read by the curator, designer, viz-builders and viz-reviewer before any explainer work. Appended by the knowledge-keeper.
Explainer short names: ch01 E1 `continuum_averaging_volume`, E2 `viscosity_momentum_diffusion`, E3 `heat_work_paths`,
E4 `parcel_stability`, E5 `buckingham_pi_machine`; ch02 E1 `rotation_of_axes`, E2 `cauchy_traction_principal_axes`,
E3 `strain_vs_rotation_split`, E4 `gauss_flux_box`, E5 `stokes_circulation_loop` (all PASS round 2, 127 parity rows).

## Patterns that worked
| Pattern | Where proven | Why it works |
|---|---|---|
| Walkthrough step that sets the parameters itself, then shows one inline slider ("Drag U and watch…") | `templates/viz_example.html` step 3 | the reader sees the change before being asked to act |
| Live equation card in a walkthrough step (`eq:` + `live:`) | `templates/viz_example.html` step 4 | symbols become numbers the reader controls |
| Tracer particles + streamlines + colour field on one stage | `templates/viz_example.html` | motion shows what a static streamline plot hides |
| Draggable probe showing the local velocity vector | `templates/viz_example.html` | direct manipulation of "the field at a point" |
| Explain tab in numbered sections: what the colours are → each quantity from the controls (boxed results) → at the current time → "Reading the current setting" | `templates/viz_example.html` (after `forced_damped_vibrations.html`); ch01 E1–E5 (7–8 sections each) | every number on screen is accounted for, then interpreted for the regime the reader chose |
| Derivation step that moves the picture (`set`) and shows the line with the reader's numbers (`live`) | `templates/viz_example.html` derivations `energy` step 4, `omegad` step 2 | the symbols in the derivation become the numbers and curves on screen |
| Walkthrough step quoting one derivation step with an "All steps →" button | `templates/viz_example.html` tour step 4; ch01 E3 → D14 step 10, E4 → D18 step 11, E5 → D28 step 12 | the story stays short while the full derivation is one click away |
| **Two-convention status badge**: verdict word first (with emoji), then the library's exact criterion text in the chosen convention, then the other convention's bare relation ("🌊 stable · stable ⇔ dT/dz > Γa: −6.5 > −9.8 K/km · met: 6.5 < 9.8 K/km"); a selftest row rebuilds the badge from `cfg.status` and compares it with the fluidpy text | ch01 E4 (48/48 browser cases match `lapse_rate_stability`, 13 exact-text parity rows) | both sign conventions are on screen at every size and every step; the verdict never reads as the opposite of the inequality text |
| **Exact-text parity rows**: when the Python function returns a sentence (a criterion with numbers), the JS mirror's string is compared byte for byte, not only its numbers | ch01 E4 | catches rounding, minus-sign glyph and inequality-direction drift that a numeric row misses |
| **Separate y-title strip on narrow views**: below 420 px shift the plot rect 18 px right and draw the rotated y title in that strip | ch01 E4 (`vplot`) | the library clamps the left pad to 40 px there, and a rotated title otherwise covers 4-character ticks ("−150" read as "150") |
| **Put a hidden view's key number in a visible title**: on portrait phones the Kn strip is hidden, so the curve title reads "Kn(body) 0.064 → slip flow" and the step adds a Kn readout; on landscape phones the strip gets its own row | ch01 E1 step 6 | "turn your phone" hints are not an explanation; the number the step is about is always visible |
| **A trace step sets the global parameter to the sample's own value**: the step that follows one sample's arithmetic sets L to that box (21 nm), so badge, Eq. card and inspector all describe the same box | ch01 E1 step 4 | two different boxes on one screen (the round-1 bug) make the arithmetic untraceable |
| Steady-state ghost + fading trail of earlier profiles + a thin dashed exact series on top of the numerical profile | ch01 E2 | shows the transient, its end state and the correctness check in one picture |
| Modes that change only the diffusivity (momentum ν / heat κ / species κ_m) on the same stage | ch01 E2 | "same equation, different D" is the idea; switching modes proves it |
| Linked p–v and T–s diagrams with shaded work and heat areas (negative areas hatched), plus term bars q, w, Δe and Δs vs ∫δq/T | ch01 E3 | path functions change area while the state functions' bars stay put; irreversible mode separates Δs from ∫δq/T |
| Draw only what is defined: no moving state dot on a non-equilibrium (stirring, free expansion) route | ch01 E3 | a dot would claim equilibrium states that do not exist |
| On narrow views move the legend into the view title | ch01 E3 p–v view | frees the plot area; nothing overlaps the route |
| Labels that cross curves sit on solid card-coloured boxes (now `Viz.text(…, {bg: true})`) | ch01 E4 ΔT/Δρ labels | readable over any curve without shrinking the font |
| A derivation step that moves the picture to the special case it treats (double root: Γ = Γa, flat ζ(t)) | ch01 E4 D18 step 14 | the degenerate case is seen, not just stated |
| Unit-system switch (SI → cgs → imperial) on a table where every raw value changes and every Π value stays | ch01 E5 Equations tab | invariance is watched, not asserted |
| Exact rational arithmetic in JS (small `Frac` class mirroring `fractions.Fraction`) for the exponent solve | ch01 E5 | parity with `pi_groups` is exact, fractions print as ½, not 0.49999 |
| Presets at the classic cases: dry adiabatic (neutral), standard atmosphere, nocturnal inversion, superadiabatic surface layer, ocean thermocline, salt tank; fluids at their temperatures; six Π problems | ch01 E4, E2, E5 | each preset is a regime the text talks about |
| Pointer interaction through `stage.onPointer` (drag a gas state, click a diagram to probe) | ch01 E3 round 2 | no local capture-phase workaround; the engine routes `ev.viewId` |
| **Trace steps hold still**: `autoplay: false` in the app config, and `play: false` on every walkthrough step that quotes numbers; only a step that says "press ▶" plays; a derivation's result page gets `play: false` too | ch02 E1 (round-1 Must), E2 (D17 result page), E4 | a running transport overwrites the value the step just set, so the text and the picture disagree (1.866 in the text vs 2.215 on screen) |
| Long live equation lines split with `\begin{aligned}` (one relation per row) | all ch02 explainers (Equations tab, Explain) | fits the 1000×700 notebook frame without `equation-too-wide` |
| **One view row on landscape phones**: scoped `<style data-chapter>` hides `.viz-views > .viz-row:nth-child(2)` and the preset strip under `[data-layout="landscape"]` (and under `@media (max-height: 700–720px)` for short portrait phones), with `:not([data-tab="derive"])` so the Derivation tab keeps its own view; the hidden row's key number moves into the kept view's title | ch02 E3, E4, E5 | the engine has no per-row hide flag yet; the scoped CSS keeps every view ≥ 60 px without touching the library |
| **Shrink transport**: the transport runs a parameter σ = −log₁₀h, so "play" shrinks a box or loop geometrically toward a point; the limit view uses a log-h axis | ch02 E4 (`hOf = 10^−shrink`), E5 | a linear sweep spends almost all its time at large h; the limit is only visible on a log scale |
| **Shortened chip labels** for fields and presets ("(x², 0)", "b×x", "vortex") with the full name in the title/tooltip; `.viz-seg` chips wrap (`flex-wrap: wrap`) when there are six | ch02 E4, E5 | six field chips fit on a phone without overflow |
| **Diverging colormap with a white zero** (blue–white–orange for divergence, purple–white–teal for curl), limits symmetric ±max | ch02 E4 (`Viz.CMAPS.bwo`), E5 (`curlImage`) | the sign of ∇·u or ∇×u reads at a glance; "zero" is visibly empty |
| **Counter-example table**: four candidate triples, each with its (2.8) residual under the same C, plus a "fixed array (not a tensor)" preset whose (2.12) residual is in the status badge | ch02 E1 | the definition is taught by what fails it, by order one, not by round-off |
| **Passive/active table with the current row lit** and a toggle labelled "Wikipedia R" that computes Cx | ch02 E1 | the classic confusion is settled with a picture; `activeRotate(θ) = transformVector(−θ)` pinned by a parity row |
| Clickable matrix cells with an inspector that draws the two unit vectors and their angle (cos of the angle = the cell) | ch02 E1 | "what does one number C_ij measure?" is answered by pointing |
| A mirror toggle wired to a derivation step (det C = −1 while CᵀC − I stays 0; status turns amber) | ch02 E1, D02 step 8 | the step's special case (reflection) is seen |
| **Element + decomposition + term bars**: the removed half greyed, n, f and its blue normal / rose shear parts; σn = τ₁₁n₁² + (τ₁₂+τ₂₁)n₁n₂ + τ₂₂n₂² as bars that add exactly | ch02 E2 | the stress tensor is visible as "what pushes on a cut" |
| Three linked views on one angle φ: element · σn/τs curves with λ bounds and special-cut ticks · Mohr circle with the double angle | ch02 E2 | principal and max-shear cuts are the same event in three pictures |
| A 2-D sketch that shrinks with the derivation step (h = 1 → 0.5 → 0.25) while face arrows stay and the dashed volume term shrinks | ch02 E2 D05 | the orders-of-smallness argument is seen, not stated |
| Live numbers on the derivation steps where numbers help (6 of 15 in D17: discriminant, b·b = 0, τ′ = diag, c_k, σn = Σλc², shear bound) | ch02 E2 D17, E4 D22 (5 of 9) | a ★★★ derivation stays concrete |
| **Contrast toggle for an index convention**: non-symmetric τ draws τ·n in amber next to n·τ; Mohr says "not a circle" | ch02 E2 | the first-index rule becomes a visible difference (√3) |
| Adaptive status length (phone: verdict + one number; landscape: medium; desktop: full) | ch02 E2 | no three-line badges on phones |
| Mode-dependent axis names (stress: σn, τs [Pa]; strain rate: n·S·n, s·S·n [1/s]) | ch02 E2 (round-1 Should) | the same stage serves two tensors without mislabelled units |
| **One clock, three squares**: G leans, S (teal) stretches along ±45°, A (orange) spins; walkthrough G → split → S → A → all → blind → your turn | ch02 E3 | the cleanest story of ch02: the decomposition is watched, not asserted |
| Closed-form matrix exponential `expm2` (cosh / cos / shear-limit branches) matched to scipy `expm` at three regimes (1e-10) + tracer inspector showing e^{Mt} and x(t) under G, S, A | ch02 E3 | exact trajectories of linear flows in JS; reusable for Ch. 3 pathlines |
| Invariant pinned as a bar: S:A always 0 while S:S + 2S:A + A:A = G:G | ch02 E3 | Eq. (2.29) as a live invariant; name the factor 2 in the bars title |
| Conventions pinned by parity rows ("vorticity ω₃ = vector(R)", "R = 2A") and by one readout text ("Spin ½ω₃ = A₂₁") everywhere | ch02 E3 | the factor-2 trap cannot reappear silently |
| **Waterfall bars** left → bottom → right → top → ∮ next to ∬, with an interior-faces bar when tiled; face arrows coloured by sign (out orange, in blue) | ch02 E4 | Gauss' theorem as bookkeeping the reader can watch |
| Tiling with shared faces drawn as opposite pairs + a shared-face inspector ("+0.250 … from the left tile, −… from the right, sum = 0 — D25 step 8") | ch02 E4, E5 (edges) | the cancellation step of the proof is clickable |
| Limit view with three regimes in one picture: exact at every h (polynomial), slope-2 log–log inset (smooth), 1/h² blow-up (singular source) | ch02 E4 | "limit", "order" and "hypothesis" in one plot |
| **A theorem failing on purpose**: vortex core inside the loop → Γ = 2πK vs curl flux 0, bars "core inside: sides differ by 6.283", status ⚠️, `hypothesis_ok = False` parity row | ch02 E5, E4 (point source) | hypotheses are taught by breaking them |
| "Unrolled loop" view: u·t against arc length, shaded area = Γ, split into Ex. 2.6's four sides as (hatched when negative) term bars; loop-point inspector for one of 400 terms | ch02 E5 | a line integral becomes an area |
| Paddle wheel turning at ½(∇×u)₃ on a curl heat image, plus RK4 tracers; shear preset shows "curl without curved streamlines" | ch02 E5 | the curl's meaning in one glance |
| Staircase → circle refinement in a derivation (8×8 tiles, gap printed live) | ch02 E5 D26 step 9 | the tiling limit is measured, not claimed |
| Flip-n toggle that changes both sides of the theorem together, with the status saying so | ch02 E5 | orientation is a convention, the equality is not |

## Failures and fixes
| Problem | Where | Fix |
|---|---|---|
| App shell not 100 % high → panels clipped on phones but the audit saw no overflow | engine, before chapter 1 | `#app` height 100 %; audit fails `app-does-not-fill-window` |
| A walkthrough step with text + controls + readouts + equation overflowed on 360×640 | example step 5 | ≤ 2 extras per step |
| Explore intro + callouts pushed the controls off a phone screen | example | callouts optional by default; intro hidden at density 3 on phones |
| 7 tabs squeezed the title into a one-letter column on a 768 px tablet (no overflow, but views shrank below 60 px) | example, Explain + Derivation added | header falls back to short labels, then to a tab row of its own |
| Derivation tab on phones left 19 px views for a 3-view stage | example | phones keep one view (`derivation.view`) and hide the preset strip and transport on that tab |
| A derivation line "r = … = …" and a long live line overflowed the side panel in the 1000×700 notebook frame | example `omegad` step 4, `energy` step 4 | one relation per line; definitions move to their own step or into *why* |
| Term bars in a walkthrough step showed "–" | engine | term rows index their own item (bars from several blocks) and stale rows are dropped |
| **Every `onPointer` handler threw** (`ev.view = v` assigns a read-only `UIEvent` getter in strict mode); `tools/shot.py` never clicks, so all explainers passed | engine, found during ch01 E1–E3 | lib sets `ev.vizView` and shadows `view` with `defineProperty` (b45cbef); the viz-reviewer's 2 141-action click/drag pass is now the model; machinery TODO: shot.py must click/drag every view and fail on page errors |
| A view the walkthrough talks about is hidden on portrait phones ("turn the phone sideways") | ch01 E1 Kn strip (round-1 Must) | own row on landscape phones; key number in a visible title on portrait phones |
| One step showed two different sampling boxes (badge vs inspector) | ch01 E1 step 4 (round 1) | the step sets the global L to the traced box |
| Badge and inspector use different bases for the same quantity (mean count vs crest count: ±6.7 % vs 0.0607) | ch01 E1 (open) | compute both from the same count, or label the badge "at the mean density" |
| The other sign convention was not always on screen | ch01 E4 (round-1 Must) | badge carries both conventions at every size and step; wraps to two lines on phones |
| Unstable badge led with "stable ⇔ …" (the criterion text) | ch01 E4 (round-1 Must) | verdict word first, criterion second |
| Rotated y title covered minus signs of ticks on portrait phones | ch01 E4 (round-1 Must) | separate y-title strip (pattern above) |
| Tick labels collide with a marker in a 45 px ζ(t) view at 360×640 when the badge wraps | ch01 E4 (open) | below 60 px draw only the 0 tick or move ±ζ_max into the title |
| Library 3-significant-figure tick labels repeat on narrow ranges (301.9 … 302.3 K) | ch01 E4 | local `xTicks` with decimals from the tick step (promotion candidate below) |
| Derivation tab step counts drifted from the notebook after the notebook split steps (D14 10 vs 11, D18 13 vs 14, D28 12 vs 13) | ch01 E3, E4, E5 (round 1) | reviewer compares `app.cfg.derivations` step counts and move titles with `notebooks/build_chNN.py`; renumber walkthrough deep links in the same edit |
| *Why* text too long (46 and 53 words; limit 35) | ch01 E3 D10 step 7 (fixed), E5 D28 step 9 (open) | keep the rule in *why*; move the consequence into *in words* or *watch* |
| Raw TeX ("e^{±λt}") in a plain-text *why* | ch01 E4 D18 step 13 (open) | plain text uses "e^(±λt)"; TeX only inside `$…$` |
| WATCH line pointed at a view hidden on phones; "switch the medium" hint where the chips are hidden | ch01 E4 (round 1) | name the visible substitute ("the coloured bracket in T(z)", "use the thermocline preset") |
| Code tab used names it never defined (`eps`, `L_flow`, `T`, `p`) | ch01 E1 (round 1) | every name in a snippet is assigned in the snippet, with its live value in the comment |
| Legend overlapped the route on a phone p–v view | ch01 E3 (round 1) | legend in the view title |
| One-line end-of-run card ran past the right edge of a portrait-phone view and covered a state label | ch01 E3 (open) | `Viz.card` (wraps within the canvas) or a shorter card below 420 px |
| KaTeX warnings from Unicode superscripts ("m³") inside TeX | ch01 E5 (round 1) | map units to TeX powers (`m^{3}`) before rendering |
| SI / cgs / imperial columns dropped in the 1000×700 notebook frame | ch01 E5 (round 1) | narrower column labels + footer "columns: SI · cgs · imperial" |
| Clicking a starred (solution) header added it to the custom repeating set | ch01 E5 (open, optional) | ignore header clicks on the solution variable |
| Eq. card value lags the canvas while playing | ch01 E2 step 3 (open) | refresh the live equation from the same frame state as the canvas |
| `Viz.fmt`/`Viz.tnum` printed molecular masses (~5e-26 kg) as 0 | ch01 E1 | local `fmtK`/`T`; library now has `keepTiny` in both |
| **Transport autoplayed under trace steps**: the step text quoted 1.866 at θ = 30° while the plane showed θ = 55° | ch02 E1 (round-1 Must) | `autoplay: false` + `play: false` on every trace step (pattern above); the inspector re-renders or pauses on cell click |
| Labels collided at the default angles (arc label on shadow feet, θ label on x₁, primed labels on one face) | ch02 E1 (round-1 Must) | inspector arc label across the origin, θ label beyond the arrowhead, shadow labels offset perpendicular to their axis, primed normal/shear labels on different faces |
| Unprimed "τ₁₂ = 1.000" label crossed by its own rose arrow in tensor mode | ch02 E1 (open, cosmetic) | offset along the face or put it on a `bg: true` box |
| Third matrix-assembly line clipped at the canvas bottom on 360×640 and 844×390 | ch02 E1 | below ~180 px keep one fitted line (the rest is in the badge) |
| The bars gap label gave the wrong reason in the singular regime ("midpoint rule ∝ 1/n²" for a delta-function discrepancy) | ch02 E4 (round-1 Must) | regime-dependent label: "the delta at the origin carries the flux m — Gauss needs a smooth Q" |
| A derivation step's key number was cut off in a truncated phone title (and the limit view is hidden there) | ch02 E4 D21 step 4 (round-1 Must) | lead the narrow title with the key number; add a `live` line to the step |
| A walkthrough step said "the box shrinks" while the sweep first grew it (start at h = 2.5) | ch02 E4 | start the sweep at the stated value |
| Portrait field views tiny (square ~40 px, box ~70 px) or with a wide x-range (±5/±4) because equal-aspect plots in a short row | ch02 E3, E4, E5 (partly open) | `rows: [2, 1]` on portrait or clamp the range; let the heat image fill the width |
| End-of-run card covered axis labels of a panel | ch02 E3 | place it in the empty title band or a free corner |
| Single-panel mode used a third of a wide view | ch02 E3 | `pw = min(v.w, v.h·1.6)` |
| A step referred to a view hidden on portrait phones ("the dot on the flat curve") | ch02 E5, E1 step 7 | add the readouts to the step ("watch the Γ/A readout …"); repeated ch01 lesson |
| D26 step 9 showed a 13 % staircase-vs-circle gap without saying why | ch02 E5 | use the finest tiling and print the gap live |
| *Why* over 35 words (D01 36, D02 40, D17 39–42) | ch02 E1, E2 | move the consequence into *in words* (repeated ch01 lesson) |
| **Two builders sharing one scratchpad clobbered each other's `splice.py`** | ch02 viz phase (parallel builders) | every builder uses a private subfolder (`<scratchpad>/<slug>/`) for helper scripts |
| `shot.py` `py:` parity expressions have **no builtins** (`abs`, `len`, `float` fail) | ch02 all | index dicts/tuples down to one float in the expression; put `abs` on the JS side or return the magnitude from fluidpy |
| Library 3-significant-figure ticks and manual log ticks written twice | ch02 E4 `drawLimit`, E5 `drawBars` | promotion candidate `P.axes({xlog: true})` |
| `test_viz_library_inlined_and_template_lints` failed while builders were mid-build ("STALE viz/ch02/…") | ch02 verify phase | expected while the viz phase runs; re-run `tools/viz_inline.py --all` before the merge gate |

## Promotion candidates (helpers duplicated across explainers)
| Helper | Found in | Proposed library name | Status |
|---|---|---|---|
| TeX number with mantissa renormalised (10×10⁻⁵ → 1×10⁻⁴) | ch01 E1 `T`, E2 `tn`, E5 `tn` | `Viz.tnum` fix | **promoted** (ch01 knowledge pass); + `keepTiny` option for E1 |
| canvas font string in the page font | ch01 E1 `font`, E5 `font` (inline in E2, E3, E4) | `Viz.font(px, weight, mono)` | **promoted** |
| rounded-rectangle path | ch01 E5 `rrect`, E2 inline `roundRect` fallback | `Viz.roundRect(ctx, x, y, w, h, r)` | **promoted** |
| pixel text with halo / background box | ch01 E3 `canvasText`, E4 `drawText` | `Viz.text(ctx, str, x, y, opt)` | **promoted** |
| end-of-run summary card | ch01 E1 (inline rect card), E2 `card`; E3 needs it (clipped card) | `Viz.card(ctx, cx, cy, lines, opt)` | **promoted** (wraps, stays in canvas) |
| duration in everyday units | ch01 E2 `humanTime`, E4 `fmtT` (s/min) | `Viz.fmtTime(seconds)` | **promoted** |
| exact rationals | ch01 E5 `Frac`/`Q` | `Viz.Frac` | candidate — promote when Ch. 4 §4.11 (similarity) needs it in a second explainer |
| Poisson / normal / lgamma samplers matching numpy | ch01 E1 | `Viz.num.poisson`, `Viz.num.normal`, `Viz.num.lgamma` | candidate — one explainer so far (Ch. 12 turbulence statistics may reuse) |
| plot with a y-title gutter on narrow views | ch01 E4 `vplot` | `Plot` option `{ylabelGutter: true}` (or default below 420 px) | candidate — really a library defect; fix in the engine when a second explainer hits it |
| x ticks with decimals from the tick step | ch01 E4 `xTicks` | `P.axes({xdecimals: 'auto'})` | candidate — engine defect on narrow ranges |
| length in everyday units (pm … km) | ch01 E1 `lenLabel` | `Viz.fmtLength(m)` | candidate |
| number keeping trailing zeros (1.00) | ch01 E2 `fm`, `tn` | `keepZeros` option for `Viz.fmt`/`Viz.tnum` | candidate |
| diagonal hatch pattern for negative areas | ch01 E3 `hatchPattern` | `Viz.hatch(ctx, color)` | candidate |
| regime colour/word maps (stable teal / neutral amber / unstable rose) | ch01 E4 | keep local; record the colour meaning in the chapter | not a library helper |

The ch01 explainers still carry their local copies; the next rebuild of any of them should switch to the library
helpers (and E3's end card to `Viz.card`, which also closes its open Should-fix).

### ch02 candidates (listed, **not promoted**: the ch02 knowledge pass ran while the site-publisher was reading `viz/`)
Promote in a pass where nothing else reads `viz/`: edit `assets/viz_lib.js` (documented, backwards compatible) →
`tools/viz_inline.py --all` → `tools/shot.py --chapter ch01 --quick`, `--chapter ch02 --quick` and
`tools/shot.py templates/viz_example.html --quick` must PASS, otherwise revert.

| Helper | Found in (file · local name) | Proposed library name | Priority |
|---|---|---|---|
| fixed-decimal formatters, HTML minus vs TeX minus, tiny → 0 | ch02 E2 `fx`/`tx`, E3 `f2z`/`tz`, E4 `fz`/`tz`, E5 `fz`/`tz` (4 copies) | `Viz.fixed(v, d, {tex})` | **high** (4 explainers) |
| 2×2 symmetric principal values/angle | ch02 E2 and E3 `principal2d` | `Viz.num.principal2d(S)` | high (2 copies) |
| matrix cells with row/column highlight and fitted assembly lines | ch02 E1 `drawMatrix`, E3 `matrixLayout`/`drawMatrix`; ch01 E5 matrix view | `Viz.matrix(ctx, rect, M, {hot, colors, fmt})` | high (3 explainers) |
| log-x axis with manual ticks | ch02 E4 `drawLimit`, E5 `drawBars` | `P.axes({xlog: true})` | medium |
| hatched (negative) bars | ch02 E5; ch01 E3 `hatchPattern` | `Viz.hatch(ctx, color)` | medium (2 chapters) |
| tiling helpers (rectangle or disc into k×k tiles, per-tile sums) | ch02 E4 `tiled`, E5 `tilesOf`/`tileSums` | `Viz.num.tiles(shape, k)` | medium |
| diverging colormap with white zero | ch02 E4 `Viz.CMAPS.bwo` (local extension), E5 `curlImage` | `Viz.CMAPS.bwo` + `Viz.field.heatmap({diverging: true})` (vmin = −vmax) | medium |
| text clamped inside the plot rectangle | ch02 E4 `labIn` | `Viz.textIn(P, str, x, y, opt)` | low (1 copy) |
| angle arc with a label | ch02 E1 `drawArc` | `P.arc(a0, a1, r, {color, label, labelRadius})` | low |
| label pushed past an arrow tip | ch02 E1 `tipLabel`, E2 inline `tipLabel` | `P.tipLabel(from, to, str, {color, offset})` (needs collision-avoiding offset) | low |
| try progressively shorter strings until one fits | ch02 E1 (inside `drawMatrix`) | `Viz.fitText(ctx, alternatives, maxW)` | low |
| filled polygon in plot coordinates | ch02 E2 `fillPoly` (+ `halfSquare`) | `P.fillPoly(pts, color, alpha)` | low |
| closed-form 2×2 matrix exponential; ε array | ch02 E3 `expm2`, `eps3` | `Viz.num.expm2(M, t)`, `Viz.num.eps3` | low now; **Ch. 3 pathlines will want it** |
| `Viz.card` with `align: 'left'` | ch02 E3 | option on `Viz.card` | low |
| adaptive status (short / mid / full by width) | ch02 E2 | engine: `status: s => ({short, mid, text, tone})` | engine request |
| `terms.unit` as a function of state | ch02 E2 (stress Pa vs strain 1/s) | engine: `terms.unit: s => …` | engine request |
| readout label wrap CSS | ch02 E2 | engine CSS | engine request |
| `.viz-seg` wrap when there are many chips | ch02 E5 | engine CSS default | engine request |
| hide a view row on landscape phones (and on short portrait phones), Derivation tab exempt | ch02 E1 (request), E3/E4/E5 scoped CSS | engine: `views[i].hideLandscape` / `rows` flag | engine request |
| `Viz.num.poisson/normal/lgamma`, `Viz.Frac`, y-title gutter, `xdecimals: 'auto'`, `fmtLength`, `keepZeros` | ch01 (see table above) | unchanged | carried over |

Python promotion candidate (same rule, `fluidpy/core/`, then `pytest -q`): `fluidpy.core.anim` could save animations
inside `matplotlib.rc_context({"figure.constrained_layout.use": False})`. matplotlib 3.11's `print_figure`
re-installs constrained layout after the first saved frame, so builders currently toggle the rcParam by hand around
`subplots_adjust` animations (`notebooks/build_ch02.py` l. 3482/3507). The ch02 kinematics helpers
(`velocity_gradient_preset`, `linear_flow_map`, `deform_square`, `material_line_angle`) move to `core/` when Ch. 3
calls them.

### Lesson candidates for the skills (not promoted in this pass; one line each for the next pass that may edit skills)
- `interactive-viz` Lessons: (ch02) `autoplay: false` + `play: false` on every step that quotes numbers ·
  `\begin{aligned}` for long live lines · scoped `<style data-chapter>` to drop a view row on phones with
  `:not([data-tab="derive"])` · `py:` parity expressions have no builtins · parallel builders use private scratch
  subfolders · a regime-dependent label must give the regime's reason (singular vs quadrature).
- `math-to-python` §7: (ch02) a test asserting f(B) == g(B) between two code functions is not evidence for the book's
  convention; pin to a field with a known answer · build sympy expressions from the library's own `coordinates`
  (assumptions make different symbols) · normalise symbolic vectors with `b/sqrt(b·b)`, not `.norm()` · theorems whose
  hypotheses can fail return a flag instead of asserting · preset names are physics: test them by their defining
  property.
- `verify-implementation`: (ch02) patch every pinned convention with its wrong variant and show a test fails
  (8/8 caught in ch02) · report `hypothesis_ok` rather than asserting through a singularity.
- `teaching-style` Lessons: (ch02) reorder the book when a derivation needs it, and say so ((2.15) before (2.12)) ·
  counter-example residual tables for "is it a vector/tensor?" · assert every symbolic twin against the numeric
  result in the same cell · look at every animation's last frame (clipped readouts appear when numbers grow) · demo a
  Python idiom on an object where a wrong statement fails · the reusable derivation moves listed in
  `knowledge/concept_map.md` → "Reusable derivation moves".
- `colab-notebook` / nbkit: (ch02) `nb.derivation` must not prefix "Eq." to labels like "Exercise 2.8"; design Part E
  reminder rows need a marker so `tools/coverage_check.py` stops raising 28 false name-matching warnings.

## Ideas for later chapters
(seeded from `book.yaml → viz_seeds`; the curator decides)
- **Ch. 2 tensors** (done): the ch01 E5 matrix layout became ch02 E1's clickable direction-cosine matrix; the
  arrow-on-a-plane idea became ch02 E2's element + Mohr stage.
- **Ch. 3 kinematics**: streamline/pathline/streakline stage = template field stage + ch01 E2's "ghost + trail" idea;
  for linear flows use ch02 E3's `expm2` closed form for exact pathlines. **Deformation of a fluid element (§3.4)** =
  ch02 E3's "one clock, three squares" with the book's R = G − Gᵀ and ω = ∇×u. Keep the factor-2 readout
  ("spin ½ω₃") and the parity rows "vorticity = vector(R)", "R = 2A". Principal strain rates = ch02 E2 in strain mode.
  Reynolds transport theorem = ch02 E4's box + waterfall bars with a moving boundary.
- **Ch. 4 conservation laws**: dynamic similarity explainer = E5 with presets for Re, Fr, Ro and a model-vs-prototype
  table (promote `Viz.Frac` then); Boussinesq buoyancy reuses E4's column + profile and `brunt_vaisala_sq`.
- **Ch. 4 stress and Cauchy's equation**: ch02 E2 (element + traction split + Mohr), now with the symmetry of τ proved;
  control volumes reuse ch02 E4's face bookkeeping.
- **Ch. 5 circulation and Kelvin's theorem / Ch. 6 lift**: ch02 E5 (curl heat image, paddle wheel, unrolled u·t loop,
  "Stokes failing on purpose" for a line vortex).
- **Ch. 5 vorticity diffusion / Ch. 8 Stokes' first and second problems**: E2's gap + profile + wall-stress stage with
  a transport and a steady (or periodic) ghost; E2's modes show "same PDE, different D".
- **Ch. 7 internal waves**: E4's column + profile + time series on one clock, ω ≤ N marker, atmosphere/ocean modes and
  the two-convention badge where lapse rates appear.
- **Ch. 11 instability (Richardson number)**: E4's badge pattern for Ri < ¼ with the exact criterion text from fluidpy.
- **Ch. 13 GFD**: E4 for stratification (θ view, inversion preset), E5 for Ro/Ri/Rossby-radius scaling.
- **Ch. 15 compressible**: E3's linked state diagrams with term bars for Fanno/Rayleigh lines and shock entropy rise.

## Reference explainers (the depth to match)
See skill `interactive-viz` §4–§5: Shammunul's preferred MIT-mathlet re-implementations (forced damped vibrations — the
explanation panel; amplitude and phase, second order II — a numbered live derivation; angular frequency explorer —
linked views and modes) and the fast.ai labs (FID formula lab, forward noising lab, pixels as parameters, random copy,
overfitting curves, stride & padding playground, 3-D U-Net and ResNet). Passing in-repo templates that use every engine
feature: `templates/viz_example.html`, `templates/viz_example_field.html`, `templates/viz_example_3d.html`. In-chapter
models after ch01: E4 `parcel_stability` (linked views, conventions, 4 derivations) and E3 `heat_work_paths` (term bars,
modes, linked diagrams). After ch02: `cauchy_traction_principal_axes` (the reference-depth Explain panel of the
chapter: 9 sections, a ★★★ derivation with live numbers) and `strain_vs_rotation_split` (the cleanest one-clock story).
