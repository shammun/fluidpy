# Explainer patterns — what worked, what failed, what to promote

Read by the curator, designer, viz-builders and viz-reviewer before any explainer work. Appended by the knowledge-keeper.
Explainer short names: ch01 E1 `continuum_averaging_volume`, E2 `viscosity_momentum_diffusion`, E3 `heat_work_paths`,
E4 `parcel_stability`, E5 `buckingham_pi_machine`.

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

## Ideas for later chapters
(seeded from `book.yaml → viz_seeds`; the curator decides)
- **Ch. 2 tensors**: reuse E5's "matrix view + highlighted block + live determinant" layout for rotation matrices and
  principal stresses; the Cauchy-traction idea is E4's arrow-on-a-plane annotation in 2-D.
- **Ch. 3 kinematics**: streamline/pathline/streakline stage = template field stage + E2's "ghost + trail" idea.
- **Ch. 4 conservation laws**: dynamic similarity explainer = E5 with presets for Re, Fr, Ro and a model-vs-prototype
  table (promote `Viz.Frac` then); Boussinesq buoyancy reuses E4's column + profile and `brunt_vaisala_sq`.
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
modes, linked diagrams).
