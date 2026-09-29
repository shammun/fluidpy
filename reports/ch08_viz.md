# Chapter 8 — explainer review                                  2026-09-29

Fresh-eyes gate for phase 7. `tools/viz_lint.py --chapter ch08` was re-run (9 × `ok`), `tools/shot.py` was re-run on each file one at a time
(8 sizes, every tab, every walkthrough step, every derivation page), and `tools/eq_refs.py --chapter ch08` was re-run (0 bare citations).
I read the screenshots in `reports/viz/ch08/<slug>/` (desktop and phone-tall walkthrough steps, phone and phone-land Explore, notebook Equations,
laptop Check, desktop and phone-tall Explain, desktop and phone-tall derivation pages). I compared every explainer's `cfg.derivations` step
titles with `notebooks/build_ch08.py --dump`, which is now the reference, using a Playwright script. I also read the JS physics against
`fluidpy`, including `fluidpy/core/lubrication.py::slider_bearing_state`.

| slug | lint | shot (sizes/views) | parity rows ok | explain | derivations (steps ok) | teaching | polish | verdict |
|---|---|---|---|---|---|---|---|---|
| couette_poiseuille_backflow (E1) | ok | PASS 8 / 376 | 28/28 (20 py) | 6 § + interp | D02 7/7 · D03 5/5 · D04 9/9 · **D05 7/9** | good | good | **FAIL** |
| lubrication_scaling (E2) | ok | PASS 8 / 288 | 19/19 (14 py) | 7 § + interp | D10 14/14 · D11 7/7 | very good | good | PASS |
| slider_bearing (E3) | ok | PASS 8 / 432 | 39/39 (31 py) but **missing rows** | 10 § + interp, **§7 wrong** | D12 8/8 · D13 12/12 · **D14 15/16** | very good | good | **FAIL** |
| viscous_gravity_current (E4) | ok | **FAIL** 8 / 272 (full-hd equations clip) | 24/24 (15 py) | 7 § + interp | D15 10/10 · D23 9/9 | good | raw TeX in symbols | **FAIL** |
| stokes_first_problem (E5) | ok | PASS 8 / 488 | 32/32 (23 py) | 8 § + interp | D16 6/6 · D17 8/8 · D18 9/9 · **D19 11/12** · D20 5/5 | very good | ghost labels | **FAIL** |
| similarity_exponents (E6) | ok | PASS 8 / 384 | 28/28 (19 py) | 7 § + interp | D21 8/8 · D22 14/14 · D23 9/9 | very good | good | PASS |
| oscillating_plate (E7) | ok | PASS 8 / 256 | 29/29 (23 py) | 8 § + interp | D24 10/10 · D25 6/6 | excellent | excellent | PASS |
| stokes_sphere_flow (E8) | ok | PASS 8 / 456 | 26/26 (20 py) | 6 § + interp | D28 12/12 · **D29 10/11** · D32 7/7 · D33 7/7 | very good | good | **FAIL** |
| stokes_drag_settling (E9) | ok | PASS 8 / 320 | 44/44 (33 py) | 10 § + interp | D30 11/11 · D31 13/13 | excellent | good | PASS |

The contract holds for all nine: live Explain tab, Derivation tab with the assigned D rows, synced Code tab, ≥ 5 depth features each, 6–8
walkthrough steps, 4–5 quiz questions, and every equation written out next to its number (`eq_refs` = 0). The only failed fit check was E4's
full-HD Equations pager. Parity tolerances are honest: 1e-12 for closed forms; 1e-10…1e-11 for series, quadrature and Crank–Nicolson; 1e-6 only
for the golden-section vs Brent taper optimum; 2e-5 only for E8's inertia/viscous ratio, where fluidpy uses finite differences.

Checked against the verified physics, and the explainers agree on each point:

- Couette–Poiseuille (E1): extremum y* = h/2 − μU/(h·dp/dx), JS `state()` plus a parity row at dp/dx = −4 → 7.5 mm. Backflow sets in at
  dp/dx > 2μU/h², and Q = 0 at 6μU/h².
- Engine film (E2): ρ = 870, μ = 0.05, ε²Re_L = 4.35e-3.
- Stokes' first problem (E5): ∫ω dy = +U (invariant row). 2 erfc⁻¹(0.05) = 2.772 (E5 row; E6 text "2.772, the page rounds to 2.76").
- Sphere (E8): the axis ratio → ½ Re_a r/a (invariant row). Crossover r_c = 41.5a at Re = 0.1, which is ≈ 2a/Re_a.
- Settling (E9): rear pressure −3μU/2a (parity row at θ = 0). Oseen C_D = (24/Re)(1 + 3Re/16) with Re on the diameter.
- Glacier preset (E4): μ = 10¹³ Pa s.

## couette_poiseuille_backflow

**Must fix**
1. **D05 is 2 steps short of the notebook (7 vs 9).** Append the notebook's two moves after "Solve slope < 0 for dp/dx" (desktop__derive-d4p7/8
   end at step 7):
   - step 8 "Set the slope to zero for the extremum": `\frac{du}{dy}=\frac Uh-\frac1{2\mu}\frac{dp}{dx}(h-2y^*)=0`.
     *Why:* step 3 divided by μ gives the slope at any height; the extremum is where it vanishes.
   - step 9 "Solve for the extremum's height": `y^*=\frac h2-\frac{\mu U}{h\,dp/dx}`.
     *Why:* multiply by 2μ/(dp/dx) and rearrange; it is interior only when 0 < y* < h (mid-gap for U = 0).
     Give it a `live` with the reader's numbers, e.g. the favourable preset gives 7.5 mm.
   - Add the notebook's result line "extremum at y* = h/2 − μU/(h dp/dx) (0 < y* < h)" to D05's result.
   - The walkthrough link `derive: {D05, step 7}` stays valid.

**Should fix**
- The explainer computes y* and u_max (`state()` has `y_umax`, `u_max` with parity rows) but never shows them. Add a u_max dot on the u(y) curve
  and one Explain line: y* = h/2 − μU/(h dp/dx) = …, u_max = …. That way D05 step 9 has a picture to point at.
- notebook__equations: in the narrow τ/Q view the "Q 3.33" label sits on top of the bars. Drop the value labels below ~200 px width.

**What works**
- The channel shows an amber reversed layer with tracers and a bending dye line.
- The profile view shows the rose Couette and orange Poiseuille parts with the blue sum, plus a black floor tangent that stands vertical at
  onset (phone-tall__tour-step4 shows it at 2 Pa/m).
- The inspector traces u = 0.025 − 0.0375 = −0.0125 m/s. The Q waterfall bars (Qc + Qp = Q) make "backflow with Q > 0" obvious.

## lubrication_scaling

**Must fix** — none.

**Should fix**
- phone-land__explore: the x-bars keep only their coefficients ("1", "ε²") and lose the term names. Keep a short name ("inert.", "press.",
  "fric.") or put the colour legend in the view title.
- D10 step 8's "we had" line on phones is only the multiplier "× ρε²L²/(μU)" (phone-tall__derive-d1p9). Give step 8 a full equation line (the
  x-equation after the multiplication) so "we had → now" reads as a move.

**What works**
- The log term bars with a grey "dropped" band. The p-scale chips show that Λ only rescales p*.
- On phones the y-momentum verdict "∂p*/∂y* = O(ε²) = 1×10⁻⁶" moves into the title.
- The "thick gap ✗" preset breaks the theory visibly.

## slider_bearing

**Must fix**
1. **The recirculation flag mirrors the old fluidpy.** JS `sliderState()` tests only x = L (`hin = h0*(1 + a)`). fluidpy's
   `slider_bearing_state` now reports recirculation at the **wide end**: `inlet_backflow == backflow_any`, plus `backflow_x`, true for α > 1
   **or α < −½**.
   - At α = −0.6 fluidpy returns `inlet_backflow: True, backflow_x: 0.0`; the JS returns False. The α slider reaches −0.9, so readers can hit
     this case.
   - Evaluate the wide end (x = L for α > 0, x = 0 for α < 0) and return `backflow_x` and `backflow_any` too.
   - Add parity rows: `ch08.slider_bearing_state(5e-5, -0.6, 0.05, 5.0, mu=0.05)["inlet_backflow"]` (1),
     `(…, -0.49, …)` (0), `["backflow_x"]` at α = 1.2 (0.05) and at α = −0.6 (0.0), and `["backflow_any"]`.
   - The local `recirc()` already handles both ends correctly.
2. **Explain §7 text is wrong** (`slider_bearing.html` ~l.2823).
   - "The flag of `slider_bearing_state` tests x = L only" must go. Say the flag reports recirculation at the wide end, wherever it is.
   - The "Off" branch says it changes sign "only when 3α/(2+α) > 1, i.e. α > 1". It must say "gap ratio > 2: α > 1, or α < −½ for a pad
     leading with its narrow end".
   - Rename every "inlet backflow" to **"recirculation at the wide end"**. That covers the Code tab comment `print(st["inlet_backflow"])  #`,
     the physics header comment "inlet backflow flag (pad frame, x = L)" and the parity-row names.
   - Quiz item 5 ("Above which taper…") should mention α < −½ as well.
3. **D14 has 15 steps; the notebook has 16.**
   - Insert the notebook's step 12 "Substitute C₁ and C₂ into step 7":
     `p-p_e=\frac{6\mu LU}{\alpha h_o^2}\Big[-\frac{1+\alpha}{(2+\alpha)(1+s)^2}+\frac{1}{1+s}-\frac{1}{2+\alpha}\Big],\ s=\frac{\alpha x}{L}`.
     *Why:* insert C₁ and C₂ and subtract p_e; 6μLU/(αh₀²) comes out of all three terms. Today step 12's "we had" is C₂ and "now" is the
     finished fraction (desktop__derive-d3p12), which skips this move.
   - Old step 12 → 13 "Put p − p_e over one denominator". Its *why* becomes the expanded numerator:
     −1 − α + 2 + 2s + α + αs − 1 − 2s − s² = s(α − s).
   - Old 13 → 14 "Write s back in terms of x".
   - Old 14 → 15 "Keep the linear term in α". Its *why* must change to the notebook's: "An extra approximation, not required by lubrication
     (the wall slope αh₀/L is tiny for any α of order one): for α ≪ 1, (2 + α)(1 + αx/L)² ≈ 2; it costs accuracy — the linear load is 15.6 %
     high at α = 0.1". The current *why*, "the lubrication profile of D12 is only valid for small slopes anyway", is the retracted reasoning.
   - Old 15 → 16 "Integrate for the load".
   - **Renumber the walkthrough link**: step 6 "The printed formula" uses `derive: {id: 'D14', step: 13}`, which must become 14. Shift any
     per-step `set`/`live`/`highlight` that is keyed by index (the `compare` model set and the W live line).

**Should fix**
- phone-tall__derive-d3p14: page 1 of step 14 shows only "we had" and the move title above ~450 px of blank space; "now" is on page 2. Shorten
  the "we had" block, for example by showing only the right-hand side on phones.
- With U < 0 (desktop__tour-step7) the W(α) view labels the minimum "optimum 1 + α = 2.189". Call it "strongest suction" when αU < 0.

**What works**
- The pad-frame picture with u − U profiles, tracers and orange pressure arrows, and the black dp/dx tangent on the pressure curve.
- The printed-formula ghost with a measurable slope miss (invariant row 1.398).
- A golden-section optimum that matches fluidpy's Brent, and a 10-section Explain with the flux split C₁ = pressure part + drag part.

## viscous_gravity_current

**Must fix**
1. **shot FAIL:**
   > `[full-hd] equations: pager page 1/1: viz-eq 'Front constantHuppert (1982)ηN=[15(310)1' clipped 19px (cut through a line)`

   `CLIP__full-hd__equations__pg1.png` shows the η_N card cut at the bottom. The pager does not split the page at 1920×1080, apparently because
   its height estimate misses the live line. Possible fixes: move the η_N Gamma-function card to its own page or into its item's `note`; drop
   the live "η_N = 1.41124" row (η_N already appears in the symbols line); or split the Huppert card into two items.
2. **Raw TeX shown to the reader.** In the Equations symbols list, `['\\beta', '\\rho g/3\\mu', '1/(m s)']` (`viscous_gravity_current.html`
   l.2991) renders as the literal text "\rho g/3\mu" (visible in the CLIP screenshot). Write `'ρg/(3μ): thin-film coefficient'` or wrap it in
   `$…$`.

**Should fix**
- The view titles and axis labels use caret text: "h τ^m vs x/τ^m", "h τ^⅕ / x₀", "x / (x₀ τ^⅕)", and the page title "t^(1/5)". Use
  superscript glyphs as the other explainers do.
- phone__explore: the layer view is only ~90 px of plot at 360×640, and the dome with its arrows is hard to read. Consider `hidePortrait` on the
  front view during Explore, or giving layer the larger row weight on portrait.

**What works**
- A conservative finite-volume march with volume conserved to 10⁻¹⁵ (shown in the status).
- The front on log–log bending onto slope ⅕ next to diffusion's ½ ghost. Two humps forget their start.
- The march is proven against fluidpy's implicit solver (three parity rows) and against Huppert's η_N from the Gamma function.

## stokes_first_problem

**Must fix**
1. **D19 has 11 steps; the notebook has 12.**
   - Insert the notebook's step 10 "Insert A and B into (8.29)": `F(\eta)=1-\frac{1}{\sqrt\pi}\int_0^\eta e^{-\xi^2/4}\,d\xi`.
     *Why:* put B = 1 (step 6) and A = −1/√π (step 9) back into $F(\eta)=A\int_0^\eta e^{-\xi^2/4}d\xi+B$ (8.29).
   - Old 10 → 11 "Substitute again inside F". Its *why* becomes the notebook's: "Substitute ξ = 2ζ again (dξ = 2dζ, ξ²/4 = ζ²): the factor 2
     turns 1/√π into 2/√π, and the upper limit ξ = η becomes ζ = η/2". The current "The same ξ = 2ζ; the upper limit becomes η/2" hides where
     the 2 comes from.
   - Old 11 → 12 "Recognise the error function". Move its `live` ("u/U at your probe") with it.
   - No walkthrough step links to D19, so there is nothing else to renumber.

**Should fix**
- δ(t) view (desktop__tour-step5, crop checked): the ghost labels sit at fixed offsets, so "water" sits just above the **air** line and "air"
  sits on the 1 m gridline between the air and honey lines. The water ghost is hidden under the orange current line. Anchor each label at its
  own line's end, or put "ghosts: air (upper) · honey (top)" in the title.
- On phones the rescaled view has no label naming the grey band as 1 − erf(η/2) (the desktop has "the one curve 1 − erf(η/2)"). Add it to the
  view title on portrait.

**What works**
- A log clock from 0.1 s to 1000 s, a raw/rescaled mode switch, and Crank–Nicolson dots on the erfc curve.
- The "stop at T" and "second wall" presets break similarity visibly (spread 1). The δ₉₉ coefficient 3.643 is computed rather than typed.
- The ∫ω dy = +U correction appears in an Equations note.

## similarity_exponents

**Must fix** — none.

**Should fix**
- Step titles replace the notebook's "(8.20)"/"ansatz" with words ("the diffusion equation", "the guess"). That is fine under the house rule;
  keep it consistent with the notebook's glossary ("ansatz" is used in the notebook's D22 steps 2 and 8).

**What works**
- The exponent plane with the 2m = 1 and n = m lines and a star at the answer. Bracket-power chips turn green when the powers match.
- The conserved-integral inset "∫ω dy ∝ t^0.1" shows *why* n is wrong even when the brackets match (desktop__tour-step4).
- Four systems as modes on one control set. Every displayed spread equals fluidpy's `similarity_collapse_error` (display = parity rows).

## oscillating_plate

**Must fix** — none.

**Should fix**
- Phone view axes read "y/δ_e" with a literal underscore. Use "y/δₑ".

**What works** (candidate for `knowledge/viz_patterns.md`)
- D24 step 7 draws the discarded growing root as a rose dotted curve that explodes (desktop__derive-d1p7): the reason for B = 0 is *seen*.
- The "start from rest" Crank–Nicolson march lands on the periodic solution within a period, answering "where is the initial condition?".
- The book depth 4√(ν/ω) with its 5.9 % is drawn on both views, and the lag arrow appears on u(t). Every preset is a real Stokes layer (tide,
  molecular vs eddy ν, 1 kHz air, Ekman look-alike).

## stokes_sphere_flow

**Must fix**
1. **D29 has 10 steps; the notebook has 11.** Step 4 "Apply E² again" (desktop__derive-d2p4) does two moves in one: it states g″ − 2g/r² = 0
   and expands in f. Split it as the notebook does:
   - 4 "Apply E² again, to g": `E^2(g\sin^2\theta)=\Big(g''-\frac{2g}{r^2}\Big)\sin^2\theta=0,\ g=f''-\frac{2f}{r^2}`.
     *Why:* step 3 showed E²(f sin²θ) = g sin²θ; the same rule applied to g sin²θ gives (8.44) as an ODE for g.
   - 5 "Expand the ODE in f": `f^{iv}-\frac{4f''}{r^2}+\frac{8f'}{r^3}-\frac{8f}{r^4}=0`.
     *Why:* insert g = f″ − 2f/r² and use (f/r²)″ = f″/r² − 4f′/r³ + 6f/r⁴. I checked the expansion by hand:
     g″ − 2g/r² = f⁗ − 4f″/r² + 8f′/r³ − 8f/r⁴ ✓.
   - Steps 5–10 → 6–11. Carry the per-step `set`/`watch`: "A = 0 … zoom 100" goes to step 9, and the live "C = −3Ua/4, D = Ua³/4" to step 11.
   - Rename 9 → "Apply the far field" and 10 → "Apply the two surface conditions", and give each its equation in TeX
     (ψ → ½Ur² sin²θ (8.47); ψ(a) = 0 (8.45), ∂ψ/∂r(a) = 0 (8.46)), so the steps stay aligned with the notebook's titles.
   - No walkthrough step links to D29.

**Should fix**
- phone-tall__tour-step4 ("Seen from the still fluid"): the captured frame has the sphere at the right edge of the view with only a sliver
  visible, and portrait hides the title and axis labels. Start the fluid-frame transport so the sphere is in view, or centre the window on it.

**What works**
- Three models on one control set: Stokes, Oseen and ideal.
- A dashed 1 % disturbance contour that shows the Oseen wake. The amber inertia = friction circle is computed by brentq on the exact axis ratio
  (41.5a at Re = 0.1, ≈ 2a/Re_a) with the ½Re_a r/a asymptote drawn.
- The side-line speed at 10a: 0.9248U for Stokes vs 1.0005U for ideal.

## stokes_drag_settling

**Must fix** — none.

**Should fix**
- phone-tall__tour-step2: the portrait "Stresses" view shows two curves identified by colour alone. Add "p − p∞ orange · σ_rθ rose" to its
  title on portrait.

**What works** (candidate for `knowledge/viz_patterns.md`)
- The traction sweep: orange pressure plus rose friction arrows whose sum is the same black x-push everywhere.
- Running-drag curves that reach 2πμaU and 4πμaU, and an end card "2πμaU (⅓) + 4πμaU (⅔) = 6πμaU".
- The C_D·Re/24 mode with Stokes, Oseen and Morrison curves at Re = 15.9 (×3.99 and ×1.96).
- A real-particle table with the current row highlighted, and a ring-sum vs Gauss–Legendre parity row.

## Promote to `assets/viz_lib.js` (recurring across the ch08 explainers — for the knowledge-keeper)

| helper | where | proposal |
|---|---|---|
| `erfcHP(x)`, `erfcinvHP(p)` (double-precision erfc: continued fraction above 3, series below; Newton inverse) | E5, E6 (+ ch04 E4) | `Viz.num.erfcHP / erfcinv` — `Viz.num.erfc` (1.2e-7) is too coarse for parity rows |
| `arrowPx` / `pxArrow` (pixel arrow with head) | E1 + 6 others | `Viz.arrow(ctx, x0, y0, x1, y1, {color, width, head})` |
| `hatch` / `hatchShift` (moving hatched wall) | E1, E3, E5, E7 | `Viz.hatch(ctx, x, y, w, h, {step, shift, color})` |
| `gutterPlot` (plot with a reserved y-title strip, the ch01 narrow-view lesson) | E3, E5, E7 | fold into `v.plot({ytitleStrip: true})` |
| `pn(v, k)` (Python-literal number for Code-tab `live`) | E3, E5, E7, E8, E9 (three variants) | `Viz.pyNum(v, sig)` — one rule for exponents |
| `phone()` / `lay()` (layout probe) | all 9 | `app.layout` / `Viz.isPortrait()` |
| `f3`/`t3`/`t4` short formatters | all 9 | already `Viz.fmt`/`Viz.tnum` — document the `sig` shorthand, or add `Viz.f3` |
| `fmtLen`/`fmtLenM` (µm/mm/cm/m unit switch) | E4, E6 (E5 inline) | `Viz.fmtLen(m)` beside `Viz.fmtTime` |
| `lab` (`Viz.text` with `bg: true, size: 12`) and `logPlot` (log–log axes with decade ticks) | E3, E8 and others | a `logx/logy` option on `v.plot` |

### Round-1 verdict: FAIL

Four explainers PASS: E2, E6, E7, E9. Five FAIL:

| explainer | reason |
|---|---|
| E1 couette_poiseuille_backflow | D05 is missing the notebook's 2 extremum steps |
| E3 slider_bearing | its recirculation flag disagrees with fluidpy for α < −½ (no parity row catches it); Explain §7 and the code say "x = L only / inlet backflow"; D14 is missing the inserted step 12, keeps the retracted step-15 reason, and its walkthrough link must become step 14 |
| E4 viscous_gravity_current | shot FAIL (full-HD Equations clip) and raw TeX in the symbols list |
| E5 stokes_first_problem | D19 is missing step 10 |
| E8 stokes_sphere_flow | D29 step 4 must be split |

Every fix is local. The physics, parity rows, depth and teaching are at the reference level throughout.

---

# Round 2 — re-review of the five fixed explainers                          2026-09-29

What I re-ran:
- `viz_lint --chapter ch08`: 9 × ok.
- `eq_refs --chapter ch08`: 0 bare citations.
- `shot.py` on each fixed file, one at a time: E1, E3, E4 and E8 PASS. E5 FAILED once (below); it passed a second full run and a phone-land-only run.
- Step titles compared again against `notebooks/build_ch08.py --dump` using the Playwright dump of `cfg.derivations`. **Every D row now has the notebook's step count.** The remaining title differences are all words replacing a bare "(8.xx)", for example "Use the cross-gap balance" for "Use (8.17b)" and "Insert A and B into the first integral" for "Insert A and B into (8.29)". Where such a step uses the equation, its *why* writes the equation out, which is acceptable under the house rule.

| slug | shot (round 2) | parity rows | round-1 must-fix | verdict |
|---|---|---|---|---|
| E1 couette_poiseuille_backflow | PASS 8 / 392 | 30/30 (22 py) | D05 9/9 ✓ | **PASS** |
| E3 slider_bearing | PASS 8 / 440 | 46/46 (38 py) | wide-end flag + 7 rows ✓ · §7/quiz/code wording ✓ · D14 16/16 ✓ · tour → D14 step 14 ✓ | **PASS** |
| E4 viscous_gravity_current | PASS 8 / 272 | 24/24 (15 py) | full-HD clip gone ✓ · symbols plain text ✓ | **PASS** |
| E5 stokes_first_problem | **FAIL**, then PASS on re-run | 32/32 (23 py) | D19 12/12 ✓ | **FAIL** (intermittent clip) |
| E8 stokes_sphere_flow | PASS 8 / 464 | 26/26 (20 py) | D29 11/11 ✓, per-step settings moved ✓ | **PASS** |

## Round-2 findings

### E1 couette_poiseuille_backflow — PASS
- D05 steps 8 "Set the slope to zero for the extremum" and 9 "Solve for the extremum's height" match the notebook.
- Step 9 has live numbers: y* = 5 mm + 2.5 mm = 7.5 mm, u(y*) = 11.3 cm/s.
- The should-fix is done: the profile now marks "max 11.3 cm/s, y* 7.5 mm" (desktop__derive-d4p9).
- The walkthrough link to D05 step 7 is still valid.

### E3 slider_bearing — PASS
- `sliderState()` evaluates the wide end (x = L for α > 0, x = 0 for α < 0) and returns `backflow_x` and `backflow_any`.
- New parity rows cover α = −0.6 (1), α = −0.49 (0), `backflow_x` at 1.2 (0.05) and at −0.6 (0), and `backflow_any` at −0.6 and 0.5. All pass.
- Explain §7:
  - The "On" branch now says the flag reports recirculation at the wide end and gives x.
  - The "Off" branch says "gap ratio exceeds 2: α > 1, or α < −½ for a pad whose gap narrows along x".
- Quiz 5 names both α > 1 and α < −½ ("Try α = 1.5, then α = −0.6").
- The physics header now reads "recirculation at the wide end". The only "inlet" left is fluidpy's API key `st["inlet_backflow"]`, which is correct.
- D14 now has 16 steps:
  - Step 12 "Substitute C₁ and C₂ into step 7" (desktop__derive-d3p12). I checked the three-fraction line numerically at α = 0.5, x/L = 0.4: 0.03333, equal to the exact pressure.
  - Step 13's *why* expands the numerator to s(α − s).
  - Step 15's *why* is the notebook's "extra approximation … 15.6 % high at α = 0.1".
- Walkthrough step 6 links to D14 step 14. The phone-tall blank space on D14 step 14 is the known library-pager limit, so it is not counted.

### E4 viscous_gravity_current — PASS
- full-hd__equations now pages (1/2) and nothing is clipped.
- The β symbol reads "thin-film coefficient ρg/(3μ)".
- Should-fixes done:
  - The titles use superscripts (t^(1/5) → t¹ᐟ⁵).
  - phone__explore gives the layer view most of the stage: a readable dome with its arrows at 360×640.

### E5 stokes_first_problem — FAIL (one remaining must-fix)
- D19 now has 12 steps:
  - Step 10 "Insert A and B into the first integral" writes $F(\eta)=A\int_0^\eta e^{-\xi^2/4}d\xi+B$ (8.29) in its *why*.
  - Step 11's *why* explains where 2/√π comes from.
  - Step 12 carries the erf definition and the live value.
- Should-fix done: the δ(t) ghost labels now sit on their own lines (honey / air, with "water (yours)" on the orange line).

**Must fix**
1. The first round-2 full audit failed at phone-land:
   > `[phone-land] equations: pager page 6/6: viz-eq 'Vorticity content§8.4∫0∞ωz dy=∫0∞(−∂u/∂y' clipped 18px (cut through a line)`

   `CLIP__phone-land__equations__pg6.png` shows the note under the "Vorticity content" card cut off at the bottom of the panel. A phone-land-only re-run and a second full run both passed, so the pager's height estimate is racing the render of the tall live integral. It is intermittent, but readers can hit it.

   Make the card smaller so it fits with margin:
   - render the live row inline ($\int_0^\infty\omega\,dy = 0.1$ m/s in `\textstyle`, or plain text) instead of a display-size ∫;
   - or put the one-line note ("the page prints −U …") in *why* or on its own item.

   Then run shot twice to confirm.

### E8 stokes_sphere_flow — PASS
- D29 step 4 "Apply E² again, to g" and step 5 "Expand the ODE in f" are now separate. Step 5's *why* shows the whole expansion, which I checked by hand.
- The far-field step (9) carries the r⁴ reason and (8.47) in TeX. The surface-conditions step (10) writes out (8.45) and (8.46). Step 11 solves for C and D.
- phone-tall__derive-d2p5 shows only "we had" and the move on page 1, with "now" on page 2. This is the same known pager limit, not counted.

## Per-explainer verdicts after round 2

| # | slug | verdict |
|---|---|---|
| E1 | couette_poiseuille_backflow | PASS |
| E2 | lubrication_scaling | PASS (round 1) |
| E3 | slider_bearing | PASS |
| E4 | viscous_gravity_current | PASS |
| E5 | stokes_first_problem | FAIL — intermittent phone-land Equations clip (1 must-fix) |
| E6 | similarity_exponents | PASS (round 1) |
| E7 | oscillating_plate | PASS (round 1) |
| E8 | stokes_sphere_flow | PASS |
| E9 | stokes_drag_settling | PASS (round 1) |

### Round-2 verdict: FAIL
8 of 9 PASS. E5 needs its "Vorticity content" equation card shrunk until `shot.py` passes consistently. It is a one-item fix; the physics, derivations and teaching are all at PASS level.

---

# Round 3 — E5 stokes_first_problem (phone-land clip)                        2026-09-29

What I checked:
- `viz_lint`: ok.
- `shot.py viz/ch08/stokes_first_problem.html`, run twice in a row, full 8 sizes each time:
  - run 1: **PASS** (8 sizes, 496 views, 7 steps, 32 selftest rows, 0 failures);
  - run 2: **PASS** (same numbers, 0 failures).
- Selftest: 32/32 ok (23 fluidpy parity rows). The invariant "integral of omega = +U" still passes (0.99999999999999).
- All 5 Equations pages at phone-land (844×390), captured with my own Playwright pass:
  - Page 4: the (8.29) card $F(\eta)=A\int_0^\eta\exp(-\xi^2/4)\,d\xi+B$ now has a text-size ∫ and fits.
  - Page 5: the "Vorticity content" card $\int_0^\infty\omega_z\,dy=\int_0^\infty(-\partial u/\partial y)\,dy=+U$ fits with room below its note.
- The physics is still right. The live row reads $u(0)-u(\infty)=0.1$ m/s at U = 0.1 m/s, which is +U. It is computed as `uS(0) − uS(far)`. The note says "The page prints −U; u falls from U to 0, so it is +U at every time".

E5 verdict: **PASS**. The round-2 must-fix is closed.

## Final per-explainer verdicts

| # | slug | verdict |
|---|---|---|
| E1 | couette_poiseuille_backflow | PASS (round 2) |
| E2 | lubrication_scaling | PASS (round 1) |
| E3 | slider_bearing | PASS (round 2) |
| E4 | viscous_gravity_current | PASS (round 2) |
| E5 | stokes_first_problem | PASS (round 3) |
| E6 | similarity_exponents | PASS (round 1) |
| E7 | oscillating_plate | PASS (round 1) |
| E8 | stokes_sphere_flow | PASS (round 2) |
| E9 | stokes_drag_settling | PASS (round 1) |

Known, not counted: the library pager leaves blank space on phone-tall derivation pages when one "now" box is placed as a unit (E3 D14 step 14, E8 D29 step 5). This is a candidate for a `viz_lib.js` pager improvement by the knowledge-keeper.

## Verdict: PASS
