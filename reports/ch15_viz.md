# Chapter 15 — explainer review                                  2026-10-14

**Verdict: FAIL.** 8 of 9 explainers pass. `nozzle_back_pressure` has 2 Must-fix items. Both are small text or config edits, with no physics change. Every mechanical check passes.

| Check (run in this session) | Result |
|---|---|
| `tools/viz_lint.py --chapter ch15` | 9 / 9 ok |
| `tools/shot.py --chapter ch15` (full: 8 sizes, every tab, step, derivation page and pager page) | **9 / 9 PASS, 2 360 views, 0 failures, 0 CLIP screenshots, 0 console messages, KaTeX loaded at every size**, exit 0 |
| selftest | **295 rows, 295 ok, 245 of them JS ↔ fluidpy**. Closed-form parity rows use rtol 1e-12 (root finders 1e-9 … 1e-7). Largest error/tolerance: 0.67 (`plateau flux at p_B 0.3`, rtol 1e-6, a 6-digit constant), 0.36 (NACA 1135 row, rtol 1e-4, a table value), 0.34 (`diamond thickness drag`, rtol 1e-3, see E9 S2) |
| `tools/eq_refs.py --chapter ch15` | 3 hits, all `ref:` labels on Code-tab panels of `normal_shock_jump` (`Eqs. (15.39)–(15.44a)`, `Eqs. (15.35)–(15.37)` ×2). These are panel labels; the equations are written out in the Equations tab. Accepted, as in ch14. No bare number in tour, Explain, Derivation or quiz text |
| `tools/check_public.py viz/ch15/*.html` | OK (9 files) |
| Derivation tabs against `notebooks/build_ch15.py --dump` | 18 / 18 D ids present, with step counts equal to the notebook's (D03 8, D04 5, D05 5, D10 9, D12 7, D13 13, D17 11, D19 7, D20 9, D21 13, D23 8, D25 4, D26 8, D27 8, D28 8, D29 9, D30 10, D31 7). Every walkthrough `derive: {id, step}` link points at an existing step |
| Equation fidelity (rendered pages p778, p789, p795, p800: (15.42), (15.46), (15.50)–(15.53), (15.56), (15.57)) | match. The corrected signs of (15.32), (15.47) and (15.48) are labelled "corrected" wherever they are shown |
| Numbers quoted in the text, re-computed with fluidpy | all correct. Examples: M₁ for a 10 % p₀ loss is 1.587 (air) and 1.627 (γ = 5/3). Half of p₀ is lost at M₁ = 2.498. Laval nozzle at p_B 0.7: M₁ 1.863, A/A_t 1.510. Paradox band T₂/T₁ 0.9947. Linear-rule error 2.2 / 10.5 / 20.1 % at 1° / 5° / 10°. Flat plate at 10°: lift −1.1 %, suction side +24.8 %, compression side −20.1 %. ν(2) = 26.38° |
| Walkthrough step length | every step ≤ 40 words (limit 45) |

| slug | lint | shot (sizes/views) | parity rows ok | explain | derivations (steps ok) | teaching | polish | verdict |
|---|---|---|---|---|---|---|---|---|
| sound_pulse_dalembert | ok | PASS 8 / 288 | 24/24 (21 py) | 8 sections + reading, 8 live | D03 8/8 · D04 5/5 · D05 5/5 | very good | good | **PASS** |
| area_mach_nozzle | ok | PASS 8 / 256 | 38/38 (31 py) | 9 sections + reading | D10 9/9 · D12 7/7 | very good | good | **PASS** |
| normal_shock_jump | ok | PASS 8 / 328 | 31/31 (25 py) | 10 sections + reading | D13 13/13 · D17 11/11 | very good | very good | **PASS** |
| nozzle_back_pressure | ok | PASS 8 / 264 | 35/35 (30 py) | 8 sections + reading | D19 7/7 · D20 8/9 (step 4 *why* has a gap) | very good | preset chips wrong | **FAIL** |
| duct_friction_heating | ok | PASS 8 / 304 | 41/41 (34 py) | 9 sections + reading | D21 13/13 · D23 8/8 | very good | very good | **PASS** |
| mach_cone | ok | PASS 8 / 136 | 23/23 (19 py) | 7 sections + reading | D25 4/4 | very good | good (phone overlap) | **PASS** |
| oblique_shock_wedge | ok | PASS 8 / 256 | 34/34 (28 py) | 8 sections + reading | D26 8/8 · D27 8/8 | very good | good | **PASS** |
| prandtl_meyer_fan | ok | PASS 8 / 280 | 36/36 (30 py) | 8 sections + reading | D29 9/9 · D30 10/10 | very good | good | **PASS** |
| supersonic_flat_plate | ok | PASS 8 / 248 | 33/33 (27 py) | 9 sections + reading | D28 8/8 · D31 7/7 | very good | good | **PASS** |

Depth compared with the references (`forced_damped_vibrations.html`, and AFE / APS as each storyboard names). Every explainer has a numbered "Explanation & interpretation" panel with:

- a "what the views show" section, with the colours named;
- every displayed quantity worked out with the reader's values (formula → numbers → boxed result);
- a hint for the view hidden on phones;
- a regime-dependent "Reading the current setting" section.

Each also uses 6–8 depth features (linked views, presets, status, plus terms, inspector, modes, notes or transport). None is shallower than its reference.

---

## sound_pulse_dalembert (E1)
**Must fix:** none.

**Should fix**
1. `phone__explore.png`: the label "stream U = 102 m/s →" sits on a solid box over the f/g curves and the probe dot. When `v.w < 420`, move it into the view title or to the top-right corner.
2. `phone-tall__tour-step4.png`, step 4 "The air only shuffles": the u′ arrows in the air strip are barely visible among the dots on phones, yet this step is about them. Draw them thicker (or longer) on small views.

**What works:** the erf closed form for the parcel shift, checked against a trapezoid sum of u′ along the path. "shifts drawn ×500 (true ≤ 0.886 mm)" states the exaggeration honestly. The x–t bands make U − c change sign visibly at M = 1. The Newton/Laplace medium table.

## area_mach_nozzle (E2)
**Must fix:** none.

**Should fix**
1. `phone__explore.png`: at x/L = 1 the cursor label "M 2.20" runs into the round exit handle. On narrow views, right-align the label so it ends before `X1 − 12`.
2. Explain §7 hint (`explainHtml`): "Hidden on phones — turn the phone sideways? The bars stay hidden there too" reads as a confused question. Rewrite it, e.g.: "On phones (upright or sideways) the bars are hidden; their numbers are in sections 4–5 and in the curve's title."

**What works:** three bars that add to zero exactly, with the 0/0 at the sonic throat replaced by its one-sided limit (and a selftest row that checks the limit against the neighbouring stations). The ghost "same area" root on a log-M U-curve. The printed-sign line in D12 step 3, which shows that the book's + sign would forbid a supersonic nozzle.

## normal_shock_jump (E3)
**Must fix:** none.

**Should fix**
1. Explain §6, line "cubic … (15.44a)": `gapPct = 100*(cub − S.dS)/|S.dS|` is relative to the exact value, but the text says "the exact value is X % below this [the cubic]". At M₁ = 1.1 it prints "40 % below"; the exact value is 29 % below the cubic, and the cubic is 40 % above the exact value. Either divide by `cub`, or reword to "the cubic is X % above the exact value".
2. The Code-tab `ref:` labels `Eqs. (15.39)–(15.44a)` and `Eqs. (15.35)–(15.37)` are the only eq_refs hits in the chapter. They are acceptable as panel labels, but `Jump laws (15.35)–(15.37)` would read better.

**What works:** the forbidden half (M₁ < 1) is there to experiment with. The three paired flux bars stay at "1 = 1" while the signed entropy bar turns red. That is the curation's "aha", found by experiment. The series branch of ΔS near M₁ = 1 avoids cancellation (parity at M₁ = 1.001). The still-gas frame shows that T₀ depends on the frame (576 K). The water-jump mode with its head-loss bars.

## nozzle_back_pressure (E4)
**Must fix**
1. **All seven Laval preset chips light up at once** (`desktop__tour-step6.png`, `laptop__check.png`, `phone__explore.png`: "subsonic … under-expanded" are all orange/pressed). Line 2979 builds every preset with `set: { kind: 'cd' }` and puts the back pressure in `action`. The library (`viz_lib.js` l.1600) marks a chip as pressed when every key of its `set` matches the state, so all seven match. The reader cannot tell which condition is active. Fix in the explainer: `const PRE = name => ({ label: name, set: { kind: 'cd', dn: DN(name) }, action: a => a.set({ dn: DN(name, a.state.ar, …) }) })`. This highlights the right chip at the default A_e/A_t = 2 and keeps the area-aware action. Or move `kind` into the action and give `set` the `dn` key.
2. **D20 step 4, *why* (line 3032) skips a move.** It says "The jump exceeds 1, so $p_f\lt p_d\lt p_b$." Only $p_d>p_f$ follows from the jump. That $p_d<p_b$ needs its own argument, which the notebook supplies (`build_ch15.py`, `WHY_SET[("D20", 4)]`). The argument: behind the shock $p_{02}<p_{01}$ and $A_2^*=A_t\,p_{01}/p_{02}$ grows, so $A_e/A_2^*$ falls, the subsonic exit Mach number rises and $p/p_{02}$ falls. Both factors are below those of $p_b/p_0$. Port it, within 35 words, e.g.: "The jump exceeds 1, so $p_d\gt p_f$. And $p_d\lt p_b$: the shock lowers $p_0$ and enlarges $A^*$ (step 6), so the exit runs faster and at lower pressure than on the $p_b$ curve." Keep the step text identical in meaning to the notebook's.

**Should fix**
1. On phones the 7 + 1 preset chips take three rows (`phone__explore.png`). Shorter labels ("sub", "choke", "bell", "exit", "over", "design", "under") would give the pressure view more room.

**What works:** the mass-flow plateau, with the other nozzle as a grey ghost. The shock walking down the bell as the tank is pumped down (transport on `dn`). The exit-versus-tank mismatch is drawn as a coloured tick (fans or oblique shocks outside). The "choked means / does not mean" note with $\dot m\propto p_0/\sqrt{T_0}$. The parity of the shock station is closed form against `laval_nozzle`.

## duct_friction_heating (E5)
**Must fix:** none.

**Should fix**
1. Tour step 3 "One number, one quadratic" quotes the final (15.46) but links `derive: { id: 'D21', step: 11 }` ("Simplify the discriminant"). Step 12 ("Write the result") is the matching move.

**What works:** the book's sign slips in (15.47) and (15.48) are corrected and shown beside the printed versions, with the sum rule that exposes them ("the two corrected heat terms add up to dq"). The M₂(M₁) map with the ΔS < 0 region. The choke station is found by applying f and q bit by bit. The paradox band is shown live (T₂/T₁ = 0.9947 with heat added). The T–s lines carry the label "usual names, not the book's".

## mach_cone (E6)
**Must fix:** none.

**Should fix**
1. `phone__explore.png` and `phone-land__explore.png`, M = 3, t = 0: the observer icon is about 1 km tall at this scale, and its head overlaps the source's flight-path arrow. On small views, scale the person to the view or draw the U arrow above the path.

**What works:** each circle shows its emission point and age in the inspector. The arrival-time view, whose U-shaped minimum *is* the boom ("two sounds heard at once inside the cone" is a selftest row). The triangle drawn to scale in the angle view. The real-world table with the current row highlighted.

## oblique_shock_wedge (E7)
**Must fix:** none.

**Should fix**
1. `phone-tall__tour-step4.png`: the canvas play button covers the value of "② M₂ = 0.604 · p₂/p₁ = 4.44". On narrow views, drop the p₂/p₁ part of the label or anchor it left of the button.
2. D26 step titles 4 and 6 ("Momentum across; the normal Mach numbers", "Recognise the normal shock: pressure") differ in wording from the notebook ("Momentum across the shock, and the normal Mach numbers", "Recognise the normal-shock problem; pressure"). Align them (ch01 lesson: same move titles).

**What works:** the velocity triangles at a clicked streamline (blue across, amber along, unchanged). The two roots merging on the δ–σ chart, with the sonic ring and faint curves for other M₁. Exact against linear shown as terms ("the rest" in amber). Both σ, δ and β, θ are named everywhere.

## prandtl_meyer_fan (E8)
**Must fix:** none.

**Should fix**
1. `phone-tall__tour-step6.png` (compression, −8°): in the small ν(M) view the labels "−8.0° turn", "ν(M₁) = 26.38°" and "ν(M₂) = 18.38° → M₂ = 1.719" sit on the curve and the dashed guides. On views shorter than ~200 px, keep one label (the ν(M₂) line); the title already carries M₂.
2. D29 step 1 "Write the kept component with speeds and angles" and step 7 "Temperature from a fixed T₀" differ from the notebook's titles ("Start from the unchanged tangential speed", "Temperature from constant T₀"). Align them.

**What works:** the staircase of 1 → 3 → 7 → 21 → ∞ waves, with an invariant that the error halves when N doubles. A Simpson check of the integral form against the closed form (an independent route to (15.56)). The compression case with the first crossing of the Mach lines marked. Both sign conventions (our turn and the book's δ) appear in the status, readouts and terms.

## supersonic_flat_plate (E9)
**Must fix:** none.

**Should fix**
1. `phone-tall__tour-step3.png`, "Square to the plate": the text names "the drag (rose)", but at phone size neither the rose drag arrow nor its label is visible. Draw the drag component, and its label, at every size (even if exaggerated, then say so).
2. Parity row `diamond thickness drag (5 %, M 2)` uses rtol 1e-3. It compares the JS closed form with `ackeret_airfoil`'s numerical quadrature across the slope jump (actual error 3.3e-4). Either name it "closed form vs quadrature" in the row name, or give the quadrature a breakpoint at x = 0.5 so the row can be 1e-9 like the biconvex rows.
3. D31 step 1 *why* is 36 words (limit 35).

**What works:** the Cp view, where the area between the two flat levels *is* the normal-force coefficient, with Ch. 14's slow-flow loop as a ghost. The status reports the lift error (−1.1 % at 10°) while the explain panel shows that each surface is 20–25 % off. The "errors nearly cancel" lesson, backed by a second-order invariant. The diamond preset shows drag without lift.

---

### Candidates for knowledge/viz_patterns.md
- (ch15) Presets whose `set` holds only a mode key (`{kind: 'cd'}`) and do their work in `action` all light up together: the library marks a chip as pressed when every key of `set` matches. Put the distinguishing value in `set`.
- (ch15) When the book has a sign slip, show the printed form beside the corrected one with a live number and the identity it breaks (D12 step 3, E5 §5). That works better than a footnote.
- (ch15) To show "conservation laws allow it, the second law forbids it", draw paired budget bars that stay "1 = 1" next to a signed entropy bar that turns red (E3).
- (ch15) An honest exaggeration label ("shifts drawn ×500, true ≤ 0.886 mm") when the physical motion is invisible at true scale (E1).
- (ch15) Floating canvas buttons (play) sit on top of the labels in the view's top-right corner on phones. Reserve that corner (E7).

## Verdict: FAIL
8 / 9 PASS. `nozzle_back_pressure` FAIL: (1) all seven Laval preset chips light up at once (`viz/ch15/nozzle_back_pressure.html` l.2979, `PRE` sets only `{kind: 'cd'}`); (2) the D20 step 4 *why* (l.3032) derives $p_d\lt p_b$ from "the jump exceeds 1", which does not follow; port the notebook's argument (`WHY_SET[("D20", 4)]`). A re-audit of that one file is enough after the fix: `tools/viz_lint.py viz/ch15/nozzle_back_pressure.html` and `tools/shot.py viz/ch15/nozzle_back_pressure.html`.
