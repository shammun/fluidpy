# Chapter 4 — explainer review (round 2)                                  2026-09-23

History:
- Round 1: 7 PASS; E5 and E7 FAIL, one Must fix each.
- Round 2: both fixed and verified; **all nine PASS**.
- The round-1 text below is kept for the record. Its E5 and E7 Must fixes are closed (see the round-2 section).

Reviewer: viz-reviewer (fresh eyes, Phase 7 gate). What I ran and checked:
- `tools/viz_lint.py --chapter ch04`: all nine files `ok`.
- `tools/shot.py --chapter ch04`, full run (8 sizes, every tab, step and pager page): **all nine PASS with 0 failures**.
  - 2744 views in total. KaTeX rendered at every size.
  - **250/250 selftest rows ok.**
- `tools/eq_refs.py --chapter ch04`: one hit, in `which_bernoulli.html` L12. It is the `<meta name="viz:concept">` string, which is metadata and never rendered. Should fix (below).
  - No bare equation number appears in any tour, Explain, Derivation or quiz text.
- Derivations: I compared every assigned D id with the executed notebook by script. **All 30 match in step count and move titles.**
  - Counts: D01 9 · D04 8 · D05 8 · D06 11 · D08 8 · D09 13 · D10 4 · D11 7 · D12 7 · D13 9 · D14 7 · D15 12 · D16 6 · D17 5 · D18 5 · D21 6 · D22 8 · D23 8 · D24 10 · D25 5 · D26 8 · D27 9 · D28 8 · D30 8.
  - The only differences are cosmetic: "ρ_s" vs "ρs", and the notebook's "(4.56) from (4.55)" vs the explainer's words.
  - Every walkthrough `derive: {id, step}` link points to the right move.
- Book pages read:
  - p125 (4.5), p126 (4.12), p129 (4.17), p132–133 (Ex. 4.2 and (4.19)), p137 (Ex. 4.4 rocket), p141 (4.37), p144 ((4.43)–(4.45)), p151 (4.58), p157 ((4.73)–(4.78)), p163 (4.86), p173 ((4.100)–(4.104)).
  - Every Equations card and derivation line I checked matches its page.
- Book-specific checks the gate asked for:
  - **E5**: the (4.43) +2Ω×u′ and (4.45) −2Ω×u′ signs are right in the physics, the bars, the arrows and the derivations. The side-panel term labels are not (Must fix below).
  - **E1**: the jet split is A(1 ± cos θ)/2. The larger sheet goes along the side of the plate that points downstream. Along-plate momentum closes to 1e-12.
  - **E9**: Ro = U/(2Ωl) is labelled a Ch. 13 preview in the tour, the Equations tab, the regime map and the JS header.
  - **E6**: the (4.74) gauge is taught as φ − ∫B dt′. The book prints + (checked on p157). Both signs are shown with live numbers in D26 step 6 and quiz 4.
  - **E8**: the `SCEN` table is identical to `fluidpy.ch04_conservation_laws.BOUSSINESQ_SCENARIOS`. The 20 scenario rows are parity-checked.
  - **D09 (E3)** is complete: 13 steps, including the γ = μ naming step and the Newton-law check.
  - **D15 (E5)** is complete: 12 steps. The two Ω × u′ come from steps 6 and 9, and the triple product is step 10.

## Summary

| slug | lint | shot (sizes/views) | parity rows ok | explain | derivations (steps ok) | teaching | polish | verdict |
|---|---|---|---|---|---|---|---|---|
| control_volume_budgets (E1) | ok | PASS · 8 · 264 | 37/37 (1e-12…1e-8; erf-based quadrature parity 1e-8; measured ◇ 1e-6) | 8 sections + interpretation, 4 live, face table | D01 9 ✓ · D05 8 ✓ | very good | Explore 10 pp and step 4 6 pp at 360×640 | **PASS** |
| stream_function_spacing (E2) | ok | PASS · 8 · 176 | 20/20 (1e-12 closed forms; 1e-8/1e-9 Gauss–Legendre and FD) | 7 sections + interpretation | D04 8 ✓ | very good | step 4 6 pp at 360×640 | **PASS** |
| newtonian_stress_lab (E3) | ok | PASS · 8 · 344 | 17/17 (all 1e-12) | 9 sections + interpretation, 3 live | D08 8 ✓ · D09 13 ✓ · D10 4 ✓ | excellent | Explain 28 pp at 360×640 | **PASS** |
| navier_stokes_term_balance (E4) | ok | PASS · 8 · 328 | 36/36 (1e-12 closed forms; 1e-5/1e-6 vs fluidpy stencils, justified) | 7 sections + interpretation, 5 live | D11 7 ✓ · D12 7 ✓ · D13 9 ✓ | very good | Explore 14 pp, Explain 14 pp at 360×640 | **PASS** |
| rotating_frame_coriolis (E5) | ok | PASS · 8 · 456 | 31/31 (1e-12; DOP853 path 1e-7; 4 exact-text label rows) | 8 sections × 3 mode variants + interpretation | D14 7 ✓ · D15 12 ✓ · D16 6 ✓ · D17 5 ✓ · D18 5 ✓ | excellent | labels signed per side ✓ (round 2) | **PASS** (round 2) |
| which_bernoulli (E6) | ok | PASS · 8 · 432 | 30/30 (1e-12; exact-text `which` rows) | 7 sections + interpretation | D06 11 ✓ · D24 10 ✓ · D25 5 ✓ · D26 8 ✓ | excellent | step 2 6 pp at 360×640 | **PASS** |
| viscous_dissipation_heating (E7) | ok | PASS · 8 · 312 | 24/24 (1e-12; series 1e-8) | 8 sections + interpretation, 4 live | D21 6 ✓ · D22 8 ✓ · D23 8 ✓ | very good | D22 TeX fixed ✓; status and units ✓ (round 2) | **PASS** (round 2) |
| boussinesq_buoyancy (E8) | ok | PASS · 8 · 256 | 36/36 (1e-12) | 9 sections + interpretation, 4 live | D27 9 ✓ · D28 8 ✓ | very good | Explain 24 pp at 360×640 | **PASS** |
| dynamic_similarity_models (E9) | ok | PASS · 8 · 176 | 24/24 (1e-12/1e-10; one 2e-3 asymptotic invariant) | 8 sections + interpretation, 3 live | D30 8 ✓ | very good | minor | **PASS** |

"pp" = pager pages. At 390×844 every Explain is 5–9 pages, and in the 1000×700 notebook frame 4–6.

**Checks that apply to all nine.**
- Tabs present: Walkthrough · Explore · Explain · Derivation · Equations · Code · Check yourself.
- Each has 6–8 depth features: linked views, transport, presets, status, inspector and modes everywhere; term bars in seven of the nine.
- 4 check questions each, all answerable with a preset or a `set`.
- Step 1 asks a plain-words question in every explainer.
- Audit: no overflow, no text below 12 px and no view below 60 px at any size.
- Every `py:` row calls `ch04.*` at non-trivial inputs against the explainer's own JS function (one exception in E5, Should fix).
- No book prose, figures or tables are reproduced. The ship, sphere data, blob and parcel models are synthetic and labelled "ours".

**About the long pagers at 360×640 (the builders asked).** E1 Explore is 10 pages and E4 Explain is 14 (it was 20 in the builder's run). E3 Explain is 28 and E8 Explain 24.
- These are acceptable. They are the library's designed safety net, and they drop to 5–9 pages at 390×844.
- None is a Must fix.
- A walkthrough card that needs 5–6 pages at 360×640 is too long, though (E1 step 4, E2 step 4, E6 step 2, E7 steps 2 and 5, E8 step 3). Those are listed as Should fixes.
- The cheapest fix is to move the inline equation to the step's `eq:` card and keep the text to about 30 words.

## control_volume_budgets (E1)
**Must fix** — none.

**Should fix**
1. `phone__tour-step4` is 6 pages at 360×640. The step text plus the inline (4.17) run to about 55 words with a display-width formula. Shorten the text to about 30 words and move (4.17) into `eq: 'mom'`.
2. In `phone__explore.png` (360×640) the Budget view's "±120" / "±3.65" scale labels sit on the storage ◇. Offset the label right of the first bar, or drop it when the view is shorter than 70 px.
3. The Explore pager is 10 pages at 360×640. It would shrink if the `budget` mass/momentum chips joined the mode strip (the same pattern as the scene chips).

**What works** (keep these; candidates for knowledge/viz_patterns.md)
- Five scenes on one budget object: the rule changes the box, not the law.
- The measured ◇ (finite-difference storage) lands on the formula bar in the bore, rocket and balloon.
- The residual bar is at round-off.
- The bore's "drag b until storage vanishes" is a genuine discovery interaction.
- The jet derives the (1 ± cos θ)/2 split from along-plate momentum in Explain §4.
- D05 step 3 flags the book's stray "= 0" in (4.15).

## stream_function_spacing (E2)
**Must fix** — none.

**Should fix**
1. Tour step 4 ("Velocity is a slope") is 6 pages at 360×640 because the inline ρu, ρv and the derive quote stack up. Keep the text and let the derive quote carry the formula.
2. On a phone in vortex mode (`phone__explore.png`) the x axis shows only the 0 tick. Add ±2 so the reader can judge where the contours crowd.

**What works**
- The continuous-θ unwrapping across the source cut.
- The "Spacing = speed" view: amber bars of width Δn and height Δψ/Δn against |u|, with a note that the width×height of every bar is Δψ.
- The axisymmetric 2πΔψ variant, with the correct "crowd away from the axis" reading. This silently fixes the storyboard's "toward the axis".

## newtonian_stress_lab (E3)
**Must fix** — none.

**Should fix**
1. The Explain pager is 28 pages at 360×640, and page 1 holds only the "0 What the windows show" heading (`phone__explain.png`). Shorten §0 to two sentences. In cube mode, open Explain at §7 (the cube section) instead of the stress-lab sections the reader cannot see.
2. At 390×844 the D08/D09/D10 derivation chips wrap into a tall three-line oval (`phone-tall__derive-d2p11.png`). Use the short labels (`D08 · τ₁₂=τ₂₁` → `D08`) on phones.

**What works**
- The source-coloured τ grid (−pδ orange, 2μ(S − ⅓S_mmδ) rose, μ_vS_mmδ purple), which switches to λ/μ/γ colours during D09.
- D09's rose/amber "only μ + γ acts" watch.
- The spinning cube's log–log α ∝ 1/h² with its end card.
- The terms bar that splits σn by source and sums exactly (selftest invariant).

## navier_stokes_term_balance (E4)
**Must fix** — none.

**Should fix**
1. Explore is 14 pages at 360×640. Mark `visc3` and `Re` as optional readouts on phones, since the status and the bars already carry them.
2. In `desktop__tour-step4.png` (Stokes, playing), the view title says "t = 126 ms" while the transport and the layer label use t = 210 ms (2√(νt) = 0.916 mm).
   - Check whether view titles are throttled during play. If so, drop the time from the title and keep it in the transport.
3. Tour steps 5 and 6 are 5 pages each at 360×640. Trim the text and let the `derive` quote carry (4.40).

**What works**
- Seven exact solutions with analytic derivatives, cross-checked against fluidpy's stencils.
- The two-column "ρDu/Dt = forces" bar next to log-scaled term rows.
- The regime word in the status (pressure ↔ viscous, local ↔ viscous, inertial, inviscid-idle).
- The three viscous forms drawn as coincident arrows.
- The compressible toy field is used only in D12, so the purple ∇(∇·u) term has something to show before vanishing.

## rotating_frame_coriolis (E5)
**Must fix**
1. **The side-panel term labels contradict their signed values.** This is the (4.43)/(4.45) trap the chapter must teach.
   - In `desktop__tour-step5.png` (side = force) the terms panel reads "Coriolis $2\boldsymbol\Omega\times\mathbf u'$ = 3.35 m/s²" and "centrifugal $\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')$ = −0.335 m/s²".
   - The bars beside them correctly read "−2Ω×u′ Coriolis 3.35" and "−Ω×(Ω×x′) centrifugal −0.335".
   - The values are the (4.45) forces. The labels are the (4.43) accelerations, whose across-components are −3.35 and +0.335.
   - Source: `CFG.terms.items`, L2905–2913. The labels are static; the values come from `itemVal(s, …)`, which follows `s.side`.
   - Fix: make the labels follow the side, for example "−2Ω×u′ Coriolis" / "−Ω×(Ω×x′) centrifugal" / "−dU/dt frame" / "−Ω̇×x′ angular" on the force side and "+2Ω×u′ …" on the acceleration side, as `barItems()` already does.
   - If the terms API only takes static labels, use sign-neutral labels plus a panel title that names the side, e.g. "(4.45) forces, −2Ω×u′ …" or "(4.43) accelerations, +2Ω×u′ …".
   - Pin the fix with an exact-text selftest row (the ch01 E4 pattern).

**Should fix**
1. The selftest row "centrifugal potential" tests a JS literal, `-0.5 * 0.5 * 0.5 * 3 * 3`, instead of an explainer function.
   - Add `centPot(R, Om) = −½Ω²R²` and use it in Explain `explainEff` §2 ("its potential") and D18.
   - Then parity-test that function against `ch04.centrifugal_potential`.
2. In D15 step 6 (the "first Ω × u′"), the half Coriolis arrow is barely visible at t = 0.5 s on phones (`phone-tall__derive-d2p6.png`). Set τ to about 0.5 (t = 1 s), or draw the half arrow with a minimum length of 18 px.
3. Tour steps 2, 4 and 5 are 5 pages each at 360×640. The inline (4.43)/(4.45) in steps 4–5 is display-wide. Put it in `eq:` and keep the words.

**What works**
- Two observers on one clock.
- The "chalk mark" curve in the inertial view.
- The side toggle "forces (−) / accelerations (+)", which flips every bar and arrow.
- D15 draws the Coriolis arrow in two halves (steps 6 and 9): the "why the 2" moment really clicks.
- Pole mode reproduces 9.34 km exact vs 9.45 km small-angle.
- The f-plane high/low parcels.
- The effective-gravity cross-section with an honest ×-exaggeration note.

## which_bernoulli (E6)
**Must fix** — none.

**Should fix**
1. Line 12, `<meta name="viz:concept">`, cites "(4.19) … (4.71) … (4.72) … (4.75) … (4.78)" bare. It is not rendered, but it is the only `eq_refs` hit in the chapter. Reword it to "…along a streamline, on Lamb surfaces, everywhere, unsteady potential flow, energy form" without the numbers, or add the TeX.
2. The U-tube (`desktop__tour-step6.png`) has two labelling problems:
   - The stacked label above bar B reads "−0.4905", which is B + ∂φ/∂t. The view title says "at B 0.4905", which is B alone.
   - The black total marker line crosses the "∂φ/∂t" text.
   - Fix: label the stack "B + ∂φ/∂t = −0.4905" and move the in-bar label clear of the marker.
3. Tour step 2 is 6 pages at 360×640. Shorten the text; the derive quote D06 step 10 already carries the numbers.

**What works**
- The TRUTH-vs-your-hypotheses decision table, with teal "holds here" and amber "you claimed something the flow does not have".
- The status names the valid form with its TeX.
- The whirlpool dip read directly from B(0).
- The U-tube's amber ∂φ/∂t term makes B + ∫∂u/∂t ds flat.
- D26 shows the gauge sign slip with live numbers (0 vs 2B).
- D06 fills in the book's skipped side-pressure move.

## viscous_dissipation_heating (E7)
**Must fix**
1. **Garbled KaTeX in D22 (goal page and step 1).** `start.plain` (L2867) and step 1 `why` (L2872) hold the inline (4.55) with single backslashes inside single-quoted JS strings:
   - `'…$\rho D(e+\tfrac12u_j^2)/Dt=\rho g_iu_i+\tau_{ij}\partial_iu_j+u_j\partial_i\tau_{ij}-\partial_iq_i$…'`
   - `\r` becomes a carriage return, `\t` a tab, and `\p` a plain "p".
   - The page renders "hoD(e+frac12u_j^2)/Dt=hog_iu_i+au_{ij}partial_iu_j+…" (`desktop__derive-d2p0.png`, `desktop__derive-d2p1.png`).
   - Fix: double every backslash in both strings, or reuse the `EQ55` constant. A scan of all nine files found no other single-backslash TeX.

**Should fix**
1. The steady status reads "steady: work in = heat out = 0.998 W/m²" while work in = 1.000 W/m² and 0.2 % is still being stored at κt/h² = 0.6 (`phone-tall__tour-step2.png`). Write "work in 1.000 ≈ heat out 0.998 W/m² (0.2 % still stored)", or run the clock to κt/h² ≈ 1.
2. The wall heat-flux labels in the Temperature view ("0.499") have no unit. Add W/m².
3. Tour steps 2 and 5 are 6 pages each at 360×640. Move (4.54) to `eq:` in step 2, and shorten step 5.

**What works**
- The energy budget in = stored + out, with the share-of-work curve over κt/h².
- The "μ < 0 ⛔" preset, which turns ε and the entropy production negative and trips the status.
- The 2-D dissipation keeps the missing z deviator (−S_mm/3), so the compressible parity rows are exact.
- D23 lights the Poiseuille ρε ∝ (du/dy)² picture on its "square" step.

## boussinesq_buoyancy (E8)
**Must fix** — none.

**Should fix**
1. The Explain pager is 24 pages at 360×640. §6 and §7 are long. Split §7 (C_p) into the D28 link, and on phones drop §8's five-row list (the status already names the verdict).
2. Tour step 3 is 6 pages at 360×640: the terms panel plus the derive quote. Drop `terms: true` on phones, since the bars view shows the same numbers.
3. The per-setting rows compute `p.alpha * p.dT`, `p.L * 9.81 / p.c²` and so on inline. Route them through `validity(...)` so the rows exercise the app's own function. The fluidpy comparison is already correct and is what the gate asked for.

**What works**
- The honest "our parcel model" labelling everywhere.
- The grey "dropped" bar with its ×αδT ratio.
- The four-number validity chart with the five-setting table and the current row lit.
- The buoyancy-off ghost circle.
- D28's live "pressure work folds into C_p" bar story.

## dynamic_similarity_models (E9)
**Must fix** — none.

**Should fix**
1. In ship mode with Fr unmatched (`desktop__tour-step1.png`), the model panel shows no waves while its label says "wave length 2πU²/g = 64 m". The model's waves are 16 l long, off the drawing, which is the point. Say so: "wave length 64 m = 16 l — off the picture".
2. Tour step 2 is 5 pages at 360×640 because the inline (4.101) is display-wide. Use `eq: 'nd'` and keep the words.

**What works**
- Prototype and model drawn in their own l, one t* clock.
- The paired log bars with "✔ matched / ×5 / ÷125" badges.
- Forty synthetic spheres collapsing onto Morrison's curve, with a toggle back to F vs U.
- Ship drag split into friction and wave, with the naive bar next to it.
- The Ro map, clearly flagged as a Ch. 13 preview.

## Round 2 (E5, E7)

Re-run on these two files only; the other seven are unchanged since round 1.
- `viz_lint`: both `ok`.
- `shot.py` (full, 8 sizes): both PASS with 0 failures.
  - E5: 456 views, **31/31** selftest rows.
  - E7: 312 views, **24/24** selftest rows.
  - KaTeX rendered at every size.
- Derivation step counts and titles still match the notebook: D14–D18 and D21–D23.
- The single-backslash TeX scan over all nine files finds nothing.

**E5 rotating_frame_coriolis — the Must fix is closed.**
- `desktop__tour-step5.png` (force side): the side panel now reads "frame −dU/dt 0", "Coriolis −2Ω×u′ 3.35", "angular −Ω̇×x′ 0", "centrifugal −Ω×(Ω×x′) −0.335", "sum 3.02". These agree with the bars.
- The acceleration side (tour step 4 bars) reads +2Ω×u′ −3.09, +Ω×(Ω×x′) 0.091.
- Four new exact-text rows pass:
  - "term label text: force side"
  - "…: acceleration side"
  - "term label sign = sign convention of the value"
  - "rendered side-panel label follows the side"
- `phone-tall__tour-step5.png` shows the signed labels too. Walkthrough pages at 360×640 are now 1/4/4/4/4/2/2 (were up to 5).
- Open Should fixes: round-1 items 1 (the centrifugal-potential row still tests a literal) and 2 (the D15 step 6 half arrow on phones).

**E7 viscous_dissipation_heating — the Must fix is closed.**
- `desktop__derive-d2p0.png` and `desktop__derive-d2p1.png` now render (4.55) cleanly in the goal's *in words* and in step 1's *why*: $ho D(e+	frac12u_j^2)/Dt=ho g_iu_i+	au_{ij}\partial_iu_j+u_j\partial_i	au_{ij}-\partial_iq_i$.
- Round-1 Should fixes 1 and 2 are also done:
  - The status reads "nearly steady: 99.8 % of the work in (1 W/m²) leaves as heat".
  - The wall fluxes read "0.499 W/m²".
- Should fix 3 is done: walkthrough pages at 360×640 are now 1/4/4/4/2/4 (were up to 6), and step 2 moves (4.54) to its card (`phone-tall__tour-step2.png`).
- No open Should fixes.

**Open items (Should fix, not blocking).** All round-1 Should fixes of E1, E2, E3, E4, E6, E8 and E9 stay open as listed above, including:
- the E6 `viz:concept` meta line, still the only `eq_refs` hit;
- the long 360×640 walkthrough cards in E1, E2, E4, E6, E8 and E9.

## Verdict: PASS (all explainers PASS)
Round 2: all nine explainers pass.
- 2744 views audited, 0 failures.
- **255/255** selftest rows ok (250 in round 1 plus the 5 new E5 rows).
- All 30 derivations match the notebook.
