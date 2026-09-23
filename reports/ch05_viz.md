# Chapter 5 — explainer review (round 2)                                  2026-09-23

History:
- Round 1: six PASS; E2, E3 and E7 FAIL, one Must fix each.
- Round 2: all three fixes verified; **all nine PASS** (see the round-2 section at the end).
- The round-1 text below is kept for the record. Its E2, E3 and E7 Must fixes are closed.

Reviewer: viz-reviewer (fresh eyes, Phase 7 gate).

## What I ran and checked

**Mechanical checks**
- `tools/viz_lint.py --chapter ch05`: all nine files `ok`.
- `tools/shot.py --chapter ch05`, full run (8 sizes, every tab, step and pager page). I re-ran it and waited for it to finish; all nine audits are dated 17:19–17:30 today.
  - **All nine PASS, 0 failures, 2776 views.**
  - KaTeX rendered at every size.
  - **268/268 selftest rows ok.**
- `tools/eq_refs.py --chapter ch05`: **0 hits**. No bare equation number appears in any tour, Explain, Derivation or quiz text.

**Derivations**
- I compared every assigned D id with the executed notebook `notebooks/ch05_vorticity_dynamics.ipynb` by script. **All 23 match in step count and move titles**:
  - D01 7 · D02 10 · D03 9 · D04 8 · D05 9 · D06 10 · D07 7 · D08 9 · D09 12
  - D10 12 · D11 12 · D12 6 · D13 8 · D14 8 · D15 15 · D16 6 · D17 5 · D18 9
  - D19 10 · D20 6 · D21 6 · D22 5 · D23 5
- The only differences are wording:
  - "Pull kernel and direction out." vs "…out of the integral".
  - D14 step 5: the notebook's "(5.22)" vs the explainer's "the ε form".
  - D20 step 6: the notebook's "against N39".
- Every walkthrough `derive: {id, step}` link lands on the right move. For example, E1 goes to D01 step 6 (the lids) and E6 goes to D10 step 11 (the sign).

**Book pages read**
- p201: (5.6), and Fig. 5.2 labelled "2ω".
- p208: (5.13), and (5.14) printed with −1/(4π).
- p209: (5.15)–(5.17), and the second sign slip in the product-rule line.
- p215: the Fig. 5.11 caption.
- p218: §5.8 dΓ = (u₂ − u₁)ds, and the Fig. 5.16 caption u₁ − u₂.

**Chapter-specific checks, all satisfied**
- **E6 / D10 / D11: (5.14) is shown with +1/(4π) and the printed −1/(4π) is flagged.** This appears in five places:
  - the Equations card note;
  - D10 step 11, with live printed and corrected u_θ;
  - tour step 2;
  - quiz 4;
  - the sign toggle, with a "⚠ printed" status.
- **E6 / D11 step 5: the book's second slip is flagged.** The book prints +(x − x′)/|x − x′|³ × ω; the page shows it in rose next to the correct line, and says the two slips cancel. I checked the algebra on paper against p209.
- **E6 / D12 step 3 has the corrected estimate.** It reads 2a/|x − x′|, 20 % at 10 core radii, with a live line 200a/r %.
- **D09 (E5) is complete, 12 steps.** It includes the ∇(u·u) → ∇(½u·u) remark.
- **D15 (E7) is complete, 15 steps.** It includes the dropped u_j,j(ω_n + 2Ω_n) term (step 9) and the Π→Φ slip (step 3). I checked the δ contractions in steps 4, 5 and 12.
- **D23 step 2 (E9) is corrected.** It says the short sides leave a remainder of size dn·ds, and step 4 says they vanish as dn → 0. This is right for a finite sheet: v is odd in x, so the two sides add, but only to O(ds·dn).
- **E4 baroclinic sign convention is correct**, in the physics, status, badge text, quiz 4 and the selftest invariant.
  - ∇p = (0, −ρ₀g) and ∇ρ = |∇ρ|(sin θ, −cos θ).
  - θ = 0 is stable; 0 < θ < 180° spins clockwise ("crossed 90°: clockwise", −0.0981 s⁻²); θ = 270° spins counterclockwise.
  - Lock exchange with heavy fluid on the left: (ρ₁ − ρ₂)/δ·e_x × (−ρ̄g)e_y > 0, so counterclockwise, 2.42 s⁻².
- **E2, Fig. 5.2:** ω = 2Ω. The "2ω" slip is named in Explain §1 and in the D02 start text. D02 step 7 warns that using ω as the rate would make the pressure four times too big.
- **E3, Kelvin needs ρ = ρ(p):** the hypothesis table, D05 step 6, Explain §5 ("single-valued is not enough") and tour step 4 all say so.
- **E8, Fig. 5.11:** G is the centre of vorticity. The fluid at G moves at −1.70 m/s for Γ₂ = 3Γ₁ (I checked this). This appears in an Explain §3 hint, the D21 step 4 watch line and an equation note.
- **E9, sheet sign:** γ = u₂ − u₁ with counterclockwise positive. The caption's u₁ − u₂ is shown as the clockwise count, with a toggle, a status line, Explain §6 and D23 step 5.
- **Ring speed and pair rates are form checks, not benchmarks.** `ringSelf` is labelled "a V1 form cross-check" and the word "benchmark" appears in no explainer.
- **E8 guards against singular drags:**
  - dropped and dragged vortices are kept ≥ 5 cm apart (`separate`);
  - positions are clamped inside the wall or bucket;
  - the run halts at a 4 mm gap (`TOUCH`);
  - the step is sub-divided when vortices are close.

**Physics re-derived by hand from the JS**
- E1: ∇·ω = 0 for the narrowing Gaussian tube, including ω_R.
- E2: the tornado edge and axis deficits, 609.9 and 1220 Pa.
- E4: the D06 closed forms x_G = R²∇ρ/4ρ₀ and M_G/I_G.
- E5: the Burgers balance, adv = str + dif.
- E6: the D13 substitution l = −d cot θ, giving 0.1125 m/s.
- E7: the column (−0.2f = −2.06×10⁻⁵ s⁻¹) and the ring (−4.19×10⁷ m²/s).
- E8: h₁ = Γ₂h/(Γ₁ + Γ₂).

## Summary

| slug | lint | shot (sizes/views) | parity rows ok | explain | derivations (steps ok) | teaching | polish | verdict |
|---|---|---|---|---|---|---|---|---|
| vortex_tubes_cannot_end (E1) | ok | PASS · 8 · 168 | 27/27 (1e-12 closed forms; RK4 vs DOP853 1e-7/1e-6) | 8 sections + interpretation, 4 live, 5 field variants | D01 7 ✓ | very good | step 5 is 6 pp at 360×640 | **PASS** |
| vortex_pressure_funnel (E2) | ok | PASS · 8 · 280 | 34/34 (1e-12; FD slopes 1e-7; volume 1e-15) | 10 sections + interpretation, table of real vortices | D02 10 ✓ · D03 9 ✓ | excellent | tornado y-axis hides the minus sign (Must fix); steps 3/5/6 are 5 pp at 360×640 | **PASS** (round 2) |
| kelvin_material_loop (E3) | ok | PASS · 8 · 360 | 29/29 (1e-12…1e-6; Rankine kink 1e-5, justified) | 7 sections + interpretation, 14 live | D04 8 ✓ · D05 9 ✓ · D08 9 ✓ | excellent | baroclinic step overlaps at 390×844 (Must fix) | **PASS** (round 2) |
| baroclinic_torque (E4) | ok | PASS · 8 · 264 | 24/24 (1e-12; numeric disc 1e-6/1e-9; stencil 1e-5) | 8 sections + interpretation, element and lock variants, real-case table | D06 10 ✓ · D07 7 ✓ | excellent | minor | **PASS** |
| vorticity_stretching_tilting (E5) | ok | PASS · 8 · 416 | 47/47 (1e-12; stencil 1e-5/1e-4; exact-text status rows) | 8 sections + interpretation, 14 live | D09 12 ✓ · D16 6 ✓ · D17 5 ✓ · D18 9 ✓ | very good | one clipped label | **PASS** |
| biot_savart_filament (E6) | ok | PASS · 8 · 464 | 23/23 (1e-12/1e-10; tube quadrature 1e-7; discretisation invariants 2e-4/1e-3) | 8 sections + interpretation | D10 12 ✓ · D11 12 ✓ · D12 6 ✓ · D13 8 ✓ | excellent | minor | **PASS** |
| vorticity_equation_rotating (E7) | ok | PASS · 8 · 472 | 31/31 (1e-12; stencils 1e-5…1e-7) | 5–6 sections per mode + interpretation, highlighted f-table | D14 8 ✓ · D15 15 ✓ · D19 10 ✓ · D20 6 ✓ | very good | Terms view illegible at 390×844 (Must fix); D15 result page 25 pp at 360×640 | **PASS** (round 2) |
| point_vortex_lab (E8) | ok | PASS · 8 · 208 | 32/32 (1e-12; RK4 vs DOP853 1e-6/1e-8) | 8 sections + interpretation, 11 live | D21 6 ✓ · D22 5 ✓ | excellent | invariants legend on phones | **PASS** |
| vortex_sheet_rollup (E9) | ok | PASS · 8 · 144 | 21/21 (1e-12; RK4 vs DOP853 1e-4/1e-5; L1 rows 1e-3, loose) | 8 sections + interpretation, jump and roll-up variants | D23 5 ✓ | very good | minor | **PASS** |

"pp" = pager pages. At 390×844 every Explain is 4–7 pages, and in the 1000×700 notebook frame 3–5.

**Checks that apply to all nine**
- Tabs present: Walkthrough · Explore · Explain · Derivation · Equations · Code · Check yourself.
- Each has 7–8 depth features: linked views, transport, presets, status, notes and modes everywhere, plus an inspector in eight and term bars in seven.
- 4–5 check questions each, all with a `set`, all answerable on the picture.
- Step 1 asks a plain-words question in every explainer.
- Every `py:` row calls `ch05.*` at non-trivial inputs against the explainer's own JS function.
- No book prose, figures or tables are reproduced. E2 quotes four words ("hyperboloids of the second degree"), with attribution; that is acceptable.

**About the long pagers at 360×640**
- Explain pagers of 14–19 pages (E2 19, E6 16, E1 15, E8 15, E4 14) drop to 4–7 pages at 390×844. That is acceptable, as in ch04.
- Walkthrough cards of 5–6 pages at 360×640 are too long: E1 step 5 (6 pp) and E2 steps 3, 5 and 6 (5 pp each). They are listed as Should fixes.
- **E1 step 5, judged as asked:** its 22-word text is fine. The pages come from the D01 quote plus a 3-line code excerpt. It is 4 pages at 390×844.

## vortex_tubes_cannot_end (E1)
**Must fix** — none.

**Should fix**
1. `phone` tour-step 5 ("Why: Gauss") is 6 pages at 360×640 (4 at 390×844). Drop `code: {lines: [5, 7]}` from that step; the D01 quote already carries the idea. The code highlight could move to step 4.
2. The ω̄ curve label is a combining macron ("ω̄"). At 12 px it renders as "ω¯" with the bar offset (`phone-land__explore`). Use "mean ω" or draw the bar on the canvas.

**What works** (keep these)
- The five-field chip set (narrowing · twisted · Burgers · ring · broken ⚠) makes "cannot end" a *contrast*:
  - the broken field's sum −0.40 m²/s equals ∫∇·ω dV (selftest row);
  - the theorem fails exactly where its hypothesis fails.
- The Gauss numbers are always visible in the scene caption (lower / wall / upper / sum) in the bar colours, so phones lose nothing when the budget view is hidden.
- The inspector computes ω_R − R′ω_z = 0 on the wall for a clicked line.
- The Burgers mode shows the stream tube narrowing while the vortex tube does not (N45).

## vortex_pressure_funnel (E2)
**Must fix**
1. **At 390×844 the tornado section's rotated "z [m]" title covers the minus sign of the "−100" tick**, which then reads "100" below 0 (`phone-tall__tour-step6.png`, left axis).
   - This is the ch01 lesson on narrow views.
   - Fix: when the view is narrower than 420 px, draw the rotated y title in its own strip left of the ticks, or move the unit into the view title ("z [m]") and drop the rotated title.
   - Check every mode at 360 and 390 px: the line and Rankine sections have negative z ticks.

**Should fix**
1. Walkthrough steps 3, 5 and 6 are 5 pages each at 360×640.
   - Step 3 quotes D02 step 10, whose two-line `aligned` result is tall. Quote step 9 (one line) or drop the `readouts`.
   - Steps 5 and 6: keep one extra each.
2. At 360×640 (`phone__explore.png`) the three isobar labels "+200 / +400 / +600 Pa" stack and touch. Label only the middle isobar when the view is under 200 px tall, or stagger them in r.

**What works**
- Four vortices on one pair of axes, the "Forces on the probe" bars (needed = supplied) and the finite-difference slope of p(r) shown next to ρu²/r ("the pressure field really does the pushing").
- D03 does the net viscous force three ways, with the r² lever-arm trap named.
- The real-vortex table has the current row highlighted.
- The tornado 610 / 1220 Pa numbers are consistent everywhere.

## kelvin_material_loop (E3)
**Must fix**
1. **Walkthrough step 6 "Break it: density" at 390×844** (`phone-tall__tour-step6.png`) has two overlaps.
   - (a) In the loop view, the "light 1000 kg/m³" label is drawn over the "heavy 1025 kg/m³" label. It reads "heavy light 1000 kg/m³", which suggests the heavy side is 1000.
   - (b) In the Γ(t) view (about 45 px of plot), the y ticks "0.04 / 0.02 / 0" overlap each other.
   - Fixes:
     - In the baroclinic flow on narrow views, put the heavy label top-left and the light label top-right (`Viz.text` with `bg: true`), or shorten both to "ρ₂ heavy" / "ρ₁ light".
     - In the Γ(t) view, draw only the 0 tick and the current value (the library lesson for views under 60 px), or give that step's `gamma` row more height on portrait.

**Should fix**
1. The Helmholtz-mode Γ(t) title prints "Γ = −1.036×10⁻¹⁶" (`desktop__tour-step7.png`). Print "Γ ≈ 0 (round-off 10⁻¹⁶)" with `Viz.tnum`, so round-off is not read as a result.
2. The Γ(t) title is truncated on phones ("purple = …"). Drop the legend words there; the colours are named in Explain §0.

**What works**
- Six flows and a fixed loop on one clock. The hypothesis table has a ✓/✗ per restriction and the surviving (5.10) term as a bar.
- The measured ◇ dΓ/dt lands on the formula, and in the fixed-loop mode visibly misses it ("the bars are for material loops").
- The spectral (FFT) loop integral keeps Γ at 1.155 m²/s to round-off while the loop grows ×6.7.
- The "single-valued is not enough, ρ = ρ(p) is" correction in D05 step 6.

## baroclinic_torque (E4)
**Must fix** — none.

**Should fix**
1. The Explore pager is 10 pages at 360×640. Mark `grho` (log slider) or `fluid` `optional` in element mode.
2. Tour step 4 ends "…the baroclinic term of (5.28) in the card". The formula is written earlier in the sentence, so this passes, but "…the baroclinic term, card (5.28)" would be clearer.

**What works**
- The torque route (◇, 2M_G/I_G) rides the formula's sine, and the third view shows the gap falling as R² (slope 2). A reader *sees* why the point-wise formula is exact in the limit.
- The θ = 270° "order matters" step and quiz pin the ∇ρ × ∇p order.
- A rim-point inspector gives each push's moment about G.
- The lock exchange's tank circulation rate ρ̄gH(1/ρ₁ − 1/ρ₂) is independent of δ, leading to the vortex sheet (C14).
- A real-case table (lab, sea breeze, ocean front) has "your setting" highlighted.

## vorticity_stretching_tilting (E5)
**Must fix** — none.

**Should fix**
1. In Burgers mode the scene label "outflow u_z = αz" is clipped at the bottom edge (`desktop__explain.png`). Keep it inside the view, or drop it below 300 px of height.
2. The Code pager is 15 pages at 360×640. It is acceptable, but the `burg` snippet's two-line call could become one line.

**What works**
- Stretching (purple) and tilting (blue) arrows on a 3-D line, with term bars split along e_s and across e_n.
- Linear flows use closed-form e^{Gt}, and the custom G has a scaling-and-squaring expm with parity to 1e-11.
- The transient Burgers core relaxing 4 mm → 2 mm shows that stretching and diffusion *settle* to a balance, rather than just asserting it.
- Exact-text status rows (stretching / tilting / 2-D).

## biot_savart_filament (E6)
**Must fix** — none.

**Should fix**
1. While the build transport plays, the unrolled-view title's "sum" and the scene's "so far" disagree by one piece (`phone-tall__tour-step4.png`: 0.001843 vs 0.003868 m/s). Compute both from the same `builtCount`.
2. At 844×390 (`phone-land__explore.png`) the unrolled view's x title "arc length s along the filament [m]" is clipped at the right edge. Use "arc length s [m]" when narrow.

**What works**
- "Build the sum" transport: purple pieces laid tip to tail onto the closed form.
- The unrolled integral has its area equal to the speed.
- The tube-mode sign toggle reverses the whole u_θ(r) curve (rose), which makes the (5.14) slip impossible to miss.
- D11 step 5 shows the book's second slip in rose beside the correct line.
- The D12 kernel-freeze estimate is live.

## vorticity_equation_rotating (E7)
**Must fix**
1. **Tour step 3 "Stretch the planet's spin" at 390×844** (`phone-tall__tour-step3.png`): the Terms view has about 120 px for seven rows.
   - All seven labels (local, advective, relative, planetary, baroclinic, diffusion, residual) are drawn on top of each other.
   - The values "0.2 / 0 / 0 / 0.2 / 0 …" collide.
   - This is the view the step points at ("the amber bar"), and it is illegible.
   - Fix options, any one:
     - in budget mode on portrait, hide the schematic `layer` view (it carries no number) so Terms gets the stage; the Derivation tab already does this and is legible (`phone-tall__derive-d2p9.png`);
     - or, below about 25 px per row, draw only the nonzero bars and the residual, with a "zero: advective, baroclinic, diffusion" line.

**Should fix**
1. The D15 result page ("the whole chain") is 25 pages at 360×640, 9 at 390×844 and 8 even in the 1000×700 notebook frame. Show the chain as the 5 key lines (steps 1, 5, 10, 13, 15) and link "all steps".
2. The walkthrough steps 2, 3 and 6 are 4 pages at 360×640. That is acceptable. Step 2's 34-word text could lose its last clause ("— each a bar here").

**What works**
- Three modes (column over a ridge, ring moved poleward, term budget) with the conserved quantity drawn flat in amber against the changing ζ or Γ.
- The f-by-latitude table has the current row lit.
- The budget scenes each isolate one new term: planetary in the stretched column, baroclinic in the lock, the inertial equation in Burgers.
- D15's index work is complete and correct, including the dropped u_j,j term.

## point_vortex_lab (E8)
**Must fix** — none.

**Should fix**
1. On phones the Invariants view has no legend and only the "1" tick (`phone-tall__tour-step5.png`, `phone-tall__explain.png`). Put "P_y, H flat" in the view title when narrow, as the desktop title does.
2. The dotted "x₁ (not conserved)" ghost is clipped at the top of the y range and shows flat plateaus (`desktop__tour-step3.png`). Rescale it into the axis, or drop it when it would clip.

**What works**
- Click-to-place, drag and delete vortices, with singular configurations guarded.
- The invariants plotted as Q/Q(0) flat lines give a live check on the integrator.
- The measured orbit rate (◇) against (Γ₁ + Γ₂)/2πh² in the view title ("rate 0.637 · formula 0.637").
- The Fig. 5.11 correction is taught with a number (−1.70 m/s at G) and an orange arrow.

## vortex_sheet_rollup (E9)
**Must fix** — none.

**Should fix**
1. The L1 parity rows use rtol 1e-3, but the JS Simpson and Python `quad` agree to 5e-11 (N = 10) and 1.4e-9 (N = 100). Tighten them to 1e-8 for N = 10 and 100, and to 1e-4 for N = 1000 (Simpson with 2000 panels vs a spacing of 1 mm). A loose tolerance hides regressions.
2. D23 step 2's line "+v dn − v dn" reads as identically zero. Write it as $+v(s{+}ds)\,dn-v(s)\,dn$, so the *why* ("a remainder of size dn·ds") matches the line.

**What works**
- u(y) of the row against the continuous sheet, with the jump marked γ.
- A draggable box whose side bars trade exactly while the total stays γ ds, so "the vorticity inside" becomes visible.
- The L1 error shows slope −1 against N.
- The roll-up has the linear-theory ghost labelled as a Ch. 11 preview.
- The caption-convention toggle flips only the sign.

## Round-1 verdict (superseded): FAIL
- Six explainers pass: E1, E4, E5, E6, E8 and E9.
- **Three fail, each on one phone-legibility Must fix:**
  - **E2:** the minus sign hidden on the tornado z-axis at 390×844.
  - **E3:** overlapping heavy/light labels and Γ ticks in the baroclinic step at 390×844.
  - **E7:** the Terms view is illegible in tour step 3 at 390×844.
- All physics, parity (268/268), equation fidelity, derivations (23/23) and book-slip flags are correct. Each fix is local to one drawing routine.

## Round 2 — E2, E3, E7 (2026-09-23)

What I re-ran:
- `tools/viz_lint.py` on the three files: all `ok`.
- `tools/shot.py` on each, full run (8 sizes, every tab, step and pager page). I waited for all three runs; the audits are dated 17:45–17:49.
  - E2: PASS, 280 views, 34/34 selftest rows.
  - E3: PASS, 360 views, 29/29.
  - E7: PASS, 472 views, 31/31.
  - 0 failures; KaTeX rendered at every size.

| slug | round-1 Must fix | round-2 evidence | verdict |
|---|---|---|---|
| vortex_pressure_funnel (E2) | the rotated "z [m]" title hid the minus sign of "−100" | `phone-tall__tour-step6.png`: the title has its own strip and "−100" / "−200" read correctly. Walkthrough pages at 360×640 are now 1, 3, 4, 2, 4, 4, 2 (were 1, 3, 5, 2, 5, 5, 2). Also looked at `phone__explore.png`: clean. | **PASS** |
| kelvin_material_loop (E3) | heavy/light labels overlapped; Γ(t) ticks overlapped | `phone-tall__tour-step6.png`: "ρ₂ heavy" top-left and "ρ₁ light" top-right; Γ(t) shows only the 0 tick and "Γ 0.02335 m²/s" at the dot. Also looked at `phone-tall__tour-step5.png` (viscous): the same short-plot treatment is readable. | **PASS** |
| vorticity_equation_rotating (E7) | Terms view illegible in tour step 3 | `phone-tall__tour-step3.png`: with the sketch row hidden, all seven rows are readable (local 0.2 dashed, planetary 0.2 amber, the others 0). The D15 result page is 5 pp at 360×640 (was 25), 2 at 390×844 and 2 in the notebook frame. `phone-tall__derive-d2p16.png` shows the boxed (5.30) result with the Terms view in "zero terms hidden" mode. | **PASS** |

**Still open (Should fix, not blocking):**
- All round-1 Should fixes of E1, E4, E5, E6, E8 and E9.
- The round-1 Should fixes of E2, E3 and E7, except these, which the builders closed:
  - E2 walkthrough length;
  - E7 D15 result page.
- In `phone__explore.png` (360×640), the E2 isobar labels "+200 / +400 / +600 Pa" still stack.

## Verdict: PASS (all nine explainers PASS)
