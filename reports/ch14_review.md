# Chapter 14 — Aerodynamics: derivation review (phase 6)

Independent, read-only review by the `derivation-reviewer` (fresh context), written up by the orchestrator from its reply.

**Counts: Must fix 1 · Should fix 15 · Verified 24 groups.**
**Verdict: no wrong equation found in the three modules; one false comparison/label on CORE C03 (Wagner's lift function
presented as the growth of the bound circulation). With M1 fixed: PASS.**

Reviewer's own test run: `tests/test_ch14.py` — 95 passed, 0 failed, 0 skipped in 278.9 s (run alongside its own checks).
The full suite and the mutation check were not re-run by the reviewer.

Pages read as images: p729, p731, p732, p734, p738, p739, p740, p741, p742, p749, p753–p754.

---

## Must fix (1)

### M1 — Wagner's function is a lift ratio; code, script and design present it as the growth of the bound circulation

**Where**
- `fluidpy/ch14_aerodynamics.py:568-572` — docstring of `impulsive_start_lumped_vortex` ("the first value is 1/3 of the
  steady one where Wagner's exact thin-airfoil result is 1/2"; "V5 growth against `wagner_function`").
- `scripts/ch14_starting_vortex.py:5, 31, 36, 50-52` — Wagner's curve drawn on the axis "circulation / πUc sin α"; prints
  "90 % of the steady circulation after 6.2 chords; first step … (model limit 1/3; Wagner 1/2)".
- `analysis/ch14_design.md:589` (and curation line 636) — the notebook figure and explainer E3 overlay
  `wagner_function(2*s_chords)` on the bound / shed / sum circulation curves.
- `reports/ch14_verification.md` row C03 and open item 6 repeat the same sentence.

**What is wrong.** Wagner's φ(s) is circulatory lift / steady lift (the `wagner_function` docstring says so correctly).
The exact bound circulation is a different function: by Kelvin, $\Gamma_{\text{foil}}=-\Gamma_{\text{shed}}$, and nothing
has been shed at $t=0^+$, so it starts at 0, not ½. "Wagner: one half at once" on a circulation axis contradicts the
statement C03 teaches, $\Gamma_{\text{foil}}+\Gamma_{\text{starting}}=0$.

**Numbers (reviewer's own computation).** Exact thin-airfoil start-up with a piecewise-constant wake sheet; the two
relations used are recalled from the von Kármán–Sears formulation (not cited); the scheme reproduces Wagner's lift
function to within 0.4–1 % of the cited Jones fit, which is its check.

| chords travelled | exact circulation Γ/Γ∞ | exact lift φ | Jones fit | model `ratio` (Γ_bound/Γ_steady) |
|---|---|---|---|---|
| 0⁺ | 0 | 0.50 | 0.500 | 1/3 |
| 0.5 | 0.416 | 0.600 | 0.594 | 0.553 |
| 1 | 0.550 | 0.669 | 0.666 | 0.653 |
| 2 | 0.694 | 0.758 | 0.762 | 0.760 |
| 4 | 0.823 | 0.849 | 0.855 | 0.855 |
| 10 | 0.931 | 0.937 | 0.933 | 0.939 |

- Over 2–20 chords the model's ratio is within 0.012 of the exact lift, but up to 0.066 above the exact circulation
  (0.10 above at one chord).
- 90 % is reached at 6.4 chords by the exact lift, 7.1 chords by the exact circulation, 6.2 chords by the model.
- The model tracks the lift because its bound vortex is set by the normal velocity at the three-quarter-chord point,
  which fixes the quasi-steady circulatory lift $\rho U\Gamma_b$ — so the V5 comparison is legitimate once it is named a
  lift comparison.

**Correction (no number in the code changes).**
1. Docstring 568-572: "a model, not a simulation of the boundary layer. Its first value of Γ_bound is 1/3 of the steady
   one. In the exact thin-airfoil solution the circulation starts from zero (Kelvin: nothing has been shed yet) and it is
   the circulatory lift that starts at one half (Wagner). ρUΓ_bound of this model follows Wagner's lift function (about
   1 % after two chords) and runs above the exact circulation (0.65 against 0.55 after one chord)."
2. Same docstring, Validation: "V5 `lift_circulatory`/(ρUΓ_steady) against `wagner_function` (a lift ratio)".
3. Script: draw Wagner's curve against `lift_circulatory / (ρ U Γ_steady)`, labelled "Wagner's function — circulatory
   lift / steady lift (Jones' fit)"; line 36 becomes "90 % of the steady value after … chords (model); exact lift 6.4,
   exact circulation 7.1", or the Wagner clause is dropped.
4. Design 589 and E3: the ghost curve belongs on a lift axis, or carries the caption "lift, not circulation; the exact
   circulation starts at 0". Line 581 ("half the lift arrives almost at once") is right as written.
5. Verification report row C03 and open item 6: "circulatory lift ρUΓ_bound against the cited fit".
6. Optional: add the exact circulation curve as a second reference.

**Test gap.** `test_starting_vortex_V5_growth_against_the_cited_fits_of_wagners_function` (`tests/test_ch14.py:699-716`)
compares `ratio` with φ and passes, so it cannot catch this. Its tolerance is 0.03 where the measured deviation is 0.007
(dt = 0.05 c/U); tighten to 0.012 once the quantity is named.

---

## Should fix (15)

**The ten open items of the verification report — each confirmed, with the fix**

1. `wagner_function` docstring stale (`ch14_aerodynamics.py:630-635`): Source → "S. T. M. Dawson & S. L. Brunton,
   arXiv:2104.15122v2 (2021), Eq. (1) and Table 1 (coefficients, φ(0) = ½, absolute error below 1 %); the same
   coefficients on the BYU FLOW Lab Aeroelasticity.jl page 'Wagner's Function'. Both read first-hand on 2026-10-12; the
   primary report (R. T. Jones, NACA, 1938) was not read." Validation → "V5 coefficients equal the cited ones (1e-15);
   within 0.022 of Garrick's fit 1 − 2/(4 + s); V7 φ(0) = 0.5, φ(∞) = 1, monotone. Label: benchmark (cited fit)." Add
   "φ is a lift ratio, not a circulation ratio".
2. `trailing_edge_speed` approach error (`core/airfoils.py:529`): measured relative differences 2.2e-3, 2.2e-4, 7.0e-5,
   2.2e-5 at d = 1e-4, 1e-6, 1e-7, 1e-8 m → "V1 limit of |dw/dz| along the cusp's bisector; the difference falls as √d
   (2.2e-4 of the limit at d = 1e-6 m, 2.2e-5 at 1e-8 m for the chapter's section)."
3. `staircase_downwash` "first order" (`core/lifting_line.py:801`): measured orders 1.54→1.50 (uniform panels),
   2.05→2.01 (cosine panels) → "V3 → (14.24) at panel mid-points for the elliptic loading: order 1.5 with uniform panels,
   2.0 with cosine panels (20–320 panels)"; also the comment at `tests/test_ch14.py:1580`.
4. `section_geometry` (`core/airfoils.py:181-182`): for a rounded trailing edge between vertices the camber is first
   order (ellipse 2 × 0.3 from an offset start: −2.3e-3, −5.7e-4, −1.9e-4 at n = 500, 2000, 6000). Add: "Second order
   only when the trailing edge is a contour vertex (true for every Zhukhovsky contour, index 0); for a rounded trailing
   edge that falls between two vertices the camber is first order — pass `te_index` or refine."
5. `_area` across a kink (`core/lifting_line.py:428-430`): tapered wing 14.998728 against 14.9988 (−4.8e-6); C_L differs
   by the same relative amount between the dict and callable routes. Split the Gauss rule at mid-span, or accept `A=`.
6. Qualitative / approximate / model labels — present in `suction_side_separation`, `generic_lift_curve`,
   `generic_drag_curve`, `stall_regime`, `wake_rollup`, `profile_drag_estimate`, `impulsive_start_lumped_vortex`. No
   change beyond M1.
7. Rectangular wing e_span = 0.9539, slope 4.5304 at Λ = 6, K = 2π — reproduced by the reviewer's from-scratch
   collocation (identical to 1e-15, N = 41). Keep "ours by convergence".
8. `analytic (ours)` items (leading-edge suction, rolled-up pair, course limit, pair energy) — re-derived by hand, right.
9. Slips #9 and #10 — #9 confirmed on p753–p754; #10 seen only in the extracted text (page image not rendered).
10. Validation lines in words — append the test names from the verification table to each "Validation:" line, as
    chapters 6–13 do.

**New items**

11. Sibling argument order (teaching risk): `aspect_ratio(s, A)` against `mean_chord(A, s)`
    (`core/lifting_line.py:63, 73`) — `mean_chord(12, 15)` returns a plausible 0.8 m where the answer is 1.25 m. Same
    hazard in `level_flight` / `stall_speed` / `drag_vs_speed` and `flat_plate_forces` / `zhukhovsky_force`. Signatures
    are pinned by design Part C, so notebook and explainer code call these **by keyword**.
12. `wake_kinetic_energy_per_length` label (`core/lifting_line.py:841`): it evaluates ρ∫Γw dy, i.e. (14.15) again, so by
    itself it is V1; name the verifier's log-kernel double-integral test in the Validation line.
13. `downwash` label "analytic, converged" (`core/lifting_line.py:234`) — no order is asserted; make it "analytic".
14. Docstring contract gaps (rule 5): `zhukhovsky_inverse` (no Assumptions); `ZhukhovskyFlow` class and methods (no Book
    / Label lines); `lifting_line_modes`, `lifting_line_rhs`, `horseshoe_downwash`, `staircase_from_loading` (no
    Validation / Assumptions); `chord_reynolds`, `sail_drive_curve`, `velocity_made_good` (no Validation).
15. Clarity: `zhukhovsky_lift_coefficient(…, exact=False)` silently ignores `e`; `lift_circulatory` is the quasi-steady
    part only (the ρ·dΓ/dt term is omitted); `zhukhovsky_force["xcp_over_c"]` is measured along the x-axis from the
    leading edge's x; `lifting_line_matrix` accepts a chord array but `solve_lifting_line` then fails in `_area`;
    `book_slips()["10"]` should say the page image was not re-read.

---

## Verified (against the page images and by independent calculation)

**Sections and the Kutta condition**
1. (14.3)–(14.4), p731: `zhukhovsky_airfoil`, `map_point` — symbol by symbol.
2. (14.5)–(14.6) and the arc circle, p731: `circular_arc`.
3. (14.7), (14.8), p731–p732: the circle, both printed cases; c ≃ 4b, camber ≃ ½βc, 3√3/4 for the printed 1.3.
4. (14.9), p734: `dwdz_from_zeta`; dw/dz = u − iv; the l'Hôpital limit at the edge re-derived by hand.
5. (14.10)–(14.11), p734: θ from the stream-aligned diameter through Q; Γ clockwise; radius a, not b; stagnation roots.
6. (14.12), p734, and the exact 8π(a/c) sin(α+β).

**Independent numerical checks of the section flow** (b = 0.25 m, e = 0.08, β = 4°, α = 5°, U = 30 m/s, ρ = 1.2)
7. Γ_K = 15.962 m²/s; counterclockwise circulation on a circle of radius 3 m = −15.96197; far field correct.
8. Surface-pressure integral: L = 574.63097 N/m against ρUΓ = 574.63098, D = 8.5e-7; normal velocity 1.3e-8.
9. Signs for β < 0, α < 0 and β = 0.3; D = 0 throughout.
10. `zhukhovsky_inverse` stays continuous across the slit under a cambered section.
11. Force conventions, M_le sign, leading-edge suction re-derived.
12. `thin_airfoil` coefficients — the standard forms.

**Lifting line**
13. (14.13), p738: sign and principal value; direct quadrature agrees with `downwash` and `downwash_series` to 3e-5.
14. Slip #1, p741: printed denominator confirmed wrong-signed for y = −(s/2) cos γ; last member of (14.19) right.
15. Slip #2, p739: "ε = tan(w/U)" confirmed; (14.14) right.
16. (14.15)–(14.17), p740, and (14.20)–(14.22), p741: re-derived by hand.
17. (14.23)–(14.27), p741–p742; K/(1 + K/πΛ) re-derived.
18. Solver: own collocation gives C_L = 0.3953538, C_Di = 0.00869263, e_span = 0.9539366 (Λ = 6, K = 2π, α = 5°) —
    identical.
19. Horseshoe: three components against brute-force Biot–Savart to 1e-7; centre downwash Γ/(πs).
20. Far-wake factor 2: `downwash_along_x` 7.5e-5, 0.0401, 0.0833, 0.1266, 0.1666 at x = −50, −1, 0, 1, 50.
21. Ground effect, rolling moment, rolled-up pair spacing and pair energy re-derived.

**Other**
22. Lumped-vortex start: both equations per step, steady limit πUc sin α, first value 1/3 — numbers right, only their
    interpretation not (M1).
23. Projectile formula, p749; vacuum limit.
24. Drag polar, best glide, drag against speed; apparent wind, drive, lateral force, polygon closure β_a = ε_A + ε_H.

**Tests.** Discriminating where probed (rolling moment, along-x downwash, staircase each carry an independent route).
The only test that shares the code's misreading is the Wagner one (M1).

---

## Resolution (2026-10-12)

**Verdict after fixes: PASS.** `tests/test_ch14.py`: 95 passed, 0 failed, 0 skipped (79 s, verifier's run after the last
test edit); 16/16 planted mutants killed for the three changed tests. The full suite is run at the merge gate.

- **M1 — fixed.** `impulsive_start_lumped_vortex` and `wagner_function` docstrings say Wagner's function is a lift
  ratio; `scripts/ch14_starting_vortex.py` draws it on a separate lift panel (figure now 1×4 panels) and the circulation
  panel carries only the steady Kutta value; design lines 120, 315, 589 and curation §10.6 corrected; the notebook-builder
  and the `starting_vortex` explainer builder were told before they built the overlay (the explainer shows Wagner only in
  its lift view). The V5 test now compares `lift_circulatory`/(ρUΓ_steady) with the fit at 0.012 (was 0.03; measured
  0.0060) and is renamed `test_starting_vortex_V5_circulatory_lift_against_the_cited_fits_of_wagners_function`.
  Limit, stated honestly: for a steady stream the model's lift ratio and circulation ratio are the same number, so no
  code mutant separates the two — M1 was a mislabelled comparison, not a wrong number. The exact-circulation values
  (0 at the start, 0.55 after one chord, 90 % at 7.1 chords) are this review's computation, quoted and attributed, not
  re-derived or tested elsewhere.
- **90 % crossing: 6.1 chords, not 6.2**, at dt = 0.05 c/U (6.12 by interpolation; 6.2 belongs to dt = 0.1 c/U; the
  crossing is first order in dt and tends to about 6.04). Test tightened to ±0.01 about 6.13 (α = 5°).
- **Should fix 1–5, 10–15 — done** in the three modules (docstrings; `_area` now uses one Gauss rule per half-span, so a
  tapered callable chord's area is exact: 5.10000000000001 against 5.1, was −2.2e-6 relative; the test that pinned the
  old defect now asserts exactness at 1e-12). Items 6–9 needed no change. Argument orders were not changed (item 11):
  the docstrings say "call by keyword", and the notebook and explainer briefs require keyword calls.
- **Staircase order** test now asserts the measured orders 1.5 (uniform) and 2.0 (cosine) within ±0.05.
- `tests/ch14_verify_figures.py`: the Wagner curve on the verification figure is relabelled as a lift ratio.
