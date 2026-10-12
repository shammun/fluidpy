# Chapter 14 verification — Aerodynamics                       2026-10-12 · code at commit `4823a8a` (tests, reference and this report uncommitted)

## Verdict: PASS — `tests/test_ch14.py`: 95 passed, 0 failed, 0 skipped (90.4 s with the slow script test; 94 in 36.7 s without it); 82 of 82 planted mutants killed; full suite: 1979 passed, 0 failed, 0 skipped (2446.9 s)

`tests/test_ch14.py`: **95 passed, 0 failed, 0 skipped in 90.4 s**. Full suite: **1979 passed, 0 failed, 0 skipped, 1 warning in 2446.9 s (40 min 47 s), slow tests included**. Mutation check: **82 of 82 planted wrong variants killed (81 in the first pass; see "Mutation check")**.

Every script ran (the nine `scripts/ch14_*.py` and `reference/ch14/make_refs.py` exit 0), every test passes, each of the
16 CORE rows has at least two independent evidence levels (one of V1 / V2 / V3 / V5 in each), every coded NOTE item and
every one of the 127 names of design Part C (C.0–C.3), the 14 methods of `ZhukhovskyFlow` and the 18 figure helpers of
C.4 is called by at least one test (a self-audit test enforces it), all 18 two- and three-star derivations are
re-derived with sympy independently of the chapter's own engine, the computable printed slips (#1, #2, #4, #8) fail
their independent checks as printed while the coded forms pass, and nothing is `unverified` without an Open item.
**No physics defect was found in `fluidpy/`.** Five documentation-level findings are listed under "Open items"; none
blocks. `wagner_function` **was confirmed first-hand** (two public sources read today; the primary 1938 report was not).

## Environment
python 3.11.5 · numpy 2.4.6 · scipy 1.17.1 · sympy 1.14.0 · pint 0.25.3 · matplotlib 3.11.2 (Agg) · pytest 9.1.1 · Windows 11.
`.venv/Scripts/python.exe -m pytest tests/test_ch14.py -q -p no:cacheprovider` → **95 passed, 0 failed, 0 skipped in 90.4 s** (with `-m "not slow"`:
94 passed, 1 deselected in 36.7 s; without the private `tests/book_values_ch14.json`: 91 passed, 3 skipped — checked in a scratch copy). `tools/check_public.py tests/test_ch14.py tests/ch14_verify_figures.py reference/ch14/SOURCES.md
reference/ch14/make_refs.py reference/ch14/benchmarks.json` → "public-repo check OK (5 files)".
Pages re-read as images by the verifier before any equation was tested: printed pp. 705 ((14.8)), 707 ((14.9)–(14.12)),
711 ((14.13) and the line before it), 712 (the line for ε, (14.14)), 714 ((14.18)–(14.24)), 715 ((14.25)–(14.27)).
Every coded equation agrees with its page.

## What the orchestrator asked about — rulings

| # | Item | Ruling | Measured |
|---|---|---|---|
| 1 | `circle_flow_zeta(center=0j)` must reproduce chapter 6 | **Confirmed.** With and without the keyword the arrays are identical bit for bit; both equal $w=U\big(\zeta+\frac{a^2}{\zeta}\big)+\frac{i\Gamma}{2\pi}\ln\frac{\zeta}{a}$ (6.68) to 1e-13 on 377 points; the surface speed equals $u_\theta=-2U\sin\theta-\frac{\Gamma}{2\pi a}$ (6.37) = (14.10) | all of chapter 6's tests pass in the full suite |
| 2 | `illustrative_inputs()` switched to the design's set | **Confirmed.** The values of header convention 11 are returned (b = 0.25 m, e = 0.08, β = 4°, α = 5°, U = 30 m/s; W = 9000 N, s = 12 m, A = 15 m², Λ = 9.6; demo wing Λ = 6); 32 of the design's worked numbers (Part F "Check" lines, C.6 item 4) are reproduced to every printed digit | Γ_K = 15.96 m²/s, exact C_L = 1.058, rule 0.987; loop area 1.054; lumped-vortex ratios 0.6463, 0.8156, 0.8792, 0.9388; 90 % at 6.2 chords |
| 3 | contour parameter t starts at −β (trailing edge at index 0) | **Confirmed** for `zhukhovsky_airfoil`, `zhukhovsky_surface`, `surface_point`: `t[0] = −β`, `z[0] = 2b` to 1e-14, counterclockwise, upper surface first; every consumer (`section_geometry(te_index=0)`, the force by surface pressure, the separation estimate) agrees with a route that does not use the index | pressure-route force = Blasius force, order 2.000 |
| 4 | `wagner_function` — first-hand source? | **Confirmed first-hand, as a cited fit (V5).** $\phi(s)=1-0.165e^{-0.0455s}-0.335e^{-0.3s}$, s in semi-chords, read in (a) Dawson & Brunton, arXiv:2104.15122v2, Eq. (1) and Table 1 (PDF text extracted), which also states φ(0) = 0.5, φ → 1 and an absolute error "less than 1 %"; (b) the BYU FLOW Lab `Aeroelasticity.jl` page (C₁ = 0.165, C₂ = 0.335, ε₁ = 0.0455, ε₂ = 0.3, exponent (u/b)t with b the semi-chord). **Not read:** the primary report (R. T. Jones, NACA TN 667, 1938, as Dawson & Brunton cite it; other pages cite NACA Rep. 681 / TN 682) | the coded function equals the cited coefficients to 1e-15; it differs from Garrick's independent fit $1-\frac{2}{4+s}$ (same paper, Eq. (2), error "2 %") by at most 0.0218 < 0.01 + 0.02 |
| 5 | `thin_airfoil` (not in the book) | **Verified independently, three ways.** (i) the Fourier integrals done exactly with sympy for a cubic (reflexed) camber line: A₀, A₁, A₂, C_L, α_{L=0}, C_{M,c/4} agree to 1e-12; (ii) an independent numerical method — N lumped vortices on the chord with tangency at the three-quarter points — converges to the same C_L and C_{M,c/4} at order 2.00 (3.6e-6 and 6.3e-6 at N = 400); (iii) the cited public results: slope 2π per radian, centre of pressure of a symmetric section at the quarter chord. Parabolic arc: α_{L=0} = −2h/c = −β of the Zhukhovsky arc | label `analytic (ours)`, now with evidence |
| 6 | `vortex_pair_energy_per_length` (not in the book) | **Verified independently.** Cores: $\int_0^\sigma\tfrac12\rho\big(\frac{\Gamma r}{2\pi\sigma^2}\big)^2 2\pi r\,dr=\frac{\rho\Gamma^2}{16\pi}$ each (sympy). Outside: $\tfrac12\rho\lvert u\rvert^2$ of the pair integrated over the plane outside the two discs by adaptive quadrature equals $\frac{\rho\Gamma^2}{2\pi}\ln\frac{s}{\sigma}$ to 1e-3 at σ/s = 0.02 (the formula is leading order in σ/s) | label `analytic (ours)`, leading order |
| 7 | qualitative and approximate functions | Tested for what is honestly testable and **kept as labelled**: see "Open items" 6 | — |

## Validation table

Evidence labels: `analytic` V1 · `symbolic` V2 · `converged` V3 · `conserved` V4 · `benchmark` V5 · `book-value` V6 (private) · limits/symmetry V7.

| Concept / Eq. | Tier | fluidpy target | Evidence | Numbers | Label | Tests |
|---|---|---|---|---|---|---|
| C01 the wing in numbers, $\Lambda\equiv s^2/A=s/\bar c$ (14.1), $\bar c=A/s$ (14.2); level flight $C_L=\frac{W}{\frac12\rho U^2A}$ | CORE | `aspect_ratio`, `mean_chord`, `planform`, `level_flight`, `stall_speed`, `compressibility_check`, `section_geometry`, `resolve_lift_drag`, `section_coefficient` | V1, V2 (D01 + dimensions), V3, V7 | three planforms against closed forms 1e-13; ellipse area by trapezoid: order 1.50; ½Ma² against the isentropic ratio: 6.1e-4 at Ma = 0.3, next term (2 − γ)/8 Ma⁴ | analytic, symbolic, converged | `test_aspect_ratio_V1_*`, `test_aspect_ratio_V2_*`, `test_planform_V3_*`, `test_level_flight_V1_*`, `test_compressibility_V1_*`, `test_section_geometry_V1_*`, `test_section_geometry_V3_*`, `test_resolve_lift_drag_V1_*` |
| C02 the Kutta condition, $\left.\frac{dw}{d\zeta}\right\rvert_{\zeta=b}=0$ | CORE | `kutta_circulation`, `kutta_circulation_root`, `stagnation_points`, `rear_stagnation_point` | V1, V2 (D09), V5, V7 | root-finder = closed form to 1e-12 on 18 (e, β, α) cases; dw/dζ(b) < 1e-13 U; both stagnation points are zeros of dw/dζ for 8 circulations (on and off the surface) and agree with chapter 6; three planted variants fail | analytic, symbolic, benchmark | `test_kutta_condition_V1_*`, `test_kutta_condition_V5_*`, `test_kutta_circulation_V2_derivation`, `test_stagnation_points_V1_*` |
| C03 the starting vortex, $\Gamma_{ABCD}=0=\Gamma_{\text{foil}}+\Gamma_{\text{starting}}$ | CORE | `impulsive_start_lumped_vortex`, `wagner_function` | V4, V1 (limits), V3, V5 | total circulation ≤ 1.8e-15 at every step (rigid and free wake, steps in speed and incidence); steady limit πUc sin α: deficit 1.3e-3 at 400 chords, decaying as 1/(2 chords); first value → 1/3; order in dt 1.00; against the cited fit: 0.6 % over 2–20 chords | conserved, analytic (limits), converged, benchmark; the transient is a **model** | `test_starting_vortex_V4_*`, `test_starting_vortex_V1_*`, `test_starting_vortex_V3_*`, `test_starting_vortex_V5_*` |
| C04 the Zhukhovsky airfoil, $z=\zeta+\frac{b^2}{\zeta}$ (14.3), $x=(R+b^2/R)\cos\theta$, $y=(R-b^2/R)\sin\theta$ (14.4), $c\simeq4b$, camber $\simeq\frac12\beta c$, $t_{\max}/c\simeq1.3e$ (14.8) | CORE | `zhukhovsky_circle`, `zhukhovsky_circle_from_centre`, `zhukhovsky_airfoil`, `zhukhovsky_shape`, `circular_arc`, `trailing_edge_angle`, `trailing_edge_wedge`, `map_point`, `circular_arc_sym`, `camber_from_beta` | V1, V2 (D03, D04, D05), V3, V7 | arc: every image point on $x^2+(y+2b\cot2\beta)^2=(2b\csc2\beta)^2$ to 1e-12; exact chord $b\,[3+2e+1/(1+2e)]$ to 1e-9; thickness against our own optimiser to 2e-7; `section_geometry` order 1.98 (thickness), 2.12 (camber); rule error first order in e (0.98) | analytic, symbolic, converged | `test_zhukhovsky_circle_V1_*`, `test_zhukhovsky_airfoil_V1_*`, `test_zhukhovsky_shape_V7_*`, `test_trailing_edge_V7_*`, `test_circular_arc_V2_derivation`, `test_thickness_rule_V2_derivation`, `test_cusp_V2_derivation` |
| C05 the flow round the airfoil, $\frac{dw}{dz}=\frac{dw}{d\zeta}\frac{\zeta^2}{\zeta^2-b^2}$ (14.9) | CORE | `zhukhovsky_flow` (all 14 methods), `flat_plate_flow`, `zhukhovsky_inverse`, `circle_flow_zeta(center=)` | V1, V2 (D06, D07, D08), V4, V7 | section is a streamline (ψ spread < 1e-12 Ua); far field $Ue^{-i\alpha}$; dw/dz equals centred differences of w(z) in two directions (Cauchy–Riemann) to 2e-6; circulation −Γ on three z-contours and on three pairs of partner contours to 1e-9; e = β = 0 equals the chapter-6 route to 1e-11 | analytic, symbolic, conserved | `test_zhukhovsky_flow_V1_*`, `test_zhukhovsky_flow_V4_*`, `test_zhukhovsky_flow_V7_*`, `test_zhukhovsky_inverse_V1_*`, `test_offcentre_circle_flow_V2_derivation`, `test_circulation_invariance_V2_derivation`, `test_circle_flow_zeta_V7_*`, `test_circle_flow_zeta_V1_*` |
| C06 lift, $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11), $C_L\simeq2\pi(\alpha+\beta)$ (14.12) | CORE | `zhukhovsky_force`, `zhukhovsky_lift_coefficient`, `trailing_edge_speed`, `edge_speed`, `zero_lift_angle`, `thin_airfoil`, `lift_per_span` | V1, V2 (D09), V3, V5, V7 | Blasius: L = ρUΓ and D = 0 to 1e-10 for 18 sections × 4 circulations, force normal to the stream (lift is not F_y); our own contour integral agrees to 1e-9; surface-pressure route converges to it at order 2.000 (3.7e-7 at n = 4000); flat plate $2\pi\sin\alpha$ exactly, centre of pressure at 0.25 c | analytic, symbolic, converged, benchmark | `test_zhukhovsky_lift_V1_*`, `test_zhukhovsky_lift_V3_*`, `test_zhukhovsky_lift_V5_*`, `test_zhukhovsky_lift_V7_*`, `test_trailing_edge_speed_V1_*`, `test_thin_airfoil_V1_*`, `test_thin_airfoil_V3_*` |
| C07 surface speed and pressure, $C_p=1-\lvert dw/dz\rvert^2/U^2$ (ours) | CORE | `zhukhovsky_surface`, `surface_point`, `surface_circulation`, `normal_force_coefficient`, `lift_from_cp`, `flat_plate_forces`, `edge_flow_suction`, `pressure_coefficient` | V1, V2 (D10), V3, V4 | node-by-node equality with the complex chain rule 1e-9; C_p = 1 at the front stagnation point 1e-12; loop area = Blasius F_y, order 2.000 (3.9e-7 at n = 4000); −∮u_s ds = Γ: polygon order 2.000 (1.7e-7), circle route 1e-12; leading-edge singularity $\lvert dw/dz\rvert^2r\to U^2c\sin^2\alpha$ to 1e-4 | analytic, symbolic, converged, conserved | `test_surface_pressure_V1_*`, `test_surface_pressure_V4_*`, `test_pressure_loop_V2_derivation`, `test_flat_plate_V1_*` |
| C08 bound vortex and trailing filaments, $d\Gamma_{\text{trail}}=-\frac{d\Gamma}{dy}dy$ | CORE | `trailing_sheet_strength`, `circulation_series`, `span_angle`, `span_position`, `elliptic_circulation`, `staircase_from_loading`, `rolled_up_pair`, `tip_vortex_estimate`, `wake_rollup` | V1, V4, V7 | −dΓ/dy against centred differences 1e-7; the starboard half sheds exactly Γ(0), the whole sheet nothing (1e-12); spacing b′ = twice the centroid of the half-sheet (1e-10); roll-up: impulse 3e-16, blob Hamiltonian 7.5e-8, centroid spacing 6e-5 from b′, sink rate within 2 % of Γ₀/(2πb′) | analytic, conserved; the spiral `qualitative` | `test_shed_filament_V1_*`, `test_tip_vortices_V1_*`, `test_wake_rollup_V4_*` |
| C09 downwash, $w(y_1)=\frac{1}{4\pi}\int_{-s/2}^{s/2}\frac{d\Gamma}{dy}\frac{dy}{(y_1-y)}$ (14.13) | CORE | `downwash`, `staircase_downwash`, `horseshoe_velocity`, `horseshoe_downwash`, `horseshoe_wing`, `ground_effect_factor`, `downwash_along_x`, `segment_speed`, `segment_induced_velocity` | V1, V2 (D12), V3, V4, V7 | parabolic loading against our closed form $\frac{2\Gamma_0}{\pi s^2}\big[s+y_1\ln\big\lvert\frac{s/2-y_1}{s/2+y_1}\big\rvert\big]$: 2e-13; ellipse: Γ₁/2s at 8 stations 1e-9; principal value by symmetric exclusion 1e-6; horseshoe against our own quadrature of (5.17), all three components 1e-8; field solenoidal and irrotational 1e-6; no flow through the ground plane 1e-13; staircase orders 1.52 / 1.47 / 2.02 | analytic, symbolic, converged | `test_downwash_V1_*`, `test_downwash_V3_*`, `test_downwash_integral_V2_derivation`, `test_horseshoe_V1_*`, `test_ground_effect_V1_*`, `test_downwash_along_x_V7_*` |
| C10 induced drag, $\alpha_e=\alpha-\frac wU$ (14.14), $L=\int\rho U\Gamma\,dy$, $D_i=\int\rho w\Gamma\,dy$ (14.15) | CORE | `induced_angle`, `effective_angle`, `strip_forces`, `lift_and_induced_drag`, `wake_kinetic_energy_per_length`, `vortex_pair_energy_per_length` | V1, V2 (D13), V3, V4, V7 | quadrature of (14.15) with w from the direct quadrature of (14.13) = series results to 1e-8; trapezoid order 2.00 (6.4e-7 at 1601 stations); wake energy = D_i to 2e-16, and by an independent log-kernel double integral (the Kirchhoff energy of the sheet) to 2.7e-4 at 1500 nodes, halving with the step; elliptic case = ½ (added mass of a plate) w∞² to 1e-12 | analytic, symbolic, converged, conserved | `test_induced_drag_V1_strip_*`, `test_induced_drag_V1_quadrature_*`, `test_strip_forces_V2_derivation`, `test_induced_drag_V4_*`, `test_vortex_pair_energy_V1_*` |
| C11 the sine series and its downwash, $\Gamma=\sum\Gamma_n\sin(n\gamma)$ (14.18), $w(y_1)=\frac{1}{2s}\sum n\Gamma_n\frac{\sin(n\gamma_1)}{\sin\gamma_1}$ (14.19, last member) | CORE | `downwash_series`, `glauert_integral` | V1, V2 (D14, D15), V7 | Glauert's integral three ways (closed form, singularity subtracted 1e-8, symmetric exclusion 2e-5) for n = 0…6 at six stations; recursion 1e-12; last member = principal-value quadrature with the corrected denominator; tips n and n(−1)^{n+1} | analytic, symbolic | `test_glauert_integral_V1_*`, `test_downwash_series_V1_*`, `test_series_in_downwash_integral_V2_derivation`, `test_glauert_integral_V2_derivation` |
| C12 the lifting-line equation, $\frac K2Uc(\alpha+\beta)=\sum\big(1+\frac{nKc}{4s\sin\gamma}\big)\Gamma_n\sin(n\gamma)$ (14.20) | CORE | `lifting_line_modes`, `lifting_line_matrix`, `lifting_line_rhs`, `lifting_line_residual`, `solve_lifting_line`, `section_lift_coefficient`, `elliptic_twist_rectangular`, `twist_for_elliptic_loading`, `wake_vortex_rolling_moment` | V1, V2 (D16), V3, V7 | matrix and right-hand side entry by entry 1e-13; elliptic planform: higher coefficients 7e-17 Γ₁ and C_L exact; rectangular wing Λ = 6: Cauchy order 3.77 (C_L), 3.77 (C_Di), residual of (14.17) between stations 0.38 → 1.1e-5 for N = 4 → 128; (14.17) with w by direct quadrature of (14.13): 2.2e-5 at N = 32; twisted wings: e_span = 1 to 1e-11 | analytic, symbolic, converged | `test_lifting_line_matrix_V1_*`, `test_lifting_line_V1_elliptic_*`, `test_lifting_line_V3_*`, `test_lifting_line_V1_twist_*`, `test_lifting_line_equation_V2_derivation`, `test_rolling_moment_V1_*` |
| C13 $L=\frac{\pi s}{4}\rho U\Gamma_1$ (14.21), $D_i=\frac\pi8\rho\sum n\Gamma_n^2$ (14.22) | CORE | `lift_from_series`, `induced_drag_from_series`, `span_efficiency` | V1 (quadrature, with C10), V2 (D17, D18), V4 (wake energy), V7 | 40 random loadings: e_span < 1, D_i = (elliptic D_i)/e_span to 1e-13, lift unchanged | analytic, symbolic, conserved | `test_series_forces_V7_*`, `test_lift_from_series_V2_derivation`, `test_induced_drag_from_series_V2_derivation`, `test_induced_drag_V1_quadrature_*` |
| C14 the elliptic loading, $\Gamma=\Gamma_1\sqrt{1-(2y/s)^2}$ (14.23), $w=\Gamma_1/2s$ (14.24), $D_i=\frac\pi8\rho\Gamma_1^2=\frac{2L^2}{\pi\rho U^2s^2}$ (14.25), $C_{D_i}=\frac{C_L^2}{\pi\Lambda}$ (14.26) | CORE | `elliptic_circulation`, `elliptic_downwash`, `elliptic_induced_drag`, `induced_drag_coefficient`, `finite_wing_lift_slope`, `elliptic_wing` | V1, V2 (D19, D20), V5, V7 | both members of (14.25) and all three of (14.26) 1e-13; equals the cited forms 1e-14; solver on nine elliptic wings reproduces K/(1 + K/πΛ) to 1e-12 | analytic, symbolic, benchmark | `test_elliptic_loading_V1_*`, `test_elliptic_optimum_V2_derivation`, `test_finite_wing_lift_slope_V2_derivation`, `test_lifting_line_V1_elliptic_*` |
| C15 the drag breakdown and the polar, $C_D=C_{D0}+\frac{C_L^2}{\pi\Lambda e_{\text{span}}}$ (ours) | CORE | `drag_polar`, `drag_breakdown`, `drag_share`, `best_glide`, `drag_vs_speed`, `speed_min_drag`, `stall_speed`, `profile_drag_estimate`, `generic_lift_curve`, `generic_drag_curve`, `stall_regime`, `chord_reynolds` | V1, V2 (D21 + dimensions), V7 | parts sum 1e-14; best glide against numerical minimisation 1e-6 (C_L) and 1e-10 (L/D); least-drag speed 34.48 m/s where the two parts are equal; induced share 0.64 at 30 m/s, 0.13 at 55 m/s (trap T9) | analytic, symbolic; the curves `qualitative`, the profile-drag estimate `approximate` | `test_drag_polar_V1_*`, `test_best_glide_V2_derivation`, `test_stall_models_V7_*`, `test_profile_drag_estimate_V7_*` |
| C16 thrust from lift, $\mathbf U_r=\mathbf U-\mathbf V$, drive $=L\sin\beta_a-D\cos\beta_a$ (ours) | CORE | `apparent_wind`, `sail_forces`, `sail_drive_curve`, `velocity_made_good`, `sailboat_course_limit`, `sailboat_equilibrium`, `oscillating_foil` | V1, V2 (D22), V4, V7 | drive and lateral force against our own vector construction 1e-12; zero drive exactly at tan β_a = D/L; steady boat speed against our own root 1e-9; force polygon closes 1e-9, β_a = ε_A + ε_H to 1e-10; fish tail: efficiency 1 for D = 0, 0 at tan θ = D/L | analytic, symbolic, conserved | `test_apparent_wind_V1_*`, `test_sailboat_V1_*`, `test_oscillating_foil_V1_*`, `test_sailing_V2_derivation` |
| NOTE: suction-side separation (N13) | NOTE | `suction_side_separation` | V7 only | moves forward monotonically 0.50 c → 0.003 c for α = 0 → 10°; same point on both sides of a symmetric section at α = 0 (1e-9); position independent of Re, momentum thickness ∝ Re^{−1/2} | **qualitative** | `test_suction_side_separation_V7_*` |
| NOTE: projectile with drag (N83) | NOTE | `projectile_range_drag` | V7, V3 | C_D → 0: vacuum range 1e-6; `solve_ivp` approaches the formula as θ₀ → 0 (7 % → 1.2 %, order 0.86) | analytic | `test_projectile_V7_*` |
| NOTE: tables | NOTE | `illustrative_inputs`, `conventions_table`, `book_slips`, `traps` | V7 | 10 slips with `kind` ("slip" for #1, #2, #5), 18 traps, 15-row letter table | — | `test_tables_V7_*`, `test_contract_V7_*` (4), `test_public_V7_*`, `test_design_numbers_V7_*` |

Evidence counts. Tests carrying each level (a test can carry two): **V1 41 · V2 21 · V3 11 · V4 12 · V5 5 · V6 3 · V7 28**
(by the level in the test name: V1 37 · V2 19 · V3 7 · V4 5 · V5 4 · V6 3 · V7 20; 95 tests, 18 of them
`…_V2_derivation`). CORE rows per validation label: **analytic 16 · symbolic 14** (all but C03 and C08, whose
derivations are one-star) **· converged 8** (C01, C03, C04, C06, C07, C09, C10, C12) **· conserved 7** (C03, C05, C07,
C08, C10, C13, C16) **· benchmark 4** (C02, C03, C06, C14) · book-value: 3 private tests · **qualitative: 5 functions**
(`suction_side_separation`, `generic_lift_curve`, `generic_drag_curve`, `stall_regime`, the spiral of `wake_rollup`)
· approximate 1 (`profile_drag_estimate`) · model 1 (the transient of `impulsive_start_lumped_vortex`) ·
**unverified 0**.

## Functions used by the notebook and explainers (design Part C)
`test_contract_V7_every_part_c_name_exists_and_is_exercised_here` holds the 127 names of C.0–C.3 (30 of
`core.airfoils`, 40 of `core.lifting_line`, 39 of `ch14`, 18 re-exported), asserts that each exists on `ch14` and appears as
`ch14.<name>` / `AF.<name>` / `LL.<name>` in the body of a test, and that each of the 14 methods of `ZhukhovskyFlow` is
called. The three wrappers added for explainer parity rows — `surface_point`, `edge_speed`, `map_point` — are tested
against `zhukhovsky_surface` at a node, `speed_near_edge` and `joukowski`, and at the trailing edge (the finite limit
with the Kutta circulation, `inf` otherwise). `test_contract_V7_signatures_match_design_part_c` compares 22 signatures
(argument order and keyword-only markers; K has no default anywhere — four calls without it raise `TypeError`);
`test_contract_V7_scalars_in_floats_out` checks 37 scalar calls and 6 dicts; `test_contract_V7_reexports_are_the_same_objects`
checks the 71 + 18 re-exports. The 18 figure helpers of `scripts/ch14_drawings.py` (C.4) each return a `Figure`
(`test_drawings_V7_*`). **Untested Part C functions: none.**

## Derivations (curation §4b, design Part F) — sympy, independent of the chapter's engine

| D | Result | Test | What is checked |
|---|---|---|---|
| D03 ★★★ | (14.4) → (14.5) → cosine rule → (14.6) → $x^2+(y+2b\cot2\beta)^2=(2b\csc2\beta)^2$, camber ≃ ½βc | `test_circular_arc_V2_derivation` | steps 1, 2, 4, 5, 7, 8, 9, 11, 12, 13, 14, 15; the cosine rule is shown to be the circle $\lvert\zeta-ib\tan\beta\rvert=b/\cos\beta$ written out; `circular_arc_sym` tested separately |
| D04 ★★★ | $a\equiv b(1+e)$ (14.7) → c ≃ 4b, $t_{\max}/c\simeq\frac{3\sqrt3}{4}e$ at x = −b | `test_thickness_rule_V2_derivation` | steps 1–11; a maximum (second derivative < 0); exact chord and its expansion 1 + e²; the cusp scaling y ~ beθ³, 2b − x ~ bθ² |
| D15 ★★★ | Glauert's integral; (14.19) last member | `test_glauert_integral_V2_derivation` | step 1 (antiderivative, both sides of γ₁), 2–3 (I₀ = 0), 4, 5, 6, 7, 8, 10, 11; direct evaluation for n = 1…6 through Chebyshev polynomials |
| D18 ★★★ | $D_i=\frac\pi8\rho\sum n\Gamma_n^2$ (14.22) | `test_induced_drag_from_series_V2_derivation` | steps 1–9 with four coefficients and two indices; the two traps (sin γ cancelled twice; factor n lost) give different results |
| D05, D06, D08, D09, D10, D12, D13, D14, D16, D17, D19, D20, D21, D22 ★★ | as in curation §4b | `test_cusp_`, `test_offcentre_circle_flow_`, `test_circulation_invariance_`, `test_kutta_circulation_`, `test_pressure_loop_`, `test_downwash_integral_`, `test_strip_forces_`, `test_series_in_downwash_integral_`, `test_lifting_line_equation_`, `test_lift_from_series_`, `test_elliptic_optimum_`, `test_finite_wing_lift_slope_`, `test_best_glide_`, `test_sailing_` `…_V2_derivation` | each result from its stated start and assumptions; D10 also by Gauss's theorem; D09 also by l'Hôpital for the finite edge speed $U\frac ba\cos(\alpha+\beta)e^{2i\beta}$; D07 (★) inside D06's test |

Every intermediate line of the 14 two-star derivations was read against the algebra as well. **No wrong line was
found in Part F.**

## Convergence studies

| Scheme | Grids | Observed order | Design order |
|---|---|---|---|
| Surface-pressure force and moment (`zhukhovsky_force(method="pressure")`) against the Blasius integral | n = 250 … 4000 | 2.000 (L), 2.000 (D), 2.000 (M_le); 3.7e-7, 3.7e-8, 3.9e-8 at n = 4000 | 2 |
| Loop area of C_p (`normal_force_coefficient`) against the Blasius F_y | n = 500 … 4000 | 2.000 | 2 |
| `surface_circulation` (polygon) | n = 500 … 4000 | 2.000; 6.8e-7 at n = 2000 (docstring: 1e-6) | 2 |
| `section_geometry` (thickness of the symmetric section; camber of the arc) | n = 250 … 4000 | 1.98; 2.12 | 2 |
| `trailing_edge_angle` of the cusp | n = 100 … 1600 | 1.00 (→ 0) | 1 |
| `planform("elliptic")` trapezoid area | n = 41 … 641 | 1.50 | 1.5 (square-root ends) |
| `impulsive_start_lumped_vortex`, time step | dt = 0.2 … 0.0125 c/U at 1 and 4 chords; free wake at 2 chords | 1.00, 1.01, 1.01 (rigid); 1.01, 1.01 (free) | 1 |
| `staircase_downwash` → (14.13), elliptic loading | 20 … 320 panels | 1.52 (uniform panels, mid-span), 1.47 (inner half), 2.02 (cosine panels, mid-span) | docstring: "first order" (Open item 3) |
| `lift_and_induced_drag` (trapezoid, stations uniform in γ) | 201 … 1601 | 2.00 (L), 2.00 (D_i) | 2 |
| `solve_lifting_line`, rectangular wing Λ = 6 (differences between N and 2N) | N = 4 … 128 | 3.77 (C_L), 3.77 (C_Di); 5.2e-8 between N = 64 and 128; condition number 1.8 → 605 | none stated |
| `lifting_line_residual` between stations | N = 4 … 128 | 0.38, 0.12, 1.0e-2, 2.1e-3, 1.4e-4, 1.1e-5 of Γ₁ | → 0 |
| Speed approaching the cusp, `edge_speed` → $U\frac ba\cos(\alpha+\beta)$ | d = 1e-3 … 1e-8 m | 0.497 (error ∝ √d); 2.2e-4 at 1e-6 | docstring: "1e-4 at 1e-6" (Open item 2) |
| Rule $t_{\max}/c\simeq\frac{3\sqrt3}4e$ against the contour | e = 0.1 … 0.005 | 0.98 (first-order rule; 9.3 % high at e = 0.1) | 1 |
| Discrete-vortex method against `thin_airfoil` (test-only) | N = 25 … 400 | 2.00 (C_L), 2.00 (C_M) | — |
| `projectile_range_drag`: numerical range → formula as θ₀ → 0 | 0.4 … 0.05 rad | 0.86 | trend only |

## Conservation / invariant residuals

| Quantity | Residual |
|---|---|
| Kelvin: Γ_bound + ΣΓ_wake in the lumped-vortex march (rigid, free, steps in U and α) | ≤ 1.8e-15 m²/s of Γ ≈ 16 m²/s |
| Circulation of the off-centre circle flow on three contours; flux through them | 1e-10 Γ; 1e-10 Γ |
| Circulation in the z-plane equals that in the ζ-plane (partner contours) | 1e-11; the velocities differ by more than 1e-3 U |
| Blasius force: L − ρUΓ, D | 1.1e-13 N/m, 3.4e-13 N/m of 587 N/m |
| Control-volume force round a plate edge, independent of the circle's radius (1e-5 … 5) | 1e-12 |
| Helmholtz: circulation shed by the starboard half − Γ(mid-span); by the whole sheet | 1e-12; 1e-12 |
| Wake kinetic energy per length − D_i of (14.22) | 2e-16 (relative) |
| Roll-up: impulse; blob Hamiltonian; centroid of each half | 3e-16; 7.5e-8; 2.5e-16 s |
| Horseshoe field: divergence, curl (three generic points) | < 1e-6 of the largest velocity gradient |
| Ground plane: normal velocity of horseshoe + image | < 1e-13 Γ/s |
| Sailing: aerodynamic + hydrodynamic force; β_a − ε_A − ε_H | 1e-9 of the force; 1e-16 rad |
| D_i U = power identity of the fish tail (D = 0): power out / power in | 1 to 1e-14 |

## Benchmarks used (all read first-hand on 2026-10-12; `reference/ch14/benchmarks.json`, `SOURCES.md`)

| Value / form | Ours | Source | Status |
|---|---|---|---|
| c_l = 2πα; cambered c_l = c_l0 + 2πα (slope 2π per radian) | slope of (14.12) = 2π to 1e-12; exact flat-plate slope at α = 1e-4 within 1e-6; `thin_airfoil` | Wikipedia, "Airfoil", section "Thin airfoil theory" | secondary source (form) |
| centre of pressure of a thin symmetric section: 0.25 c | 0.25 to 1e-10 (Blasius moment of the flat plate, three incidences) | same page | secondary (value) |
| Γ = 4πV∞R sin(α + sin⁻¹(μ_y/R)), R = √((1 − μ_x)² + μ_y²) | `kutta_circulation` and `kutta_circulation_root` to 1e-12 on 18 cases | Wikipedia, "Joukowsky transform" | secondary (form) |
| C_D,i = C_L²/(π AR e), e = 1.0 for an elliptic distribution | `induced_drag_coefficient` to 1e-14 | NASA Glenn, Beginner's Guide to Aeronautics, "Induced Drag Coefficient"; Wikipedia, "Lift-induced drag" | two sources (form) |
| D_i = L²/(½ρV²πb²) | `elliptic_induced_drag` to 1e-14 | Wikipedia, "Lift-induced drag" | secondary (form) |
| Wagner's function, Jones' fit: 0.165, 0.0455, 0.335, 0.3; φ(0) = 0.5; error < 1 % | coded coefficients identical; lumped-vortex model within 0.6 % of it over 2–20 chords | Dawson & Brunton, arXiv:2104.15122v2, Eq. (1), Table 1; BYU FLOW Lab, Aeroelasticity.jl documentation | two sources, first-hand; primary report not read |
| Garrick's fit 1 − 2/(4 + s), error 2 % | Jones' fit within 0.022 of it | Dawson & Brunton, Eq. (2) | cross-check |
| induced share of a wide-body's cruise drag ≈ 37 % | evidence for trap T9 only (ours: 13 % at cruise, 64 % at 30 m/s, 50 % at best glide) | Wikipedia, "Lift-induced drag" | not a test of our code |

Read and **not used**: NASA Glenn "Downwash Effects on Lift" writes Cl = Cl0/(1 + Cl0/(π AR)) with a lift *coefficient*
where lifting-line theory has the section *slope* K (the two agree only if Cl0 is read as the slope); our
$\frac{dC_L}{d\alpha}=\frac{K}{1+K/(\pi\Lambda)}$ is verified by sympy (D20) and by the solver instead. No readable table of
Glauert's rectangular-wing factor was found.

## Numbers from the text (private values redacted to relative differences)
Three V6 tests read `tests/book_values_ch14.json` (skipped when it is absent; none of its numbers is typed in the public
part of the test file — `test_public_V7_no_private_book_number_is_typed_in_this_file` and `tools/check_public.py` enforce it).

| What | Result |
|---|---|
| Thickness coefficient of (14.8) as printed against 3√3/4 | 0.07 % (the printed value is ours rounded to two digits) |
| Analyst's measured thickness ratios at three values of e, exact chord at e = 0.1, exact lift slope at e = 0.1 | 5e-5, 1e-6, 5e-5 (four- to seven-digit values) |
| Lift-slope per degree quoted in §14.7, converted to per radian; "most airfoils" slope of §14.6 | exact conversion; within 9 % and 5 % of 2π |
| Low-Mach limit of §14.1 | the classifier flips at that value |
| Exercise with two piecewise-linear pressure curves | loop area exact; C_L within 0.24 % of the stated answer (which is the area projected by cos α) |
| Flat-plate exercise: velocities at a named point in both planes, circulation | exact (1e-12) |
| Exercise that gives a, b, β of a cambered section | the (a, b, β) parametrisation reproduces a closed section with its edge at 2b; Kutta circulation from the same a |
| Lift coefficient at a second incidence | 1e-4 |
| Trailing-edge exponents of the wedge exercise | included angles 2π(1 − 1/n), stagnation point for each |
| Untwisted elliptic wing in level flight: induced drag, circulation, circulation half-way along a wing | 1e-4, 2e-4, 2e-4 |
| Rectangular wing with elliptic twist: flight angle, Γ₁, lift, induced drag; solver with that twist | 1e-4, 1e-5, 1e-5, 1e-4; higher coefficients < 1e-10 |
| Rolling moment in a wake vortex | 1e-5 |
| Core size and induced drag (a doubled core) | within 1.5 percentage points of the analyst's "about" value |
| Rule 9: our inputs against the exercise inputs and the quoted stall values | all different |

## Printed slips — confirmed or refuted

| # | Printed | Verdict | Test |
|---|---|---|---|
| 1 | (14.19), first two members: denominator cos γ₁ − cos γ | **Slip, confirmed on the page image (p. 714).** From (14.13) with y = −(s/2) cos γ, y₁ − y = (s/2)(cos γ − cos γ₁) (sympy); the printed middle member, evaluated by an independent principal-value quadrature, equals **minus** the last member (2e-5); `downwash_series(printed_middle=True)` returns exactly that; the direct quadrature of (14.13) and Biot–Savart (a horseshoe with Γ > 0 gives a downward velocity at its centre) both side with the last member. Everything from the last member of (14.19) on is right | `test_downwash_series_V1_*`, `test_series_in_downwash_integral_V2_derivation`, `test_glauert_integral_V1_*` |
| 2 | "ε = tan(w/U)" (p. 712) | **Slip, confirmed on the page.** The velocity triangle gives ε = tan⁻¹(w/U); the printed form is larger by ⅔(w/U)³ (5° at w/U = ½); `induced_angle(small=False)` is the arctangent | `test_induced_drag_V1_strip_*`, `test_strip_forces_V2_derivation` |
| 3 | "at location y (< y₁)" | **Loose, confirmed on the page (p. 711):** the caption of the figure on the same page says y > y₁; the formula and (14.13) are consistent with the caption (D12) | `test_downwash_integral_V2_derivation` |
| 4 | (14.8): "camber = ≅"; coefficient 1.3 | **Loose, confirmed (p. 705).** The first-order coefficient is 3√3/4 = 1.29904 (sympy, D04); the rule is first order in e: measured t_max/(c e) = 1.2734, 1.2365, 1.1785 at e = 0.02, 0.05, 0.1 (9.3 % low at 0.1); the planted radius b(1 + e) cos β puts ζ = −b outside the circle | `test_zhukhovsky_shape_V7_*`, `test_thickness_rule_V2_derivation`, `test_zhukhovsky_circle_V1_*` |
| 5 | Lanchester's year of birth | **Slip.** Born 23 October 1868 (Wikipedia, "Frederick W. Lanchester", read 2026-10-12); the printed dates are in the private file and differ | `test_book_V6_section_shape_rules_and_lift_slope` |
| 6 | "(4.107) and (4.108) in Section 4.3" | **Loose, confirmed on the page (p. 715)** — the coefficients are defined in §4.11 | — |
| 7 | upper limit N in (14.22), ∞ in (14.18)–(14.20) | **Loose, confirmed on the page (p. 714)**; the code sums whatever coefficients it is given | `test_series_forces_V7_*` |
| 8 | "chord length c (≃ 4a)", "camber/chord ratio β/2" | **Loose, confirmed (p. 707).** Exactly C_L = 8π(a/c) sin(α + β); the slope exceeds 2π by a factor 4(1 + e)/(3 + 2e + 1/(1 + 2e)) = 1.0909 at e = 0.1 (first order in e) | `test_zhukhovsky_lift_V7_*`, `test_zhukhovsky_lift_V1_*` |
| 9, 10 | span written b in two exercises; the "Λ-pattern" | not re-read by the verifier; nothing coded depends on them | — |

## Mutation check (do the tests discriminate?)
The package was copied to the scratchpad and one likely-wrong variant planted at a time (**82** in all, in
`core/airfoils.py`, `core/lifting_line.py`, `ch14_aerodynamics.py` and `core/conformal.py`): signs and factors of (14.4),
(14.9), (14.10), (14.11), (14.12), (14.13), (14.14), (14.16), (14.19), (14.20), (14.21), (14.22), (14.24), (14.25),
(14.26); Γ read counterclockwise; the circle's radius and centre; the inverse map's root chosen by |ζ|; lift resolved
with the wrong sign; loop area sign; Blasius moment sign; leading-edge suction; thin-airfoil A₁ and moment sign; the
half-line's direction; both horseshoe legs in the same sense; the image at the wrong depth; the far-wake factor; the
staircase's jump sign; roll-up spacing; wake energy with w for 2w; twist factors; Kelvin's condition dropped from the
lumped-vortex march; its control point moved; a digit swapped in Wagner's fit; drag polar, best glide and level-flight
factors; drive and lateral force swapped; apparent wind U + V; the course limit as a difference; the smoothstep
replaced by a ramp; the Mach limit; the stream turned the wrong way in `circle_flow_zeta`. **81 of 82 were killed in the first pass.** The survivor — the sign of the closing value in the
polygon rule of `surface_circulation` — changes the result only at second order (the last segment at the cusp has length
O(h²)); the test now asserts the docstring's own tolerance (1e-6 at n = 2000; measured 6.8e-7, the variant gives 4.8e-6)
and kills it. One variant (the Mach limit moved from 0.3) was killed only by a private V6 test; a public check on both
sides of the limit was added and kills it with the private file removed. **82 of 82 killed.** (Variants 0–61 were run
against the test file as it stood before the last additions — two derivation tests and four tightened tolerances —
variants 56–81 and the two re-runs against the later file; nothing was loosened in between.)

## Figures reproduced with our code (`tests/ch14_verify_figures.py` → `outputs/ch14/verify/`, local)

| Figure | Visual verdict |
|---|---|
| `fig_shapes.png` | Four circles through ζ = b map to a slit from −2 to 2, an arc that lies exactly on the grey circle with its top at 0.511 b, a symmetric section and a cambered one, all with a cusped tail at x = 2b; the measured thickness ratio leaves the line (3√3/4)e from below, the thickest station stays just behind the quarter chord, chord/4b rises slowly from 1, and the measured camber rises above β/2 as tan β/2. |
| `fig_kutta.png` | With Γ = 0 the rear stagnation point sits on the upper surface ahead of the edge and the flow turns the edge; with Γ_K it is at the edge and the streamlines leave smoothly; the edge speed is flat at U(b/a)cos(α + β) for Γ_K and rises as d^{−1/2} for 0.8 and 1.2 Γ_K; the lift curves cross zero at α = −β = −4°, the exact one steeper than 2π(α + β), the model curve drops at 14° to a plateau; the bound circulation rises from 1/3 toward 1, the shed circulation mirrors it, their sum is flat at zero, and Wagner's fit (from 1/2) is joined within the line width after about two chords. |
| `fig_pressure.png` | Suction peak at the nose growing from −0.6 to −7.5 between α = 0 and 10°, C_p = 1 at the front stagnation point, both curves meeting slightly above 0 at the trailing edge, loop areas 0.472, 1.054, 1.612; the signed surface speed changes sign once (front stagnation point) and is finite at both ends; the laminar separation point moves forward with α and jumps to the nose between 5.5° and 6° (trend only). |
| `fig_downwash.png` | One horseshoe: downwash Γ/(πs) at the centre growing without bound toward the tips, upwash outboard; staircases of 8, 32, 128 horseshoes lie on Γ₁/2s over a widening part of the span and the quadrature of (14.13) lies on it at every station; along the flight direction the downwash rises from 0 through w at the line to 2w, antisymmetric about the line. |
| `fig_lifting_line.png` | Rectangular, tapered and elliptic wings of the same span and area: the rectangular loading is the fullest toward the tips (e_span = 0.9539), the tapered one nearly elliptic (0.9902), the elliptic exact (1.0000); the differences between successive solutions fall along N⁻⁴; the polars are parabolas starting at C_D0 with the tangent from the origin touching at C_D = 2C_D0. |
| `fig_drag_and_sail.png` | Profile drag rising and induced drag falling cross at the minimum of the total, 34.5 m/s, above the stall speed; the sail's drive is negative below about 14° off the true wind, largest near 80°, while the lateral force is largest close to the wind; the model lift curve is the line up to 14°, then drops to 0.9, the drag curve has its minimum at zero lift and a step at stall. |
| `fig_wake.png` | The flat sheet curls up at both tips into two spirals while its middle sinks; impulse, Hamiltonian and centroid spacing stay on 1.000; the centroid of each half sinks along the rolled-up pair's line Γ₀t/(2πb′), five times slower than the flat sheet's own 2w. |

## Deviations & justifications
No `DEVIATION` marker appears in the three chapter modules. Choices that are ours and are stated in the docstrings: the
general circle $\mu=b-ae^{-i\beta}$, $a=b(1+e)/\cos\beta$ (reduces to both printed cases; verified); the exact lift
coefficient 8π(a/c) sin(α + β); the finite trailing-edge speed; every numerical method (collocation at γ_j = jπ/(N + 1),
singularity subtraction in `downwash`, the lumped-vortex march, vortex blobs in `wake_rollup`); the drag polar, best
glide, course limit and the thin-airfoil and Rankine-pair results (not in the book — each verified above). The corrected
forms of slips #1, #2, #4, #8 are implemented; the printed variant of #1 is kept as `downwash_series(printed_middle=True)`
and shown to have the wrong sign. Γ is clockwise everywhere (`Gamma_cw=` keyword); w is positive downward;
`horseshoe_velocity` returns uz upward and `horseshoe_downwash` its negative (both tested).

## Open items
1. **`wagner_function` docstring is stale** (`fluidpy/ch14_aerodynamics.py`): it still says "the verifier should confirm
   first-hand" and "verifier to confirm the source". Confirmed today in two public sources (ruling 4); the docstring
   should cite Dawson & Brunton, arXiv:2104.15122 (Eq. (1), Table 1) beside the BYU page and say that the primary report
   was not read. Label: `benchmark` (cited fit). Documentation only.
2. **`trailing_edge_speed` docstring understates the approach error**: "1e-4 at distance 1e-6"; measured for our section
   2.2e-4 at d = 1e-6 m, 7.0e-5 at 1e-7, falling as √d (order 0.497). Documentation only.
3. **`staircase_downwash` docstring says "first order"**; measured for the elliptic loading 1.5 (uniform panels) and 2.0
   (cosine panels) at mid-span. A conservative statement; documentation only.
4. **`section_geometry` is first order, not second, for a contour whose *rounded* trailing edge falls between two
   vertices** (the trailing edge is the vertex of largest x and is not refined; the leading edge is): an ellipse
   sampled from an offset start returns a camber of 1e-4 to 4e-4 of the semi-axis for 500–6000 points, exact when the tip
   is a vertex. Every section of this chapter has its sharp edge at index 0, where the measured order is 2. Not blocking;
   say so in the docstring or refine the trailing edge as the leading edge is.
5. **Area of a chord given as a callable with a kink** (`_area` in `core/lifting_line.py`, used by `solve_lifting_line`
   and `twist_for_elliptic_loading`): a 400-point Gauss rule across the kink of a tapered wing is off by 2.2e-6, so C_L,
   C_Di and Λ carry that error; exact with a `planform` dict (closed-form area) or a smooth chord. Not blocking; split
   the rule at mid-span or accept `A=`.
6. **Qualitative / approximate / model — tested for limits, invariants and monotonicity only, and they stay labelled so:**
   `suction_side_separation` (laminar Thwaites on the attached ideal pressure: a trend, not a maximum lift; the jump of
   the separation point to the nose near 6° is the model's, not a measurement); `generic_lift_curve`,
   `generic_drag_curve`, `stall_regime` (drawn models, every constant a required keyword); the spiral of `wake_rollup`
   past the first turn (its invariants are `conserved`); `profile_drag_estimate` (`approximate`: twice the plate
   friction times an input form factor; its turbulent coefficient is still the secondary-sourced one of chapter 9);
   `impulsive_start_lumped_vortex` (a model: first value 1/3 where Wagner's is 1/2; up to 0.17 below the fit before two
   chords, within 0.6 % after).
7. **Ours by convergence only** (no public table found): the rectangular wing's span efficiency 0.95393 and slope
   4.5304 per radian at Λ = 6, K = 2π (order 3.8, five digits by N = 64).
8. **`analytic (ours)` without a first-hand source, verified independently instead:** `thin_airfoil`,
   `vortex_pair_energy_per_length`, the rolled-up pair (spacing from the conserved impulse; sink rate confirmed by the
   blob run to 2 %), the sailing course limit, the leading-edge suction. All cited benchmarks are **secondary** sources
   (Wikipedia, NASA's Beginner's Guide, one arXiv paper, one documentation page): forms, not tables.
9. Slips #9 and #10 (notation in two exercises, the "Λ-pattern") were not re-read by the verifier.
10. Docstring "Validation:" lines of the three modules name evidence in words; the implementer may now replace them by
    the test names of this report (as chapters 6–13 do).

## Found wrong in design / analysis
Nothing material. The design's 32 worked numbers of Part F and C.6 that the tests assert are right, and so are more than 40
"expect" numbers of the Part A storyboard spot-checked by hand or by a one-off run (circle centre and radius, Mach
numbers, stall speeds, the tapered area, arc geometry, the rule-against-measured table, far-field velocity, edge speed,
suction peak −2.24 at 0.45 % chord, loop areas 0.472 / 1.054 / 1.612, force components and centre of pressure at
35.7 % chord, plate forces, tip-vortex numbers, horseshoe and ground-effect values, the N = 3 matrix, right-hand side and
solution). One rounding: the storyboard expects a thickness ratio of 0.0965 for our section where the code gives
0.09645 (prints 0.0964). The reused functions of C.0 that `ch14` does not re-export (`CONF.angle_preservation`,
`POT.power_branch`, `POT.cylinder`, `BS.wall_image_system`, `BL.thwaites`, `BL.separation_point`,
`BL.blasius_drag_coefficient`, `fluidpy.core.vortices.rankine_vortex`) exist with the signatures the design lists. D12 step 1 describes the half-line by the angles 0 and
90° (the leg seen from its far end) where the code's own comment uses 90° and 180° (seen from the wing) — the same
number, Γ/(4πd), tested both ways. The analysis's expected "3 % against the fit over 2–20 chords" is met with a wide
margin (0.6 %).

## Rest of the suite
Full suite, `.venv/Scripts/python.exe -m pytest -q -p no:cacheprovider` (slow tests included), run once after the last
change to `tests/test_ch14.py`: **1979 passed, 0 failed, 0 skipped, 1 warning in 2446.9 s (40 min 47 s)** — the 95 tests
of this chapter and 1884 of chapters 1–13 and the machinery. The one warning is an `IntegrationWarning` inside a
chapter-13 test (round-off in an adaptive quadrature at a 1e-13 tolerance), unrelated to this chapter.
`fluidpy/core/conformal.py` (`circle_flow_zeta` gained the keyword `center=0j`): the default path is unchanged — identical
arrays with and without the keyword, equal to (6.68) to 1e-13 — and all of chapter 6's tests pass. The switch of
`illustrative_inputs()` and of the contour parameter to t = −β broke nothing: all nine chapter scripts run in fast mode
(inside the slow test, 48 s) and in full mode (3–8 s each, 20 figures, run separately), and no other chapter's test
depends on them.
