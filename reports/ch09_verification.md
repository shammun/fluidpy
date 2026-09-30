# Chapter 9 verification — Boundary Layers and Related Topics
2026-09-30 (re-verification after fixes F1-F4) · repository state `e5beebd` + uncommitted implementer fixes (`fluidpy/core/{boundary_layer,bluff_body,jets}.py`, `scripts/ch09_drawings.py`) + uncommitted verifier files (`tests/test_ch09.py`, `tests/ch09_verify_figures.py`, `reference/ch09/`) · verifier: `math-verifier`

## Verdict: PASS

All five earlier failures (F1, F2, F3a, F3b, F4) are fixed in `fluidpy/`; the previously failing tests now pass with **unchanged assertions and tolerances**. Provenance, stated honestly: `tests/test_ch09.py` is untracked (never committed), so `git diff HEAD` cannot show it. The only edits made to the five formerly failing tests in this round replace their "FAILS on the current code" comment lines by "(F#, fixed)" comments (checked by reading the file); every `assert` and tolerance (order ±0.25, 1e-3, 1e-7, 1e-6, 1e-8, the 8 % bound) is as in the FAIL report. Nine new tests were added for the repaired behaviours (below).

Suite results (2026-09-30): **ch01-ch09 test files: 1109 passed, 0 failed (10 min 24 s); `tests/test_machinery.py`: 9 passed (19 s) => 1118 passed, 0 failed.** ch09 alone: 134 collected items (128 test functions plus parametrisation), all pass. A single-command `pytest -q` over the whole tree was run three times: run 1 (before the nine new tests were appended) 1109 passed in 12 min 08 s; runs 2 and 3 (with the new tests) were killed by the harness time limit at about 96 %, during the last stage, i.e. the same browser-based machinery stall already noted in the FAIL report (that file alone passes in 19 s when nothing else runs). The final count is therefore given as the two pieces above. ch01-ch08 and machinery show no regression from the fixes.

Verdict rules: all 22 `ch09_*.py` demo scripts run headless (`ch09_drawings.py` is a helper, imported and smoke-tested) · every computable CORE row has >= 2 independent levels · every coded NOTE row >= 1 · every design Part C function (including C.5 drawings) is at least smoke-tested · 22/22 derivations re-derived in sympy · 25 planted wrong variants killed · nothing `unverified` without an Open item.

### Re-check of the five earlier failures
| # | test (unchanged assertion) | earlier | now |
|---|---|---|---|
| F1 | `test_falkner_skan_shoot_V3_hiemenz_by_shooting` | ValueError (no sign change) for m = 1, 4 | passes: Hiemenz 1.2325876568 and m = 4: 2.4057248594 to 1e-7; new test: shoot vs bvp agree to 1e-7 for m = -0.05, 0, 0.3, 0.5, 1, 2, 4 |
| F2 | `test_thwaites_bug_V1_thwaites_named_wedge_fails_for_n_above_1p4` | ValueError for every n >= 1.5 | passes: theta error -7.57 % at n = 4 (-0.0757 +- 1e-3); new test: n = 1.5, 2, 3, 4 give lambda = 0.45n/(5n+1) to 1e-6 and theta within 8 % of exact Falkner-Skan; the guard still rejects theta0 != 0 at a true 0/0 start |
| F3a | `test_march_V3_falkner_skan_tau0_is_second_order_in_dsigma` (n = 0.5, ny = 100...800) | order **0.98** | order **2.24** (pairwise 2.14, 2.21, 2.38; errors 1.36e-4, 3.07e-5, 6.6e-6, 1.27e-6); plate n = 0: 2.03 (1.82, 1.98, 2.29) |
| F3b | `test_march_V5_falkner_skan_tau0_within_1e_3_..._at_default_settings` | 0.48 %, 0.71 %, 0.88 % (n = 0.2, 0.5, 1) | 1.7e-4, 6.0e-6, 3.7e-5 (< 1e-3) |
| F4 | `test_wall_jet_ode_V3_design_expect_row_at_the_default_eta_max` | f_inf = 0.99961, err 5.5e-4 | f_inf = 1 - 1.2e-12, err_vs_9_83 = 6.5e-12; new test: f_inf^3 = 72 f''(0) to 1e-8 at f''(0) = 0.005, 1/72, 0.2, 3 with default arguments |

### New tests (9), each names the fixed item
`test_falkner_skan_shoot_V3_agrees_with_bvp_across_the_documented_range` · `test_thwaites_named_wedge_V1_power_law_lambda_and_theta_for_every_n_above_1p4` · `test_thwaites_V7_stagnation_guard_still_rejects_a_wrong_theta0` · `test_wall_jet_ode_V7_default_eta_max_scales_with_the_free_scale` · `test_march_V3_inlet_station_wall_shear_falls_with_refinement` (inlet tau0 error 2.29e-4, 2.87e-5, 3.6e-6, 4.5e-7 for ny = 100...800: monotone, order about 3; it used to grow 2.8e-3 to 8.7e-3) · `test_march_V7_below_the_fold_inlet_raises_a_clear_error_and_accepts_a_supplied_profile` (ValueError naming the separation value; a supplied Blasius profile runs) · `test_cylinder_cd_schematic_V1_is_continuous_at_re_1_and_equals_lamb_there` (C_D(1) = 8 pi/2.002 = 12.5538; jump across Re = 1 < 1e-6 relative; no step on a 400-point log sweep) · `test_drawings_V1_ch09_drawings_imports_and_every_helper_returns_a_figure` · `test_falkner_skan_V5_default_eta_max_blasius_within_1e_10_of_toepfer` (default eta_max now 12 for m = 0).

Counts per validation label (a test counts once for each label in its definition comment; 128 test functions, 134 collected items): **V1 37 · V2 34 · V3 17 · V4 9 · V5 16 · V6 5 · V7 24** (V2 includes the 22 derivation tests). `qualitative` items: C11's regime table, angles and drag-crisis magnitudes (Open items). `unverified`: none.

Observation (not a failure): wedge n = -0.05 wall shear at fixed 161 x-stations does not converge in Delta-sigma below 3e-4 because the x-discretisation and the slowly decaying inlet memory then dominate (errors 3.5e-5, 2.2e-4, 2.7e-4, 2.9e-4 for ny = 100...800); the order claim is asserted only where the sigma error dominates (n = 0.5 and the plate).


## Environment
Python 3.11.5 · numpy 2.4.6 · scipy 1.17.1 · sympy 1.14.0 · pint 0.25.3 · matplotlib 3.11.2 (Windows 11, `.venv`; `progress.json → environment` checked 2026-09-12 — unchanged). Commands: `.venv/Scripts/python.exe -m pytest tests/test_ch09.py -q -p no:cacheprovider` (non-slow part ≈ 90 s; the 3 `slow` tests — 23 scripts, two reversed-branch continuations — add ≈ 2 min); figures: `.venv/Scripts/python.exe tests/ch09_verify_figures.py`; references: `reference/ch09/make_refs.py`.
Public-repo rule: `tests/book_values_ch09.json` (git-ignored) was **not valid JSON** (`2.0/3.0` and `1.0/3.0` literals on the "sec_9_10_free_jet" line, so `json.load` raised at line 75); the verifier replaced them by decimal numbers (file is in `tests/`, git-ignored, no content change). V6 tests read only that file; this report quotes relative errors only.

## Validation table (CORE rows; every row has ≥ 2 independent levels)
| Concept / Eq. | Tier | fluidpy target | Evidence (V-levels) | Numbers (error, order, residual) | Label | Notes |
|---|---|---|---|---|---|---|
| C01 scaling, BL equations (9.4)–(9.11) | CORE | `boundary_layer_scales`, `to/from_bl_variables`, `bl_x_momentum_residual`, `bl_pressure_variation`, `outer_flow`, `ch09.bl_nondim_sympy`, `march_boundary_layer` | V1 (δ̄/L = Re^{−1/2} to 1e-13, adv = visc, (9.6) round trip) · V2 (coefficient table exactly Re⁰, Re⁻¹, Re⁻², D01 all 14 lines) · V3 ((9.9) residual on the Blasius field: order 1.77 over n = 20, 40, 80; continuity 1.78) · V7 (pressure variation ∝ Re^{−1.00}, ratio 7.6e-3 → 7.6e-5 for Re 1e3 → 1e5) | order 1.77 (design 2, span 1.1 decades ⇒ ±0.25) | analytic, symbolic, converged | R1 printed (9.7) fails the dimension check and the sympy residual |
| C02 δ*, θ, δ₉₉ (9.16)–(9.17) | CORE | `thicknesses`, `displacement/momentum_thickness`, `delta_level`, `profile_shape` | V1 (1 − e^{−y/a}: δ* = a, θ = a/2 to 1e-8; sine, cubic, linear, power closed forms via sympy to 1e-13) · V4 (ρU²θ = ∫τ₀dx to 1e-6 at 60 stations; v∞ = U dδ*/dx to 1e-6) · V2 (D03, D04) · V7 (H limits) | — | analytic, conserved | |
| C03 Blasius reduction (9.19)–(9.27) | CORE | `ch09.similarity_reduce_sympy("blasius")`, `blasius_fields`, `similarity_collapse_error` | V2 (residual 0, ηf′f″ terms cancel, brackets U²δ′/δ and νU/δ², D05 all 13 lines) · V1 (data collapse < 1e-12 at n = 0, m = ½; > 1e-2 for wrong exponents) · V4 (u = ψ_y, v = −ψ_x to 1e-6/1e-5) · V3 ((9.18) residual order 1.77) | collapse 1e-12 | symbolic, analytic | |
| C04 Blasius numbers (9.28)–(9.33) | CORE | `blasius_constants`, `falkner_skan(0.0,…)`, `blasius_delta99/…/drag` | V5 (f″(0): 5.7e-15 vs Wikipedia; Belden et al. κ√½: 7.7e-14) · V3 (Töpfer IVP vs `solve_bvp` η_max = 16: 2e-14; vs brentq shoot 1e-9; truncation 8 → 10 → 12: 1.9e-6, 9.9e-10, 7e-12) · V2 (ODE residual < 2e-6·max f″, BCs) · V7 (θ = 2f″(0) to 1e-10, H, v∞ = δ*/2) · V6 (0.02–0.4 %) | see Convergence | benchmark, converged | R2 (4.93 vs root 4.910, 0.41 %), R5 (sides = 2 doubles) planted variants fail |
| C05 Falkner–Skan (9.34)–(9.36) | CORE | `falkner_skan`, `falkner_skan_state`, `falkner_skan_separation`, `falkner_skan_fields`, `falkner_skan_table` | V5 (Belden κ at β = 0.5, 0, −0.12: ≤ 1.5e-12; separation β: 2.4e-10; reversed branch β = −0.12, −0.02: 5.5e-12, 6.4e-13; Hiemenz 1.232588: 2.8e-7; δ*/δ 0.6479: 7.3e-7) · V2 (residual of (9.36) < 3e-6, D07 all 12 lines with symbolic n) · V1 (momentum integral identity f″(0) = ((3m+1)/2)I_θ + mI_δ to 1e-8 for 7 members) · V3 (bvp vs shoot 5e-9) · V7 (f‴(0) = −m) | fold m = −0.090429 | benchmark, symbolic, analytic | F1 fixed: shoot ≡ bvp to 1e-7 for m ∈ [−0.05, 4] |
| C06 momentum integral (9.43) | CORE | `momentum_integral_residual`, `karman_pohlhausen`, `ch09.momentum_integral_sympy` | V1/V4 ((9.43) closes to 1e-7 for Blasius and 3 FS members; marched flat and diffuser layers close it to 7e-4 and 1.3e-3 of τ₀) · V2 (D08: exact sympy integrals on a manufactured flow, τ₀/ρ = d(U²θ)/dx + Uδ*U′ with residual 0) · V5 (Pohlhausen cubic/sine/quartic: δ, θ, C_f√Re exact to 1e-5; θ within 6 % of Blasius) | — | analytic, conserved | |
| C07 Thwaites (9.44)–(9.50) | CORE | `thwaites`, `thwaites_named`, `thwaites_l/H/L`, `thwaites_closure_table`, `thwaites_cylinder*`, `example_9_1/2` | V1 (flat plate θ² = 0.45νx/U to 1e-12; diffuser λ = −(0.45/4)[(1+x/L)⁴−1] to 1e-12; power law λ = 0.45n/(5n+1) to 1e-9; stagnation limit 0.075 exact; x_sep = 1.8^{1/4}−1) · V5 (θ error vs exact FS: +4.4 % (n = −0.09), +3.1 (−0.05), +1.0 (0), −1.6 (0.1), −4.2 (1/3), −6.3 (1), −7.6 (4): **all ≤ 8 %**) · V2 (D09, D10, D11; `thwaites_sympy`) · V6 · V7 | see Numbers from the text | analytic, benchmark, book-value | F2 fixed; τ₀ errors: ≤ 2.4 % for n ≥ 0 (FS closure), 5.5 % at n = −0.05, **−56 % at n = −0.085** (near the fold) |
| C08 separation, (9.51)–(9.52) | CORE | `wall_curvature`, `profile_inflection`, `separation_point`, `march_boundary_layer` | V1 (μ u_yy = dp/dx from the solved field, 2e-3, three exponents) · V7 (favourable none, Blasius at the wall, adverse finite η = 1.650 at n = −0.05, moving out towards the fold) · V5/V3 (PDE marching: n = −0.085, −0.089 stay attached to x = 50; n = −0.095, −0.1, −0.12 separate at x = 3.43, 0.71, 0.27 (Blasius inlet at x₀ = 0.1; after fix F3 the n = −0.095 position moved from 2.03 to 3.43 — near the fold the separation point is very sensitive to the inlet profile; the fold n = −0.09043 is still bracketed) — the Belden fold n = −0.09043 seen by an independent solver) | — | analytic, converged | |
| C09 form drag | CORE | `separated_pressure_drag`, `separated_cp`, `pressure_drag_from_cp`, `cp_ideal_cylinder` | V1 (closed form C_D,p = sin φ_s(1 − (4/3) sin²φ_s − C_b) vs Gauss–Legendre: 1e-12; vs sympy integral exact) · V2 (D13 all 7 lines) · V3 (sampled trapezoid with a jump: order 0.98) · V7 (d'Alembert 1e-14; C_b monotone) | order 0.98 (design 1) | analytic, qualitative (model) | see G5: at fixed C_b the drag is not monotone in φ_s |
| C10 Kármán street | CORE | `karman_street_ratio/velocity/spectrum/spectrum_periodic/growth/growth_closed/positions`, `ch09.karman_street_sympy`, `shedding_frequency` | V1 (arccosh√2/π; growth (πΓ/2a²)\|½ − sech²(πb/a)\|; street speed by direct 80 001-vortex Biot–Savart sum: 2e-4) · V2 (D14 all 15 lines incl. facing rows; matrix and its factorised polynomial by sympy) · V3 (Bloch 4×4 vs finite periodic cell vs closed form: 1e-12; lattice truncation order 1.99) · V5 (0.281 corroboration, 1.6e-3 = 3 printed digits) · V7 (facing rows π/4 for every b/a) | order 1.99 | analytic, symbolic, converged | |
| C11 drag crisis, sports balls | CORE | `cylinder_flow_regime`, `sphere_flow_regime`, `cylinder_state`, `drag_crisis_state/pair`, `cylinder_cd_schematic`, `ball_swing_deflection`, `magnus_sign` | V1 (swing kinematics, model closed forms) · V7 (regime lookup at every threshold, φ_sep 82° → 125°, roughness triggers at Re_cr/3, Magnus truth table incl. R10) · V5 (reused Morrison sphere correlation: crisis dip C_D 0.09 vs 0.43 at 1e5; Stokes/Oseen bracket) · V6 (thresholds, angles, wire 1000 Hz, cricket 0.2 %) | model ratio 0.65 | qualitative (flagged), analytic | conceptual/illustrative by design — see Open items |
| C12 free jet (9.53)–(9.76) | CORE | `free_jet*`, `free_jet_ode_solve`, `free_jet_constants`, `ch09.jet_momentum_sympy` | V2 (sympy: tanh solves 3f‴ + ff″ + f′² = 0, first integrals, C = 4√6/3, D15–D18) · V3 (`solve_bvp` vs √6 tanh: 4e-13; tolerance 1e-4 → 1e-10 monotone; Dirichlet truncation e12 ≫ e18 ≫ e24; (9.18) and continuity on the field: order 1.87/1.91) · V4 (ρ∫u²dy = J to 1e-10 at 6 stations 0.01…10 m; ṁ slope 1/3, u₀ −1/3, δ 2/3) · V5 (Bickley 0.4543, 0.2752, 3.3019, 0.5503: ≤ 1.4e-4 = printed digits) | 1e-10 | symbolic, converged, conserved, benchmark | R6: printed 5.6152 puts u at 4 % not 1 % of u₀ — planted variant fails |
| C13 wall jet (9.77)–(9.85) | CORE | `wall_jet*`, `wall_jet_ode_solve`, `wall_jet_profile`, `wall_jet_invariant`, `ch09.wall_jet_sympy` | V2 (sympy: 4f‴ + ff″ + 2f′² from ψ; first integrals; partial fractions; D19–D21 numerics on the exact solution) · V3 (IVP vs the implicit (9.83) by brentq: 6.5e-12 at η_max = 120; PDE/continuity order 1.96/1.97) · V4 (invariant constant in x to 1e-5 with a different quadrature per station; ρ∫u²dy ∝ x^{−1/4}; ṁ ∝ x^{1/4}) · V7 (f_∞³ = 72 f″(0) to 1e-9 at three scales; scaling f → λf(λη) exact to 1e-12) | 1e-5 (quadrature) | symbolic, converged, conserved | R3 (coefficient 1 gives f_∞ = 0.397 and violates (9.83) by 0.31), R4 (printed integrand fails); F4 fixed (default η_max) |
| C14 teacup | CORE | `secondary_flow_radial_force`, `secondary_flow_layer_profile` | V1 (1000, 750, 437.5, 0 N/m³ exact) · V2 (D22: r-momentum of a pure swirl by sympy; dimension check) · V7 (demonstration: four illustrative profiles vanish on the floor, force ≥ 0, largest on the floor, zero in the core) | exact | analytic, qualitative (profile) | conceptual demonstration |

NOTE rows that are coded (all ≥ 1 level): N04/N06/N07 (scales, (9.6)), N08/N09 (sympy table), N10 (`bl_pressure_variation`), N11 ((9.11), (9.35)), N14/N24 (marching, memory: erasing of the inlet 1e-5 for n = 0.5, 2 % plate, 7.6 % for n = −0.05), N16 (2/√Re vs 0.664: factor 3.0), N19/N21, N33 (far field within 0.3 %), N34, N35 (R2), N36–N40, N41 (0.664 < 1.128 < 1.328), N42, N45–N47 (R11), N55, N56–N67, N68/N69 (plate curve; V5 secondary), N70–N73, N76–N86, N87–N124 (all jet notes). RECAP rows R01–R09 are tested in their own chapters.

## Functions used by the notebook and explainers (design Part C) — each is called by at least one test
| Part C rows | tests |
|---|---|
| 1.1–1.5, 1.20, 1.22–1.24 | `test_bl_scales_*`, `test_bl_variables_*`, `test_bl_x_momentum_V3_*`, `test_bl_pressure_V7_*`, `test_outer_flow_V1_*`, `test_thwaites_*`, `test_inflection_V7_*`, `test_separation_point_V1_*`, `test_march_*`, `test_transition_state_V7_*`, `test_plate_drag_*` |
| 1.6–1.7 | `test_thicknesses_*`, `test_profile_shape_V1_*`, `test_delta_level_V1_*` |
| 1.8–1.16 | `test_blasius_*`, `test_falkner_skan_*`, `test_hiemenz_V5_*`, `test_falkner_skan_table_V1_*` |
| 1.17–1.19, 1.20b, 1.21 | `test_momentum_integral_*`, `test_karman_pohlhausen_*`, `test_thwaites_named_*`, `test_thwaites_closure_*`, `test_thwaites_cylinder_*` |
| 2.1–2.12 | `test_free_jet_*`, `test_wall_jet_*` (incl. `wall_jet_constants`, `wall_jet_integrals`, `wall_jet_first_integral_residual`) |
| 3.1–3.8 | `test_cylinder_regimes_V7_*`, `test_drag_crisis_V7_*`, `test_separated_drag_*`, `test_karman_*`, `test_shedding_V1_*`, `test_ball_dynamics_V1_*` |
| 4.1–4.10 | `test_bl_nondim_V2_*`, `test_sympy_engines_V2_*` (thwaites, jets, cylinder, `derive_all`), `test_karman_street_V2_derivation_D14`, `test_thwaites_V1_flat_plate_and_example_9_1`, `test_thwaites_V1_diffuser_*`, `test_teacup_*`, `test_book_slips_V1_*` |
| all of Part C | `test_part_c_V1_every_contract_function_exists_and_is_scalar_callable`: every public name of BL/JET/BB is re-exported by `ch09`; 38 scalar-callable functions return finite floats; every dict-returning function has its documented keys |
| reused C.0 | called with ch09 inputs in the same smoke test (`reynolds_number`, `similarity_variable`, `diffusion_thickness`, `pressure_coefficient`, Oseen/Stokes, `temporal_bl_wall_stress`, ch06 sphere/cylinder C_p, `potential.cylinder`) |
| C.5 scripts | `test_scripts_V1_every_ch09_script_runs` (slow): all 22 `ch09_*.py` demos + `ch09_common.py` helper: exit 0 with `--no-show` in 97 s. `scripts/ch09_drawings.py` (design C.5) now exists and `test_drawings_V1_*` imports it and draws all five helpers on Agg |

## Convergence studies (observed order, `tools/convergence.observed_order`)
| scheme | grids | observed order (pairwise) | design order |
|---|---|---|---|
| central differences of (9.18) on `blasius_fields` (interior, two boundary nodes excluded) | n = 20, 40, 80 | 1.77 (1.69, 1.85) | 2 (±0.25) |
| continuity u_x + v_y, same field | n = 20, 40, 80 | 1.78 (1.70, 1.85) | 2 |
| (9.18) on `free_jet` field | n = 40, 80, 160 | 1.87 (1.82, 1.93); continuity 1.91 | 2 |
| (9.18) on `wall_jet` field | n = 40, 80, 160 | 1.96 (1.95, 1.97); continuity 1.97 | 2 |
| marching, backward Euler in Δx (plate, ny = 800) | nx = 10, 20, 40 | 1.006 (1.009, 1.004) | 1 |
| marching, central in Δσ (**plate**) | ny = 100…800 | 2.03 (1.82, 1.98, 2.29) | 2 |
| marching inlet-station τ₀ (Blasius inlet, x₀ = 0.1) | ny = 100…800 | ≈ 3.0 (monotone: 2.3e-4, 2.9e-5, 3.6e-6, 4.5e-7) | ≥ 2 |
| marching, central in Δσ (**wedge n = 0.5**) | ny = 100…800 | **2.24** (2.14, 2.21, 2.38); was 0.98 before fix F3 | 2 (±0.25) |
| marching, BDF2 in Δx (wedge n = 0.5, self-reference) | nx = 20, 40, 80 | 2.20 (2.24, 2.15) | 2 (±0.25, span 1.3 decades) |
| marching, BDF2 in Δx (wedge n = −0.05) | nx = 40, 80, 160 | 1.87 (1.82, 1.92) | 2 (±0.25) |
| marching, backward Euler (wedge n = 0.5; −0.05) | as above | 1.04; 1.06 | 1 |
| separated-drag trapezoid with a jump at a node | Δφ = 2, 1, 0.5, 0.25° | 0.984 | 1 |
| Kármán lattice sums truncated at ±N | N = 500…4000 | 1.99 (errors 1.6e-6 → 2.5e-8) | 2 |
| Blasius truncation η_max = 8, 10, 12 | f″(0) error 1.9e-6, 9.9e-10, 7e-12 | Gaussian tail e^{−η²/4} (not a power law) | — |

## Conservation / invariant residuals
- Momentum flux of the free jet: ρ∫u²dy / J − 1 < 1e-10 at x = 0.01, 0.05, 0.1, 0.5, 2, 10 m (adaptive quadrature); trapezoid helper 1e-8; ṁ = (36Jρ²νx)^{1/3} by quadrature 1e-10.
- Wall-jet invariant ∫u(∫_y u²)dy / Ψ − 1 < 1e-5 at five stations with a different grid at each (measured 7.8e-6 → 8.6e-7 with resolution); ordinary momentum flux slope −0.25 (± 2e-5), mass flux slope +0.25.
- Momentum thickness vs wall friction: ρU²θ/∫τ₀dx − 1 < 1e-6 (Blasius, 60 stations); (9.43) residual/τ₀ < 1e-7 for Blasius and Falkner–Skan; marched layers 7e-4 (plate, sine inlet), 1.3e-3 (weak diffuser).
- D19 numerics on the exact wall-jet solution: T₁ + T₂ = 0 and 2T₁ = ∫u²v dy to 1e-4 of the integrand scale; d/dx of the double integral < 2e-4 |A₁|.
- Exact impossibility checks: d'Alembert (∮ ideal C_p cos φ dφ) 1e-14; Kármán marginal growth 1e-16.

## Benchmarks used (public sources; see `reference/ch09/SOURCES.md`, all read 2026-09-30)
| value | ours | source | relative difference |
|---|---|---|---|
| Blasius f″(0) = 0.332057336215196 | 0.3320573362151941 | Wikipedia "Blasius boundary layer" | 5.7e-15 |
| Blasius κ = 0.469599988361 (β = 0, their normalisation) | 0.469599988361036 | Belden et al., arXiv:1907.09912 | 7.7e-14 |
| FS β = 0.5: 0.927680039836653; β = −0.12: 0.28176052424 | …652; …0411 | same | 1.1e-15; 1.5e-12 |
| FS separation β = −0.198837735 | −0.198837735047 | same | 2.4e-10 |
| FS reversed branch β = −0.12: −0.1429351943576; β = −0.02: −0.065168585542904 | −0.1429351943568; −0.0651685855429 | same | 5.5e-12; 6.4e-13 |
| Hiemenz F″(0) = 1.232588 | 1.2325876568 | Weidman & Turner (Surrey preprint) | 2.8e-7 |
| Hiemenz δ*/δ = 0.6479 | 0.6479005 | Wikipedia "Stagnation point flow" | 7.3e-7 |
| Bickley 0.4543 / 0.2752 / 3.3019 / 0.5503 | 0.45428 / 0.27516 / 3.30193 / 0.55032 | Wikipedia "Bickley jet" (Bickley 1937) | 4.4e-5 / 1.4e-4 / 8.3e-6 / 3.9e-5 |
| Thwaites 0.45 + 6m, θ² form, m_sep ≈ 0.09 | exact (forms) | Agrawal et al., arXiv:2310.16337 | form check |
| turbulent flat-plate law | 0.074 Re^{−1/5} vs 5/4 × 0.0576 = 0.0720 | Wikipedia "Skin friction drag" (secondary) | +2.8 % (inside the 3 % secondary band) |
| Kármán b/a ≈ 0.281 | 0.28055 | Horváth et al. 2020 (search text) | 1.6e-3 (3 digits) |
Not found today (recorded, nothing asserted): Glauert (1956) wall-jet constants, Lienhard (1966) cylinder C_D, Roshko Strouhal numbers, Howarth separation angles, a primary source for the turbulent-plate coefficient. Note: the Blasius Wikipedia page prints a second, different number 0.332043934904293 for the boundary-condition conversion; it is 4e-5 from the converged constant and is not used.

## Numbers from the text (book vs ours) — private values redacted to relative differences
| item | difference | comment |
|---|---|---|
| η₉₉ (9.30) | 0.41 % | slip R2: figure reading; the root of f′ = 0.99 is 4.910 |
| δ*, θ, τ₀, C_f, C_D, v∞ coefficients of §9.3 | 0.05 %, 0.02 %, 0.02 %, 0.02 %, 0.13 %, 0.05 % | printed with 3 digits |
| air example: δ₉₉ for Re_x = 6e4 | 0.25 % (ν implied by Re_x); ν = 1.5e-5 gives 1.9 cm | the book's ν ≈ 1.67e-5 |
| C_f ≈ 2/√Re (crude) vs Blasius | factor 3.0 | as the book says |
| n_sep = −0.0904 | 0.03 % | |
| stagnation C_f√Re_x (exercise value) | 1.0e-5 | 5 printed digits |
| Fig. 9.7 members (7 exponents) | all attached solutions | |
| Example 9.1: θ, l(0), C_f, θ error, C_f error | 0.03 %, 0.2 %, 0.2 %, 0.0 pp, 0.2 pp | |
| Example 9.1: H(0), δ* | 0.7 %, 0.7 % | the table's H(0) is 2.61, exact Blasius 2.591 (ours from the FS closure) |
| Example 9.2: λ at x/L = 0.05…0.2 | ≤ 1.8e-4 | closed form |
| Example 9.2: x_sep/L | 1.1 % | the book prints 0.16 (2 digits); closed form 0.158 |
| Table 9.1 vs the "white" closure (l, −0.06 ≤ λ ≤ 0.25) | ≤ 2.5 % | |
| Table 9.1 vs our exact-FS closure | l: ≤ 6.0 % (0 ≤ λ ≤ 0.1), ≤ 3.6 % (−0.016 ≤ λ < 0), 6.4 % at −0.024, 14.5 % at −0.04, 43 % at −0.06; H: ≤ 4.7 % (0…0.1), 12 % at −0.06 | **the design/analysis statement "differs by < 5 % where it overlaps" is false** (G1) |
| Thwaites accuracy claim ±3 % (favourable) / ±10 % (adverse) | θ errors −1.6 % … −7.6 % for n = 0.1 … 4 exceed 3 %; τ₀ errors ≤ 2.4 % (FS closure) meet it | (9.50) θ does not depend on the closure; the fit L = 0.45 − 6λ is a 6–8 % effect on θ at strong acceleration |
| regime thresholds, separation angles 82°/125°, St = 0.2, wire 1000 Hz | equal (exactly the rounded values); Kármán 0.28 vs 0.28055: 0.2 % | rounded experimental values |
| drag crisis C_D 1.2 → 0.33 (× 0.28) | model pressure drag × 0.65 | illustrative model, same direction, weaker magnitude (Open item) |
| sphere critical Re | code 3e5 vs book 5e5 (40 %) | design contract value; not asserted (Open item) |
| cricket ball 0.4 W, 30 m/s: swing over the flight time | 0.2 % | |
| potential-flow sphere C_p,min = −1.25 | exact | |
| exercise-type Blasius point (x = 0.15 m, U = 6 m/s, air) | y-η pair 0.02 %, v 0.9 %, ∂u/∂y 0.1 % | 2–3 printed digits |
| jets: C = 4√6/3, ṁ coefficient 36, printed h₉₉ argument/coefficient | exact / exact / 1.6e-5 (printed 5.6152 = √6 arccosh 5 to 4 digits) | slip R6 explained |

## Figures reproduced with our code (`outputs/ch09/verify/`, local, git-ignored)
| Fig. | PNG | visual verdict (each PNG was opened and read) |
|---|---|---|
| 9.5, 9.6 | `fig_9_5_9_6_blasius.png` | f′ rises smoothly from 0 to 1 with f″(0) = 0.332 at the wall, the 99 % line at η = 4.910 meets f′ = 0.99 and the printed 4.93 line sits just above; δ* and θ lines at 1.72 and 0.66; dimensional profiles at 0.1, 0.5, 2 m stretch as √x; v√Re_x/U saturates at 0.8604 |
| 9.7, 9.8 | `fig_9_7_9_8_falkner_skan.png` | seven FS profiles thicken monotonically from n = 4 to n = −0.0904 whose wall slope is zero (S-shaped); f″(0) vs n is monotone, the Belden points lie on the curve, the fold is at f″(0) = 0 with the reversed-branch triangles below it; L(λ) is a straight line parallel to 0.45 − 6λ, 0.01 below at λ = 0 and 0.06 above it at λ = 0.085 |
| 9.9, D11 | `fig_thwaites.png` | diffuser λ(x) overlays the closed form and crosses −0.0681 at x/L = 0.126 and −0.09 at 0.158; θ error is +4.4 % near the fold, +1.0 % at Blasius and falls to −7.6 % at n = 4 (outside the book's ±3 %, inside our 8 % bound); cylinder λ(φ) starts at 0.075, passes 0 at 90° and crosses −0.09 at 103.11° with the quadrature dots on the curve |
| — | `fig_marching.png` | wall-shear error vs Δσ: the plate and the wedge n = 0.5 both fall with slope 2 (F3 fixed); Δx errors: BDF2 slope ≈ 2, backward Euler slope 1; the τ₀(x) of wedge flows n = −0.095 and −0.12 goes to zero and stops (separation; x = 3.43 and 0.27) while n = −0.089 stays attached to x = 50 |
| 9.16–9.22, 9.11 | `fig_bluff_and_street.png` | Kármán growth is a V with its zero at b/a = 0.2805 and the three numerical routes on the curve, facing rows flat at π/4; the spectrum at 0.2805 is two pure-imaginary pairs; model C_p follows 1 − 4sin²φ to the separation angle then a plateau (jump at φ_s); sphere Morrison curve shows the crisis dip to 0.09 at Re ≈ 4e5 and the Stokes/Oseen asymptotes; cylinder schematic is now continuous at Re = 1 (anchor = Lamb's 12.55); model φ_sep steps 82° → 125° (rough curve earlier); the plate curves: laminar 1/√Re line, turbulent 0.074Re^{-1/5}, mixed curve leaving the laminar line at 5e5 and joining the turbulent one near 1e7 (the model cylinder line in that panel is a plotting artefact of the script, not a claim) |
| 9.29 etc. | `fig_jets.png` | free-jet profile: `solve_bvp` (dashed) is indistinguishable from sech²(η/√6) (4e-13); exponents −1/3, 2/3, 1/3 are straight lines on log axes; the wall-jet f′ peaks near η = 8 and f rises to 1, the brentq inversion of (9.83) lies on the IVP, the printed ODE saturates at 0.397; the invariant is flat while ∫u²dy falls like x^{−1/4} |

## Deviations & justifications
Every `DEVIATION` in `fluidpy/`: `bl_pressure_variation` (denominator = the O(1) size of the terms ∂p*/∂x* would balance, because ∂p*/∂x* = 0 on a Blasius field — accepted, exactly the order-of-magnitude reading; measured order 1.00), `falkner_skan` (BVP/Töpfer/shooting are our choices — see F1 and the η_max = 10 truncation note under "Smaller code observations"), `karman_pohlhausen` (shape integrals by sympy: verified exactly), `thwaites_closure_table` (the exact-FS closure replaces the private Table 9.1 — see G1: it differs by up to 6 % in l over 0…0.1 and much more near the fold; the "white" fit follows the table to 2.5 %), `thwaites` (Gauss–Legendre for ∫U_e⁵: verified 1e-12 against the closed forms).
Verifier choices (tolerances not loosened): the (9.18) residual orders exclude two boundary nodes because one-sided edge stencils cost an order in u_yy (the compact interior stencil converges at 2.0); the BDF2 marching order is measured by self-reference at fixed Δσ (a flat-plate solution is super-convergent, 3.4, and is not used for that claim).

## Open items (verbatim for the orchestrator)
**Defects: none open.** F1-F4 are fixed and verified (table at the top); the smaller code observations of the FAIL report are resolved:
- march_boundary_layer inlet-station tau0 now falls with refinement (2.3e-4 ... 4.5e-7, `test_march_V3_inlet_station_*`); a local m0 below the fold raises a clear ValueError that names the separation value and the u_inlet way out (`test_march_V7_below_the_fold_*`).
- falkner_skan_state(-0.05) now uses eta_max = 16 (switch is m <= -0.05); falkner_skan(0.0) default eta_max = 12 agrees with the Toepfer value to < 1e-10 (`test_falkner_skan_V5_default_eta_max_*`).
- free_jet_ode_solve(bc="dirichlet") docstring now says 2.4e-3; thwaites_closure_table docstring now says 2e-8; separated_pressure_drag docstring now says 0.8838 and carries the non-monotone-in-phi_s note (G5); cylinder_cd_schematic is continuous at Re = 1; scripts/ch09_drawings.py exists (imported and drawn by a test).

**Remaining observations (no test fails):**
- Wedge n = -0.05 wall shear at fixed x-stations stalls at about 3e-4 as ny grows (x-discretisation and inlet memory dominate); the marched separation position for n = -0.095 depends strongly on the inlet profile (2.03 before the F3 fix, 3.43 after) because that exponent is only 5 % below the fold. Do not quote a marched separation x near the fold as a number; quote the bracket (n = -0.089 attached to x = 50, n = -0.12 separated at 0.27).
- D10 `thwaites_closure_table` fast=True (24 members) was not separately compared to the 60-member table beyond the smoke test.
- The whole-tree single-command `pytest -q` stalls near 96 % when the browser-based machinery test runs after the long ch09 script test (harness time limit); run the machinery file separately (19 s) if this recurs. Not a ch09 defect.

## Documentation corrections required (for the orchestrator / notebook builder; the design and analysis documents are not edited by the verifier)
G1-G9 — the designer's lines that were checked and are wrong or imprecise (the notebook and explainers must use the corrected statements):
- G1 D10 step 6/"Check" and analysis §6: "Book's Table 9.1 differs by < 5 % where it overlaps" is false for l (6.0 % over 0 ≤ λ ≤ 0.1; 6.4 % at −0.024, 14.5 % at −0.04, 43 % at −0.06); true for H over 0…0.1 (4.7 %). "reproduces to 1e-9" → 2e-8.
- G2 D10 "Result": "L(λ) ≈ 0.45 − 6λ within 0.06" — measured 0.063 at λ = 0.0855 (0.05 on −0.068 ≤ λ ≤ 0.05).
- G3 D06 step 6 line: "f″(0) = ½[I_θ(η) − f(1 − f′)]" is only true as η → ∞; the exact finite-η relation is f″(0) = ½I_θ(η) − ½f(1 − f′) + f″(η) (sympy: E′ = 0 by (9.27); numeric check at η = 1, 2, 4 to 1e-11). The end result θ = 2f″(0) is right.
- G4 D14 step 15: "det(M − μI) = (μ² − π²/16)²" — in the eigenvalue μ of M/(π/2) used in steps 11–14 the polynomial is (μ² − ¼)²; (λ² − π²/16)² holds for λ = (π/2)μ. The eigenvalues ±π/4 and the conclusion are right (sympy for both).
- G5 Design V7 (analysis §6, Part C 3.3 note) "later separation ⇒ smaller drag at fixed base pressure": false on this model (D(82°) = 0.884, D(90°) = 0.867, D(125°) = 1.07 for C_b = −1.2); only the base-pressure rise makes the drag fall. The docstring of `separated_pressure_drag` repeats it.
- G6 D11 step 6 "Thwaites is 12 % low at the stagnation point": 1 − 0.075/0.0855 = 12.3 % (fine), but the reverse ratio is 14 %; state which.
- G7 Every other displayed line of D01–D22 (D01 all 14, D02–D05, D07–D09, D11–D13, D15–D22) was re-derived in sympy and matches; the only wrong lines are G3 and G4 above (and the imprecise G2, G6).
- G8 Analysis §8 "Hiemenz value unverified" — now verified (Weidman & Turner 1.232588; Wikipedia δ*/δ 0.6479); the Belden et al. table also supplies the *reversed branch* values that verify R11.
- G9 Curation §8 item 5 / design C.3 3.1 quote 3e5 for the sphere critical Re; the book text (private) says 5e5 — one of the two should be changed or flagged in the regime table.

**Qualitative / non-numerical (must stay labelled qualitative in the notebook and explainers):** regime thresholds (Re = 1, 4, 40, 200, 3000, 3e5, 6e5; sphere 130, 3e5, 8e5), separation angles (82°, 125°, 80°, 120°), St = 0.2 band, base pressures (−1.2, −0.6), `cylinder_cd_schematic`, the model drag-crisis ratio (0.65 vs the book's 0.28) and `transition_state` thresholds (transitional band one decade wide; the book's "fully turbulent above 1e7" is 2× the code's 5e6). Turbulent flat-plate coefficient 0.074 is secondary-sourced (5/4 × 0.0576 = 0.0720; +2.8 %) — cite a primary source before the notebook prints it. Wall-jet constants have no independent published number (Glauert 1956 not accessed): verified by V2 + V3 (two routes) + V4 + V7 only.

## Discrimination (planted wrong variants, run on a scratch copy of `fluidpy/`; the selected tests pass on the unmutated copy and fail on each mutant)
25 mutants, **all killed**: θ = f″(0) (factor 2) · Thwaites 0.45 → 0.5 · Thwaites U⁶ → U⁵ · Falkner–Skan I_θ computed as ∫f′ (not f′(1−f′)) · free-jet ṁ coefficient 36 → 32 · half-width defaulting to the printed 5.6152 · wall-jet ODE default coefficient 4 → 1 · Kármán growth sech² → sech · (9.43) without U_eδ*dU_e/dx · pressure drag with sin in place of cos · η₉₉ = 4.93 by default · diffuser sign flipped · wall-jet invariant with u instead of u² · Holstein–Bohlen λ = mI_θ (not squared) · wall-jet similarity engine with coefficient 1 · f = St·U·d · δ̄ = √(νLU) · cylinder λ with sin · entrainment √6/2 · wall curvature sign · C_f = f″(0) · drag on both faces by default · free-jet u ∝ sech (not sech²) · u₀ ∝ x^{−1/2} · δ* without the exponential tail. One further mutant (sign of the Falkner–Skan forcing term in the ODE right-hand side) made `solve_bvp` loop and the selected tests did not terminate within 900 s (a hang, not a clean failure; the suite would be stopped by the timeout).
