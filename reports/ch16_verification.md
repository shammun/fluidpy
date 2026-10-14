# Chapter 16 verification — Introduction to Biofluid Mechanics                2026-10-14 · loop 3 (after the derivation review), working tree above `76eb2ae`

## Verdict: FAIL — `tests/test_ch16.py -W error::RuntimeWarning`: 98 passed, 1 failed, 0 skipped in 333 s (slow script test and private V6 test included). The one failure is a new robustness defect found in this pass (`collapsible_frictionless_alpha` raises `ValueError` on valid wide upstream sections, finding L3-1 below); a one-line fix makes all 9 collapsible tests pass under `-W error` in a scratch copy. Everything the review asked for is in place and passes. `tests/test_machinery.py` 9 passed; `tools/check_public.py` OK (818 files)

### Loop 3 (after `reports/ch16_review.md`) — what changed and how it was checked
The implementer changed `fluidpy/ch16_biofluid_mechanics.py`, `fluidpy/core/womersley.py` (used only by ch16; grep finds
no other importer) and `scripts/ch16_collapsible.py`. The verifier changed only `tests/test_ch16.py` (91 → 99 tests),
`tests/ch16_verify_figures.py` and this report.

| Review item | Test(s) | Result |
|---|---|---|
| Must-fix 1: range guard for (16.188), α ≤ 1 | `test_collapsible_V1_tube_law_wave_speed_and_M` changed: α = 1.7 removed (grid 0.15, 0.4, 0.8, 0.97, 1.0); ∂P/∂α by a second-order **backward** stencil, so α = 1 never steps above 1; `shapiro_M`'s central differences used only at α < 1 | passes under `-W error::RuntimeWarning`; c² = (α/ρ)K_p∂P/∂α (16.198) and Bramwell–Hill (16.185) to 1e-8, round trip to 1e-13, both forms |
| | new `test_collapsible_V7_tube_law_range_guard_alpha_le_1` | scalar α = 1.7: RuntimeWarning naming (16.188) and the extrapolated value (1e-15) for `collapsible_tube_law` (both forms), `_derivatives`, `collapsible_wave_speed`, `alpha_from_tube_law(P = 0.3)`, `collapsible_terms(1.3, …)`; array [0.5, 1, 1.7, 0.9, 2.5]: NaN **exactly** at the two α > 1 entries, the others unchanged (1e-15); silent at α = 1, α = 1 + 1e-13 (round-off allowance), α ≤ 1 grids, `alpha_from_tube_law(0) = 1` |
| | new `test_collapsible_V2_throat_condition_alpha_t_le_0397` | sympy: H = K_pP(α_t) + ½ρc_t² = K_p[1 − (1 − n/2)α_t⁻ⁿ], H = 0 at α_t = (1 − n/2)^{1/n}, dH/dα_t = K_pn(1 − n/2)α_t^{−n−1} > 0 (n < 2); at n = 3/2 the bound is 0.25^{2/3} = 0.39685. Numerically (H ≤ 0) ⇔ (α_t ≤ 0.39685) on α_t = 0.2 … 0.6; the demo's numbers reproduced: H = −41.47 Pa, Q* = 47.646 cm³/s, α = 0.8758 at x = −12 cm (5e-5), the subcritical branch rises monotonically towards α_far = 0.88194 < 1 with **no warning** (warnings → errors); α_t = 0.45 (H > 0) gives an upstream root > 1 **with** the RuntimeWarning |
| (2) frictionless throat test | `test_collapsible_V1_frictionless_throat_chokes_at_the_throat` moved to α_t = 0.35 (asserts H < 0) and marked `filterwarnings("error::RuntimeWarning")` | passes, no warning; march = Bernoulli to 1e-6, driver 0 at the throat (1e-9) |
| Should-fix 1: scalar G | new `test_womersley_flow_rate_V7_scalar_gradient_has_one_meaning` | `womersley_flow_rate_harmonics` raises ValueError for 100.0, `np.float64`, 0-d array, int and complex scalars; `womersley_flow_rate(G)` = `…_harmonics([0, G])` (1e-14); `[G₀]` = Poiseuille πG₀a⁴/(8μ) in both functions (1e-14) |
| Should-fix 5: V3 names without an order | renamed to V1 (independent-route cross-check): `test_distensible_flux_V1_ode_route_and_rigid_limit`, `test_windkessel_V1_periodic_state_reached_by_marching`, `test_elastic_tube_V1_independent_4x4_determinant`, `test_murray_V1_numerical_branch_point_reproduces_angles`, `test_junction_V1_grid_pulse_matches_theory` (now also asserts nx = 200 and 400 give the same R, T to 1e-12: C = 1 characteristic transport is exact, there is nothing to refine); `test_tube1d_V3_printed_16_54…` → `test_tube1d_V7_…` (structure: a simple wave leaves rest behind); the `collapsible_march` test's header comment "V1 + V3" → "V1" | all pass |
| | **new genuine order studies** | `test_distensible_flux_V3_quadrature_of_16_44_orders_two_and_four`: trapezoid and Simpson on the corrected (16.44) with the module's a(p), p₁ = 30 kPa (β to 0.34), n = 8 … 64 → orders **1.999** (pairwise 1.998, 1.999, 2.000) and **3.989** (3.973, 3.993, 3.998), Simpson error 8.8e-9. `test_windkessel_V3_crank_nicolson_march_order_two`: trapezoidal march of (16.49) from rest over three beats, Δt = T/40 … T/320 (T_s on a node), against the closed form (16.52) → order **2.000** (pairwise 2.001, 2.000, 2.000), finest error 6.5e-5 of the peak. `test_elastic_tube_V3_ode_determinant_series_start_order_six`: residual of `elastic_tube_determinant_numeric` at an exact root of (16.136) against the series start s₀ = 0.2 … 0.025 → order **6.000** (α = 2; pairwise 6.000, 6.000, 5.999) and **6.023** (α = 8; 6.069, 6.005, 6.000), the design order of the s⁴-truncated start (first dropped term β̂⁶s⁶/2304); residual 1.6e-14 / 6.5e-13, ≥ 1e8 smaller than off the root. `test_dean_V3_spectral_resolution_study_and_fit_truncation`: flux ratio at x = K/576 = 0.35 against N = 25: N = 9, 13, 17, 21 → 1.1e-4, 4.2e-7, 1.0e-8, 3.2e-11 (geometric, spectral); M = 8, 12 against 16 → 2.8e-10, 7.9e-14; the c₄ fit against the exact 6471982981/542442700800 with 3, 4, 5, 6 fitted terms → 4.7e-3, 2.3e-4, 9.4e-6, 3.3e-7 (7 terms over-fits, 8.2e-7); c₂ to 1.4e-10 |
| `tests/ch16_verify_figures.py` | tube-law panel limited to α ∈ [0.15, 1] with a dotted line at α = 1 | runs with `-W error::RuntimeWarning`; panel verdict below |
| `scripts/ch16_collapsible.py` | run headless with `-W error::RuntimeWarning` | exit 0, no warning (throat mode α_t = 0.35; Starling Q_choke 51.772 cm³/s) |

**The new tests discriminate.** On the pre-fix module (`git show HEAD:` of both changed files, scratch tree) the guard
test, the throat-condition test (no warning for α_t = 0.45) and the scalar-G test fail; the rewritten
`tube_law_wave_speed_and_M` and the moved throat test pass there too (they only stopped needing α > 1).

**Finding L3-1 (implementer; the reason for FAIL) — `collapsible_frictionless_alpha` raises on wide upstream sections.**
It brackets the critical area ratio with a fixed `brentq(S2 − 1, 1e-6, 1e6)`. For the power law,
α* = (ρQ²/(A₀²K_pn))^{1/(2−n)} = α_t(A₀,t/A₀)⁴ at n = 3/2, so with the demo's throat (α_t = 0.35) α* < 1e-6 once
A₀/A₀,t ≳ 24.3 and the call raises `ValueError: f(a) and f(b) must have different signs` instead of returning the
subcritical root (0.88 there, near α_far = 0.88194), contrary to its docstring ("NaN if H is below the minimum").
`test_collapsible_V1_frictionless_alpha_on_wide_upstream_sections` fails with: *A₀/A₀,t = 25 (α* = 8.96e-07), 30 (4.32e-07),
100 (3.50e-09): ValueError*. The demo (|x| ≤ 12 cm, A₀/A₀,t ≤ 6.8) and the notebook range are not affected; a user
sliding the upstream section out (or an explainer that does) is. **Fix (one line, verified):** replace the bracketed
search by the closed form `a_star = (rho*Q**2/(A0**2*Kp*n))**(1/(2-n))` (valid for both forms — they share ∂P/∂α), or
widen the bracket adaptively. With that line in a scratch copy, all 9 collapsible tests pass under `-W error`, including
this one (Bernoulli residual < 1e-9 K_p, α increasing with A₀ towards α_far, within 1e-4 of it at A₀/A₀,t = 100).

**Docstring labels after loop 3 (for the implementer, not verdict-blocking).** `elastic_tube_determinant_numeric` now
has a genuine order study (order 6 in s₀) — its label may say `converged` again (V3 order 6 in the series start);
`elastic_tube_flow_rate` is the reference of the new quadrature V3 study (its own label "symbolic, analytic" is right).
`test_docstrings_V7_no_stale_validation_text` passes.

### Loop 2 — what changed and how it was checked (kept for the record)
The implementer applied Must-fix 1 and Should-fix 2 of loop 1 (below, kept for the record). Only
`fluidpy/ch16_biofluid_mechanics.py` changed (47 + / 20 − lines), so the 44-minute full suite was not re-run; loop 1's
full run was **2176 passed, 1 failed** (the failure was this Dean test) with chapters 1–15 untouched.

| Item | Check | Result |
|---|---|---|
| **S34 — (16.174) advection sign** (slip, changes physics) | `test_dean_flux_V3_order4_coefficient_is_the_navier_stokes_value` (unchanged since loop 1) | **passes**: `DEAN_C4_DERIVED` = 6471982981/542442700800 = 0.01193118; spectral Navier–Stokes c₄ = 0.0119311 (rel. 9e-6) |
| `printed_16_174=True` reproduces the wrong variant | new `test_dean_flux_V2_V3_printed_16_174_option_reproduces_the_wrong_variant` | `dean_series_sympy(4, printed_16_174=True)["c4"]` = 821328579941/67805337600000 = 0.0121130 exactly, and it **fails** the spectral NS value (+1.53 %, test requires > 1 %); default series c₄ equals `DEAN_C4_DERIVED` exactly; ψ₀, w₁ and c₂ = 1541/50400 identical in both variants, ψ₁ differs (where the sign first enters); `dean_flux_ratio(exact=True, exact_c4=True)` uses the corrected value (1e-15) |
| slips table | `test_tables_V7_structure_and_slips_with_printed_options` updated | 34 rows S1–S34, 12 of kind "slip"; S34 is "slip", `changes_physics` True, cites (16.174), `printed_16_174` and `DEAN_C4_DERIVED`; all six `printed_*` options listed in `book_slips()` **and** the module docstring ("34 rows") |
| Windkessel overflow (Should-fix 2) | new `test_windkessel_V7_periodic_state_warning_free_for_tiny_RK` (warnings → errors inside the test) | K = 1e-14, 1e-16, 1e-20: finite, non-negative, diastolic p < 1e-200 (no storage), scalar stays scalar, `windkessel_summary` finite; array shape kept |
| The new tests discriminate | pre-fix module (`git show HEAD:…`) in a scratch tree | all four fail on the pre-fix code: K⁴ test (+1.53 %), printed-option test (`TypeError`, no option), slips table (33 rows), Windkessel (`RuntimeWarning: overflow encountered in exp`) |
| Option wiring mutants | post-fix module, (a) sign always printed, (b) option ignored (always corrected) | both **killed** by the printed-option test (the constant-based K⁴ test alone would not see (a)) |
| Fix changes nothing else | pre-fix vs post-fix on calls the fix should not change | `windkessel_periodic_state` over t ∈ [−1, 3] s, 6 compliances × 3 heart rates, arrays and scalars: max relative difference **0.0**; `dean_flux_ratio` (book and exact-c₂ coefficients, K = 1 … 576): **0.0**; order-2 series c₂ identical |

**Mutation tally: 56 of 56 planted wrong variants killed** (loop 1's 52, with the Dean K⁴ test now included and killing
the pre-fix "printed-sign" module; plus the two option-wiring mutants and the pre-fix two-branch Windkessel evaluation;
the pre-fix Dean module counted once).

All 13 `scripts/ch16_*.py` exit 0 headless inside the slow test; `reference/ch16/make_refs.py` is idempotent (unchanged).

### What was checked, in one paragraph
Every computable CORE item (C01–C17) has at least two independent evidence levels, at least one of V1/V2/V3/V5; every ★★
and ★★★ derivation row of curation §4b (★★ D05, D06, D10, D11, D12, D14, D15, D16, D17, D18, D22, D26, D29, D32;
★★★ D13, D19, D20, D21, D24, D27, D30, D33) has a sympy test that re-derives it from the stated starting equations,
checking the intermediate lines of the ★★★ rows one by one (and the ★ rows D01–D04, D07–D09, D23, D25, D31, D34; D28 by finite differences); all 165
functions and constants of design Part C (C.1–C.10) are called by at least one test (an automatic contract test parses
Part C and fails if a name is missing); the six `printed_*` options (five from loop 1, `printed_16_174` from loop 2) give the wrong answers and a test shows each;
four independent numerical routes were written in the test file (Dean equations in Navier–Stokes form, a 4 × 4 Bessel
determinant for the viscous elastic tube, a closed-form Starling-resistor choke, a 30-digit Bessel evaluation of
Womersley's profile). **Mutation check (loop 1): 52 of 52 planted wrong variants killed; 56 of 56 after loop 2** (scratch copy of `fluidpy/`, one
variant at a time; in loop 1 the Dean K⁴ test was deselected because it failed on the unmutated code): factors of 2 and signs in every CORE
formula, the i^{3/2} of the Bessel argument, the half-cell axis row of the Crank–Nicolson solver, the junction node
update of the grid solver, Dean's w₁ sign and source factor, the Starling α* exponent, each printed slip swapped in for
its correction. The one variant first caught only by a V3 test (the Windkessel p₀ closure) is now also caught by the V2
test, after its periodicity check was tightened to p(T⁻) = p(0).

## Environment
Python 3.11.5 · numpy 2.4.6 · scipy 1.17.1 · sympy 1.14.0 · pint 0.25.3 · mpmath 1.3.0 (Windows 11, `.venv`).
Run: `.venv/Scripts/python.exe -m pytest tests/test_ch16.py -q -p no:cacheprovider` (≈ 330 s with the slow script test in loop 3; `-m "not slow"` skips the
script runs). Figures: `.venv/Scripts/python.exe tests/ch16_verify_figures.py`. References:
`.venv/Scripts/python.exe reference/ch16/make_refs.py` (idempotent).

## The implementer's four points — settled
| # | Point | Independent route used here | Result |
|---|---|---|---|
| 1 | Starling choking by d(ln α) integration: Q_choke = 1.99899 × 10⁻⁵ m³/s, p₂ − p_e = −946.13 Pa | Closed form: Bernoulli with laminar friction d(p + ½ρu²) = −R′Q dx with the power-law tube law integrates in α to x(α) = G(α_in) − G(α), G(α) = nK_pA₀²α^{2−n}/((2−n)8πμQ) − ρQ ln α/(8πμ); Q_choke solves L = G(1) − G(α*(Q)) with α* where S = 1 | **Confirmed.** Closed form 1.998988800752 × 10⁻⁵ m³/s, α* = 0.312271, p₂ − p_e = −946.12802 Pa; `flow_limitation_curve` agrees to 2e-16 (Q) and 1e-15 (p₂). Subcritical outlet states (`starling_state`) agree to < 1e-9; `collapsible_march` α(x) agrees with x(α) to 1e-8 m. The design's earlier x-march value (2.00290 × 10⁻⁵) was the error; E9's parity row can be tightened from rtol 5e-3 |
| 2 | Dean order-4 flux coefficient 0.0121130 vs the printed 0.0120 (+0.94 %) | (a) Chebyshev–Fourier Newton solution of the Dean equations in the **Navier–Stokes form** of the Imperial College M4A33 notes (independent of the book's (16.173)–(16.174)), 25 Dean numbers, fit in (K/576)²; (b) sympy: eliminating p₁ from the book's own (16.167)–(16.168) with (16.172) | **0.0121130 is the right expansion of the printed equations but not the physics.** (b) shows the printed (16.174) has the nonlinear term −(1/r)(ψ_r∂_α − ψ_α∂_r)∇²ψ where the momentum equations give **+**. (a) gives c₂ = 0.0305754 (= 1541/50400 to 1e-7, so the K mapping and the leading orders are right) and **c₄ = 0.0119311**; with the book's sign the same solver gives c₄ = 0.0121130 (agrees with the module to 2e-5). The module's own series engine with the sign flipped gives the exact **c₄ = 6471982981/542442700800 = 0.01193118**, matching (a) to 9e-6 relative. The printed 0.0120 is 0.58 % above the physical value (it rounds 0.01193 up — possibly Dean's own figure; Dean (1928) was not reached). ψ₀ and w₁ (sign and size) are unaffected and are confirmed by (a) to 1e-9 / 3e-6 |
| 3 | Elastic-tube characteristic (16.136) confirmed by an ODE determinant | (a) sympy: (16.128)–(16.129) from (16.126)–(16.127); (16.130)–(16.133); the rows (16.134), (16.135) from the boundary conditions; their determinant = (const/X)·(16.136); (b) a 4 × 4 determinant of the **un-eliminated** system (16.122)–(16.125), (16.128)–(16.129) in (A, η̂, ξ̂, p̂), roots by `fsolve` | **Confirmed.** Both roots of (16.136) agree with the 4 × 4 roots to ≤ 1e-15 at α = 2, 5, 8, 20 (and ν̂ = 0.3); at α = 8 the pressure root gives c/c₀ = 0.954064, k₂/k₁ = 0.062169, 0.67664 of the amplitude kept per wavelength |
| 4 | Womersley Crank–Nicolson order 2.00 | grids (n_r, n_t) = (20, 40) … (160, 320) against the Bessel profile at α = 3 and 8 | **Confirmed.** Observed order 1.999 (α = 3; pairwise 1.997, 1.999, 2.000) and 1.995 (α = 8; 1.990, 1.996, 1.999); finest error 2.5e-5 / 5.1e-5 of the peak speed |

## Validation table
| Concept / Eq. | Tier | fluidpy target | Evidence (V-levels) | Numbers | Label | Notes |
|---|---|---|---|---|---|---|
| PV loop, W = −∮p dV (ours) | C01 | `pv_loop_work`, `pv_loop_schematic(_work)` | V1 rectangle + orientation; V2 sympy ∮ over the legs; V3 polygon order | exact 1e-13; polygon error order 2.00 ± 0.15 | analytic, symbolic, converged | D01 |
| SV, EF, CO, CI, (16.9)–(16.11), units | N06–N14 | bookkeeping functions | V1 hand values; V2 pint | ≤ 1e-14 | analytic | ⅓ rule = time mean for a ⅓/⅔ square wave (1e-5) |
| Casson law (16.6)/(16.207) | C02 | `casson_stress`, `casson_shear_rate`, `casson_apparent_viscosity`, `differential_viscosity` | V1 inverse pair, secant, FD tangent; V2 sympy D02; V7 limits | 1e-12; μ_app − K_c² ∝ γ̇^{−0.5} (slope −0.5000) | analytic, symbolic | T14, T27 |
| Whitmore (16.7), yield stress (16.5) | N22, N19 | `whitmore_casson`, `yield_stress_correlation` | V1 | 1e-13 | analytic | |
| Womersley number (16.12) | C03 | `womersley_number`, `womersley_from_re_st`, `stokes_number` | V1 three routes (a√(ω/ν), ½√(Re St), √2a/δ_e via ch08); V2 sympy + pint | 1e-14 | analytic | T6 |
| Resistance (16.20), Poiseuille recap | C04, R04–R07 | `vascular_resistance`, `poiseuille_*`, `network_resistance` | V1 quadrature; V4 force balance τ_w2πaL = πa²Δp; V2 sympy D04 + units | 1e-12 / 1e-14 | analytic | S1 sign: τ > 0 |
| Entrance length (16.37) | N27 | `entrance_length` | V1; V5 Wikipedia (Bergman & Incropera, Çengel) | laminar 0.06 vs 0.0575: +4.3 %; turbulent 0.693 vs 1.359: factor 1.96 | analytic, benchmark | S31 reported |
| Two-layer model (16.21)–(16.34) | C05 | `two_layer_velocity`, `two_layer_flow_rates`, `apparent_viscosity_two_layer`, `apparent_viscosity_thin_layer` | V2 sympy D05–D07 (profiles from the stress law, fluxes, μ_app, thin-layer series); V1 quadrature; V7 limits | 1e-11; thin-layer error order 2.0 ± 0.15 | symbolic, analytic | **printed_16_30 fails** the integral (> 10 %) |
| FL closure + Charm–Kurland (16.35)–(16.36), (16.8) | N30, N36 | `fahraeus_lindqvist_model`, `charm_kurland_core_viscosity`, `cell_conservation` | V1 residual, cell balance; V7 trend | residual < 1e-12 | analytic (model, ours) | magnitude vs Pries: Open item 3 |
| Pries in-vitro viscosity | N29 | `pries_relative_viscosity` | V5 Secomb & Pries (2013) Eqs. (3)–(5) and their "minimum near 7 µm" | 1e-13; minimum at D = 6.85 µm | benchmark | read first-hand |
| Tube law (16.39)–(16.41), (16.55), (16.58) | C06 | `hoop_stress`, `tube_radius_hookean`, `tube_law_hookean`, `tube_law_linear`, `wall_bulge` | V2 sympy D08; V1 quadrature of (16.38), Hooke residual, round trips; V7 guard (T16), E → ∞ | 1e-12 | symbolic, analytic | |
| Distensible flux (16.42)–(16.46), compliance | N39–N45 | `elastic_tube_flow_rate`, `elastic_tube_pressure`, `elastic_tube_resistance`, `compliance_from_tube` | V2 sympy + pint (S4: printed a₀ is not a pressure); V1 solve_ivp of (16.43) (independent-route cross-check); V3 trapezoid/Simpson quadrature of the corrected (16.44) (loop 3); V1 FD of volume; V7 rigid limit | p(L) to 1e-6 p₁; profile 1e-9; rigid 1e-8; orders 1.999 / 3.989 | symbolic, analytic, converged | |
| Windkessel (16.49)–(16.52) | C07 | `windkessel_ode`, `windkessel_pressure`, `windkessel_periodic_state/p0`, `windkessel_summary`, `windkessel_terms` | V1 constant/zero inflow; V2 sympy D09–D10 (integrating factor = closed systolic form) + ODE residual of the periodic state; V1 ODE and quadrature from rest → periodic (cross-check); V3 Crank–Nicolson march, order 2.000 (loop 3); V4 ∫Q dt = ∫p/R dt; V5 Wikipedia diastolic form | 1e-10–1e-13 | analytic, converged, conserved | |
| Moens–Korteweg (16.56)–(16.63) | C08 | `moens_korteweg_speed`, `bramwell_hill_speed`, `youngs_modulus_from_pwv` | V1 Bramwell–Hill at A₀, inverse; V2 sympy linearisation D11 + cross-differentiation D12; V3 1-D model pulse speed; V5 Wikipedia form | 1e-9; speed within 0.5 % | analytic, symbolic, benchmark | S5: both momentum forms linearise alike (sympy), differ at 2nd order |
| d'Alembert, water hammer (16.65)–(16.68) | R09, R10 | `pulse_wave`, `velocity_from_pulse(_shapes)`, `pulse_state` | V1 wave-equation and ρu_t = −p_x residuals | 1e-4 / 1e-5 (FD) | analytic | |
| Nonlinear 1-D tube (16.53)–(16.55) | N50 | `tube1d_run` | V3 order 2.01 (grids 200/400/800 against translation); V3 pulse speed c(A₀) to 2.7e-4; V4 Σ A Δx to 1e-13 | | converged, conserved | **printed_16_54** sheds a reflected wave 4.0e-3 vs 1.7e-6 (corrected) |
| Womersley profile (16.72)–(16.88) | C09 | `womersley_velocity`, `womersley_F1`, `womersley_ode_residual`, `womersley_pde_solve`, `_small/_large_alpha` | V1 vs 30-digit mpmath (α = 0.7 … 45); V2 sympy D13 every line, D14, D15; V3 Crank–Nicolson order 2; V5 Wikipedia forms; V7 α² and 1/α error slopes, ch08 Stokes layer | 1e-11; orders 2.00, 2.00, 1.0 | analytic, symbolic, converged, benchmark | S6: residual = −G/μ |
| Flow, amplitude, phase (16.89)–(16.100) | C10 | `womersley_flow_rate(_harmonics)`, `womersley_wall_shear`, `womersley_XY`, `womersley_flow_ratio`, `womersley_phase_lag`, `womersley_F2(_series/_asymptote)` | V1 quadrature, FD wall shear, peak time = φ/ω; V2 sympy D16–D17 (+ mpmath ∫); V5 Wikipedia Q and τ_w; V7 series/asymptote slopes | 1e-10; slopes 4, 6, 3 (corrected) | analytic, symbolic, benchmark | S8 (τ_w sign) and **printed_16_96** (error slope 2, worse than the leading term) |
| Viscous elastic tube (16.116)–(16.137) | C11 | `elastic_tube_characteristic`, `elastic_tube_roots`, `elastic_tube_wave(_curve)`, `elastic_tube_determinant_numeric`, `_elastic_tube_modes` | V2 sympy D18–D21; V1 independent 4 × 4 determinant (cross-check); V3 ODE determinant order 6 in the series start (loop 3); V1 mode residuals; V7 1 − c/c₀ ∝ α⁻¹, k₂/k₁ ∝ α⁻¹ | roots 1e-15; slopes 1.0 ± 0.1; order 6.000 / 6.023 | symbolic, analytic, converged | T7: the other root is a fast wall mode |
| Murray (16.138)–(16.147) | C12 | `murray_*`, `bifurcation_cost/_point` | V2 sympy D22–D24 (stationarity, second derivative, cost split, Fermat balance); V1 grid search; V1 numeric branch point (cross-check); V7 symmetric 37.4673°; V5 Wikipedia cost | angles 1e-6 rad | symbolic, analytic, benchmark | **printed_16_140** fails stationarity and units; 70/30 split 25.706°/49.734° (= design) |
| Junction (16.148)–(16.163) | C13 | `junction_coefficients`, `junction_energy`, `junction_waveforms`, `junction_from_geometry`, `network_pulse_run`, admittance/impedance | V1 linear solve; V2 sympy D25–D26; V1 grid pulse (C = 1 transport exact; nx = 200 and 400 identical to 1e-12); V4 energy R² + (ΣY/Y₁)T² = 1 | 1e-14; grid ≤ 1e-3 | analytic, symbolic, conserved | Murray fork with equal c: ΣY/Y₁ = 2^{1/3}, R < 0 |
| Dean (16.164)–(16.182) | C14 | `dean_*`, `curved_resistance_ratio`, `DEAN_C4_DERIVED` | V2 sympy D27 steps 1–14 (1541/1045094400 → 1541/50400); V2 (16.174) sign from (16.167)–(16.168); V3 Navier–Stokes spectral solver (ψ₀, w₁, c₂, c₄) with a resolution study in N and M and a fit-truncation study (loop 3); V1 conventions | c₂ 2e-6; N = 21 vs 25: 3.2e-11; c₄ = 6471982981/542442700800 = 0.0119312 vs spectral 0.0119311 (9e-6); printed variant 0.0121130 fails (+1.53 %) | symbolic, converged | slip S34 (16.174) coded corrected, `printed_16_174` kept (loop 2) |
| Shapiro (16.185)–(16.206), Starling | C15 | `collapsible_*`, `shapiro_M`, `choking_flux`, `starling_state`, `flow_limitation_curve`, `wall_friction_term`, `friction_R_per_length` | V2 own sympy derivation of (16.204), (16.205) and (16.201); V1 closed-form choke, x(α) of the march, (16.205) along the march, frictionless throat (Bernoulli, α_t = 0.35); V2 + V1 throat condition α_t ≤ 0.25^{2/3} ⇔ H ≤ 0 (loop 3); V7 range guard α ≤ 1 (loop 3) | 1e-11 / 1e-8 | symbolic, analytic | **printed_16_201** non-zero when P depends on x; at S = 1 RHS(16.205) = −M·RHS(16.204); **L3-1: `collapsible_frictionless_alpha` raises for A₀/A₀,t ≳ 24 (failing test)** |
| Casson pipe (16.207)–(16.220) | C16 | `casson_pipe_velocity`, `casson_core_velocity`, `casson_pipe_flow_rate`, `casson_flux_ratio` | V2 sympy D31–D33 (r = s² substitution; (16.212)–(16.220) every line); V1 quadrature, annulus law by FD; V4 force balance; V7 limits | 1e-11 | symbolic, analytic, conserved | Q ∝ (1 − ξ)³ at the yield end; S26 |
| Water potential (16.221), xylem, phloem (16.222) | C17, N110–N113 | `water_potential(_terms)`, `osmotic_pressure`, `xylem_column`, `sieve_tube_pressure_drop`, … | V1 terms, column, still water ψ uniform, sieve limits; V2 pint + book-unit form | 1e-12 | analytic | S30: printed ρgz not ÷10⁶ |
| Power-law, Herschel–Bulkley | S01 pointer | `power_law_*`, `herschel_bulkley_pipe_flow_rate` | V1 quadrature, Buckingham–Reiner | 1e-10 | analytic | |

Counts of test functions (loop 3, 99 tests, one marked slow) — by the **first** level in the name: V1 36 · V2 27 · V3 9 ·
V4 4 · V5 4 · V6 1 · V7 15; by **any** level named: V1 36 · V2 27 · V3 10 · V4 4 · V5 5 · V6 1 · V7 17. Every test
named V3 now asserts a refinement order or a resolution study (loop 2 counted V3 12 by any level, 6 of which had none).
Loop 2's counts: V1 33 · V2 26 · V3 11 · V4 4 · V5 4 · V6 1 · V7 12 (first level).
(Loop 1's line "V1 30 · V2 25 · … · V7 13, 89 tests" used a different counting rule; these are recounted by script.)

## Functions used by the notebook and explainers (design Part C)
All 165 names of C.1–C.10 exist and each is called in `tests/test_ch16.py` — enforced by
`test_contract_V7_every_design_part_c_function_is_exercised_here` (parses Part C; fails on any missing or uncalled name).

## Convergence studies
| scheme | grids | observed order | design order |
|---|---|---|---|
| PV-loop polygon (trapezoid) | n = 25 … 400 | 2.00 | 2 |
| Womersley Crank–Nicolson (Δr, Δt halved together) | (20, 40) … (160, 320) | 1.999 (α = 3), 1.995 (α = 8) | 2 |
| 1-D tube MacCormack + 4th-difference damping | 200, 400, 800 | 2.008 (pairwise 2.009, 2.007) | 2 |
| quasi-steady limit (16.87) error | α = 0.05 … 0.4 | 2.0 | 2 (α²) |
| large-α limit (16.88) error | α = 20 … 160 | 1.0 | 1 (α⁻¹) |
| F₂ series (16.93), order 2 / 4 | α = 0.02 … 0.16 | 4.0 / 6.0 | 4 / 6 |
| F₂ asymptote corrected / printed (16.96) | α = 40 … 320 | 3.0 / 2.0 | 3 / (2) |
| thin-layer (16.34) relative error | δ/a = 0.0025 … 0.02 | 2.0 | 2 |
| elastic tube 1 − c/c₀, k₂/k₁ | α = 25 … 200 | 1.0, 1.0 | 1 |
| Dean spectral solver (ψ₀ against κψ₀) | K_lec = 0.25 … 1 | 2.0 | 2 (κ²) |
| Dean spectral solver, resolution (loop 3) | N = 9 … 25 (M = 16), M = 8 … 16 | geometric: 1.1e-4, 4.2e-7, 1.0e-8, 3.2e-11 | spectral |
| Dean c₄ fit, number of fitted terms (loop 3) | 3 … 6 terms | 4.7e-3, 2.3e-4, 9.4e-6, 3.3e-7 relative to the exact rational | (series truncation) |
| corrected (16.44) quadrature, trapezoid / Simpson (loop 3) | n = 8 … 64 | 1.999 / 3.989 | 2 / 4 |
| Windkessel Crank–Nicolson march (loop 3) | Δt = T/40 … T/320 | 2.000 | 2 |
| ODE determinant, series start s₀ (loop 3) | s₀ = 0.2 … 0.025 | 6.000 (α = 2), 6.023 (α = 8) | 6 |

## Conservation / invariant residuals
Windkessel ∫Q dt − ∫p/R dt over a beat: 1e-10 relative · 1-D tube Σ A Δx: < 1e-13 · junction energy R² + (ΣY/Y₁)T² − 1:
< 1e-14 (200 random networks) · Poiseuille and Casson force balance τ_w 2πaL − πa²Δp: < 1e-14 · Starling plateau equals the
closed-form choke to 2e-16.

## Benchmarks used
| value / form | our value | source | verified |
|---|---|---|---|
| Pries in-vitro viscosity, Eqs. (3)–(5); minimum "down to about 7 µm" | identical form (1e-13); minimum at 6.85 µm | Secomb & Pries, C. R. Physique 14, 470 (2013), doi:10.1016/j.crhy.2013.04.002 | 2026-10-14 |
| Womersley u, Q, τ_w (Λ = α n^{1/2} i^{3/2}) | identical (1e-11) with P′ = −G; their τ_w = +μ∂u/∂r, ours −μ∂u/∂r | Wikipedia "Pulsatile flow" | 2026-10-14 |
| Moens–Korteweg √(Eh/(2rρ)), r = radius | identical | Wikipedia "Moens–Korteweg equation" | 2026-10-14 |
| Windkessel I = P/R + C dP/dt; diastolic decay | identical (1e-12) | Wikipedia "Windkessel effect" | 2026-10-14 |
| Murray cost and cube law | identical | Wikipedia "Murray's law" | 2026-10-14 |
| entrance length 0.0575 Re D, 1.359 D Re^{1/4} | (16.37): +4.3 % laminar, ×0.51 turbulent | Wikipedia "Entrance length (fluid dynamics)" | 2026-10-14 |
| Dean equations, Navier–Stokes form (no number) | solved here: c₂ = 0.0305754, c₄ = 0.0119311 | Imperial College M4A33 notes | 2026-10-14 |
Full citations: `reference/ch16/SOURCES.md`; values: `reference/ch16/benchmarks.json` (written by `make_refs.py`).

## Numbers from the text (book vs ours) — private values redacted to relative errors
All in `test_book_V6_worked_numbers` (passes): the cardiac-output example 0 %; the PRU example — the book rounds the exact value to
one decimal (3.8 % from the exact value, rounding); the Womersley-number estimate 0.5 % (the book's ω is rounded; with
2π·HR/60 it rounds to the printed integer); the aortic Reynolds number 0.37 %; the α = 8 example: Y within 0.5 %, X 3.2 %
below the print (slip S11, confirmed), |Q|/|Q₀| within rounding; the Casson minimum pressure drop: the printed mmHg value
is right, the printed dyn/cm² value is 1000× too small (S26, confirmed); the xylem Reynolds number needs water's ν
(S28, confirmed); the yield-stress range at H = 0.45 lies below the correlation's (S29, confirmed); the book's R in
cm³ bar/(mol K) is 4e-5 below CODATA; 133 Pa/mmHg is 0.24 % below the exact 133.322; the Dean coefficient 0.0306 is
0.08 % above 1541/50400; the printed resistance coefficient of (16.182) is consistent with inverting (16.178) (6.5e-5).

## Figures reproduced with our code (`tests/ch16_verify_figures.py` → `outputs/ch16/verify/`, local)
| figure | png | verdict |
|---|---|---|
| Fig. 16.6 analogue (PV loop) | `ch16_heart_blood.png` (a) | counter-clockwise loop through A (ESV, low p) → B → C → D, work 1.049 J > 0, isovolumic legs vertical |
| Figs. 16.7–16.8 analogues (Casson) | `ch16_heart_blood.png` (b, c) | √(τ/μ_p) is a straight line in √γ̇ with intercept 2.0; μ_app falls monotonically onto the K_c² asymptote, the tangent viscosity lies below the secant |
| Fig. 16.21 analogue (Casson plug) | `ch16_heart_blood.png` (d) | flat plug of radius a_c at the centre, zero at the wall, below the Poiseuille profile with μ = K_c² everywhere |
| Fig. 16.12 analogue + FL effect | `ch16_heart_blood.png` (e, f) | two-layer profile has its slope kink at a − δ; both μ_app(D) curves rise monotonically from 20 to 500 µm, the closure 20–35 % below Pries (Open item 3) |
| Fig. 16.14 analogue (Womersley) | `ch16_womersley_windkessel_waves.png` (a–c) | α = 2.5 near-parabolic with reversal, α = 5 near-wall overshoot, α = 12 flat core with thin wall layers; profiles vanish at r = ±a and are symmetric |
| amplitude and lag against α | `…waves.png` (d) | |Q|/|Q₀| = 1 at small α and follows 8/α² at large α; φ rises from 0 to 90° through α ≈ 3 |
| Windkessel beat | `…waves.png` (e) | periodic pressure 73/109 mmHg (ours), the run from rest converges onto it, the true mean R·CO lies above the ⅓-rule estimate |
| viscous elastic tube | `…waves.png` (f) | c/c₀ and the amplitude kept per wavelength both rise to 1 as α grows |
| Fig. 16.15 cost | `ch16_branches_bends_collapse_plants.png` (a) | total cost minimum at the dashed optimum radius, where pumping is half the upkeep |
| Figs. 16.17–16.18 analogues (Dean) | `…plants.png` (b, c) | two antisymmetric counter-rotating cells; the axial maximum sits 0.187a toward the outer bend at K = 576 |
| Fig. 16.19 analogue (tube law) | `…plants.png` (d) | (loop 3, α ≤ 1 only) both forms rise steeply as α → 0.15 (−P ≈ 16–17), the shapiro curve reaches −P = 0 and the pure curve −P = 1 exactly at the dotted α = 1 line where the fitted range ends |
| Starling resistor | `…plants.png` (e) | Q rises as p₂ falls and plateaus at 20.0 mL/s for −(p₂ − p_e) > 946 Pa (dashed) |
| §16.4 xylem column | `…plants.png` (f) | p falls linearly with height faster than ψ (gravity + drag), ψ falls by the viscous drop only, gravity term grows linearly |

## Derivation storyboards (design Part F) — intermediate lines
Every ★★★ storyboard's lines were re-run as sympy (tests named `*_V2_derivation*`): D13 (each of (16.76)–(16.82), with
the corrected right side −G/μ), D19–D21 (the rows (16.134)–(16.135) come out as −(16.124) and −(i/ω)(16.125); their
determinant is (const/X)·(16.136)), D24 (the force balance and the three cosines), D27 (steps 1–14, including the
intermediate ∫r f g dr = 1541/1045094400), D30 (16.204)–(16.205) from the differentials, D33 (16.212)–(16.220). The
"expect" numbers of D24 (25.706°, 49.734° for a 70/30 split) and D27 (ψ₀ max 0.00757 at r = 0.429, axial maximum 0.187a
toward the outer bend and flux ratio 0.9814 at K = 576) are reproduced. No wrong intermediate line was found in Part F;
the D27 "Result" line keeps the printed 0.0120 for the K⁴ coefficient; the code now derives 0.0119312 (S34) — **for the
designer/notebook**: if D27 or the explainer shows the K⁴ coefficient, show 0.0120 as printed and 0.0119312 as the
consistent value, with (16.174) written with the corrected + sign (and a note that the printed − sign gives 0.0121130).

## Label audit (docstring `Validation:` levels against the tests that call each function)
After tightening (sympy expressions now also evaluated against `two_layer_velocity`, `casson_pipe_velocity`,
`casson_core_velocity`, `tube_law_linear`, `elastic_tube_flow_rate`, `elastic_tube_characteristic`, `dean_flux_ratio`,
`casson_constants`, `casson_stress_plasma_form`, `womersley_profile_complex`), 21 functions still name a level that differs
from the level in the name of the test that calls them, e.g. a "V1" claim checked inside a test named V3/V4 or a "V7"
limit checked inside a V1 test (`HalfSineEjection`, `ejection_half_sine`, `windkessel_periodic_p0(_general)`,
`windkessel_terms`, `pv_loop_schematic`, `network_resistance`, `murray_theta`, `murray_symmetric_angle`,
`casson_min_pressure_drop`, `casson_core_velocity`, `power_law_*`, `herschel_bulkley_pipe_flow_rate`,
`yield_stress_correlation`, `elastic_tube_wave_curve`, `water_potential_book_units`, `womersley_ode_residual`); the evidence
itself exists in each case. `tube1d_run` is exercised through the helper `_tube1d_case` (V3, V4). No over-claim of
substance was found.

## Deviations & justifications (`# DEVIATION` in the code)
All checked against the page images and an independent test: S3 (16.30), S4 (16.44)–(16.45), S5 (16.54), S6 (16.76),
S8 (16.92), S10 (16.96), S15 (16.140), S21 (f′_L constant), S23 (16.201), S26, S27, S30; the 1-D tube scheme and damping
(ours); the Womersley periodic-map solve (ours); the Starling α-integration (ours; confirmed by the closed form);
**S34 (16.174)** — `# DEVIATION` in `dean_series_sympy` (the advection term) and `dean_flux_ratio` (c₄), loop 2; confirmed
by two independent routes (sympy elimination of p₁ from (16.167)–(16.168) with (16.172); a spectral Newton solve of the
Navier–Stokes form of the Dean equations) and, through `printed_16_174=True`, shown to be the only difference.

## Recorded for the analyst/designer (no code change; not verdict-blocking)
- **S34** should be added to analysis §9 (row: (16.174), printed "(2/κ)∇⁴ψ − (1/r)(ψ_r∂_α − ψ_α∂_r)∇²ψ = −2w(…)",
  correct "+ (1/r)(…)", slip, changes physics: K⁴ flux coefficient only) and the §8 B6 note updated (our 0.0119312 is
  0.16 % below the often-quoted 0.01195, which was not read first-hand). `analysis/ch16.md` currently has no S34 row.
- **Analysis §9 row S3 sign note (correction).** The row says "printed minus correct = πΔp a²δ(δ − 2a)/(4μ_pL)". Sympy
  (`test_two_layer_V2_derivation_profiles_fluxes_viscosity`) gives **printed − correct = +πΔp a²δ(2a − δ)/(4μ_pL)**,
  i.e. the printed Q_c of (16.30) is the larger by that amount; the analysis row has the sign reversed (its
  "correct" expression for Q_c is right). `analysis/ch16.md` still carries the old sign as of this pass.
- Minor wording (scripts, implementer): `scripts/ch16_dean.py` still prints "(printed 0.0120: our exact solution of
  (16.173)–(16.176) is … % off — reported, not resolved)"; with S34 the printed number is now −0.57 % off and the cause
  is resolved (the line should name S34). Output is correct; the comment is stale.

## Findings for the implementer (fix loop 1 — all resolved in loop 2)
1. **Must-fix (physics) — new slip in (16.174).** Eliminating p₁ from (16.167)–(16.168) with (16.172) gives
   $\dfrac2\kappa\nabla_1^4\psi+\dfrac1r\left(\dfrac{\partial\psi}{\partial r}\dfrac{\partial}{\partial\alpha}-\dfrac{\partial\psi}{\partial\alpha}\dfrac{\partial}{\partial r}\right)\nabla_1^2\psi=-2w\left(\sin\alpha\dfrac{\partial w}{\partial r}+\dfrac{\cos\alpha}r\dfrac{\partial w}{\partial\alpha}\right)$ — the book prints a minus before the bracket (test
   `test_dean_vorticity_V2_printed_16_174_nonlinear_sign_is_inconsistent`). Only the order-K⁴ flux coefficient is affected.
   Fix: in `dean_series_sympy` subtract (not add) the ψ-advection term of the vorticity equation (keep `printed_16_174=True`
   for the printed sign), set `DEAN_C4_DERIVED = Fraction(6471982981, 542442700800)` (0.0119312), update the docstrings of
   `dean_flux_ratio` (exact_c4: "0.6 % below the printed 0.0120") and `dean_series_sympy`, add S34 to `book_slips()` (kind
   "slip", changes_physics True) and to the module docstring. The failing test
   `test_dean_flux_V3_order4_coefficient_is_the_navier_stokes_value` then passes unchanged.
2. **Should-fix (robustness).** `windkessel_periodic_state` evaluates both branches of `np.where`; for very small RK the
   diastolic exponential overflows and warns (`RuntimeWarning: overflow encountered in exp`, seen with K = 1e-14).
   Evaluate the diastolic branch only where s ≥ T_s.
3. **For the analyst/designer (no code change).** Analysis §9 row S3 states "printed minus correct = πΔp a²δ(δ − 2a)/(4μ_pL)";
   sympy gives **+πΔp a²δ(2a − δ)/(4μ_pL)** (the printed Q_c is the larger) — the sign in the analysis row is wrong.
   Design C.13 difference (1) is settled in favour of the α-integration; E9's parity tolerance can drop to rtol 1e-10.

## Open items
0. **L3-1 (blocks PASS):** `collapsible_frictionless_alpha` raises `ValueError` when the critical area ratio α* < 1e-6
   (A₀/A₀,throat ≳ 24 for the demo throat); failing test `test_collapsible_V1_frictionless_alpha_on_wide_upstream_sections`;
   one-line closed-form fix verified in a scratch copy (see Loop 3).
1. Dean (1928) itself was not reached; the K⁴ coefficient therefore rests on two independent routes of ours (sympy series
   of the corrected equations and a spectral Navier–Stokes solution, agreeing to 9e-6), not on a cited value. The
   often-quoted 0.01195 (analysis §8 B6, not read first-hand) would be 0.16 % from ours. Label for c₄: `symbolic` +
   `converged` (no V5).
2. No public numeric table was reached for the viscous elastic-tube wave (Womersley 1955/1957; analysis B11): C11 rests
   on sympy (D18–D21) and an independent 4 × 4 determinant (V1 cross-check), both of which agree to 1e-15, plus the order-6 study of the ODE determinant (V3, loop 3).
3. The Fåhræus–Lindqvist closure (two layers + Charm–Kurland, δ = 2.5 µm, ours — trap T19) reproduces the direction of the
   effect but is 20 % (D = 20 µm) to 35 % (D = 300 µm) below the cited Pries in-vitro curve, because Charm–Kurland at
   310 K gives a whole-blood viscosity ratio near 2 where Pries's high-shear data give 3.2; label `qualitative` for the
   magnitude. It also has no root below D ≈ 15 µm at δ = 2.5 µm (documented `ValueError`).
4. Slip S21 (A_e in (16.193)) stays unresolved without Shapiro (1977a); coded as a constant f′_L keyword.
5. The Wikipedia Womersley-number page's human values do not state the heart rate or length scale; not used.
