# Chapter 10 verification — Computational Fluid Dynamics
2026-09-30 · repository state `675fcfc` (implement phase) + uncommitted verifier files (`tests/test_ch10.py`,
`tests/ch10_verify_figures.py`, `reference/ch10/{make_refs.py (cross-check added), SOURCES.md, verify_runs.py,
cavity_mck_re100_n128.csv}`) · verifier: `math-verifier` (loop 1)

## Verdict: PASS (after fix loop 1 and the derivation-review fixes)

### Review fixes (after `reports/ch10_review.md`, 2026-09-30)
The implementer fixed the review's three code Must-fixes; tests were added, none loosened.
- **`tests/test_ch10.py`: 340 items, 340 passed, 0 failed (103 s, slow tests included); 165 test functions.** Labels now
  V1 76 · V2 36 · V3 28 · V4 7 · V5 5 · V6 2 · V7 20 (+ 3 smoke).
- **Must-fix 1, `mac.primary_vortex_centre` sign:** the fitted ψ_min was pushed *above* the discrete minimum (the tables were
  written by the same code, so the old tests could not see it). New tests:
  `test_primary_vortex_centre_V1_exact_on_a_separable_paraboloid[20, 33]` (vertex (0.6137, 0.7311) and ψ_min = −0.1 returned
  to 1e-12/1e-14; fitted ≤ grid min) and `test_primary_vortex_centre_V7_fitted_minimum_below_the_discrete_one_on_cavity_runs`
  (live 16², 32²; the regenerated tables ψ_min 32/64/128² = −0.101363, −0.102986, −0.103392 → **p = 2.00, f_ext −0.103528,
  GCI_fine 0.164 %**; was p = 1.94). Vortex positions were unaffected (the x, y corrections were right).
- **Must-fix 2, `fd.rod_heating_exact` (10.199) stopping rule:** it stopped where sin((2m − 1)πx/L) = 0 (x = 0.2, t = 0.001 gave
  −0.110). Now it stops on the term envelope. New tests against the independent erfc image solution:
  `test_rod_heating_V1_series_agrees_with_the_image_solution_where_sines_vanish` (x = 0.2, t = 0.001: 7.7442e-06; x = 1/3,
  t = 0.01: 0.0184246 — both equal to the image sum to 1e-11; the "0.018422" quoted to me was a one-term short-time
  approximation, the full image sum agrees exactly) and `test_rod_heating_V7_odd_cell_grids_nonnegative_and_equal_to_images`
  (n = 5, 15, 21, 45 at t = 0.001, 0.01, 0.1: ≥ 0 and = images to 1e-11).
- **Must-fix 3, `maccormack.cavity_maccormack` time step:** new default `dt_rule="additive"` (`maccormack_dt_additive`, a
  documented DEVIATION from (10.110); the book's rule stays as `dt_rule="tap"`) and a `check_stability=True` guard. New tests:
  `test_cavity_maccormack_V7_low_Re_runs_stay_finite_with_the_additive_step[5.0-64, 2.0-32]` (finite, |u| ≤ U to t = 1; the
  old inviscid bracket produces NaN within 50 steps in the same runs) and
  `test_cavity_maccormack_V7_too_large_step_raises_and_additive_limits` (dt = 0.1 raises; unknown rule raises; → inviscid bracket
  as μ → 0 and → the Heun diffusion limit as c, u, v → 0). Regenerated MacCormack tables: max deviation from Ghia **4.14 %,
  1.38 %, 0.441 %** (32/64/128², Re = 100), **5.86 %** (Re = 400, 64²); mass drift −0.68/−0.56/−0.49 %.
- Follow-ups: `mac.run`/`step` lid default is now None (grid walls) — the test passes `lid=1.0` so the convective limit
  0.02 < 0.05 raises; stale comments refreshed (MacCormack percentages, F2/F3 marked fixed, block p ≈ 0.43).
  `stability_verdict` with diffusion: `test_stability_verdict_V1_reasons_and_flags_with_diffusion_on_a_C_beta_grid` (26 × 15 grid;
  stable ⇔ |C| + 2β ≤ 1 (upwind), C² + 2β ≤ 1 (Lax–Wendroff) — closed forms derived by the verifier: for upwind |G|² is quadratic in q = 1 − cos θ with a positive leading coefficient, so the ends q → 0, 2 decide;
  for LW, a²−C²−2β = a(a − 1) with a = C² + 2β; reason strings exact). `numerical_diffusivity` without C for "upwind"/"ftcs":
  `test_numerical_diffusivity_V7_unsteady_schemes_require_the_courant_number`.
- **Mutants on the changed functions: 9 planted, 9 killed** — M34 old ψ_min sign (killed by the C15 ψ_min test, and the
  paraboloid test), M35 old stopping rule (image test), M36 additive Δt without the viscous rate (low-Re run), M37 guard
  disabled (raise test), M38 old verdict reason (C–β grid), plus M18, M31, M32, M33 on the loop-1 fixes. Total over the chapter:
  38 planted, 38 killed.
- `fd`, `mac`, `maccormack` are used by no earlier chapter (grep), so the full-tree result stands. The design/analysis wording
  items stay open for the notebook-builder (Open item 1).

### Loop 1 (re-verification after the implementer's fixes, 2026-09-30)
The implementer changed `fluidpy/core/mac.py` (`_lap1d` built link by link from real neighbours), `fluidpy/core/fd.py`
(`advect_periodic` reference = integer `np.roll` shift when k·C is whole, else a rounded argument; `grid_convergence_index`
returns NaN f_ext/GCIs when rᵖ − 1 = 0, p = the observed 0.0), wording in `fd.py`, `maccormack.py`, `ch10.book_slips()` (R3)
and a non-asymptotic line in `scripts/ch10_block.py`. `tests/test_ch10.py` and every tolerance are unchanged.
- **`tests/test_ch10.py`: 327 passed, 0 failed (120 s, slow tests included)** — re-run by the verifier.
- Spot checks: `pressure_poisson_matrix(MacGrid(3, 1, 3, 1))` = [[−1, 1, 0], [1, −2, 1], [0, 1, −1]] (the D20 pipe); square pulse
  at C = 0.8 now reports rms_error 0.1442 against the corrected reference; `grid_convergence_index(1.0, 1.1, 1.0)` →
  p = 0.0, f_ext = GCIs = NaN, monotone False (no exception).
- **Mutants re-run on the fixed code: 33 planted, 33 killed, none hung** (baseline 324 passed; M18 re-targeted at the new
  `_lap1d` link weight, and three new mutants on the fixed code: M31 reference shifted one node too far → killed by
  `test_advect_periodic_V1_reported_error_is_zero_for_an_exact_shift`; M32 p = 0 guard removed → killed by
  `test_gci_V7_equal_and_opposite_differences_return_nan_not_an_exception`; M33 wall link reflected onto the inner neighbour →
  killed by `test_projection_V4_divergence_removed_to_round_off`).
- Docstring corrections confirmed in the code: R3 downwind wording (`fd.py` module docstring, `steady_cd_fd`, `book_slips()` R3
  evaluator); the stretched grid "shrinks, does not remove" the wiggles (`steady_cd_fd`); `lax_demo` (|G| = 1.04, step 690);
  MacCormack cavity "not mass-conserving" and the low-Mach caveat (`maccormack.py` module docstring, `cavity_maccormack`,
  `wc_taylor_green`); block study labelled *qualitative / non-asymptotic* (`block_channel`), and `ch10_block.py` prints it.
- `fd`, `mac`, `maccormack` are used by no earlier chapter (grep), so the loop-0 full-tree run (ch01–ch09 all green) stands;
  `tools/check_public.py` OK.

Loop-0 record follows (the three defects below are fixed).

Suite results at loop 0 (2026-09-30): **`tests/test_ch10.py`: 327 items — 316 passed, 11 failed** (non-slow part 78 s; 3 slow tests ≈ 2 min,
all pass). **Full tree `pytest tests --ignore=tests/test_machinery.py` (ch01–ch10, slow tests included): 1427 passed, 11 failed in
18 min 06 s** — the 11 failures are exactly the ch10 items F1–F3 below; ch01–ch09 show no regression. **`tests/test_machinery.py`:
9 passed (20 s).** All 11 runnable `scripts/ch10_*.py` run headless, exit 0 (`ch10_common.py` is a helper). Evidence gates: every computable CORE row has ≥ 2
independent levels (15/15); every coded NOTE row ≥ 1; every design Part C name is exported, documented and exercised; 23/23
derivations re-derived; 30/30 planted wrong variants killed; nothing `unverified`. The only blocker is the failing tests.

**Loop 0 was FAIL — 3 defects in `fluidpy/`, 11 failing test items (all fixed in loop 1); every other test passed with unchanged tolerances.** No physics
transcription error was found; all three defects are local and cheap to fix.

| # | failing test(s) | metric | hypothesis | fix (for the implementer) |
|---|---|---|---|---|
| F1 | `test_pressure_poisson_V1_single_cell_direction_rows_sum_to_zero_D20_pipe`; `test_pressure_poisson_V1_constant_null_vector_on_thin_grids[4-1-per0, 1-4-per1, 4-1-per2, 1-4-per3, 2-4-per4, 4-2-per5]` (7 items) | `MAC.pressure_poisson_matrix(MacGrid(3, 1, 3.0, 1.0))` = [[−2, 1, 0], [1, −3, 1], [0, 1, −2]] instead of D20's [[−1, 1, 0], [1, −2, 1], [0, 1, −1]]; max \|A·1\| = 1/h² (must be 0) on every grid with one cell in a direction (walls or periodic) and on a periodic direction with two cells | `_lap1d(n, h, periodic)` builds tridiag(1, −2, 1) and then *assigns* the end entries: with n = 1 both wall assignments hit the same entry (−1 instead of 0 — a Dirichlet ghost), periodic n = 1 gives −1, periodic n = 2 overwrites the wrap entry with 1 where the two neighbours coincide (needs 2) | build the 1-D operator from the neighbour count (diag = −number of real neighbours, off-diagonals *added* with `+=`); for walls and n = 1 the entry is 0. The D20 tiny example ("3-cell pipe", notebook and E6) is exactly this grid, and `solve_pressure` on it returns a wrong p |
| F2 | `test_advect_periodic_V1_reported_error_is_zero_for_an_exact_shift[1.0, 0.3, 2.0]` (3 items) | `advect_periodic("square", C = 1, N = 50, n_rev, "upwind")` reports `rms_error` 0.1414 (n_rev 1), 0.1414 (0.3), 0.2 (2) although the computed T equals the exactly shifted T0 to 1e-15 (`test_upwind_V1_design_example_and_exact_shift_at_C_equal_one` passes) | the reference profile is evaluated as `f(x − k·C/N)` in floating point: (0.2 − 1.0) % 1 = 0.19999999999999996 < 0.2 moves the square's edge by one node — the *reference* is wrong, not the scheme. Every square-pulse `rms_error` (E3 end-of-run card, F3 figure, A2) carries this edge bias | evaluate the reference by an integer node shift when k·C is an integer, or round the argument (e.g. `np.round((x − s) % 1, 12)`) before the edge comparison |
| F3 | `test_gci_V7_equal_and_opposite_differences_return_nan_not_an_exception` (1 item) | `grid_convergence_index(1.0, 1.1, 1.0)` raises `ZeroDivisionError` | \|(f₃ − f₂)/(f₂ − f₁)\| = 1 ⇒ p = 0 ⇒ rᵖ − 1 = 0 is not guarded (d = 0 is); reachable from C15/E7 when a user picks such grids | return NaN for f_ext and the GCIs when rᵖ − 1 = 0 |

## Environment
Python 3.11.5 · numpy 2.4.6 · scipy 1.17.1 · sympy 1.14.0 · pint 0.25.3 · matplotlib 3.11.2 · pytest 9.1.1 (Windows 11, `.venv`;
`progress.json → environment` checked 2026-09-12 — unchanged). Commands: `.venv/Scripts/python.exe -m pytest tests/test_ch10.py -q
-p no:cacheprovider` (non-slow part ≈ 80 s; the 3 `slow` tests add ≈ 2 min: the 11 scripts headless, a live Δx = 1/8 block run
reproducing the committed table to 6 s.f., a live 32² MacCormack cavity); figures `.venv/Scripts/python.exe
tests/ch10_verify_figures.py`; references `reference/ch10/make_refs.py` (Ghia/Hou, cross-check now asserted) and
`reference/ch10/verify_runs.py` (MacCormack 128², ≈ 6.5 min). The 11 scripts were also run by hand with `--no-show`: all exit 0
(block 2 s, cavity 2 s, cylinder 6 s, drawings 2 s, fem1d 4 s, lax 3 s, projection 21 s, stability 5 s, steady_cd 3 s, stencils
3 s, tables 2 s — the heavy runs hit their caches). `tests/book_values_ch10.json` is valid JSON and git-ignored
(`.gitignore:25`); this report quotes relative differences only.

Test inventory: 157 test functions, 327 collected items (153 of them the parametrised Part C smoke test). Labels (a test counts
once per level in its definition comment): **V1 72 · V2 36 · V3 28 · V4 7 · V5 5 · V6 2 · V7 15**, plus 3 smoke tests. 23/23
derivations re-derived. `qualitative`: the block C_D grid study and the confined-cylinder Strouhal number (Open items).
`unverified`: none.

## Validation table (CORE rows; every computable row has ≥ 2 independent levels)
| Concept / Eq. | Tier | fluidpy target | Evidence (V-levels) | Numbers | Label | Notes |
|---|---|---|---|---|---|---|
| C01 stencils (10.4)–(10.8) | CORE | `fd_weights`, `fd_leading_error`, `stencil_taylor_coefficients`, `stencil_derivative/error`, `fd_derivative`, `taylor_table_sympy`, `mixed_derivative_onesided` | V1 exact rational weights (10.6), (10.7), p. 450, (10.149); exact on polynomials; (10.150) exact on bi-quadratics · V2 leading coefficient of 6 stencils from the series of e^{ax} (independent of the Taylor table); D01 all 9 lines · V3 orders on sin · V6 | orders 1.004 / 0.996 / 2.000 / 2.000 / 1.977 (forward, backward, centred, second derivative, one-sided); design numbers (−4.294e-2 …) to 5e-6 | analytic, symbolic, converged | round-off floor ε/h and ε/h² visible (figure) |
| C02 FTCS / BTCS (10.9)–(10.13) | CORE | `transport_1d_step`, `solve_transport_1d`, `ftcs_coefficients`, `cfl_number`, `diffusion_number` | V1 five-point example (1e-15); from-scratch `np.roll` update (1e-14); u = 0 ≡ ch01 `ftcs_diffusion_1d` (1e-14); BTCS residual of (10.13) 1e-13; Neumann ghost exact on linear T (1e-12) · V4 periodic sum · V3 · V2 D02 + dimensions · V7 guard/cap | FTCS 1.99 [1.85, 2.09, 2.00]; CN 2.03; BTCS 1.98 [1.95, 1.99, 2.00] on n = 80…640 | analytic, conserved, converged | BTCS default grids pre-asymptotic (1.80) — doc correction |
| C03 consistency, truncation error (10.14)–(10.17) | CORE | `truncation_error_sympy`, `truncation_terms`, `convergence_study`, `error_norm` | V2 code returns (10.17) exactly, BTCS (−Δt/2 T_tt), upwind, steady modified equations; D03 re-derived on e^{ax+bt} incl. step 11 and the 0.0406 check · V3 remainder (measured − (10.17)) · V1 norms | remainder order 4.00 [3.99, 4.00, 4.00] with Δt ∝ Δx²; one-step residual within 5 % of (10.17) | symbolic, converged | |
| C04 von Neumann (10.18)–(10.28) | CORE | `amplification_factor`, `amplification_modulus`, `ftcs_amplification_modulus2`, `amplification_curve`, `max_amplification`, `worst_theta`, `is_von_neumann_stable`, `ftcs_stable`, `stability_verdict`, `propagate_error`, `fourier_mode` | V1 \|G(10.24)\|² = (10.26) (1e-12, 20 random (α, β)); one Fourier mode grows by \|G\|ⁿ (1e-9) for FTCS, BTCS, CN, upwind ± u, Lax–Wendroff; Noye closed form = brute force (0 mismatches on 31 × 29) · V7 both edges of (10.27) · V2 D04–D07 | \|G\|max = 1.04 at θ = π for β = 0.51; (0.4, 0.2) fails at θ < π/2 | analytic | |
| C05 upwind, CFL, Lax (10.29)–(10.31), (10.199) | CORE | `transport_1d_step("upwind")`, `advect_periodic`, `numerical_diffusivity`, `phase_error`, `lax_demo`, `rod_heating_exact`, `cfl_time_step` | V1 C = 1 exact shift (1e-15, upwind/MacCormack/LW); D08 modulus (1e-13); upwind Gaussian 10× closer to the D_num = \|u\|Δx(1 − C)/2 solution than to pure translation; rod series solves the heat equation · V7 CFL edges, R11 printed stencil unstable for u < 0, u → −u mirror (1e-15), Lax demo · V3 · V2 D08 | upwind 0.947 [0.91, 0.95, 0.976] on N = 160…1280; FTCS → rod series 1.98; β = 0.51 blows at step 690 | analytic, converged | **F2** (reported rms_error of the square pulse) |
| C06 weak form (10.32)–(10.38) | CORE | `weak_form_sympy`, `weak_to_strong_sympy`, `weak_residual`, `bilinear_form` | V2 engines vanish; D09, D10 re-derived · V1 exact steady solution with a flux datum passes three test functions (< 1e-9), a wrong flux fails (> 1e-3) · V3 natural slope approached, T(0) exact | slope order 0.95 [0.92, 0.96, 0.98] | symbolic, analytic | |
| C07 Galerkin (10.39)–(10.63) | CORE | `galerkin_equations_sympy`, `hat`, `hat_basis`, `interpolate`, `solve_steady`, `solve_transport`, `interior_stencil` | V2 check_interior 0, M tridiagonal, F row; D11 (my element integrals = engine = `assemble_1d` numbers), D12 · V1 steady FE = centred FD (≤ 1e-12, three (n, R)); partition of unity; θ-scheme = method of lines (1e-4) · V3 | Neumann steady 2.00 [2.005, 2.001, 2.000]; transient half rod vs series 2.00 | analytic, symbolic, converged | |
| C08 element matrices and assembly (10.64)–(10.78) | CORE | `element_matrices_linear`, `element_integrals`, `element_force`, `connectivity`, `assemble_1d`, `assembly_trace`, `shape_slopes`, `map_to_element/parent`, `parent_shapes` | V1 closed form = Gauss (1e-13); assembled interior rows = (10.63)·h (1e-14); trace rebuilds `assemble_1d`; a(N_A, N_B) = K_AB · V2 D13, D14 with numeric parity · R1, R2 printed variants fail | design numbers exact | analytic, symbolic | |
| C09 cell Péclet (10.84)–(10.94) | CORE | `steady_cd_exact`, `cd_layer_thickness`, `steady_cd_fd`, `steady_cd_discrete_exact`, `discrete_root`, `wiggle_indicator`, `stretched_grid`, `numerical_diffusivity` | V2 (10.86) solves (10.84), R → 0 gives x; D15–D17 · V1 scaled = printed (1e-13, \|R\| ≤ 50), finite at R = 1e4; tridiagonal = closed form (1e-11, 36 cases incl. R_cell = 2 and R < 0); e⁻¹, e⁻² levels · V7 wiggles iff R_cell > 2 · V3 | central 2.04 [2.13, 2.03, 2.007, 2.002]; upwind 0.89 [0.79, 0.88, 0.93, 0.96]; gap upwind ↔ exact of (10.94) 0.0134 → 1.1e-4 (order 1.73, 1 % of the error at n = 160) | analytic, symbolic, converged | the implementer's "gap shrinks, not zero": confirmed |
| C10 MacCormack (10.95)–(10.110) | CORE | `maccormack_step`, `maccormack_advection_1d`, `lax_wendroff_advection_1d`, `maccormack_linear_sympy`, `ns_fluxes`, `ns_coefficients`, `cavity_coefficients`, `coefficients_from_ns`, `ns_predictor`, `ns_corrector`, `weakly_compressible_step`, `maccormack_dt(_asymptotic)`, `wc_taylor_green`, `cavity_maccormack` | V1 = Lax–Wendroff every step (1e-11, both arrangements); generic flux form = scalar; corrector = six-substep periodic step (1e-14); ρ′ storage invariant (1e-13) · V2 engine (LW difference 0, local error −(uΔx²/6)(1 − C²)T_xxx, G² identity); D18 all 14 steps · V3 Gaussian; predictor consistency with (10.96)–(10.98) on a manufactured field (independent of the coefficient code); FF/BB − BB/FF · V7 C = 1.05 blows, 0.95 bounded · V4 periodic mass exact · V5 cavity (C14) | Gaussian 1.96 [1.89, 1.99, 2.00]; one-sided predictor 0.97, averaged 1.99; arrangement difference 1.98; TG 2.85 [2.82, 2.88]; periodic mass drift 0.0 | analytic, symbolic, converged, conserved | cavity mass drift (Open item 3), low-Mach error ∝ 1/Ma (Open item 4) |
| C11 splitting, Θ-scheme, projection (10.111)–(10.118), (10.129)–(10.133) | CORE | `split_linear_system`, `marchuk_yanenko`, `theta_scheme_linear`, `theta_scheme_amplification_sympy`, `splitting_order`, `projection_sympy`, `project`, `projection_stages` | V3 orders · V2 Θ Δt² coefficient vanishes only at 1 − 1/√2 (any λ₁, λ₂), z³ = (106 − 75√2)/6; projection engine; D19 incl. tiny example · V4 ∇·u after projection · V1 curl of the correction 0; boundary u* irrelevant | MY 0.99; Θ 2.00; Θ(¼) 0.995 on fine Δt (0.86 on the default list); div after/before ≤ 1e-12 | converged, symbolic, conserved | |
| C12 staggered grid (10.119)–(10.128) | CORE | `MacGrid`, `divergence`, `gradient`, `pressure_poisson_matrix`, `solve_pressure`, `pin_pressure`, `correct`, `checkerboard`, `collocated_gradient`, `gradient_null_space`, `dt_limit`, `ftcs2d_max_amplification`, `channel_poiseuille`, `taylor_green` | V1 checkerboard (collocated 0 exactly, staggered ±2/Δx); null spaces 4 vs 1; Poisson matrix symmetric, rank N − 1 (8 × 8: 63), pinned full rank; Poiseuille exact (quadratic ghosts); von Neumann scan of (10.127)–(10.128) · V2 D20 (cell budget), D21 · V3 · V4 | TG 1.99; linear-ghost Poiseuille 2.00 [2.00, 2.00]; face Laplacian/convection 1.98; Poiseuille 1.7e-13 | analytic, converged, conserved | **F1** (grids with a 1-cell direction) |
| C13 mixed FE, LBB, cylinder (10.134)–(10.137), (10.156)–(10.198) | CORE | `p2_shape(_grad)`, `p1_shape`, `iso_map`, `jacobian`, `tri_quad_7pt`, `integrate_element`, meshes, `p2_node_counts`, `euler_check`, `poiseuille_test`, `kovasznay_test`, `infsup_constant/table`, `stokes_cavity`, `assemble_saddle/newton_system`, `apply_dirichlet`, `newton_solve`, `cylinder_steady`, `march_unsteady`, `cylinder_forces`, `weak_ns_*` | V1 δ_ab, Σφ = 1, gradients (1e-8); quadrature exact to degree 5 (1e-15), not 6; J = 2 × area; Poiseuille exact; counts, Euler 1 and 0; steady state a fixed point of `march_unsteady` · V3 Kovasznay; Newton quadratic; β_h mesh-independent · V2 Kovasznay solves NS; D22; engines · V7 P1–P1 spurious modes; symmetric mesh ⇒ C_L = C_M = 0 · R6 printed fails | Kovasznay u 3.00, p 2.23 (16 → 32); Poiseuille 1e-14 (replace), 2e-11 (penalty); β_h P2–P1 0.3666/0.3677/0.3662 (n = 2, 4, 8), P1–P1 7 spurious modes, β_nonspurious 0.10 → 0.072 (n = 4 → 8) | analytic, converged, symbolic | the implementer's "P2–P1 β ≈ 0.366 mesh-independent" confirmed |
| C14 lid-driven cavity (Figs. 10.6–10.8) | CORE | `cavity`, `cavity_centreline(_v)`, `streamfunction`, `primary_vortex_centre`, `ghia_centreline`, `ghia_vortex_centre`, `hou_centres`, `cavity_error_vs_ghia` | V5 live MAC 64² vs Ghia; centre vs Ghia and Hou; our 128² tables; MacCormack 128² (verifier run) · V3 live 16/24/32; MacCormack 32/64/128 · V4 ∇·u, ψ on walls · V7 reversed lid mirrors the flow | MAC 64² max dev 0.00257 (0.26 % of U, the implementer's 0.0026 confirmed live); 128²: 0.46 % (Re 100), 0.32 % (Re 400); MacCormack 4.14 %, 1.38 %, **0.44 %** at 32/64/128² (regenerated tables); deviation order 2.36 (MAC 16–32), 1.58 (MacCormack) | benchmark, converged, conserved | MacCormack at 64² is outside the 1 % tolerance, 128² inside |
| C15 grid convergence, Richardson, GCI (Fig. 10.13, §10.6) | CORE | `grid_convergence_index`, `richardson_three`, `richardson_extrapolate`, `block_channel`, `body_forces`, `block_wall_density` | V1 GCI exact on f₀ + Khᵖ; D23; block wall density exact on manufactured fields (4 faces); body forces (uniform p → 0; linear p → exact C_D, C_L) · V3 cavity ψ_min · V7 block C_D monotone, \|C_L\| → 0 | ψ_min p = 2.00 (32/64/128; was 1.94 before the review fix), GCI_fine 0.164 %, f_ext −0.103528; block C_D p = 0.43, GCI_fine 8.9 % | analytic, converged, qualitative (block) | **F3**; block not asymptotic (Open item 2) |

NOTE rows that are coded (each ≥ 1 level): N05 (`advected_gaussian` solves (10.1), mass conserved) · N06 (ghost) · N07 (orders) ·
N08 (BTCS) · N10 (norms) · N11–N15 · N16 (upwind) · N17 (BTCS at β = 100) · N18 (Lax demo) · N19 (R_cell = 2) · N21–N28 (Galerkin
sympy, hats) · N29 (θ vs method of lines) · N30–N33 · N34–N38 (incl. R1, R2) · N43–N49 (incl. R3) · N52 (artificial
compressibility: Poiseuille 1.2e-12, ∇·u 1.7e-10 after 17 315 pseudo-steps) · N53 (weakly compressible TG) · N54–N57 · N58
(`maccormack_dt`) · N60–N61 · N62–N63 (predictor) · N64–N72 · N73 (MAC stability scan) · N74 (Θ) · N75 · N76 (saddle structure
[[A, B], [Bᵀ, 0]]) · N79–N84 (cavity wall densities exact on quadratics) · N85 (R5) · N86–N89 (block) · N90 (ρ′ storage) · N91
(R12) · N93–N95 · N97 (α, β orders 0.99, 2.00) · N98 (Newton) · N99 (through Kovasznay/Poiseuille) · N102/N103 (assembled blocks) ·
N104, N105, N108, N109 · N110 (Strouhal) · N113 (rod series). C-depth NOTE rows (N01, N03, N04, N20, N24, N39, N40, N50, N59, N77,
N92, N96, N100, N101, N103, N107, N112) are named only; the coded ones among them (N101 → R6, N103, N107) are covered above.
RECAP rows are tested in their own chapters; new instances here: R04, R05, R08, R09, R11, R12, R13 (see the tests of the same
name). SKIP rows S01, S02 are not coded.

## Derivations (curation §4b) — all 23 D rows re-derived with sympy, independently of the chapter's own engines
| D | ★ | test | what is checked |
|---|---|---|---|
| D01 | ★ | `test_stencils_D01_V2_derivation` | steps 2–8 coefficient by coefficient; the check numbers −0.0421, −9.005e-4 |
| D02 | ★ | `test_ftcs_D02_V2_derivation` | (10.9) solved for T^{n+1} = (10.10); weights sum 1; α, β dimensionless |
| D03 | ★★ | `test_truncation_D03_V2_derivation_on_an_exponential_mode` | steps 3–7 on e^{ax+bt} (independent of the code's Taylor table); step 11; the 0.0406 check |
| D04 | ★★ | `test_von_neumann_D04_V2_derivation` | steps 3–4 collection, (10.24), G(0) = 1, 0.6 − 0.2i, G(−θ) = conj G |
| D05 | ★★ | `test_von_neumann_D05_V2_derivation` | Re G, Im G, half angle, (10.26) |
| D06 | ★★ | `test_noye_D06_V2_derivation` | steps 2–5 factorisation, both ends (steps 9–10), step 13 |
| D07 | ★★ | `test_btcs_D07_V2_derivation` | bracket = 1 + 4β sin²(θ/2) + 2iα sin θ; \|den\|² − (Re)² = 4α² sin²θ; 1/401 |
| D08 | ★★ | `test_upwind_D08_V2_derivation` | steps 2–9 incl. C = 1 shift, \|G(π)\|² = (2C − 1)², 1.2 at C = 1.1, printed stencil for C < 0 |
| D09 | ★★ | `test_weak_form_D09_V2_derivation` | integration by parts on cubic T, w(0) = 0; the tiny example 1 = 2 − 1 |
| D10 | ★★ | `test_weak_form_D10_V2_derivation` | (10.37) identity; w = x/L gives D[T_x(L) − q] |
| D11 | ★★ | `test_galerkin_D11_V2_derivation` | my element integrals of 3 hats = engine M, K, F = `assemble_1d` numbers; K symmetric iff u = 0 |
| D12 | ★★ | `test_galerkin_D12_V2_derivation` | 2h/3, h/6, ±½, 0, 2/h, −1/h; row sums |
| D13 | ★ | `test_element_D13_V2_derivation` | (10.65), (10.66), dξ/dx = 2/h, slopes ∓1/h |
| D14 | ★★ | `test_element_D14_V2_derivation` | m^e, k^e on the parent element; the scatter-added row (step 10–11); numeric parity |
| D15 | ★ | `test_steady_cd_D15_V2_derivation` | `dsolve` with the BCs = (10.86); roots {0, u/D}; (10.88) limit; e⁻¹, e⁻² |
| D16 | ★★ | `test_steady_cd_D16_V2_derivation` | quadratic, roots {1, (1 + P/2)/(1 − P/2)}, closed form satisfies (10.91) for every j; tiny example; upwind root 1 + P |
| D17 | ★★ | `test_upwind_steady_D17_V2_derivation` | balance per Δx² = (P/Δx)a − (1 + P/2)a² at fixed P; D_num = uΔx/2 |
| D18 | ★★★ | `test_maccormack_D18_V2_derivation_step_by_step` + `test_maccormack_V2_chapter_engine_lax_wendroff_error_and_G` | every line: steps 1–5 (step 3 and 4 separately), 6–8 exact series, 9–12 local error C(C − 1)(C + 1)Δx³T_xxx/6, 13 G from the stencil, 14 both forms of \|G\|², C = 1 shift |
| D19 | ★★ | `test_projection_D19_V2_derivation` | steps 1–3, 8; tiny example u* = (x + y, 0) → p = x²/2, u = (y, 0), vorticity −1 before and after |
| D20 | ★★ | `test_pressure_poisson_D20_V2_derivation` | cell budget → (10.124) with symbolic Δx, Δy, Δt; wall cell has no p outside; 3-cell pipe numbers (0, 1, 2) |
| D21 | ★★ | `test_checkerboard_D21_V2_derivation` | (10.125) = 0 and (10.126) = ∓2/Δx for symbolic i, j; zigzag velocity |
| D22 | ★★ | `test_weak_ns_D22_V2_derivation` + `test_weak_ns_V2_chapter_engines_D22_N95` | steps 2, 4 (divergence of ∇u + ∇uᵀ), 8 (symmetric : antisymmetric = 0) |
| D23 | ★★ | `test_richardson_D23_V2_derivation` | ratio rᵖ, Richardson f₀ exact |

Every intermediate line of Part F that was checked follows from the previous one — **no wrong intermediate line found**
(including the D18 local-error sign, D16's discriminant P², D17's R_eff = 16/3 example, D20's rank 63, D23's 1 + h² example).

## Functions used by the notebook and explainers (design Part C) — each is called by at least one test
| Part C rows | tests |
|---|---|
| C.0 reused (`observed_order`, `pairwise_orders`, `richardson`, `grid_convergence_index`, `core.diffusion.ftcs_diffusion_1d`) | `test_gci_*`, `test_ftcs_V1_u_zero_reproduces_the_ch01_diffusion_solver`, every order test |
| C.1 1.1–1.8 | `test_stencils_*`, `test_stencil_taylor_coefficients_*`, `test_fd_derivative_*`, `test_mixed_derivative_*`, `test_advected_gaussian_*`, `test_part_c_V1_direct_calls_*` |
| C.1 1.9–1.18 | `test_ftcs_*`, `test_btcs_*`, `test_neumann_ghost_*`, `test_solve_transport_*`, `test_truncation_*`, `test_error_norms_*`, `test_von_neumann_*`, `test_fourier_mode_*` |
| C.1 1.19–1.26 | `test_von_neumann_*`, `test_noye_*`, `test_stability_verdict_*`, `test_phase_error_*`, `test_upwind_*`, `test_mac_stability_*` |
| C.1 1.27–1.32 | `test_steady_cd_*`, `test_layer_thickness_*`, `test_cell_peclet_*`, `test_stretched_grid_*`, `test_forward_scheme_R3_*`, `test_upwind_*`, `test_advect_periodic_*`, `test_time_differences_*`, `test_lax_equivalence_*` |
| C.2 (FEM1, 14 rows) | `test_hat_basis_*`, `test_parent_map_*`, `test_shape_slopes_*`, `test_element_*`, `test_connectivity_*`, `test_bilinear_form_*`, `test_interior_stencil_*`, `test_galerkin_*`, `test_fe_steady_*`, `test_fe_transient_*`, `test_solve_transport_V1_theta_*`, `test_weak_form_V1_*` |
| C.3 (MAC, 11 rows) | `test_projection_*`, `test_checkerboard_*`, `test_pressure_poisson_*`, `test_solve_pressure_*`, `test_mac_*`, `test_step_and_run_*`, `test_cavity_mac_*`, `test_part_c_V1_direct_calls_*` |
| C.4 (MCK, 8 rows) | `test_maccormack_*`, `test_ns_*`, `test_arrangements_*`, `test_weakly_compressible_*`, `test_cavity_density_bc_*`, `test_cavity_maccormack_*`, `test_block_*`, `test_body_forces_*` |
| C.5 (FEM2, 10 rows) | `test_p2_shapes_*`, `test_quadrature_*`, `test_iso_map_*`, `test_meshes_*`, `test_mixed_fe_*`, `test_kovasznay_*`, `test_lbb_*`, `test_stokes_cavity_*`, `test_newton_*`, `test_fe_time_*`, `test_fe_unsteady_*`, `test_fe_forces_*`, `test_apply_dirichlet_*` |
| C.6 (chapter module, 15 rows) | `test_weak_form_V2_*`, `test_compressible_ns_*`, `test_maccormack_V2_*`, `test_projection_V2_*`, `test_weak_ns_*`, `test_splitting_*`, `test_theta_scheme_*`, `test_checkerboard_*`, `test_artificial_compressibility_*`, `test_ghia_reference_*`, `test_cavity_*`, `test_fe_forces_*`, `test_rod_heating_*`, `test_book_slips_*`, `test_cfl_time_step_*`, `test_maccormack_D18_supports_*` (`derive_all`), `test_gci_*` (`richardson_three`) |
| C.7 `tools/convergence.grid_convergence_index` | `test_gci_V1_*`, `test_grid_convergence_*` |
| C.8 scripts and `ch10_drawings` helpers | `test_scripts_V7_every_ch10_script_runs_headless` (slow), `test_scripts_drawings_V1_*` (10 helpers draw) |
| every contract name (153) | `test_part_c_smoke_every_contract_name_is_exported_and_callable` (exported from `ch10`, callable, docstring present) |

## Convergence studies (observed order, `tools/convergence.observed_order`; pairwise orders in brackets)
| scheme / quantity | grids | observed | design |
|---|---|---|---|
| forward / backward / centred / second-derivative / one-sided stencils on sin | h = 0.1 … 0.00625 | 1.004 / 0.996 / 2.000 / 2.000 / 1.977 | 1 / 1 / 2 / 2 / 2 |
| time differences (10.8) forward, backward, leapfrog, trapezoidal on y′ = −y | Δt = 1/20 … 1/160 | 1.01, 0.99, 2.06, 2.00 | 1, 1, 2, 2 |
| FTCS, advected Gaussian, Δt = 0.25Δx²/D | n = 20 … 160 | 1.99 [1.85, 2.09, 2.00] | 2 |
| Crank–Nicolson, same | n = 20 … 160 | 2.03 | 2 |
| BTCS, same | n = 20 … 160 / 80 … 640 | 1.80 [1.53, 1.88, 1.95] / **1.98** [1.95, 1.99, 2.00] | 2 (pre-asymptotic on the default grids) |
| upwind, Δt = 0.5Δx/u | n = 20 … 160 / 160 … 1280 | 0.76 [0.75, 0.71, 0.83] / **0.947** [0.91, 0.95, 0.976] | 1 (pre-asymptotic on the default grids) |
| Lax–Wendroff (advective rule, D = 0.01) | n = 20 … 160 | 2.00 | 2 |
| MacCormack, Gaussian one revolution, C = 0.5 | N = 100 … 800 | 1.96 [1.89, 1.99, 2.00] | 2 |
| FTCS truncation remainder (measured − (10.17)) | Δx = 0.01/s, Δt = 0.001/s² | 4.00 [3.99, 4.00, 4.00] | 4 |
| FTCS → rod series (10.199) | n = 10 … 80 | 1.98 | 2 |
| steady centred / upwind (R = 10) | n = 10 … 160 | 2.04 [2.13, 2.03, 2.007, 2.002] / 0.89 [0.79, 0.88, 0.93, 0.96] | 2 / 1 |
| upwind vs exact solution of (10.94) | n = 10 … 160 | 1.73 [1.46, 1.69, 1.83, 1.91] (gap 0.0134 → 1.1e-4, not zero) | → 0 faster than the error |
| FE steady with Dirichlet + natural Neumann | n = 8 … 64 | 2.00 [2.005, 2.001, 2.000] (slope at L: 0.95) | 2 |
| FE transient (θ = ½) vs rod series | n = 8, 16, 32 | 2.00 | 2 |
| NS predictor consistency (one-sided / FF+BB average) | N = 32, 64, 128 | 0.97 / 1.99 | 1 / 2 |
| FF/BB − BB/FF after one step | N = 32, 64, 128 | 1.98 | 2 |
| weakly compressible TG, Ma = 0.05 | n = 16, 32, 64 | 2.85 [2.82, 2.88] | ≥ 2 |
| Marchuk–Yanenko; Θ (θ = 1 − 1/√2); Θ (θ = ¼) | Δt = 0.1 … 0.0125 (¼: 0.004 … 0.0005) | 0.99; 2.00; 0.995 (0.86 on the default list) | 1; 2; 1 |
| MAC Taylor–Green (Δt ∝ Δx²) | n = 16, 32 | 1.99 | 2 |
| MAC Poiseuille, linear wall ghosts | n_y = 8, 16, 32 | 2.00 [2.00, 2.00] (quadratic ghosts: exact, 1.7e-13) | 2 |
| MAC face Laplacian / advective = conservative convection | n = 16, 32, 64 | 1.98 / 1.98 | 2 |
| MAC cavity max deviation from Ghia, Re = 100 | 16, 24, 32 (live) | 2.36 [2.65, 1.99] | 2 |
| MAC cavity ψ_min (our tables) | 32, 64, 128 | 2.00 (GCI_fine 0.164 %, f_ext −0.103528; after the review fix of `primary_vortex_centre`, was 1.94) | 2 |
| MacCormack cavity max deviation from Ghia | 32, 64, 128 | 1.58 [1.56, 1.58] | 2 (wall closures, see Open items) |
| Kovasznay P2–P1 (L² u, p) | n = 8 → 16 → 32 | u 3.05 → 3.00; p 2.74 → 2.23 | 3; 2 |
| Newton residuals (Kovasznay n = 16) | iterations | 0.27, 0.069, 1.6e-4, 8.7e-9, 9.7e-17 — quadratic | quadratic |
| block C_D, Re = 20 (mean over t > 20) | Δx = 1/8, 1/16, 1/32 | **0.43** (4.272, 4.423, 4.534; GCI_fine 8.9 %) | 2 — not asymptotic (qualitative) |

## Conservation / invariant residuals
| quantity | residual |
|---|---|
| FTCS periodic node sum, 500 steps | < 1e-12 relative |
| advected Gaussian ∫T dx (infinite line), t = 0, 0.5, 2 | 1e-9 relative (trapezoid on ±10) |
| MAC projection ∇·u (20 × 14 walled, random u*) | ≤ 1e-12 × ∇·u* ; periodic FFT 1e-11 |
| discrete curl of the correction −Δt∇p | 0 to round-off (≤ 1e-10 relative) |
| MAC cavity 64² final ∇·u; ψ on the four walls | < 1e-11; < 1e-12 |
| MAC run 30 steps (10²) max ∇·u per step | < 1e-12 |
| compatibility Σ rhs ΔxΔy | ≤ 1e-9 (5.6e-17 for the E6 field) |
| weakly compressible MacCormack, periodic TG mass | 0.0 (exactly conserved) |
| MacCormack cavity mass (walls, continuity closures) | **−0.68 % (32²), −0.56 % (64²), −0.49 % (128²)** (tables regenerated with the additive Δt) at Ma = 0.08; −0.22 % (Ma 0.04), −1.9 % (Ma 0.16) at 64² — not conserved (Open item 3) |
| artificial compressibility channel ∇·u | 1.7e-10 at steady state |
| FE: symmetric mesh ⇒ C_L, C_M of the steady cylinder | < 1e-10 |
| block Re = 20: \|C_L\| (symmetric flow) | 0.082, 0.037, 0.018 for Δx = 1/8, 1/16, 1/32 → 0 at order ≈ 1.1 |
| mirror symmetry of the cavity under U → −U (MAC, 16², t = 2) | 1.7e-16 |

## Benchmarks used (public sources; all re-read 2026-09-30, see `reference/ch10/SOURCES.md`)
| benchmark | our value | source | verified |
|---|---|---|---|
| Ghia, Ghia & Shin (1982) Table I, u on x = ½, Re = 100 | live MAC 64²: max deviation 0.00257 = 0.26 % of U; cached 128²: 0.46 %; MacCormack 128² (verifier run): 0.48 % | J. Comput. Phys. 48, 387 (1982), doi:10.1016/0021-9991(82)90058-4; transcription: ivan-pi gist | three values re-read today in CMC 36, 1–21 (2013) Table 1 (u_min at Re 100, 400, 1000: −0.21090, −0.32726, −0.38289 — agree); all 17 Re = 100 values in the independent transcription cfdgasman/lid-driven-cavity `cavity/ghia.py` (agree); `make_refs.py` now asserts 5 of them |
| Ghia (1982) Table I, Re = 400 | cached MAC 128²: 0.32 %; MacCormack 64²: 5.86 % (Open item) | same | same |
| Ghia (1982) primary vortex centres | MAC 64² (0.6158, 0.7371) vs (0.6172, 0.7344): Δ = 0.0014, 0.0027 (< one cell) | via Hajabdollahi & Premnath, arXiv:1202.6351 | PDF re-read today |
| Hou et al. (1995) primary centres (lattice Boltzmann 256²) | Re 100: Δ = (0.0038, 0.0002) at 64²; Re 400 (our 128² table): (0.5553, 0.6061) vs (0.5608, 0.6078): Δ = 0.0055, 0.0017 | J. Comput. Phys. 118, 329 (1995), arXiv:comp-gas/9401003 Fig. 1 | PDF re-read today |
| Kovasznay (1948) exact steady NS flow | sympy residual of the coded field < 1e-12 (Re 10, 40); P2–P1 orders 3.00 / 2.23 | Kovasznay, Proc. Camb. Phil. Soc. 44, 58 (1948) (formula as coded; the test proves it solves NS) | sympy |

## Numbers from the text (book vs ours) — private values redacted to relative differences
| item | result |
|---|---|
| FTCS truncation coefficients (10.17) | identical (0 %) |
| stability limits: β ≤ ½, C ≤ 1, R_cell ≤ 2, FTCS orders (1, 2) | reproduced exactly as the boundaries of the computed \|G\| ≤ 1 / r > 0 regions; measured orders 1.01 and 1.99 |
| consistent-mass weights (10.63) | identical |
| layer levels at one and two layer thicknesses (percent, as printed) | ours e⁻¹ = 36.79 %, e⁻² = 13.53 %: 0.57 % and 0.25 % from the printed values — the book rounds to a whole percent and to 0.1 % (inside half a unit) |
| upwind numerical diffusivity factor (10.94) | identical |
| Θ-scheme θ (printed to five digits) | 1 − 1/√2 differs by 3.2e-6 absolute (1.1e-5 relative) — rounding |
| c₅ denominator (10.109) | identical |
| primary-eddy centres of Hou et al. quoted by the book | identical to the public data |
| the book's MacCormack centres (its much finer grid; Re 100, 400) vs our MAC 128² | inside the book's stated uncertainty (\|Δ\| ≤ 0.005 at Re 100, ≤ 0.015 at Re 400) |
| FE mesh counts (both meshes) | V − E + T = 0 (one hole) — consistent |
| Strouhal number from the book's period | consistent with the printed rounded value to 0.25 % (values private: tests/book_values_ch10.json) |
| our confined cylinder (W = 5d, ours) vs the book's Strouhal number | 2.2 % lower — different confinement (the book's geometry is private); qualitative (R9) |
| block C_D at Re = 20 | not compared: our channel (H = 4, 8 ahead, 20 behind) differs from the book's |
| 7-point quadrature | point count and degree identical |

## Figures reproduced with our code (`outputs/ch10/verify/`, local, git-ignored; `tests/ch10_verify_figures.py`)
| figure | png | visual verdict |
|---|---|---|
| (10.6)–(10.7) stencil errors | `c01_stencil_orders.png` | slopes 1 (forward, backward) and 2 (centred, second derivative, one-sided) over h = 1e-4…0.3, and the round-off floor rising as ε/h (first derivatives) and ε/h² (second derivative) below h ≈ 1e-5 / 1e-3 — as C01 teaches. |
| (10.24)–(10.27) von Neumann | `c04_von_neumann.png` | β = 0.51 exceeds 1 only near θ = π (1.04), (0.4, 0.2) exceeds 1 at mid θ, pure convection peaks at θ = π/2 (1.166); the brute-force stable set is exactly the lens between 4α² = 2β and 2β = 1. |
| C05/C10 advection at C = 0.8 | `c05_c10_upwind_maccormack.png` | upwind smears (Gaussian amplitude 0.745, square rounded), MacCormack keeps amplitude (0.978) with dispersive ripples trailing the square (upstream side), overshoot 1.17 — diffusion vs dispersion as E3 shows. |
| Fig. 10.x layer (C09) | `c09_cell_peclet.png` | at R_cell = 4 the centred nodes alternate (−0.33 at x = 0.9), upwind is monotone and sits on the exact solution of (10.94); the centred root passes through ±∞ at R_cell = 2 and is negative beyond, upwind r = 1 + R_cell stays positive. |
| FE = FD (N33) | `c07_fe_equals_fd.png` | FE circles and FD crosses coincide (≤ 8e-17) for n = 4, 8, 16 and approach the exact layer. |
| projection (C11/C12) | `c12_projection.png` | ∇·u* up to 9.3, the pressure saddle from (10.124), ∇·u after ≤ 9.5e-13 everywhere. |
| LBB (C13, F7) | `c13_infsup.png` | P2–P1 flat at 0.366, P1–P1 smallest nonzero constant falling like h (plus 7 exact zeros). |
| Figs. 10.7–10.8 cavity | `c14_cavity.png` | centreline u converges onto Ghia's points at Re 100 (MAC 16² → 128², MacCormack 128² on top); at Re 400 MAC 128² on the points while MacCormack 64² misses the minimum (−0.27 vs −0.327); 64² streamlines with the primary vortex at Ghia's/Hou's markers and the two bottom corner eddies (the right one larger). |
| Fig. 10.13 grid convergence | `c15_grid_convergence.png` | block C_D histories keep an acoustic oscillation (±0.05) to t = 30 and their means rise with refinement at p ≈ 0.43 (not asymptotic, Richardson 4.86 far off); cavity ψ_min converges at p = 1.94 to −0.10352. |
| Fig. 10.21 FE cylinder forces | `c13_cylinder_re100_forces.png` | after the rotation kick (t < 2) C_L settles into a periodic oscillation (amplitude ≈ 0.4, St = 0.2054), C_D ≈ 1.8 steady-periodic — a Hopf-type periodic wake in the confined channel. |

## Deviations & justifications (every `# DEVIATION` in the code)
| where | deviation | verifier check |
|---|---|---|
| `fd.stretched_grid` | the design's map clusters nodes at x = 0; the code uses tanh(β_s j/n)/tanh β_s, clustering at the layer (x = L) | correct: last spacing < first (test); the design formula indeed has dx/dj largest at j = n. The stretched grid shrinks but does not remove the wiggles (doc correction) |
| `mac._pad_u` | tangential no-slip by ghost values (linear 2U − u₀ or quadratic) — the book leaves it open | linear ghosts give an exact O(Δy²) error (Poiseuille 2.00), quadratic ones reproduce Poiseuille to 1.7e-13 |
| `maccormack.maccormack_dt` | zero speeds left out of Re_Δ = min(…) | the formula → (10.155) × 1/(1 + √2 Ma) as Re_Δ → ∞ (test to 1e-9) |
| `maccormack.make_cavity_bc` | top-corner density extrapolated instead of (10.139)/(10.141) at the corners | runs stable; cavity converges to Ghia (0.48 % at 128²); mass not conserved (−0.5 %, Open item 3) |
| `maccormack.block_channel` | inflow density by zero-gradient extrapolation instead of continuity | runs stable to t = 30 at three grids; C_D not asymptotic (Open item 2) |

## Open items (verbatim for the orchestrator)
1. ~~F1–F3~~ **closed in loop 1** (327/327 pass; 33/33 mutants killed). Remaining for the **notebook-builder** (analysis/design
   wording, not code): design C.1 1.14 default-grid expectations (upwind 1.0, BTCS 2.0 hold only on n = 160…1280 / 80…640);
   design C.6 6.6 (Θ at θ = ¼ gives 0.86 on the default Δt list, 0.995 on fine Δt); cavity prose (MAC 64² 0.26 %, MacCormack
   64² 1.38 %, 128² 0.44 %; MAC 128² vs 64² at Ghia's own error level — cite ψ_min p = 2.00); D20 note (the 3-cell pipe is now
   exactly `MacGrid(3, 1)`); N47 "removes the wiggles" → "shrinks them below 1e-4"; `kovasznay_test` orders — quote 3.00/2.23
   (16 → 32), not the 8 → 16 pair.
2. **Block C_D is not grid-converged (qualitative):** Re = 20, mean over t > 20: 4.272, 4.423, 4.534 for Δx = 1/8, 1/16, 1/32 —
   observed p = 0.43 (the implementer's 0.39 with another window), GCI_fine 8.9 %, Richardson 4.86 not trustworthy; the lift does
   go to zero (0.082 → 0.018). Hypothesis: the heuristic block-wall density closures (10.151)–(10.154) with averaged corners (and
   the corner singularities) are first-order-or-worse, plus a persistent acoustic oscillation (±0.05) of the weakly compressible
   run. The notebook's C15 must call this a *non-asymptotic* study (it cannot say "grid independent"); a Δx = 1/64 run and a longer
   average would be the next evidence.
3. **MacCormack cavity mass is not conserved:** relative drift −0.68 %, −0.56 %, −0.49 % (32², 64², 128², Ma 0.08); −0.22 % and
   −1.9 % at Ma 0.04 and 0.16 (64²). The wall-density updates (10.139)–(10.146) are not in conservation form. Report it with the
   runs; do not label the cavity MacCormack solver "conserved".
4. **Low-Mach accuracy of MacCormack (N53 wording):** on a fixed grid the weakly compressible error *grows* as Ma falls (TG 32²:
   0.0062, 0.0099, 0.0197 for Ma = 0.1, 0.05, 0.025; cavity 64²: 1.64 %, 1.43 %, 1.41 % for Ma = 0.04, 0.08, 0.16) — Δt ∝ Ma·Δx
   and the scheme's acoustic truncation error ∝ cΔx² accumulate. The O(Ma²) compressibility error is real but hidden behind this;
   the notebook should say so (it is why the book's cavity uses a much finer grid).
5. **MacCormack cavity at 64² is 1.38 % of U from Ghia (outside a 1 % tolerance);** within tolerance at 128² (0.441 %, verifier
   run cached in `outputs/ch10`, table `reference/ch10/cavity_mck_re100_n128.csv`). At Re = 400, 64²: 5.86 % — no 128²
   MacCormack Re = 400 run exists (≈ 7 min; optional).
6. **Cached-data evidence:** the MAC 128² tables, the block tables and the unsteady cylinder history are checked as data (V5/V3/V7
   on the committed CSVs); the default suite recomputes MAC 16²–64², MacCormack 16² (and 32² in `slow`), the Δx = 1/8 block run
   (`slow`, reproduces its table to 6 s.f.) and the steady cylinder, but not the heaviest runs. Mutants in code only exercised by the
   heavy runs are therefore invisible to the default suite.
7. **Confined-geometry numbers are qualitative:** cylinder St = 0.2054 (W = 5d, ours), block C_D — never compare them with
   unbounded or the book's data as if equal (R9).

## Documentation corrections required (analysis / design / docstrings; the verifier does not edit them)
Loop 1: the code-side items (R3 wording, stretched-grid docstring, `lax_demo`, MacCormack conservation/low-Mach, block label)
are fixed in `fluidpy/`; the design/analysis items stay open for the notebook-builder (Open item 1). `kovasznay_test`'s
docstring still says "≈ 3 and ≈ 2" — acceptable, but the notebook should quote 3.00/2.23.
- **R3 wording (code text shown to learners):** `fd.py` module docstring (line 17), `steady_cd_fd` docstring and
  `ch10.book_slips()` R3 evaluator say the downwind ("forward") difference "wiggles for every R_cell". It does not: its root is
  r = 1/(1 − R_cell) > 0 for R_cell < 1 (monotone, anti-diffusive D(1 − 0.5R_cell)); it wiggles only for R_cell > 1
  (`test_forward_scheme_R3_V1_downwind_root_wiggles_above_Rcell_one`).
- **N47 / stretched grid:** "removes the wiggles" → "shrinks them from 0.35 to below 1e-4 (visibly removes them)"; the sign
  alternation remains in the coarse upstream cells (local R_cell ≈ 8) — `wiggle_indicator` still reports it.
- **Design C.1 1.14 expectations:** "upwind advective 1.0 ± 0.15" and "BTCS diffusive 2.0 ± 0.15" fail on the *default* n_list
  (0.76, 1.80 — pre-asymptotic, Gaussian width 2.5–20 cells); they hold on n = 160…1280 (0.947) and 80…640 (1.98). State the grids
  or change the defaults.
- **Design C.6 6.6:** "Θ-scheme 1.0 at θ = 0.25" on the default dt_list gives 0.86; Δt = 0.004…0.0005 gives 0.995.
- **`lax_demo` docstring / design C.1 1.32:** "≈ 1.03 per step", "≈ 900–1000 steps" → \|G(π)\| = \|1 − 4β\| = 1.04 and the blow-up
  (growth cap 1e6) happens at step 690.
- **Cavity comparisons in prose:** say "MAC 64²: 0.26 % of U; MacCormack 64²: 1.38 %, 128²: 0.44 %"; note that MAC 128² (0.46 %)
  is not closer to Ghia than 64² (0.26 %) because both are at the level of Ghia's own 129² discretisation error — the
  self-convergence (ψ_min, p = 2.00) is the evidence of convergence.
- **`kovasznay_test` docstring:** "orders ≈ 3 and ≈ 2" — measured 3.05/2.74 (8 → 16) and 3.00/2.23 (16 → 32): quote the finer pair.
- **Design D20 check:** the 3-cell pipe rows hold for `MacGrid(3, 1)` only after F1 is fixed.

## Discrimination (planted wrong variants on a scratch copy of `fluidpy/`)
Harness: copy `fluidpy/`, `tools/`, `scripts/`, `reference/ch10/`, `book.yaml` and the ch10 tests to the scratchpad (no
`outputs/` caches, so every cached function recomputes), plant one text change, run `pytest tests/test_ch10.py -x -m "not slow"`
minus the 11 known F1–F3 items. Baseline on the unmutated copy: 312 passed. **30 mutants planted, 30 killed, none hung** (the
first failing test is named):

| # | planted variant | killed by |
|---|---|---|
| M01 | α = uΔt/Δx (the 2 of the centred difference lost, (10.11)) | `test_ftcs_V1_design_five_point_example_and_coefficients` |
| M02 | FTCS convective sign flipped (10.10) | same |
| M03 | G of (10.24) with 1 − β | `test_von_neumann_V1_G_modulus_equals_10_26_and_design_values` |
| M04 | Noye 2β ≤ 1 → β ≤ 1 (10.27) | `test_solve_transport_V7_stability_guard_and_growth_cap` |
| M05 | upwind side not taken from sign(u) (slip R11) | `test_von_neumann_V1_a_single_fourier_mode_…[upwind--0.3-0.1]` |
| M06 | stencil divided by Δx instead of Δx^m | `test_fd_derivative_V1_exact_on_polynomials_and_V3_periodic_order` |
| M07 | exact layer (10.86) at the wrong wall | `test_steady_cd_exact_V2_solves_10_84_and_V1_overflow_safe_form` |
| M08 | discrete root without the ½ | `test_steady_cd_fd_V1_tridiagonal_solve_equals_the_discrete_closed_form` |
| M09 | numerical diffusivity 1.0 R_cell D instead of 0.5 (10.94) | `test_upwind_steady_V1_nodes_approach_the_modified_equation_solution` |
| M10 | Richardson without the −1 | `test_gci_V1_exact_on_power_law_models` |
| M11 | element mass (h/6)[[4, 1], [1, 4]] | `test_interior_stencil_V1_assembled_rows_are_10_63_times_h` |
| M12 | convective element block transposed | `test_galerkin_V1_steady_FE_equals_centred_FD_N33` |
| M13 | Dirichlet column enters F with + sign (10.57) | `test_fe_steady_V3_dirichlet_and_natural_neumann_converge` |
| M14 | natural Neumann term D q dropped (10.77) | same |
| M15 | printed slope labels by default (slip R1) | `test_shape_slopes_V1_chain_rule_and_R1_printed_labels_fail` |
| M16 | discrete divergence with −∂v/∂y (10.123) | `test_projection_V4_divergence_removed_to_round_off` |
| M17 | correction u + Δt∇p (10.121) | same |
| M18 | wall rows of (10.124) Dirichlet-like (diag −3) | same |
| M19 | viscous term × Re in the predictor (10.119) | `test_mac_V1_channel_poiseuille_exact_with_quadratic_ghosts_and_order_two_linear` |
| M20 | convective term u u_x − v u_y | `test_mac_V3_taylor_green_decays_at_the_exact_rate` |
| M21 | MAC limit (10.127) without the factor 2 | `test_mac_stability_V1_von_neumann_scan_of_10_127_and_10_128_N73` |
| M22 | MacCormack corrector forward instead of backward | `test_maccormack_V1_equals_lax_wendroff_every_step` |
| M23 | cross viscous term a₉ with the wrong sign | `test_ns_predictor_V3_consistent_with_the_compressible_equations_10_96_to_10_98` |
| M24 | a₅ without the 4/3 | `test_ns_coefficients_V1_10_109_and_cavity_coefficients_are_the_same_numbers` |
| M25 | pressure coefficient a₃ with Ma instead of Ma² | same |
| M26 | block front face 8/9 → 4/9 (10.151) | `test_block_wall_density_V1_exact_on_manufactured_fields_10_147_to_10_154` |
| M27 | P2 shape φ₆ = 4ηξ (10.185) | `test_p2_shapes_V1_kronecker_partition_and_gradients_N105` |
| M28 | 7-point rule centroid weight ¼ (10.198) | `test_quadrature_V1_seven_point_rule_exact_to_degree_five_not_six_N108` |
| M29 | pressure sign in the FE momentum residual (10.195) | `test_mixed_fe_V1_poiseuille_exact_R6_printed_fails_N109` |
| M30 | checkerboard (−1)^i only | `test_checkerboard_V1_collocated_gradient_zero_staggered_two_over_dx_D21` |

The printed slips kept as named options are, in addition, each asserted to FAIL their check: R1 `shape_slopes(printed=True)`,
R2 `connectivity(printed=True)`, R3 `steady_cd_fd(scheme="forward")` (wiggles for R_cell > 1), R5
`cavity_maccormack(printed_step5=True)` (NaN mass within 1 time unit on 16²), R6 `poiseuille_test(printed=True)` (error 1.7e55),
R11 `"upwind_printed"` with u < 0 (\|G\|max 1.5, explodes), R12 asymptotic Δt 1.71× the full limit at Ma = 0.5. Limitation
(Open item 6): the mutants did not touch code reachable only through the cached heavy runs.
