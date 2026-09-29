# Chapter 8 verification — Laminar Flow                     2026-09-28, commit 641a3dd + working tree (loop 2)

Verifier: math-verifier. Suite: `tests/test_ch08.py` (101 test functions = 101 collected items; 1 marked `slow`),
figures/metrics: `tests/ch08_verify_figures.py` → `outputs/ch08/verify/` (git-ignored), cited data: `reference/ch08/`
(`make_refs.py`, `benchmarks.json`, `SOURCES.md`). Nothing in `fluidpy/` or `scripts/` was edited by the verifier.

## Environment
Python 3.11.5 · numpy 2.4.6 · scipy 1.17.1 · sympy 1.14.0 · pint 0.25.3 · matplotlib 3.11.2 (Windows, `.venv`).

Runs (serial, no `-n`), loop 2: `pytest tests/test_ch08.py -q -p no:cacheprovider` → **101 passed in 169 s** (with the
private book file; without it the three V6 tests skip cleanly — checked in loop 1). Full suite
`pytest tests -q -p no:cacheprovider` → **983 passed in 660 s** (ch01–ch07 + machinery 882 + the 101 ch08 items; no
warnings reported). Loop 1 for reference: 97 passed in 101 s, full suite 979 in 606 s; the sympy engines dominate
(`stokes_sphere_sympy` 8 s, D28 6 s, `lubrication_nondim_sympy` 6 s, `slider_bearing_sympy` 5 s), and the loop-2
wall time also reflects a busier machine (other agents were running) — every test kept, none marked slow except the scripts.
All 12 `scripts/ch08_*.py` exit 0 headless with `--no-show` (2–5 s each; `ch08_common.py`, `ch08_drawings.py` are
helpers). `tests/ch08_verify_figures.py` runs and writes 10 PNGs.

Collected items per evidence tag (first tag of the `def` line): **V1 46 · V2 31 (22 of them `V2 derivation`) · V3 5 ·
V4 4 · V5 6 · V6 3 · V7 6** (= 101; loop 2 added three V1 and one V2, two of which also carry a second level inside). Most tests carry further levels inside (a V1 test that also measures an order, a V2
derivation that also checks the design's numbers); the table lists every level exercised per item.

## Loop 1 — what the first runs found (all test-side; no physics defect)
First full run: 86 of 96 passed. All ten failures were in the tests, fixed before this report, none by loosening a
derived tolerance: finite-difference steps chosen too small (round-off: the circular-Couette ODE with Δ = 1e-7R₂; the
`inertia_viscous_ratio` ∝ U check asserted 1e-8 where its own nested stencils carry ε·r/h_rel² ≈ 2e-6 — now 2e-5 with
the bound written in the test) or too tight for the stencil truncation they measure (pressure round-off of p₀ = 1e5 Pa;
the circular-Couette (4.39b) residual and the zero net viscous force now asserted by their O(h²) *order* plus a
(h/R)²-sized bound; `stokes_first_stopped` residual now an order-2 study in Δy); linear interpolation (h²/8·F″ ≈ 3e-6)
used to compare two `solve_bvp` grids (now cubic splines, 2e-7 = twice the per-solution 1e-7); a sympy `erf` vs `erfc`
simplification; a design-D32 reading (see O3: the θ-component of u·∇u has no 1/r² term — the radial one carries it); a
wrong-variant check normalised by the largest C_D (now point-wise); and the "Re_a r/a ~ 1 crossover" asserted with
prefactor 1 (loop 1 then quoted a radial-component-only prefactor (8cos²θ − 4sin²θ)/(16cos θ); that was wrong — the
magnitudes give ½ on the axis and at π/2, see O3 as corrected in loop 2).

## Loop 2 — after the derivation review (`reports/ch08_review.md`) and the implementer's fixes
Four tests added (97 → 101), one test comment and one test comment-block corrected, nothing loosened:

| review item | new test | what it proves | wrong variant planted (scratch copy) | caught |
|---|---|---|---|---|
| M1 extremum of (8.5) | `test_couette_poiseuille_state_V1_extrema_match_dense_grid` | for h = 0.01 m, μ = 1e-3 Pa s, (U, dp/dx) ∈ {(0.01, ±0.5), (0.01, ±1), (0, ±0.5)}: u_max, u_min agree with the max/min of `channel_flow` on 20 001 points (walls included) within the grid's own bound 0 ≤ u* − u_grid ≤ \|dp/dx\|/(2μ)(Δy/2)² (= 3.1e-11 m/s at 1 Pa/m); y_umax, y_umin within Δy/2 = 2.5e-7 m of the grid arg-extrema (ties at both walls allowed); interior extremum at y* = h/2 − μU/(h dp/dx) to 1e-12 (U = 1 cm/s, dp/dx = −0.5 Pa/m ⇒ 7 mm); backflow False for U = 0 at dp/dx = ±0.5, ±1 (new semantics); threshold 2μU/h² = 0.2 Pa/m unchanged: 0.19/0.21 → False/True, mirrored for U = −0.01 | `ys = h/2 + μU/(h g)` (the old sign) | yes — u_min at (0.01, 0.5) reported 0.0 at a wall against the grid's −2.25e-3 m/s |
| SF1 Oseen wake direction | `test_oseen_V2_streamfunction_solves_oseen_equation_and_wake_downstream` | (1) coded ψ = the (8.53) closed form on a 30 × 31 field for Re = 0.3, 1, 3 (max rel. 1e-11, so the sympy statements are about the coded function); (2) sympy, exact rational points (30 digits): E⁴ψ − (U/ν)∂ₓ(E²ψ) = 0 with ν = 2aU/Re, ∂ₓ = cos θ ∂_r − (sin θ/r)∂_θ, residual < 1e-20 of the term size at 4 points; the opposite sign leaves ≥ 1e-3 of the term size (e.g. −0.24 at r = 2.3a, θ = 0.7, Re = 0.8); (3) fluid-frame \|ψ\|(θ) > 1.05 \|ψ\|(π − θ) for 14 θ ∈ [0.1, 1.4] at r = 1.5a, 3a, 8a (Re = 1), and the ratio grows with r (0.36 rad: 1.49 → 1.92 → 3.77) | θ ↔ π − θ in `oseen_streamfunction` (s with 1 + cos θ) | yes — parity 2.1e-2 (and the older `test_oseen_V1_velocity…` also fails, since ψ and u then disagree) |
| O3 (my loop-1 finding, corrected) | `test_inertia_viscous_ratio_V1_half_Re_a_r_asymptote` | ratio/(Re_a r/a) at θ = 0, π/2, π for r/a = 100 … 800: 0.4981 (axis) and 0.4991 (π/2) at 400; deviation from ½ ∝ a/r (observed order 1.00 ± 0.05); Richardson 2f(2r) − f(r) = ½ within (a/r_min)² + 3ε r_max/(h_rel² a) = 1.5e-4 (observed ≤ 5.2e-5); V2 sympy on (8.49): the r → ∞ limit of \|u·∇u\|/\|ν∇²u\| ÷ (Ua/ν · r/a) is exactly ½ at θ = 0 and π/2 | ratio × ¼ (the loop-1 "≈ 1/8" prefactor) | yes |
| coordinator: slider recirculation at either end | `test_slider_backflow_both_ends` | brute force from `slider_gap_velocity(frame="pad")` on 401 × 799 (x, η) points, independent of `slider_bearing_state`: recirculation (pad-frame u against the flux sign) iff h > 1.5h_m, h_m = 2(1 + α)h₀/(2 + α) — re-derived here as u/U = −s + 3(1 − h_m/h)s(1 − s), s from the pad; α = 0.99/1.01 → False/True, −0.49/−0.51 → False/True, 0.5 False, 2.0 and −0.7 True, each for U = ±5 m/s; `backflow_any` = `inlet_backflow` = grid verdict; `backflow_x` = the station of the strongest reversal (L for α > 0, 0 for α < 0), NaN otherwise; the set of recirculating stations = {h > 1.5h_m} except within one x-cell of the edge; h₀ = 50 µm, α = −0.7, U = 5 m/s, μ = 0.05 Pa s: pad-frame u at x = 0 spans −4.98 … +0.2930 m/s on this grid (analytic peak 0.29304 m/s, grid agrees to 1e-5; the floor value itself is −U = −5 m/s, so the coordinator's "about −4.74" is the last interior point of a coarser y grid — the lower end depends on the grid, the upper end does not) | (a) `x_wide = L` always (the old x = L-only check); (b) threshold h > h_m instead of 1.5h_m | both yes |
| test comment at the temporal-BL check | — | the comment said "stokes_first(2U)", the assertion (correct) uses U: U − u(y) = U erfc(y/2√(νt)) = `stokes_first_problem(y, t, U)`; the jump across the sheet is 2U, the half-profile deficit is U | — | — |

The first run of the new ratio test failed at its a-priori Richardson bound 5e-5 (observed 5.18e-5 at θ = π/2, pair
400/800): that bound ignored the finite-difference round-off of `inertia_viscous_ratio` itself, ε r/(h_rel² a) ≈ 1.8e-5
per value at r = 800a, tripled by 2f(2r) − f(r). The bound is now written from that budget in the test (1.5e-4 total);
the physics claim (limit ½, O(a/r) approach, order 1.00) was never in question and is additionally pinned by sympy.

`slider_bearing_state["inlet_backflow"]` changed meaning with the implementer's generalisation: it is now True for
α < −½ as well (recirculation at the wide end, wherever it is), and equals `backflow_any`. The old assertion
`inlet_backflow is (α > 1)` for α ∈ {0.1, 0.8, 1.5} still holds; the explainer/notebook text that calls it "inlet"
should be read as "wide end" (O9).

One earlier probe (not a test) looked like a physics failure and was mine: `thin_film_spread`'s clock starts at 0, so a
run started on Huppert's profile at t₀ = 10 s must be compared with the similarity solution at t₀ + T, not at T.

## Validation table
| Concept / Eq. | Tier | fluidpy target | Evidence (V-levels) | Numbers (error, order, residual) | Label | Notes |
|---|---|---|---|---|---|---|
| C01 Re = Ud/ν, regimes; ν = μ/ρ; L²/ν (§8.1) | CORE | `pipe_flow_regime`, `momentum_diffusivity`, `diffusion_time`, `inertia_viscous_scales` | V1, V5, V7, V6 | Re exact; band edges exact at 2000/3000 (1e-12 either side); ν_air/ν_w(293 K) = 15.007; Sutherland from USSA constants 1e-12; scale/units invariance 1e-14 | analytic, benchmark | wrong: radius Re; μ instead of ν |
| R01, R04, R05 (inertia/viscous, (8.2), (8.3)) | RECAP | `inertia_viscous_scales`, `wall_bc_residuals` | V1, V7 | ratio = Re (1e-14); normal/tangential residual split exact | analytic | wrong: t·u ignored |
| C02 (8.5) Couette–Poiseuille | CORE | `channel_flow` (+ `_rate`, `_shear_stress`, `_backflow_threshold`, `couette_poiseuille_state`) | V1, V2, V7 | (4.39b) residual ≤ round-off bound (stencil exact for quadratics); ch04 parity 1e-14; sympy residual 0; D04 re-derived; Q by quad 1e-12; backflow threshold by bisection 2e-3 (y-grid) of 2μU/h² = 2 Pa/m; y_rev = 5 mm; loop 2: state extrema = dense-grid max/min within (dp/dx)/(2μ)(Δy/2)², locations within Δy/2, y* = h/2 − μU/(h dp/dx) (M1) | analytic, symbolic | wrong: G passed as dp/dx; printed V (units); old +μU/(h dp/dx) extremum |
| N04–N11 (v ≡ 0, (8.4a,b), constant dp/dx, Q, V, cases) | NOTE | `parallel_flow_sympy("channel")`, `channel_flow_rate` | V2, V1 | all engine residuals 0; non-constant dp/dx leaves −f′(x)/ρ ≠ 0 | symbolic | R9 printed V = ∫u dy has m²/s (pint) |
| C03 (8.6) pipe Poiseuille, (8.7), (8.8), Q, f = 64/Re | CORE | `pipe_poiseuille`, `pipe_shear_stress`, `pipe_wall_stress`, `pipe_flow_rate`, `pipe_friction_factor` | V1, V2, V4 | ch04 parity 1e-14; ch03 `pipe_profile` (z→∞) 1e-12; 3-D (4.39b) residual at round-off; D06, D07 re-derived; Q quad 1e-12; u_max = 2V; f·Re = 64 (1e-12); CV force balance πa²Δp + 2πaLτ₀ = 0 (1e-15) | analytic, symbolic, conserved | wrong: 1/(2μ) (Laplacian gives 2G/μ); A ln R unbounded |
| C04 (8.9)–(8.12) circular Couette | CORE | `circular_couette` (+ `_pressure`, `_shear_stress`, `_power`, `_state`) | V1, V2, V4, V7 | walls 1e-13; ω_z = 2A by FD; Euler ODE residual ≤ 1e-5 scale; (4.39b) residual → 0 at order 2.00; D08, D09 re-derived; Wikipedia A, B form identical; R₂ = 1e6R₁ vs ∞ branch 1e-10; ch05 `rotating_cylinder_flow(ω = 2Ω₁)` 1e-14 (ω = Ω₁ off by 2) | analytic, symbolic, conserved | Rayleigh criterion μ > η² reproduced |
| N20, R10, R11, R12, N21, N22 | NOTE / RECAP | `circular_couette_pressure`, `_shear_stress`, `_power`, `advective_acceleration_check` | V1, V4 | dp/dR = ρu²/R (FD 1e-7); torque 2πR²σ constant (1e-13); power in = dissipation (closed 1e-12, quad 1e-10); R₂ = ∞: 4πμΩ₁²R₁²; ch05 dissipation parity 1e-10; zero net viscous force → 0 at order 2 (vs O(1) for a Lamb–Oseen vortex) | analytic, conserved | R15: printed (2πR₁)σu < 0 = −power |
| C05 (8.14)–(8.17) lubrication scaling and balance | CORE | `lubrication_nondim_sympy`, `lubrication_scales`, `lubrication_term_magnitudes` | V2, V1 | coefficient sets {ε²Re_L, 1/Λ, ε², 1}, {ε⁴Re_L, 1/Λ, ε⁴, ε²}; D10 re-derived line by line with O(1) test fields; D11 limit + re-dimensionalisation; engine film ε = 1e-3, Re_L = 4350, ε²Re_L = 4.35e-3, Λ = 49.35; numeric magnitudes = sympy coefficients (1e-12) | symbolic, analytic | R6: printed ∂p/∂x in (8.13b) gives ε/Λ; R7: (8.17a) without ν fails pint |
| C06 (8.18)–(8.19), gap flux, 1-D Reynolds equation | CORE | `lubrication_velocity`, `lubrication_flux`, `reynolds_pressure_1d`, `reynolds_equation_sympy` | V1, V2, V3 | u(0) = U₀, u(h) = U_h; flux = quad 1e-12; μu″ = p_x; D12, D13 re-derived (Leibniz, kinematic cancellation, h/6 and h/2 pieces); callable-h route = exact slider 1e-10; array-h trapezoid order **2.000** | analytic, symbolic, converged | R8b: `form="book"` gives u(h) = U_h + U₀ |
| N36 Hele-Shaw | NOTE | `hele_shaw_velocity`, `_potential`, `_mean_velocity`, `hele_shaw_cylinder`, `hele_shaw_streamfunction_grid` | V1, V2, V3 | u = ∂φ/∂x 1e-12; no slip; ū = −(h²/12μ)∇p (quad 1e-12); gap average = ch04 ideal cylinder 1e-12; ∇²φ = 0 for harmonic p (sympy); grid solve order **1.11** (staircase disk) | analytic, symbolic, converged | wrong: φ sign; see O5 on the order |
| C07 Example 8.1 slider bearing, load, optimum | CORE | `slider_bearing`, `_load`, `_optimum_taper`, `_state`, `slider_gap_velocity`, `slider_bearing_sympy` | V2, V1, V5, V7 | D14 re-derived through all 15 steps (step 12 numerator identity, step 13 squared denominator, steps 14–15); engine residuals 0; quad ∫p dx = W (1e-9) for α = 1e-4…3; series/closed branches continuous (1e-8); San Andrés K_opt 2.18870 (−9e-5), W* 0.026707, P_max(K) identity 1e-12, K = 2.4142, 0.0429; pad-frame flux = C₁ at every x (1e-9); loop 2: recirculation iff h > 1.5h_m at either end (α > 1 or α < −½), brute-force pad-frame sign test on 401 × 799 points for α = ±(threshold ± 0.01), U = ±5 m/s | symbolic, analytic, benchmark | R8: printed first power fails the ODE and carries a > 5 % wrong load; (1 − αx/L) integrands miss p(L) = p_e |
| C08 Example 8.3 thin film | CORE | `thin_film_velocity`, `thin_film_flux`, `thin_film_spread`, `viscous_current_similarity`, `thin_film_state`, `viscous_current_eta_N` | V1, V2, V4, V5 | BCs u(0) = 0, u_y(h) = 0; flux = quad 1e-12; D15 re-derived; Huppert form solves the PDE exactly (sympy) and holds ∫h = A (1e-8); volume drift **2.2e-16**; front slope **0.1965** (target 0.2 ± 0.005); x_front/x_N = 1.019 at 1000 s; from the similarity profile L2 < 2 % at N = 200, 400, 800; η_N = 1.411245 (Huppert 1.411) | analytic, symbolic, conserved, benchmark (`thin_film_spread` itself: conserved, benchmark — no grid-refinement order is claimed; the N = 200/400/800 row below is a benchmark distance, not a convergence study) | wrong: h³/(2μ) (1.5× flux) |
| C09 (8.20)–(8.31) Stokes' first problem | CORE | `stokes_first_problem`, `similarity_variable`, `stokes_first_vorticity`, `vorticity_content`, `diffusion_thickness`, `stokes_first_state`, `stokes_first_stopped`, `similarity_ode_solve`, `crank_nicolson_1d` | V1, V2, V3 | collapse over three (U, ν, t) sets 1e-14; ch04 parity 1e-14 (and identity re-export); sympy residual 0; D17–D19 re-derived; `solve_bvp` 3.7e-13 (Stokes), 2.5e-9 (line vortex), η_max 10…14 agree 2e-7; CN time order **2.000**, space **1.995**; BE **0.996** (wrong variant); impulsive start with two BE steps **1.99**, pure CN stalls; stable at 50× FTCS; ∫ω dy = +U (1e-10); η₉₉ = 3.6428 | analytic, symbolic, converged (the "converged" is Crank–Nicolson's; `similarity_ode_solve` is analytic — closed-form agreement and η_max-insensitivity, no mesh order measured, although its test name carries "V3") | R10 (−U) rejected; wrong erfc(y/√νt) fails (8.20) |
| N56 stopped plate (ours) | NOTE | `stokes_first_stopped` | V1 | = (8.30) for t ≤ T (1e-15); u(0, t > T) = 0; diffusion residual → 0 at order 2; ∫u dy = 2U√(νt/π) (1e-9) | analytic | |
| C10 (8.32a,b), Examples 8.4–8.7 | CORE | `similarity_reduce_sympy`, `similarity_collapse_error`, `vortex_sheet_diffusion`, `transition_width`, `temporal_bl_wall_stress`, `line_vortex_decay`, `line_vortex_spinup` | V2, V1, V4 | D21, D22 (★★★, 14 steps), D23 re-derived; collapse spread < 1e-12 at (0, ½), (½, ½), (1, ½), (⅕, ⅕) and > 1e-2 off by 0.1; −∫ω dy = 2U (1e-10) at three t; ch05 parity with γ = −2U 1e-14; τ_w by FD 1e-7; C_f = 1.1284 Re_x^{−1/2}; Lamb–Oseen = Gaussian vortex σ = 2√(νt) (1e-14) and ch04 preset (1e-12) | symbolic, analytic, conserved | R11: 2.76 is not the 95 % point (2.7718); wrong γ = +2U, σ² = 2νt, n = ¼ |
| C11 (8.33)–(8.38) Stokes' second problem | CORE | `stokes_second_problem`, `stokes_layer`, `stokes_layer_state`, `stokes_layer_envelope`, `stokes_second_sympy` | V2, V1, V3 | D24 re-derived; wall U cos ωt 1e-15; envelope over a period 2e-5; crest speed by zero-crossing tracking = √(2νω) (1e-6); e^{−2√2} = 0.05911; CN space order **1.99**, time **2.00**; from rest after 10 periods 7.7e-5U within 6δ_e | symbolic, analytic, converged | wrong: growing root unbounded |
| C12 (8.39)–(8.43) Stokes equations | CORE | `low_re_scaling_sympy`, `stokes_residual` | V2, V1, V3 | D26 re-derived (dynamic {1, 1, 1/Re}, ×Re → pressure lost; viscous {Re, 1, 1}); Stokes sphere residual < 1e-6·3μUa/r³ on 50 random points, → 0 at order **2.0**; ideal-flow sphere fails (> 0.1 scale) | symbolic, analytic, converged | |
| C13 (8.44)–(8.49) stream function | CORE | `E2`, `E4_residual`, `stokes_sphere_sympy`, `stokes_sphere_streamfunction`, `stokes_sphere_velocity(_xyz)`, `side_line_speed` | V2, V1, V7 | D27, D28 (★★★: ∇×∇×(Ae_φ) = −E²(r sinθ A)/(r sinθ)e_φ for a generic A), D29 re-derived; roots {−1, 1, 2, 4}; constants; no slip 1e-15; far field 1e-5 at 1e6a; ∇·u (FD) < 1e-8; (6.83) by FD 1e-7; fluid frame fore–aft symmetric 1e-13; ideal side speed = `core.potential.sphere` 1e-12 | symbolic, analytic | R17: ∇⁴ψ ≠ 0 |
| C14 (8.50)–(8.52) drag, settling, Millikan | CORE | `stokes_sphere_pressure`, `_surface_stresses`, `stokes_drag`, `stokes_drag_running`, `sphere_drag_quadrature`, `stokes_drag_coefficient`, `terminal_velocity`, `radius_from_terminal_velocity`, `settling_state`, `millikan_charge`, `synthetic_millikan` | V2, V1, V5 | D30 (★★★, with `core.curvilinear.vector_laplacian`) and D31 (★★★, with `strain_rate`) re-derived: 2πμaU + 4πμaU; Gauss–Legendre exact for n ≥ 1 (1e-12); running integrals' derivative = integrand (1e-7); σ_rθ by FD 1e-4 (one-sided); C_D parity with ch04 1e-15; drag = weight 1e-14; 10 µm droplet U_t = 1.21 cm/s, Re = 0.016, D = 4.10e-11 N; synthetic e error **5.1e-4** (CODATA) | symbolic, analytic, benchmark (`sphere_drag_quadrature`: analytic — Gauss–Legendre exact from n = 1, nothing to converge) | R12: printed +3μU/2a minimum rejected; slip-sphere 4πμaU rejected |
| C15 far-field breakdown, Oseen (8.53) | CORE | `inertia_viscous_ratio`, `oseen_streamfunction`, `oseen_velocity`, `oseen_limit_sympy`, `oseen_linearisation_sympy`, drag laws | V7, V2, V1, V5 | log–log slope **0.995** (r/a 50–500, three θ); ∝ 1/ν exact, ∝ U 2e-5; loop 2: ratio/(Re_a r/a) → **½** on θ = 0, π/2, π (0.4981/0.4991 at r = 400a, O(a/r) approach, order 1.00; sympy limit exactly ½); D32 (asymptotics), D33 re-derived; Re → 0 limit 0 (sympy) and 1e-9 numerically; ψ = 0 on the axis; (6.83) by FD 1e-6; wall slip ∝ Re (order **0.996**); loop 2: coded ψ = (8.53) (1e-11) and (8.53) solves E⁴ψ − (U/ν)∂ₓE²ψ = 0 exactly (sympy; the opposite sign fails), wake downstream (\|ψ\| ratio 1.05–3.8 at θ < π/2 vs π − θ); Oseen/PP = radius forms 1e-14; Morrison between Stokes and Oseen for 0.1–5 | analytic, symbolic, benchmark (`inertia_viscous_ratio`: analytic; `oseen_streamfunction`: symbolic, analytic; `oseen_drag_coefficient`, `proudman_pearson_drag_coefficient`: analytic (form cross-check) — the V5 row checks the published *form* and a correlation bracket, no drag data are compared) | R13: printed +∂p/∂x differs by 2∂p/∂x; wrong: 3/8 on diameter Re; θ ↔ π − θ; ratio × ¼ |
| NOTE figures N103–N107 | NOTE | 12 scripts | V1 smoke | all exit 0 | — | |

C-depth NOTEs (N02, N11, N37, N64, N74–N76, N80, N94, N102) are named only and not coded; SKIP S01 not coded.

## Derivations (curation §4b) — every ★★/★★★ D row has a V2 test
| D | ★ | test | intermediate lines checked |
|---|---|---|---|
| D04 | ★★ | `test_couette_poiseuille_V2_derivation` | steps 1–9, book A = −c₁, check number −0.0125 m/s |
| D06 | ★★ | `test_pipe_poiseuille_V2_derivation_engine_and_D06` | steps 4, 6–10; check 0.25 m/s |
| D07 | ★★ | `test_pipe_V2_derivation_D07_flux_friction_factor` | steps 2–10 |
| D08 (+D09 ★) | ★★ | `test_circular_couette_V2_engine_and_D08_derivation` | steps 5–10; limits; Γ; Wikipedia form |
| D10 | ★★★ | `test_lubrication_scaling_V2_derivation` | steps 3–14 each as a coefficient ratio; printed (8.13b) |
| D11 | ★★ | `test_lubrication_balance_V2_derivation_and_units` | steps 2–7, units, Λ = 49.35 |
| D12 | ★★ | `test_lubrication_profile_V2_derivation` | steps 2–8 |
| D13 | ★★★ | `test_reynolds_equation_V2_derivation` | steps 1–12 with generic h(x, t), p(x, t) |
| D14 | ★★★ | `test_slider_bearing_V2_derivation` | steps 3–15 incl. step 12's numerator identity; exact-load series |
| D15 | ★★ | `test_thin_film_equation_V2_derivation` | steps 3–10, units |
| D17, D18, D19 | ★★ | `test_stokes_first_V2_derivation_D17_D18_D19` | Π groups, ∂η/∂t = −η/2t, (8.26), A = −1/√π, erfc, check 0.4795 |
| D21 | ★★ | `test_similarity_example_8_4_V2_derivation` | steps 2–8 |
| D22 | ★★★ | `test_vortex_sheet_V2_derivation` | steps 1–14; checks 2.772, 5.544 mm, −5.64 s⁻¹ |
| D23 | ★★ | `test_similarity_example_8_7_V2_derivation` | step 4, steps 7–9; wrong exponents leave t |
| D24 | ★★ | `test_stokes_second_V2_engine_and_D24_derivation` | steps 1–10 |
| D26 | ★★ | `test_low_re_scaling_V2_engine_and_D26_derivation` | steps 2–9; droplet 0.022 Pa vs 1.7e-4 Pa |
| D27, D29 | ★★ | `test_stokes_stream_function_V2_derivation_D29_D27` | D29 steps 3–10; D27 curl∇p = 0, curl∇² = ∇²curl |
| D28 | ★★★ | `test_stokes_stream_function_V2_derivation_D28` | steps 2–12 (identity for generic A, then (8.48)) |
| D30 | ★★★ | `test_stokes_pressure_V2_derivation_D30` | steps 1–11 (g′ = 0 checked) |
| D31 | ★★★ | `test_stokes_drag_V2_derivation_D31` | steps 2–13 and the running integral |
| D32 | ★★ | `test_far_field_V2_derivation_D32` + `test_inertia_viscous_ratio_V1_half_Re_a_r_asymptote` (sympy prefactor ½) | steps 1–7 (see O3) |
| D33 | ★★ | `test_oseen_V2_limit_linearisation_and_D33` | steps 1–7 |
The ★ rows of CORE items are checked numerically: D01 (O2), D02–D03 (engine), D05 (5 mm, 2 Pa/m, u_max = 1.5V), D09,
D16 (via the stokes1 reduction), D20 (3.643, 3.64 cm, 21.9 cm, 2.772), D25 (0.0591, 0.56 mm, 1.6 mm, 3.5 mm/s).
**No wrong intermediate line in the ★★★ derivations.** Design-text discrepancies: O2, O3, O4.

## Functions used by the notebook and explainers (design Part C) — test name each
C.1 `core.laminar`: `channel_flow` → `test_channel_flow_V1_navier_stokes_residual_and_walls`, `…_parity_with_ch04…`,
`…_V7_superposition…` · `channel_flow_rate` → `test_channel_flow_rate_V1_quadrature_and_printed_V_units` ·
`channel_shear_stress` → `test_channel_shear_stress_V1_derivative_and_wall_values` · `channel_backflow_threshold`,
`couette_poiseuille_state` → `test_channel_backflow_V1_threshold_by_bisection_and_state`,
`test_couette_poiseuille_state_V1_extrema_match_dense_grid` · `pipe_poiseuille` →
`test_pipe_poiseuille_V1_parity_ns_residual_and_ch03` · `pipe_shear_stress`, `pipe_wall_stress` →
`test_pipe_wall_stress_V4_control_volume_force_balance` · `pipe_flow_rate`, `pipe_friction_factor` →
`test_pipe_flow_rate_V1_quadrature_umax_and_friction_factor` · `circular_couette` → `test_circular_couette_V1_walls_ode_vorticity`,
`…_V7_limits_and_ch05_ch03_parity` · `circular_couette_pressure` → `test_circular_couette_pressure_V1_radial_balance` ·
`circular_couette_shear_stress`, `circular_couette_power` → `test_circular_couette_V4_torque_and_power_equal_dissipation` ·
`circular_couette_state` → `test_circular_couette_state_V1_rayleigh_criterion` · `similarity_variable`,
`stokes_first_problem` → `test_stokes_first_V1_boundaries_collapse_and_parity` · `stokes_first_vorticity`,
`stokes_first_state` → `test_vorticity_content_V1_plus_U_and_vorticity_field` · `stokes_first_stopped` →
`test_stokes_first_stopped_V1_superposition` · `diffusion_thickness` → `test_diffusion_thickness_V1_level_inversion_and_D20_numbers` ·
`transition_width`, `vortex_sheet_diffusion` → `test_vortex_sheet_V4_conserved_jump_and_ch05_parity` ·
`temporal_bl_wall_stress` → `test_temporal_bl_V1_wall_stress_and_cf` · `line_vortex_decay`, `line_vortex_spinup` →
`test_line_vortex_V1_parity_circulation_and_axis`, `test_line_vortex_V2_residuals` · `stokes_second_problem`,
`stokes_layer`, `stokes_layer_state`, `stokes_layer_envelope` → `test_stokes_second_V1_wall_envelope_phase_speed_and_D25` ·
`crank_nicolson_1d` → `test_crank_nicolson_V3_second_order_space_and_time`, `…_impulsive_start_and_stability`,
`test_stokes_second_V3_crank_nicolson_orders_and_transients`.
C.2 `core.lubrication`: `lubrication_scales`, `lubrication_term_magnitudes` → `test_lubrication_scales_V1_definitions_and_term_magnitudes` ·
`lubrication_velocity`, `lubrication_flux` → `test_lubrication_velocity_V1_walls_flux_and_printed_form` ·
`reynolds_pressure_1d` → `test_reynolds_pressure_1d_V1_exact_slider_and_V3_order_two` · `slider_bearing`,
`slider_bearing_load` → `test_slider_bearing_V1_quadrature_load_and_numbers`, `…_V7_reversal_and_ends` ·
`slider_optimum_taper` → `test_slider_optimum_taper_V5_san_andres` · `slider_bearing_state`, `slider_gap_velocity` →
`test_slider_bearing_state_V1_explainer_numbers`, `test_slider_backflow_both_ends` · `hele_shaw_velocity`, `hele_shaw_potential` (+ `hele_shaw_mean_velocity`) →
`test_hele_shaw_V1_velocity_potential_and_mean` · `thin_film_flux` (+ `thin_film_velocity`) →
`test_thin_film_V1_profile_walls_and_flux` · `thin_film_spread` → `test_thin_film_spread_V4_volume_conserved`,
`…_V5_huppert_shape_and_t_one_fifth` · `viscous_current_similarity`, `thin_film_state` →
`test_viscous_current_similarity_V1_pde_volume_front` (+ `viscous_current_eta_N` → `test_viscous_current_eta_N_V5_huppert`;
`hele_shaw_streamfunction_grid` → `test_hele_shaw_grid_V3_staircase_first_order`).
C.3 `core.creeping`: `stokes_residual` → `test_stokes_residual_V1_sphere_field_and_ideal_flow_fails` · `E2`,
`E4_residual`, `stokes_sphere_sympy` → `test_stokes_sphere_sympy_V2_engine`, `…_derivation_D28` ·
`stokes_sphere_streamfunction`, `stokes_sphere_velocity`, `stokes_sphere_velocity_xyz` →
`test_stokes_sphere_velocity_V1_walls_far_field_divergence`, `test_stokes_sphere_V7_frames_and_fore_aft_symmetry` ·
`stokes_sphere_pressure` → `test_stokes_pressure_V1_extremes_and_printed_minimum` · `stokes_sphere_surface_stresses`,
`stokes_drag`, `stokes_drag_running`, `sphere_drag_quadrature` → `test_stokes_drag_V1_quadrature_parts_and_slip_sphere` ·
`stokes_drag_coefficient` → `test_stokes_law_V1_form_and_drag_coefficient` · `oseen_drag_coefficient`,
`proudman_pearson_drag_coefficient` → `test_drag_laws_V5_oseen_proudman_pearson_morrison` · `terminal_velocity`,
`radius_from_terminal_velocity`, `settling_state` → `test_terminal_velocity_V1_force_balance_round_trip_and_warning` ·
`millikan_charge` → `test_millikan_V5_synthetic_experiment_recovers_e` · `inertia_viscous_ratio` →
`test_inertia_viscous_ratio_V7_linear_growth_in_r_and_Re`, `test_inertia_viscous_ratio_V1_half_Re_a_r_asymptote` ·
`oseen_streamfunction`, `oseen_velocity` → `test_oseen_V1_velocity_axis_wake_and_no_slip_order`,
`test_oseen_V2_streamfunction_solves_oseen_equation_and_wake_downstream` · `side_line_speed` → `test_stokes_sphere_V7_frames_and_fore_aft_symmetry`.
C.4 `ch08`: `pipe_flow_regime` → `test_pipe_flow_regime_V1_definition_and_band_edges`, `…_V7_scale_invariance` ·
`inertia_viscous_scales`, `diffusion_time` → `test_diffusion_time_V1_D01_numbers` · `momentum_diffusivity` →
`test_momentum_diffusivity_V5_air_water_from_published_property_laws` · `wall_bc_residuals` →
`test_wall_bc_residuals_V1_no_through_flow_and_no_slip` · `parallel_flow_sympy` → the three C02–C04 V2 tests ·
`advective_acceleration_check` → `test_circular_couette_V1_navier_stokes_residual_with_pressure` ·
`lubrication_nondim_sympy` → `test_lubrication_nondim_V2_coefficient_sets_and_printed_8_13b` · `reynolds_equation_sympy` →
`test_reynolds_equation_V2_derivation` · `slider_bearing_sympy` → `test_slider_bearing_V2_engine_and_printed_slips` ·
`hele_shaw_cylinder` → `test_hele_shaw_cylinder_V1_equals_ideal_cylinder` · `similarity_reduce_sympy` →
`test_similarity_reduce_V2_all_cases` · `similarity_ode_solve` → `test_similarity_ode_solve_V3_bvp_matches_closed_forms` ·
`vorticity_content` → `test_vorticity_content_V1_plus_U_and_vorticity_field` · `similarity_collapse_error` →
`test_similarity_collapse_V1_right_exponents_only` · `stokes_second_sympy` → `test_stokes_second_V2_engine_and_D24_derivation` ·
`low_re_scaling_sympy` → `test_low_re_scaling_V2_engine_and_D26_derivation` · `oseen_linearisation_sympy`,
`oseen_limit_sympy` → `test_oseen_V2_limit_linearisation_and_D33` · `synthetic_millikan` (+ `estimate_elementary_charge`) →
`test_millikan_V5_synthetic_experiment_recovers_e` · `stokes_first_pi_groups` → `test_stokes_first_V2_derivation_D17_D18_D19` ·
`G0`, `G_BOOK`, `P_ATM` and scalar-callability of every `*_state` dict and 12 explainer-mirrored functions →
`test_part_c_V1_every_contract_function_exists_and_is_scalar_callable` (98 names).
C.0 reused, called with ch08 parameters: `ns_incompressible_terms`, `exact_solution(_fields)` (couette, poiseuille,
pipe_poiseuille, stokes_first, lamb_oseen, solid_body, cylinder), `viscous_force_forms`, `reynolds_number`,
`nondimensional_ns_coefficients` (through `low_re_scaling_sympy`), `sphere_drag_coefficient` ("stokes", "morrison"),
`ftcs_diffusion_1d`, `stable_time_step`, `couette_startup_profile` (→ `test_reuse_V1_couette_startup_ftcs_and_navier_stokes_presets`),
`gaussian_vortex`, `solid_body_rotation`, `line_vortex`, `core.curvilinear` (`curl`, `gradient`, `vector_laplacian`,
`advective_acceleration`, `strain_rate`), `core.potential.sphere`, `core.laplace_solvers` (through
`hele_shaw_streamfunction_grid`), `ch01.pi_groups`, `water_viscosity`, `water_density`, `sutherland_viscosity`,
`ch03.pipe_profile`, `ch04.stokes_first_problem`, `plane_poiseuille`, `ch05.rotating_cylinder_flow`,
`diffusing_vortex_sheet`, `dissipation_outside_cylinder`, `tools.convergence.observed_order`. The machinery rows 0.1–0.4
(`style`, `anim`, `interact`, `embed`) carry no physics and are covered by `tests/test_machinery.py`.

## Convergence studies (scheme | grids | observed order | design order)
| scheme | grids | observed | design |
|---|---|---|---|
| Crank–Nicolson (8.20), time, smooth start at t = 100 s | Δt = 30, 15, 7.5, 3.75 s (801 nodes) | 2.000 (pairwise 2.000, 2.000, 2.000) | 2 |
| backward Euler (θ = 1), same (wrong variant) | same | 0.996 | 1 |
| Crank–Nicolson (8.20), space | 41, 81, 161, 321 nodes on 12√(νt) | 1.995 | 2 |
| CN impulsive start with two BE steps (Rannacher) | Δt = T/50 … T/400, Δy ∝ Δt | 1.99 (pure CN stalls: last pairwise < 1) | 2 |
| Crank–Nicolson (8.33) oscillating wall, space | 61 … 481 nodes on 30δ_e | 1.986 | 2 |
| Crank–Nicolson (8.33), time | T/25 … T/200 | 2.00 | 2 |
| steady Reynolds pressure, trapezoid (array h) | 21 … 321 nodes | 2.000 | 2 |
| Hele-Shaw masked five-point Laplace (staircase disk) | 33, 65, 129, 257 per side | 1.11 (pairwise 1.05, 1.12, 1.16) | 1 (staircase; O5) |
| Stokes residual stencil (8.43) on the sphere field | h = 1e-2, 5e-3, 2.5e-3 | 2.0 | 2 |
| (4.39b) residual of circular Couette (Cartesian stencils) | h = 4, 2, 1 × 10⁻⁴ m | 2.0 | 2 |
| net viscous force of the viscous vortex (R10) | same | 2.0 (→ 0) | 2 |
| `stokes_first_stopped` diffusion residual | Δy = 4, 2, 1 × 10⁻⁴ m | 2.0 | 2 |
| Oseen wall slip vs Re (not a scheme: O(Re) no slip, R18) | Re = 1e-3, 1e-2, 1e-1 | 0.996 | 1 |
| `solve_bvp` (8.26) / Example 8.6 ODE | η_max = 10, 12, 14 | max error 3.7e-13 / 2.5e-9 | ≤ 1e-7 |
| thin film vs Huppert (backward Euler + Newton, precursor) | N = 200, 400, 800 | L2 0.4–1.1 % (not monotone: front singularity) | < 2 % |

## Conservation / invariant residuals
Volume of `thin_film_spread`: 2.2e-16 relative over 13 outputs to 1000 s · vortex-sheet jump −∫ω dy − 2U: 1e-10
relative at t = 0.5, 50, 5000 s · Stokes' first problem ∫ω dy − U: 1e-10 · pipe CV force balance: 1e-15 N ·
circular-Couette torque 2πR²σ_Rφ across the gap: 1e-13 relative; power in − dissipation: 1e-12 (closed) and 1e-10
(quad); inner + outer torque: 0 to 1e-18 · Millikan: synthetic drops recover every integer n.

## Benchmarks used (value | our value | source + DOI | date verified)
| value | ours | source | verified |
|---|---|---|---|
| slider K_opt = 2.1889 | 2.18870 | San Andrés, Modern Lubrication Notes 2 Appendix (TAMU 2009/2012) | 2026-09-28 |
| slider W(K_opt) = 0.0267 | 0.026707 | same | 2026-09-28 |
| peak-pressure K = 2.414, P_max = 0.043 | 2.4142, 0.04289 | same | 2026-09-28 |
| planar gravity current η_N = 1.411 | 1.411245 | Ball & Huppert preprint App. (a); Huppert, JFM 121, 43 (1982) | 2026-09-28 |
| e = 1.602176634e-19 C (exact) | synthetic estimate within 5.1e-4 | NIST CODATA 2022 | 2026-09-28 |
| Morrison sphere C_D (correlation) | lies between Stokes and Oseen for 0.1 ≤ Re ≤ 5 | Morrison (2016), reused from `reference/ch04` | 2026-09-23 |
| Sutherland β, S (USSA-1976) | ν_air reproduced to 1e-12 | NASA-TM-X-74335, reused from `reference/ch01` | 2026-09-12 |
| forms: Oseen/Proudman–Pearson, Stokes' law, Stokes' second problem, Hagen–Poiseuille, Taylor–Couette A, B | identical (1e-14 or symbolic) | Wikipedia pages listed in `reference/ch08/SOURCES.md` | 2026-09-28 |

## Numbers from the text (book vs ours)  [private values redacted to percentages]
Three V6 tests (`test_book_V6_*`, skipped without `tests/book_values_ch08.json`) reproduce every printed closed form of
§8.2–§8.6 in the JSON to 1e-12 at random parameters, and: the §8.1 air/water diffusivity ratio to 0.05 %, ν_air and
ν_water to < 1 %; the §8.3 oil-film reduced Reynolds number exactly (0.0 %); the (8.31) thickness factor to 0.08 %
(rounds to the printed value); Example 8.5's width factor to 0.07 %; the §8.5 amplitude to its one printed significant
figure (the exact value is 1.5 % below the rounded print); the §8.6 front-stagnation pressure exactly, the drag and its
fractions exactly. Printed slips confirmed as slips (asserted as wrong variants): Example 8.5's η₉₅ (0.43 % off — its
own width is 2 × 2.772); ∫ω dy printed −U; the rear pressure minimum printed with +; Example 8.1's first-power
denominator (differs from the exact pressure by > 0.1 % at α = 0.37 and violates the ODE); (8.19)'s U₀ placement;
(8.13b)'s ∂p/∂x; (8.17a)'s missing ν; the channel V without 1/h; Oseen's +∂p/∂x_i.

## Figures reproduced with our code (figure | outputs/ch08/verify/<png> (local) | one-sentence visual verdict)
| figure | png | verdict |
|---|---|---|
| Fig. 8.4 | `fig8_4_couette_poiseuille.png` | favourable profile fuller than Couette, adverse profile reverses in the lower half and crosses zero at y = 5 mm (dp/dx = 2× threshold), Poiseuille parabola peaks at mid-gap, and its \|τ\| is a V with zero at mid-gap — as expected. |
| Figs. 8.5–8.6 | `fig8_5_8_6_pipe_circular_couette.png` | paraboloid with u(0) = 0.25 m/s and u(a) = 0, τ linear from 0 to −0.5 Pa at the wall; every circular-Couette curve starts at Ω₁R₁ = 10 mm/s and ends at Ω₂R₂, rigid rotation is a straight line. |
| Fig. 8.9 / Ex. 8.1 | `fig8_9_slider_bearing.png` | exact humps vanish at both ends with the peak shifted toward the exit (x = L/(2 + α)), the O(α) parabolas are symmetric and overshoot, the printed first-power ghost differs; W(K) peaks at K = 2.1887, W* = 0.0267. |
| Fig. 8.11 / Exs. 8.3, 8.7 | `fig8_11_spreading.png` | numerical bead overlays Huppert's (1 − x²/x_N²)^{1/3} cap at 32 s and 1000 s, and the front approaches the t^{1/5} line from above once the initial box is forgotten. |
| Figs. 8.12–8.13 | `fig8_13_stokes_first.png` | three dimensional profiles thicken like √t and collapse onto one curve in y/(2√(νt)); CN errors fall with slope 2, backward Euler with slope 1. |
| Figs. 8.14–8.15 | `fig8_14_8_15_vortex_diffusion.png` | the Gaussian vorticity layer halves its peak and doubles its width from 5 s to 20 s; the line-vortex profile has a solid-body core and a 1/r tail whose peak moves out and down as νt grows. |
| Fig. 8.16 | `fig8_16_oscillating_plate.png` | four phase profiles start at ±U or 0 at the wall, stay inside the ±e^{−y/δ_e} envelope and are negligible beyond the book depth 2.83δ_e. |
| Fig. 8.17 | `fig8_17_sphere_surface.png` | p a/(μU) runs from −1.5 at the rear (θ = 0) to +1.5 at the front, σ_rθ = −1.5 sin θ, and the x-traction is the uniform 1.5. |
| Figs. 8.19–8.20 | `fig8_19_8_20_stokes_oseen.png` | the Stokes fluid-frame streamlines are fore–aft mirror images, the Oseen (Re = 1) ones are asymmetric (a wake). |
| C_D(Re) | `fig_drag_curve.png` | all laws merge on 24/Re below Re ≈ 0.1; the Morrison correlation lies between Stokes and Oseen up to Re ≈ 5; Proudman–Pearson turns up beyond its range (Re ≳ 3). |

## Deviations & justifications (every `# DEVIATION` in the code)
- `core.laminar.channel_flow_rate`: V = Q/h (the printed middle form omits 1/h, R9) — verified by pint and by quad.
- `core.laminar.circular_couette_power`: power into the fluid = −2πR₁σ_Rφu_φ > 0 (R15) — verified equal to ∫ε dA.
- `core.lubrication.lubrication_velocity`: U₀(1 − y/h) (R8b) — the printed form fails u(h) = U_h when U₀ ≠ 0.
- `core.lubrication.slider_bearing`: squared denominator (R8) — verified by D14 step by step and by the Reynolds solver.
- `core.lubrication.thin_film_spread`: Newton instead of Picard (same fixed point) — volume conserved to 2e-16, Huppert
  shape reached to < 2 %.
Other conventions (not marked DEVIATION, verified): (8.50) coded as printed, only the prose minimum is corrected (R12);
Oseen evaluated in the −expm1 form (exact at Re = 0).

## Discrimination proofs (wrong variants planted in a scratch copy; the repository untouched)
Method: a scratch copy of `fluidpy/`, `tools/`, `tests/test_ch08.py` (+ conftest, private JSON, `reference/ch01`, `ch04`,
`ch08`); one wrong variant planted at a time by exact string replacement; `pytest -x -m "not slow"`; the first failing
test recorded. **36/36 caught in loop 1; loop 2 adds rows 37–41 (run with `-k` on the new test): 41/41 caught.**

| # | wrong variant | where | first test that fails |
|---|---|---|---|
| 1 | G passed as dp/dx (sign of the parabola) | `core.laminar.channel_flow` | `test_channel_flow_V1_navier_stokes_residual_and_walls` |
| 2 | Q with h³/(6μ) instead of h³/(12μ) | `channel_flow_rate` | `test_channel_flow_rate_V1_quadrature_and_printed_V_units` |
| 3 | pipe profile with 1/(2μ) | `pipe_poiseuille` | `test_pipe_poiseuille_V1_parity_ns_residual_and_ch03` |
| 4 | Hagen–Poiseuille πa⁴/(4μ) | `pipe_flow_rate` | `test_pipe_flow_rate_V1_quadrature_umax_and_friction_factor` |
| 5 | circular-Couette B with the wrong sign | `_cc_coeffs` | `test_circular_couette_V1_walls_ode_vorticity` |
| 6 | power with the printed sign (R15) | `circular_couette_power` | `test_circular_couette_V4_torque_and_power_equal_dissipation` |
| 7 | erfc(y/√(νt)) (missing 2) | `stokes_first_problem` | `test_stokes_first_V1_boundaries_collapse_and_parity` |
| 8 | vortex-sheet ω with + (γ = +2U) | `vortex_sheet_diffusion` | `test_vortex_sheet_V2_derivation` |
| 9 | Lamb–Oseen with σ² = 2νt | `line_vortex_decay` | `test_line_vortex_V1_parity_circulation_and_axis` |
| 10 | Stokes layer with phase ωt + y/δ_e | `stokes_second_problem` | `test_stokes_second_V1_wall_envelope_phase_speed_and_D25` |
| 11 | δ₉₉ on the figures' η/2 axis (R22) | `diffusion_thickness` | `test_diffusion_thickness_V1_level_inversion_and_D20_numbers` |
| 12 | wall vorticity without √π | `stokes_first_vorticity` | `test_vorticity_content_V1_plus_U_and_vorticity_field` |
| 13 | C_f without the ½ | `temporal_bl_wall_stress` | `test_temporal_bl_V1_wall_stress_and_cf` |
| 14 | circular-Couette pressure 2AB → AB | `circular_couette_pressure` | `test_circular_couette_V1_navier_stokes_residual_with_pressure` |
| 15 | (8.19) with the printed U₀ placement (R8b) | `lubrication_velocity` | `test_lubrication_velocity_V1_walls_flux_and_printed_form` |
| 16 | slider pressure with the printed first-power denominator (R8) | `slider_bearing` | `test_reynolds_pressure_1d_V1_exact_slider_and_V3_order_two` |
| 17 | thin-film flux h³/(2μ) | `thin_film_flux` | `test_thin_film_V1_profile_walls_and_flux` |
| 18 | solver with β = ρg/(2μ) | `thin_film_spread` | `test_thin_film_spread_V5_huppert_shape_and_t_one_fifth` |
| 19 | η_N with exponent −1/2 | `viscous_current_eta_N` | `test_viscous_current_similarity_V1_pde_volume_front` |
| 20 | Re_L with the gap h (analysis wrong variant) | `lubrication_scales` | `test_lubrication_balance_V2_derivation_and_units` |
| 21 | Reynolds pressure with half the Couette flux | `reynolds_pressure_1d` | `test_reynolds_pressure_1d_V1_exact_slider_and_V3_order_two` |
| 22 | Hele-Shaw p = −6μφ̄/h² | `hele_shaw_cylinder` | `test_hele_shaw_cylinder_V1_equals_ideal_cylinder` |
| 23 | Hele-Shaw φ with + sign | `hele_shaw_potential` | `test_hele_shaw_V1_velocity_potential_and_mean` |
| 24 | Stokes pressure +3μaU cos θ/2r² (R12) | `stokes_sphere_pressure` | `test_stokes_residual_V1_sphere_field_and_ideal_flow_fails` |
| 25 | u_θ with the wrong sign | `stokes_sphere_velocity` | `test_stokes_residual_V1_sphere_field_and_ideal_flow_fails` |
| 26 | slip-sphere drag 4πμaU | `stokes_drag` | `test_stokes_drag_V1_quadrature_parts_and_slip_sphere` |
| 27 | Oseen 3/8 on the diameter Re | `oseen_drag_coefficient` | `test_drag_laws_V5_oseen_proudman_pearson_morrison` |
| 28 | surface shear +(3μU/2a) sin θ | `stokes_sphere_surface_stresses` | `test_stokes_drag_V1_quadrature_parts_and_slip_sphere` |
| 29 | settling speed with ρ′ + ρ | `terminal_velocity` | `test_stokes_law_V1_form_and_drag_coefficient` |
| 30 | Millikan balance with the weight subtracted | `millikan_charge` | `test_millikan_V5_synthetic_experiment_recovers_e` |
| 31 | Oseen ψ with e^{−s} for (1 − e^{−s})/s | `oseen_streamfunction` | `test_oseen_V1_velocity_axis_wake_and_no_slip_order` |
| 32 | inertia/viscous ratio with ν² | `inertia_viscous_ratio` | `test_inertia_viscous_ratio_V7_linear_growth_in_r_and_Re` |
| 33 | pipe regime upper edge 2300 | `ch08.pipe_flow_regime` | `test_pipe_flow_regime_V1_definition_and_band_edges` |
| 34 | air μ instead of ν | `ch08.momentum_diffusivity` | `test_momentum_diffusivity_V5_air_water_from_published_property_laws` |
| 35 | θ-scheme without its explicit half | `core.diffusion.crank_nicolson_1d` | `test_crank_nicolson_V3_second_order_space_and_time` |
| 36 | Stokes-layer envelope √(ω/ν) (missing 2) | `stokes_layer_envelope` | `test_stokes_second_V1_wall_envelope_phase_speed_and_D25` |
| 37 | extremum at y* = h/2 + μU/(h dp/dx) (review M1, the pre-fix code) | `couette_poiseuille_state` | `test_couette_poiseuille_state_V1_extrema_match_dense_grid` |
| 38 | Oseen ψ mirrored θ ↔ π − θ (upstream wake; review SF1) | `oseen_streamfunction` | `test_oseen_V2_streamfunction_solves_oseen_equation_and_wake_downstream` (and `test_oseen_V1_velocity_axis_wake_and_no_slip_order`) |
| 39 | inertia/viscous ratio × ¼ (the loop-1 "≈ 1/8" prefactor) | `inertia_viscous_ratio` | `test_inertia_viscous_ratio_V1_half_Re_a_r_asymptote` |
| 40 | recirculation checked at x = L only | `slider_bearing_state` | `test_slider_backflow_both_ends` |
| 41 | recirculation threshold h > h_m instead of 1.5h_m | `slider_bearing_state` | `test_slider_backflow_both_ends` |

The printed slips that live only in the sympy engines (R6 (8.13b), R7 (8.17a), R8's (1 − αx/L) integrands, R13 Oseen's
sign, R17 the biharmonic) are asserted as failing variants directly inside the V2 tests (the engines return them as
separate keys), and the independent re-derivations do not use the engines.

## Open items
O1 — **closed (loop 2).** `grep "Validation (planned)" fluidpy/` returns 0; every ch08 docstring names its tests
(`Validation — tests/test_ch08.py: …`). Labels now read as in the validation table: `thin_film_spread` "conserved,
benchmark", `sphere_drag_quadrature` "analytic", `oseen_drag_coefficient` / `proudman_pearson_drag_coefficient`
"analytic (form cross-check)", `similarity_ode_solve` "analytic (closed-form agreement, truncation-insensitive)".

O2 (designer, Part F D01 check) — "The property functions give 15.1": `momentum_diffusivity("air")/("water")` at
293.15 K gives **15.007** (15.0). The step-4 hand value (0.018 × 833 = 15) is fine; only the check line's 15.1 is off.

O3 (designer, Part F D32) — step 3 sizes the inertia with u_r ∂u_θ/∂r; in (8.49) the θ-component of (u·∇)u has no
O(U²a/r²) term (its leading term is (9/32)U²a² sin 2θ/r³), the O(U²a/r²) inertia is in the **radial** component:
(u·∇u)_r ≈ (3U²a/16r²)(8cos²θ − 4sin²θ) = (3U²a/4r²)(2 − 3sin²θ). That θ-part of the finding stands.
**Corrected in loop 2:** loop 1 compared this with the radial viscous component only and concluded "prefactor ≈ 1/8,
crossover 10–20 a/Re_a" — wrong. The viscous term is a vector too: |ν∇²u| = |∇p|/ρ = (3νUa/2r³)(4cos²θ + sin²θ)^{1/2},
so on the axis ratio = (1.5U²a/r²)/(3νUa/r³) = ½ Re_a (r/a), at θ = π/2 (0.75U²a/r²)/(1.5νUa/r³) = ½ Re_a (r/a) as
well, and the crossover there is **r/a ≈ 2/Re_a** (smaller in between: the leading inertia vanishes near
sin²θ = 2/3). Pinned by `test_inertia_viscous_ratio_V1_half_Re_a_r_asymptote` (0.4981 / 0.4991 at r = 400a, O(a/r)
approach, sympy limit ½). The design's check line "Re_a = 0.01 → r ≈ 100a" is therefore low by a factor 2 on the axis
(≈ 200a); still an order-of-magnitude statement and should be worded as one. Suggested wording "u·∇u ~ U ∂u′/∂x" (the
stream advecting the Stokeslet), which is what Oseen keeps (N95).

O4 (designer, consistency) — the engine-film ν: design §0 item 10 says ν = 10⁻⁴ m²/s with μ = 0.05 Pa s, Part A/B use
ν = 5.75 × 10⁻⁵ (ρ = 870) → ε²Re_L = 4.35 × 10⁻³, Λ = 49; the curation N31 uses ν = 10⁻⁴ → ε²Re_L = 2.5 × 10⁻³, Λ ≈ 10³
(with a longer L). Both are correct for their inputs; the notebook should quote one set. D27's check cites
`stokes_sphere_sympy()["curlcurl_identity"]` for −∇×∇×ω = 0; that key holds ∇p + μ∇×ω (the D30 identity). −∇×∇×ω = 0
itself is verified here with `core.curvilinear.curl` (0 for Stokes' ω).

O5 (analyst plan vs implementation, low) — analysis §6 V15 planned the Hele-Shaw grid solve "order ≈ 2 away from the
boundary (±0.25)"; `hele_shaw_streamfunction_grid` measures the error where r ≥ 2a and its staircase disk makes the whole
field first order: observed **1.11** (its docstring says "between 1 and 2"). The test asserts order 1 ± 0.15. This is our
numerical cross-check route, not the book's method (N36 is verified exactly by the analytic cylinder), so it is not a
physics defect; a cut-cell or body-fitted boundary would be needed for order 2.

O6 (runtime, low) — `tests/test_ch08.py` without the slow script test takes ≈ 75 s (target < 60 s): four sympy
engines take 25 s. Kept rather than marked slow, because they are the ★★★ evidence.

O7 (book value) — none outstanding: every printed number of the private file is reproduced (to printed precision where
the book rounds) or rejected as a documented slip (R8–R13). The §8.5 "0.06" is the 1-significant-figure rounding of
0.0591.

O8 (implementer, docstrings, low) — the four loop-2 tests are not yet named in the docstrings they verify:
`couette_poiseuille_state` (→ `test_couette_poiseuille_state_V1_extrema_match_dense_grid`), `oseen_streamfunction`
(→ `test_oseen_V2_streamfunction_solves_oseen_equation_and_wake_downstream`; its "Checks" line could now say "V2 solves
Oseen's equation"), `inertia_viscous_ratio` (→ `test_inertia_viscous_ratio_V1_half_Re_a_r_asymptote`), and
`slider_bearing_state`, whose docstring still says "test_slider_backflow_both_ends (to be added by the verifier …)" —
it exists now. Text only.

O9 (implementer / notebook / explainer E3, low) — `slider_bearing_state["inlet_backflow"]` now means "recirculation at
the wide end" (x = L for α > 0, x = 0 for α < 0) and equals `backflow_any`; for α < 0 with U > 0 the wide end is the
exit, so the key's name is misleading there. Existing callers that assume α > 1 only still get the right answer for
α > 0. Any text that says "inlet backflow" should say "recirculation at the wide end". The docstring's η in
"(u − U)/U = −(1 − η) + 3(1 − h_m/h)η(1 − η)" is measured from the floor, while `slider_gap_velocity`'s y = 0 is the
moving pad side (u(0) = U in the ground frame): consistent physics, but the two descriptions put the pad on
opposite walls — worth one clarifying phrase.

## Verdict: PASS (loop 2) — 101/101 ch08 tests pass (169 s; full suite 983 passed in 660 s); all 15 CORE items have ≥ 2 independent levels including V1 or V2; every coded NOTE ≥ 1; all 25 ★★/★★★ derivations re-derived with sympy (no wrong intermediate line in the ★★★ rows; design-text notes O2–O4, O3 corrected: prefactor ½, crossover ≈ 2a/Re_a); review M1 (Couette–Poiseuille extremum sign), SF1 (Oseen wake direction) and the slider recirculation at either end are pinned by new tests; every design Part C function and all 12 scripts exercised; 41/41 planted wrong variants caught; no tolerance loosened below a derived bound; O1 closed; open items are docstring/naming text (O8, O9), design-text notes, a numerical-route order and runtime.
