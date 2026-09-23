# Chapter 4 verification — Conservation Laws                     2026-09-23, commit ca36949 + working tree (loop 2: implementer fix of F1, O2, O4, O6; designer fix of O3)

Verifier: math-verifier. Suite: `tests/test_ch04.py` (140 test functions = 147 collected items; the exact-solution
test is parametrised over 8 solutions), figures/metrics: `tests/ch04_verify_figures.py` → `outputs/ch04/verify/`
(git-ignored), cited data: `reference/ch04/` (`make_refs.py`, `benchmarks.json`, `SOURCES.md`). Nothing in `fluidpy/`
or `scripts/` was edited by the verifier.

## Environment
Python 3.11.5 · numpy 2.4.6 · scipy 1.17.1 · sympy 1.14.0 · pint 0.25.3 · matplotlib 3.11.2 (Windows, `.venv`).

Runs (loop 2, final): `pytest tests/test_ch04.py -q -p no:cacheprovider` → **147 passed in 57 s**, no warnings (with
the private book file; without it **144 passed, 3 skipped** — the three V6 tests). Full suite `pytest -q` → **474
passed in 115 s** (ch01–ch03 + machinery 327 items). All 8 `scripts/ch04_*.py` exit 0 (3–5 s each).
`tools/check_public.py` → OK.

## Loop history
| loop | result | what changed |
|---|---|---|
| 1 | FAIL: 145/146 pass | F1: `couette_heating_transient` truncated its sine series at the first (symmetry-)zero coefficient and aliased the high modes with a fixed 128-node quadrature — t = 0 error 2.98e-3 ΔT_max (dp/dx = 0) and 1.80 ΔT_max (dp/dx = −2000 Pa/m). O2–O6 opened. |
| 2 | **PASS: 147/147** | Implementer: b_n in closed form (recursion for ∫y^m sin(nπy/h) dy), all modes summed; O2 `lamb_vector` docstring sign; O4 FTCS claim removed; O6 quad epsrel 1e-11 (no warning). Designer: O3 contract numbers 3.3661 m/s, 2.7269 mm. Verifier: new independent test `test_couette_heating_V1_transient_equals_independent_quad_series` and 5 new discrimination variants of the series. |

F1 re-checked on our own terms (loop 2): b_n for n = 1…40 by adaptive `quad` of (2/h)∫(T_s − T₀) sin(nπy/h) dy
(independent of the implementer's recursion) rebuild T(y, t) that equals the function to < 1e-9 ΔT_max at
t = 0.002, 0.02, 0.2 h²/κ for dp/dx = 0, −2000, +5000 Pa/m; even modes vanish (< 1e-12) for pure Couette; t = 100h²/κ
gives the steady profile (1e-12); t = 0 error 2.1e-7 (dp/dx = 0) and 6.4e-7 (−2000) of ΔT_max at nterms = 200, falling
at least as nterms⁻² (4.3e-5, 1.3e-5, 1.7e-6, 2.1e-7 for 25, 50, 100, 200 — pairwise orders 1.8, 3.0, 3.0; at the
grid points the n⁻³ coefficients' alternating tail converges faster than the n⁻² bound); interior PDE residual < 1e-4 of φ at t = 0.2h²/κ.

Collected items per evidence tag (tag of the `def` line; many tests carry more levels inside): **V1 79 · V2 41 · V3 9 ·
V4 6 · V5 6 · V6 3 · V7 3** (= 147).

**Discrimination proofs** (skill lesson ch01): twenty-two wrong variants (17 in loop 1, 5 in loop 2) were patched one at a time into a scratch copy of
`fluidpy/` (scratchpad `discriminate.py`; the repository was not touched) and the ch04 suite re-run — every variant
fails at least one test (loop-1 counts exclude the then-failing F1 test):

| wrong variant patched in | tests that fail |
|---|---|
| Coriolis without the factor 2 (kinematic term and force) | 2 |
| Cauchy with the second index of τ (∂τ_ij/∂x_j, the book's prose) | 1 |
| 2-D stream function with the opposite sign (u = −∂ψ/∂y) | 3 |
| viscous term with the wrong sign in (4.39b) | 8 |
| centrifugal = +Ω × (Ω × x) (the centripetal sign) | 2 |
| (4.74) gauge with the printed + sign | 1 |
| plane-flow dissipation without the (3, 3) deviator entry | 1 |
| Lamb vector ω × u instead of u × ω | 1 |
| CV outflux with u·n instead of (u − b)·n | 6 |
| wake drag ρ∫U∞(U∞ − U) (deficit weighted by U∞) | 3 |
| Stokes assumption as λ = 0 (instead of μ_v = 0) | 5 |
| Bernoulli head U² instead of ½U² | 3 |
| kinematic condition ∂η/∂t − u·∇η | 2 |
| Boussinesq buoyancy with the wrong sign | 1 |
| meniscus height with cos θ instead of sin θ | 2 |
| (4.29) third term γδ_imδ_jn (duplicating the μ term) | 1 (after adding the component check; the first run had 0 — every other test contracts K with a symmetric S, where only μ + γ acts) |
| energy budget heat flux with the wrong sign | 1 |
| *loop 2:* Couette series truncated after 5 modes (the loop-1 defect) | 2 |
| *loop 2:* I_m recursion without (−1)ⁿ | 2 |
| *loop 2:* J_m recursion with the wrong sign | 2 |
| *loop 2:* decay rate κ(nπ/h)² with n not squared | 2 |
| *loop 2:* b_n with 1/h instead of 2/h | 2 |

## F1 (loop 1, resolved in loop 2) — kept for the record
**F1 — `ch04.couette_heating_transient` did not start from T₀ and was corrupted by aliasing** (`fluidpy/ch04_conservation_laws.py`, lines 676–705).
- Test: `test_couette_heating_V1_transient_starts_from_T0_at_documented_accuracy` (its docstring claims "t = 0 is
  reproduced to ~1/nterms² relative", i.e. ≈ 2.5e-5 for the default nterms = 200; the test allows 10/nterms² = 2.5e-4).
- Input: `y = linspace(0, 1e-3, 11)`, `t = 0`, U = 1 m/s, h = 1e-3 m, μ = 1e-3 Pa s, k = 0.6 W/(m K), ρ = 1000, C_p = 4182.
- Expected: T(y, 0) − T₀ = 0 (to ~2.5e-5 ΔT_max). Got: max|T − T₀|/ΔT_max = **2.98e-3** for dp/dx = 0 and **1.80**
  (i.e. 180 %) for dp/dx = −2000 Pa/m. At t = 0.001 s the dp/dx ≠ 0 profile is still off by 3e-3 of ΔT_max; by
  t = 0.2h²/κ both cases are fine (the PDE-residual and long-time tests pass).
- Evidence and cause (two defects, confirmed by calling the function with `tol` changed):
  1. *Early break on a zero coefficient.* `if n > 4 and abs(term) < tol*scale: break` fires at n = 6 when dp/dx = 0,
     because T_s − T₀ ∝ y(h − y) is symmetric about h/2 and every even-n coefficient b_n is exactly 0. The series is
     truncated after n = 5 (error 3e-3, visible as wiggles in `outputs/ch04/verify/couette_heating_verify.png`, left).
  2. *Aliasing of the projection.* b_n is computed with a fixed 128-node Gauss–Legendre rule, which cannot resolve
     sin(nπy/h) for n ≳ 100. With the break disabled (`tol=0`) the error at t = 0 is 1.3e-5 (nterms = 50), 1.7e-6 (100),
     but **1.95 ΔT_max** (200, the default) — high-n coefficients are garbage. With dp/dx ≠ 0 the break never fires, so
     the default already sums the aliased terms (the 180 % error; right panel of the same figure).
- Suggested fix (implementer's choice): compute b_n in closed form (T_s − T₀ is a quartic polynomial in y, so
  b_n = (2/h)∫(T_s − T₀) sin(nπy/h) dy has an exact expression), or scale the quadrature with nterms (≥ 2·nterms nodes);
  and stop the series on |b_n| ≤ tol over **two consecutive** n (or use the n⁻³ bound), not on one zero coefficient.
  The docstring's "V3 FTCS agrees" has no FTCS comparison behind it anywhere (O4).

## Validation table — A items (CORE, ≥ 2 independent levels, one of V1/V2/V3/V5)
| Concept / Eq. | Tier | fluidpy target | Evidence (V-levels) | Numbers (error, order, residual) | Label | Notes |
|---|---|---|---|---|---|---|
| C01 mass for a moving CV (4.5); N01, N03–N06, R01; D01 | CORE | `interval_mass_budget`, `mass_budget`, `MassBudget`, `material_mass`, `material_interval`, `expanding_flow(_fields)` | V1 contract numbers (−1, 1, 2, residual 0; material ends all 0), two arbitrary moving intervals = Leibniz (1e-12); callable route (quad + FD) = preset (1e-8), a non-uniform compressible pair closes (1e-7), a non-conserving pair leaves exactly ∫(ρ_t + (ρu)_x); V4 fixed box: storage = analytic dM/dt (1e-7), outflux = −storage (1e-12); material CV (b = u): outflux 0, storage 0 (box and sphere); material mass constant for t ∈ [0, 5] (ptp < 1e-12) and ends on `pathline` (1e-9); V3 moving + growing sphere, residual vs dt; V2 D01 in 1-D (generic density, arbitrary moving ends; (4.2)–(4.5); the u·n slip rejected) | order 2.023 (2.055, 2.014, 2.003) | analytic, conserved, converged, symbolic | |
| C02 continuity (4.7)–(4.10); N07–N12, R02; D02, D03 | CORE | `continuity_residual(_sym)`, `continuity_terms`, `continuity_material_terms`, `density_material_rate`, `divergence_free_check`, `stratified_shear_flow(_fields)`, `is_incompressible_regime`, `ball_integral` | V2 expanding flow 1-, 2-, 3-D ≡ 0, wrong exponent ≠ 0; V1 40 points × 3 times (1e-8), non-uniform 1-D field (1e-6), wrong ρ gives O(1); V3 flux-divergence stencil; V1 (4.8) = (4.7)/ρ pointwise; stratified shear: Dρ/Dt = 0, ∇·u = 0, ∂ρ/∂z = −0.4; V7 Galilean shift of U₀; V5 c(288.15 K) vs USSA-1976 340.29 m/s (1e-4), M(100) = 0.2939; V1 ball mean → point value O(r²); V2 D02 (Gauss on a box, polynomial flux), D03 (product rule, generic ρ, u) | order 2.001; ball 2.000 | analytic, symbolic, converged, benchmark | |
| C03 stream function; N13–N18; D04 | CORE | `velocity_from_streamfunction_2d(_sym)`, `…_axisym(_sym)`, `flux_between_streamlines`, `flux_along_path`, `streamfunction_preset`, `velocity_preset`, `STREAMFUNCTION_PRESETS`, `mass_flux_from_vector_potential`, `mass_flux_from_stream_functions`, `stream_function_pair`, `stream_surface_check`, `stream_tube_mass_flux`, `stream_surface_example` | V1 six presets × 200 points = closed-form velocities (1e-7), cylinder = ch03 field and ψ (1e-13), sign pinned (ψ = Uy → +U; ψ = kxy → (kx, −ky)); V2 ∇·u ≡ 0 (generic ψ, also ρ(x, y)), axisymmetric div ≡ 0; V1/V4 flux across 16 random gates and a bent gate = Δψ (1e-8); V3 stencil; V1 axisymmetric uniform stream and sphere equator 3U/2; V2 3-D: ∇·(∇χ × ∇ψ) ≡ 0, u·∇χ = u·∇ψ = 0 (generic χ, ψ), ∇×(χ∇ψ) = ∇χ × ∇ψ, parabolic pair = (1, 0, 2x); patch flux = (b − a)(d − c) (1e-9, two pairs); V2 D04 (χ = −z, flux across a curved gate = ψ(end) − ψ(start) symbolically) | order 2.000 | analytic, symbolic, converged | |
| C04 momentum for a moving CV (4.17); N19–N29, N31, N32, N82–N84, R03, R04; D05 | CORE | `momentum_budget`, `MomentumBudget`, `angular_momentum_budget/flux`, `body_force_from_potential`, `gravity_potential`, wake, bore, rocket, `jet_plate_force`, `cv_scenario`, `sprinkler_*` | V1 steady solid-body rotation through an off-centre fixed box and a translating, growing box (residual < 1e-9 / 1e-7 of the flux), traction route = τ route; uniform flow; Archimedes = ch01 `net_pressure_force_on_box` (1e-12); V3 unsteady Taylor–Green through a moving, growing box, residual vs dt; V4 material CV (4.13); V1 wake closed form (1e-10) and mass leak Δ√πb, contract 3.6523 N/m, 0.35449 m²/s; V7 Δ → 0, H-independence (1e-10), thrust sign; V2 wake closed form by sympy; bore: V1 U → √(gh) at O(Δh/h), p_o cancels, Bélanger form (1e-13), V2 mass + momentum ⇒ U; rocket: Tsiolkovsky form (1e-8), closed-form trajectory (1e-8), rtol sweep; V2 ODE; E1 scenarios all close (1e-8), bore with b ≠ U closes; sprinkler numeric disc flux = 2aρAU²cos α (1e-10), zero-swirl free spin; V4 angular-momentum budget of rigid rotation closes (1e-9); body force: −∇Φ (1e-7), zero work round a loop; V2 D05 in 1-D with arbitrary moving ends ((4.15) is an identity and not 0 by itself) | order 1.994 (dt); bore 0.998 | analytic, symbolic, converged, conserved | book (4.15) "= 0" typo confirmed by the D05 test |
| C05 Bernoulli (4.19); N30, N103–N105; D06 | CORE | `bernoulli_head`, `bernoulli_solve`, `stream_tube_element_balance_sym`, `pitot_speed(_from_heads)`, `stagnation_pressure`, `dynamic_pressure`, `torricelli_speed`, `orifice_mass_flow`, `tank_drain(_time)` | V1 cylinder surface C_p = 1 − 4 sin²θ from (4.19) (1e-12); V4 B constant along 4 traced streamlines (ptp/B < 1e-8) and a hydrostatic column; V1 solve ↔ head round trips (1e-14), raises; V2 dimensions (pint) of ½U² + gz + p/ρ and the example statement's (½)ρU² + gz does not add (analysis §9 item 7); pitot/heads/Torricelli/drain: closed-form level h(t) (1e-8), t_empty (1e-6), monotone in tank area; V5 C_c = 0.611; V2 D06 re-run independently (mass O(ds) = ρ(UA_s + AU_s); dropping the conical side force gives a wrong result; d(½U² + gz + p/ρ) = 0) | contract 28.868, 4.4294 | analytic, conserved, symbolic, benchmark | |
| C06 Cauchy (4.20)–(4.24); N02, N33–N37; D07 | CORE | `cauchy_terms`, `cauchy_residual`, `CauchyTerms`, `momentum_conservative_residual`, `conservative_to_advective_sym`, `stress_divergence`, `divergence_first_index_demo`, `closure_count` | V1 rigid rotation with τ = −pδ: both forms 0 (1e-9 of ρΩ²R); V1 **first index pinned**: τ₂₁ = x₁ at rest satisfies Cauchy with ∂τ_ij/∂x_i, the second-index divergence is 1; symmetric τ: both agree; V1 (4.23) numerically: flux form − advective form = u_j × continuity (1e-6) on a non-conserving field; V3 stencil; V2 (4.23) for generic ρ, u; V2 D07 (Gauss per component on a box pins the first index; the prose's ∂τ_ij/∂x_j differs; (4.20a)); V1 closure ledger 6 < 13, 4 < 5, 5 = 5, 7 = 7 | order 1.999 | analytic, symbolic, converged | book prose ∂τ_ij/∂x_j (analysis §9 item 8) discriminated |
| C07 Newtonian stress (4.25)–(4.37); N38–N51, R05, R06; D08, D09 ★★★, D10 | CORE | `newtonian_stress`, `viscous_stress`, `static_stress`, `total_stress`, `isotropic_fourth_order`, `linear_stress`, `mean_pressure`, `thermodynamic_pressure_from_stress`, `pressure_difference`, `bulk_viscosity`, `lam_from_bulk`, `stokes_assumption_holds`, `stress_on_plane`, `stress_lab_gradient`, `shear_stress_parallel_flow`, `deviatoric_part`, `cube_spin_acceleration` | V1 parallel shear τ₁₂ = μγ = ch01 `newton_shear_stress` (1e-14), rigid rotation σ = 0, rest = −pδ, traction −pn; V7 K isotropic under 50 random rotations (< 1e-12); V1 K:S = σ(μ_v = λ + ⅔μ) on 20 random G, only μ + γ acts, (4.31) ≡ (4.37), rotation covariance σ(CGCᵀ) = Cσ(G)Cᵀ; V1 (4.32)–(4.36) on 10 random states (1e-12); (4.35) and the incompressible guard; stress_on_plane contract (0.010, 0.0), rotation τ_s = 0, expansion σ_n = −p + 3μ_v r, = `tensors.normal_shear_stress`; cube spin slope −2 exactly, 0 when symmetric, 60 rad/s; V2 D08 (I = ρh⁵/6, α = 6(τ₁₂ − τ₂₁)/ρh²); V2 **D09 ★★★** re-running the design's cell (K from KroneckerDelta, symmetric S, only λ and μ + γ survive, σ symmetric for any γ, γ = μ gives 2μS + λS_mmδ, K symmetric in i, j ⇔ γ = μ, symbolic rotation about z leaves 20 random K_ijmn unchanged, parallel shear gives μγ̇); V2 D10 ((4.31) → (4.37), p̄ = p − μ_v∇·u, Stokes p = p̄) | slope −2.000 | analytic, symbolic | |
| C08 Navier–Stokes (4.38)–(4.41); N52–N55; D11–D13 | CORE | `ns_incompressible_terms`, `NSTerms`, `navier_stokes_residual`, `navier_stokes_sym`, `newtonian_stress_field`, `newtonian_viscous_stress_field`, `exact_solution(_fields)`, `exact_field`, `EXACT_SOLUTIONS`, `ns_terms_preset`, `viscous_force_forms`, `stokes_first_problem`, `plane_poiseuille` | V1 all 8 exact solutions (Couette–Poiseuille, plane Poiseuille with gravity, pipe Poiseuille, Stokes' first problem, Taylor–Green, Lamb–Oseen, ideal cylinder (Euler), solid body with gravity): (4.39b) residual < 1e-6 of the largest term, ∇·u = 0, and (4.38) form too; V2 five solutions satisfy (4.39b) and (4.38) symbolically (independent sympy fields), wrong viscous sign fails; V2 (4.38) ≡ (4.39a) for constant μ, μ_v, (4.39a) − (4.39b) = (μ_v + μ/3)∇(∇·u), variable μ(x) differs (D12), (4.38) = (4.24) + (4.37) (D11); V3 residual order ((4.39b) and nested (4.38)); V1 Taylor–Green = Wikipedia form (1e-12, form cross-check); Lamb–Oseen pressure = −∫ρu_θ²/r dr by quad (1e-9) and = ch03 Gaussian vortex with σ² = 4νt; ns_terms_preset Poiseuille contract (±100 N/m³), three viscous forms; gravity absorbed in p (N107); V1 (4.38) with μ_v and with μ(x) = sympy (1e-5); V2 D13 (ε–δ, generic u), stencil forms agree on a solenoidal field, compressible difference ±μ∇(∇·u), solid body 0 | order 1.999 / 1.999 | analytic, symbolic, converged | |
| C09 noninertial frame (4.42)–(4.45); N56–N65; D14–D18 | CORE | `rotating_basis`, `basis_rate(_exact)`, `inertial_velocity`, `frame_acceleration_terms`, `apparent_body_forces`, `coriolis_acceleration`, `coriolis_force`, `centrifugal_acceleration`, `centrifugal_potential`, `effective_gravity`, `earth_oblateness_diameter`, `coriolis_parameter`, `projectile_paths`, `coriolis_deflection`, `coriolis_projectile`, `rotating_pump_terms/equations`, `high_low_flow`, `core.curvilinear` | V3 de′/dt = Ω × e′; V1 12 random inertial paths seen from frames with Ω(t), U(t) (sympy-exact derivatives): the five terms reassemble the inertial acceleration (1e-10), "Coriolis without the 2" off by > 1e-3; (4.45) bracket = g − (terms) (1e-12); V2 **D15 ★★★** (design cell: x = X + R(t)x′ with generic X_i, x′_i, θ(t): D14 velocity, the five terms of (4.43), factor 1 leaves exactly Ω × ẋ′, triple product); V1 projectile ODE = closed form (7e-11 m), contract numbers, NH right / SH left, V4 speed conserved without centrifugal, D17 ratio → 1 at (Ωt)²/6 (order 2.000); V2 D17 series; V2 D18 centrifugal = −∇(−½Ω²R²), = +Ω²R e_R; V1 effective gravity = independent vector sum (1e-14), peak at 45.00°; V5 WGS-84 Ω, a, b, 2(a − b) = 42.77 km, Ω²a = 0.033916 m/s²; V2 D16 ∇²(U + Ω × x′) = 0, forces = −accelerations; V1 curvilinear operators vs Cartesian (cylindrical Laplacian, solid body div/curl/vector Laplacian/advective/strain, spherical potential flow, ∇²(1/r) = 0); V2 Example 4.5 rotation terms ρ(2Ω_zu_φ + Ω_z²R), −2ρΩ_zu_R and the curvature terms −u_φ²/R, u_Ru_φ/R | order 2.000 (basis), 2.000 (deflection) | analytic, symbolic, converged, conserved, benchmark | |
| C10 internal energy (4.57); N66–N81, R07, R08; D19–D23 | CORE | `energy_budget`, `EnergyBudget`, `stress_work_split`, `kinetic_energy_budget`, `internal_energy_terms`, `internal_energy_residual`, `total_energy_residual_sym`, `energy_identities_sym`, `energy_forms_sym`, `dissipation_rate`, `entropy_terms`, `entropy_production`, `couette_heating`, `couette_heating_transient` | V2 the implementer's chain (4.53)→(4.55), u·(4.24)→(4.56), (4.55)−(4.56)→(4.57), (4.54), (4.60)⇔(4.112) all ≡ 0; V2 **independent** 1-D re-derivation of D19–D22 (RTT + Gauss for an interval, E × continuity, u × Cauchy, (p/ρ) × continuity, force-work cancellation) and a 3-D box Gauss check of (4.49)–(4.51) (with dV, the book prints dA); V1 (4.53) residual ≡ 0 on the exact Couette-heating field, a doubled heating is not a solution; V1 steady Couette T satisfies (4.60) (2e-6 of μγ²), ΔT_max = μU²/8k, contract 2.0833e-4 K, 1.0 W/m²; V4 work in = heat out = ∫ρε dy (1e-10) with and without dp/dx; `energy_budget` box across the gap: shear work μU²/h = heat out, residual < 1e-9; V1 dissipation two routes on 1000 random G (1e-12), ≥ 0, μ < 0 gives ε < 0, νγ², 2νS:S, rotation invariance, plane 2 × 2 G = padded 3 × 3 (1e-13); V2 D23 (σ:S = 2μ dev² + μ_vS_mm², δ_ijδ_ij = 3, σ:R = 0); V1 stress-work split, kinetic-energy budget (Taylor–Green, Couette), (4.57) terms on Couette; V1/V2 entropy split and production sign, (4.62) forms equal; transient: long-time limit, interior PDE residual, independent quad series (1e-9) and t = 0 start (2.1e-7 / 6.4e-7 ΔT_max) pass | see text | analytic, symbolic, conserved | F1 fixed in loop 2 |
| C11 Bernoulli function (4.66)–(4.72), (4.76)–(4.78); N85–N96, N102; D24, D25 | CORE | `lamb_vector`, `lamb_identity_terms/sym`, `pressure_function`, `bernoulli_function`, `bernoulli_along_line`, `lamb_surface_check/terms`, `rankine_vortex_pressure`, `rankine_bernoulli`, `stagnation_enthalpy`, `stagnation_temperature`, `BERNOULLI_FORMS`, `which_bernoulli(_text)`, `bernoulli_scenario`, `BERNOULLI_SCENARIOS` | V2 Lamb identity ≡ 0 for generic u, wrong sign fails; V3 stencil; V1 solid body: u × ω = +2Ω²(x, y, 0), (u·∇)u = −Ω²(x, y); V1 pressure function: constant, isothermal, isentropic = quad of ρ(p) (1e-10), isentropic = C_p(T − T_o); V1/V4 Rankine: B = −1, 0 (contract), inside = closed form, varies inside, flat outside (1e-14), dp/dr = ρu²/r (1e-7), ∇B = u × ω and u·∇B = 0 (1e-6); cylinder B uniform at 1000 points (ptp/B < 1e-12); V1 T₀/T = 1 + (γ−1)M²/2 (1e-13), 304.98 K; V1 hypothesis table (six cases) and text form; six E6 scenarios consistent and flat where claimed; V2 D24 ((4.66) ≡ (4.69) with a barotropic P(p) and generic Φ), D25 (u·(u × ω) = ω·(u × ω) = 0) | order 2.000 | analytic, symbolic, converged, conserved | `lamb_vector` docstring sign (O2) |
| C12 unsteady Bernoulli (4.73)–(4.83); N97–N101; D26 | CORE | `unsteady_bernoulli_B`, `unsteady_bernoulli_pressure`, `gauge_absorbed_bracket`, `viscous_irrotational_residual`, `unsteady_streamline_bernoulli`, `u_tube_column`, `accelerating_sphere_pressure`, `accelerating_sphere_fields` | V1 accelerating sphere: the (4.74) bracket uniform at 30 points (1e-6), and the same (φ, p) satisfy Euler's equation by stencils (1e-4 of the terms — the independent route); V1 added-mass form cross-check: force = −½ρV dU/dt (1e-10), surface-pressure formula = field (1e-12); V1 gauge −1 → 0, printed + → 2B; U-tube residual 0 and L dU/dt + 2gh = 0; pipe p₁ − p₂ = ρL dU/dt; viscous force of three potential flows ≈ 0 (1e-5), a curved shear flow's is not; V2 D26 (Euler ≡ ∇[bracket] for a generic φ; the φ − ∫B sign absorbs B, + doubles it) | — | analytic, symbolic | book (4.74) sign slip confirmed |
| C13 Boussinesq (4.84)–(4.89); N106–N114, R10; D27, D28 | CORE | `perturbation_fields`, `buoyancy`, `boussinesq_momentum_terms`, `heat_equation_terms`, `temperature_equation_residual`, `boussinesq_density`, `boussinesq_validity`, `boussinesq_scenario`, `gaussian_blob_advection_diffusion`, `gaussian_blob_fields`, `blob_rise` | V1 hydrostatic column: p′ = ρ′ = 0 (1e-8); offsets recovered; buoyancy sign; N² = ch01 `brunt_vaisala_sq` (1e-9); linear EOS = ch01 `seawater_density_linear` (1e-14); rest state all 0, hydrostatic anomaly balanced (1e-9), unbalanced anomaly sinks; V1 blob residual (4.89) < 1e-6 of κ∇²T, (4.88) with ε = 0; V3 residual order; V4 ∫T′ dA = 2πσ₀² at 3 times (1e-8), centre moves at U; V1 validity contract numbers, c²/g for air = 11.8 km, four scenarios valid, deep atmosphere not; blob-rise ODE and terminal speed; V2 D27 (−∇p + ρg = −∇p′ + ρ′g with a hydrostatic base), D28 (perfect gas: ρC_vDT/Dt + p∇·u = ρC_pDT/Dt, α = 1/T) | order 2.010 | analytic, symbolic, converged, conserved | |
| C14 kinematic BC (4.90)–(4.93); N115–N129, R11; D29 | CORE | `kinematic_bc_residual`, `surface_normal_speed`, `relative_normal_velocity`, `interface_mass_flux`, `surface_preset`, `SURFACE_PRESETS`, `pillbox_limit`, `two_layer_conduction`, `two_fluid_couette`, `navier_slip_couette`, `linear_wave_surface/fields`, `wave_kinematic_residual`, `cap_pressure_force`, `cap_surface_tension_force`, `laplace_jump_from_balance`, `capillary_length`, `helmholtz_free_energy`, `spheroid_area`, `meniscus_height`, `meniscus_profile_x`, `meniscus_profile_ode` | V1 moving wall: Dη/Dt = 0, normal speed V; plane front: (u_rel)_n = U − V, flux ρ(U − V) (constant and field ρ); oblique, scaled η; V3 full (4.91) on the linear wave ∝ a² (residual/(aω) ≈ ka), linearised exactly 0 (deep and finite depth); V2 D29 (chain rule along a surface particle; wave residual starts at a²); V1/V3 pillbox residual ∝ l (slope 1); two-layer T and flux continuous; two-fluid Couette τ equal, u continuous; slip u = b du/dy; V1 cap pressure force dblquad = closed form (1e-10), x, y parts 0; V3 exact rim → small-ζ form at O(ζ); Laplace = ch01 (1e-12), 145.6 Pa; V2 the book's quarter-path integral → πab(1/R₁ + 1/R₂) → (1.5); V5 capillary length of water at 20 °C (IAPWS σ) vs 2.71 mm; Bo = (l/ℓ_c)²; V1 meniscus: h² = 2δ²(1 − sin θ), ODE (arc length) vs closed form (1.8e-7 δ, 4 angles); V2 closed form satisfies the first integral and has a negative slope (the book's dropped minus); V2 (∂f/∂v)_T = −p; spheroid area = quadrature (1e-10), minimum at aspect 1 | order 2.024 (wave), 0.995 (cap), 1.000 (pillbox) | analytic, symbolic, converged, benchmark | book curve-C sign and Ex. 4.7 slips handled |
| C15 dimensionless NS (4.99)–(4.119); N130–N156, R12; D30 | CORE | `Scales`, `nondimensional_ns_coefficients/sym`, `nondimensional_energy_coefficients`, `nondimensional_continuity_coefficient`, named numbers, `pressure/drag/lift_coefficient`, `reference_area`, `sphere_drag_coefficient`, `synthetic_sphere_drag_data`, `sphere_drag_pi_groups`, `froude_scaled_speed`, `model_prototype`, `ship_drag_extrapolation`, `prandtl_of` | V2 coefficients Ωl/U, 1, 1, gl/U², μ/ρUl; viscous and hydrostatic pressure scalings; Ec, Ec/Re, 1/(Pr Re); U²/c²; V2 **D30 independent**: a concrete scaled field substituted into (4.39b) and divided by ρU²/l equals the starred equation; V1 collapse of Stokes-layer and Poiseuille solutions (1e-12, 1e-14) and cylinder C_p for 3 (U, a); V2 every group dimensionless (pint); identities Ca = We/Re, Ri = 1/Fr′², Fr′(N) = Fr′(g′ = N²l); V1 coefficients and areas; V5 Morrison formula reproduced at 200 Re (1e-12), 24/Re within 1 % at Re = 0.1, synthetic data collapse; V1 Π groups Fρ/μ² = (F/ρU²D²)Re² exactly; V1 ship: wave-drag ratio (ρ_p/ρ_m)λ⁻³, Re ratio λ^−3/2 = 125, Re-matching; Scales round trip (1e-15), oscillating body St = 1; V5 Pr of air (0.70–0.73 at 250/280/300 K), water at 300 K, monatomic 2/3 | — | symbolic, analytic, benchmark, book-value | |

Every computable A item has ≥ 2 independent levels with at least one of V1/V2/V3/V5.

## Validation table — coded B/C items (NOTE, ≥ 1 level) and recaps
| Concept | fluidpy target | Evidence (test) | Label |
|---|---|---|---|
| N02 equations vs unknowns | `closure_count` | V1 `test_closure_count_V1_…` | analytic |
| N03–N06 (4.1)–(4.4) | `mass_budget(material=True)`, `material_mass` | V4, V2 D01 | conserved, symbolic |
| N07 localisation | `ball_integral` | V1 (order 2.000) | analytic |
| N08–N12 | continuity family, `mach_number`, `is_incompressible_regime`, `incompressible_speed_limit` | V1, V2, V5, V7 | analytic, symbolic, benchmark |
| N13–N18 (4.11)–(4.12), axisymmetric ψ | `core.streamfunction` | V1, V2, V3 | analytic, symbolic, converged |
| N19–N25 (4.13)–(4.16) | `momentum_budget` | V1, V3, V4, V2 D05 | analytic, converged, conserved, symbolic |
| N26–N28 body/surface forces, (4.18), drag sign | `body_force_from_potential`, `gravity_potential`, `cv_scenario`, `wake_drag_per_span` | V1 (−∇Φ, zero loop work, thrust < 0, surface = −F_D) | analytic |
| N29 Ex. 4.1 | wake family | V1, V2, V7, V6 | analytic, symbolic, book-value |
| N30 Ex. 4.2 | `stream_tube_element_balance_sym` | V2 | symbolic |
| N31 Ex. 4.3 | `bore_speed`, `bore_outlet_velocity`, `bore_pressure_force` | V1, V2, form cross-check, V6 | analytic, symbolic |
| N32 Ex. 4.4 | `rocket_*` | V1, V2 | analytic, symbolic |
| N33–N37 (4.20)–(4.23) | Cauchy family | V1, V2 | analytic, symbolic |
| N38 (4.25) | `cube_spin_acceleration` | V1, V2 D08 | analytic, symbolic |
| N40–N51 (4.26)–(4.37) | `core.constitutive` | V1, V2, V7 | analytic, symbolic |
| N52–N55 (4.38)–(4.41) | NS family | V1, V2, V3 | analytic, symbolic, converged |
| N57–N65 (4.42)–(4.45), Ex. 4.5 | `core.rotating`, `rotating_pump_*`, `high_low_flow`, `core.curvilinear` | V1, V2, V3, V4, V5 | as C09 |
| N66–N81 (4.46)–(4.63) | energy family | V1, V2, V4 | as C10 |
| N82–N84 (4.64)–(4.65), Ex. 4.6 | `angular_momentum_*`, `sprinkler_*` | V1, V4, V6 | analytic, conserved |
| N86–N96, N102 | Bernoulli-function family | V1, V2, V3 | as C11 |
| N97–N101 | unsteady family | V1, V2 | as C12 |
| N103–N105 pitot, orifice | pitot/orifice/drain | V1, V5 | analytic, benchmark |
| N106–N114 | Boussinesq family | V1, V2, V3, V4 | as C13 |
| N116–N129 | interfaces, surface tension, meniscus | V1, V2, V3, V5 | as C14 |
| N131 Fig. 4.21 | `sphere_drag_coefficient`, `synthetic_sphere_drag_data` | V5, V1 | benchmark |
| N133–N156 | `core.similarity`, `ship_drag_extrapolation`, `model_prototype`, `prandtl_of` | V1, V2, V5, V6 | as C15 |
| R01–R12 | recaps (ch01–ch03) | reused functions are called and asserted here where Ch. 4 relies on them (ch01 `newton_shear_stress`, `laplace_pressure_jump`, `seawater_density_linear`, `net_pressure_force_on_box`, `brunt_vaisala_sq`; ch03 `cylinder_flow`, `cylinder_streamfunction`, `gaussian_vortex`, `pathline`, `streamline`) | analytic |

Conceptual NOTE items without computable output (N01, N39, N56, N85, N92, N115, N124, N130, N132, N138, N140, N144)
are words only (curation §2); SKIP S01, S02 not coded.

## Derivations (curation §4b, design Part F) — every ★★ and ★★★ re-derived with sympy
| D | ★ | test | what is re-derived (intermediate lines checked) | finding |
|---|---|---|---|---|
| D01 | ★★ | `test_mass_cv_V2_derivation` | 1-D: generic ρ, u from continuity, arbitrary x₀(t), x₁(t); storage + [ρ(u − b)] = 0 (4.5); Leibniz = (4.3); material (4.2); the u·n slip rejected | ✓ |
| D02 | ★★ | `test_continuity_V2_derivation` | Gauss on a box for a polynomial flux; localisation numerically (`test_localisation_V1_…`) | ✓ |
| D03 | ★ | `test_continuity_V2_derivation` | ∇·(ρu) = u·∇ρ + ρ∇·u ⇒ (4.8) | ✓ |
| D04 | ★★ | `test_streamfunction_V2_derivation` | χ = −z ⇒ ρu = ∂ψ/∂y, ρv = −∂ψ/∂x; flux across a curved gate = Δψ | ✓ |
| D05 | ★★ | `test_momentum_cv_V2_derivation` | 1-D with a stress built to satisfy Cauchy: (4.17) closes; (4.15) is an identity, not 0 | ✓ (book "= 0" typo) |
| D06 | ★★ | `test_bernoulli_V2_derivation` | independent Taylor re-run; conical side force needed; d(½U² + gz + p/ρ) | ✓ |
| D07 | ★★ | `test_cauchy_V2_derivation`, `test_cauchy_V2_conservative_to_advective_symbolic` | Gauss per component pins the first index; (4.20a); (4.23) = u_j × continuity | ✓ |
| D08 | ★★ | `test_stress_symmetry_V2_derivation` | I = ρh⁵/6, torque, α = 6(τ₁₂ − τ₂₁)/ρh² | ✓ |
| D09 | ★★★ | `test_newtonian_law_V2_derivation` | the design's check cell step by step (steps 4–12), plus the book's K-symmetry argument | ✓ |
| D10 | ★ | `test_bulk_viscosity_V2_derivation` | (4.31) → (4.37), (4.33)–(4.34), Stokes | ✓ |
| D11 | ★★ | `test_navier_stokes_V2_forms_4_38_4_39a_4_39b_agree` | (4.24) + (4.37) = (4.38) for generic fields | ✓ |
| D12 | ★★ | same | (4.38) ≡ (4.39a) (constant μ, μ_v), − (4.39b) = (μ_v + μ/3)∇(∇·u), variable μ differs | ✓ |
| D13 | ★★ | `test_viscous_force_V2_derivation`, `test_viscous_force_V2_three_forms_…` | ε–δ for generic u: 2∂S_ij/∂x_i = ∇²u + ∇(∇·u), −∇×ω = ∇²u − ∇(∇·u) | ✓ |
| D14 | ★★ | `test_rotating_frame_V2_derivation` | Rᵀẋ = RᵀẊ + ẋ′ + Ω × x′ | ✓ |
| D15 | ★★★ | `test_rotating_frame_V2_derivation` | the design's cell: five terms of (4.43) for generic X_i(t), x′_i(t), θ(t); factor 1 leaves Ω × ẋ′; triple product | ✓ |
| D16 | ★★ | `test_rotating_ns_V2_derivation` | ∇²(U + Ω × x′) = 0; forces = −accelerations | ✓ |
| D17 | ★ | `test_coriolis_V2_derivation_small_time_deflection` | series ut sin Ωt → Ωut², angle Ωt | ✓ |
| D18 | ★ | `test_centrifugal_V2_potential_and_effective_gravity` | −Ω × (Ω × x) = −∇(−½Ω²R²) | ✓ |
| D19 | ★★ | `test_energy_V2_derivation_one_dimensional_chain`, `test_energy_budget_V2_derivation_box_gauss` | interval RTT + Gauss ⇒ (4.53); 3-D box Gauss of the flux terms (dV, not the printed dA) | ✓ (book (4.51) dA typo) |
| D20 | ★★ | same (+ `test_energy_V2_identity_chain`) | (4.53) − E × continuity = (4.55) | ✓ |
| D21 | ★★ | same | u × Cauchy (4.24) = (4.56) | ✓ |
| D22 | ★★ | same | (4.55) − (4.56) = ρ(4.57) + (p/ρ) × continuity; force-work cancels | ✓ |
| D23 | ★★ | `test_dissipation_V2_derivation` | σ:S = 2μ dev:dev + μ_vS_mm², δ_ijδ_ij = 3, σ:R = 0 | ✓ |
| D24 | ★★ | `test_bernoulli_function_V2_derivation`, `test_lamb_identity_V2_symbolic_and_sign` | barotropic P(p) (FTC + chain rule), (4.66) ≡ (4.69), Lamb identity, sign | ✓ |
| D25 | ★ | `test_bernoulli_function_V2_derivation` | u·(u × ω) = ω·(u × ω) = 0 | ✓ |
| D26 | ★★ | `test_unsteady_bernoulli_V2_derivation` | Euler ≡ ∇[bracket]; φ − ∫B dt absorbs B, + doubles it | ✓ (book (4.74) sign) |
| D27 | ★★ | `test_boussinesq_V2_derivation` | −∇p + ρg = −∇p′ + ρ′g with a hydrostatic base | ✓ |
| D28 | ★★ | same | ρC_vDT/Dt + p∇·u = ρ(C_v + R)DT/Dt, α = 1/T | ✓ |
| D29 | ★★ | `test_kinematic_bc_V2_derivation` | chain rule along a surface particle; wave residual O(a²), linearised exact | ✓ |
| D30 | ★★ | `test_similarity_V2_derivation`, `test_similarity_V2_ns_coefficients` | independent concrete substitution ÷ ρU²/l = starred equation; the implementer's coefficient routine | ✓ |

Part F's intermediate lines were checked where they carry a formula (D01, D05, D07, D09 steps 4–12, D15 steps 1–12,
D19–D23, D26, D29, D30); no wrong intermediate line was found.

## Functions used by the notebook and explainers (design Part C) — test name each
`test_contract_V1_every_part_c_name_exists_and_is_reexported` checks **all 203 Part C names** (195 from C.1–C.11 parsed from the
design plus the 8 C.6 operators reached as `ch04.CU.*`) exist on `ch04`, that they are the core modules' objects (re-exports, not copies),
and the contracted leading parameters. `test_contract_V1_every_part_c_name_is_exercised_in_this_file` asserts that
every one is called as `ch04.<name>` (or `CU.<name>`) somewhere in the suite with a physical assertion (the result
classes and constants in `test_result_types_V1_named_fields_and_constants`). `test_scalar_callable_V1_…` checks that the
22 functions the explainers' parity rows call return plain floats. C.12 drawing helpers: `test_scripts_V1_drawing_helpers_run`.

Contract numbers of design Part C, each asserted: interval budgets, `material_mass` = 1.0, wake 3.6523 N/m and
0.35449 m²/s, sprinkler 0.8660 N m, pitot 28.868 m/s, Torricelli 4.4294 m/s, τ₁₂ = 0.010 Pa, `stress_on_plane`
(0.010, 0.0), cube spin 60, Poiseuille ±100 N/m³, projectile 36 000 / 9450.6 / 9342.4 m / 0.26252 rad, equatorial
centrifugal 0.033916 m/s², ε = 1.0 W/kg, Couette 2.0833e-4 K and 1.0 W/m², Rankine B = −1, 0, `which_bernoulli` lists,
304.98 K, Boussinesq 2.0e-3 and 0.01962 m/s², Laplace 145.6 Pa, Re = 1e4, Morrison 0.3926, Froude 2.0 m/s, Ri = 98.1,
Ro = 0.06857, `bore_speed(1.0, 1.1)` = 3.3661 m/s and `capillary_length(0.0728, 998.0)` = 2.7269 mm (the last two
corrected in the design in loop 2, O3).

## Convergence studies (scheme | steps | observed order (pairwise) | design order)
| scheme | steps | observed | design |
|---|---|---|---|
| continuity ∇·(ρu) stencil | h = 0.1 … 0.0125 | 2.001 (2.002, 2.001, 2.000) | 2 |
| ψ → (u, v) central differences | h = 0.08 … 0.01 | 2.000 (2.001, 2.000, 2.000) | 2 |
| mass budget (4.5), moving growing sphere | dt = 0.2 … 0.025 | 2.023 (2.055, 2.014, 2.003) | 2 |
| momentum budget (4.17), moving growing box, unsteady flow | dt = 0.4 … 0.05 | 1.994 (1.985, 1.996, 1.999) | 2 |
| Cauchy (4.24) residual | h = 0.1 … 0.0125 | 1.999 (1.997, 1.999, 2.000) | 2 |
| NS (4.39b) residual, Taylor–Green | h = 0.1 … 0.0125 | 1.999 (1.997, 1.999, 2.000) | 2 |
| NS (4.38) residual (nested stress stencil) | h = 0.1 … 0.0125 | 1.999 (1.997, 1.999, 2.000) | 2 |
| Lamb identity (4.68) residual | h = 0.1 … 0.0125 | 2.000 (2.000, 2.000, 2.000) | 2 |
| de′/dt = Ω × e′ | h = 0.1 … 0.0125 | 2.000 (1.999, 2.000, 2.000) | 2 |
| (4.89) residual, advected Gaussian | h = 0.02 … 0.0025 | 2.010 (2.023, 2.006, 2.001) | 2 |
| ball mean − point value (localisation) | r = 0.4 … 0.05 | 2.000 | 2 |
| full (4.91) on a linear wave | a = 0.08 … 0.01 | 2.024 (2.042, 2.021, 2.010) | 2 (in a) |
| cap (4.98) exact / small-ζ − 1 | ζ = 1e-5 … 1.25e-6 | 0.995 (0.991, 0.996, 0.998) | 1 |
| pillbox side + volume terms | l = 0.1 … 0.0125 | 1.000 | 1 |
| bore U/√(gh) − 1 | Δh = 0.2 … 0.025 | 0.998 (0.997, 0.999, 0.999) | 1 |
| Coriolis ut sin Ωt vs Ωut² | t = 0.4 … 0.05 | 2.000 (1.999, 2.000, 2.000) | 2 |
| cube spin acceleration vs h | h = 0.1 … 0.0125 | −2.000 | −2 (divergence) |
| rocket rtol sweep | 1e-4 … 1e-10 | monotone, final 1e-7 relative | — |
| Couette transient at t = 0 vs nterms (loop 2) | 25 … 200 | pairwise 1.8, 3.0, 3.0 (4.3e-5 → 2.1e-7 ΔT_max) | ≥ 2 (b_n ∝ n⁻³) |

## Conservation / invariant residuals
Fixed-box mass: storage + outflux < 1e-7 of |dM/dt|; material volume: outflux 0, storage < 1e-8 of the mass; material
interval mass ptp < 1e-12 over t ∈ [0, 5]; steady rotating flow through fixed/moving boxes: momentum residual < 1e-9 /
1e-7 of the flux; angular-momentum budget < 1e-9; Couette heating: work in = heat out = ∫ρε dy (1e-10) and the
`energy_budget` box residual < 1e-9 of μU²/h; Bernoulli constant along traced cylinder streamlines (ptp/B < 1e-8);
cylinder B uniform at 1000 points (1e-12); Rankine B flat on circles (1e-12) and outside the core (1e-14); projectile
speed constant without the centrifugal term (1e-9 of u₀); ∫T′ dA of the advected blob = 2πσ₀² (1e-8); E1 scenario
budgets (1e-8); ψ-flux across any gate = Δψ (1e-8).

## Benchmarks used (value | our value | source + URL | date verified)
Only published numbers are V5; Wikipedia forms (Bélanger, Tsiolkovsky, added mass, Taylor–Green) are V1 "form
cross-checks (not V5)".
| value | ours | source | verified |
|---|---|---|---|
| Morrison C_D(Re) formula, 0.1 ≤ Re ≤ 1e6; 24/Re creeping limit | identical at 200 Re (1e-12); 240.22 vs 240 at Re = 0.1 (0.09 %) | Morrison (2016), Michigan Tech, Eq. (1), https://pages.mtu.edu/~fmorriso/DataCorrelationForSphereDrag2016.pdf (PDF text extracted) | 2026-09-23 |
| WGS-84 a = 6 378 137 m, 1/f = 298.257223563, b = 6 356 752.314245 m, ω = 7.292115e-5 rad/s | constants identical; 2(a − b) = 42.769 km; Ω²a = 0.0339157 m/s² | Wikipedia "World Geodetic System" | 2026-09-23 |
| USSA-1976 sea-level c = 340.29 m/s | √(γRT) with ch01's R_air: agrees to < 1e-4 | NASA-TM-X-74335 (reference/ch01) | 2026-09-12 |
| C_c = 0.611 (sharp orifice) | used as `Cc`; jet area and ṁ ratio exact | Wikipedia "Vena contracta" | 2026-09-23 |
| capillary length of water at 20 °C = 2.71 mm | 2.7254 mm (IAPWS σ = 72.74 mN/m, ρ = 998.2, g₀) — 0.57 %, inside 1 % + the 3-digit rounding | Wikipedia "Capillary length" | 2026-09-23 |
| Pr of air 0.70–0.73 (250–1000 K) | 0.707–0.710 at 250, 280, 300 K | Wikipedia "Prandtl number" | 2026-09-23 |
| Pr of water 5.9 at 300 K | 5.811 (1.5 % below; inside 1 % + the 2-digit rounding band) | Wikipedia "Prandtl number" | 2026-09-23 |
| Pr = 2/3 monatomic; Eucken 4γ/(9γ − 5) | 0.666… exactly | Wikipedia "Prandtl number" | 2026-09-23 |
| Bélanger h₂/h₁ = (√(1 + 8Fr₁²) − 1)/2 (form) | bore speed satisfies it (1e-13) | Wikipedia "Hydraulic jump" | 2026-09-23 |
| Tsiolkovsky Δv = v_e ln(m₀/m_f) (form) | 1e-8 | Wikipedia "Tsiolkovsky rocket equation" | 2026-09-23 |
| added mass of a sphere ½ρV (form) | force = −½ρV dU/dt (1e-10) | Wikipedia "Added mass" | 2026-09-23 |
| Taylor–Green u, v, p (form) | identical (1e-12) | Wikipedia "Taylor–Green vortex" | 2026-09-23 |
Morrison's remark "plateaus at C_D ≈ 0.14 at the highest Re" is a verbal description of the curve; the formula itself
gives 0.1296 at Re = 1e6 (we reproduce the formula exactly), so it is not used as a numeric benchmark.

## Numbers from the text (book vs ours)  [private values redacted to relative errors]
| book item | our route | relative error | comment |
|---|---|---|---|
| Ex. 4.8 model speed, Re_m, Re_p, model friction and wave drag, prototype friction | `ship_drag_extrapolation(rho_p = 1000)` | 0 (to 2e-16) | |
| Ex. 4.8 prototype wave drag, uncorrected estimate | same | +5.6e-4, +5.3e-4 | book rounding |
| Ex. 4.8 prototype total | same | +1.1e-3 | the book adds its rounded wave drag (analysis §9 item 12) |
| Ex. 4.1, 4.3, 4.6, 4.7 closed forms; cap forces; projectile deflection (parsed from the private JSON) | wake, bore, sprinkler, meniscus, cap, projectile functions | ≤ 1e-7 | formulas identical |
| Ex. 4.5 radial rotation terms | `rotating_pump_terms` | 0 (sympy) | |
| capillary scale at 288 K | `capillary_length` | −1.0e-3 | |
| Earth oblateness | `earth_oblateness_diameter` (WGS-84) | +1.8 % | book truncates 42.77 km (pre-registered 2 % in analysis §6) |
| §4.2/§4.11 incompressibility speed and the bound on M² | `mach_number`, `incompressible_speed_limit`, `compressibility_parameter` | both bounds hold; the speed at M = 0.3 is within 2.1 % of the book's round number | asserted as the book's inequalities |
| kinetic-theory (monatomic) Prandtl number | `eucken_prandtl(5/3)` | −5.0e-3 | book rounds 2/3 |
| c²/g for air (§4.9) | 11.8 km | +18 %; equal at the book's 1-significant-figure precision ✓ | order of magnitude |
| sharp-orifice contraction (§4.9) | Wikipedia 0.611 | −1.5 % | two sources; pre-registered 2 % |
| Pr of air and of water at 20 °C (§4.11) | `prandtl_of` at 293.15 K: 0.709, 6.947 | −1.5 %, −2.2 % | **not asserted** (property-table difference: our μ(T), k, C_p come from ch01's Sutherland/Vogel + Incropera rows; the V5 checks against Wikipedia pass) — Open item O5 |

## Figures reproduced with our code (outputs/ch04/verify/, local)
| figure | file | one-sentence visual verdict |
|---|---|---|
| Fig. 4.8 | `fig4_8_projectile_verify.png` | The inertial path is the straight segment 0 → 12 m on y = 0, and the rotating-frame path curls clockwise (to the right) from the origin through (3.7, −5) to (−2.7, −11.7) m, with the 61 ODE markers on the closed-form curve (max gap 7e-11 m). |
| Fig. 4.9 | `fig4_9_effective_gravity_verify.png` | The deflection of g_e from the radial is 0 at the equator and pole and peaks at 45.00° with 5.96 arcmin, while |g_e| rises monotonically from 9.766 to 9.800 m/s². |
| Fig. 4.21 | `fig4_21_sphere_drag_verify.png` | The Morrison curve follows 24/Re below Re ≈ 1, flattens near 0.4 for 10³–10⁵, drops to ≈ 0.09 near Re ≈ 4e5 and rises to 0.13 at 1e6; synthetic air, water and glycerine points of many sizes collapse on it. |
| Ex. 4.7 | `ex4_7_meniscus_verify.png` | Meniscus profiles for θ = 10°, 30°, 60°, 85° start at h/δ = 1.28, 1.00, 0.52, 0.09 and decay monotonically, with the closed-form x(ζ) dots on the ODE curves (max 1.8e-7 δ). |
| (4.60) Couette heating | `couette_heating_verify.png` | The steady profile is the parabola peaking at y = h/2 (ΔT_max = μU²/8k); loop 2: the t = 0 curve lies on T₀ in both panels and the transients grow smoothly toward the steady profile (peak at y = h/2 without, shifted toward the fixed wall with dp/dx = −2000 Pa/m); in loop 1 it wiggled (3e-3) and was noise of ±2 ΔT_max (F1). |
| (4.71) Rankine | `rankine_bernoulli_verify.png` | B rises from −1 on the axis as r² to 0 at r = σ and is flat outside, while p/ρ and ½u_θ² vary everywhere (the kink of u_θ at σ). |
| convergence | `convergence_verify.png` | Eleven straight log–log lines, nine of slope 2.00–2.02 and two of slope 1.00, none flattening at round-off. |

## Deviations & justifications
One `# DEVIATION` in the code: `meniscus_profile_ode` integrates the curvature balance in arc length instead of the
book's dζ/dx separation (the wall slope −cot θ is infinite as θ → 0); verified against the closed form to 1.8e-7 δ
(`test_meniscus_V1_…`). Book typos implemented as corrected physics and tested: (4.15)'s trailing "= 0" (D05 test shows
the integral is not 0 by itself); (4.51) dA → dV (box Gauss test); (4.74) gauge sign (D26 test, `gauge_absorbed_bracket`);
Cauchy prose ∂τ_ij/∂x_j (first-index discrimination); curve C sign (quarter-path integral with +); Ex. 4.7's dropped
minus (negative slope asserted); Ex. 4.2's (½)ρU² (pint). Test-side choices recorded honestly: stencil steps for the
exact-solution residuals are 2e-4 × the solution's length scale (the ns_terms_preset default); the steady Couette
residual uses h = 5e-5 m in a 1 mm gap (T is quadratic, so the stencil is exact and h only controls round-off).

## Open items
- **Resolved in loop 2:** F1 (closed-form b_n), O2 (`lamb_vector` docstring now +2Ω²(x, y, 0)), O3 (design contract
  numbers 3.3661 m/s, 2.7269 mm), O4 (FTCS claim removed), O6 (quad epsrel 1e-11; no IntegrationWarning in the suite).
- **O5 (knowledge, still standing — a note, not a defect):** the book's Prandtl numbers of air and of water at 20 °C
  differ from `prandtl_of` at 293.15 K by −1.5 % and −2.2 %: a property-table difference (our μ(T), k, C_p come from
  ch01's Sutherland/Vogel fits and Incropera rows; they agree with Wikipedia's air range and 300 K water value). If the
  notebook quotes the book's values next to ours it should say so. The V6 tests do not assert these two values.
- Nothing is labelled `qualitative` or `unverified`.

## Verdict: PASS
Loop 2 of 3. All 8 scripts ran; `tests/test_ch04.py` 147/147 pass (144 + 3 skipped without the private file); the full
suite 474/474; every computable CORE row has ≥ 2 independent levels (one of V1/V2/V3/V5), every coded NOTE row ≥ 1,
all 203 Part C names are exercised, every ★★/★★★ derivation is re-derived symbolically, and 22 wrong variants are each
caught. F1 is fixed and confirmed by an independent quadrature of the series coefficients; the only open item (O5) is
a documented property-table difference.

## Review follow-up (reports/ch04_review.md) — 2026-09-23
The code review found four paths the suite did not reach; the implementer fixed them and the verifier added tests (no
`fluidpy/` edits, no tolerances loosened). New tests, all V1:

| test | what it pins | numbers |
|---|---|---|
| `test_jet_on_plate_V1_split_from_along_plate_momentum` | E1 "jet": Q₁,₂ = Q(1 ± cos θ)/2 from the along-plate momentum balance, derived by hand in the test (in ρV²A cos θ, out ρV²(A₁ − A₂), no plate shear) and compared with the function's `residual_momentum_along` and per-face `momentum_flux_along`; θ = 5°…170° sweep; normal force ρV²A sin θ | θ = 45°, ρ = 1000, V = 10, A = 1e-3: sheets 8.5355 / 1.4645 kg/s, force 70.711 N; along-plate residual < 1e-13 ρV²A for all θ; the old 50/50 split leaves −ρV²A cos θ (≠ 0 except at 90°); other scenarios report 0 |
| `test_mean_pressure_V1_plane_tensor_needs_tau33` | (4.33)–(4.34) for 2 × 2 plane-flow stress: ValueError without `tau33` (both `mean_pressure` and `thermodynamic_pressure_from_stress`); with τ₃₃ = −p + (μ_v − ⅔μ)∇·u = the 3 × 3 result; 3 × 3 input plus `tau33` raises | G = diag(1, 0), p = 100, μ = 1, μ_v = 0.5: p − p̄ = 0.5 = μ_v∇·u (1e-13); p recovered (1e-14) |
| `test_rotating_frame_V1_translation_vector_and_plane_inputs` | a nonzero scalar dU/dt raises (both frame functions), scalar 0 accepted; 2-component inputs are padded to 3 components and equal the 3-vector call with z = 0 | `apparent_body_forces([1, 0], [0.5, 0.2], 0.3)`: Coriolis (0, −0.6, 0), centrifugal (0.045, 0.018, 0) |
| `test_boussinesq_validity_V1_mach_rule_agrees_with_incompressible_regime` | the Mach flag of `boussinesq_validity` uses the §4.2 rule M < 0.3, the same as `is_incompressible_regime` | M = 0.1, 0.2, 0.29 not flagged, 0.31, 0.5 flagged, identical to `is_incompressible_regime` at the same c; M = 0.2 "Boussinesq valid" |

Discrimination (the pre-fix behaviour planted back into a scratch copy, one at a time; scratchpad `disc_review.py`):
| wrong variant | tests that fail |
|---|---|
| jet: the old 50/50 split | 1 |
| Boussinesq Mach flagged at the generic threshold 0.1 | 1 |
| scalar dU/dt read silently as a z-component | 1 |
| `mean_pressure` as −tr τ/d for a 2 × 2 tensor (the old behaviour) | 1 |

Runs: `tests/test_ch04.py` → **151 passed** (114 s); full suite `pytest -q` → **478 passed** (209 s).
Collected items per evidence tag: **V1 83 · V2 41 · V3 9 · V4 6 · V5 6 · V6 3 · V7 3** (= 151; 144 test functions,
the exact-solution test parametrised over 8). 26 wrong variants caught in total (17 + 5 + 4).

### Verdict after review follow-up: PASS
Open items unchanged: O5 (property-table note on the book's Prandtl numbers) only.

