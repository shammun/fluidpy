# Chapter 6 verification — Ideal Flow                     2026-09-23, commit 02a3d52 + working tree (loop 3, review follow-up)

Verifier: math-verifier. Suite: `tests/test_ch06.py` (131 test functions = 137 collected items; two stagnation-point tests are
parametrised over four source strengths each), figures/metrics: `tests/ch06_verify_figures.py` → `outputs/ch06/verify/`
(git-ignored), cited data: `reference/ch06/` (`make_refs.py`, `benchmarks.json`, `SOURCES.md`). Nothing in `fluidpy/`
or `scripts/` was edited by the verifier.

## Environment
Python 3.11.5 · numpy 2.4.6 · scipy 1.17.1 · sympy 1.14.0 · pint 0.25.3 · matplotlib 3.11.2 (Windows, `.venv`).

Runs (loop 3, final): `pytest tests/test_ch06.py -q -p no:cacheprovider` → **136 passed in 108 s** (with the private
book file; without it the four V6 tests skip). Full suite `pytest -q` → **749 passed in 542 s** (613 earlier + 136 ch06).
Loop 2: 133 passed, full suite 746 passed.
All 12 `scripts/ch06_*.py` exit 0 (asserted by `test_scripts_V1_every_ch06_script_runs`). All five modules import
(`python -c "import fluidpy.ch06_ideal_flow"`). `tools/check_public.py` → OK. Loop 1: 123 of 128 passed, full suite
736 passed / 5 failed (F1 × 4, F2).

## Loop history
| loop | result | what changed |
|---|---|---|
| 1 | FAIL: 123/128 (full suite 736 passed, 5 failed) | F1 (stagnation search misses weak-source half-bodies), F2 (three Part C 6.1 drawing helpers missing); O1–O6 opened; 28/28 planted wrong variants caught. |
| 2 | **PASS: 133/133** (full suite 746 passed) | Implementer: F1 — strength-scaled Newton seeds (rings at ℓ·(¼, ½, 1, 2, 4) round every source, vortex and doublet; ℓ = m/2πU, Γ/2πU or √(d/2πU)), on-axis roots snapped to +0 imaginary part, a zero-residual start accepted as a root, a speed scale for stream-free flows; F2 — `draw_body`, `pressure_arrows`, `separated_cp_band` added, `flow_net` given the Part C signature; O4/O5 — docstrings; design Part G errata (G1–G5) for O1–O3. Verifier: 5 new robustness items (below); the drawing test now checks the band follows (6.35) ahead of separation and the arrows point along −p n; the 28 wrong variants re-planted into the loop-2 code: **28/28 caught**. |

Collected items per evidence tag (tag of the `def` line), after the odd-N addendum: **V1 86 · V2 27 · V3 10 · V4 6 · V5 2 ·
V6 4 · V7 2** (= 137); loop 3 was V1 85 (= 136); loop 2 was V1 82 · V2 27 · V3 10 · V4 6 · V5 2 · V6 4 · V7 2 (= 133).

## Loop 3 — review follow-up (reports/ch06_review.md), verifier's test updates
| item | change in `tests/test_ch06.py` / `reference/ch06/` | result |
|---|---|---|
| M1 (6.82) is not a slip | removed from the slip list; `test_axisym_velocity_V1_spherical_components_from_psi_and_phi` now asserts with sympy that the book's printed form (1/r)∂(r²u_r)/∂r + (1/sin θ)∂(u_θ sin θ)/∂θ vanishes for the sphere flow (6.90), that it equals r × the App. B divergence for a generic field (`core.curvilinear`), and numerically that it equals r × `spherical_continuity_residual` on a non-solenoidal test field (1e-8). The loop-1 "hybrid" assertion (a form the book never prints) is gone. | pass |
| S2 flux check is a tautology | the old test is relabelled `test_example_6_2_V1_flux_identity_and_maximum_principle` (V1); new `test_example_6_2_V4_discrete_circulation_vanishes`: `cell_circulation` recomputed independently as −(ψ_E + ψ_W + ψ_N + ψ_S − 4ψ_P) (1e-15), NaN off the unknowns; converged max 3.3e-10 ≤ 4 tol (Γ_cell = −4 × the average-rule defect, tol = 1e-10); after 1, 5, 20, 50 sweeps it falls monotonically, ≈ 0.0999 after 5 sweeps; the direct solve at refine 4 < 1e-12 | pass |
| S3 geometry/BCs public | new `test_example_6_2_V1_geometry_and_boundary_values` (Q = 1, no book values): 24 unknowns, solid = {x > 5, y < 2}, inlet ψ = Qy/5, outlet ψ = Q(y − 2)/3 with uniform steps Q/3, ψ = Q on top, ψ = 0 on the lower wall, the step face and the outlet's lower wall; 121 unknowns at Δ = ½. Planted non-uniform outlet (ψ ∝ ((y − 2)/3)²) is caught by this public test. | pass; WV caught |
| S6 `sphere_potential_vector(a=)` | the test uses `a=a`; both d_vec and a → ValueError, neither → ValueError. The implementation keeps a positional scalar as a **deprecated** radius (DeprecationWarning) rather than rejecting it as the review proposed; the test pins the warning and the correct value. | pass (note) |
| S1 tilted ellipse | `test_blasius_V4_three_routes_agree`: F_perp = ρUΓ_cw (1e-10), F_par ≈ 0 (1e-9), stream_angle 15°, c₀ = Ue^{−iα} (1e-12), D − iL = −iρUΓe^{−iα} (1e-9); also α = 0°, 30°, −20° explicitly (α = 0 stays 0) and the cylinder (F_perp = L) | pass |
| S5 net sink | `test_superposition_state_V1_explainer_terms`: stream + sink 2π → "open body (net sink): extends upstream", net −2π, stagnation at x = +m/2πU = +1 m; planted "no body" status caught | pass; WV caught |
| S4 ρ defaults | audited every force call: all loop-1/2 assertions pass ρ explicitly or already expected 1.2; new `test_lift_V1_shared_air_density_default` pins ρ = 1.2 kg/m³ in `lift_per_span`, `surface_pressure_force`, `blasius_force`, `cv_force_on_body`, `cylinder_circulation_state`, `laurent_contributions`, `blasius_state` (all 24 N/m), `force_on_held_singularity` (−7.2 N/m) and `cylinder_surface_pressure` (60 Pa) | pass |
| S7 primary citations | `reference/ch06/SOURCES.md` and `make_refs.py`: added Lamb, *Hydrodynamics*, 6th ed. (1932), Ch. V §92 (section number as cited by the review; the text was not fetched — only the edition confirmed) and Rayleigh, Phil. Mag. (6) 34, 94–98 (1917), doi:10.1080/14786440808635681 (confirmed by search); Wikipedia kept as the fetched secondary | done |
| S8 Exercise 6.42 | the V6 test now asserts the one correct cubic (z/a)²(1 − z/a) = +Q/4πUa² at the nose and −Q/4πUa² at the tail (1e-9), and only records that the book prints a ± form | pass |

No tolerance was relaxed in loop 3; all new bounds are round-off or derived (4 tol for Γ_cell). Two wrong variants planted
for S3/S5 are caught (2/2); the 28 earlier variants are unaffected (code paths unchanged).

## Loop-2 checks (on our own terms)
- **F1 re-check and robustness** (`test_stagnation_points_V1_default_search_finds_the_half_body_nose[m]`, m = 0.5, 1,
  2, 3; `test_stagnation_points_V1_extreme_source_strengths[m]`, m = 1e-4, 0.01, 20, 200 m²/s; and
  `test_stagnation_points_V1_robustness_cases`): half-body x_S = −m/2πU to 1e-12 relative from 1.6e-5 m to 32 m, with
  ψ_S = m/2 (the +0 snap reads ψ above the source's cut); two equal sources without a stream → the midpoint (1e-12);
  source + sink without a stream → none (correct: dw/dz ≠ 0 at every finite point); cylinder Γ_cw = 6π (U = a = 1) →
  one free point −i(6 + √20)/4 (1e-12), both roots with `keep_inside` and r₊r₋ = a²; Γ_cw = 0, 2, 20 equal the closed
  form (1e-12); the merged case Γ_cw = 4π is a double root, found to 1e-8 (asserted < 1e-6: a double root's position
  is determined only to √(residual tolerance) — a property of the problem, not a loosened tolerance); corner n = 4
  (triple root at the tip) → 0 (1e-12); a source outside the default box at 7 + 5i, a vortex in a stream, a tiny
  cylinder a = 0.05, a Rankine oval with m = 0.05 through `superposition_state`, and U = 0.01 (x_S = −15.9 m) all to
  1e-12.
- **F2 re-check**: the three helpers exist; `separated_cp_band` returns lo ≤ hi, contains the ideal 1 − 4 sin²β ahead of
  separation (above its documented −1.05 floor) and a negative wake plateau behind; `pressure_arrows` draws exactly
  −p n (quiver components compared); `draw_body` fills one patch; `scripts/ch06_cylinder_cp.py` uses the band with the
  label "qualitative sketch band, not data".
- **Docstring audit (O4, O5)**: no function or class in the five modules still says "Validation (planned)"; none lost its
  docstring; all modules import. A script parsed the 42 scripted "Validation: V… — tests/test_ch06.py: … Label: …"
  lines: every cited test exists and mentions the function; every cited test's V-tag is among the listed levels (one
  apparent exception, `sphere_motion` citing `test_sphere_motion_V4_…` whose `def` comment reads "V1/V4/V7" — the
  listed V4 is carried, so it under-claims, which is honest); every listed label word matches its levels. Read by hand
  (11): `sphere_motion`, `example_6_2`, `source_panels` (now "exact on the circle, order ≈ 2 on the ellipse", off-body
  first order), `joukowski_inverse`, `solve_laplace`, `blasius_force`, `added_mass_sphere`,
  `ComplexFlow.stagnation_points`, `Uniform`, `cylinder`, `half_body` — all honest. Note: 116 public callables carry a
  `Label:` but no "Validation:" line (as before loop 1; not a regression) and four small helpers (`as_complex_point`,
  `ComplexFlow`, `BlasiusForce`, `polygon_signed_area`) have no `Label:` — a docstring nicety for the knowledge phase.
- **Tolerances**: `tests/test_ch06.py` was edited by no one but the verifier (its timestamp precedes the loop-2 code
  changes). No assertion was relaxed in loop 2. Two loop-1 bounds were set (not relaxed) against measured error sources:
  the cylinder `ideal_flow_residuals` bounds are the O(h²) stencil level (1e-5, with the h² rate itself asserted), and
  the late-time Example 6.1 check allows 1e-8 Pa because the exact closed form is 4.0e-9 Pa at t = 10⁶ s.

## F1 (loop 1, resolved in loop 2) — `Flow.stagnation_points` (default search) misses the half-body's stagnation point; `superposition_state` then reports "no stagnation point"
- Functions: `core.potential.ComplexFlow.stagnation_points` (Part C 1.6) and, through it, `ch06.superposition_state`
  (Part C 5.13, the Python reference of explainer E1 `superposition_sandbox`, its parity rows and status).
- Test: `test_stagnation_points_V1_default_search_finds_the_half_body_nose[m]`, m = 0.5, 1, 2, 3 m²/s (all four fail).
- Input / expected / got: `superposition_state([{"kind": "uniform", "U": 1.0}, {"kind": "source", "m": m}], 1.0, 1.0)`
  — expected `stagnation = [[−m/2πU, 0]]` = [[−0.0796, 0]], [[−0.159, 0]], [[−0.318, 0]], [[−0.477, 0]] and
  `psi_dividing = m/2`; **got `stagnation = []` and `psi_dividing = NaN`** for every one. The same call finds the point
  for m = 2π (the E1 default), which is why the implementer's checks passed.
- Cause (measured): the Newton iteration z ← z − f/f′ on dw/dz = U + m/2πz is e_{k+1} = (U/c) e_k² with c = m/2π, so
  its basin is the disc |z + c/U| < c/U = a (the stagnation distance itself). The default start grid is 9 × 9 points
  over (−3, 3)² (spacing 0.75 m); no start point lies within a of −a when a ≲ 0.5 m, and every start diverges. The
  design's "complex Newton on dw/dz **with deflation**" is not implemented, and the result is not checked against a
  second route. E1's source slider runs 0 … 20 m²/s (design Part B), so most of its range is affected; C05's slider over
  m/U likewise.
- Suggested fix (implementer's choice): seed Newton also from the singular points' neighbourhoods (e.g. starts at
  z_k ± ε·e^{iθ_j} around every source/vortex/doublet centre, ε a few percent of the box) and/or a finer default grid
  scaled by the element strengths; or solve dw/dz = 0 as a polynomial (clear the denominators of the rational dw/dz and
  use `np.roots`), which is exact for every Flow of sources, sinks, vortices, doublets and a stream.

## F2 (loop 1, resolved in loop 2) — design Part C row 6.1: three drawing helpers are missing from `scripts/ch06_drawings.py`
- Contract (design Part C C.6, row 6.1): `flow_net(...)`, `draw_body(ax, pts, **kw)`, `pressure_arrows(ax, pts, p,
  normals, scale)`, `separated_cp_band(theta_front_deg, sep_deg=80.0) -> (lo, hi)` (the **qualitative** band of N31 /
  Fig. 6.10, "labelled so in every figure").
- Test: `test_drawing_helpers_V1_part_c_6_1_exist_and_draw` — `flow_net` and `circle` exist and draw;
  **`draw_body`, `pressure_arrows`, `separated_cp_band` are absent** (AttributeError list printed by the test).
  The notebook's C06 figure (ideal vs qualitative real band) and E2's "real band" toggle depend on
  `separated_cp_band`; nothing else in the repository defines it.
- Suggested fix: add the three helpers (no physics; the band must carry a "qualitative, not measured" label).

## Discrimination proofs (wrong variants planted in a scratch copy; the repository untouched)
Each variant was written into a fresh copy of `fluidpy/` and `tests/test_ch06.py` was run with `-x` (scripts, drawing
and F1 tests excluded). **28 of 28 planted variants are caught** in loop 1, and again 28/28 on the loop-2 code.

| # | planted wrong variant (file) | first failing test |
|---|---|---|
| 1 | book Γ read as counterclockwise | `test_polar_velocity_V1_cylinder_components` |
| 2 | half-body principal θ branch | `test_no_through_flow_V1_bodies_are_streamlines` |
| 3 | Zhukhovsky inverse: numpy principal root | `test_blasius_V3_polygon_contours_converge_at_second_order` |
| 4 | (6.104) printed bracket sign | `test_moving_sphere_V1_kinematic_condition_and_surface_velocity` |
| 5 | cylinder doublet pointing downstream | `test_polar_velocity_V1_cylinder_components` |
| 6 | vortex image with the same sign | `test_images_V1_wall_conditions_and_wrong_variant` |
| 7 | 5-point Laplacian: Δy instead of Δy² | `test_laplace_5pt_V3_second_order_and_exact_on_quadratics` |
| 8 | Blasius factor iρ instead of iρ/2 | `test_lift_V1_form_matches_published_kutta_joukowski` |
| 9 | added mass = whole displaced mass | `test_added_mass_V1_force_quadrature_three_directions` |
| 10 | axial influence matrix sign | `test_axial_method_V3_ellipsoid_strengths_converge` |
| 11 | source panel self term −½ | `test_source_panels_V1_circle_is_exact_and_closed` |
| 12 | Example 6.1 fraction sign | `test_example_6_1_V1_closed_form_vs_numeric_route` |
| 13 | Stokes operator = Laplacian sign | `test_stokes_stream_function_V1_field_equation_is_not_laplace` |
| 14 | half-body C_p: −k² sign | `test_half_body_V1_surface_pressure_and_its_zero` |
| 15 | (6.49) scalar d as +d e_x | `test_doublet_limit_V3_pair_approaches_doublet_at_second_order` |
| 16 | flux without 2π | `test_axisym_flux_V1_two_pi_dpsi_and_source_strength` |
| 17 | circle flow: vortex sign flipped in w only | `test_elliptic_cylinder_V1_boundary_far_field_and_lift` |
| 18 | stagnation sin θ = +Γ/4πaU | `test_lift_V1_stagnation_points_closed_form_and_newton` |
| 19 | (6.105) 9/4 instead of 9/8 | `test_moving_sphere_V1_dphidt_and_pressure` |
| 20 | far-field doublet term sign in (6.61) | `test_kutta_zhukhovsky_V2_series_residue_and_book_slip` |
| 21 | Gauss–Seidel matrix (I + ωL) | `test_laplace_solver_V1_methods_agree_and_sweep_counts` |
| 22 | Example 6.2 outlet ψ over 2 m | `test_example_6_2_V4_flux_and_maximum_principle` |
| 23 | (6.41) curve: −x² in the quadratic | `test_two_sources_V1_streamline_equation_and_complex_form` |
| 24 | Laurent coefficients R^(−k) | `test_blasius_V1_laurent_coefficients_and_contributions` |
| 25 | source velocity m/4π | `test_delta_flux_V1_vortex_and_source_strengths_on_every_circle` |
| 26 | (6.101) ∂φ/∂x_s = +∇φ | `test_moving_sphere_V1_dphidt_and_pressure` |
| 27 | rankine oval half-length without the ma/πU term | `test_rankine_oval_V1_closed_body` |
| 28 | sphere doublet d = πa³U | `test_sphere_V1_stream_surface_velocity_and_cp` |

In-test discrimination (always run): printed (6.104) bracket adds exactly u_s; numpy's principal √(z² − 4b²) lands
inside the circle at every left-half point off the slit; the printed (6.61) 1/z² coefficient differs from the correct
−(Ud/π + Γ²/4π²); ψ = +(Γ/2π) ln r gives flux
+Γ; a source written m/4π gives m/2; a same-sign vortex image leaks > 0.1 m/s through the wall; `Gamma_ccw = +2`
gives L = −24 N/m; numpy's principal half-body branch puts ψ = −m/2 on the lower body; dropping the 2π in (6.78)
misses by 6.28×; the reversed dipole misses the pair's far field by > 0.1; a Poiseuille field keeps μ∇²u = −G.

## Validation table — A items (CORE, ≥ 2 independent levels, one of V1/V2/V3/V5)
| Concept / Eq. | Tier | fluidpy target | Evidence (V-levels) | Numbers (error, order, residual) | Label | Notes |
|---|---|---|---|---|---|---|
| C01 ideal-flow equations (6.1); N01–N04, R01; D01 | CORE | `ideal_flow_residuals`, `ideal_flow_applicability`, `flow_field_callables` | V1 cylinder field: ∇·u ≤ 1e-5 U/a, Euler/inertia ≤ 1e-5, μ∇²u ≤ 1e-5 μU/a², ω ≤ 1e-5 U/a, all falling as h² (order 2.0 ± 0.15); corner exact; Poiseuille control μ∇²u = −G = Euler residual, ω = Gy/μ; half-body Euler residual ≤ 1e-7 of inertia; applicability table (Re = 50, M = 0.5, baroclinic → "Kelvin fails", wake…); V2 D01 sympy (curl-of-curl identity for a generic field; harmonic φ ⇒ ∇·u = ω = μ∇²u = 0; Bernoulli p balances (u·∇)u; Poiseuille keeps −G); pint N/m³ | stencil order 2.0 | analytic, symbolic | N01 number Re = 6.7e5, δ ≈ 1.2 mm |
| C02 ψ, φ, ω = −∇²ψ (6.4)–(6.15), (6.21)–(6.22); N05–N15, N18–N19, R03, R08–R09; D02–D04 | CORE | `vorticity_from_psi(_sym)`, `rankine_vortex_psi`, `laplacian_residual`, `delta_flux_check`, `orthogonality_check/report`, `polar_velocity(_sym)`, `Uniform`, `Source`, `Vortex` | V1 ω of sin x sin y (1e-8), Rankine core Γ/πa² (1e-9) and 0 outside, ψ continuous with ch03 u_θ; every element harmonic (≤ 2e-6 with h = 1e-3), NaN at centres; flux −Γ and m on radii 0.01…100 (1e-13, V4 radius independence), 0 when missed; ∇φ·∇ψ ≤ 1e-7 |∇φ|², |∇φ|/|∇ψ| − 1 ≤ 1e-7; polar = rotated Cartesian (1e-9); element forms (6.7), (6.8), (6.14), (6.15) and velocities vs ch05 Biot–Savart (1e-13); V3 9-point residual order 3.999; V2 sympy twins, D02/D03/D04 (∇² ln r = 0, flux 2π, ∇ψ = e_z×∇φ and the analysis' reversed form fails) | order 3.999 | analytic, converged, symbolic | N10 Poisson solver order 2 (manufactured) |
| C03 superposition, no through-flow (6.16)–(6.18); N16, N17, R05; D05 | CORE | `Flow`, `normal_velocity_on/values`, `far_field_check`, `Flow.pressure`, `pressure_coefficient`, `stream_function`, `velocity_potential` | V1 w, dw/dz, u of a 5-element sum = sums (1e-13), contributions; u·n on r = a ≤ 1e-12 for 9 Γ and ψ = 0 there; half-body u·∇ψ ≤ 1e-8 and ψ = m/2 on both halves (1e-12); V7 far-field decay slopes −2.00 (cylinder, oval) and −1.00 (half-body, circulation) ± 0.05; V4 one Bernoulli constant on the field (ptp ≤ 1e-9), C_p = 1 at every Newton stagnation point; V2 D05 sympy | slopes −2.000, −1.000 | analytic, conserved, symbolic | F1 affects the default stagnation search |
| C04 doublet limit (6.28)–(6.29); N20–N26, R10; D06 | CORE | `source_sink_pair`, `doublet_limit_error/frames`, `Doublet`, `Doublet.from_book_scalar`, `Corner`, `Constant`, `harmonic_polynomials` | V3 pair → doublet error order 2.011 (pairwise 2.03, 2.01, 2.00); V1 D06 numbers (0.318352 vs 0.318310, gap ε²/3 within 2 %), (6.29)/(6.49) forms (1e-13), dipole from sink to source (reversed vector misses by > 0.1); corner family (6.24)–(6.27) incl. rotate and complex A; V2 D06 sympy (Taylor, ε² terms cancel, limit with 2mε fixed); harmonic polynomials 2 per degree | order 2.011 | converged, analytic, symbolic | |
| C05 half-body (6.30)–(6.32); N27–N29; D07, D08 | CORE | `half_body`, `half_body_numbers/shape/surface_cp/cp_zero_angle/net_force`, `rankine_oval`, `superposition_state`, `flow_from_spec` | V1 stagnation −m/2πU (Newton from a nearby start, 1e-12), ψ = m/2 on the body (1e-12·m), h(θ) (1e-12), (6.30)/(6.50) forms; surface C_p reduction = field C_p (1e-12), zero at 113.218°, +1 at the nose, −0.4053 at 90°, min −0.587 at 63.0°; V4 flux inside the body = m at x = 0.5, 5, 200 m (1e-9), h → m/2U; V3 net force → 0 as (1/x_end)^1.97; V2 D07/D08 sympy step by step; V1 Rankine oval (ψ = 0 body, root equation by sympy); V6 Fig. 6.8 angle | order 1.97 (force) | analytic, conserved, converged, symbolic, book-value | F1 (default search) |
| C06 cylinder, d'Alembert (6.33)–(6.35); R12–R14, N30, N31; D09 | CORE | `cylinder`, `cylinder_surface_cp/pressure/speed`, `surface_pressure_force` | V1 parity with `ch03.cylinder_flow` (body and fluid frames) and `cylinder_streamfunction` (1e-13), (6.33)–(6.34) (1e-13), C_p = 1 − 4 sin²θ, symmetric; D = L = 0 by the spectral pressure integral (≤ 1e-11); N30 +60/−180/+60 Pa; shoulder speed 2U; V1 form cross-check (Wikipedia, cited); V2 D09 sympy | — | analytic, symbolic | N31 real band qualitative (F2: helper missing) |
| C07 lift with circulation (6.36)–(6.40); N32–N38; D10, D11 | CORE | `cylinder(Gamma_cw/ccw)`, `cylinder_stagnation_points`, `lift_per_span`, `circulation_family_check`, `cylinder_circulation_state` | V1 L = ρUΓ_cw by the pressure integral for 6 Γ incl. > 4πaU (≤ 1e-10), D = 0; Γ_ccw = +2 → L = −24 (WV); stagnation closed form = Newton (1e-12), −9.158°/−170.842°, merged at −a, off-body r₊ = 0.2618, r₋ = 0.0382, r₊r₋ = a² (1e-12); (6.37), (6.39) vs the field; V4 loop circulation −Γ_cw on every radius (1e-10); regimes of the E2 state; V2 D10/D11 sympy (Vieta product, term-by-term integrals); pint N/m; V1 form cross-check (Wikipedia K–J, clockwise contour) | — | analytic, conserved, symbolic | |
| C08 images, Example 6.1; N39, N40, R15–R17; D12, D13 | CORE | `mirror`, `two_sources`, `two_source_streamline`, `circle_theorem`, `example_6_1(_times)`, `example_6_1_wall_pressure` | V1 wall conditions on 1001 points ≤ 1e-14 (vortex, source, doublet; walls y = 0, x = 0), same-sign image leaks (WV), parity with ch05 wall images; (6.41) curve on ψ mod m/2 (1e-10), (6.53) Im w = ψ mod m (1e-12); circle theorem u·n ≤ 1e-10 and parity with ch05 circle images (1e-12); Example 6.1 closed vs numeric route ≤ 1e-6·ρΓ²/4π², path vs ch05 evolver 1e-8 m, −25.330 Pa, zero at 4π s, max 3.1663 Pa at 21.766 s; V3 ∂φ/∂t differences order 1.9997; V4 ξ_x conserved (ptp < 1e-9); V2 D12, D13 sympy (incl. maximum at s = √3h) ; V6 printed closed forms | order 1.9997 | analytic, converged, conserved, symbolic, book-value | |
| C09 complex potential (6.42)–(6.53); N41–N52; D14, D15 | CORE | `Flow.w/dwdz`, `cauchy_riemann_residual(_sym)`, `complex_potential_probe`, `complex_potential_family`, `corner_info`, `corner_speed_exponent`, `Corner` | V1 every element's w vs (6.46)–(6.52) on 2000 points (1e-12); CR residual ≤ 1e-6 for 7 families × 3 points, z* → (2, 0) (discriminates); dw/dz = u − iv vs differences of ψ (1e-7); N43 number u = 2, v = −2; corner walls ψ = 0 (1e-12), log–log slope n − 1 (± 0.01) for n = 4…½, no jump inside the wedge; V2 D14/D15 sympy | slopes exact to 0.01 | analytic, symbolic | |
| C10 Blasius and Kutta–Zhukhovsky (6.54)–(6.62); N53–N63, R18; D16–D18 | CORE | `blasius_force(_report)`, `laurent_coefficients`, `laurent_contributions`, `blasius_state`, `cv_force_on_body`, `contour_force`, `complex_force_from_pressure`, `kutta_zhukhovsky_sym`, `force_on_held_singularity` | V1 (0, 24) N/m on R = a…100a (1e-12), clockwise contour rejected, crossing flagged; c₀ = U, c₋₁ = iΓ_cw/2π, c₋₂ = −Ua² (1e-12), shape-free c₋₁ for the ellipse; only the 1/z bar is nonzero; V4 three routes (CV on R = 5a…50a, surface pressure, Blasius) agree to 1e-10 for cylinder and ellipse; tilted ellipse force ⟂ stream, size ρUΓ; Gauss checks p = y → L = −A (1e-12); held source D = −ρmU; V3 polygon Blasius order 1.9965 (circle) and 1.9994 (square round the ellipse); V2 D16, D17, D18 sympy incl. the printed (6.61) slip; V1 form cross-check (Wikipedia Blasius) | orders 1.997, 1.999 | analytic, conserved, converged, symbolic | |
| C11 conformal mapping, Zhukhovsky (6.63)–(6.69); N64–N70; D19–D21 | CORE | `angle_preservation`, `map_elements(_info)`, `grid_image`, `map_grid_lines`, `conformal_map`, `joukowski(_derivative)`, `joukowski_inverse(_derivative)`, `circle_flow_zeta`, `mapped_flow`, `joukowski_ellipse`, `elliptic_cylinder_flow`, `joukowski_state`, `ellipse_surface_speed`, `cot_flow` | V1 angles kept (1e-8) for 6 maps, doubled (z²) and tripled (z³) at 0; slit, ellipse (1e-13), foci 2b; inverse round trip on 10⁴ random ζ in all quadrants (1e-12), |ζ| ≥ b on a grid, principal root wrong exactly for Re z < 0 (WV); u·n on the ellipse ≤ 1e-12, far field 1/R, Blasius ρUΓ, b → 0 limit (1e-9); exact ellipse speed = mapped-flow speed (1e-11, independent); V2 D19–D21 sympy (Vieta, branch check, chain rule); V1 form cross-check (Wikipedia Joukowsky, W = W̃/(1 − 1/ζ²)) | — | analytic, symbolic | D21 check line: see O1 |
| C12 finite-difference Laplace (6.70)–(6.73), Example 6.2; N71–N74, R19–R21; D22, D23 | CORE | `laplacian_5pt`, `node_update`, `jacobi/gauss_seidel/sor_sweep`, `residual_norm/field`, `sweep_order`, `assemble_laplace`, `solve_laplace`, `solve_poisson`, `optimal_sor_omega`, `jacobi_spectral_radius`, `four_point_system`, `example_6_2(_geometry)`, `relaxation_state` | V3 5-point order 2.0000 (fixed point, truncation = (Δx²/12)·2π⁴ψ to 1e-3), solver order 1.995 on sin πx sinh πy; V1 exact on xy, x² − y²; 4-point (1, 2, 2, 4) by every method; sweeps 34/19/12 (design D23 check), ω_opt 1.0718; first GS sweep 0, .75, .75, 3.375; V1 flux identity Q through every section (1e-9; relabelled from V4 in loop 3); V4 cell circulation ≤ 4 tol = 4e-10 converged (3.3e-10), 0.0999 after 5 sweeps, < 1e-12 direct at Δ = ¼; V1 public geometry/BCs (24 unknowns); discrete max principle and average rule on every node; V3 Example 6.2 refinement order 1.36/1.39/1.30 at three points → **4/3 (corner pollution)**; V2 D22/D23 sympy (Δx²/12 ψ_xxxx, ρ_J = ½, ρ_GS = ¼); V6 Fig. 6.25 values | orders 2.000, 1.995, 1.36 (→ 4/3) | converged, analytic, conserved, symbolic, book-value | O2 (design said "≈ 2 away from the corner") |
| C13 axisymmetric flow, sphere (6.74)–(6.92); N75–N85, R22–R30; D24, D25 | CORE | `AxisymUniform`, `PointSource3D`, `Doublet3D`, `LineSource3D`, `AxisymFlow`, `sphere`, `axisym_velocity_spherical(_sym)`, `sphere_potential_vector`, `stokes_operator_residual(_sym)`, `stokes_operator_sym`, `axisym_laplacian_residual/sym`, `spherical_continuity_residual`, `axisym_flux_between`, `sphere_surface_cp`, `cylinder_vs_sphere`, `perturbation_radius`, `axisym_half_body` | V1 ψ = 0 on r = a (1e-11) and the axis, (6.89)–(6.91) (1e-13), (6.92) vector form; parity with ch05's Hill exterior (1e-12) and ch04's axisymmetric velocity tool; ψ and φ routes of (6.83) agree; (6.77) residual ≤ 1e-5 for 5 fields, R z → −z/R², R² z a solution; (6.80) for φ; App. B continuity 0 and the printed (6.82) = r × App. B (sympy, loop 3; not a slip — review M1); flux = 2πΔψ (1e-12, WV without 2π); V4 source flux Q through every sphere (1e-12); V3 3-D pair → doublet order 2.002; V7 relief table; V2 D24/D25 sympy; V1 form cross-check McDonald (2015) | order 2.002 | analytic, conserved, converged, symbolic | |
| C14 airship, axial method; panels (N86–N90); D26, D27 | CORE | `line_sink_stream_function`, `airship`, `axisym_body_target`, `axisym_body_fit`, `axial_state`, `ellipsoid_linear_source_strength`, `axial_influence_matrix`, `axial_singularity_solve/psi/velocity`, `source_panels`, `panel_geometry`, `panel_induced_velocity`, `panel_velocity`, `panel_cp_error` | V1 (6.94) = quad of (6.93) (1e-10), LineSource3D ψ, φ, u consistent (1e-7); airship ψ = 0 on the body (1e-10), closure 0, nose −0.2521, tail 1.0696, L = 1.3217 m, our cubic (z/a)²(z/a − 1) = ±Q/4πUa² (1e-9); V1+V3 ellipsoid k_n → Kξ at order 1.064 (N = 16…64), body ψ error order 2.97; Rankine oval: dipole moment error 1.4e-3, 1.6e-4, 2.9e-6 (N = 10, 20, 40), Σk_nΔξ ≤ 1e-5, cond ↑; V1/V4 panels: circle C_p exact (≤ 1e-12, every N), Σλ_jS_j ≤ 1e-12, rotated stream; V3 ellipse C_p order 2.057; off-body velocity order 0.93 (→ 1); V2 D26 sympy (substitution, (6.94), far-field cancellation, axis cubic) | orders 1.064, 2.057, 0.93 | analytic, converged, conserved, symbolic | sphere/airship targets and N > 40 fenced qualitative (O3) |
| C15 accelerating sphere, added mass (6.96)–(6.109); N91–N106, R31, R32; D28–D31 | CORE | `moving_sphere_potential/velocity/surface_velocity/dphidt/surface_pressure`, `sphere_force_quadrature`, `added_mass_sphere`, `added_mass_by_energy`, `cylinder_added_mass`, `moving_sphere_state`, `added_mass_state`, `moving_sphere_dphidt_sym`, `sphere_motion`, `rayleigh_collapse_time` | V1 n·∇φ = n·u_s on |ξ| = a (1e-13), ∇φ = differences of φ (1e-8), (6.104) = (6.103) on the surface, printed bracket off by u_s (WV); ∂φ/∂t = time differences for a turning motion (1e-8); (6.105) = unsteady Bernoulli (ch04 tool, 1e-10) = ch04's accelerating sphere (1e-10); steady part = (6.91) with θ ↔ π − θ; force quadrature: speed part 0, acceleration part −M du/dt for 3 × 2 cases (1e-10); V1 energy route (surface and volume) = M (1e-10/1e-9); cylinder ρπa² by energy and force; V5 coefficient ½ (Wikipedia), Rayleigh 0.91468 (Wikipedia, 5 s.f.); V4 work = kinetic energy (1e-8), u = Ft/(m + M), bubble 2g, steel ball 34.77 kg; V2 D28–D31 sympy (design's concrete motion, surface and volume integrals); pint | — | analytic, benchmark, conserved, symbolic | |

## Validation table — coded B/C items (NOTE, ≥ 1 level) and recaps
| Items | fluidpy target | Evidence | Label |
|---|---|---|---|
| N01–N04, R01 (applicability, BCs, exclusions) | `ideal_flow_applicability` | V1 decision table incl. Kelvin (WV), N01 numbers; V6 threshold 10³ is the default | analytic, book-value |
| N05, N06, N10–N12, N26 (Laplace, deltas, Poisson) | `laplacian_residual`, `delta_flux_check`, `solve_poisson` | V1 residual table, fluxes; V3 Poisson order 2 | analytic, converged |
| N07, N08, N13, N14, N24, R11 (elements, velocities) | `Uniform`, `Source`, `Vortex` | V1 forms and ch05 parity (1e-13) | analytic |
| N09, N18, N19, R03–R09 (orthogonality, polar forms, recaps) | `orthogonality_check`, `polar_velocity(_sym)`, `core.streamfunction`, `core.curvilinear` | V1/V2 | analytic, symbolic |
| N15 (named) | — | pointer only | — |
| N16, N17, R05 | `normal_velocity_on`, `far_field_check`, `Flow.pressure` | V1, V7, V4 | analytic, conserved |
| N20–N23, N25, R10 | `harmonic_polynomials`, `Corner`, `source_sink_pair` | V2, V1 | symbolic, analytic |
| N27–N29 | `half_body`, `half_body_surface_cp`, `half_body_cp_zero_angle`, `half_body_net_force` | V1, V3, V6 | analytic, converged, book-value |
| N30, R12–R14 | `cylinder_surface_cp/pressure`, `cylinder`, ch03 parity | V1 | analytic |
| N31 (qualitative band) | `separated_cp_band` (loop 2) | V1 smoke: lo ≤ hi, contains (6.35) ahead of separation, negative wake plateau | qualitative (Open O6) |
| N32–N38 | `cylinder(Gamma_cw)`, `cylinder_surface_speed/pressure`, `cylinder_stagnation_points`, `circulation_family_check` | V1, V4, V2 | analytic, conserved, symbolic |
| N39, N40, R15–R17, N51 | `mirror`, `two_sources`, `two_source_streamline`, `circle_theorem` | V1, V2 | analytic, symbolic |
| N41–N50, N52 | element `w`, `Flow.w`, `cauchy_riemann_residual(_sym)`, `complex_potential_probe` | V1, V2 | analytic, symbolic |
| N53–N63, R18 | `contour_force`, `complex_force_from_pressure`, `blasius_force`, `laurent_coefficients`, `cv_force_on_body`, `kutta_zhukhovsky_sym` | V1, V4, V3, V2 | analytic, conserved, converged, symbolic |
| N64–N70 | `map_elements`, `angle_preservation`, `grid_image`, `joukowski*`, `mapped_flow`, `elliptic_cylinder_flow`, `cot_flow` | V1, V2 | analytic, symbolic |
| N71 (named); N72–N74; R19–R21 | `four_point_system`, `solve_laplace`, `example_6_2` | V1, V3, V4, V6 | analytic, converged, conserved, book-value |
| N75–N85, R22–R30 | Stokes operator, axisymmetric elements, `sphere`, `axisym_flux_between`, `sphere_potential_vector` | V1, V2, V3, V4 | analytic, symbolic, converged, conserved |
| N86–N90 | `line_sink_stream_function`, `airship`, axial method, `source_panels` | V1, V3, V4 | analytic, converged, conserved |
| N91–N105, R31, R32 | moving-sphere family, `added_mass_*`, `sphere_motion` | V1, V2, V4, V5 | analytic, symbolic, conserved, benchmark |
| N106, N37, N57, N76 (named) | — | pointers only | — |
| Exercises (optional functions) | `rankine_oval`, `force_on_held_singularity`, `perturbation_radius`, `axisym_half_body`, `cylinder_added_mass`, `rayleigh_collapse_time`, `ellipsoid_linear_source_strength` | V1, V2, V5, V6 | analytic, benchmark, book-value |

## Derivations (curation §4b, design Part F) — every ★★ and ★★★ re-derived with sympy
| D | Star | Test | What is re-run |
|---|---|---|---|
| D01 | ★★ | `test_ideal_flow_V2_derivation` | curl-of-curl for a generic field (steps 4–5), harmonic φ ⇒ μ∇²u = 0 (6), Bernoulli balances Euler (7), Poiseuille control |
| D03 | ★★ | `test_delta_flux_V2_derivation` | steps 1–4 (Cartesian and polar), 7–9 constants |
| D06 | ★★ | `test_doublet_limit_V2_derivation` | step 1 factoring, series: ε¹ = (mε/π)x/r², ε² = 0, ε³ ≠ 0, limit with 2mε fixed, polar form, dipole vector |
| D07, D08 | ★★ | `test_half_body_V2_derivation` | every step (velocity, x_S, ψ_S, h(θ), limits, mass balance; k, |u|²/U², C_p, factorisation) |
| D10, D11 | ★★ | `test_lift_V2_derivation` | clockwise ψ, (6.36)–(6.38), Vieta r₊r₋ = a²; square expansion, lift and drag integrals |
| D13 | ★★ | `test_example_6_1_V2_derivation` | ∂φ/∂t and v at the origin, one fraction, the maximum at s = √3h, drift |
| D14, D15 | ★★ | `test_complex_potential_V2_derivation` | quotients along x and iy, 1/i = −i, CR, u − iv, Laplace, orthogonality; polar zⁿ, walls, n = 2, derivative with α |
| D16 | ★★ | `test_force_V2_derivation_D16` | outward normal, −∮p n ds = (D, L), circle consistency with D11 |
| D17 | ★★★ | `test_blasius_V2_derivation` | step 2 identity, step 9 on the body, ∮dz* = 0, (6.60) = −iρUΓ on r = a and r = 2a, pressure route = Blasius |
| D18 | ★★★ | `test_kutta_zhukhovsky_V2_derivation` (+ `…_series_residue_and_book_slip`) | w terms, square, true vs printed 1/z² coefficient, ∮z⁻ⁿdz, residue, D = 0, L = ρUΓ, open-body thrust |
| D19, D20 | ★★ | `test_conformal_V2_derivation` | turn cancels, doubling at 0; slit, ellipse, foci, critical points, two-to-one |
| D21 | ★★★ | `test_joukowski_inverse_V2_derivation` | roots, Vieta, the design's branch grid, far-field ζ ≈ z, chain-rule velocity, zeros of 1 − b²/ζ² |
| D23 (+ D22 ★) | ★★ | `test_laplace_V2_derivation` | Taylor Δx²/12 ψ_xxxx; A, b, diagonal dominance, eigenvalues of N, ρ_J = ½, ρ_GS = ¼, k = 27/14, first sweep |
| D24, D25 | ★★ | `test_axisym_V2_derivation` | cross products, (6.75), ω_φ = −(6.77), R·(6.77) vs ∇²ψ; source flux, φ, ψ, 3-D doublet limit, d = 2πa³U, (6.90), (6.91), (6.85) |
| D26 | ★★ | `test_airship_V2_derivation` | (6.93) integrated = (6.94), dξ = R dα/sin²α, the substituted integral, axis cubic, far field |
| D29 (+ D28 ★) | ★★★ | `test_moving_sphere_V2_derivation` | `moving_sphere_dphidt_sym`; the design's motion x_s = (t²/2, 0, t): direct ∂φ/∂t = chain rule, (6.104) and the printed bracket, (6.105) = −∂φ/∂t − ½|∇φ|², steps 10–12, (6.106) |
| D30, D31 | ★★ / ★★★ | `test_added_mass_V2_derivation` | (6.107)–(6.108) integrals, speed part zero for any u_s; volume = surface = 2πU²a³/3, M, ∇²φ = 0, far surface ~ R⁻³, surface values |
★ rows also checked: D02 (`test_psi_vorticity_V2_sympy_twin`), D04, D05, D09, D12, D22, D28; D27 numerically.
No intermediate line of Part F was found wrong; one **check line** is: D21's "print(bad.mean()) (→ 1.0)" gives 0.998 (O1).

## Functions used by the notebook and explainers (design Part C) — test name each
All 100 Part C rows are exercised (loop 2: the three drawing helpers of F2 now exist and are tested). C.1 `core.potential`: elements,
`Flow`, bodies, checks, forces, axisymmetric and moving-sphere families — C02–C15 tests above plus
`test_primitives_V1_branches_points_and_function_flows` (`log_branch`, `power_branch`, `d2wdz2`, `FunctionFlow`,
`singular_points`); C.2 `core.conformal` — C11 tests; C.3 `core.laplace_solvers` — C12 tests; C.4 `core.panels` — C14
tests; C.5 chapter module — every name appears in a named test (`test_exports_V1_every_part_c_name_is_reachable`
asserts the API); Part B parity rows — `test_parity_rows_V1_design_part_b_expressions_run` (26 rows evaluated with the
design's literal arguments); C.0 reuse rows called with chapter-6 inputs — `test_reused_tools_V1_smoke_with_chapter_6_inputs`
(`ch01.fluid_properties`, `ch04.is_incompressible_regime`, `ch05.kelvin_hypotheses_text`, `gauss_legendre_nodes`,
`viscous_irrotational_residual`, `interact.slider_figure`, `point_vortex_evolve`); `anim`, `embed`, `style` are covered
by `tests/test_machinery.py`; C.6 scripts — `test_scripts_V1_every_ch06_script_runs`; C.6 drawing helpers —
`test_drawing_helpers_V1_part_c_6_1_exist_and_draw` (passes in loop 2).

## Convergence studies (scheme | steps | observed order (pairwise) | design order)
| Scheme | Steps | Observed | Design |
|---|---|---|---|
| source–sink pair → doublet (6.28)→(6.29) | ε = 0.2 … 0.025 | 2.011 (2.03, 2.01, 2.00) | 2 |
| 3-D pair → doublet (6.88) | ε = 0.1 … 0.0125 | 2.002 (2.004, 2.001, 2.000) | 2 |
| 9-point `laplacian_residual` | h = 0.1 … 0.0125 | 3.999 (3.999, 4.000, 3.999) | 4 |
| 5-point `laplacian_5pt` (fixed point) | Δ = 1/16 … 1/128 | 2.000 (2.000, 2.000, 2.000) | 2 |
| `solve_laplace`, harmonic sin πx sinh πy | Δ = 1/8 … 1/64 | 1.995 (1.988, 1.997, 1.999) | 2 |
| Poisson (manufactured) | Δ = 1/8 … 1/64 | ≈ 2 (asserted ± 0.15) | 2 |
| **Example 6.2 refinement** (probe far (1, 4), mid (2, 2), near (5.5, 2.5)) | Δ = ½ … 1/64 | 1.36 (1.41, 1.38, 1.37, 1.35), 1.39 (1.50, 1.44, 1.40, 1.37), 1.30 (1.12, 1.24, 1.29, 1.31) | **4/3** = 2n for the 270° corner (n = ⅔ in (6.46)); the design's "≈ 2 away from the corner" does not hold (O2) |
| Blasius on a polygon (circle R = 0.3, N = 16 … 128) | segment trapezoid | 1.997 (1.992, 1.998, 1.999) | 2 |
| Blasius on a square round the Zhukhovsky ellipse | N = 8 … 64 per side | 1.999 (1.998, 2.000, 2.000) | 2 |
| Example 6.1 numeric ∂φ/∂t | Δt = 0.4 … 0.05 s | 1.9997 | 2 |
| constant-source panels, ellipse C_p (A = 1, B = 0.5) | N = 32 … 256 | 2.057 (2.05, 2.08, 2.03) | 2 |
| constant-source panels, circle C_p | N = 8 … 128 | exact (≤ 1e-14 at every N) | — (the curation's "order 2 on the circle" does not apply) |
| panels, off-body velocity (circle) | N = 32 … 256 | 0.93 (0.88, 0.94, 0.97) → 1 | not stated (O4) |
| axial method, ellipsoid k_n (A/B = 5) | N = 16, 32, 64 | 1.064 (1.117, 1.011) | 1 (end effects of constant segments) |
| axial method, ellipsoid body ψ error | N = 16, 32, 64 | 2.97 (2.74, 3.20) | — (reported) |
| half-body net force vs 1/x_end | x_end = 10, 100, 1000 m | 1.97 (1.93, 2.00) | 2 (tail C_p ~ 1/x) |
| far-field decay (V7) | R = 10 … 80 | −2.000 closed, −1.000 open/circulation | −2, −1 |
| `ideal_flow_residuals` stencil | h = 4e-4 … 1e-4 | 2.0 | 2 |

Condition numbers of the axial method (1-norm): ellipsoid 3.1e3, 9.2e5, 4.1e10 (N = 16, 32, 64) and 5.9e17 at N = 128
(k_error 15 — breakdown); Rankine oval 2.7e3, 3.3e6, 2.3e12, 3.5e17 (N = 10, 20, 40, 80; Σk_nΔξ = 0.134 at N = 80);
sphere 3.8e6, 4.1e13, 1.8e18 (N = 10, 20, 40).

## Conservation / invariant residuals
| Invariant | Residual |
|---|---|
| flux of ∇ψ (vortex) / ∇φ (source) through circles r = 0.01 … 100 m | ≤ 1e-13 relative (radius independent) |
| mass balance inside the half-body at x = 0.5, 5, 200 m | ≤ 1e-9 relative |
| Bernoulli constant over a field (cylinder with Γ) | ptp ≤ 1e-9 relative |
| loop circulation −Γ_cw on r = 1.5a, 3a, 10a (4 Γ) | ≤ 1e-10 m²/s |
| CV force = pressure force = Blasius (R = 5a … 50a) | ≤ 1e-10 relative |
| Example 6.2: discrete circulation round every control cell (converged, tol 1e-10) | 3.3e-10 ≤ 4 tol (the Q-flux through each section is an identity of the BCs, V1) |
| 3-D source: Q through spheres r = 0.1, 1, 10 m | ≤ 1e-12 relative |
| Example 6.1: ξ_x = h (impulse) | ptp ≤ 1e-9 m |
| sphere motion: work of F_E = ½(m + M)u² | ≤ 1e-8 J |
| panels Σλ_jS_j (closed body) | ≤ 1e-12 |

## Benchmarks used (value | our value | source + URL | date verified)
| Value | Ours | Source | Verified | Label |
|---|---|---|---|---|
| added mass of a sphere = ½ displaced mass, (2/3)πr³ρ | 0.5 (1e-14); energy route 0.5 (1e-10) | primary: Lamb, *Hydrodynamics*, 6th ed. (1932), Ch. V §92 (as cited by the review; text not fetched); secondary (fetched): Wikipedia "Added mass", https://en.wikipedia.org/wiki/Added_mass | 2026-09-23 | V5 |
| Rayleigh/Besant empty-cavity collapse 0.91468 R₀√(ρ/P∞) | 0.914681 (Beta-function closed form = quad, 1e-10) | primary: Rayleigh, Phil. Mag. (6) 34, 94–98 (1917), doi:10.1080/14786440808635681; secondary (fetched): Wikipedia "Rayleigh–Plesset equation" | 2026-09-23 | V5 |
| cylinder ϕ, V_r, V_θ, C_p ∈ [−3, +1], zero drag | identical (1e-12) | Wikipedia "Potential flow around a circular cylinder" | 2026-09-23 | V1 form |
| Kutta–Joukowski L′ = ρVΓ (clockwise contour) | 24 N/m = ρUΓ_cw (1e-12) | Wikipedia "Kutta–Joukowski theorem" | 2026-09-23 | V1 form |
| Blasius F_x − iF_y = (iρ/2)∮(dw/dz)²dz | our D − iL (1e-12) | Wikipedia "Blasius theorem" | 2026-09-23 | V1 form |
| Joukowsky: unit circle → [−2, 2], W = W̃/(1 − 1/ζ²) | identical | Wikipedia "Joukowsky transform" | 2026-09-23 | V1 form |
| sphere Φ, u_r, u_θ, cylindrical components, Stokes ψ sign | identical (1e-12) | K. T. McDonald, Princeton (2015), http://kirkmcd.princeton.edu/examples/sphereflow.pdf | 2026-09-23 | V1 form |
| Rankine half-body ψ, b = m/2πU, body r(θ) | identical | Wikipedia "Rankine half body" | 2026-09-23 | V1 form (in `reference/`; the tests use the same forms via D07) |
| constant-strength source panels | method only | Hess & Smith, Prog. Aerosp. Sci. 8, 1–138 (1967), doi:10.1016/0376-0421(67)90003-6 | 2026-09-23 | pedigree |

## Numbers from the text (book vs ours)  [private values redacted to relative errors]
| Book item | Agreement |
|---|---|
| §6.1 Re threshold | used as the default `Re_min` (exact) |
| §6.3 closed forms (half-body a, ψ_body, h(θ), h_max; cylinder d, C_p, u_θ(a), sin θ_s, r₊, critical Γ, L; (6.41) curve) | identical (≤ 1e-12 relative) |
| Fig. 6.8 zero-C_p angle (printed to the degree) | 0.19 % (113.218° from tan θ = −2(π − θ)) |
| Fig. 6.10 ideal C_p range | exact |
| Example 6.1 closed forms (trajectory, ∂φ/∂t, v, p) | identical (≤ 1e-12) at three (Γ, h, t) |
| Example 6.2: 24 printed grid values of Fig. 6.25 | all within the printed 2-decimal rounding (max |Δ| 0.0049 < 0.005) at the text's sweep count and when converged; at the code's sweep count 3 of 24 differ in the last digit (max |Δ| 0.0059) — as the analysis noted |
| §6.5 residue and the printed (6.61) coefficient | residue identical; printed coefficient ≠ correct (the slip) |
| §6.6 ellipse semi-axes, foci, slit | identical |
| §6.8 sphere ψ, φ, u_r, u_θ, C_p, d; line sink, airship ψ | identical (≤ 1e-12) |
| §6.9 added mass and fraction ½; steady pressure part | identical |
| Exercises 6.12, 6.18, 6.19, 6.34, 6.42, 6.43, 6.45, 6.49 answers | identical (≤ 1e-9) |

## Figures reproduced with our code (outputs/ch06/verify/, local)
| Figure | PNG | Visual verdict |
|---|---|---|
| Figs. 6.7–6.8 | `fig6_7_6_8_half_body_verify.png` | Blunt open body through the stagnation point at x = −1 m, widening toward the dotted h_max = π line; surface C_p runs from −0 far downstream through a −0.59 minimum near 63° to 0 at 113.2° and +1 at the nose. |
| Fig. 6.10 (ideal only) | `fig6_10_cylinder_cp_verify.png` | Symmetric curve from +1 at 0° to −3 at 90° and back to +1 at 180°, crossing zero at 30° and 150°. |
| Fig. 6.12 | `fig6_12_cylinder_circulation_verify.png` | Stagnation points at 0°/180°, then at −30°/−150°, merged at the bottom for Γ = 4πaU and free below the body for 1.5 × 4πaU; the printed L equals ρUΓ in every panel. |
| Example 6.1 | `example6_1_wall_pressure_verify.png` | Suction −25.3 Pa at t = 0 rising through zero at 12.57 s to a +3.17 Pa maximum at 21.77 s and decaying; the numeric route's dots sit on the closed form. |
| Fig. 6.25 | `fig6_25_contraction_verify.png` | Streamlines crowd round the 270° step corner and are straight and evenly spaced in the inlet and outlet; the book grid and the ×8 grid agree away from the corner. |
| Fig. 6.21 | `fig6_21_joukowski_verify.png` | Flow round the flat ellipse with both stagnation points slightly below the ends (clockwise Γ); the principal-root failure region is exactly the left half-plane. |
| Figs. 6.27–6.28 | `fig6_27_6_28_sphere_airship_verify.png` | Fore–aft symmetric sphere streamlines; the airship hull closes from −0.25 m to 1.07 m with a blunt nose at the source and a tapered tail over the line sink. |
| (6.105) | `eq6_105_sphere_pressure_verify.png` | Speed part fore–aft symmetric (500 → −625 → 500 Pa), acceleration part ±100 Pa, total 600 Pa in front, 400 Pa behind. |

## Deviations & justifications
The implementation carries no `# DEVIATION` marker; the documented corrections of book slips are verified as follows:
(6.61) coefficient (code returns the correct −(Ud/π + Γ²/4π²), the printed one kept separately and shown to differ);
(6.104) bracket (`printed_bracket=True` reproduces the slip; the default is the gradient); (6.108) (M = 2πρa³/3 by
the pressure integral and by energy); Example 6.2 (our Q = 1 by default; with the book's Q every printed value is reproduced); §6.7 "first-order"
differences are second order (order 2.000 measured).

## Open items
- **O1 (designer, check line; carried as design Part G1)** D21's sympy-check intent says `print(bad.mean())` → 1.0; it gives 0.998 because the
  three grid points on the slit (y = 0, |x| ≤ 2b) have |ζ| = b for both roots. Say "every left-half point off the
  slit" or exclude the slit in the cell.
- **O2 (designer, notebook text; Part G2)** C12/N74: the Example 6.2 refinement order is 4/3 **everywhere** (pollution by the
  corner singularity ψ ~ r^{2/3} sin(2θ/3), n = ⅔ in (6.46)), not "≈ 2 away from the corner"; the second-order claim
  holds for the smooth test problem (2.000, 1.995). The notebook's refinement cell should state 4/3 with this reason.
- **O3 (designer / viz-builder, E8 text and ranges; Part G3)** (a) "the segments converge to two narrow spikes" (Rankine oval):
  the k_n alternate in sign (± 9.6 at N = 20, ± 18.6 at N = 40); only the moments converge (dipole error 1.4e-3 →
  2.9e-6); (b) "the fitted bars reproduce a point source and a flat line sink" (airship): they do not (k from −4.9 to
  14.1 at N = 20, Σk_nΔξ = −0.094); (c) "panels: C_p converges like 1/N²" is true on the ellipse, not on the circle
  (exact at every N); (d) the N slider reaches 80, where the Rankine oval breaks down (cond 3.5e17, Σk_nΔξ = 0.13):
  cap N at 40 for the axial mode or show a conditioning warning; the ellipsoid breaks at N = 128.
- **O4 (resolved in loop 2: panel docstrings updated)** `panel_velocity` off the body converges at first order (0.93 → 1; the
  panel dipole moment has relative error ≈ 1.4/N) while the control-point C_p is exact (circle) or second order
  (ellipse). Our hypothesis: inherent to constant-strength panels on an inscribed polygon (the classical Hess–Smith
  circle case also over-predicts the strengths at coarse N); not asserted beyond the observed rate. `source_panels`'
  docstring still says "V3 circle C_p → 1 − 4 sin²θ at order 2" (it is exact) — update.
- **O5 (resolved in loop 2: no "Validation (planned)" left; scripted labels audited)** 47 docstrings still read "Validation (planned): …" (19 in the chapter module, 22 in
  `core.potential`, 2 each in `core.conformal`, `core.laplace_solvers`, `core.panels`); replace by the evidence now in
  `tests/test_ch06.py` (skill lesson: stale validation text).
- **O6 (qualitative, labelled)** Fig. 6.10's measured high-Re curve has no fetched dataset (candidates Roshko 1961,
  Achenbach 1968); N31's band stays *qualitative*. The sphere and airship targets of the axial method are fenced
  qualitative beyond the moments (conditioning).

## Addendum (2026-09-24) — odd N in the axial singularity method
The implementer added an odd-N branch to `core.panels.axial_singularity_solve`: for a fore–aft symmetric target the
reflection z → −z gives A = −PAP (P the flip), so for odd N the symmetric strength subspace has one dimension more than
the antisymmetric ψ subspace and A is exactly singular; the solver now warns (RuntimeWarning), returns the minimum-norm
least-squares (antisymmetric) k and sets `odd_symmetric = True`. New test
`test_axial_method_V1_fore_aft_symmetry_makes_odd_n_singular`:
(a) A = −PAP to 1e-13·max|A| for the Rankine oval and the ellipsoid at N = 10, and a symmetric k maps onto an odd ψ;
(b) Rankine oval N = 9, 11: the warning, `odd_symmetric` True, rank A = N − 1, k antisymmetric (1e-9 relative),
Σk_nΔξ < 1e-9, body error 4.8e-3 m and 1.5e-3 m (< 1e-2 m; 0.3–0.5 m before the branch);
(c) even N unchanged: Rankine N = 20 cond₁ = 3337914.99387 (1e-9) and body error 1.66106771e-4 m (1e-6), no warning,
`odd_symmetric` False; (d) the airship (not symmetric) at N = 9, 11 takes `np.linalg.solve` without a warning (residual
< 1e-8). Observation for the explainer (O3): the airship's body error at odd N is larger (0.28 m at N = 9, 0.013 m at
N = 11, vs 3.1e-3 m at N = 10) — conditioning, not the new branch. Runs after the addendum: `tests/test_ch06.py`
**137 passed in 71 s**; full suite **750 passed in 390 s**.

## Verdict: PASS
Loop 3 (review follow-up). All 12 scripts ran; `tests/test_ch06.py` **137/137 pass** after the odd-N addendum (133 + 4 skipped without
the private file); the full suite **750 passed**; every computable CORE row has ≥ 2 independent levels (one of
V1/V2/V3/V5), every coded NOTE row ≥ 1, all 100 Part C rows exercised, all 23 ★★/★★★ derivations re-derived
symbolically, all planted wrong variants caught (28 in loops 1–2, 2 more for S3/S5 in loop 3). Review items M1 and
S1–S8 are reflected in the tests and references; (6.82) is no longer called a slip (it is r × the App. B divergence).
O1–O3 remain carried in the design's Part G errata; O6 (the qualitative real-cylinder band, the fenced sphere/airship
axial targets) remains a labelled qualitative item. Nothing is `unverified`.
