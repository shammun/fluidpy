# Chapter 5 verification — Vorticity Dynamics                     2026-09-23, commit a37d713 + working tree (loop 2)

Verifier: math-verifier. Suite: `tests/test_ch05.py` (113 test functions = 133 collected items; the scenario and
exact-solution tests are parametrised), figures/metrics: `tests/ch05_verify_figures.py` → `outputs/ch05/verify/`
(git-ignored), cited data: `reference/ch05/` (`make_refs.py`, `benchmarks.json`, `SOURCES.md`). Nothing in `fluidpy/`
or `scripts/` was edited by the verifier.

## Environment
Python 3.11.5 · numpy 2.4.6 · scipy 1.17.1 · sympy 1.14.0 · pint 0.25.3 · matplotlib 3.11.2 (Windows, `.venv`).

Runs (loop 2, final): `pytest tests/test_ch05.py -q -p no:cacheprovider` → **133 passed in 125 s**, no warnings (with the
private book file; without it the four V6 tests skip). Full suite `pytest -q` → **611 passed in 329 s** (478 earlier + 133 ch05). All 9
`scripts/ch05_*.py` exit 0 (3–23 s each; also asserted by `test_scripts_V1_every_ch05_script_runs`).
`tools/check_public.py` → OK.

## Loop history
| loop | result | what changed |
|---|---|---|
| 1 | FAIL: 130/133 (full suite 608 passed, 3 failed) | F1 (`vortex_pressure_scenario` at r = a: half the wall torque, a step-dependent "net force" under a "no net viscous force" status) and F2 (Hill exterior u_R sign). O1–O6 opened. 40 wrong variants all caught. |
| 2 | **PASS: 133/133** | Implementer: F1 — at r = a the fluid-side outer profile Γ/2πr is differentiated (never a stencil across the kink), `net_viscous_force` = NaN there, new key `edge_line_force` = the stress jump, edge status texts ("Rankine core edge …", "cylinder wall …"); F2 — `uR_out = +1.5·U·a³·R·z/r⁵`; O1 — `biot_savart_2d` docstring states the interior accuracy. Designer (coordinator): O5 (D07 step 5 exact at the interface centre), O6 (D12 step 3: 30 % at 10a, 3 % at 100a). Verifier: the core-edge test now checks `edge_line_force` independently (below); F1 and F2 re-planted into a scratch copy — both caught. |

Collected items per evidence tag (tag of the `def` line; most tests carry more levels inside): **V1 86 · V2 23 · V3 10 ·
V4 8 · V5 2 · V6 4** (= 133).

## F1 (loop 1, resolved in loop 2) — `vortex_pressure_scenario` at r = a (Rankine core edge, rotating-cylinder wall): half the torque and a step-dependent "net force"
- Function: `ch05.vortex_pressure_scenario` (`fluidpy/ch05_vorticity_dynamics.py`, the block after
  `hh = 1e-5 * max(r_, 1e-9)`): at r = a the point is treated as outside (`rigid = r_ < a` is False) and
  `polar_net_viscous_force` differentiates u_θ with a central step h = 1e-4 r that **straddles the kink** of u_θ.
- Tests: `test_vortex_pressure_scenario_V1_core_edge_is_consistent[cylinder-kw1]`, `[rankine-kw0]`.
- Input / expected / got:
  - kind "cylinder", ω = 2 s⁻¹, a = 0.1 m, μ = 1e-3 Pa s, r = a (the cylinder's surface, where the torque is applied):
    expected σ_rθ = −μΓ/πa² = −2.000e-3 Pa and torque −2μΓ = −1.2566e-4 N m/m (`torque_per_length`, N08);
    **got σ_rθ = −9.9985e-4 Pa and torque −6.2822e-5 N m/m (half)**, net_viscous_force **−100.0 N/m³** with status
    "irrotational: viscous stress ≠ 0 but no net viscous force". At r = 1.001a the scenario is correct (−1.2566e-4).
  - kind "rankine", Γ = 1 m²/s, a = 0.1 m: torque −9.9985e-4 vs −2.000e-3; net force −1591.5 N/m³.
  - The E2 tornado preset evaluated at its core edge (ρ = 1.2, Γ = 1e4, a = 50 m) reports −0.127 N/m³.
- Evidence: σ_rθ at r = a is the average of the rigid value 0 and the outside value (a central difference across the
  jump of r d(u_θ/r)/dr); the force is that jump divided by 2h, so it grows as 1/h (`outputs/ch05/verify/core_edge_F1_verify.png`).
  Physically the net viscous force at the edge is concentrated (a jump of ω makes −μ∇×ω a line force); it has no finite
  pointwise value, and the stress and torque on the fluid side are the outside values.
- Suggested fix (implementer's choice): at r = a evaluate σ_rθ and the torque one-sided from the fluid side (the closed
  form −μΓ/πr² holds for r ≥ a), and give r = a the design's third status "Rankine core edge" (design Part B E2) with the
  force reported as concentrated/undefined rather than a number — or keep a finite number only with a status that says
  so. `radial_residual`, p, B and u_θ at r = a are correct.

## F2 (loop 1, resolved in loop 2) — `hill_spherical_vortex` outside the sphere: u_R has the wrong sign (fluid crosses the sphere)
- Function: `core.vortices._hill_parts` (`fluidpy/core/vortices.py`), line `uR_out = -1.5 * U * a_ ** 3 * R_ * z_ / r ** 5`;
  reached by `ch05.hill_spherical_vortex(..., outside=True)` (the default) and `hill_spherical_vortex_field` outside
  r = a. The stream function `hill_stream_function` (ψ_out = −(U/2)R²(1 − a³/r³)) and everything inside the sphere are
  correct.
- Test: `test_hill_vortex_V1_outside_velocity_is_the_stream_functions`.
- Input / expected / got: A = 1.3, a = 0.7, (R, z) = (0.8, 0.9): expected u_R = −R⁻¹∂ψ/∂z = +0.0124273 m/s (central
  differences of the coded ψ; closed form +(3/2)Ua³Rz/r⁵); **got −0.0124273**. u_z agrees (1e-7). On r = a⁺ (A = a = 1)
  the normal velocity u·e_r reaches **±0.154 m/s** instead of 0 (inside: < 4e-7), and the tangential velocity is not
  continuous — `outputs/ch05/verify/hill_outside_F2_verify.png`.
- Suggested fix: `uR_out = +1.5 * U * a_ ** 3 * R_ * z_ / r ** 5`. (Hill's vortex is N23, a C-depth NOTE "named", and
  the exterior is our addition — but the notebook draws the flow outside the sphere, so the field must be right.)

## Loop-2 re-check of F1 and F2 (on our own terms)
- **F1:** `edge_line_force` is checked against an independent computation in
  `test_vortex_pressure_scenario_V1_core_edge_is_consistent`. On the core side u_θ/r = Γ/2πa² is constant, so
  σ_rθ(a⁻) = 0; on the fluid side μ r d(Γ/2πr²)/dr = −μΓ/πa², so the jump σ_rθ(a⁺) − σ_rθ(a⁻) = −μΓ/πa². The code gives
  Rankine (Γ = 1, a = 0.1, μ = 1e-3) **−0.0318310 Pa** (σ, jump), torque **−2.000e-3 N m/m**; cylinder (ω = 2, a = 0.1)
  σ = jump = **−2.0000e-3 Pa**, torque **−1.25664e-4 N m/m** (= −2μΓ, 1e-7). Also asserted: the edge force per unit
  length balances the torque carried outward, 2πa·jump·a = −2μΓ (1e-7); `net_viscous_force` is NaN at r = a and the
  status no longer claims "no net viscous force"; at 0.999a σ = jump = 0, at 1.001a the jump is 0 and the force ≈ 0; the
  line vortex reports no jump. The tornado preset's edge now reports jump −1.2732e-3 Pa and NaN force.
- **F2:** (0.8, 0.9) gives u_R = +0.0124273 = −R⁻¹∂ψ/∂z; the normal velocity 1e-6a outside the sphere is ≤ 4.0e-7
  (the O(Δr) offset of the probe; the tangential velocity is continuous to 1e-4) — `hill_outside_F2_verify.png` now
  shows arrows along the ψ contours.
- **Re-planted in a scratch copy** (the repository untouched): the loop-1 kink-straddling stencil → 2 tests fail (both
  core-edge cases); `edge_line_force` with the wrong sign → 2; a finite net force kept at r = a → 2; the loop-1 Hill
  sign → 1 (`test_hill_vortex_V1_outside_velocity_is_the_stream_functions`).

## Discrimination proofs (wrong variants planted in a scratch copy; the repository untouched)
Forty wrong variants were patched one at a time into a scratch copy of `fluidpy/` (scratchpad `discriminate_ch05.py`,
the ch05 suite without the slow script-smoke test; baseline = the 3 F1/F2 failures) — **every one of the 40 fails at least one test** (median 2, max 10). The first run caught 39; the planetary-transpose variant
slipped through the rotating-column test (whose velocity gradient is symmetric in its z row and column), so a planetary-tilting
check was added and the variant re-run (fails 1).

| wrong variant patched in | tests that fail (beyond the baseline) |
|---|---|
| (5.14) with the printed −1/(4π) as the default sign | 1 — test_velocity_from_curl_omega_V1_corrected_sign_and_printed_sign |
| ellipk/ellipe given the modulus k instead of the parameter m = k² | 2 — test_ring_dynamics_V4_leapfrog_and_ring_toward_a_wall, test_ring_ring_velocity_V3_elliptic_form_equals_polygon |
| sheet strength with the caption's sign as the default | 2 — test_book_V6_printed_forms_of_section_5_1_and_5_7, test_vortex_sheet_strength_V1_text_convention_and_caption_variant |
| vorticity budget without the tilting (off-diagonal) part of (ω·∇)u | 1 — test_vorticity_terms_V1_exact_solutions_balance[abc] |
| stretching/tilting split with the tilting dropped | 2 — test_stretching_tilting_scenario_V1_statuses_and_numbers, test_stretching_tilting_split_V1_decomposition_V7_invariance |
| (5.1) with ω passed as ch03's rotation rate | 5 — test_bernoulli_across_vortex_V1_rotational_vs_irrotational, test_solid_body_V1_omega_is_the_vorticity_not_the_rate, test_solid_body_pressure_V1_gradients_isobars_contract, test_vortex_pressure_scenario_V1_full_balance_every_kind[solid-kw0] … |
| baroclinic term ∇p × ∇ρ (order reversed) | 5 — test_baroclinic_term_V1_barotropic_zero_V3_order_V2_symbolic, test_lock_exchange_V1_field_route_sense_and_contract, test_pressure_torque_V1_contract_barotropic_and_callables, test_vorticity_budget_V1_rotating_column_and_absolute_vorticity … |
| lock-exchange rate without the factor 2 | 4 — test_baroclinic_V2_dimensions, test_book_V6_exercise_closed_forms, test_lock_exchange_V1_field_route_sense_and_contract, test_lock_exchange_V2_derivation |
| net viscous force without the polar metric r² (naive dσ/dr) | 6 — test_stress_vs_net_force_table_V1_three_rows, test_viscous_force_V1_naive_divergence_misses_the_metric, test_viscous_stress_V2_derivation, test_vortex_pressure_scenario_V1_full_balance_every_kind[cylinder-kw3] … |
| Kelvin pressure term with a constant (mean) density | 1 — test_kelvin_force_terms_V1_barotropic_zero_baroclinic_stokes |
| line-vortex pressure with 4π² instead of 8π² | 9 — test_bernoulli_across_vortex_V1_rotational_vs_irrotational, test_book_V6_printed_forms_of_section_5_1_and_5_7, test_line_vortex_pressure_V1_radial_balance_and_funnel, test_solid_body_pressure_V1_gradients_isobars_contract … |
| Biot–Savart kernel (x − x′) × ω (cross product reversed) | 2 — test_biot_savart_volume_V3_spectral_convergence_and_long_tube_limit, test_velocity_from_curl_omega_V1_corrected_sign_and_printed_sign |
| point-vortex velocity r × e_z (clockwise-positive) | 10 — test_book_V6_exercise_closed_forms, test_circle_images_V1_circle_is_a_streamline_and_knife_pair, test_cylinder_as_vortex_sheet_V1_ring_of_point_vortices, test_discrete_sheet_V3_first_order_in_N_and_continuous_form … |
| wall image with the same sign | 2 — test_book_V6_exercise_closed_forms, test_wall_image_V1_no_penetration_and_drift |
| fluid column (ζ + f)h = const (inverted) | 2 — test_fluid_column_V1_ridge_trough_and_ring_of_air, test_fluid_column_V2_derivation |
| absolute circulation Γ − 2Ω·A | 1 — test_absolute_circulation_V4_material_loop_in_rotating_frame |
| pressure torque taken about the geometric centre O instead of G | 5 — test_baroclinic_V2_dimensions, test_baroclinic_element_scenario_V1_closed_form_equals_numeric, test_pressure_torque_V1_contract_barotropic_and_callables, test_pressure_torque_V2_derivation … |
| helix e_n toward the centre of curvature (Frenet N) | 2 — test_frenet_frame_V1_numeric_frame_on_the_helix, test_natural_coordinates_V2_derivation |
| Kelvin ring speed with the hollow-core constant for 'uniform' | 1 — test_ring_self_velocity_V5_kelvin_thin_ring |
| vortex near a wall at Γ/2πh (image distance h) | 2 — test_book_V6_printed_forms_of_section_5_1_and_5_7, test_wall_image_V1_no_penetration_and_drift |
| circle image at the projection on the circle, not the inverse point | 1 — test_circle_images_V1_circle_is_a_streamline_and_knife_pair |
| Burgers vortex with Wikipedia's α in the book's formula (a factor 2) | 2 — test_book_V6_exercise_closed_forms, test_burgers_vortex_V1_wikipedia_form_pressure_and_ns |
| uniform-strain vorticity e^{Gᵀt}ω₀ | 2 — test_stretching_tilting_scenario_V1_statuses_and_numbers, test_uniform_strain_V1_presets_V2_expm_solves_the_equation |
| segment law with the end terms swapped (sign) | 4 — test_filament_V3_polygon_ring_converges_order_two_and_contracts, test_ring_ring_velocity_V3_elliptic_form_equals_polygon, test_segment_V1_infinite_line_recovers_5_2_and_wikipedia_form, test_segment_induced_velocity_V1_equals_quadrature_of_5_17 |
| Kelvin contour term du = Gᵀdx | 2 — test_kelvin_rate_V3_total_equals_dGamma_dt_of_material_loop, test_kelvin_scenario_rate_V1_every_scenario_full_split[rotating-3.0] |
| loop vector area without the ½ | 3 — test_absolute_circulation_V2_derivation, test_absolute_circulation_V4_material_loop_in_rotating_frame, test_loop_helpers_V1_spectral_quadrature_on_closed_loops |
| diffusing sheet u = +(γ/2) erf (sign) | 1 — test_diffusing_vortex_sheet_V2_symbolic_V1_invariants |
| Rankine pressure: core constant not matched at r = a | 3 — test_bernoulli_across_vortex_V1_rotational_vs_irrotational, test_book_V6_exercise_5_2_tornado, test_rankine_pressure_V1_matched_core_and_tornado_contract |
| planetary term 2Gᵀ·Ω in (5.30) | 1 — test_vorticity_budget_V1_rotating_column_and_absolute_vorticity (**0 in the first run**: the rotating column's G has equal z-row and z-column; the planetary-tilting check u = (sz, 0, 0) and a generic-field comparison were added, then the variant fails) |
| tube budget: lower end counted with the inward normal | 2 — test_tube_flux_budget_V1_broken_field_fails_exactly_as_computed, test_tube_flux_budget_V4_gauss_on_tube_pieces |
| cylinder wall power without the factor 2 of the torque | 2 — test_rotating_cylinder_V1_no_slip_torque_and_energy, test_vortex_pressure_V2_dimensions |
| E3 Coriolis loop term with Ω instead of 2Ω | 2 — test_kelvin_scenario_rate_V1_every_scenario_full_split[rotating-0.0], test_kelvin_scenario_rate_V1_every_scenario_full_split[rotating-3.0] |
| frozen-in element D(δx)/Dt = Gᵀδx | 1 — test_uniform_strain_V1_presets_V2_expm_solves_the_equation |
| sheet roll-up kernel u with the wrong sign | 1 — test_sheet_rollup_V1_periodic_kernel_V4_centroid |
| channel image series: a fixed 10 000 terms, no tail | 2 — test_book_V6_exercise_closed_forms, test_channel_images_V1_series_closed_form_and_rhs |
| rotating NS residual with Ω × u (factor 2 lost) | 2 — test_rotating_lamb_form_V1_equals_5_20_residual_on_random_fields, test_rotating_ns_V1_corotating_rest_and_inertial_oscillation |
| (5.25) Lamb term with ω + Ω | 1 — test_rotating_lamb_form_V1_equals_5_20_residual_on_random_fields |
| Gaussian tube ω_R with the wrong sign (not solenoidal) | 5 — test_gaussian_tube_V1_cartesian_field_and_sections, test_gaussian_tube_V2_flux_function_construction_is_solenoidal, test_tube_flux_budget_V4_gauss_on_tube_pieces, test_tube_flux_budget_traced_V3_closed_polyhedron_converges … |
| centre of vorticity nearer the weaker vortex | 1 — test_point_vortex_evolve_V4_invariants_and_V1_pair_orbits |
| Lamb–Oseen circulation rate with the wrong sign | 3 — test_kelvin_scenario_rate_V1_every_scenario_full_split[lamb_oseen-0.0], test_kelvin_scenario_rate_V1_every_scenario_full_split[lamb_oseen-20.0], test_material_circulation_V1_lamb_oseen_decays_as_closed_form |

## Validation table — A items (CORE, ≥ 2 independent levels, one of V1/V2/V3/V5)
| Concept / Eq. | Tier | fluidpy target | Evidence (V-levels) | Numbers (error, order, residual) | Label | Notes |
|---|---|---|---|---|---|---|
| C01 tubes cannot end (5.4); R01, N01–N03, N45; D01 | CORE | `vorticity_field`, `vortex_line`, `vortex_tube_strength`, `Tube`, `vortex_tube`, `tube_flux_budget(_traced)`, `gaussian_tube_*`, `broken_tube_field`, `vortex_ring_vorticity`, `tube_core_radius` | V1 u = b×x → 2b (1e-10); ABC vortex lines = streamlines (1e-6 m) and tangent to ω; helical swirl lines keep zR² (ptp < 1e-8); Lamb–Oseen tube strength both routes = Γ(1 − e⁻¹) (circle 1e-12, disc route order 2, 2e-6 at nr = 256), orientation flip; V2 sympy: the flux-function tube is solenoidal, its walls R/a(z) = const are flux surfaces, code = formulas (1e-12); V4 Gauss budget: total ≤ 2e-12, side ≤ 2e-12, four tube pieces; "broken" total = flux₀(e⁻¹ − 1) = −0.399577 (1e-9, fails as it must); traced lines stay on the wall (ptp 2e-10); V3 traced polyhedron end flux → Γ(1 − e⁻¹) at order 2 in 1/n_lines, closed total ≤ 5e-9; ring meridional flux = Γ (1e-9) | curl order 2.000; traced 2.0 | analytic, symbolic, conserved, converged | contract numbers 0.632121, 20.121, 54.695, 86.526 asserted |
| C02 pressure of the two basic vortices (5.5)–(5.7); R02–R08, N04–N09, N46, N47; D02, D03 | CORE | `solid_body_from_vorticity`, `line_vortex_gamma`, `solid_body_pressure(_gradients)`, `line_vortex_pressure`, `rankine_pressure`, `isobar_height`, `rotating_tank_free_surface`, `bernoulli_across_vortex`, `polar_viscous_stress`, `line_vortex_viscous_stress`, `polar_net_viscous_force`, `vortex_stress_force`, `stress_vs_net_force_table`, `rotating_cylinder_flow`, `torque_per_length`, `dissipation_outside_cylinder`, `vortex_pressure_scenario`, `tornado_circulation`, `rankine_radius_at_pressure` | V1 curl of (5.1) = ω (WV: ω as rate → 2ω); ∮ = Γ on circles; (5.5a) FD balance 1e-7/1e-6; isobars; tank volume by independent quad (1e-9, dry bottom, lid); Rankine = ch04 `rankine_vortex_pressure` (1e-12), continuous p and slope; tornado 31.831 m/s, −607.93, −1215.85 Pa; B = ω²r²/4, 0, −1 (ch04); σ_rθ, zero force three routes, naive dσ/dr ≠ 0; torque −2μΓ at 7 radii (1e-13); dissipation by our polar quad = power in − out (1e-10), → power in as R_out → ∞; E2 full balance for 4 kinds × 5 radii (radial residual ≤ 1e-6, σ, force, torque, B, surface) and at r = a (edge jump −μΓ/πa², torque −2μΓ, NaN force); V2 D02/D03 sympy (below); pint (5.6), (5.7), σ_rθ, torque, power | — | analytic, symbolic, conserved | F1 fixed (loop 2): edge jump checked |
| C03 Kelvin (5.8)–(5.11); N10–N14, N16, R09, N48; D04, D05 | CORE | `material_loop`, `material_circulation`, `loop_circulation`, `loop_length`, `loop_tangent`, `loop_line_integral`, `circle/square_loop_points`, `kelvin_rate_terms`, `kelvin_force_terms`, `kelvin_scenario(_circulation, _rate, _gamma)`, `lamb_oseen_circulation`, `lamb_oseen_viscous_loop_integral`, `kelvin_hypotheses(_text)`, `cellular_flow` | V4 cellular loop stretched ×6.68 over 6 turnovers, Γ = 1.155130 constant to 8.2e-13 (flux route by dblquad 1e-11; points on pathlines 1e-8; ψ conserved); Rankine straddle = 2 × lens area (2e-5, kink), Gaussian 1e-9; V1 Lamb–Oseen circle = Γ₀(1 − e^{−r²/4ν(t+t₀)}) (1e-11), viscous loop integral = ∂Γ/∂t (1e-6); V3 (5.9) total = central dΓ/dt of a material loop in a manufactured unsteady flow, order 1.996; contour term ≤ 1e-8 of the acceleration term; −∮dp/ρ = ∫(∇ρ×∇p/ρ²)dA by dblquad (1e-8), ρ(p) → 0, ρ₀ constant → 0 (WV), body terms; all 16 hypothesis combinations; every E3 scenario's full split at 11 (scenario, t) pairs (acceleration = pressure + viscous + Coriolis, each term as predicted); V2 D04/D05 | order 1.996 | conserved, analytic, converged, symbolic | |
| C04 baroclinic pressure torque (Fig. 5.6), lock exchange; N15, N49, N50; D06, D07 | CORE | `pressure_torque_on_element`, `baroclinic_term(_sym)`, `baroclinic_rate_2d`, `lock_exchange_initial_vorticity_rate`, `lock_exchange_fields`, `baroclinic_element_scenario` | V2 D06 exact disc integrals: F = −πR²∇p, M_O = 0, x_G = R²∇ρ/4ρ₀, M_G = πR⁴(∇ρ×∇p)/4ρ₀, I_G exact, spin = (∇ρ×∇p)/ρ₀²/(1 − ∣∇ρ∣²R²/8ρ₀²); code = exact (1e-9), ratio 1 + 1.25e-9 (1e-11); V3 ratio − 1 order 2.000 and equal to the exact factor (1e-6); V1 barotropic → 0, callables = presets (1e-9), lock field route 2.422222 (1e-7), sense (WV ∇p×∇ρ), interface-integrated rate ρ̄gΔρ/(ρ₁ρ₂) independent of δ; `baroclinic_term` = sympy (order 2.00), ρ(p) → 0; E4 closed form = quadrature; pint | order 2.000 | symbolic, analytic, converged | |
| C05 Helmholtz's theorems; N17, N51; D08 | CORE | `frozen_in_check`, `kelvin_scenario("helmholtz", "helmholtz_abc")`, `abc_flow`, `stretched_gaussian_vortex_field` | V1 ABC ×2 and the inviscid stretched Gaussian: angle(δx, ω) < 1e-8, ∣ω∣/∣δx∣ − 1 < 1e-8 (integrated) and < 1e-7 (field ω); on the axis ∣δx∣ grows as e^{αt}; viscous Burgers: ratio → < 0.5 (WV); V4 patch on a tube wall: flux stays < 1e-10 (Gaussian), ABC patch Γ = 8.016e-10 constant (1e-6 rel) with normal ⟂ ω (≤ 2e-5), a patch facing ω carries ∣ω∣A; V2 D08: Dχ/Dt = 0 for the inviscid tube walls, ≠ 0 with ν; ABC ω = u, Lamb vector 0 | — | analytic, conserved, symbolic | |
| C06 vorticity equation (5.12)–(5.13); N18–N20, N22, N23, R10, R11; D09 ★★★ | CORE | `vorticity_equation_sym`, `vorticity_terms`, `vorticity_terms_sym`, `vorticity_budget_preset`, `diffusing_vortex_sheet`, `hill_spherical_vortex(_field)`, `burgers_vortex(_field)`, `lamb_oseen_field`, `planar_field_3d` | V2 D09 re-run for a generic u = ∇×A (every step; (B.3.10) for unrelated a, b; curl(NS residual) ≡ (5.13) residual; the ω∇·u term when compressible); function outputs on generic and test fields (Hill residual 0); V1 residual ≤ 1e-5 of the terms on Lamb–Oseen, Taylor–Green, Burgers, Hill, ABC and the viscous stretched Gaussian (every term nonzero); term relations per flow; V3 residual order 1.999 (Burgers) and 2.0 (Lamb–Oseen, fully nested); WV wrong ν → residual O(terms); presets contract + burgers_balance closed forms; sheet: sympy diffusion eq, ∫ω = γ, u(∞) = −γ/2, FTCS order 2.0; Hill: Wikipedia form (sign-mapped), ψ inside = code; Burgers: Wikipedia form (α_w = α/2), NS residual with `burgers_pressure` ≤ 1e-5, radial balance | order 1.999 | symbolic, analytic, converged | F2 fixed (loop 2) |
| C07 Biot–Savart (5.14)–(5.16); N24–N28, N52; D10 ★★★, D11 ★★★ | CORE | `poisson_green_3d`, `velocity_from_curl_omega`, `gaussian_tube_fields`, `cylinder_quadrature`, `biot_savart_volume`, `biot_savart_2d`, `velocity_from_vorticity_fft`, `curl_theorem_box` | V1 Gaussian tube (L = 4 m) at r = 0.5 m: +1/(4π) route = 0.308838 = 1-D reference (5e-6, three directions), printed −1/(4π) → −0.308838 (WV), (5.14) route = (5.16) route (1e-9); V3 quadrature error 4.1e-3 → 1.2e-6 (faster than any power); long tube → Γ(1 − e^{−r²/σ²})/2πr at (0.5/L)² (V7); result solenoidal (stencil order 2); FFT: Taylor–Green 2-D and ABC 3-D to 1e-13, periodic images fade with the box; plane kernel: outside spectral (< 1e-8), dual-grid interior order > 3.5, = from-scratch double loop (1e-12); (5.15) exact for polynomials (1e-12), b×x → 2b, midpoint order 2.0; V2 D10, D11 | orders 2.0 | analytic, converged, symbolic | O1 resolved (docstring) |
| C08 filament law (5.17); N29; D12, D13 | CORE | `segment_induced_velocity`, `segment_speed`, `filament_velocity`, `filament_contributions`, `filament_preset`, `filament_velocity_preset`, `ring_axis_velocity`, `ring_ring_velocity`, `ring_self_velocity`, `RING_CORES` | V1 closed form = quad of (5.17) at 6 random geometries (1e-9); long segment → Γ/2πd (1e-10), semi-infinite half, Wikipedia (cos A − cos B) form, 0.1125395; square 2√2/π (1e-12); V3 polygon ring order 2.007 (1.054786, 1.013052, 1.000804), on-axis off-plane order 2; elliptic ring = polygon at order 1.997, axis limit (1e-12), symmetry; WV modulus k for m misses by > 5 %; V5 Kelvin ring speed (3 cases, 1e-14), 0.328816; V2 D12 (ring axis ΓR²/2(R²+z²)^{3/2}), D13 | orders 2.007, 1.997 | analytic, converged, benchmark, symbolic | |
| C09 rotating/baroclinic vorticity equation (5.18)–(5.30); R12–R17, N30–N37; D14, D15 ★★★ | CORE | `vorticity_budget`, `vorticity_budget_sym`, `vorticity_budget_preset`, `rotating_lamb_form_terms`, `rotating_ns_residual`, `vorticity_divergence`, `absolute_vorticity`, `planetary_vorticity_terms`, `baroclinic_term` | V2 D14 ((5.21)–(5.24) for generic u; (5.25) ≡ (5.20) for u = ∇×A) and D15 with `sp.LeviCivita` (Lamb term = the full expression incl. the dropped u_{j,j}(ω_n + 2Ω_n); gradient term 0; pressure = (5.28); viscous = νω_{n,jj}; curl(5.25) ≡ (5.30) on a divergence-free polynomial field; `vorticity_budget_sym` residual = that curl, reduces_to_513); V1 co-rotating rest and inertial oscillation residual < 1e-8 (WV: no centrifugal, Coriolis ×1); (5.25) = (5.20) residual numerically on random fields (1e-5 of the terms) and = sympy residual; rotating column: local = (ζ + 2Ω)α, planetary 2Ωα, residual < 1e-7 (WV relative only); ∇·ω round-off; 2GΩ = 2Ω∂u/∂z, 1.458423e-9; absolute vorticity round trip | — | symbolic, analytic | |
| C10 stretching and tilting (5.31)–(5.32); N21, N38, N53; D16–D18 | CORE | `stretching_tilting_split`, `strain_preset`, `uniform_strain_vorticity`, `stretched_tube`, `stretching_tilting_scenario`, `helix`, `helix_frame`, `frenet_frame`, `helical_vortex_line`, `burgers_balance`, `burgers_core_radius` | V1 split = Gω, ∥/⟂ (10 random, 1e-12), 2-D → 0; V7 rotation invariance; presets e² = 7.389056, s t ω_x = 2, planar 1; V2 expm solves Dω/Dt = Gω for every preset and G² symmetric (steady Euler); frozen-in element = e^{Gt}ω₀ (1e-8, independent integration); tube ωA, AL conserved; statuses at the right instants (shear tilt: tilting at t = 0, rate 0.4/0.2 at t = 2); V2 D16 (frame, projections), helix κ = a/(a²+c²), τ = c/(a²+c²), book e_n = −N (1e-14); numeric Frenet frame; V2 D18 (ω ∝ L; Burgers ODE, total derivative, ∫ω dA = Γ, core √(4ν/α)) | — | analytic, symbolic | |
| C11 absolute circulation (5.33), fluid column; R18, N39, N54; D19, D20 | CORE | `absolute_circulation`, `loop_vector_area`, `column_relative_vorticity`, `column_over_slope`, `relative_circulation_after_move`, `coriolis_parameter` | V4 material loop advected in the rotating scenario: Γ(t) = 2ΩA₀(1 − e^{−αt}) (1e-9), Γ_a = π constant (≤ 1e-9); fluid at rest in the inertial frame: Γ = −2ΩA, Γ_a = 0; tilted loop 2Ω·A_vec; V2 D19 (triple product; DA_vec/Dt = ∮u×dx for a generic material loop; ∮(Ω×x)·dx = 2Ω·A_vec on a non-planar Fourier loop, = `loop_vector_area`), D20 ((ζ+f)/h, Dζ/Dt = (ζ+f)/h Dh/Dt); V1 ridge anticyclonic, trough cyclonic, ratio constant (1e-12); f(30°, 45°, 60°); ring of air −4.19261e7 m²/s, −5.33820e-5 s⁻¹, Γ + fA conserved | — | conserved, analytic, symbolic | |
| C12 point vortices; N40, N41, N55, N56; D21 | CORE | `point_vortex_velocity`, `point_vortex_rhs`, `point_vortex_evolve`, `point_vortex_invariants`, `centre_of_vorticity`, `point_vortex_preset`, `vortex_pair` | V1 single vortex Γ/2πr counterclockwise, self excluded; = from-scratch double loop (1e-13); V4 three vortices: P, I, H flat to 1e-10; V1 equal pair back after 10 periods (1e-8), quarter turn; unequal pair G fixed at 0.25 m (1e-10), rate 4/2π (1e-9); opposite pair +y at Γ/2πh; caption slip: fluid at G moves at −1.70 m/s; V5 García & Haziot (2023) Ω₀ = 1/(4πl²) and V₀ = −i/(4πl) (1e-14 / 1e-9) | — | analytic, conserved, benchmark | |
| C13 images, rings; N42, N43, N57–N59; D22 | CORE | `wall_image_system`, `circle_image_system`, `image_vortices`, `channel_image_velocity(_exact)`, `vortex_near_wall_speed`, `ring_dynamics`, `ring_self_velocity` | V1 wall: v = 0 at 200 points (1e-15), wall speed Γh/π(x²+h²) (1e-12), drift Γ/4πh by integration (1e-6) (WV same-sign image); circle inside/outside: normal velocity 0 at 400 points (1e-13), cylinder circulation 0.7 recovered (1e-9); knife pair: −x first, stays inside, separates ×3; channel series = (Γ/4H)cot(πh/H) at 4 heights (1e-9), 0 at H/2, → Γ/4πh at the wall, = ±20 000-pair direct sum (1e-4), coth kernel in the rhs (1e-12) and its drift; V4 leap-frog impulse ΣΓπR² drift 1.7e-9, ≥ 4 pass-throughs, a²R conserved; single ring at Kelvin speed; ring toward a wall: R monotone up, approach speed monotone down | — | analytic, conserved | O3 (thin-core ring near a wall) |
| C14 vortex sheet; N44, N60; D23 | CORE | `vortex_sheet_strength`, `vortex_sheet_velocity`, `continuous_sheet_velocity`, `discrete_sheet_u`, `discrete_sheet_convergence`, `sheet_rollup` | V1 text convention γ = u_below − u_above (WV caption −2); a counterclockwise row gives u_above = −γ/2 = continuous-sheet value (1e-6); circulation of a thin box = γ ds (5e-3); continuous sheet closed form = quad of the kernel (1e-10); V3 L1 error 0.043968, 0.0044131, 0.00044128, order 0.999; roll-up kernel = the cot form (1e-6) and the ±3000-period sum (O(1/M)); V4 mean y conserved (1e-15); Exercise 5.19 geometry: a sheet ring gives Ωa²/r outside (1e-9), rest inside | order 0.999 | analytic, converged, conserved | roll-up shape qualitative (O4) |

Every computable A item has ≥ 2 independent levels with at least one of V1/V2/V3/V5.

## Validation table — coded B/C items (NOTE, ≥ 1 level) and recaps
| Concept | fluidpy target | Evidence (test) | Label |
|---|---|---|---|
| N02 (5.3) vortex lines | `vortex_line` | V1 ABC and helical-swirl tests | analytic |
| N03 tube strength | `vortex_tube_strength` | V1 two routes, V3 disc order | analytic, converged |
| N05, R06 (5.5a, b) | `solid_body_pressure_gradients` | V1 FD balance, V2 D02 | analytic, symbolic |
| N06, N09 σ_rθ, stress vs force | `polar_viscous_stress`, `vortex_stress_force`, `stress_vs_net_force_table`, `polar_net_viscous_force` | V2 D03, V1 three routes + naive WV | symbolic, analytic |
| N07 (5.7), Rankine tornado | `line_vortex_pressure`, `rankine_pressure`, `tornado_circulation`, `rankine_radius_at_pressure` | V1, V2, V6 | analytic, book-value |
| N08 torque, dissipation | `torque_per_length`, `dissipation_outside_cylinder` | V1 / V4 energy balance | analytic, conserved |
| N10–N13 (5.9)–(5.11) | `kelvin_rate_terms`, `kelvin_force_terms`, `lamb_oseen_circulation`, `lamb_oseen_viscous_loop_integral` | V3, V1, V2 | converged, analytic, symbolic |
| N14, N16 three sources, four restrictions | `kelvin_hypotheses(_text)`, `kelvin_scenario_rate` | V1 16 combinations; per-scenario split | analytic |
| N15 lock exchange | `lock_exchange_*` | V1, V2 D07, V6 | analytic, symbolic, book-value |
| N17 Helmholtz proof | `kelvin_scenario("helmholtz", "helmholtz_abc")` | V4 | conserved |
| N19–N20 (5.12), (B.3.10) | `vorticity_equation_sym` | V2 | symbolic |
| N21 Burgers | `burgers_vortex(_field)`, `burgers_core_radius`, `burgers_balance`, `VX.burgers_pressure` | V1 (form, NS residual), V2 D18, V6 | analytic, symbolic, book-value |
| N22 diffusing sheet | `diffusing_vortex_sheet` | V2, V1, V3 FTCS | symbolic, analytic, converged |
| N23 Hill | `hill_spherical_vortex(_field)`, `VX.hill_stream_function`, `VX.hill_translation_speed` | V1 form, exterior = ψ's velocity (F2 fixed), V2 (5.13) residual 0, V6 | analytic, symbolic |
| N24 FFT Poisson | `velocity_from_vorticity_fft` | V1, V3 | analytic, converged |
| N25–N28 (5.14)–(5.15), V′ | `velocity_from_curl_omega`, `poisson_green_3d`, `curl_theorem_box` | V1, V2 D10–D11, V3 | analytic, symbolic, converged |
| N29 infinite line | `segment_induced_velocity`, `segment_speed` | V1, V2 D13 | analytic, symbolic |
| N30 (5.18) | `vorticity_divergence` | V1 round-off, V2 identity | analytic, symbolic |
| N31–N36 (5.24)–(5.29) | sympy in D14/D15; `rotating_lamb_form_terms` | V2, V1 | symbolic, analytic |
| N37 reading (5.30) | `vorticity_budget_preset` | V1 six scenes, residual ≤ 1e-6 | analytic |
| N38 (5.31) | `stretching_tilting_split` | V1, V7, V2 D16 | analytic, symbolic |
| N39 planetary terms | `planetary_vorticity_terms` | V1 | analytic |
| N40–N41 pairs | `vortex_pair`, `point_vortex_evolve` | V1, V5, V6 | analytic, benchmark, book-value |
| N42–N43 rings | `ring_dynamics`, `ring_ring_velocity`, `ring_self_velocity` | V4, V3, V5 | conserved, converged, benchmark |
| N44, N60 sheet | sheet functions | V1, V3 | analytic, converged |
| N45–N60 figures | `tests/ch05_verify_figures.py` | visual verdicts below | — |
| R01–R18 recaps | reused ch01–ch04 functions | called and asserted where ch05 relies on them (`K.streamline`, `K.pathline`, `K.acceleration`, `K.vorticity_in_rotating_frame`, `NS.exact_solution("lamb_oseen")`, `NS.ns_incompressible_terms`, `core.bernoulli.rankine_vortex_pressure`, `core.diffusion.ftcs_diffusion_1d`, `core.curvilinear` operators, `coriolis_parameter`) | analytic |

Conceptual NOTE items without computable output (N01, N04, N18, N51) are words only (curation §2); SKIP S01, S02 not
coded.

## Derivations (curation §4b, design Part F) — every ★★ and ★★★ re-derived with sympy
| D | ★ | test | what is re-derived (intermediate lines checked) | finding |
|---|---|---|---|---|
| D01 | ★ | (tube tests) | numeric Gauss budget, sympy solenoidal tube | ✓ |
| D02 | ★★ | `test_solid_body_pressure_V2_derivation` | cylindrical Euler via `core.curvilinear` → −u_θ²/r e_r (steps 3–4); ∂p/∂r = ρω²r/4; partial integration + `dsolve` for f(z) (steps 8–9); (5.6); isobar (step 10); trap ×4; (5.7) by the same moves; B grows as ω²r²/4, uniform for (5.2); code = derived forms | ✓ |
| D03 | ★★ | `test_viscous_stress_V2_derivation` | steps 1–9: σ_rθ = −μΓ/πr², σ_rr = σ_θθ = 0, metric divergence 0, naive ≠ 0, Laplacian form 0, `CU.vector_laplacian` 0; Gaussian contrast μ∂ω/∂r ≠ 0 | ✓ |
| D04 | ★★ | `test_kelvin_rate_V2_derivation` | generic labelled loop X(s, t) moving with a generic unsteady field: d/dt of the integrand = Du/Dt·X_s + ∂_s(½∣u∣²) (steps 3–7); ∮dθ = 2π trap | ✓ |
| D05 | ★★ | `test_kelvin_V2_derivation` | ρ(p) ⇒ curl(∇p/ρ) = 0; baroclinic curl = (∇ρ×∇p)_z/ρ²; conservative Φ; Lamb–Oseen ∂Γ/∂t = 2πrν(∇²u)_θ; the ∇(∇·u) caveat | ✓ |
| D06 | ★★ | `test_pressure_torque_V2_derivation` | all 10 steps with exact disc integrals (F, M_O, M, x_G, M_G, I_G, spin, limit) | ✓ (step 8's "≃" is exact up to the stated factor) |
| D07 | ★★ | `test_lock_exchange_V2_derivation` | tanh step + hydrostatic mean-density p: rate at x = 0 = 2Δρg/((ρ₁+ρ₂)δ) exactly; slope; sense; δ-independent integral ρ̄gΔρ/(ρ₁ρ₂) | ✓ (step 5's "ρ² ≈ ρ̄²" is exact at the interface centre — note O5) |
| D08 | ★★ | `test_helmholtz_V2_derivation` | tube walls χ = const material (inviscid), not with ν; ABC Beltrami and Euler | ✓ |
| D09 | ★★★ | `test_vorticity_equation_V2_derivation` | the design's cell step by step on u = ∇×A (steps 2–12), (B.3.10) for generic a, b, the compressible ω∇·u difference | ✓ |
| D10 | ★★★ | `test_poisson_green_V2_derivation` | curl of curl; ∇²(1/r) = 0; flux −4π; G unit source; constructive check ∇²φ = +q for φ = ∫Gq (radial source); the printed −1/(4π) solves ∇²u = +∇×ω | ✓ (book sign slip confirmed) |
| D11 | ★★★ | `test_biot_savart_V2_derivation` | product rule; ∇′(1/r) = +(x − x′)/r³; correct step 5; the book's step-5 variant differs by 2ω×(x − x′)/r³; each curl component is a divergence (Gauss per component) | ✓ (book's compensating slip confirmed) |
| D12 | ★★ | `test_filament_V2_derivation` | (5.17) round a circle: ΓR²/2(R²+z²)^{3/2} e_z, centre Γ/2R; frozen-kernel change 3a/D | ✓ |
| D13 | ★★ | `test_segment_V2_derivation` | e_z × (d e_d − l e_z) = d e_φ; ∫ = (cos θ_a − cos θ_b)/d; substitution l = −d cot θ (0 < θ < π); limits 2/d and 1/d | ✓ |
| D14 | ★★ | `test_rotating_lamb_form_V2_derivation` | (5.21), (5.22) generic; (5.23) with the ν∇(∇·u) remainder; (5.24); (5.25) ≡ (5.20) for u = ∇×A | ✓ |
| D15 | ★★★ | `test_rotating_vorticity_equation_V2_derivation` | the design's LeviCivita cell (steps 1–13, nothing dropped); curl(5.25) ≡ (5.30) on a divergence-free polynomial field with variable ρ, p; `vorticity_budget_sym` | ✓ (the book's silent drop of u_{j,j}(ω_n + 2Ω_n) confirmed) |
| D16 | ★★ | `test_natural_coordinates_V2_derivation` | (ω·∇)u = ω G e_s; projections; helix T, κ, τ, book e_n = −N | ✓ |
| D17 | ★ | `test_uniform_strain_V1_…` | numeric (5.32) numbers | ✓ |
| D18 | ★★ | `test_burgers_balance_V2_derivation` | ω ∝ L; steady ODE; total derivative; bracket 0; ∫ω dA = Γ; core; curl of u_φ; = `burgers_balance` | ✓ |
| D19 | ★★ | `test_absolute_circulation_V2_derivation` | triple product; DA_vec/Dt = ∮u × dx for a generic material loop (steps 6–7); ∮(Ω×x)·dx = 2Ω·A_vec on a non-planar loop | ✓ |
| D20 | ★★ | `test_fluid_column_V2_derivation` | (ζ + f)A, Ah ⇒ (ζ + f)/h; Dζ/Dt = (ζ + f)/h Dh/Dt | ✓ |
| D21–D23 | ★ | pair, wall, sheet tests | numeric | ✓ (Fig. 5.11 caption slip: G is not a stagnation point unless Γ₁ = Γ₂ — asserted) |

## Functions used by the notebook and explainers (design Part C) — test name each
`test_contract_V1_every_part_c_name_exists_and_is_reexported` checks **all 119 Part C names** (C.1–C.5) plus 30 public
additions the scripts and notebook use exist on `ch05`, are the core objects (re-exports, not copies) and keep the
contracted leading parameters; `test_contract_V1_every_part_c_name_is_exercised_in_this_file` asserts that every one is
called as `ch05.<name>` somewhere in the suite with a physical assertion; `test_result_types_V1_named_fields` pins the
NamedTuples; `test_scalar_callable_V1_…` checks that 34 functions the explainers' parity rows call return plain floats;
C.6 drawing helpers: `test_scripts_V1_drawing_helpers_run` (all nine); every script: `test_scripts_V1_every_ch05_script_runs`.
All design contract numbers (Part C, "Contract numbers …") are asserted, e.g. 0.632121, 125 Pa, 12.742 mm,
0.093629/0.106371, −1266.515, −0.129104, −2533.030, 0.25, −0.0318310, −0.002, 0.0315127/0.0318310/3.18310e-4, 31.831 m/s,
−607.93/−1215.85 Pa, 0.00464739/−3.34538e-4, 1.155130, 1.098212, 1.985865, 2.422222, 2.5e-7/−7.7047e-7/−0.0981/1 + 1.25e-9,
0.002 m, 79.5775, 282.095, −0.260250, 7.389056, 2.0, 1.0, (20, 5e-5, 1e-3), κ/τ, 0.308838, 0.112540, 1/2π, 1/4π, 0.900316,
1.0, 0.353553, 1.054786/1.013052/1.000804, 0.328816, 1.45842e-9, 1.0e-5, −4.19261e7/−5.33820e-5, f(30/45/60°),
pair values, 0.25, H = 0, ±2.0, 0.043968/0.0044131/0.00044128, −0.936549/−0.936551.

## Convergence studies (scheme | steps | observed order (pairwise) | design order)
| scheme | steps | observed | design |
|---|---|---|---|
| curl stencil (`vorticity_field`, ABC) | h = 0.1 … 0.0125 | 2.000 (1.999, 2.000, 2.000) | 2 |
| (5.13) residual, nested stencils (Burgers) | h = 4e-5 … 5e-6 m | 1.999 (1.998, 2.000, 2.001) | 2 |
| (5.13) residual, fully nested (Lamb–Oseen) | h = 8e-4 … 1e-4 m | 1.99 (1.96, 1.99, 2.00) | 2 |
| polygon ring at the centre (5.17) | M = 16 … 128 | 2.007 (2.017, 2.004, 2.001) | 2 |
| elliptic ring vs polygon | M = 32 … 256 | 1.997 (1.993, 1.998, 2.000) | 2 |
| pressure torque → (5.28) | R = 0.4 … 0.05 m | 2.000 (2.001, 2.000, 2.000) | 2 |
| (5.9) vs central dΓ/dt of a material loop | dt = 0.2 … 0.025 s | 1.996 (1.992, 1.998, 1.999) | 2 |
| tube-strength flux route (ch02 midpoint disc) | nr = 32 … 256 | 2.000 | 2 |
| traced-tube end flux vs polygon lines | n_lines = 8, 16, 32 | ≈ 2 (±0.25) | 2 |
| `baroclinic_term` vs sympy | h = 0.1 … 0.0125 | 2 (±0.15) | 2 |
| FTCS vs diffusing sheet | dy (n = 101 … 401) | 2 (±0.15) | 2 |
| (5.15) midpoint rule | n = 8 … 32 | 2 (±0.15) | 2 |
| Biot–Savart solenoidality (stencil) | h = 4e-3 … 1e-3 | 2 (±0.15) | 2 |
| `biot_savart_2d` interior on the dual grid | n = 100 … 400 | > 3.5 | ≥ 2 |
| discrete sheet L1 | N = 10, 100, 1000 | 0.999 (0.998, 1.000) | 1 |
| Biot–Savart tube quadrature | 16³ … 40³ nodes | error 4.1e-3, 2.1e-4, 1.4e-5, 1.2e-6 (spectral) | spectral |

## Conservation / invariant residuals
Kelvin: cellular loop Γ constant to 8.2e-13 (relative 7e-13) while stretched ×6.68; Gaussian straddle 1e-9; rotating
frame Γ_a = π to 1e-9; Helmholtz patches: flux < 1e-10, ABC Γ = 8.016e-10 constant to 1e-6 relative; tube Gauss budget
total 2.2e-12 (axisymmetric), ≤ 5e-9 (traced polyhedron); point vortices: P, I, H over 20 s ≤ 1e-10; equal pair back
to 1e-8 after 10 periods; leap-frogging rings ΣΓπR² 1.7e-9 relative over 30 s, a²R 1e-12; sheet roll-up centroid 1e-17;
cylinder energy budget dissipation = power in − out (residual 7e-18 W/m); tank volume residual < 1e-12.

## Benchmarks used (value | our value | source + URL | date verified)
| value | ours | source | verified |
|---|---|---|---|
| Kelvin thin-ring speed U = Γ/4πR[ln(8R/a) − ¼] (uniform core) | identical at three (R, a, Γ) (1e-14); 0.328816 m/s | Wikipedia "Vortex ring", https://en.wikipedia.org/wiki/Vortex_ring | 2026-09-23 |
| corotating pair Ω₀ = 1/(4πl²); counter-rotating V₀ = −i/(4πl) (vortices at ±l) | `vortex_pair` identical (1e-14); integrated angle and drift 1e-9 | García & Haziot, Commun. Math. Phys. (2023), doi:10.1007/s00220-023-04741-6, arXiv:2204.11327 §2.1, https://arxiv.org/html/2204.11327 | 2026-09-23 |
| Burgers vortex (form; α_w = α/2) | identical (1e-12) | Wikipedia "Burgers vortex" | 2026-09-23 (V1 form) |
| Hill's spherical vortex ψ_in, ψ_out, ω_φ (form; U_w = −U, A = 15U/2a²) | ψ identical inside and outside (1e-12), ω_φ (1e-12) | Wikipedia "Hill's spherical vortex" | 2026-09-23 (V1 form) |
| Biot–Savart segment Γ/(4πr)(cos A − cos B), line Γ/2πr | identical | Wikipedia "Biot–Savart law" | 2026-09-23 (V1 form) |
| Kelvin / Poincaré–Bjerknes statement ∫(ω + 2Ω)·n dS | (5.33) as coded | Wikipedia "Kelvin's circulation theorem" | 2026-09-23 (V1 form) |
Reused: WGS-84 Ω (reference/ch04). Wikipedia forms are V1 "form cross-checks (not V5)", as in ch04.

## Numbers from the text (book vs ours)  [private values redacted to relative errors]
| book item | our route | relative error | comment |
|---|---|---|---|
| §5.1 printed pressure forms (tank, line vortex), σ_rθ, rotating-cylinder u_θ and Γ | `solid_body_pressure`, `line_vortex_pressure`, `line_vortex_viscous_stress`, `rotating_cylinder_flow` | ≤ 1e-14 | forms identical |
| §5.7 V₁, V₂, opposite-pair and wall speeds | `vortex_pair`, `vortex_near_wall_speed` | ≤ 1e-15 | |
| §5.8 text and caption strengths | `vortex_sheet_strength` ccw / caption | 0 | both conventions pinned |
| Exercise 5.1 uncovered area | `rotating_tank_free_surface(closed=True)` | −0.23 % | book rounds ("nearly") |
| Exercise 5.2 tornado circulation | `tornado_circulation` with ρ(25 °C) = 1.184 kg/m³ | −0.13 % | property value |
| Exercises 5.3, 5.5, 5.6, 5.11 (inside), 5.12, 5.13, 5.19, 5.20 closed forms | helical swirl, lock rate, sheet, Hill, Burgers, channel series, cylinder sheet, wall image | ≤ 1e-7 | identical |

## Figures reproduced with our code (outputs/ch05/verify/, local)
| figure | file | one-sentence visual verdict |
|---|---|---|
| Figs. 5.2, 5.3 | `fig5_2_5_3_isobars_verify.png` | Three tank isobars are parallel paraboloids 20 mm apart rising 0.115 m over 0.3 m, with the volume-kept free surface above them; the line-vortex funnels plunge to −∞ at the axis while the Rankine isobars flatten to finite bottoms (−0.05, −0.20 m) inside a = 0.1 m and coincide with the funnel outside. |
| Fig. 5.1 / (5.4) | `fig5_1_tube_budget_verify.png` | The flux line is flat at 0.632 while the area falls and the mean vorticity rises exponentially along z; the Gauss bars are −0.632/0/+0.632/0 for the vorticity field and −0.632/0/+0.233/−0.400 for the broken field. |
| Fig. 5.4 | `fig5_4_kelvin_verify.png` | The circle is wound into a long spiral (×6.7 length) while its relative change of Γ stays within ±7e-13, and the viscous Lamb–Oseen circle's Γ falls from 4.65e-3 to 0.85e-3 m²/s exactly on the closed-form curve. |
| Figs. 5.5, 5.6 | `fig5_5_5_6_baroclinic_verify.png` | The torque-route error is a straight log–log line of slope 2.00; the spin-up is a sine of the isopycnal tilt (±0.049 s⁻² at ∓90°, 0 at 0°); the lock-exchange source is a vertical band peaked at 2.42 s⁻² on x = 0. |
| C07 (Fig. 5.8 context) | `biot_savart_tube_verify.png` | Biot–Savart, (5.14) with +1/(4π) and the 1-D reference coincide (peak 0.86 m/s near r = 0.13 m, below the infinite-tube 1.0 because the tube is 4 m long), while the printed −1/(4π) gives the exact mirror image below zero. |
| Figs. 5.11–5.14 | `fig5_11_5_14_point_vortices_verify.png` | The equal pair circles its midpoint on r = 0.5, the unequal pair circles G at x = 0.25 (radii 0.75 and 0.25); the opposite pair climbs two parallel straight lines; the knife pair inside the bucket runs to the wall and splits along it; wall arrows lie flat (max ∣v∣ = 0). |
| Fig. 5.15 | `fig5_15_rings_verify.png` | The ring rises to the wall at z = 3, turns and spreads along it with its image mirrored above; two coaxial rings alternate radii 0.77–1.19 m ten times while advancing 16 m. |
| Fig. 5.16 | `fig5_16_sheet_verify.png` | Rows of 10, 100, 1000 filaments approach the continuous jump from −1 to +1 m/s at y = 0 (N = 1000 indistinguishable), L1 order 0.999. |
| Fig. 5.10 + Burgers | `fig5_10_column_burgers_verify.png` | Ridge ζ dips to −2.1e-5 s⁻¹ and trough ζ peaks at +2.1e-5 s⁻¹ at x = 0, the slope goes monotonically from −2 to +2e-5; Burgers' stretching (79.6 at R = 0) is balanced by diffusion inside 2 mm and by advection outside, residual flat at 0. |
| F2 evidence | `hill_outside_F2_verify.png` | Loop 1: the exterior arrows cut across the ψ contours and u·e_r on r = a⁺ swung ±0.154 m/s; loop 2: the arrows follow the ψ contours and u·e_r is ≤ 4e-7 on both sides (the title still says "F2: fluid crosses" — it is the evidence plot, now flat). |
| F1 evidence | `core_edge_F1_verify.png` | Loop 1: a single spike to half the torque at r = a with a −100 N/m³ force; loop 2: the torque is −1.2566e-4 N m/m for every r ≥ a with no spike, and the force is 0 outside with a gap (NaN) at r = a. |
| convergence | `convergence_verify.png` | Eight straight log–log lines, seven of slope 2.00–2.01 and the sheet's 1.00, none flattening at round-off. |

## Deviations & justifications
Two `# DEVIATION` markers in the code, both tested: `velocity_from_curl_omega` uses +1/(4π) (the printed (5.14) sign is
kept as `sign=−1` only; `test_velocity_from_curl_omega_V1_corrected_sign_and_printed_sign`, `test_poisson_green_V2_derivation`)
and `vortex_sheet_strength` follows the text's u₂ − u₁ with the caption's sign as `convention="caption"`
(`test_vortex_sheet_strength_V1_…`, V6). Book slips handled and pinned by tests: (5.14) sign and its compensating
step-5 slip (D10/D11 tests); the dropped u_{j,j}(ω_n + 2Ω_n) in (5.27) (D15 test); Fig. 5.2's "2ω" (ω = vorticity, WV);
"single valued" is not why ∮dp/ρ = 0 (baroclinic counterexample); Fig. 5.11's caption (G is not a stagnation point);
Fig. 5.16's caption sign. Test-side choices recorded honestly: stencil-residual tolerances are relative to the size of
the terms (for Lamb–Oseen to ∣u∣∣ω∣/σ, the scale of the cancelling advective product; the Lamb–Oseen case supplies the
exact ω); the smoothed square loop's *length* converges algebraically (its parametrisation is C² at the corners),
asserted as order > 2.5 and 1e-6 at 512 points.

## Open items
- **F1, F2** — resolved in loop 2 (above).
- **O1 (docstring) — resolved in loop 2:** `biot_savart_2d` said u_θ is recovered "outside and inside the core". Outside it is spectral
  (< 1e-8); inside, only at points on the dual (cell-corner) grid (order > 3.5); at arbitrary interior points with
  eps = 0 the singular kernel leaves 1–7 % errors that do not converge with n. The docstring now says so and advises
  eps ≈ the spacing.
- **O2 (notebook/explainer note):** the flux route of `vortex_tube_strength` inherits ch02's midpoint disc (order 2): with
  the default `planar_disc` (nr = 32) the Lamb–Oseen flux is 0.632232 against 0.632121 (1.8e-4); the contract's
  six-digit "both routes" needs nr ≥ 256 (2e-6). Quote four digits or pass a finer disc.
- **O3 (qualitative, as labelled in `ring_dynamics`):** in the thin-core model a ring near a wall keeps widening
  (R = 1 → 75.7 m in 20 s, 2433 m in 60 s at Γ = 1 m²/s, a₀ = 0.1 m, wall 3 m away) while hugging the wall with a
  shrinking core (a²R conserved); the asserted properties (R monotone up, approach speed monotone down, z < wall) hold,
  but explainers should stop the run while the gap to the wall exceeds a few core radii.
- **O4 (qualitative, as labelled in `sheet_rollup`):** the roll-up shape is qualitative; its kernel (cot form), the
  centroid invariant and the Kelvin–Helmholtz growth (> ×10 by t = 4L/γ) are asserted.
- **O5 (designer, wording) — resolved in loop 2 (design Part F D07 step 5 edited):** D07 step 5 says ρ² ≈ ρ̄² "because the contrast is small"; at the interface centre of the
  tanh profile ρ = ρ̄ exactly, so the result 2(ρ₂ − ρ₁)g/((ρ₂ + ρ₁)δ) is exact there (the approximation only matters off
  the centre).
- **O6 (designer) — resolved in loop 2 (D12 step 3 now reads "30 % at 10 core radii, 3 % at 100"):** D12 step 3 says the frozen-kernel change is "about
  3a/|x − x′| — 3 % at 10 core radii". The formula is right (sympy series in `test_filament_V2_derivation`: 3a/D), but at
  D = 10a it is **30 %**; 3 % needs D = 100a. Every other Part F "live"/check number was recomputed and holds (D02
  2500 Pa/m and 12.7 mm, D06 −7.70e-7, D07 2.42, D13 0.1125, D15 1.458e-9, D19 −4.19e7, D20 1.0e-5, D21 −1.70 m/s).
- Nothing is labelled `unverified`.

## Verdict: PASS
Loop 2 of 3. All 9 scripts ran; `tests/test_ch05.py` 133/133 pass (129 + 4 skipped without the private file); the full
suite **611 passed in 329 s** (478 earlier + 133 ch05); every computable CORE row has ≥ 2 independent levels (one of V1/V2/V3/V5), every coded NOTE row ≥ 1,
all 119 Part C names (+ 30 extras) are exercised, all 18 ★★/★★★ derivations are re-derived symbolically, all 40 planted
wrong variants and the 4 re-planted F1/F2 variants are caught. F1 and F2 are fixed and re-checked independently
(the edge stress jump −μΓ/πa² from both sides' profiles); O1, O5, O6 resolved. Remaining open items O2–O4 are notes
(tolerance of a default disc, two models labelled qualitative), not defects.
