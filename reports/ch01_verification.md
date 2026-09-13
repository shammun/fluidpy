# Chapter 1 verification — Introduction                     2026-09-13 (update after review), commit a50bef3 + uncommitted review fixes

Verifier: math-verifier. Suite: `tests/test_ch01.py` (107 tests), figures/metrics: `tests/ch01_verify_figures.py`,
cited data: `reference/ch01/` (`make_refs.py`, `SOURCES.md`). Nothing in `fluidpy/` or `scripts/` was edited by the
verifier. This update answers `reports/ch01_review.md` (Must fix 1 test side, Should fix test items, labels).

## Environment
Python 3.11.5 · numpy 2.4.6 · scipy 1.17.1 · sympy 1.14.0 · pint 0.25.3 · matplotlib 3.11.2 · plotly 7.0.0 (Windows, `.venv`).

Runs: `pytest tests/test_ch01.py` → **107 passed, 0 failed, 0 skipped in 11.5 s** (with the private book file;
without it 102 passed, 5 skipped — the skip guard works). Full suite `pytest tests` → **115 passed in 15.6 s**.
All 14 `scripts/ch01_*.py --no-show` exit 0 (2–4 s each). `tools/check_public.py` → OK (127 files).

Tests per evidence tag (a test carrying two tags counts for both): V1 58 · V2 19 · V3 8 · V4 7 · V5 13 · V6 5 · V7 19 ·
smoke 2.

Changes since the first report (review response):
- `test_synthetic_column_V1_layers_consistent`: default slope is `ch01.MIXED_LAYER_DT_DZ` = `adiabatic_lapse_rate()`;
  every level below 790 m is `"neutral"` with N² = 0 (max |N²| < 1e-15 s⁻², θ ptp < 1e-9 K); inversion and upper
  layers `"stable"`; the old −9.8e-3 slope is shown to be `"unstable"`. Finding F1 closed.
- New `test_parcel_equation_V4_energy_first_integral_divides_by_parcel_density` pins the division by the parcel density
  (proof that it fails on an environment-density version below).
- New `test_seawater_eos_V5_unesco_eos80_check_values` from the primary UNESCO Tech. Pap. 44 table. O2/F4 closed.
- The book's worked scale-height example is no longer asserted in the public test; public test uses 288.15 K and
  H = RT/g analytically; the book value stays in the private V6 test.
- New evidence for the new keyword options: `standard_atmosphere(geometric=True)` (V5), `blast_*(geometry="hemisphere")`
  (V1), `ftcs_diffusion_1d` unknown `bc` (V1), `couette_startup_profile` at t = 0 for h = 1e-12…1e3 (V1),
  `adiabatic_lapse_rate(T=…)` without `alpha` warns (V1).
- Re-labels: wedge (C18) V3 → V1; nonlinear parcel ζ₀ → 0 (C50) V3 → V7; N² three routes (C51) V4 → V1 consistency;
  lapse-rate sweep (C54) V4 → V7 (invariance under the sign convention); Buckingham property tests (C69) V4 → V1;
  unit-change invariance (C64) V4 → V7; the two statistical C06 checks V4 → V1 (statistical); θρ_θ = p_ref/R V4 → V1
  (an identity). V4 now counts only genuine conservation laws / state-function invariants.

## Validation table — A items (CORE, ≥ 2 independent levels, one of V1/V2/V3/V5)
| Concept / Eq. | Tier | fluidpy target | Evidence (V-levels) | Numbers | Label | Notes |
|---|---|---|---|---|---|---|
| C06 continuum hypothesis (§1.4), D34, D35 | CORE | `sample_density`, `density_noise_expected`, `box_average_density`, `molecular_pressure`, `wall_impact_pressure` | V1 Poisson law on the sampler; V1 from-scratch uniform positions (binomial); V1 box average vs quadrature; V2 sympy D35 and D34; V1 (statistical) sampled molecules → n k_B T; V1 (statistical) sample mean → box average | noise slope −1.503 (design −1.5); ratio to (nL³)^−1/2 within 5 SE; kinetic p within 1 % (N = 10⁶, SE 0.08 %); box average vs quadrature < 1e-8 | analytic, symbolic | statistical tolerances are 5 standard errors |
| C12 Newton viscosity (1.3), (1.4) | CORE | `shear_stress_profile`, `newton_shear_stress`, `ftcs_diffusion_1d`, `couette_startup_profile`, `wall_shear_history` | V1 exact stress of a parabola; V3 stencil order; V3 FTCS → Couette series; V2 pint dims of (1.1)–(1.4); V4 zero-flux integral; V7 stability limit; V1 t = 0 initial state for any gap; V1 bc names + periodic discrete decay; V5 USSA μ, ν | stencils 1.996; FTCS vs series 2.002 (pairwise 2.013, 1.998, 1.999, 1.999); steady wall stress μU/h to 1e-6; integral drift 2e-16; periodic sine vs discrete amplification < 1e-13 | analytic, converged, conserved, benchmark | |
| C20 hydrostatic law (1.8), D05, D37 | CORE | `integrate_hydrostatic`, `hydrostatic_pressure_uniform`, `atmosphere_from_temperature`, `standard_atmosphere`, `net_pressure_force_on_box` | V2 sympy residual + pint; V1 integrator vs uniform/isothermal/polytropic closed forms; V3 from-scratch Euler order 1; V5 USSA-1976 Table 1 T, p, ρ (geopotential input and the new `geometric=True`); V1+V4 box face integral = ρgV, horizontal = 0 | closed forms < 1e-8; Euler 1.002; USSA max \|ΔT\| 4.6e-4 K, p 1.2e-4, ρ 8.5e-5 rel (identical with `geometric=True`; 10 km geometric → 223.2521 K vs 223.252); box Fz/ρgV − 1 = 2e-16 | analytic, symbolic, converged, benchmark | |
| C25 first law (1.10), (1.11) | CORE | `process_heat_work`, `process_path`, `path_heat_work_totals` | V1 isothermal/isochoric closed forms; V3 trapezoid order; V4 Δe, Δs path-independent while q, w differ | w = −RT ln 2 to 1e-6; order 1.999; ptp(Δe) < 1e-9, ptp(Δs) < 1e-12, ptp(q) > 1e4 J/kg | analytic, converged, conserved | |
| C35 Gibbs relations (1.18), D10 | CORE | `perfect_gas_entropy_change(_p)`, `entropy_change_reversible` | V2 sympy D10; V1 both forms agree; V4 closed reversible cycle Δs → 0 | cycle Δs = −8.4e-6 J/(kg K) with net w = 4896 J/kg | symbolic, analytic, conserved | |
| C36 speed of sound (1.19) | CORE | `sound_speed_from_eos`, `perfect_gas_sound_speed`, `tait_pressure` | V1 p = Kρ^γ exact; V5 USSA c column; V7 Tait c ∝ √K0; V3 central-difference order | c(1.225) = 340.29 m/s; USSA c rel 2.1e-5; order 2.000 | analytic, benchmark | |
| C40 p = ρRT (1.22), D11 | CORE | `perfect_gas_density/pressure/state`, `gas_constant`, constants | V5 CODATA k_B, N_A, R; V5 USSA sea level, p/(ρT), Table 2 n, speed, H_p; V5 Table 3 M; V1 D11 chain; V2 pint + kmol trap | ρ(101325, 288.15) rel −1.8e-5; M0 3.7e-6; n 3.8e-5, speed 1.9e-5, H_p 8.2e-5 | benchmark, analytic, symbolic | R_AIR uses CODATA R_u |
| C45 isentropic law (1.25), D14 | CORE | `isentropic_pressure`, `isentropic_ratios` | V2 sympy D14; V3 from-scratch Euler; V1 ratio consistency; V7 γ → 1 | Euler 0.998; pump 351.3 K, ρ ratio 1.641 | symbolic, converged, analytic | |
| C50 parcel equation, D18 | CORE | `parcel_displacement`, `parcel_ode`, `parcel_ode_from_gradients`, `parcel_ode_atmosphere` | V2 sympy D18 (steps 3–12); V1 code solution satisfies ζ″ + N²ζ = 0; **V4 exact energy first integral of the nonlinear parcel (÷ρ_p)**; **V1 turning point and period vs quadrature**; V7 nonlinear → linear as ζ₀ → 0 (ocean, atmosphere); V7 runaway cap → NaN | ρ0 = 1, dρ/dz = −0.01, dρ_a/dz = +0.005, ζ0 = 40 m: ptp(E)/\|V(ζ0)\| = 3.6e-8 (tol 1e-6); turning point −35.2876 m rel err 3.4e-11 (÷ρ_e would give −54.74 m); period 16.40641 s rel err 1.6e-12 (linear 16.382 s); limit order 2.000 (ocean), 2.002 (atmosphere) | symbolic, analytic, conserved | review test gap closed |
| C51 N² (1.29), D36 | CORE | `brunt_vaisala_sq(_from_lapse/_theta)`, `stability_timescale`, `classify_stability` | V1 thermocline number, period vs zero crossing; V2 sympy isothermal N² = g²/(C_pT); V1 consistency of three independent routes (density, lapse, θ) | route differences 1.9e-7, 1.8e-7 of max\|N²\|; period 642 s | analytic, symbolic | |
| C54 adiabatic lapse rate (1.30), D19, conventions | CORE | `adiabatic_lapse_rate`, `lapse_rate_convention`, `lapse_rate_stability`, `lapse_rate`, `parcel_temperature` | V2 sympy D19 incl. every intermediate line; V1 −g/C_p, general form, from-scratch lifted-parcel slope, T-without-α warning; V5 AMS 9.8 °C/km; V6 book value; V1 exact texts; V7 invariance of the verdict under the convention (sweep) | Γa = −9.7607 K/km (AMS rel −0.4 %); release slope rel 2e-11; 2504 sweep points identical codes | symbolic, analytic, benchmark | F3 still open (Nice) |
| C55 potential temperature (1.31), (1.32), D20, D36 | CORE | `potential_temperature`, `potential_temperature_gradient`, `potential_density`, `temperature_from_potential` | V1 θ = T at p_ref, 250 K/500 hPa → 304.8 K; V2 sympy D20, (1.32), D36; V4 θ constant along an independently integrated dry adiabat; V1 identity θρ_θ = p_ref/R; V1 (1.32), (1.34) vs FD | θ ptp/θ = 3.8e-12; θρ_θ residual 4.4e-16; (1.32) vs FD < 1e-6 | analytic, symbolic, conserved | |
| C64 dimensional homogeneity | CORE | `units.dimensional_check`, `check_dimensions`, `dimension_vector`, `rescale_units`, `group_value` | V2 pint checks of chapter laws + a wrong law fails; V1 dimension vectors vs pint dicts; V7 group values invariant under SI → cgs → imperial | Π₁ = 10 in all systems to 1e-12 | symbolic, analytic | |
| C67 dimensional matrix (1.39) | CORE | `dimensional_matrix`, `minor_determinant`, `rank_by_minors` | V1 exact (1.39), all 35 minors by cofactor expansion; V2 rank vs numpy and sympy on 200 random matrices; V6 book matrices | 200/200 ranks agree | analytic, symbolic, book-value | |
| C69 Buckingham Π (1.37), D28 | CORE | `pi_groups`, `solve_exponents`, `groups_independent`, `group_*`, presets | V2 sympy D28 nullspace, rank–nullity, span, rescaling; V1 exponent solve vs `np.linalg.solve`; V1 properties on 7 presets × 2 repeating sets; V1 Poiseuille collapse | collapse 1e-12; rescaled repeating variables = 1 to 1e-12 | symbolic, analytic | |

Every computable A item has ≥ 2 independent levels with at least one of V1/V2/V3/V5.

## Validation table — coded B/C items (NOTE, ≥ 1 level)
| Concept | fluidpy target | Evidence | Numbers | Label |
|---|---|---|---|---|
| C02/N01 fluid vs solid | `shear_deformation_history` | V1 + V7 | exact closed forms; μ → ∞ no flow | analytic |
| C03 normal/shear stress | `traction_components` | V1 | < 1e-12 | analytic |
| C04 molecular spacing | `number_density`, `mean_molecular_spacing` | V1, V5 USSA n | n rel 3.8e-5 | analytic, benchmark |
| C07/N03 Knudsen, mean free path | `knudsen_number`, `mean_free_path_jennings`, `mean_free_path_air` | V5 Jennings, V2 pint, V7, V6 | 67.18 vs 67.3 nm (−0.18 %) | benchmark |
| C09 mass fraction | `mass_fractions`, `DRY_AIR` | V1, V5 USSA Table 3 | ΣY = 1 | analytic, benchmark |
| C10/C11 Fick, Fourier | `fick_mass_flux`, `fourier_heat_flux` | V1, V2 pint | | analytic |
| C13 μ(T) | `sutherland_viscosity`, `viscosity_power_law`, `water_viscosity` | V5 USSA Table 2; V5 IAPWS R12-08 (one point); V7 | μ rel 3.1e-5; water 298.15 K +0.079 % | benchmark |
| C14 ν, κ, h²/ν | `kinematic_viscosity`, `thermal_diffusivity`, `diffusion_time`, `FLUIDS`, `fluid_properties` | V5 USSA ν; V2; V1 | ν rel 5.3e-5 | benchmark |
| C15/C16/N07 surface tension, Laplace | `surface_tension_water`, `laplace_pressure_jump` | V5 IAPWS Table 1, V1, V2 | max \|Δσ\| 0.0045 mN/m | benchmark, analytic |
| C17/C21 gauge, layered tank | `gauge_pressure`, `absolute_pressure`, `layered_pressure` | V1 | | analytic |
| C18 isotropy (1.6) | `wedge_pressure_difference`, `wedge_face_forces` | V1 closed-form balance ∝ dz (re-labelled from V3: an exact expression, not a discretisation) | fitted power 1.000; residuals < 1e-8 N/m | analytic |
| C22 capillary rise | `capillary_rise(_deg)`, `alpha_from_contact_angle` | V1, V7 | residual < 1e-14 | analytic |
| C23/C24 collision time | `mean_molecular_speed`, `collision_time` | V1, V1 (statistical), V5 USSA speed | speed rel 1.9e-5 | benchmark |
| C28–C31, N13 | `perfect_gas_state`, `specific_volume`, `enthalpy`, …, `specific_heat_cp/cv` | V1, V3 order 2 | | analytic, converged |
| C34 Clausius–Duhem | `irreversible_process`, `free_expansion`, `stirred_isochoric_process` | V1 + V4 (entropy produced > ∫δq/T) | Δs = R ln 2 > 0 | analytic |
| C37/C48 α | `thermal_expansion_coefficient`, `perfect_gas_expansion_coefficient`, `water_density` | V1, V7 | ρ_max at 277.13 K ± 0.05 | analytic |
| C41 van der Waals | `van_der_waals_*` | V1 Gibbs identity numerically | rel < 1e-5 | analytic |
| C42/C43/N16 | `cv_from_cp`, `gamma_from_cp`, `cp_from_gamma` | V1 | | analytic |
| C46/C47 | `isentropic_ratios`, `perfect_gas_sound_speed` | V1, V5 | | analytic, benchmark |
| C49 static medium | `atmosphere_from_temperature`, `linear_lapse_pressure`, `standard_atmosphere` (+ `geometric=True`) | V1, V5 | geometric Table 1 match as in C20 | benchmark |
| C52 classification | `classify_stability` | V7 | | analytic |
| C56–C59 | `potential_temperature_gradient`, `potential_density` | V1, V2, V4 | | symbolic |
| C60/N21/C61 ocean | `seawater_density_linear`, `seawater_linear_coefficients`, `seawater_density_eos80`, `isentropic_density_gradient`, `ocean_potential_density_gradient` | **V5 UNESCO EOS-80 check values (p = 0)**; V1 sign vs −N², c → ∞, perfect-gas identity; linear EOS within 0.5 kg/m³ of EOS-80 | ρ at (S, t68) = (0, 5), (0, 25), (35, 5), (35, 25): Δρ = +7.9e-7, −2.5e-6, −4.7e-6, −1.5e-6 kg/m³ (≤ 5e-6 = half the last printed digit; rel ≤ 4.6e-9); V column rel ≤ 4.6e-9; feeding t68 as ITS-90 misses by ≥ 1.8e-5 | benchmark, analytic |
| C62/C63 isothermal, scale height | `isothermal_pressure/density`, `scale_height` | V2 sympy ODE, V1 H = RT/g at 288.15 K and p(H) = p0/e, V5 USSA H_p, V6 (private book example) | USSA H_p rel 8.2e-5 | symbolic, benchmark |
| N18 synthetic profile | `synthetic_boundary_layer_profile/column`, `MIXED_LAYER_DT_DZ` | V1 layers, neutral mixed layer (N² = 0, θ constant), N² per layer, θρ_θ, pressures vs integration | below 790 m all "neutral", max \|N²\| < 1e-15; p vs integration 1e-7 | analytic |
| C72–C76 examples | `poiseuille_pressure_drop`, `blast_energy/radius` (+ `geometry="hemisphere"`), `TAYLOR_K_GAMMA14`, `rayleigh_scattering_ratio`, `pythagoras_phi`, `wavelength_to_rgb` | V1, V5 Taylor K, V7; V1 hemisphere ≡ free sphere of 2E | hemisphere E = Kρ D⁵/(2t²) and D(E) = D_sphere(2E) to 1e-14; radius ratio 2^(1/5) | analytic, benchmark (`wavelength_to_rgb` qualitative, O1) |
| R03 °C ↔ K | `celsius_to_kelvin`, `kelvin_to_celsius` | V1 | | analytic |

## Functions used by the notebook and explainers (design Part C — 144 callables)
Every callable of Part C is called in at least one test (unchanged set; the new keyword options are covered by the new
tests listed above). Mapping (renamed tests in bold):

| Part C rows | Tests |
|---|---|
| 0.1–0.4 `style`, `anim`, `interact`, `embed.show_viz` | `test_machinery_smoke_style_anim_interact_embed` (smoke) |
| 0.5–0.6 units | `test_dimensional_homogeneity_V2_*`, `test_units_V1_*` |
| 1–12 dimensional | `test_dimensional_matrix_*`, `test_buckingham_*` (**`test_buckingham_V1_groups_count_…`**), `test_dimensional_homogeneity_*` (**`…_V7_values_change_groups_do_not`**), `test_examples_V1_*` |
| 13–35 thermo | `test_perfect_gas_*`, `test_first_law_*`, `test_gibbs_*`, `test_clausius_duhem_*`, `test_specific_heats_*`, `test_sound_speed_*`, `test_thermal_expansion_*`, `test_isentropic_*`, `test_specific_heat_relations_*`, `test_equation_of_state_*` |
| 36–45 statics | `test_hydrostatics_*`, **`test_standard_atmosphere_V5_geometric_height_option`**, `test_buoyancy_*`, `test_isothermal_atmosphere_*` |
| 46–63 stratification | `test_brunt_vaisala_*` (**`…_V1_density_form_equals_lapse_and_theta_forms`**), `test_parcel_*` (**`…_V7_nonlinear_parcel_converges_to_linear`**, **`…_V7_atmospheric_parcel_converges_to_linear`**, **`…_V4_energy_first_integral_divides_by_parcel_density`**), `test_stability_classification_*`, `test_adiabatic_lapse_rate_*`, `test_lapse_rate_*` (**`test_lapse_rate_stability_V7_same_code_…`**), `test_potential_temperature_*`, `test_ocean_criterion_*` |
| 64–67 diffusion | `test_couette_startup_*` (**`…_V1_initial_state_any_gap_width`**), `test_ftcs_*` (**`test_ftcs_V1_boundary_names_validated`**), `test_newton_viscosity_V3_*` |
| 68–88 chapter module | `test_continuum_*` (**`…_V1_sample_mean_tracks_box_average`**), `test_kinetic_pressure_*` (**`…_V1_sampled_molecules_give_nkT`**), `test_knudsen_*`, `test_newton_viscosity_*`, `test_viscosity_temperature_*`, `test_water_viscosity_*`, `test_fluid_properties_*`, `test_fluid_vs_solid_*`, `test_normal_shear_stress_*`, `test_mass_fraction_*`, `test_surface_tension_*`, `test_laplace_jump_*`, `test_capillary_rise_*`, **`test_pressure_isotropy_V1_wedge_closed_form_vanishes_linearly`**, `test_seawater_eos_*` (**`…_V5_unesco_eos80_check_values`**), `test_synthetic_column_*`, `test_examples_V1_*`, `test_blast_*` (**`test_blast_V1_hemisphere_equals_free_sphere_of_twice_the_energy`**) |
| 89–90 drawing helpers | `test_drawing_helpers_smoke` (smoke) |

Functions from Part C without evidence: **none**.

## Convergence studies (discretisations only)
| scheme | grids | observed order | design order |
|---|---|---|---|
| FTCS vs Couette start-up series (r = 0.4, t = 0.05 h²/ν) | N = 11, 21, 41, 81 (+161 in the figure) | 2.002 (pairwise 2.013, 1.998, 1.999, 1.999) | 2 |
| FTCS vs spreading Gaussian, zero-flux ends (r = 0.25) | N = 61…481 | 2.007 | 2 |
| `shear_stress_profile` / `derivative_2nd_order` incl. end stencils | N = 21…161 | 1.996 | 2 |
| `partial_derivative` central difference | h = 0.1…0.0125 | 2.000 | 2 |
| `process_heat_work` cumulative trapezoid (isentrope) | n = 11…81 | 1.999 | 2 |
| `lapse_rate` (np.gradient edge_order 2, non-uniform grid) | N = 21…161 | within 2 ± 0.25 | 2 |
| from-scratch Euler hydrostatic march | n = 50…400 | 1.002 | 1 |
| from-scratch Euler isentrope dp/dρ = γp/ρ | n = 100…800 | 0.998 | 1 |
| FTCS stability limit | r = 0.49 decays; r = 0.51 grows > 10³×; r > 0.5 raises | — | r ≤ ½ |

Limits (V7, not convergence): nonlinear parcel → linear as ζ₀ → 0, error ∝ ζ₀^2.000 (ocean, ζ₀ = 8…1 m) and
ζ₀^2.002 (atmosphere, 200…25 m); wedge p₂ − p₁ ∝ dz^1.000 (closed form, V1).

## Conservation / invariant residuals (V4)
| invariant | residual |
|---|---|
| FTCS trapezoid integral, zero-flux ends, 2000 steps | 2.0e-16 relative |
| energy first integral ½ζ′² + V(ζ) of the nonlinear parcel (÷ρ_p), 1.25 periods | ptp(E)/\|V(ζ₀)\| = 3.6e-8 (O(dt²) of the central-difference velocity) |
| θ along an integrated dry adiabat (0–10 km) | ptp/θ = 3.8e-12 |
| closed reversible 4-leg cycle | Δe = 0, Δs = −8.4e-6 J/(kg K) (trapezoid), q + w = 0 to 1e-9 |
| Δe, Δs over 4 two-leg paths (state functions) | ptp 1e-9, 1e-12 |
| box face integral of (1.9) | F_x = F_y = 0; F_z/ρgV − 1 = 2.2e-16 |

Consistency / invariance checks (not V4): N² by three routes 1.9e-7 of max\|N²\| (V1); θρ_θ = p_ref/R 4.4e-16 (V1
identity); lapse-rate verdict identical in both conventions on 2504 points (V7); Π-group values under SI → cgs →
imperial < 1e-12 (V7).

## Proof that the parcel-density test discriminates (review test gap)
Scratch copy of `fluidpy/`, `tools/`, `tests/`, `reference/`, `scripts/` in the session scratchpad (never in the repo),
with one line of `fluidpy/core/stratification.py::parcel_ode.rhs` changed:
```
-        return [y[1], -g * (rp - rho_env_fn(z0 + y[0])) / rp]  # Newton II with weight and buoyancy
+        return [y[1], -g * (rp - rho_env_fn(z0 + y[0])) / rho_env_fn(z0 + y[0])]  # PATCHED: divides by environment density
```
`pytest -k parcel` on the patched copy: **1 failed, 6 passed** — the six earlier parcel tests (D18 sympy, linear ODE,
both ζ₀ → 0 limits, runaway cap) all still pass, confirming the gap; the new test fails with
`ptp(E)/|V(ζ₀)| = 1.607` (tolerance 1e-6). Full chapter file on the patched copy: that test fails, plus
`test_drawing_helpers_smoke`, which fails only because `book.yaml` was not copied to the scratch tree (unrelated).
On the real code the same test passes with 3.6e-8.

## Benchmarks used
| value | our value | source | date verified |
|---|---|---|---|
| k_B, N_A, R (exact) | identical | NIST CODATA CUU pages | 2026-09-12 |
| USSA-1976 Table 1: T, p, ρ, c at Z = 0…50 km | max \|ΔT\| 4.6e-4 K; p 1.2e-4, ρ 8.5e-5, c 2.1e-5 rel; same with `geometric=True` | PDAS bigtables (from NASA-TM-X-74335) | 2026-09-12 |
| USSA-1976 Table 2: μ, ν, H_p, n, mean speed | 3.1e-5, 5.3e-5, 8.2e-5, 3.8e-5, 1.9e-5 rel | PDAS bigtables | 2026-09-12 |
| USSA-1976 constants, Tables 3–4, r0 = 6356.766 km (errata) | identical; M0 3.7e-6 rel; `USSA_R0` = r0 | NASA-TM-X-74335, NTRS 19770009539 | 2026-09-12 |
| IAPWS R1-76(2014) σ(t) Table 1 | ≤ 0.0045 mN/m | iapws.org Surf-H2O-2014.pdf | 2026-09-12 |
| IAPWS R12-08 μ(298.15 K, 998 kg/m³) | +0.079 % | iapws.org viscosity release, Table 4 | 2026-09-13 |
| Jennings mean free path 67.3 nm | −0.18 % | Tsalikis et al. 2024, doi:10.1080/02786826.2024.2333859 | 2026-09-13 |
| Taylor S(1.4)^−5 = 0.856 | identical | Díaz, arXiv:2009.05674 | 2026-09-13 |
| AMS dry-adiabatic lapse rate ≈ 9.8 °C/km | −0.40 % | AMS Glossary | 2026-09-13 |
| **UNESCO EOS-80 check values ρ(S, t68, 0): 999.96675, 997.04796, 1027.67547, 1023.34306 kg/m³ and V** | rel ≤ 4.6e-9 (≤ 5e-6 kg/m³ absolute) | Fofonoff & Millard, Unesco Tech. Pap. Mar. Sci. 44 (1983), p. 19 (from Unesco Report 38, p. 191); scanned PDF (WHOI/MBL darchive; also JODC), page image read | 2026-09-13 |

The reviewer's cross-check values (1027.67547, 1023.34306, 999.96675) agree with the primary table. The p = 10000 dbar
rows are stored but not used (`seawater_density_eos80` has no pressure argument). Rejected: Engineering ToolBox c_p,
CFD-forum Sutherland constants, Kell (1975) table (primary not fetched), PDAS mean free path (different formula).

## Numbers from the text (book vs ours) — private values redacted to relative differences
| item | relative difference | comment |
|---|---|---|
| k_B, Avogadro per kmol, R_u, M_w of dry air, R of air, γ, C_p | 0.025 %, 0.014 %, 0.006 %, 0.006 %, 0.020 %, 0 %, 0.070 % | book rounds pre-2019 constants |
| atmospheric pressure (kPa, bar) | 0.025 % | rounding |
| adiabatic lapse rate (book's sign) | 2.4 % | book prints a 1-significant-figure value |
| lab-scale temperature change over 1 m | same at 1 s.f. | |
| scale height at the book's temperature | 0.25 % | now asserted only in the private V6 test |
| isothermal band of the lower atmosphere | USSA max deviation 15.26 % vs the book's rounded band | marginal (O3) |
| mean free path of air | ours 30 % larger | order of magnitude in the book |
| pipe matrix (1.39), minors, rank, exponents; Exs. 1.2–1.5 | exact | |

## Figures reproduced with our code
| figure | file (local) | visual verdict |
|---|---|---|
| Fig. 1.9 (a, b) | `outputs/ch01/verify/fig1_9_T_theta_synthetic.png` | the 0–800 m mixed layer is now orange (neutral) and exactly parallel to the grey dry adiabats; the inversion (800–1000 m) and the layer above are teal (stable); θ(z) is a vertical line in the mixed layer, then rises steeply through the inversion and slowly above — matches the book's picture. |
| §1.10 nonlinear parcel (new) | `parcel_energy_first_integral.png` | the dashed `parcel_ode` trajectory lies exactly on the teal closed curve ½ζ′² + V_ρp = const (turning points +40 m and −35.3 m) and well inside the pink ÷ρ_e curve (left turning point ≈ −54.7 m): the code divides by the parcel density. |
| Fig. 1.8 idea (parcel) | `parcel_regimes.png` | stable parcel oscillates at ≈ 642 s, nonlinear on linear; neutral stays; unstable grows like cosh to the cap. |
| Figs. 1.2–1.3 idea | `couette_ftcs_vs_series.png` | FTCS dots on the series curves; error parallel to slope 2. |
| §1.4 continuum | `continuum_noise_and_drift.png` | noise on the (nL³)^−1/2 line; plateau then drift and sinc ringing. |
| lapse-rate criterion | `lapse_rate_sweep.png` | N² linear in dT/dz, zero at Γa = −9.76 K/km in both axes. |
| USSA-1976 | `ussa1976_vs_table.png` | T(H) through every PDAS point, kinks at 11, 20, 32, 47 km. |
| §1.11 pipe collapse | `pi_collapse_pipe.png` | 80 laminar cases on Π₁ = 32 Π₂Π₄. |
| §1.6 σ(T) | `iapws_surface_tension.png` | monotone decrease through all IAPWS points. |

## Derivations re-derived (curation §4b, design Part F)
| D | ★ | test | intermediate lines checked | result |
|---|---|---|---|---|
| D14 isentropic law | ★★ | `test_isentropic_law_V2_derivation_D14` | steps 2, 3, 6, 7, 10 | holds |
| D18 parcel equation | ★★ | `test_parcel_equation_V2_derivation_D18` | steps 3, 6, 7–9, 10, 12; the un-linearised step 3 (÷ρ_p) is now also pinned numerically by the V4 first-integral test | holds |
| D19 adiabatic lapse rate | ★★★ | `test_adiabatic_lapse_rate_V2_derivation_D19` | steps 4, 7, 10, 12, 15, 16 | holds |
| D28 Buckingham Π | ★★★ | `test_buckingham_V2_derivation_D28_nullspace_rank_nullity` | steps 3, 4, 5, 7; 9–11 numerically | holds |
| D34 kinetic pressure | ★★ | `test_kinetic_pressure_V2_derivation_D34` | steps 5, 7, 8 | holds |
| D36 N² from θ | ★★ | `test_potential_temperature_V2_derivation_D36` | steps 2, 3, 5, 6–7, 8 | holds |
| D10, D20, D35, D37 (★) | ★ | `test_gibbs_V2_*`, `test_potential_temperature_V2_derivation_D20_and_eq_1_32`, `test_continuum_V2_*`, `test_buoyancy_V2_*` | as in the tests | holds |

## Deviations & justifications
| `# DEVIATION` in code | assessment |
|---|---|
| `seawater_density_eos80`: salinity in g/kg used as practical salinity (≈ 0.5 %) | acceptable; the UNESCO check values are in practical salinity and are reproduced when S is passed as that number |
| `FLUIDS` solute diffusivity of viscous liquids from Stokes–Einstein scaling | order of magnitude only; positivity and ν = μ/ρ tested |
Also documented: `standard_atmosphere` uses CODATA `R_AIR` instead of USSA R* (1.7e-5); `seawater_density_eos80` takes
ITS-90 kelvin and converts to IPTS-68 (t68 = 1.00024 t90) — the V5 test shows the conversion is required.

## Findings
- F1 (boundary-layer default not neutral) — **closed**: default is Γa, test pins "neutral" below 790 m.
- F2 ("Validation (planned) … Label: pending") — **closed**: no such strings remain in `fluidpy/`.
- F4 (EOS-80 check values unverified) — **closed** by the V5 test.
- F3 · Nice (unchanged) — neutral tolerances of `lapse_rate_stability` (1e-9 K/m) and `classify_stability` (1e-12 s⁻²)
  differ; physically irrelevant, the sweep test excludes the band.
- F5 · Should fix (docstrings only, for the implementer; no test impact) — validation text now stale relative to the
  tests: `stratification.py:287` ("division by the parcel density is not yet pinned by a test" → pinned by
  `test_parcel_equation_V4_energy_first_integral_divides_by_parcel_density`; label can gain "conserved"; ζ₀ → 0 is V7);
  `ch01_introduction.py:1281` ("benchmark pending", O2 → V5 UNESCO Tech. Pap. 44, label benchmark);
  `ch01_introduction.py:1505` and `:1554` ("awaits its test update", "the test still pins the old −9.8e-3 slope");
  `ch01_introduction.py:1686` ("hemisphere factor 1/2 is not yet tested"); "V4 (statistical) … conserved
  (statistical)" on `molecular_pressure`, `wall_impact_pressure`, `box_average_density`, `sample_density`,
  `mean_molecular_speed` (lines 281–283, 321–323, 396, 480–482, 1195–1196) should read "V1 (statistical) … analytic";
  "V4 … conserved" on `dimensional.py:568–570, 613` (unit invariance → V7) and `stratification.py:882–884`
  (θρ_θ = p_ref/R is an identity → V1), per the re-labelling above. (`potential_temperature`'s V4 — θ constant along
  an independently integrated adiabat — stays V4.)

## Open items
- O1 (qualitative, by design): `wavelength_to_rgb` is a visual aid; smoke-tested only (blue-dominant at 450 nm, red at
  650 nm, black outside 380–780 nm).
- O3 (book statement): the book's isothermal-band claim is met by USSA-1976 only at the book's rounding (15.26 % vs a
  2-figure band); the notebook should present it as an order-of-magnitude statement.

(O2 closed 2026-09-13: EOS-80 now has a cited V5 benchmark.)

## Verdict: PASS
