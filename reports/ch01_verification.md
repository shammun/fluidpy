# Chapter 1 verification — Introduction                     2026-09-13, commit a3409d1 (tests uncommitted)

Verifier: math-verifier. Suite: `tests/test_ch01.py` (101 tests), figures/metrics: `tests/ch01_verify_figures.py`,
cited data: `reference/ch01/` (`make_refs.py`, `SOURCES.md`). Nothing in `fluidpy/` or `scripts/` was edited.

## Environment
Python 3.11.5 · numpy 2.4.6 · scipy 1.17.1 · sympy 1.14.0 · pint 0.25.3 · matplotlib 3.11.2 · plotly 7.0.0 (Windows, `.venv`).

Runs: `pytest tests/test_ch01.py` → **101 passed, 0 failed, 0 skipped in 9.7 s** (with the private book file;
without it 96 passed, 5 skipped — the skip guard works). Full suite `pytest tests` → **109 passed in 13.3 s**.
All 14 `scripts/ch01_*.py --no-show` exit 0 (2–5 s each). `tools/check_public.py` → OK.

Tests per evidence tag (a test carrying two tags counts for both): V1 50 · V2 19 · V3 11 · V4 13 · V5 11 · V6 5 · V7 15 ·
smoke 2.

## Validation table — A items (CORE, ≥ 2 independent levels, one of V1/V2/V3/V5)
| Concept / Eq. | Tier | fluidpy target | Evidence (V-levels) | Numbers | Label | Notes |
|---|---|---|---|---|---|---|
| C06 continuum hypothesis (§1.4), D34, D35 | CORE | `sample_density`, `density_noise_expected`, `box_average_density`, `molecular_pressure`, `wall_impact_pressure` | V1 Poisson law on the sampler; V1 from-scratch uniform positions (independent binomial mechanism); V1 box average vs quadrature; V2 sympy D35 and D34; V4 sampled molecules → n k_B T; V4 sample mean → box average | noise slope −1.503 (design −1.5); ratio to (nL³)^−1/2 within 5 SE; kinetic p within 1 % (N = 10⁶); box average vs quadrature < 1e-8 | analytic, symbolic, conserved | statistical tolerances are 5 standard errors of the estimator |
| C12 Newton viscosity (1.3), (1.4) | CORE | `shear_stress_profile`, `newton_shear_stress`, `ftcs_diffusion_1d`, `couette_startup_profile`, `wall_shear_history` | V1 exact stress of a parabola; V3 stencil order; V3 FTCS → Couette series; V2 pint dims of (1.1)–(1.4); V4 zero-flux integral; V7 stability limit; V5 USSA μ, ν | stencils order 1.996; FTCS vs series order 2.002 (pairwise 2.013, 1.998, 1.999, 1.999); steady wall stress = μU/h to 1e-6; integral drift 2e-16 | analytic, converged, benchmark | |
| C20 hydrostatic law (1.8), D05, D37 | CORE | `integrate_hydrostatic`, `hydrostatic_pressure_uniform`, `atmosphere_from_temperature`, `standard_atmosphere`, `net_pressure_force_on_box` | V2 sympy residual + pint; V1 integrator vs uniform/isothermal/polytropic closed forms; V3 from-scratch Euler order 1; V5 USSA-1976 Table 1 T, p, ρ at 0–50 km; V1+V4 box face integral = ρgV, horizontal = 0 | integrator vs closed forms < 1e-8; Euler order 1.002; USSA max \|ΔT\| 4.6e-4 K, p 1.2e-4, ρ 8.5e-5 rel; box Fz/ρgV − 1 = 2e-16 | analytic, symbolic, converged, benchmark | geometric Z converted to geopotential with USSA r0 = 6356.766 km |
| C25 first law (1.10), (1.11) | CORE | `process_heat_work`, `process_path`, `path_heat_work_totals` | V1 isothermal/isochoric closed forms; V3 trapezoid order; V4 Δe, Δs path-independent while q, w differ (4 two-leg paths) | w = −RT ln 2 to 1e-6 (n = 4001); order 1.999; ptp(Δe) < 1e-9, ptp(Δs) < 1e-12, ptp(q) > 1e4 J/kg | analytic, converged, conserved | |
| C35 Gibbs relations (1.18), D10 | CORE | `perfect_gas_entropy_change(_p)`, `entropy_change_reversible` | V2 sympy D10 steps + residuals of both forms; V1 both forms agree (Δs = C_p ln 2 = 696.4); V4 closed reversible cycle Δs → 0 | cycle Δs = −8.4e-6 J/(kg K) (trapezoid, n = 2001 per leg) with net w = 4896 J/kg | symbolic, analytic, conserved | |
| C36 speed of sound (1.19) | CORE | `sound_speed_from_eos`, `perfect_gas_sound_speed`, `tait_pressure` | V1 p = Kρ^γ exact; V5 USSA c column (0–50 km); V7 Tait c ∝ √K0 → ∞; V3 central-difference order | c(1.225) = 340.29 m/s; USSA c max rel 2.1e-5; order 2.000 | analytic, benchmark | |
| C40 p = ρRT (1.22), D11 | CORE | `perfect_gas_density/pressure/state`, `gas_constant`, constants | V5 CODATA exact k_B, N_A, R; V5 USSA sea-level ρ, p/(ρT) on every row, Table 2 n, mean speed, H_p; V5 Table 3 mean M; V1 D11 constants chain; V2 pint + kmol trap | ρ(101325, 288.15) rel −1.8e-5; Table 3 M0 = 28.96451 vs 28.9644 (3.7e-6); n 3.8e-5, speed 1.9e-5, H_p 8.2e-5 | benchmark, analytic, symbolic | R_AIR uses CODATA R_u; USSA R* differs by 1.7e-5 |
| C45 isentropic law (1.25), D14 | CORE | `isentropic_pressure`, `isentropic_ratios` | V2 sympy D14 step by step; V3 from-scratch Euler of dp/dρ = γp/ρ; V1 ratio consistency; V7 γ → 1 | Euler order 0.998; pump T = 351.3 K, ρ ratio 1.641 | symbolic, converged, analytic | |
| C50 parcel equation, D18 | CORE | `parcel_displacement`, `parcel_ode_from_gradients`, `parcel_ode_atmosphere`, `parcel_ode` | V2 sympy D18 (steps 3–12); V1 code solution satisfies ζ″ + N²ζ = 0; V3 nonlinear → linear as ζ0 → 0 (ocean and atmosphere); V7 runaway cap → NaN | order 2.000 (ocean), 2.002 (atmosphere); neutral atmosphere drift < 1e-8 m | symbolic, analytic, converged | |
| C51 N² (1.29), D36 | CORE | `brunt_vaisala_sq(_from_lapse/_theta)`, `stability_timescale`, `classify_stability` | V1 thermocline number, period vs zero crossing; V2 sympy isothermal N² = g²/(C_pT); V4 three independent routes (density form with −ρg/c², lapse form, θ form) on a non-trivial T(z) | route differences 1.9e-7, 1.8e-7 of max\|N²\|; period 642 s | analytic, symbolic, conserved | |
| C54 adiabatic lapse rate (1.30), D19, conventions | CORE | `adiabatic_lapse_rate`, `lapse_rate_convention`, `lapse_rate_stability`, `lapse_rate`, `parcel_temperature` | V2 sympy D19 for an arbitrary G(T, p) incl. every intermediate line (steps 4, 7, 10, 12, 15, 16); V1 −g/C_p, general form, from-scratch lifted-parcel slope; V5 AMS 9.8 °C/km; V6 book value; V1 exact texts + auto-decimals; V4 sweep −15…+10 K/km both conventions vs sign N² | Γa = −9.7607 K/km (meteorology +9.7607, AMS rel −0.4 %); release slope rel 2e-11; elsewhere slope = −(g/C_p)T_p/T_e to 1e-6; 2504 sweep points identical codes | symbolic, analytic, benchmark, conserved | see Findings F1, F3 |
| C55 potential temperature (1.31), (1.32), D20, D36 | CORE | `potential_temperature`, `potential_temperature_gradient`, `potential_density`, `temperature_from_potential` | V1 θ = T at p_ref, 250 K/500 hPa → 304.8 K; V2 sympy D20, (1.32) residual, D36 all steps; V4 θ constant along independently integrated dry adiabat; V4 θρθ = p_ref/R; V1 (1.32) and (1.34) vs finite differences | θ ptp/θ = 3.8e-12; θρθ residual 4.4e-16; (1.32) vs FD < 1e-6 | analytic, symbolic, conserved | |
| C64 dimensional homogeneity | CORE | `units.dimensional_check`, `check_dimensions`, `dimension_vector`, `rescale_units`, `group_value` | V2 pint checks of (1.5), (1.9), (1.22), (1.29), (1.30), Ex. 1.1, H, Ex. 1.4, c + wrong law fails; V1 dimension vectors vs pint dimensionality dicts; V4 group values invariant (cgs, imperial; rescale vs pint conversion) | Π₁ = 10 in SI, cgs, imperial to 1e-12 | symbolic, analytic, conserved | |
| C67 dimensional matrix (1.39) | CORE | `dimensional_matrix`, `minor_determinant`, `rank_by_minors` | V1 exact (1.39), minors 0 and −1, from-scratch cofactor expansion of all 35 3×3 minors; V2 rank vs numpy and sympy on 200 random integer matrices; V6 book matrices of Exs. 1.2–1.5 | 200/200 ranks agree; witness minors nonzero | analytic, symbolic, book-value | |
| C69 Buckingham Π (1.37), D28 | CORE | `pi_groups`, `solve_exponents`, `groups_independent`, `group_*`, presets | V2 sympy D28: nullspace, rank–nullity, Π₁ in span, code groups span null(A), steps 9–11 unit rescaling numerically; V1 exponent solve vs `np.linalg.solve`; V4 properties on all 7 presets × 2 repeating sets (count n − r, zero dimension by code and by pint, independence); V1 Poiseuille collapse Π₁ = 32 Π₂Π₄ | collapse 1e-12; rescaled repeating variables = 1 to 1e-12 | symbolic, analytic, conserved | groups tested as properties, not exact matching (except the book's repeating set) |

## Validation table — coded B/C items (NOTE, ≥ 1 level)
| Concept | fluidpy target | Evidence | Numbers | Label |
|---|---|---|---|---|
| C02/N01 fluid vs solid | `shear_deformation_history` | V1 + V7 | exact closed forms; μ → ∞ no flow; Maxwell recovers τ/G | analytic |
| C03 normal/shear stress | `traction_components` | V1 | τ·n̂ < 1e-12, recomposition < 1e-12 | analytic |
| C04 molecular spacing | `number_density`, `mean_molecular_spacing` | V1, V5 (USSA n) | water 0.31 nm, air 3.4 nm; n rel 3.8e-5 | analytic, benchmark |
| C07/N03 Knudsen, mean free path | `knudsen_number`, `mean_free_path_jennings`, `mean_free_path_air` | V5 Jennings, V2 pint, V7 l ∝ 1/p, V6 | 67.18 nm vs 67.3 nm (−0.18 %) | benchmark |
| C09 mass fraction | `mass_fractions`, `DRY_AIR` | V1, V5 USSA Table 3 | ΣY = 1; M_mix 28.96 | analytic, benchmark |
| C10/C11 Fick, Fourier | `fick_mass_flux`, `fourier_heat_flux` | V1 down-gradient, V2 pint | | analytic |
| C13 μ(T) | `sutherland_viscosity`, `viscosity_power_law`, `water_viscosity` | V5 USSA Table 2; V5 IAPWS R12-08 check point; V7 monotonic | μ rel 3.1e-5; water 298.15 K +0.079 % | benchmark |
| C14 ν, κ, h²/ν | `kinematic_viscosity`, `thermal_diffusivity`, `diffusion_time`, `FLUIDS`, `fluid_properties` | V5 USSA ν; V2 pint; V1 | ν rel 5.3e-5; ν_air/ν_water = 15.0 | benchmark |
| C15/C16/N07 surface tension, Laplace | `surface_tension_water`, `laplace_pressure_jump` | V5 IAPWS Table 1 (calc. column and exp. ± unc.), V1 special cases, V2 pint | max \|Δσ\| = 0.0045 mN/m vs 2-decimal column | benchmark, analytic |
| C17/C21 gauge, layered tank | `gauge_pressure`, `absolute_pressure`, `layered_pressure` | V1 | | analytic |
| C18 isotropy (1.6) | `wedge_pressure_difference`, `wedge_face_forces` | V3 order 1 in dz, V1 | order 1.000 | converged |
| C22 capillary rise | `capillary_rise(_deg)`, `alpha_from_contact_angle` | V1 force balance + Laplace/(1.9) consistency, V7 | residual < 1e-14; R = 1 mm → 14.9 mm | analytic |
| C23/C24 collision time | `mean_molecular_speed`, `collision_time` | V1, V4 Maxwellian mean speed, V5 USSA speed | speed rel 1.9e-5 | benchmark |
| C28–C31, N13 state, h, C_p, C_v | `perfect_gas_state`, `specific_volume`, `enthalpy`, `perfect_gas_internal_energy/enthalpy`, `partial_derivative`, `specific_heat_cp/cv` | V1, V3 order 2 | | analytic, converged |
| C34 Clausius–Duhem | `irreversible_process`, `free_expansion`, `stirred_isochoric_process` | V1 + V4 | Δs = R ln 2 > 0 = ∫δq/T | analytic |
| C37/C48 α | `thermal_expansion_coefficient`, `perfect_gas_expansion_coefficient`, `water_density` | V1 (1/T, linear EOS), V7 sign change near 4 °C | ρ_max at 277.13 K ± 0.05 | analytic |
| C41 van der Waals contrast | `van_der_waals_pressure`, `van_der_waals_internal_energy`, `van_der_waals_constants_per_mass` | V1 Gibbs identity (∂e/∂v)_T = T(∂p/∂T)_v − p numerically; a = b = 0 → perfect gas | rel < 1e-5 | analytic |
| C42/C43/N16 | `cv_from_cp`, `gamma_from_cp`, `cp_from_gamma` | V1 round trips | | analytic |
| C46/C47 | `isentropic_ratios`, `perfect_gas_sound_speed` | V1, V5 | | analytic, benchmark |
| C49 static medium | `atmosphere_from_temperature`, `linear_lapse_pressure`, `standard_atmosphere` | V1, V5 | | benchmark |
| C52 classification | `classify_stability` | V7 | | analytic |
| C56–C59 | `potential_temperature_gradient`, `potential_density` | V1, V2, V4 | | symbolic |
| C60/N21/C61 ocean | `seawater_density_linear`, `seawater_linear_coefficients`, `seawater_density_eos80`, `isentropic_density_gradient`, `ocean_potential_density_gradient` | V1 sign vs −N² (200 random gradients), c → ∞ limit, perfect-gas gradient identity; linear EOS within 0.5 kg/m³ of EOS-80 | ρg/c² = 4.47e-3 kg m⁻⁴ | analytic (EOS-80 itself: Open item O2) |
| C62/C63 isothermal, scale height | `isothermal_pressure/density`, `scale_height` | V2 sympy ODE, V1 p(H) = p0/e, V5 USSA H_p, V6 | H(250 K) = 7318 m | symbolic, benchmark |
| N18 synthetic profile | `synthetic_boundary_layer_profile/column` | V1 layers, N² per layer, θρθ, pressures vs direct integration (1e-7) | | analytic (Finding F1) |
| C72–C76 examples | `poiseuille_pressure_drop`, `blast_energy/radius`, `TAYLOR_K_GAMMA14`, `rayleigh_scattering_ratio`, `pythagoras_phi`, `wavelength_to_rgb` | V1, V5 Taylor K, V7 | blue/red = 5.855; D ∝ t^0.4 exact; Ex. 1.2 constant = 1 | analytic, benchmark (`wavelength_to_rgb`: qualitative by design) |
| R03 °C ↔ K | `celsius_to_kelvin`, `kelvin_to_celsius` | V1 round trip + pint | | analytic |

## Functions used by the notebook and explainers (design Part C — 144 callables)
Every one of the 144 callables is called in at least one test (checked by script against the Part C table; the only
unmatched token is the type name `NamedTuple`). Main mapping:

| Part C rows | Tests |
|---|---|
| 0.1–0.4 `style.setup_notebook/savefig`, `anim.animate/show_animation`, `interact.slider_figure/animate_figure/live`, `embed.show_viz` | `test_machinery_smoke_style_anim_interact_embed` (smoke) |
| 0.5–0.6 units | `test_dimensional_homogeneity_V2_*`, `test_units_V1_*` |
| 1–12 dimensional | `test_dimensional_matrix_*`, `test_buckingham_*`, `test_dimensional_homogeneity_*`, `test_examples_V1_*` |
| 13–35 thermo | `test_perfect_gas_*`, `test_first_law_*`, `test_gibbs_*`, `test_clausius_duhem_*`, `test_specific_heats_*`, `test_sound_speed_*`, `test_thermal_expansion_*`, `test_isentropic_*`, `test_specific_heat_relations_*`, `test_equation_of_state_*` |
| 36–45 statics | `test_hydrostatics_*`, `test_buoyancy_*`, `test_isothermal_atmosphere_*` |
| 46–63 stratification | `test_brunt_vaisala_*`, `test_parcel_*`, `test_stability_classification_*`, `test_adiabatic_lapse_rate_*`, `test_lapse_rate_*`, `test_potential_temperature_*`, `test_ocean_criterion_*` |
| 64–67 diffusion | `test_couette_startup_*`, `test_ftcs_*`, `test_newton_viscosity_V3_*` |
| 68–88 chapter module | `test_continuum_*`, `test_kinetic_pressure_*`, `test_knudsen_*`, `test_newton_viscosity_*`, `test_viscosity_temperature_*`, `test_water_viscosity_*`, `test_fluid_properties_*`, `test_fluid_vs_solid_*`, `test_normal_shear_stress_*`, `test_mass_fraction_*`, `test_surface_tension_*`, `test_laplace_jump_*`, `test_capillary_rise_*`, `test_pressure_isotropy_*`, `test_seawater_eos_*`, `test_synthetic_column_*`, `test_examples_V1_*`, `test_blast_V5_*` |
| 89–90 drawing helpers | `test_drawing_helpers_smoke` (smoke) |

Functions from Part C without evidence: **none**.

## Convergence studies
| scheme | grids | observed order | design order |
|---|---|---|---|
| FTCS vs Couette start-up series (r = 0.4, t = 0.05 h²/ν) | N = 11, 21, 41, 81 (+161 in the figure) | 2.002 (pairwise 2.013, 1.998, 1.999, 1.999) | 2 |
| FTCS vs spreading Gaussian, zero-flux ends (r = 0.25) | N = 61…481 | 2.007 | 2 |
| `shear_stress_profile` / `derivative_2nd_order` incl. end stencils | N = 21…161 | 1.996 | 2 |
| `partial_derivative` central difference | h = 0.1…0.0125 | 2.000 | 2 |
| `process_heat_work` cumulative trapezoid (isentrope) | n = 11…81 | 1.999 | 2 |
| `lapse_rate` (np.gradient edge_order 2, non-uniform grid) | N = 21…161 | within 2 ± 0.25 | 2 |
| nonlinear parcel → linear, ocean mode (error vs ζ0) | ζ0 = 8, 4, 2, 1 m | 2.000 | 2 (O(ζ0²)) |
| nonlinear parcel → linear, dry atmosphere (error vs ζ0) | ζ0 = 200…25 m | 2.002 | 2 |
| from-scratch Euler hydrostatic march | n = 50…400 | 1.002 | 1 |
| from-scratch Euler isentrope dp/dρ = γp/ρ | n = 100…800 | 0.998 | 1 |
| wedge p₂ − p₁ → 0 | dz = 0.1…0.0125 m | 1.000 | 1 |
| FTCS stability limit | r = 0.49 decays; r = 0.51 grows > 10³×; r > 0.5 raises | — | r ≤ ½ |

## Conservation / invariant residuals
| invariant | residual |
|---|---|
| FTCS trapezoid integral, zero-flux ends, 2000 steps | 2.0e-16 relative |
| θ along an integrated dry adiabat (0–10 km) | ptp/θ = 3.8e-12 |
| θ ρ_θ = p_ref/R (USSA 0–11 km) | 4.4e-16 |
| closed reversible 4-leg cycle | Δe = 0, Δs = −8.4e-6 J/(kg K) (trapezoid), q + w = 0 to 1e-9 |
| box face integral of (1.9) | F_x = F_y = 0; F_z/ρgV − 1 = 2.2e-16 |
| N² by three routes (density/(1.29), lapse/(1.32), θ) | 1.9e-7 of max\|N²\| |
| lapse-rate code, 2504-point sweep, both conventions | identical; equals sign of N² everywhere outside the tolerance band of Finding F3 |
| Π-group values under SI → cgs → imperial | < 1e-12 |

## Benchmarks used
| value | our value | source | date verified |
|---|---|---|---|
| k_B = 1.380649e-23 J/K, N_A = 6.02214076e23, R = 8.314462618 J/(mol K) (exact) | identical (R_U/1000 to 1e-10) | NIST CODATA CUU pages | 2026-09-12 |
| USSA-1976 Table 1: T, p, ρ, c at Z = 0, 5, 10, 15, 20, 25, 30, 40, 45, 50 km | max \|ΔT\| 4.6e-4 K; p 1.2e-4, ρ 8.5e-5, c 2.1e-5 rel | PDAS bigtables (from NASA-TM-X-74335) | 2026-09-12 |
| USSA-1976 Table 2: μ, ν, H_p, n, mean particle speed at 0–50 km | 3.1e-5, 5.3e-5, 8.2e-5, 3.8e-5, 1.9e-5 rel | PDAS bigtables | 2026-09-12 |
| USSA-1976 Table 2 constants g0, P0, T0, β = 1.458e-6, γ = 1.40; S = 110.4 K and r0 = 6356.766 km (with NASA errata); R* = 8314.32; Table 3 composition; Table 4 layers | constants identical; layers identical; M0 3.7e-6 rel | NASA-TM-X-74335, NTRS 19770009539 (scanned pages read; errata sheet p. 242) | 2026-09-12 |
| IAPWS R1-76(2014) σ(t), Table 1, 0.01–200 °C | ≤ 0.0045 mN/m vs calculated column; inside experimental uncertainty at every row | iapws.org Surf-H2O-2014.pdf | 2026-09-12 |
| IAPWS R12-08 μ(298.15 K, 998 kg/m³) = 889.735100 µPa s | +0.079 % (Vogel fit; claims ±1 %) | iapws.org viscosity release, Table 4 | 2026-09-13 |
| Jennings mean free path of air 67.3 nm at 300 K, 1 atm | 67.18 nm (−0.18 %) | Tsalikis et al., Aerosol Sci. Technol. 58(8) 2024, doi:10.1080/02786826.2024.2333859 | 2026-09-13 |
| Taylor S(1.4)^−5 = 0.856 | `TAYLOR_K_GAMMA14` = 0.856; 1.032^−5 = 0.8543 (rounding of S) | Díaz, arXiv:2009.05674 | 2026-09-13 |
| AMS dry-adiabatic lapse rate g/c_pd ≈ 9.8 °C/km | 9.7607 K/km (−0.40 %) in the meteorology sign | AMS Glossary (HTML) | 2026-09-13 |

Rejected: Engineering ToolBox c_p (snippet only; C_p is derived from USSA γ = 1.40 instead); CFD-forum Sutherland
constants (typo β = 1.458e-5; replaced by the NASA primary); UNESCO EOS-80 check value (snippet only, primary not
fetched — Open item O2); Kell 1975 table (primary not fetched); PDAS mean free path (hard-sphere formula, not Jennings).

## Numbers from the text (book vs ours) — private values redacted to relative differences
| item | relative difference | comment |
|---|---|---|
| k_B, Avogadro per kmol, R_u, M_w of dry air, R of air, γ, C_p | 0.025 %, 0.014 %, 0.006 %, 0.006 %, 0.020 %, 0 %, 0.070 % | book rounds pre-2019 constants; all ≤ 0.5 % |
| atmospheric pressure (kPa, bar) | 0.025 % | rounding |
| adiabatic lapse rate (book's sign Γ ≡ dT/dz) | 2.4 % | book prints a 1-significant-figure value; ours rounds to it, same sign |
| lab-scale temperature change over 1 m | same at 1 s.f. | |
| scale height at the book's temperature | 0.25 % | |
| isothermal band of the lower atmosphere | USSA max deviation 15.26 % vs the book's rounded band | marginal: holds only at the book's 2-figure rounding |
| mean free path of air | ours 30 % larger | order-of-magnitude statement in the book; Jennings benchmark agrees with ours |
| pipe matrix (1.39), minors, rank, exponents; Exs. 1.2–1.5 matrices, ranks, group counts, Ex. 1.2 exponents; rank example | exact | |
| ocean density change per km (book ≈ value) | not tested | depends on a seawater sound speed that is not a coded constant |

## Figures reproduced with our code
| figure | file (local) | visual verdict |
|---|---|---|
| Fig. 1.9 (a, b) | `outputs/ch01/verify/fig1_9_T_theta_synthetic.png` | T(z) cools along the mixed layer parallel to the grey dry adiabats, reverses in the 800–1000 m inversion and cools slowly above; θ(z) is vertical in the mixed layer, jumps through the inversion and rises above — correct shape, but the mixed layer is coloured *unstable* (Finding F1). |
| Fig. 1.8 idea (parcel) | `parcel_regimes.png` | stable parcel oscillates at period ≈ 642 s with amplitude 5 m and the nonlinear dashed curve lies on the linear one; neutral stays at 5 m; unstable grows like cosh and leaves the frame. |
| Figs. 1.2–1.3 idea (momentum diffusion) | `couette_ftcs_vs_series.png` | FTCS dots sit on the series curves from the early boundary-layer shape to the linear Couette profile; error line parallel to the slope-2 guide over 1.2 decades of Δy. |
| §1.4 continuum | `continuum_noise_and_drift.png` | sampled relative noise on the (nL³)^−1/2 line from ≈ 1 at 3 nm to 6e-6 at 10 µm; noiseless box density shows the plateau at ρ(x0) for L ≪ L_flow, falls to the mean at L = L_flow and rings around it (sinc) beyond. |
| lapse-rate criterion (both conventions) | `lapse_rate_sweep.png` | N² is linear in dT/dz, crosses zero exactly at the dashed Γa = −9.76 K/km, and the stability code steps from −1 to +1 at the same place; the top axis shows the meteorology Γ with the same crossing at +9.76. |
| USSA-1976 | `ussa1976_vs_table.png` | our T(H) broken line passes through every PDAS point (kinks at 11, 20, 32, 47 km); log p(H) points on the curve. |
| §1.11 pipe collapse | `pi_collapse_pipe.png` | 80 random laminar cases fall on Π₁ = 32 Π₂Π₄ over four decades. |
| §1.6 σ(T) | `iapws_surface_tension.png` | monotone decrease from 75.6 mN/m at 0 °C to 0 at 374 °C through all IAPWS points. |

## Derivations re-derived (curation §4b, design Part F)
| D | ★ | test | intermediate lines checked | result |
|---|---|---|---|---|
| D14 isentropic law | ★★ | `test_isentropic_law_V2_derivation_D14` | steps 2, 3, 6, 7, 10 | holds |
| D18 parcel equation | ★★ | `test_parcel_equation_V2_derivation_D18` | steps 3, 6, 7–9, 10, 12 (+ code parity of N²) | holds |
| D19 adiabatic lapse rate | ★★★ | `test_adiabatic_lapse_rate_V2_derivation_D19` | steps 4, 7, 10, 12, 15, 16 for arbitrary G(T, p); α = 1/T and C_p for the perfect-gas G | holds |
| D28 Buckingham Π | ★★★ | `test_buckingham_V2_derivation_D28_nullspace_rank_nullity` | steps 3, 4, 5, 7; steps 9–11 numerically with `rescale_units` | holds |
| D34 kinetic pressure | ★★ | `test_kinetic_pressure_V2_derivation_D34` | steps 5, 7, 8 | holds |
| D36 N² from θ | ★★ | `test_potential_temperature_V2_derivation_D36` | steps 2, 3, 5, 6–7, 8 | holds |
| D10, D20, D35, D37 (★) | ★ | `test_gibbs_V2_*`, `test_potential_temperature_V2_derivation_D20_and_eq_1_32`, `test_continuum_V2_*`, `test_buoyancy_V2_*` | as listed in the tests | holds |

No wrong intermediate line was found in Part F. The Part F numbers checked (696.4 J/(kg K), 287.06, 1.225 kg/m³,
351.3 K, 1.641, 9.57e-5 s⁻², 642 s, −9.76 K/km, −0.10 K/km for water, 304.8 K, 3.83e-4 s⁻², 2.55e10 molecules,
6.3e-6, 9.807 N, Π₁ = 10) all reproduce.

## Deviations & justifications
| `# DEVIATION` in code | assessment |
|---|---|
| `seawater_density_eos80`: salinity in g/kg used as practical salinity (≈ 0.5 %) | acceptable for teaching; stated in docstring |
| `FLUIDS` solute diffusivity of glycerine/honey/oil from Stokes–Einstein scaling | order of magnitude only; used for E2 presets, tested only for positivity and ν = μ/ρ |
Also documented (not marked DEVIATION): `standard_atmosphere` uses CODATA-based `R_AIR` instead of USSA R* (1.7e-5
relative; table agreement 1.2e-4 in p at 50 km).

## Findings (for the implementer/designer; none blocks this verdict)
- **F1 · Should fix — `synthetic_boundary_layer_profile` default mixed layer is not near-neutral.** Default
  `mixed_dT_dz = −9.8e-3` K/m is steeper than Γa = −9.7607e-3 K/m, so `synthetic_boundary_layer_column` labels the whole
  0–800 m layer **unstable** (N² = −1.35e-6 s⁻²) and `lapse_rate_stability` prints "−9.80 < −9.76 K/km". The docstring
  calls it "marginally steeper … near-neutral" and Fig. 1.9 intends a neutral layer; the C55 slider and E4 inversion
  preset will show a red unstable layer. Suggested: default `mixed_dT_dz = adiabatic_lapse_rate()` (exactly neutral) or a
  value slightly *less* steep (e.g. −9.7e-3, marginally stable). Not encoded in any test.
- **F2 · Nice — docstrings still say "Validation (planned) … Label: pending".** With this report the labels can be
  filled (CLAUDE.md rule 5), e.g. from the Label column above.
- **F3 · Nice — neutral tolerances of `lapse_rate_stability` (|margin| ≤ 1e-9 K/m) and `classify_stability`
  (|N²| ≤ 1e-12 s⁻², i.e. |margin| ≤ ≈2.7e-11 K/m) differ.** For 2.7e-11 < |dT/dz − Γa| ≤ 1e-9 K/m one helper says
  neutral and the other stable/unstable. Physically irrelevant; the V4 sweep test excludes that band explicitly.
- **F4 · Nice — `seawater_density_eos80` docstring quotes UNESCO check values that could not be verified** from a primary
  source today (O2).

## Open items
- O1 (qualitative, by design): `wavelength_to_rgb` is a visual aid; smoke-tested only (blue-dominant at 450 nm, red at
  650 nm, black outside 380–780 nm).
- O2 (missing benchmark): `seawater_density_eos80` has self-consistency evidence only (linear EOS agrees within
  0.5 kg/m³; monotone in S). The UNESCO Tech. Pap. Mar. Sci. 44 check values were not fetched; add them to
  `reference/ch01/` when the primary PDF is obtained.
- O3 (book statement): the book's isothermal-band claim is met by USSA-1976 only at the book's rounding (15.26 % vs a
  2-figure band); the notebook should present it as an order-of-magnitude statement.

## Verdict: PASS
