# Chapter 7 verification — Gravity Waves                     2026-09-24, commit 6b69363 + working tree (loop 1)

Verifier: math-verifier. Suite: `tests/test_ch07.py` (130 test functions = 130 collected items), figures/metrics:
`tests/ch07_verify_figures.py` → `outputs/ch07/verify/` (git-ignored), cited data: `reference/ch07/` (`make_refs.py`,
`benchmarks.json`, `SOURCES.md`). `tests/conftest.py` gained a `pytest_configure` that registers the `slow` marker (the
script smoke test). Nothing in `fluidpy/` or `scripts/` was edited by the verifier.

## Environment
Python 3.11.5 · numpy 2.4.6 · scipy 1.17.1 · sympy 1.14.0 · pint 0.25.3 · matplotlib 3.11.2 (Windows, `.venv`).

Runs: `pytest tests/test_ch07.py -q -p no:cacheprovider` → **130 passed in 125 s** (with the private book file; without
it the three V6 tests skip; `-m "not slow"` drops the 44 s script test → 129 in ≈ 81 s). Full suite `pytest -q` →
**880 passed in 579 s** (ch01–ch06 + machinery + the 130 ch07 items; one warning, from `core.waves` — O9). Re-run of the final `tests/test_ch07.py` after the last (V6 comment/assert) edits: 130 passed in 160 s under load. All 11 `scripts/ch07_*.py` exit 0 headless (asserted by `test_scripts_V1_every_ch07_script_runs`,
44 s in total). `tests/ch07_verify_figures.py` runs and writes 11 PNGs.

Collected items per evidence tag (tag of the `def` line): **V1 58 · V2 36 · V3 10 · V4 5 · V5 5 · V6 3 · V7 13** (= 130).
Many tests carry a second level inside (e.g. a V1 test that also measures an order); the table below lists every level
exercised per item.

## Loop 1 — what the first run found (all test-side; no physics defect)
First run: 119 of 129 passed. All ten failures were in the tests, fixed before this report:
ch04's `linear_wave_surface` overflows at kH = 900 (plain cosh/sinh — a ch04 limitation, see O8; the parity is now
restricted to kH < 300 and ch07's overflow-safe fields are asserted finite there); the shallow pressure limit needed the
(kH)²/2 correction in its tolerance; a sympy check of D23 compared two different `Subs` representations of the same
derivative (rewritten with the concrete ω(k, H(x))); D25's quadratic was solved with a positive symbol (lost the negative
root that Vieta's step needs); the KdV right-hand-side check had the soliton's tail on the periodic seam; `dyed_line`
releases particles *on* the line (their mean depth is lower by aS), so it is now compared with
`particle_path(start="mean")`; a size-1 array → float; a missing sympy symbol; two tolerances I had set tighter than the
design's own check numbers (ray ω-drift at machine precision instead of the solver's rtol 1e-10; the chirped-train
(7.74) residual at 1e-10 instead of the design's 1e-9, whose nested-FD round-off bound is 1.8e-9). No tolerance was
loosened below a documented or derived value.

## Discrimination proofs (wrong variants planted in a scratch copy; the repository untouched)
Method: a scratch copy of `fluidpy/`, `tools/`, `tests/test_ch07.py` (+ private JSON, `reference/ch07/`); one wrong variant planted at a time; `pytest -x -m "not slow"`; the first failing test recorded. **35/35 caught.**

| # | wrong variant | where | first test that fails |
|---|---|---|---|
| 1 | ω = √(gk)·tanh kH instead of √(gk tanh kH) | `core.waves.omega_gravity` | `test_linear_evolve_V1_single_modes_and_superposition` |
| 2 | c_g missing the ½ of (7.69) | `core.waves.group_velocity` | `test_dispersion_state_V1_explainer_numbers` |
| 3 | σk²/ρ instead of σk/ρ inside c² | `core.waves.phase_speed` | `test_capillary_minimum_V1_numerical_minimisation` |
| 4 | c_min exponent ½ instead of ¼ | `core.waves.capillary_minimum` | `test_capillary_minimum_V1_numerical_minimisation` |
| 5 | ε² = (ρ₂ − ρ₁)/ρ₂ (analysis' wrong variant) | `core.waves.interface_omega` | `test_interface_fields_V1_residuals_vortex_sheet_and_parity` |
| 6 | g′ with ρ₁ in the book form (7.117) | `core.waves.reduced_gravity_book` | `test_two_layer_modes_V1_amplitude_ratios_and_pressure` |
| 7 | printed (7.145) for every sign of k | `core.waves.internal_wave_velocities` | `test_internal_wave_velocities_V1_gradient_parity_both_signs` |
| 8 | one-way evolution with the wrong sign | `core.waves.linear_evolve` | `test_linear_evolve_V1_single_modes_and_superposition` |
| 9 | ray dk/dt = +∂ω/∂x | `core.waves.ray_trace` | `test_ray_trace_V4_frequency_conserved_while_k_grows` |
| 10 | (7.66) with the printed ½Δω x by default | `core.waves.beat_wave` | `test_beat_wave_V1_equals_sum_printed_fails_and_chord_to_tangent` |
| 11 | deep error 1 − tanh kH (no square root) | `core.waves.depth_regime` | `test_depth_regime_V1_book_thresholds_and_errors` |
| 12 | Doppler ω − U·K | `core.waves.doppler_frequency` | `test_doppler_V1_probe_frequency_of_a_translated_pattern` |
| 13 | ψ with the GFD (Ch. 13) sign | `ch07.wave_fields` | `test_wave_fields_V1_closed_forms_and_streamfunction` |
| 14 | ξ with + sign (anticlockwise orbits) | `ch07.orbit_linear` | `test_orbit_V1_ellipses_clockwise_and_constant_foci` |
| 15 | Stokes drift decaying like e^{kz} | `ch07.stokes_drift` | `test_orbit_state_V1_explainer_numbers` |
| 16 | Bélanger with 4Fr₁² instead of 8Fr₁² | `ch07.hydraulic_jump` | `test_hydraulic_jump_V1_momentum_mass_and_numbers` |
| 17 | jump energy change with the dropped minus (text extraction) | `ch07.hydraulic_jump` | `test_hydraulic_jump_V1_momentum_mass_and_numbers` |
| 18 | KdV nonlinear coefficient 3 instead of 3/2 (solver) | `ch07.kdv_solve` | `test_kdv_solve_V1_soliton_translates_at_c0_1_plus_a_over_2H` |
| 19 | KdV dispersive coefficient 1/3 instead of 1/6 (solver) | `ch07.kdv_solve` | `test_kdv_solve_V1_soliton_translates_at_c0_1_plus_a_over_2H` |
| 20 | energy flux at c instead of c/2(1 + …) | `ch07.energy_flux` | `test_energy_flux_V1_quadrature_equals_closed_form` |
| 21 | seiche tanh without π | `ch07.seiche_modes` | `test_seiche_V1_walls_modes_and_dispersion` |
| 22 | exact kinematic residual +η_xφ_x | `ch07.free_surface_residuals` | `test_free_surface_V1_exact_residuals_equal_independent_closed_forms` |
| 23 | interfacial E_p over λ/2 (printed middle form, ⅛) | `ch07.interface_energy` | `test_interface_energy_V1_quarter_from_direct_integration` |
| 24 | ρ′ in phase with w (no factor i) | `ch07.internal_wave_fields` | `test_internal_wave_fields_V3_residuals_second_order` |
| 25 | two-layer b with swapped exponentials | `ch07._two_layer_constants` | `test_two_layer_residuals_V1_numeric_and_printed_variant` |
| 26 | solitary width √(3a/H³) | `ch07.solitary_wave` | `test_kdv_solve_V1_soliton_translates_at_c0_1_plus_a_over_2H` |
| 27 | Stokes third harmonic ⅛ instead of 3/8 | `ch07.stokes_wave_profile` | `test_stokes_wave_profile_V1_harmonics_and_permanence` |
| 28 | pressure factor cosh/sinh instead of cosh/cosh | `ch07.pressure_response` | `test_pressure_response_V1_surface_bottom_and_limits` |
| 29 | standing ψ without the minus of ψ₂ (cos kx) | `ch07.standing_wave_fields` | `test_standing_wave_V1_sum_of_opposite_waves` |
| 30 | internal flux F_z with the wrong sign | `ch07.internal_wave_energy` | `test_internal_wave_energy_V1_closed_forms_equal_averages` |
| 31 | E_k = ½ρga² (double) | `ch07.wave_energy` | `test_wave_energy_V1_quadrature_equals_closed_form` |
| 32 | rigid lid with h₁, h₂ swapped | `ch07.two_layer_rigid_lid_omega` | `test_two_layer_V7_limits_and_reduced_gravity` |
| 33 | cnoidal trough with the wrong sign | `ch07.cnoidal_wave` | `test_cnoidal_V1_form_and_soliton_limit` |
| 34 | PE interface limit without ⟨cos²⟩ = ½ | `ch07.internal_pe_interface_limit` | `test_internal_pe_interface_limit_V3_epsilon_order_one` |
| 35 | Stokes speed with ½(ka)² inside c² | `ch07.stokes_wave_speed` | `test_stokes_wave_V5_speed_and_limiting_steepness` |

## Validation table — A items (CORE, ≥ 2 independent levels, one of V1/V2/V3/V5)
| Concept / Eq. | Tier | fluidpy target | Evidence (V-levels) | Numbers (error, order, residual) | Label | Notes |
|---|---|---|---|---|---|---|
| C01 sinusoid, k, ω, c, K (7.1)–(7.9) | CORE | `W.sinusoid`, `crest_positions`, `wave_parameters`, `plane_wave`, `phase_velocity_vector`, `trace_velocities`, `doppler_frequency` | V1 · V2 (D01) · V7 | crests at ω/k to 1e-12; round trips 1e-12; 1/c² = Σ1/c_i² for 50 random K (1e-12); probe frequency ω + Uk to 1e-6 (Hilbert phase); rotation invariance 1e-13 | analytic, symbolic | the "c = (ω/k, ω/l, ω/m)" reading is shown to overstate c |
| C02 linearised free-surface problem (7.11)–(7.21) | CORE | `free_surface_residuals`, `free_surface_residual_scan`, `surface_normal`, `surface_velocity`, `linear_bernoulli_pressure` | V1 · V2 (D02, D03, D04) · V3 | linear conditions (7.18), (7.21), (7.55) ≤ 1e-14 aω / 1e-13 ag; exact residuals equal our independent closed forms to 1e-12; **exact-condition residuals: absolute slope 2.010 (kin.) / 2.003 (dyn.) at kH = 1, 2.014 / 2.022 deep; relative (÷ aω, ÷ ag) slope 1.010 / 1.003 and 1.014 / 1.022**; parity with ch04 `kinematic_bc_residual(z − η)` 1e-8 aω | analytic, symbolic, converged | implementer's note 2 confirmed: the dropped terms are O(ka) *relative* (slope 1), O((ka)²) absolute at fixed k (slope 2) — analysis §6 V25's "(ka)² × aω" is wrong (asserted: relative slope ≠ 2); passing the elevation as the level set (R1 trap) fails |
| C03 dispersion relation (7.28) | CORE | `W.omega_gravity`, `wavenumber_from_omega`, `period_from_wavelength`, `wavelength_from_period`, `surface_wave_sympy` | V1 · V2 (D05, D06) · V5 · V7 | ω² = gk tanh kH 1e-14; inverse ∘ forward per point 1e-12 over ω ∈ [1e-4, 1e3] rad/s × H ∈ [1e-3, 1e4] m and ∞; parity with ch04's test field 1e-12; Fenton & McKee (1990): max λ error 1.658 % ≤ 1.7 % (rounds to 1.7), optimal ν 1.483 vs 1.49 (0.46 %); Guo (2002, β = 2.4908): 0.753 % vs 0.75 % | analytic, symbolic, benchmark | see "Benchmarks" for the 2006 note's 1.5 % / 0.7 % (not reproducible; implementer's note 3) |
| C04 phase speed (7.29), deep/shallow (7.45)–(7.52) | CORE | `W.phase_speed`, `depth_regime`, `pressure_response`, `wave_fields`, `dispersion_state`, `wave_regime_label` | V1 · V2 (D07, D09 step 3; potential relations) · V3 · V7 | deep 1e-15 at kH = 40; shallow √(gH)(1 − (kH)²/6) to (kH)⁴; dc/dk < 0 everywhere; c < √(gH); deep error at kH = 2: 1.815 %; shallow error at H = 0.07λ: 3.039 %; e^{−π} at −λ/2 1e-14; 5-point ∇²φ residual order 2.00 (± 0.15); no overflow at kH = 500 | analytic, symbolic, converged | |
| C05 particle orbits (7.32)–(7.37) | CORE | `orbit_linear`, `orbit_semi_axes`, `particle_path`, `orbit_state` | V1 · V2 (D10, D11) | ellipse (7.36) = 1 to 1e-12 at 4 depths; foci a/sinh kH at every depth (spread 1e-13); signed area < 0 (clockwise); "linear" path = (7.35) to 1e-9 a; "exact" path = `core.kinematics.pathline` (ch03) to 1e-8 a; ψ (7.37) → (7.27) through ch04's `velocity_from_streamfunction_2d` 1e-8 | analytic, symbolic | |
| C06 wave energy and flux (7.38)–(7.44) | CORE | `wave_energy`, `energy_flux`, `W.wave_energy_density` | V1 · V2 (D12, D13, D14) · V4 · V3 | quad = closed 1e-8 for kH ∈ {0.1, 1, 10, ∞}; E_k = E_p; F quad = (7.44) 1e-8 at kH = 0.35, 2.8, 14, ∞; F = E c_g 1e-13; ∫η² of a one-way packet conserved to 1e-12 (V4); its centroid speed → c_g with error 1.9e-3, 4.7e-4, 1.2e-4 (order 2.0) | analytic, symbolic, conserved | |
| C07 capillary–gravity waves (7.53)–(7.60) | CORE | `W.omega_capillary_gravity`, `phase_speed(σ)`, `capillary_minimum`, `min_group_velocity`, `curvature`, `capillary_surface_pressure`, `capillary_state` | V1 · V2 (D15, D16) · V5 · V7 | minimum = `minimize_scalar` (c 1e-12) for water, an oil-like and a mercury-like fluid; c_g = c at k_m 1e-13; c_g,min at σk²/ρg = 2/√3 − 1 (1e-14); IAPWS σ(20 °C), ρ = 998.2: 0.2312 m/s at 1.712 cm vs published 0.23 / 1.7; circle curvature 1/R 1e-14; parity with ch04 Laplace jump 1e-12 | analytic, symbolic, benchmark | |
| C08 standing waves and seiches (7.61)–(7.65) | CORE | `standing_wave_fields`, `seiche_modes`, `basin_modes`, `seiche_state` | V1 · V2 (D17) · V3 · V7 | standing = sum of ± waves 1e-13; (7.62)–(7.63) 1e-13; fixed nodes; zero mean flux (Ex. 7.8); u = 0 on both walls; **our 2nd-order FD sloshing (Steklov) eigenproblem converges to (7.65) at order 1.98 (pairwise 1.86, 1.97, 1.99; max error 2.4e-4 at Δx = L/160)**; shallow-limit error ∝ (kH)² (order 2.0) | analytic, symbolic, converged | the FD check is independent of (7.28): it solves Laplace with the free-surface condition directly |
| C09 group velocity (7.66)–(7.71) | CORE | `W.group_velocity`, `group_velocity_numeric`, `beat_wave`, `gaussian_packet`, `linear_evolve`, `envelope`, `pond_ripples`, `packet_state`, `viscous_decay` | V1 · V2 (D19, D20, D21, a-D68) · V3 · V7 | analytic = complex-step c_g 1e-12 (gravity, finite depth, capillary–gravity, Rossby-type, internal); central fallback order 4.00; beats = sum of cosines 1e-13, the printed ½Δω x fails; chord → tangent order 2.0; order-1 packet rigid at c_g 1e-14; order-2 packet = FFT evolution for quadratic ω 1e-12; pond cosine sum → Cauchy–Poisson integral order 2.0 | analytic, symbolic, converged | |
| C10 kinematic waves and rays (7.72)–(7.79) | CORE | `local_wavenumber_frequency`, `crest_conservation_residual`, `W.ray_trace`, `snell_ray_plane_beach`, `refraction_state` | V1 · V2 (D22, D23, D24) · V3 · V4 | local k order 4.00; (7.74) 1.3e-10 (≤ 1e-9); (7.79) 0 and local dispersion 1e-12 over a slope, (7.75) ≠ 0 there (homogeneous media only); homogeneous ray straight at c_g 1e-8; ω drift 1.9e-9 over a 1:50 beach to H = 2 cm while k grows 25×; k_y conserved 2e-12; ray ≡ Snell closed form 4.5e-11 (relative); T = 8 s, 30° at 20 m → 11.27° at 2 m | analytic, symbolic, converged, conserved | |
| C11 hydraulic jump (7.80)–(7.81); KdV (7.87)–(7.88) | CORE | `hydraulic_jump`, `jump_momentum_residual`, `jump_state`, `kdv_rhs`, `kdv_solve`, `kdv_invariants`, `kdv_linear_phase_speed`, `solitary_wave`, `cnoidal_wave`, `kdv_residual_sympy`, `simple_wave_evolve`, `nonlinear_wavelet_speed`, `ursell_number` | V1 · V2 (D25, D26, Ex. 7.15) · V3 · V4 | momentum residual ≤ 1e-13; mass; Fr₂ ≤ 1; E₂ − E₁ = −g(H₂ − H₁)³/4H₁H₂ = direct energies 1e-10; dE ≤ 0 and decreasing for Fr₁ ≥ 1, > 0 (forbidden) below; bore speed = ch04 `bore_speed` 1e-13; KdV soliton after 30 m: shape error 2.2e-7·a; IF-RK4 order 4.00 in dt; spectral in N (×470 from 64 → 128); invariants: mass 2e-16, momentum 7e-10, Hamiltonian 9e-10 (drift ∝ dt^4.9); linear KdV speed error ∝ (kH)⁴ (order 4.0); cnoidal = published formulas 1e-13, m → 1 (trough datum) → soliton 3e-9; breaking time = 1/max(−∂c/∂x) 1e-6 | analytic, symbolic, converged, conserved | implementer's note 5 confirmed (the mean datum converges only like 1/K(m)); note 6: the default dt drift is ~1e-9 relative for a splitting hump (≈1e-11 was for the soliton) |
| C12 Stokes waves and drift (7.82)–(7.86) | CORE | `stokes_drift`, `stokes_drift_numeric`, `eulerian_mean_u`, `particle_path(model="taylor1")`, `dyed_line`, `stokes_wave_profile`, `stokes_wave_speed`, `stokes_expansion_sympy`, `STOKES_LIMIT_STEEPNESS` | V1 · V2 (D27, a-D40, a-D42) · V3 · V5 | (7.85)/(7.86) closed forms 1e-13; Eulerian mean ≤ 1e-10 aω; numeric drift from exact path lines → (7.86) with relative error 0.61 %, 1.21 %, 2.39 %, 4.62 % at ka = 0.01…0.08 (**order 0.976**; pairwise 0.99, 0.98, 0.95); rtol halving 1e-7; independent third-order expansion: η₂ = ½, η₃ = 3/8, ω² = 1 + ε², φ one harmonic (a₁₃ = −1/8 = the coded β + ½); published speed (1 + ½(ka)²)√(g/k) agrees to O((ka)⁴) (order 4.0); z = −λ/4 fraction 4.32 % ("about 4 %") | analytic, symbolic, converged, benchmark | implementer's note 1 confirmed independently: the consistent third-order solution gives γ = 1, the literal truncated set-up γ = 3/8 |
| C13 interfacial waves (7.89)–(7.96) | CORE | `W.interface_omega`, `interface_fields`, `interface_residuals`, `interface_energy`, `W.eps2_density` | V1 · V2 (D28) · V7 | residuals of (7.90)–(7.94) ≤ 1e-12 scale; A = −B = iωa/k; u₁ = −u₂ at z = 0; sheet strength 2ωa cos θ = ch05 `vortex_sheet_strength` 1e-13; E_k = E_p = ¼Δρga² by dblquad 1e-8 (**¼ from direct integration; the printed middle form over ∫₀^{λ/2} gives ⅛**, implementer's note 4); ρ₁ → 0 → √(gk) 1e-15; ρ₁ > ρ₂ → NaN + RuntimeWarning; published two-fluid form (σ = 0) identical | analytic, symbolic | |
| C14 two-layer modes (7.97)–(7.119) | CORE | `two_layer_modes`, `two_layer_residuals`, `two_layer_sympy`, `two_layer_rigid_lid_omega`, `two_layer_state`, `W.two_layer_free_surface_omega`, `W.two_layer_long_wave_speed`, `W.reduced_gravity_book` | V1 · V2 (D29, D30, D31) · V7 | every D29/D30/D31 line re-derived (see Derivations); residual of (7.101) = (ag²k/ω²)·(7.110) exactly; numeric residuals ≤ 1e-12 scale for both modes; printed (7.105) fails (Laplace and interface residuals ≠ 0); η/ζ = (7.114) 1e-9; long-wave c = √(g′H) and η/ζ → −Δρ/ρ₁ (1e-3 at kH = 4e-4); kH = 60 → (7.95) 1e-12; g′_lower/g′_upper = ρ₁/ρ₂; rigid lid → (7.95) and g(ρ₂ − ρ₁)h₁h₂/(ρ₁h₂ + ρ₂h₁) | analytic, symbolic | |
| C15 internal-wave dispersion (7.120)–(7.139) | CORE | `boussinesq_linear_sympy`, `W.internal_wave_omega`, `W.beam_angle`, `internal_wave_fields`, `w_equation_residual`, `layered_flow_check` | V1 · V2 (D32, D33, D34) · V3 · V7 | (7.134) eliminated exactly for arbitrary N(z); FD residuals of (4.10), (7.128), (7.130), (7.131) order 2.00; w-equation residual order 2.0, a 10 % wrong ω leaves 100× more; K·u = 0 to 1e-14; ρ′ = N²ρ₀ζ/g 1e-13; m = 0 → ω = N = ch01 parcel frequency; ω(sK) = ω(K) | analytic, symbolic, converged | |
| C16 c ⟂ c_g, beams, F = c_g E (7.140)–(7.159) | CORE | `W.internal_wave_velocities`, `W.group_velocity_vector`, `internal_wave_energy`, `internal_energy_budget_residual`, `internal_pe_interface_limit`, `st_andrews_cross`, `W.linear_evolve_2d`, `internal_wave_state` | V1 · V2 (D35–D37) · V3 · V4 · V7 | c_g = ∇_Kω (complex step) 1e-12 for all four sign quadrants, printed (7.145) fails for k < 0; c·c_g ≤ 1e-15; E_k = E_p, F = c_g E 1e-13, closed forms = period averages 1e-12; (7.147) pointwise residual order 2.00 (5.8e-5 W/m³ at h = 2.5e-3 vs terms O(7)); ε → 0 limit → ¼Δρga² at order 0.97–1.00; beam axes at arccos(ω/N) (1e-10) and field maximum within 1°; 2-D packet centroid moves at c_g within 0.6 % | analytic, symbolic, converged, conserved | St Andrew's cross is our labelled illustration (Mowbray & Rarity 1967 pedigree) |

## Validation table — coded B/C items (NOTE, ≥ 1 level) and recaps
| Item | fluidpy target | Test (level) | Result |
|---|---|---|---|
| N02 Fourier superposition | `linear_evolve` | `test_linear_evolve_V1_single_modes_and_superposition` (V1), `…_V7_nondispersive_translation` (V7) | 1e-13 |
| N03–N05, D01 | `sinusoid`, `crest_positions`, `wave_parameters` | V1, V2 | 1e-12 |
| N06–N10, N167 | `plane_wave`, `phase_velocity_vector`, `trace_velocities` | V1, V7 | 1e-12 |
| N11 Doppler | `doppler_frequency` | V1 | 1e-6 (measured), exact form |
| R02–R06, N13–N17 | `free_surface_residuals`, `surface_normal`, `surface_velocity`, `linear_bernoulli_pressure` | V1, V2, V3 | as C02 |
| N19–N25 | `surface_wave_sympy` | V2 (all residuals 0, wrong bottom sign ≠ 0) | exact |
| N26, N42, N46 | `pressure_response` | V1 | 1e-14 |
| N27–N30, N40, N41, N44, N45, N169–N171 | `particle_path`, `orbit_linear`, `orbit_semi_axes`, `wave_fields` | V1, V7 | as C05 |
| N31–N36 | `wave_energy`, `energy_flux` | V1, V2 | 1e-8 |
| N37–N39, N43 | `depth_regime`, `wave_regime_label`, `wavelength_from_period` | V1, V6 | 1.815 %, 3.039 % |
| N47, N48, N174, N175, N184 | `ray_trace`, `snell_ray_plane_beach` | V1, V4 | 4.5e-11, 1.9e-9 |
| R09, N50–N57, N176 | `curvature`, `capillary_surface_pressure`, `capillary_minimum`, `phase_speed(g = 0)` | V1, V5, V7 | as C07 |
| N58–N63, N177, N178 | `standing_wave_fields`, `seiche_modes`, `basin_modes` | V1, V3 | as C08 |
| N64–N73, N179–N182 | `beat_wave`, `gaussian_packet`, `group_velocity`, `min_group_velocity`, `pond_ripples`, `viscous_decay` | V1, V2, V3, V7 | as C09 |
| N75–N80 | `local_wavenumber_frequency`, `crest_conservation_residual` | V1, V3 | order 4.00 |
| N81, N185 (our simple wave, labelled) | `nonlinear_wavelet_speed`, `simple_wave_evolve` | V1 | t_b 1e-6; both speeds = c₀(1 + 3η/2H) + O(η²) (order 2.0) |
| N82–N85, R10, N186 | `hydraulic_jump`, `jump_state` | V1, V2, V4 | as C11 |
| N86–N92, N187, N188 | `stokes_wave_profile`, `stokes_wave_speed`, `stokes_drift*`, `eulerian_mean_u`, `dyed_line` | V1, V2, V3, V5 | as C12 |
| N93–N97, N189 | `kdv_*`, `solitary_wave`, `cnoidal_wave`, `ursell_number` | V1, V2, V3, V4 | as C11 |
| N99 complex notation | `real_field` | V1 | ½Re(AB*) 1e-13 |
| N105, N107, N190 | `interface_energy`, `interface_fields` | V1 | ¼Δρga² 1e-8; sheet 1e-13 |
| N109–N131, N193, N194 | `two_layer_*`, `reduced_gravity_book` | V1, V2, V7 | as C14 |
| R13, R18, N132–N147, N196 | `boussinesq_linear_sympy`, `internal_wave_fields`, `layered_flow_check`, `ST.parcel_displacement` parity | V1, V2, V3, V7 | as C15 |
| N148–N166, N195, N197–N199 | `internal_wave_velocities`, `internal_wave_energy`, `internal_energy_budget_residual`, `internal_pe_interface_limit`, `st_andrews_cross`, `linear_evolve_2d` | V1, V3, V4, V7 | as C16 |
| dimensional homogeneity (V49) | formulas of (7.28), (7.56), (7.58), (7.42), (7.44), (7.85), (7.117), (7.158), (7.54) | `test_dispersion_family_V2_dimensional_homogeneity` (pint) | all consistent |
| units invariance | 16 dimensional functions | `test_wave_functions_V7_change_of_units` (m → cm, s → ds, kg → g) | every output rescales by its own dimension to 1e-13 |

RECAP items R01–R20 were tested in their own chapters; the ones this chapter re-uses through a parity (R02 ch06 stencil,
R04 ch04 kinematic residual, R06 ch04 Bernoulli, R08 ch03 path lines, R09 ch04 Laplace jump, R10 ch04 bore speed, R18
ch01 parcel) are exercised again above.

## Derivations (curation §4b, design Part F) — every ★★ and ★★★ re-derived with sympy
| D | ★ | test | what is checked, line by line |
|---|---|---|---|
| D01 | ★ | `test_phase_speed_V2_derivation_crest_condition` | x_crest(t), n drops out, c = λν, cos = 1 on crests |
| D02 | ★★ | `test_kinematic_condition_V2_derivation` | steps 1–8: ∇f, n (7.14), U_s (7.15), ×\|∇f\|, the dot products, (7.16) and D(z − η)/Dt |
| D03 | ★★ | `test_linearised_kinematic_V2_derivation` | O(a) part of (7.16) at z = η is (7.18) (= 0); the O(a²) part is exactly ηφ_zz − η_xφ_x at z = 0 (steps 4–6), = ka·aω × f(θ, kH) |
| D04 | ★★ | `test_linearised_dynamic_V2_derivation` | O(a) part of Bernoulli at z = η is (7.21) and vanishes with (7.28); O(a²) part = ηφ_tz + ½\|∇φ\|² = ka·(aω²/k) × f(θ, kH) |
| D05 | ★★ | `test_potential_V2_derivation` | steps 1–11: Laplace → f″ − k²f, dsolve, B = Ae^{−2kH}, A, the cosh regrouping, 2Ae^{−kH} = aω/(k sinh kH), (7.26); wrong-sign B fails |
| D06 | ★★ | `test_dispersion_V2_derivation` | φ_t, the coth line, ω² = gk tanh kH, positive root, T–λ form = `period_from_wavelength` |
| D07, D09 | ★ | `test_phase_speed_V2_derivation_longer_is_faster` | d/dλ[λ tanh(2πH/λ)] = (½ sinh 2x − x)/cosh²x > 0; c/√(gH) = 1 − (kH)²/6 |
| D10, D11 | ★★ | `test_orbits_V2_derivation` | ∫(7.34) → (7.35), zero-mean constants, ellipse, B/A = tanh, focal distance |
| D12, D13 | ★★ | `test_wave_energy_V2_derivation` | ⟨cos²⟩, ∫cosh², ∫sinh², steps 7–10 → ½ρg⟨η²⟩; D13's ∫_{−H}^{η} − ∫_{−H}^0 = η²/2 |
| D14 | ★★ | `test_energy_flux_V2_derivation` | steps 5–10 each, then (7.44) = E dω/dk |
| D15 | ★★ | `test_capillary_V2_derivation_tension_condition` | small-slope curvature, η_xx = −k²η, (7.55) → (7.56), (7.57), p > p_a under a crest |
| D16 | ★★ | `test_capillary_minimum_V2_derivation` | k_m, equal terms, λ_m, c_min², c_min, second derivative, c_g = c at k_m |
| D17 | ★★ | `test_standing_wave_V2_derivation` | sum-to-product, ψ₂ sign, (7.62), (7.63); the missing minus gives cos kx (fails) |
| D19 | ★ | `test_beats_V2_derivation` | (7.66) with ½Δω t; the printed ½Δω x is not the sum and is frozen in time |
| D20 | ★★★ | `test_packet_envelope_V2_derivation` | step 6 phase split; steps 4, 8–10 envelope a(x − c_g t) exactly; step 11 the kept ω″ term: closed form solves ∂_t env = (iω″/2)∂²_ξ env with the right initial value, \|env\|² peaks at ξ = 0, spreads as 1/√(1 + (ω″s²t)²); = `gaussian_packet(order=2)` 1e-12 |
| D21 | ★★ | `test_group_velocity_V2_derivation` | steps 1–7, sech²/tanh = 2/sinh 2x, (7.69), limits (7.70), and the coded capillary–gravity c_g |
| D22 | ★★ | `test_crest_conservation_V2_derivation` | (7.74) from mixed partials; chain rule; constant along dx/dt = c |
| D23 | ★★ | `test_frequency_along_rays_V2_derivation` | (7.78), (7.79) = c_g × (7.74), and dk/dt = −(∂ω/∂x)_k, with ω = √(gk tanh kH(x)) |
| D24 | ★★ | `test_snell_V2_derivation` | ∂ω/∂y = 0, \|k\| sin α = const, α → 0 as H → 0 |
| D25 | ★★ | `test_hydraulic_jump_V2_derivation_belanger` | face force ½ρgH², (7.80), common factor, Q², 2Fr² = r(1 + r), Vieta, the positive root |
| D26 | ★★ | `test_hydraulic_jump_V2_derivation_energy_loss` | steps 2–7 each |
| D27 | ★★ | `test_stokes_drift_V2_derivation` | ξu_x = a²ωke^{2kz}sin², ζu_z = …cos², sum time-independent, ⟨u⟩ = 0; (7.86) by the same moves at any depth |
| D28 | ★★ | `test_interface_V2_derivation` | Laplace, A = iωa/k, A = −B, φ_t, pressure condition → (7.95), ρ₁ → 0 |
| D29 | ★★★ | `test_two_layer_constants_V2_derivation` | steps 3, 4, 5–6 (7.106)–(7.107), 7, 8 (C = A − Be^{2kH}, (7.108)), 9, 10, 11 (both pieces), 12 (7.109), H → 0 check; printed (7.105) fails; `two_layer_sympy` residuals all 0 |
| D30 | ★★★ | `test_two_layer_dispersion_V2_derivation` | step 3 = step 5 regrouping; step 6 (−iωA, −iωB); step 8 both sides in s; step 9 factor (s − 1)/s; step 11 hyperbolic form; step 12 = (2/s)·(7.110); residual of (7.101) = (ag²k/ω²)·(7.110) (also for the coded `two_layer_sympy`) |
| D31 | ★★ | `test_two_layer_modes_V2_derivation` | (7.112); step 5 brackets; step 6 b = −aρ₁e^{kH}/(ρ₂ − ρ₁); (7.114); kH → ∞ → (7.95); series → (7.115); c² = g′H; (7.119) |
| D32 | ★★ | `test_boussinesq_linear_V2_derivation` | ε-expansion of (7.125): O(ε) = (7.126), O(ε²) the dropped products; (7.131) with N²; `boussinesq_linear_sympy` r7126–r7131 = 0 |
| D33 | ★★★ | `test_w_equation_V2_derivation` | steps 1–2 (swap), 3–4 (7.132) as a residual combination, 5–7 (7.133), 8 (∇_H² with N(z)), 9–12: p′ eliminated exactly for arbitrary N(z); `boussinesq_linear_sympy` r7132–r7137 = 0 |
| D34 | ★★ | `test_internal_dispersion_V2_derivation` | ∇² → −K², (7.137), N\|k\|/K for either sign of k, cos θ = \|k\|/K |
| D36 | ★★ | `test_internal_group_velocity_V2_derivation` | ∂ω/∂k, ∂ω/∂m, (7.145), (7.144), c·c_g = 0, c_z = −c_g,z; sign-safe form for k = ±q |
| D37 | ★★ | `test_polarization_and_flux_V2_derivation` | p̂, ρ̂, û from (7.132), (7.131), (7.128); ½Re(p̂û*), ½Re(p̂ŵ*); c_gE with ω = kN/K = F (7.159); E_k = E_p |
| a-D40 (demoted) | ★★★ | `test_stokes_expansion_V2_third_order_coefficients` | independent third-order expansion (ours) + the coded one |
| a-D46 (demoted) | ★★★ | `test_kdv_V2_solitary_wave_residual` | coded tanh-polynomial residual 0 (3 → wrong coefficient ≠ 0) and our sech form at 60 random points (1e-12) |
| a-D68 (demoted) | — | `test_viscous_decay_V2_derivation` | S_ijS_ij of (7.47), ∫2νS_ijS_ij dz, dE/dt = −D → da/dt = −2νk²a |

No intermediate line of Part F was found wrong: every checked line follows from the previous one. Two statements
*outside* Part F are wrong and are reported to the designer/analyst: analysis §6 V25 ("(ka)² × aω") — the relative
scaling is (ka)¹ (C02 above), and analysis §6 V8 / the `fenton_mckee_kh`, `guo_kh` docstrings ("≤ 1.5 %", "≤ 0.7 %").

## Functions used by the notebook and explainers (design Part C) — test name each
Every C.1 and C.2 callable is asserted to exist and the nine `*_state` explainer functions to return scalar-only dicts
(`test_part_c_V1_every_contract_function_exists_and_is_scalar_callable`); each also has a physics test:
`sinusoid`, `crest_positions` → `test_sinusoid_V1…`; `wave_parameters` → `test_wave_parameters_V1…`; `plane_wave`,
`phase_velocity_vector`, `trace_velocities` → `test_plane_wave_V1…`, `…_V7…`; `doppler_frequency` → `test_doppler_V1…`;
`omega_gravity` → `test_dispersion_V1…`, `…_V7…`; `omega_capillary_gravity`, `phase_speed` → `test_capillary_V7_limits`,
`test_phase_speed_V7…`; `period_from_wavelength`, `wavelength_from_period`, `wavenumber_from_omega` →
`test_wavenumber_from_omega_V1…`, `test_dispersion_V2_derivation`; `fenton_mckee_kh`, `guo_kh` → `test_dispersion_V5_*`;
`group_velocity` → `test_group_velocity_V1…`, `…_V7…`; `group_velocity_numeric` → `…_V1_complex_step_parity`,
`test_group_velocity_numeric_V3…`; `group_velocity_vector` → `test_internal_wave_velocities_V1…`; `depth_regime`,
`wave_regime_label` → `test_depth_regime_V1…`; `capillary_minimum`, `min_group_velocity` → `test_capillary_minimum_V1…`;
`beat_wave` → `test_beat_wave_V1…`; `gaussian_packet` → `test_gaussian_packet_V1…`, `test_packet_envelope_V2…`;
`linear_evolve` → `test_linear_evolve_V1…`, `…_V7…`, `test_energy_flux_V4…`; `envelope` → `test_envelope_V1…`;
`linear_evolve_2d` → `test_linear_evolve_2d_V1…`; `ray_trace` → `test_ray_trace_V1…` (×2), `…_V4…`;
`wave_energy_density`, `viscous_decay` → `test_wave_energy_V1…`, `test_viscous_decay_V2…`; `interface_omega`,
`eps2_density` → `test_interface_omega_V7…`; `two_layer_free_surface_omega`, `two_layer_long_wave_speed`,
`reduced_gravity_book` → `test_two_layer_V7…`; `internal_wave_omega`, `beam_angle` → `test_internal_wave_omega_V7…`;
`internal_wave_velocities` → `test_internal_wave_velocities_V1…`; `real_field` → `test_real_field_V1…`;
`depth_profiles`, `cosh_over_sinh`, `sinh_over_sinh`, `cosh_over_cosh` → `test_pressure_response_V1…`;
`surface_normal`, `surface_velocity` → `test_surface_normal_V1…`; `linear_bernoulli_pressure` →
`test_linear_bernoulli_V1…`; `surface_wave_sympy` → `test_potential_V2…`, `test_wave_fields_V2…`; `wave_fields` →
`test_wave_fields_V1/V3/V7…`; `free_surface_residuals`, `free_surface_residual_scan` → `test_free_surface_V1…` (×2),
`…_V3…`; `pressure_response` → `test_pressure_response_V1…`; `dispersion_state` → `test_dispersion_state_V1…`;
`particle_path` → `test_particle_path_V1…`, `test_stokes_drift_V3…`; `orbit_linear`, `orbit_semi_axes` →
`test_orbit_V1…`; `orbit_state` → `test_orbit_state_V1…`; `dyed_line` → `test_dyed_line_V1…`; `wave_energy` →
`test_wave_energy_V1…`; `energy_flux` → `test_energy_flux_V1…`; `curvature`, `capillary_surface_pressure` →
`test_curvature_V1…`; `capillary_state` → `test_capillary_state_V1…`; `standing_wave_fields` →
`test_standing_wave_V1…`; `seiche_modes`, `basin_modes` → `test_seiche_V1…`, `…_V3…`; `seiche_state` →
`test_seiche_state_V7…`; `packet_state` → `test_packet_state_V1…`; `pond_ripples` → `test_pond_ripples_V3…`;
`local_wavenumber_frequency`, `crest_conservation_residual` → `test_local_wavenumber_V3…`; `snell_ray_plane_beach` →
`test_ray_trace_V1_snell_closed_form_parity`; `refraction_state` → `test_refraction_state_V1…`;
`nonlinear_wavelet_speed`, `simple_wave_evolve` → `test_simple_wave_V1…`; `hydraulic_jump`, `jump_momentum_residual` →
`test_hydraulic_jump_V1/V4…`; `jump_state` → `test_jump_state_V1…`; `stokes_wave_profile`, `stokes_wave_speed`,
`STOKES_LIMIT_STEEPNESS` → `test_stokes_wave_profile_V1…`, `test_stokes_wave_V5…`; `stokes_expansion_sympy` →
`test_stokes_expansion_V2…`; `stokes_drift`, `stokes_drift_numeric`, `eulerian_mean_u` → `test_stokes_drift_V1/V3/V5…`;
`kdv_rhs`, `kdv_solve`, `kdv_invariants`, `kdv_linear_phase_speed`, `ursell_number` → `test_kdv_solve_V1/V3/V4…`,
`test_kdv_linear_phase_speed_V1…`; `solitary_wave`, `solitary_wave_speed`, `cnoidal_wave`, `kdv_residual_sympy` →
`test_kdv_V2…`, `test_cnoidal_V1…`; `interface_fields`, `interface_residuals` → `test_interface_fields_V1…`;
`interface_energy` → `test_interface_energy_V1…`; `two_layer_modes`, `two_layer_residuals` →
`test_two_layer_modes_V1…`, `test_two_layer_residuals_V1…`; `two_layer_sympy` → `test_two_layer_constants_V2…`,
`…_dispersion_V2…`; `two_layer_rigid_lid_omega` → `test_two_layer_V7…`; `two_layer_state` → `test_two_layer_state_V1…`;
`boussinesq_linear_sympy` → `test_boussinesq_linear_V2…`, `test_w_equation_V2…`; `internal_wave_fields` →
`test_internal_wave_fields_V1/V3…`; `w_equation_residual` → `test_w_equation_residual_V1…`; `layered_flow_check` →
`test_layered_flow_V1…`; `internal_wave_energy`, `internal_energy_budget_residual` → `test_internal_wave_energy_V1…`,
`test_internal_energy_budget_V4…`; `internal_pe_interface_limit` → `test_internal_pe_interface_limit_V3…`;
`st_andrews_cross` → `test_st_andrews_cross_V7…`; `internal_wave_state` → `test_internal_wave_velocities_V1…`; C.0
re-used functions are exercised through parities (ch04 `linear_wave_surface`, `bore_speed`, `kinematic_bc_residual`,
`unsteady_bernoulli_pressure`, `velocity_from_streamfunction_2d`, `laplace_jump_from_balance`, `reduced_gravity`; ch03/core
`pathline`; ch05 `vortex_sheet_strength`; ch01 `surface_tension_water`, `water_density`; core `parcel_displacement`); C.3
scripts → `test_scripts_V1_every_ch07_script_runs`.

## Convergence studies (scheme | steps | observed order (pairwise) | design order)
| scheme | steps | observed order | design |
|---|---|---|---|
| exact free-surface residuals vs ka (C02), absolute kin./dyn. kH = 1 | ka = 0.005…0.08 | 2.011 / 2.002 | 2 |
| same, relative (÷ aω, ÷ ag) | same | 1.011 / 1.002 (deep 1.014 / 1.022) | 1 |
| 5-point ∇²φ of `wave_fields` | h = 0.2…0.025 | 2.0 | 2 |
| packet centroid speed → c_g (`linear_evolve`) | σ_x = 10, 20, 40 (δk halved) | 2.0 (1.9e-3, 4.7e-4, 1.2e-4) | 2 |
| symmetric chord Δω/Δk → dω/dk | Δk = 0.2…0.025 | 2.0 | 2 |
| `group_velocity_numeric` central fallback | h = 0.08…0.01 | 4.03, 4.01, 4.00 | 4 |
| `pond_ripples` midpoint rule → Cauchy–Poisson integral | n = 256, 512, 1024 | 2.0 | 2 |
| FD sloshing eigenproblem → (7.65) (ours) | nx = 20, 40, 80, 160 | 1.86, 1.97, 1.99 (fit on 40–160: 1.98) | 2 |
| seiche shallow-limit error | H/L = 0.02, 0.01, 0.005 | 2.0 | 2 ((kH)²/6) |
| local wavenumber (4th-order FD) | h = 0.4…0.05 | 4.00, 4.00, 4.00 | 4 |
| KdV IF-RK4 in time | dt = 0.04, 0.02, 0.01 | 3.97, 4.00 | 4 |
| KdV in space | N = 64 → 128 | error ÷ 470 | spectral |
| linear KdV speed vs (7.29) | kH = 0.05…0.4 | 4.0 | 4 |
| Stokes drift, exact paths → (7.86) | ka = 0.01…0.08 | 0.976 (0.99, 0.98, 0.95) | 1 |
| published Stokes speed vs (7.83) | ka = 0.02…0.16 | 4.0 | 4 |
| wavelet speeds vs c₀(1 + 3η/2H) | η/H = 5e-5…5e-3 | 2.0 | 2 |
| internal-wave FD residuals (4.10), (7.128), (7.130), (7.131) | h = 1e-2…2.5e-3 | 2.00 each | 2 |
| w-equation FD residual | h = 0.08…0.02 | 2.0 | 2 |
| internal energy budget (7.147) FD residual | h = 1e-2…2.5e-3 | 2.00 | 2 |
| smoothed density jump ε → 0 (N157) | ε = 4…0.25 m | 0.97, 0.98, 0.99, 1.00 | ≥ 1 |

## Conservation / invariant residuals
| invariant | residual |
|---|---|
| ∫η² of a one-way linear packet (`linear_evolve`) | ptp/E < 1e-12 |
| KdV mass (splitting hump, 20 s) | 2e-16 relative |
| KdV momentum ∫η², Hamiltonian | 7.0e-10, 9.1e-10 relative (default dt); ÷ 29 when dt halves (truncation, not a leak) |
| ω along a ray over a 1:50 beach (T = 8 s, to H = 2 cm) | 1.9e-9 relative (design D23 check: < 1e-8) |
| k_y along the same ray (Snell) | 2.3e-12 relative |
| hydraulic jump momentum (7.80) | ≤ 1e-13 m³/s²; mass exact to 1e-14 |
| second law: E₂ − E₁ ≤ 0 for Fr₁ ≥ 1 | 91 values of Fr₁ ∈ [1, 10], all ≤ 0, monotone |
| internal-wave energy equation (7.147), pointwise | 5.8e-5 W/m³ at h = 2.5e-3 (terms ≈ 6.9 W/m³), order 2; ⟨gρ′w⟩ = 0 to 1e-12 |
| c·c_g (internal waves) | ≤ 1e-15 |

## Benchmarks used (value | our value | source + URL | date verified)
| benchmark | published | ours | source | verified |
|---|---|---|---|---|
| Fenton & McKee (1990) Eq. (21), max wavelength error | "always better than 1.7 %" | 1.658 % (λ), 1.631 % (kd) | Coastal Eng. 14, 499–513, doi:10.1016/0378-3839(90)90032-R, https://johndfenton.com/Papers/Fenton90c+McKee-On-calculating-the-lengths-of-water-waves.pdf | 2026-09-24 (PDF read) |
| same, optimal exponent | ν = 1.49 | 1.483 | same, p. 507 | 2026-09-24 |
| Guo (2002), β = 2.4908 | 0.75 % | 0.753 % | Coastal Eng. 45, 71–74, doi:10.1016/S0378-3839(02)00039-X (**search-confirmed**: publisher pages 403) | 2026-09-24 |
| Fenton (2006) note | FM "within 1.5 %"; Guo (5/2) "about 0.7 %", "about half" | 1.66 % / 0.79 %; ratio 0.48 | https://johndfenton.com/Papers/Dispersion-Relation.pdf | 2026-09-24 |
| capillary minimum, air–water | 0.23 m/s at 1.7 cm | 0.2312 m/s at 1.712 cm (IAPWS σ, ρ = 998.2) | Wikipedia "Capillary wave" | 2026-09-24 |
| limiting Stokes wave | H/λ = 0.1410633 ± 4e-7, 120° | constant reproduced; a_max ≈ 0.0705λ | HandWiki "Physics:Stokes wave" (attributes Schwartz & Fenton); arXiv:1507.02784 (angle) | 2026-09-24 |
| deep-water Stokes drift | ωka²e^{2kz}; "about 4 %" at z = −λ/4 | identical form; 4.32 % | Wikipedia "Stokes drift" | 2026-09-24 |
| Stokes amplitude dispersion | c = (1 + ½(ka)²)√(g/k) | (7.83) equal to O((ka)⁴) | HandWiki "Physics:Stokes wave" | 2026-09-24 |
| Bélanger relations (V1 form) | y₂/y₁, ΔE = (y₂ − y₁)³/4y₁y₂ | identical 1e-12 | Wikipedia "Hydraulic jumps in rectangular channels" | 2026-09-24 |
| KdV, cnoidal, solitary, Ursell (V1 form) | published formulas | identical 1e-13 | Wikipedia "Cnoidal wave" | 2026-09-24 |
| two-fluid capillary–gravity dispersion (V1 form) | ω² = \|k\|[((ρ − ρ′)/(ρ + ρ′))g + σk²/(ρ + ρ′)] | identical 1e-14 (ρ′ = 0 and σ = 0 cases) | Wikipedia "Capillary wave" | 2026-09-24 |
| St Andrew's cross | beams at arccos(ω/N) | geometry (pedigree only) | Mowbray & Rarity, J. Fluid Mech. 28, 1–16 (1967) | 2026-09-24 (search) |

**Implementer's note 3, settled honestly.** The analysis and the `fenton_mckee_kh`/`guo_kh` docstrings quote the 2006
note's "1.5 %" and "0.7 %". Measured against our exact root, the maxima are 1.66 % (in λ) and 0.79 %. The *primary*
sources settle it: Fenton & McKee (1990) state 1.7 % for their formula (reproduced: 1.658 % ≤ 1.7 %, rounds to 1.7) and
ν_opt = 1.49 (ours 1.483); Guo (2002) states 0.75 % for *his* exponent β = 2.4908 (reproduced: 0.753 %), and the
0.79 % belongs to the 5/2 rounding that Fenton's note (and `guo_kh`) uses. The tests assert the primary numbers at their
printed precision and the note's "about half" ratio; the note's 1.5 % is recorded as not reproducible (it matches the
*period* error, 1.49 %, not the wavelength error).

## Numbers from the text (book vs ours) [private values redacted to relative errors]
Source: `tests/book_values_ch07.json` (git-ignored), tested by the three `test_book_V6_*` tests. Every printed closed
form of §§7.2–7.8 in the JSON that is an expression (49 forms: fields, dispersion relations, energies, fluxes, orbits,
drift, jump, KdV, two-layer roots and ratios, internal-wave energies, the Exercise 7.20 rigid-lid form) equals the fluidpy
value at non-trivial points of our choosing — 43 of them to ≤ 1e-12 relative, six (limits and simple products) at
pytest's default 1e-6; the shallow-water orbit forms (7.50) to 1e-7 at kH = 4.5e-5 (their own O((kH)²) truncation). The
vector and prose entries (c, c_g tuples, polarization text) are covered by the V2 derivations instead.
| number | ours vs book | comment |
|---|---|---|
| tanh 2 | 0.00 % | |
| "deep within 2 %" at kH = 2 | ours 1.82 % | consistent |
| H/λ at kH = 2 | −0.53 % | book rounds 0.318 |
| "shallow within 3 %" at H = 0.07λ | ours 3.04 % | rounding |
| λ/H ratio for the shallow threshold | +2.0 % | rounding (1/0.07) |
| wind-wave λ for the book's period | +4.1 % | rounding in the book (analysis T9) |
| deep pressure fraction at λ/2 | +8.0 % | one significant figure in the book (e^{−π} = 4.3 %) |
| c_min, λ_m (7.59) | +0.15 %, +0.23 % | |
| c_g,min (Ex. 7.10) | +0.17 % | |
| a_max/λ of the limiting wave | +0.76 % | from the published H/λ; the book rounds |
| density decrease per 10 °C | 0.26 % vs one printed digit | Kell (1975) densities; rounds to the book's value |
| Fig. 7.33 beam angle | −0.52 % | 44.77° vs the rounded figure label |
| Exercise 7.6 lake period | **+0.72 %** (+0.74 % with g = 9.80665) | analysis T8: no reasonable g or rounding reproduces the book's answer; (7.65) is verified independently (V1 + the FD eigenproblem, V3) — a book-value discrepancy, not a code error (implementer's note 7) |
| Exercise 7.16 speed, g′ | +0.07 %, +0.51 % | g′ from Kell densities vs the book's table |
| Exercise 7.17 angle | exact (60°) | |

## Figures reproduced with our code (outputs/ch07/verify/, local)
| figure | PNG | visual verdict |
|---|---|---|
| Figs. 7.3–7.4 orbits | `fig7_4_orbits.png` | deep circles shrinking like e^{kz₀}, intermediate ellipses flattening with depth, shallow ellipses of nearly equal width collapsing to a line on the bed; the arrows at the top of every orbit point forward (+x): clockwise. |
| Fig. 7.10 c(λ) | `fig7_10_phase_speed.png` | one capillary branch ∝ λ^{−1/2} and one gravity branch ∝ λ^{1/2} meeting at the red minimum 0.231 m/s at 1.71 cm; finite depths saturate at √(gH) (0.99, 3.13, 9.90 m/s); all depths coincide for λ ≪ H. |
| C02 residual scan | `c02_residual_scan.png` | straight lines on log–log axes: absolute exact-condition residual slope 2.01, relative slope 1.01, kH = 1 and deep nearly parallel. |
| Figs. 7.13–7.15 packet | `fig7_15_packet_cg.png` | the Hilbert envelope peaks exactly on the purple c_g t lines (156.6 m, 313 m) while the c t line (orange) runs ahead; the carrier shape is preserved (deep water, 200 s). |
| Fig. 7.20 jump | `fig7_20_jump.png` | H₂/H₁ = 1 at Fr₁ = 1, nearly linear growth beyond; E₂ − E₁ zero at Fr₁ = 1, negative (cubic) above, positive (forbidden) below. |
| Figs. 7.21–7.22 Stokes | `fig7_22_stokes.png` | peaked crests and flat troughs growing from order 1 to 3; the dyed line leans forward more at the top every period; drift error line slope 0.976. |
| Fig. 7.23 KdV | `fig7_23_kdv.png` | a Gaussian hump splits into a rank-ordered train, tallest (fastest) soliton in front plus a small dispersive tail; no oscillation on the right-hand face. |
| Figs. 7.27–7.28 two layers | `fig7_27_two_layer.png` | barotropic ω/k = √(g/k) (slope −½); baroclinic flat at √(g′H) for kH ≪ 1 and merging with the deep interface relation (7.95) for kH ≳ 3. |
| Figs. 7.29, 7.33 internal | `fig7_29_7_33_internal.png` | four beams exactly on the dashed lines at 53.13° from the vertical for ω/N = 0.6, phase lines along the beams; for k < 0, K and c point up-left and c_g down-left at a right angle. |
| V5 Fenton/Guo | `v5_fenton_guo_errors.png` | both errors vanish at both ends; FM peaks at +1.66 % (inside the ±1.7 % band) near ω√(d/g) ≈ 0.6, Guo(5/2) at +0.79 %, both change sign near 1. |
| C08 seiche FD | `c08_seiche_fd_convergence.png` | a straight log–log line of slope ≈ 2 from 1.4e-2 to 2.4e-4. |
No book page, crop or scan was used; every PNG is our code with our parameters.

## Deviations & justifications
No `# DEVIATION` markers in `fluidpy/core/waves.py` or `fluidpy/ch07_gravity_waves.py`. Documented choices that differ
from the printed book and are verified: (7.66) with ½Δω t (printed x fails, D19); (7.105) with e^{i(kx − ωt)} (printed kz
fails, D29); sign-safe (7.138)/(7.145) (printed fails for k < 0); g′ with ρ₂ below (7.117) alongside ch04's ρ₁ form
(ratio ρ₁/ρ₂ asserted); interfacial E_p = ¼Δρga² (the printed middle form's λ/2 upper limit gives ⅛ — shown numerically);
KdV integrated with our IF-RK4 pseudo-spectral scheme (book gives none; order 4 and invariants verified); the Riemann
simple wave and the St Andrew's cross are labelled extensions/illustrations; ω(k) inverses by `brentq` with physics
brackets (residual asserted inside the function).

## Open items
O1 (implementer, docstrings) — 81 docstrings still say "Validation (planned)" (35 in `fluidpy/core/waves.py`, 46 in
`fluidpy/ch07_gravity_waves.py`). Replace each with the test names of this report (as ch06 did in its loop 2). Not a
physics defect; it keeps the evidence traceable.

O2 (implementer, docstring numbers/citations) — (a) `fenton_mckee_kh` "≤ 1.5 %" → the primary Fenton & McKee (1990)
bound is 1.7 % (ours 1.66 % in λ); `guo_kh` "≤ 0.7 %" → 0.79 % for the coded exponent 5/2 (Guo's own β = 2.4908 gives
0.75 %); the same numbers in analysis §6 V8 and §8. (b) `STOKES_LIMIT_STEEPNESS`: HandWiki attributes 0.1410633 to
Schwartz & Fenton, not to Dyachenko, Lushnikov & Korotkevich (whose abstract gives the 120° angle only) — value correct,
citation text wrong (also analysis §8). (c) `ray_trace` claims "ω drifts < 1e-9 up to the shore": measured 1.9e-9 on a
1:50 beach to H = 2 cm (the design's D23 check, < 1e-8, holds). (d) `internal_pe_interface_limit`: "≈ 1.2 kε, defaults
give 1.2 %" → measured 1.37 % (≈ 1.37 kε). (e) `depth_regime`: "3.05 %" → 3.04 %. (f) module docstring of `core.waves`
says point sets carry components on the *first* axis while `plane_wave` documents the *last* (the code accepts both).

O3 (analyst/designer) — analysis §6 V25 says the exact-condition residuals of the linear solution scale as (ka)² × aω.
They scale as (ka)·aω (relative slope 1.01, absolute slope 2.01 at fixed k); design Part F D03 already states it
correctly — the notebook's C02 text must follow D03, not V25 (implementer's note 2 confirmed).

O4 (derivation reviewer) — Exercise 7.2 set up literally (potential amplitude fixed at aω/k, kinematic condition to (ka)¹)
gives γ = 3/8 in c² = (g/k)(1 + γk²a²); the consistent third-order Stokes solution (the coded one, and our independent
expansion) gives γ = 1, i.e. (7.83), with the potential's first-harmonic amplitude corrected by −(ka)²/8. The notebook
should present (7.82)–(7.83) with the consistent derivation and mention why the truncated set-up misses γ (implementer's
note 1 confirmed).

O5 (book value, user) — Exercise 7.6 lake period: our (7.65) is 0.72 % above the printed answer (0.74 % with g = 9.80665);
(7.65) is verified independently (V1 and the finite-difference sloshing eigenproblem, V3). Reported as a book-value
discrepancy (analysis T8), not a code error (implementer's note 7 confirmed).

O6 (benchmark provenance) — Guo (2002)'s 0.75 % and β = 2.4908 are *search-confirmed* (publisher pages returned 403);
reproduced to the printed digits. If the user has journal access, a direct read of the abstract would upgrade it.

O7 (benchmark, open since the analysis) — the minimum group velocity (Exercise 7.10) has no verified public primary
number; its evidence is V1 (`minimize_scalar`) + V2 (sympy stationarity 3x² + 6x − 1 = 0) + V6 (0.17 %).

O8 (ch04, low) — `ch04.linear_wave_surface` uses plain cosh/sinh and returns NaN for kH ≳ 710 (found through the parity
test at kH = 900). ch07's overflow-safe `wave_fields` supersedes it; suggest re-pointing ch04 to `core.waves.cosh_over_sinh`.

O9 (implementer, low) — the complex-step path of `core.waves._tanh_kH` emits "RuntimeWarning: overflow encountered in
tanh" for kH ≳ 355 (numpy's complex tanh); the returned c_g is still correct (1e-12). Wrap in `np.errstate(over="ignore")`.
**Closed in loop 2:** `_tanh_kH` now uses `np.errstate(over="ignore")`; asserted by the kH = 800 case of
`test_group_velocity_V1_negative_k_is_signed_derivative` under `warnings.simplefilter("error")`.

O10 (qualitative, by design) — `st_andrews_cross` is a labelled illustration (Mowbray & Rarity 1967 pedigree): its beam
geometry (arccos(ω/N), c ⟂ c_g, phase lines along beams, energy away from the source) is verified (V7); it is not a
solution of the forced problem. The Riemann simple wave (`simple_wave_evolve`, N81) is our labelled extension, verified
analytically (breaking time, first-order speed), without a published benchmark.

## Loop 2 (post-review) — 2026-09-24, commit 04e4606 + working tree

**The change** (review finding M1, negative k in the complex-step group velocity). `core.waves._abs_k` is now the
analytic continuation of |k| for complex input (−k where Re k < 0, k elsewhere); before, it returned a complex k
unchanged, so a complex step at k < 0 evaluated √(g k) on the wrong branch and gave c_g ≈ 1e20 m/s. `group_velocity`
now returns the signed derivative dω/dk = sgn(k)·(c/2)[… + 2|k|H/sinh 2|k|H] of ω(|k|) from (7.28)/(7.56) (identical to
(7.69) for k > 0); `phase_speed` still returns the speed |ω/k| ≥ 0. `_tanh_kH` now silences numpy's spurious complex-tanh
overflow flag (closes O9).

**Existing expectations unchanged.** `git diff HEAD -- tests/test_ch07.py` is 54 added lines (the two tests below) and
no modified or deleted line; every pre-existing ch07 test passes unchanged against the new sign convention.

**The two new tests (both V1), audited independently.**

| Test | Physics basis of the expectation | Tolerances | Honest? |
|---|---|---|---|
| `test_group_velocity_V1_negative_k_is_signed_derivative` | ω = √(\|k\|(g + σk²/ρ) tanh \|k\|H) (7.28)/(7.56) depends on \|k\| only ⇒ dω/dk is odd in k; asserts oddness exactly (`array_equal(neg, −pos)`), sign(c_g) = sign(k) per element, complex step = analytic, an **independent real ω(\|k\|) lambda** differenced by 4th-order central differences = analytic, and the closed forms at k = −1: −½√g (7.70 deep) and −(c/2)(1 + 4/sinh 4) (7.69, H = 2); shallow limit −√(gH) (7.70). 6 k × 5 (H, σ) cases incl. kH = 800 under `warnings.simplefilter("error")` | complex step 1e-12 (`maxrel`, normalised by the array max; **element-wise worst is 3.2e-16**, measured); central 1e-9 (truncation h⁴ with h = 1e-3\|k\| plus round-off ε/h ≈ 1e-13); closed forms 1e-13; shallow limit 1e-9 at kH = 1e-5 (the c_g correction ½(kH)² is 5e-11) | yes |
| `test_gaussian_packet_V1_negative_k0_is_the_mirror_image` | η₀ = a e^{−x²/2σ²} cos k₀x is even in x and k₀ and ω(\|k\|) is even ⇒ η(x, t; −k₀) = η(−x, t; k₀) (deep and H = 5, orders 1 and 2); order-1 envelope rides rigidly left at c_g = −½√(g/\|k₀\|) (7.68)/(7.70); order 2 vs the FFT evolution `linear_evolve` for a quadratic Ω(\|k\|) at k₀ = −1 (c_g = −0.8, ω″ = 0.6 exactly) | c_g 1e-13/1e-12 (closed-form arithmetic); ω″ 1e-8 (finite difference of complex-step values, step 1e-3\|k₀\|: error O(h²ω⁗) ≈ 1e-9); η 1e-12 of a = 1; envelope 1e-14 | yes |

**Planted variants** (monkeypatched in memory from a scratch script, `fluidpy/` untouched):

| Variant | signed-derivative test | mirror-packet test |
|---|---|---|
| A — old `_abs_k` (complex k returned unchanged, HEAD) | **fails** (sign / complex-step parity) | **fails** (c_g ≈ 1e20, c_g(−k₀) ≠ −c_g(k₀)) |
| B — unsigned analytic `group_velocity` (HEAD, \|c_g\|) | **fails** (oddness) | **fails** (c_g ≠ `group_velocity(−0.5)`) |
| C — full HEAD behaviour (A + B) | **fails** | **fails** |
| E — continuation that loses the sign (−Re k + i Im k) | **fails** (sign of c_g) | **fails** |
| F — `gaussian_packet` moves at \|c_g\| | passes (out of scope) | **fails** |
| D — `np.abs` in `_abs_k` (non-analytic; ω real, `auto` falls back to central differences) | passes — correctly: the fallback returns the right signed c_g to 4e-13 element-wise, not a physics error | **fails** (η 1e-12 bound: ω″ loses precision) |

Every behaviour that is physically wrong (A, B, C, E, F) is caught by at least one of the two tests; the old code fails
both. 6/6 variants caught by the pair (41/41 cumulative for ch07).

**Runs.** `pytest tests/test_ch07.py -q -p no:cacheprovider` → **132 passed in 72 s**; full suite `pytest tests` →
**882 passed in 409 s**, no warnings reported. Collected items per tag: **V1 60 · V2 36 · V3 10 · V4 5 · V5 5 · V6 3 ·
V7 13** (= 132). C09 (group velocity) now has V1 evidence for both signs of k.

## Verdict: PASS (loop 2) — 132/132 ch07 tests pass (72 s; full suite 882 passed); M1 fixed and proven by two discriminating V1 tests (old behaviour fails both); no existing test expectation changed. Loop 1 verdict: 130/130 ch07 tests pass (125 s; full suite 880 passed); all 16 CORE items have V1 and V2 evidence and 15 of them a third level (V3, V4, V5 or V7; C05 has V1 + V2); every coded NOTE ≥ 1; every ★★/★★★ derivation re-derived with sympy (no wrong intermediate line); every design Part C function and all 11 scripts exercised; 35/35 planted wrong variants caught; no tolerance loosened; open items are docstring/citation fixes, book-value and provenance notes, and two labelled qualitative illustrations.
