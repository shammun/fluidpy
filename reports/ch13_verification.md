# Chapter 13 verification — Geophysical Fluid Dynamics                       2026-10-07 · code on the working tree at `68f22a4` (ch13 files staged, uncommitted)

## Verdict: PASS (re-verified after the code review) — `tests/test_ch13.py`: 91 passed, 0 failed, 0 skipped (121 s with the two slow tests); 39 of 39 planted mutants killed; full suite `-m "not slow"`: 1865 passed, 0 failed, 19 deselected (2994 s on a machine shared with other runs)

Every script ran (13 chapter scripts and the cache writer exit 0), every test passes, each of the 17 CORE rows has at
least two independent evidence levels (one of V1 / V2 / V3 / V5 in each), every coded NOTE item and every one of the 182
names of design Part C is called by at least one test (a self-audit test enforces it), all 27 two- and three-star
derivations are re-derived with sympy independently of the chapter's own engines, the six planted printed variants fail
their independent checks while the corrected forms pass, and nothing is `unverified` without an Open item. **No physics
defect was found in `fluidpy/`.** Four documentation-level findings and two design/curation wording errors are listed
under "Open items" and "Found wrong in design / analysis"; none blocks.

## Environment
python 3.11.5 · numpy 2.4.6 · scipy 1.17.1 · sympy 1.14.0 · pint · matplotlib (Agg) · Windows 11.
**After the review (this revision):** `tests/test_ch13.py` → **91 passed in 120.6 s** (all, slow included; 89 non-slow); full suite `-m "not slow"`: 1865 passed, 0 failed, 19 deselected (2994 s on a machine shared with other runs). First pass, before the review:
`.venv/Scripts/python.exe -m pytest tests/test_ch13.py -q -p no:cacheprovider` → **89 passed in 87.55 s** (all, including the two
`slow` tests: 13 scripts 35 s, cache regeneration 13 s); with `-m "not slow"` → **87 passed, 2 deselected in 47 s**.
`tools/check_public.py tests/test_ch13.py tests/ch13_verify_figures.py reference/ch13/SOURCES.md reference/ch13/make_refs.py
reference/ch13/benchmarks.json` → "public-repo check OK (5 files)". Rest of the suite and `tests/test_machinery.py`: see
the last section.

## The ten claims the orchestrator asked about — rulings

| # | Claim | Ruling | Measured |
|---|---|---|---|
| 1 | `shallow_water_omega` solves the cubic $\omega^3-c^2\omega K^2-f_0^2\omega-c^2\beta k=0$ (13.76) | **Confirmed.** 648 points (6 latitudes of both signs × 3 speeds × 9 wavelengths × 4 directions) | relative residual ≤ 1.9e-16; sum of roots/largest ≤ 1e-12; product = $c^2\beta k$ to 1e-9; fast roots equal `numpy.roots` to 1e-9 (the slow root of `numpy.roots` carries the round-off, ours does not); NaN **exactly** where `shallow_water_discriminant` < 0 (400-point sweep; there `numpy.roots` returns a complex pair); limits → Poincaré $\omega^2=f^2+gHK^2$ (13.82) and Rossby $\omega=-\beta k/(k^2+l^2+f_0^2/c^2)$ (13.118) to 2e-3, → $cK$ to 3e-4 at $K\Lambda>50$, → ch07 long waves at $f=\beta=0$ |
| 2 | $\int\psi_m\psi_n\,dz=0$ with weight 1 for the rigid lid **and** the free surface | **Proved.** Lagrange identity for arbitrary $N^2(z)$ from $\frac{d}{dz}\big(\frac1{N^2}\frac{d\psi_n}{dz}\big)+\frac1{c_n^2}\psi_n=0$ (13.56): the boundary bracket $(\psi_n\psi_m'-\psi_m\psi_n')/N^2$ vanishes at the bottom by $d\psi_n/dz=0$ (13.64) and at the top by $\frac{d\psi_n}{dz}+\frac{N^2}{g}\psi_n=0$ (13.65), because that condition does not contain the eigenvalue | exact uniform-N free-surface modes by adaptive quadrature: overlap ≤ 1e-12; with a surface term added it is 3e-5 (refuted); finite-volume modes (uniform N and thermocline, both lids): off-diagonal ≤ 1e-12; the surface term belongs to the energy relation only (`kind="energy"`, ≤ 1e-8; dropping it there gives > 1e-6) |
| 3 | Bottom Ekman layer: largest u is 1.0670 U at z = 3πδ/4 | **Confirmed**, symbolically ($du/ds=0$ at $s=3\pi/4$, $u/U=1+e^{-3\pi/4}/\sqrt2$) and numerically (200 001 points, both hemispheres) | 1.067020 U at z/δ = 2.35619; at z = πδ: $1+e^{-\pi}$ = 1.043214 U |
| 4 | Eady numbers | **Confirmed** by our own root-finding of $c=\frac{U_0}{2}\pm\frac{U_0}{\alpha H}\sqrt{\big(\frac{\alpha H}{2}-\tanh\frac{\alpha H}{2}\big)\big(\frac{\alpha H}{2}-\coth\frac{\alpha H}{2}\big)}$ (13.141) | cut-off 2.39936 (agrees to 1e-12), fastest wave 1.60612, growth coefficient 0.30982; `eady_numeric_eigs` (Chebyshev) equals the closed form to 2e-8 U₀ at seven wavenumbers on both sides of the cut-off, both hemispheres; spectral convergence in n; published values 0.3098 and 1.606 matched to every printed digit (V5) |
| 5 | Geostrophic adjustment (ours) | **Confirmed.** ODE solved from scratch with `dsolve`; linear PV kept point by point; KE/PE released = 1/3 exactly (sympy) and to 1e-14 (code vs quadrature); jet sign = sign(f) | the 1-D model: discrete PV drift < 1e-12 (round-off), time-mean over the 6th inertial period within 0.8 % of the end state (errors fall 1.3 % → 0.6 % → 0.5 % over periods 2, 4, 6), numerical KE/PE = 0.335 |
| 6 | Westward flow over a step ends displaced by $Y_p=f_0(h_1-h_0)/(\beta h_0)$ | **Slip in the book (false as printed), confidence high (about 90 %).** Pages re-read as images (printed pp. 661–662): the figure of the westward case draws the streamline back at its upstream latitude over the changed depth, its caption says the deflection "recovers", and the text says the particle ends "at its original latitude". With $\frac{D}{Dt}\big(\frac{\zeta+f}{h}\big)=0$ (13.94) and a uniform stream far downstream (ζ = 0) the only possible latitude is $f=f_0h_1/h_0$, i.e. $Y=Y_p$ — an exact statement, no linearisation. "Back at the original latitude" misses PV conservation by $(h_0-h_1)/h_1$. The remaining 10 %: the authors may have had a ridge of finite width in mind (then the stream does return), but both the caption and the sketch say "step change in depth". The page's other statements (upstream influence, counter-clockwise curvature before the step, sudden drop of vorticity at the step, no oscillation) are all reproduced by our solution | coded streamlines conserve PV to 2e-7 on both sides for U > 0 and U < 0; eastward: meander between 0 and 2Y_p about Y_p (as the book's own eastward figure shows); a mutant that returns to Y = 0 is killed |
| 7 | Numerical models | **Confirmed.** | orders and residuals in the two tables below; the 2/3 rule is strict (a product landing exactly on 2/3 of Nyquist is removed; the `<=` mutant is killed); `wkb_error` × (H m) → 0.230 (constant to 0.01 %) |
| 8 | Cached runs reproduce | **Confirmed.** `scripts/ch13_make_caches.py --dest <tmp>` regenerates all four `.npz` files within 4.1e-4 of each array's range (stated tolerance 2e-3: float16 storage) and `explainer_constants.json` identically | 13 s |
| 9 | Rayleigh–Kuo with β (sech² jet, k = 1.4) | **Confirmed for β ≥ 0**, by an independent shooting integration written in the test file (radiation condition in the far field): 0.119126, 0.074385, 0.015074 at β = 0, 0.3, 0.6 (agree to 1e-7 in c); no growing mode at β = 0.7 > max U″ = 2/3. **k = 0.9, β = −1 characterised:** not a resolution problem but a boundary-condition one — see Open item 1 | integral identity $c_i\int(\beta-U'')\lvert\phi\rvert^2/\lvert U-c\rvert^2dy=0$ holds to 2e-6 on the shooting mode |
| 10 | f → 0 limits; (13.112) | **Confirmed.** `inertia_gravity_omega(…, f=0)` equals ch07's `internal_wave_omega` to 1e-14 (with and without l), the group velocity equals ch07's `group_velocity_vector` to 1e-6, the 1-D model with f = 0 moves a pulse at $\sqrt{gH}$. Page re-read: (13.112) is printed only as $\omega^2-f^2=\frac{k^2}{m^2}(N^2-\omega^2)$; the code's explicit form satisfies it to 1e-10 on a 40-point grid | — |

## Validation table

Evidence labels: `analytic` V1 · `symbolic` V2 · `converged` V3 · `conserved` V4 · `benchmark` V5 · `book-value` V6 (private) · limits/symmetry V7.

| Concept / Eq. | Tier | fluidpy target | Evidence | Numbers | Label | Tests |
|---|---|---|---|---|---|---|
| C01 thin-shell equations, $2\boldsymbol\Omega\times\mathbf u\cong(-fv,\ fu,\ -2\Omega u\cos\theta)$ (13.7), $f=2\Omega\sin\theta$ (13.8), $f=f_0+\beta y$ (13.10), friction (13.5)–(13.6) | CORE | `coriolis_acceleration_local`, `beta_plane*`, `thin_layer_terms`, `thin_layer_residuals`, `eddy_friction`, `anisotropic_eddy_stress` | V1, V2 (D01, D02), V2 dimensions, V7 | full cross product to 1e-13; β-plane error order 2.00; residuals of a balanced state < 1e-15 | analytic, symbolic | `test_planetary_V1_*`, `test_beta_plane_V1_*`, `test_coriolis_components_V2_derivation`, `test_friction_force_V2_derivation`, `test_friction_force_V2_dimensions_printed_form_fails`, `test_boussinesq_V2_rest_state_printed_sign_fails`, `test_thin_layer_V1_term_sizes_and_residuals` |
| C02 geostrophic balance $-fv=-\frac1{\rho_0}\frac{\partial p}{\partial x}$, $fu=-\frac1{\rho_0}\frac{\partial p}{\partial y}$ (13.11)–(13.12), Ro (13.13) | CORE | `geostrophic_velocity`, `geostrophic_from_field`, `geostrophic_from_height`, `geostrophic_streamfunction`, `pressure_centre`, `parcel_adjust`, `rossby_number`, `ekman_number` | V1, V3, V2, V7 (both hemispheres) | u·∇p = 0 to 1e-12; gridded gradient order 2.0 (edges included); parcel ODE to 1e-6 | analytic, converged | `test_geostrophic_V1_flow_along_isobars_both_hemispheres`, `test_geostrophic_V3_gridded_gradient_is_second_order`, `test_parcel_adjust_V2_ode_and_geostrophic_limit` |
| C03 thermal wind $\frac{\partial v}{\partial z}=-\frac{g}{\rho_0f}\frac{\partial\rho}{\partial x}$, $\frac{\partial u}{\partial z}=\frac{g}{\rho_0f}\frac{\partial\rho}{\partial y}$ (13.15); Taylor–Proudman $\partial\mathbf u/\partial z=0$ (13.21) | CORE | `thermal_wind_shear`, `thermal_wind_from_temperature`, `thermal_wind_integrate`, `jet_section`, `taylor_proudman_residual` | V2 (D04, D05), V1, V4, V7 | jet in balance to 1e-10; westerly shear in both hemispheres | symbolic, analytic | `test_thermal_wind_V2_derivation`, `test_thermal_wind_V1_jet_section_in_exact_balance`, `test_taylor_proudman_V2_derivation` |
| C04 surface Ekman layer $\frac{d^2V}{dz^2}=\frac{if}{\nu_v}V$ (13.27), $\delta=\sqrt{2\nu_v/f}$ (13.29) | CORE | `ekman_surface`, `ekman_surface_printed`, `ekman_depth`, `ekman_residual`, `ekman_vorticity_balance`, `ekman_solve` | V2 (D06), V1, V3, V7 | residual < 2e-5 f|V|; 45° to the right / left exactly; `ekman_solve` order 2.00 | symbolic, analytic, converged | `test_ekman_surface_V2_derivation`, `test_ekman_surface_V1_spiral_residual_stress_and_hemispheres`, `test_ekman_solve_V3_second_order_and_limits` |
| C05 Ekman transport $\int_{-\infty}^0u\,dz=0$, $\int_{-\infty}^0v\,dz=-\frac{\tau}{\rho f}$ (13.30); pumping (ours) | CORE | `ekman_transport`, `ekman_transport_partial`, `ekman_pumping`, `ekman_pumping_from_curl`, `ekman_pumping_bottom`, `sverdrup_transport`, `coastal_upwelling` | V2 (D07, D09), V1 (quadrature), V3, V4 (transport independent of K(z)), V7 | quadrature to 1e-9; pumping stencil order 2.0; variable-K transport to 2e-4 | symbolic, analytic, converged | `test_ekman_transport_V2_derivation`, `test_ekman_transport_V1_quadrature_partial_and_hemispheres`, `test_ekman_pumping_V2_derivation` |
| C06 bottom Ekman layer $u=U[1-e^{-z/\delta}\cos(z/\delta)]$, $v=Ue^{-z/\delta}\sin(z/\delta)$ (13.41) | CORE | `ekman_bottom`, `ekman_bottom_transport`, `ekman_force_balance`, `ekman_finite_depth`, `eddy_viscosity_from_depth` | V2 (D08), V1, V4 (three forces close), V7 | overshoot 1.067020 at 3π/4; ½Uδ to 1e-8; force sum < 1e-15 | symbolic, analytic | `test_ekman_bottom_V2_derivation`, `test_ekman_bottom_V1_overshoot_transport_angle_forces_hemispheres` |
| C07 shallow-water set (13.44)–(13.45) | CORE | `ShallowWater`, `continuity_tendency`, `momentum_tendencies`, `step`, `run`, `energy`, `volume`, `dt_limit` | V2 (D10), V3, V4 | space order 1.98, time order 2.98; volume to round-off; dE/dt of the space discretisation < 1e-8 | symbolic, converged, conserved | `test_shallow_water_set_V2_derivation`, `test_cgrid_V3_poincare_wave_second_order_in_space_third_in_time`, `test_cgrid_V4_volume_and_energy_conserved_by_the_space_discretisation` |
| C08 vertical modes (13.56), (13.64), (13.65), $\tan\frac{NH}{c_n}=\frac{c_nN}{g}$ (13.69), $c_n=\frac{NH}{n\pi}$ (13.71) | CORE | `vertical_modes`, `vertical_modes_shooting`, `modes_uniform_N`, `orthogonality_matrix`, `project`, `reconstruct`, `w_structure`, `rho_structure`, `rigid_lid_error`, `modal_amplitudes`, `wkb_mode_speed` | V2 (D11), V1, V3, V4 | finite volume order 2.000–2.001 per mode; three solvers agree (2.5e-4, 1e-6); orthogonality 1e-12 | symbolic, analytic, converged | `test_vertical_modes_V2_derivation`, `test_vertical_modes_V1_orthogonality_weight_one_exact_quadrature`, `test_vertical_modes_V3_finite_volume_second_order_and_cross_methods`, `test_vertical_modes_V1_structures_projection_and_amplitudes`, `test_vertical_modes_V7_first_root_is_small_printed_statement_fails` |
| C09 the cubic (13.76) from (13.75) | CORE | `shallow_water_omega`, `shallow_water_discriminant`, `shallow_water_branches`, `shallow_water_regime`, `dispersion_term_sizes` | V2 (D12, D13), V1, V4 (Vieta), V7 | see claim 1 | symbolic, analytic | `test_v_equation_V2_derivation`, `test_dispersion_cubic_V2_derivation`, `test_shallow_water_omega_V1_roots_against_numpy_and_invariants`, `test_shallow_water_omega_V7_limits_poincare_kelvin_rossby_and_slip7` |
| C10 Poincaré waves (13.80), $\omega^2=f^2+gHK^2$ (13.82) | CORE | `poincare_omega`, `poincare_amplitudes`, `poincare_fields`, `poincare_group_velocity`, `poincare_orbit`, `inertial_oscillation`, `inertial_radius` | V2 (D14), V1, V3 (C-grid run), V7 | fields satisfy the linear set to the differencing error; $c_pc_g=c^2$ to 1e-12; clockwise for f > 0 | symbolic, analytic | `test_poincare_V2_derivation`, `test_poincare_V1_fields_satisfy_the_linear_set_and_orbits` |
| C11 Kelvin wave $\eta=\eta_0e^{-fy/c}\cos k(x-ct)$ (13.87) | CORE | `kelvin_wave`, `kelvin_residuals`, `kelvin_decay_side`, `kelvin_omega`, `kelvin_state`, cache `kelvin_basin` | V2 (D15), V1, V3 (model), V7 | residuals < 1e-7; model returns after one period to 0.24 % (both hemispheres); basin pulse speed = $\sqrt{gH}$ within the grid resolution, offshore scale = Λ to 5 % | symbolic, analytic | `test_kelvin_V2_derivation`, `test_kelvin_V1_residuals_trapping_and_hemispheres`, `test_cgrid_V1_balanced_states_vorticity_and_kelvin_wave_both_hemispheres`, `test_cached_runs_V4_committed_files_show_the_physics` |
| C12 Rossby radius, geostrophic adjustment (ours) | CORE | `rossby_radius*`, `geostrophic_adjustment_1d`, `adjustment_energy`, `linear_1d_step`, `linear_1d_run` | V2 (D16), V1, V3, V4 | see claim 5 | symbolic, analytic (ours), converged, conserved | `test_geostrophic_adjustment_V2_derivation`, `test_geostrophic_adjustment_V1_end_state_energy_and_jet_sign`, `test_linear_1d_V3_first_order_in_time_second_in_space`, `test_linear_1d_V4_adjustment_keeps_pv_and_reaches_the_end_state` |
| C13 potential vorticity $\frac{D}{Dt}\big(\frac{\zeta+f}{h}\big)=0$ (13.94); flow over a step | CORE | `potential_vorticity`, `step_vorticity`, `flow_over_step`, `sw_potential_vorticity`, `relative_vorticity`, `advect_particles`, `particle_potential_vorticity`, `interpolate` | V2 (D17, D18), V1, V4 | every line (13.91)–(13.94) follows from the previous one; PV on particles: 0.66 % drift of a 44 % spread (48²), 1.5 % (32²) | symbolic, analytic, conserved (approximately) | `test_potential_vorticity_V2_derivation`, `test_flow_over_step_V2_derivation`, `test_flow_over_step_V1_coded_streamlines_and_permanent_shift`, `test_particles_V1_interpolation_advection_and_pv_on_paths` |
| C14 inertia–gravity waves (13.96), $m^2=\frac{(k^2+l^2)(N^2-\omega^2)}{\omega^2-f^2}$ (13.99), (13.112) | CORE | `inertia_gravity_*`, `wkb_vertical_structure`, `vertical_structure_solve`, `wkb_error`, `w_equation_rotating_residual`, `lee_wave_m`, `lee_wave_field` | V2 (D19, D20, D21), V1, V3 (WKB scaling), V7 (f → 0 = ch07) | see claim 10; $\mathbf c\cdot\mathbf c_g=0$ to 1e-12 | symbolic, analytic | `test_w_equation_V2_derivation`, `test_vertical_structure_and_dispersion_V2_derivation`, `test_group_velocity_inertia_gravity_V2_derivation`, `test_inertia_gravity_V1_dispersion_band_regimes_group_velocity`, `test_inertia_gravity_V1_fields_hodograph_and_w_equation_residual`, `test_wkb_V3_error_falls_as_one_over_Hm`, `test_lee_waves_V1_vertical_wavenumber_field_and_tilt` |
| C15 Rossby waves (13.117)–(13.120); Rayleigh–Kuo $\beta-\frac{d^2U}{dy^2}$ (13.124) | CORE | `rossby_*`, `stationary_rossby_wavelength`, `basin_crossing_time`, `qg_linear_evolve*`, `absolute_vorticity_gradient`, `rayleigh_kuo_criterion`, `rayleigh_kuo_eigs`, `ST.rayleigh_eigs_contour(beta=)` | V2 (D22, D23, D24), V1, V3 (two eigen-solvers), V4 (integral identity) | spectral evolution exact to 1e-12; group velocity vs differences 1e-6; see claim 9 | symbolic, analytic, converged | `test_qg_vorticity_V2_derivation`, `test_rossby_dispersion_V2_derivation`, `test_rossby_V1_frequency_group_velocity_circle_and_limits`, `test_qg_linear_evolve_V1_plane_wave_and_packet_group_velocity`, `test_rayleigh_kuo_V2_derivation`, `test_rayleigh_kuo_V3_growth_rates_against_independent_shooting`, `test_rayleigh_contour_V7_beta_default_unchanged_and_radiating_mode_limit` |
| C16 Eady problem (13.128), (13.136), (13.141), (13.142) | CORE | `eady_*`, `eady_matrix`, `eady_numeric_eigs`, `eady_basic_state` | V2 (D25, D26, D27), V1, V3, **V5** | see claim 4; 2 × 2 system derived from w′ = 0 equals the page's and the coded matrix; roots of its determinant equal (13.141) to 1e-24 | symbolic, analytic, converged, benchmark | `test_eady_equation_V2_derivation`, `test_eady_phase_speed_V2_derivation`, `test_eady_cutoff_and_fastest_wave_V2_derivation`, `test_eady_V5_published_growth_rate_and_wavenumbers`, `test_eady_V3_chebyshev_eigenvalues_against_the_closed_form`, `test_eady_V1_mode_structure_vertical_velocity_and_fluxes` |
| C17 two-dimensional turbulence (13.143)–(13.145) | CORE | `fjortoft_transfer`, `enstrophy_spectrum`, `two_d_cascade_spectrum`, `rhines_length`, `barotropic_*`, `spectral_centroids`, `random_vorticity`, caches | V2 (D28, D29), V1, V4, V3 (drift order) | energy drift 2.1e-6, enstrophy 1.0e-5 over t = 1 (order 5.0 and 4.9 in dt); ⟨ψ, rhs⟩ and ⟨ζ, rhs⟩ = 0 to 1e-17 | symbolic, analytic, conserved; cascade pictures `qualitative` | `test_fjortoft_V2_two_constraints_and_ratios`, `test_cascade_spectra_V2_derivation`, `test_barotropic_V4_energy_and_enstrophy_and_strict_dealiasing`, `test_barotropic_V1_single_modes_rossby_wave_viscous_decay_and_diagnostics` |
| NOTE: stratification (13.1), profiles, standard atmosphere, lapse rate in both conventions | NOTE | `ocean_density_gradient_budget`, `ocean_profile_idealized`, `thermocline_N2`, `buoyancy_frequency_sq`, `atmosphere_layers`, `lapse_rate_table`, `scale_height` | V1, V3, V7 | N² stencil order 2.0; identical verdict in both conventions for four gradients (`Gamma_a` is a required keyword) | analytic | `test_stratification_V1_*`, `test_lapse_rate_V7_both_conventions_same_verdict` |
| NOTE: tables and names | NOTE | `illustrative_inputs`, `conventions_table`, `book_slips`, `traps`, `wind_from_to`, `viscous_layer_thicknesses`, `sympy_engines`, `reexport_audit` | V7 | 21 inputs; 14 slips (each with `kind`: "loose" for rows 10 and 11, "slip" otherwise; row 9 lists four cross-references); 18 traps; 152 re-exports | — | `test_tables_V7_*`, `test_contract_V7_*` (4 tests) |

## Functions used by the notebook and explainers (design Part C)
`test_contract_V7_every_part_c_name_exists_and_is_exercised_here` holds the list of 182 contract names (C.1–C.4) plus 23
further public names the scripts and caches use, asserts that each exists on `ch13`, and that each appears as
`ch13.<name>` / `GFD.` / `VM.` / `SW.` in the body of a test. `test_contract_V7_signatures_match_design_part_c` compares 21
signatures and 6 keyword-only parameters with the design; `test_contract_V7_scalars_in_floats_out_and_f_zero_raises`
checks the return-shape convention (16 scalar calls return `float`, f = 0 raises `ValueError` in 8 functions,
`inertial_period(0)` is the one documented `inf`); `test_contract_V7_reexports_are_the_same_objects` checks the 152
re-exports. Figure helpers (8) and sketch helpers (12): `test_figures_V7_*`, `test_drawings_V7_*`. **Untested Part C
functions: none.**

## Convergence studies

| Scheme | Grids | Observed order | Design order |
|---|---|---|---|
| C-grid shallow water, Poincaré wave, one period (`run`) | 16², 32², 64² (both hemispheres) | 1.98 (errors 14 %, 3.6 %, 0.90 % of the amplitude) | 2 |
| SSP-RK3 time stepping (`step`) at fixed 16² grid | dt = dt_max/1, /2, /4, /8 | 2.98 (solution), 2.90 (energy drift) | 3 |
| Stability limit of `step` | dt_max = √3/ω_max | stable at 1.00 dt_max, `ValueError` at 1.01 | as stated |
| 1-D forward–backward (`linear_1d_run`), time | dt, dt/2, dt/4, dt/8 at 64 cells | 1.0 | 1 (docstring: "first order in time") |
| the same, space (time error removed by Richardson) | 16, 32, 64 cells | 2.0 | 2 |
| Finite-volume vertical modes (`vertical_modes`, "fd") | 41, 81, 161, 321 nodes | 2.000, 2.000, 2.001 (free, n = 1–3); 2.000–2.001 (rigid, n = 1–4); barotropic speed already at 6e-9 | 2 |
| `ekman_solve` | 81 … 641 nodes | 2.00 (constant K, vs the spiral); 2.00 (K(z), vs a 5121-node run) | 2 |
| `geostrophic_from_field` / `ekman_pumping` stencils | 24 … 192 / 33 … 129 | 2.0 / 2.0 (edges included) | 2 |
| `buoyancy_frequency_sq` | 41 … 321 | 2.0 | 2 |
| `geostrophic_state` balance | 32², 64², 128² | 2.0 | 2 |
| `eady_numeric_eigs` (Chebyshev) | n = 4, 6, 8, 12 | error falls faster than 30× per step, < 1e-8 U₀ at n = 12 | spectral |
| `wkb_error` against 1/(H m) | H m = 3.7, 7.5, 15, 30 | 1.0; error × H m → 0.230 | 1 |
| `barotropic_run` invariant drift against dt | dt, dt/2, dt/4 | 5.0 (energy), 4.9 (enstrophy) | ≥ 4 (docstring says dt⁴ — see Open item 3) |
| `rayleigh_kuo_eigs` vs shooting | N = 120 box 16 vs adaptive ODE | agree to 1e-7 for β = 0 … 0.6 | — |

## Conservation / invariant residuals

| Quantity | Residual |
|---|---|
| Volume, C-grid (periodic, channel, closed; linear and nonlinear) | Σ ∂η/∂t ≤ 1e-13 of Σ|∂η/∂t|; after one wave period 1e-12 of amplitude × area |
| Energy, C-grid space discretisation (random state, three boundary types, β ≠ 0, linear and nonlinear) | dE/dt ≤ 2e-10 (relative) — the Coriolis/vorticity average does no work |
| Energy drift of the time stepping, one Poincaré period | −8.4e-3, −1.2e-3, −1.6e-4 on 16², 32², 64² |
| Discrete linear PV in `linear_1d_step` (δ_x v − f η̄/H) | 3e-14 relative after six inertial periods (round-off) |
| Volume in `linear_1d_run` | ≤ 1e-12 |
| Energy and enstrophy, `barotropic_run`, inviscid, 64², t = 1 | −2.1e-6 and −1.0e-5; −6.5e-8 and −3.5e-7 with dt/2; without the 2/3 rule the run blows up |
| ⟨ψ, rhs⟩ and ⟨ζ, rhs⟩ of `barotropic_vorticity_rhs` (β = 3) | −1.3e-17 and 2.1e-18 |
| PV on marked particles (nonlinear vortex) | 1.5 % (32², one inertial period), 0.66 % (48², three periods) of a 44 % spread |
| Ekman transport with variable K(z) | equals −iτ/(ρf) to 2e-4 (trapezoid on 5121 nodes) |
| Ekman-layer force balance (three vectors) | sum ≤ 1e-15 |
| Rayleigh–Kuo integral identity on the shooting mode | 2e-6 of ∫|…| |

## Benchmarks used

| Value | Ours | Source | Verified |
|---|---|---|---|
| Eady growth coefficient 0.3098 | 0.30982 | K. A. Emanuel, MIT OpenCourseWare 12.803 (Fall 2009), Lecture 19, text after its Eq. (19.17) — **read first-hand** (PDF page 7); second source: NCL `eady_growth_rate` documentation page — read first-hand (the papers it cites, Eady 1949 and Lindzen & Farrell 1980, were **not** read) | 2026-10-07 |
| Fastest wavenumber 1.606 | 1.60612 | Emanuel, same page — read first-hand | 2026-10-07 |
| Cut-off, as the source rounds it (two digits) | 2.39936 | Emanuel, same lecture (pages 7–8) — read first-hand | 2026-10-07 |

Stored by `reference/ch13/make_refs.py` in `reference/ch13/benchmarks.json`; cited in `reference/ch13/SOURCES.md`.
Not used as benchmarks (and labelled so in SOURCES.md): the Earth's rotation rate and radius (constant consistency), the
adjustment energy split (Gill 1982 unread), β-plane jet growth rates (Kuo unread), the −3 spectrum and the Rhines length
(Kraichnan, Rhines unread).

## Numbers from the text (private values redacted to relative differences)
Four V6 tests read `tests/book_values_ch13.json` (skipped when it is absent; none of its numbers is typed in the test
file — `test_public_V7_no_private_book_number_is_typed_in_this_file` and `tools/check_public.py` enforce it).

| What | Result |
|---|---|
| Rotation rate and polar f as the page rounds them (solar day) | within 0.4 % and 0.3 %; sidereal value 0.27 % larger (trap T1) |
| Rossby-number example of §13.5 | exact |
| Ocean adiabatic density gradient, scale heights | 9 % and < 20 % (one-digit and "about" values) |
| Laminar Ekman thickness of §13.7 | the formula gives a value 37 % larger than the page prints (slip #10 confirmed); the eddy viscosity implied by the observed layer: exact |
| Uniform-N modes with the chapter's typical values | exact roots to 1e-6 of the analyst's; the first root NH/c₀ is a few hundredths, not 1 (slip #6); the page's first baroclinic speed is the formula's value rounded up by a quarter |
| Regime ratios of §13.10 | 1 %; inertial-circle example exact; deep-sea Kelvin example within 0.7 % (two-digit rounding) |
| Internal Rossby radius of §13.12 | the page's typical value is three times $NH/(\pi f)$ with the chapter's own N, H, f (slip #11) |
| Rossby-wave maximum frequency, long-wave speed | exact against the analyst's values; "shortest period more than a year" comes out 0.4 % short of a year with the page's round inputs |
| Exercises: sea-surface slope, Ekman transport, tank Ekman layer, tidal Poincaré wave, coastal Kelvin wave, baroclinic Rossby speed | 0.1 %, 0.05 %, 0.2 %, 0.2 %, 0.5 %, 0.05 %; the printed Rossby speed needs the round β (0.2 %) and is 19 % above the value with β at the stated latitude (slip #12) |
| Eady cut-off and the two wavelength multiples; Fjørtoft's worked triad | 0.03 %, 0.7 %, 0.3 % (two-digit rounding); exact |

## Printed slips — confirmed or refuted (each page re-read as an image unless marked)

| # | Printed | Verdict | Test |
|---|---|---|---|
| 1 | (13.2) with $+\frac1{\rho_0}\nabla p$ | **Slip.** A fluid at rest would accelerate at $2g\bar\rho/\rho_0$; the next two displays of the same page carry the minus | `test_boussinesq_V2_rest_state_printed_sign_fails` |
| 2 | $F_i=\partial\tau_{ij}/\partial x_j$ as a force per unit mass | **Slip** (dimensions: N/m³ against m/s²) | `test_friction_force_V2_dimensions_printed_form_fails`, `test_friction_force_V2_derivation` |
| 3 | (13.17) $-2\Omega u=-\frac1\rho\frac{\partial p}{\partial y}$ | **Slip.** Contradicts (13.12) at f = 2Ω; cross-differentiation then gives $\partial u/\partial x-\partial v/\partial y=0$ and the divergence $q_{xy}/(\Omega\rho)$ survives | `test_taylor_proudman_V2_derivation` |
| 4 | (13.26) "as z → ∞" | **Slip.** The ocean is z < 0; the root bounded upward grows without limit downward; the page itself says z → −∞ two lines later | `test_ekman_surface_V2_derivation` |
| 6 | "the first root occurs for NH/c_n = 1" | **Slip** (≪ 1): X = 1 misses (13.69) by more than 1; the root is $N\sqrt{H/g}$ = 0.056 for our inputs and tan x ≈ x holds there | `test_vertical_modes_V7_first_root_is_small_printed_statement_fails` |
| 7 | "first term … negligible for ω ≫ f" | **Slip** (≪): on the slow root ω³ is the smallest term (< 1e-3 of the β term); at ω = 3f it is the largest | `test_shallow_water_omega_V7_limits_poincare_kelvin_rossby_and_slip7` |
| — | `ekman_surface_printed` (cos / sin components) | **True as printed for f > 0 — a trap, not a slip.** Equals `ekman_surface` to 1e-15 for f > 0 and refuses f ≤ 0 | `test_ekman_surface_V2_derivation` |
| 8 | "let N be depth independent" (§13.14 opening) | Slip confirmed from the text (the sentence announces more generality than ch07, and (13.99) writes N²(z)); no code variant | — |
| 10, 11, 12 | numerical inconsistencies | confirmed in the private tests (table above) | `test_book_V6_*` |
| 5, 9, 13 | wording, cross-references, attribution | not re-read by the verifier; nothing coded depends on them | — |
| 14 | westward flow over a step "at its original latitude" (§13.13) | **Slip, confidence about 90 %** (claim 6); **now row "14" of `ch13.book_slips()`** | `test_flow_over_step_V2_derivation`, `test_flow_over_step_V1_coded_streamlines_and_permanent_shift`, `test_flow_over_step_V4_printed_westward_statement_fails_pv_conservation` |

`ch13.book_slips()` has 14 rows and `ch13.traps()` 18 (T18: the sign wording of the vertical Coriolis term); the tables test asserts both and the `kind` of every row.

## Mutation check (do the tests discriminate?)
The package was copied to the scratchpad and one likely-wrong variant planted at a time (34 in all): transport sign,
sign(f) dropped in the Ekman root, thermal-wind sign, sign of the β term of the cubic, Kelvin direction, jet without
sign(f), adjustment energy coefficient, factor ½ in the Eady speed, surface term dropped from the mode solver, a surface
term added to the weight-1 inner product, Rossby group-velocity term, non-strict de-aliasing, one-sided Coriolis average
in the 1-D step, group velocity without f², westward step returning to its latitude, β sign in the Rayleigh solver,
Poincaré amplitude sign, bottom pumping without sign(f), C-grid Coriolis sign, m² with the wrong sign, δ without the 2,
geostrophic sign, Eady matrix sign, β with sin for cos, step vorticity denominator, Fjørtoft's two energies swapped,
Rossby frequency sign, Poincaré frequency without f, surface height without g, QG evolution sign, bottom transport sign,
growth rate from c_r, inertial-circle sense, lapse-rate operators swapped. **34 of 34 were killed in the first pass; after the review the check was re-run on the changed code with five more (no projection of the start in `barotropic_run`, `poincare_orbit` path flipped in y, `inertial_oscillation` path flipped in y, north-wall `kelvin_state` travelling the wrong way, north-wall state decaying from the south wall): 39 of 39 killed.** One survived the
first run in one of its two tests — `numpy.allclose` has a default absolute tolerance of 1e-8, larger than β itself — so
every comparison in the file now goes through a helper with **no** absolute tolerance unless one is written.

## Figures reproduced with our code (`tests/ch13_verify_figures.py` → `outputs/ch13/verify/`, local)

| Figure | Visual verdict |
|---|---|
| `fig_ekman_layers.png` | Surface hodographs start 45° to the right (f > 0) and left (f < 0) of the stress and spiral inward as mirror images; the bottom hodograph leaves the origin along the 45° line and winds onto (1, 0); u overshoots to 1.067 U at z = 3πδ/4 and v first vanishes at πδ. |
| `fig_vertical_modes.png` | tan X crosses the falling curve once just above 0 and once just above each nπ; finite-volume dots sit on the exact uniform-N modes with n interior zeros; with a thermocline the baroclinic modes change sign within the upper 500 m. |
| `fig_dispersion.png` | The fast root lies on $\sqrt{f^2+c^2K^2}$, flat at ω = f for KΛ < 0.3 and joining ω = cK beyond KΛ ≈ 10; the slow root lies on the Rossby curve four decades below, with its maximum at KΛ = 1. |
| `fig_adjustment.png` | The model's inertial-period mean lies on $\eta_0\,\mathrm{sgn}(x)(1-e^{-|x|/\Lambda})$ and on the jet $e^{-|x|/\Lambda}$ peaking at gη₀/c; the f < 0 jet is its mirror image; the instantaneous field at half an inertial period shows grid-scale ripples behind the two fronts (see Open item 6). |
| `fig_eady.png` | Growth rises from 0, peaks at (1.606, 0.310) and falls steeply to 0 at 2.399, with the Chebyshev eigenvalues on the curve; c_r = U₀/2 while growing, then two branches; c_i starts at 0.289 = 1/(2√3); the mode's amplitude is symmetric and its phase rises by π/2. |
| `fig_step_and_rayleigh_kuo.png` | Eastward flow: undisturbed upstream, then a meander between 0 and 2Y_p; westward flow: a monotonic bend that begins before the step, passes Y_p/2 at the step and settles on Y_p; growth rate of the jet falls almost linearly with β and vanishes just below β = 2/3. |
| `fig_models.png` | The Kelvin pulse hugs the south wall, moves toward +x, turns up the east wall and along the north wall (coast on its right); enstrophy decays to 0.21 while energy stays at 0.73 and the energy centroid falls to 0.37 (0.41 with β); PV on particles stays within ±0.7 %. |

## Deviations & justifications
`fluidpy/core/shallow_water.py` carries the chapter's one `DEVIATION` (method only): **every numerical scheme is ours —
the book contains no numerical method** (C-grid with an energy-conserving Coriolis average and SSP-RK3; the
forward–backward 1-D line; pseudo-spectral RK4 with the 2/3 rule). Orders and invariants are measured above. The
corrected forms of slips #1–#4, #6, #7 are implemented, with the printed variants kept (`printed=True`,
`friction_force_dimensions(printed=True)`) and shown to fail. `OMEGA_EARTH` is the sidereal value; the page's solar-day
value is `OMEGA_SOLAR_DAY` (0.27 % smaller). Lapse rates are computed in Kundu's sign with `Gamma_a` a required keyword;
`lapse_rate_table` shows both conventions with the same verdict.

## Open items
1. **`rayleigh_kuo_eigs` at negative β is a box-truncation result, not a converged one** (k = 0.9, β = −1: 0.0546, 0.0599,
   0.0565 for (N, y_max) = (120, 16), (200, 24), (280, 32)). Cause, measured: the far field of this mode is a weakly
   damped radiating Rossby wave, $\kappa=\sqrt{k^2+\beta/c}$ with Re κ = 0.07, so φ = 0 at ±y_max is the wrong condition
   until y_max ≫ 14. An independent shooting integration with the decaying far-field solution gives
   **c = 0.883848 + 0.063886 i, growth 0.057498**, the same to 1e-10 from starting points 14 and 24. The docstring's
   "0.055 against 0.060" compares two box values; neither is the answer. Label for β < 0: `approximate` (5 %). Suggested
   (documentation only; `fluidpy/ch13_geophysical_fluid_dynamics.py::rayleigh_kuo_eigs`): state the cause and the value,
   and that the wrapper is validated for β ≥ 0. The notebook and explainers use β ≥ 0 only (design Part A/B) — not blocking.
2. **`vertical_modes(method="cheb")` loses the barotropic speed**: 8.5e-5 relative at `n_cheb=96` on the thermocline
   profile (baroclinic speeds 4e-8). The docstring already says the non-symmetric solve loses digits (2e-7 at 64 on
   uniform N); the measured number at 96 is larger. Use "fd" or `vertical_modes_shooting` for c₀. Not blocking.
3. **`barotropic_run` docstring understates the scheme**: it says the invariant drift is ∝ dt⁴; measured order 5.0
   (energy) and 4.9 (enstrophy). Harmless.
4. **`analytic (ours)` without a first-hand source**: the adjustment end state and its 1/3 (Gill 1982 unread); β-plane
   jet growth rates (computed twice, Kuo unread); the −3 range and the Rhines length (dimensional analysis only).
5. **`qualitative`**: the cascade pictures of `turbulence_f.npz` / `turbulence_beta.npz` (64², viscous): energy centroid
   falls 8.6 → 3.2, zonal energy fraction 0.11 without β and 0.29 with it — illustrations, not measurements of a spectral
   slope.
6. **For the notebook builder**: the instantaneous field of `linear_1d_run` started from a discontinuous step carries
   grid-scale ripples of up to about half of η₀ behind the fronts (non-dissipative scheme, visible in
   `fig_adjustment.png`); the design's "a few per cent" holds for the inertial-period mean near the step, not for a
   snapshot. Show the mean, or smooth the initial step over a few cells.
7. Slips #5, #9, #13 (wording, cross-references, attribution) were not re-read by the verifier.

## Found wrong in design / analysis (for the designer and curator; no code change)
- `analysis/ch13_design.md` Part C, row `orthogonality_matrix(modes)`: "(plus the surface term for `lid="free"`)" is
  wrong for the default `kind="psi"` and contradicts Part F D11 step 13 and the code. The surface term belongs to
  `kind="energy"` only (claim 2).
- `analysis/ch13_curation.md` §4b, D11 traps: "the free surface adds a boundary term to the orthogonality relation" —
  the same error; it adds one to the energy relation.
- (Settled after the review: the westward-step statement is now row "14".) Before: it was not in `ch13.book_slips()` (rows "1" … "13"); it should become a slip row or,
  if the user prefers caution, a trap with the 90 % confidence stated.
- The task brief lists `ekman_surface_printed` among the printed variants that must fail; it is true as printed for
  f > 0 (a trap). The test asserts exactly that.
- Design "expect" numbers checked and found right: D11 (3.61 m/s), D16 (Λ = 43.15 km, jet 0.136 m/s), D18 (Y_p = −223 km,
  wavelength 5983 km, e-folding 952 km), D27 (0.3098). D16's "0.2 % and 0.4 % over the fifth inertial period" was not
  reproduced at our resolution (0.5 % and 0.6 % over the sixth period with cells of 0.2 Λ) — a hypothesis about a
  different grid, not an error.

## After the code review (`reports/ch13_review.md`) — what changed in the tests
| Change in `fluidpy` | Test evidence added | Measured |
|---|---|---|
| `book_slips()` 14 rows with `kind`; `traps()` 18 rows | `test_tables_V7_inputs_conventions_slips_traps_and_names` asserts the 14 keys, every `kind`, the four cross-references of row 9, T1–T18; new `test_flow_over_step_V4_printed_westward_statement_fails_pv_conservation` (independent of the solver: "back at its latitude with ζ = 0" gives PV × h₀/h₁, three depth pairs; `flow_over_step` far downstream sits on Y_p with PV kept to 1e-12) | the printed statement misses (13.94) by 5 %, 7 % and 48 % for our three steps |
| `barotropic_run` projects the start onto the retained modes | new `test_barotropic_V4_white_noise_start_is_projected_and_conserved`: white noise plus a large-scale mode, and a sharp-edged top-hat (both with > 2 % of their spectral power outside the mask), 64², inviscid, t = 4, quarter step | energy −6.2e-9, enstrophy −1.8e-7 (noise case; the implementer's −3.6e-8 / −1.9e-7 is the same size for a different field); drift order 4.92 / 4.97; frames stay inside the mask to 1e-9; the mutant without the projection fails |
| reviewer's surviving mutants | `poincare_orbit` and `inertial_oscillation` now checked for dy/dt = v, the path's axis ratio and a closed circle about its centre; `kelvin_state(wall="north")`: largest on the north wall, envelope exp(−(L_y − y)/Λ) to 1e-9, travels toward −sign(f)·x (tendency and quarter-period shift), returns after a period to 0.4 %, cross-channel geostrophy to 2e-3; `kelvin_residuals` at two more points | both y-flipped mutants and both north-wall mutants die |
| public-repo hygiene | no public test types the page's typical f, β or f/ω pair any more: public tests use `F_N`, `F_S` (35°, ours) and β = 1.7e-11 or 1.9e-11; the two book-only tests read f and β from the private JSON | `tools/check_public.py` OK on the five files |
| docstring numbers | `ekman_residual`: 4001 nodes over 14 δ, asserted within 0.7e-6 … 2.1e-6 of f|V| (docstring 1.4e-6); `absolute_vorticity_gradient`: 401 points over |y| ≤ 5, asserted within 4e-4 … 1.3e-3 (docstring 8.3e-4) | both inside their bands |

Open item 3 of the first pass (the docstring of `barotropic_run` said dt⁴) is superseded by the measured order above; the
other open items stand.

## Rest of the suite
Full suite after the review, `.venv/Scripts/python.exe -m pytest tests -q -p no:cacheprovider -m "not slow"`: **1865 passed, 0 failed, 19 deselected (2994 s on a machine shared with other runs)**.
`tests/test_machinery.py` is part of that run. First pass: 1776 passed in the other files, `test_machinery.py` 9 passed.
`fluidpy/core/stability.py` (`rayleigh_eigs_contour` gained `beta=0.0`): the default path is unchanged — identical arrays
with and without the keyword, equal to the independent shooting value to 1e-7, and all of chapter 11's tests pass.
