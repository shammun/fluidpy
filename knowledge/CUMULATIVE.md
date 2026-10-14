# CUMULATIVE knowledge — fluidpy
(Rewritten by the knowledge-keeper after every chapter. Every agent reads this first. Last rewrite: after ch16 — the
book's last chapter — 2026-10-14, knowledge pass on commit a042a3c.)

## Open items after ch16 — read before a revision pass or the library pass (details `knowledge/ch16.md` §0, §9)
- **Needs the user**: ch13 slip #14 (westward flow over a step) still has no ruling.
- **Code**: ch16 sibling argument order deferred — `vascular_resistance(mu, L, a)` vs `poiseuille_flow_rate(a, dp, L, mu)`,
  `tube_radius_hookean` vs `tube_law_hookean` (call by keyword; also ch14 `aspect_ratio`/`mean_chord`, ch15
  `ackeret_pressure`/`weak_wave_pressure`); FL closure magnitude qualitative; S21 unresolved; carried ch15 items (§0 there).
- **Tooling**: `tools/shot.py` lacks `.katex-error` / control-character checks (the ch16 Must class), real clicks,
  page-then-step, Explain state change, 844×345, preset exclusivity, Code-tab runs, a per-file write retry; `eq_refs`
  skips derivation `start`/`result`; `testpaths = tests`; delete `outputs/ch14/_resume_scratch/mut/`.
- **Text**: Should-fix lists of ch14–ch16 (viz and lesson) in `knowledge/viz_patterns.md`. **Unread sources**: Dean 1928,
  Womersley 1955/57, Shapiro 1977a, Wagner/Jones, Gill, Kuo.

## Project rules in force (added during ch01–ch03; sharpened since)
- **5–10 explainers per chapter, as many as the CORE ideas need** (ch03 7, ch04–ch09 9, ch10 8, ch11 9, ch12 10, ch13 10,
  ch14 9, ch15 9, **ch16 10**).
- **Every book equation is shown in full next to its number** — notebook, explainer text, `<meta>` strings and the design
  file; **a number stands only beside the equation actually printed under it** (render the page); never start a sentence
  with a symbol. Enforced by `tools/coverage_check.py` check 8, `tools/eq_refs.py`, the builder's `self_check_near`.
- **No control characters, no backslash-less TeX** (`coverage_check` check 9 — notebooks only); **every TeX-bearing JS
  string in an explainer is `String.raw`** (ch16); **write files with Edit/Write only.**
- **Rule 9 (public repo)**: book-quoted numbers live only in `tests/book_values_chNN.json` (`_forbidden_public` strings and
  regex patterns, enforced by `tools/check_public.py` on tracked or named files). "Illustrative" constants must be visibly
  not the book's; **a book value written in words ("a quarter below") is still a book value — the tool cannot see it** (ch13).
- **Pager and audit**: `tools/shot.py` pages every pager at every size; `build_chNN.py --dump` is the reference for
  derivation steps. **Two values, the book's first** (ch09); a **DEVIATION box** wherever the code departs from the book.
- **Eigen-solver contract (ch11)**: leading mode or raise; box and N doubling; an independent route for every headline
  eigenvalue. **Contract discipline (ch12, ch13)**: **design Part C first, reconciled mechanically**; design "expect"
  values are hypotheses; convention-carrying arguments and empirical constants are required keywords; an unread source
  is not a benchmark.
- **Figure notes are written from the rendered figure, never from the storyboard** (ch13). **Slips and traps are kept
  apart**; a slip row carries `kind` ("slip" or "loose"); items not in the book are labelled "ours".

## Physics pipeline so far
```
ch01 Introduction — continuum ρ, p, T, u (D35, Kn); τ = μ du/dy, FTCS [C12]; statics, buoyancy, USSA-1976; first law →
     Gibbs → c² → p = ρRT → p/ρ^γ; N², Γ_a, θ (Kundu Γ ≡ dT/dz in code, meteorology −dT/dz shown); Π = null space
ch02 Cartesian tensors — summation, C_ij, Cauchy f_i = τ_ji n_j, τ' = CᵀτC, invariants, ε–δ, S + A, principal axes, ∇ on
     grids, Gauss, Stokes
ch03 Kinematics — DF/Dt (3.5); stream/path/streak; Galilean; du = G·dx → S, ½ω; Rankine, Gaussian vortices; RTT (3.35)
ch04 Conservation laws (the trunk) — mass, ψ; momentum, Bernoulli, Cauchy; τ = −pδ + 2μS + λS_mmδ → NS; rotating frame
     ±2Ω×u′; energy, ε ≥ 0; Lamb, unsteady Bernoulli; Boussinesq; interfaces; dimensionless NS
ch05 Vorticity — tubes, Kelvin, baroclinic torque, Helmholtz; Dω/Dt (5.13); Biot–Savart; rotating baroclinic (5.30);
     (ζ + f)/h; point vortices, images, sheets
ch06 Ideal flow — Laplace, superposition; cylinder C_p, D = 0, L = ρUΓ; w = φ + iψ; Blasius, Kutta–Zhukhovsky; conformal
     maps; Laplace on a grid (D22, D23); Stokes ψ; added mass
ch07 Gravity waves — plane wave + BCs → ω² = gk tanh kH (7.28); orbits, energy, c_g, rays, jump, KdV, Stokes drift,
     interfaces √(g′H), internal waves ω = N cos θ
ch08 Laminar flow — Couette–Poiseuille, pipe f = 64/Re, circular Couette; lubrication (Reynolds equation, slider, thin
     film); Stokes' problems η = y/√(νt), δ_e = √(2ν/ω); similarity by exponent matching; creeping flow 6πμaU, Oseen
ch09 Boundary layers — δ ~ L Re^{−1/2}; Blasius (0.332057); Falkner–Skan fold at n = −0.090429; Thwaites (9.50);
     separation; Kármán street b/a = 0.28055; drag crisis; free jet J const, wall jet Ψ const; teacup → Ekman
ch10 CFD — stencils, order, round-off floor; FTCS, von Neumann G (10.24), Noye (10.27); upwind, CFL, Lax; Galerkin; cell
     Péclet > 2 ⇒ wiggles; MacCormack (reused by ch15 `euler1d_run`); projection on a C-grid; LBB; cavity vs Ghia 0.26 %
ch11 Instability — "basic state + small disturbance → normal modes → eigenvalue problem → neutral curve → critical number;
     the theorems are one-way": KH ΔU_min 6.70 m/s; Bénard Ra_c 1707.76 rigid, 27π⁴/4 free (11.44); Miles–Howard Ri > ¼
     (11.67); Orr–Sommerfeld Re_c 5772.22, Blasius 519.06; dE/dt = P − Λ (11.88); Lorenz r_H 24.737   (D01–D25)
ch12 Turbulence — "an average fails on exactly one thing, a product": RANS; ε̄ = 15ν mean(u²)/λ_g² (12.43); η = (ν³/ε̄)^{1/4}
     (12.50); S₁₁ = C₁ε̄^{2/3}k₁^{−5/3} (12.54); log law, DNS κ = 0.384 in a stated window; mixing length, k–ε; Ri = Pr_T Rf
     (12.109), L_M; Taylor dispersion   (D01–D28)
ch13 Geophysical fluid dynamics — "one set of equations on a thin rotating shell with terms switched off; (ζ + f)/h is what
     a column remembers": f = 2Ω sin θ (13.8); geostrophy, thermal wind; Ekman δ = √(2ν_v/f) (13.29), transport −τ/(ρf)
     (13.30); vertical modes c_n² ≡ gH_e (13.62); the cubic (13.76), Poincaré, Kelvin, Rossby waves; adjustment keeps ⅓
     (ours); Rayleigh–Kuo β − U″ (13.124); Eady σ = 0.30982 fU₀/(NH)   (D01–D29)
ch14 Aerodynamics — "lift is circulation; the trailing edge chooses it; the wing tips make you pay for it" (Γ clockwise):
     Λ ≡ s²/A (14.1); Kutta condition; starting vortex; z = ζ + b²/ζ (14.3), t_max/c ≃ 1.3e (14.8) = (3√3/4)e; Γ = 4πUa
     sin(α + β) (14.11), C_L ≃ 2π(α + β) (14.12); downwash (14.13), α_e = α − w/U (14.14); Γ = ΣΓ_n sin nγ (14.18), (14.20);
     L = (πs/4)ρUΓ₁ (14.21), D_i = (π/8)ρΣnΓ_n² (14.22); C_Di = C_L²/(πΛ) (14.26); polar, thrust from lift (ours)   (D01–D22)
ch15 Compressible flow — "three conservation laws and the gas law in one number M; 1 − M² changes sign at sound": p′ ≅ c²ρ′
     (15.5) → (15.12); d'Alembert (15.15)–(15.17); stagnation ratios (15.28)–(15.30); A/A* (15.31); (15.34) [slip #4]; normal
     shock (15.39)–(15.44); choking, Laval (ours); friction/heat quadratic (15.46), (15.47)–(15.48) [slips #11, #12]; Mach
     cone (15.49); oblique shock (15.53); Prandtl–Meyer (15.56); Ackeret (15.59)   (D01–D31)
ch16 Biofluid mechanics — "odd fluid, pulsed, through elastic, branching, curved, collapsing tubes; Poiseuille's R or a
     wall-set wave speed, modified by α": W = −∮p dV (ours), CO (16.3); Casson (16.6); α = a√(ω/ν) (16.12); R = 8μL/(πa⁴)
     (16.20); two-layer μ_app (16.33) [S3]; tube law (16.41); Windkessel (16.49), (16.52); c = √(Eh/(2ρa₀)) (16.62) [S5];
     Womersley (16.82) [S6], Q and lag (16.90)–(16.100) [S8, S10]; elastic-tube quadratic (16.136); Murray (16.140) [S15],
     (16.144); junction R (16.153); Dean c₂ 1541/50400, c₄ [S34]; Shapiro (16.204), S = 1 chokes, (16.188) only α ≤ 1 [S23];
     Casson plug (16.210), flux (16.220); ψ = p − RTc + ρgz (16.221)   (D01–D34; the book's last chapter)
```
**The book is complete** (appendices A–D are reference only). What follows: the library pass, the tooling TODOs and a
revision pass over the open Should-fix items; any capstone (e.g. a climate notebook) starts from P396 and reuses the
primitives below.

## Available primitives
### Machinery (`fluidpy/core/`, ready before chapter 1)
| module | what | used for |
|---|---|---|
| `project` | repo root, `book.yaml` config, Pages/raw/Colab URLs | every tool and notebook |
| `style` | `setup_notebook()`, `COLORS`, `savefig` | every notebook figure |
| `embed` | `show_viz(chapter, slug)` — Jupyter (srcdoc) / Colab (Pages src) / page | explainer cells |
| `anim` | `animate`, `show_animation(player="video"/"frames")` | animations (matplotlib 3.11 constrained-layout trap) |
| `interact` | `slider_figure`, `animate_figure` (plotly), `live` (ipywidgets) | Python interactives |
| `units` | pint `ureg`/`Q_`, `dimensional_check`, `check_dimensions`, `DIM`, °C ↔ K | V2 tests |
| `refdata` | reference-data registry | benchmarks |
| `_util` | `as_scalar_if_0d`, `require_positive`, `require_nonnegative` | every physics function |
| `_stencil` (ch04) | private explicit 2nd-order central differences on callables | residual functions |
| `tools/convergence` | `observed_order(h, err)`, `pairwise_orders`, `richardson`, **`grid_convergence_index(f1, f2, f3, r, p, Fs)` (ch10)**, `refinement_study`, `mms` | every V3 test |

### Physics primitives (re-exported as `chNN.<name>`; details in `knowledge/chNN.md` §3)
| function family | module | book § / Eq. | label |
|---|---|---|---|
| constants (`K_B`, `R_AIR`, `GAMMA_AIR`, `G0`, `P_ATM`, `G_BOOK` = 9.81 …); perfect gas, first law, entropy, sound speed, isentropic, van der Waals | `thermo` | §1.8–1.9, (4.94) | analytic, symbolic, converged, conserved, benchmark |
| hydrostatics, USSA-1976, buoyancy; N², `classify_stability`, lapse-rate conventions, θ; FTCS, `crank_nicolson_1d` (ch08); Π groups | `statics`, `stratification`, `diffusion`, `dimensional` | ch01, ch08 | analytic, symbolic, converged, benchmark |
| index notation, `tensors`, grids `[k, j, i]` + `operators`, `fields`, `integral_theorems` | ch02 modules | (2.1)–(2.36) | analytic, symbolic, converged, conserved, benchmark |
| `kinematics`, `coords`, `transport` (Leibniz, RTT); `vortices` | ch03/ch05 modules | (3.1)–(3.35), §5.4 | analytic, symbolic, converged, conserved, benchmark |
| `conservation`, `constitutive`, `navier_stokes` (residuals, `navier_stokes_sym(mu_v=)`, `exact_solution`), `curvilinear`, `streamfunction`, `rotating`, `bernoulli`, `interfaces`, `similarity` | ch04 modules | (4.5)–(4.119) | analytic, symbolic, converged, conserved, benchmark |
| `vorticity`, `biot_savart` | ch05 modules | (5.3)–(5.33) | analytic, symbolic, converged, conserved |
| `potential`, `conformal`, `laplace_solvers`, `panels` | ch06 modules | (6.1)–(6.109) | analytic, symbolic, converged, conserved, benchmark |
| `waves` (dispersion, signed `group_velocity`, packets, rays, two-layer and internal waves) | `waves` (ch07) | (7.1)–(7.159) | analytic, symbolic, converged, conserved, benchmark |
| `laminar`, `lubrication`, `creeping` | ch08 modules | (8.1)–(8.53) | analytic, symbolic, converged, conserved, benchmark |
| `boundary_layer` (BL, 47), `bluff_body` (BB, 22), `jets` (JET, 22), `similarity_reduce` (SR) | ch09 modules | (9.1)–(9.85) | analytic, symbolic, converged, conserved, benchmark, qualitative (labelled) |
| **`fd` (FD)**: `fd_weights(offsets, m)` (exact rationals), `fd_derivative`, `fd_apply`, `mixed_derivative_onesided`; `transport_1d_step(T, alpha, beta, scheme)` (ftcs/btcs/cn/upwind/upwind_printed/lax_wendroff), `solve_transport_1d` (guard + cap), `advected_gaussian`, `advect_periodic`; `amplification_factor(theta, alpha, beta, scheme)`, `max_amplification`, `worst_theta`, `ftcs_stable`, `stability_verdict`, `phase_error`, `ftcs2d_max_amplification`, `propagate_error`, `fourier_mode(convention=)`; `truncation_error_sympy`, `truncation_terms`, `error_norm`, `convergence_study`, `ode_scheme_order`, `richardson_three`, `grid_convergence_index`; `steady_cd_exact/fd/discrete_exact`, `discrete_root`, `wiggle_indicator`, `stretched_grid`, `cell_peclet`, `numerical_diffusivity(C=)`; `rod_heating_exact` (10.199), `lax_demo` | `fd` (ch10) | (10.1)–(10.31), (10.84)–(10.94), (10.199) | analytic, symbolic, converged, conserved |
| **`fem1d` (FEM1)**: `hat`, `hat_basis`, `interpolate`, `element_matrices_linear`, `element_integrals`, `element_force`, `connectivity(printed=)`, `assemble_1d(lumped=)`, `assembly_trace`, `solve_steady`, `solve_transport(method=)`, `interior_stencil`, `bilinear_form`, `weak_residual`, `shape_slopes(printed=)`, `gauss_legendre` | `fem1d` (ch10) | (10.32)–(10.78) | analytic, symbolic, converged |
| **`mac` (MAC)**: `MacGrid` (C-grid p[ny, nx], u[ny, nx+1], v[ny+1, nx]), `divergence`, `gradient`, `curl`, `vorticity`, `predictor(body=)`, `pressure_poisson_matrix`, `solve_pressure(splu/fft/sor)`, `pin_pressure`, `correct`, `project`, `projection_stages`, `dt_limit`, `step`, `run`, `taylor_green`, `channel_poiseuille`, `cavity`, `streamfunction`, `primary_vortex_centre`, `cavity_centreline` | `mac` (ch10) | (10.111)–(10.128) | analytic, converged, conserved, benchmark |
| **`maccormack` (MCK)**: `maccormack_step(U, flux_E, flux_F, …)`, `maccormack_advection_1d`, `lax_wendroff_advection_1d`, `ns_*`, `weakly_compressible_step`, `maccormack_dt` (10.110), `maccormack_dt_additive` (default, DEVIATION), `maccormack_dt_asymptotic` (10.155), `cavity_maccormack`, `wc_taylor_green`, `block_channel`, `block_wall_density`, `body_forces` | `maccormack` (ch10) | (10.95)–(10.110), (10.138)–(10.155) | analytic, symbolic, converged, conserved (periodic only), qualitative (block) |
| **`fem2d` (FEM2)**: `p2_shape(_grad)`, `p1_shape`, `iso_map`, `jacobian`, `tri_quad_7pt`, `integrate_element`, `Mesh`, `structured_square_mesh`, `cylinder_channel_mesh`, `assemble_saddle`, `assemble_newton_system`, `apply_dirichlet`, `newton_solve`, `stokes_solve`, `stokes_cavity`, `infsup_constant/table`, `kovasznay_test`, `poiseuille_test`, `cylinder_steady`, `march_unsteady`, `cylinder_forces` | `fem2d` (ch10) | (10.134)–(10.137), (10.156)–(10.198) | analytic, converged, symbolic, qualitative (confined St) |
| **`stability` (ST)**: `cheb`, `cheb_matrices`, `clenshaw_curtis_weights`, `cheb_grid(map="linear"\|"tan"\|"algebraic")` → `SpectralGrid`; `constrained_eig` (boundary unknowns eliminated), `apply_bc_rows` + `generalized_eigs`, `converged_mask`, `converged_eigs`; `normal_mode`, `sigma_from_c`, `c_from_sigma`, `stability_class`, `stability_verdict`, `marginal_type`; `orr_sommerfeld_eigs(k, Re, U, Upp, bc="wall"\|"decay"\|"semi_infinite")`, `os_mode(index=)`; `rayleigh_eigs` (real axis), **`rayleigh_eigs_contour(k, U, Up, Upp, delta)`** (complex path; walls/unbounded only), **`rayleigh_shoot`**; `taylor_goldstein_eigs(k, U, Upp, N2)`; **`decay_box(k)`** = max(30, 12/k), **`decay_map_scale(k)`** = min(0.5, max(0.035, 0.025/k)), `decay_box_map_scale(y_max)`, `far_field_fraction`; `max_growth`, `neutral_curve`, `critical_point`; `howard_semicircle`, `in_howard_semicircle`, `inflection_points`, `squire_transform`; `disturbance_energy_budget` | `stability` (ch11) | (11.1), (11.61)–(11.62), (11.71)–(11.72), (11.77)–(11.84), (11.88) | analytic, converged, conserved, benchmark |
| **`turbstats` (TS)** — statistics of any fluctuating signal: `make_ensemble` (members on axis 0), `ensemble_average`, `moment`, `central_moment`, `statistics(normalized=)`, `standard_error`, `time_average`, `volume_average`, `reynolds_decompose`, `product_average_split`, `check_averaging_rules`; `correlation`, `correlation_coefficient`, `correlated_pair`, `autocorrelation(method="fft", unbiased=)`, `cross_correlation`, `correlation_time`, **`integral_scale(upto="first_zero"\|"all")`**, `effective_samples`, `taylor_microscale`, `spatial_correlation`, `longitudinal_transverse_correlation`; `spectrum_from_correlation` (Filon), `correlation_from_spectrum`, `spectrum_variance`, **`periodogram(segments=, window="boxcar"\|"hann")`**, `taylor_frozen`, `frequency_to_wavenumber_spectrum`, `correlation_spectrum_pair`, `shell_spectrum`; `smooth_signal`, `synthetic_solenoidal_field(dim=2\|3)` (kinematic), `white_noise_field`; `velocity_covariance`, `reynolds_stress`, `turbulent_kinetic_energy`, `anisotropy_tensor` | `turbstats` (ch12) | (12.1)–(12.26), (12.38)–(12.39), (12.45) | analytic, symbolic, converged, conserved |
| **`wall_turbulence` (WT)** — wall and surface-layer mean-flow summaries: `friction_velocity`, `viscous_length`, `friction_reynolds_number`, `wall_units`, `from_wall_units`, `viscous_sublayer`, `log_law(yplus, *, kappa, B)`, `log_law_defect`, `friction_law_from_overlap`, `log_law_crossing`, **`fit_log_law(window=)`**, `log_law_indicator`, `layer_name`, `spalding_uplus/yplus/slope`, `coles_wake`, `composite_profile(_plus)`, `velocity_defect`, `total_stress`, `channel_total_stress`, `channel_pressure_gradient`, `pipe_pressure_gradient`, `stress_partition`, `boundary_layer_stress_from_profile`, `zpg_boundary_layer`, `skin_friction_zpg`, `nagib_chauhan_kappa/B`, `pipe_friction_factor_turbulent`; **`rough_wall_log_law`, `friction_velocity_from_wind`, `drag_coefficient_neutral`, `dimensionless_shear(zeta, beta, unstable=)`, `surface_layer_wind(z, u_star, z0, L_M, *, kappa, unstable=, psi_at_z0=)`**; `LOG_LAW_CONSTANTS` (cited presets; κ, B, Π, z₀ are required keywords) | `wall_turbulence` (ch12) | (12.76)–(12.93), §12.11 log-linear profile | analytic, symbolic, benchmark (Lee & Moser DNS), approximate (Spalding, stress partition), unread source (Businger–Dyer) |
| **`gfd` (GFD)** — closed forms for a rotating, stratified thin layer (103 functions; f signed, `ValueError` at f = 0): `coriolis_parameter_deg`, `hemisphere`, `inertial_period`, `beta_parameter`, `beta_plane`, `coriolis_acceleration_local`, `thin_layer_terms`, `thin_layer_residuals(*, …)`, `eddy_friction`; `geostrophic_velocity`, `geostrophic_from_field/height`, `rossby_number(U, f, L)`, `ekman_number`, `parcel_adjust`; `thermal_wind_shear`, `thermal_wind_from_temperature(…, *, alpha)`; `ekman_depth(convention=)`, `ekman_surface`, `ekman_transport(_partial)`, `ekman_bottom`, `ekman_force_balance`, `ekman_finite_depth`, **`ekman_solve(z, K, f, *, …)`** (K(z)), `ekman_pumping*`, `sverdrup_transport`; `long_wave_speed`, `equivalent_depth`, `baroclinic_mode_speed`; **`shallow_water_omega(k, l, c, f0, beta)`** (three roots or NaN), `shallow_water_discriminant`, `dispersion_term_sizes`, `poincare_*`, `inertial_oscillation`, `kelvin_wave`, `kelvin_decay_side`; `rossby_radius*`, `geostrophic_adjustment_1d`, `adjustment_energy`; `inertia_gravity_*`, `wkb_vertical_structure`, `lee_wave_m`; `rossby_omega`, `rossby_group_velocity`, `rossby_omega_circle`, `stationary_rossby_wavelength`; `potential_vorticity`, `step_vorticity`, `rayleigh_kuo_criterion`; `eady_*` (`eady_max_growth_rate(f, N, dUdz)`: keyword calls); `fjortoft_transfer`, `two_d_cascade_spectrum`, `rhines_length`; constants `OMEGA_SOLAR_DAY`, `SIDEREAL_DAY`, `EARTH_RADIUS_MEAN` | `gfd` (ch13) | (13.5)–(13.145) | analytic, symbolic, converged, benchmark (Eady), analytic (ours: adjustment, −3 range, Rhines) |
| **`vertical_modes` (VM)** — Sturm–Liouville modes of N²(z): `vertical_modes(z, N2, g, n_modes, lid="free"\|"rigid", method="fd"\|"cheb")` → `Modes(c, psi, He, z, weights)`, `vertical_modes_shooting`, `modes_uniform_N`, `orthogonality_matrix(kind="psi"\|"energy")`, `project`, `reconstruct`, `w_structure`, `rho_structure`, `rigid_lid_error`, `modal_amplitudes`, `wkb_mode_speed` (ψ_n(0) = 1; rigid lid: index 0 is the first baroclinic mode) | `vertical_modes` (ch13) | (13.52)–(13.71) | symbolic, analytic, converged ("cheb" loses the barotropic speed) |
| **`shallow_water` (SW)** — numerical models, **every scheme ours**: `ShallowWater(nx, ny, Lx, Ly, H, f0, beta, g, bc, bottom)` (C-grid `[j, i]`, state dict `eta`/`u`/`v`), `make_state`, `gaussian_bump`, `geostrophic_state`, `kelvin_state(wall=)`, `continuity_tendency`, `momentum_tendencies`, `step` (SSP-RK3), `run`, `dt_limit`, `energy`, `volume`, `relative_vorticity`, `sw_potential_vorticity`, `interpolate`, `advect_particles`, `particle_potential_vorticity`; `linear_1d_step/run` (forward–backward); `qg_linear_evolve(_1d)` (exact spectral); **`barotropic_run(n, L, zeta0, *, t_end, dt, beta, nu, dealias)`** (pseudo-spectral RK4, strict 2/3 rule, start projected), `barotropic_vorticity_rhs`, `barotropic_velocity`, `barotropic_invariants`, `spectral_centroids`, `zonal_energy_fraction`, `barotropic_spectrum`, `random_vorticity` | `shallow_water` (ch13) | (13.44)–(13.45), (13.88)–(13.94), (13.117), (13.122), (13.143)–(13.144) | converged (space 1.98, time 2.98), conserved, analytic; turbulence caches qualitative |
| **`airfoils` (AF)** — sections from the Zhukhovsky map and their ideal flow (Γ clockwise via `Gamma_cw=`; dw/dz = u − iv; radians; section forces per unit span; contour counterclockwise from the trailing edge at index 0): `zhukhovsky_circle(b, e, beta, *, a)` (μ = b − a e^{−iβ}, a = b(1 + e)/cos β, ours), `zhukhovsky_circle_from_centre`, `zhukhovsky_airfoil`, `zhukhovsky_inverse` (root nearest μ), `section_geometry(x, y, te_index)`, `trailing_edge_angle`, `zhukhovsky_shape(exact=)`, `circular_arc`, `zero_lift_angle`, `camber_from_beta`; `kutta_circulation(U, a, alpha, beta)`, `kutta_circulation_root`, `stagnation_points`, `rear_stagnation_point`, `trailing_edge_speed`; **`zhukhovsky_flow(U, b, e, beta, alpha, *, Gamma_cw)` → `ZhukhovskyFlow`** (14 methods: `w`, `dwdz`, `velocity`, `speed`, `psi`, `phi`, `cp`, `zeta_of`, `inside`, `speed_near_edge`, …), `flat_plate_flow`; `zhukhovsky_surface`, `surface_circulation`, **`zhukhovsky_force(…, method="blasius"\|"pressure")`**, `zhukhovsky_lift_coefficient(alpha, beta, e, exact=)`, `normal_force_coefficient`, `lift_from_cp`, `flat_plate_forces(rho, U, c, alpha)`, `resolve_lift_drag`, `section_coefficient`, `circulation_from_cl`, `cl_from_circulation`, `thin_airfoil(dyc_dx, alpha, n)` (ours) | `airfoils` (ch14) | (14.3)–(14.12) | analytic, symbolic, converged, conserved, benchmark (cited forms) |
| **`lifting_line` (LL)** — finite wings (y from −s/2 port to +s/2 starboard, span angle γ from 0 to π; **w positive downward; K per radian, a required keyword**): `aspect_ratio(s, A)`, `mean_chord(A, s)` (keyword calls), `planform(kind, s, c_root, taper, n)`; `span_angle`, `span_position`, `circulation_series(Gn, gamma)`, `trailing_sheet_strength`, `elliptic_circulation`; `glauert_integral`, `downwash(y1, Gamma_fn, s, method)` (principal value), `downwash_series(Gn, gamma1, s, printed_middle=)`, `elliptic_downwash`; `induced_angle(small=)`, `effective_angle`, `strip_forces(exact=)`, `lift_and_induced_drag`, `lift_from_series`, `induced_drag_from_series`, `span_efficiency`, `section_lift_coefficient`; `lifting_line_modes/matrix/rhs/residual`, **`solve_lifting_line(c, s, U, alpha, *, K, beta, twist, N)`** (c a callable or a `planform` dict; collocation at γ_j = jπ/(N + 1); Cauchy order 3.8), `elliptic_induced_drag`, `induced_drag_coefficient(CL, AR, e_span)`, `finite_wing_lift_slope(K, AR)`, `elliptic_wing`, `twist_for_elliptic_loading`, `elliptic_twist_rectangular`; `horseshoe_velocity` (uz upward), `horseshoe_downwash`, `horseshoe_wing`, `ground_effect_factor(h, s)`, `downwash_along_x`, `staircase_from_loading(spacing="uniform"\|"cosine")`, `staircase_downwash` (order 1.5 / 2.0), `rolled_up_pair`, `wake_kinetic_energy_per_length` | `lifting_line` (ch14) | (14.1)–(14.2), (14.13)–(14.27) | analytic, symbolic, converged, conserved, benchmark (cited forms) |
| **`gasdynamics` (GD)** — one-dimensional and planar perfect-gas flow in M (82 names; γ an argument, R and C_p only via `PerfectGas`, **one gas per call**, station 1 upstream, radians): `PerfectGas(gamma, R)`, `AIR_BOOK`, `mach_from_speed`; `mass_flow`, `total_enthalpy`, `energy_residual`, `wall_force`, `nozzle_isentropic`; `stagnation_{temperature,pressure,density}_ratio`, `mach_from_pressure_ratio`, `isentropic_state`, `isentropic_table`, `sonic_ratios`, `critical_pressure_ratio`, `pitot_compressibility_factor`, `max_speed`; `area_ratio`, `mach_from_area_ratio(…, *, branch)`, `mass_flux_function`, `choked_mass_flow`; **`influence_coefficients(M, gamma, *, dA_over_A, df, dq, printed_15_32/47/48)`**, `area_velocity`, `duct_role`; `normal_shock(M1, *, allow_expansion)`, `normal_shock_entropy`, `weak_shock_entropy`, `shock_admissible` (M₁ ≥ 1), `normal_shock_stagnation_pressure_ratio`, `normal_shock_fluxes`, `hugoniot_pressure_ratio`, `strong_shock_limits`, `moving_shock`; `mass_flow_vs_back_pressure`, `convergent_nozzle`, `laval_critical_pressures`, `laval_regime`, `laval_shock_location`, `laval_nozzle`, `laval_exit_state`; `duct_parameter_A`, `duct_outlet_mach(M1, f, q, *, branch)`, `duct_state_ratios`, `duct_choking_limits`, `fanno_*`, `rayleigh_*` (ours); `mach_angle`, `oblique_shock(M1, *, delta, sigma, branch)`, `deflection_angle` (bare), `max_deflection`, `sonic_deflection`, `detachment_mach`, `weak_wave_pressure(M, delta)`, `prandtl_meyer`, `mach_from_prandtl_meyer`, `prandtl_meyer_turn(M1, turn)` (expansion positive); `riemann_exact` | `gasdynamics` (ch15) | (15.1)–(15.56) | analytic, symbolic, conserved, benchmark (NACA 1135) |
| **`acoustics` (AC)** — linear acoustics and simple waves (subscript 0 = undisturbed; f right-going, g left-going): `isentropic_pressure_perturbation(order=)`, `acoustic_density`, `linearisation_error`, `acoustic_impedance`, `sound_pressure_level`; `convected_wave_residual`, `wave_equation_step(p_now, p_prev, C, bc)` (leapfrog, ours; exact at C = 1); `dalembert(x, t, c, f, g)`, `plane_wave_velocity`, `velocity_from_pressure`, `dalembert_initial_value`, `dalembert_in_stream`, `wave_speeds`, `xt_characteristics`, `acoustic_energy`; `simple_wave_evolve(u0, x, t, c0, gamma)` (γ = 2 is ch07's), `breaking_time`; `source_wavefronts`, `in_zone_of_action`, `mach_cone_delay`, `doppler_factor` | `acoustics` (ch15) | (15.4)–(15.17), (15.49) | analytic, symbolic, converged, conserved |
| **`womersley` (WOM)** — oscillatory flow in a rigid round pipe (e^{+iωt}, G = −dp/dx, real part physical; ratios with `jve`, no overflow; scalar-callable): `womersley_number(a, omega, nu)`, `stokes_number`, `womersley_from_re_st` (diameter Re, St), `womersley_regime(alpha, *, low, high)`, `bessel_j0_series`, `bessel_j1_series`, `bessel_ratio_j0_j1`, `I12`, `I32`; `womersley_F1`, `womersley_F2`, `womersley_F2_series`, `womersley_F2_asymptote(*, printed_16_96)`, `womersley_profile_complex`, `womersley_velocity(r, t, a, omega, rho, mu, G)`, `womersley_velocity_nd`, `womersley_harmonics(…, G_n)`, `womersley_flow_rate(…, G)`, `womersley_flow_rate_harmonics(…, G_n)` (raises for a bare scalar), `womersley_wall_shear`, `womersley_XY`, `womersley_flow_ratio`, `womersley_phase_lag`, `womersley_small_alpha`, `womersley_large_alpha`, `womersley_ode_residual`, `womersley_pde_solve` (Crank–Nicolson, ours) | `womersley` (ch16) | (16.12), (16.69)–(16.100) | analytic (mpmath 1e-11), symbolic, converged (CN 2.00), benchmark (forms) |
| **`rheology` (RH)** — generalised-Newtonian laws and pipe fluxes (γ̇, τ magnitudes ≥ 0; SI with CGS converters): `newtonian_stress`, `apparent_viscosity` (secant), `casson_constants(k_c, s0, mu_p)`, `casson_stress(gamma_dot, tau_y, Kc)`, `casson_shear_rate`, `casson_apparent_viscosity`, `differential_viscosity` (tangent), `casson_stress_plasma_form`, `whitmore_casson`, `yield_stress_correlation(H, C_F, *, warn_range)` (H a fraction); `casson_wall_stress`, `casson_min_pressure_drop`, `casson_plug_radius`, `casson_pipe_velocity(r, a, dp, L, tau_y, Kc)`, `casson_core_velocity`, `casson_pipe_flow_rate`, `casson_poiseuille_flux_ratio(xi)`, `casson_flux_ratio`; `power_law_pipe_velocity/flow_rate(a, dp, L, m, n)`, `herschel_bulkley_pipe_flow_rate` (ours); `DYN_CM2`, `POISE` | `rheology` (ch16) | (16.5)–(16.7), (16.207)–(16.220) | analytic, symbolic, conserved |
`ST.rayleigh_eigs_contour` gained the keyword **`beta=0.0`** in ch13 (Rayleigh–Kuo; default path unchanged to the last bit).
`CONF.circle_flow_zeta` gained the keyword **`center=0j`** in ch14 (off-centre circle; the centred call is bit-identical).
Chapter-only code (move to `core/` on second use): **ch16** units (`MMHG`, `PRU`, `cP_to_Pa_s`, …), cardiac bookkeeping
and `pv_loop_work`, `vascular_resistance(mu, L, a)` (keyword calls), two-layer and FL closure (qualitative), tube law and
distensible tube, Windkessel (`windkessel_periodic_state` closed form), Moens–Korteweg, `tube1d_run` (MacCormack, ours),
elastic-tube quadratic, Murray, junctions + `network_pulse_run`, Dean series (`DEAN_C4_DERIVED`), collapsible tube /
Starling (`collapsible_tube_law` guarded to α ≤ 1), xylem and water potential; **ch15** `wedge_flow`, `expansion_fan`, `fan_staircase`,
`curved_wall_mach_lines`, `ackeret_pressure`, `ackeret_airfoil`, `flat_plate_supersonic(*, method=)`, `lift_slope_vs_mach`,
`becker_shock_*`, `nozzle_with_shock`, `rocket_nozzle` (warns without a gas), `blast_constant_crude` (qualitative),
`euler1d_run` (MacCormack + artificial viscosity, DEVIATION), `mach_regime`, the sympy engines; **ch14**
`impulsive_start_lumped_vortex` (a model), `wagner_function(s)` (a **lift** ratio), `wake_rollup`, `generic_*`,
`stall_regime` (qualitative), `drag_polar`, `best_glide`, `compressibility_check`, `apparent_wind`, `sail_forces`,
`oscillating_foil` (→ **Ch. 16**); **ch13** `flow_over_step`, `rayleigh_kuo_eigs` (β ≥ 0), `eady_numeric_eigs`,
`lapse_rate_table(dT_dz, *, Gamma_a)`; **ch12** `gradient_richardson_thermal(…, *, Gamma_a)`, `monin_obukhov_*`,
`k_epsilon_*`; **ch11** `os_leading_mode`, `benard_*`, `lorenz_*`; tables `illustrative_inputs`, `conventions_table`,
`book_slips` (ch15: a DataFrame with `kind`, `changes_physics`), `traps` (ch13–ch15); the rest in `knowledge/chNN.md` §3.
**Default g**: ch04, ch05, ch07, `core.waves` G_BOOK = 9.81; the rest G0 = 9.80665. **Default ρ**: ch06 1.2 vs 1000; ch07,
ch08 1000; ch09 BL 1.2; ch10, ch11 non-dimensional; ch12–ch14 arguments; ch15 none (ratios in M and γ); ch16 arguments
(`BLOOD` preset ρ 1060 kg/m³, μ 3.5 mPa s, μ_p 1.2 mPa s, cited).
**Three air-constant sets (ch15 trap T4)**: `CP_BOOK` 1004.5; `PerfectGas()` R = 287.058 → C_p 1004.70; `AIR_BOOK` R = 287.0.

### Explainer engine (`assets/viz_lib.js`)
`Viz.app` (tabs Walkthrough / Explore / Explain / Derivation / Equations / Code / Check; fit-to-window with density levels and
pagers; transport, modes, presets, status, terms, inspector, notes, selftest parity), `Plot`, `Viz.field`, `Viz.num` (RK4,
odeint, brentq, erf/erfc **~1e-7 only — a defect**, niceTicks), `Viz.work`, `Viz.three`, KaTeX with fallback; `Viz.font`,
`Viz.roundRect`, `Viz.text`, `Viz.card`, `Viz.fmtTime`, `Viz.tnum(keepTiny)`; pager slicing, `ev.viewId`/`ev.vizView`.
Quirks: Code-tab `{{placeholders}}` share one namespace; `#param=` deep links are overridden by tour step 1; `Viz.tnum(x) +
'^2'` is a KaTeX double superscript (bracket with a local `tsq`/`tpow`); `hidePortrait` keeps the view's flex share; view
titles are not re-evaluated on a tab change; no `.c-amber` class; no log axes, no linear algebra, no complex numbers;
**a preset chip is pressed when every key of its `set` matches the state** — a `set` holding only a shared mode key lights
every chip (ch15 Must; a `get set()` getter is evaluated on each refresh). **Library fixes during ch14** (commits `8798d2a`,
`dbde04f`): L1 pages settle after a state change with a tab open; L2 `Pager.setItems(list, keepPage)` keeps the page only
when asked; L3 one handler on "Next →". **No JS helper promoted after ch02–ch16** (the site-publisher was reading `viz/`
each time). Ranked candidates in `knowledge/viz_patterns.md` (ch16 top five: TeX-safety checks in shot.py/lint,
`Viz.num.cbessel`, `Viz.fitText`, `Viz.mix`/`ramp`, `Viz.arrowPx` + `v.size`; ch15: preset-exclusivity probe, `Viz.gas`;
ch14: shot.py click and Explain-state checks, `Viz.cx`, `transport.max`/`holdAt`, a y-title strip). **JS strings: a lone
backslash is silently eaten (`\f` form feed, `\r` carriage return) — `String.raw` for TeX.** Reference:
`templates/viz_example.html`.

## Explainer inventory
| chapter | slug | CORE idea | reusable stage pattern |
|---|---|---|---|
| ch01 | `continuum_averaging_volume` · `viscosity_momentum_diffusion` · `heat_work_paths` · `parcel_stability` · `buckingham_pi_machine` | C06 · C12 · C25–C45 · C50–C55 · C64–C69 | molecules vs box; gap + steady ghost; piston · p–v · T–s; column + two-convention badge; chips → matrix → groups |
| ch02 | `rotation_of_axes` · `cauchy_traction_principal_axes` · `strain_vs_rotation_split` · `gauss_flux_box` · `stokes_circulation_loop` | C02–C16 | clickable C matrix; element + Mohr; three squares; box + waterfall flux bars; unrolled u·t |
| ch03 | `flow_lines_unsteady` · `material_derivative_probe` · `galilean_frames_cylinder` · `fluid_element_deformation` · `spin_and_principal_axes` · `vortex_paddle_wheels` · `reynolds_transport_cv` | C01–C15 | streamline/trail/dye on one clock; local/advective bars; observer slider; G deforms a square; paddle wheels; moving CV + waterfall |
| ch04 | `control_volume_budgets` · `stream_function_spacing` · `newtonian_stress_lab` · `navier_stokes_term_balance` · `rotating_frame_coriolis` · `which_bernoulli` · `viscous_dissipation_heating` · `boussinesq_buoyancy` · `dynamic_similarity_models` | C01–C15 | five scenes, one budget; ψ spacing = speed; G → S → τ; term bars summing to 0; two observers; hypothesis table; kept vs dropped bars |
| ch05 | `vortex_tubes_cannot_end` · `vortex_pressure_funnel` · `kelvin_material_loop` · `baroclinic_torque` · `vorticity_stretching_tilting` · `biot_savart_filament` · `vorticity_equation_rotating` · `point_vortex_lab` · `vortex_sheet_rollup` | C01–C14 | 3-D tube + Gauss bars; torque ◇; "build the sum"; conserved quantity flat in amber; click-to-place vortices; roll-up with KH ghost |
| ch06 | `superposition_sandbox` · `cylinder_circulation_lift` · `vortex_wall_images` · `complex_potential_corners` · `blasius_kutta_contour` · `conformal_joukowski` · `laplace_relaxation` · `axial_singularity_bodies` · `added_mass_sphere` | C03–C15 | element kits; ◇ ρUΓ bars; sensor trace; Laurent bars; two planes, same particles; stencil stepping + honest order |
| ch07 | `dispersion_relation` · `particle_orbits` · `capillary_gravity_waves` · `seiche_standing_waves` · `group_velocity_packets` · `wave_rays_refraction` · `hydraulic_jump` · `two_layer_modes` · `internal_wave_beams` | C03–C16 | tank + c(λ) with regime bands; linear vs exact paths; chord vs tangent; in-page rays; budget bars + forbidden preset; square K plane |
| ch08 | `couette_poiseuille_backflow` · `lubrication_scaling` · `slider_bearing` · `viscous_gravity_current` · `stokes_first_problem` · `similarity_exponents` · `oscillating_plate` · `stokes_sphere_flow` · `stokes_drag_settling` | C02–C15 | line + parabola = sum; log bars with a "dropped" band; printed ghost with a measured miss; log clock + raw/rescaled; exponent plane; rejected root drawn; traction sweep |
| ch09 | `bl_scaling_thicknesses` · `blasius_similarity_collapse` · `falkner_skan_family` · `thwaites_marching` · `cylinder_drag_crisis` · `karman_street_stability` · `free_jet_similarity` · `wall_jet_invariant` · `teacup_secondary_flow` | C01–C14 | streamlines lifted by δ*; raw ↔ rescaled + shooting overlay; one dial through a fold; outer-flow picker → criterion, book first; log-Re dial; configuration + growth + spectrum (Ch. 11 template); invariant bars (one falls, one flat); force bars adding with height (Ekman seed) |
| ch10 | `fd_stencil_order` · `von_neumann_amplification` · `upwind_cfl_advection` · `cell_peclet_wiggles` · `fem_hat_assembly` · `mac_projection_staggered` · `lid_driven_cavity` · `mixed_fe_lbb` | C01, C03 (D01, D03) · C04, C02 (D04–D07) · C05, C10 (D08, D18 ★★★) · C09 (D15–D17) · C06–C08 (D09, D11, D12, D14) · C11, C12 (D19–D21) · C14, C15 (D23) · C13 (D22) | log–log error vs h with a round-off band and best-h ◆; complex-plane G(θ) + parameter region with a brute-force scan + live kick; x–t stencil with the characteristic's foot and interpolation weights; discrete-root plot with an r < 0 band and a sign strip of rʲ; assembly as a transport (element block lands, outlined); algorithm stages on a transport (u* → Poisson → push → uⁿ⁺¹); live solver (banded LU, background pump); unknown counting → argument |
| ch11 | `normal_mode_growth` · `kelvin_helmholtz_boundary` · `benard_neutral_curve` · `salt_fingers` · `taylor_couette_onset` · `richardson_shear_instability` · `inviscid_shear_criteria` · `orr_sommerfeld_neutral_curve` · `lorenz_attractor` | C01, C05 (D01, D12) · C02 (D02–D04) · C03, C04 (D09–D11) · C06 (D13) · C07 (D14, D15) · C08, C09 (D17, D18) · C12, C10 (D22, D19) · C11, C13, C14 (D21, D23) · C15 (D24, D25) | three systems in one state: the wave, σ(k) with the growing band, the σ-plane; term bars adding to a discriminant + c-plane where two roots collide; (K, Ra) plane + **table of fluidpy growth rates with exact zero on the live neutral curve**; parcel + regime map with a verdict text pinned by checksum rows; ring swap with energy bars, a slider whose displayed value depends on another parameter; profile on a log axis with the criterion band + growth map + three-way status; necessary conditions as bars on a computed mode, interpolation ending on the exact neutral point; click a point on a neutral-curve map → mode + production/dissipation bars (Ch. 12/13 template); two runs on one clock, canvas 3-D orbit, status computed from the run |
| ch12 | `reynolds_averaging_window` · `correlation_and_spectrum` · `reynolds_stress_parcels` · `energy_cascade_spectrum` · `turbulent_energy_budget` · `turbulent_jet_similarity` · `law_of_the_wall` · `mixing_length_closure` · `stratified_surface_layer` · `taylor_dispersion` | C01, C02 (D01, D03) · C02, C03 (D03–D05) · C04 (D06) · C07, C08 (D11, D12) · C06 (D09, D10 ★★★) · C09 (D13–D16) · C10, C11 (D17–D20) · C12 (D21) · C14, C15 (D24, D25) · C16 (D26–D28) | two estimators of one mean on one clock with an error split (noise · drift · leak, ◆ formula beside ◇ measured); drag a lag to build r(τ), area = scale, linked spectrum; click one sample → its quadrant → thousands with a ± 5 s.e. band; ladder on a log axis + model spectrum + real-case table + units inspector; **paired budgets sharing one mirrored term, bars summing to zero** (Ch. 13 balances); raw ↔ rescaled with an invariant line that turns rose; **draggable fit window on DNS data**; "κ turns it, damping slides it" with the model's failure shown first; **two-line two-convention badge, a flux as the transport, a formula ending with ✕ at the edge of its range** (Ch. 13 template); log-paced transport, puff ↔ plume modes, constant-D ghost |
| ch13 | `geostrophic_balance` · `thermal_wind` · `ekman_spiral` · `ekman_force_balance` · `vertical_modes` · `shallow_water_dispersion` · `kelvin_wave` · `geostrophic_adjustment` · `rossby_waves` · `eady_instability` | C02 (D03) · C03 (D04) · C04, C05 (D06 ★★★, D07) · C06 (D08) · C08 (D11 ★★★) · C09, C10 (D12 ★★★, D13) · C11 (D15) · C12 (D16 ★★★) · C15, C13 (D22 ★★★, D23) · C16 (D25, D26 ★★★) | parcel on a pressure map with force arrows closing tip to tail and per-axis bars that sum to zero; **sign bars that swap sides while the result bar stays**; **one cursor shared by a 3-D curve, its hodograph and a running integral**; force triangle leaning toward low pressure in either hemisphere; shape · ladder · travelling mode on one mode number with a live overlap; **log term bars "smallest / largest" for the selected root**, branch diagram with a shaded no-real-roots band; **equations as pairs of equal bars with the refused solution drawn**; **snapshot vs period mean vs prediction, conserved quantity flat in amber** (Ch. 14 starting vortex); **chord and tangent on ω(k) + wavenumber circle + a clicked column's arithmetic**; a threshold crossed by one cursor with ◆ / ▲ presets; hemisphere chips wherever a direction word appears |
| ch14 | `zhukhovsky_airfoil_builder` · `kutta_condition` · `starting_vortex` · `airfoil_pressure_lift` · `lift_curve_stall` · `downwash_induced_drag` · `lifting_line_solver` · `elliptic_loading_aspect_ratio` · `thrust_from_lift` | C04 (D03 ★★★, D04 ★★★) · C02, C05, C06 (D07, D09) · C03 (D02) · C07 (D10) · C06 (model, no derivation) · C08–C10 (D11–D13) · C11–C13 (D15 ★★★, D18 ★★★) · C14, C15 (D19, D21) · C16 (D22) | a draggable parameter point with guides naming its two directions + rule-against-measured bars; **a hypothesis run as an experiment, a quantity reading "unbounded" until the right value**; **conserved total flat in amber with a cut, a cited curve only on the axis of the quantity it is a ratio of**; three pictures of one flow with an area that is a force; a model labelled "model" everywhere; **"build the sum" 1 → 3 → 7 → 21 → ∞ with signed term bars** on a canvas line projection; **the matrix row with the reader's numbers** + log bars with shares + "make it perfect"; **a U-shaped curve from one falling and one rising term with a wall, a ghost and a computed badge** (Ch. 15 area–Mach); **a force polygon that closes, modes = one triangle in three systems** (Ch. 16) |
| ch15 | `sound_pulse_dalembert` · `area_mach_nozzle` · `normal_shock_jump` · `nozzle_back_pressure` · `duct_friction_heating` · `mach_cone` · `oblique_shock_wedge` · `prandtl_meyer_fan` · `supersonic_flat_plate` | C01, C02 (D03–D05) · C05, C06 (D10, D12) · C07, C08 (D13 ★★★, D17 ★★★) · C09, C10 (D19, D20) · C11, C12 (D21 ★★★, D23) · C13 (D25) · C14 (D26, D27) · C15 (D29, D30 ★★★) · C16 (D28, D31) | two shapes sliding apart over an air strip "drawn ×500, true ≤ 0.886 mm"; **term bars adding to zero through a removable 0/0**; **"1 = 1" conservation bars beside a signed entropy bar that turns red** (admissibility); a transport walking a shock down a nozzle, the other nozzle as a ghost; corrected slips live beside the printed ones; **an arrival-time view whose minimum is the boom**; velocity triangles at a clicked streamline, two roots merging; a staircase of waves with an error-halving invariant; an area between two levels is a force, ch14's loop as a ghost |
| ch16 | `casson_blood_flow` · `windkessel_heartbeat` · `pulse_wave_artery` · `womersley_profiles` · `womersley_flow_phase` · `murray_branching` · `junction_reflection` · `dean_curved_pipe` · `collapsible_tube_choking` · `xylem_water_potential` | C02, C16 (D02, D31–D33 ★★★) · C07, C01 (D01, D09, D10) · C08, C06 (D08, D11, D12) · C09, C03 (D03, D13 ★★★, D14, D15) · C10 (D16, D17) · C12 (D22, D23, D24 ★★★) · C13 (D25, D26) · C14 (D27 ★★★) · C15 (D28, D29, D30 ★★★) · C17 (D34) | a stress panel sharing the r axis so τ = τ_y is the plug edge; beat-by-beat closed form with storage + run-off bars; one exaggeration factor with true sizes; full-diameter profiles with both limits as ghosts and a Hankel-scaled J₀; **time traces + Bode pair + complex locus on one α (the APS pattern)**; a draggable point whose force triangle closes at the optimum, a printed slip as a mode; exact closed-form pulses with a grid-echo row and energy bars to 1; convention chips κ/K/D on one slider, a corrected slip with both numbers; **a model's validity range computed in the status line, news arcs that stop at S = 1**; a ψ ledger of term bars at stations |

## Notation
See `knowledge/notation.md` (sign traps first, then the register, then per-chapter conventions). Everywhere: SI and kelvin
inside functions; **z up** (except ch06 §6.8); **lapse rate Kundu Γ ≡ dT/dz in code, meteorology shown**; `p0` ≠ `p_ref`;
pressures absolute unless `_gauge` (ch07 gauge); signed τ_xy = μ ∂u/∂y. ch02: passive C; traction on the first index; grid
`[k, j, i]`. ch03: R = G − Gᵀ, spin ½ω. ch04: Coriolis acceleration +2Ω × u′ vs force −2Ω × u′; **2-D ψ: u = ∂ψ/∂y** (ch07:
u = ∂ψ/∂z). ch05: γ = u_below − u_above. ch06: **Γ counterclockwise in code, `Gamma_cw=`/`Gamma_ccw=`**. ch07: **ω = angular
frequency**, η = elevation, g′ = g(ρ₂ − ρ₁)/ρ₂. ch08: **dp/dx book sign, `G=` = −dp/dx**; η = y/√(νt). ch09: x along the
wall; dp/dx > 0 adverse; θ = momentum thickness; FS exponent book n = code `m`. ch10: three grids; θ = Fourier angle; C = 2α
Courant; Ma for Mach. ch11: **σ = −ikc is the growth rate**; K Bénard wavenumber; Γ three ways; Re per flow. ch12: over-bar
= **ensemble** average; **lower-case u, v, w, T′ are fluctuations**; **κ = von Kármán**; ē turbulent kinetic energy; **η =
Kolmogorov length**; **h = full channel height**; heat flux positive upward, L_M > 0 stable; `uv_plus` keys hold MINUS the
correlation; **spectra two-sided, ∫S = variance**; `gradient_richardson_thermal` REQUIRES `Gamma_a`.
**ch13:** x east, y north, z up; **f signed** (abs for scales, sign for directions, `ValueError` at f = 0); which z = 0
changes by section; **ψ with u = −∂ψ/∂y**; Ro = U/(abs(f)L); three Rossby radii; V = u + iv; §13.18 spectra one-sided.
**ch14:** **Γ clockwise** (`Gamma_cw=`); dw/dz = u − iv; **w positive downward**; Λ = aspect ratio; `e_span`; K per radian.
**ch15:** station 1 upstream; **subscript 0 undisturbed in §15.2, stagnation after**; σ, δ (β, θ elsewhere), turn = −δ_book;
μ = Mach angle, ν = Prandtl–Meyer (not viscosities); f = F/(p₁A), q = Q/h₁; three air-constant sets, one gas per call.
**ch16:** **G = −dp/dx, τ = −μ du/dr ≥ 0 (ch08 signed)**; **two time conventions — e^{+iωt} rigid tube (real part, i^{3/2}),
e^{i(kx − ωt)} elastic tube** (T1); **transmural p − p_e**, a₀ at p = p_e; **Shapiro's P < 0 in collapse, (16.188) only for
α ≤ 1**; S = u/c, M = 2 − n (not Mach); α, K, R, δ, X overloaded → distinct code names; Dean κ, K = 2κ², D = 4√2κ;
hematocrit a fraction; α on the radius, Re and St on the diameter; slips S1–S34, traps T1–T27 (all chapters: notation.md).

## Teaching lessons
- **Depth tiered, coverage exhaustive** (A / B / C, derivations, cells): ch01 15/70/21, 12, 496 · ch02 16/60/18, 15, 405
  · ch03 15/52/12, 24, 413 · ch04 15/151/19, 30, 611 · ch05 14/69/11, 23, 496 · ch06 15/116/29, 31, 564 · ch07 16/193/27,
  37 (321 steps), 617 · ch08 15/113/15, 33 (301), 501 · ch09 14/118/16, 22 (224), 472 · ch10 15/105/23, 23 (222), 602 ·
  ch11 15/131/24, 25 (258), 618 · ch12 16/210/17, 28 (266), 690 · ch13 17/154/47, 29 (253), 821 · ch14 16/77/20, 22
  (173 steps), 549 · ch15 16/121/14, 31 (243 steps), 606 · **ch16 17/107/21, 34 (≈ 266 steps), 577**.
- **A spine sentence, pointed back to in every block** (ch11 "one recipe, seven times"; ch12 "an average fails on exactly
  one thing — a product"; ch13 "one set of equations with terms switched off").
- **One-way theorems get a counter-example and the word "necessary" or "sufficient"** (ch11; ch13 Rayleigh–Kuo); "theorem
  vs observation" for look-alike criteria (ch12); a ★★★ check cell rebuilds the derivation's own steps.
- **Convention callouts** ("which p_o?", Coriolis term vs force, Γ_cw vs Γ_ccw, which ω/η/g′/dp/dx/δ/Re/θ/λ/κ, **which
  hemisphere, which z = 0, which Rossby radius**): state each convention, the size of the difference, which the code uses;
  a conventions table up front; a library string in another notation is translated at every appearance. **Both lapse-rate
  conventions on screen, always. Every direction sentence names the hemisphere and gives the mirror rule** (ch13).
- **Book slips are taught by computing both versions**: a slips table up front, printed forms as code options a test must
  fail, a slip box that never over-claims; **"loose" slips are inconsistencies, not errors; an unruled slip keeps its
  caveat** (ch13 #14). **"Ours — not in the book" opens every extension**; "Our choice" boxes mark our schemes (ch13).
- **Make an approximation or an artefact measurable**; test an explanation by removing the alleged cause (ch12); a fitted
  constant is quoted with its window (ch12); **"honest caveat" boxes with a criterion and computed cases**; a ringing run
  is compared through its period mean, not a snapshot (ch13).
- **Every key number by two independent routes; every prose number printed by a cell**; captions count what they claim;
  **a figure note is written from the rendered figure** (ch13: six Must items came from notes written from the storyboard).
- **A tool first used inside a derivation gets its primer first**; boundary conditions are signposted just before the
  derivation that needs them. `knowledge/primers.md` lists 396 (P01–P395 + P218a); **the next free number is P396**.
- **Builder hygiene**: design TeX converted, not copied; raw strings; never a heredoc; self-checks tested on planted
  damage. **Statistics are estimates** (ch12): seeded cells, 5 standard errors, synthetic fields captioned "kinematic".
- **Climate hooks with numbers** (ch04 Ro; ch07 swell, tsunami; ch10 CFL on an ocean grid; ch12 30-year normals, L_M; ch13
  weather maps, upwelling, a 43 km Rossby radius, storms in 1.6 days, zonal jets; ch15 the standard atmosphere at 10 km; ch16 the soil–plant–air water path in ψ).
- Reusable derivation moves (ch01–ch16, 270 rows; 16 from ch16) are tabulated in `knowledge/concept_map.md`.
- **ch14**: a letter table "here / elsewhere"; a "first look" primer (P351); a cited curve on the axis of its ratio. **ch15**:
  a bridge between the book's parameter and the field's, slice by slice (f against 4f_FL/D); an inequality by factoring a
  ratio; sliders open on the case the text reads; every ratio says what it is a ratio of. **ch16**: **a quantitative word in
  a reading note gets a printing cell above it**; a definitions cell for symbols the book introduces silently (p_e,
  transmural, a₀); a numerical demo separates physics from artefacts and starts from an exact state (u + 4c uniform).

## Global pitfalls confirmed in this project
**Machinery and environment**
- Extracted text garbles maths (`¼` = `=`, `ð…Þ` = parentheses, missing minus, ω → `u`, γ → g, θ → q, σ/τ → s, ν → n,
  ε → 3, ρ → r, Ω → U, ψ → j, φ → 4, η → h, ζ → z, ∂ → v, ∇ → V): read page images.
- **Windows: heredocs and `sed` mangle backslashes, quotes and `|`** — write code, explainers, design and knowledge files
  with Write/Edit or Python files; builder strings with LaTeX must be raw; single-backslash TeX in JS strings becomes
  control characters (ch04, ch06, ch09–ch12). **A scratch script named like a standard-library module (`numbers.py`)
  breaks `import numpy`** (ch13). `MPLBACKEND=Agg` strips inline figures; Anaconda's kernelspec can shadow `fluidpy-venv`.
- **`tools/shot.py` blind spots**: no click or drag, **no `.katex-error` or control-character check, Explain audited only in
  the default and step states, no 844×345, no preset-exclusivity check (ch15)**; `py:` rows have restricted builtins
  (**plain lambdas cannot see the chapter module**); a table that rounds its own abscissa needs rtol 1e-4.
- **Machine load**: intermittent pager clips (a walkthrough `eq` card with live rows, ch13); run the merge gate serially;
  notebook runtime is load-dependent — measure on an idle machine.
- **Parallel agents need private scratch subfolders and unique file names** (collisions in ch10–ch13; a stray interpreter
  in ch13); do not edit a solver module while a shot run evaluates parity rows; after a crash, audit code against Part C.
- `coverage_check` matches ledger reminders by exact primer name; `eq_refs` misses numbers inside template strings. Book
  wording can hide in docstrings, exercise answers in computed numbers, table constants in analysis, design, scripts and
  reports, **and book values in words in derivation steps** (ch13 D13, D27).
- **Design and analysis documents contain errors; downstream agents compute, not copy** (every chapter; **ch13: a boundary
  term in the wrong orthogonality relation, the height of the Ekman overshoot, the WKB error order, the sign of the
  cubic's slow root, "settles" for a frictionless parcel**).

**Mathematics → code**
- **Defaults encode physics** (ch09 `cp_base=None`; ch10 a default Δt that dropped the viscous rate; **ch12 `Gamma_a=0.0`
  called the standard atmosphere unstable**): convention-carrying and empirical arguments are required keywords; explicit
  solvers compute the full stability limit and raise.
- **Fitted constants belong to their window** (ch12: κ = 0.400 from y⁺ = 30, 0.384 from 350 on the same DNS). **Names
  carry signs and normalisations** (`uv_plus` holds −mean(uv); `C1` one-sided; ch15 `stagnation_pressure` = p + ½ρU²).
  **Estimators have artefacts of predictable size** (the biased correlation's k⁻² tail); seed, 5 standard errors.
- **Classifiers and verdict strings need tests on both sides of every boundary and at degenerate inputs** (ch11: 128
  tests passed with two classifier bugs; ch12: boundaries compared on the ratio and mirrored identically in JS; a regime
  name is not a stability verdict; NaN outside a formula's range). **A pinned number can pin a defect**: pin sets and physics.
- **Unbounded profiles (ch11)**: scale the box with 1/k, prove it by box and N doubling; **a convergence filter can make a
  false zero or a lower mode** — leading mode or raise; a complex path near a neutral point (`knowledge/ch11.md` §4).
- **A test comparing two of our own functions is not evidence**, nor is a table written by the same function; **prove
  discrimination by planting the wrong variant** (ch02 8/8 … ch10 38/38, ch11 60/61, **ch12 42/42 + 18/18**); an identity
  that holds for any input cannot test which input was used (ch12 R16); mutants must reach `@slow` tests too; heavy cached
  runs are invisible to mutants — keep one slow test that recomputes a coarse cached run.
- **A test on |x| cannot see the sign of x** (ch09): signed tests on manufactured fields, both signs of every parameter.
  **Series stop on the term envelope, never on the actual term** (ch10). **Build operators link by link** and test null
  vectors on 1- and 2-cell grids; guard every division's degenerate case; cap solver iterations.
- **Folds and bifurcations are mathematics** (ch09): parametrise by the quantity that passes through. **BVPs on infinite
  domains**: scaling symmetry or `solve_bvp` with continuation. **Complex-step derivatives need analytic continuations**.
- **Code the corrected form; keep the printed one as a named option a test must fail** (every chapter since ch07).
  **Never type a constant a function can compute** (ch01, ch09; ch12 explainers: a typed 2.240 went stale).
- **Overflow and cancellation**: cosh/sinh ratios; expm1; erfc; scaled exponentials beyond R ≈ 700; 0/0 starts with a
  finite limit. **Stiff and degenerate PDEs**: CN after a jump needs BE start-up; a √ψ wall singularity needs a mapping.
- **Observed orders only in the asymptotic range**; singular corners, staircase boundaries and heuristic closures set the
  global order (ch10 block p = 0.43); finite-difference tolerances from a written budget, never loosened by feel.
- **A telescoping sum is not a conservation test**; genuine V4 invariants: ch07 ∫η², ch09 J, Ψ, ch10 ∇·u, ch11 a mode's
  energy budget, ch12 Parseval, ch14 Kelvin's sum, **ch15 jump fluxes, Laval mass flow, acoustic energy, Riemann invariant**.
- **Sibling functions share defaults and argument meanings** or are called by keyword; convention-named keywords
  (`Gamma_cw=`, `ref=`, `G=`, `closure=`, `printed=`, `convention=`, `scheme=`, **`Gamma_a=`, `upto=`, `two_sided=`,
  `unstable=`, `window=`**).
- **Benchmarks from the primary source at its printed precision, read first-hand** (ch10 Ghia by two transcriptions;
  **ch12: an unread "Businger–Dyer" citation was relabelled V1 and flagged**). **Validation labels**: V4 only for real
  conservation laws; "converged" only with an asserted order; "benchmark" only inside the stated tolerance; **"constant
  consistency" for a constant compared with its citation (not V5)**; "approximate" for models against data with the band
  stated; "secondary" for compilations; "qualitative" for confined/heuristic studies.
- **ch13**: either sign of f (abs for scales, sign for directions), raise at f = 0; the 2/3 mask on the state too; **a new
  parameter is its own keyword**; **`numpy.allclose` has an absolute tolerance of 1e-8, larger than β** — relative
  tolerances only; measure an order before writing it; a box-truncated eigenproblem fails silently on a radiating mode.
- **Which quantity is the curve a ratio of? (ch14)** Wagner's lift ratio (½ at the start) sat on a circulation axis (0 by
  Kelvin) in five documents and **no mutant could see it**. A time-to-threshold carries its time step; split quadratures
  at kinks; test each member of a printed chain; a library fix can cause the next defect (L1 → L2).
- **ch15**: one gas per call (γ ≠ 1.4 tests); form the driver before dividing by 1 − M² (0·∞); guard bare formulas in their
  callers with a round-off allowance; a series near a branch point; a re-exported name that means something else gets a
  warning wrapper; measure an approximation's order; presets light together when `set` holds a shared key; port fixes.
- **ch16**: **a fitted law carries its validity range** ((16.188) was used at α = 1.5; now warn / NaN above α = 1, demo
  bound α_t ≤ 0.39685 computed); **closed forms beat bracketed root searches** (a fixed `brentq` bracket raised for wide
  sections); one scalar, one meaning across siblings; a PDE sign slip can hide until fourth order (S34, found by an
  independent Navier–Stokes solve); **lone backslashes in JS strings ate TeX and passed lint and shot** (`String.raw`);
  slug-prefixed scratch files; a builder's `--help` rebuilt the notebook (argparse); retry Windows `OSError [Errno 22]`.

## Benchmark inventory
| value / table | source (citation) | stored in | used by |
|---|---|---|---|
| k_B, N_A, R (exact); e; USSA-1976 (sea-level c = 340.29 m/s), Sutherland; IAPWS σ(T), μ(T); Jennings mean free path; Taylor blast 0.856; AMS Γ_d; EOS-80 | CODATA 2018/2022, NIST CUU; NASA-TM-X-74335, PDAS; IAPWS R1-76(2014), R12-08; Tsalikis et al. 2024; arXiv:2009.05674; AMS Glossary; Fofonoff & Millard 1983 | `reference/ch01/`, `reference/ch08/` | ch01, ch04, ch07, ch08, ch15 (c, mean free path) |
| Divergence theorem, Levi-Civita, Mohr, rotations, Stokes; Lamb–Oseen peak α = 1.256; sphere drag C_D(Re); WGS-84, orifice C_c 0.611, capillary length | Wikipedia; scipy docs; arXiv:2210.05223; F. A. Morrison (2016) | `reference/ch02/`–`ch04/` | ch02–ch05, ch08, ch09 |
| Sphere added mass ½; Rayleigh collapse 0.91468; source panels; Fenton & McKee, Guo, capillary minimum, Stokes drift H/λ = 0.1410633; slider K_opt 2.1889; viscous gravity current η_N 1.411 | Lamb (1932); Rayleigh 1917; Hess & Smith 1967; Coastal Eng. 1990, 2002; arXiv:1507.02784; San Andrés (Texas A&M); Huppert, JFM 121 (1982) | `reference/ch06/`–`ch08/` | ch06–ch08 |
| Blasius f″(0) = 0.332057336215196 (ours 5.7e-15); Falkner–Skan κ, separation β = −0.198837735; Hiemenz 1.232588, δ*/δ 0.6479; Bickley jet 0.4543, 0.2752, 3.3019, 0.5503; Kármán b/a ≈ 0.281 | Wikipedia; Belden et al., arXiv:1907.09912; Weidman & Turner; Horváth et al. 2020 | `reference/ch09/` | ch09 C04, C05, C10, C12 |
| **Lid-driven cavity centreline u(y) at x = ½, Re = 100, 400 (17 points)** (MAC 64² 0.26 % of U; MacCormack 128² 0.44 %) | Ghia, Ghia & Shin, J. Comput. Phys. 48, 387 (1982), Table I; transcriptions ivan-pi gist + cfdgasman/lid-driven-cavity; CMC 36 (2013) Table 1 cross-check | `reference/ch10/ghia1982_table1.csv` (`make_refs.py` asserts 5 values) | ch10 C14 (V5), E7 |
| **Primary-vortex centres** Ghia (0.6172, 0.7344) Re 100; Hou et al. lattice Boltzmann Re 100, 400 | Ghia via Hajabdollahi & Premnath arXiv:1202.6351; Hou et al., J. Comput. Phys. 118, 329 (1995), arXiv:comp-gas/9401003 | `reference/ch10/cavity_vortex_centres.csv` | ch10 C14 (V5) |
| Kovasznay flow (form; pressure ours by sympy); Radon 7-point rule (corroboration) | Kovasznay 1948 / Wikipedia; JAX-BEM issue #1 | `reference/ch10/SOURCES.md` | ch10 C13 (V1/V3) |
| Unbounded cylinder Re = 100 St 0.164–0.165 (comparison only; ours is confined) | compilation in arXiv:2303.09262 (Williamson 1996 et al.) | `reference/ch10/SOURCES.md` | ch10 R13 (qualitative) |
| **Plane Poiseuille Re_c 5772.22, α_c 1.02056, c_r 0.26400; eigenvalue 0.23752649 + 0.00373967i at Re = 10⁴, α = 1** (ours 4e-8; 1e-9) | Orszag, JFM 50, 689 (1971) | `reference/ch11/benchmarks.json`, `SOURCES.md` | ch11 C11, C13 (V5), E8 |
| **Bénard Ra_c**: rigid–rigid 1707.762 at 3.117, rigid–free 1100.65 at 2.682, free–free 657.511 at 2.2214; odd mode 17610.39 at 5.365 (source page not reachable; two routes) | Chandrasekhar (1961) via Nek5000 examples and arXiv:nlin/0302057 | `reference/ch11/benchmarks.json` | ch11 C04, C05 (V5), E3 |
| Blasius (parallel) Re_c 519.2, k 0.303, ω 0.120 (ours 519.0601, 2.7e-4); Bickley jet Re_c ≈ 4.0 at k ≈ 0.2 (ours 4.0170 at 0.1728; approximate source); tanh layer k 0.4446, kc_i 0.1897 (digits secondary) | Thomas via Gallagher, Griffiths & Stephen, Phys. Fluids 28, 074107 (2016), Table II; Tatsumi & Kakutani (1958) via a 2017 abstract; Michalke (1964) | `reference/ch11/benchmarks.json`, `critical_points.json` | ch11 C13, C12, C08 |
| Taylor narrow gap, co-rotation asymptote 1707.76[1 − 0.00761((1 − μ)/(1 + μ))²] (tertiary; ours 1.2e-6 … 1.1e-4); Lorenz r_H 24.74; Feigenbaum δ, a_n | Wikipedia "Taylor–Couette flow" (rev. 1373953105), "Lorenz system", "Feigenbaum constants" | `reference/ch11/benchmarks.json` | ch11 C07, C15 |
| Our ch11 data (public, labelled ours): Bénard neutral curves, Taylor table, TG growth map (620 points), Rayleigh spectra (8 profiles), OS grids and neutral curves (Poiseuille, Blasius, tanh, Bickley two bands), critical points, explainer tables | our runs (`reference/ch11/make_refs.py`, `scripts/ch11_tables.py`, `tools/viz_tables/ch11/`) | `reference/ch11/*.csv`, `*.json` | ch11 C04–C15, E3, E5–E8 |
| Our data (public, labelled ours): MAC cavity 16–128² (Re 100, 400), MacCormack 32–128², block C_D/C_L histories, FE cylinder steady and Re = 100 forces, inf–sup table, explainer tables | our runs (`scripts/ch10_*.py`, `reference/ch10/verify_runs.py`) | `reference/ch10/*.csv`, `explainer_tables.json` | ch10 C13–C15, E6–E8 |
| Forms (V1 cross-checks, not V5): Rankine, Lamb–Oseen, RTT; Bélanger, Taylor–Green; Kelvin ring, Burgers, Hill; cylinder, Kutta–Joukowski; KdV; Oseen, Stokes, Hagen–Poiseuille, Taylor–Couette | Wikipedia pages; K. T. McDonald | `reference/ch03/`–`reference/ch08/` | ch03–ch08 |
| **Channel DNS mean profiles U⁺(y⁺), dU⁺/dy⁺ at Re_τ = 182, 543, 1001, 1995, 5186; κ = 0.384 ± 0.004** (ours 0.3837 from y⁺ = 350, 0.3845 from 3√Re_τ; sublayer within 5e-4 for y⁺ < 1; Spalding and mixing-length channel "approximate") — cited, no licence statement; **redistributed as is by the user's decision of 2026-10-07** | Lee & Moser, J. Fluid Mech. 774, 395 (2015), doi:10.1017/jfm.2015.268 (arXiv:1410.7809); files at turbulence.oden.utexas.edu/channel2015/data | `reference/ch12/lee_moser_2015_channel_mean.csv` (281-row subset), `benchmarks.json`, `SOURCES.md`, `make_refs.py` | ch12 C10, C11, C12 (V5), E7, E8 |
| Prandtl's friction law 1/√f = 2.0 log₁₀(Re√f) − 0.8; Pope's c_L = 6.78; one-dimensional Kolmogorov constant 0.53 ± 0.055 (ours 0.491); k–ε standard constants; van Driest A⁺ = 26 — **all secondary** (the last three "constant consistency", not V5). **Not a benchmark**: the "Businger–Dyer" φ_m coefficients (unread first-hand) | Prandtl (1935) via McKeon et al., JFM 538, 429 (2005); Pope (2000) via arXiv:1705.04917 and the ATOMIX wiki; Sreenivasan, Phys. Fluids 7, 2778 (1995) via ATOMIX; Launder & Sharma (1974) via documentation pages | `reference/ch12/benchmarks.json`, `SOURCES.md` | ch12 C08, C11, C13, C15 |
| Our ch12 data (public, labelled ours): model-channel and van Driest tables (byte-reproducible) | `ch12.explainer_tables`, `write_reference_tables` | `reference/ch12/explainer_tables.json` | ch12 E5, E8 |
| **Eady problem: growth coefficient 0.3098, fastest wavenumber 1.606, cut-off (two digits)** (ours 0.30982, 1.60612, 2.39936) | K. A. Emanuel, MIT OpenCourseWare 12.803 (Fall 2009), Lecture 19 — read first-hand by the verifier, not re-read by the reviewer; NCL `eady_growth_rate` documentation as a second source | `reference/ch13/benchmarks.json`, `SOURCES.md`, `make_refs.py` | ch13 C16 (V5), E10 |
| Our ch13 data (public, labelled ours): cached runs `kelvin_basin`, `pv_particles`, `turbulence_f`, `turbulence_beta` (float16 `.npz`, regenerated within 4.1e-4 of range), explainer constants. **Not benchmarks**: Earth's rotation rate and radius (constant consistency); adjustment energy split, β-plane jet growth rates, −3 range, Rhines length (sources unread) | `scripts/ch13_make_caches.py`, `ch13.write_reference_runs` | `reference/ch13/*.npz`, `explainer_constants.json` | ch13 C11, C13, C17, E-parity rows |
| **Wagner's function, Jones' two-exponential fit** (coefficients 0.165, 0.0455, 0.335, 0.3; φ(0) = ½; error below 1 %) — a **lift** ratio; ours: coefficients identical, the lumped-vortex model's circulatory-lift ratio within 0.0060 over 2–20 chords; Garrick's fit 1 − 2/(4 + s) as a cross-check (0.022) | Dawson & Brunton, arXiv:2104.15122v2, Eq. (1)–(2), Table 1; BYU FLOW Lab Aeroelasticity.jl documentation — both read first-hand; **primary report (R. T. Jones, NACA, 1938) not read** | `reference/ch14/benchmarks.json`, `SOURCES.md`, `make_refs.py` | ch14 C03 (V5, cited fit), E3 |
| Forms (secondary sources, V5 as forms): thin-airfoil slope 2π per radian and centre of pressure at c/4; Kutta circulation of an off-centre circle; C_Di = C_L²/(π AR e) with e = 1 for an elliptic distribution; D_i in terms of lift and span. **Not benchmarks**: the cited share of induced drag of an airliner at cruise (about 37 %, evidence for trap T9); the exact start-up circulation (the reviewer's computation); the rectangular wing's e_span 0.95393 and slope 4.5304 at Λ = 6 (ours by convergence) | Wikipedia "Airfoil", "Joukowsky transform", "Lift-induced drag"; NASA Glenn Beginner's Guide to Aeronautics ("Induced Drag Coefficient") — read first-hand 2026-10-12 | `reference/ch14/benchmarks.json`, `SOURCES.md` | ch14 C02, C06, C14 |
| **NACA Report 1135, Tables I and II** (rows M = 0.50–0.60 and 1.60–2.10; isentropic, ν, μ and normal-shock columns) **and equations (105)–(171c)**, weak-shock and expansion series, wedge limit 45.6° — 109 of 111 entries to half a unit of the last digit, all to one unit; forms to 1e-10 … 1e-13 | Ames Research Staff, NACA Report 1135 (1953), public domain — read first-hand from the page images | `reference/ch15/benchmarks.json`, `SOURCES.md`, `make_refs.py` | ch15 C04, C05, C07, C08, C13–C15 (V5) |
| **Pries in-vitro relative viscosity** Eqs. (3)–(5), minimum near 7 µm (ours: identical form 1e-13, minimum 6.85 µm); FL closure 20–35 % below (qualitative) | Secomb & Pries, C. R. Physique 14, 470 (2013) — read first-hand | `reference/ch16/benchmarks.json`, `SOURCES.md`, `make_refs.py` | ch16 C05 (V5) |
| Forms (secondary, V5 as forms): Womersley u, Q, τ_w; Moens–Korteweg; Windkessel and its diastolic decay; Murray's cost; entrance length 0.0575 Re D, 1.359 D Re^{1/4} ((16.37) +4.3 % / ×0.51, S31); Dean equations in Navier–Stokes form (solved here: c₂ = 0.0305754, c₄ = 0.0119311). **Not reached**: Dean 1928, Womersley 1955/57 tables | Wikipedia "Pulsatile flow", "Moens–Korteweg equation", "Windkessel effect", "Murray's law", "Entrance length"; Imperial College M4A33 notes | `reference/ch16/` | ch16 C04, C07–C10, C12, C14 |
| Normal-shock rows M₁ = 1.5–3.0 (5e-5); p*/p₀ = 0.528 and eight gases; Fanno and Rayleigh forms (T_max at M = 0.845); Riemann test 1 p* = 0.303130, u* = 0.927453 — **secondary, forms** | Wikipedia "Normal shock tables", "Choked flow", "Fanno flow", "Rayleigh flow"; Toro (cited form) | `reference/ch15/` | ch15 C05, C07, C12, N114 |
Book-printed values (private, git-ignored): `tests/book_values_ch01.json` … `…_ch16.json` (worked examples, tables, figure constants,
exercise answers; `_forbidden_public` strings and regex patterns). **Open**: measured high-Re cylinder C_p (ch06), primary c_g,min (ch07),
primary pipe transition (ch08), Glauert wall-jet constants, Lienhard C_D, Roshko St, a primary 0.074 (ch09), Schäfer–Turek
(ch10), a counter-rotating Taylor number, primary Michalke digits (ch11), first-hand Businger–Dyer, primary k–ε, Prandtl,
van Driest (ch12), first-hand Gill, Kuo, Kraichnan, Rhines, Eady (ch13), the primary Jones / Wagner report, a public
lifting-line table, a measured lift curve (ch14), a first-hand Toro table and a similarity blast constant (ch15), **Dean
1928, Womersley 1955/57, Shapiro 1977a (ch16)**.

## Validation summary per chapter
Counts are concept-map rows carrying each label (a row can carry several); tests per evidence level in parentheses.
| chapter | analytic | symbolic | converged | conserved | benchmark | book-value | qualitative | unverified | verdict |
|---|---|---|---|---|---|---|---|---|---|
| ch01 | 39 (V1 58) | 16 (V2 19) | 5 (V3 8) | 5 (V4 7) | 18 (V5 13) | 1 (V6 5) | 1 | 0 | **PASS** — 107 tests; 5 explainers (115 rows); 496 cells |
| ch02 | 43 (V1 41) | 21 (V2 13) | 11 (V3 15) | 2 (V4 1) | 9 (V5 6) | 8 (V6 4) | 0 | 0 | **PASS** — 92 tests; 8/8 variants; 5 explainers (127 rows); 405 cells |
| ch03 | 38 (V1 58) | 23 (V2 32) | 14 (V3 14) | 5 (V4 4) | 2 (V5 1) | 8 (V6 4) | 0 | 0 | **PASS** — 120 tests; 13/13 (+7); 7 explainers (145 rows); 413 cells |
| ch04 | 41 (V1 83) | 35 (V2 41) | 17 (V3 9) | 13 (V4 6) | 11 (V5 6) | 5 (V6 3) | 0 | 0 | **PASS** — 151 tests; 26/26; 9 explainers (255 rows); 611 cells |
| ch05 | 39 (V1 90) | 24 (V2 23) | 13 (V3 10) | 13 (V4 8) | 0 (V5 0) | 5 (V6 4) | 2 (fenced) | 0 | **PASS** — 135 tests; 40/40 + 5; 9 explainers (268 rows); 496 cells |
| ch06 | 33 (V1 86) | 24 (V2 27) | 14 (V3 10) | 16 (V4 6) | 3 (V5 2) | 7 (V6 4) | 2 (labelled) | 0 | **PASS** — 137 tests; 30/30; 9 explainers (291 rows); 564 cells |
| ch07 | 38 (V1 60) | 28 (V2 36) | 17 (V3 10) | 8 (V4 5) | 5 (V5 5) | 3 (V6 3) | 0 (2 labelled) | 0 | **PASS** — 132 tests (+ V7 13); 41/41; 9 explainers (222 rows); 617 cells |
| ch08 | 29 (V1 46) | 20 (V2 31) | 6 (V3 5) | 7 (V4 4) | 7 (V5 6) | 0 rows (V6 3) | 0 | 0 | **PASS** loop 2 — 101 tests (+ V7 6); 41/41; 9 explainers (278 rows); 501 cells |
| ch09 | 26 (V1 37) | 13 (V2 34) | 8 (V3 17) | 6 (V4 9) | 9 (V5 16) | 2 (V6 5) | 5 (labelled) | 0 | **PASS** loop 1 — 136 items after review; 22/22 D; 25 variants; 9 explainers (248 rows); 472 cells |
| ch10 | 30 (V1 76) | 13 (V2 36) | 20 (V3 28) | 7 (V4 7) | 1 (V5 5) | 0 rows (V6 2) | 3 (labelled: block C_D, confined St/C_D, heuristic closures) | 0 | **PASS** loop 1 + review (4 Must: ψ_min sign, (10.199) stopping, MacCormack Δt, rule-9 leak) — 340 items (+ V7 20); 23/23 D sympy-checked; 38/38 variants; 8 explainers PASS round 2 (310 rows, 270 vs fluidpy); notebook PASS round 2 (602 cells, 72 s) |
| ch11 | 15 (V1 70) | 8 (V2 22) | 12 (V3 22) | 6 (V4 8) | 6 (V5 10) | 0 rows (V6 2) | 2 (labelled: Lorenz Lyapunov exponent; Falkner–Skan upper branch at Re ≥ 2e4) | 0 | **PASS** after 2 fix loops (TG decay box; k-dependent map scale) + review (2 Must: `salt_finger_regime`, Ri at U′ = 0) + 6 later defects (Rayleigh false zeros → complex path; OS far-field box → Bickley two bands; Blasius Re_c; Bénard wrong mode; contour `semi_infinite` refused; slip S13) — 160 tests (+ V7 41); 25/25 D sympy-checked (★★/★★★); 61 variants, 1 equivalent survivor; 9 explainers PASS round 2 (417 rows, 350 vs fluidpy); notebook PASS round 3 (618 cells, 126–239 s by load) |
| ch12 | 4 (V1 54) | 9 (V2 30) | 1 (V3 6) | 0 primary (V4 8) | 2 (V5 5) | 0 rows (V6 2) | 3 (labelled: `jet_tke_budget`, `k_epsilon_channel`, free-shear amplitude constants) | 0 | **PASS** on the first pass (117 tests, 42/42 variants) + review (no wrong equation; 4 Must: `Gamma_a` default, table constants in scripts and design, 11 stale docstring statements, an unread citation) — 173 tests (+ V7 10, statistical 10, constant consistency 2); every ★★/★★★ derivation re-derived independently (22 `_V2_derivation` tests); 42/42 + 18/18 variants; 10 explainers PASS round 2 (597 rows, 444 vs fluidpy); notebook PASS round 3 (690 cells, ≈ 240 s) |
| ch13 | 16 (V1 25) | 16 (V2 31) | 8 (V3 8) | 4 (V4 6) | 1 (V5 2) | 0 rows (V6 4) | 1 (labelled: 2-D turbulence caches) | 0 | **PASS** on the first pass (89 tests, 34/34 variants; no physics defect) + review (no wrong equation; 2 Must: `barotropic_run` start not projected, a WKB docstring order) — 91 tests (+ V7 15); all 27 ★★/★★★ derivations re-derived independently with sympy; 39/39 variants; 14 slips (7 with a planted printed form that fails), 18 traps; 10 explainers PASS round 2 (482 rows, 374 vs fluidpy); notebook PASS round 3 (821 cells, 124–165 s) |
| ch14 | 16 (V1 41) | 14 (V2 21) | 8 (V3 11) | 7 (V4 13) | 4 (V5 5) | 0 rows (V6 3) | 5 functions (labelled: separation estimate, generic lift and drag curves, stall regime, roll-up spiral) + 1 approximate, 1 model | 0 | **PASS** on the first pass (95 tests, 82/82 variants; no wrong equation) + review (1 Must: Wagner's function is a lift ratio, not the growth of the bound circulation — a mislabelled comparison no mutant could see; 15 Should; `_area` exact across a kink) — 95 tests (+ V7 28), 16/16 further variants; all 18 ★★/★★★ derivations re-derived independently with sympy; 10 slips (printed forms of #1, #2, #4, #8 fail their checks), 18 traps; 9 explainers PASS round 4 (328 rows, 240 vs fluidpy; rounds 2–4 held by library defects L1–L3, fixed); notebook PASS round 2 (549 cells, 126 s) |
| ch15 | 11 (V1 37) | 12 (V2 29) | 1 (V3 3) | 7 (V4 11) | 8 (V5 7) | 0 rows (V6 5) | 1 function (labelled: `blast_constant_crude`) + 2 sketches (`wedge_flow` detached, the plate's trailing-edge closure) | 0 | **PASS** after 1 fix loop (F1 NaN at M = 1, F2 σ < μ₁ accepted; 103/103 + 45/46 variants, 1 equivalent) + review (no wrong equation; 1 Must: one gas per call; 9 Should) — 111 tests (+ V7 19); three physics-changing slips coded corrected, the printed forms fail; 23 sympy re-derivations of every ★★/★★★ D row; NACA 1135 109/111; 17 slips, 22 traps; 9 explainers PASS (295 rows, 245 vs fluidpy; round 1: 2 Must in one file, fixed); notebook PASS round 2 (606 cells, 68 s) |
| ch16 | 16 (V1 36) | 13 (V2 27) | 7 (V3 9) | 5 (V4 4) | 7 (V5 4) | 0 rows (V6 1) | 1 (labelled: the Fåhræus–Lindqvist closure magnitude) | 0 | **PASS** after 1 fix loop (S34: (16.174) advection sign, c₄ 0.0119312; Windkessel overflow) + review (no formula, sign or unit error; 1 Must: (16.188) used above α = 1 → guard; 5 of 6 Should fixed, argument order deferred) + L3-1 (closed-form α*) — 99 tests (+ V7 15); 56/56 variants; every ★★/★★★ D row re-derived with sympy; 34 slips (9 physics slips coded corrected, six `printed_*` options fail), 27 traps; 10 explainers PASS (229 rows, 177 vs fluidpy; round 1: 4 lone-backslash Must in three files, fixed); notebook PASS round 2 (577 cells, 133 s) |
(ch13–ch16: CORE rows carrying each label. V7 tests: ch01 19 … ch14 28, ch15 19, ch16 15. Full suite: ch15 phase 5 2080
collected; **ch16 loop 1 2177 collected, 2176 passed (the 1 = S34, fixed)**, ≈ 44 min; run as `pytest -q tests`.)

## Open across chapters
- **ch16 open**: the top of this file and `knowledge/ch16.md` §0, §9. **ch15–ch11 open** (their §0, §9): docstring
  citations, crude blast constant, slip #17 (ch15); unread Wagner report (ch14); slip #14, β < 0 Rayleigh–Kuo (ch13);
  `skin_friction_zpg` 0.074 since ch09, secondary k–ε, Prandtl, van Driest (ch12); no Taylor benchmark for μ < 0 (ch11).
- **Library pass (orchestrator, when nothing reads `viz/` or `assets/`)**: ranked lists in `knowledge/viz_patterns.md` —
  **ch16**: TeX-safety checks, `Viz.num.cbessel`, `Viz.fitText`, `Viz.mix`/`ramp`, `Viz.arrowPx` + `v.size` (26 files carry
  a window-based `narrow()`), log axes, `axes2`; **ch15**: preset-exclusivity probe, `Viz.gas`; **ch14**: `Viz.cx`,
  `transport.max`/`holdAt`, `Viz.num.linalg`; **ch13**: marker glyphs, `.c-amber`; **ch12**: `Viz.tpow`, precise `erf`. Then
  `viz_inline.py --all`, shot quick on every chapter **and the click / page-then-step / preset / `.katex-error` probes**.
- **Skill pass**: ch14–ch16 lessons are in the skills (`interactive-viz`, `math-to-python`, `verify-implementation`,
  `teaching-style`, `colab-notebook`); the ch02–ch13 candidates in `knowledge/viz_patterns.md` are still not.
- **Machinery / tooling TODO**: the `shot.py` checks at the top of this file; `viz_lint` lone-backslash warning; `eq_refs`
  on derivation `start`/`result`; `testpaths = tests`; a keyword-call lint (sibling orders ch14–ch16); a derivation-step
  sync check (counts **and titles**); the label-audit script into `tools/`; a figure-note guard; word patterns for
  `check_public.py`; argparse in every `build_chNN.py`.
- **Explainer follow-ups**: per chapter in `knowledge/chNN.md` §9 and `viz_patterns.md`. Backups not built: ch04
  `kinematic_free_surface`, ch05 `vortex_rings`, ch06 `flow_net_sources_vortices`, ch07 `linearised_free_surface`, ch08
  `rotating_cylinders_couette`, ch09 `ball_swing_magnus`, ch10 `operator_splitting_theta`, ch11 `period_doubling_route`,
  ch12 `k_epsilon_decay`, ch13 `inertia_gravity_beams`, ch14 `wake_rollup`, ch15 `wave_steepening`, ch16
  `fahraeus_lindqvist_layer`.
