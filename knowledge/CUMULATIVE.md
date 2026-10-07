# CUMULATIVE knowledge — fluidpy
(Rewritten by the knowledge-keeper after every chapter. Every agent reads this first. Last rewrite: after ch12,
2026-10-07, knowledge pass on commit f0b48b0.)

## Open items carried from ch12 — read before starting ch13 (details `knowledge/ch12.md` §0, §9)
- **Needs a human**: the Businger–Dyer coefficients/attribution of φ_m(ζ) are unread first-hand (`verified_first_hand:
  false`; not a benchmark). **User decisions pending**: the Lee & Moser DNS subset is redistributed with its citation but
  no licence statement; constants coinciding with the book's table sit in two earlier LOCAL, unpushed commits of 2026-10-07.
- **Tooling**: `tools/shot.py` does not fail on `.katex-error` nodes and audits Explain only in default and step states;
  14 library findings and ranked promotion candidates in `knowledge/viz_patterns.md` (log axes first).
- **Code**: `time_average_variance_ou` has no fluidpy counterpart (E1 invariant rows only); `skin_friction_zpg` 0.074
  unverified since ch09; secondary sources only for the k–ε constants, Prandtl's law, van Driest's A⁺; core cap 0.09
  uncited; the ch01 library string "Γ < Γa" means Γ_met < Γ_d — change it across chapters?
- **Text**: slips #16–#18 unnumbered in three explainers; design Part A shows 6.75e-3 (κ = 0.4) where the notebook prints
  6.59e-3 (κ = 0.41), and D13 calls slip #18 a "trap"; Should-fix lists in `reports/ch12_lesson.md`, `reports/ch12_viz.md`.

## Project rules in force (added during ch01–ch03; sharpened since)
- **5–10 explainers per chapter, as many as the CORE ideas need** (ch03 7, ch04–ch09 9 each, ch10 8, ch11 9, **ch12 10**).
- **Every book equation is shown in full next to its number** — notebook, explainer text, `<meta>` strings **and the
  design file** (ch12: 214 number-only mentions inlined). **A number stands only beside the equation actually printed under
  it** (render the page); never start a sentence with a symbol. Enforced by `tools/coverage_check.py` check 8,
  `tools/eq_refs.py`, the builder's `self_check_near`; reviewers still scan by eye.
- **No control characters, no backslash-less TeX**: `coverage_check` check 9; the builder's `self_check_ctrl` and (ch12)
  `self_check_eaten`. **Write files with Edit/Write only — shell heredocs and sed destroy backslashes.** Explainers are
  still not linted for it.
- **Rule 9 (public repo) in practice**: book-quoted numbers live only in `tests/book_values_chNN.json` (valid JSON), with
  `"_forbidden_public"` strings and (ch12) `"_forbidden_public_regex"` patterns enforced by `tools/check_public.py`, which
  exits on an unparsable private file and **scans tracked files only — run it after `git add`**. **"Illustrative"
  constants must be visibly not the book's** (they coincided five times in ch12); a *computed* value can coincide with a
  table row too (pick another input); reports that quote findings can leak constants.
- **Pager and audit**: an item taller than a page is sliced at natural breaks; `tools/shot.py` pages through every pager at
  every size. The notebook (`build_chNN.py --dump`) or design Part F is the reference for derivation steps and titles.
- **Two values, the book's first** (ch09); **our own run parameters when the book's are private**, and a **DEVIATION box**
  wherever the code departs from the book (ch10).
- **Eigen-solver contract (ch11)**: leading mode or raise; a 0 or NaN in a table has a stated meaning; box
  max(profile box, 12/k) shown converged by box **and** N doubling; an independent route for every headline eigenvalue.
- **Contract discipline (ch12)**: design "expect" values are hypotheses computed before the code existed — builders report
  differences, they do not edit; after a crash, audit code against design Part C mechanically; **convention-carrying
  arguments (`Gamma_a`) and empirical constants (κ, B, Π, z₀) are required keywords**; "constant consistency" is not V5;
  an unread source is not a benchmark.

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
ch09 Boundary layers — δ ~ L Re^{−1/2}, (9.9)–(9.11); δ*, θ, H; Blasius f‴ + ½ff″ = 0 (0.332057); Falkner–Skan with a fold
     at n = −0.090429; momentum integral (9.43) → Thwaites (9.50); separation μu_yy = dp/dx; form drag; Kármán street
     b/a = 0.28055 (perturb–linearise–eigenvalues); drag crisis; free jet J const, sech²; wall jet Ψ const; teacup → Ekman
ch10 CFD — stencils = weighted Taylor series, order, round-off floor; FTCS (10.10), α, β, C = 2α; truncation error (10.17);
     von Neumann G (10.24), Noye 0 ≤ 4α² ≤ 2β ≤ 1 (10.27); upwind, CFL, Lax; weak form → Galerkin M ḋ + K d = F → element
     blocks; cell Péclet R_cell > 2 ⇒ wiggles, upwind = +uΔx/2 diffusion; MacCormack = Lax–Wendroff; splitting, projection
     ∇²p = ∇·u*/Δt; staggered C-grid, checkerboard; mixed FE, LBB (P2–P1 β_h ≈ 0.366); cavity vs Ghia 0.26 %; observed
     order + Richardson/GCI (ψ_min p = 2.00; a failed block study p = 0.43)   (D01–D23; details `knowledge/ch10.md`)
ch11 Instability — "basic state + small disturbance → linearise → normal modes → eigenvalue problem for the growth rate
     → neutral curve → its minimum is the critical number; the theorems are one-way; linear onset is not transition"
  normal mode (11.1), σ = −ikc (D01) [C01]; Kelvin–Helmholtz (11.18), ΔU_min 6.70 m/s (D02–D04) [C02]; Bénard: Ra (11.21),
     Ra_c 1707.76 rigid–rigid, 27π⁴/4 free–free (11.44) (D05–D12) [C03–C05]; salt fingers (11.46) (D13) [C06];
     Taylor–Couette Ta (11.52) (D14, D15) [C07]; Taylor–Goldstein (11.61), Miles–Howard Ri > ¼ ⇒ stable (11.67), Howard's
     semicircle (11.72) (D16–D19) [C08–C10]; Squire (11.78), Orr–Sommerfeld (11.79), Rayleigh (11.84), Fjørtoft (11.86),
     complex-path solver (D20–D22) [C11, C12]; Poiseuille Re_c 5772.22, Blasius 519.06 (δ*), Bickley 4.017 [C13];
     disturbance energy dE/dt = P − Λ (11.88) (D23) [C14]; Lorenz (11.91), r_H 24.737 (D24 ★★★, D25) [C15]
ch12 Turbulence (done) — "an average is a sum; it commutes with every linear operation and fails on exactly one thing,
     a product — the leftover covariance is the Reynolds stress, the production, the eddy flux; closures guess it"
  ensemble average (12.1), rules (12.4)–(12.9), mean(ũṽ) = ŪV̄ + mean(uv) (D01) [C01]; R₁₁(τ) (12.17), Λ_t (12.18),
     λ_t (12.19) (D02–D04) [C02]; S_e ↔ R₁₁ (12.20)–(12.22), S_e(0) = mean(u²)Λ_t/π (D05) [C03]; RANS (12.30), Reynolds
     stress −ρ₀mean(u_iu_j), closure problem (D06) [C04]; isotropy: g = f + (r/2)f′, ε̄ = 15ν mean(u²)/λ_g² (12.43)
     (D07, D08 ★★★) [C05]; energy budgets (12.46), (12.47): production is one term with two signs (D09, D10 ★★★) [C06];
     η = (ν³/ε̄)^{1/4}, η/L ~ Re^{−3/4} (12.50)–(12.51) (D11) [C07]; S₁₁ = C₁ε̄^{2/3}k₁^{−5/3} (12.54, slip #1) (D12) [C08];
     plane jet: J_s invariant (12.62), δ ∝ x, U_CL ∝ x^{−1/2} (12.66), seven flows by one recipe (D13–D16; D14 ★★★) [C09];
     U⁺ = f(y⁺) (12.80), U⁺ = y⁺ (12.82) (D17, D18) [C10]; log law by overlap matching (12.88), rough wall (12.93); DNS
     κ = 0.384 in a stated window (D19, D20) [C11]; eddy viscosity (12.94), mixing length → log law, B = −1.23 undamped,
     5.28 with van Driest damping (D21) [C12]; k–ε (12.103)–(12.105), κ_model = 0.433 (D22, D23) [C13]; Rf (12.107),
     Ri = Pr_T Rf (12.109), Γ_a required (D24) [C14]; L_M (12.110), Rf = z/L_M (12.111), log-linear wind (D25) [C15];
     Taylor dispersion (12.119): ballistic → diffusive, D_T → mean(u²)Λ_t for t ≫ Λ_t (12.129, slip #12) (D26–D28) [C16]

next: ch13 Geophysical Fluid Dynamics — `TS` for eddy fluxes, spectra and effective degrees of freedom; `WT` + ch12
      Monin–Obukhov functions for the surface layer and bulk drag (move to core there); `shear_flow_eddy_viscosity_solve`
      as the Ekman-with-K(z) template; `ST.rayleigh_eigs_contour` with U″ − β; `core.waves`, `core.rotating`, `core.mac`
      (C-grid); Kundu Γ computed, meteorological shown, `Gamma_a` explicit; primers start at P307.
```
The book ahead: **Ch13 GFD (uses 4, 5, 7–12)** → Ch14 aerodynamics · Ch15 compressible · Ch16 biofluids.
Where earlier chapters feed in: ch01 N², θ, lapse-rate conventions → **Ch. 13**; isentropic gas → Ch. 15. ch02 Stokes →
Ch. 13 PV. ch03 vortices → Ch. 13, 14. ch04 rotating frame, Boussinesq, Ri, Ro → **Ch. 13**; CV budgets → Ch. 14, 15.
ch05 **(5.30), (ζ + f)/h → Ch. 13**. ch06 `core.potential`, `conformal`, `panels` → **Ch. 14**; `laplace_solvers` → Ch. 13
PV inversion. ch07 `core.waves` → **Ch. 13** (√(gH), Poincaré/Kelvin/Rossby, adjustment), Ch. 15 (jump ↔ shock). ch08
`core.laminar` → **Ch. 13** (Stokes layer ≡ Ekman with ω → f, spin-up). ch09 `core.boundary_layer` → **Ch. 14**; teacup →
**Ch. 13 Ekman**; `jets` → Ch. 13 plumes. ch10 `core.fd` → **Ch. 13** (CFL for gravity/Kelvin waves, upwind tracers);
`core.fem1d` → Ekman ODE with K(z); **`core.mac` → Ch. 13 C-grid shallow water**; `core.maccormack` → **Ch. 15**. ch11
`core.stability` → **Ch. 13** (Rayleigh–Kuo = Rayleigh's equation with U″ → U″ − β; Eady/Charney by `cheb_grid` +
`constrained_eig`; Taylor–Goldstein; inertial instability = Rayleigh's criterion with absolute angular momentum; **Squire
fails with rotation or stratification**). **ch12** `core.turbstats` → **Ch. 13** (mean + eddy decomposition,
$\overline{v'T'}$, spectra, decorrelation times), Ch. 14; `core.wall_turbulence` + the Monin–Obukhov and Richardson
functions → **Ch. 13** (surface stress on the Ekman layer, C_D, stable boundary layers), Ch. 14; eddy coefficients and
Taylor dispersion → **Ch. 13** (K-closures, tracer mixing); D06 → averaging the rotating equations (the Coriolis term is
linear and averages to itself); D09–D10 → the energy cycle.

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
Chapter-only code (move to `core/` on second use): **ch12** `gradient_richardson_thermal(…, *, Gamma_a)`, `flux_richardson`,
`monin_obukhov_length`, `monin_obukhov_from_fluxes`, `surface_layer_state`, `stratified_tke_budget`, `turbulence_regime`
(→ **Ch. 13**: a `core/surface_layer.py`), `eddy_viscosity_stress`, `gradient_diffusion_flux`, `mixing_length_*`,
`shear_flow_eddy_viscosity_solve` (→ **Ch. 13** Ekman with K(z)), `taylor_dispersion*`, `eddy_diffusivity_*`,
`langevin_particles`, `richardson_diffusivity` (→ Ch. 13 tracers), `kolmogorov_scales`, `scale_table`, `model_spectrum`,
`batchelor_scale`, `convective_velocity_scale`, `mean_energy_budget`, `tke_budget`, `plane_jet_*`, `k_epsilon_*`,
`book_slips`; **ch11** `os_leading_mode`, `parallel_profile`, `gradient_richardson`, `miles_howard_stable`,
`rayleigh_number`, `benard_*`, `salt_finger_regime`, `rayleigh_circulation_criterion` (→ **Ch. 13**), `lorenz_*`; ch05
`column_relative_vorticity`, ch07 `wave_fields`, `stokes_drift`, `boussinesq_linear_sympy`, ch09
`secondary_flow_radial_force`, ch10 `theta_scheme_linear`, `cfl_time_step` (all → **Ch. 13**); ch06 `lift_per_span`
(→ Ch. 14); ch07 `simple_wave_evolve` (→ Ch. 15); ch10 `dominant_frequency` is superseded by `TS.periodogram`; the rest in
`knowledge/chNN.md` §3. **Default g**: ch04, ch05, ch07, `core.waves` G_BOOK = 9.81; the rest (incl. ch11, ch12)
G0 = 9.80665. **Default ρ**: ch06 1.2 vs 1000; ch07, ch08 1000; ch09 BL 1.2; ch10, ch11 non-dimensional; ch12 an argument.

### Explainer engine (`assets/viz_lib.js`)
`Viz.app` (tabs Walkthrough / Explore / Explain / Derivation / Equations / Code / Check; fit-to-window with density levels and
pagers; transport, modes, presets, status, terms, inspector, notes, selftest parity), `Plot`, `Viz.field`, `Viz.num` (RK4,
odeint, brentq, erf/erfc **~1e-7 only — a defect**, niceTicks), `Viz.work`, `Viz.three`, KaTeX with fallback; `Viz.font`,
`Viz.roundRect`, `Viz.text`, `Viz.card`, `Viz.fmtTime`, `Viz.tnum(keepTiny)` (values below 1e-12 print as 0 without it);
pager slicing, `ev.viewId`/`ev.vizView`. Quirks: Code-tab `{{placeholders}}` share one namespace; `#param=` deep links are
overridden by tour step 1; **`Viz.tnum(x) + '^2'` is a KaTeX double superscript when x is in scientific notation (ch12 —
bracket with a local `tsq`/`tpow`)**; `hidePortrait` keeps the view's flex share; no log axes, no linear algebra, no
complex numbers in the library. **No JS promoted after ch02–ch12** (the site-publisher was reading `viz/` each time); ch11
changed only `assets/viz_base.css`. Ranked candidate lists per chapter in `knowledge/viz_patterns.md` (ch12: log axes in
eight files, `tpow`, `fitText`, precise `erf`, Gauss–Legendre, seeded Gaussian/OU draws, log-paced transport, two-line
status, per-layout `hideOn`, two CSS fixes). Reference: `templates/viz_example.html`.

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

## Notation
See `knowledge/notation.md` (sign traps first, then the register, then per-chapter conventions). Everywhere: SI and kelvin
inside functions; **z up** (except ch06 §6.8); **lapse rate Kundu Γ ≡ dT/dz in code, meteorology shown**; `p0` ≠ `p_ref`;
pressures absolute unless `_gauge` (ch07 gauge); signed τ_xy = μ ∂u/∂y. ch02: passive C; traction on the first index; grid
`[k, j, i]`. ch03: R = G − Gᵀ, spin ½ω. ch04: Coriolis acceleration +2Ω × u′ vs force −2Ω × u′; **2-D ψ: u = ∂ψ/∂y** (ch07:
u = ∂ψ/∂z). ch05: γ = u_below − u_above. ch06: **Γ counterclockwise in code, `Gamma_cw=`/`Gamma_ccw=`**. ch07: **ω = angular
frequency**, η = elevation, g′ = g(ρ₂ − ρ₁)/ρ₂. ch08: **dp/dx book sign, `G=` = −dp/dx**; η = y/√(νt). ch09: x along the
wall; dp/dx > 0 adverse; θ = momentum thickness; FS exponent book n = code `m`; St with the cyclic frequency. ch10: three
grids (node / element / staggered cell); θ = Fourier angle per cell; α, β three ways; C = 2α Courant; Ma for Mach, M mass
matrix; slips "#1…#12" (R1–R12 in code). ch11: **σ = σ_r + iσ_i = −ikc is the growth rate** (surface tension σ_s); K Bénard
wavenumber vs k; **Γ three ways** ((11.21) −dT̄/dz, Kundu dT/dz in code, meteorology shown) and Γ = circulation in §11.6; Ra
signed in §11.5; `mu` = Ω₂/Ω₁; **Re per flow**; N² is (7.127); slips S1–S13 = "#1…#13".
**ch12 (changes):** over-bar = **ensemble** average; **lower-case u, v, w, T′ are fluctuations** (members on axis 0);
**κ = von Kármán** (`kappa`; thermal diffusivity `kappa_th`, eddy `kappa_T`); k₁, K wavenumbers, **ē** = turbulent kinetic
energy (`e`, the "k" of k–ε), ε̄ dissipation; λ_t, λ_f, λ_g Taylor microscales and Λ_t, Λ_f, Λ_g integral scales (ch11's Λ was
dissipation); **η = Kolmogorov length**; τ = time lag vs stress; f, g correlation functions vs wall function, Darcy f and
gravity; R_ij = correlation tensor; **h = full channel height**; θ potential temperature vs θ_m momentum thickness; σ² = 2νt
(ch03: 4νt); **heat flux positive upward, L_M > 0 stable**; **`uv_plus` keys hold MINUS the correlation** (alias
`minus_uv_plus`); **spectra two-sided, ∫S = variance**, `C1` one-sided; $\overline{u^2}$ = one component; T̄, T′ potential
temperature and **`gradient_richardson_thermal` REQUIRES `Gamma_a`** (Γ_a ≈ −9.8 K/km Kundu, Γ_d ≈ +9.8 K/km meteorology;
the library string "Γ < Γa" means Γ_met < Γ_d); `surface_layer_regime` says nothing about stability; the log-linear
unstable wind ends at z = |L_M|/5; slips "#1…#18".

## Teaching lessons
- **Depth tiered, coverage exhaustive** (A / B / C, derivations, cells): ch01 15/70/21, 12, 496 · ch02 16/60/18, 15, 405
  · ch03 15/52/12, 24, 413 · ch04 15/151/19, 30, 611 · ch05 14/69/11, 23, 496 · ch06 15/116/29, 31, 564 · ch07 16/193/27,
  37 (321 steps), 617 · ch08 15/113/15, 33 (301 steps), 501 · ch09 14/118/16, 22 (224 steps), 472 · ch10 15/105/23, 23
  (222 steps), 602 · ch11 15/131/24, 25 (258 steps), 618 · **ch12 16/210/17, 28 (266 steps), 690**.
- **A spine sentence, pointed back to in every block** (ch11 "one recipe, seven times"; ch12 "an average fails on exactly
  one thing — a product").
- **One-way theorems get a counter-example and the word "necessary" or "sufficient"** (ch11, P255); **"theorem vs
  observation" for look-alike criteria** (ch12: Ri > ¼ of (11.67) vs the observed Rf_cr ≈ ¼); "trust what does not move
  with N"; say what a table does not resolve next to its figure; a ★★★ check cell tests the derived result itself.
- **Convention callouts** ("which p_o?", Coriolis term vs force, Γ_cw vs Γ_ccw, which ω/η/g′/dp/dx/δ/Re/θ/λ/κ): state each
  convention, the size of the difference in numbers, which the code uses; **a conventions table up front for every
  overloaded letter, enforced in every later cell** (ch10, ch12); **a library string in another notation is translated at
  every appearance** (ch12: "Γ < Γa" read as Γ_met < Γ_d). **Both lapse-rate conventions on screen, always**: a two-row
  table, both inequalities in every legend, badge and slider trace; the thermometer trap run in code (+1.110 vs −2.213).
- **Book slips are taught by computing both versions**: a slips table up front ("slip #k"), printed forms as code options a
  test must fail, a slip box that says what is printed, why it fails, what is right — and never over-claims (ch12 slip #5:
  "the printed family fails", not "only power laws").
- **Make an approximation or an artefact measurable** (ch01, ch07–ch10; ch12: "Where this holds" under an order-of-magnitude
  estimate; the biased estimator's k⁻² tail beside the exact route with a one-term guide); **test an explanation by
  removing the alleged cause** (ch12 "leakage"; ch11 the Bickley "kink"). Pre-asymptotic honesty (ch10): grids with every
  observed order, a failed convergence study taught as failed; **a fitted constant is quoted with its window** (ch12).
- **Every key number by two independent routes; every prose number printed by a cell** — a printed gap before any
  comparative sentence; captions count what they claim; drawings draw what the reading note describes.
- **A tool first used inside a derivation gets its primer first**, and so does a named test case.
  `knowledge/primers.md` lists 307 (P01–P306 + P218a); **ch13 starts at P307**.
- **Builder rendering hygiene**: design TeX converted, not copied; raw strings; never a heredoc; `self_check_ctrl` /
  `self_check_eaten` / `self_check_near` tested on planted damage; templated comments (113 lines in ch12) and stray
  `<matplotlib.legend.Legend …>` outputs still appear — fix the template.
- **Statistics are estimates** (ch12): seeded cells, assertions at 5 standard errors (P281), synthetic fields captioned
  "kinematic — prescribed spectrum, no cascade".
- **Climate hooks with numbers** (ch04 Ro; ch07 swell, tsunami, √(g′H)/f; ch08 Ekman look-alike; ch09 teacup; ch10 CFL
  126 s ocean / 78 s atmosphere at 25 km, the Arakawa C-grid; ch11 billows, salt fingers, Rayleigh–Kuo, Lorenz; **ch12
  30-year normals as an averaging window, L_M = +83 m on a clear night and −17 m on a sunny afternoon, a 1 km model grid
  against η, a climate model's surface scheme, eddy diffusivity near a source**).
- Reusable derivation moves (ch01–ch12, 210 rows; 20 from ch12) are tabulated in `knowledge/concept_map.md`.

## Global pitfalls confirmed in this project
**Machinery and environment**
- Extracted text garbles maths (`¼` = `=`, `ð…Þ` = parentheses, missing minus, ω → `u`, γ → g, θ → q, σ/τ → s, ν → n,
  ε → 3, ρ → r, Ω → U, ψ → j, φ → 4, η → h, ζ → z, ∂ → v, ∇ → V): read page images.
- **Windows: heredocs and `sed` mangle backslashes, quotes and `|`** — write code, explainers, design and knowledge files
  with Write/Edit or Python files; builder strings with LaTeX must be raw; single-backslash TeX in JS strings turns `\f`,
  `\t`, `\r`, `\n`, `\b` into control characters (ch04, ch06, ch09–ch11; **ch12: form feeds in the design file broke an
  equation in the notebook**).
- `MPLBACKEND=Agg` strips inline figures; scripts take `--no-show`; Anaconda's kernelspec can shadow `fluidpy-venv`.
- **`tools/shot.py` blind spots**: it does not click or drag (stale Explain text after a slider drag); **it does not fail on
  `.katex-error` nodes or control characters, and audits Explain only in the default and step states** (ch11 E2; ch12 E1,
  E4 — builders and reviewers run their own DOM probe over presets × slider ends × tabs × pager pages); `py:` expressions
  have restricted builtins (no ndarray methods); tolerances near the achieved agreement; a table that rounds its own
  abscissa needs rtol 1e-4 (ch12).
- **Machine load**: intermittent pager clips (1–26 px), state-dependent clips; **run the merge gate serially** (full pytest
  ≈ 43 min at 1737 tests alone); **notebook runtime is load-dependent** (ch11 126–239 s; ch12 budget ≈ 170 s, committed
  run 241 s) — measure on an idle machine; animations are the main cost.
- **Parallel agents need private scratch subfolders and unique file names** (collisions in ch10, ch11, ch12); do not edit
  a solver module while a shot run evaluates parity rows; after an API overload or crash, audit code against the design
  contract mechanically (ch12).
- `coverage_check` matches ledger reminders by exact primer name; `eq_refs` misses numbers inside template strings and
  flags `ref:` labels with a suffix. Book wording can hide in docstrings, exercise answers in computed numbers, **run
  parameters and table constants in analysis, design, scripts and reports** (ch10, ch12).
- **Design and analysis documents contain errors; downstream agents compute, not copy** (every chapter; **ch12: the
  log-law fit expectation, an analysis row with the wrong sign of the jet stress, κ = 0.4 vs 0.41 in a worked number**).

**Mathematics → code**
- **Defaults encode physics** (ch09 `cp_base=None`; ch10 a default Δt that dropped the viscous rate; **ch12 `Gamma_a=0.0`
  called the standard atmosphere unstable**): convention-carrying and empirical arguments are required keywords; explicit
  solvers compute the full stability limit and raise.
- **Fitted constants belong to their window** (ch12: a correct fit returns κ = 0.383 on a composite profile from y⁺ = 30
  because the profile reaches the logarithm from below; DNS 0.400 from 30, 0.384 from 350).
- **Names carry signs and normalisations** (ch12: `uv_plus` holds −mean(uv) → alias `minus_uv_plus`; `C1` is one-sided
  and the two-sided default halves it; "to the first zero" is not an integral scale when r has a negative lobe — 26 %).
- **Estimators have artefacts of predictable size** (ch12: the biased correlation estimate's k⁻² tail); sampled assertions
  need a seed and 5 standard errors; an OU signal has no Taylor microscale; random-phase fields are kinematic; 2-D
  isotropy differs from 3-D.
- **Classifiers and verdict strings need tests on both sides of every boundary and at degenerate inputs** (ch11: 128
  tests passed with two classifier bugs; ch12: boundaries compared on the ratio and mirrored identically in JS; a regime
  name is not a stability verdict; NaN outside a formula's range). **A pinned number can pin a defect**: pin sets and physics.
- **Unbounded profiles (ch11)**: scale the truncation box with 1/k (`ST.decay_box`), prove it by box and N doubling; **a
  convergence filter can make a false zero or return a lower mode** — leading mode or raise; near a neutral point use a
  complex path and an independent route; spectral convergence is digits, with a round-off rise at large N; track mode
  families by continuation, never join across NaN (`knowledge/ch11.md` §4).
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
- **A telescoping sum is not a conservation test**; genuine V4 invariants: ch07 ∫η², KdV; ch08 thin-film volume; ch09 J,
  Ψ; ch10 projection ∇·u, periodic mass; ch11 the energy budget of a mode; **ch12 Parseval, ∇·u of synthetic fields,
  2ν∫K²E dK = ε̄, J_s and the slot-fluid flux, channel work = dissipation + production**.
- **Sibling functions share defaults and argument meanings** or are called by keyword; convention-named keywords
  (`Gamma_cw=`, `ref=`, `G=`, `closure=`, `printed=`, `convention=`, `scheme=`, **`Gamma_a=`, `upto=`, `two_sided=`,
  `unstable=`, `window=`**).
- **Benchmarks from the primary source at its printed precision, read first-hand** (ch10 Ghia by two transcriptions;
  **ch12: an unread "Businger–Dyer" citation was relabelled V1 and flagged**). **Validation labels**: V4 only for real
  conservation laws; "converged" only with an asserted order; "benchmark" only inside the stated tolerance; **"constant
  consistency" for a constant compared with its citation (not V5)**; "approximate" for models against data with the band
  stated; "secondary" for compilations; "qualitative" for confined/heuristic studies.

## Benchmark inventory
| value / table | source (citation) | stored in | used by |
|---|---|---|---|
| k_B, N_A, R (exact); e | CODATA 2018/2022, NIST CUU | `reference/ch01/`, `reference/ch08/` | ch01 `thermo`; ch08 Millikan |
| USSA-1976, sea-level c = 340.29 m/s; Sutherland | NASA-TM-X-74335; PDAS tables | `reference/ch01/` | ch01, ch04, ch08 ν_air |
| IAPWS σ(T), μ(T) | IAPWS R1-76(2014), R12-08 | `reference/ch01/` | ch01, ch04, ch07, ch08 |
| Jennings mean free path; Taylor blast 0.856; AMS Γ_d; EOS-80 | Tsalikis et al. 2024; arXiv:2009.05674; AMS Glossary; Fofonoff & Millard 1983 | `reference/ch01/benchmarks.json` | ch01 |
| Divergence theorem, Levi-Civita, Mohr, rotations, Stokes | Wikipedia; scipy docs | `reference/ch02/` | ch02 |
| Lamb–Oseen peak α = 1.256 | arXiv:2210.05223 | `reference/ch03/` | ch03 C14 |
| Sphere drag C_D(Re) | F. A. Morrison (2016) | `reference/ch04/` | ch04 C15; ch08; ch09 C11 |
| WGS-84; orifice C_c 0.611; capillary length | Wikipedia | `reference/ch04/` | ch04, ch05 |
| Sphere added mass ½; Rayleigh collapse 0.91468; source panels | Lamb (1932); Rayleigh 1917; Hess & Smith 1967 | `reference/ch06/` | ch06 |
| Fenton & McKee; Guo; capillary minimum; Stokes drift H/λ = 0.1410633 | Coastal Eng. 1990, 2002; Wikipedia; arXiv:1507.02784 | `reference/ch07/` | ch07 |
| Slider K_opt 2.1889; viscous gravity current η_N 1.411 | San Andrés (Texas A&M); Huppert, JFM 121 (1982) | `reference/ch08/` | ch08 C07, C08 |
| Blasius f″(0) = 0.332057336215196 (ours 5.7e-15) | Wikipedia "Blasius boundary layer" | `reference/ch09/` | ch09 C04 |
| Falkner–Skan κ, separation β = −0.198837735, reversed branch | Belden et al., arXiv:1907.09912 | `reference/ch09/` | ch09 C05 |
| Hiemenz 1.232588, δ*/δ 0.6479; Bickley jet 0.4543, 0.2752, 3.3019, 0.5503; Kármán b/a ≈ 0.281 | Weidman & Turner; Wikipedia; Horváth et al. 2020 | `reference/ch09/` | ch09 C05, C12, C10 |
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
| **Channel DNS mean profiles U⁺(y⁺), dU⁺/dy⁺ at Re_τ = 182, 543, 1001, 1995, 5186; κ = 0.384 ± 0.004** (ours 0.3837 from y⁺ = 350, 0.3845 from 3√Re_τ; sublayer within 5e-4 for y⁺ < 1; Spalding and mixing-length channel "approximate") — **cited, no licence statement: user decision pending** | Lee & Moser, J. Fluid Mech. 774, 395 (2015), doi:10.1017/jfm.2015.268 (arXiv:1410.7809); files at turbulence.oden.utexas.edu/channel2015/data | `reference/ch12/lee_moser_2015_channel_mean.csv` (281-row subset), `benchmarks.json`, `SOURCES.md`, `make_refs.py` | ch12 C10, C11, C12 (V5), E7, E8 |
| Prandtl's friction law 1/√f = 2.0 log₁₀(Re√f) − 0.8; Pope's c_L = 6.78; one-dimensional Kolmogorov constant 0.53 ± 0.055 (ours 0.491); k–ε standard constants; van Driest A⁺ = 26 — **all secondary** (the last three "constant consistency", not V5). **Not a benchmark**: the "Businger–Dyer" φ_m coefficients (unread first-hand) | Prandtl (1935) via McKeon et al., JFM 538, 429 (2005); Pope (2000) via arXiv:1705.04917 and the ATOMIX wiki; Sreenivasan, Phys. Fluids 7, 2778 (1995) via ATOMIX; Launder & Sharma (1974) via documentation pages | `reference/ch12/benchmarks.json`, `SOURCES.md` | ch12 C08, C11, C13, C15 |
| Our ch12 data (public, labelled ours): model-channel and van Driest tables (byte-reproducible) | `ch12.explainer_tables`, `write_reference_tables` | `reference/ch12/explainer_tables.json` | ch12 E5, E8 |
Book-printed values (private, git-ignored): `tests/book_values_ch01.json` … `…_ch12.json` (ch12: the free-shear table's
constants and half-widths, per-flow (κ, B) pairs, wake strength, worked-example inputs and answers; `_forbidden_public`
strings and `_forbidden_public_regex` patterns). **Open**: measured high-Re cylinder C_p (ch06), primary c_g,min (ch07),
primary pipe transition (ch08), Glauert wall-jet constants, Lienhard C_D, Roshko St, **a primary turbulent-plate
coefficient 0.074** (ch09, still open in ch12), Schäfer–Turek DFG cylinder benchmark, Dennis & Chang Re = 40 (ch10), a
citable counter-rotating Taylor number and primary Michalke digits (ch11), **first-hand Businger et al. 1971 / Dyer 1974,
primary k–ε, Prandtl and van Driest papers, a public table of free-shear constants (ch12)**.

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
(Counts are CORE rows by primary label for ch12. V7 tests: ch01 19, ch02 12, ch03 7, ch04 3, ch06 2, ch07 13, ch08 6,
ch09 24, ch10 20, ch11 41, ch12 10. Full suite at the ch12 gate: **1793 tests**; 1737 took 2564 s serially.)

## Open across chapters
- **ch12 open**: the list at the top of this file and `knowledge/ch12.md` §0, §9 (also: `@slow` tests without mutants;
  animations trimmed to A2 45, A3 60, A5 60 frames at dpi 60; the committed notebook carries no outputs until publish).
- **ch11 open** (`knowledge/ch11.md` §9): no citable Taylor benchmark for μ < 0; Bickley lower edges below k = 0.02;
  `benard_growth_rate` explicit N ≥ 60 round-off; complex-path limits; `tools/viz_tables/ch11/` absolute paths.
- **Library pass (orchestrator, when nothing reads `viz/` or `assets/`)**: ranked lists in `knowledge/viz_patterns.md` —
  **ch12**: log axes for `v.plot()` (eight files), `Viz.tpow`, `fitText`, precise `erf`, Gauss–Legendre, seeded
  Gaussian/OU draws, log-paced transport, two-line status, `hideOn`, two CSS fixes; **ch11**: complex numbers,
  `bisect`/`goldenMax`, table interpolation, cubic solver, `arrowPx`; **ch10**: `Viz.num.linalg`, `gutterPlot`. Then
  `tools/viz_inline.py --all`, `tools/shot.py --chapter ch01…ch12 --quick` and `templates/viz_example.html --quick`.
- **Skill pass**: lesson candidates for the five skills in `knowledge/viz_patterns.md` (ch02–ch12) are not yet in the
  skills (ch12 first: bracket powers and probe `.katex-error`; required keywords for conventions and empirical constants;
  constant consistency ≠ V5; a number only beside its own equation; test an explanation by removing its cause).
- **Machinery / tooling TODO**: **`shot.py` failing on `.katex-error` and control characters and walking presets with
  Explain open** (three chapters); hash deep links after `goStep(0)`; code-line length lint ≈ 44 chars; `eq_refs`
  allow-marker and `ref:` suffixes; `viz_lint` control-character rule (seven chapters); a slider-drag Explain test;
  `check_public.py` on untracked files; a contract-audit tool (design Part C vs `inspect.signature`); builder self-checks
  and a templated-comment lint in `tools/nbkit.py` (four chapters).
- **Explainer follow-ups**: per chapter in `knowledge/chNN.md` §9 and `viz_patterns.md`. Backups not built: ch04
  `kinematic_free_surface`, ch05 `vortex_rings`, ch06 `flow_net_sources_vortices`, ch07 `linearised_free_surface`, ch08
  `rotating_cylinders_couette`, ch09 `ball_swing_magnus`, ch10 `operator_splitting_theta`, ch11 `period_doubling_route`,
  ch12 `k_epsilon_decay`. Older verification notes: ch10 block C_D non-asymptotic; ch09 sphere Re_cr, turbulent 0.074;
  ch08 O5, O8, O9; ch07 O5–O7; ch06 O4, O6; ch05 O2.
