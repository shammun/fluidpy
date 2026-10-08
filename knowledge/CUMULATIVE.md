# CUMULATIVE knowledge — fluidpy
(Rewritten by the knowledge-keeper after every chapter. Every agent reads this first. Last rewrite: after ch13,
2026-10-08, knowledge pass on commit 0a7c012.)

## Open items carried from ch13 — read before starting ch14 (details `knowledge/ch13.md` §0, §9; ch12 §0)
- **Needs the user**: **ch13 slip #14** (westward flow over a step, §13.13) was promoted from trap to slip after four
  independent confirmations and keeps its "true for a ridge of finite width" caveat — no ruling yet. From ch12: the
  Businger–Dyer coefficients/attribution of φ_m(ζ) are unread first-hand (not a benchmark).
- **Sources unread first-hand (ch13, all labelled "ours")**: Gill 1982 (adjustment, ⅓), Kuo (β-plane jet), Kraichnan (−3
  range), Rhines; the Eady benchmark (Emanuel, MIT OCW) was read by the verifier but not re-read by the reviewer.
- **Code**: `rayleigh_kuo_eigs` approximate (≈ 5 %) for β < 0; `vertical_modes(method="cheb")` loses digits of the
  barotropic speed; `scripts/ch13_instability.py` uses other ocean inputs than the notebook and explainer; argument order
  of `eady_max_growth_rate(f, N, …)` and `poincare_omega(K, f, c)` differs from their siblings (call by keyword). From
  ch12: `time_average_variance_ou` has no fluidpy counterpart; plotly title/legend collision in `core/interact.py`; the
  ch01 library string "Γ < Γa" means Γ_met < Γ_d.
- **Tooling**: `tools/shot.py` audits neither 844×345 nor `.katex-error` (fourth chapter); the pager leaves headings
  alone on a page and measures `eq` cards before their live rows are filled; candidates in `knowledge/viz_patterns.md`.
- **Text**: ch13 lesson Should-fix (3) and five figures never compared with their notes; ch13 viz Should-fix (velocity has
  no single colour; sparse phone pages); ch12 Derivation-tab clipping at 844×345; ch12 slips #16–#18 unnumbered.

## Project rules in force (added during ch01–ch03; sharpened since)
- **5–10 explainers per chapter, as many as the CORE ideas need** (ch03 7, ch04–ch09 9, ch10 8, ch11 9, ch12 10, **ch13 10**).
- **Every book equation is shown in full next to its number** — notebook, explainer text, `<meta>` strings and the design
  file; **a number stands only beside the equation actually printed under it** (render the page); never start a sentence
  with a symbol. Enforced by `tools/coverage_check.py` check 8, `tools/eq_refs.py`, the builder's `self_check_near`.
- **No control characters, no backslash-less TeX** (`coverage_check` check 9); **write files with Edit/Write only.**
- **Rule 9 (public repo)**: book-quoted numbers live only in `tests/book_values_chNN.json` (`_forbidden_public` strings and
  regex patterns, enforced by `tools/check_public.py` on tracked or named files). "Illustrative" constants must be visibly
  not the book's; **a book value written in words ("a quarter below") is still a book value — the tool cannot see it** (ch13).
- **Pager and audit**: `tools/shot.py` pages through every pager at every size; the notebook (`build_chNN.py --dump`) or
  design Part F is the reference for derivation steps. **Two values, the book's first** (ch09); a **DEVIATION box**
  wherever the code departs from the book (ch10; ch13: every numerical scheme is ours).
- **Eigen-solver contract (ch11)**: leading mode or raise; a 0 or NaN has a stated meaning; box and N doubling; an
  independent route for every headline eigenvalue (ch13: closed form + Chebyshev; collocation + shooting).
- **Contract discipline (ch12, ch13)**: **design Part C is written first and reconciled mechanically by the implementer**;
  design "expect" values are hypotheses (builders report differences; a correction pass follows the implementation);
  convention-carrying arguments (`Gamma_a`) and empirical constants are required keywords; "constant consistency" is not
  V5; an unread source is not a benchmark.
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
ch12 Turbulence — "an average is a sum; it commutes with every linear operation and fails on exactly one thing, a
     product — the leftover covariance is the Reynolds stress, the production, the eddy flux; closures guess it"
  ensemble average (12.1), rules (12.4)–(12.9) (D01) [C01]; R₁₁(τ), Λ_t, λ_t (12.17)–(12.19) [C02]; S_e ↔ R₁₁
     (12.20)–(12.22) (D05) [C03]; RANS (12.30) (D06) [C04]; isotropy, ε̄ = 15ν mean(u²)/λ_g² (12.43) [C05]; energy budgets
     (12.46), (12.47): production is one term with two signs (D09, D10) [C06]; η = (ν³/ε̄)^{1/4} (12.50) [C07]; S₁₁ =
     C₁ε̄^{2/3}k₁^{−5/3} (12.54, slip #1) [C08]; plane jet, seven flows by one recipe [C09]; U⁺ = f(y⁺) (12.80), log law
     (12.88), DNS κ = 0.384 in a stated window [C10, C11]; mixing length, k–ε [C12, C13]; Rf (12.107), Ri = Pr_T Rf
     (12.109), L_M (12.110) [C14, C15]; Taylor dispersion (12.119) [C16]   (details `knowledge/ch12.md`)
ch13 Geophysical fluid dynamics (done) — "one set of equations on a thin rotating shell with terms switched off: every
     balance keeps two or three terms, every wave is a disturbed balance, c/|f| is where rotation matters as much as
     gravity, and (ζ + f)/h is what a column remembers"
  thin-shell equations (13.9), f = 2Ω sin θ (13.8), β-plane (13.10) (D01, D02) [C01]; geostrophy (13.11)–(13.12), Ro (13.13)
     (D03) [C02]; thermal wind (13.15), Taylor–Proudman (13.21) (D04, D05) [C03]; Ekman spiral d²V/dz² = (if/ν_v)V (13.27),
     δ = √(2ν_v/f) (13.29) (D06 ★★★) [C04]; transport −τ/(ρf) (13.30), pumping (D07, D09) [C05]; bottom layer (13.41),
     maximum 1.067 U at 3πδ/4 (D08) [C06]; shallow water (13.45) on a C-grid (D10) [C07]; vertical modes (13.56), c_n² ≡ gH_e
     (13.62), weight-1 orthogonality for both lids (D11 ★★★) [C08]; the cubic ω³ − c²ωK² − f₀²ω − c²βk = 0 (13.76) (D12
     ★★★, D13) [C09]; Poincaré ω² = f² + gHK² (13.82) (D14) [C10]; Kelvin wave (13.87) (D15) [C11]; Rossby radius,
     adjustment with ⅓ of the energy kept (D16 ★★★, ours) [C12]; PV (13.94), flow over a step (slip #14) (D17 ★★★, D18)
     [C13]; inertia–gravity waves (13.96), (13.112) (D19 ★★★, D20, D21) [C14]; Rossby waves (13.117)–(13.120), Rayleigh–Kuo
     β − U″ (13.124) (D22 ★★★, D23, D24) [C15]; Eady (13.136), (13.141): cut-off 2.39936, fastest 1.60612, σ = 0.30982
     fU₀/(NH) (D25, D26 ★★★, D27) [C16]; Fjørtoft (13.145), −3 range, Rhines length (D28, D29) [C17]

next: ch14 Aerodynamics — ch06 `core.potential`, `conformal`, `panels`, `lift_per_span`; ch05 `biot_savart` (lifting line);
      ch09 `core.boundary_layer`, ch12 `core.wall_turbulence` (skin friction); ch13 complex velocity (Ekman's V = u + iv;
      potential flow uses u − iv), `SW.barotropic_run` (β = 0) for 2-D vortices, zero-sum term budgets; primers from P336.
```
The book ahead: **Ch14 aerodynamics (uses 3–6, 9, 12)** → Ch15 compressible · Ch16 biofluids.
Where earlier chapters feed in: ch01 isentropic gas, sound speed → **Ch. 15**. ch03 vortices, ch05 Kelvin's theorem,
Biot–Savart, sheets → **Ch. 14** (starting vortex, lifting line, wake roll-up). ch04 CV budgets, Bernoulli → **Ch. 14**, 15.
ch06 `core.potential`, `conformal`, `panels` → **Ch. 14**. ch07 `core.waves`, jump ↔ shock → Ch. 15. ch08 `core.laminar` →
Ch. 16. ch09 `core.boundary_layer`, ch12 `core.wall_turbulence` → **Ch. 14**. ch10 `core.maccormack`, ch11 `core.stability`,
ch13 `core.shallow_water` (hydraulic analogy) → **Ch. 15**; `core.vertical_modes` → duct and pulse-wave modes.

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
`ST.rayleigh_eigs_contour` gained the keyword **`beta=0.0`** in ch13 (Rayleigh–Kuo; default path unchanged to the last bit).
Chapter-only code (move to `core/` on second use): **ch13** `flow_over_step`, `rayleigh_kuo_eigs` (β ≥ 0 validated),
`eady_numeric_eigs`, `eady_matrix`, `eady_basic_state`, `jet_section`, `pressure_centre`, `coastal_upwelling`,
`lee_wave_field`, `vertical_structure_solve`, `wkb_error`, `thermocline_N2`, `ocean_profile_idealized`,
`lapse_rate_table(dT_dz, *, Gamma_a)`, `wind_from_to`, `illustrative_inputs`, `conventions_table`, `book_slips`, `traps`,
`reference_run` / `load_reference_run` (float16 caches), twelve sympy engines; **ch12** `gradient_richardson_thermal(…, *,
Gamma_a)`, `flux_richardson`, `monin_obukhov_*`, `surface_layer_state`, `eddy_viscosity_stress`, `mixing_length_*`,
`taylor_dispersion*`, `kolmogorov_scales`, `plane_jet_*`, `k_epsilon_*`; **ch11** `os_leading_mode`, `parallel_profile`,
`gradient_richardson`, `benard_*`, `lorenz_*`; ch06 `lift_per_span` (→ **Ch. 14**); ch07 `simple_wave_evolve` (→ Ch. 15);
the rest in `knowledge/chNN.md` §3. **Default g**: ch04, ch05, ch07, `core.waves` G_BOOK = 9.81; the rest (incl. ch11–ch13)
G0 = 9.80665. **Default ρ**: ch06 1.2 vs 1000; ch07, ch08 1000; ch09 BL 1.2; ch10, ch11 non-dimensional; ch12, ch13 arguments.

### Explainer engine (`assets/viz_lib.js`)
`Viz.app` (tabs Walkthrough / Explore / Explain / Derivation / Equations / Code / Check; fit-to-window with density levels and
pagers; transport, modes, presets, status, terms, inspector, notes, selftest parity), `Plot`, `Viz.field`, `Viz.num` (RK4,
odeint, brentq, erf/erfc **~1e-7 only — a defect**, niceTicks), `Viz.work`, `Viz.three`, KaTeX with fallback; `Viz.font`,
`Viz.roundRect`, `Viz.text`, `Viz.card`, `Viz.fmtTime`, `Viz.tnum(keepTiny)` (values below 1e-12 print as 0 without it);
pager slicing, `ev.viewId`/`ev.vizView`. Quirks: Code-tab `{{placeholders}}` share one namespace; `#param=` deep links are
overridden by tour step 1; `Viz.tnum(x) + '^2'` is a KaTeX double superscript in scientific notation (bracket with a local
`tsq`/`tpow`); `hidePortrait` keeps the view's flex share; view titles are not re-evaluated on a tab change; a walkthrough
`eq` card is measured before its live rows are filled (intermittent phone clip); no `.c-amber` class; no log axes, no
linear algebra, no complex numbers in the library. **No JS promoted after ch02–ch13** (the site-publisher was reading `viz/`
each time). Ranked candidates in `knowledge/viz_patterns.md` (ch13 top five: shot.py 844×345 + `.katex-error`,
`Viz.arrowPx`, a y-title strip, log axes, pager "keep a heading with the next block"). Reference: `templates/viz_example.html`.

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
**ch13 (changes):** x east, y north, z up; **f signed** (north +, south −; abs(f) for scales, sign(f) for directions,
`ValueError` at f = 0); every "to the right / clockwise / coast on the right" of the book assumes f > 0; **which z = 0
changes by section** (sea surface §13.6, §13.9, §13.14; solid surface §13.7; flat bottom §13.8; lower lid §13.17); p and ρ
are perturbations from §13.4 on; **ψ with u = −∂ψ/∂y** (opposite to ch04, ch06, ch11); **Ro = U/(abs(f)L)** (ch04:
U/(2ΩL)); **three Rossby radii** (external c/abs(f), internal with π, Eady NH/abs(f) without π); δ = Ekman e-folding
thickness (Ekman depth πδ); **V = u + iv** (potential flow: u − iv); η = surface displacement again; c = long-wave speed,
sound speed or complex phase speed by section; α three ways; θ latitude vs wavevector angle; **Rossby functions return
signed ω for signed k**; the cubic returns NaN where its discriminant is negative; modes ψ_n(0) = 1, rigid lid: index 0 is
the first baroclinic mode; **§13.18 spectra one-sided, no ½**; sidereal Ω; slips "#1…#14" with `kind`, traps T1–T18.

## Teaching lessons
- **Depth tiered, coverage exhaustive** (A / B / C, derivations, cells): ch01 15/70/21, 12, 496 · ch02 16/60/18, 15, 405
  · ch03 15/52/12, 24, 413 · ch04 15/151/19, 30, 611 · ch05 14/69/11, 23, 496 · ch06 15/116/29, 31, 564 · ch07 16/193/27,
  37 (321 steps), 617 · ch08 15/113/15, 33 (301), 501 · ch09 14/118/16, 22 (224), 472 · ch10 15/105/23, 23 (222), 602 ·
  ch11 15/131/24, 25 (258), 618 · ch12 16/210/17, 28 (266), 690 · **ch13 17/154/47, 29 (253 steps), 821**.
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
  derivation that needs them. `knowledge/primers.md` lists 336 (P01–P335 + P218a); **ch14 starts at P336**.
- **Builder rendering hygiene**: design TeX converted, not copied; raw strings; never a heredoc; self-checks tested on
  planted damage; templated comments and stray `<matplotlib.legend.Legend …>` outputs still appear — fix the template.
- **Statistics are estimates** (ch12): seeded cells, assertions at 5 standard errors, synthetic fields captioned "kinematic".
- **Climate hooks with numbers** (ch04 Ro; ch07 swell, tsunami, √(g′H)/f; ch10 CFL on an ocean grid; ch12 30-year normals,
  L_M by day and night; **ch13: reading a weather map, the jet above the strongest temperature gradient, upwelling from
  τ/(ρf), a first baroclinic Rossby radius of 43 km against a model grid, a year for a long Rossby wave to cross 10 000 km
  at 12°, storms growing in 1.6 days, the Rhines scale and zonal jets**).
- Reusable derivation moves (ch01–ch13, 225 rows; 16 from ch13) are tabulated in `knowledge/concept_map.md`.

## Global pitfalls confirmed in this project
**Machinery and environment**
- Extracted text garbles maths (`¼` = `=`, `ð…Þ` = parentheses, missing minus, ω → `u`, γ → g, θ → q, σ/τ → s, ν → n,
  ε → 3, ρ → r, Ω → U, ψ → j, φ → 4, η → h, ζ → z, ∂ → v, ∇ → V): read page images.
- **Windows: heredocs and `sed` mangle backslashes, quotes and `|`** — write code, explainers, design and knowledge files
  with Write/Edit or Python files; builder strings with LaTeX must be raw; single-backslash TeX in JS strings becomes
  control characters (ch04, ch06, ch09–ch12). **A scratch script named like a standard-library module (`numbers.py`)
  breaks `import numpy`** (ch13). `MPLBACKEND=Agg` strips inline figures; Anaconda's kernelspec can shadow `fluidpy-venv`.
- **`tools/shot.py` blind spots**: it does not click or drag; **it does not fail on `.katex-error` nodes or control
  characters, audits Explain only in the default and step states, and does not audit the in-page frame 844×345** (ch11–
  ch13: builders and reviewers run their own probes); `py:` expressions have restricted builtins (no ndarray methods;
  **plain lambdas cannot see the chapter module**, ch13); a table that rounds its own abscissa needs rtol 1e-4.
- **Machine load**: intermittent and state-dependent pager clips (**a walkthrough `eq` card with live rows**, ch13); run
  the merge gate serially (the not-slow suite took 2994 s at 1865 tests on a shared machine); notebook runtime is
  load-dependent (ch13 124–165 s unloaded) — measure on an idle machine.
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
- **Rotation and sign (ch13)**: every function takes either sign of f (abs for scales, sign for directions) and raises at
  f = 0; plant mutants that drop sign(f), flip a path in y or swap the north and south wall (two survived at first).
- **De-aliasing (ch13)**: apply the 2/3 mask to the state as well as to the product (white noise drifted +18 % in energy
  without it), with a strict inequality. **A new parameter is its own keyword**: β folded into U″ moved the solver's path.
- **Small roots and tolerances (ch13)**: compute a cubic's slow root from its own balance and check Vieta's product; return
  NaN and expose the discriminant where "all roots real" fails; **`numpy.allclose` has an absolute tolerance of 1e-8,
  larger than β** — relative tolerances only. **Measure before writing**: WKB error ∝ 1/(Hm), drift order 5, the Ekman
  maximum at 3πδ/4; a box-truncated eigenproblem fails silently on a radiating mode; non-dissipative schemes ring.

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
| **Channel DNS mean profiles U⁺(y⁺), dU⁺/dy⁺ at Re_τ = 182, 543, 1001, 1995, 5186; κ = 0.384 ± 0.004** (ours 0.3837 from y⁺ = 350, 0.3845 from 3√Re_τ; sublayer within 5e-4 for y⁺ < 1; Spalding and mixing-length channel "approximate") — cited, no licence statement; **redistributed as is by the user's decision of 2026-10-07** | Lee & Moser, J. Fluid Mech. 774, 395 (2015), doi:10.1017/jfm.2015.268 (arXiv:1410.7809); files at turbulence.oden.utexas.edu/channel2015/data | `reference/ch12/lee_moser_2015_channel_mean.csv` (281-row subset), `benchmarks.json`, `SOURCES.md`, `make_refs.py` | ch12 C10, C11, C12 (V5), E7, E8 |
| Prandtl's friction law 1/√f = 2.0 log₁₀(Re√f) − 0.8; Pope's c_L = 6.78; one-dimensional Kolmogorov constant 0.53 ± 0.055 (ours 0.491); k–ε standard constants; van Driest A⁺ = 26 — **all secondary** (the last three "constant consistency", not V5). **Not a benchmark**: the "Businger–Dyer" φ_m coefficients (unread first-hand) | Prandtl (1935) via McKeon et al., JFM 538, 429 (2005); Pope (2000) via arXiv:1705.04917 and the ATOMIX wiki; Sreenivasan, Phys. Fluids 7, 2778 (1995) via ATOMIX; Launder & Sharma (1974) via documentation pages | `reference/ch12/benchmarks.json`, `SOURCES.md` | ch12 C08, C11, C13, C15 |
| Our ch12 data (public, labelled ours): model-channel and van Driest tables (byte-reproducible) | `ch12.explainer_tables`, `write_reference_tables` | `reference/ch12/explainer_tables.json` | ch12 E5, E8 |
| **Eady problem: growth coefficient 0.3098, fastest wavenumber 1.606, cut-off (two digits)** (ours 0.30982, 1.60612, 2.39936) | K. A. Emanuel, MIT OpenCourseWare 12.803 (Fall 2009), Lecture 19 — read first-hand by the verifier, not re-read by the reviewer; NCL `eady_growth_rate` documentation as a second source | `reference/ch13/benchmarks.json`, `SOURCES.md`, `make_refs.py` | ch13 C16 (V5), E10 |
| Our ch13 data (public, labelled ours): cached runs `kelvin_basin`, `pv_particles`, `turbulence_f`, `turbulence_beta` (float16 `.npz`, regenerated within 4.1e-4 of range), explainer constants. **Not benchmarks**: Earth's rotation rate and radius (constant consistency); adjustment energy split, β-plane jet growth rates, −3 range, Rhines length (sources unread) | `scripts/ch13_make_caches.py`, `ch13.write_reference_runs` | `reference/ch13/*.npz`, `explainer_constants.json` | ch13 C11, C13, C17, E-parity rows |
Book-printed values (private, git-ignored): `tests/book_values_ch01.json` … `…_ch13.json` (ch13: rounded rotation rate,
typical f, β, N, H and radii, worked-example inputs, exercise answers; `_forbidden_public` strings and regex patterns).
**Open**: measured high-Re cylinder C_p (ch06), primary c_g,min (ch07), primary pipe transition (ch08), Glauert wall-jet
constants, Lienhard C_D, Roshko St, a primary turbulent-plate coefficient 0.074 (ch09), Schäfer–Turek DFG cylinder
benchmark (ch10), a counter-rotating Taylor number, primary Michalke digits (ch11), first-hand Businger–Dyer, primary k–ε,
Prandtl and van Driest papers (ch12), **first-hand Gill 1982, Kuo, Kraichnan 1967, Rhines 1975, Eady 1949 (ch13)**.

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
(ch12: CORE rows by primary label; ch13: CORE rows carrying each label. V7 tests: ch01 19, ch02 12, ch03 7, ch04 3, ch06 2,
ch07 13, ch08 6, ch09 24, ch10 20, ch11 41, ch12 10, ch13 15. Full suite at the ch13 gate: **1884 tests** (19 slow).)

## Open across chapters
- **ch13 open**: the list at the top of this file and `knowledge/ch13.md` §0, §9 (also: slips #5, #9, #13 not re-read by
  the verifier; two late equations of §13.16 not re-read by the reviewer; no explainer for C01, C07, C14, C17).
- **ch12 open** (`knowledge/ch12.md` §0, §9): `skin_friction_zpg` 0.074 unverified since ch09; secondary sources for the
  k–ε constants, Prandtl's law, van Driest's A⁺; core cap 0.09 uncited; `@slow` tests without mutants. **ch11 open**
  (`knowledge/ch11.md` §9): no citable Taylor benchmark for μ < 0; Bickley lower edges below k = 0.02; complex-path limits.
- **Library pass (orchestrator, when nothing reads `viz/` or `assets/`)**: ranked lists in `knowledge/viz_patterns.md` —
  **ch13**: `Viz.arrowPx` (sixth chapter), a y-title strip in `v.plot()`, log axes with decade ticks and a symmetric-log
  bar scale, pager "keep a heading with the next block" and "re-measure after the live fill", marker glyphs, `.c-amber`;
  **ch12**: `Viz.tpow`, `fitText`, precise `erf`, Gauss–Legendre, seeded draws, log-paced transport, two-line status,
  `hideOn`; **ch11**: complex numbers, cubic solver; **ch10**: `Viz.num.linalg`. Then `viz_inline.py --all` and shot quick.
- **Skill pass**: lesson candidates for the five skills in `knowledge/viz_patterns.md` (ch02–ch13) are not yet in the
  skills (ch13 first: figure notes from the rendered figure; no live `eq` card in a phone walkthrough step; project the
  state when de-aliasing; relative tolerances for small quantities; a solver-independent invariant for a disputed claim).
- **Machinery / tooling TODO**: **`shot.py`: 844×345, `.katex-error`, control characters, presets with Explain open**; a
  figure-note guard in the builder; word patterns for `check_public.py`; the chapter module in the scope of `py:` lambdas;
  code-line length lint ≈ 46 chars; `viz_lint` control-character rule; a contract-audit tool; per-agent scratch subfolders.
- **Explainer follow-ups**: per chapter in `knowledge/chNN.md` §9 and `viz_patterns.md`. Backups not built: ch04
  `kinematic_free_surface`, ch05 `vortex_rings`, ch06 `flow_net_sources_vortices`, ch07 `linearised_free_surface`, ch08
  `rotating_cylinders_couette`, ch09 `ball_swing_magnus`, ch10 `operator_splitting_theta`, ch11 `period_doubling_route`,
  ch12 `k_epsilon_decay`, ch13 `inertia_gravity_beams`.
