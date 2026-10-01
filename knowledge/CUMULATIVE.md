# CUMULATIVE knowledge — fluidpy
(Rewritten by the knowledge-keeper after every chapter. Every agent reads this first. Last rewrite: after ch11,
2026-10-01, knowledge pass on commit 4bc7060.)

## Project rules in force (added during ch01–ch03; sharpened since)
- **5–10 explainers per chapter, as many as the CORE ideas need** (`book.yaml → project.min/max_explainers_per_chapter`);
  ch03 has 7, ch04–ch09 9 each, ch10 8, ch11 9.
- **Every book equation is shown in full next to its number** — notebook prose, derivation steps, traps, recaps and
  explainer tour/Explain/Derivation/quiz/notes/status text **and `<meta>` strings**. Enforced by `tools/coverage_check.py`
  check 8 (markdown), `tools/eq_refs.py` (explainers; misses numbers inside template strings — ch10 E8) and, since ch10, the
  builder's own `self_check_near` (a bare number with no equation nearby); reviewers still scan by eye.
- **No control characters**: `coverage_check` check 9 for notebooks; ch10's builder adds `self_check_ctrl`. **Explainers are
  still not linted for it** (five chapters shipped one; machinery TODO).
- **Rule 9 (public repo) in practice**: a computed value can round to an exercise's printed answer — build it from the
  function and print 4 s.f.; the book's value lives only in `tests/book_values_chNN.json` (valid JSON). **Since ch10 the
  analyst writes book run parameters (Mach, grids, geometry, σ, periods) only there, under `"_forbidden_public"`, and
  `tools/check_public.py` fails if any tracked text file contains one of those strings** (ch10 leaked them into analysis,
  curation and a report; the unpushed commits were squashed). check_public skips an invalid JSON silently — keep it valid.
- **Pager and audit**: an item taller than a page is sliced at natural breaks; `tools/shot.py` pages through every page of
  every pager at every size. **The notebook (`build_chNN.py --dump`) is the reference for derivation step counts and
  titles** — re-check after the lesson review too (ch09 drift; ch10 matched 20/20).
- **Two values, the book's first** (ch09): when the book's number and the exact/computed one differ, show both everywhere.
- **Our own run parameters when the book's are private** (ch10: cavity Ma 0.08, block Ma 0.06, H = 4, cylinder W = 5d,
  σ = 0.8), and a **DEVIATION box** wherever the code departs from the book (formula, reason, both numbers).
- **Eigen-solver contract (ch11)**: return the **leading mode or raise**; **a 0 or NaN in a table has a stated meaning**
  (docstring + file header); unbounded profiles use the box max(profile box, 12/k) (`ST.decay_box`), shown converged by box
  **and** N doubling; every headline eigenvalue has an independent route (shooting, determinant, compound matrix, Galerkin).
- **`tools/check_public.py` runs after `git add`** (it sees tracked/staged files only — ch11); a chained `… | tail` hides a
  failing exit code.

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
ch11 Instability (done) — "basic state + small disturbance → linearise → normal modes → eigenvalue problem for the growth
     rate → neutral curve → its minimum is the critical number; the theorems (Ri > ¼, semicircle, inflection point) are
     one-way; linear onset is not transition; past it, three exact ODEs are unpredictable"
  normal mode (11.1), σ = −ikc, stable = for every k (D01) [C01]; Kelvin–Helmholtz (11.18), ΔU_min 6.70 m/s with depth and
     tension (D02–D04) [C02]; Bénard: Ra (11.21), amplitude problem (11.36)–(11.37), σ real (D05–D08) [C03]; rigid–rigid
     determinant → Ra_c 1707.76 at K_c 3.1163 (D09, D10 ★★★) [C04]; free–free (11.44), 27π⁴/4, growth quadratic (D11, D12)
     [C05]; salt fingers Rs − Ra = 27π⁴/4 (11.46) (D13) [C06]; Taylor–Couette: ring interchange, (11.51), Ta (11.52),
     Ta_c(μ = 0) 3389.9, Bénard limit (D14 ★★★, D15) [C07]; Taylor–Goldstein (11.61), neutral J = k(1 − k) (D16, D17) [C08];
     Miles–Howard Ri > ¼ ⇒ stable (11.67) (D18 ★★★) [C09]; Howard's semicircle (11.72) (D19 ★★★) [C10]; Squire (11.78),
     Orr–Sommerfeld (11.79) by Chebyshev collocation (D20, D21) [C11]; Rayleigh (11.84), Fjørtoft (11.86), critical layer,
     complex-path solver (D22) [C12]; Poiseuille Re_c 5772.22, Blasius 519.06 (δ*), Bickley 4.017 (two bands + gap) [C13];
     disturbance energy dE/dt = P − Λ (11.88) (D23) [C14]; Lorenz (11.91), r_H 24.737, period doubling (D24 ★★★, D25) [C15]

next: ch12 Turbulence — (11.88)'s −⟨uv⟩U′ is the turbulence production (`ST.disturbance_energy_budget`); transition after
      linear onset (`ch11.table_11_1`, TS waves `ts_mode`); `MAC.run`/`taylor_green` for 2-D decaying turbulence;
      `ch10.dominant_frequency` for spectra (move to core); `core.laminar` f = 64/Re, `core.boundary_layer` momentum
      integral and H, `core.jets` (turbulent jets); Lorenz/logistic tools for sensitivity; turbulent 0.074 needs a primary source.
```
The book ahead: Ch12 turbulence → **Ch13 GFD (uses 4, 5, 7–12)** → Ch14 aerodynamics · Ch15 compressible · Ch16 biofluids.
Where earlier chapters feed in: ch01 N², θ → Ch. 11, 13; isentropic gas → Ch. 15. ch02 Gauss → every CV; Stokes → Ch. 13 PV;
invariants → Ch. 12. ch03 vortices → Ch. 13, 14. ch04 NS residual tools → every exact solution; CV budgets → Ch. 12, 14, 15;
rotating frame, Boussinesq, Ri, Ro → **Ch. 13**. ch05 sheet roll-up → Ch. 11; **(5.30), (ζ + f)/h → Ch. 13**. ch06
`core.potential`, `conformal`, `panels` → **Ch. 14**; `laplace_solvers` → Ch. 13 PV inversion. ch07 `core.waves` → **Ch. 13**
(√(gH), Poincaré/Kelvin/Rossby, adjustment), **Ch. 11** (KH, RT), **Ch. 15** (jump ↔ shock). ch08 `core.laminar` → **Ch. 11**
base states, **Ch. 12** (f = 64/Re), **Ch. 13** (Stokes layer ≡ Ekman with ω → f, spin-up); `lubrication`, `creeping` → Ch. 13,
16. ch09 `core.boundary_layer` → **Ch. 11** (Blasius/FS base profiles), **Ch. 12** (momentum integral, H), **Ch. 14** (Thwaites
on airfoils); `bluff_body` → Ch. 11, 13, 14; `jets` → **Ch. 12**, Ch. 11 (Bickley), Ch. 13 (plumes); teacup → **Ch. 13 Ekman**.
**ch10** `core.fd` → **Ch. 11** (G(θ), FD operators), **Ch. 12** (numerical vs eddy diffusion), **Ch. 13** (CFL for
gravity/Kelvin waves, upwind tracer advection), Ch. 15; `core.fem1d` → Ch. 13 (Ekman ODE with K(z)), Ch. 11 (Galerkin
eigenproblems); **`core.mac` → Ch. 13 C-grid shallow water and Boussinesq pressure solve**, Ch. 11 saturation, Ch. 12 2-D
turbulence; `core.maccormack` → **Ch. 15** shocks; `core.fem2d` → Ch. 14 low-Re bodies, Ch. 11 wake base flows, Ch. 16.
**ch11** `core.stability` → **Ch. 13** (Rayleigh–Kuo = Rayleigh's equation with U″ → U″ − β via `rayleigh_eigs_contour`;
Eady/Charney by `cheb_grid` + `constrained_eig`; `taylor_goldstein_eigs`, Ri for mixing; inertial instability = Rayleigh's
circulation criterion with absolute angular momentum; **Squire fails with rotation or stratification**), **Ch. 12** (energy
budget → TKE production; TS waves → transition), Ch. 14 (laminar-flow airfoils), Ch. 15 (compressible vortex sheet).

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
Chapter-only code (move to `core/` on second use): **ch11** `os_box_numerics`, `os_leading_mode` (discrete mode vs
continuum), `parallel_profile`, `gradient_richardson`, `miles_howard_stable`, `rayleigh_number`, `benard_*`, `taylor_*`,
`salt_finger_regime`, `rayleigh_circulation_criterion` (→ **Ch. 13**), `lorenz_*`, `logistic_*`, `table_11_1` (→ Ch. 12); ch03 `pipe_profile`; ch04 `couette_heating`, `two_fluid_couette`,
`wake_drag_per_span`, `linear_wave_surface` (overflows for kH ≳ 710; use `ch07.wave_fields`); ch05 `rotating_cylinder_flow`,
`lamb_oseen_circulation`, `sheet_rollup` (→ Ch. 11), `column_relative_vorticity` (→ Ch. 13), `ring_dynamics`; ch06
`elliptic_cylinder_flow`, `lift_per_span` (→ Ch. 14), `sphere_motion`; ch07 `wave_fields`, `particle_path`/`stokes_drift`
(→ Ch. 13), `hydraulic_jump`, `simple_wave_evolve` (→ Ch. 15), `kdv_solve`, `boussinesq_linear_sympy` (→ Ch. 13); ch08 and
ch09 sympy engines; ch09 `secondary_flow_radial_force` (→ **Ch. 13**); **ch10** sympy engines, `split_linear_system`,
`marchuk_yanenko`, `theta_scheme_linear` (→ Ch. 13 splitting), `checkerboard`, `collocated_gradient`, `gradient_null_space`,
`artificial_compressibility_channel`, `ghia_centreline`, `hou_centres`, `cavity_error_vs_ghia`, `dominant_frequency`
(→ **Ch. 12** spectra: move to core there), `strouhal_from_period`, `cfl_time_step` (→ **Ch. 13**: move to core there),
`book_slips`. **Default g**: ch04, ch05, ch07, `core.waves` G_BOOK = 9.81; the rest (incl. ch11) G0 = 9.80665. **Default ρ**:
ch06 1.2 vs 1000; ch07, ch08 1000; ch09 BL 1.2, teacup 1000; ch10 non-dimensional (ρ = 1) except §10.2–10.3 (no ρ); ch11
non-dimensional per flow (KH dimensional, air over water in the examples).

### Explainer engine (`assets/viz_lib.js`)
`Viz.app` (tabs Walkthrough / Explore / Explain / Derivation / Equations / Code / Check; fit-to-window with density levels and
pagers; transport, modes, presets, status, terms, inspector, notes, selftest parity), `Plot`, `Viz.field`, `Viz.num` (RK4,
odeint, brentq, erf/erfc **~1.2e-7 only — a defect**, niceTicks), `Viz.work`, `Viz.three`, KaTeX with fallback; `Viz.font`,
`Viz.roundRect`, `Viz.text`, `Viz.card`, `Viz.fmtTime`, `Viz.tnum(keepTiny)`; pager slicing, `--viz-stage-need`,
`ev.viewId`/`ev.vizView`; a user-driven `set()` always rebuilds Explain (ch09 fix). **No linear algebra in the library**: ch10's
builders wrote dense LU, min-norm lstsq, banded LU, two tridiagonal solvers, Cholesky + Jacobi eigen, a background solver
pump — `Viz.num.linalg` is the top ch10 promotion candidate. Quirk (ch10): **Code-tab `{{placeholders}}` share one namespace
across code blocks**. **No JS promoted after ch02–ch11** (the site-publisher was reading `viz/` each time). ch11 changed only
`assets/viz_base.css` (orchestrator, re-inlined everywhere): coloured `<b class="c-…">` legends in view titles keep their
colour; an empty stage badge is hidden. ch11 builders wrote complex-pair formatters, `bisect`/`goldenMax`, rectangular-table
interpolation, a real cubic solver, a text-checksum parity helper and a canvas orthographic 3-D projector. Ranked lists in
`knowledge/viz_patterns.md`. Reference: `templates/viz_example.html`.

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

## Notation
See `knowledge/notation.md` (sign traps first, then the register, then per-chapter conventions). Everywhere: SI and kelvin
inside functions; **z up** (except ch06 §6.8); **lapse rate Kundu Γ ≡ dT/dz in code, meteorology shown**; `p0` ≠ `p_ref`;
pressures absolute unless `_gauge` (ch07 gauge); signed τ_xy = μ ∂u/∂y. ch02: passive C; traction on the first index; grid
`[k, j, i]`. ch03: R = G − Gᵀ, spin ½ω. ch04: Coriolis acceleration +2Ω × u′ vs force −2Ω × u′; **2-D ψ: u = ∂ψ/∂y** (ch07:
u = ∂ψ/∂z). ch05: γ = u_below − u_above; `ellipk(m)`, m = k². ch06: **Γ counterclockwise in code, `Gamma_cw=`/`Gamma_ccw=`**.
ch07: **ω = angular frequency**, η = elevation, g′ = g(ρ₂ − ρ₁)/ρ₂. ch08: **dp/dx book sign, `G=` = −dp/dx**; η = y/√(νt).
ch09: x along the wall; dp/dx > 0 adverse; η = y/δ(x) per flow; θ = momentum thickness; FS exponent book n = code `m`; body
angles from the forward stagnation point; λ = Thwaites parameter and street eigenvalue; St with the cyclic frequency.
**ch10 (changes):** three grids (node / element / staggered cell, named in every docstring); **θ = Fourier angle per cell**
(book e^{iπkx_i}, θ = kπΔx); **i = grid index and √−1** (upright i, index j in our lines); **α, β three ways** (FTCS, Θ-scheme
α_Θ, β_Θ, FE time α_t, β_t); **C = 2α** Courant; **Ma** for Mach (book M), **M** mass matrix; **K** stiffness vs K_e error
constant; **𝒮** trial space vs St; **R** global Péclet vs R_cell; **δ** layer thickness vs δ_AB; **F** FE load vs y-flux;
**G** amplification vs G_A; **p** pressure vs order p; **D** diffusivity vs **D**[u]; **g** Dirichlet value vs body force;
**λ** eigenvalue; **p′** Newton correction; printed slips "#1…#12" in the notebook (R1–R12 in code) vs recap IDs R01–R13.
**ch11 (changes):** **σ = σ_r + iσ_i = −ikc is the growth rate** (surface tension σ_s); **K** Bénard wavenumber vs k; **Γ three
ways** ((11.21) −dT̄/dz, Kundu dT/dz in code, meteorology shown) and Γ = circulation in §11.6; **Ra signed in §11.5**; `mu` =
Ω₂/Ω₁ (not viscosity); β haline; **Re per flow** (half-width, δ*, L); J bulk Richardson number; r, b, Pr of Lorenz; φ and ψ
three conventions each; N² is (7.127); slips S1–S13 = "#1…#13".

## Teaching lessons
- **Depth tiered, coverage exhaustive** (A / B / C, derivations, cells): ch01 15/70/21, 12, 496 · ch02 16/60/18, 15, 405
  · ch03 15/52/12, 24, 413 · ch04 15/151/19, 30, 611 · ch05 14/69/11, 23, 496 · ch06 15/116/29, 31, 564 · ch07 16/193/27,
  37 (321 steps), 617 · ch08 15/113/15, 33 (301 steps), 501 · ch09 14/118/16, 22 (224 steps), 472 · **ch10 15/105/23, 23
  (222 steps), 602** · **ch11 15/131/24, 25 (258 steps), 618**.
- **One-way theorems get a counter-example and the word "necessary" or "sufficient"** (ch11, P255): Ri < ¼ only allows
  instability (R = 3, J = 0.4 grows with Ri_centre 0.4); Rayleigh and Fjørtoft met and nothing grows (sin y, 2b < π);
  top-heavy yet stable below 27π⁴/4; linear onset ≠ transition. **"Trust what does not move with N"** (raw spectra at two
  resolutions overlaid); say what a table does not resolve next to its figure; a ★★★ check cell must test the derived
  result itself (ch11 D24 checked two steps of fifteen).
- **Convention callouts** ("which p_o?", Coriolis term vs force, Γ_cw vs Γ_ccw, which ω/η/g′/p′/dp/dx/δ/Re/θ/λ): state each
  convention, the size of the difference in numbers, which the code uses. **ch10: a conventions table up front for every
  overloaded letter, then enforced in every later cell** (the lesson review found (10.129)–(10.133) shown with bare α, β, θ
  after the table had renamed them).
- **Book slips are taught by computing both versions**: a slips table up front, printed forms as code options a test must
  fail; check a suspected slip on the page first. **ch10: name slips "#k" so they cannot collide with recap IDs "Rkk"**.
- **Pre-asymptotic honesty** (ch10): state the grids with every observed order (BTCS 1.80 → 1.98, upwind 0.76 → 0.947,
  Θ(¼) 0.86 → 0.995); teach a failed convergence study as failed (block C_D p = 0.43); a benchmark saturates at its own error
  (MAC 128² not closer to Ghia than 64²) — self-convergence is the convergence evidence.
- **Make an approximation measurable** (ch01, ch07, ch08, ch09 Thwaites; ch10 upwind vs its modified equation: gap 0.0134 →
  1.1e-4; numerical diffusivity 500 m²/s vs model eddy diffusivities); a result sensitive to its inlet is a **bracket**.
- **Every key number by two independent routes** (ch06 lift; ch07 FD eigenproblem; ch09 Blasius three ways; **ch10: every
  stability region by its closed form and a brute-force scan; Ghia by two transcriptions; Kovasznay by sympy and by orders**).
- **Every prose number is printed by a cell**, including numbers written during a review round; captions count what they
  claim; **drawings must draw what the reading note describes** (ch10 round 1: three drawings did not).
- **DEVIATION boxes** (ch10): when the code departs from the book (MacCormack additive Δt, corner and inflow densities) say
  so where the book's formula is taught, with the reason and both numbers computed (2.79e-4 vs 5.62e-4).
- **A tool first used inside a derivation gets its primer first** (ch09 P218a); **so does a named test case** (ch10:
  Taylor–Green P254 was used before being explained; ch11: P263, P269 placed after first use). `knowledge/primers.md` lists
  280 (P01–P279 + P218a); **ch12 starts at P280**.
- **A sign-checking tiny example** for every transcribed equation whose sign matters (ch09 (9.8)).
- **Builder rendering hygiene**: design TeX converted, not copied; raw strings for LaTeX; regex tidies never inside maths;
  headings keep their equations; **ch10: `self_check_ctrl` + `self_check_near` in the builder** (a heredoc turned `\to` into
  TAB + "o"); templated code comments ("# show the numbers computed above") still appear — fix in the builder template.
- **Climate hooks with numbers** (ch04 Ro; ch05 tornado; ch06 ψ inversion; ch07 swell, tsunami, √(g′H)/f; ch08 Ekman look-alike;
  ch09 teacup, island wakes; **ch10 CFL steps 126 s ocean / 78 s atmosphere at 25 km, the Arakawa C-grid, upwind numerical
  diffusivity 500 m²/s at 10 km vs tens–100 m²/s used, dynamics/physics splitting, "reduced speed of sound"**).
- Climate hooks, ch11: billows and Ri < ¼ mixing, salt fingers, cloud streets, inertial instability, Rayleigh–Kuo, Lorenz.
- Reusable derivation moves (ch01–ch11, 190 rows; 14 from ch11) are tabulated in `knowledge/concept_map.md`.

## Global pitfalls confirmed in this project
**Machinery and environment**
- Extracted text garbles maths (`¼` = `=`, `ð…Þ` = parentheses, missing minus, ω → `u`, γ → g, θ → q, σ/τ → s, ν → n,
  ε → 3, ρ → r, Ω → U, ψ → j, φ → 4, η → h, ζ → z, ∂ → v, ∇ → V): read page images.
- **Windows: heredocs and `sed` mangle backslashes, quotes and `|`** — write code, explainers and knowledge files with
  Write/Edit or Python files; builder strings with LaTeX must be raw; single-backslash TeX in JS strings turns `\f`, `\t`,
  `\r`, `\n` into control characters (ch04, ch06, ch09 ×2, **ch10 `\to` in the notebook and `\approx` in an explainer**).
- `MPLBACKEND=Agg` strips inline figures; every chapter script takes `--no-show`; Anaconda's kernelspec can shadow
  `fluidpy-venv`; JS `UIEvent.view` is read-only.
- `tools/shot.py` does not click/drag (cannot see stale Explain text after a slider drag); its `py:` expressions have
  restricted builtins — **no ndarray methods (`.min()`, `.sum()`): use `np.min`, `np.add.reduce`** (ch10); set parity
  tolerances near the achieved agreement.
- **Intermittent pager clips (1–26 px) under machine load**: leave margin, run shot 3×; **clips can be state-dependent** —
  the audit reaches a tab with whatever preset the last step left (ch10: Equations at β = 100); KaTeX-late clips need one line
  of slack and text-size integrals in live rows.
- **Run the merge gate serially**: the full pytest (≈ 10 min alone at 1460 tests) ran past an hour in parallel with notebook
  execution and shot audits; `tests/test_machinery.py` separately (19 s) avoids the 96 % stall.
- **Parallel agents need their own scratchpad subfolders** (ch10 collision; **ch11 again — builders overwrote each other's
  scripts**). Do not edit a solver module while a shot run evaluates parity rows (ch11: `AttributeError` crash).
- **`tools/shot.py` does not fail on `.katex-error` nodes or control characters in rendered text** (ch11 E2 shipped two KaTeX
  errors); reviewers sweep the DOM by hand. **API overload interruptions** (ch11): short steps, consistent disk state.
- matplotlib 3.11: animations using `subplots_adjust`/`fig.text` disable constrained layout (ch10 cell 159 warning with a
  local path in the output — drop `subplots_adjust`).
- `coverage_check` matches ledger reminders by exact primer name; `eq_refs` misses bare numbers inside template strings and
  rejects en-dash ranges on `ref:` badges.
- Book wording can hide in docstrings and exercise answers in computed numbers; **book run parameters hide in analysis and
  reports** (ch10) — `_forbidden_public` + check_public.
- **Design and analysis documents contain errors; downstream agents compute, not copy** (every chapter; ch10: default-grid
  order expectations, a design stretched-grid map clustering at the wrong end, "stretched grid removes the wiggles",
  Θ(¼) order on the default list, "R3 forward wiggles for every R_cell"; ch11: g = 9.81 vs G0, Ra(6) 2672 vs 2680.85,
  Poiseuille band 0.80–1.07 vs 0.797–1.095, a root for τ = 0.01 quoted at τ = 0.0107, U″ for U′ in a signature, wrong quiz claims).

**Mathematics → code**
- **Unbounded profiles (ch11): a fixed truncation box is a different problem at long wavelength** — scale it with 1/k
  (`ST.decay_box`) and prove it by box **and** N doubling: TG false growth above the neutral curve at k = 0.05; Bickley lower
  branch 0.0676 → 0.0211 at Re = 7.2; Blasius Re_c 519.0765 → 519.0601, lower branch +9 % at Re = 6000.
- **A convergence filter can make a false zero or return a lower mode (ch11)**: TG tongue at k ≥ 0.65, Rayleigh jet at
  k = 1.8, 1.9, Bénard −39.6 instead of +7.12 (the finer solve was the worse one). Leading mode or raise.
- **Near a neutral point more real-axis points do not help (ch11)**: complex path round the critical layer
  (`rayleigh_eigs_contour`), nodes where the wave is (`decay_map_scale`), an independent route always. Spectral convergence
  is digits, not an order: report the error ladder and the round-off rise at large N.
- **Mode families (ch11)**: least-damped ≠ the mode to interpolate; neutral curves can have several bands (Bickley gap);
  track by continuation, never join across NaN; tell discrete modes from the continuum (`far_field_fraction`).
- **Classifiers and verdict strings need tests on both sides of every boundary and at degenerate inputs (ch11)**: top-heavy
  below the (11.46) margin, Ri at U′ = 0, KH with ΔU = 0 — 128 tests passed with both bugs. **A pinned number can pin a
  defect** (299 positives, 0.011389): pin sets and physics, re-derive when the numerics change.
- **A test comparing two of our own functions is not evidence** — and **a table written by the same function is not a
  reference** (ch10 ψ_min: the tests compared live runs with CSVs made by the same buggy fit). Pin conventions and fits to
  fields with known answers (a paraboloid), assert must-hold inequalities (fitted min ≤ grid min), and **prove discrimination
  by planting the wrong variant** (ch02 8/8 … ch09 25, **ch10 38/38**).
- **A test on |x| cannot see the sign of x** (ch09 (9.8)); every signed transcription gets a signed test on a manufactured
  field; tests exercise every parameter at non-zero values with both signs; **both signs of u for upwind (ch10 R11)**.
- **Series stop on the term envelope, never on the actual term** (ch10 (10.199): a term vanishes where the sine does);
  test against an independent solution at those points.
- **Defaults encode physics and stability**: ch09 `cp_base=None` gave 2.59; **ch10 a default Δt that dropped the viscous rate
  gave NaN — explicit solvers take an explicit rule, compute the full limit and raise** (`check_stability=True`, growth cap
  1e6, never NaN plots); scheme-dependent quantities take required arguments (`numerical_diffusivity(C=)`).
- **Build operators link by link** (`+=` over real neighbours), and test null vectors on 1-cell and 2-cell grids (ch10 F1).
- **Floating-point references can be the error** (ch10 F2: `(0.2 − 1.0) % 1`); shift by integers, round before edges.
- **Guard every division's degenerate case** (ch10 GCI with p = 0 → NaN).
- **A verifier's correction is a claim too**; **a mutant that hangs is not killed** — cap solver iterations (ch10: none hung).
- **Folds and bifurcations are mathematics** (ch09 Falkner–Skan): parametrise by the quantity that passes through.
- **BVPs on infinite domains**: scaling symmetry (Töpfer) or `solve_bvp` with continuation; η_max per member.
- **Complex-step derivatives need analytic continuations** (ch07); keep |c| and signed dω/dk apart.
- **Code the corrected form; keep the printed one as a named option a test must fail** (every chapter since ch07).
- **Never type a constant a function can compute** (ch01, ch09); keep the book's value as a separate named constant.
- **Overflow and cancellation**: cosh/sinh ratios; expm1; erfc; small-α series; **ch10 scaled (e^{Rx/L} − 1)/(e^R − 1) beyond
  R ≈ 700**; 0/0 starts with a finite limit (Thwaites at a stagnation point).
- **Stiff and degenerate PDEs**: CN after a jump needs BE start-up; parabolic marching with a √ψ wall singularity needs a
  mapping and BDF2 (ch09).
- **Observed orders only in the asymptotic range**; singular corners and staircase boundaries set the global order (ch06 4/3,
  ch08 1.11, ch09 0.98); **heuristic boundary closures can wreck it (ch10 block p = 0.43)**; a benchmark saturates at its own
  accuracy.
- **Weak compressibility on a fixed grid degrades as Ma falls** (Δt ∝ Ma·Δx); non-conservative boundary closures leak mass
  (ch10 MacCormack cavity −0.5 %) — label it.
- **Finite-difference tolerances from a written budget**, never loosened by feel.
- **A telescoping sum is not a conservation test**; genuine V4 invariants so far: ch07 ∫η², KdV; ch08 thin-film volume, ∫ω dy;
  ch09 J, Ψ, ρU²θ = ∫τ₀dx; **ch10 projection ∇·u ≤ 1e-12, periodic MacCormack mass (exact), FTCS periodic sum, ∫T of the
  advected Gaussian**.
- **Sibling functions share defaults and argument meanings** or are called by keyword; convention-named keywords (`Gamma_cw=`,
  `ref=`, `G=`, `sides=`, `closure=`, `branch=`, `printed=`, **`convention=`, `dt_rule=`, `scheme=`, `method=`, `lid=None`**).
- **Two quantities with one name** — the code names the reference in the keyword or the function name (ch10: θ, α, β, M, K,
  S, R, δ, F, G, p, D, g, λ; see notation).
- **Benchmarks from the primary source at its printed precision**, cross-checked against a second transcription (ch10 Ghia).
- **Validation labels**: V4 only for real conservation laws; "converged" only with an asserted order (ch10 review removed it
  from six functions); "benchmark" only inside the stated tolerance (ch10 MacCormack 64² 1.38 % is not); "approximate" for
  few-percent agreement; "secondary" for compilations; "qualitative" for confined/heuristic studies.
- **Heavy cached runs are invisible to mutants**; keep one slow test that recomputes a coarse cached run to 6 s.f.

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
Book-printed values (private, git-ignored): `tests/book_values_ch01.json` … `…_ch10.json` (ch10: §10.2 numbers, §10.5 run
parameters and results, Θ printed digits, periods, St, exercise inputs; `_forbidden_public` strings). **Open**: measured
high-Re cylinder C_p (ch06), primary c_g,min (ch07), primary pipe transition (ch08), Glauert wall-jet constants, Lienhard
C_D, Roshko St, Howarth angles, a primary turbulent-plate coefficient (ch09), Schäfer–Turek DFG cylinder benchmark (would
make the FE cylinder V5), Dennis & Chang Re = 40 primary (ch10), **a citable counter-rotating Taylor number (μ < 0;
Chandrasekhar 1961 on paper), primary Michalke digits and the odd Bénard mode page (ch11)**. `tests/book_values_ch11.json`
holds the book's rounded critical numbers, Tollmien coefficients, figure parameters and exercise inputs.

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
(V7 tests: ch01 19, ch02 12, ch03 7, ch04 3, ch06 2, ch07 13, ch08 6, ch09 24,
ch10 20, ch11 41. Full suite at the ch11 gate: **1620 tests**, 1594 s serially.)

## Open across chapters
- **ch11 open** (`knowledge/ch11.md` §9): no citable Taylor benchmark for μ < 0; Bickley lower edges below k = 0.02;
  `benard_growth_rate` explicit N ≥ 60 round-off; complex-path limits; notebook runtime 126–239 s (animations); E4 EQTAG
  workaround; E8 clipped axis title at 360×640; `tools/viz_tables/ch11/` absolute paths and the manual E8 splice.
- **Library pass (orchestrator, when nothing reads `viz/`)**: ranked lists in `knowledge/viz_patterns.md` — **ch11**: complex
  numbers (`Viz.cx` + pair formatting; asked since ch09), `bisect`/`goldenMax`, rectangular-table interpolation, cubic solver,
  `arrowPx`, checksum-row helper, canvas 3-D projector; **ch10**: `Viz.num.linalg`, log axes + `gutterPlot`, `polyFill`, matrix
  drawer, `pyRound`, `Run`, background pump, `views[i].hideOn`; ch09 (lint rule, formatters, `Viz.profiles`, hatch); ch08
  (pager-label CSS, erfc precision); ch07 (`Viz.waves`, per-mode hiding before Ch. 13). Then `tools/viz_inline.py --all` and
  `tools/shot.py --chapter ch01…ch11 --quick` + `templates/viz_example.html --quick`.
- **Skill pass**: lesson candidates for `interactive-viz`, `math-to-python` §7, `verify-implementation`, `teaching-style` and
  `colab-notebook` in `knowledge/viz_patterns.md` (ch02–ch11) are not yet in the skills (ch11 first: box ∝ 1/k with box and N
  doubling; leading mode or raise; a 0 has a meaning; independent route; both sides of every classifier boundary; pin sets
  not counts; builders recompute design numbers; one-way theorems with counter-examples).
- **Machinery / tooling TODO**: **`shot.py` failing on `.katex-error` and control characters (L3)**, hash deep links applied
  after `goStep(0)` (L4), code-line length lint ≈ 44 chars (L5), `eq_refs` allow-marker for mirrored fluidpy strings (L6);
  `viz_lint` control-character rule (**six chapters**); `eq_refs` en-dash ranges and template strings; a slider-drag Explain
  test; `check_public.py` on untracked files and **failing on an invalid `book_values_*.json`**; Derivation-tab vs `--dump`
  diff in `embed_check`; `shot.py` exposing `np` to `py:` rows, resetting presets per tab, a 3-run flake mode; nbkit —
  templated/echo-comment lint (three chapters), `self_check_ctrl`/`self_check_near` in `tools/nbkit.py`.
- **Explainer follow-ups**: ch10 E1–E7 items in `knowledge/ch10.md` §9; earlier chapters in their `chNN.md` §9. Backups not
  built: ch04 `kinematic_free_surface`, ch05 `vortex_rings`, ch06 `flow_net_sources_vortices`, ch07
  `linearised_free_surface`, ch08 `rotating_cylinders_couette`, ch09 `ball_swing_magnus`, ch10 `operator_splitting_theta`,
  ch11 `period_doubling_route`.
- Verification notes: ch10 block C_D non-asymptotic (Δx = 1/64 next), MacCormack cavity mass and low-Mach degradation, no
  Re = 400 MacCormack 128², cached-run mutants invisible; ch09 sphere Re_cr 5e5 vs 3e5 (G9), `karman_street_spectrum` σ → λ,
  turbulent 0.074 needs a primary source before Ch. 12; ch08 O5, O8, O9; ch07 O5–O7; ch06 O4, O6; ch05 O2.
