# CUMULATIVE knowledge — fluidpy
(Rewritten by the knowledge-keeper after every chapter. Every agent reads this first. Last rewrite: after ch08,
2026-09-29, knowledge pass on commit 61ba691.)

## Project rules in force (added during ch01–ch03; sharpened since)
- **5–10 explainers per chapter, as many as the CORE ideas need** (`book.yaml → project.min/max_explainers_per_chapter`);
  ch03 has 7, ch04–ch08 9 each.
- **Every book equation is shown in full next to its number** — notebook prose, derivation steps, traps, recaps and
  explainer tour/Explain/Derivation/quiz/notes/status text **and `<meta>` strings** (ch04 E6, ch06 E6). Enforced by
  `tools/coverage_check.py` check 8 (markdown) and `tools/eq_refs.py` (explainers); reviewers still scan by eye.
- **`coverage_check` check 9**: no control characters in cells (a LaTeX backslash eaten by a non-raw string).
- **Rule 9 (public repo) in practice (ch07)**: a value we compute can round to an exercise's printed answer. Numbers
  that coincide with a book value are built from the function (never typed) and printed with **4 significant figures**;
  the book's value lives only in `tests/book_values_chNN.json`.
- **Pager and audit**: an item taller than a page is sliced at natural breaks; `tools/shot.py` pages through every page
  of every pager at every size. **The notebook is the reference for derivation step counts and titles** (ch08): explainer
  Derivation tabs are checked against `notebooks/build_chNN.py --dump`, not against design Part F.

## Physics pipeline so far
```
ch01 Introduction — molecules → continuum ρ, p, T, u (D35, Kn) [C06]; transport τ = μ du/dy, ν = μ/ρ, FTCS [C12];
     statics dp/dz = −ρg (D05), buoyancy (D37), USSA-1976 [C20]; thermodynamics first law → Gibbs (D10) → c² → p = ρRT
     (D11) → p/ρ^γ (D14) [C25–C45]; stratification ζ″ + N²ζ = 0 (D18), N², Γ_a = −gαT/C_p (D19), θ (D20, D36)
     [C50–C55] (Kundu Γ ≡ dT/dz in code, meteorology −dT/dz shown); dimensional analysis → Π = null space (D28) [C64–C69]
ch02 Cartesian tensors — summation [C01] → C_ij, CCᵀ = I (D02, D03), x' = Cᵀx (D01) [C02, C03]; stress, Cauchy
     f_i = τ_ji n_j (D05) [C04, C05] → τ' = CᵀτC (D06), invariants (D18) [C06, C07]; ε–δ (D09) [C08]; S + A, R ↔ ω (D14,
     D15) [C12]; principal axes (D17 ★★★) [C13]; ∇φ, ∇·u, ∇×u on grids [C09–C11]; Gauss (D25) [C14, C15]; Stokes (D26) [C16]
ch03 Kinematics — particle ⇄ field (D01), DF/Dt = ∂F/∂t + u·∇F (3.5) (D02) [C01, C02]; stream/path/streak lines (D03–D05)
     [C03, C04]; Galilean (3.9) (D06) [C05]; du = G·dx → S, ½ω (D07–D15) [C06–C11]; principal axes (D16) [C12]; Rankine,
     Gaussian vortices (D17–D20) [C13, C14]; Leibniz (D21) → RTT (3.35) (D22 ★★★) [C15]
ch04 Conservation laws (the trunk) — mass (D01–D03) [C01, C02] → ψ (D04) [C03]; momentum (D05), Bernoulli (4.19) (D06),
     Cauchy (4.24) (D07) [C04–C06]; τ = −pδ + 2μS + λS_mmδ (D08, D09 ★★★) [C07] → NS (D11–D13) [C08]; rotating frame
     +2Ω×u′ (4.43) / −2Ω×u′ (4.45) (D14–D18) [C09]; energy, ε ≥ 0 (D19–D23) [C10]; Lamb, unsteady Bernoulli (4.75)
     (D24–D26) [C11, C12]; Boussinesq (D27, D28) [C13]; Dη/Dt = 0, jumps, Laplace (D29) [C14]; dimensionless NS (D30) [C15]
ch05 Vorticity — tubes (D01) [C01]; vortex pressure (D02, D03) [C02]; Kelvin (D04, D05) [C03]; baroclinic torque (D06,
     D07) [C04]; Helmholtz (D08) [C05]; Dω/Dt (5.13) (D09 ★★★) [C06]; Biot–Savart (D10, D11 ★★★), filaments (D12, D13)
     [C07, C08]; rotating baroclinic (5.30) (D14, D15 ★★★) [C09]; stretching/tilting (D16–D18) [C10]; (ζ + f)/h (D19,
     D20) [C11] ← PV seed; point vortices, images, sheets γ = u₂ − u₁ (D21–D23) [C12–C14]
ch06 Ideal flow — Laplace; solutions add; any streamline can be a wall (6.1) (D01) [C01]; ω_z = −∇²ψ, φ ⟂ ψ (D02–D04)
     [C02]; superposition (D05) [C03]; doublet, half-body, cylinder C_p = 1 − 4 sin²θ, D = 0, L = ρUΓ, images (D06–D13)
     [C04–C08]; w = φ + iψ, corners Azⁿ (D14, D15) [C09]; Blasius (6.60), Kutta–Zhukhovsky (D16–D18 ★★★) [C10]; conformal
     maps, outside root (D19–D21 ★★★) [C11]; Laplace on a grid, order 4/3 at a 270° corner (D22, D23) [C12]; Stokes ψ,
     sphere (D24, D25) [C13]; airship, axial method (D26, D27) [C14]; added mass ½ displaced mass (D28–D31 ★★★) [C15]
ch07 Gravity waves — "plane wave in, boundary conditions → algebra, ω(k) out; read everything from ω(k)"
     c = ω/k, Doppler (D01) [C01]; linearised free surface, relative error ka (D02–D04) [C02] → ω² = gk tanh kH (7.28)
     (D05, D06) [C03] → deep/shallow, pressure (D07–D09) [C04]; orbits (D10, D11) [C05]; E = ½ρga², flux (D12–D14) [C06];
     capillary c_min (D15, D16) [C07]; seiches (D17, D18) [C08]; c_g, F = E c_g (D19–D21, D20 ★★★) [C09]; rays, Snell
     (D22–D24) [C10]; jump, KdV (D25, D26) [C11]; Stokes drift (D27) [C12]; interface ε√(gk) (D28) [C13]; barotropic/
     baroclinic √(g′H), g′ = gΔρ/ρ₂ (D29–D31 ★★★) [C14]; w-equation, ω = N cos θ (D32–D34 ★★★) [C15]; c ⟂ c_g, F = c_gE
     (D35–D37) [C16]
ch08 Laminar flow (done) — "symmetry kills u·∇u → linear ODE + wall conditions; else scale until one balance
     dominates; no imposed length or time → a similarity variable"
  ν = μ/ρ a diffusivity, L²/ν (D01) [C01]; parallel-flow reduction v ≡ 0, "f(x) = g(y) ⇒ constant" (D02, D03) →
     Couette–Poiseuille (8.5), backflow dp/dx > 2μU/h², extremum y* = h/2 − μU/(h dp/dx) (D04, D05) [C02]; pipe (8.6),
     τ₀ = (a/2)dp/dz as a CV balance, Q ∝ a⁴, f = 64/Re (D06, D07) [C03]; circular Couette AR + B/R (8.9)–(8.12) (D08, D09) [C04]
  lubrication: two-length scaling, inertia × ε²Re_L (D10 ★★★), 0 ≅ −p_x/ρ + νu_yy (8.17a) (D11) [C05] → local Couette +
     Poiseuille (8.19), Reynolds equation h_t + q_x = 0 (D12, D13 ★★★) [C06] → slider hump and load, optimum 1 + α = 2.189,
     recirculation iff gap ratio > 2 (D14 ★★★) [C07] → thin film h_t = (ρg/3μ)(h³h_x)_x (D15) [C08]
  similarity: u_t = νu_yy, η = y/√(νt), erfc, δ₉₉ = 3.643√(νt) (D16–D20) [C09]; ansatz At^{−n}F(ξ/δ) + exponent matching:
     Ex. 8.4, vortex sheet n = ½ (D22 ★★★), bead n = m = ⅕ (D21–D23) [C10]; Stokes layer e^{−y/δ}cos(ωt − y/δ),
     δ = √(2ν/ω) (D24, D25) [C11]
  creeping: viscous pressure scale → ∇p = μ∇²u (D26) [C12]; ∇²ω = 0, (E²)²ψ = 0, Stokes ψ (8.48) (D27–D29, D28 ★★★)
     [C13]; p (8.50), D = 6πμaU = ⅓ pressure + ⅔ friction (D30, D31 ★★★) [C14]; inertia/viscous → ½Re_a r/a, Oseen
     C_D = (24/Re)(1 + 3Re/16) (D32, D33) [C15]

next: ch09 boundary layers — the (8.14) two-length scaling is the boundary-layer approximation; √(νt) ↔ √(νx/U);
      similarity engine (`similarity_reduce_sympy`, `similarity_collapse_error`) for Blasius and Falkner–Skan; ch06
      `core.potential` outer flows U_e(x); `crank_nicolson_1d` for parabolic marching; `core.laminar` for code tests.
```
The book ahead: Ch9 boundary layers → Ch10 CFD → Ch11 instability → Ch12 turbulence → **Ch13 GFD (uses 4, 5, 7, 8, 11,
12)** → Ch14 aerodynamics · Ch15 compressible · Ch16 biofluids.
Where earlier chapters feed in: ch01 N², θ → Ch. 11, 13; isentropic gas → Ch. 15; FTCS → Ch. 10. ch02 Gauss → every CV;
Stokes → Ch. 13 PV; invariants → Ch. 12. ch03 vortices → Ch. 13, 14. ch04 **NS residual/term tools → every exact
solution**; CV budgets → Ch. 9, 14, 15; ε → Ch. 12; rotating frame, Boussinesq, Ri, Ro → **Ch. 13**. ch05 vorticity
budget → Ch. 10, 12; sheet roll-up → Ch. 11; **(5.30), (ζ + f)/h → Ch. 13**; filaments → Ch. 14. ch06 `core.potential`
→ **Ch. 9** outer flows (U_e(x), wedge Azⁿ → Falkner–Skan), **Ch. 14**; `core.conformal`, `core.panels` → Ch. 14;
`core.laplace_solvers` → **Ch. 10**, Ch. 13 PV inversion. ch07 `core.waves` → **Ch. 13** (√(gH), Poincaré/Kelvin/Rossby
c_g by complex step, √(g′H)/f, inertia–gravity waves, WKB rays, geostrophic adjustment), **Ch. 11** (KH, σk²/ρ, RT flag),
**Ch. 15** (jump ↔ shock), Ch. 10 (`kdv_solve`). **ch08** `core.laminar` → **Ch. 11** (Couette/Poiseuille/Taylor–Couette
base states, Rayleigh in A, B), **Ch. 12** (f = 64/Re, τ₀ → u_τ, linear total stress), **Ch. 9** (√(νt), temporal layer),
**Ch. 13** (Stokes layer ≡ Ekman layer with ω → f, spin-up L²/ν, cyclostrophic balance), Ch. 10 (exact test solutions);
`core.lubrication` → **Ch. 9** (thin-layer scaling), **Ch. 13** (hydrostatic thin layers, gravity currents), Ch. 16;
`core.creeping` → **Ch. 13** (droplet/sediment settling), Ch. 16 (micro-swimmers); `crank_nicolson_1d` → **Ch. 10**.

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
| `tools/convergence.observed_order(h, err)` | log–log slope | every V3 test |

### Physics primitives (re-exported as `chNN.<name>`; details in `knowledge/chNN.md` §3)
| function family | module | book § / Eq. | label |
|---|---|---|---|
| constants `K_B, N_A, R_U, R_AIR, GAMMA_AIR, CP_AIR, G0, P_ATM, P_REF`, `G_BOOK` = 9.81; perfect gas, first law, entropy, sound speed, isentropic, van der Waals | `thermo` | §1.8–1.9, (4.94) | analytic, symbolic, converged, conserved, benchmark |
| hydrostatics, USSA-1976, buoyancy; N² family, `classify_stability`, `parcel_displacement`, lapse-rate conventions, θ; FTCS diffusion; Π groups | `statics`, `stratification`, `diffusion`, `dimensional` | ch01 | analytic, symbolic, converged, benchmark |
| index notation, `tensors`, grids `[k, j, i]` + `operators`, `fields`, `integral_theorems` | ch02 modules | (2.1)–(2.36) | analytic, symbolic, converged, conserved, benchmark |
| `kinematics`, `coords`, `transport` (Leibniz, RTT); `vortices` | ch03/ch05 modules | (3.1)–(3.35), §5.4 | analytic, symbolic, converged, conserved, benchmark |
| `conservation` budgets on any `ControlVolume`; `constitutive`; `navier_stokes` residuals and term splits (8 exact solutions); `curvilinear` (cylindrical/spherical grad, div, curl, vector Laplacian, strain rate); `streamfunction`; `rotating`; `bernoulli`; `interfaces`; `similarity` (`Scales`, 20+ groups, `nondimensional_ns_coefficients(pressure_scale=)`, `sphere_drag_coefficient`) | ch04 modules | (4.5)–(4.119) | analytic, symbolic, converged, conserved, benchmark |
| `vorticity`, `biot_savart` | ch05 modules | (5.3)–(5.33) | analytic, symbolic, converged, conserved |
| `potential`, `conformal`, `laplace_solvers`, `panels` | ch06 modules | (6.1)–(6.109) | analytic, symbolic, converged, conserved, benchmark |
| **`waves`**: vocabulary; dispersion (`omega_gravity`, `phase_speed` \|c\|, `group_velocity` **signed**, inverse by brentq, `depth_regime`, capillary minima); any-ω (`group_velocity_numeric` complex step, `group_velocity_vector`, packets, `linear_evolve(_2d)`, `ray_trace`); layers and stratification (`interface_omega`, `two_layer_*`, `reduced_gravity_book(ref=)`, `internal_wave_omega`, `beam_angle`, `internal_wave_velocities`); overflow-safe `cosh_over_sinh` family | `waves` (ch07) | (7.1)–(7.159) | analytic, symbolic, converged, conserved, benchmark |
| **`laminar`** parallel flows: `channel_flow(y, h, U, dpdx, mu, G=)`, `_rate`, `_shear_stress`, `channel_backflow_threshold`, `couette_poiseuille_state` (extremum y*); `pipe_poiseuille`, `pipe_wall_stress`, `pipe_flow_rate`, `pipe_friction_factor`; `circular_couette` (exact R₂ = ∞, R₁ = 0), `_pressure`, `_shear_stress`, `_power`, `_state` (Rayleigh) | `laminar` (ch08) | (8.4)–(8.12) | analytic, symbolic, conserved |
| **`laminar`** unsteady: `similarity_variable(half=)`, `stokes_first_problem` (ch04 re-exports it), `stokes_first_vorticity`, `diffusion_thickness(level=)`, `stokes_first_stopped`, `vortex_sheet_diffusion(U)`, `temporal_bl_wall_stress`, `line_vortex_decay`/`spinup`; `stokes_second_problem`, `stokes_layer` (δ_e **and** δ_book), `stokes_layer_state`; `*_state` dicts | `laminar` (ch08) | (8.20)–(8.38) | analytic, symbolic, converged, conserved |
| **`lubrication`**: `lubrication_scales` (ε, Re_L, ε²Re_L, Λ, μUL/h²), `lubrication_term_magnitudes`, `lubrication_velocity(form=)`, `lubrication_flux`, `reynolds_pressure_1d` (any h(x)); `slider_bearing(model=)`, `_load`, `slider_optimum_taper`, `slider_bearing_state` (wide-end recirculation), `slider_gap_velocity(frame=)`; `hele_shaw_*`; `thin_film_flux`, **`thin_film_spread`** (implicit Newton, precursor, volume 2e-16), `viscous_current_similarity`, `viscous_current_eta_N` | `lubrication` (ch08) | (8.13)–(8.19), Ex. 8.1–8.3, 8.7 | analytic, symbolic, converged, conserved, benchmark |
| **`creeping`**: `stokes_residual`, `E2`/`E4_residual`/`stokes_sphere_sympy`, `stokes_sphere_streamfunction/velocity(frame=)`, `stokes_sphere_pressure`, `_surface_stresses`, `stokes_drag(parts=)`, `stokes_drag_running`, `sphere_drag_quadrature`, `side_line_speed(model=)`; drag laws (Stokes, Oseen diameter-Re, Proudman–Pearson); `terminal_velocity` (g = G0, warns), `settling_state`, `millikan_charge`; `inertia_viscous_ratio`, `oseen_streamfunction`, `oseen_velocity` | `creeping` (ch08) | (8.39)–(8.53) | analytic, symbolic, converged, benchmark |
| **`diffusion.crank_nicolson_1d`**(u0, y, dt, nsteps, D, bc_left, bc_right, startup_be=2, theta=0.5, …) — θ-scheme, two BE start-up steps, order 2 | `diffusion` (ch08 addition) | (8.20), (8.33) (our scheme) | converged |
Chapter-only code (move to `core/` on second use): ch03 `pipe_profile`; ch04 `couette_heating`, `two_fluid_couette`,
`wake_drag_per_span` (→ Ch. 9), `linear_wave_surface` (**overflows for kH ≳ 710; use `ch07.wave_fields`**); ch05
`rotating_cylinder_flow` (ω = 2Ω), `lamb_oseen_circulation`, `sheet_rollup` (→ Ch. 11), `column_relative_vorticity`
(→ Ch. 13), `ring_dynamics` (→ Ch. 14); ch06 `elliptic_cylinder_flow`, `lift_per_span` (→ Ch. 14), `half_body_surface_cp`
(→ Ch. 9), `sphere_motion`; ch07 `wave_fields`, `free_surface_residuals`, `particle_path`/`stokes_drift` (→ Ch. 13),
`hydraulic_jump`, `simple_wave_evolve` (→ Ch. 15), `kdv_solve` (→ Ch. 10), `seiche_modes`, `internal_wave_fields`,
`boussinesq_linear_sympy` (→ Ch. 13); **ch08** sympy engines (`parallel_flow_sympy`, `lubrication_nondim_sympy`,
`reynolds_equation_sympy`, `slider_bearing_sympy`, `stokes_second_sympy`, `low_re_scaling_sympy`, `oseen_*_sympy`),
**`similarity_reduce_sympy` + `similarity_collapse_error` + `similarity_ode_solve`** (→ `core` when Ch. 9 Blasius/
Falkner–Skan calls them), `pipe_flow_regime`, `momentum_diffusivity`, `synthetic_millikan`. **Default g**: ch04, ch05,
ch07 and `core.waves` use `G_BOOK` = 9.81; the rest of `core` (incl. `lubrication`, `creeping`) G0 = 9.80665.
**Default ρ**: ch06 1.2 (2-D forces) vs 1000 (§6.9); ch07 1000; ch08 1000 (water), μ = 1e-3.

### Explainer engine (`assets/viz_lib.js`)
`Viz.app` (tabs Walkthrough / Explore / Explain / Derivation / Equations / Code / Check; fit-to-window with density
levels and pagers; transport, modes, presets, status, terms, inspector, notes, selftest parity), `Plot`, `Viz.field`,
`Viz.num` (RK4, odeint, brentq, erf/erfc **~1.2e-7 only — a defect; ch08 E5/E6 carry `erfcHP`**, niceTicks), `Viz.work`,
`Viz.three`, KaTeX with fallback; `Viz.font`, `Viz.roundRect`, `Viz.text`, `Viz.card`, `Viz.fmtTime`,
`Viz.tnum(keepTiny)`; pager slicing, `--viz-stage-need`, `ev.viewId`/`ev.vizView`. **Nothing promoted after ch02–ch08**
(the site-publisher was reading `viz/` each time). Ranked ch08 list in `knowledge/viz_patterns.md` (pager-label nowrap CSS
12 files, erfc precision, `arrowPx` 32, formatters incl. `pyNum`/`fmtLen`, layout probes → engine, gutter + log axes,
hatch, gamma, golden, derivation-result chain option, pager "now"-box slicing) on top of the ch07 list (`Viz.waves` and
per-mode hiding before Ch. 13). Reference: `templates/viz_example.html`.

## Explainer inventory
| chapter | slug | CORE idea | reusable stage pattern |
|---|---|---|---|
| ch01 | `continuum_averaging_volume` · `viscosity_momentum_diffusion` · `heat_work_paths` · `parcel_stability` · `buckingham_pi_machine` | C06 · C12 · C25/C35/C45 · C50–C55 · C64–C69 | molecules vs box size; gap with tracers + steady ghost; piston · p–v · T–s on one path; column + ζ(t) with a two-convention badge; chips → matrix → groups (exact rationals) |
| ch02 | `rotation_of_axes` · `cauchy_traction_principal_axes` · `strain_vs_rotation_split` · `gauss_flux_box` · `stokes_circulation_loop` | C02–C06 · C04/C05/C13 · C12 · C14/C15 · C11/C16 | clickable C matrix; element with a cut + Mohr; one clock, three squares; draggable box + waterfall flux bars; unrolled u·t along a loop |
| ch03 | `flow_lines_unsteady` · `material_derivative_probe` · `galilean_frames_cylinder` · `fluid_element_deformation` · `spin_and_principal_axes` · `vortex_paddle_wheels` · `reynolds_transport_cv` | C03/C04 · C01/C02 · C05 · C06–C09 · C10–C12 · C13/C14 · C15 | streamline/trail/dye on one clock; local/advective/total bars with measured ◇; observer slider; G sliders deform a square; paddle wheels; moving CV + budget waterfall |
| ch04 | `control_volume_budgets` · `stream_function_spacing` · `newtonian_stress_lab` · `navier_stokes_term_balance` · `rotating_frame_coriolis` · `which_bernoulli` · `viscous_dissipation_heating` · `boussinesq_buoyancy` · `dynamic_similarity_models` | C01/C04 · C03 · C07 · C08 · C09 · C05/C11/C12 · C10 · C13 · C15 | five scenes on one budget; ψ spacing = speed; G → S, R → τ → traction; 7 exact solutions with term bars summing to 0; two observers; hypothesis table; in = stored + out; kept vs dropped bars; paired group bars |
| ch05 | `vortex_tubes_cannot_end` · `vortex_pressure_funnel` · `kelvin_material_loop` · `baroclinic_torque` · `vorticity_stretching_tilting` · `biot_savart_filament` · `vorticity_equation_rotating` · `point_vortex_lab` · `vortex_sheet_rollup` | C01–C14 | 3-D tube + Gauss bars; needed = supplied force; ✓/✗ hypothesis table; torque ◇; Burgers balance; "build the sum"; conserved quantity flat in amber; click-to-place vortices; roll-up with KH ghost |
| ch06 | `superposition_sandbox` · `cylinder_circulation_lift` · `vortex_wall_images` · `complex_potential_corners` · `blasius_kutta_contour` · `conformal_joukowski` · `laplace_relaxation` · `axial_singularity_bodies` · `added_mass_sphere` | C03–C15 | element kits; term bars with ◇ ρUΓ; sensor trace; branch-cut status; Laurent bars; two planes, same particles; stencil stepping + honest order; moments converge, bars don't; speed + acceleration parts |
| ch07 | `dispersion_relation` | C03, C04 | tank + c(λ) with ghost limits and shaded regime bands + depth profile; sea-bed banner on phones |
| ch07 | `particle_orbits` · `capillary_gravity_waves` · `seiche_standing_waves` · `group_velocity_packets` · `wave_rays_refraction` · `hydraulic_jump` · `two_layer_modes` · `internal_wave_beams` | C05/C12 · C07 · C08 · C09 · C10 · C11 · C13/C14 · C15/C16 | linear vs exact paths + measured drift; restoring pushes + c² floor; right/left/sum ghosts; chord vs tangent + pond; in-page Hamiltonian rays; momentum/energy bars + forbidden preset; two magnifications; tank + square K plane on phones |
| ch08 | `couette_poiseuille_backflow` | C02 (D02–D05) | rose line + orange parabola = blue sum; floor tangent vertical at onset; amber reversed layer; Q waterfall; u_max dot at y* |
| ch08 | `lubrication_scaling` | C05 (D10, D11) | log term bars with a grey "dropped" band; pressure-scale chips; "thick gap ✗" preset |
| ch08 | `slider_bearing` | C06, C07 (D12–D14) | pad-frame profiles + tracers + pressure arrows; printed-formula ghost with a measured miss; wide-end recirculation; golden-section optimum |
| ch08 | `viscous_gravity_current` | C08, C10 (D15, D23) | conservative march with volume in the status; log–log front onto slope ⅕ beside a ½ ghost; glacier preset 10¹³ Pa s |
| ch08 | `stokes_first_problem` | C09 (D16–D20) | log clock + raw/rescaled mode; CN dots on erfc; presets that break similarity (stop at T, second wall) |
| ch08 | `similarity_exponents` | C10 (D21–D23) | exponent plane (n, m) with constraint lines; bracket-power chips; conserved-integral inset; four systems as modes |
| ch08 | `oscillating_plate` | C11 (D24, D25) | the rejected growing root drawn exploding; start-from-rest CN march; book depth with its 5.9 %; tide/eddy/Ekman presets |
| ch08 | `stokes_sphere_flow` | C13, C15 (D28, D29, D32, D33) | Stokes/Oseen/ideal on one control set; 1 % disturbance contour; measured inertia = friction circle |
| ch08 | `stokes_drag_settling` | C14 (D30, D31) | traction sweep (pressure + friction arrows = uniform push); running drag to 2πμaU + 4πμaU; particle table |

## Notation
See `knowledge/notation.md` (sign traps first, then the register). Conventions that matter everywhere: SI and kelvin
inside functions; **z up** (except ch06 §6.8, z horizontal); **lapse rate Kundu Γ ≡ dT/dz in code, meteorology shown**;
`p0` ≠ `p_ref`; pressures absolute unless `_gauge` (ideal flow: from hydrostatic; ch07: **gauge**); signed τ_xy = μ ∂u/∂y.
ch02: passive C; traction on the first index; A:B = A_ij B_ji; G[i, j] = ∂u_i/∂x_j; grid `[k, j, i]`. ch03: R = G − Gᵀ,
spin ½ω; u′ = u − U; ω′ = ω − 2Ω; RTT outward n. ch04: Coriolis acceleration +2Ω × u′ vs force −2Ω × u′; μ_v = 0;
**2-D ψ: u = ∂ψ/∂y (ch07: u = ∂ψ/∂z; GFD often opposite)**; four primes; two g defaults. ch05: ω vorticity; Γ circulation,
γ = u_below − u_above; ∇ρ × ∇p order; `ellipk(m)`, m = k². ch06: **Γ counterclockwise in code, `Gamma_cw=`/`Gamma_ccw=`**;
w = φ + iψ, dw/dz = u − iv; outside Zhukhovsky root; Stokes ψ [m³/s]. ch07: **ω = angular frequency**, ω ≥ 0 with
direction in sgn k; η = elevation; p′ per section; **g′ = g(ρ₂ − ρ₁)/ρ₂** (ch04 ρ₁); complex amplitudes with ½Re(AB*).
**ch08 (changes):** **dp/dx book sign (favourable < 0) with `G=` = −dp/dx alias (ch04)**; channel walls y = 0 fixed, y = h
moving; capital R cylindrical radius; **η = y/√(νt)** (ch04 and the book's figures: y/(2√(νt)); `half=True`); **ω =
vorticity (§8.1, §8.6) and oscillation frequency (§8.5)**; δ_e = √(2ν/ω) vs the book's 4√(ν/ω); **§8.6 θ from the
downstream axis** (rear stagnation θ = 0), `frame="body"|"fluid"`; four Re (pipe diameter; Re_L with ε²Re_L; generic;
**sphere 2aU/ν diameter vs Re_a radius**: Oseen 3/16 ↔ 3/8, far field ½Re_a r/a); lubrication p* = p/P_a (Λ) vs μUL/h²;
ch05 `rotating_cylinder_flow` takes ω = 2Ω₁ and `diffusing_vortex_sheet` γ = −2U for Example 8.5; `inlet_backflow` means
"recirculation at the wide end". Multi-meaning symbols (register): through ch07 as before (α, θ, T, h, R, σ, c, D, n, k,
q, γ, λ, τ, ε, H, Γ, ω, φ, A, C, G, K, m, s, t, a, b, primes, Φ, Ψ, Ω, M, f, F, B, C_p, η, ζ, κ, N, w, ξ, d, Q, ρ, g′, ν,
E); **ch08 adds η, ε (fineness), δ (three), θ, φ (Hele-Shaw potential), γ (ansatz), A–D (constants, D drag), E (field,
E²), n, h, U, f (friction factor), Λ**.

## Teaching lessons
- **Depth tiered, coverage exhaustive** (A / B / C, derivations, cells): ch01 15/70/21, 12, 496 · ch02 16/60/18, 15, 405
  · ch03 15/52/12, 24, 413 · ch04 15/151/19, 30, 611 · ch05 14/69/11, 23, 496 · ch06 15/116/29, 31, 564 · ch07
  16/193/27, 37 (321 steps), 617 · **ch08 15/113/15, 33 (301 steps), 501**.
- **Convention callouts** ("which p_o?", Coriolis term vs force, Γ_cw vs Γ_ccw, which ω/η/g′/p′, **which dp/dx sign, which
  η, which δ, which Re, which θ origin**): state each convention, give the size of the difference in numbers, say which
  the code uses — and carry **one numeric case through table → tiny example → code → explainer**.
- **Book slips are taught by computing both versions**: a **slips table up front** (printed | corrected | where), printed
  forms as ghosts that visibly fail and as code options a test must fail (ch07 frozen envelope; ch08 nine options);
  **check a suspected slip before teaching it** ((6.82) was not one; Stokes γ = 1 is right). The notebook's automatic
  equation insertion must not "correct" a sentence that names the printed form (ch08 lesson M2).
- **Make an approximation measurable** (ch01 linear vs full; ch07 residual scan; **ch08 measure every "~" claim**: "r ~ a/Re"
  → prefactor ½, crossover 2a/Re_a, direction dip; collapse spread 1e-16 vs 1.02 for wrong exponents) and label an extra
  approximation with its cost (ch08 slider linear load +15.6 % at α = 0.1 is not lubrication's assumption).
- **Every key number by two independent routes** (ch06 lift, added mass; ch07 (7.65) by FD eigenproblem; **ch08 three-way
  parity cells**: ch03 developing profile far downstream = ch04 exact solution = ch08 derived profile).
- **Every prose number is printed by a cell; captions count what they claim**; thresholds and crossovers computed, not
  quoted (ch07 Ursell 4π²/9; ch08 pipe radius for Re = 2000 by brentq, 2.00 mm, not a typed 2.5).
- **Flag inserted derivation steps** ("the book skips this move; we add it") — ch08 inserted four; the explainers must
  follow (the notebook is the step-count reference).
- **Primers come before first use** (tool reminders moved into the "Tools from earlier chapters" cell before the first
  use, ch08 S5); `knowledge/primers.md` lists 199 (P01–P199); ch09 starts at **P200**.
- **Builder rendering hygiene**: design TeX must be converted, not copied; **regex tidies never run inside maths they just
  created** (ch08 `$e^(…)$` ×20); headings keep their equations (`_title_no_numbers` ate "(8.29)").
- **Climate hooks with numbers** (ch04 Ro; ch05 tornado; ch06 ψ inversion; ch07 38 kW/m swell, 198 m/s tsunami, √(g′H)/f;
  **ch08 Stokes layer = Ekman layer with ω → f, spin-up L²/ν, cloud droplet 1.21 cm/s, sediment 0.36 mm/s, lava and ice
  as viscous gravity currents (μ to 10¹³ Pa s), cyclostrophic circular Couette**).
- Reusable derivation moves (ch01–ch08, 140; 20 from ch08) are tabulated in `knowledge/concept_map.md`.

## Global pitfalls confirmed in this project
**Machinery and environment**
- Extracted text garbles maths (`¼` = `=`, `ð…Þ` = parentheses, missing minus, ω → `u`, γ → g, θ → q, σ/τ → s, ν → n,
  ε → 3, ρ → r, Ω → U, ψ → j, φ → 4, η → h, ζ → z, ∂ → v, ∇ → V): read page images.
- **Windows: shell heredocs and `sed` mangle backslashes, quotes and `|`** — write code, explainers and knowledge files
  with Write/Edit or Python files; builder strings with LaTeX must be raw; single-backslash TeX in JS strings turns `\f`,
  `\t`, `\r` into control characters (not linted yet).
- `MPLBACKEND=Agg` while executing a notebook strips inline figures; **scripts calling `plt.show()` block headless runs —
  every chapter script takes `--no-show`** (ch07, ch08). Anaconda's kernelspec can shadow `fluidpy-venv`. JS `UIEvent.view`
  is read-only.
- `tools/shot.py` does not click/drag; its `py:` expressions have no builtins and no `np` (bind via lambda defaults); set
  parity tolerances near the achieved agreement; **run shot serially when memory is short and confirm a clip twice**
  (ch08); runtimes under parallel audits are ~5× slower (ch08 notebook ~196 s alone).
- **KaTeX-late pager clips are intermittent**: keep one line of slack in cards near the page height (ch07 E7) and use
  **text-size integrals in live rows** (ch08 E5: a display-size ∫ clipped one run in three).
- `test_viz_library_inlined_and_template_lints` fails whenever `assets/viz_lib.js` changed or an explainer is mid-build;
  run `tools/viz_inline.py --all` before the merge gate.
- matplotlib 3.11: animations using `subplots_adjust`/`fig.text` disable constrained layout for the save.
- `coverage_check` matches ledger reminders by exact primer name; it does not see literal TeX outside `$…$`.
- Book wording can hide in docstrings and **exercise answers can hide in computed numbers** (`check_public.py` catches
  neither): reviewers scan against the private JSON.
- Library: pager packs before the final `fit()`; **pager places a tall "now" box as a unit (blank phone pages, ch08)**;
  `Viz.num.erf/erfc` ~1.2e-7; fixed transport range; no per-mode controls (Q4); `onChange` after `set` (Q5). **Phone
  legibility is the usual viz failure** (rotated y titles, labels over labels, colour-only legends, squashed strips).
- **Design and analysis documents contain errors; downstream agents compute, not copy** (ch05 three; ch06 four; ch07
  three; **ch08: an "8 Pa/m" preset already past zero net flow, ice at 10⁴ instead of 10¹³ Pa s, two oil densities, a
  "Re = 2" bead that is not creeping, a D01 check "15.1" (15.007), a D32 check low by 2**). **Explainers built from design
  Part F drift when the notebook adds moves** — build Derivation tabs from the notebook dump.

**Mathematics → code**
- **A test comparing two of our own functions is not evidence for the book's convention.** Pin conventions to fields
  with known answers and **prove discrimination by planting the wrong variant** (ch02 8/8 … ch07 41/41, **ch08 41/41**,
  incl. the old code of every review fix).
- **Tests must exercise every parameter at a non-zero value, with both signs** (ch08 M1: the Couette–Poiseuille extremum
  sign hid behind U = 0 tests); **derive a regime flag from the field and brute-force it over the whole range the UI
  allows** (ch08: slider recirculation checked at x = L only; true rule gap ratio > 2, α > 1 or α < −½ — found by an
  explainer builder's independent JS); **every asymmetric solution gets a direction test** (ch08 Oseen wake).
- **A verifier's correction is a claim too**: compare magnitudes of vector terms, not one component (ch08 loop-1 "⅛,
  crossover 10–20 a/Re_a" overturned by the implementer: ½Re_a r/a, crossover 2a/Re_a; pinned by sympy + asymptote).
- **Complex-step derivatives need analytic continuations of every non-analytic piece** (ch07 M1); keep a speed (|c|) and a
  signed derivative (dω/dk) as different functions.
- **Printed formulas can assume a sign** ((7.138), (7.145) need k > 0): code the gradient form, test all quadrants.
- **Code the corrected form; keep the printed one as a named option a test must fail** (`form="book"`, `model="book"`,
  `printed_8_13b=True`, `printed=True`); **keep the book's sign in the argument (`dpdx`) with the earlier chapter's as a
  keyword alias (`G=`)**.
- **Before calling an equation a book slip, compare the whole printed form with the reference form times every plausible
  factor** ((6.82)); a literal exercise set-up may give a different coefficient (Stokes γ).
- **Overflow and cancellation**: cosh/sinh ratios in the e^{kz}(1 ± e^{−2k(z+H)})/(1 ∓ e^{−2kH}) form; x/sinh x via expm1;
  erfc, never 1 − erf, at large η; (1 − e^{−s})/s via −expm1; small-α series for the slider load; exact branches for
  R₂ = ∞ and R₁ = 0 (never 1e300); NaN inside bodies; t ≤ 0 guards.
- **Stiff and degenerate PDEs**: Crank–Nicolson after a sudden jump oscillates — two backward-Euler start-up steps restore
  order 2 (pure CN stalls); degenerate nonlinear diffusion (D ∝ h³) needs implicit Newton/Picard, a precursor film and
  conservative faces; KdV needs an integrating factor and 2/3 dealiasing.
- **Inverse relations by a physics bracket** (brentq, residual asserted) and iterative searches seeded at the problem's
  own length scale; a double root is located only to √tol.
- **Complex branches are physics** (Zhukhovsky outside root, `arctan2`); **complex amplitudes need real parts before
  products** (½Re(AB*)); keep only the bounded root (ch08 Stokes layer).
- **Singular corners and staircase boundaries set the global convergence order** (ch06 4/3; ch08 Hele-Shaw 1.11);
  asymptotic orders are stated relative or absolute; condition numbers reported with every inverse method.
- **Finite-difference tolerances from a written budget** (stencil truncation + ε·r/h_rel² round-off, tripled by
  Richardson) — never loosened by feel (ch08 ratio test).
- **A telescoping sum is not a conservation test**; genuine V4 invariants so far: ch07 ∫η², KdV invariants, ω along rays;
  **ch08 thin-film volume, vortex-sheet jump, ∫ω dy = +U, pipe CV force balance, circular-Couette torque and power =
  dissipation**.
- **Sibling functions share defaults and argument meanings** — or are documented and called by keyword; convention-named
  keywords (`Gamma_cw=`, `ref="lower"`, `G=`, `frame=`, `half=`, `p_scale=`).
- **Two quantities with one name**: g′ (ρ₁ vs ρ₂), η (level set, elevation, y/√(νt) vs y/(2√(νt))), p′ (three bases), δ
  (e-folding vs 4√(ν/ω)), Re (radius vs diameter) — the code names the reference in the keyword or the function name.
- **Benchmarks from the primary source at its printed precision** (Fenton & McKee 1.7 %; San Andrés K_opt 2.1889).
- Signed and ordered quantities stay signed: b·n, wall stress (τ₀ < 0 for forward flow), lapse rate, ∇ρ × ∇p, u₂ − u₁,
  Γ_cw vs Γ_ccw, dw/dz = u − iv, dω/dk, dp/dx, power into the fluid −2πR₁σ_Rφu_φ. Probe piecewise functions at their
  joins; check a field against its own stream function; test every scenario's full balance.
- Library conventions: `ellipk(m)` with m = k²; `np.gradient` first order at edges; principal branches; pint °C offset;
  kmol vs mol; geopotential vs geometric altitude; FFT evolution is periodic.
- **Validation labels**: V4 only for real conservation laws; identities V1; published closed forms reproduced identically
  are V1 form cross-checks, not V5; "converged" only with an asserted order (ch08: Gauss–Legendre exact from n = 1 is
  analytic; `thin_film_spread` is conserved + benchmark, not converged); a correlation bracket is V5 for the bracket only.

## Benchmark inventory
| value / table | source (citation) | stored in | used by |
|---|---|---|---|
| k_B, N_A, R (exact); e = 1.602176634e-19 C | CODATA 2018/2022, NIST CUU | `reference/ch01/constants.json`, `reference/ch08/` | ch01 `thermo`; ch08 synthetic Millikan (5.1e-4) |
| USSA-1976 constants, Tables 1–2; sea-level c = 340.29 m/s; Sutherland β, S | NASA-TM-X-74335; PDAS tables | `reference/ch01/ussa1976_*` | ch01, ch04 C02, ch08 ν_air |
| IAPWS σ(T), μ(T) | IAPWS R1-76(2014), R12-08 | `reference/ch01/` | ch01, ch04, ch07 C07 (72.74 mN/m), ch08 ν_water |
| Jennings mean free path; Taylor blast K = 0.856; AMS Γ_d; EOS-80 | Tsalikis et al. 2024; Díaz arXiv:2009.05674; AMS Glossary; Fofonoff & Millard 1983 | `reference/ch01/benchmarks.json` | ch01 |
| Divergence theorem 8π/3; Levi-Civita; Mohr; rotations; Stokes | Wikipedia; scipy docs | `reference/ch02/` | ch02 |
| Lamb–Oseen vortex peak α = 1.256 | Canivete Cuissa & Steiner 2022, arXiv:2210.05223 | `reference/ch03/` | ch03 C14 |
| Sphere drag C_D(Re), 0.1–1e6 | F. A. Morrison (2016), Michigan Tech | `reference/ch04/` | ch04 C15; **ch08 C15 bracket (between Stokes and Oseen for 0.1 ≤ Re ≤ 5)** |
| WGS-84 a, 1/f, b, ω | Wikipedia "World Geodetic System" | `reference/ch04/` | ch04 C09, ch05 |
| Orifice C_c = 0.611; capillary length 2.71 mm; Prandtl ranges | Wikipedia | `reference/ch04/` | ch04 |
| Sphere added mass = ½ displaced mass; Rayleigh collapse 0.91468 | Lamb (1932) §92; Rayleigh, Phil. Mag. 34 (1917) 94 | `reference/ch06/` | ch06 |
| Source panels (pedigree) | Hess & Smith (1967), doi:10.1016/0376-0421(67)90003-6 | `reference/ch06/SOURCES.md` | ch06 C14 |
| Explicit inverse dispersion: Fenton & McKee 1.7 % (ours 1.658 %); Guo β = 2.4908 → 0.75 % (0.753 %) | Coastal Eng. 14, 499–513 (1990); 45, 71–74 (2002) | `reference/ch07/` | ch07 C03 |
| Capillary–gravity minimum 0.23 m/s at 1.7 cm (ours 0.2312 at 1.712 cm); Stokes drift, limiting H/λ = 0.1410633 | Wikipedia "Capillary wave", "Stokes drift"; HandWiki "Stokes wave"; arXiv:1507.02784 | `reference/ch07/` | ch07 C07, C12 |
| **Inclined slider bearing: K_opt = 2.1889 (ours 2.18870), W(K_opt) = 0.0267 (0.026707), peak-pressure K = 2.414, P_max = 0.043** | L. San Andrés, *Modern Lubrication*, Notes 2 Appendix (Texas A&M 2009, rev. 2012), rotorlab.tamu.edu | `reference/ch08/` | ch08 C07 (V5) |
| **Planar viscous gravity current η_N = 1.411 (ours 1.411245), similarity shape, t^{1/5}** | Ball & Huppert preprint App. (a); Huppert, J. Fluid Mech. 121, 43–58 (1982) | `reference/ch08/` | ch08 C08 (V5) |
| Forms (V1 cross-checks, not V5): Rankine, Lamb–Oseen, RTT; Bélanger, Taylor–Green; Kelvin ring, Burgers, Hill; cylinder, Kutta–Joukowski, Blasius; KdV/cnoidal/solitary; **Oseen/Proudman–Pearson, Stokes' law, Stokes' second problem, Hagen–Poiseuille, Taylor–Couette A, B (ch08)** | Wikipedia pages; K. T. McDonald | `reference/ch03/`–`reference/ch08/` (`make_refs.py`, `SOURCES.md`) | ch03–ch08 |
Book-printed values (private, git-ignored): `tests/book_values_ch01.json` … `…_ch08.json` (ch08: every §8.2–§8.6 closed
form ≤ 1e-12; printed slips confirmed as slips). **Open**: a measured high-Re cylinder C_p (ch06 N31); a primary number
for the minimum group velocity (ch07 O7); a primary pipe-transition number (ch08 keeps the book's 2000–3000, V6).

## Validation summary per chapter
Counts are concept-map rows carrying each label (a row can carry several); tests per evidence level in parentheses.
| chapter | analytic | symbolic | converged | conserved | benchmark | book-value | qualitative | unverified | verdict |
|---|---|---|---|---|---|---|---|---|---|
| ch01 | 39 (V1 58) | 16 (V2 19) | 5 (V3 8) | 5 (V4 7) | 18 (V5 13) | 1 (V6 5) | 1 (`wavelength_to_rgb`) | 0 | **PASS** — 107 tests; 5 explainers (115 rows); 496 cells |
| ch02 | 43 (V1 41) | 21 (V2 13) | 11 (V3 15) | 2 (V4 1) | 9 (V5 6) | 8 (V6 4) | 0 | 0 | **PASS** — 92 tests; 8/8 wrong variants; 5 explainers (127 rows); 405 cells |
| ch03 | 38 (V1 58) | 23 (V2 32) | 14 (V3 14) | 5 (V4 4) | 2 (V5 1) | 8 (V6 4) | 0 | 0 | **PASS** — 120 tests; 13/13 (+7); 7 explainers (145 rows); 413 cells |
| ch04 | 41 (V1 83) | 35 (V2 41) | 17 (V3 9) | 13 (V4 6) | 11 (V5 6) | 5 (V6 3) | 0 | 0 | **PASS** — 151 tests; 26/26; 9 explainers (255 rows); 611 cells |
| ch05 | 39 (V1 90) | 24 (V2 23) | 13 (V3 10) | 13 (V4 8) | 0 (V5 0) | 5 (V6 4) | 2 (fenced) | 0 | **PASS** — 135 tests; 40/40 + 5; 9 explainers (268 rows); 496 cells |
| ch06 | 33 (V1 86) | 24 (V2 27) | 14 (V3 10) | 16 (V4 6) | 3 (V5 2) | 7 (V6 4) | 2 (labelled) | 0 | **PASS** — 137 tests; 30/30; 9 explainers (291 rows); 564 cells |
| ch07 | 38 (V1 60) | 28 (V2 36) | 17 (V3 10) | 8 (V4 5) | 5 (V5 5) | 3 (V6 3) | 0 (2 labelled extensions) | 0 | **PASS** — 132 tests (+ V7 13); 41/41; 9 explainers (222 rows); 617 cells |
| ch08 | 29 (V1 46) | 20 (V2 31) | 6 (V3 5) | 7 (V4 4) | 7 (V5 6) | 0 rows (V6 3 tests) | 0 | 0 | **PASS** loop 2 — 101 tests (+ V7 6); 41/41 wrong variants; review 1 Must (extremum sign) + 7 Should resolved, slider rule generalised; all 25 ★★/★★★ D re-derived in sympy; 9 explainers PASS by round 3 (278 rows, 207 vs fluidpy, 3 312 views, 26 D matched); notebook PASS round 2 (501 cells, ~196 s) |
(V7 tests: ch01 19, ch02 12, ch03 7, ch04 3, ch06 2, ch07 13, ch08 6. Full suite at the ch08 merge gate: **983 tests**.)

## Open across chapters
- **Library pass (orchestrator, when nothing reads `viz/`)**: the ranked ch08 list in `knowledge/viz_patterns.md`
  (pager-label nowrap CSS, **erfc precision**, `Viz.arrowPx`, formatters incl. `pyNum`/`fmtLen`, layout probes → engine,
  gutter + log axes, pager "now"-box slicing, derivation-result chain option, hatch, gamma, golden) plus the ch07 list
  (**`Viz.waves` and per-mode hiding before Ch. 13**, `Viz.cx`, Gauss–Legendre); then `tools/viz_inline.py --all` and
  `tools/shot.py --chapter ch01…ch08 --quick` + `templates/viz_example.html --quick` (serially if memory is short).
- **Skill pass**: lesson candidates for `interactive-viz`, `math-to-python` §7, `verify-implementation`, `teaching-style`
  and `colab-notebook` in `knowledge/viz_patterns.md` (ch02–ch08) are not yet in the skills (ch08 first: notebook = step
  reference, non-zero both-sign parameter tests, brute-forced regime flags, direction tests, magnitude comparisons,
  measure every "~").
- **Machinery TODO**: nbkit — `show_eqs` skips "printed" mentions, regex tidies only outside `$…$` with a build guard,
  headings keep equations, `tidy_raw_tex` into nbkit; a scripted Derivation-tab vs notebook-dump check in `embed_check`;
  `viz_lint` single-backslash TeX rule; `run_notebook` refuses `MPLBACKEND=Agg`; `coverage_check` literal TeX; `shot.py`
  clicks/drags and exposes `np` to `py:` rows; pager re-pack after KaTeX; rule-9 scan against the private JSON; `eq_refs`
  exempts `<meta>` and `ref:` badges; re-point ch04 `linear_wave_surface` to `core.waves.cosh_over_sinh`; an alias
  `wide_end_backflow` for `slider_bearing_state["inlet_backflow"]`.
- Non-blocking explainer follow-ups: ch01–ch07 (their `knowledge/chNN.md` §9); **ch08** E1 (narrow "Q" label), E2 (phone-
  land names, D10 step 8), E3 ("strongest suction", phone D14 s14 blank), E5 (phone rescaled label), E6 ("ansatz"), E7
  (δₑ glyph), E8 (phone fluid-frame view, D29 s5 blank), E9 (phone Stresses legend) — `knowledge/ch08.md` §9. Backups not
  built: ch04 `kinematic_free_surface`, ch05 `vortex_rings`, ch06 `flow_net_sources_vortices`, ch07
  `linearised_free_surface`, ch08 `rotating_cylinders_couette`.
- Verification notes: ch08 O5 Hele-Shaw staircase order 1.11, O6 test runtime ~75 s, O8 four docstrings to name the
  loop-2 tests, O9 `inlet_backflow` naming, settling "valid" Re < 0.1 (ours); ch07 O5 Exercise 7.6 +0.72 % (book value),
  O6 Guo search-confirmed, O7 no primary c_g,min; ch06 O4, O6; ch04 O5 Prandtl numbers; ch05 O2 tube-flux route 1.8e-4.
