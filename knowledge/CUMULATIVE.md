# CUMULATIVE knowledge — fluidpy
(Rewritten by the knowledge-keeper after every chapter. Every agent reads this first. Last rewrite: after ch06,
2026-09-24, commit b1bd4c5.)

## Project rules in force (added during ch01–ch03; unchanged in ch04–ch06)
- **5–10 explainers per chapter, as many as the CORE ideas need** (`book.yaml → project.min/max_explainers_per_chapter`).
  ch03 has 7, ch04–ch06 9 each.
- **Every book equation is shown in full next to its number** — notebook prose, derivation steps and headings, traps,
  recaps and explainer tour/Explain/Derivation/quiz/notes/status text **and `<meta>` strings** (ch04 E6, ch06 E6).
  Enforced by `tools/coverage_check.py` check 8 (markdown) and `tools/eq_refs.py` (explainers). The lesson reviewer also
  scans by eye: a cell with *other* maths passes check 8 while citing a bare number (ch05, ch06 round 1: ~20 hits).
- **`coverage_check` check 9**: no control characters in cells (a LaTeX backslash eaten by a non-raw string). The same
  bug in explainer JS strings is not linted yet (ch06 E2) — see Open.
- **Pager and audit**: an item taller than a page is sliced at natural breaks; `tools/shot.py` pages through every page
  of every pager at every size.

## Physics pipeline so far
```
ch01 Introduction (done)
  molecules ──(box average, D35 noise N^-1/2, Kn = l/L)──► continuum fields ρ, p, T, u at a point      [C06]
     ├─ molecular transport: flux = −D·gradient; Newton τ = μ du/dy; ν = μ/ρ; FTCS diffusion      [C12]
     ├─ statics: dp/dz = −ρg (D05), buoyancy ρgV (D37), USSA-1976, H = RT/g                        [C20]
     ├─ thermodynamics: first law [C25] → Gibbs (D10) [C35] → c² [C36] → p = ρRT (D11) [C40] → p/ρ^γ (D14) [C45]
     ├─ stratification: ζ″ + N²ζ = 0 (D18) [C50], N² [C51], Γ_a = −gαT/C_p (D19) [C54], θ (D20, D36) [C55]
     │     (Kundu Γ ≡ dT/dz in code; meteorology −dT/dz always shown)
     └─ dimensional analysis: homogeneity [C64] → matrix, rank [C67] → Π groups = null space (D28) [C69]

ch02 Cartesian tensors (done) — the mathematical language
  summation [C01] → C_ij, CCᵀ = I (D02, D03) [C02] → x' = Cᵀx (D01) [C03]; stress τ_ij [C04] → Cauchy f_i = τ_ji n_j
  (D05) [C05] → τ' = CᵀτC (D06) [C06] → invariants (D18) [C07]; ε_ijk, ε–δ (D09) [C08]; S + A, R = G − Gᵀ ↔ ω (D14,
  D15) [C12]; principal axes (D17 ★★★) [C13]; ∇φ, ∇·u, ∇×u on the grid [C09–C11]; Gauss (D25) [C14] → integral
  definitions (D21, D22) [C15]; Stokes (D26) [C16]

ch03 Kinematics (done) — how fluid moves, before any force
  particle ⇄ field (3.2) (D01) [C01] → DF/Dt = ∂F/∂t + u·∇F (3.5) (D02) [C02]; streamlines, path and streak lines
  (D03–D05) [C03, C04]; Galilean invariance (3.9) (D06) [C05]; du = G·dx (D07) [C06] → n·S·n, S₁₂, ∇·u, spin ½ω,
  du = S·dx + ½ω × dx (D08–D15) [C07–C11]; principal axes (D16) [C12]; polar ω_z, Rankine, Gaussian (D17–D20)
  [C13, C14]; Leibniz (D21) → RTT (3.35) (D22 ★★★) [C15] → (3.14) (D23), (3.5) + F∇·u (D24)

ch04 Conservation laws (done) — the trunk
  material volume + coincident CV + RTT (3.35): mass (4.5) (D01) [C01] → continuity (4.7)–(4.10) (D02, D03) [C02] → ψ
  (D04) [C03]; momentum (4.17) (D05) [C04], Bernoulli (4.19) (D06) [C05] → Cauchy (4.24) (D07) [C06]
  constitutive τ = −pδ + 2μS + λS_mmδ (4.31) (D08, D09 ★★★, D10) [C07] → NS (4.38)–(4.39b), −μ∇×ω (D11–D13) [C08]
  rotating frame +2Ω×u′ (4.43) / −2Ω×u′ (4.45) (D14–D18, D15 ★★★) [C09]; energy, ε ≥ 0 (D19–D23) [C10]
  Lamb identity → ∇B = u × ω (4.69) (D24, D25) [C11]; unsteady Bernoulli (4.75) (D26) [C12]; Boussinesq (4.86), (4.89)
  (D27, D28) [C13]; Dη/Dt = 0 (4.91), jumps, Laplace (D29) [C14]; dimensionless NS (4.101) (D30) [C15]

ch05 Vorticity dynamics (done) — vorticity's life
  tubes cannot end (5.4) (D01) [C01]; pressure of the basic vortices (5.6), (5.7) (D02, D03) [C02]; Kelvin DΓ/Dt = 0
  (5.8) (D04, D05; needs ρ = ρ(p)) [C03]; baroclinic torque ∇ρ×∇p/ρ² (D06, D07) [C04]; Helmholtz (D08) [C05];
  Dω/Dt = (ω·∇)u + ν∇²ω (5.13) (D09 ★★★) [C06]; Biot–Savart (5.16) with +1/(4π) in (5.14) (D10, D11 ★★★) [C07];
  filament law (5.17) (D12, D13) [C08]; rotating baroclinic (5.30) (D14, D15 ★★★) [C09]; stretching/tilting (5.32)
  (D16–D18) [C10]; Γ_a and (ζ + f)/h (5.33) (D19, D20) [C11] ← seed of PV (Ch. 13); point vortices (D21) [C12];
  images (D22) [C13]; vortex sheet γ = u₂ − u₁ (D23) [C14]

ch06 Ideal flow (done) — "irrotational + incompressible ⇒ Laplace; solutions add; any streamline can be a wall"
  ω = 0 kills μ∇²u = −μ∇×ω ⇒ (6.1), only no-through-flow kept (D01) [C01]
  ω_z = −∇²ψ (6.4) (D02); ln r is harmonic but has flux 2π ⇒ vortex/source = δ sources (D03); φ ⟂ ψ (D04) [C02]
  superposition + ∂ψ/∂s = 0 on a wall (6.16) (D05) + Bernoulli everywhere (6.18) [C03]
  element kit → doublet = pair limit, d from sink to source (6.29) (D06) [C04] → half-body (6.31) (D07, D08) [C05]
     → cylinder (6.33), C_p = 1 − 4 sin²θ, D = 0 (D09) [C06] → + clockwise Γ: stagnation (6.38), L = ρUΓ (6.40)
     (D10, D11) [C07] → images, Example 6.1 wall pressure by unsteady Bernoulli (D12, D13) [C08]
  w = φ + iψ, Cauchy–Riemann, dw/dz = u − iv (D14); corner w = Azⁿ, speed ∝ r^{n−1} (D15) [C09]
     → Blasius D − iL = (iρ/2)∮(dw/dz)²dz (6.60) (D16, D17 ★★★) → Kutta–Zhukhovsky L = ρUΓ for any body (D18 ★★★) [C10]
     → conformal maps keep angles (D19); Zhukhovsky z = ζ + b²/ζ: circle → ellipse, outside root (D20, D21 ★★★) [C11]
  on a grid: average rule (6.72) (D22), Aψ = b, Jacobi/GS/SOR (D23); Example 6.2 order 4/3 (270° corner) [C12]
  axisymmetric: Stokes ψ (6.75), (6.77) ≠ Laplace (D24); 3-D elements, sphere C_p = 1 − (9/4)sin²θ (D25) [C13]
     → airship (6.95) (D26), axial singularity method (D27), source panels (ours) [C14]
  accelerating sphere: (6.97) (D28) → pressure (6.105) (D29 ★★★) → M = ½ displaced mass (6.108) (D30), by energy
     (D31 ★★★) [C15]

next: ch07 gravity waves — φ with ∇²φ = 0 (ch06 C02, `core.potential`), unsteady Bernoulli (ch04 C12, ch06 D13) as the
      dynamic free-surface condition, Dη/Dt = 0 (ch04 C14) as the kinematic one, images/bottom BC (ch06 D05, D12),
      N² and hydrostatics (ch01) for internal waves, added mass (ch06 C15) for wave forces; `Viz.cx` would help.
```
The book ahead: Ch7 gravity waves · Ch8 laminar flow → Ch9 boundary layers → Ch10 CFD → Ch11 instability → Ch12
turbulence → **Ch13 GFD (uses 4, 5, 7, 8, 11, 12)** → Ch14 aerodynamics · Ch15 compressible · Ch16.
Where earlier chapters feed in: ch01 N², θ → Ch. 7 §7.8, 11, 13; τ, ν → Ch. 8, 9; hydrostatic base → Ch. 7, 13;
isentropic gas → Ch. 15; FTCS → Ch. 8, 10. ch02 Gauss → every CV; Stokes/circulation → Ch. 13 PV; invariants → Ch. 12.
ch03 vortices → Ch. 8, 13, 14; path lines → Ch. 7, 10; principal strain → Ch. 12, 13. ch04 **NS residual/term tools →
every exact solution**; unsteady Bernoulli + kinematic BC → Ch. 7; CV budgets → Ch. 9, 14, 15; ε → Ch. 12; rotating
frame, Boussinesq, Ri, Ro → **Ch. 13**. ch05 vorticity budget → Ch. 8, 10, 12; sheet roll-up → Ch. 11; **(5.30), (5.33),
(ζ + f)/h → Ch. 13**; filaments, rings, images → Ch. 14. **ch06**: `core.potential` outer flows → **Ch. 9** (U_e(x),
dp/dx; wedge Azⁿ → Falkner–Skan, n = 2 → Hiemenz) and **Ch. 14** (lift, Kutta, circle theorem); `core.conformal` →
Ch. 14 airfoils; `core.laplace_solvers` → **Ch. 10** (Poisson baseline, ω–ψ) and Ch. 13 (ψ/PV inversion in a basin);
`core.panels` → Ch. 14 panel methods; moving sphere/added mass → Ch. 7 wave forces, Ch. 13 parcels, Ch. 16.

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

### Physics primitives from ch01–ch03 (re-exported as `chNN.<name>`)
| function family | module | book § / Eq. | label |
|---|---|---|---|
| constants `K_B, N_A, R_U, R_AIR, GAMMA_AIR, CP_AIR, G0, P_ATM, P_REF`, `G_BOOK` = 9.81; perfect gas, first law, entropy, Gibbs, sound speed, isentropic ratios, van der Waals | `thermo` | §1.8–1.9, (4.94) | analytic, symbolic, converged, conserved, benchmark |
| `integrate_hydrostatic`, `standard_atmosphere`, `scale_height`, `buoyancy_force`, `net_pressure_force_on_box` | `statics` | (1.8), (1.9) | analytic, converged, benchmark |
| N² family, `classify_stability`, `parcel_displacement`, `adiabatic_lapse_rate` (negative), `lapse_rate_convention`, `lapse_rate_stability`, θ, ρ_θ | `stratification` | (1.29)–(1.35) | analytic, symbolic, benchmark |
| `ftcs_diffusion_1d`, `stable_time_step`, `couette_startup_profile`, `gaussian_spreading` | `diffusion` | §1.5 | converged, conserved, analytic |
| `dimensional_matrix`, `rank_by_minors`, `pi_groups`, `rescale_units` | `dimensional` | (1.36)–(1.40) | symbolic, analytic |
| index notation; `tensors` (passive C, `traction` first index, `double_dot(convention=)`, `invariants`, `rotation_tensor`, `principal_axes`); grids `[k, j, i]` + `operators` (`vector_gradient` G[i, j] = ∂u_i/∂x_j, `curl`, `laplacian`); `fields`; `integral_theorems` (Gauss, Stokes, `circulation`, `gauss_legendre_nodes`, `curl_theorem_box`) | ch02 modules | (2.1)–(2.36) | analytic, symbolic, converged, conserved, benchmark |
| `kinematics` (field callables `u(x, t)`, material derivative, flow lines, strain/spin, presets); `coords`; `transport` (Leibniz, CV shapes, `reynolds_transport`) | ch03 modules | (3.1)–(3.35) | analytic, symbolic, converged, conserved |
| `vortices`: solid body (rate), line, Rankine, Gaussian; ch05 Burgers, Hill, Lamb–Oseen, Gaussian tubes, rings, ABC | `vortices` | (3.22)–(3.29), §5.4 | analytic, symbolic, converged, benchmark |

### Conservation-law primitives from ch04 (`ch04.<name>`; details `knowledge/ch04.md` §3)
| function family | module | book § / Eq. | label |
|---|---|---|---|
| `mass_budget`, `momentum_budget`, `energy_budget`, `angular_momentum_budget` on any `ControlVolume`; flux (u − b)·n | `conservation` | (4.5), (4.17), (4.48), (4.65) | analytic, conserved, converged, symbolic |
| `newtonian_stress`, `viscous_stress`, `mean_pressure(tau, tau33)`, `bulk_viscosity`, `stress_on_plane`, `dissipation_rate` | `constitutive` | (4.25)–(4.37), (4.58) | analytic, symbolic |
| residuals and term splits (numeric + sympy): continuity, Cauchy, NS, `viscous_force_forms`, `exact_solution` (8), `lamb_vector`, energy/entropy, Boussinesq | `navier_stokes` | (4.7)–(4.89) | analytic, symbolic, converged, conserved |
| Appendix-B operators (cylindrical, spherical), sympy | `curvilinear` | App. B | analytic, symbolic |
| ψ ↔ u (planar u = ∂ψ/∂y, axisymmetric), fluxes, presets | `streamfunction` | (4.11), (4.12) | analytic, symbolic, converged |
| `frame_acceleration_terms` (+2Ω×u′), `apparent_body_forces` (−2Ω×u′), `coriolis_parameter`, `effective_gravity`, `OMEGA_EARTH` | `rotating` | (4.42)–(4.45) | analytic, symbolic, converged, benchmark |
| `bernoulli_function`, `pressure_function`, `which_bernoulli(_text)`, unsteady family, stagnation, pitot | `bernoulli` | (4.19), (4.66)–(4.83) | analytic, symbolic, conserved, benchmark |
| `kinematic_bc_residual`, `interface_mass_flux`, `laplace_jump_from_balance`, `capillary_length` | `interfaces` | (4.90)–(4.98) | analytic, converged, benchmark |
| `Scales`, 20+ named numbers (Re … Ro), C_p/C_D/C_L, `sphere_drag_coefficient` (Morrison) | `similarity` | (4.99)–(4.119) | analytic, symbolic, benchmark |

### Vorticity primitives from ch05 (`ch05.<name>`; details `knowledge/ch05.md` §3)
| function family | module | book § / Eq. | label |
|---|---|---|---|
| loops (`circle_loop_points`, `loop_circulation` spectral, `loop_vector_area`), `material_loop`, `material_circulation`, `kelvin_rate_terms`, `kelvin_force_terms`, `absolute_circulation` → (Γ, Γ_a) | `vorticity` | (5.8)–(5.11), (5.33) | conserved, analytic, converged, symbolic |
| `vorticity_field`, `vortex_line`, `vortex_tube_strength`, `vortex_tube`, `tube_flux_budget`, `frozen_in_check` | `vorticity` | (5.3)–(5.4), §5.3 | analytic, symbolic, conserved, converged |
| **`vorticity_budget(u, x, t, nu, Omega, rho, p)`** → local, advective, relative, planetary, baroclinic, diffusion, residual; `_sym`; presets; `baroclinic_term`, `planetary_vorticity_terms`, `stretching_tilting_split` | `vorticity` | (5.12)–(5.32) | symbolic, analytic, converged |
| `poisson_green_3d`, `velocity_from_curl_omega(sign=+1)`, `biot_savart_volume`, `biot_savart_2d`, `velocity_from_vorticity_fft` (periodic) | `biot_savart` | (5.14)–(5.16) | analytic, converged, symbolic |
| segments, filaments, rings (`ring_ring_velocity` m = k², `ring_self_velocity`), point vortices (`point_vortex_evolve(boundary=wall/circle/channel)`, invariants, images), sheets | `biot_savart` | (5.17), §5.7–5.8 | analytic, conserved, converged |

### Ideal-flow primitives from ch06 (`ch06.<name>`; details `knowledge/ch06.md` §3)
| function family | module | book § / Eq. | label |
|---|---|---|---|
| elements `Uniform(U, V)`, `Source(m, z0, cut_angle)`, `Vortex(Gamma, z0, cut_angle)` (ccw), `Doublet(d_vec)` / `.from_book_scalar`, `Corner(A, n, cut_angle, rotate)`, `Constant`; `Flow([...])`, `FunctionFlow(w, dwdz)`: `.w`, `.dwdz`, `.velocity`, `.psi`, `.phi`, `.cp`, `.pressure(rho=1.2)`, `.stagnation_points()`, `.contributions`; `log_branch`, `power_branch` | `potential` | (6.4)–(6.53) | analytic, symbolic, conserved, converged |
| bodies `half_body`, `cylinder(Gamma_cw=/Gamma_ccw=)`, `gamma_ccw_from`; images `mirror(wall=)`, `circle_theorem`; forces `blasius_force(flow, R \| contour, rho)`, `laurent_coefficients` (FFT), `contour_crosses_body`; checks `laplacian_residual`, `normal_velocity_on`, `far_field_check`, `polar_velocity(_sym)`, `pressure_coefficient` | `potential` | (6.16)–(6.62) | analytic, conserved, converged, symbolic |
| axisymmetric `AxisymUniform`, `PointSource3D`, `Doublet3D`, `LineSource3D`, `AxisymFlow` (Stokes ψ [m³/s]), `sphere`, `sphere_potential_vector(a=)`; moving sphere `moving_sphere_potential/velocity/surface_velocity/dphidt/surface_pressure`, `sphere_force_quadrature`, `added_mass_sphere` | `potential` | (6.74)–(6.109) | analytic, symbolic, conserved, converged, benchmark |
| `joukowski(_inverse)(_derivative)` (outside root), `circle_flow_zeta(U, a, Gamma_cw, alpha)`, `mapped_flow`, `grid_image`, `map_grid_lines`, `map_elements`, `angle_preservation`, `conformal_map` | `conformal` | (6.63)–(6.69) | analytic, symbolic |
| masked-grid `laplacian_5pt`, `jacobi_sweep`, `gauss_seidel_sweep`, `sor_sweep`, `residual_norm/field`, `assemble_laplace` (sparse), `solve_laplace(mask, bc, method, tol)`, `solve_poisson`, `optimal_sor_omega`, `jacobi_spectral_radius` | `laplace_solvers` | (6.70)–(6.73) | converged, analytic, conserved |
| `axial_influence_matrix`, `axial_singularity_solve` (odd-N symmetric ⇒ min-norm LSQ + warning), `axial_singularity_psi/velocity`; 2-D `source_panels`, `panel_geometry`, `panel_induced_velocity`, `panel_velocity` | `panels` | (6.93)–(6.95), §6.8 | analytic, converged, conserved |
Chapter-only code stays in chapter modules (move to `core/` when a second chapter calls it): ch03 `pipe_profile` (→ Ch.
8); ch04 `linear_wave_surface`, `bore_speed` (→ Ch. 7), `couette_heating`, `two_fluid_couette`, `stokes_first_problem`
(→ Ch. 8), `wake_drag_per_span` (→ Ch. 9); ch05 `rotating_cylinder_flow`, `torque_per_length`, `lamb_oseen_circulation`
(→ Ch. 8), `diffusing_vortex_sheet`, `sheet_rollup` (→ Ch. 11), `column_relative_vorticity`, `column_over_slope`,
`relative_circulation_after_move` (→ Ch. 13), `ring_dynamics` (→ Ch. 14), `uniform_strain_vorticity` (→ Ch. 12); **ch06**
`elliptic_cylinder_flow(alpha=)`, `lift_per_span`, `joukowski_state` (→ Ch. 14), `cylinder_surface_speed`,
`half_body_surface_cp`, `sphere_surface_cp` (→ Ch. 9 outer flows), `sphere_motion` (added-mass dynamics → Ch. 7, 13,
16), `example_6_1` (→ Ch. 13 coastal eddy), `rayleigh_collapse_time` (→ Ch. 15/16). ch03 `cylinder_flow` and ch04
`accelerating_sphere_*` are now matched by `core.potential` (parity 1e-13, 1e-10); the old names stay. Note `ch04.G`,
`ch05.G` = `G_BOOK` = 9.81 while `core` defaults to `G0` = 9.80665; **ch06 ρ defaults: 1.2 (2-D forces) vs 1000 (§6.9,
Ex. 6.1)**.

### Explainer engine (`assets/viz_lib.js`)
`Viz.app` (tabs Walkthrough / Explore / Explain / Derivation / Equations / Code / Check; fit-to-window with density
levels and pagers; transport, modes, presets, status, terms, inspector, notes, selftest parity), `Plot`, `Viz.field`,
`Viz.num` (RK4, odeint, brentq, erf/erfc **~1.2e-7 only**, niceTicks), `Viz.work`, `Viz.three`, KaTeX with fallback;
`Viz.font`, `Viz.roundRect`, `Viz.text`, `Viz.card`, `Viz.fmtTime`, `Viz.tnum(keepTiny)`; pager slicing,
`--viz-stage-need`, `ev.viewId`/`ev.vizView`. **Nothing promoted after ch02–ch06** (the site-publisher was reading
`viz/` each time). Ranked ch06 candidates in `knowledge/viz_patterns.md`: **`Viz.cx` complex module** (4 files, needed by
Ch. 7, 11, 14), `arrowPx` (15 files), the narrow-view y-title gutter (7 local `gutterPlot`s), formatters (44/44),
per-mode hiding (Q4), contour `maxJump`/mask for branch cuts, field image, log axes, bars with faint/solid pairs,
Gauss–Legendre, LU/cond1; plus four open ch03 library bugs and quirks Q1–Q8. Reference: `templates/viz_example.html`.

## Explainer inventory
| chapter | slug | CORE idea | reusable stage pattern |
|---|---|---|---|
| ch01 | `continuum_averaging_volume` | C06 continuum (+Kn, D35) | molecules · value vs log box size with band · regime strip; click-a-sample inspector |
| ch01 | `viscosity_momentum_diffusion` | C12 Newton viscosity, ν | gap with tracers · profile with steady ghost and trail · wall stress |
| ch01 | `heat_work_paths` | C25, C35, C45 | piston · p–v · T–s linked on a path parameter; path vs state term bars |
| ch01 | `parcel_stability` | C50, C51, C54, C55 | column · profile with adiabat · ζ(t) with linear ghost; two-convention badge |
| ch01 | `buckingham_pi_machine` | C64, C67, C69 | variable chips → matrix with minor → groups; exact rational JS |
| ch02 | `rotation_of_axes` | C02, C03, C06 | plane with fixed arrow and turning axes · clickable C matrix · τ′(θ) |
| ch02 | `cauchy_traction_principal_axes` | C04, C05, C13 (D17 ★★★) | element with a cut, f split blue/rose · σn/τs vs φ · Mohr; term bars |
| ch02 | `strain_vs_rotation_split` | C12 | one clock, three squares (G, S, A); `expm2`; Frobenius bars |
| ch02 | `gauss_flux_box` | C14, C15 | heatmap + draggable box, signed face arrows · waterfall flux bars · limit view vs log h |
| ch02 | `stokes_circulation_loop` | C11, C16 | curl image + paddle wheel + draggable loop · unrolled u·t · Stokes failing on purpose |
| ch03 | `flow_lines_unsteady` | C03, C04 | streamline, particle trail and dye on one clock · answer-sheet view |
| ch03 | `material_derivative_probe` | C01, C02 | fixed probe + carried float · local/advective/total bars with measured ◇ |
| ch03 | `galilean_frames_cylinder` | C05 | observer slider lake → body · bars that trade places · u = U + u′ triangle |
| ch03 | `fluid_element_deformation` | C06–C09 | four G sliders deform a square + ring · rate curves with measured dots |
| ch03 | `spin_and_principal_axes` | C10–C12 | threads + paddle wheel · single-thread rates vs pair average · co-rotating observer |
| ch03 | `vortex_paddle_wheels` | C13, C14 | vortex with orbiting wheels and a draggable loop · u_θ, ω_z, Γ(r) · real-vortex table |
| ch03 | `reynolds_transport_cv` | C15 (D22 ★★★) | moving CV with signed swept band · budget waterfall + measured ◇ · dropped-term log–log |
| ch04 | `control_volume_budgets` | C01, C04 | five scenes on one budget object · face flux bars + measured ◇ |
| ch04 | `stream_function_spacing` | C03 | ψ at equal Δψ over a speed map · spacing = speed bars · axisymmetric mode |
| ch04 | `newtonian_stress_lab` | C07 (D09 ★★★) | G → S, R → τ grid coloured by source → rotatable plane traction · cube mode |
| ch04 | `navier_stokes_term_balance` | C08 | 7 exact solutions, click a point, term bars summing to 0, regime word |
| ch04 | `rotating_frame_coriolis` | C09 (D15 ★★★) | two observers on one clock · signed side toggle · split Coriolis arrow |
| ch04 | `which_bernoulli` | C05, C11, C12 | six flows with probes along/across streamlines · hypothesis decision table |
| ch04 | `viscous_dissipation_heating` | C10 | Couette gap: profile, ε, T(y, t) · in = stored + out bars |
| ch04 | `boussinesq_buoyancy` | C13 | rising blob with buoyancy-off ghost · kept vs dropped bars · validity chart |
| ch04 | `dynamic_similarity_models` | C15 | prototype and model on one t* clock · paired group bars · sphere collapse |
| ch05 | `vortex_tubes_cannot_end` | C01 (D01) | 3-D tube with a sliding section · flux/area/mean ω along the tube · Gauss bars; "broken ⚠" field |
| ch05 | `vortex_pressure_funnel` | C02 (D02, D03) | four vortices on one (r, z) section · needed = supplied force bars with FD slope |
| ch05 | `kelvin_material_loop` | C03, C05 (D04, D05, D08) | six flows + material vs fixed loop · ✓/✗ hypothesis table · measured ◇ dΓ/dt |
| ch05 | `baroclinic_torque` | C04 (D06, D07) | disc with tilted isopycnals · torque ◇ on the formula's sine · R² gap panel |
| ch05 | `vorticity_stretching_tilting` | C06, C10 (D09, D16–D18) | 3-D vortex line with stretching/tilting arrows · Burgers balance settling |
| ch05 | `biot_savart_filament` | C07, C08 (D10–D13) | "build the sum" transport · unrolled integrand · printed-sign toggle |
| ch05 | `vorticity_equation_rotating` | C09, C11 (D14, D15, D19, D20) | column over a ridge / ring moved poleward, conserved quantity flat in amber · 7-term budget |
| ch05 | `point_vortex_lab` | C12, C13 (D21, D22) | click-to-place vortices (guarded) · trajectories · invariants Q/Q(0); images |
| ch05 | `vortex_sheet_rollup` | C14 (D23) | row of filaments vs continuous jump · draggable circuit · roll-up with KH ghost |
| ch06 | `superposition_sandbox` | C03–C05 (D05–D07) | element kits · ψ-bars tip-to-tail onto the velocity arrow · labelled source-fluid tracers · "squeeze" limit button |
| ch06 | `cylinder_circulation_lift` | C06, C07 (D09–D11) | term bars of expanded (6.39) with ◇ ρUΓ · cross term in C_p(θ) vs "without it" · both Γ conventions in the status |
| ch06 | `vortex_wall_images` | C08 (D12, D13) | vortex + image + wall washed by p − p∞ · sensor trace with special-time markers · same-sign image leak |
| ch06 | `complex_potential_corners` | C09 (D14, D15) | dw/dz vs u + iv mirror arrows · two difference quotients · branch-cut status · log–log tip speed |
| ch06 | `blasius_kutta_contour` | C10 (D16–D18, two ★★★) | Laurent bars "size vs what survives" · force vs contour radius, rose "cuts the body" band · tilted-stream mode |
| ch06 | `conformal_joukowski` | C11 (D19–D21 ★★★) | ζ and z planes with the same particles · live β = α cross · ghost second root · wrong-branch streamlines in rose |
| ch06 | `laplace_relaxation` | C12 (D22, D23) | stencil stepping · matrix with b − Aψ · residual plot with method ghosts · honest refinement order |
| ch06 | `axial_singularity_bodies` | C14 (D26, D27) | bars vs converging moments · fit vs N with a cond band · matrix-row inspector · panel mode |
| ch06 | `added_mass_sphere` | C15 (D28–D31, two ★★★) | speed + acceleration parts = total · force bars · added mass three ways · release mode with a density table |

## Notation
See `knowledge/notation.md` (symbol register and sign traps). Conventions that matter everywhere: SI and kelvin inside
functions; **z up** (except ch06 §6.8, where z is the horizontal symmetry axis); **lapse rate Kundu Γ ≡ dT/dz in code,
meteorology −dT/dz shown alongside**; `p0` ≠ `p_ref`; q to / w on the system; gas constants per kmol; pressures absolute
unless `_gauge` (ideal flow: measured from hydrostatic); signed τ_xy = μ ∂u/∂y.
From ch02: passive C (x′ = Cᵀx); traction on the first index; book A:B = A_ij B_ji; G[i, j] = ∂u_i/∂x_j; Stokes n_c
into A; grid `[k, j, i]` with `h` in (x, y, z).
From ch03: R = G − Gᵀ (no ½), ω = ∇×u, spin ½ω; Galilean u′ = u − U; rotating ω′ = ω − 2Ω; RTT outward n, signed b·n;
field callables `u(x, t)`; spherical θ from +z.
From ch04: Cauchy's divergence on the first index; flux (u − b)·n; **Coriolis acceleration +2Ω × u′ (4.43) vs force
−2Ω × u′ (4.45)**; Stokes assumption μ_v = 0; 2-D ψ: u = ∂ψ/∂y (GFD often opposite); four primes; two g defaults.
From ch05: **ω is the vorticity** (tank at ω/2); **Γ = circulation [m²/s], γ = sheet strength [m/s] = u_below − u_above**;
counterclockwise positive; ∇ρ × ∇p in that order; (5.14) with +1/(4π); e_n = −Frenet N; `ellipk(m)` with m = k².
From ch06:
- **Γ counterclockwise in code; the book's (6.36)–(6.40), (6.52), (6.61)–(6.62), (6.68), Ex. 6.1 use a clockwise Γ**
  (circulation −Γ, L = +ρUΓ): functions take `Gamma_cw=`/`Gamma_ccw=`; explainers show both.
- **2-D doublet vector from sink to source** (cylinder d = −2πUa² e_x; (6.49) scalar d = dipole −d e_x); 3-D sphere
  d = 2πa³U; moving sphere d(t) = +2πa³u_s.
- **Three angle origins**: θ from +x (downstream) in polar formulas; from the upstream stagnation point in Fig. 6.10-style
  comparisons; θ_s from the sphere's velocity in (6.106).
- **w = φ + iψ, dw/dz = u − iv** (velocity = conjugate); w is the z-velocity in §6.9; ζ = Zhukhovsky plane, outside root.
- **Stokes ψ in m³/s**, flux 2πΔψ, (6.77) ≠ Laplace; (6.82) = r × App. B divergence (**not** a slip).
- ρ defaults 1.2 (2-D forces) / 1000 (§6.9, Ex. 6.1); pass ρ explicitly.
Symbols with several meanings so far: α, θ, T, h, R, σ/S, c, D, n, k, q, γ, λ, τ, ε, H (ch01); Γ, ω, φ, A, C, G, K, m, s,
t, a, b (ch02); Γ, ω, σ, φ, θ, α, a, b, h, primes, U (ch03); primes, σ, ε, γ, Φ, Ψ, Ω, h, H, M, R, f, F, B, C_p, λ, b, η,
ζ, G (ch04); Γ/γ, ω vs Ω, σ, ε, ζ, h, α, s, τ, κ, m, H, N, G, A (ch05); **Γ (cw/ccw), w (potential vs velocity), ζ
(Zhukhovsky plane), z (horizontal in §6.8), ξ (axial coordinate vs x − x_s), d (doublet), m (source), Q, M (added
mass), k (half-body ratio, line-sink and segment strengths), a, b, n (corner exponent), N (segments), ρ (spectral
radius), θ (three origins) (ch06)**.

## Teaching lessons
- **Depth tiered, coverage exhaustive** (A / B / C, derivations, cells): ch01 15/70/21, 12, 496 · ch02 16/60/18, 15, 405
  · ch03 15/52/12, 24, 413 · ch04 15/151/19, 30, 611 · ch05 14/69/11, 23, 496 · **ch06 15/116/29, 31 (257 steps), 564**.
- **Convention callouts** (lapse rate, "Which p_o?", passive/active, first/second index, A:B, R vs A, Coriolis term vs
  force, four primes, ψ sign, ω vs rate, Γ vs γ, ∇ρ × ∇p order, **Γ_cw vs Γ_ccw, doublet direction, three angle origins,
  w vs w**): state each convention, give the size of the difference in numbers, say which one the code uses and why —
  and (ch06) **carry one numeric case through the conventions table, the tiny example, the code and the explainer**
  (stagnation −9.16°/−170.84°, L = 24 N/m), with convention-named keywords (`Gamma_cw=`).
- **Book typos are taught, not hidden — by computing both versions** (ch02–ch05; **ch06 (6.61) coefficient, (6.104)
  bracket, (6.108) dφ, printed forms shown beside the correct ones in sympy cells**). **And a suspected slip is checked
  before it is taught**: (6.82) was flagged against a hybrid form the book never prints; it is r × App. B (retracted).
- **State what converges and what doesn't** (ch06): Example 6.2 at order 4/3 with its 270°-corner reason; axial
  strengths alternate while moments converge; circle panels exact, ellipse O(1/N²), off-body first order.
- **Running threads** make a chapter's logic visible (ch04 equations-vs-unknowns ledger; ch05 "life of vorticity"; ch06
  "solutions of Laplace add; any streamline can be a wall"). **Decision tables** sort results by hypotheses (ch04
  Bernoulli; ch05 Kelvin; ch06 `ideal_flow_applicability`).
- **Every key number by two or three independent routes** (ch06: lift by surface integral, far-field CV, Blasius on many
  contours and the residue; added mass by force, by energy (surface and volume) and by formula; Example 6.1 in closed form
  and by moving vortices) — also as a measured ◇ in the explainer.
- **★★★ sympy checks re-run the derivation's own construction** and keep what the book drops (ch05 D15; ch06 D29 direct
  ∂φ/∂t vs the chain rule for a concrete motion; D21 checks the step-11 velocity, added in round 1).
- **Every prose number is the printed number; every caption matches the plotted data** (ch06 M6: "dashed curves on three
  teal lines" that were not drawn; S4 halo "round-off" that was stencil error); *why* texts are checked for direction and
  sign claims on the whole body (D17: parallel **or antiparallel**) and for limit orders (D26).
- **Builder rendering hygiene** (ch06 round 1: 98 literal `\lvert` in prose, 43 `*…*` inside `\text{}`, glued `|x|word`):
  design TeX conventions must be converted, not copied, into non-maths text.
- **Primers come before the first use**; `knowledge/primers.md` lists 164 (P01–P164); ch07 starts at P165.
- **Climate hooks with numbers** (ch04 pole projectile, Ro; ch05 tornado, lock exchange 2.42 s⁻², ring of air; ch06 QG ψ
  inversion, frontogenesis deformation field ψ = 2Axy, bottom-pressure recorder under an eddy −25.3 Pa, overturning
  streamfunction geometry, bubbles/thermals ≤ 2g).
- Reusable derivation moves (ch01–ch06, 104; 20 from ch06) are tabulated in `knowledge/concept_map.md`.

## Global pitfalls confirmed in this project
**Machinery and environment**
- Extracted text garbles maths (`¼` = `=`, `ð…Þ` = parentheses, missing minus, ω → `u`, γ → g, θ → q, σ/τ → s, ν → n,
  ε → 3, ρ → r, Ω → U, ψ → j, χ → c, φ → 4, Φ → F, η → h, ζ → z, ∂ → v, ∇ → V) — read equations from rendered pages.
- **Windows: shell heredocs and `sed` mangle backslashes, quotes and `|`** (GNU sed treats `\|` as alternation) — write
  code, design, explainer and knowledge files with Write/Edit or Python; Python builder strings with LaTeX must be raw.
- **Single-backslash TeX in explainer JS strings** turns `\f`, `\t`, `\r` into control characters (ch06 E2 round-1 Must;
  not linted yet).
- **`MPLBACKEND=Agg` while executing a notebook strips inline figures** (coverage fails "no visual output"); scripts that
  call `plt.show()` block a headless run unless they switch to Agg (`--no-show` in ch06 scripts).
- Anaconda's `python3` kernelspec can shadow the venv's (use `fluidpy-venv`); JS `UIEvent.view` is read-only;
  `tools/shot.py` does not click/drag (a reviewer runs a click pass).
- `test_viz_library_inlined_and_template_lints` fails whenever `assets/viz_lib.js` changed or an explainer is mid-build;
  run `tools/viz_inline.py --all` before the merge gate.
- `tools/shot.py` `py:` parity expressions have **no builtins**; set tolerances near the achieved agreement.
- Parallel viz-builders use private scratch subfolders; runtimes measured while builders run browser audits are ~5×
  slower (ch06 notebook ~170 s alone).
- matplotlib 3.11: animations using `subplots_adjust`/`fig.text` disable `figure.constrained_layout.use` for the save.
- `coverage_check` matches ledger reminder rows by exact primer name → false warnings; it does not see literal TeX
  commands outside `$…$` (ch06 M1) or unbuilt explainers (warnings only).
- Book wording can hide in docstrings (`check_public.py` does not catch paraphrase): reviewers scan for it.
- Library: pager packs before the final `fit()`; tall ∫ clip at slice breaks; `Viz.num.erf/erfc` only ~1.2e-7;
  transport range fixed per app (ch06 E7 `syncRange`); no per-mode controls (Q4); `onChange` after `set` (Q5).
- **Phone legibility is the usual viz failure** (rotated y titles over minus signs, labels over labels, > 7 term rows in
  ~120 px, tick crowding in < 60 px views); ch06 copied a `gutterPlot` helper six times to avoid it.
- **Design documents contain errors; downstream agents compute, not copy** (ch05 three; ch06: "≈ 2 away from the
  corner" (really 4/3), "segments converge to two spikes" (they alternate), "panels O(1/N²) on the circle" (exact), the
  D21 fraction 1.0 (0.998 with the slit); the design's Part G errata override Parts A, B, F). Fix the design; builder
  `pf_sub` patches drift.
- Docstring contract numbers, validation lines and explainer texts go stale after fixes — re-grep them (ch06: 47
  "Validation (planned)" docstrings; E8 odd-N text).

**Mathematics → code**
- **A test comparing two of our own functions is not evidence for the book's convention.** Pin conventions to fields
  with known answers and **prove discrimination by patching in the wrong variant** (ch02 8/8; ch03 13/13 + 7; ch04 26/26;
  ch05 40/40 + 5; **ch06 30/30**) — and choose a field where the wrong variant differs.
- **Before calling an equation a book slip, compare the whole printed form with the reference form times every
  plausible factor** ((6.82) = r × App. B divergence).
- **A telescoping sum is not a conservation test** (ch06 S2: Σ Δψ across a section equals the BCs by construction; the
  real V4 invariant is the discrete cell circulation). Public tests pin a book example's geometry and BCs (S3).
- **Iterative searches need seeds at the problem's own length scale** (ch06 F1: Newton for the half-body's stagnation
  point failed for m ≲ 3 m²/s; seeds at ℓ·(¼…4), ℓ = m/2πU, plus deflation); test over the whole slider range. A double
  root is located only to √tol.
- **Complex branches are physics**: Zhukhovsky's outside root √(z − 2b)√(z + 2b), never numpy's principal
  √(z² − 4b²); cuts rotated into bodies (`cut_angle`); `arctan2`, not `arctan`.
- **Singular corners set the global convergence order** (Example 6.2: 4/3 everywhere, n = ⅔ in Azⁿ); report the order
  with its reason. **Stop iterative solvers on a residual far below the wanted accuracy** (residual and change ≈
  (1 − ρ) × error when ρ → 1).
- **Symmetric collocation at odd N is exactly singular** (A = −PAP); constant-source panels are exact on a circle;
  condition numbers are reported with every inverse method (axial method up to 10¹⁷).
- **Sibling functions share defaults and argument meanings**: one ρ per medium (ch06 S4: 20 vs 24 N/m), keyword-only
  when ambiguous (`sphere_potential_vector(a=)`, `Gamma_cw=`); return components *and* invariants in rotated problems
  (F_perp, F_par, stream angle).
- **Probe piecewise functions exactly at their joins** (ch05 F1); **check a field against its own stream function**
  (ch05 F2); test every scenario's full balance per direction (ch04).
- **Signed and ordered quantities stay signed/ordered**: b·n, (u − b)·n, wall stress, lapse-rate conventions, Stokes
  orientation, ∇ρ × ∇p, the sheet's u₂ − u₁, the (5.14) factor, **Γ_cw vs Γ_ccw, the doublet direction, dw/dz = u − iv**.
- **Book sign slips can cancel** ((5.14)); **or be harmless** ((6.61)'s 1/z² term drops out) — code the correct form, keep
  the printed one as an option (`printed_bracket=True`), test that it differs.
- Library parameter conventions: `ellipk(m)` with m = k²; `np.gradient` first order at edges; `np.sqrt`/`np.log`
  principal branches; pint °C offset; kmol vs mol; geopotential vs geometric altitude.
- Models have ranges; users have fingers: thin-core rings stop at a 3a gap; dragged vortices guarded; axial N capped at
  40 and even; separated-cylinder C_p band labelled qualitative.
- Rotating ≠ Galilean (ω′ = ω − 2Ω; Γ_a = Γ + 2Ω·A_vec); traction and Cauchy on the first index; sympy assumptions
  make different symbols, `.norm()` adds `Abs`, cache slow curvilinear simplifications.
- **Validation labels**: V4 only for real conservation laws; identities V1; **published closed forms reproduced
  identically are V1 form cross-checks, not V5**; "converged" only with an asserted order; property-table differences
  reported, not asserted.

## Benchmark inventory
| value / table | source (citation) | stored in | used by |
|---|---|---|---|
| k_B, N_A, R (exact) | CODATA 2018/2022, NIST CUU | `reference/ch01/constants.json` | ch01 `thermo` |
| USSA-1976 constants, Tables 1–2; sea-level c = 340.29 m/s | NASA-TM-X-74335; PDAS tables | `reference/ch01/ussa1976_*` | ch01, ch04 C02 |
| IAPWS σ(T), μ(T) | IAPWS R1-76(2014), R12-08 | `reference/ch01/` | ch01, ch04 capillary length |
| Jennings mean free path; Taylor blast K = 0.856; AMS Γ_d; EOS-80 check values | Tsalikis et al. 2024; Díaz arXiv:2009.05674; AMS Glossary; Fofonoff & Millard 1983 | `reference/ch01/benchmarks.json` | ch01 |
| Divergence theorem 8π/3; Levi-Civita identities; Mohr; rotation matrices; Stokes | Wikipedia; scipy docs | `reference/ch02/` | ch02 |
| Gaussian (Lamb–Oseen) vortex peak α = 1.256 | Canivete Cuissa & Steiner 2022, arXiv:2210.05223 | `reference/ch03/` | ch03 C14 |
| Sphere drag C_D(Re), 0.1–1e6 | F. A. Morrison (2016), Michigan Tech | `reference/ch04/` | ch04 C15 |
| WGS-84 a, 1/f, b, ω | Wikipedia "World Geodetic System" | `reference/ch04/` | ch04 C09, ch05 C09/C11 |
| Sharp-orifice C_c = 0.611; capillary length 2.71 mm; Prandtl ranges | Wikipedia | `reference/ch04/` | ch04 C05, C14, C15 |
| **Sphere added mass = ½ displaced mass** | Lamb, *Hydrodynamics* (1932) Ch. V §92 (primary, as cited); Wikipedia "Added mass" | `reference/ch06/` | ch06 C15 (V5) |
| **Rayleigh cavity collapse 0.91468 R₀√(ρ/Δp)** | Rayleigh, Phil. Mag. 34 (1917) 94, doi:10.1080/14786440808635681; Wikipedia "Rayleigh–Plesset" | `reference/ch06/` | ch06 exercise (V5) |
| Constant-strength source panels (method pedigree) | Hess & Smith (1967), doi:10.1016/0376-0421(67)90003-6 | `reference/ch06/SOURCES.md` | ch06 C14 |
| Forms (V1 cross-checks, not V5): Rankine, Lamb–Oseen, RTT, Leibniz (ch03); Bélanger, Tsiolkovsky, added mass, Taylor–Green (ch04); Kelvin ring, García–Haziot pairs, Burgers, Hill, segment law (ch05); **cylinder potential flow, Kutta–Joukowski, Blasius, Joukowsky W = W̃/(1 − 1/ζ²), McDonald (2015) sphere, Rankine half-body (ch06)** | Wikipedia pages; Commun. Math. Phys.; K. T. McDonald (Princeton) | `reference/ch03/`–`reference/ch06/` (`make_refs.py`, `SOURCES.md`) | ch03–ch06 |
Book-printed values (private, git-ignored): `tests/book_values_ch01.json` … `…_ch06.json` (ch06: §6.3 closed forms,
Fig. 6.8 angle 0.19 %, Example 6.1, Fig. 6.25's 24 values within the printed rounding, eight exercise answers ≤ 1e-9).
**Open**: a measured high-Re cylinder C_p (Roshko 1961, Achenbach 1968) would turn ch06 N31 from qualitative into V5.

## Validation summary per chapter
Counts are concept-map rows carrying each label (a row can carry several); tests per evidence level in parentheses.
| chapter | analytic | symbolic | converged | conserved | benchmark | book-value | qualitative | unverified | verdict |
|---|---|---|---|---|---|---|---|---|---|
| ch01 | 39 (V1 58) | 16 (V2 19) | 5 (V3 8) | 5 (V4 7) | 18 (V5 13) | 1 (V6 5) | 1 (`wavelength_to_rgb`) | 0 | **PASS** — 107 tests; 5 explainers (115 rows); 496 cells |
| ch02 | 43 (V1 41) | 21 (V2 13) | 11 (V3 15) | 2 (V4 1) | 9 (V5 6) | 8 (V6 4) | 0 | 0 | **PASS** — 92 tests (+ V7 12); 8/8 wrong variants; 5 explainers (127 rows); 405 cells |
| ch03 | 38 (V1 58) | 23 (V2 32) | 14 (V3 14) | 5 (V4 4) | 2 (V5 1) | 8 (V6 4) | 0 | 0 | **PASS** — 120 tests (+ V7 7); 13/13 (+7) wrong variants; 7 explainers (145 rows); 413 cells |
| ch04 | 41 (V1 83) | 35 (V2 41) | 17 (V3 9) | 13 (V4 6) | 11 (V5 6) | 5 (V6 3) | 0 | 0 | **PASS** — 151 tests (+ V7 3); 26/26 wrong variants; 9 explainers (255 rows); 611 cells |
| ch05 | 39 (V1 90) | 24 (V2 23) | 13 (V3 10) | 13 (V4 8) | 0 (V5 0) | 5 (V6 4) | 2 (ring past a 3a gap, roll-up shape; fenced) | 0 | **PASS** — 135 tests; 40/40 + 5 wrong variants; 9 explainers (268 rows); 496 cells |
| ch06 | 33 (V1 86) | 24 (V2 27) | 14 (V3 10) | 16 (V4 6) | 3 (V5 2) | 7 (V6 4) | 2 (separated-cylinder C_p band; axial-method sphere/airship targets beyond moments; both labelled) | 0 | **PASS** — 137 tests (+ V7 2); 30/30 wrong variants; review 1 Must (M1 retracted slip) + 8 Should resolved; every ★★/★★★ derivation re-derived in sympy; orders 2.01, 2.00, 4.00, 1.997, 2.057, 4/3 (corner), 1.06; 9 explainers PASS round 2 (291 rows, 2 864 views, 31 derivations matched); notebook PASS round 2 (564 cells, ~170 s) |
(ch01 V7: 19 tests; ch02: 12; ch03: 7; ch04: 3; ch05: inside V1; ch06: 2. Full suite at the ch06 merge gate: **750 tests**.)

## Open across chapters
- **Library pass (orchestrator, when nothing reads `viz/`)**: fix the four ch03 library bugs and quirks Q1–Q8 (Q4
  per-mode controls, Q5 `onChange` after `set`, fixed transport range); make the narrow-view y-title gutter an engine
  default; add `Viz.cx` (complex arithmetic with numpy branches) before Ch. 7/14; promote the ranked JS candidates
  (`knowledge/viz_patterns.md` ch06 list); then `tools/viz_inline.py --all` and `tools/shot.py --chapter ch01…ch06
  --quick` + `templates/viz_example.html --quick`.
- **Skill pass**: the lesson candidates for `interactive-viz`, `math-to-python` §7, `verify-implementation`,
  `teaching-style` and `colab-notebook` listed in `knowledge/viz_patterns.md` (ch02–ch06) are not yet in the skills
  (ch06 first: single-backslash TeX, "check a slip before teaching it", telescoping sums ≠ V4, seeds at the problem's
  length scale, one ρ per medium, one numeric case through every convention view).
- **Machinery TODO**: `viz_lint` rule for single-backslash TeX in JS strings; `run_notebook` refuses/unsets
  `MPLBACKEND=Agg`; script template switches to Agg when headless; `coverage_check` scans literal TeX outside `$…$`
  and understands ledger reminder rows; nbkit converter for design TeX → prose (spaces kept, no `*` inside `\text{}`);
  `shot.py` clicks/drags every view; `eq_refs` exempts `<meta>` and selftest names; `core.anim` saves inside a
  constrained-layout-off `rc_context`.
- Non-blocking explainer follow-ups: ch01–ch05 (their `knowledge/chNN.md` §9); ch06 E1 (Γ conventions, KaTeX Unicode,
  label overlap, 16 Explore pages), E2 (Explain hint), E3 (step 1 paused, units), E4 (step 2 with z³), E5 (Γ
  conventions, quiz 4), E7 (D23 tolerance, labels), E8 (axis hook), E9 (subscript, paused step 6) — `knowledge/ch06.md`
  §9. Backups not built: ch04 `kinematic_free_surface`, ch05 `vortex_rings`, ch06 `flow_net_sources_vortices`.
- ch06 notebook: glued "|x|word" spaces, one inset over a curve; ch06 O4 (off-body panel order, hypothesis), O6 (no
  measured cylinder C_p); ch04 O5 Prandtl numbers (−1.5 %, −2.2 %); ch05 O2 tube-flux route 1.8e-4.
