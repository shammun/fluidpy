# CUMULATIVE knowledge — fluidpy
(Rewritten by the knowledge-keeper after every chapter. Every agent reads this first. Last rewrite: after ch07,
2026-09-24, commit 3e32618.)

## Project rules in force (added during ch01–ch03; sharpened since)
- **5–10 explainers per chapter, as many as the CORE ideas need** (`book.yaml → project.min/max_explainers_per_chapter`);
  ch03 has 7, ch04–ch07 9 each.
- **Every book equation is shown in full next to its number** — notebook prose, derivation steps, traps, recaps and
  explainer tour/Explain/Derivation/quiz/notes/status text **and `<meta>` strings** (ch04 E6, ch06 E6). Enforced by
  `tools/coverage_check.py` check 8 (markdown) and `tools/eq_refs.py` (explainers); reviewers still scan by eye.
- **`coverage_check` check 9**: no control characters in cells (a LaTeX backslash eaten by a non-raw string).
- **Rule 9 (public repo) in practice (ch07)**: a value we compute can round to an exercise's printed answer. Numbers
  that coincide with a book value are built from the function (never typed) and printed with **4 significant figures**
  (ch07: the minimum group speed as 17.76 cm/s); the book's value lives only in `tests/book_values_chNN.json`.
- **Pager and audit**: an item taller than a page is sliced at natural breaks; `tools/shot.py` pages through every page
  of every pager at every size.

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
ch06 Ideal flow — "irrotational + incompressible ⇒ Laplace; solutions add; any streamline can be a wall"
     (6.1) (D01) [C01]; ω_z = −∇²ψ, ln r as a delta, φ ⟂ ψ (D02–D04) [C02]; superposition, wall = streamline (D05) [C03];
     doublet (D06), half-body (D07, D08), cylinder C_p = 1 − 4 sin²θ, D = 0 (D09), L = ρUΓ (D10, D11), images and
     Example 6.1 (D12, D13) [C04–C08]; w = φ + iψ, dw/dz = u − iv, corners Azⁿ (D14, D15) [C09]; Blasius (6.60) and
     Kutta–Zhukhovsky (D16–D18 ★★★) [C10]; conformal maps, outside root (D19–D21 ★★★) [C11]; Laplace on a grid, order 4/3
     at a 270° corner (D22, D23) [C12]; Stokes ψ, sphere (D24, D25) [C13]; airship, axial method (D26, D27) [C14];
     accelerating sphere, added mass ½ displaced mass (D28–D31 ★★★) [C15]
ch07 Gravity waves (done) — "plane wave in, boundary conditions → algebra, ω(k) out; read everything from ω(k)"
  vocabulary: phase, c = ω/k (D01), K, trace speeds ≥ c, Doppler ω₀ = ω + U·K (7.9) [C01]
  linear free surface: ∇²φ = 0, kinematic (7.16) with f = z − η (D02), Taylor-moved to z = 0 at relative cost ka
     (7.18), (7.21) (D03, D04) [C02] → separable φ (7.26) (D05) → ω² = gk tanh kH (7.28) (D06) [C03]
     → c(λ) dispersive, deep √(g/k), shallow √(gH), p′ ∝ cosh k(z+H)/cosh kH, hydrostatic when shallow (D07–D09) [C04]
     → closed clockwise ellipses (D10, D11) [C05]; E = ½ρga², E_k = E_p, flux (7.44) (D12–D14) [C06]
  σ: g → g + σk²/ρ, c_min = (4gσ/ρ)^¼ (D15, D16) [C07]; standing waves and seiches (D17, D18) [C08]
  groups: beats (D19), packet envelope at c_g (D20 ★★★), c_g = (c/2)(1 + 2kH/sinh 2kH), F = E c_g (D21) [C09]
     → crest conservation, ω constant along rays, Snell refraction (D22–D24) [C10]
  nonlinear: jump H₂/H₁ = ½(−1 + √(1 + 8Fr₁²)), energy loss ⇒ Fr₁ ≥ 1 (D25, D26); KdV soliton [C11];
     Stokes γ = 1 (consistent 3rd order), drift ū_L = a²ωke^{2kz₀} (D27) [C12]
  layers: ε√(gk) and the vortex sheet (D28) [C13]; barotropic/baroclinic, √(g′H), g′ = gΔρ/ρ₂ (D29, D30 ★★★, D31) [C14]
  stratified: linear Boussinesq (D32) → w-equation (7.134) (D33 ★★★) → ω = N cos θ (D34) [C15];
     K·u = 0, c ⟂ c_g, beams at arccos(ω/N), F = c_g E (D35–D37) [C16]

next: ch08 laminar flows — Couette/Poiseuille/Stokes problems from ch04 `exact_solution` + `navier_stokes` residuals,
      ch03 `pipe_profile`, ch05 `rotating_cylinder_flow`/`lamb_oseen_circulation`, ch06 moving sphere vs Stokes drag,
      ch07 `viscous_decay` (a₀e^{−2νk²t}) and the overflow-safe `core.waves.cosh_over_sinh` for oscillating layers.
```
The book ahead: Ch8 laminar → Ch9 boundary layers → Ch10 CFD → Ch11 instability → Ch12 turbulence → **Ch13 GFD (uses
4, 5, 7, 8, 11, 12)** → Ch14 aerodynamics · Ch15 compressible · Ch16.
Where earlier chapters feed in: ch01 N², θ → Ch. 11, 13; τ, ν → Ch. 8, 9; isentropic gas → Ch. 15; FTCS → Ch. 8, 10.
ch02 Gauss → every CV; Stokes → Ch. 13 PV; invariants → Ch. 12. ch03 vortices → Ch. 8, 13, 14. ch04 **NS residual/term
tools → every exact solution**; CV budgets → Ch. 9, 14, 15; ε → Ch. 12; rotating frame, Boussinesq, Ri, Ro → **Ch. 13**.
ch05 vorticity budget → Ch. 8, 10, 12; sheet roll-up → Ch. 11; **(5.30), (ζ + f)/h → Ch. 13**; filaments → Ch. 14. ch06
`core.potential` → **Ch. 9** outer flows (U_e(x), wedge Azⁿ → Falkner–Skan) and **Ch. 14** (lift, Kutta, circle theorem);
`core.conformal` → Ch. 14; `core.laplace_solvers` → **Ch. 10**, Ch. 13 PV inversion; `core.panels` → Ch. 14. **ch07**
`core.waves` → **Ch. 13** (√(gH), Poincaré/Kelvin/Rossby dispersion and c_g by complex step, reduced gravity and the
Rossby radius √(g′H)/f, inertia–gravity waves from `internal_wave_omega` + f, WKB rays, geostrophic adjustment by
`linear_evolve`), **Ch. 11** (interface ε√(gk) + shear = Kelvin–Helmholtz, σk²/ρ cut-off, Rayleigh–Taylor flag),
**Ch. 15** (jump ↔ normal shock, simple wave ↔ characteristics, acoustic rays), Ch. 10 (`kdv_solve` IF-RK4 spectral).

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
| hydrostatics, USSA-1976, scale height, buoyancy; N² family, `classify_stability`, `parcel_displacement`, lapse-rate conventions, θ; FTCS diffusion; Π groups | `statics`, `stratification`, `diffusion`, `dimensional` | ch01 | analytic, symbolic, converged, benchmark |
| index notation, `tensors`, grids `[k, j, i]` + `operators`, `fields`, `integral_theorems` | ch02 modules | (2.1)–(2.36) | analytic, symbolic, converged, conserved, benchmark |
| `kinematics` (field callables `u(x, t)`, material derivative, flow lines, `pathline`), `coords`, `transport` (Leibniz, RTT); `vortices` | ch03/ch05 modules | (3.1)–(3.35), §5.4 | analytic, symbolic, converged, conserved, benchmark |
| `conservation` budgets on any `ControlVolume`; `constitutive`; `navier_stokes` residuals and term splits (8 exact solutions); `curvilinear`; `streamfunction`; `rotating`; `bernoulli`; `interfaces` (kinematic residual, Laplace jump, capillary length); `similarity` (`Scales`, 20+ groups, `reduced_gravity` with ρ₁) | ch04 modules | (4.5)–(4.119) | analytic, symbolic, converged, conserved, benchmark |
| `vorticity` (loops, material circulation, Kelvin terms, tubes, **`vorticity_budget`**), `biot_savart` (Green, filaments, rings, point vortices, sheets) | ch05 modules | (5.3)–(5.33) | analytic, symbolic, converged, conserved |
| `potential` (elements, `Flow`, `FunctionFlow`, bodies, images, `blasius_force`, axisymmetric, moving sphere, added mass), `conformal` (Zhukhovsky, outside root, `mapped_flow`), `laplace_solvers` (masked 5-point, Jacobi/GS/SOR, sparse, `solve_poisson`), `panels` (axial singularities, source panels) | ch06 modules | (6.1)–(6.109) | analytic, symbolic, converged, conserved, benchmark |
| **`waves`** vocabulary: `sinusoid`, `wave_parameters`, `crest_positions`, `plane_wave`, `phase_velocity_vector`, `trace_velocities`, `doppler_frequency`, `real_field` | `waves` (ch07) | (7.1)–(7.9) | analytic, symbolic |
| **`waves`** dispersion: `omega_gravity`, `omega_capillary_gravity`, `phase_speed` (\|c\|), `group_velocity` (**signed**), `period_from_wavelength`, `wavenumber_from_omega` (brentq), `wavelength_from_period`, `fenton_mckee_kh`, `guo_kh`, `depth_regime` (label + errors), `capillary_minimum`, `min_group_velocity`; overflow-safe `cosh_over_sinh`, `sinh_over_sinh`, `cosh_over_cosh`, `depth_profiles` | `waves` | (7.28)–(7.31), (7.45)–(7.60), (7.69) | analytic, symbolic, benchmark |
| **`waves`** any ω: `group_velocity_numeric` (complex step, central fallback), `group_velocity_vector` (∇_K ω), `beat_wave`, `gaussian_packet(order=1\|2)`, `linear_evolve` (FFT), `envelope`, `linear_evolve_2d`, `ray_trace` (Hamiltonian), `wave_energy_density`, `viscous_decay` | `waves` | (7.66)–(7.79), (7.42) | analytic, symbolic, converged, conserved |
| **`waves`** layers and stratification: `eps2_density`, `interface_omega` (NaN + warning if ρ₁ > ρ₂), `two_layer_free_surface_omega`, `reduced_gravity_book(ref="lower"\|"upper")`, `two_layer_long_wave_speed`, `internal_wave_omega` (N\|k\|/K), `beam_angle`, `internal_wave_velocities` (sign-safe) | `waves` | (7.95), (7.110)–(7.117), (7.137)–(7.146) | analytic, symbolic |
Chapter-only code (move to `core/` on second use): ch03 `pipe_profile` (→ Ch. 8); ch04 `couette_heating`,
`two_fluid_couette`, `stokes_first_problem` (→ Ch. 8), `wake_drag_per_span` (→ Ch. 9), `linear_wave_surface` (**overflows
for kH ≳ 710; superseded by `ch07.wave_fields`**); ch05 `rotating_cylinder_flow`, `lamb_oseen_circulation` (→ Ch. 8),
`sheet_rollup` (→ Ch. 11), `column_relative_vorticity` (→ Ch. 13), `ring_dynamics` (→ Ch. 14); ch06
`elliptic_cylinder_flow`, `lift_per_span` (→ Ch. 14), `cylinder_surface_speed`, `half_body_surface_cp` (→ Ch. 9),
`sphere_motion` (→ Ch. 13, 16); **ch07** `wave_fields`, `free_surface_residuals`, `particle_path`/`stokes_drift`/`dyed_line`
(→ Ch. 13), `hydraulic_jump`, `simple_wave_evolve` (→ Ch. 15), `kdv_solve` (IF-RK4 pseudo-spectral → Ch. 10, 15),
`seiche_modes`/`basin_modes`, `two_layer_*`, `internal_wave_fields`, `internal_wave_energy`, `boussinesq_linear_sympy`
(→ Ch. 13). **Default g**: `ch04.G`, `ch05.G`, `ch07` and `core.waves` use `G_BOOK` = 9.81, the rest of `core` G0 =
9.80665. **Default ρ**: ch06 1.2 (2-D forces) vs 1000 (§6.9); ch07 1000; teaching water σ = 0.07274 N/m, ρ = 998.2.

### Explainer engine (`assets/viz_lib.js`)
`Viz.app` (tabs Walkthrough / Explore / Explain / Derivation / Equations / Code / Check; fit-to-window with density
levels and pagers; transport, modes, presets, status, terms, inspector, notes, selftest parity), `Plot`, `Viz.field`,
`Viz.num` (RK4, odeint, brentq, erf/erfc **~1.2e-7 only**, niceTicks), `Viz.work`, `Viz.three`, KaTeX with fallback;
`Viz.font`, `Viz.roundRect`, `Viz.text`, `Viz.card`, `Viz.fmtTime`, `Viz.tnum(keepTiny)`; pager slicing,
`--viz-stage-need`, `ev.viewId`/`ev.vizView`. **Nothing promoted after ch02–ch07** (the site-publisher was reading `viz/`
each time). Ranked ch07 list in `knowledge/viz_patterns.md`: pager-label nowrap CSS (5 files), `arrowPx`/`dotPx` (21 + 5
files), gutter + square plots (12 + 1), formatters incl. length/speed (all 53), **`Viz.waves` dispersion mirror before
Ch. 13**, per-mode hiding (19 workarounds), log axes, `Viz.cx`, Gauss–Legendre. Reference: `templates/viz_example.html`.

## Explainer inventory
| chapter | slug | CORE idea | reusable stage pattern |
|---|---|---|---|
| ch01 | `continuum_averaging_volume` | C06 continuum | molecules · value vs log box size · regime strip; sample inspector |
| ch01 | `viscosity_momentum_diffusion` | C12 | gap with tracers · profile with steady ghost · wall stress |
| ch01 | `heat_work_paths` | C25, C35, C45 | piston · p–v · T–s on one path parameter; path vs state bars |
| ch01 | `parcel_stability` | C50–C55 | column · profile with adiabat · ζ(t) with linear ghost; two-convention badge |
| ch01 | `buckingham_pi_machine` | C64–C69 | chips → matrix with minor → groups; exact rational JS |
| ch02 | `rotation_of_axes` | C02, C03, C06 | fixed arrow and turning axes · clickable C matrix · τ′(θ) |
| ch02 | `cauchy_traction_principal_axes` | C04, C05, C13 | element with a cut · σn/τs vs φ · Mohr; term bars |
| ch02 | `strain_vs_rotation_split` | C12 | one clock, three squares (G, S, A); `expm2` |
| ch02 | `gauss_flux_box` | C14, C15 | heatmap + draggable box · waterfall flux bars · limit vs log h |
| ch02 | `stokes_circulation_loop` | C11, C16 | curl image + paddle wheel + loop · unrolled u·t · Stokes failing on purpose |
| ch03 | `flow_lines_unsteady` | C03, C04 | streamline, trail and dye on one clock · answer sheet |
| ch03 | `material_derivative_probe` | C01, C02 | probe + float · local/advective/total bars with measured ◇ |
| ch03 | `galilean_frames_cylinder` | C05 | observer slider · bars that trade places |
| ch03 | `fluid_element_deformation` | C06–C09 | G sliders deform a square + ring · measured rates |
| ch03 | `spin_and_principal_axes` | C10–C12 | threads + paddle wheel · co-rotating observer |
| ch03 | `vortex_paddle_wheels` | C13, C14 | vortex with wheels and a loop · u_θ, ω_z, Γ(r) |
| ch03 | `reynolds_transport_cv` | C15 (D22 ★★★) | moving CV with swept band · budget waterfall · dropped-term log–log |
| ch04 | `control_volume_budgets` | C01, C04 | five scenes on one budget object · face flux bars + ◇ |
| ch04 | `stream_function_spacing` | C03 | ψ at equal Δψ · spacing = speed · axisymmetric mode |
| ch04 | `newtonian_stress_lab` | C07 (D09 ★★★) | G → S, R → τ grid → rotatable plane traction · cube mode |
| ch04 | `navier_stokes_term_balance` | C08 | 7 exact solutions, click a point, term bars summing to 0 |
| ch04 | `rotating_frame_coriolis` | C09 (D15 ★★★) | two observers · signed side toggle · split Coriolis arrow |
| ch04 | `which_bernoulli` | C05, C11, C12 | six flows with probes · hypothesis decision table |
| ch04 | `viscous_dissipation_heating` | C10 | Couette: profile, ε, T(y, t) · in = stored + out bars |
| ch04 | `boussinesq_buoyancy` | C13 | rising blob with ghost · kept vs dropped bars · validity chart |
| ch04 | `dynamic_similarity_models` | C15 | prototype and model on one t* clock · paired group bars |
| ch05 | `vortex_tubes_cannot_end` | C01 | 3-D tube with sliding section · Gauss bars; "broken" field |
| ch05 | `vortex_pressure_funnel` | C02 | four vortices on one section · needed = supplied force bars |
| ch05 | `kelvin_material_loop` | C03, C05 | material vs fixed loop · ✓/✗ hypothesis table · measured ◇ |
| ch05 | `baroclinic_torque` | C04 | tilted isopycnals · torque ◇ on the sine · R² gap panel |
| ch05 | `vorticity_stretching_tilting` | C06, C10 | 3-D vortex line with arrows · Burgers balance settling |
| ch05 | `biot_savart_filament` | C07, C08 | "build the sum" transport · unrolled integrand · printed-sign toggle |
| ch05 | `vorticity_equation_rotating` | C09, C11 | conserved quantity flat in amber · 7-term budget |
| ch05 | `point_vortex_lab` | C12, C13 | click-to-place vortices (guarded) · invariants Q/Q(0); images |
| ch05 | `vortex_sheet_rollup` | C14 | filaments vs jump · draggable circuit · roll-up with KH ghost |
| ch06 | `superposition_sandbox` | C03–C05 | element kits · ψ-bars tip-to-tail · labelled tracers · squeeze limit |
| ch06 | `cylinder_circulation_lift` | C06, C07 | term bars of (6.39) with ◇ ρUΓ · both Γ conventions in the status |
| ch06 | `vortex_wall_images` | C08 | vortex + image + wall washed by p − p∞ · sensor trace with markers |
| ch06 | `complex_potential_corners` | C09 | dw/dz vs u + iv mirror arrows · branch-cut status · log–log tip speed |
| ch06 | `blasius_kutta_contour` | C10 (two ★★★) | Laurent bars "size vs what survives" · force vs contour radius |
| ch06 | `conformal_joukowski` | C11 (D21 ★★★) | two planes, same particles · ghost root · wrong-branch lines in rose |
| ch06 | `laplace_relaxation` | C12 | stencil stepping · matrix with b − Aψ · residual ghosts · honest order |
| ch06 | `axial_singularity_bodies` | C14 | moments converge, bars don't · cond band · matrix-row inspector |
| ch06 | `added_mass_sphere` | C15 (two ★★★) | speed + acceleration parts = total · added mass three ways |
| ch07 | `dispersion_relation` | C03, C04 (D01, D05–D09) | tank + c(λ) with ghost deep/shallow limits and shaded regime bands + depth profile; sea-bed banner on phones |
| ch07 | `particle_orbits` | C05, C12 (D10, D11, D27) | linear vs exact path lines · measured vs formula drift · dyed line |
| ch07 | `capillary_gravity_waves` | C07 (D15, D16) | restoring pushes at the crest · term bars over a c²_min floor · liquids as modes · computed 4-s.f. constants |
| ch07 | `seiche_standing_waves` | C08 (D17, D18) | right/left/sum ghosts · period-vs-mode ladder with deep/shallow curves |
| ch07 | `group_velocity_packets` | C09 (D19–D21, D20 ★★★) | chord vs tangent on ω(k) · "group at c" ghost · ∫E dx strip · pond calm disc |
| ch07 | `wave_rays_refraction` | C10 (D22–D24) | in-page Hamiltonian rays with ω drift · Snell inspector · one code block per geometry |
| ch07 | `hydraulic_jump` | C11 (D25, D26) | momentum/energy term bars · "forbidden" preset · bore and solitary modes |
| ch07 | `two_layer_modes` | C13, C14 (D28–D31, two ★★★) | surface and interface magnified separately · both roots with limits · "which g′" |
| ch07 | `internal_wave_beams` | C15, C16 (D33 ★★★, D34–D36) | tank + square K plane side by side on phones · sign-safe k < 0 · c ⟂ c_g |

## Notation
See `knowledge/notation.md` (sign traps first, then the register). Conventions that matter everywhere: SI and kelvin
inside functions; **z up** (except ch06 §6.8, z horizontal); **lapse rate Kundu Γ ≡ dT/dz in code, meteorology shown**;
`p0` ≠ `p_ref`; pressures absolute unless `_gauge` (ideal flow: from hydrostatic; ch07: **gauge**); signed τ_xy = μ ∂u/∂y.
ch02: passive C; traction on the first index; A:B = A_ij B_ji; G[i, j] = ∂u_i/∂x_j; grid `[k, j, i]`. ch03: R = G − Gᵀ,
spin ½ω; u′ = u − U; ω′ = ω − 2Ω; RTT outward n. ch04: Coriolis acceleration +2Ω × u′ vs force −2Ω × u′; μ_v = 0;
**2-D ψ: u = ∂ψ/∂y (ch07: u = ∂ψ/∂z; GFD often opposite)**; four primes; two g defaults. ch05: ω vorticity; Γ circulation,
γ = u_below − u_above; ∇ρ × ∇p order; `ellipk(m)`, m = k². ch06: **Γ counterclockwise in code, `Gamma_cw=`/`Gamma_ccw=`**
where the book is clockwise; doublet from sink to source; three angle origins; w = φ + iψ, dw/dz = u − iv; outside
Zhukhovsky root; Stokes ψ [m³/s].
**ch07 (changes):** **ω = angular frequency** (was vorticity), ω ≥ 0 with direction in sgn k, `phase_speed` |c| vs
signed `group_velocity`; **η = elevation** (ch04: level-set function; call `kinematic_bc_residual` with z − η); p′ = p + ρgz
(§7.2) vs p − p̄(z) (§7.8); ζ = particle excursion **and** interface displacement; θ = phase **and** K's angle above the
horizontal (= the beam's angle from the vertical); **g′ = g(ρ₂ − ρ₁)/ρ₂** (ch04 ρ₁; `reduced_gravity_book(ref=)`);
printed (7.138)/(7.145) assume k > 0 (code ∇_K ω); **complex amplitudes with Re dropped, ½Re(AB*) for means** (in ch06
the imaginary part was ψ, here a phase shift); E per horizontal area (surface), per volume (internal), per mass (jump).
Symbols with several meanings (see the register): α, θ, T, h, R, σ, c, D, n, k, q, γ, λ, τ, ε, H (ch01); Γ, ω, φ, A, C,
G, K, m, s, t, a, b (ch02–ch03); primes, Φ, Ψ, Ω, M, f, F, B, C_p, η, ζ (ch04); Γ/γ, κ, N (ch05); w, ζ, z, ξ, d, m, Q, M,
k, n, ρ, θ (ch06); **ω, η, ζ, θ, g′, k sign, m, ν, T, σ/σ_x, γ, ε², b, E, F, f, H, Q, w₀/ŵ (ch07)**.

## Teaching lessons
- **Depth tiered, coverage exhaustive** (A / B / C, derivations, cells): ch01 15/70/21, 12, 496 · ch02 16/60/18, 15, 405
  · ch03 15/52/12, 24, 413 · ch04 15/151/19, 30, 611 · ch05 14/69/11, 23, 496 · ch06 15/116/29, 31, 564 · **ch07
  16/193/27, 37 (321 steps), 617**.
- **Convention callouts** ("which p_o?", passive/active, Coriolis term vs force, Γ_cw vs Γ_ccw, **which ω, which η, which
  g′, which p′**): state each convention, give the size of the difference in numbers, say which the code uses — and
  carry **one numeric case through table → tiny example → code → explainer** (ch06 L = 24 N/m; ch07 g′ 0.01910 vs 0.01914).
- **Book slips are taught by computing both versions**, shown as ghosts that visibly fail (ch07: the (7.66) envelope that
  does not move, (7.105) residuals that do not vanish); **check a suspected slip before teaching it** ((6.82) was not one;
  Stokes γ = 1 is right — the literal truncated exercise set-up gives 3/8).
- **Make an approximation measurable** (ch01 linear vs full; ch07 residual scan of the linearised free-surface conditions:
  relative slope 1, absolute 2) and **state what converges and what doesn't** (ch06 4/3; ch07 Stokes drift order 0.98).
- **Every key number by two independent routes** (ch06 lift, added mass; ch07 (7.65) by closed form and an FD sloshing
  eigenproblem, the drift by formula and exact path lines, c_g analytic and by complex step).
- **Every prose number is printed by a cell; captions count what they claim** (ch06 M6; ch07 island rays reaching the
  lee printed by the cell); **thresholds and crossovers are computed, not quoted** (ch07 Ursell 4π²/9, regime errors
  1.815 % / 3.039 %).
- **Honest forward use**: "take it on trust for now (derived in C03–C04)" (ch07 §7.1); scaling from the forcing (the
  surface rises at aω ⇒ velocities aω).
- **Primers come before first use**; idiom primers whose two-line demo output *is* the explanation (ch07 P179–P184);
  `knowledge/primers.md` lists 184 (P01–P184); ch08 starts at **P185**.
- **Builder rendering hygiene**: design TeX must be converted, not copied (ch06 `\lvert`; ch07 raw `e^{kz}` in *why* lines
  → `tidy_raw_tex` + a build-time guard).
- **Climate hooks with numbers** (ch04 Ro; ch05 tornado, lock exchange; ch06 ψ inversion, bottom-pressure recorders;
  **ch07 38 kW/m swell, 198 m/s tsunami, 0.99 m/s thermocline wave vs 22 m/s surface wave, Rossby radius √(g′H)/f, Stokes
  drift and Langmuir cells, internal tides and lee waves**).
- Reusable derivation moves (ch01–ch07, 120; 16 from ch07) are tabulated in `knowledge/concept_map.md`.

## Global pitfalls confirmed in this project
**Machinery and environment**
- Extracted text garbles maths (`¼` = `=`, `ð…Þ` = parentheses, missing minus — ch07's jump energy E₂ − E₁ and e^{−2kH} —,
  ω → `u`, γ → g, θ → q, σ/τ → s, ν → n, ε → 3, ρ → r, Ω → U, ψ → j, φ → 4, η → h, ζ → z, ∂ → v, ∇ → V): read page images.
- **Windows: shell heredocs and `sed` mangle backslashes, quotes and `|`** — write code, explainers and knowledge files
  with Write/Edit or Python files; builder strings with LaTeX must be raw; **single-backslash TeX in JS strings** turns
  `\f`, `\t`, `\r` into control characters (ch06 E2; not linted yet).
- `MPLBACKEND=Agg` while executing a notebook strips inline figures; scripts calling `plt.show()` block headless runs
  (`--no-show`). Anaconda's kernelspec can shadow the venv's (`fluidpy-venv`). JS `UIEvent.view` is read-only.
- `tools/shot.py` does not click/drag (a reviewer runs a click pass); **its `py:` expressions have no builtins and no
  `np` — bind functions via lambda default args** (`lambda X, th=ch07.np.tanh: …`); set parity tolerances near the
  achieved agreement.
- `test_viz_library_inlined_and_template_lints` fails whenever `assets/viz_lib.js` changed or an explainer is mid-build;
  run `tools/viz_inline.py --all` before the merge gate. Runtimes measured under parallel browser audits are ~5× slower
  (ch07 notebook ~189 s alone, KaTeX-late pager clips appear intermittently: keep one line of slack).
- matplotlib 3.11: animations using `subplots_adjust`/`fig.text` disable constrained layout for the save.
- `coverage_check` matches ledger reminders by exact primer name; it does not see literal TeX outside `$…$`.
- Book wording can hide in docstrings and **exercise answers can hide in computed numbers** (`check_public.py` catches
  neither): reviewers scan against the private JSON.
- Library: pager packs before the final `fit()`; `Viz.num.erf` ~1.2e-7; fixed transport range; no per-mode controls (Q4,
  19 local workarounds); `onChange` after `set` (Q5). **Phone legibility is the usual viz failure** (rotated y titles,
  labels over labels, > 7 term rows in ~120 px, squashed strips — ch07 E9 needed two square views side by side).
- **Design and analysis documents contain errors; downstream agents compute, not copy** (ch05 three; ch06 four; ch07:
  "(ka)² × aω" was relative O(ka), the 2006 note's 1.5 %/0.7 % were not the primary bounds, γ = 1 was attributed to a set-up
  that gives 3/8). Docstring numbers and "Validation (planned)" lines go stale — re-grep them (ch06 47, ch07 108).

**Mathematics → code**
- **A test comparing two of our own functions is not evidence for the book's convention.** Pin conventions to fields
  with known answers and **prove discrimination by planting the wrong variant** (ch02 8/8 … ch06 30/30, **ch07 41/41**) —
  on a case where the variant differs (ch07: k < 0 for the signed c_g).
- **Complex-step derivatives need analytic continuations of every non-analytic piece** (ch07 M1: `_abs_k` returned a
  complex k unchanged, c_g = 3e20 m/s for k < 0); keep a speed (|c|) and a signed derivative (dω/dk) as different
  functions; `np.abs` is not a continuation.
- **Printed formulas can assume a sign** ((7.138), (7.145) need k > 0; the book's own figure has k < 0): code the gradient
  form ∇_K ω, test all quadrants, angles from cos θ = |k|/K.
- **Before calling an equation a book slip, compare the whole printed form with the reference form times every
  plausible factor** ((6.82)); a literal exercise set-up may give a different coefficient than the published result
  (Stokes γ) — check the expansion's consistency order.
- **Overflow**: cosh/sinh ratios in the e^{kz}(1 ± e^{−2k(z+H)})/(1 ∓ e^{−2kH}) form (finite at kH = 500); x/sinh x via
  expm1; numpy's complex tanh raises a spurious overflow flag (`np.errstate`); ch04's plain form NaNs at kH ≳ 710.
- **Inverse relations by a physics bracket** (deep and shallow roots bound k; brentq; residual asserted) and **iterative
  searches seeded at the problem's own length scale** (ch06 F1); a double root is located only to √tol.
- **Complex branches are physics** (Zhukhovsky outside root, cuts inside bodies, `arctan2`); **complex amplitudes need
  real parts before products** (½Re(AB*); ρ′ 90° from w).
- **Singular corners set the global convergence order** (Example 6.2: 4/3); **asymptotic orders are stated relative or
  absolute** (ch07 O(ka) relative = slope 2 absolute at fixed k); stop iterative solvers on a residual far below the wanted
  accuracy; symmetric collocation at odd N is singular; condition numbers are reported with every inverse method.
- **A telescoping sum is not a conservation test** (ch06 S2); the genuine V4 invariants in ch07 are ∫η² of a packet, KdV
  mass/momentum/Hamiltonian, ω along a ray, jump momentum and the internal-wave energy budget.
- **Sibling functions share defaults and argument meanings** — or are documented and called by keyword (ch06 one ρ per
  medium; ch07 `omega_capillary_gravity(k, H, sigma, rho, g)` vs `phase_speed(k, H, g, sigma, rho)` with different σ
  defaults, kept for API stability); convention-named keywords (`Gamma_cw=`, `ref="lower"`).
- **Two quantities with one name**: ch04 vs ch07 g′ (ρ₁ vs ρ₂, factor ρ₂/ρ₁), ch04 vs ch07 η, three p′ bases — the code
  names the reference in the keyword or the function name.
- **Benchmarks from the primary source at its printed precision** (Fenton & McKee 1.7 %, Guo 0.75 % for β = 2.4908) —
  a later summary's numbers (1.5 %, 0.7 %) were not reproducible.
- Signed and ordered quantities stay signed/ordered: b·n, (u − b)·n, wall stress, lapse rate, ∇ρ × ∇p, u₂ − u₁, Γ_cw vs
  Γ_ccw, dw/dz = u − iv, **dω/dk, clockwise orbits, K up-left**. Probe piecewise functions at their joins; check a field
  against its own stream function; test every scenario's full balance.
- Library conventions: `ellipk(m)` with m = k²; `np.gradient` first order at edges; `np.sqrt`/`np.log` principal
  branches; pint °C offset; kmol vs mol; geopotential vs geometric altitude; FFT evolution is periodic (no wrap of the
  fastest group; ω(0)).
- **Validation labels**: V4 only for real conservation laws; identities V1; published closed forms reproduced identically
  are V1 form cross-checks, not V5; "converged" only with an asserted order (ch07 review: ka-scaling checks are analytic +
  V7, not "converged"); labelled extensions (simple wave, St Andrew's cross) say so in code, figures and explainers.

## Benchmark inventory
| value / table | source (citation) | stored in | used by |
|---|---|---|---|
| k_B, N_A, R (exact) | CODATA 2018/2022, NIST CUU | `reference/ch01/constants.json` | ch01 `thermo` |
| USSA-1976 constants, Tables 1–2; sea-level c = 340.29 m/s | NASA-TM-X-74335; PDAS tables | `reference/ch01/ussa1976_*` | ch01, ch04 C02 |
| IAPWS σ(T), μ(T) | IAPWS R1-76(2014), R12-08 | `reference/ch01/` | ch01, ch04, **ch07 C07 (σ(20 °C) = 72.74 mN/m)** |
| Jennings mean free path; Taylor blast K = 0.856; AMS Γ_d; EOS-80 | Tsalikis et al. 2024; Díaz arXiv:2009.05674; AMS Glossary; Fofonoff & Millard 1983 | `reference/ch01/benchmarks.json` | ch01 |
| Divergence theorem 8π/3; Levi-Civita; Mohr; rotations; Stokes | Wikipedia; scipy docs | `reference/ch02/` | ch02 |
| Lamb–Oseen vortex peak α = 1.256 | Canivete Cuissa & Steiner 2022, arXiv:2210.05223 | `reference/ch03/` | ch03 C14 |
| Sphere drag C_D(Re), 0.1–1e6 | F. A. Morrison (2016), Michigan Tech | `reference/ch04/` | ch04 C15 |
| WGS-84 a, 1/f, b, ω | Wikipedia "World Geodetic System" | `reference/ch04/` | ch04 C09, ch05 |
| Orifice C_c = 0.611; capillary length 2.71 mm; Prandtl ranges | Wikipedia | `reference/ch04/` | ch04 C05, C14, C15 |
| Sphere added mass = ½ displaced mass | Lamb, *Hydrodynamics* (1932) §92; Wikipedia "Added mass" | `reference/ch06/` | ch06 C15 |
| Rayleigh cavity collapse 0.91468 R₀√(ρ/Δp) | Rayleigh, Phil. Mag. 34 (1917) 94, doi:10.1080/14786440808635681 | `reference/ch06/` | ch06 exercise |
| Source panels (pedigree) | Hess & Smith (1967), doi:10.1016/0376-0421(67)90003-6 | `reference/ch06/SOURCES.md` | ch06 C14 |
| **Explicit inverse dispersion: Fenton & McKee max λ error 1.7 % (ours 1.658 %), ν_opt 1.49 (1.483); Guo β = 2.4908 → 0.75 % (0.753 %)** | Coastal Eng. 14, 499–513 (1990), doi:10.1016/0378-3839(90)90032-R; Coastal Eng. 45, 71–74 (2002), doi:10.1016/S0378-3839(02)00039-X (search-confirmed); Fenton (2006) note | `reference/ch07/` | ch07 C03 (V5) |
| **Capillary–gravity minimum 0.23 m/s at 1.7 cm** (ours 0.2312 at 1.712 cm) | Wikipedia "Capillary wave" + IAPWS σ | `reference/ch07/` | ch07 C07 (V5) |
| **Stokes drift ωka²e^{2kz}; Stokes speed (1 + ½(ka)²)√(g/k); limiting H/λ = 0.1410633, 120° crest** | Wikipedia "Stokes drift"; HandWiki "Physics:Stokes wave" (digits via Dyachenko, Lushnikov & Korotkevich 2016; Schwartz & Fenton 1982); arXiv:1507.02784 | `reference/ch07/` | ch07 C12 (V5) |
| St Andrew's cross at arccos(ω/N) (pedigree) | Mowbray & Rarity, J. Fluid Mech. 28, 1–16 (1967) | `reference/ch07/SOURCES.md` | ch07 C16 |
| Forms (V1 cross-checks, not V5): Rankine, Lamb–Oseen, RTT (ch03); Bélanger, Tsiolkovsky, Taylor–Green (ch04); Kelvin ring, Burgers, Hill (ch05); cylinder, Kutta–Joukowski, Blasius, Joukowsky, McDonald sphere (ch06); **Bélanger + head loss, KdV/cnoidal/solitary/Ursell, two-fluid capillary–gravity dispersion (ch07)** | Wikipedia pages; K. T. McDonald | `reference/ch03/`–`reference/ch07/` (`make_refs.py`, `SOURCES.md`) | ch03–ch07 |
Book-printed values (private, git-ignored): `tests/book_values_ch01.json` … `…_ch07.json` (ch07: 49 printed closed forms
≤ 1e-12 or 1e-6; thresholds and ocean numbers within the book's rounding; Exercise 7.6 lake period +0.72 %).
**Open**: a measured high-Re cylinder C_p (ch06 N31); a primary number for the minimum group velocity (ch07 O7).

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
| ch07 | 38 (V1 60) | 28 (V2 36) | 17 (V3 10) | 8 (V4 5) | 5 (V5 5) | 3 (V6 3) | 0 (2 labelled extensions: simple wave, St Andrew's cross) | 0 | **PASS** — 132 tests (+ V7 13); 41/41 wrong variants; review 2 Must (signed c_g for k < 0, Stokes γ attribution) + 8 Should resolved; every ★★/★★★ D re-derived in sympy; 9 explainers PASS round 2 (222 rows, 3 312 views, 29 D matched); notebook PASS round 2 (617 cells, ~189 s) |
(V7 tests: ch01 19, ch02 12, ch03 7, ch04 3, ch06 2, ch07 13. Full suite at the ch07 merge gate: **882 tests**.)

## Open across chapters
- **Library pass (orchestrator, when nothing reads `viz/`)**: the ranked ch07 list in `knowledge/viz_patterns.md`
  (pager-label nowrap CSS, `Viz.arrowPx`/`dotPx`, y-title gutter + square plots, formatters incl. length/speed,
  **`Viz.waves` before Ch. 13**, per-mode hiding Q4, log axes, `Viz.cx`, Gauss–Legendre), the four ch03 library bugs and
  quirks Q1–Q8; then `tools/viz_inline.py --all` and `tools/shot.py --chapter ch01…ch07 --quick` +
  `templates/viz_example.html --quick`.
- **Skill pass**: the lesson candidates for `interactive-viz`, `math-to-python` §7, `verify-implementation`,
  `teaching-style` and `colab-notebook` in `knowledge/viz_patterns.md` (ch02–ch07) are not yet in the skills (ch07 first:
  never type a number a function can print + 4 s.f. rule, complex-step continuations, sign-safe gradient forms, primary
  benchmark bounds, measurable linearisations, square K plane on phones).
- **Machinery TODO**: `viz_lint` single-backslash TeX rule; `run_notebook` refuses `MPLBACKEND=Agg`; nbkit converter for
  design TeX + ch07's `tidy_raw_tex` guard into nbkit; `coverage_check` literal TeX; `shot.py` clicks/drags and exposes
  `np` to `py:` rows; pager re-pack after KaTeX; rule-9 scan of explainer strings and notebook markdown against the
  private JSON; `eq_refs` exempts `<meta>` and Derivation/Code `ref:` badges; re-point ch04 `linear_wave_surface` to
  `core.waves.cosh_over_sinh`; keyword-only sibling signatures in `core.waves` (review Should-fix 6, documented only).
- Non-blocking explainer follow-ups: ch01–ch06 (their `knowledge/chNN.md` §9); **ch07** E1 (0.970, phone c(λ)), E2 (phone
  pages), E3 (4-s.f. live label near λ_m), E4 (phone intro), E5 (`lam0`, axis unit), E6 (lee criterion, phone labels),
  E7 (intermittent 3 px clip), E8 (whole (7.110), (7.96) note), E9 (ray label, k-axis title) — `knowledge/ch07.md` §9.
  Backups not built: ch04 `kinematic_free_surface`, ch05 `vortex_rings`, ch06 `flow_net_sources_vortices`, ch07
  `linearised_free_surface`.
- ch07 open verification notes: O5 Exercise 7.6 +0.72 % (book value), O6 Guo search-confirmed, O7 no primary c_g,min; ch06
  O4, O6; ch04 O5 Prandtl numbers; ch05 O2 tube-flux route 1.8e-4.
