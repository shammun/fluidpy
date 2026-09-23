# CUMULATIVE knowledge — fluidpy
(Rewritten by the knowledge-keeper after every chapter. Every agent reads this first. Last rewrite: after ch05,
2026-09-23, commit 33bb48a.)

## Project rules in force (added during ch01–ch03; unchanged in ch04–ch05)
- **5–10 explainers per chapter, as many as the CORE ideas need** (`book.yaml → project.min/max_explainers_per_chapter`).
  ch03 has 7, ch04 9, ch05 9.
- **Every book equation is shown in full next to its number** — notebook prose, derivation steps, traps, recaps and
  explainer tour/Explain/Derivation/quiz/notes/status text. Enforced by `tools/coverage_check.py` check 8 (markdown) and
  `tools/eq_refs.py` (explainers; still flags selftest names and `<meta>` content, quirk Q8). Builders use a
  `show_eqs`/`finalize_equations` helper (ch03–ch05 builders). The lesson reviewer also scans by eye: a cell with *other*
  maths passes check 8 while citing a bare number (ch05 round 2).
- **`coverage_check` check 9**: no control characters in cells (a LaTeX backslash eaten by a non-raw string).
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
  (D04) [C03]; momentum (4.17) (D05) [C04], Bernoulli (4.19) (D06) [C05] → Cauchy (4.24) (D07) [C06]   6 < 13
  constitutive τ = −pδ + 2μS + λS_mmδ (4.31) (D08, D09 ★★★, D10) [C07] → NS (4.38)–(4.39b), −μ∇×ω (D11–D13) [C08]  4 < 5
  rotating frame +2Ω×u′ (4.43) / −2Ω×u′ (4.45) (D14–D18, D15 ★★★) [C09]; energy, ε ≥ 0 (D19–D23) [C10]   7 = 7
  Lamb identity → ∇B = u × ω (4.69) (D24, D25) [C11]; unsteady Bernoulli (4.75) (D26) [C12]; Boussinesq (4.86), (4.89)
  (D27, D28) [C13]; Dη/Dt = 0 (4.91), jumps, Laplace (D29) [C14]; dimensionless NS (4.101) (D30) [C15]

ch05 Vorticity dynamics (done) — vorticity's life
  where it lives: vortex lines (5.3), tubes; Gauss on a tube piece −Γ_lower + Γ_upper = 0 (5.4) (D01): cannot end   [C01]
  what it does to pressure: tank p − p_o = ⅛ρω²r² − ρgz (5.6) (D02; ω = vorticity, tank at ω/2); funnel (5.7);
     line vortex: stress −μΓ/πr² but no net viscous force (D03)                                      [C02]
  when it is conserved: material loop (5.9) (D04) + momentum (5.10) → (5.11) → Kelvin DΓ/Dt = 0 (5.8) (D05; needs
     ρ = ρ(p)) [C03]; pressure torque on a disc = ∇ρ×∇p/ρ² (D06), lock exchange (D07) [C04]; Helmholtz (D08) [C05]
  how it changes: curl of (4.39b) (5.12) → Dω/Dt = (ω·∇)u + ν∇²ω (5.13) (D09 ★★★) [C06];
     rotating Lamb form (5.25) (D14) → curl → Dω/Dt = (ω+2Ω)·∇u + ∇ρ×∇p/ρ² + ν∇²ω (5.30) (D15 ★★★) [C09];
     (ω·∇)u = ω∂u/∂s (5.31) (D16) → stretching/tilting (5.32) (D17), stretched tube + Burgers (D18) [C10];
     absolute circulation DΓ_a/Dt = 0 (5.33) (D19) → fluid column (ζ + f)/h = const (D20) [C11]   ← seed of PV (Ch. 13)
  how it makes velocity: ∇×ω = −∇²u → Green's function, (5.14) with +1/(4π) (D10 ★★★) → Biot–Savart (5.16)
     (D11 ★★★) [C07] → filament law (5.17) (D12), segment and infinite line Γ/2πd (D13) [C08]
  how vortices move each other: point vortices, centre of vorticity (D21) [C12]; wall image Γ/4πh (D22), circle and
     channel images, rings [C13]; vortex sheet γ = u₂ − u₁ (D23), roll-up [C14]

next: ch06 ideal flow — Kelvin (ch05 C03, the four restrictions) says irrotational flow stays irrotational ⇒ u = ∇φ,
      ∇²φ = 0; ψ tools (ch04 C03), images (ch05 C13: wall, circle `circle_image_system`), cylinder with circulation
      (point vortex + ch03 cylinder), Bernoulli for pressure (ch04 C05/C11), Biot–Savart superposition (ch05 C07–C08),
      Hill's exterior = sphere flow; stress without net force (ch05 D03) for "viscous potential flow".
```
The book ahead: Ch6 ideal flow · Ch7 gravity waves · Ch8 laminar flow → Ch9 boundary layers → Ch10 CFD → Ch11
instability → Ch12 turbulence → **Ch13 GFD (uses 4, 5, 7, 8, 11, 12)** → Ch14 aerodynamics · Ch15 compressible · Ch16.
Where ch01 feeds in: N², θ → Ch. 7 §7.8, Ch. 11 §11.7, Ch. 13; Newtonian τ, ν → Ch. 8, 9; hydrostatic base → Ch. 7, 13;
isentropic gas → Ch. 15; FTCS → Ch. 8, 10. ch02: Gauss → every CV; Stokes/circulation → Ch. 6 lift, Ch. 13 PV;
invariants → Ch. 12. ch03: vortices → Ch. 8, 13, 14; streamlines/ψ/cylinder → Ch. 6; path lines → Ch. 7, 10; principal
strain → Ch. 12, 13 frontogenesis. ch04: **NS residual/term tools → every exact solution of Ch. 6–13**; ψ tools → Ch.
6, 8, 9, 13; unsteady Bernoulli + kinematic BC → Ch. 7; CV budgets → Ch. 9, 14, 15; ε → Ch. 12; rotating frame,
Boussinesq, Ri, Fr′, Ro → **Ch. 13**. **ch05**: Kelvin + images → Ch. 6 (lift, starting vortex, cylinders); vorticity
budget → Ch. 8 (decay), Ch. 10 (ω–ψ, vortex methods), Ch. 12 (stretching, Burgers); sheet roll-up → Ch. 11 (KH);
**(5.30), (5.33), (ζ + f)/h, baroclinic torque → Ch. 13 (PV, QG, Rossby waves, fronts, thermal wind)**; filaments,
rings, images → Ch. 14 (lifting line, induced drag, ground effect).

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
| constants `K_B, N_A, R_U, R_AIR, GAMMA_AIR, CP_AIR, G0, P_ATM, P_REF`; **`G_BOOK` = 9.81 (ch05)**; perfect gas, first law, entropy, Gibbs, sound speed, isentropic ratios, van der Waals, `helmholtz_free_energy` | `thermo` | §1.8–1.9, (4.94) | analytic, symbolic, converged, conserved, benchmark |
| `integrate_hydrostatic`, `standard_atmosphere`, `scale_height`, `buoyancy_force`, `net_pressure_force_on_box` | `statics` | (1.8), (1.9) | analytic, converged, benchmark |
| N² family, `classify_stability`, `parcel_displacement`, `adiabatic_lapse_rate` (negative), `lapse_rate_convention`, `lapse_rate_stability`, θ, ρ_θ | `stratification` | (1.29)–(1.35) | analytic, symbolic, benchmark |
| `ftcs_diffusion_1d`, `stable_time_step`, `couette_startup_profile`, `gaussian_spreading` | `diffusion` | §1.5 | converged, conserved, analytic |
| `dimensional_matrix`, `rank_by_minors`, `pi_groups`, `rescale_units` | `dimensional` | (1.36)–(1.40) | symbolic, analytic |
| index notation; `tensors` (passive C, `transform_*`, `traction` first index, `double_dot(convention=)`, `invariants`, `levi_civita`, `rotation_tensor` = G − Gᵀ, `principal_axes`); grids `[k, j, i]` + `operators` (`vector_gradient` G[i, j] = ∂u_i/∂x_j, `tensor_divergence(index=…)`, `curl`, `laplacian`); `fields`; `integral_theorems` (Gauss, Stokes, `circulation`, `gauss_legendre_nodes`, **`curl_theorem_box` (ch05)**) | ch02 modules | (2.1)–(2.36) | analytic, symbolic, converged, conserved, benchmark |
| `kinematics` (field callables `u(x, t)`, material derivative, flow lines, `vorticity`, strain/spin, `vorticity_in_rotating_frame`, linear-flow presets); `coords`; `transport` (Leibniz, CV shapes, `reynolds_transport`) | ch03 modules | (3.1)–(3.35) | analytic, symbolic, converged, conserved |
| `vortices`: ch03 `solid_body_rotation(r, omega0 = rate)`, `line_vortex`, `rankine_vortex`, `gaussian_vortex`, `circulation_circle`; **ch05 + `burgers_vortex(_field)`, `burgers_pressure`, `burgers_core_radius`, `hill_spherical_vortex(_field)`, `hill_stream_function`, `lamb_oseen_field`, `gaussian_tube_*`, `broken_tube_field`, `vortex_ring_vorticity`, `stretched_gaussian_vortex_field`, `abc_flow_field`** | `vortices` | (3.22)–(3.29), §5.4 | analytic, symbolic, converged, benchmark |

### Conservation-law primitives from ch04 (`ch04.<name>`; details `knowledge/ch04.md` §3)
| function family | module | book § / Eq. | label |
|---|---|---|---|
| `mass_budget`, `momentum_budget`, `energy_budget`, `angular_momentum_budget` on any ch03 `ControlVolume`; flux (u − b)·n | `conservation` | (4.5), (4.17), (4.48), (4.65) | analytic, conserved, converged, symbolic |
| `newtonian_stress`, `viscous_stress`, `mean_pressure(tau, tau33)`, `bulk_viscosity`, `stress_on_plane`, `deviatoric_part`, `dissipation_rate` | `constitutive` | (4.25)–(4.37), (4.58) | analytic, symbolic |
| residuals and term splits (numeric + sympy): continuity, Cauchy (first index), NS, `viscous_force_forms`, `exact_solution` (8), `lamb_vector`, energy/entropy, Boussinesq families | `navier_stokes` | (4.7)–(4.89) | analytic, symbolic, converged, conserved |
| Appendix-B operators (cylindrical, spherical), sympy | `curvilinear` | App. B | analytic, symbolic |
| ψ ↔ u (planar u = ∂ψ/∂y, axisymmetric), fluxes, presets | `streamfunction` | (4.11), (4.12) | analytic, symbolic, converged |
| `frame_acceleration_terms` (+2Ω×u′), `apparent_body_forces` (−2Ω×u′), `coriolis_parameter`, `effective_gravity`, `OMEGA_EARTH` | `rotating` | (4.42)–(4.45) | analytic, symbolic, converged, benchmark |
| `bernoulli_function`, `pressure_function`, `which_bernoulli(_text)`, unsteady family, `rankine_vortex_pressure`, stagnation, pitot | `bernoulli` | (4.19), (4.66)–(4.83) | analytic, symbolic, conserved, benchmark |
| `kinematic_bc_residual`, `interface_mass_flux`, `laplace_jump_from_balance`, `capillary_length` | `interfaces` | (4.90)–(4.98) | analytic, converged, benchmark |
| `Scales`, 20+ named numbers (Re … Ro), C_p/C_D/C_L, `sphere_drag_coefficient` (Morrison) | `similarity` | (4.99)–(4.119) | analytic, symbolic, benchmark |

### Vorticity primitives from ch05 (`ch05.<name>`; details `knowledge/ch05.md` §3)
| function family | module | book § / Eq. | label |
|---|---|---|---|
| loops (`circle_loop_points`, `loop_circulation` spectral, `loop_vector_area`), `material_loop`, `material_circulation`, `kelvin_rate_terms`, `kelvin_force_terms`, `absolute_circulation` → (Γ, Γ_a), Kelvin scenarios | `vorticity` | (5.8)–(5.11), (5.33) | conserved, analytic, converged, symbolic |
| `vorticity_field`, `vortex_line`, `vortex_tube_strength`, `vortex_tube`, `tube_flux_budget(_traced)`, `frozen_in_check` | `vorticity` | (5.3)–(5.4), §5.3 | analytic, symbolic, conserved, converged |
| **`vorticity_budget(u, x, t, nu, Omega, rho, p)`** → local, advective, relative, planetary, baroclinic, diffusion, residual; `_sym`; `vorticity_terms`, `vorticity_equation_sym`, presets; `baroclinic_term` (∇ρ×∇p/ρ²), `planetary_vorticity_terms`, `stretching_tilting_split`, `rotating_lamb_form_terms`, `pressure_torque_on_element` | `vorticity` | (5.12)–(5.32) | symbolic, analytic, converged |
| `poisson_green_3d`, `velocity_from_curl_omega(sign=+1)`, `biot_savart_volume`, `biot_savart_2d`, `velocity_from_vorticity_fft` (periodic) | `biot_savart` | (5.14)–(5.16) | analytic, converged, symbolic |
| `segment_induced_velocity`, `filament_velocity`, `filament_preset`, `ring_axis_velocity`, `ring_ring_velocity` (m = k²), `ring_self_velocity` (`RING_CORES`) | `biot_savart` | (5.17), §5.7 | analytic, converged |
| `point_vortex_velocity`, `point_vortex_rhs/evolve(boundary=wall/circle/channel)`, `point_vortex_invariants`, `centre_of_vorticity`, image systems, `vortex_sheet_velocity`, `continuous_sheet_velocity` | `biot_savart` | §5.7–5.8 | analytic, conserved, converged |
Chapter-only code stays in chapter modules (move to `core/` when a second chapter calls it): ch03 `cylinder_flow` (→ Ch.
6), `pipe_profile` (→ Ch. 8); ch04 `linear_wave_surface`, `bore_speed` (→ Ch. 7), `couette_heating`, `two_fluid_couette`,
`stokes_first_problem` (→ Ch. 8), `wake_drag_per_span` (→ Ch. 9), `accelerating_sphere_*` (→ Ch. 6); **ch05**
`rotating_cylinder_flow`, `torque_per_length`, `lamb_oseen_circulation` (→ Ch. 8), `diffusing_vortex_sheet`,
`sheet_rollup` (→ Ch. 11), `column_relative_vorticity`, `column_over_slope`, `relative_circulation_after_move` (→ Ch.
13), `kelvin_hypotheses(_text)` (→ Ch. 6), `ring_dynamics` (→ Ch. 14), `uniform_strain_vorticity` (→ Ch. 12). Note `ch04.G`
and `ch05.G` = `G_BOOK` = 9.81 while `core` defaults to `G0` = 9.80665.

### Explainer engine (`assets/viz_lib.js`)
`Viz.app` (tabs Walkthrough / Explore / Explain / Derivation / Equations / Code / Check; fit-to-window with density
levels and pagers; transport, modes, presets, status, terms, inspector, notes, selftest parity), `Plot`, `Viz.field`,
`Viz.num` (RK4, odeint, brentq, erf/erfc **~1.2e-7 only**, niceTicks), `Viz.work`, `Viz.three`, KaTeX with fallback;
`Viz.font`, `Viz.roundRect`, `Viz.text`, `Viz.card`, `Viz.fmtTime`, `Viz.tnum(keepTiny)`; pager slicing,
`--viz-stage-need`, `ev.viewId`/`ev.vizView`. **Nothing promoted after ch02, ch03, ch04 or ch05** (the site-publisher
was reading `viz/` each time). Ranked candidates (ch05 list: formatters in all 35 explainers, per-mode control hiding
16, `arrowPx` 7 (+2 variants), `vec3` 5, log axes 5, bars 11, field image 7, Gauss–Legendre 5, `rk4` for arrays 3, 3-D
canvas projector 3), four open library bugs (ch03) and eight quirks Q1–Q8 (ch04) are in `knowledge/viz_patterns.md`.
Reference: `templates/viz_example.html`.

## Explainer inventory
| chapter | slug | CORE idea | reusable stage pattern |
|---|---|---|---|
| ch01 | `continuum_averaging_volume` | C06 continuum (+Kn, D35) | molecules · value vs log box size with band · regime strip; click-a-sample inspector |
| ch01 | `viscosity_momentum_diffusion` | C12 Newton viscosity, ν | gap with tracers · profile with steady ghost and trail · wall stress; modes = different diffusivity |
| ch01 | `heat_work_paths` | C25, C35, C45 | piston · p–v · T–s linked on a path parameter; hatched areas; path vs state term bars |
| ch01 | `parcel_stability` | C50, C51, C54, C55 | column · profile with adiabat · ζ(t) with linear ghost; two-convention badge with exact-text parity |
| ch01 | `buckingham_pi_machine` | C64, C67, C69 | variable chips → matrix with minor → groups; exact rational JS; unit-system switch |
| ch02 | `rotation_of_axes` | C02, C03, C06 | plane with fixed arrow and turning axes · clickable C matrix · τ′(θ); passive/active table |
| ch02 | `cauchy_traction_principal_axes` | C04, C05, C13 (D17 ★★★) | element with a cut, f split blue/rose · σn/τs vs φ · Mohr; term bars |
| ch02 | `strain_vs_rotation_split` | C12 | one clock, three squares (G, S, A); `expm2`; Frobenius bars |
| ch02 | `gauss_flux_box` | C14, C15 | heatmap + draggable box, signed face arrows · waterfall flux bars · limit view vs log h |
| ch02 | `stokes_circulation_loop` | C11, C16 | curl image + paddle wheel + draggable loop · unrolled u·t · Stokes failing on purpose |
| ch03 | `flow_lines_unsteady` | C03, C04 | unsteady flow with streamline, particle trail and dye on one clock · answer-sheet view |
| ch03 | `material_derivative_probe` | C01, C02 | fixed probe + carried float · local/advective/total bars with measured ◇ |
| ch03 | `galilean_frames_cylinder` | C05 | observer slider lake → body · bars that trade places · u = U + u′ triangle |
| ch03 | `fluid_element_deformation` | C06–C09 | four G sliders deform a square + ring · rate curves with measured dots |
| ch03 | `spin_and_principal_axes` | C10–C12 | threads + paddle wheel · single-thread rates vs flat pair average · co-rotating observer |
| ch03 | `vortex_paddle_wheels` | C13, C14 | vortex with orbiting wheels and a draggable loop · u_θ, ω_z, Γ(r) · real-vortex table |
| ch03 | `reynolds_transport_cv` | C15 (D22 ★★★) | moving CV with signed swept band · budget waterfall + measured ◇ · dropped-term log–log |
| ch04 | `control_volume_budgets` | C01, C04 | five scenes on one budget object · face flux bars + measured ◇ · "drag b until storage vanishes" |
| ch04 | `stream_function_spacing` | C03 | ψ at equal Δψ over a speed map · spacing = speed bars · bendable gate · axisymmetric mode |
| ch04 | `newtonian_stress_lab` | C07 (D09 ★★★) | G → S, R → τ grid coloured by source → rotatable plane traction · cube spin-up mode |
| ch04 | `navier_stokes_term_balance` | C08 | 7 exact solutions, click a point, term bars summing to 0, regime word |
| ch04 | `rotating_frame_coriolis` | C09 (D15 ★★★) | two observers on one clock · signed side toggle · split Coriolis arrow · f-plane highs/lows |
| ch04 | `which_bernoulli` | C05, C11, C12 | six flows with probes along/across streamlines · hypothesis decision table · stacked B terms |
| ch04 | `viscous_dissipation_heating` | C10 | Couette gap: profile, ε, T(y, t) · in = stored + out bars · "μ < 0 ⛔" preset |
| ch04 | `boussinesq_buoyancy` | C13 | rising blob with buoyancy-off ghost · kept vs dropped bars · validity chart with lit settings table |
| ch04 | `dynamic_similarity_models` | C15 | prototype and model on one t* clock · paired group bars · sphere collapse · Ro preview |
| ch05 | `vortex_tubes_cannot_end` | C01 (D01) | 3-D tube with a sliding section · flux/area/mean ω along the tube · Gauss bars; "broken ⚠" field fails by ∫∇·ω dV |
| ch05 | `vortex_pressure_funnel` | C02 (D02, D03) | four vortices on one (r, z) section · u_θ, B, σ_rθ profiles · needed = supplied force bars with FD slope; real-vortex table |
| ch05 | `kelvin_material_loop` | C03, C05 (D04, D05, D08) | six flows + material vs fixed loop on one clock · ✓/✗ hypothesis table with surviving term · measured ◇ dΓ/dt |
| ch05 | `baroclinic_torque` | C04 (D06, D07) | disc with tilted isopycnals, rim forces, G · torque ◇ on the formula's sine · R² gap panel; lock-exchange mode |
| ch05 | `vorticity_stretching_tilting` | C06, C10 (D09, D16–D18) | 3-D vortex line with purple stretching / blue tilting arrows · \|ω\|(t) · Burgers balance settling |
| ch05 | `biot_savart_filament` | C07, C08 (D10–D13, two ★★★) | filament + field point · "build the sum" transport · unrolled integrand · speed vs distance; printed-sign toggle |
| ch05 | `vorticity_equation_rotating` | C09, C11 (D14, D15 ★★★, D19, D20) | column over a ridge / ring moved poleward with the conserved quantity flat in amber · 7-term budget · f-table |
| ch05 | `point_vortex_lab` | C12, C13 (D21, D22) | click-to-place vortices (guarded) · trajectories · invariants Q/Q(0); wall and bucket images |
| ch05 | `vortex_sheet_rollup` | C14 (D23) | row of filaments vs continuous jump · draggable circuit with trading side bars · L1 vs N · roll-up with KH ghost |

## Notation
See `knowledge/notation.md` (symbol register and sign traps). Conventions that matter everywhere: SI and kelvin inside
functions; **z up**; **lapse rate Kundu Γ ≡ dT/dz in code, meteorology −dT/dz shown alongside**; `p0` ≠ `p_ref`; q to /
w on the system; gas constants per kmol; pressures absolute unless `_gauge`; signed τ_xy = μ ∂u/∂y.
From ch02: passive C (x′ = Cᵀx); traction on the first index; book A:B = A_ij B_ji; G[i, j] = ∂u_i/∂x_j; Stokes n_c
into A; grid `[k, j, i]` with `h` in (x, y, z).
From ch03: R = G − Gᵀ (no ½), ω = ∇×u, spin ½ω; Galilean u′ = u − U; rotating ω′ = ω − 2Ω; RTT outward n, signed b·n;
field callables `u(x, t)`; spherical θ from +z.
From ch04: Cauchy's divergence on the first index; flux (u − b)·n; **Coriolis acceleration +2Ω × u′ (4.43) vs force
−2Ω × u′ (4.45)**; Stokes assumption μ_v = 0; 2-D ψ: u = ∂ψ/∂y (GFD often opposite); the prime's four meanings with
one code name each; g defaults 9.80665 (core) vs 9.81 (`G_BOOK`, chapter modules).
From ch05:
- **ω is the vorticity**; a tank turning at Ω has ω = 2Ω (ch03's `omega0` is the rate; Fig. 5.2's "2ω" is a slip).
- **Γ = circulation [m²/s] always; γ = sheet strength [m/s] = u_below − u_above**; counterclockwise positive for ω_z, Γ,
  point vortices and sheets; Γ_a = Γ + 2Ω·A_vec. Γ/γ now has five-plus meanings project-wide (notation table).
- **Baroclinic source ∇ρ × ∇p/ρ² in that order**; E4 angle ∇ρ = |∇ρ|(sin θ, −cos θ) with ∇p = (0, −ρ₀g): 0 < θ < 180°
  spins clockwise; lock exchange with heavy fluid on the left spins counterclockwise.
- **(5.14) with +1/(4π)** (book prints −); ∇′(1/|x − x′|) = +(x − x′)/|x − x′|³.
- **e_n = −(Frenet N)** in the natural frame (5.31); τ = torsion there.
- **scipy `ellipk`/`ellipe` take m = k²**; ring cores ¼ / ½ / 0.558; Burgers' book α = 2 × Wikipedia's.
Symbols with several meanings so far: α, θ, T, h, R, σ/S, c, D, n, k, q, γ, λ, τ, ε, H (ch01); Γ, ω, φ, A, C, G, K, m, s,
t, a, b (ch02); Γ, ω, σ, φ, θ, α, a, b, h, primes, U (ch03); primes (4), σ, ε, γ, Φ, Ψ, Ω, h, H, M, R, f, F, B, C_p, λ,
b, η, ζ, G (ch04); **Γ/γ, ω vs Ω, σ (viscous stress vs core radius), ε (Levi-Civita, dissipation, kernel smoothing),
ζ, h (column height, spacing, wall distance), α (strain rate), s (arc length, shear rate), τ (torsion), κ (curvature), m
(elliptic parameter), H (channel width, Kirchhoff energy), N (filaments), G (centre of vorticity), A (Hill constant) (ch05)**.

## Teaching lessons
- **Depth is tiered; coverage is exhaustive.** ch01: 15 A / 70 B / 21 C, 12 derivations, 496 cells. ch02: 16 A / 60 B /
  18 C, 15 derivations, 405 cells. ch03: 15 A / 52 B / 12 C, 24 derivations (171 steps), 413 cells. ch04: 15 A / 151 B /
  19 C, 30 derivations (235 steps), 611 cells. ch05: 14 A / 69 B / 11 C, 23 derivations (194 steps), 496 cells.
- **Convention callouts** (lapse rate, "Which p_o?", passive/active, first/second index, A:B, R vs A, γ vs Γ vs
  circulation, translating vs rotating frame, signed b·n, Coriolis term vs force, four primes, ψ sign, λ vs μ_v, **ω vs
  rate, Γ vs γ, ∇ρ × ∇p order, e_n vs Frenet N, m vs k**): state the book's, the field's and the code's convention,
  give the size of the difference in numbers, say which one the code uses and why.
- **Book typos are taught, not hidden — by computing both versions**: "book prints X; Y is right" with a test or a
  curve that shows the printed form wrong (ch02 Ex. 2.3/2.4/2.6; ch03 (3.6); ch04 (4.15), (4.51), (4.74); **ch05 (5.14)
  sign as a dashed reversed curve, Fig. 5.11's G moving at −1.70 m/s, the (5.27) dropped term kept in sympy, Fig. 5.2
  "2ω", Fig. 5.16 caption, Kelvin needs ρ = ρ(p), exercise pointers 5.8 → 5.9 and 5.11 → 5.10**).
- **Running threads** make a chapter's logic visible: ch04's equations-vs-unknowns ledger; ch05's "life of vorticity"
  (where it lives → pressure → conservation → change → velocity → interaction). **Decision tables** sort results by
  hypotheses (ch04 four Bernoulli forms; ch05 Kelvin's 16 combinations naming the surviving term).
- **Every key number by two or three independent routes** (formula vs measured/finite difference; ch05: net viscous
  force three ways, torque route vs formula with an R² gap, loop vs area to 12 digits) — also as a measured ◇ in the
  explainer.
- **★★★ sympy checks re-run the derivation's own construction** and **keep what the book drops** (ch05 D15 carries
  u_{j,j}(ω_n + 2Ω_n) for a generic u); a derivation check must describe exactly what its cell runs (failed in ch03,
  ch04 and ch05 round 1).
- **Every prose number is the printed number** (ch05 round 1: four mismatches; numpy print options hid 2.9e-9 as 0 —
  use format strings); every "Dxx step N" pointer is checked by script after a split.
- **Physical analogies and "what would change if" predictions are checked like equations** (ch05: a rotating tank is
  effective gravity, not geostrophy; viscous cellular flow decays as e^{−2νt} whatever the loop shape).
- **Make an approximation measurable; fence qualitative models** (ch03 ellipse ≈ γ/4; ch04 Boussinesq dropped bar;
  ch05 frozen kernel 2a/distance, ring model stopped at a 3a gap with `valid`/`stop_time`).
- **Primers come before the first use** (ch05 round 1: meshgrid `indexing="ij"`, Gauss–Legendre, `np.deg2rad`, polar
  vector Laplacian, Taylor–Green, moment transfer were late or missing); `knowledge/primers.md` lists 148 (P01–P148);
  ch06 starts at P149.
- **Climate hooks with numbers** (ch04 pole projectile, geopotential, Ro preview; ch05 tornado deficit, lock
  exchange/sea breeze 2.42 s⁻², planetary stretching 1.46e-9 s⁻², ridge ζ = −0.2f, ring of air −4.19e7 m²/s, Fujiwhara).
- Reusable derivation moves (ch01–ch05, 84 so far; 21 from ch05) are tabulated in `knowledge/concept_map.md` →
  "Reusable derivation moves".

## Global pitfalls confirmed in this project
**Machinery and environment**
- Extracted text garbles maths (`¼` = `=`, `ð…Þ` = parentheses, missing minus, ω → `u`, γ → g, θ → q, σ/τ → s, ν → n,
  ε → 3, ρ → r, Ω → U, ψ → j, χ → c, φ → 4, Φ → F, η → h, ζ → z, ∂ → v, ∇ → V) — read equations from rendered pages.
- **Windows: shell heredocs mangle backslashes and quotes** — write code, design and explainer files with Write/Edit;
  Python builder strings with LaTeX must be raw (coverage check 9); scan explainers for single-backslash TeX.
- Anaconda's `python3` kernelspec can shadow the venv's — notebooks execute on `fluidpy-venv`.
- JS: `UIEvent.view` is read-only. `tools/shot.py` does not click/drag; a reviewer runs a click pass.
- `test_viz_library_inlined_and_template_lints` fails whenever `assets/viz_lib.js` changed or an explainer is mid-build;
  run `tools/viz_inline.py --all` before the merge gate.
- `tools/shot.py` `py:` parity expressions have **no builtins**; parity rows call the explainer's own functions; set the
  tolerance near the achieved agreement (a 1e-3 row hid 5e-11 agreement in ch05 E9).
- Parallel viz-builders use **private scratch subfolders**; **runtimes measured while builders run browser audits are
  5× slower** (ch05 notebook 67 s alone vs 345 s; full suite 329 s vs 1156 s).
- **matplotlib 3.11**: animations using `subplots_adjust`/`fig.text` disable `figure.constrained_layout.use` for the save.
- `tools/coverage_check.py` matches ledger reminder rows by exact primer name → false warnings; `nb.derivation`
  prefixes "Eq." to its label.
- Book wording can hide in docstrings (`check_public.py` does not catch paraphrase): reviewers scan for it.
- Library: pager packs before the final `fit()`; tall ∫ clip at slice breaks; `Viz.num.erf/erfc` only ~1.2e-7;
  transport range fixed per app; no per-mode controls (Q4, hit in all 9 ch05 explainers); `onChange` after `set` (Q5).
- **Phone legibility is the usual viz failure** (ch01 E4, ch05 E2/E3/E7 round-1 Musts): rotated y titles over tick minus
  signs, labels over labels, > 7 term rows in ~120 px, tick crowding in < 60 px views, 25-page derivation results.
- **Design documents contain errors; downstream agents compute, not copy** (ch05: baroclinic angle reversed, D12 3a →
  2a, D23 "cancel" → vanish as dn → 0). Fix the design itself; builder text patches (`pf_sub`) break when the design is
  corrected upstream.
- Docstring contract numbers and derivation IDs go stale after fixes — re-grep them.

**Mathematics → code**
- **A test comparing two of our own functions is not evidence for the book's convention.** Pin conventions to fields
  with known answers and **prove discrimination by patching in the wrong variant** (ch01; ch02 8/8; ch03 13/13 + 7;
  ch04 26/26; ch05 40/40 + 5 re-planted) — and **choose a field where the wrong variant differs** (a symmetric G hid a
  transposed planetary term in ch05).
- **Probe piecewise functions exactly at their joins** (ch05 F1: a central difference straddling the Rankine kink gave
  half the torque and a 1/h "force"); **check a field against its own potential or stream function** (ch05 F2: Hill's
  exterior u_R sign).
- **Test every scenario's full balance per direction, not just totals** (ch04 oblique jet, 2 × 2 mean pressure).
- **Reduced-dimension inputs hide physics** (plane stress needs τ₃₃; 2-vectors padded explicitly; scalars raise).
- **Series solutions**: no early stop at a zero coefficient, no fixed-node projection (ch04 F1).
- **Preset names and letters are physics** (γ = 2S₁₂ vs Γ ≡ S₁₂ vs Γ circulation vs γ sheet strength; λ vs μ_v; Coriolis
  + vs −; ω vorticity vs rate) — each pinned by a test.
- **Signed and ordered quantities stay signed/ordered**: b·n, (u − b)·n, Leibniz lower term, wall stress, lapse-rate
  conventions, Stokes orientation, ω₃ = −γ, the (4.74) gauge, **∇ρ × ∇p, the sheet's u₂ − u₁, the (5.14) factor**.
- **Book sign slips can cancel** ((5.14) and its integrand rewrite): code the corrected intermediate, keep the printed
  form as an option, test that it reverses the physics.
- **Library parameter conventions**: `scipy.special.ellipk(m)` with m = k²; `np.gradient` first order at edges; pint °C
  offset; kmol vs mol; geopotential vs geometric altitude.
- **Models have ranges; users have fingers**: thin-core rings stop at a 3a gap (`stop_time`, `valid`); a vortex dragged
  to a circle's centre must not NaN the system; explainers guard separations and sub-step near approaches.
- **Rotating ≠ Galilean**: ω′ = ω − 2Ω; four apparent accelerations with the factor 2; Γ_a = Γ + 2Ω·A_vec.
- Index conventions: traction and Cauchy on the first index, ch02 `tensor_divergence` default on the second; passive C
  vs active R; book vs Frobenius double dot; R = 2A; grid `[k, j, i]`; natural frame e_n = −N.
- Units of numerical defaults; two g defaults; vectorised inputs via `as_scalar_if_0d`; closed loops wrap; tolerances
  relative to the size of the cancelling terms (ch05 stencil residuals vs |u||ω|/σ).
- sympy: assumptions make different symbols (use `coordinates`); `.norm()` adds `Abs`; `Subs(Derivative)` does not
  simplify; cache slow curvilinear simplifications; `sp.LeviCivita` for index-notation curls.
- Cancellation: `-np.expm1(-x)`; stiff geometry in arc length; unbounded ODEs stop with an event.
- **Validation labels**: V4 only for real conservation laws; unit/frame invariance is V7; identities V1; **published
  closed forms reproduced identically are V1 form cross-checks, not V5** (ch03 Rankine/Lamb–Oseen; ch04 Bélanger,
  Tsiolkovsky; ch05 Kelvin ring speed, García–Haziot pairs, Burgers, Hill, segment law); "converged" only with an
  asserted order; property-table differences are reported, not asserted.

## Benchmark inventory
| value / table | source (citation) | stored in | used by |
|---|---|---|---|
| k_B, N_A, R (exact) | CODATA 2018/2022, NIST CUU | `reference/ch01/constants.json` | ch01 `thermo` |
| USSA-1976 constants, Tables 1–2; sea-level c = 340.29 m/s | NASA-TM-X-74335; PDAS tables | `reference/ch01/ussa1976_*` | ch01, ch04 C02 |
| IAPWS σ(T), μ(T) | IAPWS R1-76(2014), R12-08 | `reference/ch01/iapws_sigma.csv`, `benchmarks.json` | ch01, ch04 capillary length |
| Jennings mean free path; Taylor blast K = 0.856; AMS Γ_d; EOS-80 check values | Tsalikis et al. 2024; Díaz arXiv:2009.05674; AMS Glossary; Fofonoff & Millard 1983 | `reference/ch01/benchmarks.json` | ch01 |
| Divergence theorem 8π/3; Levi-Civita identities; Mohr; rotation matrices; Stokes | Wikipedia; scipy docs | `reference/ch02/` | ch02 |
| Gaussian (Lamb–Oseen) vortex peak α = 1.256 | Canivete Cuissa & Steiner 2022, arXiv:2210.05223 | `reference/ch03/` | ch03 C14 |
| Sphere drag C_D(Re), 0.1–1e6 | F. A. Morrison (2016), Michigan Tech, Eq. (1) | `reference/ch04/` | ch04 C15 |
| WGS-84 a, 1/f, b, ω (reused by ch05 for Ω_E) | Wikipedia "World Geodetic System" | `reference/ch04/` | ch04 C09, ch05 C09/C11 |
| Sharp-orifice C_c = 0.611; capillary length 2.71 mm; Prandtl ranges | Wikipedia | `reference/ch04/` | ch04 C05, C14, C15 |
| Forms (V1 cross-checks, not V5): Rankine, Lamb–Oseen, RTT, Leibniz (ch03); Bélanger, Tsiolkovsky, added mass, Taylor–Green (ch04); **Kelvin thin-ring speed, García & Haziot (2023, arXiv:2204.11327) pair rates, Burgers, Hill, Biot–Savart segment, Kelvin/Poincaré–Bjerknes (ch05)** | Wikipedia pages; Commun. Math. Phys. | `reference/ch03/`–`reference/ch05/` (`make_refs.py`, `SOURCES.md`) | ch03–ch05 |
Book-printed values (private, git-ignored): `tests/book_values_ch01.json` … `…_ch05.json` (ch05: §5.1/5.7/5.8 forms,
Ex. 5.1 −0.23 %, Ex. 5.2 −0.13 %, exercise closed forms ≤ 1e-7). **ch05 has no V5 evidence** — a cited lab/DNS number
(ring speed, lock-exchange front speed) would be welcome when Ch. 7/13 revisit these flows.

## Validation summary per chapter
Counts are concept-map rows carrying each label (a row can carry several); tests per evidence level in parentheses.
| chapter | analytic | symbolic | converged | conserved | benchmark | book-value | qualitative | unverified | verdict |
|---|---|---|---|---|---|---|---|---|---|
| ch01 | 39 (V1 58) | 16 (V2 19) | 5 (V3 8) | 5 (V4 7) | 18 (V5 13) | 1 (V6 5) | 1 (`wavelength_to_rgb`) | 0 | **PASS** — 107 tests; 5 explainers (115 rows); notebook 496 cells |
| ch02 | 43 (V1 41) | 21 (V2 13) | 11 (V3 15) | 2 (V4 1) | 9 (V5 6) | 8 (V6 4) | 0 | 0 | **PASS** — 92 tests (+ V7 12); 8/8 wrong variants; 5 explainers (127 rows); 405 cells |
| ch03 | 38 (V1 58) | 23 (V2 32) | 14 (V3 14) | 5 (V4 4) | 2 (V5 1) | 8 (V6 4) | 0 | 0 | **PASS** — 120 tests (+ V7 7); 13/13 (+7) wrong variants; 7 explainers (145 rows); 413 cells |
| ch04 | 41 (V1 83) | 35 (V2 41) | 17 (V3 9) | 13 (V4 6) | 11 (V5 6) | 5 (V6 3) | 0 | 0 | **PASS** — 151 tests (+ V7 3); 26/26 wrong variants; 9 explainers (255 rows); 611 cells |
| ch05 | 39 (V1 90) | 24 (V2 23) | 13 (V3 10) | 13 (V4 8) | 0 (V5 0) | 5 (V6 4) | 2 (ring beyond a 3a gap, roll-up shape; both fenced) | 0 | **PASS** — 135 tests; 40/40 wrong variants + 5 re-planted fixes; review no Must-fix; every ★★/★★★ derivation re-derived in sympy; orders 1.99–2.01 (sheet 0.999 by design); 9 explainers PASS round 2 (268 rows, 2 776 views, 23 derivations matched); notebook PASS round 2 (496 cells, 67 s) |
(ch01 V7: 19 tests; ch02: 12; ch03: 7; ch04: 3; ch05: V7 checks inside V1 tests. Full suite at the ch05 merge gate:
**613 tests**.)

## Open across chapters
- **Library pass (orchestrator, when nothing reads `viz/`)**: fix the four ch03 library bugs and quirks Q1–Q8 (Q4
  per-mode controls and Q5 `onChange` after `set` first — every ch05 explainer works around Q4); make the narrow-view
  y-title strip an engine default; promote the ranked JS candidates (`knowledge/viz_patterns.md`: formatters (35/35),
  per-mode hiding, `arrowPx`, `vec3` before Ch. 13, log axes, `Viz.bars` with zero-row collapse, field image,
  Gauss–Legendre, double-precision erf); then `tools/viz_inline.py --all` and `tools/shot.py --chapter ch01…ch05 --quick`
  + `templates/viz_example.html --quick`.
- **Skill pass**: the lesson candidates for `interactive-viz` (ch05 phone-legibility checklist first), `math-to-python`
  §7, `verify-implementation` (ch04 "full balance per direction", ch05 "probe piecewise joins"), `teaching-style` and
  `colab-notebook` listed in `knowledge/viz_patterns.md` (ch02–ch05) are not yet in the skills.
- Machinery TODO: `shot.py` clicks/drags every view; `coverage_check` understands ledger reminder rows and flags a bare
  number even in a cell with other maths; `nb.derivation` skips "Eq." for non-numeric labels; `show_eqs` → `tools/nbkit.py`;
  `core.anim` saves inside a constrained-layout-off `rc_context`; `eq_refs` exempts selftest names and `<meta>`.
- Non-blocking explainer follow-ups: ch01–ch04 (their `knowledge/chNN.md` §9); ch05 E1 step 5 and "ω̄", E2 isobar labels,
  E3 round-off title, E4 Explore 10 pages, E5 clipped label, E6 "sum" vs "so far", E7 4-page steps, E8 invariants
  legend, E9 parity tolerances (`knowledge/ch05.md` §9). Backups not built: ch04 `kinematic_free_surface`, ch05
  `vortex_rings`.
- ch04 O5: the book's Prandtl numbers vs our property tables (−1.5 %, −2.2 %), reported not asserted. ch05 O2: tube-flux
  route at the default disc 1.8e-4 off (quote four digits).
