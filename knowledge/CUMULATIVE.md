# CUMULATIVE knowledge — fluidpy
(Rewritten by the knowledge-keeper after every chapter. Every agent reads this first. Last rewrite: after ch04,
2026-09-23, commit b4619bb.)

## Project rules in force (added during ch01–ch03; unchanged in ch04)
- **5–10 explainers per chapter, as many as the CORE ideas need** (`book.yaml → project.min/max_explainers_per_chapter`).
  ch03 has 7, ch04 has 9.
- **Every book equation is shown in full next to its number** — notebook prose, derivation steps, traps, recaps and
  explainer tour/Explain/Derivation/quiz/notes/status text. Enforced by `tools/coverage_check.py` check 8 (markdown) and
  `tools/eq_refs.py` (explainers; `ref:` labels, metadata and selftest names are meant to be exempt — ch04 found it
  still flags selftest names and `<meta>` content, quirk Q8). Builders use a `show_eqs(text)` helper (ch03/ch04 builders).
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
  material volume + coincident CV + RTT (3.35):
     mass (4.5) (D01) [C01] ─Gauss+localisation─► continuity (4.7) (D02) → (4.8)–(4.10) ∇·u = 0 (D03) [C02]
        └─ steady continuity solved by ψ: u = ∂ψ/∂y, v = −∂ψ/∂x, flux = Δψ (D04) [C03]
     momentum (4.17) (D05) [C04]: wake drag, bore → √(gh), rocket, jet, sprinkler; stream tube → Bernoulli (4.19) (D06) [C05]
        └─ Gauss (first index!) + localisation ─► Cauchy ρDu/Dt = ρg + ∂τ_ij/∂x_i (4.24) (D07) [C06]   6 < 13
  constitutive: τ symmetric (D08); linear + isotropic ⇒ τ = −pδ + 2μS + λS_mmδ (4.31) (D09 ★★★); μ_v, Stokes (D10) [C07]
     └─► Navier–Stokes (4.38) (D11) → (4.39a,b) (D12); viscous force −μ∇×ω (D13); Euler (4.41) [C08]      4 < 5
  noninertial frame: u = U + u′ + Ω×x′ (D14) → +2Ω×u′ … (4.43) (D15 ★★★) → forces −2Ω×u′, Ω²R e_R (4.45) (D16);
     pole projectile Ωut² (D17); effective gravity (D18) [C09]
  energy: total (4.53) (D19, D20) − mechanical (4.56) (D21) = internal (4.57) (D22); ε ≥ 0 (4.58) (D23); entropy (4.63) [C10]  7 = 7
  special forms: Lamb identity → ∇B = u × ω (4.69) (D24) → B const on streamlines and vortex lines (4.71) (D25) [C11];
     unsteady potential Bernoulli (4.75) with φ − ∫B dt′ (D26) [C12]; Boussinesq (4.86), (4.89) (D27, D28) [C13]
  edges: Dη/Dt = 0 (4.91) (D29), pillbox jumps, no-slip, Laplace jump derived [C14]
  scaling: dimensionless NS (4.101) (D30): St, Re, Fr, Fr′, Ri, M, Ec, Pr, We, Bo, Ca; similarity, C_D(Re) [C15]

next: ch05 vorticity dynamics — curl of (4.39b) using D13 (−μ∇×ω) and D24 (Lamb identity); Kelvin's circulation
      theorem (barotropic, from D24's pressure function); vorticity diffusion = ch03 Lamb–Oseen σ² = 4νt (exact solution
      already in `exact_solution`); rotating frames (D14–D16) → absolute vorticity ω + 2Ω; Boussinesq → baroclinic torque.
```
The book ahead: Ch5 vorticity · Ch6 ideal flow · Ch7 gravity waves · Ch8 laminar flow → Ch9 boundary layers → Ch10 CFD →
Ch11 instability → Ch12 turbulence → **Ch13 GFD (uses 4, 5, 7, 8, 11, 12)** → Ch14 aerodynamics · Ch15 compressible ·
Ch16 biofluids.
Where ch01 feeds in: N², θ → Ch. 7 §7.8, Ch. 11 §11.7, Ch. 13 (and ch04 Fr′, Ri); Newtonian τ, ν → Ch. 8, 9; hydrostatic
base → Ch. 7, 13; Π groups → ch04 §4.11 done; isentropic gas → Ch. 15; FTCS → Ch. 5, 8, 10.
Where ch02 feeds in: principal axes → Ch. 12; Gauss → every CV; Stokes/circulation → Ch. 5 Kelvin, Ch. 6 lift, Ch. 13 PV;
invariants → Ch. 12; ε → every cross product.
Where ch03 feeds in: vortices (Lamb–Oseen) → Ch. 5, 8, 13, 14; streamlines/ψ/cylinder → Ch. 6; path lines → Ch. 7 orbits,
Ch. 10; principal strain → Ch. 12, Ch. 13 frontogenesis; ω′ = ω − 2Ω → Ch. 13 absolute vorticity; RTT → Ch. 15 CVs.
Where ch04 feeds in: **NS residual/term tools → verify every exact solution of Ch. 5–13**; D13 + D24 → Ch. 5 vorticity
equation and Kelvin; ψ tools → Ch. 6, 8, 9, 13; unsteady Bernoulli + kinematic BC → Ch. 7 free surface; exact
solutions → Ch. 8; CV budgets → Ch. 9 momentum integral, Ch. 14, 15; ε and energy budgets → Ch. 12; rotating frame,
Boussinesq, Ri, Fr′, Ro → **Ch. 13**; stagnation quantities, μ_v, M → Ch. 15.

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
| `_stencil` (ch04) | private explicit 2nd-order central differences on callables (`ddx`, `d2dx2`, `ddt`, `grad`, `div`, `laplacian`) | residual functions |
| `tools/convergence.observed_order(h, err)` | log–log slope | every V3 test |

### Physics primitives from ch01 (re-exported as `ch01.<name>`)
| function family | module | book § / Eq. | label |
|---|---|---|---|
| constants `K_B, N_A(_KMOL), R_U, M_W_AIR, R_AIR, GAMMA_AIR, CP_AIR, CV_AIR, G0, P_ATM, P_REF` | `thermo` | §1.9 | benchmark |
| perfect gas, specific heats, first law paths, entropy, Gibbs, sound speed, isentropic ratios, van der Waals; `helmholtz_free_energy` (ch04) | `thermo` | (1.10)–(1.28), (4.94) | analytic, symbolic, converged, conserved, benchmark |
| `integrate_hydrostatic`, `standard_atmosphere(geometric=)`, `scale_height`, `buoyancy_force`, `net_pressure_force_on_box` | `statics` | (1.8), (1.9), §1.10 | analytic, converged, benchmark |
| `brunt_vaisala_sq(_from_lapse/_from_theta)`, `classify_stability`, `parcel_displacement`, `adiabatic_lapse_rate` (negative), `lapse_rate_convention`, `lapse_rate_stability`, `potential_temperature`, `potential_density` | `stratification` | (1.29)–(1.35) | analytic, symbolic, benchmark |
| `ftcs_diffusion_1d`, `stable_time_step`, `couette_startup_profile`, `gaussian_spreading` | `diffusion` | §1.5 | converged, conserved, analytic |
| `dimensional_matrix`, `rank_by_minors`, `pi_groups`, `rescale_units`, presets | `dimensional` | (1.36)–(1.40) | symbolic, analytic |

### Mathematical and field primitives from ch02 (`ch02.<name>`)
| function family | module | book § / Eq. | label |
|---|---|---|---|
| `expand_indices`, `classify_indices`, `rename_dummy`, `comma_to_partial`; `coordinates(dim)` | `index_notation` | (2.2), (2.36) | analytic, symbolic |
| `direction_cosines`, `rotation_matrix_2d/3d` (passive C), `transform_vector` (Cᵀu), `transform_tensor`, `traction(tau, n)` (first index), `normal_shear_stress`, `mohr_circle_2d`, `double_dot(convention=)`, `invariants`, `levi_civita`, `cross`, `strain_rate_tensor(G)`, **`rotation_tensor(G)` = G − Gᵀ**, `principal_axes` | `tensors` | (2.1)–(2.29), §2.11 | analytic, symbolic, benchmark |
| `grid`/`grid2d` (`[k, j, i]`), `partial`, `gradient`, `divergence`, `vector_gradient` (G[i, j] = ∂u_i/∂x_j), `tensor_divergence(index=1 default; **index=0 for Cauchy**)`, `curl`, `laplacian`, `is_solenoidal`, `is_irrotational` | `grids`, `operators` | (2.22)–(2.25) | analytic, converged |
| `VectorField`, `ScalarField` | `fields` | Ex. 2.3 | analytic, symbolic |
| `divergence_theorem_box/_sphere/_tiled`, `integral_*`, `circulation`, `curl_flux`, `stokes_theorem_check`, `gauss_legendre_nodes` | `integral_theorems` | (2.30)–(2.35) | analytic, converged, conserved, benchmark |

### Kinematics primitives from ch03 (`ch03.<name>`)
| function family | module | book § / Eq. | label |
|---|---|---|---|
| field callables `u(x, t)` (x `(d,)`/`(d, N)`, components on axis 0); `material_derivative(_terms, _sym)`, `acceleration`, `lagrangian_*`, `streamline`, `pathline`, `streakline`, `galilean_transform`, `velocity_gradient_at`, `vorticity`, strain/spin rates, `vorticity_in_rotating_frame`, `principal_strain_rates`, `strain_ellipse_axes`, linear-flow presets | `kinematics` | (3.1)–(3.21) | analytic, symbolic, converged |
| polar/cylindrical/spherical conversions (θ from +z), unit vectors, velocity components | `coords` | §3.1 | analytic |
| `solid_body_rotation`, `line_vortex`, `rankine_vortex`, `gaussian_vortex`, `gaussian_vortex_max_radius`, `polar_vorticity_z(_sym)`, `circulation_circle`, `vortex_velocity_field` | `vortices` | (3.22)–(3.29) | analytic, symbolic, converged, benchmark |
| `leibniz_terms`; CVs `MovingBox`, `GrowingSphere`, `GrowingCylinder`, `GrowingCone`, `MovingEllipse2D`; `reynolds_transport`, `rtt_check`, `swept_terms`, `material_volume_rate` | `transport` | (3.30)–(3.35) | analytic, symbolic, converged, conserved |

### Conservation-law primitives from ch04 (`ch04.<name>`; details `knowledge/ch04.md` §3)
| function family | module | book § / Eq. | label |
|---|---|---|---|
| `mass_budget`, `interval_mass_budget`, `momentum_budget`, `energy_budget`, `angular_momentum_budget` → dataclasses with every term and the residual, on any ch03 `ControlVolume`; flux (u − b)·n | `conservation` | (4.5), (4.17), (4.48), (4.65) | analytic, conserved, converged, symbolic |
| `newtonian_stress(G, p, mu, lam|mu_v)`, `viscous_stress`, `isotropic_fourth_order`, `linear_stress`, `mean_pressure(tau, tau33)`, `bulk_viscosity`, `stress_on_plane`, `deviatoric_part`, `dissipation_rate` | `constitutive` | (4.25)–(4.37), (4.58) | analytic, symbolic |
| **residuals and term splits (numeric + sympy)**: continuity, Cauchy (first index), NS (4.38)/(4.39b), `viscous_force_forms`, `exact_solution` (8 solutions), `ns_terms_preset`, Lamb identity, energy/entropy family, Boussinesq family | `navier_stokes` | (4.7)–(4.89) | analytic, symbolic, converged, conserved |
| Appendix-B operators (cylindrical, spherical), sympy | `curvilinear` | App. B | analytic, symbolic |
| ψ ↔ u (planar u = ∂ψ/∂y, axisymmetric), fluxes, χψ, presets | `streamfunction` | (4.11), (4.12) | analytic, symbolic, converged |
| `frame_acceleration_terms` (+2Ω×u′), `apparent_body_forces` (−2Ω×u′), `coriolis_force`, `centrifugal_*`, `effective_gravity`, `coriolis_parameter`, `projectile_paths`, `OMEGA_EARTH` | `rotating` | (4.42)–(4.45) | analytic, symbolic, converged, benchmark |
| `bernoulli_head/solve`, `pressure_function`, `bernoulli_function`, `which_bernoulli(_text)`, unsteady family, `gauge_absorbed_bracket`, stagnation, pitot, Torricelli | `bernoulli` | (4.19), (4.66)–(4.83) | analytic, symbolic, conserved, benchmark |
| `kinematic_bc_residual`, `surface_normal_speed`, `interface_mass_flux`, `pillbox_limit`, cap forces, `laplace_jump_from_balance`, `capillary_length` | `interfaces` | (4.90)–(4.98) | analytic, converged, benchmark |
| `Scales`, 20+ named numbers (Re, Fr, St, Fr′, Ri, Ri_g, M, Ec, Pr, We, Bo, Ca, Ro), C_p/C_D/C_L, `sphere_drag_coefficient` (Morrison), `model_prototype`, `nondimensional_*_coefficients` | `similarity` | (4.99)–(4.119) | analytic, symbolic, benchmark |
Chapter-only code stays in chapter modules (move to `core/` when a second chapter calls it): ch01 (continuum sampling,
Kn, μ(T), seawater EOS, examples); ch02 (polar components, examples, test fields); ch03 (`cylinder_flow`/
`cylinder_streamfunction` → Ch. 6, `pipe_profile` → Ch. 8, `thermal_front`, examples); **ch04** (`linear_wave_surface`,
`bore_speed` → Ch. 7; `couette_heating(_transient)`, `two_fluid_couette`, `navier_slip_couette`, `stokes_first_problem`
→ Ch. 8; `wake_drag_per_span` → Ch. 9; `accelerating_sphere_*`, `u_tube_column` → Ch. 6/7; Boussinesq scenarios, meniscus,
ship model). Note `ch04.G` = 9.81 while `core` defaults to `G0` = 9.80665.

### Explainer engine (`assets/viz_lib.js`)
`Viz.app` (tabs Walkthrough / Explore / Explain / Derivation / Equations / Code / Check; fit-to-window with density
levels and pagers; transport, modes, presets, status, terms, inspector, notes, selftest parity), `Plot`, `Viz.field`,
`Viz.num` (RK4, odeint, brentq, erf/erfc **~1.2e-7 only**, niceTicks), `Viz.work`, `Viz.three`, KaTeX with fallback;
`Viz.font`, `Viz.roundRect`, `Viz.text`, `Viz.card`, `Viz.fmtTime`, `Viz.tnum(keepTiny)` (ch01 pass); pager slicing,
`--viz-stage-need`, `ev.viewId`/`ev.vizView` (ch03). **Nothing promoted after ch02, ch03 or ch04** (the site-publisher
was reading `viz/` each time). Ranked candidates (ch04 list: formatters in all 26 explainers, waterfall bars 9,
per-mode control hiding 7, field image 5, principal2d 5, erf/erfc precision, `vec3`, …), four open library bugs (ch03)
and eight quirks Q1–Q8 (ch04: shared Code-tab placeholder namespace, badge over the first view, fixed transport max,
no per-mode controls, `onChange` overwriting `set()`, chips not wrapping, erf precision, `eq_refs` false hits) are in
`knowledge/viz_patterns.md`. Reference: `templates/viz_example.html`.

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
| ch04 | `control_volume_budgets` | C01, C04 (D01, D05) | five scenes on one budget object (wake, bore, jet, rocket, balloon) · face flux bars + measured ◇ · "drag b until storage vanishes" |
| ch04 | `stream_function_spacing` | C03 (D04) | ψ at equal Δψ over a speed map · spacing = speed bars · bendable gate with flux Δψ · axisymmetric mode |
| ch04 | `newtonian_stress_lab` | C07 (D08, D09 ★★★, D10) | G → S, R → τ grid coloured by source → rotatable plane traction · cube spin-up mode α ∝ 1/h² |
| ch04 | `navier_stokes_term_balance` | C08 (D11–D13) | 7 exact solutions, click a point, term bars summing to 0, regime word, time scrub |
| ch04 | `rotating_frame_coriolis` | C09 (D14–D18, D15 ★★★) | two observers on one clock · signed side toggle forces/accelerations · split Coriolis arrow · f-plane highs/lows |
| ch04 | `which_bernoulli` | C05, C11, C12 (D06, D24–D26) | six flows with probes along/across streamlines · hypothesis decision table · stacked B terms |
| ch04 | `viscous_dissipation_heating` | C10 (D21–D23) | Couette gap: profile, ε, T(y, t) · in = stored + out bars · "μ < 0 ⛔" preset |
| ch04 | `boussinesq_buoyancy` | C13 (D27, D28) | rising blob with buoyancy-off ghost · kept vs dropped bars · validity chart with lit settings table |
| ch04 | `dynamic_similarity_models` | C15 (D30) | prototype and model on one t* clock · paired group bars with badges · sphere collapse · Ro preview |

## Notation
See `knowledge/notation.md` (symbol register and sign traps). Conventions that matter everywhere: SI and kelvin inside
functions; **z up**; **lapse rate Kundu Γ ≡ dT/dz in code, meteorology −dT/dz shown alongside**; `p0` ≠ `p_ref`; q to /
w on the system; gas constants per kmol; pressures absolute unless `_gauge`; signed τ_xy = μ ∂u/∂y.
From ch02: passive C (x′ = Cᵀx); traction on the first index; book A:B = A_ij B_ji; G[i, j] = ∂u_i/∂x_j; Stokes n_c
into A; grid `[k, j, i]` with `h` in (x, y, z).
From ch03: R = G − Gᵀ (no ½), ω = ∇×u, spin ½ω; γ = 2S₁₂ vs ch02 Γ ≡ S₁₂ vs Γ circulation; Galilean u′ = u − U;
rotating ω′ = ω − 2Ω; RTT outward n, signed b·n; field callables `u(x, t)`; spherical θ from +z; line-vortex strength B.
From ch04:
- **Cauchy's divergence contracts the first index** ∂τ_ij/∂x_i (`tensor_divergence(index=0)`); the book's prose writes
  ∂τ_ij/∂x_j (harmless only because τ is symmetric).
- **Flux through a moving wall uses (u − b)·n**; drag is on the body, the fluid feels −F_D.
- **Coriolis acceleration +2Ω × u′ (4.43) vs force −2Ω × u′ (4.45)**; centripetal Ω × (Ω × x′) vs centrifugal
  −Ω × (Ω × x′) = Ω²R e_R; NH deflection to the right; Ro = U/(2Ωl) (forward pointer; Ch. 13 U/(fL)).
- **Stokes assumption μ_v = 0, not λ = 0**; p − p̄ = μ_v∇·u; plane stress needs τ₃₃.
- 2-D ψ: u = ∂ψ/∂y, v = −∂ψ/∂x (GFD texts often the opposite); axisymmetric u_R = −(1/R)∂ψ/∂z.
- **The prime has four meanings** (ch02 rotated axes, ch03 translating frame, ch04 rotating frame, ch04 perturbation;
  plus a dummy variable in (4.67)): one code name per meaning (`u_prime`/`u_rot`, `p_pert`).
- Two g defaults: `core` 9.80665, `ch04` 9.81 — pass `g=` where they meet.
Symbols with several meanings so far: α, θ, T, h, R, σ/S, c, D, n, k, q, γ, λ, τ, ε, H (ch01); Γ, ω, φ, A, C, G, K, m, s,
t, a, b (ch02); Γ (4 meanings), ω, σ, φ, θ, α, a, b, h, primes, U (ch03); primes (4), σ, ε, γ, Φ, Ψ, Ω (frame rotation vs
imposed frequency), h, H, M, R, f, F, B, C_p (pressure coefficient vs specific heat), λ, b, η, ζ, G (ch04).

## Teaching lessons
- **Depth is tiered; coverage is exhaustive.** ch01: 15 A / 70 B / 21 C, 12 derivations, 496 cells. ch02: 16 A / 60 B /
  18 C, 15 derivations, 405 cells. ch03: 15 A / 52 B / 12 C, 24 derivations (171 steps), 413 cells. ch04: 15 A / 151 B /
  19 C, 30 derivations (235 steps), 611 cells.
- **Convention callouts** (lapse rate, "Which p_o?", passive/active, first/second index, A:B, R vs A, γ vs Γ vs
  circulation, translating vs rotating frame, signed b·n, **Coriolis term vs force, four primes, ψ sign, λ vs μ_v**):
  state the book's, the field's and the code's convention, give the size of the difference in numbers, say which one
  the code uses and why.
- **Book typos are taught, not hidden**: "book prints X; Y is right" with a test that shows the printed form wrong (ch02
  Ex. 2.3/2.4/2.6; ch03 (3.6), Ex. 3.2; ch04 (4.15) "= 0", (4.51) dA, (4.74) gauge sign, curve C, Ex. 4.7, Ex. 4.2).
- **A running ledger across a chapter** (ch04 equations vs unknowns 0 → 6/13 → 4/5 → 5/5 → 7/7) makes the chapter's
  logic visible; **decision tables** sort results with different hypotheses (ch04's four Bernoulli forms).
- **Every key number by two independent routes** (formula vs measured/finite difference; brentq vs Lambert W; loops vs
  einsum vs library) — in the notebook and as a measured ◇ in the explainer.
- **★★★ sympy checks re-run the derivation's own construction** (ch01 D28, D19; ch02 D17; ch03 D22; ch04 D09 builds K
  from δ's and rotates it, D15 builds x = X + Rx′ and shows the factor 2 is necessary).
- **A derivation check must cite a cell that runs** (ch04 round 1: five did not); **every reading note of a figure is
  checked against printed numbers** (ch04 Lamb–Oseen "arrows forward outside" was false); explainer IDs ("E1") are
  jargon in the notebook — name the explainer; keep one scenario per name (ch04's three "lakes").
- **Make an approximation measurable** (ch03 ellipse ≈ γ/4; ch04 the Boussinesq "dropped" bar ×αδT and validity chart).
- **Every prose number is computed; every "Dxx step N" pointer is checked by script after a split.**
- **Primers come before the first derivation that uses them**; `knowledge/primers.md` lists 133 (P01–P133); ch05 starts
  at P134.
- **Climate hooks with numbers** (ch03 thermal advection, ω′ = ω − 2Ω; ch04 pole projectile 9.34 km, geopotential and
  effective gravity, Boussinesq ocean/atmosphere, Ri, Fr′, Ro preview, highs and lows).
- Reusable derivation moves (ch01–ch04, 63 so far; 21 from ch04) are tabulated in `knowledge/concept_map.md` → "Reusable derivation
  moves".

## Global pitfalls confirmed in this project
**Machinery and environment**
- Extracted text garbles maths (`¼` = `=`, `ð…Þ` = parentheses, missing minus, ω → `u`, γ → g, θ → q, σ/τ → s, ν → n,
  ε → 3, ρ → r, Ω → U, ψ → j, χ → c, φ → 4, Φ → F, η → h, ζ → z, ∂ → v, ∇ → V) — read equations from rendered pages.
- **Windows: shell heredocs mangle backslashes** — write code, design and explainer files with Write/Edit (ch04: a heredoc
  corrupted the design file; single-backslash TeX in E7's JS strings rendered `\rho` as CR + "ho", `\tfrac` as TAB).
  Python builder strings with LaTeX must be raw (coverage check 9). Scan explainers for single-backslash TeX.
- Anaconda's `python3` kernelspec can shadow the venv's — notebooks execute on `fluidpy-venv`.
- JS: `UIEvent.view` is read-only. `tools/shot.py` does not click/drag; a reviewer runs a click pass.
- `test_viz_library_inlined_and_template_lints` fails whenever `assets/viz_lib.js` changed or an explainer is mid-build;
  run `tools/viz_inline.py --all` before the merge gate.
- `tools/shot.py` `py:` parity expressions have **no builtins**; parity rows must call the explainer's own functions,
  never a JS literal.
- Parallel viz-builders use **private scratch subfolders**.
- **matplotlib 3.11**: animations using `subplots_adjust`/`fig.text` disable `figure.constrained_layout.use` for the save.
- `tools/coverage_check.py` matches ledger reminder rows by exact primer name → false warnings.
- `nb.derivation` prefixes "Eq." to its label — pass non-equation labels in the title.
- Book wording can hide in docstrings (`check_public.py` does not catch paraphrase): reviewers scan for it.
- Library: pager packs before the final `fit()`; tall ∫ clip at slice breaks (`\textstyle` rows); `Viz.num.erf/erfc`
  only ~1.2e-7 (local double-precision versions in ch04 E1, E4); transport range fixed per app; no per-mode controls.
- Docstring contract numbers and derivation IDs go stale after fixes — re-grep them (ch04 review Should-fix 3–4).

**Mathematics → code**
- **A test comparing two of our own functions is not evidence for the book's convention.** Pin conventions to fields
  with known answers and **prove discrimination by patching in the wrong variant** (ch01; ch02 8/8; ch03 13/13 + 7;
  ch04 26/26).
- **Test every scenario's full balance per direction, not just totals** (ch04: the oblique jet split 50/50 and the 2 × 2
  mean pressure passed 147 tests; the derivation review caught both).
- **Reduced-dimension inputs hide physics**: plane stress needs τ₃₃ for p̄ and the (3, 3) deviator for ε; a scalar where a
  vector is expected raises; 2-vectors are padded explicitly.
- **Series solutions**: no early stop at a zero coefficient (symmetry), no fixed-node projection (aliasing); closed-form
  coefficients (ch04 F1).
- **Preset names and letters are physics** — γ = 2S₁₂ vs Γ ≡ S₁₂ vs Γ circulation; λ vs μ_v; Coriolis + vs −; each pinned.
- **Signed quantities stay signed**: RTT sliver b·n, (u − b)·n, Leibniz lower term, wall stress, lapse-rate conventions,
  Stokes orientation, clockwise spin ω₃ = −γ, heat flux out of a CV, the (4.74) gauge sign.
- **Rotating ≠ Galilean**: ω′ = ω − 2Ω; a rotating frame adds four apparent accelerations with the factor 2.
- Index conventions: traction and Cauchy on the first index, ch02 `tensor_divergence` default on the second; passive C
  vs active R; book vs Frobenius double dot; G vs the transposed integral gradient; R = 2A; grid `[k, j, i]`.
- **Units of numerical defaults**: a time-step default is never the space step; two g defaults (9.80665 vs 9.81).
- Vectorised inputs: trace over axes (0, 1) for (d, d, N) arrays; `as_scalar_if_0d`. Closed loops wrap.
- sympy: assumptions make different symbols (use `coordinates`); `.norm()` adds `Abs`; `Subs(Derivative)` does not
  simplify — substitute generic polynomials; cache slow curvilinear simplifications.
- A default at a regime boundary is computed from the boundary function. kmol vs mol; pint °C offset; geopotential vs
  geometric altitude; `np.gradient` first order at edges (project stencils are explicit second order).
- Cancellation: `-np.expm1(-x)`; stiff geometry: integrate in arc length (ch04 meniscus). Unbounded ODEs: `solve_ivp`
  event, NaN beyond; burn-out events (rocket).
- **Validation labels**: V4 only for real conservation laws; unit/frame invariance is V7; identities V1; encyclopedia
  *form* cross-checks are V1, not V5 (ch03 Rankine/Lamb–Oseen/RTT; ch04 Bélanger, Tsiolkovsky, added mass,
  Taylor–Green); "converged" only with an asserted order; closed-form identities are V1, not V4 (ch04 `material_mass`);
  book-rounded constants only in the private JSON; property-table differences are reported, not asserted (ch04 O5).

## Benchmark inventory
| value / table | source (citation) | stored in | used by |
|---|---|---|---|
| k_B, N_A, R (exact) | CODATA 2018/2022, NIST CUU | `reference/ch01/constants.json` | ch01 `thermo` |
| USSA-1976 constants, Tables 1–2; sea-level c = 340.29 m/s | NASA-TM-X-74335; PDAS tables | `reference/ch01/ussa1976_*` | ch01, ch04 C02 |
| IAPWS σ(T), μ(T) | IAPWS R1-76(2014), R12-08 | `reference/ch01/iapws_sigma.csv`, `benchmarks.json` | ch01 C13, C15; ch04 capillary length |
| Jennings mean free path; Taylor blast K = 0.856; AMS Γ_d; EOS-80 check values | Tsalikis et al. 2024; Díaz arXiv:2009.05674; AMS Glossary; Fofonoff & Millard 1983 | `reference/ch01/benchmarks.json` | ch01 |
| Divergence theorem 8π/3; Levi-Civita identities; Mohr; rotation matrices; Stokes statement | Wikipedia; scipy docs | `reference/ch02/` | ch02 |
| Gaussian (Lamb–Oseen) vortex peak α = 1.256 | Canivete Cuissa & Steiner 2022, SWIRL I, arXiv:2210.05223 | `reference/ch03/` | ch03 C14 |
| **Sphere drag C_D(Re), 0.1–1e6** | F. A. Morrison (2016), Michigan Tech, Eq. (1) | `reference/ch04/benchmarks.json`, `SOURCES.md` | ch04 C15 (`sphere_drag_coefficient`, E9) |
| **WGS-84 a, 1/f, b, ω** (2(a − b) = 42.77 km) | Wikipedia "World Geodetic System" | `reference/ch04/` | ch04 C09 (`OMEGA_EARTH`, effective gravity) |
| **Sharp-orifice C_c = 0.611** | Wikipedia "Vena contracta" | `reference/ch04/` | ch04 C05 |
| **Capillary length of water 2.71 mm (20 °C)** | Wikipedia "Capillary length" | `reference/ch04/` | ch04 C14 |
| **Prandtl: air 0.70–0.73, water 5.9 (300 K), monatomic 2/3** | Wikipedia "Prandtl number" | `reference/ch04/` | ch04 C15 (`prandtl_of`) |
| Forms (V1 cross-checks, not V5): Rankine, Lamb–Oseen, RTT, Leibniz (ch03); Bélanger jump, Tsiolkovsky, added mass, Taylor–Green (ch04) | Wikipedia pages | `reference/ch03/`, `reference/ch04/` | ch03, ch04 |
Book-printed values (private, git-ignored): `tests/book_values_ch01.json` … `…_ch04.json` (Ex. 4.1, 4.3, 4.5–4.8,
cap forces, projectile, capillary scale, Mach threshold, Prandtl values, Earth bulge).

## Validation summary per chapter
Counts are concept-map rows carrying each label (a row can carry several); tests per evidence level in parentheses.
| chapter | analytic | symbolic | converged | conserved | benchmark | book-value | qualitative | unverified | verdict |
|---|---|---|---|---|---|---|---|---|---|
| ch01 | 39 (V1 58) | 16 (V2 19) | 5 (V3 8) | 5 (V4 7) | 18 (V5 13) | 1 (V6 5) | 1 (`wavelength_to_rgb`) | 0 | **PASS** — 107 tests; 5 explainers (115 rows); notebook 496 cells |
| ch02 | 43 (V1 41) | 21 (V2 13) | 11 (V3 15) | 2 (V4 1) | 9 (V5 6) | 8 (V6 4) | 0 | 0 | **PASS** — 92 tests (+ V7 12); 8/8 wrong variants; 5 explainers (127 rows); notebook 405 cells |
| ch03 | 38 (V1 58) | 23 (V2 32) | 14 (V3 14) | 5 (V4 4) | 2 (V5 1) | 8 (V6 4) | 0 | 0 | **PASS** — 120 tests (+ V7 7); 13/13 (+7) wrong variants; 7 explainers (145 rows); notebook 413 cells |
| ch04 | 41 (V1 83) | 35 (V2 41) | 17 (V3 9) | 13 (V4 6) | 11 (V5 6) | 5 (V6 3) | 0 | 0 | **PASS** — 151 tests (+ V7 3); 26/26 wrong variants; review 2 Must + 7 Should applied; every ★★/★★★ derivation re-derived in sympy; orders 1.99–2.02 (1.00 where first order by design); 9 explainers PASS round 2 (255 rows, 2 744 views, 30 derivations matched); notebook PASS round 2 (611 cells, 113 s) |
(ch01 V7: 19 tests; ch02: 12; ch03: 7; ch04: 3. Full suite at the ch04 merge gate: **478 tests**.)

## Open across chapters
- **Library pass (orchestrator, when nothing reads `viz/`)**: fix the four ch03 library bugs and the ch04 quirks Q1–Q8;
  promote the ranked JS candidates (`knowledge/viz_patterns.md`: formatters in all 26 explainers first; waterfall bars,
  per-mode controls, field image, principal2d, double-precision erf/erfc, `vec3` before Ch. 13); then
  `tools/viz_inline.py --all` and `tools/shot.py --chapter ch01…ch04 --quick` + `templates/viz_example.html --quick`.
- **Skill pass**: the lesson candidates for `interactive-viz`, `math-to-python` §7, `verify-implementation`,
  `teaching-style` and `colab-notebook` listed in `knowledge/viz_patterns.md` (ch02 + ch03 + ch04) are not yet in the
  skills — the ch04 "test every scenario's full balance" lesson is the most important.
- Machinery TODO: `shot.py` clicks/drags every view; `coverage_check` understands ledger reminder rows; `nb.derivation`
  skips "Eq." for non-numeric labels; `show_eqs` → `tools/nbkit.py`; `core.anim` saves inside a constrained-layout-off
  `rc_context`; `eq_refs` exempts selftest names and `<meta>`.
- Non-blocking explainer follow-ups: ch01 E1–E5, ch02 E1/E3/E5, ch03 E1–E3/E6/E7 (their `knowledge/chNN.md` §9); ch04
  long 360×640 walkthrough cards in 6 explainers, E3/E8 Explain 28/24 pages, E5 literal parity row, E6 meta, E9 ship
  label (`knowledge/ch04.md` §9).
- ch04 O5: the book's Prandtl numbers vs our property tables (−1.5 %, −2.2 %), reported not asserted.
