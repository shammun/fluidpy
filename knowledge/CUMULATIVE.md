# CUMULATIVE knowledge — fluidpy
(Rewritten by the knowledge-keeper after every chapter. Every agent reads this first. Last rewrite: after ch03,
2026-09-23, commit 1f166f6.)

## Project rules added during ch03 (they apply to every chapter from now on)
- **5–10 explainers per chapter, as many as the CORE ideas need** (was 4–5; `book.yaml → project.min/max_explainers_
  per_chapter`, CLAUDE.md, skills, agents and `coverage_check` check 4 updated; commit 0792849). ch03 has 7.
- **Every book equation is shown in full next to its number** — in notebook prose, derivation steps, traps, recaps and
  in explainer tour/Explain/Derivation/quiz/notes/status text (commit 33df5c8). Enforced by `tools/coverage_check.py`
  **check 8** (a markdown cell naming "(3.5)" or "Eq. 3.5" must contain LaTeX; the reviewer checks it is the RIGHT
  equation) and by **`tools/eq_refs.py`** for explainers (lists JS strings citing a number without TeX; `ref:` labels next
  to shown TeX, metadata and selftest names are exempt). Builders use a `show_eqs(text)` helper
  (`notebooks/build_ch03.py`; JS `showEqs` in ch03 E5) that writes the equation after its first bare mention.
- **`coverage_check` check 9**: no control characters (`\f`, `\b`, `\v`, `\a`) in any cell — the sign of a LaTeX
  backslash eaten by a non-raw Python string (`"\frac"` → form feed). Write LaTeX in raw strings.
- **ch01 and ch02 were retrofitted** to the equation rule: explainers (a54f4c3; review `reports/eq_retrofit_viz.md`,
  fixes 7a9f7ca, cb68262, 266a1ca) and notebooks (d212904: 58 + 77 cells). All 10 ch01–ch02 explainers PASS the
  every-page audit.
- **Pager and audit (68f214e)**: an item taller than a page is sliced at natural breaks with a "continued" cue;
  `tools/shot.py` pages through **every page of every pager** at every size (`CLIP__<size>__<view>__pgK.png`).

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
  summation convention [C01] → direction cosines, CCᵀ = I (D02, D03) [C02] → x' = Cᵀx (D01) [C03]
  stress τ_ij [C04] → Cauchy f_i = τ_ji n_j (D05) [C05] → τ' = CᵀτC (D06) [C06] → invariants (D18) [C07]
  ε_ijk, epsilon–delta (D09) [C08]; S + A split, R = G − Gᵀ ↔ ω (D14, D15) [C12]; principal axes (D17 ★★★) [C13]
  fields on the grid [k, j, i]: ∇φ [C09], ∇·u [C10], ∇×u (D12) [C11]
  Gauss (D25) [C14] → integral definitions (D21, D22) [C15]; Stokes, circulation per area (D26) [C16]

ch03 Kinematics (done) — how fluid moves, before any force
  particle r(t; r_o, t_o) ⇄ field F(x, t) via (3.2) (D01) [C01]
     └─► material derivative DF/Dt = ∂F/∂t + u·∇F (3.5) (D02) [C02]  ── Ch. 4 writes every law with it
  flow lines: streamline (3.7) (D03) [C03] · path line (3.8), streak line, Ex. 3.1 (D04, D05) [C04]
  frames: Galilean invariance of Du/Dt (3.9) (D06) [C05]; rotating ω′ = ω − 2Ω (D13)
  small element: du = G·dx (3.10) (D07) [C06] ─► n·S·n (D08) [C07] · S₁₂ = ½D(α+β)/Dt, rigid ⇒ S = 0 (D09, D10) [C08]
     · ∇·u = volume growth rate (3.14) (D11) [C09] · spin = ½ω (D12, D14) [C10] · du = S·dx + ½ω × dx (3.19) (D15) [C11]
     · principal axes, sphere → ellipsoid (D16) [C12]
  plane vortices: ω_z in polar form (3.23) (D17) [C13]; solid body vs line vortex; Rankine (D18), Gaussian (D19, D20) [C14]
  Leibniz (3.30) (D21) → Reynolds transport theorem (3.35) (D22 ★★★) [C15] → (3.14) (D23), (3.5) + F∇·u (D24)

next: ch04 conservation laws / Navier–Stokes (the trunk): continuity from RTT + (3.14); Cauchy's equation from RTT +
      ch02 Cauchy traction; Newtonian stress from S (only S, (3.19)); rotating frames (Coriolis); Bernoulli; Boussinesq
      with ch01's N²; similarity with ch01's Π groups.
```
The book ahead: **Ch4 conservation laws (the trunk)** → Ch5 vorticity · Ch6 ideal flow · Ch7 gravity waves · Ch8 laminar
flow → Ch9 boundary layers → Ch10 CFD → Ch11 instability → Ch12 turbulence → **Ch13 GFD (uses 4, 5, 7, 8, 11, 12)** →
Ch14 aerodynamics · Ch15 compressible · Ch16 biofluids.
Where ch01 feeds in: N², θ → Ch. 4 Boussinesq, Ch. 7 §7.8, Ch. 11 §11.7, Ch. 13; Newtonian stress, ν → Ch. 4 §4.5, 8, 9;
hydrostatic base state → Ch. 4, 7, 13; Π groups → Ch. 4 §4.11; isentropic gas → Ch. 15; FTCS → Ch. 5, 8, 10.
Where ch02 feeds in: tensors, Cauchy, principal axes → Ch. 4 §4.3–4.5; tensor divergence → Cauchy's equation; Gauss →
every control volume; Stokes/circulation → Ch. 5 Kelvin, Ch. 6 lift, Ch. 13 PV; invariants → Ch. 12; ε → Coriolis.
Where ch03 feeds in: D/Dt and RTT → **Ch. 4 §4.1–4.4 (continuity, momentum, energy for fixed/moving/material volumes)**;
(3.19) "only S deforms" → Ch. 4 §4.5 Newtonian stress; ω′ = ω − 2Ω and the Galilean/rotating contrast → Ch. 4 §4.7,
Ch. 13 (relative vs absolute vorticity, PV); vortices (Lamb–Oseen σ² = 4νt) → Ch. 5, 8, 13 (cyclones), 14 (tip
vortices); streamlines, stream function, cylinder flow → Ch. 6; path lines → Ch. 7 orbits, Ch. 10 particle tracking;
principal strain axes → Ch. 12, Ch. 13 frontogenesis; thermal advection −u·∇T → Ch. 4 energy, Ch. 13.

## Available primitives
### Machinery (`fluidpy/core/`, ready before chapter 1)
| module | what | used for |
|---|---|---|
| `project` | repo root, `book.yaml` config, Pages/raw/Colab URLs | every tool and notebook |
| `style` | `setup_notebook()`, `COLORS`, `savefig` | every notebook figure |
| `embed` | `show_viz(chapter, slug)` — Jupyter (srcdoc) / Colab (Pages src) / page | explainer cells |
| `anim` | `animate`, `show_animation(player="video"/"frames")` | animations (matplotlib 3.11 constrained-layout trap, see pitfalls) |
| `interact` | `slider_figure`, `animate_figure` (plotly, work on the page), `live` (ipywidgets) | Python interactives |
| `units` | pint `ureg`/`Q_`, `dimensional_check`, `check_dimensions`, `DIM`, °C ↔ K | V2 tests |
| `refdata` | reference-data registry | benchmarks |
| `_util` | `as_scalar_if_0d` (scalar in → float out, for parity rows), `require_positive`, `require_nonnegative` | every physics function |
| `tools/convergence.observed_order(h, err)` | log–log slope | every V3 test |

### Physics primitives from ch01 (re-exported as `ch01.<name>`)
| function family | module | book § / Eq. | label |
|---|---|---|---|
| constants `K_B, N_A(_KMOL), R_U, M_W_AIR, R_AIR, GAMMA_AIR, CP_AIR, CV_AIR, G0, P_ATM, P_REF` | `thermo` | §1.9 (CODATA, USSA-1976) | benchmark |
| perfect gas state, specific heats, first law paths, entropy, Gibbs, sound speed, isentropic ratios, van der Waals | `thermo` | (1.10)–(1.28) | analytic, symbolic, converged, conserved, benchmark |
| hydrostatics `integrate_hydrostatic`, `standard_atmosphere(geometric=)`, `scale_height`, `buoyancy_force` | `statics` | (1.8), (1.9), §1.10 | analytic, converged, benchmark |
| `brunt_vaisala_sq(_from_lapse/_from_theta)`, `classify_stability`, `parcel_displacement`, `adiabatic_lapse_rate` (negative), `lapse_rate_convention`, `lapse_rate_stability`, `potential_temperature`, `potential_density` | `stratification` | (1.29)–(1.35) | analytic, symbolic, benchmark |
| `ftcs_diffusion_1d`, `stable_time_step`, `couette_startup_profile`, `gaussian_spreading` | `diffusion` | §1.5 | converged, conserved, analytic |
| `dimensional_matrix`, `rank_by_minors`, `pi_groups`, `rescale_units`, presets `PIPE, PENDULUM, …` | `dimensional` | (1.36)–(1.40) | symbolic, analytic |

### Mathematical and field primitives from ch02 (re-exported as `ch02.<name>`)
| function family | module | book § / Eq. | label |
|---|---|---|---|
| `expand_indices`, `classify_indices`, `rename_dummy`, `comma_to_partial`; `coordinates(dim)` (real symbols — use these) | `index_notation` | (2.2), (2.36) | analytic, symbolic |
| `direction_cosines`, `rotation_matrix_2d/3d` (passive C), `transform_vector` (Cᵀu), `transform_tensor`, `traction(tau, n)` (first index), `normal_shear_stress`, `mohr_circle_2d`, `contract`, `double_dot(convention=)`, `invariants`, `levi_civita`, `cross`, `symmetric_part`, `antisymmetric_part`, `strain_rate_tensor(G)`, **`rotation_tensor(G)` = G − Gᵀ**, `vector_from_antisymmetric`, `principal_axes` (λ ascending, det +1) | `tensors` | (2.1)–(2.29), §2.11 | analytic, symbolic, benchmark |
| `grid`/`grid2d` (`[k, j, i]`, `h` in (x, y, z)), `partial`, `gradient`, `divergence`, `vector_gradient` (G[i, j] = ∂u_i/∂x_j), `tensor_divergence` (second index), `curl`, `laplacian`, `is_solenoidal`, `is_irrotational` | `grids`, `operators` | (2.22)–(2.25) | analytic, converged (orders 1.97–2.03) |
| `VectorField(exprs, coords, params, singular_at)`, `ScalarField` | `fields` | Ex. 2.3 | analytic, symbolic |
| `divergence_theorem_box/_sphere/_tiled`, `integral_gradient/divergence/curl`, `circulation`, `curl_flux`, `stokes_theorem_check` → `hypothesis_ok`, loops/surfaces, `gauss_legendre_nodes`, `midpoint_nodes` | `integral_theorems` | (2.30)–(2.35) | analytic, converged, conserved, benchmark |

### Kinematics primitives from ch03 (re-exported as `ch03.<name>`)
| function family | module | book § / Eq. | label |
|---|---|---|---|
| field callables `u(x, t)` (x `(d,)` or `(d, N)`, components on axis 0); `as_coord_field`, `from_coord_field` adapters to ch02's `fn(X, Y[, Z])` | `kinematics` | project convention | analytic |
| `lagrangian_velocity_acceleration(_sym)`, `lagrangian_to_eulerian` | `kinematics` | (3.1), (3.2) | analytic, symbolic, converged |
| `material_derivative`, `material_derivative_terms` → (local, advective, total), `material_derivative_sym`, `streamwise_derivative` (\|u\| ∂F/∂s) | `kinematics` | (3.3)–(3.6) | analytic, symbolic, converged |
| `streamline(u, x0, t_frozen, s_max)`, `pathline(u, r0, t0, t_eval)`, `streakline(u, x0, t, t_release)` | `kinematics` | (3.7), (3.8) | analytic, conserved |
| `galilean_transform(u, U, x0p)`, `acceleration(u, x, t)` → (a, local, advective) | `kinematics` | (3.9) | analytic, symbolic |
| `velocity_gradient_at`, `relative_velocity`, `vorticity(_from_gradient)`, `linear_strain_rate(G, n)`, `shear_strain_rate(G, n1, n2)` (= n₁·S·n₂), `volumetric_strain_rate`, `material_volume_ratio` (e^{t tr G}), `measured_strain_rates` | `kinematics` | (3.10)–(3.14) | analytic, converged, symbolic |
| `element_rotation_rate` (½ω₃), `material_line_rotation_rate`, `perpendicular_pair_rotation_rate`, `vorticity_in_rotating_frame` (ω − 2Ω), `velocity_potential_2d`, `relative_velocity_split` | `kinematics` | (3.15)–(3.19) | analytic, symbolic |
| `principal_strain_rates` (ascending), `strain_velocity_principal`, `strain_ellipse_axes(method="first_order"/"exact")`, `deform_circle`, `deform_sphere`; **promoted from ch02**: `velocity_gradient_preset`, `linear_flow_map`, `deform_square`, `material_line_angle` | `kinematics` | (3.20), (3.21), §2.10 | analytic, converged |
| polar/cylindrical/spherical conversions (θ from +z), `unit_vectors_*` (rows), `velocity_components`, `cartesian_components` | `coords` | §3.1, Fig. 3.3 | analytic |
| `solid_body_rotation`, `line_vortex(r, B)`, `rankine_vortex(r, Γ, σ)`, `gaussian_vortex(r, Γ, σ)` → (u_θ, ω_z), `gaussian_vortex_max_radius`, `vortex_profile(kind, …)`, `polar_vorticity_z(_sym)`, `circulation_circle`, `mean_vorticity_in_disc`, `vortex_velocity_field` | `vortices` | (3.22)–(3.29) | analytic, symbolic, converged, benchmark |
| `leibniz_terms`, `leibniz_check`; CVs `MovingBox`, `GrowingSphere`, `GrowingCylinder`, `GrowingCone`, `MovingEllipse2D` (`volume_nodes`, `surface_nodes` → X, outward n, dA, b; `volume`, `volume_rate`); `reynolds_transport`, `rtt_check`, `volume_integral(_rate_fd)`, `surface_flux_term`, `volume_rate_term`, `swept_terms(_sphere)`, `material_volume_rate`, `rtt_ellipse_2d` | `transport` | (3.30)–(3.35) | analytic, symbolic, converged, conserved |
Chapter-only code stays in the chapter modules: ch01 (continuum sampling, Kn, surface tension, μ(T), seawater EOS,
Examples 1.2–1.5); ch02 (polar components, Examples 2.2–2.4, test fields, `operator_convergence`); ch03
(`cylinder_flow`/`cylinder_streamfunction`/`cylinder_velocity_field`, `frame_acceleration_terms`, `thermal_front(_terms)`,
`unsteady_flow_preset`, `example_3_1`, `example_3_2`, `parallel_shear_kinematics`, `annular_sector_circulation`,
`rtt_field`, `leibniz_example`, `cross_section_average`, `pipe_profile`, `flux_through_disc`, `rigid_body_velocity`,
`potential_velocity`). **Move to `core/` when a second chapter calls it** — Ch. 6 will call the cylinder flow, Ch. 8 the
pipe profile.

### Explainer engine (`assets/viz_lib.js`)
`Viz.app` (tabs Walkthrough / Explore / Explain / Derivation / Equations / Code / Check; fit-to-window with density
levels and pagers; transport, modes, presets, status, terms, inspector, notes, selftest parity), `Plot`, `Viz.field`,
`Viz.num` (RK4, odeint, brentq, erf, niceTicks), `Viz.work` (Explain builders), `Viz.three`, KaTeX with fallback.
After ch01: `Viz.font`, `Viz.roundRect`, `Viz.text` (halo/bg), `Viz.card`, `Viz.fmtTime`, `Viz.tnum` + `keepTiny`.
During ch03 (68f214e): pager slices oversized items at natural breaks; `P.tall`/`data-room`; portrait stage floor
`--viz-stage-need` with `data-paged`; live notes/inspector re-pack; shot.py audits every pager page. Pointer events:
`ev.viewId`, `ev.vizView` (never assign `ev.view`). **Nothing promoted after ch02 or ch03** (the publisher was reading
`viz/`). The ranked candidates (fixed formatters in 11 explainers, per-layout row hide in 10, waterfall bars in 5,
`principal2d` 4, `expm2` 3, …) and **four open library bugs** (no keep-with-next; page packed at a different height than
shown; tall ∫ clipped at a slice break; `optional` controls hidden on phones while driving the picture) are in
`knowledge/viz_patterns.md`. Reference: `templates/viz_example.html`.

## Explainer inventory
| chapter | slug | CORE idea | reusable stage pattern |
|---|---|---|---|
| ch01 | `continuum_averaging_volume` | C06 continuum (+Kn, D35) | molecules · value vs log box size with band · regime strip; click-a-sample inspector |
| ch01 | `viscosity_momentum_diffusion` | C12 Newton viscosity, ν | gap with tracers · profile with steady ghost and trail · wall stress; modes = different diffusivity |
| ch01 | `heat_work_paths` | C25, C35, C45 | piston · p–v · T–s linked on a path parameter; hatched areas; path vs state term bars |
| ch01 | `parcel_stability` | C50, C51, C54, C55 | column · profile with adiabat · ζ(t) with linear ghost; two-convention badge with exact-text parity |
| ch01 | `buckingham_pi_machine` | C64, C67, C69 | variable chips → matrix with minor → groups; exact rational JS; unit-system switch |
| ch02 | `rotation_of_axes` | C02, C03, C06 (D01, D02, D06) | plane with fixed arrow and turning axes · clickable C matrix · τ′(θ); passive/active table; counter-example table |
| ch02 | `cauchy_traction_principal_axes` | C04, C05, C13 (D05, D17 ★★★) | element with a cut, f split blue/rose · σn/τs vs φ · Mohr; term bars; non-symmetric contrast |
| ch02 | `strain_vs_rotation_split` | C12 (D14, D15) | one clock, three squares (G, S, A); `expm2`; Frobenius bars with S:A ≡ 0 |
| ch02 | `gauss_flux_box` | C14, C15 (D25, D22, D21) | heatmap + draggable box, signed face arrows · waterfall flux bars · limit view vs log h; shrink transport |
| ch02 | `stokes_circulation_loop` | C11, C16 (D26, D12) | curl image + paddle wheel + draggable loop · unrolled u·t · side bars; Stokes failing on purpose |
| ch03 | `flow_lines_unsteady` | C03, C04 (D03, D04, D05) | unsteady flow with streamline pattern, particle trail and dye on one clock · answer-sheet view of the closed forms · port-velocity signal; threshold presets |
| ch03 | `material_derivative_probe` | C01, C02 (D01, D02) | fixed probe + carried float on one clock · local/advective/total bars with the float's measured ◇ · thermometer series; front / stretching-map modes |
| ch03 | `galilean_frames_cylinder` | C05 (D06) | cylinder flow in the observer's frame (observer slider lake → body) · bars that trade places with a dashed body-frame line · u = U + u′ triangle |
| ch03 | `fluid_element_deformation` | C06–C09 (D07, D08, D09, D11) | four G sliders deform a square + tracer ring · rate curves with measured dots · area vs e^{t tr G} |
| ch03 | `spin_and_principal_axes` | C10–C12 (D12–D16) | threads + paddle wheel in a chosen flow · single-thread rates vs flat pair average · split view with ellipse table; co-rotating observer; "Rate k" |
| ch03 | `vortex_paddle_wheels` | C13, C14 (D17–D20) | vortex with orbiting wheels and a draggable loop · u_θ, ω_z, Γ(r) profiles · log–log mean vorticity; real-vortex table |
| ch03 | `reynolds_transport_cv` | C15 (D21, D22 ★★★, D23) | moving CV with swept band by sign of b·n · budget waterfall + measured ◇ · dropped-term log–log; 1-D / ellipse / cone modes |

## Notation
See `knowledge/notation.md` (symbol register and sign traps). Conventions that matter everywhere: SI and kelvin inside
functions; **z up** (ch01); **lapse rate Kundu Γ ≡ dT/dz in code, meteorology −dT/dz shown alongside**; `p0` ≠ `p_ref`;
q to / w on the system; gas constants per kmol; pressures absolute unless `_gauge`; signed τ_xy = μ ∂u/∂y.
From ch02: passive C (x′ = Cᵀx); traction on the first index, tensor divergence on the second; book A:B = A_ij B_ji
(pass `convention=`); G[i, j] = ∂u_i/∂x_j (integral gradient is the transpose); Stokes n_c into A, t counterclockwise
about n; grid `[k, j, i]` with `h` and components in (x, y, z); λ^k a label.
From ch03:
- **R = G − Gᵀ (no ½), ω = vector(R) = ∇×u, an element spins at ½ω**; ω₃ = R₂₁ in 2-D; solid body ω = 2ω₀.
- **Γ has four meanings**: ch01 lapse rate; `velocity_gradient_preset` rate (simple shear Γ = du₁/dx₂ = ch03 γ = 2S₁₂);
  ch02 Ex. 2.4 Γ ≡ S₁₂; **Γ = circulation [m²/s] already in ch03** (`core/vortices.py`). γ = du₁/dx₂ = 2S₁₂.
- **Galilean (3.9)**: x = x′ + Ut + x′_o, u′ = u − U; primes = translating frame (ch02 primes = rotated axes). Towed
  cylinder: body frame steady (fluid at +U), fluid frame unsteady (cylinder at −U). Rotating: ω′ = ω − 2Ω.
- **RTT**: outward n, surface velocity b, **signed** b·n; Leibniz lower term subtracted.
- **(3.6) is DF/Dt = ∂F/∂t + |u| ∂F/∂s** (book drops F). Ex. 3.2: b·n = 0 on the base (not b = 0); "[?]" = z.
- Field callables `u(x, t)` with components on axis 0; spherical θ from +z; plane polar θ from +x; counterclockwise
  positive (shear spins clockwise, ω₃ = −γ); line-vortex strength B (book's Fig. 3.16 writes C; ch02 wrote K).
- ω also means Ex. 3.1's angular frequency (and Ch. 7's); σ = vortex core radius (ch01: surface tension); φ = azimuth
  and velocity potential; explainers name a mode-dependent rate "Rate k" with a meaning line.
Symbols with several meanings so far: α, θ, T, h, R, σ/S, c, D, n, k, q, γ, λ, τ, ε, H (ch01); Γ, ω, φ, A, C, G, K, m, s,
t, a, b (ch02); Γ (4 meanings), ω, σ, φ, θ, α, a, b, h, primes, U (ch03).

## Teaching lessons
- **Depth is tiered; coverage is exhaustive.** ch01: 15 A / 70 B / 21 C, 12 derivations, 496 cells. ch02: 16 A / 60 B /
  18 C, 15 derivations (122 steps), 405 cells. ch03: 15 A / 52 B / 12 C, 24 derivations (171 steps), 413 cells.
- **Convention callouts** (lapse rate, "Which p_o?", passive/active, first/second index, A:B, R vs A, **γ vs Γ vs
  circulation, translating vs rotating frame, signed b·n**): state the book's, the field's and the code's convention,
  give the size of the difference in numbers, say which one the code uses and why.
- **Book typos are taught, not hidden**: "book prints X; dimensions require Y; implemented Y" with a test that shows the
  printed form is wrong (ch03 (3.6), Ex. 3.2; ch02 Ex. 2.3/2.4/2.6).
- **Every key number by two independent routes** (formula vs tracked/measured/finite difference; brentq vs Lambert W;
  stencil vs `solve_ivp` path) — in the notebook and as a measured ◇ in the explainer.
- **★★★ sympy checks re-run the derivation's own construction** (ch01 D28, D19; ch02 D17; ch03 D22 builds the swept
  shell of a growing sphere). When sympy leaves `Subs(Derivative(...))`, substitute generic polynomials first.
- **Make an approximation measurable** and explain a measured deviation with a tiny example (ch03: the strain ellipse
  turns at ≈ γ/4, not γ/2; orders of smallness as fitted slopes).
- **Every prose number is computed; every "Dxx step N" pointer is checked by script after a split** (ch01 and again
  ch03 round 2: D08/D24 renumbering left three stale pointers). Symbolic twins are asserted against the numeric value.
- **Primers come before the first derivation that uses them**; `knowledge/primers.md` lists 110 (P01–P110); ch04 starts
  at P111.
- **The observer-dependence thread** (steady/unsteady, local/advective, relative/absolute vorticity) with climate hooks
  (thermal advection, planetary vorticity f = 2Ω sin φ) is ch03's model for threading one idea through a chapter.
- **Reorder the book when a derivation needs it, and say so**; **counter-examples teach definitions**; **theorems with
  hypotheses report, they do not assert**; **look at every animation's last frame**; **demo an idiom where a wrong
  statement would fail**; **every cited equation is written out** (rule above).
- Reusable derivation moves (ch01–ch03, 42 so far) are tabulated in `knowledge/concept_map.md` → "Reusable derivation
  moves".

## Global pitfalls confirmed in this project
**Machinery and environment**
- Extracted text garbles maths (`¼` = `=`, `ð…Þ` = parentheses, missing minus, ω → `u`, γ → g, θ → q, σ/τ → s, ν → n,
  ε → 3) — read equations from rendered pages.
- Windows: shell heredocs mangle backslashes — write code files with Write/Edit; very long shell commands fail; console
  needs UTF-8. **Python builder strings containing LaTeX must be raw** (`"\frac"` becomes a form feed; coverage check 9).
- Anaconda's `python3` kernelspec can shadow the venv's — notebooks execute on `fluidpy-venv`.
- JS: `UIEvent.view` is read-only. `tools/shot.py` does not click/drag; a reviewer runs a click pass.
- `test_viz_library_inlined_and_template_lints` fails whenever `assets/viz_lib.js` changed or an explainer is mid-build;
  run `tools/viz_inline.py --all` before the merge gate (ch03 O5: uncommitted library edits made 14 copies STALE).
- `tools/shot.py` `py:` parity expressions have **no builtins** — index dicts/tuples down to one float.
- Parallel viz-builders use **private scratch subfolders**.
- **matplotlib 3.11**: animations using `subplots_adjust`/`fig.text` disable `figure.constrained_layout.use` for the save.
- `tools/coverage_check.py` matches ledger reminder rows by exact primer name → false warnings (28 in ch02, 63 in ch03).
- `nb.derivation` prefixes "Eq." to its label — pass non-equation labels in the title.
- Book wording can hide in docstrings (`check_public.py` does not catch paraphrase): reviewers scan for it.
- A pager page is packed before the final `fit()` (library bug 2); tall display integrals clip at slice breaks (bug 3):
  keep derivation rows `\textstyle` on phones.

**Mathematics → code**
- **A test comparing two of our own functions is not evidence for the book's convention** (ch02 M1). Pin conventions to
  fields with known answers and **prove discrimination by patching in the wrong variant** (ch01 parcel ÷ρ_env; ch02 8/8;
  ch03 13/13 + 7 reviewer variants).
- **Preset names are physics**; **letters are physics too** — γ = 2S₁₂ vs Γ ≡ S₁₂ vs Γ circulation each pinned by a test.
- **Signed quantities stay signed**: RTT sliver b·n (|b·n| fails 5 tests), Leibniz lower term, wall stress, principal
  radii, lapse-rate conventions, Stokes orientation, clockwise spin ω₃ = −γ.
- **Change of variable back**: an extremum found in x = r²/σ² gives r = σ√x, not σx.
- **Rotating ≠ Galilean**: ω′ = ω − 2Ω; a translating frame changes the split, not the acceleration.
- **Units of numerical defaults**: a time-step default is never the space step (`ht` vs `h`).
- **Vectorised inputs**: trace over axes (0, 1) for (d, d, N) arrays; wrap scalars with `as_scalar_if_0d`.
- **Closed loops wrap** (index modulo n) or the last segment goes missing.
- Index-order conventions (traction first index vs divergence second; passive C vs active R; book vs Frobenius double
  dot; G vs the transposed integral gradient; R = 2A); grid layout `[k, j, i]`.
- sympy: symbols with assumptions are different symbols (use `core.index_notation.coordinates`); `.norm()` adds `Abs`;
  `Subs(Derivative)` does not simplify — substitute concrete generic functions.
- A default at a regime boundary is computed from the boundary function (ch01 −9.8e-3 vs Γ_a = −9.7607e-3).
- kmol vs mol; pint's [substance]; pint °C offset; `p_ref` vs `p0`; geopotential vs geometric altitude.
- `np.gradient` is first order at the edges; the project stencils are explicit second order. `bc="periodic"` only on
  periodic grids. Discrete ∇×∇φ and ∇·∇×u vanish to round-off.
- Cancellation: `-np.expm1(-x)` for 1 − e^{−x} near 0 (Gaussian vortex at r → 0).
- Unbounded ODE growth: cap with a `solve_ivp` event, return NaN beyond it. Statistical checks: seed, 5 standard errors.
- **Validation labels**: V4 only for real conservation laws; unit/frame invariance is V7; identities are V1; an
  encyclopedia *form* cross-check is V1, not V5 (ch03 relabelled Rankine, Lamb–Oseen, RTT, Leibniz); "converged" only
  when an order is asserted (ch03 `pathline`, `rtt_ellipse_2d` relabelled analytic); book-rounded constants only in the
  private JSON.

## Benchmark inventory
| value / table | source (citation) | stored in | used by |
|---|---|---|---|
| k_B, N_A, R (exact) | CODATA 2018/2022, NIST CUU | `reference/ch01/constants.json` | ch01 `thermo` |
| USSA-1976 constants, Tables 1–2 (0–50 km) | NASA-TM-X-74335; PDAS tables | `reference/ch01/ussa1976_*.{json,csv}` | ch01 C13–C63 |
| IAPWS σ(T), μ(T) | IAPWS R1-76(2014), R12-08 | `reference/ch01/iapws_sigma.csv`, `benchmarks.json` | ch01 C13, C15 |
| Jennings mean free path; Taylor blast K = 0.856; AMS Γ_d ≈ 9.8 K/km; EOS-80 check values | Tsalikis et al. 2024; Díaz arXiv:2009.05674; AMS Glossary; Fofonoff & Millard 1983 | `reference/ch01/benchmarks.json` | ch01 |
| Divergence theorem 8π/3; Levi-Civita identities; Cauchy index placement, Mohr; rotation matrices (scipy); Stokes statement | Wikipedia; scipy docs | `reference/ch02/benchmarks.json`, `SOURCES.md` | ch02 C02–C16 |
| **Gaussian (Lamb–Oseen) vortex peak α = r_max²/σ² = 1.256** | Canivete Cuissa & Steiner 2022, SWIRL I, arXiv:2210.05223, Sect. 2.4 Eq. (10) | `reference/ch03/benchmarks.json`, `SOURCES.md` | ch03 C14 (x* = 1.2564312086 inside the 3-digit interval; profiles identical to 1e-12) |
| Rankine and Lamb–Oseen profile forms; RTT and Leibniz sign forms (form cross-checks, V1) | Wikipedia "Rankine vortex", "Lamb–Oseen vortex", "Reynolds transport theorem", "Leibniz integral rule" | same | ch03 C14, C15 |
Book-printed values (private, git-ignored): `tests/book_values_ch01.json`, `…_ch02.json` (Examples 2.1–2.4, Ex. 2.1a),
`…_ch03.json` (Ex. 3.1 circles, §3.5 shear flow, (3.22)–(3.27), Gaussian radius, Ex. 3.2).

## Validation summary per chapter
Counts are concept-map rows carrying each label (a row can carry several); tests per evidence level in parentheses.
| chapter | analytic | symbolic | converged | conserved | benchmark | book-value | qualitative | unverified | verdict |
|---|---|---|---|---|---|---|---|---|---|
| ch01 | 39 (V1 58) | 16 (V2 19) | 5 (V3 8) | 5 (V4 7) | 18 (V5 13) | 1 (V6 5) | 1 (`wavelength_to_rgb`) | 0 | **PASS** — 107 tests; 5 explainers (115 parity rows); notebook PASS (496 cells) |
| ch02 | 43 (V1 41) | 21 (V2 13) | 11 (V3 15) | 2 (V4 1) | 9 (V5 6) | 8 (V6 4) | 0 | 0 | **PASS** — 92 tests (+ V7 12); 44 rows; 8/8 wrong variants caught; 5 explainers (127 rows); notebook PASS (405 cells) |
| ch03 | 38 (V1 58) | 23 (V2 32) | 14 (V3 14) | 5 (V4 4) | 2 (V5 1) | 8 (V6 4) | 0 | 0 | **PASS** — 120 tests (+ V7 7); 44 rows + 24 D rows; every CORE ≥ 2 levels; 24 derivations sympy-checked; 13/13 wrong variants caught (+7 by the reviewer); orders 2.00 (stencils, RTT, sector) and 1.00 (tracked); review 0 Must / 9 Should applied; 7 explainers PASS round 2 (145 rows, 2 312 views, every pager page); notebook PASS round 3 (413 cells, 129 s) |
(ch01 V7: 19 tests, smoke 2; ch02 V7: 12; ch03 V7: 7. Full suite at the ch03 merge gate: 327 tests.)

## Open across chapters
- **Library pass (orchestrator, when nothing reads `viz/`)**: fix the four open library bugs; promote the ranked JS
  candidates (ch02 + ch03 lists in `knowledge/viz_patterns.md`); then `tools/viz_inline.py --all` and
  `tools/shot.py --chapter ch01/ch02/ch03 --quick` + `templates/viz_example.html --quick`.
- **Skill pass**: the lesson candidates for `interactive-viz`, `math-to-python` §7, `verify-implementation`,
  `teaching-style` and `colab-notebook` listed in `knowledge/viz_patterns.md` (ch02 + ch03) are not yet in the skills.
- Machinery TODO: `shot.py` clicks/drags every view; `coverage_check` understands ledger reminder rows; `nb.derivation`
  skips "Eq." for non-numeric labels; `show_eqs` → `tools/nbkit.py`; `core.anim` saves inside a constrained-layout-off
  `rc_context`; `scripts/ch02_drawings.py` `sys.path` bootstrap.
- Non-blocking explainer follow-ups: ch01 E1–E5 (`knowledge/ch01.md` §9); ch02 E1 τ₁₂ label, E3/E5 portrait x-range
  (`knowledge/ch02.md` §9); ch03 E1 D04 label and ticks, E2 negative bar direction, E3 portrait legend and 21-page
  Explain, E6 km, E7 `\dot b`, narrow bar labels, 19-page D22 result (`knowledge/ch03.md` §9).
