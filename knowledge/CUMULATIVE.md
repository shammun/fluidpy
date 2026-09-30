# CUMULATIVE knowledge — fluidpy
(Rewritten by the knowledge-keeper after every chapter. Every agent reads this first. Last rewrite: after ch09,
2026-09-30, knowledge pass on commit aa34f17.)

## Project rules in force (added during ch01–ch03; sharpened since)
- **5–10 explainers per chapter, as many as the CORE ideas need** (`book.yaml → project.min/max_explainers_per_chapter`);
  ch03 has 7, ch04–ch09 9 each.
- **Every book equation is shown in full next to its number** — notebook prose, derivation steps, traps, recaps and
  explainer tour/Explain/Derivation/quiz/notes/status text **and `<meta>` strings** (ch04 E6, ch06 E6). Enforced by
  `tools/coverage_check.py` check 8 (markdown) and `tools/eq_refs.py` (explainers); reviewers still scan by eye.
- **`coverage_check` check 9**: no control characters in cells (a LaTeX backslash eaten by a non-raw string). Explainers
  have **no such lint yet** (ch09 shipped two, found by the reviewer — machinery TODO).
- **Rule 9 (public repo) in practice (ch07)**: a value we compute can round to an exercise's printed answer. Numbers
  that coincide with a book value are built from the function (never typed) and printed with **4 significant figures**;
  the book's value lives only in `tests/book_values_chNN.json` (check it is valid JSON — ch09's was not).
- **Pager and audit**: an item taller than a page is sliced at natural breaks; `tools/shot.py` pages through every page
  of every pager at every size. **The notebook is the reference for derivation step counts and titles** (ch08): explainer
  Derivation tabs are checked against `notebooks/build_chNN.py --dump` — **again after the lesson review** (ch09: D06 and
  D14 drifted after the viz review had passed).
- **Two values, the book's first** (ch09): when the book's number and the exact/computed one differ, both are shown
  everywhere, the book's first (Example 9.2: λ = −0.090 → x/L 0.158, then the exact Falkner–Skan fold −0.068148 → 0.126).

## Physics pipeline so far
```
ch01 Introduction — continuum ρ, p, T, u (D35, Kn) [C06]; τ = μ du/dy, ν = μ/ρ, FTCS [C12]; statics dp/dz = −ρg (D05),
     buoyancy (D37), USSA-1976 [C20]; first law → Gibbs (D10) → c² → p = ρRT (D11) → p/ρ^γ (D14) [C25–C45]; N², Γ_a,
     θ (D18–D20, D36) [C50–C55] (Kundu Γ ≡ dT/dz in code, meteorology −dT/dz shown); Π = null space (D28) [C64–C69]
ch02 Cartesian tensors — summation, C_ij (D01–D03), Cauchy f_i = τ_ji n_j (D05), τ' = CᵀτC (D06), invariants (D18),
     ε–δ (D09), S + A, R ↔ ω (D14, D15), principal axes (D17 ★★★), ∇ on grids, Gauss (D25), Stokes (D26)
ch03 Kinematics — DF/Dt (3.5) (D02); stream/path/streak (D03–D05); Galilean (D06); du = G·dx → S, ½ω (D07–D15);
     Rankine, Gaussian vortices (D17–D20); Leibniz (D21) → RTT (3.35) (D22 ★★★)
ch04 Conservation laws (the trunk) — mass (D01–D04, ψ); momentum, Bernoulli (4.19), Cauchy (4.24) (D05–D07);
     τ = −pδ + 2μS + λS_mmδ (D08, D09 ★★★) → NS (D11–D13); rotating frame ±2Ω×u′ (D14–D18); energy, ε ≥ 0 (D19–D23);
     Lamb, unsteady Bernoulli (D24–D26); Boussinesq (D27, D28); interfaces (D29); dimensionless NS (D30)
ch05 Vorticity — tubes, vortex pressure, Kelvin, baroclinic torque, Helmholtz (D01–D08); Dω/Dt (5.13) (D09 ★★★);
     Biot–Savart (D10, D11 ★★★); rotating baroclinic (5.30) (D14, D15 ★★★); stretching/tilting; (ζ + f)/h (D19, D20);
     point vortices, images, sheets γ = u₂ − u₁ (D21–D23)
ch06 Ideal flow — Laplace, superposition, walls (D01–D05); cylinder C_p = 1 − 4 sin²θ, D = 0, L = ρUΓ (D06–D13);
     w = φ + iψ (D14, D15); Blasius (6.60), Kutta–Zhukhovsky (D16–D18 ★★★); conformal maps (D19–D21 ★★★); Laplace on a
     grid (D22, D23); Stokes ψ, sphere, airship (D24–D27); added mass (D28–D31 ★★★)
ch07 Gravity waves — "plane wave in, boundary conditions → algebra, ω(k) out": linearised surface (D02–D04) →
     ω² = gk tanh kH (7.28) (D05, D06) → limits, orbits, energy (D07–D14); capillary (D15, D16); seiches; c_g, F = Ec_g
     (D19–D21); rays, Snell (D22–D24); jump, KdV (D25, D26); Stokes drift (D27); interfaces, √(g′H) (D28–D31 ★★★);
     internal waves ω = N cos θ, c ⟂ c_g (D32–D37)
ch08 Laminar flow — "symmetry kills u·∇u → linear ODE + wall conditions; else scale until one balance dominates; no
     imposed length or time → a similarity variable": Couette–Poiseuille, backflow, y* (D02–D05); pipe, τ₀ as a CV
     balance, f = 64/Re (D06, D07); circular Couette AR + B/R (D08, D09); lubrication ε²Re_L (D10 ★★★, D11) → Reynolds
     equation (D13 ★★★), slider (D14 ★★★), thin film (D15); Stokes' first problem η = y/√(νt), erfc (D16–D20); ansatz +
     exponent matching (D21–D23); Stokes layer δ_e = √(2ν/ω) (D24, D25); creeping flow, (E²)²ψ = 0, 6πμaU (D26–D31);
     far field ½Re_a r/a, Oseen (D32, D33)
ch09 Boundary layers (done) — "viscosity lives in a layer δ ~ L Re^{−1/2} whose pressure is imposed from outside;
     solve it by similarity, by an integral balance, or by marching; where U_e decelerates too much τ₀ → 0 and the
     layer separates; drag = the pressure the wake fails to recover"
  scaling (8.14 with ε = Re^{−1/2}) → (9.7)–(9.10), −(1/ρ)dp/dx = U_eU_e′ (9.11) (D01, D02) [C01]; δ*, θ, H; ρU²θ = ∫τ₀dx
     by a streamline-roofed CV (D03, D04) [C02]; ψ = Uδf(η) → f‴ + ½ff″ = 0 (9.27) (D05) [C03]; Töpfer, f″(0) = 0.332057,
     η₉₉ 4.910, θ = 2f″(0), C_f 0.664/√Re_x, C_D 1.328/√Re_L one side (D06) [C04]; Falkner–Skan (9.36), fold n = −0.090429
     (D07) [C05]; momentum integral (9.43) (D08) [C06] → Thwaites θ²U_e⁶/ν = 0.45∫U_e⁵dx (9.50), exact FS closure,
     cylinder 103.1° (D09–D11) [C07]; μu_yy(wall) = dp/dx, inflection, τ₀ = 0 (D12) [C08]; form drag
     sin φ_s(1 − (4/3)sin²φ_s − C_b) (D13) [C09]; regimes, St 0.2, Kármán street b/a = 0.28055 (D14 ★★★) [C10]; drag
     crisis, balls [C11]; free jet J const, sech², ṁ ∝ x^{1/3} (D15–D18, D16 ★★★) [C12]; wall jet Ψ const,
     4f‴ + ff″ + 2f′² = 0, (9.83) (D19–D21 ★★★) [C13]; teacup ρ(u_e² − u²)/R (D22) [C14] ← Ekman seed

next: ch10 CFD — the parabolic marcher (`march_boundary_layer`, von Mises + σ² mapping + BDF2, orders 1 and 2),
      Blasius/FS as reference solutions, continuation through a fold, `crank_nicolson_1d`, `core.laplace_solvers`,
      ch04 NS residuals + ch08 exact solutions as code-verification tests, `kdv_solve`.
```
The book ahead: Ch10 CFD → Ch11 instability → Ch12 turbulence → **Ch13 GFD (uses 4, 5, 7, 8, 9, 11, 12)** → Ch14
aerodynamics · Ch15 compressible · Ch16 biofluids.
Where earlier chapters feed in: ch01 N², θ → Ch. 11, 13; isentropic gas → Ch. 15; FTCS → Ch. 10. ch02 Gauss → every CV;
Stokes → Ch. 13 PV; invariants → Ch. 12. ch03 vortices → Ch. 13, 14. ch04 **NS residual/term tools → every exact
solution**; CV budgets → Ch. 12, 14, 15; ε → Ch. 12; rotating frame, Boussinesq, Ri, Ro → **Ch. 13**. ch05 vorticity
budget → Ch. 10, 12; sheet roll-up → Ch. 11; **(5.30), (ζ + f)/h → Ch. 13**; filaments → Ch. 14. ch06 `core.potential`,
`core.conformal`, `core.panels` → **Ch. 14** (and U_e(x) for ch09's Thwaites); `core.laplace_solvers` → **Ch. 10**, Ch. 13 PV
inversion. ch07 `core.waves` → **Ch. 13** (√(gH), Poincaré/Kelvin/Rossby, √(g′H)/f, WKB, adjustment), **Ch. 11** (KH, RT),
**Ch. 15** (jump ↔ shock), Ch. 10 (`kdv_solve`). ch08 `core.laminar` → **Ch. 11** (Couette/Poiseuille/Taylor–Couette base
states), **Ch. 12** (f = 64/Re, u_τ), **Ch. 13** (Stokes layer ≡ Ekman with ω → f, spin-up), Ch. 10; `core.lubrication` →
**Ch. 13** (hydrostatic thin layers, gravity currents); `core.creeping` → **Ch. 13** (settling), Ch. 16; `crank_nicolson_1d`
→ **Ch. 10**. **ch09** `core.boundary_layer` → **Ch. 10** (benchmarks, marching), **Ch. 11** (Blasius/FS base profiles,
inflection), **Ch. 12** (momentum integral, H, plate drag), **Ch. 14** (Thwaites on airfoils, stall, δ*); `core.bluff_body` →
**Ch. 11** (perturb–linearise–eigenvalues), Ch. 13 (island wakes), Ch. 14 (form drag, Magnus); `core.jets` → **Ch. 12**
(self-similar jets), Ch. 11 (Bickley jet), Ch. 13 (plumes, entrainment); teacup → **Ch. 13 Ekman pumping / spin-down**.

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
| constants (`K_B`, `R_AIR`, `GAMMA_AIR`, `G0`, `P_ATM`, `G_BOOK` = 9.81 …); perfect gas, first law, entropy, sound speed, isentropic, van der Waals | `thermo` | §1.8–1.9, (4.94) | analytic, symbolic, converged, conserved, benchmark |
| hydrostatics, USSA-1976, buoyancy; N², `classify_stability`, lapse-rate conventions, θ; FTCS; Π groups | `statics`, `stratification`, `diffusion`, `dimensional` | ch01 | analytic, symbolic, converged, benchmark |
| index notation, `tensors`, grids `[k, j, i]` + `operators`, `fields`, `integral_theorems` | ch02 modules | (2.1)–(2.36) | analytic, symbolic, converged, conserved, benchmark |
| `kinematics`, `coords`, `transport` (Leibniz, RTT); `vortices` | ch03/ch05 modules | (3.1)–(3.35), §5.4 | analytic, symbolic, converged, conserved, benchmark |
| `conservation` (any `ControlVolume`); `constitutive`; `navier_stokes` residuals + term splits; `curvilinear`; `streamfunction`; `rotating`; `bernoulli`; `interfaces`; `similarity` (`Scales`, groups, `sphere_drag_coefficient`) | ch04 modules | (4.5)–(4.119) | analytic, symbolic, converged, conserved, benchmark |
| `vorticity`, `biot_savart` | ch05 modules | (5.3)–(5.33) | analytic, symbolic, converged, conserved |
| `potential`, `conformal`, `laplace_solvers`, `panels` | ch06 modules | (6.1)–(6.109) | analytic, symbolic, converged, conserved, benchmark |
| `waves`: dispersion, `group_velocity` (signed, complex step), packets, `ray_trace`, two-layer and internal waves, `cosh_over_sinh` family | `waves` (ch07) | (7.1)–(7.159) | analytic, symbolic, converged, conserved, benchmark |
| `laminar`: channel/pipe/circular Couette (+ `_state`), Stokes' first/second problems, `similarity_variable(half=)`, `diffusion_thickness`, `vortex_sheet_diffusion`, `temporal_bl_wall_stress`, `stokes_layer` (δ_e, δ_book) | `laminar` (ch08) | (8.4)–(8.38) | analytic, symbolic, converged, conserved |
| `lubrication` (scales, profile, Reynolds equation, slider, Hele-Shaw, `thin_film_spread`, Huppert); `creeping` (Stokes sphere, drag laws, settling, Oseen) | ch08 modules | (8.13)–(8.53) | analytic, symbolic, converged, conserved, benchmark |
| `diffusion.crank_nicolson_1d` (θ-scheme, BE start-up, order 2) | `diffusion` (ch08) | (8.20), (8.33) | converged |
| **`boundary_layer`** (BL, 47 public): `boundary_layer_scales`, `to/from_bl_variables`, `bl_x_momentum_residual`, `bl_dpdy_scaled`, `bl_pressure_variation`, `outer_flow(kind)`; `thicknesses`, `displacement/momentum_thickness`, `delta_level`, `profile_shape`; `blasius_constants`, `blasius_profile`, `blasius_fields`, `blasius_delta99(printed=)`, `blasius_drag(sides=1)` …; `falkner_skan(m, branch=, method=)` (m = book n), `falkner_skan_state`, `_separation`, `_fields`, `_table`; `momentum_integral_residual`, `karman_pohlhausen`, `thwaites(closure=, lam_sep=)`, `thwaites_named`, `thwaites_l/H/L`, `thwaites_closure_table`, `thwaites_cylinder*`, `LAMBDA_SEP_BOOK` −0.09, `LAMBDA_SEP_FS` −0.068148 (lazy); `wall_curvature`, `profile_inflection`, `separation_point`, **`march_boundary_layer`** (von Mises, σ², BDF2); `transition_state`, `plate_drag_coefficient` | `boundary_layer` (ch09) | (9.1)–(9.52), Ex. 9.1–9.2 | analytic, symbolic, converged, conserved, benchmark |
| **`bluff_body`** (BB, 22): `cylinder_flow_regime`, `sphere_flow_regime` (rounded, qualitative), `cylinder_state`, `drag_crisis_state/pair`, `cylinder_cd_schematic`; `cp_ideal_cylinder`, `separated_cp`, `pressure_drag_from_cp`, `separated_pressure_drag(cp_base=-1.2)`; `karman_street_ratio/velocity/spectrum/spectrum_periodic/growth/growth_closed/positions`; `shedding_frequency` (cyclic St), `strouhal_of_re`; `ball_swing_deflection`, `magnus_sign` | `bluff_body` (ch09) | §9.7–9.9, (4.102) | analytic, symbolic, converged, qualitative (labelled) |
| **`jets`** (JET, 22): `free_jet_constants`, `free_jet_profile(_table)`, `free_jet`, `_centreline`, `_thickness`, `_mass_flux`, `_entrainment_velocity`, `_at_level`, `_halfwidth(printed=)`, `_reynolds`, `free_jet_ode_solve(coeff=)`, `jet_momentum_flux`; `wall_jet_profile`, `wall_jet`, `_mass_flux`, `wall_jet_ode_solve(coeff=4, printed=)`, `wall_jet_invariant`, `_integrals`, `_K1`, `_constants` | `jets` (ch09) | (9.53)–(9.85) | analytic, symbolic, converged, conserved, benchmark |
| **`similarity_reduce`** (moved from ch08 in ch09): `similarity_reduce_sympy(case, printed=)` (stokes1, sheet, vortex, bead, **blasius, falkner_skan, free_jet, wall_jet**), `similarity_ode_solve`, `similarity_collapse_error` (ch08 re-exports them) | `similarity_reduce` (ch08 → core) | (8.32), (9.19)–(9.82) | symbolic, analytic |
Chapter-only code (move to `core/` on second use): ch03 `pipe_profile`; ch04 `couette_heating`, `two_fluid_couette`,
`wake_drag_per_span`, `linear_wave_surface` (**overflows for kH ≳ 710; use `ch07.wave_fields`**); ch05
`rotating_cylinder_flow` (ω = 2Ω), `lamb_oseen_circulation`, `sheet_rollup` (→ Ch. 11), `column_relative_vorticity`
(→ Ch. 13), `ring_dynamics` (→ Ch. 14); ch06 `elliptic_cylinder_flow`, `lift_per_span` (→ Ch. 14), `sphere_motion`; ch07
`wave_fields`, `particle_path`/`stokes_drift` (→ Ch. 13), `hydraulic_jump`, `simple_wave_evolve` (→ Ch. 15), `kdv_solve`
(→ Ch. 10), `boussinesq_linear_sympy` (→ Ch. 13); ch08 sympy engines (`parallel_flow_sympy`, `lubrication_nondim_sympy`,
`reynolds_equation_sympy`, `stokes_second_sympy`, `low_re_scaling_sympy`, `oseen_*_sympy`); **ch09** sympy engines
(`bl_nondim_sympy`, `momentum_integral_sympy`, `thwaites_sympy`, `jet_momentum_sympy`, `wall_jet_invariant_sympy`,
`wall_jet_sympy`, `karman_street_sympy`, `cylinder_thwaites_sympy`, `derive_all`), `example_9_1/2`,
`secondary_flow_radial_force` (→ **Ch. 13** Ekman: move to core there), `book_slips`. **Default g**: ch04, ch05, ch07,
`core.waves` G_BOOK = 9.81; the rest of `core` G0 = 9.80665. **Default ρ**: ch06 1.2 vs 1000 (§6.9); ch07, ch08 1000;
ch09 BL helpers 1.2 (air), teacup 1000.

### Explainer engine (`assets/viz_lib.js`)
`Viz.app` (tabs Walkthrough / Explore / Explain / Derivation / Equations / Code / Check; fit-to-window with density
levels and pagers; transport, modes, presets, status, terms, inspector, notes, selftest parity), `Plot`, `Viz.field`,
`Viz.num` (RK4, odeint, brentq, erf/erfc **~1.2e-7 only — a defect; local `erfcHP` in ch08**, niceTicks), `Viz.work`,
`Viz.three`, KaTeX with fallback; `Viz.font`, `Viz.roundRect`, `Viz.text`, `Viz.card`, `Viz.fmtTime`, `Viz.tnum(keepTiny)`;
pager slicing, `--viz-stage-need`, `ev.viewId`/`ev.vizView`. **ch09 fix (aa34f17): a user-driven `set()` always rebuilds
the Explain HTML** (before, a hand-dragged transport parameter refreshed only `explain.live`, so Explain went stale when the
transport parameter was physical — x, n, Re; `shot.py` cannot see this). **Nothing else promoted after ch02–ch09** (the
site-publisher was reading `viz/` each time). Ranked lists in `knowledge/viz_patterns.md` (ch09: lint rule for control
characters, `views[i].hideOn`, `arrowPx`, formatters, gutter + log axes, `Viz.profiles`, hatch, `Viz.cx` before Ch. 11;
ch08: pager-label CSS, erfc precision; ch07: `Viz.waves`, per-mode hiding before Ch. 13). Reference: `templates/viz_example.html`.

## Explainer inventory
| chapter | slug | CORE idea | reusable stage pattern |
|---|---|---|---|
| ch01 | `continuum_averaging_volume` · `viscosity_momentum_diffusion` · `heat_work_paths` · `parcel_stability` · `buckingham_pi_machine` | C06 · C12 · C25–C45 · C50–C55 · C64–C69 | molecules vs box; gap + steady ghost; piston · p–v · T–s; column + two-convention badge; chips → matrix → groups |
| ch02 | `rotation_of_axes` · `cauchy_traction_principal_axes` · `strain_vs_rotation_split` · `gauss_flux_box` · `stokes_circulation_loop` | C02–C16 | clickable C matrix; element + Mohr; three squares; box + waterfall flux bars; unrolled u·t |
| ch03 | `flow_lines_unsteady` · `material_derivative_probe` · `galilean_frames_cylinder` · `fluid_element_deformation` · `spin_and_principal_axes` · `vortex_paddle_wheels` · `reynolds_transport_cv` | C01–C15 | streamline/trail/dye on one clock; local/advective bars with ◇; observer slider; G deforms a square; paddle wheels; moving CV + waterfall |
| ch04 | `control_volume_budgets` · `stream_function_spacing` · `newtonian_stress_lab` · `navier_stokes_term_balance` · `rotating_frame_coriolis` · `which_bernoulli` · `viscous_dissipation_heating` · `boussinesq_buoyancy` · `dynamic_similarity_models` | C01–C15 | five scenes, one budget; ψ spacing = speed; G → S → τ; term bars summing to 0; two observers; hypothesis table; kept vs dropped bars |
| ch05 | `vortex_tubes_cannot_end` · `vortex_pressure_funnel` · `kelvin_material_loop` · `baroclinic_torque` · `vorticity_stretching_tilting` · `biot_savart_filament` · `vorticity_equation_rotating` · `point_vortex_lab` · `vortex_sheet_rollup` | C01–C14 | 3-D tube + Gauss bars; ✓/✗ table; torque ◇; "build the sum"; conserved quantity flat in amber; click-to-place vortices; roll-up with KH ghost |
| ch06 | `superposition_sandbox` · `cylinder_circulation_lift` · `vortex_wall_images` · `complex_potential_corners` · `blasius_kutta_contour` · `conformal_joukowski` · `laplace_relaxation` · `axial_singularity_bodies` · `added_mass_sphere` | C03–C15 | element kits; ◇ ρUΓ bars; sensor trace; Laurent bars; two planes, same particles; stencil stepping + honest order |
| ch07 | `dispersion_relation` · `particle_orbits` · `capillary_gravity_waves` · `seiche_standing_waves` · `group_velocity_packets` · `wave_rays_refraction` · `hydraulic_jump` · `two_layer_modes` · `internal_wave_beams` | C03–C16 | tank + c(λ) with regime bands; linear vs exact paths; chord vs tangent + pond; in-page rays; budget bars + forbidden preset; two magnifications; square K plane |
| ch08 | `couette_poiseuille_backflow` · `lubrication_scaling` · `slider_bearing` · `viscous_gravity_current` · `stokes_first_problem` · `similarity_exponents` · `oscillating_plate` · `stokes_sphere_flow` · `stokes_drag_settling` | C02–C15 | line + parabola = sum, floor tangent at onset; log bars with a "dropped" band; printed ghost with a measured miss; log clock + raw/rescaled; exponent plane; rejected root drawn; three models + measured circle; traction sweep |
| ch09 | `bl_scaling_thicknesses` | C01, C02 (D01, D03, D04) | streamlines lifted by exactly δ*; shape chips at fixed δ₉₉; D04 closed with the reader's numbers |
| ch09 | `blasius_similarity_collapse` | C03, C04 (D05, D06) | raw ↔ rescaled collapse; shooting overlay (too small / too large / converged); τ₀(x) with "area = F_D" |
| ch09 | `falkner_skan_family` | C05, C08 (D07, D12) | one dial through a fold + dashed second branch; equal wall term bars μu_yy = dp/dx |
| ch09 | `thwaites_marching` | C06–C08 (D08–D11) | outer-flow picker → ∫U_e⁵ → λ crossing; θ² bar (integral vs inlet memory); inspector at a cursor; book/exact criteria, book first |
| ch09 | `cylinder_drag_crisis` | C09–C11 (D13) | log-Re dial moving cartoon + C_p + C_D + regime table; "set C_b by hand" check |
| ch09 | `karman_street_stability` | C10 (D14 ★★★) | configuration + growth curve + complex spectrum reaching the axis at the marginal value (the Ch. 11 template) |
| ch09 | `free_jet_similarity` | C12 (D15–D18) | wrong-exponent toggle breaks the invariant; level slider exposes the printed half-width |
| ch09 | `wall_jet_invariant` | C13 (D19–D21 ★★★) | two invariant bars (one falls, one flat to 4e-7); printed-ODE toggle; gauge preset |
| ch09 | `teacup_secondary_flow` | C14 (D22) | force bars that add with height (rose + teal = orange); river-bend mode; Ekman preview (Ch. 13 seed) |

## Notation
See `knowledge/notation.md` (sign traps first, then the register, then per-chapter conventions). Everywhere: SI and kelvin
inside functions; **z up** (except ch06 §6.8); **lapse rate Kundu Γ ≡ dT/dz in code, meteorology shown**; `p0` ≠ `p_ref`;
pressures absolute unless `_gauge` (ideal flow: from hydrostatic; ch07: **gauge**); signed τ_xy = μ ∂u/∂y. ch02: passive C;
traction on the first index; grid `[k, j, i]`. ch03: R = G − Gᵀ, spin ½ω; RTT outward n. ch04: Coriolis acceleration
+2Ω × u′ vs force −2Ω × u′; **2-D ψ: u = ∂ψ/∂y** (ch07: u = ∂ψ/∂z). ch05: γ = u_below − u_above; `ellipk(m)`, m = k². ch06:
**Γ counterclockwise in code, `Gamma_cw=`/`Gamma_ccw=`**; w = φ + iψ. ch07: **ω = angular frequency**, η = elevation,
**g′ = g(ρ₂ − ρ₁)/ρ₂**. ch08: **dp/dx book sign (favourable < 0), `G=` = −dp/dx**; walls y = 0 fixed, y = h moving; **η =
y/√(νt)** (figures y/(2√(νt))); ω vorticity and frequency; δ_e vs the book's 4√(ν/ω); §8.6 θ from the downstream axis; Re on
the diameter vs Re_a. **ch09 (changes):** x along the wall, y normal, u = ψ_y, v = −ψ_x; **dp/dx > 0 adverse**; **η = y/δ(x)**
with the flow's own δ (Blasius √(νx/U) = ch08's √(νt) with t = x/U; FS √(νx/U_e); free jet ∝ x^{2/3}; wall jet ∝ x^{3/4});
δ = similarity length (not δ₉₉); **θ = momentum thickness** (a length); **Falkner–Skan exponent: book n, code `m`**, β =
2m/(m + 1); **body angles φ from the FORWARD stagnation point**; **λ = Thwaites parameter (C07) and the street eigenvalue
(C10)**; **γ, σ in D14 = street coefficients** ½ − sech²(πb/a), sinh/cosh²; St with the **cyclic** frequency; (9.33) one
side (`sides=`); six Re (P200); C four times (Blasius, `C_fj` = 4√6/3, wall-jet dimensional C, constants); Ψ in m⁵/s³; J [N/m].
Multi-meaning symbols (register): through ch08 as before; **ch09 adds η, δ, θ, λ, γ, σ, n/m, a, C, L/l, J, Ψ, f, u₀, φ**.

## Teaching lessons
- **Depth tiered, coverage exhaustive** (A / B / C, derivations, cells): ch01 15/70/21, 12, 496 · ch02 16/60/18, 15, 405
  · ch03 15/52/12, 24, 413 · ch04 15/151/19, 30, 611 · ch05 14/69/11, 23, 496 · ch06 15/116/29, 31, 564 · ch07
  16/193/27, 37 (321 steps), 617 · ch08 15/113/15, 33 (301 steps), 501 · **ch09 14/118/16, 22 (224 steps), 472**.
- **Convention callouts** ("which p_o?", Coriolis term vs force, Γ_cw vs Γ_ccw, which ω/η/g′/p′/dp/dx/δ/Re/θ, **which λ,
  which η, which angle origin**): state each convention, give the size of the difference in numbers, say which the code uses —
  one numeric case carried through table → tiny example → code → explainer.
- **Book slips are taught by computing both versions**: a slips table up front (printed | corrected | where), printed forms
  as ghosts and as code options a test must fail (ch09: `printed=True`, `printed_9_7=True`, `coeff=1.0`); **check a suspected
  slip on the page before teaching it** (ch06 (6.82) was not one; **ch09 R5 "one side" is the book's own statement — a ⚠️
  trap, not a slip**; ch09 R17 was an unrecorded slip found by reading a sentence against (9.9)).
- **Two values, the book's first** (ch09 Example 9.2, cylinder 103.1° then 100.9°): both computed, same order everywhere.
- **Make an approximation measurable** (ch01, ch07, ch08 "measure every ~"; **ch09 Thwaites θ −7.6 … +4.4 %, τ₀ −56 % near
  the fold; Pohlhausen within 6 %; the drag-crisis model 0.65 of the measured collapse**); a result sensitive to its inlet is
  reported as a **bracket, not a number**.
- **Every key number by two independent routes** (ch06 lift; ch07 FD eigenproblem; ch08 three-way parity; **ch09 Blasius
  Töpfer = solve_bvp = shooting; Kármán closed form = Bloch 4×4 = periodic cell; wall jet IVP = implicit (9.83)**).
- **Every prose number is printed by a cell** — **including numbers in reading notes written during a review round** (ch09:
  a new note said 110°, the computed crossing is 132°); captions count what they claim.
- **Flag inserted derivation steps**; **a tool first used inside a derivation gets its primer first** (ch09 P218a integration
  by parts, moved before D06, D06 split 10 → 14 steps); `knowledge/primers.md` lists 221 (P01–P220 + P218a); ch10 starts at
  **P221**.
- **A sign-checking tiny example** for every transcribed equation whose sign matters (ch09 (9.8) solved for ∂p*/∂y*).
- **Builder rendering hygiene**: design TeX converted, not copied; regex tidies never run inside maths they created;
  headings keep their equations.
- **Climate hooks with numbers** (ch04 Ro; ch05 tornado; ch06 ψ inversion; ch07 swell, tsunami, √(g′H)/f; ch08 Stokes layer
  = Ekman with ω → f, spin-up, droplets, lava/ice; **ch09 teacup = Ekman pumping seed (1000 N/m³ at the floor), Kármán
  streets behind islands, entrainment in jets and plumes**).
- Reusable derivation moves (ch01–ch09, 157 rows counted; 19 from ch09) are tabulated in `knowledge/concept_map.md`.

## Global pitfalls confirmed in this project
**Machinery and environment**
- Extracted text garbles maths (`¼` = `=`, `ð…Þ` = parentheses, missing minus, ω → `u`, γ → g, θ → q, σ/τ → s, ν → n,
  ε → 3, ρ → r, Ω → U, ψ → j, φ → 4, η → h, ζ → z, ∂ → v, ∇ → V): read page images.
- **Windows: heredocs and `sed` mangle backslashes, quotes and `|`** — write code, explainers and knowledge files with
  Write/Edit or Python files; builder strings with LaTeX must be raw; **single-backslash TeX in JS strings turns `\f`, `\t`,
  `\r`, `\n` into control characters — four chapters now (ch04, ch06, ch09 ×2, one inside `String.raw` because the bytes were
  written into the file); explainers are not linted for it yet**.
- `MPLBACKEND=Agg` strips inline figures; **every chapter script takes `--no-show`**; Anaconda's kernelspec can shadow
  `fluidpy-venv`; JS `UIEvent.view` is read-only.
- `tools/shot.py` does not click/drag (**so it cannot see stale Explain text after a slider drag — ch09 R2-M1**); its `py:`
  expressions have no builtins and no `np`; set parity tolerances near the achieved agreement; run shot serially when memory
  is short and confirm a clip twice.
- **KaTeX-late pager clips are intermittent**: one line of slack; text-size integrals in live rows.
- `test_viz_library_inlined_and_template_lints` fails whenever `assets/viz_lib.js` changed or an explainer is mid-build;
  run `tools/viz_inline.py --all` before the merge gate. **Whole-tree `pytest -q` can stall near 96 %** when the browser
  machinery test follows the long ch09 script test — run `tests/test_machinery.py` separately (19 s).
- matplotlib 3.11: animations using `subplots_adjust`/`fig.text` disable constrained layout for the save.
- `coverage_check` matches ledger reminders by exact primer name; it does not see literal TeX outside `$…$`.
- Book wording can hide in docstrings and **exercise answers can hide in computed numbers** (ch09: Exercise 9.8 and 9.19
  values come out of `blasius_fields` and `falkner_skan(1)`; print ≤ 4 s.f.); reviewers scan against the private JSON.
- Library: pager packs before the final `fit()`; a tall "now" box is placed as a unit; `Viz.num.erf/erfc` ~1.2e-7; fixed
  transport range; no per-mode controls or per-layout view hiding (`views[i].hideOn` requested since ch03). **Phone legibility
  is the usual viz failure.**
- **Design and analysis documents contain errors; downstream agents compute, not copy** (ch05–ch08 several each; **ch09:
  "Table 9.1 differs < 5 %" false (G1), D06 and D14 lines wrong (G3, G4), "later separation ⇒ less drag at fixed C_b" false
  (G5), 60 vs 64 closure rows, a K with the wrong dimension**). **Explainer Derivation tabs drift when the notebook changes
  after the viz review** — re-diff against `--dump` after every lesson round.

**Mathematics → code**
- **A test comparing two of our own functions is not evidence for the book's convention.** Pin conventions to fields with
  known answers and **prove discrimination by planting the wrong variant** (ch02 8/8 … ch08 41/41, **ch09 25**).
- **A test on |x| cannot see the sign of x** (ch09: (9.8) coded negated passed an |δp| ∝ 1/Re test): every signed
  transcription gets a signed test on a manufactured field. **Tests must exercise every parameter at a non-zero value, with
  both signs** (ch08 M1); **derive a regime flag from the field and brute-force it** (ch08 slider); **every asymmetric solution
  gets a direction test** (ch08 Oseen wake).
- **A verifier's correction is a claim too** (ch08 far-field prefactor ½, not ⅛); **a mutant that hangs is not killed**, and a
  reviewer re-running a sample of mutants can find a survivor (ch09: 6 of 25 re-run, the (9.8) sign survived).
- **Folds and bifurcations are mathematics** (ch09 Falkner–Skan n = −0.090429): parametrise by the quantity that passes
  through (f″(0)); raise a clear error beyond; follow the second branch separately.
- **Boundary-value problems on infinite domains**: prefer a scaling symmetry (Töpfer) or `solve_bvp` with continuation to
  long shooting (sensitivity e^{η²/4}); brackets that scale with the solution; η_max per member (thick layers near a fold) and
  show insensitivity.
- **Complex-step derivatives need analytic continuations of every non-analytic piece** (ch07); keep |c| and signed dω/dk apart.
- **Printed formulas can assume a sign** ((7.138), (7.145)): code the gradient form, test all quadrants.
- **Code the corrected form; keep the printed one as a named option a test must fail**; keep the book's sign in the argument
  (`dpdx`) with the earlier chapter's as a keyword alias (`G=`).
- **Never type a constant a function can compute** (ch01 neutral lapse rate; **ch09 `LAMBDA_SEP_FS`**); keep the book's value
  as a separate named constant. **Defaults encode physics** (ch09 `cp_base=None` gave C_D,p 2.59).
- **Overflow and cancellation**: cosh/sinh ratios; x/sinh x via expm1; erfc, never 1 − erf; −expm1; small-α series; exact
  branches for R₂ = ∞ and R₁ = 0; NaN inside bodies; t ≤ 0 guards; **0/0 starts with a finite limit (Thwaites at a
  stagnation point θ² = 0.45ν/(6U_e′(0)))**.
- **Stiff and degenerate PDEs**: CN after a jump needs BE start-up steps; degenerate nonlinear diffusion needs implicit
  Newton + precursor + conservative faces; **parabolic marching with a √ψ wall singularity needs a mapping (ψ = ψ_max σ²) and
  BDF2, not CN** (ch09: order 0.98 → 2.24); measure each direction's order where its error dominates.
- **Inverse relations by a physics bracket** (brentq, residual asserted); implicit solutions inverted by brentq (ch09 wall jet).
- **Complex branches are physics** (Zhukhovsky outside root, `arctan2`); ½Re(AB*); keep the bounded root; **lattice sums
  converge only in ±n pairs** (ch09 street, order 1.99 asserted).
- **Singular corners and staircase boundaries set the global convergence order** (ch06 4/3; ch08 1.11; ch09 a pressure jump
  on a sampled trapezoid 0.98).
- **Finite-difference tolerances from a written budget**, never loosened by feel.
- **A telescoping sum is not a conservation test**; genuine V4 invariants so far: ch07 ∫η², KdV, ω along rays; ch08 thin-film
  volume, sheet jump, ∫ω dy, pipe CV, circular-Couette torque; **ch09 free-jet J (1e-10), wall-jet Ψ (1e-5), ρU²θ = ∫τ₀dx
  (1e-6), (9.43) residual**.
- **Sibling functions share defaults and argument meanings** or are called by keyword; convention-named keywords
  (`Gamma_cw=`, `ref=`, `G=`, `frame=`, `half=`, `p_scale=`, **`sides=`, `closure=`, `branch=`, `printed=`**).
- **Two quantities with one name**: g′, η, p′, δ, Re, **λ (Thwaites vs eigenvalue), θ (thickness vs angle), φ origin, n vs m**
  — the code names the reference in the keyword or the function name.
- **Benchmarks from the primary source at its printed precision**; an encyclopedia correlation is secondary (ch09 0.074).
- Signed and ordered quantities stay signed: b·n, wall stress, lapse rate, ∇ρ × ∇p, u₂ − u₁, Γ_cw vs Γ_ccw, dω/dk, dp/dx,
  **∂p*/∂y* in (9.8), wall curvature u_yy (R17)**.
- Library conventions: `ellipk(m)` with m = k²; `np.gradient` first order at edges; **numpy 2 has no `np.trapz`**; pint °C
  offset; kmol vs mol; geopotential vs geometric altitude; FFT evolution is periodic.
- **Validation labels**: V4 only for real conservation laws; identities V1; published closed forms reproduced are V1 form
  cross-checks; "converged" only with an asserted order; **"approximate" for few-percent agreement (Pohlhausen); a course
  handout is not V5**; a correlation bracket is V5 for the bracket only.

## Benchmark inventory
| value / table | source (citation) | stored in | used by |
|---|---|---|---|
| k_B, N_A, R (exact); e = 1.602176634e-19 C | CODATA 2018/2022, NIST CUU | `reference/ch01/`, `reference/ch08/` | ch01 `thermo`; ch08 Millikan |
| USSA-1976, sea-level c = 340.29 m/s; Sutherland β, S | NASA-TM-X-74335; PDAS tables | `reference/ch01/ussa1976_*` | ch01, ch04, ch08 ν_air |
| IAPWS σ(T), μ(T) | IAPWS R1-76(2014), R12-08 | `reference/ch01/` | ch01, ch04, ch07, ch08 |
| Jennings mean free path; Taylor blast 0.856; AMS Γ_d; EOS-80 | Tsalikis et al. 2024; arXiv:2009.05674; AMS Glossary; Fofonoff & Millard 1983 | `reference/ch01/benchmarks.json` | ch01 |
| Divergence theorem, Levi-Civita, Mohr, rotations, Stokes | Wikipedia; scipy docs | `reference/ch02/` | ch02 |
| Lamb–Oseen peak α = 1.256 | arXiv:2210.05223 | `reference/ch03/` | ch03 C14 |
| Sphere drag C_D(Re), 0.1–1e6 | F. A. Morrison (2016), Michigan Tech | `reference/ch04/` | ch04 C15; ch08 C15 bracket; **ch09 C11 crisis (C_D 0.09 near 4e5)** |
| WGS-84; orifice C_c 0.611; capillary length | Wikipedia | `reference/ch04/` | ch04, ch05 |
| Sphere added mass ½; Rayleigh collapse 0.91468; source panels | Lamb (1932) §92; Rayleigh 1917; Hess & Smith 1967 | `reference/ch06/` | ch06 |
| Fenton & McKee 1.7 %; Guo 0.75 %; capillary minimum; Stokes drift, H/λ = 0.1410633 | Coastal Eng. 1990, 2002; Wikipedia; arXiv:1507.02784 | `reference/ch07/` | ch07 |
| Slider K_opt 2.1889, W 0.0267; viscous gravity current η_N 1.411 | San Andrés (Texas A&M 2009/2012); Huppert, JFM 121 (1982) | `reference/ch08/` | ch08 C07, C08 |
| **Blasius f″(0) = 0.332057336215196** (ours 5.7e-15) | Wikipedia "Blasius boundary layer" | `reference/ch09/benchmarks.json` | ch09 C04 (V5) |
| **Falkner–Skan κ, separation β = −0.198837735, reversed branch β = −0.12, −0.02** (≤ 1.5e-12; 2.4e-10) | Belden et al., arXiv:1907.09912 | `reference/ch09/` | ch09 C05 (V5), R11 |
| **Hiemenz 1.232588, δ*/δ 0.6479** (2.8e-7, 7.3e-7) | Weidman & Turner (preprint); Wikipedia "Stagnation point flow" | `reference/ch09/` | ch09 C05 |
| **Bickley jet 0.4543, 0.2752, 3.3019, 0.5503** (≤ 1.4e-4) | Wikipedia "Bickley jet" (Bickley 1937) | `reference/ch09/` | ch09 C12 (V5) |
| **Kármán b/a ≈ 0.281** (0.28055, 3 digits) | Horváth et al. 2020 | `reference/ch09/` | ch09 C10 (corroboration) |
| Thwaites 0.45 + 6m and θ² forms (form check); turbulent plate 0.074Re^{−1/5} (**secondary**, +2.8 %) | Agrawal et al., arXiv:2310.16337; Wikipedia "Skin friction drag" | `reference/ch09/SOURCES.md` | ch09 C07, N69 |
| Forms (V1 cross-checks, not V5): Rankine, Lamb–Oseen, RTT; Bélanger, Taylor–Green; Kelvin ring, Burgers, Hill; cylinder, Kutta–Joukowski, Blasius force; KdV/cnoidal/solitary; Oseen, Stokes' law and second problem, Hagen–Poiseuille, Taylor–Couette | Wikipedia pages; K. T. McDonald | `reference/ch03/`–`reference/ch08/` | ch03–ch08 |
Book-printed values (private, git-ignored): `tests/book_values_ch01.json` … `…_ch09.json` (ch09: §9.3 coefficients, Table 9.1,
Examples 9.1–9.2, regime thresholds and angles, jet constants, exercise answers). **Open**: measured high-Re cylinder C_p
(ch06 N31), primary c_g,min (ch07 O7), primary pipe transition (ch08), **Glauert wall-jet constants, Lienhard cylinder C_D,
Roshko Strouhal, Howarth separation angles, a primary turbulent-plate coefficient (ch09)**.

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
| ch07 | 38 (V1 60) | 28 (V2 36) | 17 (V3 10) | 8 (V4 5) | 5 (V5 5) | 3 (V6 3) | 0 (2 labelled extensions) | 0 | **PASS** — 132 tests (+ V7 13); 41/41; 9 explainers (222 rows); 617 cells |
| ch08 | 29 (V1 46) | 20 (V2 31) | 6 (V3 5) | 7 (V4 4) | 7 (V5 6) | 0 rows (V6 3) | 0 | 0 | **PASS** loop 2 — 101 tests (+ V7 6); 41/41; 9 explainers (278 rows, 207 vs fluidpy); 501 cells |
| ch09 | 26 (V1 37) | 13 (V2 34) | 8 (V3 17) | 6 (V4 9) | 9 (V5 16) | 2 (V6 5) | 5 (labelled: model drag, regimes, crisis, teacup profile, thresholds) | 0 | **PASS** loop 1 (F1–F4 fixed) — 134 items at verify (+ V7 24), **136 after review** (2 Must: (9.8) sign, R5 relabelled; 9 Should applied); 22/22 D sympy-checked; 25 variants planted (1 hung; reviewer re-ran 6); 9 explainers PASS round 2 after the app.set fix (248 rows, 178 vs fluidpy); notebook PASS round 3 (472 cells, ~94 s) |
(V7 tests: ch01 19, ch02 12, ch03 7, ch04 3, ch06 2, ch07 13, ch08 6, ch09 24. Full suite at the ch09 merge gate: **1120 tests**.)

## Open across chapters
- **Library pass (orchestrator, when nothing reads `viz/`)**: ranked lists in `knowledge/viz_patterns.md` — ch09 (lint rule
  for control characters first; `views[i].hideOn`; `arrowPx`; formatters; plot gutter + log axes; `Viz.profiles` after Ch. 12;
  hatch; `Viz.cx` before Ch. 11), ch08 (pager-label CSS, **erfc precision**, gamma, golden, pager "now"-box slicing,
  derivation-result chain option), ch07 (**`Viz.waves` and per-mode hiding before Ch. 13**, Gauss–Legendre); then
  `tools/viz_inline.py --all` and `tools/shot.py --chapter ch01…ch09 --quick` + `templates/viz_example.html --quick`.
- **Skill pass**: lesson candidates for `interactive-viz`, `math-to-python` §7, `verify-implementation`, `teaching-style` and
  `colab-notebook` in `knowledge/viz_patterns.md` (ch02–ch09) are not yet in the skills (ch09 first: signed tests on
  manufactured fields, folds, book-first criteria, computed numbers in reading notes, primer before a derivation's first use,
  re-dump after lesson rounds, transport ≠ time).
- **Machinery TODO**: `viz_lint` control-character rule (**four chapters hit it**); a machinery test that drags a transport
  slider and asserts the Explain text changes; a Derivation-tab vs `--dump` diff **after the lesson review** in
  `embed_check`; solver iteration caps so mutants fail instead of hanging; JSON validity of `tests/book_values_*.json`;
  nbkit — `show_eqs` skips "printed" mentions, regex tidies only outside `$…$`, headings keep equations; `run_notebook`
  refuses `MPLBACKEND=Agg`; `shot.py` clicks/drags and exposes `np` to `py:` rows; pager re-pack after KaTeX; rule-9 scan
  against the private JSON; `eq_refs` exempts `<meta>` and `ref:` badges; ch08 `wide_end_backflow` alias.
- **ch09 explainer follow-ups** (`knowledge/ch09.md` §9, `viz_patterns.md` → ch09 open items): **D06 (10 vs 14 steps) and
  D14 (15 vs 16) must be rebuilt from `--dump`**; E1 `tau_0` backslash, δ*/θ colours, phone title; E2 collapse from data,
  phone quote; E3 phone view height; E5 mode-aware §0, label overlap; E6 `$+\Gamma$`, 31 vs 61, ±Γ convention; E7 label
  collision, exponent wrap; E8 legend overprint; E9 u(z) scale, labels. Earlier chapters' follow-ups in their `chNN.md` §9.
  Backups not built: ch04 `kinematic_free_surface`, ch05 `vortex_rings`, ch06 `flow_net_sources_vortices`, ch07
  `linearised_free_surface`, ch08 `rotating_cylinders_couette`, ch09 `ball_swing_magnus`.
- Verification notes: ch09 sphere Re_cr 5e5 (book) vs 3e5 (code) — decide (G9); marched separation near the fold quoted as a
  bracket; `thwaites_closure_table(fast=True)` not separately compared; `karman_street_spectrum` docstring σ → λ; turbulent
  0.074 needs a primary source before Ch. 12 prints it; ch08 O5 staircase order 1.11, O8 docstrings, O9 `inlet_backflow`
  naming; ch07 O5–O7; ch06 O4, O6; ch05 O2.
