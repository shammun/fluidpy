# CUMULATIVE knowledge — fluidpy
(Rewritten by the knowledge-keeper after every chapter. Every agent reads this first. Last rewrite: after ch02,
2026-09-16, commit b8c8837.)

## Physics pipeline so far
```
ch01 Introduction (done)
  molecules ──(box average, D35 noise N^-1/2, Kn = l/L)──► continuum fields ρ, p, T, u at a point      [C06]
     │                                                        │
     ├─ molecular transport: flux = −D·gradient                ├─ statics: dp/dz = −ρg (D05), buoyancy ρgV (D37)   [C20]
     │   Fick (1.1) · Fourier (1.2) · Newton τ = μ du/dy (1.3)  │     └─ static atmosphere from T(z): USSA-1976, H = RT/g
     │   ν = μ/ρ sets the clock h²/ν  → FTCS diffusion [C12]   │
     │                                                        ├─ thermodynamics: first law (1.10) [C25] → Gibbs (1.18, D10) [C35]
     │                                                        │     → c² = (∂p/∂ρ)_s (1.19) [C36] → p = ρRT (1.22, D11) [C40]
     │                                                        │     → isentropic p/ρ^γ (1.25, D14) [C45]
     │                                                        └─ stratification: parcel ζ″ + N²ζ = 0 (D18) [C50] → N² (1.29) [C51]
     │                                                              → Γ_a = −gαT/C_p (1.30, D19) [C54] → θ, N² = (g/θ)dθ/dz (D20, D36) [C55]
     │                                                              (Kundu Γ ≡ dT/dz in code; meteorology −dT/dz always shown)
     └─ dimensional analysis: homogeneity [C64] → dimensional matrix, rank (1.39) [C67] → Π groups = null space (1.37, D28) [C69]

ch02 Cartesian tensors (done) — the mathematical language of everything after
  summation convention (2.2) [C01] ──► direction cosines C_ij = e_i·e'_j, CCᵀ = I (D02, D03) [C02]
     └─► vectors: x' = Cᵀx (2.5, D01), "vector = transforms like this" (2.8) [C03]
  stress τ_ij (first index = face) [C04] ──► Cauchy f_i = τ_ji n_j (2.15, D05) [C05] ──► tensors: τ' = CᵀτC (2.12, D06) [C06]
     ├─► contraction, invariants I₁ I₂ I₃, book A:B = A_ij B_ji (D18) [C07]
     ├─► ε_ijk, epsilon–delta (2.19, D09), cross product (2.21) [C08]
     ├─► S + A split (D14); antisymmetric ↔ vector R_ij = −ε_ijk ω_k (D15); book R = G − Gᵀ ↔ ω = ∇×u [C12]
     └─► symmetric eigenproblem: real λ, orthogonal axes, bounds (D17 ★★★) [C13]
  fields on the project grid [k, j, i]: ∇φ (2.22) [C09] · ∇·u (2.23) [C10] · ∇×u (2.24–2.25, D12) [C11]
     └─► integral theorems: Gauss (2.30, D25) [C14] → integral definitions (2.31–2.33, D21, D22) [C15]
                            Stokes (2.34–2.35, D26), circulation per area [C16]

next: ch03 kinematics (S_ij, R_ij, vorticity, Reynolds transport — built on ch02's tensors, grids, operators, theorems)
      → ch04 conservation laws / Navier–Stokes (the trunk)
```
The book's structure ahead: Ch3 kinematics → **Ch4 conservation laws / Navier–Stokes (the trunk)** → Ch5 vorticity · Ch6
ideal flow · Ch7 gravity waves · Ch8 laminar flow → Ch9 boundary layers → Ch10 CFD → Ch11 instability → Ch12 turbulence
→ **Ch13 GFD (uses 4, 5, 7, 8, 11, 12)** → Ch14 aerodynamics · Ch15 compressible · Ch16 biofluids.
Where ch01 feeds in: N² and θ → Ch. 4 Boussinesq, Ch. 7 §7.8, Ch. 11 §11.7, Ch. 13; Newtonian stress and ν → Ch. 4 §4.5,
Ch. 8, 9; hydrostatic base state → Ch. 4, 7, 13; Π groups → Ch. 4 §4.11 and every scaling; isentropic gas → Ch. 15;
FTCS diffusion → Ch. 5, 8, 10.
Where ch02 feeds in: `strain_rate_tensor`, `rotation_tensor`, `vector_gradient`, grids, operators and `VectorField`
are **the kinematics spine of Ch. 3** (§3.4 deformation, vorticity, principal strain rates); Cauchy + principal axes +
τ symmetry → Ch. 4 §4.3–4.5; tensor divergence → Cauchy's equation; Gauss → every control-volume argument (Ch. 3
Reynolds transport, Ch. 4); Stokes/circulation → Ch. 5 Kelvin, Ch. 6 lift, Ch. 13 potential vorticity; invariants →
Ch. 12; ε identities → every vector identity (Coriolis 2Ω × u in Ch. 4/13).

## Available primitives
### Machinery (`fluidpy/core/`, ready before chapter 1)
| module | what | used for |
|---|---|---|
| `project` | repo root, `book.yaml` config, Pages/raw/Colab URLs | every tool and notebook |
| `style` | `setup_notebook()` (matplotlib style, plotly renderer, FAST), `COLORS`, `savefig` | every notebook figure |
| `embed` | `show_viz(chapter, slug)` — explainer in Jupyter (srcdoc) / Colab (Pages src) / page | explainer cells |
| `anim` | `animate`, `show_animation(player="video"/"frames")` | animations (matplotlib 3.11: disable constrained layout for `subplots_adjust` animations, see pitfalls) |
| `interact` | `slider_figure`, `animate_figure` (plotly, work on the page), `live` (ipywidgets) | Python interactives |
| `units` | pint registry `ureg`/`Q_`, `dimensional_check`, `check_dimensions`, `DIM`, `celsius_to_kelvin`, `kelvin_to_celsius`, `nondimensional` | V2 tests, temperatures |
| `refdata` | reference-data resolution + registry | benchmarks |
| `_util` (private) | `as_scalar_if_0d` (scalar in → float out, for parity rows), `require_positive`, `require_nonnegative` | every physics function |
| `tools/convergence.observed_order(h, err)` | log–log slope | every V3 test and convergence figure |

### Physics primitives from ch01 (re-exported as `ch01.<name>`)
| function | module | book § / Eq. | label |
|---|---|---|---|
| `K_B, N_A, N_A_KMOL, R_U, M_W_AIR, R_AIR, GAMMA_AIR, CP_AIR, CV_AIR, G0, P_ATM, P_REF` | `thermo` | §1.9 (CODATA, USSA-1976) | benchmark |
| `perfect_gas_density/pressure/state`, `gas_constant`, `molecular_gas_pressure`, `molecule_mass` | `thermo` | (1.21), (1.22) | benchmark, analytic |
| `specific_volume`, `enthalpy`, `perfect_gas_internal_energy/enthalpy`, `specific_heat_cp/cv`, `partial_derivative` | `thermo` | (1.12)–(1.15) | analytic, converged |
| `process_path`, `path_heat_work_totals`, `process_heat_work`, `join_paths` | `thermo` | (1.10), (1.11) | analytic, converged, conserved |
| `perfect_gas_entropy_change(_p)`, `entropy_change_reversible`, `irreversible_process`, `free_expansion`, `stirred_isochoric_process` | `thermo` | (1.16)–(1.18) | symbolic, analytic, conserved |
| `sound_speed_from_eos`, `perfect_gas_sound_speed`, `tait_pressure`; `thermal_expansion_coefficient`; `cv_from_cp`, `gamma_from_cp`, … | `thermo` | (1.19)–(1.24), (1.27), (1.28) | analytic, benchmark |
| `isentropic_pressure`, `isentropic_ratios`; `van_der_waals_*` | `thermo` | (1.25), (1.26) | symbolic, converged, analytic |
| `gauge_pressure`, `absolute_pressure`, `hydrostatic_pressure_uniform`, `layered_pressure`, `integrate_hydrostatic(z, rho_fn, p0)` | `statics` | (1.8), (1.9) | analytic, converged |
| `atmosphere_from_temperature`, `linear_lapse_pressure`, `standard_atmosphere(z, geometric=False)`, `isothermal_pressure/density`, `scale_height` | `statics` | §1.10, USSA-1976 | benchmark, symbolic |
| `buoyancy_force`, `net_pressure_force_on_box` | `statics` | D37 | analytic, conserved |
| `brunt_vaisala_sq(_from_lapse/_from_theta)`, `classify_stability`, `stability_timescale`, `parcel_displacement`, `parcel_ode*` | `stratification` | (1.29), D18 | analytic, symbolic, conserved |
| `adiabatic_lapse_rate` (negative), `lapse_rate`, `parcel_temperature`, `lapse_rate_convention`, `lapse_rate_stability` | `stratification` | (1.30), (1.32) | symbolic, analytic, benchmark |
| `potential_temperature`, `temperature_from_potential`, `potential_temperature_gradient`, `potential_density`, `isentropic_density_gradient`, `ocean_potential_density_gradient` | `stratification` | (1.31)–(1.35) | analytic, symbolic, conserved |
| `ftcs_diffusion_1d`, `stable_time_step`, `couette_startup_profile`, `gaussian_spreading`, `derivative_2nd_order` | `diffusion` | §1.5 (model PDE ours) | converged, conserved, analytic |
| `dimension_vector`, `dimensional_matrix`, `minor_determinant`, `rank_by_minors`, `solve_exponents`, `pi_groups`, `groups_independent`, `group_value`, `rescale_units`, `group_latex`; presets `PIPE, PENDULUM, SPHERE_DRAG, SCALE_HEIGHT, BLAST, RAYLEIGH, PYTHAGORAS` | `dimensional` | (1.36)–(1.40) | symbolic, analytic |

### Mathematical and field primitives from ch02 (re-exported as `ch02.<name>`)
| function | module | book § / Eq. | label |
|---|---|---|---|
| `expand_indices(term)`, `expand_indices_str`, `classify_indices`, `tensor_order`, `rename_dummy`, `comma_to_partial`; `coordinates(dim)` (real symbols — use these), `symbol_vector/matrix/tensor`, `field_functions` | `index_notation` | (2.2), (2.9), (2.17), (2.19), (2.21), (2.36) | analytic, symbolic |
| `unit_vectors`, `vector_from_components`, `dot`, `inner`, `outer`, `tensor_product`, `contract(A, B, pattern)`, `dot_tensor_vector`, `double_dot(A, B, convention=)`, `trace`, `invariants`, `characteristic_polynomial` | `tensors` | (2.1)–(2.3), (2.9), (2.14), Ex. 2.9 | analytic, symbolic |
| `direction_cosines`, `rotation_matrix_2d(theta)` (passive C), `rotation_matrix_3d(axis, angle)`, `rotation_angle`, `random_rotation(rng)`, `is_orthogonal`, `is_proper_rotation`, `orthogonality_residual` | `tensors` | §2.2, Ex. 2.8 | analytic, symbolic, benchmark |
| `transform_vector(u, C)` = Cᵀu, `inverse_transform_vector` = Cu', `transform_tensor(T, C)` (any order), `transforms_as_vector/tensor` (residual) | `tensors` | (2.5)–(2.8), (2.12), (2.13) | analytic, symbolic |
| `STRESS_CUBE_FACES`, `cube_face_tractions`, `stress_component_meaning`, `tetrahedron_face_areas`, `traction(tau, n)` (first index), `traction_components` (moved from ch01), `normal_shear_stress`, `mohr_circle_2d` | `tensors` | §2.4, (2.15) | analytic, symbolic, benchmark |
| `kronecker_delta`, `levi_civita`, `permutation_sign`, `epsilon_delta_residual`, `triple_product`, `is_isotropic`, `cross`, `cross_einsum`, `angle_between` | `tensors` | (2.16)–(2.21) | analytic, symbolic, benchmark |
| `is_symmetric`, `is_antisymmetric`, `independent_components`, `symmetric_part`, `antisymmetric_part`, `strain_rate_tensor(G)`, **`rotation_tensor(G)` = G − Gᵀ**, `antisymmetric_from_vector`, `vector_from_antisymmetric`, `symmetric_double_contraction` | `tensors` | (2.26)–(2.29) | analytic, symbolic |
| `principal_axes(tau)` (λ ascending, det B = +1, raises if non-symmetric), `diagonalize`, `normal_stress_bounds` | `tensors` | §2.11 | analytic, symbolic, benchmark |
| `grid(bounds, n, periodic=False)` → `Grid3D`, `grid2d` → `Grid2D`, `axis_of_direction`, `spacing`, `evaluate_field` | `grids` | §2.9 (project layout) | analytic |
| `partial(f, direction, h, bc)`, `second_partial`, `gradient`, `directional_derivative`, `divergence`, `vector_gradient` (G[i, j] = ∂u_i/∂x_j), `tensor_divergence(index=1)`, `curl`, `curl_components`, `laplacian`, `is_solenoidal`, `is_irrotational` | `operators` | (2.22)–(2.25) | analytic, converged (orders 1.97–2.03) |
| `VectorField(exprs, coords, params, singular_at)` (callable; `.div_fn`, `.curl_fn`, `.grad_fn`, sympy exprs), `ScalarField` | `fields` | Ex. 2.3 | analytic, symbolic |
| `volume_integral_box`, `divergence_theorem_box(rule=)`, `gauss_gradient_box` ([i, j] = ∂Q_j/∂x_i), `flux_through_box_faces`, `divergence_theorem_sphere`, `flux_through_sphere`, `flux_through_faces`, `divergence_theorem_rect2d`, `divergence_theorem_tiled`; `midpoint_nodes`, `simpson_nodes`, `gauss_legendre_nodes` | `integral_theorems` | (2.30) | analytic, converged, conserved, benchmark |
| `integral_gradient`, `integral_divergence`, `integral_curl`, `integral_divergence_2d` | `integral_theorems` | (2.31)–(2.33) | analytic, symbolic, converged |
| `boundary_tangent(n_c, n)`, `planar_loop`, `rectangle_loop`, `planar_disc`, `planar_rectangle`, `circulation`, `curl_flux`, `stokes_theorem_check` → `StokesCheck(lhs, rhs, hypothesis_ok, note)`, `integral_curl_component` | `integral_theorems` | (2.34), (2.35) | analytic, symbolic, converged |
Chapter-only code stays in the chapter modules: ch01 (continuum sampling, Kn, surface tension, μ(T), seawater EOS,
Examples 1.2–1.5); ch02 (`polar_components`, `example_2_2/2_3/2_4`, `traction_2d`, `stress_vs_angle`,
`principal_angle_2d`, test fields `radial_field`, `solid_body_rotation_field`, `shear_field`,
`irrotational_vortex_field`, `point_source_field`, `potential_field`, `smooth_test_field`, `periodic_*`,
`exact_div_curl`; linear-flow kinematics `velocity_gradient_preset`, `linear_flow_map`, `deform_square`,
`material_line_angle`; `operator_convergence`, `integral_definition_convergence`). **Move a function to `core/` when a
second chapter calls it** — Ch. 3 will call the linear-flow kinematics and the test fields.

### Explainer engine (`assets/viz_lib.js`)
`Viz.app` (tabs: Walkthrough / Explore / Explain / Derivation / Equations / Code / Check; fit-to-window with density
levels and pagers; transport, modes, presets, status, terms, inspector, notes, selftest parity), `Plot` (world
coordinates), `Viz.field` (grid, contours, heatmap, streamlines, quiver, particles), `Viz.num` (RK4, odeint, brentq,
erf, niceTicks), `Viz.work` (Explain builders), `Viz.three` (3-D), KaTeX with fallback. Added after ch01: `Viz.font`,
`Viz.roundRect`, `Viz.text` (halo/bg), `Viz.card`, `Viz.fmtTime`, `Viz.tnum` fix + `keepTiny`. Pointer events:
`ev.viewId`, `ev.vizView` (never assign `ev.view`). **Nothing promoted after ch02** (the publisher was reading `viz/`);
the ch02 candidates (`Viz.fixed`, `Viz.num.principal2d`, `Viz.matrix`, `P.axes({xlog})`, `Viz.hatch`, `Viz.num.tiles`,
diverging `bwo` map, `Viz.num.expm2`, engine row-hide flag, adaptive status) are listed in `knowledge/viz_patterns.md`.
Reference: `templates/viz_example.html`.

## Explainer inventory
| chapter | slug | CORE idea | reusable stage pattern |
|---|---|---|---|
| ch01 | `continuum_averaging_volume` | C06 continuum (+Kn, D35) | molecules view · measured value vs log box size with expected band · regime strip; resample + click-a-sample inspector; hidden-view number in a visible title |
| ch01 | `viscosity_momentum_diffusion` | C12 Newton viscosity, ν | gap with shearing tracers · profile with steady ghost and fading trail · wall stress vs time; modes = different diffusivity; draggable plate; end-of-run card |
| ch01 | `heat_work_paths` | C25, C35, C45 | piston · p–v · T–s linked on a path parameter; shaded/hatched areas; term bars (path vs state functions); reversible/irreversible modes |
| ch01 | `parcel_stability` | C50, C51, C54, C55 | column coloured by θ/ρ · profile with parcel adiabat · ζ(t) with linear ghost, one clock; two-convention badge with exact-text parity; atmosphere/ocean/lab modes; y-title strip |
| ch01 | `buckingham_pi_machine` | C64, C67, C69 | variable chips → matrix with highlighted minor and determinant → groups; exact rational JS; unit-system switch; problem presets |
| ch02 | `rotation_of_axes` | C02, C03, C06 (D01, D02, D06) | plane with fixed arrow and turning axes · clickable direction-cosine matrix with angle inspector · τ′(θ) curves; vector/tensor/polar modes; passive/active table; "is it a vector?" counter-example table; mirror toggle (det −1) |
| ch02 | `cauchy_traction_principal_axes` | C04, C05, C13 (D05, D17 ★★★) | stressed element with a cut, n, f and blue/rose split · σn/τs vs φ with λ bounds · Mohr circle, one angle; σn term bars; stress/strain modes; non-symmetric contrast toggle; shrinking tetrahedron in D05 |
| ch02 | `strain_vs_rotation_split` | C12 (D14, D15) | one clock, three squares (G, S, A) + matrix view + material-line angles; `expm2` exact trajectories; Frobenius term bars with S:A ≡ 0; kinematic presets |
| ch02 | `gauss_flux_box` | C14, C15 (D25, D22, D21) | field heatmap (white-zero diverging) with a draggable box and signed face arrows · waterfall flux bars (∮ vs ∬, interior faces) · limit view (1/A)∮ vs log h; tiling with shared-face inspector; shrink transport; singular-source preset |
| ch02 | `stokes_circulation_loop` | C11, C16 (D26, D12) | curl heat image + streamlines + tracers + paddle wheel with a draggable loop · unrolled u·t vs arc length (area = Γ) · side/term bars; circle/rectangle modes; flip-n; vortex-core preset where Stokes fails on purpose; tiles in D26 |

## Notation
See `knowledge/notation.md` (symbol register and sign traps). Conventions that matter everywhere: SI and kelvin inside
functions; **z up** (ch01); **lapse rate Kundu Γ ≡ dT/dz in code, meteorology −dT/dz always shown alongside**; `p0`
(pressure at z = 0) ≠ `p_ref` (θ reference, 1000 hPa); q heat added to / w work done on the system; gas constants per
kmol; pressures absolute unless `_gauge`; signed τ_xy = μ ∂u/∂y.
From ch02 (project-wide):
- **passive direction cosines** C_ij = e_i·e'_j, x' = Cᵀx, τ' = CᵀτC; scipy/Wikipedia R is active (= C applied as Cᵀ);
- **traction contracts the first index** f_i = τ_ji n_j; **tensor divergence the second** ∂τ_ij/∂x_j;
- **book A:B = A_ij B_ji** vs Frobenius — always pass `convention=`;
- **velocity gradient G[i, j] = ∂u_i/∂x_j**; `integral_gradient`/`gauss_gradient_box` return the transpose;
- **book R = G − Gᵀ with vector ω = ∇×u; A = ½R with vector ½∇×u — never call the vector of A "ω"**;
- **Stokes: n_c points INTO A, t = n_c × n counterclockwise about n**; closed surfaces use the outward n;
- **grid arrays `[k, j, i]` = (z, y, x), x on the last axis, components on axis 0; `Grid.h`, directions and
  components ordered (x, y, z)**; nodes; periodic grids drop the end node;
- λ^k, b^k: k is a label.
Symbols with several meanings so far: α, θ, T, h, R, σ/S, c, D, n, k, q, γ, λ, τ, ε, H (ch01); plus Γ (lapse rate
ch01 / shear rate ch02 / circulation ch05 — ch02 writes Γ_circ), ω (vorticity from ch02 / angular frequency ch07), φ,
A, C, G, K, m, s, t, a, b (ch02).

## Teaching lessons
- **Depth is tiered; coverage is exhaustive** (adopted 2026-09-12). Every inventory row gets an ID and a depth; only
  12–18 load-bearing ideas get the full treatment (A: picture → question → derivation → worked number → code →
  figure); B stated inside the nearest A block, C named with a pointer; derivations only for A items or results the
  book never writes out. ch01: 15 A / 70 B / 21 C, 12 derivations, 496 cells, 2 review rounds. ch02: 16 A / 60 B / 18 C,
  15 derivations (122 steps), 405 cells, 2 review rounds.
- **Convention callouts** ("Which p_o?", the lapse-rate table, passive/active, first/second index, A:B, R vs A): state
  the book's, the field's and the code's convention, give the size of the difference in numbers, say which one the
  code uses and why.
- **★★★ sympy checks re-run the derivation's own construction** (ch01 D28, D19; ch02 D17 builds C from unit
  eigenvectors and checks CᵀτC is diagonal), not just the final identity.
- **Make an approximation measurable** (ch01 C50 linear-vs-full table; ch02 C15 integral definitions converging at
  order 2, C16 staircase → circle).
- **Every prose number is computed, not typed**; step pointers ("step 12 of D28") are checked by script after any
  renumbering. **Symbolic twins are asserted against the numeric result in the same cell** (ch02: a symbol-assumption
  mismatch printed 0 under a comment claiming 3a).
- **Primers come before the first derivation that uses them**; `knowledge/primers.md` lists 86 (P01–P86). P numbers
  continue across chapters.
- **Reorder the book when a derivation needs it, and say so** (ch02: Cauchy (2.15) before the tensor rule (2.12)).
- **Counter-examples teach definitions**: a residual table where non-vectors fail by order one (ch02 C03, E1).
- **Theorems with hypotheses report, they do not assert** (ch02 `stokes_theorem_check.hypothesis_ok`; E5 breaks
  Stokes on purpose).
- **Look at every animation's last frame** (clipped readouts appear when numbers grow) and **demo a Python idiom on an
  object where a wrong statement would fail** (ch02 `np.transpose` primer).
- Reusable derivation moves (ch01 and ch02) are tabulated in `knowledge/concept_map.md` → "Reusable derivation moves".

## Global pitfalls confirmed in this project
**Machinery and environment**
- Extracted text garbles maths (`¼` = `=`, `ð…Þ` = parentheses, missing minus, ω → `u`, γ → g, θ → q, σ/τ → s, ν → n,
  ε → 3); matrices and minus signs vanish — read equations from rendered pages.
- Windows: shell heredocs mangle backslashes — write code files with Write/Edit; very long shell commands fail
  (`ENAMETOOLONG`) — use the Write tool for large files; console needs UTF-8.
- Anaconda's `python3` kernelspec can shadow the venv's — notebooks execute on `fluidpy-venv`.
- JS: `UIEvent.view` is read-only (assigning it broke every pointer handler). `tools/shot.py` still does not
  click/drag; a reviewer runs a click pass.
- `tests/test_machinery.py::test_viz_library_inlined_and_template_lints` fails whenever `assets/viz_lib.js` changed or
  an explainer is mid-build; run `tools/viz_inline.py --all` before the merge gate.
- `tools/shot.py` `py:` parity expressions are evaluated **without builtins** — index dicts/tuples down to one float.
- Parallel viz-builders must use **private scratch subfolders** (two builders clobbered one `splice.py` in ch02).
- **matplotlib 3.11**: `print_figure` re-installs constrained layout after the first saved frame; animations that use
  `subplots_adjust`/`fig.text` set `plt.rcParams["figure.constrained_layout.use"] = False` for the save and restore
  it afterwards.
- `tools/coverage_check.py` matches design-ledger reminder rows by exact primer name → false warnings (28 in ch02).
- `nb.derivation` prefixes "Eq." to its label — pass non-equation labels ("Exercise 2.8", "§2.10") in the title.
- Book wording can hide in docstrings (`tools/check_public.py` does not catch paraphrase): reviewers scan for it.

**Mathematics → code**
- **A test comparing two of our own functions (f(B) == g(B)) is not evidence for the book's convention** (ch02 M1:
  `rotation_tensor` was ½ of the book's R and passed). Pin conventions to fields with known answers
  (∇×(b × x) = 2b), and **prove discrimination by patching in the wrong variant** (ch01 parcel ÷ρ_env; ch02 8/8
  convention variants caught).
- **Preset names are physics** (ch02 M2: "uniaxial extension" was divergence-free): test a preset by its defining
  property (trace, spin, area growth).
- **Index-order conventions** (ch02): traction first index vs tensor divergence second; passive C vs active R;
  book vs Frobenius double dot; G[i, j] = ∂u_i/∂x_j vs the transposed integral gradient; R = 2A.
- **Grid layout** `[k, j, i]` with (x, y, z) component/spacing order — a deviation silently transposes curls.
- **sympy symbols with assumptions are different symbols** (`Symbol('x1', real=True)` ≠ `Symbol('x1')`); build
  expressions from `core.index_notation.coordinates`. sympy 1.14 `.norm()` adds `Abs` — normalise with `sqrt(b·b)`.
- A default at a regime boundary must be computed from the boundary function (ch01 −9.8e-3 vs Γ_a = −9.7607e-3).
- Sign conventions: lapse rate (Kundu vs meteorology), first-law q/w, signed wall stress, signed principal radii,
  Stokes orientation (n_c into A).
- kmol vs mol (factor 1000); pint's [substance] dimension dropped consistently; pint °C is an offset unit.
- Reference pressure `p_ref` vs surface pressure `p0`; geopotential vs geometric altitude; old EOS temperature scales.
- `np.gradient` is first order at the edges unless `edge_order=2`; the project stencils are explicit second-order
  one-sided (patching in `np.gradient` failed 4 of 5 order tests). `bc="periodic"` only on periodic fields and grids.
- Discrete ∇×∇φ and ∇·(∇×u) vanish to round-off on the project grid (stencils along different axes commute).
- `np.transpose(a, axes)` puts old axis `axes[p]` in position p.
- Unbounded ODE growth (N² < 0): cap with a `solve_ivp` event, return NaN beyond it.
- Neutral tolerance bands differ between sibling functions; sibling functions with different argument orders — call
  by keyword.
- Statistical checks: seed, 5 standard errors.
- Validation labels: V4 only for real conservation laws or state-function invariants (ch02: zero net outflux of a
  solenoidal field); unit/convention/rotation invariance is V7; identities are V1; a scipy cross-check counts as V5
  only when cited as a library benchmark.

## Benchmark inventory
| value / table | source (citation) | stored in | used by |
|---|---|---|---|
| k_B, N_A, R (exact) | CODATA 2018/2022, NIST CUU | `reference/ch01/constants.json` | ch01 `thermo` constants |
| USSA-1976 constants, composition, layer gradients, errata | NASA-TM-X-74335 (NTRS 19770009539) | `reference/ch01/ussa1976_constants.json` | ch01 `statics.standard_atmosphere`, `sutherland_viscosity` |
| USSA-1976 Tables 1–2 (T, p, ρ, c, g; μ, ν, H_p, n, speed, M; 0–50 km) | PDAS tables computed from NASA-TM-X-74335 | `reference/ch01/ussa1976_table{1,2}.csv` | ch01 C13, C14, C20, C23, C36, C40, C49, C63 |
| IAPWS R1-76(2014) σ(T); R12-08 μ(T) | iapws.org | `reference/ch01/iapws_sigma.csv`, `benchmarks.json` | ch01 C13, C15 |
| Jennings mean free path 67.3 nm | Tsalikis et al., Aerosol Sci. Technol. 58 (2024) | `reference/ch01/benchmarks.json` | ch01 C07/N03 |
| Taylor blast K(γ = 1.4) = 0.856 | Díaz, arXiv:2009.05674 | same | ch01 C75 |
| AMS dry-adiabatic lapse rate ≈ 9.8 K/km | AMS Glossary of Meteorology | same | ch01 C54 |
| UNESCO EOS-80 check values | Fofonoff & Millard, Unesco Tech. Pap. Mar. Sci. 44 (1983) | same | ch01 C60 |
| Divergence theorem, F = (2x, y², z²) through the unit sphere = 8π/3 | Wikipedia "Divergence theorem" | `reference/ch02/benchmarks.json`, `SOURCES.md` | ch02 C14 (rel 4e-16 / 8e-16) |
| Levi-Civita identities (ε–δ, 2δ, 6, cross product, det, pseudotensor) | Wikipedia "Levi-Civita symbol" | same | ch02 C08 |
| Cauchy index placement T_j = σ_ij n_i; Mohr principal stresses and τ_max | Wikipedia "Cauchy stress tensor" | same | ch02 C05, C13 |
| Active R(θ), passive = Rᵀ, Rodrigues; scipy `Rotation.from_rotvec` | Wikipedia "Rotation matrix"; scipy docs | same | ch02 C02 (library cross-check) |
| Stokes' theorem statement and right-hand rule | Wikipedia "Stokes' theorem" | same | ch02 C16 |
Book-printed values (private, git-ignored): `tests/book_values_ch01.json` (constants, Γ_a, scale height, mean free
path, pipe matrix, Examples 1.2–1.5), `tests/book_values_ch02.json` (Examples 2.1–2.4, Exercise 2.1a).

## Validation summary per chapter
Counts are concept-map rows carrying each label (a row can carry several); tests per evidence level in parentheses.
| chapter | analytic | symbolic | converged | conserved | benchmark | book-value | qualitative | unverified | verdict |
|---|---|---|---|---|---|---|---|---|---|
| ch01 | 39 (V1 58) | 16 (V2 19) | 5 (V3 8) | 5 (V4 7) | 18 (V5 13) | 1 (V6 5) | 1 (`wavelength_to_rgb`) | 0 | **PASS** — 107 tests; every computable A item ≥ 2 independent levels; 144 Part C callables tested; FTCS order 2.002; review Must-fix 2/2 closed; 5 explainers PASS (115 parity rows); notebook PASS |
| ch02 | 43 (V1 41) | 21 (V2 13) | 11 (V3 15) | 2 (V4 1) | 9 (V5 6) | 8 (V6 4) | 0 | 0 | **PASS** — 92 tests (+ V7 12); 44 concept rows; every A item ≥ 2 independent levels; 143 Part C names contract-tested; stencil orders 1.97–2.03, integral definitions 1.997–1.998, Simpson 4.00; 8/8 wrong-convention variants caught; review Must-fix 2/2 closed (M1 R = G − Gᵀ, M2 presets); 5 explainers PASS round 2 (127 parity rows, 1 496 views); notebook PASS round 2 (405 cells, 59 s) |
(ch01 V7 limits/invariance: 19 tests; smoke: 2. ch02 V7: 12 tests. Full suite at the ch02 merge gate: 207 tests.)

## Open across chapters
- Machinery TODO (orchestrator): `tools/shot.py` clicks/drags every view and fails on page errors; `tools/coverage_check.py`
  fails loudly on a blank tier cell and understands ledger reminder rows (28 false warnings in ch02); `nb.derivation`
  skips the "Eq." prefix for non-numeric labels; `scripts/ch02_drawings.py` needs a `sys.path` bootstrap (or helper
  modules are excluded from the "every script exits 0" check).
- Pending promotions (a pass where nothing reads `viz/`): the JS candidates in `knowledge/viz_patterns.md`; the
  `core.anim` constrained-layout wrapper; skill lesson lines listed there under "Lesson candidates for the skills".
- Non-blocking follow-ups: ch01 E1–E5 and the book-map figure (`knowledge/ch01.md` §9); ch02 E1 τ₁₂ label, E3/E5
  portrait x-range, E5 label overlap (`knowledge/ch02.md` §9).
