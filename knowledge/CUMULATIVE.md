# CUMULATIVE knowledge — fluidpy
(Rewritten by the knowledge-keeper after every chapter. Every agent reads this first. Last rewrite: after ch01,
2026-09-13, commit 19507fe.)

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
     │                                                        │
     │                                                        └─ stratification: parcel ζ″ + N²ζ = 0 (D18) [C50] → N² (1.29) [C51]
     │                                                              → Γ_a = −gαT/C_p (1.30, D19) [C54] → θ, N² = (g/θ)dθ/dz (D20, D36) [C55]
     │                                                              (Kundu Γ ≡ dT/dz in code; meteorology −dT/dz always shown)
     └─ dimensional analysis: homogeneity [C64] → dimensional matrix, rank (1.39) [C67] → Π groups = null space (1.37, D28) [C69]

next: ch02 Cartesian tensors → ch03 kinematics → ch04 conservation laws / Navier–Stokes (the trunk)
```
The book's structure ahead: Ch2 tensors → Ch3 kinematics → **Ch4 conservation laws / Navier–Stokes (the trunk)** → Ch5
vorticity · Ch6 ideal flow · Ch7 gravity waves · Ch8 laminar flow → Ch9 boundary layers → Ch10 CFD → Ch11 instability →
Ch12 turbulence → **Ch13 GFD (uses 4, 5, 7, 8, 11, 12)** → Ch14 aerodynamics · Ch15 compressible · Ch16 biofluids.
Where ch01 feeds in: N² and θ → Ch. 4 Boussinesq, Ch. 7 §7.8, Ch. 11 §11.7, Ch. 13; Newtonian stress and ν → Ch. 4 §4.5,
Ch. 8, 9; hydrostatic base state → Ch. 4, 7, 13; Π groups → Ch. 4 §4.11 and every scaling; isentropic gas → Ch. 15;
FTCS diffusion → Ch. 5, 8, 10.

## Available primitives
### Machinery (`fluidpy/core/`, ready before chapter 1)
| module | what | used for |
|---|---|---|
| `project` | repo root, `book.yaml` config, Pages/raw/Colab URLs | every tool and notebook |
| `style` | `setup_notebook()` (matplotlib style, plotly renderer, FAST), `COLORS`, `savefig` | every notebook figure |
| `embed` | `show_viz(chapter, slug)` — explainer in Jupyter (srcdoc) / Colab (Pages src) / page | explainer cells |
| `anim` | `animate`, `show_animation(player="video"/"frames")` | animations |
| `interact` | `slider_figure`, `animate_figure` (plotly, work on the page), `live` (ipywidgets) | Python interactives |
| `units` | pint registry `ureg`/`Q_`, `dimensional_check`, `check_dimensions`, `DIM`; ch01 added `celsius_to_kelvin`, `kelvin_to_celsius`, `nondimensional` | V2 tests, temperatures |
| `refdata` | reference-data resolution + registry | benchmarks |
| `_util` (private, ch01) | `as_scalar_if_0d` (scalar in → float out, needed by explainer parity rows), `require_positive`, `require_nonnegative` | every physics function |

### Physics primitives (all created in ch01; re-exported as `ch01.<name>`; labels from `reports/ch01_verification.md`)
| function | module | book § / Eq. | used by chapters | validation label |
|---|---|---|---|---|
| `K_B, N_A, N_A_KMOL, R_U, M_W_AIR, R_AIR, GAMMA_AIR, CP_AIR, CV_AIR, G0, P_ATM, P_REF` | `thermo` | §1.9 (CODATA, USSA-1976) | ch01 | benchmark |
| `perfect_gas_density/pressure/state`, `gas_constant`, `molecular_gas_pressure`, `molecule_mass` | `thermo` | (1.21), (1.22) | ch01 | benchmark, analytic |
| `specific_volume`, `enthalpy`, `perfect_gas_internal_energy/enthalpy`, `specific_heat_cp/cv`, `partial_derivative` | `thermo` | (1.12)–(1.15) | ch01 | analytic, converged |
| `process_path`, `path_heat_work_totals`, `process_heat_work`, `join_paths` | `thermo` | (1.10), (1.11) | ch01 | analytic, converged, conserved |
| `perfect_gas_entropy_change(_p)`, `entropy_change_reversible`, `irreversible_process`, `free_expansion`, `stirred_isochoric_process` | `thermo` | (1.16)–(1.18), Clausius–Duhem | ch01 | symbolic, analytic, conserved |
| `sound_speed_from_eos`, `perfect_gas_sound_speed`, `tait_pressure` | `thermo` | (1.19), (1.27) | ch01 | analytic, benchmark |
| `thermal_expansion_coefficient`, `perfect_gas_expansion_coefficient` | `thermo` | (1.20), (1.28) | ch01 | analytic |
| `cv_from_cp`, `gamma_from_cp`, `cp_from_gamma`, `cv_from_gamma` | `thermo` | (1.23), (1.24) | ch01 | analytic |
| `isentropic_pressure`, `isentropic_ratios` | `thermo` | (1.25), (1.26) | ch01 | symbolic, converged, analytic |
| `van_der_waals_pressure/internal_energy/constants_per_mass` | `thermo` | §1.9 (counter-example) | ch01 | analytic |
| `gauge_pressure`, `absolute_pressure`, `hydrostatic_pressure_uniform`, `layered_pressure` | `statics` | (1.9) | ch01 | analytic |
| `integrate_hydrostatic(z, rho_fn, p0)` | `statics` | (1.8) | ch01 | analytic, converged |
| `atmosphere_from_temperature`, `linear_lapse_pressure`, `standard_atmosphere(z, geometric=False)` | `statics` | §1.10, USSA-1976 | ch01 | benchmark |
| `isothermal_pressure/density`, `scale_height` | `statics` | §1.10 | ch01 | symbolic, benchmark |
| `buoyancy_force`, `net_pressure_force_on_box` | `statics` | D37 | ch01 | analytic, conserved |
| `brunt_vaisala_sq`, `…_from_lapse`, `…_from_theta`, `classify_stability`, `stability_timescale` | `stratification` | (1.29) | ch01 | analytic, symbolic |
| `parcel_displacement`, `parcel_ode`, `parcel_ode_from_gradients`, `parcel_ode_atmosphere` | `stratification` | §1.10, D18 | ch01 | symbolic, analytic, conserved |
| `adiabatic_lapse_rate` (negative), `lapse_rate`, `parcel_temperature`, `lapse_rate_convention`, `lapse_rate_stability` | `stratification` | (1.30), (1.32) | ch01 | symbolic, analytic, benchmark |
| `potential_temperature`, `temperature_from_potential`, `potential_temperature_gradient`, `potential_density` | `stratification` | (1.31)–(1.34) | ch01 | analytic, symbolic, conserved |
| `isentropic_density_gradient`, `ocean_potential_density_gradient` | `stratification` | (1.35) | ch01 | analytic |
| `ftcs_diffusion_1d`, `stable_time_step`, `couette_startup_profile`, `gaussian_spreading`, `derivative_2nd_order` | `diffusion` | §1.5 (model PDE ours) | ch01 | converged, conserved, analytic |
| `dimension_vector`, `dimensional_matrix`, `minor_determinant`, `rank_by_minors` | `dimensional` | (1.39) | ch01 | analytic, symbolic |
| `solve_exponents`, `pi_groups`, `groups_independent`, `group_value`, `rescale_units`, `group_latex`; presets `PIPE, PENDULUM, SPHERE_DRAG, SCALE_HEIGHT, BLAST, RAYLEIGH, PYTHAGORAS` | `dimensional` | (1.36)–(1.40) | ch01 | symbolic, analytic |
Chapter-only physics (continuum sampling, Kn, surface tension, capillary rise, μ(T), seawater EOS, Examples 1.2–1.5)
stays in `fluidpy/ch01_introduction.py`; move a function to `core/` when a second chapter calls it.

### Explainer engine (`assets/viz_lib.js`)
`Viz.app` (tabs: Walkthrough / Explore / Explain / Derivation / Equations / Code / Check; fit-to-window with density
levels and pagers; transport, modes, presets, status, terms, inspector, notes, selftest parity), `Plot` (world
coordinates), `Viz.field` (grid, contours, heatmap, streamlines, quiver, particles), `Viz.num` (RK4, odeint, brentq,
erf, niceTicks), `Viz.work` (Explain builders), `Viz.three` (3-D), KaTeX with fallback. **Added after ch01:**
`Viz.font`, `Viz.roundRect`, `Viz.text` (pixel label with halo/bg), `Viz.card` (wrapped end-of-run card),
`Viz.fmtTime`; `Viz.tnum` mantissa fix + `keepTiny`. Pointer events: `ev.viewId`, `ev.vizView` (never assign
`ev.view`). Reference: `templates/viz_example.html`.

## Explainer inventory
| chapter | slug | CORE idea | reusable stage pattern |
|---|---|---|---|
| ch01 | `continuum_averaging_volume` | C06 continuum (+Kn, D35) | molecules view · measured value vs log box size with expected band · regime strip; resample + click-a-sample inspector; hidden-view number in a visible title |
| ch01 | `viscosity_momentum_diffusion` | C12 Newton viscosity, ν | gap with shearing tracers · profile with steady ghost and fading trail · wall stress vs time; modes = different diffusivity; draggable plate; end-of-run card |
| ch01 | `heat_work_paths` | C25, C35, C45 | piston · p–v · T–s linked on a path parameter; shaded/hatched areas; term bars (path vs state functions); reversible/irreversible modes |
| ch01 | `parcel_stability` | C50, C51, C54, C55 | column coloured by θ/ρ · profile with parcel adiabat · ζ(t) with linear ghost, one clock; two-convention badge with exact-text parity; atmosphere/ocean/lab modes; y-title strip |
| ch01 | `buckingham_pi_machine` | C64, C67, C69 | variable chips → matrix with highlighted minor and determinant → groups; exact rational JS; unit-system switch; problem presets |

## Notation
See `knowledge/notation.md` (symbol register). Conventions that matter everywhere: SI and kelvin inside functions;
**z up**; **lapse rate Kundu Γ ≡ dT/dz in code, meteorology −dT/dz always shown alongside** (negate and flip the
inequality); `p0` (pressure at z = 0) ≠ `p_ref` (θ reference, 1000 hPa); q heat added to / w work done on the system;
gas constants per kmol; pressures absolute unless `_gauge`; signed τ_xy = μ ∂u/∂y. Symbols with several meanings
already in ch01: α, θ, T, h, R, σ/S, c, D, n, k, q, γ, λ, τ, ε, H.

## Teaching lessons
- **Depth is tiered; coverage is exhaustive** (adopted 2026-09-12, during ch01). The first ch01 curation made every new
  idea CORE at full depth: 76 items, 37 derivations. The design and implement agents then each ran ~50 min and ~500k
  tokens without finishing, and everything would have arrived at the same volume, so the ideas that carry the rest of
  the book did not stand out. Now every inventory row still gets an ID and a depth in an auditable chapter map, but
  only 12–18 load-bearing ideas get the full treatment (A: picture → question → derivation → worked number → code →
  figure). The rest are stated and explained inside the nearest A block (B), or named with a pointer to where they
  are used later (C). A derivation is written out only if it belongs to an A item or the book never writes it out.
  Policy: `book.yaml → policy`; rule: `.claude/agents/concept-curator.md`, `/do-chapter` §2, CLAUDE.md rule 6.
  ch01 result: 15 A / 70 B / 21 C, 12 derivations, 496-cell notebook — passed review in 2 rounds.
- **Convention callouts** ("Which p_o?", the lapse-rate table): state the book's, the field's and the code's
  convention, give the size of the difference in numbers, say which one the code uses and why (ch01 C54, D20).
- **★★★ sympy checks re-run the derivation's own construction** (D28 builds the groups exactly as step 9 says and
  checks they form a basis; D19 tests step 12 for an arbitrary Gibbs function), not just the final identity.
- **Make an approximation measurable**: a linear-vs-unlinearised table whose gap grows with amplitude (ch01 C50).
- **Every prose number is computed, not typed** (ch01 round 2: "40 km" should have been 17.9 km from USSA); step
  pointers ("step 12 of D28") are checked by script after any renumbering.
- **Primers come before the first derivation that uses them** (∂, ∇, complex roots, Maxwell relations were round-1
  Must-fixes in ch01); `knowledge/primers.md` lists the 61 already written.

## Global pitfalls confirmed in this project
**Machinery and environment**
- Extracted text garbles maths (`¼` = `=`, `ð…Þ` = parentheses, missing minus, ω → `u`, γ → g, θ → q, σ/τ → s, ν → n,
  ε → 3); matrices and minus signs vanish — read equations from rendered pages.
- Extracted text contains ligatures in the PDF; `tools/split_pdf.py` expands them (grep "fluid" works).
- Windows: shell heredocs mangle backslashes — write code files with Write/Edit; console needs UTF-8 (settings.json env).
- Anaconda's `python3` kernelspec can shadow the venv's — notebooks execute on `fluidpy-venv`.
- JS: `UIEvent.view` is a read-only getter; assigning it threw inside every pointer handler and no audit clicked.
  Until `tools/shot.py` clicks and drags every view, a reviewer must run a click/drag pass and count page errors.
- `tests/test_machinery.py::test_viz_library_inlined_and_template_lints` fails whenever `assets/viz_lib.js` changed and
  `tools/viz_inline.py --all` has not been run.
- Book wording can hide in docstrings (`tools/check_public.py` does not catch paraphrase): reviewers scan for it.

**Mathematics → code**
- A default at a regime boundary must be computed from the boundary function (−9.8e-3 K/m typed as "neutral" was
  unstable against Γ_a = −9.7607e-3).
- A test suite can pass on wrong physics: prove discrimination by patching a scratch copy (parcel ÷ρ_environment passed
  six tests; an exact energy first integral caught it).
- Sign conventions: lapse rate (Kundu vs meteorology), first-law q/w, signed wall stress, signed principal radii.
- kmol vs mol (factor 1000) in every gas constant; pint's [substance] dimension must be dropped consistently.
- pint °C is an offset unit — kelvin inside functions.
- Reference pressure `p_ref` vs surface pressure `p0`; geopotential vs geometric altitude (USSA tables are geopotential).
- Old equations of state use old temperature scales (EOS-80 in IPTS-68: t68 = 1.00024 t90).
- `np.gradient` is first order at the edges unless `edge_order=2`; never measure a scheme's order with it.
- Unbounded ODE growth (N² < 0): cap with a `solve_ivp` event, return NaN beyond it.
- Neutral tolerance bands differ between sibling functions — exclude the band when sweeping.
- Sibling functions with different argument orders — call by keyword.
- Statistical checks: seed, 5 standard errors.
- Validation labels: V4 only for real conservation laws or state-function invariants; unit/convention invariance is V7;
  identities are V1.

## Benchmark inventory
| value / table | source (citation) | stored in | used by |
|---|---|---|---|
| k_B, N_A, R (exact) | CODATA 2018/2022, NIST CUU | `reference/ch01/constants.json` | ch01 `thermo` constants |
| USSA-1976 constants, composition (M0 = 28.9644), layer gradients, errata (S = 110.4 K, r0) | NASA-TM-X-74335 (NTRS 19770009539) | `reference/ch01/ussa1976_constants.json` | ch01 `statics.standard_atmosphere`, `sutherland_viscosity` |
| USSA-1976 Table 1 T, p, ρ, c, g (0–50 km) | PDAS tables computed from NASA-TM-X-74335 | `reference/ch01/ussa1976_table1.csv` | ch01 C20, C36, C49 |
| USSA-1976 Table 2 μ, ν, H_p, n, mean speed, M | same | `reference/ch01/ussa1976_table2.csv` | ch01 C13, C14, C23, C40, C63 |
| IAPWS R1-76(2014) σ(T) of water | iapws.org Surf-H2O-2014 | `reference/ch01/iapws_sigma.csv`, `benchmarks.json` | ch01 C15 |
| IAPWS R12-08 μ(298.15 K), μ(373.15 K) | iapws.org viscosity release Table 4 | `benchmarks.json` | ch01 C13 |
| Jennings mean free path 67.3 nm (air, 300 K, 1 atm) | Tsalikis et al., Aerosol Sci. Technol. 58 (2024), doi:10.1080/02786826.2024.2333859 | `benchmarks.json` | ch01 C07/N03 |
| Taylor blast K(γ = 1.4) = 0.856 | Díaz, arXiv:2009.05674 (quoting Taylor 1950) | `benchmarks.json` | ch01 C75 |
| AMS dry-adiabatic lapse rate ≈ 9.8 K/km | AMS Glossary of Meteorology | `benchmarks.json` | ch01 C54 |
| UNESCO EOS-80 check values ρ(S, t68, 0) | Fofonoff & Millard, Unesco Tech. Pap. Mar. Sci. 44 (1983) | `benchmarks.json` | ch01 C60 |
Book-printed values (private, git-ignored): `tests/book_values_ch01.json` (constants, Γ_a, scale height, mean free
path, pipe matrix and Examples 1.2–1.5).

## Validation summary per chapter
Counts are concept-map rows carrying each label (a row can carry several); tests per evidence level in parentheses.
| chapter | analytic | symbolic | converged | conserved | benchmark | book-value | qualitative | unverified | verdict |
|---|---|---|---|---|---|---|---|---|---|
| ch01 | 39 (V1 58) | 16 (V2 19) | 5 (V3 8) | 5 (V4 7) | 18 (V5 13) | 1 (V6 5) | 1 (`wavelength_to_rgb`) | 0 | **PASS** — 107 tests; every computable A item ≥ 2 independent levels; 144 Part C callables all tested; FTCS order 2.002; review Must-fix 2/2 closed; 5 explainers PASS (115 parity rows); notebook PASS |
(V7 limits/invariance: 19 tests; smoke: 2.)

## Open across chapters
- Machinery TODO (orchestrator): `tools/shot.py` clicks/drags every view and fails on page errors; `tools/coverage_check.py`
  fails loudly on a blank tier cell.
- Non-blocking ch01 follow-ups (explainers E1–E5, book-map figure): `knowledge/ch01.md` §9.
