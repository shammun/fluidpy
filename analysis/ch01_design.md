# Chapter 1 — Introduction: lesson design
(from the user-approved, **re-curated** `analysis/ch01_curation.md` (A 15 · B 70 · C 21 = 106 rows; CORE 15 · RECAP 3 ·
NOTE 86 · SKIP 2; 12 derivations; E1–E5 + backup B1) and `analysis/ch01.md`; equations re-read on the rendered pages
p034–p052 (`tools/render_pages.py ch01 --eq 1.8 / 1.10 / 1.11 / 1.18 / 1.22 / 1.25 / 1.29 / 1.30 / 1.31 / 1.37 / 1.39`);
signatures checked against the committed WIP `fluidpy/` by introspection; 2026-09-12, lesson-designer. Rewrites the
earlier partial design (Part C only, written against the 76-CORE curation).)

**Binding conventions for every builder.**
1. **Lapse rate (user decision).** Code computes with Kundu's $\Gamma \equiv dT/dz$: `adiabatic_lapse_rate()` returns
   $\Gamma_a = -g/C_p = -9.76\times10^{-3}$ K/m (our constants), arguments are named `dT_dz`, stability is
   $dT/dz > \Gamma_a$. The **meteorology convention** $\Gamma \equiv -dT/dz$ ($\Gamma_a = +9.76$ K/km, stable when
   $\Gamma < \Gamma_a$) is shown **next to it every time** a lapse rate is printed, plotted or explained: a two-row
   convention table in C54, both numbers in every E4 readout, a Kundu ↔ meteorology toggle in E4 that changes only the
   bookkeeping. Conversion is taught as "negate the number **and** flip the inequality". Helpers:
   `lapse_rate_convention`, `lapse_rate_stability` (Part C, NEW).
2. **z is upward** in §1.7 and §1.10; $p_0$ = pressure at $z = 0$; depth $h = -z$. The reference pressure of θ and
   $\rho_\theta$ is `p_ref` $= p_o = 1.0\times10^5$ Pa and is never confused with $p_0$.
3. **Constants are ours** (CODATA $k_B$, $N_A$; USSA-1976 $M_w = 28.9644$ kg/kmol; γ = 1.4 ⇒ $R = 287.06$,
   $C_p = 1004.70$, $C_v = 717.64$ J kg⁻¹ K⁻¹; $g = 9.80665$ m s⁻²). Tiny examples may round ($g \approx 9.81$) and
   say so. Numbers quoted by the book never appear in cells, explainers or this file; every number below was computed
   with the WIP `fluidpy` (the "expected output" values).
4. **Suspected typos taught as analysed (analysis §9):** Clausius–Duhem with the actual heat
   $s_2 - s_1 \ge \int \delta q/T$; (1.35) as "has the same sign as"; Ex. 1.1's α is the complement of the usual contact
   angle; the Step-5 "(1.32)" reference means (1.39).
5. **Colours (text, figures, explainers):** pressure / heat `orange` · viscous / friction / irreversible `rose` ·
   inertia / motion / parcel `teal` · work `blue` · state functions (Δe, Δs, θ) `accent` purple · stable `teal`,
   neutral `amber`, unstable `rose` · molecules and ghosts `muted` grey · environment profile `text` (bold black).
6. **Depth.** A items (CORE) get the full nine-part block. B items are `nb.note` paragraphs (plain words + equation +
   one number) inside the A block named in the curation (or, in §1.1/1.2/1.3/1.6, in their own section tagged with
   the A parent). C items are one `nb.note` sentence with the pointer. RECAP → `nb.recap`; SKIP → `nb.pointer`.
   The 25 demoted derivations (§4c) are **stated, never derived**. Only the 12 D rows are written out (Part F).
7. **Parsing.** Part E is the only part with table rows after the Part E heading; Part F has no lines that start
   with `|`, so `tools/coverage_check.py` reads only the ledger. Explainer headings are exactly `### E1 · <slug>` …
   `### B1 · <slug>` (`tools/embed_check.py`).

Order of parts in this file: C (the implementer's contract, written first) · A (notebook storyboard) · B (explainers) ·
D (runtime) · E (prerequisite ledger) · F (derivations).

---

## Part C — functions the builders will call (the implementer's contract)

**Status column.** **WIP** = already in the committed WIP with exactly this signature (keep it; test it at ≥ 1 evidence
level, ≥ 2 for A items). **NEW** = must be added. **§4 F·n** = planned in `analysis/ch01.md` §4; *not in §4* = added by
the curation or this design (flagged, as the phase asks). Every callable is reachable as `ch01.<name>` (the chapter module
re-exports `core`), which is what `tools/shot.py` parity rows need: its `py:` expressions are evaluated with **no
builtins** (no `len`, `float`, `list`) and the result must convert with `float(...)`, so parity rows index dicts/tuples
down to a number (`ch01.path_heat_work_totals(…)["w"]`, `ch01.stability_timescale(1e-4)[1]`).

Defaults everywhere: `g=G0` (9.80665 m s⁻²), `R=R_AIR`, `gamma=GAMMA_AIR` (1.4), `cp=CP_AIR`, `cv=CV_AIR`,
`p_ref=P_REF` (1.0e5 Pa), `p_atm=P_ATM` (101325 Pa). SI in and out, kelvin inside, numpy-vectorised where an argument
can be an array.

### C.0 Machinery (existing, `fluidpy/core/`)
| # | Callable (signature) | Returns | Used by | Status |
|---|---|---|---|---|
| 0.1 | `style.setup_notebook() -> bool` (FAST flag), `style.COLORS`, `style.savefig(fig, chapter, name)` | FAST · palette dict · path | setup cell, every figure | WIP |
| 0.2 | `anim.animate(update, frames=60, fig=None, interval=50, blit=False, init=None, repeat=True)` · `anim.show_animation(anim, player="auto"\|"video"\|"frames")` | FuncAnimation · HTML | C06, C12, C25, C50 | WIP |
| 0.3 | `interact.slider_figure(fn, name, values, *, unit, xlabel, ylabel, title, xrange, yrange, active, height, modes)` · `interact.animate_figure(frame_fn, times, …)` · `interact.live(fn, **widgets)` | plotly Figure · widget | C06, C20, C40, C54, C55 · C51 live | WIP |
| 0.4 | `embed.show_viz(chapter, slug)` | display | 5 explainer cells | WIP |
| 0.5 | `units.Q_`, `units.ureg`, `units.DIM`, `units.check_dimensions(value, expected)`, `units.dimensional_check(fn, expected, **kw)` | pint objects · bool · result | R01, R02, C64 | WIP |
| 0.6 | `units.celsius_to_kelvin(T_C)` · `units.kelvin_to_celsius(T_K)` (⟲ `ch01.`) | K · °C | R03, C13, C36 | WIP · §4 F1 |

### C.1 `fluidpy/core/dimensional.py` (⟲ all re-exported)
| # | Callable (signature) | Returns [unit] | Book | Used by | Status |
|---|---|---|---|---|---|
| 1 | `BASIS = ("M", "L", "T", "Θ")` · `LATEX_SYMBOLS` · `UNIT_SYSTEMS` · `PRESET_INFO` | constants | §1.11 Step 2 | C66, E5 | WIP · not in §4 |
| 2 | `dimension_vector(q, basis=BASIS, drop_substance=True) -> ndarray[int]` | exponents of M, L, T, Θ (warns when [substance] is dropped) | §1.11 Step 2 | C64, C66, E5 | WIP · §4 F2 |
| 3 | `dimensional_matrix(variables: Mapping[str, str], basis=BASIS, drop_zero_rows=True) -> (A: ndarray[int], names: list[str], rows: list[str])` | integer matrix, columns in dict order | Eq. (1.39) | C67, C68, C69, C73–C76, E5 | WIP · §4 F3 (numpy, not sympy) |
| 4 | `minor_determinant(A, rows: Sequence[int], cols: Sequence[int]) -> int \| Fraction` | exact determinant (cofactor expansion) | §1.11 Step 3 | C67 worked number, C68, E5 | WIP · not in §4 |
| 5 | `rank_by_minors(A) -> (r: int, rows: tuple, cols: tuple)` | rank + first nonzero r×r witness (lexicographic; PIPE → cols (0, 1, 4)) | §1.11 Step 3 | C67, C68, C73, C76, E5 | WIP · §4 F4 |
| 6 | `solve_exponents(target: str, repeating: Sequence[str], variables) -> dict[str, Fraction]` | `{target: 1, rep: a, …}` (zero exponents kept); `ValueError` if the repeating minor is singular | §1.11 Step 5 | C69 (C70), E5 | WIP · §4 F5 |
| 7 | `pi_groups(variables, solution: str \| None = None, repeating: Sequence[str] \| None = None) -> list[dict[str, Fraction]]` | n − r groups; group 1 holds `solution` to power 1 | Eq. (1.37), (1.40) | C69, C72–C76, E5 | WIP · §4 F6 |
| 8 | `group_latex(group, symbols=None) -> str` · `group_dimension(group, variables) -> ndarray` · `group_expression(group, labels=None) -> sympy.Expr` | LaTeX · (M, L, T, Θ) exponents (zero for a Π) · sympy | §1.11 | C69–C76 printouts, D28 check | WIP · §4 F6 |
| 9 | `groups_independent(groups, names=None) -> bool` · `exponent_matrix(groups, names) -> ndarray` | rank(E) == len(groups) · rows = groups | §1.11 combining | C71 | WIP · §4 F7 |
| 10 | `group_value(group, values: Mapping[str, float]) -> float` | ∏ valueᵏ [–] | §1.11 | C64, C69, E5 | WIP · not in §4 |
| 11 | `rescale_units(values, variables, system="cgs" \| "SI" \| "imperial" \| Mapping) -> dict[str, float]` | the same physical values as numbers in another unit system | §1.11 opening | C64, C69 invariance figure, E5 | WIP · not in §4 |
| 12 | presets `PIPE`, `PENDULUM`, `SPHERE_DRAG`, `SCALE_HEIGHT`, `PYTHAGORAS`, `BLAST`, `RAYLEIGH` (`dict[str, str]` name → SI unit) | — | (1.38), Ex. 1.2–1.5 | C65–C76, E5 | WIP · not in §4 |

Preset keys used by every cell and explainer: `PIPE = {dp: Pa, dx: m, d: m, eps: m, U: m/s, rho: kg/m**3, mu: Pa*s}` ·
`SCALE_HEIGHT = {H: m, T0: K, Mw: kg/kmol, g: m/s**2, Ru: J/(kmol*K)}` · `BLAST = {E: J, D: m, rho: kg/m**3, t: s}` ·
`RAYLEIGH = {S: W/m**2, I: W/m**2, lam: m, V: m**3, n_s: dimensionless, d: m}` · `PENDULUM = {tau: s, Lp: m, m: kg,
g: m/s**2, theta0: rad}` · `SPHERE_DRAG = {F: N, D: m, U: m/s, rho: kg/m**3, mu: Pa*s}` · `PYTHAGORAS = {a: m**2,
beta: rad, C: m}`.

### C.2 `fluidpy/core/thermo.py` (⟲ all re-exported)
| # | Callable (signature) | Returns [unit] | Book | Used by | Status |
|---|---|---|---|---|---|
| 13 | constants `K_B` [J/K], `N_A` [1/mol], `N_A_KMOL` [1/kmol], `R_U` [J kmol⁻¹ K⁻¹], `M_W_AIR` = 28.9644, `R_AIR` = 287.058, `GAMMA_AIR` = 1.4, `CP_AIR` = 1004.70, `CV_AIR` = 717.64, `G0`, `P_ATM`, `P_REF`, `MOLAR_MASS`, `GAMMA_BY_ATOMICITY`, `VDW_CO2` | — | §1.9 | everywhere | WIP · §4 F8 |
| 14 | `molecular_mass(M_w)` (alias `molecule_mass`) | kg | §1.9 | C06 (D34), C38 | WIP · not in §4 |
| 15 | `molecular_gas_pressure(n, V, T)` | Pa | Eq. (1.21) | C38 | WIP · §4 F9 |
| 16 | `gas_constant(M_w)` | J kg⁻¹ K⁻¹ | (1.22) | C39, C40, D11 | WIP · §4 F10 |
| 17 | `perfect_gas_pressure(rho, T, R=R_AIR)` · `perfect_gas_density(p, T, R=R_AIR)` | Pa · kg/m³ | Eq. (1.22) | C40, C49, E3 | WIP · §4 F10 |
| 18 | `perfect_gas_state(p=None, rho=None, T=None, R=R_AIR) -> (p, rho, T)` | exactly one argument None | (1.12), (1.22) | C28, C40 3-D figure | WIP · §4 F11 |
| 19 | `specific_volume(rho)` · `enthalpy(e, p, v)` · `perfect_gas_internal_energy(T, cv=CV_AIR, T_ref=0, e_ref=0)` · `perfect_gas_enthalpy(T, cp=CP_AIR, T_ref=0, h_ref=0)` | m³/kg · J/kg · J/kg · J/kg | §1.8, (1.13) | N13, C29, C35 | WIP · §4 F12 (+ not in §4) |
| 20 | `partial_derivative(f, var, at: Mapping, rel_step=1e-6, step=None)` · `specific_heat_cp(h_of_Tp, T, p)` · `specific_heat_cv(e_of_Tv, T, v)` | ∂f/∂var · J kg⁻¹ K⁻¹ | (1.14), (1.15) | C30, C31, C36 primer | WIP · §4 F13, F14 |
| 21 | `cv_from_cp(cp, R)` · `gamma_from_cp(cp, R)` · `cp_from_gamma(gamma, R)` | J kg⁻¹ K⁻¹ · – | (1.23), (1.24) | C42, C43, D14 | WIP · §4 F20 |
| 22 | `process_path(kind, state1: (v, T), state2: (v, T), n=401, gamma, R) -> dict(v, T, p, corner)` — kinds `isothermal`, `isochoric`, `isobaric`, `isentropic`, `isochoric-isobaric`, `isobaric-isochoric`, `isentropic-isochoric`, `isothermal-isochoric` | arrays [m³/kg, K, Pa] | §1.8 | C25, C26, E3 | WIP · not in §4 |
| 23 | `process_heat_work(v_path, T_path, cv, R) -> dict(w, de, q, p)` | cumulative arrays [J/kg]; `w = −∫p dv` (trapezoid) | (1.10), (1.11) | C25 (library side of from-scratch), E3 | WIP · §4 F15 |
| 24 | `path_heat_work_totals(kind, v1, T1, v2, T2, gamma, R) -> dict(q, w, de, ds, int_dq_over_T)` | exact leg-by-leg totals [J/kg, J kg⁻¹ K⁻¹] | (1.10), (1.11), (1.16), (1.18) | C25, C26, C33, E3 mirror | WIP · not in §4 |
| 25 | `irreversible_process(kind: "stirring" \| "free_expansion", T1, v1, cv, R, T2=None, v2=None) -> dict(q, w, de, ds, int_dq_over_T)` | totals; `ds > int_dq_over_T` | Clausius–Duhem (§1.8 ii) | C34, C44, E3 irreversible mode | WIP · not in §4 |
| 26 | `entropy_change_reversible(q_path, T_path, cumulative=False)` | ∫dq/T [J kg⁻¹ K⁻¹] | (1.16), (1.17) | C33, N15, E3 | WIP · §4 F16 |
| 27 | `perfect_gas_entropy_change(T1, v1, T2, v2, cv, R)` · `perfect_gas_entropy_change_p(T1, p1, T2, p2, cp, R)` | J kg⁻¹ K⁻¹ | integrated (1.18) | C35, C34, E3 | WIP · §4 F17 |
| 28 | `sound_speed_from_eos(p_of_rho_s: Callable[[rho, s], p], rho, s=0.0, rel_step=1e-6)` | m/s (raises if ∂p/∂ρ ≤ 0) | Eq. (1.19) | C36 | WIP · §4 F18 |
| 29 | `tait_pressure(rho, rho0=1000.0, K0=2.2e9, n=7.15, p0=P_ATM)` | Pa | model for (1.19) | C36 | WIP · not in §4 |
| 30 | `thermal_expansion_coefficient(rho_of_Tp: Callable[[T, p], rho], T, p)` | 1/K | Eq. (1.20) | C37 | WIP · §4 F19 |
| 31 | `isentropic_pressure(rho, p0, rho0, gamma)` | Pa | Eq. (1.25) | C45, E3 | WIP · §4 F21 |
| 32 | `isentropic_ratios(p_over_p0, gamma) -> (T_over_T0, rho_over_rho0)` | – | Eq. (1.26) | C46, C55 (D20) | WIP · §4 F22 |
| 33 | `perfect_gas_sound_speed(T, gamma, R)` | m/s | Eq. (1.27) | C47, C36 check | WIP · §4 F23 |
| 34 | `perfect_gas_expansion_coefficient(T)` | 1/K | Eq. (1.28) | C48 | WIP · §4 F24 |
| 35 | `van_der_waals_pressure(T, v, a, b, R)` · `van_der_waals_internal_energy(T, v, a, cv, e_ref=0)` | Pa · J/kg | contrast for C41 | C41 | WIP · not in §4 |

### C.3 `fluidpy/core/statics.py` (⟲ all re-exported)
| # | Callable (signature) | Returns [unit] | Book | Used by | Status |
|---|---|---|---|---|---|
| 36 | `gauge_pressure(p, p_atm=P_ATM)` · `absolute_pressure(p_gauge, p_atm=P_ATM)` | Pa | §1.7 | C17, B1 | WIP · §4 F25 |
| 37 | `hydrostatic_pressure_uniform(z, p0, rho, g)` | Pa (absolute when p0 is) | Eq. (1.9) | C20, C21, D37, B1 | WIP · §4 F26 |
| 38 | `layered_pressure(z, thicknesses, densities, p_surface=P_ATM, g)` (z = 0 at the free surface, z < 0 below) | Pa | (1.9) layer by layer | C20 slider, B1 | WIP · not in §4 |
| 39 | `integrate_hydrostatic(z, rho_fn: Callable[[z, p], rho], p0, g, z0=None, rtol=1e-10, atol=None, method="DOP853")` | p(z) [Pa] | Eq. (1.8) | C20, C40 (C62 check) | WIP · §4 F27 |
| 40 | `net_pressure_force_on_box(p_fn: Callable[[x, y, z], p], box: (x0, x1, y0, y1, z0, z1), n=41) -> ndarray(3)` | N | D37 | C20 (D37 check), B1 | WIP · not in §4 |
| 41 | `buoyancy_force(rho_fluid, volume, g)` | N (up) | D37 | C20, C50, B1 | WIP · not in §4 |
| 42 | `atmosphere_from_temperature(z, T_fn, p0=P_ATM, R, g, z0=None, rtol=1e-10) -> (p, rho, T)` | Pa, kg/m³, K | (1.8) + (1.22) | C49 (in C54), C54 from-scratch | WIP · §4 F28 |
| 43 | `standard_atmosphere(z, g, R) -> (T, p, rho)` (USSA-1976, 0–84.852 km) | K, Pa, kg/m³ | benchmark | C07 slider, C40, C54, E4 preset | WIP · not in §4 |
| 44 | `isothermal_pressure(z, p0, T, R, g)` · `isothermal_density(z, rho0, T, R, g)` · `scale_height(T, R, g)` | Pa · kg/m³ · m | §1.10 (unnumbered) | C62, C63, C73 | WIP · §4 F29, F30 |
| 45 | `linear_lapse_pressure(z, p0, T0, dT_dz, R, g)` | Pa | (1.8)+(1.22), constant Γ | C54 figure, E4 environment | WIP · not in §4 |

### C.4 `fluidpy/core/stratification.py` (⟲ all re-exported)
| # | Callable (signature) | Returns [unit] | Book | Used by | Status |
|---|---|---|---|---|---|
| 46 | `brunt_vaisala_sq(rho0, drho_dz, drho_a_dz, g)` | N² [1/s²] | Eq. (1.29) | C50, C51, C61, E4 ocean/lab | WIP · §4 F31 |
| 47 | `brunt_vaisala_sq_from_lapse(T, dT_dz, cp, g)` = (g/T)(dT/dz + g/C_p) | 1/s² | (1.32) + D36 | C54, E4 | WIP · not in §4 |
| 48 | `brunt_vaisala_sq_from_theta(theta, dtheta_dz, g)` = (g/θ)dθ/dz | 1/s² | D36 | C55, E4 | WIP · not in §4 |
| 49 | `classify_stability(N2, tol=1e-12)` | `"stable"`/`"neutral"`/`"unstable"` (array of str) | §1.10 | C52, C55 figure | WIP · §4 F32 |
| 50 | `stability_timescale(N2, tol=1e-12) -> (kind: "period"\|"efold"\|"none", seconds: float)` | s | §1.10 | C51, E4 status | WIP · not in §4 |
| 51 | `parcel_displacement(t, zeta0, N2, w0=0.0)` | ζ [m]: cos / linear / cosh | D18 solution | C50, C51, E4 | WIP · §4 F33 |
| 52 | `parcel_ode_from_gradients(t, zeta0, rho0, drho_dz, drho_a_dz, g, w0=0.0, zeta_max=None)` | ζ [m] (NaN after the cap) | D18 before linearising | C50 (nonlinear check), C51 live | WIP · §4 F33 |
| 53 | `parcel_acceleration_atmosphere(zeta, T0, dT_dz, cp, g)` · `parcel_ode_atmosphere(t, zeta0, T0, dT_dz, cp, g, w0=0.0, zeta_max=None, rtol=1e-10, atol=1e-12)` | m/s² · ζ [m] | D18, perfect gas | C54 animation, E4 nonlinear path | WIP · not in §4 |
| 54 | `lapse_rate(T, z)` | dT/dz [K/m] (Kundu sign), 2nd order | §1.10 | C54 from-scratch | WIP · not in §4 |
| 55 | `adiabatic_lapse_rate(T=None, cp=CP_AIR, alpha=None, g)` = −gαT/C_p; `alpha=None` → −g/C_p | K/m, **negative** | Eq. (1.30) | C54, N17, D19 check, E4 | WIP · §4 F34 |
| 56 | `parcel_temperature(T0, z, z0=0.0, cp, g)` = T0 + Γ_a(z − z0) | K | (1.30) integrated | C54 figure, E4 adiabat | WIP · not in §4 |
| 57 | `potential_temperature(T, p, p_ref=P_REF, gamma)` · `temperature_from_potential(theta, p, p_ref, gamma)` | K | Eq. (1.31) | C55, E4 | WIP · §4 F35 |
| 58 | `potential_temperature_gradient(T, dT_dz, theta=None, cp, g, p=None, p_ref, gamma)` = (θ/T)(dT/dz + g/C_p) | K/m | Eq. (1.32) | C56, E4 Explain | WIP · §4 F36 |
| 59 | `potential_density(rho, p, p_ref, gamma)` | kg/m³ | Eq. (1.33) | C58 | WIP · §4 F37 |
| 60 | `isentropic_density_gradient(rho, c, g)` = −ρg/c² | kg/m⁴ | §1.10 (N21) | C51 (N21), E4 ocean mode | WIP · §4 F38 |
| 61 | `ocean_potential_density_gradient(drho_dz, rho, c, g)` = dρ/dz + ρg/c² | kg/m⁴ (negative = stable; **same sign** as dρ_θ/dz) | Eq. (1.35) | C61, E4 ocean mode | WIP · §4 F39 |
| 62 | **`lapse_rate_convention(dT_dz, convention: Literal["kundu", "meteorology"] = "kundu") -> float \| ndarray`** — the lapse rate Γ expressed in the requested convention: `"kundu"` returns `dT_dz` unchanged (Γ ≡ dT/dz), `"meteorology"` returns `−dT_dz` (Γ ≡ −dT/dz). Negation is its own inverse, so `lapse_rate_convention(Γ_met, "meteorology")` converts a meteorology number back to dT/dz. Raises `ValueError` for any other convention string. First docstring line: "Convention switch: Γ_met = −Γ_Kundu; convert the number AND flip the stability inequality." | K/m | §1.10 (Γ ≡ dT/dz) + standard meteorology | C53, N17, C54 table and slider legend, E4 every readout | **NEW in this contract** → now in WIP · not in §4 (curation §8 item 2) |
| 63 | **`lapse_rate_stability(dT_dz, Gamma_a=None, convention="kundu", tol=1e-9, per_km=True, prefix=True, digits=None, ascii_minus=False, cp=CP_AIR, g=G0) -> LapseStability`** — `LapseStability` is a `NamedTuple(verdict: str, text: str, code: int, Gamma: float, Gamma_a: float, margin: float)`; first two fields = the curation's `(verdict, inequality_text)`. `Gamma_a=None` → `adiabatic_lapse_rate()` (Kundu sign). Kundu: `Gamma = dT_dz`, `Gamma_a` negative, `margin = Gamma − Gamma_a`, text `"stable ⇔ dT/dz > Γa: −6.5 > −9.8 K/km"`; meteorology: `Gamma = −dT_dz`, `Gamma_a = −Gamma_a_kundu` (positive), `margin = Gamma_a − Gamma`, text `"stable ⇔ Γ < Γa: 6.5 < 9.8 K/km"`. `margin` is identical in both conventions (so `verdict` and `code` never depend on the convention); `code` = +1 stable (margin > tol), 0 neutral (\|margin\| ≤ tol K/m), −1 unstable; for a neutral or unstable case the text shows the actual relation (`"=" `, or the reversed sign, e.g. `"−12.0 < −9.8 K/km"` → unstable). Numbers in the text: K/km with one decimal when `per_km`, Unicode minus "−". Scalar input only. `prefix=False` gives the bare relation ("−12.0 < −9.8 K/km") used by the E4 status line; `digits` rises automatically when two different numbers would print alike. | tuple (str, str, int, K/m, K/m, K/m) | §1.10, (1.30), (1.32) | C54 (code, figure table, slider titles), E4 status + Explain §5 | **NEW in this contract** → now in WIP · not in §4 (curation §8 item 2) |

Both helpers were added to the WIP by the implementer while this design was written (checked by introspection: the texts for −6.5 K/km are exactly the strings above). Both helpers were added to the WIP by the implementer while this design was written (checked by introspection: the texts for −6.5 K/km are exactly the strings above). Both helpers were added to the WIP by the implementer while this design was written (checked by introspection: the texts for −6.5 K/km are exactly the strings above). Both helpers were added to the WIP by the implementer while this design was written (checked by introspection: the texts for −6.5 K/km are exactly the strings above). Tests the verifier adds for 62–63 (curation §8.2): `lapse_rate_convention(adiabatic_lapse_rate(), "meteorology") ==
+g/C_p`; round trip; `lapse_rate_stability(x, convention="kundu").code == lapse_rate_stability(x,
convention="meteorology").code` for a sweep x ∈ [−15, +10] K/km; the code agrees with `classify_stability(
brunt_vaisala_sq_from_lapse(T, x))` mapped to ±1/0; the texts for −6.5 K/km are exactly the two strings above.

### C.5 `fluidpy/core/diffusion.py` (⟲ all re-exported)
| # | Callable (signature) | Returns [unit] | Book | Used by | Status |
|---|---|---|---|---|---|
| 64 | `ftcs_diffusion_1d(f0, D, dy, dt, nsteps, bc=("dirichlet", "dirichlet"), values=None, save_every=1, check_stability=True)` | ndarray (nsave, N) | model PDE (ours; Ch. 4 derives it) | C12 animation, E2 mirror | WIP · §4 F40 |
| 65 | `stable_time_step(D, dy, safety=0.9)` (alias `ftcs_stable_time_step`) | s | FTCS limit | C12 | WIP · not in §4 |
| 66 | `couette_startup_profile(y, t, U, h, nu, nterms=200, tol=1e-12, nmax=20000)` | m/s ("result derived in Ch. 8; used only as a check") | check of (1.3) diffusion picture | C12 check, E2 ghost + parity | WIP · not in §4 |
| 67 | `derivative_2nd_order(f, y)` | df/dy, 2nd order incl. ends | tool | C12 from-scratch reference | WIP · not in §4 |

### C.6 `fluidpy/ch01_introduction.py` (chapter module)
| # | Callable (signature) | Returns [unit] | Book | Used by | Status |
|---|---|---|---|---|---|
| 68 | `shear_deformation_history(t, tau, kind="fluid" \| "solid" \| "bingham" \| "maxwell", G=None, mu=None, t_off=None, tau_y=None)` | γ(t) [rad] | §1.3, (1.3) | C02 sketch (§1.3), N01, E2 tracer lines | WIP · §4 F41 |
| 69 | `traction_components(traction, normal) -> (sigma_n, tau_vec, tau_mag)` (σ_n > 0 = tension) | Pa | §1.3 | C03 | WIP · not in §4 |
| 70 | `number_density(rho, M_w)` · `mean_molecular_spacing(rho, M_w)` | 1/m³ · m | §1.3–1.4 | C04, C06 | WIP · not in §4 |
| 71 | `maxwellian_velocities(n, T, m, rng=None, seed=0) -> ndarray(n, 3)` · `molecular_pressure(number_density, m, velocities)` · `wall_impact_pressure(velocities, m, number_density, dt=1e-12, area=1.0)` | m/s · Pa · Pa | §1.4, D34 | C05, D34 check | WIP · §4 F42 |
| 72 | `box_average_density(L, rho0, variation=0.0, L_flow=1.0, x0=None)` = ρ₀[1 + ε sin(2πx₀/L_flow)·sinc(πL/L_flow)], x₀ default L_flow/4 (a crest, so the point value is ρ₀(1 + ε)) · `density_expected(L, number_density, m, gradient=0.0, L_flow=1.0)` (linear alternative) | kg/m³ (noiseless box average) | §1.4 | C06 code and figure (drift), E1 mirror | WIP · not in §4 |
| 73 | `density_noise_expected(L, number_density)` = (nL³)^(−1/2) | – | D35 | C06, D35 check, E1 mirror | WIP · not in §4 |
| 74 | `sample_density(L_box, number_density, m, rng=None, n_samples=200, gradient=0.0, L_flow=1.0, seed=0, variation=0.0, x0=None) -> (mean, std)` | kg/m³ | §1.4 | C06 (tested side of from-scratch), E1 | WIP · §4 F43 |
| 75 | `knudsen_number(l, L)` · `mean_free_path_jennings(mu, rho, p)` · `mean_free_path_air(T, p=P_ATM)` | – · m · m | §1.4 | C07, N03, E1 | WIP · §4 F44 |
| 76 | `mean_molecular_speed(T, m)` · `collision_time(l, T, m)` | m/s · s | §1.8 | C23, C24 | WIP · not in §4 |
| 77 | `DRY_AIR` (mole fraction, kg/kmol) · `mass_fractions(mole_fractions, molar_masses)` | – | §1.5 | C09 | WIP · not in §4 |
| 78 | `fick_mass_flux(rho, kappa_m, grad_Y)` · `fourier_heat_flux(k, grad_T)` | kg m⁻² s⁻¹ · W/m² | (1.1), (1.2) | C10, C11, E2 modes | WIP · §4 F45, F46 |
| 79 | `newton_shear_stress(mu, du_dy)` · `shear_stress_profile(u, y, mu)` · `wall_shear_history(u_hist, dy, mu) -> (tau_bottom, tau_top)` | Pa | Eq. (1.3) | C12, E2 | WIP · §4 F47 |
| 80 | `viscosity_power_law(T, mu_ref, T_ref, n=0.5)` · `sutherland_viscosity(T, beta=1.458e-6, S=110.4)` · `water_viscosity(T)` | Pa s | §1.5 | C13 | WIP · §4 F48 |
| 81 | `kinematic_viscosity(mu, rho)` · `thermal_diffusivity(k, rho, cp)` · `diffusion_time(L, D)` · `FLUIDS` (`air`, `water`, `glycerine`, `honey`, `engine_oil`: `rho, mu, k, cp, kappa_m, nu`) · `fluid_properties(name, T=293.15, p=P_ATM) -> dict(rho, mu, nu[, sigma])` | m²/s · m²/s · s · dict | Eq. (1.4) | C14, C12 worked number, E2 presets | WIP · §4 F49 |
| 82 | `surface_tension_water(T)` · `laplace_pressure_jump(sigma, R1, R2=None)` | N/m · Pa | §1.6, (1.5) | C15, N07, C16, C22, B1 | WIP · §4 F50, F51 |
| 83 | `wedge_pressure_difference(rho, dz, theta, g) -> dict(p2_minus_p1, p3_minus_p1)` | Pa | (1.6) | C18 number | WIP · §4 F52 |
| 84 | `capillary_rise(sigma, alpha, rho, R, g)` · `alpha_from_contact_angle(theta_c)` | m · rad | Ex. 1.1 | C22, B1 | WIP · §4 F53 |
| 85 | `water_density(T)` · `seawater_density_linear(T, S, rho0=1027.0, alpha_T=1.67e-4, beta_S=7.6e-4, T0=283.15, S0=35.0)` | kg/m³ | data for C37, C60 | C37, C60, E4 ocean mode | WIP · not in §4 |
| 86 | `synthetic_boundary_layer_profile(z, T_surface=288.15, mixed_top=800.0, inversion_depth=200.0, inversion_dT_dz=0.01, mixed_dT_dz=MIXED_LAYER_DT_DZ (= adiabatic_lapse_rate(), neutral mixed layer), upper_dT_dz=-0.0045)` · `synthetic_boundary_layer_column(z, …same…, p_surface, p_ref, cp, R, g, gamma) -> dict(z, T, dT_dz, p, rho, theta, dtheta_dz, rho_theta, N2, stability)` | K · dict | §1.10 (N18) | C55 slider, E4 inversion preset | WIP · not in §4 |
| 87 | `poiseuille_pressure_drop(mu, U, dx, d)` ("result derived in Ch. 8; data generator only") | Pa | data for C72 | C69 collapse figure | WIP · not in §4 |
| 88 | `pythagoras_phi(beta)` · `blast_energy(D, t, rho, K=1.0)` · `blast_radius(E, t, rho, K=1.0)` · `TAYLOR_K_GAMMA14` · `rayleigh_scattering_ratio(V, d, lam, phi3=1.0)` · `wavelength_to_rgb(lam_nm)` | – · J · m · – · – · RGB | Ex. 1.3–1.5 | C74, C75, C76 | WIP · §4 F54, F55 (+ not in §4) |

### C.7 `scripts/` (drawing helpers, pure matplotlib, no physics; imported as `from scripts.ch01_drawings import …`)
| # | Callable (signature) | Returns | Used by | Status |
|---|---|---|---|---|
| 89 | `ch01_book_map.draw_book_map(ax=None, highlight=("ch13",))` | Axes | C01 (§1.1) | WIP · not in §4 |
| 90 | `ch01_drawings.draw_cube_forces(ax)` · `draw_wedge(ax, dz=1.0, theta=0.61)` · `draw_capillary(ax, R=1.0, alpha=1.22, h=3.0)` · `draw_parcel_column(ax, z=None, T_env=None, z0=500.0, zeta=300.0)` · `draw_triangle_split(ax, beta=0.52)` | Axes | C20 (D05), C18, C22, C50, C74 | WIP · not in §4 |

**Count.** 90 contract rows naming **144 callables** (plus constants and presets): 10 machinery (C.0), 6 drawing helpers (C.7) and 128 physics/teaching functions reachable as `ch01.<name>`. **NEW in this contract: 2** (`lapse_rate_convention`, `lapse_rate_stability`) — the implementer added both to the WIP while this design was written, so today **144 of 144 exist** in the WIP with the signatures above. Not planned in `analysis/ch01.md` §4 (flagged): 63 physics/teaching functions (`absolute_pressure`, `alpha_from_contact_angle`, `box_average_density`, `brunt_vaisala_sq_from_lapse`, `brunt_vaisala_sq_from_theta`, `buoyancy_force`, `collision_time`, `couette_startup_profile`, `cp_from_gamma`, `density_expected`, `density_noise_expected`, `derivative_2nd_order`, `diffusion_time`, `exponent_matrix`, `fluid_properties`, `group_dimension`, `group_expression`, `group_latex`, `group_value`, `irreversible_process`, `isothermal_density`, `lapse_rate`, `lapse_rate_convention`, `lapse_rate_stability`, `layered_pressure`, `linear_lapse_pressure`, `mass_fractions`, `mean_free_path_air`, `mean_molecular_spacing`, `mean_molecular_speed`, `minor_determinant`, `molecular_mass`, `net_pressure_force_on_box`, `number_density`, `parcel_acceleration_atmosphere`, `parcel_ode_atmosphere`, `parcel_ode_from_gradients`, `parcel_temperature`, `path_heat_work_totals`, `perfect_gas_enthalpy`, `perfect_gas_entropy_change_p`, `perfect_gas_internal_energy`, `poiseuille_pressure_drop`, `process_path`, `pythagoras_phi`, `rescale_units`, `seawater_density_linear`, `stability_timescale`, `stable_time_step`, `standard_atmosphere`, `synthetic_boundary_layer_column`, `synthetic_boundary_layer_profile`, `tait_pressure`, `temperature_from_potential`, `thermal_diffusivity`, `traction_components`, `van_der_waals_internal_energy`, `van_der_waals_pressure`, `wall_impact_pressure`, `wall_shear_history`, `water_density`, `water_viscosity`, `wavelength_to_rgb`), all curation or design additions. Superseded names of the earlier partial design: none dropped — `stirred_isochoric_path` survives as an alias, `box_average_density` is joined by the WIP's `density_expected`.

## Part A — notebook storyboard (`notebooks/build_ch01.py` → `notebooks/ch01_introduction.ipynb`)

**How to read this part.** Each numbered row is **one** `nbkit` call, in order. Flags in brackets say which optional
arguments are used: `[code]` a primer with a runnable demo, `[explain]` a "What does the code above do?" cell,
`[check]` a derivation with `check_src`. Cell numbers below are computed from these rows (section cell = the
`nb.section` heading; primer[code] = 2 cells, code[explain] = 2, figure = 2 (+1 explain), animation/plotly[explain] = 2,
live[explain] = 3, derivation[check] = 2, explainer = 2, everything else 1). Our words only; equations with book
numbers; every number is ours (computed with the WIP `fluidpy`, constants of convention 3).

**Section check (book section → IDs → notebook cells).**
- front matter (title, explainer index, setup) — cells 1–4
- §1.1 Fluid Mechanics → no A item; B/C, recap and skip items C01 — cells 5–10
- §1.2 Units of Measurement → no A item; B/C, recap and skip items R01 R02 R03 — cells 11–22
- §1.3 Solids, Liquids, and Gases → no A item; B/C, recap and skip items C02 N01 C03 C04 N02 — cells 23–37
- §1.4 Continuum Hypothesis → C06; B/C, recap and skip items C05 C07 N03 (also stated here from other sections: C24); derivations D34 D35; explainer E1 — cells 38–88
- §1.5 Molecular Transport Phenomena → C12; B/C, recap and skip items C08 C09 C10 N04 C11 N05 C13 C14 N06; explainer E2 — cells 89–125
- §1.6 Surface Tension → no A item; B/C, recap and skip items C15 N07 C16 N08 N09 — cells 126–133
- §1.7 Fluid Statics → C20; B/C, recap and skip items C17 C18 N10 C19 C21 N11 C22 N12; derivations D05 D37 — cells 134–174
- §1.8 Classical Thermodynamics → C25 C35 C36; B/C, recap and skip items C23 C26 N13 C27 C28 C29 C30 C31 C32 C33 C34 N14 N15 C37 (stated elsewhere: C24 in §1.4); derivations D10 — cells 175–243
- §1.9 Perfect Gas → C40 C45; B/C, recap and skip items C38 C39 C42 C43 N16 C41 C44 C46 C47 C48 (also stated here from other sections: C62 C63); derivations D11 D14; explainer E3 — cells 244–291
- §1.10 Stability of Stratified Fluid Media → C50 C51 C54 C55; B/C, recap and skip items C49 C52 C53 N17 N18 N19 C56 C57 C58 N20 C59 C60 N21 C61 (stated elsewhere: C62 in §1.9, C63 in §1.9); derivations D18 D19 D20 D36; explainer E4 — cells 292–381
- §1.11 Dimensional Analysis → C64 C67 C69; B/C, recap and skip items N22 C65 C66 C68 N23 C70 N24 C71 C72 N25 C73 C74 C75 C76 S01 S02; derivations D28; explainer E5 — cells 382–457
- summary — cell 457 (last); **total ≈ 457 cells**. CORE block cell ranges: C06 39–88 · C12 90–125 · C20 135–174 · C25 176–206 · C35 207–231 · C36 232–243 · C40 245–268 · C45 269–291 · C50 293–307 · C51 308–330 · C54 331–357 · C55 358–381 · C64 383–393 · C67 394–415 · C69 416–457.

**Primer register for this chapter** (all new: `knowledge/primers.md` is empty). 61 📎 primers (P01–P61, each with a
2–4-line demo except where marked *text*), plus 10 glosses (one plain sentence inside a note, no demo; curation §4
"gloss" rows). The term string in each `nb.primer(term, …)` call is **exactly** the Part E "Concept" text before any
parenthesis, so `tools/coverage_check.py` can match them.

### A.0 Front matter
1. `nb.title` — *big idea:* "A fluid is matter that keeps deforming under any shear, described as a **continuum** whose
   density, velocity, pressure and temperature exist at every point. Chapter 1 builds the vocabulary everything else
   rests on: molecules averaged into fields (§1.4), molecules carrying momentum and heat down gradients (§1.5), pressure
   in a fluid at rest (§1.7), the thermodynamics of a fluid particle (§1.8–1.9), when a stratified column is stable
   (§1.10) and how dimensional analysis compresses a problem (§1.11)." *Roadmap* (the 15 A ideas, in order): continuum
   hypothesis (C06) · Newton's viscosity law (C12) · hydrostatic law (C20) · first law (C25) · entropy and Gibbs (C35)
   · speed of sound (C36) · p = ρRT (C40) · isentropic law (C45) · the displaced parcel (C50) · N² (C51) · adiabatic
   lapse rate in both sign conventions (C54) · potential temperature (C55) · dimensional homogeneity (C64) ·
   dimensional matrix (C67) · Buckingham Π (C69). *Prerequisites:* "school physics (force, energy, SI units),
   calculus of one variable, a little linear algebra, Python with numpy — everything else is primed where it is used".
2. `nb.explainer_index` — rows: (`continuum_averaging_volume`, "When does 'density at a point' make sense?", "a value
   at a point is the plateau of an average") · (`viscosity_momentum_diffusion`, "How does the fluid learn that a plate
   moved?", "viscosity diffuses momentum; ν sets the clock") · (`heat_work_paths`, "Same two states: what depends on
   the path?", "q and w depend on the route, Δe and Δs do not") · (`parcel_stability`, "Push a parcel up: does it come
   back?", "stability compares two lapse rates — in both sign conventions") · (`buckingham_pi_machine`, "Why can 7
   pipe variables shrink to 4 numbers?", "Π groups are the null space of a matrix").
3. `nb.setup` — standard setup cell (FAST flag).

### A.1 §1.1 Fluid Mechanics — no A item (placement rule: tagged → C06)
1. `nb.section("1.1", "Fluid Mechanics")` — intro "**What is this section about?** What the subject covers and the three
   ways it makes progress; a map of where this chapter's ideas are used later in the book."
2. `nb.note` — **C01 [B → C06]**: "Fluid mechanics studies liquids and gases at rest and in motion — from blood in a
   capillary to the jet stream. It advances by three routes that check each other: **analysis** (equations solved by
   hand, most of this book), **computation** (Ch. 10) and **experiment** (laboratory and field data). All three need
   the same starting point, the continuum picture developed in the C06 block of §1.4." No equation; the "number" is
   the map itself.
3. `nb.primer[code]("matplotlib figures", …)` — P01: "`fig, ax = plt.subplots()` makes a figure and an axis; `ax.plot`
   draws, `ax.set_xlabel` labels (always with units). Every figure in this notebook follows the same three reading
   notes." *demo:* plot y = x² for x = 0…3 with labelled axes.
4. `nb.figure` — the book map: `from scripts.ch01_book_map import draw_book_map` → `draw_book_map(highlight=("ch13",))`.
   *see:* 16 chapters as nodes with arrows "needs", Chapter 13 (GFD) highlighted. *read:* follow arrows backwards from
   Ch. 13 to see that it rests on Chs. 1, 4, 5, 7, 8, 11, 12. *change:* highlight `ch15` instead and the thermodynamics
   of §1.8–1.9 lights up as its root.

### A.2 §1.2 Units of Measurement — no A item (recaps; tagged → C64 and C40)
1. `nb.section("1.2", "Units of Measurement")` — intro "SI units are pre-book knowledge; we recap them once and let the
   `pint` library carry units through calculations."
2. `nb.recap("R01", "SI base quantities and derived units", …, where="pre-book (school physics); used in C64")` —
   [B → C64]: "Four base quantities are enough for this book: length (m), mass (kg), time (s), temperature (K).
   Everything else is built from them: newton N = kg m s⁻², pascal Pa = N m⁻², joule J = N m, watt W = J s⁻¹, hertz
   Hz = s⁻¹. **A unit is a product of powers of base units** — exactly the idea §1.11 turns into dimensional analysis
   (C64)." Equation $[\mathrm{Pa}] = \mathrm{kg\,m^{-1}\,s^{-2}}$.
3. `nb.primer[code]("pint quantities", …)` — P02: "`Q_(value, unit)` attaches a unit; `.to('Pa')` converts;
   `.dimensionality` shows the base-dimension powers; temperatures in °C are *offset* units (use kelvin inside
   formulas)." *demo:* `Q_(1, "bar").to("Pa")` → 100000 Pa; `(Q_(1000, "kg/m**3") * Q_(9.81, "m/s**2") * Q_(10,
   "m")).to("Pa")` → 98100 Pa.
4. `nb.primer[code]("numpy arrays", …)` — P03: "A numpy array holds many numbers and arithmetic acts on all at once
   (vectorised) — no loops." *demo:* `T_C = np.array([0.0, 15.0, 100.0]); T_C + 273.15` → `[273.15 288.15 373.15]`.
5. `nb.primer[code]("f-strings", …)` — P04: "`f\"{x:.3g} Pa\"` prints a number with 3 significant figures and its
   unit." *demo:* `p = 101325.0; print(f"p = {p/1e3:.1f} kPa")` → `p = 101.3 kPa`.
6. `nb.recap("R02", "SI prefixes", …, where="pre-book; used wherever a number is printed")` — [C → C64]: "Prefixes scale
   by powers of ten (n 10⁻⁹, µ 10⁻⁶, m 10⁻³, k 10³, M 10⁶, G 10⁹); `pint` handles them, so hPa in §1.10, mN/m in §1.6
   and nm in §1.4 need no hand conversion."
7. `nb.recap("R03", "Celsius and kelvin", …, where="pre-book; used in C40")` — [B → C40]: "$T_{[^\circ\mathrm C]} =
   T_{[\mathrm K]} - 273.15$. Laws such as p = ρRT need the **absolute** temperature: 15 °C is 288.15 K, and doubling
   the kelvin temperature (not the Celsius one) doubles p at fixed ρ."
8. `nb.code[explain]` — *code:* `ch01.celsius_to_kelvin(np.array([-40.0, 0.0, 15.0, 20.0]))`,
   `ch01.kelvin_to_celsius(288.15)`, and the pint offset trap `Q_(15, "degC").to("K")`. *expect:* `[233.15 273.15 288.15
   293.15]`, `15.0`, `288.15 kelvin`. *explain:* 1. converts four Celsius readings at once; 2. converts back; 3. pint
   treats °C as an offset unit, so it converts correctly but refuses to multiply it.

### A.3 §1.3 Solids, Liquids, and Gases — no A item (tagged → C12, C20, C06)
1. `nb.section("1.3", "Solids, Liquids, and Gases")` — intro "**What is this section about?** What makes something a
   fluid (it cannot resist shear at rest), how fluids respond to squeezing and pulling, and why gases and liquids differ
   (molecular spacing)."
2. `nb.primer[code]("stress", …)` — P05: "**Stress = force per area** [Pa]. On a surface it splits into a **normal**
   part (pushing or pulling across the surface) and a **shear** part (sliding along it). Shear deforms a small cube into
   a leaning one; the lean angle γ [rad] is the **shear strain**, and dγ/dt [1/s] the **strain rate**." *demo:* 2 N
   sliding force on a 0.01 m² lid → `2/0.01` = 200 Pa shear stress.
3. `nb.note` — **C02 [B → C12]**: "A solid under a small, steady shear stress leans to a fixed angle and stops
   (γ = τ/G, where G [Pa] is its shear modulus — *gloss:* G is how stiff the solid is in shear). A fluid **never stops**:
   any shear stress, however small, makes γ grow for as long as it acts, γ = τt/μ for the simple fluids of this book.
   The rate law behind the second formula is Newton's viscosity law, developed in the C12 block of §1.5."
   Equation $\gamma_{\rm solid} = \tau/G \qquad \gamma_{\rm fluid} = \tau t/\mu$.
4. `nb.primer[code]("np.linspace and np.logspace", …)` — P06: "`np.linspace(a, b, n)` gives n evenly spaced numbers;
   `np.logspace(p, q, n)` gives n numbers from 10ᵖ to 10ᵠ, evenly spaced on a log axis." *demo:* `np.linspace(0, 1,
   5)` → `[0. 0.25 0.5 0.75 1.]`; `np.logspace(-9, -6, 4)` → `[1e-09 1e-08 1e-07 1e-06]`.
5. `nb.figure` — **N02 [C → C12]** sketch: *code:* `t = np.linspace(0, 3, 301)`; `ch01.shear_deformation_history(t,
   10.0, "solid", G=100.0, t_off=1.5)` and `(…, "fluid", mu=1.0, t_off=1.5)`; plot γ(t) for both, shaded loading
   interval 0–1.5 s. *expect:* solid 0.1 rad flat then 0; fluid rises to 15 rad and stays. *see:* the solid (grey) jumps
   to a fixed lean and springs back; the fluid (teal) keeps leaning while loaded and keeps its deformation after.
   *read:* a flat line while loaded = a solid; a rising line = a fluid. *change:* halve μ and the fluid line rises
   twice as steeply; the solid line would not change. (One-sentence pointer in the markdown: "our drawing of the idea of
   the book's Fig. 1.1; Newton's law behind the slope is C12.")
6. `nb.note` — **N01 [C → C12]**: "Some materials sit in between: **plastics** behave as solids until the stress passes
   a yield value; **viscoelastic** materials (paint, egg white) partly spring back. These non-Newtonian fluids are where
   Newton's law fails (Ch. 4 §4.5; blood in Ch. 16 §16.3)."
7. `nb.note` — **C03 [B → C20]**: "Normal stress can **push** (compression) or **pull** (tension). Both solids and
   fluids resist compression, but a liquid pulled hard enough breaks: when its pressure falls below the **vapour
   pressure** (*gloss:* the pressure at which the liquid boils at that temperature, ≈ 2.3 kPa for water at 20 °C) vapour
   cavities form — **cavitation**, seen on ship propellers. On a surface with unit normal n̂ (*gloss:* the dot product
   t·n̂ picks out the part of a vector along n̂) the normal stress is σ_n = t·n̂; pressure (C20) is the compressive normal
   stress of a fluid at rest." Equation $\sigma_n = \mathbf t\cdot\hat{\mathbf n}$. *Number (same cell as C04 below):*
   `ch01.traction_components([3.0, 4.0, 0.0], [1.0, 0.0, 0.0])` → σ_n = 3 Pa (tension), |τ| = 4 Pa.
8. `nb.primer[code]("mole, kilomole, molecular weight and Avogadro's number", …)` — P07: "A kilomole (kmol) is
   6.022×10²⁶ molecules (Avogadro's number per kmol, A_o). The molecular weight M_w is the mass of one kmol in kg (air
   28.96 kg/kmol, water 18.02). So one molecule weighs m = M_w/A_o, and a density ρ holds n = ρA_o/M_w molecules per
   m³. **This book uses kmol, not mol — a factor 1000 hides here.**" *demo:* `28.9644/6.02214076e26` → 4.81e-26 kg;
   `1.225*6.02214076e26/28.9644` → 2.55e25 m⁻³.
9. `nb.note` — **C04 [B → C06]**: "In a liquid the molecules touch their neighbours; in a gas they are about ten times
   farther apart than their size, so a gas is mostly empty space and easy to compress. A liquid poured into a container
   forms a **free surface**; a gas fills all the space it is given. The spacing is n^{-1/3}: water 0.31 nm, sea-level
   air 3.4 nm. How many molecules a sampling box holds is the question of C06." Equation $\ell_{\rm spacing} =
   (n)^{-1/3},\; n = \rho A_o/M_w$.
10. `nb.code[explain]` — *code:* `ch01.traction_components(...)` (C03 number above); `ch01.mean_molecular_spacing(998.0,
    18.015)`, `ch01.mean_molecular_spacing(1.225, ch01.M_W_AIR)`. *expect:* `(3.0, [0, 4, 0], 4.0)`; 3.11e-10 m;
    3.40e-9 m (ratio ≈ 11).

### A.4 §1.4 Continuum Hypothesis — C06 (+C05, C07, N03, C24 · D34, D35 · E1)
1. `nb.section("1.4", "Continuum Hypothesis")` — intro "**What is this section about?** Why we may speak of density,
   velocity and temperature *at a point* although matter is made of molecules, and the test (Kn) that says when this
   fails."
2. `nb.core("C06", "The continuum hypothesis: a value at a point is the plateau of an average", question="Matter is
   mostly empty space between molecules — so what can 'the density at this point' possibly mean?")`
3. `nb.md` — **The problem in plain words:** "A weather model stores one temperature per grid box; a pipe-flow formula
   uses the velocity *at* the wall. Zoom into either and you find molecules flying about with nothing in between. Yet
   every equation from Chapter 2 on treats ρ, **u**, p and T as smooth functions of position. We need to know why that
   works — and when it stops working (microchannels, the upper atmosphere, aerosols)."
4. `nb.md` — **The idea:** "Measure density by counting: put a box of side L around the point, weigh the molecules
   inside, divide by the volume. Grow the box and watch the reading." ASCII sketch:
   ```
   box side L     1 nm      10 nm     100 nm     10 µm        1 mm        10 cm
   molecules N    0.03      25        2.5e4      2.5e10       2.5e16      2.5e22
   relative noise  —        20 %      0.6 %      6e-6         6e-9        —
   reading        0 or huge  noisy    settling   ──── plateau ────       drifts with the flow
                  └ molecular noise ┘            └ continuum window ┘    └ flow's own variation ┘
   ```
   "**The density at a point is the plateau value** — it exists because a wide range of box sizes lies between the
   molecular scale and the scale on which the flow itself changes."
5. `nb.primer[code]("temperature as molecular kinetic energy", …)` — P08: "Temperature measures the mean random
   kinetic energy of molecules: ½m⟨|u|²⟩ = (3/2)k_BT, with Boltzmann's constant k_B = 1.380649×10⁻²³ J/K." *demo:*
   `np.sqrt(3*ch01.K_B*288.15/ch01.molecular_mass(ch01.M_W_AIR))` → 498 m/s (rms speed of air molecules).
6. `nb.primer[code]("Newton's second law and momentum", …)` — P09: "Momentum = mass × velocity [kg m/s]. Newton's second
   law: the net force on a body equals the rate of change of its momentum (F = ma for fixed mass). Draw every force on a
   body (a free-body diagram), add them with signs, set the sum equal to m a." *demo:* a 0.1 kg ball hits a wall at
   2 m/s and bounces back: `0.1*2 - 0.1*(-2)` → momentum change 0.4 kg m/s; delivered in 0.01 s → force `0.4/0.01` =
   40 N.
7. `nb.primer[code]("np.random.default_rng", …)` — P10: "`rng = np.random.default_rng(0)` makes a seeded random
   generator (same numbers every run); `rng.normal(0, s, n)` draws Gaussian numbers, `rng.poisson(lam, n)` counts."
   *demo:* `rng.normal(0, 2, 100000).std()` → ≈ 2.0; `rng.poisson(100, 100000).std()` → ≈ 10.
8. `nb.primer[code]("Gaussian velocity components and mean square speed", …)` — P11: "In a gas at rest each velocity
   component is Gaussian with zero mean and variance k_BT/m. The mean square speed adds the three components:
   ⟨|u|²⟩ = ⟨u_x²⟩ + ⟨u_y²⟩ + ⟨u_z²⟩ = 3⟨u_x²⟩." *demo:* `u = ch01.maxwellian_velocities(100000, 288.15, m_air)`;
   `(u**2).sum(1).mean() / (u[:,0]**2).mean()` → ≈ 3.0.
9. `nb.note` — **C05 [B → C06]**: "**Pressure is momentum delivered by molecular impacts.** Molecules hitting a wall
   bounce back; each bounce pushes the wall. Averaged over the ~10²⁷ impacts per second on every square metre, the push
   per area is steady — that average is the pressure. The derivation below makes this exact." Equation
   $p = \tfrac13\,n\,m\,\langle|\mathbf u|^2\rangle$ (D34).
10. `nb.derivation("D34", "Pressure from molecular impacts", ref="kinetic theory, §1.4 idea", …)` — Part F · D34
    (8 steps).
11. `nb.code[explain]` — *code:* `n_air = ch01.number_density(1.225, ch01.M_W_AIR)`; `m_air =
    ch01.molecular_mass(ch01.M_W_AIR)`; `u = ch01.maxwellian_velocities(200_000 if not FAST else 50_000, 288.15,
    m_air, seed=0)`; `p_kin = ch01.molecular_pressure(n_air, m_air, u)`; `p_wall = ch01.wall_impact_pressure(u, m_air,
    n_air)`; `p_nkT = n_air*ch01.K_B*288.15`. *expect:* p_nkT = 101 327 Pa; p_kin ≈ 1.015×10⁵ and p_wall ≈ 1.03×10⁵
    (sampling scatter < 2 %); printed ratio to p_nkT. *explain:* 1. molecules per m³ and mass of one molecule (P07);
    2. 200 000 random velocities (P10, P11); 3. step 7 of D34 as a formula; 4. step 5 of D34 as a count of wall hits;
    5. step 8 (p = nk_BT).
12. `nb.primer[code]("Poisson counting statistics", …)` — P12: "When many independent things each have a small chance
    to fall in a box, the count N is Poisson: its variance equals its mean N̄, so its standard deviation is √N̄ and the
    **relative** scatter is 1/√N̄." *demo:* `c = rng.poisson(25, 100000); c.std()/c.mean()` → ≈ 0.20 = 1/√25.
13. `nb.primer[code]("power laws and log–log plots", …)` — P13: "y = a·x^k is a straight line of slope k on log–log axes
    (`ax.loglog`), because log y = log a + k log x." *demo:* `x = np.logspace(0, 3, 4); np.polyfit(np.log10(x),
    np.log10(5*x**-1.5), 1)[0]` → −1.5.
14. `nb.md` — **The maths, step by step:** "1. Around a point **x** take a small volume δV [m³] holding mass δm [kg].
   2. Define the continuum density as the ratio, $\rho(\mathbf x) = \delta m/\delta V$ (§1.4). 3. The ratio is useful
   only if δV holds so many molecules that the count's scatter is negligible (lower bound, from D35 below) **and** is so
   small that the flow's density does not change across it (upper bound, set by the flow's length scale L_flow [m]).
   4. Every field variable of the book — **u**, p, T — is defined the same way, as a plateau average." Symbols: n
   [m⁻³] molecules per volume, m [kg] mass of one molecule, N̄ = nδV expected count, ε [–] fractional amplitude of the
   flow's own density variation.
15. `nb.derivation("D35", "How noisy is a density measured in a small box?", ref="§1.4, our addition", …)` — Part F ·
    D35 (5 steps).
16. `nb.worked_example("a 10 µm cube of sea-level air", …)` — "1. n = ρA_o/M_w = 1.225 × 6.022×10²⁶ / 28.96 =
    2.55×10²⁵ m⁻³. 2. δV = (10⁻⁵ m)³ = 10⁻¹⁵ m³. 3. N̄ = nδV = 2.55×10¹⁰ molecules. 4. δm = N̄m = 2.55×10¹⁰ × 4.81×10⁻²⁶
   = 1.225×10⁻¹⁵ kg, so δm/δV = 1.225 kg/m³. 5. Relative noise 1/√N̄ = 1/(1.6×10⁵) = 6.3×10⁻⁶. 6. Same air, a 10 nm
   cube: N̄ = 25, noise = 1/5 = 20 % — no useful density."
17. `nb.primer[code]("tuple unpacking", …)` — P14: "A function may return several values as a tuple; `a, b = f()`
    unpacks them into names." *demo:* `mean, std = (1.225, 0.01); print(mean, std)`.
18. `nb.code[explain]` — *code:* `L = np.logspace(-9, 0, 37 if not FAST else 19)`; `mean, std = ch01.sample_density(L,
    n_air, m_air, n_samples=200 if not FAST else 60, variation=0.2, L_flow=1.0, seed=0)`; `rel = std/mean`;
    `noise = ch01.density_noise_expected(L, n_air)`; `rho_pt = 1.225*(1 + 0.2)` (the density at the point: the box sits on a crest of the variation); `drift =
    np.abs(ch01.box_average_density(L, 1.225, 0.2, 1.0)/rho_pt - 1)`; print a small table at L = 1e-8, 1e-5, 1e-1. *expect:* rel(10 nm) ≈ 0.2, rel(10 µm) ≈ 6×10⁻⁶, drift(10 µm) ≈
    3×10⁻¹¹, drift(10 cm) ≈ 2.7×10⁻³; mean(1 nm…10 µm) ≈ 1.47 kg/m³. *explain:* 1. 37 box sizes from 1 nm to 1 m (P06); 2. 200 random boxes per size,
    counts drawn as in D35 around the box average of a flow whose density varies ±20 % over 1 m; 3. measured relative
    scatter; 4. the D35 prediction; 5. how far the noiseless box average has drifted from the value at the point.
19. `nb.primer[code]("assert np.allclose", …)` — P15: "`assert np.allclose(a, b, rtol=…)` stops the notebook if two
    arrays differ by more than the tolerance — our proof that a hand-written version matches the tested library." *demo:*
    `assert np.allclose([1.0, 2.0], [1.0, 2.0 + 1e-9])`.
20. `nb.check_agree` — **from scratch** (curation §7): place molecules and count them. *code:* `rng =
    np.random.default_rng(1)`; `n_mol = 100_000 if not FAST else 40_000`; `pos = rng.random((n_mol, 3))` (unit cube,
    so n = n_mol per unit volume); for side b in `[0.05, 0.1, 0.2]`: choose 300 (FAST 120) random corners in
    `[0, 1-b]³`, `counts = np.array([np.all((pos >= c) & (pos < c + b), axis=1).sum() for c in corners])`; `rel_mine =
    counts.std(ddof=1)/counts.mean()`. Then `m_lib, s_lib = ch01.sample_density(b, n_mol, 1.0, n_samples=300)`;
    `assert np.allclose(rel_mine, ch01.density_noise_expected(b, n_mol), rtol=0.25)`; `assert np.allclose(rel_mine,
    s_lib/m_lib, rtol=0.25)`. *expect:* rel ≈ 0.28, 0.10, 0.035 (N̄ = 12.5, 100, 800). Markdown after: "Statistical
    agreement (25 %) is the right test here: both sides are random; the deterministic law is D35."
21. `nb.figure` — the noise → plateau → drift figure (`ax.loglog`): teal dots `rel`, grey dashed `noise` (slope −3/2),
    orange `drift`, green band where both < 10⁻³ (between L ≈ 0.3 µm and ≈ 6 cm), vertical markers at the molecular
    spacing (3.4 nm) and at L_flow. *see:* scatter falling on a line of slope −3/2, then an orange curve rising; a
    shaded window about five decades wide. *read:* inside the window any box gives the same density to better than
    0.1 % — that shared number is ρ(**x**). *change:* water (n 1300× larger) moves the grey line down by √1300 ≈ 36, so
    the window starts ≈ 11× smaller; a thin boundary layer (L_flow = 1 mm) pulls the orange curve three decades left.
22. `nb.primer[code]("animate and show_animation", …)` — P16: "`animate(update, frames)` calls `update(i)` once per
    frame to move artists already drawn; `show_animation(anim, player='frames')` shows ◀ ▮▮ ▶ buttons." *demo:* a dot
    moving along a line, 10 frames.
23. `nb.animation[explain]` — zoom animation (`player="frames"`, 37 frames, FAST 19): left panel a 2-D slice of the box
    region — random dots when the expected count in view is < 2000, otherwise a smooth colour field of the ±20 %
    variation; the sampling square in accent purple; right panel the reading trace so far on the axes of item 21 with
    the current point enlarged. *explain:* frame i sets L = 10^(−9 + i/4); dots are redrawn from the seeded generator;
    the trace appends one sample. Followed by `nb.md` with see/read/change: *see* dots → noisy reading → smooth field →
    drift; *read* the moment the dots become too many to draw is roughly where the reading settles; *change* with
    variation = 0 the right end stays flat forever.
24. `nb.primer[code]("slider_figure", …)` — P17: "`slider_figure(fn, name, values)` precomputes curves for every slider
    position so the slider works on the web page without Python." *demo:* curves y = a·x for a ∈ {1, 2, 3}.
25. `nb.plotly[explain]` — **C07/N03 slider**: altitude z = 0…80 km (17 steps); for each z, `T, p, rho =
    ch01.standard_atmosphere(z)`, `l = ch01.mean_free_path_air(T, p)`, traces "Kn = l/L" vs body size L = 1 nm…10 m
    (log x, log y) and two horizontal guides "Kn = 0.01 (continuum limit)" and "Kn = 0.1". *explain:* 1. the USSA-1976
    state at each altitude; 2. Jennings mean free path; 3. Kn for every body size. Reading notes in the markdown: at sea
    level only bodies < 7 µm are outside the continuum; at 80 km (l ≈ 4.4 mm) even a 40 cm probe is.
26. `nb.note` — **C07 [B → C06]**: "The same scale question for *motion*: molecules travel a mean free path l between
    collisions. If the body or channel size L is not much larger, the gas cannot smooth out velocity differences and
    the continuum description fails. The test is the **Knudsen number**: continuum for Kn ≲ 0.01, slip flow for
    0.01–0.1, transition up to ~10, free-molecular beyond." Equation $Kn = l/L$ (§1.4). Number in item 29.
27. `nb.note` — **N03 [B → C06]**: "For air at 300 K and 1 atm, l ≈ 67 nm (Jennings formula with Sutherland viscosity —
    our value, a published order of magnitude). In liquid water molecules are in contact, so l is comparable with the
    0.3 nm molecular size and Kn is negligible for any practical L."
28. `nb.note` — **C24 [B → C06]** (book §1.8, taught here): "A **fluid particle** is the plateau box followed as the
    flow carries it: 1) it always contains the same molecules, 2) it is large enough that its properties are well
    defined (the lower end of the window), 3) it is small enough — and molecular collisions frequent enough — that it
    **relaxes** to equilibrium (settles after any disturbance; the relaxation time is defined with C23 in §1.8) much
    faster than the flow changes it. A 10 µm particle of air holds 2.5×10¹⁰ molecules (noise 6×10⁻⁶) and its molecules
    collide every ≈ 1.4×10⁻¹⁰ s. Ch. 3 turns this into the material volume."
29. `nb.code[explain]` — numbers for C07, N03, C24: `l = ch01.mean_free_path_air(300.0)`; `ch01.knudsen_number(l,
    1e-6)`; `ch01.collision_time(l, 300.0, m_air)`; `ch01.mean_free_path_air(*ch01.standard_atmosphere(80e3)[:2])`.
    *expect:* 6.72×10⁻⁸ m; 0.067; 1.43×10⁻¹⁰ s; 4.4×10⁻³ m.
30. `nb.primer[code]("show_viz", …)` — P18: "`show_viz('ch01', slug)` embeds an interactive explainer: in Jupyter from
    disk, in Colab from the web; on the published page it fills the window." *demo:* none needed beyond the call — the
    primer text only, plus `print(show_viz.__doc__.splitlines()[0])`.
31. `nb.explainer("continuum_averaging_volume", heading="When does 'density at a point' make sense?", why="The
    continuum window is a sweep across nine decades with randomness in it. Dragging the box, resampling and switching
    the medium shows noise, plateau and drift in a way a single static curve cannot.", tries=[…])` — tries: "Press
    ▶ in the Walkthrough and watch the teal samples squeeze into the grey band as the box grows."; "Switch the medium to
    liquid water and find where the reading first stays within 1 %."; "Set L_flow = 1 mm (a thin boundary layer): does a
    continuum window still exist?"; "Choose '80 km air' and read Kn for a 1 cm sensor."
32. `nb.md` — **What would change if…** "…the gas were ten times thinner (40 km up)? n drops 10×, so the noise line rises
    by √10 and l grows 10×: the window shrinks from both ends. **Next:** once properties exist at points they can differ
    from point to point — and molecules carry them across. That is §1.5 (C12)."

### A.5 §1.5 Molecular Transport Phenomena — C12 (+C08, C09, C10, N04, C11, N05, C13, C14, N06 · E2)
1. `nb.section("1.5", "Molecular Transport Phenomena")` — intro "**What is this section about?** Random molecular
   motion carries salt, heat and momentum from where there is more to where there is less. Three look-alike laws
   describe it; the one for momentum defines viscosity."
2. `nb.core("C12", "Newton's law of viscosity: shear stress is a flux of momentum", question="Hold one plate still and
   drag the other: how does the fluid in between find out, and what force does it exert on the plates?")`
3. `nb.md` — **The problem in plain words:** "Stir honey and it resists; stir water and it barely does. Oil in a
   bearing, wind dragging the ocean surface along, the atmospheric boundary layer — all are fluid layers sliding over
   each other. We want the force per area between layers and how fast motion spreads into still fluid."
4. `nb.md` — **The idea:** molecules hop between neighbouring layers, carrying their layer's property:
   ```
   y ↑   fast layer  u + du   ● → ● → ●     a molecule hopping DOWN brings extra x-momentum
         ─ ─ ─ ─ ─ ─ ─ ─ A B ─ ─ ─ ─ ─      a molecule hopping UP brings a deficit
         slow layer  u        ● → ●         net: x-momentum flows from fast to slow = a drag τ on AB
   ```
   Table: | what spreads | its gradient | law | diffusivity [m²/s] | · species Y | ∇Y | Fick (1.1) | κ_m | · heat | ∇T |
   Fourier (1.2) | κ = k/ρC_p | · x-momentum | du/dy | Newton (1.3) | ν = μ/ρ |. "**Viscosity is diffusion of
   momentum.**"
5. `nb.md` — reminder: "In §1.3 (C02) a fluid kept deforming under any shear stress; the law below says how fast."
6. `nb.primer[code]("ordinary derivative as a slope", …)` — P19: "du/dy is the slope of the profile u(y): how much u
   changes per metre of y [1/s for a velocity]." *demo:* `y = np.array([0, 1e-3]); u = np.array([0, 1.0]);
   (u[1]-u[0])/(y[1]-y[0])` → 1000 s⁻¹.
7. `nb.note` — **C08 [B → C12]**: "Where a property (salt, heat, momentum) is unevenly spread, random molecular motion
   carries more of it from high to low than back: a **flux** down the gradient, proportional to the gradient for small
   gradients. A profile that is not straight therefore relaxes — its bumps fill in. For a profile f(y) this is
   the model diffusion equation (ours in Ch. 1; derived in Ch. 4 §4.8 for heat), with D the diffusivity."
   Equation $\partial f/\partial t = D\,\partial^2 f/\partial y^2$.
8. `nb.note` — **C09 [C → C12]**: "The amount of a constituent is measured by its mass fraction Y (kg of it per kg of
   mixture; its partial density is ρY) — salinity S in §1.10 is one; salt and heat diffusing at different rates drive
   double-diffusive instability (Ch. 11 §11.5)."
9. `nb.note` — **C10 [B → C12]**: "**Fick's law**: the mass flux of a constituent points down its mass-fraction gradient
   (*gloss:* ∇Y is the vector of slopes ∂Y/∂x, ∂Y/∂y, ∂Y/∂z, pointing where Y increases fastest; Ch. 2 §2.9). κ_m [m²/s]
   is the mass diffusivity. Number: water vapour in air, κ_m = 2.6×10⁻⁵ m²/s, ρ = 1.20 kg/m³, Y increasing upward by
   0.01 per metre → J = −3.1×10⁻⁷ kg m⁻² s⁻¹ (downward)." Equation $\mathbf J_m = -\rho\,\kappa_m\,\nabla Y$ (1.1).
10. `nb.note` — **N04 [C → C12]**: "Drawn as a profile with a downward flux arrow, this is the 'species' mode of the
    explainer at the end of this block."
11. `nb.note` — **C11 [B → C12]**: "**Fourier's law**: heat flux points down the temperature gradient, with thermal
    conductivity k [W m⁻¹ K⁻¹]. Divided by the heat stored per kelvin per volume, ρC_p, it gives the thermal
    diffusivity κ = k/(ρC_p) (*gloss:* C_p, the heat per kg per kelvin, is defined properly in §1.8). Number: water,
    k = 0.60, temperature falling upward by 100 K/m → q = +60 W/m² upward; κ = 1.44×10⁻⁷ m²/s." Equation
    $\mathbf q = -k\,\nabla T$ (1.2).
12. `nb.primer[code]("boundary conditions", …)` — P20: "What the solution must do at the edges. At a solid wall a
    viscous fluid moves with the wall (**no-slip**). **Steady** means nothing changes in time; **transient** is the
    approach to it." *demo:* `u0 = np.zeros(5); u0[-1] = 1.0; u0` → the initial profile with the top plate moving.
13. `nb.md` — **The maths, step by step:** "1. Layers slide in x with speed u(y) [m/s]; y [m] across the gap.
    2. Across a surface AB of constant y the molecules carry x-momentum per area per time — a momentum flux [kg m⁻¹ s⁻²
    = N/m² = Pa], i.e. a stress. 3. Following (1.1)–(1.2), the flux is proportional to the gradient du/dy; the constant
    μ [Pa s = kg m⁻¹ s⁻¹] is the **dynamic viscosity**: $\tau = \mu\,du/dy$ (1.3), the drag the faster fluid above AB
    exerts on the fluid below. 4. Units check: Pa s × s⁻¹ = Pa ✓. 5. Divide by density to get a diffusivity with the
    same units as κ_m and κ: $\nu \equiv \mu/\rho$ (1.4) [m²/s]. 6. A diffusivity D spreads a disturbance across a
    distance h in a time of order h²/D." Symbols: τ [Pa], μ [Pa s], ν [m²/s], h [m] gap, U [m/s] plate speed.
14. `nb.note` — **C13 [B → C12]**: "μ depends on temperature, in opposite ways: a **gas** gets more viscous when hot
    (faster molecules carry momentum farther; roughly μ ∝ T^{1/2}), a **liquid** gets less viscous (the cohesive forces
    that resist sliding weaken). Numbers (item 22): air 1.79×10⁻⁵ Pa s at 288 K falls to 1.42×10⁻⁵ at 217 K (Sutherland);
    the √T rule gives 1.55×10⁻⁵, 9 % high. Water 1.00×10⁻³ at 20 °C falls to 3.5×10⁻⁴ at 80 °C."
15. `nb.note` — **C14 [B → C12]**: "**Kinematic viscosity** ν = μ/ρ is the diffusivity of momentum. It decides how
    fast motion spreads, so the ranking can flip: air's μ is 55× smaller than water's, but its ν is 15× larger — motion
    diffuses through air 15× faster." Equation $\nu \equiv \mu/\rho$ (1.4).
16. `nb.worked_example("water between plates 1 mm apart", …)` — "1. Top plate U = 1 m/s, bottom fixed, gap h = 1 mm:
    steady slope du/dy = U/h = 1000 s⁻¹. 2. Water μ = 1.0×10⁻³ Pa s → τ = μU/h = 1.0 Pa. 3. On a 10 cm × 10 cm plate:
    F = τA = 1.0 × 0.01 = 0.01 N. 4. ν = μ/ρ = 1.0×10⁻³/998 = 1.0×10⁻⁶ m²/s → h²/ν = 10⁻⁶/10⁻⁶ = 1.0 s. 5. Air:
    μ = 1.81×10⁻⁵ → τ = 0.018 Pa; ν = 1.51×10⁻⁵ → h²/ν = 0.066 s."
17. `nb.primer[code]("Python dictionaries", …)` — P23: "A dict maps names to values: `d['mu']`. fluidpy returns
    several named results this way." *demo:* `fp = {'rho': 998.2, 'mu': 1.0e-3}; fp['mu']/fp['rho']` → 1.0e-6.
18. `nb.code[explain]` — *code:* for name in ("water", "air"): `fp = ch01.fluid_properties(name)`; `tau =
    ch01.newton_shear_stress(fp["mu"], U/h)`; `nu = ch01.kinematic_viscosity(fp["mu"], fp["rho"])`; `t_d =
    ch01.diffusion_time(h, nu)`; print. *expect:* water τ = 1.002 Pa, ν = 1.004×10⁻⁶ m²/s, t = 0.996 s; air τ = 0.0181
    Pa, ν = 1.506×10⁻⁵, t = 0.0664 s. *explain:* 1. properties at 20 °C from cited correlations; 2. Eq. (1.3) for the
    steady Couette slope; 3. Eq. (1.4); 4. h²/ν.
19. `nb.primer[code]("finite differences", …)` — P21: "A derivative from samples: the **central difference**
    (f[i+1] − f[i−1])/(2Δy) is second-order accurate (error ∝ Δy²). Stepping ∂f/∂t = D∂²f/∂y² forward with
    f_i^{new} = f_i + r(f_{i+1} − 2f_i + f_{i−1}), r = DΔt/Δy², is the **FTCS** scheme; it is stable only for r ≤ ½."
    *demo:* central difference of sin at 0.5 with Δ = 0.01 vs cos(0.5): error 4×10⁻⁶.
20. `nb.code[explain]` — the transient: *code:* `N = 101 if not FAST else 51`; `y = np.linspace(0, h, N)`; `dy = y[1] -
    y[0]`; `nu = ch01.kinematic_viscosity(1.0e-3, 998.2)`; `dt = ch01.stable_time_step(nu, dy)`; `nsteps =
    int(1.5*h**2/nu/dt)`; `u0 = np.zeros(N); u0[-1] = U`; `F = ch01.ftcs_diffusion_1d(u0, nu, dy, dt, nsteps, values=(0.0,
    U), save_every=max(nsteps//60, 1))`; `times = …`; `tau_b, tau_t = ch01.wall_shear_history(F, dy, 1.0e-3)`;
    `u_exact = ch01.couette_startup_profile(y, times[-1], U, h, nu)`; `print(np.abs(F[-1] - u_exact).max())`. *expect:*
    ≈ 33 000 steps (N = 101), 61 saved profiles, τ_b and τ_t → 1.00 Pa, max error vs series < 10⁻³ m/s. *explain:* grid,
    stable step (P21), initial profile with the top plate moving (P20), the march, wall stresses at every saved time, the
    analytic series ("derived in Ch. 8; used here only as a check").
21. `nb.primer[code]("np.gradient", …)` — P22: "`np.gradient(f, y, edge_order=2)` returns central differences inside
    and second-order one-sided differences at the ends (the default `edge_order=1` is only first order there)." *demo:*
    `np.gradient(np.array([0., 1., 4., 9.]), 1.0, edge_order=2)` → `[0. 2. 4. 6.]` (exact for y²).
22. `nb.check_agree` — **from scratch** (curation §7): *code:* `u = ch01.couette_startup_profile(y, 0.1*h**2/nu, U, h,
    nu)`; `tau_mine = np.empty_like(u)`; interior loop `tau_mine[i] = mu*(u[i+1] - u[i-1])/(2*dy)`; ends
    `mu*(-3u[0] + 4u[1] - u[2])/(2dy)` and `mu*(3u[-1] - 4u[-2] + u[-3])/(2dy)`; `assert np.allclose(tau_mine,
    ch01.shear_stress_profile(u, y, mu))`; `assert np.allclose(tau_mine, mu*np.gradient(u, y, edge_order=2))`. Also the
    C13 numbers: `ch01.sutherland_viscosity(np.array([288.15, 216.65]))`, `ch01.viscosity_power_law(216.65,
    ch01.sutherland_viscosity(288.15), 288.15)`, `ch01.water_viscosity(np.array([293.15, 353.15]))`, and C10/C11:
    `ch01.fick_mass_flux(1.20, 2.6e-5, np.array([0, 0.01, 0]))`, `ch01.fourier_heat_flux(0.60, np.array([0, -100.0,
    0]))`, `ch01.thermal_diffusivity(0.60, 998.0, 4182.0)`. *expect:* the numbers quoted in the notes.
23. `nb.animation[explain]` — `player="video"`, 60 frames (FAST 30): left u(y, t)/U (teal bold) with the steady line
    (grey dashed) and the series solution (orange thin); right τ(y, t) from `shear_stress_profile` (rose) with the steady
    value μU/h (grey dashed); title "t = 0.12 h²/ν = 0.12 s". Fixed axes. Then `nb.md` **What you see / How to read it /
    What would change if…**: *see* a boundary layer creeping down from the moving plate; the stress spike at the top
    wall decays and spreads until τ is the same everywhere. *read* uniform τ means every layer passes on exactly the
    momentum it receives — steady Couette flow; by t ≈ 0.5 h²/ν the profile is within 1 % of straight. *change* double
    h and everything happens 4× slower in seconds but identically in units of h²/ν.
24. `nb.note` — **N05 [B → C12]**: "The animation is our version of the book's relaxation picture: a sheared profile
    relaxing toward the straight one, with the stress on a surface AB equal to μ times the local slope."
25. `nb.note` — **N06 [C → C12]**: "All three laws are **linear** in the gradient and use **first** derivatives only —
    accurate because molecular steps are tiny compared with the scale of the gradients; Ch. 4 §4.5 generalises (1.3) to
    a tensor law, and the diffusion equations appear in Ch. 4 §4.8 and Ch. 8 §8.4."
26. `nb.explainer("viscosity_momentum_diffusion", heading="How does the fluid learn that a plate moved?", why="The
    phenomenon is a transient: a profile creeping into the gap and settling. Playing, scrubbing and swapping fluids shows
    h²/ν at work and the μ-versus-ν surprise, which no single frame can.", tries=[…])` — tries: "Play with water, then
    with air, and compare the settling times on the τ_w(t) panel."; "Pick honey with a 5 mm gap and predict: faster or
    slower than water?"; "Switch to heat mode and read the Prandtl-number comparison in Explain."; "Click the profile to
    see the τ(y) arithmetic."
27. `nb.md` — **What would change if…** "…the top plate stopped again? Momentum would diffuse out through both walls and
    the profile would decay on the same clock h²/ν. **Next:** when nothing moves at all, every shear stress is zero and
    only the normal stress — pressure — remains: fluid statics (§1.7, C20). First, one more force that acts at
    interfaces (§1.6)."

### A.6 §1.6 Surface Tension — no A item (placement rule: tagged → C20)
1. `nb.section("1.6", "Surface Tension")` — intro "**What is this section about?** An interface between two fluids pulls
   like a stretched sheet; when it is curved, the pressure on its concave side is higher. The C20 block (§1.7) uses this
   for capillary rise."
2. `nb.note` — **C15 [B → C20]**: "**Surface tension σ** is the pull per unit length along any line drawn in the
   interface [N/m] — equivalently the energy needed to make one extra square metre of interface [J/m² = N/m]. It depends
   on the pair of fluids, the temperature and traces of surfactants (soap lowers it). Water against air at 20 °C: σ =
   72.7 mN/m (IAPWS)."
3. `nb.note` — **N07 [B → C20]**: "A spherical drop or bubble of radius R: the pull σ × 2πR around a cut through its
   middle balances the pressure excess acting on the cut's area πR² (stated, not derived). A 1 mm drop of water has
   145 Pa extra inside; a 1 µm fog droplet 1.45×10⁵ Pa — about 1.4 atmospheres." Equation $p_i - p_o = 2\sigma/R$.
4. `nb.note` — **C16 [B → C20]**: "**Laplace's pressure jump** for any gently curved interface uses its two principal
   radii of curvature (*gloss:* at a point on a surface, the circles that best fit it in two perpendicular directions;
   their radii R₁, R₂, taken negative when the centre is on the outer side, as on a saddle). The concave side has the
   higher pressure. Cylinder (R₂ → ∞): σ/R. Sphere (R₁ = R₂): 2σ/R. Symmetric saddle (R₂ = −R₁): no jump." Equation
   $p_i - p_o = \sigma\left(\dfrac1{R_1} + \dfrac1{R_2}\right)$ (1.5).
5. `nb.code[explain]` — *code:* `sigma = ch01.surface_tension_water(293.15)`; `ch01.laplace_pressure_jump(sigma,
   np.array([1e-3, 1e-6]))`; `ch01.laplace_pressure_jump(sigma, 1e-3, np.inf)`; `ch01.laplace_pressure_jump(sigma, 1e-3,
   -1e-3)`. *expect:* 0.07274 N/m; [145.5, 1.455e5] Pa; 72.7 Pa; 0.0 Pa.
6. `nb.note` — **N08 [C → C20]**: "Narrow tubes where these jumps matter are **capillary** tubes and the phenomena
   **capillarity**; surface tension becomes a boundary condition in Ch. 4 §4.10 and drives capillary waves in Ch. 7
   §7.3."
7. `nb.note` — **N09 [C → C20]**: "The book sketches the hemisphere force balance and a patch with its two principal
   radii; those radii return with capillary waves (Ch. 7 §7.3)."

### A.7 §1.7 Fluid Statics — C20 (+C17, C18, N10, C19, C21, N11, C22, N12; uses C03, C15, C16 · D05, D37)
1. `nb.section("1.7", "Fluid Statics")` — intro "**What is this section about?** Pressure in a fluid at rest: it acts
   equally in all directions, does not change sideways, and grows downward by exactly the weight of the fluid above."
2. `nb.core("C20", "The hydrostatic law: pressure grows with depth by the weight above", question="A cube of lake water
   sits perfectly still. What must the pressure do between its top and bottom faces so that nothing moves?")`
3. `nb.md` — **The problem in plain words:** "Your ears hurt at the bottom of a pool; a dam is thicker at its base; a
   barometer weighs the atmosphere; the deep ocean crushes submarines. The base state of every stratified flow, wave and
   weather system in this book is this one balance."
4. `nb.md` — **The idea:** ASCII cube:
   ```
            (p + dp)·dx dy   ↓        top face pushed down
          ┌────────────┐
          │  weight    │  ρ g dx dy dz ↓
          └────────────┘
            p·dx dy          ↑        bottom face pushed up
   at rest:  (bottom push) − (top push) = weight   ⇒   pressure must fall going up
   ```
5. `nb.figure` — **N11 [B]** drawing: `from scripts.ch01_drawings import draw_cube_forces`; `draw_cube_forces(ax)`.
   *see:* the fluid cube with the two vertical pressure forces and z upward. *read:* the arrows' difference is what holds
   the cube's weight. *change:* a denser fluid needs a larger difference for the same cube.
6. `nb.note` — **N11 [B → C20]**: "The pressure difference between a fluid element's top and bottom faces balances its
   weight — the picture above, made exact in D05."
7. `nb.primer[code]("weight and gravitational acceleration", …)` — P24: "Weight W = mg points down; g = 9.80665 m/s²
   (standard). For fluid of density ρ in volume V: W = ρVg." *demo:* 1 L of water `1000*1e-3*9.80665` → 9.81 N.
   One-line reminder: "Newton's second law was primed in C06 (P09)."
8. `nb.primer[code]("partial derivative", …)` — P25: "For p(x, y, z), ∂p/∂x is the slope in x with y and z held fixed.
   If p depends on z only, ∂p/∂z = dp/dz." *demo:* `p = lambda x, z: 1e5 - 9807*z + 0*x`; `(p(0.01, 0) - p(-0.01,
   0))/0.02` → 0.0; `(p(0, 0.01) - p(0, -0.01))/0.02` → −9807.
9. `nb.primer[code]("first-order Taylor expansion", …)` — P26: "Close to a point, a smooth function is its value plus
   slope × step: p(z + dz) ≈ p(z) + (dp/dz)dz. The neglected terms shrink like dz², so after dividing by dz they vanish
   as dz → 0 (**orders of smallness**)." *demo:* `np.exp(0.1)` vs `1 + 0.1` → error 0.005 ≈ 0.1²/2.
10. `nb.note` — **C17 [B → C20]**: "In a fluid at rest the normal stress is the **pressure** (the compressive normal stress
    of C03, §1.3). **Absolute** pressure is measured from vacuum; **gauge** pressure from the local atmosphere, p_gauge =
    p − p_atm. Standard atmosphere 101 325 Pa; 1 bar = 10⁵ Pa. A tyre at 2.0 bar gauge holds 3.01 bar absolute."
    Equation $p_{\rm gauge} = p - p_{\rm atm}$.
11. `nb.note` — **C18 [B → C20]**: "**Pressure is the same in every direction** at a point (stated, not derived): take a
    tiny wedge of fluid at rest. The forces on its faces scale with face area (∝ size²) but its weight with volume
    (∝ size³), so as the wedge shrinks the weight becomes negligible and the face pressures must match for any
    orientation. With a 1 mm tall wedge of water the difference is only ½ρg dz = 4.9 Pa, and it vanishes with dz
    (`ch01.wedge_pressure_difference(1000.0, 1e-3, 0.6)`). So pressure is a scalar." Equation $p_1 = p_2 = p_3$ (1.6).
12. `nb.note` — **N10 [C → C20]**: "This direction-independent normal stress becomes the −pδᵢⱼ part of the stress tensor
    (Ch. 2 §2.6, Ch. 4 §4.5)."
13. `nb.note` — **C19 [B → C20]**: "**Pascal's law** (stated): on a fluid cube at rest the two faces normal to x feel
    p(x) and p(x + dx); nothing else acts sideways (gravity is vertical), so they must be equal and ∂p/∂x = 0; likewise
    in y. Points at the same height in the same connected, resting fluid have the same pressure — why water levels in
    connected vessels line up." Equation $\partial p/\partial x = \partial p/\partial y = 0$ (1.7).
14. `nb.derivation("D05", "The hydrostatic law", ref="1.8", …)` — Part F · D05 (8 steps).
15. `nb.primer[code]("definite integral", …)` — P27: "∫ₐᵇ f dz adds up f dz from a to b; for a constant, ∫ₐᵇ c dz =
    c(b − a)." *demo:* `np.trapezoid(np.full(11, 9807.0), np.linspace(-10, 0, 11))` → 98 070.
16. `nb.note` — **C21 [B → C20]**: "For uniform density, integrating (1.8) from z = 0 (pressure p₀) to z gives a straight
    line (stated): pressure grows by ρgh at depth h = −z. Ten metres of water add 9.81×10⁴ Pa ≈ one atmosphere. And
    because pressure pushes harder on the bottom of a submerged body than on its top, the body feels an upward net force —
    **buoyancy**, derived next." Equation $p = p_0 - \rho g z$ (1.9).
17. `nb.primer[code]("net force from pressure", …)` — P28: "Pressure p pushes on each bit of surface dA along the inward
    normal, force −p n̂ dA. The net force on a closed body is the sum over its whole surface, F = −∮ p n̂ dA — for a box
    just six faces." *demo:* `ch01.net_pressure_force_on_box(lambda x, y, z: 1e5 + 0*z, (0, 1, 0, 1, 0, 1))` → `[0. 0.
    0.]` (uniform pressure pushes equally from all sides).
18. `nb.derivation("D37", "Buoyancy: the net pressure force on a submerged body", ref="from 1.9", …)` — Part F · D37
    (7 steps).
19. `nb.note` — **C03 reminder** (stated in §1.3): "A liquid cannot be pulled below its vapour pressure: a sealed water
    column hanging from a pump cannot be lifted more than ≈ 10 m, because by (1.9) the pressure at its top would have to
    fall below ≈ 2.3 kPa."
20. `nb.note` — **C22 [B → C20]** (Ex. 1.1, stated): "In a thin tube dipped into water the curved meniscus (C15, C16 in
    §1.6) lowers the pressure just under it, so water rises until its weight balances the surface-tension pull around the
    rim: σ·2πR·sin α = ρgh·πR². **The book's α is measured from the horizontal** — the complement of the usual contact
    angle θ_c (α = 90° − θ_c; α = 90° is perfect wetting; `ch01.alpha_from_contact_angle`). Clean glass, R = 1 mm:
    h = 14.9 mm; R = 0.1 mm: 149 mm." Equation $h = \dfrac{2\sigma\sin\alpha}{\rho g R}$.
21. `nb.note` — **N12 [C → C20]**: "Along the axis the pressure falls linearly from atmospheric at the outside level to
    p_atm − ρgh just under the meniscus (drawn in item 27); pressure jumps at free surfaces return in Ch. 7 §7.3."
22. `nb.worked_example("10 m of water and a 1-litre block", …)` — "(g ≈ 9.81) 1. p(−10 m) = 101 325 + 1000 × 9.81 × 10 =
    199 425 Pa absolute; gauge 98 100 Pa. 2. 1-litre block: buoyancy ρ_water g V = 1000 × 9.81 × 0.001 = 9.81 N,
    whatever it is made of. 3. Steel (7850 kg/m³) weighs 77.0 N → net 67.2 N down. 4. Pine (600 kg/m³) weighs 5.89 N →
    net 3.92 N up: it floats. 5. Oil (800 kg/m³, 1 m) over 9 m of water: p = 101 325 + 800×9.81×1 + 1000×9.81×9 =
    197 463 Pa (with g = 9.80665: 197 430 Pa)."
23. `nb.primer[code]("functions as arguments and lambda", …)` — P29: "Python functions are values: pass one to another
    function. `lambda z, p: 1000.0` is a one-line function (here: density that ignores z and p)." *demo:* `apply = lambda
    f, x: f(x); apply(lambda z: 2*z, 3.0)` → 6.0.
24. `nb.primer[code]("explicit stepping", …)` — P30: "To follow dy/dx = f(x, y), take small steps: y_{k+1} = y_k +
    f(x_k, y_k)Δx (**Euler**; error ∝ Δx). For oscillators, update the velocity first and use the new velocity for the
    position (**Euler–Cromer**; keeps energy bounded)." *demo:* Euler for dy/dx = −y from 1 to x = 1 with 1000 steps →
    0.3677 vs e⁻¹ = 0.3679.
25. `nb.primer("scipy.integrate.solve_ivp", …)` — P31 *text + one line*: "`integrate_hydrostatic` hands dp/dz = −ρ(z, p)g
    to scipy's adaptive Runge–Kutta solver `solve_ivp`, which chooses its own steps to reach a tolerance (here 10⁻¹⁰)."
    *demo:* `from scipy.integrate import solve_ivp; solve_ivp(lambda x, y: -y, (0, 1), [1.0], rtol=1e-10).y[0, -1]` →
    0.367879.
26. `nb.code[explain]` — *code:* `z = np.linspace(0, -10, 201)`; `p_w = ch01.hydrostatic_pressure_uniform(z, ch01.P_ATM,
    1000.0)`; `p_layers = ch01.layered_pressure(z, [1.0, 9.0], [800.0, 1000.0])`; lake with salt stratification
    `rho_fn = lambda z, p: 1000.0 - 0.5*z` (z < 0, so 0.5 kg/m³ heavier per metre down); `p_strat =
    ch01.integrate_hydrostatic(z, rho_fn, ch01.P_ATM)`; closed form `p_exact = ch01.P_ATM + ch01.G0*(1000.0*(-z) +
    0.25*z**2)`; `assert np.allclose(p_strat, p_exact, rtol=1e-10)`; buoyancy `F = ch01.net_pressure_force_on_box(lambda
    x, y, zz: ch01.hydrostatic_pressure_uniform(zz, ch01.P_ATM, 1000.0), (0, 0.1, 0, 0.1, -1.1, -1.0))`;
    `ch01.buoyancy_force(1000.0, 1e-3)`; `ch01.gauge_pressure(p_w[-1])`; capillary `ch01.capillary_rise(
    ch01.surface_tension_water(293.15), np.pi/2, 998.2, np.array([1e-3, 1e-4]))`. *expect:* p_w(−10) = 199 391.5 Pa;
    p_layers(−10) = 197 430 Pa; p_strat(−10) = 199 636.7 Pa; F = [0, 0, 9.807] N = buoyancy_force; gauge 98 066.5 Pa;
    h = [0.0149, 0.149] m. *explain:* numbered, one item per call.
27. `nb.check_agree` — **from scratch** (curation §7): Euler march p[k+1] = p[k] − ρ(z_k)gΔz downward with Δz = −1 mm
    on the stratified lake: loop over 10 000 steps; `assert np.allclose(p_mine[::50], p_strat, rtol=1e-6)`. Markdown:
    "Euler's error is ½ρ′gΔz per metre — about 0.02 Pa here, far inside 10⁻⁶."
28. `nb.figure` — p(z) for three still waters: uniform (teal), oil over water (orange, kink at z = −1 m where the slope
    changes from 800g to 1000g), stratified lake (blue, slightly curving); top x-axis in gauge kPa; inset: the capillary
    drawing `draw_capillary(ax_in, R=1.0, alpha=np.radians(90), h=3.0)` with the N12 pressure line. *see:* straight lines
    whose slopes are −ρg. *read:* slope = −ρg, so the kink marks the density jump; gauge and absolute differ by a constant.
    *change:* mercury (13 600 kg/m³) would make the line 13.6× steeper — why barometers use it.
29. `nb.plotly[explain]` — slider over the upper-layer density 600…1000 kg/m³ (25 steps) in a 2 m oil-over-8 m-water
    tank: traces "p(z) [kPa gauge]" (`layered_pressure`) and "block faces" (markers: gauge pressure on the top and bottom
    faces of a 10 cm cube centred 1.5 m deep, inside the upper layer) and "buoyancy × 10³/A" hover. *explain:* one
    `layered_pressure` call per slider position; the markers are p at z = −1.45 and −1.55 m; their difference × area
    equals ρ_top g V (D37). Reading note: the face-pressure gap grows with the upper density, and so does the buoyancy.
30. `nb.md` — **What would change if…** "…the fluid were air over 1 km? Its density falls with height (it is
    compressible), so (1.8) still holds but p(z) is no longer a straight line — that needs the equation of state (§1.9,
    C40). **Next:** a fluid particle also has temperature and energy; §1.8 is its thermodynamics (C25)."

### A.8 §1.8 Classical Thermodynamics — C25, C35, C36
1. `nb.section("1.8", "Classical Thermodynamics")` — intro "**What is this section about?** The energy bookkeeping of a
   fluid particle: what heat and work do to it (first law), which quantities describe its state, entropy and the Gibbs
   relations, and two material properties the rest of the book needs — the speed of sound and the thermal expansion
   coefficient. (The fluid particle itself, C24, was already met in §1.4.)"

#### Block C25
2. `nb.core("C25", "The first law: heat in plus work done changes the internal energy", question="Take a kilogram of air
   from one state to another by two different routes. What is the same at the end, and what depends on the route?")`
3. `nb.md` — **The problem in plain words:** "A bicycle pump gets warm; an engine turns heat into work cycle after cycle;
   air rising over a mountain cools although nobody removes heat. To follow a fluid particle's temperature we need a
   bookkeeping rule for energy — Chapter 4 turns it into the energy equation."
4. `nb.md` — **The idea:** a bank account: "the balance (internal energy e) depends only on where you are; deposits
   arrive through two channels (heat q, work w), and how much came through each channel depends on the route."
   ASCII p–v sketch of two routes A (isothermal curve) and B (down then across) between the same states 1 and 2.
5. `nb.note` — **C23 [B → C25]**: "A **thermodynamic system** is a fixed amount of matter exchanging heat and work, but
   no mass, with its surroundings. After a disturbance it **relaxes** back to **equilibrium**, where p, T and the rest
   are well defined; the time this takes is the **relaxation time** — in air about the time between molecular
   collisions, ≈ 10⁻¹⁰ s, which is why the fluid particle of C24 can be treated as always in equilibrium."
6. `nb.primer[code]("internal and kinetic energy per unit mass", …)` — P32: "Internal energy e [J/kg] is the energy of
   the molecules' random motion (and their interactions); the particle's bulk motion adds kinetic energy u²/2 [J/kg].
   Thermodynamics here uses only e." *demo:* `0.5*10.0**2` → 50 J/kg for a 10 m/s breeze vs `717.6*288.15` ≈ 2.07×10⁵
   J/kg stored as e (perfect-gas working model, P33).
7. `nb.primer[code]("perfect-gas law as a working model", …)` — P33: "Until §1.9 derives it, we use air as a perfect gas
   so the examples have numbers: p = ρRT (equivalently pv = RT) with R = 287 J kg⁻¹ K⁻¹, and e = C_vT with C_v = 718 J
   kg⁻¹ K⁻¹ (C40, C41 show why)." *demo:* `v1 = ch01.R_AIR*300/1e5` → 0.861 m³/kg.
8. `nb.primer[code]("differentials", …)` — P34: "dx is a tiny change of a quantity that has a value in every state (an
   **exact** differential: adding the pieces gives x₂ − x₁ whatever the route). δq and δw are tiny amounts *transferred*
   along a route: they are not changes of anything (**inexact**), so their sums depend on the route." *demo:* sum of
   dT over two different lists of steps from 300 to 450 K → both 150.
9. `nb.md` — **The maths, step by step:** "1. Per unit mass, energy is conserved: heat added δq plus work done **on** the
   particle δw equals the rise of internal energy, $\delta q + \delta w = \Delta e$ (1.10) [J/kg]. 2. Sign convention:
   work done *on* the gas is positive (compression), so an expanding gas has w < 0. 3. For a slow, frictionless
   (reversible) change the only work is the boundary moving: δw = −p dv (stated as C27 below). 4. Along a route from 1
   to 2: w = −∫p dv, the area under the path on a p–v diagram, with sign. 5. Then q = Δe − w."
10. `nb.note` — **N13 [C → C25]**: "v = 1/ρ [m³/kg] is the specific volume; it appears in (1.11)–(1.18) and in the
    energy equation of Ch. 4 §4.8."
11. `nb.note` — **C26 [B → C25]**: "**Heat and work are path functions** — energy in transit across the boundary; the
    internal energy is a **state function** — a property with one value per state, so Δe = e₂ − e₁ ignores the route.
    A **reversible** process is one carried out so slowly and without friction that the system is always in
    equilibrium; only then is the path a line on a p–v diagram. Number: routes A, B, C below all have Δe = 0 but w =
    −59.7, −43.1 and −86.1 kJ/kg."
12. `nb.note` — **C27 [B → C25]**: "For a reversible change with only boundary work, the work done on unit mass is −p dv
    (a piston face of area A moved outward by dx does work pA dx on the surroundings, and A dx is the volume increase),
    so (stated) the first law reads" Equation $de = dq - p\,dv$ (1.11). "Friction work is excluded."
13. `nb.note` — **C28 [B → C25]**: "An **equation of state** ties state functions together; for a simple one-component
    fluid **two** properties fix the state, e.g. p = p(v, T) (thermal) or e = e(p, T) (caloric). Seawater needs a third,
    salinity (§1.10). Air at v = 0.861 m³/kg and 300 K has p = 1.00×10⁵ Pa (`ch01.perfect_gas_state(rho=1/0.861,
    T=300.0)`)." Equation $p = p(v, T)\ \text{or}\ e = e(p, T)$ (1.12).
14. `nb.primer[code]("line integral along a path", …)` — P35: "∫p dv along a route adds p × (tiny volume change) step by
    step along that route; on a p–v diagram it is the area under the route (negative if v decreases)." *demo:* for a
    horizontal route p = 5×10⁴ Pa from v = 0.861 to 1.722: `5e4*(1.722-0.861)` → 43 050 J/kg.
15. `nb.primer[code]("natural logarithm and exponential", …)` — P36: "ln x undoes eˣ; ∫dv/v = ln(v₂/v₁); ln 2 = 0.693.
    e⁻¹ ≈ 0.368 — a quantity that shrinks by that factor has decayed one **e-folding**." *demo:* `np.log(2),
    np.exp(-1)`.
16. `nb.worked_example("doubling the volume of 1 kg of air at 300 K, three ways", …)` — "State 1: T = 300 K, p = 10⁵ Pa
    → v₁ = RT/p = 287.06 × 300/10⁵ = 0.861 m³/kg. State 2: v₂ = 1.722 m³/kg, T = 300 K. Δe = C_v(300 − 300) = 0 on every
    route. **A (isothermal):** p = RT/v, w = −RT ln(v₂/v₁) = −287.06 × 300 × 0.693 = −59.7 kJ/kg; q = Δe − w = +59.7.
    **B (cool at fixed v to 150 K, then expand at fixed p):** leg 1 w = 0, q = C_vΔT = 717.6 × (−150) = −107.6 kJ/kg;
    leg 2 p = R·150/v₁ = 5.0×10⁴ Pa, w = −p(v₂ − v₁) = −43.1 kJ/kg, Δe = +107.6, q = 150.7; totals w = −43.1, q = +43.1.
    **C (expand at fixed p to 600 K, then cool at fixed v):** w = −p₁(v₂ − v₁) = −86.1, q = +86.1 kJ/kg."
17. `nb.primer[code]("trapezoid rule", …)` — P37: "Approximates ∫f dx by adding trapezoids between samples; error ∝ Δx².
    `np.trapezoid(f, x)` returns the total, `scipy.integrate.cumulative_trapezoid` the running sum." *demo:*
    `np.trapezoid(1/np.linspace(1, 2, 101), np.linspace(1, 2, 101))` → 0.693153 vs ln 2 = 0.693147.
18. `nb.code[explain]` — *code:* `v1 = ch01.R_AIR*300.0/1e5`; `for kind in ("isothermal", "isochoric-isobaric",
    "isobaric-isochoric"): tot = ch01.path_heat_work_totals(kind, v1, 300.0, 2*v1, 300.0)`; print q, w, Δe in kJ/kg;
    `path = ch01.process_path("isochoric-isobaric", (v1, 300.0), (2*v1, 300.0))`; `run = ch01.process_heat_work(path["v"],
    path["T"])`; `run["w"][-1]`. *expect:* (q, w, Δe) = (59.69, −59.69, 0), (43.06, −43.06, 0), (86.12, −86.12, 0)
    kJ/kg; run["w"][-1] = −43.06 kJ/kg (trapezoid, exact here because the legs are straight on p–v). *explain:*
    1. the start state; 2. exact leg-by-leg totals; 3. the sampled route with its corner; 4. running sums along it.
19. `nb.check_agree` — **from scratch** (curation §7): *code:* for each route's samples `v, T`: `p = ch01.R_AIR*T/v`;
    loop `w -= 0.5*(p[k] + p[k+1])*(v[k+1] - v[k])`; `de = ch01.CV_AIR*(T[-1] - T[0])`; `q = de - w`; `assert
    np.allclose([q, w, de], [tot["q"], tot["w"], tot["de"]], rtol=1e-4, atol=1e-6)` and against `run["w"][-1]`.
    *expect:* passes for the three routes (isothermal trapezoid with 401 samples: relative error ≈ 10⁻⁶).
20. `nb.figure` — p–v diagram: faint isotherms 150, 300, 600 K; routes A (blue), B (orange), C (rose) with arrows; the
    work area under route B's isobaric leg shaded blue; state dots; legend "A: w = −59.7 kJ/kg" etc. (`ax.fill_between`
    explained in the code comment). *see:* three routes, one start, one end, three different areas. *read:* area under a
    route = work done by the gas (−w); since Δe = 0, the heat taken in equals that area. *change:* reverse route C (go
    2 → 1) and every w and q flips sign — the loop A-forward + C-back encloses an area: a heat engine's net work.
21. `nb.animation[explain]` — `player="frames"` (40 frames per route, FAST 20): a piston–cylinder above the p–v plot;
    route A then route B; the piston moves with v, the gas colour follows T, a heater/cooler glyph shows the sign of
    δq, running counters q, w, Δe. Then `nb.md` see/read/change: *see* both routes end with the same piston position and
    colour; *read* the counters for q and w end differently, Δe returns to zero on both; *change* a route through 600 K
    needs more heat, because the expansion happens at higher pressure.
22. `nb.md` — **What would change if…** "…we divided each little δq by the temperature at which it arrives? The next
    block shows that this sum is the same on every route — a new state function, entropy (C35)."

#### Block C35
23. `nb.core("C35", "Entropy and the Gibbs relations: a heat-like quantity that is a state function", question="Heat
    depends on the route. Can we divide it by something so that the result no longer does?")`
24. `nb.md` — **The problem in plain words:** "Why does a hot cup cool but never warm itself? Why can the potential
    temperature of rising air stay constant (§1.10)? Both need a property that tracks heat but belongs to the state."
25. `nb.md` — **The idea:** table for the three routes of C25: | route | q [kJ/kg] | ∫dq/T [J kg⁻¹ K⁻¹] | · A | 59.7 |
    199.0 | · B | 43.1 | 199.0 | · C | 86.1 | 199.0 |. "**Heat depends on the route; heat divided by temperature,
    summed along a reversible route, does not.** That sum is the entropy change."
26. `nb.primer[code]("partial derivative with a variable held fixed", …)` — P39: "In thermodynamics every property
    depends on two others, so a partial derivative must say which one is held fixed: $(\partial h/\partial T)_p$ is the
    rate of change of h with T **at constant p**. Changing what is held fixed changes the answer." *demo:* `h = lambda T,
    p: 1004.7*T; ch01.partial_derivative(h, "T", {"T": 300.0, "p": 1e5})` → 1004.7.
27. `nb.note` — **C29 [B → C35]**: "**Enthalpy** h = e + pv [J/kg]: internal energy plus the 'flow work' pv needed to
    push the particle into place; it is the natural energy for constant-pressure processes and flows. Air at 300 K:
    h − e = RT = 86.1 kJ/kg." Equation $h \equiv e + pv$ (1.13).
28. `nb.note` — **C30, C31 [B → C35]**: "**Specific heats**: how fast h and e grow with temperature — C_p at constant
    pressure, C_v at constant volume [J kg⁻¹ K⁻¹]. They are properties (defined from other properties, not from heat).
    Air: C_p = 1004.7, C_v = 717.6." Equation $C_p \equiv (\partial h/\partial T)_p,\quad C_v \equiv (\partial
    e/\partial T)_v$ (1.14, 1.15).
29. `nb.note` — **C32 [C → C35]**: "Only along a reversible constant-volume (constant-pressure) route with nothing but
    p dv work is the heat per kelvin equal to C_v (C_p); stirring heats a gas with no heat at all — used when h is
    written C_pT in Ch. 4 §4.8 and in stagnation enthalpy, Ch. 15 §15.4."
30. `nb.note` — **C33 [B → C35]**: "**Second law (i): entropy exists.** There is a state function s [J kg⁻¹ K⁻¹] whose
    change is the sum of dq/T along **any reversible** route between the two states. Routes A, B, C give 199.0 J kg⁻¹
    K⁻¹ = R ln 2." Equation $s_2 - s_1 = \displaystyle\int_1^2 \frac{dq_{\rm rev}}{T}$ (1.16).
31. `nb.note` — **N15 [B → C35]**: "In differential form, for a reversible change:" Equation $T\,ds = dq$ (1.17).
32. `nb.primer[code]("product rule for differentials", …)` — P38: "d(pv) = p dv + v dp — both factors may change; forgetting
    one term is the classic slip." *demo:* p, v from (1e5, 0.861) to (1.001e5, 0.8615): exact Δ(pv) vs p Δv + v Δp →
    differ only by Δp·Δv = 0.05 J/kg (tiny).
33. `nb.derivation("D10", "The Gibbs relations", ref="1.18", …)` — Part F · D10 (7 steps).
34. `nb.note` — **C34 [B → C35]**: "**Second law (ii), Clausius–Duhem:** for *any* process the entropy change is at
    least the actual heat received divided by the temperature at which it arrives — equality only when reversible. (We
    write the actual heat δq; with dq_rev the right side would just equal Δs by (1.16).) Free expansion of air into
    vacuum to twice the volume: q = w = Δe = 0 but Δs = R ln 2 = 199 J kg⁻¹ K⁻¹ > 0. Stirring 300 → 310 K at fixed v:
    δq = 0, Δs = C_v ln(310/300) = 23.5 J kg⁻¹ K⁻¹ > 0." Equation $s_2 - s_1 \ge \displaystyle\int_1^2 \frac{\delta
    q}{T}$.
35. `nb.note` — **N14 [C → C35]**: "**Second law (iii):** μ and k must be positive, otherwise momentum or heat would
    un-mix by itself; so viscous dissipation is never negative (Ch. 4 §4.8)."
36. `nb.worked_example("heating air from 300 K to 600 K at 1 bar", …)` — "Constant p, so v doubles too. First Gibbs form
    integrated (perfect gas): Δs = C_v ln(T₂/T₁) + R ln(v₂/v₁) = 717.6 × 0.693 + 287.1 × 0.693 = 497.4 + 199.0 = 696.4
    J kg⁻¹ K⁻¹. Second form: Δs = C_p ln(T₂/T₁) − R ln(p₂/p₁) = 1004.7 × 0.693 − 0 = 696.4 J kg⁻¹ K⁻¹. Same number,
    two routes through the algebra."
37. `nb.code[explain]` — *code:* `ch01.perfect_gas_entropy_change(300.0, v1, 600.0, 2*v1)`;
    `ch01.perfect_gas_entropy_change_p(300.0, 1e5, 600.0, 1e5)`; along the sampled isobaric route
    `path = ch01.process_path("isobaric", (v1, 300.0), (2*v1, 600.0))`, `run = ch01.process_heat_work(path["v"],
    path["T"])`, `ch01.entropy_change_reversible(run["q"], path["T"])`; `ch01.irreversible_process("free_expansion",
    300.0, v1, v2=2*v1)`; `ch01.irreversible_process("stirring", 300.0, v1, T2=310.0)`. *expect:* 696.41, 696.41, 696.4
    (trapezoid), {q: 0, w: 0, de: 0, ds: 198.97}, {q: 0, w: 7176, de: 7176, ds: 23.53}.
38. `nb.primer[code]("sympy", …)` — P40: "sympy does algebra with symbols: `sp.symbols`, `sp.diff`, `sp.simplify`; a
    result that simplifies to 0 is an identity." *demo:* `x = sp.symbols('x'); sp.simplify(sp.diff(x**3, x) - 3*x**2)` → 0.
39. `nb.code[explain]` — sympy residual of both Gibbs forms for the perfect gas (curation §2b): `T, v, cv, R =
    sp.symbols('T v c_v R', positive=True)`; `s = cv*sp.log(T) + R*sp.log(v)`; `e = cv*T`; `p = R*T/v`; `h = e + p*v`;
    for a path T(t), v(t): residual `T*ds - de - p*dv` and `T*ds - dh + v*dp` with `dp` from p(T, v) → both `0`.
40. `nb.figure` — T–s diagram: faint isochores T ∝ e^{(s−s₁)/C_v} and isobars T ∝ e^{(s−s₁)/C_p} through state 1; the
    three C25 routes drawn; the area under route B shaded orange = its heat. *see:* the same three routes as on the p–v
    diagram, now ending at the same s. *read:* on a T–s diagram heat is the area under the route (∫T ds), so different
    areas = different q, one Δs. *change:* an adiabatic reversible route would be a vertical line (ds = 0) — the
    isentrope of C45.
41. `nb.md` — **What would change if…** "…the process were not reversible (friction)? The Gibbs relations still give Δs
    because they only link state functions; only T ds = dq fails. **Next:** how stiff a fluid is when squeezed at
    constant entropy sets the speed of sound (C36)."

#### Block C36
42. `nb.core("C36", "The speed of sound: stiffness at constant entropy", question="Why does a sound pulse cross 340 m of
    air in one second but about 1480 m of water?")`
43. `nb.md` — **The problem in plain words:** "Thunder arrives after lightning; sonar measures ocean depth; weather
    models filter sound out. More importantly for this book: comparing a flow's speed with c decides whether density
    changes can be ignored (Ch. 4 §4.9) — the incompressible approximation — and c is the start of Ch. 15."
44. `nb.md` — **The idea:** "A pressure push squeezes the next layer, which pushes the next… The signal speed grows with
    **stiffness** (pressure change per density change) and falls with inertia. The squeeze is too fast for heat to flow,
    so the stiffness is taken at constant entropy." Table: squeeze by Δp = 1 kPa → air Δρ/ρ ≈ 0.7 %, water ≈ 4.5×10⁻⁷.
45. `nb.md` — **The maths, step by step:** "1. Speed of sound (proof in Ch. 15 §15.2): $c^2 = (\partial p/\partial
    \rho)_s$ (1.19) [m²/s²]. 2. Define the bulk modulus K ≡ ρ(∂p/∂ρ)_s [Pa] (how many pascals per fractional density
    change), so c² = K/ρ. 3. For air squeezed isentropically p = p₀(ρ/ρ₀)^γ (derived in C45): ∂p/∂ρ = γp₀ρ^{γ−1}/ρ₀^γ =
    γp/ρ, so K = γp. 4. For water K ≈ 2.2 GPa (modified Tait model, `tait_pressure`). 5. **Incompressible limit:** if
    ∂ρ/∂p → 0 then ∂p/∂ρ → ∞ and c → ∞ — pressure changes are felt everywhere at once." Primer reminders in the text:
    P39 (held fixed), P21 (central difference).
46. `nb.worked_example("air and water", …)` — "Air at 101 325 Pa, 1.225 kg/m³, γ = 1.4: c² = γp/ρ = 1.4 × 101 325/1.225
    = 115 800 m²/s² → c = 340.3 m/s. Water: c = √(K/ρ) = √(2.2×10⁹/1000) = 1483 m/s. To compress water by 1 % needs
    Δp = 0.01K = 22 MPa — the pressure 2.2 km down in the ocean."
47. `nb.code[explain]` — *code:* `ch01.sound_speed_from_eos(lambda rho, s: ch01.isentropic_pressure(rho, ch01.P_ATM,
    1.225), 1.225)`; `ch01.sound_speed_from_eos(lambda rho, s: ch01.tait_pressure(rho), 1000.0)`;
    `ch01.perfect_gas_sound_speed(288.15)`. *expect:* 340.29, 1483.2, 340.30 m/s. *explain:* the EOS is passed as a
    function (P29); the library differentiates it at fixed s.
48. `nb.check_agree` — **from scratch** (curation §7): `drho = 1e-6*rho`; `c_mine = np.sqrt((P(rho + drho) - P(rho -
    drho))/(2*drho))` for both EOS; `assert np.allclose(c_mine, [c_air_lib, c_water_lib], rtol=1e-6)`; `assert
    np.isclose(c_mine[0], ch01.perfect_gas_sound_speed(288.15), rtol=1e-4)`.
49. `nb.figure` — c vs bulk modulus K (log–log, K = 10⁴…10¹¹ Pa): lines c = √(K/ρ) for ρ = 1.2 and 1000 kg/m³; air
    (K = γp = 1.42×10⁵ Pa) and water (2.2×10⁹) marked; an arrow "K → ∞: incompressible, c → ∞". *see:* two parallel
    lines of slope ½. *read:* for the same density, 100× stiffer means 10× faster. *change:* heating air raises p at
    fixed ρ, hence K and c (c ∝ √T, C47).
50. `nb.note` — **C37 [B → C36]**: "The **thermal expansion coefficient** α [1/K] is the same kind of held-fixed
    derivative, at constant pressure: the fractional density drop per kelvin. Water is strange: α(20 °C) = +2.07×10⁻⁴
    K⁻¹ but α(2 °C) = −3.3×10⁻⁵ K⁻¹ — below ≈ 4 °C it contracts on warming, which is why lakes freeze from the top. α
    returns in the adiabatic lapse rate (C54) and the Boussinesq approximation (Ch. 4)." Equation $\alpha \equiv
    -\dfrac1\rho\left(\dfrac{\partial\rho}{\partial T}\right)_p$ (1.20). *Number cell:* `ch01.thermal_expansion_coefficient(
    lambda T, p: ch01.water_density(T), np.array([275.15, 293.15]), 1e5)` (loop if not vectorised) → −3.27e-5, 2.07e-4.
51. `nb.md` — **What would change if…** "…the fluid were a perfect gas? Then K = γp = γρRT and c = √(γRT) depends only
    on temperature — §1.9 makes this and α = 1/T explicit. **Next:** the perfect-gas law itself (C40)."

### A.9 §1.9 Perfect Gas — C40 (+R03, C38, C39, C41, C48, C62, C63 · D11), C45 (+C42, C43, N16, C44, C46, C47 · D14 · E3)
1. `nb.section("1.9", "Perfect Gas")` — intro "**What is this section about?** The equation of state of air and most
   gases at ordinary conditions, and what it implies: R from molecular weight, C_p − C_v = R, the isentropic law, c and
   α of a gas. (The isothermal atmosphere and scale height of book §1.10 are stated here, in C40, because they are
   p = ρRT plus hydrostatics.)"

#### Block C40
2. `nb.core("C40", "The perfect-gas law in continuum form: p = ρRT", question="The molecular gas law counts molecules
   in a box. How does it become a law for the density at a point?")`
3. `nb.md` — **The problem in plain words:** "Aircraft performance, weather-balloon soundings, air density for drag,
   altimeters, every atmospheric model: all need density from pressure and temperature. The molecular law is written for
   a container; we need it for a fluid particle."
4. `nb.md` — **The idea:** the chain
   ```
   p = (n/V) k_B T  →  multiply and divide by m  →  ρ (k_B/m) T  →  m = M_w/A_o  →  ρ (R_u/M_w) T  =  ρ R T
        molecules per volume          continuum density (C06)        per kilomole         one gas constant per gas
   ```
5. `nb.md` — R03 reminder: "T is always in kelvin here (recap in §1.2)."
6. `nb.note` — **C38 [B → C40]**: "For n non-interacting molecules in a volume V at absolute temperature T — the result
   D34 (C06) derived from impacts:" Equation $pV = nk_BT$ (1.21). "It holds when attractions between molecules are
   negligible and each molecule's own volume is tiny compared with V/n."
7. `nb.note` — **C39 [B → C40]**: "Constants (CODATA, exact since 2019): k_B = 1.380 649×10⁻²³ J/K; Avogadro's number per
   kilomole A_o = 6.022 140 76×10²⁶ kmol⁻¹; universal gas constant R_u = k_BA_o = 8314.46 J kmol⁻¹ K⁻¹; dry air M_w =
   28.9644 kg/kmol (USSA-1976), so R = R_u/M_w = 287.06 J kg⁻¹ K⁻¹. (Kilomoles: primer P07.)"
8. `nb.derivation("D11", "From the molecular gas law to p = ρRT", ref="1.22", …)` — Part F · D11 (6 steps).
9. `nb.worked_example("air density on three days", …)` — "Sea level, 15 °C: ρ = p/(RT) = 101 325/(287.06 × 288.15) =
   1.225 kg/m³. Hot day, 35 °C: 101 325/(287.06 × 308.15) = 1.145 kg/m³ (−6.5 %). High plateau, 84 kPa and 15 °C:
   84 000/(287.06 × 288.15) = 1.016 kg/m³ — aircraft need longer runways."
10. `nb.code[explain]` — *code:* `R = ch01.gas_constant(ch01.M_W_AIR)`; `ch01.perfect_gas_density(np.array([101325.0,
    101325.0, 84000.0]), np.array([288.15, 308.15, 288.15]))`; `ch01.perfect_gas_state(p=101325.0, T=288.15)`;
    `ch01.perfect_gas_pressure(1.225, 288.15)`. *expect:* 287.058; [1.2250, 1.1455, 1.0155]; (101325, 1.2250, 288.15);
    101 327 Pa.
11. `nb.check_agree` — **from scratch** (curation §7): `rho_mine = p/((ch01.K_B*ch01.N_A_KMOL/ch01.M_W_AIR)*T)`; `assert
    np.allclose(rho_mine, ch01.perfect_gas_density(p, T))`; and the molecular route `n = p/(ch01.K_B*T)`, `rho2 =
    n*ch01.molecular_mass(ch01.M_W_AIR)`; `assert np.allclose(rho2, rho_mine)`.
12. `nb.figure` — isobars ρ(T) for p = 50, 70, 101.325 kPa, T = 200…320 K, with the three worked-example states marked.
    *see:* hyperbola-like falling curves. *read:* at fixed p, ρ ∝ 1/T; at fixed T, ρ ∝ p (vertical spacing). *change:*
    humid air has smaller M_w (water 18 < 28.96), larger R and so lower density at the same p and T.
13. `nb.primer[code]("plotly 3-D surface", …)` — P41: "`go.Figure(go.Surface(x=…, y=…, z=…))` draws an interactive
    surface you can rotate; the page keeps it interactive without Python." *demo:* a 10 × 10 paraboloid.
14. `nb.plotly[explain]` — **C28 illustrated:** surface p(v, T) = RT/v over v = 0.5…2 m³/kg, T = 200…350 K, with a red
    dot at (0.861, 300 K, 10⁵ Pa) and the C25 routes drawn as 3-D lines on the surface (height 520). *explain:* "two
    properties fix the state" means every state is one point on this surface; routes are curves on it.
15. `nb.note` — **C41 [B → C40]**: "For a perfect gas the internal energy and enthalpy depend on **temperature only**,
    e = e(T), h = h(T) — and conversely only a perfect gas has this property (stated; the proof uses the Gibbs relations
    with p = ρRT to show (∂e/∂v)_T = 0). A real gas with attractions does not: for van der Waals CO₂, e = C_vT − a/v, so
    squeezing from v = 0.1 to 0.01 m³/kg at fixed 300 K lowers e by 16.9 kJ/kg (`ch01.van_der_waals_internal_energy`)."
    Equation $e = e(T),\ h = h(T)$.
16. `nb.note` — **C48 [B → C40]**: "For a perfect gas at fixed p, ρ = p/(RT) ∝ 1/T, so (stated) the expansion
    coefficient is simply 1/T: 3.3×10⁻³ K⁻¹ for air at 300 K — sixteen times water's at 20 °C." Equation $\alpha = 1/T$
    (1.28).
17. `nb.note` — **C62 [B → C40]** (book §1.10): "If a layer of atmosphere had one temperature T, hydrostatics (1.8) with
    ρ = p/(RT) gives dp/p = −(g/RT)dz, and (stated) pressure falls exponentially with height." Equation $p(z) =
    p_0\,e^{-gz/RT}$. "At 5 km with T = 250 K: p = 51.2 kPa; the USSA-1976 value is 54.0 kPa (5 % higher, because the
    real lower atmosphere is warmer than 250 K)."
18. `nb.note` — **C63 [B → C40]** (book §1.10): "The height over which p falls by the factor e is the **scale height**
    H = RT/g: 287.06 × 250/9.807 = 7.32 km — a good one-number thickness of the atmosphere (P36 e-folding)." Equation $H =
    RT/g$.
19. `nb.plotly[explain]` — slider T = 200…300 K (21 steps): traces "isothermal p(z)" (`isothermal_pressure(z, P_ATM,
    T)`), "USSA-1976 p(z)" (`standard_atmosphere(z)[1]`, fixed), "p₀/e at z = H" (marker at `scale_height(T)`) over z =
    0…40 km, p in kPa (log y optional). *explain:* the model is one line per T; H moves 0.29 km per 10 K. The same cell
    first prints: `ch01.scale_height(250.0)`, `ch01.isothermal_pressure(5000.0, ch01.P_ATM, 250.0)`,
    `ch01.standard_atmosphere(5000.0)`, `ch01.perfect_gas_expansion_coefficient(300.0)`, `ch01.van_der_waals_internal_energy(
    300.0, 0.1, ch01.VDW_CO2["a"], 655.0) - ch01.van_der_waals_internal_energy(300.0, 0.01, ch01.VDW_CO2["a"], 655.0)`
    → 7318 m; 51 166 Pa; (255.65 K, 54 020 Pa, 0.736); 3.33e-3; 16 914 J/kg.
20. `nb.md` — **What would change if…** "…the gas were squeezed too quickly for heat to escape? Its temperature would rise
    and p would climb faster than ρ. **Next:** the isentropic law (C45)."

#### Block C45
21. `nb.core("C45", "The isentropic law: p/ρ^γ = const", question="Squeeze a gas quickly — no time for heat to leak and no
    friction. How does its pressure follow its density?")`
22. `nb.md` — **The problem in plain words:** "A bicycle pump's barrel heats up; a diesel engine ignites fuel by
    compression alone; a rising air parcel cools as it expands (§1.10); sound is a rapid squeeze. All follow the same
    power law."
23. `nb.md` — **The idea:** table: | squeeze to 2× density | isothermal (heat leaks out) | isentropic (no heat, no
    friction) | · pressure ratio | 2 | 2^1.4 = 2.64 | · temperature | unchanged | ×2^0.4 = 1.32 |. "**The work done stays
    in the gas as heat-like motion, so the pressure climbs faster.** On log–log axes: slope 1 vs slope γ."
24. `nb.note` — **C42 [B → C45]**: "For a perfect gas h = e + pv = e + RT; with e and h functions of T only (C41) the
    partial derivatives in (1.14)–(1.15) become ordinary ones (*gloss:* a total derivative d/dT, because nothing else
    varies), and differentiating gives (stated) the gap between the specific heats: 1004.7 − 717.6 = 287.1 J kg⁻¹ K⁻¹."
    Equation $R = C_p - C_v$ (1.23).
25. `nb.note` — **C43 [B → C45]**: "The ratio of specific heats:" Equation $\gamma \equiv C_p/C_v$ (1.24). "Kinetic theory
    gives 5/3 for monatomic gases (He, Ar), 7/5 for diatomic (N₂, O₂, air) and about 4/3 for triatomic (CO₂, H₂O)
    (`ch01.GAMMA_BY_ATOMICITY`)."
26. `nb.note` — **N16 [B → C45]**: "Air at ordinary temperatures: γ = 1.40 and C_p ≈ 1005 J kg⁻¹ K⁻¹ (our constant 1004.7
    from γ = 1.4 and R). Both C_p and C_v rise slowly with temperature; we treat them as constant, which is what (1.25)
    needs."
27. `nb.note` — **C44 [B → C45]**: "**Adiabatic** = no heat crosses the boundary. **Isentropic** = adiabatic *and*
    frictionless, so s stays constant. Stirring an insulated gas is adiabatic but not isentropic (Δs = 23.5 J kg⁻¹ K⁻¹ in
    the C34 example)."
28. `nb.primer[code]("separation of variables", …)` — P42: "For dy/y = k dx, put each variable on its own side and
    integrate both: ln y = kx + C, so y = y₀e^{kx}." *demo:* `np.log(2.0**1.4)` vs `1.4*np.log(2.0)` → both 0.9704.
29. `nb.primer[code]("exponent rules", …)` — P43: "a^m a^n = a^{m+n}; (a^m)^n = a^{mn}; a^{−m} = 1/a^m; (a/b)^m =
    a^m/b^m." *demo:* `2.0**0.4 * 2.0**1.0`, `2.0**1.4` → 2.639 both.
30. `nb.derivation("D14", "The isentropic law for a perfect gas", ref="1.25", …)` — Part F · D14 (10 steps).
31. `nb.worked_example("a bicycle pump from 1 bar to 2 bar", …)` — "Start 288.15 K. Density ratio (1.25): ρ₂/ρ₁ =
    (p₂/p₁)^{1/γ} = 2^{0.714} = 1.641. Temperature via p = ρRT: T₂/T₁ = (p₂/p₁)(ρ₁/ρ₂) = 2/1.641 = 1.219 → T₂ = 351.3 K,
    63 K hotter. Isothermal comparison: ρ ratio 2, no warming."
32. `nb.code[explain]` — *code:* `ch01.isentropic_ratios(2.0)`; `288.15*ch01.isentropic_ratios(2.0)[0]`; `rho1 =
    ch01.perfect_gas_density(1e5, 288.15)`; `ch01.isentropic_pressure(1.641*rho1, 1e5, rho1)`;
    `ch01.cv_from_cp(ch01.CP_AIR)`, `ch01.gamma_from_cp(ch01.CP_AIR)`; `ch01.perfect_gas_sound_speed(np.array([250.0,
    288.15, 300.0]))`. *expect:* (1.2190, 1.6407); 351.26 K; ≈ 2.000×10⁵ Pa; 717.64; 1.4000; [317.0, 340.3, 347.2] m/s.
33. `nb.check_agree` — **from scratch** (curation §7): integrate dp/dρ = γp/ρ with 20 000 Euler steps (P30) from (ρ₀,
    p₀) to 1.641ρ₀; `assert np.allclose(p_mine, ch01.isentropic_pressure(rho_grid, 1e5, rho1), rtol=1e-4)`.
34. `nb.figure` — log–log p/p₀ vs ρ/ρ₀ from 0.5 to 3: isentrope (slope 1.4, accent) and isotherm (slope 1, grey), the
    pump states marked, slope triangles annotated. *see:* two straight lines through one point. *read:* the slope is the
    exponent: 1.4 for no heat, 1 for heat leaking freely. *change:* helium (γ = 5/3) is steeper still — it heats more
    when pumped.
35. `nb.note` — **C46 [B → C45]**: "Combining (1.25) with p = ρRT (stated) gives temperature and density along an
    isentrope directly from the pressure ratio — the pump numbers above; θ in §1.10 is built from the first." Equation
    $T/T_0 = (p/p_0)^{(\gamma-1)/\gamma},\quad \rho/\rho_0 = (p/p_0)^{1/\gamma}$ (1.26).
36. `nb.note` — **C47 [B → C45]**: "Differentiating p = Kρ^γ at constant s gives (∂p/∂ρ)_s = γp/ρ = γRT, so (stated,
    closing C36) the speed of sound depends on temperature only: 340.3 m/s at 288.15 K, 317 m/s at 250 K. **Trap:**
    Newton used the isothermal slope RT and got √(RT) = 287.6 m/s, 15 % too slow." Equation $c = \sqrt{\gamma RT}$
    (1.27).
37. `nb.explainer("heat_work_paths", heading="Same two states: what depends on the path?", why="Path dependence clicks
    only when you drag a route between the same two states and watch q and w change while Δe and Δs stay fixed; the
    irreversible modes then show where T ds = δq breaks. It ties together C25, C35 and this block.", tries=[…])` —
    tries: "Run 'isothermal ×2', then 'cool then expand', and compare the q and w bars."; "Open the Derivation tab and
    step through D10 with the stirring mode on at step 7."; "Choose 'pump (isentropic)' and check that the bold curve
    lies on the faint isentrope."; "Free expansion: which bar stays at zero, and what does Δs do?"
38. `nb.md` — **What would change if…** "…a parcel of air were carried upward instead of pumped? Its pressure falls, so
    it expands and cools along exactly this isentrope — while the air around it follows hydrostatics. Whether it ends up
    heavier or lighter than its surroundings is §1.10."

### A.10 §1.10 Stability of Stratified Fluid Media — C50 (D18), C51 (+C52, C60, N21, C61), C54 (+C49, C53, N17 · D19), C55 (+C56, N19, C57, N18, C58, N20, C59 · D20, D36 · E4)
1. `nb.section("1.10", "Stability of Stratified Fluid Media")` — intro "**What is this section about?** Whether a fluid
   whose density changes with height is stable: displace a parcel and see whether buoyancy pushes it back. For air this
   becomes a comparison of two lapse rates — written with **two opposite sign conventions** in the literature, both used
   here — and it is captured in one line by the potential temperature. (This book section also contains the isothermal
   atmosphere and scale height, already stated as C62 and C63 in §1.9.)"

**Lapse-rate conventions used in every cell of this section** (convention 1): compute with Kundu's Γ ≡ dT/dz; print the
meteorology value Γ_met = −dT/dz next to it; state stability with the inequality in the form being used.

#### Block C50
2. `nb.core("C50", "The displaced parcel: buoyancy as a restoring or runaway force", question="Nudge a small blob of
   fluid upward and let go. What force acts on it, and does it come back?")`
3. `nb.md` — **The problem in plain words:** "Morning fog stays trapped in a valley; a chimney plume spreads flat under
   an inversion; the ocean's thermocline blocks mixing between warm surface water and the cold deep; on a hot afternoon
   air parcels shoot up into cumulus towers. The same test decides all of them."
4. `nb.md` — **The idea:** ASCII:
   ```
   z_o + ζ ── parcel density  ρ(z_o) + (dρ_a/dz)·ζ      environment  ρ(z_o) + (dρ/dz)·ζ
      ↑ ζ      the parcel changes along ITS OWN path     the surroundings keep THEIR profile
   z_o ────── parcel = environment (at rest)
   parcel heavier than its new surroundings  → pushed back down (stable)
   parcel lighter                            → pushed further up (unstable)
   ```
5. `nb.md` — reminders: "Newton's second law (P09, C06), first-order Taylor (P26, C20) and buoyancy = weight of the
   displaced fluid (D37, C20) are the tools."
6. `nb.primer[code]("linear second-order ODE", …)` — P44: "ζ″ + N²ζ = 0 with constant N². Try ζ = e^{λt}: λ² + N² = 0.
   If N² > 0, λ = ±iN and ζ oscillates like cos Nt; if N² < 0, λ = ±√(−N²) and ζ grows like cosh; if N² = 0, ζ moves in a
   straight line." *demo:* check numerically that ζ = 5 cos(0.01t) satisfies the equation: `t = 100.0; -5*0.01**2*np.cos(
   0.01*t) + 0.01**2*5*np.cos(0.01*t)` → 0.0.
7. `nb.derivation("D18", "The parcel equation and N²", ref="§1.10 parcel equation, 1.29", …)` — Part F · D18 (12 steps).
8. `nb.worked_example("a parcel in the ocean thermocline", …)` — "ρ = 1025 kg/m³; density falls upward by 1 kg/m³ per
   100 m, dρ/dz = −0.01 kg m⁻⁴; treat the parcel as incompressible, dρ_a/dz = 0. N² = −(g/ρ)(dρ/dz − dρ_a/dz) = −(9.81/
   1025)(−0.01 − 0) = 9.57×10⁻⁵ s⁻². N = 9.78×10⁻³ s⁻¹; period 2π/N = 642 s ≈ 10.7 min. Released 5 m up: ζ(t) = 5 cos(
   0.00978 t) m."
9. `nb.code[explain]` — *code:* `N2 = ch01.brunt_vaisala_sq(1025.0, -0.01, 0.0)`; `t = np.linspace(0, 1800, 601)`; `zeta =
   ch01.parcel_displacement(t, 5.0, N2)`; `zeta_nl = ch01.parcel_ode_from_gradients(t, 5.0, 1025.0, -0.01, 0.0)`;
   `ch01.stability_timescale(N2)`; `np.abs(zeta - zeta_nl).max()`. *expect:* 9.567×10⁻⁵; ('period', 642.4); max
   difference < 0.01 m (the nonlinear version keeps ρ_p in the denominator exactly). *explain:* 1. Eq. (1.29);
   2. the D18 solution; 3. the unlinearised Newton law solved by `solve_ivp` (P31); 4. period.
10. `nb.check_agree` — **from scratch** (curation §7): Euler–Cromer (P30) with dt = 1 s: `w += -N2*z*dt; z += w*dt` for
    1800 steps; `assert np.allclose(z_mine[::3], zeta, atol=0.05)` (1 % of the amplitude).
11. `nb.animation[explain]` — `player="video"`, 60 frames (FAST 30) over 600 s: three columns side by side, coloured by
    density — stable N² = +10⁻⁴ s⁻² (teal parcel), neutral 0 (amber), unstable −10⁻⁴ (rose; stops at the column top) —
    and below them ζ(t) so far for each, with the linear solutions as faint ghosts; released 20 m up from rest
    (`parcel_displacement`, unstable branch capped with `np.minimum`). Then `nb.md` see/read/change: *see* one parcel
    bobbing (period 10.5 min), one resting where released, one accelerating away. *read* the curve shapes are cos, a
    flat line and cosh; nothing but the sign of N² differs. *change* doubling N² shortens the period by √2 (to 7.4 min).
12. `nb.md` — **What would change if…** "…the parcel felt friction or exchanged heat? The oscillation would die out
    (the model has neither, so it rings forever). **Next:** everything hinges on the number N² — read it (C51)."

#### Block C51
13. `nb.core("C51", "The Brunt–Väisälä frequency N²: one number that classifies a column", question="What do the sign and
    size of N² tell you about a fluid column before you push anything?")`
14. `nb.md` — **The problem in plain words:** "Internal waves in the ocean, lee waves behind mountains, the Richardson
    number of Ch. 11, the convection switch in climate models — all read N². It is the natural frequency of a
    stratified fluid."
15. `nb.md` — **The idea:** table: | N² | motion after release | time scale | · > 0 | oscillation cos Nt | period 2π/N |
    · = 0 | stays (or drifts at its initial speed) | none | · < 0 | runaway cosh √(−N²)t | e-folding 1/√(−N²) |.
16. `nb.primer[code]("square root of a negative number", …)` — P45: "√(−a) = i√a. An imaginary 'frequency' iσ turns
    cos(iσt) into cosh(σt): growth at rate σ." *demo:* `np.cos(1j*2.0).real, np.cosh(2.0)` → 3.762 both.
17. `nb.md` — **The maths, step by step:** "1. From D18: $N^2 = -\dfrac{g}{\rho(z_o)}\left(\dfrac{d\rho}{dz} -
    \dfrac{d\rho_a}{dz}\right)$ (1.29) [1/s²]; g [m/s²], ρ [kg/m³], gradients [kg/m⁴]. 2. N² > 0 ⇔ dρ/dz < dρ_a/dz: the
    environment's density falls with height **faster** than a displaced parcel's does. 3. N = √N² [rad/s], period T_N =
    2π/N [s]. 4. N² < 0: growth rate σ = √(−N²) [1/s], e-folding time 1/σ; after 3 e-folds the displacement is ×10."
18. `nb.note` — **C52 [B → C51]**: "**Stable** (N² > 0): a displaced parcel returns and oscillates. **Unstable** (N² < 0):
    it accelerates away. **Neutral** (N² = 0): it stays where it is put — reached either when density does not vary at
    all (dρ/dz = dρ_a/dz = 0) or when the environment's density falls exactly as fast as a parcel's (dρ/dz = dρ_a/dz), an
    atmosphere of uniform entropy."
19. `nb.worked_example("reading three values of N²", …)` — "N² = 10⁻⁴ s⁻² → N = 0.01 s⁻¹ → period 628 s = 10.5 min.
    N² = 4×10⁻⁴ → period 314 s. N² = −10⁻⁴ → σ = 0.01 s⁻¹ → e-folds every 100 s; after 5 minutes (3 e-folds) cosh(3) =
    10.1 times the release height."
20. `nb.primer[code]("np.where and np.select", …)` — P46: "`np.where(cond, a, b)` picks elementwise; `np.select([c1, c2],
    [a1, a2], default)` handles several regimes — how labels like 'stable' are assigned to arrays." *demo:*
    `np.select([x > 0, x < 0], ["stable", "unstable"], "neutral")` for `x = np.array([1e-4, 0, -1e-4])`.
21. `nb.code[explain]` — *code:* `N2 = np.array([4e-4, 1e-4, 0.0, -1e-4])`; `ch01.classify_stability(N2)`;
    `[ch01.stability_timescale(x) for x in N2]`. *expect:* ['stable' 'stable' 'neutral' 'unstable']; [('period', 314.2),
    ('period', 628.3), ('none', inf), ('efold', 100.0)].
22. `nb.figure` — ζ(t)/ζ₀ for N² ∈ {4×10⁻⁴, 10⁻⁴, 0, −2.5×10⁻⁵, −10⁻⁴} s⁻² over 0–900 s (`parcel_displacement`), y
    clipped to ±5, teal for stable, amber neutral, rose unstable; dashed verticals at each period, dots at each e-folding
    time. *see:* two cosines of different periods, a flat line, two exponentials leaving the frame. *read:* the sign of N²
    picks the shape, its size the time scale. *change:* halve |N²| and every time scale grows by √2.
23. `nb.primer("live widgets", …)` — P47 *text*: "`live(fn, a=(lo, hi, step))` builds ipywidgets sliders that re-run
    `fn` — only while Python is running (Colab or local); the web page shows a note instead."
24. `nb.live[explain]` — sliders N² (−2×10⁻⁴…4×10⁻⁴) and ζ₀ (1…200 m): plots `parcel_displacement` and
    `parcel_ode_from_gradients` (ρ = 1025, dρ/dz = −N²ρ/g, dρ_a/dz = 0) together; the gap appears only for large ζ₀ with
    N² < 0 (cap reached). *explain:* 1. the slider callback; 2. the two solutions; 3. why they differ only for runaway.
25. `nb.note` — **C60 [B → C51]**: "Seawater density depends on temperature, pressure **and salinity** S (grams of salt
    per kilogram, ≈ 35 g/kg on average). A displaced parcel keeps its salt, so potential density is defined at constant
    S and oceanographers refer densities to sea-level pressure. Linear model: 10 °C, 35 g/kg → 1027.0 kg/m³; 20 °C →
    1025.3 kg/m³ (`ch01.seawater_density_linear`)." Equation $\rho = \rho(T, p, S)$.
26. `nb.note` — **N21 [B → C51]**: "A parcel moved vertically in the ocean is squeezed or relaxed by the hydrostatic
    pressure of its surroundings at constant entropy and salinity, so (stated) its density changes at a rate set by the
    speed of sound: with c = 1500 m/s and ρ = 1025 kg/m³, dρ_a/dz = −4.47×10⁻³ kg m⁻⁴." Equation $\dfrac{d\rho_a}{dz}
    \cong -\dfrac{\rho g}{c^2}$.
27. `nb.note` — **C61 [B → C51]**: "So (1.29) turns into the ocean's stability test: the column is stable when the
    compressibility-corrected density gradient is negative. **Taught as 'has the same sign as' dρ_θ/dz** (analysis §9:
    the book's equality holds only up to a positive factor near the reference pressure). Thermocline with dρ/dz = −0.01:
    −0.01 + 0.00447 = −5.53×10⁻³ kg m⁻⁴ < 0 → stable, but N² drops from 9.57×10⁻⁵ to 5.29×10⁻⁵ s⁻² (period 10.7 → 14.4
    min) because the parcel also compresses as it sinks." Equation $\dfrac{d\rho_\theta}{dz} \;\sim\; \dfrac{d\rho}{dz} +
    \dfrac{\rho g}{c^2}$ (1.35). *Number cell:* `ch01.seawater_density_linear(np.array([283.15, 293.15]), 35.0)`;
    `ch01.isentropic_density_gradient(1025.0, 1500.0)`; `ch01.ocean_potential_density_gradient(-0.01, 1025.0, 1500.0)`;
    `ch01.brunt_vaisala_sq(1025.0, -0.01, ch01.isentropic_density_gradient(1025.0, 1500.0))`.
28. `nb.code` — the number cell of item 27 (counts as the block's second code cell). *expect:* [1027.0, 1025.28];
    −4.467e-3; −5.533e-3; 5.293e-5.
29. `nb.md` — **What would change if…** "…the fluid were air, whose parcels expand a lot as they rise? Then dρ_a/dz is
    large and it is easier to compare temperatures than densities. **Next:** the adiabatic lapse rate (C54)."

#### Block C54
30. `nb.core("C54", "The adiabatic lapse rate — and the two ways its sign is written", question="Air gets colder as you go
    up. Does that make the atmosphere unstable?")`
31. `nb.md` — **The problem in plain words:** "Mountain tops are cold, yet the air over them is usually stable; smog
    sits under inversions; thunderstorms build on hot afternoons. In climate science the **lapse-rate feedback** — how
    the vertical temperature profile changes as the planet warms — is one of the main feedbacks, and those papers write
    the lapse rate as a *positive* number for cooling with height. You need both sign habits."
32. `nb.md` — **The idea:** ASCII:
   ```
   z ↑      environment T(z): cools at 6.5 K/km
     |   \   parcel lifted dry: cools at 9.8 K/km (it expands as pressure drops)
     |    \ \
     |     \  \          at the same height the parcel is COLDER than its surroundings
     |______\___\___ T   → denser → sinks back → STABLE, although T falls with height
   ```
   "**Stability compares two cooling rates: the environment's and a rising parcel's own.**"
33. `nb.note` — **C49 [B → C54]**: "In a fluid at rest hydrostatics (1.8) and the equation of state (1.12) link p, ρ and
    T, so one profile fixes the other two: give T(z), integrate dp/dz = −pg/(RT(z)), then ρ = p/RT. For T(z) = 288.15 −
    6.5 z/km: p(5 km) = 54.0 kPa, ρ(5 km) = 0.736 kg/m³ (`ch01.atmosphere_from_temperature`) — the USSA-1976 values."
34. `nb.note` — **C53 [B → C54]**: "The **lapse rate** is the environment's vertical temperature gradient. **Two sign
    conventions are in daily use; both are shown everywhere in this section:**" followed by the convention table (a
    markdown table in the note):
    | Convention | Definition | Standard troposphere | Dry adiabatic Γ_a | Stable when |
    | Kundu (this book, fluidpy code) | Γ ≡ dT/dz | −6.5 K/km | −9.8 K/km | Γ > Γ_a, e.g. −6.5 > −9.8 |
    | Meteorology (climate literature) | Γ ≡ −dT/dz | +6.5 K/km | +9.8 K/km | Γ < Γ_a, e.g. 6.5 < 9.8 |
    "**Converting: negate the number and flip the inequality** — Γ_met = −Γ_Kundu. The physics never changes."
35. `nb.primer[code]("inequalities under a sign change", …)` — P48: "Multiplying both sides of an inequality by a
    negative number reverses it: a > b ⇔ −a < −b." *demo:* `a, b = -6.5, -9.8; (a > b) == (-a < -b)` → True.
36. `nb.code[explain]` — *code:* `Ga = ch01.adiabatic_lapse_rate()`; `Ga_met = ch01.lapse_rate_convention(Ga,
    "meteorology")`; `for dTdz in (-6.5e-3, -12.0e-3, 0.0, Ga): k = ch01.lapse_rate_stability(dTdz); m =
    ch01.lapse_rate_stability(dTdz, convention="meteorology"); print(f"{k.verdict:9s} | Kundu: {k.text} | meteorology:
    {m.text}")`; `assert k.code == m.code`. *expect:* Ga = −9.761×10⁻³ K/m, Ga_met = +9.761×10⁻³; lines "stable | stable ⇔
    dT/dz > Γa: −6.5 > −9.8 K/km | stable ⇔ Γ < Γa: 6.5 < 9.8 K/km", "unstable | −12.0 < −9.8 | 12.0 > 9.8", "stable |
    0.0 > −9.8 | 0.0 < 9.8", "neutral | = | =". *explain:* the helper computes in Kundu's convention and only formats the
    other; the assert proves the verdict is convention-free.
37. `nb.primer[code]("chain rule", …)` — P49: "If f depends on x and y, and both change, df = (∂f/∂x)_y dx + (∂f/∂y)_x dy.
    Divide by dz to get a rate along a path." *demo:* f = x²y at (1, 2) with dx = 10⁻⁶, dy = 2×10⁻⁶: `((1+1e-6)**2*(2+2e-6)
    - 2)` vs `2*1*2*1e-6 + 1*2e-6` → 6.0×10⁻⁶ both.
38. `nb.primer[code]("Gibbs free energy", …)` — P50: "g ≡ h − Ts [J/kg]. We use it only as a device: from the Gibbs
    relation dh = T ds + v dp it has the neat differential dg = v dp − s dT." *demo:* sympy: `sp.expand(sp.symbols('T')*
    sp.symbols('ds') + sp.symbols('v')*sp.symbols('dp') - (sp.symbols('T')*sp.symbols('ds') + sp.symbols('s')*
    sp.symbols('dT')))` → `dp*v - dT*s`.
39. `nb.primer[code]("exact differentials and Maxwell relations", …)` — P51: "If dg = A dp + B dT for a smooth state
    function g, then A = (∂g/∂p)_T, B = (∂g/∂T)_p, and because mixed second derivatives are equal, (∂A/∂T)_p =
    (∂B/∂p)_T. Applied to dg = v dp − s dT this gives the Maxwell relation (∂v/∂T)_p = −(∂s/∂p)_T." *demo:* sympy
    `f = x**2*sp.sin(y)`: `sp.diff(f, x, y) == sp.diff(f, y, x)` → True.
40. `nb.derivation[check]("D19", "The adiabatic lapse rate without perfect-gas relations", ref="1.30", …)` — Part F ·
    D19 (16 steps, sympy check).
41. `nb.note` — **N17 [B → C54]**: "Earth's dry air: Γ_a = −g/C_p = −9.81/1004.7 = **−9.76 K/km in Kundu's convention =
    +9.76 K/km in the meteorology convention** (usually quoted as 9.8 K/km). Water is barely affected: α = 1.5×10⁻⁴ K⁻¹,
    C_p = 4190 J kg⁻¹ K⁻¹ at 283 K give Γ_a = −0.10 K/km."
42. `nb.worked_example("is the standard troposphere stable?", …)` — "Environment dT/dz = −6.5 K/km; Γ_a = −9.76 K/km.
    **Kundu:** −6.5 > −9.76 ✓ stable. **Meteorology:** Γ = +6.5, Γ_a = +9.76, 6.5 < 9.76 ✓ stable. Same margin 3.26 K/km
    in both. Lift a parcel 1 km from 288.15 K: it cools to 278.39 K, the environment there is 281.65 K — the parcel is
    3.26 K colder, so it sinks back. N² = (g/T)(dT/dz − Γ_a) = (9.81/288.15)(3.26×10⁻³) = 1.11×10⁻⁴ s⁻² → period 9.9 min.
    Over hot ground dT/dz = −12 K/km: −12 < −9.76 (meteorology 12 > 9.76) → unstable, e-folding 117 s at 300 K."
43. `nb.code[explain]` — *code:* `N2 = ch01.brunt_vaisala_sq_from_lapse(288.15, -6.5e-3)`; `ch01.stability_timescale(N2)`;
    `ch01.parcel_temperature(288.15, 1000.0)`; `288.15 - 6.5`; `ch01.brunt_vaisala_sq_from_lapse(300.0, -12e-3)`;
    `ch01.adiabatic_lapse_rate(T=283.15, cp=4190.0, alpha=1.5e-4)`. *expect:* 1.110×10⁻⁴; ('period', 596.4); 278.39 K;
    281.65; −7.32×10⁻⁵ → ('efold', 116.9); −9.94×10⁻⁵ K/m.
44. `nb.check_agree` — **from scratch** (curation §7): environment = USSA-1976 `T, p, _ =
    ch01.standard_atmosphere(z)` on z = 0…200 m (Δz = 1 m); a dry parcel lifted from the ground follows (1.26) with the
    environment's pressure: `T_parcel = 288.15*(p/p[0])**((ch01.GAMMA_AIR - 1)/ch01.GAMMA_AIR)`; `dTdz =
    np.gradient(T_parcel, z, edge_order=2)`; `assert np.isclose(dTdz[0], ch01.adiabatic_lapse_rate(), rtol=1e-3)`; `assert
    np.isclose(-dTdz[0], ch01.lapse_rate_convention(ch01.adiabatic_lapse_rate(), "meteorology"), rtol=1e-3)`. Markdown
    after: "At the release height parcel and environment coincide, so the slope is exactly (1.30). Higher up the parcel
    is colder than the air around it and its slope is −(g/C_p)(T_parcel/T_env) — (1.30) is the *local* rate where the two
    still match, the point D19 makes at step 13."
45. `nb.figure` — T(z) 0–3 km: environment profiles −6.5 K/km (black), −12 K/km (rose), +5 K/km inversion (teal) from
    288.15 K; dry adiabat through the ground point (dashed purple); a side table in the axes with both conventions for
    each profile and the verdict with its inequality (from `lapse_rate_stability`). *see:* three environment lines and one
    parcel line. *read:* where the environment line leans less than the dashed adiabat (cools more slowly), a lifted
    parcel ends up colder → stable; where it leans more → unstable. *change:* make the environment cool at exactly 9.76
    K/km and it lies on the adiabat: neutral, a parcel is at home at every height.
46. `nb.plotly[explain]` — slider over the environment dT/dz = −15…+10 K/km (26 steps): traces "environment T(z)" and
    "dry adiabat" (fixed) for z = 0…3 km. After `slider_figure(...)` the cell appends to every slider step a title update
    (`fig.layout.sliders[0].steps[i].args.append({"title.text": …})` — the helper's steps use method "update"), so the
    title shows e.g. "dT/dz = −6.5 > Γa = −9.8 K/km  ⇔  Γ_met = 6.5 < 9.8 K/km: **stable**, N = 0.0105 s⁻¹". *explain:*
    one call to `lapse_rate_stability` per slider value in each convention; the two texts are joined with ⇔. Reading
    notes: the two inequalities always agree; the flip happens at −9.8 K/km (+9.8 in meteorology).
47. `nb.md` — **What would change if…** "…you read a paper that says 'the lapse rate increases under warming'? In its
    convention Γ = −dT/dz grows, i.e. the air cools *faster* with height and moves *toward* Γ_a — less stable. In our
    code that is dT/dz becoming more negative. **Next:** a label that removes the parcel's own cooling altogether (C55)."

#### Block C55
48. `nb.core("C55", "Potential temperature θ: a temperature label that survives lifting", question="Can we label a parcel of
    air with a temperature that does not change when it rises or sinks without heating?")`
49. `nb.md` — **The problem in plain words:** "Is air at 500 hPa and −23 °C 'colder' than air at the surface at 15 °C?
    Bring it down and it warms by compression. Weather maps, isentropic analysis and all of dry atmospheric dynamics
    (Ch. 13) use a temperature with that effect removed."
50. `nb.md` — **The idea:** "Bring every parcel adiabatically to a reference pressure p_o = 1000 hPa and read its
    thermometer there." ASCII:
    ```
    500 hPa   T = 250 K  ──(compress adiabatically down to 1000 hPa)──►  θ = 305 K
    1000 hPa  T = 288 K  ──(already there)────────────────────────────►  θ = 288 K
    θ increases upward  ⇒  a lifted parcel is always cooler (in θ) than air above it  ⇒  stable
    ```
51. `nb.derivation("D20", "Potential temperature", ref="1.31", …)` — Part F · D20 (5 steps). Reminder before it: exponent
    rules (P43, C45) and (1.26) (C46, C45).
52. `nb.worked_example("two parcels", …)` — "θ = T(p_o/p)^{R/C_p}, R/C_p = 287.06/1004.7 = 0.2857 = 2/7. Parcel at 500 hPa,
    250 K: θ = 250 × 2^{0.2857} = 250 × 1.2190 = 304.8 K. Surface air 1000 hPa, 288 K: θ = 288 K. The upper parcel is 38 K
    colder in T but 17 K warmer in θ — the column is stable."
53. `nb.code[explain]` — *code:* `ch01.potential_temperature(250.0, 5.0e4)`; `z = np.linspace(0, 11e3, 111)`; `T, p, rho =
    ch01.standard_atmosphere(z)`; `theta = ch01.potential_temperature(T, p)`; `dth = np.gradient(theta, z,
    edge_order=2)`; `N2_theta = ch01.brunt_vaisala_sq_from_theta(theta, dth)`; `N2_lapse =
    ch01.brunt_vaisala_sq_from_lapse(T, -6.5e-3)`; `ch01.potential_temperature_gradient(T[0], -6.5e-3, p=p[0])`.
    *expect:* 304.75 K; θ rises from 287.06 K at the ground (p_o = 1000 hPa, not 1013) to ≈ 330 K at 11 km;
    N2_theta ≈ N2_lapse (1.11×10⁻⁴ at the ground, rtol 10⁻³); dθ/dz = 3.25 K/km.
54. `nb.check_agree` — **from scratch** (curation §7): `theta_mine = T*(ch01.P_REF/p)**((ch01.GAMMA_AIR - 1)/
    ch01.GAMMA_AIR)`; `assert np.allclose(theta_mine, ch01.potential_temperature(T, p))`; `assert
    np.allclose(ch01.temperature_from_potential(theta_mine, p), T)`.
55. `nb.note` — **C56 [B → C55]**: "Taking logs of (1.31), differentiating in z, and using dp/dz = −ρg, p = ρRT and
    α = 1/T gives (stated) the link between θ's gradient and the two lapse rates. In both conventions:" Equation
    $\dfrac{T}{\theta}\dfrac{d\theta}{dz} = \dfrac{dT}{dz} + \dfrac{g}{C_p} = \Gamma - \Gamma_a\ \text{(Kundu)} =
    \Gamma_{a,\rm met} - \Gamma_{\rm met}\ \text{(meteorology)}$ (1.32). "Standard troposphere at the ground: (θ/T)(−6.5 +
    9.76) K/km = 3.25 K/km."
56. `nb.note` — **N19 [C → C55]**: "The logarithmic derivative of (1.31) is step 5 of the next derivation; the same move
    gives the atmosphere's θ profiles in Ch. 13 §13.2."
57. `nb.primer[code]("logarithmic differentiation", …)` — P52: "d(ln f)/dz = (1/f) df/dz, and logs turn products and
    powers into sums: ln(ab^k) = ln a + k ln b. So a relative rate of a product is the sum of relative rates." *demo:*
    f = z³ at z = 2: `(np.log((2+1e-6)**3) - np.log(2**3))/1e-6` vs `3/2` → 1.5.
58. `nb.derivation("D36", "N² from the potential-temperature gradient", ref="links 1.29 and 1.32", …)` — Part F · D36
    (8 steps).
59. `nb.note` — **C57 [B → C55]**: "So an atmosphere is **stable, neutral or unstable as θ increases, stays constant or
    decreases with height** — the gradient of θ, not of T, decides. In a laboratory tank the difference hardly matters:
    over 1 m the adiabatic cooling is only 9.8×10⁻³ K, so T ≈ θ." Equation $N^2 = \dfrac{g}{\theta}\dfrac{d\theta}{dz}$.
60. `nb.note` — **N18 [B → C55]**: "A typical lower atmosphere (our synthetic profile, not digitised from the book): a
    mixed layer near the ground cooling at almost exactly Γ_a (near-neutral, stirred by turbulence), an inversion where T
    *rises* with height (very stable: smog and fog stay below it), and a stable layer above cooling more slowly than Γ_a.
    In θ: vertical near the ground, a sharp increase through the inversion, a steady increase above."
61. `nb.plotly[explain]` — slider over the inversion strength dT/dz_inv = −5…+20 K/km (26 steps) with
    `ch01.synthetic_boundary_layer_column(z, inversion_dT_dz=x/1000, mixed_dT_dz=ch01.adiabatic_lapse_rate())` on z =
    0…2 km (≤ 400 points; the mixed layer set exactly to Γ_a so it reads neutral): traces "T(z)" (black), "θ(z)" (purple),
    "θ where stable (N² > 0)" (teal markers), "θ where neutral/unstable (N² ≤ 0)" (rose markers); per-step title with
    N² in the inversion and its stability period, both lapse-rate conventions for the inversion. *explain:* the column
    function builds T, p, ρ, θ, N² layer by layer; markers are θ masked by `stability`. Reading notes: T falls then rises
    then falls; θ never decreases unless the inversion strength is below −9.8 K/km (then the "inversion" layer itself is
    unstable).
62. `nb.note` — **C58 [B → C55]**: "The same trick for density: the **potential density** ρ_θ is the density a parcel
    would have after an isentropic trip to p_o (stated): ρ = 0.70 kg/m³ at 500 hPa → ρ_θ = 1.148 kg/m³." Equation $\rho(z)
    = \rho_\theta(z)\,\big(p(z)/p_o\big)^{1/\gamma}$ (1.33). *(number via `ch01.potential_density(0.7, 5e4)` printed in
    item 65)*
63. `nb.note` — **N20 [C → C55]**: "Multiplying (1.31) and (1.33) and using p = ρRT, the exponents add to 1 and θρ_θ =
    p_o/R is the same for every parcel — one line of algebra, not used again."
64. `nb.note` — **C59 [C → C55]**: "Taking logs of θρ_θ = const gives (1.34), −(1/ρ_θ)dρ_θ/dz = (1/θ)dθ/dz: potential
    density must **decrease** upward for stability — the form the ocean uses (C61) and Ch. 13 §13.2 uses for vertical
    density variation."
65. `nb.code` — *code:* `ch01.potential_density(0.7, 5.0e4)`; `rho_th = ch01.potential_density(rho, p)`; `assert
    np.allclose(theta*rho_th, ch01.P_REF/ch01.R_AIR)` (N20 as a check). *expect:* 1.1485; the assert passes.
66. `nb.explainer("parcel_stability", heading="Push a parcel up: does it come back?", why="Stability links three things at
    once — the environment's profile, the parcel's own adiabat and ζ(t) — and the sign convention adds a fourth. Reshaping
    the profile while the parcel moves, and flipping Kundu ↔ meteorology while nothing physical changes, makes the
    criterion stick.", tries=[…])` — tries: "Start from 'standard atmosphere', press play, then switch the convention chip
    and read the status line in both forms."; "Drag dT/dz below −9.8 K/km and watch cos turn into cosh."; "Turn on the
    nocturnal inversion and release the parcel inside it."; "Ocean mode: toggle compressibility and compare the periods
    (Check question 4)."
67. `nb.md` — **What would change if…** "…the air were moist and clouds formed? Condensation releases heat, so a
    saturated parcel cools more slowly than Γ_a (≈ 5–6 K/km instead of 9.8) and a column stable for dry air can be
    unstable for cloudy air — beyond this chapter. **Next:** §1.11 asks what dimensional reasoning alone can tell us."

### A.11 §1.11 Dimensional Analysis — C64 (R01, R02 reminders), C67 (+C65, C66, C68), C69 (+N22, N23, C70, N24, C71, C72, N25, C73, C74, C75, C76 · D28 · E5), S01, S02
1. `nb.section("1.11", "Dimensional Analysis")` — intro "**What is this section about?** Because nature does not care
   about our units, every correct equation has matching dimensions, and every law can be written with dimensionless
   groups — exactly n − r of them. That turns a seven-variable experiment into a four-number one."

#### Block C64
2. `nb.core("C64", "Dimensional homogeneity: a correct law cannot care about units", question="Could a true law of nature
   give a different answer in feet than in metres?")`
3. `nb.md` — **The problem in plain words:** "A spacecraft was lost because one team worked in pound-seconds and another
   in newton-seconds. Checking that every term of an equation has the same dimensions catches such mistakes — and every
   test in this project does it with `pint`. The same principle lets us predict laws before solving anything."
4. `nb.md` — **The idea:** "**Units** are our choice (m, ft); **dimensions** are what is measured (length). Write a
   dimension as powers of the base dimensions, [q] = M^a L^b T^c Θ^d, i.e. an exponent vector (a, b, c, d). Adding two
   terms is only meaningful when their vectors are equal." Small table: | quantity | M | L | T | Θ | · p | 1 | −1 | −2 | 0 |
   · ρ | 1 | −3 | 0 | 0 | · g | 0 | 1 | −2 | 0 | · z | 0 | 1 | 0 | 0 |.
5. `nb.md` — reminders: "R01 (SI units, §1.2): a unit is a product of powers of base units. R02: prefixes only rescale
   numbers — they never change a dimension."
6. `nb.md` — **The maths, step by step:** "1. [ρgz] = (1, −3, 0, 0) + (0, 1, −2, 0) + (0, 1, 0, 0) = (1, −1, −2, 0) —
   exponents **add** when quantities multiply (P43). 2. [p₀] = (1, −1, −2, 0): the same, so p = p₀ − ρgz can be right.
   3. A wrong formula p = p₀ − ρgz² has (1, 0, −2, 0) ≠ [p]: rejected without any experiment. 4. Divide the correct law
   by p₀: p/p₀ = 1 − ρgz/p₀ — every term is now a pure number, the same in any unit system (conclusion 1: laws can be
   written dimensionlessly; conclusion 2: terms must match — dimensional homogeneity)."
7. `nb.worked_example("10 m of water in SI and in cgs", …)` — "SI: ρgz/p₀ = 1000 × 9.81 × 10 / 101 325 = 0.968. cgs: ρ =
   1 g/cm³, g = 981 cm/s², z = 1000 cm, p₀ = 1 013 250 dyn/cm² → 1 × 981 × 1000 / 1 013 250 = 0.968. The group is
   identical; the individual numbers changed by factors of 1000, 100 and 10."
8. `nb.code[explain]` — *code:* `ch01.dimension_vector("kg/m**3") + ch01.dimension_vector("m/s**2") +
   ch01.dimension_vector("m")`; `ch01.dimension_vector("Pa")`; `(Q_(1000, "kg/m**3")*Q_(9.81, "m/s**2")*Q_(10,
   "m")).to("Pa")`; a `try: (Q_(1000, "kg/m**3")*Q_(9.81, "m/s**2")*Q_(10, "m")**2).to("Pa") except
   pint.DimensionalityError as err: print("rejected:", err)`; `V = {"rho": "kg/m**3", "g": "m/s**2", "z": "m", "p0":
   "Pa"}`; `vals = {"rho": 1000.0, "g": 9.81, "z": 10.0, "p0": 101325.0}`; `group = {"rho": 1, "g": 1, "z": 1, "p0":
   -1}`; `ch01.group_value(group, vals)`, `ch01.group_value(group, ch01.rescale_units(vals, V, "cgs"))`,
   `ch01.group_value(group, ch01.rescale_units(vals, V, "imperial"))`. *expect:* [1 −1 −2 0] twice; 98100 pascal;
   "rejected: Cannot convert from 'kilogram / second ** 2' …"; 0.9682 three times.
9. `nb.figure` — exponent bars: two panels (left "p = p₀ − ρgz ✓", right "p = p₀ − ρgz² ✗"); grouped bars of the M, L,
   T, Θ exponents for each term (colours: M orange, L blue, T teal, Θ amber); mismatching bars outlined in rose. *see:*
   in the left panel every term's bars are identical; on the right one term sticks out in L. *read:* equal bars for every
   term = dimensionally homogeneous. *change:* add a term ½ρu² (dynamic pressure) and its bars match p's too — Bernoulli
   (Ch. 4) passes this test.
10. `nb.md` — **What would change if…** "…a problem had seven variables? Checking one equation is easy; finding *all*
    unit-free combinations needs bookkeeping. **Next:** put the exponent vectors side by side as the columns of a matrix
    (C67)."

#### Block C67
11. `nb.core("C67", "The dimensional matrix: dimensions as columns of numbers", question="How do we write down what
    dimensions each variable has, so that linear algebra can reason about them?")`
12. `nb.md` — **The problem in plain words:** "An engineer wants the pressure drop along a pipe for any fluid, diameter,
    roughness and speed. Seven variables at ten values each would be 10⁷ experiments. Dimensional analysis promises far
    fewer — and it starts by tabulating dimensions."
13. `nb.md` — **The idea:** "Stack the exponent vectors of C64 as columns: rows = base dimensions, columns = variables.
    Everything dimensional analysis needs is now in one integer matrix."
14. `nb.note` — **C65 [B → C67]**: "**Step 1 — choose the variables.** The most important step: one solution variable
    (Δp) plus everything it can depend on — geometry (Δx, d, ε), flow (U) and material (ρ, μ). Leave one out and no
    amount of algebra recovers it (dropping μ is tested in the C69 figure); add unnecessary ones and the result gets
    weaker." Equation $f(\Delta p, \Delta x, d, \varepsilon, U, \rho, \mu) = 0$ (1.38).
15. `nb.note` — **C66 [B → C67]**: "**Step 2 — dimensions.** [q] means 'the dimension of q', written with M, L, T and Θ:
    [U] = L/T, [p] = M/(LT²), [C_p] = L²/(ΘT²). Temperature often appears only as k_BΘ, RΘ or C_pΘ, whose dimension is
    L²/T², so three base dimensions suffice for most flows. (The kilomole in R_u is treated as a pure number, as the
    book's Example 1.2 does; `dimension_vector` says so with a warning.)"
16. `nb.primer[code]("matrices, determinants and minors", …)` — P53: "A matrix is a table of numbers; the determinant of
    a square one is a single number that is zero exactly when its columns are dependent. A **minor** is the determinant
    of a square piece cut from a bigger matrix. 3×3 determinants by **cofactor expansion** along the first row:
    a₁₁(a₂₂a₃₃ − a₂₃a₃₂) − a₁₂(a₂₁a₃₃ − a₂₃a₃₁) + a₁₃(a₂₁a₃₂ − a₂₂a₃₁)." *demo:* `np.linalg.det(np.array([[2, 0, 0],
    [0, 3, 0], [0, 0, 4]]))` → 24.0.
17. `nb.primer[code]("linear independence and rank", …)` — P54: "Vectors are independent if none is a combination of the
    others. The **rank** of a matrix is the number of independent columns (= independent rows) = the size of its largest
    nonzero minor." *demo:* `np.linalg.matrix_rank(np.array([[1, 2], [2, 4]]))` → 1 (second column = 2 × first).
18. `nb.md` — **The maths, step by step:** "1. Columns in the order of (1.38): Δp, Δx, d, ε, U, ρ, μ. 2. [Δp] = M L⁻¹ T⁻²
    → column (1, −1, −2); [Δx] = [d] = [ε] = L → (0, 1, 0); [U] → (0, 1, −1); [ρ] → (1, −3, 0); [μ] = M L⁻¹ T⁻¹ →
    (1, −1, −1). 3. The matrix:" Equation $\begin{array}{c|ccccccc} & \Delta p & \Delta x & d & \varepsilon & U & \rho &
    \mu\\\hline \mathrm M & 1&0&0&0&0&1&1\\ \mathrm L & -1&1&1&1&1&-3&-1\\ \mathrm T & -2&0&0&0&-1&0&-1\end{array}$ (1.39)
    "4. Its rank r is at most 3 (three rows)."
19. `nb.note` — **C68 [B → C67]**: "**Step 3 — rank by minors.** r is the size of the largest square sub-matrix with a
    nonzero determinant. The first three columns (Δp, Δx, d) give 0 — Δx and d are both pure lengths, so those columns
    are dependent — but the last three (U, ρ, μ) give −1 ≠ 0, so r = 3 (stated; worked below). A rank below the number of
    rows happens when one row is a combination of others, e.g. in statics problems where mass enters only through force."
20. `nb.worked_example("the (U, ρ, μ) minor by hand", …)` — "Rows M, L, T of columns U, ρ, μ: [[0, 1, 1], [1, −3, −1],
    [−1, 0, −1]]. Cofactor expansion along row 1: 0·((−3)(−1) − (−1)(0)) − 1·((1)(−1) − (−1)(−1)) + 1·((1)(0) − (−3)(−1))
    = 0 − 1·(−1 − 1) + (0 − 3) = 2 − 3 = **−1** ≠ 0 → r = 3. The first-three-columns minor [[1, 0, 0], [−1, 1, 1],
    [−2, 0, 0]]: expansion gives 1·(1·0 − 1·0) − 0 + 0 = 0."
21. `nb.primer[code]("itertools.combinations", …)` — P55: "`combinations(range(7), 3)` lists every way to choose 3 of 7
    column indices, without repeats — 35 of them." *demo:* `len(list(itertools.combinations(range(7), 3)))` → 35.
22. `nb.code[explain]` — *code:* `A, names, rows = ch01.dimensional_matrix(ch01.PIPE)`; `print(rows, names); print(A)`;
    `r, ri, ci = ch01.rank_by_minors(A)`; `ch01.minor_determinant(A, (0, 1, 2), (4, 5, 6))`; `ch01.minor_determinant(A,
    (0, 1, 2), (0, 1, 2))`. *expect:* ['M', 'L', 'T'] ['dp', 'dx', 'd', 'eps', 'U', 'rho', 'mu']; the matrix (1.39);
    r = 3 with witness rows (0, 1, 2), columns (0, 1, 4) = (Δp, Δx, U) (first nonzero in lexicographic order); −1; 0.
23. `nb.primer[code]("np.linalg.det and np.linalg.matrix_rank", …)` — P56: "numpy's floating-point versions: `det`
    returns e.g. −0.9999999999999998 for an integer −1 (round-off), `matrix_rank` counts singular values above a
    tolerance." *demo:* `round(np.linalg.det(A[:, 4:7]))` → −1.
24. `nb.check_agree` — **from scratch** (curation §7): hand-written `det3(M)` by cofactor expansion; `nonzero =
    [c for c in itertools.combinations(range(7), 3) if det3(A[:, c]) != 0]`; `assert len(nonzero) > 0` → r = 3; `assert
    3 == ch01.rank_by_minors(A)[0] == np.linalg.matrix_rank(A)`; `assert all(det3(A[:, c]) == round(np.linalg.det(A[:,
    c])) for c in itertools.combinations(range(7), 3))`. *expect:* 26 of the 35 minors are nonzero (the zero ones contain
    two of the three pure lengths Δx, d, ε, or all three).
25. `nb.figure` — annotated heatmap of (1.39) (`ax.imshow`, diverging colours, integers printed in cells), columns
    labelled with the symbols, the (U, ρ, μ) minor outlined in accent purple with "det = −1", the (Δp, Δx, d) minor
    outlined in grey with "det = 0"; a small bar at the right: "nonzero 3×3 minors: 26 / 35". *see:* a 3 × 7 integer
    grid with two outlined blocks. *read:* one nonzero 3×3 block is enough to prove r = 3; the three identical L-only
    columns are what make many minors vanish. *change:* add temperature T as an 8th variable and a Θ row appears — r
    would become 4 (Example 1.2 has r = 4).
26. `nb.md` — **What would change if…** "…we asked which exponent combinations make a column combination dimensionless?
    That is solving A·k = 0 — the null space of this matrix, and the Π theorem (C69)."

#### Block C69
27. `nb.core("C69", "Buckingham's Π theorem: a law depends on n − r numbers", question="Why can seven pipe-flow variables be
    squeezed into four dimensionless numbers — and why exactly four?")`
28. `nb.md` — **The problem in plain words:** "Wind-tunnel models of aircraft, scale models of harbours, rotating-tank
    experiments for the atmosphere: a small experiment stands for a big system when the right dimensionless numbers
    match. Reynolds, Froude and Rossby numbers (Ch. 4 §4.11, Ch. 13) all come from this theorem."
29. `nb.md` — **The idea:** "A product of powers q₁^{k₁}⋯q_n^{k_n} is dimensionless exactly when its exponent vector k
    solves A·k = 0. The solutions form a subspace — the **null space** — whose dimension is n − r. Each basis vector is
    one Π group." ASCII: `7 variables  →  matrix of rank 3  →  null space of dimension 7 − 3 = 4  →  Π₁, Π₂, Π₃, Π₄`.
30. `nb.note` — **N22 [C → C69]**: "The starting point is a relation among the n variables, written as the input
    dictionary of `pi_groups`; Ch. 4 §4.11 starts from the same form." Equation $f(q_1, q_2, \dots, q_n) = 0$ (1.36).
31. `nb.note` — **N23 [B → C69]**: "**Step 4 — count the groups:** n − r. Pipe: 7 − 3 = 4."
32. `nb.primer[code]("null space and rank–nullity theorem", …)` — P58: "The null space of A is the set of all vectors k
    with A·k = 0. **Rank–nullity:** for a matrix with n columns, rank + (dimension of the null space) = n." *demo:*
    `sp.Matrix([[1, 1, 0], [0, 0, 1]]).nullspace()` → one vector (1 = 3 − 2).
33. `nb.primer[code]("scaling the base units", …)` — P59: "Change the unit of length by a factor λ_L (m → cm: λ_L = 100),
    of mass by λ_M, of time by λ_T. A quantity with dimension M^a L^b T^c gets its **number** multiplied by λ_M^a λ_L^b
    λ_T^c; the physical quantity is unchanged." *demo:* density 1000 kg/m³ in g/cm³: `1000 * 1000**1 * 100**-3` → 1.0.
34. `nb.primer[code]("sympy Matrix and nullspace", …)` — P61: "`sp.Matrix(...)` keeps integers exact; `.rank()` and
    `.nullspace()` return exact results (fractions, no round-off)." *demo:* `sp.Matrix([[1, 0, 0, 0, 0, 1, 1], [-1, 1, 1,
    1, 1, -3, -1], [-2, 0, 0, 0, -1, 0, -1]]).rank()` → 3.
35. `nb.derivation[check]("D28", "Buckingham's Π theorem", ref="1.37", …)` — Part F · D28 (12 steps, sympy check).
36. `nb.note` — **C70 [B → C69]**: "**Step 5 — build the groups (exponent algebra).** Choose r = 3 repeating variables
    whose minor is nonzero (U, d, ρ). Each group is one other variable times powers of them. For Δp: Π₁ = Δp U^a d^b ρ^c;
    matching M: c + 1 = 0, L: a + b − 3c − 1 = 0, T: −a − 2 = 0 gives a = −2, b = 0, c = −1 (stated; solved in code
    below). The same with Δx, ε and μ:" Equation $\Pi_1 = \dfrac{\Delta p}{\rho U^2},\ \Pi_2 = \dfrac{\Delta x}{d},\ \Pi_3
    = \dfrac{\varepsilon}{d},\ \Pi_4 = \dfrac{\mu}{\rho U d}$.
37. `nb.note` — **N24 [C → C69]**: "The same groups can be found **by inspection**, cancelling M, L and T one ratio at a
    time — how Navier–Stokes is made dimensionless in Ch. 4 §4.11."
38. `nb.note` — **C71 [B → C69]**: "Products and powers of groups are groups too: Δp d²ρ/μ² = Π₁/Π₄² and ε/Δx = Π₃/Π₂.
    Any such set may replace the original, but **only n − r of them are independent** — a fifth is always a combination
    of four (`groups_independent` → False)."
39. `nb.note` — **C72 [B → C69]**: "**Step 6 — the dimensionless law:**" Equation $\dfrac{\Delta p}{\rho U^2} =
    \varphi\!\left(\dfrac{\Delta x}{d}, \dfrac{\varepsilon}{d}, \dfrac{\mu}{\rho U d}\right)$ (1.40). "Π₄ is the inverse
    Reynolds number, 1/Re."
40. `nb.note` — **N25 [B → C69]**: "**Step 7 — add physics.** Far from the inlet, pressure drop is proportional to
    length, so Δp/ρU² = (Δx/d)·φ₂(ε/d, Re): three numbers become two. For slow laminar flow in a smooth pipe the full
    solution (Ch. 8) gives φ₂ = 32/Re — the line the collapse figure below follows."
41. `nb.worked_example("Π₁ in two unit systems", …)` — "Δp = 100 Pa, ρ = 1000 kg/m³, U = 0.1 m/s: Π₁ = 100/(1000 × 0.01) =
    10. In cgs: Δp = 1000 dyn/cm², ρ = 1 g/cm³, U = 10 cm/s: Π₁ = 1000/(1 × 100) = 10. With d = 1 cm and μ = 10⁻³ Pa s:
    Π₄ = 10⁻³/(1000 × 0.1 × 0.01) = 10⁻³, i.e. Re = 1000."
42. `nb.primer[code]("fractions.Fraction", …)` — P60: "`Fraction(-1, 2)` is an exact rational number; fluidpy returns
    exponents this way so ½ never becomes 0.49999." *demo:* `Fraction(1, 3) + Fraction(1, 6)` → `Fraction(1, 2)`.
43. `nb.code[explain]` — *code:* `groups = ch01.pi_groups(ch01.PIPE, solution="dp", repeating=("U", "d", "rho"))`;
    `[ch01.group_latex(g) for g in groups]`; `ch01.solve_exponents("dp", ("U", "d", "rho"), ch01.PIPE)`; `vals = {"dp":
    100.0, "dx": 1.0, "d": 0.01, "eps": 1e-5, "U": 0.1, "rho": 1000.0, "mu": 1e-3}`; `[ch01.group_value(g, vals) for g
    in groups]`; same with `ch01.rescale_units(vals, ch01.PIPE, "cgs")` and `"imperial"`; `extra = {"dp": 1, "d": 2,
    "rho": 1, "mu": -2}`; `ch01.groups_independent(groups + [extra])`; `ch01.groups_independent(groups)`. *expect:*
    ['\\frac{\\Delta p}{U^{2} \\rho}', '\\frac{\\Delta x}{d}', '\\frac{\\varepsilon}{d}', '\\frac{\\mu}{U d \\rho}'];
    {dp: 1, U: −2, d: 0, rho: −1}; [10.0, 100.0, 0.001, 0.001] in all three systems; False; True.
44. `nb.primer[code]("np.linalg.solve", …)` — P57: "`np.linalg.solve(B, b)` finds x with B·x = b for a square, nonsingular
    B." *demo:* `np.linalg.solve(np.array([[2., 0.], [0., 4.]]), np.array([2., 8.]))` → `[1. 2.]`.
45. `nb.check_agree` — **from scratch** (curation §7): `B = A[:, [names.index(k) for k in ("U", "d", "rho")]]`; `b =
    -A[:, names.index("dp")]`; `x = np.linalg.solve(B.astype(float), b.astype(float))`; `assert np.allclose(x, [-2, 0,
    -1])`; `lib = ch01.solve_exponents("dp", ("U", "d", "rho"), ch01.PIPE)`; `assert np.allclose(x, [float(lib[k]) for k
    in ("U", "d", "rho")])`.
46. `nb.figure` — **units invariance**: 60 random pipe states (seeded, log-uniform ranges); x = each Π in SI, y = the
    same Π recomputed after `rescale_units` to cgs (teal) and imperial (orange), all on the 1:1 line, log–log. *see:*
    every point on the diagonal. *read:* a group's value does not depend on units, while the raw numbers (inset histogram
    of Δp in Pa vs dyn/cm²) do. *change:* a group that is *not* dimensionless, e.g. Δp/ρU, would scatter off the line by
    the length factor.
47. `nb.figure` — **data collapse** (1×2): synthetic laminar data from `ch01.poiseuille_pressure_drop(mu, U, dx, d)`
    ("result derived in Ch. 8; data generator only") for water, glycerine and air × 4 diameters × U spanning Re = 1…2000:
    left raw Δp/Δx vs U (scattered families); right Π₁·(d/Δx) vs Re (all points on 32/Re, grey line); the markdown adds that without μ the only groups left are Δx/d and ε/d, identical for every point of one pipe,
    so nothing can collapse the three decades of scatter in Π₁. *see:* many curves → one line. *read:* the
    dimensionless law needs only Re here; that is why engineers plot friction factors against Re. *change:* turbulent
    rough-pipe data would spread by ε/d — Π₃ would matter.
48. `nb.note` — **C73 [B → C69]** (Ex. 1.2, stated): "Scale height from dimensions alone: variables H, T₀, M_w, g, R_u
    (n = 5); M, L, T, Θ all appear and r = 4, so there is one group, HgM_w/(R_uT₀), and a lone group must be a constant:
    H = const · R_uT₀/(gM_w). The isothermal atmosphere (C63) fixes const = 1: 7.32 km at 250 K."
49. `nb.note` — **C74 [C → C69]**: "Even Pythagoras follows: a right triangle's area is C²φ(β) for hypotenuse C and angle
    β, and the two smaller similar triangles (*gloss:* same angles, sides in one ratio) add up to it, giving A² + B² = C²
    (φ(β) = ¼ sin 2β as a check) — a curiosity not used later."
50. `nb.note` — **C75 [B → C69]** (Ex. 1.4, stated): "G. I. Taylor estimated a nuclear blast's energy from photos: E, ρ, D
    and t give r = 3 and one group, so E = KρD⁵/t² and D ∝ t^{2/5}; K comes only from the full similarity solution
    (*gloss:* a solution whose shape stays the same as it grows), K ≈ 0.86 for air (Taylor 1950). A front 100 m wide at
    25 ms in air of 1.2 kg/m³: E = 1.9×10¹³ K joules ≈ 4.6K kilotons of TNT."
51. `nb.note` — **C76 [B → C69]** (Ex. 1.5, stated): "Why the sky is blue: scattered intensity S from a particle of volume
    V at distance d, light of wavelength λ (*gloss:* colour; intensity ∝ wave amplitude²). Only L and the intensity
    dimension appear (the T row is −3 times the M row, so r = 2); energy spreading forces d⁻², a dipole forces V², and
    then the groups force λ⁻⁴. Blue (450 nm) scatters (700/450)⁴ = 5.9 times more than red (700 nm)." Equation
    $\dfrac{S}{I} = \dfrac{V^2}{d^2\lambda^4}\,\varphi_3(n_s)$.
52. `nb.code[explain]` — numbers for C73–C76: `ch01.rank_by_minors(ch01.dimensional_matrix(ch01.SCALE_HEIGHT)[0])[0]`,
    `[ch01.group_latex(g) for g in ch01.pi_groups(ch01.SCALE_HEIGHT)]`, `ch01.R_U*250/(ch01.G0*ch01.M_W_AIR)` vs
    `ch01.scale_height(250.0)`; `ch01.pythagoras_phi(np.radians(30))` vs `0.25*np.sin(np.radians(60))`;
    `ch01.blast_energy(100.0, 0.025, 1.2)/4.184e12`, `ch01.TAYLOR_K_GAMMA14`; `ch01.rayleigh_scattering_ratio(1, 1,
    450e-9)/ch01.rayleigh_scattering_ratio(1, 1, 700e-9)`; `ch01.rank_by_minors(ch01.dimensional_matrix(ch01.RAYLEIGH)[0])[0]`.
    *expect:* 4; ['\\frac{H M_w g}{R_u T_0}']; 7317.9 m twice; 0.2165 twice; 4.59 kt; 0.856; 5.855; 2.
53. `nb.explainer("buckingham_pi_machine", heading="Why can 7 pipe variables shrink to 4 numbers?", why="Toggling
    variables, swapping the repeating set and switching units changes the matrix, its rank and the groups at once — you
    experiment with the method instead of watching it done once.", tries=[…])` — tries: "Choose the repeating set (d, ε,
    ρ) and read why the machine refuses."; "Remove μ: how many groups are left, and which physics is lost?"; "Switch to
    the blast preset: what power of t does the radius follow?"; "Change SI → imperial and watch every Π value stay put."
54. `nb.md` — **What would change if…** "…the problem included rotation (the Earth's Ω) as an eighth variable? One more
    column, the same rank, one more group — the Rossby number of Ch. 13."
55. `nb.pointer` — **S01 [C, SKIP]**: "The chapter's 30 exercises are in the book and are not reproduced here. Of the
    derivations the book leaves to Exercises 1.10, 1.11, 1.13 and 1.14, D14 (isentropic law), D18 (parcel equation) and
    D19 (adiabatic lapse rate) are written out above; the result of Exercise 1.10 (e = e(T) for a perfect gas) is stated
    in the C40 block."
56. `nb.pointer` — **S02 [C, SKIP]**: "The book's reading list is at the end of the printed chapter; the public data we
    used (CODATA, USSA-1976, IAPWS, Jennings, Taylor) are cited in reference/ch01/SOURCES.md."
57. `nb.summary` — *clicked* (one per A item): "A density 'at a point' is the plateau of box averages between molecular
    noise and the flow's own variation (C06)." · "Viscosity is momentum diffusion: τ = μ du/dy, and ν = μ/ρ sets how fast
    motion spreads (C12)." · "In still fluid pressure falls with height at exactly ρg per metre; buoyancy is the resulting
    face-pressure difference (C20)." · "Heat plus work changes internal energy; q and w depend on the route, Δe does not
    (C25)." · "Entropy is the route-independent sum of dq/T, and the Gibbs relations tie state functions for any process
    (C35)." · "Sound speed is √((∂p/∂ρ)_s); an incompressible fluid has c → ∞ (C36)." · "p = ρRT is the molecular gas law
    per unit mass, R = R_u/M_w (C40)." · "With no heat and no friction a perfect gas follows p ∝ ρ^γ (C45)." · "A displaced
    parcel obeys ζ″ + N²ζ = 0 (C50)." · "N² > 0 oscillates, = 0 stays, < 0 runs away (C51)." · "Air is stable when it
    cools more slowly than 9.8 K/km: dT/dz > −9.8 K/km (Kundu) ⇔ Γ_met < +9.8 K/km (meteorology) (C54)." · "Potential
    temperature removes adiabatic cooling; stability is θ increasing upward, N² = (g/θ)dθ/dz (C55)." · "Every term of a
    correct law has the same dimensions (C64)." · "Dimensions of all variables form an integer matrix whose rank counts
    independent dimensions (C67)." · "Dimensionless groups are the null space of that matrix: n − r of them (C69)."
    *feeds_forward:* "Ch. 2–3: fields and fluid particles on the continuum" · "Ch. 4: tensor viscosity law, energy
    equation, Boussinesq (α), dimensionless Navier–Stokes (Π)" · "Ch. 7 and 13: N² for internal waves; θ and hydrostatics
    for GFD" · "Ch. 8: Couette start-up and Poiseuille results used here as checks" · "Ch. 15: c, isentropic relations".
    *left_out:* "the 30 exercises (book)" · "moist thermodynamics (beyond Kundu Ch. 1)" · "the proof of (1.19) (Ch. 15 §15.2)".

### A.12 Placement table (ID → notebook section → host block → call) — all 106 curation rows
| ID | Depth | Section | Host block | Call (Part A row) |
|---|---|---|---|---|
| C01 | B | §1.1 | tagged → C06 | A.1 #2 note |
| R01 | B | §1.2 | recap (→ C64; reminded A.11 #5) | A.2 #2 recap |
| R02 | C | §1.2 | recap (→ C64; reminded A.11 #5) | A.2 #6 recap |
| R03 | B | §1.2 | recap (→ C40; reminded A.9 #5) | A.2 #7 recap |
| C02 | B | §1.3 | tagged → C12 (reminded A.5 #5) | A.3 #3 note |
| N01 | C | §1.3 | tagged → C12 | A.3 #6 note |
| C03 | B | §1.3 | tagged → C20 (reminded A.7 #19) | A.3 #7 note |
| C04 | B | §1.3 | tagged → C06 | A.3 #9 note |
| N02 | C | §1.3 | tagged → C12 | A.3 #5 figure |
| C05 | B | §1.4 | C06 | A.4 #9 note |
| C06 | A | §1.4 | CORE | A.4 #2 core |
| C07 | B | §1.4 | C06 | A.4 #26 note |
| N03 | B | §1.4 | C06 | A.4 #27 note |
| C08 | B | §1.5 | C12 | A.5 #7 note |
| C09 | C | §1.5 | C12 | A.5 #8 note |
| C10 | B | §1.5 | C12 | A.5 #9 note |
| N04 | C | §1.5 | C12 | A.5 #10 note |
| C11 | B | §1.5 | C12 | A.5 #11 note |
| C12 | A | §1.5 | CORE | A.5 #2 core |
| N05 | B | §1.5 | C12 | A.5 #24 note |
| C13 | B | §1.5 | C12 | A.5 #14 note |
| C14 | B | §1.5 | C12 | A.5 #15 note |
| N06 | C | §1.5 | C12 | A.5 #25 note |
| C15 | B | §1.6 | tagged → C20 | A.6 #2 note |
| N07 | B | §1.6 | tagged → C20 | A.6 #3 note |
| C16 | B | §1.6 | tagged → C20 | A.6 #4 note |
| N08 | C | §1.6 | tagged → C20 | A.6 #6 note |
| N09 | C | §1.6 | tagged → C20 | A.6 #7 note |
| C17 | B | §1.7 | C20 | A.7 #10 note |
| C18 | B | §1.7 | C20 | A.7 #11 note |
| N10 | C | §1.7 | C20 | A.7 #12 note |
| C19 | B | §1.7 | C20 | A.7 #13 note |
| C20 | A | §1.7 | CORE | A.7 #2 core |
| C21 | B | §1.7 | C20 | A.7 #16 note |
| N11 | B | §1.7 | C20 | A.7 #6 note (+ #5 figure) |
| C22 | B | §1.7 | C20 | A.7 #20 note |
| N12 | C | §1.7 | C20 | A.7 #21 note |
| C23 | B | §1.8 | C25 | A.8 #5 note |
| C24 | B | §1.8 → §1.4 | C06 | A.4 #28 note |
| C25 | A | §1.8 | CORE | A.8 #2 core |
| C26 | B | §1.8 | C25 | A.8 #11 note |
| N13 | C | §1.8 | C25 | A.8 #10 note |
| C27 | B | §1.8 | C25 | A.8 #12 note |
| C28 | B | §1.8 | C25 | A.8 #13 note (+ A.9 #14 3-D surface) |
| C29 | B | §1.8 | C35 | A.8 #27 note |
| C30 | B | §1.8 | C35 | A.8 #28 note |
| C31 | B | §1.8 | C35 | A.8 #28 note |
| C32 | C | §1.8 | C35 | A.8 #29 note |
| C33 | B | §1.8 | C35 | A.8 #30 note |
| C34 | B | §1.8 | C35 | A.8 #34 note |
| N14 | C | §1.8 | C35 | A.8 #35 note |
| N15 | B | §1.8 | C35 | A.8 #31 note |
| C35 | A | §1.8 | CORE | A.8 #23 core |
| C36 | A | §1.8 | CORE | A.8 #42 core |
| C37 | B | §1.8 | C36 | A.8 #50 note |
| C38 | B | §1.9 | C40 | A.9 #6 note |
| C39 | B | §1.9 | C40 | A.9 #7 note |
| C40 | A | §1.9 | CORE | A.9 #2 core |
| C42 | B | §1.9 | C45 | A.9 #24 note |
| C43 | B | §1.9 | C45 | A.9 #25 note |
| N16 | B | §1.9 | C45 | A.9 #26 note |
| C41 | B | §1.9 | C40 | A.9 #15 note |
| C44 | B | §1.9 | C45 | A.9 #27 note |
| C45 | A | §1.9 | CORE | A.9 #21 core |
| C46 | B | §1.9 | C45 | A.9 #35 note |
| C47 | B | §1.9 | C45 | A.9 #36 note |
| C48 | B | §1.9 | C40 | A.9 #16 note |
| C49 | B | §1.10 | C54 | A.10 #33 note |
| C50 | A | §1.10 | CORE | A.10 #2 core |
| C51 | A | §1.10 | CORE | A.10 #13 core |
| C52 | B | §1.10 | C51 | A.10 #18 note |
| C53 | B | §1.10 | C54 | A.10 #34 note (convention table) |
| C54 | A | §1.10 | CORE | A.10 #30 core |
| N17 | B | §1.10 | C54 | A.10 #41 note |
| N18 | B | §1.10 | C55 | A.10 #60 note (+ #61 slider) |
| C55 | A | §1.10 | CORE | A.10 #48 core |
| N19 | C | §1.10 | C55 | A.10 #56 note |
| C56 | B | §1.10 | C55 | A.10 #55 note |
| C57 | B | §1.10 | C55 | A.10 #59 note |
| C58 | B | §1.10 | C55 | A.10 #62 note |
| N20 | C | §1.10 | C55 | A.10 #63 note |
| C59 | C | §1.10 | C55 | A.10 #64 note |
| C60 | B | §1.10 | C51 | A.10 #25 note |
| N21 | B | §1.10 | C51 | A.10 #26 note |
| C61 | B | §1.10 | C51 | A.10 #27 note |
| C62 | B | §1.10 → §1.9 | C40 | A.9 #17 note (+ #19 slider) |
| C63 | B | §1.10 → §1.9 | C40 | A.9 #18 note |
| C64 | A | §1.11 | CORE | A.11 #2 core |
| N22 | C | §1.11 | C69 | A.11 #30 note |
| C69 | A | §1.11 | CORE | A.11 #27 core |
| C65 | B | §1.11 | C67 | A.11 #14 note |
| C66 | B | §1.11 | C67 | A.11 #15 note |
| C67 | A | §1.11 | CORE | A.11 #11 core |
| C68 | B | §1.11 | C67 | A.11 #19 note |
| N23 | B | §1.11 | C69 | A.11 #31 note |
| C70 | B | §1.11 | C69 | A.11 #36 note |
| N24 | C | §1.11 | C69 | A.11 #37 note |
| C71 | B | §1.11 | C69 | A.11 #38 note |
| C72 | B | §1.11 | C69 | A.11 #39 note |
| N25 | B | §1.11 | C69 | A.11 #40 note |
| C73 | B | §1.11 | C69 | A.11 #48 note |
| C74 | C | §1.11 | C69 | A.11 #49 note |
| C75 | B | §1.11 | C69 | A.11 #50 note |
| C76 | B | §1.11 | C69 | A.11 #51 note |
| S01 | C | Ex | end of §1.11 | A.11 #55 pointer |
| S02 | C | Ref | end of §1.11 | A.11 #56 pointer |

Derivation placement (all 12 inside the curation's CORE block): D34, D35 → C06 · D05, D37 → C20 · D10 → C35 · D11 → C40
· D14 → C45 · D18 → C50 · D19 ★★★ (check) → C54 · D20, D36 → C55 · D28 ★★★ (check) → C69. The 25 demoted derivations
(§4c) appear only as "(stated)" sentences in the notes named in the curation.

Visual per A block (≥ 1 each): C06 figure + animation + plotly + E1 · C12 animation + E2 · C20 figure ×2 + plotly · C25
figure + animation · C35 figure · C36 figure · C40 figure + plotly 3-D + plotly slider · C45 figure + E3 · C50 animation ·
C51 figure + live · C54 figure + plotly slider · C55 plotly slider + E4 · C64 figure · C67 figure · C69 figure ×2 + E5.

## Part B — explainer storyboards

Common to all six: built with `tools/new_viz.py`, `viz:chapter ch01`, colours of convention 5, KaTeX with double
backslashes in JS strings, constants identical to fluidpy (`G0 = 9.80665`, `R_AIR = 287.0579959589441`, `CP_AIR =
1004.7029858563045`, `CV_AIR = 717.6449898973605`, `K_B = 1.380649e-23`, `N_A_KMOL = 6.02214076e26`, `P_REF = 1e5`,
`P_ATM = 101325`) so parity rows can use `rtol ≤ 1e-9`. Tabs: Walkthrough · Explore · Explain · Derivation (when listed)
· Equations · Code · Check. ≤ 5 Explore controls (the rest `optional`), readout labels ≤ 22 chars, step texts ≤ 45 words,
derivation lines one relation each. Every parity `py:` expression below is evaluated by `tools/shot.py` with no builtins
(dict/tuple literals and indexing are fine).

### E1 · continuum_averaging_volume
- **Title:** "When does 'density at a point' make sense?" · **Summary:** "Grow a sampling box from 1 nm to 1 m through a gas
  or a liquid and watch its density reading go from wild molecular noise to a steady plateau to the flow's own
  variation." · **CORE:** C06 (also C04, C05, C07, N03, C24) · **Reference:** `angular_frequency_explorer_1.html` (linked
  views on one clock, modes, a highlighted real-world table, "Right now" notes).
- **meta:** `viz:sections 1.4` · `viz:equations §1.4 ρ=δm/δV, Kn=l/L` · `viz:fluidpy ch01.sample_density
  ch01.density_noise_expected ch01.box_average_density ch01.mean_free_path_air ch01.knudsen_number` · `viz:derivations
  D35`.
- **Physics (JS ↔ Python):** `numberDensity(rho, Mw)` ↔ `ch01.number_density(rho, M_w)` · `molecularMass(Mw)` ↔
  `ch01.molecular_mass(M_w)` · `noiseExpected(L, n)` ↔ `ch01.density_noise_expected(L, number_density)` ·
  `boxAverage(L, rho0, eps, Lflow, x0)` ↔ `ch01.box_average_density(L, rho0, variation, L_flow, x0)` ·
  `sutherland(T)` ↔ `ch01.sutherland_viscosity(T)` · `meanFreePathAir(T, p)` ↔ `ch01.mean_free_path_air(T, p)` (Jennings
  with Sutherland μ and ρ = p/RT) · `knudsen(l, L)` ↔ `ch01.knudsen_number(l, L)` · `sampleReading(L, n, m, eps, Lflow,
  rng)` mirrors the method of `ch01.sample_density` (Poisson draw via `Viz.rng` for λ ≤ 10⁷, normal approximation above;
  random, so checked by invariants, not parity).
- **State:** `logL` (box side, log₁₀ m), `medium` ('air' | 'air80' | 'water'), `eps` (flow's variation), `logLflow`,
  `logLbody`, `seed`, `samples` (array of {logL, reading} so far), `probe`.
  Media: air — ρ 1.225 kg/m³, M_w 28.9644, T 288.15 K, p 101 325 Pa; air80 — from `ch01.standard_atmosphere(80e3)` (T
  196.65 K, p 0.886 Pa, ρ 1.57×10⁻⁵ kg/m³; the builder copies the three numbers the Python call prints); water — ρ 998.2,
  M_w 18.015, no mean free path (Kn strip shows the molecular spacing 0.31 nm as "l").
- **Views** (one clock = the transport over `logL`):
  1. `scene` "The sample" — a square window of side 3L around the point. When the expected number of molecules in the
     drawn slice (n·(3L)²·L) is < 3000: seeded random dots (muted grey), the sampling square (accent outline) and the count
     N inside; otherwise a smooth heat-map of the macroscopic density ρ₀[1 + ε sin(2πx/L_flow)] (teal shades) with the
     square on it. Scale bar with nm/µm/mm/m label. Animates: dots jitter while playing; redraw on `seed` change. Pointer:
     click → `seed += 1` (resample). Portrait: kept (top).
  2. `curve` "Reading vs box size" — x: L [m], log, 10⁻⁹…10⁰; y: ρ_reading/ρ_point, linear 0…2 (clipped, with arrows
     for clipped samples). Faint grey band 1 ± N^{-1/2} (expected noise), orange dashed noiseless box average
     `boxAverage/ρ_point`, green shaded continuum window (noise < 10⁻³ and |drift| < 10⁻³), bold teal samples "so far",
     a vertical marker at the current L, dotted verticals at the molecular spacing and at L_flow. Pointer: click → set
     `probe` = nearest sample (inspector) and `logL`. Portrait: kept (bottom).
  3. `kn` "Knudsen strip" (`hidePortrait`) — horizontal log axis Kn 10⁻⁴…10³ with bands continuum (< 0.01, teal),
     slip (0.01–0.1, amber), transition (0.1–10, orange), free molecular (> 10, rose); markers: Kn of the body (L_body)
     and of the box (L). Portrait: hidden, hint in Explain §4.
- **Controls:** `logL` "Box side $\\log_{10} L$" −9…0 step 0.05, default −5, unit "log m", help "the side of the sampling
  cube" · `medium` select chips (sea-level air · air at 80 km · liquid water), default air · `eps` "Flow variation
  $\\epsilon$" 0…0.5 step 0.01, default 0.2, help "fractional density change over the flow scale" · `logLflow` "Flow scale
  $\\log_{10} L_{\\rm flow}$" −4…5 step 0.1, default 0, *optional* · `logLbody` "Body size $\\log_{10} L_{\\rm body}$" −9…1
  step 0.1, default −6, *optional* · `seed` 0…99 step 1, *optional*.
- **Depth features:** Explain + Code (required) · **transport** (param `logL`, −9 → 0, rate 1.2 decades/s, end 'hold',
  hold 2 s; end-of-run card: "noise fell below 0.1 % at L = 0.32 µm · the flow's variation exceeded 0.1 % at L = 6 cm ·
  continuum window ≈ 5.3 decades") · **presets** ("sea-level air" {medium: air, eps: 0.2, logLflow: 0, logL: −5} ·
  "80 km air" {medium: air80, logL: −3, logLbody: −2} · "liquid water" {medium: water, logL: −8} · "microchannel"
  {medium: air, logLbody: −5.3, logL: −6} · "weather-model grid" {logLflow: 4, logL: 0}) · **status** (noise > 10⁻³:
  "🎲 molecular noise ±x %"; else |drift| > 10⁻³: "📈 the box sees the flow's own variation (x %)"; else "✅ continuum
  plateau: noise x, drift y") · **inspector** (click a sample: "N = nL³ = 2.55×10²⁵ × (10⁻⁵)³ = 2.55×10¹⁰ · δm = Nm =
  2.55×10¹⁰ × 4.81×10⁻²⁶ = 1.225×10⁻¹⁵ kg · δV = L³ = 10⁻¹⁵ m³ · δm/δV = **1.225 kg/m³** · noise N^{-1/2} = 6.3×10⁻⁶")
  · **notes** (regime text + highlighted table of real sizes: virus 100 nm, fog/aerosol 1 µm, cloud droplet 10 µm,
  raindrop 2 mm, pipe 5 cm, weather-model cell 10 km — the row nearest L highlighted, with its N and noise) · **linked
  views** (3).
- **Readouts:** "Molecules N" · "Relative noise" (%) · "Drift" (%) · "Kn (body)".
- **Explain** (`Viz.work`):
  0. *What the three windows show* — "The **sample** window zooms on the point: grey dots are molecules (drawn only while
     there are few enough), the purple square is your box. The **reading** graph keeps every measurement: teal dots are
     samples, the grey band is the scatter statistics predicts, the orange dashed curve is what the box would read with no
     molecules at all (just the flow's smooth variation), green marks the continuum window. The **Kn strip** compares the
     mean free path with a body size."
  1. *Molecules in your box* — `n = ρA_o/M_w = 1.225 × 6.022e26 / 28.96 = 2.55e25 m⁻³` · `N = nL³ = … ` → box
     **N = …** ("why: every molecule counts once; this is the δm of ρ = δm/δV").
  2. *How noisy the reading is* — `σ/ρ = N^{-1/2} = …` → box ("D35: Poisson counting; the grey band's half-width").
  3. *How much the flow itself changes across the box* — `ρ_box/ρ_point = [1 + ε sin(2πx₀/L_flow) sinc(πL/L_flow)]/(1 +
     ε) = …` → box drift ("sinc(u) = sin u / u; a box much smaller than L_flow reads the point value").
  4. *Mean free path and Kn* — `μ = βT^{3/2}/(T + S) = …`, `ρ = p/RT = …`, `l = √(π/8) μ/(0.4987445 √(ρp)) = …`,
     `Kn = l/L_body = …` → regime name (hint on phones: "the Kn strip is hidden; turn the phone sideways").
  5. *Your continuum window* — `L_min = (10⁶/n)^{1/3} = …` (noise = 10⁻³), `L_max` from drift = 10⁻³ (solved with
     `Viz.num.brentq`) `= …`, width `log₁₀(L_max/L_min) = … decades` (box; "no window" if L_max < L_min).
  6. *At the current sweep position* — `Viz.live('logL')`, `Viz.live('N')`, last reading.
  7. *Reading the current setting* (interpret) — noise regime: "Your box holds only N molecules; one more or one fewer
     changes the reading by N^{-1/2}. No continuum value exists at this size." · plateau: "Any box in the green range
     gives the same density to better than 0.1 %: this shared number is what the equations of fluid mechanics call
     ρ(**x**). A fluid particle (C24) is a box of this size followed through the flow." · drift: "The box is now so large
     that it averages over the flow's own changes; the reading is a regional mean, not a point value." Medium add-ons:
     water ("1300× more molecules per volume: the window starts ≈ 11× smaller"), 80 km ("l ≈ 4 mm: a sensor of a few
     centimetres already sees transition-regime Kn").
- **Derivation tab:** D35 (5 steps, copied from Part F, phone-shortened). `view: 'curve'`; goal page `set: {medium: 'air',
  eps: 0, logL: -8}`; step 1 `set {logL: -8}`, `live: N̄ = nL³ = 25.5`; step 2 live `σ_N = √N̄ = 5.05`; step 3 live
  `σ_ρ = mσ_N/L³`; step 4 **live** `σ_ρ/ρ = N̄^{-1/2} = 0.198` and `watch: "the grey band's half-width at the marker is
  this number"`; step 5 `set {logL: -6}`, `watch: "one decade of L shrinks the band 10^{1.5} ≈ 32×"`. `interpret: s =>
  "With your box (L = …) the noise is … — " + (noise < 1e-3 ? "inside the continuum window." : "too noisy for a
  continuum value.")`.
- **Code** (`code`, ≤ 20 lines, live placeholders):
  ```python
  n = ch01.number_density(rho, M_w)             # molecules per m³ = {{n}}
  m = ch01.molecular_mass(M_w)                  # one molecule = {{m}} kg
  N = n * L**3                                  # expected count in the box = {{N}}
  noise = ch01.density_noise_expected(L, n)     # N**-0.5 = {{noise}}
  rho_box = ch01.box_average_density(L, rho, eps, L_flow)  # {{rhob}} kg/m³
  mean, std = ch01.sample_density(L, n, m, n_samples=200,
                                  variation=eps, L_flow=L_flow)
  l = ch01.mean_free_path_air(T, p)             # mean free path = {{l}} m
  Kn = ch01.knudsen_number(l, L_body)           # {{Kn}} -> {{regime}}
  ```
- **Walkthrough (6 steps):**
  1. "What is the density at a point?" — "Air is mostly empty space. This box is 1 nm wide: it usually holds no molecule
     at all, so its 'density' is zero — or enormous when one wanders in." `set {medium: 'air', logL: -9, eps: 0.2}`,
     highlight `view:scene`.
  2. "Grow the box" — "Press play. As the box grows the reading jumps less and less; watch the teal samples squeeze into
     the grey band." `play: true`, `controls: ['logL']`.
  3. "How noisy? One line of statistics" — "A box holding N molecules on average reads density with relative scatter
     N^{-1/2}: 0.6 % at 100 nm, six parts per million at 10 µm." `set {logL: -7}`, `eq: 'noise'`, `derive: {id: 'D35',
     step: 4}`, `readouts: ['noise']`.
  4. "The plateau ends" — "Keep growing. Near centimetres the box starts averaging over the flow's own change (orange
     dashed): the reading drifts away from the point value." `set {logL: -1, eps: 0.3}`, `notes: true`.
  5. "Kn: the same test for motion" — "A 1 µm particle in air: the mean free path is 67 nm, Kn ≈ 0.07 — slip flow, not a
     clean continuum." `set {logLbody: -6}`, `code: {id: 'e1code', lines: [7, 8]}`, `readouts: ['Kn']`.
  6. "Your turn" — "Predict first: at 80 km altitude, where does the window start, and is a 1 cm sensor in the continuum?
     Then pick the preset and sweep." `set {medium: 'air80'}`, `controls: ['logL', 'logLbody']`.
- **Equations:** `rho` "Continuum density" ref '§1.4' `\\rho(\\mathbf x) = \\delta m/\\delta V` live `δm/δV = …` ·
  `noise` "Counting noise" ref 'D35' `\\sigma_\\rho/\\rho = \\bar N^{-1/2} = (nL^3)^{-1/2}` live · `kn` "Knudsen number"
  ref '§1.4' `Kn = l/L` live · `pkin` "Pressure from impacts" ref 'D34' `p = \\tfrac13 n m\\langle|\\mathbf u|^2\\rangle
  = nk_BT` live with the medium's T (Equations tab only). Symbols: ρ kg/m³, δm kg, δV m³, n m⁻³, N̄ –, L m, l m.
- **Check yourself:** (1) "Which box side makes the noise exactly 1 % in sea-level air?" — "N̄ = 10⁴ ⇒ L = (10⁴/n)^{1/3}
  ≈ 74 nm." `set {medium: 'air', logL: -7.13}` · (2) "Switch to water at L = 10 nm. Larger or smaller noise, and by how
  much?" — "Water has ≈ 1300× more molecules per m³, so the noise is √1300 ≈ 36× smaller (0.20 → 5.5×10⁻³)." ·
  (3) "Set L_flow = 1 mm. Is there still a continuum window?" — "Yes, but narrow: from 0.3 µm to ≈ 60 µm, about 2.3
  decades." `set {logLflow: -3}` · (4) "At 80 km, what is Kn for a 1 cm sensor?" — "l ≈ 4.4 mm ⇒ Kn ≈ 0.44: transition
  regime — no continuum." `set {medium: 'air80', logLbody: -2}`.
- **Selftest parity rows:** `{name: 'noise 10 µm air', js: noiseExpected(1e-5, numberDensity(1.225, 28.9644)), py:
  'ch01.density_noise_expected(1e-5, ch01.number_density(1.225, 28.9644))', rtol: 1e-10}` · `{name: 'box average 5 cm',
  js: boxAverage(0.05, 1.225, 0.2, 1.0, 0.25), py: 'ch01.box_average_density(0.05, 1.225, 0.2, 1.0)', rtol: 1e-10}` ·
  `{name: 'mean free path 300 K', js: meanFreePathAir(300, 101325), py: 'ch01.mean_free_path_air(300.0, 101325.0)',
  rtol: 1e-8}` · invariants `{name: 'Kn', js: knudsen(6.7e-8, 1e-6), expect: 0.067, rtol: 1e-12}`, `{name: 'Poisson
  mean', js: mean of 5000 draws at λ = 25, expect: 25, rtol: 0.03}`.
- **Fit plan:** 360×640 portrait shows `scene` (top, 45 %) + `curve` (bottom) with the status line; `kn` hidden (Explain
  §4 hint); Explore shows `logL`, `medium`, `eps` (optional ones behind "more"); the real-size table lives in the notes
  pager; Explain paged in 5-line chunks; walkthrough card paged. 844×390 landscape: `curve` + `kn` side by side, scene
  small left.

### E2 · viscosity_momentum_diffusion
- **Title:** "How does the fluid learn that a plate moved?" · **Summary:** "The top plate starts moving: watch momentum
  diffuse down through the gap, the profile relax to a straight line and the wall stress settle at μU/h — with ν, not μ,
  setting the clock." · **CORE:** C12 (also C02, C08, C10, C11, C13, C14, N04, N05) · **Reference:**
  `forced_damped_vibrations.html` (system + graph on one clock, transient → steady state, "at the current time" values,
  regime-dependent interpretation).
- **meta:** `viz:sections 1.3 1.5` · `viz:equations 1.1 1.2 1.3 1.4` · `viz:fluidpy ch01.ftcs_diffusion_1d
  ch01.couette_startup_profile ch01.newton_shear_stress ch01.kinematic_viscosity ch01.thermal_diffusivity
  ch01.wall_shear_history` · `viz:derivations none`.
- **Physics:** `ftcsStep(f, r)` (Dirichlet ends) ↔ `ch01.ftcs_diffusion_1d(f0, D, dy, dt, nsteps, values=…)` ·
  `couetteStartup(y, t, U, h, nu, nterms=200)` ↔ `ch01.couette_startup_profile(y, t, U, h, nu, nterms=200)` ·
  `newtonShear(mu, dudy)` ↔ `ch01.newton_shear_stress(mu, du_dy)` · `wallShear(u, dy, mu)` (one-sided 2nd order) ↔
  `ch01.wall_shear_history(u_hist, dy, mu)` · `kinematicViscosity(mu, rho)` ↔ `ch01.kinematic_viscosity` ·
  `thermalDiffusivity(k, rho, cp)` ↔ `ch01.thermal_diffusivity` · `diffusionTime(L, D)` ↔ `ch01.diffusion_time`. Fluid
  table copied from `ch01.FLUIDS` (air, water, glycerine, honey, engine oil: ρ, μ, k, C_p, κ_m at 20 °C).
- **Grid:** 41 nodes on y/h ∈ [0, 1]; the solver runs in nondimensional time τ = t/(h²/D) with r = 0.4 (Δτ = 0.4/40²),
  so every fluid uses the same number of steps; dimensional times are τ·h²/D.
- **Modes** (`modes: {param: 'mode'}`): momentum (f = u/U, D = ν, flux τ = μ ∂u/∂y) · heat (f = (T − T_cold)/ΔT, D = κ =
  k/ρC_p, flux q = −k ∂T/∂y) · species (f = Y, D = κ_m, flux J = −ρκ_m ∂Y/∂y). The top wall steps to 1 at τ = 0.
- **Views:**
  1. `gap` "The gap" — plates (bottom fixed grey; top with an arrow U, or a hot/salty boundary in the other modes); fluid
     band shaded by f(y, t); five dyed vertical tracer lines advected by u(y, t) (momentum mode: they lean and keep
     leaning, C02); a few molecules hopping between layers carrying the colour of their layer; a dashed line AB at mid-gap
     with a flux arrow (rose) whose length ∝ |τ_AB|. Pointer: drag the top plate horizontally to change U.
  2. `profile` "u(y, t)" — x: f from 0 to 1, y: y/h from 0 to 1; ghost steady straight line (grey dashed), analytic
     series (orange thin), FTCS profile (teal bold); AB marker. Pointer: click at a height → probe (inspector).
  3. `stress` "Wall stress τ_w(t)" (`hidePortrait`) — x: τ = t/(h²/D) 0…1.5 (top axis in seconds), y: τ_w/(μU/h) 0…3;
     bottom wall teal, top wall rose (starts off-scale, arrow), faint line at 1, vertical marker at τ = 0.5 ("within 1 %
     of steady"), bold "so far".
- **Controls:** `fluid` select chips air · water · glycerine · honey · engine oil (default water) · `h` "Gap $h$" 0.1…10
  step 0.1 mm (default 1) · `U` "Plate speed $U$" 0.1…5 step 0.1 m/s (default 1) · `logmu` "Override $\\log_{10}\\mu$"
  −5…2 step 0.05 *optional* (sets a custom fluid) · `rho` "Density $\\rho$" 0.5…2000 kg/m³ *optional*.
- **Depth features:** Explain + Code · **transport** (param `tau`, 0 → 1.5, rate 0.1/s, end 'hold', hold 2; end card:
  "steady Couette flow reached at t ≈ 0.5 h²/ν = … s · τ_w = μU/h = … Pa · force on a 10 cm plate = … N") · **linked
  views** (3) · **modes** (momentum/heat/species) · **presets** (air · water · honey · "same μ as water, air's ρ"
  {logmu: −3, rho: 1.2} · "wide gap" {h: 10}) · **status** (τ < 0.5: "⏳ momentum still spreading: t/(h²/ν) = 0.12 (t =
  0.12 s)"; else "✅ steady Couette flow: τ uniform = 1.00 Pa") · **inspector** ("τ(y) = μ ∂u/∂y ≈ μ(u_{i+1} −
  u_{i−1})/(2Δy) = 1.0×10⁻³ × (0.512 − 0.488)/(2 × 2.5×10⁻⁵ m) = **0.48 Pa**") · **notes** (fluid table with the
  current row highlighted: μ, ρ, ν, h²/ν for h = 1 mm — air 1.81×10⁻⁵, 1.20, 1.51×10⁻⁵, 0.066 s; water 1.00×10⁻³, 998,
  1.00×10⁻⁶, 1.0 s; glycerine 1.42, 1264, 1.12×10⁻³, 8.9×10⁻⁴ s; honey 19.0, 1421, 1.33×10⁻², 7.5×10⁻⁵ s; engine oil
  0.80, 888, 8.96×10⁻⁴, 1.1×10⁻³ s).
- **Readouts:** "ν = μ/ρ" · "Clock h²/ν" · "t now" · "τ top wall" · "τ bottom wall".
- **Explain:**
  0. *What you see* — "The **gap** is the physical system: the top plate moves, dye lines show how the fluid deforms.
     The **u(y, t)** graph shows the velocity profile (teal) creeping toward the straight steady line (grey dashed); the
     thin orange curve is the exact series, a check on the numerics. The **τ_w(t)** graph (rose top wall, teal bottom
     wall) shows the stress each plate feels."
  1. *From μ to ν* — `ν = μ/ρ = 1.00×10⁻³/998 = ` **1.00×10⁻⁶ m²/s** ("momentum per volume is ρu, so dividing the flux
     law by ρ gives a diffusivity").
  2. *The clock* — `h²/ν = (1.0×10⁻³)²/1.00×10⁻⁶ =` **1.00 s**; `t = τ·h²/ν = …` (live).
  3. *The steady state* — `u = Uy/h`, `τ = μU/h = 1.0×10⁻³ × 1/10⁻³ =` **1.00 Pa**; plate 10 cm × 10 cm: F = 0.010 N.
  4. *The transient* — slowest mode `(2U/π) e^{−π²t/(h²/ν)} = …` m/s; within 1 % of steady when `t > (h²/ν)
     ln(200/π)/π² = 0.42 h²/ν`; half-life `ln 2/π² × h²/ν = 0.070 h²/ν`.
  5. *At the current time* — `τ_top = …`, `τ_bottom = …`, max |u − u_exact| = … (hint: "on phones the τ_w graph is
     hidden").
  6. *The other two diffusions* — heat: `κ = k/ρC_p = …`, `Pr = ν/κ = …`, `q = −k ∂T/∂y`; species: `κ_m = …`, `Sc =
     ν/κ_m = …`, `J = −ρκ_m ∂Y/∂y` (current mode's numbers boxed).
  7. *Reading the current setting* — spreading: "The plate has moved for only t/(h²/ν) = … of the clock: fluid far from it
     does not yet know. The top wall feels a large stress because the profile is steep there." · steady: "Every layer now
     passes on exactly the momentum it receives: τ is the same at every height and equals μU/h." · μ-vs-ν comparison for
     the current fluid versus water ("air: 55× smaller μ but 15× larger ν — it settles 15× faster").
- **Derivation tab:** none (the curation assigns no D row: §1.5 derives nothing). `viz:derivations none`.
- **Code:**
  ```python
  nu = ch01.kinematic_viscosity(mu, rho)          # (1.4): {{nu}} m²/s
  t_clock = ch01.diffusion_time(h, nu)            # h²/ν = {{tc}} s
  y = np.linspace(0.0, h, 41); dy = y[1] - y[0]
  dt = ch01.stable_time_step(nu, dy)              # FTCS limit: {{dt}} s
  u0 = np.zeros_like(y); u0[-1] = U               # top plate moves at t = 0
  nsteps = int(t / dt)                            # {{nsteps}} steps to t = {{t}} s
  F = ch01.ftcs_diffusion_1d(u0, nu, dy, dt, nsteps, values=(0.0, U))
  tau_b, tau_t = ch01.wall_shear_history(F, dy, mu)   # {{taub}}, {{taut}} Pa
  u_ex = ch01.couette_startup_profile(y, t, U, h, nu) # series check
  tau_steady = ch01.newton_shear_stress(mu, U / h)    # (1.3): {{tau}} Pa
  ```
- **Walkthrough (6 steps):** 1. "A plate starts moving" — "The top plate jumps to 1 m/s. Water touching it moves at once
  (no slip). How does the rest of the gap find out?" `set {fluid: 'water', tau: 0, mode: 'momentum'}` · 2. "Momentum
  diffuses" — "Play. The profile creeps down; the dye lines lean, then keep leaning — a fluid never stops deforming under
  stress." `play: true` · 3. "Stress = μ × slope" — "Early the slope at the top wall is steep, so τ there is huge. At
  steady state the line is straight and τ = μU/h everywhere." `eq: 'newton'`, `readouts: ['tauTop', 'tauBot']`,
  `inspect: true` · 4. "ν sets the clock" — "Now air: μ 55× smaller than water's, yet it settles 15× faster. The clock is
  h²/ν, and ν = μ/ρ." `set {fluid: 'air', tau: 0}`, `play: true`, `readouts: ['nu', 'clock']` · 5. "Same law, three
  things" — "Switch to heat: the profile is temperature and the flux is Fourier's q. Species gives Fick. Same maths,
  different D." `set {mode: 'heat'}`, `code: {id: 'e2code', lines: [1, 3]}` · 6. "Predict" — "Honey has a million times
  air's viscosity. In a 5 mm gap, does it settle faster or slower than water? Predict, then play." `set {fluid: 'honey',
  h: 5, mode: 'momentum', tau: 0}`, `controls: ['fluid', 'h']`.
- **Equations:** `newton` "Newton's law of viscosity" ref 'Eq. (1.3)' `\\tau = \\mu\\,\\dfrac{du}{dy}` live · `nu`
  "Kinematic viscosity" ref 'Eq. (1.4)' `\\nu \\equiv \\mu/\\rho` live · `fick` ref 'Eq. (1.1)' `\\mathbf J_m =
  -\\rho\\kappa_m\\nabla Y` · `fourier` ref 'Eq. (1.2)' `\\mathbf q = -k\\nabla T` · `diff` "Model diffusion equation
  (ours; derived in Ch. 4)" `\\partial f/\\partial t = D\\,\\partial^2 f/\\partial y^2`.
- **Check yourself:** (1) "Water, h = 1 mm vs 2 mm: how much longer to reach steady flow?" — "4× (the clock is h²/ν)." ·
  (2) "Air or water in the same gap — which settles first, and why?" — "Air: ν_air ≈ 15 ν_water although μ_air ≪ μ_water."
  · (3) "Double U. What happens to τ_w and to the settling time?" — "τ doubles (linear law); the time is unchanged
  (diffusion does not depend on U)." · (4) "Heat mode, water, 1 mm: how long, and why longer than momentum?" — "κ =
  1.44×10⁻⁷ m²/s ⇒ h²/κ ≈ 7 s: heat diffuses 7× slower than momentum in water (Prandtl number ≈ 7)." `set {mode:
  'heat', fluid: 'water', h: 1}`. (5) "Honey in 5 mm vs water in 5 mm?" — "Honey: h²/ν = 1.9×10⁻³ s; water 25 s. Honey's
  huge μ also means a huge ν."
- **Selftest parity rows:** `{name: 'ν water', js: kinematicViscosity(1.0017e-3, 998.2), py:
  'ch01.kinematic_viscosity(1.0017e-3, 998.2)', rtol: 1e-12}` · `{name: 'Couette series mid-gap', js:
  couetteStartup(0.5e-3, 0.1, 1.0, 1e-3, 1e-6, 200), py: 'ch01.couette_startup_profile(0.5e-3, 0.1, 1.0, 1e-3, 1e-6)',
  rtol: 1e-9}` · `{name: 'one FTCS step', js: ftcsStep([0, 0, 0, 1], 0.25)[2], py: 'ch01.ftcs_diffusion_1d(np.array([0.0,
  0.0, 0.0, 1.0]), 1.0, 1.0, 0.25, 1)[-1][2]', rtol: 1e-12}` · `{name: 'τ steady', js: newtonShear(1e-3, 1000), py:
  'ch01.newton_shear_stress(1e-3, 1000.0)', rtol: 1e-12}` · `{name: 'κ water', js: thermalDiffusivity(0.60304, 998.204,
  4182.11), py: 'ch01.thermal_diffusivity(0.60304, 998.204, 4182.11)', rtol: 1e-12}` · invariant `{name: 'steady τ
  uniform', js: max|τ(y) − μU/h| at τ = 1.5, expect: 0, atol: 1e-3}`.
- **Fit plan:** portrait: `gap` (top) + `profile` (bottom); `stress` hidden (Explain §5 hint); transport step/speed hidden;
  fluid table in notes pager; Explore: fluid, h, U (+ optional). Landscape phone: three views in a row with titles off.

### E3 · heat_work_paths
- **Title:** "Same two states: what depends on the path?" · **Summary:** "Take 1 kg of air between the same two states along
  different routes: heat and work change, internal energy and entropy do not — and friction breaks ds = δq/T." · **CORE:**
  C25, C35, C45 (also C26, C27, C33, C34, C44, C46) · **Reference:** `amplitude_phase_second_order_II_3.html` (several
  windows linked by one state; a numbered live derivation in the explanation) with term bars after
  `fid_formula_lab.html`.
- **meta:** `viz:sections 1.8 1.9` · `viz:equations 1.10 1.11 1.16 1.17 1.18 1.25` · `viz:fluidpy ch01.process_path
  ch01.path_heat_work_totals ch01.process_heat_work ch01.entropy_change_reversible ch01.perfect_gas_entropy_change
  ch01.isentropic_pressure ch01.irreversible_process` · `viz:derivations D10 D14`.
- **Physics:** `pathSamples(kind, v1, T1, v2, T2, n, gamma)` ↔ `ch01.process_path(kind, (v1, T1), (v2, T2), n, gamma)` ·
  `pathTotals(kind, v1, T1, v2, T2, gamma, R)` ↔ `ch01.path_heat_work_totals(kind, v1, T1, v2, T2, gamma, R)` ·
  `runningTotals(v, T, cv, R)` (cumulative trapezoid) ↔ `ch01.process_heat_work(v_path, T_path, cv, R)` ·
  `irreversible(kind, T1, v1, cv, R, T2, v2)` ↔ `ch01.irreversible_process(...)` · `entropyChange(T1, v1, T2, v2, cv, R)`
  ↔ `ch01.perfect_gas_entropy_change(...)` · `isentropicPressure(rho, p0, rho0, gamma)` ↔ `ch01.isentropic_pressure(...)`.
  Start state fixed: T₁ = 300 K, p₁ = 10⁵ Pa, v₁ = R T₁/p₁ = 0.8612 m³/kg.
- **State:** `path` · `vr` = v₂/v₁ · `Tr` = T₂/T₁ · `gas` (γ) · `proc` (mode) · `s` (progress 0…1 along the route) ·
  `probe`.
- **Modes** (`proc`): reversible · stirring (v fixed, T₁ → T₂, q = 0) · free expansion (T fixed, v₁ → v₂, q = w = 0).
- **Views:**
  1. `piston` "The gas" — cylinder with the piston at a height ∝ v; gas colour blue→orange by T; heater glyph glows orange
     while δq > 0 and turns blue (cooling) while δq < 0; hatched insulation jacket on adiabatic legs; a spinning paddle in
     stirring mode; a burst membrane and vacuum chamber in free-expansion mode; running counters q, w, Δe. Pointer: drag
     the piston to scrub `s`.
  2. `pv` "p–v diagram" — v 0.2…4 m³/kg (log x), p 10³…10⁶ Pa (log y); faint isotherms at T₁ and T₂, faint isentrope
     through state 1, all five routes as grey ghosts; the chosen route bold blue "so far", work area under the traversed
     part shaded (blue, alpha; hatched when w > 0), states 1 and 2 dots, corner marker. Irreversible modes: only the two
     states and a dashed "no equilibrium path" connector. Pointer: click on the route → `probe`.
  3. `ts` "T–s diagram" (`hidePortrait`) — s − s₁ [J kg⁻¹ K⁻¹] vs T [K]; the same route; area under it (∫T ds) shaded
     orange = q for reversible routes; isentropes are verticals.
- **Controls:** `path` select chips (isothermal · cool→expand = 'isochoric-isobaric' · expand→cool = 'isobaric-isochoric' ·
  isentrope + isochore = 'isentropic-isochoric' · isotherm + isochore = 'isothermal-isochoric') · `vr` "Volume ratio
  $v_2/v_1$" 0.25…4 step 0.05 (default 2) · `Tr` "Temperature ratio $T_2/T_1$" 0.5…2 step 0.05 (default 1) · `gas`
  select chips monatomic γ = 5/3 · diatomic γ = 7/5 (default 7/5) · `proc` via modes. Rule: 'isothermal' needs Tr = 1 (the
  chip is disabled otherwise, with a hint).
- **Depth features:** Explain + Code · **linked views** (3) · **transport** (param `s`, 0 → 1, rate 0.2/s, end 'hold'; end
  card: "route done: q = …, w = …, Δe = … kJ/kg · Δs = … J/(kg K) — compare with the ghost routes in Explain §3") · **terms**
  ("Energy so far", unit kJ/kg: q orange, w blue, Δe accent; total label "q + w = Δe") · **presets** ("isothermal ×2"
  {path: 'isothermal', vr: 2, Tr: 1, proc: 'rev'} · "cool then expand" {path: 'isochoric-isobaric', vr: 2, Tr: 1} ·
  "expand then cool" {path: 'isobaric-isochoric', vr: 2, Tr: 1} · "pump (isentropic)" {path: 'isentropic-isochoric', vr:
  0.5, Tr: 1.3195} · "stir it" {proc: 'stirring', vr: 1, Tr: 1.2} · "free expansion" {proc: 'free', vr: 2, Tr: 1}) ·
  **status** (reversible: "🔁 same Δe = … and Δs = … on every route — this one needs q = … kJ/kg"; irreversible: "⚠️
  irreversible: Δs = … > ∫δq/T = 0") · **inspector** (click the route: "p = RT/v = 287.06 × 245.3/1.290 = 5.46×10⁴ Pa ·
  w so far = −∫p dv = … kJ/kg (trapezoid over k samples) · Δe = C_v(T − T₁) = … · q = Δe − w = …").
- **Readouts:** "Heat q" · "Work w" · "Δe" · "Δs (Gibbs)" · "∫δq/T".
- **Explain:**
  0. *The three windows* — "The **gas** window is the physical system; heater colour shows the sign of δq. On the **p–v
     diagram** (blue) the area under the route is the work done by the gas; on the **T–s diagram** (orange) the area
     under the route is the heat. Grey ghosts are the other routes between the same states."
  1. *The two end states* — `v₁ = RT₁/p₁ = 287.06 × 300/10⁵ =` **0.861 m³/kg**; `v₂ = vr·v₁ = …`, `T₂ = …`, `p₂ = RT₂/v₂
     = …` (boxed).
  2. *Internal energy — a state function* — `Δe = C_v(T₂ − T₁) = … ` (boxed) ("C_v = R/(γ − 1) = …").
  3. *Work, leg by leg* — per leg the formula with numbers: isothermal `w = −RT ln(v_b/v_a)`; isochoric `w = 0`;
     isobaric `w = −p(v_b − v_a)`; isentropic `w = Δe_leg = C_v(T_b − T_a)`; total boxed; then a `Viz.work.table` of all
     five routes' (q, w, Δe, Δs) with the current row highlighted.
  4. *Heat from the first law* — `q = Δe − w = …` (1.10) per leg and total.
  5. *Entropy two ways* — Gibbs (1.18): `Δs = C_v ln(T₂/T₁) + R ln(v₂/v₁) = …`; along the route `∫δq/T` leg by leg
     (isothermal q/T; isochoric C_v ln(T_b/T_a); isobaric C_p ln(T_b/T_a); isentropic 0) = … → equal for reversible routes
     (boxed "state function"); irreversible modes: `∫δq/T = 0 < Δs` (Clausius–Duhem, taught with the actual heat).
  6. *The T–s window* — "∫T ds under the route equals q for reversible routes" (hint on phones).
  7. *At the current point* — live v, T, p, w so far, q so far.
  8. *Reading the current setting* — reversible expansion at constant T: "all the heat taken in leaves as work"; route
     with a cooling leg first: "expansion happens at low pressure, so less work and less heat"; isentropic compression:
     "no heat — the work you do stays in the gas as internal energy: T rises by …"; stirring: "adiabatic but not
     isentropic: the paddle's work warms the gas and creates entropy"; free expansion: "nothing crosses the boundary and
     yet entropy rises by R ln(v₂/v₁)".
- **Derivation tab:** **D10** (7 steps) `view: 'pv'`, goal `set {path: 'isochoric-isobaric', vr: 2, Tr: 1, proc: 'rev'}`;
  step 2 **live** "T ds = de + p dv with this route's totals: Δs = … = Δe/… (per leg)"; step 4 `highlight
  ['readout:de']`; step 7 `set {proc: 'stirring', vr: 1, Tr: 1.2}` and `watch: "Gibbs still gives Δs = 23.5 J/(kg K)
  although δq = 0"`; interpret: `s => "For your route Δs = … J/(kg K); the two Gibbs forms agree to …"`. **D14** (10 steps)
  `view: 'pv'`, goal `set {path: 'isentropic-isochoric', vr: 0.5, Tr: 1.3195, gas: '7/5', proc: 'rev'}`; step 6 **live**
  "dp/p = −γ dv/v: at the marker −1.4 × (−0.01) = +0.014"; step 9 `watch: "the bold curve lies on the faint isentrope"`;
  step 10 live "p/ρ^γ = … at state 1 and … at the marker"; interpret: `s => "Compressing to v₂/v₁ = … raises T by …
  K — the pump example of C45."`.
- **Code:**
  ```python
  cv = ch01.cv_from_gamma(gamma)                      # C_v = {{cv}} J/(kg K)
  v1 = ch01.R_AIR * T1 / p1                           # {{v1}} m³/kg
  tot = ch01.path_heat_work_totals(kind, v1, T1, v2, T2, gamma)
  # q = {{q}} kJ/kg   w = {{w}} kJ/kg   Δe = {{de}} kJ/kg
  ds = ch01.perfect_gas_entropy_change(T1, v1, T2, v2, cv=cv)   # {{ds}}
  path = ch01.process_path(kind, (v1, T1), (v2, T2), gamma=gamma)
  run = ch01.process_heat_work(path["v"], path["T"], cv=cv)     # running sums
  s_int = ch01.entropy_change_reversible(run["q"], path["T"])   # {{sint}}
  irr = ch01.irreversible_process("free_expansion", T1, v1, cv=cv, v2=v2)
  ```
- **Walkthrough (7 steps):** 1. "One kilogram, two states" — "Air at 300 K and 1 bar will end at twice the volume and the
  same temperature. Is the heat needed the same whichever way we go?" `set preset 'isothermal ×2', s: 0` · 2. "Route A:
  isothermal" — "Play. The gas expands slowly while the heater keeps it at 300 K. It takes in 59.7 kJ/kg of heat and does
  59.7 kJ/kg of work." `play: true`, `terms: true` · 3. "Route B: cool, then expand" — "Same start, same end. Now q =
  43.1 kJ/kg and w = −43.1 kJ/kg; Δe is still zero." `set preset 'cool then expand', s: 0`, `play: true`, `terms: true` ·
  4. "Something that does not depend on the route" — "Divide each bit of heat by its temperature and add: 199 J/(kg K) on
  both routes. That is entropy." `readouts: ['ds', 'sint']`, `derive: {id: 'D10', step: 2}` · 5. "Friction breaks it" —
  "Stir an insulated gas: no heat, yet T and s rise. Δs > ∫δq/T — the second law." `set preset 'stir it'`, `notes: true`
  · 6. "No heat, no friction: the isentrope" — "A fast, clean compression follows p ∝ ρ^γ, steeper than the isotherm."
  `set preset 'pump (isentropic)'`, `eq: 'isentropic'`, `derive: {id: 'D14', step: 9}` · 7. "Your turn" — "For v₂/v₁ = 3
  at constant temperature, predict which route needs the most heat, then check the table in Explain."  `set {vr: 3, Tr:
  1, proc: 'rev'}`, `controls: ['path', 'vr']`.
- **Equations:** `first` ref 'Eq. (1.10)' `\\delta q + \\delta w = \\Delta e` · `rev` ref 'Eq. (1.11)' `de = dq - p\\,dv`
  · `entropy` ref 'Eq. (1.16)' `s_2 - s_1 = \\int_1^2 dq_{\\rm rev}/T` · `gibbs` ref 'Eq. (1.18)' `T\\,ds = de + p\\,dv =
  dh - v\\,dp` · `isentropic` ref 'Eq. (1.25)' `p/\\rho^\\gamma = \\text{const}` · `cd` "Clausius–Duhem (actual heat)" ref
  '§1.8 (ii)' `s_2 - s_1 \\ge \\int_1^2 \\delta q/T`. Each with live numbers for the current route.
- **Check yourself:** (1) "Between (v₁, 300 K) and (2v₁, 300 K), which route does the least work on the surroundings?" —
  "Cool then expand: the expansion happens at half the pressure, so w = −43.1 kJ/kg." · (2) "Why is Δs = 199 J/(kg K) on
  every reversible route although q differs?" — "Heat received at low T counts more in ∫dq/T; the weighting exactly
  compensates — s is a state function." · (3) "Free expansion to 2v₁: q, w, Δe, Δs?" — "0, 0, 0 and R ln 2 = 199
  J/(kg K)." `set preset 'free expansion'` · (4) "Switch the gas to monatomic. Does the isothermal work change? Does route
  B's heat?" — "Isothermal w = −RT ln 2 has no γ: unchanged. Route B's legs involve C_v, so their separate q values
  change, but the totals q = −w = 43.1 kJ/kg do not (Δe = 0)."
- **Selftest parity rows:** `{name: 'w isothermal', js: pathTotals('isothermal', 0.8612, 300, 1.7224, 300, 1.4,
  R_AIR).w, py: 'ch01.path_heat_work_totals("isothermal", 0.8612, 300.0, 1.7224, 300.0)["w"]', rtol: 1e-9}` · `{name: 'q
  cool→expand', js: pathTotals('isochoric-isobaric', 0.8612, 300, 1.7224, 300, 1.4, R_AIR).q, py:
  'ch01.path_heat_work_totals("isochoric-isobaric", 0.8612, 300.0, 1.7224, 300.0)["q"]', rtol: 1e-9}` · `{name: 'Δs
  Gibbs', js: entropyChange(300, 0.8612, 600, 1.7224, CV_AIR, R_AIR), py: 'ch01.perfect_gas_entropy_change(300.0,
  0.8612, 600.0, 1.7224)', rtol: 1e-10}` · `{name: 'free expansion Δs', js: irreversible('free_expansion', 300, 0.8612,
  CV_AIR, R_AIR, null, 1.7224).ds, py: 'ch01.irreversible_process("free_expansion", 300.0, 0.8612, v2=1.7224)["ds"]',
  rtol: 1e-10}` · `{name: 'isentropic p', js: isentropicPressure(2.0, 1e5, 1.0, 1.4), py: 'ch01.isentropic_pressure(2.0,
  1.0e5, 1.0)', rtol: 1e-12}` · invariant `{name: 'q + w − Δe', js: …, expect: 0, atol: 1e-6}`.
- **Fit plan:** portrait: `piston` (small, top) + `pv` (bottom); `ts` hidden (Explain §6 hint); term bars in the
  walkthrough card on steps 2–3; the five-route table in the Explain pager.

### E4 · parcel_stability
- **Title:** "Push a parcel up: does it come back?" · **Summary:** "Shape the environment's temperature (or density)
  profile, lift a parcel and release it: it oscillates at N, stays put or runs away — with the stability criterion shown
  in Kundu's and in the meteorology sign convention." · **CORE:** C50, C51, C54, C55 (also C49, C52, C53, N17, C56, C57,
  C60, N21, C61) · **Reference:** `forced_damped_vibrations.html` (the system + x(t) on one clock and a regime-dependent
  "reading the current setting").
- **meta:** `viz:sections 1.10` · `viz:equations 1.29 1.30 1.31 1.32 1.35` · `viz:fluidpy ch01.adiabatic_lapse_rate
  ch01.lapse_rate_convention ch01.lapse_rate_stability ch01.brunt_vaisala_sq_from_lapse ch01.brunt_vaisala_sq
  ch01.brunt_vaisala_sq_from_theta ch01.parcel_displacement ch01.parcel_ode_atmosphere ch01.potential_temperature
  ch01.stability_timescale ch01.ocean_potential_density_gradient` · `viz:derivations D18 D19 D20 D36`.
- **Lapse-rate conventions (binding, user decision).** All physics in JS uses Kundu's `dTdz` (K/m). The `conv` chip
  (Kundu dT/dz · meteorology −dT/dz) changes **only** which form is primary in the status line, readouts, Explain §5
  and the Equations tab; both forms are always visible (notes table). The criterion is always printed with its
  inequality: "stable ⇔ dT/dz > Γa: −6.5 > −9.8 K/km ✓" (Kundu) or "stable ⇔ Γ < Γa: 6.5 < 9.8 K/km ✓" (meteorology),
  from `lapseStability` = mirror of `ch01.lapse_rate_stability` (text format identical, Unicode minus; the status line
  uses the bare relation, `prefix=False`).
- **Physics:** `adiabaticLapseRate(cp, g)` ↔ `ch01.adiabatic_lapse_rate(cp=…, g=…)` · `lapseConvention(dTdz, conv)` ↔
  `ch01.lapse_rate_convention(dT_dz, convention)` · `lapseStability(dTdz, Ga, conv)` → {verdict, text, code, Gamma,
  Gamma_a, margin} ↔ `ch01.lapse_rate_stability(dT_dz, Gamma_a, convention)` · `N2FromLapse(T, dTdz, cp, g)` ↔
  `ch01.brunt_vaisala_sq_from_lapse` · `N2FromTheta(theta, dthdz, g)` ↔ `ch01.brunt_vaisala_sq_from_theta` ·
  `N2Density(rho0, drhodz, drhoadz, g)` ↔ `ch01.brunt_vaisala_sq` · `stabilityTimescale(N2)` ↔
  `ch01.stability_timescale` · `parcelDisplacement(t, z0, N2, w0)` ↔ `ch01.parcel_displacement` · `parcelAccelAtm(zeta,
  T0, dTdz, cp, g)` + `Viz.num.rk4Step` integration ↔ `ch01.parcel_acceleration_atmosphere` / `ch01.parcel_ode_atmosphere`
  · `linearLapsePressure(z, p0, T0, dTdz, R, g)` ↔ `ch01.linear_lapse_pressure` · `potentialTemperature(T, p, pref,
  gamma)` ↔ `ch01.potential_temperature` · `isentropicDensityGradient(rho, c, g)` ↔ `ch01.isentropic_density_gradient`
  · `oceanCriterion(drhodz, rho, c, g)` ↔ `ch01.ocean_potential_density_gradient` · `seawaterLinear(T, S)` ↔
  `ch01.seawater_density_linear`.
- **Modes** (`medium`): **atmosphere** (T profile, dry parcel, 0–3 km) · **ocean** (ρ profile from a thermocline dT/dz at
  S = 35 g/kg with the linear EOS, 0–300 m, optional parcel compressibility with c = 1500 m/s) · **lab tank** (salt
  stratification dρ/dz, incompressible parcel, 0–0.5 m).
- **Views:**
  1. `column` "The column" (`hidePortrait`) — vertical strip coloured by θ (atmosphere) or ρ (ocean/lab); the parcel (teal
     disc) at z_o + ζ(t); a force arrow on it (teal pointing back = restoring, rose pointing away = runaway); dashed rest
     level z_o; inversion layer outlined when on. Pointer: drag the parcel vertically to set ζ₀ (pauses the clock).
  2. `profile` "T(z) and the parcel's adiabat" — atmosphere: environment T(z) bold black (with the inversion kink when
     on), dry adiabat through (z_o, T(z_o)) dashed purple, θ(z) faint purple when `showTheta`, the parcel as a teal dot at
     (T_parcel, z_o + ζ) with a horizontal bracket to the environment at that height labelled "ΔT = −3.3 K → sinks";
     ocean/lab: ρ(z) and the parcel's ρ_a(z). Pointer: click a height → inspector; vertical drag on the environment line
     changes `dTdz`.
  3. `zeta` "ζ(t)" — x: t 0…t_end (min), y: ζ (m); linear solution ghost (grey: cos / line / cosh), nonlinear RK4 path
     bold so far (teal stable, amber neutral, rose unstable), ticks at each period or e-folding time; clipped at ±ζ cap
     with an arrow. Portrait: kept.
- **Controls:** `dTdz` "Environment $dT/dz$" −15…+10 step 0.1 K/km (default −6.5; help "Kundu's Γ; the meteorology lapse
  rate is its negative"; label shows the live meteorology value) · `conv` select chips "Kundu Γ = dT/dz" · "meteorology Γ
  = −dT/dz" · `zeta0` "Release height $\\zeta_0$" −300…300 step 10 m (default 100) · `inv` toggle "Inversion layer 800–1000
  m" · `invRate` "Inversion $dT/dz$" 0…+30 K/km *optional* · `showTheta` toggle *optional* · ocean/lab: `drhodz`
  "Density gradient $d\\rho/dz$" −0.05…+0.02 kg m⁻⁴ (lab: −50…+5) · `compress` toggle "Parcel compressible (c = 1500
  m/s)" *optional*.
- **Depth features:** Explain + Code · **linked views** (3) · **transport** (param `t`, 0 → t_end = 3 periods or 5
  e-foldings (≤ 60 min), rate scaled so a run lasts 12 s, end 'hold'; end card: "3 oscillations, period 9.9 min —
  undamped in this model (real air damps it)" or "ran away: ζ reached 1 km after 6.2 min") · **presets** ("dry adiabatic
  (neutral)" {medium: 'atm', dTdz: −9.76, inv: false} · "standard atmosphere" {dTdz: −6.5} · "isothermal layer" {dTdz: 0}
  · "nocturnal inversion" {inv: true, invRate: 10, dTdz: −6.5, zeta0: 50, zo: 900} · "superadiabatic ground layer"
  {dTdz: −12} · "ocean thermocline" {medium: 'ocean', drhodz: −0.01, compress: false} · "lab salt tank" {medium: 'lab',
  drhodz: −10}) · **status** ("🌊 stable: −6.5 > −9.8 K/km · period 9.9 min" / "⚖️ neutral: dT/dz = Γa" / "🚀 unstable:
  −12.0 < −9.8 K/km · e-folds in 1.9 min"; the form follows `conv`) · **modes** (atmosphere / ocean / lab) · **notes**
  (always-visible convention table, current convention row highlighted: | Convention | Γ now | Γa | stable when | verdict
  |, e.g. | Kundu Γ = dT/dz | −6.5 | −9.8 | Γ > Γa | −6.5 > −9.8 ✓ | · | meteorology Γ = −dT/dz | +6.5 | +9.8 | Γ < Γa |
  6.5 < 9.8 ✓ |) · **inspector** (click a height z: "T_env = T₀ + (dT/dz)(z − z_o) = 288.15 − 6.5 × 0.1 = 287.50 K · T_p
  = T₀ + Γa(z − z_o) = 288.15 − 9.76 × 0.1 = 287.17 K · a = g(T_p − T_e)/T_e = 9.81 × (−0.33)/287.5 = **−0.0111 m/s²**
  (restoring)").
- **Readouts** (both conventions in every lapse-rate readout): "Γ env" → "−6.5 K/km (met +6.5)" · "Γa" → "−9.8 K/km (met
  +9.8)" · "N²" · "Period / e-fold" · "dθ/dz".
- **Explain:**
  0. *What the windows show* — "The **column** is the air (or water) coloured by potential temperature (or density); the
     teal disc is your parcel. In **T(z)** the black line is the environment and the dashed purple line is the parcel's
     own dry adiabat through its starting point; the teal dot is the parcel now. **ζ(t)** shows its height: grey is the
     linear solution, the coloured curve the full nonlinear motion."
  1. *The parcel's own cooling rate* — `Γa = −g/C_p = −9.807/1004.7 = ` **−9.76×10⁻³ K/m = −9.76 K/km (Kundu)**; `Γa,met =
     −Γa =` **+9.76 K/km** ("a rising parcel expands as the pressure drops and cools at this fixed rate; α = 1/T for a
     perfect gas turns (1.30) into −g/C_p").
  2. *Your environment's lapse rate* — `dT/dz = …` (Kundu) · `Γ_met = −dT/dz = …` (meteorology) (both boxed).
  3. *N² and its time scale* — `N² = (g/T)(dT/dz − Γa) = (9.807/288.15)(−6.5 + 9.76)×10⁻³ = ` **1.11×10⁻⁴ s⁻²** →
     `N = 0.0105 s⁻¹`, `period 2π/N = 596 s = 9.9 min` (or `σ = √(−N²)`, e-folding `1/σ`) (C51, D36).
  4. *Potential temperature* — `θ(z_o) = T(p_o/p)^{R/C_p} = …`, `dθ/dz = (θ/T)(dT/dz − Γa) = …` K/km (1.32): positive ⇒
     θ increases upward.
  5. *The criterion in both conventions* — `Viz.work.table(['convention', 'Γ', 'Γa', 'stable when', 'now'], [['Kundu (dT/dz)',
     '−6.5', '−9.8', 'Γ > Γa', '−6.5 > −9.8 ✓'], ['meteorology (−dT/dz)', '+6.5', '+9.8', 'Γ < Γa', '6.5 < 9.8 ✓']],
     current)` + say "Negating both numbers reverses the inequality; the margin Γ − Γa (Kundu) = Γa,met − Γ_met = 3.26
     K/km is the same, so the verdict never depends on the convention."
  6. *The parcel now* — live t, ζ (nonlinear), ζ (linear), T_parcel, T_env, acceleration.
  7. *Ocean and lab mode* (or hint "switch the mode to see the density form") — `dρ_a/dz = −ρg/c² = −1025 × 9.81/1500²
     = −4.47×10⁻³ kg m⁻⁴`; `dρ/dz + ρg/c² = …` "has the same sign as dρ_θ/dz" (1.35); `N²` with and without
     compressibility.
  8. *Reading the current setting* — stable: "Your environment cools more slowly than a rising parcel (… K/km vs 9.76), so
     a lifted parcel ends up colder and denser than its surroundings and sinks back, overshoots and oscillates with period
     … — a buoyancy (Brunt–Väisälä) oscillation. Real air damps it by mixing." · neutral: "The environment is on the
     adiabat: a parcel is at home at every height — the well-mixed boundary layer." · unstable: "The environment cools
     faster than the parcel: a lifted parcel is warmer and lighter than its surroundings and accelerates away, e-folding
     every … s — convection, thermals, towering cumulus." · inversion: "Inside the inversion T rises with height: very
     stable (period …): smoke and pollution below cannot mix up through it." · climate note (always): "Lapse-rate-feedback
     papers write Γ = −dT/dz; a warming that makes their Γ smaller (closer to moist adiabatic) makes the upper troposphere
     warm more than the surface."
- **Derivation tab:**
  - **D18** (12 steps) `view: 'zeta'`; goal `set {medium: 'ocean', drhodz: -0.01, compress: false, zeta0: 5}`; step 2
    `set {t: 0}` highlight `view:column`; step 6 **live** "ρ_p − ρ_env = (dρ_a/dz − dρ/dz)ζ = (0 + 0.01) × 5 = 0.05 kg/m³";
    step 10 live "N² = −(9.81/1025)(−0.01 − 0) = 9.57×10⁻⁵ s⁻²"; step 12 `set {t: 0}` + play, `watch: "the grey cos curve
    has period 2π/N = 10.7 min"`; interpret `s => "With your gradient N² = …: " + regimeText`.
  - **D19** (16 steps, ★★★) `view: 'profile'`; goal `set {medium: 'atm', dTdz: -6.5, conv: 'kundu'}`; steps 1–12 no set;
    step 13 `highlight ['readout:Ga']`; step 15 **live** "Γa = −gαT/C_p = −9.807 × (1/288.15) × 288.15/1004.7 = −9.76×10⁻³
    K/m"; step 16 `watch: "the dashed adiabat's slope is exactly this"`; result page shows both conventions; interpret `s
    => "Your environment: dT/dz = … vs Γa = −9.76 (Kundu) ⇔ Γ_met = … vs +9.76 (meteorology) → …"`.
  - **D20** (5 steps) `view: 'profile'`, `set {showTheta: true}`; step 2 live "T = θ (p/p_o)^{0.2857} at z_o"; step 5
    live "θ = …".
  - **D36** (8 steps) `view: 'profile'`, `set {showTheta: true, dTdz: -6.5}`; step 7 **live** "N² = (g/θ)dθ/dz = (9.807/…)
    × … = …" and `watch: "equals the N² readout computed from the lapse rate"`; step 8 live "(g/T)(Γ − Γa) = …".
- **Code:**
  ```python
  Ga = ch01.adiabatic_lapse_rate()                        # Kundu: {{Ga}} K/m
  Ga_met = ch01.lapse_rate_convention(Ga, "meteorology")  # {{Gam}} K/m
  res = ch01.lapse_rate_stability(dT_dz, convention=conv)
  print(res.verdict, res.text)            # {{verdict}}: {{text}}
  N2 = ch01.brunt_vaisala_sq_from_lapse(T0, dT_dz)        # {{N2}} 1/s²
  kind, tau = ch01.stability_timescale(N2)                # {{kind}}, {{tau}} s
  zeta = ch01.parcel_ode_atmosphere(t, zeta0, T0, dT_dz)  # nonlinear: {{zeta}} m
  zeta_lin = ch01.parcel_displacement(t, zeta0, N2)       # linear: {{zlin}} m
  theta = ch01.potential_temperature(T0, p0)              # θ = {{theta}} K
  ```
- **Walkthrough (8 steps):** 1. "Colder aloft — unstable?" — "The standard atmosphere cools 6.5 K per km. Colder air above
  warmer air sounds unstable. Is it?" `set preset 'standard atmosphere'`, highlight `view:profile` · 2. "The parcel cools
  too" — "Lift a parcel: it expands and cools at its own fixed rate, 9.76 K/km — the dashed line." `eq: 'Ga'`, `derive:
  {id: 'D19', step: 15}` · 3. "Compare the two slopes" — "At 100 m up the parcel is 0.33 K colder than the air around it,
  so it sinks back. Play." `play: true`, `inspect: true` · 4. "N² and the period" — "The restoring force grows with ζ: an
  oscillator with N² = 1.11×10⁻⁴ s⁻², period 9.9 min." `eq: 'N2'`, `readouts: ['N2', 'tscale']`, `derive: {id: 'D18',
  step: 10}` · 5. "Two sign conventions, one verdict" — "Flip to meteorology: Γ = +6.5 < Γa = +9.8. Kundu wrote −6.5 >
  −9.8. Negate both numbers, flip the inequality; nothing physical changes." `set {conv: 'met'}`, `notes: true` · 6. "Runaway"
  — "Drag dT/dz to −12 K/km: the environment now cools faster than the parcel. The cosine becomes a cosh." `set {dTdz:
  -12, conv: 'kundu', t: 0}`, `play: true` · 7. "θ says it in one look" — "Potential temperature removes the parcel's
  cooling: stable exactly where θ increases upward, and N² = (g/θ)dθ/dz." `set {dTdz: -6.5, showTheta: true}`, `derive:
  {id: 'D36', step: 7}` · 8. "Your turn: the ocean" — "Predict: does letting the parcel feel seawater's compressibility
  make the thermocline more or less stable? Toggle and read the period." `set preset 'ocean thermocline'`, `controls:
  ['compress', 'drhodz']`.
- **Equations:** `parcel` "Parcel equation" ref '§1.10' `\\ddot\\zeta + N^2\\zeta = 0` · `N2` ref 'Eq. (1.29)' `N^2 =
  -\\dfrac{g}{\\rho(z_o)}\\Big(\\dfrac{d\\rho}{dz} - \\dfrac{d\\rho_a}{dz}\\Big)` · `Ga` ref 'Eq. (1.30)' `\\Gamma_a =
  -g\\alpha T/C_p` with note "Kundu: Γ ≡ dT/dz, Γa = −9.76 K/km, stable ⇔ Γ > Γa · meteorology: Γ ≡ −dT/dz, Γa = +9.76
  K/km, stable ⇔ Γ < Γa" · `theta` ref 'Eq. (1.31)' `T = \\theta\\,(p/p_o)^{(\\gamma-1)/\\gamma}` · `dtheta` ref 'Eq.
  (1.32)' `\\dfrac{T}{\\theta}\\dfrac{d\\theta}{dz} = \\Gamma - \\Gamma_a` (+ meteorology line `= \\Gamma_{a,\\rm met} -
  \\Gamma_{\\rm met}`) · `N2theta` ref 'D36' `N^2 = \\dfrac{g}{\\theta}\\dfrac{d\\theta}{dz}` · `ocean` ref 'Eq. (1.35)'
  (ocean mode) `\\dfrac{d\\rho_\\theta}{dz} \\sim \\dfrac{d\\rho}{dz} + \\dfrac{\\rho g}{c^2}`.
- **Check yourself:** (1) "Set dT/dz = 0 (isothermal) at 250 K. Stable? What period?" — "0 > −9.76 ⇒ stable; N² = g²/(C_pT)
  = 3.83×10⁻⁴ s⁻², period 5.4 min." `set {dTdz: 0}` · (2) "A meteorology paper reports Γ = 11 K/km. Enter it. Verdict?"
  — "dT/dz = −11 K/km < −9.76 ⇒ unstable (in the paper's form 11 > 9.76)." · (3) "Which dT/dz gives a 20-minute period at
  288 K?" — "N = 2π/1200 = 5.24×10⁻³ s⁻¹, N² = 2.74×10⁻⁵ = (g/T)(dT/dz + 9.76×10⁻³) ⇒ dT/dz ≈ −8.96 K/km (Γ_met ≈ 8.96)."
  · (4) "Ocean thermocline: toggle compressibility. What happens to N²?" — "It drops from 9.57×10⁻⁵ to 5.29×10⁻⁵ s⁻²
  (period 10.7 → 14.4 min): the parcel also compresses as it sinks, so it is less heavy than an incompressible one."
- **Selftest parity rows:** `{name: 'Γa Kundu', js: adiabaticLapseRate(CP_AIR, G0), py: 'ch01.adiabatic_lapse_rate()',
  rtol: 1e-12}` · `{name: 'Γ meteorology', js: lapseConvention(-0.0065, 'meteorology'), py:
  'ch01.lapse_rate_convention(-0.0065, "meteorology")', rtol: 1e-12}` · `{name: 'stability code (met)', js:
  lapseStability(-0.012, null, 'meteorology').code, py: 'ch01.lapse_rate_stability(-0.012, convention="meteorology")[2]',
  atol: 0.5}` · `{name: 'margin', js: lapseStability(-0.0065, null, 'kundu').margin, py:
  'ch01.lapse_rate_stability(-0.0065)[5]', rtol: 1e-9}` · `{name: 'N² from lapse', js: N2FromLapse(288.15, -0.0065,
  CP_AIR, G0), py: 'ch01.brunt_vaisala_sq_from_lapse(288.15, -0.0065)', rtol: 1e-12}` · `{name: 'ζ linear', js:
  parcelDisplacement(300, 100, 1.1097e-4, 0), py: 'ch01.parcel_displacement(300.0, 100.0, 1.1097e-4)', rtol: 1e-12}` ·
  `{name: 'ζ nonlinear', js: rk4 run to 300 s, py: 'ch01.parcel_ode_atmosphere(300.0, 50.0, 288.15, -0.0065)', rtol:
  1e-5}` · `{name: 'θ', js: potentialTemperature(250, 5e4, 1e5, 1.4), py: 'ch01.potential_temperature(250.0, 50000.0)',
  rtol: 1e-12}` · `{name: 'period', js: stabilityTimescale(1e-4).seconds, py: 'ch01.stability_timescale(1e-4)[1]', rtol:
  1e-12}` · `{name: 'ocean criterion', js: oceanCriterion(-0.01, 1025, 1500, G0), py:
  'ch01.ocean_potential_density_gradient(-0.01, 1025.0, 1500.0)', rtol: 1e-12}` · invariant `{name: 'verdict
  convention-free', js: codes equal for dTdz in −15…10, expect: 1, atol: 0}`.
- **Fit plan:** portrait: `profile` (top) + `zeta` (bottom) with the status line (which carries the inequality); `column`
  hidden (its parcel is also drawn on `profile`; Explain §0 hint); the convention table is the notes panel (always one
  tap away); Explore shows dTdz, conv, zeta0, inv (+ optional). Derivation tab on phones keeps the view named per D.
  Landscape: three views in a row.

### E5 · buckingham_pi_machine
- **Title:** "Why can 7 pipe variables shrink to 4 numbers?" · **Summary:** "Switch variables on and off, pick the
  repeating set and watch the dimensional matrix, its rank, the exponent solve and the Π groups change — then change units
  and see every Π keep its value." · **CORE:** C64, C67, C69 (also C65, C66, C68, N23, C70, C71, C72, C73, C75, C76) ·
  **Reference:** `amplitude_phase_second_order_II_3.html` (numbered live derivation) with presets and "formula with
  numbers plugged in" after `stride_padding_playground.html`.
- **meta:** `viz:sections 1.11` · `viz:equations 1.36 1.37 1.38 1.39 1.40` · `viz:fluidpy ch01.dimensional_matrix
  ch01.rank_by_minors ch01.minor_determinant ch01.solve_exponents ch01.pi_groups ch01.groups_independent ch01.group_value
  ch01.rescale_units` · `viz:derivations D28`.
- **Physics:** a small exact `Frac` class (add, sub, mul, div, normalise) — promotion candidate for `viz_lib.js` ·
  `DIM` table of the preset units → exponent vectors (M, L, T, Θ) mirroring `ch01.dimension_vector` · `dimMatrix(vars)`
  ↔ `ch01.dimensional_matrix` · `minorDet(A, rows, cols)` ↔ `ch01.minor_determinant` · `rankByMinors(A)` ↔
  `ch01.rank_by_minors` (same lexicographic witness) · `solveExponents(target, rep, vars)` (Gauss–Jordan on Frac) ↔
  `ch01.solve_exponents` · `piGroups(vars, solution, rep)` ↔ `ch01.pi_groups` · `groupValue(g, vals)` ↔
  `ch01.group_value` · `rescale(vals, vars, system)` ↔ `ch01.rescale_units` (SI → cgs: M ×10³, L ×10²; imperial: lb, ft,
  s, °R) · `independent(groups, names)` ↔ `ch01.groups_independent`.
- **Presets (modes `problem`):** pipe Δp (`PIPE`, repeating U, d, ρ) · pendulum period (`PENDULUM`) · sphere drag
  (`SPHERE_DRAG`) · scale height Ex. 1.2 (`SCALE_HEIGHT`) · blast Ex. 1.4 (`BLAST`) · Rayleigh Ex. 1.5 (`RAYLEIGH`). Each
  with default SI values for the value checks (pipe: Δp 100 Pa, Δx 1 m, d 0.01 m, ε 10⁻⁵ m, U 0.1 m/s, ρ 1000, μ 10⁻³).
- **Views:**
  1. `vars` "Variables" (`hidePortrait`) — one chip per variable: symbol, SI unit, a tiny M/L/T/Θ exponent bar group
     (orange M, blue L, teal T, amber Θ); on/off (click toggles), repeating ones ringed purple, the solution variable
     starred.
  2. `matrix` "Dimensional matrix" — rows M, L, T (Θ when present) × active columns; integers in cells (blue negative,
     orange positive, grey zero); repeating columns highlighted; the r×r repeating minor outlined with a determinant badge
     (green nonzero, red zero); big "r = 3 · n − r = 4". Pointer: click a cell → inspector; click a column header → toggle
     it as repeating.
  3. `groups` "Π groups" — each group as a KaTeX fraction; beneath it its dimension vector (0, 0, 0, 0 ✓) and its value in
     SI | cgs | imperial (identical numbers); independence badge; "Π₁ = φ(Π₂, Π₃, Π₄)" line.
- **Controls:** `problem` select (6 presets) · `rep` "Repeating set" select (options per problem; pipe: "U, d, ρ (book)",
  "U, d, μ", "Δx, U, ρ", "U, ρ, μ", "d, ε, ρ (singular)") · `drop` "Remove a variable" select (none, μ, ε, Δx) · `units`
  chips SI · cgs · imperial.
- **Depth features:** Explain + Code · **presets** (the six problems) · **inspector** (matrix cell: "[μ] = kg m⁻¹ s⁻¹ =
  M¹L⁻¹T⁻¹ → row L, column μ = −1"; group: "Π₁ = Δp U^a d^b ρ^c: M: c + 1 = 0 → c = −1 · T: −a − 2 = 0 → a = −2 · L: a + b −
  3c − 1 = 0 → b = 0") · **status** ("✅ n = 7, r = 3 → 4 independent groups" / "⚠️ singular repeating set: det = 0 —
  the chosen variables cannot cancel T" / "⚠️ n_s is dimensionless by itself → it is its own group" / "⚠️ μ removed: no
  group can describe viscous friction") · **linked views** (3) · **notes** (the physics step 7 per problem: pipe Δp ∝ Δx
  → Re; blast → D ∝ t^{2/5}; Rayleigh → λ⁻⁴).
- **Readouts:** "Variables n" · "Rank r" · "Groups n − r" · "Minor det".
- **Explain:**
  0. *The three windows* — "Chips are the variables with their dimension exponents; the **matrix** stacks them as columns;
     the **groups** are what the matrix allows. Orange = positive exponent, blue = negative; purple = repeating."
  1. *Your variable list* — n = …; a table of [q] for each.
  2. *Rank by minors* — the witness minor printed as a 3×3 array with the cofactor expansion `0·(…) − 1·(…) + 1·(…) = −1`
     → **r = 3** ("the largest square block with a nonzero determinant").
  3. *How many groups* — `n − r = 7 − 3 =` **4** (rank–nullity, D28 step 5).
  4. *Each group* — for every non-repeating variable: the exponent equations M, L, T with numbers, their solution, the
     group (boxed).
  5. *The check* — each group's dimension vector (0, 0, 0) and its value in SI, cgs and imperial (identical) — dimensional
     homogeneity (C64).
  6. *The dimensionless law* — Π₁ = φ(Π₂, …) for the current set; step 7 physics note.
  7. *Reading the current setting* — pipe: "Π₄ = μ/ρUd = 1/Re: viscosity matters only through the Reynolds number" ·
     other repeating set: "a different but equivalent basis — each new group is a product of powers of the old ones" ·
     singular: "d and ε are both pure lengths, so no power of them can cancel the T in Δp" · dropped μ: "3 groups remain;
     the law can no longer know about viscosity — laminar data would not collapse" · blast: "one group, so it is a
     constant: E ∝ ρD⁵/t²" · Rayleigh: "the T row is −3 × the M row, so r = 2 and 4 groups; physics fixes the powers of
     V and d, leaving λ⁻⁴".
- **Derivation tab:** **D28** (12 steps, ★★★) `view: 'matrix'`; goal `set {problem: 'pipe', rep: 'Udrho', drop: 'none',
  units: 'SI'}`; step 3 **live** "A·k for Π₁'s k = (1, 0, 0, 0, −2, −1, 0): M 1 − 1 = 0, L −1 − 2 + 3 = 0, T −2 + 2 = 0";
  step 5 live "7 − 3 = 4"; step 8 `highlight ['view:matrix']` (the repeating minor); step 10 `set {units: 'cgs'}` `watch:
  "every Π value in the groups window is unchanged"`; step 12 live the law for the current set; interpret `s => "Your
  problem: n = …, r = …, so … groups: …"`; `check: 'sympy nullspace of (1.39) has 4 vectors; parity rows below'`.
- **Code:**
  ```python
  A, names, rows = ch01.dimensional_matrix(ch01.PIPE)    # {{rows}} × {{names}}
  r, ri, ci = ch01.rank_by_minors(A)                     # r = {{r}}
  det = ch01.minor_determinant(A, (0, 1, 2), cols)       # {{det}}
  groups = ch01.pi_groups(ch01.PIPE, solution="dp",
                          repeating=rep)                  # {{n_groups}} groups
  tex = [ch01.group_latex(g) for g in groups]             # {{tex}}
  new = ch01.rescale_units(values, ch01.PIPE, units)
  vals = [ch01.group_value(g, new) for g in groups]      # {{vals}}
  ok = ch01.groups_independent(groups)                   # {{ok}}
  ```
- **Walkthrough (7 steps):** 1. "Seven variables" — "The pressure drop along a pipe depends on Δx, d, ε, U, ρ and μ. Seven
  variables. Could a law really need all seven separately?" `set preset pipe` · 2. "Dimensions as columns" — "Each chip's
  exponents become a column. Seven columns, three rows: the dimensional matrix (1.39)." `eq: 'matrix'`, highlight
  `view:matrix` · 3. "Rank 3" — "The (U, ρ, μ) block has determinant −1 ≠ 0, so three dimensions are independent: r = 3."
  `inspect: true`, `readouts: ['rank', 'det']` · 4. "Null space = groups" — "A product of powers is dimensionless when its
  exponents solve A·k = 0. There are n − r = 4 independent solutions." `derive: {id: 'D28', step: 5}` · 5. "Units don't
  matter" — "Switch to imperial: every number in the table changes, every Π stays." `set {units: 'imperial'}` · 6. "A bad
  choice" — "Pick (d, ε, ρ) as repeating: two pure lengths cannot cancel time. The machine refuses." `set {rep:
  'depsrho'}`, `notes: true` · 7. "Your turn" — "Remove μ. Predict how many groups remain and what physics is lost, then
  check. Try the blast preset next." `set {rep: 'Udrho', drop: 'mu'}`, `controls: ['drop', 'problem']`.
- **Equations:** `f` ref 'Eq. (1.36)' `f(q_1,\\dots,q_n) = 0` · `pi` ref 'Eq. (1.37)' `\\phi(\\Pi_1,\\dots,\\Pi_{n-r}) = 0`
  · `pipe` ref 'Eq. (1.38)' `f(\\Delta p, \\Delta x, d, \\varepsilon, U, \\rho, \\mu) = 0` · `matrix` ref 'Eq. (1.39)' the
  matrix as `\\begin{array}` · `law` ref 'Eq. (1.40)' `\\Delta p/\\rho U^2 = \\varphi(\\Delta x/d, \\varepsilon/d,
  \\mu/\\rho U d)` · `expo` "Exponents of Π₁" `M^{c+1}L^{a+b-3c-1}T^{-a-2} = M^0L^0T^0` live.
- **Check yourself:** (1) "Choose repeating (d, ε, ρ). Why does it fail?" — "d and ε are both lengths; the 3×3 minor has an
  all-zero T row, det = 0, so time in Δp and μ cannot be cancelled." · (2) "Remove μ: how many groups, and what is lost?" —
  "6 − 3 = 3; no Reynolds number, so viscous (e.g. laminar) behaviour cannot be described." · (3) "Blast preset: how does
  the radius grow with time?" — "One group Et²/(ρD⁵) = K ⇒ D ∝ t^{2/5}." · (4) "Scale-height preset: why r = 4 and what is
  the constant?" — "Θ appears (T₀, R_u) so four dimensions; one group ⇒ H = const·R_uT₀/(gM_w), const = 1 from the
  isothermal atmosphere."
- **Selftest parity rows:** `{name: 'rank pipe', js: rankByMinors(dimMatrix(PIPE)).r, py:
  'ch01.rank_by_minors(ch01.dimensional_matrix(ch01.PIPE)[0])[0]', atol: 1e-12}` · `{name: 'minor (U, ρ, μ)', js:
  minorDet(dimMatrix(PIPE), [0, 1, 2], [4, 5, 6]).valueOf(), py: 'ch01.minor_determinant(ch01.dimensional_matrix(
  ch01.PIPE)[0], (0, 1, 2), (4, 5, 6))', atol: 1e-12}` · `{name: 'exponent of U in Π₁', js: solveExponents('dp', ['U',
  'd', 'rho'], PIPE).U.valueOf(), py: 'ch01.solve_exponents("dp", ("U", "d", "rho"), ch01.PIPE)["U"]', atol: 1e-12}` ·
  `{name: 'groups n − r', js: 4, py: 'ch01.dimensional_matrix(ch01.PIPE)[0].shape[1] -
  ch01.rank_by_minors(ch01.dimensional_matrix(ch01.PIPE)[0])[0]', atol: 1e-12}` · `{name: 'Π₁ in cgs', js:
  groupValue({dp: 1, U: -2, rho: -1}, rescale(PIPE_VALUES, PIPE, 'cgs')), py: 'ch01.group_value({"dp": 1, "U": -2, "rho":
  -1}, ch01.rescale_units({"dp": 100.0, "dx": 1.0, "d": 0.01, "eps": 1e-05, "U": 0.1, "rho": 1000.0, "mu": 0.001},
  ch01.PIPE, "cgs"))', rtol: 1e-12}` · `{name: 'rank Rayleigh', py: 'ch01.rank_by_minors(ch01.dimensional_matrix(
  ch01.RAYLEIGH)[0])[0]', atol: 1e-12}` (expected 2).
- **Fit plan:** portrait: `matrix` (top) + `groups` (bottom); `vars` hidden (the Explore select + matrix header clicks
  replace it; Explain §1 lists the variables); matrix cell text ≥ 12 px means at most 8 columns — every preset fits;
  landscape: three views.

### B1 · pressure_in_still_fluid
*(Backup — built only if one of E1–E5 fails review; then add D05, D37 to those rows' "Shown in".)*
- **Title:** "Why does still water push up on a block?" · **Summary:** "Drag a block through a layered tank and watch the
  face pressures: sideways forces cancel, top and bottom differ by exactly the weight of displaced fluid; in a thin tube
  the curved meniscus lifts the water." · **CORE:** C20 (also C16, C17, C19, C21, C22, N07) · **Reference:**
  `angular_frequency_explorer_1.html` (modes, presets, Right-now notes).
- **meta:** `viz:sections 1.6 1.7` · `viz:equations 1.5 1.7 1.8 1.9` · `viz:fluidpy ch01.hydrostatic_pressure_uniform
  ch01.layered_pressure ch01.buoyancy_force ch01.net_pressure_force_on_box ch01.gauge_pressure ch01.capillary_rise
  ch01.laplace_pressure_jump` · `viz:derivations D05 D37`.
- **Physics:** `hydroUniform(z, p0, rho, g)` ↔ `ch01.hydrostatic_pressure_uniform` · `layered(z, h, rho, pS, g)` ↔
  `ch01.layered_pressure` · `buoyancy(rho, V, g)` ↔ `ch01.buoyancy_force` · `boxFaceForces(pfn, box)` (exact for linear p)
  ↔ `ch01.net_pressure_force_on_box` · `gauge(p)` ↔ `ch01.gauge_pressure` · `capillaryRise(sigma, alpha, rho, R, g)` ↔
  `ch01.capillary_rise` · `laplaceJump(sigma, R1, R2)` ↔ `ch01.laplace_pressure_jump`.
- **Modes:** tank · capillary tube.
- **Views:** 1. `tank` "The tank" — two layers (top density adjustable) with a draggable 10 cm block; orange pressure
  arrows on all four faces (length ∝ gauge pressure), rose weight arrow, teal net arrow; capillary mode: tube, meniscus of
  radius R/sin α, rise h. 2. `pz` "p(z)" — absolute (solid) and gauge (dashed) pressure vs depth, probe line at the
  block's faces. 3. `hR` "h vs R" (`hidePortrait`, capillary mode) — h = 2σ sin α/(ρgR) log–log with the current tube.
- **Controls:** `rhoTop` 500…1000 kg/m³ · `hTop` 0…5 m · `depth` 0.1…9 m (block centre) · `rhoBlock` 100…8000 kg/m³ ·
  capillary `R` 0.05…2 mm, `alphaDeg` 0…90° (book's α), *optional* `sigma`.
- **Depth features:** Explain + Code · **presets** (fresh water · oil over water · mercury manometer · glass capillary ·
  waxy tube α = 0) · **terms** ("Vertical forces", N: top face orange, bottom face orange, weight rose, net teal) ·
  **inspector** (click a face: "p = p_S + ρ₁gh₁ + ρ₂g(z_s − h₁) = … Pa · F = pA = … N") · **status** ("🟢 floats: buoyancy
  9.81 N > weight 5.89 N" / "🔴 sinks" / "capillary rise h = 14.9 mm") · **modes** · **linked views**.
- **Explain:** 0 views/colours · 1 pressure at the top and bottom faces (formula = numbers = result) · 2 the difference ×
  area = ρgV (D37) · 3 weight and net force, verdict · 4 gauge vs absolute · 5 capillary: meniscus radius R/sin α, Laplace
  jump 2σ sin α/R, p just below the meniscus, rise h (hint: h–R view hidden on phones) · 6 Reading the current setting
  (floats/sinks, why the body's density does not enter buoyancy, why thinner tubes lift higher).
- **Derivation tab:** D05 (8 steps, `view: 'tank'`, step 6 live with the current density), D37 (7 steps, `view: 'tank'`,
  step 5 live "p(z₁) − p(z₂) = ρg(z₂ − z₁) = …", watch "the teal net arrow's length").
- **Code:** `ch01.layered_pressure(z, [h1, 20], [rho1, 1000])` · `ch01.net_pressure_force_on_box(...)` ·
  `ch01.buoyancy_force(rho, V)` · `ch01.capillary_rise(sigma, alpha, rho, R)` with live comments (≤ 12 lines).
- **Walkthrough (5 steps):** 1 "Which way does pressure push?" · 2 "Sideways forces cancel" (Pascal) · 3 "Bottom minus top =
  weight of displaced water" (`derive: {id: 'D37', step: 5}`) · 4 "Oil on top" (kink in p(z)) · 5 "Your turn: the thin tube"
  (predict h for R = 0.1 mm).
- **Equations:** (1.5), (1.7), (1.8), (1.9), capillary rise (Ex. 1.1), buoyancy (D37).
- **Check yourself:** (1) "Why doesn't the block's density change the buoyancy?" — "Buoyancy is the face-pressure
  difference, set by the fluid only." · (2) "Oil (800) over water: buoyancy on a block half in each?" — "g(800·V/2 +
  1000·V/2)." · (3) "Halve R in water: h?" — "doubles (h ∝ 1/R)."
- **Selftest parity rows:** `{name: 'layered p', py: 'ch01.layered_pressure(-3.0, [1.0, 2.0], [800.0, 1000.0])', rtol:
  1e-12}` · `{name: 'capillary rise', py: 'ch01.capillary_rise(0.0728, 1.5707963267948966, 998.2, 0.001)', rtol: 1e-12}` ·
  `{name: 'buoyancy', py: 'ch01.buoyancy_force(1000.0, 0.001)', rtol: 1e-12}`.
- **Fit plan:** portrait `tank` + `pz`; `hR` hidden; terms in walkthrough steps 2–3.

## Part D — runtime budget (full run < 5 min on a laptop / Colab CPU)

Estimates from timing the WIP functions (2026-09-12): FTCS N = 201 over 88 889 steps took 0.63 s; `sample_density` for
25 sizes × 200 boxes 2 ms; both ★★★ sympy checks < 2 s together. Matplotlib animations cost ≈ 0.12 s per frame at dpi 80.

| Section | Heaviest cells | Full | FAST (`FLUIDPY_FAST=1`) |
|---|---|---|---|
| setup + imports | kernel, numpy/scipy/sympy/plotly/pint imports | 8 s | 8 s |
| §1.1–1.3 | book map, γ(t) sketch | 2 s | 2 s |
| §1.4 C06 | from-scratch counting (300 boxes × 10⁵ molecules), zoom animation 37 frames (frames player), Kn slider 17 steps, 2×10⁵ Maxwellian velocities | 12 s | 5 s (120 boxes × 4×10⁴, 19 frames, 5×10⁴ velocities, 60 samples per size) |
| §1.5 C12 | FTCS N = 101 (≈ 33 000 steps, 61 saved profiles), video animation 60 frames | 9 s | 5 s (N = 51, 30 frames) |
| §1.6 | σ and Laplace numbers | < 1 s | < 1 s |
| §1.7 C20 | `integrate_hydrostatic` (DOP853, rtol 10⁻¹⁰), 10 000-step Euler loop, p(z) figure, slider 25 steps | 3 s | 3 s |
| §1.8 C25/C35/C36 | piston animation 80 frames (frames player), sympy Gibbs residual, p–v and T–s figures | 14 s | 8 s (40 frames) |
| §1.9 C40/C45 | 3-D surface (40 × 40), isothermal slider 21 steps × 400 points, 20 000-step Euler loop | 4 s | 4 s |
| §1.10 C50–C55 | parcel video 60 frames, `parcel_ode_from_gradients` (solve_ivp rtol 10⁻¹⁰), D19 sympy, two sliders (26 steps each; `synthetic_boundary_layer_column` 26 × 400 points), live cell (not executed interactively) | 16 s | 10 s (30 frames) |
| §1.11 C64–C69 | 35 minors, D28 sympy, 60-state invariance scatter, collapse figure (≈ 250 synthetic points) | 5 s | 5 s |
| explainer cells (5 × `show_viz`) | read the HTML files | 1 s | 1 s |
| **Total** | | **≈ 75 s** | **≈ 52 s** |

**FAST plan.** Every size-dependent choice is written `a if not FAST else b` in the cell (listed above); animation frame
counts ≤ 80 (≤ 40 in FAST); plotly sliders ≤ 26 steps × ≤ 4 traces × ≤ 400 points (each figure < 250 kB); no arrays are
cached to disk — the heaviest arrays (C06 sweep, C12 FTCS history, C25 route samples) are computed once in their block's
first code cell and reused by the figure, animation and check cells of the same block. Outputs: 3 video animations
(< 2 MB each), 3 frame players (< 4 MB each), page well under 15 MB.

## Part E — prerequisite ledger
Every concept, symbol, maths tool and Python function or idiom the notebook or its explainers use, with where it is
explained. "primer (in Cxx)" = a 📎 primer placed in that block before first use; the primer term is the Concept text
before any parenthesis. IDs are curation IDs (A items explain themselves; B/C items are explained by their note).

| Concept | First used in | Explained by |
|---|---|---|
| matplotlib figures (plt.subplots, ax.plot, labels) | §1.1 (C01 book map) | primer (in §1.1, before C01's figure) |
| importing a drawing helper from scripts/ | §1.1 | C01 (code comment in the book-map cell) |
| fluid mechanics and its three routes | §1.1 | C01 |
| SI base and derived units (N, Pa, J, W, Hz) | §1.2 | R01 |
| SI prefixes | §1.2 | R02 |
| Celsius–kelvin relation, absolute temperature | §1.2 | R03 |
| pint quantities (Q_, .to, .dimensionality, offset units) | §1.2 | primer (in §1.2, R01) |
| numpy arrays (vectorised arithmetic) | §1.2 | primer (in §1.2, R03) |
| f-strings | §1.2 | primer (in §1.2) |
| stress (normal and shear, strain angle γ, strain rate) | §1.3 (C02) | primer (in §1.3, before C02) |
| fluid versus solid | §1.3 | C02 |
| elastic shear modulus G | §1.3 (C02) | C02 (gloss in the note) |
| np.linspace and np.logspace | §1.3 (N02 sketch) | primer (in §1.3) |
| plastic and viscoelastic materials | §1.3 | N01 |
| normal stress, compression, tension, cavitation | §1.3 | C03 |
| vapour pressure | §1.3 (C03) | C03 (gloss in the note) |
| dot product and unit normal vector | §1.3 (C03) | C03 (gloss in the note) |
| mole, kilomole, molecular weight and Avogadro's number (A_o, M_w, m = M_w/A_o) | §1.3 (C04) | primer (in §1.3, before C04) |
| molecular spacing n^{-1/3}; liquid vs gas; free surface | §1.3 | C04 |
| continuum hypothesis, averaging volume δV, ρ = δm/δV | C06 | C06 |
| number density n and molecular mass m | C06 | C04 (uses P07 of §1.3) |
| temperature as molecular kinetic energy (k_B) | C06 (D34) | primer (in C06) |
| Newton's second law and momentum (free-body diagram) | C06 (D34) | primer (in C06) |
| np.random.default_rng (normal, poisson) | C06 | primer (in C06) |
| Gaussian velocity components and mean square speed | C06 (D34) | primer (in C06) |
| pressure as molecular impacts | C06 | C05 |
| kinetic pressure p = ⅓nm⟨u²⟩ | C06 | C06 (D34) |
| Poisson counting statistics (relative noise N^{-1/2}) | C06 (D35) | primer (in C06) |
| power laws and log–log plots (ax.loglog, np.polyfit slope) | C06 | primer (in C06) |
| relative density noise (nL³)^{-1/2} | C06 | C06 (D35) |
| tuple unpacking | C06 | primer (in C06) |
| assert np.allclose | C06 | primer (in C06) |
| animate and show_animation | C06 | primer (in C06) |
| slider_figure | C06 | primer (in C06) |
| show_viz | C06 | primer (in C06) |
| sinc(u) = sin u/u (box average of a sinusoid) | C06 (E1 Explain §3) | C06 (stated in the figure's code comment and E1 Explain §3) |
| Knudsen number Kn = l/L; continuum / slip / transition / free-molecular | C06 | C07 |
| mean free path l (Jennings, Sutherland) | C06 | N03 |
| U.S. Standard Atmosphere 1976 (standard_atmosphere) | C06 (Kn slider) | C06 (code comment; cited in reference/ch01/SOURCES.md) |
| fluid particle | C06 | C24 |
| relaxation time, equilibrium, thermodynamic system | C06 (C24) | C23 (meaning given inline in C24, defined in C25 block) |
| collision time l/c̄, mean molecular speed | C06 (C24) | C24 |
| ordinary derivative as a slope | C12 | primer (in C12) |
| diffusion down a gradient; flux | C12 | C08 |
| model diffusion equation ∂f/∂t = D∂²f/∂y² | C12 | C08 |
| mass fraction Y, partial density ρY | C12 | C09 |
| gradient vector ∇ | C12 (C10) | C10 (gloss in the note) |
| Fick's law J_m = −ρκ_m∇Y, mass diffusivity κ_m | C12 | C10 |
| Fourier's law q = −k∇T, conductivity k | C12 | C11 |
| thermal diffusivity κ = k/ρC_p (informal C_p) | C12 (C11) | C11 (gloss in the note) |
| boundary conditions (no-slip, steady vs transient) | C12 | primer (in C12) |
| Newton's law of viscosity τ = μ du/dy, dynamic viscosity μ | C12 | C12 |
| shear stress as momentum flux | C12 | C12 |
| viscosity temperature dependence (gas ↑, liquid ↓; Sutherland) | C12 | C13 |
| kinematic viscosity ν = μ/ρ, diffusion time h²/ν | C12 | C14 |
| Couette flow (steady linear profile) | C12 | C12 |
| Couette start-up series solution (check only, Ch. 8) | C12 | C12 (code comment "derived in Ch. 8; used as a check") |
| Python dictionaries | C12 | primer (in C12) |
| finite differences (central difference, FTCS, stability limit r ≤ ½) | C12 | primer (in C12) |
| np.gradient (edge_order=2) | C12 | primer (in C12) |
| matplotlib subplots with two panels, fixed axis limits | C12 | C12 (code comment) |
| relaxation of a sheared profile | C12 | N05 |
| linearity and first-derivative form of transport laws | C12 | N06 |
| Prandtl number ν/κ | E2 (Explain §6) | C12 (E2 Explain §6 defines it with numbers) |
| surface tension σ (force per length, energy per area) | §1.6 | C15 |
| spherical pressure jump 2σ/R | §1.6 | N07 |
| radius of curvature, principal radii, signed radii | §1.6 (C16) | C16 (gloss in the note) |
| Laplace pressure jump (1.5) | §1.6 | C16 |
| capillarity, capillary tube | §1.6 | N08 |
| drawing of principal radii | §1.6 | N09 |
| hydrostatic law dp/dz = −ρg | C20 | C20 |
| weight and gravitational acceleration (W = mg, g = 9.80665) | C20 | primer (in C20) |
| partial derivative ∂/∂x | C20 | primer (in C20) |
| first-order Taylor expansion (orders of smallness) | C20 (D05) | primer (in C20) |
| absolute and gauge pressure, bar, standard atmosphere | C20 | C17 |
| pressure is isotropic (1.6) | C20 | C18 |
| stress tensor −pδ_ij (pointer) | C20 | N10 |
| Pascal's law (1.7) | C20 | C19 |
| definite integral | C20 (C21, D37) | primer (in C20) |
| uniform-density hydrostatics p = p₀ − ρgz (1.9) | C20 | C21 |
| pressure difference balances weight (cube) | C20 | N11 |
| net force from pressure (surface integral −∮p n dA) | C20 (D37) | primer (in C20) |
| buoyancy F = ρ_fluid g V | C20 | C20 (D37) |
| capillary rise h = 2σ sin α/(ρgR); book's α vs contact angle | C20 | C22 |
| pressure along the tube axis (meniscus) | C20 | N12 |
| functions as arguments and lambda | C20 | primer (in C20) |
| explicit stepping (Euler, Euler–Cromer) | C20 (from scratch) | primer (in C20) |
| scipy.integrate.solve_ivp | C20 (inside integrate_hydrostatic) | primer (in C20) |
| layered and stratified still water (layered_pressure) | C20 | C20 |
| plotly slider step title update (method "update") | C54 | C54 (code comment in the slider cell) |
| thermodynamic system, equilibrium, relaxation time | C25 | C23 |
| internal and kinetic energy per unit mass | C25 | primer (in C25) |
| perfect-gas law as a working model | C25 | primer (in C25) |
| differentials: exact d and inexact δ | C25 | primer (in C25) |
| first law δq + δw = Δe (sign of work) | C25 | C25 |
| specific volume v = 1/ρ | C25 | N13 |
| path functions and state functions; reversible process | C25 | C26 |
| reversible first law de = dq − p dv | C25 | C27 |
| equations of state; two properties fix the state | C25 | C28 |
| line integral along a path (∫p dv as an area) | C25 | primer (in C25) |
| natural logarithm and exponential (e-folding) | C25 | primer (in C25) |
| trapezoid rule (np.trapezoid, cumulative_trapezoid) | C25 | primer (in C25) |
| ax.fill_between (shaded work area) | C25 | C25 (code comment) |
| p–v diagram | C25 | C25 |
| entropy s; second law (i) | C35 | C33 |
| T ds = dq (1.17) | C35 | N15 |
| partial derivative with a variable held fixed (( )_p notation) | C35 | primer (in C35) |
| enthalpy h = e + pv | C35 | C29 |
| specific heats C_p, C_v | C35 | C30 |
| C_v definition (1.15) | C35 | C31 |
| heat per degree only along reversible paths | C35 | C32 |
| product rule for differentials | C35 (D10) | primer (in C35) |
| Gibbs relations T ds = de + p dv = dh − v dp | C35 | C35 (D10) |
| Clausius–Duhem inequality (actual heat) | C35 | C34 |
| positivity of μ and k | C35 | N14 |
| sympy (symbols, diff, simplify) | C35 | primer (in C35) |
| T–s diagram | C35 | C35 |
| stirring and free expansion (irreversible processes) | C35 | C34 |
| speed of sound c² = (∂p/∂ρ)_s | C36 | C36 |
| bulk modulus K = ρ(∂p/∂ρ)_s; incompressible limit c → ∞ | C36 | C36 |
| modified Tait equation of state (tait_pressure) | C36 | C36 (code comment, cited model) |
| thermal expansion coefficient α (1.20); water's anomaly | C36 | C37 |
| molecular perfect-gas law pV = nk_BT | C40 | C38 |
| k_B, A_o, R_u, M_w, R | C40 | C39 |
| perfect-gas law p = ρRT | C40 | C40 (D11) |
| plotly 3-D surface (go.Surface) | C40 | primer (in C40) |
| e = e(T), h = h(T); van der Waals contrast | C40 | C41 |
| α = 1/T for a perfect gas | C40 | C48 |
| isothermal atmosphere p = p₀e^{−gz/RT} | C40 | C62 |
| scale height H = RT/g | C40 | C63 |
| R = C_p − C_v | C45 | C42 |
| total versus partial derivative | C45 (C42) | C42 (gloss in the note) |
| ratio of specific heats γ; γ by atomicity | C45 | C43 |
| air values of γ and C_p | C45 | N16 |
| adiabatic versus isentropic | C45 | C44 |
| separation of variables | C45 (D14) | primer (in C45) |
| exponent rules | C45 (D14) | primer (in C45) |
| isentropic law p/ρ^γ = const | C45 | C45 (D14) |
| isentropic ratios (1.26) | C45 | C46 |
| c = √(γRT); Newton's isothermal error | C45 | C47 |
| displaced parcel, rest height z_o, displacement ζ | C50 | C50 |
| isentropic (parcel) density gradient dρ_a/dz | C50 (D18) | C50 (D18 step 4) |
| linear second-order ODE (cos, cosh solutions) | C50 | primer (in C50) |
| parcel equation ζ″ + N²ζ = 0 | C50 | C50 (D18) |
| Brunt–Väisälä frequency N² (1.29) | C50 (D18 step 10) | C51 |
| square root of a negative number (growth rate) | C51 | primer (in C51) |
| period 2π/N and e-folding time 1/√(−N²) | C51 | C51 |
| stable, neutral, unstable | C51 | C52 |
| np.where and np.select | C51 | primer (in C51) |
| live widgets (ipywidgets) | C51 | primer (in C51) |
| seawater salinity S, ρ(T, p, S), linear EOS | C51 | C60 |
| ocean parcel density gradient −ρg/c² | C51 | N21 |
| ocean stability criterion (1.35), "same sign as" | C51 | C61 |
| static medium: one profile fixes p, ρ, T | C54 | C49 |
| lapse rate Γ, Kundu Γ ≡ dT/dz and meteorology Γ ≡ −dT/dz | C54 | C53 |
| inequalities under a sign change | C54 | primer (in C54) |
| lapse_rate_convention and lapse_rate_stability (NamedTuple fields .verdict, .text, .code) | C54 | C54 (code comments) |
| chain rule | C54 (D19) | primer (in C54) |
| Gibbs free energy g = h − Ts | C54 (D19) | primer (in C54) |
| exact differentials and Maxwell relations (equality of mixed partials) | C54 (D19) | primer (in C54) |
| adiabatic lapse rate Γ_a = −gαT/C_p | C54 | C54 (D19) |
| value of Γ_a in both conventions | C54 | N17 |
| inversion, superadiabatic layer | C54 | C54 (figure notes) |
| lapse-rate feedback (climate literature convention) | C54 | C54 (problem text) |
| dry adiabat through a point (parcel_temperature) | C54 | C54 |
| potential temperature θ, reference pressure p_o (≠ p₀) | C55 | C55 (D20) |
| (T/θ)dθ/dz = Γ − Γ_a (1.32) in both conventions | C55 | C56 |
| log-derivative of (1.31) | C55 | N19 |
| logarithmic differentiation | C55 (D36) | primer (in C55) |
| N² = (g/θ)dθ/dz | C55 | C55 (D36) |
| stability from the sign of dθ/dz; lab-scale T ≈ θ | C55 | C57 |
| boundary-layer profile: mixed layer, inversion, stable layer | C55 | N18 |
| potential density ρ_θ (1.33) | C55 | C58 |
| θρ_θ = p_o/R | C55 | N20 |
| −(1/ρ_θ)dρ_θ/dz = (1/θ)dθ/dz (1.34) | C55 | C59 |
| moist adiabat (pointer only) | C55 | C55 ("What would change if…" text) |
| units versus dimensions; exponent vector [q] | C64 | C64 |
| dimensional homogeneity; dimensionless form of a law | C64 | C64 |
| pint.DimensionalityError | C64 | C64 (code comment) |
| rescale_units and group_value | C64 | C64 (code comments) |
| choosing the variable list (1.38) | C67 | C65 |
| [q] notation, base dimensions M, L, T, Θ; kmol dropped | C67 | C66 |
| matrices, determinants and minors (cofactor expansion) | C67 | primer (in C67) |
| linear independence and rank | C67 | primer (in C67) |
| dimensional matrix (1.39) | C67 | C67 |
| rank by minors | C67 | C68 |
| itertools.combinations | C67 | primer (in C67) |
| np.linalg.det and np.linalg.matrix_rank | C67 | primer (in C67) |
| ax.imshow heatmap with cell annotations | C67 | C67 (code comment) |
| relation among n variables (1.36) | C69 | N22 |
| number of groups n − r | C69 | N23 |
| null space and rank–nullity theorem | C69 | primer (in C69) |
| scaling the base units (unit-change invariance) | C69 (D28) | primer (in C69) |
| sympy Matrix and nullspace | C69 | primer (in C69) |
| Buckingham Π theorem (1.37) | C69 | C69 (D28) |
| repeating variables; exponent algebra | C69 | C70 |
| groups by inspection | C69 | N24 |
| combining groups; independence | C69 | C71 |
| dimensionless pipe law (1.40); Reynolds number | C69 | C72 |
| using physics to reduce groups (Δp ∝ Δx; 32/Re) | C69 | N25 |
| fractions.Fraction | C69 | primer (in C69) |
| np.linalg.solve | C69 | primer (in C69) |
| Hagen–Poiseuille pressure drop (data generator, Ch. 8) | C69 | C69 (code comment "derived in Ch. 8") |
| scale height by dimensional analysis (Ex. 1.2) | C69 | C73 |
| Pythagoras by dimensional analysis; similar triangles | C69 | C74 |
| blast-wave energy E = KρD⁵/t²; similarity solution | C69 | C75 |
| Rayleigh scattering λ⁻⁴; light as waves, intensity | C69 | C76 |
| exercises and reading list | end of §1.11 | pointer (S01, S02 lines) |
| symbol ρ density [kg/m³] | C06 | C06 |
| symbol p pressure [Pa] | §1.2 (R01) | R01 |
| symbol T temperature [K] | R03 | R03 |
| symbol μ dynamic viscosity [Pa s] | C12 | C12 |
| symbol ν kinematic viscosity [m²/s] | C12 | C14 |
| symbol σ surface tension [N/m] | §1.6 | C15 |
| symbol g gravitational acceleration [m/s²] | C20 | C20 (weight and g, P24) |
| symbol z height, upward; depth h = −z | C20 | C21 |
| symbols e, h, s, v per unit mass | C25 | C25 |
| symbols C_p, C_v, γ, R | C35 | C30 |
| symbol α thermal expansion [1/K] (and α contact angle in C22) | C36 | C37 |
| symbol c speed of sound [m/s] | C36 | C36 |
| symbol ζ parcel displacement [m] | C50 | C50 |
| symbol N² [1/s²] | C50 | C51 |
| symbol Γ, Γ_a lapse rates [K/m] in both conventions | C54 | C53 |
| symbol θ potential temperature [K] (≠ Θ dimension, ≠ wedge angle) | C55 | C55 |
| symbol Π dimensionless group; n, r | C69 | C69 |

### Summary of counts (Parts A–F)
Part A: 15 CORE blocks (C06, C12, C20, C25, C35, C36, C40, C45, C50, C51, C54, C55, C64, C67, C69); 86 NOTE IDs (68 at
depth B, 18 at depth C); 3 RECAP (R01, R03 at B; R02 at C); 2 SKIP pointers (S01, S02) — all 106 curation IDs in the
placement table (A.12), each exactly once. 61 primers + 10 glosses. 12 derivations (Part F). 5 explainers embedded (E1 in
C06, E2 in C12, E3 in C45, E4 in C55, E5 in C69) + backup B1 storyboarded. Part C: 144 callables. Part E: 225 ledger
rows.

## Part F — derivation storyboards

Builders copy these word for word into `nb.derivation(key, title, goal=…, start=(tex, plain), plan=[…], uses=[…],
steps=[dict(did, tex, why, plain)], result=(tex, plain), interpret=…, check=…, check_src=…)` and into the explainer's
`derivations: [...]` (phones may shorten *why* to its first sentence). Every step is one move; *why* names the rule and
says why we make the move. The book's own moves were read on the rendered pages (p038 for D05, p043 for D10, p044 for
D11/D14, p045–p046 for D18/D19, p047 for D20/D36, p049–p052 for D28); the rest fills the gaps listed in
`analysis/ch01.md` §2b. Colours: pressure orange, weight/buoyancy teal, work blue, heat orange, state functions purple.
(No line in this part starts with a table bar.)

### D05 · The hydrostatic law, Eq. (1.8) — ★, 8 steps, in C20 (notebook; B1 if built)
- **Goal.** Find how pressure must change with height inside a fluid that is not moving, so that every little piece of
  it stays at rest. This is the base state of everything stratified in the book.
- **Start.** $\sum F_z = m\,a_z = 0$ — *in words:* a small cube of fluid at rest has zero acceleration, so the vertical
  forces on it add to zero (Newton's second law, P09).
- **Plan.** (1) List the vertical forces on a cube dx × dy × dz: pressure on the bottom, pressure on the top, weight.
  (2) Relate the top pressure to the bottom one with a Taylor step. (3) Add the forces, cancel, divide by the volume.
- **Tools.** Newton's second law (P09, C06) · weight W = mg (P24, C20) · pressure acts normal to a surface and is the same
  in every direction (C18) · Pascal's law, p independent of x and y at rest (C19) · partial derivative (P25) · first-order
  Taylor expansion (P26).
- **Assumptions.** Fluid at rest (step 1: a = 0; step 4: p = p(z) only) · z points up (steps 2–3: signs of the face
  forces) · g uniform over the tiny cube (step 5) · dz small enough for a first-order Taylor step (step 3).
- **Steps.**
  1. *did:* Write the push on the bottom face · *tex:* $F_{\rm bottom} = +\,p(z)\,dx\,dy$ · *why:* Pressure pushes into the
     fluid normal to a face; on the bottom face that is upward (positive z). We start the force list with it because the
     cube is held up from below. · *plain:* The fluid underneath pushes the cube up with pressure times area.
  2. *did:* Write the push on the top face · *tex:* $F_{\rm top} = -\,p(z+dz)\,dx\,dy$ · *why:* On the top face the fluid
     above pushes down, so the force is negative with z up; it uses the pressure at the top's height z + dz. · *plain:* The
     fluid above pushes the cube down with the pressure found one step higher.
  3. *did:* Expand the top pressure to first order · *tex:* $p(z+dz) = p(z) + \dfrac{\partial p}{\partial z}\,dz$ · *why:*
     First-order Taylor expansion (P26); the dropped terms are proportional to dz² and will vanish when we divide by the
     volume and let the cube shrink. We need the two face pressures in terms of one value. · *plain:* The top pressure is
     the bottom pressure plus slope times height.
  4. *did:* Use Pascal's law to make it an ordinary derivative · *tex:* $\dfrac{\partial p}{\partial z} = \dfrac{dp}{dz}$ ·
     *why:* At rest p does not change with x or y (C19, Eq. 1.7), so p depends on z alone and its partial derivative in z
     is the ordinary one. The book writes dp directly; this is why it may. · *plain:* In still fluid pressure depends only
     on height.
  5. *did:* Write the weight · *tex:* $W = -\,\rho g\,dx\,dy\,dz$ · *why:* Mass = density × volume, weight = mass × g, pointing
     down (P24); ρ and g are treated as constant across the tiny cube. · *plain:* The cube's weight pulls it down.
  6. *did:* Add the three forces and set the sum to zero · *tex:* $p\,dx\,dy - \Big(p + \dfrac{dp}{dz}dz\Big)dx\,dy - \rho
     g\,dx\,dy\,dz = 0$ · *why:* Newton's second law with zero acceleration (the fluid is at rest), using steps 1–5. The
     horizontal forces already cancel by Pascal's law, so only this balance is left. · *plain:* Push up, push down and
     weight exactly cancel.
  7. *did:* Cancel p dx dy and divide by dx dy dz · *tex:* $-\dfrac{dp}{dz} - \rho g = 0$ · *why:* The terms p dx dy appear
     with opposite signs; the remaining terms all contain the volume dx dy dz ≠ 0, so we may divide by it. This removes the
     arbitrary size of the cube. · *plain:* Per unit volume, the upward pressure-gradient force balances the weight.
  8. *did:* Move the weight term across · *tex:* $\dfrac{dp}{dz} = -\rho g$ · *why:* Add ρg to both sides, then multiply both sides by −1 (an equation, not an inequality, so nothing else changes). We isolate dp/dz because it is what we want to know. ·
     *plain:* Pressure falls with height by ρg per metre.
- **Result.** $\dfrac{dp}{dz} = -\rho g$ (1.8) — *in words:* going up one metre in still fluid, pressure drops by the weight
  of a one-metre column of unit cross-section.
- **Check.** Units: Pa/m = kg m⁻² s⁻² and ρg = kg m⁻³ · m s⁻² = kg m⁻² s⁻² ✓. Limit: ρ → 0 (vacuum) gives constant p ✓.
  Number: water, ρg = 1000 × 9.81 = 9810 Pa/m, so 10 m of depth adds 9.81×10⁴ Pa ≈ 1 atm ✓ (notebook code prints
  `hydrostatic_pressure_uniform(-10, P_ATM, 1000)`).
- **What it means.** Pressure in still fluid is the weight of everything above per unit area. It holds for liquids and
  gases alike (ρ may vary with z or p — then integrate numerically, as `integrate_hydrostatic` does). It fails when the
  fluid accelerates (Ch. 4 adds ρ Du/Dt) or when g varies over the column (planetary scales).
- **Traps.** With z pointing *down* (depth) the sign flips: dp/dh = +ρg. Forgetting that the top-face force points down
  gives dp/dz = +ρg. Writing dp/dz before knowing p depends on z only.

### D37 · Buoyancy: the net pressure force on a submerged body — ★, 7 steps, in C20 (notebook; B1 if built)
- **Goal.** Show that still fluid pushes up on a submerged body with a force equal to the weight of the fluid the body
  displaces — whatever the body is made of. The parcel argument (D18) needs this.
- **Start.** $p(z) = p_0 - \rho g z$ (1.9) — *in words:* in fluid of uniform density ρ the pressure grows linearly with
  depth.
- **Plan.** (1) Put a box of horizontal area A between heights z₁ (bottom) and z₂ (top) into the fluid. (2) Show the side
  forces cancel. (3) Subtract top from bottom force and use (1.9). (4) Build any shape from thin boxes.
- **Tools.** (1.9) (C21) · Pascal's law (C19) · net force from pressure as a surface sum (P28) · definite integral (P27).
- **Assumptions.** Fluid at rest with uniform density near the body (step 5) · the body is fully submerged and does not
  disturb the fluid's pressure (steps 2–6) · g uniform.
- **Steps.**
  1. *did:* Pair up the side faces · *tex:* $F_{x,\rm left} + F_{x,\rm right} = p\,A_x - p\,A_x = 0$ · *why:* Opposite side
     faces sit at the same heights, so by Pascal's law (C19) they feel the same pressure at each height, pushing in opposite
     directions. We clear the horizontal forces first to be left with a one-dimensional balance. · *plain:* Sideways pushes
     cancel.
  2. *did:* Write the upward force on the bottom face · *tex:* $F_{\rm bottom} = p(z_1)\,A$ · *why:* Pressure acts normal to
     the face, upward on a bottom face; area A. · *plain:* The fluid below pushes up.
  3. *did:* Write the downward force on the top face · *tex:* $F_{\rm top} = p(z_2)\,A$ · *why:* Same rule on the top face,
     pointing down; z₂ > z₁. · *plain:* The fluid above pushes down, less hard because it is higher.
  4. *did:* Take the net upward force · *tex:* $F_{\rm net} = \big[p(z_1) - p(z_2)\big]\,A$ · *why:* Up minus down, the
     vertical part of the surface sum of P28. · *plain:* What is left is the difference of the two face pressures.
  5. *did:* Substitute the hydrostatic pressure · *tex:* $p(z_1) - p(z_2) = \rho g\,(z_2 - z_1)$ · *why:* From (1.9):
     (p₀ − ρgz₁) − (p₀ − ρgz₂); p₀ cancels. This is where the fluid's density enters — and only the fluid's. · *plain:* The
     pressure difference is ρg times the box's height.
  6. *did:* Recognise the volume · *tex:* $F_{\rm net} = \rho g\,A\,(z_2 - z_1) = \rho\,g\,V$ · *why:* A (z₂ − z₁) is the box
     volume V. · *plain:* The net push up equals the weight of fluid that would fill the box.
  7. *did:* Extend to any shape with thin vertical columns · *tex:* $F_{\rm net} = \displaystyle\int_A \rho g\,h(x,y)\,dA =
     \rho g V$ · *why:* Slice the body into vertical columns of cross-section dA and height h(x, y); each is a thin box (its
     side forces cancel against its neighbours', and a slanted cap's vertical pressure force is p times its horizontal
     projection). Adding the columns is the definite integral of P27; ∫h dA is the body's volume. · *plain:* Any submerged
     body feels an upward force equal to the weight of the displaced fluid.
- **Result.** $F_b = \rho_{\rm fluid}\,g\,V$ — *in words:* Archimedes' principle, derived from the hydrostatic pressure.
- **Check.** Units: kg m⁻³ · m s⁻² · m³ = N ✓. Special cases: ρ_fluid → 0 gives no buoyancy ✓; a body of the fluid's own
  density feels weight = buoyancy, it is neutral ✓. Number: 1 L in water → 9.81 N; `net_pressure_force_on_box` integrates
  the six faces and returns (0, 0, 9.807) N ✓.
- **What it means.** Buoyancy is nothing but the pressure difference between a body's bottom and top. It does not care what
  the body is made of; the body's own weight decides whether it rises or sinks. In a stratified fluid (D18) the relevant
  density is the environment's at the body's current height. It fails if the body touches the bottom (no pressure
  underneath) or is only partly submerged (use the submerged volume).
- **Traps.** Thinking side forces add up; putting the body's density into the buoyancy; forgetting that the net force
  points up.

### D34 · Pressure from molecular impacts — ★★, 8 steps, in C06 (notebook)
- **Goal.** Turn "pressure is the average push of molecules hitting a wall" into a formula, and see that it gives the gas
  law pV = nk_BT used in §1.9.
- **Start.** $p = \dfrac{\text{momentum delivered to the wall}}{\text{area} \times \text{time}}$ — *in words:* force is
  momentum per time (Newton), pressure is force per area.
- **Plan.** (1) Momentum given to the wall by one bounce. (2) How many molecules of one velocity hit in a time dt. (3) Add
  over all velocities. (4) Use isotropy and the kinetic meaning of temperature.
- **Tools.** Newton's second law and momentum (P09) · Gaussian velocity components and mean square speed (P11) ·
  temperature as molecular kinetic energy (P08).
- **Assumptions.** Ideal gas: molecules do not interact except by brief collisions (step 2: straight flight to the wall)
  · elastic, smooth wall (step 1) · isotropic velocities, gas at rest (steps 5–6) · equilibrium, so the kinetic definition
  of T applies (step 8).
- **Steps.**
  1. *did:* Momentum given by one bounce · *tex:* $\Delta P_{\rm wall} = 2\,m\,u_x$ · *why:* An elastic bounce off a smooth
     wall normal to x reverses u_x and keeps u_y, u_z; the molecule's momentum changes by −2mu_x, so by Newton's third law
     the wall gains +2mu_x. This is the "push" we want to count. · *plain:* Each hit hands the wall twice the molecule's
     normal momentum.
  2. *did:* Count hits of one velocity class in time dt · *tex:* $dN = n_{u_x}\,A\,u_x\,dt \quad (u_x > 0)$ · *why:* Only
     molecules moving toward the wall (u_x > 0) and within a distance u_x dt of it arrive in dt: they fill a slab of volume
     A u_x dt, and n_{u_x} is the number per volume with that velocity. · *plain:* Faster molecules reach the wall from
     farther away.
  3. *did:* Multiply hits by momentum per hit · *tex:* $dP = 2\,m\,n_{u_x}\,u_x^2\,A\,dt$ · *why:* Total momentum = number
     of hits × momentum per hit (steps 1 and 2); u_x appears twice, once from how many arrive and once from how hard each
     hits. · *plain:* The delivered momentum grows with the square of the speed.
  4. *did:* Divide by area and time · *tex:* $p_{u_x} = 2\,m\,n_{u_x}\,u_x^2$ · *why:* Force = momentum per time (P09),
     pressure = force per area; dividing by A dt gives this class's contribution. · *plain:* One velocity class adds this
     much pressure.
  5. *did:* Add all classes moving toward the wall · *tex:* $p = 2m\sum_{u_x>0} n_{u_x}u_x^2 = 2m\cdot\tfrac{n}{2}\langle u_x^2
     \rangle$ · *why:* By symmetry half of the n molecules per volume move toward the wall, and those have the same mean
     square u_x as the whole gas. We sum because pressure is the total push. · *plain:* The half of the molecules heading for
     the wall carry the average square speed.
  6. *did:* Cancel the 2 and the ½ · *tex:* $p = n\,m\,\langle u_x^2\rangle$ · *why:* 2 × ½ = 1. · *plain:* Pressure is
     number density × mass × mean square normal velocity.
  7. *did:* Replace ⟨u_x²⟩ using isotropy · *tex:* $p = \tfrac13\,n\,m\,\langle|\mathbf u|^2\rangle$ · *why:* In a gas at
     rest no direction is special, so ⟨u_x²⟩ = ⟨u_y²⟩ = ⟨u_z²⟩ and their sum is ⟨|u|²⟩ (P11); hence ⟨u_x²⟩ = ⅓⟨|u|²⟩.
     We want the speed, not one component. · *plain:* Pressure is one third of density times mean square speed.
  8. *did:* Use the kinetic meaning of temperature · *tex:* $p = n\,k_B\,T$ · *why:* ½m⟨|u|²⟩ = (3/2)k_BT (P08), so
     m⟨|u|²⟩ = 3k_BT and the ⅓ cancels the 3. This links the molecular picture to a measurable temperature. · *plain:*
     Pressure = molecules per volume × k_B × temperature.
- **Result.** $p = \tfrac13\,n\,m\,\langle|\mathbf u|^2\rangle = n\,k_B T$ — *in words:* multiplying by the volume V with n =
  N/V gives pV = Nk_BT, the molecular gas law (1.21) stated in C38.
- **Check.** Units: m⁻³ · kg · m² s⁻² = kg m⁻¹ s⁻² = Pa ✓. Number: sea-level air, n = 2.55×10²⁵ m⁻³, m = 4.81×10⁻²⁶ kg,
  ⟨|u|²⟩ = 3k_BT/m = 2.48×10⁵ m²/s² → p = ⅓ × 2.55×10²⁵ × 4.81×10⁻²⁶ × 2.48×10⁵ = 1.013×10⁵ Pa ✓. Code: 200 000 sampled
  molecules give `molecular_pressure` and `wall_impact_pressure` within 2 % of nk_BT.
- **What it means.** A steady pressure is the average of an enormous number of tiny, random pushes — which is why it is well
  defined only for boxes holding many molecules (D35). Heating a gas at fixed n raises p because molecules hit harder and
  more often. It fails for dense gases and liquids, where molecules interact over their whole flight.
- **Traps.** Losing the factor 2 of the bounce; forgetting that only half the molecules move toward the wall; using ⟨u_x²⟩ =
  ⟨|u|²⟩ instead of ⅓.

### D35 · How noisy is a density measured in a small box? — ★, 5 steps, in C06 (notebook · `continuum_averaging_volume`)
- **Goal.** Find how the scatter of a box's density reading depends on the box size — the lower end of the continuum window.
- **Start.** $\bar N = n\,\delta V = n\,L^3$ — *in words:* a box of side L in a gas with n molecules per m³ holds N̄
  molecules on average.
- **Plan.** (1) The molecule count is Poisson. (2) Scale its scatter into a density scatter. (3) Divide by the mean.
  (4) Write the result in terms of L.
- **Tools.** Poisson counting statistics (P12) · power laws and log–log plots (P13) · exponent rules (P43, primed later in
  C45; here only (L³)^{−1/2} = L^{−3/2}, spelled out in the step).
- **Assumptions.** Molecules placed independently (ideal gas; step 1) · the box holds a tiny fraction of all molecules
  (step 1: Poisson limit) · the macroscopic density does not vary across the box (step 3; the flow's variation is the
  separate upper bound).
- **Steps.**
  1. *did:* Use Poisson statistics for the count · *tex:* $\operatorname{Var}(N) = \bar N$ · *why:* Each of very many
     independent molecules has a tiny chance to be inside the box; such counts are Poisson, whose variance equals the mean
     (P12). This is where the randomness enters. · *plain:* The count scatters by about √N̄. · *live (E1):* N̄ = nL³ with your
     box.
  2. *did:* Take the square root · *tex:* $\sigma_N = \sqrt{\bar N}$ · *why:* Standard deviation is the square root of the
     variance. · *plain:* A box expecting 100 molecules typically sees 90 to 110.
  3. *did:* Turn counts into densities · *tex:* $\sigma_\rho = \dfrac{m\,\sigma_N}{L^3},\quad \rho = \dfrac{m\bar N}{L^3}$ ·
     *why:* The reading is δm/δV = mN/L³ (C06); multiplying a random number by the constant m/L³ multiplies its mean and its
     standard deviation by that constant. · *plain:* The density reading inherits the count's scatter.
  4. *did:* Divide by the mean · *tex:* $\dfrac{\sigma_\rho}{\rho} = \dfrac{\sigma_N}{\bar N} = \bar N^{-1/2}$ · *why:* m/L³
     cancels; √N̄/N̄ = N̄^{−1/2}. A *relative* scatter is what decides whether a value is usable. · *plain:* The relative noise
     is one over the square root of the expected count. · *live (E1):* 0.198 at L = 10 nm; *watch:* "the grey band's
     half-width at the marker equals this number".
  5. *did:* Write it in terms of the box side · *tex:* $\dfrac{\sigma_\rho}{\rho} = (nL^3)^{-1/2} = n^{-1/2}L^{-3/2}$ · *why:*
     Substitute N̄ = nL³ and use (ab)^k = a^k b^k and (L³)^{−1/2} = L^{−3/2}. The exponent tells the slope on a log–log plot
     (P13). · *plain:* Every factor 10 in box side cuts the noise by 10^{1.5} ≈ 32. · *set (E1):* logL −6; *watch:* "one
     decade of L shrinks the band 32×".
- **Result.** $\dfrac{\sigma_\rho}{\rho} = (n\,L^3)^{-1/2}$ — *in words:* the relative noise of a box density falls as the
  box side to the power −3/2.
- **Check.** Dimensionless ✓ (n L³ is a count). Limits: L → 0 gives infinite noise (a box with no molecules), L → ∞ zero ✓.
  Number: 10 µm cube of air, N̄ = 2.55×10¹⁰ → 6.3×10⁻⁶ ✓; the notebook's measured scatter follows the line of slope −1.5 and
  the from-scratch molecule count agrees within 25 %.
- **What it means.** The continuum's lower limit is statistical: a density exists once a box holds, say, 10⁶ molecules
  (noise 0.1 %), i.e. L ≈ 0.3 µm in sea-level air. In a liquid, molecules are not independent (they touch), so the true
  scatter is smaller than Poisson; the estimate is then an upper bound.
- **Traps.** Writing L^{−1/2} (the volume, not the side, is under the root); confusing absolute scatter σ_ρ with relative
  scatter σ_ρ/ρ.

### D10 · The Gibbs relations, Eq. (1.18) — ★, 7 steps, in C35 (notebook · `heat_work_paths`)
- **Goal.** Eliminate heat — which depends on the route — from the reversible first law, leaving relations between state
  functions only. Those relations then hold for *any* process.
- **Start.** $T\,ds = dq$ (1.17) and $de = dq - p\,dv$ (1.11) — *in words:* for a reversible change the heat received is T ds,
  and the first law with boundary work −p dv.
- **Plan.** (1) Solve the first law for dq and substitute into T ds = dq. (2) Rewrite with enthalpy using the product rule.
  (3) Explain why a route-independent relation holds for every process.
- **Tools.** (1.17) (N15, C35) · (1.11) (C27, C25) · h = e + pv (C29, C35) · differentials (P34) · product rule for
  differentials (P38).
- **Assumptions.** A reversible route is used as a device only (steps 1–2) · simple compressible substance: only p dv work
  (step 1) · states 1 and 2 are equilibrium states (step 7).
- **Steps.**
  1. *did:* Solve the reversible first law for dq · *tex:* $dq = de + p\,dv$ · *why:* Add p dv to both sides of (1.11). We
     isolate dq because it is the route-dependent quantity we want to replace. · *plain:* Heat in = rise of internal energy
     + work done by the gas.
  2. *did:* Substitute into T ds = dq · *tex:* $T\,ds = de + p\,dv$ · *why:* (1.17) says the same dq equals T ds on a
     reversible route; equal things can replace each other. · *plain:* The first Gibbs relation — no heat left in it. ·
     *live (E3):* for your route Δs, Δe and ∫p dv/T leg by leg.
  3. *did:* Differentiate the definition of enthalpy · *tex:* $dh = de + d(pv)$ · *why:* h = e + pv (1.13); the differential
     of a sum is the sum of differentials. We want a second form written with dh. · *plain:* A change in h is a change in e
     plus a change in pv.
  4. *did:* Apply the product rule · *tex:* $dh = de + p\,dv + v\,dp$ · *why:* d(pv) = p dv + v dp (P38): both p and v may
     change. · *plain:* Both pressure and volume changes feed pv.
  5. *did:* Isolate de + p dv · *tex:* $de + p\,dv = dh - v\,dp$ · *why:* Subtract v dp from both sides of step 4; the left
     side is exactly the right side of step 2. · *plain:* The same combination written with enthalpy.
  6. *did:* Substitute into step 2 · *tex:* $T\,ds = dh - v\,dp$ · *why:* Replace de + p dv by its equal from step 5. ·
     *plain:* The second Gibbs relation.
  7. *did:* Drop the reversible restriction · *tex:* $T\,ds = de + p\,dv = dh - v\,dp\quad\text{(any process)}$ · *why:* T, s,
     e, h, p, v are all state functions, so the relation between their changes for neighbouring equilibrium states cannot
     know which route took the fluid there; any irreversible change can be replaced by a reversible route between the same
     states. Only the *interpretation* T ds = dq is lost. · *plain:* The relations hold for friction and stirring too — but
     then T ds is no longer the heat. · *set (E3):* proc 'stirring'; *watch:* "Gibbs still gives Δs = 23.5 J/(kg K) although
     δq = 0".
- **Result.** $T\,ds = de + p\,dv, \qquad T\,ds = dh - v\,dp$ (1.18) — *in words:* entropy changes are fixed by changes of
  energy (or enthalpy), volume and pressure, for any process.
- **Check.** Units: J/kg on every term ✓. Perfect gas, 300 → 600 K at constant p: form 1 C_v ln 2 + R ln 2 = 696.4 J/(kg K),
  form 2 C_p ln 2 − R ln 1 = 696.4 ✓. The notebook's sympy cell shows T ds − de − p dv = 0 and T ds − dh + v dp = 0
  symbolically for s = C_v ln T + R ln v.
- **What it means.** Entropy is computable from measurable state properties alone. For stirring at constant volume
  (dv = 0) the first form gives T ds = de: the stirring work shows up as entropy although no heat entered — the
  Clausius–Duhem inequality (C34). The isentropic law (D14) and the lapse rate (D19) both start here.
- **Traps.** Concluding the relations hold only for reversible processes; dropping v dp in the product rule.

### D11 · From the molecular gas law to p = ρRT, Eq. (1.22) — ★, 6 steps, in C40 (notebook)
- **Goal.** Rewrite the molecule-counting gas law as a law for the continuum density, with one constant R per gas.
- **Start.** $pV = n\,k_B T$ (1.21) — *in words:* n non-interacting molecules in a volume V at temperature T exert pressure p
  (here n counts molecules, as in the book; n/V is the number density of §1.4).
- **Plan.** (1) Divide by V. (2) Create the density ρ = nm/V. (3) Replace one molecule's mass by M_w/A_o. (4) Name the
  constants.
- **Tools.** (1.21) (C38, derived in D34) · continuum density (C06) · kilomole, molecular weight and Avogadro's number (P07)
  · constants k_B, A_o, R_u (C39).
- **Assumptions.** Perfect gas (non-interacting molecules; the start) · Kn ≪ 1 so nm/V is a continuum density (step 3) · one
  average molecular weight for a mixture (step 4).
- **Steps.**
  1. *did:* Divide by the volume · *tex:* $p = \dfrac{n}{V}\,k_B T$ · *why:* Divide both sides by V ≠ 0; pressure is a
     local quantity and should not depend on how big a container we imagine. · *plain:* Pressure = molecules per volume ×
     k_B × T.
  2. *did:* Multiply and divide by the molecular mass · *tex:* $p = \dfrac{n\,m}{V}\,\dfrac{k_B}{m}\,T$ · *why:* m/m = 1
     changes nothing; we do it so that nm/V — mass per volume — appears. · *plain:* The same law, rearranged to expose mass.
  3. *did:* Recognise the density · *tex:* $p = \rho\,\dfrac{k_B}{m}\,T$ · *why:* nm/V is the total molecular mass per volume,
     the continuum density of C06 (valid for boxes in the continuum window, Kn ≪ 1). · *plain:* Density replaces the
     molecule count.
  4. *did:* Write one molecule's mass per kilomole · *tex:* $\dfrac{k_B}{m} = \dfrac{k_B A_o}{M_w}$ · *why:* One kilomole of
     molecules (A_o of them) has mass M_w kg, so m = M_w/A_o (P07); dividing by a fraction multiplies by its inverse. ·
     *plain:* Mass of one molecule = molecular weight ÷ Avogadro's number.
  5. *did:* Name the universal gas constant · *tex:* $p = \rho\,\dfrac{R_u}{M_w}\,T,\quad R_u \equiv k_B A_o$ · *why:*
     k_BA_o = 1.380649×10⁻²³ J/K × 6.02214076×10²⁶ kmol⁻¹ = 8314.46 J kmol⁻¹ K⁻¹ is the same for every gas; naming it
     shortens the law. · *plain:* One universal constant divided by the gas's molecular weight.
  6. *did:* Name the gas constant of this gas · *tex:* $p = \rho R T,\quad R \equiv R_u/M_w$ · *why:* For a given gas R_u/M_w
     is one number (air: 8314.46/28.9644 = 287.06 J kg⁻¹ K⁻¹). · *plain:* The perfect-gas law per unit mass.
- **Result.** $p = \rho R T$ (1.22) — *in words:* pressure = density × gas constant × absolute temperature.
- **Check.** Units: kg m⁻³ · J kg⁻¹ K⁻¹ · K = J m⁻³ = Pa ✓. Number: ρ = 101 325/(287.06 × 288.15) = 1.225 kg/m³ ✓
  (`perfect_gas_density`, and the from-scratch constants chain). Mixture: dry air as one gas with M_w = 28.96 ✓.
- **What it means.** Every atmospheric calculation in the book closes with this law. It fails when molecules attract each
  other or fill a noticeable part of the volume (dense gases, near condensation — van der Waals in C41).
- **Traps.** Mixing mol and kmol (R_u = 8.314 J mol⁻¹ K⁻¹ vs 8314 J kmol⁻¹ K⁻¹, a factor 1000); confusing R_u with R; using
  °C for T.

### D14 · The isentropic law for a perfect gas, Eq. (1.25) — ★★, 10 steps, in C45 (notebook · `heat_work_paths`)
- **Goal.** Find how pressure follows density when a perfect gas is compressed or expanded with no heat and no friction.
- **Start.** $T\,ds = de + p\,dv = dh - v\,dp$ (1.18) — *in words:* the Gibbs relations, valid for any process.
- **Plan.** (1) Set ds = 0 in both forms. (2) Use e = e(T) and h = h(T) to bring in C_v and C_p. (3) Divide the two
  equations so T drops out. (4) Switch from v to ρ and integrate.
- **Tools.** Gibbs relations (C35, D10) · e = e(T), h = h(T) for a perfect gas (C41) · C_p, C_v (C30, C31) · γ = C_p/C_v
  (C43) · differentials (P34) · separation of variables (P42) · natural logarithm (P36) · exponent rules (P43).
- **Assumptions.** Isentropic: adiabatic and frictionless, ds = 0 (steps 1, 3) · perfect gas (steps 2, 3) · constant C_p and
  C_v, hence constant γ (step 9).
- **Steps.**
  1. *did:* Set ds = 0 in the first Gibbs form · *tex:* $0 = de + p\,dv$ · *why:* Isentropic means s does not change, so
     T ds = 0 (T > 0). We start with the energy form because the work p dv appears in it. · *plain:* All the work done on
     the gas goes into internal energy.
  2. *did:* Use e = e(T) · *tex:* $C_v\,dT = -p\,dv$ · *why:* For a perfect gas e depends on T only (C41), so the partial
     derivative in (1.15) is an ordinary one and de = C_v dT. · *plain:* Compression (dv < 0) warms the gas.
  3. *did:* Set ds = 0 in the second Gibbs form and use h = h(T) · *tex:* $C_p\,dT = v\,dp$ · *why:* Same reasoning with
     (1.14): h depends on T only, so dh = C_p dT; ds = 0 gives dh = v dp. We need a second equation containing dp. ·
     *plain:* Raising the pressure raises the enthalpy.
  4. *did:* Divide step 3 by step 2 · *tex:* $\dfrac{C_p}{C_v} = \dfrac{v\,dp}{-\,p\,dv}$ · *why:* Along a compression dT ≠ 0,
     so both sides of step 2 are nonzero and we may divide; dT cancels, removing temperature from the problem. · *plain:*
     The ratio of specific heats links the pressure change to the volume change.
  5. *did:* Name the ratio γ · *tex:* $\gamma = -\dfrac{v\,dp}{p\,dv}$ · *why:* γ ≡ C_p/C_v (1.24). · *plain:* γ measures how
     much faster p changes than v, in relative terms.
  6. *did:* Separate the variables · *tex:* $\dfrac{dp}{p} = -\gamma\,\dfrac{dv}{v}$ · *why:* Multiply both sides by
     −dv/v; now each side contains one variable (P42). · *plain:* A 1 % decrease in volume raises the pressure by γ %. ·
     *live (E3):* at the marker, dp/p = −1.4 × dv/v.
  7. *did:* Switch from specific volume to density · *tex:* $\dfrac{dv}{v} = -\dfrac{d\rho}{\rho}$ · *why:* v = 1/ρ gives
     dv = −dρ/ρ² (derivative of 1/ρ); dividing by v = 1/ρ gives −dρ/ρ. Watch the sign: this is the classic slip. · *plain:*
     A 1 % increase in density is a 1 % decrease in volume.
  8. *did:* Substitute into step 6 · *tex:* $\dfrac{dp}{p} = \gamma\,\dfrac{d\rho}{\rho}$ · *why:* Replace dv/v by −dρ/ρ;
     the two minus signs cancel. · *plain:* Relative pressure change = γ × relative density change.
  9. *did:* Integrate both sides with constant γ · *tex:* $\ln p = \gamma\,\ln\rho + \text{const}$ · *why:* ∫dx/x = ln x
     (P36); γ constant because C_p and C_v are (assumption). · *plain:* On log–log axes p against ρ is a straight line of
     slope γ. · *watch (E3):* "the bold curve lies on the faint isentrope".
  10. *did:* Combine the logarithms and exponentiate · *tex:* $\dfrac{p}{\rho^\gamma} = \text{const}$ · *why:* ln p − γ ln ρ
      = ln(p/ρ^γ) by the log and exponent rules (P43); if a logarithm is constant, so is its argument. · *plain:* Along an
      isentrope p/ρ^γ never changes. · *live (E3):* p/ρ^γ at state 1 and at the marker.
- **Result.** $p/\rho^\gamma = \text{const}$ (1.25) — *in words:* with no heat and no friction, a perfect gas with constant
  specific heats keeps p ∝ ρ^γ.
- **Check.** Units: the constant carries Pa (m³/kg)^γ — fine, it is fixed by the starting state ✓. Limit γ → 1: p/ρ = const,
  the isothermal law (RT constant) ✓. Number: pump 1 → 2 bar, ρ₂/ρ₁ = 2^{1/1.4} = 1.641 and T₂ = 288.15 × 2/1.641 = 351.3 K ✓;
  the notebook's step-by-step integration of dp/dρ = γp/ρ agrees with `isentropic_pressure` to 10⁻⁴.
- **What it means.** Work done on the gas has nowhere to go but internal energy, so the gas heats and its pressure rises
  faster than for an isothermal squeeze (slope γ instead of 1). With p = ρRT it gives the ratios (1.26) and hence θ (D20)
  and c = √(γRT) (C47). It fails with friction or heat exchange (stirring, slow compression in a conducting cylinder) and
  when C_p varies strongly with T.
- **Traps.** The sign of dv/v = −dρ/ρ; assuming constant C_p, C_v silently; dividing step 2 by step 3 and getting 1/γ.

### D18 · The parcel equation and N², §1.10 and Eq. (1.29) — ★★, 12 steps, in C50 (notebook · `parcel_stability`)
- **Goal.** Find the equation of motion of a small blob of fluid pushed a little way up or down in a stratified fluid at
  rest, and the number that decides whether it comes back.
- **Start.** $\rho_p V\,\dfrac{d^2\zeta}{dt^2} = \sum F_z$ — *in words:* Newton's second law for a parcel of volume V and
  density ρ_p displaced by ζ from its rest height z_o.
- **Plan.** (1) Forces: weight and buoyancy. (2) Densities of the parcel and of its new surroundings to first order in ζ.
  (3) Keep only first-order terms. (4) Name the coefficient N² and solve.
- **Tools.** Newton's second law (P09) · buoyancy (D37, C20) · first-order Taylor expansion (P26) · linear second-order ODE
  (P44) · square root of a negative number (P45, in C51 — only needed to read the N² < 0 case, restated in step 12).
- **Assumptions.** Small displacement ζ (steps 4–7) · the parcel moves frictionlessly and adiabatically, no mixing (step 4:
  it follows its own isentropic density) · its pressure instantly equals that of its surroundings (step 2: buoyancy from the
  environment) · the background is static (step 5) · released from rest (step 12).
- **Steps.**
  1. *did:* Name the parcel's two vertical forces · *tex:* $\sum F_z = -\rho_p V g + F_b$ · *why:* Its weight points down;
     the surrounding fluid's pressure gives a net force F_b. Nothing else acts on a frictionless parcel. We list forces
     before using Newton's law. · *plain:* Weight down, buoyancy up.
  2. *did:* Use the buoyancy of D37 at the new height · *tex:* $F_b = \rho(z_o + \zeta)\,V g$ · *why:* The parcel's pressure
     matches its surroundings, so the pressure around it is the environment's hydrostatic pressure at z_o + ζ; D37 says the
     net push is the weight of displaced *environment* fluid there. · *plain:* The parcel is buoyed by the fluid it pushes
     aside at its new height. · *set (E4):* t = 0, parcel lifted.
  3. *did:* Insert both forces and divide by ρ_pV · *tex:* $\dfrac{d^2\zeta}{dt^2} = -\,g\,\dfrac{\rho_p - \rho(z_o+\zeta)}
     {\rho_p}$ · *why:* Newton's second law with steps 1–2; dividing by the parcel's mass leaves its acceleration. · *plain:*
     A parcel heavier than its surroundings accelerates down.
  4. *did:* Expand the parcel's density along its own path · *tex:* $\rho_p = \rho(z_o) + \dfrac{d\rho_a}{dz}\,\zeta$ · *why:*
     First-order Taylor (P26). The parcel started with the environment's density ρ(z_o) and changes adiabatically as it
     meets new pressure; dρ_a/dz is that isentropic rate of change. · *plain:* The parcel's density changes at its own rate.
  5. *did:* Expand the environment's density at the new height · *tex:* $\rho(z_o+\zeta) = \rho(z_o) + \dfrac{d\rho}{dz}\,\zeta$
     · *why:* First-order Taylor of the static background profile; dρ/dz is the environment's gradient. · *plain:* The
     surroundings have their own, different rate.
  6. *did:* Subtract the two · *tex:* $\rho_p - \rho(z_o+\zeta) = \Big(\dfrac{d\rho_a}{dz} - \dfrac{d\rho}{dz}\Big)\zeta$ ·
     *why:* ρ(z_o) cancels because parcel and environment matched at rest; what remains is first order in ζ. · *plain:* The
     density difference grows in proportion to the displacement. · *live (E4):* (0 − (−0.01)) × 5 = 0.05 kg/m³.
  7. *did:* Keep only the leading term in the denominator · *tex:* $\rho_p \approx \rho(z_o)$ · *why:* The numerator is already
     proportional to ζ; the correction (dρ_a/dz)ζ in the denominator would add only terms of order ζ², which we drop for a
     small displacement. · *plain:* For the inertia, the parcel's density is its starting density.
  8. *did:* Substitute steps 6 and 7 into step 3 · *tex:* $\dfrac{d^2\zeta}{dt^2} = -\,\dfrac{g}{\rho(z_o)}\Big(\dfrac{d\rho_a}
     {dz} - \dfrac{d\rho}{dz}\Big)\zeta$ · *why:* Replace the density difference and the denominator. · *plain:* The
     acceleration is proportional to the displacement.
  9. *did:* Swap the order in the bracket · *tex:* $\dfrac{d^2\zeta}{dt^2} - \dfrac{g}{\rho(z_o)}\Big(\dfrac{d\rho}{dz} -
     \dfrac{d\rho_a}{dz}\Big)\zeta = 0$ · *why:* −(a − b) = (b − a); then move the right side to the left. This is the book's
     form. · *plain:* The parcel equation of §1.10.
  10. *did:* Name the coefficient N² · *tex:* $N^2 \equiv -\dfrac{g}{\rho(z_o)}\Big(\dfrac{d\rho}{dz} - \dfrac{d\rho_a}{dz}\Big)$
      · *why:* Definition (1.29), chosen so that the equation reads ζ″ + N²ζ = 0 — the oscillator form of P44. · *plain:* One
      number collects gravity and the two density gradients. · *live (E4):* N² = −(9.81/1025)(−0.01 − 0) = 9.57×10⁻⁵ s⁻².
  11. *did:* Try an exponential solution · *tex:* $\zeta = e^{\lambda t}\ \Rightarrow\ \lambda^2 + N^2 = 0$ · *why:* For a
      linear equation with constant coefficients, e^{λt} turns derivatives into powers of λ (P44); e^{λt} ≠ 0 divides out.
      · *plain:* The motion is set by the roots λ = ±√(−N²).
  12. *did:* Apply release from rest · *tex:* $\zeta = \zeta_0\cos Nt\ (N^2>0),\quad \zeta_0\ (N^2=0),\quad \zeta_0\cosh\!
      \sqrt{-N^2}\,t\ (N^2<0)$ · *why:* ζ(0) = ζ₀ and ζ′(0) = 0 pick the even combination of the two roots: cos for
      imaginary roots, cosh for real ones (P45), a constant when both vanish. · *plain:* Oscillate, stay, or run away. ·
      *watch (E4):* "the grey curve has period 2π/N = 10.7 min".
- **Result.** $\dfrac{d^2\zeta}{dt^2} + N^2\zeta = 0,\quad N^2 = -\dfrac{g}{\rho(z_o)}\Big(\dfrac{d\rho}{dz} -
  \dfrac{d\rho_a}{dz}\Big)$ (1.29) — *in words:* a displaced parcel is a harmonic oscillator with squared frequency N²,
  positive when the environment's density falls with height faster than the parcel's own.
- **Check.** Units: (m s⁻²/kg m⁻³) × kg m⁻⁴ = s⁻² ✓. Limits: uniform incompressible fluid with an incompressible parcel,
  dρ/dz = dρ_a/dz = 0 → N² = 0, neutral ✓; environment exactly isentropic, dρ/dz = dρ_a/dz → neutral ✓. Number: thermocline
  dρ/dz = −0.01 kg m⁻⁴, dρ_a/dz = 0, ρ = 1025 → N² = 9.57×10⁻⁵ s⁻², period 642 s ✓ (`brunt_vaisala_sq`,
  `parcel_displacement`, the Euler–Cromer loop and the nonlinear `parcel_ode_from_gradients` agree).
- **What it means.** Stability is a race between two density gradients: the environment's and the parcel's own adiabatic
  one. It is the book's first stability analysis — the same "linearise, look at the sign of the growth rate" logic returns
  in Ch. 11. N is the highest frequency internal waves can have (Ch. 7 §7.8). The model has no friction, so a stable parcel
  rings forever; real fluids damp it. It fails for large displacements (the Taylor steps) and when the parcel mixes.
- **Traps.** The sign of buoyancy; using the environment's gradient for the parcel; keeping O(ζ²) terms inconsistently;
  replacing ρ_p by ρ(z_o) *before* subtracting (that would lose the whole effect).

### D19 · The adiabatic lapse rate without perfect-gas relations, Eq. (1.30) — ★★★, 16 steps + sympy, in C54 (notebook · `parcel_stability`)
- **Goal.** Find the rate at which the temperature of a parcel changes as it moves isentropically up through a fluid in
  hydrostatic balance — for *any* fluid, not just a perfect gas — and then read it in both lapse-rate sign conventions.
- **Start.** $T\,ds = dh - v\,dp$ (1.18) with $ds = 0$ and $\dfrac{dp}{dz} = -\rho g$ (1.8) — *in words:* the parcel's
  entropy stays constant, and the pressure it feels is the hydrostatic pressure of its surroundings.
- **Plan.** (1) Write dh for the isentropic parcel. (2) Expand dh in T and p, which brings in C_p and an unknown
  (∂h/∂p)_T. (3) Find (∂h/∂p)_T with the Gibbs free energy, a Maxwell relation and α. (4) Divide by dz and use hydrostatics.
- **Tools.** Gibbs relations (C35, D10) · C_p (1.14) (C30) · α (1.20) (C37) · two properties fix the state (C28) ·
  hydrostatic law (C20, D05) · chain rule (P49) · partial derivative with a variable held fixed (P39) · Gibbs free energy
  (P50) · exact differentials and Maxwell relations (P51) · inequalities under a sign change (P48, for the result).
- **Assumptions.** Isentropic parcel, ds = 0 (step 1) · single-component fluid whose state is fixed by T and p (step 2) ·
  smooth state functions, so mixed partials commute (step 7) · parcel pressure = environment pressure and the environment
  is hydrostatic (step 14) · evaluated where the parcel still has the environment's density (step 14; see Traps) · perfect
  gas only in the last step (step 16).
- **Steps.**
  1. *did:* Set ds = 0 in the second Gibbs form · *tex:* $dh = v\,dp$ · *why:* The parcel moves without heat or friction, so
     its entropy is constant (T ds = 0). We use the enthalpy form because it contains dp, which hydrostatics will supply. ·
     *plain:* The parcel's enthalpy changes only through the pressure it meets.
  2. *did:* Expand dh in T and p · *tex:* $dh = \Big(\dfrac{\partial h}{\partial T}\Big)_p dT + \Big(\dfrac{\partial h}
     {\partial p}\Big)_T dp$ · *why:* Two properties fix the state (C28), so h = h(T, p); the chain rule for a function of two
     variables (P49) gives its differential. We want dT to appear. · *plain:* h changes because T changes and because p
     changes.
  3. *did:* Recognise C_p · *tex:* $dh = C_p\,dT + \Big(\dfrac{\partial h}{\partial p}\Big)_T dp$ · *why:* (1.14) defines
     C_p = (∂h/∂T)_p. One coefficient is now a known property; the other still needs work. · *plain:* The first part is the
     specific heat times the temperature change.
  4. *did:* Get (∂h/∂p)_T from Gibbs at fixed T · *tex:* $\Big(\dfrac{\partial h}{\partial p}\Big)_T = T\Big(\dfrac{\partial s}
     {\partial p}\Big)_T + v$ · *why:* Gibbs dh = T ds + v dp holds for any change (D10); take a change at constant T and
     divide by dp (P39). This trades the unknown for an entropy derivative. · *plain:* How h depends on p at fixed T involves
     how s does.
  5. *did:* Introduce the Gibbs free energy · *tex:* $dg = v\,dp - s\,dT,\quad g \equiv h - Ts$ · *why:* dg = dh − T ds − s dT
     (product rule on Ts, P38) and dh = T ds + v dp; the T ds terms cancel (P50). We choose g because its natural variables
     are exactly p and T. · *plain:* A helper state function whose changes involve only dp and dT.
  6. *did:* Read off the partial derivatives of g · *tex:* $\Big(\dfrac{\partial g}{\partial p}\Big)_T = v,\quad
     \Big(\dfrac{\partial g}{\partial T}\Big)_p = -s$ · *why:* Comparing dg = v dp − s dT with the chain rule dg = (∂g/∂p)_T
     dp + (∂g/∂T)_p dT term by term (P51). · *plain:* v and −s are the slopes of g.
  7. *did:* Equate the mixed second derivatives · *tex:* $\Big(\dfrac{\partial v}{\partial T}\Big)_p = -\Big(\dfrac{\partial s}
     {\partial p}\Big)_T$ · *why:* g is a smooth state function, so ∂²g/∂T∂p = ∂²g/∂p∂T (equality of mixed partials, P51);
     differentiating step 6's first result in T and the second in p gives this **Maxwell relation**. It converts an entropy
     derivative into a measurable one. · *plain:* How entropy changes with pressure is set by thermal expansion.
  8. *did:* Express (∂v/∂T)_p with α · *tex:* $\Big(\dfrac{\partial v}{\partial T}\Big)_p = v\,\alpha$ · *why:* v = 1/ρ gives
     (∂v/∂T)_p = −(1/ρ²)(∂ρ/∂T)_p; by (1.20) (∂ρ/∂T)_p = −ρα, so the result is α/ρ = vα. · *plain:* Heating at fixed pressure
     swells unit mass by vα per kelvin.
  9. *did:* Insert into the Maxwell relation · *tex:* $\Big(\dfrac{\partial s}{\partial p}\Big)_T = -\,v\,\alpha$ · *why:*
     Steps 7 and 8. The sign matters: squeezing at fixed T lowers entropy when α > 0. · *plain:* Entropy falls with pressure
     at a rate vα.
  10. *did:* Complete (∂h/∂p)_T · *tex:* $\Big(\dfrac{\partial h}{\partial p}\Big)_T = v\,(1 - \alpha T)$ · *why:* Substitute
      step 9 into step 4: −Tvα + v. · *plain:* The pressure part of dh, now in measurable properties.
  11. *did:* Equate the two expressions for dh · *tex:* $C_p\,dT + v(1 - \alpha T)\,dp = v\,dp$ · *why:* Step 3 (with step 10)
      and step 1 describe the same change of the same parcel. · *plain:* Two bookkeepings of one enthalpy change.
  12. *did:* Subtract v dp from both sides · *tex:* $C_p\,dT = v\,\alpha\,T\,dp$ · *why:* v(1 − αT)dp − v dp = −vαT dp; move it
      across. Temperature change is now tied to pressure change alone. · *plain:* A parcel warms when compressed and cools
      when the pressure drops (for α > 0).
  13. *did:* Divide by dz along the parcel's path · *tex:* $C_p\,\dfrac{dT_a}{dz} = v\,\alpha\,T\,\dfrac{dp}{dz}$ · *why:*
      Both sides are changes along the vertical path; dividing by dz turns them into rates. The subscript a marks the
      adiabatic parcel. · *plain:* Temperature rate per metre ↔ pressure rate per metre.
  14. *did:* Use hydrostatics with v = 1/ρ · *tex:* $C_p\,\dfrac{dT_a}{dz} = -\,\alpha\,T\,g$ · *why:* The parcel's pressure
      is the environment's, dp/dz = −ρg (1.8); at the release height the parcel's v equals the environment's 1/ρ, so vρ = 1.
      · *plain:* The pressure drop with height sets the cooling. · *highlight (E4):* readout Γa.
  15. *did:* Divide by C_p · *tex:* $\dfrac{dT_a}{dz} \equiv \Gamma_a = -\,\dfrac{g\,\alpha\,T}{C_p}$ · *why:* C_p > 0. This is
      (1.30), derived without any perfect-gas relation. · *plain:* The adiabatic lapse rate of any fluid. · *live (E4):*
      −9.807 × (1/288.15) × 288.15/1004.7 = −9.76×10⁻³ K/m.
  16. *did:* Specialise to a perfect gas · *tex:* $\Gamma_a = -\,\dfrac{g}{C_p} = -9.76\ \text{K/km}$ · *why:* α = 1/T for a
      perfect gas (1.28, C48), so αT = 1. With g = 9.807 m/s² and C_p = 1004.7 J kg⁻¹ K⁻¹. · *plain:* Dry air cools 9.76 K
      for every kilometre it rises. · *watch (E4):* "the dashed adiabat's slope is exactly this".
- **Result.** $\Gamma_a = \dfrac{dT_a}{dz} = -\dfrac{g\alpha T}{C_p}$ (1.30); dry air $\Gamma_a = -g/C_p = -9.76$ K/km.
  *In words, in both conventions:* **Kundu (Γ ≡ dT/dz):** Γ_a = −9.76 K/km and the environment is stable when
  **Γ > Γ_a** (e.g. −6.5 > −9.76). **Meteorology (Γ ≡ −dT/dz):** Γ_a,met = +gαT/C_p = +9.76 K/km and stable when
  **Γ < Γ_a** (6.5 < 9.76). Converting multiplies both sides by −1, which reverses the inequality (P48); the verdict is
  identical.
- **Check.** Units: m s⁻² · K⁻¹ · K / (J kg⁻¹ K⁻¹) = K/m ✓. Limits: α = 0 (a fluid that does not expand) → Γ_a = 0, no
  adiabatic temperature change ✓; perfect gas → −g/C_p ✓. Numbers: air −9.76 K/km; water at 283 K (α = 1.5×10⁻⁴ K⁻¹,
  C_p = 4190) → −0.10 K/km ✓ (`adiabatic_lapse_rate(T=283.15, cp=4190.0, alpha=1.5e-4)`). The from-scratch cell
  differentiates a lifted parcel's (1.26) temperature and finds −g/C_p at the release height ✓.
- **sympy check intent** (`check_src`, required ★★★): verify steps 7–15 for an arbitrary fluid by starting from a general
  Gibbs free energy G(T, p): v = ∂G/∂p, s = −∂G/∂T, h = G + Ts, C_p = ∂h/∂T, α = (∂v/∂T)/v; the slope of an isentrope
  dT/dp at constant s = −(∂s/∂p)/(∂s/∂T); assert it equals vαT/C_p (step 12), multiply by dp/dz = −g/v and assert −gαT/C_p
  (step 15); then substitute the perfect-gas G = C_p(T − T ln T) + RT ln p and assert −g/C_p (step 16). Sketch (tested
  2026-09-12, returns 0, 0, −g/c_p):
  ```python
  import sympy as sp                                         # symbolic algebra
  T, p, g, R, cp = sp.symbols("T p g R c_p", positive=True)  # temperature, pressure, gravity, gas constant, C_p
  G = sp.Function("G")(T, p)                                 # ANY smooth Gibbs free energy per unit mass
  v = sp.diff(G, p)                                          # step 6: v = (∂g/∂p)_T
  s = -sp.diff(G, T)                                         # step 6: s = −(∂g/∂T)_p
  h = G + T*s                                                # g = h − Ts rearranged
  Cp = sp.diff(h, T)                                         # (1.14): C_p = (∂h/∂T)_p
  alpha = sp.diff(v, T)/v                                    # step 8: vα = (∂v/∂T)_p
  dTdp_s = -sp.diff(s, p)/sp.diff(s, T)                      # isentrope slope (implicit-function rule)
  print(sp.simplify(dTdp_s - v*alpha*T/Cp))                  # 0 → step 12 holds for any fluid
  print(sp.simplify(dTdp_s*(-g/v) + g*alpha*T/Cp))           # 0 → (1.30) after dp/dz = −g/v (step 14)
  Gpg = cp*(T - T*sp.log(T)) + R*T*sp.log(p)                 # a perfect gas with constant C_p
  print(sp.simplify((-g*alpha*T/Cp).subs(G, Gpg).doit()))    # −g/c_p → step 16
  ```
- **What it means.** Γ_a is the slope a well-mixed (isentropic) layer settles to — the reference against which every
  observed profile is judged. A column cooling more slowly than 9.76 K/km (dT/dz > −9.76, or Γ_met < 9.76) is stable even
  though it cools with height. It fails for moist air (latent heat makes the saturated rate ≈ 4–7 K/km), for large
  excursions where the parcel's temperature departs from the environment's (then the local rate is −(g/C_p)T_p/T_env), and
  for mixtures whose composition changes (seawater salinity).
- **Traps.** The Maxwell-relation sign; using dh = C_p dT (true only for a perfect gas); using the parcel's v with the
  environment's ρ away from the release point; **negating Γ without flipping the inequality** when switching to the
  meteorology convention.

### D20 · Potential temperature, Eq. (1.31) — ★, 5 steps, in C55 (notebook · `parcel_stability`)
- **Goal.** Define a temperature label that a parcel keeps when it moves up or down without heating: the temperature it
  would have if brought adiabatically to a standard pressure p_o.
- **Start.** $\dfrac{T}{T_0} = \Big(\dfrac{p}{p_0}\Big)^{(\gamma-1)/\gamma}$ (1.26) — *in words:* along one isentrope of a
  perfect gas, the temperature ratio of any two states is a power of their pressure ratio.
- **Plan.** (1) Choose the two states as "now" and "at the reference pressure". (2) Solve for the reference temperature.
  (3) Rewrite the exponent with R/C_p.
- **Tools.** Isentropic ratios (C46, C45) · exponent rules (P43) · R = C_p − C_v and γ = C_p/C_v (C42, C43).
- **Assumptions.** Perfect gas with constant γ (start) · the imagined trip to p_o is adiabatic and reversible (step 1) ·
  p_o is a fixed reference, 1000 hPa (step 1; not the surface pressure p₀).
- **Steps.**
  1. *did:* Take the reference state (θ, p_o) and the current state (T, p) · *tex:* $\dfrac{T}{\theta} = \Big(\dfrac{p}
     {p_o}\Big)^{(\gamma-1)/\gamma}$ · *why:* (1.26) links any two states on one isentrope; by definition the parcel taken
     adiabatically to p_o has temperature θ, so rename T₀ → θ and p₀ → p_o. · *plain:* The parcel now and the parcel at
     1000 hPa lie on the same isentrope.
  2. *did:* Multiply both sides by θ · *tex:* $T = \theta\,\Big(\dfrac{p}{p_o}\Big)^{(\gamma-1)/\gamma}$ · *why:* θ ≠ 0; this
     is the book's form (1.31). · *plain:* Actual temperature = potential temperature scaled by a pressure factor. · *live
     (E4):* T at z_o from θ and p.
  3. *did:* Solve for θ · *tex:* $\theta = T\,\Big(\dfrac{p_o}{p}\Big)^{(\gamma-1)/\gamma}$ · *why:* Divide by the power and use
     1/(a/b)^k = (b/a)^k (P43). This is the form we compute. · *plain:* Correct the thermometer reading for the pressure.
  4. *did:* Rewrite the exponent with C_v and C_p · *tex:* $\dfrac{\gamma-1}{\gamma} = 1 - \dfrac{C_v}{C_p}$ · *why:* Split the
     fraction: (γ − 1)/γ = 1 − 1/γ, and 1/γ = C_v/C_p by (1.24). · *plain:* The exponent is one minus the specific-heat ratio's
     inverse.
  5. *did:* Use R = C_p − C_v · *tex:* $\theta = T\,\Big(\dfrac{p_o}{p}\Big)^{R/C_p}$ · *why:* 1 − C_v/C_p = (C_p − C_v)/C_p
     = R/C_p by (1.23); for air R/C_p = 287.06/1004.7 = 2/7. · *plain:* The form used in meteorology, exponent 0.286. ·
     *live (E4):* θ = … K.
- **Result.** $T = \theta\,(p/p_o)^{(\gamma-1)/\gamma}$ (1.31), equivalently $\theta = T\,(p_o/p)^{R/C_p}$ — *in words:*
  potential temperature is the temperature a parcel would have at 1000 hPa.
- **Check.** The exponent is dimensionless and so is p_o/p ✓. At p = p_o, θ = T ✓. Adiabatic motion keeps θ constant, since
  it only moves the parcel along its isentrope ✓. Number: 500 hPa, 250 K → θ = 250 × 2^{2/7} = 304.8 K ✓
  (`potential_temperature(250.0, 5e4)` and the from-scratch formula agree).
- **What it means.** θ is the conserved temperature of dry adiabatic motion: two parcels can be compared regardless of height.
  Its vertical gradient is the stability test (C57, D36). It fails when heat is added (radiation, condensation — then θ
  changes, which is how diabatic heating is measured).
- **Traps.** Confusing p_o (1000 hPa, fixed) with p₀ (the surface pressure); inverting the ratio (p/p_o instead of p_o/p)
  in the computed form.

### D36 · N² from the potential-temperature gradient — ★★, 8 steps, in C55 (notebook · `parcel_stability`)
- **Goal.** Show that for a perfect-gas atmosphere the density-based N² of (1.29) is simply g/θ times the vertical gradient of
  θ — so stability is "θ increases upward" — and connect it to the two lapse rates.
- **Start.** $N^2 = -\dfrac{g}{\rho}\Big(\dfrac{d\rho}{dz} - \dfrac{d\rho_a}{dz}\Big)$ (1.29) — *in words:* the parcel
  result of D18, evaluated at the height where parcel and environment share the density ρ.
- **Plan.** (1) Write both density gradients as relative (logarithmic) rates. (2) Environment from p = ρRT; parcel from the
  isentropic ratio. (3) Subtract: the pressure terms combine into θ. (4) Use (1.32) for the lapse-rate form.
- **Tools.** (1.29) (C51, D18) · logarithmic differentiation (P52) · p = ρRT (C40) · isentropic ratios (1.26) (C46) · θ
  (1.31) (D20) · (1.32) (C56) · γ and R/C_p (C42, C43).
- **Assumptions.** Perfect gas with constant γ (steps 2–3) · the parcel's pressure equals the environment's (step 3) · at the
  parcel's rest height ρ_a = ρ (step 1).
- **Steps.**
  1. *did:* Take ρ inside the bracket · *tex:* $N^2 = -g\Big(\dfrac1\rho\dfrac{d\rho}{dz} - \dfrac1\rho\dfrac{d\rho_a}{dz}\Big)$
     · *why:* Distribute 1/ρ over both terms; at the rest height the parcel's density equals ρ, so each term is a relative
     (logarithmic) rate. Relative rates are what logs make easy. · *plain:* N² compares relative density changes.
  2. *did:* Log-differentiate the environment's p = ρRT · *tex:* $\dfrac1\rho\dfrac{d\rho}{dz} = \dfrac1p\dfrac{dp}{dz} -
     \dfrac1T\dfrac{dT}{dz}$ · *why:* ρ = p/(RT) so ln ρ = ln p − ln R − ln T; differentiate in z with R constant (P52). ·
     *plain:* Environment density falls because pressure falls, and rises where it gets colder.
  3. *did:* Log-differentiate the parcel's isentropic density · *tex:* $\dfrac1{\rho_a}\dfrac{d\rho_a}{dz} = \dfrac1\gamma\,
     \dfrac1p\dfrac{dp}{dz}$ · *why:* Along its isentrope ρ_a ∝ p^{1/γ} (1.26), and the parcel's pressure is the
     environment's p(z); ln ρ_a = (1/γ) ln p + const. · *plain:* The parcel's density follows the pressure only, more weakly.
  4. *did:* Subtract step 3 from step 2 · *tex:* $\dfrac1\rho\dfrac{d\rho}{dz} - \dfrac1\rho\dfrac{d\rho_a}{dz} = \Big(1 -
     \dfrac1\gamma\Big)\dfrac1p\dfrac{dp}{dz} - \dfrac1T\dfrac{dT}{dz}$ · *why:* Collect the two pressure terms; ρ_a = ρ at the
     rest height. · *plain:* What decides stability: a pressure part and a temperature part.
  5. *did:* Log-differentiate θ = T(p_o/p)^{(γ−1)/γ} · *tex:* $\dfrac1\theta\dfrac{d\theta}{dz} = \dfrac1T\dfrac{dT}{dz} -
     \dfrac{\gamma-1}{\gamma}\,\dfrac1p\dfrac{dp}{dz}$ · *why:* ln θ = ln T + ((γ − 1)/γ)(ln p_o − ln p) (P52); p_o is constant.
     This is the book's log-derivative of (1.31) (N19). · *plain:* θ rises with T and with falling pressure.
  6. *did:* Recognise the bracket of step 4 · *tex:* $\Big(1 - \dfrac1\gamma\Big)\dfrac1p\dfrac{dp}{dz} - \dfrac1T\dfrac{dT}{dz}
     = -\,\dfrac1\theta\dfrac{d\theta}{dz}$ · *why:* 1 − 1/γ = (γ − 1)/γ, so step 4 is exactly minus step 5. · *plain:* The
     density race is the θ gradient in disguise.
  7. *did:* Substitute into step 1 · *tex:* $N^2 = \dfrac{g}{\theta}\,\dfrac{d\theta}{dz}$ · *why:* −g × (−(1/θ)dθ/dz); the
     minus signs cancel. · *plain:* N² is g over θ times θ's gradient. · *live (E4):* (9.807/θ) × dθ/dz = …; *watch:* "equals
     the N² readout computed from the lapse rate".
  8. *did:* Use (1.32) for the lapse-rate form · *tex:* $N^2 = \dfrac{g}{T}\Big(\dfrac{dT}{dz} + \dfrac{g}{C_p}\Big) =
     \dfrac{g}{T}\,(\Gamma - \Gamma_a)$ · *why:* (1.32) gives (1/θ)dθ/dz = (1/T)(dT/dz + g/C_p) (C56); with Γ_a = −g/C_p
     (D19) the bracket is Γ − Γ_a in Kundu's convention. · *plain:* N² is proportional to how much slower the environment
     cools than a parcel. · *live (E4):* (g/T)(Γ − Γ_a) = ….
- **Result.** $N^2 = \dfrac{g}{\theta}\dfrac{d\theta}{dz} = \dfrac{g}{T}\,(\Gamma - \Gamma_a)$ (Kundu) $= \dfrac{g}{T}\,
  (\Gamma_{a,\rm met} - \Gamma_{\rm met})$ (meteorology) — *in words:* a perfect-gas atmosphere is stable exactly where θ
  increases upward; in lapse rates, where the environment cools more slowly than 9.76 K/km.
- **Check.** Units: (m s⁻²/K) × K/m = s⁻² ✓. Neutral: dθ/dz = 0 ⇔ Γ = Γ_a ✓. Isothermal layer (dT/dz = 0) at 250 K:
  N² = g²/(C_pT) = 3.83×10⁻⁴ s⁻², period 5.4 min ✓. Standard troposphere at 288.15 K: 1.11×10⁻⁴ s⁻², period 9.9 min ✓ — the
  notebook computes N² from θ(z) of USSA-1976 with `np.gradient` and from the lapse rate and they agree to 10⁻³.
- **What it means.** One profile — θ(z) — tells you everything about dry static stability; this is how Ch. 13 and weather
  analysis read soundings, and why isentropic surfaces are the natural coordinates of dry dynamics. It fails for moist air
  (use equivalent potential temperature) and for the ocean, where salinity requires potential density (C61).
- **Traps.** Using the environment's density gradient for the parcel; forgetting (γ − 1)/γ = R/C_p when comparing with
  (1.32); trying to substitute dp/dz = −ρg too early (it is not needed until step 8).

### D28 · Buckingham's Π theorem, Eq. (1.37) — ★★★, 12 steps + sympy, in C69 (notebook · `buckingham_pi_machine`)
- **Goal.** Show that any correct relation among n dimensional variables can be rewritten as a relation among exactly n − r
  independent dimensionless groups, where r is the rank of the dimensional matrix — and why the groups are not unique.
- **Start.** $f(q_1, q_2, \dots, q_n) = 0$ (1.36) with dimensional matrix $A$ (rows = base dimensions, columns = variables,
  rank $r$) — *in words:* a complete, dimensionally homogeneous relation among n variables whose dimensions we have
  tabulated.
- **Plan.** (1) A product of powers is dimensionless exactly when its exponent vector solves A·k = 0. (2) Rank–nullity
  counts the independent solutions: n − r. (3) Choosing units cleverly shows that the relation can depend only on those
  groups.
- **Tools.** Dimensional homogeneity (C64) · the dimensional matrix and its rank (C67, C68) · exponent rules (P43) ·
  linear independence and rank (P54) · null space and rank–nullity (P58) · scaling the base units (P59) · matrices and
  minors (P53).
- **Assumptions.** The variable list is complete (step 12: nothing else can enter f) · f is dimensionally homogeneous, i.e.
  true in every unit system (step 12) · all dimensions are products of powers of the base dimensions (step 2).
- **Steps.**
  1. *did:* Form a general product of powers · *tex:* $\Pi = q_1^{k_1}\,q_2^{k_2}\cdots q_n^{k_n}$ · *why:* Only products of
     powers of the variables have dimensions that are again products of powers of M, L, T, Θ, so they are the candidates for
     unit-free combinations. The exponents k = (k₁, …, k_n) are what we must find. · *plain:* A group is some variables
     multiplied together, each raised to a power.
  2. *did:* Write its dimension with the matrix entries · *tex:* $[\Pi] = \prod_d d^{\,\sum_j A_{dj} k_j}$ · *why:* [q_j] =
     M^{A_Mj} L^{A_Lj} T^{A_Tj} Θ^{A_Θj}; raising to k_j multiplies exponents and multiplying quantities adds them (P43). ·
     *plain:* The exponent of each base dimension is a weighted sum of the k's.
  3. *did:* Require every exponent to vanish · *tex:* $A\,\mathbf k = \mathbf 0$ · *why:* Dimensionless means M⁰L⁰T⁰Θ⁰;
     each row d of A gives one equation Σ_j A_dj k_j = 0. · *plain:* A group is dimensionless exactly when its exponent vector
     solves this linear system. · *live (E5):* for Π₁, k = (1, 0, 0, 0, −2, −1, 0): M 1 − 1 = 0, L −1 − 2 + 3 = 0, T −2 + 2 = 0.
  4. *did:* Name the solution set · *tex:* $\{\mathbf k : A\mathbf k = \mathbf 0\} = \operatorname{null}(A)$ · *why:*
     Definition of the null space (P58); sums and multiples of solutions are solutions, so it is a vector space. · *plain:*
     All dimensionless groups live in the null space of A.
  5. *did:* Count its dimension · *tex:* $\dim\operatorname{null}(A) = n - r$ · *why:* Rank–nullity theorem for a matrix with n
     columns: rank + nullity = n (P58). · *plain:* There are exactly n − r independent exponent vectors. · *live (E5):* 7 − 3
     = 4.
  6. *did:* Pick a basis and call the groups Π₁…Π_{n−r} · *tex:* $\mathbf k^{(1)},\dots,\mathbf k^{(n-r)} \;\to\; \Pi_1,\dots,
     \Pi_{n-r}$ · *why:* Any basis of the null space will do; its vectors are linearly independent (P54), so no group is a
     product of powers of the others (C71). Different bases give different, equivalent sets — groups are not unique. ·
     *plain:* n − r independent dimensionless numbers.
  7. *did:* Show every group is built from these · *tex:* $\mathbf k = \sum_i c_i\mathbf k^{(i)} \;\Rightarrow\; \Pi = \prod_i
     \Pi_i^{\,c_i}$ · *why:* A basis spans the space; adding exponent vectors multiplies the corresponding powers (P43). ·
     *plain:* Any other dimensionless combination is a product of powers of Π₁…Π_{n−r}.
  8. *did:* Choose r variables with independent columns · *tex:* $\det A_{\rm rep} \neq 0$ · *why:* Rank r means some r×r minor
     is nonzero (C68); its r variables are the repeating set. If A has more rows than r, drop the dependent rows first. · *plain:*
     Pick r variables that between them carry all the independent dimensions. · *highlight (E5):* the outlined minor.
  9. *did:* Rescale the base units · *tex:* $q_j \to q_j\prod_d \lambda_d^{A_{dj}}$ · *why:* Changing the unit of dimension d
     by a factor λ_d multiplies the *number* of every quantity by λ_d to its exponent (P59); physics is unchanged. · *plain:*
     New units, new numbers, same world.
  10. *did:* Choose the scales that make the repeating variables equal to 1 · *tex:* $\sum_d A_{dj}\ln\lambda_d = -\ln q_j
      \quad (j \in \text{rep})$ · *why:* Taking logs turns step 9 into r linear equations for the ln λ_d whose matrix is the
      (transposed) repeating block — nonsingular by step 8 — so a solution exists. · *plain:* We can always pick units in
      which the r repeating variables are exactly one. · *set (E5):* units cgs; *watch:* "every Π value is unchanged".
  11. *did:* Read the other variables in those units · *tex:* $q_j' = \Pi_j \quad (j \notin \text{rep})$ · *why:* Π_j = q_j
      × (powers of the repeating variables) is dimensionless, so its number is the same in every unit system (C64); in the new
      units the repeating factors are 1, so q_j′ equals Π_j. · *plain:* In these units each remaining variable *is* its group.
  12. *did:* Use homogeneity to drop the ones · *tex:* $f(1,\dots,1,\Pi_1,\dots,\Pi_{n-r}) = 0 \;\Rightarrow\; \phi(\Pi_1,\dots,
      \Pi_{n-r}) = 0$ · *why:* A correct law holds in every unit system (C64), in particular in these; the r constant ones
      are absorbed into a new function φ. Completeness of the variable list guarantees nothing else hides in f. · *plain:*
      The law depends only on the n − r groups. · *live (E5):* the law for your problem.
- **Result.** $\phi(\Pi_1, \Pi_2, \dots, \Pi_{n-r}) = 0$ or $\Pi_1 = \varphi(\Pi_2, \dots, \Pi_{n-r})$ (1.37) — *in
  words:* n variables with r independent dimensions can always be traded for n − r independent dimensionless groups.
- **Check.** Pipe: n = 7, r = 3 → 4 groups (1.40) ✓; Π₁ = Δp/ρU² = 10 in SI and in cgs ✓; `groups_independent` rejects a
  fifth group ✓; Example 1.4 (n = 4, r = 3) → one group, hence a constant ✓.
- **sympy check intent** (`check_src`, required ★★★): with the exact (1.39) matrix, compute `A.rank()` and `A.nullspace()`;
  assert the nullspace has n − r = 4 vectors, each with `A*k` zero; assert Π₁'s exponent vector (1, 0, 0, 0, −2, −1, 0) lies
  in the null space and in the span of the basis (rank does not grow when it is appended). Sketch (tested 2026-09-12: rank 3,
  nullity 4, True, True, True):
  ```python
  import sympy as sp                                          # exact linear algebra
  A = sp.Matrix([[1, 0, 0, 0, 0, 1, 1],                       # M row of (1.39): dp dx d eps U rho mu
                 [-1, 1, 1, 1, 1, -3, -1],                    # L row
                 [-2, 0, 0, 0, -1, 0, -1]])                   # T row
  ns = A.nullspace()                                          # basis of all dimensionless exponent vectors (step 4)
  print(A.rank(), len(ns), A.shape[1] - A.rank())             # 3 4 4 → rank–nullity (step 5)
  print(all((A*k).is_zero_matrix for k in ns))                # True → every basis vector is dimensionless (step 3)
  k1 = sp.Matrix([1, 0, 0, 0, -2, -1, 0])                     # exponents of Π₁ = Δp U⁻² ρ⁻¹
  print((A*k1).is_zero_matrix)                                # True → Π₁ is dimensionless
  print(sp.Matrix.hstack(*ns, k1).rank() == len(ns))          # True → Π₁ is a combination of the basis (step 7)
  ```
- **What it means.** Dimensional analysis is linear algebra: groups are null-space vectors, their number is rank–nullity, a
  different repeating set is a change of basis. It is why experiments are plotted against Re, Fr, Ro and why a small model
  can stand for a big system (Ch. 4 §4.11). It fails — silently — if a relevant variable is missing (drop μ from the pipe:
  3 groups, no Reynolds number, and laminar data will not collapse) or if the "law" mixes quantities that are only
  numerically related in one unit system.
- **Traps.** Believing the groups are unique; forgetting that the relation must be complete and homogeneous; choosing a
  singular repeating set (two pure lengths such as d and ε cannot cancel time).
