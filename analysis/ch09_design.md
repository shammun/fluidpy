# Chapter 9 — Boundary Layers and Related Topics: lesson design
(from `analysis/ch09_curation.md` (A 14 · B 118 · C 16 = 148 rows; CORE 14 · NOTE 124 · RECAP 9 · SKIP 1; 22 derivations
D01–D22 written out (★ 5 · ★★ 12 · ★★★ 5), 219 steps; E1–E9 + backup B1) and `analysis/ch09.md` (§2b derivations, §4
implementation rows, §5 core modules `core/boundary_layer.py`, `core/jets.py`, `core/bluff_body.py`, §9 conventions and
slips); every equation re-read on the rendered pages `chapters/pages/ch09/` (printed = pdf − 27): (9.7)–(9.12) p392
(printed 365, the missing squares of slip R1 confirmed on the image), (9.36) and the Fig. 9.7 labels p401, (9.42)–(9.43)
p403, Fig. 9.8 and (9.50) p406, (9.76) and the Reynolds numbers p431, (9.79)–(9.80) p432, (9.81)–(9.84) p433; the rest
from the analyst's page-image transcription (analysis §2). 2026-09-29, lesson-designer. The implementer works in
parallel from analysis §4 + curation §9; **Part C is written first and is the contract both sides keep.** New numbers in this
design were computed with scratch scripts (scipy `solve_bvp`, Töpfer IVP, sympy); each is marked *expect* so the builder
can compare an executed cell against it.)

**Binding conventions for every builder (analysis §9, curation decisions 1–6 and §8).**
1. **Imports and aliases.** `from fluidpy import ch09_boundary_layers as ch09`; `from fluidpy.core import boundary_layer as
   BL, jets as JET, bluff_body as BB`. **`ch09` re-exports every public name of `BL`, `JET`, `BB`, of
   `core.similarity_reduce`, of `core.similarity` (`sphere_drag_coefficient`, `strouhal_number`) and of `core.creeping` (`stokes_drag_coefficient`, `oseen_drag_coefficient`)**, so explainer parity
   rows write `ch09.<name>` only. Units SI: x, y, δ, θ, L, d, a, b [m]; U, u, v [m/s]; ν [m²/s]; μ [Pa s]; p, τ [Pa];
   dp/dx [Pa/m]; J [N/m] (momentum flux per unit span); ṁ [kg/(m s)]; Γ [m²/s]; Ψ [m⁴/s³]; frequencies f [Hz],
   Ω [rad/s]. Angles enter and leave as `phi_deg` measured **from the forward stagnation point** (ch06 measured θ from
   the downstream axis; ch08's sphere likewise) — the docstring and every axis label say so.
2. **Defaults.** Air: ν = 1.5 × 10⁻⁵ m²/s, ρ = 1.2 kg/m³ (μ = 1.8 × 10⁻⁵ Pa s); water: ν = 1.0 × 10⁻⁶ m²/s, ρ = 1000 kg/m³;
   g = `G_BOOK` = 9.81 (ball swing only). Blasius constants are **computed**, never typed (f″(0) = 0.3320573362 from the
   Töpfer IVP, agreeing with `solve_bvp` to < 1e-9); the notebook prints 4 significant figures (0.3321, 4.910, 1.721,
   0.6641, 2.591, 0.8604, 1.328). `np.trapezoid`, never `np.trapz`.
3. **Scalar-callable and parity-friendly.** Every public function accepts Python floats and returns a float, a tuple of
   floats or a `dict` of floats (arrays only when the input is an array). Parity `py:` expressions use `ch09.…`, `np.pi`,
   numbers, strings, lists, dicts and keyword arguments only, indexed down to one float. **`Viz.num.erf` is not needed in
   this chapter** (no erf), but JS explainers that print 1e-10 parity numbers write their own RK4 (`Viz.num.rk4Step`)
   with a step ≤ 0.005 in η.
4. **Signs and coordinates** (⚠️ callouts where first used, with numbers): x along the surface, y normal; **u = ∂ψ/∂y,
   v = −∂ψ/∂x**; the book's **dp/dx** (favourable < 0, adverse > 0; ch08's `dpdx` sign); wedge exponent **n** in
   U_e = a xⁿ (the FS `m` argument of code is the same n; β = 2n/(n+1) is the Hartree parameter); Thwaites' **λ**
   (dimensionless) versus the wavelength λ of Ch. 7; θ is the momentum thickness here but an angle elsewhere in the
   book; L(λ) the Thwaites correlation versus L the body length; **six Reynolds numbers** (Re = U∞L/ν of (9.6) ·
   Re_x = Ux/ν · Re_L = UL/ν · jet Re_x = x u₀/ν and Re_h99 · cylinder/sphere Re = U∞d/ν on the **diameter** · FS
   Re_x = a x^{n+1}/ν) — every figure axis says which. η = y/δ with δ = √(νx/U) for Blasius (a factor 1 from ch08's
   η = y/√(νt) under t ↔ x/U; ch08's figures plotted y/(2√(νt))).
5. **Book slips taught in corrected form** (each gets "the book prints X; the correct form is Y" where it is used, a
   discriminating wrong-variant option in code, never asserted as physics). The table (analysis §9 labels R1–R16 = curation
   §8 labels S-1 … S-16):

| Slip | Printed | Correct | Taught in | Code option that must fail |
|---|---|---|---|---|
| R1 (9.7) | denominators $\partial x^*$, $\partial y^*$ without squares | $\frac1{Re}\frac{\partial^2u^*}{\partial x^{*2}}+\frac{\partial^2u^*}{\partial y^{*2}}$ | C01 (N08, D01 step 9) | `ch09.bl_nondim_sympy(printed_9_7=True)` fails the dimension check |
| R2 (9.30) | η₉₉ = 4.93 | root of f′ = 0.99 is 4.910 (0.4 % smaller) | C04 (N35, D06 step 4) | `BL.blasius_delta99(…, printed=True)` |
| R3 wall-jet ODE below (9.82) | $f'''+ff''+2f'^2=0$ | $4f'''+ff''+2f'^2=0$ | C13 (N120, D20 step 11) | `JET.wall_jet_ode_solve(coeff=1.0)`, `similarity_reduce_sympy("wall_jet", printed=True)` |
| R4 wall-jet separation of variables | $\int\frac{df}{f_\infty^{3/2}f-f^2}=\frac16\int d\eta$ | $\int\frac{df}{f_\infty^{3/2}f^{1/2}-f^2}=\frac16\int d\eta$ (and the tail constant 4.29 the book omits) | C13 (N121, D21 step 8) | residual of the printed integral in the check cell |
| R5 (9.33) | drag stated for "the plate" | one side of the plate only | C04 (N40) | `BL.blasius_drag(…, sides=2)` doubles |
| R6 (9.76) | h₉₉ = 5.6152[…]^{1/3} (argument 2.2924) | 2.2924 = arccosh 5 is the 4 % point; 1 % gives arccosh 10 = 2.9932, h₉₉ = 7.3319[…]^{1/3} | C12 (N112, D18) | `JET.free_jet_halfwidth(…, printed=True)` |
| R7 (9.56) | ∂τ/∂y without 1/ρ | kinematic stress ν∂u/∂y, or τ/ρ | C12 (N91, D15 step 3) | unit check in the derivation |
| R8 (9.85) | Ψ written as if a force per length | Ψ has units m⁴/s³ | C13 (N123) | ⚠️ confusion callout |
| R9 Fig. 9.11 text | transition Re quoted as 10⁶, 5×10⁵, 10⁷ | `Re_tr` is an argument (default 5×10⁵) | C09 (N68, N69) | ⚠️ confusion callout |
| R10 §9.9 Magnus sentence | "Re < Re_cr" twice | second one is Re > Re_cr | C11 (N85) | `BB.magnus_sign` truth table |
| R11 §9.4 | solutions "exist for n < −0.0904 with reverse flow" | none with f′ → 1 below the fold; the reversed profiles are the second branch for −0.0904 < n < 0 | C05 (N47) | `BL.falkner_skan(-0.095).success is False` |
| R12 | "Blasius" | boundary-layer solution here, force theorem (6.60) in Ch. 6 — same person | C03 | ⚠️ confusion callout |
| R13 §9.9 | sphere eddy range borrowed from the cylinder | sphere separation starts near Re ≈ 20 (literature) | C11 (N83) | ⚠️ confusion callout |
| R14 §9.10 | turbulent jets "see Chapter 13" | Chapter 12 (turbulence) | C12 (N113) | text correction |
| R15 | λ, L, l, m, θ₀ | Holstein–Bohlen λ; correlation L(λ), l(λ); Thwaites' m = −λ; initial momentum thickness θ₀ | C07 | ⚠️ confusion callout |
| R16 | exercise answers are computable | print at most 4 significant figures, never the private book value | S01, all blocks | rule for the writers |

6. **Colours (text, figures, derivation terms, explainers — one meaning each; curation §5):** outer / ideal flow `teal` ·
   boundary-layer profile u/U `blue` · displacement δ* `amber`, momentum θ `orange`, δ₉₉ `muted` dashed · viscous
   (wall-shear, ν u_yy) `rose` · pressure and pressure-gradient terms `orange` · separation, reverse flow and warnings
   `amber`/red dot · similarity / rescaled curves `accent` purple dashed · exact or analytic ghost muted dashed · numerical
   dots teal · printed-slip ghosts muted dotted, labelled "as printed" · wake and jet `blue`, entrained fluid `teal`.
7. **Every book equation is shown in full next to its number** (CLAUDE.md rule 3) — in markdown, derivation steps
   ("substitute (9.19), $\psi=U\delta f(\eta)$"), traps, recaps, notes and every explainer text. Builders reuse
   `show_eqs(text, EQ)` from `notebooks/build_ch07.py` with an `EQ` dict for all 85 labels of ch09 plus the earlier labels
   cited here ((4.10), (4.19), (4.102), (6.2), (6.60), (6.91), (8.1), (8.14), (8.15), (8.32), (8.42)); JS explainers use a
   local `showEqs`. `tools/eq_refs.py` must list 0 offenders. Exercises are cited as "Exercise 9.6", never as bare equations.
8. **nbkit behaviour** (ch04–ch08 lesson): `nb.recap(...)` and `nb.section(...)` end the current CORE block, so every RECAP
   sits **before** the `nb.core` call of the block that uses it (R01–R03 before C01; R06, R07 before C03; R08, R09 before
   C13) except R04, R05, which follow the completed C01 block as "the boundary conditions every solution below uses".
   Every `nb.derivation` sits inside its curation CORE block with `ref=` a bare equation label ("9.7").
9. **Parsing.** Part E is the only part with table rows after its heading; Part F has no line starting with a table bar
   (absolute values are written \lvert … \rvert or in words). Explainer headings are exactly `### E1 · bl_scaling_thicknesses`
   … `### E9 · teacup_secondary_flow` and `### B1 · ball_swing_magnus`. Primer terms in Part A are the exact Concept text of
   the Part E rows marked "primer" (coverage_check matches the first 18 characters). "Explained by" never names an
   earlier chapter's C-number or an analysis slip label (it cites "Ch. 4 §4.9" or a primers.md entry).
10. **Public repo (rule 9).** No exercise text; Table 9.1, the regime thresholds and angles (82°, 125°, Re ≈ 40, 200,
    3×10⁵, 5×10⁵), the sports-ball figures, the printed jet constants and every exercise answer stay in
    `tests/book_values_ch09.json`. The notebook shows **our** numbers from our inputs (air plate 1 m at 1 m/s; water; a
    slot jet with J = 1 N/m; a 20 m cricket flight with F/W = 0.2; a wire d = 2 mm at 10 m/s). Regime thresholds and
    angles appear only as "the book's rounded values, experimental" in a labelled regime table, never as tested results.
    Numbers that coincide with a printed value (0.332, 0.664, 1.72, 1.328, 4.91) are built from the function and printed
    at 4 significant figures. Figures are drawn by our code (every "Fig. 9.n" is remade, none copied).
11. **Verifier notes found while designing (must reach the verifier).** (i) **Thwaites' accuracy is not ±3 %/±10 % on the
    Falkner–Skan family with the fit L = 0.45 − 6λ**: the notebook prints the measured θ errors (n = 0: +1.0 %, n = 1/3:
    −4.2 %, n = 1: −6.3 %, n = 4: −7.6 %, n = −0.05: +3.1 %) instead of quoting the book's percentages, and `test_ch09`
    must assert |error| ≤ 8 % on n ∈ [−0.09, 4], not 3 %. (ii) **Wall-jet constants:** because f′(0) = 0 the wall-jet
    similarity has a free scale (f → λf(λη) with C → C/λ²): only C f_∞² is physical, so Ψ and ṁ(x) are the same datum
    (Ψ = ṁ⁴/(40 ρ⁴ ν x)); see `JET.wall_jet_constants` below. (iii) **Kármán street:** the linearised growth rate has a
    closed form at the most dangerous wavenumber k = π/a, $\sigma=\frac{\pi\Gamma}{2a^2}\lvert\tfrac12-\mathrm{sech}^2\frac{\pi b}a\rvert$
    (numerics over k ∈ (0, π/a] confirm k = π/a is the maximum for every b/a tried), zero exactly at cosh(πb/a) = √2 — a
    V1 + V2 pair for `karman_street_growth`. (iv) The crude default of `separated_pressure_drag` (base pressure = ideal
    value at the separation point) gives C_D,p = 2.59 at 82° — the notebook always passes `cp_base` explicitly (−1.2
    subcritical, −0.6 supercritical, both *illustrative, ours*). (v) The exact-FS separation member has λ = −0.0681 with
    H ≈ 4.0 (n = −0.0904: l = 0.0041, H = 3.970); the analysis note "H = 3.81" belongs to n = −0.09.

Order of parts: C (contract) · A (notebook storyboard) · B (explainers) · D (runtime) · E (prerequisite ledger) · F
(derivations).

---
## Part C — functions the builders will call (the implementer's contract)

**Status column.** **A#** = planned in `analysis/ch09.md` §4 row # (signature kept or refined here). **§8/§9** = added by the
curation's §8/§9. **NEW** = added by this design (in neither) — flagged as the phase asks. **CHG** = analysis signature
changed here (reason given). **reuse** = exists (ch01–ch08) and is only called. Docstrings cite § and Eq. (from the rendered
page), list symbols with units and assumptions, and carry the validation label. Return types: `float`/`ndarray` for single
quantities, `dict` (keys listed) for several named results. Every function is scalar-callable unless "arrays" is stated.

### C.0 Reused (existing; called, not changed)
| # | Callable (signature) | Returns | Used by | Status |
|---|---|---|---|---|
| 0.1 | `style.setup_notebook() -> bool` (FAST) · `style.savefig(fig, "ch09", name)` | FAST · path | setup, every figure | reuse |
| 0.2 | `anim.animate(update, frames, fig, interval)` · `anim.show_animation(anim, player="video" or "frames")` | HTML | A1–A6 | reuse |
| 0.3 | `interact.slider_figure(fn, name, values, *, unit, xlabel, ylabel, title, xrange, yrange, height)` · `interact.animate_figure(frame_fn, times, …)` · `interact.live(fn, **widgets)` | plotly figure / widget | F1–F7, live cells | reuse |
| 0.4 | `embed.show_viz("ch09", slug)` | display | 9 explainer cells | reuse |
| 0.5 | `tools.convergence.observed_order(h, err)` | log–log slope | C01 (marching), C04 (truncation) | reuse |
| 0.6 | `core.laminar.similarity_variable(y, t, nu, half=False)` · `diffusion_thickness(t, nu, level=0.01)` · `stokes_first_problem(y, t, U, nu)` · `temporal_bl_wall_stress(t, U, nu, rho) -> dict(Cf_coefficient=1.1284, …)` | η · δ · u · dict | C01 (√(νt) bridge), C03, C04 (N41) | reuse |
| 0.7 | `core.similarity.reynolds_number(U, l, nu=None, rho=None, mu=None)` · `strouhal_number(Omega, l, U)` · `sphere_drag_coefficient(Re, model="morrison")` · `pressure_coefficient(p, p_inf, rho, U)` | [–] | C01, C10, C11 | reuse |
| 0.8 | `core.creeping.stokes_drag_coefficient(Re)` · `oseen_drag_coefficient(Re)` (Re on the diameter) | C_D | C11 (Fig. 9.22 asymptotes) | reuse |
| 0.9 | `core.potential.cylinder(U, a)` (ideal surface speed 2U sin φ) · `ch06.cylinder_surface_cp(theta, U, a)` (θ from +x!) · `core.potential.pressure_coefficient(speed, U)` · `ch06.sphere_surface_cp(theta)` | U_e, C_p | C07 (D11), C09 | reuse |
| 0.10 | `core.biot_savart.point_vortex_velocity(x, xv, Gamma, eps=0)` | (u, v) | C10 (from-scratch cross-check of one vortex) | reuse |
| 0.11 | `ch04.wake_drag_per_span(U_of_y, U_inf, rho, H)` | drag per span | C02 (θ-drag parity) | reuse |
| 0.12 | `core.diffusion.crank_nicolson_1d(u0, y, dt, nsteps, D, bc_left, …)` | array | C01 (parabolic contrast) | reuse |
| 0.13 | `core.navier_stokes.navier_stokes_residual`, `ns_incompressible_terms` | residual dict | C01 (full-NS residual of the Blasius field) | reuse |
| 0.14 | `core.thermo.G_BOOK`, `core.rotating` (Ekman preview constants) | constants | C11, C14 | reuse |
| 0.15 | `ch08.similarity_reduce_sympy(case)`, `similarity_ode_solve`, `similarity_collapse_error` → **moved to `core/similarity_reduce.py`** (ch08 re-exports; ch08 tests unchanged) and given the cases `"blasius"`, `"falkner_skan"`, `"free_jet"`, `"wall_jet"` | dict | C03, C05, C12, C13 | A#16 (move flagged by the analyst) |

### C.1 `fluidpy/core/boundary_layer.py` (NEW module, `BL`; re-exported by `ch09`) — attached laminar layers for Ch. 9, 10, 11, 12, 14, 15
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 1.1 | `boundary_layer_scales(U, L, nu, rho=1.2)` | (9.2)–(9.4), v ~ U Re^{−1/2}, τ₀ ~ μU/δ̄, C_f ~ 2/√Re | dict(Re, delta_over_L, delta [m], v_scale [m/s], tau0_scale [Pa], cf_estimate, adv [m/s²] = U²/L, visc [m/s²] = νU/δ̄², visc_x = νU/L²) | C01, E1 | A#1, **CHG** (`rho` added for τ₀ in Pa) |
| 1.2 | `to_bl_variables(x, y, u, v, p, L, U, rho, Re=None)` · `from_bl_variables(xs, ys, us, vs, ps, L, U, rho, Re=None)` | (9.6) | (x*, y*, u*, v*, p*) / inverse | C01 | A#3 |
| 1.3 | `bl_x_momentum_residual(x, y, u, v, dpdx, nu, rho=1.2)` (arrays, shape (ny, nx)) | (9.9) with 2nd-order central differences (`np.gradient(edge_order=2)`) | ndarray residual [m/s²] | C01, C03, C05 | A#4 |
| 1.4 | `bl_pressure_variation(Re, n_x=41, n_y=121)` | (9.8) evaluated on a Blasius field | dict(Re, ratio = max\|∂p*/∂y*\| / max\|∂p*/∂x*\|, order) — ratio → 0 as Re grows | C01 (N10) | A#5 |
| 1.5 | `outer_flow(kind, **params)` — kinds `"flat"` (U), `"wedge"` (n, a), `"stagnation"` (a), `"diffuser"` (U1, L: U_e = U1/(1 + x/L)), `"cylinder"` (U, a: U_e = 2U sin(x/a)), `"retarded"` (U0, L, c: U_e = U0(1 − c x/L), Howarth), `"custom"` (Ue, dUe=None) | (9.11), (9.35) | object with `Ue(x)`, `dUe(x)`, `dpdx(x, rho=1.0)` = −ρ U_e U_e′, `x_stag` | C01, C05, C07, C08, E4 | A#6, **CHG** (`"retarded"`, `"diffuser"` params fixed) |
| 1.6 | `displacement_thickness(y, u, Ue)` · `momentum_thickness(y, u, Ue)` · `delta_level(y, u, Ue, level=0.99)` · `thicknesses(y, u, Ue, level=0.99)` | (9.16), (9.17), δ₉₉ | float [m] / dict(delta99, delta_star, theta, H) | C02, C04, C06, E1 | A#7 (composite `scipy.integrate.simpson` + exponential tail when u(y_max) < U_e; `delta_level` by PCHIP + `brentq`) |
| 1.7 | `profile_shape(name, y, delta, p=None)` — names `"linear"`, `"sine"`, `"cubic"`, `"exponential"` (δ = a), `"power"` (u/U = (y/δ)^{1/p}, p ≥ 1; δ* = δ/(p+1), θ = pδ/((p+1)(p+2)), H = (p+2)/p, δ₉₉ = 0.99^p δ), `"blasius"` (δ = √(νx/U)) | closed forms | dict(u_over_Ue, delta_star, theta, H, delta99) exact; *expect* linear (0.5, 0.1667, 3, 0.99), sine (0.3634, 0.1366, 2.660, 0.910), cubic (0.375, 0.1393, 2.692, 0.918), exponential (1, 0.5, 2, 4.605) — in units of δ (a) | C02, E1, F1 | §8.2 |
| 1.8 | `falkner_skan(m, eta_max=10.0, method="bvp", tol=1e-10, branch="attached", n=400)` — methods `"bvp"`, `"toepfer"` (m = 0 only), `"shoot"` (m ≥ −0.05) | (9.27), (9.36), (9.28), (9.29) | dict(eta, f, fp, fpp, fppp, fpp0, success); `success is False` below the fold (R11); `branch="reversed"` returns the second branch for −0.0904 < m < 0 | C04, C05, E2, E3 | A#8 |
| 1.9 | `blasius_constants()` | (9.30)–(9.33) | dict(fpp0 = 0.332057, eta99 = 4.90999, delta_star = 1.720788, theta = 0.664115, H = 2.591100, v_inf = 0.860394, tau_coeff = fpp0, cf_coeff = 2 fpp0, cd_coeff = 4 fpp0) computed, cached | C04, E1, E2 | A#9 |
| 1.10 | `blasius_fields(x, y, U, nu)` (arrays, broadcast) | (9.19), (9.23), (9.24), (9.26) | dict(u, v, psi, eta, delta, tau0_over_rho) (f from `solve_ivp` dense output, 1e-10) | C01–C04, E1, E2 | A#10 |
| 1.11 | `blasius_delta99(x, U, nu, printed=False)` · `blasius_delta_star(x, U, nu)` · `blasius_theta(x, U, nu)` · `blasius_wall_shear(x, U, rho, nu)` · `blasius_skin_friction(Re_x)` · `blasius_drag(L, U, rho, nu, sides=1)` · `blasius_drag_coefficient(Re_L, sides=1)` | (9.30)–(9.33); R2 (`printed=True` → 4.93), R5 (`sides=2` doubles) | floats | C04, C06, E1, E2 | A#11, **CHG** (`blasius_delta_star`, `blasius_theta`, `sides` added) |
| 1.12 | `blasius_far_field(eta)` | f′ − 1 ≈ −A√π erfc((η − δ*)/2) (linearised about f ≈ η − δ*; A = f″(5)e^{(5−δ*)²/4} = 0.2339 taken once from the solved profile; within 0.3 % of the solved f′ for 4 ≤ η ≤ 8; its large-η form −(2A/s)e^{−s²/4}, s = η − δ*, is the book's (1/η)e^{−η²/4} with the shift restored) | ndarray | C04 (N33) | A#12 |
| 1.13 | `falkner_skan_state(m)` | (9.36), integrals I_δ = ∫(1 − f′), I_θ = ∫f′(1 − f′) | dict(fpp0, fppp0 = −m, I_delta, I_theta, H, lam = m I_θ², l = I_θ f″(0), cf_sqrtRex = 2 f″(0), inflection_eta or None, separated) — *expect* m = 1: (1.2326, 2.2162, 0.08546, 0.36034); m = 1/3: (0.75745, 2.2969, 0.06134, 0.32494); m = 0: (0.33206, 2.5911, 0, 0.22052); m = −0.05: (0.21348, 2.8182, −0.028235, 0.16042, η_infl = 1.650); m = 4: (2.4057, 2.1725, 0.10033, 0.38101); m = −0.0904: (0.00477, 3.970, −0.06811, 0.00414) | C05, C07, C08, E3, E4 | A#13 |
| 1.14 | `falkner_skan_separation()` | fold of the attached branch (f″(0; m) = 0) | dict(m_sep = −0.09043, beta_sep = −0.19884, fpp0_at_sep = 0) (parametrise by f″(0) to cross the fold; ±1e-5) | C05, C08, E3 | A#14 |
| 1.15 | `falkner_skan_fields(x, y, m, a, nu)` · `falkner_skan_thickness(x, m, a, nu)` | (9.34), (9.35), δ = √(νx^{1−n}/a) | dict(u, v, psi, eta, delta, Ue) / float | C05, E3 | A#15, A#56 |
| 1.16 | `falkner_skan_table(m_grid, n_eta=161, eta_max=8.0)` | tabulated (9.36) solutions for JS explainers | dict(m, eta, fp[m, η], fpp[m, η], fpp0, I_delta, I_theta, H, lam, l, inflection_eta) at 6 s.f.; m_grid 41 nodes dense near the fold (default `np.r_[−0.09043 + 1e-4 … 4]`); also `fpp0_reversed` (second branch, nan where absent, for −0.0904 < m < 0) | E3 (table), E4 (closure) | **NEW** (curation §8.1) |
| 1.17 | `momentum_integral_residual(x, Ue, theta, delta_star, tau0, rho=1.0)` | (9.43) | ndarray residual (4th-order differences) — ≈ 0 for Blasius and every FS member | C06 (from scratch), E4 | A#17 |
| 1.18 | `karman_pohlhausen(Ue, x, nu, profile="cubic", theta0=0.0)` — profiles `"cubic"`, `"sine"`, `"quartic"` | (9.43) with an assumed shape | dict(delta, theta, delta_star, tau0) — *expect* flat plate cubic: δ = 4.641√(νx/U), θ = 0.6465 (−2.7 %), δ* = 1.740 (+1.1 %), τ₀ coefficient 0.3232 (−2.7 %); sine: δ = 4.795, θ = 0.6551, τ₀ 0.3276 (−1.3 %) | C06 (N55) | A#19 |
| 1.19 | `holstein_bohlen(theta, dUe_dx, nu)` · `thwaites_l(lam, closure="falkner_skan")` · `thwaites_H(lam, closure)` · `thwaites_L(lam, closure)` · `thwaites_closure_table(n_points=60)` — closures `"falkner_skan"` (exact FS, PCHIP, λ ∈ [−0.0681, 0.1065]) and `"white"` (l = (λ + 0.09)^0.62, zero at −0.09, the book's criterion); the FS closure returns l = 0 for λ < −0.0681 and continues above 0.1065 with the white fit scaled to match | (9.44)–(9.46), (9.48) | float / dict(m, lam, l, H, L) — *expect* L(0) = 0.4410, L(0.0855) = 0.0001, L(−0.0675) = 0.818 | C07, E4, F4 | A#20 |
| 1.20 | `thwaites(x, Ue, nu, theta0=0.0, closure="falkner_skan", rho=1.0, lam_sep=None)` (`Ue` callable or array) | (9.50), (9.45), (9.46) | dict(theta, delta_star, tau0, cf, lam, H, l, separated, x_sep) — first crossing of `lam_sep` (interpolated); stagnation-point limit handled analytically (θ² ∝ x); stops the marching after separation | C07, C08, E4 | A#21, **CHG** (`lam_sep` keyword; `None` → −0.0681 for the exact-FS closure, −0.09 for `closure="white"`, the book's criterion; the FS closure returns l = 0 below λ = −0.0681, so τ₀ and the flag agree) |
| 1.20b | `thwaites_named(kind, x, nu, theta0=0.0, closure="falkner_skan", **params)` — kinds `"flat"` (U), `"diffuser"` (U1, L), `"wedge"` (n, a), `"cylinder"` (U, a), `"retarded"` (U0, L) — one-station scalars for explainer parity rows | dict(theta, delta_star, tau0, cf, lam, H, l, x_sep) at station x (θ from the marched integral of U_e⁵ with 4000 graded points) | E4 | **NEW** |
| 1.21 | `thwaites_cylinder_closed_form(phi)` = 0.45 cos φ F(φ)/sin⁶φ with F(φ) = 8/15 − cos φ + (2/3)cos³φ − (1/5)cos⁵φ · `thwaites_cylinder(phi_deg, closure=None)` (numeric, parity) · `thwaites_cylinder_separation(lam_sep=-0.09)` → φ_sep in degrees | D11 | float — *expect* λ(82°) = 0.02630, λ(90°) = 0, φ_sep = 103.11° (λ = −0.09), 100.89° (λ = −0.0681) | C07, E4 | §8.4 (**NEW** `thwaites_cylinder_separation`) |
| 1.22 | `wall_curvature(dpdx, mu)` · `profile_inflection(y, u)` (→ y_i or None) · `separation_point(x, tau0)` (first sign change, linear interpolation → x_sep or None) | μ u_yy(wall) = dp/dx; (9.51), (9.52); τ₀ = 0 | float / None | C08, E3, E4 | A#24 |
| 1.23 | `march_boundary_layer(Ue, x_grid, nu, u_inlet, ny=400, y_max_factor=12)` | (9.9) in von Mises variables (x, ψ), backward Euler in x, central in ψ; stops at τ₀ ≤ 0 | dict(x, y, u, v, tau0, x_sep or None) — τ₀ within 1e-3 of Blasius (U_e = const); forgets two different inlet profiles when U_e′ ≥ 0 | C01 (N14), C03 (N24), E1 | A#25 (**required**, not optional: C01 uses it) |
| 1.24 | `transition_state(Re_x, Re_cr=5e5)` → "laminar" / "transitional" / "turbulent" (transitional between Re_cr and 10 Re_cr) · `plate_drag_coefficient(Re_L, regime="laminar" or "turbulent" or "mixed", Re_tr=5e5, sides=1)` | laminar 1.328/√Re_L; turbulent 0.074 Re_L^{−1/5} (Prandtl one-seventh-law form; **cite Schlichting, *Boundary-Layer Theory*, and White, *Viscous Fluid Flow*, in the docstring before use**); mixed = turbulent − A/Re_L with A = Re_tr(C_turb − C_lam) ≈ 1743 at 5×10⁵ | float | C09 (N68, N69), F5 | A#26 |

### C.2 `fluidpy/core/jets.py` (NEW module, `JET`; re-exported by `ch09`) — free and wall jets for Ch. 9, 11, 12, 13
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 2.1 | `free_jet_constants()` | (9.72), (9.73), (9.76) | dict(C = 4√6/3 = 3.26599, mdot_coeff = 36^{1/3} = 3.30193, h99_arg = arccosh 10 = 2.99322, h99_coeff = √6 arccosh 10 = 7.33187, h99_printed = 5.61529, u0_coeff = (1/C²)^{1/3} = 0.45428, f_inf = √6) computed via sympy/`quad`, not typed | C12, E7 | A#32 |
| 2.2 | `free_jet(x, y, J, rho, nu)` (arrays) | (9.62)–(9.64), (9.71), (9.74) | dict(u, v, psi, eta, u0, delta) — *expect* air J = 1 N/m, x = 0.1 m: u₀ = 35.14 m/s, δ = 0.2066 mm | C12, E7 | A#33 |
| 2.3 | `free_jet_ode_solve(eta_max=12.0, n=400, coeff=3.0)` | 3f‴ + ff″ + f′² = 0 with f(0) = 0, f′(0) = 1, f′(η_max) = 0 (`solve_bvp`) | dict(eta, f, fp, max_err vs √6 tanh(η/√6) ≤ 1e-8) | C12 (from scratch), E7 | A#34 |
| 2.4 | `free_jet_centreline(x, J, rho, nu)` · `free_jet_thickness(x, J, rho, nu)` · `free_jet_mass_flux(x, J, rho, nu)` · `free_jet_halfwidth(x, J, rho, nu, level=0.01, printed=False)` · `free_jet_reynolds(x, J, rho, nu)` → dict(Re_x, Re_h99, Re_h99_printed) · `free_jet_entrainment_velocity(Re_x)` → v/u₀ limit √6/(3√Re_x) · `jet_momentum_flux(y, u, rho)` | (9.62), (9.63), (9.73), (9.76) corrected (R6; `printed=True` → 5.6152), (9.75), (9.58) | floats — *expect* air J = 1, x = 0.1: ṁ = 0.04268 kg/(m s), h₉₉ = 1.515 mm, Re_x = 2.343×10⁵, edge v = −0.0593 m/s; water J = 1: u₀ = 0.9787 m/s, δ = 0.3196 mm, ṁ = 1.533 kg/(m s), h₉₉ = 2.344 mm | C12, E7 | A#35 |
| 2.5 | `free_jet_profile(eta)` → dict(f = √6 tanh(η/√6), fp = sech²(η/√6), fpp) · `free_jet_at_level(level)` → dict(z = arccosh(1/√level), coeff = √6 z) (level 0.01 → 7.3319; 0.04 → 5.6153; 0.5 → 2.1589) | (9.70), (9.71), (9.76) | dict | E7 (level control) | **NEW** |
| 2.6 | `wall_jet_ode_solve(fpp0=1.0, eta_max=40.0, coeff=4.0)` (`coeff=1.0` is the printed wrong variant, R3) | 4f‴ + ff″ + 2f′² = 0, IVP with f(0) = f′(0) = 0 | dict(eta, f, fp, f_inf, err_vs_9_83, fpp0_over_finf_cubed = 1/72) — *expect* f″(0) = 1/72 ⇒ f_∞ = 1.0000000, f′ peak 0.078745 at η = 8.11 | C13 (from scratch), E8 | A#38 |
| 2.7 | `wall_jet_profile(eta, f_inf=1.0)` | (9.83) inverted by `brentq` on g ∈ [0, 1): f = f_∞ g², f′ = f_∞² g (1 − g³)/6 | dict(g, f, fp) | C13, E8 | A#36 |
| 2.8 | `wall_jet(x, y, C, f_inf, nu, rho=1.0)` (arrays) | (9.82), (9.84), v = −√(νC)(f − 3ηf′)/(4x^{3/4}) | dict(u, v, psi, delta, u0 = C x^{−1/2}) | C13, E8 | A#37 |
| 2.9 | `wall_jet_invariant(y, u)` | (9.80): ∫₀^∞ u(∫_y^∞ u² dy′) dy by reverse `cumulative_trapezoid` | float | C13, E8 | A#39 |
| 2.10 | `wall_jet_integrals(f_inf=1.0)` | scaling f → λf(λη): ∫f′dη = f_∞, ∫f′²dη = f_∞³/18, ∫f′∫f′² dη dη = f_∞⁴/40, f″(0) = f_∞³/72 | dict(int_fp, int_fp2, invariant, fpp0) — *expect* f_∞ = 1: (1, 0.055556, 0.025, 0.013889) | C13, E8 | **NEW** |
| 2.11 | `wall_jet_constants(rho, nu, *, Psi=None, mdot=None, x=None, f_inf=1.0)` | (9.85), (9.84): **one physical datum** (R8/§11 note (ii)): given `Psi` or (`mdot`, `x`) and the gauge `f_inf` returns C; the other is a consequence: Ψ = ṁ⁴/(40 ρ⁴ ν x); giving both returns their consistency residual | dict(C, C_f_inf_sq = C f_∞² (gauge-free), Psi, mdot_at_x, residual) | C13 (N123), E8 | A#39, **CHG** (was `(Psi, mdot, x, rho, nu) -> C, f_inf`; that pair is degenerate) |
| 2.12 | `wall_jet_mass_flux(x, C, f_inf, rho, nu)` | (9.84) ρ√(νC) f_∞ x^{1/4} | float | C13 | A#39 |

### C.3 `fluidpy/core/bluff_body.py` (NEW module, `BB`; re-exported by `ch09`) — cylinder, sphere and street physics for Ch. 9, 11, 12, 14
| # | Callable (signature) | Implements | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 3.1 | `cylinder_flow_regime(Re)` · `sphere_flow_regime(Re)` | regime tables (the book's rounded thresholds, flagged "book's rounded values, experimental" in the docstring). **Exact ASCII labels the explainers copy** — cylinder: Re < 1 "creeping flow: symmetric, no wake"; 1 ≤ Re < 4 "creeping to attached eddies"; 4 ≤ Re < 40 "two steady attached eddies"; 40 ≤ Re < 200 "laminar Karman street"; 200 ≤ Re < 3000 "irregular vortices, St stays near 0.2"; 3000 ≤ Re < 3e5 "subcritical: laminar separation, wide wake"; 3e5 ≤ Re < 6e5 "critical: drag crisis"; Re ≥ 6e5 "supercritical: turbulent separation, narrow wake". Sphere: Re < 1 "creeping flow: symmetric, no wake"; 1 ≤ Re < 130 "steady wake with an attached doughnut eddy"; 130 ≤ Re < 3e5 "unsteady wake, loops shed (no regular street)"; 3e5 ≤ Re < 8e5 "critical: drag crisis"; Re ≥ 8e5 "supercritical: turbulent separation, narrow wake" | dict(label, thresholds, separation_deg (82.0, 125.0 or None), St (0.2 for 40 ≤ Re < 3000 else None)) | C10, C11, E5 | A#27 |
| 3.2 | `cylinder_state(Re, rough=False)` · `drag_crisis_state(Re, rough=False, Re_cr=3e5)` · `cylinder_cd_schematic(Re)` (labelled *qualitative*, no dataset) | regime → φ_sep, St, model C_D; roughness shifts Re_cr down | dict(label, phi_sep_deg, St, cb, cd_model, qualitative=True): `phi_sep_deg`, `cb`, `cd_model` are `None` for Re < 3000 (the separated model applies to high-Re separated flow only); otherwise φ_s = 82° (subcritical), 125° (supercritical, or rough with Re ≥ Re_cr/3), blended across the critical band, and C_b = −1.2, −0.6 (illustrative, ours) | C11, E5 | **NEW** (curation §8.5) |
| 3.3 | `separated_cp(phi_deg, phi_sep_deg, cp_base=None)` · `separated_pressure_drag(phi_sep_deg, cp_base=None)` · `pressure_drag_from_cp(phi_deg, cp)` (Gauss–Legendre / periodic trapezoid of ½∮C_p cos φ dφ) | D13: C_p = 1 − 4 sin²φ for φ ≤ φ_s, C_b behind; C_D,p = sin φ_s (1 − (4/3) sin²φ_s − C_b); `cp_base=None` → the crude rule C_b = ideal value at φ_s | float — *expect* (82°, −1.2) = 0.8840, (125°, −0.6) = 0.5778, (90°, ideal C_b = −3) = 2.667, φ_s → 180° with C_b → 1: → 0 | C09, C11, E5 | A#30, **NEW** (`pressure_drag_from_cp`) |
| 3.4 | `karman_street_ratio()` = arccosh(√2)/π = 0.280550 | cosh(πb/a) = √2 | float | C10, E6 | A#28 |
| 3.5 | `karman_street_spectrum(b_over_a, offset=0.5, k=np.pi, Gamma=1.0, a=1.0, n_terms=2000)` · `karman_street_growth(b_over_a, offset=0.5, Gamma=1.0, a=1.0, k=None, n_terms=2000)` (`k=None`: maximum Re λ over 61 wavenumbers in (0, π/a]) · `karman_street_growth_closed(b_over_a, Gamma=1.0, a=1.0)` = (πΓ/2a²)\|½ − sech²(πb/a)\| | linearised point-vortex street (D14); lattice sums exact for k = 0 and k = π/a, else truncated at ±n_terms with the k = 0 tail 2/N added | complex ndarray(4) [1/s], sorted by (real part, imaginary part) / float [1/s] — *expect* (Γ = a = 1): growth 0.6400 (b/a = 0.1), 0.2982 (0.2), 0.1099 (0.25), 0 (0.28055, < 1e-8), 0.0663 (0.3), 0.2208 (0.35), 0.5359 (0.5), 0.7737 (1.0); non-staggered (`offset=0`) π/4 = 0.7854 for every b/a; at b/a = 0.28055 the spectrum is ±0.7854i | C10 (from scratch), E6 | A#28, **CHG** (was `(b_over_a, n_pairs, Gamma)`: an infinite periodic row is exact) |
| 3.6 | `karman_street_velocity(a, b, Gamma)` = (Γ/2a) tanh(πb/a) · `karman_street_positions(t, b_over_a, eps, mode="stable" or "unstable", n_cells=8)` (linear evolution of a perturbed street for A5, E6) | D14 step 4 | float [m/s] / dict(zA, zB) | C10, E6 | A#28, §8.6 |
| 3.7 | `shedding_frequency(U, d, St=0.2)` → dict(f [Hz] = St U/d, omega_rad [rad/s] = 2πf; the book's Ω in St = Ωd/U ≈ 0.2 is the frequency f in Hz) · `strouhal_of_re(Re)` (Roshko plateau 0.21 above Re ≈ 300; below only with a cited fit — else `nan` and a warning) | (4.102) recalled | dict / float — *expect* U = 10 m/s, d = 2 mm: f = 1000 Hz, Ω = 6283 rad/s | C10, E5 | A#29 |
| 3.8 | `ball_swing_deflection(F_over_W, distance, U, g=9.81)` = ½(F/W)g(d/U)² · `magnus_sign(Re_slow, Re_fast, Re_cr)` → "+" (both sides on the same side of Re_cr), "−" (Re_slow < Re_cr ≤ Re_fast), "none" (equal) | N84, N85 | float [m] / str — *expect* (0.2, 18 m, 35 m/s) = 0.2595 m | C11, B1 | A#31 |

### C.4 `fluidpy/ch09_boundary_layers.py` (chapter module; re-exports C.1–C.3, `core.similarity_reduce`, `sphere_drag_coefficient`, `strouhal_number`)
| # | Callable (signature) | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 4.1 | `bl_nondim_sympy(printed_9_7=False)` | D01: substitute (9.6) into (9.1), the y-momentum equation and (6.2), divide by the largest term | dict(x_momentum, y_momentum, continuity, coefficients (powers of Re: x-mom d²u/dx*² → −1, d²u/dy*² → 0; y-mom inertia → −1, pressure 0, d²v/dx*² → −2, d²v/dy*² → −1; continuity 0), limit); `printed_9_7=True` leaves un-squared denominators and fails the dimension check | C01 | A#2 |
| 4.2 | `momentum_integral_sympy()` | D08 nine displayed lines with a generic u = U_e F(y/δ(x)); Blasius residual | dict of steps | C06 | A#18 |
| 4.3 | `thwaites_sympy()` | D09 (9.47) → (9.50), `dsolve` of (9.49) | dict of steps | C07 | A#22 |
| 4.4 | `jet_momentum_sympy()` · `wall_jet_invariant_sympy()` · `wall_jet_sympy()` | D15; D19 (9.79)–(9.81); D21 first integrals with `apart` | dicts of steps | C12, C13 | A#40 |
| 4.5 | `karman_street_sympy()` | D14: builds the k = π/a matrix from the lattice sums, factorises the characteristic polynomial ((μ − γ)² + σ²)((μ + γ)² + σ²), solves γ = 0 | dict(matrix, char_poly, factors, growth, b_over_a_marginal) | C10 (check cell) | **NEW** |
| 4.6 | `cylinder_thwaites_sympy()` | D11: F(φ) by substitution c = cos φ, λ(φ), the limit φ → 0 (0.45/6) | dict | C07 | **NEW** |
| 4.7 | `example_9_1()` · `example_9_2(theta0=0.0)` | Thwaites for Blasius (θ = 0.6708 √(νx/U), +1.0 %; δ* = 1.738; C_f√Re_x = 0.6575, −1.0 %) · diffuser λ(x/L) = −(0.45/4)[(1 + x/L)⁴ − 1] (θ₀ = 0), x_sep/L = 1.8^{1/4} − 1 = 0.15829 (0.12563 with λ_sep = −0.0681) | `example_9_1` → dict(theta_coef = 0.6708, theta_err = +0.010, delta_star_coef = 1.738, cf_sqrtRex = 0.6575, cf_err = −0.010); `example_9_2` → dict(lam (callable of x/L), x_sep_over_L = 0.15829 (λ_sep = −0.09), x_sep_over_L_fs = 0.12563 (−0.0681), theta0 used) — 4 s.f. | C07, E4 | A#23 |
| 4.8 | `secondary_flow_radial_force(u_inviscid, u_layer, R, rho)` = ρ(u_e² − u²)/R [N/m³] · `secondary_flow_layer_profile(z, delta, u_e, shape="exponential" or "linear" or "power", n=1/7)` (illustrative, docstring says "illustrative, not a solution") | D22 | float / ndarray | C14, E9 | A#41, §8.8 |
| 4.9 | `book_slips()` | the slips table of convention 5 (id, printed, correct, evaluator) | list of dicts | callouts, tests | A#42, §8.9 |
| 4.10 | `derive_all()` | runs every sympy engine and returns their residual summaries | dict | tests | A#42 |

### C.5 `scripts/ch09_*.py` (runnable demos, each `--no-show` for headless runs; drawing helpers carry no physics)
The 19 scripts of analysis §4 row 43 (`boundary_layer_picture`, `displacement_streamlines`, `blasius`, `blasius_vs_temporal`,
`parabolic_marching`, `falkner_skan`, `thwaites_L`, `thwaites_accuracy`, `pressure_gradient_profiles`, `plate_regimes`,
`plate_drag_curve`, `cylinder_drag_model`, `cylinder_regimes`, `cylinder_cp_cd`, `sphere_drag_curve`, `ball_swing`,
`magnus_sign`, `jets`, `teacup`) plus `ch09_karman_street.py` (A5), `ch09_thickness_shapes.py` (F1), `ch09_free_vs_wall_jet.py`
(A6) and `ch09_drawings.py` (helpers: `plate_layer`, `cylinder_sketch`, `jet_sketch`, `cup_section`, `profile_arrows`). Each
figure function returns the `Figure` and takes `ax=None`.

### C.6 Not in `analysis/ch09.md` §4 (flag list for the implementer)
`BL.boundary_layer_scales(…, rho)` (CHG) · `BL.profile_shape` · `BL.falkner_skan_table` · `BL.thwaites_cylinder_*` ·
`BL.blasius_delta_star`, `blasius_theta`, `sides=` · `BL.thwaites(…, lam_sep=)` · `BL.march_boundary_layer` now **required** ·
`JET.free_jet_profile`, `free_jet_at_level` · `JET.wall_jet_integrals` · `JET.wall_jet_constants` (CHG: degenerate pair) ·
`BB.cylinder_state`, `drag_crisis_state`, `cylinder_cd_schematic`, `pressure_drag_from_cp` · `BB.karman_street_spectrum`,
`karman_street_growth_closed`, `karman_street_positions` (CHG signature of `karman_street_growth`) · `ch09.karman_street_sympy`,
`cylinder_thwaites_sympy` · the ch08 → `core/similarity_reduce.py` move. Nothing else in Parts A, B or F calls an unlisted name.

---
## Part A — notebook storyboard (`notebooks/build_ch09.py` → `notebooks/ch09_boundary_layers.ipynb`)

**One line per book section** (cell numbers are estimates for the builder's budget; ≈ 640 cells in all):
- §9.1 → R01 R02 R03, C01 (N01–N14, N16 · N15 N17 pointers · D01 D02 · P200 P201), R04 R05 — cells ≈ 8–92
- §9.2 → C02 (N18–N21 · D03 D04 · P202 P203 P204 P205 · **A1** · F1 · **E1**) — cells ≈ 93–150
- §9.3 → R06 R07, C03 (N22–N32, N24 · D05 · P206 · R12 callout), C04 (N33–N42 · D06 · P207 P208 · **A2** · F2 · **E2**) — cells ≈ 151–265
- §9.4 → C05 (N43–N47, N48 · D07 · P209 · **A3** · F3 · live · **E3**) — cells ≈ 266–320
- §9.5 → C06 (N49–N55 · D08 · optional Kármán–Pohlhausen) — cells ≈ 321–365
- §9.6 → C07 (N56–N67 · D09 D10 D11 · P210 P211 · **A4** · F4 · **E4**) — cells ≈ 366–440
- §9.7 → C08 (N70–N74 · D12 · P212 · F7), C09 (N68 N69 N75 · D13 · F5 · **E5** shared) — cells ≈ 441–500
- §9.8 → C10 (N76–N80 · D14 · P213 P214 · **A5** · **E6**), C11 (N81 N82 N86 · **E5**) — cells ≈ 501–560
- §9.9 → covered inside C11 (N83–N86 · B1 pointer) — cells ≈ 561–580 (same block, second half)
- §9.10 → C12 (N87–N113 · D15 D16 D17 D18 · P215 P216 P217 · **A6** · F6 · **E7**), R08 R09, C13 (N114–N124 · D19 D20 D21 · P218 P219 · **E8**) — cells ≈ 581–655
- §9.11 → C14 (D22 · P220 · **E9**), S01 pointer, summary — cells ≈ 656–690

Every CORE block below follows the order: problem in plain words → idea → primers → maths (notes and derivations, Part F) →
tiny example → code (fluidpy) + "What does the code above do?" → from-scratch check → visual(s) → notes and "What would
change if…". Code drafts are intent + exact calls; the builder writes every line commented (novice grade, units, the
equation written out next to its number). *expect* gives the numbers the executed cell must print (sanity values computed
for this design; air ν = 1.5 × 10⁻⁵ m²/s, ρ = 1.2 kg/m³, water ν = 10⁻⁶ m²/s, ρ = 1000 kg/m³ unless stated). *see / read /
change* are the three figure notes. Note ids open every note text in bold (**N09 [B]**) so the lesson-reviewer can find them.

### A.0 Front matter
1. `nb.title(big_idea=…, roadmap=[…14…], prerequisites=[…])`. **Big idea (draft):** "Air flowing over a wing, water past a
   bridge pier, wind round a chimney: at high Reynolds number the flow is almost ideal everywhere except in a razor-thin
   layer at the wall, where friction drags the speed from zero (no slip) up to the outer value. Prandtl saw in 1904 that
   this layer can be studied on its own — it is thin, so its equations are simpler than Navier–Stokes (the streamwise
   diffusion and the cross-stream pressure change drop out), and the ideal outer flow simply hands it the pressure. This
   chapter derives those equations, solves them exactly for a flat plate (Blasius) and for wedges (Falkner–Skan), turns them
   into one integral law (von Kármán) and a one-line predictor (Thwaites), and uses them to explain what a thin layer does
   to a body: it thickens, it can let go (separation), and where it lets go decides the drag of cylinders, spheres and
   cricket balls. The same equations describe layers with no wall at all — the jet — and hint at the flow a layer sets up
   sideways in a teacup, the doorway to the Ekman layers of Ch. 13." **Roadmap (one line per CORE):** C01 how thin the
   layer is and which terms survive · C02 three thicknesses δ₉₉, δ*, θ · C03 the flat plate reduces to one ODE · C04 the
   Blasius solution and its numbers · C05 wedge flows (Falkner–Skan) · C06 the von Kármán momentum integral · C07 Thwaites'
   method · C08 wall curvature and separation · C09 streamlining and form drag · C10 the wake of a cylinder and the Kármán
   street · C11 the drag crisis, spheres and sports balls · C12 the free laminar jet · C13 the wall jet · C14 secondary flow
   in a teacup. **Prerequisites:** Navier–Stokes, no slip and Bernoulli (Ch. 4 §4.9–4.10) · Reynolds number and
   non-dimensionalisation (Ch. 4 §4.11) · stream function (Ch. 4 §4.3) · ideal flow past a cylinder and sphere, d'Alembert's
   paradox (Ch. 6 §6.5, §6.8) · point vortices (Ch. 5 §5.7) · the two-length scaling, similarity variable and √(νt) (Ch. 8
   §8.3–8.4).
2. `nb.explainer_index([...])` — 9 rows: ("bl_scaling_thicknesses", "How thin is a boundary layer — and what is its
   thickness?", "δ₉₉, δ* and θ are three integrals of one profile; the layer grows like √(νx/U) = x Re^(−1/2)") ·
   ("blasius_similarity_collapse", "Why does one curve describe the whole plate?", "profiles at every x fall on f′(η) once
   y is scaled by √(νx/U); f″(0) = 0.332 then gives wall shear, drag and momentum loss") · ("falkner_skan_family", "How can
   favourable and adverse gradients be one family?", "one dial n: n > 0 fuller profile, n = 0 Blasius, n < 0 an inflection
   and, at n = −0.0904, zero wall shear") · ("thwaites_marching", "Can one integral predict separation?", "θ² ∝ ∫U_e⁵dx and
   λ = (θ²/ν)dU_e/dx crossing −0.09 tell where the wall shear vanishes") · ("cylinder_drag_crisis", "Why does drag fall
   when the flow gets faster?", "where the layer lets go sets the wake pressure; a turbulent layer holds on longer, 82° →
   125°") · ("karman_street_stability", "Why is the vortex street's shape fixed at b/a ≈ 0.28?", "a staggered double row of
   point vortices grows a wiggle at every spacing except cosh(πb/a) = √2") · ("free_jet_similarity", "A jet keeps its
   momentum, spreads and slows — where does the extra mass come from?", "constant J forces u₀ ∝ x^(−1/3), δ ∝ x^(2/3);
   the jet entrains ambient fluid") · ("wall_jet_invariant", "In a wall jet the wall eats momentum — what is conserved?",
   "the flux of exterior momentum flux; it fixes x u₀² = const and the exponent 3/4") · ("teacup_secondary_flow", "Why do
   tea leaves collect at the centre?", "a slowed layer under an unchanged pressure gradient is pushed inward").
3. `nb.setup()`.
4. `nb.code` — **chapter imports** (outside any block): `import numpy as np` · `import matplotlib.pyplot as plt` ·
   `import sympy as sp` · `from scipy import integrate, optimize, interpolate` · `from fluidpy import ch09_boundary_layers
   as ch09` · `from fluidpy.core import boundary_layer as BL, jets as JET, bluff_body as BB` · `from fluidpy.core import
   laminar as LAM, similarity as SIM, creeping as CRP, potential as PF` · `from fluidpy import ch04_conservation_laws as
   ch04, ch06_ideal_flow as ch06, ch08_laminar_flow as ch08` · `from fluidpy.core.interact import slider_figure,
   animate_figure, live` · `from fluidpy.core.anim import animate` · `from fluidpy.core.style import savefig` · `import
   sys; sys.path.insert(0, "scripts"); from ch09_drawings import plate_layer, cylinder_sketch, jet_sketch, cup_section,
   profile_arrows` · `NU_AIR, RHO_AIR = 1.5e-5, 1.2` · `NU_W, RHO_W = 1.0e-6, 1000.0`. *explain:* one line per import ("`ch09`
   re-exports the three new toolkits `BL`, `JET`, `BB`, so `ch09.blasius_constants` and `BL.blasius_constants` are the same
   function").
5. `nb.md` — **⚠️ Conventions in this chapter** (a two-column table "symbol | meaning here (and before)", then numbers): x
   along the surface, y normal; u = ∂ψ/∂y; dp/dx > 0 adverse; U_e(x) the edge speed, U the plate's constant speed;
   δ̄ (order of magnitude), δ = √(νx/U) (similarity length — *not* δ₉₉), δ₉₉, δ*, θ (a **length** here, an angle in Ch. 3, 6,
   8), h₉₉ (jet half-width); six Reynolds numbers (P200); φ from the **forward** stagnation point (82°, 125°); λ (Holstein–
   Bohlen, dimensionless) versus the wavelength of Ch. 7; L(λ) the correlation versus L the length; "Blasius" the plate
   solution here, the force theorem (6.60) in Ch. 6; C the constant of three different jobs (`C_fj`, `C_wj` in code). "We
   compute with the book's conventions and name the other one wherever it differs." A second table lists the sixteen
   printed slips of convention 5 with where each is taught ("the book prints X, correct is Y").
6. `nb.md` — **🔁 Tools from earlier chapters used in this one** (one line each, primer number): partial derivative (P25),
   Taylor (P26), definite integral (P27), `solve_ivp` (P31/P94), trapezoid rule (P37), product rule (P38), sympy (P40),
   separation of variables (P42), exponent rules (P43), integrating an ODE twice (P44), `np.where`/`np.select` (P46),
   chain rule (P49/P91), orders of smallness (P68), `meshgrid`/contour/streamplot (P76/P78), fundamental theorem of
   calculus (P84), `quad` (P87), substitution in an integral (P106), `brentq` (P108), differentiation under the integral
   (P109), log–log axes (P13), `np.gradient` and finite differences (P21/P22), `np.sign` crossings (P180), `np.interp`
   (P182), diffusion time (P185), two-length scaling (P188), Leibniz with a moving limit (P189), `cumulative_trapezoid` (P190),
   implicit stepping and Crank–Nicolson (P192/P193), `solve_bvp` (P196), exponent matching (P197), dominant balance (P198),
   eigenvalues (P80), complex conjugate (P81), hyperbolic functions (P168), `assert np.allclose` (P15), `animate` (P16),
   `slider_figure` (P17), `show_viz` (P18), live widgets (P47). Each block repeats the ones it uses in a one-line
   reminder at first use.

---

### A.1 §9.1 Introduction — R01 R02 R03, C01 (+N01–N14, N16 · N15 N17 · D01 D02 · P200 P201), R04 R05
1. `nb.section("9.1", "Introduction", intro="**What is this section about?** Ideal flow says a body in a steady stream feels
   no drag (d'Alembert, Ch. 6) — yet every real body does. Prandtl's resolution: viscosity matters only in a very thin layer
   next to the body, where the fluid must stick to the wall; outside it the flow is the ideal flow you already know. This
   section works out how thin the layer is, which terms of Navier–Stokes survive in it, and what the outer flow must tell it.
   Everything in the rest of the chapter solves the equations found here.")`
2. `nb.recap("R01", "Steady two-dimensional momentum along a surface", "Ch. 4 (4.10) and Ch. 8 (8.1) gave, for steady flow in
   the plane, $u\\frac{\\partial u}{\\partial x}+v\\frac{\\partial u}{\\partial y}=-\\frac1\\rho\\frac{\\partial p}{\\partial x}+\\nu
   \\big(\\frac{\\partial^2u}{\\partial x^2}+\\frac{\\partial^2u}{\\partial y^2}\\big)$ (9.1); here x runs along the surface and y
   away from it, valid while the layer is thin compared with the radius of curvature R of the surface (δ/R ≪ 1). It is the
   starting line of D01.", where="Ch. 4 §4.10, Ch. 8 §8.1")`
3. `nb.recap("R02", "Scaled continuity", "The two-length scaling of Ch. 8 turned continuity into
   $\\frac{\\partial u^*}{\\partial x^*}+\\frac{\\partial v^*}{\\partial y^*}=0$ (8.15) with no coefficient: v is scaled so that
   both terms have the same size. The same move returns in D01 step 6.", where="Ch. 8 §8.3")`
4. `nb.recap("R03", "Continuity", "Incompressible mass conservation $\\frac{\\partial u}{\\partial x}+\\frac{\\partial v}
   {\\partial y}=0$ (6.2); it guarantees a stream function with u = ∂ψ/∂y, v = −∂ψ/∂x (Ch. 4).", where="Ch. 6 §6.2")`
5. `nb.core("C01", "The boundary-layer equations $u u_x+v u_y=-\\frac1\\rho p_x+\\nu u_{yy}$ (9.9), $p_y=0$ (9.10) and the outer
   matching (9.11)", question="How thin is the layer at a wall, and which terms of the Navier–Stokes equations survive in it?")`
6. `nb.md` — **The problem in plain words:** "Blow air over a 1 m board at 1 m/s. Far from the board the air is nearly
   untouched; at the board it is at rest. Between them lies a layer only a few millimetres thick in which the speed climbs
   from 0 to 1 m/s. How thick? And can we write simpler equations for it than the full Navier–Stokes equations? If we can,
   we can solve them — and the drag, the heat transfer and the separation of the flow all follow from that layer. (Ideal-flow
   theory has no layer at all: it predicts zero drag, d'Alembert's paradox of Ch. 6.)"
7. `nb.md` — **The idea** (ASCII, to scale in the ratio 1 : 10 vertically):
   ```
   outer flow: ideal, irrotational, U_e(x)  ──►  hands the layer its pressure p(x)
   ─────────────────────────────────────────────────────────────
            δ(x) ~ √(νx/U)  (grows like √x, thin: δ/L ~ Re^(−1/2))
   ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─
     u(y):  0 ─▶ ─▶─▶─▶──▶ (climbs to U_e inside δ)      viscosity matters here only
   ═════════════════════════ wall, u = v = 0 ═══════════════════
   ```
8. `nb.note` — **N01 [B]** "**Prandtl's hypothesis.** At large Re viscosity acts only in a thin layer next to a solid
   surface (inner problem, no slip); outside it the flow is irrotational and inviscid (outer problem). The layer is not a
   separate fluid: it is where vorticity created at the wall (Ch. 5, Ch. 8) has diffused to. The zero drag of Ch. 6 fails
   because the ν → 0 limit is singular: a tiny viscosity still enforces no slip." + `nb.figure` ("Fig. 9.1 remade",
   `plate_layer(U=1, nu=NU_AIR, L=1, exaggeration=1)` and a second copy with the layer drawn to scale — thickness 3.9 mm on a
   1 m plate): *see* "a thin sliver along the plate, growing downstream"; *read* "the edge is the δ₉₉ curve; above it streamlines
   are the ideal ones"; *change* "…the plate were 4 m long: the layer at the end is only 2× thicker (√4)".
9. `nb.primer("six Reynolds numbers (which length?)", "The Reynolds number is a ratio of inertia to viscosity, but which length
   and speed? This chapter uses six flavours: the overall Re = U∞L/ν of (9.6) (body length L), the local Re_x = Ux/ν (distance
   from the leading edge), the plate's Re_L = UL/ν, the jet's Re_x = x u₀/ν and Re_h99, the cylinder and sphere Re = U∞d/ν
   on the **diameter**, and the Falkner–Skan Re_x = a x^{n+1}/ν. A statement 'Re ≫ 1' is only meaningful with its length.",
   code="U, L, nu = 1.0, 1.0, 1.5e-5\nprint(U*L/nu)              # 66667: overall Re of a 1 m board at 1 m/s in air\nprint(U*0.1/nu)            # 6667: local Re_x at 10 cm from the leading edge (smaller: the layer is younger)\nprint(1.0*0.02/nu)         # 1333: a 2 cm cylinder in the same wind (diameter based)")` (**P200**)
10. `nb.note` — **N02 [B], N03 [B]** "**Sizes of the two terms that compete** (U = 1 m/s, L = 1 m, ν = 1.5 × 10⁻⁵ m²/s). The
   advective term $u\\,\\partial u/\\partial x\\sim U^2/L$ (9.2) = 1 m/s². The viscous term $\\nu\\,\\partial^2u/\\partial
   y^2\\sim\\nu U/\\bar\\delta^2$ (9.3) equals 15 m/s² if the layer is 1 mm thick, 0.15 m/s² if 1 cm: the layer thickness δ̄
   is whatever makes them equal." equations (9.2), (9.3).
11. `nb.note` — **N04 [B]** "**Thickness estimate (9.4):** $\\bar\\delta\\sim\\sqrt{\\nu L/U}$ or $\\bar\\delta/L\\sim Re^{-1/2}$.
   Same square root as ch08's $\\sqrt{\\nu t}$ with t = L/U: the fluid spends a time L/U passing the board and momentum diffuses
   √(νL/U) in that time (bridge to Ch. 8 §8.4). Air over 1 m at 1 m/s: 3.9 mm; water: 1.0 mm." equation (9.4).
12. `nb.note` — **N05 [B], N06 [B]** "**Derivative sizes and v:** $\\partial/\\partial x\\sim1/L$, $\\partial/\\partial y\\sim1/\\bar\\delta$
   (9.5); continuity (6.2) then says $U/L\\sim v/\\bar\\delta$, so $v\\sim U\\,Re^{-1/2}$ — 3.9 mm/s, 250 times smaller than u."
   equation (9.5).
13. `nb.derivation("D01", …)` — Part F D01 (14 steps), ref "9.7". (Uses P188 two-length scaling, P49 chain rule, P198 dominant
    balance, P200.)
14. `nb.note` — **N07 [B]** "**The stretched variables (9.6)** $x^*=x/L,\\ y^*=\\frac yLRe^{1/2},\\ u^*=u/U,\\ v^*=\\frac vURe^{1/2},\\
   p^*=\\frac{p-p_\\infty}{\\rho U^2}$ — the ε = Re^{−1/2} of ch08's (8.14), but with the pressure scaled by ρU² (advective) rather
   than μU/δ (viscous, which gave the lubrication and Stokes limits (8.42))." equation (9.6).
15. `nb.note` — **N08 [B], N09 [B]** "**The scaled equations:** $u^*\\frac{\\partial u^*}{\\partial x^*}+v^*\\frac{\\partial u^*}{\\partial y^*}=
   -\\frac{\\partial p^*}{\\partial x^*}+\\frac1{Re}\\frac{\\partial^2u^*}{\\partial x^{*2}}+\\frac{\\partial^2u^*}{\\partial y^{*2}}$ (9.7) and
   $\\frac1{Re}\\big(u^*\\frac{\\partial v^*}{\\partial x^*}+v^*\\frac{\\partial v^*}{\\partial y^*}\\big)=-\\frac{\\partial p^*}{\\partial y^*}+
   \\frac1{Re^2}\\frac{\\partial^2v^*}{\\partial x^{*2}}+\\frac1{Re}\\frac{\\partial^2v^*}{\\partial y^{*2}}$ (9.8). > ⚠️ **The book prints
   (9.7) with denominators ∂x* and ∂y* (no squares); dimensionally and by D01 step 9 the correct ones are ∂x*² and ∂y*².**
   The second-derivative term in y keeps coefficient 1 — the difference from (4.101) and (8.42), where the two directions were
   scaled alike." equations (9.7), (9.8).
16. `nb.worked_example("air over a 1 m board at 1 m/s", "1. Re = UL/ν = 1 × 1/1.5 × 10⁻⁵ = 6.67 × 10⁴. 2. δ̄/L ~ Re^{−1/2} =
    1/258 = 3.87 × 10⁻³, so δ̄ ≈ 3.9 mm. 3. v ~ U Re^{−1/2} = 3.9 mm/s. 4. The dropped x-diffusion term ν U/L² = 1.5 ×
    10⁻⁵ m/s² against the kept ν U/δ̄² = 1.0 m/s²: ratio 1/Re = 1.5 × 10⁻⁵. 5. Crude wall stress τ₀ ~ μU/δ̄ = 1.8 × 10⁻⁵ ×
    1/3.87 × 10⁻³ = 4.6 × 10⁻³ Pa and C_f ~ 2/√Re = 7.7 × 10⁻³ (N16 — Blasius will give 0.664/√Re, a factor 3 smaller).")`
17. `nb.code` — `sc = BL.boundary_layer_scales(U=1.0, L=1.0, nu=NU_AIR, rho=RHO_AIR)`; print the dict. *expect:* Re = 66667,
    delta_over_L = 3.873e-3, delta = 3.873e-3 m, v_scale = 3.873e-3 m/s, tau0_scale = 4.648e-3 Pa, cf_estimate = 7.746e-3,
    adv = 1.0, visc = 1.0, visc_x = 1.5e-5. *explain:* 1. the sizes (9.2)–(9.4); 2. `adv ≈ visc` is the definition of δ̄;
    3. `visc_x` is the term (9.7) drops.
18. `nb.code` — sympy engine: `res = ch09.bl_nondim_sympy()`; print `res["coefficients"]` and `res["limit"]`; then
    `bad = ch09.bl_nondim_sympy(printed_9_7=True)`; print `bad["dimension_check"]`. *expect:* coefficients (powers of Re):
    x-momentum d²u/dx*² → −1, d²u/dy*² → 0; y-momentum inertia → −1, pressure 0, d²v/dx*² → −2, d²v/dy*² → −1; continuity 0;
    the printed form fails ("second derivative not dimensionless"). *explain:* the engine repeats D01 by machine.
19. `nb.note` — **N10 [B]** "**(9.10) $\\partial p/\\partial y=0$:** pressure is uniform across the layer, so the wall pressure is
   the edge pressure, which the outer ideal flow supplies — the experimental fact behind 'surface pressure of an attached layer
   equals ideal-flow theory'. It fails at separation and where the surface curves as sharply as the layer is thick." + `nb.code`
    `[BL.bl_pressure_variation(Re)["ratio"] for Re in (1e3, 1e4, 1e5, 1e6)]`. *expect:* falling roughly like 1/Re (a
    log–log slope printed with `observed_order`, between −0.5 and −1.1; the notebook asserts only ratio(1e6) < 1e-2 and monotone
    decrease). equation (9.10).
20. `nb.note` — **N11 [B]** "**(9.11) matching:** the pressure gradient inside the layer comes from Bernoulli applied to the
   outer flow, $-\\frac1\\rho\\frac{dp}{dx}=U_e\\frac{dU_e}{dx}$ (9.11); the outer speed U_e(x) is the only forcing of every
   layer in this chapter. Decelerating outer flow (dU_e/dx < 0) means dp/dx > 0: **adverse**." equation (9.11).
21. `nb.derivation("D02", …)` — Part F D02 (5 steps), ref "9.11". + `nb.code`: `of = BL.outer_flow("wedge", n=1.0, a=10.0)`;
    `print(of.Ue(0.1), of.dUe(0.1), of.dpdx(0.1, rho=RHO_AIR))`. *expect:* 1.0 m/s, 10 s⁻¹, −12 Pa/m (favourable).
22. `nb.note` — **N12 [B], N13 [B]** "**Conditions on the layer:** $u(x,y\\to\\infty)=U_e(x)$ (9.14), matching to the outer flow
   (y → ∞ means y ≫ δ); and an inlet profile $u(x_0,y)=u_{in}(y)$ (9.15) at some x₀ — needed because (9.9) is parabolic. The
   wall conditions are the recap below." equations (9.14), (9.15).
23. `nb.primer("parabolic, elliptic and marching in x", "A PDE is *elliptic* if conditions all round the domain are needed (full
   Navier–Stokes: what happens downstream affects what happens upstream); *parabolic* if it has a time-like direction and a
   start (the heat equation, and (9.9) with x in the role of time: u ∂u/∂x = ν ∂²u/∂y² + …). A parabolic problem is solved
   by *marching* from an inlet profile, one x-step at a time — information only flows downstream.", code="import numpy as np\nu = np.zeros(7); u[3] = 1.0           # a spike at the middle node of seven\nfor step in range(3):                  # march the heat equation 'in time' (or in x)\n    u[1:-1] += 0.25*(u[:-2] - 2*u[1:-1] + u[2:])\n    print(u)                           # the spike spreads outwards and flattens; nothing comes back")` (**P201**)
24. `nb.note` — **N14 [B]** "**Parabolic character.** In (9.9) the pair $u\\,\\partial u/\\partial x$ and $\\nu\\,\\partial^2u/\\partial y^2$ makes x
   a time and y a space: like Stokes' first problem (Ch. 8) the layer at x is determined by the layers upstream and the outer
   pressure, never by anything downstream. This is why we can march (and why marching stops at separation — C08)." + `nb.code`:
    `x = np.linspace(0.02, 1.0, 60)`; `march = BL.march_boundary_layer(BL.outer_flow("flat", U=1.0).Ue, x, NU_AIR,
    u_inlet=<Blasius profile at x[0] from BL.blasius_fields>)`; compare `march["tau0"]` with `BL.blasius_wall_shear(x, 1.0,
    RHO_AIR, NU_AIR)`. *expect:* relative difference < 1e-3 for x ≥ 0.1 m; the observed order in Δx printed (≈ 1, backward
    Euler). *explain:* 1. inlet profile; 2. `march_boundary_layer` steps in x with the von Mises variable; 3. τ₀ agrees with the
    exact similarity solution that C03–C04 will derive.
25. `nb.check_agree` — **from scratch (curation §7):** on the Blasius field (`BL.blasius_fields` on a 200 × 120 grid,
    x ∈ [0.5, 1] m, y ∈ [0, 3δ₉₉]) evaluate the terms of (9.1) by `np.gradient` (second order): advective $uu_x+vu_y$, viscous
    $\nu u_{yy}$, dropped $\nu u_{xx}$ at U = 0.1, 1, 10 m/s (Re_L = 6.7 × 10³, 6.7 × 10⁴, 6.7 × 10⁵). Print the ratios
    (viscous/advective, dropped/kept). `assert np.allclose(residual_9_9, 0, atol=1e-3 * scale)` with `BL.bl_x_momentum_residual`
    (the similarity solution satisfies (9.9)); `assert` dropped/kept × Re_x stays within [0.1, 10]. *expect:* viscous/advective
    ≈ 1 (0.5–2) at all three speeds; dropped/kept ≈ 1/Re_x·(0.3–3).
26. `nb.figure` — **term sizes and layer to scale** (9 × 3.4 in, three panels): (a) the layer to scale (δ₉₉(x), x ∈ [0, 1] m,
    air at 1 m/s; 20 mm scale bar) with the outer streamline dashed teal; (b) grouped bars (log axis) of advective, kept
    viscous ν u_yy, dropped ν u_xx and pressure difference at Re_L = 6.7 × 10³ / 10⁵ / 10⁷ (blue, rose, muted, orange) from the
    Blasius field; (c) `march` τ₀ (teal dots) on the Blasius line (muted dashed) with the ratio in the title. Title "Only the
    x-diffusion term is small: the rest all matter inside the layer". *see:* "three bar groups where the muted bar shrinks
    tenfold per hundredfold Re"; *read:* "advective ≈ viscous in every group: that is the definition of δ̄; dropped/kept falls
    like 1/Re"; *change:* "…Re were 10 (a short board in oil): all four bars comparable — the layer is not thin and (9.9) fails
    (N17)".
27. `nb.note` — **N16 [B]** "**A crude wall stress:** $\\tau_0\\sim\\mu U/\\bar\\delta$ gives $C_f=\\tau_0/(\\tfrac12\\rho U^2)\\sim2/\\sqrt{Re}$ (order and
   Re-dependence right; Blasius' exact 0.664 in C04 is a factor 3 smaller)."
28. `nb.pointer` — **N15 [C]** "The two-step viscous–inviscid procedure (ideal flow ⇒ pressure ⇒ layer ⇒ displaced body ⇒ repeat) is used
    in C02 (δ* correction) and in Chs. 10, 14." **N17 [C]** "Where the approximation fails: near the leading edge (Re_x ≲ 1), where δ/R is
    not small, and after separation — shown in C04 (leading-edge τ₀ ∝ x^{−1/2}) and C08."
29. `nb.md` — **What would change if…** "…we asked for the drag: we need the wall stress, which needs the profile, not only its
    thickness. C02 defines what 'thickness' means when the profile fades smoothly; C03 solves the equations for the simplest
    outer flow, a constant one."
30. `nb.recap("R04", "No slip", "Eq. (9.12) $u(x,0)=0$: the fluid sticks to a solid wall (Ch. 4 §4.10, Ch. 8 §8.2) — the reason
    a layer exists.", where="Ch. 4 §4.10, Ch. 8 §8.2")` · `nb.recap("R05", "No through-flow", "Eq. (9.13) $v(x,0)=0$: no fluid
    crosses a solid wall (Ch. 8 (8.2)); suction through the wall (Exercise 9.26) would change this and is only named.",
    where="Ch. 8 §8.2")`

### A.2 §9.2 Boundary-Layer Thickness Definitions — C02 (+N18–N21 · D03 D04 · P202 P203 P204 P205 · A1 · F1 · E1)
1. `nb.section("9.2", "Boundary-Layer Thickness Definitions", intro="**What is this section about?** The profile fades smoothly
   into the outer flow, so 'where does the layer end?' has no unique answer. There are three useful answers: where the speed
   is 99 % of the outer speed (δ₉₉, a position), how far the outer streamlines are pushed away from the wall (δ*, the
   displacement thickness) and how much momentum the wall has stolen (θ, the momentum thickness). All later results are
   quoted in these three lengths.")`
2. `nb.core("C02", "Displacement thickness $\\delta^*=\\int_0^\\infty(1-u/U_e)dy$ (9.16), with $\\delta_{99}$ and $\\theta$ (9.17)",
   question="What is the thickness of a layer that has no sharp edge — and what does the layer do to the flow outside it?")`
3. `nb.md` — **The problem in plain words:** "A wind-tunnel wall carries a slow layer. The fast air in the middle is squeezed as if
   the tunnel were narrower; the designer widens the walls by the layer's 'thickness' to keep the test speed. But which
   thickness? The 99 % point ignores how slow the fluid inside is. We want the thickness of a *fictitious* layer that would
   have the same effect on the outer flow."
4. `nb.md` — **The idea** (ASCII): "real profile ≈ 'stagnant slab of thickness δ*' + full-speed flow above; and 'slab of
   thickness θ' for momentum":
   ```
    y                        y                          y
    │      ──────U_e         │  ░░░░ deficit U_e−u      │  ░░░ loss u(U_e−u)
    │    ╭──                 │ ░░░░                     │  ▒▒▒
    │  ╭─╯  u(y)             │░░░░  area = U_e·δ*       │▒▒▒  area = U_e²·θ
    └──╨──────────  u        └─────────────  u          └─────────── u
   ```
5. `nb.primer("improper integral of a deficit (truncating the tail)", "An integral to infinity such as ∫₀^∞ (1 − u/U_e) dy
   converges when the integrand dies fast enough (for the layer it dies like a Gaussian or an exponential). Numerically
   we integrate to a y_max where the integrand is negligible and add the tail (an exponential tail is integrated exactly).",
   code="import numpy as np\ny = np.linspace(0, 30, 3001)\nd = np.exp(-y)                  # the deficit 1 - u/U for u/U = 1 - exp(-y/a), a = 1 m\nprint(np.trapezoid(d, y))        # 1.0000: the integral to 30 m is already the whole answer\nprint(np.trapezoid(d[y <= 5], y[y <= 5]))   # 0.9933: cutting at 5 m loses 0.7 %, the tail")` (**P202**)
6. `nb.primer("scipy.integrate.simpson and np.trapezoid", "Both integrate a sampled profile. `np.trapezoid(f, y)` joins samples by
   straight lines (error ∝ Δy²); `scipy.integrate.simpson(f, x=y)` fits parabolas through triples (error ∝ Δy⁴), better for
   smooth profiles. numpy 2 has no `np.trapz`; use `np.trapezoid`.", code="import numpy as np\nfrom scipy.integrate import simpson\ny = np.linspace(0, 1, 11)\nprint(np.trapezoid(y**3, y), simpson(y**3, x=y))   # 0.2525 vs 0.25: exact answer 1/4")` (**P203**)
7. `nb.primer("monotone interpolation PchipInterpolator", "To read the height where u/U = 0.99 from samples we interpolate the profile
   *without overshoot*: `PchipInterpolator` keeps monotone data monotone (a cubic spline may wiggle). Then `brentq` finds the
   crossing.", code="import numpy as np\nfrom scipy.interpolate import PchipInterpolator\nfrom scipy.optimize import brentq\ny = np.linspace(0, 5, 11); u = 1 - np.exp(-y)\nP = PchipInterpolator(y, u)\nprint(brentq(lambda s: P(s) - 0.99, 0, 5), np.log(100))   # 4.6 vs exact 4.605: the 99 % height")` (**P204**)
8. `nb.primer("control volume with a streamline as a side", "In a control volume (Ch. 4 §4.4) mass and momentum cross faces. A
   *streamline* can be one face because by definition no fluid crosses it: it carries no mass flux and no momentum flux, and
   the pressure on it is the ambient one. Use it as the 'roof' of a box over a wall layer.", code="import numpy as np\n# inflow through x=0 (height h0, speed U) must equal outflow through x=L (height h(L), profile u):\nU, h0 = 2.0, 0.10\nflux_in = U*h0                          # m^2/s per unit width\nprint(flux_in)                          # 0.2: the streamline rises to wherever u(y) integrates to this")` (**P205**)
9. `nb.note` — **N18 [B]** "**δ₉₉** $u(x,\\delta_{99})=0.99\\,U_e(x)$: convenient, arbitrary (95 %, 99.9 % are used too), and
   ignores what the profile does inside. `BL.delta_level(y, u, Ue, level=0.99)` returns it by PCHIP and `brentq`."
10. `nb.derivation("D03", …)` — Part F D03 (6 steps), ref "9.16". Then `nb.note` — **N19 [B]** "**δ* moves the outer flow.** Because
    ∂u/∂x < 0 in a growing layer, continuity gives v(y) = −∫₀^y u_x dy′ > 0: the outer streamlines are lifted by about δ*(x),
    so the outer flow 'sees' the body plus a slab of thickness δ*. Airfoil, duct and inlet designers correct the ideal
    calculation by δ*. (Numbers in the figure below: Blasius, v_∞ = 0.86 U/√Re_x = U dδ*/dx.)"
11. `nb.note` — **N20 [B]** "**θ:** $\\theta=\\int_0^\\infty\\frac u{U_e}\\big(1-\\frac u{U_e}\\big)dy$ (9.17) — the momentum flux missing compared with an
    ideal layer is ρU_e²θ. It is a *length* (in Chs. 3, 6, 8 θ is an angle!). The shape factor **H = δ*/θ** measures how 'full'
    the profile is: 2.59 for Blasius, ≈ 3.5–4 near separation, ≈ 1.3 for a turbulent layer (Ch. 12)." equation (9.17).
12. `nb.derivation("D04", …)` — Part F D04 (9 steps), ref "9.17". Then `nb.note` — **N21 [B]** "**Plate drag from θ:** the derivation
    gives $\\rho U^2\\theta(x)=\\int_0^x\\tau_0\\,dx$: measure the momentum defect far behind a plate and you have its drag (the wake
    method of Ch. 4, `ch04.wake_drag_per_span`)." equation $\\rho U^2\\theta=\\int_0^x\\tau_0dx$.
13. `nb.worked_example("two profiles by hand", "1. u/U = 1 − e^{−y/a}: δ* = ∫e^{−y/a}dy = a; θ = ∫(1 − e^{−y/a})e^{−y/a}dy = a −
    a/2 = a/2; H = 2; the 99 % height is a ln 100 = 4.605a. 2. Linear u/U = y/δ (cut at δ): δ* = δ/2, θ = ∫(y/δ)(1 − y/δ)dy = δ/6,
    H = 3, δ₉₉ = 0.99δ. Same δ₉₉ to 1 %-ish, yet δ*/δ₉₉ = 1 versus 0.5: δ₉₉ tells you nothing about the displacement.")`
14. `nb.code` — `y = np.linspace(0, 40e-3, 4001)`; `bl = BL.blasius_fields(1.0, y, 1.0, NU_AIR)`; `th = BL.thicknesses(y, bl["u"], 1.0)`;
    print th; also `BL.blasius_constants()`; for the four shapes of `BL.profile_shape`. *expect (air, U = 1 m/s, x = 1 m):* δ₉₉ =
    19.02 mm, δ* = 6.665 mm, θ = 2.572 mm, H = 2.591 (= 4.910, 1.7208, 0.6641 in √(νx/U) = 3.873 mm); shapes in units of δ:
    linear (0.5, 0.1667, 3.00), sine (0.3634, 0.1366, 2.660), cubic (0.3750, 0.1393, 2.692), exponential (1.000, 0.5000, 2.000)
    in units of a. *explain:* 1. `blasius_fields` gives u(y) (C03–C04 explain it); 2. `thicknesses` integrates with Simpson and
    a tail; 3. `profile_shape` are the exact closed forms the numerical routine must reproduce.
15. `nb.check_agree` — **from scratch (curation §7):** trapezoid rule by hand for δ* and θ on 1 − e^{−y/a} (a = 1 mm) and on
    the sampled Blasius profile; `assert np.isclose(delta_star_mine, 1e-3, rtol=1e-5)`, `np.isclose(theta_mine, 0.5e-3, rtol=1e-5)`;
    `assert np.allclose(mine, BL.thicknesses(...)["delta_star"], rtol=1e-6)`; plus the streamline lift: `dδ*/dx` from
    finite differences of `BL.blasius_delta_star(x, U, nu)` equals `v_inf = 0.8604 U/√Re_x` (`assert np.isclose(…, rtol=1e-6)`).
16. `nb.figure` — **three thicknesses on one profile** (9 × 3.4 in, three panels): (a) Blasius u/U against y [mm] (blue) with the
    three heights marked (δ₉₉ muted dashed, δ* amber, θ orange); (b) the deficit 1 − u/U shaded amber (area = δ*) and the
    momentum loss (u/U)(1 − u/U) shaded orange (area = θ); (c) streamlines above the plate: ideal (dashed teal) versus real
    (solid), lifted by δ*(x) with the amber line y = δ*(x) drawn. Title "One profile, three thicknesses". *see:* "the amber and
    orange areas; the real streamlines rising above the ideal ones"; *read:* "the amber area equals U_e δ*: δ* = 6.7 mm is 1/3
    of δ₉₉ = 19 mm; H = 2.59"; *change:* "…the profile were fuller (turbulent-like): δ* and θ shrink relative to δ₉₉ and H → 1.3".
17. `nb.animation` — **A1** (video, 60 frames, FAST 30): the layer growing along a 2 m plate; the profile u(y) at moving station x,
    δ₉₉ curve, the displaced streamline line y = δ*, tracers lifting, and the marching-solver dots (C01) riding the Blasius
    profile. Figure notes: *see* "the profile fattens as the station moves; the amber line lifts like √x"; *read* "at x = 1 m
    the marker sits at 19 mm; at 4 m at 38 mm (×2)"; *change* "…water: the same movie at 1/3.9 the height".
18. `nb.plotly` — **F1** `slider_figure` (≤ 25 steps): the power-law family u/U = (y/δ)^{1/p} for p = 1 … 10: profile plus a trace of
    δ*/δ, θ/δ, and H = (p + 2)/p (p = 1: H = 3; p = 7: H = 1.29, a turbulent-like profile; Blasius sits at p ≈ 1.26). *expect:*
    δ*/δ = 1/(p + 1), θ/δ = p/((p + 1)(p + 2)).
19. `nb.explainer("bl_scaling_thicknesses", heading="How thin is a boundary layer — and what is its thickness?", why="Sliding U, ν
    and x shows δ ∝ √(νx/U) and Re^(−1/2) in one picture, and switching the profile shape shows the three thicknesses change by
    different ratios — a static figure shows one case.", tries=["Set air, U = 1 m/s and drag x from 0.1 to 2 m: watch δ₉₉ grow like √x.",
    "Switch the profile shape from Blasius to linear: which of δ*, θ, H changes most?", "Choose the preset 'small Re near the leading
    edge': the status warns that the layer is not thin.", "Open Derivation and step through D01: the picture shows the stretched y."])`
20. `nb.md` — **⚠️ Common confusion:** "δ₉₉, δ* and θ are three answers to three questions; none is 'the' thickness. δ = √(νx/U)
    of C03 is a *scale*, not any of them (δ₉₉ = 4.91 δ)."
21. `nb.md` — **What would change if…** "…the outer speed changes along x? The same three integrals apply with U_e(x); C06 turns
    them into one ODE, (9.43), that any layer must obey. First, the simplest case: a constant U_e over a flat plate."

---
### A.3 §9.3 Boundary Layer on a Flat Plate: Blasius Solution — R06 R07, C03 (+N22–N32 · N24 · D05 · P206), C04 (+N33–N42 · D06 · P207 P208 · A2 · F2 · E2)
1. `nb.section("9.3", "Boundary Layer on a Flat Plate: Blasius Solution", intro="**What is this section about?** The simplest
   outer flow is no flow change at all: a thin plate lined up with a steady stream U. Then the pressure gradient vanishes and
   the layer equations (9.9)–(9.10) have an exact solution. Because a semi-infinite plate has no length of its own, the
   profiles at different distances are stretched copies of one curve; that turns the PDE into an ODE (C03), which a computer
   solves in milliseconds and which gives every number we need: thickness, wall shear, drag (C04). Heinrich Blasius did this
   in 1908; the same trick will work for wedges and jets.")`
2. `nb.recap("R06", "Wall conditions", "On the plate both velocity components vanish, $u=v=0$ on $y=0$ (9.20), i.e. no slip
   (9.12) $u(x,0)=0$ and no through-flow (9.13) $v(x,0)=0$ (R04, R05 above).", where="§9.1")` · `nb.recap("R07", "Edge condition",
   "Far from the wall the velocity joins the outer stream: $u\\to U$ as $y/\\delta\\to\\infty$ (9.21), the constant-speed case of the
   matching condition (9.14).", where="§9.1")`
3. `nb.core("C03", "Blasius equation $f'''+\\frac12ff''=0$ (9.27) by similarity reduction of (9.18)", question="Can the whole
   flat-plate layer — a PDE in x and y — collapse to one ordinary differential equation?")`
4. `nb.md` — **The problem in plain words:** "At x = 0.1 m the layer on the plate is 1.9 mm thick, at 1 m it is 19 mm. Are the
   profiles at the two places different curves — or the same curve drawn at two magnifications? A semi-infinite plate has no
   length scale except the distance x itself, so the second answer is the only one that makes sense. If it holds, one function
   f(η) describes the entire plate."
5. `nb.md` — **The idea** (ASCII): "constant-η lines are parabolas y ∝ √x; each station is the same picture zoomed in y only:"
   ```
   y     η = 4    ╱──────────────                 profile at x₁:  ▏╱‾‾   ← same shape,
   │   η = 2  ╱──╯                                profile at x₂:  ▏╱‾‾‾‾  stretched by √(x₂/x₁)
   │  η = 1 ╱─╯  ...   δ(x) = √(νx/U)             ψ = Uδ(x) f(y/δ(x))   →   u = U f′(η)
   └─────────────────────── x
   ```
6. `nb.note` — **N22 [B]** "**The equation to be solved:** with dp/dx = 0 (U_e = U) the layer equation (9.9) becomes
   $u\\frac{\\partial u}{\\partial x}+v\\frac{\\partial u}{\\partial y}=\\nu\\frac{\\partial^2u}{\\partial y^2}$ (9.18) together with
   continuity (6.2). Momentum enters only through viscosity: (9.18) says the fluid element at height y loses speed to its slower
   neighbour below." equation (9.18). > ⚠️ **Common confusion (name clash):** "This Blasius solution (1908, boundary layers) is not
   the *Blasius force theorem* (6.60) of Ch. 6 (1910, forces on cylinders by contour integration) — same person, different
   result; both return in Ch. 14."
7. `nb.note` — **N23 [B]** "**The similarity ansatz.** Nothing in the problem has a length except x, so we look for
   $\\psi=U\\delta(x)f(\\eta)$, $\\eta=y/\\delta(x)$ (9.19): ψ has dimensions m²/s = (speed) × (length), so δ(x) must be a length and f
   dimensionless. It specialises ch08's (8.32) $u=U F(y/\\delta(t))$ with time replaced by x/U; the price is that ψ (not u) carries the
   δ so that continuity holds automatically." equation (9.19).
8. `nb.primer("chain rule when the similarity variable moves with x and y", "η = y/δ(x) depends on *both* variables. For a function
   F(η): ∂F/∂y = F′(η)·(1/δ) and ∂F/∂x = F′(η)·∂η/∂x = F′(η)·(−y δ′/δ²) = −F′(η)·η δ′/δ. The chain rule applies each partial
   derivative separately; a product like δ(x) f(η) also needs the product rule in x.", code="import sympy as sp\nx, y = sp.symbols('x y', positive=True)\ndelta = sp.Function('delta')(x)\nf = sp.Function('f')\neta = y/delta\nprint(sp.simplify(sp.diff(f(eta), y)))     # f'(eta)/delta  (P49 chain rule, one variable)\nprint(sp.simplify(sp.diff(f(eta), x)))     # -y f'(eta) delta'/delta^2 = -eta f' delta'/delta")` (**P206**)
9. `nb.derivation("D05", …)` — Part F D05 (13 steps), ref "9.27". (Contains N25 δ → 0 at x → 0 (9.22), N26 u = ψ_y = U f′ (9.23), N27 v = −ψ_x = Uδ′(ηf′ − f) (9.24), N28 the ηf′f″ cancellation, N29 the bracket equation (9.25), N30 δ = √(νx/U) (9.26); boundary conditions N31 (9.28) as its last step.)
10. `nb.note` — **N25–N30 [B] summary** (after the derivation, a compact box the reader can reuse): "$u=U f'(\\eta)$ (9.23), $v=U\\delta'(\\eta f'-f)$
    (9.24), $-\\big[\\frac{U^2}{\\delta}\\delta'\\big]ff''=\\big[\\frac{\\nu U}{\\delta^2}\\big]f'''$ (9.25), $\\delta(x)=[\\nu x/U]^{1/2}$ (9.26),
    $\\frac{d^3f}{d\\eta^3}+\\frac12f\\frac{d^2f}{d\\eta^2}=0$ (9.27) — third order, nonlinear, with conditions $f=f'=0$ at η = 0 (9.28)
    and $f'\\to1$ as η → ∞ (9.29): two at the wall, one at infinity." equations (9.23)–(9.29).
11. `nb.worked_example("is a candidate δ(x) allowed?", "The bracket test: the ODE keeps its x-independence only if $U\\delta\\delta'/\\nu$
    is a constant. (a) δ = √(νx/U): δδ' = ν/(2U), so the test number is 1/2 at every x ✓. (b) δ = kx with k = 0.01 (a wedge-like
    layer): U k² x/ν = 1 × 10⁻⁴ x/1.5 × 10⁻⁵ = 3.33 at x = 0.5 m but 6.67 at x = 1 m ✗ — the ODE would change with x, no similarity. (c) δ at x
    = 0.5 m for air at U = 1 m/s: √(1.5 × 10⁻⁵ × 0.5) = 2.74 mm, and y = 2δ = 5.48 mm has η = 2.")`
12. `nb.code` — `red = ch09.similarity_reduce_sympy("blasius")`; print `red["brackets"]`, `red["cancelled_terms"]`, `red["ode"]`, `red["delta"]`,
    `red["residual"]`. *expect:* brackets $-U^2\delta'/\delta$ and $\nu U/\delta^2$; cancelled terms η f′ f″ (coefficient $\mp U^2\delta'/\delta$);
    ode `f''' + f*f''/2 = 0`; delta `sqrt(nu*x/U)`; residual 0. *explain:* 1. sympy differentiates the ansatz for u and v exactly as
    steps 3–6 of D05; 2. it lists the two terms that cancel; 3. it demands proportional brackets and solves for δ; 4. the residual
    of the PDE with the ansatz inserted and the ODE substituted is zero.
13. `nb.pointer` — **N24 [C]** "Favourable gradients (U_e′ > 0) forget the inlet profile, adverse ones never (Serrin, Peletier); Blasius sits
    on the border. Below, the marching solver of C01 started from two different inlets ends on the same profile — qualitative
    only; Ch. 10 has the numerics." + `nb.code`: `march` from a linear and from a sine inlet profile at x₀ = 0.02 m (U_e = const) to
    x = 1 m; print `max |u₁ − u₂|/U` at x = 0.1, 0.5, 1. *expect:* falling below 1e-3 by x = 1 m (qualitative; printed, not asserted).
14. `nb.figure` — **the similarity test and the zoom map** (8 × 3.4 in, two panels): (a) the test number $U\delta\delta'/\nu$ against x for
    δ ∝ x^{1/2} (flat at 0.5, teal), δ ∝ x (rising, rose), δ ∝ x^{1/3} (falling, amber) on x ∈ [0.05, 1] m; (b) the x–y plane with the
    lines η = 1, 2, 3, 4, 5 (parabolas, purple dashed) and the δ₉₉ line η = 4.91 (muted), three stations x = 0.1, 0.5, 1 m marked with vertical
    ticks. Title "Only δ ∝ √x makes the ODE independent of x". *see:* "one flat line, two drifting ones; parabolas fanning out"; *read:*
    "on the fan every station meets the same η at heights that scale like √x"; *change:* "…U doubled: the fan flattens by √2; the
    flat line stays at 1/2".
15. `nb.md` — **What would change if…** "…the outer speed varied as a power of x? The same steps with U → U_e(x) give C05. First, let's
    solve (9.27) and read off the numbers."

16. `nb.core("C04", "Blasius solution: $f''(0)=0.332$ and the numbers $\\delta_{99},\\delta^*,\\theta,\\tau_0$ (9.30)–(9.33)",
    question="What is the profile — and what wall shear and drag does it give?")`
17. `nb.md` — **The problem in plain words:** "A glider wing, a ship's hull, the floor of a wind tunnel: the friction drag of a smooth
    surface at moderate Re is set by one number, the slope of the velocity profile at the wall. Blasius' ODE has that slope as its
    only unknown: f″(0) = 0.332. Get it once and the skin friction, the drag and all three thicknesses follow."
18. `nb.md` — **The idea** (ASCII): "a two-point boundary problem is awkward (we know f, f′ at η = 0 and f′ at ∞), but the equation has a
    scaling symmetry, so we shoot once and rescale":
    ```
    solve  f''' + ½ f f'' = 0,  f(0) = f'(0) = 0,  f''(0) = 1     (initial-value problem)
    read   f'(∞) = 2.0854 ≠ 1
    rescale  f(η) = λ f₁(λη), λ² · 2.0854 = 1  ⇒  f''(0) = λ³ = 0.33206
    ```
19. `nb.pointer` — **R06/R07 reminder** in one line: the conditions (9.28) and (9.29) are (9.20) and (9.21) written in η.
20. `nb.primer("scaling symmetry of an ODE and the Töpfer trick", "If f₁(η) solves the ODE, sometimes λ f₁(λη) does too (each term
    scales with the same power of λ). Then a solution with one wrong end condition can be *rescaled* into one with the right
    condition: solve with f″(0) = 1, see f′(∞) = 2.085 instead of 1, and multiply the height by λ² = 1/2.085. (Töpfer 1912.)",
    code="import numpy as np\nfrom scipy.integrate import solve_ivp\nrhs = lambda e, y: [y[1], y[2], -0.5*y[0]*y[2]]      # f''' = -f f''/2\ns = solve_ivp(rhs, [0, 12], [0, 0, 1], rtol=1e-10)     # f''(0) = 1 guess\nprint(s.y[1, -1])                                    # 2.0854: f'(infinity) with the guess\nprint(s.y[1, -1]**-1.5)                              # 0.3321: f''(0) after rescaling, lambda^3 with lambda^2 f'(inf) = 1")` (**P207**)
21. `nb.primer("shooting versus boundary-value solving", "A *boundary-value problem* (conditions at two ends) can be solved by
    *shooting*: guess the missing initial slope, integrate, and adjust the guess until the far end is right (`brentq` on the
    mismatch f′(η_max) − 1), or by `solve_bvp`, which treats the whole interval at once. For Blasius shooting is easy because
    the far value depends smoothly on the guess (here ∝ s^{2/3}); near the separation member of C05 (f″(0) → 0, two solutions
    merging) the root becomes ill-conditioned, so we shoot only for n ≥ −0.05 and use `solve_bvp` with continuation elsewhere.",
    code="import numpy as np
from scipy.integrate import solve_ivp
rhs = lambda e, y: [y[1], y[2], -0.5*y[0]*y[2]]
for s2 in (0.3320, 0.3321):                            # two guesses for f''(0)
    print(s2, solve_ivp(rhs, [0, 10], [0, 0, s2], rtol=1e-12).y[1, -1])   # f'(10) = 0.99988 and 1.00009: smooth, so brentq converges fast")` (**P208**)
22. `nb.note` — **N31 [B], N32 [B]** "**Conditions:** $f=0$ and $f'=0$ at η = 0 (9.28), $f'\\to1$ as η → ∞ (9.29). 'Infinity' is truncated
    at η_max; the solution should not care (truncation study below: η_max = 8, 12, 16 change f″(0) by < 1e-9)." equations (9.28), (9.29).
23. `nb.derivation("D06", …)` — Part F D06 (10 steps), ref "9.32".
24. `nb.note` — **N35 [B]** "**δ₉₉ (9.30):** the 99 % height is $\\delta_{99}=4.91\\sqrt{\\nu x/U}$, i.e. $\\delta_{99}/x=4.91/Re_x^{1/2}$.
    > ⚠️ The book prints 4.93 (read off a figure); the root of f′ = 0.99 is η₉₉ = 4.910, 0.4 % smaller — `BL.blasius_delta99(…, printed=True)`
    gives the printed one." equation $\\delta_{99}=4.93\\sqrt{\\nu x/U}$ (book), ref "9.30".
25. `nb.note` — **N36 [B], N37 [B], N38 [B]** "$\\delta^*=1.721\\sqrt{\\nu x/U}$, $\\theta=0.6641\\sqrt{\\nu x/U}=2f''(0)\\sqrt{\\nu x/U}$; wall
    shear $\\tau_0=\\mu\\big(\\frac{\\partial u}{\\partial y}\\big)_0=0.332\\rho U^2/\\sqrt{Re_x}$ (9.31), which blows up like x^{−1/2} at the leading
    edge (integrable); skin-friction coefficient $C_f\\equiv\\frac{\\tau_0}{\\frac12\\rho U^2}=\\frac{0.664}{\\sqrt{Re_x}}$ (9.32)." equations (9.31), (9.32).
26. `nb.note` — **N39 [B], N40 [B]** "**Drag per unit width (one side):** $F_D=\\int_0^L\\tau_0dx=0.664\\,\\rho U^2L/\\sqrt{Re_L}$, growing like
    $U^{3/2}$ (Stokes $U$, bluff bodies $U^2$); $C_D=\\frac{F_D}{\\frac12\\rho U^2L}=\\frac{1.33}{\\sqrt{Re_L}}$ (9.33), twice the local
    $C_f(L)$. > ⚠️ **One side only:** a plate with two wetted faces has twice this drag (`BL.blasius_drag(…, sides=2)`)." equation (9.33).
27. `nb.note` — **N33 [B], N34 [B]** "**Far field and cross-flow.** The approach to the free stream is Gaussian: with f ≈ η − δ*
    the linearised equation g″ + ½(η − δ*)g′ = 0 for g = f′ − 1 gives $f'-1\\approx-A\\sqrt\\pi\\,\\mathrm{erfc}\\frac{\\eta-\\delta^*}2$ with A = 0.234 (within 0.3 % of the solved profile for 4 ≤ η ≤ 8); far out this is a Gaussian, which the book quotes as (1/η)e^{−η²/4} — the shift by δ* is what makes the numbers fit. The cross-flow
    $\\frac vU=\\frac1{2\\sqrt{Re_x}}(\\eta f'-f)\\to\\frac{0.860}{\\sqrt{Re_x}}$ lifts the outer streamlines (= dδ*/dx, C02)."
28. `nb.worked_example("a 1 m plate in air", "U = 1 m/s, ν = 1.5 × 10⁻⁵: Re_x(1 m) = 6.67 × 10⁴, √Re = 258.2, √(νx/U) = 3.873 mm.
    δ₉₉ = 4.910 × 3.873 = 19.02 mm; δ* = 1.7208 × 3.873 = 6.665 mm; θ = 0.6641 × 3.873 = 2.572 mm; τ₀ = 0.3321 × 1.2 × 1/258.2 =
    1.543 × 10⁻³ Pa; C_f = 0.6641/258.2 = 2.572 × 10⁻³; F_D (L = 1 m, one side) = τ₀ integrated = 2 τ₀(L) L = 3.086 × 10⁻³ N/m;
    C_D = 1.328/258.2 = 5.14 × 10⁻³. Doubling U raises F_D by 2^{3/2} = 2.83.")`
29. `nb.code` — constants and sweeps: `bc = BL.blasius_constants()`; print; `for em in (8, 12, 16): print(em, BL.falkner_skan(0.0, eta_max=em)["fpp0"])`;
    `t = BL.falkner_skan(0.0, method="toepfer")["fpp0"]`; `b = BL.falkner_skan(0.0, method="bvp")["fpp0"]`;
    `assert abs(t - b) < 1e-8`. *expect:* fpp0 = 0.3320573362, eta99 = 4.909990, delta_star = 1.720788, theta = 0.664115, H = 2.591100,
    v_inf = 0.860394, cf_coeff = 0.664115, cd_coeff = 1.328229; f″(0) identical to 1e-9 for the three η_max; Töpfer − bvp < 1e-9;
    profile f′(η) at η = 1, 2, 3, 4, 5, 6 = 0.3298, 0.6298, 0.8460, 0.9555, 0.9915, 0.9990; (v/U)√Re_x = ½(ηf′ − f) at η = 2, 5, 7 =
    0.3048, 0.8372, 0.8601. *explain:* 1. `blasius_constants` solves (9.27) once by the Töpfer IVP and integrates the profile;
    2. the truncation loop shows η_max does not matter; 3. the two independent routes agree.
30. `nb.code` — the plate in air: `x = np.array([0.1, 0.5, 1.0])`; `print(BL.blasius_delta99(x, 1.0, NU_AIR), BL.blasius_delta_star(x, 1.0, NU_AIR),
    BL.blasius_theta(x, 1.0, NU_AIR), BL.blasius_wall_shear(x, 1.0, RHO_AIR, NU_AIR))`; `print(BL.blasius_skin_friction(1.0/NU_AIR),
    BL.blasius_drag(1.0, 1.0, RHO_AIR, NU_AIR), BL.blasius_drag_coefficient(1.0/NU_AIR))`; `print(BL.blasius_delta99(1.0, 1.0, NU_AIR, printed=True) /
    BL.blasius_delta99(1.0, 1.0, NU_AIR))`. *expect:* δ₉₉(1 m) = 19.02 mm (6.01 mm at 0.1 m); τ₀(1 m) = 1.543e-3 Pa; C_f = 2.572e-3;
    F_D = 3.086e-3 N/m; C_D = 5.14e-3; printed/true = 1.004. *explain:* the closed forms are the computed constants times powers of x.
31. `nb.check_agree` — **from scratch (curation §7):** an RK4 loop written by hand (step 0.01, η to 10) for f‴ = −½ff″ from f″(0) = 1, then
    the Töpfer rescale; `assert np.isclose(fpp0_mine, 0.3320573362, rtol=1e-6)`; `assert np.allclose(fp_mine_rescaled(np.arange(0, 6.1, 0.5)),
    bc_profile, atol=1e-6)`; identity checks θ = 2f″(0) and δ* = lim(η − f) to 1e-6.
32. `nb.figure` — **Figs. 9.5 and 9.6 remade** (9 × 3.4 in, three panels): (a) f′(η) = u/U against η ∈ [0, 7] (blue) with η₉₉ = 4.91 (muted
    dashed), δ*-arrow 1.721 (amber) and θ-arrow 0.664 (orange) on the abscissa, the zero curvature at the wall (inflection) noted, the
    Gaussian tail (dotted); (b) (v/U)√Re_x = ½(ηf′ − f) against η ∈ [0, 6] rising to 0.860 (dashed asymptote); (c) log₁₀|1 − f′| against
    η with the erfc form of N33 (purple dashed) — a downward parabola in η. Title "The Blasius profile". *see:* "an S-shaped
    profile flat at the wall, and a cross-flow that saturates"; *read:* "at η = 2 the speed is 63 % of U; the layer edge (99 %) sits at
    4.91; the slope at the wall is 0.332"; *change:* "…the truncation η_max were 5: the tail is cut and f″(0) would be off by 1e-3".
33. `nb.figure` — **the collapse** (8 × 3.4 in, two panels; N42): (a) dimensional u(y) [m/s] at x = 0.1, 0.5, 1 m for air at 1 m/s (blue, darker
    with x) — three different curves; (b) the same data against η = y√(U/νx) — one curve, with the residual `similarity_collapse_error("blasius")` < 1e-12
    in the title. *see* "three stretched curves collapse into one"; *read* "they differ only by the √x stretch in y (×1, ×2.2, ×3.2)"; *change*
    "…water: every y shrinks by √15 but panel (b) is unchanged".
34. `nb.figure` — **Blasius versus the temporal layer (N41)** (bar chart): C_f√Re_x = 0.664 (local Blasius, teal), 1.128 (ch08's Stokes-layer Galilean
    map `LAM.temporal_bl_wall_stress(...)["Cf_coefficient"]`, rose), 1.328 (plate-averaged C_D√Re_L, orange). *see* "three bars"; *read* "the
    Galilean map t = x/U overstates the wall stress by 70 %: in the plate layer the fluid near the wall moves slower than U and is
    advected less"; *change* "…a plate that was sucked (Exercise 9.26): the layer thins and the first bar rises".
35. `nb.animation` — **A2** (video, 60 frames, FAST 30): stations x = 0.1 → 2 m sweep; left the dimensional profile stretching in y, right the
    rescaled f′(η) fixed with a dot riding; δ₉₉ ticks. *see* "the left curve stretches while the right one never moves"; *read* "y scales with √x,
    u/U is a function of η only"; *change* "…U doubled: everything shrinks by √2 in y; the right panel is unchanged".
36. `nb.plotly` — **F2** `slider_figure` over x (≤ 25 values, 0.05–2 m): dimensional and rescaled Blasius profiles with the δ₉₉ marker.
37. `nb.explainer("blasius_similarity_collapse", heading="Why does one curve describe the whole plate?", why="Toggling raw ↔ rescaled makes the
    collapse happen; the truncation and shooting sliders show what 'f′ → 1 at infinity' demands; a static plot shows three profiles, not the
    mechanism.", tries=["Drag the three station sliders and watch raw profiles stretch but the rescaled ones stay together.", "Use the shooting
    slider: set f″(0) 0.30 or 0.36 and read the status (never reaches 1 / overshoots).", "Click a point on the profile: the inspector shows η, f′, u and y.",
    "Open Derivation D05: the picture highlights where the two ηf′f″ terms cancel."])`
38. `nb.md` — **⚠️ Common confusion (leading-edge):** "τ₀ ∝ x^{−1/2} is infinite at x = 0: the layer equations fail where Re_x ≲ 1 (N17), but the
    singularity is integrable, so the drag (9.33) is still right for Re_L ≳ 10³."
39. `nb.md` — **What would change if…** "…the free stream accelerates or decelerates? δ(x) then carries U_e(x) and the ODE gains two terms: C05."

### A.4 §9.4 Falkner–Skan Similarity Solutions — C05 (+N43–N47, N48 · D07 · P209 · A3 · F3 · live · E3)
1. `nb.section("9.4", "Falkner–Skan Similarity Solutions", intro="**What is this section about?** Blasius has a constant outer speed. If the outer
   speed is a power of the distance, U_e = a xⁿ, the layer is still similar — and the single family of ODEs contains the flat plate (n = 0),
   the flow towards a wall (n = 1) and, at n = −0.0904, the last attached profile before the wall shear vanishes. It is our first laboratory
   for pressure-gradient effects and for separation.")`
2. `nb.core("C05", "Falkner–Skan equation $f'''+\\frac{n+1}2ff''-nf'^2+n=0$ (9.36): wedge flows, Blasius at n = 0, stagnation flow at n = 1",
   question="What happens to the layer when the outer flow accelerates or decelerates as a power of x?")`
3. `nb.md` — **The problem in plain words:** "Wind hits a roof ridge, a bow, a wing's leading edge: the outer speed changes along the wall. If it
   rises as x^n, the pressure falls along the wall (favourable), the layer is squeezed and its wall shear is large; if it falls (n < 0) the
   pressure rises along the wall (adverse), the layer thickens and — for n = −0.0904 — the wall shear reaches zero, which is separation."
4. `nb.md` — **The idea** (a small table, the ch06 corner flow $Az^n$ is the ideal flow past a wedge of half-angle πn/(n+1)):
   | n | β = 2n/(n+1) | outer flow | pressure | profile |
   |---|---|---|---|---|
   | 4 | 1.60 | strongly accelerating | strongly favourable | very full |
   | 1 | 1.00 | stagnation point (wall ⊥ stream) | favourable | δ constant |
   | 1/3 | 0.50 | wedge | favourable | full |
   | 0 | 0 | flat plate (Blasius) | none | inflection at the wall |
   | −0.05 | −0.105 | slightly decelerating | adverse | inflection at η ≈ 1.65 |
   | −0.0904 | −0.199 | deceleration | adverse, strong | f″(0) = 0: separation |
5. `nb.note` — **N43 [B], N44 [B]** "**Ansatz and pressure gradient:** $\\psi=\\sqrt{\\nu xU_e}\\,f(\\eta)$, $\\eta=\\frac yx\\sqrt{Re_x}=y\\sqrt{\\frac a\\nu}x^{(n-1)/2}$ (9.34) with
   Re_x = a xⁿ⁺¹/ν; the outer flow imposes $-\\frac{dp}{dx}=U_e\\frac{dU_e}{dx}=na^2x^{2n-1}$ (9.35) (n > 0 favourable, n < 0 adverse)." equations (9.34), (9.35).
6. `nb.note` — **N45 [B]** "**Thickness:** $\\delta(x)=\\sqrt{\\nu x/U_e}=\\sqrt{\\nu x^{1-n}/a}$ grows for n < 1 and is *constant* at n = 1 (a layer of fixed thickness
   under a stagnation point). Numbers below."
7. `nb.primer("continuation in a parameter and a fold (saddle-node)", "To solve a family f(η; n), start at an easy n (Blasius) and move n in small
   steps, using each solution as the guess for the next: *continuation*. A solution branch can *turn back* at a fold: two solutions merge
   and vanish (a saddle-node). Beyond the fold `solve_bvp` fails — not a bug, but the mathematics saying there is no attached solution.",
   code="import numpy as np\nn = np.linspace(-0.2, 0.2, 5)\nprint(n + 0.01)              # the branch x = +sqrt(n + 0.01) exists only where this is >= 0 (a fold at n = -0.01)\nprint(np.sqrt(np.maximum(n + 0.01, 0)))  # the upper branch exists only where n + 0.01 >= 0")` (**P209**)
8. `nb.derivation("D07", …)` — Part F D07 (12 steps), ref "9.36". (Contains N43's algebra; ends with the boundary conditions (9.28)–(9.29) and
   the special cases.)
9. `nb.note` — **N46 [B]** "**The family (Fig. 9.7 remade).** The shear f″(0) rises monotonically with n; the sign of the wall curvature is the sign of n: setting η = 0
   in (9.36) (f = f′ = 0) gives $f'''(0)=-n$, so $(\\partial^2u/\\partial y^2)_{wall}\\propto-n$: negative (no inflection) for n > 0, zero at n = 0, positive for n < 0 (an
   inflection appears at finite y, C08)." equation $f'''(0)=-n$.
10. `nb.note` — **N47 [B]** "**The separation member.** As n decreases from 0 the wall shear f″(0) falls; it reaches zero at $n=-0.0904$ ($\\beta=-0.1988$) — the fold of the
    attached branch. > ⚠️ **Common confusion (slip):** the book says solutions 'exist for n < −0.0904 with reverse flow'. With f′ → 1 there is *no* solution below the fold;
    the reversed-flow profiles (f″(0) < 0) form the second branch for −0.0904 < n < 0 (Stewartson 1954)." equation $f''(0)=0$ at n = -0.0904.
11. `nb.worked_example("stagnation flow, a = 10 s⁻¹, air", "U_e = 10x, so U_e(0.1 m) = 1 m/s, dU_e/dx = 10 s⁻¹; dp/dx = −ρU_eU_e′ = −1.2 × 1 × 10 = −12 Pa/m
    (favourable). δ = √(ν/a) = √(1.5 × 10⁻⁶) = 1.225 mm at every x. f″(0) = 1.2326 (from the table below): τ₀ = μU_e f″(0)/δ = 1.8 × 10⁻⁵ × 1 ×
    1.2326/1.225 × 10⁻³ = 1.81 × 10⁻² Pa — 3.7 times the flat-plate value at the same x and U (4.9 × 10⁻³ Pa at x = 0.1 m); and τ₀ ∝ U_e at constant δ, so
    τ₀(0.2 m) = 3.62 × 10⁻² Pa.")` (the builder prints both from `BL.falkner_skan_fields` and `BL.blasius_wall_shear`).
12. `nb.code` — `for n in (4, 1, 1/3, 1/9, 0, -0.05, -0.0654, -0.0904): st = BL.falkner_skan_state(n); print(n, st["fpp0"], st["H"], st["inflection_eta"])`;
    `print(BL.falkner_skan_separation())`; `print(BL.falkner_skan(-0.095)["success"], BL.falkner_skan(-0.05, branch="reversed")["fpp0"])`;
    `sol = BL.falkner_skan(1.0)`; `assert np.allclose(sol["fppp"][0], -1.0, atol=1e-6)` (f‴(0) = −n). *expect:* f″(0) = 2.4057, 1.2326, 0.7574, 0.5118, 0.3321,
    0.2135, 0.1640, 0.0048; H = 2.172, 2.216, 2.297, 2.410, 2.591, 2.818, 2.963, 3.970; inflection η = None for n ≥ 0 (at the wall for n = 0), 1.650 for n = −0.05;
    m_sep = −0.09043, β_sep = −0.19884; `success` False for n = −0.095; the reversed branch has f″(0) < 0. *explain:* 1. `falkner_skan_state` solves (9.36) by
    `solve_bvp` with continuation and integrates the thicknesses; 2. the fold is located by parametrising with f″(0); 3. below the fold nothing attached exists.
13. `nb.check_agree` — **from scratch (curation §7):** `brentq` on the initial slope for n = 1/3 and n = 1 in η_max = 8 with `solve_ivp` (rtol 1e-11): the root of
    f′(8; s) − 1 = 0 bracketed in [0.2, 3]; `assert np.isclose(s13, 0.757448, atol=1e-5)`, `np.isclose(s1, 1.232588, atol=1e-5)`.
14. `nb.figure` — **Fig. 9.7 remade and the fold** (9 × 3.4 in, two panels): (a) f′ against the scaled variable ½√(n+1)·η of the book's figure, for n =
    4, 1, 1/3, 1/9, 0, −0.0654, −0.0904 (viridis, labelled), the inflection points as dots, the wall tangent (slope f″(0)) as a short black segment; (b) f″(0) against n
    ∈ [−0.0904, 1]: the attached branch (blue) ending at the fold (red dot at −0.09043, f″(0) = 0) and the second branch (amber dashed, f″(0) < 0) continuing back
    toward n = 0. Title "One dial n: from a full profile to zero wall shear". *see* "profiles fanning from steep to lazy; a curve that turns back at the red dot";
    *read* "steeper wall slope = larger friction; the last blue curve has zero slope at the wall"; *change* "…n = −0.1: nothing on the blue branch; only if the profile may
    have reverse flow (amber) and then f′ → 1 fails".
15. `nb.animation` — **A3** (frames player, stops at n = 0 and n = −0.0904): n from 4 → 0 → −0.0904; the profile changes, the wall tangent rotates to zero slope. Notes:
    *see* "the tangent at the wall lying down"; *read* "the last frame: zero slope = zero shear"; *change* "…run back: the layer refills".
16. `nb.plotly` — **F3** `slider_figure` (≤ 30 values of n between −0.0904 and 4, denser near the fold): f′(η) profile plus the f″(0) marker. `nb.live` — free choice of n and
    Re_x: profile, δ(x) and τ₀ (kernel-only; the page shows a note).
17. `nb.explainer("falkner_skan_family", heading="How can favourable and adverse gradients be one family?", why="The profile, its curvature at the wall and the shear f″(0)
    change continuously with n; the inflection slides in from infinity and the wall shear reaches zero at exactly one n — seven static curves cannot show the approach.",
    tries=["Drag n from 1 to 0 to −0.05: where does the inflection appear?", "Press the 'separation' preset: what is f″(0) and the status?", "Switch the scaled variable to the book's
    ½√(n+1)η and compare curves.", "Compare the terms bars: u_yy and the pressure gradient are equal at the wall."])`
18. `nb.pointer` — **N48 [C]** "Real flows are rarely of power-law form, which is why the approximate methods of the next two sections exist; the full numerical solution is Ch. 10."
19. `nb.md` — **What would change if…** "…we do not want to solve an ODE at all? The momentum integral (C06) removes y from the equation, at the price of guessing the profile."

### A.5 §9.5 Von Karman Momentum Integral Equation — C06 (+N49–N55 · D08)
1. `nb.section("9.5", "Von Karman Momentum Integral Equation", intro="**What is this section about?** Instead of solving for the whole profile, integrate the
   momentum equation across the layer. What remains is one ordinary differential equation in x linking the three quantities we already met — the momentum
   thickness θ, the displacement thickness δ* and the wall stress τ₀ — for any pressure gradient, laminar or (time-averaged) turbulent. It has three
   unknowns and one equation, so we need a closure; Thwaites' choice is the next section.")`
2. `nb.core("C06", "The von Kármán momentum integral equation $\\frac1\\rho\\tau_0=\\frac d{dx}[U_e^2\\theta]+U_e\\delta^*\\frac{dU_e}{dx}$ (9.43)", question="Can we get
   the wall shear from an ordinary differential equation instead of the PDE?")`
3. `nb.md` — **The problem in plain words:** "An engineer wants τ₀(x) along a wing or a duct for a given pressure distribution and has an afternoon, not a supercomputer. The
   layer equation is a PDE in x and y; if we integrate it across the layer, y disappears and all we need is the size of the profile (θ, δ*), not its shape."
4. `nb.md` — **The idea:** "**momentum balance on a slab** of the layer: what leaves minus what enters = wall friction + pressure" — the integral form of Ch. 4's control volume, done on
   the differential equation (that is why it holds for turbulent time-averages too):
   ```
   d/dx (momentum flux U_e²θ)  +  (pressure-gradient term U_e δ* U_e′)  =  wall friction τ₀/ρ
   ```
5. `nb.note` — **N49 [B]** "**Start:** with the pressure gradient replaced by (9.11), the layer equation reads $u\\frac{\\partial u}{\\partial x}+v\\frac{\\partial u}{\\partial y}=U_e\\frac{dU_e}{dx}
   +\\frac1\\rho\\frac{\\partial\\tau}{\\partial y}$ (9.37) with τ = μ∂u/∂y." equation (9.37).
6. `nb.derivation("D08", …)` — Part F D08 (11 steps), ref "9.43". (Displays N50 (9.38), N51 (9.39), N52 (9.40), N53 (9.41), N54 (9.42) at steps 2, 3, 4, 5–7, 8; result (9.43).)
7. `nb.note` — **N55 [B]** "**Closure problem.** (9.43) has three unknowns (θ, δ*, τ₀) and one equation. Pohlhausen (1921) assumes a profile shape with one free parameter (a cubic in y/δ); Thwaites
   (1949, C07) assumes instead that shear and shape depend on a single number λ."
8. `nb.worked_example("does Blasius satisfy (9.43)?", "U_e = U constant: $\\tau_0/\\rho=U^2\\,d\\theta/dx$. With θ = 0.6641√(νx/U): dθ/dx = θ/(2x); so U²·0.6641·½√(ν/(Ux)) = 0.3321U²/√Re_x and
   τ₀/ρ = 0.3321U²/√Re_x from (9.31) ✓ — the momentum integral is exact, so it must hold for every exact solution. At x = 1 m in air: τ₀/ρ = 1.286 × 10⁻³ m²/s².")`
9. `nb.code` — residual test: `for n in (0.0, 1/3, -0.05): st = BL.falkner_skan_state(n); x = np.linspace(0.05, 1, 40); Ue = a*x**n; θ = st["I_theta"]*np.sqrt(NU_AIR*x/Ue); ds = st["I_delta"]*...;
   tau0 = RHO_AIR*NU_AIR*Ue*st["fpp0"]/np.sqrt(NU_AIR*x/Ue)` (with ρ ν = μ); `r = BL.momentum_integral_residual(x, Ue, θ, ds, tau0, RHO_AIR)`; print `abs(r).max()/abs(tau0).max()`.
   *expect:* < 1e-7 for each n. *explain:* 1. build θ, δ*, τ₀ of an exact Falkner–Skan layer from the solved f; 2. the residual is the difference of the two sides of (9.43) using 4th-order finite differences.
10. `nb.check_agree` — **from scratch (curation §7):** θ(x), δ*(x) from Blasius, `np.gradient` for d(U²θ)/dx, `assert np.allclose(rho*dθdx*U**2, tau0, rtol=1e-4)`; plus the same with
    `BL.momentum_integral_residual`.
11. `nb.code` — optional Kármán–Pohlhausen: `for p in ("cubic", "sine"): print(BL.karman_pohlhausen(BL.outer_flow("flat", U=1.0).Ue, np.array([1.0]), NU_AIR, profile=p))`. *expect (units √(νx/U)):*
    cubic δ = 4.641, θ = 0.6465 (−2.7 %), δ* = 1.740 (+1.1 %), τ₀ coefficient 0.3232 (−2.7 %); sine δ = 4.795, θ = 0.6551 (−1.4 %), τ₀ 0.3276 (−1.3 %). *explain:* an assumed shape closes (9.43) and gets within 3 %
    of Blasius — the idea Thwaites improves on.
12. `nb.figure` — **the momentum budget** (9 × 3.4 in, two panels): (a) stacked bars along x for a decelerating layer n = −0.05 (a = 1): τ₀/ρ (rose) = d(U_e²θ)/dx (teal) + U_eδ*U_e′ (orange, negative); the dashed black
    line is the sum from the code — bars and line coincide; (b) profiles: exact Blasius (black), cubic (blue dashed), sine (teal dashed) with their θ and τ₀ errors. Title "(9.43) is exact; the closure is the only
    approximation". *see:* "bar heights adding to the black line; three profiles close to each other"; *read:* "for the adverse case the pressure-gradient term subtracts, so τ₀ falls faster than dθ/dx";
    *change:* "…n > 0: the orange term turns positive and adds to the friction".
13. `nb.md` — **What would change if…** "…we choose the closure from the exact Falkner–Skan family instead of a profile? Then θ, δ* and τ₀ depend on the single number λ and (9.43) integrates in closed form: Thwaites, next."

---
### A.6 §9.6 Thwaites' Method — C07 (+N56–N67 · D09 D10 D11 · P210 P211 · A4 · F4 · E4)
1. `nb.section("9.6", "Thwaites' Method", intro="**What is this section about?** The momentum integral (9.43) has three unknowns. Thwaites noticed that
   for the exact solutions the wall shear and the shape of the profile both depend on a single dimensionless number, the pressure-gradient parameter
   λ = (θ²/ν) dU_e/dx. With that closure (9.43) becomes a first-order linear equation whose solution is one integral of U_e⁵. Give the outer speed of any
   body and you get θ(x), the wall shear and the separation point in a minute of arithmetic. We build the closure from our own Falkner–Skan solutions,
   so the notebook does not need the book's table.")`
2. `nb.core("C07", "Thwaites' method: $\\frac{\\theta^2U_e^6}{\\nu}=0.45\\int_0^xU_e^5dx'+\\frac{\\theta_0^2U_0^6}{\\nu}$ (9.50) and the separation test",
   question="Can one integral of the outer speed predict wall shear and separation without solving a PDE?")`
3. `nb.md` — **The problem in plain words:** "You have a pressure distribution on a duct wall or an airfoil — U_e(x) from ideal-flow theory — and want to know: how thick does the layer get, how
   much friction does it cause, and does it separate, and where? The answer in one pass: integrate U_e⁵, read off λ, compare with −0.09."
4. `nb.md` — **The idea** (ASCII flow chart):
   ```
   U_e(x) ──► ∫U_e⁵dx ──► θ(x) ──► λ = (θ²/ν)U_e′ ──► l(λ), H(λ) ──► τ₀ = μ(U_e/θ) l(λ),  δ* = Hθ
                                          │
                                 λ falls to −0.09 (l = 0)  ⇒  τ₀ = 0  ⇒  separation predicted
   ```
5. `nb.note` — **N56 [B], N57 [B], N58 [B]** "**Definitions (Holstein–Bohlen):** the pressure-gradient parameter $\\lambda\\equiv\\frac{\\theta^2}{\\nu}\\frac{dU_e}{dx}$ (9.44)
   (λ < 0: decelerating, adverse); the shear correlation $\\tau_0\\equiv\\mu\\frac{U_e}{\\theta}\\,l(\\lambda)$ (9.45); the shape factor $\\frac{\\delta^*}{\\theta}\\equiv H(\\lambda)$ (9.46).
   The claim is that l and H depend on λ alone. > ⚠️ **Names:** λ here has nothing to do with wavelength (Ch. 7); L(λ) below is a correlation, L elsewhere a length; Thwaites' own m is −λ; θ₀ is the initial
   momentum thickness." equations (9.44), (9.45), (9.46).
6. `nb.note` — **N59 [B]** "**The closure, from our own solutions (the book's Table 9.1 is not reproduced).** For a Falkner–Skan layer everything depends on n only: with I_δ = ∫(1 − f′)dη and
   I_θ = ∫f′(1 − f′)dη, λ = n I_θ², l = I_θ f″(0), H = I_δ/I_θ (D10). Solving (9.36) for 60 values of n gives l(λ) and H(λ) for −0.0681 ≤ λ ≤ 0.1065 (`closure='falkner_skan'`); the attached family ends at
   λ = −0.0681 (below it the exact closure returns l = 0: separated) and above 0.1065 a fit scaled to match continues it. The book's own criterion is reproduced by `closure='white'`, the fit l ≈ (λ + 0.09)^0.62 (ANSYS lesson handout), which reaches zero at λ = −0.09. Exact separation of the family: λ = −0.0681 (l = 0.004, H = 3.97 at n = −0.0904);
   Thwaites' fit puts l = 0 at λ = −0.09 — two criteria, both reported."
7. `nb.code` — closure table: `tab = BL.thwaites_closure_table(n_points=60)`; print l, H at λ = −0.06, 0, 0.0855: `BL.thwaites_l(lam)`, `BL.thwaites_H(lam)`. *expect:* l(0) = 0.2205, H(0) = 2.591; at
   λ = 0.0855 (n = 1): l = 0.3603, H = 2.216; at λ = −0.0675: l = 0.0163, H = 3.813; `BL.thwaites_l(-0.0681, closure="falkner_skan")` ≈ 0.004. *explain:* 1. each row is one solved Falkner–Skan
   profile; 2. PCHIP interpolation between rows; 3. `holstein_bohlen` is (9.44).
8. `nb.figure` — **l(λ), H(λ)** (8 × 3.4 in, two panels): l and H against λ (blue dots at the 60 Falkner–Skan members; 'white' fit as a dashed curve; vertical lines at λ = −0.0681 and −0.09; Blasius
   λ = 0 and Hiemenz λ = 0.0855 labelled). Title "The closure: shear and shape depend on λ alone". *see:* "l rises and H falls as λ increases"; *read:* "at λ = 0 (Blasius) l = 0.22, H = 2.59; near −0.07 the
   shear is almost zero and H ≈ 4"; *change:* "…a stronger favourable gradient: H tends to 2.1 and l saturates".
9. `nb.primer("first-order linear ODE and the integrating factor", "An equation y′ + p(x) y = q(x) is solved by multiplying by a *factor* μ(x) chosen so that the left side is one derivative:
   with μ′ = pμ (μ = e^{∫p}) we get (μ y)′ = μ q, which integrates directly. Here p = 6U_e′/U_e, so μ = U_e⁶.", code="import sympy as sp\nx = sp.symbols('x', positive=True)\nz = sp.Function('z')\nprint(sp.dsolve(sp.Eq(z(x).diff(x) + 6/x*z(x), 1), z(x)))   # z = x/7 + C/x**6: the factor x**6 makes (x**6 z)' = x**6")` (**P210**)
10. `nb.derivation("D09", …)` — Part F D09 (10 steps), ref "9.50". (Displays N60 (9.47), N61 l = (2 + H)λ + (U_e/2)(θ²/ν)′, N62 (9.48), N64 (9.49), and (9.50).)
11. `nb.note` — **N63 [B]** "**Fig. 9.8 remade: the fit L ≈ 0.45 − 6.0λ.** On the Falkner–Skan family L(λ) = 2l − 2(2 + H)λ is nearly linear (0.441 at Blasius; 0.818 at λ = −0.0675; ≈ 0 at n = 1); the line
    0.45 − 6λ is Thwaites' fit to several *families* (Falkner–Skan, Schubauer's ellipse, Howarth's linear retardation, Iglisch). It is a fit, not exact."
12. `nb.derivation("D10", …)` — Part F D10 (9 steps), ref "9.48". `nb.plotly` — **F4** `slider_figure` over the family (closure family): L(λ) computed on the FS members against the line 0.45 − 6λ, the slider moving the highlighted member n (30 steps).
13. `nb.note` — **N65 [B]** "**Accuracy.** The book quotes ±3 % (favourable) and ±10 % (adverse) and warns that Thwaites predicts *whether* a layer separates better than *where*. We measure it on
    the exact Falkner–Skan family: θ_Thwaites/θ_exact = √(0.45/((5n + 1)I_θ²)); table below." + `nb.code`: for n in (−0.05, 0, 1/3, 1, 4): `st = BL.falkner_skan_state(n)`; `err = np.sqrt(0.45/((5*n+1)*st["I_theta"]**2)) - 1`;
    print. *expect:* +3.1 % (n = −0.05), +1.0 % (0), −4.2 % (1/3), −6.3 % (1), −7.6 % (4). *explain:* the fit L = 0.45 − 6λ is best near Blasius; for strong acceleration it undershoots θ by up to 8 %.
14. `nb.worked_example("Example 9.1 — Thwaites for a flat plate", "U_e = U constant, so U_e′ = 0, λ = 0, and (9.50) with θ₀ = 0 gives θ² = 0.45νx/U, θ = 0.6708√(νx/U): 1.0 % above the exact Blasius 0.6641.
    δ* = H(0)θ = 2.591 × 0.6708 = 1.738√(νx/U) (+1.0 % over 1.721). τ₀ = μ(U/θ)l(0) gives C_f√Re_x = 2 l(0)/0.6708 = 0.6575 (−1.0 % against 0.6641). One integral of a constant reproduces Blasius to 1 %.")`
15. `nb.code` — Example 9.1: `ex = ch09.example_9_1()`; print. *expect:* theta_coef = 0.6708, +1.0 %; delta_star_coef = 1.738; cf_sqrtRex = 0.6575, −1.0 %. And `BL.thwaites(x, BL.outer_flow("flat", U=1.0).Ue, NU_AIR)` at x = 1 m agrees:
    θ = 2.598 mm (0.6708 × 3.873 mm). *explain:* `example_9_1` is the closed form; `thwaites` is the numerical route; they must agree to 1e-8.
16. `nb.worked_example("Example 9.2 — a straight diffuser", "Area grows A = A₁(1 + x/L) so U_e = U₁/(1 + x/L). With θ₀ = 0: ∫U_e⁵dx = U₁⁵L[1 − (1 + x/L)⁻⁴]/4 and dU_e/dx = −U₁/(L(1 + x/L)²), so
    λ = −(0.45/4)[(1 + x/L)⁴ − 1]: at x/L = 0.05, 0.10, 0.15, 0.20 that is −0.0242, −0.0522, −0.0843, −0.1208. λ = −0.09 at (1 + x/L)⁴ = 1.8, x/L = 0.1583. With an initial θ₀ = 0.1 mm
    (U₁ = 10 m/s, L = 0.5 m, ν = 1.5 × 10⁻⁵) the extra term −0.01333(1 + x/L)⁴ moves the crossing to x/L = 0.1263: a thicker inlet layer separates sooner. With the exact-FS criterion −0.0681 the crossing is at 0.1256 (θ₀ = 0).")`
17. `nb.code` — Example 9.2: `ex2 = ch09.example_9_2()`; `ex2b = ch09.example_9_2(theta0=1e-4)` (with U₁ = 10, L = 0.5, ν = 1.5 × 10⁻⁵ defaults documented); `x = np.linspace(0, 0.3, 601)*0.5`;
    `th = BL.thwaites(x, BL.outer_flow("diffuser", U1=10.0, L=0.5).Ue, NU_AIR)`; print `th["x_sep"]/0.5` for the default (exact-FS closure, λ_sep = −0.0681) and for `closure="white"` (the book's −0.09); `assert np.isclose(th["lam"][60], ex2["lam"](0.05), atol=2e-4)`.
    *expect:* λ(0.05, 0.10, 0.15, 0.20) = −0.02424, −0.05221, −0.08426, −0.12078; x_sep/L = 0.15829 (`closure="white"`, λ = −0.09), 0.12563 (exact-FS closure, λ = −0.0681); with θ₀ = 0.1 mm: 0.1263. *explain:* the closed form and the marching
    agree; separation is where l(λ) → 0.
18. `nb.check_agree` — **from scratch (curation §7):** Thwaites in six lines with `cumulative_trapezoid` of U_e⁵ (θ² = 0.45 ν/U_e⁶ × I, λ = θ²U_e′/ν, l from the linear-L route?) — the hand version uses the closed form
    for the diffuser: `lam_mine = -(0.45/4)*((1 + s)**4 - 1)` and `np.allclose(lam_mine, th["lam"], atol=2e-4)`; and a hand `np.sign` search for the first crossing of −0.09 using
    `np.interp` (P180/P182): `assert np.isclose(xs_mine, 0.1583*0.5, rtol=2e-3)`.
19. `nb.primer("integrals of powers of sine (∫sin⁵ by c = cos φ)", "∫sin⁵φ dφ is done by writing sin⁵φ dφ = (1 − cos²φ)² sin φ dφ and substituting c = cos φ, dc = −sin φ dφ: ∫(1 − c²)²(−dc) = −c + (2/3)c³ − c⁵/5. The
    definite integral from 0 to φ is that minus its value at c = 1: 8/15 − c + (2/3)c³ − c⁵/5.", code="import numpy as np\nfrom scipy.integrate import quad\nphi = np.radians(60); c = np.cos(phi)\nprint(8/15 - c + 2/3*c**3 - c**5/5, quad(lambda p: np.sin(p)**5, 0, phi)[0])   # 0.0958 twice")` (**P211**)
20. `nb.derivation("D11", …)` — Part F D11 (8 steps), ref "9.50". Then `nb.code`: `print(BL.thwaites_cylinder_separation(), BL.thwaites_cylinder_separation(lam_sep=-0.0681))`; `phi = np.radians([30, 60, 82, 90, 100, 103])`;
    `print(BL.thwaites_cylinder_closed_form(phi))`; `assert np.allclose(BL.thwaites_cylinder_closed_form(phi), BL.thwaites_cylinder(np.degrees(phi)), atol=1e-4)`. *expect:* 103.11° and 100.89°; λ(30°, 60°, 82°, 90°, 100°) =
    0.0722, 0.0589, 0.0263, 0.0000, −0.0603; the limit φ → 0 is 0.075. *explain:* Thwaites on the *ideal-flow* speed 2U sin φ predicts separation near 103°, but real cylinders separate near 82° (the book's rounded value, C09).
21. `nb.figure` — **three Thwaites pictures** (9 × 3.4 in, three panels): (a) Fig. 9.8 remade: L(λ) on the FS members (blue dots) with the line 0.45 − 6λ (black) and the residual in the title; (b) the diffuser (Example 9.2) λ(x/L) for θ₀ = 0 and 0.1 mm with the
    −0.09 and −0.0681 lines and the crossings marked; (c) λ(φ) on the cylinder with the crossings at 103.1° / 100.9° and the ideal-flow speed sketched. Title "One integral predicts where the wall shear dies". *see:* "λ falling through the
    lines"; *read:* "a bigger θ₀ moves the crossing upstream; in (c) λ crosses zero at 90° where the ideal speed peaks"; *change:* "…a longer diffuser (smaller angle) is the same curve stretched in x: separation is set by x/L".
22. `nb.animation` — **A4** (frames player, 30 frames FAST 16): marching through the diffuser: the profile thickening (Thwaites θ, H), the marker λ sliding toward −0.09 and τ₀(x) bar shrinking to zero. Notes: *see* "the red marker touches the line"; *read* "τ₀ = 0
    where λ = −0.09"; *change* "…θ₀ larger: the marker starts closer to the line".
23. `nb.explainer("thwaites_marching", heading="Can one integral predict separation?", why="The reader chooses the body's outer flow, the θ, λ, H, τ₀ curves update along x, and the separation point jumps as the diffuser angle, θ₀ or closure changes;
    static curves for one body cannot show that.", tries=["Choose 'diffuser' and raise θ₀: how far upstream does separation move?", "Switch the closure to 'white' and compare the −0.09 and −0.0681 criteria.", "Choose 'cylinder' and read the predicted separation angle; then compare with 82°.",
    "Click on the θ² bars: which term is larger downstream, ∫U_e⁵ or the θ₀ memory?"])`
24. `nb.md` — **What would change if…** "…the outer flow keeps decelerating? λ keeps falling and l(λ) reaches zero: the wall shear vanishes. Next section: what that means physically, and what it does to a body."

### A.7 §9.7 Transition, Pressure Gradients, and Boundary-Layer Separation — C08 (+N70–N73, N74 · D12 · P212 · F7), C09 (+N68 N69 N75 · D13 · F5)
1. `nb.section("9.7", "Transition, Pressure Gradients, and Boundary-Layer Separation", intro="**What is this section about?** Three things can happen to a laminar layer as it runs along a body: it can turn
   turbulent (transition), it can be squeezed and stay attached under a favourable pressure gradient, or, under an adverse gradient, it can stop clinging to the wall (separation). Separation is the boundary layer's way
   of failing: the stream leaves the wall, a wake forms and a drag appears that ideal-flow theory never saw. We first see why an adverse gradient must end in zero wall shear (C08), then what separation does to a
   body's drag (C09).")`
2. `nb.core("C08", "Boundary-layer separation: wall curvature $\\mu u_{yy}|_{wall}=dp/dx$, (9.51)–(9.52), the inflection point and $\\tau_0=0$", question="Why does an adverse pressure gradient make the wall shear vanish?")`
3. `nb.md` — **The problem in plain words:** "Roll a marble up a hill: it slows and, if the hill is steep or long enough, turns back. The fluid nearest the wall has almost no momentum (it is slowed by friction), so a rising pressure along the
   wall — an adverse gradient — reverses it first. When the wall shear reaches zero the layer separates: the stream detaches, a wake forms, the wing stalls."
4. `nb.md` — **The idea** (ASCII, three profiles at the same x):
   ```
   favourable dp/dx<0      zero dp/dx = 0         adverse dp/dx>0           separation
     ▏ ╭─────               ▏ ╭─────               ▏  ╭───╮ inflection      ▏   ╭──►
     ▏╭╯   u_yy(0)<0        ▏╭╯   u_yy(0)=0        ▏ ╭╯   ╰──  u_yy(0)>0     ▏  ╭╯       τ₀ = 0, u_y(0) = 0
     ▏╯   no inflection     ▏╯ (inflection at wall) ▏ ╯                       ▏◄─╯ reverse flow near wall
   ```
5. `nb.primer("inflection point", "An *inflection point* of a curve is where its second derivative changes sign (the curve stops bending one way and starts bending the other). For a velocity profile u(y) it is where u_yy = 0
   and changes sign. Profiles with an inflection are the ones Ch. 11 shows to be unstable to small disturbances (Rayleigh's criterion).", code="import numpy as np\ny = np.linspace(0, 1, 201)\nu = 3*y**2 - 2*y**3 - 0.3*y            # a made-up profile\nuyy = np.gradient(np.gradient(u, y), y)     # second derivative\ni = np.nonzero(np.diff(np.sign(uyy[2:-2])))[0][0] + 2\nprint(y[i])                          # 0.5: u_yy = 6 - 12 y changes sign at y = 1/2")` (**P212**)
6. `nb.note` — **N70 [B]** "**Read (9.9) at the wall.** There u = v = 0, so the two advective terms vanish and $0=-\\frac1\\rho\\frac{dp}{dx}+\\nu\\big(\\frac{\\partial^2u}{\\partial y^2}\\big)_{wall}$, i.e.
   $\\mu\\big(\\frac{\\partial^2u}{\\partial y^2}\\big)_{wall}=\\frac{dp}{dx}$ — the wall curvature *is* the pressure gradient (divided by μ)."
7. `nb.derivation("D12", …)` — Part F D12 (6 steps), ref "9.52". (Displays N71 (9.51) accelerating: $(\\partial^2u/\\partial y^2)_{wall}<0$; N72 (9.52) decelerating: $>0$ hence an inflection point.)
8. `nb.worked_example("sign and size", "Air, adverse dp/dx = +20 Pa/m: u_yy(wall) = 20/1.8 × 10⁻⁵ = 1.1 × 10⁶ (m s)⁻¹ > 0; near the edge u_yy < 0 (the profile rounds off onto U_e), so u_yy must change sign somewhere in between: an
   inflection. For a Falkner–Skan layer f‴(0) = −n: n = −0.05 gives +0.05 (positive, inflection at η = 1.65); n = +1 gives −1 (negative, no inflection); n = 0 gives 0 (the inflection sits at the wall).")`
9. `nb.code` — `print(BL.wall_curvature(20.0, 1.8e-5))`; `for n in (1.0, 0.0, -0.05): sol = BL.falkner_skan(n); print(n, sol["fppp"][0], BL.profile_inflection(sol["eta"], sol["fp"]))`; the diffuser wall shear:
   `th = BL.thwaites(x, BL.outer_flow("diffuser", U1=10.0, L=0.5).Ue, NU_AIR)`; `print(BL.separation_point(x, th["tau0"]))`. *expect:* 1.11e6; (1, −1, None), (0, 0, 0), (−0.05, +0.05, 1.650); x_sep = 0.0628 m (x/L = 0.1256) with the default exact-FS closure, whose l reaches zero at λ = −0.0681; with `closure="white"`
   (the book's criterion λ = −0.09) it is 0.0791 m (x/L = 0.1583). *explain:* 1. `fppp[0]` = f‴(0) = −n exactly; 2. `profile_inflection` finds the sign change of u_yy; 3. `separation_point` interpolates τ₀'s first zero.
10. `nb.note` — **N73 [B]** "**Adverse gradient thickens the layer.** v(y) = −∫₀^y u_x dy′: in a decelerating stream u_x < 0 also *outside* the layer, so there is more fluid to push outward; the layer thickens by diffusion *and* by
    being lifted off the wall (Fig. 9.12 remade below)."
11. `nb.check_agree` — **from scratch (curation §7):** the sign-change search written by hand with `np.sign` and linear interpolation (P180/P182) on the τ₀ array of the diffuser; `assert np.isclose(x_sep_mine, BL.separation_point(x, th["tau0"]), rtol=1e-9)`.
12. `nb.figure` — **Fig. 9.12 and the wall relation** (9 × 3.4 in, three panels): (a) profiles u/U_e against y for n = 1, 0, −0.05 (teal, blue, rose) at the same η range, the wall tangent, inflection dots; (b) u_yy against η with the sign at the wall labelled (negative,
    zero, positive: μ u_yy = dp/dx); (c) τ₀(x) along the diffuser with the region dp/dx > 0 shaded and x_sep marked. Title "Adverse pressure gradient: inflection, then zero shear". *see:* "the rose profile gets a second bend"; *read:* "wall curvature =
    pressure gradient; the inflection appears as soon as the gradient is adverse"; *change:* "…dp/dx stronger: the tangent at the wall lies flat, τ₀ = 0".
13. `nb.plotly` — **F7** `slider_figure`: pressure gradient (n from 1 to −0.0904, 25 steps) with the profile, the wall tangent and the sign of the curvature in the trace name.
14. `nb.pointer` — **N74 [C]** "After separation the boundary-layer equations fail (Goldstein singularity: marching stops at τ₀ = 0; the reversed flow is not thin, the pressure is no longer the ideal one, the flow is unsteady). Chs. 10 and 14."
15. `nb.md` — **What would change if…** "…the layer were turbulent? Its wall-nearest fluid has more momentum, so it survives a stronger adverse gradient: separation is delayed. That is the whole of the drag crisis (C11). First: what does separation do to the drag?"

16. `nb.core("C09", "Streamlining and form drag: separated wake pressure, laminar and turbulent separation angles", question="Why does letting go of the wall create drag, and how does streamlining reduce it?")`
17. `nb.md` — **The problem in plain words:** "Ideal flow past a cylinder is symmetrical front to back: the pressure that pushes on the front is recovered on the rear and the net force is zero (d'Alembert). A real cylinder has a wide wake at low, nearly uniform pressure —
    the rear pressure is not recovered — so the front pushes and the rear does not push back. That missing recovery is the *form drag* (pressure drag). Streamlining is the art of keeping the layer on the wall long enough to recover pressure."
18. `nb.md` — **The idea** (ASCII, C_p around a cylinder against angle φ from the forward stagnation point):
    ```
    C_p   1 ●╲                        ideal (potential) flow: 1 − 4 sin²φ, recovers to 1 at φ = 180°
          0 ──╲──────────╱────         real: follows the ideal curve up to separation φ_s, then stays
         −2    ╲__    __╱              flat at the low wake pressure C_b  ⇒  drag ∝ area between the curves
       0°     φ_s ────────── 180°      later separation (turbulent, 125° vs laminar 82°) ⇒ smaller drag
    ```
19. `nb.note` — **N68 [B]** "**Transition.** A laminar layer becomes unstable beyond a critical Reynolds number (free-stream turbulence, roughness, curvature and pressure gradient move it); on a flat plate, roughly Re_x ≈ 10⁶ within a factor of five.
    The regimes along a plate: leading edge (Re_x ~ 1), laminar similarity region, first instability (Tollmien–Schlichting, Ch. 11), nonlinear breakdown, fully turbulent. > ⚠️ **Common confusion:** the text uses Re_L = UL/ν of the whole plate and the local
    Re_x = 5 × 10⁵ for the start of transition and 10⁶ or 10⁷ elsewhere — `BL.transition_state(Re_x, Re_cr=5e5)` keeps Re_cr an argument." + `nb.code`: `print([BL.transition_state(r) for r in (1e5, 5e5, 2e6, 1e7)])` (*expect:* laminar, transitional, transitional, turbulent).
20. `nb.note` — **N69 [B]** "**Plate drag curve (Fig. 9.11 remade).** Laminar: $C_D=1.33/\\sqrt{Re_L}$ (9.33); turbulent (Prandtl's one-seventh-law form, Ch. 12): 0.074 Re_L^{−1/5}; mixed: the turbulent line minus a patch A/Re_L (A = 1743 at Re_tr = 5 × 10⁵). A turbulent
    layer has larger wall shear (fuller profile) but at high Re it keeps a plate's drag lower than an extended laminar one would be — laminar flow that long is not stable." + `nb.code`: `for r in (1e5, 1e6, 1e7): print(BL.plate_drag_coefficient(r, "laminar"), BL.plate_drag_coefficient(r, "turbulent"), BL.plate_drag_coefficient(r, "mixed"))`.
    *expect:* laminar 4.20e-3, 1.328e-3, 4.20e-4; turbulent (Re_L = 1e6: 4.669e-3; 1e7: 2.946e-3); mixed at 1e6 2.926e-3, at 1e7 2.772e-3 (and equal to laminar at Re_tr = 5e5: 1.878e-3). *(the turbulent correlation must be cited in its docstring — Part C 1.24).*
21. `nb.figure` — **plate regimes and drag curve** (8 × 3.4 in, two panels): (a) a schematic plate with the five regimes (`ch09_drawings.plate_layer`, labelled "schematic, ours"); (b) log–log C_D against Re_L (10⁵–10⁹): laminar line slope −½ (blue), turbulent (rose), the mixed curve (black) leaving the laminar
    line near 5 × 10⁵ and joining the turbulent line by 10⁷. *see:* "two lines and a bridge"; *read:* "at Re_L = 10⁶ a fully laminar plate would have C_D = 1.3 × 10⁻³, a mixed one 2.9 × 10⁻³"; *change:* "…a rougher leading edge trips earlier: the bridge moves left".
22. `nb.derivation("D13", …)` — Part F D13 (7 steps), `ref=""` (ours: the book has no number for it). Result $C_{D,p}=\\sin\\varphi_s\\big(1-\\frac43\\sin^2\\varphi_s-C_b\\big)$.
23. `nb.worked_example("the drag of the separated model by hand", "Take φ_s = 90° and the ideal value of C_p there, C_b = 1 − 4 = −3: C_D = 1·(1 − 4/3 + 3) = 2.667. A realistic wake pressure C_b = −1 gives 1·(1 − 4/3 + 1) = 0.667; C_b = −1.2 and φ_s = 82° gives 0.9903 × (1 − 1.3075 + 1.2) = 0.884
    (sin 82° = 0.9903, sin² = 0.9806); φ_s = 125° with C_b = −0.6: 0.8192 × (1 − 0.8947 + 0.6) = 0.578. And φ_s → 180°, C_b → 1 gives 0: d'Alembert. (The two C_b values are illustrative, ours.)")`
24. `nb.code` — `phi = np.linspace(0, 180, 361)`; `cp_ideal = 1 - 4*np.sin(np.radians(phi))**2`; `cp_sub = BB.separated_cp(phi, 82.0, cp_base=-1.2)`; `cp_sup = BB.separated_cp(phi, 125.0, cp_base=-0.6)`; `for ps, cb in ((82., -1.2), (125., -0.6), (90., None), (179.9, None)): print(ps, cb,
    BB.separated_pressure_drag(ps, cb))`; the Gauss–Legendre cross-check `BB.pressure_drag_from_cp(phi_full, cp_full)`. *expect:* 0.8840, 0.5778, 2.667, 1.4e-8; the integral of the piecewise C_p agrees with the closed form to 1e-12. *explain:* 1. C_p ideal from (6.x) 1 − 4 sin²φ; 2. `separated_cp` freezes the
    pressure at C_b after φ_s; 3. `separated_pressure_drag` is D13's closed form; 4. the quadrature is an independent route.
25. `nb.check_agree` — **from scratch (curation §7):** hand trapezoid of ½∮C_p cos φ dφ over 0…2π (symmetric halves, 4001 points); `assert np.isclose(cd_mine, BB.separated_pressure_drag(82.0, -1.2), rtol=1e-6)`.
26. `nb.figure` — **form drag as an area** (8 × 3.4 in, two panels): (a) C_p(φ) ideal (black dashed), subcritical model 82° (rose, C_b = −1.2) and supercritical 125° (teal, C_b = −0.6) with the region between the curve and the ideal recovery shaded; the caption states "illustrative, ours"; (b) C_D,p against φ_s from 60° to 180° at two fixed C_b with the two
    operating points and the d'Alembert limit at 180°. Title "Later separation, smaller wake, less drag". *see:* "a red wake plateau and a teal one that sits higher"; *read:* "the shaded area between real and ideal is the missing pressure recovery = the form drag"; *change:* "…a lower C_b: bigger drag at the same φ_s (raising the wake suction)".
27. `nb.plotly` — **F5** `slider_figure` over the separation angle φ_s (60°–180°, 25 steps): C_p(φ) and C_D,p in the title.
28. `nb.pointer` — **N75 [C]** "Internal separation in diffusers, elbows and valves is the same mechanism (Example 9.2 is a diffuser); Ch. 12 pipes." + `nb.md` **⚠️ Common confusion** "Streamlined bodies have mostly friction drag, bluff bodies mostly form drag; the same body can change class by delaying separation."
29. `nb.md` — **What would change if…** "…the Reynolds number grew until the layer went turbulent before separating? Then φ_s moves back and the drag falls. That is the subject of the next two blocks, after we meet the wake at lower Re."

### A.8 §9.8 Flow Past a Circular Cylinder — C10 (+N76 N77 N78 N80 · N79 · D14 · P213 P214 · A5 · E6), C11 (+N81 · N82 · E5 — continued in §9.9)
1. `nb.section("9.8", "Flow Past a Circular Cylinder", intro="**What is this section about?** A cylinder in a stream is the standard bluff body. As the Reynolds number rises from below 1 to above 10⁶ the flow passes through a ladder of states: symmetric creeping flow, two steady eddies
   behind the body, a periodic **vortex street**, an irregular wake, and finally the layer itself turns turbulent (the drag crisis). We take the ladder in two blocks: the wake and its street (C10) and the crisis (C11).")`
2. `nb.core("C10", "Cylinder wake regimes and the Kármán vortex street with its stable spacing $b/a=\\frac1\\pi\\cosh^{-1}\\sqrt2=0.2805$", question="Why does a steady stream past a cylinder shed a periodic street — and why has the street a fixed shape?")`
3. `nb.md` — **The problem in plain words:** "Wind sings in power lines, chimneys sway, clouds line up in a zig-zag downstream of an island. Behind a cylinder the steady wake is unstable: it breaks into vortices of alternating sign. Kármán (1912) modelled them as two rows of point vortices and found that only one
   shape of the street survives — the ratio of row separation to spacing is 0.28."
4. `nb.md` — **The idea** (ASCII, the ladder and the street):
   ```
   Re:  <1        4…40             40…200          ≳200 (irregular)     several thousand+      3×10⁵ (crisis)
        creeping  two steady       laminar         irregular vortices,   turbulent wake        layer turbulent,
        symmetric attached eddies  Kármán street   St ≈ 0.2 persists     beyond a few d        wake narrows
   street:   ●  ○  ●  ○   row A (Γ), spacing a          b/a = 0.2805 ⇒ cosh(πb/a) = √2
               ○  ●  ○  ●   row B (−Γ), offset a/2
   ```
5. `nb.note` — **N76 [B], N77 [B], N80 [B]** "**Regimes (the book's rounded thresholds — experimental, not computed):** Re < 1 symmetric creeping flow (vorticity diffuses, then is advected, Ch. 8 Oseen); 4 < Re < 40 two steady eddies whose length grows with Re; Re ≳ 40 the wake
   oscillates and sheds a street; Re < 200 laminar street; above 200 irregular vortices but the frequency stays; several thousand: periodic only near the cylinder, a turbulent wake beyond." + `nb.code`: `for Re in (0.5, 10, 100, 1000, 1e5, 1e6): print(Re, BB.cylinder_flow_regime(Re)["label"])`.
   *expect:* "creeping flow: symmetric, no wake" · "two steady attached eddies" · "laminar Karman street" · "irregular vortices, St stays near 0.2" · "subcritical: laminar separation, wide wake" · "supercritical: turbulent separation, narrow wake". (The table is marked "book's rounded values".)
6. `nb.note` — **N78 [B]** "**Strouhal number** $St=\\Omega d/U_\\infty\\approx0.2$ (4.102) (recalled from Ch. 4): here the book's Ω is the shedding frequency in cycles per second, so f = St U/d — a wire of d = 2 mm in a 10 m/s wind sings at
   f = 0.2 × 10/0.002 = 1000 Hz. > ⚠️ **Common confusion (Ω vs f):** read Ω as an *angular* frequency (rad/s) and the same formula gives St = 2π × 0.2 = 1.26; the angular frequency of the 1000 Hz tone is 6283 rad/s. Say which one your Ω is." + `nb.code`: `print(BB.shedding_frequency(10.0, 0.002))`;
   `print(SIM.strouhal_number(1000.0, 0.002, 10.0), SIM.strouhal_number(6283.19, 0.002, 10.0))`. *expect:* f = 1000 Hz, omega_rad = 6283 rad/s; St = 0.2 with Ω = f, and 1.2566 (= 2π × 0.2) with the angular frequency — the trap.
7. `nb.pointer` — **N79 [C]** "Vortex-induced vibration (spiral strakes break the spanwise coherence) and atmospheric Kármán streets behind mountains — stratification makes the flow two-dimensional; Ch. 13."
8. `nb.primer("a row of vortices: the cotangent sum", "A point vortex of strength Γ at z₀ induces the conjugate velocity w = u − iv = (Γ/2πi)/(z − z₀) (Ch. 5). An infinite row at z₀ + na (n = 0, ±1, ±2, …) adds up, pairing +n with −n so the sum converges, to
   (Γ/2ia) cot(π(z − z₀)/a). Similar lattice sums: Σ 1/(z − na)² = (π/a)²/sin²(πz/a) and Σ(−1)ⁿ/(z − na)² = (π/a)² cos(πz/a)/sin²(πz/a).", code="import numpy as np\nz = 0.3 + 0.2j                              # a point off the row (a = 1)\nn = np.arange(-2000, 2001)\nprint(np.sum(1/(z - n)), np.pi/np.tan(np.pi*z))   # both about 1.353-2.297j: the row sum equals pi cot(pi z)\nprint(np.sum(1/(z - n)**2), (np.pi/np.sin(np.pi*z))**2)   # squares: pi^2/sin^2")` (**P213**) — the builder prints the actual complex numbers; the two numbers of each print agree to about 1e-3 (truncation of the row at ±2000).
9. `nb.primer("linear stability of a steady configuration (perturb, linearise, eigenvalues)", "To test whether a steady arrangement is stable: displace it slightly, keep only terms linear in the displacement to get d(displacement)/dt = M × displacement, and find the eigenvalues λ of M
   (P80): a solution grows like e^{λt}, so Re λ > 0 for any eigenvalue means instability; if all Re λ ≤ 0 the configuration is stable (or neutral when Re λ = 0).", code="import numpy as np\nM = np.array([[0., 1.], [1., 0.]])       # d/dt (x, v) = (v, x): a saddle\nprint(np.linalg.eigvals(M))              # [ 1. -1.]: one growing mode (Re lambda = 1 > 0): unstable\nprint(np.linalg.eigvals(np.array([[0., 1.], [-1., 0.]])))   # +-i: neutral oscillation")` (**P214**)
10. `nb.derivation("D14", …)` — Part F D14 (15 steps), ref "" (the book states the ratio only); `check_src` (★★★) in Part F.
11. `nb.worked_example("the street speed and the wire", "At the stable spacing tanh(πb/a) = tanh(arccosh √2) = 1/√2 (since sech² = 1/2), so the street moves at U_s = (Γ/2a)/√2. For a = 5 cm and Γ = 0.5 m²/s that is 5 × 0.7071 = 3.54 m/s — slower than the 10 m/s wind that created it (the vortices are left behind). The row separation
    is b = 0.2805 a = 1.40 cm. For d = 2 mm in a 10 m/s wind, f = 1000 Hz.")`
12. `nb.code` — `print(BB.karman_street_ratio())`; `for r in (0.1, 0.2, 0.25, BB.karman_street_ratio(), 0.3, 0.5, 1.0): print(r, BB.karman_street_growth(r), BB.karman_street_growth_closed(r))`; `print(BB.karman_street_growth(0.3, offset=0.0), BB.karman_street_spectrum(BB.karman_street_ratio()))`;
    `print(BB.karman_street_velocity(0.05, 0.05*BB.karman_street_ratio(), 0.5))`. *expect:* 0.280550; growth (Γ = a = 1): 0.6400, 0.2982, 0.1099, 0 (< 1e-8), 0.0663, 0.5359, 0.7737 — identical to the closed form; non-staggered 0.7854; spectrum ±0.7854 i at the marginal spacing; street speed 3.536 m/s.
    *explain:* 1. `karman_street_growth` builds the linearised velocity of the periodic street (lattice sums) and takes the largest growth over the perturbation wavenumber; 2. the closed form (D14 step 13) is its value at k = π/a, the most dangerous mode; 3. a non-staggered row is unstable for every spacing.
13. `nb.check_agree` — **from scratch (curation §7):** numerical Jacobian of the induced velocity of a periodic cell (two vortices, images summed ±2000), central differences of step 1e-6, 4×4 real matrix at the alternating wavenumber (δ(−1)ⁿ), `np.linalg.eigvals`; `assert np.isclose(max(ev.real), BB.karman_street_growth_closed(b), rtol=1e-5)` for b/a = 0.2, 0.5; the sympy engine `ch09.karman_street_sympy()` gives ((μ − γ)² + σ²)((μ + γ)² + σ²). Also a one-vortex parity: `core.biot_savart.point_vortex_velocity` for a single vortex equals (Γ/2π r).
14. `nb.figure` — **the street and its stability** (9 × 3.4 in, three panels): (a) the staggered double row at b/a = 0.2805 (rows A blue Γ, B rose −Γ) with the street velocity arrow and streamlines in the street frame; (b) growth rate σ (Γ/a²) against b/a: numerical (dots), closed form (curve), zero at 0.2805 (green dot) and the
    non-staggered value π/4 (grey dashed); (c) eigenvalue spectrum in the complex plane for b/a = 0.15, 0.2805, 0.5 (three colours): quartets ±γ ± iσ on the imaginary axis at the marginal spacing. Title "Only one spacing does not grow". *see:* "a V-shaped growth curve touching zero"; *read:* "at 0.2805 all four eigenvalues sit on the imaginary axis (pure oscillation)"; *change:* "…offset 0 (non-staggered): σ = π/4 at every b/a".
15. `nb.animation` — **A5** (video, 60 frames FAST 30): the staggered street at b/a = 0.2805 advecting as an unchanged pattern, then at b/a = 0.15 with a small random kick that grows into a clump; the growth rate in the corner. Notes: *see* "an intact street versus a tearing one"; *read* "growth e-folds every 1/σ ≈ 2 time units (b/a = 0.15: σ = 0.48)"; *change* "…b/a = 0.5: growth again".
16. `nb.explainer("karman_street_stability", heading="Why is the vortex street's shape fixed at b/a ≈ 0.28?", why="Dragging b/a and the stagger changes the growth live, and a small kick visibly grows or not; the eigenvalue spectrum collapses onto the imaginary axis at exactly cosh(πb/a) = √2.", tries=["Drag b/a from 0.1 to 0.5: where does the growth touch zero?", "Set the stagger to 0: what happens?", "Press 'kick' at the stable spacing and at 0.15.", "Open Derivation D14 and watch the factorised polynomial."])`
17. `nb.md` — **What would change if…** "…viscosity and finite cores are added? The point-vortex street is only marginally stable at 0.2805; real streets are helped by viscous spreading. Kármán's result explains the shape, not the onset at Re ≈ 40 — that is Ch. 11."

18. `nb.core("C11", "Drag crisis: the layer turns turbulent, separation moves aft (82° → 125°) and $C_D$ falls — the base of sports-ball dynamics", question="Why does the drag fall when the flow gets faster?")`
19. `nb.md` — **The problem in plain words:** "Roughen a golf ball with dimples and it flies farther; a rough cylinder can have less drag than a smooth one; a cricket ball swings when its seam side is rough. In each case the same thing: transition to a turbulent layer happens at a lower speed on the rough side. Faster is not always more drag."
20. `nb.md` — **The idea** (ASCII): "laminar layer separates near 82° (from the front stagnation point) → wide wake, low pressure, C_D ~ 1; turbulent layer, fuller profile, holds on to ~125° → narrower wake, higher wake pressure, C_D ~ 0.3–0.6."
21. `nb.note` — **N81 [B]** "**Subcritical flow (Re below the crisis).** Laminar separation near 82°, a nearly flat wake pressure below the front value, and drag that is mostly form drag: C_D near unity and constant over decades of Re. The measured C_p distribution is drawn qualitatively by the separated model of C09 (Fig. 9.20 remade)."
22. `nb.code` — `for Re, rough in ((1e5, False), (3e5, False), (1e6, False), (1e5, True)): print(Re, rough, BB.drag_crisis_state(Re, rough=rough))`; `print(BB.cylinder_cd_schematic(np.logspace(-1, 7, 5)))`. *expect:* states labelled subcritical (φ_s ≈ 82°) / critical / supercritical (≈ 125°); the rough cylinder at 1e5 is already supercritical (roughness lowers Re_cr) — all marked qualitative; the numbers 82°, 125°, 3×10⁵ are the book's rounded experimental values. *explain:* the state table is a look-up of the book's rounded thresholds, not a calculation.
23. `nb.figure` — **Figs. 9.20 and 9.21 remade (qualitative)** (9 × 3.4 in, two panels): (a) C_p(φ): ideal (black dashed), subcritical (rose) and supercritical (teal) from `BB.separated_cp` with the illustrative C_b; (b) log–log C_D of a cylinder against Re (0.1–10⁷) as a **schematic labelled 'qualitative, no dataset'** with the four regimes shaded and the crisis dip; the caption says "the sphere curve (tested) is in §9.9". *see:* "two plateaus and a dip"; *read:* "before the crisis C_D ≈ 1, after ≈ 0.3–0.6 (book's rounded)"; *change:* "…rougher surface: the dip moves to smaller Re".
24. `nb.explainer("cylinder_drag_crisis", heading="Why does drag fall when the flow gets faster?", why="Re is a log-scale dial from 0.1 to 10⁷ that moves the schematic flow, the separation angle, the C_p distribution and C_D together, so the reader sees the causal chain, not four disconnected figures; the roughness toggle shows the crisis move.",
    tries=["Slide Re from 10 to 10⁷: name each regime as the status changes.", "Set 'rough' at Re = 10⁵ (before the crisis): what happened to C_D?", "Switch to the sphere mode: how does the dip compare?", "Change the wake pressure C_b: how sensitive is the drag?"])`
25. `nb.md` — **⚠️ Common confusion / the three counter-intuitive points (N82 [C]):** "(i) A tiny viscosity changes everything (ν → 0 is singular: d'Alembert's paradox is resolved by the boundary layer and separation, C01 and C09); (ii) a symmetric problem has an asymmetric solution (the street); (iii) roughness can *reduce* the drag of a blunt body. Ch. 11 and Ch. 14 return to them."

### A.9 §9.9 Flow Past a Sphere and the Dynamics of Sports Balls — C11 (continued: N83–N85 · N86) — no A item of its own
1. `nb.section("9.9", "Flow Past a Sphere and the Dynamics of Sports Balls", intro="**What is this section about?** The sphere shows the same drag crisis as the cylinder (at a higher Reynolds number, near 5 × 10⁵ for a smooth sphere) — and the crisis explains
   why a cricket ball swings, why a tennis ball can curve the 'wrong' way, and why a baseball knuckles. These are all consequences of C11, so the notes below continue that block; each has a figure and a number from our functions, none from the book.")`
2. `nb.md` — **C11 (continued): the sphere.** `nb.note` — **N83 [B]** "**Sphere regimes.** At low Re a doughnut-shaped attached eddy; above Re ≈ 130 the wake oscillates and sheds loops, not a regular street; transition at Re_cr ≈ 5 × 10⁵ with a sudden dip of C_D. > ⚠️ **Common confusion:** the book borrows the cylinder's 4 < Re < 40 for the sphere's steady
   eddy; sphere separation starts near Re ≈ 20." + `nb.code`: `Re = np.logspace(-1, 7, 200)`; `cd = SIM.sphere_drag_coefficient(Re, "morrison")`; `print(SIM.sphere_drag_coefficient(np.array([1e2, 1e5, 2e5, 4e5, 1e6])))`; `print(CRP.stokes_drag_coefficient(0.1), CRP.oseen_drag_coefficient(0.1))`.
   *expect:* C_D = 1.038 (Re = 100), 0.426 (1e5), 0.416 (2e5), 0.093 (4e5), 0.130 (1e6); Stokes 240.0 and Oseen 244.5 at Re = 0.1 (Morrison 240.2); Morrison's dip sits near 4 × 10⁵ (the book's rounded Re_cr ≈ 5 × 10⁵ — same crisis, different rounding).
   `nb.figure` — **Fig. 9.22 remade** (log–log C_D against Re, Morrison curve tested (black), Stokes 24/Re (blue dashed, valid < 1), Oseen 24/Re(1 + 3Re/16) (teal dashed), points A (before) and B (after) the crisis marked; badge "tested (V5)"). *see* "a −1 slope, a plateau near 0.4, a sudden dip";
   *read* "the dip at 3–5 × 10⁵ is the layer turning turbulent"; *change* "…a rough sphere: the dip appears earlier".
3. `nb.note` — **N84 [B]** "**Cricket swing.** A new ball with a smooth side and a seam-tripped side, at a speed just below the crisis: the smooth side stays laminar and separates early (≈ 85°), the seam side turns turbulent and separates late (≈ 120°): the pressure on the seam side is nearer the ideal minimum (−5/4 at the sphere's shoulder, (6.91) $C_{p,min}=-\frac54$) —
   a net side force F. A constant force gives a parabolic path $y=\frac12\frac Fm t^2=\frac12\frac FWg\Big(\frac dU\Big)^2$." equation $y=\tfrac12 a t^2$. `nb.code`: `print(BB.ball_swing_deflection(0.2, 18.0, 35.0))`; `for U in (25, 30, 35, 40): print(U, BB.ball_swing_deflection(0.2, 18.0, U))`. *expect:* 0.2595 m (t = 0.514 s); 0.5085 m at 25 m/s, 0.3532 at 30, 0.2595 at 35, 0.1986 at 40 — our F/W = 0.2, d = 18 m, not the book's.
   `nb.figure` — **swing schematic** (two panels): sketch of the ball with separation points 85° / 120° and C_p bars from the separated model (rose laminar side, teal seam side), and the parabolic trajectory y(x) for three speeds. *see/read/change* as usual ("…too slow: both sides laminar, no side force; too fast: both turbulent").
4. `nb.note` — **N85 [B]** "**Spin and the sign of Magnus.** Ch. 6's ideal Magnus lift L = ρUΓ (6.40) pushes a backspinning ball up. With a boundary layer the sign can flip: on a smooth ball at a speed below the crisis the side moving against the stream has the higher *relative* speed, its layer turns turbulent first, separates later and the force reverses (**negative Magnus effect**); a rough tennis ball has a lower Re_cr, both sides are past the crisis, and the faster side now separates *earlier*
   (**positive**). > ⚠️ **The book prints 'Re < Re_cr' twice; the second must be Re > Re_cr.**" + `nb.code`: `for lo, hi in ((0.8e5, 1.2e5), (2e5, 4e5), (4e5, 6e5)): print(lo, hi, BB.magnus_sign(lo, hi, 3e5))`. *expect:* '+' (both laminar), '−' (only the fast side past the crisis), '+' (both turbulent). `nb.figure` — **Fig. 9.25 remade**: a 3 × 2 truth table drawn as coloured tiles (sign, side Re vs Re_cr).
5. `nb.pointer` — **N86 [C]** "Baseball: a curveball has sidespin like the tennis ball; a knuckleball, with almost no spin, tumbles so the seam's position varies and the side force is irregular. See Ch. 6 (Magnus) and Ch. 14."
6. `nb.pointer` — the backup explainer `ball_swing_magnus` (B1) is built only if an explainer above fails review.
7. `nb.md` — **What would change if…** "…the surface has no wall at all? The same layer equations (9.18) describe a jet, where a wall's friction is replaced by entrainment. Next section."

---
### A.10 §9.10 Two-Dimensional Jets — C12 (+N87–N113 · D15 D16 D17 D18 · P215 P216 P217 · A6 · F6 · E7), R08 R09, C13 (+N114–N124 · D19 D20 D21 · P218 P219 · E8)
1. `nb.section("9.10", "Two-Dimensional Jets", intro="**What is this section about?** The layer equations (9.18) do not need a wall. A slot that blows fluid into still fluid makes a **free jet**: it is a boundary layer without a boundary, and (9.18) has an exact similarity solution, the
   sech² profile. Blown along a wall it makes a **wall jet**: a boundary layer under a free jet. Two conservation laws decide the shapes: in the free jet the momentum flux is constant and the mass flux grows (entrainment); in the wall jet the ordinary momentum flux dies at the wall, but a hidden 'flux of exterior momentum flux'
   survives. Jets and their entrainment recur in turbulent jets (Ch. 12) and plumes (Ch. 13).")`
2. `nb.core("C12", "The free two-dimensional laminar jet: $u=u_0\\,\\mathrm{sech}^2(\\eta/\\sqrt6)$ (9.71), $\\delta\\propto x^{2/3}$, entrainment $\\dot m=(36J\\rho^2\\nu x)^{1/3}$ (9.73)", question="A jet keeps its momentum but spreads and slows down — how do these fit, and where does the extra mass come from?")`
3. `nb.md` — **The problem in plain words:** "Air leaves a narrow slot at 30 m/s into a quiet room. A metre downstream the stream is wider, slower — and carries more air than left the slot. The jet drags in the room's air by viscous friction (entrainment). Momentum has nothing to push on: there is no wall and, with the pressure uniform, no net force — so the
   *momentum flux* stays the same all the way, while the *mass flux* grows. That single fact fixes how fast the jet slows and widens."
4. `nb.md` — **The idea** (ASCII):
   ```
   slot ══►   ──── x ────►      J = ρ∫u²dy  constant  (momentum flux, N/m)
   ══►     ╲   u₀(x) ∝ x^(−1/3) (slower)         ṁ = ρ∫u dy  grows ∝ x^(1/3)
   ══►  ───►  δ(x) ∝ x^(2/3)   (wider)           arrows in from both sides: entrainment
   ══►     ╱   profile: sech²(η/√6) at every x, η = y/δ(x)
   ```
5. `nb.note` — **N87 [B]** "**Set-up.** Slot in the x-direction into still fluid; boundary-layer approximation (∂/∂y ≫ ∂/∂x, v ≪ u) and no external pressure gradient (dp/dx = 0), so (9.18) $u\\frac{\\partial u}{\\partial x}+v\\frac{\\partial u}{\\partial y}=\\nu\\frac{\\partial^2u}{\\partial y^2}$ and (6.2) again. The jet is symmetric about y = 0.
   Wakes and shear layers obey the same equations." **N88–N90 [B]** "**Conditions:** $u=0$ for $y\\to\\pm\\infty$, $x>0$ (9.53); $v=0$ on the axis $y=0$ (9.54); an inlet profile $u=\\tilde u(y)$ on $x=x_0$ (9.55), forgotten far downstream." equations (9.53)–(9.55).
6. `nb.primer("momentum flux and mass flux through a cross-section", "Through a vertical line at x, per unit span, the *mass flux* is ṁ = ρ∫u dy [kg/(m s)] and the *momentum flux* is J = ρ∫u² dy [N/m] (the x-momentum carried across). A jet in still fluid at constant pressure has no force on it, so J cannot change with x (control volume, Ch. 4); ṁ can, because fluid enters through the sides.",
   code="import numpy as np\ny = np.linspace(-0.05, 0.05, 2001)                 # m, across a 10 cm wide stream\nu = 20*np.exp(-(y/0.01)**2)                        # a Gaussian jet, 20 m/s on the axis\nrho = 1.2\nprint(rho*np.trapezoid(u, y), rho*np.trapezoid(u**2, y))   # mdot = 0.4254 kg/(m s), J = 6.02 N/m")` (**P217**)
7. `nb.derivation("D15", …)` — Part F D15 (7 steps), ref "9.57". (Displays N91 (9.56) $\\int 2u\\frac{\\partial u}{\\partial x}dy+\\int\\big[u\\frac{\\partial v}{\\partial y}+v\\frac{\\partial u}{\\partial y}\\big]dy=\\int\\frac{\\partial\\tau}{\\partial y}dy$ (⚠️ printed τ lacks 1/ρ: kinematic stress), N92 (9.57), N93 (9.58) $\\int u^2dy=\\text{const}=J/\\rho$.)
8. `nb.md` — **Similarity from the invariant.** `nb.note` — **N94 [B], N95 [B], N96 [B]** "Look for $\\psi=u_0(x)\\delta(x)f(\\eta)$, $\\eta=y/\\delta(x)$, $\\delta=[\\nu x/u_0(x)]^{1/2}$ (9.59) — the ansatz of C03 with the outer speed replaced by the unknown centre-line speed u₀; then $u=\\partial\\psi/\\partial y=u_0f'(\\eta)$ (9.60) and
   $\\frac J\\rho=u_0^2\\delta\\int f'^2d\\eta$ (9.61). The number $C\\equiv\\int f'^2d\\eta$ is a constant: the x-dependence of J/ρ must cancel, which fixes u₀ and δ (next)." equations (9.59)–(9.61).
9. `nb.note` — **N97 [B], N98 [B], N99 [B]** "**Exponents:** $u_0(x)=\\big[J^2/(C^2\\rho^2\\nu x)\\big]^{1/3}\\propto x^{-1/3}$ (9.62), $\\delta(x)=\\big[C\\rho\\nu^2x^2/J\\big]^{1/3}\\propto x^{2/3}$ (9.63), $\\psi=[J\\nu x/C\\rho]^{1/3}f(\\eta)$, $\\eta=y\\big/[C\\rho\\nu^2x^2/J]^{1/3}$ (9.64). With J = 1 N/m in air at x = 0.1 m: u₀ = 35.1 m/s, δ = 0.207 mm; at x = 0.4 m: 22.1 m/s and 0.521 mm (× 4^{−1/3} = 0.63 and × 4^{2/3} = 2.52)." equations (9.62)–(9.64).
10. `nb.derivation("D16", …)` — Part F D16 (14 steps), ref "9.65" — the reduction of (9.65) $\\frac{\\partial\\psi}{\\partial y}\\frac{\\partial}{\\partial x}\\Big(\\frac{\\partial\\psi}{\\partial y}\\Big)-\\frac{\\partial\\psi}{\\partial x}\\frac{\\partial}{\\partial y}\\Big(\\frac{\\partial\\psi}{\\partial y}\\Big)=\\nu\\frac{\\partial^2}{\\partial y^2}\\Big(\\frac{\\partial\\psi}{\\partial y}\\Big)$ (N100) to $3f'''+ff''+f'^2=0$ (N101, the book skips the
    differentiation; sympy-checked) with $f'\\to0$ as η → ±∞ (9.66), $f'=1$ on η = 0 (9.67), $f=0$ on η = 0 (9.68) (N102–N104).
11. `nb.primer("sech, arccosh and (tanh)′ = sech²", "sech x = 1/cosh x, a bell-shaped curve equal to 1 at 0 that dies like 2e^{−|x|}; d(tanh x)/dx = sech² x = 1 − tanh² x. Its inverse: arccosh y = ln(y + √(y² − 1)) for y ≥ 1 answers 'where does sech reach 0.1?'. (Extends the tanh primer of Ch. 7.)",
    code="import numpy as np\nprint(1/np.cosh(0.0), 1/np.cosh(2.0)**2)         # 1.0 and 0.0707: sech^2 at x = 0 and 2\nprint(np.arccosh(10.0), np.log(10 + np.sqrt(99)))   # 2.9932 twice: the point where sech = 0.1")` (**P215**)
12. `nb.primer("the total-derivative move (look for the derivative of a product)", "Before integrating an ODE, ask whether the left side is already the derivative of something: (f f′)′ = f′² + f f″ and (3f′ + f²/2)′ = 3f″ + f f′ (product rule read backwards, Ch. 3 P38). Then integrating is free: the ODE becomes 'that something = constant'.",
    code="import sympy as sp\nx = sp.symbols('x')\nf = sp.Function('f')\nprint(sp.expand(sp.diff(f(x)*f(x).diff(x), x)))           # f'^2 + f f''\nprint(sp.expand(sp.diff(3*f(x).diff(x) + f(x)**2/2, x)))  # 3 f'' + f f'")` (**P216**)
13. `nb.derivation("D17", …)` — Part F D17 (11 steps), ref "9.71". (Displays N105 (9.69) $3f''+ff'=0$, then $3f'+\\frac{f^2}2=3$; N106 (9.70) $\\tanh^{-1}\\frac f{\\sqrt6}=\\frac\\eta{\\sqrt6}$; N107 (9.72) $C=\\frac{4\\sqrt6}3$; N108, N109 (9.73); N110 (9.74); N111 (9.75).)
14. `nb.note` — **N108–N111 [B] summary** "$u(x,y)=u_0(x)\\,\\mathrm{sech}^2\\Big(\\frac y{\\sqrt6}\\Big[\\frac J{C\\rho\\nu^2x^2}\\Big]^{1/3}\\Big)$ (9.71) (Bickley's jet, 1937: u₀ = 0.4543(J²/(ρ²νx))^{1/3}); mass flux $\\dot m=\\rho u_0\\delta\\,2\\sqrt6=(36J\\rho^2\\nu x)^{1/3}$ (9.73), growing like x^{1/3}: entrainment; cross-flow
    $v=-\\frac13\\big(\\frac{J\\nu}{C\\rho x^2}\\big)^{1/3}[f-2\\eta f']$ (9.74) and $\\frac v{u_0}\\to\\mp\\frac{\\sqrt6}{3\\sqrt{Re_x}}$ as η → ±∞ (9.75): ambient fluid flows toward the jet from both sides, with $Re_x=xu_0/\\nu$." equations (9.71), (9.73)–(9.75).
15. `nb.derivation("D18", …)` — Part F D18 (5 steps), ref "9.76". `nb.note` — **N112 [B]** "> ⚠️ **The book prints h₉₉ = 5.6152[Cρν²x²/J]^{1/3} (9.76) from sech² = 0.01 → 2.2924. The 1 % point is arccosh 10 = 2.9932 (h₉₉ = 7.3319[…]^{1/3}); 2.2924 = arccosh 5 is the 4 % point.**" and **N113 [B]** "$Re_x=\\big(\\frac{3Jx}{4\\sqrt6\\rho\\nu^2}\\big)^{2/3}$ and $Re_{h_{99}}=h_{99}u_0/\\nu$ (with the corrected coefficient 7.3319);
    a laminar jet at Re ≫ 1 is unstable (the sech² profile has inflection points: Ch. 11); the turbulent jet spreads like x, faster than the laminar x^{2/3} (Ch. 12 — the book's cross-reference to 'Chapter 13' is a slip)." equation (9.76).
16. `nb.worked_example("a slot jet in air and in water", "J = 1 N/m. Air (ρ = 1.2, ν = 1.5 × 10⁻⁵), x = 0.1 m: u₀ = 35.14 m/s, δ = 0.2066 mm, ṁ = 0.04268 kg/(m s), h₉₉ = 7.3319 δ = 1.515 mm, Re_x = 2.34 × 10⁵, edge cross-flow v = −0.0593 m/s (toward the jet). Water (ρ = 1000, ν = 10⁻⁶): u₀ = 0.979 m/s, δ = 0.320 mm,
    ṁ = 1.533 kg/(m s), h₉₉ = 2.344 mm. Four times further downstream in air: u₀ down to 63 %, width up 2.52 ×, ṁ up 1.59 ×, but J still 1 N/m.")`
17. `nb.code` — `k = JET.free_jet_constants()`; print. `x = np.array([0.05, 0.1, 0.2, 0.4])`; `for xi in x: print(xi, JET.free_jet_centreline(xi, 1.0, RHO_AIR, NU_AIR), JET.free_jet_thickness(xi, 1.0, RHO_AIR, NU_AIR), JET.free_jet_mass_flux(xi, 1.0, RHO_AIR, NU_AIR), JET.free_jet_halfwidth(xi, 1.0, RHO_AIR, NU_AIR))`;
    `print(JET.free_jet_halfwidth(0.1, 1.0, RHO_AIR, NU_AIR, printed=True), JET.free_jet_reynolds(0.1, 1.0, RHO_AIR, NU_AIR))`. *expect:* C = 3.26599, mdot_coeff 3.30193, arccosh10 = 2.99322, h99_coeff = 7.33187, printed 5.61529, u0_coeff = 0.45428, f_inf = 2.44949; the table of the worked example; printed h₉₉ = 1.160 mm (23 % too small) vs 1.515 mm.
    *explain:* 1. `free_jet_constants` computes C by `quad` and confirms 4√6/3; 2. the four functions are (9.62), (9.63), (9.73), (9.76) corrected; 3. `printed=True` reproduces the slip.
18. `nb.code` — momentum conservation and Bickley: `y = np.linspace(-40, 40, 8001)*JET.free_jet_thickness(x0, ...)`; for x in (0.05, 0.1, 0.2, 0.4): `u = JET.free_jet(x, y, 1.0, RHO_AIR, NU_AIR)["u"]`; `print(JET.jet_momentum_flux(y, u, RHO_AIR))`; `print(k["u0_coeff"], 1/(np.sqrt(6)*k["C"]**(1/3)), k["mdot_coeff"])`.
    *expect:* J = 1.0000 N/m at every x (agreement 1e-8); 0.45428, 0.27520, 3.30193 (Bickley's 0.4543, 0.2752, 3.3019). *explain:* the integral of u² is the same at four stations while u₀ and δ move: the invariant of D15.
19. `nb.check_agree` — **from scratch (curation §7):** `solve_bvp` on 3f‴ + ff″ + f′² = 0 with f(0) = 0, f′(0) = 1, f′(η_max) = 0 (η_max = 12), from a guess tanh; `JET.free_jet_ode_solve()`; `assert np.allclose(f_bvp, np.sqrt(6)*np.tanh(eta/np.sqrt(6)), atol=1e-7)`; `assert np.isclose(np.trapezoid(fp**2, eta)*2, 4*np.sqrt(6)/3, rtol=1e-6)` (∫f′² over ±∞ is twice the half-line); the sympy engine `ch09.similarity_reduce_sympy("free_jet")["residual"] == 0`.
20. `nb.figure` — **the free jet** (9 × 3.4 in, three panels): (a) the jet in the x–y plane (air, J = 1 N/m, x ∈ [0, 0.5] m, y in mm) coloured by u/u₀ with streamlines bending inward (entrainment) and the h₉₉ envelope; (b) u(y) at x = 0.05, 0.1, 0.2, 0.4 m (raw, blue) and the same divided by u₀ against η (purple dashed, sech²(η/√6) exact); (c) log–log u₀, h₉₉ and ṁ against x with fitted slopes −1/3, 2/3, 1/3 printed.
    Title "A jet keeps J, loses speed, gains mass". *see:* "stretching sech² curves that collapse; arrows into the jet"; *read:* "the area under u² is the same at every x; the area under u grows"; *change:* "…doubling J: u₀ ×1.59, δ ×0.79, ṁ ×1.26 (2^{2/3}, 2^{−1/3}, 2^{1/3})".
21. `nb.plotly` — **F6** `slider_figure` over x (25 values 0.03–0.6 m): free-jet profile with the wall-jet profile (C13) for comparison; the exponents in the title.
22. `nb.animation` — **A6** (video, 60 frames FAST 30): the free jet and the wall jet spreading side by side over x with their profiles at successive stations, δ ∝ x^{2/3} versus x^{3/4} and ṁ ∝ x^{1/3} versus x^{1/4} counters (shared with C13).
23. `nb.explainer("free_jet_similarity", heading="A jet keeps its momentum, spreads and slows — where does the extra mass come from?", why="Moving x and J changes u₀, δ and ṁ together while the momentum-flux integral stays fixed; entrainment arrows lengthen toward the axis; static profiles hide the invariant.",
    tries=["Drag x from 0.05 to 0.6 m: read the slopes −1/3, 2/3, 1/3.", "Change the half-width level from 1 % to 4 %: which coefficient does the book's 5.6152 belong to?", "Click a point in the jet: v and u with the arithmetic.", "Open Derivation D16 and watch the powers of x cancel."])`
24. `nb.md` — **What would change if…** "…a wall is added under the jet? Friction removes momentum flux, so J is no longer conserved; a subtler invariant appears (C13)."

25. `nb.recap("R08", "Wall conditions of the wall jet", "Eq. (9.77) $u=v=0$ on $y=0$, $x>0$: no slip (9.12) and no through-flow (9.13) again.", where="§9.1, R04 R05")` · `nb.recap("R09", "Far-field condition", "Eq. (9.78) $u(x,y)\\to0$ as $y\\to\\infty$: the still ambient fluid, the one-sided version of (9.53) and of (9.14) with U_e = 0.", where="§9.1, §9.10")`
26. `nb.core("C13", "The wall jet: the invariant $\\frac d{dx}\\int_0^\\infty u\\Big(\\int_y^\\infty u^2dy'\\Big)dy=0$ (9.80), the ODE $4f'''+ff''+2f'^2=0$ and the implicit solution (9.83)", question="In a wall jet the wall eats momentum — so what is conserved, and why does the jet spread as x^{3/4} instead of x^{2/3}?")`
27. `nb.md` — **The problem in plain words:** "A jet of air runs along a ceiling; a pressure-washer's spray slides along a wall; a cooling jet hugs a turbine blade. Near the wall the fluid sticks (a boundary layer); far from it the fluid moves like a free jet. Because the wall pulls back, the jet's ordinary momentum flux ρ∫u²dy *decreases* downstream. Yet there is still a conserved
    quantity, and it fixes how the jet spreads."
28. `nb.md` — **The idea** (ASCII):
    ```
      slot ═►  free-jet-like outer part ─────────►       ordinary momentum flux ∫u²dy  ↓ (wall shear)
      ═══════ wall ═══ inner boundary layer, u = 0 at y = 0    but  ∫₀^∞ u ( ∫_y^∞ u² dy′ ) dy  = constant
    ```
29. `nb.note` — **N114 [B]** "**Laminar wall jet (Glauert 1956).** Slot along a plane wall; near the wall a boundary layer, far away a free jet; p ≈ const; (6.2) and (9.18) again, with the conditions u = v = 0 on y = 0 (9.77) and u → 0 as y → ∞ (9.78) (recaps above)." equations (9.77), (9.78).
30. `nb.primer("integration by parts with a variable lower limit", "∫₀^∞ u(y) G(y) dy with G(y) = ∫_y^∞ g dy′ can be integrated by parts: G′ = −g, so ∫₀^∞ u G dy = [W G]₀^∞ + ∫₀^∞ W g dy with W = ∫₀^y u. Differentiating a double integral with respect to x adds only the x-derivative under the integral signs, since the limits 0, ∞, y do not depend on x.",
    code="import numpy as np\nfrom scipy.integrate import quad\nG = lambda y: np.exp(-2*y)/2                       # G(y) = int_y^inf u^2 dy for u = exp(-y)\nprint(quad(lambda y: np.exp(-y)*G(y), 0, np.inf)[0])                 # 0.16667: int u G dy\nprint(quad(lambda y: (1 - np.exp(-y))*np.exp(-2*y), 0, np.inf)[0])   # 0.16667 again: by parts, W = int_0^y u = 1 - exp(-y) times u^2")` (**P218**) — both prints give 1/6 = 0.16667: integrating by parts trades the double integral for a single one.
31. `nb.derivation("D19", …)` — Part F D19 (13 steps), ref "9.80". (Displays N115 the chain, N116 (9.79), N117 (9.80), N118 (9.81); ends with $xu_0^2=C^2$, $u_0=Cx^{-1/2}$.) `nb.note` — **N119 [B]** "**Exponents:** $u_0=Cx^{-1/2}$, $\\psi(x,y)=[\\nu Cx^{1/2}]^{1/2}f(\\eta)$, $\\eta=y/\\delta(x)$, $\\delta(x)=[\\nu x^{3/2}/C]^{1/2}\\propto x^{3/4}$ (9.82) — thicker than the free jet's x^{2/3}: the wall slows the fluid so the layer must spread more." equation (9.82).
32. `nb.derivation("D20", …)` — Part F D20 (12 steps), ref "9.82". `nb.note` — **N120 [B]** "> ⚠️ **The book prints the reduced equation as f‴ + ff″ + 2f′² = 0. The correct one is $4f'''+ff''+2f'^2=0$** (sympy: residual −(4f‴ + ff″ + 2f′²)C²/(4x²)); the book's own next line 4ff″ − 2f′² + f²f′ = 0 needs the 4.
    Wall conditions $f(0)=0,\\ f'(0)=0$ and $f'(\\infty)=0$." + `nb.code`: `print(ch09.similarity_reduce_sympy("wall_jet")["ode"], ch09.similarity_reduce_sympy("wall_jet", printed=True)["residual"])`. *expect:* `4*f''' + f*f'' + 2*f'^2 = 0`; the printed variant leaves the non-zero residual −(3/4)C²f‴/x² (coefficient difference 3).
33. `nb.primer("partial fractions and sympy apart; inverting an implicit solution", "A rational function such as 1/(1 − g³) is split into simple pieces that can be integrated (partial fractions); `sympy.apart` does it. When an integral gives η as an explicit function of g (η(g)) but we want g(η), invert numerically: for each η, `brentq` finds the g with η(g) − η = 0.",
    code="import sympy as sp\nfrom scipy.optimize import brentq\ng = sp.symbols('g')\nprint(sp.apart(1/(1 - g**3), g))                        # 1/(3*(1 - g)) + (g + 2)/(3*(g**2 + g + 1))\nprint(brentq(lambda t: t**3 + t - 2, 0, 2))             # 1.0: invert t^3 + t = 2 numerically")` (**P219**)
34. `nb.derivation("D21", …)` — Part F D21 (15 steps), ref "9.83". `nb.note` — **N121 [B]** "> ⚠️ **The book's separation of variables reads ∫df/(f_∞^{3/2}f − f²) = ⅙∫dη; the correct integrand is 1/(f_∞^{3/2}f^{1/2} − f²)** (from f′ = f^{1/2}(f_∞^{3/2} − f^{3/2})/6); the far-field statement 1 − g ≈ e^{−f_∞η/4} also omits the factor 4.29 = √3 e^{√3π/6}."
35. `nb.note` — **N122 [B]** "**Entrainment:** $\\dot m=\\int_0^\\infty\\rho u\\,dy=\\rho u_0\\delta\\int f'd\\eta=\\rho\\sqrt{\\nu C}\\,f_\\infty x^{1/4}$ (9.84): ∝ x^{1/4}, less than the free jet's x^{1/3} — the wall entrains from one side only." equation (9.84).
    **N123 [B]** "**The constants (9.85).** $u_0^2(x)\\,\\nu x\\int_0^\\infty\\Big(f'\\int_\\eta^\\infty f'^2d\\eta'\\Big)d\\eta=C^2\\nu\\int_0^\\infty\\Big(f'\\int_\\eta^\\infty f'^2d\\eta'\\Big)d\\eta=\\Psi$: the invariant fixes C; ṁ at one x fixes f_∞. **Our reading:** because f′(0) = 0 (unlike the free jet's f′(0) = 1) the ODE has a free scale, f → λf(λη), which sends
    C → C/λ² and f_∞ → λf_∞; the only physical constant is the combination $C f_\\infty^2$. With the gauge f_∞ = 1 the invariant integral is 1/40 (in general f_∞⁴/40) and Ψ = ṁ⁴/(40ρ⁴νx) — so Ψ and ṁ are the same datum, not two. > ⚠️ Ψ has units m⁴/s³ (u³L²), not a force per length." equation (9.85).
36. `nb.worked_example("a wall jet in air, gauge f_∞ = 1", "For f_∞ = 1 the IVP from f″(0) = f_∞³/72 = 1/72 = 0.013889 gives ∫f′dη = 1, ∫f′²dη = 1/18 = 0.0556 and ∫f′(∫_η^∞f′²) = 1/40 = 0.025; the profile peaks at f′ = 0.0787 at η = 8.11. If ṁ(1 m) = 0.05 kg/(m s) in air: √(νC) = ṁ/(ρ f_∞ x^{1/4}) = 0.04167, C = 115.7 m^{3/2}/s,
    so u₀(1 m) = 115.7 m/s, the peak speed u₀ f′_max = 9.11 m/s, δ = (νx^{3/2}/C)^{1/2} = 0.360 mm and the peak sits at y = 8.11 δ = 2.9 mm. At 4 m: ṁ ×√2 = 0.0707, the peak speed ×½ = 4.56 m/s and δ × 4^{3/4} = 2.83.")`
37. `nb.code` — `sol = JET.wall_jet_ode_solve(fpp0=1/72)`; print `sol["f_inf"], sol["err_vs_9_83"], sol["fpp0_over_finf_cubed"], sol["fp"].max()`; `print(JET.wall_jet_integrals(1.0))`; `bad = JET.wall_jet_ode_solve(fpp0=1/72, coeff=1.0)`; print `bad["err_vs_9_83"]`; `pr = JET.wall_jet_profile(np.array([2., 8., 20.]), f_inf=1.0)`; print `pr["fp"]`;
    `c = JET.wall_jet_constants(RHO_AIR, NU_AIR, mdot=0.05, x=1.0, f_inf=1.0)`; print c. *expect:* f_inf = 1.000000, error vs (9.83) < 1e-9, ratio 1/72 = 0.0138889, f′ max = 0.078745 (at η = 8.11); integrals (1, 0.055556, 0.025, 0.013889); the printed ODE (same f″(0)) ends at f_∞ = 0.397 instead of 1 and misses the (9.83) relation by 0.6–3; f′(2, 8, 20) = 0.0276, 0.0787, 0.0133 (checked against the IVP to 1e-9);
    C = 115.7, Psi = 0.05⁴/(40 × 1.2⁴ × 1.5 × 10⁻⁵ × 1) = 5.0 × 10⁻³ … (units m⁴/s³). *explain:* 1. the IVP from f″(0) = f_∞³/72 lands on f_∞ = 1 without being told; 2. the closed forms of the three integrals follow from the scaling; 3. the wrong ODE fails the implicit relation; 4. `wall_jet_constants` shows Ψ and ṁ are one datum.
38. `nb.check_agree` — **from scratch (curation §7):** `solve_ivp` of the ODE from f″(0) = f_∞³/72 (rtol 1e-12) for f_∞ = 1 and 2 (f″(0) = 1/72, 8/72); check the implicit relation (9.83) at η = 2, 5, 10 with a hand-written `brentq` inverse; `assert np.isclose(sol["y"][0, -1], f_inf, rtol=1e-8)`; scaling test: f₂(η) = 2 f₁(2η) reproduces the f_∞ = 2 solution (`assert np.allclose(...)`, rtol 1e-7).
39. `nb.figure` — **Fig. 9.29 remade and the two invariants** (9 × 3.4 in, three panels): (a) f(η) (blue) and f′(η) (rose) for f_∞ = 1, η ∈ [0, 30], with the exponential tail 4.29-corrected (purple dashed) and the peak marked; (b) log–log growth of δ (free jet slope 2/3, wall jet 3/4) and ṁ (1/3 vs 1/4) with J and Ψ fixed to the same reference at x = 0.1 m; (c) bars at
    x = 0.1, 0.2, 0.4, 0.8 m for the wall jet: ∫u²dy (rose, falling like x^{−1/4}) and the double integral (teal, constant), with the free jet's flat ∫u²dy as a grey ghost. Title "When the wall eats the obvious invariant, the next one survives". *see:* "a falling rose bar and a flat teal one"; *read:* "wall jet: ∫u²dy decays like x^{−1/4} while the double integral is constant"; *change:* "…f_∞ doubled: C quarters, but every dimensional profile stays put (gauge freedom)".
40. `nb.explainer("wall_jet_invariant", heading="In a wall jet the wall eats momentum — what is conserved?", why="In free-jet mode ∫u²dy is constant with x; in wall-jet mode it decreases but the double integral remains constant — both bars at the same instant; the f_∞ slider shows the scaling symmetry that makes f_∞ free.", tries=["Switch between 'free' and 'wall': which bar stays flat?", "Move the f_∞ slider: the dimensional profile does not change (C compensates).", "Turn on 'ODE as printed' and read the residual.", "Compare δ ∝ x^{3/4} with x^{2/3}."])`
41. `nb.md` — **What would change if…** "…the jet were round, or turbulent? The same two lessons (constant momentum flux, growing mass flux) hold; only the exponents change (Ch. 12)."

### A.11 §9.11 Secondary Flows — C14 (+D22 · P220 · E9), S01
1. `nb.section("9.11", "Secondary Flows", intro="**What is this section about?** So far a boundary layer only slows the flow. When the streamlines outside the layer are curved, the layer does more: it creates a flow *across* the main flow. The teacup is the everyday example — stir it and the tea leaves go to the middle, not to the rim. The same mechanism (a friction layer in a rotating fluid) makes the Ekman layers of the ocean and atmosphere in Ch. 13.")`
2. `nb.core("C14", "Secondary flow: the teacup — a slowed layer under an unchanged pressure gradient", question="Why do tea leaves collect at the centre when the water is spinning outward?")`
3. `nb.md` — **The problem in plain words:** "Stir tea and let it spin. The leaves pile up in the middle of the cup's floor. Spinning water is flung outward, so why do heavy leaves move inward? Because the water next to the floor is slowed by friction, but the pressure difference that holds the fast water on its circle is not reduced. In the slow layer the pressure gradient wins: it pushes the water — and the leaves — inward."
4. `nb.md` — **The idea** (ASCII cup section):
   ```
     axis                         free surface (slightly dipped: the pressure gradient ∂p/∂R = ρu²/R)
      │ ▲ up      ◄── out ──                 fast core: circular streamlines, ∂p/∂R = ρ u_e²/R
      │ │        ╱                            slowed floor layer: u < u_e, centrifugal ρu²/R too small
      │ │  in ───►  ●leaves  ─ floor layer ─   net force ρ(u_e² − u²)/R pointing inward
   ```
5. `nb.primer("radial force balance in a swirl (and a thin layer)", "A fluid parcel going round a circle of radius R at speed u needs an inward force per volume ρu²/R (the centripetal acceleration, Ch. 4); in the fast core that is supplied by the pressure gradient, ∂p/∂R = ρu_e²/R. In a thin layer the pressure hardly changes across the layer (∂p/∂z ≈ 0, the boundary-layer argument (9.10) with z for y),
   so the same ∂p/∂R acts on the slower fluid, which needs less: the surplus pushes it inward.", code="rho, ue, R = 1000.0, 0.2, 0.04          # water, core swirl 0.2 m/s, radius 4 cm\ndpdR = rho*ue**2/R                       # 1000 N/m^3: the radial pressure gradient set by the core\nfor u in (0.2, 0.1, 0.0):                # core, half speed, on the floor\n    print(u, dpdR - rho*u**2/R)          # net inward force per volume: 0, 750, 1000 N/m^3")` (**P220**)
6. `nb.derivation("D22", …)` — Part F D22 (7 steps), ref "".
7. `nb.worked_example("water in a 4 cm-radius cup", "ρ = 1000 kg/m³, core swirl u_e = 0.2 m/s, R = 4 cm. The core needs ρu_e²/R = 1000 N/m³. On the floor u → 0, so the net inward force is 1000 N/m³ = 10 % of the weight density ρg (9810 N/m³): small but unopposed by any centrifugal force; at half speed 750, at 0.75 u_e 437.5, at u_e zero.")`
8. `nb.code` — `z = np.linspace(0, 5e-3, 101)`; `u = ch09.secondary_flow_layer_profile(z, delta=1.5e-3, u_e=0.2, shape="exponential")` (illustrative: u = u_e(1 − e^{−z/δ})); `F = ch09.secondary_flow_radial_force(0.2, u, 0.04, 1000.0)`; print `F[0], F[50], F[-1]` and the layer depth where F falls to 10 % of its floor value; `print(ch09.secondary_flow_radial_force(0.2, 0.0, 0.04, 1000.0))`.
   *expect:* F(0) = 1000, F(z = 2.5 mm) = 1000(1 − (1 − e^{−1.667})²) = 342, F(5 mm) = 1000(1 − (1 − e^{−3.33})²) = 70 N/m³; the height where F falls to 10 % of the floor value is z = 2.97 δ = 4.46 mm (exponential profile, illustrative); floor value 1000.0. *explain:* 1. the profile is *illustrative*, not a solution of anything (docstring says so); 2. the force is ρ(u_e² − u²)/R; 3. it is largest on the floor and vanishes in the core.
9. `nb.check_agree` — **from scratch (curation §7):** a two-line balance `dpdR = rho*ue**2/R; F_mine = dpdR - rho*u**2/R` on the sample profile; `assert np.allclose(F_mine, F)`; sign test `assert (F >= 0).all()`.
10. `nb.figure` — **the teacup** (9 × 3.4 in, three panels): (a) a cup in cross-section (`cup_section`) with the swirl arrows, the meridional loop (in along the floor, up the axis, out along the top, down the side wall) and leaves drifting inward (schematic, labelled); (b) profiles versus height: u(z)/u_e (blue), pressure-gradient force ρu_e²/R (orange, uniform in z), centrifugal ρu²/R (rose) with the difference shaded (the net inward force);
    (c) F(z) alone. Title "A slowed layer under a fixed pressure gradient is pushed inward". *see:* "the shaded wedge between orange and rose vanishing at the top of the layer"; *read:* "the net force is 1000 N/m³ on the floor, 70 N/m³ at 5 mm"; *change:* "…a thicker layer (viscous oil): the wedge is taller and the inflow reaches higher; a spinning-up cup (u_e falling): the wedge shrinks".
11. `nb.md` — **The Ekman hook (preview, not in this chapter's book text):** "The inflow speed follows from balancing this force with viscosity: the layer thickness is δ ~ √(ν/Ω) with Ω = u_e/R = 5 s⁻¹, i.e. 0.45 mm for water, and the water in the cup spins down (or up) in a time ~ H/√(νΩ) = 0.05/√(10⁻⁶ × 5) = 22 s for a 5 cm depth — far faster than the pure-diffusion time H²/ν = 2500 s. This is the mechanism of spin-up
    and of the Ekman layers of the ocean and atmosphere (Ch. 13); a river bend (Exercise 9.28) has the same balance across the channel."
12. `nb.explainer("teacup_secondary_flow", heading="Why do tea leaves collect at the centre?", why="Dragging the height z shows the imbalance ρ(u_e² − u²)/R growing toward the floor; changing the swirl profile shape or the layer thickness shows how the inflow depends on friction; a static picture cannot show the balance profile by profile.", tries=["Drag the height slider from the core to the floor.", "Change the profile shape and the layer thickness.", "Switch to 'river bend' mode: the same balance across a channel.", "Watch the force bars add to the net inward force."])`
13. `nb.pointer` — **S01 [SKIP]** "The 28 exercises (plate in two orientations, RK solution of the Blasius equation, pipe entry length, Thwaites for a horn, a wavy wall, perpetual separation, a round jet, a river bend, …), the Literature Cited and the Supplemental Reading are not reproduced here; the exercises' numerical answers are computable from the functions above and are kept out of this repository. Chapter 10 solves the boundary-layer equations numerically."
14. `nb.summary(clicked=[…14…], feeds_forward=[…], left_out=[…])`. **clicked (one sentence per CORE):** C01 "Where advection and viscosity balance the layer is δ ~ L Re^{−1/2} thick, and the equations (9.9)–(9.10) with the pressure from (9.11) are all that is left." · C02 "δ₉₉, δ* and θ are three integrals of one profile; H = δ*/θ tells how full it is." · C03 "A plate has no length scale, so ψ = Uδ(x)f(η) turns the PDE into f‴ + ½ff″ = 0." · C04 "One shot with f″(0) = 1 and a rescale gives 0.3321; every Blasius number is an integral of that curve." · C05 "A power-law outer flow keeps the similarity form; n is the dial from stagnation flow through Blasius to separation at n = −0.0904." · C06 "Integrating the momentum equation across the layer gives (9.43), exact for any layer, at the price of a closure." · C07 "If shear and shape depend on λ alone, (9.43) integrates to θ² ∝ ∫U_e⁵ and separation is λ crossing −0.09 (or −0.0681)." · C08 "At the wall μ u_yy = dp/dx: an adverse gradient bends the profile, creates an inflection and drives τ₀ to zero." · C09 "Separation replaces the pressure recovery by a flat wake pressure: C_D,p = sin φ_s (1 − 4/3 sin²φ_s − C_b), so later separation means less drag." · C10 "The wake of a cylinder sheds a street of alternating vortices, and only b/a = 0.2805 (cosh πb/a = √2) does not grow." · C11 "A turbulent layer separates later (82° → 125°) and the drag falls; spin, seam and roughness decide which side of a ball is past the crisis." · C12 "A free jet keeps J = ρ∫u²dy, so u₀ ∝ x^{−1/3}, δ ∝ x^{2/3}, sech² profile, and ṁ ∝ x^{1/3} grows by entrainment." · C13 "A wall jet loses ordinary momentum flux to the wall, but ∫u(∫_y^∞u²) is invariant: u₀ ∝ x^{−1/2}, δ ∝ x^{3/4}, and f_∞³ = 72 f″(0)." · C14 "A friction layer under an unchanged pressure gradient is pushed inward: the teacup, and the door to Ekman layers.".
    **feeds_forward:** Ch. 10 (marching solver and Blasius/Falkner–Skan as test problems; `BL` as the reference), Ch. 11 (inflection points of (9.52), jet and wake profiles sech² as unstable shear flows, Blasius profile for Orr–Sommerfeld, Tollmien–Schlichting), Ch. 12 (θ, H, C_f laws for turbulent layers, turbulent jets spreading ∝ x, the plate drag curve), Ch. 13 (Ekman layer from C14, stratified Kármán streets, spin-up), Ch. 14 (separation and stall, θ for airfoils, Thwaites on an airfoil, Blasius force theorem again), Ch. 15 (compressible boundary layers). **left_out:** numerical solution of the layer equations (Ch. 10), transition mechanism (Ch. 11), turbulent layers and jets (Ch. 12), Ekman layers (Ch. 13), Table 9.1 of Thwaites (replaced by our exact-FS closure), all exercises.

### A.12 Placement check (every curation id has exactly one home)
- CORE C01–C14: C01 §9.1 · C02 §9.2 · C03, C04 §9.3 · C05 §9.4 · C06 §9.5 · C07 §9.6 · C08, C09 §9.7 · C10, C11 §9.8 (C11 continued in §9.9) · C12, C13 §9.10 · C14 §9.11 — **14 blocks, each with code and at least one visual.**
- RECAP: R01, R02, R03 before C01 · R04, R05 after C01 · R06, R07 before C03 · R08, R09 before C13 (9 recaps).
- NOTE ids: N01–N17 in C01 (N15, N17 pointers) · N18–N21 in C02 · N22–N32 in C03 (N24 pointer) · N33–N42 in C04 · N43–N48 in C05 (N48 pointer) · N49–N55 in C06 · N56–N67 in C07 (N59 closure, N63 Fig. 9.8, N65 accuracy, N66 Example 9.1, N67 Example 9.2) · N70–N74 in C08 (N74 pointer) · N68, N69, N75 in C09 (N75 pointer) ·
  N76–N80 in C10 (N79 pointer) · N81–N82 in C11 (N82 pointer) and N83–N86 in the §9.9 continuation (N86 pointer) · N87–N113 in C12 · N114–N124 in C13.
- DERIVATION D01–D22: D01 D02 in C01 · D03 D04 in C02 · D05 in C03 · D06 in C04 · D07 in C05 · D08 in C06 · D09 D10 D11 in C07 · D12 in C08 · D13 in C09 · D14 in C10 · D15 D16 D17 D18 in C12 · D19 D20 D21 in C13 · D22 in C14 (22, matching the curation's CORE column).
- Explainers: E1 after C02 · E2 after C04 · E3 after C05 · E4 after C07 · E5 after C11 (§9.8) · E6 after C10 · E7 after C12 · E8 after C13 · E9 after C14 (9 embeds; B1 not embedded). SKIP S01 in §9.11.
- Primers P200–P220 (21): P200 six Reynolds numbers, P201 parabolic/elliptic/marching (C01) · P202 improper integral of a deficit, P203 simpson and trapezoid, P204 PCHIP, P205 CV with a streamline side (C02) · P206 chain rule with η(x, y) (C03) · P207 Töpfer scaling symmetry, P208 shooting vs BVP (C04) · P209 continuation and fold (C05) · P210 integrating factor, P211 ∫sin⁵ (C07) ·
  P212 inflection point (C08) · P213 row of vortices, P214 linear stability (C10) · P215 sech and arccosh, P216 total-derivative move, P217 momentum and mass flux (C12) · P218 integration by parts with a variable limit, P219 partial fractions and implicit inverse (C13) · P220 radial force balance in a swirl (C14).

---
## Part B — explainer storyboards

Common to all ten: created with `tools/new_viz.py`; `<meta name="viz:chapter" content="ch09">`; tabs Walkthrough · Explore · Explain · Derivation · Equations · Code · Check; every displayed number is computed by a JS function that mirrors a `ch09` callable and is proved by `selftest()` parity rows (`py:` expressions use only `ch09.…`, `np.pi`, numbers, strings, lists, dicts and keywords, indexed down to one float — Part C convention 3). Explain is
"Explanation & interpretation" in numbered sections built with `Viz.work.step / line / box / table / hint / interpret`, modelled on `forced_damped_vibrations.html` and `amplitude_phase_second_order_II_3.html`: **0** what the views show and what each colour means (two sentences on phones) · **1…n** every displayed quantity from the controls ("formula = substituted = result — why", results boxed) · a section per view hidden on
phones, or a hint · the values at the current time (live) · **Reading the current setting** (regime-dependent). Derivation steps are copied from Part F (same `did` titles, same step count; phones shorten *why* to its first sentence; plain-text *why* and *watch* never contain raw TeX). Every tour, Explain, Derivation, notes, status, equation and quiz text that names a book equation **writes it out** next to its number
(`tools/eq_refs.py` → 0; no bare numbers in `<meta>` strings). Drafts below write equations in Unicode for readability; builders set each in TeX (backslashes doubled in JS strings). Colours as in convention 6 (outer flow teal, profile blue, δ* amber, θ orange, δ₉₉ muted dashed, viscous rose, pressure orange, similarity/rescaled purple dashed, analytic ghost muted dashed, numerical dots teal, separation red/amber, printed-slip
muted dotted). Walkthrough texts ≤ 45 words, step 1 ≤ 24 words, ≤ 2 extras per step, `play: false` on steps that quote numbers, `autoplay: false` in the app config. A view hidden on portrait phones never carries a step's key number (repeat it in a visible title or readout). Lengths print with a unit switch (µm/mm/cm/m) at 3–4 s.f.; numbers that coincide with printed book values (0.332, 0.664, 1.72, 1.328, 4.91) are shown from the function at 4 s.f.; regime
thresholds and angles that are the book's rounded experimental values (Re ≈ 40, 200, 3×10⁵, 5×10⁵, 82°, 125°) are tagged "book's rounded values" wherever displayed and never enter a parity row. Each explainer fits 360×640 … 1920×1080 and the 1000×700 notebook frame with no scrolling (fit plan per explainer). **Data embedded in JS** (labelled "ours, computed by fluidpy"): `FS_TAB` (E3, E4; from `BL.falkner_skan_table`, 41
exponents × 161 η-points, 6 s.f.), `CLOSURE` (E4; from `BL.thwaites_closure_table`, 60 rows). Parallel builders use private scratch subfolders (`<scratchpad>/<slug>/`). **Library note:** no explainer needs `Viz.num.erf`; E2 needs a fixed-step RK4 (`Viz.num.rk4Step`, h = 0.005) and E6 a complex 4 × 4 eigen-solver (characteristic polynomial by Faddeev–LeVerrier + Durand–Kerner; ~40 lines, local).

### E1 · bl_scaling_thicknesses
- **Title:** "How thin is a boundary layer — and what is its thickness?" · **Summary:** "Slide U, ν and x to watch a layer grow like √(νx/U) = x·Re^(−1/2), then switch the profile shape: δ₉₉, δ* and θ are three integrals of one profile and change by different ratios." · **CORE:** C01, C02 (also N04–N08, N10, N14, N16, N18–N21) · **Reference:** `angular_frequency_explorer_1.html` (linked views on one clock, modes, presets, "right now" notes) and `amplitude_phase_second_order_II_3.html` (numbered Explain with live numbers).
- **meta:** `viz:order 1` · `viz:sections 9.1 9.2` · `viz:equations 9.4 9.6 9.9 9.10 9.11 9.16 9.17` · `viz:fluidpy ch09.boundary_layer_scales ch09.thicknesses ch09.blasius_fields ch09.blasius_constants ch09.profile_shape ch09.blasius_delta99 ch09.blasius_delta_star ch09.blasius_theta` · `viz:derivations D01 D03 D04`.
- **Physics (JS ↔ Python):** `scales(s)` → {Re, delta_over_L, delta, v_scale, tau0_scale, cf_estimate} ↔ `ch09.boundary_layer_scales(U, L=x, nu, rho)`; `blasius()` — Töpfer RK4 (f‴ = −½ff″ from f″(0) = 1, h = 0.005, η ≤ 12, rescale by λ = f′(∞)^{−1/2}) giving f, f′, f″ tables ↔ `ch09.blasius_constants` and `ch09.blasius_fields`; `shape(name, eta)` closed-form profiles (`linear`, `sine`, `cubic`, `exponential`, `power` with p) ↔ `ch09.profile_shape`; `thick(u, y, Ue)` composite Simpson for δ*, θ and interpolated δ₉₉ (with exponential tail) ↔ `ch09.thicknesses`; `bl99(x, U, ν)`, `blstar`, `bltheta` = coefficient × √(νx/U) ↔ `ch09.blasius_delta99` etc.
- **Views** (rows [1.15, 1]): 1. `plate` "The plate" (row 0, flex 1.4) — x ∈ [0, 2] m, y in mm auto-scaled with the exaggeration printed in the title ("y ×0.2 of x"); the δ₉₉ curve (muted dashed), the amber line y = δ*(x), ideal streamlines (teal dashed) and real streamlines (solid) lifted by δ*, profile arrows (blue) at the station x, the plate as a hatched strip; title "Re_x 6.7e4 · δ₉₉ 19.0 mm · δ*/δ₉₉ 0.35". 2. `profile` "u/U at the station" (row 0, flex 1) — u/U (x) against y (mm, or y/√(νx/U) in *scaled* mode) with the profile (blue), the deficit area 1 − u/U shaded amber (its area is δ*), the momentum-loss area (u/U)(1 − u/U) shaded orange (area θ), horizontal ticks for δ₉₉/δ*/θ; pointer: click → probe height (inspector). 3. `growth` "δ₉₉, δ*, θ against x" (row 1, `hidePortrait: true`) — log–log, x ∈ [0.01, 2] m; three straight lines of slope ½ (muted slope-½ ghost); a dot on each at the station; on portrait the numbers δ₉₉ and Re_x are repeated in the `plate` title.
- **Controls (≤ 5 visible):** `U` "Free-stream speed $U$" 0.1 … 50 m/s, step 0.1, default 1 · `fluid` chips air (ν = 1.5×10⁻⁵) / water (10⁻⁶) / oil (10⁻⁴) · `x` "Station $x$" 0.01 … 2 m, step 0.01, default 1 (the transport parameter) · `shape` chips Blasius / linear / sine / cubic / exponential / power-law p = 7 (each scaled to the same δ₉₉ as Blasius, so only the shape differs) · `scale` toggle "scaled by √(νx/U)" (mode raw / scaled; optional). Help lines: "smaller ν or bigger U thins the layer as √(ν/U)".
- **Transport:** `x` 0.01 → 2 m, rate 0.2 m/s, `end: 'loop'` (the station sweeps downstream; profile arrows and the lifted streamlines follow).
- **Presets:** "air 1 m/s at 1 m" {fluid: air, U: 1, x: 1} · "water 1 m/s at 1 m" {fluid: water, U: 1, x: 1} · "layer not thin (Re_x = 10)" {fluid: oil, U: 0.1, x: 0.01} · "fast plate (Re_x ≈ 10⁶)" {fluid: air, U: 15, x: 1} · "turbulent-like profile" {shape: power-law p = 7}.
- **Status:** "🟢 thin attached layer: Re_x = 6.67e4, δ₉₉/x = 0.0190" (Re_x ≥ 100 and < 5×10⁵) · "⚠️ layer not thin: Re_x = 10, δ₉₉/x = 1.55 — (9.9) needs Re_x ≫ 1" (Re_x < 100) · "🌀 Re_x = 1.0e6 ≥ 5×10⁵: expect transition (the laminar solution is only the start; Chs. 11–12)" (rule of thumb, Re_cr is an argument in fluidpy).
- **Readouts:** "Re_x" · "δ₉₉" (mm) · "δ*" (mm) · "θ" (mm) · "H = δ*/θ".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (raw / scaled) · **presets** · **status** · **transport** · **inspector** (click a height: "u = U f′(η) = 1 × 0.6298 = 0.630 m/s at y = 5.5 mm, η = y/√(νx/U) = 2.00") · **notes** (regime text below).
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**The plate**: blue arrows are the velocity profile at your station; the grey dashed curve is δ₉₉ (99 % of the outer speed); the amber line is δ*, the height by which the layer lifts the outer streamlines (teal dashed = ideal, solid = real). **u/U**: blue is the profile; the amber area is δ* (missing speed), the orange area is θ (missing momentum). **Growth** (hidden on phones): the three thicknesses on log axes — parallel lines of slope ½."
  1. *How big is the layer?* — "Re_x = Ux/ν = **U** × **x** / **ν** = **Re**; (9.4) δ̄/x ~ Re^(−1/2) = **1/√Re** ⇒ δ̄ ~ **δ̄** mm; cross-flow v ~ U Re^(−1/2) = **v** mm/s; dropped/kept diffusion = 1/Re = **1/Re**" boxed; why: (9.4) $\bar\delta/L\sim Re^{-1/2}$ balances u ∂u/∂x ~ U²/x against ν u_yy ~ νU/δ̄²; (9.6)–(9.7) say the x-diffusion term is 1/Re smaller than the y-diffusion term.
  2. *Three thicknesses of your profile* — "√(νx/U) = **s** mm. δ₉₉ = **k99** × s = **d99** mm · δ* = ∫(1 − u/U)dy (9.16) = **k*** × s = **dstar** mm · θ = ∫(u/U)(1 − u/U)dy (9.17) = **kθ** × s = **theta** mm · H = δ*/θ = **H**" boxed (for the Blasius shape the coefficients are 4.910, 1.721, 0.6641; for other shapes they come from Simpson integration of the drawn curve); why: three different integrals of one profile.
  3. *Where the outer streamlines go* — "v_∞ = U dδ*/dx = 0.8604 U/√Re_x = **v∞** mm/s: the outer streamlines rise with slope dδ*/dx = **slope**; over the whole plate they rise by δ*(x) = **dstar** mm" (Blasius shape only; other shapes: 'the layer here is only sketched, the lift uses your δ*'). why: continuity (6.2) with ∂u/∂x < 0 in the layer (D03).
  4. *A crude wall stress against Blasius* — "μU/δ̄ = **tau0crude** Pa (N16); Blasius τ₀ = 0.332 ρU²/√Re_x (9.31) = **tau0** Pa; ratio **r** ≈ 3" boxed; why: the crude estimate gets the Re-dependence right and the factor wrong.
  5. *At the current station* — "x = **x** m: δ₉₉ = **live** mm, δ* = **live** mm, θ = **live** mm (scaled: they are 4.910, 1.721, 0.6641 in units of √(νx/U))" (live).
  6. *Reading the current setting* — thin: "The layer is a sliver: δ₉₉/x = **ratio**. Moving x doubles δ by √2 only; the layer is thin *because* Re is large. All three thicknesses grow ∝ √x and their ratios never change for Blasius (H = 2.59)." · not thin: "At Re_x = **Re** the layer is as thick as the distance from the leading edge: the layer equations (9.9)–(9.10) are not valid here (N17), nor is the outer-flow pressure." · transition: "Above Re_x ≈ 5×10⁵ a real layer starts to become turbulent (C09): this laminar picture stops being the whole story." · shape-specific: "**linear**: H = 3, δ* is half of δ₉₉; **exponential**: H = 2; **power-law p = 7** (like a turbulent layer): H = 1.29 — a fuller profile has δ* and θ much smaller than δ₉₉."
- **Derivation tab:** **D01** (14 steps) `view: 'plate'`; goal `set` {fluid: 'air', U: 1, x: 1}; step 3 (δ̄/L ~ Re^{−1/2}) live "δ̄/L ~ 1/√(6.67e4) = 3.87e-3"; step 9 (9.7) `watch` "the dropped term's coefficient is 1/Re = 1.5e-5"; step 13 (9.10) `watch` "pressure is uniform across the layer". **D03** (6 steps) `view: 'profile'`; step 3 `highlight: ['area:deficit']`, step 6 live "v∞ = U dδ*/dx = 3.33 mm/s". **D04** (9 steps) `view: 'profile'`; step 8 `highlight: ['area:loss']`, step 9 live "ρU²θ = 1.2 × 2.572e-3 = 3.09e-3 N/m = ∫τ₀dx" (Blasius). Interpret: `s => "With your numbers Re_x = " + Re + ": the layer is " + ratio + " of x thick; the dropped x-diffusion term is " + inv + " times the kept y-diffusion term."`.
- **Code:**
  ```python
  U, nu, x = {{U}}, {{nu}}, {{x}}                 # speed [m/s], viscosity [m^2/s], station [m]
  sc = ch09.boundary_layer_scales(U, x, nu)       # (9.4): Re = {{Re}}, delta/L = {{dl}}
  y = np.linspace(0, 8*sc["delta"], 801)          # heights [m], well beyond the layer
  u = ch09.blasius_fields(x, y, U, nu)["u"]       # Blasius profile (C03): u(y) [m/s]
  th = ch09.thicknesses(y, u, U)                  # (9.16), (9.17): three integrals of that profile
  print(th["delta99"], th["delta_star"])          # {{d99}} m, {{dstar}} m
  print(th["theta"], th["H"])                     # {{theta}} m, H = {{H}}
  ```
- **Walkthrough (6 steps):** 1. "How thin?" — "Blow air over a board. Where does the wall's influence end? Press ▶ and watch the layer grow along the plate." `set` air preset, `play: true` · 2. "Balance decides" — "Advection U²/x against viscosity νU/δ²: equating gives δ ~ x/√Re (9.4). Slide U: doubling it thins the layer by √2." `controls: ['U']`, `readouts: ['Re']`, `derive: {id: 'D01', step: 3}` · 3. "Which thickness?" — "The profile fades smoothly. δ₉₉ is the 99 % height; δ* is the amber area; θ the orange area. Click the profile to probe it." `inspect: true` · 4. "δ* lifts the flow" — "Real streamlines rise by δ* above the ideal ones: the layer acts like a slightly thicker body." `derive: {id: 'D03', step: 6}`, `eq: 'dstar'` · 5. "Shape matters" — "Switch to 'linear': δ₉₉ stays the same but δ* and θ change by different ratios: H = 3." `set` {shape: 'linear'}, `controls: ['shape']`, `readouts: ['H']` · 6. "Your turn" — "Predict H for the exponential profile, then check. Then press 'layer not thin'." `controls: ['shape', 'x']`.
- **Equations:** `bl` "Layer thickness" ref 'Eq. (9.4)' `\bar\delta/L\sim Re^{-1/2}` live "δ̄/x = …" · `dstar` "Displacement thickness" ref 'Eq. (9.16)' `\delta^*=\int_0^\infty\Big(1-\frac u{U_e}\Big)dy` live "δ* = … mm" · `theta` "Momentum thickness" ref 'Eq. (9.17)' `\theta=\int_0^\infty\frac u{U_e}\Big(1-\frac u{U_e}\Big)dy` live · `d99` "Blasius δ₉₉" ref 'Eq. (9.30)' `\delta_{99}=4.91\sqrt{\nu x/U}` (note: the book prints 4.93; the root of f′ = 0.99 is 4.910) · `match` "Outer matching" ref 'Eq. (9.11)' `-\frac1\rho\frac{dp}{dx}=U_e\frac{dU_e}{dx}` · symbols U, ν, x, δ̄, δ*, θ, H with units.
- **Check yourself:** (1) "Quadruple x at fixed U and ν: what happens to δ₉₉ and to δ₉₉/x?" — "δ₉₉ doubles (√4); δ₉₉/x halves (Re_x quadruples)." `set {x: 0.25}` then x: 1 · (2) "Switch the shape to exponential and read H." — "2.0: δ* = a and θ = a/2." `set {shape: 'exponential'}` · (3) "Water instead of air at the same U and x: how does δ change?" — "Thinner by √15 = 3.9: ν is 15 times smaller." `set {fluid: 'water'}` · (4) "Where is the layer *not* thin?" — "At Re_x ≲ 100 near the leading edge or for slow viscous fluids: press 'layer not thin'."
- **Selftest parity rows:** `{name: 'Re', js: scales(S0).Re, py: 'ch09.boundary_layer_scales(1.0, 1.0, 1.5e-5)["Re"]', rtol: 1e-12}` · `{name: 'delta scale', js: scales(S0).delta, py: 'ch09.boundary_layer_scales(1.0, 1.0, 1.5e-5)["delta"]', rtol: 1e-12}` · `{name: 'delta99 Blasius', js: bl99(1, 1, 1.5e-5), py: 'ch09.blasius_delta99(1.0, 1.0, 1.5e-5)', rtol: 1e-7}` · `{name: 'delta* Blasius', js: blstar(1, 1, 1.5e-5), py: 'ch09.blasius_delta_star(1.0, 1.0, 1.5e-5)', rtol: 1e-7}` · `{name: 'theta Blasius', js: bltheta(...), py: 'ch09.blasius_theta(1.0, 1.0, 1.5e-5)', rtol: 1e-7}` · `{name: 'H sine', js: thick(sine).H, py: 'ch09.profile_shape("sine", [0.0, 1.0], 1.0)["H"]', rtol: 1e-9}` · `{name: 'H power p=7', js: (7+2)/7, py: 'ch09.profile_shape("power", [0.0, 1.0], 1.0, p=7.0)["H"]', rtol: 1e-12}` (S0 = {U: 1, x: 1, ν: 1.5e-5}).
- **Fit plan:** 360×640: status one line, `plate` (45 %) above `profile` (55 %), `growth` hidden (δ₉₉ and Re_x in the `plate` title), presets as a wrapping chip row, walkthrough card paged. 844×390: `plate` | `profile` side by side, second row hidden by scoped CSS. 1000×700 and desktop: rows [1.15, 1] with `growth` under `profile`.

### E2 · blasius_similarity_collapse
- **Title:** "Why does one curve describe the whole plate?" · **Summary:** "Profiles at different distances are one curve stretched in y: rescale by √(νx/U) and they collapse onto f′(η), whose wall slope 0.332 fixes the shear, drag and momentum loss." · **CORE:** C03, C04 (also N22–N32, N33–N42, R06, R07) · **Reference:** `amplitude_phase_second_order_II_3.html` (Explain as numbered derivation with live numbers; crosshair readouts) and `forced_damped_vibrations.html` (regime status + interpretation).
- **meta:** `viz:order 2` · `viz:sections 9.3` · `viz:equations 9.19 9.26 9.27 9.28 9.29 9.31 9.32 9.33` · `viz:fluidpy ch09.blasius_constants ch09.blasius_fields ch09.falkner_skan ch09.blasius_wall_shear ch09.blasius_skin_friction ch09.blasius_drag ch09.blasius_drag_coefficient` · `viz:derivations D05 D06`.
- **Physics (JS ↔ Python):** `ivp(s, etaMax)` RK4 of f‴ = −½ff″ with f(0) = f′(0) = 0, f″(0) = s (h = 0.005) returning f, f′, f″ tables; `topfer()` → {fpp0, eta99, delta_star, theta, H, vinf} ↔ `ch09.blasius_constants`; `fields(x, y, U, ν)` ↔ `ch09.blasius_fields`; `tau0(x)`, `Cf(Re_x)`, `FD(L)`, `CD(Re_L, sides)` ↔ `ch09.blasius_wall_shear`, `blasius_skin_friction`, `blasius_drag`, `blasius_drag_coefficient`; `shootStatus(s, etaMax)` → 'small' / 'large' / 'ok' from f′(η_max) − 1 (and `maxfp`) ↔ `ch09.falkner_skan(0.0, method="shoot")` branch logic.
- **Views** (rows [1.1, 1]): 1. `raw` "u(y) at three stations" (row 0, flex 1.3) — u [m/s] (x) against y [mm]; three curves (blue, darker with x) at x₁ = x/10, x₂ = x/2, x₃ = x, each with a small ghost of the outer flow (u = U line) and its δ₉₉ tick; hidden in mode *rescaled*. 2. `eta` "u/U against η = y√(U/νx)" (row 0, flex 1) — the collapsed f′(η) (blue) with three coloured dots sets (one per station, riding the curve), η₉₉ = 4.910 (muted), δ*/√(νx/U) = 1.721 (amber arrow), θ/√(νx/U) = 0.6641 (orange arrow); the *shoot* overlay: the IVP profile for the guessed f″(0) (rose dashed) with its f′(η_max) printed; pointer: click → probe (inspector). 3. `shear` "f″(η) and τ₀(x)" (row 1, `hidePortrait: true`) — left half f″(η) with the wall tangent slope 0.332; right half τ₀(x) = 0.332ρU²/√Re_x (Pa) with the three stations as dots.
- **Controls (≤ 5 visible):** `x` "Far station $x$" 0.05 … 2 m, step 0.01, default 1 (stations at x/10, x/2, x) · `U` "Speed $U$" 0.1 … 20 m/s, default 1 · `mode` chips raw / rescaled / both · `etaMax` "Truncation $\eta_{max}$" 3 … 16, step 0.5, default 12 (optional) · `guess` "Shooting guess for $f''(0)$" 0.20 … 0.50, step 0.001, default 0.332 (optional; the preset 'converged' sets it to the solved 0.33206) · `fluid` chips air / water (optional).
- **Transport:** `x` 0.05 → 2 m, rate 0.25 m/s, `end: 'loop'`.
- **Presets:** "air 1 m/s at 1 m" {fluid: air, U: 1, x: 1} · "thin water layer" {fluid: water, U: 5, x: 0.2} · "truncate too early" {etaMax: 3.5} · "guess too small" {guess: 0.30} · "guess too large" {guess: 0.36} · "converged" {guess: 0.33206}.
- **Status:** "✅ converged: f′(η_max) = 1.0000, f″(0) = 0.3321" · "🔻 guess 0.300 too small: f′ levels off at 0.93 — the layer never reaches the free stream" · "🔺 guess 0.360 too large: f′ overshoots and levels off at 1.06" · "⚠️ η_max = 3.5 truncates the layer: f′(3.5) = 0.91, f″(0) off by the error printed live".
- **Readouts:** "f″(0)" · "η₉₉" · "τ₀" (Pa) · "C_f at x" · "F_D per width" (mN/m).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (raw / rescaled / both) · **presets** · **status** · **inspector** (click the profile: η = y√(U/νx) = 2.00, f′(η) = 0.6298, u = U f′ = 0.630 m/s, y = 5.48 mm at x = 0.5 m) · **notes** ("Right now" table).
- **Explain:**
  0. *What the views show* — "**u(y)**: three profiles of the same layer at three distances, in metres — different curves. **u/U against η**: the same data with y divided by √(νx/U) — all on one curve f′(η) (blue). **f″ and τ₀**: the slope of the profile at the wall and the wall stress that follows from it."
  1. *The similarity variable* — "δ(x) = √(νx/U) (9.26): at x₁ **x1** m it is **d1** mm, at x₂ **d2** mm, at x₃ **d3** mm; η = y/δ (9.19) $\psi=U\delta f(\eta)$. A point at y = 5 mm sits at η = **e1**, **e2**, **e3** on the three curves" — why: no length in the problem except x, so y can appear only as y/δ(x) (D05).
  2. *Solving for f* — "Shoot with f″(0) = 1: f′(∞) = **fpinf**. Rescale f = λf₁(λη) with λ² f′(∞) = 1: λ = **lam**, f″(0) = λ³ = **fpp0**" boxed (Töpfer: f‴ + ½ff″ = 0 (9.27) is invariant under this stretch; D06 steps 1–3). Then your guess **guess** gives f′(η_max) = **fpEnd** → **verdict**.
  3. *Thicknesses in units of δ and in mm* — "η₉₉ = **4.910** (root of f′ = 0.99; the book prints 4.93), δ* = 1.721, θ = 0.6641 = 2f″(0). At x₃: δ₉₉ = **d99** mm, δ* = **dstar** mm, θ = **theta** mm" boxed.
  4. *Wall shear and drag* — "τ₀ = μU f″(0)/δ = 0.332 ρU²/√Re_x (9.31) = **tau0** Pa; C_f = 0.664/√Re_x (9.32) = **cf**; F_D = ∫τ₀dx = 0.664ρU²L/√Re_L = **fd** N/m (one side); C_D = 1.328/√Re_L (9.33) = **cd** = 2 C_f(L)" boxed; hint "F_D ∝ U^{3/2}: double U and the drag rises by 2.83".
  5. *Truncation and shooting* — "Cutting the layer at η_max = **em** changes f″(0) by **err** (from f′ → 1 at ∞); the guess **guess** misses f′(∞) = 1 by **miss**: the solution is only right for f″(0) = 0.33206."
  6. *At the current time (station x)* — "x = **x** m: Re_x = **Rex**, δ₉₉ = **d99** mm, τ₀ = **tau0** Pa (live)."
  7. *Reading the current setting* — converged/rescaled: "All three profiles lie on one curve: the flow is *similar* — profiles differ only by the √x stretch (× **s2** and × **s3** relative to the first). One number, f″(0) = 0.3321, gives the wall shear at every x." · raw: "In metres the profiles look unrelated; the layer is **r** times thicker at x₃ than at x₁." · too small/large: "A wrong f″(0) does not satisfy f′ → 1: the third-order ODE has exactly one attached solution (C04)." · truncated: "Too short a domain forces f′ = 1 at η_max where the true profile is still 0.91."
- **Derivation tab:** **D05** (13 steps) `view: 'eta'`; goal `set` {mode: 'both'}; step 2 (ansatz) `watch` "the dashed lines η = const are parabolas y ∝ √x"; step 8 (cancellation) `watch` "the two ηf′f″ terms cancel — only f f″ survives"; step 11 live "δ = √(νx/U) = **d3** mm at x₃"; step 13 `set` {mode: 'rescaled'}. **D06** (10 steps) `view: 'eta'`; step 3 live "f″(0) = λ³ = 0.6925³ = 0.3321", `set` preset 'converged'; step 4 `highlight: ['mark:eta99']`; step 9 live "F_D = 2τ₀(L)L = **fd** N/m". Interpret: `s => "With your U and x: τ₀ = " + tau0 + " Pa, δ₉₉ = " + d99 + " mm."`.
- **Code:**
  ```python
  s = {{guess}}                                      # guess for f''(0)
  rhs = lambda e, y: [y[1], y[2], -0.5*y[0]*y[2]]    # f''' = -f f''/2, (9.27)
  sol = solve_ivp(rhs, [0, 12], [0, 0, 1], rtol=1e-10)    # shoot with f''(0) = 1 (Töpfer)
  lam = sol.y[1, -1]**-0.5                           # rescale factor: lam^2 f'(inf) = 1, lam = {{lam}}
  print(lam**3)                                      # f''(0) = {{fpp0}}: the one number
  bc = ch09.blasius_constants()                      # tested library route
  print(bc["fpp0"], bc["eta99"], bc["theta"])        # {{fpp0}}, {{eta99}}, {{theta}}
  tau0 = ch09.blasius_wall_shear({{x}}, {{U}}, {{rho}}, {{nu}})   # (9.31): {{tau0}} Pa
  ```
- **Walkthrough (7 steps):** 1. "Three stations, three curves" — "Air at 1 m/s. The layer is thicker further downstream. Are these three different profiles?" `set` air preset, mode raw · 2. "No length but x" — "Nothing in the problem has a length except x, so y can only appear as η = y/δ(x). Switch to 'rescaled'." `set` {mode: 'rescaled'}, `derive: {id: 'D05', step: 2}` · 3. "One ODE" — "Insert ψ = Uδ f(η) into (9.18): two terms cancel, the rest is f‴ + ½ff″ = 0 (9.27)." `eq: 'blasius'`, `derive: {id: 'D05', step: 12}` · 4. "Shoot" — "f″(0) is unknown. Try 0.30 and 0.36: the profile never reaches 1 or overshoots. Press 'converged'." `set` guess too small, `controls: ['guess']` · 5. "One number" — "f″(0) = 0.3321 gives δ*, θ = 2f″(0), τ₀ and C_D — read them on the right." `readouts: ['fpp0', 'tau0']`, `derive: {id: 'D06', step: 6}` · 6. "Truncation" — "Set η_max = 3.5: the layer is cut, the shear is wrong; η_max ≥ 8 is safe." `set` {etaMax: 3.5}, `controls: ['etaMax']` · 7. "Your turn" — "Predict τ₀ at x = 0.2 m in water at 5 m/s, then press 'thin water layer'." `controls: ['x', 'U']`.
- **Equations:** `psi` "Similarity ansatz" ref 'Eq. (9.19)' `\psi=U\delta(x)f(\eta),\ \eta=y/\delta(x)` · `delta` "Layer scale" ref 'Eq. (9.26)' `\delta(x)=\sqrt{\nu x/U}` live "δ = … mm" · `blasius` "Blasius equation" ref 'Eq. (9.27)' `f'''+\tfrac12ff''=0` with 'Eq. (9.28)' `f(0)=f'(0)=0` and 'Eq. (9.29)' `f'(\infty)=1` in the note · `tau` "Wall shear and friction" ref 'Eq. (9.31)' `\tau_0=0.332\rho U^2/\sqrt{Re_x}` · `cf` "Skin friction" ref 'Eq. (9.32)' `C_f=0.664/\sqrt{Re_x}`, and 'Eq. (9.33)' `C_D=1.33/\sqrt{Re_L}` in the note (one side) · symbols.
- **Check yourself:** (1) "Double U at fixed x: what happens to δ, to η at fixed y, and to τ₀?" — "δ shrinks by √2; η at fixed y grows by √2; τ₀ grows by 2^{3/2} (∝ U²/√U)." `set {U: 2}` · (2) "Why is f′(∞) = 2.085 after shooting with f″(0) = 1 not a problem?" — "Because f → λf(λη) is a symmetry: rescale by λ² = 1/2.085." · (3) "Set η_max = 3.5. Which quantity is wrong first?" — "f″(0): the truncation forces f′ = 1 too early." `set {etaMax: 3.5}` · (4) "Where is the wall shear infinite?" — "At the leading edge x → 0 (τ₀ ∝ x^{−1/2}); the drag stays finite because ∫x^{−1/2}dx = 2√L."
- **Selftest parity rows:** `{name: 'fpp0', js: topfer().fpp0, py: 'ch09.blasius_constants()["fpp0"]', rtol: 1e-8}` · `{name: 'eta99', js: topfer().eta99, py: 'ch09.blasius_constants()["eta99"]', rtol: 1e-7}` · `{name: 'delta*', js: topfer().delta_star, py: 'ch09.blasius_constants()["delta_star"]', rtol: 1e-7}` · `{name: 'u at eta=2', js: fields(0.5, 5.477225575e-3, 1, 1.5e-5).u, py: 'ch09.blasius_fields(0.5, 0.005477225575, 1.0, 1.5e-5)["u"]', rtol: 1e-7}` · `{name: 'tau0', js: tau0(1, 1, 1.2, 1.5e-5), py: 'ch09.blasius_wall_shear(1.0, 1.0, 1.2, 1.5e-5)', rtol: 1e-8}` · `{name: 'CD', js: CD(66666.6667, 1), py: 'ch09.blasius_drag_coefficient(66666.6667)', rtol: 1e-8}` · `{name: 'CD two sides', js: CD(66666.6667, 2), py: 'ch09.blasius_drag_coefficient(66666.6667, sides=2)', rtol: 1e-8}` · invariant `{name: 'shoot small', js: shootStatus(0.30, 12).f1 < 1, expect: 1}`.
- **Fit plan:** 360×640: status, `eta` (55 %) above `raw` (45 %) (mode *both* shows `raw` first), `shear` hidden (f″(0) and τ₀ in readouts), the guess slider visible when the shoot presets are on. 844×390: `raw` | `eta`. Desktop: rows [1.1, 1].

### E3 · falkner_skan_family
- **Title:** "How can favourable and adverse gradients be one family?" · **Summary:** "One dial n: a fuller profile and more wall shear for n > 0, Blasius at n = 0, an inflection for n < 0 and zero wall shear at n = −0.0904 — at the wall the viscous curvature equals the pressure gradient." · **CORE:** C05, C08 (also N43–N47, N70–N73, N74) · **Reference:** `forced_damped_vibrations.html` (regime status + interpretation panel) and `fid_formula_lab.html` (term bars).
- **meta:** `viz:order 3` · `viz:sections 9.4 9.7` · `viz:equations 9.34 9.35 9.36 9.51 9.52` · `viz:fluidpy ch09.falkner_skan_state ch09.falkner_skan_separation ch09.falkner_skan ch09.wall_curvature ch09.profile_inflection ch09.falkner_skan_table` · `viz:derivations D07 D12`.
- **Physics (JS ↔ Python):** `FS_TAB` (m grid 41 nodes dense near the fold, η 0 … 8 step 0.05: f′, f″, plus f″(0), I_δ, I_θ, H, λ, l, inflection η, and `fpp0_reversed` for the second branch) ↔ `ch09.falkner_skan_table`; `state(n)` monotone-cubic interpolation of the node values ↔ `ch09.falkner_skan_state`; `profile(n)` interpolates the arrays in n; `shoot(n)` RK4 shooting for n ≥ −0.05 (independent route, parity) ↔ `ch09.falkner_skan(n, method="shoot")`; `beta(n) = 2n/(n+1)`, `halfAngle(n) = π n/(n+1)`; `wallCurv(n) = −n`; `infl(prof)` sign change of f‴; `dpdxSign(n)` ↔ `ch09.wall_curvature`.
- **Views** (rows [1.2, 1]): 1. `prof` "u/U_e and its slope at the wall" (row 0, flex 1.3) — f′ against η (or the book's ½√(n+1)η in *scale* mode) with the current profile (blue), a ghost Blasius (muted), the wall tangent (black segment of slope f″(0)), the inflection dot (amber) or 'at the wall' for n = 0, optional f″ (rose) and f‴ (teal) curves by `show` chips; pointer: click → probe. 2. `curve` "f″(0) against n" (row 0, flex 1) — attached branch (blue) from n = 4 to the fold, the second branch (amber dashed, f″(0) < 0) for −0.0904 < n < 0, the fold dot (red, n = −0.09043), the current n as a dot with its f″(0) label, n = 0 and n = 1 marked. 3. `outer` "U_e = a xⁿ and dp/dx" (row 1, `hidePortrait: true`) — U_e(x)/U_e(1) on x ∈ [0.1, 2] and the pressure gradient shading (blue-green favourable for n > 0, rose adverse for n < 0).
- **Controls (≤ 5 visible):** `n` "Exponent $n$ in $U_e=ax^n$" −0.0904 … 4, step 0.001, default 1/3 · `show` chips f′ / f′ f″ / f′ f″ f‴ · `scale` toggle "book's scaled η (½√(n+1)η)" · `branch` chips attached / second (optional) · `mu` no. Readout β = 2n/(n+1) in the `curve` title.
- **Presets:** "Hiemenz (n = 1)" {n: 1} · "Blasius (n = 0)" {n: 0} · "wedge (n = 1/3)" {n: 0.3333} · "mild adverse (n = −0.05)" {n: −0.05} · "separation (n = −0.0904)" {n: −0.09043}.
- **Status:** "🟢 favourable: n = 0.333 > 0, f″(0) = 0.7574, u_yy(wall) < 0, no inflection" · "⚪ zero pressure gradient (Blasius): f″(0) = 0.3321, inflection at the wall" · "🟠 adverse: n = −0.050, f″(0) = 0.2135, inflection at η = 1.650" · "🛑 separated at the wall: f″(0) = 0.0000 at n = −0.09043 (β = −0.1988)" · second branch: "↩️ reversed flow at the wall: f″(0) = … < 0 (second branch)".
- **Terms:** "At the wall: viscous curvature $u_{yy}\delta^2/U_e=f'''(0)=-n$ (rose) against pressure $\frac{dp}{dx}\frac{\delta^2}{\mu U_e}=-n$ (orange)" — two bars of equal length (their difference, 0, printed); click a bar to isolate it in the profile view (draws the curvature arrow).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **presets** · **status** · **terms** · **notes** ("Right now" table, below) · **modes** (scale) · **inspector**.
- **Explain:**
  0. *What the views show* — "**Profile**: blue = f′(η) = u/U_e of the layer under U_e = a xⁿ; the black segment is the wall slope (∝ shear); the amber dot is the inflection point. **f″(0) against n**: the shear as a function of the dial; it reaches zero at the red fold. **Outer flow** (hidden on phones): U_e(x) and whether the pressure gradient is favourable or adverse."
  1. *The dial* — "n = **n**, β = 2n/(n+1) = **beta**, wedge half-angle πn/(n+1) = **ang**°; U_e = a xⁿ so −dp/dx = U_eU_e′ = n a²x^{2n−1} (9.35): **verdict**" boxed.
  2. *The similarity solution* — "f‴ + ((n+1)/2) f f″ − n f′² + n = 0 (9.36) with f(0) = f′(0) = 0, f′(∞) = 1 gives f″(0) = **fpp0**, I_δ = **Id**, I_θ = **It**, H = **H**" boxed (table interpolation; parity with a JS shoot at your n when n ≥ −0.05).
  3. *The wall balance* — "Set η = 0 in (9.36): f‴(0) = −n = **fppp**. So u_yy(wall) = −n U_e/δ² and μ u_yy = dp/dx: bars **rose** = **orange** (difference 0)" boxed; hint "this is the wall curvature relation of C08."
  4. *Inflection* — "u_yy changes sign at η_i = **infl** (n < 0); none for n > 0; at the wall for n = 0 (f‴(0) = 0)" boxed.
  5. *Thickness* — "δ(x) = √(νx/U_e) = √(ν x^{1−n}/a): at x = 0.1, 0.5, 1 m: **d1**, **d2**, **d3** mm (constant for n = 1)" (table with the current row highlighted).
  6. *Right now* — table of the six presets with current row lit (n, f″(0), H, inflection, verdict).
  7. *Reading the current setting* — favourable: "The outer flow accelerates, the pressure falls along the wall and squeezes the layer: a full profile, large shear, no inflection." · Blasius: "No pressure gradient: the profile is curved only by its own diffusion; u_yy(wall) = 0." · adverse: "The rising pressure decelerates the fluid nearest the wall first: the profile flattens at the wall (smaller f″(0)) and bulges (inflection). Adverse gradients thicken the layer." · separation: "f″(0) = 0: the wall shear vanishes at every x along the wedge — separation imminent all along the surface. Below this n there is no attached solution (the book's remark that solutions exist with reverse flow refers to the second branch)." · second branch: "This solution has reverse flow at the wall; f′ → 1 still, but it is not the attached layer."
- **Derivation tab:** **D07** (12 steps) `view: 'prof'`; goal `set` {n: 1}; step 9 (v u_y) `watch` "the η terms in steps 8 and 9 have opposite signs and cancel"; step 12 result: `set` {n: 0.3333} live "f″(0) = **fpp0**". **D12** (6 steps) `view: 'prof'`; step 4 `set` {n: 1} `watch` "u_yy < 0 everywhere: no inflection"; step 5 `set` {n: −0.05} `highlight: ['mark:inflection']`; step 6 `set` preset separation, live "f″(0) = 0.0048 ≈ 0". Interpret: `s => "At n = " + n + " the wall curvature is " + (-n) + " in units of U_e/δ² — " + verdict`.
- **Code:**
  ```python
  n = {{n}}                                       # wedge exponent, U_e = a x^n
  st = ch09.falkner_skan_state(n)                 # (9.36) solved by solve_bvp with continuation
  print(st["fpp0"], st["H"])                      # wall shear f''(0) = {{fpp0}}, shape factor {{H}}
  print(st["fppp0"])                              # f'''(0) = -n = {{fppp}}: the wall curvature
  print(st["inflection_eta"], st["separated"])    # inflection at eta = {{infl}}, separated: {{sep}}
  sep = ch09.falkner_skan_separation()            # the fold of the attached branch
  print(sep["m_sep"], sep["beta_sep"])            # -0.09043, -0.19884
  ```
- **Walkthrough (6 steps):** 1. "A dial for the pressure gradient" — "Under U_e = a xⁿ one ODE describes every layer. Drag n and watch the profile." `set` {n: 1}, `controls: ['n']` · 2. "Blasius is n = 0" — "No pressure gradient: press Blasius and see the wall tangent." `set` {n: 0} · 3. "Adverse" — "Make n negative. The rising pressure flattens the profile at the wall and an amber inflection appears." `set` {n: −0.05}, `derive: {id: 'D12', step: 5}` · 4. "The wall balance" — "Set η = 0 in (9.36): f‴(0) = −n, so μ u_yy = dp/dx. The two bars are always equal." `terms: true`, `derive: {id: 'D12', step: 2}`, `eq: 'wall'` · 5. "Separation" — "At n = −0.0904 the slope at the wall is zero: no shear, separation. Below it nothing attached exists." `set` preset separation, `readouts: ['fpp0']` · 6. "Your turn" — "Predict the inflection height for n = −0.08, then check." `controls: ['n']`, `inspect: true`.
- **Equations:** `ansatz` "Falkner–Skan ansatz" ref 'Eq. (9.34)' `\psi=\sqrt{\nu xU_e}\,f(\eta),\ \eta=\frac yx\sqrt{Re_x}` · `dpdx` "Pressure gradient" ref 'Eq. (9.35)' `-\frac{dp}{dx}=U_e\frac{dU_e}{dx}=na^2x^{2n-1}` · `fs` "Falkner–Skan equation" ref 'Eq. (9.36)' `f'''+\frac{n+1}2ff''-nf'^2+n=0` live "f″(0) = …" · `acc` "Accelerating stream" ref 'Eq. (9.51)' `\Big(\frac{\partial^2u}{\partial y^2}\Big)_{wall}<0` · `dec` "Decelerating stream" ref 'Eq. (9.52)' `\Big(\frac{\partial^2u}{\partial y^2}\Big)_{wall}>0` · `wall` (note) `\mu\big(\partial^2u/\partial y^2\big)_{wall}=dp/dx`.
- **Check yourself:** (1) "For which sign of n is there an inflection point?" — "n < 0 (decelerating stream); none for n > 0." `set {n: -0.05}` · (2) "What is f‴(0) for n = 1/3?" — "−1/3: read it from (9.36) at η = 0." · (3) "Why can't you set n = −0.1?" — "The attached branch ends at the fold n = −0.0904; the second branch has reverse flow." · (4) "What does the layer thickness δ(x) do at n = 1?" — "Stays constant: δ = √(ν/a) — the stagnation-point layer."
- **Selftest parity rows:** `{name: 'fpp0 n=1', js: shoot(1).fpp0, py: 'ch09.falkner_skan_state(1.0)["fpp0"]', rtol: 1e-7}` · `{name: 'fpp0 n=1/3', js: state(0.3333333333).fpp0, py: 'ch09.falkner_skan_state(0.3333333333)["fpp0"]', rtol: 1e-6}` · `{name: 'fppp0', js: -(-0.05), py: 'ch09.falkner_skan_state(-0.05)["fppp0"]', rtol: 1e-6}` · `{name: 'H n=-0.05', js: state(-0.05).H, py: 'ch09.falkner_skan_state(-0.05)["H"]', rtol: 1e-6}` · `{name: 'infl n=-0.05', js: state(-0.05).infl, py: 'ch09.falkner_skan_state(-0.05)["inflection_eta"]', rtol: 1e-3}` · `{name: 'm_sep', js: FOLD, py: 'ch09.falkner_skan_separation()["m_sep"]', rtol: 1e-4}` · `{name: 'lam n=1', js: state(1).lam, py: 'ch09.falkner_skan_state(1.0)["lam"]', rtol: 1e-6}`.
- **Fit plan:** 360×640: status, `prof` (55 %) above `curve` (45 %), `outer` hidden (verdict in the status), `show` chips wrap. 844×390: `prof` | `curve`. Desktop: rows [1.2, 1].

### E4 · thwaites_marching
- **Title:** "Can one integral predict separation?" · **Summary:** "Choose a body's outer speed: θ² follows from ∫U_e⁵dx, λ = (θ²/ν)dU_e/dx gives wall shear and shape, and separation is λ crossing its threshold — one integral, no PDE." · **CORE:** C06, C07, C08 (also N49–N55, N56–N67, N70–N73) · **Reference:** `amplitude_phase_second_order_II_3.html` (linked windows, numbered live explanation) and `fid_formula_lab.html` (term bars).
- **meta:** `viz:order 4` · `viz:sections 9.5 9.6 9.7` · `viz:equations 9.11 9.43 9.44 9.45 9.46 9.48 9.50` · `viz:fluidpy ch09.thwaites ch09.thwaites_named ch09.thwaites_l ch09.thwaites_H ch09.thwaites_L ch09.thwaites_cylinder_closed_form ch09.thwaites_cylinder_separation ch09.example_9_2` · `viz:derivations D08 D09 D10 D11`.
- **Physics (JS ↔ Python):** `outer(kind, s, x)` → U_e, U_e′ for `flat`, `diffuser` (U₁/(1 + x/L)), `cylinder` (2U sin(x/a)), `stagnation` (a x), `retarded` (U₀(1 − 0.1x/L)) ↔ `ch09.outer_flow(kind, **p)`; `march(s)` — graded grid (600 points, x = x_max (i/N)² so stagnation points are resolved), cumulative trapezoid of U_e⁵, θ² = ν(0.45 ∫U_e⁵ + θ₀²U₀⁶/ν)/U_e⁶ (9.50), λ (9.44), l and H from `closure(λ)` (PCHIP over the embedded `CLOSURE` rows or the white fit; l = 0 below λ = −0.0681 for the exact closure), τ₀ = μU_e l/θ (9.45), δ* = Hθ, first crossing x_sep by interpolation ↔ `ch09.thwaites`; `closure(λ)` ↔ `ch09.thwaites_l/H/L`; `named(kind, x, s)` scalars at one station ↔ `ch09.thwaites_named`; `cylPhiSep(lamSep)` brentq on 0.45 cos φ F(φ)/sin⁶φ ↔ `ch09.thwaites_cylinder_separation`; `momentumResidual` finite-difference check of (9.43) on the marched arrays ↔ `ch09.momentum_integral_residual`.
- **Views** (rows [1, 1.3]): 1. `body` "Body and layer" (row 0, flex 1, `hidePortrait: true`) — the diffuser (walls diverging over L), the cylinder (with φ ticks every 30° from the forward stagnation point) or the plate, the layer band (thickness ×20, amber for δ*), the separated region in red beyond x_sep, the cursor marker; 2. `ue` "U_e(x) and dp/dx" (row 0, flex 1) — U_e (teal) and the pressure-gradient sign shading (green favourable, rose adverse), red vertical line at x_sep; 3. `lam` "λ(x), θ(x), τ₀(x)" (row 1) — three stacked mini-plots sharing the x axis and the cursor line: λ with the lines −0.0681 (exact FS, red) and −0.09 (Thwaites, amber); θ (orange, mm); τ₀ (rose, Pa) with its zero; mode `lview: closure` swaps the top plot for L(λ) with the Falkner–Skan dots and the line 0.45 − 6λ (D10). Pointer: click or drag in `lam`/`ue` moves the cursor.
- **Controls (≤ 5 visible):** `flow` chips Blasius / diffuser / cylinder / stagnation / retarded · `L` "Length scale $L$" (diffuser and retarded: length over which the area doubles / the speed drops 10 %) 0.2 … 2 m, default 0.5 · `theta0` "Initial θ₀" 0 … 0.5 mm, step 0.01, default 0 · `closure` chips exact FS / Thwaites' (white) · `x` "Cursor $x$" (transport parameter) 0 … x_max · optional: `U1` "Inlet speed" 1 … 50 m/s, default 10 (only τ₀ in Pa depends on it), `lview` chips x / closure.
- **Transport:** `x` 0 → x_max, rate = x_max/12 per s, `end: 'hold'` with an end card ("separation predicted at x = 0.063 m" or "attached to the end").
- **Presets:** "Blasius (constant U_e)" {flow: flat, theta0: 0} · "diffuser (Example 9.2)" {flow: diffuser, L: 0.5, theta0: 0} · "thick inlet layer" {flow: diffuser, theta0: 0.1} · "cylinder (ideal speed)" {flow: cylinder} · "stagnation flow" {flow: stagnation} · "Howarth-type retardation" {flow: retarded}.
- **Status:** "🟢 attached: λ = −0.052 at x = 0.038 m, above −0.0681 (exact FS)" · "🛑 separation predicted at x = 0.063 m (x/L = 0.126, λ = −0.0681); beyond it the layer equations no longer apply" · with the white closure "…(λ = −0.09) at x = 0.079 m (x/L = 0.158)" · "🟢 accelerating everywhere: λ > 0" (stagnation, flat).
- **Terms:** stacked bar of θ²U_e⁶/ν at the cursor = 0.45∫₀ˣU_e⁵dx′ (teal) + θ₀²U₀⁶/ν (rose) with each share in % (θ₀ = 0 leaves one bar); click a bar to isolate its effect on θ(x).
- **Inspector:** click x in `lam`: "θ²U_e⁶/ν = 0.45 × ∫U_e⁵dx + θ₀²U₀⁶/ν = **rhs** ⇒ θ = **theta** mm; λ = (θ²/ν)U_e′ = **lam**; l(λ) = **l**, H = **H** ⇒ τ₀ = μU_e l/θ = **tau0** Pa" (every number from the current state, none typed).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **terms** · **presets** · **status** · **inspector** · **notes** (regime text).
- **Explain:**
  0. *What the views show* — "**Body and layer**: the geometry with the layer thickness (×20) and, in red, where the layer has 'let go'. **U_e(x)**: the ideal-flow speed at the edge; where it falls (adverse, rose) the pressure rises. **λ, θ, τ₀**: the pressure-gradient parameter, the momentum thickness and the wall shear along x. Orange = θ, rose = τ₀ (viscous), teal = outer flow."
  1. *The outer flow* — "U_e = **Ue** m/s, dU_e/dx = **dUe** 1/s, so dp/dx = −ρU_eU_e′ (9.11) = **dpdx** Pa/m: **favourable/adverse**" boxed.
  2. *The integral of U_e⁵* — "∫₀ˣU_e⁵dx′ = **I** (m⁶/s⁵); θ²U_e⁶/ν = 0.45 × **I** + θ₀²U₀⁶/ν = **rhs**; θ = **theta** mm (9.50)" boxed; why: (9.43) with the closure L ≈ 0.45 − 6λ is a linear ODE with integrating factor U_e⁶ (D09).
  3. *The pressure-gradient parameter* — "λ = (θ²/ν) dU_e/dx (9.44) = **theta2/nu** × **dUe** = **lam**" boxed; hint "for θ₀ = 0 the layer 'forgets' U₁, ν and L: λ depends on the shape of U_e only".
  4. *Shear, shape, thickness* — "l(**lam**) = **l**, H(**lam**) = **H** (closure **closure**); τ₀ = μ(U_e/θ)l (9.45) = **tau0** Pa; δ* = Hθ (9.46) = **dstar** mm; C_f = 2τ₀/(ρU_e²) = **cf**" boxed.
  5. *The separation test* — "λ = **lam** against the thresholds −0.0681 (the exact Falkner–Skan member where l = 0) and −0.09 (Thwaites' fit, the book's criterion): **verdict**; x_sep = **xsep** m (x/L = **xsepL**)" boxed; hint "the closure table is ours, from solved Falkner–Skan profiles; the book's Table 9.1 is not used".
  6. *How good is it?* — "Against the exact Falkner–Skan layers Thwaites' θ is +1.0 % (n = 0), −4.2 % (1/3), −6.3 % (1), −7.6 % (4) and +3.1 % (n = −0.05); it predicts *whether* a layer separates better than *where*." (static table from the notebook, our numbers).
  7. *At the cursor* — "x = **x** m: U_e = **Ue**, θ = **theta** mm, λ = **lam**, τ₀ = **tau0** Pa (live)."
  8. *Reading the current setting* — diffuser: "Straight diffuser: U_e falls like 1/(1 + x/L); the pressure rises along the wall, λ = −(0.45/4)[(1 + x/L)⁴ − 1] falls through the threshold at x/L = 0.158 (Thwaites) or 0.126 (exact). A thicker inlet layer (θ₀ > 0) separates sooner." · cylinder: "On the ideal cylinder speed the threshold is crossed at φ = **phi**° (or 100.9°), yet real cylinders separate near 82° (the book's rounded experimental value): the ideal pressure ceases to be the real one once the wake forms (C09)." · Blasius: "Constant U_e: λ = 0 exactly, θ = 0.6708√(νx/U), 1 % above Blasius." · stagnation: "Accelerating: λ = 0.075 at the stagnation point (exact Hiemenz 0.0855), τ₀ large, no separation." · retarded: "A linearly retarded outer flow separates at an x set by its slope (qualitative: the book quotes a value we have not verified)."
- **Derivation tab:** **D08** (11 steps) `view: 'ue'`; goal `set` {flow: 'flat'}; step 2 `watch` "adding u times continuity turns the advective terms into flux form"; step 11 live "τ₀/ρ = U² dθ/dx = **tau0rho** m²/s² on both sides". **D09** (10 steps) `view: 'lam'`; step 5 live "λ = θ²U_e′/ν = **lam**"; step 8 `set` {lview: 'closure'}; step 10 live "θ²U_e⁶/ν = **rhs**", `set` {lview: 'x'}. **D10** (9 steps) `view: 'lam'`, `set` {lview: 'closure'}; step 3 live "λ = n I_θ² = 1 × 0.2923² = 0.0855"; step 8 `watch` "the dots follow the line to within 0.06". **D11** (8 steps) `view: 'ue'`, `set` {flow: 'cylinder'}; step 5 live "λ(φ) = 0.45 cos φ F(φ)/sin⁶φ = **lamphi**"; step 7 live "φ_sep = 103.11° (λ = −0.09), 100.89° (−0.0681)". Interpret: `s => "With your body: separation " + (xsep ? "predicted at x = " + xsep + " m" : "not predicted along the body") + "."`.
- **Code:**
  ```python
  x = np.linspace(0, {{xmax}}, 601)                # stations along the body [m]
  Ue = ch09.outer_flow("{{flow}}", {{params}})     # outer flow object: Ue(x), dUe(x), dpdx (9.11)
  th = ch09.thwaites(x, Ue.Ue, {{nu}}, theta0={{theta0m}}, closure="{{closure}}")   # (9.50) then (9.44)-(9.46)
  print(th["theta"][-1], th["lam"][-1])            # theta = {{theta}} m, lambda = {{lam}} at the cursor
  print(th["tau0"][-1], th["H"][-1])               # wall shear {{tau0}} Pa, H = {{H}}
  print(th["separated"], th["x_sep"])              # {{sep}}, x_sep = {{xsep}} m
  print(ch09.thwaites_l({{lam}}), ch09.thwaites_H({{lam}}))   # the closure at the cursor
  ```
- **Walkthrough (7 steps):** 1. "The question" — "You know the outer speed U_e(x) of a body. Does the layer separate, and where? Press ▶ to march along a diffuser." `set` diffuser preset, `play: true` · 2. "Pressure from the outer flow" — "Where U_e falls the pressure rises (9.11): $-\frac1\rho\frac{dp}{dx}=U_e\frac{dU_e}{dx}$. Rose = adverse." `eq: 'match'`, `derive: {id: 'D08', step: 10}` · 3. "One integral" — "Thwaites' closure turns (9.43) into θ²U_e⁶/ν = 0.45∫U_e⁵dx (9.50). Drag the cursor: the teal bar is that integral." `terms: true`, `derive: {id: 'D09', step: 10}`, `eq: 'thw'` · 4. "λ tells the story" — "λ = (θ²/ν)U_e′ (9.44) falls as the flow decelerates; l(λ) and H(λ) give shear and shape." `readouts: ['lam', 'l']`, `derive: {id: 'D09', step: 6}` · 5. "Threshold" — "τ₀ = 0 when λ reaches −0.0681 (exact) or −0.09 (Thwaites): the layer separates. Switch the closure." `controls: ['closure']`, `readouts: ['xsep']` · 6. "Memory" — "Set θ₀ = 0.1 mm: a thicker inlet layer separates earlier. The rose bar is that memory." `set` thick inlet layer, `controls: ['theta0']` · 7. "Your turn" — "Choose 'cylinder' and read the predicted angle; compare with 82° for a real cylinder." `set` cylinder preset, `inspect: true`.
- **Equations:** `match` "Outer matching" ref 'Eq. (9.11)' `-\frac1\rho\frac{dp}{dx}=U_e\frac{dU_e}{dx}` · `vk` "Von Kármán momentum integral" ref 'Eq. (9.43)' `\frac1\rho\tau_0=\frac d{dx}\big[U_e^2\theta\big]+U_e\delta^*\frac{dU_e}{dx}` live "τ₀/ρ = …" · `lam` "Pressure-gradient parameter" ref 'Eq. (9.44)' `\lambda\equiv\frac{\theta^2}\nu\frac{dU_e}{dx}` and 'Eq. (9.45)' `\tau_0=\mu\frac{U_e}\theta\,l(\lambda)` in the note · `L` "Universal function" ref 'Eq. (9.48)' `L(\lambda)=2l-2(2+H)\lambda\approx0.45-6\lambda` · `thw` "Thwaites' closed form" ref 'Eq. (9.50)' `\frac{\theta^2U_e^6}\nu=0.45\int_0^xU_e^5dx'+\frac{\theta_0^2U_0^6}\nu` live "θ = … mm".
- **Check yourself:** (1) "Double the inlet speed U₁ of the diffuser (θ₀ = 0): does the separation point move?" — "No: λ(x/L) depends only on the shape of U_e; x_sep is fixed by x/L." `set {U1: 20}` · (2) "Which closure separates first, exact FS or white?" — "The exact one (λ = −0.0681 is reached before −0.09)." `set {closure: 'white'}` · (3) "What does λ equal at Blasius' plate?" — "Zero: U_e′ = 0." · (4) "Why does Thwaites on the cylinder separate at 103° when real ones separate near 82°?" — "It uses the ideal pressure; once the wake forms the real pressure differs (C09)."
- **Selftest parity rows:** `{name: 'L(0)', js: closure(0).L, py: 'ch09.thwaites_L(0.0)', rtol: 1e-6}` · `{name: 'l(0)', js: closure(0).l, py: 'ch09.thwaites_l(0.0)', rtol: 1e-6}` · `{name: 'H(-0.05)', js: closure(-0.05).H, py: 'ch09.thwaites_H(-0.05)', rtol: 1e-6}` · `{name: 'cyl lambda 82', js: cylLam(1.43117), py: 'ch09.thwaites_cylinder_closed_form(1.43117)', rtol: 1e-12}` · `{name: 'cyl phi_sep', js: cylPhiSep(-0.09), py: 'ch09.thwaites_cylinder_separation()', rtol: 1e-8}` · `{name: 'diffuser lam 0.1', js: named('diffuser', 0.05, S0).lam, py: 'ch09.thwaites_named("diffuser", 0.05, 1.5e-5, U1=10.0, L=0.5)["lam"]', rtol: 1e-4}` · `{name: 'x_sep/L', js: xsepDiffuser(), py: 'ch09.example_9_2()["x_sep_over_L"]', rtol: 1e-9}` (S0 = {U1: 10, L: 0.5, ν: 1.5e-5, θ₀: 0}).
- **Fit plan:** 360×640: status, `ue` (35 %) above `lam` (65 %), `body` hidden (the red x_sep line is drawn in `ue` and `lam`); the terms bar sits in the `lam` title; presets wrap. 844×390: `ue` | `lam`. 1000×700 and desktop: rows [1, 1.3] with `body` left of `ue`.

### E5 · cylinder_drag_crisis
- **Title:** "Why does drag fall when the flow gets faster?" · **Summary:** "A log-scale Re moves a schematic flow, the separation angle, C_p and C_D together: a laminar layer lets go near 82°, a turbulent one holds on to 125°, the wake narrows and the drag collapses." · **CORE:** C09, C10, C11 (also N68, N69, N76–N78, N80, N81, N83–N85) · **Reference:** `angular_frequency_explorer_1.html` (modes, presets, "right now" table).
- **meta:** `viz:order 5` · `viz:sections 9.7 9.8 9.9` · `viz:equations 4.102 9.33` (the pressure-drag integral of D13 is ours and has no number) · `viz:fluidpy ch09.separated_cp ch09.separated_pressure_drag ch09.cylinder_flow_regime ch09.cylinder_state ch09.shedding_frequency ch09.sphere_drag_coefficient ch09.stokes_drag_coefficient ch09.oseen_drag_coefficient` · `viz:derivations D13`.
- **Physics (JS ↔ Python):** `regime(Re, mode)` → {label, phi_sep, St} ↔ `ch09.cylinder_flow_regime`, `sphere_flow_regime`; `state(Re, rough)` → {phi_sep_deg, cb, label} (roughness lowers Re_cr by a factor 3 — schematic) ↔ `ch09.cylinder_state`; `cpSep(phi, phis, cb)` ↔ `ch09.separated_cp`; `cdp(phis, cb)` = sin φ_s (1 − (4/3) sin²φ_s − C_b) ↔ `ch09.separated_pressure_drag`; `cdQual(Re)` hand-drawn polyline (labelled qualitative) ↔ `ch09.cylinder_cd_schematic`; `sphereCD(Re)` Morrison correlation (port, 1e-10) ↔ `ch09.sphere_drag_coefficient`; `stokes(Re)`, `oseen(Re)` ↔ `ch09.stokes_drag_coefficient`, `oseen_drag_coefficient`; `shed(U, d)` ↔ `ch09.shedding_frequency`.
- **Views** (rows [1.1, 1]): 1. `flow` "The flow past the body" (row 0, flex 1.3) — schematic streamlines and wake by regime: creeping (symmetric), twin attached eddies (length grows with Re), a street of alternating vortices advected on the transport clock (60 particles), subcritical wake (wide, separation dots at 82°), supercritical (narrow wake, dots at 125°); the sphere mode draws a doughnut eddy / loops; separation points marked red. 2. `cp` "C_p(φ)" (row 0, flex 1) — ideal 1 − 4 sin²φ (black dashed), separated model (blue) with the wake plateau C_b, the region between the model and the ideal recovery shaded (the missing pressure recovery); title carries C_D,p. 3. `cd` "C_D against Re" (row 1, `hidePortrait: true`) — log–log; cylinder mode: qualitative schematic (dashed, badge "qualitative, hand-drawn"); sphere mode: Morrison (solid, badge "tested") with Stokes (blue dashed, Re < 1) and Oseen (teal dashed) lines; the current Re as a dot; the four regimes shaded.
- **Controls (≤ 5 visible):** `logRe` "Reynolds number (log₁₀ Re)" −1 … 7, step 0.02, default 2 (Re on the diameter) · `mode` chips cylinder / sphere · `surface` chips smooth / rough · `cb` "Wake pressure $C_b$" −2 … 0.5, step 0.05, default auto (optional; auto = −1.2 subcritical, −0.6 supercritical, blended across the crisis — illustrative, ours) · `d` "Diameter" 1 mm … 1 m (optional), `U` "Speed" (optional) for the frequency.
- **Transport:** `t` 0 → 20 s, rate 1 s/s, `end: 'loop'` (shedding clock, street advection).
- **Presets:** "creeping (Re 0.5)" {logRe: −0.3} · "twin eddies (Re 20)" {logRe: 1.3} · "street (Re 100)" {logRe: 2} · "subcritical (Re 10⁵)" {logRe: 5} · "critical (Re 3×10⁵)" {logRe: 5.48} · "supercritical (Re 10⁶)" {logRe: 6} · "rough cylinder at 10⁵" {logRe: 5, surface: rough}.
- **Status:** "🌀 subcritical: laminar separation at φ ≈ 82° (book's rounded value), C_b = −1.2, C_D,p = 0.884 (model)" · "⚡ critical: drag crisis — the layer turns turbulent, C_D collapsing" · "🌬️ supercritical: separation at φ ≈ 125°, C_D,p = 0.578 (model)" · "🌊 laminar street: St ≈ 0.2, f = 1000 Hz for 2 mm at 10 m/s" · "◯ creeping: symmetric, no wake".
- **Notes:** a regime table (Re range · label · separation · St) with the current row lit; header "the book's rounded thresholds, experimental — not computed".
- **Depth features:** Explain + Code + Derivation · **modes** (cylinder / sphere) · **linked views** (3) · **presets** · **status** · **transport** · **notes** (regime table) · **inspector** (click `cp`: "at φ = 60°: C_p = 1 − 4 sin²60° = −2.00 (ideal); the model follows it up to φ_s").
- **Explain:**
  0. *What the views show* — "**Flow**: a cartoon of the wake for the Reynolds number you chose (red dots = where the layer separates). **C_p(φ)**: pressure around the body from the front stagnation point (φ = 0): the black dashed line is ideal flow, blue is the separated model — flat behind the separation point; the shaded gap is the pressure the wake does not recover. **C_D**: drag against Re; the dot is you."
  1. *Reynolds number* — "Re = U d/ν = **U** × **d** / **ν** = **Re** (on the diameter, not the length or the radius); regime: **label** (book's rounded values)" boxed.
  2. *Shedding* — "St = Ωd/U ≈ 0.2 (4.102) ⇒ f = St U/d = 0.2 × **U** / **d** = **f** Hz; the *angular* frequency is 2πf = **omega** rad/s (the book's Ω is the frequency in Hz)" boxed; hint "using rad/s for Ω would give St = 1.26".
  3. *Where the layer lets go* — (for Re < 3000 the `cp` view shows only the ideal curve and this section reads "the separated model applies to high-Reynolds-number flow (Re ≥ 3000); here the wake is the regime of step 1") "φ_s = **phis**° from the forward stagnation point (laminar ≈ 82°, turbulent ≈ 125°: the book's rounded values); ideal C_p(φ_s) = 1 − 4 sin²φ_s = **cpid**; the wake stays at C_b = **cb** (illustrative)" boxed.
  4. *Pressure drag of the model* — "C_D,p = ½∮C_p cos φ dφ = sin φ_s (1 − (4/3) sin²φ_s − C_b) = **s** × (1 − **4/3 s²** − **cb**) = **cdp**" boxed; ideal limit: φ_s → 180°, C_b → 1 ⇒ 0 (d'Alembert). Why: D13.
  5. *Sphere mode* — "Morrison's C_D(Re) = **cds** (tested, V5); Stokes 24/Re = **stokes**, Oseen = 24/Re(1 + 3Re/16) = **oseen** at low Re; the dip near Re ≈ 4×10⁵ is the crisis (the book's rounded Re_cr ≈ 5×10⁵)."
  6. *Right now* — table of the seven presets with the current row lit.
  7. *Reading the current setting* — creeping: "Viscosity dominates: fore–aft symmetric flow, no wake, C_D ∝ 1/Re (Stokes/Oseen)." · twin eddies: "Two steady attached eddies grow with Re; the flow is still laminar and steady." · street: "The wake is unstable: alternating vortices shed at St ≈ 0.2 (C10); the drag is mostly form drag." · subcritical: "A laminar layer separates near 82°; the wake pressure stays low, C_D ≈ 1 over decades of Re." · critical: "The layer turns turbulent before separating: fuller profile, separation moves aft, the wake shrinks and C_D drops — the drag crisis (faster is *less* drag)." · supercritical: "Separation near 125°, narrower wake, higher wake pressure: C_D,p ≈ 0.58 (model); slowly rising again as Re grows." · rough: "Roughness (or free-stream turbulence) trips the layer early: the crisis moves to lower Re — the reason for golf-ball dimples and cricket seams."
- **Derivation tab:** **D13** (7 steps) `view: 'cp'`; goal `set` {logRe: 5}; step 3 (drag coefficient as ½∮C_p cos φ) `watch` "cos φ weights the front (push) and rear (pull) pressures with opposite signs"; step 4 `highlight: ['curve:model']`; step 7 live "C_D,p = sin φ_s(1 − 4/3 sin²φ_s − C_b) = **cdp**", `set` preset supercritical. Interpret: `s => "At Re = " + Re + " the model gives C_D,p = " + cdp + "; the ideal (φ_s = 180°) would give 0."`.
- **Code:**
  ```python
  Re, d, U = {{Re}}, {{d}}, {{U}}                    # Re = U d / nu on the diameter
  st = ch09.cylinder_state(Re, rough={{rough}})      # regime table look-up (book's rounded thresholds)
  print(st["label"], st["phi_sep_deg"], st["St"])    # {{label}}, phi_sep = {{phis}} deg from the front stagnation point
  cb = {{cb}}                                        # wake pressure (illustrative, ours)
  print(ch09.separated_pressure_drag(st["phi_sep_deg"], cb))   # C_D,p = sin(phi_s)(1 - 4/3 sin^2 phi_s - C_b) = {{cdp}}
  print(ch09.shedding_frequency(U, d))               # f = St U/d = {{f}} Hz, omega = {{omega}} rad/s
  print(ch09.sphere_drag_coefficient(Re))            # Morrison (tested): {{cds}}
  ```
- **Walkthrough (7 steps):** 1. "Faster, less drag?" — "Slide the Reynolds number from 10 to 10⁷ and name each wake." `set` {logRe: 1}, `controls: ['logRe']`, `play: true` · 2. "The street" — "At Re ≈ 100 vortices peel off alternately; St ≈ 0.2 (4.102) gives the frequency." `set` street preset, `readouts: ['f']` · 3. "Where the layer lets go" — "Laminar layers separate near 82°. Behind that the pressure stays low." `set` subcritical preset, `derive: {id: 'D13', step: 4}` · 4. "Drag is missing recovery" — "The shaded gap between real and ideal pressure is the form drag: C_D,p = sin φ_s(1 − 4/3 sin²φ_s − C_b)." `eq: 'cdp'`, `derive: {id: 'D13', step: 7}` · 5. "The crisis" — "At Re ≈ 3×10⁵ the layer turns turbulent, holds on to 125° and C_D collapses." `set` supercritical preset · 6. "Roughness" — "A rough cylinder does this at 10⁵. Toggle 'rough'." `set` rough preset, `controls: ['surface']` · 7. "Your turn" — "Switch to the sphere: find the dip and compare with Stokes at small Re." `set` {mode: sphere}, `controls: ['mode', 'logRe']`.
- **Equations:** `st` "Strouhal number" ref 'Eq. (4.102)' `St=\frac{\Omega d}{U_\infty}\approx0.2` live "f = … Hz" · `cp` "Ideal surface pressure" (ch06) `C_p=1-4\sin^2\varphi` · `cdp` "Pressure drag of the separated model (ours)" `C_{D,p}=\tfrac12\oint C_p\cos\varphi\,d\varphi=\sin\varphi_s\big(1-\tfrac43\sin^2\varphi_s-C_b\big)` live · `cd` "Sphere drag" — Stokes `C_D=24/Re`, Oseen `C_D=\frac{24}{Re}\big(1+\frac{3Re}{16}\big)` (Ch. 8) · `plate` "Flat-plate drag (laminar), one side" ref 'Eq. (9.33)' `C_D=1.33/\sqrt{Re_L}` (linked from C09).
- **Check yourself:** (1) "At φ_s = 90° with the ideal C_b = −3, what is C_D,p? With C_b = −1?" — "2.667; 0.667." `set` via presets · (2) "Why does the rough cylinder at Re = 10⁵ have lower drag than the smooth one?" — "Roughness trips the layer: separation moves from 82° to 125°." · (3) "Which is higher: f in Hz or Ω in rad/s for the same wire?" — "Ω = 2πf ≈ 6.3 f." · (4) "What is the limit φ_s → 180°?" — "C_D,p → 0 (d'Alembert's paradox recovered)."
- **Selftest parity rows:** `{name: 'cdp 82', js: cdp(82, -1.2), py: 'ch09.separated_pressure_drag(82.0, -1.2)', rtol: 1e-12}` · `{name: 'cdp 125', js: cdp(125, -0.6), py: 'ch09.separated_pressure_drag(125.0, -0.6)', rtol: 1e-12}` · `{name: 'cp model', js: cpSep(60, 82, -1.2), py: 'ch09.separated_cp(60.0, 82.0, -1.2)', rtol: 1e-12}` · `{name: 'cp wake', js: cpSep(120, 82, -1.2), py: 'ch09.separated_cp(120.0, 82.0, -1.2)', rtol: 1e-12}` · `{name: 'f shed', js: shed(10, 0.002).f, py: 'ch09.shedding_frequency(10.0, 0.002)["f"]', rtol: 1e-12}` · `{name: 'sphere 4e5', js: sphereCD(4e5), py: 'ch09.sphere_drag_coefficient(4e5)', rtol: 1e-10}` · `{name: 'Oseen', js: oseen(0.1), py: 'ch09.oseen_drag_coefficient(0.1)', rtol: 1e-12}` · exact-text `{name: 'regime label', js: regime(100).label, py: 'ch09.cylinder_flow_regime(100.0)["label"]'}` · invariant `{name: 'ideal limit', js: cdp(179.9, 1 - 4*Math.pow(Math.sin(179.9*Math.PI/180), 2)), expect: 0, atol: 1e-6}`.
- **Fit plan:** 360×640: status one line, `flow` (50 %) above `cp` (50 %), `cd` hidden (C_D in the `cp` title), presets wrap; 844×390: `flow` | `cp`; desktop: rows [1.1, 1] with `cd` under `cp`.

### E6 · karman_street_stability
- **Title:** "Why is the vortex street's shape fixed at b/a ≈ 0.28?" · **Summary:** "A staggered double row of point vortices grows a wiggle at every spacing except cosh(πb/a) = √2; drag b/a and the stagger to watch the growth rate and the four eigenvalues." · **CORE:** C10 (also N78, N80) · **Reference:** `forced_damped_vibrations.html` (a state view next to a response curve, regime interpretation).
- **meta:** `viz:order 6` · `viz:sections 9.8` · `viz:equations 4.102` (the row velocity and the marginal condition are shown as ours, unnumbered) · `viz:fluidpy ch09.karman_street_ratio ch09.karman_street_growth ch09.karman_street_growth_closed ch09.karman_street_spectrum ch09.karman_street_velocity ch09.karman_street_positions` · `viz:derivations D14`.
- **Physics (JS ↔ Python):** `matrixPi(b, off)` — the 2 × 2 complex system for the alternating mode from the closed-form lattice sums P = π²/sin²(πζ), Q = π² cos(πζ)/sin²(πζ), ζ = −o + ib (as D14 steps 7–10), returned as the 4 × 4 real matrix ↔ the k = π/a case of `ch09.karman_street_spectrum`; `spectrum(b, off, k)` general k: lattice sums by direct summation (±2000, k = 0 tail 2/N) + 4 × 4 complex eigenvalues (Faddeev–LeVerrier + Durand–Kerner) ↔ `ch09.karman_street_spectrum`; `growthClosed(b)` = (π/2)\|½ − sech²(πb)\| ↔ `ch09.karman_street_growth_closed`; `growthMax(b, off)` ↔ `ch09.karman_street_growth`; `positions(t, b, off, eps)` RK4 of the linear system for the alternating mode ↔ `ch09.karman_street_positions`; `speed(a, b, Γ)` = (Γ/2a) tanh(πb/a) ↔ `ch09.karman_street_velocity`; `ratio()` ↔ `ch09.karman_street_ratio`.
- **Views** (rows [1.1, 1]): 1. `street` "The double row" (row 0, flex 1.4) — eight cells of vortices in the street frame: row A (blue, Γ) at y = +b/2, row B (rose, −Γ) at y = −b/2 offset by o·a, counter-clockwise/clockwise arrows on each; displacements from the linear evolution of a small kick (amplitude exaggerated ×20, stated in the title); the street velocity arrow; a dashed ghost of the undisturbed street; animates with the clock. 2. `growth` "Growth rate against b/a" (row 0, flex 1) — σ (units Γ/a²) against b/a: closed-form curve (blue), the numerical max over k as dots at 20 values, the marginal spacing (green dot at 0.2805), the non-staggered value π/4 (grey dashed) and the current b/a as a moving dot. 3. `spec` "Eigenvalues in the complex plane" (row 1, `hidePortrait: true`) — the four λ (× Γ/a²) for the chosen b/a, offset and k with the axes; imaginary axis highlighted; the growth is the largest real part.
- **Controls (≤ 5 visible):** `b` "Row separation ratio $b/a$" 0.05 … 1, step 0.005, default 0.2 · `off` "Stagger offset (0 = facing rows, 0.5 = staggered)" 0 … 0.5, step 0.01, default 0.5 · `k` "Perturbation wavenumber $k\,a/\pi$" 0.05 … 1, step 0.01, default 1 (1 = alternating, the most dangerous; optional) · `kick` "Kick" button (increments a seed; restarts the growth) · `amp` "Kick size" (optional).
- **Transport:** `t` 0 → 12 (units a²/Γ), rate 1, `end: 'hold'` (end card: "grew ×… in t = 12" or "stayed at ×1.0").
- **Presets:** "stable spacing b/a = 0.2805" {b: 0.2805, off: 0.5} · "too tight (0.15)" {b: 0.15} · "too loose (0.5)" {b: 0.5} · "non-staggered" {off: 0, b: 0.28}.
- **Status:** "⚖️ marginal: growth 6e-9 Γ/a² ≈ 0 (cosh(πb/a) = √2)" · "⚠️ unstable: growth σ = 0.2982 Γ/a² (e-fold every 3.4 a²/Γ), b/a = 0.20" · "⚠️ non-staggered rows: σ = π/4 = 0.7854 Γ/a² at every b/a" · verdict for the chosen k: "this k grows at 0.23 (the maximum, at k = π/a, is 0.2982)".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** · **status** · **notes** · **inspector** (click a vortex: "vortex A₃: position (3.00, 0.14) + kick (0.03, −0.02); induced speed from row B: (Γ/2a)tanh(πb/a) = 0.3536").
- **Explain:**
  0. *What the views show* — "**Street**: two rows of point vortices (blue +Γ, rose −Γ) carried along at the street speed; the ghost is the perfect street, the dots are the kicked one. **Growth**: the fastest growth rate over all disturbances, against b/a. **Eigenvalues** (hidden on phones): the four exponents λ of the linearised motion; a positive real part is growth."
  1. *The velocity of a row* — "One vortex: w = u − iv = (Γ/2πi)/(z − z₀). A row of spacing a: Σ 1/(z − na) = (π/a)cot(πz/a) ⇒ w = (Γ/2ia)cot(π(z − z₀)/a). At your b/a = **b** the street moves at U_s = (Γ/2a)tanh(πb/a) = **Us** (Γ = a = 1)" boxed.
  2. *The lattice sums* — "For the alternating disturbance (each vortex displaced by d(−1)ⁿ): P = π²/sin²(πζ) = **P**, Q = π² cos(πζ)/sin²(πζ) = **Q** with ζ = −o + ib = **zeta**; the rows' own contributions are π²/3 and −π²/6" boxed.
  3. *The linear system* — "γ = ½ − 1/cosh²(πb) = **gamma**, σ = sinh(πb)/cosh²(πb) = **sigma**: d/dt (x_A, y_A, x_B, y_B) = (π/2)[[0, −γ, −σ, 0], [−γ, 0, 0, σ], [σ, 0, 0, γ], [0, −σ, γ, 0]] (x_A, y_A, x_B, y_B)" (staggered case; offsets ≠ 0.5 use the general P, Q).
  4. *Eigenvalues* — "λ = (π/2)(±γ ± iσ) = ±**re** ± **im** i (units Γ/a²); growth = (π/2)\|γ\| = **growth**" boxed; the factorised polynomial ((μ − γ)² + σ²)((μ + γ)² + σ²).
  5. *The marginal spacing* — "γ = 0 ⇔ cosh²(πb/a) = 2 ⇔ b/a = arccosh(√2)/π = **ratio**; there σ = ½ and λ = ±i(π/4)Γ/a²: pure oscillation, no growth; your b/a is **d** from it" boxed.
  6. *Other disturbances* — "For wavenumber ka/π = **k** the sums are truncated at ±2000: the largest real part is **gk** (the alternating mode k = π/a is the maximum for every b/a tried). For rows facing each other (o = 0) every b/a grows at π/4 = 0.7854."
  7. *At the current time* — "t = **t**: the kick has grown by e^{σt} = **amp** (largest displacement / initial) (live)."
  8. *Reading the current setting* — stable: "At cosh(πb/a) = √2 the staggered street neither grows nor decays (marginal): this is the shape nature shows behind cylinders (b/a ≈ 0.28). Real streets are also helped by viscous spreading of the cores." · too tight/too loose: "Kicks grow: the vortices clump and the street tears; only one spacing survives." · non-staggered: "Facing rows are unstable at every spacing." · low k: "Long-wavelength wiggles grow more slowly than the alternating mode."
- **Derivation tab:** **D14** (15 steps) `view: 'street'`; goal `set` {b: 0.2805, off: 0.5}; step 4 (street speed) live "U_s = (Γ/2a)tanh(π b/a) = 0.3536 Γ/a"; step 6 (alternate disturbance) `watch` "adjacent vortices move in opposite directions"; step 9 live "sin(πζ) = −cosh(πb/a) = −1.4142"; step 11 (matrix) live γ, σ; step 12 (factorisation) `watch` "the quartic splits into two conjugate pairs"; step 14 `set` {b: 0.2805} `highlight: ['dot:marginal']`; step 15 `set` {off: 0} `watch` "growth π/4 at every b/a". Interpret: `s => "At b/a = " + b + " the largest growth rate is " + growth + " Γ/a²; the marginal spacing is " + ratio + "."`.
- **Code:**
  ```python
  b, off = {{b}}, {{off}}                          # b/a and stagger (0.5 = staggered)
  print(ch09.karman_street_ratio())                # arccosh(sqrt 2)/pi = 0.2805
  ev = ch09.karman_street_spectrum(b, offset=off)  # four eigenvalues (k = pi/a), units Gamma/a^2
  print(ev)                                        # {{ev}}
  g = ch09.karman_street_growth(b, offset=off)     # largest real part over all wavenumbers = {{growth}}
  print(g, ch09.karman_street_growth_closed(b))    # closed form (pi/2)|1/2 - sech^2(pi b)| = {{growthc}}
  print(ch09.karman_street_velocity(1.0, b, 1.0))  # street speed (Gamma/2a) tanh(pi b/a) = {{Us}}
  ```
- **Walkthrough (7 steps):** 1. "Two rows of vortices" — "Behind a cylinder alternating vortices form a staggered double row. Is its shape arbitrary? Press ▶ with b/a = 0.28." `set` stable preset, `play: true` · 2. "A row moves a vortex" — "A row of vortices has the velocity (Γ/2ia)cot(π(z − z₀)/a): the street travels at (Γ/2a)tanh(πb/a)." `derive: {id: 'D14', step: 4}` · 3. "Kick it" — "Press 'kick' at b/a = 0.15: the wiggle grows and the street tears." `set` too tight, `controls: ['kick', 'b']` · 4. "Linearise" — "For small displacements the motion is linear: a 4 × 4 matrix with γ = ½ − sech²(πb/a) and σ." `derive: {id: 'D14', step: 11}`, `readouts: ['gamma', 'sigma']` · 5. "Eigenvalues" — "Growth = (π/2)|γ|: zero exactly when cosh(πb/a) = √2." `eq: 'marg'`, `derive: {id: 'D14', step: 14}` · 6. "Facing rows" — "Set the stagger to 0: growth π/4 at every b/a." `set` non-staggered · 7. "Your turn" — "Find the two b/a values where the growth is 0.2." `controls: ['b']`.
- **Equations:** `row` "Velocity of a row of vortices (ours)" `w=u-iv=\frac{\Gamma}{2ia}\cot\frac{\pi(z-z_0)}a` live · `speed` "Street speed" `U_s=\frac\Gamma{2a}\tanh\frac{\pi b}a` live "U_s = …" · `growth` "Growth rate (alternating mode)" `\sigma=\frac{\pi\Gamma}{2a^2}\Big|\frac12-\mathrm{sech}^2\frac{\pi b}a\Big|` live · `marg` "Marginal spacing" `\cosh\frac{\pi b}a=\sqrt2\ \Leftrightarrow\ \frac ba=\frac1\pi\cosh^{-1}\sqrt2=0.2805` · `st` "Strouhal number" ref 'Eq. (4.102)' `St=\frac{\Omega d}{U_\infty}\approx0.2`.
- **Check yourself:** (1) "At b/a = 0.2 is the street stable? What is the growth?" — "No; σ = 0.298 Γ/a²." `set {b: 0.2}` · (2) "What are the two b/a values with σ = 0.1099?" — "0.25 and 0.313; the curve is V-shaped around 0.2805." · (3) "What happens to the growth when the rows face each other?" — "π/4 at every b/a: unstable." `set {off: 0}` · (4) "Why does the street move slower than the wind?" — "It moves at (Γ/2a)tanh(πb/a) relative to the fluid at rest at infinity; the cylinder carries the stream."
- **Selftest parity rows:** `{name: 'ratio', js: ratio(), py: 'ch09.karman_street_ratio()', rtol: 1e-14}` · `{name: 'growth closed 0.2', js: growthClosed(0.2), py: 'ch09.karman_street_growth_closed(0.2)', rtol: 1e-12}` · `{name: 'growth num 0.2', js: growthMax(0.2, 0.5), py: 'ch09.karman_street_growth(0.2)', rtol: 1e-5}` · `{name: 'growth nonstag', js: growthMax(0.28, 0), py: 'ch09.karman_street_growth(0.28, offset=0.0)', rtol: 1e-5}` · `{name: 'spectrum re', js: spectrum(0.5, 0.5, Math.PI)[0].re, py: 'ch09.karman_street_spectrum(0.5).real[0]', rtol: 1e-8}` (eigenvalues sorted by real part, then imaginary part in both) · `{name: 'speed', js: speed(1, 0.2805, 1), py: 'ch09.karman_street_velocity(1.0, 0.2805, 1.0)', rtol: 1e-12}` · invariant `{name: 'marginal', js: growthClosed(0.28054992616959), expect: 0, atol: 1e-9}`.
- **Fit plan:** 360×640: status, `street` (55 %) above `growth` (45 %), `spec` hidden (the four λ in the `growth` title: "λ = ±0.48 ± 0.62i"); 844×390: `street` | `growth`; desktop: rows [1.1, 1] with `spec` under `growth`.

### E7 · free_jet_similarity
- **Title:** "Where does a jet's extra mass come from?" · **Summary:** "Constant momentum flux J forces u₀ ∝ x^(−1/3) and δ ∝ x^(2/3); the sech² profile keeps ∫u²dy fixed while ṁ ∝ x^(1/3) grows as the jet entrains ambient fluid." · **CORE:** C12 (also N87–N113) · **Reference:** `amplitude_phase_second_order_II_3.html` (numbered derivation with live numbers; crosshair readouts).
- **meta:** `viz:order 7` · `viz:sections 9.10` · `viz:equations 9.57 9.59 9.62 9.63 9.64 9.70 9.71 9.72 9.73 9.75 9.76` · `viz:fluidpy ch09.free_jet ch09.free_jet_centreline ch09.free_jet_thickness ch09.free_jet_mass_flux ch09.free_jet_halfwidth ch09.free_jet_at_level ch09.free_jet_reynolds ch09.free_jet_constants ch09.jet_momentum_flux` · `viz:derivations D15 D16 D17 D18`.
- **Physics (JS ↔ Python):** `C = 4√6/3` (computed by Simpson of sech⁴ at start, parity with the closed form); `u0(x)`, `delta(x)`, `mdot(x)` ↔ `ch09.free_jet_centreline/thickness/mass_flux`; `jet(x, y)` → {u, v, eta} with f = √6 tanh(η/√6), f′ = sech² ↔ `ch09.free_jet`; `atLevel(level)` → {z = arccosh(1/√level), coeff = √6 z} ↔ `ch09.free_jet_at_level`; `h(x, level)` ↔ `ch09.free_jet_halfwidth`; `Rex(x)`, `vEdge` ↔ `ch09.free_jet_reynolds`, `free_jet_entrainment_velocity`; `momentum(y, u)` Simpson of ρ∫u²dy ↔ `ch09.jet_momentum_flux`.
- **Views** (rows [1.1, 1]): 1. `jet` "The jet" (row 0, flex 1.3) — x ∈ [0, 0.6] m, y ∈ ±h_max (mm, auto), colour u/u₀ (light to dark blue), streamlines from `psi` bending toward the axis (teal arrows at the edges show the entrainment velocity), the h-level envelope (muted dashed) and the station marker; title "u₀ 35.1 m/s · δ 0.207 mm · ṁ 0.0427 kg/(m s)". 2. `prof` "Profiles" (row 0, flex 1) — raw u(y) at x/4, x/2, x (blue, darker with x) or, in *rescaled* mode, u/u₀ against η with the exact sech²(η/√6) ghost (muted dashed); half-width marks at the chosen level (amber); click → probe. 3. `slopes` "Centre-line speed, half-width, mass flux against x" (row 1, `hidePortrait: true`) — log–log lines u₀ (blue, slope −1/3), h (amber, 2/3), ṁ (teal, 1/3) with fitted slopes printed; the *wrong exponent* toggle overlays u₀ ∝ x^{−1/2} (rose dotted) and shows that ρ∫u²dy then falls like x^{−1/4}.
- **Controls (≤ 5 visible):** `logJ` "Momentum flux $J$ (log₁₀ N/m)" −2 … 2, step 0.05, default 0 (J = 1 N/m) · `fluid` chips air / water · `x` "Station $x$" 0.02 … 1 m, step 0.01, default 0.1 (transport parameter) · `level` "Half-width level" 1 % … 50 %, step 1, default 1 · `wrong` toggle "test a wrong exponent" (optional).
- **Transport:** `x` 0.02 → 1 m, rate 0.1 m/s, `end: 'loop'`.
- **Presets:** "air slot jet (J = 1 N/m)" {fluid: air, logJ: 0, x: 0.1} · "water jet" {fluid: water, logJ: 0, x: 0.1} · "half-width at 1 % (corrected 7.3319)" {level: 1} · "half-width at 4 % (the book's printed 5.6152)" {level: 4}.
- **Status (rule-of-thumb bands, ours, flagged):** "🟢 Re_x = 40 < 100: a laminar jet is observable" · "🟡 100 ≤ Re_x < 1000: marginal" · "⚠️ Re_x = 2.3×10⁵ ≫ 1: the sech² profile has inflection points and is unstable (Ch. 11); real jets turn turbulent and spread ∝ x (Ch. 12)"; always followed by "momentum check ρ∫u²dy = 1.0000 N/m at x, x/2, x/4".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** · **status** · **inspector** (click a point in `jet`: "η = y/δ = 2.00, f′ = sech²(2/√6) = 0.5967, u = u₀f′ = 20.97 m/s; v = −(f − 2ηf′)/(3√Re_x)·u₀ = …") · **notes**.
- **Explain:**
  0. *What the views show* — "**The jet**: blue shade = speed relative to the centre line; teal arrows = ambient fluid being drawn in; the dashed lines mark where the speed has fallen to your level. **Profiles**: u(y) at three stations (or rescaled, on the exact sech² curve). **Slopes** (hidden on phones): the three power laws on log axes."
  1. *The invariant and its constant* — "J = ρ∫u²dy = **J** N/m at every x (9.57) $\frac d{dx}\int u^2dy=0$; with u = u₀f′(η), J/ρ = u₀²δC (9.61) and C = ∫f′²dη = ∫sech⁴(η/√6)dη = 4√6/3 = **C** (9.72)" boxed.
  2. *Speed and width* — "u₀ = [J²/(C²ρ²νx)]^{1/3} = **u0** m/s (9.62), δ = [Cρν²x²/J]^{1/3} = **delta** mm (9.63); doubling x: u₀ × **0.794**, δ × **1.587**" boxed.
  3. *The profile* — "u = u₀ sech²(η/√6), η = y/δ (9.71): at your probe η = **eta**, u = **u** m/s" boxed.
  4. *The mass that comes from outside* — "ṁ = ρu₀δ·2√6 = (36Jρ²νx)^{1/3} (9.73) = **mdot** kg/(m s); at x/2 it was **mdot2**: entrained between = **ent** kg/(m s) (× **ratio**)" boxed; hint "36^{1/3} = 3.302 is Bickley's coefficient."
  5. *Cross-flow* — "v/u₀ → ∓√6/(3√Re_x) (9.75) with Re_x = xu₀/ν = **Rex**: v = **vedge** m/s toward the jet" boxed.
  6. *Half-width at your level* — "sech²(z) = **level** ⇒ z = arccosh(1/√**level**) = **z**; h = √6 z δ = **coeff** δ = **h** mm (9.76). For 1 %: 7.3319 δ. The book's 5.6152 is arccosh 5 = 2.2924, the **4 %** point" boxed.
  7. *At the current station* — "x = **x** m: u₀ = **u0**, δ = **delta**, ṁ = **mdot**, momentum check ρ∫u²dy = **Jnum** N/m (Simpson on the drawn profile) (live)."
  8. *Reading the current setting* — "**Momentum stays, mass grows**: the jet slows because the same momentum is shared with more fluid: u₀ × 0.79 and ṁ × 1.26 per doubling of x. In air at this J the Reynolds number is huge, so a real jet would already be turbulent — the laminar solution shows two lessons (constant J, growing ṁ) that survive in round and turbulent jets. **Wrong exponent** (if on): u₀ ∝ x^{−1/2} (with δ = √(νx/u₀)) would make ρ∫u²dy fall like x^{−1/4}: momentum would vanish with no wall to take it."
- **Derivation tab:** **D15** (7 steps) `view: 'prof'`; goal `set` {wrong: false}; step 4 `watch` "the two advective terms are one derivative each"; step 7 live "ρ∫u²dy = **Jnum** = J". **D16** (14 steps) `view: 'slopes'`; step 4 live "u₀ ∝ x^{−1/3}: u₀(0.1) = 35.14, u₀(0.2) = 27.89"; step 12 `watch` "every power of x is x^{−5/3}: they cancel"; step 14 `set` {level: 1} live "3f‴ + ff″ + f′² = 0". **D17** (11 steps) `view: 'prof'`; `set` {rescaled}; step 8 `highlight: ['curve:sech2']`; step 9 live "C = 3.2660"; step 11 `view` note: "entrainment arrows on `jet`". **D18** (5 steps) `view: 'prof'`; step 4 `set` {level: 1} live "h = 7.3319 δ = **h** mm"; step 5 `set` {level: 4} live "h = 5.6153 δ". Interpret: `s => "At x = " + x + " m the jet is " + h + " mm half-wide at the " + level + " level and carries " + mdot + " kg/(m s), " + ratio + " × the slot's."`.
- **Code:**
  ```python
  J, rho, nu, x = {{J}}, {{rho}}, {{nu}}, {{x}}     # N/m, kg/m^3, m^2/s, m
  u0 = ch09.free_jet_centreline(x, J, rho, nu)      # (9.62): centre-line speed = {{u0}} m/s
  d  = ch09.free_jet_thickness(x, J, rho, nu)       # (9.63): delta = {{delta}} m
  md = ch09.free_jet_mass_flux(x, J, rho, nu)       # (9.73): (36 J rho^2 nu x)^(1/3) = {{mdot}} kg/(m s)
  h  = ch09.free_jet_halfwidth(x, J, rho, nu, level={{level}})   # corrected (9.76): {{h}} m
  print(ch09.free_jet_constants()["C"])             # 4 sqrt(6)/3 = 3.2660
  print(ch09.free_jet_reynolds(x, J, rho, nu))      # Re_x = {{Rex}}: how unstable the sech^2 profile is
  ```
- **Walkthrough (7 steps):** 1. "A jet in still air" — "Fluid leaves a slot. Downstream it is wider, slower — and carries more mass. Press ▶." `set` air preset, `play: true` · 2. "Momentum is conserved" — "No force acts on the jet (dp/dx = 0): d/dx ∫u²dy = 0 (9.57). The mass flux is not conserved: it grows." `derive: {id: 'D15', step: 7}`, `readouts: ['J']` · 3. "The exponents are forced" — "With u = u₀f′(η): J/ρ = u₀²δC. Constant J forces u₀ ∝ x^{−1/3}, δ ∝ x^{2/3} (9.62)–(9.63). Toggle 'wrong exponent'." `set` {wrong: true}, `derive: {id: 'D16', step: 4}` · 4. "sech²" — "The reduced ODE 3f‴ + ff″ + f′² = 0 integrates twice to f = √6 tanh(η/√6) (9.70), so u = u₀ sech²(η/√6) (9.71)." `set` {wrong: false}, `derive: {id: 'D17', step: 8}`, `eq: 'sech'` · 5. "Entrainment" — "The arrows show fluid drawn in from both sides: ṁ = (36Jρ²νx)^{1/3} (9.73)." `readouts: ['mdot']`, `derive: {id: 'D17', step: 10}` · 6. "Half-width" — "Choose 1 %: h = 7.3319 δ. The book's 5.6152 belongs to 4 %." `set` preset 1 %, `controls: ['level']`, `derive: {id: 'D18', step: 4}` · 7. "Your turn" — "Double J: predict u₀, δ, ṁ, then check." `controls: ['logJ']`.
- **Equations:** `mom` "Momentum flux" ref 'Eq. (9.57)' `\frac d{dx}\int_{-\infty}^{\infty}u^2\,dy=0` and 'Eq. (9.58)' `\int u^2dy=J/\rho` · `ans` "Similarity ansatz" ref 'Eq. (9.59)' `\psi=u_0\delta f(\eta),\ \delta=\sqrt{\nu x/u_0}` · `sech` "Bickley's jet" ref 'Eq. (9.71)' `u=u_0\,\mathrm{sech}^2(\eta/\sqrt6),\ \ f=\sqrt6\tanh(\eta/\sqrt6)` with 'Eq. (9.70)' in the note · `C` "The constant" ref 'Eq. (9.72)' `C=\int f'^2d\eta=\frac{4\sqrt6}{3}` and 'Eq. (9.73)' `\dot m=(36J\rho^2\nu x)^{1/3}` · `h99` "Half-width (corrected)" ref 'Eq. (9.76)' `h_{99}=7.3319\Big[\frac{C\rho\nu^2x^2}J\Big]^{1/3}` (note: the book prints 5.6152 = the 4 % point) · symbols.
- **Check yourself:** (1) "Double x: by what factors do u₀, δ, ṁ change?" — "2^{−1/3} = 0.794; 2^{2/3} = 1.587; 2^{1/3} = 1.260." · (2) "Where does the mass carried by the jet come from?" — "Ambient fluid drawn in from the sides by viscous friction (entrainment)." · (3) "Set the level to 4 %: which coefficient do you get, and whose is it?" — "5.6153: the book's printed 5.6152 (it is the 4 % point)." `set {level: 4}` · (4) "Why is a laminar jet at Re_x = 10⁵ not observed?" — "The sech² profile has inflection points and is unstable (Ch. 11)."
- **Selftest parity rows:** `{name: 'u0', js: u0(0.1, 1, 1.2, 1.5e-5), py: 'ch09.free_jet_centreline(0.1, 1.0, 1.2, 1.5e-5)', rtol: 1e-12}` · `{name: 'delta', js: delta(...), py: 'ch09.free_jet_thickness(0.1, 1.0, 1.2, 1.5e-5)', rtol: 1e-12}` · `{name: 'mdot', js: mdot(...), py: 'ch09.free_jet_mass_flux(0.1, 1.0, 1.2, 1.5e-5)', rtol: 1e-12}` · `{name: 'h99', js: h(..., 0.01), py: 'ch09.free_jet_halfwidth(0.1, 1.0, 1.2, 1.5e-5, level=0.01)', rtol: 1e-12}` · `{name: 'h printed', js: 5.61528786 * delta(...), py: 'ch09.free_jet_halfwidth(0.1, 1.0, 1.2, 1.5e-5, printed=True)', rtol: 1e-8}` · `{name: 'coeff 4%', js: atLevel(0.04).coeff, py: 'ch09.free_jet_at_level(0.04)["coeff"]', rtol: 1e-12}` · `{name: 'C', js: Csimpson(), py: 'ch09.free_jet_constants()["C"]', rtol: 1e-9}` · `{name: 'Re_x', js: Rex(...), py: 'ch09.free_jet_reynolds(0.1, 1.0, 1.2, 1.5e-5)["Re_x"]', rtol: 1e-12}` · invariant `{name: 'J numeric', js: momentum(profile(0.2)), expect: 1.0, rtol: 1e-6}`.
- **Fit plan:** 360×640: status, `jet` (45 %) above `prof` (55 %), `slopes` hidden (the three exponents in the `jet` title: "u₀ ∝ x^(−1/3) · h ∝ x^(2/3) · ṁ ∝ x^(1/3)"); 844×390: `jet` | `prof`; desktop rows [1.1, 1].

### E8 · wall_jet_invariant
- **Title:** "In a wall jet the wall eats momentum — what is conserved?" · **Summary:** "A wall removes ordinary momentum flux, but the flux of exterior momentum flux ∫u(∫_y^∞u²)dy survives: it fixes x u₀² = const, δ ∝ x^(3/4), and a free scale f_∞ that does not change the flow." · **CORE:** C13 (also N114–N124; the free jet of C12 as contrast) · **Reference:** `angular_frequency_explorer_1.html` (modes) and `amplitude_phase_second_order_II_3.html` (derivation with live numbers).
- **meta:** `viz:order 8` · `viz:sections 9.10` · `viz:equations 9.80 9.82 9.83 9.84 9.85` · `viz:fluidpy ch09.wall_jet_profile ch09.wall_jet ch09.wall_jet_ode_solve ch09.wall_jet_invariant ch09.wall_jet_integrals ch09.wall_jet_constants ch09.wall_jet_mass_flux ch09.jet_momentum_flux` · `viz:derivations D19 D20 D21`.
- **Physics (JS ↔ Python):** `odeSolve(fpp0, coeff)` RK4 (h = 0.01, η ≤ 40) of coeff·f‴ + ff″ + 2f′² = 0 ↔ `ch09.wall_jet_ode_solve` (coeff = 4 correct, 1 printed); `profile(η, finf)` = inverse of (9.83) by bisection on g ∈ [0, 1) ↔ `ch09.wall_jet_profile`; `integrals(finf)` = {∫f′, ∫f′², ∫f′∫f′²} by cumulative trapezoid of the IVP solution ↔ `ch09.wall_jet_integrals`; `wj(x, y, C, finf)` ↔ `ch09.wall_jet`; `invariant(y, u)` ↔ `ch09.wall_jet_invariant`; `constants(mdot, x, finf)` ↔ `ch09.wall_jet_constants`; `mdot(x, C, finf)` ↔ `ch09.wall_jet_mass_flux`; free-jet mode reuses E7's functions.
- **Views** (rows [1.1, 1]): 1. `jetxy` "The jet in x–y" (row 0, flex 1.3) — *wall* mode: the wall (hatched) at y = 0, the inner layer and the outer free-jet-like region, colour u/u_max, envelope δ(x) ∝ x^{3/4}, streamlines; *free* mode: the symmetric jet of C12 for comparison. 2. `prof` "f′(η) and f(η)" (row 0, flex 1) — f′ (rose) and f (blue) against η ∈ [0, 30/f_∞]: the implicit solution (9.83) (solid), the IVP integration from f″(0) = f_∞³/72 (teal dots), the peak marked, the far-field tail (purple dashed, with the 4.29 factor), f_∞ as a dotted level; in *free* mode sech² instead. 3. `bars` "Two integrals against x" (row 1, `hidePortrait: true`) — at four stations x = 0.1, 0.2, 0.4, 0.8 m grouped bars: ∫u²dy (rose), the double integral (teal), ṁ (blue, in kg/(m s)); the numbers printed above; *free* mode: ∫u²dy constant, the double integral not defined (greyed).
- **Controls (≤ 5 visible):** `mode` chips wall / free · `finf` "Gauge $f_\infty$" 0.5 … 3, step 0.05, default 1 (wall mode) · `x` "Station $x$" 0.05 … 1 m, default 0.3 (transport) · `mdot` "Mass flux ṁ at 1 m" 0.01 … 0.2 kg/(m s), default 0.05 · `printed` toggle "ODE as printed (coefficient 1)" (optional).
- **Transport:** `x` 0.05 → 1 m, rate 0.1 m/s, `end: 'loop'`.
- **Presets:** "gauge f_∞ = 1" {mode: wall, finf: 1} · "gauge f_∞ = 2 (C quarters)" {finf: 2} · "printed ODE (coefficient 1)" {printed: true} · "free jet for comparison" {mode: free}.
- **Terms:** two bars for the wall jet at the current x: ordinary momentum flux ρ∫u²dy (rose; falls like x^{−1/4}) and the invariant Ψ = ν x u₀² ∫f′∫f′² = C²ν f_∞⁴/40 (teal; constant), with their ratios to the values at x = 0.1 m; click a bar to isolate its curve in `bars`.
- **Status:** "✅ invariant constant: Ψ varies by 2e-9 over x = 0.1 … 0.8 m; ρ∫u²dy falls 19 %" · "🟠 free jet: ρ∫u²dy = J at every x (no wall)" · "⚠️ printed ODE: f_∞ = 0.397 instead of 1 and the (9.83) relation misses by 0.6–3" · "gauge f_∞ = 2: C = 28.9 (a quarter), dimensional profile unchanged (max difference 3e-12)".
- **Depth features:** Explain + Code + Derivation · **modes** (wall / free) · **terms** · **linked views** (3) · **presets** · **status** · **transport** · **notes**.
- **Explain:**
  0. *What the views show* — "**x–y**: a wall jet (slot at the left, wall at the bottom) with its thin inner layer and outer jet. **f′ and f**: the similarity profile; solid = the closed implicit solution (9.83), teal dots = an independent integration of 4f‴ + ff″ + 2f′² = 0. **Bars** (hidden on phones): the ordinary momentum flux (rose) and the exterior-momentum invariant (teal)."
  1. *The gauge* — "f′(0) = 0 here (unlike the free jet's f′(0) = 1), so the ODE has a free scale: f → λf(λη) sends C → C/λ² and f_∞ → λf_∞. Only K = C f_∞² is physical: K = (ṁ/ρ)²/(νx^{1/2}) = **K** m^{3/2}/s; with your f_∞ = **finf**: C = **C** m^{3/2}/s" boxed.
  2. *Centre speed and width* — "u₀ = Cx^{−1/2} = **u0** m/s at x = **x** (9.82); δ = [νx^{3/2}/C]^{1/2} = **delta** mm ∝ x^{3/4}; the peak speed u₀ f′_max = **umax** m/s at y = **ypk** mm" boxed.
  3. *Mass flux* — "ṁ = ρ√(νC) f_∞ x^{1/4} (9.84) = **mdot** kg/(m s): ∝ x^{1/4}, slower than the free jet's x^{1/3}" boxed.
  4. *The invariant* — "Ψ = C²ν ∫f′(∫_η^∞f′²)dη dη = C²ν f_∞⁴/40 (9.85) = **Psi** m⁴/s³ and Ψ = ṁ⁴/(40ρ⁴νx) — the same datum as ṁ (Ψ has units of u³L², not a force per length)" boxed.
  5. *What the wall takes* — "ρ∫u²dy = ρ u₀²δ ∫f′²dη = ρ C^{3/2}√ν f_∞³ x^{−1/4}/18 = **mom** N/m; at x = 0.1 m it was **mom1**: **frac** lost to the wall, while Ψ did not change" boxed.
  6. *Solving the ODE* — "f″(0) = f_∞³/72 = **fpp0** from f^{−1/2}f′ + f^{3/2}/6 = f_∞^{3/2}/6 at η → 0; the IVP lands on f_∞ = **finfnum**; the far field 1 − g ≈ 4.29 e^{−f_∞η/4}, g = √(f/f_∞) (the book omits 4.29)" boxed.
  7. *The printed ODE* — "With coefficient 1 instead of 4 the same f″(0) ends at f_∞ = **finfp** and misses (9.83) by **miss**: the 4 comes from the algebra of D20."
  8. *Reading the current setting* — wall: "The wall has removed momentum flux (rose bar falls) but the double integral is unchanged: 'the flux of exterior momentum flux' is conserved. It gives u₀ ∝ x^{−1/2}, δ ∝ x^{3/4} — thicker than the free jet's x^{2/3} because the wall slows the fluid." · free: "No wall: ρ∫u²dy is conserved and the exponents are −1/3, 2/3." · gauge: "Changing f_∞ changes C but not one dimensional profile — the gauge freedom of a scale-invariant ODE." · printed: "The printed equation is not the one that reduces (9.65): the residual shows it."
- **Derivation tab:** **D19** (13 steps) `view: 'bars'`; goal `set` {mode: 'wall'}; step 5 `watch` "the boundary terms vanish: no slip at the wall and u → 0 far away"; step 9 (9.79) `watch` "two integrals equal"; step 12 live "d/dx∫u(∫_y^∞u²) = 0: Ψ = **Psi**". **D20** (12 steps) `view: 'prof'`; step 9 `watch` "the ηf′f″ terms cancel again"; step 11 `set` {printed: true} `watch` "the printed coefficient leaves 3f‴ over"; step 12 `set` {printed: false}. **D21** (15 steps) `view: 'prof'`; step 4 live "4ff″ − 2f′² + f²f′ = 0 at η = 0 ✓"; step 8 `set` {printed: false}; step 12 (9.83) live "residual of (9.83) at η = 5: 3e-13"; step 14 live "f″(0) = f_∞³/72 = **fpp0**"; step 15 `set` {finf: 2} `watch` "C quarters, profile identical". Interpret: `s => "With f_∞ = " + finf + " and ṁ = " + mdot + " kg/(m s) at 1 m: C = " + C + ", u₀(1 m) = " + u0 + " m/s, peak speed " + umax + " m/s."`.
- **Code:**
  ```python
  finf, mdot1, x = {{finf}}, {{mdot}}, {{x}}         # gauge, mass flux at 1 m [kg/(m s)], station [m]
  sol = ch09.wall_jet_ode_solve(fpp0=finf**3/72)     # 4f''' + f f'' + 2 f'^2 = 0 from f''(0) = f_inf^3/72
  print(sol["f_inf"], sol["err_vs_9_83"])            # lands on f_inf = {{finfnum}} without being told
  print(ch09.wall_jet_integrals(finf))               # {{ints}}: int f', int f'^2 = f^3/18, invariant f^4/40
  c = ch09.wall_jet_constants(1.2, 1.5e-5, mdot=mdot1, x=1.0, f_inf=finf)   # one datum: Psi and mdot are the same
  print(c["C"], c["Psi"])                            # C = {{C}}, Psi = {{Psi}}
  wj = ch09.wall_jet(x, y, c["C"], finf, 1.5e-5, rho=1.2)   # dimensional u(y) at the station
  print(ch09.wall_jet_invariant(y, wj["u"]))         # constant in x: {{inv}}
  ```
- **Walkthrough (7 steps):** 1. "A jet along a wall" — "Fluid leaves a slot along a wall. The wall pulls back, so ordinary momentum flux must fall. Press ▶ on 'free' vs 'wall'." `set` wall preset, `play: true` · 2. "What falls and what stays" — "In the bars, ρ∫u²dy (rose) falls with x; the double integral ∫u(∫_y^∞u²)dy (teal) does not (9.80)." `terms: true`, `derive: {id: 'D19', step: 12}` · 3. "The exponents" — "Insert u = u₀f′: x u₀² = const ⇒ u₀ ∝ x^{−1/2}, δ ∝ x^{3/4} (9.82)." `derive: {id: 'D19', step: 13}`, `eq: 'exp'` · 4. "The ODE" — "The reduced equation is 4f‴ + ff″ + 2f′² = 0. The printed coefficient 1 is wrong: press 'printed ODE'." `set` printed preset, `derive: {id: 'D20', step: 11}` · 5. "Two integrations" — "An integrating factor f and an exact derivative give f^{−1/2}f′ + f^{3/2}/6 = f_∞^{3/2}/6, then partial fractions give (9.83)." `derive: {id: 'D21', step: 12}` · 6. "A free scale" — "Set f_∞ = 2: C quarters but the dimensional profile does not move." `set` f_∞ = 2, `controls: ['finf']` · 7. "Your turn" — "Predict ṁ and the peak speed at 4 m, then check." `controls: ['x', 'mdot']`.
- **Equations:** `inv` "Flux of exterior momentum flux" ref 'Eq. (9.80)' `\frac d{dx}\int_0^\infty u\Big(\int_y^\infty u^2dy'\Big)dy=0` · `exp` "Exponents" ref 'Eq. (9.82)' `u_0=Cx^{-1/2},\ \delta=[\nu x^{3/2}/C]^{1/2}` · `ode` "Reduced equation" `4f'''+ff''+2f'^2=0` (the book prints coefficient 1) · `imp` "Implicit solution" ref 'Eq. (9.83)' `-\ln(1-g)+\sqrt3\tan^{-1}\frac{2g+1}{\sqrt3}+\ln(1+g+g^2)^{1/2}=\frac{f_\infty}4\eta+\sqrt3\tan^{-1}\frac1{\sqrt3}` (g² = f/f_∞) · `mdot` "Entrainment" ref 'Eq. (9.84)' `\dot m=\rho\sqrt{\nu C}\,f_\infty x^{1/4}` with 'Eq. (9.85)' `\Psi` in the note.
- **Check yourself:** (1) "Which bar stays flat in wall mode?" — "The teal double integral." · (2) "Double f_∞: what happens to C, to u₀, to δ and to the dimensional profile?" — "C → C/4, u₀ → u₀/4, δ → 2δ, while f′ peaks four times higher (f′ ∝ f_∞²): u₀f′ and hence the dimensional profile are unchanged." `set {finf: 2}` · (3) "Why does the wall jet spread as x^{3/4} instead of x^{2/3}?" — "The wall slows the fluid, so the layer must be thicker for the same invariant." · (4) "What is f″(0) for f_∞ = 2?" — "8/72 = 0.111 (f_∞³/72)."
- **Selftest parity rows:** `{name: 'integrals', js: integrals(1).inv, py: 'ch09.wall_jet_integrals(1.0)["invariant"]', rtol: 1e-8}` · `{name: 'int f2', js: integrals(1).int_fp2, py: 'ch09.wall_jet_integrals(1.0)["int_fp2"]', rtol: 1e-8}` · `{name: 'f_inf', js: odeSolve(1/72, 4).finf, py: 'ch09.wall_jet_ode_solve(fpp0=0.013888888888)["f_inf"]', rtol: 1e-7}` · `{name: 'fp at 8', js: profile(8, 1).fp, py: 'ch09.wall_jet_profile(8.0)["fp"]', rtol: 1e-8}` · `{name: 'printed f_inf', js: odeSolve(1/72, 1).finf, py: 'ch09.wall_jet_ode_solve(fpp0=0.013888888888, coeff=1.0)["f_inf"]', rtol: 1e-7}` · `{name: 'mdot', js: mdotFn(1, 115.7, 1), py: 'ch09.wall_jet_mass_flux(1.0, 115.7, 1.0, 1.2, 1.5e-5)', rtol: 1e-12}` · `{name: 'constants C', js: constants(0.05, 1, 1).C, py: 'ch09.wall_jet_constants(1.2, 1.5e-5, mdot=0.05, x=1.0, f_inf=1.0)["C"]', rtol: 1e-10}` · invariant `{name: 'Psi vs x', js: PsiSpread(), expect: 0, atol: 1e-7}`.
- **Fit plan:** 360×640: status, `prof` (55 %) above `bars` (45 %) — `bars` carries the point of the explainer, so `jetxy` is the view hidden on portrait (`hidePortrait`) and the invariant values are repeated in the `bars` title; 844×390: `prof` | `bars`; desktop rows [1.1, 1] with `jetxy` above-left.

### E9 · teacup_secondary_flow
- **Title:** "Why do tea leaves collect at the centre?" · **Summary:** "A thin layer slows the swirl at the floor but not the pressure gradient set by the fast water above: the unbalanced pressure ρ(u_e² − u²)/R drives fluid inward — the seed of Ekman layers." · **CORE:** C14 (also the Ekman hook) · **Reference:** `angular_frequency_explorer_1.html` (modes) and `fid_formula_lab.html` (term bars that add up).
- **meta:** `viz:order 9` · `viz:sections 9.11` · `viz:equations` none numbered (`∂p/∂R = ρu²/R` is ch04's; the net-force formula is ours) · `viz:fluidpy ch09.secondary_flow_radial_force ch09.secondary_flow_layer_profile` · `viz:derivations D22`.
- **Physics (JS ↔ Python):** `layer(z, δ, ue, shape)` illustrative swirl profile (exponential u_e(1 − e^{−z/δ}), linear, 1/7-power) ↔ `ch09.secondary_flow_layer_profile`; `force(ue, u, R, ρ)` = ρ(u_e² − u²)/R ↔ `ch09.secondary_flow_radial_force`; `parts(z)` → {pressure gradient force ρu_e²/R, centrifugal ρu²/R, net}; `height10()` where net falls to 10 % of the floor value; `Omega = ue/R`, `ekman = √(ν/Ω)`, `spinup = H/√(νΩ)` (preview constants, closed forms).
- **Views** (rows [1.1, 1]): 1. `cup` "The cup" (row 0, flex 1.3) — cross-section: water depth H = 5 cm, radius R, swirl arrows (blue, longer where faster), the meridional circulation (teal loop: in along the floor, up the axis, out along the top, down the side wall) with speeds ∝ the net force, leaves (brown dots) drifting inward on the floor on the transport clock; *river bend* mode: plan view of a bend, the same balance across the channel (fast outside/slow near the bed). 2. `prof` "u(z) and the two forces" (row 0, flex 1) — u(z)/u_e against z (blue), the pressure-gradient force ρu_e²/R (orange, uniform in z) and the centrifugal force ρu²/R (rose) on a force axis [N/m³], the difference shaded teal (the net inward force); a draggable height marker (`z`). 3. `net` "Net inward force against z" (row 1, `hidePortrait: true`) — F(z) with the 10 % height and the floor value.
- **Controls (≤ 5 visible):** `ue` "Core swirl speed $u_e$" 0.05 … 0.5 m/s, step 0.01, default 0.2 · `R` "Radius $R$" 2 … 10 cm, default 4 · `delta` "Layer thickness $\delta$" 0.5 … 5 mm, default 1.5 · `shape` chips exponential / linear / power 1/7 (illustrative) · `z` "Height above the floor" 0 … 10 mm, default 0 (probe) · optional `mode` chips teacup / river bend.
- **Transport:** `t` 0 → 30 s, rate 3 s/s, `end: 'loop'` (leaves drift; the loop turns).
- **Presets:** "tea cup (u_e 0.2 m/s, R 4 cm)" · "thick layer (viscous oil, δ 4 mm)" {delta: 4} · "slow stirring" {ue: 0.05} · "river bend (u_e 1 m/s, R 50 m)" {mode: river bend, ue: 1, R: 50}.
- **Terms:** bars at height z: pressure-gradient force ρu_e²/R (orange) − centrifugal force ρu²/R (rose) = net inward force (teal); the bars add: orange = rose + teal (click a bar to isolate it in `prof`).
- **Status:** "⬅️ inflow: net inward force 1000 N/m³ at the floor (10 % of ρg)" · "⚖️ at z = 8 mm the layer has recovered: net force 9.7 N/m³ (< 1 % of the floor value)" · "✓ core: 0 (u = u_e: the pressure gradient exactly supplies the centripetal acceleration)".
- **Depth features:** Explain + Code + Derivation · **modes** (teacup / river bend) · **terms** · **linked views** (3) · **presets** · **status** · **transport** · **notes** (Ch. 13 hook).
- **Explain:**
  0. *What the views show* — "**The cup**: swirl arrows (blue), the meridional loop (teal) and leaves drifting inward. **u(z) and forces**: orange = the pressure-gradient force per volume that the fast core needs and sets; rose = the centrifugal force the slowed layer needs; the teal gap is the net inward push. **Net force** (hidden on phones): the same gap against height."
  1. *The pressure gradient set by the core* — "∂p/∂R = ρu_e²/R = 1000 × **ue**² / **R** = **dpdR** N/m³ (∂p/∂z ≈ 0 across the thin layer, (9.10) with z for y)" boxed.
  2. *What the slowed fluid needs* — "at z = **z** mm: u = **u** m/s, ρu²/R = **cf** N/m³" boxed.
  3. *The net force* — "F_in = ρ(u_e² − u²)/R = **dpdR** − **cf** = **F** N/m³ (inward = positive); on the floor F = **Ffloor** N/m³ = **pg** % of ρg" boxed.
  4. *How far it reaches* — "F falls to 10 % of its floor value at z = **z10** mm (illustrative profile, thickness δ = **delta** mm)".
  5. *Ekman preview (Ch. 13, not derived here)* — "Ω = u_e/R = **Om** s⁻¹, δ ~ √(ν/Ω) = **ek** mm, spin-up time H/√(νΩ) = **spin** s against the diffusion time H²/ν = **diff** s" boxed.
  6. *Values at the current time* — "leaves at r = **r** cm after t = **t** s (live)."
  7. *Reading the current setting* — teacup: "The water at the floor is slowed by friction but the pressure gradient is unchanged, so the net force is inward: the floor flow carries leaves to the centre; continuity sends the water up the axis and out at the top." · thick layer: "A thicker layer means more of the depth is slowed: a taller shaded wedge and a stronger, deeper inflow." · slow stirring: "F scales like u_e²: stirring half as fast gives a quarter of the push." · river bend: "The same balance across a channel: the slow bed layer is pushed toward the inside of the bend — sand collects on the inner bank (Exercise 9.28)." · core: "No layer, no imbalance: the pressure gradient exactly supplies the centripetal acceleration."
- **Derivation tab:** **D22** (7 steps) `view: 'prof'`; goal `set` preset tea cup; step 2 `watch` "the orange bar is the same at every height"; step 3 live "F = ρ(u_e² − u²)/R at z: **F** N/m³"; step 5 live "floor: 1000 N/m³ (10 % of ρg)"; step 6 `view` note: "the loop on `cup`"; step 7 live Ekman scales. Interpret: `s => "At z = " + z + " mm the net inward force is " + F + " N/m³; on the floor " + Ffloor + " N/m³."`.
- **Code:**
  ```python
  rho, ue, R, delta = 1000.0, {{ue}}, {{R}}, {{delta}}      # water, core swirl [m/s], radius [m], layer thickness [m]
  z = np.linspace(0, 5*delta, 101)                           # heights above the floor [m]
  u = ch09.secondary_flow_layer_profile(z, delta, ue, shape="{{shape}}")   # illustrative swirl profile
  F = ch09.secondary_flow_radial_force(ue, u, R, rho)        # rho (ue^2 - u^2)/R [N/m^3], inward > 0
  print(F[0], F[-1])                                         # floor {{Ffloor}}, top {{Ftop}} N/m^3
  print(rho*9.81)                                            # weight density 9810 N/m^3: the push is {{pg}} % of it
  Om = ue/R; print(np.sqrt(1.0e-6/Om))                       # Ekman-like layer depth sqrt(nu/Omega) = {{ek}} m
  ```
- **Walkthrough (6 steps):** 1. "Where do the leaves go?" — "Stir tea, let it spin: the leaves gather in the middle. Press ▶." `set` tea cup preset, `play: true` · 2. "The core sets the pressure" — "Fast water on circles needs ∂p/∂R = ρu_e²/R. The thin layer feels the same gradient (orange)." `derive: {id: 'D22', step: 2}` · 3. "A slow layer needs less" — "Friction slows the water near the floor: the centrifugal need ρu²/R (rose) falls." `terms: true` · 4. "The gap pushes inward" — "Net force ρ(u_e² − u²)/R on the floor: 1000 N/m³, a tenth of the weight." `readouts: ['F']`, `derive: {id: 'D22', step: 5}` · 5. "Thicker layer, stronger inflow" — "Try the oil preset." `set` thick layer, `controls: ['delta']` · 6. "Your turn" — "Predict the inflow force at half the stirring speed, then check; then try the river bend." `controls: ['ue']`.
- **Equations:** `rad` "Radial balance in the core (Ch. 4)" `\frac{\partial p}{\partial R}=\frac{\rho u_e^2}R` · `thin` "Thin layer" `\frac{\partial p}{\partial z}\approx0` (compare (9.10) $0=-\frac{\partial p}{\partial y}$) · `net` "Net inward force per volume (ours)" `F_{in}=\frac{\rho\,(u_e^2-u^2)}R` live "F = … N/m³" · `ek` "Ekman-layer scale (Ch. 13 preview)" `\delta_E\sim\sqrt{\nu/\Omega}`.
- **Check yourself:** (1) "What is the net force in the core (u = u_e)?" — "Zero: the pressure gradient exactly supplies the centripetal acceleration." · (2) "Halve u_e: by what factor does the floor force change?" — "1/4 (∝ u_e²)." `set {ue: 0.1}` · (3) "Why does the inflow happen only near the floor?" — "Only there is u < u_e (no slip)." · (4) "Which way does sand move in a river bend?" — "Toward the inside bank along the bed." `set` river-bend preset.
- **Selftest parity rows:** `{name: 'force mid', js: force(0.2, 0.1, 0.04, 1000), py: 'ch09.secondary_flow_radial_force(0.2, 0.1, 0.04, 1000.0)', rtol: 1e-12}` · `{name: 'force floor', js: force(0.2, 0, 0.04, 1000), py: 'ch09.secondary_flow_radial_force(0.2, 0.0, 0.04, 1000.0)', rtol: 1e-12}` · `{name: 'profile exp', js: layer(1.5e-3, 1.5e-3, 0.2, 'exponential'), py: 'ch09.secondary_flow_layer_profile(0.0015, 0.0015, 0.2, shape="exponential")', rtol: 1e-12}` · `{name: 'profile 1/7', js: layer(1e-3, 1.5e-3, 0.2, 'power'), py: 'ch09.secondary_flow_layer_profile(0.001, 0.0015, 0.2, shape="power")', rtol: 1e-12}` · invariant `{name: 'core net zero', js: force(0.2, 0.2, 0.04, 1000), expect: 0, atol: 1e-9}`.
- **Fit plan:** 360×640: status, `cup` (45 %) above `prof` (55 %), `net` hidden (the floor value in the `prof` title); 844×390: `cup` | `prof`; desktop rows [1.1, 1].

### B1 · ball_swing_magnus
- **Title:** "Why does the seam or the spin change which way the ball curves?" · **Summary:** "Which side of the ball has its layer past the crisis decides the side force — swing, negative and positive Magnus effect on one control set." · **CORE:** C11 (also N84, N85, N86) · **Reference:** `angular_frequency_explorer_1.html` (modes, presets). Built only if an explainer above fails review (the notebook does not embed it).
- **Physics:** `magnusSign(Re_slow, Re_fast, Re_cr)` ↔ `ch09.magnus_sign`; `swing(F/W, d, U)` ↔ `ch09.ball_swing_deflection`; `cpSide(phis, cb)` ↔ `ch09.separated_cp` per side; speed ↔ Re = U d/ν.
- **Views:** `ball` (ball in a stream with the two separation points and C_p bars, the force arrow) · `table` (truth table: side Re against Re_cr with the current cell lit) · `path` (trajectory y = ½ a t², hidePortrait). **Controls:** speed 10…50 m/s · spin (rev/s, sign) · surface roughness · seam side · Re_cr slider (optional). **Modes:** cricket / tennis. **Presets:** "too slow", "swing window", "both turbulent", "rough tennis ball". **Status:** "− (negative Magnus)" etc. **Depth features:** modes, presets, status, terms (C_p bars). **Selftest:** `ch09.magnus_sign(0.8e5, 1.2e5, 3e5)` = "+", `(2e5, 4e5, 3e5)` = "−", `ch09.ball_swing_deflection(0.2, 18.0, 35.0)` = 0.2595 m. No derivations. Walkthrough 5 steps, ≥ 3 check questions, fit plan as E5.

## Part D — runtime budget (full run < 5 min on a laptop / Colab CPU)

Chapter 9 is mostly closed forms (µs per point) plus five numerical workhorses (the Falkner–Skan `solve_bvp` with continuation ≈ 60 ms per exponent, the von Mises marching solver, the closure table of 60 exponents, the wall-jet IVP, the point-vortex lattice sums) and nine sympy engines. The costs: `falkner_skan_state` × 8 in C05 ≈ 0.5 s (each solution `lru_cache`d, key = m rounded to 1e-9, so C06–C08 reuse them), the fold search ≈ 1 s,
`thwaites_closure_table(60)` ≈ 3 s (built once per kernel and reused by C08 and F4), `march_boundary_layer` (400 ψ-nodes × 60 steps with a Picard–backward-Euler step) ≈ 2 s (its solution is computed once in C01 and reused by A1 and C08), `bl_pressure_variation` × 4 ≈ 1 s, sympy (`bl_nondim_sympy` ≈ 1 s, `similarity_reduce_sympy` × 4 cases ≈ 6 s, `momentum_integral_sympy`, `thwaites_sympy`, `cylinder_thwaites_sympy`, `jet_momentum_sympy`,
`wall_jet_invariant_sympy`, `wall_jet_sympy`, `karman_street_sympy` ≈ 1–3 s each — all `lru_cache`d), the Kármán scan (61 wavenumbers × 4 lattice sums × 4001 terms ≈ 0.1 s per b/a; 8 values in C10 ≈ 0.8 s, the 20-point growth curve ≈ 2 s), 6 animations and 7 plotly figures. ch08 ran in ≈ 196 s with 15 CORE blocks; ch09 has 14.

| Section | Heaviest cells | Full | FAST (`FLUIDPY_FAST=1`) |
|---|---|---|---|
| setup + imports | numpy/scipy/sympy/plotly, three new core modules | 9 s | 9 s |
| §9.1 C01 | `bl_nondim_sympy`, Blasius field 200 × 120 with `np.gradient` at three Re, marching 60 stations, `bl_pressure_variation` × 4, three-panel figure | 8 s | 5 s (ny 200, 30 stations) |
| §9.2 C02 | thicknesses on 4001 points, four shapes, trapezoid from scratch, three-panel figure, **A1 video 60**, F1 (25 steps) | 14 s | 8 s (30 frames, 12 steps) |
| §9.3 C03 | `similarity_reduce_sympy("blasius")`, two marching runs, two-panel figure | 5 s | 4 s |
| §9.3 C04 | `blasius_constants`, three truncations, RK4 by hand (1000 steps), three figures, **A2 video 60**, F2 (25 × 2 traces) | 16 s | 9 s |
| §9.4 C05 | 8 × `falkner_skan_state`, `falkner_skan_separation`, two `brentq` shoots, family + fold figure, **A3 frames 30**, F3 (30 steps), live (not run on the page) | 11 s | 7 s (16 frames, 15 steps) |
| §9.5 C06 | 3 FS states, residual tests, two Kármán–Pohlhausen runs, `momentum_integral_sympy`, budget figure | 7 s | 4 s |
| §9.6 C07 | `thwaites_closure_table(60)`, `thwaites_sympy`, `cylinder_thwaites_sympy`, diffuser + cylinder marches, accuracy table (5 FS states), three figures, **A4 frames 30**, F4 (30 steps) | 21 s | 12 s (30 rows, 16 frames) |
| §9.7 C08 | FS profiles, τ₀ of the diffuser, sign-change search, figure, F7 (25 steps) | 5 s | 3 s |
| §9.7 C09 | plate drag curves, separated-model drag, Gauss–Legendre, two figures, F5 (25 steps) | 6 s | 4 s |
| §9.8 C10 | Kármán growth table (7 values), Jacobian eigenvalues, `karman_street_sympy`, three-panel figure (20-point growth curve), **A5 video 60** | 13 s | 8 s (k-scan 21, n_terms 1000, 30 frames) |
| §9.8 C11 + §9.9 | drag-crisis states, schematic and Morrison curves, swing and Magnus figures | 7 s | 5 s |
| §9.10 C12 | `similarity_reduce_sympy("free_jet")`, `free_jet_ode_solve`, momentum flux at 4 stations (8001 points), 200 × 200 jet field, F6 (25 steps), **A6 video 60** (shared with C13) | 15 s | 9 s (100², 30 frames) |
| §9.10 C13 | `similarity_reduce_sympy("wall_jet")` × 2, `wall_jet_sympy`, `wall_jet_invariant_sympy`, three IVPs, three-panel figure | 9 s | 6 s |
| §9.11 C14 | closed forms, cup figure | 2 s | 2 s |
| explainers (9 × `show_viz`) | read the HTML files | 1 s | 1 s |
| **Total** | | **≈ 149 s** | **≈ 96 s** |

**FAST plan.** Every size-dependent choice is written `a if not FAST else b` in the cell: marching ny 400 → 200 and stations 60 → 30; closure table 60 → 30 exponents; Kármán scan 61 → 21 wavenumbers and n_terms 2000 → 1000; contour grids 200² → 100²; videos 60 → 30 frames, frame players 30 → 16; plotly sliders ≤ 25 → ≤ 12 steps (≤ 4 traces × ≤ 400 points, < 250 kB each); the wall-jet IVPs use `rtol=1e-10` instead of 1e-12. **Cached arrays / results:** `falkner_skan` solutions (`lru_cache`, key m rounded to 1e-9; C05–C08 share them), the Blasius profile (`blasius_constants` cached; C01–C04 and A1–A2 share it), `thwaites_closure_table` (once per kernel), every
`*_sympy` engine (`functools.lru_cache`), the marching solution of C01 (reused by A1 and by C08's diffuser comparison), and the wall-jet IVP (reused by C13's three figures). Sympy cells never call `simplify` on expressions with more than one generic function: `similarity_reduce_sympy` substitutes the ansatz for ψ, differentiates once and simplifies only the final bracketed residual. **Outputs:** 4 videos (A1, A2, A5, A6) + 2 frames players (A3, A4) at dpi 80 (each < 3 MB), 7 plotly figures, ≈ 36 static figures — the page stays under 15 MB.

---

## Part E — prerequisite ledger
Every concept, symbol, maths tool and Python function or idiom the notebook or its explainers use, with where it is explained. "primer (in Cxx)" = a 📎 primer placed in that block before first use (the Concept text is the primer term, used verbatim in `nb.primer`); "knowledge/primers.md: <term> (chNN Pnn) — reminder" = a one-line reminder naming the earlier primer; a CORE/RECAP id alone = taught there; "Cxx (Nnn)" = the NOTE placed in that block;
"Cxx (D0n)" = the derivation step where it is used; "Cxx (gloss …)" = one sentence where it is used. Earlier chapters' material that is neither a primer nor a ch09 RECAP is named by section ("Ch. 4 §4.9"). New primers P200–P220 in first-use order: six Reynolds numbers, parabolic/elliptic and marching (C01) · improper integral of a deficit, simpson and trapezoid, PCHIP, control volume with a streamline side (C02) · chain rule with η(x, y) (C03) ·
Töpfer scaling symmetry, shooting versus BVP (C04) · continuation and fold (C05) · integrating factor, ∫sin⁵ (C07) · inflection point (C08) · row of vortices, linear stability (C10) · sech and arccosh, total-derivative move, momentum and mass flux (C12) · integration by parts with a variable limit, partial fractions and implicit inverse (C13) · radial force balance in a swirl (C14) — 21 primers.

| Concept | First used in | Explained by |
|---|---|---|
| Prandtl's boundary-layer hypothesis, inner and outer problem | C01 | C01 (N01) |
| d'Alembert's paradox (zero drag in ideal flow) | C01 | Ch. 6 §6.5 (recalled in C01, N01) |
| irrotational outer flow and Bernoulli along a streamline | C01 | Ch. 4 §4.9 (recalled in C01, D02) |
| steady 2-D momentum along a surface (9.1) | C01 | R01 |
| continuity (6.2) and the stream function u = ∂ψ/∂y, v = −∂ψ/∂x | C01 | R03 |
| scaled continuity (8.15) | C01 | R02 |
| no-slip condition (9.12) | C01 | R04 |
| no through-flow condition (9.13) | C01 | R05 |
| six Reynolds numbers (which length?) | C01 | primer (in C01) |
| kinematic viscosity ν = μ/ρ and its diffusion time | C01 | knowledge/primers.md: diffusivity and the diffusion time L²/ν (ch08 P185) — reminder |
| order-of-magnitude scaling, the symbol ~ | C01 | knowledge/primers.md: order-of-magnitude scaling (ch04 P130) — reminder |
| limits and orders of smallness | C01 | knowledge/primers.md: limits and orders of smallness (ch02 P68) — reminder |
| anisotropic scaling with two length scales | C01 | knowledge/primers.md: anisotropic scaling with two length scales (ch08 P188) — reminder |
| scaled variables and the chain rule | C01 | knowledge/primers.md: scaled variables and the chain rule (ch04 P133) — reminder |
| chain rule for derivatives | C01 | knowledge/primers.md: chain rule (ch01 P49) — reminder |
| dominant balance | C01 | knowledge/primers.md: dominant balance (ch08 P198) — reminder |
| partial derivative | C01 | knowledge/primers.md: partial derivative (ch01 P25) — reminder |
| sympy symbols, subs, diff, simplify | C01 | knowledge/primers.md: sympy (ch01 P40) — reminder |
| boundary-layer thickness δ̄/L ~ Re^{−1/2} (9.4), the √(νt) bridge | C01 | C01 (N04, D01) |
| sizes of advective and viscous terms (9.2), (9.3), derivative sizes (9.5), size of v | C01 | C01 (N02, N03, N05, N06) |
| stretched variables (9.6) | C01 | C01 (N07, D01) |
| scaled x- and y-momentum equations (9.7), (9.8) | C01 | C01 (N08, N09, D01) |
| boundary-layer equations (9.9), (9.10) and their limit Re → ∞ | C01 | C01 (D01) |
| ∂p/∂y = 0: pressure imposed by the outer flow | C01 | C01 (N10) |
| matching to the outer flow (9.11), favourable and adverse gradient | C01 | C01 (N11, D02) |
| edge condition (9.14), inlet profile (9.15) | C01 | C01 (N12, N13) |
| parabolic, elliptic and marching in x | C01 | primer (in C01) |
| crude wall-stress estimate τ₀ ~ μU/δ̄, C_f ~ 2/√Re | C01 | C01 (N16) |
| viscous–inviscid iteration, where the approximation fails | C01 | C01 (N15, N17) |
| finite-difference residual with np.gradient | C01 | knowledge/primers.md: np.gradient (ch01 P22) — reminder |
| power laws and log–log plots | C01 | knowledge/primers.md: power laws and log–log plots (ch01 P13) — reminder |
| observed order of a scheme (tools.convergence.observed_order) | C01 | Ch. 8 §8.4 (used there for the Crank–Nicolson order) |
| BL.boundary_layer_scales, BL.bl_pressure_variation, BL.bl_x_momentum_residual, ch09.bl_nondim_sympy | C01 | C01 (code explain) |
| BL.march_boundary_layer, BL.outer_flow | C01 | C01 (N14, D02, code explain) |
| Python dictionaries (fluidpy results) | C01 | knowledge/primers.md: Python dictionaries (ch01 P23) — reminder |
| matplotlib figures | C01 | knowledge/primers.md: matplotlib figures (ch01 P01) — reminder |
| f-strings | C01 | knowledge/primers.md: f-strings (ch01 P04) — reminder |
| tuple unpacking | C01 | knowledge/primers.md: tuple unpacking (ch01 P14) — reminder |
| functions as arguments and lambda | C01 | knowledge/primers.md: functions as arguments and lambda (ch01 P29) — reminder |
| assert np.allclose | C01 | knowledge/primers.md: assert np.allclose (ch01 P15) — reminder |
| delta_99, the 99 % thickness | C02 | C02 (N18) |
| displacement thickness δ* (9.16) | C02 | C02 (D03) |
| momentum thickness θ (9.17), shape factor H = δ*/θ | C02 | C02 (N20, D04) |
| streamline lift, v(y) = −∫u_x dy′ | C02 | C02 (N19, D03) |
| plate drag from θ, ρU²θ = ∫τ₀dx | C02 | C02 (N21, D04) |
| definite integral | C02 | knowledge/primers.md: definite integral (ch01 P27) — reminder |
| improper integral of a deficit (truncating the tail) | C02 | primer (in C02) |
| scipy.integrate.simpson and np.trapezoid | C02 | primer (in C02) |
| trapezoid rule | C02 | knowledge/primers.md: trapezoid rule (ch01 P37) — reminder |
| monotone interpolation PchipInterpolator | C02 | primer (in C02) |
| root finding with brentq | C02 | knowledge/primers.md: scipy.optimize.brentq (ch03 P108) — reminder |
| control volume with a streamline as a side | C02 | primer (in C02) |
| control-volume momentum balance | C02 | Ch. 4 §4.4 (recalled in D04) |
| BL.thicknesses, BL.profile_shape, BL.blasius_fields | C02 | C02 (code explain) |
| animate and show_animation | C02 | knowledge/primers.md: animate and show_animation (ch01 P16) — reminder |
| slider_figure | C02 | knowledge/primers.md: slider_figure (ch01 P17) — reminder |
| show_viz (embedded explainers) | C02 | knowledge/primers.md: show_viz (ch01 P18) — reminder |
| flat-plate layer equation (9.18) | C03 | C03 (N22) |
| similarity variable η = y/δ(x) and similarity solutions | C03 | C03 (N23); Ch. 8 §8.4 (recalled) |
| dimensional argument for ψ = Uδf(η) (9.19) | C03 | C03 (N23) |
| chain rule when the similarity variable moves with x and y | C03 | primer (in C03) |
| product rule | C03 | knowledge/primers.md: product rule for differentials (ch01 P38) — reminder |
| 'a function of x times a function of η ⇒ the bracket is a constant' | C03 | C03 (D05 step 10); Ch. 8 §8.2 (recalled) |
| wall conditions (9.20), edge condition (9.21) | C03 | R06, R07 |
| δ → 0 at the leading edge (9.22) | C03 | C03 (N25, D05 step 11) |
| u = ψ_y (9.23), v = −ψ_x (9.24) for the similarity ψ | C03 | C03 (N26, N27, D05) |
| cancellation of the ηf′f″ terms, reduced equation (9.25) | C03 | C03 (N28, N29, D05) |
| δ = √(νx/U) (9.26) | C03 | C03 (N30, D05) |
| Blasius equation (9.27) and its conditions (9.28), (9.29) | C03 | C03 (D05 step 13) |
| memory of the inlet profile (Serrin, Peletier) | C03 | C03 (N24) |
| Blasius' boundary-layer solution versus Blasius' force theorem (6.60) | C03 | C03 (⚠️ callout); Ch. 6 (force theorem) |
| ch09.similarity_reduce_sympy | C03 | C03 (code explain); Ch. 8 §8.4 |
| boundary-value problem with a condition at infinity, truncation at η_max | C04 | C04 (N31, N32) |
| scaling symmetry of an ODE and the Töpfer trick | C04 | primer (in C04) |
| shooting versus boundary-value solving | C04 | primer (in C04) |
| solve_ivp | C04 | knowledge/primers.md: scipy.integrate.solve_ivp (ch01 P31) — reminder |
| solve_bvp | C04 | knowledge/primers.md: scipy.integrate.solve_bvp (ch08 P196) — reminder |
| explicit stepping (RK4 written by hand) | C04 | knowledge/primers.md: explicit stepping (ch01 P30) — reminder |
| η₉₉ as the root of f′ = 0.99 | C04 | C04 (N35, D06 step 4) |
| δ*, θ = 2f″(0), the identities in terms of f | C04 | C04 (N36, D06) |
| wall shear (9.31), skin friction (9.32), drag F_D and C_D (9.33), one side only | C04 | C04 (N37–N40, D06) |
| far-field Gaussian tail and cross-flow (v/U)√Re_x → 0.86 | C04 | C04 (N33, N34) |
| complementary error function erfc | C04 | knowledge/primers.md: complementary error function erfc (ch04 P123) — reminder |
| temporal boundary layer C_f = 1.128/√Re_x (Stokes' first problem) | C04 | Ch. 8 §8.4 (N41 compares) |
| data collapse and the similarity spread | C04 | C04 (N42); Ch. 8 §8.4 |
| BL.blasius_constants, BL.falkner_skan (method toepfer, bvp), blasius_delta99 and friends | C04 | C04 (code explain) |
| power-law outer flow U_e = a xⁿ and the wedge (corner flow Azⁿ) | C05 | C05 (N43); Ch. 6 §6.4 (corner flow, recalled) |
| Falkner–Skan ansatz (9.34) and pressure gradient (9.35) | C05 | C05 (N43, N44, D07) |
| Falkner–Skan equation (9.36), Hartree parameter β = 2n/(n + 1) | C05 | C05 (D07) |
| stagnation-point (Hiemenz) flow, n = 1 | C05 | C05 (N45, worked example) |
| generic thickness δ(x) = √(νx/U_e) | C05 | C05 (N45) |
| continuation in a parameter and a fold (saddle-node) | C05 | primer (in C05) |
| wall curvature f‴(0) = −n | C05 | C05 (N46) |
| separation member n = −0.0904 and the second branch | C05 | C05 (N47) |
| exponent rules | C05 | knowledge/primers.md: exponent rules (ch01 P43) — reminder |
| live widgets | C05 | knowledge/primers.md: live widgets (ch01 P47) — reminder |
| BL.falkner_skan_state, BL.falkner_skan_separation, BL.falkner_skan_fields | C05 | C05 (code explain) |
| von Kármán momentum integral (9.43) | C06 | C06 (D08) |
| conservative (flux) form (9.38), integral of continuity (9.39) | C06 | C06 (N50, N51, D08) |
| Leibniz rule and differentiation under the integral sign | C06 | knowledge/primers.md: differentiation under the integral sign (ch03 P109) — reminder |
| Leibniz rule with a moving limit | C06 | knowledge/primers.md: Leibniz rule with a moving upper limit (ch08 P189) — reminder |
| fundamental theorem of calculus | C06 | knowledge/primers.md: fundamental theorem of calculus (ch02 P84) — reminder |
| integrals that diverge separately but combine (work on [0, h]) | C06 | C06 (D08 step 4) |
| closure problem, assumed profiles (Pohlhausen) | C06 | C06 (N55) |
| BL.momentum_integral_residual, BL.karman_pohlhausen, ch09.momentum_integral_sympy | C06 | C06 (code explain) |
| Holstein–Bohlen parameter λ (9.44), shear correlation l(λ) (9.45), shape factor H(λ) (9.46) | C07 | C07 (N56–N58) |
| closure table from Falkner–Skan solutions (Table 9.1 not reproduced) | C07 | C07 (N59, D10) |
| universal function L(λ) (9.48) and the linear fit (9.49) | C07 | C07 (N62–N64, D09) |
| first-order linear ODE and the integrating factor | C07 | primer (in C07) |
| Thwaites' closed form (9.50), θ₀ and the stagnation-point limit | C07 | C07 (D09) |
| accuracy of Thwaites' method | C07 | C07 (N65) |
| Examples 9.1 (Blasius) and 9.2 (diffuser) | C07 | C07 (N66, N67) |
| integrals of powers of sine (∫sin⁵ by c = cos φ) | C07 | primer (in C07) |
| running integral with cumulative_trapezoid | C07 | knowledge/primers.md: scipy.integrate.cumulative_trapezoid (ch08 P190) — reminder |
| sign-change search with np.sign | C07 | knowledge/primers.md: np.sign and np.nonzero: where a curve crosses zero (ch07 P180) — reminder |
| np.interp | C07 | knowledge/primers.md: np.interp: reading a curve between samples (ch07 P182) — reminder |
| pressure-gradient sign convention (adverse dp/dx > 0) | C07 | C01 (N11) |
| BL.thwaites, BL.thwaites_l, thwaites_H, thwaites_L, thwaites_closure_table, ch09.example_9_1, example_9_2 | C07 | C07 (code explain) |
| wall relation μ u_yy = dp/dx at the wall | C08 | C08 (N70, D12) |
| inflection point | C08 | primer (in C08) |
| accelerating and decelerating streams (9.51), (9.52) | C08 | C08 (N71, N72, D12) |
| Rayleigh inflection criterion (named only; taught in Ch. 11) | C08 | C08 (N72, named) |
| separation as τ₀ = 0 and reverse flow | C08 | C08 (D12 step 6) |
| adverse gradient thickens the layer | C08 | C08 (N73) |
| after separation the layer equations fail (Goldstein singularity, named) | C08 | C08 (N74) |
| BL.wall_curvature, BL.profile_inflection, BL.separation_point | C08 | C08 (code explain) |
| transition of a laminar layer (named; mechanism in Ch. 11) | C09 | C09 (N68) |
| flat-plate drag curve, turbulent branch (named; Ch. 12) | C09 | C09 (N69) |
| form (pressure) drag and friction drag, streamlining | C09 | C09 |
| pressure coefficient C_p = (p − p∞)/(½ρU²) and the ideal value 1 − 4 sin²φ | C09 | Ch. 6 §6.5 (recalled in C09) |
| pressure force on a body, ∮p n dS | C09 | C09 (D13 steps 1–3) |
| Gauss–Legendre quadrature | C09 | knowledge/primers.md: Gauss–Legendre quadrature in 3-D and a smoothed kernel (ch05 P143) — reminder |
| BB.separated_cp, BB.separated_pressure_drag, BB.pressure_drag_from_cp, BL.plate_drag_coefficient, BL.transition_state | C09 | C09 (code explain) |
| internal separation in diffusers (named) | C09 | C09 (N75) |
| cylinder wake regimes (the book's rounded thresholds) | C10 | C10 (N76, N77, N80) |
| Strouhal number (4.102) and shedding frequency | C10 | Ch. 4 §4.11 (recalled in C10, N78) |
| frequency in Hz versus angular frequency in rad/s | C10 | C10 (N78, ⚠️ callout) |
| complex velocity of a point vortex | C10 | Ch. 5 §5.7 (recalled in the primer below) |
| a row of vortices: the cotangent sum | C10 | primer (in C10) |
| linear stability of a steady configuration (perturb, linearise, eigenvalues) | C10 | primer (in C10) |
| eigenvalues and eigenvectors | C10 | knowledge/primers.md: eigenvalues and eigenvectors (ch02 P80) — reminder |
| complex conjugate | C10 | knowledge/primers.md: complex conjugate (ch02 P81) — reminder |
| the complex plane in numpy | C10 | knowledge/primers.md: the complex plane in numpy (ch06 P153) — reminder |
| hyperbolic functions cosh, sinh, tanh | C10 | knowledge/primers.md: hyperbolic functions cosh, sinh, tanh (ch07 P168) — reminder |
| finite-difference Jacobian | C10 | knowledge/primers.md: finite differences (ch01 P21) — reminder |
| Kármán street stability, cosh(πb/a) = √2 | C10 | C10 (D14) |
| vortex-induced vibration and atmospheric streets (named; Ch. 13) | C10 | C10 (N79) |
| BB.karman_street_ratio, karman_street_growth, karman_street_growth_closed, karman_street_spectrum, karman_street_velocity, shedding_frequency, cylinder_flow_regime | C10 | C10 (code explain) |
| drag crisis, critical Reynolds number, roughness | C11 | C11 (N81) |
| sphere flow regimes, Stokes and Oseen asymptotes, Morrison's correlation | C11 | Ch. 8 §8.6 and Ch. 4 §4.11 (recalled in C11, N83) |
| minimum pressure coefficient −5/4 on a sphere | C11 | Ch. 6 §6.8 (recalled in N84) |
| cricket-ball swing, constant lateral force ⇒ parabolic path | C11 | C11 (N84) |
| Magnus lift L = ρUΓ | C11 | Ch. 6 §6.5 (recalled in N85) |
| negative and positive Magnus effect | C11 | C11 (N85) |
| baseball curveball and knuckleball (named) | C11 | C11 (N86) |
| BB.drag_crisis_state, BB.cylinder_state, BB.ball_swing_deflection, BB.magnus_sign | C11 | C11 (code explain) |
| free jet, entrainment, dp/dx = 0 | C12 | C12 (N87) |
| conditions (9.53)–(9.55) | C12 | C12 (N88–N90) |
| momentum flux and mass flux through a cross-section | C12 | primer (in C12) |
| conserved momentum flux (9.56)–(9.58) | C12 | C12 (N91–N93, D15) |
| jet similarity ansatz (9.59)–(9.64) and the exponents −1/3, 2/3 | C12 | C12 (N94–N99, D16) |
| exponent matching from a conserved integral | C12 | knowledge/primers.md: exponent matching for similarity forms (ch08 P197) — reminder |
| sech, arccosh and (tanh)′ = sech² | C12 | primer (in C12) |
| the total-derivative move (look for the derivative of a product) | C12 | primer (in C12) |
| jet ODE 3f‴ + ff″ + f′² = 0 and the first integrals (9.69) | C12 | C12 (N100–N105, D16, D17) |
| the sech² profile (9.71), C = 4√6/3 (9.72), mass flux (9.73), cross-flow (9.74), entrainment (9.75) | C12 | C12 (N106–N111, D17) |
| jet half-width (9.76) and the corrected coefficient | C12 | C12 (N112, D18) |
| jet Reynolds numbers and the instability of the laminar jet (named; Ch. 11, Ch. 12) | C12 | C12 (N113) |
| JET.free_jet, free_jet_constants, free_jet_ode_solve, free_jet_halfwidth, jet_momentum_flux | C12 | C12 (code explain) |
| wall jet, boundary layer under a free jet | C13 | C13 (N114) |
| wall conditions (9.77), far-field condition (9.78) | C13 | R08, R09 |
| integration by parts with a variable lower limit | C13 | primer (in C13) |
| the invariant (9.79)–(9.81): flux of exterior momentum flux | C13 | C13 (N115–N118, D19) |
| exponents u₀ ∝ x^{−1/2}, δ ∝ x^{3/4} (9.82) | C13 | C13 (N119, D19) |
| reduced ODE 4f‴ + ff″ + 2f′² = 0 | C13 | C13 (N120, D20) |
| partial fractions and sympy apart; inverting an implicit solution | C13 | primer (in C13) |
| first integrals of the wall-jet ODE and the implicit solution (9.83) | C13 | C13 (N121, D21) |
| free scale f → λf(λη) and the gauge C f_∞² | C13 | C13 (N123); C04 (scaling symmetry, primer) |
| mass flux (9.84) and the constants (9.85) | C13 | C13 (N122, N123) |
| JET.wall_jet_ode_solve, wall_jet_profile, wall_jet_invariant, wall_jet_integrals, wall_jet_constants | C13 | C13 (code explain) |
| secondary flow through a boundary layer | C14 | C14 |
| radial force balance in a swirl (and a thin layer) | C14 | primer (in C14) |
| centripetal acceleration | C14 | Ch. 4 §4.7 (recalled in the primer) |
| Ekman layer and spin-up (preview) | C14 | C14 (Ch. 13 hook, not derived) |
| river bend (named) | C14 | C14 (gloss) |
| ch09.secondary_flow_radial_force, secondary_flow_layer_profile | C14 | C14 (code explain) |
| exercises, literature (not reproduced) | C14 | S01 |

---
## Part F — derivation storyboards

Builders copy these word for word into `nb.derivation(key, title, goal=…, start=(tex, plain), plan=[…], uses=[…], steps=[dict(did, tex, why, plain)], result=(tex, plain), interpret=…, check=…, check_src=…)` and into the explainer's `derivations: [...]` (phones may shorten *why* to its first sentence; `live`, `set` and `watch` are the explainer's and are listed in Part B). Every step is one move; *why* names the rule and says why we make it (≤ 35 words, ≥ 6); *did* ≤ 8 words. The book's own moves were read on the rendered pages (p392 for D01–D02; p393–p394 for D03; D04 is left to Exercise 9.6;
p396–p398 for D05–D06; p401 for D07; p403 for D08; p405–p406 for D09; D10 and D11 are ours; p410–p411 for D12; D13 and D14 are ours (the book states the street's ratio only); p429 for D15; p429–p430 for D16–D17; p431 for D18; p432–p433 for D19–D21; p435 for D22); the gaps listed in `analysis/ch09.md` §2b are filled and the notebook says so ("the book skips this move; we add it"). Every equation named by number is written out. Colours: pressure orange, viscous rose, inertia teal, velocity blue, momentum
orange, displacement amber. No line of this part starts with a table bar; absolute values are written with \lvert \rvert or in words. Sympy checks use `check_src`, every line commented; ★★★ derivations require it.

### D01 · Scaled equations $u^*u^*_{x^*}+v^*u^*_{y^*}=-p^*_{x^*}+\frac1{Re}u^*_{x^*x^*}+u^*_{y^*y^*}$ (9.7), (9.8), the layer equations (9.9), (9.10) and $\bar\delta/L\sim Re^{-1/2}$ (9.4) — ★★, 14 steps, in C01 (notebook · `bl_scaling_thicknesses`)
- **Goal.** Show how thick a boundary layer is and which terms of the Navier–Stokes equations survive inside it when Re is large.
- **Start.** The x-momentum equation (9.1) $u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=-\frac1\rho\frac{\partial p}{\partial x}+\nu\big(\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}\big)$, continuity (6.2) $\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0$ and the y-momentum equation $u\frac{\partial v}{\partial x}+v\frac{\partial v}{\partial y}=-\frac1\rho\frac{\partial p}{\partial y}+\nu\big(\frac{\partial^2v}{\partial x^2}+\frac{\partial^2v}{\partial y^2}\big)$ (the book never writes the last one; we do) — *in words:* steady plane flow along a thin layer.
- **Plan.** (1) Size the terms to find the thickness. (2) Introduce variables in which everything is of order one. (3) Substitute into all three equations. (4) Divide by the biggest term and read the powers of Re. (5) Drop what vanishes as Re → ∞.
- **Tools.** order-of-magnitude scaling (ch04 P130) · anisotropic scaling with two length scales (ch08 P188) · chain rule (ch01 P49) · dominant balance (ch08 P198) · six Reynolds numbers (primer P200) · sympy check (ch01 P40).
- **Assumptions.** Re ≫ 1; δ̄ ≪ L; the surface radius of curvature ≫ δ̄ (Cartesian axes along the surface); steady flow; every starred quantity O(1) (step 4). Each is used where marked.
- **Steps.**
  1. *did:* Size the derivatives · *tex:* $\partial_x\sim\frac1L,\ \ \partial_y\sim\frac1{\bar\delta},\ \ u\sim U$ (9.5) · *why:* The flow changes over the body length L along the surface but over the thin layer δ̄ across it — the two-length picture of ch08 (8.14); we write the sizes down so terms can be compared. · *plain:* across the layer things change far faster than along it.
  2. *did:* Size v from continuity · *tex:* $\frac UL\sim\frac v{\bar\delta}\ \Rightarrow\ v\sim\frac{\bar\delta}L\,U\ll U$ · *why:* Continuity (6.2) says $u_x=-v_y$, so the two terms have equal size; that fixes v once δ̄ is known (δ̄ ≪ L makes it small). · *plain:* the cross-stream velocity is much smaller than the streamwise one.
  3. *did:* Balance advection against viscosity · *tex:* $\frac{U^2}L\sim\frac{\nu U}{\bar\delta^2}\ \Rightarrow\ \frac{\bar\delta}L\sim Re^{-1/2}$ (9.4) · *why:* Dominant balance (P198): a layer exists where viscosity competes with advection, $uu_x\sim U^2/L$ (9.2) against $\nu u_{yy}\sim\nu U/\bar\delta^2$ (9.3); this equality *defines* δ̄. · *plain:* at high Reynolds number the layer is thin.
  4. *did:* Introduce stretched variables · *tex:* $x^*=\frac xL,\ y^*=\frac{\sqrt{Re}\,y}L,\ u^*=\frac uU,\ v^*=\frac{\sqrt{Re}\,v}U,\ p^*=\frac{p-p_\infty}{\rho U^2}$ (9.6) · *why:* Steps 1–3 make every starred quantity O(1). The pressure is scaled by ρU² because it must balance the advective size U²/L (contrast μU/δ̄ for slow flow, (8.42)). · *plain:* each variable is measured in its natural unit.
  5. *did:* Convert the derivatives · *tex:* $\partial_x=\frac1L\partial_{x^*},\ \ \partial_y=\frac{\sqrt{Re}}L\partial_{y^*},\ \ \partial_{yy}=\frac{Re}{L^2}\partial_{y^*y^*}$ · *why:* Chain rule (P49): $x^*=x/L$, $y^*=y\sqrt{Re}/L$, so each derivative gains the factor of its variable. · *plain:* y-derivatives are magnified by √Re, and second ones by Re.
  6. *did:* Substitute into continuity · *tex:* $\frac UL\Big(\frac{\partial u^*}{\partial x^*}+\frac{\partial v^*}{\partial y^*}\Big)=0$ (8.15) · *why:* $u_x=\frac UL u^*_{x^*}$ and $v_y=\frac U{\sqrt{Re}}\frac{\sqrt{Re}}Lv^*_{y^*}=\frac UL v^*_{y^*}$: the same factor, so no coefficient — the reason v was scaled by Re^{−1/2}. · *plain:* continuity keeps its form.
  7. *did:* Substitute into (9.1): inertia, pressure · *tex:* $\frac{U^2}L\big(u^*u^*_{x^*}+v^*u^*_{y^*}\big)=-\frac{U^2}Lp^*_{x^*}+\cdots$ · *why:* $vu_y=\frac U{\sqrt{Re}}\frac{U\sqrt{Re}}Lv^*u^*_{y^*}$ and $-\frac1\rho p_x=-\frac1\rho\frac{\rho U^2}Lp^*_{x^*}$: all three have the size U²/L. · *plain:* three terms of one size.
  8. *did:* Substitute the viscous terms · *tex:* $\nu u_{xx}+\nu u_{yy}=\frac{\nu U}{L^2}u^*_{x^*x^*}+\frac{\nu U\,Re}{L^2}u^*_{y^*y^*}$ · *why:* Step 5's chain-rule factors: the y-derivative carries an extra Re because y was stretched. · *plain:* y-diffusion is Re times stronger than x-diffusion.
  9. *did:* Divide by U²/L, use ν = UL/Re · *tex:* $u^*u^*_{x^*}+v^*u^*_{y^*}=-p^*_{x^*}+\frac1{Re}u^*_{x^*x^*}+u^*_{y^*y^*}$ (9.7) · *why:* $\frac{\nu U}{L^2}\big/\frac{U^2}L=\frac1{Re}$ and $\frac{\nu U Re}{L^2}\big/\frac{U^2}L=1$. The book prints the denominators as ∂x* and ∂y* without the squares; dimensionally the second derivatives need them. · *plain:* only the x-diffusion term carries 1/Re.
  10. *did:* Substitute into y-momentum: inertia · *tex:* $uv_x+vv_y=\frac{U^2}{L\sqrt{Re}}\big(u^*v^*_{x^*}+v^*v^*_{y^*}\big)$ · *why:* $v_x=\frac U{\sqrt{Re}\,L}v^*_{x^*}$ and $v_y=\frac ULv^*_{y^*}$ (step 6); multiplied by u = Uu* and v = Uv*/√Re each product has size U²/(L√Re). · *plain:* inertia in y is √Re weaker than in x.
  11. *did:* Substitute pressure and viscous terms · *tex:* $-\frac1\rho p_y=-\frac{U^2\sqrt{Re}}Lp^*_{y^*},\ \ \nu v_{xx}=\frac{U^2}{L\,Re^{3/2}}v^*_{x^*x^*},\ \ \nu v_{yy}=\frac{U^2}{L\sqrt{Re}}v^*_{y^*y^*}$ · *why:* $p_y=\rho U^2\frac{\sqrt{Re}}Lp^*_{y^*}$; ν = UL/Re turns the viscous terms into powers of Re. · *plain:* the pressure term is the largest.
  12. *did:* Divide by $U^2\sqrt{Re}/L$ · *tex:* $\frac1{Re}\big(u^*v^*_{x^*}+v^*v^*_{y^*}\big)=-p^*_{y^*}+\frac1{Re^2}v^*_{x^*x^*}+\frac1{Re}v^*_{y^*y^*}$ (9.8) · *why:* Dividing by the largest term gives it coefficient 1 and every other term its power of 1/Re. · *plain:* the pressure derivative stands alone at order one.
  13. *did:* Let Re → ∞ in (9.8) · *tex:* $\frac{\partial p^*}{\partial y^*}=0\ \Rightarrow\ \frac{\partial p}{\partial y}=0$ (9.10) · *why:* Every other term has a coefficient that vanishes as Re → ∞ while its starred factors stay O(1) (assumption of step 4). · *plain:* the pressure does not change across the layer.
  14. *did:* Let Re → ∞ in (9.7), restore units · *tex:* $u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=-\frac1\rho\frac{dp}{dx}+\nu\frac{\partial^2u}{\partial y^2}$ (9.9) · *why:* Drop the 1/Re x-diffusion term; the surviving $u^*_{y^*y^*}$ becomes $\nu u_{yy}$ because $\nu U\,Re/L^2=U^2/L$; p depends on x only by step 13. · *plain:* streamwise diffusion is gone, cross-stream diffusion stays.
- **Result.** $\bar\delta/L\sim Re^{-1/2}$ (9.4), and for Re ≫ 1 the layer obeys (6.2), $uu_x+vu_y=-\frac1\rho\frac{dp}{dx}+\nu u_{yy}$ (9.9) and $\partial p/\partial y=0$ (9.10) — *in words:* thin, with cross-stream friction and a pressure the outer flow imposes.
- **Check.** Units: $\nu U/\bar\delta^2$ = (m²/s)(m/s)/m² = m/s² = U²/L ✓. Air, U = 1 m/s, L = 1 m, ν = 1.5 × 10⁻⁵ m²/s: Re = 6.67 × 10⁴, δ̄/L = 3.87 × 10⁻³, dropped/kept diffusion = 1/Re = 1.5 × 10⁻⁵ ✓. `ch09.bl_nondim_sympy()` returns the coefficient table of steps 9 and 12; the printed (9.7) fails its dimension check.
- **What it means.** Everything in Ch. 9 solves (9.9) with (9.11). It fails when Re_x ≲ 1 (leading edge), when δ is not small against the radius of curvature, and after separation.
- **Traps.** Taking the printed (9.7) at face value (missing squares); dropping 1/Re terms without checking that their factors are O(1); scaling p with μU/δ̄ (the viscous scale of Ch. 8) instead of ρU²; forgetting that the y-equation exists at all.

### D02 · The pressure from the outer flow: $-\frac1\rho\frac{dp}{dx}=U_e\frac{dU_e}{dx}$ (9.11) — ★, 5 steps, in C01 (notebook)
- **Goal.** Get the pressure gradient inside the layer from the ideal outer flow, so the layer equations have no unknown pressure.
- **Start.** Bernoulli's equation along a streamline of the irrotational outer flow $p+\frac12\rho U_e^2=\text{const}$ (Ch. 4 §4.9) and (9.10) $\frac{\partial p}{\partial y}=0$ — *in words:* outside the layer the flow is ideal; inside, the pressure is the same as at the edge.
- **Plan.** (1) Differentiate Bernoulli along the edge. (2) Carry the pressure through the layer with (9.10). (3) Insert it into (9.9).
- **Tools.** Bernoulli's equation (Ch. 4 §4.9) · chain rule (ch01 P49).
- **Assumptions.** Steady, incompressible, irrotational outer flow (used in step 1); (9.10) for the constancy across the layer (step 2).
- **Steps.**
  1. *did:* Write Bernoulli at the layer's edge · *tex:* $p_e(x)+\tfrac12\rho\,U_e^2(x)=\text{const}$ · *why:* Just above the layer the flow is steady, irrotational and inviscid (Prandtl's picture), so Bernoulli holds along the edge streamline. · *plain:* where the outer flow speeds up, the pressure falls.
  2. *did:* Differentiate along x · *tex:* $\frac{dp_e}{dx}+\rho\,U_e\frac{dU_e}{dx}=0$ · *why:* The derivative of a constant is zero; the chain rule (P49) gives $\frac d{dx}\big(\tfrac12U_e^2\big)=U_eU_e'$. · *plain:* the pressure gradient is minus ρ U_e U_e′.
  3. *did:* Rearrange · *tex:* $-\frac1\rho\frac{dp_e}{dx}=U_e\frac{dU_e}{dx}$ (9.11) · *why:* Divide by −ρ; nothing else changes. · *plain:* an accelerating outer flow (U_e′ > 0) has a falling pressure — favourable.
  4. *did:* Carry the pressure across the layer · *tex:* $p(x,y)=p_e(x)\ \Rightarrow\ \frac{\partial p}{\partial x}=\frac{dp_e}{dx}$ · *why:* (9.10) says p does not depend on y inside the layer, so it equals its edge value. · *plain:* the wall feels the edge pressure.
  5. *did:* Insert into (9.9) · *tex:* $u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=U_e\frac{dU_e}{dx}+\nu\frac{\partial^2u}{\partial y^2}$ (9.37, with τ = μ∂u/∂y) · *why:* Replace $-\frac1\rho\frac{dp}{dx}$ by step 3; this is the form every later section uses. · *plain:* the outer speed is the only forcing.
- **Result.** $-\frac1\rho\frac{dp}{dx}=U_eU_e'$ (9.11) — *in words:* the layer is driven only by the outer flow's speed distribution.
- **Check.** Units (m²/s²)/m = m/s² ✓. U_e = 10x s⁻¹ (n = 1, a = 10): at x = 0.1 m U_e = 1 m/s, U_e′ = 10 s⁻¹, dp/dx = −1.2 × 1 × 10 = −12 Pa/m (favourable) ✓ (`outer_flow("wedge", n=1, a=10)`).
- **What it means.** Decelerating outer flow (U_e′ < 0) means dp/dx > 0: adverse — the seed of separation (C08).
- **Traps.** Applying Bernoulli *inside* the layer (it does not hold there: viscous forces act); the sign of dp/dx for a decelerating flow.

### D03 · Displacement thickness $\delta^*=\int_0^\infty\big(1-\frac u{U_e}\big)dy$ (9.16) — ★, 6 steps, in C02 (notebook · `bl_scaling_thicknesses`)
- **Goal.** Define one number that says how far a boundary layer pushes the outer flow away from the wall.
- **Start.** Mass flux through a vertical line from the wall to a height h above the layer, and the same flux if there were no layer — *in words:* the layer carries less mass than an ideal stream.
- **Plan.** (1) Compute the missing mass flux. (2) Replace the layer by a stagnant slab with the same deficit. (3) Read off the thickness. (4) See what it does to the streamlines.
- **Tools.** improper integral of a deficit (primer P202) · definite integral (ch01 P27) · continuity (6.2).
- **Assumptions.** Incompressible flow; the layer has an edge beyond which u = U_e (used in step 5).
- **Steps.**
  1. *did:* Write the real mass flux · *tex:* $\dot m_{real}=\rho\int_0^h u\,dy$ · *why:* Definition of the mass flux per unit span through a vertical line of height h > layer thickness; the wall is at y = 0. · *plain:* what actually flows.
  2. *did:* Subtract it from the ideal flux · *tex:* $\rho\int_0^hU_e\,dy-\rho\int_0^hu\,dy=\rho\int_0^h(U_e-u)\,dy$ · *why:* Without a layer the flux would be ρU_eh; the difference is the mass the layer cannot carry, which the outer flow must divert. · *plain:* the missing mass flux.
  3. *did:* Model the deficit as a stagnant slab · *tex:* $\int_0^h(U_e-u)\,dy\equiv U_e\,\delta^*$ · *why:* Define δ* as the thickness of fluid *at rest* that removes the same flux from a stream at U_e; the outer flow cannot tell the slab from the real layer. · *plain:* an equivalent layer of thickness δ*.
  4. *did:* Solve for δ* · *tex:* $\delta^*=\int_0^h\Big(1-\frac u{U_e}\Big)dy$ · *why:* U_e does not depend on y (outer flow), so divide by it and bring it under the integral. · *plain:* δ* is the area of the deficit 1 − u/U_e.
  5. *did:* Let h → ∞ · *tex:* $\delta^*=\int_0^\infty\Big(1-\frac u{U_e}\Big)dy$ (9.16) · *why:* Above the layer u = U_e, so the integrand is zero there (it decays like a Gaussian for Blasius): enlarging h adds nothing (P202). · *plain:* the integral converges.
  6. *did:* Read the lift of the outer streamlines · *tex:* $v(h)=U\frac{d\delta^*}{dx}\quad(U_e=U)$ · *why:* Continuity: $-v(h)=\int_0^hu_x\,dy=\frac d{dx}\int_0^hu\,dy=\frac d{dx}\big[U(h-\delta^*)\big]=-U\delta^{*\prime}$ (h fixed). · *plain:* the layer acts like a body of thickness δ*: streamlines above rise with slope δ*′.
- **Result.** $\delta^*=\int_0^\infty(1-u/U_e)dy$ (9.16) and, on a flat plate, $v_\infty=U\,d\delta^*/dx$ — *in words:* the layer's slowness displaces the outer flow by δ*.
- **Check.** Units m ✓. u/U = 1 − e^{−y/a}: δ* = a exactly. Blasius: δ* = 1.7208√(νx/U); v_∞ = U d/dx(1.7208√(νx/U)) = 0.8604U/√Re_x ✓ (equal to (9.39)'s v_∞ of C06).
- **What it means.** Wind-tunnel and duct designers enlarge the walls by δ* to keep the core speed; ideal-flow calculations of airfoils are corrected by displacing the body by δ* (the viscous–inviscid iteration, N15).
- **Traps.** Confusing δ* with δ₉₉ (δ* = 0.35 δ₉₉ for Blasius); stopping the integral at δ₉₉ (loses the tail); the sign of v_∞ (outward, positive).

### D04 · Momentum thickness θ (9.17) and $\rho U^2\theta=\int_0^x\tau_0\,dx'$ on a flat plate — ★★, 9 steps, in C02 (notebook · `bl_scaling_thicknesses`)
- **Goal.** Show that θ measures the momentum the wall has removed, so a plate's drag is ρU²θ.
- **Start.** A control volume over the plate: wall below, inflow at x = 0, outflow at x, a streamline as roof — *in words:* the book leaves this to Exercise 9.6; we do it.
- **Plan.** (1) Mass balance fixes the roof height. (2) Momentum balance. (3) Combine the two integrals and recognise θ.
- **Tools.** control volume with a streamline as a side (primer P205) · CV momentum balance (Ch. 4 §4.4) · improper integral of a deficit (primer P202).
- **Assumptions.** Steady, dp/dx = 0 (flat plate) and ∂p/∂y = 0 (9.10), so the pressure force on the volume vanishes (step 4).
- **Steps.**
  1. *did:* Draw the control volume · *tex:* $0\le y\le h(x)$ with the roof the streamline through $(0,h_0)$ · *why:* A streamline carries no mass or momentum flux (P205) and has the ambient pressure on it, so it is a free side of the box. · *plain:* a box whose roof moves with the fluid.
  2. *did:* Balance the mass · *tex:* $\rho\int_0^hu\,dy=\rho\,U\,h_0$ · *why:* Steady flow, nothing through the wall or the roof: what enters at x = 0 (speed U over height h₀) leaves at x. · *plain:* the roof rises as the layer slows the flow.
  3. *did:* Write the net momentum outflow · *tex:* $\rho\int_0^hu^2\,dy-\rho\,U^2h_0$ · *why:* x-momentum leaves through the outlet and enters through the inlet; the roof carries none (streamline). · *plain:* momentum out minus momentum in.
  4. *did:* Write the force on the fluid · *tex:* $F_x=-\int_0^x\tau_0\,dx'$ · *why:* Uniform pressure exerts no net force on a closed surface; the wall pulls the fluid back with the stress τ₀ (Newton's third law). · *plain:* friction is the only force.
  5. *did:* Apply the momentum theorem · *tex:* $\rho\int_0^hu^2dy-\rho U^2h_0=-\int_0^x\tau_0\,dx'$ · *why:* Steady control-volume momentum balance (Ch. 4 §4.4): net outflow of momentum = net force. · *plain:* the fluid loses the momentum the wall takes.
  6. *did:* Eliminate h₀ with step 2 · *tex:* $\rho U^2h_0=\rho U\int_0^hu\,dy$ · *why:* Multiply the mass balance by U so that the inlet flux uses the same integral as the outlet. · *plain:* both fluxes in terms of u(y).
  7. *did:* Combine the integrals · *tex:* $\int_0^x\tau_0\,dx'=\rho\int_0^hu\,(U-u)\,dy$ · *why:* Move to one side: $\rho U\int u-\rho\int u^2=\rho\int u(U-u)$ over the same range. · *plain:* the deficit weighted by u.
  8. *did:* Let the roof rise above the layer · *tex:* $\int_0^x\tau_0\,dx'=\rho U^2\int_0^\infty\frac uU\Big(1-\frac uU\Big)dy$ · *why:* Above the layer u = U, so u(U − u) = 0 and the upper limit may be ∞ (P202); factor out U². · *plain:* only the layer contributes.
  9. *did:* Name the integral · *tex:* $\int_0^x\tau_0\,dx'=\rho U^2\,\theta(x)$ (9.17) · *why:* The dimensionless integral is θ (9.17) $\theta=\int_0^\infty\frac u{U_e}\big(1-\frac u{U_e}\big)dy$: a *length*. · *plain:* the wall's total friction equals the momentum defect ρU²θ.
- **Result.** $\rho U^2\theta(x)=\int_0^x\tau_0\,dx'$ — *in words:* θ is the wall's friction integrated, expressed as a thickness.
- **Check.** Units: ρU²θ = N/m ✓ = ∫τ₀dx. Blasius: ρU²θ = 1.2 × 1 × 2.572 × 10⁻³ = 3.09 × 10⁻³ N/m at x = 1 m and ∫₀¹τ₀dx = 2 × 1.543 × 10⁻³ = 3.09 × 10⁻³ ✓ (τ₀ ∝ x^{−1/2} integrates to 2τ₀(1)).
- **What it means.** Measure the momentum defect far downstream of a body and you have its drag (`ch04.wake_drag_per_span`); (9.43) generalises this to a pressure gradient (C06).
- **Traps.** Forgetting that the roof is a streamline, not a horizontal line, so mass flows out through a horizontal roof; θ is a length, not an angle (Chs. 3, 6, 8).

### D05 · Blasius reduction: $\psi=U\delta f(\eta)$ (9.19) turns (9.18) into $f'''+\frac12ff''=0$ (9.27) with $\delta=\sqrt{\nu x/U}$ (9.26) — ★★, 13 steps, in C03 (notebook · `blasius_similarity_collapse`)
- **Goal.** Show that the flat-plate layer equations collapse to one ODE and find the length δ(x) that makes it so.
- **Start.** $u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=\nu\frac{\partial^2u}{\partial y^2}$ (9.18) with (6.2), conditions $u=v=0$ on y = 0 (9.20) and $u\to U$ (9.21) — *in words:* (9.9) with no pressure gradient.
- **Plan.** (1) Use a stream function. (2) Assume the profile is a stretched copy of one curve. (3) Differentiate. (4) See two terms cancel. (5) Demand that x drops out.
- **Tools.** stream function (Ch. 4 §4.3) · chain rule when η depends on x and y (primer P206) · product rule (ch01 P38) · function of x times function of η (Ch. 8 §8.2) · sympy (ch01 P40).
- **Assumptions.** No imposed length scale; δ → 0 at the leading edge (9.22), used in step 11.
- **Steps.**
  1. *did:* Introduce the stream function · *tex:* $u=\frac{\partial\psi}{\partial y},\qquad v=-\frac{\partial\psi}{\partial x}$ · *why:* Continuity (6.2) then holds automatically (mixed partial derivatives commute), so (9.18) is the only equation left, for one unknown ψ. · *plain:* one function instead of two.
  2. *did:* Choose the similarity form · *tex:* $\psi=U\,\delta(x)\,f(\eta),\qquad\eta=\frac y{\delta(x)}$ (9.19) · *why:* The plate has no length except x, so y can only enter as y/δ(x); ψ has dimensions velocity × length, hence Uδ f. · *plain:* each station is the same shape stretched.
  3. *did:* Differentiate η · *tex:* $\frac{\partial\eta}{\partial y}=\frac1\delta,\qquad\frac{\partial\eta}{\partial x}=-\frac{\eta\,\delta'}{\delta}$ · *why:* Quotient rule on η = y/δ(x): δ depends on x, hence η does too (P206). · *plain:* moving downstream at fixed y changes η.
  4. *did:* Compute u · *tex:* $u=\psi_y=U\delta f'\cdot\frac1\delta=U f'(\eta)$ (9.23) · *why:* Chain rule: ∂/∂y acts on ψ only through η; the δ cancels. · *plain:* u/U = f′(η).
  5. *did:* Compute v · *tex:* $v=-\psi_x=U\delta'\,(\eta f'-f)$ (9.24) · *why:* Product rule: ψ_x = Uδ′f + Uδ f′·η_x = Uδ′f − Uδ′ηf′ (step 3); then change the sign. · *plain:* v is small and comes from the layer growing.
  6. *did:* Differentiate u · *tex:* $u_x=-U\frac{\delta'}\delta\,\eta f'',\quad u_y=\frac U\delta f'',\quad u_{yy}=\frac U{\delta^2}f'''$ · *why:* Chain rule with step 3 ($\eta_x=-\eta\delta'/\delta$, $\eta_y=1/\delta$), applied twice for u_yy. · *plain:* the three derivatives (9.18) needs.
  7. *did:* Form the advective terms · *tex:* $uu_x=-\frac{U^2\delta'}\delta\,\eta f'f'',\quad vu_y=\frac{U^2\delta'}\delta\big(\eta f'f''-ff''\big)$ · *why:* Multiply steps 4–6: $Uf'\cdot(-U\frac{\delta'}\delta\eta f'')$ and $U\delta'(\eta f'-f)\cdot\frac U\delta f''$. · *plain:* two products, each with an η f′ f″ piece.
  8. *did:* Add them · *tex:* $uu_x+vu_y=-\frac{U^2\delta'}{\delta}\,ff''$ · *why:* The two $\eta f'f''$ terms are equal and opposite and cancel (the first and third of the book's three terms). · *plain:* only f f″ survives.
  9. *did:* Insert into (9.18) · *tex:* $-\Big[\frac{U^2\delta'}\delta\Big]ff''=\Big[\frac{\nu U}{\delta^2}\Big]f'''$ (9.25) · *why:* Left side from step 8, right side ν u_yy from step 6; each bracket depends on x only, the rest on η only. · *plain:* an ODE in η whose coefficients depend on x.
  10. *did:* Demand similarity · *tex:* $f'''+\Big[\frac{U\delta\delta'}{\nu}\Big]ff''=0,\qquad\frac{U\delta\delta'}\nu\equiv\frac12$ · *why:* Divide by νU/δ². If the bracket varied with x the ODE would change from station to station; f depends on η alone only if it is constant — we choose ½, which merely defines δ. · *plain:* the coefficient must not depend on x.
  11. *did:* Solve for δ · *tex:* $\frac{d(\delta^2)}{dx}=\frac\nu U\ \Rightarrow\ \delta^2=\frac{\nu x}U+D=\frac{\nu x}U$ (9.26) · *why:* Integrate δδ′ = ν/(2U); δ → 0 as x → 0 (9.22) forces the constant D = 0. · *plain:* δ = √(νx/U), the same root as √(νt) with t = x/U.
  12. *did:* State the Blasius equation · *tex:* $\frac{d^3f}{d\eta^3}+\frac12f\frac{d^2f}{d\eta^2}=0$ (9.27) · *why:* Insert the bracket ½ into step 10. · *plain:* third order, nonlinear, no parameters.
  13. *did:* Translate the conditions · *tex:* $f(0)=f'(0)=0\ \ (9.28),\qquad f'(\infty)=1\ \ (9.29)$ · *why:* No slip u = Uf′(0) = 0 and no through-flow v = Uδ′(0·f′ − f(0)) = 0 from (9.20); u → U (9.21) gives f′ → 1; the leading-edge and edge conditions collapse into one at η = ∞, the consistency test. · *plain:* two conditions at the wall, one at infinity.
- **Result.** $f'''+\frac12ff''=0$ (9.27), $f(0)=f'(0)=0$, $f'(\infty)=1$, $\delta=\sqrt{\nu x/U}$ (9.26) — *in words:* a single curve f′(η) describes the layer at every x.
- **Check.** Units: δ = √((m²/s)(m)/(m/s)) = m ✓; the bracket Uδδ′/ν is dimensionless ✓. sympy (`ch09.similarity_reduce_sympy("blasius")`): residual of (9.18) with the ansatz and (9.27) is 0. Easy numbers: x = 0.5 m, U = 1 m/s, air: δ = 2.74 mm; y = 5.48 mm is η = 2.
- **What it means.** Profiles at different x are the same curve stretched by √x (the collapse figure); the same move works for wedges (D07) and jets (D16, D20).
- **Traps.** Forgetting that η_x ≠ 0 (v would come out wrong); the constant ½ is a *choice* that fixes the meaning of δ (not δ₉₉); demanding proportional brackets is what makes x drop out.

### D06 · Blasius numbers: f″(0) = 0.3321 by the Töpfer trick; η₉₉, δ*, θ = 2f″(0), τ₀ (9.31), C_f (9.32), C_D (9.33) — ★★, 10 steps, in C04 (notebook · `blasius_similarity_collapse`)
- **Goal.** Get every Blasius number from one numerical solution of (9.27) and read off the wall shear and drag.
- **Start.** $f'''+\frac12ff''=0$ (9.27) with $f(0)=f'(0)=0$, $f'(\infty)=1$ (9.28), (9.29) — *in words:* a two-point boundary problem with one unknown, f″(0).
- **Plan.** (1) Turn the problem into an initial-value problem by a scaling symmetry. (2) Read thicknesses from the profile. (3) Convert f″(0) to stress and drag.
- **Tools.** scaling symmetry of an ODE and the Töpfer trick (primer P207) · shooting versus boundary-value solving (primer P208) · `solve_ivp` (ch01 P31) · `brentq` (ch03 P108) · integration by parts.
- **Assumptions.** As D05; drag counted on one side of the plate (used in step 9).
- **Steps.**
  1. *did:* Find the scaling symmetry · *tex:* $f_\lambda(\eta)=\lambda f_1(\lambda\eta)\ \Rightarrow\ f_\lambda'''+\tfrac12f_\lambda f_\lambda''=\lambda^4\big(f_1'''+\tfrac12f_1f_1''\big)$ · *why:* Each term is a product of total weight λ⁴: f′″ brings λ⁴, f f″ brings λ·λ³; so λf₁(λη) solves (9.27) whenever f₁ does. · *plain:* a stretched solution is a solution.
  2. *did:* Integrate with a guessed slope · *tex:* $f_1(0)=f_1'(0)=0,\ f_1''(0)=1\ \Rightarrow\ f_1'(\infty)=2.0854$ · *why:* With all three initial values known this is an initial-value problem (`solve_ivp`); read the far value from the solution. · *plain:* the guess gives the wrong far speed.
  3. *did:* Rescale to satisfy (9.29) · *tex:* $\lambda^2f_1'(\infty)=1\ \Rightarrow\ \lambda=0.6925,\quad f''(0)=\lambda^3=0.33206$ · *why:* f_λ′(∞) = λ²f₁′(∞); choose λ to make it 1; f_λ″(0) = λ³f₁″(0). · *plain:* one number, 0.3321.
  4. *did:* Read the 99 % point · *tex:* $f'(\eta_{99})=0.99\ \Rightarrow\ \eta_{99}=4.910$ · *why:* Root of a monotone function (`brentq`, P108); δ₉₉ = η₉₉√(νx/U). The book prints 4.93, read from a figure: the root is 0.4 % smaller. · *plain:* the edge sits at η = 4.91.
  5. *did:* Displacement thickness · *tex:* $\frac{\delta^*}{\delta}=\int_0^\infty(1-f')\,d\eta=\lim_{\eta\to\infty}(\eta-f)=1.7208$ · *why:* With y = δη, (9.16) is an integral in η; since f′ = df/dη, ∫(1 − f′)dη = η − f. · *plain:* read it off the tail of f.
  6. *did:* Momentum thickness · *tex:* $\frac\theta\delta=\int_0^\infty f'(1-f')\,d\eta=2f''(0)=0.6641$ · *why:* Integrate (9.27) from 0 to η and by parts: f″(0) = ½[I_θ(η) − f(1 − f′)]; the last term → 0 because 1 − f′ decays like a Gaussian. · *plain:* θ is twice the wall slope.
  7. *did:* Wall shear · *tex:* $\tau_0=\mu\Big(\frac{\partial u}{\partial y}\Big)_0=\mu\frac U\delta f''(0)=0.332\,\frac{\rho U^2}{\sqrt{Re_x}}$ (9.31) · *why:* u = Uf′(η), ∂η/∂y = 1/δ, and μU/δ = ρU²/√Re_x since ν = μ/ρ and δ = √(νx/U). · *plain:* wall stress falls like x^{−1/2}.
  8. *did:* Skin friction · *tex:* $C_f=\frac{\tau_0}{\frac12\rho U^2}=\frac{2f''(0)}{\sqrt{Re_x}}=\frac{0.664}{\sqrt{Re_x}}$ (9.32) · *why:* Divide (9.31) by the dynamic pressure ½ρU²; 2f″(0) is the same number as θ/δ. · *plain:* C_f√Re_x = θ/δ.
  9. *did:* Total friction on one side · *tex:* $F_D=\int_0^L\tau_0\,dx=\frac{0.664\,\rho U^2L}{\sqrt{Re_L}}$ · *why:* ∫₀ᴸ x^{−1/2}dx = 2√L, so the integrable singularity at the leading edge costs nothing: the mean is twice the end value. · *plain:* drag ∝ U^{3/2}.
  10. *did:* Drag coefficient · *tex:* $C_D=\frac{F_D}{\frac12\rho U^2L}=\frac{1.328}{\sqrt{Re_L}}$ (9.33) · *why:* Divide by the dynamic pressure and the plate length; it equals 2C_f(L). One side of the plate only. · *plain:* a two-sided plate has twice this.
- **Result.** f″(0) = 0.33206, η₉₉ = 4.910, δ* = 1.7208, θ = 0.6641 (units √(νx/U)); τ₀ = 0.332ρU²/√Re_x, C_f = 0.664/√Re_x, C_D = 1.328/√Re_L — *in words:* one curve gives all of them.
- **Check.** θ = 2f″(0) exactly (0.66411 = 2 × 0.33206). Air, U = 1 m/s, L = 1 m: τ₀(1 m) = 1.543 × 10⁻³ Pa, F_D = 3.09 × 10⁻³ N/m, C_D = 5.14 × 10⁻³; `solve_bvp` and the Töpfer IVP agree to 1e-9; F_D ∝ U^{3/2}.
- **What it means.** The yardstick for all approximate methods (C06–C07) and for every Ch. 10 boundary-layer code.
- **Traps.** Reading η₉₉ as 4.93; a factor 2 between θ/δ and f″(0) (they are equal only in this combination); one side of the plate; naive shooting with η_max large.

### D07 · Falkner–Skan reduction: $\psi=\sqrt{\nu xU_e}f(\eta)$ (9.34) gives $f'''+\frac{n+1}2ff''-nf'^2+n=0$ (9.36) — ★★, 12 steps, in C05 (notebook · `falkner_skan_family`)
- **Goal.** Show that a power-law outer flow keeps the similarity form and find the ODE.
- **Start.** $uu_x+vu_y=U_eU_e'+\nu u_{yy}$ (9.9) with (9.11), $U_e=ax^n$ so $-\frac1\rho\frac{dp}{dx}=U_eU_e'=na^2x^{2n-1}$ (9.35) — *in words:* Blasius' problem with a pressure gradient.
- **Plan.** (1) Generalise the ansatz. (2) Differentiate with two x-dependences. (3) Watch the η terms cancel. (4) Collect.
- **Tools.** chain rule with η(x, y) (primer P206) · power laws (ch01 P43) · sympy (ch01 P40).
- **Assumptions.** U_e = axⁿ exactly; Re_x ≫ 1.
- **Steps.**
  1. *did:* Note the forcing · *tex:* $U_eU_e'=\frac{n\,U_e^2}{x}$ · *why:* U_e′ = naxⁿ⁻¹ = nU_e/x. · *plain:* the pressure term is n times the natural scale U_e²/x.
  2. *did:* Choose the ansatz · *tex:* $\psi=\sqrt{\nu xU_e}\,f(\eta),\quad\eta=\frac y\delta,\quad\delta=\sqrt{\frac{\nu x}{U_e}}$ (9.34) · *why:* Blasius' (9.19) with U replaced by U_e(x): ψ = U_eδ f has the same form; no length other than x is imposed. · *plain:* δ(x) is the local Blasius scale.
  3. *did:* Write powers of x · *tex:* $g\equiv\sqrt{\nu xU_e}\propto x^{\frac{n+1}2},\quad\delta\propto x^{\frac{1-n}2},\quad\frac g\delta=U_e$ · *why:* Exponent rules (P43) with U_e = axⁿ; then $g'=\frac{n+1}{2x}g$ and $\eta_x=-\frac{1-n}{2x}\eta$. · *plain:* both scales are power laws.
  4. *did:* Compute u · *tex:* $u=\psi_y=\frac g\delta f'=U_ef'(\eta)$ · *why:* Chain rule; g/δ = U_e. · *plain:* u/U_e = f′(η).
  5. *did:* Compute v · *tex:* $v=-\psi_x=-\frac g{2x}\big[(n+1)f+(n-1)\eta f'\big]$ · *why:* Product rule on g f(η) with step 3: $-g'f-gf'\eta_x=-\frac{(n+1)g}{2x}f+\frac{(1-n)g}{2x}\eta f'$. · *plain:* v is small, from the growth of the layer.
  6. *did:* Differentiate u in x · *tex:* $u_x=\frac{U_e}x\Big[nf'-\frac{1-n}2\eta f''\Big]$ · *why:* Product rule: U_e′f′ + U_ef″η_x with U_e′ = nU_e/x. · *plain:* two contributions: the outer speed changes and η moves.
  7. *did:* Differentiate u in y · *tex:* $u_y=\frac{U_e}\delta f'',\qquad u_{yy}=\frac{U_e^2}{\nu x}f'''$ · *why:* Chain rule with η_y = 1/δ, and $1/\delta^2=U_e/(\nu x)$. · *plain:* the viscous term in the same units.
  8. *did:* Form u u_x · *tex:* $uu_x=\frac{U_e^2}x\Big[nf'^2-\frac{1-n}2\eta f'f''\Big]$ · *why:* Multiply steps 4 and 6. · *plain:* an n f′² term and an ηf′f″ term.
  9. *did:* Form v u_y · *tex:* $vu_y=-\frac{U_e^2}{2x}\big[(n+1)ff''+(n-1)\eta f'f''\big]$ · *why:* Multiply steps 5 and 7 with g/δ = U_e. · *plain:* an f f″ term and another ηf′f″ term.
  10. *did:* Add the advective terms · *tex:* $uu_x+vu_y=\frac{U_e^2}x\Big[nf'^2-\frac{n+1}2ff''\Big]$ · *why:* The η f′ f″ coefficients are $-\frac{1-n}2$ and $-\frac{n-1}2$, which sum to zero: they cancel. · *plain:* again the ηf′f″ terms vanish.
  11. *did:* Form the right side · *tex:* $U_eU_e'+\nu u_{yy}=\frac{U_e^2}x\big[n+f'''\big]$ · *why:* Step 1 for the pressure term and step 7 for the viscous term. · *plain:* n plus the viscous curvature.
  12. *did:* Equate and divide by U_e²/x · *tex:* $f'''+\frac{n+1}2ff''-nf'^2+n=0$ (9.36) · *why:* Move the advective terms to the right: n f′² − ((n+1)/2) f f″ = n + f‴; all x-dependence has cancelled — similarity holds with (9.28), (9.29). · *plain:* n = 0 is Blasius (9.27); n = 1 is stagnation flow.
- **Result.** $f'''+\frac{n+1}2ff''-nf'^2+n=0$ (9.36) — *in words:* one family of ODEs for every wedge-like outer flow, parameter n.
- **Check.** n = 0: (9.27) ✓. sympy (`ch09.similarity_reduce_sympy("falkner_skan")`) residual 0. n = 1, a = 10 s⁻¹: δ = √(ν/a) = 1.22 mm at every x ✓ (constant).
- **What it means.** The pressure gradient enters only through the constant n; separation shows up as f″(0) → 0 at n = −0.0904 (C05, C08).
- **Traps.** Differentiating η with respect to x once too few times; forgetting that U_e also depends on x; the forcing appears as +n because it is the pressure term moved to the left.

### D08 · The von Kármán momentum integral: $\frac1\rho\tau_0=\frac d{dx}\big[U_e^2\theta\big]+U_e\delta^*\frac{dU_e}{dx}$ (9.43) — ★★, 11 steps, in C06 (notebook · `thwaites_marching`)
- **Goal.** Integrate the layer equation across the layer to get one ODE in x that any layer must obey.
- **Start.** (9.37) $u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=U_e\frac{dU_e}{dx}+\frac1\rho\frac{\partial\tau}{\partial y}$ and (6.2) — *in words:* (9.9) with the pressure from (9.11) and τ = μ∂u/∂y.
- **Plan.** (1) Put the equation in flux form. (2) Integrate over y. (3) Remove v with continuity. (4) Rearrange with the product rule.
- **Tools.** product rule (ch01 P38) · differentiation under the integral sign (ch03 P109) · Leibniz with a moving limit (ch08 P189) · fundamental theorem of calculus (ch02 P84).
- **Assumptions.** Steady, incompressible; τ → 0 outside the layer; u = v = 0 at the wall.
- **Steps.**
  1. *did:* Restate (9.37) · *tex:* $uu_x+vu_y=U_eU_e'+\frac1\rho\frac{\partial\tau}{\partial y}$ (9.37) · *why:* Replace −(1/ρ)dp/dx by U_eU_e′ using (9.11) and νu_yy by (1/ρ)∂τ/∂y with τ = μu_y. · *plain:* the layer equation with the outer forcing.
  2. *did:* Add u times continuity · *tex:* $\frac{\partial(u^2)}{\partial x}+\frac{\partial(uv)}{\partial y}=U_eU_e'+\frac1\rho\frac{\partial\tau}{\partial y}$ (9.38) · *why:* (6.2) says u(u_x + v_y) = 0; adding gives 2uu_x + uv_y + vu_y, which the product rule reads as ∂(u²)/∂x + ∂(uv)/∂y. · *plain:* the flux form.
  3. *did:* Integrate continuity across the layer · *tex:* $\int_0^h\frac{\partial u}{\partial x}dy=-v(h)$, so $-v_\infty$ as h → ∞ (9.39) · *why:* u_x = −v_y and the fundamental theorem (P84): $-[v]_0^h$ with v(0) = 0 (no through-flow). · *plain:* the layer's slowing pushes fluid outward.
  4. *did:* Integrate (9.38) up to a height h · *tex:* $\int_0^h\partial_x(u^2)dy+U_ev(h)-U_eU_e'h=-\frac{\tau_0}\rho$ (9.40) · *why:* Fundamental theorem: [uv]₀ʰ = U_ev(h); τ(h) = 0 outside, τ(0) = τ₀. Keep h finite: ∫u²dy and ∫U_eU_e′dy diverge separately as h → ∞. · *plain:* every term is now a number per unit x.
  5. *did:* Pull d/dx out, eliminate v(h) · *tex:* $\frac d{dx}\int_0^hu^2dy-U_eU_e'h-U_e\int_0^h\frac{\partial u}{\partial x}dy=-\frac{\tau_0}\rho$ (9.41) · *why:* h is a constant, so d/dx passes through the integral (P109); step 3 replaces v(h) by −∫u_xdy. · *plain:* only integrals of u remain.
  6. *did:* Rewrite the last integral · *tex:* $U_e\int_0^hu_x\,dy=\frac d{dx}\Big(U_e\int_0^hu\,dy\Big)-U_e'\int_0^hu\,dy$ · *why:* Product rule (P38): d(U_eG)/dx = U_e′G + U_eG′ with G = ∫u dy. · *plain:* U_e moved inside the derivative.
  7. *did:* Substitute and collect · *tex:* $\frac d{dx}\int_0^h(u^2-U_eu)\,dy+U_e'\int_0^h(u-U_e)\,dy=-\frac{\tau_0}\rho$ (9.42) · *why:* Combine the two d/dx terms; the leftover $-U_eU_e'h+U_e'\int u\,dy$ is $U_e'\int_0^h(u-U_e)dy$. · *plain:* both integrands now vanish above the layer.
  8. *did:* Let h → ∞ · *tex:* the same with upper limits ∞ · *why:* Both integrands are zero where u = U_e, so each integral converges and the limit is allowed (the reason we combined them first). · *plain:* h has dropped out.
  9. *did:* Factor out U_e², flip the sign · *tex:* $\frac{\tau_0}\rho=\frac d{dx}\Big[U_e^2\!\int_0^\infty\!\frac u{U_e}\Big(1-\frac u{U_e}\Big)dy\Big]+U_e'U_e\!\int_0^\infty\!\Big(1-\frac u{U_e}\Big)dy$ · *why:* u² − U_eu = −U_e²(u/U_e)(1 − u/U_e) and u − U_e = −U_e(1 − u/U_e); U_e(x) is independent of y, so it moves in or out of the y-integrals. · *plain:* dimensionless integrals.
  10. *did:* Name θ and δ* · *tex:* $\frac1\rho\tau_0=\frac d{dx}\big[U_e^2\theta\big]+U_e\delta^*\frac{dU_e}{dx}$ (9.43) · *why:* The integrals are the definitions (9.17) and (9.16). · *plain:* wall friction = rate of change of momentum thickness + pressure-gradient term.
  11. *did:* Test on Blasius · *tex:* $\frac{\tau_0}\rho=U^2\frac{d\theta}{dx}=\frac{0.332\,U^2}{\sqrt{Re_x}}$ · *why:* U_e = U constant, θ = 0.6641√(νx/U), dθ/dx = θ/(2x); the result equals (9.31). · *plain:* (9.43) is exact, so it holds for every exact solution.
- **Result.** $\frac1\rho\tau_0=\frac d{dx}[U_e^2\theta]+U_e\delta^*\frac{dU_e}{dx}$ (9.43) — *in words:* one ODE linking θ, δ* and τ₀ for any layer, laminar or time-averaged turbulent; three unknowns need a closure.
- **Check.** Units: τ₀/ρ = m²/s² and U²·dθ/dx = m²/s² ✓. `BL.momentum_integral_residual` ≈ 1e-8 for Blasius and every Falkner–Skan member. At x = 1 m in air τ₀/ρ = 1.286 × 10⁻³ m²/s².
- **What it means.** Thwaites (C07), Ch. 12 (turbulent layers) and Ch. 14 all start here.
- **Traps.** Letting h → ∞ before combining the integrals (each diverges); dropping the −v_∞ term; forgetting τ(∞) = 0.

### D09 · Thwaites' closed form: $\frac{\theta^2U_e^6}\nu=0.45\int_0^xU_e^5dx'+\frac{\theta_0^2U_0^6}\nu$ (9.50) — ★★, 10 steps, in C07 (notebook · `thwaites_marching`)
- **Goal.** Turn the momentum integral (9.43), with an empirical closure, into a formula that gives θ(x) from one integral of the outer speed.
- **Start.** (9.43) $\frac1\rho\tau_0=\frac d{dx}[U_e^2\theta]+U_e\delta^*\frac{dU_e}{dx}$ with the definitions $\lambda=\frac{\theta^2}\nu\frac{dU_e}{dx}$ (9.44), $\tau_0=\mu\frac{U_e}\theta l(\lambda)$ (9.45), $\frac{\delta^*}\theta=H(\lambda)$ (9.46) — *in words:* three unknowns, but shear and shape are assumed to depend on λ alone.
- **Plan.** (1) Make (9.43) dimensionless. (2) Express everything through λ, l and H. (3) Eliminate θ to get an equation for λ. (4) Use Thwaites' linear fit and an integrating factor.
- **Tools.** product rule (ch01 P38) · first-order linear ODE and the integrating factor (primer P210) · sympy `dsolve` (ch01 P40).
- **Assumptions.** l and H are functions of λ only (empirical, checked in D10); the linear fit L ≈ 0.45 − 6λ (step 9).
- **Steps.**
  1. *did:* Restate (9.43) · *tex:* $\frac{\tau_0}\rho=\frac d{dx}\big[U_e^2\theta\big]+U_e\delta^*U_e'$ · *why:* Starting point: one equation, three unknowns (θ, δ*, τ₀). · *plain:* we must close it.
  2. *did:* Multiply by $\rho\theta/(\mu U_e)$ · *tex:* $\frac{\theta\tau_0}{\mu U_e}=\frac{\theta}{\nu U_e}\frac{d(U_e^2\theta)}{dx}+\frac{\theta\delta^*U_e'}\nu$ (9.47) · *why:* With ν = μ/ρ every term becomes dimensionless and θ²/ν appears, the combination in λ. · *plain:* a dimensionless form.
  3. *did:* Expand the derivative · *tex:* $\frac{d(U_e^2\theta)}{dx}=2U_eU_e'\,\theta+U_e^2\,\theta'$ · *why:* Product rule (P38) on U_e²θ. · *plain:* two pieces.
  4. *did:* Collect the terms · *tex:* $\frac{\theta\tau_0}{\mu U_e}=2\frac{\theta^2}\nu U_e'+\frac{U_e\theta\theta'}\nu+\frac{\theta^2}\nu\frac{\delta^*}\theta U_e'$ · *why:* Multiply out $\frac\theta{\nu U_e}(2U_eU_e'\theta+U_e^2\theta')$ and rewrite the last term with δ* = (δ*/θ)θ. · *plain:* every piece has θ²/ν or θθ′.
  5. *did:* Recognise λ, H and θθ′ · *tex:* $\frac{\theta\tau_0}{\mu U_e}=(2+H)\lambda+\frac{U_e}2\frac d{dx}\Big(\frac{\theta^2}\nu\Big)$ · *why:* λ = θ²U_e′/ν (9.44), δ*/θ = H (9.46), and θθ′ = ½(θ²)′ by the chain rule. · *plain:* two terms: pressure gradient and growth.
  6. *did:* Use the shear correlation · *tex:* $l(\lambda)=(2+H)\lambda+\frac{U_e}2\frac d{dx}\Big(\frac{\theta^2}\nu\Big)$ · *why:* (9.45) makes the left side exactly l(λ) (the unnumbered relation the book states without the algebra). · *plain:* l is fixed by λ and the growth of θ².
  7. *did:* Eliminate θ with (9.44) · *tex:* $\frac{\theta^2}\nu=\frac\lambda{U_e'}$ · *why:* Invert the definition of λ, so the only unknown function is λ(x). · *plain:* θ² is λ over the pressure gradient.
  8. *did:* Define the universal function · *tex:* $U_e\frac d{dx}\Big(\frac\lambda{U_e'}\Big)=2l-2(2+H)\lambda\equiv L(\lambda)$ (9.48) · *why:* Multiply step 6 by 2 and insert step 7; l and H depend on λ alone, so the right side is a function of λ only. · *plain:* one equation for λ(x).
  9. *did:* Use Thwaites' linear fit · *tex:* $L\approx0.45-6\lambda\ \Rightarrow\ \frac d{dx}\Big(\frac{\theta^2}\nu\Big)+\frac{6U_e'}{U_e}\frac{\theta^2}\nu=\frac{0.45}{U_e}$ (9.49) · *why:* With Z = θ²/ν = λ/U_e′, (9.48) reads U_eZ′ = 0.45 − 6ZU_e′; divide by U_e. The fit is empirical (D10). · *plain:* a linear ODE for θ²/ν.
  10. *did:* Multiply by the integrating factor $U_e^6$ · *tex:* $\frac d{dx}\Big(U_e^6\frac{\theta^2}\nu\Big)=0.45\,U_e^5\ \Rightarrow\ \frac{\theta^2U_e^6}\nu=0.45\int_0^xU_e^5dx'+\frac{\theta_0^2U_0^6}\nu$ (9.50) · *why:* μ′ = (6U_e′/U_e)μ gives μ = U_e⁶ (P210), turning the left side into one derivative; integrate from x = 0, the constant coming from the initial θ₀ (zero at a stagnation point). · *plain:* θ² is one integral of U_e⁵.
- **Result.** $\theta^2=\frac{0.45\,\nu}{U_e^6(x)}\int_0^xU_e^5dx'+\frac{\theta_0^2U_0^6}{U_e^6(x)}$ (9.50); then λ (9.44), τ₀ (9.45), δ* (9.46) — *in words:* the layer remembers the whole pressure history through ∫U_e⁵.
- **Check.** Units: θ²U_e⁶/ν = m⁶/s⁵ and ∫U_e⁵dx = m⁶/s⁵ ✓. U_e = U: θ² = 0.45νx/U ⇒ θ = 0.6708√(νx/U), 1.0 % above Blasius 0.6641 ✓ (`ch09.example_9_1()`). Diffuser U_e = U₁/(1 + x/L): λ = −(0.45/4)[(1 + x/L)⁴ − 1] ✓.
- **What it means.** One integral, a few lines of arithmetic: predicts shear and separation for any body given its ideal-flow speed.
- **Traps.** Where 0.45 and 6.0 come from (an empirical fit, not derivation); the integrating factor U_e⁶ blows up at a stagnation point, so θ₀ = 0 and the limit θ² ∝ x is needed; l and H are *closures*, not consequences of (9.43).

### D10 · The exact closure on the Falkner–Skan family: λ = nI_θ², l = I_θ f″(0), H = I_δ/I_θ, and L(λ) against 0.45 − 6λ — ★★, 9 steps, in C07 (notebook · `thwaites_marching`)
- **Goal.** Build the functions l(λ) and H(λ) of (9.45), (9.46) from our own solved layers, and see how good Thwaites' linear fit of L(λ) is.
- **Start.** The Falkner–Skan solutions of (9.36) for U_e = axⁿ — *in words:* for these flows every quantity is a function of n only, so it is a function of λ only.
- **Plan.** (1) Express θ, δ*, τ₀ and λ through integrals of f. (2) Tabulate over n. (3) Evaluate L and compare with the line.
- **Tools.** Falkner–Skan solver (C05) · PCHIP interpolation (primer P204) · sympy for the identities (ch01 P40).
- **Assumptions.** Similar flows only; the closure is exact for them and a fit across families.
- **Steps.**
  1. *did:* Use the local scale · *tex:* $\theta=I_\theta\sqrt{\tfrac{\nu x}{U_e}},\ \ \delta^*=I_\delta\sqrt{\tfrac{\nu x}{U_e}}$ · *why:* For Falkner–Skan y = δη with δ = √(νx/U_e), so (9.17) and (9.16) become integrals in η: $I_\theta=\int f'(1-f')d\eta$, $I_\delta=\int(1-f')d\eta$ (as in D06 steps 5–6). · *plain:* thicknesses in units of the local scale.
  2. *did:* Form the shape factor · *tex:* $H=\frac{\delta^*}\theta=\frac{I_\delta}{I_\theta}$ (9.46) · *why:* Ratio of the two integrals; the common scale cancels. · *plain:* H depends on n only.
  3. *did:* Form λ · *tex:* $\lambda=\frac{\theta^2}\nu U_e'=n\,I_\theta^2$ (9.44) · *why:* θ² = I_θ²νx/U_e and U_e′ = nU_e/x, so x and U_e cancel. · *plain:* λ depends on n only.
  4. *did:* Form the wall shear · *tex:* $\tau_0=\mu\frac{U_e}{\delta}f''(0)=\mu\frac{U_e}\theta\,I_\theta f''(0)$ · *why:* The wall slope of u = U_ef′(η) is U_ef″(0)/δ; multiply and divide by θ = I_θδ. · *plain:* the shear in the form of (9.45).
  5. *did:* Read off l · *tex:* $l=\frac{\theta\tau_0}{\mu U_e}=I_\theta\,f''(0)$ (9.45) · *why:* Compare with τ₀ ≡ μ(U_e/θ)l(λ). · *plain:* l depends on n only.
  6. *did:* Tabulate over n · *tex:* $n\ \mapsto\ (\lambda,\ l,\ H)$ · *why:* Solve (9.36) for 60 values from just above the fold (−0.0904) to large n; PCHIP through the rows (P204) gives l(λ), H(λ) on −0.0681 ≤ λ ≤ 0.1065. · *plain:* our own closure table.
  7. *did:* Evaluate the universal function · *tex:* $L=2l-2(2+H)\lambda$ (9.48) · *why:* Apply the definition to each row: Blasius (λ = 0) 2 × 0.2205 = 0.441; n = 1 (λ = 0.0855) 0.7207 − 0.7206 ≈ 0; λ = −0.0675 gives 0.818. · *plain:* a nearly straight function of λ.
  8. *did:* Compare with the line · *tex:* $0.45-6\lambda:\ 0.450,\ -0.063,\ 0.855$ at $\lambda=0,\ 0.0855,\ -0.0675$ · *why:* Thwaites fitted several families (Falkner–Skan, Schubauer's ellipse, Howarth); along Falkner–Skan alone the line is 0.01 too high at Blasius, 0.06 too low for strong acceleration, 0.04 too high near separation. · *plain:* a good fit, not exact.
  9. *did:* Note the two separation criteria · *tex:* $l(\lambda_{sep})=0$: exact family $\lambda=-0.0681$ ($H\approx4.0$); Thwaites' fit $\lambda=-0.09$ · *why:* The attached family ends at the fold n = −0.0904 where f″(0) = 0 (l = 0); the empirical table puts l = 0 at −0.09. We report both. · *plain:* two thresholds, 0.0681 and 0.09.
- **Result.** $\lambda=nI_\theta^2,\ l=I_\theta f''(0),\ H=I_\delta/I_\theta$; $L(\lambda)\approx0.45-6\lambda$ within 0.06 — *in words:* Thwaites' closure is exact for similar flows and a fit for the rest.
- **Check.** Blasius: I_θ = 0.6641, f″(0) = 0.33206 ⇒ l = 0.2205, H = 2.5911; n = 1: λ = 0.0855, l = 0.3603, H = 2.216; `BL.thwaites_closure_table` reproduces them to 1e-9. Book's Table 9.1 (private) differs by < 5 % where it overlaps.
- **What it means.** The closure is why Thwaites' method is good to a few per cent for smooth outer flows and poor for strong acceleration (measured: θ −6 % at n = 1, −8 % at n = 4).
- **Traps.** The fit is *not* exact for Falkner–Skan; two separation thresholds (−0.0681 exact, −0.09 Thwaites); l(λ) is undefined beyond the fold.

### D11 · Thwaites on the cylinder: λ(φ) = 0.45 cos φ F(φ)/sin⁶φ and separation at φ ≈ 103° for the ideal pressure — ★★, 8 steps, in C07 (notebook · `thwaites_marching`)
- **Goal.** Predict where a laminar layer on a cylinder would separate if the outer flow were the ideal one.
- **Start.** (9.50) with θ₀ = 0 and the ideal surface speed $U_e=2U\sin\varphi$ (Ch. 6) — *in words:* the layer starts at the forward stagnation point and is pushed along by the ideal-flow speed.
- **Plan.** (1) Write U_e, U_e′ in φ. (2) Do the integral of U_e⁵. (3) Form λ. (4) Find where it reaches the threshold.
- **Tools.** integrals of powers of sine (primer P211) · substitution in an integral (ch03 P106) · `brentq` (ch03 P108).
- **Assumptions.** Ideal-flow pressure all the way (wrong beyond separation); laminar attached layer; θ₀ = 0 since U_e(0) = 0.
- **Steps.**
  1. *did:* Write the outer flow along the surface · *tex:* $U_e=2U\sin\varphi,\quad x=a\varphi,\quad U_e'=\frac{2U}a\cos\varphi$ · *why:* Ideal-flow surface speed of a cylinder, φ measured from the forward stagnation point, arc length x = aφ. · *plain:* the speed rises to 2U at 90° and falls again.
  2. *did:* Insert into Thwaites · *tex:* $\theta^2=\frac{0.45\,\nu}{U_e^6}\int_0^xU_e^5dx'$ · *why:* (9.50) with θ₀ = 0 because x = 0 is a stagnation point. · *plain:* the layer starts from nothing.
  3. *did:* Change the integration variable · *tex:* $\int_0^xU_e^5dx'=a(2U)^5\!\int_0^\varphi\sin^5\varphi'\,d\varphi'$ · *why:* x′ = aφ′ so dx′ = a dφ′; (2U)⁵ is constant. · *plain:* an integral of sin⁵.
  4. *did:* Evaluate the integral · *tex:* $F(\varphi)=\frac8{15}-\cos\varphi+\frac23\cos^3\varphi-\frac15\cos^5\varphi$ · *why:* sin⁵φ dφ = −(1 − c²)²dc with c = cos φ (P211); integrate from c = 1 to cos φ. · *plain:* a polynomial in cos φ.
  5. *did:* Form λ · *tex:* $\lambda=\frac{\theta^2}\nu U_e'=\frac{0.45\cos\varphi\,F(\varphi)}{\sin^6\varphi}$ · *why:* Insert steps 1–4: 0.45·(2U/a)cos φ · a(2U)⁵F/(2U sin φ)⁶ — U, a and ν all cancel. · *plain:* the same curve for every cylinder, speed and fluid.
  6. *did:* Take the limit at the stagnation point · *tex:* $F\approx\frac{\varphi^6}6,\ \sin\varphi\approx\varphi\ \Rightarrow\ \lambda\to\frac{0.45}6=0.075$ · *why:* For small φ, sin⁵ ≈ φ⁵ so F ≈ φ⁶/6; compare the exact Hiemenz value 0.0855 (D10). · *plain:* Thwaites is 12 % low at the stagnation point.
  7. *did:* Find the crossing · *tex:* $\lambda(\varphi_{sep})=-0.09\ \Rightarrow\ \varphi_{sep}=103.1^\circ\quad(100.9^\circ\text{ for }-0.0681)$ · *why:* `brentq` on λ(φ) between 60° and 120°; λ passes through 0 at 90°, where cos φ = 0 (the speed peaks). · *plain:* it separates about 13° behind the shoulder.
  8. *did:* Compare with observation · *tex:* $\varphi_{sep}^{obs}\approx82^\circ<103^\circ$ · *why:* The calculation uses the ideal pressure; once the wake exists the real pressure differs and the layer lets go earlier (C09). 82° is the book's rounded experimental value. · *plain:* an attached-flow prediction that reality overrules.
- **Result.** $\lambda(\varphi)=\frac{0.45\cos\varphi\,F(\varphi)}{\sin^6\varphi}$; φ_sep = 103.1° (Thwaites) or 100.9° (exact FS) — *in words:* independent of size, speed and viscosity.
- **Check.** λ(30°, 60°, 82°, 90°, 100°) = 0.0722, 0.0589, 0.0263, 0.0000, −0.0603; `BL.thwaites_cylinder` (numerical) equals the closed form to 1e-4.
- **What it means.** The method flags separation but the location from ideal pressure is too late: the reason for the drag-crisis story of C09–C11.
- **Traps.** φ is from the *forward* stagnation point (ch06's θ is from the downstream axis); ideal pressure, not the real one; U_e → 0 at φ = 0 needs the limit.

### D12 · Wall curvature μ u_yy = dp/dx, the inflection points (9.51)–(9.52) and separation as τ₀ = 0 — ★, 6 steps, in C08 (notebook · `falkner_skan_family`)
- **Goal.** Show why a favourable pressure gradient keeps a profile convex-up and an adverse one gives an inflection and finally zero wall shear.
- **Start.** The layer equation (9.9) $u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=-\frac1\rho\frac{dp}{dx}+\nu\frac{\partial^2u}{\partial y^2}$ at the wall — *in words:* what the equation says exactly at y = 0.
- **Plan.** (1) Evaluate at the wall. (2) Compare the curvature at the wall and near the edge. (3) Conclude.
- **Tools.** evaluating a PDE on a boundary · inflection point (primer P212) · intermediate value idea (a continuous function that changes sign has a zero).
- **Assumptions.** Attached laminar layer with a single edge; no suction (v = 0 at the wall).
- **Steps.**
  1. *did:* Evaluate (9.9) at the wall · *tex:* $0=-\frac1\rho\frac{dp}{dx}+\nu\Big(\frac{\partial^2u}{\partial y^2}\Big)_{wall}$ · *why:* No slip (9.12) and no through-flow (9.13) make u = v = 0 at y = 0, so both advective terms vanish. · *plain:* only pressure and friction are left at the wall.
  2. *did:* Rearrange · *tex:* $\mu\Big(\frac{\partial^2u}{\partial y^2}\Big)_{wall}=\frac{dp}{dx}$ · *why:* Multiply by ρ (μ = ρν). · *plain:* the wall curvature is the pressure gradient (over μ).
  3. *did:* Look near the edge · *tex:* $\frac{\partial u}{\partial y}>0\ \text{falls to}\ 0\ \Rightarrow\ \frac{\partial^2u}{\partial y^2}<0\ \ (y\to\delta)$ · *why:* The profile rises to U_e and joins it smoothly, so its slope decreases from positive to zero. · *plain:* near the edge the profile always bends over.
  4. *did:* Accelerating stream · *tex:* $\Big(\frac{\partial^2u}{\partial y^2}\Big)_{wall}<0\quad\text{for }\frac{dp}{dx}<0$ (9.51) · *why:* Step 2 with dp/dx < 0: the same sign as near the edge, so no sign change is required — no inflection point. · *plain:* favourable gradient: one smooth bend.
  5. *did:* Decelerating stream · *tex:* $\Big(\frac{\partial^2u}{\partial y^2}\Big)_{wall}>0\quad\text{for }\frac{dp}{dx}>0$ (9.52) · *why:* Positive at the wall (step 2), negative near the edge (step 3); a continuous u_yy must pass through zero in between: an inflection point at finite y (Ch. 11 links it to instability). · *plain:* adverse gradient: the profile gets an S-bend.
  6. *did:* Follow the wall shear · *tex:* $\tau_0=\mu\Big(\frac{\partial u}{\partial y}\Big)_{wall}\to0$ (separation) · *why:* A stronger or longer adverse gradient keeps flattening the wall slope; when it reaches zero the wall fluid stops, and beyond it reverses. · *plain:* separation is zero shear.
- **Result.** $\mu u_{yy}|_{wall}=dp/dx$ (9.51)–(9.52); adverse ⇒ inflection; separation at τ₀ = 0 — *in words:* the pressure gradient is the wall curvature.
- **Check.** Falkner–Skan: f‴(0) = −n ⇒ u_yy(wall) = −nU_e²/(νx) = (1/μ)dp/dx since dp/dx = −ρnU_e²/x ✓. n = −0.05: f‴(0) = +0.05 and the inflection at η = 1.650 ✓.
- **What it means.** Every adverse-gradient flow (diffusers, airfoils past the suction peak) is at risk; the mechanism returns in Ch. 14 (stall).
- **Traps.** Assuming the wall curvature is zero (only true for Blasius); mixing the sign of dp/dx (adverse > 0).

### D13 · Pressure drag of the separated model: $C_{D,p}=\sin\varphi_s\big(1-\frac43\sin^2\varphi_s-C_b\big)$ — ★, 7 steps, in C09 (notebook · `cylinder_drag_crisis`)
- **Goal.** Show how the position of separation and the wake pressure set a cylinder's form drag.
- **Start.** Ideal pressure $C_p=1-4\sin^2\varphi$ (Ch. 6) up to the separation angle φ_s, and a constant wake pressure C_b behind — *in words:* our simple model; the ideal case has zero drag (d'Alembert).
- **Plan.** (1) The pressure force in the flow direction. (2) A coefficient. (3) Integrate the model.
- **Tools.** force from pressure on a surface · substitution in an integral (ch03 P106) · Gauss–Legendre check (ch05 P143).
- **Assumptions.** Pressure drag only (friction is small on blunt bodies); symmetric separation; constant C_b; φ from the forward stagnation point.
- **Steps.**
  1. *did:* Write the pressure force element · *tex:* $dF_x=-p\,n_x\,dS=p\cos\varphi\;a\,d\varphi$ (per unit span) · *why:* The force on the body is −p n dS; with φ from the forward stagnation point the outward normal has n_x = −cos φ. · *plain:* front pressure pushes downstream, rear pressure upstream.
  2. *did:* Subtract the ambient pressure · *tex:* $dF_x=(p-p_\infty)\cos\varphi\;a\,d\varphi$ · *why:* ∮cos φ dφ = 0, so a uniform p∞ exerts no net force. · *plain:* only differences matter.
  3. *did:* Form a coefficient · *tex:* $C_{D,p}=\frac{F_x}{\tfrac12\rho U^2\,(2a)}=\frac12\oint C_p\cos\varphi\,d\varphi$ · *why:* C_p = (p − p∞)/(½ρU²) and the frontal area per unit span is 2a. · *plain:* an integral of C_p weighted by cos φ.
  4. *did:* State the model · *tex:* $C_p=1-4\sin^2\varphi\ \ (0\le\varphi\le\varphi_s),\qquad C_p=C_b\ \ (\varphi_s\le\varphi\le2\pi-\varphi_s)$ · *why:* Ideal-flow pressure until the layer lets go (C08), then a flat wake pressure; symmetric about the axis (our model, qualitative). · *plain:* the wake does not recover pressure.
  5. *did:* Integrate the front part · *tex:* $2\int_0^{\varphi_s}(1-4\sin^2\varphi)\cos\varphi\,d\varphi=2\Big(\sin\varphi_s-\tfrac43\sin^3\varphi_s\Big)$ · *why:* Substitute s = sin φ (P106): ∫(1 − 4s²)ds = s − 4s³/3; the two halves are equal. · *plain:* the push on the front.
  6. *did:* Integrate the wake part · *tex:* $C_b\int_{\varphi_s}^{2\pi-\varphi_s}\cos\varphi\,d\varphi=-2C_b\sin\varphi_s$ · *why:* ∫cos φ dφ = sin φ, and sin(2π − φ_s) − sin φ_s = −2 sin φ_s. · *plain:* the wake sucks back.
  7. *did:* Add and halve · *tex:* $C_{D,p}=\sin\varphi_s\Big(1-\tfrac43\sin^2\varphi_s-C_b\Big)$ · *why:* ½ × (step 5 + step 6). Ideal limit: φ_s → π with C_b → 1 gives sin φ_s → 0, so C_D,p → 0 (d'Alembert). · *plain:* later separation or a higher wake pressure means less drag.
- **Result.** $C_{D,p}=\sin\varphi_s\big(1-\frac43\sin^2\varphi_s-C_b\big)$ — *in words:* form drag = missing pressure recovery.
- **Check.** φ_s = 90°, C_b = −3 (ideal at 90°): 2.667; C_b = −1: 0.667. φ_s = 82°, C_b = −1.2: 0.884; φ_s = 125°, C_b = −0.6: 0.578 (illustrative pair, ours). The Gauss–Legendre integral of the piecewise C_p equals the formula to 1e-12.
- **What it means.** Explains why streamlining (delaying separation) and the drag crisis (separation 82° → 125°) reduce drag; the numbers are qualitative.
- **Traps.** φ from the *forward* stagnation point (ch06 measured from downstream, D09 there); the crude default C_b = ideal value at φ_s gives 2.59 at 82° — the wake pressure is a *model input*.

### D14 · Kármán street stability: growth rate $\frac{\pi\Gamma}{2a^2}\lvert\frac12-\mathrm{sech}^2\frac{\pi b}a\rvert$, marginal at $\cosh\frac{\pi b}a=\sqrt2$, $b/a=0.2805$ — ★★★, 15 steps, in C10 (notebook · `karman_street_stability`)
- **Goal.** Explain why only one spacing of a staggered double row of point vortices does not grow.
- **Start.** The velocity of a point vortex $w=u-iv=\frac{\Gamma}{2\pi i}\frac1{z-z_0}$ (Ch. 5; Ch. 6 complex potential) — *in words:* the book only states Kármán's result 0.28; we derive it (we call it the closed form of the linear stability).
- **Plan.** (1) Sum the rows. (2) Find the street's own speed. (3) Perturb and linearise. (4) Take the alternating disturbance and reduce to a 4 × 4 matrix. (5) Read its eigenvalues.
- **Tools.** a row of vortices: the cotangent sum (primer P213) · linear stability (primer P214) · eigenvalues (ch02 P80) · complex conjugate (ch02 P81) · hyperbolic functions (ch07 P168) · sympy (ch01 P40).
- **Assumptions.** Inviscid point vortices; infinite rows; small displacements; the alternating disturbance (wavenumber π/a) is the fastest-growing one — confirmed by scanning all k numerically (notebook).
- **Steps.**
  1. *did:* Set up the street · *tex:* $z_A=na+\tfrac{ib}2\ (\Gamma),\qquad z_B=(n+\tfrac12)a-\tfrac{ib}2\ (-\Gamma)$ · *why:* Two rows of counter-rotating vortices, spacing a, separation b, the lower row shifted by a/2 (the staggered street); the facing rows are checked in step 15. · *plain:* two staggered rows.
  2. *did:* Recall one vortex · *tex:* $w=\frac{\Gamma}{2\pi i}\,\frac1{z-z_0},\qquad\frac{dz}{dt}=\bar w$ · *why:* The conjugate velocity of a point vortex (Ch. 5, Ch. 6); the vortex is carried by the local velocity of all the others. · *plain:* each vortex moves with the flow the others make.
  3. *did:* Sum a whole row · *tex:* $w_{row}=\frac{\Gamma}{2ia}\cot\frac{\pi(z-z_0)}a$ · *why:* Σₙ 1/(z − na) = (π/a)cot(πz/a): pairing +n with −n makes the sum converge (primer P213). · *plain:* an infinite row gives a cotangent.
  4. *did:* Find the street's own speed · *tex:* $U_s=\frac\Gamma{2a}\tanh\frac{\pi b}a$ · *why:* Vortex A feels only row B (its own row cancels by symmetry): w = −(Γ/2ia)cot(πζ/a) with ζ = z_A − z_B = −a/2 + ib, and cot(−π/2 + iπb/a) = −i tanh(πb/a); w is real and the same for A and B, so the street translates rigidly. · *plain:* the pattern moves as a whole.
  5. *did:* Perturb and linearise · *tex:* $\dot\delta_i=\sum_{j,n}\frac{\Gamma_j}{2\pi i}\,\frac{\bar\delta_i-\bar\delta_{j,n}}{\overline{(z_i-z_{j,n})}^{\,2}}$ · *why:* Displace each vortex z → z + δ; differentiate w_i = Σ(Γ_j/2πi)/(z_i − z_j) to first order: dw = −Σ(Γ_j/2πi)(δ_i − δ_j)/(z_i − z_j)², and dδ/dt = conj(dw). · *plain:* a linear equation for the displacements.
  6. *did:* Take the alternating disturbance · *tex:* $\delta_{A,n}=d_A(-1)^n,\qquad\delta_{B,n}=d_B(-1)^n$ · *why:* The wavenumber π/a is the most dangerous one (the k-scan in the notebook confirms it for every b/a tried); it reduces the infinite system to two complex unknowns. · *plain:* neighbours move in opposite directions.
  7. *did:* Reduce the sums to lattice sums · *tex:* $P_{ij}=\sum_n\frac1{(\zeta_{ij}-na)^2},\qquad Q_{ij}=\sum_n\frac{(-1)^n}{(\zeta_{ij}-na)^2}$ · *why:* Step 6 in step 5 gives $\dot d_i=\sum_j\frac{\Gamma_j}{2\pi i}\big(\bar d_i\bar P_{ij}-\bar d_j\bar Q_{ij}\big)$; P = π²/(a² sin²(πζ/a)), Q = π² cos/(a² sin²) (P213); for i = j (n ≠ 0): π²/3a² and −π²/6a². · *plain:* closed forms.
  8. *did:* Write the two equations · *tex:* $\dot d_A=\frac\Gamma{2\pi i}\Big[\frac{\pi^2}{2a^2}\bar d_A-\bar P\bar d_A+\bar Q\bar d_B\Big];\ \ \dot d_B=-\frac\Gamma{2\pi i}\Big[\frac{\pi^2}{2a^2}\bar d_B-\bar P\bar d_B+\bar Q\bar d_A\Big]$ · *why:* Own row: P_ii − Q_ii = π²/2a²; the other row has the opposite circulation; P_BA = P_AB and Q_BA = Q_AB because sin² and cos are even. · *plain:* two coupled equations.
  9. *did:* Evaluate at the staggered offset · *tex:* $\sin\frac{\pi\zeta}a=-\cosh\frac{\pi b}a,\ \ \cos\frac{\pi\zeta}a=i\sinh\frac{\pi b}a$ · *why:* ζ = −a/2 + ib: sin(−π/2 + iπb/a) = −cos(iπb/a) = −cosh and cos(−π/2 + iπb/a) = sin(iπb/a) = i sinh. Hence P = π²/(a²cosh²) is real and Q = iπ² sinh/(a²cosh²) imaginary. · *plain:* the lattice sums in hyperbolic functions.
  10. *did:* Insert and simplify · *tex:* $\dot d_A=-\frac{i\pi}2\big[\gamma\bar d_A-i\sigma\bar d_B\big],\ \ \dot d_B=\frac{i\pi}2\big[\gamma\bar d_B-i\sigma\bar d_A\big]$ · *why:* With Γ = a = 1: Γ/2πi = −i/2π; γ = ½ − sech²(πb) and σ = sinh(πb)/cosh²(πb) collect the coefficients (units Γ/a² restored at the end). · *plain:* two numbers γ, σ control everything.
  11. *did:* Split into real and imaginary parts · *tex:* $\dot{\mathbf x}=\frac\pi2\begin{pmatrix}0&-\gamma&-\sigma&0\\-\gamma&0&0&\sigma\\\sigma&0&0&\gamma\\0&-\sigma&\gamma&0\end{pmatrix}\mathbf x,\ \ \mathbf x=(x_A,y_A,x_B,y_B)$ · *why:* d = x + iy, d̄ = x − iy: collect real and imaginary parts of step 10. A real linear system: its eigenvalues decide stability (P214). · *plain:* a 4 × 4 matrix.
  12. *did:* Factorise the characteristic polynomial · *tex:* $\det(M-\mu I)=\big[(\mu-\gamma)^2+\sigma^2\big]\big[(\mu+\gamma)^2+\sigma^2\big]$ · *why:* Direct expansion of the determinant (sympy check): two conjugate quadratics; the ODE eigenvalues are λ = (πΓ/2a²)μ. · *plain:* μ = ±γ ± iσ.
  13. *did:* Read the growth rate · *tex:* $\mathrm{Re}\,\lambda=\frac{\pi\Gamma}{2a^2}\Big\lvert\tfrac12-\mathrm{sech}^2\frac{\pi b}a\Big\rvert$ · *why:* The real part of μ = ±γ ± iσ is ±γ; growth is the positive one, |γ|. · *plain:* it grows unless γ = 0.
  14. *did:* Find the marginal spacing · *tex:* $\gamma=0\ \Leftrightarrow\ \cosh\frac{\pi b}a=\sqrt2\ \Leftrightarrow\ \frac ba=\frac1\pi\cosh^{-1}\sqrt2=0.2805$ · *why:* γ = ½ − 1/cosh² vanishes when cosh² = 2; there σ = ½, so λ = ±i(π/4)Γ/a²: pure oscillation. · *plain:* one neutral spacing.
  15. *did:* Check the facing rows · *tex:* $\det(M-\mu I)=\Big(\mu^2-\frac{\pi^2}{16}\Big)^2\ \Rightarrow\ \mathrm{Re}\,\lambda=\frac{\pi\Gamma}{4a^2}$ · *why:* Repeat steps 7–12 with ζ = ib: P = −π²/(a² sinh²), Q = −π² cosh/(a² sinh²) (real); the polynomial no longer depends on b (sympy check): unstable at every spacing. · *plain:* facing rows never settle.
- **Result.** growth $=\frac{\pi\Gamma}{2a^2}\lvert\frac12-\mathrm{sech}^2\frac{\pi b}a\rvert$, zero only at $b/a=\frac1\pi\cosh^{-1}\sqrt2=0.2805$; facing rows grow at π/4 Γ/a² for every b — *in words:* the street's shape is fixed.
- **Check.** Units 1/s ✓ (Γ/a²). Γ = a = 1: b/a = 0.1 → 0.6400, 0.2 → 0.2982, 0.25 → 0.1099, 0.3 → 0.0663, 0.5 → 0.5359 (all equal to the numerical maximum over k of the full linearisation, 1e-6). At 0.2805 the eigenvalues are ±0.7854 i. Street speed at the stable spacing: tanh(arccosh √2) = 1/√2, U_s = 0.3536 Γ/a.
- **What it means.** Nature's streets have b/a ≈ 0.28 (and St ≈ 0.2, C10); viscosity and finite cores are not in the model, so it explains the shape, not the onset at Re ≈ 40.
- **Traps.** The sum over images converges only in symmetric pairs; both rows are counter-rotating; the growth is zero only at the marginal spacing (a V-shaped curve); the alternating mode is the worst, not the only one; overall factors Γ/a².
- **sympy check intent (★★★).** Rebuild step 11's matrix from the lattice sums (step 9), expand its characteristic polynomial and verify the factorisation (step 12) symbolically, solve γ = 0 (step 14), and repeat with ζ = ib (step 15); cross-check numbers with the notebook's Jacobian of the full periodic velocity. `check_src` sketch:
  ```python
  import sympy as sp                                    # symbolic algebra
  b, mu = sp.symbols('b mu', positive=True)             # b/a (a = Gamma = 1) and the eigenvalue
  c, s = sp.cosh(sp.pi*b), sp.sinh(sp.pi*b)             # cosh, sinh of pi b/a
  gam, sig = sp.Rational(1, 2) - 1/c**2, s/c**2         # step 10: gamma and sigma
  M = sp.Matrix([[0, -gam, -sig, 0], [-gam, 0, 0, sig],
                 [sig, 0, 0, gam], [0, -sig, gam, 0]])  # step 11 (the factor pi/2 is outside)
  p = sp.simplify((M - mu*sp.eye(4)).det())             # characteristic polynomial
  print(sp.simplify(p - ((mu - gam)**2 + sig**2)*((mu + gam)**2 + sig**2)))   # 0: step 12
  print(sp.solve(sp.cosh(sp.pi*b)**2 - 2, b))           # step 14: b = acosh(sqrt 2)/pi
  # facing rows (step 15): P = -pi^2/s^2, Q = -pi^2 c/s^2 rebuilt from the same equations -> (mu^2 - pi^2/16)^2
  ```

### D15 · Jet momentum flux is conserved: $\frac d{dx}\int_{-\infty}^\infty u^2dy=0$ (9.57), $\int u^2dy=J/\rho$ (9.58) — ★★, 7 steps, in C12 (notebook · `free_jet_similarity`)
- **Goal.** Show that a free jet keeps its momentum flux J constant, so we may use it to fix the exponents.
- **Start.** (9.18) $u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=\nu\frac{\partial^2u}{\partial y^2}$ and (6.2), with dp/dx = 0 — *in words:* a layer with no pressure gradient and no wall.
- **Plan.** (1) Make the equation a flux balance. (2) Integrate over all y. (3) Show the boundary terms vanish.
- **Tools.** product rule read backwards (ch01 P38) · fundamental theorem of calculus (ch02 P84) · momentum flux (primer P217).
- **Assumptions.** dp/dx = 0; u and ∂u/∂y → 0 as y → ±∞ (9.53), v bounded.
- **Steps.**
  1. *did:* Multiply continuity by u · *tex:* $u\frac{\partial u}{\partial x}+u\frac{\partial v}{\partial y}=0$ · *why:* (6.2) equals zero, so any multiple is zero; we want a flux form. · *plain:* a harmless zero.
  2. *did:* Add it to (9.18) · *tex:* $2u\frac{\partial u}{\partial x}+u\frac{\partial v}{\partial y}+v\frac{\partial u}{\partial y}=\nu\frac{\partial^2u}{\partial y^2}$ · *why:* Sum of the two equations; the left side now contains total derivatives. · *plain:* three advective pieces.
  3. *did:* Read the left side as derivatives · *tex:* $\frac{\partial(u^2)}{\partial x}+\frac{\partial(uv)}{\partial y}=\frac{\partial}{\partial y}\Big(\nu\frac{\partial u}{\partial y}\Big)$ · *why:* Product rule backwards: 2uu_x = ∂(u²)/∂x and uv_y + vu_y = ∂(uv)/∂y. The right side is the derivative of the kinematic stress ν∂u/∂y (the book's (9.56) writes τ without the 1/ρ, a unit slip). · *plain:* every term is a derivative.
  4. *did:* Integrate over all y · *tex:* $\frac d{dx}\int_{-\infty}^\infty u^2dy+\big[uv\big]_{-\infty}^{\infty}=\Big[\nu\frac{\partial u}{\partial y}\Big]_{-\infty}^{\infty}$ (9.56) · *why:* Fundamental theorem for the y-derivatives; the limits do not depend on x, so d/dx passes through the first integral. · *plain:* boundary terms remain.
  5. *did:* Check the boundary terms · *tex:* $\big[uv\big]_{-\infty}^\infty=0,\qquad\Big[\nu\frac{\partial u}{\partial y}\Big]_{-\infty}^{\infty}=0$ · *why:* u → 0 and ∂u/∂y → 0 far away (9.53) while v stays bounded. · *plain:* nothing leaks out of the sides.
  6. *did:* Read the result · *tex:* $\frac d{dx}\int_{-\infty}^\infty u^2\,dy=0$ (9.57) · *why:* Steps 4 and 5. · *plain:* the momentum flux is the same at every x.
  7. *did:* Name the constant · *tex:* $\int_{-\infty}^\infty u^2dy=\int\tilde u^2dy=\frac J\rho$ (9.58) · *why:* Evaluate at the inlet profile ũ (9.55); J is the momentum flux per unit span, in N/m. · *plain:* J is set by the slot.
- **Result.** $\rho\int u^2dy=J$ at every x — *in words:* the jet keeps its momentum.
- **Check.** Units: ρ∫u²dy = kg/m³ · m²/s² · m = N/m ✓. The notebook integrates u = u₀sech² at four stations: J = 1.0000 N/m each.
- **What it means.** With the conserved J the similarity exponents follow (D16); the mass flux ρ∫u dy is *not* conserved: the jet entrains fluid.
- **Traps.** The printed (9.56) writes τ without 1/ρ; the boundary terms are dropped because u and u_y decay, not because they are zero at the wall (no wall here).

### D16 · Free-jet similarity: exponents from constant J, then $3f'''+ff''+f'^2=0$ with $f'\to0$ (9.66), $f'(0)=1$ (9.67), $f(0)=0$ (9.68) — ★★★, 14 steps, in C12 (notebook · `free_jet_similarity`)
- **Goal.** Show that the free jet is similar, find how fast it slows and widens, and reduce (9.65) to an ODE (the book skips the differentiation).
- **Start.** (9.65) $\frac{\partial\psi}{\partial y}\frac\partial{\partial x}\Big(\frac{\partial\psi}{\partial y}\Big)-\frac{\partial\psi}{\partial x}\frac\partial{\partial y}\Big(\frac{\partial\psi}{\partial y}\Big)=\nu\frac{\partial^2}{\partial y^2}\Big(\frac{\partial\psi}{\partial y}\Big)$ — (9.18) with u = ψ_y, v = −ψ_x — and $\rho\int u^2dy=J$ (9.58) from D15.
- **Plan.** (1) Use the conserved J to fix u₀(x) and δ(x). (2) Differentiate ψ four ways. (3) Watch every power of x cancel. (4) Read the ODE.
- **Tools.** chain rule with x^{−2/3} inside η (primer P206) · exponent rules (ch01 P43) · momentum flux (primer P217) · sympy (ch01 P40).
- **Assumptions.** dp/dx = 0; similarity (no imposed length); J constant (D15).
- **Steps.**
  1. *did:* State the ansatz · *tex:* $\psi=u_0(x)\,\delta(x)\,f(\eta),\quad\eta=\frac y\delta,\quad\delta=\Big[\frac{\nu x}{u_0}\Big]^{1/2}$ (9.59) · *why:* Blasius' (9.19) with the outer speed replaced by the unknown centre-line speed u₀(x); no length except x. · *plain:* one shape, stretched and rescaled.
  2. *did:* Compute u · *tex:* $u=\frac{\partial\psi}{\partial y}=u_0\,f'(\eta)$ (9.60) · *why:* Chain rule; the δ cancels. With f′(0) = 1 (9.67), u₀ is the centre-line speed. · *plain:* u/u₀ = f′(η).
  3. *did:* Write the momentum flux · *tex:* $\frac J\rho=\int u^2dy=u_0^2\,\delta\int f'^2d\eta\equiv u_0^2\,\delta\,C$ (9.61) · *why:* Change variable y = δη (dy = δ dη at fixed x); C = ∫f′²dη is a pure number. · *plain:* J/ρ = u₀²δC.
  4. *did:* Solve for u₀ · *tex:* $u_0^{3/2}\sqrt{\nu x}\,C=\frac J\rho\ \Rightarrow\ u_0=\Big[\frac{J^2}{C^2\rho^2\nu x}\Big]^{1/3}$ (9.62) · *why:* δ = (νx/u₀)^{1/2} makes u₀²δ = u₀^{3/2}(νx)^{1/2}; J is constant (D15), so u₀^{3/2} ∝ x^{−1/2}. · *plain:* u₀ ∝ x^{−1/3}.
  5. *did:* Solve for δ · *tex:* $\delta=\Big[\frac{\nu x}{u_0}\Big]^{1/2}=\Big[\frac{C\rho\nu^2x^2}{J}\Big]^{1/3}$ (9.63) · *why:* Insert step 4: (νx)^{1/2}(νx)^{1/6}(C²ρ²/J²)^{1/6} = (νx)^{2/3}(Cρ/J)^{1/3}. · *plain:* δ ∝ x^{2/3}.
  6. *did:* Write ψ in pure powers of x · *tex:* $\psi=A\,x^{1/3}f(\eta),\ \ \eta=\frac y{Bx^{2/3}},\ \ AB=\nu$ (9.64) · *why:* u₀δ = (Jνx/Cρ)^{1/3}, so A = (Jν/Cρ)^{1/3} and B = (Cρν²/J)^{1/3}; their product is ν exactly. · *plain:* now only x-powers remain.
  7. *did:* Differentiate η · *tex:* $\eta_x=-\frac{2\eta}{3x},\qquad\eta_y=\frac1{Bx^{2/3}}$ · *why:* Power rule on η = y B⁻¹x^{−2/3}. · *plain:* η moves as x^{−2/3}.
  8. *did:* Compute u and v · *tex:* $u=\frac AB\,x^{-1/3}f',\qquad v=-\frac A3\,x^{-2/3}\big(f-2\eta f'\big)$ (9.74) · *why:* u by the chain rule (A x^{1/3}·1/(Bx^{2/3})); v = −ψ_x by the product rule, A[⅓x^{−2/3}f + x^{1/3}f′η_x], with step 7. · *plain:* the second is the entrainment flow.
  9. *did:* Differentiate u in x · *tex:* $u_x=\frac AB\,x^{-4/3}\Big[-\frac{f'}3-\frac{2\eta f''}3\Big]$ · *why:* Product rule on x^{−1/3}f′(η) with η_x from step 7. · *plain:* both the prefactor and η change.
  10. *did:* Differentiate u in y · *tex:* $u_y=\frac A{B^2}\,x^{-1}f'',\qquad u_{yy}=\frac A{B^3}\,x^{-5/3}f'''$ · *why:* Chain rule with η_y (each y-derivative brings 1/(Bx^{2/3})). · *plain:* the viscous term.
  11. *did:* Form the advective terms · *tex:* $uu_x=\frac{A^2}{B^2}x^{-5/3}f'\Big[-\frac{f'}3-\frac{2\eta f''}3\Big],\ \ vu_y=-\frac{A^2}{3B^2}x^{-5/3}\big[ff''-2\eta f'f''\big]$ · *why:* Multiply steps 8–10; both come out ∝ x^{−5/3}. · *plain:* the same power of x in both.
  12. *did:* Add them · *tex:* $uu_x+vu_y=-\frac{A^2}{3B^2}\,x^{-5/3}\big[f'^2+ff''\big]$ · *why:* The ηf′f″ terms have coefficients −⅔ and +⅔ and cancel — as in Blasius and Falkner–Skan. · *plain:* only f′² and f f″ survive.
  13. *did:* Insert into (9.65) and equate · *tex:* $-\frac{A^2}{3B^2}\big[f'^2+ff''\big]=\frac{\nu A}{B^3}f'''$ · *why:* The right side is νu_yy (step 10). Every term carries x^{−5/3}: it cancels, which is the proof that the exponents of steps 4–5 are right. · *plain:* an ODE with constant coefficients.
  14. *did:* Simplify with AB = ν · *tex:* $3f'''+ff''+f'^2=0$ · *why:* Multiply by −3B³/A: 3νf‴ + AB(f′² + ff″) = 0 and AB = ν (step 6). Conditions: f′ → 0 as η → ±∞ (9.66), f′(0) = 1 (9.67), f(0) = 0 (9.68). · *plain:* C dropped out: it is fixed by the solution (D17).
- **Result.** $u_0\propto x^{-1/3}$, $\delta\propto x^{2/3}$ and $3f'''+ff''+f'^2=0$ — *in words:* one conservation law forces the exponents; the ODE has integer coefficients.
- **Check.** Units: J/ρ is m³/s², so A = (Jν/Cρ)^{1/3} is m^{5/3}/s and ψ = A x^{1/3} f is m²/s ✓. `ch09.similarity_reduce_sympy("free_jet")` returns a zero residual and the ODE with coefficients 3, 1, 1. Air, J = 1 N/m: u₀(0.1) = 35.14 m/s, u₀(0.4) = 22.14 m/s (ratio 4^{−1/3} = 0.63) ✓.
- **What it means.** Any conserved integral of the form ∫u² dy forces exponents through dimensional bookkeeping; the ODE then determines the profile (D17).
- **Traps.** ψ has four derivatives to compute (ψ_x, ψ_y, ψ_xy, ψ_yy — here folded into steps 7–11); C is not a free constant of the ODE; AB = ν is what makes the coefficient 3 come out an integer.
- **sympy check intent (★★★).** Build ψ, u, v symbolically with f an undefined function of η(x, y), compute u u_x + v u_y − ν u_yy, divide by the common factor and check that it equals −(A²/3B²)x^{−5/3}(3f‴ + ff″ + f′²)/… i.e. that the reduced ODE is exactly the residual. `check_src` sketch:
  ```python
  import sympy as sp                                         # symbolic algebra
  x, y, J, rho, nu, C = sp.symbols('x y J rho nu C', positive=True)   # positive symbols
  f = sp.Function('f')                                       # the unknown similarity function
  A = (J*nu/(C*rho))**sp.Rational(1, 3)                      # step 6
  B = (C*rho*nu**2/J)**sp.Rational(1, 3)                     # step 6
  eta = y/(B*x**sp.Rational(2, 3))                           # similarity variable
  psi = A*x**sp.Rational(1, 3)*f(eta)                        # stream function (9.64)
  u, v = sp.diff(psi, y), -sp.diff(psi, x)                   # u = psi_y, v = -psi_x
  res = u*sp.diff(u, x) + v*sp.diff(u, y) - nu*sp.diff(u, y, 2)   # (9.18) residual
  print(sp.simplify(res*3*B**2*x**sp.Rational(5, 3)/A**2))   # -(f'^2 + f f'' + 3 f''') evaluated at eta (uses A*B = nu)
  ```

### D17 · Free-jet solution: $f=\sqrt6\tanh(\eta/\sqrt6)$ (9.70), $u=u_0\,\mathrm{sech}^2(\eta/\sqrt6)$ (9.71), $C=\frac{4\sqrt6}3$ (9.72), $\dot m=(36J\rho^2\nu x)^{1/3}$ (9.73), entrainment (9.74)–(9.75) — ★★, 11 steps, in C12 (notebook · `free_jet_similarity`)
- **Goal.** Solve the jet ODE in closed form and read off the profile, the constant C, the mass flux and the entrainment flow.
- **Start.** $3f'''+ff''+f'^2=0$, $f(0)=0$, $f'(0)=1$, $f'\to0$ as η → ±∞ (D16) — *in words:* a third-order ODE that can be integrated twice.
- **Plan.** (1) Spot exact derivatives and integrate twice. (2) Solve the first-order equation. (3) Compute C and the fluxes.
- **Tools.** the total-derivative move (primer P216) · sech, arccosh and (tanh)′ = sech² (primer P215) · separation of variables (ch01 P42) · sympy (ch01 P40).
- **Assumptions.** Symmetric jet; f bounded at infinity.
- **Steps.**
  1. *did:* Restate the ODE and conditions · *tex:* $3f'''+ff''+f'^2=0,\ \ f(0)=0,\ f'(0)=1,\ f'(\pm\infty)=0$ · *why:* The result of D16. · *plain:* what we integrate.
  2. *did:* Spot a total derivative · *tex:* $3f'''+\big(ff'\big)'=0$ · *why:* Product rule backwards (P216): (ff′)′ = f′² + ff″. · *plain:* two pieces are one derivative.
  3. *did:* Integrate once · *tex:* $3f''+ff'=C_1=0$ · *why:* Integrate a derivative equal to zero; far away f″, f′ → 0 with f bounded, so C₁ = 0. · *plain:* a second-order equation.
  4. *did:* Spot another derivative · *tex:* $3f''+ff'=\Big(3f'+\frac{f^2}2\Big)'$ · *why:* d(f²/2)/dη = ff′ and d(3f′)/dη = 3f″. · *plain:* again a total derivative.
  5. *did:* Integrate again · *tex:* $3f'+\frac{f^2}2=3$ (9.69) · *why:* The constant follows from η = 0: f′(0) = 1 and f(0) = 0 give C₂ = 3. · *plain:* a first-order equation.
  6. *did:* Find f at infinity · *tex:* $f'\to0\ \Rightarrow\ f_\infty^2=6,\ \ f_\infty=\sqrt6$ · *why:* Put f′ = 0 in step 5 far from the axis. · *plain:* the stream function saturates at √6.
  7. *did:* Separate variables · *tex:* $f=\sqrt6\,t\ \Rightarrow\ \frac{dt}{1-t^2}=\frac{d\eta}{\sqrt6}$ · *why:* Step 5 gives f′ = 1 − f²/6 = 1 − t²; with f′ = √6 t′ the variables separate (P42). · *plain:* an integrable equation.
  8. *did:* Integrate and apply f(0) = 0 · *tex:* $\tanh^{-1}t=\frac\eta{\sqrt6}\ \Rightarrow\ f=\sqrt6\tanh\frac\eta{\sqrt6},\ \ f'=\mathrm{sech}^2\frac\eta{\sqrt6}$ (9.70), (9.71) · *why:* ∫dt/(1 − t²) = artanh t; t(0) = 0 fixes the constant to 0; (tanh)′ = sech² (P215). · *plain:* the jet profile u = u₀sech²(η/√6).
  9. *did:* Compute the constant C · *tex:* $C=\int f'^2d\eta=\sqrt6\int_{-\infty}^\infty\mathrm{sech}^4s\,ds=\sqrt6\cdot\frac43=\frac{4\sqrt6}3$ (9.72) · *why:* Substitute s = η/√6; ∫sech⁴s ds = ∫(1 − t²)dt over t = tanh s ∈ (−1, 1) = 4/3. · *plain:* C = 3.266.
  10. *did:* Compute the mass flux · *tex:* $\dot m=\rho\int u\,dy=\rho u_0\delta\,[f]_{-\infty}^{\infty}=2\sqrt6\,\rho u_0\delta=(36J\rho^2\nu x)^{1/3}$ (9.73) · *why:* ∫f′dη = f(∞) − f(−∞) = 2√6; u₀δ = (Jνx/Cρ)^{1/3} (D16) so ṁ³ = (2√6)³ρ³Jνx/(Cρ) = 36 Jρ²νx since (2√6)³/C = 117.6/3.266 = 36. · *plain:* ṁ ∝ x^{1/3}: entrainment.
  11. *did:* Compute the cross-flow · *tex:* $\frac v{u_0}=-\frac{f-2\eta f'}{3\sqrt{Re_x}}\ \to\ \mp\frac{\sqrt6}{3\sqrt{Re_x}}\quad(\eta\to\pm\infty)$ (9.74), (9.75) · *why:* v from D16 step 8 divided by u₀, with δ/x = 1/√Re_x (Re_x = xu₀/ν); as η → ±∞, f → ±√6 and ηf′ → 0. · *plain:* ambient fluid flows toward the jet from both sides.
- **Result.** $u=u_0\,\mathrm{sech}^2(\eta/\sqrt6)$, $C=\frac{4\sqrt6}3$, $\dot m=(36J\rho^2\nu x)^{1/3}$, $v/u_0\to\mp\sqrt6/(3\sqrt{Re_x})$ — *in words:* Bickley's jet, growing by entrainment.
- **Check.** Bickley's coefficients: u₀ = 0.4543(J²/ρ²νx)^{1/3} (= (1/C²)^{1/3}), η-scaling 0.2752, ṁ coefficient 36^{1/3} = 3.3019 ✓. sympy: f = √6 tanh(η/√6) solves 3f‴ + ff″ + f′² = 0 identically; `solve_bvp` matches to 1e-8. Air J = 1, x = 0.1 m: ṁ = 0.04268 kg/(m s).
- **What it means.** Momentum flux is conserved, mass flux is not; the same two facts hold for round and turbulent jets (Ch. 12).
- **Traps.** C₁ = 0 needs f″, f′ → 0 with f bounded; ∫sech⁴ = 4/3 (not 2); √6 appears from f(∞).

### D18 · The jet half-width $h_{99}=7.3319\big[\frac{C\rho\nu^2x^2}J\big]^{1/3}$ (9.76 corrected) — ★, 5 steps, in C12 (notebook · `free_jet_similarity`)
- **Goal.** Find where the jet speed has fallen to 1 % of the centre-line value, and correct the book's number.
- **Start.** $u=u_0\,\mathrm{sech}^2(\eta/\sqrt6)$ (9.71) and δ from (9.63) — *in words:* the profile decays like a bell.
- **Plan.** (1) Set sech² equal to 0.01. (2) Invert with arccosh. (3) Multiply by √6 δ.
- **Tools.** sech, arccosh and (tanh)′ = sech² (primer P215).
- **Assumptions.** The half-width is defined at the 1 % level (any level works the same way).
- **Steps.**
  1. *did:* Impose the level · *tex:* $\mathrm{sech}^2\Big(\frac{h_{99}}{\sqrt6\,\delta}\Big)=0.01$ · *why:* Definition of h₉₉ (u = 1 % of u₀), analogous to δ₉₉. · *plain:* where the speed is 1 % of the peak.
  2. *did:* Take the square root · *tex:* $\mathrm{sech}\,z=0.1\ \Rightarrow\ \cosh z=10,\quad z=\frac{h_{99}}{\sqrt6\,\delta}$ · *why:* sech² = 0.01 gives sech = 0.1 (positive branch), and cosh = 1/sech. · *plain:* cosh z = 10.
  3. *did:* Invert · *tex:* $z=\cosh^{-1}10=\ln\big(10+\sqrt{99}\big)=2.9932$ · *why:* arccosh y = ln(y + √(y² − 1)) for y ≥ 1 (P215). · *plain:* z = 2.993.
  4. *did:* Restore units · *tex:* $h_{99}=\sqrt6\,z\,\delta=7.3319\Big[\frac{C\rho\nu^2x^2}J\Big]^{1/3}$ · *why:* Multiply by √6δ with δ from (9.63): √6 × 2.9932 = 7.3319. · *plain:* width ∝ x^{2/3}.
  5. *did:* Explain the printed value · *tex:* $\cosh^{-1}5=2.2924,\ \ \sqrt6\times2.2924=5.6152$ · *why:* The book's 2.2924 satisfies sech²z = 1/25 = 0.04: it is the 4 % point, not the 1 % point. · *plain:* 5.6152 belongs to 4 %.
- **Result.** $h_{99}=7.3319\,\delta$; the book's 5.6152 δ is $h_{04}$ — *in words:* corrected half-width; Re_h99 inherits the correction.
- **Check.** 1/cosh²(2.2924) = 0.0400 ✓; 1/cosh²(2.9932) = 0.0100 ✓. Air, J = 1, x = 0.1 m: h₉₉ = 7.3319 × 0.2066 mm = 1.515 mm (printed would give 1.160 mm, 23 % too small).
- **What it means.** Level matters: state it (1 %, 4 %, 50 %): the coefficient is √6 arccosh(1/√level) (0.5 → 2.159).
- **Traps.** sech vs sech²; mixing the 1 % and 4 % levels.

### D19 · The wall-jet invariant: $\frac d{dx}\int_0^\infty u\Big(\int_y^\infty u^2dy'\Big)dy=0$ (9.80), hence $x\,u_0^2=$ const (9.81) — ★★★, 13 steps, in C13 (notebook · `wall_jet_invariant`)
- **Goal.** Find the quantity a wall jet conserves, now that the wall removes ordinary momentum flux.
- **Start.** (9.18) $u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=\nu\frac{\partial^2u}{\partial y^2}$, (6.2), $u=v=0$ on y = 0 (9.77) and $u\to0$ as y → ∞ (9.78) — *in words:* the book gives four displayed integrals; we fill every move.
- **Plan.** (1) Integrate (9.18) from y to ∞ so the viscous term becomes a stress. (2) Multiply by u and integrate again. (3) Use continuity and integration by parts to turn the result into a total derivative.
- **Tools.** integration by parts with a variable lower limit (primer P218) · fundamental theorem of calculus (ch02 P84) · product rule (ch01 P38) · differentiation under the integral sign (ch03 P109).
- **Assumptions.** dp/dx = 0; no slip and no through-flow at the wall; u decays fast enough at infinity (exponentially, D21 step 13).
- **Steps.**
  1. *did:* Collect the equations · *tex:* $uu_x+vu_y=\nu u_{yy},\quad u_x+v_y=0,\quad u=v=0\ (y=0),\ u\to0\ (y\to\infty)$ · *why:* (9.18), (6.2), (9.77), (9.78). · *plain:* the whole problem.
  2. *did:* Integrate (9.18) from y to ∞ · *tex:* $\int_y^\infty(uu_x+vu_y)\,dy'=-\nu\,u_y(y)$ · *why:* Fundamental theorem: ν[u_y]_y^∞ = −νu_y(y) since u_y → 0 at infinity. · *plain:* the viscous term becomes a stress at height y.
  3. *did:* Write uu_x as a derivative · *tex:* $\frac\partial{\partial x}\int_y^\infty\frac{u^2}2\,dy'+\int_y^\infty vu_y\,dy'=-\nu u_y$ · *why:* uu_x = ∂(u²/2)/∂x; the limits y, ∞ do not depend on x, so ∂/∂x passes outside. · *plain:* a flux term and a cross-flow term.
  4. *did:* Multiply by u and integrate over the wall layer · *tex:* $T_1+T_2=-\nu\int_0^\infty uu_y\,dy$, with $T_1=\int_0^\infty u\,\frac\partial{\partial x}\Big(\int_y^\infty\frac{u^2}2dy'\Big)dy,\ T_2=\int_0^\infty u\Big(\int_y^\infty vu_y\,dy'\Big)dy$ · *why:* The equation of step 3 holds at every y; multiplying by u and integrating over 0…∞ keeps it true. · *plain:* two double integrals.
  5. *did:* Evaluate the right side · *tex:* $-\nu\int_0^\infty uu_y\,dy=-\frac\nu2\big[u^2\big]_0^\infty=0$ · *why:* uu_y = ∂(u²/2)/∂y; u = 0 at the wall (9.77) and at infinity (9.78). · *plain:* friction drops out entirely.
  6. *did:* Rewrite vu_y with continuity · *tex:* $vu_y=\frac{\partial(uv)}{\partial y}-uv_y=\frac{\partial(uv)}{\partial y}+uu_x$ · *why:* Product rule backwards, then v_y = −u_x from (6.2). · *plain:* a derivative plus a term we already have.
  7. *did:* Integrate the inner integral of T₂ · *tex:* $\int_y^\infty vu_y\,dy'=-u\,v(y)+\frac\partial{\partial x}\int_y^\infty\frac{u^2}2dy'$ · *why:* Fundamental theorem on ∂(uv)/∂y gives [uv]_y^∞ = −u(y)v(y) (uv → 0 at infinity); the other piece is the flux term of step 3. · *plain:* T₂ contains T₁'s integrand.
  8. *did:* Insert into T₂ · *tex:* $T_2=-\int_0^\infty u^2v\,dy+T_1$ · *why:* Multiply step 7 by u and integrate over y. · *plain:* T₂ = T₁ minus a cross-flow term.
  9. *did:* Use T₁ + T₂ = 0 · *tex:* $\int_0^\infty u\,\frac\partial{\partial x}\Big(\int_y^\infty u^2dy'\Big)dy-\int_0^\infty u^2v\,dy=0$ (9.79) · *why:* Steps 4–5 give T₁ + T₂ = 0; with T₂ = T₁ − ∫u²v dy (step 8) this is 2T₁ − ∫u²v dy = 0, and 2T₁ is the first integral above (the ½ of u²/2 is absorbed). · *plain:* the book's (9.79).
  10. *did:* Differentiate the double integral · *tex:* $\frac d{dx}\int_0^\infty uG\,dy=\int_0^\infty u_xG\,dy+\int_0^\infty u\,G_x\,dy,\quad G(y)\equiv\int_y^\infty u^2dy'$ · *why:* Product rule under the integral sign (P109); limits 0, ∞ are constants. · *plain:* two pieces.
  11. *did:* Transform the first piece · *tex:* $\int_0^\infty u_xG\,dy=-\int_0^\infty v_yG\,dy=-\big[vG\big]_0^\infty+\int_0^\infty vG_y\,dy=-\int_0^\infty u^2v\,dy$ · *why:* u_x = −v_y (6.2); integrate by parts (P218): [vG]₀^∞ = 0 since v(0) = 0 (9.77) and G(∞) = 0; G_y = −u². · *plain:* the cross-flow term of step 9.
  12. *did:* Add the pieces · *tex:* $\frac d{dx}\int_0^\infty u\Big(\int_y^\infty u^2dy'\Big)dy=-\int_0^\infty u^2v\,dy+\int_0^\infty u\,G_x\,dy=0$ (9.80) · *why:* Step 11 for the first piece; the sum vanishes by (9.79). · *plain:* the flux of exterior momentum flux is conserved.
  13. *did:* Insert the similarity form · *tex:* $\nu x\,u_0^2\!\int_0^\infty\!\Big(f'\!\int_\eta^\infty\! f'^2d\eta'\Big)d\eta=\text{const}\ \Rightarrow\ u_0=C\,x^{-1/2}$ (9.81), (9.82) · *why:* With u = u₀f′(η), δ² = νx/u₀ (9.59): u∫u² dy′ = u₀³δ²·I = νxu₀²·I; the double integral I of f is a number, so x u₀² = C². · *plain:* u₀ ∝ x^{−1/2}.
- **Result.** $\int_0^\infty u(\int_y^\infty u^2dy')dy$ is independent of x and $xu_0^2=C^2$ — *in words:* a conserved 'momentum flux of the momentum flux' fixes the exponent.
- **Check.** Units: u³L² = m⁴/s³ (Ψ, not a force per length). f_∞ = 1: the double integral is 1/40 (notebook, 1e-9); `JET.wall_jet_invariant` on the solution at four stations is constant to 1e-9 while ρ∫u²dy falls like x^{−1/4}.
- **What it means.** When the obvious invariant is lost to the wall, the next-order one survives; the exponent 3/4 (D20) follows.
- **Traps.** The ordinary momentum flux is *not* conserved (wall shear); the boundary terms in steps 5, 7, 11 need the wall conditions and the decay at infinity; the factor 2 between T₁ and the integrand of (9.79).
- **sympy check intent (★★★).** Verify the identity of step 6 and the by-parts identity of step 11 on a generic decaying pair (u, v) satisfying continuity via a stream function ψ(x, y) = e^{−y}·(something) (u = ψ_y, v = −ψ_x), and check that the double integral's x-derivative vanishes for the exact similarity solution built from `JET.wall_jet_profile`. `check_src` sketch:
  ```python
  import sympy as sp                                      # symbolic algebra
  x, y, yp = sp.symbols('x y yp', positive=True)          # streamwise, wall-normal and dummy variable
  psi = x**sp.Rational(1, 4)*sp.exp(-y/x**sp.Rational(3, 4))*(1 - sp.exp(-y/x**sp.Rational(3, 4)))   # a smooth test stream function, wall at 0
  u, v = sp.diff(psi, y), -sp.diff(psi, x)                # u = psi_y, v = -psi_x: continuity holds automatically
  G = sp.integrate((u**2).subs(y, yp), (yp, y, sp.oo))   # G(y) = int_y^inf u^2 dy'
  lhs = sp.integrate(sp.diff(u, x)*G + u*sp.diff(G, x), (y, 0, sp.oo))   # d/dx of the double integral, term by term (step 10)
  rhs = -sp.integrate(u**2*v, (y, 0, sp.oo)) + sp.integrate(u*sp.diff(G, x), (y, 0, sp.oo))   # step 11 substituted
  print(sp.simplify(lhs - rhs))                           # 0: the by-parts identity (needs v = 0 at the wall and G -> 0 at infinity; if sympy cannot integrate symbolically, evaluate at x = 1 with quad)
  ```

### D20 · Wall-jet reduction: $\psi=[\nu Cx^{1/2}]^{1/2}f(\eta)$ (9.82) turns (9.65) into $4f'''+ff''+2f'^2=0$ (the book prints coefficient 1) — ★★★, 12 steps, in C13 (notebook · `wall_jet_invariant`)
- **Goal.** Derive the ODE of the wall jet, including the coefficient the book prints wrongly.
- **Start.** (9.65) with u = ψ_y, v = −ψ_x and the similarity form (9.59) with u₀ = Cx^{−1/2} from D19 — *in words:* the same algebra as D16 with different exponents.
- **Plan.** (1) Write ψ in powers of x. (2) Differentiate four ways. (3) See the ηf′f″ terms cancel. (4) Compare coefficients.
- **Tools.** chain rule with x^{−3/4} in η (primer P206) · exponent rules (ch01 P43) · sympy (ch01 P40).
- **Assumptions.** Boundary-layer equations, dp/dx = 0; u₀ = Cx^{−1/2} (D19).
- **Steps.**
  1. *did:* Write ψ in powers of x · *tex:* $\psi=A\,x^{1/4}f(\eta),\ \ \eta=\frac y{Bx^{3/4}},\ \ A=\sqrt{\nu C},\ B=\sqrt{\nu/C}$ (9.82) · *why:* u₀δ = Cx^{−1/2}(νx^{3/2}/C)^{1/2}: δ = (νx/u₀)^{1/2} = Bx^{3/4}; note A/B = C, AB = ν. · *plain:* pure powers of x again.
  2. *did:* Differentiate η · *tex:* $\eta_x=-\frac{3\eta}{4x},\qquad\eta_y=\frac1{Bx^{3/4}}$ · *why:* Power rule on η = yB⁻¹x^{−3/4}. · *plain:* η moves as x^{−3/4}.
  3. *did:* Compute u · *tex:* $u=\frac AB\,x^{-1/2}f'=Cx^{-1/2}f'(\eta)$ · *why:* Chain rule: A x^{1/4}·1/(Bx^{3/4}); A/B = C. · *plain:* u₀ = Cx^{−1/2}, as D19 demanded.
  4. *did:* Compute v · *tex:* $v=-\frac A4\,x^{-3/4}\big(f-3\eta f'\big)$ · *why:* v = −ψ_x by the product rule: −A[¼x^{−3/4}f + x^{1/4}f′η_x] with step 2. · *plain:* cross-flow, ∝ x^{−3/4}.
  5. *did:* Differentiate u in x · *tex:* $u_x=Cx^{-3/2}\Big[-\frac{f'}2-\frac{3\eta f''}4\Big]$ · *why:* Product rule on x^{−1/2}f′(η) with η_x from step 2. · *plain:* both the prefactor and η change.
  6. *did:* Differentiate u in y · *tex:* $u_y=\frac C B\,x^{-5/4}f'',\qquad u_{yy}=\frac C{B^2}\,x^{-2}f'''=\frac{C^2}\nu x^{-2}f'''$ · *why:* Chain rule with η_y; 1/B² = C/ν. · *plain:* the viscous term.
  7. *did:* Form uu_x · *tex:* $uu_x=C^2x^{-2}f'\Big[-\frac{f'}2-\frac{3\eta f''}4\Big]$ · *why:* Multiply steps 3 and 5. · *plain:* f′² and ηf′f″.
  8. *did:* Form vu_y · *tex:* $vu_y=-\frac{C^2}4\,x^{-2}\big[ff''-3\eta f'f''\big]$ · *why:* Multiply steps 4 and 6; AC/B = C². · *plain:* ff″ and ηf′f″.
  9. *did:* Add them · *tex:* $uu_x+vu_y=-\frac{C^2}4\,x^{-2}\big[2f'^2+ff''\big]$ · *why:* The ηf′f″ coefficients −¾ and +¾ cancel; the f′² coefficient is −½ = −2/4. · *plain:* both sides scale as x^{−2}.
  10. *did:* Form the viscous term · *tex:* $\nu u_{yy}=C^2x^{-2}f'''$ · *why:* Step 6. · *plain:* same prefactor C²x^{−2}.
  11. *did:* Equate and compare coefficients · *tex:* $-\frac14\big[2f'^2+ff''\big]=f'''\ \Rightarrow\ 4f'''+ff''+2f'^2=0$ · *why:* Divide by C²x^{−2} (all powers of x cancel: the check that exponents 1/2 and 3/4 are right) and multiply by −4. The book prints f‴ + ff″ + 2f′² = 0: the coefficient of f‴ is 4, as its own next line 4ff″ − 2f′² + f²f′ = 0 requires. · *plain:* the factor 4 matters.
  12. *did:* State the conditions · *tex:* $f(0)=0,\quad f'(0)=0,\quad f'(\infty)=0$ · *why:* (9.77): u = v = 0 at y = 0 gives f′(0) = 0 and f(0) = 0; (9.78): u → 0 far away. All three are homogeneous, so f has a free scale (D21). · *plain:* no normalisation condition.
- **Result.** $4f'''+ff''+2f'^2=0$ with $f(0)=f'(0)=f'(\infty)=0$; δ = (νx^{3/2}/C)^{1/2} ∝ x^{3/4} — *in words:* thicker than the free jet (x^{2/3}).
- **Check.** sympy: substituting ψ leaves residual −(C²/4x²)(4f‴ + ff″ + 2f′²); the printed coefficient leaves the extra −(3/4)C²f‴/x². `JET.wall_jet_ode_solve(coeff=1.0)` ends at f_∞ = 0.397 instead of 1 for the same f″(0) — the printed ODE fails the (9.83) relation.
- **What it means.** The invariant and the ODE together fix exponents and profile; the wall jet spreads faster (3/4) than the free jet (2/3) because the wall slows the fluid.
- **Traps.** The coefficient 4 (not 1); differentiating η with x^{−3/4}; the y-dependence must cancel: only the x-powers may cancel, not the f-terms.
- **sympy check intent (★★★).** Substitute ψ symbolically, form the residual of (9.18), divide by C²x^{−2} and show it equals −(4f‴ + ff″ + 2f′²)/4; then repeat with coefficient 1 to see the leftover. `check_src` sketch:
  ```python
  import sympy as sp                                       # symbolic algebra
  x, y, nu, C = sp.symbols('x y nu C', positive=True)      # positive symbols
  f = sp.Function('f')                                     # the similarity function
  A, B = sp.sqrt(nu*C), sp.sqrt(nu/C)                      # step 1
  eta = y/(B*x**sp.Rational(3, 4))                         # similarity variable
  psi = A*x**sp.Rational(1, 4)*f(eta)                      # stream function (9.82)
  u, v = sp.diff(psi, y), -sp.diff(psi, x)                 # velocities
  res = u*sp.diff(u, x) + v*sp.diff(u, y) - nu*sp.diff(u, y, 2)   # residual of (9.18)
  print(sp.simplify(res*4*x**2/C**2))                      # -(4 f''' + f f'' + 2 f'^2) at eta: the coefficient of f''' is 4
  ```

### D21 · Wall-jet integration: $4ff''-2f'^2+f^2f'=0$, $f^{-1/2}f'+\frac{f^{3/2}}6=\frac{f_\infty^{3/2}}6$, the implicit solution (9.83) and $f''(0)=f_\infty^3/72$ — ★★★, 15 steps, in C13 (notebook · `wall_jet_invariant`)
- **Goal.** Integrate the wall-jet ODE in closed form, including the constants the book leaves out.
- **Start.** $4f'''+ff''+2f'^2=0$ (D20) with $f(0)=f'(0)=0$, $f'(\infty)=0$ — *in words:* three integrations' worth of structure, if we find the right integrating factors.
- **Plan.** (1) Multiply by f and spot total derivatives. (2) Divide by f^{3/2} and integrate again. (3) Separate variables with g² = f/f_∞ and use partial fractions. (4) Read off f″(0) and the free scale.
- **Tools.** the total-derivative move (primer P216) · partial fractions, `apart` and implicit inversion (primer P219) · exponent rules (ch01 P43) · scaling symmetry (primer P207) · sympy (ch01 P40).
- **Assumptions.** f > 0 for η > 0; f → f_∞ (finite) at infinity.
- **Steps.**
  1. *did:* Restate the problem · *tex:* $4f'''+ff''+2f'^2=0,\quad f(0)=f'(0)=0,\ f'(\infty)=0$ · *why:* The result of D20 (coefficient 4). · *plain:* what we integrate.
  2. *did:* Multiply by f · *tex:* $4ff'''+f^2f''+2ff'^2=0$ · *why:* An integrating factor: we hope the result is a sum of total derivatives. · *plain:* every term has one more f.
  3. *did:* Spot the derivatives · *tex:* $4ff'''=4(ff'')'-2(f'^2)',\qquad f^2f''+2ff'^2=(f^2f')'$ · *why:* Product rule backwards (P216): (ff″)′ = f′f″ + ff‴ and (f′²)′ = 2f′f″, so 4ff‴ = 4(ff″)′ − 4f′f″ = 4(ff″)′ − 2(f′²)′; and (f²f′)′ = 2ff′² + f²f″. · *plain:* everything is a derivative.
  4. *did:* Integrate once · *tex:* $4ff''-2f'^2+f^2f'=0$ · *why:* The constant is zero because f = f′ = 0 at η = 0 (the book evaluates it at η = 0). · *plain:* a second-order equation.
  5. *did:* Divide by $4f^{3/2}$ · *tex:* $f^{-1/2}f''-\tfrac12f^{-3/2}f'^2+\tfrac14f^{1/2}f'=0$ · *why:* A second integrating factor; f > 0 for η > 0. · *plain:* again look for derivatives.
  6. *did:* Recognise them · *tex:* $\Big(f^{-1/2}f'+\tfrac16f^{3/2}\Big)'=0$ · *why:* (f^{−1/2}f′)′ = f^{−1/2}f″ − ½f^{−3/2}f′², and (f^{3/2}/6)′ = ¼f^{1/2}f′. · *plain:* another total derivative.
  7. *did:* Integrate · *tex:* $f^{-1/2}f'+\tfrac16f^{3/2}=\tfrac16f_\infty^{3/2}$ · *why:* As η → ∞, f′ → 0 (faster than f^{1/2}) and f → f_∞, which fixes the constant. · *plain:* a first-order equation.
  8. *did:* Solve for f′ and separate variables · *tex:* $\frac{df}{f_\infty^{3/2}f^{1/2}-f^2}=\frac{d\eta}6$ · *why:* f′ = f^{1/2}(f_∞^{3/2} − f^{3/2})/6 and f^{1/2}f^{3/2} = f². The book prints f in place of f^{1/2} in the first term of the denominator. · *plain:* an integrable equation.
  9. *did:* Substitute g² = f/f_∞ · *tex:* $f=f_\infty g^2\ \Rightarrow\ \frac{dg}{1-g^3}=\frac{f_\infty}{12}\,d\eta$ · *why:* df = 2f_∞g dg, f^{1/2} = f_∞^{1/2}g, f^{3/2} = f_∞^{3/2}g³ turn the left side into 2dg/(f_∞(1 − g³)). · *plain:* g runs from 0 to 1.
  10. *did:* Split into partial fractions · *tex:* $\frac1{1-g^3}=\frac13\Big[\frac1{1-g}+\frac{g+2}{1+g+g^2}\Big]$ · *why:* 1 − g³ = (1 − g)(1 + g + g²); putting the right side over the common denominator gives 1 (sympy `apart`, P219). · *plain:* two simple integrals.
  11. *did:* Integrate · *tex:* $\tfrac13\Big[-\ln(1-g)+\tfrac12\ln(1+g+g^2)+\sqrt3\tan^{-1}\tfrac{2g+1}{\sqrt3}\Big]=\tfrac{f_\infty}{12}\eta+c$ · *why:* ∫dg/(1 − g) = −ln(1 − g) and (g + 2)/(1 + g + g²) = ½(2g + 1)/(1 + g + g²) + (3/2)/(1 + g + g²), whose integrals are ½ln(1 + g + g²) and √3 arctan((2g + 1)/√3). · *plain:* logs and an arctangent.
  12. *did:* Fix the constant · *tex:* $-\ln(1-g)+\sqrt3\tan^{-1}\frac{2g+1}{\sqrt3}+\ln(1+g+g^2)^{1/2}=\frac{f_\infty}4\eta+\sqrt3\tan^{-1}\frac1{\sqrt3}$ (9.83) · *why:* Multiply by 3; g(0) = 0 sets c so that the right side carries √3 arctan(1/√3) = √3π/6. To get g(η) invert numerically with `brentq` (P219). · *plain:* left-implicit closed form.
  13. *did:* Read the far field · *tex:* $1-g\approx\sqrt3\,e^{\sqrt3\pi/6}\,e^{-f_\infty\eta/4}=4.29\,e^{-f_\infty\eta/4}$ · *why:* As g → 1 the −ln(1 − g) term dominates; solve for 1 − g with the constants of step 12 (√3π/6 − √3π/3 − ½ln3). The book omits the factor 4.29. · *plain:* exponential decay with rate f_∞/4.
  14. *did:* Find f″(0) · *tex:* $\sqrt{2f''(0)}=\frac{f_\infty^{3/2}}6\ \Rightarrow\ f''(0)=\frac{f_\infty^3}{72}$ · *why:* Step 7 at η → 0: f ≈ ½f″(0)η² and f′ ≈ f″(0)η, so f^{−1/2}f′ → √(2f″(0)) while the f^{3/2} term vanishes. · *plain:* f_∞ is tied to the initial curvature.
  15. *did:* Identify the free scale · *tex:* $f\to\lambda f(\lambda\eta):\ \ f_\infty\to\lambda f_\infty,\ \ C\to C/\lambda^2$ · *why:* Every term of 4f‴ + ff″ + 2f′² scales like λ⁴ (D06 step 1) and the conditions are homogeneous, so f_∞ is not fixed by the problem; the physical constant is C f_∞². · *plain:* f_∞ is a gauge; only C f_∞² matters.
- **Result.** f^{−1/2}f′ + f^{3/2}/6 = f_∞^{3/2}/6, the implicit solution (9.83), f″(0) = f_∞³/72, f′ ≈ 2.14 f_∞² e^{−f_∞η/4} far out, and f_∞ a free scale — *in words:* a closed-form wall-jet profile with one gauge freedom.
- **Check.** f_∞ = 1: an IVP from f″(0) = 1/72 ends at f_∞ = 1.000000, satisfies (9.83) to 1e-11 at η = 2, 5, 10, 20; ∫f′ = 1, ∫f′² = 1/18, ∫f′∫f′² = 1/40; peak f′ = 0.0787 at η = 8.11; 4.29 confirmed numerically (4.2896).
- **What it means.** Together with D19–D20: u₀ ∝ x^{−1/2}, δ ∝ x^{3/4}, ṁ ∝ x^{1/4}; one physical parameter (Ψ ~ ṁ⁴/x).
- **Traps.** f^{1/2}, not f, in the separated integral; arctan branch and the constant √3π/6; g(0) = 0 fixes c; the free scale means f_∞ and C are not both determined by data.
- **sympy check intent (★★★).** Verify by differentiation that the two first integrals (steps 4 and 7) satisfy the ODE, that the partial-fraction split (step 10) is exact, and that the implicit relation (9.83) reproduces the IVP solution. `check_src` sketch:
  ```python
  import sympy as sp                                       # symbolic algebra
  g = sp.symbols('g', positive=True)                       # g = sqrt(f / f_inf), 0 < g < 1
  print(sp.simplify(sp.apart(1/(1 - g**3), g) - sp.Rational(1, 3)*(1/(1 - g) + (g + 2)/(1 + g + g**2))))   # 0: step 10
  eta_of_g = 4*(-sp.log(1 - g) + sp.sqrt(3)*sp.atan((2*g + 1)/sp.sqrt(3)) + sp.log(1 + g + g**2)/2 - sp.sqrt(3)*sp.atan(1/sp.sqrt(3)))   # (9.83) solved for eta (f_inf = 1)
  print(sp.simplify(sp.diff(eta_of_g, g) - 12/(1 - g**3)))  # 0: d(eta)/dg = 12/(f_inf (1 - g^3)) from step 9
  f = sp.Function('f'); e = sp.symbols('e')                # check the first integral against the ODE
  first = 4*f(e)*f(e).diff(e, 2) - 2*f(e).diff(e)**2 + f(e)**2*f(e).diff(e)
  print(sp.simplify(first.diff(e) - f(e)*(4*f(e).diff(e, 3) + f(e)*f(e).diff(e, 2) + 2*f(e).diff(e)**2)))   # 0: step 4 differentiated equals f times the ODE
  ```

### D22 · The teacup: net inward force per volume $\rho(u_e^2-u^2)/R$ in the slowed bottom layer — ★★, 7 steps, in C14 (notebook · `teacup_secondary_flow`)
- **Goal.** Show why a friction layer beneath a swirling fluid drives a flow toward the axis.
- **Start.** Circular motion of the core at speed u_e(R) with radial balance ∂p/∂R = ρu_e²/R (Ch. 4) — *in words:* the book states the argument in words; we write the force balance.
- **Plan.** (1) The pressure gradient set by the core. (2) It is the same inside a thin layer. (3) Compare it with what circular motion needs there.
- **Tools.** radial force balance in a swirl (primer P220) · the thin-layer argument ∂p/∂z ≈ 0 (as (9.10)).
- **Assumptions.** Thin layer; steady, axisymmetric swirl; u_φ(z) illustrative (no viscous solution here — Ch. 13).
- **Steps.**
  1. *did:* Balance in the inviscid core · *tex:* $\frac{\partial p}{\partial R}=\frac{\rho u_e^2}{R}$ · *why:* A parcel on a circle of radius R needs the inward force ρu_e²/R (Ch. 4 centripetal acceleration); in the core only the pressure gradient supplies it. · *plain:* pressure holds the fast water on its circle.
  2. *did:* Carry the pressure into the layer · *tex:* $\frac{\partial p}{\partial z}\approx0\ \Rightarrow\ p(R,z)=p_e(R)$ · *why:* Thin-layer argument as (9.10): the pressure does not change across the layer, so the bottom layer feels the core's ∂p/∂R. · *plain:* the layer inherits the pressure field.
  3. *did:* Write the net force per volume · *tex:* $F_{in}=\frac{\partial p}{\partial R}-\frac{\rho u^2}R=\frac{\rho\,(u_e^2-u^2)}R$ · *why:* Force pushing inward minus the inward force a parcel of speed u needs (its centripetal requirement ρu²/R). Positive = net inward. · *plain:* the surplus of pressure over need.
  4. *did:* Find the sign · *tex:* $u<u_e\ \Rightarrow\ F_{in}>0$ · *why:* Friction slows the water near the floor (no slip, u → 0); u_e² − u² > 0 there. · *plain:* inward on the floor, zero in the core.
  5. *did:* Put numbers · *tex:* $F_{in}(0)=\frac{\rho u_e^2}R=\frac{1000\times0.04}{0.04}=1000\ \mathrm{N/m^3}\approx0.10\,\rho g$ · *why:* Water, u_e = 0.2 m/s, R = 4 cm; the weight density is ρg = 9810 N/m³. · *plain:* a small but unopposed push.
  6. *did:* Close the loop with continuity · *tex:* $\text{inflow on the floor}\ \Rightarrow\ \text{up near the axis}\to\text{out at the top}\to\text{down the side wall}$ · *why:* Mass conservation (6.2): fluid moved inward must go somewhere; the side wall and free surface complete the meridional loop; leaves ride the floor inflow. · *plain:* the leaves collect in the middle.
  7. *did:* State the limit and the hook · *tex:* $u=u_e\ \Rightarrow\ F_{in}=0;\quad\delta_E\sim\sqrt{\nu/\Omega},\ \Omega=u_e/R$ · *why:* In the core pressure and centripetal need balance exactly; the layer depth follows from viscous diffusion over one rotation time (Ekman scale, Ch. 13: 0.45 mm here). · *plain:* the same mechanism drives Ekman layers and spin-up.
- **Result.** $F_{in}=\rho(u_e^2-u^2)/R>0$ inside the friction layer, zero outside — *in words:* a slowed layer under an unchanged pressure gradient is pushed toward the axis.
- **Check.** Units: (kg/m³)(m²/s²)/m = N/m³ ✓. u_e = 0.2, R = 0.04, ρ = 1000: F at u = 0, 0.5u_e, 0.75u_e, u_e = 1000, 750, 437.5, 0 N/m³ ✓ (`secondary_flow_radial_force`).
- **What it means.** Secondary flows arise whenever a boundary layer sits under curved streamlines (river bends, Ekman layers, cyclone spin-up).
- **Traps.** The inward flow is *not* caused by centrifugal force alone; the pressure is set by the fast core; sign of F (inward positive); the layer profile used is illustrative, not a solution.

