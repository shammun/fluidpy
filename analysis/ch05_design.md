# Chapter 5 — Vorticity Dynamics: lesson design
(from `analysis/ch05_curation.md` (A 14 · B 69 · C 11 = 94 rows; CORE 14 · NOTE 60 · RECAP 18 · SKIP 2; 23 derivations
D01–D23 written out (★ 5 · ★★ 14 · ★★★ 4); E1–E9 + backup B1) and `analysis/ch05.md` (§2b derivations a-D01…a-D31, §4
implementation rows 1–39, §5 core modules, §9 slips and conventions); every equation placed below was re-read on the
rendered pages `chapters/pages/ch05/p198–p218.png` (printed pp. 171–191: (5.1)–(5.3) p199, (5.4), (5.5a, b) and
Fig. 5.1 p200, (5.6), the isobar line, −½u_θ² + gz + p/ρ, σ_rθ and Fig. 5.2 p201, (5.7), the funnel, u_θ = ωa²/2r and
Fig. 5.3 p202, Γ = πa²ω, (5.8)–(5.10) p203, Fig. 5.4, D(dx)/Dt = du, ∮u_i du_i = 0, (5.11) p204, Fig. 5.5 and the four
restrictions p205, Fig. 5.6, restriction (4) and Helmholtz's theorems p206, Fig. 5.7, the proof and (5.12) p207, the
Lamb step, (B.3.10), (5.13), ∇×ω = −∇²u and (5.14) p208, the integrand rewrite, (5.15)–(5.17) p209, Fig. 5.8, (5.18),
(5.19, 5.20) p210, (5.21)–(5.27) p211, (5.28)–(5.30) p212, Fig. 5.9, (5.31), (5.32) and the planetary component
equations p213, Fig. 5.10 and (5.33) p214, Figs. 5.11–5.12 and V₁, V₂ p215, Fig. 5.13 and Γ/(4πh) p216, Figs. 5.14–5.15
p217, leap-frogging, §5.8, dΓ = (u₂ − u₁)ds and Fig. 5.16 p218); 2026-09-23, lesson-designer. The implementer works in
parallel from analysis §4 + curation §8; **Part C is written first and is the contract both sides keep.**)

**Binding conventions for every builder (analysis §9, curation decisions 4–6, ch04 conventions carried over).**
1. **Fields are callables** exactly as in Ch. 3–4: `u(x, t) -> ndarray` with `x` of shape `(d,)` or `(d, N)` (components
   on axis 0, d = 2 or 3); scalar fields `F(x, t)`; a density `rho` may be a float **or** a callable. Every public
   function is scalar-callable (scalar in → float out). `py:` parity expressions have **no builtins and no lambdas**:
   functions that take a field also accept a **preset name** plus keywords, lists stand for vectors, and results are
   indexed down to one number (`["flux"]`, `[1]`, `.total`).
2. **ω means the vorticity in every ch05 function** — never a rotation rate. The solid-body rotation of (5.1),
   $u_\theta=\omega r/2$, turns at **ω/2**; ch03's `core.vortices.solid_body_rotation(r, omega0)` takes the rotation rate
   ω₀ = ω/2, and `ch05.solid_body_from_vorticity` converts. A tank turning at `Omega_tank` [rad/s] has ω = 2Ω_tank. ⚠️
   Fig. 5.2 labels the tank's rotation "2ω" — a slip: with (5.1) the tank turns at ω/2 (taught at R02, D02).
3. **Axes, gravity, orientation.** z **up**, $\mathbf g=-g\mathbf e_z$ (2-D sketches: y up, $\mathbf g=-g\mathbf e_y$),
   Φ = gz, **g = 9.81 m/s²** default. Plane vorticity, circulation and point-vortex strengths are **counterclockwise
   positive** about +z; a loop's circulation uses the counterclockwise (right-hand, n = +e_z) orientation unless the
   function says otherwise. Water ρ = 1000 kg/m³, μ = 1.0×10⁻³ Pa s, ν = 1.0×10⁻⁶ m²/s in worked numbers; air ρ = 1.2.
4. **Five meanings of Γ and two of γ** (⚠️ callout in the notation cell and again at C14): `Gamma` [m²/s] is **always a
   circulation** (of a loop, a vortex, a filament, a tube); `gamma` [m/s] is the **vortex-sheet strength** (circulation per
   unit length; the book prints Γ ≡ dΓ/ds in §5.8, the exercises write γ). Earlier meanings (lapse rate, shear-rate
   preset) never appear in ch05 code. γ is not C_p/C_v here.
5. **σ, ε, Π, the prime.** σ_ij = viscous stress (C02, C03); the Gaussian core radius of Ch. 3 is `sigma` only inside
   `core.vortices` calls and is written σ_c in prose (⚠️ at N06). ε_ijk = Levi-Civita (C09); the dissipation rate is
   `eps_diss`, written ε only in N08 with a ⚠️. The primes of §5.5 (x′, ∇′, V′, A′) mark the **source** point and
   integration variables — not a moving frame. (5.26)'s "Π" is a slip for Φ (N33).
6. **Two derivations of the vorticity equation.** (5.13) (§5.4: ρ constant, inertial frame, vector route, D09) and (5.30)
   (§5.6: Boussinesq, frame rotating at constant Ω, index route, D15). (5.30) with Ω = 0 and ∇ρ ∥ ∇p must give (5.13):
   asserted in D15 step 15, in `vorticity_budget_sym` and in E7. In (5.30) u and ω are **relative** (measured in the
   rotating frame); the stretching term acts on the **absolute** vorticity ω + 2Ω; the pressure term keeps the full
   1/ρ² although the continuity equation is Boussinesq (⚠️ D15 trap: to leading order ∇ρ×∇p/ρ² ≈ ∇ρ′×g/ρ₀ with
   g = −g e_z, cf. (4.86)).
7. **Book slips taught in corrected form** (each gets "the book prints X; the correct form is Y" where it is used):
   (5.14)'s −1/(4π) → **+1/(4π)** (N25, D10 steps 10–11, E6 sign toggle) and the compensating slip in the integrand
   rewrite, printed "+((x − x′)/|x − x′|³) × ω", correct "−((x − x′)/|x − x′|³) × ω = +ω × (x − x′)/|x − x′|³" (N26,
   D11 steps 3–5) · "(see Exercise 5.8)" after (5.14) → Exercise 5.9 (N25) · "(see Exercise 5.11)" after (5.33) →
   Exercise 5.10 (C11, D19) · Fig. 5.2's "2ω" → the tank turns at ω/2 (R02, D02 step 1) · the Kelvin proof's reason
   "ρ and p are single valued" → the integral vanishes because ρ = ρ(p) makes dp/ρ exact (R09, D05 step 6) · "C lies
   in irrotational fluid ⇒ no net viscous force" holds for incompressible constant-μ flow only (N13, D05 step 9) · the
   Lamb step's "∇(u·u)" → ∇(½u·u), harmless (R11, D09 step 5) · (5.26)'s Π → Φ (N33, D15 step 3) · (5.27) silently drops
   u_{j,j}(ω_n + 2Ω_n) (N34, D15 step 8) · Fig. 5.16's caption dΓ/ds = u₁ − u₂ vs the text's u₂ − u₁ (C14, D23 step 5) ·
   "hyperboloids of revolution of the second degree" → the funnel isobars (c − z)r² = const are cubic surfaces (N07) ·
   **`scipy.special.ellipk(m)` takes the parameter m = k², not the modulus k** (a code trap, P148, N42) · (5.10) writes
   ∂σ_ij/∂x_j where Ch. 4's Cauchy contracts the first index — harmless, σ is symmetric (D05 step 1) · **new
   (designer, p215):** Fig. 5.11's caption calls G "the point where the combined velocity induced by the two vortices
   is zero" — true only for Γ₁ = Γ₂; G is the centre of vorticity, the fixed centre of the pair's rotation (N40, D21
   check, E8 Explain §3).
8. **Colours (text, figures, explainers — one meaning each; curation §5):** relative vorticity `teal` · planetary
   vorticity 2Ω `amber` · baroclinic `orange` · diffusion `rose` · stretching `accent` purple · tilting `blue` ·
   advective (u·∇)ω `muted` grey · local ∂ω/∂t `blue` dashed (bars) · residual black · pressure `orange`, density `blue`
   (E4, isobars orange, isopycnals blue) · **positive (counterclockwise) vortices `teal`, negative `rose`** · centripetal
   acceleration `teal`, pressure-gradient force `orange`, viscous force `rose` (E2) · the tube's lower end `blue`, side
   `muted`, upper end `orange` (E1, D01) · ghosts and references `muted` grey dashed.
9. **Every book equation is shown in full next to its number** (CLAUDE.md rule 3; memory rule) — in markdown, derivation
   steps ("substitute (5.10), $\oint_C\frac{Du_i}{Dt}dx_i=-\oint_C\frac1\rho dp-\dots$"), traps, recaps, and explainer tour,
   Explain, Derivation, notes, status and quiz text. Builders reuse `show_eqs(text, EQ)` (from `notebooks/build_ch04.py`)
   with an `EQ` dict covering all 34 labels of ch05 plus the Ch. 3–4 labels cited here ((3.5), (3.7), (3.18), (3.22),
   (3.25), (3.28), (3.29), (4.18), (4.19), (4.24), (4.37), (4.39b), (4.40), (4.41), (4.45), (4.59), (4.67), (4.68),
   (4.69), (4.86), (2.19), (2.30), (2.34)); JS explainers use a local `showEqs`. `tools/eq_refs.py` must list 0
   offenders. (B.3.10), (B.3.13) are Appendix-B identity labels: always written out too.
10. **nbkit behaviour** (ch04 lesson): `nb.recap(...)` and `nb.section(...)` end the current CORE block and `nb.core(id)`
    is called once — so every RECAP sits **before** the `nb.core` call of the block that uses it (R01 before C01; R02–R08
    before C02; R09 before C03; R10 R11 before C06; R12–R17 before C09; R18 before C11); every `nb.derivation` sits inside
    its curation CORE block; `nb.derivation(ref=…)` gets an equation number only ("5.4"), never "Exercise 5.9". Items
    taught in another section's block get a one-line `nb.pointer` in their own section (N21 Burgers in §5.4 → C10;
    C05's field-equation check in §5.3 → end of C06).
11. **Parsing.** Part E is the only part with table rows after its heading; Part F has no line starting with `|`.
    Explainer headings are exactly `### E1 · vortex_tubes_cannot_end` … `### E9 · vortex_sheet_rollup` and
    `### B1 · vortex_rings`. Primer terms in Part A are the exact Concept text of the Part E rows marked "primer"
    (coverage_check matches the first 18 characters). "Explained by" never names an earlier chapter's C-number (cite
    "Ch. 3 §3.5" instead).
12. **Public repo.** Exercise text is never reproduced; the numbers of Exercises 5.1, 5.2 and the printed closed forms of
    Exercises 5.5, 5.6, 5.11–5.13, 5.19, 5.20 stay in `tests/book_values_ch05.json`. The notebook shows **our**
    derivations (D03, D07, D10, D18–D21) and **our** numbers.

Order of parts: C (contract) · A (notebook storyboard) · B (explainers) · D (runtime) · E (prerequisite ledger) · F
(derivations).

---

## Part C — functions the builders will call (the implementer's contract)

**Status column.** **§4 #n** = planned in `analysis/ch05.md` §4 row n (signature kept or refined here). **§8** = added by
the curation's notes for the implementer (`analysis/ch05_curation.md` §8). **NEW** = added by this design (in neither) —
flagged as the phase asks. **reuse** = exists (ch01–ch04) and is only called. Every callable is reachable as
**`ch05.<name>`**: the chapter module `fluidpy/ch05_vorticity_dynamics.py` re-exports every public name of
`core/vorticity.py`, `core/biot_savart.py` and the ch05 additions to `core/vortices.py` and `core/integral_theorems.py`
(so the notebook imports one module and parity rows write `ch05.…`). Notebook alias: `from fluidpy import
ch05_vorticity_dynamics as ch05`. Units SI: x [m], t [s], u [m/s], ω [1/s], Γ [m²/s], γ [m/s], p [Pa], ρ [kg/m³],
μ [Pa s], ν [m²/s], Ω [rad/s], G [1/s]. Docstrings cite § and Eq. (from the rendered page), list symbols with units and
assumptions, and carry the validation label. Return types: `NamedTuple` for term splits (unpackable, attribute access),
`dict` for scenario/geometry results, `@dataclass(frozen=True)` for traced objects (`Tube`).

### C.0 Reused (existing; called, not changed)
| # | Callable (signature) | Returns | Used by | Status |
|---|---|---|---|---|
| 0.1 | `style.setup_notebook() -> bool` (FAST) · `style.COLORS` · `style.savefig(fig, "ch05", name)` | FAST · palette · path | setup, every figure | reuse |
| 0.2 | `anim.animate(update, frames, fig, interval)` · `anim.show_animation(anim, player="video"\|"frames")` | HTML | C03, C10, C12, C13, C14 | reuse |
| 0.3 | `interact.slider_figure(fn, name, values, *, unit, xlabel, ylabel, title, xrange, yrange, height)` · `interact.animate_figure` | plotly Figure | C02, C04, C06, C07, C08, C11, C14 | reuse |
| 0.4 | `embed.show_viz("ch05", slug)` | display | 9 explainer cells | reuse |
| 0.5 | `tools.convergence.observed_order(h, err) -> float` | slope | C04 (O(R²)), C07 (quadrature), C08 (polygon ring), C14 (1/N) | reuse |
| 0.6 | `core.vortices`: `solid_body_rotation(r, omega0)` (ω₀ = rotation rate!), `line_vortex(r, B)` (Γ = 2πB), `rankine_vortex(r, Gamma, sigma) -> (u_theta, omega_z)`, `gaussian_vortex(r, Gamma, sigma) -> (u_theta, omega_z)`, `vortex_velocity_field(profile, **params) -> u(x, t)`, `circulation_circle(u_theta, r, center, n, **kw)`, `polar_vorticity_z(_sym)` | arrays | R01–R03, R08, C01 (N03), C03, C07 | reuse |
| 0.7 | `core.integral_theorems`: `Loop`, `Surface`, `planar_loop(center, normal, radius, n)`, `rectangle_loop`, `planar_disc(center, normal, radius, nr, ntheta)`, `circulation(u_fn, loop)`, `curl_flux(u_fn, surface, curl_fn, fd_step)`, `stokes_theorem_check`, `divergence_theorem_box(Q_fn, bounds, n, div_fn, rule)`, `gauss_legendre_nodes(lo, hi, n)`, `fd_partials` | floats | C01 (N03), C03, C07 (N27) | reuse |
| 0.8 | `core.kinematics`: `streamline(u, x0, t_frozen, s_max, both, n)`, `pathline(u, r0, t0, t_eval, rtol, atol)`, `vorticity(u, x, t, h)`, `velocity_gradient_at(u, x, t, h)`, `vorticity_from_gradient(G)`, `vorticity_in_rotating_frame(omega, Omega)`, `acceleration(u, x, t, h, ht)`, `material_derivative_terms`, `linear_flow_map(G, t)` | arrays | C01, C03, C05, C06, C09, C10, R16, R18 | reuse |
| 0.9 | `core.navier_stokes`: `exact_solution(name, x, t, **p)` / `exact_field(name, **p)` (`"lamb_oseen"`, `"taylor_green"`, `"solid_body"`), `ns_incompressible_terms`, `viscous_force_forms(u, x, t, mu, h) -> (lap, div2S, minus_curl_omega)`, `lamb_vector`, `lamb_identity_terms(u, x, t, h)`, `lamb_identity_sym(u_exprs, coords)`, `navier_stokes_sym`, `continuity_residual` | arrays, sympy | R05, R11, R13, R15, R17, C02 (N09), C03 (N13), C06 | reuse |
| 0.10 | `core.constitutive`: `newtonian_stress(G, p, mu, …)`, `viscous_stress(G, mu, mu_v)`, `dissipation_rate(G, rho, mu, mu_v, form)` · `core.tensors`: `levi_civita()`, `strain_rate_tensor(G)`, `rotation_tensor(G)` | arrays | R05, N06, N08, R16 | reuse |
| 0.11 | `core.curvilinear`: `coordinates("cylindrical")`, `curl(A, system)`, `strain_rate(u, system)`, `vector_laplacian(A, system)`, `advective_acceleration(u, system)`, `laplacian(f, system)` (sympy) | sympy | C02 (D02 check, D03), C10 (D18 check) | reuse |
| 0.12 | `core.bernoulli.bernoulli_function(speed, p, z, rho, g, kind, p_o)` · `core.rotating`: `OMEGA_EARTH`, `coriolis_parameter(lat_rad, Omega)`, `frame_acceleration_terms`, `apparent_body_forces`, `coriolis_acceleration(Omega, u)` | floats | R07, C09, C11 | reuse |
| 0.13 | `core.diffusion.ftcs_diffusion_1d(f0, D, dy, dt, nsteps, …)` · `core.statics.hydrostatic_pressure_uniform(z, p0, rho, g)` · `ch01.fluid_properties(name, T)` | arrays | C06 (N22 check), R06, E2 fluids | reuse |

### C.1 `fluidpy/core/vorticity.py` (NEW module, `core.VD`) — vortex kinematics, circulation theorems, vorticity-equation term splits
Stencil functions take `h` (space, default **1e-3·L** where the docstring states the field's length scale; nested
stencils — ω from u, then ∇²ω — need h ≈ 1e-3 L, not 1e-6) and `ht` (time, default `core.kinematics.DEFAULT_HT`);
second-order central differences through `core._stencil`. Loops are point arrays `pts` of shape `(d, N)` (d = 2 or 3),
**counterclockwise, not closed** (the last point is not repeated), uniformly parametrised by a label s_k = k/N, so the
periodic trapezoid rule and spectral differentiation apply.
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 1.1 | `vorticity_field(u, h=1e-4) -> Callable` (ω(x, t) by stencils; wraps `core.kinematics.vorticity`; 2-D fields give ω of shape (1,) or (1, N)) | R01, ω = ∇×u | 1/s | C01, C05, C06 | §4 #2 |
| 1.2 | `vortex_line(omega, x0, t=0.0, s_max=1.0, n=400, both=True) -> (s, X)` (arc-length ODE dx/ds = ω/\|ω\| via `core.kinematics.streamline(omega, x0, t_frozen=t, …)`; `X` shape (3, n); `omega` a callable or a preset name of 3.4/3.6 with its keywords) | (5.3), N02 | m | C01 (N02), E1 | §4 #2 |
| 1.3 | `TubeStrength(circulation, flux, diff)` · `vortex_tube_strength(u, loop, surface, omega=None, t=0.0) -> TubeStrength` (`circulation` = ∮u·dx by `integral_theorems.circulation`; `flux` = ∫ω·n dA by `curl_flux` (exact `omega` if given); both with the loop's orientation (n_c into A, ch02)) | N03, Stokes (2.34) | m²/s | C01 (N03) | §4 #3 (signature tightened) |
| 1.4 | `Tube` frozen dataclass (`lines`: list of (3, n) arrays, `seed_center`, `seed_normal`, `seed_radius`, `s`: arc lengths) · `vortex_tube(omega, seed_center, seed_normal, seed_radius, n_lines=24, s_max=1.0, n=200, t=0.0, **p) -> Tube` (vortex lines from n_lines points on the seed circle, both directions) | Fig. 5.1, N45 | m | C01 figure (N45 plotly 3-D) | §4 #4 |
| 1.5 | `TubeFlux(lower, side, upper, total)` · `tube_flux_budget(field="gaussian_tube", z_a=0.0, z_b=1.0, R0=0.1, n=48, R_of_z=None, **p) -> TubeFlux` — for an **axisymmetric** tube whose wall is the surface of revolution R = R_t(z) through (R0, z = 0): `lower` = −∫ω_z dA over the disc at z_a (outward normal −e_z), `upper` = +∫ω_z dA at z_b, `side` = ∫ω·n dA over the wall (n outward, from R_t′(z)), `total` = sum (= ∫∇·ω dV by Gauss); `field` = `"gaussian_tube"` (3.4, R_t from the flux function), `"lamb_oseen"` (straight tube), `"broken"` (3.6, wall R_t = R0, straight) or a callable ω with `R_of_z` required; Gauss–Legendre in R and z | (5.4), D01 | m²/s | C01, E1 parity | §4 #4 (axisymmetric form fixed here) |
| 1.6 | `material_loop(u, pts0, t_eval, rtol=1e-10, atol=1e-12) -> ndarray (T, d, N)` (all N points advected in **one** `solve_ivp` call, state flattened `y.reshape(d, N)`; DOP853) | C03, N48 | m | C03, C05, E3 | §4 #11 |
| 1.7 | `loop_circulation(u, pts, t=0.0) -> float` (∮u·dx by the periodic trapezoid rule with dx/ds from FFT differentiation of the periodic parametrisation: Γ = Σ_k u(x_k)·x′(s_k)/N) · `loop_length(pts) -> float` · `material_circulation(u, pts0, t_eval, rtol=1e-10) -> ndarray (T,)` | (3.18), (5.8) | m²/s, m | C03, C11, E3 | §4 #11 (+ `loop_circulation`, `loop_length` **NEW**: the from-scratch check and the "×6.7 longer" readout need them) |
| 1.8 | `KelvinRate(acceleration, contour, total)` · `kelvin_rate_terms(u, pts, t=0.0, h=1e-4, ht=None) -> KelvinRate` (`acceleration` = ∮(Du/Dt)·dx with `core.kinematics.acceleration`; `contour` = ∮u·du with du = (dx·∇)u along the loop; `total` = sum = DΓ/Dt) | (5.9), N10, N11 | m²/s² | C03 (D04 check), E3 | §4 #12 |
| 1.9 | `KelvinForces(pressure, body, viscous)` · `kelvin_force_terms(pts, grad_p, rho, grad_Phi=None, visc_force=None, t=0.0) -> KelvinForces` (`pressure` = −∮(1/ρ)∇p·dx, `body` = −∮∇Φ·dx, `viscous` = ∮f_v·dx with `visc_force(x, t)` the viscous force **per unit mass** (1/ρ)∂σ_ij/∂x_j [m/s²]; callables of (x, t); `rho` float or callable) | (5.10), (5.11), N12, N14 | m²/s² | C03 (D05 check), C04 (N15), E3 | §4 #12 |
| 1.10 | `kelvin_scenario(name, **p) -> dict(u, p, rho, grad_p, visc_force, nu, Omega, pts0, t_end, label, broken)` — `"cellular"` (ψ = U sin(x/ℓ) sin(y/ℓ) with U = 1 m/s, ℓ = 1 m: u = sin x cos y, v = −cos x sin y; loop = circle r = 0.5 m about (π/2, 0.9); t_end = 6 turnover times 6·2π ≈ 37.7 s; the loop's length grows 3.14 → 10.9 m at 18.8 s → 21.0 m at 37.7 s (×6.7), Γ constant to 1e-12 with 2048 points — point insertion keeps the spacing < 0.05 m), `"rankine_straddle"` (Γ = 2π, σ_c = 1 m; circle r = 0.5 about (0.8, 0)), `"lamb_oseen"` (Γ₀ = 0.01, ν = 1e-6, t₀ = 10 s: the clock starts at t = 0 with the vortex already 10 s old; circle r = 5 mm about 0; material because u_r = 0), `"baroclinic"` (fluid at rest at t = 0 in the lock-exchange field of 5.9 with ρ₁ = 1000, ρ₂ = 1025, δ = 0.1 m; square loop of side 0.2 m centred on the interface; Γ(t) ≈ t·dΓ/dt(0)), `"rotating"` (frame Ω = 0.5 rad/s, axisymmetric convergence α = 0.2 s⁻¹: u = −(α/2)x + (ζ(t)/2)e_z×x with ζ(t) = 2Ω(e^{αt} − 1), loop = unit circle; Γ(t) = 2ΩA₀(1 − e^{−αt}), Γ_a = 2ΩA₀), `"helmholtz"` (as the WIP defines it: an inviscid stretched Gaussian vortex tube, 3-D, Γ = 1, σ₀ = 1, α = 0.5, with a circle of radius 0.3 drawn on the tube wall R = 1 — its flux is 0 and stays 0; E3's Helmholtz picture), `"helmholtz_abc"` (**NEW**: the ABC flow of 5.16, 3-D, a circle of radius 0.01 m about the origin in the plane spanned by ω(0)/\|ω(0)\| and (1, −1, 0)/√2, i.e. lying on the local tube wall of a flow whose tubes are not axisymmetric; t_end = 2 s; Γ = 8.016e-10 m²/s constant, cos(∠(A_vec, ω)) ≤ 2e-5 — C05's quantitative demo); `broken` names the violated restriction (None, "viscous", "baroclinic", "inertial") | E3 flows, C03 animation | mixed SI | C03, E3 | §8 (definitions fixed here) |
| 1.11 | `kelvin_scenario_circulation(name, t, n=256, **p) -> float` (Γ(t) of the scenario's material loop; `"lamb_oseen"`, `"rotating"` closed form, others by 1.6–1.7) · `kelvin_scenario_rate(name, t=0.0, **p) -> dict(acceleration, contour, pressure, body, viscous, coriolis, total)` (the terms of (5.9)/(5.10) at time t, scalar floats) | E3 parity | m²/s, m²/s² | E3 | **NEW** (E3's parity rows need scalar mirrors without lambdas) |
| 1.12 | `pressure_torque_on_element(p="linear", rho="linear", center=(0.0, 0.0), radius=0.01, n=4000, nr=400, grad_rho=(10.0, 0.0), grad_p=(0.0, -9810.0), rho0=1000.0, p0=1.0e5) -> dict(force, x_G, torque, I_G, spin_up, baroclinic, ratio)` — a disc in the (x, y) plane; `"linear"` builds ρ = ρ₀ + ∇ρ·x, p = p₀ + ∇p·x (else callables of (x, y)); `force` = −∮p n ds (2,), `x_G` centre of mass (2,), `torque` = z-torque of the pressure force about G, `I_G` = ∫ρ\|x − x_G\|² dA, `spin_up` = 2·torque/I_G, `baroclinic` = (∇ρ×∇p)_z/ρ(center)², `ratio` = spin_up/baroclinic | D06, Fig. 5.6, N50 | N/m, m, N, kg m²/m, 1/s² | C04, E4 parity | §4 #16 (keyword form fixed here) |
| 1.13 | `baroclinic_term(rho, p, x, t=0.0, h=1e-4) -> (3,)` ((∇ρ×∇p)/ρ² by stencils; 2-D inputs return (0, 0, value)) · `baroclinic_term_sym(rho_expr, p_expr, coords) -> sympy.Matrix` | (5.28), N35 | 1/s² | C04 (D07 check), C09, E4, E7 | §4 #29 |
| 1.14 | `frozen_in_check(u, x0, delta0, t_span, t_eval, h=1e-5, nu=0.0) -> dict(t, x, delta, omega, angle, ratio)` (one `solve_ivp` for the particle x, a material element D(δx)/Dt = G δx and the vorticity Dω/Dt = G ω (+ ν∇²ω from stencils when nu > 0) along the path; `angle` between δx and ω [rad]; `ratio` = \|ω\|/\|δx\| normalised to 1 at t = 0) | N17, Helmholtz 1 | –, rad | C06 end (C05's field-equation route), E5 | §4 #17 |
| 1.15 | `vorticity_equation_sym(u_exprs, coords, t, nu, p_expr=None, Phi_expr=None, rho=None) -> dict(curl_local, curl_advective, curl_pressure, curl_gravity, curl_viscous, lamb_curl, identity_B310, residual_513)` (sympy curls of each term of (4.39b); `curl_pressure` and `curl_gravity` simplify to 0 for constant ρ; `identity_B310` = ∇×(ω×u) − [(u·∇)ω − (ω·∇)u + ω∇·u − u∇·ω] (→ 0 for any fields); `residual_513` = the difference between the curl of NS and (5.13)) | (5.12), (5.13), (B.3.10), D09 | sympy | C06 (D09 check) | §4 #18 (outputs fixed here) |
| 1.16 | `VorticityTerms(local, advective, stretching_tilting, diffusion, residual)` · `vorticity_terms(u, x, t=0.0, nu=0.0, h=1e-3, ht=None) -> VorticityTerms` ((3,) vectors; residual = local + advective − stretching_tilting − diffusion) | (5.13) | 1/s² | C06, C10, E5 | §4 #19 |
| 1.17 | `VorticityBudget(local, advective, stretching_tilting, planetary, baroclinic, diffusion, residual)` · `vorticity_budget(u, x, t=0.0, nu=0.0, Omega=(0.0, 0.0, 0.0), rho=None, p=None, h=1e-3, ht=None) -> VorticityBudget` (relative ω; `stretching_tilting` = (ω·∇)u, `planetary` = (2Ω·∇)u, `baroclinic` = (∇ρ×∇p)/ρ², residual = local + advective − the other four) · `vorticity_budget_sym(u_exprs, coords, t, nu, Omega, rho_expr=None, p_expr=None) -> dict(terms, residual_530, reduces_to_513)` | (5.30), N37 | 1/s² | C09, C11 (N39), E7 | §4 #19 |
| 1.18 | `vorticity_budget_preset(name, x=0.0, y=0.0, z=0.0, t=0.0, component=2, **p) -> dict(local, advective, stretching_tilting, planetary, baroclinic, diffusion, residual)` (scalar floats of one component) — `"burgers"` (Γ = 1e-3, α = 1, ν = 1e-6; probe at R = 1 mm), `"lamb_oseen"` (Γ = 0.01, ν = 1e-6, t₀ = 10 s), `"taylor_green"`, `"rotating_column"` (Ω = 0.5 rad/s, α = 0.2 s⁻¹ axial stretching of a column at rest in the rotating frame at t = 0: planetary = 2Ωα), `"lock_exchange"` (u = 0, t = 0: baroclinic = 2.4222 s⁻² at the interface, all else 0), `"hill"` (A = 1, a = 1, ν = 0) | E5, E7 scenes, N37 | 1/s² | C09 figure, E5, E7 parity | **NEW** (parity-friendly wrapper of 1.16–1.17, the ch04 `ns_terms_preset` pattern) |
| 1.19 | `vorticity_divergence(u, x, t=0.0, h=1e-3) -> float` (∇·(∇×u) by nested stencils; round-off level) | (5.18), N30 | 1/(m s) | C09 (N30) | §4 #26 |
| 1.20 | `rotating_lamb_form_terms(u, p, rho, Phi, x, t=0.0, nu=0.0, Omega=(0.0, 0.0, 0.0), h=1e-4, ht=None) -> dict(local, grad_B, lamb_abs, pressure, viscous_curl, residual)` (the five terms of (5.25): ∂u/∂t, ∇(½u² + Φ), −u×(ω + 2Ω), −(1/ρ)∇p, −ν∇×ω; residual = local + grad_B + lamb_abs − pressure − viscous_curl, equal to the (5.20) residual) | (5.25), N32 | m/s² | C09 (D14 check) | §4 #28 |
| 1.21 | `stretching_tilting_split(omega, G) -> dict(e_s, rate, stretching, tilting, total)` (`e_s` = ω/\|ω\|; `rate` = e_s·G e_s [1/s]; `stretching` = rate·ω (∥ ω); `tilting` = Gω − stretching (⟂ ω); `total` = Gω = (ω·∇)u) | (5.31), (5.32), N38 | 1/s, 1/s² | C10, E5 parity | §4 #31 |
| 1.22 | `planetary_vorticity_terms(G, Omega) -> (3,)` (2(Ω·∇)u = 2GΩ; for Ω = Ωe_z this is 2Ω∂u/∂z) | N39 | 1/s² | C11 (N39) | §4 #33 |
| 1.23 | `loop_vector_area(pts) -> (3,)` (½∮x×dx by the periodic trapezoid rule; a planar counterclockwise loop in the (x, y) plane gives (0, 0, A)) · `absolute_circulation(u, pts, Omega, t=0.0) -> (Gamma, Gamma_a)` (Γ_a = Γ + 2Ω·A_vec) | (5.33), D19 | m², m²/s | C11, E7 | §4 #34 |

### C.2 `fluidpy/core/biot_savart.py` (NEW module, `core.BS`) — velocity from vorticity
Positions are arrays `(3,)`/`(3, N)` (3-D) or `(2,)`/`(2, N)` (planar point vortices, strengths counterclockwise
positive). Every kernel excludes the self term (\|x − x′\| < 1e-14 contributes nothing) and accepts a smoothing length
`eps` (Rosenhead–Moore: \|r\|² → \|r\|² + ε²) for points inside the vorticity.
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 2.1 | `poisson_green_3d(x, xp) -> float \| ndarray` (G = −1/(4π\|x − x′\|): ∇²G = δ) | D10, P138 | 1/m | C07 (D10 check) | §4 #22 |
| 2.2 | `velocity_from_curl_omega(curl_omega, x, nodes, weights, sign=+1.0) -> (3,)` (u = sign·(1/4π)∫(∇′×ω)(x′)/\|x − x′\| d³x′; **sign = +1 is correct, sign = −1 reproduces the printed (5.14)**) · `gaussian_tube_fields(Gamma=1.0, sigma=0.1, L=4.0) -> dict(omega, curl_omega, bounds)` (a straight Gaussian tube along z of length L with a smooth end taper, for the quadrature tests) | (5.14) corrected, N25 | m/s | C07 (the discriminating test), E6 | §4 #22 (+ `gaussian_tube_fields` **NEW**) |
| 2.3 | `velocity_from_vorticity_fft(omega, Lbox) -> ndarray` (periodic box of side Lbox, n points per side; 2-D: ω (n, n) → (u, v) of shape (2, n, n) via ψ̂ = ω̂/k², u = ∂ψ/∂y, v = −∂ψ/∂x; 3-D: ω (3, n, n, n) → û = i k × ω̂/k²; k = 0 mode set to 0) | N24, Poisson | m/s | C07 (N24) | §4 #21 |
| 2.4 | `biot_savart_volume(omega, x, nodes, weights, eps=0.0) -> (3,)` ((5.16) as a vectorised quadrature sum; `omega` callable of x′) · `biot_savart_2d(omega_z, X, Y, x, eps=0.0) -> (2,)` (planar kernel e_z×r/(2π\|r\|²) summed over grid cells of area ΔxΔy; `X, Y` from `np.meshgrid(indexing="xy")`) | (5.16) | m/s | C07 (code + from-scratch), E6 | §4 #24 |
| 2.5 | `segment_induced_velocity(x, xa, xb, Gamma, eps=0.0) -> (3,) \| (3, N)` (closed form of (5.17) along a straight segment: (Γ/4πd)(cos θ_a − cos θ_b) along e_ω × e_d, vectorised over x) · `segment_speed(d, theta_a, theta_b, Gamma) -> float` (the scalar formula itself) | (5.17), D13, N29 | m/s | C08, E6 | §4 #25 (+ `segment_speed` **NEW**: E6's closed-form readout) |
| 2.6 | `filament_velocity(x, polyline, Gamma, closed=True, eps=0.0) -> (3,) \| (3, N)` (sum of 2.5 over the polyline's segments) · `filament_contributions(x, polyline, Gamma, closed=True) -> (3, M)` (one row per segment, for E6's per-segment arrows and inspector) | (5.17) | m/s | C08, E6 | §4 #25 (+ `filament_contributions` **NEW**) |
| 2.7 | `filament_preset(name, M=64, **p) -> ndarray (3, M+1 or M)` — `"straight"` (half_length=1.0: from (0, 0, −ℓ) to (0, 0, ℓ), open), `"semi_infinite"` (length=1000.0: from the origin up the z-axis, open), `"square"` (side=1.0, in z = 0, closed, M points on the perimeter), `"ring"` (R=0.5, in z = 0, closed, counterclockwise about +z), `"helix"` (R=0.3, pitch=0.4, turns=3, open) · `filament_velocity_preset(name, x, y, z, Gamma=1.0, M=64, component=None, **p) -> float \| (3,)` | E6 shapes | m, m/s | C08, E6 parity | §8 (+ the `_velocity_preset` wrapper **NEW**) |
| 2.8 | `ring_axis_velocity(z, R, Gamma) -> float` (ΓR²/(2(R² + z²)^{3/2}), along +z for a counterclockwise ring about +z) | C08 | m/s | C08, E6, B1 | §4 #25 |
| 2.9 | `ring_ring_velocity(R, z, R0, z0, Gamma0) -> (u_R, u_z)` (velocity at a point (R, z) of the meridional plane induced by a coaxial circular filament of radius R0 at z0 with circulation Γ0: Stokes stream function with `scipy.special.ellipk(m)`, `ellipe(m)`, **m = k² = 4RR0/((R + R0)² + (z − z0)²)**; guard m → 1) · `ring_self_velocity(R, a, Gamma, core="uniform") -> float` (Γ/(4πR)[ln(8R/a) − ¼]; `core="hollow"` uses ½ instead of ¼) | N42, N43 | m/s | C13, B1 | §4 #38 |
| 2.10 | `point_vortex_velocity(x, xv, Gamma, eps=0.0) -> (2,) \| (2, N)` (Σ_j (Γ_j/2π) e_z×(x − x_j)/\|x − x_j\|², self term excluded) · `point_vortex_rhs(xv, Gamma, boundary=None, **bp) -> (2, M)` (velocities of the vortices themselves, images included) | (5.2) superposed, C12 | m/s | C12, C13, E8 | §4 #35 (+ `point_vortex_rhs` **NEW**) |
| 2.11 | `point_vortex_evolve(xv0, Gamma, t_eval, boundary=None, rtol=1e-11, atol=1e-13, **bp) -> ndarray (T, 2, M)` (`solve_ivp` DOP853; `boundary` None, `"wall"` (`wall_y`=0.0), `"circle"` (`a`, `center`), `"channel"` (`H`)) | C12, C13 | m | C12, C13, E8 | §4 #35 |
| 2.12 | `point_vortex_invariants(xv, Gamma) -> dict(P, I, H)` (P = (ΣΓx, ΣΓy) (2,), I = ΣΓ\|x\|², H = −(1/2π)Σ_{j<k}Γ_jΓ_k ln\|x_j − x_k\|; unbounded plane only) · `centre_of_vorticity(xv, Gamma) -> (2,)` (P/ΣΓ; NaN if ΣΓ = 0) | D21, N40 | m³/s, m⁴/s, m⁴/s² | C12, E8 | §4 #35 (+ `centre_of_vorticity` **NEW**) |
| 2.13 | `point_vortex_preset(name) -> dict(xv, Gamma, boundary, bp, t_end, label)` — `"equal_pair"` (Γ = (1, 1), at (−½, 0), (½, 0)), `"unequal_pair"` (Γ = (1, 3), same places), `"opposite_pair"` (Γ = (1, −1), at (−½, 0), (½, 0): both move in +y at 0.159155 m/s, D21 step 6), `"near_wall"` (Γ = 1 at (0, 0.5), wall y = 0), `"knife_bucket"` (Γ = −1 at (0, 0.15) and +1 at (0, −0.15): the pair marches in −x; inside a circle a = 1 about (0.3, 0)), `"three_vortices"` (Γ = (1, 1, −0.5)) | E8 | m, m²/s | C12, C13, E8 | §8 (positions fixed here) |
| 2.14 | `wall_image_system(xv, Gamma, wall_y=0.0) -> (xv_all, Gamma_all)` (mirror points, strength −Γ) · `circle_image_system(xv, Gamma, a, center=(0.0, 0.0), inside=True, cylinder_circulation=0.0) -> (xv_all, Gamma_all)` (inverse points c + a²(x − c)/\|x − c\|² with −Γ; `inside=False` adds a centre vortex ΣΓ + cylinder_circulation — the exterior problem of Ch. 6) · `channel_image_velocity(h, H, Gamma, tol=1e-12, n_max=10**6) -> float` (drift of a vortex at height h between walls y = 0 and y = H by the image series summed until the increment < tol; closed form (Γ/4H)cot(πh/H) in the docstring as the test) | C13, N57, N58 | m, m/s | C13, E8 | §4 #37 (the `inside` switch **NEW**: the knife pair in a bucket is an interior problem, where no centre vortex is allowed) |
| 2.15 | `vortex_sheet_velocity(x, sheet_points, gamma, ds, eps=0.0) -> (2,) \| (2, N)` (N point vortices of strength γ·ds, counterclockwise positive) · `continuous_sheet_velocity(x, y, gamma, L=1.0) -> (u, v)` (flat finite sheet on −L/2 < x′ < L/2: u = −(γ/2π)[arctan((L/2 − x)/y) + arctan((L/2 + x)/y)], v = (γ/4π) ln(((x + L/2)² + y²)/((x − L/2)² + y²))) | N44, C14 | m/s | C14, E9 | §4 #39 (+ `continuous_sheet_velocity` **NEW**: the reference for the 1/N convergence) |

### C.3 `fluidpy/core/vortices.py` (additions)
| # | Callable (signature) | Implements (Eq.) | Returns | Used by | Status |
|---|---|---|---|---|---|
| 3.1 | `burgers_vortex(R, z, Gamma, alpha, nu) -> (u_R, u_phi, u_z, omega_z)` (u_R = −αR/2, u_z = αz, u_φ = (Γ/2πR)(1 − e^{−αR²/4ν}), ω_z = (αΓ/4πν)e^{−αR²/4ν}; book's α = 2 × Wikipedia's) · `burgers_vortex_field(Gamma, alpha, nu) -> u(x, t)` (Cartesian 3-D callable) · `burgers_core_radius(alpha, nu) -> float` (√(4ν/α)) | Ex. 5.12 form, N21, D18 | m/s, 1/s, m | C10, E5, E1 (straight tube) | §4 #20 (+ `burgers_core_radius` **NEW**) |
| 3.2 | `hill_spherical_vortex(R, z, A, a) -> (u_R, u_z, omega_phi)` (inside r ≤ a, frame moving with the vortex; ψ = (Aa⁴/10)(R²/a²)(1 − R²/a² − z²/a²), u_R = −R⁻¹∂ψ/∂z, u_z = R⁻¹∂ψ/∂R, ω_φ = AR) · `hill_spherical_vortex_field(A, a) -> u(x, t)` | N23 | m/s, 1/s | C06 (N23) | §4 #20 |
| 3.3 | `lamb_oseen_field(Gamma, nu, t0=0.0) -> u(x, t)` (thin wrapper of `core.navier_stokes.exact_field("lamb_oseen", …)` returning the velocity callable only) | (3.29) with σ_c² = 4ν(t + t₀) | m/s | C03, C06 | **NEW** (convenience; avoids unpacking in notebook cells) |
| 3.4 | `gaussian_tube_vorticity(R, z, Gamma=1.0, a0=0.1, L=1.0, twist=0.0) -> (omega_R, omega_phi, omega_z)` — a **narrowing Gaussian vortex tube**, exactly solenoidal: flux function χ(R, z) = (Γ/2π)(1 − e^{−R²/a(z)²}) with a(z)² = a0²e^{−z/L}; ω_z = (1/R)∂χ/∂R = (Γ/πa²)e^{−R²/a²}, ω_R = −(1/R)∂χ/∂z = −(Γ/2π)(R/(a²L))e^{−R²/a²}, ω_φ = twist·R·ω_z (twist [1/m] makes the lines helices) · `gaussian_tube_field(Gamma=1.0, a0=0.1, L=1.0, twist=0.0) -> omega(x, t)` (Cartesian callable) | E1, C01, D01 | 1/s | C01, E1 | §8 ("twisted_gaussian_tube_field", renamed and made exact here) |
| 3.5 | `gaussian_tube_section(z, R0=0.1, Gamma=1.0, a0=0.1, L=1.0) -> dict(radius, area, flux, mean_omega, peak_omega)` (the tube through (R0, 0): radius R0e^{−z/2L}, area πR0²e^{−z/L}, flux Γ(1 − e^{−R0²/a0²}) (independent of z), mean = flux/area, peak = Γ/(πa(z)²)) | (5.4), N03 | m, m², m²/s, 1/s | C01, E1 parity | **NEW** |
| 3.6 | `broken_tube_field(Gamma=1.0, a0=0.1, L=1.0) -> omega(x, t)` (ω = (Γ/πa0²)e^{−R²/a0²}e^{−z/L} e_z — **not a vorticity field** (∇·ω ≠ 0), flagged in its docstring; exists to show (5.4) failing) · `vortex_ring_vorticity(R, z, Gamma=1.0, R0=0.5, a=0.05) -> omega_phi` (Gaussian-core ring, for E1's closed-tube preset) | E1 presets | 1/s | C01, E1 | §8 |

### C.4 `fluidpy/core/integral_theorems.py` (addition)
| # | Callable (signature) | Implements (Eq.) | Returns | Used by | Status |
|---|---|---|---|---|---|
| 4.1 | `CurlCheck(volume, surface, diff)` · `curl_theorem_box(F_fn, bounds, n=24, rule="gauss", curl_fn=None) -> CurlCheck` (∫_V∇×F dV vs ∮_A n×F dA on a box, same field convention as `divergence_theorem_box`: `F_fn(X, Y, Z) -> (3, …)`; `curl_fn` exact or `fd_partials`) | (5.15), N27, N28 | (3,) each | C07 (N27, N28) | §4 #23 |

### C.5 `fluidpy/ch05_vorticity_dynamics.py` (chapter module; re-exports C.1–C.4)
| # | Callable (signature) | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 5.1 | `solid_body_from_vorticity(r, omega) -> (u_theta, omega_z)` (u_θ = ωr/2 via `solid_body_rotation(r, omega/2)`, ω_z = ω) · `line_vortex_gamma(r, Gamma) -> (u_theta, omega_z)` (Γ/2πr via `line_vortex(r, Gamma/2π)`; ω_z = 0 for r > 0, NaN at 0) | (5.1), (5.2) | m/s, 1/s | R02, R03, C02, E2 | §4 #1 |
| 5.2 | `solid_body_pressure_gradients(r, z, omega, rho=1000.0, g=9.81) -> (dp_dr, dp_dz)` (ρu_θ²/r = ρω²r/4 and −ρg) · `solid_body_pressure(r, z, omega, rho=1000.0, g=9.81, p_o=0.0) -> p` (p_o + ⅛ρω²r² − ρgz) · `line_vortex_pressure(r, z, Gamma, rho=1000.0, g=9.81, p_inf=0.0) -> p` (p_∞ − ρΓ²/(8π²r²) − ρgz) · `rankine_pressure(r, z, Gamma, a, rho=1000.0, g=9.81, p_inf=0.0) -> p` (funnel outside a, paraboloid inside with core vorticity Γ/πa², matched at a; p(0) − p_∞ = −ρΓ²/(4π²a²)) | (5.5a, b), (5.6), (5.7), R08 | Pa/m, Pa | C02, E2 | §4 #5, #8 |
| 5.3 | `isobar_height(r, param, dp_over_rho_g=0.0, kind="solid", g=9.81) -> z` (`"solid"`: ω²r²/8g − Δp/ρg with param = ω; `"line"`: −Γ²/(8π²r²g) − Δp/ρg with param = Γ) · `rotating_tank_free_surface(R, depth, Omega_tank, g=9.81, H=None) -> dict(z_vertex, z_rim, rim_rise, centre_drop, r_dry, uncovered_area, spills)` (volume-conserving paraboloid z = z_v + Ω_tank²r²/2g; if z_v < 0 the bottom uncovers: `brentq` for the dry radius with the volume constraint; `spills` if H given and z_rim > H) | (5.6), Fig. 5.2 | m, m² | C02, E2 | §4 #5 (`H` optional, `rim_rise`, `centre_drop` **NEW**) |
| 5.4 | `bernoulli_across_vortex(kind, r, z=0.0, param=1.0, rho=1000.0, g=9.81, r_ref=0.0, a=None) -> float` (B(r) − B(r_ref) with B = ½u_θ² + gz + p/ρ from 5.2: `"solid"` ω²(r² − r_ref²)/4, `"line"` 0 (r_ref > 0 required), `"rankine"` jumps only inside the core) | R07, (4.69) | m²/s² | C02, E2 | §4 #6 |
| 5.5 | `polar_viscous_stress(u_r_expr, u_theta_expr, r, theta, mu) -> dict(sigma_rr, sigma_rtheta, sigma_thetatheta)` (sympy, 2μS in polar coordinates via `core.curvilinear.strain_rate`) · `line_vortex_viscous_stress(r, Gamma, mu) -> float` (−μΓ/πr²) · `vortex_stress_force(kind, r, mu=1e-3, **p) -> dict(sigma_rtheta, force_divergence, force_laplacian, force_curl)` (θ-components [N/m³] of the net viscous force by the three routes of D03: (1/r²)∂(r²σ_rθ)/∂r, μ(u_θ″ + u_θ′/r − u_θ/r²), −μ(∇×ω)_θ = μ∂ω_z/∂r; kinds `"solid"` (omega), `"line"` (Gamma), `"gaussian"` (Gamma, sigma_c)) | N06, N09, D03 | Pa, N/m³ | C02, E2 | §4 #7 (+ `vortex_stress_force` **NEW**: the N09 table and D03's three routes in one call) |
| 5.6 | `rotating_cylinder_flow(r, a, omega) -> (u_theta, omega_z)` (Rankine with Γ = πa²ω; cylinder turns at ω/2) · `torque_per_length(r, Gamma, mu) -> float` (2πr²σ_rθ = −2μΓ, independent of r) · `dissipation_outside_cylinder(a, R_out, Gamma, mu) -> dict(dissipated, power_in, power_out)` (∫_a^{R_out}2πrρε dr by `quad` with `core.constitutive.dissipation_rate`; power = \|torque\| × angular speed u_θ/r = μΓ²/πr²) | R08, N08 | m/s, N m/m, W/m | C02 | §4 #9, #10 |
| 5.7 | `vortex_pressure_scenario(kind, r, z=0.0, rho=1000.0, mu=1e-3, g=9.81, **p) -> dict(u_theta, p, B, sigma_rtheta, net_viscous_force, centripetal, dp_dr, surface_z)` — kinds `"solid"` (omega), `"line"` (Gamma), `"rankine"` (Gamma, a), `"cylinder"` (omega, a); `centripetal` = u_θ²/r, `dp_dr`/ρ equals it (the (5.5a) balance), `surface_z` = the free-surface height with p = 0 gauge at the reference point (solid: r = 0; line/Rankine: r → ∞) | E2 | mixed SI | C02, E2 parity | §8 (keys fixed here) |
| 5.8 | `lamb_oseen_circulation(r, t, Gamma, nu, t0=0.0) -> (Gamma_r, dGamma_dt)` (Γ(1 − e^{−r²/4ν(t+t₀)}) and its time derivative) · `lamb_oseen_viscous_loop_integral(r, t, Gamma, nu, t0=0.0) -> float` (∮ν∇²u·dx on the circle, = dΓ/dt: the (5.11) check) | N13 | m²/s, m²/s² | C03, E3 | §4 #13 (+ the loop integral **NEW**) |
| 5.9 | `kelvin_hypotheses(inviscid, barotropic, conservative, inertial) -> dict(holds, surviving_terms, verdict, text)` (`surviving_terms` ⊂ {"viscous", "baroclinic", "body", "coriolis"}; `text` names each surviving term with its TeX) · `kelvin_hypotheses_text(...) -> str` (the verdict sentence, for exact-text parity) | N16, (5.10), (5.11) | – | C03 table, E3 status | §4 #15 (+ `_text` **NEW**) |
| 5.10 | `lock_exchange_initial_vorticity_rate(rho1, rho2, delta, g=9.81) -> float` (2(ρ₂ − ρ₁)g/((ρ₂ + ρ₁)δ), counterclockwise positive with ρ₂ on the **left**, as in Fig. 5.5) · `lock_exchange_fields(rho1, rho2, delta, g=9.81, H=1.0, p_top=0.0) -> dict(rho, p, grad_rho, grad_p)` (callables of (x, y): ρ = ρ̄ − (Δρ/2)tanh(2x/δ), hydrostatic p with the mean density ρ̄ = (ρ₁ + ρ₂)/2 at t = 0⁺; y up, interface at x = 0) · `baroclinic_rate_2d(grad_rho, grad_p, rho) -> float` ((∇ρ×∇p)_z/ρ² for 2-vectors) | N15, D07, (5.28) | 1/s² | C04, E4 | §4 #14 (+ `baroclinic_rate_2d` **NEW**: E4's formula readout) |
| 5.11 | `rotating_ns_residual(u, p, rho, x, t=0.0, nu=1e-6, Omega=(0.0, 0.0, 0.5), g_eff=(0.0, 0.0, -9.81), h=1e-4, ht=None) -> (3,)` ((5.20) residual) · `absolute_vorticity(omega_rel, Omega) -> (3,)` (ω + 2Ω) | R14, R18 | m/s², 1/s | C09, C11 | §4 #27, #30 |
| 5.12 | `diffusing_vortex_sheet(y, t, gamma, nu) -> (u, omega_z)` (u = −(γ/2)erf(y/2√(νt)), ω_z = (γ/(2√(πνt)))e^{−y²/4νt}) | N22 | m/s, 1/s | C06, E5 (optional) | §4 #20 |
| 5.13 | `strain_preset(name, rate=1.0) -> (3, 3)` (`"axial_stretch"` diag(−r/2, −r/2, r); `"shear_tilt"` G[2, 0] = r (w = r x); `"planar"` diag(r, −r, 0); `"burgers"` = axial_stretch) · `uniform_strain_vorticity(omega0, G, t, rate=1.0) -> (3,)` (Cauchy's frozen-in solution ω(t) = e^{Gt}ω₀ for a weak vortex line in a steady linear flow, inviscid; `G` an array or a preset name) · `stretched_tube(L, L0, omega0, A0) -> dict(omega, A, Gamma)` (inviscid: A = A₀L₀/L, ω = ω₀L/L₀, Γ = ω₀A₀) | (5.32), D17, D18 | 1/s | C10, E5 parity | §4 #32 (+ `strain_preset` **NEW**) |
| 5.14 | `helix(s, a=1.0, c=0.3) -> (3,)` (arc-length parametrised helix of radius a and pitch 2πc) · `helix_frame(s, a=1.0, c=0.3) -> dict(e_s, e_n, e_m, curvature, torsion)` (**book convention: e_n points away from the centre of curvature** (−Frenet N), e_m = e_s × e_n; κ = a/(a² + c²), τ = c/(a² + c²)) | Fig. 5.9, N53, D16 | m, 1/m | C10, E5 | §8 (renamed from `helical_vortex_line`, `frenet_frame`) |
| 5.15 | `column_relative_vorticity(h, h0, zeta0=0.0, f=1.0e-4) -> float` (ζ = (ζ₀ + f)h/h₀ − f, f = 2Ω_z the local planetary vorticity) · `column_over_slope(x, depth="ridge", lat_deg=45.0, zeta0=0.0, h0=1000.0, bump=0.2, width=1.0e5, Omega=OMEGA_EARTH) -> dict(h, zeta, f)` (depth profile preset `"ridge"` / `"trough"` / `"slope"`, or a callable h(x)) · `relative_circulation_after_move(Gamma0, A0, lat0_deg, A1, lat1_deg, Omega=OMEGA_EARTH) -> float` (Γ₁ = Γ₀ + 2Ω(A₀ sin φ₀ − A₁ sin φ₁)) | D20, N54, (5.33) | 1/s, m²/s | C11, E7 parity | §4 #33, #34, §8 (signature: f instead of Ω; latitude in degrees at the interface, `_deg`) |
| 5.16 | `helical_swirl_field(a=1.0) -> u(x, t)` (u_φ = aRz: ω_R = −aR, ω_z = 2az, vortex lines zR² = const — the N02 test field) · `cellular_flow(x, t=0.0, U=1.0, ell=1.0) -> (2,)` (ψ = Uℓ sin(x/ℓ) sin(y/ℓ)) · `abc_flow(A=1.0, B=1.0, C=1.0) -> u(x, t)` (u = (A sin z + C cos y, B sin x + A cos z, C sin y + B cos x): a steady Euler (Beltrami) flow with ω = u) | N02, C03, C05 | m/s | C01, C03, C05, C06, C09 | **NEW** (named test fields for the notebook; ABC is the inviscid field of the Helmholtz demos) |
| 5.17 | `vortex_pair(Gamma1, Gamma2, h) -> dict(V1, V2, centre_from_1, rotation_rate, period, translation_speed)` (V₁ = Γ₁/2πh at vortex 2, V₂ = Γ₂/2πh at vortex 1, centre h₁ = Γ₂h/(Γ₁ + Γ₂), rate (Γ₁ + Γ₂)/2πh², period 2π/rate; Γ₁ + Γ₂ = 0: centre `inf`, rate 0, translation Γ₁/2πh) · `vortex_near_wall_speed(Gamma, h) -> float` (Γ/4πh) | D21, D22 | m/s, m, rad/s, s | C12, C13, E8 parity | §4 #36, #37 |
| 5.18 | `ring_dynamics(rings0, t_eval, wall_z=None, core="uniform") -> dict(t, R, z, a, impulse)` (coaxial thin-core rings, each dict(R, z, Gamma, a); self speed 2.9 + mutual 2.9 + images of every ring in the wall z = wall_z; core volume a²R conserved; `impulse` = ΣΓπR² per time) | N42, N43, Fig. 5.15 | m, m⁴/s | C13, B1 | §4 #38 |
| 5.19 | `vortex_sheet_strength(u_above, u_below, convention="ccw") -> float` (`"ccw"`: γ = u_below − u_above = u₂ − u₁ (the book's text); `"caption"`: u₁ − u₂ (Fig. 5.16's caption, the clockwise magnitude)) · `discrete_sheet_u(x, y, gamma, N, L=1.0) -> float` (u of a finite row of N filaments at the cell midpoints of [−L/2, L/2]) · `discrete_sheet_convergence(gamma, N_list=(10, 100, 1000), L=1.0, y_max=0.1) -> dict(N, l1_error, u_above, u_below)` (L1 error ∫\|u_N(0, y) − u_cont(0, y)\|dy over \|y\| < y_max; ∝ 1/N) | D23, N44, N60 | m/s | C14, E9 parity | §4 #39 (+ `discrete_sheet_u` **NEW**) |
| 5.20 | `sheet_rollup(N=100, gamma=1.0, amplitude=0.01, delta=0.05, t_eval=None, L=1.0) -> dict(t, x, y)` (periodic vortex sheet of period L, N points, Krasny δ-smoothed periodic kernel u − iv = (Γ_j/2iL)cot(π(z − z_j)/L) with δ; initial y = amplitude·sin(2πx/L); DOP853; arrays (T, N)) | C14 animation, E9 roll-up | m | C14, E9 | §8 |

### C.6 `scripts/ch05_drawings.py` (drawing helpers, no physics; `from scripts.ch05_drawings import …`)
`tube_mesh(ax3d, tube, color)` (the Fig. 5.1 analogue from a `Tube`), `tank_section(ax, R, H)` (walls of a tank in the
(r, z) half-plane), `disc_element(ax, center, radius, isobars, isopycnals)` (the Fig. 5.6 analogue: rim pressure arrows,
G, the net-force line, a curved torque arrow), `loop_with_element(ax, pts, k)` (a loop with dx and du arrows at point k,
Fig. 5.4 analogue), `biot_savart_geometry(ax3d, polyline, x)` (Fig. 5.8 analogue), `helix_with_frame(fig, s0)` (plotly,
Fig. 5.9 analogue), `column_sketch(ax, x, h, zeta)` (Fig. 5.10 analogue), `vortex_marks(ax, xv, Gamma)` (teal/rose
circular arrows), `sheet_circuit(ax, ds, dn)` (Fig. 5.16 analogue). Figures → `outputs/ch05/`.

**Count:** C.1 23 rows (30 functions + 7 result classes) · C.2 15 rows (20) · C.3 6 rows (11) · C.4 1 row (1 + 1
class) · C.5 20 rows (46) · C.6 9 drawing helpers — **108 public functions and 8 result classes** in 2 new core modules,
2 extended core modules and the chapter module. **Beyond analysis §4 and curation §8 (NEW, flagged):** `loop_circulation`, `loop_length`,
`kelvin_scenario_circulation`, `kelvin_scenario_rate`, `vorticity_budget_preset` (C.1) · `gaussian_tube_fields`,
`segment_speed`, `filament_contributions`, `filament_velocity_preset`, `point_vortex_rhs`, `centre_of_vorticity`, the
`inside` switch of `circle_image_system`, `continuous_sheet_velocity` (C.2) · `burgers_core_radius`, `lamb_oseen_field`,
`gaussian_tube_section` (C.3) · `rim_rise`/`centre_drop`, `vortex_stress_force`, `lamb_oseen_viscous_loop_integral`,
`kelvin_hypotheses_text`, `baroclinic_rate_2d`, `strain_preset`, `helical_swirl_field`, `cellular_flow`, `abc_flow`,
`discrete_sheet_u` (C.5). **Signatures fixed here where the plans left them open:** the axisymmetric
`tube_flux_budget`, the Gaussian tube's flux function, `pressure_torque_on_element`'s keyword form and keys,
`kelvin_scenario` geometry and the rotating scenario, `point_vortex_preset` positions, `vortex_pressure_scenario` keys,
`helix_frame`'s e_n sign (book convention), `column_relative_vorticity` taking f, latitudes in degrees at the interface.

**Contract numbers the notebook and explainers rely on (computed independently for this design, 2026-09-23 — please
keep):**
- C01 · `gaussian_tube_section(0.0)` (R0 = a0 = 0.1, Γ = 1, L = 1) → radius 0.1, area 0.031416 m², flux 0.632121 m²/s,
  mean 20.121 s⁻¹, peak 31.831 s⁻¹ · `gaussian_tube_section(1.0)` → radius 0.060653, area 0.011557 m², flux 0.632121,
  mean 54.695 s⁻¹, peak 86.526 s⁻¹ · `tube_flux_budget("gaussian_tube", 0.0, 1.0, 0.1)` → lower −0.632121, side ≈ 0,
  upper +0.632121, total ≤ 1e-8 · `tube_flux_budget("broken", 0.0, 1.0, 0.1)` → lower −0.632121, side 0, upper
  +0.232544 (= 0.632121·e⁻¹), total −0.399577 · Lamb–Oseen tube strength at r = σ_c: 1 − e⁻¹ = 0.632121 of Γ by both
  routes.
- C02 · `solid_body_pressure(0.1, 0.0, 10.0)` = 125.0 Pa · rim–centre height ω²R²/8g at ω = 10 s⁻¹, R = 0.1 m = 0.012742 m
  · `rotating_tank_free_surface(0.1, 0.1, 5.0)` → z_vertex 0.093629 m, z_rim 0.106371 m (rise = drop = 6.371 mm) ·
  `line_vortex_pressure(0.1, 0.0, 1.0)` = −1266.515 Pa · `isobar_height(0.1, 1.0, kind="line")` = −0.129104 m ·
  `rankine_pressure(0.1, 0.0, 1.0, 0.1)` = −1266.515 Pa, `rankine_pressure(0.0, 0.0, 1.0, 0.1)` = −2533.030 Pa ·
  `bernoulli_across_vortex("solid", 0.1, param=10.0)` = 0.25 m²/s² · `line_vortex_viscous_stress(0.1, 1.0, 1e-3)` =
  −0.0318310 Pa · `torque_per_length(r, 1.0, 1e-3)` = −0.002 N m/m for every r · `dissipation_outside_cylinder(0.1, 1.0,
  1.0, 1e-3)` → dissipated 0.0315127, power_in 0.0318310, power_out 3.18310e-4 W/m (in − out = dissipated) · tornado
  preset (air ρ = 1.2, Γ = 10⁴ m²/s, a = 50 m): u_max 31.831 m/s, core-edge deficit 607.93 Pa, centre deficit 1215.85 Pa.
- C03 · `lamb_oseen_circulation(0.005, 10.0, 0.01, 1e-6)` → (0.00464739, −3.34538e-4) · `kelvin_scenario_circulation(
  "cellular", 0.0)` = 1.155130 m²/s (flux route: ∫2 sin x sin y dA over the disc, same number) and constant to 1e-9 up to
  t_end · `kelvin_scenario_circulation("rankine_straddle", 0.0)` = 1.098212 m²/s (2 × the lens area 0.549106 m²; the kink at the core edge makes a 256-point loop sum accurate to 1e-5 only) ·
  `kelvin_scenario_circulation("rotating", 5.0)` = 1.985865 m²/s while Γ_a = π = 3.141593 m²/s ·
  `kelvin_hypotheses(True, True, True, True)["holds"]` = True, `(True, False, True, True)["surviving_terms"]` =
  ["baroclinic"].
- C04 · `lock_exchange_initial_vorticity_rate(1000.0, 1025.0, 0.1)` = 2.422222 s⁻² · `pressure_torque_on_element()`
  (defaults: R = 0.01 m, ∇ρ = (10, 0) kg/m⁴, ∇p = (0, −9810) Pa/m, ρ₀ = 1000) → x_G = (2.5000e-7, 0) m, I_G = 1.57079e-5
  kg m²/m, torque −7.7047e-7 N m/m, spin_up −0.0981000001 s⁻², baroclinic −0.0981 s⁻², ratio 1 + 1.25e-9 (exact:
  1/(1 − \|∇ρ\|²R²/8ρ₀²)) · `baroclinic_rate_2d([10, 0], [0, -9810], 1000)` = −0.0981.
- C06/C10 · `burgers_core_radius(1.0, 1e-6)` = 0.002 m · `burgers_vortex(0.0, 0.0, 1e-3, 1.0, 1e-6)[3]` = 79.5775 s⁻¹ ·
  `diffusing_vortex_sheet(0.0, 1.0, 1.0, 1e-6)[1]` = 282.095 s⁻¹ and `(1e-3, 1.0, 1.0, 1e-6)[0]` = −0.260250 m/s ·
  `uniform_strain_vorticity([0, 0, 1], "axial_stretch", 2.0)[2]` = 7.389056 · `uniform_strain_vorticity([1, 0, 0],
  "shear_tilt", 2.0)[2]` = 2.0 · `uniform_strain_vorticity([0, 0, 1], "planar", 2.0)[2]` = 1.0 · `stretched_tube(2.0,
  1.0, 10.0, 1e-4)` → ω 20.0, A 5e-5, Γ 1e-3 · `helix_frame(0.0)` → κ 0.917431, τ 0.275229 (a = 1, c = 0.3).
- C07/C08 · Gaussian tube (Γ = 1, σ_c = 0.1) at r = 0.5: u_θ = 0.308838 m/s for the contract's finite tube L = 4 (`gaussian_tube_fields(...)['u_theta_reference']`; the infinite-line value is 0.318310, `u_theta_infinite`) (+1/(4π)); the printed −1/(4π) gives
  −0.318310 · `segment_induced_velocity([1, 0, 0], [0, 0, -1], [0, 0, 1], 1.0)[1]` = 0.112540 m/s · infinite line
  1/2π = 0.159155, semi-infinite 0.0795775 · square loop of side 1 at its centre 2√2/π = 0.900316 · `ring_axis_velocity(
  0.0, 0.5, 1.0)` = 1.0, `(0.5, 0.5, 1.0)` = 0.353553 · polygon ring (R = 0.5, M segments) at its centre: M = 8 → 1.054786,
  16 → 1.013052, 64 → 1.000804 (order 2 in 1/M) · `ring_self_velocity(1.0, 0.1, 1.0)` = 0.328816 m/s.
- C09/C11 · 2Ω∂w/∂z with Ω = 7.292115e-5, ∂w/∂z = 1e-5 → 1.45842e-9 s⁻² · `column_relative_vorticity(1100.0, 1000.0,
  0.0, 1e-4)` = 1.0e-5 s⁻¹ · `relative_circulation_after_move(0.0, π(5e5)², 30.0, π(5e5)², 60.0)` = −4.19261e7 m²/s
  (mean ζ = −5.33820e-5 s⁻¹) · f(30°) = 7.29212e-5, f(45°) = 1.03126e-4, f(60°) = 1.26303e-4 s⁻¹.
- C12/C13 · `vortex_pair(1, 1, 1)` → V1 = V2 = 0.159155, centre_from_1 0.5, rotation_rate 0.318310 rad/s, period
  19.7392 s · `vortex_pair(1, 3, 1)` → centre_from_1 0.75, rate 0.636620 · `vortex_pair(1, -1, 1)` → translation 0.159155,
  rate 0 · `vortex_near_wall_speed(1.0, 0.5)` = 0.159155 · `channel_image_velocity(0.25, 1.0, 1.0)` = 0.25 m/s (direct
  sum over ±20 000 image pairs 0.250002) · `point_vortex_invariants([[0, 1], [0, 0]], [1, 1])["H"]` = 0.
- C14 · `vortex_sheet_strength(-1.0, 1.0)` = 2.0 (ccw), `convention="caption"` → −2.0 · `discrete_sheet_convergence(2.0)`
  → l1_error 0.043968 (N = 10), 0.0044131 (100), 0.00044128 (1000) — slope −1 · `continuous_sheet_velocity(0.0, 0.05,
  2.0)[0]` = −0.936549 and `discrete_sheet_u(0.0, 0.05, 2.0, 100)` = −0.936551.


**Implementation cross-check (probe of the implementer's WIP, 2026-09-23).** `fluidpy/core/vorticity.py` (1282 lines),
`fluidpy/core/biot_savart.py` (795) and `fluidpy/ch05_vorticity_dynamics.py` (1119), with the additions to
`core/vortices.py` and `core/integral_theorems.py`, already define 112 of the 113 Part C names, and all 35 contract
numbers probed above evaluate as written (tube sections and budgets, tank and funnel pressures, stress and dissipation,
Lamb–Oseen and the four Kelvin scenarios, the lock-exchange rate, the element spin-up, Burgers, strain presets, the
helix frame, segments, square and polygon rings, ring self-speed, column and ring-of-air circulation, pairs, wall and
channel drifts, sheet strength, discrete sheet and its L1 convergence, and two budget presets). **Two gaps for the
implementer:** `abc_flow(A=1.0, B=1.0, C=1.0)` (C.5 row 5.16) and the `"helmholtz_abc"` scenario of `kelvin_scenario`
(C.1 row 1.10) were added by this design after the WIP was started — C05, C06 (frozen-in) and C09 (N30) call them (the
WIP's own `"helmholtz"` scenario is kept for E3); and the `"cellular"` scenario's default `t_end` is **6** turnover
times (37.7 s) here, 3 in the WIP — the loop grows only ×3.5 in 3 turnovers, ×6.7 in 6, and the notebook and E3 quote
the 6-turnover numbers.
The probe also corrected one design number: the Rankine-straddle loop's Γ(0) is 1.098212 m²/s (exact lens area), which
the WIP already returns.

---

## Part A — notebook storyboard (`notebooks/build_ch05.py` → `notebooks/ch05_vorticity_dynamics.ipynb`)

**Section check (book order; one `nb.section` per book section; cell numbers are estimates ± 5).**
- §5.1 → R01, C01 (N01 N02 N03 N45 · D01 · **E1**), R02 R03 R04 R05 R06 R07 R08, C02 (N04 N05 N06 N07 N08 N09 N46 N47 · D02 D03 · **E2**) — cells 8–92
- §5.2 → R09, C03 (N10 N11 N12 N13 N14 N16 N48 · D04 D05 · **E3**), C04 (N15 N49 N50 · D06 D07 · **E4**) — cells 93–170
- §5.3 → C05 (N17 N51 · D08; pointer to E3's Helmholtz mode) — cells 171–190
- §5.4 → R10 R11, C06 (N18 N19 N20 N22 N23 · D09 ★★★; C05's field-equation check at its end; pointer N21 → C10) — cells 191–240
- §5.5 → C07 (N24 N25 N26 N27 N28 N52 · D10 ★★★ D11 ★★★), C08 (N29 · D12 D13 · **E6**) — cells 241–318
- §5.6 → R12 R13 R14 R15 R16 R17, C09 (N30–N37 · D14 D15 ★★★), C10 (N38 N21 N53 · D16 D17 D18 · **E5**), R18, C11 (N39 N54 · D19 D20 · **E7**) — cells 319–440
- §5.7 → C12 (N40 N41 N55 N56 · D21), C13 (N42 N43 N57 N58 N59 · D22 · **E8**) — cells 441–500
- §5.8 → C14 (N44 N60 · D23 · **E9**) — cells 501–530
- end → S01, S02, summary — cells 531–534

Every CORE block below has the nine parts in order (question · idea · primers · maths with its D rows · tiny example ·
code + "What does the code above do?" · from-scratch check where curation §7 asks · ≥ 1 visual · notes + "What would
change if…"). *code:* is the intent (exact `ch05.` calls and arguments; every line gets a novice comment in the builder);
*expect:* the numbers the executed cell must print (computed for this design — the builder re-checks them against the
implemented `fluidpy`); *explain:* the numbered "What does the code above do?" list. Markdown drafts are in our words;
every equation is written out in LaTeX next to its number (convention 9). Derivation cells are copied from Part F word
for word. B1 (`vortex_rings`) is the backup and is **not** embedded (9 explainers embedded; 5–10 allowed). Primer titles
are the exact Part E Concept texts. **Equations in drafts:** where a draft cites a book equation by number in a markdown cell
that would otherwise show no maths (section intros, see/read/change texts, explainer `why`/`tries`, pointers), the builder
runs it through `show_eqs(text, EQ)`, which writes the equation in LaTeX after its number; ASCII idea sketches in code
fences get one line **below** the fence writing each cited equation in LaTeX ("Here (5.4) is $\int_V\nabla\cdot\boldsymbol\omega\,dV=-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0$."),
because TeX inside a fence does not render. `tools/coverage_check.py` rule 8 and `tools/eq_refs.py` must both pass.

### A.0 Front matter
1. `nb.title(big_idea="Vorticity — twice the local spin of the fluid — is the quantity that makes a flow swirl, and this
   chapter follows it: where it lives (vortex lines and tubes, which cannot end: $\int_V\nabla\cdot\boldsymbol\omega\,dV=-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0$
   *(Eq. 5.4)*), what it does to pressure (the bowl of a spinning bucket, the funnel of a drain), when it is conserved
   (Kelvin: $D\Gamma/Dt=0$ *(Eq. 5.8)*; Helmholtz: vortex lines move with the fluid), how it changes (it is carried,
   stretched, tilted and diffused: $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$ *(Eq. 5.13)*; on a
   rotating planet with density contrasts it gains the planet's vorticity and a baroclinic source,
   $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u+\frac1{\rho^2}\nabla\rho\times\nabla p+\nu\nabla^2\boldsymbol\omega$
   *(Eq. 5.30)*), and how it makes velocity everywhere (Biot–Savart,
   $\mathbf u=\frac1{4\pi}\int\frac{\boldsymbol\omega\times(\mathbf x-\mathbf x')}{\lvert\mathbf x-\mathbf x'\rvert^3}d^3x'$ *(Eq. 5.16)*) — so that
   vortices push each other, slide along walls and roll up.", roadmap=["§5.1 vortex lines and tubes; the two basic
   vortices, their pressure (5.6), (5.7) and their viscous stress", "§5.2 Kelvin's circulation theorem (5.8)–(5.11); the
   three sources of vorticity; barotropic vs baroclinic", "§5.3 Helmholtz's four vortex theorems", "§5.4 the vorticity
   equation (5.13)", "§5.5 velocity from vorticity: Poisson, Green, Biot–Savart (5.14)–(5.17)", "§5.6 the vorticity
   equation in a rotating frame (5.30); stretching and tilting (5.32); absolute circulation (5.33)", "§5.7 point vortices,
   images, rings", "§5.8 vortex sheets"], prerequisites=["vorticity ω = ∇×u, spin ½ω, circulation Γ = ∮u·dx and the
   solid-body, line, Rankine and Gaussian vortices (Ch. 3 §3.4–3.5)", "Gauss' and Stokes' theorems, the ε–δ identity,
   the curl of a gradient (Ch. 2)", "Navier–Stokes, Euler, the viscous force −μ∇×ω, the Lamb identity, the Bernoulli
   function, rotating frames and Boussinesq (Ch. 4)", "hydrostatics (Ch. 1)"])`
2. `nb.explainer_index([("vortex_tubes_cannot_end", "Can a vortex just stop in the middle of the water?", "C01: vortex
   lines and tubes; the strength is the same at every section (5.4)"), ("vortex_pressure_funnel", "Why is a spinning
   bucket a bowl and a drain a funnel?", "C02: pressure (5.6), (5.7); stress without net force"), ("kelvin_material_loop",
   "Stretch and tangle a loop of dye — what stays the same?", "C03 C05: Kelvin's theorem and its hypotheses; Helmholtz"),
   ("baroclinic_torque", "How can density make fluid spin?", "C04: the pressure torque = ∇ρ×∇p/ρ²; the lock exchange"),
   ("vorticity_stretching_tilting", "How can a flow spin fluid faster without a torque?", "C06 C10: (5.13), stretching and
   tilting (5.32), Burgers' vortex"), ("biot_savart_filament", "How does a vortex push water far away?", "C07 C08:
   Biot–Savart (5.16), the filament law (5.17), the sign of (5.14)"), ("vorticity_equation_rotating", "What does a spinning
   planet add to the vorticity equation?", "C09 C11: (5.30), the fluid column, absolute circulation (5.33)"),
   ("point_vortex_lab", "A vortex cannot push itself — so how do vortices move?", "C12 C13: pairs, the centre of
   vorticity, images"), ("vortex_sheet_rollup", "What is a velocity jump made of?", "C14: sheet strength = jump; roll-up")])`
   (the index text is plain; the explainers show every equation in full).
3. `nb.setup()`, then `nb.code` (no CORE yet) — *code:* `import numpy as np; import sympy as sp; import
   matplotlib.pyplot as plt; import plotly.graph_objects as go` · `from fluidpy import ch05_vorticity_dynamics as ch05` ·
   `from fluidpy import ch01_introduction as ch01` · `from fluidpy.core import vortices, kinematics, navier_stokes as ns,
   integral_theorems as itg, curvilinear as cu, rotating` · `from fluidpy.core.interact import slider_figure` ·
   `from fluidpy.core.anim import animate` · `from fluidpy.core.style import COLORS, savefig` · `from tools.convergence
   import observed_order`. *expect:* no output. *explain:* 1. numerical, symbolic and plotting libraries; 2. the chapter
   module — every function of this notebook lives there (or in the core modules it re-exports) and is tested in
   `tests/test_ch05.py`; 3. Ch. 1–4 machinery reused (vortex profiles, stencils, loops, Appendix-B operators, rotating
   frames); 4. house helpers for sliders, animations and colours.
4. `nb.md` — "**Notation, conventions and colours used in this notebook.**" Table | symbol | meaning | unit |:
   ω vorticity [1/s] (never a rotation rate in this chapter); u_θ, u_r polar velocity [m/s]; Γ circulation [m²/s]; γ
   vortex-sheet strength [m/s]; ω_z plane vorticity, counterclockwise positive; p pressure [Pa]; ρ density [kg/m³];
   μ, ν viscosities; σ_ij viscous stress [Pa]; Ω frame rotation [rad/s]; 2Ω planetary vorticity, ω + 2Ω absolute
   vorticity; f = 2Ω sin φ; x′ (with ∇′, V′, A′) the source point in §5.5; (s, n, m) natural coordinates on a vortex line.
   Then four ⚠️ lines: "**ω is the vorticity.** In (5.1), $u_\theta=\omega r/2$, the fluid turns at ω/2: a tank spinning at
   1 rad/s has ω = 2 s⁻¹. Ch. 3's `solid_body_rotation(r, omega0)` takes the turning rate ω₀ = ω/2." · "**Γ and γ.** Γ
   (`Gamma`, m²/s) is always a circulation; the strength of a vortex sheet is γ (`gamma`, m/s) — the book writes Γ for
   it in §5.8." · "**σ and ε.** σ is the viscous stress here (Ch. 3's Gaussian core radius is written σ_c); ε is the
   Levi-Civita symbol except in N08, where ε is the dissipation rate." · "**Primes in §5.5** mark the source point x′
   being integrated over, not a moving frame." Colours (convention 8) as a one-line legend.
5. `nb.md` — "**Where this chapter is used later**" (two columns): (5.4) and Helmholtz → trailing vortices of wings
   (Ch. 14) · Kelvin (5.8) → irrotational flow stays irrotational (Ch. 6), the starting vortex and lift (Ch. 6, 14) ·
   $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$ (5.13) → vortex decay (Ch. 8), vorticity–stream-function
   CFD (Ch. 10), the turbulent cascade (Ch. 12) · Biot–Savart (5.16), (5.17) → vortex methods (Ch. 10), lifting line and
   induced drag (Ch. 14) · (5.30), (5.33), the column (ζ + f)/h → potential vorticity, Rossby waves, baroclinic
   instability (Ch. 13) · vortex sheets → Kelvin–Helmholtz instability (Ch. 11), wakes (Ch. 14). Climate hook: "Ch. 13's
   potential-vorticity conservation is the fluid column of C11 plus stratification; the sea breeze and every front start
   with the baroclinic torque of C04."
6. `nb.md` — "**Book slips we correct in this notebook** (each is explained where it is used):" a two-column table
   | the book prints | we use |: (5.14) $\mathbf u=-\frac1{4\pi}\int\frac{\nabla'\times\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}d^3x'$ | the
   factor +1/(4π) (D10) · the rewrite's $+\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}\times\boldsymbol\omega$ | $-\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}\times\boldsymbol\omega$
   (D11; the two slips cancel, so (5.16) is right) · Fig. 5.2's rotation label "2ω" | the tank turns at ω/2 (R02) · "ρ
   and p single valued" as the reason ∮dp/ρ = 0 | ρ = ρ(p) makes dp/ρ exact (R09) · "C in irrotational fluid ⇒ no
   viscous term" | only for incompressible, constant-μ flow (N13) · the Lamb step's ∇(u·u) | ∇(½u·u) (R11) · (5.26)'s Π
   | Φ (N33) · (5.27) drops $u_{j,j}(\omega_n+2\Omega_n)$ silently | zero by continuity, said aloud (D15) · "Exercise 5.8"
   after (5.14), "Exercise 5.11" after (5.33) | Exercises 5.9 and 5.10 · Fig. 5.16's caption $d\Gamma/ds=u_1-u_2$ | the text's
   $u_2-u_1$, counterclockwise positive (D23) · Fig. 5.11's caption: G is where the induced velocity vanishes | G is the
   centre of vorticity; the fluid there is at rest only for equal strengths (D21) · "hyperboloids of the second degree"
   (Fig. 5.3) | cubic surfaces (c − z)r² = const (N07).

### A.1 §5.1 Introduction — R01, C01 (+N01 N02 N03 N45 · D01 · E1), R02–R08, C02 (+N04–N09 N46 N47 · D02 D03 · E2)
1. `nb.section("5.1", "Introduction", intro="**What is this section about?** Vorticity, ω = ∇×u, is twice the spin of a
   fluid particle. This section gives the words for where it sits — vortex lines, vortex tubes and their strength —
   proves that a tube cannot end inside the fluid, and then looks hard at the two simplest vortices: the rigidly
   spinning bucket and the ideal line vortex. Their pressure fields explain a bowl-shaped and a funnel-shaped surface,
   and their viscous stresses teach that irrotational does not mean stress-free.")`
2. `nb.recap("R01", "Vorticity and spin", "Vorticity is the curl of the velocity, $\boldsymbol\omega=\nabla\times\mathbf u$, and a small
   fluid element spins at ½ω (Ch. 3: the average turning rate of any two perpendicular material lines). A
   concentration of nearly parallel vorticity is a *vortex*; motion on nearly circular streamlines is *vortex motion*.",
   where="Ch. 3 §3.4")` — followed by `nb.code`: *code:* `u_sb = lambda x, t: np.array([-1.0*x[1], 1.0*x[0], 0*x[0]])`
   (a solid body turning at 1 rad/s) · `u_lv = vortices.vortex_velocity_field("line", Gamma=2*np.pi)` (the ideal line
   vortex, Γ = 2π) · `print(kinematics.vorticity(u_sb, np.array([0.3, 0.2, 0.0]))[2], kinematics.vorticity(u_lv,
   np.array([1.0, 0.5]))[-1])`. *expect:* `2.0` and ≈ `0.0` (1e-8). *explain:* "the turning tank has ω = 2 × its rotation
   rate everywhere; the line vortex has none off its axis."

**C01 — Vortex lines and tubes; a tube has one strength and cannot end (5.4)**
3. `nb.core("C01", "Vortex tubes cannot end: $\\int_V\\nabla\\cdot\\boldsymbol\\omega\\,dV=-\\Gamma_{\\text{lower end}}+\\Gamma_{\\text{upper end}}=0$ (5.4)",
   question="A smoke ring, a tornado, the swirl down a drain — each is a tube of spinning fluid. Can such a tube simply
   fade out in the middle of the air or water?")`
4. `nb.note` — **N01 [C]** "**The chapter's programme.** Think of vorticity as carried by fluid elements: an element's
   vorticity can be turned, concentrated or spread by how the element moves and deforms and by the torques of its
   neighbours. The rest of the chapter makes each of these precise — carried (C06), stretched and tilted (C10),
   twisted by baroclinic torques (C04, C09), spread by viscosity (N22)."
5. `nb.md` — **The problem in plain words:** "When a weather forecaster tracks a cyclone, an engineer a wing-tip vortex
   or an oceanographer an eddy, they follow a tube of concentrated vorticity. What rules does such a tube obey whatever
   the flow does? The first rule is purely geometric, and it is strict: the tube's strength is the same all along it —
   so it can thin, fatten, bend and twist, but it cannot stop."
6. `nb.md` — **The idea:**
   ```
   stream tube (Ch. 3)                      vortex tube (here)
   walls made of streamlines                walls made of vortex lines
   no flow crosses the wall                 no vorticity crosses the wall
   ∇·u = 0  ⇒ same volume flux Q            ∇·ω = 0 ALWAYS ⇒ same strength Γ at every section
   thinner ⇒ faster flow                    thinner ⇒ stronger spin (ω̄ = Γ/A)
   ```
   "**A vortex tube carries a fixed amount of spin, like a hose carrying water — squeeze it and the spin speeds up.**"
7. `nb.note` — **N02 [B]** "**Vortex lines** are drawn tangent to ω everywhere, exactly as streamlines (3.7) are drawn
   tangent to u:" equation `dx/\omega_x=dy/\omega_y=dz/\omega_z`, ref "5.3". "The ratio form breaks when a component is zero;
   we trace lines with the arc-length form $d\mathbf x/ds=\boldsymbol\omega/\lvert\boldsymbol\omega\rvert$ (ch03 P92, P93 reminders). In
   solid-body rotation every vertical line is a vortex line; in the line vortex only the axis is. Irrotational flow has
   no vortex lines at all."
8. `nb.code` — *code:* `u_h = ch05.helical_swirl_field(a=1.0)` (u_φ = aRz: a swirl that grows with height and radius) ·
   `om = ch05.vorticity_field(u_h)` (ω by stencils) · `s, X = ch05.vortex_line(om, np.array([1.0, 0.0, 1.0]), s_max=1.5)`
   (one vortex line through (1, 0, 1)) · `R = np.hypot(X[0], X[1])` · `print(np.ptp(X[2]*R**2))`. *expect:* spread of zR²
   along the line ≤ 1e-6 (the exact lines of this field are zR² = const). *explain:* 1. a test field whose vorticity is
   ω_R = −aR, ω_z = 2az; 2. its vorticity from central differences; 3. the arc-length ODE dx/ds = ω/\|ω\| solved by
   `solve_ivp`; 4. zR² stays constant along the traced line — the analytic answer dR/(−aR) = dz/(2az) of (5.3).
9. `nb.note` — **N03 [B]** "**Tube and strength.** The vortex lines through a closed curve form a *vortex tube*. Its
   *strength* is the circulation round a loop on the tube that encircles it once — by Stokes (2.34), the vorticity flux
   through a cross-section:" equation `\Gamma=\oint_C\mathbf u\cdot d\mathbf x=\int_A\boldsymbol\omega\cdot\mathbf n\,dA`, no ref. "Table |
   stream tube | vortex tube |: dQ = u·n dA | dΓ = ω·n dA. Number: the Gaussian (Lamb–Oseen) vortex (3.29) with
   Γ = 1 m²/s, σ_c = 0.1 m, measured on the circle r = σ_c: both routes give Γ(1 − e⁻¹) = 0.632 m²/s."
10. `nb.code` — *code:* `u_g = vortices.vortex_velocity_field("gaussian", Gamma=1.0, sigma=0.1)` · `loop =
    itg.planar_loop((0, 0, 0), (0, 0, 1), radius=0.1, n=512)`; `disc = itg.planar_disc((0, 0, 0), (0, 0, 1), radius=0.1)` ·
    `ts = ch05.vortex_tube_strength(lambda x, t=0: u_g(x, t), loop, disc)` · `print(ts)` · flipped: `loop_f =
    itg.planar_loop((0, 0, 0), (0, 0, -1), radius=0.1)`; `print(itg.circulation(u_g, loop_f))`. *expect:*
    `TubeStrength(circulation≈0.632121, flux≈0.632121, diff≤1e-6)` and −0.632121. *explain:* 1. the vortex as a
    velocity callable; 2. a circle and the disc it bounds; 3. the circulation route and the flux route of N03 agree;
    4. turning the loop over flips the sign — strength is measured with the right-hand rule.
11. `nb.derivation("D01", …)` — Part F D01 (7 steps), ref "5.4".
12. `nb.worked_example("a tube squeezed to a quarter of its area", "A thin vortex tube of strength Γ = 0.01 m²/s has a
    cross-section of 1 cm² = 10⁻⁴ m² at one place. 1. Mean vorticity there: ω̄ = Γ/A = 0.01/10⁻⁴ = 100 s⁻¹. 2. Further along
    the section is 0.25 cm² = 2.5×10⁻⁵ m². By (5.4), $-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0$, the strength is still
    0.01 m²/s. 3. So ω̄ = 0.01/2.5×10⁻⁵ = 400 s⁻¹ — four times the spin for a quarter of the area. 4. The fluid spins at
    ½ω̄: from 50 to 200 rad/s, about 8 to 32 turns per second.")`
13. `nb.code` — *code:* `for z in (0.0, 0.5, 1.0): print(z, ch05.gaussian_tube_section(z, R0=0.1, Gamma=1.0, a0=0.1,
    L=1.0))` · `print(ch05.tube_flux_budget("gaussian_tube", 0.0, 1.0, 0.1, Gamma=1.0, a0=0.1, L=1.0))` ·
    `print(ch05.tube_flux_budget("broken", 0.0, 1.0, 0.1, Gamma=1.0, a0=0.1, L=1.0))`. *expect:* flux 0.632121 at all three
    z; area 0.031416 → 0.019055 → 0.011557 m²; mean ω 20.12 → 33.17 → 54.69 s⁻¹ · `TubeFlux(lower=-0.632121, side≈0,
    upper=0.632121, total≈0)` · `TubeFlux(lower=-0.632121, side=0.0, upper=0.232544, total=-0.399577)`. *explain:* 1. the
    narrowing Gaussian tube (ω_z = (Γ/πa²)e^{−R²/a²} with a² = a₀²e^{−z/L}) — the tube through R₀ = 0.1 m keeps its
    flux while its area shrinks; 2. Gauss on the piece between z = 0 and 1 m: lower + side + upper = 0 — (5.4); 3. a field
    that is not a curl (∇·ω ≠ 0) breaks the budget: vorticity would be "created" inside.
14. `nb.check_agree` — **from scratch (curation §7):** midpoint sums of ω·n over the two lids: `om = ch05.gaussian_tube_field
    (Gamma=1.0, a0=0.1, L=1.0)`; for each lid z ∈ {0, 1}: radius `Rt = 0.1*np.exp(-z/2)`, polar midpoint grid `r =
    (np.arange(200)+0.5)/200*Rt`, `th = (np.arange(400)+0.5)/400*2*np.pi`, points `X = r cos θ, Y = r sin θ, Z = z`,
    `flux = np.sum(om(np.array([X, Y, Z]).reshape(3, -1), 0)[2] * (r*(Rt/200)*(2*np.pi/400)).ravel())` · `assert
    np.allclose(flux_top, flux_bottom, rtol=1e-4)` and `assert np.allclose(flux_top, ch05.tube_flux_budget("gaussian_tube",
    0.0, 1.0, 0.1).upper, rtol=1e-4)`. Markdown: "Our own sums over the two lids agree with the library and with each
    other."
15. `nb.plotly` — **N45 [B] (Fig. 5.1 analogue)**: plotly 3-D, two tubes side by side traced in the Burgers vortex
    (`vortices.burgers_vortex_field(Gamma=1e-3, alpha=1.0, nu=1e-6)`, lengths in mm): left, a **stream tube** — 16
    streamlines (`kinematics.streamline`) from a circle of radius 3 mm at z = 1 mm, spiralling inward and upward and
    narrowing as u_z = αz grows; right, a **vortex tube** of the narrowing Gaussian field — 16 vortex lines
    (`ch05.vortex_tube(ch05.gaussian_tube_field(), (0, 0, 0), (0, 0, 1), 0.1, n_lines=16, s_max=1.2)`), with the lids at
    z = 0 and 1 coloured blue and orange. Followed by see/read/change: *see:* "left, lines that twist and crowd together
    as they rise; right, lines that converge toward the axis." *read:* "each tube's wall is made of its own lines, so
    nothing crosses it: the stream tube carries a fixed volume flux, the vortex tube a fixed strength." *change:* "…the
    twist of the vortex field turned on (`twist=5.0`): the vortex lines become helices, the lids' fluxes do not change."
16. `nb.figure` — flux, area and mean vorticity against height z ∈ [0, 1.5] m for the Gaussian tube through R₀ = 0.1 m
    (`gaussian_tube_section`, three curves: flux teal flat, area blue falling, mean ω purple rising; log-y), plus the broken
    field's "flux" dashed rose, falling. Title "The strength is pinned; the spin pays for the squeeze". *see:* "a flat teal
    line, a falling blue curve and a rising purple one; the dashed rose line falls." *read:* "flux = area × mean vorticity
    stays constant (5.4); only a field that is not a curl lets the 'strength' change." *change:* "…L = 0.5 m (a faster
    narrowing): the area and mean-ω curves steepen twice as fast; the flux line does not move."
17. `nb.explainer("vortex_tubes_cannot_end", heading="Can a vortex just stop in the middle of the water?", why="A static
    figure shows one section; here you slide the section along a narrowing, twisting tube and watch the flux stay pinned
    while the area and the spin trade off — then break the field and watch Gauss' budget fail.", tries=["Drag the section
    from z = 0 to z = 1 m and watch the lid-flux bars: they stay equal while the mean ω grows by e ≈ 2.7.", "Switch to the
    twisted tube: the lines spiral, the wall bar stays at zero.", "Choose the broken field: the upper lid carries less than
    the lower — the badge says why this cannot be a vorticity field.", "Open the Derivation tab and step through D01 with
    the three coloured pieces."])`
18. `nb.md` — "> ⚠️ **Common confusion:** 'the side of the tube carries no flux because the vorticity is zero there.' The
    vorticity on the wall can be large; it simply points *along* the wall (ω·n = 0)."
19. `nb.md` — **What would change if…** "…the fluid were viscous? Nothing in (5.4): it is geometry, true for every flow,
    because ∇·(∇×u) = 0 always. Viscosity will change *which* lines form a tube as time goes on (C05) — but at every
    instant each tube has one strength."

**Recaps for C02 (they end C01's block and sit before the `nb.core` call of C02; one short code cell each where it helps)**
20. `nb.recap("R02", "Solid-body rotation from uniform vorticity", "A uniform plane-normal vorticity ω makes the fluid turn
    rigidly: $u_\theta=\omega r/2$ *(Eq. 5.1)* — Ch. 3's (3.22) $u_\theta=\omega_0r$ with the rotation rate ω₀ = ω/2. ⚠️ Here ω is the
    vorticity: a tank turning at Ω_tank = 5 rad/s has ω = 10 s⁻¹. (Fig. 5.2 labels the tank's rotation '2ω' — with (5.1)
    the tank turns at ω/2.)", where="Ch. 3 §3.5")` + `nb.code`: `print(ch05.solid_body_from_vorticity(0.1, 10.0),
    vortices.solid_body_rotation(0.1, 5.0))` → *expect:* `(0.5, 10.0) 0.5`.
21. `nb.recap("R03", "The ideal line vortex", "A perfect concentration of vorticity on the axis with circulation Γ gives
    $u_\theta=\Gamma/2\pi r$ *(Eq. 5.2)* — Ch. 3's $u_\theta=B/r$ with Γ = 2πB; irrotational for r > 0, and ∮u·dx = Γ on
    every circle round the axis.", where="Ch. 3 §3.5")` + `nb.code`: `print(ch05.line_vortex_gamma(0.1, 1.0),
    [vortices.circulation_circle(lambda r: ch05.line_vortex_gamma(r, 1.0)[0], r) for r in (0.1, 1.0)])` → *expect:*
    `(1.5915, 0.0) [1.0, 1.0]`.
22. `nb.recap("R04", "Paddle wheels in the two vortices", "Both vortices are steady with circular streamlines, but only
    the first turns a paddle wheel; in the line vortex a small wheel translates round the axis without spinning. The Ch. 3
    explainer `vortex_paddle_wheels` shows it (open it from the Ch. 3 page; not embedded again).", where="Ch. 3 §3.5")`
23. `nb.recap("R05", "No deformation, no viscous stress", "Solid-body rotation does not deform elements, S_ij = 0, so the
    Newtonian stress (4.37) is pure pressure, τ_ij = −pδ_ij, and Cauchy's equation (4.24) reduces to Euler's
    $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g$ (4.41).", where="Ch. 4 §4.5–4.6")` + `nb.code`: `G_sb = np.array([[0, -5.0,
    0], [5.0, 0, 0], [0, 0, 0]])` (the velocity gradient of a tank turning at 5 rad/s) · `from fluidpy.core import
    constitutive` · `print(constitutive.newtonian_stress(G_sb, p=1.0e5, mu=1e-3))` → *expect:* diag(−1e5, −1e5, −1e5), zero
    off-diagonal: pure pressure.
24. `nb.recap("R06", "Hydrostatics", "With no vertical motion the vertical momentum balance is hydrostatic:
    $0=-\partial p/\partial z-\rho g$ *(Eq. 5.5b)* — Ch. 1's (1.8).", where="Ch. 1 §1.7")`
25. `nb.recap("R07", "The Bernoulli function across streamlines", "B = ½u² + gz + p/ρ (4.69) is constant along streamlines
    in steady inviscid barotropic flow and constant everywhere only if the flow is also irrotational (4.71)–(4.72) —
    Ch. 4's `which_bernoulli` explainer.", where="Ch. 4 §4.9")`
26. `nb.recap("R08", "The Rankine vortex", "Uniform vorticity inside r ≤ a, a line vortex outside: $u_\theta=\frac{\Gamma}{2\pi a^2}r$
    inside, $\frac{\Gamma}{2\pi r}$ outside (3.28). A solid cylinder of radius a spinning at ω/2 in viscous fluid drives exactly the
    outer part, $u_\theta=\omega a^2/2r$, with Γ = πa²ω (Ch. 8 solves it).", where="Ch. 3 §3.5")` + `nb.code`:
    `print(ch05.rotating_cylinder_flow(np.array([0.05, 0.1, 0.2]), 0.1, 10.0))` → *expect:* u_θ (0.25, 0.5, 0.25) m/s,
    ω_z (10, 10, 0) s⁻¹ (Γ = π·0.01·10 = 0.314 m²/s).

**C02 — Pressure and viscous stress in the two basic vortices (5.6), (5.7)**
27. `nb.core("C02", "The bowl and the funnel: pressure in the two basic vortices (5.6), (5.7)", question="Spin a bucket
    of water and its surface becomes a bowl; pull the plug of a bath and a narrow funnel appears. Why the different
    shapes — and which of the two vortices is the viscous one?")`
28. `nb.md` — **The problem in plain words:** "Every vortex in nature, from a stirred cup to a tornado to a hurricane's
    eye wall, keeps its fluid on curved paths, and something must push that fluid inward: pressure, lower at the centre.
    How much lower depends on how the speed varies with radius. The same two basic profiles, (5.1) and (5.2), also
    settle a subtle question about viscosity: the vortex whose fluid does *not* rotate is the one with viscous stress."
29. `nb.md` — **The idea:**
    ```
    every parcel on a circle needs an inward push  ρu_θ²/r  →  pressure rises outward: ∂p/∂r = ρu_θ²/r   (5.5a)
    solid body  u_θ = ωr/2   : push grows with r   →  p ∝ r²   → a BOWL (paraboloid)            no viscous stress
    line vortex u_θ = Γ/2πr  : push ∝ 1/r³ near the axis → p ∝ −1/r² → a deep FUNNEL   viscous stress, zero net force
    ```
30. `nb.primer("partial integration with an unknown function", "When you integrate a partial derivative ∂p/∂r with respect
    to r, the 'constant' of integration may depend on the other variable z, because ∂/∂r of any f(z) is zero. So
    ∂p/∂r = 2r gives p = r² + f(z), and a second equation (for ∂p/∂z) fixes f.", code="r, z = sp.symbols('r z')
    # two variables\nf = sp.Function('f')                      # an unknown function of z alone\np = sp.integrate(2*r, r) + f(z)            # ∂p/∂r = 2r integrated in r\nprint(p, sp.diff(p, r))                      # r**2 + f(z), 2*r: the f(z) is invisible to ∂/∂r")`
    — *expect:* `r**2 + f(z) 2*r`.
31. `nb.note` — **N04 [C]** "Viscosity diffuses vorticity and lets vortex lines *reconnect* — cut and rejoin — which an
    ideal fluid forbids (C05). The book only names reconnection; it returns with turbulence (Ch. 12). Here we look at the
    two basic vortices with viscosity switched on."
32. `nb.derivation("D02", …)` — Part F D02 (10 steps), ref "5.6".
33. `nb.note` — **N05 [B]** "D02's step 5 is (5.5a): pressure supplies the centripetal acceleration —" equation
    `-\rho u_\theta^2/r=-\partial p/\partial r`, ref "5.5a". "At r = 0.1 m in a tank with ω = 10 s⁻¹ (u_θ = 0.5 m/s): ∂p/∂r =
    1000 × 0.25/0.1 = 2500 Pa/m."
34. `nb.code` — *code:* `print(ch05.solid_body_pressure_gradients(0.1, 0.0, 10.0))` · `print(ch05.solid_body_pressure(0.1,
    0.0, 10.0), ch05.isobar_height(0.1, 10.0, kind="solid"))` · `print(ch05.rotating_tank_free_surface(0.1, 0.1, 5.0))` ·
    sympy: `R, Z, w, rho, g = sp.symbols('r z omega rho g', positive=True)`; `p = rho*w**2*R**2/8 - rho*g*Z`;
    `print(sp.simplify(sp.diff(p, R) - rho*(w*R/2)**2/R), sp.simplify(sp.diff(p, Z) + rho*g))`. *expect:* `(2500.0,
    -9810.0)` · `125.0 0.012742` · `z_vertex 0.093629, z_rim 0.106371, rim_rise 0.006371, centre_drop 0.006371, r_dry 0,
    spills False` · `0 0`. *explain:* 1. the two gradients of (5.5a, b); 2. the pressure of (5.6) and the height of the
    isobar through the origin at r = 0.1 m; 3. a real tank of radius 10 cm filled 10 cm deep, spun at 5 rad/s: the water
    keeps its volume, so the centre drops by exactly as much as the rim rises; 4. sympy confirms (5.6) satisfies both
    (5.5a) and (5.5b).
35. `nb.check_agree` — **from scratch (curation §7):** integrate (5.5a) outward numerically: `r = np.linspace(0, 0.1,
    2001)`; `dpdr = 1000*(10*r/2)**2/np.where(r > 0, r, 1)`; `p_num = scipy.integrate.cumulative_trapezoid(dpdr, r,
    initial=0.0)`; `assert np.allclose(p_num, ch05.solid_body_pressure(r, 0.0, 10.0), atol=1e-6)`; the funnel from r = 1 m
    inward: `rr = np.linspace(1.0, 0.1, 4001)`; `dp = 1000*(1/(2*np.pi*rr))**2/rr`; `p_f = ch05.line_vortex_pressure(1.0,
    0.0, 1.0) + cumulative_trapezoid(dp, rr, initial=0.0)`; `assert np.allclose(p_f[-1], ch05.line_vortex_pressure(0.1,
    0.0, 1.0), rtol=1e-6)`. Markdown: "Adding up ρu_θ²/r dr ourselves gives the library's bowl and funnel."
36. `nb.note` — **R07 used: B across streamlines.** (the recap was given above; this note states the tank result) "With
    ρω²r²/8 = ρu_θ²/2, (5.6) reads $-\tfrac12u_\theta^2+gz+\frac{p(r,z)}{\rho}=\text{const.}$ — so the Bernoulli function $B=u_\theta^2/2+gz+p/\rho$
    **grows** outward, B − B(0) = ω²r²/4, as it must in a rotational flow (4.71)–(4.72)." (builder: an `nb.md`, not a
    second `nb.recap`) + `nb.code`: `print(ch05.bernoulli_across_vortex("solid", 0.1, param=10.0),
    ch05.bernoulli_across_vortex("line", 0.1, param=1.0, r_ref=1.0))` → *expect:* `0.25 0.0`.
37. `nb.note` — **N06 [B]** "**The line vortex is sheared.** Its viscous stress is" equation
    `\sigma_{r\theta}=\mu\Big[\frac1r\frac{\partial u_r}{\partial\theta}+r\frac{\partial}{\partial r}\Big(\frac{u_\theta}{r}\Big)\Big]=-\frac{\mu\Gamma}{\pi r^2}`, no ref, "nonzero
    everywhere because elements deform (Ch. 3's Fig. 3.16 picture) — yet the net viscous force on every element is zero
    (the book leaves this to Exercise 5.4; D03 below does it three ways). ⚠️ σ here is the viscous stress, not Ch. 3's
    core radius."
38. `nb.derivation("D03", …)` — Part F D03 (9 steps), no ref (the result is stated after (5.7): builder passes ref "5.7").
39. `nb.note` — **N07 [B]** "**The funnel.** The same moves as D02 with u_θ = Γ/2πr (∫ρΓ²/(4π²r³)dr = −ρΓ²/(8π²r²), and
    p → p_∞ far away at z = 0) give" equation `p(r,z)-p_\infty=-\frac{\rho\Gamma^2}{8\pi^2r^2}-\rho gz`, ref "5.7". "Isobars are
    $z=-\frac{\Gamma^2}{8\pi^2r^2g}-\frac{p-p_\infty}{\rho g}$, and $\tfrac12u_\theta^2+gz+p/\rho$ is the same everywhere (irrotational). ⚠️ The book
    calls them 'hyperboloids of revolution of the second degree'; the formula says (c − z)r² = const, a cubic surface
    — read the shape from the formula. Number: Γ = 1 m²/s, water, r = 0.1 m → deficit 1266.5 Pa and the free surface
    12.9 cm below its far level. The ideal funnel has no bottom; a real core (Rankine) gives it one: a tornado."
40. `nb.code` — *code:* `print(ch05.line_vortex_pressure(0.1, 0.0, 1.0), ch05.isobar_height(0.1, 1.0, kind="line"))` ·
    `print(ch05.line_vortex_viscous_stress(0.1, 1.0, 1e-3))` · `for kind, p in [("solid", dict(omega=10.0)), ("line",
    dict(Gamma=1.0)), ("gaussian", dict(Gamma=1.0, sigma=0.05))]: print(kind, ch05.vortex_stress_force(kind, 0.1, mu=1e-3,
    **p))` · `print(ch05.rankine_pressure(np.array([0.0, 0.05, 0.1, 0.2]), 0.0, 1.0, 0.1))`. *expect:* `-1266.515 -0.129104`
    · `-0.031831` · solid: σ 0, forces 0, 0, 0 · line: σ −0.031831, forces 0, 0, 0 (≤ 1e-12) · gaussian: σ ≠ 0, the three
    force routes agree and are nonzero · Rankine p: −2533.03, −2216.40, −1266.51, −316.63 Pa. *explain:* 1. (5.7) and its
    isobar; 2. the stress of N06; 3. the N09 table in numbers — three routes to the net viscous force for three vortices;
    4. the Rankine composite: continuous pressure, twice as deep at the centre as at the core edge.
41. `nb.note` — **N08 [B]** "**Torque to infinity.** Round a spinning cylinder of radius a, the viscous torque per unit
    length transmitted across every circle of radius r is the same:" equation `T'=2\pi r^2\sigma_{r\theta}=-2\mu\Gamma`, no ref.
    "The power it carries, |T′| × (u_θ/r) = μΓ²/πr², falls outward; the difference between two radii is dissipated as
    heat: $\int_a^R2\pi r\,\rho\varepsilon\,dr=\frac{\mu\Gamma^2}{\pi}\big(\frac1{a^2}-\frac1{R^2}\big)$. As R → ∞ all the power the
    cylinder puts in is dissipated — and a container at any radius would feel the full torque. ⚠️ ε here is the
    dissipation rate (Ch. 4 (4.58)), not the Levi-Civita symbol. Number: a = 0.1 m, R = 1 m, Γ = 1 m²/s, water: in
    0.031831 W/m, out 0.000318 W/m, dissipated 0.031513 W/m."
42. `nb.code` — *code:* `print([ch05.torque_per_length(r, 1.0, 1e-3) for r in (0.1, 1.0, 10.0)])` · `d =
    ch05.dissipation_outside_cylinder(0.1, 1.0, 1.0, 1e-3); print(d, d["power_in"] - d["power_out"] - d["dissipated"])`.
    *expect:* `[-0.002, -0.002, -0.002]` · dissipated 0.0315127, power_in 0.0318310, power_out 3.1831e-4, balance ≤ 1e-12.
43. `nb.note` — **N09 [B]** "**The principle.** Irrotationality does not mean *no viscous stress*; it means *no net
    viscous force*: for incompressible flow $\mu\nabla^2\mathbf u=-\mu\nabla\times\boldsymbol\omega$ (4.40) vanishes where ω = 0. Solid-body
    rotation is the only motion with no viscous stress at all." Table (from cell 40) | vortex | σ_rθ | net viscous force |:
    solid body | 0 | 0 · line vortex | −μΓ/πr² ≠ 0 | 0 · Gaussian (Lamb–Oseen) | ≠ 0 | ≠ 0 (it decays, Ch. 8).
44. `nb.worked_example("a spinning bucket and a bathtub drain", "1. **Bucket**: radius 10 cm, spun at 5 rad/s, so
    ω = 10 s⁻¹ (⚠️ not 5). Rim–centre height of the free surface: ω²R²/8g = 100 × 0.01/(8 × 9.81) = 0.0127 m = 12.7 mm.
    Volume is conserved, so the centre drops 6.4 mm and the rim rises 6.4 mm. 2. **Drain**: Γ = 1 m²/s (strong). At
    r = 10 cm: u_θ = Γ/2πr = 1.59 m/s; pressure deficit ρΓ²/8π²r² = 1000/(8 × 9.87 × 0.01) = 1266 Pa; the surface there
    lies 1266/(1000 × 9.81) = 0.129 m below the far level. At r = 5 cm the depth quadruples to 0.52 m — a funnel. 3.
    **Tornado** (air, ρ = 1.2, Γ = 10⁴ m²/s, core a = 50 m): peak wind Γ/2πa = 31.8 m/s, pressure deficit 608 Pa at the
    core edge and 1216 Pa at the centre (Rankine) — about 1 % of the atmospheric pressure.")`
45. `nb.plotly` — **N46, N47 [B] (Figs. 5.2 and 5.3 analogues)** `slider_figure` over the tank rate Ω_tank = 0…10 rad/s
    (21 steps, FAST 11): (r, z) half-plane r ∈ [−0.15, 0.15] m; traces "free surface (5.6)" (teal), "isobars p = 500,
    1000 Pa (5.6)" (orange dashed), "funnel free surface (5.7), Γ = 0.05 m²/s" (rose, fixed), "Rankine surface a = 1 cm"
    (rose dashed, fixed). Title "A spinning bucket is a bowl; a drain is a funnel". *see:* "the teal bowl deepens as the
    slider moves; the rose funnel plunges near the axis, the Rankine curve gives it a flat bottom." *read:* "both shapes
    come from ∂p/∂r = ρu_θ²/r (5.5a); the bowl because u_θ grows with r, the funnel because it grows toward the axis."
    *change:* "…the tank rate doubled: the bowl's depth quadruples (ω² in (5.6)), its shape stays a paraboloid."
46. `nb.figure` — two panels on shared r axes (0–0.3 m): left, u_θ(r) teal for solid body (ω = 10), line vortex (Γ = 1),
    Rankine (a = 0.1); right, σ_rθ(r) rose and the net viscous force (black, zero for the first two, nonzero for the
    Gaussian of σ_c = 0.05 m). *see / read / change* written from N09: "the line vortex has the largest stress near the
    axis and zero force; the Gaussian has both." *change:* "…μ doubled: every stress and force doubles; the zeros stay
    zeros."
47. `nb.explainer("vortex_pressure_funnel", heading="Why is a spinning bucket a bowl and a drain a funnel?", why="The
    bowl, the funnel, the Rankine tornado and the rotating cylinder share one radial balance; dragging the rotation, the
    circulation and the core size while a probe element shows its pressure, stress and net force makes 'stress without
    force' visible.", tries=["Preset 'bucket 5 rad/s': read the rim–centre height (12.7 mm) and the zero stress bar.",
    "Switch to the line vortex: the stress bar fills, the net-force bar stays at zero.", "Try 'tornado': the same
    formula in air, 600 Pa at the core edge.", "Open Derivation → D03 step 6 and watch the two faces of the probe
    element."])`
48. `nb.md` — **What would change if…** "…the fluid were spinning but the frame of reference turned with it? The water
    would be at rest in that frame, yet the free surface would still be a paraboloid: in the rotating frame the
    centrifugal force (Ch. 4 §4.7) does the job that the centripetal acceleration does here. That is the laboratory
    analogue of geostrophic balance (Ch. 13)."

### A.2 §5.2 Kelvin's Circulation Theorem — R09, C03 (+N10–N14 N16 N48 · D04 D05 · E3), C04 (+N15 N49 N50 · D06 D07 · E4)
1. `nb.section("5.2", "Kelvin's Circulation Theorem", intro="**What is this section about?** Circulation round a loop
   that moves with the fluid is conserved — if the fluid is inviscid, barotropic, pushed only by conservative forces and
   watched from a non-rotating frame. The proof shows exactly which term wakes up when each hypothesis fails: those are
   the only ways vorticity can be created. The last of them, baroclinicity, has a simple mechanical picture: a torque.")`
2. `nb.recap("R09", "Barotropic fluid and the pressure function", "A fluid is barotropic when its density depends on
   pressure alone, ρ = ρ(p) (constant density, isothermal or isentropic gas). Then dp/ρ is the exact differential of
   the pressure function $\mathcal P(p)=\int_{p_o}^{p}dp'/\rho(p')$ (4.67), just as g·dx = −dΦ is exact for a conservative
   body force (4.18).", where="Ch. 4 §4.9")` + `nb.code`: `from fluidpy.core import bernoulli`;
   `print(bernoulli.pressure_function(2.0e5, 1.0e5, kind="isothermal", T=288.15))` → *expect:* ≈ 57 330 J/kg (RT ln 2).

**C03 — Kelvin's circulation theorem (5.8)**
3. `nb.core("C03", "Kelvin's circulation theorem: $D\\Gamma/Dt=0$ (5.8)", question="Put a ring of dye into water and let
   the flow stretch it seven times longer and wind it into a spiral. What has stayed exactly the same?")`
4. `nb.md` — **The problem in plain words:** "Circulation measures the net swirl inside a loop. If we could prove it
   never changes for loops that move with the fluid, we would know that a fluid starting from rest (zero circulation
   everywhere) can never start swirling on its own — the basis of all 'irrotational' flow in Ch. 6, of how wings get
   lift, and, on a rotating planet, of how cyclones spin up (C11). The proof also tells us precisely how swirl *can* be
   made."
5. `nb.md` — **The idea:**
   ```
   Γ(t) = ∮ u·dx round a loop of dye         d/dt:  (acceleration along the loop)  +  (stretching of the loop)
                                                        = ∮ forces·dx / ρ                 = ∮ d(½u²) = 0 always
   forces: pressure −∇p/ρ  → −∮dp/ρ  = 0 if ρ = ρ(p)         (barotropic)
           gravity −∇Φ     → −∮dΦ    = 0 if Φ exists          (conservative)
           viscous         → ∮...    = 0 if inviscid          (ideal fluid)
           Coriolis        → ≠ 0 in a rotating frame          (inertial frame needed)
   ```
6. `nb.primer("closed-loop integral of an exact differential", "If F is a single-valued function, going once round a
   closed loop brings F back to its starting value, so ∮dF = 0. It fails for a multi-valued 'function' such as the polar
   angle θ, which grows by 2π on a loop round the origin. Kelvin's proof uses this twice: for ½u² and for the pressure
   function.", code="s = np.linspace(0, 2*np.pi, 2001)                  # a loop: the unit circle\nx, y = np.cos(s), np.sin(s)\nF = x**2*y + 3*x                                   # single-valued F(x, y)\nprint(np.sum(np.diff(F)))                           # ∮dF = F(end) − F(start) = 0\ntheta = np.unwrap(np.arctan2(y, x))                  # the polar angle, followed continuously\nprint(theta[-1] - theta[0])                          # 2π: θ is not single-valued round the origin")`
   — *expect:* `0.0` (round-off) and `6.283185`.
7. `nb.primer("periodic trapezoid rule on a closed loop", "For a smooth closed loop parametrised by s ∈ [0, 1) with N
   equally spaced points, the plain average Σf(s_k)/N is extraordinarily accurate — the errors of the trapezoid rule
   cancel round a periodic curve (error falls faster than any power of 1/N). We compute ∮u·dx this way, with dx/ds from
   an FFT derivative (numpy's `np.fft` splits the periodic coordinates into waves e^{2πiks}; differentiating a wave just
   multiplies it by 2πik — C07's primer says more).", code="for N in (8, 16, 32):                                   # three loop resolutions\n    s = np.arange(N)/N                                  # equally spaced tags, end point not repeated\n    x, y = np.cos(2*np.pi*s), np.sin(2*np.pi*s)          # the unit circle\n    dxds, dyds = -2*np.pi*y, 2*np.pi*x                    # tangent dx/ds\n    u, v = -y, x                                          # solid body: Γ = 2 × area = 2π\n    print(N, np.mean(u*dxds + v*dyds) - 2*np.pi)          # error ~ 1e-15 already at N = 8")`
   — *expect:* errors ≤ 1e-14 for every N.
8. `nb.primer("advecting many points in one solve_ivp call", "To move a whole loop of N particles we give `solve_ivp`
   one long state vector (all x's, then all y's) and a right-hand side that reshapes it to (2, N), evaluates the velocity
   at every point at once (vectorised) and flattens it back. One adaptive integrator then moves the whole loop
   consistently.", code="from scipy.integrate import solve_ivp\nX0 = np.array([[1.0, 0.0, -1.0], [0.0, 1.0, 0.0]])    # three particles, shape (2, 3)\nrhs = lambda t, y: np.array([-y.reshape(2, -1)[1], y.reshape(2, -1)[0]]).ravel()   # solid body u = (−y, x)\nsol = solve_ivp(rhs, (0, np.pi/2), X0.ravel(), rtol=1e-10)\nprint(sol.y[:, -1].reshape(2, -1).round(6))             # each point turned by 90°")`
   — *expect:* `[[0, -1, 0], [1, 0, -1]]` (to 1e-6).
9. `nb.md` — gloss: "**Material label.** Tag each particle of the loop once with a number s ∈ [0, 1); the tags travel
   with the particles, so the loop is always 'the particles with tags 0…1' — fixed limits, even though the loop moves
   (Ch. 3 P89)."
10. `nb.derivation("D04", …)` — Part F D04 (8 steps), ref "5.9".
11. `nb.note` — **N10 [B]** "D04's step 4 is (5.9):" equation
    `\frac{D\Gamma}{Dt}=\frac{D}{Dt}\oint_Cu_i\,dx_i=\oint_C\frac{Du_i}{Dt}dx_i+\oint_Cu_i\frac{D}{Dt}(dx_i)`, ref "5.9". "Two causes of change: the
    fluid accelerating along the loop, and the loop's elements changing."
12. `nb.note` — **N11 [B]** "**The element's rate is the velocity difference** (Fig. 5.4): $\mathbf u+d\mathbf u=\frac{D}{Dt}(\mathbf x+d\mathbf x)$,
    so $D(d\mathbf x)/Dt=d\mathbf u$, and the second integral of (5.9) is $\oint_Cu_i\,du_i=\oint_Cd(\tfrac12u_i^2)=0$." equation
    `\frac{D}{Dt}(d\mathbf x)=d\mathbf u`, no ref.
13. `nb.derivation("D05", …)` — Part F D05 (9 steps), ref "5.11".
14. `nb.note` — **N12 [B]** "D05's step 4 is (5.10):" equation
    `\oint_C\frac{Du_i}{Dt}dx_i=-\oint_C\frac1\rho dp-\oint_Cd\Phi+\oint_C\Big(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j}\Big)dx_i`, ref "5.10".
15. `nb.note` — **N13 [B]** "What survives is (5.11):" equation
    `\frac{D\Gamma}{Dt}=\oint_C\Big(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j}\Big)dx_i`, ref "5.11". "Kelvin holds if the fluid is
    inviscid or the net viscous force vanishes along C. Number, the Lamb–Oseen vortex (Γ₀ = 0.01 m²/s, ν = 10⁻⁶ m²/s) on the
    circle r = 5 mm at t = 10 s: Γ = Γ₀(1 − e^{−r²/4νt}) = 0.00465 m²/s and falling at 3.35×10⁻⁴ m²/s² — exactly the
    (5.11) integral. ⚠️ The book adds 'this occurs when C lies entirely in irrotational fluid' — true for incompressible,
    constant-μ flow, where the viscous force is −μ∇×ω (4.40); a compressible flow keeps a ∇(∇·u) force."
16. `nb.worked_example("circulation that survives a stretching — and one that leaks away", "1. **Cellular flow**
    ψ = sin x sin y (u = sin x cos y, v = −cos x sin y, in m and m/s): a steady inviscid flow whose vorticity
    ω = 2 sin x sin y is constant on streamlines. A circle of radius 0.5 m round (π/2, 0.9 m) has Γ(0) = ∬2 sin x sin y dA
    = 1.155 m²/s. Kelvin predicts Γ(t) = 1.155 m²/s forever. 2. **Lamb–Oseen** (viscous): the circle r = 5 mm with
    Γ₀ = 0.01 m²/s, ν = 10⁻⁶ m²/s at t = 10 s: r²/4νt = 25×10⁻⁶/(4×10⁻⁵) = 0.625; Γ = 0.01(1 − e^{−0.625}) = 0.01 ×
    0.4647 = 0.00465 m²/s; at t = 20 s: r²/4νt = 0.3125, Γ = 0.00268 m²/s — the vorticity diffuses out through the loop.")`
17. `nb.code` — *code:* `sc = ch05.kelvin_scenario("cellular")` · `t = np.linspace(0, sc["t_end"], 7)` · `G =
    ch05.material_circulation(sc["u"], sc["pts0"], t)` · `L = [ch05.loop_length(P) for P in ch05.material_loop(sc["u"],
    sc["pts0"], t)]` · `print(np.round(G, 10), np.round(L, 3))` · `print(ch05.kelvin_rate_terms(sc["u"], sc["pts0"], 0.0))` ·
    `print(ch05.lamb_oseen_circulation(0.005, 10.0, 0.01, 1e-6), ch05.lamb_oseen_viscous_loop_integral(0.005, 10.0, 0.01,
    1e-6))` · `print(ch05.kelvin_scenario_circulation("rotating", 5.0))`. *expect:* Γ = 1.155130 at every time (spread ≤
    1e-9) while the length grows from 3.14 m to 21.0 m (×6.7) · `KelvinRate(acceleration≈0, contour≈1e-13, total≈0)` ·
    `(0.00464739, -3.34538e-04) -3.34538e-04` · `1.985865`. *explain:* 1. the scenario's velocity and initial loop (256
    points); 2. the loop advected by `solve_ivp` in one call (primer) and its circulation by the periodic trapezoid rule;
    3. its length, to show it really is stretched; 4. the two terms of (5.9); 5. the Lamb–Oseen circulation and its rate
    equal the (5.11) viscous integral; 6. in a rotating frame the relative circulation grows (C11 explains).
18. `nb.check_agree` — **from scratch (curation §7):** RK4 by hand (ch03 P95) on the 256 loop points of the cellular flow,
    `dt = 0.01`, to t = 6 s; `Gamma = np.mean(np.sum(u(P)*dPds, axis=0))` with `dPds` from `np.fft` differentiation of the
    periodic coordinates; `assert np.allclose(Gamma, ch05.material_circulation(sc["u"], sc["pts0"], [0, 6.0])[-1],
    rtol=1e-7)` and `assert abs(Gamma - 1.155130) < 1e-6`. Markdown: "Our own particles and our own loop sum agree: the
    circulation did not move."
19. `nb.animation` — **N48 [B] (Fig. 5.4 analogue)** (`player="video"`, 60 frames, FAST 30): left, the cellular flow's
    streamlines (grey) with the material loop (teal, 256 points) stretched into a spiral over 6 turnover times (37.7 s), one
    element dx (purple arrow) with u and u + du at its ends; right, Γ(t) (teal, flat) with the Lamb–Oseen Γ(t) of a 5 mm
    circle (rose, falling) and the dashed ghost of Γ(0). *see:* "the loop winds into a thin spiral; its curve on the right
    does not move." *read:* "however long the loop gets, the acceleration along it integrates to zero (D05): Kelvin. The
    rose curve is the viscous leak of (5.11)." *change:* "…ν switched on in the cellular flow: the teal curve would start
    to fall like the rose one, fastest where the spiral is thinnest."
20. `nb.note` — **N14 [B]** "**Three ways to make or destroy vorticity** — one per broken hypothesis." Table | source |
    term of (5.10) | example | number |: non-conservative body force | ∮g·dx ≠ 0 (Coriolis) | a drain vortex in a tank on
    the rotating Earth | the rotating scenario gains 1.99 m²/s in 5 s · non-barotropic (baroclinic) | −∮dp/ρ ≠ 0 | lock
    exchange (C04), sea breeze | the lock-exchange loop's initial rate `ch05.kelvin_scenario_rate("baroclinic")["pressure"]`
    = 0.0467 m²/s² (the baroclinic source ∇ρ×∇p/ρ² integrated over the 0.2 m × 0.2 m square, by Stokes) · net
    viscous force | ∮(1/ρ)∂σ_ij/∂x_j dx_i ≠ 0 | boundary layers at walls (Ch. 9) | −3.35×10⁻⁴ m²/s² (Lamb–Oseen).
21. `nb.note` — **N16 [B]** "**Four restrictions** keep irrotational flow irrotational: no net viscous force along C
    (boundary layers diffuse vorticity in), conservative body forces (they act through the centre of mass), a barotropic
    fluid (else baroclinic: temperature, salinity, composition), an inertial frame (§4.7's extra terms)." + `nb.code`:
    `import itertools; for combo in itertools.product([True, False], repeat=4): r = ch05.kelvin_hypotheses(*combo);
    print(combo, r["holds"], r["surviving_terms"])` → *expect:* 16 rows; only (True, True, True, True) holds; each False
    adds exactly its term ("viscous", "baroclinic", "body", "coriolis").
22. `nb.md` — "> ⚠️ **Common confusion:** 'circulation is conserved round any loop.' Only round a **material** loop. A
    loop fixed in space in the same inviscid flow sees Γ change as different fluid (with different vorticity) passes
    through it — E3's 'fixed loop' mode shows it."
23. `nb.explainer("kelvin_material_loop", heading="Stretch and tangle a loop of dye — what stays the same?", why="Kelvin's
    theorem is about time: watching the loop tangle while Γ(t) stays flat, then breaking one hypothesis at a time and
    seeing which term of (5.10) wakes up, teaches both the theorem and its limits; its Helmholtz mode shows C05.",
    tries=["Press ▶ on the cellular flow: the loop grows ×6.7, Γ stays at 1.155 m²/s.", "Switch to 'fixed loop': Γ
    now changes — the loop is not material.", "Choose the baroclinic flow: the orange ∮dp/ρ curve lifts off zero.",
    "Try the rotating frame: Γ grows while Γ_a stays at π — a preview of C11."])`
24. `nb.md` — **What would change if…** "…the fluid were stratified (density depending on temperature as well as
    pressure)? Then ∮dp/ρ need not vanish: where surfaces of constant density cross surfaces of constant pressure,
    circulation — and vorticity — is born. C04 shows the torque that does it."

**C04 — Barotropic vs baroclinic: the pressure torque (Fig. 5.6)**
25. `nb.core("C04", "Baroclinic torque: crossed isobars and isopycnals spin fluid up", question="On a sunny afternoon
    the land warms, the sea stays cool, and a sea breeze starts to blow. Why does a density difference *sideways* make air
    turn over, while a density difference up–down leaves it at rest?")`
26. `nb.md` — **The problem in plain words:** "Warm and cold air side by side, fresh river water meeting salty sea water,
    a front between two air masses: in each, density changes horizontally while gravity pulls vertically. Kelvin's
    theorem says the barotropic hypothesis is broken there, so circulation must change. But *why*, mechanically? A small
    blob of fluid gives the answer: pressure pushes through its middle, its weight hangs off-centre, and it is spun."
27. `nb.md` — **The idea:**
    ```
    barotropic element                 baroclinic element
    isobars ∥ isopycnals               isobars cross isopycnals
    pressure force through centre O    pressure force through centre O   (a circle: every push points at O)
    centre of mass G = O (by symmetry  G shifted toward the heavy side
      along the isobars)
    no lever arm → no torque           lever arm O–G → torque ∝ ∇ρ × ∇p  →  spin-up  Dω/Dt = ∇ρ×∇p/ρ²
    ```
28. `nb.md` — glosses: "**Centre of mass** of a body with varying density: x_G = ∫ρx dA/∫ρ dA — pulled toward the heavy
    side. **Disc's moment of inertia** about its centre, I = ½MR² (P116 gave the cube's). **Smoothed step**:
    ρ(x) = ρ̄ − (Δρ/2)tanh(2x/δ) goes from ρ̄ + Δρ/2 to ρ̄ − Δρ/2 across a layer of width ≈ δ, with slope −Δρ/δ at its
    centre."
29. `nb.derivation("D06", …)` — Part F D06 (10 steps), ref "5.28".
30. `nb.worked_example("a 1 cm blob in water with a sideways density gradient", "Disc radius R = 1 cm; ρ₀ = 1000 kg/m³
    increasing to the right at 10 kg/m⁴ (∇ρ = (10, 0)); hydrostatic ∇p = (0, −ρ₀g) = (0, −9810) Pa/m. 1. Centre of mass
    offset: R²|∇ρ|/4ρ₀ = 10⁻⁴ × 10/4000 = 2.5×10⁻⁷ m to the right. 2. Pressure force: πR²|∇p| = 3.1416×10⁻⁴ × 9810 =
    3.08 N/m upward, through the centre. 3. Torque about G: 2.5×10⁻⁷ × 3.08 = 7.70×10⁻⁷ N m/m, clockwise (the upward push
    acts left of G). 4. I_G = ρ₀πR⁴/2 = 1.571×10⁻⁵ kg m. 5. Spin-up: 2 × 7.70×10⁻⁷/1.571×10⁻⁵ = 0.0981 s⁻², clockwise —
    the heavy right side sinks. 6. Formula: (∇ρ×∇p)_z/ρ₀² = 10 × (−9810)/10⁶ = −0.0981 s⁻² ✓.")`
31. `nb.code` — *code:* `r = ch05.pressure_torque_on_element()` (defaults = the example) · `print({k: r[k] for k in
    ("x_G", "torque", "I_G", "spin_up", "baroclinic", "ratio")})` · `print(ch05.pressure_torque_on_element(grad_rho=(0.0,
    -10.0))["spin_up"])` (stable stratification: ∇ρ ∥ ∇p) · `R = np.array([0.1, 0.03, 0.01]); err = [abs(ch05.
    pressure_torque_on_element(radius=Ri, grad_rho=(100.0, 0.0))["ratio"] - 1) for Ri in R]; print(err,
    observed_order(R, err))`. *expect:* x_G (2.5e-7, 0), torque −7.70e-7, I_G 1.5708e-5, spin_up −0.0981, baroclinic −0.0981,
    ratio 1.0000000013 · `0.0` · errors ≈ (1.25e-5, 1.1e-6, 1.25e-7), order 2.0. *explain:* 1. the disc of the example:
    the torque route and the formula agree; 2. density increasing downward (stable, barotropic-like here): no torque; 3.
    the small difference between the two routes falls as R² — in the limit of a point they are the same thing.
32. `nb.check_agree` — **from scratch (curation §7):** 4000 rim points: `th = (np.arange(4000)+0.5)/4000*2*np.pi`; `xb, yb
    = R*np.cos(th), R*np.sin(th)`; `p = 1e5 + 0*xb - 9810*yb`; `F = [np.sum(-p*np.cos(th))*ds, np.sum(-p*np.sin(th))*ds]`;
    the area sums for M, x_G, I_G on a 400 × 800 polar grid; `tau = np.sum(((xb - xG)*(-p*np.sin(th)) - (yb - yG)*(-p*
    np.cos(th)))*ds)`; `assert np.allclose(2*tau/IG, r["spin_up"], rtol=1e-8)`. Markdown: "Summing the pushes ourselves
    gives the same spin-up."
33. `nb.figure` — **N50 [B] (Fig. 5.6 analogue)**: two discs side by side (`scripts.ch05_drawings.disc_element`): left
    barotropic (isobars and isopycnals both horizontal), right baroclinic (isopycnals tilted 45°); inward pressure arrows
    round the rim (orange, ∝ p − p_min), isobars orange solid, isopycnals blue dashed, the geometric centre O, the centre
    of mass G (offset exaggerated ×10⁴ and labelled so), the net-force line through O, a curved torque arrow on the right
    disc. *see:* "on the left, force line and G coincide; on the right, the force line misses G." *read:* "the pressure
    force always passes through O; only a density gradient across the isobars moves G off that line — that lever arm is
    the torque." *change:* "…the isopycnals turned to vertical (90° to the isobars): the offset, the lever arm and the
    torque are largest (sin 90° = 1)."
34. `nb.plotly` — `slider_figure` over the tilt angle between ∇ρ and ∇p, 0…180° (19 steps, FAST 10): spin-up rate vs angle
    (sine curve, orange) with the current angle's dot and the formula ghost; title "The torque switches on as the isolines
    cross". *see / read / change:* "zero at 0° and 180° (barotropic), largest at 90°; the sign flips past 180° (the other
    sense of rotation)."
35. `nb.note` — **N15 [B]** "**Lock exchange** (Fig. 5.5): fresh water (ρ₁) and salt water (ρ₂) side by side, a gate
    pulled out: the heavy fluid slumps under, the light rides over, the interface tilts — vorticity has been made. Its
    first-instant rate (derived below; the book leaves it to Exercise 5.5):" equation
    `\frac{D\omega_z}{Dt}=\frac{2(\rho_2-\rho_1)g}{(\rho_2+\rho_1)\delta}`, no ref. "Number: ρ₁ = 1000, ρ₂ = 1025 kg/m³, interface 10 cm thick:
    2.42 s⁻². Sea breezes, gravity currents and oceanic fronts start the same way (Ch. 13)."
36. `nb.derivation("D07", …)` — Part F D07 (7 steps), no ref (builder passes ref "5.28" — the baroclinic term it evaluates).
37. `nb.code` — *code:* `print(ch05.lock_exchange_initial_vorticity_rate(1000.0, 1025.0, 0.1))` · `F =
    ch05.lock_exchange_fields(1000.0, 1025.0, 0.1)` · `print(ch05.baroclinic_term(F["rho"], F["p"], np.array([0.0, 0.5]))
    [2])` · `print([ch05.lock_exchange_initial_vorticity_rate(1000.0, 1025.0, d) for d in (0.2, 0.05, 0.01)])`. *expect:*
    `2.422222` · `2.4222` (to 1e-6; the field route) · `[1.2111, 4.8444, 24.222]`. *explain:* 1. the closed form of D07;
    2. the field route: `baroclinic_term` differentiates the tanh-step density and the hydrostatic pressure at the
    interface — the same number; 3. a thinner interface spins faster: in the limit it is a vortex sheet (C14).
38. `nb.figure` — **N49 [B] (Fig. 5.5 analogue)**: left, the tank before (two colours, gate dashed); right, the colour map
    of (∇ρ × ∇p)_z/ρ² on the smoothed interface (orange, from `baroclinic_term` on a 201 × 101 grid) with counterclockwise
    tumbling arrows. *see:* "a bright vertical strip at the interface." *read:* "vorticity is created only where density
    changes sideways — at the interface — in the counterclockwise sense: heavy under, light over." *change:* "…the gate
    between two fluids of equal density: the strip vanishes; nothing happens."
39. `nb.explainer("baroclinic_torque", heading="How can density make fluid spin?", why="Turning the isopycnals against the
    isobars and watching the force line leave the centre of mass, the torque grow as the sine of the angle and two
    independent routes agree, makes (5.28) a mechanism instead of a formula; the lock-exchange mode shows it on a whole
    interface.", tries=["Set the tilt to 0°: no torque however large ∇ρ.", "Drag the tilt to 90°: the badge reads the
    spin-up and its sense.", "Shrink the disc: the two routes converge (the ratio → 1).", "Lock exchange: thin the
    interface and watch the rate grow as 1/δ."])`
40. `nb.md` — **What would change if…** "…the fluid were also rotating, like the atmosphere? The baroclinic torque would
    still act, but the new spin would be turned and stretched by the planet's own vorticity — both appear together in
    the rotating-frame vorticity equation (C09)."

### A.3 §5.3 Helmholtz's Vortex Theorems — C05 (+N17 N51 · D08; E3's Helmholtz mode)
1. `nb.section("5.3", "Helmholtz's Vortex Theorems", intro="**What is this section about?** Under Kelvin's four
   restrictions, vorticity is glued to the fluid: vortex lines move with it, and a tube keeps its strength along its
   length and in time. Two of the four theorems are (5.4) again; the other two follow from Kelvin in a few lines.")`

**C05 — Helmholtz's vortex theorems**
2. `nb.core("C05", "Helmholtz's theorems: vortex lines move with the fluid", question="A smoke ring drifts across a room.
   Is the ring of spin always made of the same smoke — or does the swirl slide through the air?")` (the loop-vector-area
   primer below is this chapter's first use of A_vec; C11 reuses it)
3. `nb.md` — **The problem in plain words:** "If vorticity stays attached to the same fluid, we can follow a vortex
   simply by following its fluid: the whole of §5.7 (point vortices, rings) and the potential-vorticity thinking of
   Ch. 13 rest on it. Helmholtz found four rules; we prove the one that needs work."
4. `nb.md` — **The idea:**
   ```
   (1) vortex lines move with the fluid       ← Kelvin on patches of a tube wall (D08)
   (2) a tube has one strength along it       ← (5.4), geometry (D01)
   (3) a tube cannot end in the fluid          ← (5.4)
   (4) a tube keeps its strength in time        ← Kelvin on a loop round the tube
   (1), (4) need Kelvin's four restrictions; (2), (3) hold always
   ```
5. `nb.note` — **N51 [C] (Fig. 5.7)** "Our 3-D picture below shows a patch S on the wall of a tube and its material image
   S′ some time later; we make no separate computation of it."
6. `nb.derivation("D08", …)` — Part F D08 (9 steps), no ref (builder passes ref "5.8").
7. `nb.note` — **N17 [B]** "**Why 'for every patch' matters.** One patch with zero circulation proves only that the total
   flux through it vanishes; because *every* patch on the wall keeps zero circulation, the flux vanishes locally (D08
   step 6) — the carried wall is again a tube wall. A second proof (the field-equation route: ω and a material element
   obey the same equation when ν = 0) needs the vorticity equation and is checked at the end of C06."
8. `nb.primer("loop vector area ½∮x × dx", "For a closed loop, A_vec = ½∮x × dx is a vector whose length is the area
   enclosed (for a flat loop) and whose direction is the loop's normal by the right-hand rule; for a curved loop it is
   ∫n dA over any surface spanning it. The 2-D version is the surveyor's shoelace formula. We use it to ask whether a
   small loop still lies on a tube wall (its normal ⟂ ω) and, in C11, for the planetary circulation.", code="s =
   2*np.pi*np.arange(400)/400                          # tags round a loop\nP = np.array([2*np.cos(s), np.sin(s), 0*s])                   # an ellipse, semi-axes 2 and 1\ndP = np.roll(P, -1, axis=1) - P                                # the chord to the next point (np.roll shifts by one)\nprint(0.5*np.sum(np.cross(P.T, dP.T), axis=0), 2*np.pi)         # ≈ (0, 0, 6.28293) vs πab = 6.28319")`
   — *expect:* `[0. 0. 6.28293] 6.28319` (400 chords cut the corners by 4×10⁻⁵); the library's spectral
   `ch05.loop_vector_area(P)` → (0, 0, 6.283185).
9. `nb.worked_example("a small loop painted on a tube wall", "The ABC flow u = (sin z + cos y, sin x + cos z, sin y +
   cos x) (in m and m/s) is a steady solution of Euler's equation whose vorticity equals its velocity, ω = u — Kelvin
   holds exactly. 1. At the origin ω = (1, 1, 1) s⁻¹, \|ω\| = 1.732 s⁻¹. 2. Paint a circle of radius 1 cm there in the plane
   spanned by ω and (1, −1, 0)/√2: its normal (1, 1, −2)/√6 is perpendicular to ω, so the patch lies on the local tube
   wall. 3. Its flux ω·n A is zero to first order; the curvature of the field leaves Γ ≈ 8×10⁻¹⁰ m²/s, against
   \|ω\|A = 5.4×10⁻⁴ m²/s for a patch facing ω. 4. Let the flow carry it for 2 s: Kelvin keeps Γ at 8×10⁻¹⁰, so the carried
   patch still faces across ω — it is still on a tube wall.")`
10. `nb.code` — *code:* `sc = ch05.kelvin_scenario("helmholtz_abc")` (the ABC flow and the painted loop) · `t = np.linspace(0,
   2.0, 5)` · `P = ch05.material_loop(sc["u"], sc["pts0"], t)` · `G = ch05.material_circulation(sc["u"], sc["pts0"], t)` ·
   `cosang = [ch05.loop_vector_area(Pk) @ sc["u"](Pk.mean(axis=1), 0.0) / (np.linalg.norm(ch05.loop_vector_area(Pk)) *
   np.linalg.norm(sc["u"](Pk.mean(axis=1), 0.0))) for Pk in P]` · `print(G, np.max(np.abs(cosang)))`. *expect:* Γ =
   8.0159e-10 m²/s at every time (spread ≤ 1e-15), max \|cos\| ≤ 2×10⁻⁵ (the loop's normal stays at 90° to ω) while the
   loop's length changes by ≈ 10 %. *explain:* 1. an inviscid steady flow and a small loop lying on a tube wall; 2. the
   loop advected; 3. its circulation stays at its tiny starting value (Kelvin); 4. the angle between the loop's vector
   area and ω stays 90°: the carried patch is still on a tube wall — Helmholtz 1 seen in numbers.
11. `nb.check_agree` — **from scratch (curation §7):** the vector area by the shoelace chords of the primer, `A = 0.5*np.sum(
    np.cross(Pk.T, (np.roll(Pk, -1, axis=1) - Pk).T), axis=0)`, for the last loop; `assert np.allclose(A,
    ch05.loop_vector_area(Pk), rtol=1e-3)`; the flux through it `A @ sc["u"](Pk.mean(axis=1), 0.0)` ≤ 1e-8 m²/s. Markdown:
    "Our own chord sum says the same: no vorticity pierces the carried patch."
12. `nb.plotly` — plotly 3-D schematic (N51): the narrowing Gaussian tube (translucent teal surface of revolution from
    `gaussian_tube_section`), a patch S on its wall at t = 0 (orange) and the same particles' patch S′ after the tube has
    stretched (orange, darker), vortex lines drawn on the wall through both. *see:* "the patch is stretched along the tube
    and narrowed round it, but it stays on the wall." *read:* "Helmholtz 1: the patch's edge has zero circulation for
    ever, so it never tilts off the wall." *change:* "…viscosity: the patch's circulation would leak and the lines could
    reconnect — the patch would drift off the wall on the diffusion time scale."
13. `nb.md` — pointer to the explainer (already embedded in C03): "🎮 Open **Stretch and tangle a loop of dye — what stays
    the same?** above and switch on its *Helmholtz* mode: a small loop on a tube wall keeps zero flux while it is
    carried." (builder: `nb.pointer`, not a second `nb.explainer` — each explainer is embedded once).
14. `nb.md` — **What would change if…** "…the vortex were a straight line vortex in a flow that carries it along? Then
    Helmholtz 1 says the vortex goes wherever the local flow (made by everything else) takes it — the rule of §5.7."

### A.4 §5.4 Vorticity Equation in a Nonrotating Frame — R10 R11, C06 (+N18 N19 N20 N22 N23 · D09 ★★★; pointer N21 → C10)
1. `nb.section("5.4", "Vorticity Equation in a Nonrotating Frame", intro="**What is this section about?** Taking the curl
   of Navier–Stokes removes pressure and gravity and leaves an equation for the vorticity alone: vorticity is carried
   by the flow, stretched and tilted by velocity gradients along it, and diffused by viscosity. This section derives it
   for a fluid of constant density in an inertial frame; §5.6 adds rotation and density variations.")`
2. `nb.recap("R10", "Vorticity has no divergence", "$\nabla\cdot\boldsymbol\omega=\nabla\cdot(\nabla\times\mathbf u)=0$ for every smooth flow —
   the divergence of a curl vanishes (Ch. 2). §5.6 proves it again in index form as (5.18).", where="Ch. 2 §2.13")`
3. `nb.recap("R11", "The Lamb identity", "$(\mathbf u\cdot\nabla)\mathbf u=\nabla(\tfrac12\mathbf u\cdot\mathbf u)+\boldsymbol\omega\times\mathbf u$ (4.68): the advective
   acceleration is the gradient of the kinetic energy plus the Lamb vector ω × u. ⚠️ The book's §5.4 line writes
   ∇(u·u) without the ½ — harmless there, because the next move takes a curl, which kills any gradient.",
   where="Ch. 4 §4.9")` + `nb.code`: `print(ns.lamb_identity_terms(vortices.vortex_velocity_field("gaussian", Gamma=1.0,
   sigma=0.1), np.array([0.07, 0.02]))["residual"])` → *expect:* ≤ 1e-8.

**C06 — The vorticity equation (5.13)**
4. `nb.core("C06", "The vorticity equation: $\\frac{D\\boldsymbol\\omega}{Dt}=(\\boldsymbol\\omega\\cdot\\nabla)\\mathbf u+\\nu\\nabla^2\\boldsymbol\\omega$ (5.13)",
   question="Stir a cup of tea and the swirl slowly spreads and fades; a tornado tightens as its column is stretched.
   What can change a fluid particle's vorticity — and why can pressure never do it?")`
5. `nb.md` — **The problem in plain words:** "Navier–Stokes tells how a particle's velocity changes. For swirl, we want the
   rule for its vorticity. Pressure is usually the hardest unknown in a flow; an equation without it is a gift — the
   basis of vorticity-based computer models (Ch. 10), of vortex decay (Ch. 8), of the turbulence cascade (Ch. 12) and of
   the quasi-geostrophic equations of Ch. 13."
6. `nb.md` — **The idea:**
   ```
   momentum:   Du/Dt = −∇p/ρ + g + ν∇²u         take the curl  ─────►   vorticity:  Dω/Dt = (ω·∇)u + ν∇²ω
                         ▲     ▲                                              carried   stretched/tilted   diffused
                  gradients: curl = 0  (pressure and gravity push through the centre of mass: no twist)
   ```
7. `nb.note` — **N18 [C]** "§5.4's hypotheses: constant density (so the fluid is barotropic), constant viscosity,
   conservative body force, inertial frame. §5.6 (C09) relaxes the first and the last."
8. `nb.note` — **N19 [B]** "The first move of the derivation is (5.12), the curl of incompressible Navier–Stokes (4.39b):"
   equation `\nabla\times\Big\{\frac{D\mathbf u}{Dt}=-\frac1\rho\nabla p+\mathbf g+\nu\nabla^2\mathbf u\Big\}`, ref "5.12". "The curls of −∇p/ρ and of
   g = −∇Φ vanish because both are gradients."
9. `nb.note` — **N20 [B]** "**A vector identity we need** (the book's (B.3.10)), in its general form:" equation
   `\nabla\times(\boldsymbol\omega\times\mathbf u)=(\mathbf u\cdot\nabla)\boldsymbol\omega-(\boldsymbol\omega\cdot\nabla)\mathbf u+\boldsymbol\omega(\nabla\cdot\mathbf u)-\mathbf u(\nabla\cdot\boldsymbol\omega)`, no ref. "It
   follows from the ε–δ identity (2.19) and the product rule (Ch. 2 D09's moves); the sympy cell of D09 checks it for
   arbitrary fields. For incompressible flow the last two terms drop."
10. `nb.derivation("D09", …)` — Part F D09 (12 steps + sympy `check_src`), ref "5.13".
11. `nb.worked_example("the three terms at one point of Burgers' vortex", "Burgers' vortex (derived in C10):
    u_R = −αR/2, u_z = αz, ω_z = (αΓ/4πν)e^{−αR²/4ν} with α = 1 s⁻¹, Γ = 10⁻³ m²/s, ν = 10⁻⁶ m²/s. At R = 1 mm:
    1. ω_z = 79.58 e^{−0.25} = 61.97 s⁻¹. 2. Its radial slope: ω_z′ = −(αR/2ν)ω_z = −500 × 61.97 = −3.10×10⁴ s⁻¹/m.
    3. Advection: u_R ω_z′ = (−0.5×10⁻³)(−3.10×10⁴) = +15.49 s⁻². 4. Stretching: ω_z ∂u_z/∂z = 61.97 × 1 = +61.97 s⁻².
    5. Steady, so (5.13) says diffusion = advection − stretching = 15.49 − 61.97 = −46.48 s⁻²: viscosity carries
    vorticity outward, stretching makes more, the inflow brings it back.")`
12. `nb.code` — *code:* `for name, R in (("burgers", 1e-3), ("lamb_oseen", 5e-3), ("taylor_green", 0.3)):` `print(name,
    ch05.vorticity_budget_preset(name, x=R, y=0.0, z=0.0, component=2))` · `u_b = vortices.burgers_vortex_field(1e-3, 1.0,
    1e-6)`; `print(ch05.vorticity_terms(u_b, np.array([1e-3, 0.0, 0.0]), nu=1e-6, h=1e-5))`. *expect:* Burgers: local 0,
    advective +15.49, stretching_tilting +61.97, diffusion −46.48, residual ≤ 1e-4 (planetary, baroclinic 0) ·
    Lamb–Oseen: stretching 0, local = diffusion, advective 0 · Taylor–Green (2-D): stretching 0 · the `VorticityTerms`
    vectors with the same z-numbers. *explain:* 1. the parity-friendly preset (one component at one point) for three
    classic flows; 2. the general call: vorticity by stencils of u, then each term of (5.13) by stencils of ω (h = 10⁻⁵ m
    here, ≈ R/100); 3. the residual is the check that (5.13) holds.
13. `nb.check_agree` — **from scratch (curation §7):** central differences by hand at x₀ = (1 mm, 0, 0): `h = 1e-5`;
    `w = lambda x: kinematics.vorticity(u_b, x, h=h)` (ω by stencils); `G = kinematics.velocity_gradient_at(u_b, x0,
    h=h)`; `stretch = G @ w(x0)` ((ω·∇)u = Gω); `lap = sum((w(x0 + h*e) - 2*w(x0) + w(x0 - h*e))/h**2 for e in
    np.eye(3))`; `assert np.allclose(stretch[2], 61.97, rtol=1e-3)`; `assert np.allclose(1e-6*lap[2], -46.48, rtol=1e-3)`;
    `assert np.allclose([stretch[2], 1e-6*lap[2]], [T.stretching_tilting[2], T.diffusion[2]], rtol=1e-6)`. Markdown: "(ω·∇)u
    is just G times ω; ∇²ω is three second differences."
14. `nb.figure` — term bars (local blue dashed, advective grey, stretching purple, diffusion rose, residual black) at
    probe points of three flows side by side: Lamb–Oseen at R = 5 mm (t₀ = 10 s), Taylor–Green, Burgers at R = 1 mm; y
    axis symmetric-log. Title "Three flows, three balances". *see:* "Lamb–Oseen: only local and diffusion bars (equal);
    Taylor–Green: no stretching bar; Burgers: advective + stretching balance diffusion, no local bar." *read:* "in 2-D
    flows vorticity is only carried and diffused; stretching needs a velocity gradient along ω (C10)." *change:*
    "…ν halved in Burgers: the core narrows by √2 and every bar at a fixed R changes, but the three still add to zero."
15. `nb.note` — **N22 [B]** "**A sheet of vorticity diffusing.** Start with all the vorticity on a plane, ω = γδ(y)e_z (a
    velocity jump of γ, C14). (5.13) reduces to 1-D diffusion, ∂ω_z/∂t = ν∂²ω_z/∂y², whose solution (the book leaves it to
    Exercise 5.6) is" equation `\omega_z(y,t)=\frac{\gamma}{2\sqrt{\pi\nu t}}\exp\Big\{-\frac{y^2}{4\nu t}\Big\}`, no ref, "with velocity
    $u=-\frac\gamma2\,\mathrm{erf}\big(\frac{y}{2\sqrt{\nu t}}\big)$ (gloss: erf η = 1 − erfc η, Ch. 4 P123). The total ∫ω dy = γ never changes;
    the width 2√(νt) grows. Number: γ = 1 m/s in water, after 1 s the layer is 2 mm thick and the peak vorticity 282 s⁻¹.
    → Stokes' first problem (Ch. 8) is half of this; shear layers (Ch. 11)."
16. `nb.code` — *code:* `y = np.linspace(-0.01, 0.01, 2001)` · `for t in (0.1, 1.0, 10.0): u, w = ch05.diffusing_vortex_sheet(y,
    t, 1.0, 1e-6); print(t, np.trapezoid(w, y), w.max(), u[0], u[-1])` · FTCS cross-check: `from fluidpy.core import
    diffusion; w0 = ch05.diffusing_vortex_sheet(y, 0.1, 1.0, 1e-6)[1]; w_num = diffusion.ftcs_diffusion_1d(w0, 1e-6, y[1]-y[0],
    dt, nsteps, …)` to t = 1 s → compare with the exact profile. *expect:* ∫ω dy = 1.0000 (to 1e-6) at each t; peak 892.1,
    282.1, 89.2 s⁻¹; u(∓1 cm) → ±0.5 m/s (0.5000 at t = 0.1 and 1 s, 0.4873 at 10 s: the layer, 2√(νt) = 6.3 mm thick, reaches the
    window's edge); FTCS vs exact ≤ 1e-3 relative. *explain:* 1. the closed form; 2. the conserved total;
    3. the far velocities ±γ/2; 4. an independent finite-difference solution of the same diffusion equation agrees.
17. `nb.plotly` — `slider_figure` over t = 0.05…10 s (20 steps, FAST 10, log-spaced): ω_z(y, t) (teal) and u(y, t)
    (orange, right axis) for γ = 1 m/s, ν = 10⁻⁶; y ∈ [−2, 2] cm; the value of ∫ω dy printed in the title ("∫ω dy = 1.000
    m/s at every t"). *see:* "a sharp spike that widens and sinks; a velocity step that softens." *read:* "vorticity
    spreads like heat; its total — the velocity jump — is conserved." *change:* "…ν ten times larger (air): the same
    shapes at t/10."
18. `nb.note` — **N23 [C]** "**Hill's spherical vortex** is a second exact solution, steady and inviscid (the book leaves
    it to Exercise 5.11): inside a sphere of radius a, $\psi=\frac{Aa^4}{10}\frac{R^2}{a^2}\Big(1-\frac{R^2}{a^2}-\frac{z^2}{a^2}\Big)$ and
    ω = AR e_φ — ring-shaped vortex lines stretched just enough to keep ω/R constant on each streamline. The flow outside
    the sphere is Ch. 6's." + `nb.code`: `print(ch05.vorticity_budget_preset("hill", x=0.3, z=0.2, component=1))` →
    *expect:* residual ≤ 1e-6, stretching = advective (steady, ν = 0).
19. `nb.pointer` — **N21 → C10**: "Burgers' vortex (an exercise attached to this section) is taught with stretching in
    C10, where its stretching–diffusion balance is the point (D18)."
20. `nb.md` — **Helmholtz's first theorem from the field equation (C05, second proof).** "With ν = 0, (5.13) reads
    Dω/Dt = (ω·∇)u = Gω — the same linear equation as a material line element, D(δx)/Dt = Gδx (Ch. 3 P99). Start δx
    along ω and they stay parallel, with \|ω\|/\|δx\| constant (Cauchy): vortex lines are material." + `nb.code`: `u_abc =
    ch05.abc_flow()` (a steady Euler flow with ω = u) · `r = ch05.frozen_in_check(u_abc, np.zeros(3), 1e-3*u_abc(np.zeros(3),
    0.0), (0, 5), np.linspace(0, 5, 11))` · `print(np.max(r["angle"]), np.ptp(r["ratio"]))` · viscous: `r2 =
    ch05.frozen_in_check(vortices.burgers_vortex_field(1e-3, 1.0, 1e-6), np.array([1e-3, 0, 1e-3]), np.array([0, 0,
    1e-6]), (0, 2), np.linspace(0, 2, 5), nu=1e-6)`; `print(np.ptp(r2["ratio"]))`. *expect:* angle ≤ 1e-8 rad, ratio spread ≤
    1e-8 (inviscid) · ratio drift ≥ 1e-2 (viscous: diffusion unhooks ω from the fluid). + `nb.figure`: angle(t) and ratio(t)
    for both runs (teal inviscid flat, rose viscous drifting). *see / read / change:* "the inviscid curves lie on 0 and 1;
    the viscous ratio drifts as diffusion moves vorticity relative to the fluid." *change:* "…ν → 10ν: the drift grows ten
    times faster."
21. `nb.md` — "> ⚠️ **Common confusion:** '(5.13) is linear in ω.' It only looks linear: u is itself determined by ω
    (Biot–Savart, C07), so the advective term (u·∇)ω and the stretching term (ω·∇)u are both quadratic in the vorticity."
22. `nb.md` — pointer to the explainer (embedded at the end of C10): "🎮 **How can a flow spin fluid faster without a
    torque?** (after C10) steps through D09 in its Derivation tab." · **What would change if…** "…the density varied or
    the frame rotated? The pressure curl would no longer vanish (a baroclinic source, C04) and the Coriolis term would add
    the planet's vorticity to the stretching term — (5.30), C09."

### A.5 §5.5 Velocity Induced by a Vortex Filament: Law of Biot and Savart — C07 (+N24–N28 N52 · D10 ★★★ D11 ★★★), C08 (+N29 · D12 D13 · E6)
1. `nb.section("5.5", "Velocity Induced by a Vortex Filament: Law of Biot and Savart", intro="**What is this section
   about?** If you know the vorticity everywhere, you know the velocity everywhere (up to a potential flow). This section
   derives the formula — Biot–Savart, the same law that gives a magnetic field from an electric current — and its
   working form for thin vortex filaments. The book's printed (5.14) carries a sign slip; we derive the correct sign.")`

**C07 — The Biot–Savart law (5.16)**
2. `nb.core("C07", "Biot–Savart: $\\mathbf u=\\frac1{4\\pi}\\int_{V'}\\frac{\\boldsymbol\\omega\\times(\\mathbf x-\\mathbf x')}{\\lvert\\mathbf x-\\mathbf x'\\rvert^3}d^3x'$ (5.16)",
   question="A vortex spins in one corner of a pond. Water on the far side moves too, although nothing touches it there.
   Given where the vorticity is, is the velocity everywhere decided?")`
3. `nb.md` — **The problem in plain words:** "Aerodynamicists replace a wing by the vortices it sheds and compute the
   downwash they induce; oceanographers invert a map of vorticity for the currents; vortex methods in CFD move blobs of
   vorticity with the velocity they induce on each other. All of them need one formula: velocity from vorticity."
4. `nb.md` — **The idea:**
   ```
   ω = ∇×u,  ∇·u = 0   ──curl again──►   ∇²u = −∇×ω          (a Poisson equation: like gravity or electrostatics)
   point-source solution of ∇²: −1/(4πr)  ──add up──►  u = (1/4π)∫ (∇'×ω)/r d³x'      (5.14), sign corrected
   move the curl off ω (product rule + Gauss)  ──►  u = (1/4π)∫ ω×(x−x')/r³ d³x'   (5.16)  ≡ current → magnetic field
   ```
5. `nb.primer("Poisson equation and Green's function", "A Poisson equation ∇²φ = q asks for the field produced by a
   source density q (gravity from mass, voltage from charge). Its building block is the response to one point source,
   the Green's function: in 3-D, G = −1/(4π\|x − x′\|) solves ∇²G = δ(x − x′). Because the equation is linear, the
   response to any q is the sum of point responses, φ(x) = ∫G(x, x′)q(x′)d³x′.", code="x, y, z = sp.symbols('x y z')
   # a point in space (source at the origin)\nr = sp.sqrt(x**2 + y**2 + z**2)                                # distance from the source\nG = -1/(4*sp.pi*r)                                              # the Green's function of ∇²\nprint(sp.simplify(sum(sp.diff(G, v, 2) for v in (x, y, z))))    # ∇²G = 0 away from the source")`
   — *expect:* `0`.
6. `nb.primer("gradient of 1/distance with respect to the source point", "With x fixed and the source point x′ moving,
   ∇′(1/\|x − x′\|) = +(x − x′)/\|x − x′\|³ — it points from the source toward the field point. Differentiating with respect
   to x instead flips the sign; the book's two sign slips in §5.5 are exactly this trap.", code="X = sp.symbols('x y z');
   Xp = sp.symbols('xp yp zp')                  # field point, source point\nd = [a - b for a, b in zip(X, Xp)]; r = sp.sqrt(sum(di**2 for di in d))\ngrad_p = [sp.diff(1/r, v) for v in Xp]                                 # ∇′(1/r)\nprint([sp.simplify(g - di/r**3) for g, di in zip(grad_p, d)])           # [0, 0, 0]")`
   — *expect:* `[0, 0, 0]`.
7. `nb.primer("curl of a product (product rule)", "For a scalar f and a vector A, ∇×(fA) = f∇×A + ∇f×A — the
   product rule, with the order of the cross product kept. D11 uses it to move the curl off ω.", code="x, y, z =
   sp.symbols('x y z'); f = sp.Function('f')(x, y, z)\nA = sp.Matrix([sp.Function(n)(x, y, z) for n in 'PQR'])       # a generic vector field\ncurl = lambda V: sp.Matrix([sp.diff(V[2], y) - sp.diff(V[1], z), sp.diff(V[0], z) - sp.diff(V[2], x), sp.diff(V[1], x) - sp.diff(V[0], y)])\ngrad = sp.Matrix([sp.diff(f, v) for v in (x, y, z)])\nprint(sp.simplify(curl(f*A) - f*curl(A) - grad.cross(A)))       # zero vector")`
   — *expect:* `Matrix([[0], [0], [0]])`. (Gloss in the text: a × b = −b × a.)
8. `nb.note` — **N24 [B]** "**Velocity obeys a Poisson equation.** For incompressible flow the curl of the vorticity is
    minus the Laplacian of the velocity:" equation `\nabla\times\boldsymbol\omega=\nabla\times(\nabla\times\mathbf u)=\nabla(\nabla\cdot\mathbf u)-\nabla^2\mathbf u=-\nabla^2\mathbf u`,
    no ref. "(the curl-of-curl identity, Ch. 4 P122). On a periodic box the Poisson equation is solved in one line with
    Fourier modes — the next primer and cell. → Vorticity–stream-function CFD (Ch. 10), PV inversion (Ch. 13)."
9. `nb.primer("Fourier modes and the FFT Poisson solver", "On a periodic box every smooth field is a sum of waves e^{i k·x};
   a derivative ∂/∂x multiplies a wave by i k_x and ∇² by −\|k\|². So ∇²ψ = −ω becomes ψ̂ = ω̂/\|k\|² for each wave
   (k = 0 is the mean, set to 0). `numpy.fft.fft2` finds the waves, `ifft2` adds them back.", code="n, L = 32, 2*np.pi;
   x = np.arange(n)*L/n; X, Y = np.meshgrid(x, x, indexing='ij')\nw = 2*np.sin(X)*np.sin(Y)                                  # vorticity of ψ = sin x sin y\nk = np.fft.fftfreq(n, L/n)*2*np.pi; KX, KY = np.meshgrid(k, k, indexing='ij'); K2 = KX**2 + KY**2; K2[0, 0] = 1\npsi = np.real(np.fft.ifft2(np.fft.fft2(w)/K2))                   # ψ̂ = ω̂/k²\nprint(np.max(abs(psi - np.sin(X)*np.sin(Y))))                    # ~1e-15")`
   — *expect:* ≤ 1e-14. + `nb.code`: `n = 32; …; U = ch05.velocity_from_vorticity_fft(w, 2*np.pi)`; compare with the
   cellular flow u = sin x cos y, v = −cos x sin y → *expect:* max error ≤ 1e-13.
10. `nb.derivation("D10", …)` — Part F D10 (12 steps + sympy `check_src`), ref "5.14".
11. `nb.note` — **N25 [B]** "**(5.14) with its sign corrected.** D10's result is" equation
    `\mathbf u(\mathbf x,t)=\frac1{4\pi}\int_{V'}\frac{\nabla'\times\boldsymbol\omega(\mathbf x',t)}{\lvert\mathbf x-\mathbf x'\rvert}d^3x'`, ref "5.14". "> ⚠️ **The book prints**
    $\mathbf u(\mathbf x,t)=-\frac{1}{4\pi}\int_{V'}\frac{1}{\lvert\mathbf x-\mathbf x'\rvert}(\nabla'\times\boldsymbol\omega(\mathbf x',t))d^3x'$ **(5.14) with −1/(4π); the
    correct factor is +1/(4π).** It also cites 'Exercise 5.8' for the proof; the Green's-function exercise is 5.9. The
    code below computes both signs on a smooth vortex tube — only one of them spins the right way."
12. `nb.md` — gloss: "**3-D quadrature.** We integrate over a box with a tensor product of Gauss–Legendre nodes (the 1-D
    rule of Ch. 2 used along each axis); points inside the vorticity use a slightly smoothed kernel \|r\|² → \|r\|² + ε²."
    (builder: this is the primer below.)
13. `nb.primer("Gauss–Legendre quadrature in 3-D and a smoothed kernel", "Gauss–Legendre picks n nodes and weights per
    axis so that polynomials up to degree 2n − 1 integrate exactly; a box uses all n³ combinations (weights multiply).
    Near a singular kernel like 1/r³ we replace r² by r² + ε² (a tiny smoothing length) so no node ever divides by
    zero.", code="from fluidpy.core.integral_theorems import gauss_legendre_nodes\nx, w = gauss_legendre_nodes(0.0, 1.0, 4)            # 4 nodes on [0, 1]\nX, Y, Z = np.meshgrid(x, x, x, indexing='ij'); W = w[:, None, None]*w[None, :, None]*w[None, None, :]\nprint(np.sum(W*X**3*Y**2*Z), 1/4*1/3*1/2)              # exact for a polynomial: 0.041667 twice")`
   — *expect:* `0.0416667 0.0416667`.
14. `nb.code` — **the sign test** — *code:* `F = ch05.gaussian_tube_fields(Gamma=1.0, sigma=0.1, L=4.0)` (a straight tube
    along z with smooth ends) · `nodes, weights = …` 24³ Gauss–Legendre nodes on `F["bounds"]` (FAST 16³) · `for s in (+1,
    -1): print(s, ch05.velocity_from_curl_omega(F["curl_omega"], np.array([0.5, 0.0, 0.0]), nodes, weights, sign=s)[1])` ·
    `print(1.0*(1 - np.exp(-25))/(2*np.pi*0.5))` · `print(ch05.biot_savart_volume(F["omega"], np.array([0.5, 0.0, 0.0]),
    nodes, weights)[1])`. *expect (orchestrator note: for the finite L = 4 tube the exact reference is 0.308838, `u_theta_reference`; 0.318310 is the infinite line):* +0.3183 (sign +1) and −0.3183 (printed sign), exact 0.318310, Biot–Savart 0.3183 (all
    to ≈ 1e-3 at 24³, the tube's finite length included). *explain:* 1. a smooth vortex tube, its vorticity and the curl
    of its vorticity as callables; 2. nodes and weights; 3. (5.14) with both signs: the book's −1/(4π) turns the swirl
    backwards; 4. the line-vortex value Γ/2πr; 5. (5.16) directly — same answer as the corrected (5.14).
15. `nb.derivation("D11", …)` — Part F D11 (12 steps + sympy `check_src`), ref "5.16".
16. `nb.note` — **N26 [B]** "**The rewrite of the integrand** (D11 steps 2–5), with the correct sign:" equation
    `\frac{\nabla'\times\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}=\nabla'\times\Big(\frac{\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}\Big)+\boldsymbol\omega\times\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}`, no ref. "⚠️ The book's second
    line prints $+\big(\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}\big)\times\boldsymbol\omega$; it should be minus that (= +ω × (x − x′)/\|x − x′\|³). This second slip
    cancels the first, so (5.16) is right."
17. `nb.note` — **N27 [B]** "**Gauss' theorem in curl form**, used in D11 steps 7–9:" equation
    `\int_{V'}\nabla'\times\Big(\frac{\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}\Big)d^3x'=\int_{A'}\frac{\mathbf n\times\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}d^2x'`, ref "5.15". "Any smooth
    field F: ∫_V∇×F dV = ∮_A n×F dA (apply (2.30) to each component ε_kij F_j)." + `nb.code`: `Fpoly = lambda X, Y, Z:
    np.array([X*Y**2, Y*Z + X**3, Z*X**2*Y])` (a polynomial field) · `print(ch05.curl_theorem_box(Fpoly, [(0, 1), (0, 2), (-1,
    1)], n=6))` → *expect:* volume = surface to 1e-12.
18. `nb.note` — **N28 [B]** "**Why the surface term vanishes** (Fig. 5.8): V′ is a piece of tube with flat ends
    perpendicular to ω (n × ω = 0 there) and a side wall outside the vortex (ω = 0 there)." + `nb.code`: the curl theorem on
    a box around a piece of the Gaussian tube, `ch05.curl_theorem_box(lambda X, Y, Z: F["omega"](np.array([X, Y, Z]), 0.0) /
    np.sqrt((X - 0.5)**2 + Y**2 + Z**2), [(-0.4, 0.4), (-0.4, 0.4), (-0.5, 0.5)], n=24)` → *expect:* surface term ≈ 0 up to
    the tube-end contributions (≤ 1e-3 of the volume term's scale). *explain:* "with a box whose sides are 4 core radii
    out and whose ends cut ω at right angles, the surface integral of (5.15) vanishes."
19. `nb.note` — **N52 [B] (Fig. 5.8 analogue)** + `nb.figure`: a curved filament (a quarter ring) in 3-D (matplotlib
    `projection='3d'`), a field point x, the vector x − x′ from one element, the element's contribution du (purple arrow,
    perpendicular to both e_ω and x − x′), and the dashed volume V′ around that element. *see:* "one element, its lever
    vector to the field point and its push." *read:* "each element pushes perpendicular to itself and to the line joining
    it to x, with strength ∝ 1/distance² (5.16)." *change:* "…the field point moved onto the filament's axis of symmetry:
    all the pushes line up along the axis."
20. `nb.worked_example("a vortex's velocity from its vorticity, by hand", "A Gaussian vortex with Γ = 1 m²/s and core
    σ_c = 0.1 m (ω_z = (Γ/πσ_c²)e^{−r²/σ_c²}, peak 31.8 s⁻¹). 1. Outside the core, r = 0.5 m, nearly all the vorticity is
    inside: the planar Biot–Savart sum must give Γ(1 − e^{−25})/(2π × 0.5) = 0.3183 m/s. 2. At the core edge r = 0.1 m:
    Γ(1 − e⁻¹)/(2π × 0.1) = 0.6321/0.6283 = 1.006 m/s. 3. Inside, r = 0.05 m: Γ(1 − e^{−0.25})/(2π × 0.05) = 0.2212/0.3142 =
    0.704 m/s. The vorticity *outside* a circle adds nothing there — the planar version of Newton's shell theorem.")`
21. `nb.code` — *code:* `x = np.linspace(-0.5, 0.5, 201); X, Y = np.meshgrid(x, x)` (grid, spacing 5 mm) · `wz =
    vortices.gaussian_vortex(np.hypot(X, Y), 1.0, 0.1)[1]` · `for r in (0.05, 0.1, 0.5): print(r, ch05.biot_savart_2d(wz, X,
    Y, np.array([r, 0.0]), eps=1e-4)[1], vortices.gaussian_vortex(r, 1.0, 0.1)[0])`. *expect:* grid sums 0.70360, 1.00558, 0.318310
    against the exact 0.704101, 1.006051, 0.318310 — within 7×10⁻⁴ inside the core (the kernel's 1/r singularity next to
    the grid node) and 4×10⁻⁸ outside, where ω ≈ 0. *explain:* 1. the
    vorticity sampled on a 201 × 201 grid; 2. the planar Biot–Savart sum (kernel e_z × r/2π\|r\|²) at three radii;
    3. agreement with the exact profile inside and outside the core.
22. `nb.check_agree` — **from scratch (curation §7):** a direct double loop over grid cells: `u = v = 0.0`; `for j, i in
    np.ndindex(X.shape): dx, dy = 0.1 - X[j, i], 0.0 - Y[j, i]; r2 = dx*dx + dy*dy + 1e-8; u += -wz[j, i]*dy/r2*dA/(2*np.pi);
    v += wz[j, i]*dx/r2*dA/(2*np.pi)` (FAST: every 2nd cell with 4 × dA) · `assert np.allclose(v, ch05.biot_savart_2d(wz, X, Y,
    np.array([0.1, 0.0]), eps=1e-4)[1], rtol=1e-6)`. Markdown: "Forty thousand little pushes, added by hand."
23. `nb.plotly` — `slider_figure` over the grid size n = 21…201 (10 steps): u_θ(r) from `biot_savart_2d` (teal dots) vs the
    exact Gaussian profile (black), and the printed-sign curve (rose, mirror image below zero), r ∈ [0, 0.5] m. Title
    "Biot–Savart rebuilds the vortex — with the right sign". *see:* "teal dots settle onto the black curve as the grid
    refines; the rose curve is its upside-down twin." *read:* "(5.16) recovers the swirl inside and outside the core;
    the book's −1/(4π) would spin every vortex backwards." *change:* "…a Rankine vortex instead: the dots converge more
    slowly (its vorticity jumps at the core edge), to the kink at r = σ_c."
24. `nb.md` — **What would change if…** "…the vortex were thin and we stood far from it? Then the details of the core stop
    mattering: only its strength Γ and its path count — the filament law (C08)."

**C08 — The filament law (5.17)**
25. `nb.core("C08", "The filament law: $d\\mathbf u=\\frac{\\Gamma\\,dl}{4\\pi}\\mathbf e_\\omega\\times\\frac{\\mathbf x-\\mathbf x'}{\\lvert\\mathbf x-\\mathbf x'\\rvert^3}$ (5.17)",
    question="A smoke ring pushes itself forward; a wing-tip vortex drags air down behind a plane. How much velocity does
    one short piece of a thin vortex make?")`
26. `nb.md` — **The problem in plain words:** "Real vortices are often thin compared with the distances that matter: the
    trailing vortices of an aircraft, a tornado seen from a kilometre away, a smoke ring. Then we can replace the tube by a
    line carrying a circulation Γ and add up its pieces — the lifting-line theory of wings (Ch. 14) is exactly this."
27. `nb.md` — **The idea:**
    ```
    thin tube of strength Γ, seen from far away            one piece of length dl pushes the fluid at x:
    ─────►══════════════════►  (all vorticity on a line)     du = (Γ dl/4π) e_ω × (x − x')/|x − x'|³
                                                            sum along a straight line → Γ/4πd (cos θa − cos θb)
                                                            infinite line → Γ/2πd  (the line vortex (5.2) again)
    ```
28. `nb.md` — glosses: "**Far-field ('frozen kernel').** If a function barely changes across a small region, pull it
    out of the integral: moving across a core of radius a changes 1/\|x − x′\|³ by about 3a/\|x − x′\| — 3 % at ten core
    radii. **1 + cot²θ = csc²θ** (divide sin² + cos² = 1 by sin²)."
29. `nb.derivation("D12", …)` — Part F D12 (6 steps), ref "5.17".
30. `nb.primer("improper integral as a limit", "An integral to infinity means: integrate to a finite end, then let the end
    run away. ∫₀^∞ e^{−x} dx is the limit of 1 − e^{−b} as b → ∞, i.e. 1. D13 lets a straight segment grow into an
    infinite line this way.", code="from scipy.integrate import quad\nf = lambda l: 1.0/(1.0 + l**2)**1.5                     # d = 1: the segment integrand\nfor b in (1, 10, 100):\n    print(b, quad(f, -b, b)[0])                              # → 2 as b → ∞ (so Γ·2/(4π) = Γ/2π)\nprint(quad(f, -np.inf, np.inf)[0])                           # scipy does the limit for you: 2.0")`
    — *expect:* `1 1.41421`, `10 1.99007`, `100 1.99990`, `2.0`.
31. `nb.derivation("D13", …)` — Part F D13 (8 steps), no ref (builder passes ref "5.2" — the result it recovers).
32. `nb.note` — **N29 [B]** "**The loop closes.** An infinite straight filament from (5.17) gives back $u_\theta=\Gamma/2\pi r$
    (5.2); a finite segment gives $\frac{\Gamma}{4\pi d}(\cos\theta_a-\cos\theta_b)$ — the book never writes this. Numbers (Γ = 1 m²/s,
    d = 1 m): infinite line 0.159 m/s, semi-infinite 0.080 m/s, the segment from −1 to +1 m 0.113 m/s; a square loop of
    side 1 m at its centre 2√2Γ/πL = 0.900 m/s; a ring of radius 0.5 m at its centre Γ/2R = 1.0 m/s."
33. `nb.worked_example("the push of a 2 m segment on a point 1 m away", "Γ = 1 m²/s, segment from z = −1 m to z = +1 m, field
    point at x = 1 m, z = 0 (d = 1 m). 1. The ends are seen at θ_a = 45° and θ_b = 135° (cos = ±0.707). 2. u = Γ/(4πd) ×
    (0.707 + 0.707) = 0.0796 × 1.414 = 0.1125 m/s. 3. Direction e_ω × e_d = e_z × e_x = e_y. 4. Stretch the segment to
    ±10 m: cosines ±0.995, u = 0.0796 × 1.990 = 0.1584 m/s — within 0.5 % of the infinite line's 0.1592 m/s.")`
34. `nb.code` — *code:* `print(ch05.segment_induced_velocity(np.array([1.0, 0, 0]), np.array([0, 0, -1.0]), np.array([0,
    0, 1.0]), 1.0))` · `for L in (1, 10, 100, 1000): print(L, ch05.segment_induced_velocity(np.array([1.0, 0, 0]), np.array([0,
    0, -L]), np.array([0, 0, L]), 1.0)[1])` · `print(ch05.filament_velocity_preset("square", 0.0, 0.0, 0.0, side=1.0,
    component=2))` · `for M in (8, 16, 64): print(M, ch05.filament_velocity(np.zeros(3), ch05.filament_preset("ring", M=M,
    R=0.5), 1.0)[2])` · `print(ch05.ring_axis_velocity(0.0, 0.5, 1.0), ch05.ring_axis_velocity(0.5, 0.5, 1.0))`. *expect:*
    `[0, 0.112540, 0]` · 0.112540, 0.158367, 0.159147, 0.159155 · 0.900316 · 1.054786, 1.013052, 1.000804 · 1.0 0.353553.
    *explain:* 1. the closed form of D13 for the worked example; 2. longer segments approach Γ/2πd; 3. a square loop
    (four segments); 4. a ring as an M-gon: the error falls as 1/M² (8 → 16 divides it by ≈ 4); 5. the exact ring-axis
    formula ΓR²/2(R² + z²)^{3/2}.
35. `nb.check_agree` — **from scratch (curation §7):** a midpoint sum of (5.17) along the straight segment: `l =
    (np.arange(20000) + 0.5)/20000*2 - 1` (midpoints on [−1, 1] m); `r = np.stack([1.0 - 0*l, 0*l, -l])` (x − x′);
    `du = np.cross(np.array([0, 0, 1.0]), r.T)/np.linalg.norm(r, axis=0)[:, None]**3 * (2/20000)/(4*np.pi)`;
    `assert np.allclose(du.sum(axis=0), ch05.segment_induced_velocity(np.array([1.0, 0, 0]), np.array([0, 0, -1.0]),
    np.array([0, 0, 1.0]), 1.0), atol=1e-10)`. Markdown: "Twenty thousand pieces of (5.17), summed, give the closed form."
36. `nb.plotly` — `slider_figure` over the number of segments M = 4…128 (10 steps, log-spaced): left trace, the M-gon ring
    (R = 0.5 m) seen from above; right trace (second subplot or offset axis), the on-axis velocity u_z(z) of the M-gon
    (teal dots) vs ΓR²/2(R² + z²)^{3/2} (black), z ∈ [−1.5, 1.5] m. Title "(5.17) summed round a ring". *see:* "a square,
    then a hexagon, then a circle; the teal dots close onto the black curve." *read:* "a ring of vorticity pushes fluid
    through its middle, strongest at the centre, Γ/2R; the polygon error falls like 1/M²." *change:* "…the ring radius
    doubled: the centre speed halves (Γ/2R), the curve widens."
37. `nb.figure` — speed at distance d = 1 m from a straight segment of half-length ℓ against ℓ (log axis, 0.1…100 m),
    with the infinite-line ghost Γ/2πd (grey dashed) and the semi-infinite half (dotted). *see / read / change:* "the
    curve climbs to the ghost once ℓ ≫ d; a semi-infinite line gives exactly half — the reason a wing's trailing vortex
    makes half the downwash at the wing than far behind it (Ch. 14)." *change:* "…d halved: the whole curve doubles and
    reaches the ghost at half the length."
38. `nb.explainer("biot_savart_filament", heading="How does a vortex push water far away?", why="Dragging a field point
    around straight, bent, square, ring and helical filaments while each segment's contribution and their vector sum are
    drawn makes the 1/distance² law and the perpendicular direction tangible; the sign toggle shows what the printed
    (5.14) would do.", tries=["Straight segment: lengthen it and watch the speed approach Γ/2πd.", "Ring: put the point on
    the axis and read Γ/2R at the centre.", "Click a segment: the inspector shows its dl, r and du.", "Flip the sign to the
    printed −1/(4π): every arrow reverses — the badge says why that is wrong."])`
39. `nb.md` — **What would change if…** "…the filament were curved and we asked for the velocity *on* it? The sum
    diverges logarithmically near the point itself: a curved vortex moves itself only because its core has a finite size
    — the smoke ring's self-speed Γ/4πR[ln(8R/a) − ¼] (N42, C13)."

### A.6 §5.6 Vorticity Equation in a Rotating Frame — R12–R17, C09 (+N30–N37 · D14 D15 ★★★), C10 (+N38 N21 N53 · D16 D17 D18 · E5), R18, C11 (+N39 N54 · D19 D20 · E7)
1. `nb.section("5.6", "Vorticity Equation in a Rotating Frame", intro="**What is this section about?** The ocean and the
   atmosphere are watched from a rotating Earth and have density contrasts. Redoing the vorticity equation in a frame
   turning at Ω, for a nearly incompressible (Boussinesq) fluid of variable density, adds two things: the planet's own
   vorticity 2Ω is stretched and tilted like any other, and density gradients across pressure gradients create vorticity.
   Two consequences follow: stretching and tilting of vortex lines (the engine of 3-D turbulence), and the conservation
   of absolute circulation — the seed of potential vorticity.")`
2. `nb.recap("R12", "Hypotheses and comma notation of §5.6", "A frame rotating at constant Ω; a non-barotropic but nearly
   incompressible (Boussinesq) fluid, so ∇·u ≈ 0; and the comma notation of Ch. 2 §2.14: $u_{i,j}\equiv\partial u_i/\partial x_j$,
   $u_{i,jk}\equiv\partial^2u_i/\partial x_j\partial x_k$. Gloss: at the end of a derivation we may rename the free index (n → i).",
   where="Ch. 4 §4.9 (Boussinesq), Ch. 2 §2.14")`
3. `nb.recap("R13", "Continuity", "Nearly incompressible flow: $u_{i,i}=0$ *(Eq. 5.19)* — Ch. 4's (4.10).",
   where="Ch. 4 §4.2")`
4. `nb.recap("R14", "Momentum in a steadily rotating frame", "$\frac{\partial u_i}{\partial t}+u_ju_{i,j}+2\varepsilon_{ijk}\Omega_ju_k=-\frac1\rho p_{,i}+g_i+\nu u_{i,jj}$
   *(Eq. 5.20)* — Ch. 4's (4.45) with Ω constant; g is the effective gravity (true gravity plus the centrifugal term).
   Ch. 4's explainer `rotating_frame_coriolis` shows the frame.", where="Ch. 4 §4.7")` + `nb.code`: `u_rest = lambda x,
   t: np.zeros(3)` (fluid at rest in a frame turning at 0.5 rad/s) · `p_h = lambda x, t: 1e5 - 1000*9.81*x[2]` ·
   `print(ch05.rotating_ns_residual(u_rest, p_h, 1000.0, np.array([0.3, 0.1, -0.2]), Omega=(0, 0, 0.5)))` → *expect:*
   `[0, 0, 0]` (1e-8): rest in the rotating frame with effective gravity is a solution.
5. `nb.recap("R15", "The Lamb identity in index form", "$u_ju_{i,j}=-(\mathbf u\times\boldsymbol\omega)_i+\tfrac12(u_j^2)_{,i}$ *(Eq. 5.21)* — the index
   form of (4.68) (R11).", where="Ch. 4 §4.9")`
6. `nb.recap("R16", "The antisymmetric gradient is the vorticity", "$\varepsilon_{ijk}\omega_k=u_{j,i}-u_{i,j}$ *(Eq. 5.22)* — Ch. 3's
   R_ij = −ε_ijkω_k (3.15)–(3.16), from ε–δ (2.19).", where="Ch. 3 §3.4")` + `nb.code`: `G = np.random.default_rng(5).
   normal(size=(3, 3))` (a random velocity gradient, G[i, j] = u_{i,j}) · `w = kinematics.vorticity_from_gradient(G)` ·
   `from fluidpy.core import tensors; E = tensors.levi_civita()` · `print(np.max(abs(np.einsum('ijk,k->ij', E, w) - (G.T -
   G))))` → *expect:* ≤ 1e-15 (the (i, j) entry of G is u_{i,j}, so u_{j,i} − u_{i,j} = (Gᵀ − G)_{ij}).
7. `nb.recap("R17", "The viscous force as a curl of the vorticity", "$\nu u_{i,jj}=-\nu\varepsilon_{ijk}\omega_{k,j}$ *(Eq. 5.23)*:
   ν∇²u = −ν∇×ω for ∇·u = 0 — Ch. 4's (4.40).", where="Ch. 4 §4.6")`

**C09 — The vorticity equation in a rotating frame with baroclinic generation (5.30)**
8. `nb.core("C09", "The full vorticity equation: $\\frac{D\\boldsymbol\\omega}{Dt}=(\\boldsymbol\\omega+2\\boldsymbol\\Omega)\\cdot\\nabla\\mathbf u+\\frac1{\\rho^2}\\nabla\\rho\\times\\nabla p+\\nu\\nabla^2\\boldsymbol\\omega$ (5.30)",
   question="The air over a warm coast on a rotating planet: where does its spin come from — the winds, the planet, or the
   warm land?")`
9. `nb.md` — **The problem in plain words:** "Weather maps are maps of vorticity: cyclones, troughs, fronts. The equation
   that governs it for the atmosphere and ocean must include the planet's rotation and the fact that warm and cold, fresh
   and salty fluid lie side by side. It is the same idea as (5.13) with two new sources — and every Ch. 13 theory
   (quasi-geostrophy, potential vorticity, baroclinic instability) starts from it."
10. `nb.md` — **The idea:**
    ```
    (5.13)  Dω/Dt = (ω·∇)u                       + ν∇²ω         constant ρ, inertial frame
    (5.30)  Dω/Dt = (ω + 2Ω)·∇u   + ∇ρ×∇p/ρ²    + ν∇²ω         Boussinesq, frame turning at Ω
                     ▲ the planet's vorticity   ▲ baroclinic torque (C04)
                       is stretched and tilted
                       like the fluid's own
    ```
11. `nb.note` — **N30 [B]** "**∇·ω = 0 in index form.** ω_i = ε_inq u_{q,n}, so" equation
    `\omega_{i,i}=(\varepsilon_{inq}u_{q,n})_{,i}=\varepsilon_{inq}u_{q,ni}=0\quad\text{or}\quad\nabla\cdot\boldsymbol\omega=0`, ref "5.18". "ε_inq is antisymmetric in i, n
    and u_{q,ni} symmetric (Schwarz): the contraction vanishes (Ch. 2). True for compressible and unsteady flow too." +
    `nb.code`: `print(ch05.vorticity_divergence(ch05.abc_flow(), np.array([0.3, 1.1, -0.4])))` → *expect:* ≤ 1e-8
    (round-off of nested stencils).
12. `nb.md` — gloss: "**Quotient rule for a gradient**: ∇(1/ρ) = −∇ρ/ρ², so ∇×(∇p/ρ) = ∇(1/ρ) × ∇p = −∇ρ×∇p/ρ²."
13. `nb.derivation("D14", …)` — Part F D14 (8 steps), ref "5.25".
14. `nb.note` — **N31 [B]** "D14's step 6 is (5.24), the Coriolis term with its dummy indices swapped:" equation
    `2\varepsilon_{ijk}\Omega_ju_k=-2\varepsilon_{ijk}\Omega_ku_j`, ref "5.24". "(Ω × u = −u × Ω.) Check with numbers: Ω = (0, 0, 1),
    u = (1, 0, 0): both sides give (0, 2, 0)." + one-line `nb.code`: `print(2*np.cross([0, 0, 1.0], [1.0, 0, 0]),
    -2*np.cross([1.0, 0, 0], [0, 0, 1.0]))` → `[0 2 0] [0 2 0]`.
15. `nb.note` — **N32 [B]** "D14's result is (5.25), rotating Navier–Stokes in Lamb form:" equation
    `\partial u_i/\partial t+\Big(\tfrac12u_j^2+\Phi\Big)_{,i}-\varepsilon_{ijk}u_j(\omega_k+2\Omega_k)=-(1/\rho)p_{,i}-\nu\varepsilon_{ijk}\omega_{k,j}`, ref "5.25".
    + `nb.code`: on the inertial-oscillation field u = U(cos ft, −sin ft, 0) (a solution of (5.20) with f = 2Ω, p
    hydrostatic): `print(ch05.rotating_lamb_form_terms(u_io, p_h, 1000.0, Phi, x0, t=1.0, Omega=(0, 0, 0.5))["residual"])`
    → *expect:* ≤ 1e-8; equal to `rotating_ns_residual` at the same point.
16. `nb.derivation("D15", …)` — Part F D15 (15 steps + sympy `check_src`), ref "5.30".
17. `nb.note` — **N33 [B]** "D15's step 1 is (5.26):" equation
    `\frac{\partial}{\partial t}\big(\varepsilon_{nqi}u_{i,q}\big)+\varepsilon_{nqi}\Big(\tfrac12u_j^2+\Phi\Big)_{,iq}-\varepsilon_{nqi}\varepsilon_{ijk}\big[u_j(\omega_k+2\Omega_k)\big]_{,q}=-\varepsilon_{nqi}\Big(\frac1\rho p_{,i}\Big)_{,q}-\nu\varepsilon_{nqi}\varepsilon_{ijk}\omega_{k,jq}`,
    ref "5.26". "⚠️ The book's sentence after it calls the scalar Π; it is Φ."
18. `nb.note` — **N34 [B]** "D15's steps 5–9 are (5.27):" equation
    `-\varepsilon_{nqi}\varepsilon_{ijk}\big[u_j(\omega_k+2\Omega_k)\big]_{,q}=-u_{n,j}(\omega_j+2\Omega_j)+u_j\omega_{n,j}`, ref "5.27". "⚠️ Between its
    second and third lines the book drops $u_{j,j}(\omega_n+2\Omega_n)$ without comment; it is zero because $u_{j,j}=0$ (5.19)."
19. `nb.note` — **N35 [B]** "D15's step 11 is (5.28), the baroclinic vector:" equation
    `-\varepsilon_{nqi}\Big(\frac1\rho p_{,i}\Big)_{,q}=-\frac1\rho\varepsilon_{nqi}p_{,iq}+\frac1{\rho^2}\varepsilon_{nqi}\rho_{,q}p_{,i}=0+\frac1{\rho^2}[\nabla\rho\times\nabla p]_n`,
    ref "5.28". "It is exactly the torque of C04 (D06). Zero when ρ = ρ(p)." + `nb.code`: `xs, ys, zs = sp.symbols('x y
    z')`; `print(ch05.baroclinic_term_sym(1000 - 2*xs, 1e5 - 9810*zs, (xs, ys, zs)).T)` (density falling with x, hydrostatic
    p) and `print(ch05.baroclinic_term_sym(sp.Function('R')(1e5 - 9810*zs), 1e5 - 9810*zs, (xs, ys, zs)).T)` (ρ = R(p)).
    *expect:* `[0, -19620/(1000 - 2x)**2, 0]` · `[0, 0, 0]`. *explain:* "a sideways density gradient under vertical
    gravity gives a horizontal vorticity source (the sea-breeze roll axis); a barotropic field gives none."
20. `nb.note` — **N36 [B]** "D15's steps 12–13 are (5.29):" equation
    `-\nu\varepsilon_{nqi}\varepsilon_{ijk}\omega_{k,jq}=-\nu(\delta_{nj}\delta_{qk}-\delta_{nk}\delta_{qj})\omega_{k,jq}=-\nu\omega_{k,nk}+\nu\omega_{n,jj}=\nu\omega_{n,jj}`, ref "5.29".
21. `nb.worked_example("a stretched column on the Earth and a lock exchange", "1. **Planetary stretching.** At the pole,
    Ω = 7.292×10⁻⁵ rad/s; a column of air at rest stretched at ∂w/∂z = 10⁻⁵ s⁻¹ (a weak convergence: 1 cm/s of updraft gained per
    kilometre of height): the planetary term 2Ω ∂w/∂z = 2 × 7.292×10⁻⁵ × 10⁻⁵ = 1.46×10⁻⁹ s⁻² — kept up for a day (86 400 s)
    it gives ζ ≈ 1.3×10⁻⁴ s⁻¹, comparable with f itself (1.46×10⁻⁴ s⁻¹ at the pole): cyclones spin up this way. 2. **Baroclinic.** The lock exchange of C04:
    2.42 s⁻² at the interface. 3. **Relative stretching** in Burgers' vortex (C06): 62 s⁻². The same equation spans
    thirteen orders of magnitude.")`
22. `nb.code` — *code:* `for name in ("burgers", "lock_exchange", "rotating_column"): print(name,
    ch05.vorticity_budget_preset(name))` (each at its scene's default probe) · `print(2*rotating.OMEGA_EARTH*1e-5)` ·
    `print(ch05.vorticity_budget_sym(u_poly, (xs, ys, zs), ts, nu, (0, 0, Om), rho_expr, p_expr)["reduces_to_513"])`.
    *expect:* Burgers: stretching 61.97, advective 15.49, diffusion −46.48, planetary 0, baroclinic 0, residual ≈ 0 ·
    lock exchange: baroclinic 2.4222, all else 0 · rotating column: planetary 0.2 (= 2Ωα with Ω = 0.5, α = 0.2), local 0.2,
    residual ≈ 0 · `1.458423e-09` · `True`. *explain:* 1. the budget of (5.30) at a probe of three scenes — one term
    dominates in each; 2. the Earth number of the worked example; 3. sympy: (5.30) with Ω = 0 and ρ = const is (5.13).
23. `nb.check_agree` — **from scratch (curation §7):** the baroclinic vector by central differences at the lock-exchange
    interface: `F = ch05.lock_exchange_fields(1000.0, 1025.0, 0.1)`; `h = 1e-4`; `x0 = np.array([0.0, 0.5])`; `drx = (F["rho"](x0 +
    [h, 0]) - F["rho"](x0 - [h, 0]))/(2*h)`, similarly `dry`, `dpx`, `dpy`; `b = (drx*dpy - dry*dpx)/F["rho"](x0)**2`;
    `assert np.allclose(b, ch05.baroclinic_term(F["rho"], F["p"], x0)[2], rtol=1e-8)`.
24. `nb.figure` — **N37 [B]** term bars of (5.30) (local blue dashed, advective grey, relative stretching/tilting purple,
    planetary amber, baroclinic orange, diffusion rose, residual black) for the three scenes, each normalised by its
    largest term (a symmetric-log inset gives the absolute sizes). Title "Same equation, three regimes". *see:* "Burgers
    — purple, grey and rose; lock exchange — one orange bar; rotating column — amber equals local." *read:* "the left side
    of (5.30) is the rate following a particle; ν∇²ω is molecular diffusion (as ν∇²u diffuses velocity); the baroclinic
    bar exists only where ∇ρ crosses ∇p; stretching matters even with Ω = 0." *change:* "…Ω switched off in the rotating
    column: the amber bar vanishes and nothing spins it up."
25. `nb.md` — "> ⚠️ **Common confusion:** 'Boussinesq means the density can be treated as constant.' In (5.30) ρ is kept
    inside the pressure term — that is where the baroclinic source lives. To leading order $\frac{1}{\rho^2}\nabla\rho\times\nabla p\approx\frac{1}{\rho_0}\nabla\rho'\times\mathbf g$
    with g = −g e_z (cf. (4.86)): only the density *perturbation* crossed with gravity matters."
26. `nb.md` — pointer: "🎮 **What does a spinning planet add to the vorticity equation?** (after C11) walks through D14 and
    D15 with term bars in its Derivation tab." · **What would change if…** "…we looked at one vortex line in an inviscid,
    barotropic, inertial flow? Only (ω·∇)u would remain — and it splits cleanly into stretching and tilting (C10)."

**C10 — Stretching and tilting of vortex lines (5.32)**
27. `nb.core("C10", "Stretching and tilting: $\\frac{D\\omega_s}{Dt}=\\omega\\frac{\\partial u_s}{\\partial s}$, $\\frac{D\\omega_n}{Dt}=\\omega\\frac{\\partial u_n}{\\partial s}$, $\\frac{D\\omega_m}{Dt}=\\omega\\frac{\\partial u_m}{\\partial s}$ (5.32)",
    question="A figure skater spins faster by pulling in her arms; a tornado tightens as its column is stretched upward.
    How can a flow make vorticity stronger without any torque at all?")`
28. `nb.md` — **The problem in plain words:** "In three dimensions a flow can amplify the vorticity it carries simply by
    stretching vortex lines, and redirect it by tilting them. This is how turbulence builds its thin, intense vortices
    (Ch. 12) and why a flat, two-dimensional flow behaves so differently. It is also, with the planet's vorticity, how a
    converging column of air becomes a cyclone (C11)."
29. `nb.md` — **The idea:**
    ```
    (ω·∇)u = |ω| ∂u/∂s        (how u changes ALONG the vortex line)
       ├─ along the line  (∂u_s/∂s > 0): the line is stretched → thinner → spins faster   (the skater)
       └─ across the line (∂u_n/∂s, ∂u_m/∂s): the line is turned → vorticity appears in new directions (tilting)
    2-D flow: ω ⟂ plane, nothing varies along it → neither
    ```
30. `nb.primer("Frenet frame of a curve", "At each point of a smooth curve: the unit tangent e_s, the principal normal
    (toward the centre of curvature) and the binormal (their cross product). The book's Fig. 5.9 takes e_n **away** from
    the centre of curvature and e_m as the second normal; for a helix of radius a and pitch 2πc the curvature is
    a/(a² + c²) and the torsion c/(a² + c²).", code="d = ch05.helix_frame(0.0, a=1.0, c=0.3)            # frame at s = 0 of a helix\nprint(np.round([d['e_s'], d['e_n'], d['e_m']], 4))   # three orthonormal vectors\nprint(round(d['curvature'], 6), round(d['torsion'], 6), np.dot(d['e_s'], d['e_n']))   # 0.917431 0.275229 0.0")`
    — *expect:* e_n = (1, 0, 0) (away from the axis at s = 0), curvature 0.917431, torsion 0.275229, dot 0.0.
31. `nb.primer("angular momentum of a spinning cylinder", "A cylinder of mass m and radius R spinning at Ω about its axis
    has angular momentum L = IΩ with I = ½mR² = mA/2π (A its cross-section). With no torque L is fixed; stretch the cylinder
    at fixed mass and volume and A shrinks, so Ω grows like 1/A — like the length.", code="m, A0, Om0 = 1.0, 1e-4, 10.0
    # mass [kg], section [m²], spin [rad/s]\nfor stretch in (1, 2, 4):                            # length × 1, 2, 4 at fixed volume\n    A = A0/stretch; Om = Om0*A0/A                     # L = (m A/2π) Ω fixed\n    print(stretch, A, Om)                             # spin doubles when the length doubles")`
    — *expect:* `1 0.0001 10.0`, `2 5e-05 20.0`, `4 2.5e-05 40.0`.
32. `nb.derivation("D16", …)` — Part F D16 (6 steps), ref "5.31".
33. `nb.note` — **N38 [B]** "D16's result is (5.31):" equation
    `(\boldsymbol\omega\cdot\nabla)\mathbf u=\Big[\boldsymbol\omega\cdot\Big(\mathbf e_s\frac{\partial}{\partial s}+\mathbf e_n\frac{\partial}{\partial n}+\mathbf e_m\frac{\partial}{\partial m}\Big)\Big]\mathbf u=\omega\frac{\partial\mathbf u}{\partial s}`,
    ref "5.31". "The stretching part is (e_s·G e_s)ω along ω; the tilting part is the rest of Gω, perpendicular to ω."
34. `nb.derivation("D17", …)` — Part F D17 (5 steps), ref "5.32".
35. `nb.worked_example("stretching and tilting of one vortex line", "ω = (1, 0, 1) s⁻¹ in the axial strain
    G = diag(−½, −½, 1) s⁻¹ (u = −x/2, v = −y/2, w = z). 1. (ω·∇)u = Gω = (−0.5, 0, 1) s⁻². 2. e_s = ω/\|ω\| = (1, 0, 1)/√2;
    stretching rate e_s·G e_s = (−0.5 + 1)/2 = 0.25 s⁻¹. 3. Stretching part 0.25 × ω = (0.25, 0, 0.25) — along ω. 4.
    Tilting part (−0.5, 0, 1) − (0.25, 0, 0.25) = (−0.75, 0, 0.75) — perpendicular to ω (dot product −0.75 + 0.75 = 0): the
    line is being turned toward the stretching axis z.")`
36. `nb.code` — *code:* `print(ch05.stretching_tilting_split(np.array([1.0, 0, 1.0]), ch05.strain_preset("axial_stretch")))`
    · `for name, w0 in (("axial_stretch", [0, 0, 1]), ("shear_tilt", [1, 0, 0]), ("planar", [0, 0, 1])): print(name,
    ch05.uniform_strain_vorticity(w0, name, 2.0))` · `print(ch05.stretched_tube(2.0, 1.0, 10.0, 1e-4))`. *expect:* e_s
    (0.7071, 0, 0.7071), rate 0.25, stretching (0.25, 0, 0.25), tilting (−0.75, 0, 0.75), total (−0.5, 0, 1) · axial
    (0, 0, 7.389) · shear (1, 0, 2) · planar (0, 0, 1) · `{omega: 20.0, A: 5e-05, Gamma: 0.001}`. *explain:* 1. the split of
    the worked example; 2. three steady linear flows for 2 s: stretching multiplies ω_z by e^{αt}, shear tilts ω_x into
    ω_z = s t ω_x, the planar flow leaves ω ⟂ plane alone; 3. the stretched tube of D18.
37. `nb.check_agree` — **from scratch (curation §7):** the projection by hand: `w = np.array([1.0, 0, 1.0]); G =
    ch05.strain_preset("axial_stretch"); es = w/np.linalg.norm(w); rate = es @ G @ es; st = rate*w; ti = G @ w - st`;
    `assert np.allclose(st, S["stretching"]) and np.allclose(ti, S["tilting"]) and abs(ti @ w) < 1e-14`.
38. `nb.plotly` — **N53 [B] (Fig. 5.9 analogue)**: plotly 3-D: a helical vortex line (`ch05.helix`, a = 1, c = 0.3, three
    turns) with, at a chosen point, its frame e_s (teal), e_n (orange, away from the axis), e_m (blue) as cones, and the
    stretching (purple) and tilting (blue) arrows of the axial-strain G. A dropdown switches G among the three presets.
    *see:* "the frame rides the helix; the purple arrow lies along the line, the blue one across it." *read:* "(5.31):
    only the change of u along the line counts; its along-component stretches, its across-components tilt." *change:*
    "…the planar preset with ω along z: both arrows vanish."
39. `nb.animation` — (`player="frames"`, 24 frames, FAST 12) a vortex tube in axial strain α = 1 s⁻¹ over 2 s: the tube
    (drawn as a cylinder with a spinning stripe) lengthens by e², its radius shrinks by e⁻¹, the stripe turns faster;
    beside it, a straight vortex segment in the shear w = s x tilting from x toward z; a small panel with \|ω\|(t) (purple,
    e^{αt}, the ghost dashed) and the tilt angle (blue). *see / read / change:* "stop at t = ln 2 = 0.69 s: length × 2,
    vorticity × 2 — (D18); the tilted segment's ω_z grows linearly." *change:* "…α negative (compression): the tube
    fattens and slows down."
40. `nb.note` — **N21 [B]** "**Burgers' vortex** (the book leaves it to Exercise 5.12): an axial strain u_R = −αR/2,
    u_z = αz with a swirl" equation `u_\varphi=\frac{\Gamma}{2\pi R}\Big[1-e^{-\alpha R^2/4\nu}\Big],\qquad\omega_z=\frac{\alpha\Gamma}{4\pi\nu}e^{-\alpha R^2/4\nu}`, no ref, "is
    steady because stretching (which would thin and intensify the core) is balanced by diffusion (which would spread it).
    D18 derives the balance. Number: α = 1 s⁻¹, water, Γ = 10⁻³ m²/s: core radius √(4ν/α) = 2 mm, peak vorticity 79.6 s⁻¹.
    → the classic model of the finest vortices of turbulence (Ch. 12)."
41. `nb.derivation("D18", …)` — Part F D18 (9 steps), no ref (builder passes ref "5.13").
42. `nb.code` — *code:* `R = np.linspace(0, 8e-3, 401)` · `uR, uphi, uz, wz = ch05.burgers_vortex(R, 0.0, 1e-3, 1.0, 1e-6)` ·
    `print(ch05.burgers_core_radius(1.0, 1e-6), wz[0], np.trapezoid(wz*2*np.pi*R, R))` · `print(ch05.vorticity_budget_preset(
    "burgers", x=2e-3))`. *expect:* `0.002 79.5775 0.001000` (Γ recovered to 1e-5) · at R = 2 mm: stretching 29.27,
    advective 29.27, diffusion 0 (the inflection of the Gaussian: ∇²ω_z = 0 at R = core radius), residual ≈ 0. *explain:*
    1. the closed-form Burgers vortex; 2. core radius, peak and total circulation; 3. at R = √(4ν/α) diffusion changes sign:
    inside it removes vorticity, outside it adds.
43. `nb.figure` — Burgers ω_z(R) (teal) with the three term curves across R ∈ [0, 8] mm: advective (grey), stretching
    (purple), diffusion (rose), their residual (black ≈ 0). *see / read / change:* "stretching is largest on the axis;
    diffusion is negative inside the core and positive outside; the three add to zero everywhere — the vortex is steady."
    *change:* "…α four times larger: the core halves (√(4ν/α)) and the peak quadruples."
44. `nb.explainer("vorticity_stretching_tilting", heading="How can a flow spin fluid faster without a torque?", why="A
    vortex line and a material element ride a chosen flow on one clock while \|ω\| grows or turns, the stretching and
    tilting arrows update and the term bars of (5.13) fight; switching on ν brings in diffusion and the Burgers balance —
    and the Derivation tab carries the ★★★ D09.", tries=["Axial strain, press ▶: \|ω\| grows as e^{αt}; the material element
    stays parallel to ω.", "Shear tilt: a new component appears — tilting.", "Planar flow: nothing happens — 2-D.",
    "Burgers: raise ν and watch the core widen to √(4ν/α)."])`
45. `nb.md` — "> ⚠️ **Common confusion:** 'stretching needs a torque.' It needs none: angular momentum is conserved while
    the moment of inertia shrinks — the spin rises on its own, exactly as the skater's does."
46. `nb.md` — **What would change if…** "…the frame rotated? Then a vertical fluid column contains the planet's vorticity
    2Ω even when it does not spin relative to the Earth — and stretching it spins it up (C11)."

**Recap for C11**
47. `nb.recap("R18", "Planetary and absolute vorticity", "Seen from a frame rotating at Ω, a fluid at rest in the frame
    has inertial vorticity 2Ω — the *planetary vorticity*; the *absolute vorticity* is ω + 2Ω (Ch. 3: ω′ = ω − 2Ω). On the
    Earth only the local vertical part matters for horizontal motion: f = 2Ω sin φ (Ch. 4 P126), and the absolute
    vertical vorticity is ζ + f.", where="Ch. 3 §3.4, Ch. 4 §4.7")` + `nb.code`: `print(ch05.absolute_vorticity(np.array([0,
    0, 1e-5]), np.array([0, 0, rotating.OMEGA_EARTH])))` · `print([rotating.coriolis_parameter(np.deg2rad(p)) for p in (0, 30,
    45, 60, 90)])` · round trip `kinematics.vorticity_in_rotating_frame(ch05.absolute_vorticity(w, Om), Om) == w`. *expect:*
    `[0, 0, 1.5584e-4]` · `[0, 7.2921e-5, 1.0313e-4, 1.2630e-4, 1.4584e-4]` · True.

**C11 — Kelvin in a rotating frame: absolute circulation (5.33) and the fluid column**
48. `nb.core("C11", "Absolute circulation: $\\frac{D\\Gamma_a}{Dt}=0,\\ \\Gamma_a=\\Gamma+2\\int_A\\boldsymbol\\Omega\\cdot\\mathbf n\\,dA$ (5.33)",
    question="Air crosses a mountain range and is squashed into a thinner layer; a ring of air drifts toward the pole. On a
    spinning planet, what is conserved — and which way must each turn?")`
49. `nb.md` — **The problem in plain words:** "Kelvin's theorem failed in a rotating frame because the Coriolis force is
    not conservative (C03). Something must replace it — otherwise the atmosphere would have no rule for how cyclones spin up
    over plains and anticyclones sit over mountains. The replacement is the circulation of the *absolute* vorticity, and
    for a column it becomes (ζ + f)/h = const: potential vorticity, the single most useful conserved quantity of Ch. 13."
50. `nb.md` — **The idea:**
    ```
    inertial frame:  Γ conserved (Kelvin)
    rotating frame:  Γ changes, because the Coriolis loop integral = −2Ω · (rate of the loop's vector area)
                     →  Γ + 2Ω·A_vec  =  ∫(ω + 2Ω)·n dA  =  Γ_a  is conserved            (5.33)
    thin column:     (ζ + f) A = const  and  A h = const   →   (ζ + f)/h = const      stretch → cyclonic, squash → anticyclonic
    ```
51. `nb.md` — reminder: "**Loop vector area** (primed in C05): A_vec = ½∮x × dx, the loop's area times its normal; for a
    horizontal loop on the Earth only its vertical part meets the vertical planetary vorticity." (not a second primer)
52. `nb.derivation("D19", …)` — Part F D19 (10 steps), ref "5.33".
53. `nb.md` — gloss: "**Column mass.** A vertical column of incompressible fluid keeps its volume: section × height =
    const, so a taller column is a thinner one."
54. `nb.derivation("D20", …)` — Part F D20 (6 steps), no ref (builder passes ref "5.33").
55. `nb.note` — **N39 [B]** "**Planetary stretching and tilting.** With Ω = Ωe_z, 2(Ω·∇)u = 2Ω∂u/∂z, so, keeping only
    this term of (5.30)," equation `\frac{D\omega_z}{Dt}=2\Omega\frac{\partial w}{\partial z},\qquad\frac{D\omega_x}{Dt}=2\Omega\frac{\partial u}{\partial z},\qquad\frac{D\omega_y}{Dt}=2\Omega\frac{\partial v}{\partial z}`,
    no ref. "Stretching vertical *fluid* lines (not vortex lines — they contain 2Ω already) raises ω_z; tilting them makes
    horizontal vorticity. By continuity ∂w/∂z = −(∂u/∂x + ∂v/∂y): horizontal convergence spins fluid up. This is D20
    step 6's small-change form." + `nb.code`: `G = np.array([[-0.1, 0, 0], [0, -0.1, 0], [0, 0, 0.2]])*1e-4` (convergence
    stretching w) · `print(ch05.planetary_vorticity_terms(G, np.array([0, 0, rotating.OMEGA_EARTH])))` → *expect:*
    `[0, 0, 2.9168e-09]` (2Ω × 2×10⁻⁵).
56. `nb.worked_example("a column over a ridge and a ring moving poleward", "1. **Column** at 45° N (f = 1.03×10⁻⁴ s⁻¹),
    at rest (ζ₀ = 0), 1000 m deep. It moves over a 200 m high ridge: h = 800 m. ζ = (0 + f)(800/1000) − f = −0.2f =
    −2.06×10⁻⁵ s⁻¹: anticyclonic (clockwise in the north). Past the ridge, h = 1000 m again and ζ = 0. Into a trough
    1200 m deep: ζ = +0.2f, cyclonic. 2. **Ring of air** of radius 500 km (A = 7.85×10¹¹ m²) at rest at 30° N, carried to
    60° N keeping its area: Γ₁ = 2ΩA(sin 30° − sin 60°) = 1.458×10⁻⁴ × 7.854×10¹¹ × (−0.366) = −4.19×10⁷ m²/s; mean
    ζ = Γ₁/A = −5.34×10⁻⁵ s⁻¹ — anticyclonic, half the local f.")`
57. `nb.code` — *code:* `f45 = rotating.coriolis_parameter(np.deg2rad(45.0))` · `print([ch05.column_relative_vorticity(h,
    1000.0, 0.0, f45) for h in (800.0, 1000.0, 1200.0)])` · `A = np.pi*5e5**2; print(ch05.relative_circulation_after_move(0.0,
    A, 30.0, A, 60.0), ch05.relative_circulation_after_move(0.0, A, 30.0, A, 60.0)/A)` · the rotating scenario of E3:
    `sc = ch05.kelvin_scenario("rotating"); t = np.linspace(0, 5, 6); P = ch05.material_loop(sc["u"], sc["pts0"], t);
    print([ch05.absolute_circulation(sc["u"], Pk, np.array([0, 0, 0.5]), t=tk) for Pk, tk in zip(P, t)])`. *expect:*
    `[-2.0625e-05, 0.0, 2.0625e-05]` · `-4.19261e+07 -5.33820e-05` · pairs (Γ, Γ_a): (0, 3.14159), …, (1.98587, 3.14159)
    — Γ_a constant to 1e-8. *explain:* 1. the column rule of D20 across a ridge and a trough; 2. the poleward ring of
    (5.33); 3. a material loop in a converging rotating flow: Γ grows, Γ_a does not move.
58. `nb.check_agree` — **from scratch (curation §7):** the vector area ½∮x × dx by the trapezoid rule on the loop's 256
    points (np.roll differences, as in the primer) vs `ch05.loop_vector_area` → `assert np.allclose(A_hand, A_lib,
    rtol=1e-4)`; then Γ_a = Γ + 2Ω·A_vec by hand from `ch05.loop_circulation` → `assert np.allclose(Ga_hand, np.pi,
    rtol=1e-6)`.
59. `nb.figure` — **N54 [B] (Fig. 5.10 analogue)**: top, a layer over a sloping and ridged bottom (depth h(x), 1000 → 800
    → 1200 m) with columns drawn at five places, each with a small curved arrow whose size and sense show ζ
    (`ch05.column_over_slope(x, "ridge", 45.0)`); bottom, ζ(x)/f (teal) and h(x)/h₀ (grey), and the conserved (ζ + f)/h (amber,
    flat). *see:* "the columns spin clockwise over the ridge and counterclockwise over the trough; the amber line is flat."
    *read:* "(ζ + f)/h is carried unchanged by each column (D20): squashed columns lose absolute vorticity, stretched ones
    gain it." *change:* "…the same layer at 10° N: f is 2.4 times smaller, so every ζ is 2.4 times smaller — the equator
    has little planetary vorticity to lend."
60. `nb.plotly` — `slider_figure` over latitude 0…90° (19 steps, FAST 10): ζ(h) for a column starting at rest with
    h₀ = 1000 m (teal line through the origin, slope f/h₀) and the relative circulation Γ₁(φ₁) of a 500 km ring starting at
    rest at 30° N (second trace, right axis). Title "The planet lends more spin toward the pole". *see / read / change:*
    "the ζ(h) line steepens with latitude (f = 2Ω sin φ); the ring's Γ₁ is zero at 30° N and grows more negative
    poleward." *change:* "…the Southern Hemisphere (φ < 0): f < 0, the signs flip — cyclones turn clockwise there."
61. `nb.explainer("vorticity_equation_rotating", heading="What does a spinning planet add to the vorticity equation?",
    why="Dragging a column over a ridge, moving a ring of air in latitude and switching the scene of a term-bar budget of
    (5.30) shows the same conservation law — and the same equation — in three regimes, with the ★★★ D15 in its Derivation
    tab.", tries=["Column mode: drag h below h₀ and watch the column turn anticyclonic while (ζ + f)/h stays put.", "Ring
    mode: move the ring from 30° N to 60° N and read Γ = −4.19×10⁷ m²/s.", "Budget mode: switch scenes and watch which bar
    dominates.", "Southern Hemisphere preset: every sign flips."])`
62. `nb.md` — "> ⚠️ **Common confusion:** '2ΩA is the planetary part.' It is 2Ω·A_vec — only the component of Ω along the
    loop's normal counts. On the Earth, for a horizontal loop, that is f A = 2Ω sin φ A."
63. `nb.md` — **What would change if…** "…the layer were also stratified? The column would be a slab between two density
    surfaces and h its thickness: (ζ + f)/h becomes Ertel's potential vorticity (Ch. 13) — the same bookkeeping."

### A.7 §5.7 Interaction of Vortices — C12 (+N40 N41 N55 N56 · D21), C13 (+N42 N43 N57 N58 N59 · D22 · E8)
1. `nb.section("5.7", "Interaction of Vortices", intro="**What is this section about?** Idealise each vortex as a line
   vortex. A straight line vortex cannot move itself, but Helmholtz says it moves with the fluid — and the fluid at its
   position is moved by all the *other* vortices (Biot–Savart). That one rule makes like vortices orbit, opposite ones
   travel together, a vortex slide along a wall (its mirror image does the pushing) and smoke rings leap-frog.")`

**C12 — Point vortices move each other: pairs and the centre of vorticity**
2. `nb.core("C12", "Point vortices: $\\frac{d\\mathbf x_k}{dt}=\\sum_{j\\ne k}\\frac{\\Gamma_j}{2\\pi}\\frac{\\mathbf e_z\\times(\\mathbf x_k-\\mathbf x_j)}{\\lvert\\mathbf x_k-\\mathbf x_j\\rvert^2}$",
   question="Two hurricanes that come within a thousand kilometres start to circle each other; two dimples made by a
   paddle stroke glide away side by side. A vortex cannot push itself — so how do vortices move?")`
3. `nb.md` — **The problem in plain words:** "Weather forecasters watch binary cyclones dance (the Fujiwhara effect);
   airport controllers space landing aircraft because each plane leaves a pair of vortices that sinks behind it;
   modellers of 2-D turbulence and of ocean eddies use swarms of point vortices. All of these rest on one simple rule."
4. `nb.md` — **The idea:**
   ```
   real vortex (core with vorticity)  ≈  line vortex of the same Γ     (outside the core the flow is irrotational)
   each vortex k moves with the velocity induced at x_k by all the OTHERS:  dx_k/dt = Σ_{j≠k} (Γ_j/2π) e_z×(x_k−x_j)/|x_k−x_j|²
   two alike  →  orbit their centre of vorticity         two opposite  →  translate together at Γ/2πh
   ```
5. `nb.md` — glosses: "**Superposition.** Outside the cores the flow is irrotational and (Ch. 6) obeys the linear Laplace
   equation, so velocities of several vortices add. **Lever rule.** The balance point of two weights w₁, w₂ a distance h
   apart is h w₂/(w₁ + w₂) from the first. **Circular motion.** A point turning at rate θ̇ at distance d from the centre
   moves at d θ̇."
6. `nb.primer("systems of ODEs for interacting bodies and their invariants", "N vortices give 2N coupled ODEs (each
   position's rate depends on all the others). `solve_ivp` integrates them as one state vector. Because the exact motion
   conserves some quantities — here ΣΓx, ΣΓ\|x\|² and the energy H — watching them stay constant tells us the numerical
   solution is trustworthy.", code="xv0 = np.array([[-0.5, 0.5], [0.0, 0.0]]); G = np.array([1.0, 1.0])   # two like vortices\nt = np.linspace(0, 20, 5)\nX = ch05.point_vortex_evolve(xv0, G, t)                               # (T, 2, M) positions\nprint([ch05.point_vortex_invariants(Xk, G)['I'] for Xk in X])        # ΣΓ|x|² = 0.5 every time")`
   — *expect:* `[0.5, 0.5, 0.5, 0.5, 0.5]` (to 1e-10).
7. `nb.derivation("D21", …)` — Part F D21 (6 steps), no ref (builder passes ref "5.2").
8. `nb.note` — **N40 [B] (Fig. 5.11)** "**Two like vortices** Γ₁, Γ₂ > 0 a distance h apart push each other sideways at"
   equation `V_1=\Gamma_1/2\pi h,\qquad V_2=\Gamma_2/2\pi h`, no ref, "in opposite directions, so the pair turns about its centre of
   vorticity G, at h₁ = Γ₂h/(Γ₁ + Γ₂) from vortex 1, at the rate (Γ₁ + Γ₂)/2πh² (D21; the book leaves G to Exercise 5.18).
   Numbers: Γ₁ = Γ₂ = 1 m²/s, h = 1 m → each moves at 0.159 m/s, period 19.7 s; Γ₂ = 3Γ₁ → G is ¾h from vortex 1. ⚠️ The
   caption of Fig. 5.11 calls G the point where the induced velocity is zero — true only for equal strengths."
9. `nb.note` — **N41 [B] (Figs. 5.12–5.13)** "**Equal and opposite pair**: each moves the other at Γ/2πh in the same
   direction, so the pair translates at" equation `V=\Gamma/(2\pi h)`, no ref, "relative to the fluid — made by a paddle stroke
   or a knife blade drawn briefly through a bucket (Lighthill's picture, our simulation in C13). The wake of an aircraft
   is such a pair, sinking behind it (Ch. 14)."
10. `nb.worked_example("two like vortices, then two opposite", "1. Γ₁ = Γ₂ = 1 m²/s at (−0.5, 0) and (0.5, 0) m: vortex 2 is
    pushed up at 1/(2π × 1) = 0.159 m/s, vortex 1 down at 0.159 m/s; G is the midpoint; the joining line turns at
    (0.159 + 0.159)/1 = 0.318 rad/s — one orbit in 2π/0.318 = 19.7 s. 2. Γ₂ = 3 m²/s instead: V₁ = 0.159, V₂ = 0.477 m/s;
    G at 0.75 m from vortex 1 (the lever rule with weights 1 and 3); rate (1 + 3)/(2π) = 0.637 rad/s. 3. Γ₂ = −1 m²/s:
    both move up at 0.159 m/s — a straight march.")`
11. `nb.code` — *code:* `for G2 in (1.0, 3.0, -1.0): print(G2, ch05.vortex_pair(1.0, G2, 1.0))` · `P =
    ch05.point_vortex_preset("equal_pair"); t = np.linspace(0, 19.7392, 200); X = ch05.point_vortex_evolve(P["xv"],
    P["Gamma"], t); print(np.abs(X[-1] - X[0]).max())` (back to the start after one period) · `inv = [ch05.
    point_vortex_invariants(Xk, P["Gamma"]) for Xk in X]; print(np.ptp([v["I"] for v in inv]), np.ptp([v["H"] for v in
    inv]))` · `print(ch05.centre_of_vorticity(ch05.point_vortex_preset("unequal_pair")["xv"], [1.0, 3.0]))`. *expect:* the
    three dicts of the worked example (V1 0.159155, V2 0.159155/0.477465/0.159155, centre_from_1 0.5/0.75/inf, rate
    0.318310/0.636620/0, period 19.739/9.870/inf, translation 0/0/0.159155) · ≤ 1e-8 · ≤ 1e-12, ≤ 1e-12 · `[0.25, 0.0]`
    (¾ of the way from −0.5 to 0.5). *explain:* 1. the closed forms of D21; 2. the full point-vortex integration returns
    to its start after the predicted period; 3. the invariants stay flat — the integrator is trustworthy; 4. the centre of
    vorticity of the unequal pair.
12. `nb.check_agree` — **from scratch (curation §7):** dx_k/dt by a double loop over pairs: `v = np.zeros((2, M)); for k
    in range(M): for j in range(M): if j != k: d = xv[:, k] - xv[:, j]; v[:, k] += G[j]/(2*np.pi)*np.array([-d[1],
    d[0]])/(d @ d)` · `assert np.allclose(v, ch05.point_vortex_rhs(xv, G), atol=1e-14)` for the three-vortex preset.
    Markdown: "e_z × (a, b) = (−b, a): each pair of vortices, one line."
13. `nb.animation` — **N55, N56 [B] (Figs. 5.11, 5.12 analogues)** (`player="video"`, 60 frames, FAST 30): three panels
    on one clock — equal pair orbiting G (teal dots, trails), Γ₂ = 3Γ₁ orbiting its off-centre G, opposite pair (teal and
    rose) marching up; G marked with a black cross; the velocity arrows of each vortex. *see:* "two circles of different
    sizes round a fixed cross; a pair of dots moving straight up." *read:* "every vortex is carried by the others: like
    ones orbit their centre of vorticity (which never moves), opposite ones travel together at Γ/2πh." *change:* "…h
    halved: the orbits are four times faster ((Γ₁ + Γ₂)/2πh²) and the march twice as fast (Γ/2πh)."
14. `nb.figure` — trajectories of the three presets over 20 s plus the three-vortex preset (Γ = 1, 1, −0.5), with the
    invariants P, I, H vs t in a small panel (flat lines). *see / read / change* in the same words as the animation;
    change: "…a fourth vortex added: the motion can become chaotic, but the invariants stay flat."
15. `nb.md` — **What would change if…** "…one of the two vortices were replaced by a wall? The wall acts like a mirror
    vortex of opposite sign — the pair becomes a vortex and its image (C13)."

**C13 — The method of images: a vortex near a wall**
16. `nb.core("C13", "Images: a vortex drifts along a wall at $\\Gamma/4\\pi h$", question="A vortex near the bottom of a
    tank slides sideways along it; a smoke ring blown at a wall grows wider as it slows. How does a wall push a vortex
    without touching it?")`
17. `nb.md` — **The problem in plain words:** "A wall forbids flow through it. For vortices this boundary condition has a
    beautiful solution: pretend the wall is a mirror and put an opposite vortex behind it. The mirror twin then moves the
    real vortex — and the same trick handles circular tanks, channels, the ground under a landing aircraft (Ch. 14) and
    cylinders in a stream (Ch. 6)."
18. `nb.md` — **The idea:**
    ```
         vortex A (+Γ) •                      wall = the line where A's push and its twin's push have
    ───────────────────────── wall            equal and opposite normal parts  →  no flow through it
         image B (−Γ)  •                      A is moved only by B (distance 2h): V = Γ/(2π·2h) = Γ/4πh, along the wall
    ```
19. `nb.md` — gloss: "**Mirror symmetry.** Reflect a counterclockwise vortex in the wall and reverse its spin: at every wall
    point the two velocity vectors are mirror images, so their normal parts cancel and their tangential parts add."
20. `nb.derivation("D22", …)` — Part F D22 (5 steps), no ref (builder passes ref "5.2").
21. `nb.worked_example("a vortex half a metre above the floor", "Γ = 1 m²/s at height h = 0.5 m. 1. Image: −1 m²/s at
    −0.5 m. 2. At the wall point directly below, both vortices are 0.5 m away; each gives 1/(2π × 0.5) = 0.318 m/s, both
    horizontal and in the same direction: 0.637 m/s along the floor, 0 through it. 3. At the wall point 0.5 m to the side,
    distance √0.5 = 0.707 m: each gives 0.225 m/s at 45°; the vertical parts cancel, the horizontal ones add to 0.318 m/s.
    4. The vortex itself moves at 1/(4π × 0.5) = 0.159 m/s along the wall. 5. In a channel 1 m high with the vortex at
    0.25 m, all the images add up to (Γ/4H)cot(πh/H) = 0.25 m/s — slower than the single-wall 0.318 m/s, because the
    ceiling's images push back.")`
22. `nb.code` — *code:* `xa, Ga = ch05.wall_image_system(np.array([[0.0], [0.5]]), np.array([1.0]))` · `xw = np.array([np.
    linspace(-3, 3, 200), np.zeros(200)])` · `vel = ch05.point_vortex_velocity(xw, xa, Ga); print(np.abs(vel[1]).max(),
    vel[0][[100]])` · `print(ch05.vortex_near_wall_speed(1.0, 0.5), ch05.point_vortex_rhs(np.array([[0.0], [0.5]]),
    np.array([1.0]), boundary="wall"))` · `print(ch05.channel_image_velocity(0.25, 1.0, 1.0), 1/(4*1.0)/np.tan(np.pi*0.25))`.
    *expect:* max \|v\| on the wall ≤ 1e-12; the wall speed at x ≈ 0 near 0.637 m/s · `0.159155 [[0.159155], [0.0]]` ·
    `0.25 0.25`. *explain:* 1. the vortex and its image; 2. the velocity at 200 wall points has no normal part — the wall
    is a streamline; 3. the vortex's own drift, from the image only; 4. the channel's infinite image series vs its closed
    form.
23. `nb.check_agree` — **from scratch (curation §7):** the image sum by hand at the vortex: `d = np.array([0.0, 0.5]) -
    np.array([0.0, -0.5])`; `v = -1.0/(2*np.pi)*np.array([-d[1], d[0]])/(d @ d)`; `assert np.allclose(v, [1/(4*np.pi*0.5),
    0.0])`; and on the wall at x = 1 m: `assert abs(v_wall_A[1] + v_wall_B[1]) < 1e-15`.
24. `nb.figure` — **N58 [B] (Fig. 5.14 analogue)**: the vortex A (teal), its image B (rose, faded below the wall), and at
    9 wall points the two velocity arrows (thin teal/rose) and their sum (black, lying along the wall). *see:* "every black
    arrow lies flat on the wall." *read:* "the image makes the wall a streamline; the black arrows are the slip velocity an
    inviscid fluid is allowed along a wall." *change:* "…an image of the same sign: the black arrows would point straight
    into the wall — no longer a wall."
25. `nb.note` — **N57 [B] (Fig. 5.13)** "**The knife-blade pair in a bucket** (after Lighthill): a knife drawn briefly
    through a bucket leaves an opposite pair that marches off; nearing the wall the two vortices separate and run along
    it in opposite directions, each pushed by its own image. For a circle of radius a the image of a vortex at x is at
    the inverse point a²x/\|x\|² (strength −Γ) — below, our simulation."
26. `nb.animation` — (`player="video"`, 60 frames, FAST 30) the knife pair `ch05.point_vortex_preset("knife_bucket")` in a
    circle of radius 1 m (images by `circle_image_system(inside=True)`), trails and the images drawn faintly outside the
    circle. *see:* "the pair crosses the bucket, splits at the wall and the two vortices run round it in opposite
    directions." *read:* "far from the wall the pair translates at Γ/2πh (C12); near it each vortex is dragged by its own
    image — the wall drift Γ/4πh." *change:* "…a stronger stroke (larger Γ): the same path, faster."
27. `nb.note` — **N42 [B] (Fig. 5.15)** "**A smoke ring meeting a wall.** A ring moves itself (its own curvature, the
    thin-core self-speed below — Kelvin's formula, cited) and is moved by its image ring behind the wall. The image pushes
    it outward and holds it back: the ring **widens** and **slows** as it approaches, never reaching the wall. We compute
    the path with the velocity one circular filament induces at another (elliptic integrals — primer below) and the
    self-speed" equation `U=\frac{\Gamma}{4\pi R}\Big[\ln\frac{8R}{a}-\frac14\Big]`, no ref. "Number: Γ = 1 m²/s, R = 1 m, core a = 0.1 m:
    U = 0.329 m/s."
28. `nb.primer("complete elliptic integrals with the parameter m", "The velocity induced by a circular filament involves
    K(m) = ∫₀^{π/2} dθ/√(1 − m sin²θ) and E(m) = ∫₀^{π/2} √(1 − m sin²θ) dθ. ⚠️ `scipy.special.ellipk(m)` and `ellipe(m)` take
    the **parameter m = k²**, not the modulus k — passing k is a classic silent bug (K → ∞ as m → 1: a point near the
    filament).", code="from scipy.special import ellipk, ellipe\nfrom scipy.integrate import quad\nm = 0.5
    # the parameter m = k²\nK = quad(lambda th: 1/np.sqrt(1 - m*np.sin(th)**2), 0, np.pi/2)[0]\nprint(ellipk(m), K)                          # 1.854075 1.854075: same definition\nprint(ellipk(0.5**0.5))                     # 2.085974 — what you get by passing k instead of m")`
    — *expect:* `1.854075 1.854075` and `2.085974` (a different number: the bug is silent).
29. `nb.code` — *code:* `print(ch05.ring_self_velocity(1.0, 0.1, 1.0))` · `print(ch05.ring_ring_velocity(0.0, 0.5, 0.5,
    0.0, 1.0), ch05.ring_axis_velocity(0.5, 0.5, 1.0))` (on the axis: the elliptic formula meets the closed form) · `res =
    ch05.ring_dynamics([dict(R=0.5, z=0.0, Gamma=1.0, a=0.05)], np.linspace(0, 30, 301), wall_z=2.0)` (a ring with Γ > 0
    moves toward +z, here toward a wall at z = 2 m) · `print(res["R"][[0, 100, 200, 300], 0], res["z"][[0, 100, 200,
    300], 0])` · `print(np.all(np.diff(res["R"][:, 0]) > 0), np.all(np.diff(np.gradient(res["z"][:, 0])) < 0))`. *expect:* `0.328816` · u_R ≈ 0, u_z =
    0.353553 twice · R increases monotonically, z climbs toward (but never reaches) the wall at 2 m, the approach speed
    decreases monotonically: `True True`. *explain:* 1. Kelvin's self-speed; 2. the ring-induced velocity on the axis equals
    ΓR²/2(R² + z²)^{3/2}; 3. a ring and its image: widening and slowing — Fig. 5.15b in numbers.
30. `nb.figure` — **N59 [B] (Fig. 5.15 analogue)**: the ring's cross-section path (R(t), z(t)) approaching the wall, with
    its image path mirrored below (faded), and the approach speed −dz/dt vs z (inset). *see / read / change:* "the path
    bends outward along the wall; the speed falls toward zero." *change:* "…a fatter core (a = 0.2 m): a slower ring
    (ln(8R/a) smaller) that widens the same way."
31. `nb.note` — **N43 [B]** "**Leap-frogging rings.** Two coaxial rings of the same sense (the book leaves it to Exercise
    5.15): the front ring is widened and slowed by the rear one, the rear ring narrowed and sped up; it passes through the
    front ring and the roles swap — for ever in an ideal fluid. The total impulse ΣΓπR² is conserved." + `nb.animation`
    (`player="video"`, 60 frames, FAST 30): the meridional plane with both rings' cross-sections (teal) and trails over two
    passes, R₁, R₂ vs t below, the impulse (flat). *see / read / change:* "the rings take turns passing through each
    other; the impulse line is flat." *change:* "…unequal strengths: the passes become irregular and may stop." + `nb.code`:
    `res2 = ch05.ring_dynamics([dict(R=1.0, z=0.0, Gamma=1.0, a=0.1), dict(R=1.0, z=0.6, Gamma=1.0, a=0.1)], np.linspace(0,
    40, 401))`; `print(np.ptp(res2["impulse"])/res2["impulse"][0], np.sum(np.diff(np.sign(res2["z"][:, 1] - res2["z"][:,
    0])) != 0))` → *expect:* ≤ 1e-9 and ≥ 2 pass-throughs.
32. `nb.explainer("point_vortex_lab", heading="A vortex cannot push itself — so how do vortices move?", why="Dropping
    vortices of either sign, dragging them and pressing play turns D21 and D22 into motion: the orbit about a fixed centre
    of vorticity, the marching pair, the image that moves in lock-step behind a wall — with the invariants as a live
    accuracy check.", tries=["Equal pair: time one orbit (19.7 s) and watch G stay put.", "Set Γ₂/Γ₁ = 3: G moves to ¾ of
    the way; the orbit speeds up.", "Opposite pair near a wall: they march, then split along it.", "Click a vortex: the
    inspector adds up the others' pushes."])`
33. `nb.md` — **What would change if…** "…instead of a few vortices we had infinitely many, packed side by side along a
    line? Their pushes would add up to a jump in velocity across the line — a vortex sheet (C14)."

### A.8 §5.8 Vortex Sheet — C14 (+N44 N60 · D23 · E9)
1. `nb.section("5.8", "Vortex Sheet", intro="**What is this section about?** Line up infinitely many line vortices side by
   side and the velocity jumps across the line: above it the fluid moves one way, below it the other. The strength of
   such a vortex sheet — its circulation per unit length — is exactly that jump. Sheets are the idealisation of every
   shear layer, wake and wing.")`

**C14 — The vortex sheet: strength = jump in tangential velocity**
2. `nb.core("C14", "A vortex sheet: $\\gamma=d\\Gamma/ds=u_2-u_1$", question="Where a fast stream runs beside a slow one — a
   river entering a lake, the jet stream over calmer air, the air leaving a wing's trailing edge — the velocity changes
   almost abruptly. What is a velocity jump made of?")`
3. `nb.md` — **The problem in plain words:** "A thin layer across which the velocity jumps looks like a boundary, not like
   vorticity. Yet all the spin of the flow is concentrated there: the jump *is* a sheet of vortices. Seeing it this way
   explains why such layers roll up into billows (Kelvin–Helmholtz, Ch. 11) and how a wing's lift is carried by vortices
   (Ch. 14)."
4. `nb.md` — **The idea:**
   ```
   u₁ = −γ/2    ← ← ← ← ← ← ← ←        just above
   ⊙ ⊙ ⊙ ⊙ ⊙ ⊙ ⊙ ⊙ ⊙ ⊙ ⊙ ⊙         a row of counterclockwise filaments, γ ds each
   u₂ = +γ/2    → → → → → → → →        just below
   circulation of a thin box ds × dn round it:  dΓ = (u₂ − u₁) ds   →   strength γ = u₂ − u₁
   ```
   "**A velocity jump is a sheet of vortices, and the jump is its strength.**"
5. `nb.note` — **N44 [B]** "**A sheet of line vortices** (Fig. 5.16): infinitely many parallel filaments on a surface. The
   tangential velocity jumps across it; the normal velocity is continuous (zero if the sheet does not move across
   itself); a real sheet has a finite thickness that spreads the jump (N22)." equation `[u_t]\neq0,\qquad[v_n]=0`, no ref.
   "Number: filaments of total strength γ = 2 m/s per metre → far above u = −1 m/s, far below +1 m/s."
6. `nb.derivation("D23", …)` — Part F D23 (5 steps), no ref (builder passes ref "3.18").
7. `nb.md` — "> ⚠️ **Common confusion — two signs and two Γ's.** The text writes $\Gamma\equiv\frac{d\Gamma}{ds}=u_2-u_1$ (counterclockwise
   circulation per length, u₁ above, u₂ below); Fig. 5.16's caption prints $d\Gamma/ds=u_1-u_2$ — the clockwise sense, which
   is positive for its clockwise filaments. And the book reuses Γ for a quantity in m/s; we write γ. Our code:
   `vortex_sheet_strength(u_above, u_below, convention="ccw")` = u₂ − u₁ (switch `"caption"` for u₁ − u₂)."
8. `nb.md` — gloss: "**From a sum to a sheet.** N filaments of strength γL/N spaced L/N apart become a continuous sheet as
   N → ∞ — a Riemann sum turning into an integral (Ch. 2); far from the sheet (more than a spacing away) the row already
   looks continuous."
9. `nb.worked_example("a row of 10 filaments, then 100", "γ = 2 m/s on a sheet 1 m long, probe at x = 0. 1. N = 10:
   each filament has Γ = γ × 0.1 = 0.2 m²/s; at y = 5 cm the sum gives −0.854 m/s (the continuous finite sheet: −0.937 m/s).
   2. N = 100 (Γ = 0.02 m²/s each): −0.936551 m/s, equal to the continuous value to 2×10⁻⁶. 3. Just above the middle of
   the continuous sheet: −γ/2 = −1 m/s; just below +1 m/s; jump u₂ − u₁ = 2 m/s = γ ✓. 4. Far above a finite sheet the
   velocity fades (the sheet looks like a single vortex of circulation γL = 2 m²/s).")`
10. `nb.code` — *code:* `print(ch05.vortex_sheet_strength(-1.0, 1.0), ch05.vortex_sheet_strength(-1.0, 1.0,
    convention="caption"))` · `for N in (10, 100): print(N, ch05.discrete_sheet_u(0.0, 0.05, 2.0, N))` ·
    `print(ch05.continuous_sheet_velocity(0.0, 0.05, 2.0), ch05.continuous_sheet_velocity(0.0, 1e-6, 2.0),
    ch05.continuous_sheet_velocity(0.0, -1e-6, 2.0))` · `c = ch05.discrete_sheet_convergence(2.0); print(c["l1_error"],
    observed_order(1/np.array(c["N"]), c["l1_error"]))`. *expect:* `2.0 -2.0` · `10 -0.853907`, `100 -0.936551` ·
    `(-0.936549, 0.0)`, `(-1.0, 0)`, `(1.0, 0)` (to 1e-6) · `[0.043968, 0.0044131, 0.00044128]` and order 1.00.
    *explain:* 1. both sign conventions; 2. the discrete row at one height; 3. the continuous finite sheet: ∓γ/2 just above
    and below the middle; 4. the L1 error of the discrete profile falls as 1/N.
11. `nb.check_agree` — **from scratch (curation §7):** `N, L, g = 100, 1.0, 2.0; xs = (np.arange(N) + 0.5)*L/N - L/2; y =
    0.05; u = np.sum(g*L/N/(2*np.pi)*(-y)/(xs**2 + y**2))` · `assert np.allclose(u, ch05.vortex_sheet_velocity(np.array([0.0,
    0.05]), np.array([xs, 0*xs]), 2.0, L/N)[0], rtol=1e-12)`. Markdown: "A hundred line vortices, one line of numpy."
12. `nb.plotly` — **N60 [B] (Fig. 5.16 analogue)** `slider_figure` over N = 4…1000 (12 steps, log-spaced): u(y) at x = 0
    (teal) for the discrete row vs the continuous sheet (black) for y ∈ [−0.2, 0.2] m, and the filaments as dots on the
    x-axis inset; the dn × ds circuit drawn in the inset. Title "From a row of vortices to a velocity jump". *see:* "a
    wiggly profile that sharpens into a clean step from +1 to −1 m/s." *read:* "the jump across the sheet is its strength
    γ; within about one filament spacing of the row the discreteness shows, beyond it the row is a sheet." *change:*
    "…γ doubled: the step doubles; the convergence with N is unchanged."
13. `nb.animation` — (`player="video"`, 60 frames, FAST 30) a periodic sheet of period 1 m, N = 200 (FAST 100) point
    vortices with Krasny smoothing δ = 0.05, initial ripple 0.01 m, rolling up over t = 0…2 s (γ = 1 m/s)
    (`ch05.sheet_rollup`): the sheet (teal line through the points) wraps into a spiral; a text line with the total
    circulation (constant). *see:* "the ripple steepens, the sheet bunches at one point and winds into a spiral." *read:*
    "the sheet's own induced velocity carries its vortices toward the places where they are already crowded — the
    Kelvin–Helmholtz instability (Ch. 11); the smoothing δ only stands in for a finite thickness." *change:* "…a smaller
    δ: the spiral winds tighter and the computation needs more points."
14. `nb.explainer("vortex_sheet_rollup", heading="What is a velocity jump made of?", why="Setting the number of filaments
    and watching u(y) sharpen into a jump whose size is the sheet's strength, then rippling the sheet and pressing play to
    watch it roll up, links the jump, the vortex row and the instability — which needs motion.", tries=["N = 10 → 1000:
    the profile becomes a step of height γ.", "Drag the circuit's height: its circulation is γ ds once it spans the
    sheet.", "Roll-up mode, ▶: the ripple winds into cat's eyes.", "Switch to the caption's convention: only the sign
    changes."])`
15. `nb.md` — **What would change if…** "…the sheet were a wall? A no-slip wall with a stream above it is a vortex sheet
    stuck to the wall — its strength the wall slip it cancels. Viscosity then spreads it into a boundary layer (Ch. 9)."

### A.9 End matter — S01, S02, summary
1. `nb.pointer` — **S01** "📝 **Exercises 5.1–5.20** in the book: practise on them. We have written out, in our own words,
   the results the text leaves to Exercises 5.4 (the line vortex's zero net viscous force, D03), 5.5 (the lock exchange,
   D07), 5.9 (the Green's function, D10), 5.10 (Kelvin in a rotating frame, D19), 5.12 (Burgers' vortex, D18) and 5.18 (the
   centre of vorticity, D21), and stated 5.6 (N22) and 5.11 (N23); 5.7–5.8 (Vazsonyi and Crocco) return with compressible
   flow in Ch. 15."
2. `nb.pointer` — **S02** "📚 **Literature**: Lighthill (1986) for the knife-blade pair, Sommerfeld (1964) for Helmholtz's
   first theorem and leap-frogging rings; for more, Saffman's *Vortex Dynamics*, and Pedlosky for the geophysical side
   (Ch. 13)."
3. `nb.summary(clicked=[…], feeds_forward=[…], left_out=[…])` — **clicked** (one per CORE): C01 "A vortex tube has one
   strength along its whole length, $\int_V\nabla\cdot\boldsymbol\omega\,dV=-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0$ (5.4): it can thin
   and spin faster but never end in the fluid." · C02 "Pressure falls toward every vortex's centre to turn its fluid:
   $p-p_o=\tfrac18\rho\omega^2r^2-\rho gz$ (5.6) makes a bowl, $p-p_\infty=-\rho\Gamma^2/8\pi^2r^2-\rho gz$ (5.7) a funnel; the irrotational
   vortex is sheared but feels no net viscous force." · C03 "Round a material loop in an ideal, barotropic,
   conservatively forced fluid seen from an inertial frame, $D\Gamma/Dt=0$ (5.8); each broken hypothesis is a source of
   vorticity." · C04 "Where isopycnals cross isobars, pressure pushes through the centre while weight hangs off it: the
   spin-up is $\frac1{\rho^2}\nabla\rho\times\nabla p$." · C05 "Vortex lines move with the fluid, and tubes keep their strength (Helmholtz)." ·
   C06 "$\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$ (5.13): vorticity is carried, stretched and
   tilted, and diffused — never made by pressure in a uniform fluid." · C07 "Vorticity decides velocity:
   $\mathbf u=\frac1{4\pi}\int\frac{\boldsymbol\omega\times(\mathbf x-\mathbf x')}{\lvert\mathbf x-\mathbf x'\rvert^3}d^3x'$ (5.16) — with +1/(4π) in (5.14)." · C08 "A thin
   filament's pieces push as $d\mathbf u=\frac{\Gamma\,dl}{4\pi}\mathbf e_\omega\times\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ (5.17); an
   infinite line gives Γ/2πd back." · C09 "On a rotating planet with density contrasts,
   $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u+\frac1{\rho^2}\nabla\rho\times\nabla p+\nu\nabla^2\boldsymbol\omega$ (5.30)." · C10
   "Stretching a vortex line spins it up, tilting turns it (5.32) — neither happens in 2-D; Burgers' vortex balances
   stretching against diffusion." · C11 "The absolute circulation is conserved, $D\Gamma_a/Dt=0$ (5.33), so a column
   keeps (ζ + f)/h: squashed columns turn anticyclonic, stretched ones cyclonic." · C12 "Vortices are moved only by each
   other: like ones orbit their centre of vorticity at (Γ₁ + Γ₂)/2πh², opposite ones march at Γ/2πh." · C13 "A wall is a
   mirror vortex of opposite sign: a vortex slides along it at Γ/4πh; a ring widens and slows near it." · C14 "A velocity
   jump is a vortex sheet whose strength is the jump, γ = u₂ − u₁; left alone, it rolls up." **feeds_forward:** "Ch. 6:
   irrotational flow stays irrotational (Kelvin); images for cylinders." · "Ch. 8: vortex decay and the rotating cylinder
   (8.11)." · "Ch. 10: vorticity–stream-function methods and vortex methods (Biot–Savart)." · "Ch. 11: Kelvin–Helmholtz
   roll-up of a vortex sheet; baroclinic instability." · "Ch. 12: vortex stretching and the cascade; Burgers' vortex." ·
   "Ch. 13: (5.30), absolute circulation, (ζ + f)/h → potential vorticity, Rossby waves, fronts." · "Ch. 14: lifting line,
   trailing vortices, induced drag, ground effect." **left_out:** "Vazsonyi and Crocco equations (Exercises 5.7–5.8) →
   Ch. 15" · "vortex reconnection (named only) → Ch. 12" · "the viscous rotating-cylinder solution → Ch. 8 (8.11)".

### A.10 Placement table (ID → notebook section → host block → call)
Every curation ID appears once. `note` = `nb.note`; `md` = a markdown cell inside the block; `recap` = `nb.recap` before
the block's `nb.core`; `deriv` = `nb.derivation`.
- §5.1: R01 recap (before C01) · C01 core · N01 note, N02 note, N03 note, N45 note+plotly (in C01) · D01 deriv (C01) ·
  E1 explainer (C01) · R02, R03, R04, R05, R06, R07, R08 recap (before C02) · C02 core · N04, N05, N06, N07, N08, N09 note,
  N46/N47 note+plotly (in C02) · D02, D03 deriv (C02) · E2 explainer (C02).
- §5.2: R09 recap (before C03) · C03 core · N10, N11, N12, N13, N14, N16 note, N48 note+animation (in C03) · D04, D05 deriv ·
  E3 explainer (C03) · C04 core · N15 note, N49 note+figure, N50 note+figure (in C04) · D06, D07 deriv · E4 explainer (C04).
- §5.3: C05 core · N17 note, N51 note (in C05) · D08 deriv · pointer to E3.
- §5.4: R10, R11 recap (before C06) · C06 core · N18, N19, N20, N22, N23 note (in C06) · D09 deriv · pointer N21 → C10 ·
  C05's field-equation check (md + code + figure at the end of C06).
- §5.5: C07 core · N24, N25, N26, N27, N28, N52 note (in C07) · D10, D11 deriv · C08 core · N29 note · D12, D13 deriv ·
  E6 explainer (C08).
- §5.6: R12, R13, R14, R15, R16, R17 recap (before C09) · C09 core · N30, N31, N32, N33, N34, N35, N36 note, N37 note+figure
  (in C09) · D14, D15 deriv · C10 core · N38, N21 note, N53 note+plotly (in C10) · D16, D17, D18 deriv · E5 explainer (C10) ·
  R18 recap (before C11) · C11 core · N39 note, N54 note+figure (in C11) · D19, D20 deriv · E7 explainer (C11).
- §5.7: C12 core · N40, N41 note, N55/N56 note+animation (in C12) · D21 deriv · C13 core · N42, N43, N57, N58, N59 note (in
  C13) · D22 deriv · E8 explainer (C13).
- §5.8: C14 core · N44 note, N60 note+plotly (in C14) · D23 deriv · E9 explainer (C14).
- End: S01, S02 pointer · summary. B1 `vortex_rings` not embedded (backup).

## Part B — explainer storyboards

Common to all ten: created with `tools/new_viz.py`; `<meta name="viz:chapter" content="ch05">`; tabs Walkthrough ·
Explore · Explain · Derivation · Equations · Code · Check; `autoplay: false` and `play: false` on every walkthrough step
that quotes numbers (ch02 lesson); every displayed number is computed by a JS function that mirrors a `ch05` callable and
is proved by `selftest()` parity rows (`py:` expressions use only `ch05.…`, `np.pi`, lists, dicts, floats and strings —
no builtins, no lambdas; results indexed down to one number; functions that take fields are called with a **preset
name**, Part C convention 1); Explain is "Explanation & interpretation" in numbered sections built with `Viz.work.step /
line / box / table / hint / interpret`, modelled on `forced_damped_vibrations.html` and
`amplitude_phase_second_order_II_3.html` (0 what the views show and what each colour means · 1…n every displayed
quantity from the controls, "formula = substituted = result — why", results boxed · a section per view hidden on phones
or a hint · the values at the current time (live) · **Reading the current setting**); derivation steps are copied from
Part F (same `did` titles, same step count; phones shorten *why* to its first sentence; plain-text *why* never contains
raw TeX — write e^(αt), ∇ρ×∇p/ρ² as plain Unicode); every tour, Explain, Derivation, notes, status, equation and quiz
text that names a book equation **writes it out** next to its number (convention 9; `tools/eq_refs.py` → 0; the
`viz:concept` meta is worded without bare numbers). Equations in the walkthrough, Explain and quiz drafts below are written in Unicode for
readability; the builders set each one in TeX (`$…$`, backslashes doubled in JS strings) exactly as drafted. **Walkthrough texts ≤ 30 words** (ch04 lesson: cards must fit a
360×640 phone without more than two pager slices). Colours as in convention 8. A view hidden on portrait phones never
carries a step's key number (repeat it in a visible title or readout). Mode-dependent controls are hidden per mode with
a class toggled by the chapter script (the ch04 `.xxx-off` pattern), never merely marked `optional`. 3-D views use
`Viz.three` with the offline fallback note; their 2-D labels come from `T.project`. Parallel builders use private
scratch subfolders. Each explainer fits 360×640 … 1920×1080 and the 1000×700 notebook frame with no scrolling (fit plan
per explainer).

### E1 · vortex_tubes_cannot_end
- **Title:** "Can a vortex just stop in the middle of the water?" · **Summary:** "Slide a cross-section along a narrowing,
  twisting vortex tube: its area and mean spin trade off while the flux through it — the tube's strength — stays pinned;
  break the field and Gauss' budget fails." · **CORE:** C01 (also R01, N02, N03, N45) · **Reference:** `np_resnet_3d.html`
  (orbitable 3-D scene, a token travelling through, inspector synced to the selection) with
  `forced_damped_vibrations.html`'s explanation panel.
- **meta:** `viz:order 1` · `viz:sections 5.1` · `viz:equations 5.3 5.4` · `viz:fluidpy ch05.gaussian_tube_section
  ch05.tube_flux_budget ch05.gaussian_tube_vorticity ch05.vortex_ring_vorticity ch05.burgers_vortex` ·
  `viz:derivations D01`.
- **Physics (JS ↔ Python):** `tubeOmega(R, z, p)` ↔ `ch05.gaussian_tube_vorticity(R, z, Gamma, a0, L, twist)` (ω_R, ω_φ,
  ω_z) · `tubeSection(z, p)` ↔ `ch05.gaussian_tube_section(z, R0, Gamma, a0, L)` (radius R₀e^{−z/2L}, area, flux
  Γ(1 − e^{−R₀²/a₀²}), mean, peak) · `tubeBudget(za, zb, p)` ↔ `ch05.tube_flux_budget("gaussian_tube" | "broken", za, zb,
  R0, …)` (closed forms: lower −flux, side 0, upper +flux; broken: upper = flux·e^{−(zb−za)/L}) · `traceLine(R0, φ0, p)`
  (vortex line by RK4 on dx/ds = ω/\|ω\| — mirrors `ch05.vortex_line`; for the Gaussian tube the exact line is R = R₀e^{−z/2L},
  φ = φ₀ + twist·(…)·z, used as a check) · `ringFlux(p)` ↔ `ch05.vortex_ring_vorticity` integrated over a meridional
  section.
- **Modes** (chips "narrow · twisted · Burgers · ring · broken"): **narrow** (Gaussian tube, twist 0) · **twisted**
  (twist = 5 m⁻¹, helical lines) · **Burgers** (a straight tube: `burgers_vortex` ω_z, area constant — the flux and area
  both flat; also shows the *stream tube* of the same flow as a faint ghost, N45) · **ring** (a closed vortex ring: the
  section slider runs round the ring; every meridional section has flux Γ) · **broken** (ω = (Γ/πa₀²)e^{−R²/a₀²}e^{−z/L}e_z,
  not solenoidal — flagged).
- **Views** (rows [1.25, 1]):
  1. `scene` "The tube" (3-D via `Viz.three`, orbitable; z ∈ [0, 1.5] m, R ≤ 0.2 m, drawn ×2 radially for visibility and
     said so in the title) — 16 vortex lines (teal) seeded on the circle R₀ at z = 0, the tube wall as a translucent
     surface, the **lower lid** (blue disc at z = 0), a moving **section disc** (orange) at the slider position z_s, the
     field's ω arrows on a few lines; in broken mode the lines are straight but fade in colour upward (ω decreasing).
     Pointer: drag vertically on the tube to move the section; click a line → inspector.
  2. `profile` "Along the tube" — against z ∈ [0, 1.5] m: flux (teal, flat), area (blue, falling), mean ω = flux/area
     (purple, rising), log-y axis; the current z_s as a vertical marker with its three dots; ghost: the broken field's
     "flux" (rose dashed) in broken mode.
  3. `budget` "Gauss on the piece 0 → z_s" (`hidePortrait: true`) — three term bars lower (blue), side (grey), upper
     (orange) and the total (black) = ∫∇·ω dV; the key number (total) also in the `profile` title.
  Portrait: `scene` + `profile`; the budget's total sits in the profile title ("lower −0.632 · upper +0.632 · sum 0").
- **Controls (≤ 5 visible):** `zs` "Section height $z_s$" 0…1.5 m, step 0.01, default 0.5, help "where the orange lid
  cuts the tube" · `R0` "Seed radius $R_0$" 0.02…0.2 m, step 0.005, default 0.1, help "which tube (the lines through this
  circle at z = 0)" · `L` "Narrowing length $L$" 0.3…5 m, log, default 1, help "area falls by e every L metres" (narrow,
  twisted, broken) · `twist` "Twist" 0…10 m⁻¹, default 5 (twisted only) · `Gamma` "Vortex strength $\Gamma$" 0.1…5 m²/s,
  default 1 (optional). Transport `zs` 0 → 1.5 m, rate 0.25 m/s, `end: 'hold'`, hold 1.5 s; end-of-run `Viz.card`:
  "The section shrank to 22 % of its area; the flux stayed 0.632 m²/s; the mean spin rose ×4.5."
- **Depth features:** Explain + Code + Derivation · **3-D view** · **linked views** (3) · **transport** (scrub the section
  along the tube) · **terms** (lower/side/upper/total, click to isolate) · **presets** "narrowing tube" {mode: 'narrow',
  zs: 0.5, L: 1}, "strong squeeze" {L: 0.4, zs: 1.5}, "twisted" {mode: 'twisted'}, "straight (Burgers)" {mode: 'burgers'},
  "vortex ring" {mode: 'ring'}, "not a vorticity field" {mode: 'broken'} · **status** (solenoidal modes: "✅ same strength
  everywhere · Γ = 0.632 m²/s at z = 0.50 m"; broken: "⚠️ ∇·ω ≠ 0 — not a curl: the 'strength' leaks (sum −0.40 m²/s)";
  ring: "⭕ a closed tube: no ends at all") · **inspector** (click a line: "line through (R₀, φ₀): at z_s it is at
  R = 0.1·e^(−z/2L) = 0.0779 m; ω there = (ω_R, ω_φ, ω_z) = (…) s⁻¹; ω·n on the wall = 0").
- **Readouts:** "Flux $\Gamma$" (m²/s) · "Area $A$" (m²) · "Mean $\bar\omega$" (1/s) · "Gauss sum" (m²/s).
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**The tube** (3-D; drag to orbit) is made of <b class="c-teal">vortex lines</b> through a
     circle of radius R₀ at the bottom. The <b class="c-blue">blue</b> lid is where we start; the <b
     class="c-orange">orange</b> lid is the section you move. **Along the tube** plots, against height, the flux through
     the section (teal), its area (blue) and the mean vorticity (purple). **Gauss on the piece** (hidden on phones) shows
     the three surfaces of the piece between the lids."
  1. *The field* — $\omega_z=\frac{\Gamma}{\pi a^2}e^{-R^2/a^2}$ with $a^2=a_0^2e^{-z/L}$; with your numbers a(z_s) =
     0.1 × e^{−z_s/2L} = `live a` m, peak ω_z = Γ/πa² = `live peak` s⁻¹. "The core narrows upward; the radial part ω_R < 0 bends
     the lines inward."
  2. *The section* — "radius R₀e^{−z_s/2L} = `live r` m, area πR² = `live A` m²" (boxed).
  3. *The flux through it* — $\Gamma=\int_A\boldsymbol\omega\cdot\mathbf n\,dA=\Gamma_{tot}(1-e^{-R_0^2/a_0^2})$ = `live flux` m²/s —
     "the same number at every z: R²/a² is the same on the tube wall at every height, and so is the fraction of the
     vortex inside it." Boxed.
  4. *The spin it forces* — $\bar\omega=\Gamma/A$ = `live flux` / `live A` = **`live mean`** s⁻¹; "a quarter of the area means
     four times the spin."
  5. *Gauss' budget* — $\int_V\nabla\cdot\boldsymbol\omega\,dV=\int_A\boldsymbol\omega\cdot\mathbf n\,dA=-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0$
     (5.4): lower `live lo` + side `live side` + upper `live up` = **`live sum`** m²/s; "the side is zero because the wall
     is made of vortex lines, not because ω vanishes there."
  6. *The budget view* (hint on phones: "turn the phone to see the three surfaces as bars").
  7. *Vortex lines* — (5.3) $dx/\omega_x=dy/\omega_y=dz/\omega_z$: "each teal line is tangent to ω; we trace it with dx/ds = ω/\|ω\|."
  8. *Reading the current setting* — narrow/twisted: "The tube has thinned to `live pct` % of its starting area, so its
     fluid spins `live gain` times faster: one strength, carried from lid to lid (Helmholtz's second and third
     theorems). It can only end on a wall or close on itself." · Burgers: "A straight tube: area and flux both constant.
     The faint stream tube of the same flow does narrow — stream tubes and vortex tubes are different objects." · ring:
     "A closed tube has no ends: every cut through it carries Γ." · broken: "This field has ∇·ω = −ω_z/L ≠ 0: it is not the
     curl of any velocity, so it cannot be vorticity — its 'strength' falls from `live lo` to `live up` m²/s. Real
     vorticity can never do that."
- **Derivation tab:** **D01** (7 steps) `view: 'scene'`; goal page `set {mode: 'narrow', zs: 0.5}`; step 1 `set {zs: 0.5}`
  highlights both lids and the wall; step 5 `set {mode: 'twisted'}` with `watch` "the grey wall bar stays at zero even
  when the lines twist"; step 6 **live** "−(`lo`) + `up` = `sum` m²/s for your tube"; step 7 `set {zs: 1.5}` with `watch`
  "mean ω rises as the area falls; the flux bar does not move"; interpret `s => "Your tube carries ${flux} m²/s at every
  height; at z = ${zs} m it spins ${gain}× faster than at the bottom."`.
- **Code:**
  ```python
  sec = ch05.gaussian_tube_section({{zs}}, R0={{R0}}, L={{L}})   # the orange lid
  print(sec["flux"], sec["area"], sec["mean_omega"])      # {{flux}} m²/s, {{A}} m², {{mean}} 1/s
  b = ch05.tube_flux_budget("{{field}}", 0.0, {{zs}}, {{R0}}, L={{L}})
  print(b.lower, b.side, b.upper, b.total)                # {{lo}} {{side}} {{up}} → {{sum}}
  # (5.4): −Γ_lower + Γ_upper = 0 for any real vorticity field
  ```
- **Walkthrough (6 steps):** 1. "A tube of spin" — "A tornado, a smoke ring, a drain: tubes of spinning fluid. Can one
  just stop in mid-water?" `set {mode: 'narrow', zs: 0}`, highlight `view:scene` · 2. "Lines and a tube" — "Teal vortex
  lines follow ω, dx/ω_x = dy/ω_y = dz/ω_z (5.3). The lines through one circle form a tube." `eq: 'lines'` · 3. "Measure
  its strength" — "The flux of ω through a section is the tube's strength Γ — by Stokes, also the circulation round it."
  `readouts: ['flux']`, highlight `readout:flux` · 4. "Slide up the tube" — "Press ▶: the section climbs. Area falls, mean
  ω rises, flux stays put." `play: true`, `controls: ['zs']` · 5. "Why: Gauss" — "Lids and wall enclose a piece: the wall
  has no flux, so −Γ_lower + Γ_upper = 0 (5.4)." `terms: true`, `derive: {id: 'D01', step: 6}`, `code: {lines: [3, 4]}` ·
  6. "Your turn" — "Predict: halve L. Does the flux change? Then pick the broken field and see what must fail." `set {L:
  0.5}`, `controls: ['L']`.
- **Equations:** `lines` "Vortex lines" ref 'Eq. (5.3)' `dx/\omega_x=dy/\omega_y=dz/\omega_z` live "at your line: ω = (…) s⁻¹" ·
  `strength` "Tube strength" (no ref) `\Gamma=\oint_C\mathbf u\cdot d\mathbf x=\int_A\boldsymbol\omega\cdot\mathbf n\,dA` live "= 0.632 m²/s" ·
  `gauss` "Gauss on a piece of tube" ref 'Eq. (5.4)'
  `\begin{aligned}\textstyle\int_V\nabla\cdot\boldsymbol\omega\,dV&=\textstyle\int_A\boldsymbol\omega\cdot\mathbf n\,dA\\&=-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0\end{aligned}`
  live "−0.632 + 0 + 0.632 = 0" · `field` "The narrowing tube" `\omega_z=\frac{\Gamma}{\pi a^2}e^{-R^2/a^2},\ a^2=a_0^2e^{-z/L}` ·
  symbols ω, Γ, A, n, a, L with units.
- **Check yourself:** (1) "Halve L. What happens to the flux at z = 1 m, and to the mean vorticity?" — "The flux is
  unchanged (5.4); the area at z = 1 m falls to e⁻² of its start, so the mean vorticity rises ×7.4 instead of ×2.7."
  `set {L: 0.5, zs: 1}` · (2) "Why is the grey wall bar zero in the twisted tube, although ω is large on the wall?" — "The
  wall is made of vortex lines: ω lies along it, ω·n = 0." `set {mode: 'twisted'}` · (3) "In the broken field, where
  does the missing strength go?" — "Nowhere real: the field has sources (∇·ω ≠ 0), which a curl can never have — it is
  not vorticity." `set {mode: 'broken'}` · (4) "Make R₀ ten times smaller. Does the strength of the new tube stay
  constant along it?" — "Yes, a smaller value (Γ(1 − e^{−R₀²/a₀²}) = 1 − e^{−0.04} = 0.039 m²/s), still the same at every height." `set
  {R0: 0.02}`.
- **Selftest parity rows:** `{name: 'flux z=0', js: tubeSection(0, P0).flux, py:
  'ch05.gaussian_tube_section(0.0)["flux"]', rtol: 1e-12}` · `{name: 'mean z=1', js: tubeSection(1, P0).mean, py:
  'ch05.gaussian_tube_section(1.0)["mean_omega"]', rtol: 1e-12}` · `{name: 'area z=1 L=0.5', js: tubeSection(1, {...P0,
  L: 0.5}).area, py: 'ch05.gaussian_tube_section(1.0, L=0.5)["area"]', rtol: 1e-12}` · `{name: 'budget total', js:
  tubeBudget(0, 1, P0).total, py: 'ch05.tube_flux_budget("gaussian_tube", 0.0, 1.0, 0.1).total', rtol: 0, atol: 1e-8}` ·
  `{name: 'broken upper', js: tubeBudget(0, 1, {...P0, broken: true}).upper, py: 'ch05.tube_flux_budget("broken", 0.0,
  1.0, 0.1).upper', rtol: 1e-8}` · `{name: 'omega_R', js: tubeOmega(0.05, 0.3, P0)[0], py:
  'ch05.gaussian_tube_vorticity(0.05, 0.3)[0]', rtol: 1e-12}` · invariant `{name: 'line stays on wall', js:
  traceWallError(P0), expect: 0, atol: 1e-6}`.
- **Fit plan:** 360×640: `scene` (3-D, 55 %) + `profile`; status one line ("✅ Γ 0.632 m²/s · z 0.50 m"); mode chips wrap
  to two rows; presets hidden on short portrait screens (scoped CSS, Derivation tab exempt). Landscape phone: `scene` +
  `profile` in one row. Notebook 1000×700 and desktop: rows [1.25, 1], `budget` beside `profile`. If three.js fails to
  load, `scene` falls back to a 2-D (R, z) cut with the lids as lines (said in the view title).

### E2 · vortex_pressure_funnel
- **Title:** "Why is a spinning bucket a bowl and a drain a funnel?" · **Summary:** "Four vortices on one pair of axes —
  a spinning tank, an ideal line vortex, a Rankine core, a rotating cylinder — with their free surfaces and isobars, the
  radial force balance at a probe, and the viscous stress that the irrotational vortex has but does not feel." ·
  **CORE:** C02 (also R02, R03, R05, R06, R07, R08, N05, N06, N07, N08, N09) · **Reference:**
  `angular_frequency_explorer_1.html` (modes on one stage, "Right now" notes, a table of real vortices with the current
  row highlighted); links (not duplicated) to ch03 `vortex_paddle_wheels` and ch04 `which_bernoulli`.
- **meta:** `viz:order 2` · `viz:sections 5.1` · `viz:equations 5.1 5.2 5.5a 5.5b 5.6 5.7` · `viz:fluidpy
  ch05.vortex_pressure_scenario ch05.solid_body_pressure ch05.line_vortex_pressure ch05.rankine_pressure
  ch05.isobar_height ch05.rotating_tank_free_surface ch05.line_vortex_viscous_stress ch05.torque_per_length` ·
  `viz:derivations D02 D03`.
- **Physics (JS ↔ Python):** `scenario(kind, r, z, p)` ↔ `ch05.vortex_pressure_scenario(kind, r, z, rho, mu, g, **p)`
  returning u_θ, p, B, σ_rθ, net force, centripetal, dp_dr, surface_z (closed forms: (5.6), (5.7), Rankine matched at a,
  cylinder = Rankine with Γ = πa²ω; σ_rθ = μ r d(u_θ/r)/dr) · `tankSurface(R, depth, Om)` ↔
  `ch05.rotating_tank_free_surface(R, depth, Omega_tank)` · `isobar(r, kind, param, dp)` ↔ `ch05.isobar_height` ·
  `stressForce(kind, r, p)` ↔ `ch05.vortex_stress_force(kind, r, mu, **p)`.
- **Modes** (chips "tank · line · Rankine · cylinder"): **tank** (solid body, ω = 2Ω_tank, the tank's walls drawn, the
  water volume conserved) · **line** (ideal line vortex, funnel without a bottom — clipped at r = 2 mm with a note) ·
  **Rankine** (core a; the tornado) · **cylinder** (a solid cylinder of radius a spinning at ω/2; fluid only outside).
- **Views** (rows [1.3, 1]):
  1. `section` "The vortex in cross-section" (`equal: false`; r ∈ [−r_max, r_max], z below and above the rest level;
     r_max = 0.15 m for tank/line, 150 m for the tornado preset — axis units switch with a note) — the free surface (teal
     thick), three isobars (orange dashed) labelled with their gauge pressure, a colour wash of p (orange ramp), the axis,
     the tank walls or the cylinder (grey), and a **probe element** (a small sector at radius r_p, drawn deformed by the
     local shear in line/Rankine-outside, undeformed in the tank). Pointer: drag the probe along r; click anywhere →
     inspector (p arithmetic at that (r, z)).
  2. `profiles` "u_θ, B and σ_rθ against r" — u_θ(r) (teal), B(r) − B_ref (amber), σ_rθ(r) (rose, right axis), the probe
     marked; ghosts: the other two vortices' u_θ faint.
  3. `balance` "Forces on the probe" (`hidePortrait: true`) — term bars per unit mass: centripetal u_θ²/r (teal) vs
     pressure gradient (1/ρ)∂p/∂r (orange) — equal (5.5a) — and the net viscous force (rose: 0 for tank and line, nonzero
     only across the Rankine core edge, where the vorticity jumps).
  Portrait: `section` + `profiles`; the balance's numbers go into the `profiles` title ("u_θ²/r = (1/ρ)∂p/∂r = 25 m/s²").
- **Controls (≤ 5 visible, per mode):** `Om` "Tank rate $\Omega_{tank}$" 0…20 rad/s, step 0.1, default 5 (tank, cylinder),
  help "the vorticity is ω = 2Ω_tank" · `Gamma` "Circulation $\Gamma$" 0.01…2 m²/s, log, default 1 (line, Rankine) · `a`
  "Core radius $a$" 0.005…0.1 m, default 0.05 (Rankine, cylinder) · `rp` "Probe radius $r$" 0.005…r_max, default 0.1 ·
  `fluid` select water/air (optional; ρ, μ from `ch01.fluid_properties` at 20 °C, hard-coded: water 998.2/1.002e-3, air
  1.204/1.82e-5). Tank geometry fixed: R = 0.1 m, depth 0.1 m, height 0.2 m.
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (4) · **presets** "bucket 5 rad/s"
  {mode: 'tank', Om: 5}, "fast bucket (bottom shows)" {Om: 20}, "bathtub drain" {mode: 'line', Gamma: 0.01}, "strong
  drain" {Gamma: 1}, "tornado (air)" {mode: 'rankine', fluid: 'air', Gamma: 1e4, a: 50}, "spinning rod" {mode:
  'cylinder', Om: 10, a: 0.01} · **status** (tank: "🥣 rotational: B grows outward, no viscous stress"; line: "🌪
  irrotational: σ_rθ = −0.0318 Pa, net viscous force 0"; Rankine: "🌀 core r < a rotates, outside irrotational"; tank with
  z_v < 0: "⚠️ the bottom is uncovered out to r = … m") · **terms** (the radial balance) · **inspector** (click (r, z):
  "p − p_o = ⅛ρω²r² − ρgz = ⅛ × 1000 × 10² × 0.08² − 1000 × 9.81 × (−0.03) = 80 + 294.3 = **374.3** Pa") · **notes**
  (regime text below).
- **Readouts:** "u_θ at probe" (m/s) · "Pressure p − p_ref" (Pa) · "Surface dip" (mm) · "Stress σ_rθ" (Pa) · "Net viscous
  force" (N/m³).
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**The vortex in cross-section**: the <b class="c-teal">free surface</b>, <b
     class="c-orange">isobars</b> (dashed) and the pressure colour; the small sector is a probe element. **u_θ, B, σ_rθ**
     plots the swirl (teal), the Bernoulli function (amber) and the viscous shear stress (rose) against radius. **Forces
     on the probe** (hidden on phones) weighs the inward push against the pressure gradient."
  1. *The swirl at the probe* — tank: $u_\theta=\omega r/2$ (5.1) with ω = 2Ω_tank = `live w` s⁻¹: `live u` m/s; line:
     $u_\theta=\Gamma/2\pi r$ (5.2) = `live u` m/s; Rankine: inside Γr/2πa², outside Γ/2πr. Boxed.
  2. *The inward push it needs* — u_θ²/r = `live c` m/s²; "every parcel on a circle is accelerating toward the axis."
  3. *The pressure gradient that supplies it* — (5.5a) $-\rho u_\theta^2/r=-\partial p/\partial r$ → ∂p/∂r = ρu_θ²/r = `live dpdr`
     Pa/m; vertically (5.5b) $0=-\partial p/\partial z-\rho g$.
  4. *The pressure field* — tank (5.6) $p-p_o=\tfrac18\rho\omega^2r^2-\rho gz$ = `live p` Pa; line (5.7)
     $p-p_\infty=-\frac{\rho\Gamma^2}{8\pi^2r^2}-\rho gz$ = `live p` Pa; Rankine: continuous at a, centre deficit
     ρΓ²/4π²a² = `live pc` Pa. Boxed.
  5. *The free surface* — p = p_atm: tank $z=\frac{\omega^2r^2}{8g}$ + vertex; rim–centre height ω²R²/8g = `live h` mm;
     volume kept: vertex `live zv` m; line: dip Γ²/(8π²r²g) = `live dip` m at the probe — "the book calls these
     'hyperboloids of the second degree'; the formula (c − z)r² = const is a cubic."
  6. *Bernoulli across the streamlines* — tank: B − B(0) = ω²r²/4 = `live B` m²/s² (grows: rotational, R07); line: B
     constant (irrotational) — "compare Ch. 4's which_bernoulli."
  7. *Stress and net force* — σ_rθ = μ r d(u_θ/r)/dr = `live s` Pa; net force (1/r²)d(r²σ_rθ)/dr = `live F` N/m³ —
     "irrotational means no net viscous force, not no stress" (N09); torque per length 2πr²σ_rθ = −2μΓ = `live T` N m/m
     at every r (N08).
  8. *Reading the current setting* — tank: "The surface is a bowl because the needed push grows with r; no element is
     deformed, so there is no viscous stress at all — only while the tank spins up does viscosity act." · line: "The
     push grows toward the axis like 1/r³, so the surface plunges: a funnel with no bottom. Elements are sheared (σ_rθ ≠
     0) but the stresses on each element's faces cancel." · Rankine: "Inside the core the fluid turns rigidly (a bowl),
     outside it swirls irrotationally (a funnel): together a funnel with a rounded bottom — a tornado's pressure dip." ·
     cylinder: "The spinning rod drives the irrotational outer flow; the torque 2μΓ it applies travels undiminished to
     any radius — a container would need to resist it."
- **Derivation tab:** **D02** (10 steps) `view: 'section'`; goal `set {mode: 'tank', Om: 5, rp: 0.1}`; step 4 draws
  e_θ turning at the probe (`watch` "the unit vector turns toward the axis"); step 5 **live** "at r = 0.1 m: ρu_θ²/r =
  1000 × 0.5²/0.1 = 2500 Pa/m"; steps 8–10 draw the bowl growing from its parts (r-part, then z-part); step 10 **live**
  "rim–centre height = 12.7 mm"; interpret with the reader's Ω. **D03** (9 steps) `view: 'section'`; goal `set {mode:
  'line', Gamma: 1, rp: 0.1}`; step 4 **live** "σ_rθ = −μΓ/πr² = −0.0318 Pa"; step 6 zooms the probe element showing the
  shear stress arrows on its inner and outer faces (inner larger, outer on a larger face and lever arm); step 7 `watch`
  "the rose force bar stays at 0 while σ_rθ is not"; interpret "your line vortex: σ_rθ = ${s} Pa, net force 0."
- **Code:**
  ```python
  sc = ch05.vortex_pressure_scenario("{{mode}}", {{rp}}, 0.0, **{{params}})
  print(sc["u_theta"], sc["centripetal"], sc["dp_dr"])   # {{u}} m/s, {{c}} m/s², {{dpdr}} Pa/m
  print(sc["p"], sc["surface_z"])                         # {{p}} Pa, {{zs}} m
  print(sc["sigma_rtheta"], sc["net_viscous_force"])     # {{s}} Pa, {{F}} N/m³
  tank = ch05.rotating_tank_free_surface(0.1, 0.1, {{Om}})
  print(tank["z_rim"] - tank["z_vertex"])                  # (5.6): ω²R²/8g = {{h}} m
  ```
- **Walkthrough (7 steps):** 1. "Bowl or funnel?" — "A spun bucket's water forms a bowl; a drain makes a narrow funnel.
  Same physics?" `set {mode: 'tank', Om: 5}` · 2. "Turning needs a push" — "A parcel on a circle accelerates inward at
  u_θ²/r. Only pressure can push it: −ρu_θ²/r = −∂p/∂r (5.5a)." `terms: true`, `eq: 'radial'` · 3. "The bowl" — "u_θ = ωr/2:
  pressure grows as r², p − p_o = ⅛ρω²r² − ρgz (5.6). Rim–centre height 12.7 mm." `readouts: ['dip']`, `derive: {id:
  'D02', step: 10}` · 4. "The funnel" — "u_θ = Γ/2πr: the push soars near the axis, p − p_∞ = −ρΓ²/8π²r² − ρgz (5.7)."
  `set {mode: 'line', Gamma: 1}`, `eq: 'funnel'` · 5. "Stress without force" — "The line vortex is sheared, σ_rθ = −μΓ/πr²,
  yet the net viscous force is zero." `readouts: ['stress', 'force']`, `derive: {id: 'D03', step: 7}`, `code: {lines: [4,
  4]}` · 6. "A real tornado" — "A Rankine core gives the funnel a bottom: 600 Pa at the core edge, twice that at the
  centre." `set {mode: 'rankine', fluid: 'air', Gamma: 1e4, a: 50}` · 7. "Your turn" — "Predict: double the tank rate —
  how much deeper is the bowl? Then drag Ω." `set {mode: 'tank', Om: 5}`, `controls: ['Om']`.
- **Equations:** `sb` "Solid-body rotation" ref 'Eq. (5.1)' `u_\theta=\omega r/2` · `lv` "Line vortex" ref 'Eq. (5.2)'
  `u_\theta=\Gamma/2\pi r` · `radial` "Radial balance" ref 'Eq. (5.5a)' `-\rho u_\theta^2/r=-\partial p/\partial r` live "2500 = 2500 Pa/m" ·
  `hydro` "Vertical balance" ref 'Eq. (5.5b)' `0=-\partial p/\partial z-\rho g` · `bowl` "Tank pressure" ref 'Eq. (5.6)'
  `p-p_o=\tfrac18\rho\omega^2r^2-\rho gz` live · `funnel` "Line-vortex pressure" ref 'Eq. (5.7)'
  `p-p_\infty=-\frac{\rho\Gamma^2}{8\pi^2r^2}-\rho gz` live · `stress` "Shear stress of the line vortex" (no ref)
  `\sigma_{r\theta}=-\mu\Gamma/\pi r^2` live.
- **Check yourself:** (1) "Double the tank rate. How much deeper is the bowl?" — "Four times: the height ω²R²/8g grows with
  the square of the rate (12.7 → 51.0 mm)." `set {mode: 'tank', Om: 10}` · (2) "In which vortex is there no viscous stress
  at all, and why?" — "The tank: rigid rotation does not deform elements, so S = 0 and σ = 0." · (3) "The line vortex has
  σ_rθ ≠ 0. Why does its velocity profile never change?" — "The net viscous force is zero: the stress on each element's
  outer face balances the inner one (r²σ_rθ is constant)." `set {mode: 'line'}` · (4) "Tornado preset: where is the
  pressure lowest, and by how much?" — "On the axis: ρΓ²/4π²a² = 1216 Pa, twice the core-edge deficit." `set {mode:
  'rankine', fluid: 'air', Gamma: 1e4, a: 50}`.
- **Selftest parity rows:** `{name: 'tank p', js: scenario('solid', 0.1, 0, {omega: 10}).p, py:
  'ch05.solid_body_pressure(0.1, 0.0, 10.0)', rtol: 1e-12}` · `{name: 'rim height', js: tankSurface(0.1, 0.1, 5).zrim -
  tankSurface(0.1, 0.1, 5).zv, py: 'ch05.rotating_tank_free_surface(0.1, 0.1, 5.0)["rim_rise"] +
  ch05.rotating_tank_free_surface(0.1, 0.1, 5.0)["centre_drop"]', rtol: 1e-10}` · `{name: 'funnel p', js:
  scenario('line', 0.1, 0, {Gamma: 1}).p, py: 'ch05.line_vortex_pressure(0.1, 0.0, 1.0)', rtol: 1e-12}` · `{name: 'funnel
  isobar', js: isobar(0.1, 'line', 1, 0), py: 'ch05.isobar_height(0.1, 1.0, kind="line")', rtol: 1e-12}` · `{name:
  'rankine centre', js: scenario('rankine', 0, 0, {Gamma: 1, a: 0.1}).p, py: 'ch05.rankine_pressure(0.0, 0.0, 1.0, 0.1)',
  rtol: 1e-12}` · `{name: 'sigma', js: scenario('line', 0.1, 0, {Gamma: 1}).sigma, py: 'ch05.line_vortex_viscous_stress(0.1,
  1.0, 0.001)', rtol: 1e-12}` · `{name: 'torque', js: torque(0.3, 1, 1e-3), py: 'ch05.torque_per_length(0.3, 1.0, 0.001)',
  rtol: 1e-12}` · invariant `{name: 'line net force', js: scenario('line', 0.07, 0, {Gamma: 1}).force, expect: 0, atol:
  1e-12}`.
- **Fit plan:** 360×640: `section` + `profiles`; status one line ("🥣 bowl · 12.7 mm"); mode chips short ("tank · line ·
  core · rod"); presets on a second strip hidden on short portrait screens. Landscape phone: one row `section` +
  `profiles`. Desktop/notebook: rows [1.3, 1], `balance` beside `profiles`. The tornado preset rescales the axes and the
  titles say "m" instead of "cm".

### E3 · kelvin_material_loop
- **Title:** "Stretch and tangle a loop of dye — what stays the same?" · **Summary:** "A loop of dye carried by five flows
  on one clock: its circulation stays flat while it is tangled — until you break one of Kelvin's four hypotheses and
  the matching term of (5.10) switches on; a Helmholtz mode keeps a loop on a tube wall at zero flux." · **CORE:** C03,
  C05 (also N10, N11, N12, N13, N14, N16, N17, R09, N48) · **Reference:** `forced_damped_vibrations.html` (one time slider
  for the phenomenon and the graph, a ghost reference, regime-dependent explanation); link to ch02
  `stokes_circulation_loop` (one instant) — E3 adds time and material loops.
- **meta:** `viz:order 3` · `viz:sections 5.2 5.3` · `viz:equations 5.8 5.9 5.10 5.11` · `viz:fluidpy
  ch05.kelvin_scenario ch05.kelvin_scenario_circulation ch05.kelvin_scenario_rate ch05.lamb_oseen_circulation
  ch05.kelvin_hypotheses_text ch05.material_circulation` · `viz:derivations D04 D05 D08`.
- **Physics (JS ↔ Python):** `flowU(name, x, y, t, p)` ↔ the `u` of `ch05.kelvin_scenario(name)` (cellular, Rankine,
  Lamb–Oseen with t₀ = 10 s, lock-exchange fluid at rest, rotating convergence) · `advectLoop(pts, t0, t1, p)` (RK4, 256
  points, dt = 0.01 s; point insertion when a gap exceeds 0.05 m, max 2048) ↔ `ch05.material_loop` · `loopGamma(pts, name,
  t)` (periodic trapezoid with FFT-free centred differences of the parametrisation) ↔ `ch05.loop_circulation` ·
  `kelvinTerms(name, t, p)` ↔ `ch05.kelvin_scenario_rate(name, t)` (acceleration, contour, pressure ∮dp/ρ, body, viscous,
  coriolis) · `lambOseen(r, t, p)` ↔ `ch05.lamb_oseen_circulation` · `hypothesesText(inv, baro, cons, inert)` ↔
  `ch05.kelvin_hypotheses_text` · `helmholtzLoop(t)` (the circle painted on the wall of the stretched Gaussian tube: it is carried inward and
  upward and stays on the narrowing wall) ↔ `ch05.kelvin_scenario("helmholtz")` with `ch05.material_circulation`.
- **Modes** (`loop` chips "material · fixed"; `flow` chips "cellular · Rankine · viscous · baroclinic · rotating ·
  Helmholtz"): the six flows below; **fixed loop** keeps the loop's points still while the flow passes (Γ(t) changes in
  inviscid flow — the "material" matters).
- **Views** (rows [1.3, 0.9]):
  1. `flow` "The loop in the flow" (`equal: true`; cellular x, y ∈ [0, π] m; Rankine [−2, 2] m; viscous [−1, 1] cm;
     baroclinic [−0.5, 0.5] × [0, 1] m with blue isopycnal shading and orange isobars; rotating [−1.5, 1.5] m with a
     turning frame arrow; Helmholtz: a side view (R, z) of the stretched tube's wall with the painted loop and the vertical vortex lines) — grey streamlines, the loop (teal,
     thickening to purple where most stretched), one element dx (purple arrow) with u and u + du (Fig. 5.4 analogue),
     tracer dots. Pointer: drag the loop's centre (at t = 0 only); click a loop point → inspector.
  2. `gamma` "Γ(t)" — Γ(t) so far (bold teal) over its faint full run, the ghost Γ(0) (grey dashed); in the viscous mode the
     Lamb–Oseen closed form (rose dashed); in rotating mode Γ_a (amber, flat); in fixed-loop mode the material Γ as ghost.
     The loop length (purple, right axis, log) climbs.
  3. `table` "Which hypothesis holds?" (`hidePortrait: true`) — the four-restriction table (inviscid · barotropic ·
     conservative · inertial) with ✓/✗ for the current flow and the surviving term of (5.10) as a coloured bar: ∮dp/ρ
     orange, ∮dΦ grey, viscous rose, Coriolis amber, contour term ∮u·du black (≈ 0 always).
  Portrait: `flow` + `gamma`; the verdict and the surviving term go into the status badge.
- **Controls (≤ 5 visible):** `t` time (transport; per flow t_end: cellular 6 turnover times ≈ 37.7 s, Rankine 20 s,
  viscous 20 s, baroclinic 5 s, rotating 5 s, Helmholtz 2 s) · `loopr` "Loop radius" (per flow range; cellular 0.2…0.8 m,
  default 0.5) · `nu` "Viscosity ν" 1e-7…1e-5 m²/s, log, default 1e-6 (viscous only) · `Om` "Frame rotation Ω" 0…1 rad/s,
  default 0.5 (rotating only) · `loop` chips material/fixed. Transport rate per flow (≈ t_end/8 per s), `end: 'hold'`,
  end-of-run cards: cellular "Loop 6.7× longer; Γ = 1.155 m²/s unchanged (to 1e-9)" · viscous "Γ fell from 0.00465 to
  0.00268 m²/s: vorticity diffused out through the loop" · rotating "Γ rose to 1.99 m²/s; Γ_a stayed π — see (5.33)".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **modes** (six flows ×
  material/fixed) · **presets** "ideal cellular flow" {flow: 'cellular', loop: 'material'}, "loop across a Rankine core"
  {flow: 'rankine'}, "viscous vortex" {flow: 'viscous'}, "lock exchange (baroclinic)" {flow: 'baroclinic'}, "rotating
  frame" {flow: 'rotating'}, "fixed loop (not material)" {flow: 'cellular', loop: 'fixed'}, "Helmholtz: loop on a tube
  wall" {flow: 'helmholtz'} · **status** ("✅ Kelvin holds: DΓ/Dt = 0" / "🍯 viscous: DΓ/Dt = ∮(1/ρ)∂σ_ij/∂x_j dx_i ≠ 0" /
  "🌡 baroclinic: ∮dp/ρ ≠ 0" / "🔄 rotating frame: Γ changes, Γ_a = Γ + 2∫Ω·n dA does not" / "📌 fixed loop: not material —
  Kelvin does not apply" / "🧵 Helmholtz: the loop stays on the tube wall, flux ≈ 0") — the badge text is
  `kelvin_hypotheses_text` with the equation written out · **terms** (the (5.10) bars) · **inspector** (click a loop
  point: "u = (…), dx/ds = (…), u·dx/ds = …; this point's share of Γ").
- **Readouts:** "Γ now" (m²/s) · "Loop length" (m) · "DΓ/Dt" (m²/s²) · "Γ_a" (rotating only).
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**The loop in the flow**: grey streamlines, the <b class="c-teal">loop of dye</b> (purple
     where most stretched), one element dx with u and u + du. **Γ(t)** plots the circulation so far against its start
     (grey dashed). **Which hypothesis holds?** (hidden on phones) ticks Kelvin's four restrictions and shows the one term
     of (5.10) that survives."
  1. *The loop now* — length `live L` m (×`live stretch` since t = 0); area `live A` m².
  2. *Its circulation* — Γ = ∮u·dx by the periodic trapezoid rule over N = `live N` points = **`live G`** m²/s (boxed);
     "flux route ∫ω·n dA agrees at t = 0: `live G0`."
  3. *The rate, term by term* — (5.9) $\frac{D\Gamma}{Dt}=\oint_C\frac{Du_i}{Dt}dx_i+\oint_Cu_i\frac{D}{Dt}(dx_i)$: acceleration term `live acc`,
     contour term ∮u·du = `live con` (≈ 0: exact differential). (5.10)
     $\oint_C\frac{Du_i}{Dt}dx_i=-\oint_C\frac1\rho dp-\oint_Cd\Phi+\oint_C\big(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j}\big)dx_i$: pressure
     `live pr`, body `live bo`, viscous `live vi`, Coriolis `live co` — sum **`live rate`** m²/s².
  4. *Which hypothesis is broken* — the table as text: "inviscid ✓ · barotropic ✓ · conservative ✓ · inertial ✓ → (5.8)
     $D\Gamma/Dt=0$" or the broken one with its term.
  5. *The viscous case* — Lamb–Oseen $\Gamma(r,t)=\Gamma_0(1-e^{-r^2/4\nu t})$ with r = `live r` m, t + t₀ = `live tt` s:
     `live GLO` m²/s; (5.11) $\frac{D\Gamma}{Dt}=\oint_C(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j})dx_i$ = `live vi` m²/s².
  6. *The table view* (hint on phones).
  7. *At the current time* — t = `live t` s, Γ = `live G`, Γ − Γ(0) = `live dG`.
  8. *Reading the current setting* — cellular: "The loop has been stretched ×`live stretch` into a spiral, yet the fluid
     on it is never accelerated along it on balance: pressure and gravity are exact differentials, so Γ is frozen —
     Kelvin." · Rankine: "Part of the loop is inside the rotating core, part outside; shear winds it up, Γ = 1.098 m²/s
     stays." · viscous: "Vorticity diffuses outward through the loop: the net viscous force along it is negative and Γ
     leaks away." · baroclinic: "Density changes sideways while pressure changes vertically: ∮dp/ρ ≠ 0, circulation is
     being made — the lock exchange starts to tumble (C04)." · rotating: "In the turning frame the Coriolis force does
     net work round the loop as its area shrinks: Γ grows, but Γ + 2Ω·A_vec does not (C11)." · fixed: "Fluid of different
     vorticity streams through a loop fixed in space; Kelvin is about loops that move with the fluid." · Helmholtz: "The
     small loop was painted on a vortex tube's wall: zero circulation then, zero now — so it is still on a tube wall.
     Vortex lines move with the fluid."
- **Derivation tab:** **D04** (8 steps) `view: 'flow'`; goal `set {flow: 'cellular', t: 0}`; step 1 numbers the particles
  (tags drawn); step 4 `set {t: 1.5}` `watch` "the acceleration curve (teal) and the contour curve (grey)"; step 5 draws
  the element with u and u + du; step 8 **live** "contour term = `con` m²/s² for your loop". **D05** (9 steps) `view:
  'gamma'`; goal `set {flow: 'cellular'}`; each substitution step lights one term bar; step 6 `set {flow: 'baroclinic'}`
  `watch` "the orange ∮dp/ρ bar is not zero here"; step 7 `set {flow: 'viscous'}` **live** "viscous integral = dΓ/dt =
  −3.35e-4 m²/s²"; step 8 `set {flow: 'cellular'}`. **D08** (9 steps) `view: 'flow'`; goal `set {flow: 'helmholtz', t:
  0}`; step 4 `set {t: 2}` `watch` "the flux through the small loop stays ≈ 0"; interpret "Γ of the painted loop: ${G}
  m²/s, compared with |ω|A = ${wA}: it is still on the wall."
- **Code:**
  ```python
  sc = ch05.kelvin_scenario("{{flow}}")          # velocity, loop, hypotheses
  t = np.linspace(0, {{t}}, 50)
  G = ch05.material_circulation(sc["u"], sc["pts0"], t)
  print(G[0], G[-1])                              # {{G0}} → {{G}} m²/s
  r = ch05.kelvin_scenario_rate("{{flow}}", {{t}})
  print(r["pressure"], r["viscous"], r["coriolis"]) # the (5.10) terms
  print(ch05.kelvin_hypotheses_text({{hyp}}))     # {{verdict}}
  ```
- **Walkthrough (7 steps):** 1. "A loop of dye" — "Put a ring of dye in a swirling flow. Press ▶ and watch it be
  stretched into a spiral." `set {flow: 'cellular', t: 0}`, `play: true` · 2. "Γ does not move" — "Its circulation
  Γ = ∮u·dx stays 1.155 m²/s while the loop grows ×6.7." `set {t: 37.7}`, `play: false`, `readouts: ['G', 'len']` ·
  3. "Why: only acceleration counts" — "DΓ/Dt = ∮(Du/Dt)·dx + ∮u·du (5.9); the second is ∮d(½u²) = 0." `derive: {id:
  'D04', step: 8}`, `eq: 'k59'` · 4. "Forces round a loop" — "Pressure and gravity give −∮dp/ρ − ∮dΦ: zero when ρ = ρ(p)
  and Φ exists. So DΓ/Dt = 0 (5.8)." `terms: true`, `derive: {id: 'D05', step: 6}` · 5. "Break it: viscosity" —
  "Lamb–Oseen vortex: Γ leaks as vorticity diffuses out, DΓ/Dt = ∮(1/ρ)∂σ_ij/∂x_j dx_i (5.11)." `set {flow: 'viscous'}`,
  `code: {lines: [5, 6]}` · 6. "Break it: density" — "Heavy beside light fluid: ∮dp/ρ ≠ 0 and circulation appears."
  `set {flow: 'baroclinic'}`, `inspect: true` · 7. "Your turn" — "Predict first: in a loop fixed in space, is Γ conserved?
  Switch to 'fixed' and check." `set {flow: 'cellular', loop: 'fixed'}`.
- **Equations:** `k58` "Kelvin" ref 'Eq. (5.8)' `D\Gamma/Dt=0` · `k59` "Rate of a material loop's circulation" ref
  'Eq. (5.9)' `\frac{D\Gamma}{Dt}=\oint_C\frac{Du_i}{Dt}dx_i+\oint_Cu_i\frac{D}{Dt}(dx_i)` live · `k510` "Forces round the loop" ref
  'Eq. (5.10)'
  `\begin{aligned}\oint_C\tfrac{Du_i}{Dt}dx_i&=-\oint_C\tfrac1\rho dp-\oint_Cd\Phi\\&\quad+\oint_C\big(\tfrac1\rho\tfrac{\partial\sigma_{ij}}{\partial x_j}\big)dx_i\end{aligned}` live ·
  `k511` "Viscous leak" ref 'Eq. (5.11)' `\frac{D\Gamma}{Dt}=\oint_C\Big(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j}\Big)dx_i` live · `lo`
  "Lamb–Oseen circulation" (no ref) `\Gamma(r,t)=\Gamma_0(1-e^{-r^2/4\nu t})` live.
- **Check yourself:** (1) "Double the loop radius in the cellular flow. Is Γ still constant?" — "Yes, at a new value (the
  loop encloses more vorticity); Kelvin holds for every material loop." `set {flow: 'cellular', loopr: 0.8}` · (2) "Why
  does Γ change in the fixed-loop mode although the flow is ideal?" — "Different fluid, carrying different vorticity,
  passes through a loop fixed in space; Kelvin is only about material loops." `set {loop: 'fixed'}` · (3) "In the
  baroclinic flow nothing moves at t = 0. Why does Γ start to change?" — "∮dp/ρ ≠ 0 because ρ is not a function of p:
  the pressure torque of C04." `set {flow: 'baroclinic', t: 0}` · (4) "Rotating frame: which quantity stays flat?" —
  "Γ_a = Γ + 2Ω·A_vec: the loop's shrinking area is paid for by relative circulation." `set {flow: 'rotating'}`.
- **Selftest parity rows:** `{name: 'cellular G0', js: loopGamma(initLoop('cellular', P), 'cellular', 0), py:
  'ch05.kelvin_scenario_circulation("cellular", 0.0)', rtol: 1e-8}` · `{name: 'cellular G(6s)', js: gammaAt('cellular',
  6), py: 'ch05.kelvin_scenario_circulation("cellular", 6.0)', rtol: 1e-6}` · `{name: 'rankine G0', js: gammaAt('rankine',
  0), py: 'ch05.kelvin_scenario_circulation("rankine_straddle", 0.0)', rtol: 1e-4} (the velocity has a kink where the loop crosses the core edge, so a 256-point JS loop sum is good to 1e-5, not spectrally)` · `{name: 'LO Gamma', js:
  lambOseen(0.005, 10, {G0: 0.01, nu: 1e-6})[0], py: 'ch05.lamb_oseen_circulation(0.005, 10.0, 0.01, 1e-6)[0]', rtol:
  1e-12}` · `{name: 'rotating G(5)', js: gammaAt('rotating', 5), py: 'ch05.kelvin_scenario_circulation("rotating", 5.0)',
  rtol: 1e-8}` · `{name: 'baroclinic rate', js: kelvinTerms('baroclinic', 0).pressure, py:
  'ch05.kelvin_scenario_rate("baroclinic", 0.0)["pressure"]', rtol: 1e-6}` · `{name: 'verdict text', js:
  hypothesesText(true, false, true, true), py: 'ch05.kelvin_hypotheses_text(True, False, True, True)', exact: true}`.
- **Fit plan:** 360×640: `flow` + `gamma`; status one line ("✅ Kelvin · Γ 1.155"); flow chips short ("cells · core · visc
  · baro · rot · helm"), wrapping; presets hidden on short portrait screens. Landscape phone: one row. Desktop/notebook:
  rows [1.3, 0.9], `table` beside `gamma`.

### E4 · baroclinic_torque
- **Title:** "How can density make fluid spin?" · **Summary:** "A fluid disc with isobars and isopycnals: turn one set
  against the other and the pressure force's line of action leaves the centre of mass; the torque, divided by the
  moment of inertia, is the baroclinic spin-up ∇ρ×∇p/ρ² — then the same torque tips over a lock exchange." · **CORE:**
  C04 (also N14, N15, N49, N50; forward link N35 (5.28)) · **Reference:** `amplitude_phase_second_order_II_3.html`
  (several windows on one state, a numbered live derivation in the explanation).
- **meta:** `viz:order 4` · `viz:sections 5.2 5.6` · `viz:equations 5.28` · `viz:fluidpy
  ch05.pressure_torque_on_element ch05.baroclinic_rate_2d ch05.lock_exchange_initial_vorticity_rate
  ch05.lock_exchange_fields ch05.baroclinic_term` · `viz:derivations D06 D07`.
- **Physics (JS ↔ Python):** `element(p)` → {F, xG, torque, IG, spin, baro, ratio} with the closed forms of D06 (x_G =
  R²∇ρ/4ρ₀, M_G = πR⁴(∇ρ×∇p)/4ρ₀, I_G = πR⁴(ρ₀ − \|∇ρ\|²R²/8ρ₀)/2) ↔ `ch05.pressure_torque_on_element(radius=R,
  grad_rho=…, grad_p=…, rho0=…)` · `baroRate(gr, gp, rho)` ↔ `ch05.baroclinic_rate_2d` · `lockRate(r1, r2, d)` ↔
  `ch05.lock_exchange_initial_vorticity_rate` · `lockField(x, y, p)` (tanh step, hydrostatic p with ρ̄) ↔
  `ch05.lock_exchange_fields` + `ch05.baroclinic_term`.
- **Modes** (chips "element · lock exchange"): **element** (the disc) · **lock** (a tank with the interface; the initial
  rate map and the tumbling arrows; a short kinematic "t = 0⁺" animation of the interface tilting at the initial rate,
  labelled as a first-instant picture).
- **Views** (rows [1.3, 1]):
  1. `disc` "The fluid element" (`equal: true`, x, y ∈ [−1.6R, 1.6R]) — orange isobars (horizontal: ∇p = −ρ₀g e_y) and blue
     dashed isopycnals (tilted by the angle θ between ∇ρ and ∇p); inward pressure arrows round the rim (orange, length ∝ p
     − p_min); the geometric centre O, the centre of mass G (offset exaggerated with a "× 10⁴" tag, and the true value in
     the title), the net-force line through O (black), the lever arm O–G (purple), a curved torque arrow (teal for
     counterclockwise, rose for clockwise). Pointer: drag round the rim to rotate the isopycnals (sets θ); click a rim point
     → inspector.
  2. `curve` "Spin-up rate vs angle" — (∇ρ×∇p)_z/ρ² against θ ∈ [0, 360°] (orange sine), the current θ as a dot; the torque
     route (teal ◇ measured, from `element`) on the curve; ghost: θ = 0 (barotropic) line at zero.
  3. `lock` "Lock exchange" (`hidePortrait: true`) — the tank (x ∈ [−0.5, 0.5] m, y ∈ [0, 1] m), heavy ρ₂ left (dark blue),
     light ρ₁ right, the interface of thickness δ, the colour map of (∇ρ×∇p)_z/ρ² (orange) and counterclockwise tumbling
     arrows; in lock mode it takes row 0.
  Portrait: `disc` + `curve` (element mode), `lock` + `curve` (lock mode: curve shows the rate vs δ, 1/δ).
- **Controls (≤ 5 visible):** `theta` "Angle between ∇ρ and ∇p" 0…360°, step 1, default 90 · `grho` "Density gradient
  |∇ρ|" 0…50 kg/m⁴, default 10 · `R` "Element radius R" 0.001…0.5 m, log, default 0.01 · (lock) `drho` "ρ₂ − ρ₁" 0…50
  kg/m³, default 25 · `delta` "Interface thickness δ" 0.01…0.5 m, log, default 0.1. ρ₀ = 1000 kg/m³, g = 9.81 fixed.
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (2) · **presets** "barotropic (0°)"
  {theta: 0}, "stable layering (180°)" {theta: 180}, "sea breeze (air)" {theta: 90, rho0: 1.2, grho: 4e-7} (ours: a 5 K contrast over 50 km, spin-up
  3.3×10⁻⁶ s⁻², about 0.012 s⁻¹ after an hour; `rho0` is a preset-only parameter, 1000 otherwise), "fresh/salt lock exchange" {mode: 'lock', drho: 25, delta: 0.1}, "sharp front" {mode: 'lock', delta:
  0.01} · **status** ("⚖ barotropic: ∇ρ ∥ ∇p, no torque" / "🌀 baroclinic: spins counterclockwise at 0.0981 s⁻²" / "🌀
  clockwise at …"; lock: "🌊 lock exchange: 2.42 s⁻² at the interface") · **inspector** (click a rim point: "p = p₀ +
  ∇p·x = …; force −p n ds = (…); its moment about G = (x − x_G) × f = …") · **terms** (in the Explain: torque from the
  upper half vs lower half of the rim).
- **Readouts:** "Offset of G" (m) · "Torque about G" (N m/m) · "Spin-up (torque)" (1/s²) · "∇ρ×∇p/ρ²" (1/s²).
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**The fluid element**: <b class="c-orange">isobars</b> (solid), <b class="c-blue">isopycnals</b>
     (dashed), pressure pushes round the rim, the centre O, the centre of mass G, the net-force line and the torque arrow.
     **Spin-up vs angle** plots the formula (orange) with the torque route (◇). **Lock exchange** (hidden on phones in
     element mode) shows the same torque on a whole interface."
  1. *The two gradients* — ∇p = (0, −ρ₀g) = (0, −9810) Pa/m; ∇ρ = \|∇ρ\|(sin θ, −cos θ)… (angle measured from ∇p) =
     (`live grx`, `live gry`) kg/m⁴.
  2. *Their cross product* — (∇ρ×∇p)_z = `live gx`×`live py` − `live gy`×`live px` = **`live cross`** Pa kg m⁻⁵ (boxed); "zero
     when the isolines are parallel."
  3. *Where the pressure force acts* — F = −πR²∇p = `live F` N/m, through O (every rim push points at O).
  4. *Where the weight hangs* — x_G = R²∇ρ/4ρ₀ = `live xG` m toward the dense side.
  5. *The torque* — M_G = −x_G × F = πR⁴(∇ρ×∇p)/4ρ₀ = **`live M`** N m/m.
  6. *Spin-up* — I_G = ½ρ₀πR⁴ = `live I` kg m; Dω/Dt = 2M_G/I_G = **`live spin`** s⁻², vs the formula (5.28)
     $\frac1{\rho^2}\nabla\rho\times\nabla p$ = `live baro` s⁻²; ratio `live ratio` (→ 1 as R → 0).
  7. *The lock exchange* (lock mode; hint in element mode) — $\frac{D\omega_z}{Dt}=\frac{2(\rho_2-\rho_1)g}{(\rho_2+\rho_1)\delta}$ =
     2 × `live drho` × 9.81/(`live sum` × `live delta`) = **`live lock`** s⁻².
  8. *Reading the current setting* — θ = 0 or 180°: "Isolines parallel: the weight hangs straight below the line of the
     pressure force; no lever arm, no spin — however strong the density contrast (the ocean's stable layers)." · 0 < θ <
     180°: "The dense side is off the force line: the element turns counterclockwise at `spin` s⁻², the heavy side going
     down. This is how sea breezes, fronts and gravity currents start." · lock: "Heavy water on the left slumps under the
     light water: the interface tumbles counterclockwise at 2.42 s⁻² at first; a sharper interface spins faster — in the
     limit, a vortex sheet (C14)."
- **Derivation tab:** **D06** (10 steps) `view: 'disc'`; goal `set {mode: 'element', theta: 90, R: 0.01}`; step 1 lights
  the rim arrows; step 3 draws the arrows meeting at O; step 5 slides G (`watch` "G moves toward the dense side"); step
  6 draws the lever arm; step 7 **live** "π(0.01)⁴ × (−98 100)/(4 × 1000) = −7.70e-7 N m/m"; step 10 **live** "−0.0981
  s⁻² by the torque, −0.0981 by the formula". **D07** (7 steps) `view: 'lock'`; goal `set {mode: 'lock', drho: 25,
  delta: 0.1}`; step 3 draws ∇p with ρ̄; step 6 **live** "2 × 25 × 9.81/(2025 × 0.1) = 2.42 s⁻²"; step 7 `watch` "the
  curved arrows turn counterclockwise".
- **Code:**
  ```python
  r = ch05.pressure_torque_on_element(radius={{R}},
          grad_rho={{grho}}, grad_p=(0.0, -9810.0))
  print(r["x_G"], r["torque"], r["I_G"])     # {{xG}} m, {{M}} N m/m, {{I}}
  print(r["spin_up"], r["baroclinic"])        # {{spin}}, {{baro}} 1/s²  (5.28)
  rate = ch05.lock_exchange_initial_vorticity_rate(
          1000.0, {{rho2}}, {{delta}})         # {{lock}} 1/s²
  ```
- **Walkthrough (6 steps):** 1. "Sideways density" — "Warm land, cool sea: density changes sideways, gravity pulls down.
  Why does the air turn over?" `set {mode: 'element', theta: 0}` · 2. "Parallel isolines" — "Isobars and isopycnals
  parallel: pressure pushes through the centre, the weight hangs there too. No torque." `readouts: ['torque']` ·
  3. "Cross them" — "Tilt the isopycnals: the centre of mass G slides toward the dense side; the push still goes through
  O." `set {theta: 90}`, `derive: {id: 'D06', step: 5}` · 4. "Torque → spin" — "Torque ÷ inertia × 2 = ∇ρ×∇p/ρ², the
  baroclinic term (5.28): −0.098 s⁻² here." `eq: 'baro'`, `readouts: ['spin', 'baro']`, `code: {lines: [4, 4]}` · 5. "A
  lock exchange" — "Heavy left, light right, gate pulled: 2(ρ₂ − ρ₁)g/((ρ₂ + ρ₁)δ) = 2.42 s⁻²." `set {mode: 'lock'}`,
  `eq: 'lock'` · 6. "Your turn" — "Predict: halve δ. Faster or slower? Then drag the thickness." `controls: ['delta']`.
- **Equations:** `baro` "Baroclinic term" ref 'Eq. (5.28)'
  `-\varepsilon_{nqi}\Big(\frac1\rho p_{,i}\Big)_{,q}=\frac1{\rho^2}[\nabla\rho\times\nabla p]_n` live "= −0.0981 s⁻²" · `torque` "Torque on
  the disc" (no ref) `M_G=\frac{\pi R^4}{4\rho_0}\nabla\rho\times\nabla p` live · `inertia` "Moment of inertia" (no ref)
  `I_G=\tfrac12\rho_0\pi R^4` · `lock` "Lock exchange" (no ref) `\frac{D\omega_z}{Dt}=\frac{2(\rho_2-\rho_1)g}{(\rho_2+\rho_1)\delta}` live.
- **Check yourself:** (1) "Set θ = 180° (density increasing downward). Torque?" — "Zero: ∇ρ ∥ ∇p; stable layering does
  not spin." `set {theta: 180}` · (2) "Shrink R from 0.1 m to 1 mm. What happens to the two spin-up routes?" — "They
  converge: the torque route differs from ∇ρ×∇p/ρ² by a fraction \|∇ρ\|²R²/8ρ₀² — a point has exactly the formula."
  `set {R: 0.001}` · (3) "Why does the pressure force pass through O, however the pressure varies?" — "On a circle every
  push points at the centre (x ∥ n); with linear p the resultant has no moment about O." · (4) "Lock exchange: what
  happens as δ → 0?" — "The rate grows as 1/δ; the circulation per unit height stays finite: a vortex sheet forms." `set
  {mode: 'lock', delta: 0.01}`.
- **Selftest parity rows:** `{name: 'spin R=1cm', js: element({R: 0.01, grho: 10, theta: 90}).spin, py:
  'ch05.pressure_torque_on_element()["spin_up"]', rtol: 1e-9}` · `{name: 'torque', js: element({R: 0.01, grho: 10, theta:
  90}).torque, py: 'ch05.pressure_torque_on_element()["torque"]', rtol: 1e-6}` · `{name: 'xG', js: element({R: 0.01, grho:
  10, theta: 90}).xG[0], py: 'ch05.pressure_torque_on_element()["x_G"][0]', rtol: 1e-6}` · `{name: 'baro formula', js:
  baroRate([10, 0], [0, -9810], 1000), py: 'ch05.baroclinic_rate_2d([10.0, 0.0], [0.0, -9810.0], 1000.0)', rtol: 1e-12}` ·
  `{name: 'lock rate', js: lockRate(1000, 1025, 0.1), py: 'ch05.lock_exchange_initial_vorticity_rate(1000.0, 1025.0, 0.1)',
  rtol: 1e-12}` · invariant `{name: 'barotropic no torque', js: element({R: 0.05, grho: 20, theta: 0}).torque, expect: 0,
  atol: 1e-15}`.
- **Fit plan:** 360×640: `disc` + `curve`; status one line ("🌀 ccw · 0.098 s⁻²"); the × 10⁴ tag on G stays inside the
  view (`Viz.text` with `bg: true`). Landscape phone: one row. Desktop/notebook: rows [1.3, 1], `lock` beside `curve`.

### E5 · vorticity_stretching_tilting
- **Title:** "How can a flow spin fluid faster without a torque?" · **Summary:** "A vortex line and a material element
  ride a chosen flow on one clock: stretching (purple) spins the line up, tilting (blue) turns it, both vanish in 2-D;
  viscosity (rose) joins in and Burgers' vortex holds stretching against diffusion — with the terms of (5.13) as bars." ·
  **CORE:** C06, C10 (also N19, R11, N20, N22, N38, N21, N53; C05's field-equation check) · **Reference:**
  `forced_damped_vibrations.html` (transient → steady state on one clock, boxed numbers) with `fid_formula_lab.html`'s
  term bars; the term-bar pattern of ch04 `navier_stokes_term_balance` (linked).
- **meta:** `viz:order 5` · `viz:sections 5.4 5.6` · `viz:equations 5.12 5.13 5.31 5.32` · `viz:fluidpy
  ch05.stretching_tilting_split ch05.uniform_strain_vorticity ch05.strain_preset ch05.stretched_tube ch05.burgers_vortex
  ch05.burgers_core_radius ch05.vorticity_budget_preset ch05.helix_frame` · `viz:derivations D09 D16 D17 D18`.
- **Physics (JS ↔ Python):** `strain(name, rate)` ↔ `ch05.strain_preset` · `omegaAt(w0, name, t, rate)` (closed forms of
  e^{Gt}ω₀: diagonal exponentials; nilpotent shear I + Gt; custom via a 3 × 3 Taylor/Padé expm) ↔
  `ch05.uniform_strain_vorticity` · `split(w, G)` ↔ `ch05.stretching_tilting_split` · `elementAt(dx0, name, t)` (the same
  e^{Gt} — a material element; the angle to ω stays 0) ↔ `ch05.frozen_in_check` (ABC row for the invariant) · `burgers(R,
  p)` ↔ `ch05.burgers_vortex` · `burgersTerms(R, p)` (closed forms: advective (αR/2)(αR/2ν)ω, stretching αω, diffusion
  their difference) ↔ `ch05.vorticity_budget_preset("burgers", x=R)` · `helixFrame(s, a, c)` ↔ `ch05.helix_frame`.
- **Modes** (flow chips "axial · shear · planar · Burgers · custom"): **axial** (G = diag(−α/2, −α/2, α), ω along z or at the
  chosen angle) · **shear** (G[2, 0] = s: w = s x; ω₀ along x) · **planar** (G = diag(s, −s, 0), ω ⟂ plane) · **Burgers**
  (axisymmetric strain + ν: the steady Gaussian core; the 3-D view shows the stretched tube and the profile) · **custom**
  (a 3 × 3 G from three sliders: stretch, shear, rotation — optional controls).
- **Views** (rows [1.3, 1]):
  1. `scene` "The vortex line" (3-D, `Viz.three`, orbitable) — a vortex line segment (a helix of a = 1, c = 0.3 in axial
     mode, a straight segment otherwise) drawn as a tube whose radius shrinks as \|ω\| grows (volume kept) with a spinning
     stripe; at a marked point its frame e_s (teal), e_n (orange, away from the centre of curvature), e_m (blue); the
     stretching arrow (purple, along e_s) and tilting arrow (blue, across); a material element δx (grey) riding along —
     parallel to ω while ν = 0. Burgers mode: the tube with inflow arrows and the core radius √(4ν/α) marked.
  2. `growth` "\|ω\|(t) and its turning" — \|ω\|(t)/\|ω₀\| (purple, log-y) with the ghost e^{αt} (grey dashed) or 1 + s t,
     and the angle of ω from its start (blue, right axis); Burgers mode: ω_z(R) profile (teal) with the core radius mark.
  3. `bars` "Terms of (5.13)" (`hidePortrait: true`) — local (blue dashed), advective (grey), stretching (purple), tilting
     (blue), diffusion (rose), residual (black), per component at the marked point (Burgers: at R = the probe radius).
  Portrait: `scene` + `growth`; the stretching rate and the growth factor appear in the growth title ("rate e_s·Ge_s =
  1.00 s⁻¹ · \|ω\| × 7.39").
- **Controls (≤ 5 visible):** `rate` "Strain rate α (or s)" 0…2 s⁻¹, default 1 · `ang` "Initial tilt of ω" 0…90°, default 0
  (angle from the stretching axis) · `nu` "Viscosity ν" 1e-7…1e-4 m²/s, log, default 1e-6 (Burgers only) · `probe`
  "Probe radius R" 0…8 mm, default 1 (Burgers only) · `t` time (transport 0 → 3 s, rate 0.5, `end: 'hold'`; end card:
  "\|ω\| grew ×20.1 = e³ with no torque — pure stretching"). Custom mode adds `g1`, `g2`, `g3` (optional).
- **Depth features:** Explain + Code + Derivation · **3-D view** · **linked views** (3) · **transport** · **terms** ·
  **presets** "skater: axial stretch" {flow: 'axial', rate: 1, ang: 0}, "stretch at 45°" {ang: 45}, "tilting by shear"
  {flow: 'shear'}, "flat world (2-D)" {flow: 'planar'}, "Burgers' vortex" {flow: 'burgers', nu: 1e-6}, "compression"
  {flow: 'axial', rate: -0.5} (allowed via the preset only) · **status** ("⬆ stretching: \|ω\| × 7.39 at t = 2 s" / "↪
  tilting: ω_z = s t ω_x" / "▭ 2-D: neither — ω ⟂ plane" / "⚖ Burgers: stretching = −diffusion + advection, core 2.0
  mm") · **inspector** (click the marked point: "G e_s = (…); e_s·G e_s = … s⁻¹; stretching = … e_s; tilting = …").
- **Readouts:** "Stretch rate e_s·Ge_s" (1/s) · "Growth \|ω\|/\|ω₀\|" · "Tilt angle" (deg) · "Core √(4ν/α)" (mm, Burgers).
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**The vortex line** (3-D; drag to orbit) with its frame e_s, e_n, e_m, the <b
     class="c-accent">stretching</b> arrow along it and the <b class="c-blue">tilting</b> arrow across it; the grey arrow
     is a material element. **|ω|(t)** tracks the growth (purple) and the turning (blue). **Terms of (5.13)** (hidden on
     phones) balances local, advective, stretching, tilting and <b class="c-rose">diffusion</b>."
  1. *The flow and the vorticity now* — G = `live G`; ω = e^{Gt}ω₀ = (`live w`) s⁻¹ (Cauchy: a weak vortex line in a steady
     linear flow evolves like a material element).
  2. *The line's direction* — e_s = ω/\|ω\| = (`live es`).
  3. *Stretching* — rate e_s·G e_s = **`live rate`** s⁻¹; stretching term (e_s·Ge_s)ω = (`live st`) s⁻² — (5.32)
     $\frac{D\omega_s}{Dt}=\omega\frac{\partial u_s}{\partial s}$.
  4. *Tilting* — Gω − stretching = (`live ti`) s⁻², perpendicular to ω — $\frac{D\omega_n}{Dt}=\omega\frac{\partial u_n}{\partial s}$,
     $\frac{D\omega_m}{Dt}=\omega\frac{\partial u_m}{\partial s}$.
  5. *Growth since t = 0* — \|ω\|/\|ω₀\| = **`live gain`** (axial: e^{αt} = e^{`live at`}); angle turned `live turn`°.
  6. *Diffusion and Burgers* (Burgers mode; hint otherwise) — core √(4ν/α) = √(4 × `live nu`/`live al`) = **`live core`** mm;
     at R = `live R` mm: advective `live adv`, stretching `live str`, diffusion `live dif`, residual ≈ 0 — (5.13)
     $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$.
  7. *At the current time* — t = `live t` s.
  8. *Reading the current setting* — axial: "The line is pulled longer along its own direction: thinner, it must spin
     faster (angular momentum, no torque). Turbulence builds its thin intense vortices this way." · at an angle: "Part of
     the stretching turns ω toward the stretching axis — stretching and tilting at once." · shear: "The flow varies along
     the line in a direction across it: the line is turned and a new component appears, ω_z = s t ω_x." · planar: "In a
     flat flow ω points out of the plane and nothing varies along it: no stretching, no tilting — why 2-D turbulence
     cannot concentrate vorticity." · Burgers: "Stretching piles vorticity into the core, diffusion spreads it out; they
     balance at a core radius √(4ν/α) = `core` mm, whatever the vortex started as."
- **Derivation tab:** **D09 ★★★** (12 steps) `view: 'bars'` (phones keep `bars` — the Derivation tab keeps its view);
  goal `set {flow: 'burgers', probe: 1}`; steps 2–3 grey out a "pressure" and a "gravity" placeholder bar (crossed:
  "curl of a gradient = 0"); step 5 shows ω × u as an arrow at the probe; step 8 **live** "Burgers at R = 1 mm: local 0 +
  advective 15.49 = stretching 61.97 + diffusion −46.48"; step 12 `set {flow: 'axial'}` lights the three bars of (5.13).
  **D16** (6 steps) `view: 'scene'`; goal `set {flow: 'shear'}`; step 1 draws e_s, e_n, e_m; step 6 `watch` "purple arrow
  along ω, blue across". **D17** (5 steps) `view: 'scene'`; step 1 `set {flow: 'axial', nu: 0}`; step 3 **live** "e^(αt) =
  7.39 at t = 2 s"; step 5 `set {flow: 'planar'}`. **D18** (9 steps) `view: 'growth'`; step 3 **live** "ω₀ = 10 s⁻¹, L 1 →
  2 m: ω = 20 s⁻¹"; step 4 `set {flow: 'burgers', rate: 1}`; step 9 **live** "core 2.0 mm, peak 79.6 s⁻¹ for Γ = 10⁻³ m²/s".
- **Code:**
  ```python
  G = ch05.strain_preset("{{flow}}", rate={{rate}})
  w = ch05.uniform_strain_vorticity({{w0}}, G, {{t}})   # e^{Gt} ω0 = {{w}}
  s = ch05.stretching_tilting_split(w, G)               # (5.31)–(5.32)
  print(s["rate"], s["stretching"], s["tilting"])      # {{rate_s}}, {{st}}, {{ti}}
  print(ch05.burgers_core_radius({{rate}}, {{nu}}))     # √(4ν/α) = {{core}} m
  ```
- **Walkthrough (7 steps):** 1. "Spin without torque" — "A skater spins faster by pulling in her arms. Can a flow do that
  to vorticity?" `set {flow: 'axial', t: 0}` · 2. "The vorticity equation" — "Curl Navier–Stokes: pressure and gravity
  vanish, leaving Dω/Dt = (ω·∇)u + ν∇²ω (5.13)." `eq: 'e513'`, `derive: {id: 'D09', step: 12}` · 3. "Stretch it" — "Press
  ▶: the flow pulls the line longer; \|ω\| grows as e^{αt}." `play: true`, `readouts: ['gain']` · 4. "Along vs across" —
  "(ω·∇)u = \|ω\|∂u/∂s (5.31): the along part stretches (purple), the across part tilts (blue)." `set {ang: 45, t: 1}`,
  `play: false`, `inspect: true` · 5. "Tilting" — "Shear w = s x turns ω_x into ω_z = s t ω_x: vorticity in a new
  direction." `set {flow: 'shear', t: 2}`, `code: {lines: [3, 4]}` · 6. "Flat world" — "In 2-D, ω ⟂ plane: neither
  stretching nor tilting. Nothing happens." `set {flow: 'planar'}` · 7. "Your turn" — "Burgers: predict the core radius
  for ν = 4×10⁻⁶ m²/s, then drag ν." `set {flow: 'burgers'}`, `controls: ['nu']`.
- **Equations:** `e512` "Curl of Navier–Stokes" ref 'Eq. (5.12)'
  `\nabla\times\Big\{\frac{D\mathbf u}{Dt}=-\frac1\rho\nabla p+\mathbf g+\nu\nabla^2\mathbf u\Big\}` · `e513` "Vorticity equation" ref 'Eq. (5.13)'
  `\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega` live · `e531` "Natural coordinates" ref
  'Eq. (5.31)' `(\boldsymbol\omega\cdot\nabla)\mathbf u=\omega\,\partial\mathbf u/\partial s` · `e532` "Stretching and tilting" ref 'Eq. (5.32)'
  `\begin{aligned}\tfrac{D\omega_s}{Dt}&=\omega\tfrac{\partial u_s}{\partial s}\\\tfrac{D\omega_n}{Dt}&=\omega\tfrac{\partial u_n}{\partial s},\ \tfrac{D\omega_m}{Dt}=\omega\tfrac{\partial u_m}{\partial s}\end{aligned}`
  live · `burg` "Burgers' vortex" (no ref) `\omega_z=\frac{\alpha\Gamma}{4\pi\nu}e^{-\alpha R^2/4\nu}` live.
- **Check yourself:** (1) "Axial strain α = 1 s⁻¹: how long until \|ω\| doubles?" — "ln 2/α = 0.69 s; stretching
  multiplies ω by e^{αt}." `set {flow: 'axial', rate: 1}` · (2) "Put ω at 90° to the stretching axis (along x). Does it
  grow?" — "No — it shrinks as e^{−αt/2}: the flow compresses along x (G_xx = −α/2)." `set {ang: 90}` · (3) "Planar flow:
  why do the bars all vanish?" — "ω is along z and nothing varies with z: (ω·∇)u = ω ∂u/∂z = 0." `set {flow: 'planar'}` ·
  (4) "Burgers: quadruple ν. The core?" — "Doubles: √(4ν/α)." `set {flow: 'burgers', nu: 4e-6}`.
- **Selftest parity rows:** `{name: 'axial growth', js: omegaAt([0, 0, 1], 'axial_stretch', 2, 1)[2], py:
  'ch05.uniform_strain_vorticity([0, 0, 1], "axial_stretch", 2.0)[2]', rtol: 1e-12}` · `{name: 'shear tilt', js: omegaAt([1,
  0, 0], 'shear_tilt', 2, 1)[2], py: 'ch05.uniform_strain_vorticity([1, 0, 0], "shear_tilt", 2.0)[2]', rtol: 1e-12}` ·
  `{name: 'split rate', js: split([1, 0, 1], strain('axial_stretch', 1)).rate, py:
  'ch05.stretching_tilting_split([1.0, 0.0, 1.0], ch05.strain_preset("axial_stretch"))["rate"]', rtol: 1e-12}` ·
  `{name: 'burgers core', js: Math.sqrt(4e-6/1), py: 'ch05.burgers_core_radius(1.0, 1e-6)', rtol: 1e-12}` · `{name:
  'burgers peak', js: burgers(0, {G: 1e-3, a: 1, nu: 1e-6}).wz, py: 'ch05.burgers_vortex(0.0, 0.0, 0.001, 1.0,
  0.000001)[3]', rtol: 1e-12}` · `{name: 'burgers stretching R=1mm', js: burgersTerms(1e-3, P).str, py:
  'ch05.vorticity_budget_preset("burgers", x=0.001)["stretching_tilting"]', rtol: 1e-5}` · `{name: 'helix curvature', js:
  helixFrame(0, 1, 0.3).kappa, py: 'ch05.helix_frame(0.0)["curvature"]', rtol: 1e-12}` · invariant `{name: 'tilting ⟂ ω',
  js: dot(split([1, 0, 1], strain('axial_stretch', 1)).tilting, [1, 0, 1]), expect: 0, atol: 1e-14}`.
- **Fit plan:** 360×640: `scene` + `growth`; status one line ("⬆ ×7.39 · 1.00 s⁻¹"); the D09 derivation keeps `bars` on
  phones (twelve steps; the four longest lines split per convention). Landscape phone: one row. Desktop/notebook: rows
  [1.3, 1], `bars` beside `growth`. three.js fallback: a 2-D projection with the frame arrows.

### E6 · biot_savart_filament
- **Title:** "How does a vortex push water far away?" · **Summary:** "Drag a point around a straight, bent, square, ring or
  helical filament: every segment's push (5.17) is drawn and summed; a long straight segment gives back Γ/2πd, and the
  printed sign of (5.14) would spin everything backwards." · **CORE:** C07, C08 (also N24, N25, N26, N27, N28, N29, N52)
  · **Reference:** `forward_noising_lab.html` (click a point to see its exact arithmetic, a sum shown as parts + total)
  with `forced_damped_vibrations.html`'s explanation panel.
- **meta:** `viz:order 6` · `viz:sections 5.5` · `viz:equations 5.14 5.15 5.16 5.17` · `viz:fluidpy
  ch05.segment_induced_velocity ch05.segment_speed ch05.filament_velocity ch05.filament_contributions
  ch05.filament_preset ch05.filament_velocity_preset ch05.ring_axis_velocity` · `viz:derivations D10 D11 D12 D13`.
- **Physics (JS ↔ Python):** `segVel(x, a, b, G, eps)` (closed form (Γ/4πd)(cos θ_a − cos θ_b) along e_ω × e_d; the
  standard vector form Γ/4π · (r₁ × r₂)/\|r₁ × r₂\|² · r₀·(r₁/\|r₁\| − r₂/\|r₂\|)) ↔ `ch05.segment_induced_velocity` ·
  `polyline(name, M, p)` ↔ `ch05.filament_preset` · `filVel(x, poly, G, closed)` and `contribs(…)` ↔
  `ch05.filament_velocity`, `ch05.filament_contributions` · `ringAxis(z, R, G)` ↔ `ch05.ring_axis_velocity` ·
  `tubeU(r, G, sigma, sign)` = sign·Γ(1 − e^{−r²/σ²})/2πr ↔ the Gaussian tube test of `ch05.velocity_from_curl_omega`
  (parity with the closed form only; the 3-D quadrature stays in Python).
- **Modes** (`shape` chips "segment · half-line · square · ring · helix · tube"): **segment** (straight, half-length ℓ) ·
  **half-line** (semi-infinite, 1000 m) · **square** (side 1 m) · **ring** (radius R) · **helix** (R = 0.3, pitch 0.4, 3
  turns) · **tube** (the volume form (5.16) for a Gaussian core: the line becomes a tube of radius σ, and the sign toggle
  shows (5.14) printed vs corrected).
- **Views** (rows [1.3, 1]):
  1. `scene` "Filament and field point" (3-D `Viz.three`, orbitable; the segment along z ∈ [−ℓ, ℓ], the ring in z = 0) —
     the filament as M segments (teal, arrowheads along e_ω), the field point x (draggable in the plane y = 0 by the
     pointer on the 2-D overlay; its (x, z) shown), each segment's contribution du as a thin arrow at x (purple, opacity ∝
     \|du\|), their sum (black thick), for the segment mode the angles θ_a, θ_b drawn at the ends and the distance d.
     Pointer: drag x; click a segment → inspector.
  2. `unrolled` "Contribution along the filament" — \|du\|/dl (purple curve) against the arc length along the filament, the
     shaded area = the speed at x (the "unrolled integral"), the element clicked marked; tube mode: u_θ(r) printed vs
     corrected (teal vs rose) against r.
  3. `distance` "Speed vs distance" (`hidePortrait: true`) — segment mode: speed at distance d (log-log) for the current ℓ
     (teal) with the infinite-line ghost Γ/2πd (grey dashed) and the semi-infinite half (dotted); ring mode: on-axis
     u_z(z) with the closed form ΓR²/2(R² + z²)^{3/2}; the current point marked.
  Portrait: `scene` + `unrolled`; the total speed and the closed form go in the `unrolled` title ("sum 0.1125 = closed
  form 0.1125 m/s").
- **Controls (≤ 5 visible):** `Gamma` "Circulation Γ" 0.1…5 m²/s, default 1 · `ell` "Half-length ℓ (or ring R)" 0.1…20 m,
  log, default 1 · `M` "Segments M" 4…200, log-int, default 32 · `px`, `pz` field point (via drag; readouts show them;
  sliders optional) · `sign` toggle "Use the printed −1/(4π)" (tube mode; default off).
- **Depth features:** Explain + Code + Derivation · **3-D view** · **linked views** (3) · **inspector** (click a segment:
  "dl = 0.0625 m at x′ = (…); x − x′ = (…), \|x − x′\| = …; e_ω × (x − x′)/\|x − x′\|³ = (…); du = Γdl/4π × that = (…) m/s") ·
  **presets** "short segment" {shape: 'segment', ell: 1}, "long segment → line" {ell: 20}, "half-line" {shape:
  'half-line'}, "square loop" {shape: 'square'}, "ring, point at centre" {shape: 'ring', ell: 0.5, px: 0, pz: 0}, "helix"
  {shape: 'helix'}, "tube: printed sign" {shape: 'tube', sign: true} · **status** ("✅ sum of M segments = closed form to
  0.1 %" / "📏 → infinite line: Γ/2πd = 0.159 m/s within 0.5 %" / "⚠️ printed (5.14): the swirl is reversed — the correct
  factor is +1/(4π)") · **terms** (the per-segment contributions as bars in the Explain, grouped by quarter of the
  filament).
- **Readouts:** "Speed at x" (m/s) · "Closed form" (m/s) · "Distance d" (m) · "Infinite-line value" (m/s).
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**Filament and field point** (3-D): the <b class="c-teal">filament</b> cut into M pieces,
     each piece's push at x (<b class="c-accent">purple</b>) and their sum (black). **Contribution along the filament**
     unrolls the sum: the shaded area is the speed. **Speed vs distance** (hidden on phones) compares with the infinite
     line (grey)."
  1. *One piece* — (5.17) $d\mathbf u=\frac{\Gamma\,dl}{4\pi}\mathbf e_\omega\times\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$: for the clicked
     piece `live du` m/s — "perpendicular to the filament and to the line joining them, weakening as 1/distance²."
  2. *The sum* — Σ over M = `live M` pieces = **`live u`** m/s (boxed).
  3. *The closed form* (segment) — d = `live d` m, cos θ_a = `live ca`, cos θ_b = `live cb`:
     $\frac{\Gamma}{4\pi d}(\cos\theta_a-\cos\theta_b)$ = `live G`/(4π × `live d`) × (`live ca` − (`live cb`)) = **`live uc`** m/s; ring
     centre: Γ/2R = `live ur` m/s; on the axis ΓR²/2(R² + z²)^{3/2}.
  4. *The infinite line* — $\Gamma/2\pi d$ = `live ul` m/s: "(5.2) $u_\theta=\Gamma/2\pi r$ again — the book never closes this
     loop; D13 does."
  5. *Where the law comes from* — ∇²u = −∇×ω, Green's function −1/(4π\|x − x′\|) → (5.14) with +1/(4π) → Biot–Savart (5.16)
     $\mathbf u=\frac1{4\pi}\int\frac{\boldsymbol\omega\times(\mathbf x-\mathbf x')}{\lvert\mathbf x-\mathbf x'\rvert^3}d^3x'$ → thin tube (5.17).
  6. *The sign* (tube mode) — corrected u_θ(r = `live r`) = +`live ut` m/s; printed −1/(4π): `live utp` m/s — "the book's
     (5.14) prints −1/(4π); with it a vortex would turn against its own vorticity."
  7. *The distance view* (hint on phones).
  8. *Reading the current setting* — segment short: "Your point sees the segment's ends at wide angles: it gets only
     `pct` % of the infinite line's push." · long: "The ends are far: cos θ_a → 1, cos θ_b → −1, and the push is Γ/2πd." ·
     half-line: "Half of every infinite line: Γ/4πd — the downwash at a wing from its trailing vortices is half the value
     far behind." · ring: "Every piece of the ring pushes the centre the same way: Γ/2R; on the axis the push fades as
     1/z³ far away." · helix: "Nearby turns dominate; the sum wraps round the helix axis." · tube: "Inside the core the
     push is less than Γ/2πr because only the vorticity inside r counts (the planar shell theorem)."
- **Derivation tab:** **D10 ★★★** (12 steps) `view: 'unrolled'` (tube mode); goal `set {shape: 'tube', sign: false}`;
  step 3 shows the Poisson equation next to the tube's ∇×ω; step 7 **live** "flux of ∇(1/r) through a sphere: −4π for
  every radius"; step 11 `set {sign: true}` `watch` "the printed sign reverses every arrow"; step 12 `set {sign: false}`.
  **D11 ★★★** (12 steps) `view: 'scene'`; step 5 shows the book's printed line in rose beside the correct one in teal;
  steps 10–11 draw the volume V′ (dashed cylinder) round one segment with its ends ⟂ ω and its side outside the core
  (`watch` "the dashed volume V′ hugs the filament piece"); step 12 `set {shape: 'segment'}`. **D12** (6 steps) `view:
  'scene'`; step 1 `set {M: 8}`; step 3 **live** "3a/\|x − x′\| = 3 × 0.05/1 = 15 %: move the point farther"; step 6 `set
  {M: 64}`. **D13** (8 steps) `view: 'scene'`; step 1 `set {shape: 'segment', ell: 1, px: 1, pz: 0}`; step 5 draws θ at
  an element; step 7 **live** "(1/4π)(0.707 + 0.707) = 0.1125 m/s"; step 8 `set {ell: 20}`.
- **Code:**
  ```python
  poly = ch05.filament_preset("{{shape}}", M={{M}}, **{{geom}})
  x = np.array([{{px}}, 0.0, {{pz}}])               # the field point
  du = ch05.filament_contributions(x, poly, {{G}})  # one push per piece (5.17)
  print(du.sum(axis=1))                              # {{u}} m/s
  print(ch05.segment_speed({{d}}, {{ta}}, {{tb}}, {{G}}))   # closed form {{uc}}
  print({{G}}/(2*np.pi*{{d}}))                       # infinite line (5.2): {{ul}}
  ```
- **Walkthrough (7 steps):** 1. "Action at a distance" — "A vortex here, water moving over there. How much, and which
  way?" `set {shape: 'segment', ell: 1, px: 1, pz: 0}` · 2. "Velocity from vorticity" — "Curl ω = −∇²u: a Poisson equation.
  Its solution is Biot–Savart (5.16), with +1/(4π) in (5.14)." `derive: {id: 'D10', step: 11}`, `eq: 'bs'` · 3. "One
  piece" — "Each piece pushes as Γdl/4π · e_ω × (x − x′)/\|x − x′\|³ (5.17). Click one." `inspect: true`, `eq: 'fil'` ·
  4. "Add them up" — "The purple pushes sum to 0.1125 m/s — the closed form (Γ/4πd)(cos θ_a − cos θ_b)." `readouts:
  ['u', 'uc']`, `code: {lines: [4, 5]}` · 5. "An endless line" — "Stretch the segment: the speed climbs to Γ/2πd = 0.159
  m/s, the line vortex (5.2)." `set {ell: 20}`, `derive: {id: 'D13', step: 8}` · 6. "A ring" — "Bend it into a ring:
  every piece pushes the centre the same way, Γ/2R." `set {shape: 'ring', ell: 0.5, px: 0, pz: 0}` · 7. "Your turn" —
  "Predict: move the point to d = 0.5 m. Double or half? Then flip to the printed sign." `set {shape: 'segment', ell: 20}`,
  `controls: ['px', 'sign']`.
- **Equations:** `poisson` "Velocity's Poisson equation" (no ref) `\nabla\times\boldsymbol\omega=-\nabla^2\mathbf u` · `green`
  "Green's-function solution" ref 'Eq. (5.14)' `\mathbf u=+\frac{1}{4\pi}\int_{V'}\frac{\nabla'\times\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}d^3x'`
  note "the book prints −1/(4π); corrected here" · `curlgauss` "Gauss in curl form" ref 'Eq. (5.15)'
  `\int_{V'}\nabla'\times\Big(\frac{\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}\Big)d^3x'=\int_{A'}\frac{\mathbf n\times\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}d^2x'` ·
  `bs` "Biot–Savart" ref 'Eq. (5.16)' `\mathbf u=\frac1{4\pi}\int_{V'}\frac{\boldsymbol\omega\times(\mathbf x-\mathbf x')}{\lvert\mathbf x-\mathbf x'\rvert^3}d^3x'` ·
  `fil` "Filament law" ref 'Eq. (5.17)' `d\mathbf u=\frac{\Gamma\,dl}{4\pi}\mathbf e_\omega\times\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}` live ·
  `seg` "Straight segment" (no ref) `u=\frac{\Gamma}{4\pi d}(\cos\theta_a-\cos\theta_b)` live.
- **Check yourself:** (1) "Long segment: halve d. What happens to the speed?" — "It doubles: Γ/2πd." `set {shape:
  'segment', ell: 20, px: 0.5}` · (2) "Why does a semi-infinite line give exactly half?" — "Its ends are at θ = π/2 and 0:
  (Γ/4πd)(1 − 0)." `set {shape: 'half-line'}` · (3) "Ring with the point at the centre: which piece pushes most?" — "All
  equally — every piece is R away and pushes along the axis: Γ/2R." `set {shape: 'ring'}` · (4) "Printed sign: which way
  would the fluid circle a vortex whose vorticity points up?" — "Clockwise seen from above — against its own vorticity;
  impossible, so the printed −1/(4π) must be wrong." `set {shape: 'tube', sign: true}`.
- **Selftest parity rows:** `{name: 'segment', js: segVel([1, 0, 0], [0, 0, -1], [0, 0, 1], 1)[1], py:
  'ch05.segment_induced_velocity([1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 0.0, 1.0], 1.0)[1]', rtol: 1e-12}` · `{name:
  'segment speed', js: segSpeed(1, Math.PI/4, 3*Math.PI/4, 1), py: 'ch05.segment_speed(1.0, np.pi/4, 3*np.pi/4, 1.0)',
  rtol: 1e-12}` · `{name: 'square centre', js: filVel([0, 0, 0], polyline('square', 64, {side: 1}), 1, true)[2], py:
  'ch05.filament_velocity_preset("square", 0.0, 0.0, 0.0, side=1.0, component=2)', rtol: 1e-10}` · `{name: 'ring M=16', js:
  filVel([0, 0, 0], polyline('ring', 16, {R: 0.5}), 1, true)[2], py: 'ch05.filament_velocity_preset("ring", 0.0, 0.0,
  0.0, M=16, R=0.5, component=2)', rtol: 1e-10}` · `{name: 'ring axis', js: ringAxis(0.5, 0.5, 1), py:
  'ch05.ring_axis_velocity(0.5, 0.5, 1.0)', rtol: 1e-12}` · invariant `{name: 'long segment → line', js: segVel([1, 0, 0],
  [0, 0, -1e4], [0, 0, 1e4], 1)[1], expect: 1/(2*Math.PI), rtol: 1e-7}`.
- **Fit plan:** 360×640: `scene` + `unrolled`; status one line ("✅ 0.1125 = closed form"); shape chips short ("seg · half ·
  sq · ring · helix · tube"). Landscape phone: one row. Desktop/notebook: rows [1.3, 1], `distance` beside `unrolled`. The
  D10/D11 lines longer than the 1000 px side panel are split exactly as in Part F (steps 7, 10, 12 of D11 use two-row
  `aligned`). three.js fallback: an (x, z) side view.

### E7 · vorticity_equation_rotating
- **Title:** "What does a spinning planet add to the vorticity equation?" · **Summary:** "A column of fluid crossing a
  ridge, a ring of air moving poleward, and a term-bar budget of (5.30) in three scenes: on a rotating planet the relative
  vorticity changes, but (ζ + f)/h and the absolute circulation do not." · **CORE:** C09, C11 (also R12–R18, N30, N32,
  N35, N37, N39, N54) · **Reference:** `angular_frequency_explorer_1.html` (modes, "Right now" notes with a highlighted
  table) with `fid_formula_lab.html`'s term bars; link to ch04 `rotating_frame_coriolis`.
- **meta:** `viz:order 7` · `viz:sections 5.6` · `viz:equations 5.20 5.25 5.28 5.30 5.33` · `viz:fluidpy
  ch05.column_relative_vorticity ch05.column_over_slope ch05.relative_circulation_after_move ch05.vorticity_budget_preset
  ch05.planetary_vorticity_terms ch05.absolute_vorticity` · `viz:derivations D14 D15 D19 D20`.
- **Physics (JS ↔ Python):** `fOf(latDeg)` = 2Ω sin φ ↔ `core.rotating.coriolis_parameter` (via `ch05.column_over_slope`'s
  `f`) · `columnZeta(h, h0, z0, f)` ↔ `ch05.column_relative_vorticity` · `depthProfile(x, name, p)` ↔ `ch05.column_over_slope`
  (ridge / trough / slope) · `ringGamma(G0, A0, lat0, A1, lat1)` ↔ `ch05.relative_circulation_after_move` · `budget(scene,
  p)` ↔ `ch05.vorticity_budget_preset(scene, …)` (closed forms: Burgers terms; lock exchange baroclinic; rotating column
  planetary = 2Ωα).
- **Modes** (chips "column · ring · budget"): **column** (a layer over a bottom with a ridge or trough; the column moves
  along x with the clock) · **ring** (a ring of air on a globe cap moved from latitude φ₀ to φ₁ at constant area) ·
  **budget** (term bars of (5.30) for three scenes).
- **Views** (rows [1.3, 1]):
  1. `layer` "The rotating layer" — column mode: side view x ∈ [−3, 3]×10⁵ m, depth h(x) with the bottom shaded, the
     column (a vertical bar) at the clock's x, its spin shown by a turning marker on top (teal counterclockwise = cyclonic
     in the NH, rose clockwise), Ω arrow at the side; ring mode: an orthographic globe cap with the ring (teal) at its
     latitude and its relative circulation as an arrow; budget mode: a schematic of the scene (Burgers core / lock tank /
     stretched column).
  2. `ratio` "ζ vs h (or latitude)" — column: ζ against h (teal line ζ = (ζ₀ + f)h/h₀ − f) with the column's path as a dot
     moving along it, the amber horizontal line of (ζ + f)/h (flat) on a second axis; ring: Γ₁ against φ₁ (teal curve) with
     Γ_a (amber, flat); budget: the (5.30) term bars (local blue dashed, advective grey, relative stretching/tilting
     purple, planetary amber, baroclinic orange, diffusion rose, residual black).
  3. `table` "f at your latitude" (`hidePortrait: true`) — the table f(0°), f(30°), f(45°), f(60°), f(90°) with the current
     latitude's row lit and ζ/f for the column.
  Portrait: `layer` + `ratio`; the conserved ratio and f go into the `ratio` title ("(ζ + f)/h = 1.03×10⁻⁷ m⁻¹s⁻¹ · f =
  1.03×10⁻⁴ s⁻¹").
- **Controls (≤ 5 visible, per mode):** column — `lat` "Latitude φ" −90…90°, default 45 · `bump` "Ridge height" −400…400 m,
  default 200 (negative = trough) · `zeta0` "Initial ζ₀/f" −1…1, default 0 · `t` (transport: the column crosses the
  domain, 0 → 1, `end: 'hold'`); ring — `lat0` "From φ₀" 0…90°, default 30 · `lat1` "To φ₁" 0…90°, default 60 · `rad` "Ring
  radius" 100…2000 km, default 500; budget — `scene` chips burgers/lock/column · `Om` "Frame Ω" 0…1 rad/s, default 0.5
  (column scene).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (3) · **transport** (the column's
  crossing; end card "Back at h₀: ζ = ζ₀ again — nothing was lost") · **presets** "column over a ridge (45° N)" {mode:
  'column', bump: 200}, "column into a trough" {bump: -200}, "ring 30° → 60° N" {mode: 'ring', lat0: 30, lat1: 60}, "ring
  to the equator" {lat1: 0}, "Southern Hemisphere" {lat: -45}, "budget: stretched column" {mode: 'budget', scene: 'column'}
  · **terms** (budget mode) · **status** ("🌀 cyclonic spin-up: ζ = +2.06×10⁻⁵ s⁻¹" / "🔃 anticyclonic: ζ = −2.06×10⁻⁵ s⁻¹" /
  "⚖ (ζ + f)/h conserved to 1e-15" / ring "↻ anticyclonic Γ = −4.19×10⁷ m²/s") · **notes** (regime text) · the **table**
  of real values.
- **Readouts:** "f = 2Ω sin φ" (1/s) · "h" (m) · "ζ" (1/s) · "(ζ + f)/h" (1/(m s)) · ring: "Γ₁" (m²/s), "mean ζ" (1/s).
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**The rotating layer**: a column (bar) riding over a <b class="c-muted">bottom</b> with a
     ridge; its marker turns <b class="c-teal">cyclonically</b> or <b class="c-rose">anticyclonically</b>. **ζ vs h** plots
     the column's relative vorticity against its height, with the conserved (ζ + f)/h in <b class="c-amber">amber</b>. **f at
     your latitude** (hidden on phones) lists the planetary vorticity."
  1. *The planet's vorticity here* — f = 2Ω sin φ = 2 × 7.292×10⁻⁵ × sin(`live lat`°) = **`live f`** s⁻¹.
  2. *The column now* — h = h₀ − bump·shape(x) = `live h` m (h₀ = 1000 m).
  3. *Conservation* — from (5.33) $\frac{D\Gamma_a}{Dt}=0$ with A h fixed: $\frac{\zeta+f}{h}=\frac{\zeta_0+f}{h_0}$ = `live q` m⁻¹s⁻¹ (boxed).
  4. *Its relative vorticity* — ζ = (ζ₀ + f)h/h₀ − f = **`live z`** s⁻¹ = `live zf` f — sense: `live sense`.
  5. *The ring of air* (ring mode) — $\Gamma_1=\Gamma_0+2\Omega(A_0\sin\varphi_0-A_1\sin\varphi_1)$ = **`live G1`** m²/s; mean ζ = Γ₁/A = `live zr`
     s⁻¹.
  6. *The budget* (budget mode) — (5.30) $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u+\frac1{\rho^2}\nabla\rho\times\nabla p+\nu\nabla^2\boldsymbol\omega$
     at the probe: local `live lo`, advective `live ad`, relative stretching `live rs`, planetary `live pl`, baroclinic `live
     ba`, diffusion `live di`, residual `live re` s⁻².
  7. *The table* (hint on phones).
  8. *Reading the current setting* — squashed column: "Over the ridge the column is squashed: absolute vorticity must fall
     with h, so it turns anticyclonically — the reason highs sit over mountains." · stretched: "Into the trough the column
     is stretched and borrows spin from the planet: cyclonic." · ring: "Moving poleward, the ring encloses more planetary
     vorticity; to keep Γ_a it must spin the other way: anticyclonic Γ = `G1` m²/s." · Southern Hemisphere: "f < 0: every
     sense flips — cyclones turn clockwise." · budget: "Burgers: relative stretching and diffusion; lock exchange: only
     the baroclinic bar; rotating column: the planetary bar, 2Ω ∂w/∂z."
- **Derivation tab:** **D14** (8 steps) `view: 'ratio'` (budget mode); step 8 `set {mode: 'budget'}`. **D15 ★★★** (15
  steps) `view: 'ratio'`; goal `set {mode: 'budget', scene: 'column'}`; steps 5–9 light the purple+amber bar they build;
  step 11 lights the orange bar (`set {scene: 'lock'}`); step 13 the rose bar (`set {scene: 'burgers'}`); step 15 `set
  {scene: 'burgers'}` `watch` "planetary and baroclinic bars drop to zero: (5.13) again". **D19** (10 steps) `view:
  'layer'`; goal `set {mode: 'ring'}`; step 7 draws the ring's vector area; step 10 **live** "Γ = 2ΩA(sin 30° − sin 60°) =
  −4.19×10⁷ m²/s". **D20** (6 steps) `view: 'ratio'`; goal `set {mode: 'column', bump: -100}`; step 5 **live** "f =
  10⁻⁴ s⁻¹, h 1000 → 1100 m: ζ = 1.0×10⁻⁵ s⁻¹" with `set {lat: 43.3, bump: -100}` (f ≈ 10⁻⁴); step 6 shows the small-change
  tangent line.
- **Code:**
  ```python
  f = rotating.coriolis_parameter(np.deg2rad({{lat}}))   # {{f}} 1/s
  z = ch05.column_relative_vorticity({{h}}, 1000.0, {{z0}}, f)
  print(z, (z + f)/{{h}})                                # {{z}}, {{q}} (conserved)
  A = np.pi*({{rad}})**2
  print(ch05.relative_circulation_after_move(0.0, A, {{lat0}}, A, {{lat1}}))  # {{G1}}
  print(ch05.vorticity_budget_preset("{{scene}}"))       # (5.30) terms
  ```
- **Walkthrough (7 steps):** 1. "Spin from the planet" — "Air crossing mountains, converging into a low: where does its
  spin come from?" `set {mode: 'column', bump: 200, t: 0}` · 2. "The full equation" — "In a rotating, stratified fluid,
  Dω/Dt = (ω + 2Ω)·∇u + ∇ρ×∇p/ρ² + ν∇²ω (5.30)." `eq: 'e530'`, `derive: {id: 'D15', step: 15}` · 3. "Stretch the planet's
  spin" — "Keep one term: Dω_z/Dt = 2Ω ∂w/∂z. A stretched column gains ζ." `set {mode: 'budget', scene: 'column'}`,
  `terms: true` · 4. "Absolute circulation" — "Kelvin in a rotating frame: Γ_a = Γ + 2∫Ω·n dA is conserved (5.33)."
  `set {mode: 'ring'}`, `eq: 'e533'` · 5. "The column rule" — "Thin loop + fixed mass: (ζ + f)/h is constant. Press ▶: the
  column crosses the ridge." `set {mode: 'column', t: 0}`, `play: true`, `derive: {id: 'D20', step: 4}` · 6. "Read it" —
  "Squashed over the ridge: ζ = −0.2f, anticyclonic. Back at h₀: ζ = 0." `set {t: 0.5}`, `play: false`, `readouts: ['z',
  'q']`, `code: {lines: [2, 3]}` · 7. "Your turn" — "Predict: the same ridge at 15° N — bigger or smaller ζ? Then drag
  the latitude." `controls: ['lat']`.
- **Equations:** `e520` "Rotating momentum" ref 'Eq. (5.20)'
  `\frac{\partial u_i}{\partial t}+u_ju_{i,j}+2\varepsilon_{ijk}\Omega_ju_k=-\frac1\rho p_{,i}+g_i+\nu u_{i,jj}` · `e525` "Lamb form" ref 'Eq. (5.25)'
  `\begin{aligned}&\partial_tu_i+(\tfrac12u_j^2+\Phi)_{,i}-\varepsilon_{ijk}u_j(\omega_k+2\Omega_k)\\&=-(1/\rho)p_{,i}-\nu\varepsilon_{ijk}\omega_{k,j}\end{aligned}` ·
  `e528` "Baroclinic term" ref 'Eq. (5.28)' `\frac1{\rho^2}[\nabla\rho\times\nabla p]_n` · `e530` "Vorticity equation, rotating" ref
  'Eq. (5.30)'
  `\begin{aligned}\tfrac{D\boldsymbol\omega}{Dt}&=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u\\&\quad+\tfrac1{\rho^2}\nabla\rho\times\nabla p+\nu\nabla^2\boldsymbol\omega\end{aligned}` live ·
  `e533` "Absolute circulation" ref 'Eq. (5.33)'
  `\frac{D\Gamma_a}{Dt}=0,\ \Gamma_a=\Gamma+2\int_A\boldsymbol\Omega\cdot\mathbf n\,dA` live · `col` "Fluid column" (no ref)
  `\frac{\zeta+f}{h}=\text{const}` live.
- **Check yourself:** (1) "At the equator, what does the ridge do to the column?" — "Nothing: f = 0 at the equator and
  ζ₀ = 0, so ζ = 0·h/h₀ = 0." `set {mode: 'column', lat: 0}` · (2) "A column starting with ζ₀ = −f (zero absolute
  vorticity) crosses the ridge. ζ?" — "ζ = −f all the way: (ζ + f)/h = 0 stays 0." `set {zeta0: -1}` · (3) "Move the ring from 30° N to the equator. Sign of Γ₁?" — "Positive (cyclonic): it
  encloses less planetary vorticity than before." `set {mode: 'ring', lat1: 0}` · (4) "Budget, lock-exchange scene: why is
  only one bar nonzero?" — "At t = 0 the fluid is at rest: no stretching, advection or diffusion — only ∇ρ×∇p/ρ²." `set
  {mode: 'budget', scene: 'lock'}`.
- **Selftest parity rows:** `{name: 'column 1100', js: columnZeta(1100, 1000, 0, 1e-4), py:
  'ch05.column_relative_vorticity(1100.0, 1000.0, 0.0, 1e-4)', rtol: 1e-12}` · `{name: 'ring 30→60', js: ringGamma(0,
  Math.PI*25e10, 30, Math.PI*25e10, 60), py: 'ch05.relative_circulation_after_move(0.0, np.pi*2.5e11, 30.0, np.pi*2.5e11,
  60.0)', rtol: 1e-10}` · `{name: 'f45', js: fOf(45), py: 'ch05.column_over_slope(0.0, lat_deg=45.0)["f"]', rtol: 1e-12}` ·
  `{name: 'planetary budget', js: budget('rotating_column').planetary, py:
  'ch05.vorticity_budget_preset("rotating_column")["planetary"]', rtol: 1e-6}` · `{name: 'lock baroclinic', js:
  budget('lock_exchange').baroclinic, py: 'ch05.vorticity_budget_preset("lock_exchange")["baroclinic"]', rtol: 1e-6}` ·
  invariant `{name: 'ratio conserved', js: (columnZeta(800, 1000, 0, 1e-4) + 1e-4)/800 - 1e-4/1000, expect: 0, atol:
  1e-18}`.
- **Fit plan:** 360×640: `layer` + `ratio`; status one line ("🔃 anticyclonic · ζ −0.2f"); mode chips "column · ring ·
  budget". Landscape phone: one row. Desktop/notebook: rows [1.3, 1], `table` beside `ratio`. The D15 lines use two-row
  `aligned` blocks (steps 1, 5, 14) to fit the 1000×700 frame.

### E8 · point_vortex_lab
- **Title:** "A vortex cannot push itself — so how do vortices move?" · **Summary:** "Drop vortices of either sign, drag
  them and press play: like vortices orbit their fixed centre of vorticity, opposite ones march off together, a vortex
  near a wall slides along it with its mirror twin — and the invariants stay flat." · **CORE:** C12, C13 (also N40, N41,
  N55, N56, N57, N58) · **Reference:** `random_copy_lab.html` (click-to-place interaction, code comments showing the live
  values) with `angular_frequency_explorer_1.html`'s end-of-run summary.
- **meta:** `viz:order 8` · `viz:sections 5.7` · `viz:equations 5.2` · `viz:fluidpy ch05.point_vortex_velocity
  ch05.point_vortex_rhs ch05.point_vortex_evolve ch05.point_vortex_invariants ch05.centre_of_vorticity
  ch05.wall_image_system ch05.circle_image_system ch05.vortex_pair ch05.vortex_near_wall_speed ch05.point_vortex_preset` ·
  `viz:derivations D21 D22`.
- **Physics (JS ↔ Python):** `rhs(xv, G, boundary)` ↔ `ch05.point_vortex_rhs` (images: wall mirror, circle inverse points
  with `inside: true`) · `step(xv, G, dt)` (RK4 with dt = 0.01 s and 4 sub-steps when two vortices come within 0.1 m;
  mirrors the DOP853 `point_vortex_evolve` to 1e-6 over one period) · `invariants(xv, G)` ↔ `ch05.point_vortex_invariants` ·
  `pair(G1, G2, h)` ↔ `ch05.vortex_pair` · `wallSpeed(G, h)` ↔ `ch05.vortex_near_wall_speed` · `preset(name)` ↔
  `ch05.point_vortex_preset`.
- **Views** (rows [1.35, 0.9]):
  1. `plane` "The vortices" (`equal: true`, x, y ∈ [−2, 2] m; circle mode shows the bucket of radius 1 m) — vortices as
     discs (teal Γ > 0, rose Γ < 0, area ∝ \|Γ\|), their velocity arrows, trails (faint), the centre of vorticity G (black
     cross; hidden when ΣΓ = 0), the wall (grey hatched) or circle, images faint behind the boundary, a faint streamline
     field (optional toggle). Pointer: click empty space → add a vortex of the chosen sign; drag a vortex (pauses the
     clock); click a vortex → inspector; right-click/long-press → delete.
  2. `inv` "Invariants" — P_x, P_y (ΣΓx, ΣΓy), I = ΣΓ\|x\|² and H against t, each normalised to its start (flat lines at
     1); unbounded plane only (with a wall: P_x and H of the vortex + image system).
  3. `pairview` "Pair distance and angle" (`hidePortrait: true`) — for two vortices: the separation (flat) and the angle of
     the joining line vs t with the predicted slope (Γ₁ + Γ₂)/2πh² (grey dashed ghost); for a wall: the height (flat) and
     the drift distance with slope Γ/4πh.
  Portrait: `plane` + `inv`; the orbit rate or drift speed in the `inv` title ("rate 0.318 rad/s · predicted 0.318").
- **Controls (≤ 5 visible):** `preset` chips "equal · 1:3 · opposite · wall · bucket · three" · `ratio` "Γ₂/Γ₁" −3…3,
  default 1 (pair presets) · `h` "Separation h (or wall distance)" 0.2…2 m, default 1 · `boundary` chips "none · wall ·
  circle" · `t` (transport 0 → 40 s, rate 2, `end: 'hold'`; end card "G never moved: ΣΓx conserved to 1e-12" or "The pair
  marched 6.4 m at 0.159 m/s").
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** (six) ·
  **inspector** (click a vortex k: "from vortex j: Γ_j/2π · e_z × (x_k − x_j)/\|x_k − x_j\|² = (…) m/s; … sum = (…) m/s — its
  velocity") · **status** ("🔄 co-rotating about G: period 19.7 s" / "➡ translating pair: 0.159 m/s" / "🧱 sliding along
  the wall: Γ/4πh = 0.159 m/s" / "🪣 in a bucket: images outside the wall") · **notes**.
- **Readouts:** "Orbit rate" (rad/s) · "Pair speed" (m/s) · "Wall drift" (m/s) · "G" (m).
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**The vortices**: <b class="c-teal">counterclockwise</b> and <b class="c-rose">clockwise</b>
     vortices with their velocity arrows and trails; the black cross is the centre of vorticity G; faint discs behind a
     wall are images. **Invariants** must stay flat — a check on the computation. **Pair distance and angle** (hidden on
     phones) compares the motion with D21."
  1. *Each vortex's velocity* — $\frac{d\mathbf x_k}{dt}=\sum_{j\ne k}\frac{\Gamma_j}{2\pi}\frac{\mathbf e_z\times(\mathbf x_k-\mathbf x_j)}{\lvert\mathbf x_k-\mathbf x_j\rvert^2}$ — the
     line vortex (5.2) $u_\theta=\Gamma/2\pi r$ of every other vortex; "a straight vortex never pushes itself."
  2. *Two vortices* — V₁ = Γ₁/2πh = `live V1`, V₂ = Γ₂/2πh = `live V2` m/s (boxed).
  3. *The centre of vorticity* — h₁ = Γ₂h/(Γ₁ + Γ₂) = **`live h1`** m from vortex 1; "fixed, because ΣΓ_kx_k never changes;
     ⚠️ the fluid at G is at rest only if Γ₁ = Γ₂ (the book's Fig. 5.11 caption says otherwise)."
  4. *The orbit* — rate (Γ₁ + Γ₂)/2πh² = **`live rate`** rad/s, period `live T` s; opposite pair: translation Γ/2πh = `live V`
     m/s.
  5. *A wall* — image −Γ at the mirror point; drift Γ/(2π·2h) = Γ/4πh = **`live Vw`** m/s; normal velocity on the wall 0.
  6. *Invariants now* — P = (`live Px`, `live Py`), I = `live I`, H = `live H` — changes since t = 0: `live dI`, `live dH`.
  7. *The pair view* (hint on phones).
  8. *Reading the current setting* — equal: "Each vortex drags the other round: a dumbbell turning about its middle." ·
     1:3: "The stronger vortex moves less: G sits ¾ of the way toward it; the pair turns twice as fast." · opposite: "The
     pushes point the same way: no turning, a straight march — an aircraft's wake pair." · wall: "The wall's mirror twin
     does the pushing: the vortex slides along, never toward, the wall." · bucket: "The pair crosses, meets the wall and
     splits: each vortex is then carried by its own image (Lighthill's knife)." · three: "Three vortices can dance in
     complicated ways; the invariants still hold."
- **Derivation tab:** **D21** (6 steps) `view: 'plane'`; goal `set {preset: 'equal', t: 0}`; step 3 draws both velocity
  arrows ⟂ the joining line; step 4 `set {preset: 'unequal', ratio: 3}` marks G (**live** "h₁ = 0.75 m"); step 5 **live**
  "0.318 rad/s, period 19.7 s"; step 6 `set {preset: 'opposite'}`. **D22** (5 steps) `view: 'plane'`; goal `set {preset:
  'wall', h: 0.5}`; step 1 draws the image; step 3 draws wall arrows (`watch` "the wall arrows lie flat along the wall");
  step 5 **live** "Γ/4πh = 1/(4π × 0.5) = 0.159 m/s".
- **Code:**
  ```python
  P = ch05.point_vortex_preset("{{preset}}")
  t = np.linspace(0, {{t}}, 200)
  X = ch05.point_vortex_evolve(P["xv"], P["Gamma"], t, boundary={{bnd}})
  print(ch05.vortex_pair({{G1}}, {{G2}}, {{h}}))        # rate {{rate}}, centre {{h1}}
  inv = ch05.point_vortex_invariants(X[-1], P["Gamma"])
  print(inv["I"], inv["H"])                             # {{I}}, {{H}} (unchanged)
  ```
- **Walkthrough (6 steps):** 1. "Who pushes a vortex?" — "A straight vortex never moves itself. It goes where the others
  carry it." `set {preset: 'equal', t: 0}` · 2. "Like vortices orbit" — "Press ▶: each pushes the other sideways at
  Γ/2πh; they circle their centre G." `play: true`, `derive: {id: 'D21', step: 5}` · 3. "Unequal strengths" — "Γ₂ = 3Γ₁: G
  moves to ¾ of the way; the orbit is twice as fast." `set {preset: 'unequal', ratio: 3, t: 0}`, `play: false`, `readouts:
  ['rate', 'G']` · 4. "Opposite vortices" — "Pushes point the same way: the pair marches at Γ/2πh." `set {preset: 'opposite'}`,
  `code: {lines: [4, 4]}` · 5. "A wall is a mirror" — "Put an opposite twin behind the wall: no flow through it; the vortex
  slides at Γ/4πh." `set {preset: 'wall', h: 0.5}`, `derive: {id: 'D22', step: 5}` · 6. "Your turn" — "Predict: halve h in
  the equal pair — how much faster? Then click to add a third vortex." `set {preset: 'equal'}`, `controls: ['h']`.
- **Equations:** `pv` "Point-vortex motion" (no ref)
  `\frac{d\mathbf x_k}{dt}=\sum_{j\ne k}\frac{\Gamma_j}{2\pi}\frac{\mathbf e_z\times(\mathbf x_k-\mathbf x_j)}{\lvert\mathbf x_k-\mathbf x_j\rvert^2}` · `lv` "Line vortex" ref
  'Eq. (5.2)' `u_\theta=\Gamma/2\pi r` · `pair` "Pair" (no ref) `V_1=\Gamma_1/2\pi h,\ V_2=\Gamma_2/2\pi h` live · `centre` "Centre of vorticity"
  (no ref) `h_1=\Gamma_2h/(\Gamma_1+\Gamma_2),\ \dot\theta=(\Gamma_1+\Gamma_2)/2\pi h^2` live · `wall` "Wall drift" (no ref)
  `V_A=\Gamma/4\pi h` live.
- **Check yourself:** (1) "Equal pair: halve h. Orbit period?" — "Quartered: rate ∝ 1/h²." `set {preset: 'equal', h: 0.5}` ·
  (2) "Γ₂ = −Γ₁: where is G?" — "At infinity (ΣΓ = 0): the pair translates instead of turning." `set {preset: 'opposite'}`
  · (3) "Wall preset: why does the vortex not move toward the wall?" — "The image's push at A is perpendicular to AB, i.e.
  along the wall." `set {preset: 'wall'}` · (4) "Bucket: why do the two vortices part at the wall?" — "Each is carried
  along the wall by its own image, in opposite directions (opposite spins)." `set {preset: 'bucket'}`.
- **Selftest parity rows:** `{name: 'pair rate', js: pair(1, 1, 1).rate, py: 'ch05.vortex_pair(1.0, 1.0, 1.0)["rotation_rate"]',
  rtol: 1e-12}` · `{name: 'centre 1:3', js: pair(1, 3, 1).h1, py: 'ch05.vortex_pair(1.0, 3.0, 1.0)["centre_from_1"]', rtol:
  1e-12}` · `{name: 'wall drift', js: wallSpeed(1, 0.5), py: 'ch05.vortex_near_wall_speed(1.0, 0.5)', rtol: 1e-12}` · `{name:
  'rhs three', js: rhs(preset('three_vortices').xv, [1, 1, -0.5], null)[0][0], py:
  'ch05.point_vortex_rhs(ch05.point_vortex_preset("three_vortices")["xv"], [1.0, 1.0, -0.5])[0][0]', rtol: 1e-12}` · `{name:
  'H pair', js: invariants([[-0.5, 0.5], [0, 0]], [1, 1]).H, py: 'ch05.point_vortex_invariants([[-0.5, 0.5], [0.0, 0.0]],
  [1.0, 1.0])["H"]', rtol: 0, atol: 1e-15}` · invariant `{name: 'period return', js: periodReturnError(), expect: 0, atol:
  1e-6}`.
- **Fit plan:** 360×640: `plane` + `inv`; status one line ("🔄 T 19.7 s"); preset chips short, two rows; add/delete by tap
  and long-press on phones. Landscape phone: one row. Desktop/notebook: rows [1.35, 0.9], `pairview` beside `inv`.

### E9 · vortex_sheet_rollup
- **Title:** "What is a velocity jump made of?" · **Summary:** "A row of N filaments becomes a clean velocity jump whose
  size is the sheet's strength; a thin circuit round it measures γ = u₂ − u₁; ripple the sheet, press play and it rolls
  up." · **CORE:** C14 (also N44, N60; C12's point-vortex motion; forward link to Ch. 11) · **Reference:**
  `overfitting_curves.html` (a minimal two-slider figure with a verdict line) with `forced_damped_vibrations.html`'s
  explanation panel.
- **meta:** `viz:order 9` · `viz:sections 5.8` · `viz:equations 3.18` · `viz:fluidpy ch05.vortex_sheet_strength
  ch05.discrete_sheet_u ch05.continuous_sheet_velocity ch05.discrete_sheet_convergence ch05.vortex_sheet_velocity
  ch05.sheet_rollup` · `viz:derivations D23`.
- **Physics (JS ↔ Python):** `sheetU(x, y, g, N, L)` ↔ `ch05.discrete_sheet_u` · `contU(x, y, g, L)` ↔
  `ch05.continuous_sheet_velocity` · `strength(uA, uB, conv)` ↔ `ch05.vortex_sheet_strength` · `l1err(g, N)` ↔
  `ch05.discrete_sheet_convergence` (trapezoid on 4001 points of y ∈ [−0.1, 0.1] m) · `rollupStep(z, g, delta, L, dt)` (RK4
  with the Krasny-smoothed periodic kernel u − iv = Σ(Γ_j/2iL)cot(π(z − z_j)/L), smoothing δ) ↔ `ch05.sheet_rollup`
  (parity at t = 0.5 s, rtol 1e-4).
- **Modes** (chips "jump · roll-up"): **jump** (a finite sheet of length 1 m, N filaments, the circuit) · **roll-up** (a
  periodic sheet of period 1 m with a ripple, point vortices advancing in time).
- **Views** (rows [1.25, 1]):
  1. `sheet` "The sheet" — jump mode: x ∈ [−0.6, 0.6] m, y ∈ [−0.25, 0.25] m; the N filaments as small circular arrows
     (teal ccw for γ > 0), tangential velocity arrows above (u₁, pointing left) and below (u₂, pointing right) at a few
     heights, the draggable dn × ds circuit (black rectangle; its sides coloured as in D23: bottom +u₂ds, top −u₁ds, sides
     grey); roll-up mode: the sheet as a line through its points (teal), trails, one period with ghosts of the neighbours.
     Pointer: jump — drag the circuit's height and width, click a height → inspector; roll-up — none.
  2. `profile` "u(y) across the sheet" — u(0, y) of the discrete row (teal) and of the continuous sheet (black), y ∈
     [−0.2, 0.2] m, the jump ±γ/2 marked; roll-up: the sheet's centre-point displacement vs t (growth curve, log-y) with
     the linear-theory ghost e^{γkt/2} (k = 2π/L), labelled "preview of Ch. 11".
  3. `conv` "Circuit and convergence" (`hidePortrait: true`) — jump mode: the circuit's circulation vs its height dn (flat
     at γ ds once dn spans the sheet) and the L1 error vs N on log-log axes (slope −1, three reference dots).
  Portrait: `sheet` + `profile`; the circuit's Γ and the jump go in the profile title ("jump 2.000 = γ · circuit 0.200 =
  γ ds").
- **Controls (≤ 5 visible):** `N` "Filaments N" 4…1000, log-int, default 100 · `gamma` "Strength γ" 0.1…5 m/s, default 2 ·
  `ds` "Circuit width ds" 0.02…0.5 m, default 0.1 (jump) · `amp` "Ripple" 0…0.05 m, default 0.01 (roll-up) · `delta`
  "Smoothing δ" 0.02…0.2, default 0.05 (roll-up) · `conv` toggle "Caption convention (u₁ − u₂)" (optional) · transport
  `t` 0 → 3 s (roll-up), rate 0.5, `end: 'hold'`, end card "The ripple grew ×… and rolled into a spiral: the
  Kelvin–Helmholtz instability (Ch. 11)".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** (roll-up) · **modes** (2) ·
  **presets** "fine sheet N = 1000" {mode: 'jump', N: 1000}, "coarse N = 10" {N: 10}, "rippled sheet" {mode: 'rollup', amp:
  0.01}, "caption convention" {conv: true} · **status** ("✅ jump = γ to 0.1 %" / "〰 coarse row: within one spacing the
  profile wiggles" / "🌀 rolling up" / "🔁 caption convention: γ = u₁ − u₂ = −2 m/s — same physics") · **inspector** (click
  a height y: "u(0, y) = Σ_k (γ ds/2π)(−y)/(x_k² + y²) = (N terms) = … m/s; the continuous sheet −(γ/π)arctan(L/2y) = …").
- **Readouts:** "u just above" (m/s) · "u just below" (m/s) · "Jump u₂ − u₁" (m/s) · "Circuit Γ" (m²/s).
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**The sheet**: <b class="c-teal">N filaments</b> of strength γ ds each, the velocity above
     and below, and the thin circuit. **u(y)** compares the row (teal) with a continuous sheet (black). **Circuit and
     convergence** (hidden on phones) shows the circuit's circulation and how fast the row becomes a sheet."
  1. *Each filament* — Γ_k = γ ds = `live g` × `live dsN` = `live Gk` m²/s (ccw positive).
  2. *Just above and just below* — continuous sheet: u₁ = −γ/2 = `live u1`, u₂ = +γ/2 = `live u2` m/s; the row at y = ±5 cm:
     `live ua`, `live ub` m/s.
  3. *The jump* — u₂ − u₁ = **`live jump`** m/s = γ (boxed).
  4. *The circuit* — $d\Gamma=u_2\,ds+v\,dn-u_1\,ds-v\,dn=(u_2-u_1)\,ds$ = `live jump` × `live ds` = **`live Gc`** m²/s; "the two short sides
     cancel because v is continuous."
  5. *Convergence* — L1 error of the row's profile = `live err` (∝ 1/N: `live err10`, `live err100`, `live err1000` for N = 10,
     100, 1000).
  6. *The convention* — text: γ = u₂ − u₁ (counterclockwise); Fig. 5.16's caption: u₁ − u₂ = `live cap` — "the caption
     counts clockwise, the sense of its drawn filaments."
  7. *Roll-up* (roll-up mode; hint otherwise) — t = `live t` s, ripple amplitude `live A` m (×`live grow`); "vortices
     crowd where the sheet is already bunched: the sheet is unstable."
  8. *Reading the current setting* — fine: "A thousand filaments already look like a sheet: the velocity jumps from +1 to
     −1 m/s across it — the sheet's strength." · coarse: "With ten filaments the row wiggles within a spacing of itself;
     farther away it is indistinguishable from a sheet." · roll-up: "The sheet's own velocity makes the ripple grow; it
     winds into cat's eyes — how shear layers and wakes roll up into billows."
- **Derivation tab:** **D23** (5 steps) `view: 'sheet'`; goal `set {mode: 'jump', N: 100}`; steps 1–3 light the circuit's
  sides one by one; step 4 shrinks dn (`watch` "only the jump survives"); step 5 **live** "γ = 2 m/s: u₂ = +1, u₁ = −1 m/s"
  then `set {conv: true}` `watch` "the caption convention flips the sign only".
- **Code:**
  ```python
  print(ch05.discrete_sheet_u(0.0, 0.05, {{g}}, {{N}}))   # the row: {{ua}} m/s
  print(ch05.continuous_sheet_velocity(0.0, 0.05, {{g}})) # the sheet: {{uc}}
  print(ch05.vortex_sheet_strength({{u1}}, {{u2}}))       # u2 − u1 = {{jump}}
  c = ch05.discrete_sheet_convergence({{g}})
  print(c["l1_error"])                                    # ∝ 1/N
  R = ch05.sheet_rollup(N=100, gamma={{g}}, amplitude={{amp}})
  ```
- **Walkthrough (6 steps):** 1. "A jump in velocity" — "Fast stream over slow stream: the velocity jumps across a thin
  layer. What is the layer?" `set {mode: 'jump', N: 10}` · 2. "A row of vortices" — "Line up N small vortices: each
  pushes left above it, right below. Raise N." `controls: ['N']`, `set {N: 100}` · 3. "The jump" — "Just above −γ/2, just
  below +γ/2: the jump u₂ − u₁ = γ." `readouts: ['ua', 'ub', 'jump']` · 4. "Measure it" — "Circulation of a thin box:
  dΓ = u₂ds + v dn − u₁ds − v dn = (u₂ − u₁)ds." `derive: {id: 'D23', step: 4}`, `eq: 'circ'`, `code: {lines: [3, 3]}` ·
  5. "It rolls up" — "Ripple the sheet and press ▶: its own velocity wraps it into a spiral." `set {mode: 'rollup', amp:
  0.01, t: 0}`, `play: true` · 6. "Your turn" — "Predict: double γ — does the jump double, and does roll-up go faster?"
  `set {mode: 'jump'}`, `controls: ['gamma']`.
- **Equations:** `jumps` "Jumps across a sheet" (no ref) `[u_t]\neq0,\quad[v_n]=0` · `circ` "Circuit" ref 'Eq. (3.18)'
  `d\Gamma=u_2\,ds+v\,dn-u_1\,ds-v\,dn=(u_2-u_1)\,ds` live · `strength` "Sheet strength" (no ref) `\gamma=\frac{d\Gamma}{ds}=u_2-u_1` live ·
  `far` "Far from the middle" (no ref) `u=\mp\gamma/2` live.
- **Check yourself:** (1) "Coarse row, N = 10: where does the profile differ from the sheet?" — "Only within about one
  spacing (0.1 m) of the row; farther away the sum is already smooth." `set {N: 10}` · (2) "Make the circuit taller than
  the sheet's influence. Does its circulation change?" — "No — once it spans the sheet, only the jump counts: γ ds." `set
  {ds: 0.2}` · (3) "With the caption convention, is the physics different?" — "No: it counts circulation clockwise, so every
  sign flips; the flow is the same." `set {conv: true}` · (4) "Roll-up: what does a larger δ do?" — "It smooths the kernel
  like a thicker sheet: slower, fatter roll-up." `set {mode: 'rollup', delta: 0.15}`.
- **Selftest parity rows:** `{name: 'row N=100', js: sheetU(0, 0.05, 2, 100, 1), py: 'ch05.discrete_sheet_u(0.0, 0.05, 2.0,
  100)', rtol: 1e-12}` · `{name: 'row N=10', js: sheetU(0, 0.05, 2, 10, 1), py: 'ch05.discrete_sheet_u(0.0, 0.05, 2.0, 10)',
  rtol: 1e-12}` · `{name: 'continuous', js: contU(0, 0.05, 2, 1)[0], py: 'ch05.continuous_sheet_velocity(0.0, 0.05, 2.0)[0]',
  rtol: 1e-12}` · `{name: 'strength ccw', js: strength(-1, 1, 'ccw'), py: 'ch05.vortex_sheet_strength(-1.0, 1.0)', rtol:
  1e-12}` · `{name: 'strength caption', js: strength(-1, 1, 'caption'), py: 'ch05.vortex_sheet_strength(-1.0, 1.0,
  convention="caption")', rtol: 1e-12}` · `{name: 'L1 N=100', js: l1err(2, 100), py:
  'ch05.discrete_sheet_convergence(2.0)["l1_error"][1]', rtol: 1e-3}` · invariant `{name: 'jump = gamma', js:
  contU(0, -1e-7, 2, 1)[0] - contU(0, 1e-7, 2, 1)[0], expect: 2, rtol: 1e-6}`.
- **Fit plan:** 360×640: `sheet` + `profile`; status one line ("✅ jump 2.000 = γ"); in roll-up mode the `profile` view
  shows the growth curve. Landscape phone: one row. Desktop/notebook: rows [1.25, 1], `conv` beside `profile`.

### B1 · vortex_rings
- **Title:** "Why does a smoke ring widen at a wall — and play leap-frog?" · **Summary:** "Axisymmetric thin-core rings on
  one clock: a ring moved by its own curvature is widened and slowed by its image in a wall; two coaxial rings take turns
  passing through each other while their total impulse stays fixed." · **CORE:** C13, C08 (also N42, N43, N59) ·
  **Reference:** `angular_frequency_explorer_1.html` (linked views on one clock, end-of-run summary). **Backup only — not
  embedded unless one of E1–E9 fails review.**
- **meta:** `viz:order 10` · `viz:sections 5.5 5.7` · `viz:equations 5.17` · `viz:fluidpy ch05.ring_ring_velocity
  ch05.ring_self_velocity ch05.ring_axis_velocity ch05.ring_dynamics` · `viz:derivations none` (D12 is linked from the
  Explain).
- **Physics:** `ringSelf(R, a, G)` ↔ `ch05.ring_self_velocity` (Γ/4πR[ln(8R/a) − ¼]) · `ringRing(R, z, R0, z0, G0)` with
  `ellipK(m)`, `ellipE(m)` by the arithmetic–geometric mean (**parameter m = k²**) ↔ `ch05.ring_ring_velocity` · `evolve(
  rings, dt)` (RK4, core volume a²R conserved) ↔ `ch05.ring_dynamics`.
- **Views:** (1) `meridional` "Rings in the (R, z) plane" — cross-sections (teal ccw about +z), their images below the wall
  (faint), trails; (2) `paths` "R(t), z(t)" with the impulse ΣΓπR² (flat, amber); (3) `speed` (`hidePortrait`) — approach
  speed vs distance to the wall.
- **Controls:** `mode` chips "wall · leap-frog" · `Gamma` 0.1…5 m²/s · `R0` 0.2…2 m · `a` core 0.02…0.3 m · `sep` start
  separation 0.3…2 m (leap-frog) · transport t 0 → 60 s.
- **Depth features:** Explain, Code, linked views (3), transport (end card "Impulse conserved to 1e-9; 3 passes"),
  presets ("ring meets a wall", "leap-frog", "fat cores"), status ("🧱 widening, slowing: z → wall, never reaching it" /
  "🐸 pass #2 at t = … s").
- **Explain** (outline): 0 what the views show · 1 self-speed with numbers (0.329 m/s for Γ = R = 1, a = 0.1) · 2 the
  image ring's push at A (outward and backward) · 3 the other ring's push (leap-frog) · 4 impulse · 5 reading the current
  setting (wall: "the image pushes the ring outward along the wall and holds it back"; leap-frog: "the rear ring shrinks,
  speeds up and passes through").
- **Walkthrough (5 steps):** "A smoke ring moves itself" · "Its speed Γ/4πR[ln(8R/a) − ¼]" · "Near a wall: an image ring"
  · "Two rings: leap-frog" · "Your turn: change the separation".
- **Check yourself (3):** fatter core → slower ring (ln(8R/a) smaller) · a ring never touches the wall (the image's push
  grows as it approaches) · leap-frog conserves ΣΓπR² (a wider front ring means a narrower rear one).
- **Selftest parity rows:** `{name: 'self speed', js: ringSelf(1, 0.1, 1), py: 'ch05.ring_self_velocity(1.0, 0.1, 1.0)',
  rtol: 1e-12}` · `{name: 'ring-ring on axis', js: ringRing(0, 0.5, 0.5, 0, 1)[1], py: 'ch05.ring_ring_velocity(0.0, 0.5,
  0.5, 0.0, 1.0)[1]', rtol: 1e-10}` · `{name: 'ring-ring off axis', js: ringRing(0.8, 0.3, 1.0, 0.0, 1)[0], py:
  'ch05.ring_ring_velocity(0.8, 0.3, 1.0, 0.0, 1.0)[0]', rtol: 1e-9}`.
- **Fit plan:** as E8 (two views on phones, the speed view hidden).

## Part D — runtime budget (full run < 5 min on a laptop / Colab CPU)

Chapter 5 is mostly closed forms and point stencils (µs–ms); the costs are the animations (≈ 0.12 s per frame at dpi 80,
ch01–ch04 timings), material-loop and point-vortex integrations (`solve_ivp` DOP853), the vortex-sheet roll-up (O(N²)
right-hand side), the 3-D Biot–Savart and Green quadratures (24³–40³ nodes, vectorised), the planar Biot–Savart grid sums,
the plotly slider figures that precompute every step, and four ★★★ sympy checks (D09, D10, D11, D15). The ch04 notebook
ran in 89 s; this one has 14 CORE blocks, 23 derivations and 9 explainers and must stay under 5 min.

| Section | Heaviest cells | Full | FAST (`FLUIDPY_FAST=1`) |
|---|---|---|---|
| setup + imports | numpy/scipy/sympy/plotly/pint | 9 s | 9 s |
| §5.1 C01 | vortex line (`solve_ivp`), tube strength (512-point loop, 32 × 64 disc), tube sections and budgets (Gauss–Legendre 48 × 48), from-scratch lid sums (200 × 400 × 2), plotly 3-D (16 streamlines + 16 vortex lines × 200 points, two surfaces 40 × 40), profile figure | 7 s | 5 s (24 × 24 lids, 8 + 8 lines) |
| §5.1 C02 | pressure closed forms, tank `brentq`, sympy (5.5) check, `cumulative_trapezoid` (2001 + 4001 points), `vortex_stress_force` × 3 (sympy-free closed forms), dissipation `quad`, slider 21 steps × 4 traces × 300 points, two-panel figure | 5 s | 4 s (11 steps) |
| §5.2 C03 | cellular loop: 256 points (with insertion) to t = 37.7 s (DOP853, rtol 1e-10, point insertion) ≈ 2 s, `kelvin_rate_terms`, 16-row hypotheses table, from-scratch RK4 (256 points × 600 steps, vectorised), rotating scenario, animation 60 frames (loop + Γ(t)) | 16 s | 9 s (30 frames, 128 points, RK4 300 steps) |
| §5.2 C04 | element closed forms, three radii for the O(R²) order, from-scratch disc (4000 rim + 400 × 800 area cells), two-disc figure, slider 19 steps × 2 traces, lock-exchange map (`baroclinic_term` on 201 × 101) | 5 s | 4 s (10 steps, 101 × 51) |
| §5.3 C05 | ABC loop (256 points, 2 s), vector areas, plotly 3-D tube + two patches | 4 s | 3 s |
| §5.4 C06 | **D09 sympy ★★★** (curl of each NS term for generic P, Q, R; (B.3.10) check; residual on a polynomial field) ≈ 3 s, term presets (nested stencils at 3 points), from-scratch Laplacian, diffusing-sheet slider 20 steps × 2 traces × 400 points, FTCS cross-check (2001 points to t = 1 s: ≈ 0.3 s), Hill residual, frozen-in (ABC 5 s, Burgers 2 s) and its figure | 12 s | 9 s (10 steps, 1001-point FTCS) |
| §5.5 C07 | FFT Poisson (32² and 64² periodic), Green primer sympy, **D10 sympy ★★★** (∇²(1/r), sphere flux) ≈ 1 s, sign test (40³ = 64 000 Gauss nodes × 2 signs + Biot–Savart, one field point), **D11 sympy ★★★** (product rule, ∇′(1/r), step 5) ≈ 2 s, curl theorem boxes (n = 6 exact, n = 24), planar Biot–Savart 201² × 3 points, from-scratch double loop (40 401 cells, one point) ≈ 0.4 s, slider over the grid (10 steps: grids 21² … 201², 60 radii each) ≈ 2 s, Fig. 5.8 3-D figure | 14 s | 8 s (24³ nodes, loop over every 2nd cell, 6 slider steps up to 121²) |
| §5.5 C08 | segments and polygons (closed forms), `quad` limit primer, from-scratch 20 000-piece sum, M-slider 10 steps × 2 traces, speed-vs-length figure | 3 s | 3 s |
| §5.6 C09 | rotating residuals, `baroclinic_term_sym`, **D15 sympy ★★★** (LeviCivita sums over 5 indices for 3 components with generic u, then a polynomial divergence-free field; `expand` before `simplify`, one component at a time) ≈ 6 s, `vorticity_budget_preset` × 3, `vorticity_budget_sym` on a polynomial field ≈ 2 s, term-bar figure | 14 s | 12 s |
| §5.6 C10 | split, `expm` × 3, plotly 3-D helix with frame (3 dropdown states), frames player 24 frames (tube + segment + growth panel), Burgers profile + term curves (401 points, closed forms) | 8 s | 5 s (12 frames) |
| §5.6 C11 | column and ring closed forms, rotating-scenario loop (256 points, 5 s), vector areas, column-over-ridge figure, latitude slider 19 steps × 2 traces | 5 s | 4 s (10 steps) |
| §5.7 C12 | pair closed forms, one-period integration (200 outputs), invariants, from-scratch double loop (3 vortices), animation 60 frames × 3 panels, trajectories + invariants figure (4 presets × 20 s) | 12 s | 7 s (30 frames) |
| §5.7 C13 | wall images (200 wall points), channel series (until the increment < 1e-12 ≈ 10⁵ terms, vectorised in blocks), knife-in-bucket animation 60 frames (circle images), ring self/ring–ring (elliptic), ring–wall `ring_dynamics` (301 outputs), leap-frog `ring_dynamics` (401 outputs) + animation 60 frames | 22 s | 12 s (30 frames, 151/201 outputs) |
| §5.8 C14 | sheet sums (N ≤ 1000), `discrete_sheet_convergence` (3 × `quad` with 2000 subintervals or a 4001-point trapezoid) ≈ 1 s, N-slider 12 steps × 2 traces × 400 points, **roll-up** N = 200 to t = 2 s (DOP853, vectorised O(N²) kernel) ≈ 6 s + animation 60 frames | 16 s | 7 s (N = 100, 30 frames) |
| explainers (9 × `show_viz`) | read the HTML files | 1 s | 1 s |
| **Total** | | **≈ 153 s** | **≈ 102 s** |

**FAST plan.** Every size-dependent choice is written `a if not FAST else b` in the cell: animations 60 frames → 30,
frames players 24 → 12; slider figures ≤ 21 steps → ≤ 11 (≤ 4 traces × ≤ 400 points, each < 250 kB); 3-D Gauss
quadrature 40³ → 24³; material loops 256 → 128 points; roll-up N = 200 → 100; lid sums 200 × 400 → 100 × 200; the
from-scratch Biot–Savart double loop runs over every second cell (4 × dA); ring runs 301/401 → 151/201 outputs; plotly
3-D surfaces 40² → 24². **Cached arrays / results:** `ch05.kelvin_scenario("cellular")`'s advected loop (C03) is reused
by the C03 animation and the E3 notebook cell; the Gaussian-tube fields and Gauss nodes of C07 are built once and reused
by N27–N28 and the sign test; the Burgers field `vortices.burgers_vortex_field(1e-3, 1.0, 1e-6)` is built in C01 (N45)
and reused in C06, C10; `ch05.abc_flow()` is built in C05 and reused in C06 (frozen-in) and C09 (N30); the point-vortex
runs of C12 are reused by its figure; the library caches the sympy results of D09, D15 (`vorticity_equation_sym`,
`vorticity_budget_sym`) and the ring-ring elliptic evaluations with `functools.lru_cache`, so a second call is free.
**Outputs:** 6 videos (C03 loop, C12 pairs, C13 bucket and leap-frog, C14 roll-up, plus the C10 frames player), 6 slider
figures, 4 plotly 3-D figures (C01 tubes, C05 patch, C10 helix, and the C07 geometry as matplotlib 3-D), ≈ 20 static
figures — page well under 15 MB (each video < 2 MB at dpi 80). Sympy cells never `simplify` expressions with more than
three generic functions at once: D15's check expands and simplifies one component and one term at a time.

---

## Part E — prerequisite ledger
Every concept, symbol, maths tool and Python function or idiom the notebook or its explainers use, with where it is
explained. "primer (in Cxx)" = a 📎 primer placed in that block before first use (the primer term is the Concept text —
the builder uses it verbatim in `nb.primer`); "knowledge/primers.md: <term> (chNN Pnn)" = a one-line reminder naming the
earlier primer; a CORE/RECAP id alone = taught there; "Cxx (Nnn note)" = the NOTE of this chapter placed in that block;
"Cxx (gloss …)" = one sentence where it is used. Earlier chapters' material that is neither a primer nor a ch05 RECAP is
named by section ("Ch. 3 §3.5"), never by an earlier chapter's C-number (convention 11). Explainer-only items are listed
with the explainer's host block. New primers P134–P148 in first-use order: partial integration (C02), closed-loop exact
differential, periodic trapezoid, many points in one solve_ivp (C03), loop vector area (C05), Poisson/Green, gradient of
1/distance, curl of a product, FFT Poisson, 3-D Gauss–Legendre (C07), improper integral (C08), Frenet frame, spinning
cylinder (C10), ODE systems with invariants (C12), elliptic integrals (C13) — 15 primers.

| Concept | First used in | Explained by |
|---|---|---|
| vorticity ω = ∇×u and particle spin ½ω | C01 | R01 |
| vortex and vortex motion (words) | C01 | R01 |
| the chapter's programme (carried, stretched, tilted, diffused) | C01 | C01 (N01 note) |
| kinematics.vorticity (stencil curl of a callable) | C01 (R01 code) | R01 (code comment; Ch. 3 §3.4) |
| functions as arguments and lambda | C01 (R01 code) | knowledge/primers.md: functions as arguments and lambda (ch01 P29) — reminder |
| numpy arrays and vectorised arithmetic | C01 | knowledge/primers.md: numpy arrays (ch01 P03) — reminder |
| f-strings | C01 | knowledge/primers.md: f-strings (ch01 P04) — reminder |
| vortex line (5.3) | C01 | C01 (N02 note) |
| streamline (3.7) as the analogue of a vortex line | C01 | C01 (N02 note; Ch. 3 §3.3) |
| parametric curves, tangent vector and arc length | C01 | knowledge/primers.md: parametric curves, tangent vector and arc length (ch03 P92) — reminder |
| parallel vectors and the cross-product test | C01 | knowledge/primers.md: parallel vectors and the cross-product test (ch03 P93) — reminder |
| helical swirl test field u_φ = aRz (vortex lines zR² = const) | C01 | C01 (N02 note and code explain) |
| vorticity_field and vortex_line (library calls) | C01 | C01 (code explain) |
| np.ptp and np.hypot | C01 | C01 (gloss in code comments: "peak-to-peak spread", "√(x² + y²)") |
| vortex tube and tube strength | C01 | C01 (N03 note) |
| stream tube and volume flux dQ = u·n dA | C01 | C01 (N03 table; Ch. 3 §3.3) |
| circulation Γ = ∮u·dx (3.18) | C01 | C01 (N03 note; Ch. 3 §3.4) |
| Stokes' theorem (2.34) | C01 | C01 (N03 note; Ch. 2 §2.13) |
| Gaussian (Lamb–Oseen) vortex (3.29) | C01 | C01 (N03 note; Ch. 3 §3.5) |
| planar_loop, planar_disc, circulation, curl_flux | C01 | C01 (code explain; Ch. 2 §2.13) |
| right-hand rule and orientation of a loop | C01 | knowledge/primers.md: right-hand rule and orientation (ch02 P74) — reminder |
| divergence of a curl is zero, ∇·(∇×u) = 0 | C01 (D01) | C01 (D01 start line; Ch. 2 §2.13; recapped as R10 in §5.4) |
| Gauss' divergence theorem (2.30) | C01 (D01) | C01 (D01 tools; Ch. 2 §2.12) |
| outward normal on a closed surface | C01 (D01) | C01 (D01 step 6) |
| tubes cannot end (5.4) | C01 | C01 |
| mean vorticity ω̄ = Γ/A | C01 | C01 (D01 step 7, worked example) |
| narrowing Gaussian tube and its flux function | C01 | C01 (code explain) |
| field that is not solenoidal (broken tube) | C01 | C01 (code explain) |
| named results (NamedTuple TubeStrength, TubeFlux) | C01 | knowledge/primers.md: dataclasses and named results (ch04 P111) — reminder |
| volume and surface integrals as midpoint sums (polar grid) | C01 | knowledge/primers.md: volume and surface integrals as midpoint sums (ch02 P83) — reminder |
| assert np.allclose | C01 | knowledge/primers.md: assert np.allclose (ch01 P15) — reminder |
| plotly 3-D lines and surfaces | C01 | knowledge/primers.md: plotly 3-D arrows, lines and meshes (ch02 P64) — reminder |
| Burgers vortex as a test field (before C10) | C01 (N45 figure) | C01 (gloss in the N45 cell: "a vortex fed by axial stretching, derived in C10") |
| log axes | C01 | knowledge/primers.md: power laws and log–log plots (ch01 P13) — reminder |
| matplotlib figures and labels with units | C01 | knowledge/primers.md: matplotlib figures (ch01 P01) — reminder |
| show_viz and embedded explainers | C01 | knowledge/primers.md: show_viz (ch01 P18) — reminder |
| solid-body rotation from vorticity (5.1), rotation rate ω/2 | C02 | R02 |
| ω is the vorticity, not the rotation rate (Fig. 5.2's "2ω") | C02 | R02 (⚠️ callout) |
| ideal line vortex (5.2), Γ = 2πB | C02 | R03 |
| paddle-wheel test | C02 | R04 |
| no deformation ⇒ no viscous stress ⇒ Euler | C02 | R05 |
| Newtonian stress (4.37) and Euler's equation (4.41) | C02 | R05 |
| Cauchy's equation (4.24) | C02 (D02) | R05 (and D02 start; Ch. 4 §4.4) |
| hydrostatic balance (5.5b) | C02 | R06 |
| Bernoulli function (4.69)–(4.72) across streamlines | C02 | R07 |
| Rankine vortex (3.28) and the rotating cylinder | C02 | R08 |
| centripetal acceleration u_θ²/r | C02 (D02) | C02 (D02 steps 3–4) |
| cylindrical and spherical unit vectors | C02 (D02) | knowledge/primers.md: cylindrical and spherical unit vectors (ch03 P88) — reminder |
| derivative of a rotating unit vector ∂e_θ/∂θ = −e_r | C02 (D02) | knowledge/primers.md: derivative of a rotating unit vector (ch04 P124) — reminder |
| partial integration with an unknown function | C02 (D02) | primer (in C02) |
| radial pressure balance (5.5a) | C02 | C02 (N05 note) |
| paraboloid isobars and the rotating-tank pressure (5.6) | C02 | C02 (D02) |
| volume-conserving free surface of a tank | C02 | C02 (code explain) |
| scipy.optimize.brentq | C02 | knowledge/primers.md: scipy.optimize.brentq (ch03 P108) — reminder |
| sympy symbols, diff, simplify, Function | C02 | knowledge/primers.md: sympy (ch01 P40) — reminder |
| trapezoid rule and cumulative_trapezoid | C02 | knowledge/primers.md: trapezoid rule (ch01 P37) — reminder |
| polar strain rate S_rθ and polar stress divergence | C02 (D03) | C02 (gloss in D03 steps 1 and 6, sympy check) |
| viscous stress σ = 2μS (4.59) | C02 (D03) | C02 (D03 start; Ch. 4 §4.5) |
| viscous force μ∇²u = −μ∇×ω (4.40) | C02 (D03) | C02 (D03 step 9; Ch. 4 §4.6) |
| vector Laplacian in polar form | C02 (D03) | C02 (D03 step 8; Appendix B via `core.curvilinear`) |
| σ means viscous stress here, σ_c the core radius | C02 | C02 (⚠️ in N06) |
| line-vortex pressure (5.7) and the funnel isobars | C02 | C02 (N07 note) |
| isobars are cubic, not "hyperboloids of the second degree" | C02 | C02 (⚠️ in N07) |
| torque per unit length and its power | C02 | C02 (N08 note) |
| torque, moment arm and moment of inertia | C02 | knowledge/primers.md: torque, moment arm and moment of inertia (ch04 P116) — reminder |
| power of a force | C02 | knowledge/primers.md: power of a force and heat flux through a surface (ch04 P127) — reminder |
| dissipation rate ε (4.58); ε ≠ ε_ijk | C02 | C02 (N08 note, ⚠️) |
| scipy.integrate.quad and dblquad | C02 | knowledge/primers.md: scipy.integrate.quad and dblquad (ch03 P87) — reminder |
| irrotational = no net viscous force, not no stress | C02 | C02 (N09 note) |
| Lamb–Oseen vortex with σ_c² = 4ν(t + t₀) | C02 (N09) | C02 (N09 table; Ch. 3 §3.5) |
| cyclostrophic balance, tornado pressure deficit | C02 | C02 (worked example) |
| slider_figure | C02 | knowledge/primers.md: slider_figure (ch01 P17) — reminder |
| fluid properties of water and air | C02 (E2) | C02 (explainer note; Ch. 1 §1.3) |
| barotropic fluid and the pressure function (4.67) | C03 | R09 |
| conservative body force and potential Φ (4.18) | C03 | knowledge/primers.md: conservative force and its potential (ch04 P115) — reminder |
| material loop | C03 | C03 |
| material label (particle tags) | C03 (D04) | C03 (gloss before D04) |
| functions of time with parameters | C03 (D04) | knowledge/primers.md: functions of time with parameters (ch03 P89) — reminder |
| differentiation under the integral sign | C03 (D04) | knowledge/primers.md: differentiation under the integral sign (ch03 P109) — reminder |
| material derivative D/Dt (3.5) | C03 (D04) | C03 (D04 step 3; Ch. 3 §3.2) |
| material line element D(δx)/Dt = δu | C03 (D04) | knowledge/primers.md: material line element (ch03 P99) — reminder |
| multivariable first-order Taylor expansion | C03 (D04) | knowledge/primers.md: multivariable first-order Taylor expansion (ch03 P98) — reminder |
| chain rule for the kinetic energy u·du = d(½u²) | C03 (D04) | knowledge/primers.md: chain rule for the kinetic energy (ch04 P128) — reminder |
| closed-loop integral of an exact differential | C03 | primer (in C03) |
| np.arctan2 and np.unwrap | C03 (primer code) | knowledge/primers.md: np.arctan2 (ch02 P70) — reminder (np.unwrap glossed in the primer's comment) |
| multivariable chain rule along a path | C03 (D05) | knowledge/primers.md: multivariable chain rule along a path (ch03 P91) — reminder |
| periodic trapezoid rule on a closed loop | C03 | primer (in C03) |
| FFT derivative of a periodic sequence (dx/ds) | C03 (code) | C03 (gloss inside the periodic-trapezoid 📎 entry; C07 develops it) |
| advecting many points in one solve_ivp call | C03 | primer (in C03) |
| solve_ivp options (rtol, t_eval, DOP853) | C03 | knowledge/primers.md: solve_ivp options: t_eval, dense_output, events, backward integration (ch03 P94) — reminder |
| RK4 by hand | C03 | knowledge/primers.md: RK4 by hand (ch03 P95) — reminder |
| 2-D stream function ψ (cellular flow) | C03 | C03 (worked example; Ch. 4 §4.3) |
| Kelvin's circulation theorem (5.8) | C03 | C03 |
| rate of a material loop's circulation (5.9) | C03 | C03 (N10 note) |
| D(dx)/Dt = du and ∮u·du = 0 | C03 | C03 (N11 note) |
| forces round the loop (5.10) | C03 | C03 (N12 note) |
| viscous leak (5.11) and its incompressible-only remark | C03 | C03 (N13 note) |
| three sources of vorticity | C03 | C03 (N14 note) |
| four restrictions (decision table) | C03 | C03 (N16 note) |
| itertools.product | C03 | C03 (gloss in the N16 code comment; Ch. 2 used it) |
| Coriolis force is not conservative | C03 | C03 (N14 note; Ch. 4 §4.7) |
| rotating frame of reference | C03 | knowledge/primers.md: rotating frame of reference (ch03 P103) — reminder |
| animate and show_animation | C03 | knowledge/primers.md: animate and show_animation (ch01 P16) — reminder |
| barotropic vs baroclinic element (Fig. 5.6) | C04 | C04 |
| centre of mass of a body with varying density | C04 (D06) | C04 (gloss before D06) |
| moment of inertia of a disc I = ½MR² | C04 (D06) | C04 (gloss before D06) |
| ∫x_i x_j dA = (πR⁴/4)δ_ij over a disc | C04 (D06) | C04 (D06 step 5) |
| net force from pressure and the gradient theorem | C04 (D06) | knowledge/primers.md: net force from pressure (ch01 P28) — reminder |
| vorticity = twice the spin (factor 2) | C04 (D06) | R01 |
| baroclinic spin-up ∇ρ×∇p/ρ² | C04 | C04 (D06) |
| smoothed step with np.tanh | C04 | C04 (gloss before D06) |
| lock exchange and its initial rate | C04 | C04 (N15 note, D07) |
| reduced gravity | C04 | knowledge/primers.md: reduced gravity and buoyancy (ch04 P131) — reminder |
| observed order of convergence (log–log slope) | C04 | knowledge/primers.md: power laws and log–log plots (ch01 P13) — reminder |
| np.cross on arrays of vectors | C04 | C04 (gloss in code comment; Ch. 4 used it) |
| Helmholtz's four vortex theorems | C05 | C05 |
| small-patch localisation ("for every S") | C05 (D08) | knowledge/primers.md: continuity of a function and the small-ball argument (ch04 P112) — reminder |
| ABC (Beltrami) flow, ω = u | C05 | C05 (worked example) |
| loop vector area ½∮x × dx | C05 | primer (in C05) |
| np.roll (next point round a loop) | C05 | C05 (gloss in the loop-area demo's code comment) |
| Cauchy's frozen-in solution (ω and δx obey one equation) | C06 | C06 (md cell after N23) |
| divergence of a curl, recap | C06 | R10 |
| Lamb identity (4.68) | C06 | R11 |
| curl of a gradient is zero | C06 (D09) | C06 (D09 step 2; Ch. 2 §2.13) |
| Schwarz's theorem | C06 (D09) | knowledge/primers.md: Schwarz's theorem (ch04 P121) — reminder |
| identity (B.3.10) for ∇×(ω×u) | C06 | C06 (N20 note) |
| ε–δ identity (2.19) | C06 | C06 (N20 note; Ch. 2 §2.7) |
| vorticity equation (5.13) | C06 | C06 |
| curl of Navier–Stokes (5.12) | C06 | C06 (N19 note) |
| §5.4 hypotheses | C06 | C06 (N18 note) |
| velocity gradient tensor G and (ω·∇)u = Gω | C06 | C06 (from-scratch cell; Ch. 3 §3.4) |
| directional derivative (ω·∇) | C06 | knowledge/primers.md: level sets and the directional derivative (ch02 P75) — reminder |
| Cartesian vector Laplacian commuting with the curl | C06 (D09) | C06 (D09 step 7) |
| one-dimensional diffusion equation | C06 (N22) | C06 (N22 note; Ch. 1 §1.5) |
| error function erf = 1 − erfc | C06 (N22) | C06 (gloss in N22) |
| complementary error function erfc | C06 (N22) | knowledge/primers.md: complementary error function erfc (ch04 P123) — reminder |
| finite differences and FTCS | C06 | knowledge/primers.md: finite differences (ch01 P21) — reminder |
| Hill's spherical vortex | C06 | C06 (N23 note) |
| axisymmetric (Stokes) stream function | C06 (N23) | C06 (N23 note; Ch. 4 §4.3) |
| Taylor–Green vortex | C06 | C06 (code explain; Ch. 4 §4.6) |
| vorticity_terms and vorticity_budget_preset (library calls) | C06 | C06 (code explain) |
| Poisson equation and Green's function | C07 | primer (in C07) |
| Dirac delta as a point source | C07 | C07 (the Poisson-equation 📎 entry; Ch. 3 §3.5 gloss) |
| gradient of 1/distance with respect to the source point | C07 | primer (in C07) |
| curl of a product (product rule) | C07 | primer (in C07) |
| cross-product antisymmetry a × b = −b × a | C07 (D11) | C07 (gloss in the curl-of-a-product 📎 entry) |
| curl of a curl identity | C07 (N24) | knowledge/primers.md: curl of a curl identity (ch04 P122) — reminder |
| Fourier modes and the FFT Poisson solver | C07 | primer (in C07) |
| np.meshgrid and the project grid layout | C07 | knowledge/primers.md: np.meshgrid and the project grid layout (ch02 P76) — reminder |
| Gauss–Legendre quadrature in 3-D and a smoothed kernel | C07 | primer (in C07) |
| superposition for a linear equation | C07 (D10) | C07 (the Poisson-equation 📎 entry; D10 step 9) |
| corrected (5.14): +1/(4π) | C07 | C07 (N25 note, D10) |
| integrand rewrite and the second sign slip | C07 | C07 (N26 note, D11) |
| Gauss' theorem in curl form (5.15) | C07 | C07 (N27 note) |
| choice of the volume V′ (Fig. 5.8) | C07 | C07 (N28 note) |
| Biot–Savart law (5.16) | C07 | C07 |
| planar Biot–Savart kernel e_z × r/(2π r²) | C07 | C07 (code explain) |
| primes mark the source point x′ | C07 | C07 (notation cell ⚠️; D10 step 10) |
| outside vorticity adds nothing (planar shell theorem) | C07 | C07 (worked example) |
| matplotlib 3-D axes | C07 (N52 figure) | C07 (gloss in the figure code comment) |
| np.ndindex (loop over grid cells) | C07 | C07 (gloss in the from-scratch code comment) |
| filament law (5.17) | C08 | C08 |
| far-field ("frozen kernel") approximation | C08 | C08 (gloss) |
| trig identity 1 + cot²θ = csc²θ | C08 (D13) | C08 (gloss) |
| substitution in an integral | C08 (D13) | knowledge/primers.md: substitution in an integral (ch03 P106) — reminder |
| improper integral as a limit | C08 | primer (in C08) |
| straight segment (Γ/4πd)(cos θ_a − cos θ_b); infinite line Γ/2πd | C08 | C08 (N29 note, D13) |
| vortex ring and its axis velocity ΓR²/2(R² + z²)^{3/2} | C08 | C08 (N29 note) |
| polygon (M-gon) convergence ∝ 1/M² | C08 | C08 (code explain) |
| comma notation u_{i,j} and renaming a free index | C09 | R12 |
| Boussinesq (nearly incompressible, ρ kept with gravity) | C09 | R12 (and the ⚠️ after N37; Ch. 4 §4.9) |
| continuity (5.19) | C09 | R13 |
| rotating-frame momentum (5.20) and effective gravity | C09 | R14 |
| index Lamb identity (5.21) | C09 | R15 |
| antisymmetric gradient = vorticity (5.22) | C09 | R16 |
| np.einsum index strings | C09 (R16 code) | knowledge/primers.md: np.einsum index strings (ch02 P62) — reminder |
| np.random.default_rng | C09 (R16 code) | knowledge/primers.md: np.random.default_rng (ch01 P10) — reminder |
| viscous term as −ν∇×ω (5.23) | C09 | R17 |
| planetary vorticity 2Ω (first mention) | C09 | C09 (gloss in the idea cell; recapped in full as R18 before C11) |
| ∇·ω = 0 in index form (5.18) | C09 | C09 (N30 note) |
| symmetric × antisymmetric contraction = 0 | C09 | C09 (N30 note; Ch. 2 §2.7) |
| quotient rule for a gradient ∇(1/ρ) = −∇ρ/ρ² | C09 | C09 (gloss) |
| dummy-index relabelling | C09 (D14) | C09 (D14 step 6; Ch. 2 §2.1) |
| Coriolis term reordered (5.24) | C09 | C09 (N31 note) |
| rotating Navier–Stokes in Lamb form (5.25) | C09 | C09 (N32 note, D14) |
| inertial oscillation as a test field | C09 (N32 code) | C09 (gloss in N32's code comment) |
| permutations, cyclic order and parity of ε | C09 (D15) | knowledge/primers.md: permutations, cyclic order and parity (ch02 P72) — reminder |
| product rule for differentials | C09 (D15) | knowledge/primers.md: product rule for differentials (ch01 P38) — reminder |
| curl of (5.25) in index form (5.26); Π is Φ | C09 | C09 (N33 note) |
| the dropped term u_{j,j}(ω_n + 2Ω_n) (5.27) | C09 | C09 (N34 note, D15 step 9) |
| baroclinic vector (5.28) | C09 | C09 (N35 note) |
| viscous term → ν∇²ω (5.29) | C09 | C09 (N36 note) |
| vorticity equation in a rotating frame (5.30) | C09 | C09 |
| reading the terms of (5.30) | C09 | C09 (N37 note and figure) |
| symmetric-log axis | C09 (figure) | C09 (gloss in the figure code comment) |
| Frenet frame of a curve | C10 | primer (in C10) |
| e_n away from the centre of curvature (book convention) | C10 | C10 (the Frenet-frame 📎 entry) |
| angular momentum of a spinning cylinder | C10 | primer (in C10) |
| natural coordinates (5.31) | C10 | C10 (N38 note, D16) |
| stretching and tilting (5.32) | C10 | C10 |
| orthonormal basis, projection and completeness | C10 (D16) | knowledge/primers.md: orthonormal basis, projection and completeness (ch02 P65) — reminder |
| scipy.linalg.expm and e^{Gt}ω₀ | C10 | knowledge/primers.md: scipy.linalg.expm (ch02 P79) — reminder |
| strain presets (axial, shear, planar) | C10 | C10 (code explain) |
| Burgers' vortex | C10 | C10 (N21 note, D18) |
| cylindrical Laplacian of an axisymmetric scalar | C10 (D18) | C10 (gloss in D18 step 5) |
| incompressible tube: A L = const | C10 (D18) | C10 (gloss in D18 step 2) |
| separation of variables | C10 (D18) | knowledge/primers.md: separation of variables (ch01 P42) — reminder |
| natural logarithm and exponential | C10 | knowledge/primers.md: natural logarithm and exponential (ch01 P36) — reminder |
| plotly dropdown menus | C10 | C10 (gloss in the plotly code comment) |
| planetary and absolute vorticity; f = 2Ω sin φ | C11 | R18 |
| latitude, Earth's rotation rate and the local vertical | C11 | knowledge/primers.md: latitude, Earth's rotation rate and the local vertical (ch04 P126) — reminder |
| np.deg2rad | C11 | knowledge/primers.md: cosines of angles between unit vectors (ch02 P66) — reminder |
| absolute circulation (5.33) | C11 | C11 |
| product rule for a cross product | C11 (D19) | knowledge/primers.md: product rule for a cross product (ch04 P125) — reminder |
| cyclic triple product (Ω × u)·dx = Ω·(u × dx) | C11 (D19) | C11 (D19 step 4; Ch. 2 §2.7) |
| column mass: A h = const | C11 (D20) | C11 (gloss before D20) |
| fluid column (ζ + f)/h = const | C11 | C11 (D20) |
| planetary stretching and tilting | C11 | C11 (N39 note) |
| cyclonic vs anticyclonic | C11 | C11 (worked example) |
| Taylor–Proudman and potential vorticity (named, Ch. 13) | C11 | C11 (D20 assumptions; What would change if…) |
| point vortex (line-vortex idealisation) and superposition | C12 | C12 (gloss) |
| lever rule | C12 (D21) | C12 (gloss) |
| circular motion: speed = rate × radius | C12 (D21) | C12 (gloss) |
| systems of ODEs for interacting bodies and their invariants | C12 | primer (in C12) |
| centre of vorticity and the orbit rate | C12 | C12 (D21, N40 note) |
| equal and opposite pair translating | C12 | C12 (N41 note) |
| Fig. 5.11 caption slip (G is not a point of zero velocity) | C12 | C12 (⚠️ in N40, D21 check) |
| Kirchhoff energy H with a logarithm | C12 | C12 (the ODE-systems 📎 entry) |
| method of images and mirror symmetry | C13 | C13 (gloss) |
| inverse point a²x/r² (circle images) | C13 | C13 (N57 note) |
| channel image series and its cot form | C13 | C13 (worked example) |
| complete elliptic integrals with the parameter m | C13 | primer (in C13) |
| thin-core ring self-speed (Kelvin's formula) | C13 | C13 (N42 note) |
| ring impulse ΣΓπR² | C13 | C13 (N43 note) |
| vortex sheet: tangential jump, continuous normal velocity | C14 | C14 (N44 note) |
| sheet strength γ = u₂ − u₁ and the two conventions | C14 | C14 (D23, ⚠️) |
| Γ (circulation, m²/s) vs γ (sheet strength, m/s) | C14 | C14 (⚠️; notation cell) |
| Riemann sum turning into an integral | C14 | C14 (gloss) |
| L1 error between two curves | C14 | C14 (gloss in code explain) |
| Krasny smoothing and the periodic sheet kernel | C14 | C14 (gloss in the roll-up animation cell) |
| Kelvin–Helmholtz roll-up (named, Ch. 11) | C14 | C14 (animation read) |
| linear-growth ghost e^{γkt/2} (explainer only) | C14 (E9) | C14 (gloss in E9's growth view: "preview of Ch. 11") |
| rotating-frame symbols Ω [rad/s], f [1/s], ζ relative vorticity [1/s] | C09 | C09 (notation cell; R18 before C11) |
| h column height [m], δ interface thickness [m], α strain rate [1/s] | C04 | C04 (worked example; C10, C11 where each appears) |
| Φ force potential [J/kg] | C03 | C03 (D05 step 4; ch04 P115 reminder) |
| A_vec vector area [m²], Γ_a absolute circulation [m²/s] | C05 | C05 (loop vector area), C11 |
| p_o (pressure at r = 0, z = 0) vs p_∞ (far from a line vortex) | C02 | C02 (D02 step 9, N07 note) |

---

## Part F — derivation storyboards

Builders copy these word for word into `nb.derivation(key, title, goal=…, start=(tex, plain), plan=[…], uses=[…],
steps=[dict(did, tex, why, plain)], result=(tex, plain), interpret=…, check=…, check_src=…)` and into the explainer's
`derivations: [...]` (phones may shorten *why* to its first sentence; `live`, `set` and `watch` are the explainer's).
Every step is one move; *why* names the rule and says why we make the move (≤ 35 words, ≥ 6); *did* ≤ 8 words. The
book's own moves were read on the rendered pages (p200 for D01, p200–p201 for D02, p201–p202 for D03 (the book states
only σ_rθ and leaves the force to Exercise 5.4), p203–p204 for D04–D05, p205–p206 for D06–D07 (the book gives only the
figure and, for D07, an exercise), p206–p207 for D08, p207–p208 for D09, p208 for D10 (the book gives only the result,
with the wrong sign, and sends the reader to "Exercise 5.8", meaning 5.9), p209 for D11–D12, none for D13, p210–p212
for D14–D15, p213 for D16–D17, none for D18, p214 for D19–D20 (result only; "Exercise 5.11" means 5.10), p215 for D21,
p215–p216 for D22, p218 for D23); the gaps listed in `analysis/ch05.md` §2b are filled and the notebook says so ("the
book skips this move"). Every equation named by number is written out. Colours: lower end blue, side grey, upper end
orange (D01); centripetal teal, pressure orange, viscous rose (D02, D03); relative vorticity teal, planetary amber,
baroclinic orange, diffusion rose, stretching purple, tilting blue (D09, D15, D16–D18); positive vortices teal, negative
rose (D21–D23). (No line in this part starts with a table bar; absolute values are written with \lvert \rvert or in
words.)

### D01 · A vortex tube has the same strength everywhere, so it cannot end: (5.4) — ★, 7 steps, in C01 (notebook · `vortex_tubes_cannot_end`)
- **Goal.** Show that the circulation measured round a vortex tube does not depend on where along the tube we measure
  it — and conclude that a tube (or a line) cannot simply stop in the middle of the fluid.
- **Start.** $\nabla\cdot\boldsymbol\omega=\nabla\cdot(\nabla\times\mathbf u)=0$ — *in words:* vorticity, being a curl, has no
  sources or sinks anywhere (R10, Ch. 2 §2.13).
- **Plan.** (1) Cut a finite piece out of a tube. (2) Apply Gauss' theorem to ω on that piece. (3) Split the surface into
  its two ends and its side. (4) Read the end fluxes as tube strengths, minding the outward normals.
- **Tools.** Gauss' theorem $\int_V\nabla\cdot\mathbf Q\,dV=\int_A\mathbf Q\cdot\mathbf n\,dA$ (2.30) (Ch. 2 §2.12) · tube strength
  $\Gamma=\oint_C\mathbf u\cdot d\mathbf x=\int_A\boldsymbol\omega\cdot\mathbf n\,dA$ (N03, Stokes (2.34)) · right-hand rule and
  orientation (ch02 P74, reminder).
- **Assumptions.** u twice continuously differentiable so ∇·ω = 0 (step 3); the tube wall is made of vortex lines
  (step 5). Nothing about viscosity, density or time: the result is kinematic.
- **Steps.**
  1. *did:* Cut a piece out of a vortex tube · *tex:* $\partial V=A_{\text{lower}}\cup A_{\text{side}}\cup A_{\text{upper}}$ · *why:*
     Any finite stretch of tube is a closed volume V bounded by two cross-sections and the tube wall; we pick one so
     Gauss' theorem has a closed surface to work on. · *plain:* A short length of tube is a sealed can: two lids and a
     wall. · *set:* {mode: 'narrow', zs: 0.5} · *watch:* "the blue lid, grey wall and orange lid".
  2. *did:* Apply Gauss' theorem to ω · *tex:* $\int_V\nabla\cdot\boldsymbol\omega\,dV=\int_{\partial V}\boldsymbol\omega\cdot\mathbf n\,dA$ ·
     *why:* Divergence theorem (2.30) with Q = ω, n the outward normal; allowed because ω is smooth. We want to turn a
     statement about the inside into one about the surface. · *plain:* What the inside sources equals the net vorticity
     flux out through the skin.
  3. *did:* Use that ω has no divergence · *tex:* $\int_{\partial V}\boldsymbol\omega\cdot\mathbf n\,dA=0$ · *why:* ∇·(∇×u) = 0 for
     any smooth u (ch02, R10), so the volume integral is zero. The book states this after (5.4); we use it first. ·
     *plain:* As much vorticity flux leaves the can as enters it.
  4. *did:* Split the skin into three pieces · *tex:* $\Big\{\int_{\text{lower}}+\int_{\text{side}}+\int_{\text{upper}}\Big\}\boldsymbol\omega\cdot\mathbf n\,dA=0$ ·
     *why:* An integral over a surface made of pieces is the sum over the pieces (additivity). This is the first line of
     (5.4). · *plain:* Lower lid plus wall plus upper lid add to zero.
  5. *did:* Drop the wall: no flux through it · *tex:* $\int_{\text{side}}\boldsymbol\omega\cdot\mathbf n\,dA=0$ · *why:* The wall is
     built from vortex lines, so ω is tangent to it and ω·n = 0 there — not because ω vanishes on the wall. ·
     *plain:* Vorticity runs along the wall, never through it. · *set:* {mode: 'twisted'} · *watch:* "the grey wall bar
     stays at zero even when the lines twist".
  6. *did:* Read each lid as a signed strength · *tex:* $-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0$ · *why:* On the
     upper lid the outward n points along ω, giving +Γ; on the lower lid the outward n points against ω, giving
     −Γ (orientation, P74). This is (5.4). · *plain:* The strength entering at the bottom equals the strength leaving at
     the top. · *live:* "−(0.632) + 0.632 = 0 m²/s for your tube".
  7. *did:* Shrink a section to see "cannot end" · *tex:* $\bar\omega=\Gamma/A\to\infty\ \text{as}\ A\to0$ · *why:* By step 6 every
     section carries the same Γ ≠ 0; Γ/A would blow up if a section closed to a point, so the tube must end on a wall
     or a free surface, or close on itself. · *plain:* A tube can thin and spin faster, but never
     just stop. · *set:* {zs: 1.0} · *watch:* "mean ω rises as the area falls; the flux bar does not move".
- **Result.** $\int_V\nabla\cdot\boldsymbol\omega\,dV=\int_A\boldsymbol\omega\cdot\mathbf n\,dA=-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0$
  (5.4) — *in words:* a vortex tube has one strength along its whole length; tubes and lines end only on boundaries or
  close as loops.
- **Check.** Units: m²/s on both lids ✓. Numbers (narrowing Gaussian tube, Γ = 1 m²/s, a₀ = R₀ = 0.1 m, L = 1 m): flux
  0.632121 at z = 0 and at z = 1 m while the area falls from 0.031416 to 0.011557 m² and the mean vorticity rises from
  20.12 to 54.69 s⁻¹ ✓. Broken field (∇·ω ≠ 0): lower −0.632, upper +0.233, total −0.400 ≠ 0 — the theorem fails exactly
  where its hypothesis does ✓.
- **What it means.** Helmholtz's second and third theorems are this one line; they hold in viscous flow too, because
  nothing dynamical was used. Squeeze a tube and its vorticity must grow (C10); a wing's bound vortex must continue as
  trailing vortices (Ch. 14).
- **Traps.** Thinking the side has no flux because ω = 0 there (it is ω·n that vanishes). Dropping the minus sign of the
  lower lid. Believing "cannot end" needs an inviscid fluid.

### D02 · Pressure in a steadily rotating tank: Euler → (5.5a, b) → (5.6) — ★★, 10 steps, in C02 (notebook · `vortex_pressure_funnel`)
- **Goal.** Find the pressure everywhere in a tank of water turning like a solid body, and from it the bowl shape of
  every surface of constant pressure — including the free surface.
- **Start.** $\mathbf u=u_\theta\mathbf e_\theta,\ u_\theta=\omega r/2$ (5.1) and Cauchy's equation
  $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ (4.24) — *in words:* the fluid turns rigidly; ω is the vorticity,
  so the tank turns at ω/2 (⚠️ Fig. 5.2's label "2ω" is a slip).
- **Plan.** (1) Show the stress is pure pressure, so Euler applies. (2) Compute the particle acceleration in polar
  coordinates. (3) Write the r- and z-components, (5.5a, b). (4) Integrate each and match them.
- **Tools.** Newtonian stress with S = 0 (R05) · derivative of a rotating unit vector ∂e_θ/∂θ = −e_r (ch03 P88, ch04
  P124, reminders) · partial integration with an unknown function (primer in C02) · hydrostatics (R06).
- **Assumptions.** Steady (step 3); axisymmetric, ∂/∂θ = 0 (step 3); gravity −g e_z with z up (step 6); ρ and g
  constant (steps 8–9).
- **Steps.**
  1. *did:* Note there is no viscous stress · *tex:* $S_{ij}=0\ \Rightarrow\ \tau_{ij}=-p\,\delta_{ij}$ · *why:* Rigid rotation does
     not deform elements (ch03), and the Newtonian law (4.37) makes the viscous stress proportional to S. So Cauchy
     (4.24) becomes Euler (4.41). · *plain:* In a rigidly turning tank only pressure pushes on a parcel.
  2. *did:* Write Euler's equation for this flow · *tex:* $\rho\frac{D\mathbf u}{Dt}=-\nabla p-\rho g\,\mathbf e_z$ · *why:* Euler (4.41)
     with g = −g e_z. We need Du/Dt in polar components to read off each balance. · *plain:* Mass times acceleration
     equals the pressure push plus the weight.
  3. *did:* Keep only the turning of the velocity · *tex:* $\frac{D\mathbf u}{Dt}=\frac{u_\theta}{r}\frac{\partial}{\partial\theta}\big(u_\theta\mathbf e_\theta\big)$ ·
     *why:* The book skips this: steady (∂/∂t = 0) and u = u_θ e_θ, so (u·∇) = (u_θ/r)∂/∂θ. u_θ does not depend on θ,
     but the unit vector e_θ does. · *plain:* A particle keeps its speed; only its direction changes.
  4. *did:* Differentiate the turning unit vector · *tex:* $\frac{D\mathbf u}{Dt}=-\frac{u_\theta^2}{r}\,\mathbf e_r$ · *why:* ∂e_θ/∂θ = −e_r
     (the unit vector turns toward the axis; ch03 P88). This is the centripetal acceleration, pointing inward. ·
     *plain:* Going round a circle means accelerating toward its centre. · *set:* {mode: 'tank', Om: 5}.
  5. *did:* Take the radial component · *tex:* $-\rho u_\theta^2/r=-\partial p/\partial r$ · *why:* Dot step 2 with e_r using step 4.
     This is (5.5a): only a pressure rising outward can supply the inward acceleration. · *plain:* Pressure must grow
     with radius to hold every parcel on its circle. · *live:* "at r = 0.1 m: ρu_θ²/r = 1000 × 0.5²/0.1 = 2500 Pa/m".
  6. *did:* Take the vertical and azimuthal components · *tex:* $0=-\partial p/\partial z-\rho g,\qquad \partial p/\partial\theta=0$ · *why:*
     The acceleration has no z or θ part. The z-balance is (5.5b), plain hydrostatics (R06); the θ-balance says p is
     axisymmetric. · *plain:* Vertically the water just holds up its own weight.
  7. *did:* Substitute the solid-body speed · *tex:* $\frac{\partial p}{\partial r}=\frac{\rho\omega^2r}{4}$ · *why:* u_θ = ωr/2 (5.1), so
     u_θ²/r = ω²r/4. Using ω as the rotation rate here would give ω²r, four times too big — the classic trap. ·
     *plain:* The outward pressure gradient grows linearly with radius.
  8. *did:* Integrate in r at fixed z · *tex:* $p(r,z)=\frac{\rho\omega^2r^2}{8}+f(z)$ · *why:* Integrating a partial derivative
     leaves a "constant" that may depend on the variable held fixed (primer: partial integration). The book writes this
     line; it is the first of two. · *plain:* The pressure is a bowl in r plus an unknown function of height.
  9. *did:* Fix f(z) with the vertical balance · *tex:* $f'(z)=-\rho g\ \Rightarrow\ f(z)=p_o-\rho gz$ · *why:* Differentiate step 8
     in z and compare with step 6: only f depends on z. p_o is the pressure at r = 0, z = 0. · *plain:* The unknown
     function is ordinary hydrostatic pressure.
  10. *did:* Assemble and solve for an isobar · *tex:* $p-p_o=\tfrac18\rho\omega^2r^2-\rho gz\ \Rightarrow\ z=\frac{\omega^2r^2}{8g}-\frac{p-p_o}{\rho g}$ ·
      *why:* Put f into step 8: that is (5.6). Holding p fixed and solving for z gives the height of an isobar. ·
      *plain:* Every surface of constant pressure — the free surface too — is a paraboloid. · *live:* "rim–centre
      height ω²R²/8g = 10² × 0.1²/(8 × 9.81) = 12.7 mm".
- **Result.** $p(r,z)-p_o=\tfrac18\rho\omega^2r^2-\rho gz$ (5.6), isobars $z=\frac{\omega^2r^2}{8g}-\frac{p(r,z)-p_o}{\rho g}$ — *in words:*
  the pressure rises with radius just enough to turn every parcel, so constant-pressure surfaces are paraboloids.
- **Check.** Units: ρω²r² is kg m⁻³ s⁻² m² = Pa ✓. ∂p/∂r = ρω²r/4 = ρu_θ²/r ✓ (5.5a); ∂p/∂z = −ρg ✓ (5.5b). ω = 0 gives
  hydrostatics ✓. Numbers: tank at 5 rad/s (ω = 10 s⁻¹), R = 0.1 m: 12.7 mm rim–centre height ✓; volume conservation
  puts the vertex 6.37 mm below the rest level (code).
- **What it means.** The same inward pressure force turns every vortex — cyclostrophic balance in a tornado, the gradient
  wind of Ch. 13. Since u_θ²/2 = ω²r²/8, (5.6) reads −½u_θ² + gz + p/ρ = const: the Bernoulli function B = ½u_θ² + gz
  + p/ρ grows outward (R07) — as it must in a rotational flow.
- **Traps.** Using ω as the tank's rotation rate (Fig. 5.2's "2ω" invites it). Getting −u_θ²/r from ∂u_θ/∂θ (it is
  zero) instead of from the turning e_θ. Treating f as a constant.

### D03 · The line vortex has viscous stress but no net viscous force — ★★, 9 steps, in C02 (notebook · `vortex_pressure_funnel`)
- **Goal.** Compute the viscous stress in the irrotational vortex $u_\theta=\Gamma/2\pi r$ (5.2), show it is not zero, and show
  that the forces it puts on a fluid element nevertheless cancel — the book states it and leaves the force to an
  exercise; we do all of it, three ways.
- **Start.** $u_r=0,\ u_\theta=\Gamma/2\pi r$ (5.2) and the incompressible Newtonian stress $\sigma_{ij}=2\mu S_{ij}$ ((4.59) with
  ∇·u = 0) — *in words:* particles move on circles slower farther out; viscous stress is twice μ times the strain rate.
- **Plan.** (1) Evaluate the polar strain rate S_rθ. (2) Multiply by 2μ. (3) Take the divergence of the stress in polar
  coordinates. (4) Confirm with the Laplacian form and with −μ∇×ω.
- **Tools.** Polar strain rate and polar stress divergence (gloss, `core.curvilinear` sympy check) · viscous force
  $\mu\nabla^2\mathbf u=-\mu\nabla\times\boldsymbol\omega$ (4.40) (Ch. 4 §4.6) · product rule (ch01 P38, reminder).
- **Assumptions.** Incompressible, constant μ (steps 4, 8, 9); r > 0 (all steps — the axis is singular).
- **Steps.**
  1. *did:* Write the polar shear strain rate · *tex:* $S_{r\theta}=\tfrac12\Big[\frac1r\frac{\partial u_r}{\partial\theta}+r\frac{\partial}{\partial r}\Big(\frac{u_\theta}{r}\Big)\Big]$ ·
     *why:* Appendix-B form of the strain rate in polar coordinates (gloss; sympy prints it). The r ∂(u_θ/r)/∂r part
     removes the rigid rotation u_θ ∝ r, which does not deform. · *plain:* The shear rate measures how fast right angles
     at a point close up.
  2. *did:* Drop the radial-velocity term · *tex:* $S_{r\theta}=\tfrac12\,r\frac{\partial}{\partial r}\Big(\frac{\Gamma}{2\pi r^2}\Big)$ · *why:* u_r = 0
     everywhere, so its θ-derivative vanishes; u_θ/r = Γ/(2πr²). · *plain:* Only the change of the angular speed with
     radius deforms elements.
  3. *did:* Differentiate in r · *tex:* $S_{r\theta}=\tfrac12\,r\Big(-\frac{\Gamma}{\pi r^3}\Big)=-\frac{\Gamma}{2\pi r^2}$ · *why:* d(r⁻²)/dr =
     −2r⁻³. The angular speed falls outward, so neighbouring rings slide past each other. · *plain:* Every element is
     being sheared, more strongly near the axis.
  4. *did:* Multiply by 2μ for the stress · *tex:* $\sigma_{r\theta}=2\mu S_{r\theta}=-\frac{\mu\Gamma}{\pi r^2}\neq0$ · *why:* Newtonian
     incompressible law σ = 2μS ((4.59), ∇·u = 0). This is the book's σ_rθ line. · *plain:* The irrotational vortex is
     full of viscous stress. · *live:* "at r = 0.1 m, Γ = 1 m²/s, water: σ_rθ = −0.0318 Pa".
  5. *did:* Check the normal stresses vanish · *tex:* $S_{rr}=\frac{\partial u_r}{\partial r}=0,\quad S_{\theta\theta}=\frac1r\frac{\partial u_\theta}{\partial\theta}+\frac{u_r}{r}=0$ ·
     *why:* u_r = 0 and u_θ is independent of θ; so σ_rr = σ_θθ = 0 and only σ_rθ can produce a force. · *plain:* There
     is no viscous push or pull, only viscous sliding.
  6. *did:* Write the net viscous force, θ-part · *tex:* $f_\theta=\frac{1}{r^2}\frac{\partial}{\partial r}\big(r^2\sigma_{r\theta}\big)$ · *why:*
     Polar divergence of a stress with only σ_rθ(r) nonzero (gloss). The r² comes from the element's two curved faces
     having different areas and lever arms — dropping it is the classic error. · *plain:* The force is how much the
     sliding stress on the outer face beats the inner one.
  7. *did:* Evaluate: the bracket is constant · *tex:* $r^2\sigma_{r\theta}=-\frac{\mu\Gamma}{\pi}\ \Rightarrow\ f_\theta=0$ · *why:* The r²
     cancels the 1/r² of step 4, so the derivative of a constant is zero; the radial part is zero because σ_rr = σ_θθ = 0. ·
     *plain:* The stresses on opposite faces of every element cancel exactly. · *set:* {mode: 'line', Gamma: 1} ·
     *watch:* "the rose force bar stays at 0 while σ_rθ is not".
  8. *did:* Confirm with the Laplacian form · *tex:* $\mu(\nabla^2\mathbf u)_\theta=\mu\Big(u_\theta''+\frac{u_\theta'}{r}-\frac{u_\theta}{r^2}\Big)=\frac{\mu\Gamma}{2\pi}\Big(\frac{2}{r^3}-\frac{1}{r^3}-\frac{1}{r^3}\Big)=0$ ·
     *why:* Incompressible, constant μ: the net viscous force is μ∇²u; its θ-component for u_θ(r) is Appendix B's
     vector Laplacian. The three terms cancel. · *plain:* A second route, same zero.
  9. *did:* Confirm with the vorticity form · *tex:* $\mu\nabla^2\mathbf u=-\mu\nabla\times\boldsymbol\omega=\mathbf 0\quad(\boldsymbol\omega=0,\ r>0)$ ·
     *why:* (4.40): with ∇·u = 0, μ∇²u = −μ∇×ω. The line vortex is irrotational off the axis, so the force vanishes
     without any computing. · *plain:* No vorticity, no net viscous force.
- **Result.** $\sigma_{r\theta}=-\mu\Gamma/\pi r^2\neq0$ but $\mu\nabla^2\mathbf u=\mathbf 0$ for r > 0 — *in words:* the irrotational
  vortex is sheared everywhere, yet no element feels a net viscous force; momentum obeys Euler again, so its pressure
  follows from the same radial balance, $p-p_\infty=-\rho\Gamma^2/8\pi^2r^2-\rho gz$ (5.7).
- **Check.** Units: μΓ/r² is Pa s · m²/s / m² = Pa ✓. Three routes agree (0, 0, 0) ✓. Contrast: the Gaussian vortex has
  ω ≠ 0 and a nonzero force μ∂ω_z/∂r (code: `vortex_stress_force("gaussian", …)`) ✓. The torque 2πr²σ_rθ = −2μΓ is the
  same at every radius (N08) ✓.
- **What it means.** **Irrotational means no net viscous force, not no viscous stress.** Viscosity still dissipates
  energy here (N08) — the torque applied at a rotating cylinder travels outward undiminished. It fails for
  compressible flow (an extra ∇(∇·u) force) and on the axis.
- **Traps.** Dropping the r² in step 6. Confusing stress (per area, on one face) with force (per volume, the imbalance
  of faces). Calling σ here the Gaussian core radius of Ch. 3.

### D04 · The rate of change of the circulation of a material loop: (5.9) — ★★, 8 steps, in C03 (notebook · `kelvin_material_loop`)
- **Goal.** Differentiate the circulation of a loop that moves with the fluid, and show that its stretching and turning
  contribute nothing — only the acceleration of the fluid along the loop counts.
- **Start.** $\Gamma(t)=\oint_{C(t)}u_i\,dx_i$ (3.18) with C(t) a material loop — *in words:* add up the velocity along a loop
  of dye that is carried by the flow.
- **Plan.** (1) Label the loop's particles so the integral has fixed limits. (2) Differentiate under the integral with
  the product rule. (3) Show the change of a loop element is the velocity difference du. (4) Show ∮u·du vanishes.
- **Tools.** Material label (gloss) and functions with parameters (ch03 P89) · differentiation under the integral sign
  (ch03 P109, reminder) · material line element D(δx)/Dt = δu (ch03 P99) · closed-loop integral of an exact
  differential (primer in C03).
- **Assumptions.** C is material (step 1); u smooth and single-valued (steps 3, 8).
- **Steps.**
  1. *did:* Label the loop by fixed particle tags · *tex:* $\mathbf x=\mathbf X(s,t),\quad 0\le s<1,\quad \mathbf X(0,t)=\mathbf X(1,t)$ ·
     *why:* The book skips this: each particle keeps its tag s while it moves, so the loop's parameter range [0, 1)
     never changes even though the loop does (material label gloss). · *plain:* Number the dye particles once; the
     numbers travel with them. · *set:* {flow: 'cellular', t: 0}.
  2. *did:* Write Γ with fixed limits · *tex:* $\Gamma(t)=\int_0^1u_i\big(\mathbf X(s,t),t\big)\,\frac{\partial X_i}{\partial s}\,ds$ · *why:*
     Along the loop dx_i = (∂X_i/∂s)ds (chain rule on a parametrised curve, ch03 P92). Now the limits no longer depend on
     t. · *plain:* The circulation is an ordinary integral over the particle tags.
  3. *did:* Differentiate under the integral sign · *tex:* $\frac{d\Gamma}{dt}=\int_0^1\Big[\frac{Du_i}{Dt}\frac{\partial X_i}{\partial s}+u_i\frac{\partial}{\partial t}\Big(\frac{\partial X_i}{\partial s}\Big)\Big]ds$ ·
     *why:* Fixed limits allow d/dt inside (P109); product rule on the integrand. At fixed s the time derivative of
     u(X(s, t), t) follows one particle, so it is Du_i/Dt (3.5). · *plain:* Γ changes because the fluid accelerates and
     because the loop's elements change.
  4. *did:* Return to loop notation · *tex:* $\frac{D\Gamma}{Dt}=\oint_C\frac{Du_i}{Dt}dx_i+\oint_Cu_i\frac{D}{Dt}(dx_i)$ · *why:*
     (∂X_i/∂s)ds = dx_i, and ∂/∂t at fixed tag is D/Dt of the element. This is (5.9). · *plain:* Acceleration term plus
     element-change term. · *set:* {t: 1.5} · *watch:* "the acceleration curve (teal) and the contour curve (grey)".
  5. *did:* Follow both ends of one element · *tex:* $\mathbf u+d\mathbf u=\frac{D}{Dt}(\mathbf x+d\mathbf x)=\frac{D\mathbf x}{Dt}+\frac{D}{Dt}(d\mathbf x)$ ·
     *why:* The points x and x + dx are both particles (Fig. 5.4), so each moves at the local velocity: u at x, u + du at
     x + dx (first-order Taylor, ch03 P98). · *plain:* The element's tail and head each ride with the flow.
  6. *did:* Subtract the tail's velocity · *tex:* $\frac{D}{Dt}(dx_i)=du_i$ · *why:* Dx/Dt = u, so the remainder is du: a material
     element grows at the velocity difference across it (ch03 P99). · *plain:* An element stretches and turns at the
     rate its two ends move apart.
  7. *did:* Rewrite the second integral · *tex:* $\oint_Cu_i\frac{D}{Dt}(dx_i)=\oint_Cu_i\,du_i=\oint_Cd\big(\tfrac12u_iu_i\big)$ · *why:*
     u_i du_i = d(½u_iu_i) (chain rule for the kinetic energy, ch04 P128). We want an exact differential. · *plain:* The
     element-change term is the change of ½u² round the loop.
  8. *did:* Close the loop: exact differential vanishes · *tex:* $\frac{D\Gamma}{Dt}=\oint_C\frac{Du_i}{Dt}dx_i$ · *why:* ½u² is
     single-valued, so going once round the loop brings it back to its start: ∮d(½u²) = 0 (primer). · *plain:* However
     the loop is stretched or tangled, only the fluid's acceleration along it can change Γ. · *live:* "contour term
     = 3e-13 m²/s² for your loop".
- **Result.** $\frac{D\Gamma}{Dt}=\oint_C\frac{Du_i}{Dt}dx_i$, from (5.9) $\frac{D\Gamma}{Dt}=\oint_C\frac{Du_i}{Dt}dx_i+\oint_Cu_i\frac{D}{Dt}(dx_i)$ with the
  second term zero — *in words:* a loop's circulation changes only if the fluid on it is accelerated along it.
- **Check.** Units: m²/s² ✓. Numerically (`kelvin_rate_terms` on the cellular flow): contour term ≈ 1e-13, acceleration
  term = dΓ/dt from central differences of Γ(t) ✓. A fixed (non-material) loop fails the argument — step 1 needs the
  particles.
- **What it means.** Circulation is a property the fluid carries, not the loop's shape: stretching a loop seven-fold
  changes nothing by itself. Everything that follows (Kelvin, Helmholtz, Bjerknes) is about the acceleration term.
- **Traps.** Moving d/dt inside a loop integral whose limits move (only fixed tags allow it). Using a multi-valued
  function (the polar angle θ round the origin) in step 8 — ∮dθ = 2π.

### D05 · Momentum in (5.9): (5.10) → (5.11) → Kelvin's theorem (5.8), and the three sources of vorticity — ★★, 9 steps, in C03 (notebook · `kelvin_material_loop`)
- **Goal.** Put the momentum equation into $\frac{D\Gamma}{Dt}=\oint_C\frac{Du_i}{Dt}dx_i$ and find exactly which forces can change a
  material loop's circulation — then list the hypotheses under which none can.
- **Start.** $\frac{D\Gamma}{Dt}=\oint_C\frac{Du_i}{Dt}dx_i$ (D04) — *in words:* only the acceleration along the loop matters.
- **Plan.** (1) Replace Du/Dt by the forces per unit mass. (2) Turn the pressure and body-force terms into loop integrals
  of differentials. (3) Show they vanish for barotropic, conservative flow. (4) Read off what is left.
- **Tools.** Cauchy's equation with τ = −pδ + σ (Ch. 4 §4.4–4.5) · conservative force and potential g = −∇Φ (4.18) (ch04
  P115, reminder) · barotropic pressure function (R09, (4.67)) · closed-loop integral of an exact differential (primer
  in C03) · Stokes' theorem (2.34) (Ch. 2 §2.13).
- **Assumptions.** Inertial frame (step 1); conservative body force (step 5); barotropic (step 6); inviscid (step 8).
- **Steps.**
  1. *did:* Substitute the momentum equation · *tex:* $\frac{Du_i}{Dt}=-\frac1\rho\frac{\partial p}{\partial x_i}+g_i+\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j}$ ·
     *why:* Cauchy (4.24) with the stress split τ = −pδ + σ (4.59), divided by ρ, in an inertial frame. The book
     differentiates σ on its second index; harmless because σ is symmetric. · *plain:* Acceleration = pressure push +
     body force + viscous force, each per kilogram.
  2. *did:* Split into three loop integrals · *tex:* $\oint_C\frac{Du_i}{Dt}dx_i=\oint_C\Big(-\frac1\rho\frac{\partial p}{\partial x_i}+g_i+\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j}\Big)dx_i$ ·
     *why:* Linearity of the integral. We examine each force separately. · *plain:* Three forces, three contributions.
  3. *did:* Write the pressure term as a differential · *tex:* $\frac{\partial p}{\partial x_i}dx_i=dp\ \Rightarrow\ -\oint_C\frac1\rho\,dp$ ·
     *why:* The change of p along a small step dx is ∇p·dx (chain rule along a path, ch03 P91). · *plain:* Along the
     loop, pressure changes by dp per step.
  4. *did:* Write gravity through its potential · *tex:* $g_i\,dx_i=-\frac{\partial\Phi}{\partial x_i}dx_i=-d\Phi$ · *why:* g = −∇Φ (4.18)
     for a conservative body force (Φ = gz). With steps 1–3 this is (5.10):
     $\oint_C\frac{Du_i}{Dt}dx_i=-\oint_C\frac1\rho dp-\oint_Cd\Phi+\oint_C\big(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j}\big)dx_i$. · *plain:* Gravity's work
     along a step is the drop in potential.
  5. *did:* Close the loop on the body force · *tex:* $\oint_Cd\Phi=0$ · *why:* Φ is single-valued, so an exact differential
     integrates to zero round a closed loop (primer). A non-conservative force (Coriolis) has no such Φ. · *plain:*
     Gravity does no net work round a loop.
  6. *did:* Close the loop on the pressure term · *tex:* $\rho=\rho(p)\ \Rightarrow\ \frac{dp}{\rho}=d\mathcal P\ \Rightarrow\ \oint_C\frac{dp}{\rho}=0$ · *why:*
     Barotropic: the pressure function 𝒫(p) = ∫dp/ρ(p) (4.67) makes dp/ρ exact. The book's reason, "ρ and p are single
     valued", is not enough: in a baroclinic field both are single-valued and ∮dp/ρ ≠ 0. · *plain:* Only if density
     is tied to pressure does pressure do no net work. · *set:* {flow: 'baroclinic'} · *watch:* "the orange ∮dp/ρ
     curve is not zero here".
  7. *did:* Keep what survives · *tex:* $\frac{D\Gamma}{Dt}=\oint_C\Big(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j}\Big)dx_i$ · *why:* Steps
     5–6 removed two of the three integrals of (5.10). This is (5.11). · *plain:* Only the net viscous force along the
     loop is left. · *set:* {flow: 'viscous'} · *live:* "viscous integral −3.35e-4 = dΓ/dt −3.35e-4 m²/s²".
  8. *did:* Let the fluid be inviscid · *tex:* $\frac{D\Gamma}{Dt}=0$ · *why:* σ = 0 when μ = μ_v = 0 — or when the net viscous force
     vanishes all along C. This is Kelvin's theorem (5.8). · *plain:* In an ideal, barotropic fluid with conservative
     forces, a material loop keeps its circulation. · *set:* {flow: 'cellular'}.
  9. *did:* Read off every way to break it · *tex:* $\frac{D\Gamma}{Dt}\neq0\ \Leftarrow\ \oint\mathbf g\cdot d\mathbf x\ne0,\ \oint\frac{dp}{\rho}\ne0,\ \oint\frac{1}{\rho}\frac{\partial\sigma_{ij}}{\partial x_j}dx_i\ne0$ ·
     *why:* Undo one assumption at a time: non-conservative forces (Coriolis in a rotating frame), baroclinicity (by
     Stokes −∮dp/ρ = ∫(∇ρ×∇p/ρ²)·n dA), net viscous force. "C in irrotational fluid" kills the last only if incompressible, constant μ. · *plain:*
     Three sources of vorticity, one per broken hypothesis.
- **Result.** (5.10) $\oint_C\frac{Du_i}{Dt}dx_i=-\oint_C\frac1\rho dp-\oint_Cd\Phi+\oint_C\big(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j}\big)dx_i$ →
  (5.11) $\frac{D\Gamma}{Dt}=\oint_C\big(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j}\big)dx_i$ → (5.8) $D\Gamma/Dt=0$ — *in words:* in an
  inviscid, barotropic fluid with conservative body forces, seen from an inertial frame, circulation round a material
  loop is conserved.
- **Check.** Units: every term m²/s² ✓. Lamb–Oseen (Γ₀ = 0.01 m²/s, ν = 1e-6 m²/s, r = 5 mm, t = 10 s): the (5.11)
  integral equals ∂Γ/∂t = −3.345e-4 m²/s² ✓. Baroclinic lock-exchange field: ∮dp/ρ round the square equals −∫(∇ρ×∇p/ρ²)·n dA to 1e-10 ✓.
- **What it means.** Circulation — hence vorticity (Stokes) — can be made or destroyed only by viscosity, baroclinicity
  or non-conservative forces; flow that starts irrotational in an ideal fluid stays irrotational (Ch. 6). The rotating
  frame adds the Coriolis term, which C11 turns into (5.33).
- **Traps.** Accepting "single-valued" as the reason in step 6. Forgetting that C must be material. Thinking the
  viscous term vanishes whenever the fluid near C is irrotational (compressible flow keeps a ∇(∇·u) force).

### D06 · The pressure torque on a small element is the baroclinic vorticity source — ★★, 10 steps, in C04 (notebook · `baroclinic_torque`)
- **Goal.** Explain Fig. 5.6 with numbers: compute the torque that pressure exerts about the centre of mass of a small
  fluid disc when isobars and isopycnals cross, and show the resulting spin-up rate is exactly
  $\frac1{\rho^2}\nabla\rho\times\nabla p$ — the baroclinic term the book derives later as (5.28). The book never writes this out.
- **Start.** A disc of radius R about the origin (unit depth), with $\rho=\rho_0+\mathbf x\cdot\nabla\rho$ and
  $p=p_0+\mathbf x\cdot\nabla p$ (∇ρ, ∇p constant across the small disc) — *in words:* over a small enough element both
  fields are linear.
- **Plan.** (1) Find the pressure force and the point it acts through. (2) Find the centre of mass. (3) Take the torque
  about the centre of mass. (4) Divide by the moment of inertia and double.
- **Tools.** Net force from pressure (ch01 P28, reminder) · torque, moment arm and moment of inertia (ch04 P116,
  reminder) · centre of mass of a non-uniform body and I = ½MR² for a disc (gloss) · vorticity = twice the spin (ch03,
  R01) · ∫x_ix_j dA = (πR⁴/4)δ_ij over a disc (gloss).
- **Assumptions.** Small disc, linear fields (steps 2, 5); gravity acts through the centre of mass (step 6); ρ₀ ≫
  R\|∇ρ\| (step 8).
- **Steps.**
  1. *did:* Add the pressure pushes round the rim · *tex:* $\mathbf F=-\oint p\,\mathbf n\,ds$ · *why:* Pressure pushes along the
     inward normal on every piece of the rim (P28). · *plain:* The net pressure force is the sum of all the inward
     arrows of Fig. 5.6. · *set:* {mode: 'element', theta: 90}.
  2. *did:* Use the gradient theorem · *tex:* $\mathbf F=-\int\nabla p\,dA=-\pi R^2\,\nabla p$ · *why:* ∮p n ds = ∫∇p dA (Gauss for a
     scalar, ch02); ∇p is constant. · *plain:* The force points down the pressure gradient.
  3. *did:* Find where the force acts · *tex:* $\mathbf M_O=\oint\mathbf x\times(-p\,\mathbf n)\,ds=\mathbf 0\quad(\mathbf x=R\,\mathbf n)$ · *why:* On a
     circle every rim point's position is parallel to its normal, so each push has zero moment about the centre O. ·
     *plain:* The pressure force passes straight through the disc's geometric centre.
  4. *did:* Compute the disc's mass · *tex:* $M=\int\rho\,dA=\rho_0\pi R^2$ · *why:* The linear part x·∇ρ is odd across the
     centre, so it integrates to zero. · *plain:* The heavy side's extra mass is the light side's deficit.
  5. *did:* Locate the centre of mass · *tex:* $\mathbf x_G=\frac1M\int(\mathbf x\cdot\nabla\rho)\,\mathbf x\,dA=\frac{R^2}{4\rho_0}\nabla\rho$ · *why:*
     Centre of mass definition (gloss); ∫x_ix_j dA = (πR⁴/4)δ_ij over a disc. · *plain:* The centre of mass sits a
     little toward the heavy side. · *watch:* "G slides off the centre toward the dense side".
  6. *did:* Move the torque to the centre of mass · *tex:* $\mathbf M_G=\mathbf M_O+(\mathbf 0-\mathbf x_G)\times\mathbf F=-\mathbf x_G\times\mathbf F$ · *why:*
     Changing the reference point adds (r_O − r_G) × F (P116). Gravity acts at G, so it adds no torque about G. ·
     *plain:* A force through O pulls with a lever arm about G.
  7. *did:* Substitute force and offset · *tex:* $\mathbf M_G=\frac{\pi R^4}{4\rho_0}\,\nabla\rho\times\nabla p$ · *why:* −(R²∇ρ/4ρ₀) × (−πR²∇p);
     the two minus signs cancel. The torque is zero when ∇ρ ∥ ∇p (barotropic). · *plain:* Crossed isolines give a
     torque; parallel ones give none. · *live:* "π(0.01)⁴ × (−98 100)/(4 × 1000) = −7.70e-7 N m/m".
  8. *did:* Compute the moment of inertia · *tex:* $I_G\simeq\int\rho\,\lvert\mathbf x\rvert^2dA=\frac{\rho_0\pi R^4}{2}$ · *why:* Disc I = ½MR²
     (gloss); the density tilt and the offset of G change it only by a fraction \|∇ρ\|²R²/8ρ₀² (exact form in the code). ·
     *plain:* The disc's resistance to spinning, as if uniform.
  9. *did:* Divide torque by moment of inertia · *tex:* $\dot\Omega_{el}=\frac{M_G}{I_G}=\frac{\nabla\rho\times\nabla p}{2\rho_0^2}$ · *why:*
     Newton's law for rotation: angular acceleration = torque / moment of inertia (P116). The R⁴ cancels — the result
     does not depend on the element's size. · *plain:* The element's spin rate grows at a size-independent rate.
  10. *did:* Double spin to get vorticity · *tex:* $\frac{D\boldsymbol\omega}{Dt}=2\dot\Omega_{el}=\frac1{\rho^2}\nabla\rho\times\nabla p$ · *why:*
      Vorticity is twice the angular velocity of an element (R01). This is the baroclinic term of (5.28)/(5.30). ·
      *plain:* Crossed pressure and density surfaces spin fluid up at ∇ρ × ∇p/ρ². · *live:* "−0.0981 s⁻² by the
      torque, −0.0981 s⁻² by the formula".
- **Result.** $\frac{D\boldsymbol\omega}{Dt}\Big\rvert_{\text{pressure torque}}=\frac1{\rho^2}\nabla\rho\times\nabla p$ — *in words:* when surfaces of
  equal density cross surfaces of equal pressure, pressure pushes through the geometric centre while weight hangs off
  it, and the element spins up.
- **Check.** Units: (kg m⁻⁴)(Pa m⁻¹)/(kg m⁻³)² = s⁻² ✓. ∇ρ ∥ ∇p gives zero ✓. Numbers: R = 1 cm, ∇ρ = 10 kg/m⁴ along x,
  ∇p = −9810 Pa/m along y: −0.0981 s⁻² both routes; the numeric disc (4000 rim points, 400 × 800 area cells) agrees to
  1e-9 relative and the error of the torque route falls as R² ✓.
- **What it means.** The mechanism of sea breezes, fronts, the lock exchange and every baroclinic instability (Ch. 13):
  horizontal density contrasts under vertical gravity always cross the isobars.
- **Traps.** Taking torques about O, where the pressure torque is zero, instead of about G. Forgetting the factor 2
  between spin and vorticity. Writing ∇p × ∇ρ (the opposite sign).

### D07 · The initial spin-up rate of a lock exchange — ★★, 7 steps, in C04 (notebook · `baroclinic_torque`)
- **Goal.** Two fluids side by side, heavy ρ₂ on the left and light ρ₁ on the right, separated by an interface of
  thickness δ; the gate is pulled at t = 0. Find how fast vorticity appears at the interface at the first instant — the
  book leaves this to an exercise.
- **Start.** $\frac{D\boldsymbol\omega}{Dt}=\frac{1}{\rho^2}\nabla\rho\times\nabla p$ (D06, the baroclinic part of (5.30)) — *in words:* crossed
  density and pressure gradients spin the fluid.
- **Plan.** (1) Argue every other term is zero at t = 0⁺. (2) Estimate ∇ρ across the interface and ∇p from
  hydrostatics. (3) Cross them. (4) Divide by ρ².
- **Tools.** D06's result · hydrostatics (R06) · e_x × e_y = e_z (ch02 P74, reminder) · smoothed step with `np.tanh`
  (gloss).
- **Assumptions.** Fluid at rest at t = 0⁺ (step 1); thin interface δ ≪ depth (step 2); pressure still hydrostatic with
  the mean density (step 3); ρ₂ − ρ₁ ≪ ρ̄ (step 5).
- **Steps.**
  1. *did:* Note the fluid is still at rest · *tex:* $\mathbf u=\mathbf 0,\ \boldsymbol\omega=\mathbf 0\ \Rightarrow\ \frac{D\boldsymbol\omega}{Dt}=\frac{1}{\rho^2}\nabla\rho\times\nabla p$ ·
     *why:* With u = ω = 0 the stretching, tilting, planetary and diffusion terms of (5.30) all vanish; only the
     baroclinic source can act. · *plain:* At the first instant only the torque of D06 works. · *set:* {mode: 'lock'}.
  2. *did:* Estimate the density gradient · *tex:* $\nabla\rho\approx\frac{\rho_1-\rho_2}{\delta}\,\mathbf e_x$ · *why:* Density drops from ρ₂ to
     ρ₁ across a layer of width δ going in +x (x to the right, y up). A tanh profile has exactly this slope at its
     centre (gloss). · *plain:* Density falls steeply from left to right.
  3. *did:* Use hydrostatic pressure with the mean density · *tex:* $\nabla p\approx-\bar\rho\,g\,\mathbf e_y,\qquad\bar\rho=\tfrac12(\rho_1+\rho_2)$ ·
     *why:* Before anything moves, pressure is hydrostatic (R06); at the interface the weight above is on average ρ̄. ·
     *plain:* Pressure increases downward.
  4. *did:* Take the cross product · *tex:* $\nabla\rho\times\nabla p=\frac{(\rho_2-\rho_1)\bar\rho g}{\delta}\,\mathbf e_z$ · *why:* e_x × e_y = e_z
     (right-hand rule); the two minus signs (ρ₁ − ρ₂ < 0 and −ρ̄g) cancel. · *plain:* The source points out of the
     page: counterclockwise.
  5. *did:* Divide by the density squared · *tex:* $\frac{D\omega_z}{Dt}=\frac{(\rho_2-\rho_1)g}{\bar\rho\,\delta}$ · *why:* at the interface centre ρ = ρ̄ exactly, so ρ² = ρ̄²
     there (off the centre the small contrast, 2.5 % here, keeps ρ² ≈ ρ̄²); one ρ̄ cancels. · *plain:* Rate = reduced gravity divided by the interface
     thickness.
  6. *did:* Substitute the mean density · *tex:* $\frac{D\omega_z}{Dt}=\frac{2(\rho_2-\rho_1)g}{(\rho_2+\rho_1)\delta}$ · *why:* ρ̄ = (ρ₁ + ρ₂)/2. This is
     the form the exercise asks for. · *plain:* The heavier the contrast and the sharper the interface, the faster the
     spin-up. · *live:* "2 × 25 × 9.81/(2025 × 0.1) = 2.42 s⁻²".
  7. *did:* Read the sense of rotation · *tex:* $\frac{D\omega_z}{Dt}>0\ \Rightarrow\ \text{counterclockwise}$ · *why:* Positive ω_z is
     counterclockwise (convention 3): the heavy left fluid slides right along the bottom, the light fluid left along the
     top — the tumbling of Fig. 5.5. · *plain:* The interface tips over, heavy under light. · *watch:* "the curved
     arrows turn counterclockwise".
- **Result.** $\frac{D\omega_z}{Dt}\Big\rvert_{t=0^+}=\frac{2(\rho_2-\rho_1)g}{(\rho_2+\rho_1)\delta}$ — *in words:* the interface starts to tumble
  at the reduced gravity divided by its thickness.
- **Check.** Units: (m/s²)/m = s⁻² ✓. ρ₂ = ρ₁ gives 0 ✓. δ → 0 gives an infinite rate at a vortex sheet (C14) — the
  circulation per unit height, rate × δ = 2(ρ₂−ρ₁)g/(ρ₂+ρ₁), stays finite ✓. Field route (`baroclinic_term` on the tanh
  step with hydrostatic p) = 2.4222 s⁻² at x = 0 ✓.
- **What it means.** Density currents, sea breezes and oceanic fronts begin exactly like this (Ch. 13); the reduced
  gravity g′ of Ch. 4 reappears.
- **Traps.** Using ρ₂ − ρ₁ with the heavy fluid on the right (the sign flips). Letting the pressure gradient have a
  horizontal part at t = 0 (it appears only once the fluid moves).

### D08 · Helmholtz's vortex theorems from Kelvin's theorem and (5.4) — ★★, 9 steps, in C05 (notebook · `kelvin_material_loop`)
- **Goal.** Prove that vortex lines move with the fluid (Helmholtz 1) using only Kelvin's theorem and Stokes' theorem,
  and see that the other three theorems are (5.4) and (5.8) restated.
- **Start.** $D\Gamma/Dt=0$ (5.8) for every material loop, and $-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0$ (5.4) — *in words:*
  circulation is carried by the fluid; tube strength is the same along a tube.
- **Plan.** (1) Put a small patch on a tube's wall. (2) Its edge has zero circulation. (3) Carry the patch with the flow;
  Kelvin keeps the zero. (4) Conclude the carried patch still lies on a tube wall, for every patch.
- **Tools.** Kelvin's theorem (C03) · tube wall has ω·n = 0 (C01) · Stokes' theorem (2.34) (Ch. 2 §2.13) · the small-ball
  (small-patch) argument (ch04 P112, reminder).
- **Assumptions.** Kelvin's four restrictions — inviscid, barotropic, conservative body forces, inertial frame
  (step 4); smooth ω (step 6).
- **Steps.**
  1. *did:* Draw a patch on the tube wall · *tex:* $S\subset\text{tube wall},\quad \partial S=\text{its edge}$ · *why:* Fig. 5.7: a small
     surface lying in the wall that does not wrap round the tube. We will follow it in time. · *plain:* Paint a small
     patch on the side of the tube. · *set:* {flow: 'helmholtz', t: 0}.
  2. *did:* Note vorticity lies along the patch · *tex:* $\boldsymbol\omega\cdot\mathbf n=0\ \text{on}\ S$ · *why:* The wall is made of vortex
     lines (C01), so ω is tangent to it everywhere. · *plain:* No vorticity pierces the patch.
  3. *did:* Apply Stokes to the patch's edge · *tex:* $\Gamma_{\partial S}=\int_S\boldsymbol\omega\cdot\mathbf n\,dA=0$ · *why:* Stokes' theorem
     (2.34): circulation round an edge = vorticity flux through the surface it bounds. · *plain:* A loop drawn round
     the patch has zero circulation.
  4. *did:* Carry the edge with the fluid · *tex:* $\Gamma_{\partial S'}(t)=\Gamma_{\partial S}(0)=0$ · *why:* The edge is a material loop, so
     Kelvin's theorem (5.8) keeps its circulation; its particles now bound a new surface S′. · *plain:* The carried loop
     still has zero circulation. · *set:* {t: 2} · *watch:* "the flux through the small loop stays ≈ 0".
  5. *did:* Apply Stokes to the carried patch · *tex:* $\int_{S'}\boldsymbol\omega\cdot\mathbf n\,dA=0$ · *why:* Stokes again, now on S′. ·
     *plain:* No net vorticity pierces the carried patch.
  6. *did:* Use "for every patch" to localise · *tex:* $\boldsymbol\omega\cdot\mathbf n=0\ \text{at every point of the carried wall}$ · *why:* The
     book skips this: steps 1–5 hold for every patch; if ω·n ≠ 0 somewhere, a tiny patch there would have nonzero flux
     (continuity, P112) — contradiction. · *plain:* Not just on average: nowhere does vorticity pierce the carried
     wall.
  7. *did:* Conclude tube walls and lines are material · *tex:* $\text{carried tube wall}=\text{tube wall}\ \Rightarrow\ \text{vortex lines move with the fluid}$ ·
     *why:* A surface everywhere tangent to ω is a tube wall. Shrinking the tube to zero section turns the statement
     about walls into one about lines. · *plain:* Helmholtz 1: particles on a vortex line stay on a vortex line.
  8. *did:* Recognise theorems 2 and 3 · *tex:* $\Gamma_{\text{lower end}}=\Gamma_{\text{upper end}}$ · *why:* That is (5.4) (D01): kinematic,
     true even with viscosity. It forbids tubes ending in the fluid. · *plain:* Helmholtz 2 and 3: one strength per
     tube, no loose ends.
  9. *did:* Recognise theorem 4 · *tex:* $\frac{D}{Dt}\oint_{C_{\text{tube}}}\mathbf u\cdot d\mathbf x=0$ · *why:* A loop round the tube on its wall is
     material and stays on the wall (step 7), so by Kelvin its circulation — the tube strength — is constant in time. ·
     *plain:* Helmholtz 4: a tube's strength never changes.
- **Result.** Helmholtz's theorems — (1) vortex lines move with the fluid; (2) a tube's strength is the same along it;
  (3) tubes cannot end in the fluid; (4) a tube's strength is constant in time — under Kelvin's four restrictions (1, 4)
  or always (2, 3).
- **Check.** Numerically (after C06): a material element δx started along ω in the ABC flow
  (steady, inviscid, ω = u) stays parallel to ω (angle < 1e-8) and \|ω\|/\|δx\| stays 1 to 1e-8; the viscous Burgers vortex drifts ✓. The small
  ABC loop on a tube wall keeps Γ = 8.0×10⁻¹⁰ m²/s (vs \|ω\|A = 5.4×10⁻⁴) and its normal at 90° to ω ✓.
- **What it means.** Vortices are carried by the flow — the rule behind point vortices, rings and sheets (C12–C14) and
  behind "frozen-in" potential vorticity (Ch. 13). Viscosity (reconnection) or baroclinicity breaks theorems 1 and 4.
- **Traps.** Proving it for one patch only (step 6 needs every patch). Believing theorem 3 needs an ideal fluid.

### D09 · The vorticity equation by the vector route: (5.12) → (5.13) — ★★★, 12 steps, in C06 (notebook · `vorticity_stretching_tilting`)
- **Goal.** Turn the momentum equation into an equation for the vorticity alone, to see what can change a particle's
  vorticity — and why pressure and gravity cannot.
- **Start.** $\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=-\frac1\rho\nabla p+\mathbf g+\nu\nabla^2\mathbf u$ (4.39b), ∇·u = 0 — *in words:*
  incompressible Navier–Stokes; §5.4 takes ρ and ν constant, g = −∇Φ, an inertial frame (N18).
- **Plan.** (1) Take the curl of every term, (5.12). (2) Kill the gradient terms. (3) Rewrite the advective term with
  the Lamb identity. (4) Expand ∇×(ω×u) with (B.3.10) and drop what continuity and ∇·ω = 0 remove. (5) Regroup into
  D/Dt.
- **Tools.** curl of a gradient = 0 (Ch. 2 §2.13) · Schwarz's theorem (ch04 P121, reminder) · Lamb identity
  $(\mathbf u\cdot\nabla)\mathbf u=\nabla(\tfrac12\mathbf u\cdot\mathbf u)+\boldsymbol\omega\times\mathbf u$ (4.68) (R11) · identity (B.3.10) (N20;
  ε–δ (2.19), Ch. 2 §2.7) · ∇·ω = 0 (R10) · material derivative (3.5) (Ch. 3 §3.2).
- **Assumptions.** ρ constant (step 2); g conservative (step 3); ν constant (step 7); ∇·u = 0 (step 10); smooth fields
  (steps 4, 7).
- **Steps.**
  1. *did:* Take the curl of every term · *tex:* $\nabla\times\Big\{\frac{D\mathbf u}{Dt}=-\frac1\rho\nabla p+\mathbf g+\nu\nabla^2\mathbf u\Big\}$ · *why:*
     ω = ∇×u, so the curl of the momentum equation is where an equation for ω must come from. This is (5.12). ·
     *plain:* Ask how each force twists the fluid rather than pushes it. · *set:* {flow: 'burgers'}.
  2. *did:* Kill the pressure term · *tex:* $\nabla\times\Big(-\frac1\rho\nabla p\Big)=-\nabla\times\nabla\Big(\frac p\rho\Big)=\mathbf 0$ · *why:* With ρ
     constant, (1/ρ)∇p = ∇(p/ρ), a gradient, and the curl of a gradient is zero (Ch. 2). If ρ varied this would be the
     baroclinic term (5.28). · *plain:* Pressure pushes through the centre of mass: no twist.
  3. *did:* Kill the gravity term · *tex:* $\nabla\times\mathbf g=-\nabla\times\nabla\Phi=\mathbf 0$ · *why:* Conservative body force g = −∇Φ
     (4.18); again the curl of a gradient. · *plain:* Gravity cannot spin fluid of uniform density.
  4. *did:* Swap curl and time derivative · *tex:* $\nabla\times\frac{\partial\mathbf u}{\partial t}=\frac{\partial}{\partial t}(\nabla\times\mathbf u)=\frac{\partial\boldsymbol\omega}{\partial t}$ ·
     *why:* Mixed partial derivatives commute for smooth fields (Schwarz, P121). · *plain:* The curl of the velocity's
     rate is the vorticity's rate.
  5. *did:* Rewrite advection with the Lamb identity · *tex:* $(\mathbf u\cdot\nabla)\mathbf u=\nabla\big(\tfrac12\mathbf u\cdot\mathbf u\big)+\boldsymbol\omega\times\mathbf u$ ·
     *why:* (4.68) (R11). The book prints ∇(u·u); the correct factor is ½ — harmless, because the next step's curl kills
     the gradient anyway. · *plain:* Advection = gradient of kinetic energy + the Lamb vector.
  6. *did:* Take its curl · *tex:* $\nabla\times[(\mathbf u\cdot\nabla)\mathbf u]=\nabla\times(\boldsymbol\omega\times\mathbf u)$ · *why:* The curl of the
     gradient ∇(½u²) is zero (Ch. 2). · *plain:* Only the Lamb vector can twist.
  7. *did:* Move the curl through the Laplacian · *tex:* $\nabla\times(\nu\nabla^2\mathbf u)=\nu\nabla^2(\nabla\times\mathbf u)=\nu\nabla^2\boldsymbol\omega$ · *why:*
     ν constant and derivatives commute (Cartesian components, Schwarz). · *plain:* Viscosity diffuses vorticity as it
     diffuses momentum.
  8. *did:* Collect the surviving terms · *tex:* $\frac{\partial\boldsymbol\omega}{\partial t}+\nabla\times(\boldsymbol\omega\times\mathbf u)=\nu\nabla^2\boldsymbol\omega$ · *why:*
     Steps 2–7 in (5.12). This is the book's intermediate line. · *plain:* A first vorticity equation, still with a
     tangled middle term. · *live:* "Burgers at R = 1 mm: the three bars and their residual".
  9. *did:* Expand the double cross product · *tex:* $\nabla\times(\boldsymbol\omega\times\mathbf u)=(\mathbf u\cdot\nabla)\boldsymbol\omega-(\boldsymbol\omega\cdot\nabla)\mathbf u+\boldsymbol\omega(\nabla\cdot\mathbf u)-\mathbf u(\nabla\cdot\boldsymbol\omega)$ ·
     *why:* Identity (B.3.10) in its general form (N20), proved with the ε–δ identity (2.19) and the product rule;
     the sympy cell checks it for arbitrary fields. · *plain:* The twist splits into carrying, stretching and two
     divergence terms.
  10. *did:* Drop the two divergence terms · *tex:* $\nabla\times(\boldsymbol\omega\times\mathbf u)=(\mathbf u\cdot\nabla)\boldsymbol\omega-(\boldsymbol\omega\cdot\nabla)\mathbf u$ · *why:*
      Incompressible: ∇·u = 0; and ∇·ω = ∇·(∇×u) = 0 for any flow (R10). The book uses both silently. · *plain:*
      Carrying minus stretching.
  11. *did:* Substitute back · *tex:* $\frac{\partial\boldsymbol\omega}{\partial t}+(\mathbf u\cdot\nabla)\boldsymbol\omega-(\boldsymbol\omega\cdot\nabla)\mathbf u=\nu\nabla^2\boldsymbol\omega$ · *why:*
      Replace the middle term of step 8 by step 10. · *plain:* Nearly there: the first two terms belong together.
  12. *did:* Recognise D/Dt and rearrange · *tex:* $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$ · *why:*
      ∂ω/∂t + (u·∇)ω = Dω/Dt (3.5); move the stretching term to the right. This is (5.13). · *plain:* Following a
      particle, vorticity changes only by stretching/tilting and by diffusion. · *set:* {flow: 'axial'}.
- **Result.** $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$ (5.13) — *in words:* a particle's
  vorticity is changed only by the stretching and tilting of vortex lines and by viscous diffusion; pressure and
  gravity never appear.
- **Check.** Units: s⁻² for every term ✓. 2-D (ω ⟂ plane): (ω·∇)u = 0, so Dω/Dt = ν∇²ω ✓. Lamb–Oseen: stretching 0,
  local = diffusion ✓. Burgers (steady): advective + stretching = diffusion, residual ≈ 1e-6 of the largest term at
  h = 1e-3 ✓. **sympy check intent (★★★):** for a generic smooth field u = (P, Q, R) of (x, y, z, t) and symbols ν, ρ,
  re-run the construction: take the curl of each term of (4.39b) separately (step 1), show the pressure and gravity curls
  simplify to 0 (steps 2–3), verify (B.3.10) component by component for generic ω and u (step 9), then substitute a
  divergence-free polynomial field (u = (y z², x z, −x y)·t-free plus a Taylor–Green piece) and show curl(NS residual)
  − [(5.13) residual] ≡ 0.
  `check_src` sketch (every line commented in the builder):
  `x, y, z, t, nu = sp.symbols('x y z t nu')` (coordinates, time, viscosity) · `P, Q, R = [sp.Function(n)(x, y, z, t)
  for n in 'PQR']` (a generic velocity) · `terms = ch05.vorticity_equation_sym([P, Q, R], (x, y, z), t, nu,
  p_expr=sp.Function('p')(x, y, z, t), Phi_expr=sp.Function('Phi')(x, y, z))` (step 1: the curl of every term) ·
  `print([sp.simplify(c) for c in terms['curl_pressure']], [sp.simplify(c) for c in terms['curl_gravity']])` (steps 2–3:
  zeros) · `print([sp.simplify(c) for c in terms['identity_B310']])` (step 9: zeros) · `u_div0 = [y*z**2, x*z, -x*y]`
  (a divergence-free test field) · `print(sp.simplify(sum(sp.diff(ui, xi) for ui, xi in zip(u_div0, (x, y, z)))))` (0:
  it is incompressible) · `r = ch05.vorticity_equation_sym(u_div0, (x, y, z), t, nu)['residual_513']` (steps 10–12) ·
  `print([sp.simplify(c) for c in r])` (zeros: (5.13) is the curl of NS).
- **What it means.** Vorticity is carried by the flow, amplified or turned by velocity gradients along it, and
  smoothed by viscosity — never created by pressure in a uniform-density fluid. The equation looks linear in ω, but
  u is determined by ω (C07), so both the advective and the stretching terms are nonlinear.
- **Traps.** Dropping ∇·u and ∇·ω without saying why (step 10). A sign slip in ∇×(ω×u) = −∇×(u×ω). Forgetting that
  step 2 needs constant ρ — without it the baroclinic term survives (C09).

### D10 · Velocity from vorticity: ∇×ω = −∇²u and its Green's-function solution (5.14), with the correct sign — ★★★, 12 steps, in C07 (notebook · `biot_savart_filament`)
- **Goal.** Given the vorticity everywhere, find the velocity it induces — an explicit integral. The book prints the
  answer (5.14) with a wrong sign and leaves the proof to "Exercise 5.8" (meaning 5.9); we derive it.
- **Start.** $\boldsymbol\omega=\nabla\times\mathbf u$, $\nabla\cdot\mathbf u=0$, all space, u → 0 far away — *in words:* an incompressible flow whose
  vorticity we know.
- **Plan.** (1) Take the curl of ω and use the curl-of-curl identity to get a Poisson equation for u. (2) Build the
  Green's function of ∇². (3) Superpose point sources. (4) Track the signs.
- **Tools.** curl of a curl (ch04 P122, reminder) · Poisson equation and Green's function (primer in C07) · Gauss'
  theorem on a small sphere (2.30) (Ch. 2 §2.12) · superposition for a linear equation (gloss).
- **Assumptions.** Incompressible (step 3); unbounded domain, u → 0 at infinity (steps 9, 12).
- **Steps.**
  1. *did:* Take the curl of the vorticity · *tex:* $\nabla\times\boldsymbol\omega=\nabla\times(\nabla\times\mathbf u)$ · *why:* We want an equation
     that contains u and the known ω only; the curl of ω brings back second derivatives of u. · *plain:* Look at how
     the vorticity itself curls.
  2. *did:* Apply the curl-of-curl identity · *tex:* $\nabla\times(\nabla\times\mathbf u)=\nabla(\nabla\cdot\mathbf u)-\nabla^2\mathbf u$ · *why:* Identity (B.3.13),
     proved with ε–δ in Ch. 4 (P122). · *plain:* A curl of a curl is a gradient of divergence minus a Laplacian.
  3. *did:* Use incompressibility · *tex:* $\nabla^2\mathbf u=-\nabla\times\boldsymbol\omega$ · *why:* ∇·u = 0 removes the gradient term. This is
     the book's line before (5.14). · *plain:* The velocity obeys a Poisson equation whose source is the curl of the
     vorticity. · *set:* {shape: 'tube'}.
  4. *did:* Treat each component separately · *tex:* $\nabla^2u_i=q_i,\qquad q_i=-(\nabla\times\boldsymbol\omega)_i$ · *why:* In Cartesian
     components the vector Laplacian is three scalar Laplacians; three independent Poisson equations. · *plain:* Three
     copies of one scalar problem.
  5. *did:* Look for a point-source solution · *tex:* $\nabla^2G(\mathbf x,\mathbf x')=\delta(\mathbf x-\mathbf x')$ · *why:* The Green's-function
     idea (primer): solve for one concentrated unit source first, then add sources up. · *plain:* How does the field
     respond to a single pinpoint source?
  6. *did:* Check 1/r is harmonic off the source · *tex:* $\nabla^2\frac{1}{r}=0\quad(r=\lvert\mathbf x-\mathbf x'\rvert>0)$ · *why:* In spherical
     coordinates ∇²f = r⁻²(r²f′)′, and r²·(−1/r²) is constant (sympy check). · *plain:* Away from the point, 1/r has
     no sources.
  7. *did:* Measure the source strength of 1/r · *tex:* $\oint_{\lvert\mathbf x-\mathbf x'\rvert=\epsilon}\nabla\frac1r\cdot\mathbf n\,dA=-\frac{1}{\epsilon^2}\,4\pi\epsilon^2=-4\pi\ \Rightarrow\ \nabla^2\frac1r=-4\pi\,\delta$ ·
     *why:* Gauss (2.30) on a small sphere: the flux of ∇(1/r) is −4π for every radius, so all of 1/r's "source" sits
     at the centre. · *plain:* 1/r is the field of a point source of strength −4π.
  8. *did:* Normalise to a unit source · *tex:* $G(\mathbf x,\mathbf x')=-\frac{1}{4\pi\lvert\mathbf x-\mathbf x'\rvert}$ · *why:* Divide step 7 by
     −4π. · *plain:* Minus one over 4π r is the response to a unit source.
  9. *did:* Superpose all the sources · *tex:* $u_i(\mathbf x)=\int_{V'}G(\mathbf x,\mathbf x')\,q_i(\mathbf x')\,d^3x'$ · *why:* The Poisson equation is
     linear, so the response to many sources is the sum of the responses (gloss); with u → 0 at infinity the solution
     is unique. · *plain:* Add up every point's contribution.
  10. *did:* Insert G and the source · *tex:* $\mathbf u(\mathbf x)=\int_{V'}\Big(-\frac{1}{4\pi\lvert\mathbf x-\mathbf x'\rvert}\Big)\big(-\nabla'\times\boldsymbol\omega(\mathbf x')\big)d^3x'$ ·
      *why:* Steps 4 and 8. ∇′ differentiates with respect to the source point x′. · *plain:* Two minus signs appear.
  11. *did:* Cancel the two minus signs · *tex:* $\mathbf u(\mathbf x,t)=+\frac{1}{4\pi}\int_{V'}\frac{\nabla'\times\boldsymbol\omega(\mathbf x',t)}{\lvert\mathbf x-\mathbf x'\rvert}\,d^3x'$ ·
      *why:* (−)(−) = +. **The book prints −1/(4π) in (5.14): a sign slip** — with it a vortex would swirl backwards
      (code). · *plain:* Vorticity's curl, spread out as 1/distance, gives the velocity. · *set:* {sign: true} · *watch:*
      "the printed sign reverses every arrow".
  12. *did:* Note what may be added · *tex:* $\mathbf u=\mathbf u_\omega+\nabla\phi,\qquad\nabla^2\phi=0$ · *why:* Any irrotational,
      divergence-free flow adds nothing to ∇×ω, so the integral gives only the "vorticity-induced portion"; far-field
      decay removes the rest here (Ch. 6 adds it back). · *plain:* Vorticity fixes the velocity up to a potential flow.
      · *set:* {sign: false}.
- **Result.** $\mathbf u(\mathbf x,t)=\frac{1}{4\pi}\int_{V'}\frac{\nabla'\times\boldsymbol\omega(\mathbf x',t)}{\lvert\mathbf x-\mathbf x'\rvert}\,d^3x'$ — the corrected (5.14),
  which the book prints as $\mathbf u(\mathbf x,t)=-\frac{1}{4\pi}\int_{V'}\frac{1}{\lvert\mathbf x-\mathbf x'\rvert}(\nabla'\times\boldsymbol\omega(\mathbf x',t))d^3x'$ — *in words:*
  the velocity at x is a 1/distance-weighted sum of the curl of the vorticity everywhere.
- **Check.** Units: [∇′×ω] = 1/(m s) and d³x′/\|x − x′\| = m², so u is in m/s ✓. Gaussian tube (Γ = 1 m²/s,
  σ_c = 0.1 m) at r = 0.5 m: the + sign gives u_θ = +0.3183 m/s = Γ/2πr, the printed sign −0.3183 m/s ✓ (quadrature
  24³–40³ nodes, error 1e-4). **sympy check intent (★★★):** re-run steps 2, 6, 7 and 11: (a) the curl-of-curl identity
  for a generic field; (b) ∇²(1/r) = 0 for r > 0; (c) the flux of ∇(1/r) through a sphere of radius ε equals −4π for
  every ε (integrate in spherical coordinates); (d) that G = −1/(4πr) solves ∇²G = 0 off the source with unit flux; then
  numerically apply `velocity_from_curl_omega(sign=±1)` to the Gaussian tube and compare with Γ/2πr.
  `check_src` sketch: `x, y, z, eps = sp.symbols('x y z epsilon', positive=True)` (a point and a sphere radius) ·
  `r = sp.sqrt(x**2 + y**2 + z**2)` (distance from the source at the origin) · `lap = sum(sp.diff(1/r, v, 2) for v in (x,
  y, z))` (step 6) · `print(sp.simplify(lap))` (0) · `th, ph = sp.symbols('theta phi')` (sphere angles) · `flux =
  sp.integrate(sp.integrate(-1/eps**2 * eps**2*sp.sin(th), (th, 0, sp.pi)), (ph, 0, 2*sp.pi))` (step 7: ∂(1/r)/∂r = −1/ε²
  times the area element) · `print(flux)` (−4π) · `for s in (+1, -1): print(s, ch05.velocity_from_curl_omega(F['curl_omega'],
  [0.5, 0, 0], nodes, weights, sign=s)[1])` (step 11: +0.318 / −0.318).
- **What it means.** Vorticity is the "cause" of velocity in the sense a current causes a magnetic field: know ω and you
  know u (plus a potential flow). It is the basis of vortex methods (Ch. 10), lifting-line theory (Ch. 14) and
  potential-vorticity inversion (Ch. 13).
- **Traps.** Taking G = +1/(4πr) (then ∇²G = −δ). Letting ∇ act on x instead of x′ in the source. Trusting the
  printed sign.

### D11 · From (5.14) to Biot–Savart: product rule, Gauss in curl form (5.15), a tube-aligned volume → (5.16) — ★★★, 12 steps, in C07 (notebook · `biot_savart_filament`)
- **Goal.** Rewrite the corrected (5.14) so that ω itself, not its curl, appears in the integrand — the Biot–Savart law
  — and see why a volume hugging the vortex makes the leftover surface term vanish.
- **Start.** $\mathbf u(\mathbf x,t)=\frac1{4\pi}\int_{V'}\frac{\nabla'\times\boldsymbol\omega(\mathbf x',t)}{\lvert\mathbf x-\mathbf x'\rvert}d^3x'$ (D10's result, the corrected
  (5.14)) — *in words:* the velocity as a weighted sum of the curl of the vorticity.
- **Plan.** (1) Move the curl off ω with the product rule. (2) Compute ∇′(1/\|x − x′\|) carefully. (3) Turn the
  pure-curl integral into a surface integral (5.15). (4) Choose V′ so the surface term vanishes.
- **Tools.** Product rule for a curl (primer in C07) · gradient of 1/\|x − x′\| with respect to x′ (primer in C07) ·
  cross-product antisymmetry (gloss) · Gauss' theorem (2.30) applied component by component with ε_kij (Ch. 2 §2.7,
  §2.12).
- **Assumptions.** ω smooth inside V′ (steps 2, 8); V′ built as in Fig. 5.8 (steps 10–11).
- **Steps.**
  1. *did:* Name the kernel · *tex:* $\phi(\mathbf x')=\frac{1}{\lvert\mathbf x-\mathbf x'\rvert}$ · *why:* x is fixed while we integrate over x′, so the
     kernel is a scalar field of x′. · *plain:* Call the 1/distance weight φ.
  2. *did:* Apply the product rule for a curl · *tex:* $\nabla'\times(\phi\,\boldsymbol\omega)=\phi\,\nabla'\times\boldsymbol\omega+\nabla'\phi\times\boldsymbol\omega$ · *why:*
     ∇×(fA) = f∇×A + ∇f × A (primer), with ∇′. We want φ∇′×ω, which is our integrand. · *plain:* Curling a product
     gives two pieces.
  3. *did:* Solve for the integrand · *tex:* $\phi\,\nabla'\times\boldsymbol\omega=\nabla'\times(\phi\,\boldsymbol\omega)-\nabla'\phi\times\boldsymbol\omega$ · *why:* Move
     one term across. This is the book's first line of the rewrite, which it prints correctly. · *plain:* The
     integrand = a pure curl minus a leftover.
  4. *did:* Differentiate the kernel with respect to x′ · *tex:* $\nabla'\frac{1}{\lvert\mathbf x-\mathbf x'\rvert}=+\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ · *why:*
     Chain rule: ∂′_i\|x − x′\| = −(x_i − x′_i)/\|x − x′\|, and d(1/s)/ds = −1/s²; two minus signs (primer, sympy). ·
     *plain:* The kernel grows toward the field point.
  5. *did:* Substitute and reorder the cross product · *tex:* $-\nabla'\phi\times\boldsymbol\omega=-\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}\times\boldsymbol\omega=\boldsymbol\omega\times\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ ·
     *why:* a × b = −b × a (gloss). **The book prints "+((x − x′)/\|x − x′\|³) × ω" here — a second sign slip that
     cancels the one in (5.14).** · *plain:* The leftover is ω crossed into the direction of the field point.
  6. *did:* Split u into two integrals · *tex:* $\mathbf u=\frac1{4\pi}\int_{V'}\nabla'\times\Big(\frac{\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}\Big)d^3x'+\frac1{4\pi}\int_{V'}\frac{\boldsymbol\omega\times(\mathbf x-\mathbf x')}{\lvert\mathbf x-\mathbf x'\rvert^3}d^3x'$ ·
     *why:* Insert steps 3 and 5 into D10's result; linearity. (The book's first integral carries the wrong sign; it
     vanishes below, so no harm.) · *plain:* A pure-curl part and the Biot–Savart part.
  7. *did:* Write the curl in index notation · *tex:* $\Big[\nabla'\times\Big(\frac{\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}\Big)\Big]_k=\varepsilon_{kij}\frac{\partial}{\partial x'_i}\Big(\frac{\omega_j}{\lvert\mathbf x-\mathbf x'\rvert}\Big)$ ·
     *why:* (∇×A)_k = ε_kij ∂A_j/∂x_i (Ch. 2 §2.7). For each fixed k this is a divergence of the vector F_i = ε_kijA_j. ·
     *plain:* Each component of the curl is a divergence in disguise.
  8. *did:* Apply Gauss for each component · *tex:* $\int_{V'}\varepsilon_{kij}\frac{\partial}{\partial x'_i}\Big(\frac{\omega_j}{\lvert\mathbf x-\mathbf x'\rvert}\Big)d^3x'=\int_{A'}\varepsilon_{kij}\Big(\frac{\omega_j}{\lvert\mathbf x-\mathbf x'\rvert}\Big)n_i\,d^2x'$ ·
     *why:* Gauss (2.30) on F_i = ε_kij ω_j/\|x − x′\|, with n the outward normal of A′ = ∂V′. It needs the field
     point x outside V′ (the kernel smooth inside). · *plain:* The volume integral becomes a surface integral.
  9. *did:* Return to vector form · *tex:* $\int_{V'}\nabla'\times\Big(\frac{\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}\Big)d^3x'=\int_{A'}\frac{\mathbf n\times\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}d^2x'$ ·
     *why:* ε_kij n_i ω_j = (n × ω)_k. This is (5.15), Gauss' theorem in curl form (N27; `curl_theorem_box` checks it). ·
     *plain:* A volume curl equals a surface "n cross" integral.
  10. *did:* Cap V′ with ends normal to ω · *tex:* $\mathbf n\parallel\boldsymbol\omega\ \Rightarrow\ \mathbf n\times\boldsymbol\omega=\mathbf 0\ \text{on the ends}$ · *why:*
      Fig. 5.8: we may choose V′; flat ends perpendicular to the local vorticity make the cross product vanish there. ·
      *plain:* The lids contribute nothing.
  11. *did:* Put the side wall outside the vortex · *tex:* $\boldsymbol\omega=\mathbf 0\ \text{on the side}$ · *why:* The curved side is drawn
      where there is no vorticity (outside the core), so the integrand is zero there too. · *plain:* The wall
      contributes nothing. · *watch:* "the dashed volume V′ hugs the filament piece".
  12. *did:* Keep the Biot–Savart integral · *tex:* $\mathbf u(\mathbf x,t)=\frac1{4\pi}\int_{V'}\frac{\boldsymbol\omega(\mathbf x',t)\times(\mathbf x-\mathbf x')}{\lvert\mathbf x-\mathbf x'\rvert^3}d^3x'$ ·
      *why:* Steps 9–11 remove the first integral of step 6. This is (5.16). · *plain:* Each bit of vorticity pushes the
      fluid at x perpendicular to itself and to the line joining them, falling off as 1/distance². · *set:* {shape:
      'segment'}.
- **Result.** $\mathbf u(\mathbf x,t)=\frac1{4\pi}\int_{V'}\frac{\boldsymbol\omega(\mathbf x',t)\times(\mathbf x-\mathbf x')}{\lvert\mathbf x-\mathbf x'\rvert^3}d^3x'$ (5.16) — *in words:*
  the Biot–Savart law: vorticity induces velocity like an electric current induces a magnetic field.
- **Check.** Units: (1/s)(m)/m³ · m³ = m/s ✓. A straight uniform filament gives Γ/2πr (D13) ✓. Rankine and Lamb–Oseen
  planar vortices recovered inside and outside the core by `biot_savart_2d` (1e-6) ✓. (5.15) on a random polynomial
  field: volume = surface to 1e-12 ✓. **sympy check intent (★★★):** re-run steps 2, 4 and 5 symbolically: (a) the
  product-rule identity ∇′×(φω) − φ∇′×ω − ∇′φ×ω = 0 for generic φ(x′) and ω(x′); (b) ∇′(1/\|x − x′\|) − (x − x′)/\|x − x′\|³
  = 0; (c) −∇′φ × ω − ω × (x − x′)/\|x − x′\|³ = 0 (the correct step 5) while the book's "+(x − x′)/r³ × ω" differs by
  2ω×(x − x′)/r³. `check_src` sketch: `xs = sp.symbols('x y z'); xp = sp.symbols('xp yp zp')` (field and source points) ·
  `d = sp.Matrix(xs) - sp.Matrix(xp); r = sp.sqrt(d.dot(d))` · `W = sp.Matrix([sp.Function(f'w{i}')(*xp) for i in
  range(3)])` (generic ω(x′)) · `phi = 1/r` · `curlp = lambda A: sp.Matrix([sp.diff(A[2], xp[1]) - sp.diff(A[1], xp[2]),
  sp.diff(A[0], xp[2]) - sp.diff(A[2], xp[0]), sp.diff(A[1], xp[0]) - sp.diff(A[0], xp[1])])` (∇′×) · `gradp = sp.Matrix([
  sp.diff(phi, v) for v in xp])` (∇′φ) · `print(sp.simplify(curlp(phi*W) - phi*curlp(W) - gradp.cross(W)))` (a: zeros) ·
  `print(sp.simplify(gradp - d/r**3))` (b: zeros) · `print(sp.simplify(-gradp.cross(W) - W.cross(d)/r**3))` (c: zeros).
- **What it means.** The working form of "vorticity causes velocity": integrate over where the vorticity is. For thin
  tubes it becomes the filament law (5.17) (D12); for 2-D flows the planar kernel of point vortices (C12).
- **Traps.** ∇′(1/r) with the wrong sign (differentiating with respect to x). Applying Gauss when the field point lies
  inside V′ (the kernel is singular there — the smoothed kernel ε handles it numerically). Copying the book's two sign
  slips separately.

### D12 · The filament law: (5.16) → (5.17) — ★★, 6 steps, in C08 (notebook · `biot_savart_filament`)
- **Goal.** Reduce Biot–Savart to a rule for a thin vortex tube seen from far away: the velocity made by one short piece
  dl of it.
- **Start.** $\mathbf u(\mathbf x,t)=\frac1{4\pi}\int_{V'}\frac{\boldsymbol\omega(\mathbf x',t)\times(\mathbf x-\mathbf x')}{\lvert\mathbf x-\mathbf x'\rvert^3}d^3x'$ (5.16) —
  *in words:* every bit of vorticity pushes the fluid at x.
- **Plan.** (1) Take V′ as a short piece of a thin tube. (2) Freeze the kernel across the tube's cross-section.
  (3) Recognise the tube strength.
- **Tools.** Tube strength Γ = ∫ω·n dA (N03) · far-field ("frozen kernel") approximation (gloss, with a number).
- **Assumptions.** Thin tube: core radius a ≪ \|x − x′\| (step 3); e_ω constant across the section (step 2).
- **Steps.**
  1. *did:* Take a short piece of thin tube · *tex:* $V'=\Delta A'\,dl,\qquad d^3x'=d^2x'\,dl$ · *why:* A volume element of a tube
     is its cross-section times a length along it; we look at one piece at a time and add pieces later. · *plain:*
     Slice the vortex into short discs. · *set:* {shape: 'segment', M: 8}.
  2. *did:* Split ω into size and direction · *tex:* $\boldsymbol\omega=\lvert\boldsymbol\omega\rvert\,\mathbf e_\omega$ · *why:* In a thin tube all vortex
     lines are nearly parallel, so the unit vector e_ω is the same across the section. · *plain:* The vorticity points
     along the tube.
  3. *did:* Freeze the kernel across the section · *tex:* $\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}\approx\text{const over }\Delta A'$ · *why:* Far
     field (gloss): moving x′ across a core of radius a changes the kernel's size 1/\|x − x′\|² by about 2a/\|x − x′\| — 20 % at 10 core
     radii, 2 % at 100. This is the book's "sufficiently distant". · *plain:* From far away the tube looks like a line.
  4. *did:* Pull kernel and direction out · *tex:* $d\mathbf u\cong\frac{1}{4\pi}\Big[\int_{\Delta A'}\lvert\boldsymbol\omega\rvert d^2x'\Big]\mathbf e_\omega\times\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}dl$ ·
     *why:* Constants leave the integral (linearity). ≅ marks the far-field approximation. · *plain:* Only the total
     vorticity of the slice matters.
  5. *did:* Recognise the tube strength · *tex:* $\int_{\Delta A'}\lvert\boldsymbol\omega\rvert d^2x'=\Gamma$ · *why:* The vorticity flux through a
     section is the tube's strength (N03, Stokes), the same everywhere along it (D01). · *plain:* The bracket is Γ.
  6. *did:* Write the filament law · *tex:* $d\mathbf u=\frac{\Gamma\,dl}{4\pi}\,\mathbf e_\omega\times\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ · *why:*
     Substitute step 5 into step 4. This is (5.17). · *plain:* Each piece of filament pushes perpendicular to itself and
     to the line to x, as Γ dl/distance². · *set:* {M: 64}.
- **Result.** $d\mathbf u(\mathbf x,t)\cong\frac{\Gamma\,dl}{4\pi}\mathbf e_\omega\times\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ (5.17) — *in words:* the
  Biot–Savart law for a thin filament of strength Γ.
- **Check.** Units: (m²/s · m)/m² = m/s ✓. Summed round a ring of radius R at its centre: ΓR·2π/(4πR²) = Γ/2R ✓
  (`ring_axis_velocity(0, 0.5, 1)` = 1.0 m/s). A polygon of M segments converges to it with error ∝ 1/M² (1.0548, 1.0131,
  1.0008 for M = 8, 16, 64) ✓.
- **What it means.** The tool for rings (C13), horseshoe vortices and induced drag (Ch. 14). On the filament itself the
  sum diverges logarithmically: a curved filament needs a finite core to move itself (the ring self-speed, N42).
- **Traps.** Using (5.17) at points within a few core radii. Forgetting that dl runs along e_ω (the sense of Γ).

### D13 · A straight segment gives (Γ/4πd)(cos θ_a − cos θ_b); an infinite line gives Γ/2πd — ★★, 8 steps, in C08 (notebook · `biot_savart_filament`)
- **Goal.** Integrate the filament law along a straight segment in closed form, and close the loop the book leaves
  open: an infinite straight filament must give back the line vortex $u_\theta=\Gamma/2\pi r$ (5.2).
- **Start.** $d\mathbf u=\frac{\Gamma\,dl}{4\pi}\mathbf e_\omega\times\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ (5.17) — *in words:* each piece
  pushes as Γ dl/distance².
- **Plan.** (1) Place the segment on an axis and name the geometry. (2) Show every element pushes the same way.
  (3) Substitute an angle. (4) Integrate and take limits.
- **Tools.** Substitution in an integral (ch03 P106, reminder) · 1 + cot²θ = csc²θ (gloss) · improper integral as a limit
  (primer in C08) · right-hand rule (ch02 P74, reminder).
- **Assumptions.** Straight filament along e_z (steps 1–2); thin far-field form (5.17) (all steps).
- **Steps.**
  1. *did:* Place the filament on the z-axis · *tex:* $\mathbf x-\mathbf x'=d\,\mathbf e_d-l\,\mathbf e_z$ · *why:* e_ω = e_z; the field point is at
     distance d from the axis in direction e_d; l is the element's position along the axis measured from the foot of the
     perpendicular. · *plain:* The point sits a distance d from the line. · *set:* {shape: 'segment', ell: 1}.
  2. *did:* Compute the cross product · *tex:* $\mathbf e_z\times(d\,\mathbf e_d-l\,\mathbf e_z)=d\,\mathbf e_\varphi$ · *why:* e_z × e_z = 0 and
     e_z × e_d = e_φ (right-hand rule). Every element pushes along the same e_φ, so we can add magnitudes. · *plain:*
     All pieces push the same way round the line.
  3. *did:* Write the distance · *tex:* $\lvert\mathbf x-\mathbf x'\rvert^2=d^2+l^2$ · *why:* Pythagoras: e_d ⟂ e_z. · *plain:* Near
     elements are close, far ones far.
  4. *did:* Collect the scalar integral · *tex:* $u_\varphi=\frac{\Gamma d}{4\pi}\int_{l_a}^{l_b}\frac{dl}{(d^2+l^2)^{3/2}}$ · *why:* Steps 2–3 in
     (5.17); sum from the segment's start l_a to its end l_b. · *plain:* The speed is a sum of 1/distance³-weighted
     pieces.
  5. *did:* Substitute an angle · *tex:* $l=-d\cot\theta,\qquad dl=d\,\csc^2\theta\,d\theta$ · *why:* θ is the angle at the element
     between e_ω and the line to x (cos θ = −l/√(d² + l²)); substitution (P106) turns the power into a trig function. ·
     *plain:* Measure each element by the angle at which it sees the point.
  6. *did:* Simplify with 1 + cot² = csc² · *tex:* $\frac{d\,dl}{(d^2+l^2)^{3/2}}=\frac{d\cdot d\csc^2\theta\,d\theta}{d^3\csc^3\theta}=\frac{\sin\theta\,d\theta}{d}$ ·
     *why:* d² + l² = d²(1 + cot²θ) = d²csc²θ (gloss). · *plain:* Each angular slice contributes sin θ dθ/d.
  7. *did:* Integrate between the end angles · *tex:* $u_\varphi=\frac{\Gamma}{4\pi d}\int_{\theta_a}^{\theta_b}\sin\theta\,d\theta=\frac{\Gamma}{4\pi d}(\cos\theta_a-\cos\theta_b)$ ·
     *why:* The antiderivative of sin θ is −cos θ. θ_a, θ_b are the angles at the two ends. · *plain:* A finite
     segment's push depends only on the angles its ends subtend. · *live:* "d = 1 m, ends at ±1 m: (1/4π)(0.707 +
     0.707) = 0.1125 m/s".
  8. *did:* Let the segment become infinite · *tex:* $\theta_a\to0,\ \theta_b\to\pi:\ \ u_\varphi=\frac{\Gamma}{4\pi d}\big(1-(-1)\big)=\frac{\Gamma}{2\pi d}$ ·
     *why:* The improper integral is the limit of finite ones (primer). This is (5.2), $u_\theta=\Gamma/2\pi r$ — the book
     never closes this loop. A semi-infinite line (θ_b = π/2) gives half. · *plain:* An endless straight filament is
     the line vortex. · *set:* {ell: 20}.
- **Result.** $u_\varphi=\frac{\Gamma}{4\pi d}(\cos\theta_a-\cos\theta_b)$; infinite line $\frac{\Gamma}{2\pi d}$ = (5.2) — *in words:* summing the
  filament law along a straight line gives back the ideal line vortex.
- **Check.** Units m/s ✓. Γ = 1 m²/s, d = 1 m: infinite 0.159155, semi-infinite 0.079577 m/s; quad of the step-4
  integral agrees to 1e-12 ✓. Square loop of side 1 m at its centre: four segments, each seen at θ_a = 45°, θ_b = 135° from
  d = 0.5 m: 4 × Γ/(4π × 0.5) × (√2/2 + √2/2) = 2√2Γ/πL = 0.9003 m/s ✓.
- **What it means.** The law behind horseshoe vortices and downwash (Ch. 14): the velocity from a finite piece of a
  vortex is known exactly. Consistency with (5.2) confirms the + sign of (5.16).
- **Traps.** Measuring θ at the field point instead of at the element (then the signs of the cosines flip). Forgetting
  that the direction is set by e_ω × e_d.

### D14 · Navier–Stokes in a rotating frame in Lamb form: (5.20) → (5.21)–(5.24) → (5.25) — ★★, 8 steps, in C09 (notebook · `vorticity_equation_rotating`)
- **Goal.** Rewrite the rotating-frame momentum equation so that the vorticity appears explicitly, ready for the curl
  in D15.
- **Start.** $u_{i,i}=0$, $\frac{\partial u_i}{\partial t}+u_ju_{i,j}+2\varepsilon_{ijk}\Omega_ju_k=-\frac1\rho p_{,i}+g_i+\nu u_{i,jj}$ (5.19, 5.20) — *in words:*
  Boussinesq continuity and momentum in a frame turning at constant Ω, with effective gravity g (R12–R14).
- **Plan.** (1) Split the advective term into a gradient and a vorticity term. (2) Write the viscous term with ω.
  (3) Reorder the Coriolis indices. (4) Collect.
- **Tools.** Comma notation (gloss; Ch. 2 §2.14) · ε–δ identity (2.19) (Ch. 2 §2.7) · ω_k = ε_kmn u_{n,m} (R16) · dummy
  relabelling (Ch. 2 §2.1) · g = −∇Φ (4.18).
- **Assumptions.** ∇·u = 0 (step 5); Ω constant (step 6); ν constant (step 4); effective gravity conservative (step 7).
- **Steps.**
  1. *did:* Add and subtract u_j u_{j,i} · *tex:* $u_ju_{i,j}=u_j(u_{i,j}-u_{j,i})+u_ju_{j,i}$ · *why:* Adding zero; the bracket is
     antisymmetric, which is what vorticity is made of. · *plain:* Split advection into a rotation part and the rest.
  2. *did:* Express the bracket with ω · *tex:* $\varepsilon_{ijk}\omega_k=\varepsilon_{ijk}\varepsilon_{kmn}u_{n,m}=(\delta_{im}\delta_{jn}-\delta_{in}\delta_{jm})u_{n,m}=u_{j,i}-u_{i,j}$ ·
     *why:* ω_k = ε_kmn u_{n,m}, then ε–δ (2.19); ε_ijk ε_kmn = ε_kij ε_kmn by cyclic order. This is (5.22). · *plain:* The
     antisymmetric gradient is the vorticity in disguise.
  3. *did:* Substitute; write the rest as a gradient · *tex:* $u_ju_{i,j}=-u_j\varepsilon_{ijk}\omega_k+\tfrac12(u_ju_j)_{,i}=-(\mathbf u\times\boldsymbol\omega)_i+\tfrac12(u_j^2)_{,i}$ ·
     *why:* Step 2 gives the first term; u_j u_{j,i} = ½(u_ju_j)_{,i} (chain rule). This is (5.21), the index Lamb identity
     (R15). · *plain:* Advection = −u × ω + gradient of kinetic energy.
  4. *did:* Split the viscous term the same way · *tex:* $\nu u_{i,jj}=\nu(u_{i,j}-u_{j,i})_{,j}+\nu u_{j,ij}$ · *why:* Add and subtract
     νu_{j,ij}; Schwarz lets u_{j,ij} = u_{j,ji}. · *plain:* Diffusion = a rotation part + a divergence part.
  5. *did:* Use continuity and (5.22) · *tex:* $\nu u_{i,jj}=-\nu\varepsilon_{ijk}\omega_{k,j}$ · *why:* u_{j,ij} = (u_{j,j})_{,i} = 0 by (5.19); the
     bracket is −ε_ijkω_k by (5.22). This is (5.23): ν∇²u = −ν∇×ω (R17, (4.40)). · *plain:* Viscous force = −ν curl ω.
  6. *did:* Swap the Coriolis dummy indices · *tex:* $2\varepsilon_{ijk}\Omega_ju_k=-2\varepsilon_{ijk}\Omega_ku_j$ · *why:* Rename j ↔ k (dummies,
     Ch. 2 §2.1): 2ε_ikjΩ_ku_j, and ε_ikj = −ε_ijk. This is (5.24): Ω × u = −u × Ω. · *plain:* Put u first, like
     in the Lamb term.
  7. *did:* Write gravity with its potential · *tex:* $g_i=-\Phi_{,i}$ · *why:* (4.18); the effective gravity (true gravity plus
     centrifugal) is conservative (ch04 D18). · *plain:* Gravity is downhill in a potential.
  8. *did:* Substitute and collect the u × brackets · *tex:* $\frac{\partial u_i}{\partial t}+\Big(\tfrac12u_j^2+\Phi\Big)_{,i}-\varepsilon_{ijk}u_j(\omega_k+2\Omega_k)=-\frac1\rho p_{,i}-\nu\varepsilon_{ijk}\omega_{k,j}$ ·
     *why:* Steps 3, 5–7 in (5.20); −ε_ijku_jω_k − 2ε_ijku_jΩ_k = −ε_ijku_j(ω_k + 2Ω_k). This is (5.25). · *plain:* Momentum
     in the rotating frame, with the absolute vorticity ω + 2Ω inside the Lamb term. · *set:* {mode: 'budget'}.
- **Result.** $\partial u_i/\partial t+(\tfrac12u_j^2+\Phi)_{,i}-\varepsilon_{ijk}u_j(\omega_k+2\Omega_k)=-(1/\rho)p_{,i}-\nu\varepsilon_{ijk}\omega_{k,j}$ (5.25) —
  *in words:* Navier–Stokes in Lamb form, where the frame's rotation simply adds 2Ω to the vorticity.
- **Check.** Units m/s² ✓. Ω = 0 gives the Lamb form of Ch. 4 ✓. `rotating_lamb_form_terms` residual equals the (5.20)
  residual `rotating_ns_residual` to 1e-10 on random smooth fields ✓. Relabel check (5.24) with random numbers ✓.
- **What it means.** In a rotating frame the Coriolis force enters exactly like extra vorticity 2Ω — the planetary
  vorticity (R18) that D15 and C11 build on.
- **Traps.** Forgetting that step 5 needs ∇·u = 0. Losing the sign in step 6 (relabelling flips ε).

### D15 · The vorticity equation in a rotating frame with baroclinic generation: curl of (5.25) → (5.26)–(5.29) → (5.30) — ★★★, 15 steps, in C09 (notebook · `vorticity_equation_rotating`)
- **Goal.** Take the curl of (5.25) in index notation and arrive at the full vorticity equation (5.30), with every
  ε–δ move shown, including the term the book drops silently.
- **Start.** $\partial u_i/\partial t+(\tfrac12u_j^2+\Phi)_{,i}-\varepsilon_{ijk}u_j(\omega_k+2\Omega_k)=-(1/\rho)p_{,i}-\nu\varepsilon_{ijk}\omega_{k,j}$ (5.25)
  and $\omega_n=\varepsilon_{nqi}u_{i,q}$ — *in words:* rotating Navier–Stokes in Lamb form; the vorticity is the curl.
- **Plan.** (1) Apply ε_nqi( )_{,q} to every term, (5.26). (2) Treat the gradient, the Lamb, the pressure and the
  viscous terms one at a time with ε–δ, (5.27)–(5.29). (3) Collect and write in vector form.
- **Tools.** ε–δ identity (2.19) with cyclic reordering ε_nqi = ε_inq (ch02 P72, reminder) · symmetric × antisymmetric =
  0 (Ch. 2 §2.7) · ∇·ω = 0 (5.18) (N30) · product rule (ch01 P38) · quotient rule for ∇(1/ρ) (gloss) · comma
  notation (gloss).
- **Assumptions.** Boussinesq continuity u_{i,i} = 0 (step 9); Ω constant (step 7); ν constant (step 12); ρ kept inside
  the pressure term (step 10).
- **Steps.**
  1. *did:* Apply ε_nqi( )_{,q} to every term · *tex:* $\frac{\partial}{\partial t}(\varepsilon_{nqi}u_{i,q})+\varepsilon_{nqi}\big(\tfrac12u_j^2+\Phi\big)_{,iq}-\varepsilon_{nqi}\varepsilon_{ijk}\big[u_j(\omega_k+2\Omega_k)\big]_{,q}=-\varepsilon_{nqi}\big(\tfrac1\rho p_{,i}\big)_{,q}-\nu\varepsilon_{nqi}\varepsilon_{ijk}\omega_{k,jq}$ ·
     *why:* Since ω_n = ε_nqi u_{i,q}, this operation is the curl; applied to both sides it keeps the equation true.
     This is (5.26). · *plain:* Take the curl of every term, in index form. · *set:* {mode: 'budget', scene: 'column'}.
  2. *did:* Identify the local term · *tex:* $\frac{\partial}{\partial t}(\varepsilon_{nqi}u_{i,q})=\frac{\partial\omega_n}{\partial t}$ · *why:* ε is constant and
     derivatives commute (Schwarz); the bracket is ω_n by definition. · *plain:* The first term is the vorticity's
     local rate.
  3. *did:* Kill the gradient term · *tex:* $\varepsilon_{nqi}\big(\tfrac12u_j^2+\Phi\big)_{,iq}=0$ · *why:* ε_nqi is antisymmetric in q, i
     while the second derivative is symmetric in q, i; their contraction is zero (Ch. 2 §2.7). The book writes "Π" here
     — a slip for Φ. · *plain:* The curl of a gradient vanishes.
  4. *did:* Reorder ε and apply ε–δ · *tex:* $\varepsilon_{nqi}\varepsilon_{ijk}=\varepsilon_{inq}\varepsilon_{ijk}=\delta_{nj}\delta_{qk}-\delta_{nk}\delta_{qj}$ · *why:*
     Cyclic reordering (nqi → inq) does not change ε (P72); then (2.19) with the shared first index. · *plain:* Two ε's
     become four δ's.
  5. *did:* Contract the δ's in the Lamb term · *tex:* $-\varepsilon_{nqi}\varepsilon_{ijk}\big[u_j(\omega_k+2\Omega_k)\big]_{,q}=-\big[u_n(\omega_k+2\Omega_k)\big]_{,k}+\big[u_j(\omega_n+2\Omega_n)\big]_{,j}$ ·
     *why:* δ_njδ_qk sets j = n, q = k; δ_nkδ_qj sets k = n, q = j; the minus in front distributes. · *plain:* Two
     terms: one differentiates along ω + 2Ω, one along u.
  6. *did:* Expand the first bracket · *tex:* $-\big[u_n(\omega_k+2\Omega_k)\big]_{,k}=-u_{n,k}(\omega_k+2\Omega_k)-u_n(\omega_{k,k}+2\Omega_{k,k})$ ·
     *why:* Product rule (P38). · *plain:* Differentiate each factor in turn.
  7. *did:* Drop ∇·ω and ∇·Ω · *tex:* $-\big[u_n(\omega_k+2\Omega_k)\big]_{,k}=-u_{n,k}(\omega_k+2\Omega_k)$ · *why:* ω_{k,k} = 0 by (5.18) and Ω is
     constant, so Ω_{k,k} = 0. · *plain:* Only the stretching of the absolute vorticity remains.
  8. *did:* Expand the second bracket · *tex:* $\big[u_j(\omega_n+2\Omega_n)\big]_{,j}=u_{j,j}(\omega_n+2\Omega_n)+u_j\omega_{n,j}$ · *why:* Product
     rule again; Ω_{n,j} = 0. · *plain:* A divergence term and an advection term.
  9. *did:* Drop the divergence term the book skips · *tex:* $\big[u_j(\omega_n+2\Omega_n)\big]_{,j}=u_j\omega_{n,j}$ · *why:* u_{j,j} = 0 by (5.19).
     **The book's (5.27) drops u_{j,j}(ω_n + 2Ω_n) without a word**; this is why. Steps 5–9 are (5.27):
     $-u_{n,j}(\omega_j+2\Omega_j)+u_j\omega_{n,j}$. · *plain:* Incompressibility removes it.
  10. *did:* Differentiate the pressure term · *tex:* $-\varepsilon_{nqi}\Big(\frac1\rho p_{,i}\Big)_{,q}=-\frac1\rho\varepsilon_{nqi}p_{,iq}+\frac1{\rho^2}\varepsilon_{nqi}\rho_{,q}p_{,i}$ ·
      *why:* Product rule with (1/ρ)_{,q} = −ρ_{,q}/ρ² (quotient rule, gloss). ρ is **not** constant here. · *plain:*
      Pressure's curl splits into two pieces.
  11. *did:* Keep only the baroclinic piece · *tex:* $-\varepsilon_{nqi}\Big(\frac1\rho p_{,i}\Big)_{,q}=\frac1{\rho^2}\big[\nabla\rho\times\nabla p\big]_n$ · *why:*
      ε_nqi p_{,iq} = 0 (antisymmetric × symmetric); ε_nqi a_q b_i = (a × b)_n. This is (5.28) — D06's torque. · *plain:*
      Pressure twists fluid only where density gradients cross it.
  12. *did:* Apply ε–δ to the viscous term · *tex:* $-\nu\varepsilon_{nqi}\varepsilon_{ijk}\omega_{k,jq}=-\nu(\delta_{nj}\delta_{qk}-\delta_{nk}\delta_{qj})\omega_{k,jq}=-\nu\omega_{k,nk}+\nu\omega_{n,jj}$ ·
      *why:* Same reordering and identity as step 4. · *plain:* Two terms: a divergence and a Laplacian.
  13. *did:* Drop the divergence of ω · *tex:* $-\nu\varepsilon_{nqi}\varepsilon_{ijk}\omega_{k,jq}=\nu\omega_{n,jj}$ · *why:* ω_{k,nk} = (ω_{k,k})_{,n} = 0 by
      (5.18). This is (5.29). · *plain:* Viscosity diffuses vorticity.
  14. *did:* Assemble and move terms right · *tex:* $\frac{\partial\omega_n}{\partial t}+u_j\omega_{n,j}=u_{n,j}(\omega_j+2\Omega_j)+\frac1{\rho^2}[\nabla\rho\times\nabla p]_n+\nu\omega_{n,jj}$ ·
      *why:* Put steps 2, 3, 9, 11, 13 in (5.26) and move −u_{n,j}(ω_j + 2Ω_j) to the right (sign flips). · *plain:*
      Local + advective change = stretching + baroclinic + diffusion.
  15. *did:* Rename n → i; write vectors · *tex:* $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u+\frac1{\rho^2}\nabla\rho\times\nabla p+\nu\nabla^2\boldsymbol\omega$ ·
      *why:* ∂/∂t + u_j∂_j = D/Dt (3.5); a free index may be renamed. This is (5.30). With Ω = 0 and ∇ρ ∥ ∇p it is (5.13). ·
      *plain:* The full vorticity equation. · *set:* {scene: 'burgers'} · *watch:* "planetary and baroclinic bars drop
      to zero: (5.13) again".
- **Result.** $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u+\frac1{\rho^2}\nabla\rho\times\nabla p+\nu\nabla^2\boldsymbol\omega$ (5.30) —
  *in words:* following a particle (in the rotating frame), vorticity changes by stretching and tilting of the
  **absolute** vorticity, by baroclinic generation, and by viscous diffusion.
- **Check.** Units s⁻² every term ✓. Ω = 0, ρ = ρ(p) → (5.13) ✓ (asserted by `vorticity_budget_sym`). Numbers: a
  column at Ω = 7.292×10⁻⁵ rad/s stretched at ∂w/∂z = 10⁻⁵ s⁻¹ gains 2Ω∂w/∂z = 1.458×10⁻⁹ s⁻²; the lock exchange's
  baroclinic term 2.42 s⁻² (C04) ✓. **sympy check intent (★★★):** re-run steps 1–15 with `sp.LeviCivita` for a generic
  field: define u_i(x, y, z, t), p, ρ as generic functions and Ω as constant symbols; build each term of (5.26) by
  explicit sums over indices; show the gradient term is 0 (step 3); show the Lamb term equals
  −u_{n,j}(ω_j + 2Ω_j) + u_jω_{n,j} + u_{j,j}(ω_n + 2Ω_n) − u_n ω_{k,k} identically (steps 5–8, with the two dropped terms
  kept), then that the extra terms vanish for a divergence-free polynomial u (steps 7, 9); the pressure term equals
  [∇ρ × ∇p]_n/ρ² (step 11); the viscous term equals νω_{n,jj} for that u (step 13); finally (5.30) − curl(5.25) ≡ 0.
  `check_src` sketch: `X = sp.symbols('x0:3'); t = sp.symbols('t'); nu = sp.symbols('nu'); Om = sp.symbols('Omega0:3')` ·
  `u = [sp.Function(f'u{i}')(*X, t) for i in range(3)]` (generic velocity) · `d = lambda f, j: sp.diff(f, X[j])` (the comma)
  · `eps = sp.LeviCivita` · `w = [sum(eps(n, q, i)*d(u[i], q) for q in range(3) for i in range(3)) for n in range(3)]`
  (ω_n = ε_nqi u_{i,q}) · `lamb = [-sum(eps(n, q, i)*eps(i, j, k)*d(u[j]*(w[k] + 2*Om[k]), q) for q in range(3) for i in
  range(3) for j in range(3) for k in range(3)) for n in range(3)]` (step 1's third term) · `book = [sum(-d(u[n], j)*(w[j] +
  2*Om[j]) + u[j]*d(w[n], j) + d(u[j], j)*(w[n] + 2*Om[n]) - u[n]*d(w[j], j) for j in range(3)) for n in range(3)]` (steps
  5–8 with nothing dropped) · `print([sp.simplify(sp.expand(a - b)) for a, b in zip(lamb, book)])` (zeros) · then
  `ch05.vorticity_budget_sym(u_poly, X, t, nu, Om, rho_expr, p_expr)['residual_530']` for a divergence-free polynomial
  u_poly (zeros) and `['reduces_to_513']` (True).
- **What it means.** The backbone of geophysical vorticity dynamics: planetary vorticity 2Ω is stretched and tilted
  like any other (C11: potential vorticity), density contrasts across isobars create vorticity (fronts, sea breezes),
  and friction spreads it. u and ω are relative; the stretching acts on the absolute vorticity.
- **Traps.** Dropping u_{j,j}(ω_n + 2Ω_n) without saying why (step 9). Sign flips when moving terms right (step 14).
  Keeping ρ inside the pressure term while calling the fluid "Boussinesq" is a mixed approximation: to leading order
  ∇ρ × ∇p/ρ² ≈ ∇ρ′ × g/ρ₀ with g = −g e_z (cf. (4.86)).

### D16 · Natural coordinates on a vortex line: (ω·∇)u = ω ∂u/∂s, (5.31) — ★★, 6 steps, in C10 (notebook · `vorticity_stretching_tilting`)
- **Goal.** Rewrite the stretching–tilting term of (5.13) in coordinates that follow a vortex line, so each part has a
  picture.
- **Start.** $(\boldsymbol\omega\cdot\nabla)\mathbf u$ at a point where ω ≠ 0, with the local orthonormal frame (e_s, e_n, e_m) of
  Fig. 5.9 — *in words:* s along the vortex line, n away from its centre of curvature, m the second normal.
- **Plan.** (1) Write ∇ in the local frame. (2) Dot with ω. (3) Use that ω points along e_s. (4) Project on the frame.
- **Tools.** Frenet frame of a curve (primer in C10) · level sets and the directional derivative (ch02 P75, reminder) ·
  orthonormal basis and projection (ch02 P65, reminder).
- **Assumptions.** ω ≠ 0 at the point (step 3); the frame is frozen at the instant (step 5).
- **Steps.**
  1. *did:* Write ∇ in the line's own frame · *tex:* $\nabla=\mathbf e_s\frac{\partial}{\partial s}+\mathbf e_n\frac{\partial}{\partial n}+\mathbf e_m\frac{\partial}{\partial m}$ ·
     *why:* At one point any orthonormal basis can carry the gradient: each ∂/∂(·) is a directional derivative along its
     unit vector (P75). · *plain:* Measure changes along, across and out of the line. · *set:* {flow: 'shear'}.
  2. *did:* Dot the operator with ω · *tex:* $\boldsymbol\omega\cdot\nabla=(\boldsymbol\omega\cdot\mathbf e_s)\frac{\partial}{\partial s}+(\boldsymbol\omega\cdot\mathbf e_n)\frac{\partial}{\partial n}+(\boldsymbol\omega\cdot\mathbf e_m)\frac{\partial}{\partial m}$ ·
     *why:* Linearity of the dot product. · *plain:* Weight each directional derivative by ω's component along it.
  3. *did:* Use that ω points along the line · *tex:* $\boldsymbol\omega\cdot\mathbf e_s=\omega,\qquad\boldsymbol\omega\cdot\mathbf e_n=\boldsymbol\omega\cdot\mathbf e_m=0$ · *why:*
     A vortex line is tangent to ω, so e_s = ω/\|ω\| and ω = \|ω\| = the magnitude. · *plain:* Only the along-line
     derivative survives.
  4. *did:* Apply the operator to u · *tex:* $(\boldsymbol\omega\cdot\nabla)\mathbf u=\omega\frac{\partial\mathbf u}{\partial s}$ · *why:* Steps 2–3. This is (5.31)
     (the printed line ends in ω ∂/∂s; u is understood). · *plain:* Stretching–tilting = vorticity magnitude × how u
     changes along the line.
  5. *did:* Project on the three directions · *tex:* $\omega\frac{\partial\mathbf u}{\partial s}\ \to\ \Big(\omega\frac{\partial u_s}{\partial s},\ \omega\frac{\partial u_n}{\partial s},\ \omega\frac{\partial u_m}{\partial s}\Big)$ ·
     *why:* Components are projections on the unit vectors (P65), with the frame frozen at this point: ∂u_s/∂s means
     e_s·∂u/∂s. · *plain:* One along the line, two across it.
  6. *did:* Give each component its picture · *tex:* $\frac{\partial u_s}{\partial s}=\mathbf e_s\cdot\mathsf G\,\mathbf e_s$ · *why:* The along-line
     part is the stretching rate of a material element along e_s (Ch. 3 §3.4); ∂u_n/∂s and ∂u_m/∂s turn the element
     about m and n. · *plain:* Along = stretching; across = tilting. · *watch:* "purple arrow along ω, blue across".
- **Result.** $(\boldsymbol\omega\cdot\nabla)\mathbf u=\omega\frac{\partial\mathbf u}{\partial s}$ (5.31) — *in words:* the stretching–tilting term is
  the vorticity magnitude times the rate at which the velocity changes along the vortex line.
- **Check.** Units s⁻² ✓. `stretching_tilting_split`: stretching + tilting = Gω to 1e-12; stretching ∥ ω, tilting ⟂ ω;
  axial strain (α = 1 s⁻¹) along ω: rate = α, tilting 0 ✓.
- **What it means.** Only the variation of u along the line matters: parallel to ω it stretches, across ω it tilts.
- **Traps.** Thinking ∂u_n/∂s involves the frame's turning (the frame is frozen). Taking e_n toward the centre of
  curvature (the book's e_n points away — opposite to the usual Frenet N).

### D17 · Stretching and tilting components (5.32), and the angular-momentum reading — ★, 5 steps, in C10 (notebook · `vorticity_stretching_tilting`)
- **Goal.** Write (5.30) along the natural frame for an ideal fluid and see stretching as the spinning skater and
  tilting as turning — and why neither exists in 2-D.
- **Start.** (5.30) $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u+\frac1{\rho^2}\nabla\rho\times\nabla p+\nu\nabla^2\boldsymbol\omega$ —
  *in words:* the full vorticity equation.
- **Plan.** (1) Remove rotation, baroclinicity and viscosity. (2) Use (5.31). (3) Project. (4) Interpret.
- **Tools.** (5.31) (D16) · projection (ch02 P65, reminder) · angular momentum of a spinning cylinder (primer in C10).
- **Assumptions.** Inertial frame Ω = 0, barotropic, inviscid (step 1); frame frozen at the instant (step 3).
- **Steps.**
  1. *did:* Remove rotation, baroclinicity, viscosity · *tex:* $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u$ · *why:* Ω = 0, ∇ρ ∥ ∇p
     and ν = 0 remove three terms of (5.30). · *plain:* In an ideal fluid only stretching and tilting act. · *set:*
     {flow: 'axial', nu: 0}.
  2. *did:* Use the natural-coordinate form · *tex:* $\frac{D\boldsymbol\omega}{Dt}=\omega\frac{\partial\mathbf u}{\partial s}$ · *why:* (5.31) from D16. ·
     *plain:* The rate is set by how u changes along the line.
  3. *did:* Project on s, n, m · *tex:* $\frac{D\omega_s}{Dt}=\omega\frac{\partial u_s}{\partial s},\quad\frac{D\omega_n}{Dt}=\omega\frac{\partial u_n}{\partial s},\quad\frac{D\omega_m}{Dt}=\omega\frac{\partial u_m}{\partial s}$ ·
     *why:* Components in the frame frozen at the instant (the frame's own turning is not part of this instantaneous
     statement). This is (5.32). · *plain:* Stretching changes the along-line vorticity; tilting creates across-line
     vorticity. · *live:* "axial strain α = 1 s⁻¹: \|ω\| grows as e^{αt} = 7.39 at t = 2 s".
  4. *did:* Read stretching as angular momentum · *tex:* $L=I\,\Omega_{spin},\quad I\propto mA\ \Rightarrow\ \Omega_{spin}\propto1/A$ · *why:* A
     spinning cylinder of fixed mass m keeps L without torque; stretching it thins it (A falls), so it spins faster
     (primer). Pressure exerts no torque on it (D09). · *plain:* The skater pulling in her arms. · *set:* {flow:
     'axial'}.
  5. *did:* Check the 2-D case · *tex:* $\boldsymbol\omega=\omega\,\mathbf e_z,\ \mathbf u=(u,v,0)(x,y)\ \Rightarrow\ \frac{\partial\mathbf u}{\partial s}=\frac{\partial\mathbf u}{\partial z}=\mathbf 0$ ·
     *why:* In a planar flow ω is normal to the plane and nothing varies along it. · *plain:* In 2-D vorticity can be
     carried and diffused, never stretched or tilted. · *set:* {flow: 'planar'}.
- **Result.** $\frac{D\omega_s}{Dt}=\omega\frac{\partial u_s}{\partial s}$, $\frac{D\omega_n}{Dt}=\omega\frac{\partial u_n}{\partial s}$, $\frac{D\omega_m}{Dt}=\omega\frac{\partial u_m}{\partial s}$
  (5.32) — *in words:* stretching along a vortex line spins it up; shear across it tilts it; neither happens in 2-D.
- **Check.** Axial strain: ω_z(t) = ω₀e^{αt} (`uniform_strain_vorticity`), e² = 7.389 at αt = 2 ✓. Shear tilt w = s x:
  ω_z = s t ω_x0 (= 2.0 at s = 1 s⁻¹, t = 2 s) ✓. Planar strain with ω ⟂ plane: unchanged ✓.
- **What it means.** Stretching is how 3-D turbulence concentrates vorticity into thin intense tubes (Ch. 12); its
  absence in 2-D is why 2-D and 3-D turbulence differ; with planetary vorticity it spins up converging columns (C11).
- **Traps.** Reading (5.32) as equations for a fixed frame over finite time (the frame turns). Expecting stretching
  in a flat flow.

### D18 · A stretched tube spins faster (ω ∝ L), and Burgers' vortex balances stretching against diffusion — ★★, 9 steps, in C10 (notebook · `vorticity_stretching_tilting`)
- **Goal.** Put numbers on stretching: an ideal tube stretched to twice its length has twice the vorticity; with
  viscosity, a steady axial strain holds a vortex at a fixed core size — the book leaves this to Exercise 5.12 and
  never writes the balance.
- **Start.** Helmholtz 4 (tube strength constant in time, D08) and $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$ (5.13)
  — *in words:* an ideal tube keeps its strength; a real one also diffuses.
- **Plan.** (1) Ideal tube: strength and volume fixed, eliminate the area. (2) Burgers: write (5.13) for an axial strain
  plus swirl, integrate the steady balance once, then again.
- **Tools.** Helmholtz 4 (C05) · incompressibility A L = const (gloss) · cylindrical Laplacian of an axisymmetric scalar
  (gloss) · separation of variables (ch01 P42, reminder) · exponential solution of a first-order ODE (ch01 P36,
  reminder).
- **Assumptions.** Inviscid, uniform ω in a thin tube (steps 1–3); steady, axisymmetric, ν > 0 (steps 4–9); decay at
  R → ∞ and regularity at R = 0 (step 8).
- **Steps.**
  1. *did:* Fix the ideal tube's strength · *tex:* $\Gamma=\omega A=\text{const}$ · *why:* Helmholtz 4 (D08): an inviscid tube keeps its
     strength; for a thin tube with uniform ω the strength is ω times the section. · *plain:* Vorticity × area stays
     put.
  2. *did:* Fix its volume · *tex:* $A\,L=V=\text{const}$ · *why:* The tube is material (Helmholtz 1) and the fluid incompressible
     (gloss). · *plain:* Stretch it longer and it gets thinner.
  3. *did:* Eliminate the area · *tex:* $\omega=\frac{\Gamma}{A}=\frac{\Gamma L}{V}\ \propto\ L$ · *why:* Divide step 1 by A and use step 2. ·
     *plain:* Double the length, double the vorticity. · *live:* "ω₀ = 10 s⁻¹, L 1 → 2 m: ω = 20 s⁻¹, A halves".
  4. *did:* Take Burgers' axisymmetric strain · *tex:* $u_R=-\tfrac12\alpha R,\quad u_z=\alpha z,\quad\boldsymbol\omega=\omega_z(R)\,\mathbf e_z$ · *why:* A steady
     axial stretching flow (α > 0) that satisfies ∇·u = (1/R)∂(Ru_R)/∂R + ∂u_z/∂z = −α + α = 0; a swirl adds ω_z. ·
     *plain:* Fluid is drawn in sideways and stretched upward. · *set:* {flow: 'burgers', rate: 1}.
  5. *did:* Write the steady z-component of (5.13) · *tex:* $u_R\frac{d\omega_z}{dR}=\omega_z\frac{\partial u_z}{\partial z}+\nu\frac1R\frac{d}{dR}\Big(R\frac{d\omega_z}{dR}\Big)$ · *why:*
     Steady: Dω_z/Dt = u_R dω_z/dR; stretching ω_z∂u_z/∂z; cylindrical Laplacian of ω_z(R) (gloss). · *plain:* Inward
     carrying = stretching + diffusion.
  6. *did:* Insert u_R and ∂u_z/∂z = α · *tex:* $-\frac{\alpha R}{2}\omega_z'=\alpha\,\omega_z+\frac{\nu}{R}\big(R\,\omega_z'\big)'$ · *why:* Substitute
     step 4. · *plain:* All three processes, one ODE.
  7. *did:* Spot a total derivative · *tex:* $\frac1R\frac{d}{dR}\Big[\nu R\,\omega_z'+\frac{\alpha}{2}R^2\omega_z\Big]=0$ · *why:* (α/2R)(R²ω_z)′ = αω_z +
     (αR/2)ω_z′ (product rule), which is exactly the advection + stretching terms moved left. · *plain:* The balance is
     one derivative set to zero.
  8. *did:* Integrate once · *tex:* $\nu R\,\omega_z'+\frac{\alpha}{2}R^2\omega_z=0$ · *why:* The bracket is constant; at R = 0 it is 0
     (regularity), so it is 0 everywhere. · *plain:* Outward diffusive flux of vorticity = inward advective flux.
  9. *did:* Integrate again and normalise · *tex:* $\omega_z=\frac{\alpha\Gamma}{4\pi\nu}\,e^{-\alpha R^2/4\nu}$ · *why:* ω_z′/ω_z = −αR/2ν
     (separation, P42) gives a Gaussian; ∫ω_z 2πR dR = Γ fixes the amplitude to αΓ/4πν. · *plain:* A Gaussian core of
     radius √(4ν/α), where stretching and diffusion balance. · *live:* "α = 1 s⁻¹, ν = 10⁻⁶: core 2.0 mm, peak
     79.6 s⁻¹ for Γ = 10⁻³ m²/s".
- **Result.** Ideal tube: $\omega=\Gamma L/V\propto L$; Burgers: $\omega_z=\frac{\alpha\Gamma}{4\pi\nu}e^{-\alpha R^2/4\nu}$, core radius
  $\sqrt{4\nu/\alpha}$ — *in words:* stretching concentrates vorticity; viscosity spreads it; a steady strain holds the two in
  balance at a fixed core size.
- **Check.** Units: αΓ/ν in s⁻¹ ✓; √(ν/α) in m ✓. α → 0: the core widens without bound (pure diffusion, Lamb–Oseen) ✓.
  Term bars at R = 1 mm: advective + stretching = diffusion, residual ≈ 10⁻⁶ of the largest ✓ (`vorticity_terms`).
  The book's α is twice Wikipedia's (same flow; Wikipedia writes v_r = −αr, v_z = 2αz) ✓.
- **What it means.** Burgers' vortex is the classic model of the thin intense vortices of turbulence (Ch. 12); the core
  size is set by the strain and the viscosity, not by how the vortex started.
- **Traps.** Expecting ωA constant with viscosity. Using Wikipedia's α in the book's formula (a factor 2).

### D19 · Kelvin's theorem in a rotating frame: the absolute circulation (5.33) — ★★, 10 steps, in C11 (notebook · `vorticity_equation_rotating`)
- **Goal.** Redo Kelvin's proof in a frame turning at constant Ω, where the Coriolis force is not conservative, and
  find the quantity that is still conserved. The book states (5.33) and refers to "Exercise 5.11" (meaning 5.10).
- **Start.** $\frac{D\mathbf u}{Dt}=-\frac1\rho\nabla p-\nabla\Phi-2\boldsymbol\Omega\times\mathbf u$ ((5.20) with ν = 0) and $\Gamma=\oint_C\mathbf u\cdot d\mathbf x$ for a
  material loop — *in words:* inviscid momentum in the rotating frame; u is the velocity measured in that frame.
- **Plan.** (1) Repeat D04–D05 with the extra Coriolis term. (2) Show the Coriolis loop integral is −2Ω·(rate of the
  loop's vector area). (3) Move it to the left and use Stokes.
- **Tools.** D04–D05 moves · vector area of a closed loop (primer in C11) · triple product (Ω×u)·dx = Ω·(u×dx) (ch02
  §2.7) · closed-loop integral of an exact differential (primer in C03) · Stokes' theorem (2.34).
- **Assumptions.** Inviscid, barotropic, conservative effective gravity (step 3); Ω constant (step 8).
- **Steps.**
  1. *did:* Reuse the material-loop rate · *tex:* $\frac{D\Gamma}{Dt}=\oint_C\frac{D\mathbf u}{Dt}\cdot d\mathbf x$ · *why:* D04 did not depend on the
     frame: ∮u·du = 0 for any single-valued u. · *plain:* Only the acceleration along the loop counts. · *set:* {mode:
     'ring'}.
  2. *did:* Substitute rotating-frame momentum · *tex:* $\frac{D\Gamma}{Dt}=-\oint_C\frac{dp}{\rho}-\oint_Cd\Phi-2\oint_C(\boldsymbol\Omega\times\mathbf u)\cdot d\mathbf x$ · *why:*
     As in D05 steps 1–4, with the Coriolis force −2Ω × u added. · *plain:* A new, frame-made term appears.
  3. *did:* Drop the barotropic and conservative terms · *tex:* $\frac{D\Gamma}{Dt}=-2\oint_C(\boldsymbol\Omega\times\mathbf u)\cdot d\mathbf x$ · *why:* D05 steps
     5–6: exact differentials integrate to zero round a closed loop. · *plain:* Only the Coriolis term can change Γ.
  4. *did:* Reorder the triple product · *tex:* $(\boldsymbol\Omega\times\mathbf u)\cdot d\mathbf x=\boldsymbol\Omega\cdot(\mathbf u\times d\mathbf x)$ · *why:* Cyclic
     property of a·(b×c) (Ch. 2 §2.7); Ω is constant, so it leaves the integral. · *plain:* Pull Ω out.
  5. *did:* Define the loop's vector area · *tex:* $\mathbf A_{vec}=\tfrac12\oint_C\mathbf x\times d\mathbf x$ · *why:* For a planar loop this is area ×
     unit normal (primer); in general it is ∫n dA over any surface spanning C. · *plain:* A vector whose length is the
     loop's area.
  6. *did:* Differentiate it following the fluid · *tex:* $\frac{D\mathbf A_{vec}}{Dt}=\tfrac12\oint_C\big(\mathbf u\times d\mathbf x+\mathbf x\times d\mathbf u\big)$ · *why:*
     Product rule for a cross product (ch04 P125) with Dx/Dt = u and D(dx)/Dt = du (D04 step 6). · *plain:* The area
     changes as the loop's points move.
  7. *did:* Integrate the second term by parts · *tex:* $\oint_C\mathbf x\times d\mathbf u=-\oint_Cd\mathbf x\times\mathbf u=\oint_C\mathbf u\times d\mathbf x\ \Rightarrow\ \frac{D\mathbf A_{vec}}{Dt}=\oint_C\mathbf u\times d\mathbf x$ ·
     *why:* d(x × u) = dx × u + x × du, and ∮d(x × u) = 0 (exact differential, primer). · *plain:* The area grows at
     the loop integral of u × dx.
  8. *did:* Move the Coriolis term to the left · *tex:* $\frac{D}{Dt}\big(\Gamma+2\boldsymbol\Omega\cdot\mathbf A_{vec}\big)=0$ · *why:* Steps 3, 4 and 7 give
     DΓ/Dt = −2Ω·DA_vec/Dt; Ω is constant, so it passes inside D/Dt. · *plain:* Γ plus twice Ω dotted with the area is
     conserved.
  9. *did:* Write the planetary part as a flux · *tex:* $2\boldsymbol\Omega\cdot\mathbf A_{vec}=2\int_A\boldsymbol\Omega\cdot\mathbf n\,dA$ · *why:* A_vec = ∫n dA over any
     surface bounded by C (primer; Stokes for the uniform field Ω×x gives the same). · *plain:* It is the planetary
     vorticity 2Ω passing through the loop.
  10. *did:* Name the absolute circulation · *tex:* $\frac{D\Gamma_a}{Dt}=0,\qquad\Gamma_a\equiv\int_A(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\mathbf n\,dA=\Gamma+2\int_A\boldsymbol\Omega\cdot\mathbf n\,dA$ ·
      *why:* Γ = ∫ω·n dA (Stokes). This is (5.33); the book's pointer "Exercise 5.11" should read 5.10. · *plain:* The
      circulation of the absolute vorticity is conserved. · *live:* "ring 30° → 60° N: Γ = 2ΩA(sin 30° − sin 60°) =
      −4.19×10⁷ m²/s".
- **Result.** $\frac{D\Gamma_a}{Dt}=0$, $\Gamma_a\equiv\int_A(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\mathbf n\,dA=\Gamma+2\int_A\boldsymbol\Omega\cdot\mathbf n\,dA$ (5.33) —
  *in words:* in a rotating frame the relative circulation can change, but the absolute circulation (relative plus
  planetary vorticity through the loop) cannot.
- **Check.** Units m²/s ✓. Equivalent route: inertial velocity u + Ω × x and ∮(Ω × x)·dx = 2Ω·A_vec (Stokes) give the same
  line ✓. The rotating scenario of E3 (Ω = 0.5 rad/s, α = 0.2 s⁻¹): Γ(5 s) = 1.9859 m²/s while Γ_a = π stays fixed ✓.
  Fluid at rest in the inertial frame seen from the rotating frame: Γ = −2ΩA, Γ_a = 0 ✓.
- **What it means.** Bjerknes' circulation theorem: a ring of air moving poleward, or shrinking, must acquire
  anticyclonic relative circulation to keep Γ_a — the seed of potential vorticity (Ch. 13).
- **Traps.** Writing 2ΩA instead of 2Ω·A_vec (only equal for a horizontal loop and vertical Ω). Thinking the Coriolis
  force "does no work, so it cannot change Γ" — it does no work on a particle but its loop integral is not zero.

### D20 · A fluid column in a rotating layer keeps (ω_z + 2Ω)/h — ★★, 6 steps, in C11 (notebook · `vorticity_equation_rotating`)
- **Goal.** Explain Fig. 5.10: a column of fluid stretched or squashed in a rotating layer changes its spin in a
  predictable way. The book describes it in words only.
- **Start.** $\frac{D\Gamma_a}{Dt}=0$, $\Gamma_a=\int_A(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\mathbf n\,dA$ (5.33), applied to a thin horizontal material
  loop around a vertical column of height h — *in words:* the absolute circulation of the loop is fixed.
- **Plan.** (1) Approximate Γ_a for a thin loop. (2) Use mass conservation of the column. (3) Divide.
- **Tools.** (5.33) (D19) · column mass A h = const for an incompressible column (gloss) · f = 2Ω sin φ (ch04 P126,
  reminder).
- **Assumptions.** Inviscid, barotropic (step 2); thin loop, ω_z nearly uniform inside (step 1); the column stays
  vertical — Taylor–Proudman, Ch. 13 (step 3).
- **Steps.**
  1. *did:* Approximate Γ_a for a thin loop · *tex:* $\Gamma_a\approx(\omega_z+2\Omega)\,A$ · *why:* Over a small horizontal loop ω_z and
     the vertical planetary vorticity 2Ω (f = 2Ω sin φ on the Earth) are nearly uniform; n = e_z. · *plain:* Absolute
     vorticity times area. · *set:* {mode: 'column', bump: 0}.
  2. *did:* Hold it fixed in time · *tex:* $(\omega_z+2\Omega)\,A=\text{const}$ · *why:* (5.33): the loop is material. · *plain:*
     If the loop shrinks, the absolute vorticity must grow.
  3. *did:* Use the column's mass · *tex:* $\rho A h=\text{const}\ \Rightarrow\ A h=\text{const}$ · *why:* The column's particles stay in it
     (material) and its density is fixed (incompressible, gloss). · *plain:* A taller column is a thinner column.
  4. *did:* Divide the two constants · *tex:* $\frac{\omega_z+2\Omega}{h}=\text{const}$ · *why:* (ω_z + 2Ω)A ÷ Ah; A cancels. · *plain:*
     Absolute vorticity grows in proportion to the column's height.
  5. *did:* Solve for the new relative vorticity · *tex:* $\omega_z=(\omega_{z0}+2\Omega)\frac{h}{h_0}-2\Omega$ · *why:* Apply step 4 between the
     start (ω_z0, h₀) and now. · *plain:* Stretch a column at rest and it spins with the planet; squash it and it spins
     against. · *live:* "f = 10⁻⁴ s⁻¹, h 1000 → 1100 m: ω_z = 1.1 × 10⁻⁴ − 10⁻⁴ = 1.0×10⁻⁵ s⁻¹ (cyclonic)" · *set:*
     {bump: -100} (a 100 m trough: h = 1100 m).
  6. *did:* Check small changes against N39 · *tex:* $\frac{D\omega_z}{Dt}=\frac{\omega_z+2\Omega}{h}\frac{Dh}{Dt}\approx2\Omega\frac{\partial w}{\partial z}$ · *why:*
     Differentiate step 4; (1/h)Dh/Dt = ∂w/∂z for a uniformly stretched column, and ω_z ≪ 2Ω. This is the planetary
     stretching term Dω_z/Dt = 2Ω∂w/∂z (N39). · *plain:* The same physics, written as a rate.
- **Result.** $\frac{\omega_z+2\Omega}{h}=\text{const}$, i.e. $\omega_z=(\omega_{z0}+2\Omega)h/h_0-2\Omega$ — *in words:* a column stretched in a
  rotating layer spins up in the sense of the rotation (cyclonic); a squashed one spins down or reverses
  (anticyclonic).
- **Check.** Units s⁻¹/m on both sides ✓. h = h₀ gives ω_z = ω_z0 ✓. Ω = 0: ω_z/h = const — pure stretching of a vortex
  (D18) ✓. Number: +10 % height at f = 10⁻⁴ s⁻¹ → ζ = +1.0×10⁻⁵ s⁻¹ ✓.
- **What it means.** The first form of potential-vorticity conservation, (ζ + f)/h = const (Ch. 13): columns crossing a
  mountain ridge turn anticyclonic, air converging into a low spins up cyclonically, Taylor columns and Rossby waves
  follow.
- **Traps.** Using 2Ω instead of the local vertical component f = 2Ω sin φ on the Earth. Forgetting the column must stay
  vertical and inviscid.

### D21 · Two point vortices: the centre of vorticity, the orbit rate, and the translating pair — ★, 6 steps, in C12 (notebook · `point_vortex_lab`)
- **Goal.** Predict how two line vortices move each other: where the pair turns, how fast, and why an opposite pair
  marches off in a straight line. The book gives V₁ and V₂; the rest is left to Exercise 5.18.
- **Start.** $u_\theta=\Gamma/2\pi r$ (5.2) for each vortex, and Helmholtz 1: each vortex moves with the velocity the *other*
  induces (a straight line vortex does not move itself) — *in words:* vortices are carried by the flow.
- **Plan.** (1) Compute each vortex's velocity. (2) Show the separation is fixed, so the motion is a rigid rotation.
  (3) Find its centre and rate. (4) Take the opposite-pair limit.
- **Tools.** Line vortex (R03) · lever rule / centre of mass analogy (gloss) · circular motion speed = rate × radius
  (gloss).
- **Assumptions.** Ideal line vortices, no walls (all steps); Γ₁, Γ₂ counterclockwise positive (convention 3).
- **Steps.**
  1. *did:* Velocity of vortex 2 due to 1 · *tex:* $V_1=\frac{\Gamma_1}{2\pi h}$ · *why:* (5.2) at distance h, perpendicular to the joining
     line; for Γ₁ > 0 and 2 to the right of 1, it points up (Fig. 5.11). · *plain:* Vortex 1 pushes vortex 2 sideways.
     · *set:* {preset: 'equal'}.
  2. *did:* Velocity of vortex 1 due to 2 · *tex:* $V_2=\frac{\Gamma_2}{2\pi h}$ · *why:* The same law; vortex 1 is to the left of 2,
     so the push points down. · *plain:* Vortex 2 pushes vortex 1 the other way.
  3. *did:* Show the separation cannot change · *tex:* $\frac{d}{dt}\lvert\mathbf x_2-\mathbf x_1\rvert^2=2(\mathbf x_2-\mathbf x_1)\cdot(\mathbf v_2-\mathbf v_1)=0$ · *why:* Both
     velocities are perpendicular to x₂ − x₁. Fixed distance + perpendicular velocities = rigid rotation about a point
     G on the joining line. · *plain:* The pair turns like a rigid dumbbell.
  4. *did:* Locate the fixed point G · *tex:* $\frac{V_2}{h_1}=\frac{V_1}{h_2}\ \Rightarrow\ h_1=\frac{\Gamma_2\,h}{\Gamma_1+\Gamma_2}$ · *why:* Rigid rotation: speed = rate ×
     distance from the centre; with h₁ + h₂ = h this is the lever rule. G is the centre of vorticity
     (Γ₁x₁ + Γ₂x₂)/(Γ₁ + Γ₂), fixed since ΣΓ_kx_k is conserved. · *plain:* G sits nearer the
     stronger vortex. · *live:* "Γ₂ = 3Γ₁, h = 1 m: h₁ = 0.75 m".
  5. *did:* Find the rotation rate · *tex:* $\dot\theta=\frac{V_1+V_2}{h}=\frac{\Gamma_1+\Gamma_2}{2\pi h^2}$ · *why:* The relative speed of the
     two vortices, V₁ + V₂, divided by their distance is the angular rate of the joining line. · *plain:* Stronger or
     closer vortices orbit faster. · *live:* "Γ₁ = Γ₂ = 1 m²/s, h = 1 m: 0.318 rad/s, period 19.7 s".
  6. *did:* Take the opposite pair · *tex:* $\Gamma_2=-\Gamma_1:\ \ h_1\to\infty,\ \dot\theta=0,\ V=\frac{\Gamma}{2\pi h}$ · *why:* With Γ₁ + Γ₂ = 0
     the centre goes to infinity; both velocities are equal and parallel, so the pair translates without turning. ·
     *plain:* An opposite pair marches off together. · *set:* {preset: 'opposite'}.
- **Result.** Centre of vorticity at $h_1=\Gamma_2h/(\Gamma_1+\Gamma_2)$ from vortex 1, orbit rate $(\Gamma_1+\Gamma_2)/2\pi h^2$; an opposite
  pair translates at $\Gamma/2\pi h$ — *in words:* like vortices orbit their centre of vorticity; opposite ones travel
  together.
- **Check.** Units: m²/s/m² = rad/s ✓. Equal strengths: G at the midpoint, period 19.74 s for Γ = 1 m²/s, h = 1 m; the
  numerical integration (`point_vortex_evolve`) finds the same period and keeps ΣΓx, ΣΓ\|x\|² and H flat to 1e-10 ✓.
  ⚠️ **The caption of Fig. 5.11 says G is "where the combined velocity induced by the two vortices is zero" — true only
  for Γ₁ = Γ₂.** For Γ₂ = 3Γ₁ the fluid at G moves at 1/(2π·0.75) − 3/(2π·0.25) = −1.70 m/s; G is the fixed centre of the
  motion, not a stagnation point of the fluid.
- **What it means.** Binary cyclones orbit each other (Fujiwhara); aircraft wake vortex pairs descend together
  (Ch. 14); every point-vortex model rests on this.
- **Traps.** Letting a vortex push itself. Mistaking the centre of vorticity for a point of zero fluid velocity
  (the caption's slip). Putting G nearer the weaker vortex.

### D22 · A vortex near a wall: its image, u·n = 0 on the wall, and the drift Γ/4πh — ★, 5 steps, in C13 (notebook · `point_vortex_lab`)
- **Goal.** Show that a mirror vortex of opposite sign makes a plane wall a streamline, and find how fast the real vortex
  slides along the wall.
- **Start.** Vortex A of strength Γ at distance h from the wall y = 0, and the inviscid wall condition v = 0 at y = 0 —
  *in words:* no fluid may cross the wall.
- **Plan.** (1) Replace the wall by an image. (2) Check the wall condition by symmetry. (3) Compute A's velocity.
- **Tools.** Mirror symmetry of a vortex and its opposite image (gloss) · line vortex (R03) · superposition (gloss).
- **Assumptions.** Inviscid (slip allowed along the wall, step 3); infinite plane wall (step 1).
- **Steps.**
  1. *did:* Place an opposite image behind the wall · *tex:* $\text{B}:\ -\Gamma\ \text{at}\ (x_A,-h)$ · *why:* Method of images: in the
     region y > 0 the pair A + B is an ideal flow; if it satisfies the wall condition it is the answer (uniqueness). ·
     *plain:* Imagine a mirror twin spinning the other way. · *set:* {preset: 'wall'}.
  2. *did:* Compare the two pushes at a wall point · *tex:* $r_A=r_B\ \Rightarrow\ \lvert\mathbf V_A\rvert=\lvert\mathbf V_B\rvert=\frac{\Gamma}{2\pi r}$ · *why:*
     Every wall point P is equidistant from A and B (mirror images). · *plain:* Both vortices push equally hard at the
     wall.
  3. *did:* Show the normal parts cancel · *tex:* $v_A+v_B=0\ \text{on}\ y=0$ · *why:* Mirror symmetry: the opposite spin reflects the
     normal component and keeps the tangential one (gloss; Fig. 5.14's 90° construction). · *plain:* The wall is a
     streamline; fluid only slides along it. · *watch:* "the wall arrows lie flat along the wall".
  4. *did:* Give A the image's velocity · *tex:* $V_A=\frac{\Gamma}{2\pi(2h)}$ · *why:* A does not move itself (Helmholtz 1); only B, at
     distance 2h, moves it. · *plain:* The mirror twin drags the vortex.
  5. *did:* Simplify and read the direction · *tex:* $V_A=\frac{\Gamma}{4\pi h}\ \ \text{parallel to the wall}$ · *why:* 2π·2h = 4πh; the
     image's push at A is perpendicular to AB, i.e. along the wall (for Γ > 0 above the wall: +x). · *plain:* A vortex
     slides along a wall at Γ/4πh. · *live:* "Γ = 1 m²/s, h = 0.5 m: 0.159 m/s".
- **Result.** $V_A=\frac{\Gamma}{4\pi h}$ parallel to the wall — *in words:* a wall acts like an opposite vortex behind it; the
  vortex slides along the wall, faster the closer it is.
- **Check.** Units m/s ✓. `wall_image_system`: v = 0 at 200 wall points to 1e-12 ✓. Channel of height H: (Γ/4H)cot(πh/H)
  → Γ/4πh as h → 0 and 0 at h = H/2 ✓.
- **What it means.** Why the knife-blade pair separates along the bucket wall (Fig. 5.13), why a smoke ring widens near a
  wall (N42), and ground effect for wing-tip vortices (Ch. 14). Images enforce only no penetration — real walls also
  make boundary layers (Ch. 9).
- **Traps.** Using a same-sign image (the wall would be crossed). Using distance h instead of 2h.

### D23 · The strength of a vortex sheet is the jump in tangential velocity — ★, 5 steps, in C14 (notebook · `vortex_sheet_rollup`)
- **Goal.** Show that the circulation per unit length of a vortex sheet equals the jump in the velocity along it.
- **Start.** $\Gamma=\oint\mathbf u\cdot d\mathbf x$ (3.18) round a small rectangle ds × dn straddling the sheet (Fig. 5.16), counterclockwise;
  u₁ is the tangential velocity just above, u₂ just below, v the normal velocity — *in words:* measure the circulation of
  a thin box around the sheet.
- **Plan.** (1) Add the four sides with their directions. (2) Let the box become thin.
- **Tools.** Circulation round a rectangle (Ch. 2 §2.13) · right-hand rule and orientation (ch02 P74, reminder).
- **Assumptions.** v continuous across the sheet (step 2); the sheet thinner than dn (step 4).
- **Steps.**
  1. *did:* Walk the bottom side forward · *tex:* $+u_2\,ds$ · *why:* Counterclockwise, the lower side (below the sheet) is walked in
     +s, where the tangential velocity is u₂. · *plain:* Below the sheet we move with u₂. · *set:* {mode: 'jump'}.
  2. *did:* Walk up the right, down the left · *tex:* $v_R\,dn-v_L\,dn	o0$ · *why:* The normal velocity is continuous across the sheet, so each short side is
     just v·dn with no jump; together they are of size dn·ds and vanish as the box is flattened, dn → 0. · *plain:* The ends shrink away as the box flattens.
  3. *did:* Walk the top side backwards · *tex:* $-u_1\,ds$ · *why:* Counterclockwise, the upper side is walked in −s, where the
     velocity is u₁. · *plain:* Above the sheet we move against u₁.
  4. *did:* Add the sides; let dn → 0 · *tex:* $d\Gamma=u_2\,ds+v\,dn-u_1\,ds-v\,dn=(u_2-u_1)\,ds$ · *why:* Sum of steps 1–3; letting the
     box become thin keeps only what jumps across the sheet. This is the book's line. · *plain:* The box's circulation is
     the jump times its length.
  5. *did:* Divide by ds: the strength · *tex:* $\gamma\equiv\frac{d\Gamma}{ds}=u_2-u_1$ · *why:* Circulation per unit length; the book
     calls it Γ, we write γ [m/s] (convention 4). **Fig. 5.16's caption prints u₁ − u₂**: the clockwise sense (the
     magnitude for its clockwise filaments). · *plain:* A sheet's strength is its velocity jump. · *live:* "γ = 2 m/s:
     u₂ = +1, u₁ = −1 m/s" · *set:* {conv: true} · *watch:* "the caption convention flips the sign only".
- **Result.** $\gamma=\frac{d\Gamma}{ds}=u_2-u_1$ (counterclockwise circulation per unit length; u₁ above, u₂ below) — *in words:* a
  jump in tangential velocity *is* a sheet of vorticity, and its strength is the jump.
- **Check.** Units m/s ✓. A row of N filaments of strength γ ds gives u → ∓γ/2 just above and below (N44); the L1 error of
  u(y) against the continuous sheet falls as 1/N (0.0440, 0.0044, 0.00044 for N = 10, 100, 1000) ✓. Clockwise filaments
  (the book's figure): u₁ > 0 above, u₂ < 0 below, γ < 0 ✓.
- **What it means.** Every shear layer, mixing layer and wing wake is a (thickened) vortex sheet; it is unstable and
  rolls up (Kelvin–Helmholtz, Ch. 11); a wall's boundary layer is a sheet stuck to the wall (Ch. 9).
- **Traps.** Mixing the text's and the caption's sign conventions. Calling it Γ and confusing it with a circulation
  (units m/s, not m²/s).
