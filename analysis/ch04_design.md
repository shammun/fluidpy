# Chapter 4 — Conservation Laws: lesson design
(from `analysis/ch04_curation.md` (A 15 · B 151 · C 19 = 185 rows; CORE 15 · NOTE 156 · RECAP 12 · SKIP 2; 30 derivations
D01–D30 written out (★ 5 · ★★ 23 · ★★★ 2); E1–E9 + backup B1) and `analysis/ch04.md` (§2b derivations, §4 rows 1–70,
§5 core modules, §9 conventions and typos); the equations placed below were re-read on the rendered pages
`chapters/pages/ch04/p124–p173.png` (printed pp. 97–146: (4.1)–(4.4) p124, (4.5)–(4.8) p125, stream functions and
(4.12) p126–p127, (4.13)–(4.16d) p128, (4.17)–(4.18) p129, Ex. 4.2 and (4.19) p132–p133, rocket and (4.20a)–(4.22) p137,
(4.23)–(4.25) p138, (4.26)–(4.30) p139, (4.31)–(4.36) p140, (4.37) p141, (4.38)–(4.41) p142, (4.42) p143, (4.43)–(4.45)
p144, Coriolis projectile p145, highs/lows and Fig. 4.9 p146, (4.46)–(4.51) p149, (4.52)–(4.57) p150, (4.58)–(4.60)
p151, (4.66)–(4.71) p155, (4.73)–(4.78) p157, (4.84)–(4.85) p162, (4.86)–(4.87) p163, (4.88)–(4.89) p164, pillbox p165,
(4.90)–(4.94) p166, (4.100)–(4.104) p173); 2026-09-23, lesson-designer. The implementer works in parallel from analysis
§4 + curation §8; **Part C is the contract both sides keep.**)

**Binding conventions for every builder (analysis §9, curation decisions 5–6).**
1. **Fields are callables** exactly as in Ch. 3: `u(x, t) -> ndarray` with `x` of shape `(d,)` or `(d, N)` (components on
   axis 0), scalar fields `F(x, t)`; a density argument `rho` may be a float **or** a callable `rho(x, t)`; stress
   fields `tau(x, t) -> (3, 3)` or `(3, 3, N)`. Every public function is scalar-callable (scalar in → float out) so
   explainer parity rows can call it; `py:` parity expressions have **no builtins and no lambdas** — functions that take
   a field also accept a **preset name** plus keywords (the ch03 "kind name or callable" pattern), and results are
   indexed down to one number (`["result"]`, `[0]`, `.residual`).
2. **Index conventions.** Velocity gradient `G[i, j] = ∂u_i/∂x_j` (row = component). Stress `tau[i, j]`: **first index
   = face normal**, second = force direction; traction $f_j=n_i\tau_{ij}$ (2.15). ⚠️ **The stress divergence in
   Cauchy's equation contracts the FIRST index**, $(\nabla\cdot\boldsymbol\tau)_j=\partial\tau_{ij}/\partial x_i$
   ((4.20b), (4.24)) — on grids `core.operators.tensor_divergence(T, h, index=0)`; ch02's default is `index=1`. Every
   function below that differentiates a stress uses the first index, and a non-symmetric test τ pins it.
3. **Axes and gravity.** z (or x₃) **up**, $\mathbf g=-g\mathbf e_z$, $\Phi=gz$ so $\mathbf g=-\nabla\Phi$ (4.18). Chapter
   functions default to **g = 9.81 m/s²** (the book's value in its examples; ch01's `G0` = 9.80665 is the standard value —
   the difference is 0.03 % and no worked number changes at 3 figures). Worked examples use water μ = 1.0×10⁻³ Pa s,
   ρ = 1000 kg/m³, k = 0.6 W/(m K) exactly (easy numbers); figures that need real properties call
   `ch01.fluid_properties` / `ch01.FLUIDS`.
4. **Four meanings of the prime in one chapter** (⚠️ callout in the §4.1 notation table and again at C09, C11, C13):
   rotating (noninertial) frame x′, u′, D′/Dt (§4.7, C09); dummy integration variable p′ in
   $\int_{p_o}^{p}dp'/\rho(p')$ (4.67) (C11); perturbation from the hydrostatic state p′ = p − p_s, ρ′ = ρ − ρ_s (4.84)
   (C13); Ch. 3's translating frame. Code never uses a trailing `_prime` for two meanings: `u_rot`, `x_rot` (rotating
   frame), `p_pert`, `rho_pert` (perturbations).
5. **Coriolis sign language.** `core.rotating.frame_acceleration_terms` returns the **accelerations** of (4.43)/(4.44)
   with their + signs (`coriolis` = +2Ω × u′, `centripetal` = Ω × (Ω × x′)); `core.rotating.apparent_body_forces`,
   `coriolis_acceleration` and `centrifugal_acceleration` return the **apparent forces per unit mass** of the bracket in
   (4.45) with − signs (−2Ω × u′, −Ω × (Ω × x′)). The book itself calls both "Coriolis acceleration"; we say *Coriolis
   term* for +2Ω × u′ on the left of (4.44) and *Coriolis force per unit mass* for −2Ω × u′ on the right of (4.45), and a
   test pins that the two differ by a sign. Ω > 0 about +z in the Northern Hemisphere: the force deflects to the **right**.
6. **2-D stream function sign:** χ = −z ⇒ $\rho u=\partial\psi/\partial y$, $\rho v=-\partial\psi/\partial x$ (many GFD texts use
   the opposite sign — ⚠️ at C03). `flux_between_streamlines(psi, p1, p2)` counts flux crossing the segment p1 → p2 from
   its **left to its right** (normal n = (t_y, −t_x), t the unit tangent from p1 to p2), which makes it equal to
   ψ(p2) − ψ(p1). Axisymmetric: χ = −φ ⇒ $\rho u_R=-R^{-1}\partial\psi/\partial z$, $\rho u_z=R^{-1}\partial\psi/\partial R$.
7. **Control volumes.** Outward n on every CV; (u − b)·n signed; u and b in the same frame; the forces in (4.17) are the
   forces **on the fluid** (drag on the body is +F_D e_x, the force on the fluid −F_D e_x; ⚠️ at C04). The ch03
   `ControlVolume` shapes (`MovingBox`, `GrowingSphere`, …) with `surface_nodes(t, n) -> (x, n_hat, dA, b)` are reused;
   1-D budgets (per unit cross-section area) use `ch04.interval_mass_budget`.
8. **Symbols with several meanings** (notation table in the §4.1 section, and a ⚠️ where each second meaning appears):
   σ_ij viscous stress (C07) vs σ surface tension (C14) vs σ vortex core radius (C11, Ch. 3); ε dissipation (C10) vs
   ε_ijk (C08); **γ** in (4.29) (a Lamé-type coefficient) vs γ = C_p/C_v (C11, C15) vs Ex. 4.7's ζ/δ (C14) — so the
   **shear rate is written $\dot\gamma$** (`gamma_dot`) and the **extension rate s** in C07 and E3; Φ force potential
   (4.18) vs the unknown function in (4.99); Ψ vector potential (4.12) vs Ψ in (4.99)/(4.106); φ velocity potential vs
   azimuth; Ω frame rotation (§4.7) vs imposed frequency (§4.11); h enthalpy vs depth vs meniscus height vs gap;
   H angular momentum vs CV height; M torque vs rocket mass vs Mach number; R cylindrical radius vs R₁, R₂ vs gas
   constant; f traction vs Helmholtz free energy; κ thermal diffusivity (the book's stray κ in the (4.63) discussion
   means μ_v); λ second viscosity coefficient vs eigenvalue.
9. **Lapse rate and N²**: Kundu Γ ≡ dT/dz in code; wherever a lapse rate or N² is shown (C13, N137, E8, E9 stratified
   mode) the meteorological −dT/dz appears alongside (memory rule, ch01 E4 pattern).
10. **Colours (text, figures, explainers — one meaning each):** mass/momentum inflow `blue` · outflow `orange` · storage
    (d/dt∫) `accent` purple · force on the fluid `rose` · local ∂/∂t `blue` · advective (u·∇)u `teal` · pressure force
    `orange` · viscous force `rose` · gravity `muted` grey (NS bars) / buoyancy `blue` (E8) · Coriolis `amber` ·
    centrifugal `accent` purple · frame acceleration `muted` · kinetic (mechanical) energy `teal` · heat / dissipation
    `rose` · ½u² `teal`, p/ρ `orange`, gz `blue`, ∂φ/∂t `amber` in Bernoulli stacks · the three groups of (4.101):
    Ωl/U `blue`, gl/U² `muted`, μ/(ρUl) `rose` · ghosts and references `muted` grey dashed.
11. **Every book equation is shown in full next to its number** (CLAUDE.md rule 3) — in markdown, derivation steps
    ("substitute (4.27), $\tau_{ij}=-p\delta_{ij}+\sigma_{ij}$"), traps, recaps and in explainer tour, Explain, Derivation,
    notes, status and quiz text. A repeat mention in the same cell may use a compact inline form, never the bare number.
    Builders use `show_eqs(text, EQ)` (from `notebooks/build_ch03.py`) with an `EQ` dict covering all 124 labels; JS
    explainers use a local `showEqs`. `tools/eq_refs.py` must list 0 offenders.
12. **nbkit behaviour the builder must respect:** `nb.recap(...)` and `nb.section(...)` **end** the current CORE block
    and `nb.core(id)` can be called only once — so every RECAP sits **before** the `nb.core` call of the block that uses
    it (R01 before C01, R02 before C02, R03 R04 before C04, R05 R06 before C07, R07 R08 before C10, R09 before C12, R10
    before C13, R11 before C14, R12 before C15), and every `nb.derivation` sits inside its curation CORE block.
    `nb.derivation(ref=…)` gets an equation number only ("4.5"), never "Exercise 4.30". §4.9 rows taught earlier
    (N82–N84 in C04, N103–N105 in C05) get a one-line `nb.pointer` in the §4.9 section.
13. **Parsing.** Part E is the only part with table rows after the Part E heading; Part F has no line starting with
    `|`. Explainer headings are exactly `### E1 · control_volume_budgets` … `### E9 · dynamic_similarity_models`,
    `### B1 · kinematic_free_surface` (`tools/embed_check.py`). Primer terms in Part A are the exact Concept text of the
    Part E rows marked "primer" (coverage_check matches the first 18 characters). The "Explained by" column never names
    an earlier chapter's C-number (the checker would read it as this chapter's id): earlier chapters are cited by
    section ("Ch. 1 §1.5").
14. **Book slips taught corrected** (analysis §9 + one found here; each gets "the book prints X; the correct form is Y"
    where it is used): (4.15)'s trailing "= 0" (N21, D05 step 3) · (4.51)'s dA inside the volume integrals must be dV
    (N71, D19 step 5) · (4.74)'s gauge change must be φ → φ − ∫B dt′ (N93, D26 step 6) · the (4.63) discussion's "μ, κ,
    k > 0" means μ, μ_v, k ≥ 0 (N81) · §4.10's curve C must read ζ = x²/2R₁ **+** y²/2R₂ (N126) · Ex. 4.7 drops a minus
    sign in the separated equation and writes η for γ = h/δ (N129) · Ex. 4.2's statement "½ρU² + gz + p/ρ" is
    dimensionally inconsistent, the result (4.19) is ½U² + gz + p/ρ (D06 trap) · the Cauchy prose writes ∂τ_ij/∂x_j where
    (4.24) contracts the first index (D07 trap) · "S_mm = ∇·u … (see Section 3.6)" is §3.4, Eq. (3.14) (C07) · §4.8
    "multiplying (4.22) by u_j" — start from (4.24) (D21) · "[,]-brackets in (4.100)" are in (4.101); "(4.106), (4.107)"
    before (4.114) means (4.109), (4.112) (N150) · Ex. 4.8's total 9.14×10⁵ N is 9.15×10⁵ N unrounded (N156) · Fig. 4.9
    "budge" → bulge; "42 km" is 42.77 km with WGS-84 (N64) · **new (designer, p151):** the line after (4.57) writes
    $\sigma_{ij}(\partial u_j/\partial x_i)=\sigma_{ij}(S_{ji}+R_{ji})$; with the book's own $R=G-G^{\mathsf T}$ (3.13) it is
    $S_{ji}+\tfrac12R_{ji}$ — harmless because σ:R = 0, taught in D22 step 6.

Order of parts: C (contract) · A (notebook storyboard) · B (explainers) · D (runtime) · E (prerequisite ledger) · F
(derivations).

---

## Part C — functions the builders will call (the implementer's contract)

**Status column.** **§4 #n** = planned in `analysis/ch04.md` §4 row n with this signature (keep it). **§8** = added by
the curation's notes for the implementer (`analysis/ch04_curation.md` §8). **NEW** = added by this design (in neither) —
flagged as the phase asks. **reuse** = exists (ch01–ch03) and is only called. Every callable is reachable as
**`ch04.<name>`**: the chapter module `fluidpy/ch04_conservation_laws.py` re-exports every public name of the new core
modules below (so the notebook imports one module and parity rows write `ch04.…`). Module alias in the notebook:
`from fluidpy import ch04_conservation_laws as ch04`. Units SI: x [m], t [s], u [m/s], p [Pa], ρ [kg/m³], μ [Pa s],
ν [m²/s], G [1/s], T [K], k [W/(m K)]; docstrings cite § and Eq. (from the rendered page) and carry the validation label.
Return types: `NamedTuple` for term splits (unpackable, attribute access), `@dataclass(frozen=True)` for budgets and
`Scales` (the C01 primer teaches "results with named fields" on `MassBudget`).

### C.0 Reused (existing; called, not changed)
| # | Callable (signature) | Returns | Used by | Status |
|---|---|---|---|---|
| 0.1 | `style.setup_notebook() -> bool` (FAST) · `style.COLORS` · `style.savefig(fig, "ch04", name)` | FAST · palette · path | setup, every figure | reuse |
| 0.2 | `anim.animate(update, frames, fig, interval)` · `anim.show_animation(anim, player="video"\|"frames")` | HTML | C01, C04, C09, C12, C13, C14 | reuse |
| 0.3 | `interact.slider_figure(fn, name, values, *, unit, xlabel, ylabel, title, xrange, yrange, height)` · `interact.animate_figure` | plotly Figure | C02, C03, C04, C07, C08, C10, C11, C15 | reuse |
| 0.4 | `embed.show_viz("ch04", slug)` | display | 9 explainer cells | reuse |
| 0.5 | `tools.convergence.observed_order(h, err) -> float` | slope | C07 (cube, slope −2), C14 (cap force O(ζ)) | reuse |
| 0.6 | `core.transport`: `ControlVolume` shapes `MovingBox(lengths, rates, velocity, origin)`, `GrowingSphere(R0, Rdot, center)`; `volume_integral(F, cv, t, n)`, `volume_integral_rate_fd(F, cv, t, dt, n)`, `surface_flux_term(F, cv, t, n)`, `reynolds_transport(F, dFdt, cv, t, n)` | floats, RTTTerms | C01 (R01), C04, C10 (N68) | reuse |
| 0.7 | `core.kinematics`: `material_derivative(_terms)`, `material_derivative_sym`, `acceleration(u, x, t, h, ht)`, `velocity_gradient_at(u, x, t, h)`, `vorticity(u, x, t, h)`, `volumetric_strain_rate(G)`, `pathline(u, r0, t0, t_eval)`, `streamline(u, x0, t_frozen, s_max)`, `velocity_gradient_preset(name, Gamma)` | arrays | C01, C02, C06–C08, C11, C14 | reuse |
| 0.8 | `core.tensors`: `strain_rate_tensor(G)`, `rotation_tensor(G)`, `traction(tau, n)`, `normal_shear_stress(tau, n)`, `kronecker_delta()`, `levi_civita()`, `transform_tensor(T, C)` (any order), `random_rotation(rng)`, `is_isotropic(T, rng, n_rotations)`, `symmetric_double_contraction(tau, B)`, `double_dot(A, B, convention)`, `cross(u, v)` | arrays | C04 (R03), C07, C10 (R07) | reuse |
| 0.9 | `core.operators`: `divergence(u, h, bc)`, `tensor_divergence(T, h, bc, index=0)`, `laplacian`, `curl`, `is_solenoidal(u, h)` · `core.grids.grid2d` | grid arrays | C02 (N08, N11), C06 | reuse |
| 0.10 | `core.integral_theorems`: `divergence_theorem_box(Q_fn, bounds, n)`, `circulation(u_fn, loop)`, `planar_loop`, `rectangle_loop`, `curl_flux(u_fn, surface)`, `planar_rectangle` | floats | C02 (R02), C04 (N27), C03 (N16) | reuse |
| 0.11 | `core.statics.integrate_hydrostatic(z, rho_fn, p0, g)`, `hydrostatic_pressure_uniform(z, p0, rho, g)` · `core.stratification.brunt_vaisala_sq_from_lapse(T, dT_dz)`, `lapse_rate_convention(dT_dz, convention)` | arrays | C04 (bore sides), C13 (R10), C15 (N137) | reuse |
| 0.12 | `core.thermo`: `perfect_gas_sound_speed(T, gamma, R)`, `cp_from_gamma`, `R_AIR`, `CP_AIR`, `CV_AIR`, `GAMMA_AIR` · `core.dimensional.pi_groups(variables, solution, repeating)`, `SPHERE_DRAG` | floats, dicts | C11 (N96), C13, C15 (R12) | reuse |
| 0.13 | `ch01`: `newton_shear_stress(mu, du_dy)`, `laplace_pressure_jump(sigma, R1, R2)`, `surface_tension_water(T)`, `fluid_properties(name, T)`, `FLUIDS`, `knudsen_number(l, L)`, `seawater_density_linear(T, S)`, `fourier_heat_flux(k, grad_T)` · `ch03`: `cylinder_flow(x, y, U, a)`, `cylinder_streamfunction(x, y, U, a)`, `cylinder_velocity_field(U, a)`, `rigid_body_velocity(U, Omega, x)` · `core.vortices.rankine_vortex(r, Gamma, sigma)`, `gaussian_vortex` · `core.diffusion.couette_startup_profile` | floats, arrays | C07 (N51), C14 (N127), C15 (N152), C03, C05, C11, C08 | reuse |

### C.1 `fluidpy/core/conservation.py` (NEW module, `core.CB`) — integral budgets on a control volume
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 1.1 | `MassBudget` dataclass: `storage` (d/dt∫_{V*}ρ dV by central difference in t), `outflux` (∮ρ(u − b)·n dA), `local` (∫∂ρ/∂t dV, NaN if not given), `residual` (= storage + outflux) | (4.5) | kg/s (1-D: kg/(m² s)) | C01 | §4 #3 (field names fixed here) |
| 1.2 | `mass_budget(rho, u, cv, t, drho_dt=None, material=False, n=24, dt=1e-4) -> MassBudget` (`material=True` uses b = u on the surface instead of the CV's own b: the (4.1) check) | (4.5), (4.1), (4.2) | kg/s | C01 (3-D check), E1 parity via 1.6 | §4 #3 |
| 1.3 | `MomentumBudget` dataclass: `storage` (3,), `outflux` (3,) (∮ρu(u − b)·n dA), `body` (3,) (∫ρg dV), `surface` (3,) (∮f dA, f = traction), `residual` (3,) (= storage + outflux − body − surface) | (4.17) | N | C04 | §4 #12 |
| 1.4 | `momentum_budget(rho, u, cv, t, g=(0.0, 0.0, -9.81), traction=None, n=24, dt=1e-4) -> MomentumBudget` (`traction(x, n_hat, t) -> (3, N)`; None = no surface force) | (4.17), (4.13) | N | C04 (3-D check) | §4 #12 |
| 1.5 | `EnergyBudget` dataclass (`storage`, `outflux`, `body_work`, `surface_work`, `heat_out`, `residual`) · `energy_budget(rho, u, e, cv, t, g=(0,0,-9.81), traction=None, q=None, n=24, dt=1e-4) -> EnergyBudget` | (4.48) | W | C10 (N68 one-line check) | §4 #35 |
| 1.6 | `interval_mass_budget(x0, x1, t, dx0dt=0.0, dx1dt=0.0, flow="expanding", **p) -> dict(storage, flux_left, flux_right, net_outflux, local, residual)` — 1-D form of (4.5) per unit cross-section area for the interval [x0, x1] whose ends move at dx0dt, dx1dt; `flow` = a preset name (`"expanding"` with `a`, `rho0`) or a pair `(rho_fn, u_fn)` of callables of (x, t); `flux_left` = ρ(u − ẋ₀) at x0 (inflow when > 0), `flux_right` = ρ(u − ẋ₁) at x1, `net_outflux` = flux_right − flux_left, `storage` = d/dt∫ρ dx (exact for the preset, central difference otherwise), `local` = ∫∂ρ/∂t dx | (4.5) in 1-D, (4.1) with ẋ = u | kg/(m² s) | C01 worked number, E1 mass mode, parity | **NEW** (the C01 tiny example and a no-lambda parity row need a 1-D budget) |
| 1.7 | `angular_momentum_flux(rho, u, cv, t, origin=(0,0,0), n=24) -> ndarray (3,)` (∮(r × ρu)(u·n) dA) | (4.65) | N m | C04 (N83 check of the sprinkler) | §4 #40 |

### C.2 `fluidpy/core/streamfunction.py` (NEW module, `core.SF`)
`psi` arguments accept **either a callable ψ(x, y) or a preset name** with keywords (see 2.6).
| # | Callable (signature) | Implements (Eq.) | Returns | Used by | Status |
|---|---|---|---|---|---|
| 2.1 | `velocity_from_streamfunction_2d(psi, x, y, rho=1.0, h=1e-5, **p) -> (u, v)` (central differences; ρu = ∂ψ/∂y, ρv = −∂ψ/∂x) · `velocity_from_streamfunction_2d_sym(psi_expr, x, y, rho=1) -> (u_expr, v_expr)` | C03 2-D ψ, D04 | m/s | C03, E2 parity | §4 #10 |
| 2.2 | `flux_between_streamlines(psi, p1, p2, n=400, rho=1.0, **p) -> float` (∫ρu·n ds along the straight segment p1 → p2 by Gauss–Legendre, n = (t_y, −t_x): flux from left to right of the direction of travel; equals ψ(p2) − ψ(p1)) · `flux_along_path(psi, pts, rho=1.0, **p) -> float` (same along a polyline `pts` (2, m), for E2's bent gate) | D04 step 7, Ex. 4.8 | m²/s (per unit depth) | C03 from-scratch, E2 | §4 #10 (+ `flux_along_path` **NEW**) |
| 2.3 | `velocity_from_streamfunction_axisym(psi, R, z, rho=1.0, h=1e-5, **p) -> (u_R, u_z)` (R > 0 enforced) · `..._axisym_sym(psi_expr, R, z, rho=1) -> (uR_expr, uz_expr)` | N17 | m/s | C03 (N17), E2 axisym preset | §4 #11 |
| 2.4 | `mass_flux_from_vector_potential(Psi_exprs, coords) -> list[sympy]` (curl) · `mass_flux_from_stream_functions(chi, psi, coords) -> list[sympy]` (∇χ × ∇ψ) · `stream_function_pair(name="parabolic") -> (chi_expr, psi_expr, coords)` — `"parabolic"`: χ = y, ψ = z − x², so ρu = ∇χ × ∇ψ = (1, 0, 2x) (streamlines z = x² + c in planes y = const; ∇·(ρu) = 0) | (4.12), N14, N15 | sympy | C03 (N14, N15 plotly 3-D) | §4 #9 (+ `stream_function_pair` **NEW**: a concrete 3-D pair for Fig. 4.1) |
| 2.5 | `stream_surface_check(chi, psi, x) -> dict(u_dot_grad_chi, u_dot_grad_psi, div)` (numeric at points; chi, psi callables or `"parabolic"`) · `stream_tube_mass_flux(chi, psi, a, b, c, d, n=64) -> (numeric, closed_form)` (flux through the patch bounded by χ = a, b, ψ = c, d vs (b − a)(d − c); `"parabolic"` pair by default) | N15, N16 | –; kg/s | C03 (N16) | §4 #9 |
| 2.6 | `streamfunction_preset(name, x, y, **p) -> psi` and `velocity_preset(name, x, y, **p) -> (u, v)` (closed-form velocities of the same presets) — `"uniform"` (U, alpha: ψ = U(y cos α − x sin α)), `"stagnation"` (k: ψ = kxy), `"source_stream"` (U, m: ψ = Uy + (m/2π)θ), `"cylinder"` (U, a: = `ch03.cylinder_streamfunction`), `"vortex"` (Gamma: ψ = −(Γ/2π) ln r), `"shear"` (gamma_dot: ψ = ½γ̇y²), `"axisym_uniform"` (U: ψ = ½UR², x ↦ z, y ↦ R) · `STREAMFUNCTION_PRESETS` (names, defaults, singular points, stagnation points) | E2 flows | m²/s; m/s | C03 figures, E2 | §8 (+ `velocity_preset` **NEW**: exact velocities for parity rows) |

### C.3 `fluidpy/core/constitutive.py` (NEW module, `core.NW`)
| # | Callable (signature) | Implements (Eq.) | Returns | Used by | Status |
|---|---|---|---|---|---|
| 3.1 | `static_stress(p) -> (3,3)` (−pδ) · `total_stress(p, sigma) -> (3,3)` (−pδ + σ) | (4.26), (4.27) | Pa | C07 (N40, N41) | §4 #23 |
| 3.2 | `linear_stress(K, S) -> (3,3)` (`np.einsum("ijmn,mn->ij")`) · `isotropic_fourth_order(lam, mu, gam) -> (3,3,3,3)` | (4.28), (4.29) | Pa | C07 (D09 check), E3 | §4 #24 |
| 3.3 | `newtonian_stress(G, p=0.0, mu=1.0e-3, lam=None, mu_v=None, incompressible=False) -> (3,3)` — exactly one of `lam`, `mu_v` may be given; both None ⇒ Stokes (μ_v = 0); `incompressible=True` raises `ValueError` if \|tr G\| > 1e-9·‖G‖ | (4.31), (4.35), (4.37) | Pa | C07, E3 parity | §4 #25 |
| 3.4 | `viscous_stress(G, mu, mu_v=0.0) -> (3,3)` (σ of (4.59)) · `mean_pressure(tau) -> float` · `thermodynamic_pressure_from_stress(tau, div_u, mu, lam) -> float` · `pressure_difference(div_u, mu, lam) -> float` · `bulk_viscosity(lam, mu) -> float` · `lam_from_bulk(mu_v, mu) -> float` | (4.59), (4.32)–(4.34), N49, (4.36) | Pa, Pa s | C07, C10 | §4 #25 |
| 3.5 | `stress_on_plane(G, p=0.0, mu=1.0e-3, mu_v=0.0, theta=0.0) -> (sigma_n, tau_s)` — plane with unit normal n = (cos θ, sin θ, 0), tangent t = (−sin θ, cos θ, 0); f_j = n_iτ_ij; σ_n = f·n (tension positive), τ_s = f·t (**signed**) | traction of (4.31) | Pa | C07 figure, E3 parity | §8 (signature fixed: angle θ instead of a vector, signed shear) |
| 3.6 | `dissipation_rate(G, rho, mu, mu_v=0.0, form="sum_of_squares") -> float` (`"contraction"`: σ_ijS_ij/ρ; `"sum_of_squares"`: 2ν(S − ⅓S_mmδ)² + (μ_v/ρ)S_mm²; `"both"` returns the pair) | (4.58) | W/kg | C10, E7 parity | §4 #38 |
| 3.7 | `deviatoric_part(A) -> (3,3)` (A − ⅓ tr A δ) | D10, D23 | as A | C07 primer, C10 | **NEW** (the deviatoric primer and D23 check call it) |

### C.4 `fluidpy/core/navier_stokes.py` (NEW module, `core.NS`) — residuals and term splits (numeric stencils + sympy)
Stencil functions take `h` (space, default 1e-4 m) and `ht` (time, default `core.kinematics.DEFAULT_HT`); never reuse h
as a time step. Term splits are `NamedTuple`s whose fields are (3,) vectors unless stated.
| # | Callable (signature) | Implements (Eq.) | Returns | Used by | Status |
|---|---|---|---|---|---|
| 4.1 | `continuity_residual(rho, u, x, t, h=1e-4, ht=None) -> float` (∂ρ/∂t + ∇·(ρu)) · `continuity_residual_sym(rho_expr, u_exprs, coords, t) -> sympy.Expr` | (4.7) | kg/(m³ s) | C02 | §4 #4 |
| 4.2 | `ContinuityTerms(local, flux_divergence, residual)` · `continuity_terms(rho, u, x, t, h=1e-4, ht=None) -> ContinuityTerms` (∂ρ/∂t and ∇·(ρu) separately) | (4.7), N08 | kg/(m³ s) | C02 figure | **NEW** (the C02 slider shows the two terms of (4.7); analysis' `continuity_terms` becomes 4.3) |
| 4.3 | `continuity_material_terms(rho, u, x, t, h=1e-4, ht=None) -> (rel_rate, div_u)` ((1/ρ)Dρ/Dt and ∇·u; sum 0) · `density_material_rate(rho, u, x, t, h=1e-4, ht=None) -> float` (Dρ/Dt) · `divergence_free_check(u, x, t=0.0, h=1e-4, tol=1e-8) -> (bool, value)` | (4.8), (4.9), (4.10) | 1/s | C02 (N09–N11) | §4 #5–#7 (**renamed** from `continuity_terms`) |
| 4.4 | `momentum_conservative_residual(rho, u, tau, g, x, t, h=1e-4, ht=None) -> (3,)` · `conservative_to_advective_sym(rho_expr, u_exprs, coords, t) -> list[sympy]` (difference of the flux form and ρDu_j/Dt, which simplifies to u_j × continuity) | (4.22), (4.23) | N/m³ | C06 (N36, N37) | §4 #19, #20 |
| 4.5 | `CauchyTerms(inertia, body, stress_divergence, residual)` · `cauchy_terms(rho, u, tau, g, x, t, h=1e-4, ht=None) -> CauchyTerms` (per unit volume; stress divergence over the **first** index) · `cauchy_residual(...) -> (3,)` | (4.24) | N/m³ | C06, E4 | §4 #21 |
| 4.6 | `navier_stokes_residual(rho, u, p, x, t, mu, mu_v=0.0, g=(0,0,-9.81), h=1e-4, ht=None) -> (3,)` (variable μ allowed as a callable μ(x, t)) · `navier_stokes_sym(rho_expr, u_exprs, p_expr, coords, t, mu, mu_v=0, g=(0,0,-g)) -> list[sympy]` | (4.38) | N/m³ | C08 (D11 check) | §4 #26 |
| 4.7 | `NSTerms(local, advective, pressure, gravity, viscous, residual)` · `ns_incompressible_terms(u, p, x, t, rho, mu, g=(0,0,-9.81), h=1e-4, ht=None, per="mass") -> NSTerms` (per="mass": ∂u/∂t, (u·∇)u, −∇p/ρ, g, ν∇²u [m/s²]; per="volume": ×ρ [N/m³]); residual = local + advective − pressure − gravity − viscous | (4.39b), (4.41) (mu=0), (4.85) | m/s² or N/m³ | C08, C13 (N107), E4 | §4 #27 (`per` switch **NEW**) |
| 4.8 | `viscous_force_forms(u, x, t=0.0, mu=1.0, h=1e-4) -> (lap, div2S, minus_curl_omega)` (μ∇²u, 2μ∂S_ij/∂x_i, −μ∇×ω, each (3,)) | (4.40) | N/m³ | C08 (N54), E4 | §4 #28 |
| 4.9 | `exact_solution(name, x, t=0.0, **p) -> (u, p)` and `exact_field(name, **p) -> (u_fn, p_fn)` — `"couette"` (U, h; walls y = 0, h), `"poiseuille"` (G = −dp/dx, h, mu, rho=1000.0, g=0.0; walls y = 0, h; u = G y(h − y)/(2μ); p = −Gx − ρgy, so g ≠ 0 means gravity along −y for the N107 demo), `"stokes_first"` (U, nu: u = U erfc(y/2√(νt))), `"taylor_green"` (U0, k, nu, rho), `"lamb_oseen"` (Gamma, nu, t0: Gaussian vortex with σ² = 4ν(t + t0)), `"cylinder"` (U, a, rho, p_inf: potential flow + Bernoulli pressure, an Euler solution), `"solid_body"` (Omega, rho, g=9.81: p = ρΩ²R²/2 − ρgz) · `ns_terms_preset(name, x, y, t=0.0, component=0, per="mass", form="laplacian", **p) -> dict(local, advective, pressure, gravity, viscous, residual)` (scalar floats of one component at (x, y) — the parity-friendly wrapper of 4.7 on 4.9; `form` = `"laplacian"` (μ∇²u), `"div2S"` (2μ∂S_ij/∂x_i) or `"curl"` (−μ∇×ω) chooses how the viscous term is computed) | E4 solutions, C06, C08 | m/s, Pa | C06, C08, E4 | §8 (+ `exact_field`, `"solid_body"`, `ns_terms_preset` **NEW**) |
| 4.10 | `lamb_vector(u, x, t=0.0, h=1e-4) -> (3,)` (u × ω; a 2-D field u(x) with x of shape (2,) is treated as the plane z = 0, ω = (0, 0, ω₃)) · `lamb_identity_terms(u, x, t=0.0, h=1e-4) -> dict(advective, minus_u_cross_omega, grad_ke, residual)` · `lamb_identity_sym(u_exprs, coords) -> list[sympy]` | (4.68) | m/s² | C11 (D24 check) | §4 #43 |
| 4.11 | `stress_work_split(p, sigma, u, x, t=0.0, h=1e-4) -> dict(total, deformation_pressure, deformation_viscous, force_pressure, force_viscous, residual)` | (4.54) | W/m³ | C10 (N74) | §4 #36 |
| 4.12 | `kinetic_energy_budget(rho, u, p, sigma, g, x, t, h=1e-4, ht=None) -> dict(lhs, gravity_work, pressure_work, viscous_force_work, residual)` | (4.56) | W/m³ | C10 (D21 check) | §4 #36 |
| 4.13 | `internal_energy_terms(rho, u, p, sigma, q, x, t, h=1e-4, ht=None) -> dict(De_Dt, pressure_work, dissipation, conduction, residual)` (per unit mass, (4.57)) · `internal_energy_residual(rho, u, e, p, T, mu, mu_v, k, x, t, h=1e-4, ht=None) -> float` ((4.60)) · `total_energy_residual_sym(...)` · `energy_forms_sym() -> sympy.Expr` (0 when (4.60) ⇔ (4.112)) | (4.57), (4.60), (4.53), (4.112) | W/kg, W/m³ | C10, C15 (N148) | §4 #36, #37 |
| 4.14 | `entropy_terms(rho, T, q, eps, x, h=1e-4) -> dict(flux_divergence, conduction_production, dissipation_production)` · `entropy_production(gradT, T, k, rho, eps) -> float` (k\|∇T\|²/(ρT²) + ε/T) | (4.62), (4.63) | W/(kg K) | C10 (N80, N81), E7 | §4 #39 |
| 4.15 | `perturbation_fields(p, rho, z, rho_s, g=9.81, p_s0=None) -> (p_pert, rho_pert)` · `buoyancy(rho_pert, rho0, g=9.81) -> float` (b = −gρ′/ρ₀, upward positive) · `boussinesq_momentum_terms(u, p_pert, rho_pert, x, t, rho0, nu, g=9.81, h=1e-4, ht=None) -> dict(local, advective, pressure, buoyancy, viscous, residual, dropped)` (`dropped` = (ρ′/ρ₀)Du/Dt, the term (4.86) neglects) | (4.84), (4.86) | Pa, m/s² | C13, E8 | §4 #50, #52 |
| 4.16 | `heat_equation_terms(T, u, x, t, rho, cp, k, eps=0.0, h=1e-4, ht=None) -> dict(DT_Dt, dissipation, conduction, residual)` · `temperature_equation_residual(T, u, x, t, kappa, h=1e-4, ht=None) -> float` | (4.88), (4.89) | K/s | C13 | §4 #53, #54 |
| 4.17 | `divergence_first_index_demo(tau_fn, x, h=1e-4) -> (first, second)` (∂τ_ij/∂x_i and ∂τ_ij/∂x_j at a point — they differ for a non-symmetric τ) | ⚠️ first index | N/m³ | C06 primer, from-scratch | **NEW** (the first-index primer's numeric demo) |

### C.5 `fluidpy/core/rotating.py` (NEW module, `core.RF`)
| # | Callable (signature) | Implements (Eq.) | Returns | Used by | Status |
|---|---|---|---|---|---|
| 5.1 | `OMEGA_EARTH = 7.292115e-5` rad/s · `EARTH_RADIUS = 6.378137e6` m (WGS-84 a) | N61, N64 | – | C09 | §4 #33 |
| 5.2 | `rotating_basis(Omega_vec, t) -> (3,3)` (rows e′_i(t) for constant Ω, Rodrigues, e′_i(0) = e_i) · `basis_rate(Omega_vec, t, h=1e-6) -> (3,3)` (central difference of the rows; equals Ω × e′_i) · `inertial_velocity(U, u_rot, Omega, x_rot) -> (3,)` (U + u′ + Ω × x′) | (4.42), D14 | –, 1/s, m/s | C09 (D14 check) | §4 #29 (+ `basis_rate` **NEW**) |
| 5.3 | `frame_acceleration_terms(a_rot, u_rot, x_rot, Omega, dOmega_dt=(0,0,0), dU_dt=(0,0,0)) -> dict(frame, relative, coriolis, angular, centripetal, total)` — the **accelerations** of (4.43): dU/dt, a′, +2Ω × u′, (dΩ/dt) × x′, Ω × (Ω × x′), and their sum a | (4.43), (4.44) | m/s² | C09, E5 parity | §4 #30 |
| 5.4 | `apparent_body_forces(u_rot, x_rot, Omega, dOmega_dt=(0,0,0), dU_dt=(0,0,0), g=(0,0,-9.81)) -> dict(gravity, frame, coriolis, angular, centrifugal, total)` — the **forces per unit mass** in the bracket of (4.45): g, −dU/dt, −2Ω × u′, −(dΩ/dt) × x′, −Ω × (Ω × x′) | (4.45) | m/s² | C09, E5 | §4 #31 |
| 5.5 | `coriolis_acceleration(Omega, u) -> (3,)` (**−2Ω × u**, the Coriolis force per unit mass of (4.45)) · `coriolis_parameter(lat) -> float` (2Ω sin φ, forward pointer only) | N61, N62 | m/s², 1/s | C09, E5 | §4 #32, §8 |
| 5.6 | `centrifugal_acceleration(Omega, x) -> (3,)` (−Ω × (Ω × x) = Ω²R e_R) · `centrifugal_potential(R, Omega) -> float` (−½Ω²R²) · `effective_gravity(lat, g_n=9.80, Omega=OMEGA_EARTH, a=EARTH_RADIUS) -> (g_e, deflection)` (magnitude m/s² and angle between g_e and the radial direction [rad], spherical Earth) | N64, D18 | m/s², J/kg, rad | C09, E5 | §4 #33 |
| 5.7 | `projectile_paths(u0, Omega, t) -> (inertial, rotating)` (each (2, n): the straight inertial path (u₀t, 0) from the pole and the same path in the frame turning at Ω: (u₀t cos Ωt, −u₀t sin Ωt)) | Fig. 4.8, D17 | m | C09 animation, E5 parity | §8 |

### C.6 `fluidpy/core/curvilinear.py` (NEW module, `core.CU`) — Appendix-B operators (sympy)
| # | Callable (signature) | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 6.1 | `gradient(f, system)`, `divergence(v, system)`, `curl(v, system)`, `laplacian(f, system)`, `vector_laplacian(v, system)`, `advective_acceleration(v, system)`, `strain_rate(v, system)` — `system="cylindrical"` (R, φ, z) or `"spherical"` (r, θ from +z, φ); inputs are sympy expressions in `coordinates(system)` | Appendix B | sympy | C09 (N65) | §4 #34 |

### C.7 `fluidpy/core/bernoulli.py` (NEW module, `core.BE`)
| # | Callable (signature) | Implements (Eq.) | Returns | Used by | Status |
|---|---|---|---|---|---|
| 7.1 | `bernoulli_head(U, z, p, rho, g=9.81) -> float` (½U² + gz + p/ρ) · `bernoulli_solve(state1, state2, unknown, rho=1000.0, g=9.81) -> float` (dicts with keys `U`, `z`, `p`; the one key named by `unknown` is solved from (4.19)) | (4.19) | m²/s² | C05, E6 | §4 #16 |
| 7.2 | `pitot_speed(p_stag, p_static, rho) -> float` · `pitot_speed_from_heads(h1, h2, g=9.81, rho=1000.0, rho_atm=0.0) -> float` · `stagnation_pressure(p, U, rho) -> float` · `dynamic_pressure(U, rho) -> float` · `torricelli_speed(h, g=9.81) -> float` | N103–N105 | m/s, Pa | C05, E6 | §4 #48, #49 |
| 7.3 | `pressure_function(p, p_o, kind="constant", rho=None, T=None, R=287.058, gamma=1.4, rho_o=None) -> float` (∫_{p_o}^{p} dp′/ρ(p′): `"constant"` (p − p_o)/ρ, `"isothermal"` RT ln(p/p_o), `"isentropic"` γ/(γ − 1)(p/ρ − p_o/ρ_o) with ρ = ρ_o(p/p_o)^{1/γ}, or a callable ρ(p) by `quad`) | (4.67) | J/kg | C11 (N87) | §4 #42 |
| 7.4 | `bernoulli_function(speed, p, z, rho=1000.0, g=9.81, kind="constant", p_o=0.0, **state) -> float` (B = ½U² + ∫dp/ρ + gz, the bracket of (4.69)) · `bernoulli_along_line(u, p, points, rho=1000.0, g=9.81) -> ndarray` (B at (d, N) points of callables u, p; z = last coordinate) · `lamb_surface_check(u, B, x, h=1e-4) -> float` (\|∇B − u × ω\|) | (4.69)–(4.72) | m²/s² | C11, E6 | §4 #44 |
| 7.5 | `unsteady_bernoulli_pressure(dphi_dt, speed, z, rho=1000.0, g=9.81, C=0.0) -> float` (p = ρ(C − ∂φ/∂t − ½U² − gz)) · `unsteady_bernoulli_B(phi, p, x, t, rho=1000.0, g=9.81, h=1e-4, ht=None) -> float` (the bracket of (4.74) at a point) · `gauge_absorbed_bracket(B, sign=-1.0) -> float` (value of the bracket after φ → φ + sign·∫B dt′: 0 for sign = −1, 2B for the printed +) · `unsteady_streamline_bernoulli(dudt_along, s, state1, state2, rho=1000.0, g=9.81) -> float` (residual of (4.82): ∫₁²∂u/∂t·ds + (½U² + gz + p/ρ)₂ − (…)₁; `dudt_along` samples on `s`) · `viscous_irrotational_residual(u, x, t=0.0, mu=1.0e-3, h=1e-4) -> (3,)` (−μ∇×ω) | (4.74)–(4.75), (4.80), (4.82) | Pa, m²/s² | C12, E6 | §4 #45 (+ `gauge_absorbed_bracket` **NEW**: the (4.74) sign test in one call) |
| 7.6 | `stagnation_enthalpy(h, U, z=0.0, g=9.81) -> float` · `stagnation_temperature(T, U, cp=1004.5) -> float` | (4.78) | J/kg, K | C11 (N96), E6 | §4 #46 |
| 7.7 | `BERNOULLI_FORMS` (dict label → dict(tex, hypotheses, constant_along)) for `"4.19"`, `"4.71"`, `"4.72"`, `"4.75"`, `"4.78"`, `"4.82"` · `which_bernoulli(steady, viscous, irrotational, barotropic, isentropic=False, constant_density=False) -> list[str]` (labels that apply, in the order above) · `which_bernoulli_text(...) -> str` (the same labels comma-joined, for exact-text parity rows) | N102 | – | C11 table, E6 status and parity | §4 #47 |
| 7.8 | `bernoulli_scenario(name, **p) -> dict(s_along, B_along, s_across, B_across, valid, status, terms)` — `"pitot"` (U, rho), `"orifice"` (h), `"rankine"` (Gamma, sigma, rho), `"cylinder"` (U, a), `"u_tube"` (L, h0, t), `"hot_nozzle"` (T0, U); `terms` = dict of ½u², p/ρ, gz, ∂φ/∂t at the probe | E6 scenarios | m²/s² | E6 | §8 |

### C.8 `fluidpy/core/interfaces.py` (NEW module, `core.IF`)
`eta` arguments are callables η(x, t) or the preset names `"moving_wall"` (V: η = x − Vt) and `"linear_wave"` (a, k,
H: η = z − a cos(kx − ωt)).
| # | Callable (signature) | Implements (Eq.) | Returns | Used by | Status |
|---|---|---|---|---|---|
| 8.1 | `surface_normal_speed(eta, x, t, h=1e-4, ht=None, **p) -> float` (u_s·n = −(∂η/∂t)/\|∇η\|) · `kinematic_bc_residual(eta, u, x, t, h=1e-4, ht=None, **p) -> float` (Dη/Dt) · `relative_normal_velocity(eta, u, x, t, h=1e-4, ht=None, **p) -> float` ((1/\|∇η\|)Dη/Dt) · `interface_mass_flux(eta, u, rho, x, t, h=1e-4, ht=None, **p) -> float` | (4.90)–(4.93) | m/s, kg/(m² s) | C14, B1 | §4 #57 |
| 8.2 | `pillbox_limit(flux_top, flux_bottom, side_rate, volume_rate, l) -> dict(jump, side, volume, residual)` (the pillbox balance for a height l; side and volume terms ∝ l) | N116 | as fluxes | C14 (N116) | §4 #56 |
| 8.3 | `cap_pressure_force(dp, R1, R2, zeta, exact=True) -> float` (z-force −Δp × ellipse area; `exact=True` by `dblquad` over the cap) · `cap_surface_tension_force(sigma, R1, R2, zeta, exact=True) -> float` (σ∮t × n ds round the exact ellipse, x = a cos s, y = b sin s) · `laplace_jump_from_balance(sigma, R1, R2, zeta=1e-6) -> float` (Δp that makes the two forces cancel) | (4.97), (4.98), (1.5) | N, Pa | C14 (N125–N127) | §4 #60 |
| 8.4 | `capillary_length(sigma, rho, g=9.81, rho_other=0.0) -> float` | N128 | m | C14, C15 (N154) | §4 #61 |

### C.9 `fluidpy/core/similarity.py` (NEW module, `core.SIM`)
| # | Callable (signature) | Implements (Eq.) | Returns | Used by | Status |
|---|---|---|---|---|---|
| 9.1 | `Scales` frozen dataclass `(l, U, rho, mu, g=9.81, Omega=None, p_inf=0.0, time_scale="omega", pressure_scale="dynamic", c=None, T_o=None, T_w=None, cp=None, k=None)` with properties `St`, `Re`, `Fr`, `M`, `Ec`, `Pr` and methods `nondimensionalise(x=None, t=None, u=None, p=None) -> dict`, `redimensionalise(...) -> dict`, classmethod `from_oscillation(l, Omega, rho, mu, g=9.81)` (U = lΩ) | (4.100), (4.109), (4.113), N141 | – | C15, E9 | §4 #64 |
| 9.2 | `nondimensional_ns_coefficients(pressure_scale="dynamic") -> dict(unsteady, advective, pressure, gravity, viscous)` (sympy coefficients in symbols Ω, l, U, g, μ, ρ after dividing by ρU²/l; `"viscous"`/`"hydrostatic"` pressure scalings give the Exercise-4.59 variants) · `nondimensional_ns_sym(pressure_scale="dynamic") -> sympy.Eq` (the whole of (4.101)) · `nondimensional_energy_coefficients() -> dict` ((4.114): Ec, Ec/Re, 1/(Pr Re)) | (4.101), (4.114) | sympy | C15 (D30 check) | §4 #65 (+ `nondimensional_ns_sym` **NEW**) |
| 9.3 | `strouhal_number(Omega, l, U)` · `reynolds_number(U, l, nu)` · `froude_number(U, l, g=9.81)` · `reduced_gravity(rho1, rho2, g=9.81)` (g(ρ₂ − ρ₁)/ρ₁) · `internal_froude_number(U, g_prime, l)` · `richardson_number(g_prime, l, U)` (g′l/U² = 1/Fr′²) · `gradient_richardson_number(N2, dU_dz)` · `mach_number(U, c=None, T=None, gamma=1.4, R=287.058)` · `compressibility_parameter(U, c)` (M²) · `eckert_number(U, cp, dT)` · `prandtl_number(nu, kappa)` · `eucken_prandtl(gamma)` (4γ/(9γ − 5)) · `weber_number(rho, U, l, sigma)` · `bond_number(rho, g, l, sigma)` · `capillary_number(mu, U, sigma)` · `rossby_number(U, Omega, l)` (U/(2Ωl): the advective term U²/l over the Coriolis term 2ΩU of (4.45); Ch. 13 writes U/(fL) with f = 2Ω sin φ — equal at the pole; forward pointer) — all scalar-callable floats | (4.102)–(4.105), (4.110), (4.111), (4.115)–(4.119) | – | C15, E9, C02 (N12) | §4 #8, #66, §8 |
| 9.4 | `pressure_coefficient(p, p_inf, rho, U)` · `drag_coefficient(F, rho, U, A)` · `lift_coefficient(F, rho, U, A)` · `reference_area(shape, **dims)` (`"sphere"`: πd²/4; `"plate"`: planform b·c) | (4.106)–(4.108) | – | C15 | §4 #67, #68 |
| 9.5 | `sphere_drag_coefficient(Re, model="morrison")` (`"stokes"`: 24/Re) · `froude_scaled_speed(U_p, l_p, l_m, g_p=9.81, g_m=9.81)` · `model_prototype(l_p, U_p, scale, fluid_p="water", fluid_m="water", match="Fr", T=293.15) -> dict(U_m, Re_p, Re_m, Fr_p, Fr_m, Re_ratio, Fr_ratio, matched)` (scale = l_m/l_p; ν from `ch01.fluid_properties`) | Fig. 4.21, Ex. 4.8 | –, m/s | C15, E9 | §4 #63, #70, §8 |

### C.10 `fluidpy/core/thermo.py` (addition)
| # | Callable (signature) | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 10.1 | `helmholtz_free_energy(e, T, s) -> float` (f = e − Ts) | (4.94) | J/kg | C14 (N121) | §4 #58 |

### C.11 `fluidpy/ch04_conservation_laws.py` (chapter module; re-exports C.1–C.10)
| # | Callable (signature) | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 11.1 | `closure_count(stage) -> dict(equations, unknowns, equation_names, unknown_names)` — `"cauchy"` (6, 13), `"navier_stokes"` (4, 5), `"barotropic"` (5, 5), `"full"` (7, 7); `"empty"` returns zeros (the §4.1 table skeleton) | N02 | dict | §4.1, C06, C08, C10 | §4 #1 (+ `"empty"` **NEW**) |
| 11.2 | `expanding_flow(x, t, a=1.0, rho0=1.0, dim=1) -> (u, rho)` (u = ax/(1 + at) per component, ρ = ρ₀/(1 + at)^dim) · `material_mass(X0, X1, t, a=1.0, rho0=1.0, dim=1) -> float` (mass per unit area of the material interval that was [X0, X1] at t = 0, by `quad`) · `material_interval(X0, X1, t, a=1.0) -> (x0, x1)` (= X(1 + at)) | N03, C01 | m/s, kg/m³, kg/m² | C01, C02, E1 | §4 #2 (+ `material_interval` **NEW**) |
| 11.3 | `stratified_shear_flow(z, U0=1.0, shear=0.1, rho0=1000.0, drho_dz=-0.5) -> (u, rho)` (u = U₀ + shear·z along x, ρ = ρ₀ + (dρ/dz)z) | N10 | m/s, kg/m³ | C02 | §4 #6 |
| 11.4 | `is_incompressible_regime(U, T=288.15, threshold=0.3) -> bool` | N12 | – | C02 | §4 #8 |
| 11.5 | `ball_integral(f, x0, radius) -> float` (∫ over the ball of radius r about x0 by `tplquad` in spherical coordinates) | N07 demo | [f]·m³ | C02 (localisation demo) | **NEW** |
| 11.6 | `body_force_from_potential(Phi, x, h=1e-4) -> (3,)` (−∇Φ by stencils; Phi callable or `"gravity"`) · `gravity_potential(z, g=9.81) -> float` | (4.18) | m/s², J/kg | C04 (N27) | §4 #13 |
| 11.7 | `gaussian_wake(y, U_inf=10.0, deficit=2.0, width=0.1) -> U` (U = U∞ − Δ exp(−y²/b²)) · `wake_side_outflow(U_of_y, U_inf, H) -> float` (∫(U∞ − U)dy over [−H/2, H/2]) · `wake_drag_per_span(U_of_y, U_inf, rho, H, return_error=False) -> float` (ρ∫U(U∞ − U)dy; `U_of_y` callable or `"gaussian"` with the keywords of `gaussian_wake`) | Ex. 4.1 (N29) | m/s, m²/s, N/m | C04, E1 | §4 #14 |
| 11.8 | `bore_speed(h_in, h_out, g=9.81) -> float` · `bore_pressure_force(h_in, h_out, rho=1000.0, g=9.81, p_o=1.0e5, width=1.0) -> dict(left, right, top, net)` (p_o drops out of `net`) | Ex. 4.3 (N31) | m/s, N | C04, E1 | §4 #17 |
| 11.9 | `rocket_trajectory(M0, mdot, Ve, t_burn, g=9.81, Fs=0.0, t_eval=None) -> dict(t, z, b, M)` (stops at burn-out) · `rocket_delta_v(M0, M1, Ve) -> float` | Ex. 4.4 (N32) | s, m, m/s, kg | C04, E1 | §4 #18 |
| 11.10 | `jet_plate_force(rho, V, A, theta=np.pi/2) -> float` (ρV²A sin θ, normal force on a plate at angle θ to the jet) · `cv_scenario(name, **p) -> dict(faces, mass_out, momentum_out_x, body_x, surface_x, storage_x, residual_mass, residual_momentum, result, result_label)` — `"wake"` (U_inf, deficit, width, rho, H), `"bore"` (h_in, h_out, b = CV speed, rho), `"jet"` (rho, V, A, theta), `"rocket"` (M0, mdot, Ve, t); `faces` = list of dict(name, area, un_rel, mass_flux, momentum_flux_x); `result` = drag per span / wave speed / plate force / rocket speed | E1 scenarios | mixed SI | C04 figure, E1 parity | §8 (keys fixed here) |
| 11.11 | `stream_tube_element_balance_sym() -> dict(mass, momentum, dropped, result)` (sympy re-run of Ex. 4.2, D06) | Ex. 4.2 | sympy | C05 (D06 check) | §4 #15 |
| 11.12 | `orifice_mass_flow(h, A, rho=1000.0, Cc=1.0, g=9.81) -> float` · `tank_drain(h0, A_tank, A_orifice, Cc=1.0, g=9.81, t_eval=None) -> dict(t, h, t_empty)` | N105 | kg/s, s, m | C05 | §4 #49 |
| 11.13 | `sprinkler_torque(a, rho, A, U, alpha) -> float` (2aρAU² cos α) · `sprinkler_free_spin_rate(a, U, alpha) -> float` (U cos α / a, our extension) | Ex. 4.6 (N84) | N m, rad/s | C04 | §4 #41 (+ free spin **NEW**) |
| 11.14 | `cube_spin_acceleration(tau12, tau21, rho, h) -> float` (6(τ₁₂ − τ₂₁)/(ρh²)) | D08 | rad/s² | C07, E3 | §4 #22 |
| 11.15 | `stokes_first_problem(y, t, U, nu) -> u` (U erfc(y/(2√(νt)))) · `plane_poiseuille(y, G, h, mu) -> u` (G y(h − y)/(2μ)) | C08 figure, C15 collapse | m/s | C08, C15, E4 | §8 (+ `plane_poiseuille` **NEW**) |
| 11.16 | `coriolis_projectile(u0, Omega, t_eval) -> dict(t, forward, deflection, deflection_small, angle, rotating, inertial)` (`deflection` = u₀t sin Ωt exact, `deflection_small` = Ωu₀t² (D17), `angle` = Ωt) · `rotating_pump_equations() -> list[sympy.Eq]` (Ex. 4.5 with the rotation terms) · `high_low_flow(x, y, U_R=1.0, sense="high") -> (u, v)` (radial outflow from a high / inflow to a low, for the N62 quiver) | N61, D17, N65, N62 | m, rad | C09, E5 | §4 #32, #34 (+ `high_low_flow` **NEW**) |
| 11.17 | `couette_heating(y, U, h, mu, k, T0=293.15, dpdx=0.0) -> dict(u, dudy, eps, T, q_bottom, q_top, work_in, heat_out, dT_max)` (steady; walls y = 0 fixed, y = h moving at U, both at T0; `eps` per unit mass needs `rho` — returns ρε [W/m³] as `eps`) · `couette_heating_transient(y, t, U, h, mu, k, rho, cp, T0=293.15, nterms=200) -> T` (sine series from T = T0) | C10, N79 | m/s, W/m³, K, W/m² | C10, E7 parity | §8 (keys fixed; transient **NEW** name) |
| 11.18 | `rankine_bernoulli(r, Gamma, sigma, rho=1000.0, p_inf=0.0) -> dict(u_theta, p, B)` (pressure from the radial balance dp/dr = ρu_θ²/r, B = ½u_θ² + p/ρ) | C11 figure | m/s, Pa, m²/s² | C11, E6 | **NEW** (C11's "B along vs across circles" needs it) |
| 11.19 | `u_tube_column(t, L, h0, g=9.81, rho=1000.0) -> dict(h, U, dUdt, dp)` (frictionless column of length L: ω² = 2g/L, h = h₀ cos ωt, dp = ρL dU/dt) · `accelerating_sphere_pressure(theta, a, dUdt, rho=1000.0) -> dict(p, force, added_mass)` (unsteady part ρa cos θ (dU/dt)/2; force −½ρV dU/dt) | C12 | m, m/s, Pa, N, kg | C12, E6 | §8 |
| 11.20 | `boussinesq_validity(alpha, dT, L, U, c=340.0, g=9.81, nu=1.5e-5, cp=1004.5) -> dict(alpha_dT, H_c, L_over_Hc, heating_ratio, g_prime, valid, verdict)` · `boussinesq_scenario(name) -> dict(fluid, alpha, dT, L, U, nu, cp, c, T0)` (`"lake"`, `"thermocline"`, `"lab_tank"`, `"abl"`, `"deep_atmosphere"`) · `boussinesq_density(T, rho0, alpha, T0) -> float` · `gaussian_blob_advection_diffusion(x, t, U, kappa, sigma0, amp=1.0, dim=2) -> T′` (σ² = σ₀² + 2κt, amplitude (σ₀/σ)^dim) · `blob_rise(t, g_prime, tau_d=10.0) -> dict(w, z)` (our one-line parcel model for E8, labelled as such: dw/dt = g′ − w/τ_d, w = g′τ_d(1 − e^{−t/τ_d}), z = g′τ_d(t − τ_d(1 − e^{−t/τ_d}))) | N108, N112, N113, N114 | –, m, m/s², kg/m³, K | C13, E8 | §4 #51, #54, #55, §8 |
| 11.21 | `two_layer_conduction(y, k1, k2, L1, L2, T_hot, T_cold) -> dict(T, q, T_interface)` · `two_fluid_couette(y, mu1, mu2, h1, h2, U) -> dict(u, tau, u_interface)` · `navier_slip_couette(y, U, h, slip_length) -> u` (our extension) | N116, N117 | K, W/m², m/s, Pa | C14 | §4 #56 (+ slip **NEW**) |
| 11.22 | `linear_wave_surface(x, z, t, a, k, g=9.81, H=np.inf) -> dict(eta, u, w, omega)` (test field only; Ch. 7 derives it) · `wave_kinematic_residual(x, t, a, k, g=9.81, H=np.inf, full=True) -> float` (Dη/Dt evaluated at the surface z = η for the linear solution; `full=False` the linearised condition at z = 0) | C14, B1 | m, m/s | C14, B1 parity | §8 (+ `wave_kinematic_residual` **NEW**) |
| 11.23 | `spheroid_area(volume, aspect) -> float` · `meniscus_height(theta, sigma, rho, g=9.81) -> float` · `meniscus_profile_x(zeta, theta, sigma, rho, g=9.81) -> x` · `meniscus_profile_ode(theta, x_max, sigma, rho, g=9.81, n=200) -> (x, zeta)` | N123, N129 | m², m | C14 | §4 #59, #62 |
| 11.24 | `prandtl_of(fluid="air", T=293.15) -> float` · `ship_drag_extrapolation(L_p, U_p, S_p, scale, D_m_total, CDf_m, CDf_p, rho_m=1000.0, rho_p=1025.0, g=9.81) -> dict(U_m, D_m_friction, D_m_wave, D_p_wave, D_p_friction, D_p_total, Re_ratio)` | N152, N156 | –, m/s, N | C15, E9 | §4 #69, #70 |

### C.12 `scripts/` (drawing helpers, no physics; `from scripts.ch04_drawings import …`)
`ch04_drawings.py`: `cv_box(ax, x0, x1, y0, y1, b=0.0, label=…)` (a control volume with outward-normal ticks),
`face_flux_arrows(ax, faces, colors)` (per-face mass/momentum arrows coloured by sign: in blue, out orange),
`budget_bars(ax, items, total_label)` (waterfall bars with hatched negatives), `stream_tube_element(ax3d, …)` (Fig. 4.3
analogue with the side-pressure arrow), `rotating_frames(ax, t, Omega)` (inertial and rotating axes, Fig. 4.6
analogue), `cone_sweep(ax3d, alpha, Omega, dt)` (Fig. 4.7 analogue), `pillbox(ax, l)` (Fig. 4.18 analogue),
`curved_cap(ax3d, R1, R2, zeta)` (Fig. 4.19 analogue). Figures → `outputs/ch04/`.

**Count:** C.1 7 rows (6 functions + 3 dataclasses) · C.2 6 rows (14 functions + 1 preset table) · C.3 7 rows (16) · C.4 17 rows (36 + 3 NamedTuples) · C.5 7 rows (11 + 2 constants) · C.6 1 row (7) · C.7 8 rows (21 + 1 table) · C.8 4 rows (9) · C.9 5 rows (29 + the `Scales` dataclass) · C.10 1 · C.11 24 rows (50) · C.12 1 script (8 drawing helpers) — **≈ 200 public functions** in 9 new core modules + the chapter module, plus 7 result classes and 4 constants/tables. **Beyond analysis §4 and curation §8 (NEW, flagged):**
`interval_mass_budget` (1.6), `flux_along_path` and `stream_function_pair` and `velocity_preset` (C.2),
`deviatoric_part` (3.7), `continuity_terms` as the (4.7) split with the old D/Dt split renamed
`continuity_material_terms` (4.2–4.3), the `per` switch of `ns_incompressible_terms`, `exact_field`, `"solid_body"`,
`ns_terms_preset` (4.9), `divergence_first_index_demo` (4.17), `basis_rate` (5.2), `gauge_absorbed_bracket` (7.5),
`nondimensional_ns_sym` (9.2), `closure_count("empty")`, `material_interval`, `ball_integral`,
`sprinkler_free_spin_rate`, `plane_poiseuille`, `high_low_flow`, `couette_heating_transient`, `rankine_bernoulli`,
`navier_slip_couette`, `wave_kinematic_residual`, `blob_rise`, `which_bernoulli_text`, the `form` switch of `ns_terms_preset`. **Signatures fixed here where the plans left them open:**
`MassBudget`/`MomentumBudget` field names, `stress_on_plane` (angle θ, signed shear), `cv_scenario` keys,
`couette_heating` keys, `flux_between_streamlines` orientation, `frame_acceleration_terms` vs `apparent_body_forces`
signs, the preset-name alternative to callables everywhere a parity row needs one.

**Contract details the notebook and explainers rely on (please keep):**
`interval_mass_budget(1.0, 2.0, 0.0)` (expanding, a = 1, ρ₀ = 1) → storage −1.0, flux_left +1.0, flux_right +2.0,
net_outflux +1.0, residual 0 · with `dx0dt=1.0, dx1dt=2.0` (material ends) → all fluxes 0, storage 0 ·
`material_mass(1.0, 2.0, 1.0)` = 1.0 kg/m² · `wake_drag_per_span("gaussian", 10.0, 1.2, 2.0)` = 3.6523 N/m (H = 2 m ≫ b)
and `wake_side_outflow` = 0.35449 m²/s · `bore_speed(1.0, 1.1)` = 3.3661 m/s · `sprinkler_torque(0.2, 1000.0, 1e-4,
5.0, np.pi/6)` = 0.8660 N m · `pitot_speed(500.0, 0.0, 1.2)` = 28.868 m/s · `torricelli_speed(1.0)` = 4.4294 m/s ·
`newtonian_stress(G_shear(10), mu=1e-3)[0, 1]` = 0.010 Pa · `stress_on_plane(G, 0, 1e-3, 0, np.pi/4)` for the shear
G = [[0, 10, 0], [0, 0, 0], [0, 0, 0]] → (0.010, 0.0) Pa · `cube_spin_acceleration(1.0, 0.0, 1000.0, 0.01)` = 60.0
rad/s² · `ns_terms_preset("poiseuille", 0.0, 2.5e-4, component=0, per="volume", G=100.0, h=1e-3, mu=1e-3)` → pressure
+100.0, viscous −100.0, local 0, advective 0, residual 0 (N/m³) · `coriolis_projectile(10.0, OMEGA_EARTH, [3600.0])` →
forward 36 000 m, deflection_small 9450.6 m, deflection 9342.4 m, angle 0.26252 rad (15.04°) · `centrifugal_acceleration` at the
equator 0.033916 m/s² · `dissipation_rate(G_shear(1000), 1000.0, 1e-3)` = 1.0 W/kg · `couette_heating(…, U=1, h=1e-3,
mu=1e-3, k=0.6)["dT_max"]` = 2.0833e-4 K and `heat_out` = 1.0 W/m² · `rankine_bernoulli(0.0, 2π, 1.0)["B"]` = −1.0,
`(2.0, …)["B"]` = 0.0 m²/s² · `which_bernoulli(True, False, True, True)` = ["4.19", "4.71", "4.72", "4.75"] when
constant_density=True, else ["4.71", "4.72", "4.75"] · `stagnation_temperature(300.0, 100.0)` = 304.98 K ·
`boussinesq_validity(2e-4, 10.0, 10.0, 0.1, c=1500.0, nu=1e-6, cp=4186.0)["alpha_dT"]` = 2.0e-3 and `g_prime` =
0.01962 m/s² · `capillary_length(0.0728, 998.0)` = 2.7269 mm · `laplace_jump_from_balance(0.0728, 1e-3, 1e-3)` = 145.6
Pa · `reynolds_number(1.0, 0.01, 1e-6)` = 1.0e4 · `sphere_drag_coefficient(1e4)` = 0.3926 (Morrison) ·
`froude_scaled_speed(10.0, 100.0, 4.0)` = 2.0 m/s · `richardson_number(9.81e-3, 100.0, 0.1)` = 98.1.

**Implementation cross-check (probe of the implementer's WIP, 2026-09-23).** `fluidpy/ch04_conservation_laws.py` (1212
lines) and the nine new core modules already define 201 of the 205 Part C names, and 32 of the 33 contract numbers above
evaluate as written (interval budgets, wake, bore, sprinkler, pitot, Torricelli, τ₁₂, `stress_on_plane`, cube spin,
`ns_terms_preset` including `form="curl"`, projectile, Coriolis term vs force signs, dissipation, Couette heating,
Rankine B, stagnation temperature, Boussinesq numbers, capillary length, Laplace, Re, Morrison C_D, Froude speed, Ri,
`which_bernoulli`, the gauge sign). **Three gaps for the implementer:** (1) `rossby_number(U, Omega, l)` returns U/(Ωl);
this design fixes it to **U/(2Ωl)** — the ratio of the advective term U²/l to the Coriolis term 2ΩU of (4.45), equal to
Ch. 13's U/(fL) at the pole (the curation's "U/(Ωl)" drops the 2; E9's parity row and Explain use U/(2Ωl)); (2)
`blob_rise(t, g_prime, tau_d)` (C.11 row 11.20, E8) and (3) `which_bernoulli_text(...)` (C.7 row 7.7, E6's exact-text
parity rows) are not written yet — a few lines each.

---

## Part A — notebook storyboard (`notebooks/build_ch04.py` → `notebooks/ch04_conservation_laws.ipynb`)

**Section check (book order; one `nb.section` per book section; cell numbers are estimates ± 5).**
- §4.1 → no A item (curation §3): notation table, N01 (tagged → C01), N02 ledger skeleton (tagged → C06) — cells 7–12
- §4.2 → R01, C01 (N03 N04 N05 N06 · D01), R02, C02 (N07–N12 · D02 D03) — cells 13–72
- §4.3 → C03 (N13–N18 · D04 · **E2**) — cells 73–100
- §4.4 → R03 R04, C04 (N19–N29 N31 N32 N82 N83 N84 · D05 · **E1**), C05 (N30 N103 N104 N105 · D06), C06 (N02 N33–N37 · D07) — cells 101–186
- §4.5 → R05 R06, C07 (N38–N51 · D08 D09 D10 · **E3**) — cells 187–232
- §4.6 → C08 (N52–N55 · D11 D12 D13 · **E4**) — cells 233–268
- §4.7 → C09 (N56–N65 · D14 D15 D16 D17 D18 · **E5**) — cells 269–320
- §4.8 → R07 R08, C10 (N66–N81 · D19 D20 D21 D22 D23 · **E7**) — cells 321–372
- §4.9 → pointers to N82–N84 (C04) and N103–N105 (C05); C11 (N85–N92 N94–N96 N102 · D24 D25); R09, C12 (N93
  N97–N101 · D26 · **E6**); R10, C13 (N106–N114 · D27 D28 · **E8**) — cells 373–450
- §4.10 → R11, C14 (N115–N129 · D29) — cells 451–490
- §4.11 → R12, C15 (N130–N156 · D30 · **E9**) — cells 491–540
- end → S01, S02, summary — cells 541–544

Every CORE block below has the nine parts in order (question · idea · primers · maths with its D rows · tiny example ·
code + "What does the code above do?" · from-scratch check where curation §7 asks · ≥ 1 visual · notes + "What would
change if…"). *code:* is the intent (exact `ch04.` calls and arguments; every line gets a novice comment in the
builder); *expect:* the numbers the executed cell must print (computed by hand here — the builder re-checks them against
the implemented `fluidpy`); *explain:* the numbered "What does the code above do?" list. Markdown drafts are in our words;
every equation is written out in LaTeX next to its number (convention 11). Derivation cells are copied from Part F word
for word. B1 (`kinematic_free_surface`) is the backup and is **not** embedded (9 explainers embedded; 5–10 allowed).

### A.0 Front matter
1. `nb.title(big_idea="Every flow obeys three bookkeeping rules: mass is neither made nor lost, momentum changes only
   when a force acts (Newton), and energy is conserved (the first law). Written for a fixed set of fluid particles and
   carried to any control volume by the Reynolds transport theorem, they become integral budgets that give forces from
   fluxes; shrunk to a point, they become the continuity equation $\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0$
   *(Eq. 4.7)*, Cauchy's equation $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ *(Eq. 4.24)* and the
   energy equation. A Newtonian stress law turns Cauchy's equation into Navier–Stokes,
   $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ *(Eq. 4.39b)* — the equation of the rest of the book.
   The chapter then shows what a rotating observer adds (Coriolis), how viscosity turns motion into heat, when
   Bernoulli's equation holds, why the ocean and atmosphere can ignore density except where it meets gravity
   (Boussinesq), what happens at the edges of a flow, and which dimensionless numbers decide whether two flows are the
   same flow.", roadmap=["§4.1 integral vs differential laws; counting equations and unknowns", "§4.2 mass: a moving
   control volume (4.5) and the continuity equation (4.7)", "§4.3 stream functions: one scalar carries a 2-D flow",
   "§4.4 momentum: forces from fluxes (4.17), Bernoulli along a streamline (4.19), Cauchy's equation (4.24)", "§4.5
   the Newtonian stress law (4.31)/(4.37)", "§4.6 Navier–Stokes (4.38)–(4.41)", "§4.7 a rotating, accelerating
   observer: Coriolis and centrifugal terms (4.45)", "§4.8 energy: dissipation turns kinetic energy into heat (4.57),
   (4.58)", "§4.9 Bernoulli's four forms, unsteady Bernoulli, Boussinesq", "§4.10 boundary and interface conditions,
   the kinematic condition (4.91)", "§4.11 dimensionless Navier–Stokes (4.101) and dynamic similarity"],
   prerequisites=["the Reynolds transport theorem and material volumes (Ch. 3 §3.6)", "the material derivative
   $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$ (Ch. 3 §3.2)", "strain rate S, vorticity ω and
   ∇·u as the volume growth rate (Ch. 3 §3.4)", "stress tensor, Cauchy traction $f_j=n_i\tau_{ij}$, Gauss' theorem,
   ε–δ identity (Ch. 2)", "Newton's viscosity law, Fourier's law, hydrostatics, first law and Gibbs relation, N²,
   Π groups (Ch. 1)"])`
2. `nb.explainer_index([("control_volume_budgets", "Weigh a force by counting what flows through a box", "C01 C04:
   the moving-CV mass and momentum budgets (4.5), (4.17) — wake drag, bore, jet, rocket"), ("stream_function_spacing",
   "Can one number field hold a whole 2-D flow?", "C03: ψ contours are streamlines, their spacing is the speed, Δψ is
   the flux"), ("newtonian_stress_lab", "How does a fluid decide its stress?", "C07: G → S → τ (4.31), traction on a
   plane, bulk viscosity, the spinning cube"), ("navier_stokes_term_balance", "Which terms of Navier–Stokes are awake
   here?", "C08: the five terms of (4.39b) at a probe for six exact solutions"), ("rotating_frame_coriolis", "Why
   does a straight throw curve on a merry-go-round?", "C09: (4.43), (4.45), Coriolis deflection Ωut², effective
   gravity"), ("which_bernoulli", "Bernoulli is constant along what, exactly?", "C05 C11 C12: the four Bernoulli
   equations and their hypotheses"), ("viscous_dissipation_heating", "Where does the energy go when viscosity stops a
   flow?", "C10: ε ≥ 0 (4.58), heating of a sheared channel, entropy production"), ("boussinesq_buoyancy", "Density
   hardly changes — so why does it drive the flow?", "C13: (4.86), (4.89) and when they hold"),
   ("dynamic_similarity_models", "When does a model behave like the real thing?", "C15: (4.101), Re, Fr, St and model
   testing")])` (the index text is plain; the explainers show every equation in full).
3. `nb.setup()`, then `nb.code` (no CORE yet) — *code:* `import numpy as np; import sympy as sp; import
   matplotlib.pyplot as plt; import plotly.graph_objects as go` · `from fluidpy import ch04_conservation_laws as ch04` ·
   `from fluidpy import ch01_introduction as ch01, ch03_kinematics as ch03` · `from fluidpy.core.interact import
   slider_figure` · `from fluidpy.core.anim import animate` · `from fluidpy.core.style import COLORS, savefig` ·
   `from fluidpy.core import transport, grids, operators, tensors, kinematics, statics, integral_theorems as itg` (reused Ch. 2–3 machinery, called by module name) · `from tools.convergence import observed_order`. *expect:* no output. *explain:* 1. numerical, symbolic and plotting
   libraries; 2. the chapter module — every function of this notebook lives there (or in the core modules it
   re-exports) and is tested in `tests/test_ch04.py`; 3. Ch. 1 and Ch. 3 functions reused as test fields; 4. the house
   helpers for sliders, animations and colours.
4. `nb.md` — "**Notation, conventions and colours used in this notebook.**" Table | symbol | meaning | unit |: ρ density
   [kg/m³]; u = (u, v, w) velocity [m/s]; b velocity of a control surface [m/s]; V(t), A(t) material volume and its
   surface; V*(t), A*(t) control volume and its surface; n outward unit normal; g body force per unit mass [m/s²], z up,
   g = −g e_z, Φ = gz; τ_ij stress [Pa] (first index = face normal), p pressure, σ_ij viscous stress; S_ij strain rate,
   ω vorticity [1/s]; μ, ν, μ_v viscosities; e internal energy, h enthalpy [J/kg]; q heat flux [W/m²]; k conductivity;
   ε dissipation [W/kg]; Ω rotation of a frame [rad/s]. Then two ⚠️ lines: "**Four primes**: x′ in §4.7 is the rotating
   frame, p′ in (4.67) is an integration variable, p′ and ρ′ in §4.9 are departures from the hydrostatic state, and
   Ch. 3's x′ was a translating frame." · "**Symbols reused**: σ (viscous stress / surface tension / vortex core
   radius), ε (dissipation / ε_ijk), γ (a coefficient in (4.29) / C_p/C_v / Ex. 4.7's ζ/δ) — so we write the shear rate
   $\dot\gamma$ and an extension rate s; Φ, Ψ, φ, Ω, h, H, M, R, f each have two or more meanings (convention 8)."
   Colours (convention 10) as a one-line legend.
5. `nb.md` — "**Where this chapter is used later**" (two columns): continuity (4.7) → every flow of Ch. 5–16 · moving-CV
   budgets (4.5), (4.17) → momentum integral (Ch. 9), layer budgets (Ch. 13), shocks (Ch. 15) ·
   $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ (4.39b) → Ch. 5–13 · the Coriolis term −2Ω × u′ of
   (4.45) → geostrophy, Ekman layers, Rossby waves (Ch. 13) · dissipation ε (4.58) → turbulence (Ch. 12) ·
   $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int\frac{dp}{\rho}+gz=$ const (4.75) → waves (Ch. 7) ·
   Boussinesq (4.86), (4.89) → internal waves (Ch. 7), convection (Ch. 11), ocean and atmosphere (Ch. 13) · Dη/Dt = 0
   (4.91) → every free surface · Re, Fr, Ri, Ro → the regime maps of every later chapter. Climate hook: "the equations of
   an ocean or atmosphere model are (4.7), (4.45) with the Boussinesq simplification (4.86), and (4.89) — this chapter
   writes each one down and says when it holds."

### A.1 §4.1 Introduction — no A item (N01 tagged → C01; N02 tagged → C06)
1. `nb.section("4.1", "Introduction", intro="**What is this section about?** The chapter's plan in two sentences: the
   same three principles (mass, momentum, energy) are written twice — once as a budget for a finite volume (integral
   form) and once at a single point (differential form) — and we keep a running count of equations against unknowns
   until the count closes.")`
2. `nb.note` — **N01 [C]** "**Integral vs differential form.** An *integral* law is a budget for a finite region — what
   is inside, what flows through its surface, what forces act on it. A *differential* law is the same statement at one
   point, a field equation. Gauss' theorem plus 'true for every volume' turns the first into the second (C02 does it for
   mass). Integral forms return in Ch. 9 (momentum integral), Ch. 13 (layer budgets) and Ch. 15 (shock control volumes)."
3. `nb.md` — **The equation ledger (N02, filled in at C06, C08, C10).** "A set of equations can determine a flow only if
   there are as many independent equations as unknown fields. We count as we go:"
4. `nb.code` — *code:* `import pandas as pd` · `stages = ["cauchy", "navier_stokes", "barotropic", "full"]` ·
   `ledger = pd.DataFrame({s: ch04.closure_count("empty") for s in stages}).T` (zeros: nothing derived yet) ·
   `print(ledger[["equations", "unknowns"]])`. *expect:* a 4 × 2 table of zeros with the four stage names. *explain:*
   1. `closure_count("empty")` returns the skeleton; 2. the same call with a stage name after C06, C08 and C10 fills a
   row — the table closes (7 = 7) at C10. (`pandas.DataFrame` gloss in the comment: "a labelled table".)

### A.2 §4.2 Conservation of Mass — R01, C01 (+N03–N06 · D01), R02, C02 (+N07–N12 · D02, D03)
1. `nb.section("4.2", "Conservation of Mass", intro="**What is this section about?** Mass can be neither created nor
   destroyed. Written for a sealed bag of fluid that moves with the flow, this is one line; the work is to rewrite it
   for any box we choose (moving, deforming or fixed) and then for a single point.")`
2. `nb.recap("R01", "Material volume and material surface", "A *material volume* V(t) always contains the same
   fluid particles; its surface A(t) moves with the local fluid velocity, so in the Reynolds transport theorem its
   surface velocity is b = u. Picture a perfectly sealed, infinitely flexible balloon: it drifts, stretches and
   twists, but nothing crosses its skin. The theorem we use, $\frac{d}{dt}\int_{V^*(t)}F\,dV=\int_{V^*(t)}\frac{\partial F}{\partial t}dV+\int_{A^*(t)}F\,\mathbf b\cdot\mathbf n\,dA$
   *(Eq. 3.35)*, holds for any control volume; for the balloon, b = u.", where="Ch. 3 §3.6")` — followed by `nb.code`
   (still before C01; allowed outside a core): *code:* `u = lambda x, t: ch04.expanding_flow(x, t, a=1.0)[0]` ·
   `ends = [ch03.pathline(u, [X], 0.0, [0.0, 1.0])[0, -1] for X in (1.0, 2.0)]` · `print(ends)`. *expect:* `[2.0, 4.0]`
   (to 1e-8). *explain:* "two particles that start at x = 1 m and 2 m end at 2 m and 4 m after 1 s: the material
   interval stretches as 1 + at."

**C01 — Mass conservation for an arbitrarily moving control volume (4.5)**
3. `nb.core("C01", "Mass in a box that moves the way you choose (4.5)", question="A sealed balloon of air drifts and
   swells, so its mass never changes. A fixed box sits where the balloon passes. What does the box record — and what
   would a box that slides along at its own speed record?")`
4. `nb.md` — **The problem in plain words:** "Engineers rarely follow a particular lump of fluid: they draw a box round a
   pump, a wing or a river reach and count what enters and leaves. Meteorologists do the same with a grid cell, and
   oceanographers with a layer between two density surfaces that moves up and down. We need one mass-conservation
   statement that works for *any* box — fixed, sliding, stretching or riding with the fluid — and that reduces to 'the
   balloon keeps its mass' when the box *is* the balloon."
5. `nb.md` — **The idea:**
   ```
   material volume V(t)          control volume V*(t)              what the box sees
   moves with the fluid (b = u)  moves as you choose (b = anything)  d/dt(mass inside) = −(net mass flux through
   mass fixed: d/dt ∫ρ dV = 0    mass can change                                         its walls, RELATIVE to them)
                    the trick:  make V* and V coincide at one instant  →  (4.5)
   ```
   "The flux through a wall depends on how fast the fluid moves **relative to the wall**, u − b: a wall that moves with
   the fluid lets nothing through, however fast both go."
6. `nb.primer("dataclasses and named results", "A `dataclass` is a small Python class whose only job is to hold named
   fields; fluidpy returns every budget this way, so you write `budget.residual` instead of remembering that the
   residual is the fourth number. `@dataclass(frozen=True)` makes the result read-only (a result should not be edited
   by accident).", code="from dataclasses import dataclass           # the decorator that builds the class\n@dataclass(frozen=True)\nclass Budget:                              # three named numbers\n    storage: float; outflux: float; residual: float\nb = Budget(storage=-1.0, outflux=1.0, residual=0.0)\nprint(b.residual, b)                         # 0.0 Budget(storage=-1.0, outflux=1.0, residual=0.0)")`
   — *expect:* `0.0 Budget(storage=-1.0, outflux=1.0, residual=0.0)`.
7. `nb.note` — **N03 [B]** "**Mass of a material volume is constant** — our starting line:" equation
   `\frac{d}{dt}\int_{V(t)}\rho(\mathbf x,t)\,dV=0`, ref "4.1". "Number with our test flow $u=ax/(1+at)$,
   $\rho=\rho_0/(1+at)$ (a = 1 s⁻¹, ρ₀ = 1 kg/m³, a 1-D flow, mass per unit cross-section area): the particles that
   occupy [1, 2] m at t = 0 occupy [2, 4] m at t = 1 s, where ρ = 0.5 kg/m³ — mass 1 × 1 = 0.5 × 2 = 1 kg/m² both
   times."
8. `nb.derivation("D01", …)` — Part F D01 (9 steps), ref "4.5".
9. `nb.note` — **N04 [B]** "Step 2 of the derivation is the book's (4.2) — the transport theorem with F = ρ and b = u:"
   equation `\int_{V(t)}\frac{\partial\rho(\mathbf x,t)}{\partial t}dV+\int_{A(t)}\rho(\mathbf x,t)\,\mathbf u(\mathbf x,t)\cdot\mathbf n\,dA=0`, ref "4.2".
   "Correct but awkward: V(t) is carried by the flow, which is usually what we do not know."
10. `nb.note` — **N05 [B]** "Step 3 is (4.3), the same theorem for a box of our choosing:" equation
    `\frac{d}{dt}\int_{V^*(t)}\rho\,dV-\int_{V^*(t)}\frac{\partial\rho}{\partial t}dV-\int_{A^*(t)}\rho\,\mathbf b\cdot\mathbf n\,dA=0`,
    ref "4.3".
11. `nb.note` — **N06 [B]** "Steps 4–6 are (4.4), the **instantaneously coincident control volume**: at the one instant
    when V* and V fill the same region, integrals of the same field over the same region are equal —" equation
    `\int_{V^*}\frac{\partial\rho}{\partial t}dV=\int_{V}\frac{\partial\rho}{\partial t}dV=-\int_{A}\rho\mathbf u\cdot\mathbf n\,dA=-\int_{A^*}\rho\mathbf u\cdot\mathbf n\,dA`,
    ref "4.4". "— but the rates d/dt∫ρ dV are **not** equal: a moment later the two volumes have moved apart. That
    difference is exactly the flux term."
12. `nb.md` — "> ⚠️ **Common confusion:** '$\frac{d}{dt}\int_{V^*}\rho\,dV$ and $\int_{V^*}\frac{\partial\rho}{\partial t}dV$ are
    the same thing.' Only for a box that does not move or deform. The first is how fast the mass *in the box* changes; the
    second adds up how fast the density changes *at each fixed point*. For a moving box they differ by
    $\int_{A^*}\rho\,\mathbf b\cdot\mathbf n\,dA$ — the mass the moving walls sweep in or out (Ch. 3's signed swept volume)."
13. `nb.worked_example("the expanding flow and three boxes, all at t = 0", "Flow $u=ax/(1+at)$, $\rho=\rho_0/(1+at)$ with
    a = 1 s⁻¹, ρ₀ = 1 kg/m³; box [1, 2] m. 1. **Fixed box** (b = 0): mass flux in at x = 1: ρu = 1 × 1 = 1 kg/(m² s);
    out at x = 2: 1 × 2 = 2 kg/(m² s); net outflux 2 − 1 = +1. By (4.5), $\frac{d}{dt}\int_{V^*}\rho\,dV+\int_{A^*}\rho(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=0$,
    the mass inside changes at −1 kg/(m² s). Check: ∂ρ/∂t = −aρ₀/(1 + at)² = −1 kg/(m³ s) everywhere, times the length
    1 m = −1 ✓. 2. **Material box** (ends move at the local u: ẋ₀ = 1, ẋ₁ = 2 m/s): u − b = 0 at both ends, no flux, mass
    constant — that is (4.1). 3. **Box sliding at b = 1 m/s**: at x = 1, ρ(u − b) = 0; at x = 2, ρ(u − b) = 1 × (2 − 1) = 1;
    net outflux +1, so the mass inside falls at 1 kg/(m² s) — but for a new reason: its right wall moves slower than the
    fluid there.")`
14. `nb.code` — *code:* `print(ch04.material_interval(1.0, 2.0, 1.0))` · `print(ch04.material_mass(1.0, 2.0, 0.0),
    ch04.material_mass(1.0, 2.0, 1.0))` · `for name, (d0, d1) in {"fixed": (0.0, 0.0), "sliding b=1": (1.0, 1.0),
    "material": (1.0, 2.0)}.items():` `bud = ch04.interval_mass_budget(1.0, 2.0, 0.0, dx0dt=d0, dx1dt=d1, flow="expanding",
    a=1.0, rho0=1.0)`; `print(f"{name:12s} storage {bud['storage']:+.3f}  in {bud['flux_left']:+.3f}  out
    {bud['flux_right']:+.3f}  residual {bud['residual']:+.1e}")`. *expect:* `(2.0, 4.0)` · `1.0 1.0` · fixed: storage
    −1.000, in +1.000, out +2.000 · sliding: storage −1.000, in +0.000, out +1.000 · material: storage +0.000, in +0.000,
    out +0.000 · residuals ≤ 1e-10. *explain:* 1. where the particles of [1, 2] m are after 1 s; 2. their mass per unit
    area at t = 0 and t = 1 s (unchanged: (4.1)); 3. the budget (4.5) for three boxes on the same flow at the same instant —
    storage is measured by differencing the mass in the box in time, the fluxes by ρ(u − b) at the two ends; in every case
    storage + net outflux = 0.
15. `nb.check_agree` — **from scratch (curation §7):** midpoint sums: `x = np.linspace(1, 2, 2001)`;
    `M = lambda t, x0, x1: np.trapezoid(ch04.expanding_flow(np.linspace(x0, x1, 2001), t)[1], np.linspace(x0, x1, 2001))`;
    `dt = 1e-5`; `storage_fixed = (M(dt, 1, 2) - M(-dt, 1, 2))/(2*dt)`; `assert np.allclose(storage_fixed, -1.0,
    rtol=1e-6)` · sliding box: `(M(dt, 1+dt, 2+dt) - M(-dt, 1-dt, 2-dt))/(2*dt)` → −1 · material box:
    `(M(dt, *ch04.material_interval(1, 2, dt)) - M(-dt, *ch04.material_interval(1, 2, -dt)))/(2*dt)` → 0 · then the 3-D
    library check on a `ch03`-style box: `cv = transport.MovingBox((1.0, 1.0, 1.0), origin=(1.0, 1.0, 1.0))` (Ch. 3's box), `rho3 = lambda x, t: ch04.expanding_flow(x, t, dim=3)[1]`, `u3 = lambda x, t:
    ch04.expanding_flow(x, t, dim=3)[0]`; `B = ch04.mass_budget(rho3, u3, cv, 0.0)`; `assert abs(B.residual) < 1e-6` and
    `ch04.mass_budget(rho3, u3, cv, 0.0, material=True).outflux` ≈ 0. Markdown: "Two routes, one budget: our own
    trapezoid sums of the mass in each box, and the library's 3-D quadrature."
16. `nb.animation` — (`player="frames"`, 24 frames, FAST 12) two rows sharing t ∈ [0, 1.2] s: top, the x-axis 0–5 m with
    the density as a shaded strip (darker = denser, fading as 1/(1 + t)), the **material interval** (purple bracket) whose
    ends ride on particles, and the **fixed box** [1, 2] m (grey outline) with blue inflow and orange outflow arrows sized
    by ρu; bottom, bars for the fixed box: storage (purple), net outflux (orange), sum (black ≈ 0); a text line "mass in
    the material interval = 1.000 kg/m²". Followed by `nb.md` see/read/change: *see:* "the bracket stretches to the
    right while the strip pales; the grey box keeps its place and its bars keep adding to zero." *read:* "the balloon
    (bracket) keeps its mass because it grows exactly as fast as the density falls; the fixed box loses mass because
    more leaves through its right face than enters on the left." *change:* "…a = −0.5 s⁻¹ (a converging flow): the
    bracket shrinks, the strip darkens, and the fixed box *gains* mass — its bars flip sign together."
17. `nb.md` — **What would change if…** "…the box shrank to a point? The storage term would become ∂ρ/∂t times a tiny
    volume, the flux term a divergence times the same volume — dividing by the volume leaves one equation at every point.
    That is C02."

**C02 — The continuity equation (4.7)**
18. `nb.recap("R02", "Gauss' divergence theorem", "For a smooth vector field Q in a volume V with outward normal n on
    its surface A, $\int_V\nabla\cdot\mathbf Q\,dV=\int_A\mathbf Q\cdot\mathbf n\,dA$ *(Eq. 2.30)*: the net outflow through the
    skin equals the sum of the sources inside. Here Q = ρu, the mass flux per unit area.", where="Ch. 2 §2.12")` —
    followed by `nb.code`: `Q = lambda X, Y, Z: np.array(ch04.expanding_flow(np.array([X, Y, Z]), 0.0, dim=3)[1] *
    ch04.expanding_flow(np.array([X, Y, Z]), 0.0, dim=3)[0])` (ρu of the 3-D expanding flow on ch02's `fn(X, Y, Z)`
    convention) · `print(itg.divergence_theorem_box(Q, [(1, 2), (1, 2), (1, 2)], n=24))`. *expect:* flux = volume
    integral = 3.0 kg/s (the divergence of ρu is 3aρ₀ = 3 at t = 0; box volume 1 m³). *explain:* "both sides of (2.30) for
    the mass flux in a unit cube."
19. `nb.core("C02", "The continuity equation: mass conservation at a point (4.7)", question="A crowd pours out of a
    square faster than people arrive. The square empties. What does 'mass is conserved' say about one single point of a
    flow?")`
20. `nb.md` — **The problem in plain words:** "A numerical weather model, a pipe-flow solver and a tidal model all store
    ρ and u at points of a grid. They need a rule that ties the density's rate of change at a point to how the flow
    spreads out or converges there — a field equation, not a box budget. Every exact solution in later chapters is
    checked against it."
21. `nb.md` — **The idea:**
   ```
   box budget (4.5), fixed box      →  Gauss: surface flux = ∫ ∇·(ρu) dV     →  ∫ [∂ρ/∂t + ∇·(ρu)] dV = 0
   true for EVERY box, however small  →  the integrand itself must vanish     →  ∂ρ/∂t + ∇·(ρu) = 0   (4.7)
   ```
   "**Flux divergence ∇·(ρu)** = the net outflow of mass per unit volume from a point: where it is positive, the density
   there must fall."
22. `nb.primer("continuity of a function and the small-ball argument", "A function f is *continuous* at x₀ if its values
    near x₀ are close to f(x₀): pick any margin, and a small enough ball around x₀ keeps f within that margin. So if
    f(x₀) = 3, some small ball has f > 1.5 everywhere inside it — and the integral of f over that ball is then at least
    1.5 × (the ball's volume) > 0. That is the whole trick of the 'localisation lemma' below.", code="f = lambda x, y,
    z: 3.0 - 50*((x - 0.5)**2 + y**2 + z**2)   # f(0.5, 0, 0) = 3, continuous\nfor r in (0.2, 0.1, 0.05):              # shrink the ball\n    I = ch04.ball_integral(f, (0.5, 0.0, 0.0), r)\n    print(r, I, I/(4/3*np.pi*r**3))   # the mean over the ball → f(x0) = 3")`
    — *expect:* r = 0.2: I ≈ 0.0603, mean ≈ 1.80; r = 0.1: mean ≈ 2.70; r = 0.05: mean ≈ 2.925 (mean = 3 − 30r², since
    the ball-average of r'² is 3r²/5) — the mean tends to f(x₀) = 3 and stays positive once r < 0.316 m.
23. `nb.note` — **N07 [B]** "**Localisation lemma.** If $\int_Vf\,dV=0$ for **every** volume V and f is continuous, then
    f = 0 at every point. Proof in the derivation below (steps 5–7): a point where f ≠ 0 would have a small ball where f
    keeps its sign, and that ball's integral could not be zero." equation `\int_Vf\,dV=0\ \ \forall V\ \Rightarrow\ f\equiv0`,
    no ref. "We use it again for momentum (C06) and energy (C10)."
24. `nb.derivation("D02", …)` — Part F D02 (8 steps), ref "4.7".
25. `nb.note` — **N08 [B]** "**Flux divergence (transport) term.** ∇·(ρu) is the net mass outflow per unit volume at a
    point. Over a whole domain with no flux through its boundary it integrates to zero (Gauss): transport moves mass
    around, it cannot change the total." Followed by `nb.code`: a periodic 2-D box, `g = grids.grid2d(bounds=((0, 2*np.pi),
    (0, 2*np.pi)), n=(64, 64), periodic=True)` (Ch. 2's grid), ρ = 1 + 0.2 sin x cos y, u = (sin y, cos x)
    (a smooth divergence-free field); `flux_div = operators.divergence(rho*U, g.h, bc="periodic")`; `print(flux_div.sum() *
    g.h[0]*g.h[1], flux_div.min(), flux_div.max())`. *expect:* total ≈ 0 (≤ 1e-12) while min and max are of order ±0.1 — "local
    changes, no global change". *explain:* "∂ρ/∂t = −∇·(ρu) raises the density where the flux converges (negative
    divergence) and lowers it elsewhere; the domain total never changes."
26. `nb.primer("product rule for a divergence", "The divergence of a scalar times a vector splits like a product rule:
    $\nabla\cdot(\rho\mathbf u)=\mathbf u\cdot\nabla\rho+\rho\,\nabla\cdot\mathbf u$ — in index form
    $\partial(\rho u_i)/\partial x_i=u_i\,\partial\rho/\partial x_i+\rho\,\partial u_i/\partial x_i$ (the book's (B3.6)).
    Here: the flux of mass changes because the density varies along the flow *or* because the flow spreads out.",
    code="x, y, z = sp.symbols('x y z'); rho = sp.Function('rho')(x, y, z)\nu = [sp.Function(f'u{i}')(x, y, z) for i in (1, 2, 3)]; X = (x, y, z)\nlhs = sum(sp.diff(rho*u[i], X[i]) for i in range(3))              # ∇·(ρu)\nrhs = sum(u[i]*sp.diff(rho, X[i]) + rho*sp.diff(u[i], X[i]) for i in range(3))\nprint(sp.simplify(lhs - rhs))                                         # 0")`
    — *expect:* `0`.
27. `nb.derivation("D03", …)` — Part F D03 (5 steps), ref "4.10".
28. `nb.note` — **N09 [B]** "D03's step 3 is the book's (4.8), continuity written with the particle's own rate:"
    equation `\frac{1}{\rho(\mathbf x,t)}\frac{D}{Dt}\rho(\mathbf x,t)+\nabla\cdot\mathbf u(\mathbf x,t)=0`, ref "4.8". "Read it
    as: the fractional rate at which a particle's density rises equals minus the rate at which its volume grows (Ch. 3's
    $\frac{1}{\delta V}\frac{D(\delta V)}{Dt}=\nabla\cdot\mathbf u$ (3.14)). Number: in the expanding flow at t = 0,
    $(1/\rho)D\rho/Dt=-a/(1+at)=-1$ s⁻¹ and ∇·u = +1 s⁻¹."
29. `nb.note` — **N10 [B]** "**Incompressible flow** means each particle keeps its density:" equation
    `\frac{D\rho}{Dt}\equiv\frac{\partial\rho}{\partial t}+\mathbf u\cdot\nabla\rho=0`, ref "4.9". "Different particles may
    have different densities: a stratified lake flowing horizontally, $\mathbf u=(U(z),0,0)$, $\rho=\rho(z)$, has
    Dρ/Dt = 0 + U ∂ρ/∂x = 0 although ρ varies with depth. > ⚠️ §4.11 re-prints (4.9) as
    $\nabla\cdot\mathbf u=-\frac1\rho\frac{D\rho}{Dt}=-\frac{1}{\rho c^2}\frac{Dp}{Dt}$ — that is a general identity for any flow
    whose density changes are isentropic ($dp=c^2d\rho$), not the incompressibility condition (C15 uses it)."
30. `nb.note` — **N11 [B]** "With (4.9) in (4.8) the continuity equation of an incompressible flow is simply" equation
    `\nabla\cdot\mathbf u=0`, ref "4.10". "— zero volumetric strain rate (Ch. 3 §3.4): particles change shape but not
    volume. Number: the Ch. 3 cylinder flow has ∇·u = 0 to round-off at every point off the body (code below)."
31. `nb.note` — **N12 [B]** "**Incompressible flow vs incompressible fluid.** ρ = const everywhere is a special case of
    (4.9) (every particle has the *same* density). A fluid whose density does not respond to pressure is an incompressible
    *fluid*; a gas is not, but its *flow* is nearly incompressible when the speed is small compared with the speed of
    sound." Table | case | Dρ/Dt = 0? | ρ = const? |: lake with a thermocline — yes, no · tap water in a pipe — yes, yes
    · air at 100 m/s — nearly, M = 0.29 · air at 300 m/s — no, M ≈ 0.88. "Number: air at 288 K has c = 340.3 m/s, so
    U = 100 m/s gives M = U/c = 0.294 < 0.3 (the usual rule of thumb; C15 shows the departure from ∇·u = 0 scales with M²)."
32. `nb.md` — "> ⚠️ **Common confusion:** 'incompressible' means 'constant density'. It means *each particle* keeps its
    density; the ocean's density varies by a few kg/m³ and its flow is still incompressible to excellent accuracy."
33. `nb.worked_example("continuity at one point of the expanding flow (t = 0, x = 1.5 m)", "Flow $u=ax/(1+at)$,
    $\rho=\rho_0/(1+at)$, a = 1 s⁻¹, ρ₀ = 1 kg/m³. 1. Local term: $\partial\rho/\partial t=-a\rho_0/(1+at)^2=-1$ kg/(m³ s).
    2. Flux: $\rho u=\rho_0ax/(1+at)^2=1.5$ kg/(m² s); its x-derivative $\partial(\rho u)/\partial x=a\rho_0/(1+at)^2=+1$
    kg/(m³ s). 3. Sum: −1 + 1 = 0 ✓ — (4.7), $\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0$. 4. Particle form:
    $u\,\partial\rho/\partial x=0$ (ρ does not vary in x), so $D\rho/Dt=-1$, $(1/\rho)D\rho/Dt=-1$ s⁻¹, and
    $\nabla\cdot\mathbf u=a/(1+at)=+1$ s⁻¹: (4.8) holds. The flow is compressible (Dρ/Dt ≠ 0) — it is our test of the
    compressible form.")`
34. `nb.code` — *code:* `rho = lambda x, t: ch04.expanding_flow(x, t)[1]` · `u = lambda x, t: ch04.expanding_flow(x,
    t)[0]` (1-D fields on (1,) points) · `x0 = np.array([1.5])` · `print(ch04.continuity_terms(rho, u, x0, 0.0, ht=1e-4))` ·
    `print(ch04.continuity_material_terms(rho, u, x0, 0.0, ht=1e-4))` · stratified check:
    `print(ch04.density_material_rate(lambda x, t: ch04.stratified_shear_flow(x[2])[1], lambda x, t:
    np.array([ch04.stratified_shear_flow(x[2])[0], 0.0, 0.0]), np.array([0.3, 0.0, -2.0]), 0.0))` · cylinder:
    `print(ch04.divergence_free_check(ch03.cylinder_velocity_field(1.0, 1.0), np.array([1.5, 0.7]), 0.0))` · symbolic:
    `xs, ts, a = sp.symbols('x t a', positive=True)`; `print(sp.simplify(ch04.continuity_residual_sym(1/(1 + a*ts),
    [a*xs/(1 + a*ts)], [xs], ts)))` · `print(ch04.is_incompressible_regime(100.0), ch04.mach_number(100.0, T=288.15))`.
    *expect:* `ContinuityTerms(local≈-1.0, flux_divergence≈+1.0, residual≈0)` · `(-1.0, 1.0)` (to 1e-7) · `0.0` (Dρ/Dt
    of the stratified flow) · `(True, ~1e-12)` · `0` · `True 0.2939`. *explain:* 1. the two terms of (4.7) by
    second-order stencils in x and t; 2. the particle form (4.8) — fractional density rate and divergence cancel; 3. a
    flow with varying density that is still incompressible (N10); 4. ∇·u = 0 for the cylinder (N11); 5. sympy proves the
    residual is identically 0 for the test flow; 6. the Mach-number rule of thumb (N12).
35. `nb.check_agree` — **from scratch (curation §7):** hand-written central differences at x = 1.5, t = 0: `hx, ht = 1e-4,
    1e-4`; `drho_dt = (rho(x0, ht) - rho(x0, -ht))/(2*ht)`; `dflux_dx = (rho(x0+hx, 0)*u(x0+hx, 0) - rho(x0-hx,
    0)*u(x0-hx, 0))/(2*hx)`; `assert np.allclose(drho_dt + dflux_dx, 0.0, atol=1e-8)`; `assert
    np.allclose([drho_dt, dflux_dx], ch04.continuity_terms(rho, u, x0, 0.0, ht=ht)[:2], rtol=1e-8)` · a **wrong** density
    ρ = ρ₀ (constant) with the same u gives a residual of +1 — "the equation catches an impossible field".
36. `nb.plotly` — `slider_figure` over t = 0…2 s (21 steps, FAST 11): x ∈ [0, 5] m; traces "∂ρ/∂t (blue)",
    "∂(ρu)/∂x (orange)", "sum (black dashed)". Title "Two mirror images that always add to zero". Followed by see/read/
    change: *see:* "a flat blue line below zero and a flat orange line above it; both shrink toward 0 as t grows."
    *read:* "the density falls everywhere at the same rate because the flux diverges everywhere at the same rate; the sum
    is (4.7) and stays at zero." *change:* "…ρ₀ doubled: both lines double, the sum stays at zero — (4.7) is linear in
    ρ for a given u."
37. `nb.figure` — the stratified flow of N10: left, ρ(z) (blue, 1000 → 1005 kg/m³ over 10 m depth) and u(z) (teal
    shear); right, particles as dots coloured by density moving to the right for 60 s with their colours unchanged.
    *see:* "horizontal layers sliding over each other; every dot keeps its colour." *read:* "each particle keeps its
    density (Dρ/Dt = 0, incompressible) although density varies from layer to layer (ρ ≠ const)." *change:* "…a vertical
    velocity w pushed dense water up: ρ at a fixed height would rise (∂ρ/∂t > 0), but each particle's density would still
    be its own — Dρ/Dt = 0 is about particles, not places."
38. `nb.md` — **What would change if…** "…the flow were steady and two-dimensional? Then (4.7) reduces to
    $\partial(\rho u)/\partial x+\partial(\rho v)/\partial y=0$ — two unknown fields tied by one equation. One scalar can then
    satisfy it automatically and carry the whole flow picture: the stream function (C03)."

### A.3 §4.3 Stream Functions — C03 (+N13–N18 · D04 · E2)
1. `nb.section("4.3", "Stream Functions", intro="**What is this section about?** Steady continuity is a constraint on
   the velocity field. Writing the mass flux as a cross product of two gradients satisfies it identically; in two
   dimensions a single function ψ remains, its contours are the streamlines, and the flux between two of them is the
   difference of their labels.")`

**C03 — The 2-D stream function**
2. `nb.core("C03", "One number field for a whole 2-D flow: the stream function", question="On a weather map the
   streamlines crowd together where the wind is strong. Is that a drawing convention — or can one scalar field carry
   both the direction and the speed of a 2-D flow?")`
3. `nb.md` — **The problem in plain words:** "A 2-D flow has two velocity components tied by continuity. If we could
   describe the whole flow by *one* function whose contours are the streamlines, every picture would be a contour plot,
   and the continuity equation could never be violated by a numerical slip. Potential-flow theory (Ch. 6), the Blasius
   boundary layer (Ch. 9) and quasi-geostrophic ocean and atmosphere models (Ch. 13) are all written this way."
4. `nb.md` — **The idea:**
   ```
   ψ = const   along every streamline            u = ∂ψ/∂y ,  v = −∂ψ/∂x    (ρ = 1)
   ψ₂ − ψ₁   = volume flux between two streamlines (per metre of depth)
   equal steps Δψ between contours  →  contours close together = fast flow   (u ≈ Δψ/Δn)
   ```
5. `nb.note` — **N13 [B]** "Start from continuity (4.7) with ∂ρ/∂t = 0 — steady flow:" equation
   `\nabla\cdot(\rho\mathbf u)=0`, ref "4.11".
6. `nb.note` — **N14 [B]** "A divergence of a curl is always zero (∇·∇×A = 0, Ch. 2 §2.13), so any mass flux written
   as a curl satisfies (4.11) automatically. Choosing the vector potential as χ∇ψ and using ∇×∇ψ = 0 gives the book's"
   equation `\rho\mathbf u=\nabla\times\boldsymbol\Psi,\qquad\boldsymbol\Psi=\chi\nabla\psi\ \Rightarrow\ \rho\mathbf u=\nabla\chi\times\nabla\psi`,
   ref "4.12". "(That every steady solenoidal ρu can be written like this — locally — is asserted by the book, not
   proved; we only use the converse.) The sympy line below checks ∇·(∇χ × ∇ψ) = 0 for arbitrary χ, ψ."
   Followed by `nb.code`: `x, y, z = sp.symbols('x y z'); chi = sp.Function('chi')(x, y, z); psi =
   sp.Function('psi')(x, y, z)` · `m = ch04.mass_flux_from_stream_functions(chi, psi, (x, y, z))` · `print(sp.simplify(sum(sp.diff(m[i],
   v) for i, v in enumerate((x, y, z)))))`. *expect:* `0`.
7. `nb.derivation("D04", …)` — Part F D04 (8 steps), ref "".
8. `nb.md` — "> ⚠️ **Sign convention.** With χ = −z the book (and this notebook) uses $u=\partial\psi/\partial y$,
   $v=-\partial\psi/\partial x$; many meteorology and oceanography texts use $u=-\partial\psi/\partial y$, $v=\partial\psi/\partial x$
   (the flux is then ψ₁ − ψ₂). Same contours, opposite labels — check the convention before comparing numbers."
9. `nb.worked_example("the stagnation-point flow ψ = kxy with k = 1 s⁻¹", "1. Velocity: $u=\partial\psi/\partial y=kx$,
   $v=-\partial\psi/\partial x=-ky$; at (1, 2) m: u = 1 m/s, v = −2 m/s. 2. Streamlines: xy = const — hyperbolas.
   3. Flux between ψ = 1 and ψ = 2 m²/s: 2 − 1 = 1 m²/s per metre of depth, whichever path joins them. 4. Along the
   diagonal y = x the two contours cross at $x=1$ m and $x=\sqrt2=1.414$ m, 0.586 m apart along the diagonal, where
   the speed is about $kx\sqrt2\approx1.7$ m/s — flux ≈ speed × spacing = 1.7 × 0.586 ≈ 1.0 ✓. 5. Where the contours are
   0.1 m apart for Δψ = 1 m²/s, the speed is ≈ 10 m/s.")`
10. `nb.code` — *code:* `print(ch04.velocity_from_streamfunction_2d("stagnation", 1.0, 2.0, k=1.0))` ·
    `print(ch04.flux_between_streamlines("stagnation", (1.0, 1.0), (np.sqrt(2), np.sqrt(2)), k=1.0))` ·
    `print(ch04.flux_along_path("stagnation", np.array([[1.0, 1.2, 1.3, 1.414213562], [1.0, 0.9, 1.25, 1.414213562]]),
    k=1.0))` (a bent gate between the same two points) · the cylinder: `print(ch04.velocity_from_streamfunction_2d(lambda
    x, y: ch03.cylinder_streamfunction(x, y, 1.0, 1.0), 2.0, 0.5), ch03.cylinder_flow(2.0, 0.5, 1.0, 1.0))` · sympy:
    `X, Y, k = sp.symbols('x y k')`; `print(ch04.velocity_from_streamfunction_2d_sym(k*X*Y, X, Y))`. *expect:* `(1.0, -2.0)`
    · `1.0` · `1.0` · two equal pairs (≈ (0.7924, −0.1107)) · `(k*x, -k*y)`. *explain:* 1. velocity by central
    differences of ψ; 2. flux through a straight gate = ψ₂ − ψ₁; 3. the same through a bent gate — flux depends only on
    the end labels; 4. the cylinder's ψ from Ch. 3 gives back the Ch. 3 velocity; 5. the symbolic version.
11. `nb.check_agree` — **from scratch (curation §7):** `h = 1e-5; psi = lambda x, y: ch04.streamfunction_preset("stagnation",
    x, y, k=1.0)`; `u_fd = (psi(1, 2+h) - psi(1, 2-h))/(2*h); v_fd = -(psi(1+h, 2) - psi(1-h, 2))/(2*h)`;
    `assert np.allclose([u_fd, v_fd], [1.0, -2.0], rtol=1e-8)` · flux by `np.trapezoid` along the straight gate
    (normal to the right of travel): `s = np.linspace(0, 1, 401); P = np.outer(p1, 1-s) + np.outer(p2, s); t =
    (p2 - p1)/np.linalg.norm(p2 - p1); n = np.array([t[1], -t[0]])`; `un = ch04.velocity_preset("stagnation", *P,
    k=1.0)`; `flux = np.trapezoid(un[0]*n[0] + un[1]*n[1], s)*np.linalg.norm(p2 - p1)`; `assert np.isclose(flux, 1.0,
    rtol=1e-6)`.
12. `nb.figure` — ψ contours at **equal Δψ** (20 levels, black) over a speed heatmap (viridis, |u|/U) for the Ch. 3
    cylinder flow (U = 1 m/s, a = 1 m) on x ∈ [−3, 3], y ∈ [−2, 2] m; two contours ψ = 0.5 and 1.0 m²/s highlighted in
    teal with the band between them shaded; a gate segment (black) across the band with its label "flux = 0.500 m²/s".
    *see:* "contours bunch together over the top and bottom of the cylinder and spread out in front and behind; the heat
    map is brightest exactly where they bunch." *read:* "equal steps in ψ carry equal flux, so a narrow gap means fast
    flow: at the top of the cylinder the speed is 2U and the gap is about half its far-field value." *change:* "…U
    doubled: every contour label doubles (ψ ∝ U) — with the same levels the drawing would show twice as many contours,
    each gap carrying the same 0.05 m²/s at twice the speed."
13. `nb.plotly` — `slider_figure` over the source strength m = 0…4 m²/s (17 steps, FAST 9) of **source + uniform
    stream** (U = 1 m/s): contour lines of ψ at fixed levels drawn as traces (≤ 4 traces: the dividing streamline ψ = m/2
    in rose, three level sets), plus the stagnation point. Title "A source in a stream: the half-body grows as m/(2U)".
    see/read/change: "…the dividing streamline (rose) encloses the source's fluid; its width far downstream is m/U — the
    flux m spread over the stream speed."
14. `nb.plotly` — **N15 [B]** 3-D figure (plotly, `go.Surface` + `go.Scatter3d`): the pair χ = y, ψ = z − x² from
    `ch04.stream_function_pair("parabolic")` — two planes χ = a, b (blue, translucent) and two parabolic sheets ψ = c, d
    (orange), with four streamlines traced by `ch04.streamline` on their intersections, and the patch bounded by the four
    sheets shaded grey with its normal. Preceded by `nb.note` — **N15 [B]** "**Two stream functions.** Since ρu = ∇χ × ∇ψ
    is perpendicular to both gradients, $\rho\mathbf u\cdot\nabla\chi=\rho\mathbf u\cdot\nabla\psi=0$: the flow runs *along* every
    surface χ = const and every surface ψ = const, so a 3-D streamline is where one of each meets (our version of Fig.
    4.1)." equation `\rho\mathbf u\cdot\nabla\chi=\rho\mathbf u\cdot\nabla\psi=0`, no ref. see/read/change after the plot: *see:*
    "planes and parabolic sheets crossing along parabolas; the streamlines lie exactly on the crossings." *read:* "the
    flow u = (1, 0, 2x) climbs as it goes right; every streamline stays in its plane y = const and on its sheet
    z − x² = const." *change:* "…χ = y replaced by χ = y + x: the planes tilt, the streamlines leave the planes y = const,
    and the flux through the patch stays (b − a)(d − c)."
15. `nb.note` — **N16 [B]** "**Flux through a stream tube.** The tube bounded by χ = a, b and ψ = c, d carries
    $\dot m=\oint_C\chi\,d\psi=b(d-c)+a(c-d)=(b-a)(d-c)$ — on the ψ = const legs dψ = 0, on the χ = b leg the integral is
    b(d − c), on χ = a it is a(c − d) (Stokes' theorem, Ch. 2 §2.13)." Followed by `nb.code`:
    `print(ch04.stream_tube_mass_flux("parabolic", "parabolic", 0.0, 1.0, 0.0, 2.0))` *expect:* `(2.0000, 2.0)` (numeric
    patch flux by `curl_flux` vs (b − a)(d − c)).
16. `nb.note` — **N17 [B]** "**Axisymmetric flows** (no swirl, streamlines in planes through the z-axis) use χ = −φ:"
    equation `\rho u_R=-\frac1R\frac{\partial\psi}{\partial z},\qquad\rho u_z=\frac1R\frac{\partial\psi}{\partial R}`, no ref.
    "Test: ψ = ½UR² gives u_z = U, u_R = 0 — a uniform stream along the axis (code:
    `ch04.velocity_from_streamfunction_axisym("axisym_uniform", 0.5, 0.3, U=2.0)` → `(0.0, 2.0)`). Used for the sphere
    (Ch. 6) and pipe flow (Ch. 8); the flux between two stream surfaces is 2π(ψ₂ − ψ₁) here."
17. `nb.note` — **N18 [C]** "With constant density the same construction holds for u itself, and ψ differences are
    volume fluxes (m²/s per metre of depth in 2-D) — every function here takes `rho=1.0` by default. Ch. 6 uses it
    throughout."
18. `nb.explainer("stream_function_spacing", heading="Can one number field hold a whole 2-D flow?", why="Seven classic
    flows drawn as ψ contours at equal steps over a speed map, with tracers moving along them. Drag a gate across the
    flow: its flux stays equal to ψ₂ − ψ₁ however you tilt or bend it; click a point to read
    $u=\partial\psi/\partial y$, $v=-\partial\psi/\partial x$ from the contour slopes.", tries=["Choose 'stagnation' and drag
    the gate between the ψ = 1 and ψ = 2 contours along two different paths — the flux readout stays 1 m²/s.", "Choose
    'cylinder': where do the contours crowd? Read the speed there in Explain.", "Choose 'line vortex': the contours are
    circles packed tighter toward the centre — the speed grows as 1/r.", "Open the Derivation tab and step through D04 with
    the uniform stream."])`
19. `nb.md` — **What would change if…** "…the flow were not steady? In 2-D incompressible flow ψ still exists at every
    instant (∇·u = 0 is enough), but its contours are then the *instantaneous* streamlines, not paths (Ch. 3 §3.3).
    Knowing the flow pattern, the next question is what forces it takes to maintain it — momentum (C04)."

### A.4 §4.4 Conservation of Momentum — R03, R04, C04 (+N19–N29, N31, N32, N82–N84 · D05 · E1), C05 (+N30, N103–N105 · D06), C06 (+N02, N33–N37 · D07)
1. `nb.section("4.4", "Conservation of Momentum", intro="**What is this section about?** Newton's second law for a
   material volume, carried to any control volume (forces from fluxes: wakes, bores, rockets, sprinklers), Bernoulli's
   equation from a thin stream tube, and — shrunk to a point — Cauchy's equation of motion, Newton's law for every
   continuum.")`
2. `nb.recap("R03", "Surface force from the stress tensor", "On a surface element with unit normal n the force per unit
   area is $f_j=n_i\tau_{ij}$ *(Eq. 2.15)*: the **first** index of τ is contracted with the normal. Its normal part is
   $\mathbf n\cdot\mathbf f=n_if_i$ and its tangential part $f_k-(n_if_i)n_k$. For a fluid at rest $\tau_{ij}=-p\delta_{ij}$ and
   f = −pn, a push along the inward normal.", where="Ch. 2 §2.6")` — followed by `nb.code`: `tau = np.array([[-1e5, 3.0,
   0], [3.0, -1e5, 0], [0, 0, -1e5]])`; `n = np.array([1, 1, 0])/np.sqrt(2)`; `print(tensors.traction(tau, n),
   tensors.normal_shear_stress(tau, n))` (Ch. 2's functions). *expect:* f ≈ (−70708.6, −70708.6, 0) Pa, normal
   −99997.0 Pa, shear 0.0 (on the 45° plane the shear pair adds to the normal part). *explain:* "traction and its normal and
   shear parts for pressure plus a small shear stress."
3. `nb.recap("R04", "Surface tension acts through boundary conditions", "Surface tension is a force per unit length
   acting along lines in an interface (Ch. 1 §1.6), not a force on the bulk fluid, so it never appears in the field
   equations of this chapter; it enters through the conditions at a free surface or an interface — derived from a
   force balance on a curved cap in C14.", where="Ch. 1 §1.6")`

**C04 — Momentum conservation for an arbitrarily moving control volume (4.17)**
4. `nb.core("C04", "Forces from fluxes: momentum in a moving box (4.17)", question="A bar is held in a wind tunnel on
   a thin sting. You may not touch the bar or the sting — you may only measure the air speed behind it. Can you tell how
   hard the air pushes on the bar?")`
5. `nb.md` — **The problem in plain words:** "Measuring the force on a body directly is often impossible (a ship's hull,
   a turbine blade in a machine, the Earth's surface under a wind). But the fluid leaving a region carries momentum
   away, and Newton's law says the momentum a box loses per second must be supplied by forces. Count the momentum
   flowing in and out of a box drawn round the body and you have weighed the force. The same bookkeeping gives the speed
   of a tidal bore, the thrust of a rocket and the torque of a lawn sprinkler."
6. `nb.md` — **The idea:**
   ```
   rate of change of momentum inside  +  momentum carried OUT through the walls (relative to them)
        d/dt ∫ ρu dV                   +  ∮ ρu (u − b)·n dA
   =   body forces on the fluid       +  surface forces on the fluid
        ∫ ρg dV                        +  ∮ f dA                                   (4.17)
   steady flow, fixed box:  (momentum out − momentum in) = forces on the fluid  → force on the body = −(that)
   ```
7. `nb.note` — **N26 [B]** "**Body forces vs surface forces.**" Two-column table: body — act without contact, ∝ mass,
   specified per unit mass (an acceleration, m/s²): gravity, electromagnetic, and in accelerating frames the fictitious
   forces of §4.7 · surface — act by contact, ∝ area, specified per unit area (a stress, Pa): pressure and viscous
   stress. "In (4.17) the first is $\int\rho\mathbf g\,dV$, the second $\int\mathbf f\,dA$ with $f_j=n_i\tau_{ij}$ (R03)."
8. `nb.primer("momentum flux through a surface", "Fluid crossing a surface carries its momentum with it. Through a patch
   dA the *volume* crossing per second is (u − b)·n dA (velocity relative to the patch, component along the normal); its
   mass is ρ times that, and its momentum is u times that mass: $\rho\mathbf u\,[(\mathbf u-\mathbf b)\cdot\mathbf n]\,dA$ — a
   vector (the direction of u) times a scalar (the rate of crossing). A jet of water 1 cm² in area at 10 m/s through a
   fixed nozzle carries ρU²A = 1000 × 100 × 1e-4 = 10 N of momentum per second.", code="rho, U, A = 1000.0, 10.0, 1e-4
   # water, speed [m/s], area [m²]\nmdot = rho*U*A                       # mass flux [kg/s]\nprint(mdot, mdot*U)                  # 1.0 kg/s carries 10.0 N of x-momentum per second")`
   — *expect:* `1.0 10.0`.
9. `nb.note` — **N19 [B]** "**Newton's second law for a material volume** — momentum per volume ρu, body force per mass
   g, surface force per area f, outward n (Newton II, Ch. 1 primer P09):" equation
   `\frac{d}{dt}\int_{V(t)}\rho\mathbf u\,dV=\int_{V(t)}\rho\mathbf g\,dV+\int_{A(t)}\mathbf f(\mathbf n,\mathbf x,t)\,dA`, ref "4.13".
10. `nb.derivation("D05", …)` — Part F D05 (8 steps), ref "4.17".
11. `nb.note` — **N20 [B]** "Step 2 of D05 is (4.14), the transport theorem applied to each component of ρu with b = u:"
    equation `\int_{V(t)}\frac{\partial}{\partial t}(\rho\mathbf u)dV+\int_{A(t)}\rho\mathbf u(\mathbf u\cdot\mathbf n)dA=\int_{V(t)}\rho\mathbf g\,dV+\int_{A(t)}\mathbf f\,dA`,
    ref "4.14". "It is also the starting line of Cauchy's equation (C06)."
12. `nb.note` — **N21 [B]** "Step 3 is (4.15), the same theorem for our box, rearranged:" equation
    `\int_{V^*}\frac{\partial}{\partial t}(\rho\mathbf u)dV=\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV-\int_{A^*}\rho\mathbf u\,(\mathbf b\cdot\mathbf n)\,dA`,
    ref "4.15". "> ⚠️ **Book slip:** the printed (4.15) ends with an extra '= 0'. That would claim ∂(ρu)/∂t integrates to
    zero; (4.15) is an identity between three terms (compare (4.3), where '= 0' is right because the mass statement was
    already substituted)."
13. `nb.note` — **N22, N23, N24, N25 [B]** "Step 4 is the four coincidence equalities (4.16a–d): at the instant V* = V, the same
    integrand over the same region gives the same number —" equation
    `\int_{V}\frac{\partial(\rho\mathbf u)}{\partial t}dV=\int_{V^*}\frac{\partial(\rho\mathbf u)}{\partial t}dV,\ \int_{A}\rho\mathbf u(\mathbf u\cdot\mathbf n)dA=\int_{A^*}\rho\mathbf u(\mathbf u\cdot\mathbf n)dA,\ \int_{V}\rho\mathbf g\,dV=\int_{V^*}\rho\mathbf g\,dV,\ \int_{A}\mathbf f\,dA=\int_{A^*}\mathbf f\,dA`,
    ref "4.16a–d" (the builder writes the four as an `aligned` block, one per line). "— exactly the move of D01, now for a
    vector."
14. `nb.primer("conservative force and its potential", "A force field is *conservative* if the work it does on a particle
    moving from A to B does not depend on the route. Then it is minus the gradient of a potential energy per unit mass Φ,
    and the work round any closed loop is zero. Gravity near the ground: Φ = gz, force per mass −∇Φ = (0, 0, −g).",
    code="loop = itg.planar_loop(center=(0.0, 0.0, 0.0), normal=(1.0, 0.0, 0.0), radius=2.0)   # a vertical circle\ng = lambda x: np.tile(np.array([[0.0], [0.0], [-9.81]]), (1, np.shape(x)[1]))    # uniform gravity field\nprint(itg.circulation(g, loop))    # work of g round the loop per unit mass: 0.0")`
    — *expect:* `0.0` (to 1e-12).
15. `nb.note` — **N27 [B]** "**Conservative body force**, with z up:" equation
    `\mathbf g=-\nabla\Phi\quad\text{or}\quad g_j=-\partial\Phi/\partial x_j`, ref "4.18". "For gravity Φ = gz gives
    $\mathbf g=-g\mathbf e_z$; `ch04.body_force_from_potential("gravity", np.array([0.0, 0.0, 5.0]))` → (0, 0, −9.81) m/s².
    C11 uses Φ to fold gravity into the Bernoulli function."
16. `nb.md` — "> ⚠️ **Common confusion — whose force?** Drag F_D is the force *on the body*, positive downstream. The
    forces in (4.17) are forces *on the fluid in the box*. By Newton's third law the body pushes the fluid with −F_D e_x
    (N28). Forget the minus and the drag comes out negative."
17. `nb.worked_example("drag of a bar from its wake (Ex. 4.1, our synthetic wake)", "Far downstream the speed is
    $U(y)=U_\infty-\Delta\,e^{-y^2/b^2}$ with U∞ = 10 m/s, Δ = 2 m/s, b = 0.1 m; air ρ = 1.2 kg/m³; the box is H tall, its
    inlet far upstream (uniform U∞) and its outlet in the wake; pressure is p∞ on all faces (so the pressure forces cancel —
    Gauss on a constant). 1. **Mass:** less leaves through the outlet than enters the inlet, by
    $\int(U_\infty-U)dy=\Delta b\sqrt\pi=2\times0.1\times1.7725=0.354$ m²/s per metre of span — that must leave through the top
    and bottom (the Gaussian integral $\int e^{-y^2/b^2}dy=b\sqrt\pi$, gloss). 2. That side outflow carries x-momentum
    U∞ per unit mass: ρU∞ × 0.354. 3. **Momentum x:** out − in = $\rho\int U^2dy+\rho U_\infty(0.354)-\rho\int U_\infty^2dy$, which
    must equal the force on the fluid, −F_D/l. 4. Collect: $F_D/l=\rho\int U(U_\infty-U)dy$ (N29). 5. Numbers: $\int U(U_\infty-U)dy
    =U_\infty\Delta b\sqrt\pi-\Delta^2b\sqrt{\pi/2}=3.5449-0.5013=3.0436$ m³/s²; × 1.2 = **3.65 N/m** of span.")`
18. `nb.note` — **N29 [B]** "**Ex. 4.1 result** (given, the derivation in the worked example above):" equation
    `F_D/l=\rho\int_{-H/2}^{+H/2}U(y)\big(U_\infty-U(y)\big)dy`, no ref. "It needs the top and bottom faces far enough out
    that shear and pressure deviations there vanish; the integral of U(1 − U/U∞)/U∞ is the *momentum thickness* of Ch. 9."
19. `nb.code` — *code:* `y = np.linspace(-1, 1, 5)` · `print(ch04.gaussian_wake(y, 10.0, 2.0, 0.1))` ·
    `print(ch04.wake_side_outflow("gaussian", 10.0, 2.0))` (U∞ = 10, H = 2 m; keywords of `gaussian_wake` defaulted) ·
    `print(ch04.wake_drag_per_span("gaussian", 10.0, 1.2, 2.0, return_error=True))` · `sc = ch04.cv_scenario("wake",
    U_inf=10.0, deficit=2.0, width=0.1, rho=1.2, H=2.0)`; `for f in sc["faces"]: print(f"{f['name']:6s} mass
    {f['mass_flux']:+.4f}  x-mom {f['momentum_flux_x']:+.4f}")`; `print(sc["result"], sc["residual_mass"],
    sc["residual_momentum"])`. *expect:* the wake profile (10, 10, 8, 10, 10) at y = −1, −0.5, 0, 0.5, 1 · `0.35449` ·
    `(3.6523, <1e-10)` · faces (per metre of span): inlet mass −24.0 kg/(m s) (in), outlet +23.5746, top/bottom +0.2127 each; inlet
    x-mom −240.0 N/m, outlet +232.0938, top/bottom +2.1269 each (sum −3.6523 = force on the fluid) · `3.6523 0.0 0.0` (≤ 1e-9). *explain:* 1. our synthetic wake; 2. the
    mass leaking through the sides; 3. the drag formula by adaptive quadrature with its error estimate; 4. the same
    budget face by face (signs: out positive, in negative) — the momentum surplus of the fluid is the force on it,
    −F_D/l.
20. `nb.check_agree` — **from scratch (curation §7):** the four face integrals by `np.trapezoid`: `y = np.linspace(-1, 1,
    4001); U = ch04.gaussian_wake(y, 10, 2, 0.1)`; `m_in = -1.2*np.trapezoid(10*np.ones_like(y), y)`; `m_out =
    1.2*np.trapezoid(U, y)`; `m_side = -(m_in + m_out)` (continuity); `p_in = -1.2*np.trapezoid(100*np.ones_like(y), y)`;
    `p_out = 1.2*np.trapezoid(U**2, y)`; `p_side = m_side*10.0`; `FD = -(p_in + p_out + p_side)`; `assert
    np.isclose(FD, ch04.wake_drag_per_span("gaussian", 10.0, 1.2, 2.0), rtol=1e-8)`.
21. `nb.figure` — two panels: (a) the control volume (our Fig. 4.2 analogue, `scripts.ch04_drawings.cv_box`) with the
    bar, the U(y) profiles at inlet (uniform) and outlet (Gaussian deficit), blue inflow and orange outflow arrows on each
    face sized by mass flux, and a rose arrow "−F_D/l on the fluid"; (b) budget bars for x-momentum: inlet (blue,
    negative), outlet (orange), top + bottom (orange), sum = −F_D/l (rose). *see:* "a thin deficit in the outlet profile
    and a small but nonzero leak through the sides." *read:* "the fluid leaves with less momentum than it brought; the
    difference, including what leaks out sideways at speed U∞, is the drag on the bar." *change:* "…the box made taller:
    the side leak and every face flux change, the drag does not (it depends only on the deficit)."
22. `nb.plotly` — `slider_figure` over the box height H = 0.1…2 m (20 steps, FAST 10): x = y [m], traces "U(y) [m/s]"
    and "integrand ρU(U∞ − U) [N/m²]" with the faces ±H/2 as vertical lines; the title carries "F_D/l = … N/m (limit
    3.652)". Title "A box must enclose the whole wake". see/read/change: "…below H ≈ 3b the box cuts the wake and misses
    part of the deficit; above it F_D/l stops changing."
23. `nb.note` — **N31 [B]** "**Ex. 4.3 — a small bore (surge) with a moving control volume.** Dimensional analysis first:
    with g, h and the speed U there is one group, U²/(gh), so U ∝ √(gh). Riding with the wave (b = U e_x) the flow in the
    box is steady; hydrostatic pressure on the vertical faces, $p_o$ on the free surface (its net push cancels), mass and
    momentum give" equation `U^2\frac{h_{in}}{h_{out}}=\frac g2(h_{in}+h_{out}),\qquad U=\sqrt{\frac{gh_{out}}{2h_{in}}(h_{in}+h_{out})}\approx\sqrt{gh}`,
    no ref. "Number: still water h_in = 1 m ahead of the wave, h_out = 1.1 m behind it → U = √(9.81 × 1.1 × 2.1/2) = 3.37
    m/s, close to √(gh) = 3.13 m/s (the shallow-water wave speed of Ch. 7 and Ch. 13)." Followed by `nb.code`:
    `print(ch04.bore_speed(1.0, 1.1), np.sqrt(9.81*1.0))` · `print(ch04.bore_pressure_force(1.0, 1.1))`. *expect:*
    `3.3661 3.1321` · left/right hydrostatic forces 4905.0 and 5935.1 N per metre width plus the p_o parts; `net`
    independent of p_o (the test calls it with two p_o values).
24. `nb.animation` — (`player="video"`, 60 frames, FAST 30) a bore moving from right to left into still water (depth 1 m,
    step to 1.1 m, surface drawn with a smooth tanh front), the control volume (dashed box) riding with it, water parcels
    as dots (in the lab frame they are pushed left and stop behind the front); inset text "in the box's frame the flow is
    steady". see/read/change: *see:* "the box slides with the front; parcels ahead are at rest, parcels behind move
    slowly left." *read:* "seen from the box, water enters fast and deep 1 m, leaves slow and deep 1.1 m — a steady budget
    whose solution is the front's speed." *change:* "…h_out → h_in: the front becomes a small wave and its speed tends to
    √(gh)."
25. `nb.note` — **N32 [B]** "**Ex. 4.4 — a rocket as an accelerating control volume** (b = b(t) e_z). Exhaust leaves the
    nozzle at V_e relative to the rocket; with F_S the sum of drag and pressure thrust, mass and vertical momentum give"
    equation `\frac{dM}{dt}+\rho_eV_eA_e=0,\qquad M\frac{d^2z_R}{dt^2}=-V_e\frac{dM}{dt}-Mg+F_S`, no ref. "With g = F_S = 0
    it integrates to Tsiolkovsky's $\Delta b=V_e\ln(M_0/M_1)$ (our extension, gloss): M₀ = 2 kg, M₁ = 1 kg, V_e = 800 m/s
    → 554.5 m/s." Followed by `nb.figure`: `r = ch04.rocket_trajectory(2.0, 0.1, 800.0, 10.0, t_eval=np.linspace(0, 10,
    201))`; speed b(t) (purple) with the Tsiolkovsky ghost V_e ln(M₀/M(t)) − gt (grey dashed). *see:* "the speed grows
    faster and faster as the rocket gets lighter." *read:* "the thrust V_e|dM/dt| = 80 N is constant but it pushes less
    and less mass." *change:* "…V_e doubled: every speed doubles (minus the gravity loss gt, which does not change)."
26. `nb.primer("torque, moment arm and moment of inertia", "A force F applied at position r (from an axis point) turns a
    body with torque M = r × F; only the part of F perpendicular to r, times the distance, counts (the moment arm). The
    angular analogue of mass is the moment of inertia: for a cube of side h and density ρ spinning about an axis through
    its centre, I = ρh⁵/6. Rate of change of angular momentum = net torque.", code="r = np.array([0.2, 0.0, 0.0])
    # lever 0.2 m along x\nF = np.array([0.0, 5.0, 0.0])                  # 5 N pushing along y\nprint(np.cross(r, F))                          # torque (0, 0, 1.0) N m about z\nrho, h = 1000.0, 0.01; print(rho*h**5/6)       # I of a 1 cm water cube: 1.67e-08 kg m²")`
    — *expect:* `[0. 0. 1.]` · `1.6666666666666667e-08`. (`np.cross` gloss: the cross product of Ch. 2 on arrays.)
27. `nb.note` — **N82 [C]** "For a rigid body, the rate of change of angular momentum equals the applied torque:"
    equation `d\mathbf H/dt=\mathbf M`, ref "4.64". "Its fluid version for a fixed control volume follows."
28. `nb.note` — **N83 [B]** "**Angular momentum for a stationary control volume** (the transport theorem with
    F = r × ρu, the same move as D05):" equation
    `\frac{d}{dt}\int_{V_o}(\mathbf r\times\rho\mathbf u)dV+\int_{A_o}(\mathbf r\times\rho\mathbf u)(\mathbf u\cdot\mathbf n)dA=\int_{V_o}(\mathbf r\times\rho\mathbf g)dV+\int_{A_o}(\mathbf r\times\mathbf f)dA`,
    ref "4.65". "Internal torques cancel in pairs — because the stress is symmetric, which C07 proves."
29. `nb.note` — **N84 [B]** "**Ex. 4.6 — a lawn sprinkler held still.** Two jets of area A at radius a leave at speed U,
    tilted by α from the tangential direction; only the tangential part U cos α has a moment arm, so the torque needed to
    hold the arm is" equation `M=2a\rho AU^2\cos\alpha`, no ref. "Number: a = 0.2 m, water, A = 1 cm², U = 5 m/s, α = 30° →
    M = 2 × 0.2 × 1000 × 10⁻⁴ × 25 × 0.866 = 0.866 N m. Let go of it and (with no friction) it spins until the jets leave
    with no tangential speed, at Ω = U cos α/a = 21.7 rad/s (our extension)." Followed by `nb.code`:
    `print(ch04.sprinkler_torque(0.2, 1000.0, 1e-4, 5.0, np.pi/6), ch04.sprinkler_free_spin_rate(0.2, 5.0, np.pi/6))`.
    *expect:* `0.8660 21.65`.
30. `nb.explainer("control_volume_budgets", heading="Weigh a force by counting what flows through a box", why="One
    scene, one control volume, every face's flux as a bar: storage + outflow − inflow on the left of
    $\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV+\int_{A^*}\rho\mathbf u(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA$ (4.17),
    forces on the right, and a residual that stays at zero. Four scenarios (wake, bore, jet, rocket) and three box types
    (fixed, riding, material).", tries=["In 'wake', shrink the box height below three wake widths — the status turns amber
    and the drag readout drops: the box cuts the wake.", "In 'bore', drag the box speed b until the storage bar
    vanishes: that speed is the bore's speed (Explain shows the formula).", "Switch the budget to mass and choose the
    material box (b = u): every flux bar vanishes — that is (4.1).", "Open the Derivation tab and step through D01: the
    coincidence step freezes the balloon on top of the box."])`
31. `nb.md` — **What would change if…** "…the box were a thin stream tube in a steady, frictionless flow? Its side walls
    carry no flux, only pressure; the momentum budget along it becomes a relation between speed, pressure and height —
    Bernoulli's equation (C05)."

**C05 — Bernoulli along a streamline (4.19)**
32. `nb.core("C05", "Bernoulli along a streamline: fast flow, low pressure (4.19)", question="An airliner measures its
    airspeed with a small tube pointing into the wind and a hole in its side — no moving parts. How does a pressure
    difference become a speed?")`
33. `nb.md` — **The problem in plain words:** "Pitot tubes on aircraft and in rivers, the speed of water leaving a hole in
    a tank, the suction on the top of a wing: all trade pressure for speed along the path of the fluid. The trade follows
    from momentum conservation on a very thin stream tube in a steady, frictionless, constant-density flow."
34. `nb.md` — **The idea:**
   ```
   along ONE streamline, per unit mass:   ½U²  +  g z  +  p/ρ  = constant           (4.19)
                                          kinetic  potential  "pressure energy"
   speed up (U ↑) → pressure down (p ↓) at the same height;  stop the flow (U = 0) → pressure up by ½ρU²
   ```
35. `nb.primer("sympy expand, series, removeO, collect and subs", "Four sympy moves for approximations: `expand`
    multiplies brackets out, `series(expr, ds, 0, 2)` keeps powers of ds up to ds¹ and adds an O(ds²) marker that
    `removeO()` drops, `collect(expr, ds)` groups terms by powers of ds, and `subs({a: b})` substitutes. Here: keep only
    first-order terms of a small element of length ds.", code="ds, U, dU = sp.symbols('ds U dU')\ne = sp.expand((U + dU*ds)**2)                       # U**2 + 2*U*dU*ds + dU**2*ds**2\nprint(e.series(ds, 0, 2).removeO())                 # U**2 + 2*U*dU*ds (ds² dropped)\nprint(sp.collect(e, ds), e.subs({U: 2, dU: 1, ds: 0.1}))   # grouped; 4.41")`
    — *expect:* `2*U*dU*ds + U**2` (sympy's order) · the collected form and `4.41000000000000`.
36. `nb.note` — **N30 [B]** "**Ex. 4.2** builds (4.19) from a stream-tube element of length ds: mass with first-order
    changes of U and A, the streamwise momentum with gravity (sin θ ds = dz) and the extra pressure force on the slowly
    widening side, then drops the (ds)² terms to get" equation
    `U\frac{\partial U}{\partial s}ds=-g\,dz-\frac1\rho\frac{\partial p}{\partial s}ds`, no ref. "The derivation below fills in
    the moves the book skips (where the side pressure force comes from, why the mean area multiplies gravity)."
37. `nb.derivation("D06", …)` — Part F D06 (11 steps), ref "4.19", with `check_src` (optional, ★★): the sympy re-run
    `ch04.stream_tube_element_balance_sym()` printing `result` = `U*dU + g*dz + dp/rho` (= 0).
38. `nb.md` — "> ⚠️ **Book slip:** the statement of Ex. 4.2 writes the result as '½ρU² + gz + p/ρ = constant'; the ρ in
    the first term is a misprint (J/m³ added to J/kg). The derived (4.19), $\tfrac12U^2+gz+p/\rho=$ constant along a
    streamline, is per unit mass; multiplied by ρ it reads $p+\tfrac12\rho U^2+\rho gz=$ constant (pressures)."
39. `nb.note` — **N104 [B]** "**Stagnation and dynamic pressure.** Where a streamline meets a body head-on the fluid stops
    (a *stagnation point*, gloss). Between a far point (p, U) and the stagnation point (p₀, 0) at the same height, (4.19)
    gives" equation `p_0=p+\tfrac12\rho\lvert\mathbf u\rvert^2`, no ref. "p₀ is the *stagnation* (total) pressure and ½ρU² the
    *dynamic* pressure."
40. `nb.note` — **N103 [B]** "**Pitot tube** (our Fig. 4.15 analogue in the figure below): the tube's mouth is a stagnation
    point, the side hole reads the static pressure; a manometer between them gives" equation
    `\lvert\mathbf u\rvert_1=\sqrt{2(p_2-p_1)/\rho}=\sqrt{2g(h_2-h_1)}`, no ref. "Number: Δp = 500 Pa in air (ρ = 1.2 kg/m³)
    → 28.9 m/s (104 km/h). Only two points on one streamline are used — irrotationality is not needed (the book's
    wording says 'irrotational', a stronger assumption than required)."
41. `nb.note` — **N105 [B]** "**Orifice (Torricelli).** Water leaving a hole a depth h below the free surface of a large
    tank: at the surface U ≈ 0, p = p_atm; in the jet p = p_atm (straight, parallel streamlines have no pressure
    difference across them). So" equation `u=\sqrt{2gh},\qquad\dot m=\rho A_c\sqrt{2gh}`, no ref. "A_c is the area of
    the jet at its narrowest (the *vena contracta*, gloss) — about 0.61 of a sharp-edged hole's area (a measured value),
    close to 1 for a rounded one. Number: h = 1 m → 4.43 m/s."
42. `nb.worked_example("pitot tube and tank", "1. Pitot: Δp = 500 Pa, ρ = 1.2 kg/m³: $U=\sqrt{2\times500/1.2}=\sqrt{833.3}=28.87$
    m/s; dynamic pressure ½ × 1.2 × 28.87² = 500 Pa ✓. 2. Tank: h = 1 m, g = 9.81 m/s²: $u=\sqrt{2\times9.81\times1}=4.43$
    m/s; a 1 cm² sharp hole (C_c = 0.61) passes ρA_cu = 1000 × 0.61 × 10⁻⁴ × 4.43 = 0.270 kg/s. 3. Both from
    $\tfrac12U^2+gz+p/\rho=$ constant (4.19) between two points of one streamline.")`
43. `nb.code` — *code:* `print(ch04.bernoulli_head(28.87, 0.0, 1.0e5, 1.2), ch04.bernoulli_head(0.0, 0.0, 1.0e5 + 500, 1.2))`
    · `print(ch04.pitot_speed(1.0e5 + 500.0, 1.0e5, 1.2), ch04.dynamic_pressure(28.8675, 1.2),
    ch04.stagnation_pressure(1.0e5, 28.8675, 1.2))` · `print(ch04.torricelli_speed(1.0), ch04.orifice_mass_flow(1.0,
    1e-4, Cc=0.61))` · `print(ch04.bernoulli_solve({"U": 0.0, "z": 1.0, "p": 1.0e5}, {"z": 0.0, "p": 1.0e5}, "U"))` ·
    `print(ch04.stream_tube_element_balance_sym()["result"])`. *expect:* two equal heads 83750.1 and 83750.0 m²/s²
    (≈ 8.375e4) · `28.8675 500.0 100500.0` · `4.4294 0.2702` · `4.4294` · `U*dU + dz*g + dp/rho` (sympy's ordering).
    *explain:* 1. the Bernoulli constant at the free stream and at the stagnation point agree; 2. the pitot inversion and
    the two pressures; 3. Torricelli and the jet's mass flow; 4. the general solver (solve (4.19) for the one unknown);
    5. the sympy re-run of Ex. 4.2.
44. `nb.check_agree` — **from scratch:** drain a tank by explicit Euler steps of $dh/dt=-(C_cA_o/A_t)\sqrt{2gh}$
    (A_t = 1 m², A_o = 10⁻³ m², h₀ = 1 m, dt = 0.5 s) and compare the emptying time with the closed form
    $t_{empty}=(A_t/C_cA_o)\sqrt{2h_0/g}$ from `ch04.tank_drain(1.0, 1.0, 1e-3)["t_empty"]` = 451.5 s: `assert abs(t_euler -
    451.5) < 2.0`.
45. `nb.figure` — two panels: (a) pressure coefficient on the surface of the Ch. 3 cylinder from (4.19),
    $C_p=1-4\sin^2\theta$ (teal) vs θ, with the stagnation points (C_p = 1) and the shoulders (C_p = −3) marked, computed by
    `ch04.pressure_coefficient` on `ch03.cylinder_flow` velocities; (b) tank level h(t) for C_c = 1 (dashed) and 0.611
    (solid) from `ch04.tank_drain`. *see:* "(a) the pressure is highest where the flow stops and lowest where it is
    fastest; (b) a parabola in t that reaches zero sooner without the vena contracta." *read:* "(a) C_p = 1 − (U/U∞)² is
    (4.19) in dimensionless form; (b) the outflow slows as √h, so the level falls fast at first and slowly at the end."
    *change:* "…viscosity added: behind a real cylinder the flow separates and the pressure never recovers to C_p = 1
    (Ch. 9) — Bernoulli fails where friction acts."
46. `nb.md` — **What would change if…** "…we asked for momentum at every *point* rather than along one streamline in
    frictionless flow? We would need the stresses on a tiny cube — pressure and friction together — and Newton's law for
    it: Cauchy's equation (C06). Bernoulli's *general* forms return in C11 and C12, where E6 compares all four."

**C06 — Cauchy's equation of motion (4.24)**
47. `nb.core("C06", "Newton's law for a fluid particle: Cauchy's equation (4.24)", question="A small cube of fluid is
    pushed by its neighbours on six faces and pulled by gravity. Before we know anything about what kind of fluid it is,
    what is F = ma for it?")`
48. `nb.md` — **The problem in plain words:** "The integral law (4.17) needs a box; to predict a flow everywhere we need the
    law at every point — a differential equation for u. We can get it without any knowledge of how stress depends on
    motion: that part comes later (C07). The result holds for water, air, honey, even a steel beam."
49. `nb.md` — **The idea:**
   ```
   (4.14) for a material volume  →  Gauss on BOTH surface integrals  →  one volume integral = 0 for every volume
                                 →  localise  →  flux form (4.22)  →  subtract u_j × continuity  →  ρ Du_j/Dt = ρ g_j + ∂τ_ij/∂x_i
   mass × acceleration          =   gravity   +   net push of the neighbours (divergence of stress, FIRST index)
   ```
50. `nb.primer("tensor divergence over the first index", "The net surface force per unit volume on a small cube is the
    divergence of the stress. Because the traction on a face with normal n is $f_j=n_i\tau_{ij}$ (first index = face),
    Gauss' theorem turns $\oint n_i\tau_{ij}dA$ into $\int\partial\tau_{ij}/\partial x_i\,dV$ — the derivative acts on the
    **first** index. For a symmetric τ it does not matter; for a non-symmetric one it does, and C07 has not yet proved τ
    symmetric. (Ch. 2's `tensor_divergence` differentiated the second index by default; here we pass `index=0`.)",
    code="tau = lambda x: np.array([[0.0, x[0], 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])   # τ12 = x1, not symmetric\nfirst, second = ch04.divergence_first_index_demo(tau, np.array([0.3, 0.2, 0.1]))\nprint(first, second)     # [0. 1. 0.] vs [0. 0. 0.]: different vectors")`
    — *expect:* `[0. 1. 0.] [0. 0. 0.]`.
51. `nb.derivation("D07", …)` — Part F D07 (10 steps), ref "4.24".
52. `nb.note` — **N33 [B]** "D07's step 2 is (4.20a), Gauss applied to each component j of the momentum flux, the
    normal's index contracted:" equation
    `\int_{A}\rho\mathbf u(\mathbf u\cdot\mathbf n)dA=\int_{V}\nabla\cdot(\rho\mathbf u\mathbf u)dV=\int_{V}\frac{\partial}{\partial x_i}(\rho u_iu_j)dV`, ref "4.20a".
53. `nb.note` — **N34 [B]** "Step 3 is (4.20b), the same for the surface force — the first index again:" equation
    `\int_{A}\mathbf f\,dA=\int_{A}n_i\tau_{ij}dA=\int_{V}\frac{\partial\tau_{ij}}{\partial x_i}dV`, ref "4.20b".
54. `nb.note` — **N35 [B]** "Step 4 is (4.21), everything under one integral:" equation
    `\int_{V}\Big\{\frac{\partial}{\partial t}(\rho u_j)+\frac{\partial}{\partial x_i}(\rho u_iu_j)-\rho g_j-\frac{\partial\tau_{ij}}{\partial x_i}\Big\}dV=0`, ref "4.21".
55. `nb.note` — **N36 [B]** "Step 5 (localisation) gives the **conservative (flux) form**, the form finite-volume codes
    discretise (Ch. 10):" equation
    `\frac{\partial}{\partial t}(\rho u_j)+\frac{\partial}{\partial x_i}(\rho u_iu_j)=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}`, ref "4.22".
56. `nb.note` — **N37 [B]** "Steps 6–8 expand the left side; the bracket is the continuity equation (4.7) and vanishes:"
    equation `\frac{\partial}{\partial t}(\rho u_j)+\frac{\partial}{\partial x_i}(\rho u_iu_j)=\rho\frac{\partial u_j}{\partial t}+u_j\Big[\frac{\partial\rho}{\partial t}+\frac{\partial}{\partial x_i}(\rho u_i)\Big]+\rho u_i\frac{\partial u_j}{\partial x_i}=\rho\frac{Du_j}{Dt}`,
    ref "4.23". Followed by `nb.code`: `x, y, z, t = sp.symbols('x y z t')`; `r = sp.Function('rho')(x, y, z, t)`; `U =
    [sp.Function(f'u{i}')(x, y, z, t) for i in (1, 2, 3)]`; `print([sp.simplify(e) for e in
    ch04.conservative_to_advective_sym(r, U, (x, y, z), t)])`. *expect:* each component printed as `u_j*(∂ρ/∂t +
    ∂(ρu_i)/∂x_i)` — "the difference is u_j times continuity, so it vanishes whenever (4.7) holds".
57. `nb.md` — "> ⚠️ **Book slip (index):** the sentence after (4.24) calls the net surface force ∂τ_ij/∂x_j. The equation
    itself contracts the first index, $\partial\tau_{ij}/\partial x_i$. The two agree only because τ is symmetric — proved in
    C07, *after* (4.24)."
58. `nb.note` — **N02 [B] — the ledger, first row.** "Continuity (1 equation), Cauchy (3) and two thermodynamic equations
    of state give 6 equations; the unknowns are ρ (1), u_j (3) and τ_ij (9): 13. Not solvable yet — a *constitutive law*
    must tie τ to the motion (C07)." Followed by `nb.code`: `ledger.loc["cauchy"] = ch04.closure_count("cauchy")` (the
    §4.1 table) · `print(ledger)`. *expect:* the cauchy row `6 13`, the others still 0.
59. `nb.worked_example("water spinning like a solid body (Ω = 1 rad/s) at r = 0.1 m", "Velocity $u_\varphi=\Omega r=0.1$
    m/s; every particle moves on a circle, so its acceleration is centripetal: $a=-\Omega^2r\,\mathbf e_R=-0.1$ m/s²,
    ρa = −100 N/m³. The pressure that holds it on the circle is $p=\rho\Omega^2R^2/2-\rho gz$ (plus a constant):
    $-\partial p/\partial R=-\rho\Omega^2R=-100$ N/m³ radially, and $-\partial p/\partial z=+\rho g$ vertically balancing gravity
    $\rho g_z=-\rho g$. No viscous stress (a rigid motion has S = 0, Ch. 3). So (4.24),
    $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ with $\tau_{ij}=-p\delta_{ij}$, reads −100 = 0 + (−100)
    radially and 0 = −9810 + 9810 vertically ✓.")`
60. `nb.code` — *code:* `u_fn, p_fn = ch04.exact_field("solid_body", Omega=1.0, rho=1000.0)` · `tau = lambda x, t:
    -p_fn(x, t)*np.eye(3)` · `x0 = np.array([0.1, 0.0, 0.0])` · `T = ch04.cauchy_terms(1000.0, u_fn, tau, np.array([0,
    0, -9.81]), x0, 0.0, h=1e-4)` · `print(T.inertia, T.body, T.stress_divergence, T.residual)` · `print(
    ch04.momentum_conservative_residual(1000.0, u_fn, tau, np.array([0, 0, -9.81]), x0, 0.0))`. *expect:* inertia ≈
    (−100, 0, 0), body (0, 0, −9810), stress divergence ≈ (−100, 0, 9810), residual ≈ 0 (≤ 1e-5) N/m³; the conservative
    residual ≈ 0 as well. *explain:* 1. the exact solid-body field and its pressure; 2. the stress of a fluid with no
    viscous part; 3. the three terms of (4.24) at a point — mass × acceleration equals gravity plus the net push; 4. the
    flux form (4.22) gives the same balance (continuity holds).
61. `nb.check_agree` — **from scratch (curation §7):** first-index stress divergence by central differences in a loop:
    `div = np.zeros(3)`; `for i in range(3): e = np.eye(3)[i]*1e-4; div += (tau(x0 + e, 0)[i, :] - tau(x0 - e, 0)[i,
    :])/2e-4` (row i differentiated along x_i) · `assert np.allclose(div, T.stress_divergence, rtol=1e-6)` · the same loop
    with a non-symmetric test τ (the primer's) shows the second-index version failing: `assert not np.allclose(first,
    second)`.
62. `nb.figure` — term bars at four points of the rotating water (R = 0.05, 0.1, 0.15, 0.2 m) for the radial component:
    ρDu/Dt (teal), ρg (grey, 0 radially), ∂τ_ij/∂x_i (orange), residual (black dot); a second small panel for the vertical
    component at one point. *see:* "teal and orange bars of equal length, growing with R; nothing grey radially." *read:*
    "the inward push of the pressure (orange) supplies exactly the centripetal acceleration (teal); vertically the
    pressure holds the water up against gravity." *change:* "…Ω doubled: both radial bars grow four times (Ω²); the
    vertical balance is unchanged."
63. `nb.md` — **What would change if…** "…we knew how τ depends on the motion? Then the 13 unknowns collapse to 5 (ρ, p,
    u_j) and Cauchy's equation becomes a closed PDE for the velocity. The simplest possible law — stress linear in the
    strain rate, the same in every direction — is the Newtonian fluid (C07)."

### A.5 §4.5 Constitutive Equation for a Newtonian Fluid — R05, R06, C07 (+N38–N51 · D08, D09, D10 · E3)
1. `nb.section("4.5", "Constitutive Equation for a Newtonian Fluid", intro="**What is this section about?** Cauchy's
   equation needs a rule for the stress. We show first that the stress tensor is symmetric, then build the simplest law
   allowed by physics — stress linear in the rate of strain and the same in every direction — and find it has only two
   material constants, μ and the bulk viscosity μ_v.")`
2. `nb.recap("R05", "Only the strain rate can create viscous stress", "Viscous stress must be the same for every
   observer moving at constant velocity (Galilean invariance, Ch. 3 §3.3), so it cannot depend on u itself — only on
   the velocity gradient G. And of G = S + ½R *(Eq. 3.11)* only the strain-rate tensor S deforms an element; the rotation
   part R describes a rigid spin, which deforms nothing (Ch. 3 §3.4: a rigid motion U + Ω × x has S = 0). So the viscous
   stress is a function of S alone and vanishes when S = 0.", where="Ch. 3 §3.3–3.4")` — followed by `nb.code`:
   `G_rigid = kinematics.velocity_gradient_preset("solid_body_rotation", Gamma=2.0, dim=3)`
   · `print(tensors.strain_rate_tensor(G_rigid), ch04.viscous_stress(G_rigid, mu=1e-3))`. *expect:* two 3 × 3 zero
   matrices. *explain:* "a rigid rotation has no strain rate and so no viscous stress."
3. `nb.recap("R06", "Non-Newtonian fluids exist", "Paint, ketchup, blood, polymer melts and mud do not obey a linear
   stress–strain-rate law: they shear-thin, have a yield stress or remember their past (Ch. 1 §1.3 named Bingham and
   Maxwell materials). This chapter treats the Newtonian fluid only; Ch. 16 returns to the others.", where="Ch. 1 §1.3")`

**C07 — The Newtonian stress law (4.31), (4.37)**
4. `nb.core("C07", "How a fluid decides its stress: the Newtonian law (4.31)", question="Honey between two plates: slide
   the top plate and it resists; squeeze the plates together and it resists differently; spin the whole sandwich and it
   does not resist at all. What rule turns 'how the fluid deforms' into 'what stress it carries'?")`
5. `nb.md` — **The problem in plain words:** "Newton's law τ = μ du/dy (1.3) covers one flow: parallel layers sliding. A
   real flow shears, stretches and swells in all directions at once, and Cauchy's equation needs all nine stress
   components. We want the general rule, as simple as physics allows, with as few material constants as possible —
   the one used for water and air in every later chapter, in CFD codes and in the eddy-viscosity models of turbulence."
6. `nb.md` — **The idea:**
   ```
   motion  G = ∂u_i/∂x_j  →  keep only S (R05)  →  linear: σ_ij = K_ijmn S_mn (81 numbers)
                          →  isotropic: K = λδδ + μδδ + γδδ (3 numbers)  →  symmetric: 2 numbers
   τ_ij = −p δ_ij  +  2μ S_ij  +  λ S_mm δ_ij                                      (4.31)
          static     shear/stretch  swelling
   ```
7. `nb.note` — **N39 [C]** "A **constitutive equation** is a material's rule linking stress to deformation; the
   *Newtonian* fluid has the simplest linear one. Non-Newtonian rules: Ch. 16."
8. `nb.note` — **N38 [B]** "**The stress tensor is symmetric**, so it has six independent components:" equation
   `\tau_{ij}=\tau_{ji}`, ref "4.25". "The book leaves the proof to an exercise; the derivation below writes it out: an
   unequal pair τ₁₂ ≠ τ₂₁ would spin a small cube faster and faster as it shrinks. The only exception is a fluid with
   body *couples* (torques per unit mass, e.g. polarised molecules in an electric field)."
9. `nb.derivation("D08", …)` — Part F D08 (8 steps), ref "4.25" (the torque primer of C04 is recalled in one line
   before it: "📎 torque and moment of inertia — primed in C04 above").
10. `nb.code` — *code:* `h = np.logspace(-4, -1, 7)` · `alpha = ch04.cube_spin_acceleration(1.0, 0.0, 1000.0, h)` ·
    `print(np.c_[h, alpha])` · `print(observed_order(h, alpha))`. *expect:* α = 60 rad/s² at h = 1 cm, 6000 at 1 mm, 6×10⁵
    at 0.1 mm · slope −2.0. *explain:* 1. the angular acceleration a 1 Pa stress imbalance would give a water cube of side
    h; 2. its log–log slope: it grows as 1/h² — no finite imbalance can survive at a point.
11. `nb.note` — **N40 [B]** "**At rest** a fluid's stress is the same on every plane (isotropic, Ch. 1 §1.6), and the only
    isotropic second-order tensor is δ_ij (Ch. 2 §2.5):" equation `\tau_{ij}=-p\,\delta_{ij}`, ref "4.26". "p is the
    thermodynamic pressure (e.g. p = ρRT); the minus sign because tension counts positive. Traction on any plane: −pn
    (`tensors.traction(ch04.static_stress(1e5), n)` = −10⁵ n for any unit n)."
12. `nb.note` — **N41 [B]** "**Moving fluid:** add a viscous part σ_ij that vanishes at rest:" equation
    `\tau_{ij}=-p\,\delta_{ij}+\sigma_{ij}`, ref "4.27". "p stays the thermodynamic pressure because fluid particles are
    in local equilibrium (Ch. 1 §1.8). > ⚠️ The book calls σ the 'deviatoric' stress; it is traceless (a true deviator)
    only when μ_v = 0 or ∇·u = 0, since its trace is 3μ_v∇·u (shown after D10)."
13. `nb.note` — **N42 [B]** "**Most general linear law** with σ = 0 when S = 0:" equation `\sigma_{ij}=K_{ijmn}S_{mn}`, ref
    "4.28". "K has 3⁴ = 81 components (each of the nine σ's may depend on each of the nine S's); code:
    `ch04.linear_stress(K, S)` = `np.einsum('ijmn,mn->ij', K, S)` (Ch. 2 primer P62)."
14. `nb.primer("isotropic fourth-order tensor", "A tensor is *isotropic* if its components are the same in every rotated
    frame. For second order only multiples of δ_ij qualify. For fourth order the only possibilities are combinations of
    products of two δ's — three of them: $\delta_{ij}\delta_{mn}$, $\delta_{im}\delta_{jn}$, $\delta_{in}\delta_{jm}$ (a
    classical result we cite, not prove). In a fluid with no preferred direction K must be isotropic, which cuts 81
    constants to 3. We check it numerically: rotate K with 50 random rotations and nothing changes.", code="K =
    ch04.isotropic_fourth_order(2.0, 1.0, 0.5)                   # λ, μ, γ\nC = tensors.random_rotation(np.random.default_rng(1))       # a random rotation matrix\nprint(np.abs(tensors.transform_tensor(K, C) - K).max())     # ~1e-15: unchanged\nprint(tensors.is_isotropic(K, n_rotations=50))              # True")`
    — *expect:* a number ≤ 1e-14 and `True`.
15. `nb.note` — **N43 [B]** "**Isotropy** leaves three scalars:" equation
    `K_{ijmn}=\lambda\delta_{ij}\delta_{mn}+\mu\delta_{im}\delta_{jn}+\gamma\delta_{in}\delta_{jm}`, ref "4.29". "(⚠️ γ here is a
    material coefficient — not C_p/C_v and not a shear rate.)"
16. `nb.note` — **N44 [B]** "**Symmetry** of σ in i, j:" equation `\gamma=\mu`, ref "4.30". "The subtlety (derivation steps
    6–7): on a *symmetric* S the μ and γ terms both give S_ij, so only the sum μ + γ ever acts, and σ is automatically
    symmetric. 'γ = μ' is the convention that calls the sum 2μ — which is exactly what makes μ the viscosity of Newton's
    law (1.3). (If one also demands that K itself be symmetric in i, j, γ = μ is forced — the book's argument.)"
17. `nb.derivation("D09", …)` — Part F D09 (13 steps), ref "4.31", with `check_src` (★★★, every line commented): build
    K_ijmn from δ's in sympy with symbols λ, μ, γ; contract with a symbolic symmetric S (six symbols); show σ_ij −
    [λS_mmδ_ij + (μ + γ)S_ij] simplifies to 0 for all nine components; show σ − σᵀ = 0; substitute γ = μ and compare with
    2μS_ij + λS_mmδ_ij; finally rotate about z by a symbol θ (rotation matrix) and show K′ − K = 0 component by component
    (spot-check 20 random components to keep the cell < 2 s).
18. `nb.primer("deviatoric (traceless) part of a tensor", "Any second-order tensor splits into an isotropic part (its
    average diagonal times δ) and a traceless remainder: $A_{ij}=\tfrac13A_{mm}\delta_{ij}+\big(A_{ij}-\tfrac13A_{mm}\delta_{ij}\big)$.
    For the strain rate the first part is pure swelling (volume change at rate S_mm = ∇·u) and the second is pure shape
    change at constant volume. Its trace is zero because δ_ii = 3.", code="S = np.array([[1.0, 0.5, 0], [0.5, 2.0, 0], [0,
    0, 3.0]])     # a strain rate [1/s]\nD = ch04.deviatoric_part(S)                    # S − (tr S/3) δ\nprint(np.trace(D), np.trace(S)/3)               # 0.0 and 2.0: the average stretching")`
    — *expect:* `0.0 2.0`.
19. `nb.derivation("D10", …)` — Part F D10 (4 steps), ref "4.37".
20. `nb.note` — **N45 [B]** "**Taking the trace** of (4.31) (set i = j and sum; δ_ii = 3) gives τ_ii = −3p + (2μ + 3λ)S_mm,
    so" equation `p=-\tfrac13\tau_{ii}+\big(\tfrac23\mu+\lambda\big)\nabla\cdot\mathbf u`, ref "4.32". "(S_mm = ∇·u is the
    volumetric strain rate of Ch. 3 §3.4, Eq. (3.14), $\frac{1}{\delta V}\frac{D(\delta V)}{Dt}=\frac{\partial u_i}{\partial x_i}$ —
    the book's pointer 'Section 3.6' should read §3.4.)"
21. `nb.note` — **N46 [B]** "The **mean (mechanical) pressure** is minus the average normal stress:" equation
    `\bar p\equiv-\tfrac13\tau_{ii}`, ref "4.33".
22. `nb.note` — **N47 [B]** "They differ only in an expanding or compressing flow:" equation
    `p-\bar p=\big(\tfrac23\mu+\lambda\big)\nabla\cdot\mathbf u`, ref "4.34". "> ⚠️ **Mean vs thermodynamic pressure.** p comes
    from the equation of state; p̄ is what the normal stresses average to. For an incompressible fluid there is no
    equation of state for p: only a mechanical pressure exists, and only its *gradient* matters — adding a constant to p
    changes nothing in (4.39b) (demonstrated in C08's code)."
23. `nb.note` — **N48 [B]** "**Incompressible** Newtonian stress (S_mm = 0 removes λ):" equation
    `\tau_{ij}=-p\,\delta_{ij}+2\mu S_{ij}`, ref "4.35". "`ch04.newtonian_stress(G, p, mu, incompressible=True)` raises an
    error if tr G ≠ 0 — a guard against using (4.35) on a compressible field."
24. `nb.note` — **N49 [B]** "**Bulk viscosity** $\mu_v=\lambda+\tfrac23\mu$: the resistance to pure expansion. It matters
    for sound absorption and the inside of shock waves (Ch. 15) and is nonzero in polyatomic gases (slow exchange of
    energy with molecular rotation)."
25. `nb.note` — **N50 [B]** "**Stokes' assumption** — the code's default (`mu_v=0.0`):" equation
    `\lambda+\tfrac23\mu=0`, ref "4.36". "Accurate whenever μ_v or the expansion rate is small — nearly always outside
    acoustics and shocks."
26. `nb.note` — **N51 [B]** "**The bulk-viscosity form** (D10's result) separates shape change from volume change:"
    equation `\tau_{ij}=-p\,\delta_{ij}+2\mu\Big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\Big)+\mu_vS_{mm}\delta_{ij}`, ref "4.37".
    "For a parallel flow u = (u(y), 0, 0) it gives τ₁₂ = μ du/dy — Newton's law (1.3), $\tau=\mu\,du/dy$. Number: shear rate
    $\dot\gamma=10$ s⁻¹ in water (μ = 1.0×10⁻³ Pa s) → τ₁₂ = 0.010 Pa, the same as `ch01.newton_shear_stress(1e-3, 10)`."
27. `nb.worked_example("three motions of water (μ = 1.0×10⁻³ Pa s, p = 0 gauge)", "1. **Shear** u = ($\dot\gamma$y, 0, 0),
    $\dot\gamma=10$ s⁻¹: S₁₂ = S₂₁ = 5 s⁻¹, S_mm = 0 → τ₁₂ = 2μS₁₂ = 2 × 10⁻³ × 5 = 0.010 Pa, all normal stresses 0.
    2. **Extension** u = (sx, −sy, 0), s = 1 s⁻¹: S₁₁ = 1, S₂₂ = −1, S_mm = 0 → τ₁₁ + p = +2 mPa (a pull along x),
    τ₂₂ + p = −2 mPa (a push along y); p̄ = p. 3. **Pure expansion** u = (ex, ey, ez), e = 1 s⁻¹: S = diag(1, 1, 1),
    S_mm = 3 → the shear part $S-\tfrac13S_{mm}\delta$ is zero; with Stokes (μ_v = 0) σ = 0 and p̄ = p; with μ_v = 1.0×10⁻³ Pa s,
    every normal stress gains μ_vS_mm = 3 mPa and $p-\bar p=\mu_v\nabla\cdot\mathbf u=3$ mPa. 4. **Rigid spin**: S = 0, σ = 0 —
    no resistance at all.")`
28. `nb.code` — *code:* `G_sh = np.array([[0, 10.0, 0], [0, 0, 0], [0, 0, 0]])` · `G_ext = np.diag([1.0, -1.0, 0.0])` ·
    `G_exp = np.eye(3)` · `for name, G in {"shear": G_sh, "extension": G_ext, "expansion": G_exp}.items():` `tau =
    ch04.newtonian_stress(G, p=0.0, mu=1e-3, mu_v=1e-3)`; `print(name, np.round(tau*1e3, 3), "mPa; p - p_bar =",
    round((0.0 - ch04.mean_pressure(tau))*1e3, 3), "mPa")` (p = 0 gauge) · `print(ch04.stress_on_plane(G_sh, 0.0, 1e-3, 0.0, np.pi/4))` ·
    `print(ch04.bulk_viscosity(-2e-3/3, 1e-3), ch04.lam_from_bulk(0.0, 1e-3))` · `print(ch01.newton_shear_stress(1e-3,
    10.0))`. *expect:* shear: τ₁₂ = τ₂₁ = 10 mPa, p − p̄ = 0 · extension: diag(2, −2, 0) mPa, p − p̄ = 0 · expansion:
    diag(3, 3, 3) mPa, p̄ = −⅓τ_ii = −3 mPa so
    p − p̄ = +3 mPa = μ_v∇·u, (4.34) with ⅔μ + λ = μ_v · `(0.010, 0.0)` (on the 45°
    plane the shear is all normal) · `0.0 -0.000667` · `0.01`. *explain:* 1. the three motions of the worked example
    through (4.37) with μ_v = 10⁻³ Pa s; 2. pressure minus mean pressure is nonzero only for expansion; 3. the traction on
    a 45° plane in shear flow is a pure pull — the principal direction; 4. μ_v and λ convert into each other (Stokes:
    λ = −⅔μ); 5. Ch. 1's Newton law agrees.
29. `nb.check_agree` — **from scratch (curation §7):** explicit loops `tau = np.zeros((3, 3)); S = 0.5*(G + G.T); Smm =
    np.trace(S)`; `for i in range(3): for j in range(3): tau[i, j] = -p*(i == j) + 2*mu*S[i, j] + lam*Smm*(i == j)`
    with λ = μ_v − ⅔μ · `assert np.allclose(tau, ch04.newtonian_stress(G, p, mu, lam=lam))` for a random G · and the
    81-coefficient route: `assert np.allclose(ch04.linear_stress(ch04.isotropic_fourth_order(lam, mu, mu), S), tau + p*np.eye(3))`.
30. `nb.figure` — two panels on one angle axis θ ∈ [0°, 180°] of the plane's normal: (a) shear flow $\dot\gamma=10$ s⁻¹:
    normal stress σ_n = 0.01 sin 2θ (orange) and signed shear τ_s = 0.01 cos 2θ (rose) in Pa with the principal
    directions (45°, 135°) marked; (b) extension s = 1 s⁻¹: σ_n = 0.002 cos 2θ, τ_s = −0.002 sin 2θ — the same curves
    shifted by 45°. Right inset: the parametric curve (σ_n, τ_s) — a circle (Mohr's circle, Ch. 2's gloss). *see:* "two
    sinusoids of period 180°, a quarter-period apart; the shear-flow circle is five times the extension-flow circle."
    *read:* "shear flow *is* an extension along 45° and a compression along 135°, seen on tilted planes: the shape of the
    curves is the same, only rotated by 45° — the stress follows the principal axes of S." *change:* "…the fluid spun
    rigidly on top: G changes (it gains an antisymmetric part), the curves do not — only S matters (R05)."
31. `nb.plotly` — `slider_figure` over μ_v = 0…5×10⁻³ Pa s (11 steps) for pure expansion e = 1 s⁻¹: traces "normal
    stress on a plane + p (orange, flat in θ)" and "p − p̄ (grey line)". Title "Bulk viscosity acts only on expansion".
    see/read/change: "…with μ_v = 0 the orange line sits at 0 (p̄ = p, Stokes); raising μ_v lifts it uniformly: every
    plane feels the same extra pull 3μ_v. Switch the flow to shear (mentally): μ_v would do nothing (S_mm = 0)."
32. `nb.figure` — the D08 cube: log–log α(h) = 6(τ₁₂ − τ₂₁)/(ρh²) for imbalances 0.01, 0.1, 1 Pa (three lines, slope −2)
    with a horizontal band "fastest spin ever seen in a fluid element ~10⁴ rad/s²" for scale. *see:* "straight lines of
    slope −2 climbing to absurd values as h shrinks." *read:* "a stress imbalance at a point would spin the fluid there
    infinitely fast — impossible, so τ₁₂ = τ₂₁ (4.25)." *change:* "…body couples present: a torque per unit *mass* scales
    like h⁵ and could balance the imbalance at every size — then τ need not be symmetric."
33. `nb.explainer("newtonian_stress_lab", heading="How does a fluid decide its stress?", why="Set a velocity gradient with
    presets (shear, extension, rotation, expansion) or sliders and watch the chain G → S (and R) → τ of
    $\tau_{ij}=-p\delta_{ij}+2\mu S_{ij}+\lambda S_{mm}\delta_{ij}$ (4.31) respond together: the element deforms, the traction
    on a plane you rotate swings through its normal and shear values, and the matrices light up by source. A second mode
    shrinks the spinning cube of D08.", tries=["Preset 'pure rotation': G is not zero but every stress entry is — only S
    matters.", "Preset 'expansion' and raise μ_v: p̄ separates from p; switch to 'shear' and μ_v does nothing.", "Rotate
    the plane in shear flow until the shear readout is zero: that is 45°, a principal direction.", "Open the Derivation
    tab at D09 and watch the δ-substitution steps light the matrix entries they produce."])`
34. `nb.md` — **What would change if…** "…we put (4.37) into Cauchy's equation (4.24),
    $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$? The divergence of 2μS becomes μ∇²u plus a term in
    ∇(∇·u), and we get the Navier–Stokes equation — 4 equations for 5 unknowns (C08)."

### A.6 §4.6 Navier–Stokes Momentum Equation — C08 (+N52–N55 · D11, D12, D13 · E4)
1. `nb.section("4.6", "Navier-Stokes Momentum Equation", intro="**What is this section about?** Substitute the Newtonian
   stress into Cauchy's equation: the Navier–Stokes equation. With constant viscosity and ∇·u = 0 it takes its famous
   compact form; its viscous force can be written three ways, and with the viscosity switched off it is Euler's
   equation.")`

**C08 — Incompressible Navier–Stokes (4.39b)**
2. `nb.core("C08", "The Navier–Stokes equation: five accelerations in balance (4.39b)", question="Honey oozing down a pipe
   and air racing over a wing obey the same equation. Which of its terms matter where?")`
3. `nb.md` — **The problem in plain words:** "Weather models, ocean models, aircraft design codes and blood-flow
   simulations all solve one equation. Every exact solution of Ch. 8, every boundary layer of Ch. 9 and every instability
   of Ch. 11 is a flow in which some of its terms balance and the others are asleep. Learning to read the terms — which
   is big, which is zero, and why — is the skill this block builds."
4. `nb.md` — **The idea:**
   ```
   ρ ( ∂u/∂t  +  (u·∇)u )  =  −∇p    +   ρg     +   μ∇²u                        (4.39b)
       local    advective     pressure   gravity    viscous
       blue     teal          orange     grey       rose
   Poiseuille: orange ↔ rose     Stokes' first problem: blue ↔ rose     potential flow: teal ↔ orange (rose = 0)
   ```
5. `nb.derivation("D11", …)` — Part F D11 (7 steps), ref "4.38" (the δ-substitution gloss "δ_ij picks out the term with
   i = j" is in step 2's *why*).
6. `nb.note` — **N52 [B]** "**Navier–Stokes with variable viscosity** (D11's result):" equation
   `\rho\Big(\frac{\partial u_j}{\partial t}+u_i\frac{\partial u_j}{\partial x_i}\Big)=-\frac{\partial p}{\partial x_j}+\rho g_j+\frac{\partial}{\partial x_i}\Big[\mu\Big(\frac{\partial u_j}{\partial x_i}+\frac{\partial u_i}{\partial x_j}\Big)+\Big(\mu_v-\tfrac23\mu\Big)\frac{\partial u_m}{\partial x_m}\delta_{ij}\Big]`,
   ref "4.38". "**Ledger:** continuity (1) + Navier–Stokes (3) = 4 equations; unknowns ρ, p, u_j = 5. Closed when ρ is
   constant or a known function of p alone (a *barotropic* flow): 5 = 5." Followed by `nb.code`: `for s in
   ("navier_stokes", "barotropic"): ledger.loc[s] = ch04.closure_count(s)` · `print(ledger)`. *expect:* rows `4 5` and
   `5 5`.
7. `nb.primer("Schwarz's theorem", "For a smooth function the order of partial derivatives does not matter:
   $\frac{\partial}{\partial x}\frac{\partial f}{\partial y}=\frac{\partial}{\partial y}\frac{\partial f}{\partial x}$. Here it lets us
   move ∂/∂x_j outside ∂u_i/∂x_i, turning a term into the gradient of the divergence ∇(∇·u). It needs continuous second
   derivatives (true for every viscous flow we meet).", code="x, y = sp.symbols('x y')\nf = sp.sin(x*y) + x**3*sp.exp(y)                     # any smooth function\nprint(sp.simplify(sp.diff(f, x, y) - sp.diff(f, y, x)))   # 0")`
   — *expect:* `0`.
8. `nb.derivation("D12", …)` — Part F D12 (7 steps), ref "4.39b" (the vector-Laplacian gloss "in Cartesian components
   ∇²u means ∇² of each component" is in step 5's *why*).
9. `nb.note` — **N53 [B]** "D12's middle line is the **constant-viscosity compressible** form:" equation
   `\rho\frac{Du_j}{Dt}=-\frac{\partial p}{\partial x_j}+\rho g_j+\mu\frac{\partial^2u_j}{\partial x_i^2}+\Big(\mu_v+\tfrac13\mu\Big)\frac{\partial}{\partial x_j}\frac{\partial u_m}{\partial x_m}`,
   ref "4.39a". "Its last term is the only trace of compressibility; ∇·u = 0 removes it and leaves (4.39b),
   $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$."
10. `nb.primer("curl of a curl identity", "For any smooth vector field, $\nabla\times(\nabla\times\mathbf u)=\nabla(\nabla\cdot\mathbf u)-\nabla^2\mathbf u$.
    It follows from the ε–δ identity $\varepsilon_{ijk}\varepsilon_{klm}=\delta_{il}\delta_{jm}-\delta_{im}\delta_{jl}$ (2.19) of
    Ch. 2. With ω = ∇×u and ∇·u = 0 it says ∇²u = −∇×ω: the Laplacian of an incompressible velocity is minus the curl of
    its vorticity.", code="x, y, z = sp.symbols('x y z'); X = sp.Matrix([x, y, z])\nu = sp.Matrix([y**2*z, sp.sin(x)*z, x*y**3])              # any smooth field\ncurl = lambda F: sp.Matrix([sp.diff(F[2], y) - sp.diff(F[1], z), sp.diff(F[0], z) - sp.diff(F[2], x), sp.diff(F[1], x) - sp.diff(F[0], y)])\ndiv = sum(sp.diff(u[i], X[i]) for i in range(3)); lap = sp.Matrix([sum(sp.diff(u[k], v, 2) for v in X) for k in range(3)])\nprint(sp.simplify(curl(curl(u)) - (sp.Matrix([sp.diff(div, v) for v in X]) - lap)))   # zero vector")`
    — *expect:* `Matrix([[0], [0], [0]])`.
11. `nb.derivation("D13", …)` — Part F D13 (9 steps), ref "4.40".
12. `nb.note` — **N54 [B]** "**The viscous force three ways** (incompressible flow):" equation
    `(\mu\nabla^2\mathbf u)_j=\mu\frac{\partial^2u_j}{\partial x_i^2}=2\mu\frac{\partial S_{ij}}{\partial x_i}=\mu\frac{\partial}{\partial x_i}\Big(\frac{\partial u_j}{\partial x_i}+\frac{\partial u_i}{\partial x_j}\Big)=-\mu\varepsilon_{jik}\frac{\partial\omega_k}{\partial x_i}=-\mu(\nabla\times\boldsymbol\omega)_j`,
    ref "4.40". "**The paradox resolved:** the viscous force contains vorticity although rotation was excluded from the
    stress law — but it contains its *derivative*. Solid-body rotation has uniform ω, so ∇×ω = 0 and no viscous force;
    and then S = 0 everywhere too. Viscosity acts only where vorticity (or strain) varies in space — the vorticity
    diffusion of Ch. 5."
13. `nb.note` — **N55 [B]** "**Euler's equation** — viscosity negligible (far from walls):" equation
    `\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g`, ref "4.41". "Check below: the Ch. 3 cylinder potential flow with
    its Bernoulli pressure satisfies it to round-off (and its viscous term is exactly zero anyway, C12)."
14. `nb.primer("complementary error function erfc", "erfc(η) = 1 − erf(η) = (2/√π)∫_η^∞ e^{−s²}ds falls smoothly from 1
    at η = 0 to 0 as η grows (0.48 at η = 0.5, 0.16 at 1, 0.005 at 2). It is the shape of a quantity diffusing into a
    region from a boundary held at a fixed value: here momentum diffusing up from a plate set suddenly in motion.
    `scipy.special.erfc` evaluates it.", code="from scipy.special import erfc      # the complementary error function\nprint(erfc(np.array([0.0, 0.5, 1.0, 2.0])))   # [1. 0.4795 0.1573 0.0047]")`
    — *expect:* `[1.     0.4795 0.1573 0.0047]`.
15. `nb.worked_example("plane Poiseuille flow of water (gap h = 1 mm, G = −dp/dx = 100 Pa/m)", "Walls at y = 0 and y = h,
    μ = 1.0×10⁻³ Pa s. 1. Try $u=\frac{G}{2\mu}y(h-y)$, v = w = 0 (steady, parallel). 2. Local term 0 (steady); advective
    term u ∂u/∂x = 0 (u does not change along x). 3. Pressure force per volume $-\partial p/\partial x=+G=+100$ N/m³ everywhere.
    4. Viscous force $\mu\,\partial^2u/\partial y^2=\mu\times(-G/\mu)=-100$ N/m³ everywhere. 5. Sum 0 ✓ — (4.39b) holds with
    only two awake terms. 6. Peak speed at y = h/2: $u_{max}=Gh^2/(8\mu)=100\times10^{-6}/(8\times10^{-3})=0.0125$ m/s =
    12.5 mm/s.")`
16. `nb.code` — *code:* `u_fn, p_fn = ch04.exact_field("poiseuille", G=100.0, h=1e-3, mu=1e-3)` · `for y in (0.0, 2.5e-4,
    5e-4): print(y, ch04.ns_incompressible_terms(u_fn, p_fn, np.array([0.0, y, 0.0]), 0.0, rho=1000.0, mu=1e-3,
    g=np.zeros(3), h=1e-5, per="volume"))` · `print(ch04.ns_terms_preset("poiseuille", 0.0, 2.5e-4, component=0,
    per="volume", G=100.0, h=1e-3, mu=1e-3))` · constant added to p: `p2 = lambda x, t: p_fn(x, t) + 1e5`;
    `print(ch04.ns_incompressible_terms(u_fn, p2, np.array([0, 2.5e-4, 0]), 0.0, 1000.0, 1e-3, g=np.zeros(3),
    per="volume").residual)` · symbolic (4.38) for a variable viscosity: `y_, mu0, b = sp.symbols('y mu0 b',
    positive=True)`; `print(ch04.navier_stokes_sym(1000, [y_*(1 - y_), 0, 0], -2*sp.Symbol('x'), (sp.Symbol('x'), y_,
    sp.Symbol('z')), sp.Symbol('t'), mu=mu0*(1 + b*y_))[0])` · Euler check: `print(ch04.ns_terms_preset("cylinder", 1.5,
    0.7, component=0, U=1.0, a=1.0, rho=1.0))`. *expect:* at y = 0, 2.5e-4, 5e-4: pressure +100.0, viscous −100.0,
    local/advective 0, residual ≈ 0 (≤ 1e-4 relative) · the same dict from the preset wrapper · residual unchanged ≈ 0
    (only ∇p enters, N47) · a sympy expression showing the extra term from dμ/dy (≠ 0 unless b = 0) · cylinder: local 0,
    viscous 0, advective x ≈ −0.0321 m/s² balanced by the pressure term, residual ≈ 0. *explain:* 1. the Poiseuille field; 2. the five terms at
    three heights; 3. the parity wrapper used by E4; 4. a constant pressure offset changes nothing; 5. with μ(y) the
    viscous term keeps the dμ/dy part — (4.38), not (4.39b); 6. an Euler flow: advective acceleration balanced by the
    pressure gradient, viscous force exactly zero.
17. `nb.check_agree` — **from scratch (curation §7):** the Poiseuille residual by hand: `dy = 1e-6; y0 = 2.5e-4; u =
    lambda y: 100/(2e-3)*y*(1e-3 - y)`; `visc = 1e-3*(u(y0+dy) - 2*u(y0) + u(y0-dy))/dy**2`; `press = 100.0`;
    `assert abs(press + visc) < 1e-3` and `assert np.isclose(visc, ch04.ns_terms_preset("poiseuille", 0.0, y0,
    per="volume", G=100.0, h=1e-3, mu=1e-3)["viscous"], rtol=1e-5)`.
18. `nb.figure` — two panels: (a) the Poiseuille profile u(y) (teal) across the gap; (b) term bars per unit volume at
    five heights: pressure (orange, +100) and viscous (rose, −100), local and advective (blue, teal, zero-height ticks),
    residual dot (black). *see:* "a parabola; at every height the same two equal and opposite bars." *read:* "in fully
    developed pipe/channel flow nothing accelerates: the pressure push is used up by friction at every point." *change:*
    "…G doubled: both bars double and the parabola doubles (u ∝ G) — the balance is linear."
19. `nb.plotly` — `slider_figure` over t = 0.01…10 s (log-spaced, 20 steps, FAST 10) of **Stokes' first problem** (a
    plate suddenly moving at U = 1 m/s in water, ν = 10⁻⁶ m²/s): traces "u(y, t)/U (teal)", "local ∂u/∂t (blue, scaled)"
    and "viscous ν∂²u/∂y² (rose, scaled)" for y ∈ [0, 10] mm, using `ch04.stokes_first_problem`. Title "Momentum
    diffuses from the wall: the local term is paid by viscosity". see/read/change: "…the blue and rose curves lie on top
    of each other at every t (they balance), and the profile's thickness grows like 2√(νt): 2 mm after 1 s, 6.3 mm after
    10 s."
20. `nb.figure` — the paradox of N54: two vortices side by side — solid-body rotation (left) and a Lamb–Oseen vortex
    (right, Γ = 2π m²/s, σ = 1 m, ν = 0.01 m²/s): vorticity as a heat map and the viscous-force arrows μ∇²u from
    `ch04.viscous_force_forms` (zero on the left). *see:* "a uniform vorticity disc with no arrows; a peaked vorticity
    blob with arrows pointing backwards inside the core and forwards outside." *read:* "viscosity slows the fast core and
    speeds the slow outskirts — it spreads the vorticity (Ch. 5) — while a rigid rotation, whose vorticity is uniform,
    feels nothing." *change:* "…ν doubled: the arrows double and the core spreads twice as fast (σ² = 4νt)."
21. `nb.explainer("navier_stokes_term_balance", heading="Which terms of Navier–Stokes are awake here?", why="Pick an
    exact solution — Couette, Poiseuille, Stokes' first problem, Taylor–Green, Lamb–Oseen or the potential flow round a
    cylinder — click any point and read the five terms of
    $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ (4.39b) as bars that add to zero; scrub time for the
    unsteady ones and watch the balance move.", tries=["Poiseuille: click anywhere — only orange and rose, equal and
    opposite.", "Stokes' first problem: press ▶ and watch the blue (local) bar trade with the rose (viscous) bar.",
    "Cylinder (Euler): the viscous bar is exactly zero although μ ≠ 0 — why? (Explain, section 5).", "Open the third
    view: the viscous force computed three ways of (4.40) gives the same arrow."])`
22. `nb.md` — **What would change if…** "…the observer were on a turntable, or on the rotating Earth? Newton's law holds
    only in an inertial frame; in a rotating one extra 'apparent' forces appear — Coriolis and centrifugal (C09). For the
    atmosphere and ocean they are not small."

### A.7 §4.7 Noninertial Frame of Reference — C09 (+N56–N65 · D14–D18 · E5)
1. `nb.section("4.7", "Noninertial Frame of Reference", intro="**What is this section about?** The equations so far hold
   for an observer who is not accelerating. We rewrite Navier–Stokes for an observer who accelerates and rotates — the
   frame every geophysical flow is described in — and meet the Coriolis and centrifugal terms.")`

**C09 — Navier–Stokes in a translating, rotating frame (4.45)**
2. `nb.core("C09", "Seen from a rotating frame: Coriolis and centrifugal forces (4.45)", question="Roll a ball straight
   across a spinning merry-go-round. Someone standing on the ground sees a straight line; the rider sees it curve. Who is
   right — and what force does the rider need to invent?")`
3. `nb.md` — **The problem in plain words:** "We live on a rotating planet and measure winds and currents relative to the
   ground. For a thrown ball the rotation hardly matters; for a hurricane or an ocean gyre that lasts days it dominates.
   Turbomachines and centrifuges are also analysed in frames that turn with them. We need Navier–Stokes as seen by an
   observer who translates with acceleration dU/dt and rotates at Ω(t)."
4. `nb.md` — **The idea:**
   ```
   inertial observer:  a = F/m                       rotating observer: a' = F/m + (apparent forces)
   the rotating axes e'_i turn:  de'_i/dt = Ω × e'_i   →  velocity gains Ω × x'      (4.42)
   differentiate again           →  a = dU/dt + a' + 2Ω×u' + dΩ/dt×x' + Ω×(Ω×x')      (4.43)
                                       frame   rel.  Coriolis angular  centripetal
   move them to the force side   →  −dU/dt, −2Ω×u' (amber), −dΩ/dt×x', −Ω×(Ω×x') (purple)   in (4.45)
   ```
5. `nb.note` — **N56 [C]** "An **inertial frame** is one in which Newton's law holds as written — fixed to the distant
   stars, or moving at constant velocity relative to them. A laboratory on the Earth is nearly inertial over seconds and
   metres; over hours and hundreds of kilometres it is not (Ch. 13)."
6. `nb.md` — "> ⚠️ **The prime changes meaning again.** In this section x′, u′, a′ and D′/Dt are measured in the
   **rotating** frame O′1′2′3′ (Ch. 3's primes were a translating frame; §4.9's p′, ρ′ will be perturbations)."
7. `nb.primer("derivative of a rotating unit vector", "A unit vector e′ fixed to a frame turning at angular velocity Ω
   keeps its length but changes direction: its tip moves on a circle round the Ω axis. In time dt it moves a distance
   (sin α)|Ω| dt (α the angle between e′ and Ω), perpendicular to both — exactly the cross product: $d\mathbf e'/dt=\boldsymbol\Omega\times\mathbf e'$.
   Here it is why a rotating observer's axes add terms to every time derivative.", code="Om = np.array([0.0, 0.0, 0.5])
    # 0.5 rad/s about z\nE = ch04.rotating_basis(Om, 1.0); dE = ch04.basis_rate(Om, 1.0)   # rows e'_i(t) and d e'_i/dt\nprint(np.allclose(dE, np.cross(Om, E)))                  # True: de'/dt = Ω × e'")`
   — *expect:* `True`.
8. `nb.primer("product rule for a cross product", "$\frac{d}{dt}(\mathbf a\times\mathbf b)=\dot{\mathbf a}\times\mathbf b+\mathbf a\times\dot{\mathbf b}$ —
   like the ordinary product rule, but the order of the factors must be kept (a × b = −b × a).", code="t =
   sp.symbols('t')\na = sp.Matrix([sp.cos(t), t, 1]); b = sp.Matrix([t**2, 0, sp.sin(t)])\nprint(sp.simplify(a.cross(b).diff(t) - (a.diff(t).cross(b) + a.cross(b.diff(t)))))   # zero")`
   — *expect:* `Matrix([[0], [0], [0]])`.
9. `nb.derivation("D14", …)` — Part F D14 (7 steps), ref "4.42".
10. `nb.note` — **N57 [B]** "**Velocity seen from the inertial frame** (D14's result):" equation
    `\mathbf u=\mathbf U+\mathbf u'+\boldsymbol\Omega\times\mathbf x',\qquad\frac{d\mathbf e'_i}{dt}=\boldsymbol\Omega\times\mathbf e'_i`,
    ref "4.42". "u′ is how fast the *components* x′_i change — what the rotating observer measures."
11. `nb.derivation("D15", …)` — Part F D15 (12 steps), ref "4.43", with `check_src` (★★★, every line commented): sympy
    rotation about z by an angle θ(t) (so Ω = θ̇ e_z, Ω̇ = θ̈ e_z) — build x(t) = X(t) + R(t)x′(t) with generic functions
    X_i(t), x′_i(t), θ(t); differentiate twice; rotate the result back into the primed components (Rᵀ ẍ); subtract the
    five terms of (4.43) written in primed components (dU/dt, ẍ′, 2Ω × ẋ′, Ω̇ × x′, Ω × (Ω × x′)); `simplify` each
    component to 0. A second line shows that with a coefficient 1 instead of 2 on Ω × u′ the difference is Ω × u′ ≠ 0.
12. `nb.note` — **N58 [B]** "**Acceleration in the inertial frame**, term by term:" equation
    `\mathbf a=\frac{d\mathbf U}{dt}+\mathbf a'+2\boldsymbol\Omega\times\mathbf u'+\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'+\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')`,
    ref "4.43". "Names: frame acceleration, acceleration seen in the rotating frame, **Coriolis** term, angular-acceleration
    term, **centripetal** term. The factor 2: one Ω × u′ comes from the turning axes (differentiating u′'s basis), one
    from the moving position (differentiating Ω × x′)."
13. `nb.note` — **N59 [B]** "For a fluid particle, a = Du/Dt in the inertial frame and a′ = D′u′/Dt in the rotating one:"
    equation `\Big(\frac{D\mathbf u}{Dt}\Big)_{O123}=\Big(\frac{D'\mathbf u'}{Dt}\Big)_{O'1'2'3'}+\frac{d\mathbf U}{dt}+2\boldsymbol\Omega\times\mathbf u'+\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'+\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')`,
    ref "4.44".
14. `nb.derivation("D16", …)` — Part F D16 (6 steps), ref "4.45".
15. `nb.md` — "> ⚠️ **Coriolis *term* vs Coriolis *force*.** On the acceleration side, (4.43) has **+2Ω × u′**; moved to the
    force side, (4.45) has **−2Ω × u′** per unit mass. The book calls both 'Coriolis acceleration'. We say *Coriolis term*
    for the first and *Coriolis force per unit mass* for the second; the second is the one that deflects a moving parcel
    to the **right** in the Northern Hemisphere (Ω up). Likewise Ω × (Ω × x′) is *centripetal* (points to the axis) and
    −Ω × (Ω × x′) *centrifugal* (away from it). In code: `frame_acceleration_terms` returns the + terms,
    `apparent_body_forces` the − terms."
16. `nb.note` — **N60 [B]** "**Frame acceleration** −dU/dt: the push back into your seat when a car accelerates. An
    aircraft on a parabolic arc has dU/dt = g, so g − dU/dt = 0: weightlessness." Code line:
    `print(ch04.apparent_body_forces(np.zeros(3), np.zeros(3), np.zeros(3), dU_dt=np.array([0, 0, -9.81]))["total"])` →
    `[0. 0. 0.]`.
17. `nb.primer("latitude, Earth's rotation rate and the local vertical", "The Earth turns once per sidereal day (86 164 s),
    Ω = 2π/86164 = 7.292×10⁻⁵ rad/s, about the axis through the poles. At latitude φ the local vertical makes an angle
    with that axis, so Ω splits into a vertical part Ω sin φ (full at the poles, zero at the equator) and a horizontal
    part Ω cos φ along the meridian. For horizontal winds the vertical part does the steering — twice it, f = 2Ω sin φ, is
    the *Coriolis parameter* of Ch. 13 (named here only).", code="Om = 2*np.pi/86164.0\nfor lat in (0, 30, 45, 90):  # degrees\n    print(lat, round(Om*np.sin(np.deg2rad(lat)), 9), round(ch04.coriolis_parameter(np.deg2rad(lat)), 9))")`
    — *expect:* Ω = 7.2921e-5; at 45°: Ω sin φ = 5.156e-5, f = 1.031e-4 s⁻¹.
18. `nb.derivation("D17", …)` — Part F D17 (5 steps), ref "".
19. `nb.note` — **N61 [B]** "**The Coriolis force per unit mass −2Ω × u′** depends on the velocity, not the position; it
    is perpendicular to u′, so it changes the direction but never the speed (it does no work). A projectile fired
    horizontally from the North Pole at speed u is deflected sideways by Ωut² after time t (D17), an angle Ωt — exactly
    the Earth's rotation in that time: seen from space the path is straight and the ground turned under it. Number:
    u = 10 m/s, t = 1 h: Coriolis acceleration 2Ωu = 1.46×10⁻³ m/s², forward 36 km, deflection Ωut² = 9.45 km (the exact
    offset ut sin Ωt = 9.34 km), angle Ωt = 15.0°."
20. `nb.note` — **N62 [B]** "**Highs and lows.** In cylindrical coordinates about a pressure centre, with Ω_z > 0 (NH),
    flow leaving a high (u_R > 0) feels −2Ω × u = −2Ω_z u_R e_φ: clockwise. Flow into a low (u_R < 0) is turned
    counter-clockwise. In the Southern Hemisphere Ω_z < 0 and both reverse. (The steady wind round a high or low —
    geostrophic and gradient wind — is Ch. 13.)" equation `-2\boldsymbol\Omega\times\mathbf u=-2\Omega_zu_R\,\mathbf e_\varphi`, no ref.
21. `nb.note` — **N63 [C]** "**Angular-acceleration term** −(dΩ/dt) × x′: only when the rotation rate or axis changes (a
    spinning-up centrifuge); zero for the Earth. One bar in E5."
22. `nb.derivation("D18", …)` — Part F D18 (5 steps), ref "".
23. `nb.note` — **N64 [B]** "**Centrifugal term and effective gravity.** For steady rotation about z,
    $-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')=\Omega^2R\,\mathbf e_R$ — outward, growing with the distance R from the
    axis — and it is minus the gradient of $-\tfrac12\Omega^2R^2$ (D18). It adds to gravitation: $\mathbf g_e=\mathbf g_n+\Omega^2R\,\mathbf e_R$
    is the 'gravity' a plumb line shows, and the sea surface is an equipotential of $\Phi_e=gz-\tfrac12\Omega^2R^2$ — the
    Earth bulges at the equator by 2(a − b) = 42.77 km (WGS-84; the book's '42 km' is rounded; its Fig. 4.9 caption says
    'budge' for bulge). Number: at the equator Ω²a = (7.292×10⁻⁵)² × 6.378×10⁶ = 0.0339 m/s², 0.35 % of g. Meteorologists
    fold it into g and the geopotential (Ch. 13)."
24. `nb.note` — **N65 [B]** "**Ex. 4.5 — a von Kármán viscous pump**: a disc rotating under fluid, written in the frame
    that turns with the disc, in cylindrical components (Appendix B operators, glossed: curvilinear coordinates add
    terms such as −u_φ²/R because the unit vectors turn). The rotation terms are ρ[2Ω_z u_φ + Ω_z²R] in the radial
    equation and ρ[−2Ω_z u_R] in the azimuthal one; the axial equation has none. The cell prints the three equations from
    `ch04.rotating_pump_equations()` with the rotation terms highlighted. Rotating flows between discs and cylinders:
    Ch. 8."
25. `nb.worked_example("a projectile at the North Pole", "Ω = 7.292×10⁻⁵ rad/s (up), u = 10 m/s along +x, no friction.
    1. Coriolis force per mass $-2\boldsymbol\Omega\times\mathbf u'$: Ω = Ωe_z, u′ = u e_x, e_z × e_x = e_y, so it is
    −2Ωu e_y = −1.458×10⁻³ e_y m/s² — to the right of the motion. 2. In t = 3600 s the ball goes ut = 36 000 m forward.
    3. Sideways, with a nearly constant acceleration: ½(2Ωu)t² = Ωut² = 7.292×10⁻⁵ × 10 × 3600² = 9451 m. 4. Angle
    9451/36000 = 0.2625 rad = 15.0° = Ωt. 5. Exact (straight inertial line seen from the turning frame):
    ut sin Ωt = 9342 m — the small-angle estimate is 1.2 % high after an hour. 6. Centrifugal at the equator for scale:
    Ω²a = 0.0339 m/s².")`
26. `nb.code` — *code:* `Om = np.array([0, 0, ch04.OMEGA_EARTH])` · `terms = ch04.frame_acceleration_terms(np.zeros(3),
    np.array([10.0, 0, 0]), np.array([1000.0, 0, 0]), Om)`; `print({k: np.round(v, 8) for k, v in terms.items()})` ·
    `forces = ch04.apparent_body_forces(np.array([10.0, 0, 0]), np.array([1000.0, 0, 0]), Om)`; `print(forces["coriolis"],
    forces["centrifugal"])` · `pr = ch04.coriolis_projectile(10.0, ch04.OMEGA_EARTH, np.array([600.0, 1800.0, 3600.0]))`;
    `print(pr["forward"], pr["deflection_small"], pr["deflection"], np.rad2deg(pr["angle"]))` · `print(
    ch04.centrifugal_acceleration(Om, np.array([ch04.EARTH_RADIUS, 0, 0])), ch04.effective_gravity(np.deg2rad(45.0)))`.
    *expect:* coriolis term (0, +1.458e-3, 0) (the + sign of (4.43)), centripetal (−5.3e-6, 0, 0) · force −2Ω × u′ =
    (0, −1.458e-3, 0), centrifugal (+5.3e-6, 0, 0) · forward [6000, 18000, 36000] m, small [262.5, 2362.6, 9450.6] m,
    exact [262.4, 2355.9, 9342.4] m, angle [2.51, 7.52, 15.04]° · (0.0339, 0, 0) m/s² · g_e ≈ 9.783 m/s² and deflection
    ≈ 1.73×10⁻³ rad at 45° (with g_n = 9.80). *explain:* 1. the five acceleration terms of (4.43) for a particle moving at
    10 m/s, 1 km from the axis; 2. the same as forces per mass in (4.45) — every sign flipped; 3. the projectile after 10,
    30, 60 min: the small-angle Ωut² and the exact offset; 4. the centrifugal acceleration at the equator and effective
    gravity at 45°.
27. `nb.check_agree` — **from scratch (curation §7):** the five terms with `np.cross` by hand (`2*np.cross(Om, u_rot)`,
    `np.cross(Om, np.cross(Om, x_rot))`) vs `frame_acceleration_terms` · and the independent route: sample a straight
    inertial path x(t) = (u₀t, 0, 0), express it in the turning frame x′(t) = Rᵀ(t)x(t) with `ch04.rotating_basis`,
    difference x′ twice in time (central differences, dt = 1 s) to get a′, and check
    $\mathbf a'=-2\boldsymbol\Omega\times\mathbf u'-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')$ (the inertial a = 0 in (4.43))
    to 1e-6 relative.
28. `nb.animation` — (`player="video"`, 60 frames, FAST 30) two panels on one clock, t = 0…6 h (Ω exaggerated ×10 so the
    curve is visible, stated in the title): left, the inertial view — a polar cap (disc) turning under a straight rose
    path with a meridian line painted on the disc; right, the rotating view — the same path curving to the right, the
    Coriolis arrow (amber) always perpendicular to it, the small-angle parabola Ωut² as a grey dashed ghost. see/read/
    change: *see:* "left: a straight line over a turning floor; right: a curve bending right." *read:* "the two pictures
    are one motion; the rotating observer needs the amber arrow −2Ω × u′ to explain the curve." *change:* "…the Southern
    Hemisphere (Ω pointing down relative to the local vertical): the right panel bends left."
29. `nb.figure` — **N62** quiver: radial outflow from a high (u_R = 1 m/s at R = 100 km, `ch04.high_low_flow(...,
    sense="high")`, grey arrows) and the Coriolis force per mass on each (amber arrows, `ch04.coriolis_acceleration`);
    second panel the same for a low. *see:* "grey arrows pointing out (high) or in (low); amber arrows at right angles."
    *read:* "the amber arrows turn outflow clockwise and inflow counter-clockwise in the NH — the sense of the winds round
    highs and lows on a weather map." *change:* "…Southern Hemisphere: every amber arrow reverses."
30. `nb.plotly` — **N64** 3-D (plotly): a sphere with Ω exaggerated so the effect is visible, arrows at latitudes 0°, 30°,
    60°, 90° for g_n (blue, toward the centre), Ω²R (purple, away from the axis) and g_e (black), plus the flattened
    equipotential as a translucent spheroid; with a 2-D inset of the true (unexaggerated) deflection angle vs latitude from
    `ch04.effective_gravity` (peak ≈ 0.1° at 45°). see/read/change: "…plumb lines do not point at the Earth's centre
    except at the poles and the equator."
31. `nb.explainer("rotating_frame_coriolis", heading="Why does a straight throw curve on a merry-go-round?", why="One
    ball, one clock, two observers: the inertial view (straight path, table turning) beside the rotating view (curved
    path) with the apparent forces of (4.45),
    $\rho\frac{D'\mathbf u'}{Dt}=-\nabla'p+\rho\big[\mathbf g-\frac{d\mathbf U}{dt}-2\boldsymbol\Omega\times\mathbf u'-\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')\big]+\mu\nabla'^2\mathbf u'$,
    drawn as arrows and bars. Modes: turntable, Earth's pole, flow out of a high, effective gravity.", tries=["Preset
    'Earth pole 1 h': compare the deflection readout with Ωut² = 9.45 km.", "Set Ω = 0: the rotating view becomes the
    inertial one — no amber arrow.", "Mode 'flow out of a high': switch to the Southern Hemisphere and watch the turning
    reverse.", "Derivation tab, D15 steps 6 and 9: each adds one half of the Coriolis arrow."])`
32. `nb.md` — **What would change if…** "…we asked not for forces but for energy? Every force above does work at some rate
    — except Coriolis, which is always perpendicular to the motion. The next block follows the energy of a fluid particle
    and finds where viscosity sends it (C10)."

### A.8 §4.8 Conservation of Energy — R07, R08, C10 (+N66–N81 · D19–D23 · E7)
1. `nb.section("4.8", "Conservation of Energy", intro="**What is this section about?** The first law of thermodynamics
   for a moving fluid particle. We derive the total-energy equation, split off its mechanical part (Cauchy's equation
   dotted with u) and are left with an equation for the internal energy whose viscous term — the dissipation ε — can
   only be positive: viscosity turns motion into heat and never the reverse.")`
2. `nb.recap("R07", "A symmetric tensor contracted with an antisymmetric one gives zero", "For symmetric σ and
   antisymmetric A, $\sigma_{ij}A_{ij}=0$ (Ch. 2 §2.10: swap the dummy indices, use σ_ji = σ_ij and A_ji = −A_ij, so the sum
   equals minus itself). Chapter 2 flagged that this is exactly what makes the viscous work depend only on the strain
   rate — used in step 6 of D22 below.", where="Ch. 2 §2.10")` — followed by `nb.code`: `rng = np.random.default_rng(0)`;
   `B = rng.normal(size=(3, 3)); sig = B + B.T; A = B - B.T` · `print(tensors.symmetric_double_contraction(sig, A))` →
   ~1e-16.
3. `nb.recap("R08", "The Gibbs relation", "For a substance in equilibrium $T\,ds=de+p\,dv$ *(Eq. 1.18)*, with v = 1/ρ the
   specific volume. Applied along the path of a fluid particle (each particle is in local equilibrium) it reads
   $\frac{De}{Dt}=T\frac{Ds}{Dt}-p\frac{D(1/\rho)}{Dt}$ *(Eq. 4.61)* — the bridge from energy to entropy.", where="Ch. 1
   §1.8")`

**C10 — The internal-energy equation (4.57) and viscous dissipation (4.58)**
4. `nb.core("C10", "Where the energy goes: internal energy and viscous dissipation (4.57)", question="Stir a cup of
   coffee and let go. The swirl dies away. Its kinetic energy has to go somewhere — where, and can it ever come back?")`
5. `nb.md` — **The problem in plain words:** "Energy is conserved, but not every kind of energy is equally useful. The
   swirl's kinetic energy ends up as a (tiny) warming of the coffee; the reverse — coffee spontaneously starting to swirl
   while cooling — never happens. In the atmosphere the same one-way conversion ends the kinetic-energy cycle: winds are
   driven by heating differences and die by friction, which returns the energy as heat. Lubricated bearings get hot for
   the same reason, and in turbulence (Ch. 12) the rate of this conversion, ε, is the most important single number."
6. `nb.md` — **The idea:**
   ```
   total energy e + ½u²:  rate = work of gravity + work of surface stresses − heat conducted out        (4.53)
   kinetic ½u² alone:     Cauchy (4.24) · u  = work of gravity + work of the NET pressure and viscous forces (4.56)
   subtract:  internal e:  De/Dt = −p Dv/Dt  +  ε  −  (1/ρ)∇·q                                          (4.57)
                                   compression  dissipation ≥ 0  conduction
   ```
   "Stress does two kinds of work: *force work* (the net force times the velocity) speeds particles up or slows them
   down; *deformation work* (stress times the rate of deformation) changes their shape — and the viscous part of that
   is lost to heat."
7. `nb.primer("power of a force and heat flux through a surface", "A force F acting on something moving at velocity u
   does work at the rate F·u (watts). Per unit area of a surface the stress does work at f·u; per unit volume gravity
   does work at ρg·u. Heat crossing a surface with outward normal n at flux q (W/m²) leaves at the rate q·n per unit
   area — so a *minus* sign appears when we count heat *gained*.", code="f = np.array([2.0, 0.0, 0.0]); u =
   np.array([3.0, 4.0, 0.0])       # traction [Pa], velocity [m/s]\nprint(f @ u)                                     # 6.0 W/m²: only the part of u along f counts\nq = np.array([0.0, 0.0, 50.0]); n = np.array([0, 0, 1.0])\nprint(-q @ n)                                    # −50 W/m²: heat leaving through a lid")`
   — *expect:* `6.0` · `-50.0`.
8. `nb.note` — **N66 [B]** "**First law for a material volume** (Ch. 1's first law (1.10) per unit time: heat in + work
   done = rise of internal + kinetic energy):" equation
   `\frac{d}{dt}\int_{V(t)}\rho\big(e+\tfrac12\lvert\mathbf u\rvert^2\big)dV=\int_{V(t)}\rho\mathbf g\cdot\mathbf u\,dV+\int_{A(t)}\mathbf f\cdot\mathbf u\,dA-\int_{A(t)}\mathbf q\cdot\mathbf n\,dA`,
   ref "4.46".
9. `nb.derivation("D19", …)` — Part F D19 (8 steps), ref "4.53".
10. `nb.note` — **N67 [B]** "D19's step 2 is (4.47), the transport theorem with F = ρ(e + ½|u|²), b = u:" equation
    `\int_{V}\frac{\partial}{\partial t}\Big(\rho e+\frac\rho2\lvert\mathbf u\rvert^2\Big)dV+\int_{A}\Big(\rho e+\frac\rho2\lvert\mathbf u\rvert^2\Big)(\mathbf u\cdot\mathbf n)dA=\int_{V}\rho\mathbf g\cdot\mathbf u\,dV+\int_{A}\mathbf f\cdot\mathbf u\,dA-\int_{A}\mathbf q\cdot\mathbf n\,dA`,
    ref "4.47".
11. `nb.note` — **N68 [B]** "The same coincidence move as D01 and D05 gives the **energy budget of any control volume**
    (stated, not re-derived):" equation
    `\frac{d}{dt}\int_{V^*}\rho\big(e+\tfrac12\lvert\mathbf u\rvert^2\big)dV+\int_{A^*}\Big(\rho e+\frac\rho2\lvert\mathbf u\rvert^2\Big)(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\cdot\mathbf u\,dV+\int_{A^*}\mathbf f\cdot\mathbf u\,dA-\int_{A^*}\mathbf q\cdot\mathbf n\,dA`,
    ref "4.48". "`ch04.energy_budget` evaluates it (used for nozzles and shocks in Ch. 15)."
12. `nb.note` — **N69, N70, N71, N72 [B]** "Steps 3–6 of D19 are (4.49)–(4.52): Gauss on the energy flux, on the stress work and on
    the heat flux, then everything under one integral. The heat-flux line is" equation
    `\int_A\mathbf q\cdot\mathbf n\,dA=\int_Aq_in_i\,dA=\int_V\nabla\cdot\mathbf q\,dV=\int_V\frac{\partial q_i}{\partial x_i}dV`,
    ref "4.51". "> ⚠️ **Book slip:** (4.51) prints dA inside its two volume integrals; Gauss' theorem gives dV (as (4.52)
    then correctly uses)."
13. `nb.note` — **N73 [B]** "**Total-energy equation** (D19's result), the conservative form used by compressible codes
    (Ch. 10, 15):" equation
    `\frac{\partial}{\partial t}\Big(\rho\Big[e+\tfrac12u_j^2\Big]\Big)+\frac{\partial}{\partial x_i}\Big(\rho\Big[e+\tfrac12u_j^2\Big]u_i\Big)=\rho g_iu_i+\frac{\partial}{\partial x_i}(\tau_{ij}u_j)-\frac{\partial q_i}{\partial x_i}`,
    ref "4.53". "u_j² = u₁² + u₂² + u₃² (summed)."
14. `nb.derivation("D20", …)` — Part F D20 (7 steps), ref "4.55".
15. `nb.note` — **N74 [B]** "**Stress work = deformation work + force work** (D20's step 2):" equation
    `\frac{\partial}{\partial x_i}(\tau_{ij}u_j)=\tau_{ij}\frac{\partial u_j}{\partial x_i}+u_j\frac{\partial\tau_{ij}}{\partial x_i}=\Big(-p\frac{\partial u_j}{\partial x_j}+\sigma_{ij}\frac{\partial u_j}{\partial x_i}\Big)+\Big(-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}\Big)`,
    ref "4.54". "The first bracket deforms particles (and heats them), the second changes their kinetic energy. The same
    split, applied to fluctuations, is the turbulent kinetic-energy budget of Ch. 12. Code: `ch04.stress_work_split`
    adds the four parts back to ∂(τ_ij u_j)/∂x_i on a test field."
16. `nb.note` — **N75 [B]** "**Total energy following a particle** (D20's result; the book leaves it to an exercise):"
    equation `\rho\frac{D}{Dt}\Big(e+\tfrac12u_j^2\Big)=\rho g_iu_i+\Big(-p\frac{\partial u_j}{\partial x_j}+\sigma_{ij}\frac{\partial u_j}{\partial x_i}\Big)+\Big(-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}\Big)-\frac{\partial q_i}{\partial x_i}`,
    ref "4.55".
17. `nb.primer("chain rule for the kinetic energy", "D/Dt obeys the ordinary chain rule, so $\frac{D}{Dt}\big(\tfrac12u_j^2\big)=u_j\frac{Du_j}{Dt}$
    (sum over j): the rate of change of kinetic energy per unit mass is the velocity dotted with the acceleration.",
    code="t = sp.symbols('t'); u = [sp.Function(f'u{j}')(t) for j in (1, 2, 3)]\nprint(sp.simplify(sp.diff(sum(v**2 for v in u)/2, t) - sum(v*sp.diff(v, t) for v in u)))   # 0")`
    — *expect:* `0`. (A one-line reminder of Ch. 1's chain rule P49; kept as a primer because D21 hinges on it.)
18. `nb.derivation("D21", …)` — Part F D21 (6 steps), ref "4.56".
19. `nb.note` — **N76 [B]** "**Mechanical-energy equation** — Cauchy's equation dotted with u:" equation
    `\rho\frac{D}{Dt}\Big(\tfrac12u_j^2\Big)=\rho g_ju_j-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}`,
    ref "4.56". "> ⚠️ The book says 'multiply (4.22) by u_j'; (4.22) is the flux form, from which you must also subtract
    u_j × continuity. Starting from Cauchy's (4.24), $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$, is
    one step shorter (D21 does this). Mean and turbulent kinetic energy: Ch. 12."
20. `nb.derivation("D22", …)` — Part F D22 (8 steps), ref "4.57" (the quotient-rule gloss
    $D(1/\rho)/Dt=-\rho^{-2}D\rho/Dt$ is in step 4's *why*).
21. `nb.note` — **C10's equation, restated.** "**Internal energy** — the thermal half of the energy budget:" equation
    `\frac{De}{Dt}=-p\frac{Dv}{Dt}+\frac1\rho\sigma_{ij}S_{ij}-\frac1\rho\frac{\partial q_i}{\partial x_i}`, ref "4.57". "It is the
    first law of Ch. 1 for a particle, per unit time: compression work −p dv, viscous heating, heat conduction. Only the
    first and last can change sign."
22. `nb.primer("completing the square for tensors", "The double sum $A_{ij}A_{ij}=\sum_{i,j}A_{ij}^2$ is a sum of squares,
    so it is ≥ 0 and zero only when every component is. To show an expression is never negative, rewrite it as such
    sums: here $S_{ij}S_{ij}-\tfrac13S_{mm}^2=\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)$ (expand
    the right side and use δ_ijδ_ij = 3).", code="S = np.random.default_rng(3).normal(size=(3, 3)); S = (S + S.T)/2   # a random symmetric S\nD = S - np.trace(S)/3*np.eye(3)                      # its deviatoric part\nprint(np.isclose((S*S).sum() - np.trace(S)**2/3, (D*D).sum()), (D*D).sum() >= 0)   # True True")`
    — *expect:* `True True`.
23. `nb.derivation("D23", …)` — Part F D23 (8 steps), ref "4.58".
24. `nb.note` — **N77 [B]** "**The dissipation rate** — kinetic energy turned into heat per unit mass and time:" equation
    `\varepsilon\equiv\frac1\rho\sigma_{ij}S_{ij}=2\nu\Big(S_{ij}-\tfrac13\frac{\partial u_m}{\partial x_m}\delta_{ij}\Big)^2+\frac{\mu_v}{\rho}\Big(\frac{\partial u_m}{\partial x_m}\Big)^2\ge0`,
    ref "4.58". "Squares with positive coefficients: ε ≥ 0 whenever μ ≥ 0 and μ_v ≥ 0. Incompressible: ε = 2νS_ijS_ij.
    It grows with the square of the velocity gradients — hot bearings, glowing re-entry shields, and turbulence (Ch. 12)."
25. `nb.note` — **N78 [B]** "The **Newtonian viscous stress** used there (the σ part of (4.37)):" equation
    `\sigma_{ij}=\mu\Big(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\Big)+\Big(\mu_v-\tfrac23\mu\Big)\frac{\partial u_m}{\partial x_m}\delta_{ij}`,
    ref "4.59".
26. `nb.note` — **N79 [B]** "With (4.58) and Fourier's law q = −k∇T (Ch. 1 §1.5, $\mathbf q=-k\nabla T$ (1.2)) the internal
    energy equation becomes" equation
    `\rho\frac{De}{Dt}=-p\frac{\partial u_m}{\partial x_m}+2\mu\Big(S_{ij}-\tfrac13\frac{\partial u_m}{\partial x_m}\delta_{ij}\Big)^2+\mu_v\Big(\frac{\partial u_m}{\partial x_m}\Big)^2+\frac{\partial}{\partial x_i}\Big(k\frac{\partial T}{\partial x_i}\Big)`,
    ref "4.60". "**Ledger closes:** continuity (1) + Navier–Stokes (3) + energy (1) + two equations of state = 7; unknowns
    ρ, e, p, T, u_j = 7." Followed by `nb.code`: `ledger.loc["full"] = ch04.closure_count("full")`; `print(ledger)`.
    *expect:* the full table: cauchy 6/13, navier_stokes 4/5, barotropic 5/5, full 7/7.
27. `nb.note` — **N80 [B]** "Combining (4.57) with the Gibbs relation (R08) and the quotient rule gives the **entropy
    equation** (stated; the moves are the book's):" equation
    `\frac{Ds}{Dt}=-\frac1\rho\frac{\partial}{\partial x_i}\Big(\frac{q_i}{T}\Big)-\frac{q_i}{\rho T^2}\frac{\partial T}{\partial x_i}+\frac{\varepsilon}{T}`,
    ref "4.62".
28. `nb.note` — **N81 [B]** "The last two terms are the **entropy production**; with q = −k∇T:" equation
    `\frac{k}{\rho T^2}\lvert\nabla T\rvert^2+\frac\varepsilon T\ge0`, ref "4.63". "The second law demands it for every flow,
    so **μ ≥ 0, μ_v ≥ 0, k ≥ 0**. > ⚠️ The book writes 'μ, κ, k > 0': κ there is a leftover symbol for the bulk viscosity
    μ_v (κ is the thermal diffusivity elsewhere), and the law requires ≥, not >. An inviscid, non-conducting flow has
    Ds/Dt = 0: every particle keeps its entropy (isentropic) — used for (4.78) and in Ch. 15."
29. `nb.worked_example("water sheared between plates (plane Couette, U = 1 m/s, gap h = 1 mm)", "μ = 1.0×10⁻³ Pa s,
    ρ = 1000 kg/m³, k = 0.6 W/(m K), both walls held at T₀. 1. Velocity u = Uy/h, shear rate du/dy = 1000 s⁻¹. 2. S₁₂ =
    500 s⁻¹, so ε = 2ν(S₁₂² + S₂₁²) = 2 × 10⁻⁶ × 2 × 500² = 1.0 W/kg, ρε = μ(U/h)² = 1000 W/m³. 3. Work done by the moving
    wall per unit area: τU = μ(U/h)U = 1.0 × 1 = 1.0 W/m²; heat generated in the gap: ρε × h = 1000 × 10⁻³ = 1.0 W/m² —
    the same (all the wall's work becomes heat). 4. Steady temperature: k T″ = −ρε ⇒ a parabola, peak rise at mid-gap
    $\Delta T_{max}=\rho\varepsilon h^2/(8k)=\mu U^2/(8k)=10^{-3}/4.8=2.1\times10^{-4}$ K; each wall conducts away 0.5 W/m².
    5. Entropy production at mid-gap: ε/T = 1/293 = 3.4×10⁻³ W/(kg K) > 0.")`
30. `nb.code` — *code:* `G = np.array([[0, 1000.0, 0], [0, 0, 0], [0, 0, 0]])` · `print(ch04.dissipation_rate(G, 1000.0,
    1e-3, form="both"))` · `c = ch04.couette_heating(np.linspace(0, 1e-3, 5), U=1.0, h=1e-3, mu=1e-3, k=0.6, T0=293.15)`;
    `print(c["eps"], c["T"] - 293.15, c["work_in"], c["heat_out"], c["q_bottom"], c["q_top"], c["dT_max"])` · `print(
    ch04.entropy_production(np.zeros(3), 293.15, 0.6, 1000.0, 1.0))` · the second-law demo: `print(ch04.dissipation_rate(G,
    1000.0, -1e-3))` (a "negative viscosity" gives ε < 0 — forbidden). *expect:* `(1.0, 1.0)` W/kg · ρε = 1000 W/m³ at
    every y; ΔT = [0, 1.56e-4, 2.08e-4, 1.56e-4, 0] K; work_in 1.0, heat_out 1.0 W/m²; q_bottom = −0.5 (down, out of the
    fluid), q_top = +0.5 (up, out) W/m²; dT_max 2.083e-4 K · 3.41e-3 · `-1.0`. *explain:* 1. ε by contraction and by
    sum of squares — identical; 2. the steady Couette solution with viscous heating and its energy bookkeeping; 3. the
    entropy production at mid-gap (no temperature gradient there, so all from ε); 4. with μ < 0 dissipation would be
    negative — the second law rules it out.
31. `nb.check_agree` — **from scratch (curation §7):** ε by the explicit double sum `sig = ch04.viscous_stress(G, 1e-3);
    S = 0.5*(G + G.T); eps = sum(sig[i, j]*S[i, j] for i in range(3) for j in range(3))/1000.0` · `assert
    np.isclose(eps, 1.0)` · and 1000 random velocity gradients (compressible, μ_v = 5×10⁻⁴): the sum-of-squares form equals
    the contraction to 1e-12 and is never negative: `assert all(e >= 0 for e in eps_random)`.
32. `nb.figure` — three panels across the gap for Couette flow (U = 1 m/s) and Couette + adverse pressure gradient (dashed):
    (a) u(y) (teal); (b) ρε(y) (rose strip, uniform for Couette, peaked at the walls for the pressure-driven part);
    (c) T(y) − T₀ (rose) with the energy-budget bars beside it: wall work in (teal), heat out through each wall (rose),
    total (black ≈ 0). *see:* "a straight velocity profile, a flat dissipation strip, a parabolic temperature bump of a
    fifth of a millikelvin." *read:* "the moving wall pumps in 1 W per m²; all of it becomes heat inside and leaves by
    conduction through the two walls; the bump is where the heat has farthest to go." *change:* "…oil instead of water
    (μ 100 × larger, k 4 × smaller): the bump is 400 × larger — 0.08 K — which is why bearings are cooled."
33. `nb.plotly` — `slider_figure` over the wall speed U = 0.1…10 m/s (log, 20 steps, FAST 10): T(y) − T₀ for water and
    oil (two traces) and the peak rise in the title. Title "Viscous heating grows as U²". see/read/change: "…at 10 m/s
    water warms by only 0.02 K; oil by 8 K."
34. `nb.explainer("viscous_dissipation_heating", heading="Where does the energy go when viscosity stops a flow?", why="A
    sheared channel on one clock: the velocity profile, the dissipation ε(y) of
    $\varepsilon=2\nu\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)^2+\frac{\mu_v}\rho S_{mm}^2$ (4.58), and the temperature relaxing to its
    steady parabola, with energy-budget bars (wall work in = heat out) and the entropy production. A 'negative μ'
    preset breaks the second law on purpose.", tries=["Water Couette at 1 m/s: read the peak temperature rise in
    Explain (2.1×10⁻⁴ K).", "Oil bearing preset: the same flow heats 400 times more.", "Press ▶ and watch T(y) grow to
    its steady ghost; the end card says shear work in = heat out.", "Preset 'negative μ': ε and the entropy production turn
    negative and the status flags the second law."])`
35. `nb.md` — **What would change if…** "…the flow were frictionless and heat did not flow? Then ε = 0, q = 0, each particle
    keeps its entropy, and the energy and momentum equations can be integrated along streamlines — the Bernoulli family
    (C11, C12)."

### A.9 §4.9 Special Forms of the Equations — C11 (+N85–N92, N94–N96, N102 · D24, D25), R09, C12 (+N93, N97–N101 · D26 · E6), R10, C13 (+N106–N114 · D27, D28 · E8)
1. `nb.section("4.9", "Special Forms of the Equations", intro="**What is this section about?** Useful consequences of the
   equations under special conditions: the angular-momentum principle, four Bernoulli equations (each with its own
   hypotheses), gravity absorbed into the pressure, and the Boussinesq approximation that ocean and atmosphere models
   use.")`
2. `nb.pointer("The angular-momentum principle, $d\mathbf H/dt=\mathbf M$ (4.64) and its control-volume form (4.65), with the
   lawn sprinkler (Ex. 4.6), was taught with the momentum budgets in C04 above.")`
3. `nb.pointer("The pitot tube, $\lvert\mathbf u\rvert=\sqrt{2(p_2-p_1)/\rho}$, stagnation and dynamic pressure
   $p+\tfrac12\rho\lvert\mathbf u\rvert^2$, and the orifice $u=\sqrt{2gh}$ were taught with Bernoulli's equation
   $\tfrac12U^2+gz+p/\rho=$ const (4.19) in C05 above.")`

**C11 — The Bernoulli function: constant on streamlines and vortex lines (4.71)**
4. `nb.core("C11", "Bernoulli is constant along what, exactly? (4.71)", question="In a bathtub whirlpool the water
   surface dips in the middle, although every particle goes round at a steady speed. If Bernoulli held across the
   circles the surface would be level. Along what, exactly, is Bernoulli's sum constant?")`
5. `nb.md` — **The problem in plain words:** "C05 derived Bernoulli along one streamline of a steady, frictionless,
   constant-density flow. People use it far more widely — across streamlines, in gases, in unsteady flows — sometimes
   correctly, sometimes not. We derive it again from Euler's equation, in a way that shows precisely what it is constant
   along, and why irrotational flows are special."
6. `nb.md` — **The idea:**
   ```
   Euler + Lamb identity:   ∂u/∂t + ∇B = u × ω         B = ½u² + ∫dp/ρ + gz          (4.69)
   steady:                  ∇B = u × ω   ⊥ u and ⊥ ω   →  B constant along streamlines AND vortex lines   (4.71)
   irrotational (ω = 0):    ∇B = 0                     →  B constant everywhere                          (4.72)
   ```
7. `nb.note` — **N85 [C]** "Bernoulli equations are **not new laws**: each is a consequence of Navier–Stokes (4.38) or
   the energy equation (4.60) under stated conditions — so each comes with its own list of hypotheses (the table at the
   end of this block)."
8. `nb.note` — **N86 [B]** "Start from Euler's equation (4.41), $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g$, with
   gravity from its potential Φ = gz (4.18):" equation
   `\frac{\partial u_j}{\partial t}+u_i\frac{\partial u_j}{\partial x_i}=-\frac1\rho\frac{\partial p}{\partial x_j}-\frac{\partial\Phi}{\partial x_j}`, ref "4.66".
9. `nb.note` — **N87 [B]** "**Barotropic flow** — density a function of pressure alone, ρ = ρ(p) (constant density,
   isothermal gas, isentropic gas). Then dp/ρ is an exact differential and" equation
   `\frac1\rho\frac{\partial p}{\partial x_j}=\frac{\partial}{\partial x_j}\int_{p_o}^{p}\frac{dp'}{\rho(p')}`, ref "4.67". "(the
   fundamental theorem of calculus with a variable upper limit, Ch. 2 primer P84). > ⚠️ The prime on p′ here is only the
   integration variable. Closed forms (`ch04.pressure_function`): constant ρ → (p − p_o)/ρ; isothermal → RT ln(p/p_o);
   isentropic → the enthalpy difference h − h_o. A *baroclinic* fluid (ρ depends on T too, like the real atmosphere) has
   no such single function — the root of Ch. 13's thermal wind."
10. `nb.derivation("D24", …)` — Part F D24 (10 steps), ref "4.69".
11. `nb.note` — **N88 [B]** "**The Lamb identity** (D24 steps 3–7; the book leaves it to an exercise):" equation
    `u_i\frac{\partial u_j}{\partial x_i}=-(\mathbf u\times\boldsymbol\omega)_j+\frac{\partial}{\partial x_j}\Big(\tfrac12u_i^2\Big)`, ref
    "4.68". "The advective acceleration = gradient of the kinetic energy minus the 'Lamb vector' u × ω. Ch. 5 starts the
    vorticity equation from here."
12. `nb.note` — **N89 [B]** "**Euler in Bernoulli-function form** (D24's result):" equation
    `\frac{\partial u_j}{\partial t}+\frac{\partial}{\partial x_j}\Big[\tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz\Big]=(\mathbf u\times\boldsymbol\omega)_j`,
    ref "4.69". "The bracket is the **Bernoulli function** B (`ch04.bernoulli_function`)."
13. `nb.derivation("D25", …)` — Part F D25 (5 steps), ref "4.71".
14. `nb.note` — **N90 [B]** "**Steady flow**:" equation `\nabla B=\mathbf u\times\boldsymbol\omega`, ref "4.70". "∇B is normal to
    the surfaces B = const, and u × ω is perpendicular to both u and ω, so each B-surface (a *Lamb surface*) is woven from
    streamlines and vortex lines. Picture (figure below): a Rankine vortex with an axial current — its Lamb surfaces are
    cylinders carrying helical streamlines and axial vortex lines."
15. `nb.note` — **C11's result, restated** — "**Bernoulli 1** (steady, inviscid, barotropic, conservative body force):"
    equation `\tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=\text{constant along streamlines and vortex lines}`, ref "4.71".
16. `nb.note` — **N91 [B]** "**Irrotational** as well (ω = 0 everywhere in a connected region): ∇B = 0 and" equation
    `\tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=\text{constant everywhere}`, ref "4.72". "Number (code below): B
    at 1000 random points of the Ch. 3 cylinder flow agrees to 1e-12; inside a Rankine core B varies from circle to
    circle."
17. `nb.note` — **N92 [C]** "An inviscid, barotropic flow in an inertial frame that starts irrotational stays
    irrotational (proved with Kelvin's circulation theorem in Ch. 5). That is why the flow outside a thin boundary layer
    round a body is irrotational and (4.72) applies there (Ch. 6, Ch. 9)."
18. `nb.note` — **N94–N95 [B]** "**Energy form.** Start instead from the total-energy equation (4.55) with σ = q = 0 and
    steady flow:" equation `\rho u_i\frac{\partial}{\partial x_i}\Big(e+\tfrac12u_j^2\Big)=\rho u_ig_i-\frac{\partial}{\partial x_j}\Big(\rho u_j\frac p\rho\Big)`,
    ref "4.76". "With g = −∇(gz) and steady continuity ∂(ρu_i)/∂x_i = 0 it becomes" equation
    `\rho u_i\frac{\partial}{\partial x_i}\Big(e+\frac p\rho+\tfrac12u_j^2+gz\Big)=0`, ref "4.77".
19. `nb.note` — **N96 [B]** "**Bernoulli 3** (steady, inviscid, non-conducting; h = e + p/ρ):" equation
    `h+\tfrac12\lvert\mathbf u\rvert^2+gz=\text{constant on streamlines}`, ref "4.78". "It is an *energy* statement: it holds in a
    gas where kinetic and thermal energy trade places. For isentropic flow dp/ρ = dh (Gibbs), so (4.71) and (4.78) agree.
    Number: air at 300 K moving at 100 m/s brought to rest heats to $T_0=T+U^2/2C_p=300+10^4/2009=304.98$ K (Ch. 15)."
20. `nb.note` — **N102 [B]** "**Which Bernoulli?** — a decision table (`ch04.BERNOULLI_FORMS`, `ch04.which_bernoulli`):"
    Table | form | steady? | inviscid? | irrotational? | barotropic / ρ const? | constant along |: (4.19)
    $\tfrac12U^2+gz+p/\rho$ — steady, inviscid, not needed, ρ const — one streamline · (4.71)
    $\tfrac12u^2+\int dp/\rho+gz$ — steady, inviscid, not needed, barotropic — streamlines and vortex lines · (4.72)
    same B — steady, inviscid, yes, barotropic — everywhere · (4.75)
    $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int dp/\rho+gz$ — **unsteady allowed**, inviscid, yes,
    barotropic — everywhere (a function of t absorbed into φ) · (4.78) $h+\tfrac12\lvert\mathbf u\rvert^2+gz$ — steady, inviscid
    and non-conducting, not needed, any — streamlines · (4.82)/(4.83) — unsteady, **viscous allowed** (ρ, μ const),
    yes — along a streamline at one instant / everywhere.
21. `nb.worked_example("a Rankine vortex in water (Γ = 2π m²/s, core radius σ = 1 m, p∞ = 0)", "1. Inside the core
    (r < σ) the fluid turns like a solid body: $u_\theta=\Gamma r/(2\pi\sigma^2)=r$ (s⁻¹ × m), vorticity ω = Γ/(πσ²) = 2 s⁻¹.
    Outside it is irrotational: u_θ = Γ/(2πr) = 1/r. 2. Pressure from the radial balance dp/dr = ρu_θ²/r: outside
    p/ρ = −1/(2r²); inside p/ρ = −1 + r²/2 (matched at r = 1). 3. B = ½u_θ² + p/ρ: outside ½/r² − ½/r² = **0** everywhere;
    inside r²/2 − 1 + r²/2 = **r² − 1**: −1 at the centre, −0.75 at r = 0.5, 0 at the edge. 4. Along a circle (a
    streamline) B is constant — (4.71) holds; across circles it changes inside the rotational core and not outside —
    (4.72) holds only where ω = 0. 5. Check D25: dB/dr = 2r = 1 at r = 0.5, and |u × ω| = u_θω = 0.5 × 2 = 1 ✓.")`
22. `nb.code` — *code:* `r = np.array([0.0, 0.5, 1.0, 2.0, 4.0])` · `print(ch04.rankine_bernoulli(r, 2*np.pi, 1.0)["B"])`
    · Lamb identity at a point of the Rankine vortex (Cartesian field from `core.vortices`, via `ch03`): `u_rk =
    ch03.vortex_velocity_field("rankine", Gamma=2*np.pi, sigma=1.0)`; `print(ch04.lamb_vector(u_rk, np.array([0.5, 0.0,
    0.0]))[:2], ch04.lamb_identity_terms(u_rk, np.array([0.5, 0.0, 0.0]))["residual"])` · the cylinder: `pts =
    np.random.default_rng(0).uniform(-3, 3, size=(2, 1000)); pts = pts[:, np.hypot(*pts) > 1.05]`; `B = ch04.bernoulli_along_line(
    ch03.cylinder_velocity_field(1.0, 1.0), lambda x, t: ch04.exact_solution("cylinder", x, U=1.0, a=1.0, rho=1.0)[1],
    pts, rho=1.0, g=0.0)`; `print(B.max() - B.min())` · `print(ch04.which_bernoulli(steady=True, viscous=False,
    irrotational=False, barotropic=True))` · `print(ch04.pressure_function(2e5, 1e5, "isothermal", T=300.0))` ·
    `print(ch04.stagnation_temperature(300.0, 100.0))`. *expect:* `[-1. -0.75 0. 0. 0.]` · Lamb vector (1.0, 0.0) (radial,
    outward), residual ≈ 0 · spread ≤ 1e-10 · `['4.71']` · 59 692 J/kg (= 287.058 × 300 × ln 2) · `304.98`. *explain:*
    1. B across the vortex — varies inside, constant outside; 2. the Lamb identity at a point inside the core (u × ω
    points outward, equal to dB/dr); 3. B at a thousand points of an irrotational flow — one constant; 4. the decision
    function for a steady, inviscid, barotropic rotational flow; 5. an isothermal gas's pressure function; 6. the
    energy form's stagnation temperature.
23. `nb.check_agree` — **from scratch (curation §7):** u × ω by `np.cross` with ω from central differences
    (`kinematics.vorticity`) and the Lamb identity checked pointwise: `adv = kinematics.acceleration(u_rk, x0, 0.0).advective`;
    `grad_ke` by central differences of ½|u|²; `assert np.allclose(adv, -np.cross(u, w)[:2] + grad_ke, atol=1e-7)` and
    against `ch04.lamb_vector`.
24. `nb.figure` — two panels for the Rankine vortex: (a) B(r) across circles (rose) with the core shaded, flat zero
    outside; (b) B along one circle inside the core (r = 0.5, teal, flat at −0.75) and along one outside (r = 2, teal, flat
    at 0) vs angle. Beside them a small rendering of the N102 decision table with the (4.71) and (4.72) rows lit. *see:*
    "(a) a parabola inside the core joining a flat line outside; (b) two flat lines." *read:* "along any streamline B is
    constant (4.71); across streamlines it changes only where there is vorticity — outside the core the flow is
    irrotational and one B holds everywhere (4.72). The whirlpool's dip is exactly the core's B deficit." *change:* "…a
    uniform flow instead: B would be one constant everywhere; the whirlpool needs vorticity for its dip."
25. `nb.plotly` — `slider_figure` over the core radius σ = 0.2…2 m (19 steps, FAST 10) of the Rankine vortex (Γ = 2π):
    traces B(r) (rose) and u_θ(r) (teal) on r ∈ [0, 4] m. Title "Rotational core: B varies across streamlines only where ω
    ≠ 0". see/read/change: "…shrinking σ deepens the dip (B(0) = −Γ²/(4π²σ²)) and narrows it."
26. `nb.md` — **What would change if…** "…the flow were irrotational but unsteady — water sloshing in a U-tube, a
    balloon pushed suddenly through water? Then the ∂u/∂t term of (4.69) stays, but it is a gradient too, and a new
    Bernoulli equation holds everywhere at each instant (C12)."

**C12 — Unsteady Bernoulli for potential flow (4.75)**
27. `nb.recap("R09", "Velocity potential", "An irrotational flow (ω = ∇×u = 0) in a simply connected region can be written
    as the gradient of a scalar, $\mathbf u\equiv\nabla\phi$ *(Eq. 4.73, first met as (3.17) in Ch. 3)*, because the curl of a
    gradient is zero. φ is defined up to an added constant — or an added function of time only, which changes no
    velocity.", where="Ch. 3 §3.4")`
28. `nb.core("C12", "Unsteady Bernoulli: the pressure of an accelerating flow (4.75)", question="Water sloshes back and
    forth in a U-tube. At the instant the column is momentarily still, the pressures at the two ends still differ. Why —
    and by how much?")`
29. `nb.md` — **The problem in plain words:** "Steady Bernoulli says pressure is high where flow is slow. In an
    accelerating flow that is not enough: pushing fluid to speed it up takes a pressure difference even at zero speed.
    Water waves (Ch. 7), a body starting to move (its 'added mass', Ch. 6), a bubble collapsing and sound (Ch. 15) all
    need the unsteady form — valid everywhere, but only for irrotational flow."
30. `nb.md` — **The idea:**
   ```
   u = ∇φ  →  Lamb term vanishes, ∂u/∂t = ∇(∂φ/∂t)  →  ∇[∂φ/∂t + ½|∇φ|² + ∫dp/ρ + gz] = 0
          →  the bracket is the same everywhere: a function of time B(t) only  →  absorb it into φ  →  (4.75)
   U-tube:  p₁ − p₂ = ρ L dU/dt   (the ∂φ/∂t term, amber, carries the difference when U = 0)
   ```
31. `nb.derivation("D26", …)` — Part F D26 (8 steps), ref "4.75" (the Schwarz primer of C08 is recalled in step 2's
    *why*).
32. `nb.note` — **N93 [B]** "D26 steps 3–5 are (4.74):" equation
    `\nabla\Big[\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz\Big]=0,\quad\text{or}\quad\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=B(t)`,
    ref "4.74". "> ⚠️ **Book slip (sign):** to absorb B(t) the potential must be redefined as
    $\phi\to\phi-\int_{t_o}^{t}B(t')dt'$; the printed '+' makes the bracket 2B(t) instead of 0 (step 6 and the code below).
    The result (4.75) is unaffected; the velocity ∇φ does not change either way."
33. `nb.note` — **C12's result, restated** — "**Bernoulli 2** (unsteady, inviscid, irrotational, barotropic):" equation
    `\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=\text{constant}`, ref "4.75".
    "Ch. 7 applies it at the free surface of every water wave."
34. `nb.note` — **N97 [B]** "**Viscosity too?** Put the Lamb identity (4.68) and the viscous force $-\mu\nabla\times\boldsymbol\omega$
    (4.40) into (4.39b):" equation
    `\rho\frac{\partial\mathbf u}{\partial t}+\rho\nabla\big(\tfrac12\lvert\mathbf u\rvert^2\big)-\rho\mathbf u\times\boldsymbol\omega=-\nabla p+\rho\mathbf g-\mu\nabla\times\boldsymbol\omega`,
    ref "4.79".
35. `nb.note` — **N98 [B]** "With ω = 0 and constant ρ **both** the Lamb term and the viscous term vanish — viscosity is
    present but exerts no net force on an irrotational flow:" equation
    `\rho\frac{\partial\mathbf u}{\partial t}+\nabla\big(\tfrac12\rho\lvert\mathbf u\rvert^2+\rho gz+p\big)=0`, ref "4.80". "Code:
    `ch04.viscous_irrotational_residual(ch03.cylinder_velocity_field(1.0, 1.0), np.array([1.5, 0.7]))` → 0 to round-off.
    (Viscosity still *dissipates* energy in such a flow — ε ≠ 0 — which damps water waves slowly, Ch. 7.)"
36. `nb.note` — **N99 [B]** "Dot (4.80) with a line element ds along a streamline and integrate from point 1 to point 2:"
    equation `\int_1^2\frac{\partial\mathbf u}{\partial t}\cdot d\mathbf s+\int_1^2\frac{\partial}{\partial s}\Big(\tfrac12\lvert\mathbf u\rvert^2+gz+\frac p\rho\Big)ds=0`,
    ref "4.81". "(the directional derivative along s, Ch. 2 P75)"
37. `nb.note` — **N100 [B]** "**Bernoulli 4** (constant ρ and μ, irrotational, along a streamline at one instant):" equation
    `\int_1^2\frac{\partial\mathbf u}{\partial t}\cdot d\mathbf s+\Big(\tfrac12\lvert\mathbf u\rvert^2+gz+\frac p\rho\Big)_2=\Big(\tfrac12\lvert\mathbf u\rvert^2+gz+\frac p\rho\Big)_1`,
    ref "4.82". "Number: a water column L = 1 m long accelerating at dU/dt = 1 m/s² along a straight tube needs
    p₁ − p₂ = ρL dU/dt = 1000 Pa, even at the instant U = 0."
38. `nb.note` — **N101 [B]** "In potential form (u = ∇φ, swap ∂/∂t and ∇):" equation
    `\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+gz+\frac p\rho=\text{constant}`, ref "4.83". "— the constant-ρ
    (4.75), now shown to hold with constant viscosity as well."
39. `nb.worked_example("the U-tube and an accelerating sphere", "1. **U-tube:** a water column of length L = 1 m with its
    two surfaces displaced ±h oscillates at ω = √(2g/L) = 4.43 rad/s (period 1.42 s). At the instant the column is still
    (U = 0) its acceleration is largest: with h₀ = 5 cm, dU/dt = ω²h₀ = 19.62 × 0.05 = 0.98 m/s², and (4.82) gives
    p₁ − p₂ = ρL dU/dt = 1000 × 1 × 0.98 = 981 Pa — exactly the weight of the 2h₀ = 10 cm height difference, ρg(2h₀) =
    981 Pa ✓. 2. **Sphere** of radius a = 0.1 m accelerating from rest at 1 m/s² through still water: the unsteady term
    gives a surface pressure ½ρa cos θ dU/dt (high in front, low behind) whose net force is −½ρV dU/dt, V = 4πa³/3 =
    4.19×10⁻³ m³: the sphere must push as if it carried half its volume of water with it — 'added mass' 2.09 kg,
    force 2.09 N (Ch. 6).")`
40. `nb.code` — *code:* `c = ch04.u_tube_column(np.array([0.0, 0.3547, 0.7093]), L=1.0, h0=0.05)`; `print(c["h"], c["U"],
    c["dUdt"], c["dp"])` · `print(ch04.unsteady_streamline_bernoulli(np.full(11, 1.0), np.linspace(0, 1, 11), {"U": 0.0,
    "z": 0.0, "p": 1000.0}, {"U": 0.0, "z": 0.0, "p": 0.0}))` · `s = ch04.accelerating_sphere_pressure(np.array([0.0,
    np.pi/2, np.pi]), a=0.1, dUdt=1.0)`; `print(s["p"], s["force"], s["added_mass"])` · `print(ch04.gauge_absorbed_bracket(3.0,
    sign=-1.0), ch04.gauge_absorbed_bracket(3.0, sign=+1.0))` · `print(ch04.viscous_irrotational_residual(
    ch03.cylinder_velocity_field(1.0, 1.0), np.array([1.5, 0.7, 0.0]), mu=1e-3))`. *expect:* h = [0.05, ~0, −0.05], U =
    [0, −0.2215, ~0], dU/dt = [−0.981, ~0, +0.981] (sign conventions of the function's docstring), dp = [−981, ~0, +981] Pa
    · residual `0.0` (1000 Pa is what the acceleration needs) · p = [+50, 0, −50] Pa, force −2.094 N, added mass 2.094 kg
    · `0.0 6.0` (the printed sign doubles B = 3) · ≈ (0, 0, 0). *explain:* 1. the column's level, speed, acceleration
    and the pressure difference that drives it; 2. the streamline form (4.82) balanced by ρL dU/dt; 3. the unsteady
    surface pressure on the sphere and the added-mass force; 4. the gauge sign test of (4.74); 5. viscosity exerts no net
    force on an irrotational flow.
41. `nb.check_agree` — **from scratch:** the gauge move in sympy: `t = sp.symbols('t'); phi = sp.Function('phi')(t); Bt =
    sp.Function('B')(t)`; bracket B(t) = φ′ + rest; with `phi_new = phi - sp.Integral(Bt, t)` let `rest =
    Bt - sp.diff(phi, t)` stand for everything in the bracket except ∂φ/∂t (so that ∂φ/∂t + rest = B), then the new
    bracket `new = sp.diff(phi_new, t) + rest` is `0` for the minus sign and `2*B(t)` with `phi + sp.Integral(Bt, t)`; `assert sp.simplify(new_minus) == 0` and `assert sp.simplify(new_plus -
    2*Bt) == 0`; then `assert np.isclose(ch04.gauge_absorbed_bracket(3.0, -1.0), 0.0)`.
42. `nb.animation` — (`player="video"`, 60 frames, FAST 30) a U-tube (two vertical legs joined at the bottom), the column
    oscillating h(t) = h₀ cos ωt, with a bar beside it for p₁ − p₂ (amber, the ∂φ/∂t part) and a small readout "U = …,
    dU/dt = …"; a faint ghost of the steady-Bernoulli prediction (which would say p₁ − p₂ = 0 whenever U = 0).
    see/read/change: *see:* "the amber bar is largest when the column stops and zero when it moves fastest." *read:* "the
    pressure difference tracks the acceleration, not the speed — the ∂φ/∂t term of (4.75), invisible in any steady
    Bernoulli." *change:* "…a column twice as long: ω falls by √2, and for the same h₀ the bar keeps its size (ρL dU/dt =
    2ρgh₀) — the weight of the height difference."
43. `nb.figure` — surface pressure on the accelerating sphere vs θ (amber, ½ρa cos θ dU/dt) with the steady part
    (grey, 0 at the start of the motion) and a polar inset with pressure arrows. *see:* "a cosine: pushed in front,
    pulled behind." *read:* "summed over the surface it resists the acceleration like extra mass: half the displaced
    water." *change:* "…the sphere moving at steady speed: the amber curve vanishes and the steady Bernoulli pressure
    1 − (9/4)sin²θ times ½ρU² remains, symmetric front to back — no net force (d'Alembert's paradox, Ch. 6)."
44. `nb.explainer("which_bernoulli", heading="Bernoulli is constant along what, exactly?", why="Six scenarios — pitot
    tube, tank orifice, Rankine vortex, cylinder, U-tube, hot-gas nozzle — each with its streamlines, two probes, and B
    plotted along and across streamlines as stacked bars (½u², p/ρ, gz, ∂φ/∂t). Toggle the hypotheses and the status picks
    the valid equation: (4.19), (4.71), (4.72), (4.75) or (4.78).", tries=["Rankine vortex: B is flat along each circle
    but not across the core — the status says '(4.71) along streamlines only'.", "Cylinder: B is flat everywhere
    ((4.72)); drag the probes anywhere.", "U-tube: stop the column (U = 0) — the amber ∂φ/∂t bar carries the pressure
    difference ((4.75)).", "Hot nozzle: B of (4.71) with ρ = const fails; the energy form $h+\tfrac12\lvert\mathbf u\rvert^2+gz$
    (4.78) stays flat."])`
45. `nb.md` — **What would change if…** "…the density varied a little from place to place, as in a lake warmed at the top?
    Then gravity no longer cancels into the pressure, and a small density difference, multiplied by the large g, drives
    motion — the Boussinesq approximation (C13)."

**C13 — The Boussinesq approximation (4.86), (4.89)**
46. `nb.recap("R10", "Hydrostatic balance", "A fluid at rest in gravity has $0=-\nabla p_s+\rho_s\mathbf g$ — pressure
    increases downward at the rate ρ_s g (Ch. 1, $dp/dz=-\rho g$ (1.8)). ρ_s(z) and p_s(z) may vary with height (a
    stratified reference state).", where="Ch. 1 §1.7")` — followed by `nb.code`: `z = np.linspace(0, -100, 5)`;
    `print(statics.integrate_hydrostatic(z, lambda z, p: 1025.0, 1.013e5))` (`statics` is in the setup imports).
    *expect:* ≈ [101300, 352681, 604063, 855444, 1106825] Pa (ρg = 10055.25 Pa/m: 1 atm every ≈ 10 m).
47. `nb.core("C13", "Boussinesq: density matters only where it meets gravity (4.86)", question="A lake is 2 °C warmer at
    the top than at the bottom — its density changes by less than one part in a thousand. Yet that tiny difference stops
    the wind from mixing the lake. Why does such a small density change matter so much, and when may we ignore it?")`
48. `nb.md` — **The problem in plain words:** "In the ocean and the atmosphere, density varies by a percent or less, but
    those variations drive the overturning circulation, sea breezes, convection and internal waves. Solving the full
    compressible equations for them is wasteful (and fills the solution with sound waves). Boussinesq's idea keeps the
    density variation in exactly one place — where it is multiplied by g — and treats the fluid as incompressible
    everywhere else. It is the equation set of Ch. 7's internal waves, Ch. 11's convection and Ch. 13's ocean and
    atmosphere."
49. `nb.md` — **The idea:**
   ```
   ρ = ρ₀ + ρ'   with  ρ'/ρ₀ ~ αδT ~ 10⁻³
   inertia   ρ Du/Dt  ≈ ρ₀ Du/Dt           (error 10⁻³: drop ρ')
   gravity   ρ g      = ρ_s g + ρ' g         (ρ_s g balanced by the hydrostatic p_s; ρ' g KEPT: g ≫ Du/Dt)
   continuity  ∇·u ≈ 0;   heat  DT/Dt = κ∇²T  (pressure work turns C_v into C_p)
   ```
50. `nb.primer("order-of-magnitude scaling", "To compare terms of an equation without solving it, replace each quantity by
    its typical size: a velocity varying by U over a distance L has ∂u/∂x ~ U/L and ∂²u/∂x² ~ U/L²; a temperature varying
    by δT has ∂T/∂x ~ δT/L. Ratios of such estimates tell which terms can be dropped (≪ 1) — the method of C13 and C15.",
    code="U, L, nu = 0.1, 1.0, 1e-6                     # a lake current, its scale, water's ν\nadvective = U**2/L; viscous = nu*U/L**2          # (u·∇)u ~ U²/L,  ν∇²u ~ νU/L²\nprint(advective, viscous, advective/viscous)      # 0.01, 1e-07, 1e5 = Re")`
    — *expect:* `0.01 1e-07 100000.0`.
51. `nb.primer("reduced gravity and buoyancy", "A parcel lighter than its surroundings by Δρ feels a net upward force per
    unit mass g′ = gΔρ/ρ₀ — gravity *reduced* by the density contrast. The buoyancy b = −gρ′/ρ₀ (upward positive) is the
    same idea as a field. For a density contrast of 10⁻³, g′ ≈ 0.01 m/s²: small, but it acts on the whole mass and for as
    long as the contrast lasts. N² = ∂b/∂z is the buoyancy frequency of Ch. 1.", code="rho0, d_rho = 1000.0, -0.4
    # a parcel 0.4 kg/m³ lighter\nprint(ch04.buoyancy(d_rho, rho0), 9.81*0.4/1000)   # +0.003924 m/s² upward, both ways")`
    — *expect:* `0.003924 0.003924`.
52. `nb.note` — **N106 [B]** "Subtract the hydrostatic state (R10) from (4.39b); with p′ = p − p_s and ρ′ = ρ − ρ_s:"
    equation `\rho\frac{D\mathbf u}{Dt}=-\nabla p'+\rho'\mathbf g+\mu\nabla^2\mathbf u`, ref "4.84". "> ⚠️ **Primes again**: here p′
    and ρ′ are departures from the hydrostatic state — not the rotating frame of C09, not the dummy variable of (4.67)."
53. `nb.note` — **N107 [B]** "**Constant density**: ρ′ = 0 and gravity disappears from the equation of motion:" equation
    `\rho\frac{D\mathbf u}{Dt}=-\nabla p'+\mu\nabla^2\mathbf u`, ref "4.85". "Gravity returns as soon as there is a free
    surface, an interface or any density variation. Demo: the Poiseuille velocity in a vertical-gravity channel solves
    (4.39b) with p including −ρgy and (4.85) with p′ alone (`ch04.ns_incompressible_terms` with and without g)."
54. `nb.note` — **N108 [B]** "**When Boussinesq holds** (stated; the scaling chain is the book's): low Mach number and no
    sound; a vertical scale L ≪ c²/g (the height over which hydrostatic compression changes the density — ≈ 11.8 km for
    air at 288 K); and small temperature differences, αδT ≪ 1, because then" equation
    `\frac{(1/\rho)(D\rho/Dt)}{\nabla\cdot\mathbf u}\sim\frac{(U/\rho)(\delta\rho/L)}{U/L}=\frac{\delta\rho}{\rho}=\alpha\,\delta T\ll1`,
    no ref. "so continuity reduces to ∇·u = 0. Table (from `ch04.boussinesq_validity` on the five scenarios): lake (water,
    δT = 5 K): 1×10⁻³ · ocean thermocline (δT = 10 K): 2×10⁻³ · lab tank: 4×10⁻⁴ · atmospheric boundary layer (air, 300 K,
    δT = 10 K, L = 1 km): αδT = 0.033, L/(c²/g) ≈ 0.08 · deep atmosphere (L = 10 km): L/(c²/g) ≈ 0.8 — **fails**."
55. `nb.derivation("D27", …)` — Part F D27 (9 steps), ref "4.86".
56. `nb.note` — **C13's momentum equation, restated** — equation
    `\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u`, ref "4.86". "ρ₀ is a constant
    reference density. In buoyant convection the term ρ′g/ρ₀ is as large as ∂w/∂t or ν∇²w."
57. `nb.note` — **N109 [B]** "The energy equation (4.60) in vector form, the start of D28:" equation
    `\rho\frac{De}{Dt}=-p\nabla\cdot\mathbf u+\rho\varepsilon-\nabla\cdot\mathbf q`, ref "4.87".
58. `nb.derivation("D28", …)` — Part F D28 (8 steps), ref "4.89".
59. `nb.note` — **N110 [B]** "**Why C_p and not C_v** (D28 steps 2–5): although ∇·u ≈ 0 in continuity, the pressure work
    p∇·u is *not* negligible in the energy equation, because p is large: for a perfect gas
    $-p\nabla\cdot\mathbf u\cong-p\alpha\frac{DT}{Dt}=-\rho(C_p-C_v)\frac{DT}{Dt}$. It moves to the left and turns ρC_v DT/Dt into
    ρC_p DT/Dt. > ⚠️ For liquids the argument is different: p∇·u is small, but C_p ≈ C_v anyway, so the ocean uses the same
    equation with water's C_p."
60. `nb.note` — **N111 [B]** equation `\rho C_p\frac{DT}{Dt}=\rho\varepsilon-\nabla\cdot\mathbf q`, ref "4.88".
61. `nb.note` — **N112 [B]** "**Viscous heating is negligible** under Boussinesq conditions:" equation
    `\frac{\rho\varepsilon}{\rho C_p(DT/Dt)}\sim\frac{2\mu S_{ij}S_{ij}}{\rho C_pu_i(\partial T/\partial x_i)}\sim\frac{\mu U^2/L^2}{\rho C_pU(\delta T/L)}=\frac{\nu U}{C_p\,\delta T\,L}`,
    no ref. "Water, U = 0.1 m/s, L = 1 m, δT = 1 K: 10⁻⁷/4186 ≈ 2×10⁻¹¹ (the book's 'typical' value is 10⁻⁷ — tiny either
    way). This is the low-Eckert-number condition of C15."
62. `nb.note` — **N113 [B]** "**The Boussinesq heat equation** (D28's result), κ = k/ρC_p:" equation
    `\frac{DT}{Dt}=\kappa\nabla^2T`, ref "4.89". "Exact test solution (the animation below): a warm Gaussian blob carried by a
    uniform current U and spreading, $T'=T_a(\sigma_0/\sigma)^2\exp(-\lvert\mathbf x-\mathbf Ut\rvert^2/2\sigma^2)$ in 2-D with
    $\sigma^2=\sigma_0^2+2\kappa t$."
63. `nb.note` — **N114 [B]** "**The Boussinesq set** (boxed): continuity ∇·u = 0 (4.10); momentum (4.86) with g = −g e_z;
    heat (4.89); and a linear equation of state" equation `\rho=\rho_0\big[1-\alpha(T-T_0)\big]`, no ref. "with ρ′ =
    −ρ₀α(T − T₀), so the buoyancy is b = −gρ′/ρ₀ = gα(T − T₀) and the stratification N² = ∂b/∂z = gα dT/dz (Kundu's
    Γ ≡ dT/dz; in the meteorological convention Γ_met = −dT/dz, N² = −gαΓ_met). For seawater, salinity enters too
    (`ch01.seawater_density_linear`). Omitting the Coriolis term — Ch. 13 adds it back."
64. `nb.worked_example("three Boussinesq numbers", "1. **Lake** (water, α ≈ 2×10⁻⁴ K⁻¹, δT = 10 K): αδT = 2×10⁻³ ≪ 1;
    reduced gravity g′ = gαδT = 9.81 × 2×10⁻³ = 0.0196 m/s². 2. **Air** at 300 K (perfect gas, α = 1/T): δT = 10 K →
    αδT = 0.033 — still small. 3. **Depth limit**: c²/g for air at 288 K = 340.3²/9.81 = 11.8 km — a flow 1 km deep is
    fine, the whole troposphere is not. 4. **Rise speed** of a 1 m warm patch in the lake: √(g′L) ≈ 0.14 m/s. 5. **Viscous
    heating** in the lake: νU/(C_pδT L) = 10⁻⁶ × 0.1/(4186 × 1 × 1) ≈ 2×10⁻¹¹.")`
65. `nb.code` — *code:* `for name in ("lake", "thermocline", "lab_tank", "abl", "deep_atmosphere"): p =
    ch04.boussinesq_scenario(name)`; `v = ch04.boussinesq_validity(p["alpha"], p["dT"], p["L"], p["U"], c=p["c"],
    nu=p["nu"], cp=p["cp"])`; `print(f"{name:16s} αδT={v['alpha_dT']:.1e}  L/(c²/g)={v['L_over_Hc']:.2e}
    heating={v['heating_ratio']:.1e}  g'={v['g_prime']:.3g}  → {v['verdict']}")` · `print(ch04.boussinesq_density(303.15,
    1000.0, 2e-4, 293.15), ch04.buoyancy(-2.0, 1000.0))` · blob residual: `Tb = lambda x, t:
    ch04.gaussian_blob_advection_diffusion(x, t, U=np.array([0.01, 0.0]), kappa=1e-4, sigma0=0.05)`; `ub = lambda x, t:
    np.array([0.01, 0.0])`; `print(ch04.temperature_equation_residual(Tb, ub, np.array([0.02, 0.01]), 10.0, kappa=1e-4,
    h=1e-3, ht=1e-2))`. *expect:* five lines, the last with verdict "⚠️ deep layer: L ≈ c²/g" and the others "✅
    Boussinesq valid" · `998.0 0.0196` (m/s², upward for a lighter parcel ρ′ = −2 kg/m³) · residual ≈ 0 (≤ 1e-6 of the
    terms). *explain:* 1. the three small numbers for five real situations; 2. the linear equation of state and the
    buoyancy of a 2 kg/m³ lighter parcel; 3. the blob solves (4.89) exactly (κ = 10⁻⁴ m²/s: an eddy diffusivity, to see it
    spread in seconds).
66. `nb.check_agree` — **from scratch (curation §7):** the heat-equation residual of the blob by hand (central
    differences for ∂T/∂t, u·∇T and the Laplacian) `assert abs(dTdt + adv - 1e-4*lap) < 1e-6*abs(dTdt)` and against
    `ch04.heat_equation_terms(Tb, ub, x0, 10.0, rho=1000.0, cp=4186.0, k=1e-4*1000*4186, eps=0.0)["residual"]`.
67. `nb.animation` — (`player="video"`, 60 frames, FAST 30) a warm Gaussian blob in a 2-D tank (T′ heat map, rose)
    carried to the right by a uniform current and spreading, σ² = σ₀² + 2κt shown as a growing dashed circle, with the
    peak value (∝ σ₀²/σ²) in the title. see/read/change: *see:* "the blob drifts and flattens; its total heat (area under
    the bump) stays constant." *read:* "advection moves it (u·∇T), diffusion spreads it (κ∇²T); the peak falls as
    1/(1 + 2κt/σ₀²)." *change:* "…κ halved: the circle grows half as fast; the drift is unchanged."
68. `nb.figure` — the validity map: three horizontal log-scale bars per scenario (αδT, L/(c²/g), νU/(C_pδT L)) with the
    threshold 0.1 as a dashed line, scenarios as rows; the deep atmosphere's L/(c²/g) bar crossing the line in amber.
    *see:* "most bars far to the left of the line; one crossing." *read:* "Boussinesq is excellent for lakes, oceans and
    boundary layers; it fails for flows as deep as the atmosphere's scale height (then use the anelastic or full
    compressible equations)." *change:* "…a 50 K temperature difference in air: αδT ≈ 0.17 — marginal."
69. `nb.explainer("boussinesq_buoyancy", heading="Density hardly changes — so why does it drive the flow?", why="Choose a
    fluid, a temperature contrast and a depth: three small numbers are compared with their thresholds while a warm blob
    rises and spreads under exactly the terms Boussinesq keeps,
    $\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u$ (4.86) and
    $\frac{DT}{Dt}=\kappa\nabla^2T$ (4.89). Switch the buoyancy term off and the blob stops rising, though ρ changed by only
    0.2 %.", tries=["Lake preset: read αδT and g′ in Explain.", "Toggle 'keep buoyancy' off: the blob only diffuses.",
    "Deep-atmosphere preset: the status turns amber — L is comparable to c²/g.", "Derivation tab D27, the 'why keep ρ′g'
    step: the grey dropped bar is a thousand times smaller than the blue buoyancy bar."])`
70. `nb.md` — **What would change if…** "…we asked what happens at the edges — the walls, the free surface, the interface
    between two fluids? The field equations are complete; to solve them we need conditions on the boundaries (C14)."

### A.10 §4.10 Boundary Conditions — R11, C14 (+N115–N129 · D29)
1. `nb.section("4.10", "Boundary Conditions", intro="**What is this section about?** The field equations need conditions
   at the edges of the flow. Conservation laws applied to a thin 'pillbox' straddling an interface give continuity of
   mass flux, traction and heat flux; experience adds no-slip and no temperature jump; a moving surface made of fluid
   obeys the kinematic condition Dη/Dt = 0; and surface tension adds a pressure jump across a curved interface.")`
2. `nb.recap("R11", "Surface tension", "Molecules at a liquid surface are pulled inward more than outward; the surface
   behaves like a stretched membrane with tension σ [N/m], the energy per unit area of surface. Water at 20 °C has
   σ ≈ 0.0728 N/m (Ch. 1 §1.6, `ch01.surface_tension_water`). Here we *derive* the pressure jump it causes.", where="Ch. 1
   §1.6")`

**C14 — The kinematic boundary condition (4.91), with the wall and interface conditions**
3. `nb.core("C14", "The surface is made of fluid: the kinematic condition (4.91) and other edge rules", question="A wave
   runs across the sea. The water on the surface stays on the surface — it neither sinks into the sea nor flies off.
   What does that sentence say as an equation?")`
4. `nb.md` — **The problem in plain words:** "Every solution of Navier–Stokes is fixed by what happens at its edges: the
   wall of a pipe, the sea surface, the boundary between warm and cold layers of the ocean, the skin of a raindrop.
   Some conditions follow from the conservation laws (what crosses an interface must be continuous), some from
   experiment (fluid sticks to walls), and one from geometry (a free surface moves with the fluid on it). Ch. 7 builds
   water waves on the last one; Ch. 13's layered ocean models on the same condition at each layer interface."
5. `nb.md` — **The idea:**
   ```
   pillbox across an interface, height l → 0:  volume and side terms vanish  →  flux in = flux out:
         ρu·n continuous, traction n_iτ_ij continuous, k ∂T/∂n continuous
   experiment at a solid wall:  no slip (u·t = 0 relative to the wall), no temperature jump (T₁ = T₂)
   a surface η(x, t) = 0 made of fluid:  Dη/Dt = ∂η/∂t + u·∇η = 0   (4.91)  — normal velocities match
   curved interface with tension σ:  Δp = σ(1/R₁ + 1/R₂)
   ```
6. `nb.note` — **N115 [C]** "**What must be specified**: the velocity on every bounding surface; for an external flow, the
   velocity and thermodynamic state on a distant closed surface. The solvers of Ch. 8–10 show it in action."
7. `nb.note` — **N116 [B]** "**Interface conditions from a pillbox** (stated; the book's Fig. 4.18 argument): apply (4.5),
   (4.17) and (4.48) to a flat cylinder straddling the interface and let its height l → 0 — the volume integrals and the
   side-wall fluxes (∝ l) vanish, leaving the two end faces:" equation
   `\rho_1\mathbf u_1\cdot\mathbf n=\rho_2\mathbf u_2\cdot\mathbf n,\qquad n_i\tau^{(1)}_{ij}=n_i\tau^{(2)}_{ij},\qquad k_1\frac{\partial T_1}{\partial n}=k_2\frac{\partial T_2}{\partial n}`,
   no ref. "(interface at rest; traction continuous if surface tension is neglected — with it, the Laplace jump below
   is added). **Two worked cases** (code): (a) a composite wall, k₁ = 1, k₂ = 0.1 W/(m K), each 0.1 m thick, 30 °C and
   20 °C outside: the heat flux 10/(0.1/1 + 0.1/0.1) = 9.09 W/m² is the same in both layers, the interface is at
   29.09 °C and the temperature slopes differ by the factor k₁/k₂ = 10 — the profile has a kink; (b) two-fluid Couette,
   water (μ₁ = 10⁻³ Pa s, 1 mm) under oil (μ₂ = 0.1 Pa s, 1 mm), top plate at 1 m/s: the shear stress
   τ = U/(h₁/μ₁ + h₂/μ₂) = 0.990 Pa is continuous, the velocity slopes differ by μ₂/μ₁ = 100, the interface moves at
   0.990 m/s." Followed by `nb.code`: `print(ch04.two_layer_conduction(np.array([0.05, 0.15]), 1.0, 0.1, 0.1, 0.1, 303.15,
   293.15))` · `print(ch04.two_fluid_couette(np.array([5e-4, 1.5e-3]), 1e-3, 0.1, 1e-3, 1e-3, 1.0))` · `print(
   ch04.pillbox_limit(9.09, 9.09, side_rate=2.0, volume_rate=5.0, l=np.array([1e-1, 1e-2, 1e-3]))["residual"])`.
   *expect:* q = 9.0909 W/m² in both layers, T_interface = 302.24 K · τ = 0.9901 Pa in both, u_interface = 0.9901 m/s ·
   residuals ∝ l (0.7, 0.07, 0.007) → 0.
8. `nb.note` — **N117 [B]** "**No slip and no temperature jump** at a solid wall — not consequences of conservation
   laws, but of experiment:" equation `\mathbf u_1\cdot\mathbf t=0\ \ (\text{relative to the wall}),\qquad T_1=T_2`, no ref.
   "Exceptions: superfluid helium (no viscosity), super-hydrophobic textured surfaces (apparent slip), and rarefied gases
   whose mean free path is not small compared with the flow (Knudsen number Kn ≳ 0.01, Ch. 1). Figure: a Couette profile
   with a Navier slip length b (our extension, labelled): u = U(y + b)/(h + 2b) — invisible for b ≪ h."
9. `nb.primer("moving level set and its normal speed", "A surface can be described as the set of points where a function
   is zero, η(x, t) = 0 (a level set, Ch. 2 P75). Its unit normal is n = ∇η/|∇η| (pointing toward increasing η). If it
   moves, a point riding on it keeps η = 0, and the chain rule gives its speed along n: −(∂η/∂t)/|∇η|. Only this normal
   speed is defined — sliding along the surface does not change the surface.", code="V = 1.5                                   # a wall moving at 1.5 m/s along +x\neta = lambda x, t: x[0] - V*t              # the plane x = Vt\nprint(ch04.surface_normal_speed(eta, np.array([2.0, 0.3, 0.0]), 1.0))   # 1.5: it moves at V along n = +x")`
   — *expect:* `1.5`.
10. `nb.note` — **N118 [B]** "An observer riding on the moving surface at velocity u_s sees η stay zero:" equation
    `d\eta/dt=\partial\eta/\partial t+(\mathbf u_s\cdot\nabla)\eta=0\quad\text{on}\quad\eta=0`, ref "4.90".
11. `nb.derivation("D29", …)` — Part F D29 (8 steps), ref "4.91".
12. `nb.note` — **C14's result, restated** — "**The kinematic boundary condition** — no fluid crosses the surface:"
    equation `\partial\eta/\partial t+(\mathbf u\cdot\nabla)\eta\equiv D\eta/Dt=0\quad\text{on}\quad\eta=0`, ref "4.91". "A free
    surface, an interface between immiscible fluids and a solid wall all obey it."
13. `nb.note` — **N119 [B]** "If mass *does* cross (evaporation, a shock), the normal velocity of the fluid relative to
    the surface is" equation `(u_{rel})_n=\mathbf u\cdot\mathbf n-\mathbf u_s\cdot\mathbf n=\big(\mathbf u\cdot\nabla\eta+\partial\eta/\partial t\big)/\lvert\nabla\eta\rvert=(1/\lvert\nabla\eta\rvert)D\eta/Dt`,
    ref "4.92".
14. `nb.note` — **N120 [B]** "and the mass flux per unit area through it is" equation
    `(\rho/\lvert\nabla\eta\rvert)D\eta/Dt\quad\text{on}\quad\eta=0`, ref "4.93". "Zero flux is again Dη/Dt = 0. Moving shocks:
    Ch. 15."
15. `nb.note` — **N121 [B]** "**Surface tension as free energy.** At constant temperature the useful quantity is the
    Helmholtz free energy per unit mass (a Legendre device, like Ch. 1's Gibbs g = h − Ts, P50):" equation `f=e-Ts`, ref
    "4.94".
16. `nb.note` — **N122 [B]** equation `df=de-T\,ds-s\,dT`, ref "4.95". "At constant T and reversibly (T ds = de + p dv):
    df = −p dv — the work done on the system. (sympy line: for a perfect gas at fixed T, e is constant and
    s = C_v ln T + R ln v + const, so df/dv = −RT/v = −p.)"
17. `nb.note` — **N123 [B]** "For two fluids and their interface of area A the free energy is" equation
    `F=\rho_1V_1f_1+\rho_2V_2f_2+A\sigma`, ref "4.96". "so σ is free energy per unit area; with σ > 0 (immiscible fluids) the
    system lowers F by shrinking A — the interface contracts as far as the constraints allow. Figure: the surface area of
    a spheroid of fixed volume vs its aspect ratio (`ch04.spheroid_area`) is smallest at aspect 1 — a sphere."
18. `nb.note` — **N124 [C]** "Where σ varies along an interface (with temperature or a surfactant) the surface pulls
    toward high σ and drags fluid with it — *Marangoni* flows (tears of wine). Not treated in the book; see Ch. 16 and
    interfacial-flow texts."
19. `nb.note` — **N125 [B]** "**The Laplace jump derived** (stated with numeric checks; the book writes the construction):
    cut a small cap $z=x^2/2R_1+y^2/2R_2$ at height ζ; its rim is an ellipse with semi-axes $\sqrt{2R_1\zeta}$, $\sqrt{2R_2\zeta}$.
    The pressure difference Δp pushes on it with" equation `(F_p)_z=-\pi\Delta p\sqrt{2R_1\zeta}\sqrt{2R_2\zeta}`, ref "4.97".
20. `nb.note` — **N126 [B]** "…and the surface tension pulls along the rim ($\sigma\oint_C\mathbf t\times\mathbf n\,ds$) with a
    net upward force, for a small cap," equation
    `(F_{st})_z=\pi\sigma\sqrt{2R_1\zeta}\sqrt{2R_2\zeta}\Big(\frac1{R_1}+\frac1{R_2}\Big)`, ref "4.98". "> ⚠️ **Book slip:** the
    curve C is printed as ζ = x²/2R₁ − y²/2R₂; it must be **+** (the cap equation); the book's later lines use +."
21. `nb.note` — **N127 [B]** "Both forces are ∝ ζ; they balance when" equation
    `\Delta p=\sigma\Big(\frac1{R_1}+\frac1{R_2}\Big)`, ref "1.5". "— Laplace's law (1.5), which Ch. 1 stated, now derived.
    Higher pressure on the concave side. Number: a water drop of radius 1 mm (R₁ = R₂ = 1 mm, σ = 0.0728 N/m): Δp = 146 Pa.
    Capillary waves: Ch. 7."
22. `nb.note` — **N128 [B]** "**Capillary length** — where gravity and surface tension are equally strong:" equation
    `\ell_c=\big(\sigma/\rho g\big)^{1/2}`, no ref. "Water at 20 °C: √(0.0728/(998 × 9.81)) = 2.73 mm. Smaller than ℓ_c,
    surface tension wins (a drop is round); larger, gravity wins (a puddle is flat). Its square ratio is the Bond number
    (C15)."
23. `nb.note` — **N129 [B]** "**Ex. 4.7 — the meniscus at a vertical wall** (stated; the closed form is checked against an
    ODE solution instead of re-derived): a liquid meets a wall at contact angle θ; the surface height ζ(x) obeys"
    equation `\Big(\frac{\rho g}{2\sigma}\Big)\zeta^2+\big(1+\zeta'^2\big)^{-1/2}=1,\qquad h^2=\frac{2\sigma}{\rho g}(1-\sin\theta)`,
    no ref. "(the second is the wall height). Water on clean glass (θ = 0): h = √2 ℓ_c = 3.86 mm. > ⚠️ **Book slips in Ex.
    4.7:** the separated equation drops a minus sign (the slope is negative), and 'evaluated when η = h/δ' should read
    γ = h/δ — here γ = ζ/δ, yet another γ; the final closed form is right (verified by differentiation in the test)."
24. `nb.worked_example("three edges", "1. **A moving wall** η = x − Vt, V = 1 m/s, with fluid at u = (1, 0.3, 0) m/s on
    it: ∂η/∂t = −1, ∇η = (1, 0, 0), so Dη/Dt = −1 + 1 × 1 = 0 ✓ — the fluid's normal speed (1 m/s) equals the wall's; the
    tangential 0.3 m/s is not constrained by (4.91) (no-slip is a separate condition). 2. **A small water wave**
    η = a cos(kx − ωt) with ka = 0.1: the linearised condition ∂η/∂t = w at z = 0 holds exactly for the linear solution;
    the full condition (4.91) at the true surface leaves a residual of order (ka)² ≈ 0.01 of the phase speed — a quarter of
    that at ka = 0.05. 3. **A drop**: R = 1 mm, σ = 0.0728 N/m → Δp = 2σ/R = 146 Pa; ℓ_c = 2.73 mm.")`
25. `nb.code` — *code:* `eta_w = "moving_wall"` · `print(ch04.kinematic_bc_residual(eta_w, lambda x, t: np.array([1.0,
    0.3, 0.0]), np.array([0.5, 0.0, 0.0]), 0.5, V=1.0), ch04.surface_normal_speed(eta_w, np.array([0.5, 0, 0]), 0.5,
    V=1.0))` · `for ka in (0.2, 0.1, 0.05): k = 1.0; a = ka/k; c = np.sqrt(9.81/k)`; `res = max(abs(ch04.wave_kinematic_residual(x,
    0.0, a, k)) for x in np.linspace(0, 2*np.pi, 64))`; `print(ka, res/c, ch04.wave_kinematic_residual(0.3, 0.0, a, k,
    full=False))` · `print(ch04.laplace_jump_from_balance(0.0728, 1e-3, 1e-3), ch01.laplace_pressure_jump(0.0728, 1e-3,
    1e-3))` · `print(ch04.cap_pressure_force(146.0, 1e-3, 1e-3, 1e-6), ch04.cap_surface_tension_force(0.0728, 1e-3, 1e-3,
    1e-6))` · `print(ch04.capillary_length(0.0728, 998.0), ch04.meniscus_height(0.0, 0.0728, 998.0))`. *expect:* `0.0 1.0`
    · (0.2, ≈0.04, 0.0), (0.1, ≈0.01, 0.0), (0.05, ≈0.0025, 0.0) — slope 2 · `145.6 145.6` · two forces of magnitude ≈
    9.17×10⁻⁷ N with opposite signs · `0.0027266 0.0038563`. *explain:* 1. the wall condition; 2. the full condition's
    residual for the linear wave shrinks like (ka)², while the linearised one is exact — the error made by linearising
    (Ch. 7); 3. Laplace from the cap balance equals Ch. 1's formula; 4. the two cap forces cancel; 5. capillary length
    and the meniscus height at θ = 0.
26. `nb.check_agree` — **from scratch (curation §7):** Dη/Dt by central differences for the wave at a surface point:
    `eta = lambda x, t: x[1] - a*np.cos(k*x[0] - w*t)` (η as the level-set function z − ζ(x, t)); `dt = de = 1e-6`;
    `deta_dt = (eta(p, dt) - eta(p, -dt))/(2*dt)`; `grad = [(eta(p + e, 0) - eta(p - e, 0))/(2*de) for e in
    np.eye(2)*de]`; `u = ch04.linear_wave_surface(p[0], p[1], 0.0, a, k)`; `res = deta_dt + u["u"]*grad[0] +
    u["w"]*grad[1]`; `assert np.isclose(res, ch04.kinematic_bc_residual(eta, …, p, 0.0), rtol=1e-6)`.
27. `nb.animation` — (`player="video"`, 60 frames, FAST 30) a linear deep-water wave (ka = 0.15, exaggerated in the
    title) with a row of particles on the surface moving on their orbits (orange) and a second row 0.3λ deep (grey, smaller
    orbits); the surface normal and the two normal speeds (fluid u·n teal, surface u_s·n rose) as arrows at one marked
    point — always equal. see/read/change: *see:* "surface particles bob and circle but stay on the moving surface; the
    teal and rose arrows coincide." *read:* "that coincidence *is* (4.91): the surface moves normal to itself exactly as
    fast as the fluid on it." *change:* "…an evaporating surface: the rose arrow would lag the teal one by the evaporation
    speed and (4.92) would be nonzero."
28. `nb.figure` — two panels, the pillbox cases of N116: (a) temperature through the composite wall — two straight
    segments with a kink at the interface, heat-flux arrows of equal length on both sides; (b) two-fluid Couette —
    velocity with a kink at the interface, shear-stress arrows of equal length; plus a thin third panel: Couette with a
    slip length b = 0, 0.1h, 0.5h (N117). *see:* "kinked profiles, equal arrows." *read:* "what is continuous is the
    *flux* (heat, momentum); the gradients jump in inverse proportion to k or μ." *change:* "…k₂ = k₁: no kink, one
    straight line."
29. `nb.figure` — meniscus profiles ζ(x)/ℓ_c for θ = 0°, 30°, 60°, 90° from `ch04.meniscus_profile_x` (lines) with
    `ch04.meniscus_profile_ode` markers on top, and the wall heights h/ℓ_c = √(2(1 − sin θ)). *see:* "curves rising to the
    wall, lower for larger contact angle, flat at 90°." *read:* "surface tension lifts the liquid within a few capillary
    lengths of the wall; gravity flattens it farther away." *change:* "…mercury on glass (θ ≈ 140°): 1 − sin θ < 1 still,
    but the meniscus is depressed — the formula's h² gives |h|, the sign comes from cos θ < 0."
30. `nb.md` — "> ⚠️ **Common confusion:** '(4.91) is the no-slip condition.' It only matches the **normal** velocities; the
    tangential velocity at a free surface is set by the stress condition (no shear from the air), at a solid wall by
    no-slip."
31. `nb.md` — **What would change if…** "…we measured every length in units of the body's size, every speed in units of
    the free stream, every time in units of the forcing period? The equations would lose their dimensions and only a
    few numbers would remain — and two flows with the same numbers would be the same flow (C15)."

### A.11 §4.11 Dimensionless Forms of the Equations and Dynamic Similarity — R12, C15 (+N130–N156 · D30 · E9)
1. `nb.section("4.11", "Dimensionless Forms of the Equations and Dynamic Similarity", intro="**What is this section
   about?** Scale every variable by the problem's own sizes and Navier–Stokes keeps only three numbers — St, Fr, Re. Two
   flows with equal numbers and the same shape are the same flow: model testing works, data collapse, and every later
   chapter names its regimes by these numbers.")`
2. `nb.recap("R12", "Dimensional analysis of sphere drag", "The drag F_D on a sphere of diameter d in a fluid (ρ, μ)
   moving at U depends on 5 variables with 3 dimensions, so on 5 − 3 = 2 groups (Ch. 1's Π theorem):
   $\frac{F_D}{\rho U^2d^2}=\Psi\Big(\frac{\rho Ud}{\mu}\Big)$ or equally $\frac{F_D\rho}{\mu^2}=\Phi\Big(\frac{\mu}{\rho Ud}\Big)$ *(Eq. 4.99)* —
   two choices of repeating variables; the second group is the first times Re².", where="Ch. 1 §1.11")` — followed by
   `nb.code`: `from fluidpy.core.dimensional import pi_groups, SPHERE_DRAG` · `print(pi_groups(SPHERE_DRAG,
   solution="F_D", repeating=["rho", "U", "d"]))` · `print(pi_groups(SPHERE_DRAG, solution="F_D", repeating=["rho",
   "mu", "d"]))`. *expect:* two group sets as ch01 prints them (F_D/(ρU²d²) with Re; F_Dρ/μ² with Re).
3. `nb.core("C15", "When is a model the real thing? Dimensionless Navier–Stokes (4.101)", question="A ship designer tows
   a 4-m model in a tank to predict the drag of a 100-m ship. What must be the same in the tank for the answer to be
   right — and can everything be matched at once?")`
4. `nb.md` — **The problem in plain words:** "Wind tunnels, towing tanks, rotating tanks and computer runs all study small
   or slow copies of real flows. The copy is only useful if the physics is the same. Equations tell us exactly which
   combinations of sizes must agree. The same combinations — Reynolds, Froude, Richardson, Rossby, Mach — classify the
   flows of every later chapter: creeping or boundary-layer flow (Ch. 8–9), waves (Ch. 7), stratified instability
   (Ch. 11), geostrophic flow (Ch. 13)."
5. `nb.md` — **The idea:**
   ```
   x = l x*,  t = t*/Ω,  u = U u*,  p − p∞ = ρU² p*,  g = g g*                                (4.100)
   (4.39b) ÷ (ρU²/l)  →  [Ωl/U] ∂u*/∂t* + (u*·∇*)u* = −∇*p* + [gl/U²] g* + [μ/ρUl] ∇*²u*     (4.101)
                             St (blue)                           1/Fr² (grey)   1/Re (rose)
   same St, Fr, Re + same shape + same boundary conditions  →  same dimensionless solution
   ```
6. `nb.note` — **N130 [B]** "**Two routes to the groups.** Dimensional analysis (R12) finds *that* groups exist; scaling
   the equations tells *which* groups and what each measures (a ratio of two terms). Dynamic similarity = geometric
   similarity + equal groups. The figure at the end of the block shows two dimensional Stokes-layer solutions collapsing
   onto one dimensionless curve."
7. `nb.note` — **N132 [C]** "The parameters of a general unsteady flow: a length l, a speed U, a frequency Ω (pulsating
   flow in a tube, a swimmer's stroke, a turbine's rotation), and the fluid's ρ and μ."
8. `nb.primer("scaled variables and the chain rule", "If t* = Ωt then ∂/∂t = Ω ∂/∂t* (each second of t is Ω units of t*),
   and if x* = x/l then ∂/∂x = (1/l)∂/∂x*, ∂²/∂x² = (1/l²)∂²/∂x*². Scaling a derivative brings out its scale factor —
   the move that turns each term of an equation into (a size) × (a dimensionless term of order one).", code="t, ts, Om =
   sp.symbols('t t_s Omega', positive=True); f = sp.Function('f')\nexpr = sp.diff(f(Om*t), t)                 # d/dt of f(Ωt)\nprint(sp.simplify(expr))                    # Omega*Subs(Derivative(f(_xi), _xi), _xi, Omega*t)")`
   — *expect:* sympy's `Omega*Subs(Derivative(f(_xi_1), _xi_1), _xi_1, Omega*t)` (the factor Ω appears).
9. `nb.note` — **N133 [B]** "**Scaled variables** (Ω an imposed frequency):" equation
   `x_i^*=x_i/l,\quad t^*=\Omega t,\quad u_j^*=u_j/U,\quad p^*=(p-p_\infty)/\rho U^2,\quad g_j^*=g_j/g`, ref "4.100". "> ⚠️ The
   time and pressure scales change inside this section: t* = Ωt here but Ut/l for steady boundary conditions and in (4.109);
   p scaled by ρU² here, but μU/l (slow viscous flow) or ρgl (hydrostatics) are also used. `ch04.Scales` records which
   convention each use takes."
10. `nb.derivation("D30", …)` — Part F D30 (8 steps), ref "4.101", with `check_src` (optional): `ch04.nondimensional_ns_sym()`
    printed, and `ch04.nondimensional_ns_coefficients()` = {unsteady: Ω l/U, advective: 1, pressure: 1, gravity: g l/U²,
    viscous: μ/(ρUl)}.
11. `nb.md` — "> ⚠️ **Book slip (reference):** the text after (4.101) says the groups are 'in [,]-brackets in (4.100)';
    they are in (4.101), $\big[\frac{\Omega l}{U}\big]\frac{\partial\mathbf u^*}{\partial t^*}+(\mathbf u^*\cdot\nabla^*)\mathbf u^*=-\nabla^*p^*+\big[\frac{gl}{U^2}\big]\mathbf g^*+\big[\frac{\mu}{\rho Ul}\big]\nabla^{*2}\mathbf u^*$."
12. `nb.note` — **N134 [B]** "**Strouhal number** — unsteady over advective acceleration:" equation
    `\mathrm{St}\equiv\frac{\text{unsteady acceleration}}{\text{advective acceleration}}\propto\frac{\partial u/\partial t}{u(\partial u/\partial x)}\propto\frac{\Omega U}{U^2/l}=\frac{\Omega l}{U}`,
    ref "4.102". "Vortex shedding behind a cylinder happens at St ≈ 0.2 (Ch. 9)."
13. `nb.note` — **N135 [B]** "**Reynolds number** — inertia over viscous force:" equation
    `\mathrm{Re}\equiv\frac{\text{inertia force}}{\text{viscous force}}\propto\frac{\rho u(\partial u/\partial x)}{\mu(\partial^2u/\partial x^2)}\propto\frac{\rho U^2/l}{\mu U/l^2}=\frac{\rho Ul}{\mu}`,
    ref "4.103". "Our numbers (code): a bacterium (1 µm, 30 µm/s, water) Re ≈ 3×10⁻⁵ · a 1-cm marble at 1 m/s in water
    10⁴ · a car (4 m, 30 m/s, air) 8×10⁶ · an ocean liner (300 m, 10 m/s, water) 3×10⁹ · the Gulf Stream (100 km, 1 m/s)
    10¹¹."
14. `nb.note` — **N136 [B]** "**Froude number** — the square root of inertia over gravity:" equation
    `\mathrm{Fr}\equiv\Big[\frac{\text{inertia force}}{\text{gravity force}}\Big]^{1/2}\propto\Big[\frac{\rho u(\partial u/\partial x)}{\rho g}\Big]^{1/2}\propto\Big[\frac{\rho U^2/l}{\rho g}\Big]^{1/2}=\frac{U}{\sqrt{gl}}`,
    ref "4.104". "Gravity matters dynamically only with a free surface or density differences — a submarine deep down
    does not feel Fr; a ship does (waves, Ch. 7)."
15. `nb.note` — **N137 [B]** "**Stratified flows** use the reduced gravity g′ = g(ρ₂ − ρ₁)/ρ₁ or the buoyancy frequency N:"
    equation `\mathrm{Fr}'=\frac{U}{\sqrt{g'l}}\ \text{or}\ \frac{U}{Nl},\qquad\mathrm{Ri}=\frac1{\mathrm{Fr}'^2},\qquad\mathrm{Ri}_g=\frac{N^2}{(dU/dz)^2}`,
    ref "4.105". "Numbers: an ocean thermocline Δρ/ρ = 10⁻³ (g′ = 9.81×10⁻³ m/s²), l = 100 m, U = 0.1 m/s → Ri = 98 —
    stratification dominates. The atmosphere with dT/dz = −6.5 K/km (Kundu's Γ; meteorology's lapse rate Γ_met = +6.5
    K/km) at 288 K: $N^2=\frac gT\big(\frac{dT}{dz}+\frac g{C_p}\big)=\frac{9.81}{288}(-6.5+9.76)\times10^{-3}=1.11\times10^{-4}$ s⁻²
    (stable since dT/dz = −6.5 > Γ_a = −9.8 K/km ⇔ Γ_met = 6.5 < 9.8 K/km); with a wind shear dU/dz = 0.01 s⁻¹, Ri_g = 1.1
    — above the ¼ of Ch. 11's stability criterion."
16. `nb.note` — **N138 [C]** "Under dynamic similarity **every** dimensionless output agrees too — local ones (a pressure
    coefficient at a point) and overall ones (a drag coefficient)."
17. `nb.note` — **N139 [B]** "**Pressure coefficient**:" equation
    `C_p\equiv\frac{p-p_\infty}{\tfrac12\rho U^2}=\Psi\Big(\mathrm{St},\mathrm{Fr},\mathrm{Re};\frac{\mathbf x}{l},\Omega t\Big)`, ref "4.106".
    "Demo: the ideal cylinder's surface C_p = 1 − 4 sin²θ is identical for five (U, a) pairs — (1 m/s, 1 m), (10, 0.1),
    (0.5, 3), … (`ch04.pressure_coefficient` on `ch03.cylinder_flow`)."
18. `nb.note` — **N140 [C]** "With steady boundary conditions the only time scale is l/U: t* = Ut/l and St drops out —
    but a flow can still become unsteady by itself (vortex shedding), at a frequency that scales with U/l (Ch. 9)."
19. `nb.note` — **N141 [B]** "A body **purely oscillating** with amplitude l at frequency Ω has U = lΩ, so St = 1,
    Re = Ωl²/ν and Fr = Ω(l/g)^{1/2} (`ch04.Scales.from_oscillation`) — the Stokes-layer and Womersley parameters of Ch. 8
    and Ch. 16."
20. `nb.note` — **N142 [B]** "**Drag coefficient**:" equation `C_D\equiv\frac{F_D}{\tfrac12\rho U^2A}`, ref "4.107".
21. `nb.note` — **N143 [B]** "**Lift coefficient**:" equation `C_L\equiv\frac{F_L}{\tfrac12\rho U^2A}`, ref "4.108". "A is a
    reference area: the frontal area πd²/4 for a sphere or a car, the planform (span × chord) for a plate or a wing
    (`ch04.reference_area`). Wings: Ch. 14."
22. `nb.note` — **N131 [B]** "**Sphere drag** (our version of Fig. 4.21, below): C_D = D/(½ρU²·πd²/4) against
    Re = ρUd/μ — 24/Re (Stokes, Ch. 8) at small Re, a plateau near 0.4–0.5 for 10³ ≲ Re ≲ 2×10⁵, then the *drag crisis*
    (Ch. 9). We use the published correlation of Morrison (2013) (V5 benchmark). Number: d = 1 cm, U = 1 m/s in water →
    Re = 10⁴, C_D = 0.39."
23. `nb.note` — **N144 [C]** "For a ship, C_D = C_D(Fr, Re); far from any free surface and at low Mach number,
    C_D = C_D(Re) alone — the bridge to Ex. 4.8 below."
24. `nb.note` — **N145 [B]** "**Compressible flow scalings** (t* = Ut/l now):" equation
    `t^*=Ut/l,\qquad p^*=(p-p_\infty)/\rho_oU^2,\qquad\rho^*=\rho/\rho_o`, ref "4.109".
25. `nb.note` — **N146 [B]** "Using (4.9)'s §4.11 re-display, $\nabla\cdot\mathbf u=-\frac1\rho\frac{D\rho}{Dt}=-\frac{1}{\rho c^2}\frac{Dp}{Dt}$
    (isentropic density changes dp = c²dρ — a general identity, not the incompressibility condition), the scalings give"
    equation `\nabla^*\cdot\mathbf u^*=-\Big[\frac{U^2}{c^2}\Big]\frac1{\rho^*}\frac{Dp^*}{Dt^*}`, ref "4.110". "The departure from
    ∇·u = 0 is of order M²: 0.09 at M = 0.3 (`ch04.compressibility_parameter`)."
26. `nb.note` — **N147 [B]** "**Mach number** (the square root of inertia over compressibility forces):" equation
    `M\equiv U/c`, ref "4.111". "Incompressible treatment is good below M ≈ 0.3 (C02's N12): ≈ 100 m/s in air. Ch. 14, 15."
27. `nb.note` — **N148 [B]** "**Energy equation in enthalpy form** (stated; the book calls it 'a mild revision' of (4.60)
    and does not write the steps; a sympy one-liner proves the two agree, `ch04.energy_forms_sym()` → 0):" equation
    `\rho\frac{Dh}{Dt}=\frac{Dp}{Dt}+\rho\varepsilon+\frac{\partial}{\partial x_i}\Big(k\frac{\partial T}{\partial x_i}\Big)`, ref "4.112".
28. `nb.note` — **N149 [B]** "**Thermal scalings**:" equation
    `\varepsilon^*=\frac{\rho_ol^2}{\mu_oU^2}\varepsilon,\quad\mu^*=\mu/\mu_o,\quad k^*=k/k_o,\quad T^*=\frac{T-T_o}{T_w-T_o}`, ref "4.113".
29. `nb.note` — **N150 [B]** "**Dimensionless energy equation** (dh ≅ C_p dT):" equation
    `\rho^*\frac{DT^*}{Dt^*}=\Big[\frac{U^2}{C_p(T_w-T_o)}\Big]\frac{Dp^*}{Dt^*}+\Big[\frac{\mu_oU}{\rho_oC_p(T_w-T_o)l}\Big]\varepsilon^*+\Big[\frac{k_o}{\rho_oC_pUl}\Big]\frac{\partial}{\partial x_i^*}\Big(k^*\frac{\partial T^*}{\partial x_i^*}\Big)`,
    ref "4.114". "The brackets are Ec, Ec/Re and 1/(Pr Re) (`ch04.nondimensional_energy_coefficients()`). > ⚠️ The book
    introduces it with 'those defined in (4.106), (4.107)'; it means the scalings (4.109) and the equation (4.112)."
30. `nb.note` — **N151 [B]** "**Eckert number**:" equation `\mathrm{Ec}\equiv\frac{U^2}{C_p(T_w-T_o)}`, ref "4.115". "Small Ec
    ⇒ pressure work and viscous heating drop out and (4.112) becomes the Boussinesq heat equation (4.89),
    $\frac{DT}{Dt}=\kappa\nabla^2T$ — the hidden assumption of C13."
31. `nb.note` — **N152 [B]** "**Prandtl number** — momentum diffusivity over heat diffusivity:" equation
    `\mathrm{Pr}\equiv\frac{\nu}{\kappa}=\frac{\mu_oC_p}{k_o}`, ref "4.116". "Our values from `ch01.FLUIDS` (`ch04.prandtl_of`):
    air ≈ 0.71, water at 20 °C ≈ 7.0; kinetic theory for a monatomic gas 2/3; Eucken's estimate 4γ/(9γ − 5) = 0.74 for
    γ = 1.4. Thermal boundary layers: Ch. 9; convection: Ch. 11."
32. `nb.note` — **N153 [B]** "**Weber number** (a ratio of *forces*, not per volume):" equation
    `\mathrm{We}\equiv\frac{\rho U^2l}{\sigma}`, ref "4.117". "Break-up of drops and jets (Ch. 7, 16)."
33. `nb.note` — **N154 [B]** "**Bond number**:" equation `\mathrm{Bo}\equiv\frac{\rho l^2g}{\sigma}=(l/\ell_c)^2`, ref "4.118". "(the
    capillary length of C14)."
34. `nb.note` — **N155 [B]** "**Capillary number**:" equation `\mathrm{Ca}\equiv\frac{\mu U}{\sigma}=\mathrm{We}/\mathrm{Re}`, ref
    "4.119". "Coating flows and flow in porous media (Ch. 16)."
35. `nb.note` — **N156 [B]** "**Ex. 4.8 — the ship model** (stated with our own numbers; the book's are in the private test
    data). A 1:25 model of a 100-m ship designed for 10 m/s. **Froude matching** $U_m=U_p\sqrt{l_m/l_p}=10/5=2$ m/s makes
    the wave pattern similar — but then Re_m/Re_p = (l_m/l_p)^{3/2} = 1/125: the model's boundary layers are far too
    viscous. Froude's way out: split the drag into friction (depends on Re, from a flat-plate friction line) and wave drag
    (depends on Fr, scaled by $(\rho_p/\rho_m)(l_p/l_m)^2(U_p/U_m)^2$ = λ³ for equal densities). Numbers in the worked
    example; `ch04.ship_drag_extrapolation`. > ⚠️ The book's total 9.14×10⁵ N is 9.15×10⁵ N without the intermediate
    rounding — a detail, not an error in the method."
36. `nb.worked_example("similarity by hand", "1. **Sphere**: d = 1 cm, U = 1 m/s, water ν = 10⁻⁶ m²/s: Re = Ud/ν = 10⁴;
    C_D ≈ 0.39 (plateau) → drag ½ρU²C_D πd²/4 = 0.5 × 1000 × 1 × 0.39 × 7.85×10⁻⁵ = 0.015 N. 2. **The same Re in air**
    (ν = 1.5×10⁻⁵): a 1-cm sphere needs U = 15 m/s — then its drag coefficient is the same 0.39. 3. **Ship**: λ = 1/25,
    U_p = 10 m/s → U_m = 2 m/s (Fr = 0.32 both); Re_p = 10 × 100/10⁻⁶ = 10⁹, Re_m = 2 × 4/10⁻⁶ = 8×10⁶ — off by 125.
    4. **Drag extrapolation** (model wetted area 4 m², prototype 2500 m²; friction line C_f = 0.075/(log₁₀Re − 2)², ITTC
    1957, gloss): C_f,m = 0.00312, C_f,p = 0.00153. Model measures 40 N total; friction ½ρU²SC_f = ½ × 1000 × 4 × 4 ×
    0.00312 = 25.0 N; wave 15.0 N. Wave drag scales by (1025/1000) × 25² × 5² = 16 016 → 240.9 kN. Prototype friction
    ½ × 1025 × 100 × 2500 × 0.00153 = 196.1 kN. Total ≈ 437 kN. Scaling the whole 40 N by 16 016 would give 641 kN — 47 %
    too much, because friction does not scale with Fr.")`
37. `nb.code` — *code:* `for name, (U, l, nu) in {"bacterium": (3e-5, 1e-6, 1e-6), "marble": (1.0, 0.01, 1e-6), "car":
    (30.0, 4.0, 1.5e-5), "liner": (10.0, 300.0, 1e-6), "Gulf Stream": (1.0, 1e5, 1e-6)}.items(): print(f"{name:12s} Re =
    {ch04.reynolds_number(U, l, nu):.1e}")` · `print(ch04.sphere_drag_coefficient(1e4), ch04.sphere_drag_coefficient(0.1,
    model="stokes"))` · `print(ch04.froude_scaled_speed(10.0, 100.0, 4.0), ch04.froude_number(10.0, 100.0),
    ch04.froude_number(2.0, 4.0))` · `mp = ch04.model_prototype(100.0, 10.0, 1/25)`; `print(mp["U_m"], mp["Re_ratio"],
    mp["matched"])` · `print(ch04.ship_drag_extrapolation(100.0, 10.0, 2500.0, 1/25, 40.0, 0.00312, 0.00153))` ·
    `print(ch04.richardson_number(9.81e-3, 100.0, 0.1), ch04.prandtl_of("air"), ch04.prandtl_of("water"),
    ch04.eucken_prandtl(1.4))` · `print(ch04.nondimensional_ns_coefficients())`. *expect:* Re = 3.0e-05, 1.0e+04,
    8.0e+06, 3.0e+09, 1.0e+11 · `0.3926 240.0` · `2.0 0.3193 0.3193` · `2.0 0.008 {'Fr': True, 'Re': False}` · U_m 2.0,
    D_m_friction ≈ 25.0, D_m_wave ≈ 15.0, D_p_wave ≈ 2.409e5, D_p_friction ≈ 1.961e5, D_p_total ≈ 4.37e5 N, Re_ratio 0.008
    · `98.1 ≈0.71 ≈7.0 0.7368` · the five sympy coefficients. *explain:* 1. Reynolds numbers across 16 orders of magnitude;
    2. the sphere's drag coefficient in two regimes; 3. Froude scaling keeps Fr equal; 4. the model–prototype summary: Fr
    matched, Re not (ratio λ^{3/2}); 5. Froude's extrapolation; 6. the stratified and thermal numbers; 7. the brackets
    of (4.101).
38. `nb.check_agree` — **from scratch (curation §7):** Re and Fr for the ship by hand (`U*l/nu`, `U/np.sqrt(9.81*l)`) and
    the sympy substitution of (4.100) into the x-momentum of (4.39b): `x_, t_, l, U, Om, rho, mu, g = sp.symbols(...)`;
    `us = sp.Function('u_s')`; `u = U*us(x_/l, Om*t_)`; `lhs = rho*(sp.diff(u, t_) + u*sp.diff(u, x_))`; `visc =
    mu*sp.diff(u, x_, 2)`; divide by ρU²/l and `sp.simplify` → coefficients Ωl/U on the unsteady term and μ/(ρUl) on the
    viscous one; `assert` they equal `ch04.nondimensional_ns_coefficients()["unsteady"]` and `["viscous"]`.
39. `nb.figure` — **same equation, same solution** (N130): left, two dimensional Stokes-layer profiles u(y) — (U = 1 m/s,
    water, t = 1 s) and (U = 5 m/s, ν = 10⁻⁵ m²/s, t = 0.1 s) — on different axes; right, the same data as u/U vs
    y/(2√(νt)) lying on one erfc curve. *see:* "two different curves on the left, one on the right." *read:* "after
    scaling, the equation has no parameters left, so every Stokes layer is the same curve." *change:* "…a plate
    oscillating instead of started: a new group (Ωl²/ν) appears and the collapse needs it to match too."
40. `nb.plotly` — `slider_figure` over Re = 10⁻¹…10⁶ (log, 36 steps): the Morrison C_D(Re) curve (black), the Stokes
    line 24/Re (grey dashed) and the current point (rose). Title "One curve for every sphere: C_D depends on Re alone". see/
    read/change: "…at the slider's Re the point shows C_D; any sphere, any fluid, any speed with that Re sits on it."
41. `nb.figure` — (N131) synthetic "experiments": 60 (d, U, fluid) combinations (water, air, glycerine; d = 1 mm–10 cm)
    with C_D from Morrison × (1 ± 3 % noise): left panel C_D vs U (a scattered cloud), right panel C_D vs Re (a single
    curve). *see:* "chaos on the left, order on the right." *read:* "plotting against the right group removes the
    dimensional scatter — that is dynamic similarity measured." *change:* "…a rough sphere: the drag crisis moves to lower
    Re (a new group, roughness/d)."
42. `nb.figure` — (N156) the ship-drag waterfall: model total → friction + wave → wave × 16 016 → + prototype friction =
    prototype total, beside the naive "total × 16 016" bar (grey, hatched). *see:* "the corrected bar much shorter than the
    naive one." *read:* "friction scales with Re, waves with Fr; the model cannot match both, so each part is scaled by its
    own law." *change:* "…a model in a fluid of much lower ν (none practical exists): both Fr and Re could match — the
    reason towing tanks live with the correction."
43. `nb.explainer("dynamic_similarity_models", heading="When does a model behave like the real thing?", why="A prototype
    and a model side by side, the groups of
    $\big[\frac{\Omega l}{U}\big]\frac{\partial\mathbf u^*}{\partial t^*}+(\mathbf u^*\cdot\nabla^*)\mathbf u^*=-\nabla^*p^*+\big[\frac{gl}{U^2}\big]\mathbf g^*+\big[\frac{\mu}{\rho Ul}\big]\nabla^{*2}\mathbf u^*$
    (4.101) as paired bars with 'matched' badges, and the consequences: sphere data collapsing on C_D(Re), a ship's drag
    split and scaled. Modes: sphere, ship, stratified flow (Fr′, Ri), rotating tank (Rossby number, a forward pointer to
    Ch. 13).", tries=["Ship preset 1:25: match Fr and read the Re mismatch ×125 in Explain.", "Sphere mode: toggle the
    axes from dimensional to Re — the points collapse.", "Try to match Re too by changing the model fluid: which fluid
    would you need?", "Rotating-tank mode: match the Rossby number U/(Ωl) of an atmospheric system in a 1-m tank."])`
44. `nb.md` — **What would change if…** "…the flow had rotation and stratification together, as the ocean and atmosphere
    do? Two more groups join — the Rossby number U/(Ωl) from the Coriolis term of (4.45) and the Richardson number above —
    and the regime map of Ch. 13 is drawn in their plane."

### A.12 End matter — S01, S02, summary
1. `nb.pointer("Exercises (S01): the book's Exercises 4.1–4.63 are for practice. The derivations the text defers to
   Exercises 4.7, 4.8, 4.30, 4.38, 4.42, 4.43, 4.45, 4.46, 4.47 and 4.50 are written out above in our own words (D04,
   D08, D13, D15, D18, D20, D21, D23, D24).")`
2. `nb.pointer("Literature (S02): the chapter's references — Aris (isotropic tensors, used in C07), Batchelor and Lamb
   (classical treatments), Spiegel & Veronis (the Boussinesq approximation, C13) — are listed at the end of the book's
   chapter.")`
3. `nb.summary(clicked=[…], feeds_forward=[…], left_out=[…])` — **clicked** (one per CORE, each with its equation):
   C01 "A box of any shape and motion keeps a mass budget:
   $\frac{d}{dt}\int_{V^*}\rho\,dV+\int_{A^*}\rho(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=0$ (4.5) — flux counts relative to the walls."
   · C02 "At a point, $\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0$ (4.7); incompressible means each particle keeps
   its density, so ∇·u = 0 (4.10)." · C03 "In 2-D one scalar ψ carries the flow: $u=\partial\psi/\partial y$,
   $v=-\partial\psi/\partial x$; contours are streamlines and Δψ is the flux between them." · C04 "Forces follow from
   fluxes: (4.17) weighs the drag on a bar from its wake." · C05 "Along one streamline of a steady, frictionless,
   constant-density flow $\tfrac12U^2+gz+p/\rho$ is constant (4.19)." · C06 "Newton for every continuum:
   $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ (4.24), stress divergence over the first index." · C07 "A
   Newtonian fluid feels only its strain rate: $\tau_{ij}=-p\delta_{ij}+2\mu S_{ij}+\lambda S_{mm}\delta_{ij}$ (4.31), two
   constants." · C08 "$\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ (4.39b) is a budget of five
   accelerations; each flow wakes only some of them." · C09 "A rotating observer adds −2Ω × u′ (Coriolis) and
   −Ω × (Ω × x′) (centrifugal) to the forces of (4.45); the ball flies straight, the floor turns." · C10 "Viscosity turns
   kinetic energy into heat at the rate ε ≥ 0 (4.58), a sum of squares — never the reverse." · C11 "Steady inviscid
   barotropic flow keeps B = ½u² + ∫dp/ρ + gz constant along streamlines and vortex lines (4.71), everywhere if
   irrotational (4.72)." · C12 "Irrotational unsteady flow obeys
   $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int\frac{dp}{\rho}+gz=$ const (4.75): accelerating fluid needs
   pressure even at rest." · C13 "Boussinesq keeps density only where it meets gravity: (4.86) and
   $\frac{DT}{Dt}=\kappa\nabla^2T$ (4.89)." · C14 "A surface made of fluid moves with it: Dη/Dt = 0 on η = 0 (4.91); fluxes
   are continuous across interfaces; walls impose no-slip." · C15 "Scaled, Navier–Stokes keeps St, Fr, Re (4.101); equal
   groups make a model the real thing." **feeds_forward:** Ch. 5 vorticity equation from (4.68)/(4.39b) · Ch. 6
   potential flow with (4.72), (4.75) · Ch. 7 free-surface conditions (4.75), (4.91) · Ch. 8 exact solutions checked by
   `core.navier_stokes` · Ch. 9 momentum integral from (4.17) · Ch. 10 conservative forms (4.22), (4.53) · Ch. 11 Boussinesq
   instabilities, Ri · Ch. 12 ε and the kinetic-energy budget · Ch. 13 (4.45) with Boussinesq, Ro, Ri · Ch. 15 energy
   equation, M. **left_out:** non-Newtonian fluids (Ch. 16); the proof of the isotropic-tensor form (cited, Aris);
   Kelvin's theorem (Ch. 5); Marangoni flows (interfacial texts); anelastic equations for deep atmospheres (Ch. 13
   texts).

### A.13 Placement table (ID → notebook section → host block → call)
| Block | Items placed there, in notebook order |
|---|---|
| §4.1 section | N01 (`note`, tagged C01), N02 skeleton (`md` + `code`, tagged C06) |
| before C01 | R01 (`recap` + `code`) |
| C01 | N03, D01, N04, N05, N06 (`note`s after D01), primer "dataclasses and named results" |
| before C02 | R02 (`recap` + `code`) |
| C02 | primer "continuity of a function…", N07, D02, N08, primer "product rule for a divergence", D03, N09, N10, N11, N12 |
| C03 | N13, N14, D04, N15 (`note` + `plotly` 3-D), N16, N17, N18, **E2** |
| before C04 | R03 (`recap` + `code`), R04 (`recap`) |
| C04 | N26, primer "momentum flux…", N19, D05, N20, N21, N22, N23, N24, N25, primer "conservative force…", N27, N28 (⚠️ `md`), N29, N31, N32, primer "torque…", N82, N83, N84, **E1** |
| C05 | primer "sympy expand…", N30, D06, N104, N103, N105 |
| C06 | primer "tensor divergence over the first index", D07, N33, N34, N35, N36, N37, N02 (ledger row) |
| before C07 | R05 (`recap` + `code`), R06 (`recap`) |
| C07 | N39, N38, D08, N40, N41, N42, primer "isotropic fourth-order tensor", N43, N44, D09, primer "deviatoric…", D10, N45, N46, N47, N48, N49, N50, N51, **E3** |
| C08 | D11, N52 (+ ledger rows), primer "Schwarz's theorem", D12, N53, primer "curl of a curl identity", D13, N54, N55, primer "complementary error function erfc", **E4** |
| C09 | N56, primers "derivative of a rotating unit vector", "product rule for a cross product", D14, N57, D15, N58, N59, D16, N60, primer "latitude…", D17, N61, N62, N63, D18, N64, N65, **E5** |
| before C10 | R07 (`recap` + `code`), R08 (`recap`) |
| C10 | primer "power of a force…", N66, D19, N67, N68, N69, N70, N71, N72, N73, D20, N74, N75, primer "chain rule for the kinetic energy", D21, N76, D22, primer "completing the square for tensors", D23, N77, N78, N79 (+ ledger row), N80, N81, **E7** |
| §4.9 head | pointers to N82–N84 (C04) and N103–N105 (C05) |
| C11 | N85, N86, N87, D24, N88, N89, D25, N90, N91, N92, N94, N95, N96, N102 |
| before C12 | R09 (`recap`) |
| C12 | D26, N93, N97, N98, N99, N100, N101, **E6** |
| before C13 | R10 (`recap` + `code`) |
| C13 | primers "order-of-magnitude scaling", "reduced gravity and buoyancy", N106, N107, N108, D27, N109, D28, N110, N111, N112, N113, N114, **E8** |
| before C14 | R11 (`recap`) |
| C14 | N115, N116, N117, primer "moving level set…", N118, D29, N119, N120, N121, N122, N123, N124, N125, N126, N127, N128, N129 |
| before C15 | R12 (`recap` + `code`) |
| C15 | N130, N132, primer "scaled variables…", N133, D30, N134, N135, N136, N137, N138, N139, N140, N141, N142, N143, N131, N144, N145, N146, N147, N148, N149, N150, N151, N152, N153, N154, N155, N156, **E9** |
| end | S01, S02 (`pointer`s), summary |

Check: 15 CORE blocks; 156 NOTE ids placed (N01–N156, N22–N25 and N69–N72 as grouped notes, N28 as a ⚠️ `md` with its id);
12 RECAPs each before its host block; 2 SKIP pointers; 30 derivations each in its curation block; 9 explainers embedded
(E1 in C04, E2 in C03, E3 in C07, E4 in C08, E5 in C09, E6 in C12, E7 in C10, E8 in C13, E9 in C15); B1 not embedded.

---

## Part B — explainer storyboards

Common to all ten: created with `tools/new_viz.py`; `<meta name="viz:chapter" content="ch04">`; tabs Walkthrough ·
Explore · Explain · Derivation · Equations · Code · Check; `autoplay: false` and `play: false` on every walkthrough step
that quotes numbers (ch02 lesson); every displayed number is computed by a JS function that mirrors a `ch04` callable and
is proved by `selftest()` parity rows (`py:` expressions use only `ch04.…`, `np.pi`, lists, dicts, floats and strings —
no builtins, no lambdas; index results down to one number; functions that take fields are called with a **preset name**,
Part C convention 1); Explain is "Explanation & interpretation" in numbered sections built with `Viz.work.step / line /
box / table / hint / interpret`, modelled on `forced_damped_vibrations.html` and `amplitude_phase_second_order_II_3.html`
(0 what the views show and what each colour means · 1…n every displayed quantity from the controls, "formula =
substituted = result — why", results boxed · a section per view hidden on phones or a hint · the values at the current
time (live) · **Reading the current setting**); derivation steps are copied from Part F (same `did` titles, same step
count; phones shorten *why* to its first sentence; plain-text *why* never contains raw TeX); every tour, Explain,
Derivation, notes, status, equation and quiz text that names a book equation **writes it out** next to its number
(convention 11; `tools/eq_refs.py` → 0). Colours as in convention 10. A view hidden on portrait phones never carries a
step's key number (repeat it in a visible title or readout). Mode-dependent controls are hidden per mode with a class
toggled by the chapter script (ch03 E7 pattern), never merely marked `optional`. Parallel builders use private scratch
subfolders. Each explainer fits 360×640 … 1920×1080 and the 1000×700 notebook frame with no scrolling (fit plan per
explainer).

### E1 · control_volume_budgets
- **Title:** "Weigh a force by counting what flows through a box" · **Summary:** "One control volume drawn on a real
  scene — a wake, a bore, a jet, a rocket, an expanding gas — with every face's mass and momentum flux as a bar: storage
  plus outflow equals the forces, and the residual stays at zero." · **CORE:** C01, C04 (also N03, N06, N26, N28, N29,
  N31, N32, N84) · **Reference:** `fid_formula_lab.html` (terms that add up to a total, stage chips) with
  `forced_damped_vibrations.html`'s explanation panel; ch03 E7 `reynolds_transport_cv` (swept band, budget waterfall,
  measured ◇).
- **meta:** `viz:order 1` · `viz:sections 4.2 4.4` · `viz:equations 4.1 4.5 4.17` · `viz:fluidpy ch04.cv_scenario
  ch04.interval_mass_budget ch04.wake_drag_per_span ch04.wake_side_outflow ch04.bore_speed ch04.jet_plate_force
  ch04.rocket_delta_v` · `viz:derivations D01 D05`.
- **Physics (JS ↔ Python):** `wakeU(y, p)` ↔ `ch04.gaussian_wake(y, U_inf, deficit, width)` · `wakeBudget(p)` → faces
  (inlet, outlet, top, bottom) with closed forms $\int_{-H/2}^{H/2}e^{-y^2/b^2}dy=b\sqrt\pi\,\mathrm{erf}(H/2b)$ and
  $\int e^{-2y^2/b^2}dy=b\sqrt{\pi/2}\,\mathrm{erf}(H/\sqrt2b)$ (`Viz.num.erf`) ↔ `ch04.cv_scenario("wake", …)` (quad) ·
  `boreBudget(p)` (x to the right, the wave travelling to −x at speed U into still water of depth h_in; behind it depth
  h_out moving at −u_out with u_out = U(1 − h_in/h_out) (mass in the wave frame); the box moves at −b, so b = U rides
  with the front; faces left/right (hydrostatic ½ρgh² and momentum fluxes), top (p_o, cancels), bed (no flux);
  storage_x = −ρh_out u_out(U − b), zero when b = U) ↔ `ch04.cv_scenario("bore", h_in=…, h_out=…, b=…, rho=…)` (the same
  convention, documented in its docstring) and `ch04.bore_speed` · `jetBudget(p)` ↔
  `ch04.cv_scenario("jet", …)`, `ch04.jet_plate_force(rho, V, A, theta)` · `rocketState(t, p)` (M = M₀ − ṁt,
  b = V_e ln(M₀/M) − gt) ↔ `ch04.cv_scenario("rocket", …)`, `ch04.rocket_delta_v` · `intervalBudget(x0, x1, t, ẋ0, ẋ1,
  a, ρ₀)` ↔ `ch04.interval_mass_budget(x0, x1, t, dx0dt, dx1dt, flow="expanding", a, rho0)`.
- **Modes** (`mode` chips, short labels "wake · bore · jet · rocket · balloon"): **wake** (Ex. 4.1, fixed box around a
  bar) · **bore** (Ex. 4.3, box moving at b) · **jet** (our extension: a water jet hitting a plate at angle θ, fixed box)
  · **rocket** (Ex. 4.4, box accelerating with the rocket) · **balloon** (the expanding 1-D flow of C01; box type chips
  fixed / sliding / material; mass budget only).
- **Views** (one clock `t` for bore, rocket and balloon; rows [1.3, 0.9]):
  1. `scene` "The control volume" (`equal: false`; wake: x ∈ [−0.5, 3] m, y ∈ [−1.1, 1.1] m; bore: x ∈ [−6, 6] m,
     z ∈ [0, 1.6] m; jet: a plate and nozzle; rocket: a column; balloon: x ∈ [0, 6] m) — the phenomenon (bar and its
     wake with U(y) profiles drawn at inlet and outlet; the bore's surface with parcels; the jet splitting on the plate;
     the rocket with its plume; the density strip and the material bracket), the **dashed CV** with outward-normal
     ticks, flux arrows per face (**blue in, orange out**, length ∝ |flux|), a **rose** arrow for the force on the fluid
     and (wake, jet) a dark arrow for the force on the body (Newton III). Pointer: drag the CV's top edge (wake: H), drag
     the CV (bore: sets b; balloon: slides the box), click a face → inspector.
  2. `bars` "Budget" — waterfall bars for the chosen budget (chips mass / x-momentum): storage d/dt∫ (purple), each
     face's outflux (orange +, blue −, hatched when negative), body force (grey), surface forces (rose), and the
     residual (black, ≈ 0) landing on a **measured ◇** (storage measured by differencing the box contents in time, ch03
     E7 pattern). Title carries the key number ("F_D/l = 3.652 N/m").
  3. `result` "What the budget gives" (`hidePortrait: true`) — wake: F_D/l vs H with the converged value dashed and the
     current H marked; bore: U vs h_out/h_in with the √(gh) ghost (grey) and the current ratio marked; jet: force vs θ
     (ρV²A sin θ); rocket: b(t) with the Tsiolkovsky ghost V_e ln(M₀/M); balloon: mass in each box vs t.
  Portrait: `scene` + `bars`; `result` hidden (its number is in the `bars` title and the readouts).
- **Controls (per mode, ≤ 5 visible):** `H` "Box height $H$" 0.05…2 m, step 0.01, default 1.0 (wake), help "must enclose
  the whole wake" · `deficit` "Wake deficit $\Delta$" 0…5 m/s, default 2 (wake) · `b` "Box speed $b$" −1…6 m/s, step
  0.01, default 0 (bore), help "the box rides with the wave when b equals its speed" · `ratio` "Depth ratio
  $h_{out}/h_{in}$" 1.01…2.0, default 1.10 (bore) · `theta` "Plate angle $\theta$" 10…90°, default 90 (jet) · `Ve`
  "Exhaust speed $V_e$" 200…3000 m/s, default 800 (rocket) · `box` chips "fixed · sliding · material" (balloon) ·
  `budget` chips "mass · momentum" (all but balloon). Fixed internally: U∞ = 10 m/s, b_w = 0.1 m, ρ_air = 1.2 kg/m³,
  h_in = 1 m, water ρ = 1000, jet V = 10 m/s, A = 10⁻⁴ m², rocket M₀ = 2 kg, ṁ = 0.1 kg/s, burn 10 s; balloon a = 1 s⁻¹,
  ρ₀ = 1 kg/m³. Transport `t`: bore 0 → 4 s (the front crosses the box), rocket 0 → 10 s, balloon 0 → 1.5 s; rate 1,
  `end: 'hold'`, hold 2 s; end-of-run `Viz.card`: bore "The box rode with the front: every frame the same budget —
  steady." / rocket "Burn-out at M = 1 kg: b = 554.5 − 98.1 = 456.4 m/s." / balloon "Mass in the material box: 1.000
  kg/m² from start to end."
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **terms** (waterfall, click a bar to isolate
  it) · **modes** (5) · **presets** "wake, tall box" {mode: 'wake', H: 1.0} · "wake, box too short" {H: 0.15} · "bore,
  fixed box" {mode: 'bore', b: 0} · "bore, riding box" {b: 3.366} · "jet on a plate" {mode: 'jet', theta: 90} ·
  "rocket" {mode: 'rocket', t: 0} · "balloon (b = u)" {mode: 'balloon', box: 'material'} · **status** (wake: H ≥ 3b_w
  "✅ budget closes · F_D/l = 3.652 N/m", H < 3b_w "⚠️ box too short: it cuts the wake — F_D/l = … N/m"; bore:
  |b − U| < 0.01 "✅ riding with the wave: steady budget, U = 3.366 m/s", else "⏱ box not riding: storage = … N/m"; jet
  "🚿 force on the plate ρV²A sin θ = … N"; rocket "🚀 thrust V_eṁ = 80 N, b = … m/s"; balloon material "🎈 material box:
  no flux, mass constant — (4.1)", fixed "📦 fixed box: loses mass at … kg/(m² s)") · **inspector** (click a face: "top
  face: area H… ; (u − b)·n = …; mass flux ρ∫(u − b)·n dA = 1.2 × 0.1772 = 0.2127 kg/(m s); x-momentum flux = U∞ × that =
  2.127 N/m") · **transport**.
- **Readouts:** "Result" (drag, wave speed, force or rocket speed with unit) · "Residual" · "Side leakage" (wake) ·
  "Storage".
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**The control volume** is the dashed box drawn on the flow. Arrows on its faces show fluid
     crossing: <b class="c-blue">blue</b> entering, <b class="c-orange">orange</b> leaving; the <b class="c-rose">rose</b>
     arrow is the force on the fluid inside; the dark arrow is its reaction, the force on the body. **Budget** stacks
     storage (purple), outflow per face and the forces; the black residual must be zero. **What the budget gives**
     (hidden on phones) plots the answer against the box size or the time."
  1. *The rule* — the box's momentum budget (4.17),
     $\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV+\int_{A^*}\rho\mathbf u(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA$,
     and its mass version (4.5), $\frac{d}{dt}\int_{V^*}\rho\,dV+\int_{A^*}\rho(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=0$, with the
     current box speed: "b = `live b` m/s, so the relative velocity on each face is u − b."
  2. *Face by face* (`Viz.work.table` with the current row lit): face · area · (u − b)·n · mass flux · x-momentum flux —
     wake numbers: inlet −24.000 / −240.000, outlet +23.575 / +232.094, top +0.213 / +2.127, bottom +0.213 / +2.127
     (H = 2 m); the line "mass: out − in = 23.575 + 0.425 − 24.000 = 0 ✓ (steady, fixed box, storage 0)".
  3. *Forces* — "pressure is p∞ on every face, so its net force is zero (a constant pressure on a closed surface pushes
     equally in all directions — Gauss); gravity acts vertically and does not enter the x-budget. The only x-force on the
     fluid is the bar's reaction, −F_D/l."
  4. *The result* (boxed, per mode) — wake: $F_D/l=\rho\int U(U_\infty-U)dy=1.2\times(3.5449-0.5013)=\mathbf{3.652}$ N/m;
     bore: $U=\sqrt{gh_{out}(h_{in}+h_{out})/2h_{in}}=\sqrt{9.81\times1.1\times2.1/2}=\mathbf{3.366}$ m/s (vs √(gh) = 3.132);
     jet: F = ρV²A sin θ = 1000 × 100 × 10⁻⁴ × 1 = **10.0** N; rocket: $M\,d^2z_R/dt^2=-V_e\,dM/dt-Mg+F_S$, speed at t:
     **…** m/s; balloon: storage = −(flux_out − flux_in) = **…** kg/(m² s).
  5. *Newton's third law* — "the budget gives the force on the fluid; the body feels the opposite. Forgetting the minus
     makes drag negative."
  6. *The result view* (hint on phones: "turn the phone to see the result curve") — wake: "F_D/l stops changing once
     H ≳ 3b_w = 0.3 m; below that the box cuts the wake"; bore: "U → √(gh) as the ratio → 1 (small waves travel at the
     shallow-water speed)".
  7. *At the current time* — t = `live t` s, storage `live storage`, residual `live res`.
  8. *Reading the current setting* — wake tall: "The box encloses the whole wake: its drag reading is final. The side
     leakage is not a mistake: less fluid leaves through the outlet than enters, so the rest leaves through the top and
     bottom, carrying momentum U∞ with it." · wake short: "The box cuts the wake: part of the deficit crosses the top and
     bottom where we assumed U∞; the reading is too low." · bore not riding: "The front moves through the box, so the
     momentum inside changes: storage ≠ 0. Drag b until it vanishes — that speed is the wave's." · bore riding: "In the
     wave's frame everything is steady: water arrives 1 m deep and fast, leaves 1.1 m deep and slower; the depth
     difference's hydrostatic push is what slows it." · jet: "The plate turns the jet; the momentum it removes each
     second is the force." · rocket: "The box accelerates with the rocket; only the exhaust crosses its surface, carrying
     momentum V_e per kg backwards — that is thrust." · balloon material: "b = u on both ends: nothing crosses, mass is
     constant — (4.1), $\frac{d}{dt}\int_{V(t)}\rho\,dV=0$." · balloon fixed: "The fixed box loses mass through its right face
     faster than it gains on the left: storage −1 kg/(m² s) at t = 0."
- **Derivation tab:**
  - **D01** (9 steps) `view: 'scene'`; goal `set {mode: 'balloon', box: 'material', t: 0}`; step 1 (material volume)
    highlights the bracket; step 3 `set {box: 'fixed'}` (the arbitrary CV appears beside the material one); steps 4–6
    (the coincidence) `set {t: 0.0}` with both outlines drawn on top of each other and `watch: "the two outlines coincide
    now and separate as soon as time runs"`; step 7 **live** "$-\int\rho\mathbf u\cdot\mathbf n+\int\rho\mathbf b\cdot\mathbf n$ = −(2 − 1)
    + (0 − 0) = −1 kg/(m² s)"; step 8 `set {box: 'material'}` (b = u check: bars vanish); step 9 `set {box: 'fixed'}`
    (b = 0 check); interpret `s => "With your box (${boxName}), storage ${storage} + net outflux ${out} = 0."`.
  - **D05** (8 steps) `view: 'scene'`; goal `set {mode: 'wake', H: 1.0, budget: 'momentum'}`; step 3 shows the ⚠️ about
    the printed "= 0" of (4.15) in *plain*; step 4 lights the four bar groups one after another (the four coincidence
    equalities (4.16a–d)); step 7 **live** "out − in = 232.094 + 4.254 − 240 = −3.652 N/m = −F_D/l"; interpret "the force
    on the body is minus the force on the fluid: F_D/l = 3.652 N/m".
- **Code:**
  ```python
  sc = ch04.cv_scenario("{{mode}}", **params)   # {{params}} — the four modes of E1
  for f in sc["faces"]:                          # outward n on every face; out > 0, in < 0
      print(f["name"], f["mass_flux"], f["momentum_flux_x"])
  print(sc["residual_mass"], sc["residual_momentum"])   # {{rm}}, {{rp}}  (≈ 0)
  print(sc["result"], sc["result_label"])        # {{result}} {{label}}
  FD = ch04.wake_drag_per_span("gaussian", 10.0, 1.2, {{H}})   # (Ex. 4.1) = {{FD}} N/m
  U = ch04.bore_speed(1.0, {{hout}})             # (Ex. 4.3) = {{U}} m/s
  b = ch04.interval_mass_budget(1.0, 2.0, 0.0, dx0dt={{d0}}, dx1dt={{d1}})
  print(b["storage"], b["net_outflux"])          # {{st}}, {{no}} kg/(m² s)
  ```
- **Walkthrough (7 steps):** 1. "Weigh without touching" — "A bar sits in a wind tunnel. You may only measure the air
  behind it. Can you find the force on it?" `set {mode: 'wake', H: 1.0}`, highlight `view:scene` · 2. "The balloon keeps
  its mass" — "First, mass. A material box moves with the fluid (b = u): nothing crosses, so its mass is constant —
  $\frac{d}{dt}\int_{V(t)}\rho\,dV=0$ (4.1)." `set {mode: 'balloon', box: 'material'}`, `play: true` · 3. "Any box: count
  relative to the walls" — "A fixed box on the same flow loses mass through its right face. Flux counts u − b:
  $\frac{d}{dt}\int_{V^*}\rho\,dV+\int_{A^*}\rho(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=0$ (4.5)." `set {box: 'fixed', t: 0}`, `play:
  false`, `derive: {id: 'D01', step: 7}` · 4. "Momentum in and out" — "Back to the wake. Less momentum leaves than
  enters; the difference, side leak included, is the force on the fluid: d/dt∫ρu dV + ∮ρu(u − b)·n dA = ∫ρg dV + ∮f dA
  (4.17)." `set {mode: 'wake', budget:
  'momentum'}`, `terms: true`, `eq: 'mom'` · 5. "Whose force?" — "The budget gives −F_D/l on the fluid; the bar feels
  +F_D/l = 3.65 N/m. Newton's third law." `readouts: ['result']`, highlight `readout:result`, `code: {id: 'e1', lines: [5,
  5]}` · 6. "A box that rides" — "A bore moves into still water. Drag the box speed b until storage vanishes: you have
  found the wave's speed." `set {mode: 'bore', b: 0}`, `controls: ['b']`, `inspect: true` · 7. "Your turn" — "Predict
  first: if the box is only 0.15 m tall, will the drag reading go up or down? Then try it." `set {mode: 'wake', H: 1.0}`,
  `controls: ['H']`.
- **Equations:** `mass1` "Material volume" ref 'Eq. (4.1)' `\frac{d}{dt}\int_{V(t)}\rho\,dV=0` live "mass = 1.000 kg/m²" ·
  `mass` "Any control volume" ref 'Eq. (4.5)' `\frac{d}{dt}\int_{V^*}\rho\,dV+\int_{A^*}\rho(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=0`
  live "storage −1.000 + outflux 1.000 = 0" · `mom` "Momentum" ref 'Eq. (4.17)'
  `\begin{aligned}&\tfrac{d}{dt}\textstyle\int_{V^*}\rho\mathbf u\,dV+\int_{A^*}\rho\mathbf u(\mathbf u-\mathbf b)\cdot\mathbf n\,dA\\&=\textstyle\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA\end{aligned}`
  live "x: 0 + (−3.652) = 0 + (−3.652) N/m" · `wake` "Drag from the wake" ref 'Ex. 4.1'
  `F_D/l=\rho\int U(U_\infty-U)\,dy` live "= 3.652 N/m" · `bore` "Bore speed" ref 'Ex. 4.3'
  `U=\sqrt{\frac{gh_{out}}{2h_{in}}(h_{in}+h_{out})}` live "= 3.366 m/s". Symbols: ρ, u, b, n, V*, A*, f, g with units.
- **Check yourself:** (1) "Double the wake deficit Δ. Does the drag double?" — "More than double would be wrong and so
  would exactly double: F_D/l = ρ(U∞Δb√π − Δ²b√(π/2)) — the Δ² term grows faster, so the drag rises by less than ×2
  (3.65 → 6.10 N/m)." `set {mode: 'wake', deficit: 4}` · (2) "In the bore mode with b = 0, is the storage positive or
  negative?" — "Negative (x to the right): the front moves left through the box, replacing still water with deeper water
  moving left, so the box's x-momentum decreases at ρh_out u_out U ≈ 1133 N/m; at b = U it vanishes." `set {mode: 'bore', b: 0}` · (3)
  "Why does fluid leave through the top and bottom of the wake box?" — "The outlet carries less fluid than the inlet
  (the deficit); mass conservation sends the difference out sideways, 0.354 m²/s." `set {mode: 'wake'}` · (4) "What box
  speed makes every flux bar vanish in the balloon mode?" — "The material one: each end moving at the local u (b = u)."
  `set {mode: 'balloon'}`.
- **Selftest parity rows:** `{name: 'wake drag H=2', js: wakeBudget({U: 10, D: 2, b: 0.1, rho: 1.2, H: 2}).FD, py:
  'ch04.wake_drag_per_span("gaussian", 10.0, 1.2, 2.0)', rtol: 1e-8}` · `{name: 'side outflow', js: wakeBudget({U: 10, D:
  2, b: 0.1, rho: 1.2, H: 2}).side, py: 'ch04.wake_side_outflow("gaussian", 10.0, 2.0)', rtol: 1e-8}` · `{name: 'short box
  drag H=0.15', js: wakeBudget({U: 10, D: 2, b: 0.1, rho: 1.2, H: 0.15}).FD, py: 'ch04.wake_drag_per_span("gaussian",
  10.0, 1.2, 0.15)', rtol: 1e-8}` · `{name: 'bore speed', js: boreSpeed(1.0, 1.1, 9.81), py: 'ch04.bore_speed(1.0, 1.1)',
  rtol: 1e-12}` · `{name: 'bore storage b=0', js: boreBudget({hin: 1, hout: 1.1, b: 0, rho: 1000}).storage, py:
  'ch04.cv_scenario("bore", h_in=1.0, h_out=1.1, b=0.0, rho=1000.0)["storage_x"]', rtol: 1e-8}` · `{name: 'jet force',
  js: jetForce(1000, 10, 1e-4, Math.PI/3), py: 'ch04.jet_plate_force(1000.0, 10.0, 1e-4, np.pi/3)', rtol: 1e-12}` ·
  `{name: 'rocket dv', js: rocketDV(2, 1, 800), py: 'ch04.rocket_delta_v(2.0, 1.0, 800.0)', rtol: 1e-12}` · `{name:
  'interval storage fixed', js: intervalBudget(1, 2, 0, 0, 0, 1, 1).storage, py: 'ch04.interval_mass_budget(1.0, 2.0,
  0.0)["storage"]', rtol: 1e-8}` · invariant `{name: 'material box no flux', js: intervalBudget(1, 2, 0, 1, 2, 1,
  1).net, expect: 0, atol: 1e-12}`.
- **Fit plan:** 360×640: `scene` (top) + `bars` (bottom), status one line ("✅ F_D/l 3.65 N/m"), mode chips wrap to two
  rows, presets hidden on short portrait screens (scoped CSS, Derivation tab exempt); Explore shows the ≤ 3 controls of
  the current mode. Landscape phone: `scene` + `bars` in one row, `result` dropped. Notebook 1000×700 and desktop: rows
  [1.3, 0.9], `result` beside `bars`. Explain §2's face table uses short headers ("face · A · (u−b)·n · ṁ · Mx").

### E2 · stream_function_spacing
- **Title:** "Can one number field hold a whole 2-D flow?" · **Summary:** "Seven classic flows drawn as ψ contours at equal
  steps over a speed map with tracers; drag a gate and its flux stays ψ₂ − ψ₁ however you bend it; click to read
  u = ∂ψ/∂y, v = −∂ψ/∂x." · **CORE:** C03 (also N11, N13, N16, N17, N18) · **Reference:** `stride_padding_playground.html`
  (classic presets, the formula with numbers plugged in, click a cell to move the window) with
  `forced_damped_vibrations.html`'s explanation panel.
- **meta:** `viz:order 2` · `viz:sections 4.3` · `viz:equations 4.11 4.12` · `viz:fluidpy ch04.streamfunction_preset
  ch04.velocity_preset ch04.velocity_from_streamfunction_2d ch04.flux_between_streamlines ch04.flux_along_path
  ch04.velocity_from_streamfunction_axisym` · `viz:derivations D04`.
- **Physics (JS ↔ Python):** `psi(name, x, y, p)` ↔ `ch04.streamfunction_preset(name, x, y, **p)` · `vel(name, x, y, p)`
  ↔ `ch04.velocity_preset(name, x, y, **p)` · `velFD(name, x, y, p, h)` (central differences of ψ) ↔
  `ch04.velocity_from_streamfunction_2d(name, x, y, h=h, **p)` · `gateFlux(p1, p2)` = ψ(p2) − ψ(p1) and
  `gateFluxNumeric(pts)` (16-point Gauss–Legendre per segment of u·n with n to the right of travel) ↔
  `ch04.flux_between_streamlines`, `ch04.flux_along_path` · axisymmetric `velAxi(R, z, U)` ↔
  `ch04.velocity_from_streamfunction_axisym("axisym_uniform", R, z, U=U)`.
- **Presets** (flow chips, short labels "uniform · stagnation · source+stream · cylinder · vortex · shear · axisym"):
  uniform (U, angle α) · stagnation ψ = kxy · source + stream (U = 1, m) · cylinder (U = 1, a = 1) · line vortex
  ψ = −(Γ/2π) ln r · shear ψ = ½γ̇y² · axisymmetric uniform stream ψ = ½UR² (x ↦ z, y ↦ R, with a note that equal Δψ
  steps now carry flux 2πΔψ and the contours crowd as 1/R).
- **Views** (rows [1.4, 0.8]):
  1. `flow` "ψ contours at equal steps" (`equal: true`, x, y ∈ [−3, 3] m) — speed |u| heat map (viridis, legend in the
     title), ψ contours (black, `n_levels` at equal Δψ), tracer particles (white dots moved by `Viz.num.rk4Step` on `vel`),
     stagnation points (✕) and singular points (●), the **gate** (a draggable polyline with 2–4 handles, orange) with its
     flux label, and two highlighted contours ψ₁, ψ₂ (teal) through the gate's ends with the band between them shaded.
     Pointer: drag gate handles; click elsewhere → inspector point.
  2. `gate` "Along the gate" — ψ(s) along the gate (teal) and the normal velocity u·n(s) (orange) with the area under it
     shaded; readout "∫u·n ds = ψ₂ − ψ₁ = 1.000 m²/s".
  3. `spacing` "Spacing = speed" (`hidePortrait: true`) — along a short probe line normal to the contours through the
     inspector point: contour spacing Δn (bars) and Δψ/Δn (dots) against the true |u| (line) — the rule u ≈ Δψ/Δn made
     visible.
- **Controls:** `flow` chips · `q` "Strength $q$" 0.1…3, default 1, with a meaning line (`qMeaning`: uniform → U [m/s],
  stagnation → k [1/s], source → m [m²/s], cylinder → U, vortex → Γ [m²/s], shear → $\dot\gamma$ [1/s], axisym → U) ·
  `levels` "Contours" 6…30, default 16 · `alpha` "Stream angle $\alpha$" 0…90° (uniform only; hidden otherwise) ·
  gate handles (pointer). Transport `t` 0 → 10 s, loop, for the tracers (no key numbers depend on t).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **presets** (7) · **inspector** (click a point:
  "$u=\partial\psi/\partial y\approx\frac{\psi(x,y+h)-\psi(x,y-h)}{2h}=\frac{2.0100-1.9900}{0.02}=1.000$ m/s, $v=-\partial\psi/\partial x\approx-\frac{…}{…}=-2.000$
  m/s") · **transport** (tracers) · **status** ("↔ fast where contours crowd · max speed … at …" / "⚠️ stagnation point:
  contours cross, u = 0" when the inspector sits within 0.05 m of one / "● singular point (source or vortex centre)").
- **Readouts:** "Gate flux ψ₂ − ψ₁" (m²/s) · "Numeric ∫u·n ds" · "u, v at the probe" · "Speed at the probe".
- **Explain:**
  0. *What the views show* — "Black lines are contours of ψ drawn at **equal steps** Δψ; the colours show the speed.
     The <b class="c-orange">orange</b> gate is a line you drag; the <b class="c-teal">teal</b> contours pass through
     its ends. **Along the gate** plots ψ and the velocity crossing the gate. **Spacing = speed** (hidden on phones)
     compares the gap between contours with the local speed."
  1. *The flow you chose* — "ψ = … (strength `live q` …): at the probe (x, y) = (…, …), ψ = …."
  2. *Velocity from ψ* — "$u=\partial\psi/\partial y=$ … m/s, $v=-\partial\psi/\partial x=$ … m/s (D04's result for χ = −z);
     speed … m/s." (boxed)
  3. *Flux through the gate* — "ψ at the ends: ψ₁ = …, ψ₂ = …; flux per metre of depth = ψ₂ − ψ₁ = **…** m²/s; the numeric
     integral along your gate: … m²/s (difference … ≤ 10⁻¹⁰). Bend the gate — only the end labels matter." (boxed)
  4. *Spacing and speed* — "Between neighbouring contours Δψ = … m²/s; here they are Δn = … m apart, so |u| ≈ Δψ/Δn = …
     m/s (true … m/s)." (hint on phones: "the Spacing view is hidden; its numbers are here.")
  5. *Continuity is automatic* — "(4.11), $\nabla\cdot(\rho\mathbf u)=0$, holds for any ψ because
     $\partial u/\partial x+\partial v/\partial y=\psi_{yx}-\psi_{xy}=0$ (mixed derivatives are equal)."
  6. *At the current time* — "tracers moved for t = `live t` s; each stays on its contour."
  7. *Reading the current setting* — stagnation: "Hyperbolas: fluid comes in along y, leaves along x; the contours cross
     at the stagnation point where the speed is zero." · cylinder: "The contours crowd over the top and bottom — speed 2U
     there — and spread at the front and back stagnation points." · vortex: "Circles packed tighter toward the centre:
     u = Γ/(2πr); the centre is singular." · source + stream: "The rose dividing streamline ψ = m/2 separates the
     source's fluid from the stream's: a half-body m/U wide far downstream." · shear: "Contours closer together as y
     grows: speed $\dot\gamma y$." · uniform: "Evenly spaced straight lines: uniform speed." · axisym: "Stream surfaces of
     revolution; equal Δψ steps crowd toward the axis because the flux between them is 2πΔψ through a ring of radius R."
- **Derivation tab:** **D04** (8 steps) `view: 'flow'`; goal `set {flow: 'uniform', q: 1, alpha: 0}`; step 2 (χ = −z)
  draws the planes z = const as the page; step 5 **live** "at the probe: ρu = ∂ψ/∂y = …"; step 7 `set {flow:
  'stagnation', q: 1}` and draws the gate from ψ = 1 to ψ = 2; step 8 **live** "ψ₂ − ψ₁ = 2 − 1 = 1.000 m²/s" `watch:
  "drag a handle: the flux readout does not change"`; interpret `s => "Here the flux between your two contours is … m²/s
  per metre of depth."`.
- **Code:**
  ```python
  name, p = "{{flow}}", {{params}}                    # preset and its strength
  u, v = ch04.velocity_from_streamfunction_2d(name, {{x}}, {{y}}, **p)
  print(u, v)                                         # {{u}}, {{v}} m/s  (u = ∂ψ/∂y, v = −∂ψ/∂x)
  p1, p2 = ({{x1}}, {{y1}}), ({{x2}}, {{y2}})             # the gate's ends
  print(ch04.flux_between_streamlines(name, p1, p2, **p))   # {{flux}} m²/s = ψ(p2) − ψ(p1)
  psi = lambda x, y: ch04.streamfunction_preset(name, x, y, **p)
  print(psi(*p2) - psi(*p1))                          # the same number from the labels
  ```
- **Walkthrough (6 steps):** 1. "One field, two jobs" — "Can a single number field tell both where the fluid goes and
  how fast?" `set {flow: 'cylinder'}`, highlight `view:flow` · 2. "Contours are streamlines" — "The tracers never leave
  their black line: ψ is constant along a streamline." `play: true` · 3. "Velocity is a slope" — "Click a point:
  $u=\partial\psi/\partial y$, $v=-\partial\psi/\partial x$ — read them off the neighbouring contours." `set {flow:
  'stagnation'}`, `play: false`, `inspect: true`, `derive: {id: 'D04', step: 5}` · 4. "Flux is a difference" — "Drag the
  gate between ψ = 1 and ψ = 2: 1 m²/s crosses it, whatever its shape." `readouts: ['flux', 'num']`, `derive: {id:
  'D04', step: 8}` · 5. "Crowded means fast" — "Equal flux in every band, so a narrow band must be fast: |u| ≈ Δψ/Δn."
  `set {flow: 'cylinder'}`, highlight `view:spacing`, `readouts: ['speed']` · 6. "Your turn" — "Predict: in the line
  vortex, where are the contours closest? Check with the probe." `set {flow: 'vortex'}`.
- **Equations:** `cont` "Steady continuity" ref 'Eq. (4.11)' `\nabla\cdot(\rho\mathbf u)=0` · `twopsi` "Two stream
  functions" ref 'Eq. (4.12)' `\rho\mathbf u=\nabla\chi\times\nabla\psi` · `psi2d` "2-D stream function" ref '§4.3'
  `\rho u=\partial\psi/\partial y,\quad\rho v=-\partial\psi/\partial x` live "u = 1.000, v = −2.000" · `flux` "Flux between
  streamlines" ref 'Ex. 4.8' `\int_1^2\mathbf u\cdot\mathbf n\,ds=\psi_2-\psi_1` live "= 1.000 m²/s" · `axi` "Axisymmetric"
  ref '§4.3' `\rho u_R=-\frac1R\frac{\partial\psi}{\partial z},\quad\rho u_z=\frac1R\frac{\partial\psi}{\partial R}`.
- **Check yourself:** (1) "In the stagnation flow with k = 1 s⁻¹, what flux crosses a gate from (1, 1) to (2, 1)?" —
  "ψ = kxy: 2 − 1 = 1 m²/s." `set {flow: 'stagnation', q: 1}` · (2) "Double the cylinder's U. What happens to the flux
  between two given points?" — "It doubles: ψ ∝ U." `set {flow: 'cylinder', q: 2}` · (3) "Why do contours cross at a
  stagnation point?" — "There ∇ψ = 0 (u = v = 0): ψ has a saddle, and a saddle's level set is two crossing lines." ·
  (4) "Shear flow ψ = ½γ̇y²: where is the fluid at rest?" — "On y = 0, where the contours are infinitely far apart."
  `set {flow: 'shear'}`.
- **Selftest parity rows:** `{name: 'stagnation u', js: velFD('stagnation', 1, 2, {k: 1}, 1e-5)[0], py:
  'ch04.velocity_from_streamfunction_2d("stagnation", 1.0, 2.0, k=1.0)[0]', rtol: 1e-8}` · `{name: 'cylinder v', js:
  vel('cylinder', 2, 0.5, {U: 1, a: 1})[1], py: 'ch04.velocity_preset("cylinder", 2.0, 0.5, U=1.0, a=1.0)[1]', rtol:
  1e-12}` · `{name: 'gate flux', js: gateFluxNumeric([[1, 1], [1.2, 0.9], [1.4142135623730951, 1.4142135623730951]],
  'stagnation', {k: 1}), py: 'ch04.flux_along_path("stagnation", [[1.0, 1.2, 1.4142135623730951], [1.0, 0.9,
  1.4142135623730951]], k=1.0)', rtol: 1e-9}` · `{name: 'vortex psi', js: psi('vortex', 0.5, 0.5, {Gamma: 6.283185307179586}),
  py: 'ch04.streamfunction_preset("vortex", 0.5, 0.5, Gamma=6.283185307179586)', rtol: 1e-12}` · `{name: 'axisym uz',
  js: velAxi(0.5, 0.3, 2)[1], py: 'ch04.velocity_from_streamfunction_axisym("axisym_uniform", 0.5, 0.3, U=2.0)[1]', rtol:
  1e-8}` · invariant `{name: 'flux = label difference', js: gateFluxNumeric([[1, 1], [2, 1]], 'stagnation', {k: 1}) - 1,
  expect: 0, atol: 1e-10}`.
- **Fit plan:** 360×640: `flow` (square, top) + `gate` (bottom, short); `spacing` hidden (its numbers in Explain §4 and
  the "Speed" readout); flow chips wrap (≤ 7 short labels). Landscape phone: `flow` + `gate` side by side. Desktop:
  `flow` large left, `gate` and `spacing` stacked right.

### E3 · newtonian_stress_lab
- **Title:** "How does a fluid decide its stress?" · **Summary:** "Set a velocity gradient and watch G split into S and
  R, the Newtonian law turn S into τ, and the traction on a plane you rotate swing through its normal and shear values;
  a second mode shrinks a cube with unequal shear stresses and watches it spin up." · **CORE:** C07 (also N38, N40, N41,
  R05, N42–N47, N49–N51; C06's traction) · **Reference:** `amplitude_phase_second_order_II_3.html` (several windows on one
  state, crosshair readouts, a numbered live derivation in the explanation); links to ch02 `cauchy_traction_principal_axes`.
- **meta:** `viz:order 3` · `viz:sections 4.5` · `viz:equations 4.25 4.26 4.27 4.28 4.29 4.31 4.33 4.34 4.37` ·
  `viz:fluidpy ch04.newtonian_stress ch04.stress_on_plane ch04.mean_pressure ch04.isotropic_fourth_order
  ch04.linear_stress ch04.cube_spin_acceleration ch04.bulk_viscosity` · `viz:derivations D08 D09 D10`.
- **Physics (JS ↔ Python):** `strainRate(G)`, `rotation(G)` (R = G − Gᵀ, the book's) ↔ `tensors.strain_rate_tensor`,
  `tensors.rotation_tensor` · `stress(G, p, mu, muv)` = −pδ + 2μ(S − ⅓S_mmδ) + μ_vS_mmδ ↔ `ch04.newtonian_stress(G, p, mu,
  mu_v=muv)` · `onPlane(G, p, mu, muv, th)` → (σ_n, τ_s) with f_j = n_iτ_ij ↔ `ch04.stress_on_plane(G, p, mu, muv, th)` ·
  `meanP(tau)` ↔ `ch04.mean_pressure` · `K4(lam, mu, gam)` (81 entries) and `contractK(K, S)` ↔
  `ch04.isotropic_fourth_order`, `ch04.linear_stress` · `cubeSpin(t12, t21, rho, h)` ↔ `ch04.cube_spin_acceleration` ·
  `deformSquare(G, t)` by the closed-form 2×2 `expm2` (ch02 E3, promotion candidate) for the element picture.
- **Modes:** **stress lab** · **spinning cube** (D08).
- **Views** (rows [1.2, 0.9]):
  1. `element` "The fluid element" (`equal: true`) — stress lab: a unit square deforming under the 2-D part of G for a
     short time (outline at t and ghost at t = 0), the S eigen-directions (teal ±, stretching blue / compressing rose
     arrows), a **cutting plane** through the centre at angle θ (drag to rotate) with its normal n and the traction f
     (black) split into normal σ_n n (orange) and shear τ_s t (rose); cube mode: a square of side h with the four shear
     arrows τ₁₂ (on the x-faces) and τ₂₁ (on the y-faces), its angular acceleration arrow and the log readout.
  2. `mohr` "Stress on the plane vs its angle" — σ_n(θ) (orange) and τ_s(θ) (rose) over θ ∈ [0°, 180°] with the current θ
     marked, the principal directions (τ_s = 0) as ticks, and a small inset of the (σ_n, τ_s) circle (Mohr, the ch02
     gloss); cube mode: α(h) on log–log axes with slope −2 and the current h.
  3. `matrices` "G → S → τ" (`hidePortrait: true`) — three 3×3 grids: G (grey), S (teal) and R (orange, small), τ with each
     entry coloured by its source (−pδ orange, 2μ(S − ⅓S_mmδ) rose, μ_vS_mmδ purple); clicking an entry opens the
     inspector.
- **Controls:** `flow` chips "shear · extension · rotation · expansion · custom" · `k` "Rate $k$" 0…20 s⁻¹, default 10,
  meaning line (`kMeaning`: shear → $\dot\gamma$ = du/dy; extension → s = ∂u/∂x = −∂v/∂y; rotation → ω₀ (solid body);
  expansion → e = ∂u/∂x = ∂v/∂y = ∂w/∂z; custom → scale of the four sliders) · `mu` "Viscosity $\mu$" 10⁻⁴…1 Pa s (log),
  default 10⁻³ · `muv` "Bulk viscosity $\mu_v$" 0…5×10⁻³ Pa s, default 0 (Stokes), help "resists expansion only" · `theta`
  "Plane angle $\theta$" 0…180°, default 30 (also drag) · `p` "Pressure $p$" (optional, gauge, default 0) · custom mode
  adds G₁₁, G₁₂, G₂₁, G₂₂ (hidden in other modes) · cube mode: `imb` "Imbalance τ₁₂ − τ₂₁" 0.01…1 Pa (log), `h` "Cube
  side h" 10⁻⁴…10⁻¹ m (log). Transport (cube mode only): param `shrink` = −log₁₀h from 1 to 4, rate 0.5/s, end hold —
  shrinks the cube geometrically (ch02 E4 pattern); end card "h = 0.1 mm: α = 6×10⁵ rad/s² — impossible, so τ₁₂ = τ₂₁".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **presets** (5 flows + "Stokes vs bulk" {flow:
  'expansion', muv: 0.003}) · **inspector** (click a τ entry: "τ₁₁ = −p + 2μ(S₁₁ − ⅓S_mm) + μ_vS_mm = 0 + 2×10⁻³×(1 − 0)
  + 0 = 2.00×10⁻³ Pa") · **status** ("🔄 rigid rotation: G ≠ 0 but S = 0 — no viscous stress" / "↔ shear: τ₁₂ = μ$\dot\gamma$ =
  0.010 Pa" / "⇄ extension: τ₁₁ + p = 2μs" / "🫧 expansion: p̄ ≠ p unless μ_v = 0 (now p − p̄ = … Pa)" / cube "🌀 α = …
  rad/s² at h = …") · **modes** (stress lab / spinning cube).
- **Readouts:** "σ_n on the plane" · "τ_s on the plane" · "p − p̄" · "τ₁₂" · cube: "α".
- **Explain:**
  0. *What the views show* — "**The fluid element** deforms under your velocity gradient; the
     <b class="c-teal">teal</b> axes are the principal directions of S. The black arrow is the traction on the cutting
     plane; its <b class="c-orange">orange</b> part pulls or pushes normal to the plane, its <b class="c-rose">rose</b>
     part shears along it. **Stress on the plane** shows both as the plane turns. **G → S → τ** (hidden on phones)
     colours each stress entry by where it comes from."
  1. *Your velocity gradient* — "G = [[…]] s⁻¹ (`live`): S = ½(G + Gᵀ) = [[…]], R = G − Gᵀ = [[…]]; S_mm = ∇·u = …
     s⁻¹." (why: "only S enters the stress — a rigid spin R deforms nothing")
  2. *The Newtonian law* — "(4.37), $\tau_{ij}=-p\delta_{ij}+2\mu\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)+\mu_vS_{mm}\delta_{ij}$, entry by
     entry: τ₁₁ = … , τ₁₂ = 2μS₁₂ = 2 × … × … = **…** Pa, …" (boxed τ matrix)
  3. *Traction on your plane* — "n = (cos θ, sin θ) = (…, …); $f_j=n_i\tau_{ij}$ = (…, …) Pa; σ_n = f·n = **…** Pa,
     τ_s = f·t = **…** Pa." (boxed)
  4. *Principal directions* — "τ_s = 0 at θ = … and … (the eigenvectors of S); the largest normal stress there is … Pa."
  5. *Pressures* — "mean pressure (4.33), $\bar p=-\tfrac13\tau_{ii}=$ … Pa; thermodynamic p = … Pa; their difference (4.34),
     $p-\bar p=\big(\tfrac23\mu+\lambda\big)\nabla\cdot\mathbf u=\mu_v\nabla\cdot\mathbf u=$ … Pa." (zero unless expanding with μ_v > 0)
  6. *The matrices view* (hint on phones) — "the source colours show which term made each entry: for shear only rose
     off-diagonal entries; for expansion with μ_v, purple on the diagonal."
  7. *Cube mode* (when active) — "net torque (τ₁₂ − τ₂₁)h³ = …; moment of inertia ρh⁵/6 = …; α = 6(τ₁₂ − τ₂₁)/(ρh²) = **…**
     rad/s² — ×100 for every ×10 smaller h."
  8. *Reading the current setting* — shear: "Shear flow is extension along 45° and compression along 135°; the
     traction on a 45° plane is a pure pull of μ$\dot\gamma$." · extension: "Stretching along x pulls the x-faces (τ₁₁ + p > 0) and
     pushes the y-faces." · rotation: "The element turns without deforming: every viscous stress is zero — viscosity
     resists deformation, not motion." · expansion + μ_v = 0: "Pure swelling: with Stokes' assumption no viscous stress
     at all, p̄ = p." · expansion + μ_v > 0: "Now every face feels an extra pull μ_v∇·u: the mean normal stress differs
     from −p." · cube: "The imbalance cannot survive at a point: the smaller the cube, the faster it would spin."
- **Derivation tab:**
  - **D08** (8 steps) `view: 'element'`; goal `set {mode: 'cube', imb: 1, h: 0.01}`; step 3 draws the moment arms h/2;
    step 6 **live** "α = 6 × 1/(1000 × 10⁻⁴) = 60 rad/s²"; step 7 `set {h: 0.001}` `watch: "α jumps ×100"`; interpret "the
    only finite answer as h → 0 is τ₁₂ = τ₂₁".
  - **D09** (13 steps, ★★★) `view: 'matrices'` (phones: `mohr`); goal `set {mode: 'stress', flow: 'shear', k: 10}`; steps
    1–2 light G then S; step 4 (the three δδ terms) lights three colours in the τ grid one by one; steps 5–7 (δ
    substitutions) each light the entries they produce, step 7 **live** "δ_inδ_jmS_mn = S_ji = S_ij: τ₁₂ gets (μ + γ)S₁₂ =
    …"; step 9 (γ = μ) `watch: "only μ + γ acts on a symmetric S"`; step 11 `set {flow: 'expansion'}` (the λS_mm term
    appears on the diagonal); step 12 rotates the plane (the rotation-invariance check); step 13 **live** the full τ with
    your numbers; interpret `s => "With μ = …, λ = …: τ₁₂ = … Pa, τ₁₁ + p = … Pa."`.
  - **D10** (4 steps) `view: 'matrices'`; goal `set {flow: 'expansion', muv: 0.003}`; step 2 (add and subtract ⅔μS_mmδ)
    splits the diagonal colours into rose and purple; step 4 **live** "μ_v = λ + ⅔μ = … Pa s"; interpret "p − p̄ = μ_v∇·u =
    … Pa".
- **Code:**
  ```python
  G = np.array({{G}})                              # velocity gradient [1/s]
  tau = ch04.newtonian_stress(G, p={{p}}, mu={{mu}}, mu_v={{muv}})   # (4.37)
  print(tau)                                       # τ12 = {{t12}} Pa
  sn, ts = ch04.stress_on_plane(G, {{p}}, {{mu}}, {{muv}}, {{th}})  # traction f_j = n_i τ_ij
  print(sn, ts)                                    # σn = {{sn}}, τs = {{ts}} Pa
  print({{p}} - ch04.mean_pressure(tau))           # p − p̄ = {{dp}} Pa
  print(ch04.cube_spin_acceleration({{imb}}, 0.0, 1000.0, {{h}}))   # {{alpha}} rad/s²
  ```
- **Walkthrough (7 steps):** 1. "Same fluid, three motions" — "Slide it, stretch it, spin it. How does honey decide what
  stress to carry?" `set {mode: 'stress', flow: 'shear', k: 10}`, highlight `view:element` · 2. "Only deformation
  counts" — "Choose rotation: G is not zero, yet every stress is. Viscous stress depends on S, not on the spin R."
  `set {flow: 'rotation'}` · 3. "Linear and isotropic" — "The simplest law: σ = KS with K the same in every direction —
  three δδ terms (4.29), $K_{ijmn}=\lambda\delta_{ij}\delta_{mn}+\mu\delta_{im}\delta_{jn}+\gamma\delta_{in}\delta_{jm}$." `set {flow:
  'shear'}`, `derive: {id: 'D09', step: 4}` · 4. "Two constants survive" — "On a symmetric S the μ and γ terms coincide:
  $\tau_{ij}=-p\delta_{ij}+2\mu S_{ij}+\lambda S_{mm}\delta_{ij}$ (4.31). In shear, τ₁₂ = μ$\dot\gamma$ = 0.010 Pa." `eq: 'newt'`,
  `readouts: ['t12']` · 5. "Turn the plane" — "Drag the plane: the traction swings between pull and shear. At 45° the
  shear vanishes — a principal direction." `set {theta: 45}`, `controls: ['theta']` · 6. "Bulk viscosity" — "Expand
  the fluid and raise μ_v: the mean normal stress leaves −p. Only expansion feels μ_v." `set {flow: 'expansion', muv:
  0.003}`, `readouts: ['dp']` · 7. "Why τ₁₂ = τ₂₁" — "Predict: halve the cube, how much faster does it spin? Then press ▶."
  `set {mode: 'cube', imb: 1, h: 0.01}`, `derive: {id: 'D08', step: 6}`.
- **Equations:** `sym` "Symmetry" ref 'Eq. (4.25)' `\tau_{ij}=\tau_{ji}` · `split` "Pressure + viscous" ref 'Eq. (4.27)'
  `\tau_{ij}=-p\delta_{ij}+\sigma_{ij}` · `K` "Isotropic K" ref 'Eq. (4.29)'
  `K_{ijmn}=\lambda\delta_{ij}\delta_{mn}+\mu\delta_{im}\delta_{jn}+\gamma\delta_{in}\delta_{jm}` · `newt` "Newtonian stress" ref
  'Eq. (4.31)' `\tau_{ij}=-p\delta_{ij}+2\mu S_{ij}+\lambda S_{mm}\delta_{ij}` live "τ₁₂ = 2 × 0.001 × 5 = 0.010 Pa" · `bulk`
  "Bulk-viscosity form" ref 'Eq. (4.37)'
  `\begin{aligned}\tau_{ij}=&-p\delta_{ij}+2\mu\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)\\&+\mu_vS_{mm}\delta_{ij}\end{aligned}` ·
  `pbar` "Mean pressure" ref 'Eq. (4.33)–(4.34)' `\bar p=-\tfrac13\tau_{ii},\quad p-\bar p=\big(\tfrac23\mu+\lambda\big)\nabla\cdot\mathbf u`.
- **Check yourself:** (1) "In shear flow with $\dot\gamma$ = 10 s⁻¹ and μ = 10⁻³ Pa s, what normal stress acts on a 45° plane?" —
  "σ_n = μ$\dot\gamma$ = 0.010 Pa (a pull): shear is extension along 45°." `set {flow: 'shear', k: 10, theta: 45}` · (2) "Why is
  every stress zero for 'rotation' although the fluid moves?" — "A rigid rotation has S = 0; viscous stress depends only
  on S (Galilean invariance and 'only deformation costs')." `set {flow: 'rotation'}` · (3) "Stokes' assumption on: does
  μ_v matter for shear?" — "No: S_mm = 0 in shear; μ_v multiplies ∇·u." `set {flow: 'shear', muv: 0.003}` · (4) "A cube
  with a 0.01 Pa imbalance and h = 1 mm spins up at what rate?" — "6 × 0.01/(1000 × 10⁻⁶) = 60 rad/s²; at 0.1 mm, 6000."
  `set {mode: 'cube', imb: 0.01, h: 0.001}`.
- **Selftest parity rows:** `{name: 'tau12 shear', js: stress(shearG(10), 0, 1e-3, 0)[0][1], py:
  'ch04.newtonian_stress([[0.0, 10.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], p=0.0, mu=1e-3, mu_v=0.0)[0][1]', rtol:
  1e-12}` · `{name: 'sigma_n 30deg', js: onPlane(shearG(10), 0, 1e-3, 0, Math.PI/6)[0], py:
  'ch04.stress_on_plane([[0.0, 10.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], 0.0, 1e-3, 0.0, np.pi/6)[0]', rtol: 1e-12}`
  · `{name: 'tau_s 30deg', js: onPlane(shearG(10), 0, 1e-3, 0, Math.PI/6)[1], py:
  'ch04.stress_on_plane([[0.0, 10.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], 0.0, 1e-3, 0.0, np.pi/6)[1]', rtol: 1e-12}` · `{name: 'p -
  pbar expansion', js: 0 - meanP(stress(expG(1), 0, 1e-3, 3e-3)), py: '0.0 - ch04.mean_pressure(ch04.newtonian_stress([[1.0,
  0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], p=0.0, mu=1e-3, mu_v=3e-3))', rtol: 1e-12}` (= 0.009 Pa) · `{name: 'cube
  spin', js: cubeSpin(1, 0, 1000, 0.001), py: 'ch04.cube_spin_acceleration(1.0, 0.0, 1000.0, 0.001)', rtol: 1e-12}` ·
  `{name: 'K contraction', js: contractK(K4(2, 1, 1), symS())[0][1], py: 'ch04.linear_stress(ch04.isotropic_fourth_order(2.0,
  1.0, 1.0), [[1.0, 0.5, 0.0], [0.5, 2.0, 0.0], [0.0, 0.0, 3.0]])[0][1]', rtol: 1e-12}` · invariant `{name: 'rigid
  rotation no stress', js: Math.abs(stress(rotG(2), 0, 1e-3, 0)[0][1]), expect: 0, atol: 1e-15}`.
  (Note: the expansion row with μ_v = 3×10⁻³ gives p − p̄ = μ_v∇·u = 3×10⁻³ × 3 = 0.009 Pa.)
- **Fit plan:** 360×640: `element` (top) + `mohr` (bottom); `matrices` hidden — its key number (τ₁₂) is in the status and
  a readout; flow chips 5 short labels; custom sliders hidden unless 'custom'. Landscape phone: `element` + `mohr` in
  one row. Desktop: three views, `matrices` on the right spanning both rows.

### E4 · navier_stokes_term_balance
- **Title:** "Which terms of Navier–Stokes are awake here?" · **Summary:** "Pick an exact solution, click any point and read
  the five terms of the incompressible Navier–Stokes equation as bars that add to zero; scrub time and watch the balance
  move from local-vs-viscous to advective-vs-pressure." · **CORE:** C08 (also C06, N52–N55, N107) · **Reference:**
  `fid_formula_lab.html` (clickable terms, bars adding to a total) with `forced_damped_vibrations.html`'s explanation
  panel.
- **meta:** `viz:order 4` · `viz:sections 4.4 4.6` · `viz:equations 4.24 4.38 4.39a 4.39b 4.40 4.41` · `viz:fluidpy
  ch04.exact_solution ch04.ns_terms_preset ch04.viscous_force_forms ch04.stokes_first_problem ch04.plane_poiseuille` ·
  `viz:derivations D11 D12 D13`.
- **Physics (JS ↔ Python):** `solution(name, x, y, t, p) -> {u, v, p}` (closed forms) ↔ `ch04.exact_solution(name, [x, y,
  0], t, **p)` · `terms(name, x, y, t, comp, p) -> {local, advective, pressure, gravity, viscous, residual}` per unit
  **volume** (N/m³), from analytic derivatives of the closed forms ↔ `ch04.ns_terms_preset(name, x, y, t, component=comp,
  per="volume", **p)` (stencils; parity rtol 1e-6) · `viscForms(name, x, y, t, p) -> {lap, div2S, curlw}` ↔
  `ch04.viscous_force_forms` on `ch04.exact_field(name)` (parity through
  `ch04.ns_terms_preset(..., form="curl")`; `form="laplacian" | "div2S" | "curl"` selects how the viscous term is
  computed — Part C row 4.9) · `toy(x, a)` (compressible test field u = (a x², 0, 0) for D12, JS only,
  documented).
- **Solutions** (chips, short labels "Couette · Poiseuille · Stokes 1st · Taylor–Green · Lamb–Oseen · cylinder · solid
  body"): plane Couette (U = 1 m/s, h = 1 mm) · plane Poiseuille (G = 100 Pa/m, h = 1 mm) · Stokes' first problem
  (U = 1 m/s) · Taylor–Green (U₀ = 1 m/s, k = 2π m⁻¹, decaying) · Lamb–Oseen vortex (Γ = 2π m²/s, σ² = 4ν(t + t₀)) ·
  potential flow round a cylinder (Euler, U = 1 m/s, a = 1 m; viscous term exactly zero) · solid-body rotation
  (Ω = 1 rad/s; ω uniform, viscous force zero — the (4.40) paradox). Water ρ = 1000 kg/m³ unless ν is changed; gravity
  enters only as the hydrostatic part and is shown as a grey bar (0 in the horizontal components).
- **Views** (rows [1.3, 0.9]):
  1. `field` "The flow" — 2-D solutions: speed heat map + streamlines (`Viz.field.streamline`) + tracers; 1-D solutions
     (Couette, Poiseuille, Stokes 1st): the profile u(y, t) (teal) across the gap with a faint trail of earlier profiles
     and the steady ghost; the **probe** (click to move) with its velocity arrow. Title carries "Re_local = |u|l/ν = …".
  2. `bars` "The five terms at the probe (N/m³)" — local ρ∂u/∂t (blue), advective ρ(u·∇)u (teal), pressure −∇p (orange),
     gravity ρg (grey), viscous μ∇²u (rose), each signed, left side (local + advective) and right side (pressure + gravity
     + viscous) drawn as two stacked columns that must reach the same height; residual dot (black). Component chips x/y.
  3. `forms` "The viscous force three ways" (`hidePortrait: true`) — three arrows at the probe for μ∇²u, 2μ∂S_ij/∂x_i and
     −μ∇×ω (4.40) drawn on top of each other (they coincide), with their magnitudes; for 1-D solutions a time history of
     the local and viscous terms at the probe instead.
  Portrait: `field` + `bars`; `forms` hidden (its equality is in the status and Explain §5).
- **Controls:** `sol` chips · `nu` "Viscosity $\nu$" 10⁻⁷…10⁻² m²/s (log), default 10⁻⁶ · probe (pointer) · `comp` chips
  "x · y" · transport `t` (unsteady solutions: Stokes 0.01 → 10 s log-mapped, Taylor–Green 0 → 1/(2νk²), Lamb–Oseen
  0 → 10 s), rate 1, `end: 'hold'`; end card "Stokes: the wall's momentum has diffused 2√(νt) = … mm; local and viscous
  bars stayed equal all the way."
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **terms** (click a bar to isolate it) ·
  **transport** · **presets** (seven solutions + "Euler: ν → 0" {sol: 'cylinder'} + "paradox: solid body" {sol:
  'solid_body'}) · **inspector** (click the viscous bar: "μ∇²u ≈ μ[u(y+h) − 2u(y) + u(y−h)]/h² = 10⁻³ × (… − 2 × … + …)/10⁻¹²
  = −100.0 N/m³") · **status** ("⚖️ pressure ↔ viscous" when those two dominate (> 90 % of the largest) / "⏳ local ↔
  viscous (diffusion)" / "🌀 advective ↔ pressure (inertial)" / "🌀 inviscid: viscous force exactly 0 though μ ≠ 0" for
  the cylinder and solid body).
- **Readouts:** "Largest term" · "Residual" · "Local Re" · "Viscous force (3 forms)".
- **Explain:**
  0. *What the views show* — "**The flow** is an exact solution of $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$
     (4.39b). The probe reads it at one point. **The five terms** are the accelerations per unit volume at the probe:
     <b class="c-blue">local</b>, <b class="c-teal">advective</b>, <b class="c-orange">pressure</b>, grey gravity,
     <b class="c-rose">viscous</b>. The two columns must have the same height. **The viscous force three ways** (hidden
     on phones) checks (4.40)."
  1. *The solution at the probe* — formula and value, e.g. Poiseuille "$u=\frac{G}{2\mu}y(h-y)=\frac{100}{2\times10^{-3}}\times2.5\times10^{-4}\times7.5\times10^{-4}=$
     **9.375** mm/s".
  2. *Each term with numbers* — local ρ∂u/∂t = …; advective ρ(u·∇)u = …; pressure −∂p/∂x = G = **+100** N/m³; gravity 0
     (x); viscous μ∂²u/∂y² = −G = **−100** N/m³. (boxed sums: left = …, right = …)
  3. *Residual* — "left − right = … (round-off of the closed forms, ≤ 10⁻⁹ of the largest term)."
  4. *Which terms dominate* — ratios |advective|/|viscous| = local Re = |u|l/ν with l the solution's length scale =
     **…**; the regime word.
  5. *Viscous force three ways* — "$\mu\nabla^2\mathbf u$ = …, $2\mu\partial S_{ij}/\partial x_i$ = …, $-\mu\nabla\times\boldsymbol\omega$ =
     … N/m³ — equal because ∇·u = 0 (4.40)."
  6. *At the current time* — t = `live t` s; Stokes: layer thickness 2√(νt) = `live delta` mm; Taylor–Green: amplitude
     e^{−2νk²t} = `live amp`.
  7. *Reading the current setting* — Couette: "Nothing accelerates and there is no pressure gradient: the viscous force
     is zero everywhere (u is linear) — the shear stress is uniform." · Poiseuille: "Pressure pushes, friction resists,
     nothing accelerates — the balance of every fully developed pipe or channel." · Stokes: "The wall's motion diffuses
     into the fluid: local acceleration is paid by viscosity, nothing else is awake." · Taylor–Green: "Cells of vortices:
     advection is balanced by pressure, and the whole pattern decays as viscosity drains it (local ↔ viscous)." ·
     Lamb–Oseen: "Round the vortex, pressure provides the centripetal acceleration; along the circles viscosity spreads the
     core (local ↔ viscous)." · cylinder: "An Euler flow: advective acceleration is supplied by pressure; the viscous term
     is exactly zero because ω = 0 everywhere — viscosity is present but idle." · solid body: "The paradox resolved:
     uniform vorticity, so −μ∇×ω = 0 and no viscous force; pressure alone bends the paths."
- **Derivation tab:**
  - **D11** (7 steps) `view: 'bars'`; goal `set {sol: 'poiseuille'}`; each step lights the bar that its term becomes
    (step 2 pressure, step 4 viscous); step 7 **live** "at the probe: 0 + 0 = +100 + 0 − 100 N/m³"; interpret "(4.38) at your
    probe: …".
  - **D12** (7 steps) `view: 'bars'`; goal `set {sol: 'toy'}` (the compressible toy field, hidden chip) — steps 3–4 show
    the extra bar "(μ_v + ⅓μ)∂_j(∇·u)" (purple, only in this derivation); step 6 `set {sol: 'poiseuille'}` (∇·u = 0 removes
    it) `watch: "the purple bar vanishes"`; interpret "incompressible: μ∇²u is the whole viscous force".
  - **D13** (9 steps) `view: 'forms'` (phones: `bars`); goal `set {sol: 'lamb_oseen'}`; step 5 **live** "−μ(∇×ω)_x = …";
    step 8 `set {sol: 'solid_body'}` `watch: "all three arrows shrink to zero"`; interpret "the viscous force needs
    vorticity that varies in space".
- **Code:**
  ```python
  terms = ch04.ns_terms_preset("{{sol}}", {{x}}, {{y}}, t={{t}}, component={{c}}, per="volume", **{{params}})
  for k in ("local", "advective", "pressure", "gravity", "viscous"):
      print(k, terms[k])                    # N/m³ at the probe — {{vals}}
  print(terms["residual"])                  # ≈ 0: (4.39b) holds
  u_fn, p_fn = ch04.exact_field("{{sol}}", **{{params}})
  print(ch04.viscous_force_forms(u_fn, [{{x}}, {{y}}, 0.0], {{t}}, mu={{mu}}))   # three equal vectors (4.40)
  ```
- **Walkthrough (7 steps):** 1. "One equation, many balances" — "Honey in a pipe, air round a wing: same equation. Which
  of its five terms matter where?" `set {sol: 'poiseuille'}`, highlight `view:bars` · 2. "Two terms awake" — "In a
  channel, pressure pushes (orange) and friction resists (rose): +100 and −100 N/m³ at every point." `eq: 'ns'`,
  `readouts: ['largest', 'res']` · 3. "Where the viscous term comes from" — "Put the Newtonian stress into Cauchy's
  ρDu_j/Dt = ρg_j + ∂τ_ij/∂x_i (4.24): the divergence of 2μS is μ∇²u when ∇·u = 0." `derive: {id: 'D12', step: 6}` · 4. "Time wakes up" — "Stokes'
  first problem: press ▶. Local acceleration (blue) and viscosity (rose) trade exactly." `set {sol: 'stokes_first', t:
  0.01}`, `play: true` · 5. "Viscosity present, yet idle" — "The cylinder's potential flow: the viscous bar is exactly zero
  although μ ≠ 0 — ∇²u = −∇×ω = 0." `set {sol: 'cylinder'}`, `play: false`, `derive: {id: 'D13', step: 9}` · 6. "The
  paradox" — "Solid-body rotation has vorticity everywhere but no viscous force: ω is uniform, so its curl is zero."
  `set {sol: 'solid_body'}`, highlight `view:forms` · 7. "Your turn" — "Predict: in Taylor–Green, which two bars pair up?
  Click a cell corner and check." `set {sol: 'taylor_green'}`.
- **Equations:** `cauchy` "Cauchy" ref 'Eq. (4.24)' `\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}` ·
  `ns38` "Navier–Stokes, variable μ" ref 'Eq. (4.38)'
  `\begin{aligned}\rho\frac{Du_j}{Dt}=&-\frac{\partial p}{\partial x_j}+\rho g_j\\&+\frac{\partial}{\partial x_i}\Big[\mu\Big(\frac{\partial u_j}{\partial x_i}+\frac{\partial u_i}{\partial x_j}\Big)+\Big(\mu_v-\tfrac23\mu\Big)\frac{\partial u_m}{\partial x_m}\delta_{ij}\Big]\end{aligned}` ·
  `ns` "Incompressible" ref 'Eq. (4.39b)' `\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u` live "0 = 100 +
  0 − 100" · `visc` "Viscous force" ref 'Eq. (4.40)' `\mu\nabla^2\mathbf u=2\mu\frac{\partial S_{ij}}{\partial x_i}=-\mu\nabla\times\boldsymbol\omega` ·
  `euler` "Euler" ref 'Eq. (4.41)' `\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g`.
- **Check yourself:** (1) "In plane Couette flow, which bars are nonzero?" — "None: u is linear in y (∇²u = 0) and there is
  no pressure gradient; the shear stress is uniform and exerts no *net* force on a particle." `set {sol: 'couette'}` ·
  (2) "Double ν in Stokes' first problem: after 1 s the layer is how thick?" — "2√(νt) grows by √2: 2.0 → 2.8 mm." `set
  {sol: 'stokes_first', nu: 2e-6, t: 1}` · (3) "Why is the viscous term of the cylinder flow exactly zero?" — "The flow is
  irrotational, so ∇²u = ∇(∇·u) − ∇×ω = 0 (4.40)." · (4) "In the Lamb–Oseen vortex, what balances the pressure gradient
  in the radial direction?" — "The advective (centripetal) acceleration −u_θ²/r." `set {sol: 'lamb_oseen', comp: 'x'}`.
- **Selftest parity rows:** `{name: 'Poiseuille viscous', js: terms('poiseuille', 0, 2.5e-4, 0, 0, {G: 100, h: 1e-3, mu:
  1e-3}).viscous, py: 'ch04.ns_terms_preset("poiseuille", 0.0, 2.5e-4, component=0, per="volume", G=100.0, h=1e-3,
  mu=1e-3)["viscous"]', rtol: 1e-6}` · `{name: 'Poiseuille pressure', js: terms('poiseuille', 0, 2.5e-4, 0, 0, {G: 100, h: 1e-3, mu: 1e-3}).pressure,
  py: 'ch04.ns_terms_preset("poiseuille", 0.0, 2.5e-4, component=0, per="volume", G=100.0, h=1e-3, mu=1e-3)["pressure"]',
  rtol: 1e-9}` · `{name: 'Stokes local', js: terms('stokes_first', 0, 1e-3, 1.0, 0, {U: 1,
  nu: 1e-6, rho: 1000}).local, py: 'ch04.ns_terms_preset("stokes_first", 0.0, 1e-3, t=1.0, component=0, per="volume",
  U=1.0, nu=1e-6, rho=1000.0)["local"]', rtol: 1e-5}` (≈ 219.7 N/m³) · `{name: 'TG advective y', js: terms('taylor_green',
  0.1, 0.2, 0.0, 1, {U0: 1, k: 6.283185307179586, nu: 1e-6, rho: 1000}).advective, py: 'ch04.ns_terms_preset(
  "taylor_green", 0.1, 0.2, t=0.0, component=1, per="volume", U0=1.0, k=6.283185307179586, nu=1e-6,
  rho=1000.0)["advective"]', rtol: 1e-6}` · `{name: 'cylinder viscous = 0', js: terms('cylinder', 1.5, 0.7, 0, 0, {U: 1,
  a: 1, rho: 1, mu: 1e-3}).viscous, py: 'ch04.ns_terms_preset("cylinder", 1.5, 0.7, component=0, per="volume", U=1.0,
  a=1.0, rho=1.0, mu=1e-3)["viscous"]', rtol: 0, atol: 1e-6}` · `{name: 'Lamb-Oseen curl form', js: viscForms('lamb_oseen',
  0.8, 0.3, 1.0, {Gamma: 6.283185307179586, nu: 0.01, t0: 1}).curlw[0], py: 'ch04.ns_terms_preset("lamb_oseen", 0.8, 0.3,
  t=1.0, component=0, per="volume", form="curl", Gamma=6.283185307179586, nu=0.01, t0=1.0, rho=1000.0)["viscous"]', rtol:
  1e-5}` · invariant `{name: 'solid body viscous', js: terms('solid_body', 0.3, 0.1, 0, 0, {Omega: 1, rho: 1000}).viscous,
  expect: 0, atol: 1e-9}`.
- **Fit plan:** 360×640: `field` + `bars`; `forms` hidden; solution chips wrap (7 short labels, presets hidden on short
  portrait screens); status one line with the regime emoji and word. Landscape phone: `field` + `bars` in one row.
  Desktop: `field` left spanning both rows, `bars` and `forms` stacked right.

### E5 · rotating_frame_coriolis
- **Title:** "Why does a straight throw curve on a merry-go-round?" · **Summary:** "One ball, one clock, two observers:
  the inertial view shows a straight path over a turning table, the rotating view a curve bent by the Coriolis force,
  with every apparent force of the noninertial Navier–Stokes equation as an arrow and a bar." · **CORE:** C09 (also
  N56–N65) · **Reference:** `angular_frequency_explorer_1.html` (linked views on one clock, modes, ghost reference,
  end-of-run summary); links to ch03 `galilean_frames_cylinder`.
- **meta:** `viz:order 5` · `viz:sections 4.7` · `viz:equations 4.42 4.43 4.44 4.45` · `viz:fluidpy
  ch04.coriolis_projectile ch04.projectile_paths ch04.frame_acceleration_terms ch04.apparent_body_forces
  ch04.coriolis_acceleration ch04.centrifugal_acceleration ch04.effective_gravity` · `viz:derivations D14 D15 D16 D17 D18`.
- **Physics (JS ↔ Python):** `paths(u0, Om, t) -> {inertial: [x, y], rotating: [x′, y′]}` ↔ `ch04.projectile_paths(u0,
  Omega, t)` · `projectile(u0, Om, t) -> {forward, defl, deflSmall, angle}` ↔ `ch04.coriolis_projectile(u0, Omega, [t])` ·
  `frameTerms(a′, u′, x′, Ω, Ω̇, U̇)` ↔ `ch04.frame_acceleration_terms` · `forces(u′, x′, Ω, Ω̇, U̇, g)` ↔
  `ch04.apparent_body_forces` · `coriolis(Ω, u)` = −2Ω × u ↔ `ch04.coriolis_acceleration` · `centrifugal(Ω, x)` ↔
  `ch04.centrifugal_acceleration` · `effG(lat)` ↔ `ch04.effective_gravity(lat)` · `highLow(x, y, sense)` ↔
  `ch04.high_low_flow`.
- **Modes:** **turntable** (Ω = 0.1–2 rad/s, a ball launched from the centre) · **Earth's pole** (Ω = 7.292×10⁻⁵ rad/s,
  times in hours; distances in km) · **high / low** (radial outflow or inflow with Coriolis arrows; NH/SH) · **effective
  gravity** (latitude slider; exaggerated-Ω cross-section).
- **Views** (rows [1.2, 0.8]):
  1. `inertial` "Seen from outside (inertial)" (`equal: true`) — the turntable (or polar cap) drawn with a painted
     radius line turning at Ωt, the ball's **straight** rose path so far (bold) and to come (faint), the launch point.
     Effective-gravity mode: a meridian cross-section of an exaggerated rotating Earth with g_n (blue), Ω²R (purple) and
     g_e (black) arrows at the chosen latitude and the equipotential (dashed).
  2. `rotating` "Seen by the rider (rotating)" (`equal: true`) — the same motion in the rotating frame: the **curved** rose
     path, the small-angle ghost y′ = −Ωu₀t² (grey dashed), arrows on the ball for Coriolis −2Ω × u′ (amber) and
     centrifugal −Ω × (Ω × x′) (purple), the deflection δ marked. High/low mode: radial arrows (grey) and Coriolis arrows
     (amber) on a ring, with a curling trajectory.
  3. `bars` "Apparent forces per unit mass" (`hidePortrait: true`) — bars for g (grey, when relevant), −dU/dt (muted),
     −2Ω × u′ (amber), −(dΩ/dt) × x′ (muted), −Ω × (Ω × x′) (purple) for the component chosen (across-track by
     default), and their sum; a toggle "show the acceleration side of (4.43)" flips every bar's sign and relabels it
     (+2Ω × u′ "Coriolis term").
  Portrait: `inertial` + `rotating`; `bars` hidden (the Coriolis magnitude 2Ωu is in the status).
- **Controls:** `Om` "Rotation $\Omega$" 0…2 rad/s (turntable; a signed slider — negative = clockwise, the SH analogue),
  or in pole mode a chips control "NH · SH" · `u0` "Launch speed $u_0$" 0.5…10 m/s (turntable) / 1…50 m/s (pole) · `lat`
  "Latitude φ" 0…90° (effective gravity) · `sense` chips "high · low" (high/low mode) · `terms` toggle "show forces" ·
  `dUdt` "Frame acceleration" (optional, weightless preset). Transport `t`: turntable 0 → 2 s, pole 0 → 3 h, rate
  scaled; `end: 'hold'`; end card "The ball went straight; the table turned by Ωt = …°" (pole: "Deflection 9.34 km after
  1 h — Ωut² = 9.45 km").
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** (end card) · **modes** (4) ·
  **presets** "Earth pole, 1 h" {mode: 'pole', u0: 10, t: 3600} · "merry-go-round" {mode: 'turntable', Om: 0.5, u0: 3} ·
  "Southern Hemisphere" {mode: 'pole', hemi: 'SH'} · "weightless parabola" {mode: 'turntable', Om: 0, dUdt: 9.81} ·
  "Ω = 0" {Om: 0} · **terms** (bars) · **status** ("↪ NH: deflected right · 2Ωu = 1.46×10⁻³ m/s²" / "↩ SH: deflected left"
  / "⭕ Ω = 0: straight line, no apparent forces" / high: "↻ outflow from a high turns clockwise (NH)" / eff. gravity "g_e
  = 9.783 m/s², off the radial by 0.099°").
- **Readouts:** "Forward distance" · "Deflection (exact)" · "Ωu t² (small angle)" · "Angle Ωt" · "2Ωu".
- **Explain:**
  0. *What the views show* — "Left, an observer fixed to the stars: the ball flies in a <b class="c-rose">straight</b>
     line while the table turns under it. Right, the rider on the table: the same ball curves; to explain it the rider
     adds the <b class="c-amber">Coriolis</b> force −2Ω × u′ and the <b class="c-accent">centrifugal</b> force
     −Ω × (Ω × x′). The bars (hidden on phones) list every apparent force of (4.45)."
  1. *Your rotation and throw* — "Ω = … rad/s (period 2π/Ω = … s), u₀ = … m/s."
  2. *Coriolis force per unit mass* — "$\lvert-2\boldsymbol\Omega\times\mathbf u'\rvert=2\Omega u_0=2\times…\times…=$ **…**
     m/s², perpendicular to u′, to the right when Ω points up." (boxed)
  3. *Deflection* — "small-angle: $\delta\approx\tfrac12(2\Omega u_0)t^2=\Omega u_0t^2=$ … m; exact: $u_0t\sin\Omega t=$ **…** m;
     angle Ωt = … rad = …°." (boxed)
  4. *Centrifugal* — "$\lvert-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')\rvert=\Omega^2R=$ … m/s² at the ball's distance R =
     … m from the axis — it grows as the ball moves out (second order in Ωt, hence the gap between the two deflections)."
  5. *Which way* — "Ω up (NH): right; Ω down (SH): left; Ω = 0: straight."
  6. *Effective gravity* (mode or hint) — "$\mathbf g_e=\mathbf g_n+\Omega^2R\,\mathbf e_R$: at latitude φ = … the centrifugal part is
     Ω²a cos φ = … m/s²; g_e = … m/s², tilted … ° from the radial."
  7. *At the current time* — t = `live t`, forward `live fwd` m, deflection `live defl` m.
  8. *Reading the current setting* — "Nothing pushed the ball sideways: in the inertial view no force acts at all. The
     rider's Coriolis force is the bookkeeping of a straight path seen from a turning floor — which is why it does no work
     and never changes the speed." · pole: "On the Earth the effect is tiny per second but relentless: in an hour a
     10 m/s parcel is 9 km off course — why weather systems turn." · high: "Air leaving a high is turned clockwise (NH); it
     ends up circling the high rather than flowing out of it — geostrophic balance, Ch. 13." · Ω = 0: "No rotation, no
     apparent forces: both views agree."
- **Derivation tab:**
  - **D14** (7 steps) `view: 'rotating'`; goal `set {mode: 'turntable', Om: 0.5, t: 0}`; step 3 (the cone) draws e′₁
    sweeping its cone as t runs to 1 s (`set {t: 1}`); step 5 **live** "|de′₁/dt| = sin α |Ω| = 1 × 0.5 = 0.5 s⁻¹";
    interpret "u = U + u′ + Ω × x′: here Ω × x′ = …".
  - **D15** (12 steps, ★★★) `view: 'rotating'` (phones: `rotating`); goal `set {mode: 'turntable', Om: 0.5, u0: 3, t:
    0.5}`; step 6 (the first Ω × u′, from the turning basis) draws half the Coriolis arrow; step 9 (the second, from
    d(Ω × x′)/dt) completes it `watch: "the amber arrow doubles"`; step 10 (triple product) draws the centripetal arrow;
    step 12 **live** "a = 0 + a′ + 2 × 0.5 × 3 (⊥) + 0 + 0.25R"; interpret "the Coriolis term is 2Ωu′ = 3.0 m/s² here: one
    half from each source".
  - **D16** (6 steps) `view: 'bars'` (phones: `rotating`); goal `set {mode: 'turntable'}`; step 5 (move the frame terms
    to the right) flips every bar `watch: "accelerations become forces: every sign flips"`.
  - **D17** (5 steps) `view: 'rotating'`; goal `set {mode: 'pole', u0: 10, t: 0}`; step 3 overlays the parabola ½(2Ωu)t²;
    step 5 `set {t: 3600}` **live** "δ = 7.292×10⁻⁵ × 10 × 3600² = 9451 m"; interpret "exact 9342 m: the small-angle
    formula is 1.2 % high after an hour".
  - **D18** (5 steps) `view: 'inertial'`; goal `set {mode: 'effg', lat: 45}`; step 4 (the potential) draws the
    equipotentials; step 5 **live** "Ω²a cos φ = 0.0339 × 0.707 = 0.0240 m/s²"; interpret "g_e = … m/s², tilted … °".
- **Code:**
  ```python
  Om, u0, t = {{Om}}, {{u0}}, {{t}}               # rotation [rad/s], launch speed [m/s], time [s]
  pr = ch04.coriolis_projectile(u0, Om, [t])
  print(pr["forward"], pr["deflection"], pr["deflection_small"])   # {{fwd}}, {{defl}}, {{ds}} m
  f = ch04.apparent_body_forces([u0, 0, 0], [u0*t, 0, 0], [0, 0, Om])
  print(f["coriolis"], f["centrifugal"])          # −2Ω×u′ = {{cor}}, −Ω×(Ω×x′) = {{cen}} m/s²
  a = ch04.frame_acceleration_terms([0, 0, 0], [u0, 0, 0], [u0*t, 0, 0], [0, 0, Om])
  print(a["coriolis"])                            # +2Ω×u′ = {{acor}} (the (4.43) term: opposite sign)
  print(ch04.effective_gravity(np.deg2rad({{lat}})))   # g_e, tilt = {{ge}}
  ```
- **Walkthrough (7 steps):** 1. "Two observers, one ball" — "A ball rolls across a spinning table. Outside: a straight
  line. On the table: a curve. Who is right?" `set {mode: 'turntable', Om: 0.5, u0: 3, t: 0}`, highlight
  `view:inertial` · 2. "Press play" — "Watch both views on one clock: straight on the left, curving right on the right."
  `play: true` · 3. "The axes turn" — "The rider's axes rotate: de′/dt = Ω × e′. Velocity gains Ω × x′ — (4.42),
  $\mathbf u=\mathbf U+\mathbf u'+\boldsymbol\Omega\times\mathbf x'$." `play: false`, `derive: {id: 'D14', step: 5}` · 4. "Why the 2" —
  "Differentiate again: one Ω × u′ from the turning axes, one from the moving position —
  a = dU/dt + a′ + 2Ω×u′ + Ω̇×x′ + Ω×(Ω×x′) (4.43)." `derive: {id: 'D15', step: 9}`, `eq: 'acc'` · 5. "Forces the rider needs" — "Moved to the force side they flip sign:
  ρD′u′/Dt = −∇′p + ρ[g − dU/dt − 2Ω×u′ − Ω̇×x′ − Ω×(Ω×x′)] + μ∇′²u′ (4.45) — amber Coriolis, purple centrifugal." `terms: true`, highlight `view:rotating` · 6. "The Earth"
  — "At the pole, 10 m/s for an hour: deflected 9.3 km, angle 15° — exactly the Earth's turn." `set {mode: 'pole', u0: 10,
  t: 3600}`, `readouts: ['defl', 'angle']` · 7. "Your turn" — "Predict: in the Southern Hemisphere, left or right? Then
  switch and press ▶." `set {mode: 'pole', hemi: 'SH', t: 0}`.
- **Equations:** `vel` "Velocity" ref 'Eq. (4.42)' `\mathbf u=\mathbf U+\mathbf u'+\boldsymbol\Omega\times\mathbf x'` · `acc`
  "Acceleration" ref 'Eq. (4.43)'
  `\mathbf a=\frac{d\mathbf U}{dt}+\mathbf a'+2\boldsymbol\Omega\times\mathbf u'+\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'+\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')`
  live "2Ωu′ = 3.000 m/s²" · `ns` "Navier–Stokes, rotating frame" ref 'Eq. (4.45)'
  `\begin{aligned}\rho\frac{D'\mathbf u'}{Dt}=&-\nabla'p+\mu\nabla'^2\mathbf u'\\&+\rho\Big[\mathbf g-\frac{d\mathbf U}{dt}-2\boldsymbol\Omega\times\mathbf u'-\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')\Big]\end{aligned}` ·
  `defl` "Deflection" ref '§4.7' `\delta\approx\Omega ut^2` live "= 9451 m" · `cen` "Centrifugal" ref '§4.7'
  `-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')=\Omega^2R\,\mathbf e_R` · `hl` "Highs and lows" ref '§4.7'
  `-2\boldsymbol\Omega\times\mathbf u=-2\Omega_zu_R\,\mathbf e_\varphi`.
- **Check yourself:** (1) "Double u₀ at fixed Ω. What happens to the deflection after the same time?" — "It doubles
  (Ωut² ∝ u), but the angle Ωt does not change." `set {mode: 'turntable', u0: 6}` · (2) "Does the Coriolis force change
  the ball's speed?" — "No: it is always perpendicular to u′, so it does no work; only the direction turns." · (3) "Why is
  there a small gap between Ωut² and the exact deflection?" — "The centrifugal force (second order in Ωt) and the curving
  path; at Ωt = 15° the gap is 1.2 %." `set {mode: 'pole', t: 3600}` · (4) "Air flows into a low in the NH. Which way
  does it turn?" — "Counter-clockwise: −2Ω_z u_R e_φ with u_R < 0 points along +e_φ." `set {mode: 'highlow', sense:
  'low'}`.
- **Selftest parity rows:** `{name: 'deflection exact 1 h', js: projectile(10, 7.292115e-5, 3600).defl, py:
  'ch04.coriolis_projectile(10.0, 7.292115e-5, [3600.0])["deflection"][0]', rtol: 1e-10}` · `{name: 'deflection small', js:
  projectile(10, 7.292115e-5, 3600).deflSmall, py: 'ch04.coriolis_projectile(10.0, 7.292115e-5,
  [3600.0])["deflection_small"][0]', rtol: 1e-12}` · `{name: 'rotating path y', js: paths(3, 0.5, 1.0).rotating[1], py:
  'ch04.projectile_paths(3.0, 0.5, [1.0])[1][1][0]', rtol: 1e-12}` · `{name: 'Coriolis force y', js: coriolis([0, 0, 0.5],
  [3, 0, 0])[1], py: 'ch04.coriolis_acceleration([0.0, 0.0, 0.5], [3.0, 0.0, 0.0])[1]', rtol: 1e-12}` (= −3.0) · `{name:
  'Coriolis term y', js: frameTerms([0, 0, 0], [3, 0, 0], [1, 0, 0], [0, 0, 0.5]).coriolis[1], py:
  'ch04.frame_acceleration_terms([0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.5])["coriolis"][1]', rtol:
  1e-12}` (= +3.0) · `{name: 'centrifugal equator', js: centrifugal([0, 0, 7.292115e-5], [6378137, 0, 0])[0], py:
  'ch04.centrifugal_acceleration([0.0, 0.0, 7.292115e-5], [6378137.0, 0.0, 0.0])[0]', rtol: 1e-12}` · `{name: 'g_e 45',
  js: effG(Math.PI/4).g, py: 'ch04.effective_gravity(np.pi/4)[0]', rtol: 1e-10}` · invariant `{name: 'Coriolis does no
  work', js: dot(coriolis([0, 0, 0.5], [3, 1, 0]), [3, 1, 0]), expect: 0, atol: 1e-12}`.
- **Fit plan:** 360×640: `inertial` + `rotating` stacked (two squares, each ≥ 150 px); `bars` hidden (the status shows 2Ωu);
  mode chips 4 short labels ("table · pole · high/low · g_e"). Landscape phone: the two views side by side. Desktop:
  three views in one row, `bars` narrow.

### E6 · which_bernoulli
- **Title:** "Bernoulli is constant along what, exactly?" · **Summary:** "Six flows, each with its streamlines and two
  probes; the Bernoulli sum is plotted along and across streamlines as stacked bars, and the hypotheses you toggle pick
  which of the four Bernoulli equations applies." · **CORE:** C05, C11, C12 (also N30, N88, N90, N91, N96, N97–N101,
  N102, N103, N105) · **Reference:** `overfitting_curves.html` (a verdict line that changes with the regime) with
  `angular_frequency_explorer_1.html`'s highlighted "Right now" table.
- **meta:** `viz:order 6` · `viz:sections 4.4 4.9` · `viz:equations 4.19 4.68 4.69 4.70 4.71 4.72 4.74 4.75 4.78 4.82` ·
  `viz:fluidpy ch04.bernoulli_scenario ch04.which_bernoulli ch04.pitot_speed ch04.torricelli_speed ch04.rankine_bernoulli
  ch04.u_tube_column ch04.stagnation_temperature ch04.gauge_absorbed_bracket` · `viz:derivations D06 D24 D25 D26`.
- **Physics (JS ↔ Python):** `scenario(name, p) -> {sAlong, BAlong, sAcross, BAcross, terms, valid}` ↔
  `ch04.bernoulli_scenario(name, **p)` · `which(steady, viscous, irrot, barotropic, isentropic, constRho) -> [labels]` ↔
  `ch04.which_bernoulli` (exact-text parity on the joined labels) · `pitot(dp, rho)` ↔ `ch04.pitot_speed` · `torricelli(h)`
  ↔ `ch04.torricelli_speed` · `rankineB(r, Γ, σ)` ↔ `ch04.rankine_bernoulli(r, Gamma, sigma)["B"]` · `uTube(t, L, h0)` ↔
  `ch04.u_tube_column` · `T0(T, U)` ↔ `ch04.stagnation_temperature` · `gauge(B, sign)` ↔ `ch04.gauge_absorbed_bracket`.
- **Scenarios** (chips "pitot · orifice · vortex · cylinder · U-tube · nozzle"): pitot tube in an air stream (U, ρ = 1.2) ·
  tank and orifice (head h) · Rankine vortex in water (Γ, σ) · cylinder potential flow (U, a) · U-tube column (L, h₀, t) ·
  hot gas in a nozzle (T₀ = 300 K, U up to 300 m/s; perfect gas).
- **Views** (rows [1.2, 0.9]):
  1. `scene` "The flow and two probes" — the scenario's picture with streamlines (teal), two draggable probes A and B
     (black rings); for the pitot tube the manometer; for the orifice the tank and jet; for the vortex the whirlpool
     surface dip; for the U-tube the moving column; for the nozzle a converging duct shaded by temperature.
  2. `bplot` "B along and across" — left half: B along the streamline through A (teal line) — flat when a form holds; right
     half: B across streamlines from A to B (rose); beneath them stacked bars at A and B: ½u² (teal), p/ρ or ∫dp/ρ
     (orange), gz (blue), ∂φ/∂t (amber, U-tube), h (rose, nozzle energy form). Title carries "B_A = …, B_B = …".
  3. `table` "Which Bernoulli?" (`hidePortrait: true`) — the N102 decision table (rows (4.19), (4.71), (4.72), (4.75),
     (4.78), (4.82)) with columns steady / inviscid / irrotational / barotropic or ρ const / constant along, the rows
     `which(...)` returns lit, the current hypothesis toggles shown as ticks.
  Portrait: `scene` + `bplot`; `table` hidden (the valid labels are in the status).
- **Controls:** `sc` chips · `q` "Strength $q$" with a meaning line (pitot → U [m/s]; orifice → h [m]; vortex → Γ [m²/s];
  cylinder → U; U-tube → h₀ [m]; nozzle → U [m/s]) · `sig` "Core radius σ" (vortex only) · probes (pointer) · hypothesis
  toggles `steady`, `irrot`, `inviscid`, `constRho` (optional on phones — their state is fixed by the scenario and shown in
  the status; changing them only changes the verdict, not the flow) · transport `t` (U-tube only), 0 → 2 periods, loop.
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **presets** (six scenarios) · **status** verdict
  ("✅ (4.71): constant along this streamline" / "✅ (4.72): constant everywhere" / "⚠️ rotational: B differs across
  streamlines by …" / "⏱ unsteady: use (4.75) — ∂φ/∂t carries … m²/s²" / "🔥 energy form (4.78): h + ½u² constant" —
  each written with its equation, e.g. "✅ $\tfrac12u^2+\int dp/\rho+gz$ constant everywhere (4.72)") · **terms** (stacked
  B) · **inspector** (click a point: "B = ½|u|² + p/ρ + gz = ½ × 0.25 + (−875)/1000 + 0 = −0.750 m²/s²") · **table** with the
  current row lit (the "Right now" table of the reference).
- **Readouts:** "B at A" · "B at B" · "ΔB across" · "Predicted speed or pressure" · "Valid forms".
- **Explain:**
  0. *What the views show* — "Teal lines are streamlines; A and B are probes you can drag. **B along and across** shows the
     Bernoulli sum along the streamline through A (teal) and across the streamlines from A to B (rose); the stacked bars
     split it into <b class="c-teal">½u²</b>, <b class="c-orange">pressure</b>, <b class="c-blue">gz</b> and, when the flow
     is unsteady, <b class="c-amber">∂φ/∂t</b>. **Which Bernoulli?** (hidden on phones) lights the valid equations."
  1. *The flow's properties* — steady? inviscid? irrotational? barotropic? (ticks) → valid forms (from `which`): **…**.
  2. *B at the probes* — "B_A = ½u_A² + p_A/ρ + gz_A = … = **…**; B_B = … = **…** m²/s²." (boxed)
  3. *What it predicts* — pitot: "$U=\sqrt{2(p_2-p_1)/\rho}=\sqrt{2\times500/1.2}=$ **28.87** m/s"; orifice: "$u=\sqrt{2gh}=$ **4.43**
     m/s"; vortex: "the dip of the free surface Δz = −B(0)/g = … m"; U-tube: "p₁ − p₂ = ρL dU/dt = … Pa"; nozzle: "T = T₀ −
     U²/2C_p = … K".
  4. *Where it would fail* — "across the streamlines of the vortex core B changes (4.72 fails); in the U-tube steady
     Bernoulli (4.19) would predict p₁ = p₂ whenever U = 0 — wrong by ρL dU/dt; in the hot nozzle constant-ρ Bernoulli
     errs by … % while (4.78) holds."
  5. *The table* (hint on phones: "the decision table is in the landscape view; the status names the valid forms").
  6. *At the current time* (U-tube) — t = `live t`, h = `live h`, U = `live U`, dU/dt = `live a`, ∂φ/∂t term `live
     phit`.
  7. *Reading the current setting* — pitot: "One streamline from the free stream to the tube's mouth: (4.19) needs no
     irrotationality." · orifice: "The whole tank is (nearly) irrotational, so (4.72) lets us compare the surface with the
     jet although they are on different streamlines." · vortex, probes on one circle: "Same streamline, same B — (4.71)." ·
     vortex, probes across the core: "Different circles inside the rotational core: B differs — that difference is the
     whirlpool's dip." · cylinder: "Irrotational: one B for the whole flow (4.72)." · U-tube: "Unsteady and irrotational:
     only (4.75), $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int\frac{dp}{\rho}+gz=$ const, is valid; the amber
     bar carries the difference." · nozzle: "Compressible and fast: the energy form (4.78), $h+\tfrac12\lvert\mathbf u\rvert^2+gz$
     constant on streamlines, holds; it trades temperature for speed."
- **Derivation tab:**
  - **D06** (11 steps) `view: 'scene'`; goal `set {sc: 'pitot'}` with a stream-tube element drawn along the streamline
    into the tube; step 5 (side pressure force) highlights the conical wall; step 10 **live** "½(28.87)² + 0 + p/ρ = 416.7 +
    … = const"; interpret "(4.19) along your streamline: B_A = B_B = …".
  - **D24** (10 steps) `view: 'scene'`; goal `set {sc: 'vortex'}`; steps 3–7 (Lamb identity) draw u × ω arrows (radial)
    inside the core; step 8 (barotropic) shows the pressure function; interpret "∂u/∂t + ∇B = u × ω: here ∇B = u × ω = …
    m/s² outward".
  - **D25** (5 steps) `view: 'bplot'`; goal `set {sc: 'vortex'}`; step 3 **live** "u·∇B = u·(u × ω) = 0: B flat along the
    circle"; step 5 `set {sc: 'cylinder'}` `watch: "ω = 0: the rose curve goes flat too"`.
  - **D26** (8 steps) `view: 'bplot'`; goal `set {sc: 'u_tube', t: 0}`; step 6 (the gauge) **live** "with φ − ∫B dt′ the
    bracket is 0; with the book's + it is 2B = …" (`gauge(B, ±1)`); interpret "the amber bar is ∂φ/∂t = … m²/s²".
- **Code:**
  ```python
  sc = ch04.bernoulli_scenario("{{sc}}", **{{params}})
  print(sc["B_along"][[0, -1]], sc["B_across"][[0, -1]])   # along: {{ba}}; across: {{bx}}
  print(sc["valid"], sc["status"])                # {{valid}}
  print(ch04.which_bernoulli(steady={{st}}, viscous={{vi}}, irrotational={{ir}},
                             barotropic={{ba2}}, constant_density={{cr}}))   # {{labels}}
  print(ch04.pitot_speed({{p2}}, {{p1}}, 1.2))       # {{U}} m/s  (4.19) with u₂ = 0
  print(ch04.gauge_absorbed_bracket({{B}}, -1.0))    # 0.0: φ − ∫B dt′ absorbs B(t)
  ```
- **Walkthrough (7 steps):** 1. "One name, four equations" — "Bernoulli's sum: constant along what?" `set {sc:
  'pitot'}`, highlight `view:bplot` · 2. "Along one streamline" — "Steady, frictionless, constant ρ: $\tfrac12U^2+gz+p/\rho$ is
  constant along a streamline (4.19). The pitot tube reads 28.9 m/s from 500 Pa." `eq: 'b19'`, `readouts: ['pred']`,
  `derive: {id: 'D06', step: 10}` · 3. "A whirlpool" — "Rankine vortex: B is flat along each circle but not across the
  core." `set {sc: 'vortex'}`, `inspect: true` · 4. "Why" — "Euler with the Lamb identity: ∂u/∂t + ∇B = u × ω (4.69). Steady:
  ∇B ⟂ u, so B = ½u² + ∫dp/ρ + gz is constant along streamlines (4.71)." `derive: {id: 'D25', step: 3}` · 5. "No vorticity, one constant" —
  "The cylinder flow is irrotational: ½u² + p/ρ + gz takes one value everywhere (4.72)." `set {sc: 'cylinder'}` · 6. "Unsteady" — "The U-tube:
  stop the column; the pressures still differ. The ∂φ/∂t term of ∂φ/∂t + ½|∇φ|² + ∫dp/ρ + gz = const (4.75) carries it." `set {sc: 'u_tube', t: 0}`, `play:
  false`, highlight `view:bplot` · 7. "Your turn" — "Predict: for the hot nozzle at 250 m/s, will constant-ρ Bernoulli
  hold? Toggle the hypotheses and check the table." `set {sc: 'nozzle', q: 250}`.
- **Equations:** `b19` "Along a streamline" ref 'Eq. (4.19)' `\tfrac12U^2+gz+p/\rho=\text{const}` live "416.7 + 0 + 83333.3
  = 83750.0" · `lamb` "Lamb identity" ref 'Eq. (4.68)' `u_i\frac{\partial u_j}{\partial x_i}=-(\mathbf u\times\boldsymbol\omega)_j+\frac{\partial}{\partial x_j}\big(\tfrac12u_i^2\big)` ·
  `bfun` "Bernoulli function" ref 'Eq. (4.69)' `\frac{\partial\mathbf u}{\partial t}+\nabla B=\mathbf u\times\boldsymbol\omega` · `b71`
  "Streamlines and vortex lines" ref 'Eq. (4.71)–(4.72)' `B=\tfrac12u^2+\int\frac{dp}{\rho}+gz=\text{const}` · `b75` "Unsteady"
  ref 'Eq. (4.75)' `\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int\frac{dp}{\rho}+gz=\text{const}` · `b78` "Energy
  form" ref 'Eq. (4.78)' `h+\tfrac12\lvert\mathbf u\rvert^2+gz=\text{const}` live "T₀ = 304.98 K".
- **Check yourself:** (1) "Probes A and B on two different circles outside the vortex core: equal B?" — "Yes: outside the
  core the flow is irrotational, so (4.72) makes B the same everywhere there." `set {sc: 'vortex'}` · (2) "Does a pitot
  tube need an irrotational flow?" — "No: its two points lie on one streamline; (4.19) is enough." · (3) "In the U-tube at
  the turning point (U = 0), what carries the pressure difference?" — "The unsteady term ∂φ/∂t (ρL dU/dt = 981 Pa for
  h₀ = 5 cm)." `set {sc: 'u_tube', t: 0}` · (4) "The book redefines φ → φ + ∫B dt′. What goes wrong?" — "The bracket becomes
  2B(t) instead of 0; the correct gauge is φ − ∫B dt′." (open the Derivation tab at D26 step 6).
- **Selftest parity rows:** `{name: 'pitot 500 Pa', js: pitot(500, 1.2), py: 'ch04.pitot_speed(500.0, 0.0, 1.2)', rtol:
  1e-12}` · `{name: 'torricelli 1 m', js: torricelli(1), py: 'ch04.torricelli_speed(1.0)', rtol: 1e-12}` · `{name: 'Rankine B
  at 0.5', js: rankineB(0.5, 6.283185307179586, 1), py: 'ch04.rankine_bernoulli(0.5, 6.283185307179586, 1.0)["B"]', rtol:
  1e-12}` · `{name: 'U-tube dp at t=0', js: uTube(0, 1, 0.05).dp, py: 'ch04.u_tube_column(0.0, 1.0, 0.05)["dp"]', rtol:
  1e-10}` · `{name: 'stagnation T', js: T0(300, 100), py: 'ch04.stagnation_temperature(300.0, 100.0)', rtol: 1e-12}` ·
  `{name: 'which (steady inviscid rotational barotropic)', js: which(true, false, false, true, false, false).join(','),
  py: 'ch04.which_bernoulli_text(True, False, False, True)', exact: true}` (exact-text compare; `which_bernoulli_text`
  returns the comma-joined labels, "4.71" here — Part C row 7.7) · `{name: 'which (irrotational, ρ const)', js:
  which(true, false, true, true, false, true).join(','), py: 'ch04.which_bernoulli_text(True, False, True, True,
  constant_density=True)', exact: true}` (= "4.19,4.71,4.72,4.75") ·
  `{name: 'gauge plus sign', js: gauge(3, 1), py: 'ch04.gauge_absorbed_bracket(3.0, 1.0)', rtol: 1e-12}` (= 6) ·
  invariant `{name: 'cylinder B flat', js: spread(scenario('cylinder', {U: 1, a: 1}).BAcross), expect: 0, atol: 1e-9}`.
- **Fit plan:** 360×640: `scene` + `bplot`; `table` hidden; hypothesis toggles hidden on phones (fixed by the scenario;
  status names the valid forms with their equations, two lines max: "✅ (4.72) B = ½u² + p/ρ + gz same everywhere");
  scenario chips 6 short labels. Landscape phone: `scene` + `bplot` in a row. Desktop: `scene` left, `bplot` and `table`
  stacked right.

### E7 · viscous_dissipation_heating
- **Title:** "Where does the energy go when viscosity stops a flow?" · **Summary:** "A sheared channel on one clock: the
  velocity profile, the dissipation ε(y) where the shear is, the temperature rising to its steady parabola, and energy
  bars showing the wall's work arriving as heat — with a 'negative viscosity' preset that breaks the second law." ·
  **CORE:** C10 (also N74, N76, N77, N79, N81) · **Reference:** `forced_damped_vibrations.html` (a transient settling to a
  steady state on one time slider, boxed numbers, a regime reading).
- **meta:** `viz:order 7` · `viz:sections 4.8` · `viz:equations 4.53 4.54 4.56 4.57 4.58 4.60 4.63` · `viz:fluidpy
  ch04.couette_heating ch04.couette_heating_transient ch04.dissipation_rate ch04.entropy_production` ·
  `viz:derivations D21 D22 D23`.
- **Physics (JS ↔ Python):** `steady(y, p) -> {u, dudy, rhoEps, T, qBottom, qTop, workIn, heatOut, dTmax}` (Couette plus
  an optional pressure gradient: u = Uy/h + (G/2μ)y(h − y); ρε = μ(du/dy)²; T from k T″ = −ρε with T = T₀ at both walls,
  closed form) ↔ `ch04.couette_heating(y, U, h, mu, k, T0, dpdx=-G)` · `transient(y, t, p)` (sine series, 200 terms) ↔
  `ch04.couette_heating_transient(y, t, U, h, mu, k, rho, cp, T0)` (Couette part; the pressure-gradient part uses the
  same series with its own source) · `eps(G, rho, mu, muv)` ↔ `ch04.dissipation_rate(G, rho, mu, mu_v)` · `sgen(dTdy, T,
  k, rho, eps)` ↔ `ch04.entropy_production([0, dTdy, 0], T, k, rho, eps)`.
- **Views** (rows [1.2, 0.9]):
  1. `channel` "The sheared gap" — the two walls (top moving, arrow U), the velocity profile u(y) (teal) with arrows, and
     beside it the dissipation density ρε(y) as a heat strip (rose, darker = more heating); legend in the title
     ("ρε max … W/m³").
  2. `temp` "Temperature across the gap" — T(y, t) − T₀ (rose) rising toward the steady parabola (dashed ghost), faint trail
     of earlier profiles, the heat-flux arrows at both walls (out of the fluid); title "ΔT_max = … K (steady … K)".
  3. `budget` "Energy per unit wall area" (`hidePortrait: true`) — bars: work done by the moving wall τU (teal), heat
     generated ∫ρε dy (rose), heat stored ∂/∂t∫ρC_pT dy (purple, → 0), heat conducted out through the walls (orange), and
     the balance (black); a second small bar "entropy production at mid-gap ≥ 0" (green, or red below zero).
  Portrait: `channel` + `temp`; `budget` hidden (its key number — heat out = work in — in the status).
- **Controls:** `flow` chips "Couette · Poiseuille · both" · `U` "Wall speed $U$" 0…10 m/s, default 1 (hidden for
  Poiseuille) · `G` "Pressure gradient $G=-dp/dx$" 0…10⁵ Pa/m (Poiseuille, both) · `mu` "Viscosity $\mu$" 10⁻⁵…1 Pa s (log),
  default 10⁻³ · `k` "Conductivity $k$" 0.02…1 W/(m K) (log), default 0.6 · gap h = 1 mm, ρ = 1000, C_p = 4186 fixed (fluid
  presets change them). Transport `tau` = κt/h² from 0 to 0.6 (readout in seconds), rate 0.05/s, `end: 'hold'`; end card
  "Steady: shear work in = heat out = … W/m²; ΔT_max = … K".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **terms** (energy budget) · **transport** (end
  card) · **presets** "water Couette 1 m/s" {flow: 'couette', U: 1, mu: 1e-3, k: 0.6} · "oil bearing 5 m/s" {mu: 0.1, k:
  0.15, U: 5, rho: 900, cp: 2000} · "air" {mu: 1.8e-5, k: 0.026, rho: 1.2, cp: 1005, U: 10} · "Poiseuille water" {flow:
  'poiseuille', G: 1e4} · "negative μ (forbidden)" {mu: -1e-3} · **status** ("🔥 heating: ε > 0 · shear work in = heat out =
  1.00 W/m²" / "⛔ second law violated: μ < 0 gives ε < 0 and negative entropy production").
- **Readouts:** "ρε (max)" · "ΔT_max" · "Heat out" · "Entropy production".
- **Explain:**
  0. *What the views show* — "The top wall slides at U over a fluid layer h = 1 mm thick. <b class="c-teal">Teal</b>:
     velocity. <b class="c-rose">Rose</b>: where viscosity turns motion into heat (ε) and the temperature it causes. The
     budget (hidden on phones) balances the wall's work against the heat leaving through the walls."
  1. *Shear rate* — "du/dy = U/h = … s⁻¹ (Couette; the pressure-driven part adds (G/2μ)(h − 2y))."
  2. *Dissipation* — "(4.58), $\varepsilon=2\nu\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)^2+\frac{\mu_v}{\rho}S_{mm}^2$, for a parallel
     flow is ε = ν(du/dy)²: ρε = μ(U/h)² = … × … = **…** W/m³." (boxed)
  3. *Total heating* — "∫₀ʰ ρε dy = … W/m² = the wall's work τU = μ(U/h)U = **…** W/m² — all the work becomes heat." (boxed)
  4. *Temperature* — "steady: k T″ = −ρε ⇒ $T-T_0=\frac{\rho\varepsilon}{2k}y(h-y)$, peak at mid-gap
     $\Delta T_{max}=\mu U^2/(8k)=$ **…** K; each wall conducts … W/m² away." (boxed)
  5. *Time to heat up* — "the diffusion time h²/κ = … s (κ = k/ρC_p); after about 0.5 h²/κ the profile is within 1 % of
     steady."
  6. *Entropy production* — "(4.63), $\frac{k}{\rho T^2}\lvert\nabla T\rvert^2+\frac\varepsilon T$: at mid-gap … + … = **…**
     W/(kg K) ≥ 0." (with negative μ: "< 0 — forbidden by the second law, hence μ ≥ 0")
  7. *At the current time* — t = `live t` s (τ = `live tau`), ΔT_max now `live dT` K, heat stored `live stored` W/m².
  8. *Reading the current setting* — water: "Tiny heating (a fifth of a millikelvin at 1 m/s): viscous heating is
     negligible in most water and air flows — why the Boussinesq heat equation drops it (C13)." · oil bearing: "Hundreds of
     times more heating: bearings need cooling; viscosity itself falls as the oil warms (a feedback not in this model)." ·
     Poiseuille: "The dissipation is largest at the walls, where the shear is; the temperature bump is flatter." ·
     negative μ: "ε would be negative — heat turning spontaneously into motion. The second law forbids it: μ, μ_v, k ≥ 0
     (4.63)."
- **Derivation tab:**
  - **D21** (6 steps) `view: 'budget'` (phones: `channel`); goal `set {flow: 'couette', U: 1}`; step 3 lights the
    "kinetic energy" bar (zero in steady Couette); step 5 splits the stress work into force work and deformation work;
    interpret "the mechanical budget: the wall's force work reaches the fluid and is used up by deformation work".
  - **D22** (8 steps) `view: 'temp'`; goal as above; step 6 (the symmetric-σ step) `watch: "only S contributes"`; step 7
    moves the deformation-work bar into the heat side; step 8 **live** "De/Dt = 0 + ε − (1/ρ)∂q/∂y: 1.0 = 1.0 W/kg at
    steady state"; interpret "the fluid's internal energy gains ε and loses conduction".
  - **D23** (8 steps) `view: 'channel'`; goal as above; step 5 draws ε(y) as the square of the shear; step 8 **live**
    "ε = 2ν(S₁₂² + S₂₁²) = 2 × 10⁻⁶ × 2 × 500² = 1.0 W/kg"; interpret "a sum of squares: never negative when μ, μ_v ≥ 0".
- **Code:**
  ```python
  y = np.linspace(0, 1e-3, 5)                              # across the gap [m]
  c = ch04.couette_heating(y, U={{U}}, h=1e-3, mu={{mu}}, k={{k}}, dpdx={{dpdx}})
  print(c["eps"].max(), c["work_in"], c["heat_out"])       # {{eps}} W/m³, {{w}} = {{q}} W/m²
  print(c["dT_max"])                                       # {{dT}} K at mid-gap
  G = [[0, {{S}}, 0], [0, 0, 0], [0, 0, 0]]                # du/dy at the wall [1/s]
  print(ch04.dissipation_rate(G, {{rho}}, {{mu}}, form="both"))   # ({{e1}}, {{e1}}) W/kg
  print(ch04.entropy_production([0, 0, 0], {{T0}}, {{k}}, {{rho}}, {{e1}}))   # {{sg}} ≥ 0
  ```
- **Walkthrough (6 steps):** 1. "The swirl dies" — "Stir, let go: motion stops. Where does its energy go?" `set
  {flow: 'couette', U: 1}`, highlight `view:channel` · 2. "Two kinds of stress work" — "Stress work splits:
  $\frac{\partial}{\partial x_i}(\tau_{ij}u_j)=\tau_{ij}\frac{\partial u_j}{\partial x_i}+u_j\frac{\partial\tau_{ij}}{\partial x_i}$ (4.54) —
  deformation work plus force work (the second changes kinetic energy)." `derive: {id: 'D21', step: 5}` · 3. "The heat side" — "Subtract
  the kinetic-energy equation from the total: what is left is
  $\frac{De}{Dt}=-p\frac{Dv}{Dt}+\varepsilon-\frac1\rho\frac{\partial q_i}{\partial x_i}$ (4.57) — the fluid gains ε." `derive: {id: 'D22', step: 8}`,
  `eq: 'int'` · 4. "Always positive" — "$\varepsilon=2\nu\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)^2+\frac{\mu_v}\rho S_{mm}^2$ (4.58) is a sum of
  squares: 1.0 W/kg here, never negative." `readouts: ['eps']`,
  `derive: {id: 'D23', step: 8}` · 5. "Heating up" — "Press ▶: the temperature climbs to its steady parabola; then heat out
  equals the wall's work." `play: true`, highlight `view:temp` · 6. "Your turn" — "Predict: the oil bearing at 5 m/s —
  how much warmer than water at 1 m/s? Then choose the preset." `set {U: 1}`, `play: false`.
- **Equations:** `tot` "Total energy" ref 'Eq. (4.53)'
  `\begin{aligned}&\tfrac{\partial}{\partial t}\big(\rho[e+\tfrac12u_j^2]\big)+\tfrac{\partial}{\partial x_i}\big(\rho[e+\tfrac12u_j^2]u_i\big)\\&=\rho g_iu_i+\tfrac{\partial}{\partial x_i}(\tau_{ij}u_j)-\tfrac{\partial q_i}{\partial x_i}\end{aligned}` ·
  `mech` "Mechanical energy" ref 'Eq. (4.56)'
  `\rho\frac{D}{Dt}\big(\tfrac12u_j^2\big)=\rho g_ju_j-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}` · `int`
  "Internal energy" ref 'Eq. (4.57)' `\frac{De}{Dt}=-p\frac{Dv}{Dt}+\frac1\rho\sigma_{ij}S_{ij}-\frac1\rho\frac{\partial q_i}{\partial x_i}` ·
  `eps` "Dissipation" ref 'Eq. (4.58)' `\varepsilon=2\nu\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)^2+\frac{\mu_v}{\rho}S_{mm}^2\ge0` live
  "= 1.0 W/kg" · `sgen` "Entropy production" ref 'Eq. (4.63)' `\frac{k}{\rho T^2}\lvert\nabla T\rvert^2+\frac\varepsilon T\ge0`.
- **Check yourself:** (1) "Double U. How does the peak temperature rise change?" — "×4: ΔT_max = μU²/(8k)." `set {U: 2}`
  · (2) "In steady Couette flow the kinetic energy of each particle is constant. Where does the wall's work go?" — "All
  into deformation work, i.e. heat (ρε), which is conducted out through the walls." · (3) "Why is μ < 0 impossible?" — "It
  would make ε and the entropy production negative (4.63)." `set {mu: -1e-3}` · (4) "Where is the dissipation largest in
  Poiseuille flow?" — "At the walls, where |du/dy| is largest; zero on the centreline." `set {flow: 'poiseuille'}`.
- **Selftest parity rows:** `{name: 'rho eps Couette', js: steady(5e-4, {U: 1, h: 1e-3, mu: 1e-3, k: 0.6, G: 0}).rhoEps,
  py: 'ch04.couette_heating(5e-4, U=1.0, h=1e-3, mu=1e-3, k=0.6)["eps"]', rtol: 1e-12}` · `{name: 'dT max', js: steady(5e-4,
  {U: 1, h: 1e-3, mu: 1e-3, k: 0.6, G: 0}).dTmax, py: 'ch04.couette_heating(5e-4, U=1.0, h=1e-3, mu=1e-3, k=0.6)["dT_max"]', rtol: 1e-12}` · `{name:
  'heat out', js: steady(5e-4, {U: 1, h: 1e-3, mu: 1e-3, k: 0.6, G: 0}).heatOut, py: 'ch04.couette_heating(5e-4, U=1.0, h=1e-3, mu=1e-3,
  k=0.6)["heat_out"]', rtol: 1e-12}` · `{name: 'transient mid-gap', js: transient(5e-4, 2.0, {U: 1, h: 1e-3, mu:
  1e-3, k: 0.6, rho: 1000, cp: 4186}), py: 'ch04.couette_heating_transient(5e-4, 2.0, U=1.0, h=1e-3, mu=1e-3, k=0.6,
  rho=1000.0, cp=4186.0) - 293.15', rtol: 1e-8}` · `{name: 'eps contraction', js: eps([[0, 1000, 0], [0, 0, 0], [0, 0, 0]],
  1000, 1e-3, 0), py: 'ch04.dissipation_rate([[0.0, 1000.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], 1000.0, 1e-3)', rtol:
  1e-12}` · `{name: 'entropy production', js: sgen(0, 293.15, 0.6, 1000, 1.0), py: 'ch04.entropy_production([0.0, 0.0,
  0.0], 293.15, 0.6, 1000.0, 1.0)', rtol: 1e-12}` · invariant `{name: 'work in = heat out', js: steady(5e-4, {U: 1, h: 1e-3, mu: 1e-3, k: 0.6, G: 0}).workIn -
  steady(5e-4, {U: 1, h: 1e-3, mu: 1e-3, k: 0.6, G: 0}).heatOut, expect: 0, atol: 1e-12}`.
- **Fit plan:** 360×640: `channel` + `temp`; `budget` hidden (status carries "work in = heat out = … W/m²"); ≤ 4 controls
  (flow, U or G, μ, k). Landscape phone: `channel` + `temp` side by side. Desktop: `channel` and `temp` left, `budget`
  right.

### E8 · boussinesq_buoyancy
- **Title:** "Density hardly changes — so why does it drive the flow?" · **Summary:** "Choose a fluid, a temperature
  contrast and a depth: three small numbers are checked against their thresholds while a warm blob rises and spreads under
  the terms Boussinesq keeps; switch buoyancy off and the blob only diffuses." · **CORE:** C13 (also R10, N106–N114, N137)
  · **Reference:** `angular_frequency_explorer_1.html` ("Right now" notes with a highlighted table, modes).
- **meta:** `viz:order 8` · `viz:sections 4.9` · `viz:equations 4.84 4.86 4.88 4.89` · `viz:fluidpy ch04.boussinesq_validity
  ch04.boussinesq_scenario ch04.boussinesq_density ch04.buoyancy ch04.gaussian_blob_advection_diffusion ch04.blob_rise` ·
  `viz:derivations D27 D28`.
- **Physics (JS ↔ Python):** `validity(alpha, dT, L, U, c, g, nu, cp)` ↔ `ch04.boussinesq_validity` · `scenario(name)` ↔
  `ch04.boussinesq_scenario(name)` · `blob(x, z, t, p)` (Gaussian advected by the parcel's rise z_c(t) and spreading with
  σ² = σ₀² + 2κt — an exact solution of (4.89) because the translation velocity is uniform in space at each instant) ↔
  `ch04.gaussian_blob_advection_diffusion` (with U the mean velocity over [0, t]) · `rise(t, gp, tauD)` (our one-line
  parcel model, labelled as such: dw/dt = b − w/τ_d with b = g′, so w = g′τ_d(1 − e^{−t/τ_d})) ↔ `ch04.blob_rise(t,
  g_prime, tau_d)` · `buoy(rhoPert, rho0)` ↔ `ch04.buoyancy` · `rhoLin(T, rho0, alpha, T0)` ↔ `ch04.boussinesq_density`.
- **Modes / scenarios** (chips "lake · thermocline · lab tank · ABL · deep atmosphere"): each sets fluid properties (α, ν,
  C_p, c, κ), δT, L and U from `boussinesq_scenario`.
- **Views** (rows [1.2, 0.9]):
  1. `tank` "The blob" — a vertical tank (or a slab of atmosphere) with the temperature anomaly T′ as a heat map (rose on
     white), the blob's centre path (dashed), σ as a dashed circle, a height scale in the scenario's units; title "ρ′/ρ₀ =
     −αδT = … (…%)".
  2. `bars` "Vertical momentum per unit mass (m/s²)" — inertia Dw/Dt (teal), buoyancy −gρ′/ρ₀ (blue), drag/pressure −w/τ_d
     (orange, the parcel model's stand-in), viscous ν∇²w (rose, tiny) and the **dropped** term (ρ′/ρ₀)Dw/Dt (grey, a hair
     — its size written next to it); toggle "keep buoyancy".
  3. `validity` "Is Boussinesq valid here?" (`hidePortrait: true`) — three horizontal log bars αδT, L/(c²/g),
     νU/(C_pδT L) with the threshold 0.1 dashed and the scenario table (rows = the five scenarios, columns = the three
     numbers) with the current row lit.
  Portrait: `tank` + `bars`; `validity` hidden (the verdict and the largest ratio in the status).
- **Controls:** `sc` chips · `dT` "Temperature contrast δT" 0.1…50 K (log), default by scenario · `L` "Depth scale L"
  (log; 0.1 m … 20 km by scenario) · `buoy` toggle "keep buoyancy" (default on) · transport `t` 0 → 5τ_d, rate by
  scenario, `end: 'hold'`; end card "The blob rose … m and spread to σ = … m; density differed by only …%."
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **presets** (five scenarios) · **status** verdict
  ("✅ Boussinesq valid · αδT = 2.0×10⁻³, L/(c²/g) = …" / "⚠️ deep layer: L ≈ c²/g — use the anelastic or compressible
  equations") · **terms** (bars) · **transport** · the **table** of real values with the current row lit.
- **Readouts:** "αδT" · "g′" · "L/(c²/g)" · "Heating ratio" · "Rise speed".
- **Explain:**
  0. *What the views show* — "A warm blob in a cooler fluid (<b class="c-rose">rose</b> = warmer). **Vertical momentum**
     shows the accelerations the Boussinesq equations keep — <b class="c-blue">buoyancy</b> ρ′g/ρ₀ is the engine — and, in
     grey, the one they drop. **Is Boussinesq valid?** (hidden on phones) compares three small numbers with 0.1."
  1. *Density change* — "α = … K⁻¹ (water: a measured value; gas: 1/T), δT = … K ⇒ δρ/ρ = αδT = **…**." (boxed)
  2. *Reduced gravity* — "g′ = gαδT = 9.81 × … = **…** m/s² — the blob's buoyancy per unit mass." (boxed)
  3. *Depth limit* — "sound speed c = … m/s, c²/g = … km; L/(c²/g) = … (must be ≪ 1)."
  4. *Viscous heating* — "νU/(C_pδT L) = … × …/(… × … × …) = **…** — negligible."
  5. *What is kept and what is dropped* — "(4.86), $\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u$:
     the buoyancy bar (… m/s²) is kept because g is large; the dropped term (ρ′/ρ₀)Dw/Dt = … m/s² is αδT times smaller than
     the inertia."
  6. *The blob* — "rise speed from the parcel model w → g′τ_d = … m/s (scale √(g′L) = … m/s); spread σ² = σ₀² + 2κt with
     κ = … m²/s: σ = `live sigma` m at t = `live t` s — the exact solution of $\frac{DT}{Dt}=\kappa\nabla^2T$ (4.89)."
  7. *The table* (hint on phones) — the five scenarios' numbers with the current row lit.
  8. *Reading the current setting* — valid: "Density differs by …%, yet the blob rises: a tiny density difference acts on
     the whole mass through g. Everywhere else the fluid behaves as incompressible." · buoyancy off: "Without the ρ′g term
     nothing drives the blob: it only diffuses — Boussinesq keeps exactly the term that matters." · deep atmosphere: "Over
     10 km the air's own weight compresses it by tens of percent: L is comparable to c²/g and density can no longer be
     treated as constant outside the gravity term."
- **Derivation tab:**
  - **D27** (9 steps) `view: 'bars'`; goal `set {sc: 'lake'}`; step 1 (subtract the hydrostatic state) draws p_s(z) beside
    the tank; step 6 (why keep ρ′g) **live** "ρ′g/ρ₀ = … m/s² vs (ρ′/ρ₀)Dw/Dt = … m/s²" `watch: "the grey bar is a
    thousandth of the blue one"`; interpret "(4.86) with your numbers: buoyancy … m/s²".
  - **D28** (8 steps) `view: 'bars'` (phones: `tank`); goal `set {sc: 'abl'}` (a gas, where the C_p argument applies);
    step 4 **live** "−p∇·u = −ρ(C_p − C_v)DT/Dt = −1.2 × 287 × DT/Dt"; step 5 `watch: "the pressure-work bar folds into the
    heating: C_v becomes C_p"`; interpret "(4.89) with κ = k/ρC_p = … m²/s".
- **Code:**
  ```python
  p = ch04.boussinesq_scenario("{{sc}}")          # α, δT, L, U, ν, C_p, c for this setting
  v = ch04.boussinesq_validity(p["alpha"], {{dT}}, {{L}}, p["U"], c=p["c"], nu=p["nu"], cp=p["cp"])
  print(v["alpha_dT"], v["L_over_Hc"], v["heating_ratio"])   # {{a}}, {{l}}, {{h}}
  print(v["g_prime"], v["verdict"])               # {{gp}} m/s²  {{verdict}}
  r = ch04.blob_rise({{t}}, v["g_prime"], {{tauD}})
  print(r["w"], r["z"])                           # rise speed {{w}} m/s, height {{z}} m
  ```
- **Walkthrough (6 steps):** 1. "A fraction of a percent" — "A lake 10 K warmer at the top differs in density by 0.2 %.
  Why does that matter?" `set {sc: 'lake'}`, highlight `view:tank` · 2. "Subtract the resting state" — "Hydrostatic balance
  takes care of ρ_s g; only the departure ρ′ remains — (4.84), $\rho\frac{D\mathbf u}{Dt}=-\nabla p'+\rho'\mathbf g+\mu\nabla^2\mathbf u$."
  `derive: {id: 'D27', step: 3}` · 3. "Keep it only with g" — "In the inertia ρ′ is a 0.2 % correction (grey); next to g it
  is the whole driving force (blue)." `terms: true`, `derive: {id: 'D27', step: 6}` · 4. "Watch it rise" — "Press ▶: the
  blob rises under buoyancy and spreads by diffusion, $\frac{DT}{Dt}=\kappa\nabla^2T$ (4.89)." `play: true` · 5. "Switch it off"
  — "Toggle buoyancy off: the blob only spreads." `set {buoy: false, t: 0}`, `play: true` · 6. "Your turn" — "Predict:
  will Boussinesq hold for a 10 km deep atmosphere? Choose the preset and read the status." `set {sc: 'deep_atmosphere',
  buoy: true}`, `play: false`.
- **Equations:** `pert` "Perturbation form" ref 'Eq. (4.84)' `\rho\frac{D\mathbf u}{Dt}=-\nabla p'+\rho'\mathbf g+\mu\nabla^2\mathbf u`
  · `bous` "Boussinesq momentum" ref 'Eq. (4.86)' `\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u`
  live "buoyancy = 0.0196 m/s²" · `heat88` "Heat with C_p" ref 'Eq. (4.88)' `\rho C_p\frac{DT}{Dt}=\rho\varepsilon-\nabla\cdot\mathbf q` ·
  `heat` "Boussinesq heat" ref 'Eq. (4.89)' `\frac{DT}{Dt}=\kappa\nabla^2T` · `eos` "Linear equation of state" ref '§4.9'
  `\rho=\rho_0\big[1-\alpha(T-T_0)\big]` · `valid` "Validity" ref '§4.9' `\alpha\,\delta T\ll1,\quad L\ll c^2/g`.
- **Check yourself:** (1) "Water, δT = 10 K: by what fraction does density change?" — "αδT = 2×10⁻³ (0.2 %)." `set {sc:
  'lake', dT: 10}` · (2) "Why is the grey 'dropped' bar so small?" — "It is ρ′/ρ₀ ≈ αδT times the inertia — 1000× smaller
  here." · (3) "Which scenario fails, and which condition?" — "The deep atmosphere: L ≈ c²/g (≈ 11.8 km for air)." `set
  {sc: 'deep_atmosphere'}` · (4) "Why does the heat equation have C_p, not C_v, although ∇·u ≈ 0?" — "The pressure is large,
  so p∇·u is not negligible: for a gas it equals ρ(C_p − C_v)DT/Dt and turns C_v into C_p (D28)."
- **Selftest parity rows:** `{name: 'alpha dT lake', js: validity(2e-4, 10, 10, 0.1, 1500, 9.81, 1e-6, 4186).alphaDT, py:
  'ch04.boussinesq_validity(2e-4, 10.0, 10.0, 0.1, c=1500.0, nu=1e-6, cp=4186.0)["alpha_dT"]', rtol: 1e-12}` · `{name:
  'g prime', js: validity(2e-4, 10, 10, 0.1, 1500, 9.81, 1e-6, 4186).gp, py: 'ch04.boussinesq_validity(2e-4, 10.0, 10.0,
  0.1, c=1500.0, nu=1e-6, cp=4186.0)["g_prime"]', rtol: 1e-12}` · `{name:
  'c2/g air', js: validity(1/288, 10, 1000, 1, 340.3, 9.81, 1.5e-5, 1005).Hc, py: 'ch04.boussinesq_validity(1/288, 10.0,
  1000.0, 1.0, c=340.3, nu=1.5e-5, cp=1005.0)["H_c"]', rtol: 1e-12}` · `{name: 'heating ratio', js: validity(2e-4, 1, 1, 0.1, 1500, 9.81, 1e-6,
  4186).heat, py: 'ch04.boussinesq_validity(2e-4, 1.0, 1.0, 0.1, c=1500.0, nu=1e-6, cp=4186.0)["heating_ratio"]', rtol:
  1e-12}` · `{name: 'blob T', js: blob(0.02, 0.01, 10, {U: [0.01, 0], kappa: 1e-4, sigma0: 0.05}),
  py: 'ch04.gaussian_blob_advection_diffusion([0.02, 0.01], 10.0, U=[0.01, 0.0], kappa=1e-4, sigma0=0.05)', rtol: 1e-12}` ·
  `{name: 'rise speed', js: rise(20, 0.0196, 10).w, py: 'ch04.blob_rise(20.0, 0.0196, 10.0)["w"]', rtol: 1e-12}` ·
  invariant `{name: 'buoyancy sign', js: buoy(-2, 1000), expect: 0.01962, rtol: 1e-12}`.
- **Fit plan:** 360×640: `tank` + `bars`; `validity` hidden (status shows the verdict and the largest ratio); scenario chips
  5 short labels; dT and L controls only. Landscape: `tank` + `bars` in a row. Desktop: three views.

### E9 · dynamic_similarity_models
- **Title:** "When does a model behave like the real thing?" · **Summary:** "A prototype and its model side by side, the
  dimensionless groups as paired bars with 'matched' badges, and the consequences: sphere data collapsing on one C_D(Re)
  curve, a ship's drag split and scaled, a rotating tank standing in for the atmosphere." · **CORE:** C15 (also N130,
  R12, N131, N133–N137, N139, N142, N147, N152, N156) · **Reference:** `stride_padding_playground.html` (presets that are
  the classic settings, the formula with numbers plugged in, badges explaining the result); links to ch01
  `buckingham_pi_machine`.
- **meta:** `viz:order 9` · `viz:sections 4.11` · `viz:equations 4.100 4.101 4.102 4.103 4.104 4.105 4.107 4.111` ·
  `viz:fluidpy ch04.reynolds_number ch04.froude_number ch04.strouhal_number ch04.richardson_number ch04.rossby_number
  ch04.mach_number ch04.sphere_drag_coefficient ch04.froude_scaled_speed ch04.model_prototype ch04.ship_drag_extrapolation`
  · `viz:derivations D30`.
- **Physics (JS ↔ Python):** `Re(U, l, nu)`, `Fr(U, l, g)`, `St(Om, l, U)`, `Ri(gp, l, U)`, `Ro(U, Om, l)` (= U/(2Ωl)),
  `Ma(U, c)` ↔ the `ch04` group functions · `cd(Re)` (Morrison 2013) ↔ `ch04.sphere_drag_coefficient` · `uFroude(Up, lp,
  lm)` ↔ `ch04.froude_scaled_speed` · `mp(lp, Up, scale, fp, fm, match)` ↔ `ch04.model_prototype` · `ship(Lp, Up, Sp, scale,
  Dm, cfm, cfp)` ↔ `ch04.ship_drag_extrapolation` · `ittc(Re)` = 0.075/(log₁₀Re − 2)² (our friction line, gloss).
- **Modes:** **sphere** · **ship** (Ex. 4.8) · **stratified** (Fr′, Ri: a sill flow and its lab model) · **rotating tank**
  (Rossby number U/(2Ωl), forward pointer to Ch. 13).
- **Views** (rows [1.1, 1.0]):
  1. `pair` "Prototype and model" — the two drawn to scale side by side (ship hulls, spheres, a stratified sill, the
     atmosphere vs a rotating tank) with their l, U, ν, g (and N or Ω); the model's size set by λ.
  2. `groups` "The groups" — paired bars (prototype grey outline, model filled) on a log axis for Re, Fr, St (and M, Ri or
     Ro by mode), each with a badge "✔ matched" (within 5 %) or "× 125 off"; title carries the verdict.
  3. `curve` "The consequence" (`hidePortrait: true`) — sphere: C_D vs Re (Morrison curve black, Stokes line grey dashed)
     with 40 synthetic "experiments" (dots coloured by fluid) and a toggle to plot them against U instead (the scatter
     returns); ship: the drag waterfall model total → friction + wave → wave × (ρ_p/ρ_m)λ⁻³ → + prototype friction =
     prototype total, next to the naive bar; stratified: the interface shape for the matched Fr′; rotating tank: the
     regime line Ro = 1 with both points.
  Portrait: `pair` + `groups`; `curve` hidden (its key number in the status and Explain).
- **Controls:** `mode` chips · `lam` "Scale λ = l_m/l_p" 10⁻³…1 (log), default 1/25 · `Um` "Model speed U_m" (log range by
  mode) with a "match Fr" / "match Re" / "match Ro" button that sets it · `fluid` chips "water · air · glycerine" (model
  fluid) · `axes` toggle "dimensional · dimensionless" (sphere).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (4) · **presets** "1:25 ship in water"
  {mode: 'ship', lam: 0.04, match: 'Fr'} · "1:10 car in air" {mode: 'sphere', shape: 'car', lp: 4, Up: 30, lam: 0.1, fluid: 'air',
  match: 'Re'} (a `shape` chip sphere/car changes only the drawing and the reference area) (needs 300 m/s: Mach 0.88 — the status warns that M no longer matches) · "sphere in glycerine vs air" {mode:
  'sphere', fluid: 'glycerine'} · "lab tank vs atmosphere" {mode: 'rotating', lam: 1e-6} · **status** ("✅ dynamically
  similar: Re, Fr matched" / "⚠️ Fr matched, Re off by ×125: correct the friction" / "⚠️ Re matched but M = 0.88: compressibility
  now differs") · **inspector** (click a data point: "d = 2 cm, U = 0.5 m/s, water ν = 10⁻⁶: Re = 0.5 × 0.02/10⁻⁶ = 10⁴ →
  C_D = 0.39").
- **Readouts:** "U_m" · "Re_m/Re_p" · "Fr_m/Fr_p" · "Prototype drag" (ship) · "C_D" (sphere).
- **Explain:**
  0. *What the views show* — "Left, the real thing and its model at scale λ. Middle, each dimensionless group of the two
     (grey outline = prototype, filled = model): a model is the real flow in miniature only if the groups that matter
     match. Right (hidden on phones), what follows: collapsed data or a corrected drag."
  1. *The equation behind it* — "(4.101): $\big[\frac{\Omega l}{U}\big]\frac{\partial\mathbf u^*}{\partial t^*}+(\mathbf u^*\cdot\nabla^*)\mathbf u^*=-\nabla^*p^*+\big[\frac{gl}{U^2}\big]\mathbf g^*+\big[\frac{\mu}{\rho Ul}\big]\nabla^{*2}\mathbf u^*$
     — three brackets: St, 1/Fr², 1/Re."
  2. *Prototype groups* — Re_p = U_pl_p/ν_p = … ; Fr_p = U_p/√(gl_p) = … ; (St, M, Ri or Ro as the mode needs).
  3. *Model groups* — with λ = … and U_m = … in … : Re_m = …, Fr_m = …. (boxed ratios)
  4. *What matching demands* — "Froude: $U_m=U_p\sqrt{\lambda}=$ **…** m/s; then Re_m/Re_p = λ^{3/2}(ν_p/ν_m) = **…**. To
     match Re as well you would need ν_m = ν_pλ^{3/2} = … m²/s — no practical liquid has it."
  5. *The consequence* (per mode) — ship: "friction C_f = 0.075/(log₁₀Re − 2)²: model … , prototype …; wave drag ×
     (ρ_p/ρ_m)λ⁻³ = … kN; total **…** kN (naive … kN)"; sphere: "C_D(Re = …) = **…**: every sphere at that Re has it";
     rotating: "Ro = U/(2Ωl) = … for both — the tank's eddies turn like the atmosphere's".
  6. *The consequence view* (hint on phones).
  7. *Reading the current setting* — Fr matched, Re not: "The waves are right, the friction is not — Froude's
     hypothesis: scale the two separately." · all matched: "Dynamically similar: every dimensionless result carries over
     unchanged." · Re matched in air at high speed: "Matching Re made the model fast enough that compressibility (M)
     now differs: another group woke up." · sphere, dimensionless axes: "Many sizes, speeds and fluids — one curve: the
     equations only know Re."
- **Derivation tab:** **D30** (8 steps) `view: 'groups'`; goal `set {mode: 'sphere'}`; steps 2–6 each turn one
  dimensional term into a bracket and light its bar (step 3 unsteady → St blue, step 5 gravity → 1/Fr² grey, step 6
  viscous → 1/Re rose); step 7 **live** "St = Ωl/U = …, gl/U² = …, μ/ρUl = …"; step 8 `watch: "equal brackets ⇒ same
  equation ⇒ same solution"`; interpret "your prototype's brackets: …; your model's: …".
- **Code:**
  ```python
  lp, Up, lam = {{lp}}, {{Up}}, {{lam}}           # prototype size [m], speed [m/s], scale l_m/l_p
  mp = ch04.model_prototype(lp, Up, lam, fluid_p="{{fp}}", fluid_m="{{fm}}", match="{{match}}")
  print(mp["U_m"], mp["Re_ratio"], mp["matched"])  # {{Um}} m/s, Re_m/Re_p = {{rr}}, {{matched}}
  print(ch04.reynolds_number({{Um}}, lp*lam, {{num}}), ch04.froude_number({{Um}}, lp*lam))
  print(ch04.sphere_drag_coefficient({{Re}}))    # C_D = {{cd}}
  d = ch04.ship_drag_extrapolation(lp, Up, 2500.0, lam, 40.0, {{cfm}}, {{cfp}})
  print(d["D_p_total"])                           # {{Dp}} N
  ```
- **Walkthrough (6 steps):** 1. "A small copy" — "A 4-m model of a 100-m ship in a towing tank: what must be the same?"
  `set {mode: 'ship', lam: 0.04}`, highlight `view:pair` · 2. "The equation knows only brackets" — "Scale everything by
  l, U, Ω: [Ωl/U]∂u*/∂t* + (u*·∇*)u* = −∇*p* + [gl/U²]g* + [μ/ρUl]∇*²u* (4.101) — three brackets remain." `derive: {id: 'D30', step: 7}`, `eq: 'nd'` · 3. "Match the waves" —
  "Press 'match Fr': U_m = U_p√λ = 2 m/s." `set {match: 'Fr'}`, `readouts: ['Um']` · 4. "Reynolds is lost" — "Now
  Re_m/Re_p = λ^{3/2} = 1/125: the model's boundary layers are wrong. Correct the friction separately." highlight
  `view:groups` · 5. "One curve for all spheres" — "Sphere mode: forty experiments, three fluids — one curve when plotted
  against Re." `set {mode: 'sphere', axes: 'dimensionless'}` · 6. "Your turn" — "Predict: to match Re for a 1:10 car in
  air, how fast must the model go? Is that a good idea?" `set {mode: 'sphere', lam: 0.1, fluid: 'air'}`.
- **Equations:** `scal` "Scalings" ref 'Eq. (4.100)' `x_i^*=\frac{x_i}{l},\ t^*=\Omega t,\ u_j^*=\frac{u_j}{U},\ p^*=\frac{p-p_\infty}{\rho U^2},\ g_j^*=\frac{g_j}g` ·
  `nd` "Dimensionless Navier–Stokes" ref 'Eq. (4.101)'
  `\begin{aligned}&\Big[\frac{\Omega l}{U}\Big]\frac{\partial\mathbf u^*}{\partial t^*}+(\mathbf u^*\cdot\nabla^*)\mathbf u^*\\&=-\nabla^*p^*+\Big[\frac{gl}{U^2}\Big]\mathbf g^*+\Big[\frac{\mu}{\rho Ul}\Big]\nabla^{*2}\mathbf u^*\end{aligned}` ·
  `re` "Reynolds" ref 'Eq. (4.103)' `\mathrm{Re}=\frac{\rho Ul}{\mu}` live "= 1.0×10⁹" · `fr` "Froude" ref 'Eq. (4.104)'
  `\mathrm{Fr}=\frac U{\sqrt{gl}}` live "= 0.319" · `ri` "Richardson" ref 'Eq. (4.105)' `\mathrm{Ri}=\frac1{\mathrm{Fr}'^2}=\frac{g'l}{U^2}` ·
  `cd` "Drag coefficient" ref 'Eq. (4.107)' `C_D=\frac{F_D}{\tfrac12\rho U^2A}` · `ma` "Mach" ref 'Eq. (4.111)' `M=U/c`.
- **Check yourself:** (1) "Froude scaling at λ = 1/25: model speed for a 10 m/s ship?" — "10/√25 = 2 m/s." `set {mode:
  'ship', lam: 0.04, match: 'Fr'}` · (2) "Then how far off is Re?" — "Re_m/Re_p = λ^{3/2} = 1/125 (same ν)." · (3) "A 1-cm
  sphere at 1 m/s in water and a 1-cm sphere in air: what air speed gives the same C_D?" — "Equal Re: U_air = 1 ×
  (1.5×10⁻⁵/10⁻⁶) = 15 m/s." `set {mode: 'sphere', fluid: 'air'}` · (4) "Why can't a towing tank use a fluid that fixes
  Re?" — "It would need ν_m = ν_pλ^{3/2} ≈ 8×10⁻⁹ m²/s, 125 times less viscous than water — no such liquid."
- **Selftest parity rows:** `{name: 'Re ship', js: Re(10, 100, 1e-6), py: 'ch04.reynolds_number(10.0, 100.0, 1e-6)', rtol:
  1e-12}` · `{name: 'Fr ship', js: Fr(10, 100, 9.81), py: 'ch04.froude_number(10.0, 100.0)', rtol: 1e-12}` · `{name:
  'Froude speed', js: uFroude(10, 100, 4), py: 'ch04.froude_scaled_speed(10.0, 100.0, 4.0)', rtol: 1e-12}` · `{name: 'CD
  Re=1e4', js: cd(1e4), py: 'ch04.sphere_drag_coefficient(10000.0)', rtol: 1e-12}` · `{name: 'Re ratio 1:25', js: mp(100,
  10, 0.04, 'water', 'water', 'Fr').ReRatio, py: 'ch04.model_prototype(100.0, 10.0, 0.04)["Re_ratio"]', rtol: 1e-10}` ·
  `{name: 'ship total', js: ship(100, 10, 2500, 0.04, 40, 0.00312, 0.00153).total, py: 'ch04.ship_drag_extrapolation(100.0,
  10.0, 2500.0, 0.04, 40.0, 0.00312, 0.00153)["D_p_total"]', rtol: 1e-10}` · `{name: 'Rossby atmosphere', js: Ro(10,
  7.292115e-5, 1e6), py: 'ch04.rossby_number(10.0, 7.292115e-5, 1.0e6)', rtol: 1e-12}` · `{name: 'Richardson', js:
  Ri(9.81e-3, 100, 0.1), py: 'ch04.richardson_number(9.81e-3, 100.0, 0.1)', rtol: 1e-12}`.
- **Fit plan:** 360×640: `pair` + `groups`; `curve` hidden (status names the mismatch, Explain §5 has the consequence
  numbers); mode chips 4 short labels ("sphere · ship · layers · tank"); λ and U_m only. Landscape: `pair` + `groups` in a
  row. Desktop: `pair` top-left, `groups` bottom-left, `curve` right spanning both rows.

### B1 · kinematic_free_surface
(Backup — built only if one of E1–E9 fails review; not embedded otherwise.)
- **Title:** "Does the sea surface move with the water on it?" · **Summary:** "A wave surface with particles riding on it:
  the fluid's normal velocity equals the surface's, Dη/Dt = 0; a 'leaky' interface (evaporation, a shock) makes the
  relative normal velocity (4.92) non-zero and particles cross." · **CORE:** C14 (also N116, N117, N118, N119, N120) ·
  **Reference:** `forward_noising_lab.html` (click a point to see its exact arithmetic).
- **meta:** `viz:order 10` · `viz:sections 4.10` · `viz:equations 4.90 4.91 4.92 4.93` · `viz:fluidpy
  ch04.linear_wave_surface ch04.wave_kinematic_residual ch04.surface_normal_speed ch04.kinematic_bc_residual` ·
  `viz:derivations D29`.
- **Physics:** `wave(x, z, t, a, k, H)` ↔ `ch04.linear_wave_surface` · `residual(x, t, a, k, full)` ↔
  `ch04.wave_kinematic_residual` · `normalSpeed(eta, x, t)` ↔ `ch04.surface_normal_speed("moving_wall"/"linear_wave", …)`.
- **Modes:** free surface · moving wall · leaky interface (a prescribed mass flux m″: (u_rel)_n = m″/ρ) · pillbox (N116:
  a cylinder shrinking across an interface with flux bars converging to the jump condition).
- **Views:** `surface` (surface, particles, normal arrows u·n teal and u_s·n rose at a probe) · `residual` (log–log
  residual/c vs ka for the full condition, slope 2, and 0 for the linearised one; pillbox mode: side and volume terms vs l,
  slope 1).
- **Controls:** `ka` 0.01…0.3 (log) · wavelength λ 1…100 m · mode chips · transport t (one period, loop) · `mflux` (leaky).
- **Depth features:** Explain, Code, Derivation + **linked views** (2) · **transport** · **status** ("🌊 Dη/Dt = 0: no fluid
  crosses the surface" / "💨 leaky: (u_rel)_n = … m/s") · **inspector** (click a surface point: "Dη/Dt = ∂η/∂t + u·∇η =
  … + … = …").
- **Explain:** 0 views; 1 the surface and its normal speed −(∂η/∂t)/|∇η|; 2 the fluid's u·n at the same point; 3 the
  residual (4.91) and why the linearised version is exact for the linear solution; 4 (4.92) and (4.93) for the leaky
  mode; 5 reading: "a free surface is made of fluid particles, so it moves exactly as fast as the fluid pushes it — no
  faster, no slower."
- **Derivation tab:** D29 (8 steps, `view: 'surface'`, step 6 live "u_s·n = −(∂η/∂t)/|∇η| = …").
- **Code:** `ch04.linear_wave_surface(x, z, t, a, k)`, `ch04.wave_kinematic_residual(x, t, a, k, full=True)` with live
  values.
- **Walkthrough (5 steps):** the question → particles ride the surface → the normal speeds match (4.91) → linearising
  costs (ka)² → your turn (leaky interface: predict which way particles cross).
- **Equations:** (4.90) $d\eta/dt=\partial\eta/\partial t+(\mathbf u_s\cdot\nabla)\eta=0$ · (4.91) $D\eta/Dt=0$ on η = 0 · (4.92)
  $(u_{rel})_n=(1/\lvert\nabla\eta\rvert)D\eta/Dt$ · (4.93) $(\rho/\lvert\nabla\eta\rvert)D\eta/Dt$.
- **Check yourself:** (1) "Does (4.91) fix the tangential velocity at the surface?" — "No — only the normal part; the
  tangential one comes from the stress condition." (2) "ka halved: the full-condition residual (in units of c)?" — "A
  quarter: it scales as (ka)²." (3) "A moving wall η = x − Vt: which fluid velocity satisfies (4.91)?" — "Any with u·n = V."
- **Selftest parity rows:** `{name: 'wave residual ka=0.1', js: residual(0.3, 0, 0.1, 1, true), py:
  'ch04.wave_kinematic_residual(0.3, 0.0, 0.1, 1.0)', rtol: 1e-8}` · `{name: 'moving wall speed', js:
  normalSpeed('moving_wall', [2, 0.3, 0], 1, {V: 1.5}), py: 'ch04.surface_normal_speed("moving_wall", [2.0, 0.3, 0.0], 1.0,
  V=1.5)', rtol: 1e-8}`.
- **Fit plan:** two views stacked on portrait, side by side on landscape; ≤ 4 controls.

---

## Part D — runtime budget (full run < 5 min on a laptop / Colab CPU)

Chapter 4 is mostly closed forms and point stencils (µs each); the costs are the animations (≈ 0.12 s per frame at dpi 80,
ch01–ch03 timings), the plotly slider figures that precompute every step, the quadratures (`quad` ms, `tplquad` ≈ 50 ms,
`dblquad` ≈ 10 ms), `solve_ivp` (ms), the 3-D control-volume quadrature (24³ nodes ≈ 20 ms), and the sympy cells (the two
★★★ checks D09 and D15, Lamb identity, Appendix-B pump equations, the dimensionless NS). The ch03 notebook ran in 64 s; this
one is larger (15 CORE, 30 derivations, 9 explainers) and must still stay under 5 min.

| Section | Heaviest cells | Full | FAST (`FLUIDPY_FAST=1`) |
|---|---|---|---|
| setup + imports | numpy/scipy/sympy/plotly/pint/pandas | 9 s | 9 s |
| §4.1 | ledger table (pandas) | 0.5 s | 0.5 s |
| §4.2 C01 | R01 path lines (2 × `solve_ivp`), interval budgets, trapezoid sums (2001 points × 6), 3-D `mass_budget` on a box (24³ nodes, twice), frames player 24 frames (1-D strip + bars) | 6 s | 4 s (12 frames, 16³ nodes) |
| §4.2 C02 | R02 `divergence_theorem_box` (24³), `ball_integral` × 3 (`tplquad`, 0.15 s), periodic grid 64² divergence, product-rule sympy, stencils, sympy residual, slider 21 steps × 3 traces × 200 points, stratified figure (60 s of particles, closed form) | 6 s | 4 s (11 steps, 48² grid) |
| §4.3 C03 | sympy ∇·(∇χ × ∇ψ), stencils, flux quadrature (400 points), cylinder heat map 241 × 161 + 20 contours, slider 17 steps (contours precomputed with `matplotlib`'s contour generator → traces, ≈ 0.05 s each), plotly 3-D (4 surfaces 40 × 40 + 4 streamlines × 400 points), `stream_tube_mass_flux` (`curl_flux` 64²) | 9 s | 6 s (9 steps, 24² surfaces) |
| §4.4 C04 | wake quadratures, face budgets, figure, slider 20 steps × 2 traces × 400 points, bore animation 60 frames, rocket `solve_ivp`, sprinkler | 14 s | 8 s (30 frames, 10 steps) |
| §4.4 C05 | sympy series (Ex. 4.2 re-run ≈ 1 s), tank Euler loop (900 steps), cylinder C_p, tank figure | 4 s | 4 s |
| §4.4 C06 | sympy conservative→advective (3 components, generic functions ≈ 1.5 s), Cauchy stencils at 4 points, term-bar figure | 4 s | 4 s |
| §4.5 C07 | cube log–log, isotropy test (50 rotations of a 3⁴ tensor), **D09 sympy ★★★** (K from δ's, 81 contractions, rotation spot-check of 20 components ≈ 2 s), stress-vs-angle figure, μ_v slider 11 steps | 7 s | 6 s |
| §4.6 C08 | Poiseuille stencils, `navier_stokes_sym` with μ(y) (≈ 1 s), curl-of-curl sympy, erfc, Stokes slider 20 steps × 3 traces × 300 points, Lamb–Oseen viscous-force quiver (`viscous_force_forms` on a 21² grid ≈ 0.4 s) | 7 s | 5 s (10 steps, 15² quiver) |
| §4.7 C09 | basis-rate check, cross-product sympy, **D15 sympy ★★★** (x = X + R(t)x′ differentiated twice, 3 components ≈ 3 s), frame terms, projectile, pump equations (Appendix-B sympy ≈ 3 s, cached), projectile animation 60 frames (two panels), high/low quiver, plotly 3-D Earth (sphere 40 × 40, 12 cones) | 20 s | 12 s (30 frames, 24 × 24 sphere) |
| §4.8 C10 | dissipation of 1000 random gradients (vectorised), Couette heating (closed form), three-panel figure, U slider 20 steps × 2 traces × 200 points, sympy chain-rule and square-completion lines | 5 s | 4 s (10 steps) |
| §4.9 C11 | Lamb identity stencils, B at 1000 random points (vectorised cylinder field), Rankine figure, σ slider 19 steps × 2 traces × 300 points | 4 s | 3 s (10 steps) |
| §4.9 C12 | U-tube and sphere closed forms, gauge sympy, U-tube animation 60 frames, sphere figure | 9 s | 5 s (30 frames) |
| §4.9 C13 | hydrostatic `integrate_hydrostatic` (5 points), validity table (5 scenarios), blob residual stencils, blob animation 60 frames on a 121 × 121 grid (closed form per frame), validity figure | 10 s | 6 s (30 frames, 81² grid) |
| §4.10 C14 | pillbox, stencil residuals, wave residual sweep (3 × 64 points), cap forces (`dblquad` + `quad` round the ellipse, ≈ 20 ms), meniscus closed form + `solve_ivp` (Radau, 4 angles ≈ 0.3 s), wave animation 60 frames (two particle rows), two figures | 12 s | 7 s (30 frames) |
| §4.11 C15 | `pi_groups` × 2, sympy substitution of (4.100) into (4.39b) (≈ 1 s), `nondimensional_ns_coefficients` (cached), group table, Stokes-layer collapse figure, C_D slider 36 steps × 3 traces × 400 points (Morrison vectorised), 60 synthetic experiments, ship waterfall | 7 s | 5 s (18 steps) |
| explainers (9 × `show_viz`) | read the HTML files | 1 s | 1 s |
| **Total** | | **≈ 135 s** | **≈ 94 s** |

**FAST plan.** Every size-dependent choice is written `a if not FAST else b` in the cell: animations 60 frames → 30, frames
players 24 → 12; slider figures ≤ 36 steps → ≤ 18, ≤ 4 traces × ≤ 400 points (each < 250 kB); 3-D control-volume
quadrature 24³ → 16³; heat-map grids 241 × 161 → 121 × 81; plotly 3-D surfaces 40² → 24²; quiver grids 21² → 15²; the
meniscus ODE for 4 → 2 angles. **Cached arrays / results:** `rho`, `u` of the expanding flow (C01) are reused by C02;
`u_cyl = ch03.cylinder_velocity_field(1.0, 1.0)` is built once in C03 and reused in C05 (C_p), C08 (Euler check), C11
(B at 1000 points) and C12 (viscous-irrotational residual); `ch04.exact_field("poiseuille", …)` from C08 is reused in C13
(N107 demo); the sympy results of D09, D15, `rotating_pump_equations()` and `nondimensional_ns_coefficients()` are cached
with `functools.lru_cache` in the library, so the notebook's second call is free; the Morrison curve (400 points) is
computed once in C15 for both the slider and the synthetic-experiment figure. **Outputs:** 5 videos (bore, projectile,
U-tube, blob, wave; < 2 MB each), 1 frames player (C01, < 3 MB), 9 slider figures, 2 plotly 3-D figures (stream surfaces,
Earth), ≈ 25 static figures — page well under 15 MB. Sympy cells never `simplify` expressions with more than two generic
functions at once (conservative→advective runs one component at a time; D15 simplifies each component separately after
substituting the rotation matrix).

---

## Part E — prerequisite ledger
Every concept, symbol, maths tool and Python function or idiom the notebook or its explainers use, with where it is
explained. "primer (in Cxx)" = a 📎 primer placed in that block before first use (the primer term is the Concept text
before any parenthesis — the builder uses it verbatim in `nb.primer`); "knowledge/primers.md: <term> (chNN Pnn)" = a
one-line reminder naming the earlier primer; a CORE/RECAP id alone = taught there; "Cxx (Nnn note)" = the NOTE of this
chapter placed in that block; "Cxx (gloss …)" = one sentence where it is used. Earlier chapters' material that is
neither a primer nor a ch04 RECAP is named by section ("Ch. 1 §1.5"), never by an earlier chapter's C-number (convention
13). Explainer-only items are listed with the explainer's host block.

| Concept | First used in | Explained by |
|---|---|---|
| integral (control-volume) vs differential (point) form of a law | §4.1 | C01 (N01 note in the §4.1 section) |
| equations vs unknowns; closure of the system | §4.1 | C06 (N02 note; ledger rows at C06, C08, C10) |
| pandas DataFrame as a labelled table | §4.1 | C06 (gloss in the §4.1 ledger cell's comment) |
| material volume and material surface (b = u) | C01 | R01 |
| Reynolds transport theorem (3.35) | C01 | R01 (text of the recap, Ch. 3 §3.6) |
| control volume, control surface velocity b | C01 | C01 |
| signed b·n and the swept volume of a moving wall | C01 (D01) | knowledge/primers.md: signed swept volume of a moving surface (ch03 P110) — reminder |
| relative velocity u − b through a moving wall | C01 | C01 (the idea sketch and D01 step 8) |
| instantaneously coincident control volume | C01 (D01) | C01 (N06 note, D01 steps 4–6) |
| mass of a material volume is constant (4.1) | C01 | C01 (N03 note) |
| mass budget of any control volume (4.5) | C01 | C01 |
| expanding test flow u = ax/(1 + at), ρ = ρ₀/(1 + at) | C01 | C01 (N03 note and worked example) |
| dataclasses and named results | C01 | primer (in C01) |
| path lines by solve_ivp | C01 (R01 code) | knowledge/primers.md: solve_ivp options (ch03 P94) — reminder |
| definite integral and np.trapezoid | C01 | knowledge/primers.md: trapezoid rule (ch01 P37) — reminder |
| central finite differences in time | C01 | knowledge/primers.md: finite differences (ch01 P21) — reminder |
| differentiation under the integral sign | C01 (D01) | knowledge/primers.md: differentiation under the integral sign (ch03 P109) — reminder |
| functions as arguments and lambda | C01 | knowledge/primers.md: functions as arguments and lambda (ch01 P29) — reminder |
| assert np.allclose | C01 | knowledge/primers.md: assert np.allclose (ch01 P15) — reminder |
| animate and show_animation | C01 | knowledge/primers.md: animate and show_animation (ch01 P16) — reminder |
| f-strings | C01 | knowledge/primers.md: f-strings (ch01 P04) — reminder |
| Gauss' divergence theorem (2.30) | C02 | R02 |
| continuity of a function and the small-ball argument | C02 (D02) | primer (in C02) |
| scipy.integrate.tplquad | C02 (primer code) | C02 (gloss in the primer's code comment; extends `quad`/`dblquad`) |
| scipy.integrate.quad and dblquad | C02 | knowledge/primers.md: scipy.integrate.quad and dblquad (ch03 P87) — reminder |
| localisation lemma | C02 | C02 (N07 note, D02 steps 5–7) |
| continuity equation (4.7) | C02 | C02 |
| flux divergence (transport) term | C02 | C02 (N08 note) |
| periodic grids and numpy stencils | C02 (N08 code) | knowledge/primers.md: numpy broadcasting (ch02 P77) — reminder |
| np.meshgrid and the project grid layout | C02 | knowledge/primers.md: np.meshgrid and the project grid layout (ch02 P76) — reminder |
| product rule for a divergence | C02 (D03) | primer (in C02) |
| material derivative D/Dt (3.5) | C02 (D03) | C02 (D03 step 2's *why*; Ch. 3 §3.2) |
| continuity with D/Dt (4.8) | C02 | C02 (N09 note) |
| incompressible flow Dρ/Dt = 0 (4.9) | C02 | C02 (N10 note) |
| ∇·u = 0 (4.10) and the volumetric strain rate | C02 | C02 (N11 note; Ch. 3 §3.4) |
| incompressible flow vs incompressible fluid; Mach < 0.3 | C02 | C02 (N12 note) |
| speed of sound c | C02 (N12) | C02 (N12 note; Ch. 1 §1.9) |
| sympy symbols, Function, diff, simplify | C02 | knowledge/primers.md: sympy (ch01 P40) — reminder |
| partial derivative | C02 | knowledge/primers.md: partial derivative (ch01 P25) — reminder |
| slider_figure | C02 | knowledge/primers.md: slider_figure (ch01 P17) — reminder |
| stratified flow ρ(z) | C02 (N10) | C02 (N10 note and figure) |
| steady continuity (4.11) | C03 | C03 (N13 note) |
| ∇·(∇×A) = 0 and ∇×(∇φ) = 0 | C03 (N14, D04) | C03 (gloss in N14 and D04 step 2's *why*; Ch. 2 §2.13) |
| curl of a product ∇×(χ∇ψ) | C03 (D04) | C03 (gloss in D04 step 3's *why*, sympy check in N14's cell) |
| vector potential and two stream functions (4.12) | C03 | C03 (N14 note) |
| cross product with e_z; right-hand rule | C03 (D04) | knowledge/primers.md: right-hand rule and orientation (ch02 P74) — reminder |
| 2-D stream function ψ, ρu = ∂ψ/∂y, ρv = −∂ψ/∂x | C03 | C03 |
| flux between streamlines = ψ₂ − ψ₁ | C03 | C03 (D04 steps 6–8) |
| chain rule along a curve dψ = ∇ψ·ds | C03 (D04) | knowledge/primers.md: multivariable chain rule along a path (ch03 P91) — reminder |
| line integral along a path | C03 (D04) | knowledge/primers.md: line integral of a vector field around a loop (ch02 P86) — reminder |
| stream surfaces; 3-D streamlines as intersections | C03 | C03 (N15 note) |
| plotly 3-D surfaces and lines | C03 | knowledge/primers.md: plotly 3-D arrows, lines and meshes (ch02 P64) — reminder |
| stream-tube flux (b − a)(d − c) and Stokes' theorem | C03 | C03 (N16 note; Ch. 2 §2.13) |
| axisymmetric (Stokes) stream function | C03 | C03 (N17 note) |
| cylindrical unit vectors | C03 (N17) | knowledge/primers.md: cylindrical and spherical unit vectors (ch03 P88) — reminder |
| constant-density ψ as volume flux | C03 | C03 (N18 note) |
| GFD sign convention for ψ | C03 | C03 (⚠️ callout) |
| plt.contour and streamplot | C03 | knowledge/primers.md: plt.contour, plt.quiver and plt.streamplot (ch02 P78) — reminder |
| stagnation point | C03 (worked example) | C05 (N104 note; first met as a label in C03's figure) |
| Newton's second law for a material volume (4.13) | C04 | C04 (N19 note) |
| Newton's second law and momentum | C04 | knowledge/primers.md: Newton's second law and momentum (ch01 P09) — reminder |
| body forces vs surface forces | C04 | C04 (N26 note) |
| momentum flux through a surface | C04 | primer (in C04) |
| RTT applied component by component to a vector | C04 (D05) | C04 (gloss in D05 step 2's *why*) |
| momentum budget of any control volume (4.17) | C04 | C04 |
| book slip in (4.15) | C04 | C04 (N21 note) |
| conservative force and its potential | C04 | primer (in C04) |
| gravity potential Φ = gz, g = −∇Φ (4.18) | C04 | C04 (N27 note) |
| surface force from the stress tensor f_j = n_iτ_ij | C04 | R03 |
| surface tension enters through boundary conditions | C04 | R04 |
| Newton's third law; drag on the body vs force on the fluid | C04 | C04 (⚠️ callout N28) |
| Gaussian integral ∫e^{−y²/b²}dy = b√π | C04 (worked example) | C04 (gloss in the worked example) |
| error function erf (finite box) | C04 (E1) | C04 (gloss in the E1 explainer cell's text; `Viz.num.erf`) |
| wake drag from a velocity deficit | C04 | C04 (N29 note and worked example) |
| momentum thickness (named) | C04 | C04 (N29 note, pointer to Ch. 9) |
| bore / hydraulic jump; shallow-water speed √(gh) | C04 | C04 (N31 note) |
| hydrostatic pressure on a vertical face ½ρgh² | C04 (N31) | knowledge/primers.md: net force from pressure (ch01 P28) — reminder |
| dimensional analysis giving one group | C04 (N31) | R12 (Π theorem; reminded in C15, used in C04 in one line) |
| rocket with an accelerating control volume | C04 | C04 (N32 note) |
| Tsiolkovsky rocket equation | C04 | C04 (gloss in N32) |
| natural logarithm | C04 (N32) | knowledge/primers.md: natural logarithm and exponential (ch01 P36) — reminder |
| torque, moment arm and moment of inertia | C04 | primer (in C04) |
| np.cross | C04 (primer code) | C04 (gloss in the torque primer's code comment) |
| angular momentum dH/dt = M (4.64) and its CV form (4.65) | C04 | C04 (N82, N83 notes) |
| lawn sprinkler torque | C04 | C04 (N84 note) |
| Bernoulli along a streamline (4.19) | C05 | C05 |
| sympy expand, series, removeO, collect and subs | C05 (D06) | primer (in C05) |
| first-order Taylor expansion | C05 (D06) | knowledge/primers.md: first-order Taylor expansion (ch01 P26) — reminder |
| limits and orders of smallness | C05 (D06) | knowledge/primers.md: limits and orders of smallness (ch02 P68) — reminder |
| side pressure force on a slowly widening tube | C05 (D06) | C05 (gloss in D06 step 5's *why*) |
| stagnation and dynamic pressure | C05 | C05 (N104 note) |
| pitot tube | C05 | C05 (N103 note) |
| orifice, Torricelli's speed, vena contracta | C05 | C05 (N105 note; vena contracta glossed there) |
| explicit Euler stepping (tank draining) | C05 (from scratch) | knowledge/primers.md: explicit stepping (ch01 P30) — reminder |
| pressure coefficient on the cylinder C_p = 1 − 4 sin²θ | C05 (figure) | C15 (N139 note; used first in C05's figure with the formula written) |
| Cauchy's equation of motion (4.24) | C06 | C06 |
| tensor divergence over the first index | C06 | primer (in C06) |
| Gauss applied to each component of a tensor (4.20a, b) | C06 (D07) | C06 (N33, N34 notes) |
| conservative (flux) form (4.22) | C06 | C06 (N36 note) |
| the bracket is continuity (4.23) | C06 | C06 (N37 note) |
| stress tensor τ_ij, first index = face | C06 | knowledge/primers.md: stress (ch01 P05) — reminder, with R03 |
| rigid (solid-body) rotation and its pressure | C06 (worked example) | knowledge/primers.md: rigid-body velocity Ω × x (ch03 P101) — reminder |
| centripetal acceleration −Ω²r | C06 | C06 (worked example) |
| Galilean invariance ⇒ stress depends on S only | C07 | R05 |
| non-Newtonian fluids | C07 | R06 |
| constitutive equation | C07 | C07 (N39 note) |
| stress symmetry (4.25) | C07 | C07 (N38 note, D08) |
| log–log slope and observed_order | C07 | knowledge/primers.md: power laws and log–log plots (ch01 P13) — reminder |
| np.logspace | C07 | knowledge/primers.md: np.linspace and np.logspace (ch01 P06) — reminder |
| static stress −pδ_ij (4.26); isotropic 2nd-order tensor | C07 | C07 (N40 note; Ch. 2 §2.5) |
| viscous (deviatoric) stress σ_ij (4.27) | C07 | C07 (N41 note) |
| fourth-order tensor K_ijmn (4.28) | C07 | C07 (N42 note) |
| np.einsum index strings | C07 | knowledge/primers.md: np.einsum index strings (ch02 P62) — reminder |
| isotropic fourth-order tensor | C07 | primer (in C07) |
| random rotations and transform_tensor | C07 (primer code) | knowledge/primers.md: np.linalg.norm and np.linalg.qr (ch02 P67) — reminder |
| δ substitution in products δ_imδ_jnS_mn = S_ij | C07 (D09) | C07 (gloss in D09 steps 5–7's *why*; δ from Ch. 2 §2.1) |
| γ = μ (4.30) and the μ + γ subtlety | C07 | C07 (N44 note, D09 steps 8–9) |
| Newtonian stress law (4.31) | C07 | C07 |
| deviatoric (traceless) part of a tensor | C07 (D10) | primer (in C07) |
| bulk-viscosity form (4.37); Newton's law τ = μ du/dy (1.3) | C07 | C07 (N51 note; Ch. 1 §1.5) |
| trace, δ_ii = 3; thermodynamic vs mean pressure (4.32)–(4.34) | C07 | C07 (N45–N47 notes) |
| incompressible Newtonian stress (4.35) | C07 | C07 (N48 note) |
| bulk viscosity μ_v and Stokes' assumption (4.36) | C07 | C07 (N49, N50 notes) |
| shear rate γ̇ and extension rate s (symbols) | C07 | C07 (notation table in §4.1 and the worked example) |
| traction on a plane; normal and shear parts | C07 (figure) | R03 |
| quadratic form n·τ·n; principal directions | C07 (figure) | knowledge/primers.md: quadratic form and the Rayleigh quotient (ch02 P82) — reminder |
| Mohr's circle (named) | C07 (figure) | C07 (gloss in the figure's reading notes; Ch. 2 §2.11) |
| Navier–Stokes with variable viscosity (4.38) | C08 | C08 (N52 note, D11) |
| barotropic flow | C08 (N52) | C11 (N87 note; named first in C08's ledger note in one sentence) |
| Schwarz's theorem | C08 (D12) | primer (in C08) |
| vector Laplacian in Cartesian components | C08 (D12) | C08 (gloss in D12 step 5's *why*) |
| constant-viscosity compressible form (4.39a) | C08 | C08 (N53 note) |
| incompressible Navier–Stokes (4.39b) | C08 | C08 |
| curl of a curl identity | C08 (D13) | primer (in C08) |
| ε–δ identity (2.19) | C08 (D13) | C08 (primer text; Ch. 2 §2.7) |
| vorticity ω = ∇×u | C08 (D13) | C08 (D13 step 3's *why*; Ch. 3 §3.4) |
| viscous force three ways (4.40); the paradox | C08 | C08 (N54 note) |
| Euler's equation (4.41) | C08 | C08 (N55 note) |
| plane Poiseuille and Couette flow | C08 | C08 (worked example; Couette first in Ch. 1 §1.5) |
| complementary error function erfc | C08 | primer (in C08) |
| Stokes' first problem | C08 | C08 (plotly figure text; derived in Ch. 8) |
| Lamb–Oseen (Gaussian) vortex, σ² = 4νt | C08 (figure) | C08 (figure text; the Gaussian vortex of Ch. 3 §3.5) |
| Taylor–Green vortex | C08 (E4) | C08 (E4's Explain §1 formula, named in the explainer cell) |
| local Reynolds number (speed × l/ν) | C08 (E4) | C15 (N135 note; used first in E4's Explain §4 with its formula) |
| inertial frame | C09 | C09 (N56 note) |
| four meanings of the prime | C09 | C09 (⚠️ callout; notation table in §4.1) |
| frames of reference and relative velocity | C09 | knowledge/primers.md: frames of reference and relative velocity (ch03 P96) — reminder |
| rotating frame of reference (first look) | C09 | knowledge/primers.md: rotating frame of reference (ch03 P103) — reminder |
| derivative of a rotating unit vector | C09 (D14) | primer (in C09) |
| Rodrigues rotation of a basis | C09 (primer code) | C09 (gloss in the primer's code comment; Ch. 2 §2.2) |
| product rule for a cross product | C09 (D15) | primer (in C09) |
| product rule | C09 (D14) | knowledge/primers.md: product rule for differentials (ch01 P38) — reminder |
| vector triple product a × (b × c) = b(a·c) − c(a·b) | C09 (D15, D18) | C09 (gloss in D15 step 10's *why*; Ch. 2 §2.7) |
| velocity in a rotating frame (4.42) | C09 | C09 (N57 note, D14) |
| acceleration in a noninertial frame (4.43), (4.44) | C09 | C09 (N58, N59 notes, D15) |
| Navier–Stokes in a noninertial frame (4.45) | C09 | C09 |
| Coriolis term vs Coriolis force; centripetal vs centrifugal | C09 | C09 (⚠️ callout) |
| frame acceleration and weightlessness | C09 | C09 (N60 note) |
| latitude, Earth's rotation rate and the local vertical | C09 | primer (in C09) |
| Coriolis parameter f = 2Ω sin φ (named) | C09 | C09 (primer text, forward pointer to Ch. 13) |
| constant-acceleration kinematics s = ½at² | C09 (D17) | C09 (gloss in D17 step 3's *why*) |
| small-angle approximation | C09 (D17) | knowledge/primers.md: small-angle approximation (ch03 P100) — reminder |
| Coriolis deflection Ωut² | C09 | C09 (N61 note, D17) |
| highs and lows (flow round pressure centres) | C09 | C09 (N62 note) |
| angular-acceleration term | C09 | C09 (N63 note) |
| centrifugal term, its potential and effective gravity; geopotential | C09 | C09 (N64 note, D18) |
| cylindrical gradient of −½Ω²R² | C09 (D18) | C09 (gloss in D18 step 4's *why*) |
| Appendix-B curvilinear operators | C09 (N65) | C09 (gloss in N65) |
| plt.quiver | C09 | knowledge/primers.md: plt.contour, plt.quiver and plt.streamplot (ch02 P78) — reminder |
| σ:R = 0 for symmetric σ | C10 | R07 |
| Gibbs relation along a particle (4.61) | C10 | R08 |
| first law of thermodynamics | C10 | C10 (N66 note; Ch. 1 §1.8) |
| internal and kinetic energy per unit mass | C10 | knowledge/primers.md: internal and kinetic energy per unit mass (ch01 P32) — reminder |
| power of a force and heat flux through a surface | C10 | primer (in C10) |
| energy budget of a material volume (4.46)–(4.48) | C10 | C10 (N66–N68 notes) |
| total-energy equation (4.53) and the (4.51) slip | C10 | C10 (N69–N73 notes, D19) |
| stress work = deformation work + force work (4.54) | C10 | C10 (N74 note) |
| total energy following a particle (4.55) | C10 | C10 (N75 note, D20) |
| chain rule for the kinetic energy | C10 (D21) | primer (in C10) |
| chain rule | C10 | knowledge/primers.md: chain rule (ch01 P49) — reminder |
| mechanical-energy equation (4.56) | C10 | C10 (N76 note, D21) |
| quotient rule D(1/ρ)/Dt = −ρ⁻²Dρ/Dt | C10 (D22) | C10 (gloss in D22 step 4's *why*) |
| specific volume v = 1/ρ | C10 | C10 (D22 step 4; Ch. 1 §1.8) |
| internal-energy equation (4.57) | C10 | C10 |
| completing the square for tensors | C10 (D23) | primer (in C10) |
| dissipation rate ε (4.58) and ε ≥ 0 | C10 | C10 (N77 note, D23) |
| Newtonian viscous stress (4.59) | C10 | C10 (N78 note) |
| Fourier's law q = −k∇T | C10 | C10 (N79 note; Ch. 1 §1.5) |
| energy equation with conduction (4.60); ledger 7 = 7 | C10 | C10 (N79 note) |
| entropy equation (4.62) | C10 | C10 (N80 note) |
| entropy production (4.63); second law ⇒ μ, μ_v, k ≥ 0 | C10 | C10 (N81 note) |
| Couette flow with viscous heating | C10 | C10 (worked example) |
| sine-series transient of the heated gap (E7) | C10 (E7) | C10 (gloss in the E7 cell; the method of Ch. 1's Couette start-up, `core.diffusion.couette_startup_profile`) |
| np.random.default_rng | C10 (check) | knowledge/primers.md: np.random.default_rng (ch01 P10) — reminder |
| Bernoulli equations are consequences, not laws | C11 | C11 (N85 note) |
| Euler with a potential (4.66) | C11 | C11 (N86 note) |
| barotropic pressure function ∫dp/ρ (4.67); dummy variable p′ | C11 | C11 (N87 note) |
| fundamental theorem of calculus with a variable limit | C11 (D24) | knowledge/primers.md: fundamental theorem of calculus (ch02 P84) — reminder |
| permutations and the alternating tensor ε_ijk | C11 (D24) | knowledge/primers.md: permutations, cyclic order and parity (ch02 P72) — reminder |
| Lamb identity (4.68) | C11 | C11 (N88 note, D24) |
| Bernoulli function B, (4.69) | C11 | C11 (N89 note) |
| triple product with a repeated vector a·(a × b) = 0 | C11 (D25) | C11 (gloss in D25 step 2's *why*) |
| zero gradient everywhere ⇒ constant in space | C11 (D25) | C11 (gloss in D25 step 5's *why*) |
| Lamb surfaces (4.70) | C11 | C11 (N90 note) |
| Bernoulli 1 (4.71) | C11 | C11 |
| irrotational ⇒ B constant everywhere (4.72) | C11 | C11 (N91 note) |
| persistence of irrotationality (named) | C11 | C11 (N92 note) |
| energy Bernoulli (4.76)–(4.78); enthalpy; stagnation temperature | C11 | C11 (N94–N96 notes) |
| Rankine vortex and its pressure | C11 | C11 (worked example; the vortex of Ch. 3 §3.5) |
| polar coordinates as a moving basis | C11 | knowledge/primers.md: polar coordinates as a moving basis (ch03 P105) — reminder |
| decision table of Bernoulli forms | C11 | C11 (N102 note) |
| velocity potential u = ∇φ (4.73) | C12 | R09 |
| unsteady Bernoulli (4.74), (4.75); gauge sign slip | C12 | C12 (N93 note, D26) |
| sympy Integral in the gauge check | C12 | C12 (gloss in the check cell's comment) |
| Lamb form with viscosity (4.79); irrotational, viscous (4.80)–(4.83) | C12 | C12 (N97–N101 notes) |
| directional derivative along a streamline | C12 (N99) | knowledge/primers.md: level sets and the directional derivative (ch02 P75) — reminder |
| U-tube oscillation ω² = 2g/L | C12 | C12 (worked example) |
| added mass of an accelerating sphere | C12 | C12 (worked example, pointer to Ch. 6) |
| hydrostatic reference state | C13 | R10 |
| perturbation pressure and density (4.84); primes as perturbations | C13 | C13 (N106 note) |
| gravity absorbed into p′ for constant density (4.85) | C13 | C13 (N107 note) |
| order-of-magnitude scaling | C13 | primer (in C13) |
| reduced gravity and buoyancy | C13 | primer (in C13) |
| thermal expansion coefficient α | C13 | C13 (worked example; Ch. 1 §1.9) |
| Boussinesq validity conditions | C13 | C13 (N108 note) |
| Boussinesq momentum (4.86) | C13 | C13 |
| energy equation in vector form (4.87); pressure work ⇒ C_p | C13 | C13 (N109, N110 notes, D28) |
| perfect gas p = ρRT, R = C_p − C_v, α = 1/T | C13 (D28) | knowledge/primers.md: perfect-gas law as a working model (ch01 P33) — reminder |
| heat equation with C_p (4.88) | C13 | C13 (N111 note) |
| viscous heating negligible (Eckert) | C13 | C13 (N112 note) |
| Boussinesq heat equation (4.89); κ = k/ρC_p | C13 | C13 (N113 note) |
| linear equation of state; buoyancy b and N² with both lapse-rate conventions | C13 | C13 (N114 note; N² from Ch. 1 §1.10) |
| Gaussian blob solution of the heat equation | C13 | C13 (N113 note and animation) |
| blob rise parcel model (E8, our extension) | C13 (E8) | C13 (gloss in the E8 explainer cell) |
| surface tension | C14 | R11 |
| boundary conditions | C14 | knowledge/primers.md: boundary conditions (ch01 P20) — reminder |
| what must be specified at boundaries | C14 | C14 (N115 note) |
| pillbox interface conditions | C14 | C14 (N116 note) |
| no-slip and no temperature jump; slip length; Knudsen number | C14 | C14 (N117 note) |
| moving level set and its normal speed | C14 | primer (in C14) |
| unit normal ∇η/‖∇η‖ | C14 (D29) | knowledge/primers.md: level sets and the directional derivative (ch02 P75) — reminder |
| surface moving with u_s (4.90) | C14 | C14 (N118 note) |
| kinematic boundary condition (4.91) | C14 | C14 |
| relative normal velocity (4.92) and mass flux (4.93) | C14 | C14 (N119, N120 notes) |
| Helmholtz free energy (4.94), (4.95); Legendre device | C14 | C14 (N121, N122 notes) |
| Gibbs free energy (as the model for f) | C14 | knowledge/primers.md: Gibbs free energy (ch01 P50) — reminder |
| free energy of an interface (4.96) | C14 | C14 (N123 note) |
| Marangoni flow (named) | C14 | C14 (N124 note) |
| curved cap forces (4.97), (4.98); curve C slip | C14 | C14 (N125, N126 notes) |
| Laplace's pressure jump (1.5) derived | C14 | C14 (N127 note) |
| radius of curvature and principal radii | C14 | C14 (N125 note; Ch. 1 §1.6) |
| capillary length | C14 | C14 (N128 note) |
| meniscus at a wall; contact angle; Ex. 4.7 slips | C14 | C14 (N129 note) |
| linear water wave as a test field | C14 | C14 (worked example and animation text; derived in Ch. 7) |
| sphere-drag Π groups (4.99) | C15 | R12 |
| two routes to dimensionless groups; dynamic similarity | C15 | C15 (N130 note) |
| flow parameters l, U, Ω, ρ, μ | C15 | C15 (N132 note) |
| scaled variables and the chain rule | C15 (D30) | primer (in C15) |
| dimensionless variables (4.100) | C15 | C15 (N133 note) |
| dimensionless Navier–Stokes (4.101) | C15 | C15 |
| Strouhal, Reynolds, Froude numbers (4.102)–(4.104) | C15 | C15 (N134–N136 notes) |
| internal Froude and Richardson numbers (4.105) | C15 | C15 (N137 note) |
| lapse-rate conventions (Kundu dT/dz vs meteorology −dT/dz) | C13, C15 | C13 (N114 note), C15 (N137 note) |
| similarity of all outputs | C15 | C15 (N138 note) |
| pressure coefficient (4.106) | C15 | C15 (N139 note) |
| steady boundary conditions, time scale l/U | C15 | C15 (N140 note) |
| purely oscillating body | C15 | C15 (N141 note) |
| drag and lift coefficients (4.107), (4.108); reference areas | C15 | C15 (N142, N143 notes) |
| ship drag C_D(Fr, Re) | C15 | C15 (N144 note) |
| compressible scalings (4.109), compressibility ∝ M² (4.110), Mach number (4.111) | C15 | C15 (N145–N147 notes) |
| enthalpy form of the energy equation (4.112) | C15 | C15 (N148 note) |
| thermal scalings (4.113), dimensionless energy equation (4.114) | C15 | C15 (N149, N150 notes) |
| Eckert and Prandtl numbers (4.115), (4.116) | C15 | C15 (N151, N152 notes) |
| Weber, Bond, capillary numbers (4.117)–(4.119) | C15 | C15 (N153–N155 notes) |
| ship-model Froude scaling and friction correction; ITTC friction line | C15 | C15 (N156 note, worked example; the friction line glossed there) |
| sphere drag correlation (Morrison 2013) | C15 | C15 (N131 note) |
| Rossby number U/(2Ωl) (forward pointer) | C15 (E9) | C15 (the "What would change if…" cell and E9's Explain §5) |
| show_viz | C03 | knowledge/primers.md: show_viz (ch01 P18) — reminder |
| matplotlib figures | §4.2 | knowledge/primers.md: matplotlib figures (ch01 P01) — reminder |
| numpy arrays | §4.1 | knowledge/primers.md: numpy arrays (ch01 P03) — reminder |
| tuple unpacking | C02 | knowledge/primers.md: tuple unpacking (ch01 P14) — reminder |
| Python dictionaries | C04 | knowledge/primers.md: Python dictionaries (ch01 P23) — reminder |
| np.where | C13 (verdicts) | knowledge/primers.md: np.where and np.select (ch01 P46) — reminder |
| exponent rules | C07 (D08) | knowledge/primers.md: exponent rules (ch01 P43) — reminder |
| mean-value theorem for integrals | C02 (D02) | knowledge/primers.md: mean-value theorem for integrals (ch02 P85) — reminder |
| volume and surface integrals as midpoint sums | C01 (3-D check) | knowledge/primers.md: volume and surface integrals as midpoint sums (ch02 P83) — reminder |
| eigenvalues of S (principal directions) | C07 | knowledge/primers.md: eigenvalues and eigenvectors (ch02 P80) — reminder |
| separation of variables (rocket, tank, meniscus) | C04 | knowledge/primers.md: separation of variables (ch01 P42) — reminder |

---

## Part F — derivation storyboards

Builders copy these word for word into `nb.derivation(key, title, goal=…, start=(tex, plain), plan=[…], uses=[…],
steps=[dict(did, tex, why, plain)], result=(tex, plain), interpret=…, check=…, check_src=…)` and into the explainer's
`derivations: [...]` (phones may shorten *why* to its first sentence; `live` and `set` are the explainer's). Every step is
one move; *why* names the rule and says why we make the move (≤ 35 words, ≥ 6); *did* ≤ 8 words. The book's own moves
were read on the rendered pages (p124–p125 for D01–D03, p126–p127 for D04, p128–p129 for D05, p132–p133 for D06, p137–p138
for D07, p138 for D08 (the book only states the result), p139–p141 for D09–D10, p141–p142 for D11–D13, p143–p145 for
D14–D17, p146 for D18, p149–p151 for D19–D23, p155 and p157 for D24–D26, p162–p164 for D27–D28, p166 for D29, p173 for
D30); the gaps listed in `analysis/ch04.md` §2b are filled and the notebook says so ("the book skips this move"). Every
equation named by number is written out. Colours: inflow blue, outflow orange, storage purple, force rose; local blue,
advective teal, pressure orange, viscous rose, gravity grey, buoyancy blue; Coriolis amber, centrifugal purple; kinetic
energy teal, heat rose. (No line in this part starts with a table bar; absolute values are written with \lvert \rvert or
in words.)

### D01 · Mass for an arbitrarily moving control volume: (4.1) → (4.2) → (4.3) + (4.4) → (4.5) — ★★, 9 steps, in C01 (notebook · `control_volume_budgets`)
- **Goal.** Turn "a sealed balloon keeps its mass" into a mass budget for any box we choose — fixed, sliding, stretching —
  written only with things we can measure on the box: the mass inside and what crosses its walls.
- **Start.** (4.1): $\frac{d}{dt}\int_{V(t)}\rho\,dV=0$ — *in words:* a material volume (always the same particles) keeps its
  mass.
- **Plan.** (1) Expand the balloon's rate with the transport theorem (b = u). (2) Write the same theorem for our box V*.
  (3) Let the box and the balloon fill the same region at one instant and compare term by term. (4) Substitute and merge
  the surface integrals; check two special boxes.
- **Tools.** Reynolds transport theorem $\frac{d}{dt}\int_{V^*}F\,dV=\int_{V^*}\frac{\partial F}{\partial t}dV+\int_{A^*}F\,\mathbf b\cdot\mathbf n\,dA$
  (3.35) and the material volume (R01) · signed swept volume of a moving wall (ch03 P110) · the coincident control
  volume (N06) · linearity of integrals.
- **Assumptions.** ρ and u continuous (steps 2–3); A* piecewise smooth (step 3); V* and V coincide at the instant t (steps
  4–6); u and b measured in the same frame (step 8).
- **Steps.**
  1. *did:* Start from the sealed balloon · *tex:* $\frac{d}{dt}\int_{V(t)}\rho\,dV=0$ · *why:* (4.1): a material volume holds
     the same particles for ever and mass is neither made nor destroyed; we want this rate in terms we can evaluate. ·
     *plain:* The balloon's mass never changes. · *set:* {mode: 'balloon', box: 'material', t: 0}.
  2. *did:* Expand with the transport theorem, b = u · *tex:* $\int_{V(t)}\frac{\partial\rho}{\partial t}dV+\int_{A(t)}\rho\,\mathbf u\cdot\mathbf n\,dA=0$ ·
     *why:* RTT (3.35) with F = ρ; the skin of a material volume moves with the fluid, so its surface velocity is b = u.
     This is the book's (4.2). · *plain:* Density changes inside plus mass carried by the moving skin add to zero. ·
     *live:* "for the particles in [1, 2] m at t = 0: −1 + (2 − 1) = 0 kg/(m² s)".
  3. *did:* Write the theorem for our own box · *tex:* $\frac{d}{dt}\int_{V^*}\rho\,dV-\int_{V^*}\frac{\partial\rho}{\partial t}dV-\int_{A^*}\rho\,\mathbf b\cdot\mathbf n\,dA=0$ ·
     *why:* RTT (3.35) holds for any volume whose surface moves at b; with F = ρ, moved to one side, it is the book's (4.3). V* is the box we can draw. · *plain:* The box's mass changes because density changes
     inside it and because its walls sweep fluid in or out. · *set:* {box: 'fixed'}.
  4. *did:* Let the box and balloon coincide now · *tex:* $V^*(t)=V(t),\quad A^*(t)=A(t)$ · *why:* We may choose which
     material volume to follow: the particles that fill the box at this instant. A moment later the two part company,
     because the balloon moves with u and the box with b. · *plain:* At one instant the balloon exactly fills the box. ·
     *set:* {t: 0} (both outlines drawn) · *watch:* "the outlines coincide now and separate as soon as time runs".
  5. *did:* Equate the two local integrals · *tex:* $\int_{V^*}\frac{\partial\rho}{\partial t}dV=\int_{V}\frac{\partial\rho}{\partial t}dV$ ·
     *why:* The book skips why: same integrand ∂ρ/∂t, same region, same instant give the same number. The rates
     d/dt∫ρ dV are NOT equal: they depend on how each volume moves next. · *plain:* Adding the same density rates over the
     same region gives the same total.
  6. *did:* Replace it by the flux through the skin · *tex:* $\int_{V^*}\frac{\partial\rho}{\partial t}dV=-\int_{A}\rho\,\mathbf u\cdot\mathbf n\,dA=-\int_{A^*}\rho\,\mathbf u\cdot\mathbf n\,dA$ ·
     *why:* Step 2: the local integral over V is minus the flux through A; A and A* are the same surface now, so those fluxes agree. This chain is the book's (4.4). · *plain:* The density changes inside add
     up to what flows out through the skin.
  7. *did:* Substitute into the box's theorem · *tex:* $\frac{d}{dt}\int_{V^*}\rho\,dV+\int_{A^*}\rho\,\mathbf u\cdot\mathbf n\,dA-\int_{A^*}\rho\,\mathbf b\cdot\mathbf n\,dA=0$ ·
     *why:* In step 3 the middle term is minus the local integral; by step 6 that is plus the flux of ρu. Only box
     quantities remain. · *plain:* Box mass rate + fluid carried out − fluid swept out by the walls = 0. · *live:* "fixed
     box [1, 2] m: −1 + (2 − 1) − 0 = 0".
  8. *did:* Merge the two surface integrals · *tex:* $\frac{d}{dt}\int_{V^*}\rho\,dV+\int_{A^*}\rho\,(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=0$ ·
     *why:* Same surface A*, same ρ and n: the integrals combine (linearity). u − b is the velocity relative to the wall, so u and b share one frame. This is (4.5). · *plain:* The mass inside changes
     only by what crosses the walls relative to them.
  9. *did:* Check the two special boxes · *tex:* $\mathbf b=\mathbf u:\ \frac{d}{dt}\int_{V}\rho\,dV=0;\qquad\mathbf b=0:\ \frac{d}{dt}\int_{V^*}\rho\,dV=-\oint_{A^*}\rho\,\mathbf u\cdot\mathbf n\,dA$ ·
     *why:* A rule for any box must reproduce the known cases: a box moving with the fluid has no relative flux (that is
     (4.1)); a fixed box loses exactly what flows out. · *plain:* The general rule contains both the balloon and the fixed
     box. · *set:* {box: 'material'} then {box: 'fixed'}.
- **Result.** $\frac{d}{dt}\int_{V^*(t)}\rho\,dV+\int_{A^*(t)}\rho\,(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=0$ (4.5) — *in words:* the mass
  in any box changes only by the mass crossing its walls, counted relative to the walls.
- **Check.** Units: kg/s both terms ✓. Special cases b = u and b = 0 (step 9) ✓. Numbers (expanding flow, a = 1 s⁻¹,
  ρ₀ = 1 kg/m³, box [1, 2] m at t = 0): fixed box −1 + (2 − 1) = 0; box sliding at b = 1 m/s: −1 + [(2 − 1) − (1 − 1)] = 0;
  material box 0 + 0 = 0 ✓ (the notebook's code cell).
- **What it means.** One template for every budget of the book: replace ρ by ρu and add forces (C04), by ρ(e + ½u²) and add
  work and heat (C10). It needs a continuum and no mass sources; the box may move or deform however we like.
- **Traps.** Thinking the two d/dt∫ terms are equal because the volumes coincide (only the ∂/∂t integrals are). Mixing
  frames for u and b. Using an inward normal or |b·n| (the sign of b·n says whether a wall sweeps fluid in or out).

### D02 · The continuity equation: (4.2) → (4.6) → (4.7) — ★★, 8 steps, in C02 (notebook)
- **Goal.** Say what "mass is conserved" means at a single point: a differential equation linking how fast the density
  changes to how the mass flux spreads out.
- **Start.** (4.2): $\int_{V(t)}\frac{\partial\rho}{\partial t}dV+\int_{A(t)}\rho\,\mathbf u\cdot\mathbf n\,dA=0$ — *in words:* for any material
  volume, density changes inside balance the mass carried through its skin.
- **Plan.** (1) Turn the surface integral into a volume integral with Gauss. (2) Put both under one integral. (3) Use
  that this holds for every volume, and prove (the book gives one sentence) that the integrand must then vanish.
- **Tools.** Gauss' theorem $\int_V\nabla\cdot\mathbf Q\,dV=\int_A\mathbf Q\cdot\mathbf n\,dA$ (2.30) (R02) · continuity of a function and
  the small-ball argument (primer in C02) · the localisation lemma (N07) · summation convention.
- **Assumptions.** ρu continuously differentiable (step 3); the statement holds for every material volume (step 5); the
  integrand is continuous (step 6).
- **Steps.**
  1. *did:* Start from the expanded material law · *tex:* $\int_{V}\frac{\partial\rho}{\partial t}dV+\int_{A}\rho\,\mathbf u\cdot\mathbf n\,dA=0$ ·
     *why:* This is (4.2), D01 step 2: the transport theorem with b = u applied to (4.1). We want every term as a volume
     integral. · *plain:* The balloon budget, written out.
  2. *did:* Name the mass-flux vector · *tex:* $\mathbf Q=\rho\,\mathbf u$ · *why:* The surface term is the flux of Q through A;
     Q is mass per unit area per second. Gauss' theorem applies to any continuously differentiable vector field. ·
     *plain:* ρu is how much mass crosses a unit area each second.
  3. *did:* Apply Gauss to the surface term · *tex:* $\int_A\rho\,\mathbf u\cdot\mathbf n\,dA=\int_V\nabla\cdot(\rho\mathbf u)\,dV$ · *why:*
     Divergence theorem (2.30), R02: net outflow through a closed surface equals the divergence summed inside. It needs ρu
     smooth, true away from shocks. · *plain:* What leaves through the skin is what the inside sources.
  4. *did:* Collect under one integral · *tex:* $\int_{V}\Big\{\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)\Big\}dV=0$ ·
     *why:* Both integrals are over the same V, so they add (linearity). This is the book's (4.6). · *plain:* The sum of
     local change and flux divergence integrates to zero.
  5. *did:* Note it holds for every volume · *tex:* $\int_Vf\,dV=0\ \ \text{for every }V,\qquad f\equiv\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)$ ·
     *why:* The book skips this: V was any material volume, and any region is filled by some particles at time t. So the
     integral vanishes over every region — that is what licenses the next moves. · *plain:* Wherever we draw the region,
     the integral is zero.
  6. *did:* Suppose f is positive at a point · *tex:* $f(\mathbf x_0)=c>0\ \Rightarrow\ f>\tfrac c2\ \text{on a small ball }B_r(\mathbf x_0)$ ·
     *why:* Continuity of f (primer): near x₀ its values stay within c/2 of c. We argue by contradiction; a negative value
     is handled the same way with −f. · *plain:* If f were positive somewhere, it would be positive on a whole little
     ball.
  7. *did:* Integrate over that ball · *tex:* $\int_{B_r}f\,dV>\tfrac c2\cdot\tfrac43\pi r^3>0$ · *why:* An integrand larger than
     c/2 everywhere gives an integral larger than c/2 times the volume; but step 5 says every region, this ball included,
     gives zero: contradiction. · *plain:* The ball's integral cannot be both positive and zero.
  8. *did:* Conclude the integrand vanishes · *tex:* $\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0\quad\text{or}\quad\frac{\partial\rho}{\partial t}+\frac{\partial}{\partial x_i}(\rho u_i)=0$ ·
     *why:* No point can have f ≠ 0, so f = 0 everywhere — the localisation lemma (N07). The index form is the summation
     convention. This is (4.7). · *plain:* At every point, the density falls exactly as fast as mass flows away.
- **Result.** $\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0$ (4.7) — *in words:* the continuity equation, mass conservation
  at a point.
- **Check.** Units: kg/(m³ s) for both terms ✓. Expanding flow at x = 1.5 m, t = 0: ∂ρ/∂t = −1, ∂(ρu)/∂x = +1, sum 0 ✓.
  Constant ρ: ∇·u = 0 ✓ (D03). Integrated over a closed box it gives back the fixed-box budget of D01 step 9 ✓.
- **What it means.** One scalar equation for four unknown fields (ρ and three u's) — not enough alone; it pairs with
  momentum. Where the flux converges (∇·(ρu) < 0) density rises. It fails only where ρu is not differentiable (a shock),
  where the integral form (4.5) still holds.
- **Traps.** Applying the lemma to a single volume (it needs *every* volume). Forgetting continuity of f (a function
  nonzero at isolated points can integrate to zero). Reading ∂ρ/∂t as the particle's density rate (that is Dρ/Dt, D03).

### D03 · Continuity with D/Dt; incompressible ⇒ ∇·u = 0: (4.7) → (4.8) → (4.9) → (4.10) — ★, 5 steps, in C02 (notebook)
- **Goal.** Rewrite continuity in terms of what a fluid particle experiences, and see what "incompressible" requires of
  the velocity.
- **Start.** (4.7): $\frac{\partial\rho}{\partial t}+\frac{\partial}{\partial x_i}(\rho u_i)=0$ — *in words:* the continuity equation.
- **Plan.** (1) Expand the flux divergence with the product rule. (2) Spot the material derivative. (3) Divide by ρ.
  (4) Impose Dρ/Dt = 0.
- **Tools.** Product rule for a divergence (primer in C02) · material derivative $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+u_i\frac{\partial F}{\partial x_i}$
  (3.5) (Ch. 3 §3.2).
- **Assumptions.** ρ > 0 (step 3); each particle keeps its density (step 4).
- **Steps.**
  1. *did:* Expand the flux divergence · *tex:* $\frac{\partial\rho}{\partial t}+u_i\frac{\partial\rho}{\partial x_i}+\rho\frac{\partial u_i}{\partial x_i}=0$ ·
     *why:* Product rule for a divergence (primer, the book's (B3.6)): ∂(ρu_i)/∂x_i = u_i∂ρ/∂x_i + ρ∂u_i/∂x_i. We want
     to separate "moving through a density gradient" from "spreading out". · *plain:* Density changes by being carried and
     by the flow spreading.
  2. *did:* Recognise the material derivative · *tex:* $\frac{D\rho}{Dt}+\rho\,\nabla\cdot\mathbf u=0$ · *why:* The first two terms are
     (3.5), $\frac{D\rho}{Dt}=\frac{\partial\rho}{\partial t}+u_i\frac{\partial\rho}{\partial x_i}$: the rate a particle sees. · *plain:* A particle's
     density rate plus density times its volume growth rate is zero.
  3. *did:* Divide by ρ · *tex:* $\frac1\rho\frac{D\rho}{Dt}+\nabla\cdot\mathbf u=0$ · *why:* ρ > 0 in any real fluid, so division
     is allowed; it turns the density rate into a fractional rate. This is (4.8). · *plain:* The fractional density rate is
     minus the volume growth rate (Ch. 3's (3.14)). · *live:* "expanding flow: −1 + 1 = 0 s⁻¹".
  4. *did:* Impose incompressibility · *tex:* $\frac{D\rho}{Dt}\equiv\frac{\partial\rho}{\partial t}+\mathbf u\cdot\nabla\rho=0$ · *why:* The
     definition (4.9): every particle keeps its own density. This is weaker than ρ = const — different particles may differ
     (a stratified lake). · *plain:* No particle gets denser or lighter.
  5. *did:* Substitute into (4.8) · *tex:* $\nabla\cdot\mathbf u=0$ · *why:* The first term of (4.8) vanishes by (4.9), leaving
     (4.10): no particle changes volume. · *plain:* An incompressible flow has zero divergence.
- **Result.** $\frac1\rho\frac{D\rho}{Dt}+\nabla\cdot\mathbf u=0$ (4.8), and for incompressible flow $\nabla\cdot\mathbf u=0$ (4.10) — *in words:*
  squeezing a particle raises its density; if densities cannot change, volumes cannot either.
- **Check.** Units 1/s ✓. Stratified shear flow u = (U(z), 0, 0), ρ(z): Dρ/Dt = 0 and ∇·u = 0 though ρ varies ✓.
  Expanding flow: (1/ρ)Dρ/Dt = −1 s⁻¹ = −∇·u ✓.
- **What it means.** ∇·u = 0 is the constraint of Ch. 5–13's incompressible flows; C15 shows that for a gas it holds
  to O(M²). It fails for sound, shocks and fast gas flows.
- **Traps.** Equating incompressible with ρ = const. Treating §4.11's re-display $\nabla\cdot\mathbf u=-\frac{1}{\rho c^2}\frac{Dp}{Dt}$ as the
  incompressible condition (it is a general identity). Dividing by ρ where ρ could vanish (never in a fluid).

### D04 · The 2-D stream function: ρu = ∂ψ/∂y, ρv = −∂ψ/∂x, and ψ₂ − ψ₁ = flux — ★★, 8 steps, in C03 (notebook · `stream_function_spacing`)
- **Goal.** From the two-function form of the mass flux, get the familiar 2-D stream function and prove that the flux
  between two streamlines is the difference of their ψ values — the book leaves the last part to Exercise 4.8.
- **Start.** (4.12): $\rho\mathbf u=\nabla\times\boldsymbol\Psi$ with $\boldsymbol\Psi=\chi\nabla\psi$ — *in words:* a mass flux written as a curl
  satisfies steady continuity (4.11), $\nabla\cdot(\rho\mathbf u)=0$, automatically.
- **Plan.** (1) Expand the curl of χ∇ψ. (2) Choose χ = −z for a plane flow and take the cross product. (3) Integrate the
  flux across a curve and recognise an exact derivative.
- **Tools.** ∇·∇× = 0 and ∇×∇ = 0 (gloss, Ch. 2 §2.13) · curl of a product (gloss with the sympy check in N14's cell) ·
  right-hand rule (ch02 P74) · chain rule along a path (ch03 P91) · fundamental theorem of calculus (ch02 P84).
- **Assumptions.** Steady flow (continuity is ∇·(ρu) = 0); the flow is plane (no z dependence, w = 0) from step 4; ψ smooth.
- **Steps.**
  1. *did:* Start from the vector potential · *tex:* $\rho\mathbf u=\nabla\times(\chi\nabla\psi)$ · *why:* (4.12): the divergence of a
     curl is zero (Ch. 2), so this ρu satisfies (4.11) for any χ, ψ. We choose Ψ = χ∇ψ because it leads to stream surfaces.
     · *plain:* A flux built as a curl can never violate steady continuity.
  2. *did:* Expand the curl of a product · *tex:* $\nabla\times(\chi\nabla\psi)=\nabla\chi\times\nabla\psi+\chi\,\nabla\times\nabla\psi$ · *why:*
     Product rule for the curl (gloss): in index form ε_ijk∂_j(χ∂_kψ) = ε_ijk∂_jχ∂_kψ + χε_ijk∂_j∂_kψ. · *plain:* The curl
     splits into two pieces.
  3. *did:* Drop the curl of a gradient · *tex:* $\rho\mathbf u=\nabla\chi\times\nabla\psi$ · *why:* ∇×∇ψ = 0 for a smooth ψ (mixed
     derivatives are equal and ε is antisymmetric, Ch. 2). · *plain:* The mass flux is perpendicular to both gradients: flow
     runs along the surfaces χ = const and ψ = const.
  4. *did:* Choose χ = −z for a plane flow · *tex:* $\chi=-z\ \Rightarrow\ \nabla\chi=-\mathbf e_z$ · *why:* In a plane flow every
     streamline lies in a plane z = const, so z (or −z) is a stream function; the minus sign is the usual convention. ·
     *plain:* One family of stream surfaces is simply the stack of horizontal planes. · *set:* {flow: 'uniform'}.
  5. *did:* Write ψ's gradient in the plane · *tex:* $\nabla\psi=\Big(\frac{\partial\psi}{\partial x},\frac{\partial\psi}{\partial y},0\Big)$ ·
     *why:* ψ = ψ(x, y) does not depend on z in a plane flow. · *plain:* ψ varies only across the plane. · *live:* "at the
     probe: ∂ψ/∂x = …, ∂ψ/∂y = …".
  6. *did:* Take the cross product · *tex:* $\rho\mathbf u=-\mathbf e_z\times\nabla\psi=\Big(\frac{\partial\psi}{\partial y},\,-\frac{\partial\psi}{\partial x},\,0\Big)$ ·
     *why:* e_z × e_x = e_y and e_z × e_y = −e_x (right-hand rule), so −e_z × (ψ_x e_x + ψ_y e_y) = ψ_y e_x − ψ_x e_y. ·
     *plain:* ρu = ∂ψ/∂y and ρv = −∂ψ/∂x.
  7. *did:* Integrate the flux across a curve · *tex:* $\int_1^2\rho\mathbf u\cdot\mathbf n\,ds=\int_1^2\Big(\frac{\partial\psi}{\partial y}\frac{dy}{ds}+\frac{\partial\psi}{\partial x}\frac{dx}{ds}\Big)ds$ ·
     *why:* The book skips this (Exercise 4.8): along a curve from 1 to 2 with tangent (dx/ds, dy/ds) the normal to the
     right of travel is n = (dy/ds, −dx/ds); insert step 6. · *plain:* The flux through a gate, written with ψ. · *set:*
     {flow: 'stagnation', q: 1} with the gate from ψ = 1 to ψ = 2.
  8. *did:* Recognise an exact derivative · *tex:* $\int_1^2\rho\mathbf u\cdot\mathbf n\,ds=\int_1^2\frac{d\psi}{ds}\,ds=\psi_2-\psi_1$ · *why:*
     Chain rule along the curve (ch03 P91): dψ/ds = ψ_x dx/ds + ψ_y dy/ds; then the fundamental theorem of calculus. Any
     path between the same end points gives the same flux. · *plain:* The flux between two streamlines is the difference of
     their labels. · *live:* "ψ₂ − ψ₁ = 2 − 1 = 1.000 m²/s" · *watch:* "drag a handle: the flux readout does not change".
- **Result.** $\rho u=\partial\psi/\partial y,\ \rho v=-\partial\psi/\partial x$ and $\int_1^2\rho\mathbf u\cdot\mathbf n\,ds=\psi_2-\psi_1$ — *in words:* one
  scalar carries a 2-D flow; its contours are streamlines, and the flux between two of them is the difference of the labels.
- **Check.** Units: ψ in kg/(m s) (mass flux per metre of depth; m²/s when ρ = 1) ✓. Uniform stream ψ = Uy: u = U, v = 0,
  flux between y₁ and y₂ = U(y₂ − y₁) ✓. Ch. 3's cylinder ψ returns Ch. 3's velocity (C03 code) ✓. Continuity:
  ∂(ρu)/∂x + ∂(ρv)/∂y = ψ_yx − ψ_xy = 0 ✓.
- **What it means.** Contours of ψ at equal steps show direction *and* speed (crowded = fast). Used for every 2-D flow
  in Ch. 6, 8, 9 and the quasi-geostrophic flows of Ch. 13. Needs a plane (or axisymmetric) flow; in general 3-D flow two
  stream functions are needed.
- **Traps.** The opposite sign convention of many GFD texts (u = −∂ψ/∂y). Orienting the normal to the left of travel (the
  flux changes sign). Forgetting that the flux is per unit depth.

### D05 · Momentum for an arbitrarily moving control volume: (4.13) → (4.14) → (4.15) → (4.16a–d) → (4.17) — ★★, 8 steps, in C04 (notebook · `control_volume_budgets`)
- **Goal.** Give Newton's second law the same "any box" form that D01 gave mass conservation, so forces can be read from
  fluxes through a box we choose.
- **Start.** (4.13): $\frac{d}{dt}\int_{V(t)}\rho\mathbf u\,dV=\int_{V(t)}\rho\mathbf g\,dV+\int_{A(t)}\mathbf f\,dA$ — *in words:* the momentum of a
  material volume changes at the rate of the forces on it.
- **Plan.** (1) Expand the material rate (b = u) component by component. (2) Write the theorem for our box. (3) Use the
  coincidence trick of D01 four times. (4) Substitute and merge the flux terms; check the special boxes.
- **Tools.** RTT (3.35) applied to each component (gloss) · D01's coincidence move · momentum flux through a surface
  (primer in C04) · Newton's second law (ch01 P09).
- **Assumptions.** As D01; f is the traction (force per area) exerted on the fluid by whatever is outside, pressure
  included.
- **Steps.**
  1. *did:* Start from Newton's law for the balloon · *tex:* $\frac{d}{dt}\int_{V(t)}\rho\mathbf u\,dV=\int_{V(t)}\rho\mathbf g\,dV+\int_{A(t)}\mathbf f\,dA$ ·
     *why:* (4.13): momentum per volume ρu, body force per mass g, surface force per area f, all for the same particles.
     We want it for a box. · *plain:* Momentum changes only by forces. · *set:* {mode: 'wake', H: 1.0, budget: 'momentum'}.
  2. *did:* Expand with the transport theorem, b = u · *tex:* $\int_{V}\frac{\partial(\rho\mathbf u)}{\partial t}dV+\int_{A}\rho\mathbf u\,(\mathbf u\cdot\mathbf n)\,dA=\int_{V}\rho\mathbf g\,dV+\int_{A}\mathbf f\,dA$ ·
     *why:* RTT (3.35) applied to each component F = ρu_j (a vector is three scalars); the skin moves with u. This is
     (4.14). · *plain:* Local change of momentum plus momentum carried by the moving skin equals the forces.
  3. *did:* Write the theorem for our box · *tex:* $\int_{V^*}\frac{\partial(\rho\mathbf u)}{\partial t}dV=\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV-\int_{A^*}\rho\mathbf u\,(\mathbf b\cdot\mathbf n)\,dA$ ·
     *why:* RTT (3.35) with F = ρu for V*, rearranged. This is (4.15); the book prints an extra "= 0" at its end, a slip:
     this line is an identity between three terms. · *plain:* The local momentum change in the box = its total rate minus
     what the walls sweep.
  4. *did:* Let the box and balloon coincide now · *tex:* $\int_{V}\{\cdot\}\,dV=\int_{V^*}\{\cdot\}\,dV,\qquad\int_{A}\{\cdot\}\,dA=\int_{A^*}\{\cdot\}\,dA$ ·
     *why:* Same integrands, same region, same instant (D01 steps 4–5): the four equalities (4.16a–d) for ∂(ρu)/∂t,
     ρu(u·n), ρg and f. · *plain:* Every integral over the balloon equals the one over the box, right now. · *watch:* "the
     four bar groups light one after another".
  5. *did:* Move step 2 onto the box · *tex:* $\int_{V^*}\frac{\partial(\rho\mathbf u)}{\partial t}dV+\int_{A^*}\rho\mathbf u\,(\mathbf u\cdot\mathbf n)\,dA=\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA$ ·
     *why:* Replace each integral of step 2 by its box twin from step 4 — nothing else changes. · *plain:* The balloon's
     momentum law, written on the box.
  6. *did:* Substitute the box's theorem · *tex:* $\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV-\int_{A^*}\rho\mathbf u\,(\mathbf b\cdot\mathbf n)\,dA+\int_{A^*}\rho\mathbf u\,(\mathbf u\cdot\mathbf n)\,dA=\ldots$ ·
     *why:* Step 3 expresses the local integral through the box's total rate and the swept momentum; the right side is
     unchanged. · *plain:* Now the box's own momentum rate appears.
  7. *did:* Merge the two surface integrals · *tex:* $\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV+\int_{A^*}\rho\mathbf u\,(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA$ ·
     *why:* Same surface, same ρu: combine (u·n) − (b·n) = (u − b)·n. This is (4.17): momentum carried across the walls
     relative to them. · *plain:* Box momentum rate + momentum carried out relative to the walls = forces. · *live:* "wake
     box, x: 0 + (232.094 + 4.254 − 240.000) = −3.652 N/m = force on the fluid".
  8. *did:* Check the special boxes · *tex:* $\mathbf b=\mathbf u\Rightarrow\text{(4.13)};\qquad\mathbf b=0,\ \text{steady}:\ \oint_{A^*}\rho\mathbf u\,(\mathbf u\cdot\mathbf n)\,dA=\int_{V^*}\rho\mathbf g\,dV+\oint_{A^*}\mathbf f\,dA$ ·
     *why:* A material box gives back Newton for the balloon; a fixed box in steady flow reads forces from fluxes alone —
     the wake, jet and sprinkler cases. · *plain:* Momentum out minus in equals the forces on the fluid.
- **Result.** $\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV+\int_{A^*}\rho\mathbf u(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA$ (4.17) —
  *in words:* forces from fluxes, for any box.
- **Check.** Units: N for every term ✓. b = u gives (4.13) ✓. Wake (steady, fixed box, pressure p∞ on every face): out −
  in = −3.652 N/m = −F_D/l, so F_D/l = 3.652 N/m ✓ (C04 code). Bore with b = U: storage 0 ✓ (E1).
- **What it means.** The engineer's momentum tool: wake drag (Ex. 4.1), bores (Ex. 4.3), rockets (Ex. 4.4), the momentum
  integral of Ch. 9, lift and drag of Ch. 14, shock relations of Ch. 15. It gives the force on the *fluid*; the body feels
  the opposite.
- **Traps.** The book's "= 0" in (4.15). Treating ρu(u − b)·n as a dot product of two vectors (it is a vector times a
  scalar). Forgetting that f includes pressure. Taking the force on the body with the wrong sign (Newton III).

### D06 · Bernoulli along a streamline from a stream-tube element (Ex. 4.2) — ★★, 11 steps, in C05 (notebook · `which_bernoulli`)
- **Goal.** Derive $\tfrac12U^2+gz+p/\rho=$ constant along a streamline (4.19) from the mass and momentum budgets of a thin
  piece of stream tube, filling in the steps the book skips.
- **Start.** A fixed element of a stream tube, length ds along a streamline s: inlet area A, speed U, pressure p; outlet
  $A+\frac{\partial A}{\partial s}ds$, $U+\frac{\partial U}{\partial s}ds$, $p+\frac{\partial p}{\partial s}ds$; the tube axis rises by
  sin θ ds — *in words:* a short tapering tube whose side is made of streamlines.
- **Plan.** (1) Specialise (4.5) and (4.17) to a fixed element in steady, frictionless, constant-density flow. (2) Face by
  face: mass balance, momentum flux, forces (including the side push). (3) Simplify with the mass balance, drop small
  terms, divide. (4) Recognise exact differentials and integrate.
- **Tools.** (4.5) and (4.17) (C01, C04) · first-order Taylor expansion (ch01 P26) · orders of smallness (ch02 P68) ·
  side pressure force on a slowly widening tube (gloss in step 5) · chain rule · sympy expand/series (primer in C05) for
  the check.
- **Assumptions.** Steady (step 1: storage terms vanish); inviscid (step 1: only pressure acts on surfaces); constant ρ
  (steps 1, 9, 11); fixed element b = 0 (step 1); faces small enough for uniform U, p across each (step 2).
- **Steps.**
  1. *did:* Specialise the two box laws · *tex:* $\oint\rho\mathbf u\cdot\mathbf n\,dA=0,\qquad\oint\rho\mathbf u\,(\mathbf u\cdot\mathbf n)\,dA=\int\rho\mathbf g\,dV-\oint p\,\mathbf n\,dA$ ·
     *why:* (4.5) and (4.17) with b = 0 (fixed element), no time derivative (steady), f = −pn (inviscid: no shear). ·
     *plain:* What enters leaves; momentum out minus in equals gravity plus pressure. · *set:* {sc: 'pitot'} with the
     element drawn.
  2. *did:* Evaluate u·n on each face · *tex:* $\mathbf u\cdot\mathbf n\,dA=-U\,dA\ (\text{inlet}),\quad0\ (\text{side}),\quad+\Big(U+\frac{\partial U}{\partial s}ds\Big)dA\ (\text{outlet})$ ·
     *why:* The side is made of streamlines, so u is tangent to it; inlet and outlet are perpendicular to u. First-order
     Taylor gives the outlet values. · *plain:* Fluid enters at one end, leaves at the other, never through the side.
  3. *did:* Write the mass balance · *tex:* $-\rho UA+\rho\Big(U+\frac{\partial U}{\partial s}ds\Big)\Big(A+\frac{\partial A}{\partial s}ds\Big)=0$ ·
     *why:* Insert step 2 into the first law of step 1, with uniform speed over each small face. · *plain:* The mass flux
     at the outlet equals that at the inlet.
  4. *did:* Write the streamwise momentum flux · *tex:* $-\rho U^2A+\rho\Big(U+\frac{\partial U}{\partial s}ds\Big)^2\Big(A+\frac{\partial A}{\partial s}ds\Big)$ ·
     *why:* Momentum flux = velocity × mass flux (primer in C04); along s each face carries its speed times step 2. ·
     *plain:* Momentum leaving minus momentum entering.
  5. *did:* Add up the streamwise forces · *tex:* $-\rho g\sin\theta\Big(A+\frac{\partial A}{\partial s}\frac{ds}{2}\Big)ds+pA+\Big(p+\frac{\partial p}{\partial s}\frac{ds}{2}\Big)\frac{\partial A}{\partial s}ds-\Big(p+\frac{\partial p}{\partial s}ds\Big)\Big(A+\frac{\partial A}{\partial s}ds\Big)$ ·
     *why:* Gravity on the element (mean area × ds); pressure pushes the inlet forward, the outlet back; the book skips
     why: the widening side's mean pressure times its projected area (∂A/∂s)ds pushes forward. · *plain:* Weight, the
     two end pushes, and the side push. · *watch:* "the conical wall lights up".
  6. *did:* Simplify the flux with the mass balance · *tex:* $\rho\Big(U+\frac{\partial U}{\partial s}ds\Big)UA-\rho U^2A=\rho UA\,\frac{\partial U}{\partial s}ds$ ·
     *why:* By step 3 the outlet mass flux is ρUA, so the outlet momentum flux is (U + U_s ds) × ρUA; subtract the inlet
     ρU²A. · *plain:* The momentum gained is the mass flux times the speed gained.
  7. *did:* Expand the forces, use sin θ ds = dz · *tex:* $-\rho g\Big(A+\frac{\partial A}{\partial s}\frac{ds}2\Big)dz+\frac{\partial p}{\partial s}\frac{\partial A}{\partial s}\frac{(ds)^2}2-A\frac{\partial p}{\partial s}ds-\frac{\partial p}{\partial s}\frac{\partial A}{\partial s}(ds)^2$ ·
     *why:* Multiplying out step 5, pA and p(∂A/∂s)ds each appear with both signs and cancel; sin θ ds is the height gained,
     dz. · *plain:* What is left of the forces.
  8. *did:* Drop the second-order terms · *tex:* $\rho UA\,\frac{\partial U}{\partial s}ds=-\rho gA\,dz-A\frac{\partial p}{\partial s}ds$ · *why:*
     Terms with (ds)² or ds dz shrink faster than ds as the element shrinks and vanish after dividing by ds (orders of
     smallness, P68). · *plain:* Only first-order changes remain.
  9. *did:* Divide by ρA · *tex:* $U\frac{\partial U}{\partial s}ds=-g\,dz-\frac1\rho\frac{\partial p}{\partial s}ds$ · *why:* A > 0 and ρ is
     constant; every term is now per unit mass. · *plain:* Speeding up is paid for by falling or by a pressure drop.
  10. *did:* Recognise exact differentials · *tex:* $d\big(\tfrac12U^2\big)+g\,dz+\frac{dp}{\rho}=0$ · *why:* Along the streamline
      U(∂U/∂s)ds = d(U²/2) (chain rule) and (∂p/∂s)ds = dp; with ρ constant dp/ρ = d(p/ρ). · *plain:* The changes of three
      quantities add to zero. · *live:* "½(28.87)² + 0 + p∞/ρ = 416.7 + 83333.3 = 83750.0 m²/s²".
  11. *did:* Integrate along the streamline · *tex:* $\tfrac12U^2+gz+p/\rho=\text{constant along a streamline}$ · *why:* A sum of
      differentials that is zero means the sum itself does not change along the path (fundamental theorem of calculus).
      This is (4.19). · *plain:* Speed, height and pressure trade off along a streamline.
- **Result.** $\tfrac12U^2+gz+p/\rho=$ a constant along a streamline (4.19) — *in words:* in steady, frictionless,
  constant-density flow, kinetic energy, potential energy and "pressure energy" per unit mass trade along each streamline.
- **Check.** Units m²/s² ✓. U = 0: p + ρgz constant — hydrostatics (1.8) ✓. Pitot: 500 Pa ↔ 28.87 m/s ✓. `check_src`
  (optional, ★★): `ch04.stream_tube_element_balance_sym()` expands steps 3–7 symbolically, keeps O(ds) with
  `series(...).removeO()` and prints `U*dU + g*dz + dp/rho` = 0.
- **What it means.** The engineer's Bernoulli: pitot tubes, orifices, pressure on bodies. Its constant belongs to one
  streamline; C11 shows when it is the same on all of them. It fails with friction (behind a cylinder), unsteadiness (C12)
  and compressibility (C11's energy form).
- **Traps.** Forgetting the side pressure force (the result then has the wrong area factor). Using the inlet area for the
  weight (it is the mean area — the difference is second order anyway). The book's Ex. 4.2 statement writes
  "½ρU² + gz + p/ρ" — dimensionally inconsistent; the result is (4.19). Applying the constant across streamlines.

### D07 · Cauchy's equation: (4.14) → (4.20a, b) → (4.21) → (4.22) → (4.23) → (4.24) — ★★, 10 steps, in C06 (notebook)
- **Goal.** Newton's law for a single fluid particle, valid for any continuum, before we know how stress depends on the
  motion.
- **Start.** (4.14) for component j: $\int_V\frac{\partial(\rho u_j)}{\partial t}dV+\int_A\rho u_j(u_in_i)\,dA=\int_V\rho g_j\,dV+\int_Af_j\,dA$ — *in
  words:* the momentum balance of a material volume.
- **Plan.** (1) Turn both surface integrals into volume integrals with Gauss (first index!). (2) Localise. (3) Expand the
  flux form and recognise continuity and D/Dt.
- **Tools.** Gauss' theorem per component (R02) · traction $f_j=n_i\tau_{ij}$ (2.15) (R03) · tensor divergence over the first
  index (primer in C06) · localisation lemma (N07) · product rule (ch01 P38) · material derivative (3.5).
- **Assumptions.** ρu_iu_j and τ_ij continuously differentiable (steps 2–3); every material volume (step 5); continuity
  holds (step 9).
- **Steps.**
  1. *did:* Start from the momentum law, component j · *tex:* $\int_V\frac{\partial(\rho u_j)}{\partial t}dV+\int_A\rho u_ju_in_i\,dA=\int_V\rho g_j\,dV+\int_Af_j\,dA$ ·
     *why:* (4.14) written with indices; j is a free index (one equation per direction), i is summed in u_in_i. · *plain:*
     Momentum in direction j of a material volume.
  2. *did:* Apply Gauss to the momentum flux · *tex:* $\int_A\rho u_ju_in_i\,dA=\int_V\frac{\partial}{\partial x_i}(\rho u_iu_j)\,dV$ · *why:*
     For fixed j, Q_i = ρu_iu_j is a vector field; Gauss (2.30) contracts its index i with n_i. This is (4.20a). · *plain:*
     Momentum carried out through the skin = its flux divergence inside.
  3. *did:* Apply Gauss to the surface force · *tex:* $\int_Af_j\,dA=\int_An_i\tau_{ij}\,dA=\int_V\frac{\partial\tau_{ij}}{\partial x_i}dV$ · *why:*
     Traction f_j = n_iτ_ij (R03, first index = face); for fixed j, Q_i = τ_ij is a vector and Gauss differentiates the
     first index. This is (4.20b). · *plain:* The pushes on the skin add up to the stress divergence inside.
  4. *did:* Collect under one integral · *tex:* $\int_V\Big\{\frac{\partial(\rho u_j)}{\partial t}+\frac{\partial(\rho u_iu_j)}{\partial x_i}-\rho g_j-\frac{\partial\tau_{ij}}{\partial x_i}\Big\}dV=0$ ·
     *why:* All four terms are now integrals over the same V (linearity). This is (4.21). · *plain:* One integrand, zero
     integral.
  5. *did:* Localise · *tex:* $\frac{\partial(\rho u_j)}{\partial t}+\frac{\partial(\rho u_iu_j)}{\partial x_i}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ · *why:*
     The integral vanishes for every material volume and the integrand is continuous, so it vanishes everywhere
     (localisation lemma, D02). This is the flux form (4.22). · *plain:* Momentum conservation at a point, in flux form.
  6. *did:* Expand the time derivative · *tex:* $\frac{\partial(\rho u_j)}{\partial t}=\rho\frac{\partial u_j}{\partial t}+u_j\frac{\partial\rho}{\partial t}$ · *why:*
     Product rule (P38). We are separating "velocity changes" from "density changes". · *plain:* Momentum changes because
     velocity or density changes.
  7. *did:* Expand the flux divergence · *tex:* $\frac{\partial(\rho u_iu_j)}{\partial x_i}=u_j\frac{\partial(\rho u_i)}{\partial x_i}+\rho u_i\frac{\partial u_j}{\partial x_i}$ · *why:*
     Product rule with the two factors ρu_i and u_j. · *plain:* The flux changes because the mass flux spreads or because
     the velocity varies along the flow.
  8. *did:* Group the terms multiplying u_j · *tex:* $\rho\frac{\partial u_j}{\partial t}+u_j\Big[\frac{\partial\rho}{\partial t}+\frac{\partial(\rho u_i)}{\partial x_i}\Big]+\rho u_i\frac{\partial u_j}{\partial x_i}$ ·
     *why:* Adding steps 6 and 7 and factoring u_j out of two terms. This is the middle of (4.23). · *plain:* A bracket
     appears that we have seen before.
  9. *did:* Drop the bracket: it is continuity · *tex:* $\frac{\partial(\rho u_j)}{\partial t}+\frac{\partial(\rho u_iu_j)}{\partial x_i}=\rho\Big(\frac{\partial u_j}{\partial t}+u_i\frac{\partial u_j}{\partial x_i}\Big)$ ·
     *why:* The bracket is the continuity equation (4.7), $\frac{\partial\rho}{\partial t}+\frac{\partial}{\partial x_i}(\rho u_i)=0$, which holds at every point. · *plain:* With mass
     conserved, only the velocity changes remain.
  10. *did:* Recognise the particle acceleration · *tex:* $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ · *why:* (3.5) with
      F = u_j: ∂u_j/∂t + u_i∂u_j/∂x_i = Du_j/Dt; put into (4.22). This is (4.24). · *plain:* Mass × acceleration = weight +
      net push of the neighbours.
- **Result.** $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ (4.24) — *in words:* Cauchy's equation of motion.
- **Check.** Units N/m³ ✓. Fluid at rest, τ = −pδ: 0 = ρg − ∇p — hydrostatics ✓. Solid-body rotation (C06 worked
  example): −100 = 0 − 100 N/m³ radially ✓.
- **What it means.** Newton's law for water, air, honey or steel; 6 equations for 13 unknowns until a constitutive law
  (C07) ties τ to the motion. The flux form (4.22) is what finite-volume codes conserve.
- **Traps.** Differentiating the second index of τ (the book's prose after (4.24) writes ∂τ_ij/∂x_j — harmless only
  because τ turns out symmetric, proved after (4.24)). Missing the continuity bracket (subtract u_j × continuity, not
  u_j × ∇·u).

### D08 · The stress tensor is symmetric, (4.25) — ★★, 8 steps, in C07 (notebook · `newtonian_stress_lab`)
- **Goal.** Show τ₁₂ = τ₂₁ (and its cyclic partners): an unequal pair would spin a tiny cube infinitely fast. The book
  states the result and leaves the proof to Exercise 4.30.
- **Start.** A cube of fluid of side h centred at the origin, stresses τ_ij at its centre — *in words:* the smallest
  piece of fluid we can apply the rotational form of Newton's law to.
- **Plan.** (1) Torque of the shear pairs about the x₃ axis. (2) Show every other torque is of higher order. (3) Moment of
  inertia. (4) Angular acceleration as h → 0.
- **Tools.** Torque, moment arm and moment of inertia (primer in C04) · $d\mathbf H/dt=\mathbf M$ (4.64) · limits and orders of
  smallness (ch02 P68).
- **Assumptions.** No body couples (torques per unit mass) (step 8); stresses smooth, so they vary by O(h) across the cube
  (step 3).
- **Steps.**
  1. *did:* Take a small cube and an axis · *tex:* $dV=h^3,\qquad\text{axis }x_3\text{ through the centre}$ · *why:* Symmetry is a
     statement about a point, so we study an element and let it shrink; any axis will do, we choose x₃. · *plain:* A tiny
     cube, spinning (or not) about a vertical line. · *set:* {mode: 'cube', imb: 1, h: 0.01}.
  2. *did:* Torque of the shear pairs · *tex:* $M_3=\tau_{12}h^2\cdot h-\tau_{21}h^2\cdot h=(\tau_{12}-\tau_{21})h^3$ · *why:* On the faces
     x₁ = ±h/2 (area h²) the tractions ±τ₁₂ act along x₂ with arms ±h/2: two equal counterclockwise torques τ₁₂h³/2. The
     faces x₂ = ±h/2 carry ±τ₂₁ along x₁: clockwise. · *plain:* The two shear pairs twist the cube in opposite senses.
  3. *did:* Show other torques are smaller · *tex:* $M_3=(\tau_{12}-\tau_{21})h^3+O(h^4)$ · *why:* Normal stresses and gravity act
     through the face or mass centres (no arm); their variation across the cube, O(h), adds torques O(h⁴) and O(h⁵). ·
     *plain:* Only the shear imbalance matters for a small cube.
  4. *did:* Compute the moment of inertia · *tex:* $I_{33}=\int\rho(x_1^2+x_2^2)\,dV=\rho\frac{h^5}{6}$ · *why:* Each of the two
     integrals is ρh²·∫x²dx over [−h/2, h/2] = ρh²·h³/12; added: ρh⁵/6 (primer). · *plain:* The cube's resistance to spinning
     shrinks like h⁵.
  5. *did:* Apply the rotational law · *tex:* $I_{33}\,\dot\Omega_3=M_3$ · *why:* Rate of change of angular momentum = torque
     ((4.64), $d\mathbf H/dt=\mathbf M$) about the axis through the centre of mass. · *plain:* Torque = moment of inertia ×
     angular acceleration.
  6. *did:* Solve for the angular acceleration · *tex:* $\dot\Omega_3=\frac{6(\tau_{12}-\tau_{21})}{\rho h^2}+O\Big(\frac1h\Big)$ · *why:* Divide
     step 3 by step 4 (exponent rules: h³/h⁵ = h⁻²). · *plain:* The spin-up grows like 1/h². · *live:* "6 × 1/(1000 × 10⁻⁴) =
     60 rad/s² at h = 1 cm".
  7. *did:* Let the cube shrink · *tex:* $h\to0:\quad\dot\Omega_3\to\infty\ \text{unless}\ \tau_{12}=\tau_{21}$ · *why:* A fluid element cannot
     spin up infinitely fast; the only finite limit is a zero imbalance (orders of smallness, P68). · *plain:* An unequal
     pair would make every point whirl without bound. · *set:* {h: 0.001} · *watch:* "α jumps ×100".
  8. *did:* Repeat for the other axes · *tex:* $\tau_{ij}=\tau_{ji}$ · *why:* Relabelling the axes gives τ₂₃ = τ₃₂ and τ₃₁ = τ₁₃. With body couples (torque per unit mass) a term ∝ h³ could balance — the book's exception. This is (4.25). · *plain:* The stress tensor has six independent components.
- **Result.** $\tau_{ij}=\tau_{ji}$ (4.25) — *in words:* the stress tensor is symmetric.
- **Check.** Units: Pa/((kg/m³)·m²) = (kg m⁻¹ s⁻²)/(kg m⁻¹) = 1/s² ✓. Numbers: 1 Pa imbalance, water, h = 1 cm → 60 rad/s², h = 0.1 mm → 6×10⁵
  rad/s² (C07 code, slope −2) ✓.
- **What it means.** Six stress components, not nine; the traction rule f_j = n_iτ_ij can be read with either index;
  internal torques cancel in the angular-momentum budget (N83). Broken only by body couples (ferrofluids in fields).
- **Traps.** Including normal stresses in the torque (they have no arm about the centre). Forgetting the 1/12 in the
  moment of inertia. Arguing "equilibrium" — the argument works for accelerating fluid too.

### D09 · The Newtonian constitutive law: (4.28) → (4.29) → (4.30) → (4.31) — ★★★, 13 steps + sympy, in C07 (notebook · `newtonian_stress_lab`)
- **Goal.** Find the most general viscous stress that is linear in the strain rate and the same in every direction, and
  show it has only two material constants.
- **Start.** (4.27): $\tau_{ij}=-p\,\delta_{ij}+\sigma_{ij}$ with σ = 0 at rest — *in words:* total stress = pressure part +
  viscous part.
- **Plan.** (1) σ depends on S only, linearly: 81 coefficients. (2) Isotropy: 3 coefficients. (3) Contract with a
  symmetric S: only λ and μ + γ survive. (4) Name the constants, add the pressure, and check against Newton's law.
- **Tools.** R05 (only S enters) · isotropic fourth-order tensor (primer in C07, cited form checked numerically) · δ
  substitution in products (gloss) · symmetry of S (Ch. 3 §3.4) · np.einsum (ch02 P62) for the numbers.
- **Assumptions.** Galilean invariance and no response to rigid rotation (step 1); linearity (step 2); isotropy (step
  3); no memory (the stress depends on the present S only).
- **Steps.**
  1. *did:* Keep only the strain rate · *tex:* $\sigma_{ij}=\sigma_{ij}(\mathbf S),\qquad\sigma_{ij}=0\ \text{when}\ \mathbf S=0$ · *why:* R05: the
     stress cannot depend on u (Galilean invariance) nor on the rigid spin R; S is what deforms an element. · *plain:*
     Only deformation rates can create viscous stress. · *set:* {mode: 'stress', flow: 'shear', k: 10}.
  2. *did:* Assume the most general linear law · *tex:* $\sigma_{ij}=K_{ijmn}S_{mn}$ · *why:* Linear is the simplest law that
     vanishes at S = 0 (a first-order expansion); each of 9 σ's may depend on each of 9 S's: 81 coefficients. This is (4.28).
     · *plain:* Each stress component is a weighted sum of all strain-rate components.
  3. *did:* Demand isotropy · *tex:* $K'_{ijmn}=K_{ijmn}\ \text{in every rotated frame}$ · *why:* The fluid has no preferred direction,
     so the law must read the same in every orientation (Ch. 2 §2.5). · *plain:* Turning the axes must not change the law.
  4. *did:* Use the isotropic form · *tex:* $K_{ijmn}=\lambda\delta_{ij}\delta_{mn}+\mu\delta_{im}\delta_{jn}+\gamma\delta_{in}\delta_{jm}$ · *why:* The only
     isotropic fourth-order tensors are these products of two δ's (primer; cited, checked by 50 random rotations). This is
     (4.29). · *plain:* 81 numbers collapse to three: λ, μ, γ. · *watch:* "three colours appear in the τ grid".
  5. *did:* Contract the first term · *tex:* $\lambda\,\delta_{ij}\delta_{mn}S_{mn}=\lambda S_{mm}\delta_{ij}$ · *why:* δ_mnS_mn replaces n by m and
     sums the diagonal: S_mm = ∇·u (δ substitution). · *plain:* The first part responds to volume change, equally in every
     direction.
  6. *did:* Contract the second term · *tex:* $\mu\,\delta_{im}\delta_{jn}S_{mn}=\mu S_{ij}$ · *why:* δ_im turns m into i and δ_jn turns n
     into j (δ substitution). · *plain:* The second part copies the strain rate.
  7. *did:* Contract the third term · *tex:* $\gamma\,\delta_{in}\delta_{jm}S_{mn}=\gamma S_{ji}=\gamma S_{ij}$ · *why:* δ_in: n → i; δ_jm: m → j,
     giving S_ji; the strain rate is symmetric (Ch. 3), so S_ji = S_ij. · *plain:* The third part also copies the strain
     rate. · *live:* "S₁₂ = 5 s⁻¹: τ₁₂ gets (μ + γ) × 5".
  8. *did:* Add the three parts · *tex:* $\sigma_{ij}=\lambda S_{mm}\delta_{ij}+(\mu+\gamma)S_{ij}$ · *why:* Linearity of the contraction.
     Notice σ is already symmetric, whatever γ is. · *plain:* Only two combinations act: λ and μ + γ.
  9. *did:* Name the sum 2μ · *tex:* $\gamma=\mu\quad\Rightarrow\quad\mu+\gamma=2\mu$ · *why:* The book argues that K symmetric in i, j forces γ = μ (4.30). Here only μ + γ acts, so γ = μ just names the sum 2μ — making μ Newton's viscosity. ·
     *plain:* Two constants remain. · *watch:* "only μ + γ acts on a symmetric S".
  10. *did:* Write the viscous stress · *tex:* $\sigma_{ij}=2\mu S_{ij}+\lambda S_{mm}\delta_{ij}$ · *why:* Substitute step 9 into step 8. ·
      *plain:* Viscous stress = 2μ × strain rate + λ × volume-change rate on the diagonal.
  11. *did:* Add the pressure part · *tex:* $\tau_{ij}=-p\,\delta_{ij}+2\mu S_{ij}+\lambda S_{mm}\delta_{ij}$ · *why:* (4.27),
      $\tau_{ij}=-p\delta_{ij}+\sigma_{ij}$, with p the thermodynamic pressure. This is (4.31). · *plain:* The Newtonian stress
      law. · *set:* {flow: 'expansion'} (the λS_mm term appears on the diagonal).
  12. *did:* Check against Newton's viscosity law · *tex:* $\mathbf u=(u(y),0,0):\quad\tau_{12}=2\mu S_{12}=\mu\frac{du}{dy}$ · *why:* For a
      parallel flow S₁₂ = ½du/dy and S_mm = 0; we recover τ = μ du/dy (1.3), so μ is the ordinary viscosity. · *plain:* The
      general law contains Newton's. · *watch:* the plane rotates (the isotropy check).
  13. *did:* Count what remains · *tex:* $81\ \to\ 3\ \to\ 2:\qquad\mu,\ \lambda$ · *why:* Linearity gave 81, isotropy 3, symmetry
      of S 2; both may depend on the thermodynamic state (temperature) but not on the motion. · *plain:* Two numbers
      describe a Newtonian fluid's viscosity. · *live:* "τ₁₂ = 2 × 10⁻³ × 5 = 0.010 Pa; τ₁₁ + p = 0 in shear".
- **Result.** $\tau_{ij}=-p\,\delta_{ij}+2\mu S_{ij}+\lambda S_{mm}\delta_{ij}$ (4.31) — *in words:* a fluid that is linear and isotropic
  needs only μ (shear) and λ (volume change).
- **Check.** Units: Pa s × 1/s = Pa ✓. Rigid rotation: S = 0 ⇒ σ = 0 ✓. Parallel shear gives (1.3) ✓. Numbers: water,
  γ̇ = 10 s⁻¹ → τ₁₂ = 0.010 Pa ✓. **sympy check (`check_src`, every line commented):** build K from `sp.KroneckerDelta` with
  symbols λ, μ, γ (step 4); a symmetric S from six symbols; contract (steps 5–7) and assert σ_ij − [λS_mmδ_ij + (μ + γ)S_ij]
  simplifies to 0 for all nine components (step 8) and σ − σᵀ = 0; substitute γ = μ and compare with 2μS + λS_mmδ (steps
  9–10); build a rotation about z by a symbol θ and assert K′_ijmn − K_ijmn = 0 for 20 random index sets (step 4's
  isotropy) — the cell re-runs the construction, not just the result.
- **What it means.** The law used for water and air throughout the book; C08 puts it into Cauchy's equation. It fails for
  non-Newtonian fluids (R06) and at extreme strain rates (rarefied gases, shocks' interiors).
- **Traps.** Thinking γ = μ is forced by the symmetry of σ (σ is symmetric anyway; only μ + γ is physical). Treating λ as
  the bulk viscosity (μ_v = λ + ⅔μ, D10). Forgetting σ = 0 when S = 0 — that is why there is no constant term.

### D10 · The bulk-viscosity form (4.31) → (4.37), and the Stokes assumption (4.36) — ★, 4 steps, in C07 (notebook · `newtonian_stress_lab`)
- **Goal.** Separate the part of the viscous stress that changes a particle's shape from the part that changes its
  volume — the book states the result.
- **Start.** (4.31): $\tau_{ij}=-p\,\delta_{ij}+2\mu S_{ij}+\lambda S_{mm}\delta_{ij}$.
- **Plan.** (1) Add and subtract ⅔μS_mmδ_ij. (2) Group into a traceless part and an isotropic part. (3) Name μ_v. (4) Take
  the trace to see what μ_v does to the mean pressure.
- **Tools.** Deviatoric (traceless) part of a tensor (primer in C07) · δ_ii = 3 (Ch. 2 §2.1).
- **Assumptions.** None beyond (4.31).
- **Steps.**
  1. *did:* Add and subtract ⅔μS_mmδ_ij · *tex:* $\tau_{ij}=-p\delta_{ij}+2\mu S_{ij}-\tfrac23\mu S_{mm}\delta_{ij}+\tfrac23\mu S_{mm}\delta_{ij}+\lambda S_{mm}\delta_{ij}$ ·
     *why:* Adding zero changes nothing; we do it to build the traceless part of S, which is pure shape change. ·
     *plain:* The same stress, rearranged. · *set:* {flow: 'expansion', muv: 0.003}.
  2. *did:* Group shape and volume parts · *tex:* $\tau_{ij}=-p\delta_{ij}+2\mu\Big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\Big)+\Big(\lambda+\tfrac23\mu\Big)S_{mm}\delta_{ij}$ ·
     *why:* Factor 2μ from the first pair and S_mmδ_ij from the second; the bracket is the deviatoric part of S (primer):
     its trace is S_mm − ⅓S_mm·3 = 0. · *plain:* Shear-and-stretch at fixed volume, plus pure swelling. · *watch:* "the
     diagonal colours split into rose and purple".
  3. *did:* Name the bulk viscosity · *tex:* $\mu_v\equiv\lambda+\tfrac23\mu:\qquad\tau_{ij}=-p\delta_{ij}+2\mu\Big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\Big)+\mu_vS_{mm}\delta_{ij}$ ·
     *why:* The coefficient of the volume-change rate deserves its own name. This is (4.37). · *plain:* μ resists shape
     change, μ_v resists volume change. · *live:* "μ_v = λ + ⅔μ = … Pa s".
  4. *did:* Take the trace · *tex:* $\tau_{ii}=-3p+3\mu_vS_{mm}\ \Rightarrow\ p-\bar p=\mu_v\nabla\cdot\mathbf u$ · *why:* Set i = j and sum: δ_ii = 3
     and the deviatoric part has zero trace; with p̄ = −⅓τ_ii (4.33) this is (4.34). Stokes' assumption (4.36), μ_v = 0,
     makes p̄ = p. · *plain:* The mean normal stress differs from −p only in expansion, and only through μ_v.
- **Result.** $\tau_{ij}=-p\,\delta_{ij}+2\mu\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)+\mu_vS_{mm}\delta_{ij}$ (4.37), with $\mu_v=\lambda+\tfrac23\mu$ —
  *in words:* shape change costs μ, volume change costs μ_v.
- **Check.** Units ✓. Incompressible (S_mm = 0) → τ = −pδ + 2μS (4.35) ✓. Expansion e = 1 s⁻¹ with μ_v = 3×10⁻³ Pa s: every
  normal stress gains μ_v S_mm = 9 mPa, p − p̄ = 9 mPa (E3 parity row) ✓.
- **What it means.** Stokes' assumption λ = −⅔μ (μ_v = 0) is the code default and good almost everywhere; μ_v matters for
  sound absorption and inside shocks (Ch. 15).
- **Traps.** Calling S − ⅓S_mmδ "S". Setting λ = 0 instead of μ_v = 0 for the Stokes assumption. Believing σ is always
  traceless (only if μ_v = 0 or ∇·u = 0).

### D11 · The Navier–Stokes equation: (4.24) + (4.37) → (4.38) — ★★, 7 steps, in C08 (notebook · `navier_stokes_term_balance`)
- **Goal.** Put the Newtonian stress into Cauchy's equation and get the momentum equation of a Newtonian fluid, keeping
  a viscosity that may vary from place to place.
- **Start.** (4.24): $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ — *in words:* Newton's law for a particle, any continuum.
- **Plan.** (1) Write τ with velocity gradients. (2) Differentiate the pressure part. (3) Keep the viscous part under the
  derivative. (4) Write the acceleration out.
- **Tools.** (4.37) (C07) · strain rate $S_{ij}=\tfrac12\big(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\big)$ (3.12) (Ch. 3 §3.4) · δ
  substitution (gloss) · tensor divergence over the first index (primer in C06) · material derivative (3.5).
- **Assumptions.** Newtonian fluid (step 3); μ and μ_v may depend on position through the temperature (step 5).
- **Steps.**
  1. *did:* Start from Cauchy's equation · *tex:* $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ · *why:* (4.24) holds for
     any continuum; we now know τ for a Newtonian fluid (C07). · *plain:* Mass × acceleration = weight + net push. · *set:*
     {sol: 'poiseuille'}.
  2. *did:* Write 2μS with velocity gradients · *tex:* $2\mu S_{ij}=\mu\Big(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\Big)$ · *why:*
     Definition of the strain rate (3.12). · *plain:* Twice the strain rate is the sum of a gradient and its transpose.
  3. *did:* Write τ with velocity gradients · *tex:* $\tau_{ij}=-p\,\delta_{ij}+\mu\Big(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\Big)+\Big(\mu_v-\tfrac23\mu\Big)\frac{\partial u_m}{\partial x_m}\delta_{ij}$ ·
     *why:* (4.37): 2μ(S − ⅓S_mmδ) + μ_vS_mmδ = 2μS + (μ_v − ⅔μ)S_mmδ, with step 2 and S_mm = ∂u_m/∂x_m. · *plain:* The
     Newtonian stress in terms of the velocity field.
  4. *did:* Differentiate the pressure part · *tex:* $\frac{\partial}{\partial x_i}\big(-p\,\delta_{ij}\big)=-\frac{\partial p}{\partial x_j}$ · *why:* δ_ij is
     constant and picks the term with i = j (δ substitution; the book's "(∂p/∂x_i)δ_ij = ∂p/∂x_j"). · *plain:* Pressure
     pushes from high to low. · *watch:* "the orange bar lights".
  5. *did:* Keep μ inside the viscous derivative · *tex:* $\frac{\partial\sigma_{ij}}{\partial x_i}=\frac{\partial}{\partial x_i}\Big[\mu\Big(\frac{\partial u_j}{\partial x_i}+\frac{\partial u_i}{\partial x_j}\Big)+\Big(\mu_v-\tfrac23\mu\Big)\frac{\partial u_m}{\partial x_m}\delta_{ij}\Big]$ ·
     *why:* μ depends on temperature (it falls with T in liquids, rises in gases), so it may vary in space; we must not pull
     it outside the derivative yet. · *plain:* The net viscous force, exact for variable viscosity. · *watch:* "the rose bar
     lights".
  6. *did:* Write the acceleration out · *tex:* $\rho\frac{Du_j}{Dt}=\rho\Big(\frac{\partial u_j}{\partial t}+u_i\frac{\partial u_j}{\partial x_i}\Big)$ · *why:* Material derivative (3.5) with F = u_j: the local plus the advective rate of the velocity. · *plain:* Local plus advective acceleration.
  7. *did:* Assemble the equation · *tex:* $\rho\Big(\frac{\partial u_j}{\partial t}+u_i\frac{\partial u_j}{\partial x_i}\Big)=-\frac{\partial p}{\partial x_j}+\rho g_j+\frac{\partial}{\partial x_i}\Big[\mu\Big(\frac{\partial u_j}{\partial x_i}+\frac{\partial u_i}{\partial x_j}\Big)+\Big(\mu_v-\tfrac23\mu\Big)\frac{\partial u_m}{\partial x_m}\delta_{ij}\Big]$ ·
     *why:* Insert steps 4–6 into step 1. This is (4.38), the Navier–Stokes equation. · *plain:* Acceleration = pressure +
     gravity + viscous force. · *live:* "at the probe: 0 + 0 = +100 + 0 − 100 N/m³".
- **Result.** (4.38) as in step 7 — *in words:* the Navier–Stokes momentum equation for a Newtonian fluid with variable
  viscosity.
- **Check.** Units N/m³ ✓. μ = 0 gives Euler (4.41) ✓. Poiseuille: 0 = 100 + 0 − 100 ✓. Ledger: with (4.7), 4 equations for
  ρ, p, u_j (5 unknowns) — closed by ρ = const or ρ = ρ(p) ✓. sympy: `ch04.navier_stokes_sym` with μ = μ₀(1 + by) keeps a
  dμ/dy term that vanishes for b = 0 (C08 code).
- **What it means.** The equation of the rest of the book. The next derivation simplifies it for constant viscosity and
  incompressible flow.
- **Traps.** Writing ∂(pδ_ij)/∂x_i as ∂p/∂x_i (the free index is j). Pulling μ out when it varies (hot bearings, magma).
  Forgetting the transpose term ∂u_i/∂x_j.

### D12 · Constant viscosity, then incompressible: (4.38) → (4.39a) → (4.39b) — ★★, 7 steps, in C08 (notebook · `navier_stokes_term_balance`)
- **Goal.** Simplify Navier–Stokes when the viscosities are uniform and the flow is incompressible, and see exactly which
  term incompressibility removes — the book states the results.
- **Start.** The viscous term of (4.38): $\frac{\partial}{\partial x_i}\Big[\mu\Big(\frac{\partial u_j}{\partial x_i}+\frac{\partial u_i}{\partial x_j}\Big)+\Big(\mu_v-\tfrac23\mu\Big)\frac{\partial u_m}{\partial x_m}\delta_{ij}\Big]$.
- **Plan.** (1) Take the constant viscosities outside. (2) Swap a derivative order to reveal ∇(∇·u). (3) Collect the
  coefficients. (4) Set ∇·u = 0.
- **Tools.** Schwarz's theorem (primer in C08) · vector Laplacian in Cartesian components (gloss) · (4.10) (C02).
- **Assumptions.** μ, μ_v uniform (small temperature differences) from step 1; u twice continuously differentiable (step
  3); ∇·u = 0 from step 6.
- **Steps.**
  1. *did:* Take the viscosities outside · *tex:* $\frac{\partial\sigma_{ij}}{\partial x_i}=\mu\frac{\partial}{\partial x_i}\Big(\frac{\partial u_j}{\partial x_i}+\frac{\partial u_i}{\partial x_j}\Big)+\Big(\mu_v-\tfrac23\mu\Big)\frac{\partial}{\partial x_j}\frac{\partial u_m}{\partial x_m}$ ·
     *why:* Uniform μ, μ_v are constants for the derivative; δ_ij turns ∂/∂x_i into ∂/∂x_j in the last term. · *plain:* Two
     viscous pieces remain. · *set:* {sol: 'toy'} (a compressible field, so every piece is visible).
  2. *did:* Split the first bracket · *tex:* $\mu\frac{\partial^2u_j}{\partial x_i\partial x_i}+\mu\frac{\partial}{\partial x_i}\frac{\partial u_i}{\partial x_j}$ · *why:* Linearity of the derivative: the bracket's two terms are differentiated separately. · *plain:* A Laplacian and a mixed term.
  3. *did:* Swap the derivative order · *tex:* $\frac{\partial}{\partial x_i}\frac{\partial u_i}{\partial x_j}=\frac{\partial}{\partial x_j}\frac{\partial u_i}{\partial x_i}=\frac{\partial}{\partial x_j}(\nabla\cdot\mathbf u)$ ·
     *why:* Schwarz's theorem (primer): mixed partial derivatives of a smooth field commute. The sum over i is then the
     divergence. · *plain:* The mixed term is the gradient of the divergence. · *watch:* "a purple bar (μ_v + ⅓μ)∂_j(∇·u)
     appears".
  4. *did:* Collect the ∇(∇·u) coefficients · *tex:* $\Big(\mu+\mu_v-\tfrac23\mu\Big)\frac{\partial}{\partial x_j}\frac{\partial u_m}{\partial x_m}=\Big(\mu_v+\tfrac13\mu\Big)\frac{\partial}{\partial x_j}\frac{\partial u_m}{\partial x_m}$ ·
     *why:* The book skips this bookkeeping: μ from step 3 plus μ_v − ⅔μ from step 1 (renaming the dummy index i as m). ·
     *plain:* One coefficient in front of the gradient of the divergence.
  5. *did:* Write the constant-viscosity form · *tex:* $\rho\frac{Du_j}{Dt}=-\frac{\partial p}{\partial x_j}+\rho g_j+\mu\frac{\partial^2u_j}{\partial x_i^2}+\Big(\mu_v+\tfrac13\mu\Big)\frac{\partial}{\partial x_j}\frac{\partial u_m}{\partial x_m}$ ·
     *why:* Put steps 2–4 into (4.38); ∂²u_j/∂x_i∂x_i is the Laplacian of component j (vector Laplacian, gloss). This is
     (4.39a). · *plain:* Viscous diffusion plus a compressibility term.
  6. *did:* Impose incompressibility · *tex:* $\frac{\partial u_m}{\partial x_m}=\nabla\cdot\mathbf u=0$ · *why:* (4.10) for incompressible flow
     removes the last term entirely. · *plain:* The compressibility term vanishes. · *set:* {sol: 'poiseuille'} · *watch:*
     "the purple bar vanishes".
  7. *did:* Write it in vector form · *tex:* $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ · *why:* The three components
     together; ∇²u acts on each Cartesian component. This is (4.39b). · *plain:* The incompressible Navier–Stokes equation.
- **Result.** (4.39a), and $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ (4.39b) — *in words:* with uniform viscosity
  and no volume change, the viscous force is μ times the Laplacian of the velocity.
- **Check.** Units ✓. Poiseuille satisfies (4.39b) exactly ✓. sympy (C08): for u = (y(1 − y), 0, 0) the ∇(∇·u) term is 0
  ✓. The toy field u = (ax², 0, 0) has ∇·u = 2ax and a nonzero (μ_v + ⅓μ)·2a term ✓.
- **What it means.** Viscosity acts as diffusion of momentum (ν∇²u after dividing by ρ). The compressibility term matters
  only for sound and shocks.
- **Traps.** Coefficient bookkeeping (μ + μ_v − ⅔μ = μ_v + ⅓μ, not μ_v − ⅓μ). Pulling out μ when it varies. Using
  Cartesian "∇² of each component" in cylindrical coordinates (extra terms, Appendix B).

### D13 · The viscous force three ways, (4.40) — ★★, 9 steps, in C08 (notebook · `navier_stokes_term_balance`)
- **Goal.** Show $\mu\nabla^2\mathbf u=2\mu\,\partial S_{ij}/\partial x_i=-\mu\nabla\times\boldsymbol\omega$ for incompressible flow, and resolve the
  "paradox" that vorticity appears though rotation was excluded from the stress law. The book leaves it to Exercise 4.38.
- **Start.** $(\mu\nabla^2\mathbf u)_j=\mu\frac{\partial^2u_j}{\partial x_i\partial x_i}$ — *in words:* the viscous force of (4.39b), component j.
- **Plan.** (1) Show the strain-rate form equals it when ∇·u = 0. (2) Write ∇×ω with ε's, use the ε–δ identity, and
  compare. (3) Test on solid-body rotation.
- **Tools.** Strain rate (3.12) · Schwarz's theorem (primer in C08) · vorticity ω = ∇×u (Ch. 3 §3.4) · ε–δ identity
  $\varepsilon_{kij}\varepsilon_{klm}=\delta_{il}\delta_{jm}-\delta_{im}\delta_{jl}$ (2.19) and the curl-of-a-curl identity (primer in C08).
- **Assumptions.** Incompressible, ∇·u = 0 (steps 3, 8); u smooth (steps 3, 8).
- **Steps.**
  1. *did:* Start from the Laplacian form · *tex:* $(\mu\nabla^2\mathbf u)_j=\mu\frac{\partial^2u_j}{\partial x_i\partial x_i}$ · *why:* The viscous term
     of (4.39b); we compare it with two other expressions. · *plain:* Viscous force = μ × curvature of the velocity
     profile. · *set:* {sol: 'lamb_oseen'}.
  2. *did:* Differentiate twice the strain rate · *tex:* $2\frac{\partial S_{ij}}{\partial x_i}=\frac{\partial^2u_i}{\partial x_i\partial x_j}+\frac{\partial^2u_j}{\partial x_i\partial x_i}$ ·
     *why:* 2S_ij = ∂u_i/∂x_j + ∂u_j/∂x_i (3.12); differentiate each term by x_i and sum over i. · *plain:* The divergence of
     the strain rate has two parts.
  3. *did:* Swap and use ∇·u = 0 · *tex:* $2\frac{\partial S_{ij}}{\partial x_i}=\frac{\partial}{\partial x_j}(\nabla\cdot\mathbf u)+\nabla^2u_j=\nabla^2u_j$ · *why:*
     Schwarz's theorem moves ∂/∂x_j outside; the divergence is zero for incompressible flow. So μ∇²u = 2μ∂S_ij/∂x_i (also
     the book's μ∂/∂x_i(∂u_j/∂x_i + ∂u_i/∂x_j) form). · *plain:* The first two forms agree.
  4. *did:* Write the vorticity with ε · *tex:* $\omega_k=\varepsilon_{klm}\frac{\partial u_m}{\partial x_l}$ · *why:* ω = ∇×u (Ch. 3), the curl
     written with the alternating tensor (Ch. 2 §2.7). · *plain:* Vorticity in index form.
  5. *did:* Write the curl of the vorticity · *tex:* $(\nabla\times\boldsymbol\omega)_j=\varepsilon_{jik}\frac{\partial\omega_k}{\partial x_i}=\varepsilon_{jik}\varepsilon_{klm}\frac{\partial^2u_m}{\partial x_i\partial x_l}$ ·
     *why:* The curl of any vector is ε_jik∂_i(·)_k; insert step 4. · *plain:* A double sum of second derivatives. · *live:*
     "−μ(∇×ω)_x at the probe = … N/m³".
  6. *did:* Apply the ε–δ identity · *tex:* $\varepsilon_{jik}\varepsilon_{klm}=\delta_{jl}\delta_{im}-\delta_{jm}\delta_{il}$ · *why:* Cyclic shift
     ε_jik = ε_kji, then (2.19) with the first indices matched. · *plain:* Two ε's become two pairs of δ's.
  7. *did:* Substitute the δ's · *tex:* $(\nabla\times\boldsymbol\omega)_j=\frac{\partial^2u_i}{\partial x_i\partial x_j}-\frac{\partial^2u_j}{\partial x_i\partial x_i}$ · *why:*
     δ_jlδ_im sets l = j, m = i; δ_jmδ_il sets m = j, l = i (δ substitution). · *plain:* The curl of ω is the gradient of
     the divergence minus the Laplacian.
  8. *did:* Use ∇·u = 0 again · *tex:* $(\nabla\times\boldsymbol\omega)_j=-\nabla^2u_j\ \Rightarrow\ \mu\nabla^2\mathbf u=-\mu\nabla\times\boldsymbol\omega$ · *why:* The first
     term is ∂_j(∇·u) = 0 (Schwarz, incompressible). This is the last equality of (4.40). · *plain:* The viscous force is
     minus μ times the curl of the vorticity.
  9. *did:* Test solid-body rotation · *tex:* $\mathbf u=\boldsymbol\Omega\times\mathbf x:\quad\boldsymbol\omega=2\boldsymbol\Omega\ \text{uniform}\ \Rightarrow\ \mu\nabla^2\mathbf u=-\mu\nabla\times\boldsymbol\omega=0$ ·
     *why:* A rigid rotation has uniform vorticity (Ch. 3), whose curl is zero; its strain rate is zero too. The viscous
     force depends on the *derivative* of ω — the paradox dissolves. · *plain:* A rigidly spinning fluid feels no friction.
     · *set:* {sol: 'solid_body'} · *watch:* "all three arrows shrink to zero".
- **Result.** $(\mu\nabla^2\mathbf u)_j=2\mu\frac{\partial S_{ij}}{\partial x_i}=-\mu(\nabla\times\boldsymbol\omega)_j$ (4.40), for ∇·u = 0 — *in words:*
  the same viscous force from the Laplacian, the strain-rate gradient or the vorticity's curl.
- **Check.** Units N/m³ ✓. Lamb–Oseen vortex: the three forms agree to stencil accuracy (`viscous_force_forms`) ✓.
  Irrotational flow: ω = 0 ⇒ no viscous force (C12's N98) ✓. Poiseuille: ω = −du/dy e_z, −μ(∇×ω)_x = μu″ = −G ✓.
- **What it means.** Viscosity acts only where vorticity varies in space — it diffuses vorticity (Ch. 5) and is idle in
  irrotational flow, which is why potential flow (Ch. 6) is an exact Navier–Stokes solution away from walls.
- **Traps.** Using the identity when ∇·u ≠ 0 (the ∇(∇·u) term then stays). Index order in ε_jik (a swap flips the sign).
  Concluding "no vorticity, no viscous stress" (stress is not force: Couette flow has stress but no net viscous force).

### D14 · Velocity in a translating–rotating frame, (4.42) — ★★, 7 steps, in C09 (notebook · `rotating_frame_coriolis`)
- **Goal.** Relate the velocity a rotating observer measures to the true (inertial) velocity, and find how the rotating
  axes themselves change in time.
- **Start.** $\mathbf x=\mathbf X(t)+\mathbf x'$ — *in words:* the particle's position = position of the moving origin O′ + its
  position relative to O′ (Fig. 4.6).
- **Plan.** (1) Write x′ with the rotating basis. (2) Find de′_i/dt from the cone construction. (3) Differentiate with
  the product rule and name the pieces.
- **Tools.** Derivative of a rotating unit vector (primer in C09) · product rule (ch01 P38) · cross product and
  right-hand rule (ch02 P74).
- **Assumptions.** The rotating frame is rigid: its axes turn together at Ω(t) (step 3).
- **Steps.**
  1. *did:* Split the position · *tex:* $\mathbf x=\mathbf X(t)+\mathbf x'$ · *why:* Vector addition (Fig. 4.6): O′ is at X(t); the
     particle is at x′ from O′. · *plain:* Where the particle is = where the frame is + where it is in the frame. · *set:*
     {mode: 'turntable', Om: 0.5, t: 0}.
  2. *did:* Write x′ in the rotating basis · *tex:* $\mathbf x'=x'_i\,\mathbf e'_i(t)$ · *why:* The rotating observer measures components
     x′_i along its own axes e′_i — and those axes turn, so they are functions of time. This is the whole difference from
     Ch. 3's Galilean frame. · *plain:* The components and the axes both change.
  3. *did:* Size of a basis vector's change · *tex:* $\Big\lvert\frac{d\mathbf e'_1}{dt}\Big\rvert=\sin\alpha\,\lvert\boldsymbol\Omega\rvert$ · *why:* The cone
     construction (primer, Fig. 4.7): e′₁'s tip moves on a circle of radius sin α round the Ω axis at angular speed |Ω|. ·
     *plain:* A unit vector sweeps a cone. · *set:* {t: 1} (e′₁ sweeps its cone).
  4. *did:* Identify its direction · *tex:* $\frac{d\mathbf e'_i}{dt}=\boldsymbol\Omega\times\mathbf e'_i$ · *why:* The change is perpendicular to
     both Ω and e′₁ with size sin α|Ω| — exactly the cross product (right-hand rule); the same for every axis. · *plain:*
     Turning axes change at Ω × axis. · *live:* "|de′₁/dt| = 1 × 0.5 = 0.5 s⁻¹".
  5. *did:* Differentiate the position · *tex:* $\mathbf u=\frac{d\mathbf X}{dt}+\frac{dx'_i}{dt}\mathbf e'_i+x'_i\frac{d\mathbf e'_i}{dt}$ · *why:* Product
     rule on each term x′_ie′_i (both factors depend on t). · *plain:* Three contributions to the true velocity.
  6. *did:* Name the first two pieces · *tex:* $\frac{d\mathbf X}{dt}=\mathbf U,\qquad\frac{dx'_i}{dt}\mathbf e'_i=\mathbf u'$ · *why:* U is the frame's
     velocity; u′ is what the rotating observer measures — the rates of change of its own components. · *plain:* Frame
     velocity and velocity seen in the frame.
  7. *did:* Use the turning rule · *tex:* $\mathbf u=\mathbf U+\mathbf u'+\boldsymbol\Omega\times\mathbf x'$ · *why:* x′_i(Ω × e′_i) = Ω × (x′_ie′_i) = Ω × x′
     (linearity of ×). This is (4.42). · *plain:* True velocity = frame velocity + measured velocity + the frame's rotation
     carrying the particle.
- **Result.** $\mathbf u=\mathbf U+\mathbf u'+\boldsymbol\Omega\times\mathbf x'$ with $\frac{d\mathbf e'_i}{dt}=\boldsymbol\Omega\times\mathbf e'_i$ (4.42) — *in words:*
  a rotating observer misses the part Ω × x′ of the velocity.
- **Check.** Units m/s ✓. Ω = 0: Galilean u = U + u′ (Ch. 3) ✓. A particle at rest on the turntable (u′ = 0) moves at
  Ω × x′ in the inertial frame ✓. `basis_rate` equals Ω × e′ numerically (primer code) ✓.
- **What it means.** The rotating observer's velocity differs from the true one by a solid-body rotation — why vorticity
  seen from the Earth is ζ, not ζ + f (Ch. 3's ω′ = ω − 2Ω, Ch. 13).
- **Traps.** Treating e′_i as constant (that is the Galilean case). Using the inertial position in Ω × x′. Thinking u′ is
  the true velocity "in rotating components" (it is a different vector).

### D15 · Acceleration in a noninertial frame: (4.42) → (4.43) → (4.44) — ★★★, 12 steps + sympy, in C09 (notebook · `rotating_frame_coriolis`)
- **Goal.** Differentiate (4.42) once more to find the true acceleration in terms of what the rotating observer measures —
  and see where the famous factor 2 of the Coriolis term comes from. The book states the result (Exercise 4.42).
- **Start.** (4.42): $\mathbf u=\mathbf U+\mathbf u'+\boldsymbol\Omega\times\mathbf x'$ — *in words:* the true velocity from D14.
- **Plan.** (1) Differentiate term by term. (2) The u′ term: product rule on its rotating basis — one Ω × u′. (3) The
  Ω × x′ term: product rule for a cross product, and D14 for dx′/dt — a second Ω × u′ and the centripetal term. (4) Collect,
  then follow a fluid particle.
- **Tools.** D14 (the turning rule) · product rule for a cross product (primer in C09) · vector triple product (gloss,
  Ch. 2 §2.7) · material derivative in each frame.
- **Assumptions.** U and Ω depend on time only (uniform in space); the frame is rigid.
- **Steps.**
  1. *did:* Start from the velocity · *tex:* $\mathbf u=\mathbf U+\mathbf u'+\boldsymbol\Omega\times\mathbf x'$ · *why:* (4.42), D14's result; the true
     acceleration is its time derivative. · *plain:* The velocity in three pieces. · *set:* {mode: 'turntable', Om: 0.5,
     u0: 3, t: 0.5}.
  2. *did:* Differentiate term by term · *tex:* $\mathbf a=\frac{d\mathbf U}{dt}+\frac{d\mathbf u'}{dt}+\frac{d}{dt}\big(\boldsymbol\Omega\times\mathbf x'\big)$ · *why:*
     The derivative of a sum is the sum of the derivatives. · *plain:* Three accelerations to work out.
  3. *did:* Write u′ in the rotating basis · *tex:* $\mathbf u'=u'_i\,\mathbf e'_i,\qquad u'_i=\frac{dx'_i}{dt}$ · *why:* u′ is built from the
     rotating observer's component rates (D14 step 6), attached to the turning axes. · *plain:* The measured velocity also
     sits on turning axes.
  4. *did:* Differentiate u′ by the product rule · *tex:* $\frac{d\mathbf u'}{dt}=\frac{du'_i}{dt}\mathbf e'_i+u'_i\frac{d\mathbf e'_i}{dt}$ · *why:* Both the
     components and the axes depend on time (product rule). · *plain:* The measured velocity changes, and its axes turn.
  5. *did:* Name the measured acceleration · *tex:* $\frac{du'_i}{dt}\mathbf e'_i=\mathbf a'$ · *why:* This is what the rotating observer
     calls acceleration: the rates of its own velocity components. · *plain:* The acceleration seen in the frame.
  6. *did:* Turn the axes: first Ω × u′ · *tex:* $u'_i\frac{d\mathbf e'_i}{dt}=u'_i\,\boldsymbol\Omega\times\mathbf e'_i=\boldsymbol\Omega\times\mathbf u'$ · *why:*
     D14's turning rule de′_i/dt = Ω × e′_i, then linearity of ×. · *plain:* The turning axes contribute one Ω × u′. ·
     *watch:* "half of the amber Coriolis arrow appears".
  7. *did:* Differentiate Ω × x′ · *tex:* $\frac{d}{dt}\big(\boldsymbol\Omega\times\mathbf x'\big)=\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'+\boldsymbol\Omega\times\frac{d\mathbf x'}{dt}$ ·
     *why:* Product rule for a cross product (primer), keeping the order of the factors. · *plain:* The rotation rate may
     change, and the position moves.
  8. *did:* Differentiate x′ with D14 · *tex:* $\frac{d\mathbf x'}{dt}=\mathbf u'+\boldsymbol\Omega\times\mathbf x'$ · *why:* D14 steps 5–7 applied to
     x′ = x′_ie′_i alone (no U): component rates plus the turning axes. · *plain:* The relative position changes by the
     measured velocity plus the rotation.
  9. *did:* Substitute: second Ω × u′ · *tex:* $\boldsymbol\Omega\times\frac{d\mathbf x'}{dt}=\boldsymbol\Omega\times\mathbf u'+\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')$ ·
     *why:* Distributivity of the cross product over a sum. · *plain:* The moving position contributes another Ω × u′. ·
     *watch:* "the amber arrow doubles".
  10. *did:* Expand the triple product · *tex:* $\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')=\boldsymbol\Omega(\boldsymbol\Omega\cdot\mathbf x')-\Omega^2\mathbf x'$ ·
      *why:* a × (b × c) = b(a·c) − c(a·b) (Ch. 2); for Ω along z it is −Ω²R e_R, pointing to the axis: centripetal. ·
      *plain:* This term pulls toward the axis.
  11. *did:* Collect the two Ω × u′ · *tex:* $\mathbf a=\frac{d\mathbf U}{dt}+\mathbf a'+2\boldsymbol\Omega\times\mathbf u'+\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'+\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')$ ·
      *why:* Add steps 5, 6, 7 and 9: one Ω × u′ from each source gives the factor 2. This is (4.43). · *plain:* Frame,
      measured, Coriolis, angular and centripetal accelerations. · *live:* "2Ωu′ = 2 × 0.5 × 3 = 3.0 m/s², centripetal
      Ω²R = 0.25 × …".
  12. *did:* Follow a fluid particle · *tex:* $\Big(\frac{D\mathbf u}{Dt}\Big)_{O123}=\Big(\frac{D'\mathbf u'}{Dt}\Big)_{O'1'2'3'}+\frac{d\mathbf U}{dt}+2\boldsymbol\Omega\times\mathbf u'+\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'+\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')$ ·
      *why:* For a fluid particle the acceleration is the material derivative: a = Du/Dt in the inertial frame, a′ =
      D′u′/Dt in the rotating one. This is (4.44). · *plain:* The same bookkeeping, for fluid particles.
- **Result.** (4.43) and (4.44) as in steps 11–12 — *in words:* the true acceleration = what the rotating observer sees +
  frame acceleration + Coriolis + angular + centripetal terms.
- **Check.** Units m/s² ✓. Ω = 0, U̇ = 0: a = a′ (Galilean) ✓. A particle at rest on a steady turntable: a = Ω × (Ω × x′),
  the centripetal acceleration ✓. **sympy check (`check_src`, every line commented):** rotation about z by θ(t)
  (Ω = θ̇e_z, Ω̇ = θ̈e_z) as a sympy matrix R(t); generic functions X_i(t), x′_i(t); build x = X + R x′ (step 1 of D14);
  differentiate twice; rotate back into primed components (Rᵀẍ); subtract the five terms of (4.43) written in primed
  components (Ẍ rotated, ẍ′, 2Ω × ẋ′, Ω̇ × x′, Ω × (Ω × x′)); `simplify` each component → 0. A second line replaces the 2 by 1
  and shows the difference is Ω × ẋ′ ≠ 0 — the factor 2 is necessary.
- **What it means.** The Coriolis term's 2 is not a convention: half comes from the turning axes, half from the moving
  position. Moved to the force side (D16) these terms are the apparent forces of every geophysical flow.
- **Traps.** Forgetting the second Ω × u′ (a factor-of-2 error). Differentiating x′ as if the axes were fixed. Sign
  language: +2Ω × u′ here is an acceleration term; the force per mass −2Ω × u′ appears in (4.45).

### D16 · Navier–Stokes in a noninertial frame: (4.44) + (4.39b) → (4.45) — ★★, 6 steps, in C09 (notebook · `rotating_frame_coriolis`)
- **Goal.** Rewrite the incompressible Navier–Stokes equation for an observer who translates with acceleration and
  rotates.
- **Start.** (4.39b) in the inertial frame: $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$.
- **Plan.** (1) Replace the acceleration by (4.44). (2) Show the pressure gradient and viscous term are unchanged. (3) Move
  the frame terms to the force side.
- **Tools.** (4.44) (D15) · vector Laplacian (gloss) · the gradient of a scalar is the same arrow in every frame at the
  coincidence instant (gloss).
- **Assumptions.** Incompressible, constant μ (the starting equation); U and Ω uniform in space (step 4).
- **Steps.**
  1. *did:* Start from inertial Navier–Stokes · *tex:* $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ · *why:* (4.39b) holds
     only in an inertial frame; the book says "(4.39)" but needs the incompressible, constant-μ form. · *plain:* Newton's
     law in the lab. · *set:* {mode: 'turntable'}.
  2. *did:* Replace the acceleration by (4.44) · *tex:* $\rho\Big[\frac{D'\mathbf u'}{Dt}+\frac{d\mathbf U}{dt}+2\boldsymbol\Omega\times\mathbf u'+\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'+\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')\Big]=\ldots$ ·
     *why:* D15: the inertial acceleration in terms of rotating-frame quantities. · *plain:* The lab acceleration, written
     with the rider's measurements.
  3. *did:* Keep the pressure gradient · *tex:* $\nabla p=\nabla'p$ · *why:* Pressure is a scalar; its gradient is one arrow, only
     described in different components in the two frames at the same instant (gloss). · *plain:* Pressure pushes the same
     way for both observers.
  4. *did:* Keep the viscous term · *tex:* $\nabla^2\mathbf u=\nabla'^2\mathbf u'$ · *why:* U + Ω × x′ is linear in x′ (U, Ω uniform in space), so its second derivatives vanish; ∇² is the same for any orientation of the axes. · *plain:* A rigid
     motion adds no friction.
  5. *did:* Move the frame terms to the right · *tex:* $\rho\frac{D'\mathbf u'}{Dt}=-\nabla'p+\rho\mathbf g-\rho\frac{d\mathbf U}{dt}-2\rho\boldsymbol\Omega\times\mathbf u'-\rho\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'-\rho\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')+\mu\nabla'^2\mathbf u'$ ·
     *why:* Subtract them from both sides; every sign flips. They are now apparent (fictitious) body forces. · *plain:* The
     rider's accelerations become the rider's forces. · *watch:* "every bar flips: accelerations become forces".
  6. *did:* Group the body forces · *tex:* $\rho\frac{D'\mathbf u'}{Dt}=-\nabla'p+\rho\Big[\mathbf g-\frac{d\mathbf U}{dt}-2\boldsymbol\Omega\times\mathbf u'-\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')\Big]+\mu\nabla'^2\mathbf u'$ ·
     *why:* All are per-unit-mass accelerations multiplied by ρ, like gravity. This is (4.45). · *plain:* In a rotating frame,
     gravity is joined by four apparent forces.
- **Result.** (4.45) as in step 6 — *in words:* the rotating observer uses Navier–Stokes with extra body forces: frame
  acceleration, Coriolis, angular acceleration, centrifugal.
- **Check.** Units ✓. Inertial frame (U constant, Ω = 0): the bracket is g alone ✓. dU/dt = g: g − dU/dt = 0,
  weightlessness (N60) ✓. Continuity is unchanged (∇·u = ∇′·u′ since ∇·(Ω × x′) = 0).
- **What it means.** The starting point of Ch. 13: on the Earth −2Ω × u′ turns winds and currents, −Ω × (Ω × x′) is
  absorbed into gravity, and the other two vanish.
- **Traps.** Forgetting that moving the terms flips their signs (Coriolis force −2Ω × u′, centrifugal −Ω × (Ω × x′)).
  Using (4.45) with space-dependent Ω. Treating the apparent forces as real interactions (they do no work — Coriolis — or
  are potential forces — centrifugal).

### D17 · The Coriolis deflection of a projectile from the pole: Ωut², angle Ωt — ★, 5 steps, in C09 (notebook · `rotating_frame_coriolis`)
- **Goal.** Show where the book's "deflection Ωut²" comes from (it asserts it), and why the angle equals the Earth's
  rotation.
- **Start.** The Coriolis force per unit mass of (4.45): $-2\boldsymbol\Omega\times\mathbf u'$.
- **Plan.** (1) Evaluate it for horizontal motion at the pole. (2) Treat it as a constant sideways acceleration for a
  short time. (3) Integrate twice and compare with the Earth's turn.
- **Tools.** Cross product and right-hand rule (ch02 P74) · constant-acceleration kinematics s = ½at² (gloss).
- **Assumptions.** At the pole (Ω vertical), horizontal launch, no friction, centrifugal effect neglected (second order);
  small Ωt, so the path direction barely changes (step 2).
- **Steps.**
  1. *did:* Evaluate the Coriolis force at the pole · *tex:* $-2\boldsymbol\Omega\times\mathbf u'=-2\Omega u\,\mathbf e_z\times\mathbf e_x=-2\Omega u\,\mathbf e_y$ ·
     *why:* At the pole Ω = Ωe_z points up; the ball moves along e_x; e_z × e_x = e_y (right-hand rule). · *plain:* A
     sideways push to the right of the motion. · *set:* {mode: 'pole', u0: 10, t: 0}.
  2. *did:* Note it is perpendicular to the motion · *tex:* $\big(-2\boldsymbol\Omega\times\mathbf u'\big)\cdot\mathbf u'=0,\qquad\lvert-2\boldsymbol\Omega\times\mathbf u'\rvert=2\Omega u$ ·
     *why:* A cross product is perpendicular to its factors: no work is done, so the speed u stays constant; for small
     Ωt the direction barely changes, so the push is nearly fixed. · *plain:* The force steers but never speeds up.
  3. *did:* Integrate the sideways motion · *tex:* $\delta=\tfrac12(2\Omega u)t^2=\Omega ut^2$ · *why:* Constant acceleration 2Ωu from rest
     sideways: s = ½at² (gloss). · *plain:* The sideways offset grows with the square of time. · *watch:* "the grey
     parabola overlays the path".
  4. *did:* Divide by the forward distance · *tex:* $\frac{\delta}{ut}=\Omega t$ · *why:* For a small angle, angle ≈ sideways/forward
     distance (small-angle approximation, ch03 P100). · *plain:* The path turns through the angle Ωt.
  5. *did:* Compare with the Earth's turn · *tex:* $\Omega t=\text{the angle the ground turns in time }t$ · *why:* Seen from space the path
     is straight (no real force); the ground turned by Ωt under it, so the rider sees the path turned by that angle. ·
     *plain:* The deflection is the rotation of the ground. · *set:* {t: 3600} · *live:* "δ = 7.292×10⁻⁵ × 10 × 3600² =
     9451 m; exact ut sin Ωt = 9342 m".
- **Result.** $\delta\approx\Omega ut^2$, angle Ωt — *in words:* at the pole a projectile drifts right by Ωut², which is just the
  ground turning under a straight path.
- **Check.** Units (1/s)(m/s)(s²) = m ✓. Ω = 0 → no deflection ✓. 10 m/s for 1 h: 9.45 km (small-angle) vs 9.34 km exact —
  1.2 % apart at Ωt = 15° ✓.
- **What it means.** Tiny per second, decisive over hours: why weather systems and ocean currents turn (Ch. 13). Away from
  the pole only the vertical component Ω sin φ steers horizontal motion.
- **Traps.** Using 2Ωut² (forgetting the ½). Deflecting to the left in the NH (sign of −2Ω × u′). Applying the formula for
  large Ωt (the path curves; use the exact offset).

### D18 · The centrifugal term, its potential and effective gravity — ★, 5 steps, in C09 (notebook · `rotating_frame_coriolis`)
- **Goal.** Show that the centrifugal force per unit mass points away from the axis with size Ω²R and is the gradient of a
  potential, so it can be folded into gravity (the book leaves it to Exercise 4.43).
- **Start.** The centrifugal term of (4.45): $-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')$.
- **Plan.** (1) Expand the triple product. (2) Put Ω along z and recognise the distance from the axis. (3) Find the
  potential and add it to gravity's.
- **Tools.** Vector triple product (gloss, Ch. 2 §2.7) · conservative force and its potential (primer in C04) ·
  cylindrical gradient of a function of R (gloss).
- **Assumptions.** Steady rotation about a fixed axis (z).
- **Steps.**
  1. *did:* Expand the triple product · *tex:* $-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')=-\boldsymbol\Omega(\boldsymbol\Omega\cdot\mathbf x')+\Omega^2\mathbf x'$ ·
     *why:* a × (b × c) = b(a·c) − c(a·b) with a = b = Ω, c = x′; then the overall minus sign. · *plain:* Two simpler
     vectors. · *set:* {mode: 'effg', lat: 45}.
  2. *did:* Put Ω along z · *tex:* $\boldsymbol\Omega=\Omega\,\mathbf e_z:\quad-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')=\Omega^2\big(\mathbf x'-z'\mathbf e_z\big)$ ·
     *why:* Ω·x′ = Ωz′, so Ω(Ω·x′) = Ω²z′e_z. · *plain:* x′ with its axial part removed, times Ω².
  3. *did:* Recognise the distance from the axis · *tex:* $\mathbf x'-z'\mathbf e_z=R\,\mathbf e_R\ \Rightarrow\ -\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')=\Omega^2R\,\mathbf e_R$ ·
     *why:* Removing the axial component leaves the cylindrical radial vector (Ch. 3 §3.1). · *plain:* Outward, growing with
     the distance from the axis.
  4. *did:* Find its potential · *tex:* $-\nabla\Big(-\tfrac12\Omega^2R^2\Big)=\Omega^2R\,\mathbf e_R$ · *why:* For a function of R alone the
     cylindrical gradient is (df/dR)e_R (gloss); d(−½Ω²R²)/dR = −Ω²R. A conservative force (primer). · *plain:* The
     centrifugal force has a potential energy −½Ω²R² per unit mass.
  5. *did:* Add it to gravity · *tex:* $\mathbf g_e=\mathbf g_n+\Omega^2R\,\mathbf e_R=-\nabla\Big(\Phi_n-\tfrac12\Omega^2R^2\Big)$ · *why:* Sum of two
     conservative forces is conservative; its potential is the geopotential whose level surfaces are sea level (Fig. 4.9).
     · *plain:* A plumb line and the sea surface feel this effective gravity. · *live:* "Ω²a cos φ = 0.0339 × 0.707 = 0.0240
     m/s²".
- **Result.** $-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')=\Omega^2R\,\mathbf e_R=-\nabla\big(-\tfrac12\Omega^2R^2\big)$, and
  $\mathbf g_e=\mathbf g_n+\Omega^2R\,\mathbf e_R$ — *in words:* the centrifugal force is outward, Ω²R, and potential, so it merges with
  gravity.
- **Check.** Units m/s² ✓. On the axis R = 0: no centrifugal force ✓. Equator: Ω²a = 0.0339 m/s², 0.35 % of g ✓.
- **What it means.** The Earth's equatorial bulge and the "g" of weather models already include it (the geopotential of
  Ch. 13); only Coriolis remains as a separate apparent force.
- **Traps.** R is the distance from the rotation axis, not from the Earth's centre. The sign of the potential (−½Ω²R², so
  the force points outward). Calling Ω × (Ω × x′) "centrifugal" (that is the centripetal acceleration term of (4.43)).

### D19 · The total-energy equation: (4.46) → (4.47) → (4.49)–(4.52) → (4.53) — ★★, 8 steps, in C10 (notebook)
- **Goal.** Turn the first law of thermodynamics for a material volume into a field equation for the total energy
  e + ½u² at a point.
- **Start.** (4.46): $\frac{d}{dt}\int_{V}\rho\big(e+\tfrac12\lvert\mathbf u\rvert^2\big)dV=\int_{V}\rho\mathbf g\cdot\mathbf u\,dV+\int_{A}\mathbf f\cdot\mathbf u\,dA-\int_{A}\mathbf q\cdot\mathbf n\,dA$ —
  *in words:* the internal plus kinetic energy of a material volume rises by the work of gravity and of the surface forces,
  and falls by the heat conducted out.
- **Plan.** (1) Expand the material rate. (2) Turn the three surface integrals into volume integrals with Gauss. (3)
  Collect and localise.
- **Tools.** RTT (R01) · power of a force and heat flux through a surface (primer in C10) · Gauss (R02) · traction
  $f_j=n_i\tau_{ij}$ (R03) · localisation lemma (N07).
- **Assumptions.** Heat transfer by a flux q (conduction; radiation included if it is a flux); fields smooth; every material
  volume (step 8).
- **Steps.**
  1. *did:* Start from the first law · *tex:* $\frac{d}{dt}\int_{V}\rho\big(e+\tfrac12u_ju_j\big)dV=\int_{V}\rho g_iu_i\,dV+\int_{A}f_ju_j\,dA-\int_{A}q_in_i\,dA$ ·
     *why:* (4.46) in index form: work rates are force · velocity, heat leaving counts negative (primer). · *plain:* Energy
     in = work done + heat in.
  2. *did:* Expand with the transport theorem · *tex:* $\int_{V}\frac{\partial}{\partial t}\Big(\rho\big[e+\tfrac12u_j^2\big]\Big)dV+\int_{A}\rho\big[e+\tfrac12u_j^2\big]u_in_i\,dA=\ldots$ ·
     *why:* RTT (3.35) with F = ρ(e + ½u_j²), b = u. This is (4.47); u_j² is summed over j. · *plain:* Local change plus
     energy carried by the moving skin.
  3. *did:* Gauss on the energy flux · *tex:* $\int_{A}\rho\big[e+\tfrac12u_j^2\big]u_in_i\,dA=\int_{V}\frac{\partial}{\partial x_i}\Big(\rho\big[e+\tfrac12u_j^2\big]u_i\Big)dV$ ·
     *why:* Gauss (2.30) on the vector field ρEu_i. This is (4.49). · *plain:* Energy carried out = its flux divergence.
  4. *did:* Write the stress power with τ · *tex:* $\int_{A}f_ju_j\,dA=\int_{A}n_i\tau_{ij}u_j\,dA$ · *why:* Traction f_j = n_iτ_ij (R03, first
     index = face). · *plain:* The surface forces' power per area is n·τ·u.
  5. *did:* Gauss on the stress power · *tex:* $\int_{A}n_i\tau_{ij}u_j\,dA=\int_{V}\frac{\partial}{\partial x_i}\big(\tau_{ij}u_j\big)dV$ · *why:* For the
     vector τ_iju_j (index i), Gauss contracts i with n_i. This is (4.50). · *plain:* The surface work becomes a divergence.
  6. *did:* Gauss on the heat flux · *tex:* $\int_{A}q_in_i\,dA=\int_{V}\frac{\partial q_i}{\partial x_i}dV$ · *why:* Gauss on q. This is (4.51)
     — the book prints dA inside the volume integrals, a slip: it must be dV. · *plain:* Heat leaving = divergence of the
     heat flux.
  7. *did:* Collect under one integral · *tex:* $\int_{V}\Big\{\frac{\partial(\rho E)}{\partial t}+\frac{\partial(\rho Eu_i)}{\partial x_i}-\rho g_iu_i-\frac{\partial(\tau_{ij}u_j)}{\partial x_i}+\frac{\partial q_i}{\partial x_i}\Big\}dV=0$ ·
     *why:* All terms are integrals over V; E = e + ½u_j² for short. This is (4.52). · *plain:* One integrand, zero integral.
  8. *did:* Localise · *tex:* $\frac{\partial(\rho E)}{\partial t}+\frac{\partial(\rho Eu_i)}{\partial x_i}=\rho g_iu_i+\frac{\partial(\tau_{ij}u_j)}{\partial x_i}-\frac{\partial q_i}{\partial x_i}$ ·
     *why:* True for every material volume with a continuous integrand: the localisation lemma (D02). This is (4.53). ·
     *plain:* Total energy is conserved at every point.
- **Result.** $\frac{\partial}{\partial t}\big(\rho[e+\tfrac12u_j^2]\big)+\frac{\partial}{\partial x_i}\big(\rho[e+\tfrac12u_j^2]u_i\big)=\rho g_iu_i+\frac{\partial}{\partial x_i}(\tau_{ij}u_j)-\frac{\partial q_i}{\partial x_i}$ (4.53).
- **Check.** Units W/m³ ✓. Fluid at rest, τ = −pδ: ∂(ρe)/∂t = −∇·q — pure conduction ✓. Steady Couette: the moving wall's
  work τU arrives through ∂(τ_iju_j)/∂x_i and leaves as heat ✓ (C10 worked example).
- **What it means.** The conservative energy equation used by compressible codes (Ch. 10, 15). D20 rewrites it for a
  particle and D21–D22 split it into mechanical and thermal parts.
- **Traps.** The dA/dV slip of (4.51). Forgetting that u_j² is summed. Sign of the heat term (outflow positive q·n is a
  loss).

### D20 · Total energy following a particle, (4.53) → (4.55) — ★★, 7 steps, in C10 (notebook)
- **Goal.** Rewrite the total-energy equation with D/Dt and split the stress work into its physical parts (the book
  leaves it to Exercise 4.45).
- **Start.** (4.53) with E = e + ½u_j²: $\frac{\partial(\rho E)}{\partial t}+\frac{\partial(\rho Eu_i)}{\partial x_i}=\rho g_iu_i+\frac{\partial(\tau_{ij}u_j)}{\partial x_i}-\frac{\partial q_i}{\partial x_i}$.
- **Plan.** (1) Expand the left side and remove E × continuity (the move of (4.23)). (2) Expand the stress work with the
  product rule and τ = −pδ + σ.
- **Tools.** Product rule for a divergence (primer in C02) · continuity (4.7) · (4.27) (C07) · δ substitution (gloss).
- **Assumptions.** Continuity holds (step 4).
- **Steps.**
  1. *did:* Name the specific total energy · *tex:* $E\equiv e+\tfrac12u_ju_j$ · *why:* A short name for internal plus kinetic energy
     per unit mass keeps the algebra readable. · *plain:* Energy per kilogram.
  2. *did:* Expand the left side · *tex:* $\frac{\partial(\rho E)}{\partial t}+\frac{\partial(\rho Eu_i)}{\partial x_i}=\rho\frac{\partial E}{\partial t}+E\frac{\partial\rho}{\partial t}+\rho u_i\frac{\partial E}{\partial x_i}+E\frac{\partial(\rho u_i)}{\partial x_i}$ ·
     *why:* Product rule twice (the second with the factors E and ρu_i). · *plain:* Four pieces.
  3. *did:* Group E times continuity · *tex:* $=\rho\Big(\frac{\partial E}{\partial t}+u_i\frac{\partial E}{\partial x_i}\Big)+E\Big[\frac{\partial\rho}{\partial t}+\frac{\partial(\rho u_i)}{\partial x_i}\Big]$ ·
     *why:* Collect the terms with ρ and with E — the move of (4.23) with u_j replaced by E. · *plain:* A material
     derivative and the continuity bracket.
  4. *did:* Drop the bracket, recognise D/Dt · *tex:* $\frac{\partial(\rho E)}{\partial t}+\frac{\partial(\rho Eu_i)}{\partial x_i}=\rho\frac{DE}{Dt}$ · *why:*
     The bracket is (4.7) = 0; the first group is (3.5) with F = E. · *plain:* The flux form equals ρ times the particle's
     rate.
  5. *did:* Expand the stress work · *tex:* $\frac{\partial(\tau_{ij}u_j)}{\partial x_i}=\tau_{ij}\frac{\partial u_j}{\partial x_i}+u_j\frac{\partial\tau_{ij}}{\partial x_i}$ · *why:* Product
     rule. The first term is work that deforms the particle, the second the work of the net force on it. · *plain:*
     Deformation work plus force work.
  6. *did:* Insert τ = −pδ + σ · *tex:* $\frac{\partial(\tau_{ij}u_j)}{\partial x_i}=\Big(-p\frac{\partial u_j}{\partial x_j}+\sigma_{ij}\frac{\partial u_j}{\partial x_i}\Big)+\Big(-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}\Big)$ ·
     *why:* (4.27); δ_ij∂u_j/∂x_i = ∂u_j/∂x_j and u_jδ_ij∂p/∂x_i = u_j∂p/∂x_j (δ substitution). This is (4.54). · *plain:*
     Pressure and viscous parts of each kind of work.
  7. *did:* Assemble · *tex:* $\rho\frac{D}{Dt}\Big(e+\tfrac12u_j^2\Big)=\rho g_iu_i+\Big(-p\frac{\partial u_j}{\partial x_j}+\sigma_{ij}\frac{\partial u_j}{\partial x_i}\Big)+\Big(-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}\Big)-\frac{\partial q_i}{\partial x_i}$ ·
     *why:* Steps 4 and 6 in (4.53). This is (4.55). · *plain:* A particle's total energy changes by gravity's work, the two
     kinds of stress work, and conduction.
- **Result.** (4.55) as in step 7.
- **Check.** Units W/m³ ✓. Summing the four stress-work pieces gives back ∂(τ_iju_j)/∂x_i (`ch04.stress_work_split` on a
  random field, C10) ✓. Fluid at rest: ρDe/Dt = −∇·q ✓.
- **What it means.** The split (4.54) is the key to the next two derivations: the force work belongs to the kinetic
  energy, the deformation work to the internal energy.
- **Traps.** Subtracting E × ∇·u instead of E × continuity. Losing the δ in −p∂u_j/∂x_j. Mixing up the index order in
  σ_ij∂u_j/∂x_i (it is S_ij that survives, D22).

### D21 · The mechanical-energy equation, (4.56) — ★★, 6 steps, in C10 (notebook · `viscous_dissipation_heating`)
- **Goal.** An equation for the kinetic energy alone, to separate it from the internal energy (the book says "after some
  manipulation" and leaves it to Exercise 4.46).
- **Start.** Cauchy's equation (4.24): $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$.
- **Plan.** (1) Dot Cauchy's equation with the velocity. (2) Use the chain rule on u_jDu_j/Dt. (3) Split the stress.
  (4) Read the stress term as force work.
- **Tools.** Chain rule for the kinetic energy (primer in C10) · (4.27) (C07) · (4.54) (D20).
- **Assumptions.** None beyond Cauchy's equation (any continuum; the Newtonian law is not needed yet).
- **Steps.**
  1. *did:* Start from Cauchy's equation · *tex:* $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ · *why:* (4.24). The book
     multiplies the flux form (4.22); starting from (4.24) saves subtracting u_j × continuity (D07 steps 6–9 already did).
     · *plain:* Newton's law for a particle. · *set:* {flow: 'couette', U: 1}.
  2. *did:* Dot with the velocity · *tex:* $\rho u_j\frac{Du_j}{Dt}=\rho g_ju_j+u_j\frac{\partial\tau_{ij}}{\partial x_i}$ · *why:* Multiply component j by u_j
     and sum over j: force · velocity is power (primer in C10), so a force balance becomes a power balance. · *plain:*
     The rate at which forces do work on the particle.
  3. *did:* Recognise the kinetic energy's rate · *tex:* $\rho u_j\frac{Du_j}{Dt}=\rho\frac{D}{Dt}\Big(\tfrac12u_ju_j\Big)$ · *why:* Chain rule for
     the kinetic energy (primer): D(½u_ju_j)/Dt = u_jDu_j/Dt. · *plain:* The left side is the particle's kinetic energy
     changing. · *watch:* "the kinetic-energy bar lights (zero in steady Couette)".
  4. *did:* Split the stress · *tex:* $u_j\frac{\partial\tau_{ij}}{\partial x_i}=-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}$ · *why:* (4.27),
     τ_ij = −pδ_ij + σ_ij, and δ substitution. · *plain:* Work of the net pressure force and of the net viscous force.
  5. *did:* Read it as force work · *tex:* $u_j\frac{\partial\tau_{ij}}{\partial x_i}=\frac{\partial(\tau_{ij}u_j)}{\partial x_i}-\tau_{ij}\frac{\partial u_j}{\partial x_i}$ · *why:*
     Product rule, i.e. (4.54) read backwards: the kinetic energy receives the total stress work minus the deformation work.
     · *plain:* Kinetic energy gets only the work of the net force. · *watch:* "the stress-work bar splits in two".
  6. *did:* Assemble · *tex:* $\rho\frac{D}{Dt}\Big(\tfrac12u_j^2\Big)=\rho g_ju_j-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}$ · *why:*
     Steps 3 and 4 in step 2. This is (4.56). · *plain:* Kinetic energy changes by the work of gravity, of the pressure force
     and of the viscous force.
- **Result.** $\rho\frac{D}{Dt}\big(\tfrac12u_j^2\big)=\rho g_ju_j-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}$ (4.56) — *in words:* the
  mechanical-energy equation.
- **Check.** Units W/m³ ✓. Steady Couette: no acceleration, no pressure gradient, u_j∂σ_ij/∂x_i = u μu″ = 0 ✓ (the wall's work
  all goes into deformation, D22). Hydrostatics: 0 = 0 ✓. `ch04.kinetic_energy_budget` residual ≈ 0 on a test field (C10).
- **What it means.** Viscosity enters kinetic energy only through the net viscous force, which can speed particles up or
  slow them down; its deformation work is not here — it goes to heat. The turbulent kinetic-energy budget of Ch. 12 is this
  equation for fluctuations.
- **Traps.** Multiplying the flux form (4.22) by u_j without subtracting u_j × continuity (the book's wording). Forgetting
  the sum over j. Treating u_j∂σ_ij/∂x_i as dissipation (it is not sign-definite).

### D22 · The internal-energy equation: (4.55) − (4.56) → (4.57) — ★★, 8 steps, in C10 (notebook · `viscous_dissipation_heating`)
- **Goal.** Subtract the mechanical part from the total energy and find what heats a fluid particle.
- **Start.** (4.55) and (4.56) — *in words:* the total-energy and kinetic-energy equations following a particle.
- **Plan.** (1) Subtract; watch the force-work terms cancel. (2) Divide by ρ and turn ∇·u into −Dv/Dt with continuity. (3)
  Show only S survives in the viscous work.
- **Tools.** (4.55) (D20), (4.56) (D21) · continuity with D/Dt (4.8) (C02) · quotient rule (gloss) · σ:R = 0 for symmetric
  σ (R07) · G = S + ½R (3.11) (Ch. 3 §3.4).
- **Assumptions.** σ symmetric (step 7, from (4.25)); ρ > 0.
- **Steps.**
  1. *did:* Subtract (4.56) from (4.55) · *tex:* $\rho\frac{De}{Dt}=-p\frac{\partial u_j}{\partial x_j}+\sigma_{ij}\frac{\partial u_j}{\partial x_i}-\frac{\partial q_i}{\partial x_i}$ · *why:*
     D(e + ½u²)/Dt − D(½u²)/Dt = De/Dt; gravity's work and the force work appear in both and cancel exactly; only
     deformation work and conduction remain. · *plain:* Internal energy changes by deformation work and heat. · *set:*
     {flow: 'couette', U: 1}.
  2. *did:* Divide by ρ · *tex:* $\frac{De}{Dt}=-\frac p\rho\frac{\partial u_j}{\partial x_j}+\frac1\rho\sigma_{ij}\frac{\partial u_j}{\partial x_i}-\frac1\rho\frac{\partial q_i}{\partial x_i}$ ·
     *why:* ρ > 0; per unit mass, like the first law of Ch. 1. · *plain:* The same, per kilogram.
  3. *did:* Replace ∇·u with continuity · *tex:* $-\frac p\rho\nabla\cdot\mathbf u=\frac p\rho\cdot\frac1\rho\frac{D\rho}{Dt}$ · *why:* (4.8),
     $\frac1\rho\frac{D\rho}{Dt}+\nabla\cdot\mathbf u=0$, so ∇·u = −(1/ρ)Dρ/Dt. · *plain:* Expansion rate = minus the fractional
     density rate.
  4. *did:* Use the quotient rule · *tex:* $\frac1{\rho^2}\frac{D\rho}{Dt}=-\frac{D}{Dt}\Big(\frac1\rho\Big)=-\frac{Dv}{Dt}$ · *why:* Quotient rule (gloss):
     D(1/ρ)/Dt = −ρ⁻²Dρ/Dt; v = 1/ρ is the specific volume. · *plain:* A density rise is a specific-volume fall.
  5. *did:* Get the compression work · *tex:* $-\frac p\rho\frac{\partial u_j}{\partial x_j}=-p\frac{Dv}{Dt}$ · *why:* Combine steps 3 and 4. It is
     the −p dv of Ch. 1's first law, per unit time. · *plain:* Compressing a particle heats it; expanding cools it.
  6. *did:* Split the velocity gradient · *tex:* $\frac{\partial u_j}{\partial x_i}=S_{ji}+\tfrac12R_{ji}$ · *why:* G = S + ½R (3.11) with the
     book's R = G − Gᵀ; the book prints S_ji + R_ji (the ½ is missing, harmless as the next step shows). · *plain:* A
     deformation part and a rotation part.
  7. *did:* Drop the rotation part · *tex:* $\sigma_{ij}\frac{\partial u_j}{\partial x_i}=\sigma_{ij}S_{ji}=\sigma_{ij}S_{ij}$ · *why:* R07: σ symmetric and
     R antisymmetric ⇒ σ_ijR_ji = 0; S symmetric ⇒ S_ji = S_ij. · *plain:* Spinning a particle does no viscous work; only
     deforming it does. · *watch:* "only S contributes".
  8. *did:* Assemble · *tex:* $\frac{De}{Dt}=-p\frac{Dv}{Dt}+\frac1\rho\sigma_{ij}S_{ij}-\frac1\rho\frac{\partial q_i}{\partial x_i}$ · *why:* Steps 2, 5 and 7.
     This is (4.57). · *plain:* Internal energy rises by compression, by viscous deformation work (ε) and by heat
     conducted in. · *live:* "Couette mid-gap, steady: 0 = 0 + 1.0 − 1.0 W/kg".
- **Result.** $\frac{De}{Dt}=-p\frac{Dv}{Dt}+\frac1\rho\sigma_{ij}S_{ij}-\frac1\rho\frac{\partial q_i}{\partial x_i}$ (4.57) — *in words:* the internal-energy
  equation: the first law per unit mass and time for a particle.
- **Check.** Units W/kg ✓. Fluid at rest: De/Dt = −(1/ρ)∇·q (conduction) ✓. Steady Couette: ε = 1.0 W/kg is conducted away
  at the same rate ✓. The pressure-work term matches Ch. 1's −p dv ✓.
- **What it means.** Energy moves between the two budgets through two doors: −p∇·u (reversible, either sign) and ε
  (irreversible, D23). Ch. 13's heat equation and C13's Boussinesq form start here.
- **Traps.** Losing the exact cancellation of the force work (a sign slip leaves a spurious term). Writing −p∇·u as
  −pDv/Dt without continuity. Thinking the rotation part contributes (it cannot, σ is symmetric).

### D23 · Dissipation as a sum of squares: ε ≥ 0, (4.58) — ★★, 8 steps, in C10 (notebook · `viscous_dissipation_heating`)
- **Goal.** Show the viscous heating ε = σ_ijS_ij/ρ can never be negative — viscosity only turns motion into heat. The
  book leaves the square completion to Exercise 4.47.
- **Start.** $\varepsilon\equiv\frac1\rho\sigma_{ij}S_{ij}$ with the Newtonian σ of (4.59).
- **Plan.** (1) Insert σ and multiply out. (2) Regroup by μ and μ_v. (3) Complete the square with the deviatoric part.
- **Tools.** (4.59) (C07/C10) · completing the square for tensors (primer in C10) · deviatoric part (primer in C07) ·
  δ_ijδ_ij = 3.
- **Assumptions.** Newtonian fluid.
- **Steps.**
  1. *did:* Define the dissipation · *tex:* $\varepsilon\equiv\frac1\rho\sigma_{ij}S_{ij}$ · *why:* The viscous work term of (4.57) per unit mass;
     we want to know its sign. · *plain:* Rate at which viscous stress works on deformation. · *set:* {flow: 'couette',
     U: 1}.
  2. *did:* Insert the Newtonian viscous stress · *tex:* $\sigma_{ij}=2\mu S_{ij}+\Big(\mu_v-\tfrac23\mu\Big)S_{mm}\delta_{ij}$ · *why:* (4.59) written
     with S (2μS_ij = μ(∂u_i/∂x_j + ∂u_j/∂x_i)). · *plain:* The Newtonian law of C07.
  3. *did:* Multiply out · *tex:* $\rho\varepsilon=2\mu S_{ij}S_{ij}+\Big(\mu_v-\tfrac23\mu\Big)S_{mm}^2$ · *why:* δ_ijS_ij = S_ii = S_mm (δ
     substitution). · *plain:* A sum of squares and a squared trace — but with a minus sign inside.
  4. *did:* Regroup by viscosity · *tex:* $\rho\varepsilon=2\mu\Big(S_{ij}S_{ij}-\tfrac13S_{mm}^2\Big)+\mu_vS_{mm}^2$ · *why:* Take −⅔μS_mm² into the μ
     bracket as 2μ × (−⅓S_mm²). · *plain:* Is the μ bracket positive?
  5. *did:* Square the deviatoric part · *tex:* $\Big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\Big)^2=S_{ij}S_{ij}-\tfrac23S_{mm}S_{ii}+\tfrac19S_{mm}^2\delta_{ij}\delta_{ij}$ ·
     *why:* Expand the product (A − B)(A − B) summed over i, j (primer). · *plain:* Three terms.
  6. *did:* Use δ_ijδ_ij = 3 and S_ii = S_mm · *tex:* $\Big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\Big)^2=S_{ij}S_{ij}-\tfrac13S_{mm}^2$ · *why:* δ_ijδ_ij = δ_ii = 3
     (not 1): −⅔S_mm² + ⅓S_mm² = −⅓S_mm². · *plain:* The μ bracket is exactly this square.
  7. *did:* Substitute · *tex:* $\varepsilon=2\nu\Big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\Big)^2+\frac{\mu_v}{\rho}S_{mm}^2$ · *why:* Step 6 in step 4, divided by
     ρ (ν = μ/ρ). This is (4.58). · *plain:* Dissipation = shape-change squared + volume-change squared.
  8. *did:* Read off the sign · *tex:* $\varepsilon\ge0\quad\Longleftrightarrow\quad\mu\ge0,\ \mu_v\ge0$ · *why:* Sums of squares are ≥ 0; with
     non-negative coefficients the total is. Incompressible flow: ε = 2νS_ijS_ij. · *plain:* Viscosity can only heat. ·
     *live:* "ε = 2ν(S₁₂² + S₂₁²) = 2 × 10⁻⁶ × 2 × 500² = 1.0 W/kg".
- **Result.** $\varepsilon=2\nu\big(S_{ij}-\tfrac13\frac{\partial u_m}{\partial x_m}\delta_{ij}\big)^2+\frac{\mu_v}{\rho}\big(\frac{\partial u_m}{\partial x_m}\big)^2\ge0$ (4.58) — *in
  words:* viscous dissipation is a sum of squares, never negative.
- **Check.** Units W/kg ✓. Rigid motion: S = 0, ε = 0 ✓. Couette U = 1 m/s, h = 1 mm, water: 1.0 W/kg ✓. 1000 random
  gradients: contraction form = sum of squares to 1e-12, all ≥ 0 (C10 code) ✓.
- **What it means.** The one-way valve from motion to heat; the second law demands μ, μ_v ≥ 0 (4.63). In turbulence ε is
  the rate at which the whole eddy cascade ends (Ch. 12).
- **Traps.** δ_ijδ_ij = 3, not 1. Squaring S instead of its deviatoric part (then the μ_v term would be wrong). Thinking ε
  is the viscous *force* work (that is u_j∂σ_ij/∂x_i, D21, any sign).

### D24 · The Bernoulli function: (4.66) + (4.67) + Lamb identity (4.68) → (4.69) — ★★, 10 steps, in C11 (notebook · `which_bernoulli`)
- **Goal.** Rewrite Euler's equation so that every gradient is collected into one scalar B, leaving only ∂u/∂t and the Lamb
  vector u × ω — the book leaves the Lamb identity to Exercise 4.50.
- **Start.** (4.66): $\frac{\partial u_j}{\partial t}+u_i\frac{\partial u_j}{\partial x_i}=-\frac1\rho\frac{\partial p}{\partial x_j}-\frac{\partial\Phi}{\partial x_j}$ — *in words:* Euler's equation
  with gravity from its potential Φ = gz.
- **Plan.** (1) Write the pressure term as a gradient (barotropic). (2) Prove the Lamb identity with the ε–δ identity. (3)
  Substitute and collect the gradients.
- **Tools.** FTC with a variable upper limit (ch02 P84) and the chain rule · barotropic pressure function (N87) ·
  alternating tensor and ε–δ (2.19) (ch02 P72, C08's primer) · conservative force (primer in C04).
- **Assumptions.** Inviscid (the starting equation); barotropic ρ = ρ(p) (step 2); conservative body force (Φ).
- **Steps.**
  1. *did:* Start from Euler with a potential · *tex:* $\frac{\partial u_j}{\partial t}+u_i\frac{\partial u_j}{\partial x_i}=-\frac1\rho\frac{\partial p}{\partial x_j}-\frac{\partial\Phi}{\partial x_j}$ ·
     *why:* (4.66): (4.41) divided by ρ with g = −∇Φ (4.18), Φ = gz. · *plain:* Acceleration = pressure push + gravity. ·
     *set:* {sc: 'vortex'}.
  2. *did:* Write the pressure term as a gradient · *tex:* $\frac1\rho\frac{\partial p}{\partial x_j}=\frac{\partial}{\partial x_j}\int_{p_o}^{p}\frac{dp'}{\rho(p')}$ · *why:*
     Barotropic: ρ = ρ(p). Let P(p) = ∫dp′/ρ(p′); chain rule and FTC: ∂P/∂x_j = (1/ρ(p))∂p/∂x_j. This is (4.67); p′ is only
     the integration variable. · *plain:* In a barotropic fluid the pressure force per mass is a gradient.
  3. *did:* Write u × ω with ε · *tex:* $(\mathbf u\times\boldsymbol\omega)_j=\varepsilon_{jkl}\,u_k\,\omega_l$ · *why:* Definition of the cross product
     with the alternating tensor (Ch. 2). · *plain:* The Lamb vector in index form.
  4. *did:* Insert the vorticity · *tex:* $(\mathbf u\times\boldsymbol\omega)_j=\varepsilon_{jkl}\varepsilon_{lmn}\,u_k\frac{\partial u_n}{\partial x_m}$ · *why:* ω_l =
     ε_lmn∂u_n/∂x_m (ω = ∇×u). · *plain:* Two ε's in a row.
  5. *did:* Apply the ε–δ identity · *tex:* $\varepsilon_{jkl}\varepsilon_{lmn}=\varepsilon_{ljk}\varepsilon_{lmn}=\delta_{jm}\delta_{kn}-\delta_{jn}\delta_{km}$ · *why:* Cyclic
     shift ε_jkl = ε_ljk, then (2.19) with the shared first index. · *plain:* Two ε's become δ's.
  6. *did:* Substitute the δ's · *tex:* $(\mathbf u\times\boldsymbol\omega)_j=u_k\frac{\partial u_k}{\partial x_j}-u_k\frac{\partial u_j}{\partial x_k}$ · *why:*
     δ_jmδ_kn: m = j, n = k; δ_jnδ_km: n = j, m = k (δ substitution). · *plain:* The Lamb vector is a difference of two
     products of u and its gradient. · *watch:* "u × ω arrows appear inside the core".
  7. *did:* Recognise a kinetic-energy gradient · *tex:* $u_i\frac{\partial u_j}{\partial x_i}=-(\mathbf u\times\boldsymbol\omega)_j+\frac{\partial}{\partial x_j}\Big(\tfrac12u_iu_i\Big)$ ·
     *why:* u_k∂u_k/∂x_j = ∂(½u_ku_k)/∂x_j (chain rule); rearrange step 6 and rename k as i. This is the Lamb identity
     (4.68). · *plain:* Advective acceleration = gradient of kinetic energy minus the Lamb vector.
  8. *did:* Substitute into Euler · *tex:* $\frac{\partial u_j}{\partial t}-(\mathbf u\times\boldsymbol\omega)_j+\frac{\partial}{\partial x_j}\Big(\tfrac12u_i^2\Big)=-\frac{\partial}{\partial x_j}\int\frac{dp'}{\rho}-\frac{\partial(gz)}{\partial x_j}$ ·
     *why:* Steps 2 and 7 in step 1, Φ = gz. · *plain:* Every term but two is a gradient.
  9. *did:* Collect the gradients · *tex:* $\frac{\partial u_j}{\partial t}+\frac{\partial}{\partial x_j}\Big[\tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz\Big]=(\mathbf u\times\boldsymbol\omega)_j$ ·
     *why:* Move the gradients left and u × ω right (linearity of the gradient). This is (4.69). · *plain:* Acceleration plus
     the gradient of one bracket equals the Lamb vector.
  10. *did:* Name the Bernoulli function · *tex:* $B\equiv\tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz:\qquad\frac{\partial\mathbf u}{\partial t}+\nabla B=\mathbf u\times\boldsymbol\omega$ ·
      *why:* One scalar collects kinetic energy, the pressure function and potential energy per unit mass. · *plain:* The
      Bernoulli function B. · *live:* "inside the core at r = 0.5 m: ∇B = u × ω = 1.0 m/s² outward".
- **Result.** $\frac{\partial u_j}{\partial t}+\frac{\partial}{\partial x_j}\big[\tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz\big]=(\mathbf u\times\boldsymbol\omega)_j$ (4.69) — *in
  words:* Euler's equation as "∂u/∂t + ∇B = u × ω".
- **Check.** Units m/s² ✓. Rankine core (r = 0.5 m, Γ = 2π, σ = 1): dB/dr = 2r = 1.0 and u_θω = 0.5 × 2 = 1.0 m/s² ✓. Lamb
  identity residual ≈ 0 at a point (C11 code); `ch04.lamb_identity_sym` returns zero vectors ✓.
- **What it means.** The whole Bernoulli family follows (D25, D26); Ch. 5 starts the vorticity equation from (4.68). Needs
  barotropy: in a baroclinic fluid (the real atmosphere) no single B exists — the root of thermal wind (Ch. 13).
- **Traps.** The sign of u × ω. Reading p′ as a perturbation (it is a dummy variable). Using ∫dp/ρ for a fluid whose density
  depends on temperature too.

### D25 · B is constant along streamlines and vortex lines, (4.70) → (4.71), (4.72) — ★, 5 steps, in C11 (notebook · `which_bernoulli`)
- **Goal.** Read off from (4.69) exactly along what B is constant in steady flow — and when it is the same everywhere.
- **Start.** (4.69): $\frac{\partial\mathbf u}{\partial t}+\nabla B=\mathbf u\times\boldsymbol\omega$.
- **Plan.** (1) Steady flow. (2) Dot with u, then with ω. (3) Irrotational flow.
- **Tools.** Triple product with a repeated vector (gloss) · directional derivative (ch02 P75) · zero gradient ⇒ constant
  (gloss).
- **Assumptions.** Steady (step 1); irrotational in a connected region (step 5).
- **Steps.**
  1. *did:* Make the flow steady · *tex:* $\nabla B=\mathbf u\times\boldsymbol\omega$ · *why:* ∂u/∂t = 0. This is (4.70): ∇B is normal to the surfaces
     B = const, u × ω is perpendicular to both u and ω. · *plain:* B's gradient is the Lamb vector. · *set:* {sc: 'vortex'}.
  2. *did:* Dot with the velocity · *tex:* $\mathbf u\cdot\nabla B=\mathbf u\cdot(\mathbf u\times\boldsymbol\omega)=0$ · *why:* A triple product with a repeated
     vector vanishes (gloss): u × ω is perpendicular to u. · *plain:* B does not change in the direction of the flow.
  3. *did:* Read it along a streamline · *tex:* $\frac{dB}{ds}=\frac{\mathbf u}{\lvert\mathbf u\rvert}\cdot\nabla B=0$ · *why:* The rate of change along
     the unit tangent of a streamline is the directional derivative (ch02 P75). · *plain:* B is constant along each
     streamline. · *live:* "B along the circle r = 0.5 m: −0.750 everywhere".
  4. *did:* Dot with the vorticity · *tex:* $\boldsymbol\omega\cdot\nabla B=0\ \Rightarrow\ \tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=\text{const along streamlines and vortex lines}$ ·
     *why:* The same triple-product argument with ω. This is (4.71): each Lamb surface B = const contains both families of
     lines. · *plain:* B is constant along vortex lines too — but may differ from line to line.
  5. *did:* Remove the vorticity · *tex:* $\boldsymbol\omega=0\ \Rightarrow\ \nabla B=0\ \Rightarrow\ B=\text{const everywhere}$ · *why:* A zero gradient in a
     connected region means B is the same at every point (gloss). This is (4.72). · *plain:* In irrotational flow one
     Bernoulli constant holds for the whole flow. · *set:* {sc: 'cylinder'} · *watch:* "the rose curve goes flat too".
- **Result.** (4.71) and $\tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=$ constant everywhere for ω = 0 (4.72).
- **Check.** Rankine vortex: along circles B constant ✓; across the core B = r² − 1 varies, outside B = 0 everywhere ✓.
  Cylinder flow: B identical at 1000 points ✓.
- **What it means.** "Bernoulli holds across streamlines" is true only where the flow is irrotational — outside boundary
  layers (Ch. 6, 9). The constant along each streamline of a rotational flow is a label of that streamline.
- **Traps.** Using one constant across a rotational flow. Forgetting the connected-region condition (flow round an
  obstacle with circulation still has one B if ω = 0 everywhere in the fluid, but φ is multivalued).

### D26 · Unsteady Bernoulli for potential flow: (4.73) → (4.74) → (4.75) — ★★, 8 steps, in C12 (notebook · `which_bernoulli`)
- **Goal.** In an irrotational flow keep the time derivative and still get a Bernoulli equation valid everywhere — and fix
  the sign of the book's redefinition of φ.
- **Start.** (4.69) with ω = 0: $\frac{\partial\mathbf u}{\partial t}+\nabla\Big[\tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz\Big]=0$.
- **Plan.** (1) Write u = ∇φ and make ∂u/∂t a gradient. (2) Collect everything under one gradient: the bracket depends on
  t only. (3) Absorb that function into φ with the correct sign.
- **Tools.** Velocity potential (R09) · Schwarz's theorem (primer in C08) · zero gradient ⇒ function of t only (gloss).
- **Assumptions.** Inviscid, barotropic (from D24); irrotational in a simply connected region (step 1).
- **Steps.**
  1. *did:* Use the velocity potential · *tex:* $\mathbf u=\nabla\phi$ · *why:* (4.73), R09: an irrotational flow in a simply connected
     region is the gradient of a potential. · *plain:* The velocity is the slope of φ. · *set:* {sc: 'u_tube', t: 0}.
  2. *did:* Swap time and space derivatives · *tex:* $\frac{\partial\mathbf u}{\partial t}=\frac{\partial}{\partial t}\nabla\phi=\nabla\frac{\partial\phi}{\partial t}$ · *why:* Schwarz's
     theorem (mixed derivatives commute for smooth φ). · *plain:* The acceleration term is a gradient too.
  3. *did:* Insert into the steady-less (4.69) · *tex:* $\nabla\frac{\partial\phi}{\partial t}+\nabla\Big[\tfrac12\lvert\nabla\phi\rvert^2+\int_{p_o}^{p}\frac{dp'}{\rho}+gz\Big]=0$ · *why:*
     ω = 0 removes u × ω; u_iu_i = |∇φ|². · *plain:* Two gradients add to zero.
  4. *did:* Collect under one gradient · *tex:* $\nabla\Big[\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz\Big]=0$ · *why:*
     Linearity of the gradient. This is the first form of (4.74). · *plain:* One bracket has no spatial variation.
  5. *did:* Integrate in space · *tex:* $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=B(t)$ · *why:* A zero gradient
     everywhere means the bracket is the same at every point, but it may still change with time (gloss). The second form of
     (4.74). · *plain:* The bracket is a function of time only.
  6. *did:* Shift φ by a function of time · *tex:* $\phi_{new}=\phi-\int_{t_o}^{t}B(t')\,dt'$ · *why:* The added term has no x-dependence, so
     ∇φ_new = ∇φ (same velocity); ∂φ_new/∂t = ∂φ/∂t − B. The book prints +∫B dt′, which makes the bracket 2B(t): a sign
     slip. · *plain:* Re-label the potential without changing the flow. · *live:* "bracket with −: B − B = 0; with the
     book's +: 2B = …".
  7. *did:* Substitute the shifted potential · *tex:* $\frac{\partial\phi_{new}}{\partial t}+\tfrac12\lvert\nabla\phi_{new}\rvert^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=B(t)-B(t)=0$ ·
     *why:* Step 6 in step 5: the time function cancels. · *plain:* The time function is gone.
  8. *did:* Drop the label · *tex:* $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=\text{constant}$ · *why:* Rename
     φ_new as φ; any constant may stand on the right (a constant shift of B changes nothing). This is (4.75). · *plain:*
     Unsteady Bernoulli holds everywhere at every instant. · *watch:* "the amber ∂φ/∂t bar carries the difference when the
     column stops".
- **Result.** $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=$ constant (4.75).
- **Check.** Units m²/s² ✓. Steady: reduces to (4.72) ✓. U-tube (L = 1 m): p₁ − p₂ = ρL dU/dt = 981 Pa at the turning
  point = ρg·2h₀ ✓. sympy gauge check (C12): minus sign → 0, plus sign → 2B(t) ✓.
- **What it means.** The pressure of an accelerating irrotational flow: Ch. 7's dynamic free-surface condition, added mass
  of accelerating bodies (Ch. 6), bubbles and acoustics (Ch. 15). It needs irrotationality — unlike (4.71).
- **Traps.** The gauge sign (φ − ∫B dt′). Believing the constant cannot depend on time before it is absorbed. Using it in a
  rotational flow (no φ exists).

### D27 · The Boussinesq momentum equation: (4.39b) − hydrostatics → (4.84) → (4.86) — ★★, 9 steps, in C13 (notebook · `boussinesq_buoyancy`)
- **Goal.** Show why small density differences can be ignored in the fluid's inertia but not where they are multiplied by
  gravity.
- **Start.** (4.39b): $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$.
- **Plan.** (1) Subtract the hydrostatic reference state. (2) Divide by ρ_s. (3) Compare sizes: drop ρ′ in the inertia, keep
  it with g. (4) Replace ρ_s by a constant.
- **Tools.** Hydrostatic reference state (R10) · order-of-magnitude scaling (primer in C13) · reduced gravity (primer in
  C13).
- **Assumptions.** Density changes small, ρ′/ρ_s ~ αδT ≪ 1 (step 6); g ≫ the flow's own accelerations (step 7); base-state
  variation small (step 8); constant μ.
- **Steps.**
  1. *did:* Start from Navier–Stokes · *tex:* $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ · *why:* (4.39b), the incompressible
     equation (continuity is ∇·u = 0 to order αδT, N108). · *plain:* Newton's law for the fluid. · *set:* {sc: 'lake'}.
  2. *did:* Write the resting state · *tex:* $0=-\nabla p_s+\rho_s\mathbf g$ · *why:* R10: a fluid at rest balances pressure against
     gravity; p_s(z) and ρ_s(z) may vary with height. · *plain:* The hydrostatic balance. · *watch:* "p_s(z) is drawn beside
     the tank".
  3. *did:* Subtract the resting state · *tex:* $\rho\frac{D\mathbf u}{Dt}=-\nabla(p-p_s)+(\rho-\rho_s)\mathbf g+\mu\nabla^2\mathbf u$ · *why:* Subtracting
     zero; the large hydrostatic parts cancel. · *plain:* Only departures from rest drive motion.
  4. *did:* Name the departures · *tex:* $\rho\frac{D\mathbf u}{Dt}=-\nabla p'+\rho'\mathbf g+\mu\nabla^2\mathbf u,\qquad p'=p-p_s,\ \rho'=\rho-\rho_s$ · *why:*
     This is (4.84); primes now mean perturbations (a new meaning of ′). · *plain:* Perturbation pressure and density.
  5. *did:* Divide by ρ_s · *tex:* $\frac{\rho}{\rho_s}\frac{D\mathbf u}{Dt}=-\frac1{\rho_s}\nabla p'+\frac{\rho'}{\rho_s}\mathbf g+\frac{\mu}{\rho_s}\nabla^2\mathbf u$ · *why:* To compare
     the sizes of the terms per unit mass (the book's move). · *plain:* Every term an acceleration.
  6. *did:* Estimate the density ratio · *tex:* $\frac{\rho}{\rho_s}=1+\frac{\rho'}{\rho_s},\qquad\frac{\rho'}{\rho_s}\sim\alpha\,\delta T\sim10^{-3}$ · *why:* Order-of-
     magnitude scaling (primer): δρ/ρ = αδT for thermal density changes (water, δT = 5 K: 10⁻³). · *plain:* Density differs
     from its resting value by a tenth of a percent.
  7. *did:* Keep ρ′ only with g · *tex:* $\Big(1+\frac{\rho'}{\rho_s}\Big)\frac{D\mathbf u}{Dt}\approx\frac{D\mathbf u}{Dt},\qquad\frac{\rho'}{\rho_s}\mathbf g\ \text{kept}$ · *why:* The book
     skips why: ρ′Du/Dt/ρ_s is 10⁻³ of Du/Dt — negligible; but ρ′g/ρ_s is compared with the other accelerations, and
     g ≫ |Du/Dt|, so it is as large as they are. · *plain:* Tiny density differences matter only through gravity. · *live:*
     "ρ′g/ρ₀ = 0.0098 m/s² vs (ρ′/ρ₀)Dw/Dt ≈ 10⁻³ × 0.01 = 10⁻⁵ m/s²" · *watch:* "the grey bar is a thousandth of the
     blue one".
  8. *did:* Replace ρ_s by a constant ρ₀ · *tex:* $\frac1{\rho_s}\approx\frac1{\rho_0},\qquad\frac{\mu}{\rho_0}=\nu$ · *why:* The base state varies by a
     similarly small fraction over the depth (L ≪ c²/g); ρ₀ is a reference value, and μ/ρ₀ is the kinematic viscosity. ·
     *plain:* One constant density everywhere except in the buoyancy.
  9. *did:* Write the Boussinesq momentum equation · *tex:* $\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u$ · *why:*
     Steps 7–8 in step 5. This is (4.86). · *plain:* Density appears only in the buoyancy term.
- **Result.** $\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u$ (4.86) — *in words:* the Boussinesq momentum
  equation.
- **Check.** Units m/s² ✓. ρ′ = 0: gravity disappears (4.85) ✓. Lake: buoyancy g′ = 0.0196 m/s² for δT = 10 K (C13) ✓;
  dropped term 10⁻³ of the inertia ✓.
- **What it means.** The momentum equation of the ocean and atmosphere models of Ch. 13, of internal waves (Ch. 7) and
  convection (Ch. 11). It fails when density differences are large (αδT ~ 1) or the layer is as deep as c²/g.
- **Traps.** Dropping ρ′ everywhere (no buoyancy, no convection). Keeping ρ′ in the inertia (inconsistent order). Confusing
  p′ (perturbation) with the rotating frame's primes.

### D28 · The Boussinesq heat equation: (4.87) → (4.88) → (4.89) — ★★, 8 steps, in C13 (notebook · `boussinesq_buoyancy`)
- **Goal.** Show why the Boussinesq heat equation has C_p (not C_v) although continuity is ∇·u ≈ 0, and reach
  DT/Dt = κ∇²T.
- **Start.** (4.87): $\rho\frac{De}{Dt}=-p\nabla\cdot\mathbf u+\rho\varepsilon-\nabla\cdot\mathbf q$ — *in words:* (4.60), the internal-energy equation, in
  vector form.
- **Plan.** (1) Express the pressure work through continuity and the thermal expansion. (2) Use perfect-gas relations to
  turn it into ρ(C_p − C_v)DT/Dt. (3) Move it left: C_v → C_p. (4) Drop viscous heating, use Fourier.
- **Tools.** Continuity (4.8) (C02) · thermal expansion coefficient α (Ch. 1 §1.9) · perfect gas p = ρRT, R = C_p − C_v,
  α = 1/T (ch01 P33) · order-of-magnitude scaling (primer in C13) · Fourier's law (Ch. 1 §1.5).
- **Assumptions.** Density changes caused by temperature alone (step 3); perfect gas (steps 4–5; liquids: see the trap);
  low Eckert number (step 7); constant k (step 8).
- **Steps.**
  1. *did:* Start from the energy equation · *tex:* $\rho\frac{De}{Dt}=-p\nabla\cdot\mathbf u+\rho\varepsilon-\nabla\cdot\mathbf q$ · *why:* (4.87): (4.60) with
     ε inserted. We keep −p∇·u: although ∇·u ≈ 0 in continuity, p is large. · *plain:* Internal energy: compression work,
     viscous heating, conduction. · *set:* {sc: 'abl'}.
  2. *did:* Replace ∇·u with continuity · *tex:* $-p\nabla\cdot\mathbf u=\frac p\rho\frac{D\rho}{Dt}$ · *why:* Continuity (4.8) gives ∇·u = −(1/ρ)Dρ/Dt; multiply both sides by −p. · *plain:* The
     tiny expansion, written with the density rate.
  3. *did:* Let temperature set the density · *tex:* $\frac{D\rho}{Dt}=\Big(\frac{\partial\rho}{\partial T}\Big)_p\frac{DT}{Dt}=-\rho\alpha\frac{DT}{Dt}$ · *why:* With
     pressure effects negligible (L ≪ c²/g), ρ changes with T; α = −(1/ρ)(∂ρ/∂T)_p (1.20). · *plain:* Warming a particle
     lowers its density.
  4. *did:* Use the perfect gas · *tex:* $-p\nabla\cdot\mathbf u=-p\alpha\frac{DT}{Dt}=-\rho RT\cdot\frac1T\frac{DT}{Dt}=-\rho R\frac{DT}{Dt}$ · *why:* p = ρRT and
     α = 1/T for a perfect gas (ch01 P33). · *plain:* The pressure work is ρR times the warming rate. · *live:* "air:
     −p∇·u = −1.2 × 287 × DT/Dt".
  5. *did:* Use R = C_p − C_v · *tex:* $-p\nabla\cdot\mathbf u=-\rho(C_p-C_v)\frac{DT}{Dt}$ · *why:* Mayer's relation for a perfect gas (Ch. 1
     §1.8). · *plain:* It is the difference of the two specific heats, times the warming rate.
  6. *did:* Move the pressure work left, e = C_vT · *tex:* $\rho C_v\frac{DT}{Dt}+\rho(C_p-C_v)\frac{DT}{Dt}=\rho C_p\frac{DT}{Dt}=\rho\varepsilon-\nabla\cdot\mathbf q$ · *why:*
     de = C_v dT for a perfect gas; adding ρ(C_p − C_v)DT/Dt to both sides turns C_v into C_p. This is (4.88). · *plain:*
     The small expansion's work is what makes the heat capacity C_p. · *watch:* "the pressure-work bar folds into the
     heating: C_v becomes C_p".
  7. *did:* Drop the viscous heating · *tex:* $\frac{\rho\varepsilon}{\rho C_p(DT/Dt)}\sim\frac{\nu U}{C_p\,\delta T\,L}\ll1$ · *why:* Order-of-magnitude
     scaling (N112): water, U = 0.1 m/s, L = 1 m, δT = 1 K gives ~2×10⁻¹¹ — the low-Eckert condition of C15. · *plain:*
     Friction heating is far too small to matter here.
  8. *did:* Use Fourier's law, constant k · *tex:* $\frac{DT}{Dt}=\kappa\nabla^2T,\qquad\kappa\equiv\frac{k}{\rho C_p}$ · *why:* q = −k∇T (1.2) with
     uniform k gives −∇·q = k∇²T; divide by ρC_p. This is (4.89). · *plain:* Temperature is carried by the flow and diffuses.
- **Result.** $\rho C_p\frac{DT}{Dt}=\rho\varepsilon-\nabla\cdot\mathbf q$ (4.88) and $\frac{DT}{Dt}=\kappa\nabla^2T$ (4.89).
- **Check.** Units K/s ✓. U = 0: the heat equation of Ch. 1 ✓. The Gaussian blob with σ² = σ₀² + 2κt solves it exactly (C13
  code, residual ≈ 0) ✓. sympy (C13 note N110): e = C_vT, p = ρRT, α = 1/T ⇒ ρC_v DT/Dt + p∇·u = ρC_p DT/Dt ✓.
- **What it means.** The heat (and, with salinity, density) equation of ocean and atmosphere models. C_p appears because
  a warming parcel expands a little against a large pressure.
- **Traps.** Writing C_v because "∇·u = 0". For liquids the argument is different: p∇·u is small but C_p ≈ C_v anyway — the
  ocean uses (4.89) with water's C_p. Dropping ρε in a fast lubricating flow (high Eckert number: E7's oil bearing).

### D29 · The kinematic boundary condition: (4.90) → (4.91) → (4.92) → (4.93) — ★★, 8 steps, in C14 (notebook)
- **Goal.** Turn "no fluid crosses this moving surface" into an equation, and find the mass flux when fluid does cross.
- **Start.** A surface described as a level set $\eta(\mathbf x,t)=0$ — *in words:* the free surface, an interface or a wall is where
  a function η vanishes.
- **Plan.** (1) Follow a point riding on the surface. (2) Get its normal speed. (3) Demand the fluid's normal speed be the
  same. (4) Allow a relative velocity for leaky surfaces.
- **Tools.** Moving level set and its normal speed (primer in C14) · chain rule along a path (ch03 P91) · unit normal
  ∇η/‖∇η‖ (ch02 P75).
- **Assumptions.** η smooth with ∇η ≠ 0 on the surface (step 4); no fluid crosses (steps 5–7).
- **Steps.**
  1. *did:* Describe the surface as a level set · *tex:* $\eta(\mathbf x,t)=0$ · *why:* Any moving, deforming surface can be written as
     the zero set of a function (primer); e.g. a wave η = z − a cos(kx − ωt). · *plain:* The surface is where η is zero.
  2. *did:* Follow a point riding on it · *tex:* $\eta\big(\mathbf x_s(t),t\big)=0\ \text{for all }t,\qquad\frac{d\mathbf x_s}{dt}=\mathbf u_s$ · *why:* A point that
     stays on the surface keeps η = 0 as time runs; u_s is its velocity. · *plain:* A marker glued to the surface.
  3. *did:* Differentiate along its path · *tex:* $\frac{d\eta}{dt}=\frac{\partial\eta}{\partial t}+(\mathbf u_s\cdot\nabla)\eta=0$ · *why:* Chain rule along a
     path (ch03 P91); the derivative of a constant (0) is zero. This is (4.90). · *plain:* The marker sees η stay zero.
  4. *did:* Get the surface's normal speed · *tex:* $\mathbf u_s\cdot\mathbf n=-\frac{\partial\eta/\partial t}{\lvert\nabla\eta\rvert},\qquad\mathbf n=\frac{\nabla\eta}{\lvert\nabla\eta\rvert}$ ·
     *why:* Divide (4.90) by |∇η| ≠ 0. Only the normal part of u_s is fixed — sliding along the surface does not move it. ·
     *plain:* How fast the surface moves along its own normal. · *live:* "moving wall η = x − Vt: u_s·n = −(−1)/1 = 1.0 m/s".
  5. *did:* Demand no flow across · *tex:* $\mathbf u\cdot\mathbf n=\mathbf u_s\cdot\mathbf n$ · *why:* Fluid particles on the surface stay on it: their
     velocity relative to the surface has no normal component. · *plain:* The fluid moves normal to the surface exactly as
     fast as the surface.
  6. *did:* Multiply by the gradient's size · *tex:* $\mathbf u\cdot\nabla\eta=\mathbf u_s\cdot\nabla\eta=-\frac{\partial\eta}{\partial t}$ · *why:* Multiply step 5 by the gradient's size, n times its length being ∇η; then use (4.90). · *plain:* Written with η's gradient.
  7. *did:* Rearrange · *tex:* $\frac{\partial\eta}{\partial t}+(\mathbf u\cdot\nabla)\eta\equiv\frac{D\eta}{Dt}=0\quad\text{on}\quad\eta=0$ · *why:* Move the term left;
     it is the material derivative of η (3.5). This is (4.91). · *plain:* η stays zero for every fluid particle on the
     surface.
  8. *did:* Allow a flux across · *tex:* $(u_{rel})_n=\mathbf u\cdot\mathbf n-\mathbf u_s\cdot\mathbf n=\frac{1}{\lvert\nabla\eta\rvert}\frac{D\eta}{Dt}$ · *why:* For evaporation or a
     shock the relative normal velocity is not zero; steps 4–6 give it. This is (4.92); times ρ it is the mass flux per area
     (4.93). · *plain:* Dη/Dt measures how fast fluid crosses the surface.
- **Result.** $\partial\eta/\partial t+(\mathbf u\cdot\nabla)\eta\equiv D\eta/Dt=0$ on η = 0 (4.91); with crossing, $(u_{rel})_n=(1/\lvert\nabla\eta\rvert)D\eta/Dt$
  (4.92) and mass flux $(\rho/\lvert\nabla\eta\rvert)D\eta/Dt$ (4.93).
- **Check.** Units: η/s ✓. Moving wall: Dη/Dt = −V + u·e_x = 0 when u_x = V ✓. Linear wave: the linearised condition
  ∂η/∂t = w at z = 0 is exact for the linear solution; the full (4.91) at the true surface leaves a residual ~ (ka)² in
  units of the phase speed ✓.
- **What it means.** Every free-surface problem (Ch. 7), every layer interface (Ch. 13) and every moving shock (Ch. 15)
  uses it. It fixes only the normal velocity; tangential conditions come from stresses or no-slip.
- **Traps.** Thinking u_s is fully defined (only its normal part is). The orientation of n (toward increasing η). Taking the
  linearised condition as exact at finite amplitude.

### D30 · Dimensionless Navier–Stokes: (4.100) → (4.101) — ★★, 8 steps, in C15 (notebook · `dynamic_similarity_models`)
- **Goal.** Scale every variable by the problem's own sizes and find which combinations of them the flow can depend on.
- **Start.** (4.39b) written out: $\rho\Big(\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u\Big)=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$.
- **Plan.** (1) Introduce scaled variables. (2) Scale the derivatives with the chain rule. (3) Substitute and divide by the
  size of the advective term. (4) Name the groups; note the boundary conditions scale too.
- **Tools.** Scaled variables and the chain rule (primer in C15) · order-of-magnitude scaling (primer in C13) · sympy
  subs/collect (primer in C05).
- **Assumptions.** Incompressible, constant ρ and μ (the starting equation); one length l, speed U, frequency Ω describe the
  problem.
- **Steps.**
  1. *did:* Start from Navier–Stokes · *tex:* $\rho\frac{\partial\mathbf u}{\partial t}+\rho(\mathbf u\cdot\nabla)\mathbf u=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ · *why:* (4.39b)
     with D/Dt written out (the book's §4.11 re-display). · *plain:* Five terms with dimensions. · *set:* {mode: 'sphere'}.
  2. *did:* Scale the variables · *tex:* $\mathbf x=l\mathbf x^*,\quad t=\frac{t^*}{\Omega},\quad\mathbf u=U\mathbf u^*,\quad p=p_\infty+\rho U^2p^*,\quad\mathbf g=g\,\mathbf g^*$ · *why:*
     (4.100): starred variables are pure numbers of order one when l, U, Ω, ρU² are the right sizes. · *plain:* Measure
     everything in the problem's own units.
  3. *did:* Scale the derivatives · *tex:* $\frac{\partial}{\partial t}=\Omega\frac{\partial}{\partial t^*},\qquad\nabla=\frac1l\nabla^*,\qquad\nabla^2=\frac1{l^2}\nabla^{*2}$ · *why:*
     Chain rule (primer): t* = Ωt gives ∂/∂t = Ω∂/∂t*; x* = x/l gives ∂/∂x = (1/l)∂/∂x*. · *plain:* Each derivative brings
     out its scale.
  4. *did:* Substitute term by term · *tex:* $\rho U\Omega\frac{\partial\mathbf u^*}{\partial t^*}+\frac{\rho U^2}{l}(\mathbf u^*\cdot\nabla^*)\mathbf u^*=-\frac{\rho U^2}{l}\nabla^*p^*+\rho g\,\mathbf g^*+\frac{\mu U}{l^2}\nabla^{*2}\mathbf u^*$ ·
     *why:* Steps 2–3 in step 1; p∞ is a constant, so its gradient vanishes — only pressure differences matter. · *plain:*
     Each term = (a size) × (a term of order one). · *watch:* "each term lights its bar".
  5. *did:* Divide by ρU²/l · *tex:* $\Big[\frac{\Omega l}{U}\Big]\frac{\partial\mathbf u^*}{\partial t^*}+(\mathbf u^*\cdot\nabla^*)\mathbf u^*=-\nabla^*p^*+\Big[\frac{gl}{U^2}\Big]\mathbf g^*+\Big[\frac{\mu}{\rho Ul}\Big]\nabla^{*2}\mathbf u^*$ ·
     *why:* The advective term's size ρU²/l is the natural yardstick (and the pressure scale ρU² was chosen to match it).
     This is (4.101). · *plain:* Three brackets remain.
  6. *did:* Name the groups · *tex:* $\mathrm{St}=\frac{\Omega l}{U},\qquad\frac1{\mathrm{Fr}^2}=\frac{gl}{U^2},\qquad\frac1{\mathrm{Re}}=\frac{\mu}{\rho Ul}$ · *why:* (4.102)–(4.104): each
     bracket is a ratio of two terms' sizes — unsteady/advective, gravity/inertia, viscous/inertia. · *plain:* Strouhal,
     Froude and Reynolds numbers. · *live:* "St = …, gl/U² = …, μ/ρUl = …".
  7. *did:* Scale the boundary conditions · *tex:* $\mathbf u^*\to\mathbf e_x\ \text{far away},\qquad\mathbf u^*=0\ \text{on the body}\ (r^*=1)$ · *why:* For a
     body of size l in a stream U, the conditions in starred variables contain no l, U or Ω — only the shape. · *plain:*
     The edges look the same for every size.
  8. *did:* Conclude dynamic similarity · *tex:* $\text{equal St, Fr, Re}+\text{same shape}\ \Rightarrow\ \text{same }\mathbf u^*(\mathbf x^*,t^*),\ p^*$ · *why:* The same
     dimensionless equation with the same dimensionless conditions has the same solution. · *plain:* Two flows with equal
     groups are the same flow at different scales. · *watch:* "equal brackets ⇒ same equation ⇒ same solution".
- **Result.** $\big[\frac{\Omega l}{U}\big]\frac{\partial\mathbf u^*}{\partial t^*}+(\mathbf u^*\cdot\nabla^*)\mathbf u^*=-\nabla^*p^*+\big[\frac{gl}{U^2}\big]\mathbf g^*+\big[\frac{\mu}{\rho Ul}\big]\nabla^{*2}\mathbf u^*$ (4.101)
  — *in words:* the equation knows only St, Fr and Re.
- **Check.** Every bracket dimensionless (Ω l/U: s⁻¹·m/(m/s) ✓; gl/U² ✓; μ/(ρUl): Pa s/(kg/m³·m/s·m) ✓). Stokes-layer
  collapse (C15 figure) ✓. sympy (`check_src` optional): `ch04.nondimensional_ns_coefficients()` equals the brackets.
  Other pressure scales (μU/l, ρgl) give products of the same groups (Exercise 4.59) ✓.
- **What it means.** Model testing works when the groups that matter match; Re, Fr, St (and later Ri, Ro, M) name the
  regimes of every later chapter. When two groups cannot both match (a ship model: Fr and Re), each effect is scaled by
  its own law.
- **Traps.** Forgetting that the time and pressure scales can change (Ωt vs Ut/l; ρU² vs μU/l). Believing geometric
  similarity is enough. Book slip: the text calls the brackets "in (4.100)"; they are in (4.101).
