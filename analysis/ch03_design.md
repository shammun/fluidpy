# Chapter 3 — Kinematics: lesson design
(from `analysis/ch03_curation.md` (A 15 · B 52 · C 12 = 79 rows; CORE 15 · NOTE 55 · RECAP 7 · SKIP 2; 24 derivations
D01–D24 written out; E1–E7 + backup B1) and `analysis/ch03.md` (§2b derivations, §4 rows 1–53, §9 conventions and
typos); every equation re-read on the rendered pages p096–p116 (printed pp. 69–89: (3.1)–(3.5) p097, (3.6)–(3.7) p098,
(3.8) p099, Ex. 3.1 p100, (3.9) p101, (3.10)–(3.13) p103, linear strain p104, shear strain + (3.14) p105, (3.15)–(3.18)
p106, (3.19)–(3.20) p107, (3.21) p108, §3.5 shear flow p109, (3.22)–(3.23) p110, (3.24)–(3.27) p111, (3.28)–(3.29)
p112, (3.30) p113, (3.31)–(3.33) p114, (3.34)–(3.35) + Ex. 3.2 p115–p116); 2026-09-22, lesson-designer. The implementer
works in parallel from analysis §4 + curation §8; **Part C is the contract both sides keep.**)

**Binding conventions for every builder (analysis §9, curation decisions 4–5).**
1. **Velocity fields are callables** `u(x, t) -> ndarray` with `x` of shape `(d,)` or `(d, N)` (d = 2 or 3; 2-D fields
   never pad silently — a 2-D field takes 2-D points) returning the same shape; scalar fields `F(x, t)` return a float
   or `(N,)`. Every public function is scalar-callable (scalar in → float out, lists accepted) so explainer parity rows
   can call it (`py:` expressions have **no builtins**: index tuples/dicts down to one number).
2. **Velocity gradient** `G[i, j] = ∂u_i/∂x_j` (row = velocity component). Book **R = G − Gᵀ** (no ½), its vector is the
   vorticity **ω = ∇×u** (`vector_from_antisymmetric(rotation_tensor(G))`), and the **element spin is ½ω** — never
   written "ω". G = S + ½R (3.11).
3. **Three gammas** (⚠️ where first used): ch02's Γ ≡ S₁₂ (Ex. 2.4); ch03's shear rate **γ = du₁/dx₂ = 2S₁₂** (code
   `gamma`; `velocity_gradient_preset("simple_shear", Gamma=gamma)` gives G = [[0, γ], [0, 0]]); **Γ = circulation**
   [m²/s] (code `Gamma` in `core.vortices`, the E6/C13/C14 symbol). ω₃ = −γ is clockwise.
4. **Angles.** Plane polar (r, θ) with θ counterclockwise from +x; cylindrical (R, φ, z); spherical (r, θ, φ) with θ the
   polar angle from +z and φ the azimuth (physics convention, Fig. 3.3d). Ex. 3.2's θ is the **cone half-angle**
   (tan θ = r_o/h) — a fourth θ, said in N54. Angles in radians in every function; degrees only on explainer controls.
5. **Ex. 3.1 names** (four symbols that collide in print): drawing instant **t′** (`t_prime`), release time **t_o**
   (`t_release`), integration constants **c_x, c_y** (the book's x_o, y_o), amplitude **ξ_o** (`xi0`).
6. **Material-line angles (Figs. 3.11, D09, D12):** α = clockwise tilt of the side that starts vertical, β =
   counterclockwise tilt of the side that starts horizontal; closing rate α + β (strain), counterclockwise spin
   ½(−α + β). Code measures signed counterclockwise angles and derives both.
7. **RTT sign:** ΔV and the swept sliver (bΔt)·n dA are **signed** (b·n < 0 removes volume); one formula covers growing and
   shrinking boundaries.
8. **Galilean transformation exactly as the book (p. 74):** frame O′x′y′z′ moves at constant U relative to Oxyz;
   t = t′, **x = x′ + Ut + x′_o**, **u(x, t) = U + u′(x′, t′)**. `galilean_transform(u, U, x0p)` returns
   u′(x′, t′) = u(x′ + Ut′ + x′_o, t′) − U. Galilean ≠ rotating frame (N28, Ch. 4 §4.7).
9. **Book typos taught corrected** (curation decision 4): (3.6) second term is |u| ∂F/∂s (F missing in print); Ex. 3.2:
   b·n = 0 on the base (not b = 0), "[?]" in the integrand is the factor z; §3.4's "Section 2.12" means §2.11; (3.19)
   follows from (3.10), (3.11) and **(3.15)** (the text says (3.14)); the eigenvalue −γ/2 is "compression at rate γ/2";
   Fig. 3.16's C is B.
10. **Colours (text, figures, explainers — one meaning each):** streamline `teal` · path line `orange` · streak line
    `rose` · local term ∂/∂t `blue` · advective term `amber` · total D/Dt `accent` purple · strain part `teal` ·
    rotation part `orange` · stretching `blue` / compressing `rose` · u_θ `teal`, ω_z `orange`, Γ `accent` · volume term
    `blue`, surface term `orange`, total purple · swept band b·n > 0 `blue`, b·n < 0 `rose` · ghosts and references
    `muted` grey.
11. **Every book equation is shown in full next to its number** (CLAUDE.md rule 3) — in markdown, derivation steps
    ("substitute (3.5), $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$"), traps, recaps, explainer tour,
    Explain, Derivation and quiz text. A repeat mention in the same cell may use a compact inline form, never the bare
    number. `tools/eq_refs.py` checks the explainers.
12. **nbkit behaviour the builder must respect:** `nb.recap(...)` and `nb.section(...)` **end** the current CORE block
    (they reset `current_core`), and `nb.core(id)` cannot be called twice — so every RECAP sits **before** the
    `nb.core` call of the block that uses it (R01–R03 before C06; R04–R07 before C10), and every `nb.derivation` sits
    inside its curation CORE block (D14 in C10 in §3.4, even though its flow, N34, is re-stated at the head of §3.5).
    `nb.derivation(ref=…)` gets an equation number only ("3.9"), never "Exercise 3.12" (it prefixes "Eq.").
13. **Parsing.** Part E is the only part with table rows after the Part E heading; Part F has no line starting with
    `|`. Explainer headings are exactly `### E1 · flow_lines_unsteady` … `### E7 · reynolds_transport_cv`, `### B1 ·
    material_volume_divergence` (`tools/embed_check.py`). Primer terms in Part A are the exact Concept text of the Part E
    rows marked "primer" (coverage_check matches the first 18 characters).

Order of parts: C (contract) · A (notebook storyboard) · B (explainers) · D (runtime) · E (prerequisite ledger) · F
(derivations).

---

## Part C — functions the builders will call (the implementer's contract)

**Status column.** **§4 #n** = planned in `analysis/ch03.md` §4 row n with this signature (keep it). **§8** = added by
the curation's notes for the implementer (`analysis/ch03_curation.md` §8). **NEW** = added by this design (in neither) —
flagged as the phase asks. **reuse** = exists (ch01/ch02) and is only called. **promote** = moves from ch02 to
`core.kinematics` with a re-export from `ch02_cartesian_tensors` (analysis §4, curation §8 last bullet). Every callable
is reachable as **`ch03.<name>`**: the chapter module `fluidpy/ch03_kinematics.py` re-exports `core.kinematics`,
`core.coords`, `core.vortices`, `core.transport` and the reused `core.tensors` / `core.operators` /
`core.integral_theorems` names listed in C.0 (so the notebook imports one module and parity rows write `ch03.…`).
Module aliases in the notebook: `from fluidpy import ch03_kinematics as ch03`. Units SI: x [m], t [s], u [m/s], G [1/s],
Γ [m²/s], T [K]; docstrings cite § and Eq. and carry the validation label.

### C.0 Reused (existing; called, not changed)
| # | Callable (signature) | Returns | Used by | Status |
|---|---|---|---|---|
| 0.1 | `style.setup_notebook() -> bool` (FAST) · `style.COLORS` · `style.savefig(fig, "ch03", name)` | FAST · palette · path | setup, every figure | reuse |
| 0.2 | `anim.animate(update, frames, fig, interval)` · `anim.show_animation(anim, player="video"\|"frames")` | HTML | C04, §3.5 head (Fig. 3.14), C12, C13, C15 | reuse |
| 0.3 | `interact.slider_figure(fn, name, values, *, unit, xlabel, ylabel, title, xrange, yrange, height)` · `interact.animate_figure` | plotly Figure | C02, C03, C05, C10, C14, C15 | reuse |
| 0.4 | `embed.show_viz("ch03", slug)` | display | 7 explainer cells | reuse |
| 0.5 | `tools.convergence.observed_order(h, err) -> float` | slope | C06 (Taylor remainder), C15 (Δt terms) | reuse |
| 0.6 | `core.tensors`: `strain_rate_tensor(G)`, `rotation_tensor(G)` (= G − Gᵀ), `antisymmetric_part(G)`, `antisymmetric_from_vector(omega)` (= −ε·ω), `vector_from_antisymmetric(R)`, `principal_axes(S)` (λ ascending, columns = unit eigenvectors, det +1), `trace(A)`, `levi_civita()`, `cross(u, v)`, `transform_tensor(T, C)`, `rotation_matrix_2d(theta)` | arrays | R01–R06, C06–C12 | reuse |
| 0.7 | `core.operators`: `curl(u, h)`, `divergence(u, h)`, `vector_gradient(u, h)`, `is_irrotational(u, h)` | grid arrays | C10 (N29), C13 cross-check | reuse |
| 0.8 | `core.integral_theorems`: `circulation(u_fn, loop)`, `planar_loop(center, normal=None, radius, n)`, `rectangle_loop(...)`, `planar_disc(...)`, `divergence_theorem_sphere(Q_fn, R, n)` | floats, Loop | R07, N37 off-centre circle, D23 check | reuse |
| 0.9 | `core.index_notation.expand_indices_str("u_i dF_dx_i")` (whatever spelling the ch02 parser accepts for a derivative factor; the WIP accepts `"u_i dF/dx_i"` per analysis §2 #13 — builder checks once) | str | C02 (N11) | reuse |
| 0.10 | `ch02.solid_body_rotation_field(b)`, `ch02.shear_field(Gamma)`, `ch02.irrotational_vortex_field(K)` (`VectorField`s with `singular_at`) | fields | C10, C13 cross-checks | reuse |

### C.1 `fluidpy/core/kinematics.py` (NEW module, `core.K`)
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 1.1 | `lagrangian_velocity_acceleration(r_of_t, t, h=1e-4) -> (u, a)` (2nd-order central differences of a callable r(t)) · `lagrangian_velocity_acceleration_sym(r_exprs, t) -> (u_exprs, a_exprs)` | (3.1) | m/s, m/s² | C01 (N08), E2 mode 2 | §4 #8 |
| 1.2 | `lagrangian_to_eulerian(r_exprs, labels, coords, t) -> dict(label_of_x, u, a_lagrangian, Du_Dt)` — sympy: solve x = r(t; labels) for the labels (`label_of_x`), substitute into dr/dt → `u(coords, t)`, into d²r/dt² → `a_lagrangian`, and compute `Du_Dt` = ∂u/∂t + (u·∇)u from the Eulerian u (so `simplify(Du_Dt − a_lagrangian) == 0` is the D01 check) | (3.1), (3.2), D01 | sympy dict | C01 (D01 check) | §4 #9 (dict keys fixed here) |
| 1.3 | `material_derivative(F, u, x, t, h=1e-4, ht=None) -> float` | (3.4)/(3.5) | [F]/s | C02 | §4 #10 |
| 1.4 | `material_derivative_terms(F, u, x, t, h=1e-4) -> (local, advective, total)` | (3.5), N10 | [F]/s | C02, E2 | §4 #11 |
| 1.5 | `material_derivative_sym(F_expr, u_exprs, coords, t) -> sympy.Expr` (= ∂F/∂t + u_i ∂F/∂x_i) | (3.5) | sympy | C02 (D02 check, N09) | §4 #12 |
| 1.6 | `streamwise_derivative(F, u, x, t, h=1e-4) -> float` (= \|u\| e_u·∇F; raises `ValueError` at \|u\| = 0) | (3.6) corrected | [F]/s | C02 (N12) | §4 #13 |
| 1.7 | `streamline(u, x0, t_frozen, s_max, both=True, n=400, rtol=1e-10) -> ndarray (d, m)` (arc-length ODE dx/ds = u/\|u\| at frozen t; terminal event \|u\| < 1e-12 — a stagnation point ends the curve and m < n) | (3.7), D03 | m | C03, C04 (N13), N15 | §4 #14 |
| 1.8 | `pathline(u, r0, t0, t_eval, rtol=1e-10, atol=1e-12) -> ndarray (d, nt)` (`solve_ivp` DOP853, `dense_output`; t_eval may lie before t0) | (3.8) | m | C04, C13 animation, E1 parity | §4 #16 |
| 1.9 | `streakline(u, x0, t, t_release) -> ndarray (d, n_release)` (one vector ODE for the ensemble; column k = position at time t of the particle released from x0 at t_release[k] ≤ t) | N17, D05 | m | C04, E1 parity | §4 #17 |
| 1.10 | `galilean_transform(u, U, x0p=None) -> callable u_prime(xp, tp)` (= u(xp + U tp + x0p, tp) − U; x0p defaults to zeros) | N21, (3.9) | m/s | C05 | §4 #19 |
| 1.11 | `acceleration(u, x, t, h=1e-4) -> (a, local, advective)` (∂u/∂t and (u·∇)u by stencils) | (3.9), D02 step 9 | m/s² | C05, E3 parity | §4 #20 |
| 1.12 | `velocity_gradient_at(u, x, t, h=1e-4) -> ndarray (d, d)` (G[i, j] = ∂u_i/∂x_j, 2nd-order central stencils) | (3.10) | 1/s | C06, C08 (D10 check), C13 check | §4 #21 |
| 1.13 | `relative_velocity(G, dx) -> du` (= G·dx) | (3.10) | m/s | C06, E4 | §4 #22 |
| 1.14 | `vorticity_from_gradient(G) -> ndarray (3,)` (2×2 G → (0, 0, ω₃)) · `vorticity(u, x, t, h=1e-4) -> ndarray (3,)` | (3.15), (3.16) | 1/s | C10, R05, R06 | §4 #23, #24 |
| 1.15 | `linear_strain_rate(G, n) -> float` (n·S·n, n normalised inside) · `shear_strain_rate(G, n1, n2) -> float` (n1·S·n2; raises if \|n1·n2\| > 1e-9) | linear / shear strain rate, D08, D09 | 1/s | C07, C08, E4 | §4 #25, #26 |
| 1.16 | `volumetric_strain_rate(G) -> float` (tr G) · `material_volume_ratio(G, t) -> float` (det e^{Gt}) | (3.14), D11 | 1/s; – | C09, E4, B1 | §4 #28 |
| 1.17 | `element_rotation_rate(G) -> ndarray (3,)` (= ½ω) · `material_line_rotation_rate(G, theta) -> float` (2-D: θ̇ = e_θ·G·e(θ); vectorised in theta) | D12, D14 | rad/s | C10, E5 | §4 #29 |
| 1.18 | `vorticity_in_rotating_frame(omega, Omega) -> omega - 2*Omega` (scalars or 3-vectors) | N28, D13 | 1/s | C10, E5 | §4 #30 |
| 1.19 | `velocity_potential_2d(u, grid, x_ref=(0.0, 0.0), t=0.0) -> (phi, path_dependence)` (cumulative trapezoid line integrals from x_ref along two different grid routes; `path_dependence` = max \|φ_route1 − φ_route2\| — reported, never asserted; ≈ 2πB for the line vortex) | (3.17), N29 | m²/s | C10 (N29) | §4 #31 (return pair fixed here) |
| 1.20 | `relative_velocity_split(G, dx) -> (du, du_strain, du_rot)` (S·dx and ½ω × dx; 2-D inputs return 2-D parts) | (3.19), D15 | m/s | C11, E5 | §4 #32 |
| 1.21 | `principal_strain_rates(G) -> (lam, axes)` (`principal_axes(strain_rate_tensor(G))`: λ ascending, columns unit eigenvectors) | (3.20)–(3.21) | 1/s | C12, E5 | §4 #33 |
| 1.22 | `strain_velocity_principal(G, dx) -> (du_bar, dx_bar)` (S·dx and dx in the eigenframe) | (3.21) | m/s, m | C12 (N32) | §4 #34 |
| 1.23 | `deform_circle(G, t, n=72, radius=1.0) -> ndarray (2, n)` · `deform_sphere(G, t, n_theta=24, n_phi=48, radius=1.0) -> ndarray (3, N)` · `strain_ellipse_axes(G, t, method="first_order") -> (semi_axes, directions)` — semi-axes descending (unit initial radius), direction columns; `method="first_order"`: 1 + λt along S's eigenvectors (D16) · `"strain_only"`: e^{λt} (expm(St)) · `"exact"`: singular values / left singular vectors of expm(Gt) | D16, Fig. 3.13 | – | C12, E5 | §4 #35 (method switch fixed here) |
| 1.24 | `measured_strain_rates(G, n, dt) -> dict(stretch, closing, area)` (tracked with `linear_flow_map`: (1/ℓ)Δℓ/Δt along n, −Δ(angle)/Δt between n and its perpendicular, (1/A)ΔA/Δt; → n·S·n, 2 n₁·S·n₂, tr G as dt → 0) | C07–C09 by measurement | 1/s | C07–C09 from-scratch, E4 | §8 |
| 1.25 | `velocity_gradient_preset(name, Gamma=1.0, dim=2)`, `VELOCITY_GRADIENT_PRESETS`, `linear_flow_map(G, t)`, `deform_square(G, t, n_side=10, half_width=1.0, boundary_only=False)`, `material_line_angle(G, t, theta0=0.0)` | ch02 §2.10 → Ch. 3 | – | C07–C12, E4, E5 | promote |

### C.2 `fluidpy/core/coords.py` (NEW module, `core.C`)
| # | Callable (signature) | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 2.1 | `cylindrical_from_cartesian(x, y, z) -> (R, phi, z)` · `cartesian_from_cylindrical(R, phi, z) -> (x, y, z)` | Fig. 3.3c | m, rad | §3.1 (N05) | §4 #3 |
| 2.2 | `spherical_from_cartesian(x, y, z) -> (r, theta, phi)` (θ from +z, φ = arctan2(y, x)) · `cartesian_from_spherical(r, theta, phi)` | Fig. 3.3d | m, rad | §3.1 (N05) | §4 #4 |
| 2.3 | `unit_vectors_cylindrical(phi) -> ndarray (3, 3)` (rows e_R, e_φ, e_z) · `unit_vectors_spherical(theta, phi) -> ndarray (3, 3)` (rows e_r, e_θ, e_φ) | local bases | – | §3.1 (N05), C15 (N54 cone normal) | §4 #5 |
| 2.4 | `velocity_components(u_cart, x, system="cylindrical") -> ndarray` (`"polar"` → (u_r, u_θ) for 2-D x; `"cylindrical"` → (u_R, u_φ, u_z); `"spherical"` → (u_r, u_θ, u_φ)) | projections | m/s | §3.1 (N05), C13 | §4 #6 |

### C.3 `fluidpy/core/vortices.py` (NEW module, `core.X`)
`u_theta` arguments below accept **either a callable u_θ(r) or a kind name** `"solid" | "line" | "rankine" |
"gaussian"` with keywords `Gamma`, `sigma` (this is what makes parity rows possible without lambdas).
| # | Callable (signature) | Implements (Eq.) | Returns | Used by | Status |
|---|---|---|---|---|---|
| 3.1 | `solid_body_rotation(r, omega0) -> u_theta` | (3.22) | m/s | C13, E6 | §4 #37 |
| 3.2 | `line_vortex(r, B) -> u_theta` (inf at r = 0 guarded to NaN) | (3.25) | m/s | C13, E6 | §4 #41 |
| 3.3 | `rankine_vortex(r, Gamma, sigma) -> (u_theta, omega_z)` | (3.28) | m/s, 1/s | C14, E6 | §4 #44 |
| 3.4 | `gaussian_vortex(r, Gamma, sigma) -> (u_theta, omega_z)` (`-np.expm1(-x)/r` form, 0 at r = 0) | (3.29) | m/s, 1/s | C14, E6 | §4 #45 |
| 3.5 | `gaussian_vortex_max_radius(sigma=1.0, method="brentq") -> float` (root of 1 + 2x = eˣ in (0.5, 3); `method="lambertw"` cross-check) | N44, D20 | m | C14, E6 | §4 #46 |
| 3.6 | `vortex_profile(kind, r, Gamma, sigma) -> (u_theta, omega_z)` — dispatcher; `"solid"` uses ω₀ = Γ/(2πσ²) (same speed as Rankine inside), `"line"` uses B = Γ/(2π) (ω_z = 0 for r > 0) — both stated in the docstring | (3.22), (3.25), (3.28), (3.29) | m/s, 1/s | C14 slider, E6 | §8 |
| 3.7 | `vortex_velocity_field(profile, **params) -> callable u(x, t)` (Cartesian u = u_θ(r) e_θ; `profile` a kind name or a callable) | (3.22)–(3.29) | m/s | C13 animation, C13 checks | §4 #38 |
| 3.8 | `polar_vorticity_z(u_r, u_theta, r, theta, h=1e-5) -> float` (callables u_r(r, θ), u_θ(r, θ); 2nd-order stencils) · `polar_vorticity_z_sym(u_r_expr, u_theta_expr, r, theta) -> sympy.Expr` | (3.23), D17 | 1/s | C13, E6 | §4 #39 |
| 3.9 | `circulation_circle(u_theta, r, center=(0.0, 0.0), n=512, **kw) -> float` (centred: 2πr u_θ; off-centre: `core.V.circulation` on `planar_loop` of the Cartesian field) | (3.24), (3.26) | m²/s | C13, C14, E6 | §4 #40 |
| 3.10 | `mean_vorticity_in_disc(u_theta, r, **kw) -> float` (= Γ(r)/(πr²)) | (3.27) | 1/s | C13 (N40), E6 | §4 #42 |

### C.4 `fluidpy/core/transport.py` (NEW module, `core.R`)
Scalar fields for this module: `F(x, t)` with `x` of shape `(d, N)` (d = 2 or 3) or, for the 1-D Leibniz functions, a
float/`(N,)` array of positions.
| # | Callable (signature) | Implements (Eq.) | Returns | Used by | Status |
|---|---|---|---|---|---|
| 4.1 | `leibniz_terms(F, dFdt, a, b, dadt, dbdt, t, n=64) -> (interior, upper, lower, total)` (interior = ∫_a^b ∂F/∂t dx by Gauss–Legendre; upper = ḃF(b, t); lower = ȧF(a, t); total = interior + upper − lower) | (3.30), D21 | [F]·m/s | C15 (N46), E7 | §4 #47 |
| 4.2 | `leibniz_check(F, a_fn, b_fn, t, dt=1e-4) -> float` (central difference of ∫_{a(t)}^{b(t)} F dx by `quad`) | (3.30) left side | [F]·m/s | C15 (D21 check) | §4 #47 |
| 4.3 | `ControlVolume` protocol (`volume(t)`, `volume_nodes(t, n) -> (x (3, N), w (N,))`, `surface_nodes(t, n) -> (x, n_hat, dA, b)`) and shapes with linear-in-time sizes: `GrowingSphere(R0, Rdot, center=(0, 0, 0))` · `GrowingCylinder(R0, Rdot, L)` · `MovingBox(lengths, rates=(0, 0, 0), velocity=(0, 0, 0), origin=(0, 0, 0))` · `GrowingCone(r0, rdot, h)` (apex at the origin, axis +z, base disc at z = h; b·n = 0 on the base) | N48, Fig. 3.18 | – | C15, E7 cone mode | §4 #48 (constructor arguments fixed here) |
| 4.4 | `volume_integral(F, cv, t, n=24) -> float` | ∫_{V*}F dV | [F]·m³ | C15 (N49 demo) | **NEW** (small helper the FD and the sympy-free demos need; ≈ 5 lines on `volume_nodes`) |
| 4.5 | `volume_integral_rate_fd(F, cv, t, dt=1e-4) -> float` | (3.31) left side | [F]·m³/s | C15, E7 parity | §4 #49 |
| 4.6 | `swept_volume_integral(F, cv, t, dt) -> float` · `surface_flux_term(F, cv, t) -> float` (∮F b·n dA) | (3.34) | [F]·m³, [F]·m³/s | C15 (N52 demo) | §4 #50 |
| 4.7 | `reynolds_transport(F, dFdt, cv, t) -> (volume_term, surface_term, total)` · `rtt_check(F, dFdt, cv, t, dt=1e-4) -> (lhs_fd, rhs_total)` | (3.35), D22 | [F]·m³/s | C15, E7 | §4 #51 |
| 4.8 | `material_volume_rate(u, cv, t) -> (surface_flux, volume_div)` (b = u: ∮u·n dA and ∫∇·u dV) | N53, D23 | m³/s | C15 (D23 check), B1 | §4 #52 |
| 4.9 | `rtt_ellipse_2d(F, dFdt, a, b, adot, bdot, c=(0.0, 0.0), cdot=(0.0, 0.0), t=0.0, n=256) -> (volume_term, surface_term, total)` — ellipse with semi-axes a, b (along x, y) centred at c, all moving at the given rates; b on the boundary = ċ + (ȧ cos s, ḃ sin s) at the point (c_x + a cos s, c_y + b sin s); 2-D areas | (3.35) in 2-D | [F]·m²/s | E7 RTT mode, C15 Fig. 3.18 | §8 (**signature adapted**: scalar sizes and rates instead of the curation's `a_fn, b_fn, center_fn` callables, so parity rows need no lambdas) |
| 4.10 | `swept_terms_sphere(R, Rdot, F, dFdt, t, dt) -> dict(T1, T2, T3, T4, exact, lhs)` — the four integrals of (3.32) for `GrowingSphere(R, Rdot)` (T1 = ∫_{V*}F, T2 = ∫_{V*}Δt ∂F/∂t, T3 = ∫_{ΔV}F, T4 = ∫_{ΔV}Δt ∂F/∂t), `exact` = ∫_{V*(t+Δt)}F(t+Δt), `lhs` = (exact − T1)/Δt | (3.32)–(3.33) | [F]·m³ | C15 slider, E7 limit view | §8 |

### C.5 `fluidpy/ch03_kinematics.py` (chapter module; re-exports C.0.6–C.0.9 and all of C.1–C.4)
| # | Callable (signature) | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 5.1 | `cross_section_average(u_fn, R, z=0.0) -> float` (u_fn(r, z); area average over a circular section by `quad`) | N03 | m/s | §3.1 | §4 #1 |
| 5.2 | `cylinder_flow(x, y, U, a, frame="body") -> (u, v)` (uniform stream U along +x + doublet: u = U[1 − a²(x² − y²)/r⁴], v = −2Ua²xy/r⁴; `frame="fluid"` subtracts U; NaN for r < a) | N04, Fig. 3.2 | m/s | C05, E3 | §4 #2 |
| 5.3 | `frame_acceleration_terms(x, y, U, a, U_frame, t=0.0, h=1e-5) -> dict(u, local, advective, total)` (each a (2,) array) — observer moving with velocity U_frame relative to the still fluid (0 = fluid frame, U = body frame): u_obs(X, t) = u_body(X + (U − U_frame)t, y) + (U_frame − U)e_x; local and advective by stencils; total independent of U_frame at t = 0 | N23, (3.9) | m/s, m/s² | C05, E3 | §8 (argument meaning fixed here) |
| 5.4 | `lagrangian_map_example(X, t, alpha) -> x` (x = X e^{αt}) | C01, D01 | m | C01, E2 mode 2 | §4 #7 |
| 5.5 | `streamline_slope(u, x, t) -> float` (v/u in 2-D) | (3.7) 2-D, N14 | – | C03 | §4 #15 |
| 5.6 | `example_3_1(t_prime, xi0=1.0, omega=1.0, n=200) -> dict(streamline, pathline, streakline, path_center, streak_center, radius, slope, field)` (curves (2, n) from the closed forms; `slope` = tan ωt′; `field(x, t)` the Ex. 3.1 callable) | Ex. 3.1, D04, D05 | m | C04, E1 parity | §4 #18 (keys fixed here) |
| 5.7 | `unsteady_flow_preset(name, x, y, t, **p) -> (u, v)` — `"ex31"` (u = ωξ_o cos ωt, v = ωξ_o sin ωt; p: omega, xi0), `"ex31_current"` (+U0 along x; p: U0), `"steady"` (uniform (ωξ_o cos β, ωξ_o sin β), β fixed; p: beta), `"rotating_strain"` (u = s[x cos 2Ωt + y sin 2Ωt, x sin 2Ωt − y cos 2Ωt]: pure strain whose axes turn at Ω; p: s, Omega; Ω = 0 is a steady hyperbolic flow) · `preset_field(name, **p) -> callable u(x, t)` (the same fields in the callable convention, for `streamline`/`pathline`/`streakline`) | E1 fields | m/s | C03, C04, E1 | §8 (`preset_field` **NEW**: the callable twin the ODE functions need) |
| 5.8 | `thermal_front(x, y, t, grad_K_per_m, heating_K_per_s=0.0, front_speed=0.0, width_m=None, T0=288.15) -> T` (cold to the north: linear T = T0 + Ht − G(y − ct) when `width_m` is None, else the front T = T0 + Ht − G w tanh((y − ct)/w), whose largest gradient is G) · `thermal_front_terms(x, y, t, u, v, grad_K_per_m, heating_K_per_s=0.0, front_speed=0.0, width_m=None) -> dict(T, dTdx, dTdy, local, advective, total)` (exact) | C02 worked number, E2 | K, K/s | C02, E2 | §8 (`width_m` and the `_terms` name fixed here) |
| 5.9 | `rigid_body_velocity(U, Omega, x) -> u` (U + Ω × x) | N27, D10 | m/s | C08, C11, E5 | §4 #27 |
| 5.10 | `potential_velocity(phi_expr, coords) -> (u_exprs, curl_expr)` (sympy ∇φ and its curl, = 0) | (3.17) | sympy | C10 (N29) | §4 #31 |
| 5.11 | `parallel_shear_kinematics(gamma) -> dict(G, S, R, omega3, spin, lam, axes)` (G = [[0, γ], [0, 0]]) | N34, N35 | 1/s | §3.5 head, C10, E5 | §4 #36 |
| 5.12 | `annular_sector_circulation(u_theta, r, dr, dtheta, u_r=None, **kw) -> float` (four legs by `quad`, counterclockwise: outer arc, inner arc back, two radial legs) | N41, D17 | m²/s | C13, E6 | §4 #43 (kind-name input added) |
| 5.13 | `example_3_2(h, r0, rdot) -> dict(direct, rtt_dblquad, rtt_cv)` | Ex. 3.2, N54 | m³/s | C15, E7 | §4 #53 |
| 5.14 | `leibniz_example(t, case="x2t") -> dict(interior, upper, lower, total, exact)` (`"x2t"`: F = x²t, a = t, b = t² — the case with a closed-form integral (t⁷ − t⁴)/3) | D21 check | – | C15 | §8 (return keys fixed here) |
| 5.15 | `rtt_field(name, dim=3) -> (F, dFdt)` — `"uniform"` F = 1 · `"ramp"` F = 1 + 0.5 x₁ · `"warming"` F = 1 + 0.5 x₁ + 0.3 t · `"carried"` F = 1 + 0.5 (x₁ − 0.3 t) (a ramp moving at 0.3 m/s along x₁; a box translating with it sees the two RTT terms cancel) (x₁ in m, t in s); `dim=1` returns F(x, t) of a position array, `dim=2/3` of `(d, N)` points | the fields of E7 and C15 | – | C15, E7 parity | **NEW** (so E7's three F choices and the notebook share one definition) |
| 5.16 | `flux_through_disc(u, center, normal, radius, t=0.0, nr=32, ntheta=64) -> float` (∫u·n dA over a planar disc, midpoint rule on `planar_disc` nodes) | N15 stream tube: equal flux through two sections | m³/s | C03 (N15) | **NEW** (curation named `core.V.flux_through_faces`, which is box-only) |

### C.6 `scripts/` (drawing helpers, no physics; `from scripts.ch03_drawings import …`)
`ch03_fig3_7_flow_lines.py` (`flow_lines_figure(t_prime, xi0=1, omega=1, ax=None)` — the three curves, port, common
tangent; `flow_lines_animation(frames)` — dye + one particle + turning streamline over a period) ·
`ch03_fig3_14_shear_elements.py` (`shear_elements_frames(gamma, times)` — aligned square ABCD and 45° square PQRS side
by side, side lengths and corner angles printed) · `ch03_fig3_17_leibniz.py` (`leibniz_strips_figure(case, t, dt, ax)` —
F(x, t), F(x, t + dt), the three strips) · `ch03_fig3_18_rtt.py` (`rtt_blob_figure(a, b, adot, bdot, cdot, F_name, dt)`
— deforming ellipse, b arrows, swept band coloured by sign(b·n), waterfall bars) · `ch03_drawings.py` (re-exports the
four plus `ring_arrows(ax, G, radius, parts=("total",))` — du arrows on a ring of neighbours, parts "total" / "strain"
/ "rotation", and `paddle(ax, xy, angle, size, color)`). Figures → `outputs/ch03/`.

**Count:** C.1 25 rows (≈ 40 callables) · C.2 4 rows (8) · C.3 10 rows · C.4 10 rows (≈ 16 incl. 4 shapes) · C.5 16
rows (≈ 19) · C.6 5 scripts. **Beyond analysis §4 and curation §8 (NEW):** `volume_integral` (4.4), `rtt_field` (5.15),
`flux_through_disc` (5.16), `preset_field` (in 5.7); **signatures adapted from the curation:** `rtt_ellipse_2d` (scalar
sizes and rates), kind-name inputs for the `u_theta` arguments (3.9, 3.10, 5.12), return keys of `lagrangian_to_eulerian`,
`example_3_1`, `leibniz_example`, `thermal_front_terms`, the `method` switch of `strain_ellipse_axes`, the constructor
arguments of the control-volume shapes. Everything else is analysis §4 or curation §8 as written.

**Contract details the notebook and explainers rely on (please keep):** `principal_strain_rates` returns λ ascending
(shear γ = 1 → (−0.5, 0.5), eigenvector of +0.5 = (1, 1)/√2); `strain_ellipse_axes` returns semi-axes **descending**;
`material_line_rotation_rate` is counterclockwise-positive (shear γ = 1, θ = 30° → −0.25); `vorticity_from_gradient`
always returns 3 components; `streakline` columns follow the order of `t_release`; `thermal_front_terms` uses
∂T/∂y = −G (cold north) and "local" = ∂T/∂t at the point; `frame_acceleration_terms(0, 1.5, 1, 1, 0)["local"][1]` =
−0.592593 m/s² (the C05 worked number); `rtt_field("warming")` = 1 + 0.5x₁ + 0.3t exactly (E7 parity numbers depend on
it); `leibniz_terms` returns `lower` = ȧF(a) (the term that is **subtracted**).

**Implementation cross-check (probe of the implementer's WIP, 2026-09-23).** `fluidpy/core/kinematics.py`,
`coords.py`, `vortices.py`, `transport.py` and `fluidpy/ch03_kinematics.py` already carry every Part C name with
compatible signatures and the return keys fixed above (`lagrangian_to_eulerian` keys, `example_3_1`
`path_center`/`streak_center`/`slope`, `thermal_front_terms` positional order, `frame_acceleration_terms`,
`strain_ellipse_axes(method=…)`, `measured_strain_rates` keys `stretch`/`closing`/`area` (+ `*_formula`, `spin`),
`annular_sector_circulation` kind names, `leibniz_example` keys, `example_3_2` keys). **One gap:**
`rtt_field` knows `uniform`, `ramp`, `warming` but not **`carried`** (F = 1 + 0.5(x₁ − 0.3t), needed by E7's preset
"riding with a moving pattern" and its parity row) — ≈ 3 lines to add. **Extra callables the builders may use:**
`ch03.pipe_profile(r, z, U_mean, R, n0, L_e)` (use it for the §3.1 developing-profile figure instead of the ad-hoc
family of A.1 #7), `ch03.cylinder_velocity_field`, `ch03.cylinder_streamfunction` (streamlines of Fig. 3.2 as ψ
contours), `ch03.perpendicular_pair_rotation_rate(G, theta)` (D14's pair average in one call), `ch03.unit_vectors_polar`,
`ch03.polar_from_cartesian`, `ch03.cartesian_components`, `ch03.MovingEllipse2D` (the E7 ellipse as a
`ControlVolume`), `ch03.swept_terms(F, dFdt, cv, t, dt)` (the four (3.32) terms for any shape), and the extra E1 preset
`"steady_vortex"` (u = Ω(−y, x): curved streamlines on which the three lines coincide — a better "steady" preset than
the uniform one). Defaults differ harmlessly from the storyboard (`unsteady_flow_preset(..., s=0.5, Omega=0.5, U0=0.5,
beta=0.0)`): every notebook call and parity row above passes these arguments explicitly.

---

## Part A — notebook storyboard (`notebooks/build_ch03.py` → `notebooks/ch03_kinematics.ipynb`)

**Section check (book order; one `nb.section` per book section; cell numbers are estimates ± 4).**
- §3.1 → no A item (curation decision 1): N01 N02 N03 N05 N06 (tagged → C01), pointer to N04 (taught in C05) — cells 7–21
- §3.2 → C01 (N07 N08 · D01), C02 (N09–N12 · D02 · **E2**) — cells 22–66
- §3.3 → C03 (N14 N15 · D03), C04 (N13 N16–N19 · D04 D05 · **E1**), C05 (N04 N20–N24 · D06 · **E3**) — cells 67–140
- §3.4 → R01 R02 R03, C06 (N25 N26 · D07), C07 (D08), C08 (N27 · D09 D10), C09 (D11 · **E4**), R04 R05 R06 R07,
  C10 (N28 N29 · D12 D14 D13), C11 (N30 N33 · D15), C12 (N31 N32 · D16 · **E5**) — cells 141–262
- §3.5 → N34, N35 (head, with the Fig. 3.14 animation), C13 (N36–N41 · D17), C14 (N42–N44 · D18 D19 D20 · **E6**) — cells 263–318
- §3.6 → C15 (N45–N55 · D21 D22 D23 D24 · **E7**) — cells 319–362
- end → S01, S02, summary — cells 363–366

Every CORE block below has the nine parts in order (question · idea · primers · maths with its D rows · tiny example ·
code + "What does the code above do?" · from-scratch check where §7 of the curation asks · ≥ 1 visual · notes + "What
would change if…"). *code:* is the intent (exact `ch03.` calls and arguments; every line gets a novice comment in the
builder); *expect:* the numbers the executed cell must print (computed by hand here — the builder re-checks them
against the implemented `fluidpy`); *explain:* the numbered "What does the code above do?" list. Markdown drafts are in
our words; every equation is written out in LaTeX next to its number (convention 11). Derivation cells are copied from
Part F word for word.

### A.0 Front matter
1. `nb.title(big_idea="Kinematics describes motion without asking what causes it. A fluid can be described by following
   each particle (Lagrangian) or by watching fixed points (Eulerian); the material derivative
   $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$ *(Eq. 3.5)* is the bridge. Three kinds of flow
   line show what the flow does and differ whenever it is unsteady. Near any point the motion splits into a
   deformation (the strain-rate tensor S: stretching, shearing, swelling) and a rigid spin at half the vorticity. The
   chapter ends with the Reynolds transport theorem, the rule for differentiating an integral over a moving volume —
   the tool Chapter 4 uses to write every conservation law.", roadmap=["§3.1 steady vs unsteady, 1-/2-/3-D flows,
   cylindrical and spherical coordinates", "§3.2 Lagrangian and Eulerian descriptions; the material derivative D/Dt",
   "§3.3 streamlines, path lines, streak lines; the acceleration is the same for every steadily moving observer",
   "§3.4 relative motion near a point: stretching, shearing, volume change, spin = ½ vorticity, principal strain axes",
   "§3.5 shear flow, solid-body rotation, the irrotational vortex, Rankine and Gaussian vortices", "§3.6 Leibniz's rule
   and the Reynolds transport theorem"], prerequisites=["partial derivatives, the chain rule and first-order Taylor
   expansion (Ch. 1 primers P25, P49, P26)", "index notation, the velocity gradient G, the split G = S + ½R and the
   vorticity ω = ∇×u (Ch. 2 C01, C10–C12)", "eigenvalues and principal axes of a symmetric tensor (Ch. 2 C13)", "Gauss'
   and Stokes' theorems, circulation (Ch. 2 C14, C16)", "scipy.integrate.solve_ivp and scipy.linalg.expm (Ch. 1 P31,
   Ch. 2 P79)"])`
2. `nb.explainer_index([("flow_lines_unsteady", "Three lines through one point — why do they disagree?", "C03 C04:
   streamline, path line and streak line in an unsteady flow (Ex. 3.1)"), ("material_derivative_probe", "Why does the
   station warm while the air does not?", "C01 C02: DF/Dt = ∂F/∂t + u·∇F with a probe and a float"),
   ("galilean_frames_cylinder", "Steady or not — does the acceleration care?", "C05: the local/advective split moves,
   the total does not (3.9)"), ("fluid_element_deformation", "What does each number in S measure?", "C06–C09:
   stretching, shearing and swelling measured on a moving element"), ("spin_and_principal_axes", "Can a straight flow
   make a fluid element spin?", "C10–C12: spin = ½ω for any pair, du = S·dx + ½ω×dx, principal axes"),
   ("vortex_paddle_wheels", "Going round in circles ≠ spinning", "C13 C14: solid-body, irrotational, Rankine and
   Gaussian vortices"), ("reynolds_transport_cv", "What changes inside a moving box?", "C15: Leibniz and the Reynolds
   transport theorem as a budget")])` (the index text is plain; the explainers themselves show every equation in full).
3. `nb.setup()`, then `nb.code` (no CORE yet) — *code:* `import numpy as np; import sympy as sp; import
   matplotlib.pyplot as plt; import plotly.graph_objects as go` · `from fluidpy import ch03_kinematics as ch03` ·
   `from fluidpy.core.interact import slider_figure` · `from fluidpy.core.anim import animate` · `from fluidpy.core.style
   import COLORS, savefig` · `from tools.convergence import observed_order`. *expect:* no output.
   *explain:* 1. numerical, symbolic and plotting libraries; 2. the chapter module (every function of this notebook
   lives there and is tested in `tests/test_ch03.py`); 3. the house helpers for sliders, animations and colours.
4. `nb.md` — "**Notation and colours used in this notebook.**" A table: | symbol | meaning | unit | — x = (x₁, x₂, x₃) or
   (x, y, z) position [m]; t time [s]; u = (u₁, u₂, u₃) = (u, v, w) velocity [m/s]; F any field (temperature, a velocity
   component…); G_ij = ∂u_i/∂x_j velocity gradient [1/s] (row = component); S strain rate, R = G − Gᵀ rotation tensor
   (the book's, no ½), ω = ∇×u vorticity, **spin of a fluid element = ½ω**; γ = du₁/dx₂ shear rate (= 2S₁₂; Ch. 2's
   Γ was S₁₂); Γ circulation [m²/s]; θ: polar angle in the plane (§3.5), polar angle from +z (spherical), cone
   half-angle (Ex. 3.2). Colours: streamline teal · path line orange · streak line rose · ∂/∂t blue · advective amber
   · D/Dt purple · strain teal · rotation orange · stretching blue / compressing rose · volume term blue, surface term
   orange.
5. `nb.md` — "**Where this chapter is used later**" (two columns): $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf
   u\cdot\nabla F$ *(3.5)* → every conservation law of Ch. 4, the vorticity equation (Ch. 5), the potential-vorticity
   equation (Ch. 13) · path lines $d\mathbf r/dt=\mathbf u(\mathbf r,t)$ *(3.8)* → particle orbits under waves (Ch. 7),
   Lagrangian statistics (Ch. 12), trajectories (Ch. 13) · $du_i=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf
   x)_i$ *(3.19)* → the Newtonian stress law (Ch. 4 §4.5), vortex stretching (Ch. 5) · circulation
   $\Gamma=\oint\mathbf u\cdot d\mathbf s$ *(3.18)* → Kelvin's theorem (Ch. 5), lift (Ch. 6, 14) · $\omega'_z=\omega_z-2\Omega$
   → absolute vs relative vorticity $f+\zeta$ (Ch. 4 §4.7, Ch. 13) · Reynolds transport theorem
   $\frac{d}{dt}\int_{V^*}F\,dV=\int_{V^*}\frac{\partial F}{\partial t}dV+\int_{A^*}F\,\mathbf b\cdot\mathbf n\,dA$ *(3.35)* →
   mass, momentum and energy equations (Ch. 4 §4.2–4.4), layer budgets (Ch. 13). Climate hook: the advective term
   $\mathbf u\cdot\nabla T$ is "warm/cold advection" on every weather map.

### A.1 §3.1 Introduction and Coordinate Systems — no A item (paragraphs tagged → C01; N04 → C05)
1. `nb.section("3.1", "Introduction and Coordinate Systems", intro="**What is this section about?** The words used to
   describe *any* flow before we ask what drives it: steady or unsteady, one-, two- or three-dimensional, and the
   coordinate systems whose velocity components later chapters use (Cartesian, plane polar, cylindrical,
   spherical).")`
2. `nb.note` — **N01 [C]** "**Kinematics** describes how a fluid moves — positions, velocities, accelerations,
   stretching and spinning — without asking which forces cause it. Forces enter in Ch. 4 (§4.4 onward, Newton's second
   law for a fluid)."
3. `nb.note` — **N02 [B]** "A flow is **steady** when nothing measured at a fixed point changes with time:
   $\partial(\cdot)/\partial t=0$ for every field. Otherwise it is **unsteady**. Careful: steadiness depends on the
   observer — flow past a parked car is steady, the same air seen from the pavement as a car drives by is not. C05 turns
   this into a theorem." equation `\frac{\partial(\cdot)}{\partial t}=0\quad\text{(steady)}`, no ref.
4. `nb.primer("scipy.integrate.quad and dblquad", "Adaptive numerical integration: `quad(f, a, b)` returns
   $\int_a^b f\,dx$ and an error estimate; `dblquad(f, a, b, gfun, hfun)` does a double integral, the inner variable
   first (its limits may depend on the outer one). We use them for section averages here, for the cone of Ex. 3.2 and
   for exact checks of the Reynolds transport theorem in §3.6.", code="from scipy.integrate import quad, dblquad
   \nval, err = quad(lambda r: 2*r, 0.0, 1.0)            # ∫₀¹ 2r dr = 1\nprint(val, err)                                     # 1.0 and a tiny error
   estimate\narea, _ = dblquad(lambda r, th: r, 0.0, 2*np.pi, 0.0, 1.0)   # ∫∫ r dr dθ over the unit disc (inner variable
   r first)\nprint(area, np.pi)                                  # 3.14159… twice")` — *expect:* `1.0 1.1e-14`, `3.14159 3.14159`.
5. `nb.note` — **N03 [B]** "**1-D, 2-D and 3-D flows.** A flow is 3-D when it depends on all three coordinates, plane
   (2-D) when one coordinate can be dropped, and — as an engineering approximation — 1-D when we keep only the average
   over each cross-section of a pipe or channel (an axisymmetric pipe flow is still 3-D in the book's sense). The 1-D
   description keeps the average and throws away the profile:" equation
   `\bar u(z)=\frac{1}{A}\int_A u\,dA` (our statement of Fig. 3.1d, no number). "Number: the parabolic (Poiseuille)
   profile $u=U(1-r^2/R^2)$ averages to exactly $U/2$."
6. `nb.code` — *code:* `U, R = 2.0, 0.05` (centre speed 2 m/s, pipe radius 5 cm) · `u_pois = lambda r, z: U*(1 - (r/R)**2)`
   · `ubar = ch03.cross_section_average(u_pois, R)` · `print(f"section average = {ubar:.4f} m/s (U/2 = {U/2})")` ·
   from scratch on the same line: `r = np.linspace(0, R, 2001); mine = np.trapezoid(u_pois(r, 0)*2*np.pi*r, r)/(np.pi*R**2)`;
   `assert np.allclose(mine, ubar, rtol=1e-6)`. *expect:* `section average = 1.0000 m/s (U/2 = 1.0)`. *explain:* 1. the
   profile as a function of radius r and axial position z; 2. `cross_section_average` integrates u · 2πr dr with
   `quad` and divides by the area πR²; 3. the trapezoid version with rings of area 2πr dr agrees (curation §7: the §3.1
   from-scratch moment).
7. `nb.figure` — our Fig. 3.1c/d: a developing pipe flow family of our own, $u=U_m(z)\,(1-(r/R)^{n(z)})$ with n falling
   from 10 (flat entrance profile) to 2 (parabola) and $U_m=\bar u\,(n+2)/n$ so the average stays 1 m/s; left panel:
   profiles at z = 0, 0.5, 2 m (teal shades); right panel: ū(z) = 1 m/s (the whole 1-D description) with the centre
   speed U_m(z) dashed. *see:* "three profiles that sharpen downstream, and one flat line." *read:* "the 1-D model sees
   only the right panel: the same average at every z while the shape changes completely — mass is conserved, detail is
   lost." *change:* "…the pipe narrowed to half the radius: ū would jump ×4 (same volume flow through a quarter of the
   area), the shape question stays open."
8. `nb.md` — pointer sentence: "Fig. 3.2's two views of one cylinder — fixed in a stream, or towed through still water —
   are the same flow seen by two observers; we take them up with the Galilean transformation in C05 (§3.3)."
9. `nb.primer("cylindrical and spherical unit vectors", "At every point P the curvilinear coordinates carry their own
   right-handed unit vectors, and they turn as P moves. Cylindrical $(R,\varphi,z)$: $\mathbf e_R=(\cos\varphi,\sin\varphi,0)$,
   $\mathbf e_\varphi=(-\sin\varphi,\cos\varphi,0)$, $\mathbf e_z$. Spherical $(r,\theta,\varphi)$ with θ measured from +z:
   $\mathbf e_r=(\sin\theta\cos\varphi,\sin\theta\sin\varphi,\cos\theta)$, $\mathbf e_\theta=(\cos\theta\cos\varphi,
   \cos\theta\sin\varphi,-\sin\theta)$, $\mathbf e_\varphi=(-\sin\varphi,\cos\varphi,0)$. A velocity component is the dot
   product of **u** with one of them (a projection, Ch. 2 P65).", code="phi = np.pi/4\ne_R = np.array([np.cos(phi),
   np.sin(phi), 0.0]); e_phi = np.array([-np.sin(phi), np.cos(phi), 0.0])\nprint(e_R @ e_phi, np.cross(e_R, e_phi))
   # 0.0 and (0, 0, 1): perpendicular, and e_R × e_φ = e_z (right-handed)")` — *expect:* `0.0 [0. 0. 1.]`.
10. `nb.note` — **N05 [B]** "**Coordinate systems (Fig. 3.3).** Plane: (x, y) = (x₁, x₂) or polar (r, θ) (Ch. 2 Ex. 2.1).
    Cylindrical (R, φ, z) with velocity components (u_R, u_φ, u_z); spherical (r, θ, φ), θ the angle from +z, φ the
    azimuth, components (u_r, u_θ, u_φ). ⚠️ The same letter θ means the plane polar angle in §3.5 but the angle from
    the z-axis in spherical coordinates." equation `R=\sqrt{x^2+y^2},\quad \varphi=\tan^{-1}(y/x);\qquad
    r=\sqrt{x^2+y^2+z^2},\quad \theta=\tan^{-1}\!\big(\sqrt{x^2+y^2}/z\big)` (no number; the book states these in its
    exercises).
11. `nb.code` — *code:* `P = (1.0, 1.0, 1.0)` · `print(ch03.cylindrical_from_cartesian(*P))` ·
    `print(ch03.spherical_from_cartesian(*P))` · `u = np.array([1.0, 0.0, 0.0])` (a 1 m/s wind along x at P) ·
    `print(ch03.velocity_components(u, np.array(P), "cylindrical"))` · `print(ch03.velocity_components(u, np.array(P),
    "spherical"))` · `print(ch03.unit_vectors_spherical(*ch03.spherical_from_cartesian(*P)[1:]))`. *expect:* (R, φ, z) =
    (1.4142, 0.7854, 1.0) · (r, θ, φ) = (1.7321, 0.9553, 0.7854) (θ = 54.74°) · (u_R, u_φ, u_z) = (0.7071, −0.7071, 0) ·
    (u_r, u_θ, u_φ) = (0.5774, 0.4082, −0.7071) — both have length 1. *explain:* 1. position P in the two coordinate
    systems; 2. the same wind vector projected on the local unit vectors; 3. the squares add to |u|² = 1 in every system
    (a projection onto an orthonormal basis keeps length).
12. `nb.plotly` — 3-D figure: point P = (1, 1, 1) with the Cartesian (grey), cylindrical (teal) and spherical (orange)
    unit vectors as `go.Cone` arrows (≤ 9 cones), dashed guide lines to the axes, a dropdown "Cartesian / cylindrical /
    spherical / all"; `height=460`.
13. `nb.md` — **What you see / How to read it / What would change if…** "Three triads of arrows at one point. The
    Cartesian triad is the same everywhere; the teal and orange triads are tied to P: e_R points away from the z-axis,
    e_r away from the origin. Move P to (−1, 1, 1) (edit P in the code cell) and watch e_R and e_φ swing by 90° while the
    Cartesian arrows stay put — this is why derivatives of curvilinear components need extra terms (Appendix B, N06)."
14. `nb.note` — **N06 [C]** "Appendix B of the book lists ∇, ∇², ∇·u and (u·∇)u in cylindrical and spherical coordinates;
    we meet the first of them in C13 (the vorticity in polar coordinates, $\omega_z=\frac1r\frac{\partial}{\partial r}(ru_\theta)-\frac1r\frac{\partial u_r}{\partial\theta}$
    *(Eq. 3.23)*, derived there by hand) and the rest in Ch. 4's Navier–Stokes equations."

### A.2 §3.2 Particle and Field Descriptions of Fluid Motion — C01 (+N07, N08 · D01), C02 (+N09–N12 · D02 · E2)
1. `nb.section("3.2", "Particle and Field Descriptions of Fluid Motion", intro="**What is this section about?** Two ways
   to describe one flow — ride along with every particle, or stand still and watch every point — and the single
   formula, the material derivative, that turns one into the other.")`

**C01 — Lagrangian and Eulerian descriptions and the bridge (3.2)**
2. `nb.core("C01", "Following a particle or watching a point — and the bridge between them (3.2)", question="A float
   drifts down a river past a thermometer bolted to a bridge pier. Both report the temperature of the same water. How
   do the two sets of numbers fit together?")`
3. `nb.md` — **The problem in plain words:** "Oceanographers throw Argo floats into the sea and read where each one goes
   (it *is* a water parcel); weather services keep thermometers at fixed stations. Numerical models store fields on a
   fixed grid. We need to translate freely between 'what happens to this parcel' and 'what happens at this place' —
   the whole of fluid dynamics is written in the second language but Newton's law is about the first."
4. `nb.md` — **The idea:**
   ```
   Lagrangian (follow the particle)            Eulerian (watch the point)
   label r_o at time t_o  →  r(t; r_o, t_o)    field F(x, t): four independent variables x, y, z, t
   "what does THIS parcel do?"                 "what happens HERE?"
            the bridge:  F[r(t; r_o, t_o), t] = F(x, t)   when   x = r(t; r_o, t_o)     (3.2)
   ```
   "**A label is not a variable**: r_o only names the particle (where it was at t_o), just as a float's serial number
   does. The bridge says: at time t, the field at x *is* whatever the particle that happens to be at x carries."
5. `nb.primer("functions of time with parameters", "A trajectory $x=X e^{\alpha t}$ is a function of time t alone once
   the label X is chosen; X is a *parameter* that picks one curve out of a family. Differentiating 'with the label
   held fixed' means following one member of the family. Here X is where the particle was at t = 0.", code="alpha = 0.5
   # stretching rate α [1/s]\nx_of = lambda t, X: X*np.exp(alpha*t)          # one curve per label X\nprint([round(x_of(1.0,
   X), 3) for X in (1.0, 2.0)])   # [1.649, 3.297]: two particles at t = 1 s")` — *expect:* `[1.649, 3.297]`.
6. `nb.note` — **N07 [B]** "**Lagrangian description.** Each particle is labelled by its position r_o at a reference time
   t_o; its later position is $\mathbf r(t;\mathbf r_o,t_o)$ (Fig. 3.4). Any property it carries is written
   $F[\mathbf r(t;\mathbf r_o,t_o),t]$. The labels are fixed numbers for a given particle, not coordinates." equation
   `\mathbf r=\mathbf r(t;\mathbf r_o,t_o)`, no ref.
7. `nb.note` — **N08 [B]** "**Velocity and acceleration of a particle** are ordinary time derivatives along its own path —
   exactly single-particle mechanics:" equation `\mathbf u=d\mathbf r(t;\mathbf r_o,t_o)/dt\quad\text{and}\quad\mathbf
   a=d^2\mathbf r(t;\mathbf r_o,t_o)/dt^2`, ref "3.1". "Number (the running example): $x=Xe^{\alpha t}$ with X = 2 m, α = 0.5
   s⁻¹, t = 1 s gives u = αXe^{αt} = 1.649 m/s and a = α²Xe^{αt} = 0.824 m/s²."
8. `nb.md` — **The maths — the Eulerian description and the bridge.** "Watching fixed points, a property is a field
   $F(\mathbf x,t)$ of four independent variables. The two descriptions must agree when the particle position and the
   field point coincide, in the same coordinates and on a common clock:
   $$F[\mathbf r(t;\mathbf r_o,t_o),t]=F(\mathbf x,t)\quad\text{when}\quad\mathbf x=\mathbf r(t;\mathbf r_o,t_o).\qquad(3.2)$$
   Read it as a recipe: to get the Eulerian field at (x, t), find the label of the particle that sits at x at time t,
   then read that particle's value. The derivation below does this for a concrete flow."
9. `nb.primer("inverse functions and sympy solve", "To use (3.2) we must answer 'which particle is at x at time t?' —
   that is, solve $x=r(t;X)$ for the label X. `sympy.solve(equation, unknown)` does it symbolically; it returns a list of
   solutions.", code="X, x, t, alpha = sp.symbols('X x t alpha', positive=True)\nprint(sp.solve(sp.Eq(x, X*sp.exp(alpha*t)),
   X))   # [x*exp(-alpha*t)]: the label of the particle now at x")` — *expect:* `[x*exp(-alpha*t)]`.
10. `nb.derivation("D01", …)` — Part F D01 (6 steps), ref "".
11. `nb.worked_example("the stretching flow x = X e^{αt} with X = 2 m, α = 0.5 s⁻¹, t = 1 s", "1. Position:
    $x=2e^{0.5}=2\times1.6487=3.297$ m. 2. Lagrangian velocity from (3.1), $u=d\mathbf r/dt$ with the label fixed:
    $u=\alpha Xe^{\alpha t}=0.5\times3.297=1.649$ m/s. 3. Eulerian velocity field from D01: $u(x,t)=\alpha x=0.5\times3.297=1.649$
    m/s — the same number, as (3.2), $F[\mathbf r(t;\mathbf r_o,t_o),t]=F(\mathbf x,t)$ at $\mathbf x=\mathbf r$, demands.
    4. Acceleration: $a=\alpha^2x=0.25\times3.297=0.824$ m/s². 5. A fixed probe at x = 3.297 m reads 1.649 m/s for ever (the
    field is steady) while every particle passing it speeds up.")`
12. `nb.code` — *code:* `X, alpha, t = 2.0, 0.5, 1.0` · `x = ch03.lagrangian_map_example(X, t, alpha)` ·
    `u, a = ch03.lagrangian_velocity_acceleration(lambda s: ch03.lagrangian_map_example(X, s, alpha), t)` ·
    `print(f"x = {x:.4f} m, u = {u:.4f} m/s, alpha*x = {alpha*x:.4f}, a = {a:.4f} m/s^2")` · sympy twin:
    `Xs, xs, ts, al = sp.symbols('X x t alpha', positive=True)`; `res = ch03.lagrangian_to_eulerian([Xs*sp.exp(al*ts)],
    [Xs], [xs], ts)`; `print(res["label_of_x"], res["u"], res["a_lagrangian"], sp.simplify(res["Du_Dt"][0] -
    res["a_lagrangian"][0]))`. *expect:* `x = 3.2974 m, u = 1.6487 m/s, alpha*x = 1.6487, a = 0.8244 m/s^2` ·
    `[x*exp(-alpha*t)] [alpha*x] [alpha**2*x] 0` (list/Matrix formatting as sympy prints it). *explain:* 1. the map and
    its derivatives at fixed label (3.1) by central differences; 2. the symbolic route solves for the label, substitutes
    into dr/dt and d²r/dt², and 3. checks that the material derivative of the Eulerian u — the formula C02 derives — gives
    the same acceleration (the last printed 0).
13. `nb.check_agree` — **from scratch:** `h = 1e-4`; `u_fd = (ch03.lagrangian_map_example(X, t+h, alpha) -
    ch03.lagrangian_map_example(X, t-h, alpha))/(2*h)` (central difference: slope of the particle's path) ·
    `assert np.allclose(u_fd, alpha*x, rtol=1e-7)` (the Lagrangian slope equals the Eulerian field αx at the particle's
    position) · `a_fd = (x(t+h) − 2x(t) + x(t−h))/h**2`; `assert np.allclose(a_fd, alpha**2*x, rtol=1e-5)`.
14. `nb.figure` — two panels: (a) space–time diagram: path lines x(t) = X e^{αt} of six labelled particles X = 0.25…1.5 m
    for 0 ≤ t ≤ 3 s (orange, labels at the start), a vertical grey line "fixed probe at x = 2 m", and the tangent slope
    at one crossing annotated "u = αx = 1 m/s"; (b) the Eulerian snapshot u(x) = αx at t = 0 and t = 2 s (the same line,
    teal, dashed for t = 2 s) with the particles' positions at each time as dots. *see:* "curves fanning out in (a); one
    straight line in (b) that does not move." *read:* "every particle crossing the probe line crosses it with the same
    slope 1 m/s (the field is steady), but each curve bends upward (each particle accelerates). The two panels are the
    two descriptions of one motion." *change:* "…α doubled: the fan opens twice as fast, the line in (b) doubles its
    slope — and still does not move in time."
15. `nb.md` — "> ⚠️ **Common confusion:** 'a steady velocity field means the fluid does not accelerate.' Here
    $u=\alpha x$ does not depend on t at any fixed point, yet every particle speeds up at $a=\alpha^2x$. Steadiness is a
    statement about points; acceleration is about particles. C02 shows the missing piece is the advective term."
16. `nb.md` — **What would change if…** "…the map were $x=X+Ut$ (uniform drift)? Then $u=U$ everywhere, a = 0, and both
    descriptions are trivially the same. The interesting case is when the particle moves *through* a field that varies
    in space — then a thermometer on the float and one on the pier disagree. How fast does the float's reading change?
    That is C02."

**C02 — The material derivative (3.3)–(3.6)**
17. `nb.core("C02", "The material derivative: the rate a moving parcel feels (3.4)–(3.5)", question="A weather station
    reports warming of 0.36 K per hour, yet the air blowing past it has kept exactly its own temperature. How can both
    be true?")`
18. `nb.md` — **The problem in plain words:** "A southerly wind carries warm air north over a station. The station's
    thermometer rises; a thermometer tied to a balloon drifting with the same air does not change at all. Forecasters
    call the station's warming *warm advection*. We want the formula that separates 'the field changes here' from 'the
    parcel moves to where the field is different' — the material derivative."
19. `nb.md` — **The idea:**
   ```
   thermometer on the pier      reads  ∂T/∂t          (x fixed)
   thermometer on the balloon   reads  DT/Dt          (label fixed)
   the difference                      u·∇T           (the balloon moves through a gradient)
              DT/Dt  =  ∂T/∂t  +  u·∇T         → blue + amber = purple in every picture below
   ```
20. `nb.primer("multivariable chain rule along a path", "If $f(t)=F(x(t),y(t),z(t),t)$, then
    $\frac{df}{dt}=\frac{\partial F}{\partial x}\frac{dx}{dt}+\frac{\partial F}{\partial y}\frac{dy}{dt}+\frac{\partial F}{\partial z}\frac{dz}{dt}+\frac{\partial F}{\partial t}$:
    each argument that changes contributes (sensitivity to it) × (its rate). Time appears twice — through the moving
    position and through the explicit clock. This extends the two-variable rule of Ch. 1 (P49).", code="t = sp.symbols('t')
    \nx = sp.cos(t); F = lambda x, t: x**2*t                    # a path and a field\ndirect = sp.diff(F(x, t), t)
    # differentiate after substituting the path\nX, T = sp.symbols('X T')\nchain = (sp.diff(F(X, T), X)*sp.diff(x, t) +
    sp.diff(F(X, T), T)).subs({X: x, T: t})\nprint(sp.simplify(direct - chain))                   # 0: the chain rule
    holds")` — *expect:* `0`.
21. `nb.derivation("D02", …)` — Part F D02 (9 steps), ref "3.5".
22. `nb.note` — **N09 [B]** "Step 3 of the derivation is the book's (3.3) written out — the chain rule along a particle's
    path, with the four arguments of F each contributing:" equation
    `\frac{d}{dt}F[\mathbf r(t;\mathbf r_o,t_o),t]=\frac{\partial F}{\partial r_1}\frac{dr_1}{dt}+\frac{\partial F}{\partial r_2}\frac{dr_2}{dt}+\frac{\partial F}{\partial r_3}\frac{dr_3}{dt}+\frac{\partial F}{\partial t}=\frac{d}{dt}F(\mathbf x,t)\quad\text{when}\quad\mathbf x=\mathbf r(t;\mathbf r_o,t_o)`,
    ref "3.3". And (3.4) is steps 4–7 (second `nb.note`): equation
    `\frac{d}{dt}F[\mathbf r(t;\mathbf r_o,t_o),t]=\frac{\partial F}{\partial x_1}u_1+\frac{\partial F}{\partial x_2}u_2+\frac{\partial F}{\partial x_3}u_3+\frac{\partial F}{\partial t}=(\nabla F)\cdot\mathbf u+\frac{\partial F}{\partial t}\equiv\frac{D}{Dt}F(\mathbf x,t)`, ref "3.4".
    "Names: material, substantial or particle derivative."
23. `nb.note` — **N11 [B]** "In vector and index notation (the forms used for the rest of the book):" equation
    `\frac{DF}{Dt}\equiv\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F,\quad\text{or}\quad\frac{DF}{Dt}\equiv\frac{\partial F}{\partial t}+u_i\frac{\partial F}{\partial x_i}`, ref "3.5".
    Followed by `nb.code` — *code:* `print(ch03.expand_indices_str("u_i dF/dx_i"))` (the spelling the ch02 parser
    accepts; builder checks once). *expect:* the three-term sum `u_1*dF/dx_1 + u_2*dF/dx_2 + u_3*dF/dx_3` (as the parser
    prints it). *explain:* "the repeated index i is the dot product u·∇F (Ch. 2 C01)."
24. `nb.note` — **N10 [B]** "**Local and advective parts.** $\partial F/\partial t$ (blue) is the *local* or unsteady
    rate — what a fixed probe sees; it vanishes when F does not depend on time. $\mathbf u\cdot\nabla F$ (amber) is the
    *advective* rate — the change a particle meets because it moves to where F is different; it vanishes when F is
    uniform, when u = 0, or when u ⟂ ∇F (wind blowing along the isotherms). The book says *advection* for transport by
    the flow and keeps *convection* for heat carried by fluid motion."
25. `nb.note` — **N12 [B]** "Because $\mathbf u\cdot\nabla F=|\mathbf u|\,(\mathbf e_u\cdot\nabla F)$ with
    $\mathbf e_u=\mathbf u/|\mathbf u|$, and $\mathbf e_u\cdot\nabla F=\partial F/\partial s$ is the directional
    derivative along the path (Ch. 2 P75), the material derivative can be written with the arc length s along the
    particle's path:" equation `\frac{DF}{Dt}=\frac{\partial F}{\partial t}+|\mathbf u|\frac{\partial F}{\partial s}`, ref
    "3.6". "> ⚠️ **Book typo:** the printed (3.6) ends with $|\mathbf u|\,\partial/\partial s$ — the F is missing (both
    terms must be rates of change of F). The form is undefined at a stagnation point (u = 0, no direction)."
26. `nb.md` — "> ⚠️ **Common confusion:** '$\partial T/\partial t$ is how fast the air warms.' It is how fast a *fixed
    thermometer* warms. The air's own rate is $DT/Dt$. They differ by the advective term
    $\mathbf u\cdot\nabla T$ — often the largest term on a weather map."
27. `nb.worked_example("warm advection over a station", "Temperature falls northward by 1 K per 100 km:
    $\partial T/\partial y=-1/10^5=-1\times10^{-5}$ K/m. A southerly wind v = 10 m/s blows north. 1. Advective term:
    $\mathbf u\cdot\nabla T=v\,\partial T/\partial y=10\times(-10^{-5})=-1\times10^{-4}$ K/s. 2. The air keeps its
    temperature (no heating): $DT/Dt=0$. 3. From (3.5), $\frac{DT}{Dt}=\frac{\partial T}{\partial t}+\mathbf u\cdot\nabla T$:
    $\partial T/\partial t=0-(-10^{-4})=+1\times10^{-4}$ K/s. 4. Per hour: $10^{-4}\times3600=0.36$ K/h of warming at the
    station, although no parcel warmed at all.")`
28. `nb.code` — *code:* `G, v = 1e-5, 10.0` (gradient [K/m], wind [m/s]) · linear front carried by the wind
    (`front_speed = v`): `terms = ch03.thermal_front_terms(0.0, 0.0, 0.0, 0.0, v, G, heating_K_per_s=0.0,
    front_speed=v)`; `print({k: f"{terms[k]:+.2e}" for k in ("local", "advective", "total")})` ·
    `print(f"station warming: {terms['local']*3600:.2f} K/h")` · the general stencils on the same field:
    `F = lambda x, t: ch03.thermal_front(x[0], x[1], t, G, 0.0, v)`; `u = lambda x, t: np.array([0.0, v])`;
    `print(ch03.material_derivative_terms(F, u, np.array([0.0, 0.0]), 0.0, h=10.0))` (h = 10 m on a 100 km scale) ·
    streamwise form: `print(ch03.streamwise_derivative(F, u, np.array([0.0, 0.0]), 0.0, h=10.0))` · symbolic:
    `x_, y_, t_ = sp.symbols('x y t')`; `print(sp.simplify(ch03.material_derivative_sym(288 - sp.Rational(1, 100000)*(y_ -
    10*t_), [0, 10], [x_, y_], t_)))`. *expect:* `{'local': '+1.00e-04', 'advective': '-1.00e-04', 'total': '+0.00e+00'}` ·
    `station warming: 0.36 K/h` · `(1e-04, -1e-04, ~0)` (to 1e-12) · `-1e-04` (= |u| ∂T/∂s: the advective term) · `0`.
    *explain:* 1. the exact terms of our thermal front (the field is carried by the wind, so the parcel's rate is zero);
    2. the same numbers from general second-order stencils, which work for any F and u; 3. (3.6)'s streamwise form
    equals the advective term; 4. sympy writes (3.5) symbolically and gets 0.
29. `nb.check_agree` — **from scratch (curation §7):** hand-written central differences at the point:
    `ht, hx = 1.0, 10.0`; `local = (F(p, ht) - F(p, -ht))/(2*ht)`; `dFdy = (F(p + [0, hx], 0) - F(p - [0, hx], 0))/(2*hx)`;
    `adv = v*dFdy`; `assert np.allclose([local, adv], [terms["local"], terms["advective"]], rtol=1e-8)` · the independent
    route: sample F along `ch03.pathline(u, [0, 0], 0.0, [-60, 0, 60])` and difference in time →
    `assert abs(dF_along) < 1e-9` (the parcel's own rate = DT/Dt = 0). Markdown: "Two routes, one number: stencils at a
    fixed point, and a thermometer riding the path line."
30. `nb.figure` — two panels for a *front of finite width* (w = 100 km, largest gradient 1 K/100 km, front moving north
    at c = 5 m/s, wind v = 10 m/s, no heating; `thermal_front(…, width_m=1e5)`): (a) T(y) at t = 0, 6, 12 h (blue
    shades) with the float's positions (purple dots) and the station at y = 0 (blue square); (b) T(t) read by the station
    (blue) and by the float (purple) over 24 h. *see:* "a front sliding north and a float overtaking it; the station's
    curve rises as the front passes, the float's curve falls as it crosses into colder air." *read:* "the station sees
    ∂T/∂t = +Gc sech²(…) (the pattern moving past); the float sees DT/Dt = G(c − v) sech²(…) < 0 because it outruns the
    pattern. Where they read the same slope, the advective term is zero." *change:* "…the front moved exactly with the
    wind (c = v): the float's curve would be flat (DT/Dt = 0) and the station's warming would be pure advection — the
    worked example above."
31. `nb.plotly` — `slider_figure` over time t = 0…24 h (25 steps, FAST 13): x = time [h], y = rate [K/h]; traces "∂T/∂t
    (station, blue)", "u·∇T (amber)", "DT/Dt (float, purple)" evaluated along the float's path and drawn up to the
    slider time, with a faint full-day ghost of each. Title "The three terms along the float: the balance shifts as it
    crosses the front". Followed by `nb.md` see/read/change: "…all three peak when the float is at the front's centre;
    the purple curve is the sum of the other two at every instant; drag to t ≈ 11 h (the crossing)."
32. `nb.explainer("material_derivative_probe", heading="Why does the station warm while the air does not?", why="Two
    observers read one moving temperature pattern on one clock: a fixed probe (∂T/∂t) and a float carried by the wind
    (DT/Dt). You set the wind, the gradient and any heating; the term bars show the local and advective parts adding to
    the material derivative $\frac{DT}{Dt}=\frac{\partial T}{\partial t}+\mathbf u\cdot\nabla T$ (3.5), and the float's
    measured rate lands on the purple bar.", tries=["Pick the preset 'pure advection': the float's reading stays flat
    while the probe warms — read the numbers in Explain.", "Turn the wind until it blows along the isotherms: the amber
    bar vanishes (u ⟂ ∇T).", "Switch to the stretching-map mode (x = X e^{αt}): now F is the velocity itself and the
    float's rate is its acceleration α²x.", "Open the Derivation tab and step through D02 with your wind speed."])`
33. `nb.md` — **What would change if…** "…F were a velocity component instead of a temperature? Then D/Dt of u is the
    particle's acceleration, $\frac{D\mathbf u}{Dt}=\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u$ —
    the left side of Newton's law in Ch. 4. Section 3.3 asks what these rates look like as curves in the flow, and
    whether the split into local and advective parts depends on who is watching (C05)."

### A.3 §3.3 Flow Lines, Fluid Acceleration, and Galilean Transformation — C03 (+N14, N15 · D03), C04 (+N13, N16–N19 · D04, D05 · E1), C05 (+N04, N20–N24 · D06 · E3)
1. `nb.section("3.3", "Flow Lines, Fluid Acceleration, and Galilean Transformation", intro="**What is this section
   about?** Three curves that picture a flow — the instantaneous direction (streamline), the track of one particle
   (path line) and the line of dye from a fixed port (streak line) — and the proof that a particle's acceleration is
   the same for every observer moving at constant velocity, even though 'steady' is not.")`

**C03 — Streamlines (3.7)**
2. `nb.core("C03", "Streamlines: the flow's direction at one frozen instant (3.7)", question="Freeze the flow for an
   instant and draw curves that follow the velocity arrows everywhere. What equations do those curves obey, and why do
   they change from one instant to the next?")`
3. `nb.md` — **The problem in plain words:** "A satellite image of cloud streaks or a picture of iron filings along a
   magnet's field lines shows *directions at one instant*. A weather chart's wind streamlines are the same idea. To draw
   one we need a rule: at every point of the curve, its tangent points along u at that instant."
4. `nb.md` — **The idea:**
   ```
   freeze the clock at t   →   field of arrows u(x, t)   →   curves everywhere tangent to the arrows
   ds ∥ u   ⇔   dx/u = dy/v = dz/w   ⇔   u × ds = 0            (unsteady flow: a new picture every instant)
   ```
5. `nb.primer("parametric curves, tangent vector and arc length", "A curve can be written $\mathbf x(s)$ with a
   parameter s; its tangent is $d\mathbf x/ds$. If s is the **arc length** (distance measured along the curve),
   the tangent has length 1. To draw a curve tangent to u, march along $d\mathbf x/ds=\mathbf u/|\mathbf u|$.",
   code="s = np.linspace(0, 2*np.pi, 5)\nx, y = np.cos(s), np.sin(s)               # unit circle, s = arc length\ntx, ty =
   -np.sin(s), np.cos(s)             # tangent dx/ds\nprint(np.hypot(tx, ty))                   # [1. 1. 1. 1. 1.]: unit
   length")` — *expect:* `[1. 1. 1. 1. 1.]`.
6. `nb.primer("parallel vectors and the cross-product test", "Two vectors are parallel when one is a multiple of the
   other, $\mathbf a=\lambda\mathbf b$ — equivalently when $\mathbf a\times\mathbf b=0$ (the cross product, Ch. 2 C08,
   measures the area they span). The ratio test $a_x/b_x=a_y/b_y=a_z/b_z$ says the same but fails when a component of
   b is zero; the cross product never divides.", code="a, b = np.array([2.0, 4.0, 0.0]), np.array([1.0, 2.0, 0.0])
   \nprint(np.cross(a, b), a[0]/b[0], a[1]/b[1])   # [0 0 0] 2.0 2.0: parallel")` — *expect:* `[0. 0. 0.] 2.0 2.0`.
7. `nb.derivation("D03", …)` — Part F D03 (5 steps), ref "3.7".
8. `nb.note` — **N14 [B]** "The streamline equations are the result just derived:" equation `dx/u=dy/v=dz/w`, ref "3.7".
   "In a plane flow the first equality gives the slope $dy/dx=v/u$. Integrating (3.7) from many starting points,
   upstream and downstream, fills the picture. Number: where u = (1, 2) m/s the streamline climbs with slope 2."
9. `nb.worked_example("slopes and a circle", "1. At a point where u = 1 m/s, v = 2 m/s: $dy/dx=v/u=2$ (63.4° to the
   x-axis). 2. Solid-body rotation u = −y, v = x (s⁻¹ × m): $dy/dx=v/u=-x/y$, so $y\,dy=-x\,dx$ and
   $x^2+y^2=\text{const}$ — circles. 3. Ex. 3.1 at the instant t′ = π/4 s (ω = 1 s⁻¹): $v/u=\tan(\omega t')=1$ at every
   point, so every streamline is a straight line at 45°.")`
10. `nb.primer("solve_ivp options: t_eval, dense_output, events, backward integration", "Ch. 1 used `solve_ivp` with its
    defaults (P31). Here: `t_eval` asks for the solution at chosen times; `dense_output=True` returns a continuous
    solution `sol.sol(t)`; an `events` function stops the integration where it crosses zero (we stop a streamline at a
    stagnation point, |u| = 0); a time span that runs backwards (`t_span=(0, -5)`) integrates into the past — how a
    streamline is drawn upstream and a streak particle is traced back to the port.", code="from scipy.integrate import
    solve_ivp\nstop = lambda t, y: y[0] - 0.5; stop.terminal = True        # event: stop when y reaches 0.5\nsol =
    solve_ivp(lambda t, y: -y, (0, 5), [1.0], events=stop, dense_output=True)\nprint(sol.t_events[0], sol.sol(0.3))
    # [0.693] (= ln 2) and y(0.3) = e^-0.3 ≈ 0.741")` — *expect:* `[0.69314718] [0.74081822]` (to solver tolerance).
11. `nb.code` — *code:* `u31 = ch03.preset_field("ex31", omega=1.0, xi0=1.0)` · for `tp in (0.0, np.pi/4, np.pi/2)`:
    `sl = ch03.streamline(u31, [0.0, 0.0], tp, s_max=2.0)`; `print(tp, np.round(sl[:, -1], 3))` ·
    `print(ch03.streamline_slope(u31, np.array([0.3, -0.7]), np.pi/4))` · solid body:
    `usb = lambda x, t: np.array([-x[1], x[0]])`; `c = ch03.streamline(usb, [1.0, 0.0], 0.0, s_max=2*np.pi)`;
    `print(np.ptp(np.hypot(c[0], c[1])))`. *expect:* end points (2.0, 0.0), (1.414, 1.414), (0.0, 2.0) (straight lines at
    0°, 45°, 90°) · `1.0` · radius spread < 1e-8 (a circle). *explain:* 1. the Ex. 3.1 field as a callable; 2.
    `streamline` integrates $d\mathbf x/ds=\mathbf u/|\mathbf u|$ with the clock frozen at t′, both directions from the
    seed; 3. the slope v/u of (3.7), $dx/u=dy/v$, at an arbitrary point; 4. a closed circle for solid-body rotation (the
    radius stays constant along the computed curve).
12. `nb.figure` — three small panels: Ex. 3.1's streamline pattern (teal, `plt.streamplot` on a 21 × 21 grid) at
    t′ = 0, π/4, π/2 s, the port at the origin, the velocity arrow at the port; a fourth panel: solid-body streamlines
    (circles). *see:* "parallel teal lines that turn by 45° from panel to panel; circles in the last panel." *read:*
    "Ex. 3.1's velocity is the same at every point at a given instant, so its streamlines are parallel straight lines;
    only their direction changes with time — the pattern turns once per period 2π/ω." *change:* "…the flow were steady
    (freeze ω t′): all three panels would be identical — and so, C04 shows, would path and streak lines."
13. `nb.plotly` — `slider_figure` over the drawing instant t′ ∈ [0, 2π] s (24 steps): the Ex. 3.1 streamline through the
    origin (teal line), with the path-line circle (orange) and streak-line circle (rose) **for the same t′** (from
    `example_3_1`); fixed axes ±2.5 m. Title "Streamlines turn with t′; path and streak lines are circles that roll round
    the port". `nb.md` see/read/change: "…at every t′ the teal line touches both circles at the origin (the tangency
    shown in D05); drag t′ a quarter period and all three turn by 90°."
14. `nb.note` — **N15 [B]** "**Stream tube.** The streamlines through every point of a closed curve C form a tube (Fig.
    3.6). No fluid crosses its wall, because the wall is everywhere tangent to u; in a steady flow the same volume per
    second passes every cross-section. Number below: the axisymmetric straining flow u = (−x/2, −y/2, z) s⁻¹ (∇·u = 0)
    with a tube seeded on a circle of radius 0.5 m at z = 1 m: the tube narrows as $R(z)=0.5\sqrt{1/z}$ and the flux is
    $\pi/4=0.785$ m³/s through every section."
15. `nb.plotly` — 3-D stream tube: 12 streamlines (`ch03.streamline`, teal `go.Scatter3d` lines) seeded on the circle,
    two cross-section discs at z = 1 and z = 2 m (translucent `go.Mesh3d`), and a printed check in the cell:
    `q1 = ch03.flux_through_disc(u3, [0, 0, 1], [0, 0, 1], 0.5)`; `q2 = ch03.flux_through_disc(u3, [0, 0, 2], [0, 0, 1],
    0.5/np.sqrt(2))`; `print(q1, q2)`. *expect:* `0.7854 0.7854`. Gloss in a comment: "plotly 3-D lines: one
    `go.Scatter3d(mode='lines')` per streamline (Ch. 2 P64)". `nb.md` see/read/change: "…the tube narrows as the flow
    speeds up along z; the two discs carry the same flux — the wall lets nothing through. A wider seed circle gives a
    fatter tube and a proportionally larger (πR²) flux."
16. `nb.md` — **What would change if…** "…you followed one particle instead of freezing the clock? In a steady flow it
    would run along a streamline. In Ex. 3.1 the arrows turn while the particle moves, so its track is not any of the
    straight teal lines — it is a circle (C04)."

**C04 — Path lines and streak lines; Ex. 3.1**
17. `nb.core("C04", "Path lines and streak lines: where one particle goes, where the dye sits (3.8)", question="Dye
    leaks steadily from a port on the floor of a sloshing tank, and you track one speck of it. Why does the photograph of
    the dye not show the path of the speck?")`
18. `nb.md` — **The problem in plain words:** "Every flow picture you see in a lab is one of three things: a long
    exposure of one particle (a **path line**), a snapshot of all the dye that came out of one port (a **streak line**),
    or a snapshot of directions (a **streamline**). Near a beach, surface water under long swell sloshes back and forth;
    a drop of dye released there and a floating leaf trace quite different curves. We need both as computable objects."
19. `nb.md` — **The idea:**
   ```
   | line        | what is fixed                  | what varies along it     | how you see it                      |
   | streamline  | the instant t                  | position (arc length s)  | arrows at one instant               |
   | path line   | the particle (label r_o, t_o)  | time t                   | long exposure of one speck          |
   | streak line | the port x_o and the instant t | release time t_o         | photo of continuously injected dye  |
   steady flow: all three are the same curve.   unsteady: three different curves.
   ```
20. `nb.note` — **N16 [B]** "A **path line** is the trajectory of one particle of fixed identity, the Lagrangian
    $\mathbf x=\mathbf r(t;\mathbf r_o,t_o)$ of C01 drawn in space. Given only the Eulerian velocity field, we get it by
    solving the ODE below."
21. `nb.md` — **The maths — the path-line equation.** "Each particle moves with the velocity of the field at the place it
    currently occupies:
    $$\frac{d\mathbf r}{dt}=[\mathbf u(\mathbf x,t)]_{\mathbf x=\mathbf r}=\mathbf u(\mathbf r,t),\qquad \mathbf r(t_o)=\mathbf r_o\qquad(3.8)$$
    — an initial-value problem: three coupled ODEs, one per component. (A discretised version of (3.8) is how particle
    image velocimetry, PIV, turns pairs of photos into velocities.)"
22. `nb.note` — **N17 [B]** "A **streak line** through the port $\mathbf x_o$ at time t is the set of particles that passed
    through $\mathbf x_o$ at earlier times $t_o$: solve (3.8), $d\mathbf r/dt=\mathbf u(\mathbf r,t)$, once per release
    time with $\mathbf r(t_o)=\mathbf x_o$, then plot all positions at the same t. The parameter along it is the release
    time, $x_i=r_i(t;\mathbf x_o,t_o)$. Sometimes $t_o$ can be eliminated to give one equation for the curve (D05)."
23. `nb.note` — **N18 [B]** "**Ex. 3.1** (a spatially uniform, time-periodic flow — a caricature of the back-and-forth
    motion under long surface waves, not a wave solution; real wave orbits shrink with depth, Ch. 7): $u=\omega\xi_o\cos\omega t$,
    $v=\omega\xi_o\sin\omega t$. Find the three lines through the origin at the instant t = t′. ⚠️ Four names that look
    alike: **t′** the drawing instant, **t_o** a release time, **c_x, c_y** integration constants (the book's x_o,
    y_o), **ξ_o** the amplitude [m]." The streamline comes straight from (3.7), $dx/u=dy/v$ (curation §4c): equation
    `\frac{dy}{dx}=\frac{v}{u}=\frac{\omega\xi_o\sin(\omega t')}{\omega\xi_o\cos(\omega t')}=\tan(\omega t')\ \Rightarrow\ y=x\tan(\omega t')`,
    no ref (integrated through the origin, the constant is 0).
24. `nb.derivation("D04", …)` — Part F D04 (7 steps), ref "".
25. `nb.derivation("D05", …)` — Part F D05 (8 steps), ref "".
26. `nb.worked_example("Ex. 3.1 with ω = 1 s⁻¹, ξ_o = 1 m, t′ = 0", "1. Streamline: $y=x\tan0=0$ — the x-axis. 2. Path
    line (D04): $(x+\sin0)^2+(y-\cos0)^2=1$, i.e. $x^2+(y-1)^2=1$: a unit circle centred at (0, 1). 3. Streak line (D05):
    $(x-\sin0)^2+(y+\cos0)^2=1$, i.e. $x^2+(y+1)^2=1$: a unit circle centred at (0, −1). 4. All three touch the origin
    with slope $\tan0=0$ (horizontal): the circles lie on opposite sides of the x-axis, like Fig. 3.7. 5. A quarter
    period later (t′ = π/2 s) everything has turned by 90°: the streamline is the y-axis, the centres are (−1, 0) and
    (1, 0).")`
27. `nb.primer("RK4 by hand", "The classical fourth-order Runge–Kutta step for $d\mathbf r/dt=\mathbf f(\mathbf r,t)$
    samples the slope four times per step: $k_1=f(r,t)$, $k_2=f(r+\tfrac{h}{2}k_1,t+\tfrac h2)$,
    $k_3=f(r+\tfrac h2k_2,t+\tfrac h2)$, $k_4=f(r+hk_3,t+h)$, then $r\leftarrow r+\tfrac h6(k_1+2k_2+2k_3+k_4)$. Its error
    per unit time falls like h⁴ (Euler's like h, Ch. 1 P30). The explainers use exactly this step in JavaScript.",
    code="f = lambda y, t: -y                          # dy/dt = −y, exact y = e^−t\ny, h = 1.0, 0.1\nfor n in range(10):
    \n    k1 = f(y, n*h); k2 = f(y + h/2*k1, n*h + h/2); k3 = f(y + h/2*k2, n*h + h/2); k4 = f(y + h*k3, n*h + h)\n    y +=
    h/6*(k1 + 2*k2 + 2*k3 + k4)\nprint(y, np.exp(-1.0))                        # 0.3678798 vs 0.3678794")` — *expect:*
    `0.36787977 0.36787944`.
28. `nb.code` — *code:* `d = ch03.example_3_1(0.0, xi0=1.0, omega=1.0)` · `print(d["path_center"], d["streak_center"],
    d["slope"])` · numerical path line of the particle at the origin at t′ = 0, half a period back and forward:
    `T = 2*np.pi; ts = np.linspace(-T/2, T/2, 201)`; `p = ch03.pathline(u31, [0.0, 0.0], 0.0, ts)`;
    `print(np.max(np.abs(np.hypot(p[0] - d["path_center"][0], p[1] - d["path_center"][1]) - 1.0)))` · streak line: dye
    released over the last period: `tr = np.linspace(-T, 0.0, 181)`; `s = ch03.streakline(u31, [0.0, 0.0], 0.0, tr)`;
    `print(np.max(np.abs(np.hypot(s[0] - d["streak_center"][0], s[1] - d["streak_center"][1]) - 1.0)))`. *expect:* `[-0.
    1.] [0. -1.] 0.0` · path residual < 1e-8 · streak residual < 1e-7. *explain:* 1. the closed forms of D04 and D05; 2.
    `pathline` solves (3.8), $d\mathbf r/dt=\mathbf u(\mathbf r,t)$, for one particle; every point is 1 m from the centre
    (0, 1); 3. `streakline` solves (3.8) once per release time (one vector ODE for the whole ensemble) and collects the
    positions at t′ = 0: a circle about (0, −1).
29. `nb.check_agree` — **from scratch (curation §7):** the RK4 loop of the primer applied to Ex. 3.1:
    `r = np.array([0.0, 0.0]); h = T/400`; loop 200 steps forward from t = 0 with `f = lambda r, t: u31(r, t)`; compare
    with `ch03.pathline(u31, [0, 0], 0.0, [0.0, T/2])[:, -1]` and with the circle `(sin t, 1 − cos t)` →
    `assert np.allclose(r, lib_end, atol=1e-8)`; `assert np.allclose(r, [np.sin(T/2), 1 - np.cos(T/2)], atol=1e-8)`.
    Markdown: "Twenty lines of arithmetic reproduce the library to 10⁻⁸ — no magic inside `pathline`."
30. `nb.figure` — **N19 [B]** our Fig. 3.7 (`flow_lines_figure(t_prime=np.pi/6)` from `scripts/ch03_fig3_7_flow_lines.py`):
    streamline (teal, dashed beyond the port), path-line circle (orange, with its direction arrows counterclockwise),
    streak-line circle (rose), centres marked ×, the angle ωt′ drawn at the origin, ξ_o radii labelled. *see:* "one
    straight line and two circles of equal radius on opposite sides, all touching at the origin." *read:* "same point,
    same instant, three different curves — because the flow is unsteady. The common tangent at the origin is the
    velocity there at t′ (D05's last steps)." *change:* "…t′ larger by π/2ω: the whole picture turns by 90°
    counterclockwise; the circles keep their radius ξ_o."
31. `nb.animation` — over one period (60 frames, FAST 30, `player="video"`; `flow_lines_animation`): dye particles
    emitted from the port every T/60 (rose dots, oldest faded), one tagged particle with its growing orange trail, the
    teal streamline through the port turning, a clock readout "t = … s". Followed by `nb.md` see/read/change: "…the dye
    always lies on a circle on the far side of the port from the tagged particle's circle; after one period the particle
    is back at the port and the whole pattern repeats. With a mean current (E1's second mode) the particle would not come
    back: its path becomes a looping curve drifting downstream."
32. `nb.note` — **N13 [B]** "**Steady flow ⇒ the three lines coincide.** If u does not depend on t, a particle at a point
    moves along the streamline through it and every later particle from the same port follows the same track. Check on
    the steady hyperbolic flow u = (x, −y) s⁻¹ from the port (1, 1) m:" followed by `nb.code` — *code:* `us =
    ch03.preset_field("rotating_strain", s=1.0, Omega=0.0)` · streamline, path line (t ∈ [0, 1] s) and streak line
    (releases t_o ∈ [0, 1] s, drawn at t = 1 s) through (1, 1); compare each against the exact curve xy = 1 by
    `np.max(np.abs(c[0]*c[1] - 1))`. *expect:* three residuals < 1e-7. *explain:* "In a steady flow all three curves lie on
    the same hyperbola xy = 1. Flow-visualisation atlases (Van Dyke's *Album of Fluid Motion*) are full of such
    pictures; in unsteady flow you must know which of the three a photograph shows."
33. `nb.explainer("flow_lines_unsteady", heading="Three lines through one point — why do they disagree?", why="A clock
    runs: the streamline pattern is redrawn every instant, one particle leaves its trail and a port leaks dye. You
    change the period, the amplitude and a mean current, or switch to a steady flow — and watch the three curves
    separate or collapse onto one.", tries=["Press ▶ with the Ex. 3.1 preset and stop at ωt = 30°: compare with Fig. 3.7
    above.", "Choose 'steady': the three curves become one line — the N13 statement as an experiment.", "Add a mean
    current U₀ = 0.5 m/s: the particle no longer returns to the port; the dye forms a wave.", "Click a dye particle
    (inspector) to see when it left the port and the path it took."])`
34. `nb.md` — **What would change if…** "…you watched Ex. 3.1 from a boat drifting with the mean current? The
    streamlines would look different again — whether a flow is steady depends on the observer. The next block asks
    whether anything important depends on the observer: the particle's acceleration does not."

**C05 — Galilean invariance of the acceleration (3.9)**
35. `nb.core("C05", "Galilean invariance: the acceleration does not care who is watching (3.9)", question="A cylinder is
    towed through a still lake; a pier stands in a steady river. One flow is unsteady, the other steady, yet they are
    the same flow seen from two boats. Does a water particle's acceleration depend on the boat you watch it from?")`
36. `nb.md` — **The problem in plain words:** "Wind-tunnel engineers hold a model still and blow air past it; the real
    car moves through still air. Oceanographers describe a wave in a frame moving with the crest, where the flow is
    steady. All of this is legitimate only if Newton's law — force equals mass times *acceleration* — reads the same in
    both frames. We check that the acceleration, written as local + advective, is frame-independent for observers
    moving at constant velocity."
37. `nb.md` — **The idea:**
   ```
   body frame (pier fixed):   ∂u/∂t = 0 (steady)          (u·∇)u ≠ 0   ─┐
   fluid frame (lake fixed):  ∂u'/∂t' ≠ 0 (unsteady)      (u'·∇')u'     ├─ same particle, same total  Du/Dt
   only U (constant) separates them: the split moves, the sum does not   ─┘
   ```
38. `nb.primer("frames of reference and relative velocity", "A frame of reference is a set of axes plus a clock. Two
    frames whose axes stay parallel and whose origins separate at a constant velocity U are related by a Galilean
    transformation: positions differ by Ut (+ a fixed offset), velocities by U, and both use the same clock (no
    relativity here). A frame that rotates is different: its points move with Ω × x, which is not uniform.",
    code="U = np.array([2.0, 0.0])                     # frame O' moves at 2 m/s along x\nu_lab = np.array([2.5, 0.3])
    # a particle's velocity seen from O\nprint(u_lab - U)                           # [0.5 0.3]: the same particle seen
    from O'")` — *expect:* `[0.5 0.3]`.
39. `nb.note` — **N04 [B]** "**Fig. 3.2, one flow in two frames.** Ideal flow past a circular cylinder of radius a in a
    stream U (the uniform stream plus a doublet; derived in Ch. 6 §6.3): $u=U[1-a^2(x^2-y^2)/r^4]$,
    $v=-2Ua^2xy/r^4$. Seen from the cylinder the streamlines are fixed (steady). Seen from the still fluid, u′ = u − U:
    streamlines start and end on the moving body and change as it passes (unsteady). The two velocity fields differ by
    the constant vector U:" equation `\mathbf u=\mathbf U+\mathbf u'`, no ref.
40. `nb.figure` — two panels (`ch03.cylinder_flow` on a 121 × 81 grid, `plt.streamplot`, U = 1 m/s, a = 1 m): (a) body
    frame: streamlines around the cylinder (steady); (b) fluid frame at t = 0: streamlines of u′ = u − U (loops leaving
    the front and entering the back of the body); in both, the point P = (0, 1.5) m marked with its velocity (black) and
    its acceleration (purple, identical in both panels: (0, −0.856) m/s²). *see:* "two completely different streamline
    pictures and one identical purple arrow." *read:* "(a) is steady — the arrows never change at a fixed point; (b) is
    unsteady — the body moves through it. The acceleration at P (computed below) is the same vector in both."
    *change:* "…the observer moved at half the towing speed: a third streamline picture, the same purple arrow (E3's
    slider)."
41. `nb.note` — **N21 [B]** "**The Galilean transformation** (Fig. 3.8). Frame O′x′y′z′ moves at constant velocity U
    relative to Oxyz with parallel axes; $\mathbf x'_o$ is the vector from O to O′ at t = 0:" equation
    `\mathbf u(\mathbf x,t)=\mathbf U+\mathbf u'(\mathbf x',t'),\qquad t=t',\qquad \mathbf x=\mathbf x'+\mathbf Ut+\mathbf x'_o`,
    no ref. "`ch03.galilean_transform(u, U)` returns the primed field u′(x′, t′) = u(x′ + Ut′ + x′_o, t′) − U."
42. `nb.primer("chain rule with a moving frame", "If $g(x,t)=f(x-Ut,t)$, then at fixed x,
    $\frac{\partial g}{\partial t}=\frac{\partial f}{\partial t'}-U\frac{\partial f}{\partial x'}$: holding x fixed while
    t advances means the moving coordinate x′ = x − Ut *decreases*. A time derivative 'at fixed position' depends on
    whose position is held fixed — the classic trap of this section.", code="x, t, U = sp.symbols('x t U')\nxp, tp =
    sp.symbols('xp tp')\nf = sp.sin(xp)*sp.exp(-tp)                         # any f(x', t')\ng = f.subs({xp: x - U*t,
    tp: t})                     # the same field seen in the unprimed frame\nrhs = (sp.diff(f, tp) - U*sp.diff(f,
    xp)).subs({xp: x - U*t, tp: t})\nprint(sp.simplify(sp.diff(g, t) - rhs))              # 0")` — *expect:* `0`.
43. `nb.derivation("D06", …)` — Part F D06 (9 steps), ref "3.9". **check_src** (optional, recommended — Part F D06 sketches
    it): sympy with a general 2-D field u′(x′, y′, t′) (two `sp.Function`s), U = (U₁, U₂) symbols.
44. `nb.note` ×3 — **N20 [B]** "Back to Fig. 3.2: in the body frame the unsteady term is zero but the streamlines bend,
    so particles still accelerate — all of it advective. In the fluid frame both terms are nonzero. Because the frames
    differ by a constant U, (3.9), $\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=\frac{\partial\mathbf u'}{\partial t'}+(\mathbf u'\cdot\nabla')\mathbf u'$,
    guarantees the same acceleration at the same place relative to the cylinder." · **N23 [B]** "**The split is
    observer-dependent, the sum is not** — the term bars below trade blue for amber as the observer's speed changes,
    while the purple total stays put." · **N24 [C]** "Two everyday pictures: a traffic light that switches from east–west
    to north–south traffic is an *unsteady* (local) change at a fixed crossing; a roller coaster on a fixed track gives
    its riders *advective* acceleration although the track never changes."
45. `nb.note` — **N22 [B]** "**Linear and quadratic.** The local term ∂u/∂t is linear in u; the advective term (u·∇)u is
    quadratic — doubling every velocity doubles the first and quadruples the second. This nonlinearity is the heart of
    fluid mechanics (turbulence, waves steepening). When u is tiny it can be dropped (acoustics, Ch. 15); when u = 0
    only statics is left (Ch. 1 §1.7)." followed by `nb.code` — *code:* scale the fluid-frame cylinder field by λ = 2 at
    P, t = 0 with `ch03.acceleration` → `print(loc2/loc1, adv2/adv1)`. *expect:* `2.0 4.0` (y-components). *explain:*
    "the scaling test in two numbers."
46. `nb.worked_example("a travelling sine wave, two frames", "In a frame moving with the wave the flow is steady:
    $u'=0.5\sin x'$ m/s (k = 1 m⁻¹). The lab sees the same pattern moving at U = 2 m/s: $u=2+0.5\sin(x-2t)$. At the lab
    point $x-2t=\pi/4$: 1. Wave frame: local 0; advective $u'\,du'/dx'=0.5\sin\frac\pi4\times0.5\cos\frac\pi4=0.3536\times0.3536=0.125$
    m/s². 2. Lab frame: local $\partial u/\partial t=-2\times0.5\cos\frac\pi4=-0.707$ m/s²; advective
    $u\,\partial u/\partial x=(2+0.354)\times0.354=0.832$ m/s². 3. Sum $-0.707+0.832=0.125$ m/s² — the same, as (3.9)
    says. 4. The cylinder at P = (0, 1.5) m (U = a = 1): body frame local 0, advective −0.856; fluid frame local −0.593,
    advective −0.263 (y-components, m/s²) — again −0.856 in total.")`
47. `nb.code` — *code:* wave: `U = np.array([2.0]); up = lambda xp, tp: 0.5*np.sin(xp)` (steady in the wave frame) ·
    lab field `u = lambda x, t: U + up(x - U*t, t)` · `print(ch03.acceleration(u, np.array([np.pi/4]), 0.0))` ·
    `print(ch03.acceleration(ch03.galilean_transform(u, U), np.array([np.pi/4]), 0.0))` · cylinder:
    `for Uf in (1.0, 0.5, 0.0): print(Uf, ch03.frame_acceleration_terms(0.0, 1.5, 1.0, 1.0, Uf))`. *expect:* lab (total,
    local, advective) = (0.125, −0.707, 0.832) · wave frame (0.125, 0.0, 0.125) · U_frame = 1 (body): local (0, 0),
    advective (0, −0.856), total (0, −0.856); U_frame = 0.5: local (0, −0.296), advective (0, −0.560), total
    (0, −0.856); U_frame = 0 (fluid): local (0, −0.593), advective (0, −0.263), total (0, −0.856). *explain:* 1. a
    field that is steady in the moving frame; 2. `acceleration` returns (a, local, advective) by stencils in the lab;
    3. `galilean_transform` builds the primed field and the same call gives the wave-frame split; 4.
    `frame_acceleration_terms` does it for the cylinder for three observers — the totals agree to 10⁻⁸.
48. `nb.figure` — term bars at P for five observers U_frame = 0, 0.25, 0.5, 0.75, 1 m/s: stacked local (blue) and
    advective (amber) y-components with the total (purple diamond) on each; a dashed purple line at −0.856 m/s².
    *see:* "blue bars shrinking to zero as the observer approaches the body frame, amber bars growing, purple diamonds on
    one horizontal line." *read:* "reading left to right is moving the observer from the lake to the cylinder; the
    acceleration of the particle at P never changes — only how it is booked." *change:* "…a point far from the body
    (P = (0, 5) m): all bars shrink toward zero (the flow there is almost uniform)."
49. `nb.plotly` — `slider_figure` over the observer velocity U_frame ∈ [0, 1] m/s (21 steps): the streamlines of the
    cylinder flow in that frame (precomputed with `ch03.streamline` from 14 seeds, ≤ 4 traces × ≤ 400 points: upper
    streamlines, lower streamlines, body outline, the point P) — title "The observer's speed reshapes the streamlines"
    (the local and advective values at P printed in each step's trace name). `nb.md` see/read/change: "…the loops of the
    fluid frame open into the flow-past-a-body picture as U_frame → U; the numbers in the legend trade places while their
    sum stays −0.856 m/s²."
50. `nb.explainer("galilean_frames_cylinder", heading="Steady or not — does the acceleration care?", why="One slider
    moves the observer from the still fluid to the cylinder. The streamline picture morphs from 'unsteady, loops on the
    body' to 'steady, around the body', while the acceleration arrow of a tagged particle stays exactly the same and its
    two bars trade places — (3.9), $\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=\frac{\partial\mathbf u'}{\partial t'}+(\mathbf u'\cdot\nabla')\mathbf u'$,
    as an experiment.", tries=["Drag the observer speed from 0 to U and watch the blue bar hand its share to the amber
    bar.", "Click any point: the inspector gives both frames' terms with numbers.", "Set the preset 'body frame': the
    status says steady — yet the purple arrow is not zero.", "Step through D06 in the Derivation tab; at the
    cancellation step the two U·∇′ bars light up."])`
51. `nb.md` — "> ⚠️ **Common confusion:** 'Galilean invariance means any moving frame is as good as any other.' Only
    frames moving at *constant* velocity without rotating. In a rotating frame (the Earth!) the acceleration picks up
    Coriolis and centrifugal terms (Ch. 4 §4.7) and even the vorticity changes by 2Ω (N28 in C10) — the heart of
    geophysical fluid dynamics (Ch. 13)."
52. `nb.md` — **What would change if…** "…the observer accelerated (U depends on t)? Step 3 of D06 gains −dU/dt, which
    does not cancel: an accelerating observer sees a fictitious force. Next, §3.4 zooms into a small neighbourhood of a
    point and asks what the velocity *differences* do to a fluid element."

### A.4 §3.4 Strain and Rotation Rates — R01–R03, C06 (+N25, N26 · D07), C07 (D08), C08 (+N27 · D09, D10), C09 (D11 · E4), R04–R07, C10 (+N28, N29 · D12, D14, D13), C11 (+N30, N33 · D15), C12 (+N31, N32 · D16 · E5)
1. `nb.section("3.4", "Strain and Rotation Rates", intro="**What is this section about?** Zoom into a tiny blob of fluid.
   Its neighbours move at slightly different velocities; those differences stretch it, shear it, swell it and spin it.
   Each of these motions is one piece of the velocity-gradient tensor you met in Ch. 2 — here you *measure* them on a
   moving element.")`
2. `nb.recap("R01", "The split of the velocity gradient", "Any velocity gradient splits uniquely into a symmetric and an
   antisymmetric part: $\frac{\partial u_i}{\partial x_j}=S_{ij}+\tfrac12R_{ij}$ *(Eq. 3.11)* (Ch. 2 D14 proved the
   split is unique). `ch03.strain_rate_tensor(G) + 0.5*ch03.rotation_tensor(G)` returns G exactly.", where="Ch. 2 §2.10,
   C12")`
3. `nb.recap("R02", "The strain-rate tensor", "The symmetric part $S_{ij}=\tfrac12\Big(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\Big)$
   *(Eq. 3.12)* — `ch03.strain_rate_tensor(G)`. §3.4 shows what each of its entries *measures*.", where="Ch. 2 §2.10,
   C12")`
4. `nb.recap("R03", "The rotation tensor — the book's, without ½", "$R_{ij}=\frac{\partial u_i}{\partial x_j}-\frac{\partial u_j}{\partial x_i}$
   *(Eq. 3.13)*, `ch03.rotation_tensor(G)` = G − Gᵀ. ⚠️ The book's R is **twice** the antisymmetric part A of Ch. 2 (R
   = 2A); that is why a ½ stands in front of R in (3.11), $\frac{\partial u_i}{\partial x_j}=S_{ij}+\tfrac12R_{ij}$.",
   where="Ch. 2 §2.10, C12")`

**C06 — Relative velocity near a point (3.10)**
5. `nb.core("C06", "Relative velocity of a neighbour: du = G·dx (3.10)", question="Two corks float 1 cm apart. How
   fast does one move relative to the other — and what does that single matrix of nine numbers already tell us about
   the blob of water between them?")`
6. `nb.md` — **The problem in plain words:** "A drop of ink in a stream is stretched into a thin filament, a cloud is
   sheared into a streak, air in a rising thermal swells. All of these are about *differences* of velocity between
   neighbouring points. Chapter 4's stress law needs exactly these differences (the rates of deformation), not the
   velocity itself."
7. `nb.md` — **The idea:**
   ```
   O at x, velocity u          P at x + dx, velocity u + du          (Fig. 3.9)
   zoom in far enough and every smooth field looks linear:   du = G · dx,   G_ij = ∂u_i/∂x_j   (nine numbers)
   G = S + ½R :   S deforms the blob  (stretch, shear, swell)   ½R spins it
   ```
8. `nb.primer("multivariable first-order Taylor expansion", "Ch. 1's $f(z+dz)\approx f(z)+f'(z)\,dz$ (P26) in several
   variables: $u_i(\mathbf x+d\mathbf x)\approx u_i(\mathbf x)+\frac{\partial u_i}{\partial x_j}dx_j$ — one partial
   derivative per direction, summed over j. The neglected terms shrink like $|d\mathbf x|^2$, so after dividing by
   $|d\mathbf x|$ they vanish as the neighbour comes closer.", code="f = lambda x, y: np.sin(x)*np.exp(y)       # a smooth
   function\nx0, y0, dx, dy = 0.3, 0.1, 0.01, -0.02\nlin = f(x0, y0) + np.cos(x0)*np.exp(y0)*dx + np.sin(x0)*np.exp(y0)*dy
   \nprint(f(x0 + dx, y0 + dy) - lin)                  # ~1e-4: second order in the step")` — *expect:* a number of order
   10⁻⁴ (the builder prints the exact value).
9. `nb.derivation("D07", …)` — Part F D07 (4 steps), ref "3.10".
10. `nb.note` ×2 — **N25 [C]** "Why rates of deformation? A solid resists being *deformed* (strain); a fluid resists being
    deformed *quickly* (strain rate). Ch. 4 §4.5 relates the stress to S, the Newtonian law generalising Ch. 1's
    τ = μ du/dy." · **N26 [C]** "S (deformation) is linked to the stress in a moving fluid; R (rotation) is not — a
    rigidly spinning fluid feels no viscous stress. C11 shows why only S can enter."
11. `nb.worked_example("G = [[1, 2], [0, −1]] s⁻¹ and a neighbour 2 cm away", "1. dx = (0.01, 0.02) m. 2.
    $du_1=G_{11}dx_1+G_{12}dx_2=1\times0.01+2\times0.02=0.05$ m/s. 3. $du_2=G_{21}dx_1+G_{22}dx_2=0+(-1)\times0.02=-0.02$ m/s.
    4. Split: S = [[1, 1], [1, −1]], ½R = [[0, 1], [−1, 0]] s⁻¹, and S + ½R = [[1, 2], [0, −1]] = G ✓ (R01). 5. The field
    $u=(x_1+2x_2+x_1^2,\ -x_2+x_1x_2)$ has exactly this G at the origin; its true du at that neighbour is (0.0501,
    −0.0198) m/s — the difference (1, 2)×10⁻⁴ m/s is the second-order remainder.")`
12. `nb.code` — *code:* `ufield = lambda x, t: np.array([x[0] + 2*x[1] + x[0]**2, -x[1] + x[0]*x[1]])` ·
    `G = ch03.velocity_gradient_at(ufield, np.zeros(2), 0.0)` · `dx = np.array([0.01, 0.02])` · `print(G)` ·
    `print(ch03.relative_velocity(G, dx), ufield(dx, 0) - ufield(np.zeros(2), 0))` · `S, R =
    ch03.strain_rate_tensor(G), ch03.rotation_tensor(G)`; `assert np.allclose(S + 0.5*R, G)` · Taylor remainder:
    `e = dx/np.linalg.norm(dx); hs = 0.02*0.5**np.arange(6)`; `err = [np.linalg.norm(ufield(h*e, 0) - ufield(np.zeros(2),
    0) - G @ (h*e)) for h in hs]`; `print(observed_order(hs, err))`. *expect:* G = [[1, 2], [0, −1]] ·
    `[0.05 -0.02] [0.0501 -0.0198]` · order ≈ 2.00. *explain:* 1. a nonlinear field; 2. G by second-order stencils; 3.
    du = G·dx vs the exact difference; 4. the R01 split; 5. the remainder falls like |dx|² (slope 2 on log–log, Ch. 1
    P13).
13. `nb.figure` — two panels: (a) a ring of 24 neighbours at radius 0.1 m around O with their relative-velocity arrows
    du = G·dx (purple, `ring_arrows`), for the worked G — a "relative-velocity field" that is zero at O and grows with
    distance; (b) log–log |exact du − G·dx| vs |dx| with a slope-2 guide. *see:* "arrows pointing out along one diagonal
    and in along the other, plus a swirl; a straight line of slope 2." *read:* "(a) is everything the element does in
    the next instant — (b) says the linear picture becomes exact as the ring shrinks." *change:* "…G antisymmetric (a
    pure rotation): the ring's arrows would all be tangential — no stretching at all."
14. `nb.md` — **What would change if…** "…you watched only two neighbours on the x₁-axis? Their separation changes at
    ∂u₁/∂x₁ per unit length — the diagonal of S. That is C07."

**C07 — Linear strain rate**
15. `nb.core("C07", "Linear strain rate: how fast a material line stretches", question="A short thread of dye lies along
    the flow. At what rate, per unit of its length, is it being stretched — and in which direction would it stretch
    fastest?")`
16. `nb.md` — **The problem in plain words:** "Stretching is how ocean eddies draw tracers into long filaments and how
    vortex tubes intensify (Ch. 5). A rubber band held at both ends and pulled: the stretching rate per unit length is
    what matters, not the length itself."
17. `nb.md` — **The idea:**
   ```
   A ●────────● B     ends move at u₁ and u₁ + (∂u₁/∂x₁)δx₁   (Fig. 3.10)
   B outruns A by (∂u₁/∂x₁)δx₁ per second   →   stretch per length per time = ∂u₁/∂x₁ = S₁₁
   any direction n:   (1/ℓ) Dℓ/Dt = n·S·n
   ```
18. `nb.primer("material line element", "A *material* line element δx joins two nearby fluid particles and is carried
    with them — it stretches and turns as they move. Its rate of change following the particles is the difference of
    their velocities: $\frac{D(\delta\mathbf x)}{Dt}=\delta\mathbf u=\mathbf G\cdot\delta\mathbf x$ (C06). Every rate of
    §3.4 is measured on such elements.", code="G = np.array([[2.0, 0.0], [0.0, -2.0]])      # u = (2x, −2y)\ndx0 =
    np.array([0.01, 0.0])                         # a 1 cm element along x\ndx1 = ch03.linear_flow_map(G, 0.01) @ dx0
    # carried for 0.01 s\nprint(dx1)                                        # [0.010202 0.]: it grew by 2 %")` —
    *expect:* `[0.01020201 0.        ]`.
19. `nb.derivation("D08", …)` — Part F D08 (6 steps), ref "".
20. `nb.worked_example("u = (2x, −2y) s⁻¹, a 1 cm thread along x for 0.01 s", "1. $S_{11}=\partial u_1/\partial x_1=2$
    s⁻¹. 2. First order: new length $\approx1\times(1+2\times0.01)=1.02$ cm. 3. Exact (the ends move as $e^{2t}$):
    $e^{0.02}=1.0202$ cm. 4. A thread along y: $S_{22}=-2$ s⁻¹ — it shrinks to $e^{-0.02}=0.9802$ cm. 5. A thread at 45°:
    $\mathbf n\cdot\mathbf S\cdot\mathbf n=\tfrac12(2-2)=0$ — its length does not change at first order.")`
21. `nb.code` — *code:* `G = np.array([[1.0, 2.0], [0.0, -1.0]])` (C06's G) · for `deg in (0, 22.5, 45, 90)`:
    `n = [cos, sin]`; `print(deg, ch03.linear_strain_rate(G, n))`. *expect:* 0° → 1.000 · 22.5° → 1.414 · 45° → 1.000 ·
    90° → −1.000 (n·S·n = cos 2θ + sin 2θ). *explain:* 1. the probe direction n; 2. n·S·n (the antisymmetric part
    drops out of any n·G·n); 3. the largest stretching, √2 s⁻¹ at 22.5°, is along a principal axis (C12).
22. `nb.check_agree` — **from scratch (curation §7):** track a 1 cm segment at each angle with `ch03.linear_flow_map(G,
    dt)` for dt = 1e-6 s: `mine = (np.linalg.norm(M @ dx0) - np.linalg.norm(dx0))/(np.linalg.norm(dx0)*dt)` →
    `assert np.allclose(mine, ch03.linear_strain_rate(G, n), rtol=1e-5)`; also `ch03.measured_strain_rates(G, n,
    1e-6)["stretch"]` agrees.
23. `nb.figure` — the stretching "rose": measured rate (1/ℓ)Δℓ/Δt from tracked segments at 72 angles (dots) and the
    formula n·S·n (teal curve), plotted against θ ∈ [0°, 180°), with the zero line; principal directions marked by
    dashed verticals (blue stretch 22.5°, rose squeeze 112.5°). *see:* "dots on a smooth cosine-like wave between −√2 and
    +√2." *read:* "the measurement on moving threads and the tensor formula are the same curve: each direction has its
    own stretching rate, and the extremes sit 90° apart." *change:* "…G antisymmetric (solid-body rotation): the curve
    would be flat at zero — nothing stretches."
24. `nb.md` — **What would change if…** "…two threads start perpendicular? Besides stretching they turn toward (or away
    from) each other; the rate at which the right angle closes is the off-diagonal of S (C08)."

**C08 — Shear strain rate; rigid motion ⇒ S = 0**
25. `nb.core("C08", "Shear strain rate: how fast a right angle closes", question="Two dye threads cross at a right
    angle. How fast does that angle close, and why does the book call half of that rate $S_{12}$?")`
26. `nb.md` — **The problem in plain words:** "Push the top of a deck of cards sideways: the corners of each card stop
    being right angles. Fluids flowing past a wall are sheared this way all the time (Couette and boundary-layer flows,
    Ch. 8–9), and the shear *rate* is what the viscous shear stress depends on."
27. `nb.md` — **The idea:**
   ```
        C ┐ dα                 vertical side tilts clockwise by dα   ∝ ∂u₁/∂x₂
          │ ╲                  horizontal side tilts counterclockwise by dβ ∝ ∂u₂/∂x₁
        B └───── dβ            closing rate α + β  →  S₁₂ = ½ D(α+β)/Dt        (Fig. 3.11)
   ```
28. `nb.primer("small-angle approximation", "For an angle ε measured in radians, $\tan\varepsilon\approx\varepsilon$
    and $\cos\varepsilon\approx1$ with errors of order ε³ and ε²: a tiny tilt is its own tangent. Every angle below is
    proportional to dt, so the neglected pieces vanish after dividing by dt.", code="eps = np.array([0.1, 0.01, 0.001])
    \nprint(np.tan(eps) - eps)                  # [3.3e-4 3.3e-7 3.3e-10]: falls like ε³")` — *expect:*
    `[3.35e-04 3.33e-07 3.33e-10]`.
29. `nb.derivation("D09", …)` — Part F D09 (8 steps), ref "".
30. `nb.primer("rigid-body velocity Ω × x", "A rigid body turning at angular velocity Ω about an axis through the origin
    moves each of its points with $\mathbf v=\boldsymbol\Omega\times\mathbf x$ — perpendicular to both the axis and the
    arm, of size Ω times the distance from the axis. Adding a translation U gives every rigid motion:
    $\mathbf u=\mathbf U+\boldsymbol\Omega\times\mathbf x$.", code="Omega = np.array([0.0, 0.0, 1.0])     # 1 rad/s
    about z\nprint(np.cross(Omega, [2.0, 0.0, 0.0]))        # [0. 2. 0.]: 2 m from the axis → 2 m/s, tangential")` —
    *expect:* `[0. 2. 0.]`.
31. `nb.derivation("D10", …)` — Part F D10 (5 steps), ref "".
32. `nb.note` — **N27 [B]** "**Rigid motion does not deform.** For $\mathbf u=\mathbf U+\boldsymbol\Omega\times\mathbf x$
    with U and Ω uniform, S = 0 and R carries all of G (ω = 2Ω). Hence adding a rigid motion — watching from a
    translating or rotating frame — changes G only by an antisymmetric part: **S is the same for every such observer**
    (even if U varies in time); ω is not (N28)."
33. `nb.worked_example("simple shear u = (γy, 0), γ = 1 s⁻¹, for dt = 0.01 s", "1. The top of a 1 cm vertical side moves
    $\gamma\,\delta x_2\,dt=1\times0.01\times0.01=10^{-4}$ m further than its foot: $\tan d\alpha=10^{-4}/10^{-2}=0.01$,
    $d\alpha\approx0.01$ rad. 2. The horizontal side does not tilt: dβ = 0 (u₂ = 0 everywhere). 3.
    $S_{12}=\tfrac12(0.01+0)/0.01=0.5$ s⁻¹ $=\gamma/2$. 4. So γ = 2S₁₂: ⚠️ Ch. 2's Ex. 2.4 called S₁₂ itself Γ; here
    γ is twice that.")`
34. `nb.code` — *code:* `G = ch03.velocity_gradient_preset("simple_shear", Gamma=1.0)` ·
    `print(ch03.shear_strain_rate(G, [1, 0], [0, 1]))` · `print(ch03.measured_strain_rates(G, [1, 0], 1e-6)["closing"]/2)` ·
    rigid motion: `U, Om = np.array([0.3, -0.1, 0.2]), np.array([0.5, -1.0, 2.0])`;
    `Gr = ch03.velocity_gradient_at(lambda x, t: ch03.rigid_body_velocity(U, Om, x), np.array([0.4, 0.1, -0.7]), 0.0)`;
    `print(np.round(ch03.strain_rate_tensor(Gr), 12))`; `print(ch03.vorticity_from_gradient(Gr))`. *expect:* `0.5` ·
    `0.5` (to 1e-6) · a 3 × 3 zero matrix · `[ 1. -2.  4.]` (= 2Ω). *explain:* 1. S₁₂ from the formula n₁·S·n₂; 2. half the
    measured closing rate of two tracked perpendicular threads; 3. a random rigid motion: its gradient has no symmetric
    part and its vorticity is twice its angular velocity.
35. `nb.figure` — the angle between two initially perpendicular material threads vs t ∈ [0, 1] s (from
    `ch03.material_line_angle`) for three flows: simple shear γ = 1 (purple), pure strain along ±45°
    ("irrotational_strain" preset with Γ = 0.5, teal), solid-body rotation (orange); dashed tangent lines of slope
    −2S₁₂ at t = 0. *see:* "the rotation curve stays at 90°, the two others fall — each starting along its dashed
    tangent." *read:* "the initial slope is −2S₁₂ in rad/s: shear and pure strain close the angle at the same initial
    rate (both have S₁₂ = 0.5 s⁻¹); rotation turns both threads together and closes nothing (S = 0)." *change:* "…the
    threads were turned by 45°: in the pure-strain flow their angle would stay 90° (they lie on the principal axes) —
    the off-diagonal of S depends on the axes (C12)."
36. `nb.md` — **What would change if…** "…three edges of a small box stretch at once? Then its volume changes — at the
    sum of the three linear rates (C09)."

**C09 — Volumetric strain rate (3.14)**
37. `nb.core("C09", "Volumetric strain rate: divergence is how fast a blob swells (3.14)", question="A small parcel of
    air rises and expands. At what rate does its volume grow per unit volume, and how is that read off the velocity
    field?")`
38. `nb.md` — **The problem in plain words:** "Rising air expands, sinking air is compressed; in the ocean water hardly
    changes volume at all. The fractional growth rate of a parcel's volume is the kinematic half of mass conservation
    (Ch. 4: $D\rho/Dt=-\rho\nabla\cdot\mathbf u$). We find it from the three stretching rates of C07."
39. `nb.md` — **The idea:**
   ```
   δV = δx₁ δx₂ δx₃   each edge stretches at S₁₁, S₂₂, S₃₃   →   (1/δV) D(δV)/Dt = S₁₁ + S₂₂ + S₃₃ = ∇·u   (3.14)
   shear only tilts the faces (second order);  the trace does not depend on the orientation of the box
   ```
40. `nb.derivation("D11", …)` — Part F D11 (7 steps), ref "3.14".
41. `nb.worked_example("u = (x, y, z) s⁻¹ and a 1 cm³ box for 0.01 s", "1. ∇·u = 1 + 1 + 1 = 3 s⁻¹. 2. First order:
    $\delta V\approx1\times(1+3\times0.01)=1.03$ cm³. 3. Exact for this linear flow (Jacobi's formula, D11 step 7):
    $e^{0.03}=1.03045$ cm³. 4. Simple shear u = (γy, 0): ∇·u = 0 — the sheared box keeps its volume although it
    changes shape.")`
42. `nb.code` — *code:* `G3 = np.eye(3)` · `print(ch03.volumetric_strain_rate(G3), ch03.material_volume_ratio(G3, 0.01))`
    · Jacobi check on a random G: `rng = np.random.default_rng(3)`; `Grand = rng.normal(size=(3, 3))`;
    `print(ch03.material_volume_ratio(Grand, 0.7), np.exp(0.7*np.trace(Grand)))` · rotation invariance:
    `C = ch03.rotation_matrix_2d(0.4)`; `G2 = np.array([[1.0, 2.0], [0.0, -1.0]])`;
    `print(np.trace(ch03.transform_tensor(G2, C)), np.trace(G2))`. *expect:* `3.0 1.0304545` · two equal numbers ·
    `0.0 0.0`. *explain:* 1. the trace of G is the volumetric rate (3.14); 2. the exact finite-time volume factor
    det e^{Gt} equals e^{t tr G} (Jacobi); 3. the trace does not change when the axes turn (Ch. 2 C07).
43. `nb.check_agree` — **from scratch:** a 1 cm cube of 8 corner tracers carried by `linear_flow_map(G3, t)`; volume from
    `np.linalg.det` of its three edge vectors → `assert np.allclose(vol/1e-6, ch03.material_volume_ratio(G3, t))` at
    t = 0.01 and 0.5 s.
44. `nb.figure` — volume ratio δV(t)/δV(0) vs t ∈ [0, 1] s for the expansion G = I (tracked cube, dots; e^{3t}, blue
    curve; first-order 1 + 3t, dashed grey) and for simple shear (flat line at 1, teal). *see:* "an exponential that
    leaves its tangent line; a flat line for shear." *read:* "(3.14) gives the initial slope, 3 s⁻¹ per unit volume;
    over finite time the growth compounds (e^{3t}); shear alone changes shape, never volume." *change:* "…G = −I
    (compression): the volume would decay as e^{−3t} — a sinking, compressed parcel."
45. `nb.explainer("fluid_element_deformation", heading="What does each number in S measure?", why="Drag the four entries
    of G and watch a square element and a ring of tracers deform. The measured rates — stretching per length along a
    probe direction, the closing rate of a right angle, the area growth — land on the formula rates n·S·n, 2S₁₂ and
    tr G, and the rose curve shows the stretching rate in every direction at once.", tries=["Preset 'solid-body
    rotation': every measured rate is zero — S = 0 (N27).", "Preset 'simple shear': rotate the probe to 45° and find the
    fastest stretching.", "Preset 'expansion': the area ratio follows the ghost e^{t tr G}.", "Click a tracer: the
    inspector gives its du = G·dx arithmetic (D07)."])`
46. `nb.md` — **What would change if…** "…you looked not at how the element deforms but at how it turns? The two
    perpendicular threads of C08 turn by −α and +β; their *average* turning is the spin, and it is half the vorticity
    (C10)."
47. `nb.recap("R04", "The rotation tensor is a vector in disguise", "R is antisymmetric: zero diagonal and three
    independent entries, which pair up with the three components of the vorticity $\boldsymbol\omega=\nabla\times\mathbf u$.",
    where="Ch. 2 §2.10, C12")`
48. `nb.recap("R05", "R ↔ ω", "$R_{ij}=-\varepsilon_{ijk}(\nabla\times\mathbf u)_k=-\varepsilon_{ijk}\omega_k=\begin{bmatrix}0&-\omega_3&\omega_2\\\omega_3&0&-\omega_1\\-\omega_2&\omega_1&0\end{bmatrix}$
   *(Eq. 3.15, = Eqs. 2.26–2.27)* — `ch03.antisymmetric_from_vector(omega)`; `ch03.vorticity_from_gradient(G)` returns
   ω = vector(R).", where="Ch. 2 §2.10, D15")`
49. `nb.recap("R06", "The components of the vorticity", "$\omega_1=\frac{\partial u_3}{\partial x_2}-\frac{\partial u_2}{\partial x_3},\ \omega_2=\frac{\partial u_1}{\partial x_3}-\frac{\partial u_3}{\partial x_1},\ \omega_3=\frac{\partial u_2}{\partial x_1}-\frac{\partial u_1}{\partial x_2}$
   *(Eq. 3.16, = Eq. 2.25)* — `ch03.curl` on a grid, `ch03.vorticity(u, x, t)` at a point.", where="Ch. 2 §2.9, D12")`
50. `nb.recap("R07", "Circulation", "$\Gamma\equiv\oint_C\mathbf u\cdot d\mathbf s=\int_A\boldsymbol\omega\cdot\mathbf n\,dA$
   *(Eq. 3.18)* (Stokes' theorem, Eq. 2.34): the circulation round a loop equals the flux of vorticity through any
   surface it bounds, so **vorticity is circulation per unit area** (Eq. 2.35). ⚠️ Γ here is a circulation [m²/s], not
   Ch. 2's shear rate. `ch03.circulation(u, ch03.planar_loop(...))`. Used again in C13 and C14; Kelvin's theorem (Ch. 5)
   and lift (Ch. 6, 14) are built on it.", where="Ch. 2 §2.13, C16")` (the recap text writes (2.34), $\oint_C\mathbf
   u\cdot\mathbf t\,ds=\int_A(\nabla\times\mathbf u)\cdot\mathbf n\,dA$, and (2.35) in full as Ch. 2 did.)

**C10 — Vorticity is twice the spin of a fluid element**
51. `nb.core("C10", "Vorticity is twice the element's spin — for any pair of lines, in any frame but its own",
    question="Drop a tiny paddle wheel into a river where the water near the bank is slower. Does it turn, how fast, and
    does the answer depend on which two sticks of the wheel you watch — or on whether you watch from the spinning
    Earth?")`
52. `nb.md` — **The problem in plain words:** "The vorticity ω = ∇×u is the most important derived field in atmosphere
    and ocean dynamics (cyclones, eddies, potential vorticity). Its meaning is physical: a small paddle wheel carried
    by the flow turns at half of it. We prove that from the two threads of C08 — and find that a flow in perfectly
    straight lines can spin every element, while (C13) a flow in circles need not."
53. `nb.md` — **The idea:**
   ```
   horizontal thread turns at +dβ/dt, vertical thread at −dα/dt (counterclockwise +)
   spin of the element ≡ their average  = ½(−α̇ + β̇) = ½ω₃ = R₂₁/2                (Fig. 3.11)
   individual threads may turn at different rates (shear!) — the average of ANY perpendicular pair is the same
   ```
54. `nb.primer("angular velocity of a line", "A line segment at angle θ (counterclockwise from +x) turns at
    $\dot\theta=d\theta/dt$; counterclockwise is positive. If its tip moves relative to its tail with velocity δu, only
    the part of δu perpendicular to the segment turns it: $\dot\theta=(\mathbf e_\theta\cdot\delta\mathbf u)/\ell$ with
    $\mathbf e_\theta=(-\sin\theta,\cos\theta)$.", code="theta, ell = np.pi/2, 0.01          # a vertical 1 cm
    thread\ndu = np.array([0.01, 0.0])                # tip moves 1 cm/s to the right relative to the tail\nprint(np.array(
    [-np.sin(theta), np.cos(theta)]) @ du/ell)   # -1.0 rad/s: clockwise")` — *expect:* `-1.0`.
55. `nb.derivation("D12", …)` — Part F D12 (8 steps), ref "".
56. `nb.md` — "**The running example: a parallel shear flow.** Near a river bank (or in the wind near the ground) the
    velocity is $\mathbf u=(u_1(x_2),0)$ and locally $u_1\approx\gamma x_2$ with the shear rate $\gamma\equiv du_1/dx_2$.
    By (3.16), $\omega_3=\frac{\partial u_2}{\partial x_1}-\frac{\partial u_1}{\partial x_2}=-\gamma$: straight streamlines, yet
    nonzero vorticity — clockwise for γ > 0. How do its threads turn?" (§3.5 re-states this flow as N34.)
57. `nb.derivation("D14", …)` — Part F D14 (6 steps), ref "".
58. `nb.primer("rotating frame of reference", "An observer on a turntable (or on the Earth) turning at Ω about z sees
    every point that is fixed on the turntable at rest, although it moves with $\boldsymbol\Omega\times\mathbf x$ in
    the lab. At the instant the two sets of axes coincide, velocities are related by
    $\mathbf u_{\rm lab}=\boldsymbol\Omega\times\mathbf x+\mathbf u_{\rm rot}$. Unlike a Galilean frame, this correction
    varies from point to point — so derivatives (and the vorticity) change. Ch. 4 §4.7 develops the full rotating-frame
    equations.", code="Omega = np.array([0, 0, 0.5]); x = np.array([2.0, 0.0, 0.0])\nu_lab = np.array([0.0, 1.0, 0.0])
    # a point on a disc turning at 0.5 rad/s\nprint(u_lab - np.cross(Omega, x))            # [0. 0. 0.]: at rest for
    the rotating observer")` — *expect:* `[0. 0. 0.]`.
59. `nb.derivation("D13", …)` — Part F D13 (6 steps), ref "".
60. `nb.note` — **N28 [B]** "**Vorticity depends on the frame** (unlike S, N27): an observer rotating with the element
    sees ω′ = 0; one rotating at Ω e_z sees $\omega'_z=\omega_z-2\Omega$. **Climate hook:** on the Earth the vorticity
    measured relative to the ground (the relative vorticity ζ) and the one seen from space (the absolute vorticity) differ
    by the local normal component of the planetary vorticity, 2Ω sin φ = f (1.03 × 10⁻⁴ s⁻¹ at 45°N) — larger than ζ of a
    typical mid-latitude cyclone (≈ 10⁻⁵–10⁻⁴ s⁻¹). Ch. 4 §4.7 and Ch. 13 build on this." equation
    `\omega'_z=\omega_z-2\Omega`, no ref.
61. `nb.note` — **N29 [B]** "**Irrotational flow.** A flow with no vorticity anywhere:" equation
    `\boldsymbol\omega=0,\quad\text{or equivalently}\quad R_{ij}=\partial u_i/\partial x_j-\partial u_j/\partial x_i=0`, ref
    "3.17". "Then u can be written as the gradient of a *velocity potential*, $u_i=\partial\phi/\partial x_i$, because
    $\nabla\times\nabla\phi=0$ (Ch. 2). ⚠️ The converse needs a region without holes (simply connected): the line vortex
    of C13 has ω = 0 everywhere except its axis, yet the circulation round the axis is 2πB and the 'potential' φ = Bθ
    jumps by 2πB after one turn. Potential flow is Ch. 6." Followed by `nb.code` — *code:* `xs, ys = sp.symbols('x
    y')`; `print(ch03.potential_velocity(xs**2 - ys**2, [xs, ys]))` · on a 61 × 61 grid excluding the core:
    `phi, jump = ch03.velocity_potential_2d(ch03.vortex_velocity_field("line", Gamma=2*np.pi, sigma=1.0), grid,
    x_ref=(1.5, 0.0))`; `print(jump, 2*np.pi*1.0)`. *expect:* `([2*x, -2*y], 0)` · `6.28… 6.283` (the path dependence
    equals 2πB for B = 1 m²/s). *explain:* "a potential always gives an irrotational field; an irrotational field in a
    region with a hole need not have a single-valued potential."
62. `nb.worked_example("solid-body rotation and shear", "1. Solid body, ω₀ = 1 rad/s: G = [[0, −1], [1, 0]] s⁻¹;
    $\omega_3=1-(-1)=2$ s⁻¹; every thread turns at $\dot\theta=1$ rad/s; spin = ½ω₃ = 1 rad/s — the element spins as
    fast as it revolves. 2. Shear γ = 1 s⁻¹: $\dot\theta=-\sin^2\theta$: the horizontal thread 0, the vertical −1 rad/s,
    a thread at 30° −0.25 rad/s and its partner at 120° −0.75 rad/s; each pair averages −0.5 = ω₃/2. 3. Watch the solid
    body from a turntable at Ω = 1 rad/s: $\omega'_z=2-2\times1=0$ — the tank looks at rest.")`
63. `nb.code` — *code:* `Gsb = ch03.velocity_gradient_preset("solid_body_rotation", Gamma=1.0)` ·
    `Gsh = ch03.velocity_gradient_preset("simple_shear", Gamma=1.0)` · `print(ch03.element_rotation_rate(Gsb),
    ch03.element_rotation_rate(Gsh))` · `th = np.deg2rad([0, 30, 90, 120])`; `print(ch03.material_line_rotation_rate(Gsh,
    th))` · `pairs = np.linspace(0, np.pi, 7)`; `print(0.5*(ch03.material_line_rotation_rate(Gsh, pairs) +
    ch03.material_line_rotation_rate(Gsh, pairs + np.pi/2)))` · `print(ch03.vorticity_in_rotating_frame(2.0, 1.0))` ·
    `print(ch03.parallel_shear_kinematics(1.0)["omega3"])`. *expect:* `[0. 0. 1.] [ 0.   0.  -0.5]` · `[-0. -0.25 -1.
    -0.75]` · seven values −0.5 · `0.0` · `-1.0`. *explain:* 1. spin = ½ω from G; 2. individual threads in the shear flow
    turn at −γ sin²θ; 3. every perpendicular pair averages −γ/2 (D14); 4. the rotating-observer rule (D13); 5. the
    shear flow's ω₃ = −γ.
64. `nb.plotly` — `slider_figure` over γ ∈ [−2, 2] s⁻¹ (21 steps): θ̇(θ) for a single thread (teal curve, θ from 0 to
    180°), the partner curve θ̇(θ + 90°) (teal dashed), and their average (orange, flat); the zero line; title "Single
    threads disagree, every pair agrees: spin = ω₃/2". Followed by `nb.md` see/read/change: "…the two teal curves
    mirror each other about −γ/2 and the orange line never bends; at γ = 0 everything is zero. For solid-body rotation
    (not on this slider) both teal curves would be flat at ω₀ — every thread turns alike."
65. `nb.md` — **What would change if…** "…you split the relative velocity du = G·dx of C06 into the part from S and the
    part from R? The R part is exactly a rigid rotation at ω/2 — C11."

**C11 — Relative velocity = deformation + rigid rotation (3.19)**
66. `nb.core("C11", "Deformation plus rigid rotation: du = S·dx + ½ω × dx (3.19)", question="The arrow from one fluid
    particle to its neighbour's velocity can be split into two arrows. What are they, and why does only one of them
    matter for friction?")`
67. `nb.md` — **The problem in plain words:** "A fluid that spins like a rigid body (a stirred cup after it settles)
    feels no internal friction; a fluid that is sheared does. To write the friction law of Ch. 4 we must separate, near
    every point, the motion that deforms from the motion that merely rotates."
68. `nb.md` — **The idea:**
   ```
   du   =   S·dx          +     ½ ω × dx                 (3.19)
   total    deformation         rigid rotation at angular velocity ω/2 (same as Ω × x with Ω = ω/2)
   purple   teal                orange
   ```
69. `nb.derivation("D15", …)` — Part F D15 (6 steps), ref "3.19".
70. `nb.note` — **N30 [C]** "The second term has the form of the rigid-body velocity $\mathbf v=\boldsymbol\Omega\times\mathbf x$
    (primer in C08) with $\boldsymbol\Omega=\boldsymbol\omega/2$: near any point, the fluid turns rigidly at half its
    vorticity."
71. `nb.worked_example("shear γ = 1 s⁻¹, neighbour dx = (0, 1) m", "1. $d\mathbf u=\mathbf G\cdot d\mathbf x=(\gamma\cdot1,0)=(1,0)$
    m/s. 2. Strain part $\mathbf S\cdot d\mathbf x$ with S = [[0, ½], [½, 0]]: (0.5, 0). 3. Rotation part with
    $\boldsymbol\omega=(0,0,-1)$: $\tfrac12\boldsymbol\omega\times d\mathbf x=\tfrac12(0\cdot0-(-1)\cdot1,\ (-1)\cdot0-0\cdot0,\ 0)=(0.5,0,0)$.
    4. Sum (1, 0) ✓: half of the shear is stretching along 45°, half is a clockwise spin.")`
72. `nb.code` — *code:* `du, du_s, du_r = ch03.relative_velocity_split(Gsh, [0.0, 1.0])`; `print(du, du_s, du_r)` · a
    3-D random G: `dx3 = np.array([0.2, -0.1, 0.3])`; `du, du_s, du_r = ch03.relative_velocity_split(Grand, dx3)`;
    `w = ch03.vorticity_from_gradient(Grand)`; `print(np.allclose(du_r, 0.5*np.cross(w, dx3)), np.dot(du_r, dx3))`.
    *expect:* `[1. 0.] [0.5 0.] [0.5 0.]` · `True 0.0` (the rotation part is perpendicular to dx). *explain:* 1. the
    split for the worked example; 2. in 3-D the rotation part is exactly ½ω × dx and never changes the distance to the
    neighbour.
73. `nb.check_agree` — **from scratch (curation §7):** a triple loop over i, j, k with `eps = ch03.levi_civita()`:
    `rot[i] += -0.5*eps[i, j, k]*w[k]*dx3[j]` → `assert np.allclose(rot, du_r)` (the index form
    $-\tfrac12\varepsilon_{ijk}\omega_kdx_j$ of (3.19), $du_i=\big(S_{ij}-\tfrac12\varepsilon_{ijk}\omega_k\big)dx_j$, written as
    loops); `assert np.allclose(ch03.strain_rate_tensor(Grand) @ dx3 + rot, Grand @ dx3)`.
74. `nb.figure` — the ring of 24 neighbours (radius 0.1 m) in the shear flow γ = 1 with three arrow sets (`ring_arrows`):
    total du (purple), strain part S·dx (teal, pointing out along 45° and in along 135°), rotation part ½ω × dx (orange,
    tangential clockwise); legend in the title. *see:* "purple arrows that are all horizontal; teal arrows that fan out
    and in; orange arrows that circulate." *read:* "the shear's horizontal arrows are the sum of a pure strain and a
    rigid clockwise spin at 0.5 rad/s — the ring would become an ellipse and turn." *change:* "…γ = 0 and G =
    [[0, −1], [1, 0]] (solid body): the teal arrows vanish; the purple arrows are the orange ones."
75. `nb.note` — **N33 [C]** "**Summary of §3.4.** The relative motion near a point = a rigid rotation of the element (at
    ω/2) + a deformation (S), which itself = stretching along three perpendicular principal axes (next block)."
76. `nb.md` — **What would change if…** "…you turned the axes to where S has no off-diagonal entries? Then the
    deformation is three pure stretchings — a small sphere becomes an ellipsoid (C12)."

**C12 — Principal strain axes: sphere → ellipsoid**
77. `nb.core("C12", "Principal strain axes: a small sphere becomes an ellipsoid", question="A round drop of dye is
    released into a flow. What shape is it a moment later, and which way do its long and short axes point?")`
78. `nb.md` — **The problem in plain words:** "Satellite images of ocean colour show round patches of plankton pulled into
    ellipses and then filaments; weather fronts form where the wind's deformation squeezes temperature contours together
    (frontogenesis, Ch. 13). Both are the principal axes of the strain-rate tensor at work."
79. `nb.md` — **The idea:**
   ```
   in the principal frame (overbar) S is diagonal:  dū₁ = S̄₁₁ dx̄₁,  dū₂ = S̄₂₂ dx̄₂,  dū₃ = S̄₃₃ dx̄₃     (3.21)
   each principal direction stretches in proportion to its own length  →  circle → ellipse on those axes (Fig. 3.13)
   M (on the stretching axis) moves most, N (on the other) not at all in that direction
   ```
80. `nb.note` ×2 — **N31 [B]** "In the frame of the principal axes of S (Ch. 2 C13: a symmetric tensor has three
    perpendicular eigenvectors; ⚠️ the book's '§2.12' reference means §2.11), the strain part of (3.19) is diagonal:"
    equation `d\bar{\mathbf u}=\bar{\mathbf S}\cdot d\bar{\mathbf x}=\begin{bmatrix}\bar S_{11}&0&0\\0&\bar S_{22}&0\\0&0&\bar S_{33}\end{bmatrix}\begin{bmatrix}d\bar x_1\\d\bar x_2\\d\bar x_3\end{bmatrix}`,
    ref "3.20". · **N32 [B]** "Its three components, with $\bar S_{\alpha\alpha}$ the eigenvalues of S:" equation
    `d\bar u_1=\bar S_{11}d\bar x_1,\quad d\bar u_2=\bar S_{22}d\bar x_2,\quad d\bar u_3=\bar S_{33}d\bar x_3`, ref "3.21".
81. `nb.primer("linear map of a circle is an ellipse", "Multiplying every point of a unit circle by a matrix M gives an
    ellipse. For a symmetric M its axes are M's eigenvectors and its semi-axes the eigenvalues; for a general M they are
    the *singular values* and singular vectors (`np.linalg.svd`) — a gloss we need only for the finite-time caveat.",
    code="M = np.array([[1.1, 0.0], [0.0, 0.9]])       # stretch x by 10 %, squeeze y by 10 %\ns = np.linspace(0, 2*np.pi,
    400); P = M @ np.stack([np.cos(s), np.sin(s)])\nprint(P[0].max(), P[1].max())                  # 1.1 0.9: the
    semi-axes")` — *expect:* `1.1 0.9` (to 4 decimals).
82. `nb.derivation("D16", …)` — Part F D16 (8 steps), ref "".
83. `nb.worked_example("shear γ = 1 s⁻¹, a 1 mm drop after 0.1 s", "1. Principal rates ±γ/2 = ±0.5 s⁻¹ at +45° and
    −45° (Ch. 2 Ex. 2.4, now with γ = 2S₁₂). 2. First order (D16): semi-axes $1\times(1\pm0.5\times0.1)=1.05$ and 0.95
    mm. 3. The strain acting alone: $e^{\pm0.05}=1.0513$ and 0.9512 mm. 4. Strain + rotation together (the true shear
    flow): singular values of $e^{\mathbf Gt}$ = [[1, 0.1], [0, 1]]: 1.0513 and 0.9513 mm, long axis at 43.6° — the
    rotation part (clockwise 0.5 rad/s × 0.1 s ≈ 2.9°) has already turned it a little off 45°.")`
84. `nb.code` — *code:* `lam, axes = ch03.principal_strain_rates(Gsh)`; `print(lam, axes[:, 1])` · for `m in
    ("first_order", "strain_only", "exact")`: `a, d = ch03.strain_ellipse_axes(Gsh, 0.1, method=m)`;
    `print(m, np.round(a, 4), np.round(np.degrees(np.arctan2(d[1, 0], d[0, 0])), 1))` · `du_bar, dx_bar =
    ch03.strain_velocity_principal(Gsh, [0.001, 0.0])`; `print(du_bar/dx_bar)` (componentwise where dx̄ ≠ 0).
    *expect:* `[-0.5  0.5] [0.7071 0.7071]` · first_order [1.05 0.95] 45.0 · strain_only [1.0513 0.9512] 45.0 · exact
    [1.0513 0.9513] 43.6 · `[-0.5  0.5]` (in the eigenframe ordered as λ). *explain:* 1. eigenvalues and eigenvectors of
    S; 2. three answers to "what shape after 0.1 s": the book's first-order statement, the strain alone, and the exact
    linear flow; 3. in the eigenframe each velocity component is its own eigenvalue × its own coordinate (3.21).
85. `nb.animation` — `player="video"`, 48 frames (FAST 24), t ∈ [0, 3] s in the shear flow γ = 1: a circle of 72 tracers
    (`deform_circle`) becoming an ellipse; the principal axes of S drawn fixed at ±45° (blue stretch, rose squeeze); the
    ellipse's current long axis (purple) from the SVD; readouts "semi-axes a, b" and "long axis at …°". Followed by
    `nb.md` see/read/change: "…at first the ellipse grows along the blue 45° axis, as D16 says; as time goes on the
    purple axis tips toward the flow direction (the rotation part keeps turning it) — the principal-axis statement is
    about the *first instant*. In pure strain ('irrotational_strain' preset) the purple and blue axes would stay together
    for ever."
86. `nb.explainer("spin_and_principal_axes", heading="Can a straight flow make a fluid element spin?", why="Two
    perpendicular threads and a small paddle wheel ride in a flow you choose. Drag the pair's angle: the two threads turn
    at different rates but their average stays at ½ω₃. Drag a probe round the ring: its relative velocity splits into a
    teal strain arrow and an orange rotation arrow, (3.19) $du_i=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i$.
    The circle turns into an ellipse on the principal axes.", tries=["Preset 'parallel shear γ = 1': sweep the pair
    angle from 0 to 180° and watch the average line stay flat at −0.5 rad/s.", "Preset 'rotate with the element'
    (Ω = ω/2): the paddle wheel stops — ω′ = 0.", "Preset 'pure strain': the orange arrows vanish; the ellipse axes never
    move.", "Derivation tab, D15: at the ε-swap step the orange arrow flips sign on screen."])`
87. `nb.md` — **What would change if…** "…the streamlines were circles instead of straight lines? Then 'going round'
    and 'spinning' come apart completely — §3.5's two vortices."

### A.5 §3.5 Kinematics of Simple Plane Flows — N34, N35 (head), C13 (+N36–N41 · D17), C14 (+N42–N44 · D18, D19, D20 · E6)
1. `nb.section("3.5", "Kinematics of Simple Plane Flows", intro="**What is this section about?** Two families of plane
   flow where one coordinate does all the work: the parallel shear flow (straight streamlines, spinning elements) and
   circular flows — solid-body rotation, the irrotational vortex and the realistic vortices between them.")`
2. `nb.note` — **N34 [B]** "**Parallel shear flow** $\mathbf u=(u_1(x_2),0)$, $\gamma(x_2)\equiv du_1/dx_2$ — the flow
   you used in C10. By (3.16), $\omega_3=\frac{\partial u_2}{\partial x_1}-\frac{\partial u_1}{\partial x_2}=-\gamma$. The
   vertical thread AB turns at −γ, the horizontal BC at 0, their average −γ/2 = ω₃/2 — and D14 (in C10) showed that
   *every* perpendicular pair gives −γ/2. ⚠️ ω₃ = −γ is clockwise for γ > 0. Couette flow (Ch. 8) and boundary layers
   (Ch. 9) are this flow." equation `\gamma(x_2)\equiv du_1/dx_2,\qquad \omega_3=-\gamma`, no ref.
3. `nb.note` — **N35 [B]** "**Two elements in the same shear (Fig. 3.14).** For the square ABCD with sides along the axes,
   S has only off-diagonal entries — it shears without stretching its sides. Turned by 45°, onto the principal axes
   (⚠️ Ch. 2's Ex. 2.4 had S₁₂ = Γ; here S₁₂ = γ/2), S is diagonal — the square PQRS stretches along x̄₁ and is
   compressed along x̄₂ (eigenvalue −γ/2, i.e. compression at rate γ/2) while its corners stay right angles. Both spin at
   −γ/2." equation `S_{ij}=\begin{bmatrix}0&\gamma/2\\\gamma/2&0\end{bmatrix},\qquad \bar S_{ij}=\begin{bmatrix}\gamma/2&0\\0&-\gamma/2\end{bmatrix}`,
   no ref.
4. `nb.animation` — `player="frames"` (16 frames, FAST 10; `shear_elements_frames(1.0, times)` from
   `scripts/ch03_fig3_14_shear_elements.py`): γ = 1 s⁻¹, t = 0…0.75 s; left: square ABCD with its four side lengths and
   corner angle printed per frame; right: the 45° square PQRS with its side lengths and corner angle; a small paddle in
   each, both turned by −γt/2. Followed by `nb.md` see/read/change: "…ABCD's vertical sides lengthen only at second
   order while its corner angle falls at γ per second; PQRS's corner stays 90.0° to first order while PS grows and PQ
   shrinks; both paddles show the same clockwise turn. Step frame by frame: the numbers are the linear (C07), shear (C08)
   and spin (C10) rates made visible. Over long times PQRS's angles drift from 90° too (the finite-time effect of C12's
   animation)."

**C13 — Vorticity in polar coordinates: solid-body rotation vs the irrotational vortex (3.23)**
5. `nb.core("C13", "Going round is not spinning: vorticity in polar coordinates (3.23)", question="Two tanks of water
   both go round in circles — one spun up like a merry-go-round, one draining through a plughole. Put a small paddle
   wheel in each. Which one turns?")`
6. `nb.md` — **The problem in plain words:** "Hurricanes, bathtub drains, the flow round a stirred cup: circular
   streamlines everywhere. The vorticity tells us whether each little parcel spins as it goes round. We need it in polar
   coordinates, where these flows are simple."
7. `nb.md` — **The idea:**
   ```
   solid body  u_θ = ω₀ r :   Γ grows like r²  → ω_z = 2ω₀ everywhere → the wheel spins as it orbits (Fig. 3.15)
   line vortex u_θ = B/r  :   Γ = 2πB for every circle → ω_z = 0 except on the axis → the wheel keeps its heading (Fig. 3.16)
   vorticity = circulation per unit area of a small polar sector  →  ω_z = (1/r)∂(r u_θ)/∂r − (1/r)∂u_r/∂θ   (3.23)
   ```
8. `nb.primer("polar coordinates as a moving basis", "In the plane, $\mathbf e_r=(\cos\theta,\sin\theta)$ and
   $\mathbf e_\theta=(-\sin\theta,\cos\theta)$ change direction with θ. A small polar sector has area $r\,dr\,d\theta$
   (an arc of length r dθ times a width dr); along a circle of radius r the line element is
   $d\mathbf s=r\,d\theta\,\mathbf e_\theta$, along a ray $d\mathbf s=dr\,\mathbf e_r$. Arcs at r and r + dr have
   *different* lengths — the source of the extra 1/r terms in polar formulas.", code="r, dr, dth = 2.0, 0.01, 0.02
   \ninner, outer = r*dth, (r + dr)*dth\nprint(outer - inner, dr*dth)            # 0.0002 0.0002: the arcs differ by
   dr·dθ\nprint(r*dr*dth)                          # 0.0004 m²: the sector's area")` — *expect:* `0.0002 0.0002` ·
   `0.0004`.
9. `nb.derivation("D17", …)` — Part F D17 (10 steps), ref "3.23".
10. `nb.note` — **N36 [B]** "**Solid-body rotation** — a tank spun steadily until the fluid turns with it:" equation
    `u_r=0\quad\text{and}\quad u_\theta=\omega_0r`, ref "3.22". "By (3.23), $\omega_z=\frac1r\frac{d}{dr}(\omega_0r^2)=2\omega_0$
    everywhere: each element spins about its own centre at the rate it revolves round the axis; S = 0, nothing deforms
    (Fig. 3.15)."
11. `nb.note` — **N37 [B]** "Its circulation round a centred circle:" equation
    `\Gamma=\oint_C\mathbf u\cdot d\mathbf s=\int_0^{2\pi}u_\theta r\,d\theta=2\pi ru_\theta=2\pi r^2\omega_0`, ref "3.24".
    "= vorticity 2ω₀ × area πr². It holds for *any* circuit, centred or not: with uniform ω, Stokes (3.18),
    $\Gamma=\oint_C\mathbf u\cdot d\mathbf s=\int_A\boldsymbol\omega\cdot\mathbf n\,dA$, gives Γ = ω × enclosed area. Number:
    ω₀ = 1 s⁻¹, r = 1 m → Γ = 2π = 6.283 m²/s; an off-centre circle of radius 1 m centred at (2, 0) m also gives 6.283
    m²/s (checked below)."
12. `nb.note` — **N38 [B]** "**The irrotational (line) vortex** — the ideal limit of a drain or a tornado far from its
    core:" equation `u_r=0\quad\text{and}\quad u_\theta=B/r`, ref "3.25". "(⚠️ Fig. 3.16 writes C for B.) By (3.23),
    $\omega_z=\frac1r\frac{d}{dr}(r\cdot B/r)=\frac1r\frac{dB}{dr}=0$ for every r > 0."
13. `nb.note` — **N39 [B]** "Yet the circulation round any centred circle is the same nonzero number:" equation
    `\Gamma=\int_0^{2\pi}u_\theta r\,d\theta=2\pi ru_\theta=2\pi B`, ref "3.26". "Independent of r (u_θ r = B is
    constant)."
14. `nb.note` — **N40 [B]** "So all the vorticity sits on the axis. Taking vorticity as circulation per unit area in a
    shrinking disc:" equation `[\omega_z]_{r\to0}=\lim_{r\to0}\frac1A\int_A\omega_z\,dA=\lim_{r\to0}\frac{1}{\pi r^2}\oint_C\mathbf u\cdot d\mathbf s=\lim_{r\to0}\frac{2B}{r^2}`,
    ref "3.27". "— infinite at r = 0 with a finite area integral 2πB: a *delta function* (Ch. 2 gloss: an infinitely
    concentrated amount whose total is finite). The point vortex of Ch. 6 is exactly this."
15. `nb.note` — **N41 [B]** "A loop that does **not** enclose the axis has zero circulation — the sector ABCD of Fig.
    3.16: the radial legs BC and DA contribute nothing (u ⟂ ds), and the arcs cancel because u_θ r = B:" equation
    `\Gamma_{ABCD}=-[u_\theta r]_r\Delta\theta+[u_\theta r]_{r+\Delta r}\Delta\theta=0`, no ref. "Elements in this flow
    deform but do not spin."
16. `nb.worked_example("the two tanks with easy numbers", "1. Solid body ω₀ = 1 s⁻¹: $u_\theta=r$;
    $\omega_z=\frac1r\frac{d(r\cdot r)}{dr}=\frac{2r}{r}=2$ s⁻¹; spin ½ω_z = 1 rad/s = orbit rate u_θ/r = 1 rad/s. 2. Line
    vortex B = 1 m²/s: $u_\theta=1/r$; $\omega_z=\frac1r\frac{d(r\cdot1/r)}{dr}=\frac1r\frac{d(1)}{dr}=0$; spin 0, orbit
    rate 1/r² rad/s. 3. Circulations: circle r = 1 m: solid body 2π(1)²(1) = 6.283 m²/s; line vortex 2π(1) = 6.283 m²/s
    — equal at r = 1, but at r = 2 m: 25.13 vs 6.283 m²/s. 4. Mean vorticity in the disc r = 0.1 m for the line vortex:
    2B/r² = 200 s⁻¹; at r = 0.01 m: 20 000 s⁻¹.")`
17. `nb.code` — *code:* symbolic (3.23): `r, th, w0, B = sp.symbols('r theta omega_0 B', positive=True)`;
    `print(ch03.polar_vorticity_z_sym(0, w0*r, r, th), ch03.polar_vorticity_z_sym(0, B/r, r, th))` · numerical:
    `print(ch03.polar_vorticity_z(lambda r, t: 0.0, lambda r, t: 1.0*r, 1.3, 0.4), ch03.polar_vorticity_z(lambda r, t:
    0.0, lambda r, t: 1.0/r, 1.3, 0.4))` · circulations: `print(ch03.circulation_circle("solid", 1.0, Gamma=2*np.pi,
    sigma=1.0), ch03.circulation_circle("solid", 1.0, center=(2.0, 0.0), Gamma=2*np.pi, sigma=1.0))` ·
    `print([ch03.circulation_circle("line", r_, Gamma=2*np.pi, sigma=1.0) for r_ in (0.5, 1.0, 2.0)])` ·
    `print(ch03.annular_sector_circulation("line", 1.0, 0.1, 0.2, Gamma=2*np.pi, sigma=1.0))` ·
    `print([ch03.mean_vorticity_in_disc("line", r_, Gamma=2*np.pi, sigma=1.0) for r_ in (0.1, 0.01)])`. *expect:*
    `2*omega_0 0` · `2.0 0.0` (to 1e-8) · `6.2832 6.2832` · `[6.2832, 6.2832, 6.2832]` · `0.0` (to 1e-12) · `[200.0,
    20000.0]`. (`"solid"` with Γ = 2π, σ = 1 m means ω₀ = Γ/(2πσ²) = 1 s⁻¹; `"line"` means B = Γ/2π = 1 m²/s — the
    `vortex_profile` convention.) *explain:* 1. sympy applies (3.23) to both profiles; 2. stencils in r and θ agree; 3.
    Γ = ω × area even off-centre (N37); 4. 2πB at every radius (N39); 5. zero for a sector away from the axis (N41);
    6. the δ-core blow-up (N40).
18. `nb.check_agree` — **from scratch (curation §7):** the four legs of a small polar sector (r = 1.3, dr = 0.01,
    dθ = 0.02) for the Rankine vortex inside its core (Γ = 2π, σ = 2 m): outer arc `u_θ(r+dr)·(r+dr)·dθ`, inner arc
    `−u_θ(r)·r·dθ`, radial legs 0; divide by the area r dr dθ → `assert np.allclose(mine, ch03.polar_vorticity_z(…),
    rtol=1e-3)` and vs the Cartesian curl `ch03.vorticity(ch03.vortex_velocity_field("rankine", Gamma=2*np.pi,
    sigma=2.0), p, 0.0)[2]` (both 2π/(π·4) = 0.5 s⁻¹). Markdown: "D17 with numbers: four line integrals and one division."
19. `nb.figure` — **N40** log–log: mean vorticity in a centred disc vs its radius r ∈ [10⁻³, 10] m for the line vortex
    (slope −2, orange), solid body (flat at 2ω₀, teal) and a Rankine vortex (flat inside σ = 1 m, then falling as r⁻²,
    purple). *see:* "a straight line of slope −2, a horizontal line, and a curve that switches from one to the other at
    r = σ." *read:* "an irrotational vortex keeps all its vorticity on the axis — the smaller the disc, the larger the
    average; a real vortex has a finite core, inside which it behaves like a solid body." *change:* "…σ shrank toward
    zero: the purple curve's flat part would climb (Γ/πσ²) and the Rankine vortex would approach the line vortex."
20. `nb.animation` — `player="video"`, 60 frames (FAST 30): two panels, (left) solid-body rotation ω₀ = 1 s⁻¹, (right) line
    vortex B = 1 m²/s, each with four paddle wheels at r = 0.5, 1, 1.5, 2 m carried round (`ch03.pathline` on
    `vortex_velocity_field`) and turned at ½ω_z, plus the polar element ABCD of Figs. 3.15/3.16 tracked by its corners.
    Followed by `nb.md` see/read/change: "…on the left every wheel turns once per revolution and the element keeps its
    shape; on the right the wheels orbit at very different speeds (fast near the axis) but every one keeps pointing the
    same way, while the element is sheared into a thin wedge. Add a paddle at r = 0.2 m on the right and it would race
    round — and still not turn."
21. `nb.md` — **What would change if…** "…a real drain vortex were measured? Near the axis the speed does not blow up
    like 1/r; the core turns like a solid body and the outside like a line vortex — the Rankine and Gaussian models of
    C14."

**C14 — Vortices with a core: Rankine (3.28) and Gaussian (3.29)**
22. `nb.core("C14", "Real vortices have a core: the Rankine and Gaussian models (3.28)–(3.29)", question="A tornado's
    wind is calm at its centre, fastest a little way out and weaker far away. Which simple formula reproduces that
    profile, and where exactly is the fastest wind?")`
23. `nb.md` — **The problem in plain words:** "Tornado chasers, hurricane forecasters and aircraft-wake engineers all
    quote two numbers: the peak wind and the radius where it occurs. Pure solid-body rotation grows for ever; the line
    vortex blows up at the axis. Joining them gives a vortex with a spinning core and an irrotational skirt — the model
    behind every cyclone schematic."
24. `nb.md` — **The idea:**
   ```
                  core (r < σ): solid body, ω_z = Γ/πσ²       outside: irrotational, u_θ = Γ/2πr
   Rankine:   sharp join at r = σ (speed continuous, vorticity jumps)     → peak at r = σ exactly
   Gaussian:  smooth join, ω_z = (Γ/πσ²) e^{−r²/σ²}                        → peak at r ≈ 1.1209 σ
   ```
25. `nb.note` — **N42 [C]** "Real vortices — bathtub drains, wing-tip vortices (Ch. 14), tornadoes, tropical cyclones
    (Ch. 13) — combine a nearly solid-body core with a nearly irrotational outer flow and have bounded speeds."
26. `nb.md` — **The maths — the Rankine vortex** "Uniform vorticity inside a core of radius σ, none outside; Γ is the
    total circulation:
    $$\omega_z(r)=\begin{Bmatrix}\Gamma/\pi\sigma^2=\text{const.}&r\le\sigma\\0&r>\sigma\end{Bmatrix}\quad\text{and}\quad u_\theta(r)=\begin{Bmatrix}(\Gamma/2\pi\sigma^2)r&r\le\sigma\\\Gamma/2\pi r&r>\sigma\end{Bmatrix}\qquad(3.28)$$
    A piecewise formula is evaluated with `np.where(r <= sigma, inside, outside)` (Ch. 1 P46)."
27. `nb.derivation("D18", …)` — Part F D18 (5 steps), ref "3.28".
28. `nb.primer("substitution in an integral", "Replace a messy variable by a simpler one and convert dx too: with
    $s=r^2/\sigma^2$, $ds=2r\,dr/\sigma^2$, so $2\pi r\,dr=\pi\sigma^2\,ds$ and the limits r = 0 … R become s = 0 …
    R²/σ². The value of the integral does not change.", code="from scipy.integrate import quad\nsig, R = 1.0, 1.5
    \nlhs = quad(lambda r: np.exp(-r**2/sig**2)*2*np.pi*r, 0, R)[0]      # in r\nrhs = quad(lambda s: np.exp(-s)*np.pi*sig**2,
    0, R**2/sig**2)[0]   # in s\nprint(lhs, rhs)                                # 2.8322 2.8322")` — *expect:* `2.8322
    2.8322` (= π(1 − e^{−2.25})).
29. `nb.derivation("D19", …)` — Part F D19 (7 steps), ref "3.29".
30. `nb.note` — **N43 [B]** "**The Gaussian vortex** — the smooth version:" equation
    `\omega_z(r)=\frac{\Gamma}{\pi\sigma^2}\exp\!\big(-r^2/\sigma^2\big)\quad\text{and}\quad u_\theta(r)=\frac{\Gamma}{2\pi r}\Big(1-\exp\!\big(-r^2/\sigma^2\big)\Big)`,
    ref "3.29". "It is the Lamb–Oseen vortex of viscous flow at one instant, with $\sigma^2=4\nu t$ growing in time (Ch. 5,
    Ch. 8) — our pointer, not the book's."
31. `nb.primer("np.expm1 and cancellation near zero", "Near r = 0, $1-e^{-x}$ subtracts two numbers that are almost equal
    and loses digits (catastrophic cancellation). `np.expm1(y)` computes $e^y-1$ accurately for tiny y, so
    $1-e^{-x}=$ `-np.expm1(-x)`. `gaussian_vortex` uses it and returns exactly 0 at r = 0.", code="x = 1e-12\nprint(1 -
    np.exp(-x), -np.expm1(-x))       # 1.000088900582341e-12 vs 9.999999999995e-13")` — *expect:* the first number wrong in
    the 5th digit, the second correct.
32. `nb.primer("scipy.optimize.brentq", "Finds a root of f(x) = 0 inside a bracket [a, b] where f changes sign;
    guaranteed to converge and fast. Choose the bracket so it excludes roots you do not want.", code="from scipy.optimize
    import brentq\ng = lambda x: 1 + 2*x - np.exp(x)             # the equation of D20\nprint(g(0.5), g(3.0))
    # 0.351 > 0, −13.09 < 0: a sign change\nprint(brentq(g, 0.5, 3.0))                   # 1.2564312086")` —
    *expect:* `0.3513 -13.0855` · `1.2564312086…`.
33. `nb.derivation("D20", …)` — Part F D20 (7 steps), ref "". Lambert-W gloss in its *check* (Part F).
34. `nb.note` — **N44 [B]** "**Where the wind peaks.** Rankine: at r = σ exactly, $u_{\max}=\Gamma/2\pi\sigma$. Gaussian: where"
    equation `1+2\,\frac{r^2}{\sigma^2}=\exp\!\big(r^2/\sigma^2\big)`, no ref, "i.e. $r\approx1.1209\,\sigma$, with
    $u_{\max}\approx0.6382\,\Gamma/(2\pi\sigma)$ — 36 % below the Rankine peak for the same Γ and σ."
35. `nb.worked_example("Γ = 2π m²/s, σ = 1 m", "1. Rankine inside: $u_\theta=(2\pi/2\pi\cdot1)r=r$ m/s; $\omega_z=2\pi/\pi=2$
    s⁻¹. 2. Rankine peak at r = 1 m: 1 m/s; at r = 2 m: $2\pi/(2\pi\cdot2)=0.5$ m/s. 3. Gaussian at r = 1 m:
    $u_\theta=(1/1)(1-e^{-1})=0.632$ m/s. 4. Gaussian peak: $r=1.1209$ m, $x=r^2=1.2564$,
    $u_\theta=(1-e^{-1.2564})/1.1209=0.7153/1.1209=0.638$ m/s. 5. Far away both → $1/r$: at 5 m, 0.200 m/s
    (Gaussian 0.19999998).")`
36. `nb.code` — *code:* `Gam, sig = 2*np.pi, 1.0` · `r = np.array([0.5, 1.0, 1.1209, 2.0, 5.0])` ·
    `print(ch03.rankine_vortex(r, Gam, sig))` · `print(ch03.gaussian_vortex(r, Gam, sig))` ·
    `rstar = ch03.gaussian_vortex_max_radius(sig)`; `print(rstar, ch03.gaussian_vortex(rstar, Gam, sig)[0])` ·
    `print(ch03.gaussian_vortex_max_radius(sig, method="lambertw"))` · `print(ch03.circulation_circle("gaussian", 3.0,
    Gamma=Gam, sigma=sig)/Gam)` · a small real-vortex table (our order-of-magnitude numbers; Γ from the peak wind via
    $\Gamma=2\pi\sigma u_{\max}/0.6382$): bathtub σ = 5 mm, u_max 0.3 m/s → Γ ≈ 0.015 m²/s; tornado σ = 100 m, 60 m/s →
    5.9 × 10⁴ m²/s; tropical cyclone σ = 36 km, 50 m/s → 1.8 × 10⁷ m²/s (radius of maximum wind 1.1209σ ≈ 40 km).
    *expect:* Rankine u_θ = (0.5, 1.0, 0.8921, 0.5, 0.2), ω = (2, 2, 0, 0, 0) · Gaussian u_θ = (0.4424, 0.6321, 0.6382,
    0.4908, 0.2000), ω = (1.5576, 0.7358, 0.5694, 0.0366, 2.8e-11) · `1.1209064228 0.6381726` · the same radius · 0.99988
    (= 1 − e^{−9}). *explain:* 1. both models at five radii; 2. the maximum by `brentq` (D20) and its Lambert-W twin; 3.
    Γ(r)/Γ → 1 outside the core (D19 step 3); 4. real vortices span nine orders of magnitude in Γ with the same shape.
37. `nb.figure` — two panels vs r/σ ∈ [0, 4]: (a) u_θ for Rankine (teal, kink at 1) and Gaussian (teal dashed), line
    vortex and solid-body ghosts (grey), the peak markers at 1 and 1.1209; (b) ω_z for both (orange, step vs bell).
    Units: u_θ in Γ/(2πσ), ω_z in Γ/(πσ²). *see:* "two humps with peaks at 1 and ≈ 1.12; a step and a bell." *read:*
    "the core is where the vorticity is; outside the core both profiles follow the line vortex; the smooth vorticity
    moves the peak outward and lowers it." *change:* "…σ doubled at the same Γ: in physical units the peak halves and
    moves twice as far out (u_max ∝ Γ/σ)."
38. `nb.plotly` — `slider_figure` over σ ∈ [0.25, 2] m (15 steps) at Γ = 2π m²/s: u_θ(r) for Rankine and Gaussian and
    their peak markers (3 traces); fixed axes r ∈ [0, 5] m, u ∈ [0, 4] m/s. Title "A smaller core spins faster: u_max ∝
    Γ/σ". `nb.md` see/read/change: "…the Rankine peak sits on the kink at σ, the Gaussian peak at 1.12σ, both moving
    outward and down as σ grows; outside about 2σ the two curves coincide with Γ/2πr."
39. `nb.explainer("vortex_paddle_wheels", heading="Going round in circles ≠ spinning", why="Four vortices on one stage:
    paddle wheels ride with the flow and turn at half the local vorticity, a loop you can drag measures the
    circulation, and the profiles u_θ(r) and ω_z(r) show where the core is. The same clock runs the tracers, the
    wheels and the circulation graph.", tries=["Mode 'line vortex': the wheels orbit but always point the same way;
    drag the loop off the axis — Γ drops to 0.", "Mode 'solid body': every wheel turns once per revolution.", "Mode
    'Gaussian': drag σ and watch the peak marker sit at 1.1209σ.", "Preset 'tropical cyclone': read the radius of
    maximum wind and the core vorticity in Explain."])`
40. `nb.md` — **What would change if…** "…the fluid were viscous? The Gaussian core would spread as $\sigma^2=4\nu t$ while
    Γ stays fixed (Ch. 5). And all of this has been about *points* — §3.6 asks how a quantity inside a whole moving
    volume changes."

### A.6 §3.6 Reynolds Transport Theorem — C15 (+N45–N55 · D21, D22, D23, D24 · E7)
1. `nb.section("3.6", "Reynolds Transport Theorem", intro="**What is this section about?** How fast does the amount of
   something inside a volume change when the volume itself moves and deforms? In one dimension this is Leibniz's rule;
   in three it is the Reynolds transport theorem — the bridge from 'what happens to a moving lump of fluid' to equations
   at fixed points, used for every conservation law in Ch. 4.")`
2. `nb.core("C15", "The Reynolds transport theorem: change inside a moving volume (3.35)", question="A balloon is being
   inflated in a room that is warming up. How fast does the heat content inside the balloon change — and which part of
   that comes from the warming, which from the balloon's growing skin?")`
3. `nb.md` — **The problem in plain words:** "Conservation laws are about *things* — a lump of fluid keeps its mass,
   Newton's law acts on it — but the lump moves and deforms. Engineers draw a control volume (a pipe section, a jet
   engine, a layer of ocean) whose walls may move. In both cases we need d/dt of an integral whose region moves: the
   answer has an 'inside' part and a 'swept by the walls' part."
4. `nb.md` — **The idea:**
   ```
   d/dt ∫_{V*(t)} F dV  =  ∫_{V*} ∂F/∂t dV   +   ∮_{A*} F b·n dA           (3.35)
      change of the total    change in place       what the moving wall sweeps in (b·n > 0) or out (b·n < 0)
   1-D:  d/dt ∫_a^b F dx = ∫_a^b ∂F/∂t dx + ḃ F(b) − ȧ F(a)                  (3.30)
   ```
5. `nb.note` — **N45 [C]** "Why we need it: every integral conservation law of Ch. 4 (mass, momentum, energy) is a time
   derivative of an integral over a moving, deforming volume."
6. `nb.primer("differentiation under the integral sign", "If the limits are fixed, the time derivative may pass inside:
   $\frac{d}{dt}\int_a^bF(x,t)\,dx=\int_a^b\frac{\partial F}{\partial t}dx$ (F and ∂F/∂t continuous). It is the special
   case of Leibniz's rule with ȧ = ḃ = 0.", code="t = sp.symbols('t'); x = sp.symbols('x')\nF = sp.sin(x*t)
   \nprint(sp.simplify(sp.diff(sp.integrate(F, (x, 0, 1)), t) - sp.integrate(sp.diff(F, t), (x, 0, 1))))   # 0")` —
   *expect:* `0`.
7. `nb.note` — **N46 [B]** "**Leibniz's theorem** — d/dt of an integral with moving limits:" equation
   `\frac{d}{dt}\int_{x=a(t)}^{x=b(t)}F(x,t)\,dx=\int_a^b\frac{\partial F}{\partial t}dx+\frac{db}{dt}F(b,t)-\frac{da}{dt}F(a,t)`,
   ref "3.30". "The book cites a proof it does not give; here it is."
8. `nb.derivation("D21", …)` — Part F D21 (7 steps), ref "3.30".
9. `nb.note` — **N47 [B]** "**Fig. 3.17, the picture of (3.30):** three thin strips — a band of height ∂F/∂t·dt over
   [a, b] (the interior change), a strip of width db and height F(b) gained at the upper end, a strip of width da and
   height F(a) lost at the lower end; the corner pieces are dt² small." followed by `nb.animation` — `player="frames"`
   (12 frames, FAST 8; `leibniz_strips_figure("warming", …)`), F = 1 + 0.5x + 0.3t, a(t) = 1 − 0.2t, b(t) = 3 + 0.4t:
   the area under F between the moving limits, the three strips shaded per frame (blue interior, orange at b, rose at
   a), running totals printed. `nb.md` see/read/change: "…the orange strip grows on the right as b advances, the rose
   strip is *added* on the left here because a moves left (ȧ < 0 makes −ȧF(a) positive); at t = 0 the three rates are
   0.6, 1.0 and 0.3, total 1.9 (checked in the code below)."
10. `nb.note` — **N48 [B]** "**Control volume** V*(t) with **control surface** A*(t), outward unit normal n and surface
    velocity **b** (Fig. 3.18). The surface need not follow the fluid: b = u for a *material* volume, b = 0 for a volume
    fixed in space, anything else for a piston, a balloon or a moving layer. `ch03.GrowingSphere`, `GrowingCylinder`,
    `MovingBox`, `GrowingCone` are concrete examples."
11. `nb.primer("signed swept volume of a moving surface", "In a short time Δt a surface patch dA moving with velocity b
    sweeps a thin prism of height $(\mathbf b\Delta t)\cdot\mathbf n$ — only the normal component counts (sliding along
    the surface sweeps nothing). With the outward n the volume is **signed**: positive where the wall advances outward
    (the region gains volume), negative where it retreats.", code="n = np.array([1.0, 0.0]); dA, dt = 0.01, 0.1
    \nfor b in ([2.0, 0.0], [0.0, 5.0], [-1.0, 3.0]):\n    print(b, np.dot(b, n)*dt*dA)            # 0.002, 0.0 (sliding),
    -0.001 (retreating)")` — *expect:* `0.002`, `0.0`, `-0.001`.
12. `nb.note` ×4 — **N49 [B]** "The derivation starts from the definition of a time derivative (Fig. 3.18: solid at t,
    dashed at t + Δt):" equation
    `\frac{d}{dt}\int_{V^*(t)}F(\mathbf x,t)dV=\lim_{\Delta t\to0}\frac{1}{\Delta t}\Big\{\int_{V^*(t+\Delta t)}F(\mathbf x,t+\Delta t)dV-\int_{V^*(t)}F(\mathbf x,t)dV\Big\}`,
    ref "3.31". · **N50 [B]** "Splitting the new volume, $\Delta V\equiv V^*(t+\Delta t)-V^*(t)$, and Taylor-expanding F in
    time gives four terms:" equation
    `\int_{V^*(t+\Delta t)}F(\mathbf x,t+\Delta t)dV\cong\int_{V^*(t)}F\,dV+\int_{V^*(t)}\Delta t\frac{\partial F}{\partial t}dV+\int_{\Delta V}F\,dV+\int_{\Delta V}\Delta t\frac{\partial F}{\partial t}dV`,
    ref "3.32". · **N51 [B]** "The first cancels, the last is second order:" equation
    `\frac{d}{dt}\int_{V^*(t)}F\,dV=\lim_{\Delta t\to0}\frac{1}{\Delta t}\Big\{\int_{V^*(t)}\Delta t\frac{\partial F}{\partial t}dV+\int_{\Delta V}F\,dV\Big\}`,
    ref "3.33". · **N52 [B]** "and the sliver is a surface integral:" equation
    `\int_{\Delta V}F(\mathbf x,t)dV\cong\int_{A^*(t)}F(\mathbf x,t)(\mathbf b\Delta t\cdot\mathbf n)dA\quad\text{as}\quad\Delta t\to0`,
    ref "3.34". "Each of these is a step of D22 below, where the moves between them are filled in."
13. `nb.derivation("D22", …)` — Part F D22 (12 steps), ref "3.35", **with `check_src`** (the ★★★ sympy cell of Part F D22).
14. `nb.note` — **N53 [B]** "**Two readings of (3.35).** (1) With F = 1 it says the volume changes at the rate the surface
    sweeps: $dV^*/dt=\int_{A^*}\mathbf b\cdot\mathbf n\,dA$; for a small material volume (b = u) this is (3.14),
    $\frac{1}{\delta V}\frac{D}{Dt}(\delta V)=\frac{\partial u_i}{\partial x_i}$ — D23. (2) For a small material volume it
    contains the material derivative (3.5), $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$, plus a
    term F∇·u from the volume change — D24. Where b·n > 0 the surface advances, where b·n < 0 it retreats; one signed
    term covers both. Only for a fixed volume (b = 0) may d/dt pass inside the integral."
15. `nb.derivation("D23", …)` — Part F D23 (5 steps), ref "3.14".
16. `nb.derivation("D24", …)` — Part F D24 (7 steps), ref "3.5".
17. `nb.worked_example("a balloon growing in a warming room", "Sphere of radius R = 1 m growing at Ṙ = 0.1 m/s; F = t
    (a 'heat content per volume' rising uniformly, [F] per s = 1) evaluated at t = 1 s. 1. Volume term:
    $\int_{V^*}\partial F/\partial t\,dV=1\times\tfrac43\pi R^3=4.189$. 2. Surface term: on the sphere b·n = Ṙ = 0.1 m/s,
    F = 1: $\oint F\,\mathbf b\cdot\mathbf n\,dA=1\times0.1\times4\pi R^2=1.257$. 3. Total 5.445. 4. Direct: $\int
    F\,dV=t\cdot\tfrac43\pi R(t)^3$, whose derivative is $\tfrac43\pi R^3+t\cdot4\pi R^2\dot R=4.189+1.257=5.445$ ✓.")`
18. `nb.code` — *code:* `cv = ch03.GrowingSphere(0.9, 0.1)` (R(t) = 0.9 + 0.1t, so R = 1.0 m at t = 1 s) ·
    `F = lambda x, t: t + 0*x[0]`; `dF = lambda x, t: 1.0 + 0*x[0]` · `print(ch03.reynolds_transport(F, dF, cv, 1.0))` ·
    `print(ch03.rtt_check(F, dF, cv, 1.0))` · Leibniz: `Fw, dFw = ch03.rtt_field("warming", dim=1)`;
    `print(ch03.leibniz_terms(Fw, dFw, 1.0, 3.0, -0.2, 0.4, 0.0))` · `print(ch03.leibniz_example(2.0))` · 2-D ellipse:
    `F2, dF2 = ch03.rtt_field("ramp", dim=2)`; `print(ch03.rtt_ellipse_2d(F2, dF2, 2.0, 1.0, 0.2, 0.1, (0.0, 0.0), (0.3,
    0.0), 0.0))`. *expect:* `(4.1888, 1.2566, 5.4454)` · `(5.4454, 5.4454)` (FD vs RTT to 1e-8) · `(0.6, 1.0, -0.3, 1.9)` ·
    `interior 18.667, upper 128.0, lower 8.0, total 138.667, exact 138.667` · `(0.0, 2.1991, 2.1991)`. *explain:* 1. a
    control volume with a moving surface; 2. the two terms of (3.35) by quadrature; 3. the same total from a finite
    difference of the volume integral (the left side of (3.31)); 4. Leibniz with E7's default numbers (the lower term
    −0.3 is ȧF(a); subtracting it *adds* 0.3); 5. the D21 check case F = x²t, a = t, b = t² at t = 2 s; 6. a translating,
    growing ellipse in a ramp field (E7's RTT mode): π(ȧb + aḃ)(1 + 0.5c_x) + πab·0.5ċ_x = 1.2566 + 0.9425.
19. `nb.check_agree` — **from scratch (curation §7):** midpoint sums for the growing sphere: volume term on a spherical
    grid (nr × nθ × nφ = 16 × 16 × 32, weights r² sin θ Δr Δθ Δφ) and surface term (nθ × nφ, weights R² sin θ Δθ Δφ)
    → `assert np.allclose([vol, surf], ch03.reynolds_transport(F, dF, cv, 1.0)[:2], rtol=1e-3)`; and the finite
    difference `(I(t+dt) − I(t−dt))/(2dt)` of the same midpoint volume sum → `assert np.allclose(fd, vol + surf,
    rtol=1e-3)`.
20. `nb.plotly` — **N51** `slider_figure` over Δt ∈ 10^[−4, −1] s (16 steps): for `GrowingSphere(1.0, 0.1)` with F = x₁² t (a
    field that varies in space), bars of the four terms of (3.32) (`ch03.swept_terms_sphere`) on a log scale, plus the
    ratio T4/T3 = O(Δt) printed in the title. Followed by a static log–log `nb.figure` of |T4| and |exact − T1 − T2 −
    T3| vs Δt (slopes 2 marked, `observed_order` printed). *see:* "four bars; the last one falls twice as fast as the
    others." *read:* "(3.33) drops T4 because after dividing by Δt it still vanishes: it is second order." *change:*
    "…the sphere grew ten times faster: T3 and T4 both grow ×10 (both ∝ ΔV), but T4/T3 is still ∝ Δt."
21. `nb.note` — **N54 [B]** "**Ex. 3.2 — a growing cone.** A right circular cone of fixed height h has base radius r(t)
    growing at ṙ. Directly, $V=\tfrac13\pi hr^2$ gives $dV/dt=\tfrac23\pi hr_o\dot r$. By (3.35) with F = 1 and the cone
    as V*, only the sloping side contributes: a side point at height z moves outward at $(z/h)\dot r\,\mathbf e_R$, the
    outward normal is $\mathbf n=\mathbf e_R\cos\theta-\mathbf e_z\sin\theta$ (θ = the cone's half-angle, tan θ = r_o/h —
    a fourth θ in this chapter), and the slanted area element is $dA=z\tan\theta\,d\varphi\,dz/\cos\theta$:" equation
    `\frac{dV}{dt}=\int_{z=0}^{h}\int_{\varphi=0}^{2\pi}\frac zh\dot r\,\mathbf e_R\cdot(\mathbf e_R\cos\theta-\mathbf e_z\sin\theta)\,z\tan\theta\,d\varphi\frac{dz}{\cos\theta}=\frac{2\pi\dot r\tan\theta}{h}\int_0^hz^2dz=\tfrac23\pi h^2\dot r\tan\theta=\tfrac23\pi hr_o\dot r`,
    no ref. "⚠️ Two print slips taught corrected: the base's points move radially in its own plane, so it is **b·n = 0**
    there (not b = 0), which is all the theorem needs; the '[?]' in the printed integrand is just the factor z."
    followed by `nb.code` — *code:* `print(ch03.example_3_2(1.0, 0.5, 0.1))`. *expect:* `direct 0.10472, rtt_dblquad
    0.10472, rtt_cv 0.10472` m³/s (h = 1 m, r_o = 0.5 m, ṙ = 0.1 m/s). *explain:* "the direct derivative, the book's
    surface integral by `dblquad`, and the general quadrature on `GrowingCone` agree."
22. `nb.figure` — **N55 [B]** our Fig. 3.18 in 2-D (`rtt_blob_figure(a=2, b=1, adot=0.2, bdot=-0.1, cdot=(0.3, 0), F_name=
    "ramp", dt=0.3)`): (left) the ellipse at t (solid) and t + Δt (dashed) over a heat map of F, b arrows on the
    boundary, the swept band coloured blue where b·n > 0 and rose where b·n < 0; (right) waterfall bars: volume term
    (blue) + surface term (orange) = total (purple) next to the finite-difference d/dt ∫F (grey marker). *see:* "a
    boundary that advances on the right and retreats at the top; blue and rose slivers; bars that land on the grey
    marker." *read:* "the rose slivers are *negative* volume — the same formula counts them with the right sign; the
    sum of the bars is the measured rate." *change:* "…b = 0 (a fixed box): the band and the orange bar vanish and d/dt
    passes inside the integral."
23. `nb.explainer("reynolds_transport_cv", heading="What changes inside a moving box?", why="A moving boundary sweeps a
    band whose colour shows the sign of b·n; bars for the volume and surface terms of (3.35),
    $\frac{d}{dt}\int_{V^*}F\,dV=\int_{V^*}\frac{\partial F}{\partial t}dV+\int_{A^*}F\,\mathbf b\cdot\mathbf n\,dA$, add up to the
    measured rate of change of ∫F, and a Δt slider shows the dropped second-order term vanish. Modes: Leibniz in 1-D, a
    deforming ellipse in 2-D, and the growing cone of Ex. 3.2.", tries=["Preset 'fixed volume (b = 0)': the orange bar
    disappears.", "Preset 'rigid translation in a uniform F': the two terms cancel exactly — a uniform F carried along
    does not change.", "Make one side retreat: its sliver turns rose and subtracts.", "Derivation tab, D22: at the 'drop
    second order' step, scrub Δt on the log–log view."])`
24. `nb.md` — **What would change if…** "…F were the density ρ and V* a material volume (b = u)? The left side is the
    rate of change of the lump's mass, which is zero; D24's steps then give $\frac{D\rho}{Dt}+\rho\nabla\cdot\mathbf u=0$ —
    the continuity equation that opens Chapter 4. With F = ρu you get Newton's law for the lump; with F = ρ(e + u²/2),
    the energy equation."

### A.7 End matter — S01, S02, summary
1. `nb.pointer` — **S01** "Practise on the book's Exercises 3.1–3.30 (curvilinear operators, flow lines, the Galilean
   proof, strain and rotation of simple flows, circulation, RTT checks). The derivations the text leaves to Exercises
   3.3, 3.12, 3.17–3.20, 3.23, 3.26, 3.28 and 3.30 are written out above (D03, D06, D10, D11, D13, D17, N37, D20, D23,
   D24) in our own words."
2. `nb.pointer` — **S02** "Further reading: Van Dyke's *An Album of Fluid Motion* and Samimy et al.'s *A Gallery of Fluid
   Motion* for photographs of streak and path lines; Riley, Hobson & Bence for Leibniz's rule; Thompson (1972) for the
   geometric RTT argument (all cited in the book's literature list)."
3. `nb.summary(clicked=[…15 bullets, one per CORE…], feeds_forward=[…], left_out=[…])` — clicked: C01 "A particle is a
   label; a field is a function of (x, t); (3.2), $F[\mathbf r(t;\mathbf r_o,t_o),t]=F(\mathbf x,t)$ at $\mathbf x=\mathbf r$,
   translates between them." · C02 "$\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$ (3.5): what a
   moving parcel feels = what a fixed probe sees + what it meets by moving." · C03 "A streamline is tangent to u at one
   frozen instant: $dx/u=dy/v=dz/w$ (3.7)." · C04 "A path line follows one particle, $d\mathbf r/dt=\mathbf u(\mathbf r,t)$
   (3.8); a streak line collects everything from one port; in unsteady flow the three lines differ." · C05 "The
   acceleration $\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u$ is the same for every observer
   moving at constant velocity (3.9); only its split is not." · C06 "$du_i=(\partial u_i/\partial x_j)dx_j$ (3.10): the
   neighbourhood of a point moves linearly." · C07 "The diagonal of S is the stretching rate per length; along any
   direction n it is n·S·n." · C08 "The off-diagonal S₁₂ is half the closing rate of a right angle; rigid motion has
   S = 0." · C09 "$\frac1{\delta V}\frac{D(\delta V)}{Dt}=\frac{\partial u_i}{\partial x_i}=S_{ii}$ (3.14): divergence is
   the swelling rate." · C10 "A fluid element spins at ½ω — the average of any perpendicular pair of threads — and the
   spin depends on the observer's rotation (ω′ = ω − 2Ω)." · C11 "$du_i=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i$
   (3.19): deformation plus rigid rotation." · C12 "S stretches a small sphere into an ellipsoid on its principal axes,
   $d\bar u_\alpha=\bar S_{\alpha\alpha}d\bar x_\alpha$ (3.21), at the first instant." · C13
   "$\omega_z=\frac1r\frac{\partial}{\partial r}(ru_\theta)-\frac1r\frac{\partial u_r}{\partial\theta}$ (3.23): solid-body
   rotation spins every element (2ω₀); the line vortex spins none (except its axis)." · C14 "Real vortices have a
   solid-body core and an irrotational outside; the peak wind is at σ (Rankine) or 1.1209σ (Gaussian)." · C15
   "$\frac{d}{dt}\int_{V^*}F\,dV=\int_{V^*}\frac{\partial F}{\partial t}dV+\int_{A^*}F\,\mathbf b\cdot\mathbf n\,dA$ (3.35):
   change in place + what the moving walls sweep (signed)." · feeds_forward: "Ch. 4: RTT + (3.5) → continuity, Cauchy's
   equation, energy; S → the Newtonian stress law; rotating frames (ω′ = ω − 2Ω) → Coriolis." · "Ch. 5: vorticity,
   circulation $\Gamma=\oint\mathbf u\cdot d\mathbf s$ (3.18), vortex stretching (S along ω), Lamb–Oseen vortex." · "Ch.
   6: irrotational flow (3.17) $\boldsymbol\omega=0$, point vortex $u_\theta=B/r$ (3.25), the cylinder of Fig. 3.2." · "Ch.
   7: particle orbits via (3.8), wave frames via (3.9)." · "Ch. 13: absolute vs relative vorticity, frontogenesis on
   principal strain axes, cyclone profiles." · left_out: "Curvilinear forms of ∇, ∇², (u·∇)u (Appendix B; Ch. 4)." ·
   "Solutions for the velocity potential (Ch. 6)." · "Viscous spreading of the Gaussian vortex (Ch. 5, 8)."
4. `nb.save()`.

### A.8 Placement table (ID → notebook section → host block → call)
| ID | Section | Host | Call |
|---|---|---|---|
| N01, N02, N03, N05, N06 | §3.1 | (none; tagged C01) | `nb.note` ×5 (+ code, figure, plotly) |
| N04 | §3.3 | C05 | `nb.note` + figure (pointer sentence in §3.1) |
| C01 · N07 · N08 · D01 | §3.2 | C01 | `nb.core`, `nb.note` ×2, `nb.derivation` |
| C02 · N09 · N10 · N11 · N12 · D02 | §3.2 | C02 | `nb.core`, `nb.note` ×5 (N09 as two), `nb.derivation`; **E2** |
| C03 · N14 · N15 · D03 | §3.3 | C03 | `nb.core`, `nb.note` ×2, `nb.derivation` |
| C04 · N13 · N16 · N17 · N18 · N19 · D04 · D05 | §3.3 | C04 | `nb.core`, `nb.note` ×4, figure (N19), `nb.derivation` ×2; **E1** |
| C05 · N04 · N20 · N21 · N22 · N23 · N24 · D06 | §3.3 | C05 | `nb.core`, `nb.note` ×6, `nb.derivation`; **E3** |
| R01 · R02 · R03 | §3.4 | before C06 | `nb.recap` ×3 |
| C06 · N25 · N26 · D07 | §3.4 | C06 | `nb.core`, `nb.note` ×2, `nb.derivation` |
| C07 · D08 | §3.4 | C07 | `nb.core`, `nb.derivation` |
| C08 · N27 · D09 · D10 | §3.4 | C08 | `nb.core`, `nb.note`, `nb.derivation` ×2 |
| C09 · D11 | §3.4 | C09 | `nb.core`, `nb.derivation`; **E4** |
| R04 · R05 · R06 · R07 | §3.4 | before C10 | `nb.recap` ×4 |
| C10 · N28 · N29 · D12 · D13 · D14 | §3.4 | C10 | `nb.core`, `nb.note` ×2, `nb.derivation` ×3 |
| C11 · N30 · N33 · D15 | §3.4 | C11 | `nb.core`, `nb.note` ×2, `nb.derivation` |
| C12 · N31 · N32 · D16 | §3.4 | C12 | `nb.core`, `nb.note` ×2, `nb.derivation`; **E5** |
| N34 · N35 | §3.5 head | (none; tagged C10, C12) | `nb.note` ×2 + Fig. 3.14 animation |
| C13 · N36 · N37 · N38 · N39 · N40 · N41 · D17 | §3.5 | C13 | `nb.core`, `nb.note` ×6, `nb.derivation` |
| C14 · N42 · N43 · N44 · D18 · D19 · D20 | §3.5 | C14 | `nb.core`, `nb.note` ×3, `nb.derivation` ×3; **E6** |
| C15 · N45–N55 · D21 · D22 · D23 · D24 | §3.6 | C15 | `nb.core`, `nb.note` ×11, `nb.derivation` ×4 (D22 with `check_src`); **E7** |
| S01 · S02 | end | — | `nb.pointer` ×2 |

---

## Part B — explainer storyboards

Common to all eight: created with `tools/new_viz.py`; `<meta name="viz:chapter" content="ch03">`; tabs Walkthrough ·
Explore · Explain · Derivation · Equations · Code · Check; `autoplay: false` and `play: false` on every walkthrough step
that quotes numbers (ch02 lesson); every displayed number is computed by a JS function that mirrors a `ch03` callable and
is proved by `selftest()` parity rows (`py:` expressions use only `ch03.…`, `np.pi`, lists, floats and strings — no
builtins, no lambdas; index tuples/dicts down to one number); Explain is "Explanation & interpretation" in numbered
sections built with `Viz.work.step / line / box / table / hint / interpret`; derivation steps are copied from Part F
(same `did` titles, same step count; phones shorten *why* to its first sentence; plain-text *why* never contains raw
TeX); every tour, Explain, Derivation, equation and quiz text that names a book equation **writes it out** next to its
number (convention 11; `tools/eq_refs.py` checks). Colours as in convention 10. Angles in degrees on controls, radians
inside every function. Linear flows use the closed-form `expm2` of ch02 E3 (promotion candidate `Viz.num.expm2`); path
lines of nonlinear or time-dependent fields use `Viz.num.rk4Step` with the same field function the parity rows call.
Each fits 360×640 … 1920×1080 and the 1000×700 notebook frame with no scrolling (fit plan per explainer); a view hidden
on portrait phones never carries a step's key number (its number is repeated in a visible title or readout).

### E1 · flow_lines_unsteady
- **Title:** "Three lines through one point — why do they disagree?" · **Summary:** "A clock runs on an unsteady flow:
  the streamline pattern is redrawn every instant, one particle leaves its trail and a port leaks dye — three curves
  that coincide only when nothing changes in time." · **CORE:** C03, C04 (also N13, N14, N16, N17, N18, N19) ·
  **Reference:** `angular_frequency_explorer_1.html` (linked views on one clock, ghost reference, modes, end-of-run
  summary).
- **meta:** `viz:order 1` · `viz:sections 3.3` · `viz:equations 3.7 3.8` · `viz:fluidpy ch03.example_3_1
  ch03.unsteady_flow_preset ch03.preset_field ch03.streamline ch03.pathline ch03.streakline ch03.streamline_slope` ·
  `viz:derivations D03 D04 D05`.
- **Physics (JS ↔ Python):** `field(mode, x, y, t, p) -> [u, v]` ↔ `ch03.unsteady_flow_preset(name, x, y, t, **p)` (the
  four modes of Part C 5.7) · `pathClosed(t, tRelease, p)` (uniform fields: x = U₀(t − t_o) + ξ_o(sin ωt − sin ωt_o),
  y = ξ_o(cos ωt_o − cos ωt)) and `pathRK4(x0, t0, t1, n)` (rotating strain) ↔ `ch03.pathline(ch03.preset_field(…), r0,
  t0, t_eval)` · `streak(t, releases)` (one particle per release time, closed form or carried by RK4 each frame) ↔
  `ch03.streakline(ch03.preset_field(…), x0, t, t_release)` · `streamlineThrough(x0, t, sMax)` (RK4 in arc length at
  frozen t) ↔ `ch03.streamline` · `centres(tp, p)` → path centre (−ξ_o sin ωt′, ξ_o cos ωt′), streak centre (ξ_o sin ωt′,
  −ξ_o cos ωt′) ↔ `ch03.example_3_1(t_prime, xi0, omega)["path_center"|"streak_center"]`.
- **Modes** (`mode` chips): **Ex. 3.1** (u = ωξ_o cos ωt, v = ωξ_o sin ωt) · **+ current** (u gains U₀) · **steady**
  (uniform, direction fixed at β = 30°) · **rotating strain** (u = s[x cos 2Ωt + y sin 2Ωt, x sin 2Ωt − y cos 2Ωt], port at
  (0.6, 0.4) m; Ω = 0 is steady).
- **Views** (one clock `t`; rows [1.25, 0.8]):
  1. `flow` "The flow plane" (`equal: true`; x, y ∈ [−2.5, 2.5] ξ_o, in m) — faint teal streamline family (through a 7 × 7
     seed grid at the current instant), the **bold teal streamline through the port**, the port (black ring) with its
     velocity arrow (black), **dye particles** (rose dots, one per release every T/48, faded by age; "dye since forever"
     pre-fills the last period so the Ex. 3.1 streak circle is complete at t = 0), **one tagged particle** released at
     t = 0 with its growing **orange trail**, a clock "ωt = …°". Pointer: click a dye dot → inspector; in rotating-strain
     mode drag the port.
  2. `lines` "The three lines at t′ = t" (`equal: true`) — the answer sheet: the closed-form streamline y = x tan ωt′
     (teal), path circle (orange) and streak circle (rose) with centres ×, radii ξ_o, the angle ωt′ as an arc at the
     origin, the common tangent (grey dashed); for non-Ex.-3.1 modes the three numerically computed curves (streamline
     through the port, the tagged particle's full path over [t − T, t + T], the current streak) with the label "no closed
     form". Title carries the key number: "centres 2ξ_o = 2.00 m apart · slope tan ωt′ = 0.577".
  3. `signal` "Velocity at the port" (`hidePortrait: true`) — u(t) (solid) and v(t) (dashed) in purple over one period
     [t − T, t], faint whole period ahead as a ghost, the current-time marker; in steady modes two flat lines.
  Portrait: `flow` (top) + `lines` (bottom); `signal` hidden (its numbers are in the readouts u, v).
- **Controls:** `omega` "Angular frequency $\omega$" 0.2…3, step 0.05, default 1 s⁻¹, help "how fast the flow direction
  turns; period 2π/ω" (in rotating-strain mode it is the axes' turning rate Ω) · `xi0` "Amplitude $\xi_o$" 0.2…2, step
  0.05, default 1 m, help "radius of every particle's circle" · `U0` "Mean current $U_0$" 0…2, step 0.05, default 0 m/s
  (*optional*; active in "+ current") · `s` "Strain rate $s$" 0.2…2, default 1 s⁻¹ (*optional*; rotating strain) · `hist`
  toggle "dye since forever" default on (*optional*). Transport `t`: 0 → 12.566 s (two periods at ω = 1), rate 1 s/s,
  `end: 'hold'`, hold 2 s, with an end-of-run `Viz.card`: Ex. 3.1 "After two periods: the particle is back at the port;
  the dye circle and the streamline are where they started — only the clock moved." / + current "After two periods the
  particle drifted U₀·2T = … m downstream; the dye is a wave with wavelength U₀T = … m." / steady "One line all along:
  nothing changed in time."
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **modes** (4 fields) ·
  **presets** "ωt′ = 0" {mode: 'ex31', t: 0} · "ωt′ = 30° (Fig. 3.7)" {t: 0.5236} · "ωt′ = 90°" {t: 1.5708} · "steady"
  {mode: 'steady'} · "steady strain (Ω = 0)" {mode: 'rotstrain', omega: 0} · "mean current" {mode: 'current', U0: 0.5} ·
  **status** (steady modes: "🟰 steady: the three lines coincide"; unsteady: "↔ unsteady: three different curves ·
  centres 2ξ_o = 2.00 m apart"; + current with U₀ < ωξ_o: "➰ the path loops: the current is slower than the sloshing";
  U₀ ≥ ωξ_o: "〰 no loops: the current wins") · **inspector** (click a dye dot: "released at t_o = −2.094 s (ωt_o =
  −120°); now $x=\xi_o(\sin\omega t-\sin\omega t_o)=1\times(0.500+0.866)=1.366$ m, $y=\xi_o(\cos\omega t_o-\cos\omega t)=1\times(-0.500-0.866)=-1.366$ m"
  — highlights that particle's path in faint orange).
- **Readouts:** "u, v at the port" (m/s) · "slope tan ωt′" · "path centre" · "streak centre" · "tagged particle".
- **Explain** ("Explanation & interpretation"):
  0. *What the windows show* — "**The flow plane** is the experiment: the <b class="c-teal">teal</b> lines are
     streamlines of the velocity field *at this instant* (the bold one passes through the port); the
     <b class="c-rose">rose</b> dots are dye that left the port at earlier times; the <b class="c-orange">orange</b>
     trail is the path of one particle released at t = 0. **The three lines at t′** redraws them cleanly from the closed
     forms of Ex. 3.1. **Velocity at the port** (hidden on phones) shows how u and v swing round once per period."
  1. *The field at the port now* — "$u=\omega\xi_o\cos\omega t=1\times1\times\cos30°=$ **0.866** m/s,
     $v=\omega\xi_o\sin\omega t=$ **0.500** m/s; speed ωξ_o = **1.000** m/s, direction ωt = **30°**. The field is the same at
     every point — only its direction changes in time."
  2. *The streamline through the port (3.7)* — "From $dx/u=dy/v=dz/w$ (3.7): $\frac{dy}{dx}=\frac vu=\tan\omega t'=\tan30°=$
     **0.577**, so $y=x\tan(\omega t')$ — a straight line at 30°. Frozen clock: it will be a different line a moment
     later." (boxed)
  3. *The path line (3.8)* — "Solving $d\mathbf r/dt=\mathbf u(\mathbf r,t)$ (3.8) with $\mathbf r(t')=0$ gives the circle
     $(x+\xi_o\sin\omega t')^2+(y-\xi_o\cos\omega t')^2=\xi_o^2$: centre (−0.500, **0.866**) m, radius **1.000** m, run
     counterclockwise once per period $2\pi/\omega=$ **6.283** s." (boxed)
  4. *The streak line* — "Every dot left the port at its own time t_o; at t′ they lie on
     $(x-\xi_o\sin\omega t')^2+(y+\xi_o\cos\omega t')^2=\xi_o^2$: centre (0.500, **−0.866**) m, radius 1.000 m — the path
     circle reflected through the port." (boxed)
  5. *Where they meet* — "The two centres are $2\xi_o=$ **2.000** m apart on opposite sides of the port. At the port all
     three curves have slope tan ωt′ = 0.577: they touch there, because each is tangent to the velocity at (port, t′)."
  6. *With a mean current* (or hint "choose '+ current' to see the path drift") — "The particle now drifts
     $U_0\cdot2\pi/\omega=$ **3.142** m per period; its path loops backwards while $U_0<\omega\xi_o$ (0.5 < 1.0 now); the dye
     makes a wave of wavelength $U_0T$. The streamline is still straight (the field is still uniform) but tilted:
     slope $v/u=\omega\xi_o\sin\omega t/(U_0+\omega\xi_o\cos\omega t)$."
  7. *At the current time* — "t = `live t` s, ωt = `live deg`°, `live ndye` dye particles, the tagged particle at
     (`live px`, `live py`) m."
  8. *Reading the current setting* — Ex. 3.1: "Unsteady: the streamline turns with the clock, the particle runs round
     its own circle, the dye sits on the mirror circle — a photograph of dye (streak line) or a long exposure of one
     speck (path line) would both mislead you about the instantaneous direction." · steady / Ω = 0: "Steady: nothing at a
     fixed point changes, so every particle follows the streamline it starts on and the dye lines up behind it — the
     three curves are one (N13)." · + current, U₀ < ωξ_o: "Unsteady with drift: loops, because part of each period the
     sloshing beats the current." · U₀ ≥ ωξ_o: "The current always wins: a wavy path without loops." · rotating strain:
     "The streamlines are hyperbolas whose axes turn at Ω; path and streak lines peel away from them — set Ω = 0 and they
     collapse."
- **Derivation tab:**
  - **D03** (5 steps) `view: 'lines'`; goal `set {mode: 'ex31', t: 0.5236}`; step 1 highlights the velocity arrow at
    the port and a short ds along the teal line; step 3 **live** "$dx/u=dy/v$: $\frac{dy}{dx}=\frac{0.500}{0.866}=0.577=\tan30°$";
    step 5 `watch: "the bold teal line in The flow plane is this ODE integrated at frozen t"`; interpret `s => "At
    ωt′ = … the streamline through the port is the line at … °."`.
  - **D04** (7 steps) `view: 'flow'`; goal `set {mode: 'ex31', t: 0}`; step 2 `set {t: 1.5}` (the orange trail half
    drawn); step 3 **live** "$c_x=-\xi_o\sin\omega t'=$ …, $c_y=\xi_o\cos\omega t'=$ …"; step 7 `set {t: 6.283}` `watch:
    "the orange trail closes on the dashed circle"`; interpret: "centre (…, …), radius … m — one lap per 2π/ω = … s".
  - **D05** (8 steps) `view: 'flow'`; goal `set {mode: 'ex31', t: 0.5236, hist: true}`; step 2 colours the dye dots by
    release time (older = paler); step 6 **live** "releases over one period T = … s fill the circle"; steps 7–8 **live**
    slopes "path: −(0 + 0.500)/(0 − 0.866) = 0.577; streak: 0.577"; interpret: "all three touch the port with slope
    tan ωt′ = …".
- **Code:**
  ```python
  om, xi0, tp = {{omega}}, {{xi0}}, {{t}}              # ω [1/s], ξ_o [m], drawing instant t' [s]
  d = ch03.example_3_1(tp, xi0=xi0, omega=om)           # closed forms of Ex. 3.1 (D04, D05)
  print(d["slope"])                                     # tan ωt' = {{slope}}  (streamline slope, (3.7))
  print(d["path_center"], d["streak_center"])           # ({{pcx}}, {{pcy}}) and ({{scx}}, {{scy}}) m
  u = ch03.preset_field("ex31", omega=om, xi0=xi0)      # the velocity field as a callable u(x, t)
  sl = ch03.streamline(u, [0.0, 0.0], tp, s_max=2.0)    # tangent to u at frozen t'
  p = ch03.pathline(u, [0.0, 0.0], tp, [tp, tp + 1.0])  # (3.8): dr/dt = u(r, t), r(t') = 0
  tr = np.linspace(tp - 2*np.pi/om, tp, 48)             # dye released over the last period
  s = ch03.streakline(u, [0.0, 0.0], tp, tr)            # positions at t' → circle about the streak centre
  ```
- **Walkthrough (7 steps):** 1. "One point, three questions" — "A port on the sea floor leaks dye while the water
  sloshes. Which way is the water going *now*? Where does one speck go? Where is the dye? Three different questions."
  `set {mode: 'ex31', t: 0}`, highlight `view:flow` · 2. "Freeze the clock" — "Streamlines follow the arrows at one
  instant. Here the field is the same everywhere, so they are straight: $dx/u=dy/v$ (3.7) gives slope tan ωt′." `set
  {t: 0.5236}`, `eq: 'stream'`, `derive: {id: 'D03', step: 3}` · 3. "Follow one speck" — "Press ▶. The orange particle
  obeys $d\mathbf r/dt=\mathbf u(\mathbf r,t)$ (3.8); because the arrows keep turning it runs round a circle of radius
  ξ_o." `play: true`, `eq: 'path'`, `readouts: ['pc']` · 4. "Photograph the dye" — "Each rose dot left the port at a
  different time. Right now they sit on the mirror circle, centre (ξ_o sin ωt′, −ξ_o cos ωt′)." `set {t: 0.5236}`, `play:
  false`, `inspect: true`, `derive: {id: 'D05', step: 5}` · 5. "They touch at the port" — "All three curves have slope
  tan ωt′ at the port: each is tangent to the velocity there. Compare with the answer sheet." highlight `view:lines`,
  `readouts: ['slope']` · 6. "Make it steady" — "Choose steady: the arrows stop turning, and the three curves collapse
  onto one line." `set {mode: 'steady'}`, `play: true` · 7. "Your turn" — "Predict first: with a mean current U₀ =
  0.5 m/s, will the particle come back to the port? Then press ▶." `set {mode: 'current', U0: 0.5, t: 0}`, `controls:
  ['U0', 'omega']`.
- **Equations:** `stream` "Streamlines" ref 'Eq. (3.7)' `dx/u=dy/v=dz/w` live "slope $v/u=$ 0.500/0.866 = 0.577" symbols
  [['u,v,w', 'velocity components', 'm/s'], ['dx,dy,dz', 'arc element', 'm']] · `path` "Path lines" ref 'Eq. (3.8)'
  `d\mathbf r/dt=[\mathbf u(\mathbf x,t)]_{\mathbf x=\mathbf r}=\mathbf u(\mathbf r,t),\ \mathbf r(t_o)=\mathbf r_o` live "r(t′) = 0" ·
  `field` "Ex. 3.1 field" ref 'Ex. 3.1' `u=\omega\xi_o\cos(\omega t),\quad v=\omega\xi_o\sin(\omega t)` live "(0.866,
  0.500) m/s" · `pcirc` "Path line of Ex. 3.1" ref 'Ex. 3.1' `\big(x+\xi_o\sin\omega t'\big)^2+\big(y-\xi_o\cos\omega t'\big)^2=\xi_o^2`
  live "centre (−0.500, 0.866)" · `scirc` "Streak line of Ex. 3.1" ref 'Ex. 3.1'
  `\big(x-\xi_o\sin\omega t'\big)^2+\big(y+\xi_o\cos\omega t'\big)^2=\xi_o^2` live "centre (0.500, −0.866)".
- **Check yourself:** (1) "At ωt′ = 90°, where are the centres of the path and streak circles (ξ_o = 1 m)?" — "Path
  centre (−1, 0), streak centre (1, 0); the streamline is the y-axis: the picture of ωt′ = 0 turned by 90°." `set
  {mode: 'ex31', t: 1.5708}` · (2) "Double ω at fixed ξ_o. What happens to the circles and to the particle's lap time?"
  — "The circles keep radius ξ_o (the speed ωξ_o doubles but so does the turning rate); the lap takes 2π/ω, half as
  long." `set {omega: 2}` · (3) "In the steady mode, where is the dye relative to the particle's path?" — "On it: in a
  steady flow every particle released at the port follows the same streamline, so streak = path = streamline." `set
  {mode: 'steady'}` · (4) "With U₀ = 1.5 m/s and ωξ_o = 1 m/s, does the path loop?" — "No: the current always exceeds
  the sloshing speed, so u never reverses; loops need U₀ < ωξ_o." `set {mode: 'current', U0: 1.5}`.
- **Selftest parity rows:** `{name: 'path centre y (t′=0)', js: centres(0, {xi0: 1, omega: 1}).path[1], py:
  'ch03.example_3_1(0.0, 1.0, 1.0)["path_center"][1]', rtol: 1e-12}` · `{name: 'streak centre x (t′=π/6)', js:
  centres(0.5235987755982988, {xi0: 1, omega: 1}).streak[0], py: 'ch03.example_3_1(0.5235987755982988, 1.0,
  1.0)["streak_center"][0]', rtol: 1e-12}` · `{name: 'field v', js: field('ex31', 0.3, -0.2, 0.7, {omega: 1, xi0: 1})[1],
  py: 'ch03.unsteady_flow_preset("ex31", 0.3, -0.2, 0.7, omega=1.0, xi0=1.0)[1]', rtol: 1e-12}` · `{name: 'streak
  particle x', js: streakPoint(1.0, 0.0, {mode: 'ex31', omega: 1, xi0: 1})[0], py: 'ch03.streakline(ch03.preset_field("ex31",
  omega=1.0, xi0=1.0), [0.0, 0.0], 1.0, [0.0, 0.5])[0][0]', rtol: 1e-7}` · `{name: 'current path x (t=2)', js:
  pathClosed(2.0, 0.0, {U0: 0.5, omega: 1, xi0: 1})[0], py: 'ch03.pathline(ch03.preset_field("ex31_current", U0=0.5,
  omega=1.0, xi0=1.0), [0.0, 0.0], 0.0, [0.0, 2.0])[0][1]', rtol: 1e-7}` · `{name: 'rotating strain RK4 y', js:
  pathRK4([0.6, 0.4], 0, 1.0, 400, {s: 1, omega: 0.5})[1], py: 'ch03.pathline(ch03.preset_field("rotating_strain",
  s=1.0, Omega=0.5), [0.6, 0.4], 0.0, [0.0, 1.0])[1][1]', rtol: 1e-6}` · invariants `{name: 'tangency', js:
  slopePath(0.5235987755982988) - Math.tan(0.5235987755982988), expect: 0, atol: 1e-12}`.
- **Fit plan:** portrait 360×640: `flow` (top, square) + `lines` (bottom); `signal` hidden, u and v in the readouts;
  Explore: mode, omega, xi0 (+ U0/s at density ≤ 2); transport step/speed hidden by the engine. Landscape phone: one
  row, `flow` + `lines` (scoped CSS drops the `signal` row, Derivation tab exempt). Notebook 1000×700 and desktop:
  rows [1.25, 0.8], `signal` across the bottom.

### E2 · material_derivative_probe
- **Title:** "Why does the station warm while the air does not?" · **Summary:** "A fixed probe and a float carried by
  the wind read one moving temperature pattern on one clock; the term bars split the float's rate DT/Dt into the local
  rate ∂T/∂t and the advective rate u·∇T." · **CORE:** C01, C02 (also N07–N12) · **Reference:**
  `forced_damped_vibrations.html` (the explanation panel with boxed numbers and a regime-dependent reading) with
  `fid_formula_lab.html`'s term bars.
- **meta:** `viz:order 2` · `viz:sections 3.2` · `viz:equations 3.1 3.2 3.3 3.4 3.5 3.6` · `viz:fluidpy
  ch03.thermal_front ch03.thermal_front_terms ch03.material_derivative_terms ch03.lagrangian_map_example
  ch03.lagrangian_velocity_acceleration` · `viz:derivations D01 D02`.
- **Physics (JS ↔ Python):** `T(x, y, t, p)` ↔ `ch03.thermal_front(x, y, t, G, H, c, w)` (cold north; linear when
  `front = 'linear'`, tanh front of width w otherwise) · `terms(x, y, t, p) -> {T, dTdx, dTdy, local, advective, total}` ↔
  `ch03.thermal_front_terms(x, y, t, u, v, G, H, c, w)` · `stencilTerms(x, y, t, p, h)` (the generic central
  differences, shown in the inspector) ↔ `ch03.material_derivative_terms` · `float(t)` = start + (u, v)t (uniform wind) ·
  `floatRate(t)` = [T(float(t + δ)) − T(float(t − δ))]/2δ (the float's own thermometer) · stretching mode:
  `xmap(X, t, α) = X e^{αt}` ↔ `ch03.lagrangian_map_example`; `uLag(X, t, α)`, `aLag(X, t, α)` ↔
  `ch03.lagrangian_velocity_acceleration`.
- **Modes** (`mode` chips): **thermal front** (atmosphere; F = T) · **stretching map** (D01's flow x = X e^{αt}; F = the
  velocity u itself, so DF/Dt is the particle's acceleration).
- **Views** (rows [1.2, 0.8]):
  1. `map` "Temperature map" (`equal: true`; x, y ∈ [−500, 500] km) — heat map of T (blue cold north, orange warm south),
     isotherms every 1 K (thin), the **probe** (blue square at the origin) with its reading, the **float** (purple circle)
     with its trail, the wind arrow at the float (black), ∇T arrow at the probe (grey, pointing to warmer air). Pointer:
     drag the probe; drag the float's start point (sets `fx0, fy0`). Stretching mode: a 1-D strip x ∈ [0, 8] m, six
     particles (purple dots, labelled by X), the probe (blue tick at x_p), the velocity field u = αx as a coloured band.
  2. `bars` "The budget at the float" — term bars: local ∂T/∂t (blue) + advective u·∇T (amber) = DT/Dt (purple), in K/h,
     with the **float's measured rate** as a black diamond on the purple bar and the probe's ∂T/∂t as a blue tick; title
     "DT/Dt = −0.00 K/h · station +0.36 K/h" (the key numbers, so the phone never needs view 3).
  3. `series` "What each thermometer reads" (`hidePortrait: true`) — T(t) of the probe (blue) and of the float (purple)
     over the run, bold so far + faint whole run, slope triangles labelled ∂T/∂t and DT/Dt at the current time.
  Also `terms` config (Explore panel bars, click a term to isolate it on the map as a colour overlay) with items local
  (blue), advective (amber), total (purple).
- **Controls:** `V` "Wind speed $|\mathbf u|$" 0…20, step 0.5, default 10 m/s · `dir` "Wind blows toward" 0…360°, step 5,
  default 90° (north: a southerly wind), help "angle counterclockwise from east" · `G` "Gradient $|\nabla T|$" 0…3, step
  0.1, default 1 K/100 km · `H` "Heating of the air" −1…1, step 0.05, default 0 K/h (*optional*) · `carry` toggle "front
  carried by the wind (c = v)" default on (*optional*); stretching mode: `alpha` "Stretching rate $\alpha$" 0.1…1, default
  0.5 s⁻¹, `X` "Label $X$" 0.5…3, default 2 m. Transport `t`: 0 → 24 h (front) / 0 → 3 s (stretching), rate 2 h/s,
  `end: 'hold'`, card "After 24 h the float moved … km and its thermometer changed by … K; the station changed by … K —
  the difference is ∫u·∇T dt."
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **terms** (bars) · **transport** · **modes** ·
  **presets** "pure advection" {V: 10, dir: 90, H: 0, carry: true} · "pure heating" {V: 0, H: 0.5} · "wind along the
  front" {dir: 0} · "stalled front" {carry: false, V: 10, dir: 90} · "cold advection" {dir: 270} · "stretching map"
  {mode: 'stretch'} · **status** "🌡️ warm advection: −u·∇T = +0.36 K/h" / "❄️ cold advection: −u·∇T = −0.36 K/h" /
  "⏸ no advection: u ⟂ ∇T" (|u·∇T| < 1e-3 K/h) / stretching: "🚀 steady field, accelerating particles: a = α²x" ·
  **inspector** (click the map: at that point "∂T/∂y = −G sech²η = …; u·∇T = u ∂T/∂x + v ∂T/∂y = 0 + 10 × (−1.0e-5) =
  **−1.0e-4 K/s**; ∂T/∂t = H + Gc sech²η = …; DT/Dt = …", with the stencil version beside it).
- **Readouts:** "∂T/∂t at probe" (K/h) · "u·∇T at float" (K/h) · "DT/Dt at float" (K/h) · "float measured" (K/h).
- **Explain:**
  0. *What the windows show* — "The **map** is a temperature pattern (cold north); the <b class="c-blue">blue</b> square is
     a thermometer bolted to the ground, the <b class="c-accent">purple</b> circle a thermometer riding with the air. The
     **bars** split the purple rate into a <b class="c-blue">local</b> part and an <b class="c-amber">advective</b> part;
     the black diamond is what the float's thermometer actually measured over the last instant. The **time series**
     (hidden on phones) shows the two thermometer readings."
  1. *The temperature gradient* — "Largest gradient G = **1.0** K/100 km = 1.0 × 10⁻⁵ K/m, pointing south (toward warm
     air). At the float, η = (y − ct)/w = **0.00**, so $\partial T/\partial y=-G\,\mathrm{sech}^2\eta=$ **−1.0 × 10⁻⁵** K/m."
  2. *The advective term, component by component* — "$\mathbf u\cdot\nabla T=u\,\partial T/\partial x+v\,\partial T/\partial y=0\times0+10\times(-1.0\times10^{-5})=$
     **−1.0 × 10⁻⁴** K/s = **−0.36** K/h." (boxed)
  3. *The local term* — "The pattern slides north at c = **10** m/s and the air is heated at H = **0**:
     $\partial T/\partial t=H+Gc\,\mathrm{sech}^2\eta=$ **+0.36** K/h at the float's position." (boxed)
  4. *The material derivative (3.5)* — "$\frac{DT}{Dt}=\frac{\partial T}{\partial t}+\mathbf u\cdot\nabla T=0.36-0.36=$
     **0.00** K/h: the air parcel keeps its temperature." (boxed)
  5. *What the float measured* — "Its own thermometer changed by **0.000** K over the last 0.1 h — the black diamond sits
     on the purple bar, as (3.4) $\frac{d}{dt}F[\mathbf r,t]=\frac{DF}{Dt}$ promises."
  6. *What the station sees* — "At the probe: ∂T/∂t = **+0.36** K/h — warm advection, although no parcel warmed."
  7. *Streamwise form (3.6)* — "$\frac{DT}{Dt}=\frac{\partial T}{\partial t}+|\mathbf u|\frac{\partial T}{\partial s}$: along the
     wind, ∂T/∂s = **−1.0 × 10⁻⁵** K/m, times |u| = 10 m/s gives the same −0.36 K/h. (The book's printed (3.6) drops the
     F in the last term.)"
  8. *Stretching mode* (or hint "switch to the stretching map to see D01") — "x = X e^{αt} = 2 × e^{0.5} = **3.297** m;
     u = αx = **1.649** m/s at the probe for ever (local 0); the particle there accelerates at u ∂u/∂x = αx·α = **0.824**
     m/s² (advective) — a steady field with accelerating particles."
  9. *At the current time* — "t = `live t` h; float at (`live fx`, `live fy`) km; T_probe = `live Tp` K, T_float =
     `live Tf` K."
  10. *Reading the current setting* — warm advection (−u·∇T > 0): "The wind blows from warm toward cold: the station warms
      although the air does not — forecasters' 'warm advection'. The two thermometers disagree by exactly the amber
      bar." · cold advection: "Cold air is carried over the station: it cools while the air keeps its temperature." · no
      advection: "The wind runs along the isotherms (u ⟂ ∇T): both thermometers agree; only heating changes T." · H ≠ 0
      and carried: "Now the air itself warms at H; the station sees H plus the advective part." · stalled front: "The
      pattern is fixed to the ground: the station sees nothing (∂T/∂t = H), while the float cools as it moves into colder
      air — the same numbers, the other way round." · stretching: "Steadiness is about points, acceleration about
      particles."
- **Derivation tab:**
  - **D01** (6 steps) `view: 'map'`; goal `set {mode: 'stretch', alpha: 0.5, X: 2, t: 1}`; step 1 highlights the dot with
    label X = 2 and its slope; step 3 **live** "X = x e^{−αt} = 3.297 × e^{−0.5} = 2.000 m"; step 5 **live** "u(x, t) = αx =
    0.5 × 3.297 = 1.649 m/s"; step 6 `watch: "the blue probe at x = 3.297 m reads 1.649 m/s at every t"`; interpret
    `s => "With α = …, every particle passing the probe has u = … and a = …"`.
  - **D02** (9 steps) `view: 'bars'`; goal `set {mode: 'front', V: 10, dir: 90, G: 1, H: 0, carry: true, t: 0}`; step 3
    highlights the four inputs of the chain rule (the float's two velocity components on the map, ∂T/∂x, ∂T/∂y); step 4
    **live** "dr₁/dt = u = 0.0 m/s, dr₂/dt = v = 10.0 m/s"; step 6 `set {fx0: 150, fy0: -100}` `watch: "moving the float
    to a new start: the identity holds at every point"`; step 8 **live** "u·∇T = −0.36 K/h"; interpret: "at your
    setting DT/Dt = … K/h, ∂T/∂t = … K/h, u·∇T = … K/h".
- **Code:**
  ```python
  V, d = {{V}}, np.deg2rad({{dir}})                     # wind speed [m/s] and direction it blows toward
  u, v = V*np.cos(d), V*np.sin(d)                       # components ({{u}}, {{v}}) m/s
  G = {{G}}/1e5                                         # largest gradient [K/m] (G K per 100 km)
  c = v if {{carry}} else 0.0                           # the front moves with the wind's northward part
  T = ch03.thermal_front_terms({{fx}}, {{fy}}, {{ts}}, u, v, G, heating_K_per_s={{Hs}}, front_speed=c, width_m=2e5)
  print(T["local"]*3600)                                # ∂T/∂t = {{loc}} K/h
  print(T["advective"]*3600)                            # u·∇T  = {{adv}} K/h
  print(T["total"]*3600)                                # DT/Dt = {{tot}} K/h  (3.5)
  ```
- **Walkthrough (7 steps):** 1. "Two thermometers" — "A station on the ground and a thermometer riding with the air
  watch the same warm air arriving from the south. Do they agree?" `set {mode: 'front', t: 0}`, highlight `view:map` ·
  2. "The station warms" — "Press ▶: the blue reading climbs 0.36 K/h. That is ∂T/∂t — what happens *here*." `play:
  true`, `readouts: ['loc']` · 3. "The air does not" — "The purple float keeps its temperature: DT/Dt = 0. The pattern
  and the air move together." `set {t: 6}`, `play: false`, `readouts: ['tot']` · 4. "The missing piece" — "(3.5),
  $\frac{DT}{Dt}=\frac{\partial T}{\partial t}+\mathbf u\cdot\nabla T$: the amber bar u·∇T = −0.36 K/h makes up the
  difference." `terms: true`, `eq: 'md'`, `derive: {id: 'D02', step: 8}` · 5. "Turn the wind" — "Drag the direction to 0°:
  the wind runs along the isotherms and the amber bar vanishes." `set {dir: 0}`, `controls: ['dir']` · 6. "A steady field
  that accelerates" — "Stretching map: the probe always reads u = αx, yet each particle speeds up at α²x." `set {mode:
  'stretch', t: 1}`, `derive: {id: 'D01', step: 5}` · 7. "Your turn" — "Stall the front (untick 'carried'). Predict: which
  thermometer changes now, and at what rate? Then press ▶." `set {mode: 'front', carry: false, dir: 90, t: 0}`,
  `controls: ['carry', 'V']`.
- **Equations:** `bridge` "The bridge" ref 'Eq. (3.2)' `F[\mathbf r(t;\mathbf r_o,t_o),t]=F(\mathbf x,t)\ \text{when}\ \mathbf x=\mathbf r(t;\mathbf r_o,t_o)`
  live "T_float = T(x_float) = 288.15 K" · `chain` "Chain rule along a path" ref 'Eq. (3.3)'
  `\frac{d}{dt}F[\mathbf r,t]=\frac{\partial F}{\partial r_i}\frac{dr_i}{dt}+\frac{\partial F}{\partial t}` live "= … K/h" · `md`
  "Material derivative" ref 'Eq. (3.5)' `\frac{DF}{Dt}\equiv\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F=\frac{\partial F}{\partial t}+u_i\frac{\partial F}{\partial x_i}`
  live "0.36 − 0.36 = 0.00 K/h" · `sw` "Streamwise form" ref 'Eq. (3.6)' `\frac{DF}{Dt}=\frac{\partial F}{\partial t}+|\mathbf u|\frac{\partial F}{\partial s}`
  live "10 × (−1.0e-5) × 3600 = −0.36 K/h" note "the book prints |u|∂/∂s without F" · `lag` "Particle velocity and
  acceleration" ref 'Eq. (3.1)' `\mathbf u=d\mathbf r/dt,\quad\mathbf a=d^2\mathbf r/dt^2` live (stretching) "u = 1.649, a = 0.824".
- **Check yourself:** (1) "With the preset 'pure heating' (no wind, H = 0.5 K/h), what do the two thermometers read?"
  — "Both +0.5 K/h: with u = 0 the advective term is zero and DT/Dt = ∂T/∂t." `set {V: 0, H: 0.5}` · (2) "In 'stalled
  front', which thermometer changes and why?" — "The float: the pattern is fixed, so ∂T/∂t = 0 at the probe, but the
  float moves into colder air: DT/Dt = u·∇T = −0.36 K/h." `set {carry: false, dir: 90, H: 0}` · (3) "Double the wind
  speed with the front carried along. Does DT/Dt change?" — "No — it stays 0 (the pattern moves with the air); the
  station's warming doubles to 0.72 K/h." `set {V: 20, carry: true}` · (4) "In the stretching map, what does the probe read
  at t = 0 and t = 3 s?" — "The same u = αx_p: the field is steady; only the particles passing it change." `set {mode:
  'stretch'}`.
- **Selftest parity rows:** `{name: 'local (linear, carried)', js: terms(0, 0, 0, {V: 10, dir: 90, G: 1, H: 0, carry:
  true, front: 'linear'}).local, py: 'ch03.thermal_front_terms(0.0, 0.0, 0.0, 0.0, 10.0, 1e-05, 0.0, 10.0)["local"]',
  rtol: 1e-10}` · `{name: 'advective (front w=200 km)', js: terms(0, 50e3, 3600, {…, front: 'tanh'}).advective, py:
  'ch03.thermal_front_terms(0.0, 50000.0, 3600.0, 0.0, 10.0, 1e-05, 0.0, 10.0, 200000.0)["advective"]', rtol: 1e-10}` ·
  `{name: 'T at a point', js: T(0, 1e5, 0, {…}), py: 'ch03.thermal_front(0.0, 100000.0, 0.0, 1e-05, 0.0, 10.0, 200000.0)',
  rtol: 1e-12}` · `{name: 'stencil total', js: stencilTerms(0, 0, 0, {…}, 10).total, py:
  'ch03.thermal_front_terms(0.0, 0.0, 0.0, 0.0, 10.0, 1e-05, 0.0, 10.0)["total"]', atol: 1e-9}` · `{name: 'x = X e^{αt}',
  js: xmap(2, 1, 0.5), py: 'ch03.lagrangian_map_example(2.0, 1.0, 0.5)', rtol: 1e-12}` · invariant `{name: 'float rate =
  DT/Dt', js: floatRate(2*3600) - terms(float(2*3600)…).total, expect: 0, atol: 1e-9}`.
- **Fit plan:** portrait: `map` (top) + `bars` (bottom), `series` hidden (both thermometer rates in the bars title);
  Explore: V, dir, G (+ H, carry optional); landscape phone: one row (map + bars); notebook frame: rows [1.2, 0.8] with
  `series` below.

### E3 · galilean_frames_cylinder
- **Title:** "Steady or not — does the acceleration care?" · **Summary:** "One slider moves the observer from the still
  lake to the towed cylinder: the streamline picture morphs from unsteady to steady while a tagged particle's
  acceleration stays exactly the same and its local and advective bars trade places." · **CORE:** C05 (also N02, N04,
  N20–N24) · **Reference:** `amplitude_phase_second_order_II_3.html` (windows linked by one state, a numbered live
  derivation in the explanation).
- **meta:** `viz:order 3` · `viz:sections 3.1 3.3` · `viz:equations 3.9` · `viz:fluidpy ch03.cylinder_flow
  ch03.frame_acceleration_terms ch03.galilean_transform ch03.acceleration` · `viz:derivations D06`.
- **Physics (JS ↔ Python):** `ub(x, y, U, a)` body-frame velocity ↔ `ch03.cylinder_flow(x, y, U, a, "body")` ·
  `uObs(X, y, t, s)` = ub(X + (U − U_f)t, y) + (U_f − U, 0) · `terms(X, y, s)` → {u, local, advective, total} by analytic
  derivatives of ub (∂ub/∂x, ∂ub/∂y closed form) ↔ `ch03.frame_acceleration_terms(x, y, U, a, U_frame, t)` · tracer and
  tagged-particle motion by `Viz.num.rk4Step` in the observer frame · streamlines at frozen t by RK4 in arc length from
  14 seeds on the left edge (plus 6 near the body) ↔ `ch03.streamline(…)`.
- **Views** (rows [1.3, 0.8]):
  1. `flow` "The flow as seen by the observer" (`equal: true`; X ∈ [−5, 5], y ∈ [−2.5, 2.5] m) — the cylinder (grey disc)
     at X_c = (U_f − U)t, instantaneous streamlines (teal), 60 tracer particles (faint), the **tagged particle** (black dot)
     with its velocity arrow (black) and **acceleration arrow** (purple, scale bar 1 m/s²), a label "observer moves at
     U_f = … m/s". Pointer: click anywhere to tag a new particle (sets `px, py` at the current t).
  2. `bars` "Acceleration of the tagged particle" — for the x and y components side by side: local ∂u/∂t (blue) +
     advective (u·∇)u (amber) = total (purple), a dashed purple line at the body-frame total (the invariant); title
     "total (0.000, −0.856) m/s² in every frame".
  3. `vectors` "u = U + u′" (`hidePortrait: true`) — at the tagged particle: the body-frame velocity (black), the fluid-
     frame velocity u′ (grey), the observer's velocity (thin), and the arrow U joining them; the book's Fig. 3.8 idea.
- **Controls:** `Uf` "Observer velocity $U_f$" 0…1.5 × U, step 0.01, default 0 m/s (fraction shown) help "0 = still
  lake (fluid frame), U = riding with the cylinder (body frame)" · `U` "Towing speed $U$" 0.2…2, step 0.1, default 1 m/s ·
  `a` "Cylinder radius $a$" 0.5…1.5, default 1 m (*optional*) · `lam` "Scale all velocities ×λ" 0.5…2, default 1
  (*optional*; N22) · transport `t` −4 → 4 s (the cylinder crosses the window in the fluid frame), rate 0.8 s/s,
  `end: 'loop'`.
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **terms** · **transport** · **presets** "fluid
  frame" {Uf: 0, t: 0} · "body frame" {Uf: 1, t: 0} · "halfway" {Uf: 0.5, t: 0} · "far from the body" {px: 0, py: 5} ·
  "traffic light" (caption only: N24 analogy, sets the fluid frame) · "roller coaster" (body frame) · **status**
  ("🟰 steady in this frame: all acceleration is advective" when U_f = U; "⏱ unsteady in this frame — same
  acceleration (0.000, −0.856) m/s²" otherwise) · **inspector** (click any point: "fluid frame: local (…), advective (…);
  body frame: local (0, 0), advective (…); totals equal to 1e-12").
- **Readouts:** "local ∂u/∂t" · "advective (u·∇)u" · "total Du/Dt" (m/s², y components) · "steady?".
- **Explain:**
  0. *What the windows show* — "The **flow** is the ideal flow round a cylinder (Ch. 6) drawn in the frame of an observer
     moving at U_f: teal streamlines are the velocity directions at this instant, the black dot is the particle you
     tagged, its <b class="c-accent">purple</b> arrow is its acceleration. The **bars** split that acceleration into the
     <b class="c-blue">local</b> ∂u/∂t and the <b class="c-amber">advective</b> (u·∇)u part. **u = U + u′** (hidden on
     phones) shows the particle's velocity seen from the two extreme frames."
  1. *The velocity at the particle* — "Body frame: $u=U[1-a^2(x^2-y^2)/r^4]=1\times[1-(0-2.25)/5.0625]=$ **1.444** m/s, v =
     **0**. Observer at U_f = 0: u′ = u − U = **0.444** m/s. The two differ by the constant U — the Galilean
     transformation $\mathbf u(\mathbf x,t)=\mathbf U+\mathbf u'(\mathbf x',t')$."
  2. *The local term in this frame* — "Seen from U_f, the pattern moves past at U − U_f = **1.0** m/s, so
     $\partial\mathbf u/\partial t=(U-U_f)\,\partial\mathbf u_b/\partial x=1.0\times(0,-0.593)=$ **(0, −0.593)** m/s²." (boxed)
  3. *The advective term in this frame* — "$(\mathbf u\cdot\nabla)\mathbf u=(\mathbf u_b\cdot\nabla)\mathbf u_b+(U_f-U)\,\partial\mathbf u_b/\partial x=(0,-0.856)-(0,-0.593)=$
     **(0, −0.263)** m/s²." (boxed)
  4. *The sum (3.9)* — "$\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=$ (0, −0.593) + (0, −0.263) =
     **(0, −0.856)** m/s² — identical to the body frame's $0+(0,-0.856)$. The U terms cancelled exactly (D06 step 8)."
     (boxed)
  5. *Linear vs quadratic* (N22) — "With λ = **2**: local ×2 = (0, −1.185), advective ×4 = (0, −1.053): the advective term
     is quadratic in the velocity, the source of fluid mechanics' nonlinearity."
  6. *At the current time* — "t = `live t` s; the cylinder is at X_c = `live xc` m; the particle at (`live px`,
     `live py`) m."
  7. *Reading the current setting* — U_f = U: "You ride with the cylinder: the flow is steady, streamlines are fixed,
     and every bit of acceleration is advective — like a roller coaster on a fixed track (N24)." · 0 < U_f < U: "An
     in-between observer: the flow is unsteady and both terms share the same total." · U_f = 0: "You stand in the lake:
     the body sweeps past, the water ahead is pushed aside — mostly *local* acceleration, like cars at a traffic light
     that changes (N24)." · U_f > U: "You overtake the cylinder: the local term changes sign; the total still does not."
     Always add: "⚠️ This works because U is constant; a rotating observer adds Coriolis terms (Ch. 4 §4.7)."
- **Derivation tab:** **D06** (9 steps) `view: 'bars'`; goal `set {Uf: 0, U: 1, t: 0, px: 0, py: 1.5}`; step 1
  highlights the moving frame's origin on `flow`; step 3 **live** "local = ∂u′/∂t′ − U·∇′u′ = … − … = (0, −0.593)"; step 5
  **live** "(U + u′)·∇′u′ = …"; step 8 draws the two canceling pieces as hatched blue/amber half-bars "−U·∇′u′" and
  "+U·∇′u′" `watch: "the hatched pieces are equal and opposite"`; step 9 `set {Uf: 0.5}` `watch: "the purple bar does
  not move"`; interpret: "at U_f = … the split is (…, …) and the total is (…) — the same as in the body frame".
- **Code:**
  ```python
  U, a, Uf = {{U}}, {{a}}, {{Uf}}                        # towing speed, radius, observer speed [m/s, m, m/s]
  x, y = {{px}}, {{py}}                                  # tagged particle (observer frame, t = {{t}} s)
  ub = ch03.cylinder_flow(x, y, U, a, frame="body")      # body-frame velocity ({{ubx}}, {{uby}}) m/s
  T = ch03.frame_acceleration_terms(x, y, U, a, Uf, t={{t}})
  print(T["local"])                                      # ∂u/∂t      = ({{lx}}, {{ly}}) m/s²
  print(T["advective"])                                  # (u·∇)u     = ({{ax}}, {{ay}}) m/s²
  print(T["total"])                                      # Du/Dt      = ({{tx}}, {{ty}}) m/s²  (3.9)
  body = ch03.frame_acceleration_terms(x, y, U, a, U, t={{t}})["total"]   # same total in the body frame
  ```
- **Walkthrough (7 steps):** 1. "Tow a cylinder" — "A cylinder is towed at 1 m/s through a still lake. You stand on the
  shore. Is the flow steady? Is the tagged particle accelerating?" `set {Uf: 0, t: -2}`, `play: true` · 2. "Unsteady
  here" — "The body sweeps past: velocities at fixed points change — the blue local bar is large." `set {t: 0}`, `play:
  false`, `terms: true` · 3. "Ride with the body" — "Drag the observer speed to U: the streamlines freeze around the body.
  Steady — the blue bar is gone." `set {Uf: 1}`, `controls: ['Uf']` · 4. "Same arrow" — "The purple arrow did not move. Its
  total is (0, −0.856) m/s² in both frames: only the booking changed." `readouts: ['tot']`, `eq: 'gal'` · 5. "Why" — "The
  frames differ by a constant U. In the lake frame the local term contains −U·∇′u′, the advective term +U·∇′u′; they
  cancel." `derive: {id: 'D06', step: 8}` · 6. "Quadratic" — "Scale every velocity by λ = 2: the local bar doubles, the
  advective bar quadruples (N22)." `set {Uf: 0, lam: 2}` · 7. "Your turn" — "Predict: an observer faster than the
  cylinder (U_f = 1.3 m/s). Which way does the blue bar point? Drag and check." `set {lam: 1}`, `controls: ['Uf']`.
- **Equations:** `gtr` "Galilean transformation" ref 'Fig. 3.8' `\mathbf u(\mathbf x,t)=\mathbf U+\mathbf u'(\mathbf x',t'),\ t=t',\ \mathbf x=\mathbf x'+\mathbf Ut+\mathbf x'_o`
  live "u = 1.000 + 0.444" · `gal` "Invariant acceleration" ref 'Eq. (3.9)'
  `\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=\frac{\partial\mathbf u'}{\partial t'}+(\mathbf u'\cdot\nabla')\mathbf u'`
  live "(0, −0.593) + (0, −0.263) = (0, −0.856)" · `cyl` "Flow past a cylinder (body frame)" ref 'Ch. 6 §6.3'
  `u=U\Big[1-\frac{a^2(x^2-y^2)}{r^4}\Big],\quad v=-\frac{2Ua^2xy}{r^4}` live "(1.444, 0.000) m/s" · `ch` "Chain rule in a moving
  frame" ref 'D06 step 3' `\frac{\partial u_i}{\partial t}=\frac{\partial u'_i}{\partial t'}-U_j\frac{\partial u'_i}{\partial x'_j}` live
  "0 − 1.0 × 0.593".
- **Check yourself:** (1) "At the preset 'body frame', what is the local term, and is the particle accelerating?" —
  "Local = 0 (steady), but the total is (0, −0.856) m/s²: it is all advective — the streamlines curve over the body."
  `set {Uf: 1, t: 0}` · (2) "Tag a particle far from the body (0, 5 m). What happens to all bars?" — "They shrink toward
  zero: far away the flow is almost uniform in every frame." `set {px: 0, py: 5}` · (3) "Set λ = 2 in the fluid frame. By
  what factor does each bar change?" — "Local ×2, advective ×4: linear vs quadratic (N22)." `set {Uf: 0, lam: 2}` · (4)
  "Is the flow steady for an observer at U_f = 0.99 m/s?" — "No — any mismatch with U makes the body move in that frame,
  so ∂u/∂t ≠ 0 near it; only U_f = U is steady." `set {Uf: 0.99}`.
- **Selftest parity rows:** `{name: 'body u at P', js: ub(0, 1.5, 1, 1)[0], py: 'ch03.cylinder_flow(0.0, 1.5, 1.0,
  1.0, "body")[0]', rtol: 1e-12}` · `{name: 'fluid-frame local y', js: terms(0, 1.5, {U: 1, a: 1, Uf: 0, t: 0}).local[1],
  py: 'ch03.frame_acceleration_terms(0.0, 1.5, 1.0, 1.0, 0.0)["local"][1]', rtol: 1e-6}` · `{name: 'advective y (U_f =
  0.5)', js: terms(0, 1.5, {…, Uf: 0.5}).advective[1], py: 'ch03.frame_acceleration_terms(0.0, 1.5, 1.0, 1.0,
  0.5)["advective"][1]', rtol: 1e-6}` · `{name: 'total x at a generic point', js: terms(-1.3, 0.9, {…, Uf: 0.3}).total[0],
  py: 'ch03.frame_acceleration_terms(-1.3, 0.9, 1.0, 1.0, 0.3)["total"][0]', rtol: 1e-6}` · invariants `{name: 'total
  frame-free', js: terms(-1.3, 0.9, {Uf: 0.3}).total[1] - terms(-1.3, 0.9, {Uf: 1}).total[1], expect: 0, atol: 1e-10}` ·
  `{name: 'no flow through the body', js: normalVelocityOnBody(1, 1), expect: 0, atol: 1e-12}`.
- **Fit plan:** portrait: `flow` (top) + `bars` (bottom), `vectors` hidden (not needed for the numbers); Explore: Uf, U
  (a, lam optional); landscape phone: one row flow + bars; notebook: rows [1.3, 0.8] with `vectors` beside `bars`.

### E4 · fluid_element_deformation
- **Title:** "What does each number in S measure?" · **Summary:** "Drag the four entries of the velocity gradient and
  watch a square and a ring of tracers deform; the measured stretching, closing and swelling rates land on the formula
  rates n·S·n, 2S₁₂ and tr G." · **CORE:** C06, C07, C08, C09 (also R01–R03, N25, N27) · **Reference:**
  `forward_noising_lab.html` (click a point to see its exact arithmetic) with `forced_damped_vibrations.html`'s
  explanation panel.
- **meta:** `viz:order 4` · `viz:sections 3.4` · `viz:equations 3.10 3.11 3.12 3.13 3.14` · `viz:fluidpy
  ch03.relative_velocity ch03.linear_strain_rate ch03.shear_strain_rate ch03.volumetric_strain_rate
  ch03.material_volume_ratio ch03.linear_flow_map ch03.measured_strain_rates ch03.velocity_gradient_at` ·
  `viz:derivations D07 D08 D09 D11`.
- **Physics (JS ↔ Python):** `relVel(G, dx)` ↔ `ch03.relative_velocity` · `expm2(G, t)` (closed form, ch02 E3) ↔
  `ch03.linear_flow_map` · `nSn(G, th)` ↔ `ch03.linear_strain_rate(G, [cos, sin])` · `n1Sn2(G, th)` ↔
  `ch03.shear_strain_rate(G, n, n⊥)` · `trG(G)` ↔ `ch03.volumetric_strain_rate` · `areaRatio(G, t)` = det expm2 ↔
  `ch03.material_volume_ratio` · `measured(G, th, dt)` → {stretch, closing, area} ↔ `ch03.measured_strain_rates` ·
  nonlinear field mode: `uNL(x) = (x₁ + 2x₂ + x₁², −x₂ + x₁x₂)` and `gradNL(x)` (central differences) ↔
  `ch03.velocity_gradient_at(…)`.
- **Modes** (`field` chips): **linear u = G·x** (sliders set G) · **C06 nonlinear field** (G at the origin = [[1, 2],
  [0, −1]]; the ring radius slider shows the Taylor remainder).
- **Views** (rows [1.25, 0.85]):
  1. `element` "The element" (`equal: true`; ±1.6 × ring radius) — the square (thin grey at t = 0, filled teal-grey at t),
     a ring of 24 tracers with their relative-velocity arrows du = G·dx (purple; in nonlinear mode the exact du in black
     beside), the **probe thread** at angle θ (blue if stretching, rose if shrinking) and its perpendicular partner (grey),
     the closing corner arc between them. Pointer: drag the probe tip (sets θ); click a tracer (inspector).
  2. `rose` "Stretching rate in every direction" — n·S·n vs θ ∈ [0°, 180°) (teal curve), the measured
     (1/ℓ)Δℓ/Δt of 36 tracked threads over the last Δt (dots), the probe's dot, principal directions (dashed), and a
     second small readout strip "closing rate of the probe's right angle: measured … vs 2n₁·S·n₂ = …".
  3. `area` "Area of the element" (`hidePortrait: true`) — area ratio A(t)/A(0) (black, from the tracked square), ghost
     e^{t tr G} (blue dashed), first order 1 + t tr G (grey); title carries "tr G = … s⁻¹".
- **Controls:** `G11` "$\partial u_1/\partial x_1$" −2…2, step 0.1, default 1 s⁻¹ · `G12` "$\partial u_1/\partial x_2$"
  default 2 · `G21` "$\partial u_2/\partial x_1$" default 0 · `G22` "$\partial u_2/\partial x_2$" default −1 · `theta`
  "Probe direction $\theta$" 0…180°, default 0 · `dt` "Measurement interval $\Delta t$" 10⁻⁴…0.3 s (log), default 0.01
  (*optional*: shows the first-order error) · `ring` "Ring radius" 0.01…1 m, default 0.1 (*optional*, nonlinear mode).
  Transport `t` 0 → 1.5 s, rate 0.5, `end: 'hold'`, card "After t = 1.5 s: area × …, the probe thread × …, the corner
  closed by …°."
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** "pure strain"
  {G: diag(1, −1)} · "simple shear" {G: [[0, 1], [0, 0]]} · "solid-body rotation" {[[0, −1], [1, 0]]} · "expansion"
  {diag(0.5, 0.5)} · "C06 example" {[[1, 2], [0, −1]]} · **status** ("🔄 rigid motion: S = 0 — nothing deforms" /
  "🎈 expanding: tr G = … > 0" / "🫧 compressing: tr G < 0" / "💧 area-preserving: tr G = 0") · **inspector** (click a
  tracer: "dx = (0.071, 0.071) m; $du_1=G_{11}dx_1+G_{12}dx_2=1\times0.071+2\times0.071=$ **0.212** m/s; $du_2=$ …"; in
  nonlinear mode also "exact − linear = (5.0e-3, 5.0e-3) m/s = O(|dx|²)").
- **Readouts:** "n·S·n" · "measured stretch" · "2S₁₂ (probe pair)" · "tr G".
- **Explain:**
  0. *What the windows show* — "The **element** is a small square of fluid and a ring of neighbours around its centre;
     <b class="c-accent">purple</b> arrows are their velocities relative to the centre, du = G·dx (3.10). The probe thread
     is <b class="c-blue">blue</b> while it stretches, <b class="c-rose">rose</b> while it shrinks. The **rose curve** is
     the stretching rate for every direction at once; dots are measured on moving threads. The **area** graph (hidden on
     phones) tracks the element's size."
  1. *Your G split into S and ½R (3.11)–(3.13)* — "$S_{12}=\tfrac12(G_{12}+G_{21})=\tfrac12(2+0)=$ **1.0**, S₁₁ = **1.0**,
     S₂₂ = **−1.0** s⁻¹; $\tfrac12R_{12}=\tfrac12(G_{12}-G_{21})=$ **1.0** s⁻¹ (the spin part, E5). Check: S + ½R = G ✓."
     (table of the four entries)
  2. *Stretching along the probe (linear strain rate)* — "n = (cos 0°, sin 0°): $\mathbf n\cdot\mathbf S\cdot\mathbf n=S_{11}\cos^2\theta+2S_{12}\sin\theta\cos\theta+S_{22}\sin^2\theta=$
     **1.000** s⁻¹ — the thread grows by 1 % of its length per 0.01 s." (boxed)
  3. *Measured on the moving thread* — "Over Δt = 0.01 s: (1/ℓ)Δℓ/Δt = **1.005** s⁻¹ — within first order of the formula
     (shrink Δt and the gap closes like Δt)."
  4. *Closing of the probe's right angle (shear strain rate)* — "$2\,\mathbf n_1\cdot\mathbf S\cdot\mathbf n_2=2\times1.0=$
     **2.000** rad/s; measured **1.99** rad/s. Half of it, S₁₂ = 1.0 s⁻¹, is what the book calls the shear strain rate."
     (boxed)
  5. *Swelling (3.14)* — "$\frac1{\delta V}\frac{D(\delta V)}{Dt}=\frac{\partial u_i}{\partial x_i}=S_{ii}=1+(-1)=$ **0** s⁻¹:
     the area is conserved, although the shape changes. Exact after t: e^{t tr G} = **1.000**." (boxed)
  6. *The Taylor remainder* (nonlinear mode, or hint) — "Ring radius 0.1 m: the largest gap between the exact du and
     G·dx is **7.1 × 10⁻³** m/s; halve the radius and it falls ×4 (second order)."
  7. *At the current time* — "t = `live t` s; area ratio `live A`; probe thread length ratio `live L`; corner angle
     `live ang`°."
  8. *Reading the current setting* — S = 0 (rotation preset): "Rigid motion: every thread keeps its length and the
     right angle stays right; the element only turns — S = 0, and that is why S (not ω) enters the stress law of Ch. 4."
     · tr G = 0, S ≠ 0: "Pure deformation at constant area: the element stretches one way and shrinks the other — like
     water, which barely changes volume." · tr G > 0: "The element swells at tr G per second — rising, expanding air." ·
     tr G < 0: "It shrinks — sinking, compressed air." · probe on a principal axis (|n₁·S·n₂| < 1e-3): "Your probe lies on
     a principal axis: it stretches at an extreme rate and its right angle does not close (E5, C12)."
- **Derivation tab:**
  - **D07** (4 steps) `view: 'element'`; goal `set {field: 'nonlinear', ring: 0.2}`; step 1 **live** "u(x + dx) =
    u(x) + G·dx + (dx₁², dx₁dx₂)"; step 3 `set {ring: 0.05}` `watch: "the black and purple arrows merge as the ring
    shrinks"`; interpret: "at ring radius … the remainder is … m/s, … of du".
  - **D08** (6 steps) `view: 'element'`; goal `set {field: 'linear', G11: 1, G12: 0, G21: 0, G22: -1, theta: 0}`; step 1
    animates A′B′ with the two end arrows; step 3 **live** "(A′B′ − AB)/(AB·dt) = …"; step 6 `set {theta: 45}` `watch: "at
    45° the thread's length does not change: n·S·n = 0"`.
  - **D09** (8 steps) `view: 'element'`; goal `set {G11: 0, G12: 1, G21: 0, G22: 0, theta: 0}`; steps 2–4 draw the angles
    α (vertical side, clockwise) and β (horizontal side, counterclockwise) on the corner; step 6 **live** "½(dα + dβ)/dt =
    ½(1.00 + 0.00) = 0.50 s⁻¹"; step 8 `set {theta: 30}`.
  - **D11** (7 steps) `view: 'area'`; goal `set {G11: 0.5, G12: 0, G21: 0, G22: 0.5}`; step 1 lights the three (here two)
    edges; step 4 **live** "tr G = 0.5 + 0.5 = 1.0 s⁻¹"; step 7 `watch: "the black curve follows the blue ghost e^{t tr G},
    not the grey line"`.
- **Code:**
  ```python
  G = np.array([[{{G11}}, {{G12}}], [{{G21}}, {{G22}}]])   # velocity gradient ∂u_i/∂x_j [1/s]
  S = ch03.strain_rate_tensor(G)                         # (3.12): [[{{S11}}, {{S12}}], [{{S12}}, {{S22}}]]
  n = [np.cos({{th}}), np.sin({{th}})]                    # probe direction θ = {{deg}}°
  print(ch03.linear_strain_rate(G, n))                   # n·S·n = {{nSn}} 1/s
  n2 = [-n[1], n[0]]                                     # the perpendicular thread
  print(2*ch03.shear_strain_rate(G, n, n2))              # closing rate 2 n·S·n2 = {{close}} rad/s
  print(ch03.volumetric_strain_rate(G))                  # (3.14): tr G = {{tr}} 1/s
  print(ch03.material_volume_ratio(G, {{t}}))            # area ratio after t: {{A}}
  print(ch03.measured_strain_rates(G, n, {{dt}}))        # measured: stretch {{ms}}, closing {{mc}}
  ```
- **Walkthrough (7 steps):** 1. "Nine numbers, one blob" — "The purple arrows are how fast each neighbour moves
  relative to the centre: du = G·dx (3.10). What do the entries of G do to the blob?" `set {preset: 'C06'}`, `inspect:
  true` · 2. "Stretch" — "Pure strain: the blue thread along x grows at S₁₁ = 1 s⁻¹ per unit length." `set {preset:
  'pure', theta: 0}`, `derive: {id: 'D08', step: 4}` · 3. "Every direction" — "Drag the probe round: its stretching rate
  traces n·S·n. Dots (measured) sit on the curve." `controls: ['theta']`, highlight `view:rose` · 4. "Close a corner" —
  "Simple shear: the vertical thread tilts, the horizontal does not; the angle closes at 2S₁₂ = 1 rad/s." `set {preset:
  'shear', theta: 0}`, `derive: {id: 'D09', step: 6}` · 5. "Swell" — "Expansion: area grows at tr G = 1 s⁻¹ — (3.14),
  $\frac1{\delta V}\frac{D(\delta V)}{Dt}=S_{ii}$." `set {preset: 'expand'}`, `play: true`, `readouts: ['tr']` · 6. "Spin only" —
  "Solid-body rotation: every rate is zero. S = 0: rigid motion does not deform (N27)." `set {preset: 'rotation'}`,
  `play: false` · 7. "Your turn" — "Set G₁₂ = 2, G₂₁ = 0. Predict the direction of fastest stretching, then find it with
  the probe." `set {preset: 'C06'}`, `controls: ['G12', 'theta']`.
- **Equations:** `rel` "Relative velocity" ref 'Eq. (3.10)' `du_i=(\partial u_i/\partial x_j)\,dx_j` live "(0.212, −0.071)
  m/s at the clicked tracer" · `split` "Split of G" ref 'Eqs. (3.11)–(3.13)'
  `\frac{\partial u_i}{\partial x_j}=S_{ij}+\tfrac12R_{ij},\ S_{ij}=\tfrac12\Big(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\Big),\ R_{ij}=\frac{\partial u_i}{\partial x_j}-\frac{\partial u_j}{\partial x_i}`
  (aligned over three rows) · `lin` "Linear strain rate" ref '§3.4'
  `\frac1{\delta x_1}\frac{D}{Dt}(\delta x_1)=\frac{\partial u_1}{\partial x_1}` live "S₁₁ = 1.0" · `shr` "Shear strain rate" ref
  '§3.4' `\tfrac12\frac{D(\alpha+\beta)}{Dt}=\tfrac12\Big(\frac{\partial u_1}{\partial x_2}+\frac{\partial u_2}{\partial x_1}\Big)=S_{12}`
  live "½(2 + 0) = 1.0" · `vol` "Volumetric strain rate" ref 'Eq. (3.14)'
  `\frac1{\delta V}\frac{D}{Dt}(\delta V)=\frac{\partial u_i}{\partial x_i}=S_{ii}` live "1 + (−1) = 0".
- **Check yourself:** (1) "With simple shear, which probe direction stretches fastest and at what rate?" — "45°, at
  S₁₂ = 0.5 s⁻¹ (γ/2); at 135° it shrinks at the same rate." `set {preset: 'shear'}` · (2) "In solid-body rotation, why
  do the dots all sit at zero?" — "Rigid rotation keeps every length: S = 0 (N27); the arrows are all tangential." `set
  {preset: 'rotation'}` · (3) "Expansion with tr G = 1 s⁻¹: after 1 s is the area 2 or 2.72 times larger?" — "2.72 =
  e^{1}: (3.14) is a rate per current area, so the growth compounds; 1 + t tr G is only the first-order tangent." `set
  {preset: 'expand', t: 1}` · (4) "Increase Δt to 0.3 s. Why do the measured dots drift off the curve?" — "The
  measurement uses finite displacements; the formula is the Δt → 0 limit (the D08 limit step)." `set {dt: 0.3}`.
- **Selftest parity rows:** `{name: 'n·S·n at 45°', js: nSn([[1, 2], [0, -1]], Math.PI/4), py:
  'ch03.linear_strain_rate([[1.0, 2.0], [0.0, -1.0]], [0.7071067811865476, 0.7071067811865476])', rtol: 1e-12}` ·
  `{name: 'S12', js: n1Sn2([[1, 2], [0, -1]], 0), py: 'ch03.shear_strain_rate([[1.0, 2.0], [0.0, -1.0]], [1.0, 0.0], [0.0,
  1.0])', rtol: 1e-12}` · `{name: 'area ratio', js: areaRatio([[0.5, 0.3], [-0.2, 0.5]], 1.0), py:
  'ch03.material_volume_ratio([[0.5, 0.3], [-0.2, 0.5]], 1.0)', rtol: 1e-10}` · `{name: 'expm2 entry', js: expm2([[0, 1],
  [-1, 0]], 0.7)[0][1], py: 'ch03.linear_flow_map([[0.0, 1.0], [-1.0, 0.0]], 0.7)[0][1]', rtol: 1e-10}` · `{name:
  'measured stretch', js: measured([[1, 2], [0, -1]], 0.3927, 1e-3).stretch, py: 'ch03.measured_strain_rates([[1.0, 2.0],
  [0.0, -1.0]], [0.9238795325112867, 0.3826834323650898], 0.001)["stretch"]', rtol: 1e-8}` · `{name: 'du', js:
  relVel([[1, 2], [0, -1]], [0.01, 0.02])[0], py: 'ch03.relative_velocity([[1.0, 2.0], [0.0, -1.0]], [0.01, 0.02])[0]',
  rtol: 1e-12}` · invariant `{name: 'S + ½R = G', js: maxAbs(add(S(G), half(R(G))), G), expect: 0, atol: 1e-15}`.
- **Fit plan:** portrait: `element` + `rose`; `area` hidden, tr G in the `rose` title and the status; Explore: G11, G12,
  G21, G22, theta (dt, ring optional); landscape phone: one row; notebook frame: `element` | `rose` over `area`.

### E5 · spin_and_principal_axes
- **Title:** "Can a straight flow make a fluid element spin?" · **Summary:** "Two perpendicular threads and a paddle
  wheel ride in a flow you choose: single threads turn at different rates, every pair averages ½ω₃; a probe's
  relative velocity splits into strain and rigid rotation; a circle becomes an ellipse on the principal axes." ·
  **CORE:** C10, C11, C12 (also R05, N28, N30, N31, N32, N34, N35) · **Reference:** `amplitude_phase_second_order_II_3.html`
  (system + graphs linked by one state; numbered live derivation in the explanation).
- **meta:** `viz:order 5` · `viz:sections 3.4 3.5` · `viz:equations 3.15 3.16 3.19 3.20 3.21` · `viz:fluidpy
  ch03.material_line_rotation_rate ch03.element_rotation_rate ch03.relative_velocity_split ch03.principal_strain_rates
  ch03.strain_ellipse_axes ch03.vorticity_in_rotating_frame ch03.parallel_shear_kinematics` · `viz:derivations D12 D13
  D14 D15 D16`.
- **Physics (JS ↔ Python):** `lineRate(G, th)` = e_θ·G·e(θ) ↔ `ch03.material_line_rotation_rate` · `spin(G)` = ½ω₃ ↔
  `ch03.element_rotation_rate(G)[2]` · `split(G, dx)` → {du, strain, rot} ↔ `ch03.relative_velocity_split` ·
  `principal(G)` (closed-form 2 × 2 symmetric eigen, λ ascending) ↔ `ch03.principal_strain_rates` · `ellipse(G, t,
  method)` (first order / strain only / exact via 2 × 2 SVD of expm2) ↔ `ch03.strain_ellipse_axes` · `omegaPrime(w3, Om)`
  ↔ `ch03.vorticity_in_rotating_frame` · element motion by `expm2(G, t)`; in the rotating-observer option every drawn
  vector is rotated back by −Ωt.
- **Modes** (`flow` chips): **parallel shear** G = [[0, γ], [0, 0]] · **solid body** [[0, −γ], [γ, 0]] · **pure strain**
  [[γ/2, 0], [0, −γ/2]] · **general** [[0.4γ, γ], [−0.3γ, −0.4γ]].
- **Views** (rows [1.2, 0.9]):
  1. `element` "The element" (`equal: true`) — the pair of perpendicular material threads at angle θ (teal and teal
     dashed) and a small **paddle wheel** at the centre turning at ½ω₃ (orange); the principal axes of S (blue stretch,
     rose squeeze, dashed); a circle of tracers becoming its ellipse (thin black) with the SVD long axis (purple); the
     **probe** on the ring with three arrows (du purple, S·dx teal, ½ω × dx orange). Pointer: drag a thread end (sets θ),
     drag the probe (sets φ), click a ring point (inspector).
  2. `rates` "How fast each thread turns" — θ̇(θ) for a single thread over θ ∈ [0°, 180°] (teal), the partner curve
     θ̇(θ + 90°) (teal dashed), their average (orange, flat at ½ω₃), the pair's two dots and their midpoint; title
     "pair average −0.500 rad/s = ω₃/2".
  3. `split` "du = S·dx + ½ω × dx" (`hidePortrait: true`) — the probe's three arrows enlarged with their numbers, plus a
     small table of the ellipse semi-axes after t by the three methods; title "strain (0.50, 0) + rotation (0.50, 0)".
- **Controls:** `gamma` "Rate $\gamma$" 0.1…2, step 0.05, default 1 s⁻¹ · `theta` "Pair angle $\theta$" 0…180°, default 30 ·
  `phi` "Probe direction" 0…360°, default 90 · `Om` "Observer rotation $\Omega$" −2…2, step 0.05, default 0 rad/s
  (*optional*) · `method` select first order / strain only / exact (*optional*). Transport `t` 0 → 3 s, rate 0.5,
  `end: 'hold'`, card "After 3 s the wheel turned … rad (½ω₃t); the ellipse's long axis is at …° (principal axis
  45°)."
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **modes** (4 flows) ·
  **presets** "shear γ = 1" · "solid body" · "pure strain" · "rotate with the element (Ω = ω₃/2)" {Om: ½ω₃} · "Fig. 3.14
  pair (θ = 45°)" · **status** ("🌀 spins and strains (shear): ½ω₃ = −0.50 rad/s, principal ±0.50 s⁻¹" / "🔄 spins
  only: S = 0" / "↔ strains only: ω = 0"; with Ω ≠ 0 appended "· observer at Ω: ω′₃ = ω₃ − 2Ω = …") · **inspector**
  (click a ring point: "dx = (0, 0.1) m: S·dx = (0.050, 0.000), ½ω×dx = ½(0, 0, −1)×(0, 0.1, 0) = (0.050, 0.000), du =
  (0.100, 0.000) m/s").
- **Readouts:** "½ω₃ (spin)" · "θ̇ thread 1" · "θ̇ thread 2" · "principal rates" · "ω′₃ (observer)".
- **Explain:**
  0. *What the windows show* — "In the **element**, two <b class="c-teal">teal</b> threads start perpendicular; the
     <b class="c-orange">orange</b> paddle wheel turns at the element's spin; dashed <b class="c-blue">blue</b>/<b
     class="c-rose">rose</b> lines are the principal axes of S. **How fast each thread turns** shows single threads
     (teal) and the pair average (orange). **The split** (hidden on phones) enlarges the probe's arrows."
  1. *Vorticity from your G (3.16)* — "$\omega_3=\frac{\partial u_2}{\partial x_1}-\frac{\partial u_1}{\partial x_2}=0-1=$
     **−1.0** s⁻¹ — clockwise. In matrix form (3.15), $R_{ij}=-\varepsilon_{ijk}\omega_k$: R₂₁ = ω₃ = −1.0."
  2. *The two threads* — "A thread at θ turns at $\dot\theta=\mathbf e_\theta\cdot\mathbf G\cdot\mathbf e(\theta)=-\gamma\sin^2\theta$:
     at 30°, **−0.250** rad/s; its partner at 120°, **−0.750** rad/s."
  3. *Their average* — "$\tfrac12(-0.250-0.750)=$ **−0.500** rad/s $=\omega_3/2$ — and the same for every θ (drag it):
     the element's spin is half the vorticity." (boxed)
  4. *Principal strain axes (3.20)–(3.21)* — "S = [[0, 0.5], [0.5, 0]]: eigenvalues **±0.500** s⁻¹ along 45° and 135°. In
     that frame $d\bar u_\alpha=\bar S_{\alpha\alpha}d\bar x_\alpha$: stretching along 45°, compression (eigenvalue −γ/2) along
     135°, no shear." (boxed)
  5. *The probe's split (3.19)* — "dx = (0, 0.1) m: $d\mathbf u=\mathbf G\cdot d\mathbf x=$ (0.100, 0) = S·dx (0.050, 0) +
     ½ω × dx (0.050, 0) m/s. The orange part is a rigid rotation at ω/2 — it never changes the distance to the
     neighbour." (boxed)
  6. *The ellipse after t* — "First order (D16): semi-axes 1 ± 0.5t = **1.05, 0.95** at t = 0.1 s along 45°; strain alone
     e^{±0.05} = 1.051, 0.951; exact linear flow 1.0513, 0.9513 at **43.6°** — the rotation part has turned it."
  7. *A rotating observer* (or hint "set Ω in Explore") — "$\omega'_3=\omega_3-2\Omega=-1.0-2\times(-0.5)=$ **0**: rotating
     with the element, the wheel stands still."
  8. *At the current time* — "t = `live t` s; the wheel has turned `live turned` rad; long axis at `live ang`°."
  9. *Reading the current setting* — shear: "Straight streamlines, yet every element spins clockwise at γ/2 and
     stretches along 45°: half of the shear is rotation, half is strain." · solid body: "Every thread turns at the same
     rate γ: the element spins and does not deform (S = 0)." · pure strain: "No spin: threads along the principal axes do
     not turn at all, others turn toward the stretching axis; the circle becomes an ellipse whose axes never move." ·
     Ω = ω₃/2: "You rotate with the element: ω′ = 0, the wheel stops — vorticity depends on the frame (unlike S)." · on
     Earth: "the planetary part 2Ω sin φ = f is added the same way (Ch. 13)."
- **Derivation tab:**
  - **D12** (8 steps) `view: 'element'`; goal `set {flow: 'general', theta: 0}`; steps 1–3 draw α (vertical thread) and β
    (horizontal thread) with their signs; step 4 **live** "½(−(∂u₁/∂x₂) + ∂u₂/∂x₁) = ½(−1.0 − 0.3) = −0.65 rad/s"; step 6
    `watch: "the orange wheel turns at exactly this rate"`.
  - **D13** (6 steps) `view: 'element'`; goal `set {flow: 'solid', Om: 0}`; step 5 `set {Om: 1}` `watch: "the wheel and the
    element stand still: ω′ = 0"`; step 6 **live** "ω′₃ = 2.0 − 2 × 1.0 = 0".
  - **D14** (6 steps) `view: 'rates'`; goal `set {flow: 'shear', gamma: 1, theta: 30}`; step 3 **live** "−γ sin²30° =
    −0.25"; step 6 `play` sweeps θ 0 → 180° `watch: "the orange midpoint never leaves −0.5"`.
  - **D15** (6 steps) `view: 'element'`; goal `set {flow: 'shear', phi: 90}`; step 4 (the ε swap) flips the orange arrow
    on screen for one second (drawn reversed, then corrected) `watch: "one index swap flips the sign"`; step 6 **live**
    "(0.05, 0) + (0.05, 0) = (0.10, 0)".
  - **D16** (8 steps) `view: 'element'`; goal `set {flow: 'shear', t: 0}`; step 2 rotates the drawing onto the principal
    frame (axes at 45°); step 7 **live** "a = 1 + 0.5 × 0.1 = 1.05, b = 0.95"; step 8 `set {t: 2}` `watch: "the purple
    long axis leaves the blue 45° line: first order only"`.
- **Code:**
  ```python
  G = ch03.velocity_gradient_preset("simple_shear", Gamma={{gamma}})   # G = [[0, γ], [0, 0]]
  w3 = ch03.vorticity_from_gradient(G)[2]                # (3.16): ω₃ = {{w3}} 1/s
  th = np.deg2rad({{theta}})                             # pair angle
  r1 = ch03.material_line_rotation_rate(G, th)           # thread 1: {{r1}} rad/s
  r2 = ch03.material_line_rotation_rate(G, th + np.pi/2) # thread 2: {{r2}} rad/s
  print(0.5*(r1 + r2), 0.5*w3)                           # both {{avg}}: spin = ω₃/2
  du, du_s, du_r = ch03.relative_velocity_split(G, [{{dx1}}, {{dx2}}])   # (3.19)
  lam, axes = ch03.principal_strain_rates(G)             # ±{{lam}} 1/s along 45°, 135°
  a, d = ch03.strain_ellipse_axes(G, {{t}}, method="{{method}}")        # semi-axes {{a1}}, {{a2}}
  print(ch03.vorticity_in_rotating_frame(w3, {{Om}}))    # ω′₃ = {{wp}}
  ```
- **Walkthrough (8 steps):** 1. "A paddle wheel in a river" — "Near the bank the water is slower: straight streamlines.
  Does a small paddle wheel turn?" `set {flow: 'shear', gamma: 1, t: 0}`, `play: true` · 2. "Threads disagree" — "The
  vertical thread turns at −1 rad/s, the horizontal at 0. Which one is 'the' rotation?" `set {theta: 0, t: 0}`, `play:
  false`, `readouts: ['r1', 'r2']` · 3. "Average a pair" — "Take the average of two perpendicular threads: ½(−α̇ + β̇) =
  ½ω₃ = −0.5 rad/s." `derive: {id: 'D12', step: 6}` · 4. "Any pair" — "Drag the pair angle: the teal dots move, their
  midpoint stays on the orange line." `controls: ['theta']`, highlight `view:rates` · 5. "Split one arrow" — "Drag the
  probe: du (purple) = S·dx (teal) + ½ω × dx (orange) — (3.19), $du_i=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i$."
  `controls: ['phi']`, `eq: 'split'`, `inspect: true` · 6. "Ellipse on the principal axes" — "The circle stretches along
  the blue 45° axis at +0.5 s⁻¹ and shrinks along 135°." `set {t: 0.4}`, `eq: 'princ'` · 7. "Turn with it" — "Rotate the
  observer at Ω = ω₃/2: the wheel stops. S did not change; ω did." `set {Om: -0.5}`, `derive: {id: 'D13', step: 5}` · 8.
  "Your turn" — "Choose pure strain. Predict: does the wheel turn? Which threads do not turn at all? Then press ▶." `set
  {flow: 'pure', Om: 0, t: 0}`, `controls: ['theta', 'phi']`.
- **Equations:** `spin` "Element rotation rate" ref '§3.4'
  `\tfrac12\frac{D(-\alpha+\beta)}{Dt}=\tfrac12\Big(-\frac{\partial u_1}{\partial x_2}+\frac{\partial u_2}{\partial x_1}\Big)=-\frac{R_{12}}2=\frac{R_{21}}2`
  live "½(−1 + 0) = −0.5" · `rw` "R ↔ ω" ref 'Eq. (3.15)' `R_{ij}=-\varepsilon_{ijk}\omega_k=\begin{bmatrix}0&-\omega_3&\omega_2\\\omega_3&0&-\omega_1\\-\omega_2&\omega_1&0\end{bmatrix}`
  live "ω₃ = −1" · `split` "Deformation + rotation" ref 'Eq. (3.19)'
  `du_i=\Big(S_{ij}-\tfrac12\varepsilon_{ijk}\omega_k\Big)dx_j=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i` live
  "(0.05, 0) + (0.05, 0)" · `princ` "Principal frame" ref 'Eqs. (3.20)–(3.21)'
  `d\bar{\mathbf u}=\bar{\mathbf S}\cdot d\bar{\mathbf x},\quad d\bar u_\alpha=\bar S_{\alpha\alpha}\,d\bar x_\alpha` live "S̄ = diag(0.5,
  −0.5)" · `rot` "Rotating observer" ref 'N28' `\omega'_z=\omega_z-2\Omega` live "−1 − 2(0) = −1".
- **Check yourself:** (1) "In the shear flow, which thread direction does not turn at all?" — "The one along the flow
  (θ = 0): θ̇ = −γ sin²0 = 0; its partner turns at −γ." `set {flow: 'shear', theta: 0}` · (2) "In pure strain, what is the
  pair average, and which threads keep their direction?" — "0 (no vorticity); threads on the principal axes (0° and 90°
  here) do not turn." `set {flow: 'pure'}` · (3) "Solid body with γ = 1: what Ω makes the wheel stop?" — "Ω = 1 rad/s:
  ω₃ = 2, ω′ = 2 − 2 × 1 = 0." `set {flow: 'solid', gamma: 1}` · (4) "Why is the exact ellipse's axis at 43.6°, not 45°,
  after 0.1 s of shear?" — "The rotation part (−0.5 rad/s) turns the ellipse while it forms; 45° is exact only at the
  first instant (D16 step 8)." `set {flow: 'shear', method: 'exact', t: 0.1}`.
- **Selftest parity rows:** `{name: 'θ̇ at 30°', js: lineRate([[0, 1], [0, 0]], Math.PI/6), py:
  'ch03.material_line_rotation_rate([[0.0, 1.0], [0.0, 0.0]], 0.5235987755982988)', rtol: 1e-12}` · `{name: 'spin',
  js: spin([[0.4, 1], [-0.3, -0.4]]), py: 'ch03.element_rotation_rate([[0.4, 1.0], [-0.3, -0.4]])[2]', rtol: 1e-12}` ·
  `{name: 'rotation part x', js: split([[0, 1], [0, 0]], [0, 1]).rot[0], py: 'ch03.relative_velocity_split([[0.0, 1.0],
  [0.0, 0.0]], [0.0, 1.0])[2][0]', rtol: 1e-12}` · `{name: 'principal λ₂', js: principal([[0, 1], [0, 0]]).lam[1], py:
  'ch03.principal_strain_rates([[0.0, 1.0], [0.0, 0.0]])[0][1]', rtol: 1e-12}` · `{name: 'exact semi-axis', js:
  ellipse([[0, 1], [0, 0]], 0.1, 'exact').a[0], py: 'ch03.strain_ellipse_axes([[0.0, 1.0], [0.0, 0.0]], 0.1,
  "exact")[0][0]', rtol: 1e-10}` · `{name: 'ω′', js: omegaPrime(2, 1), py: 'ch03.vorticity_in_rotating_frame(2.0, 1.0)',
  atol: 1e-15}` · `{name: 'shear ω₃', js: -1, py: 'ch03.parallel_shear_kinematics(1.0)["omega3"]', rtol: 1e-12}` ·
  invariant `{name: 'pair average', js: maxOver(th => 0.5*(lineRate(Gs, th) + lineRate(Gs, th + Math.PI/2)) + 0.5),
  expect: 0, atol: 1e-14}`.
- **Fit plan:** portrait: `element` + `rates`; `split` hidden (its numbers in the inspector and the `element` title);
  Explore: flow chips, gamma, theta, phi (Om, method optional); landscape phone one row; notebook frame three views.

### E6 · vortex_paddle_wheels
- **Title:** "Going round in circles ≠ spinning" · **Summary:** "Solid-body, irrotational, Rankine and Gaussian
  vortices on one clock: paddle wheels orbit and turn at half the local vorticity, a draggable loop measures the
  circulation, and the u_θ(r), ω_z(r) profiles show where the core is." · **CORE:** C13, C14 (also R07, N29, N36–N44) ·
  **Reference:** `angular_frequency_explorer_1.html` (modes, linked views on one clock, "Right now" notes with a
  highlighted table).
- **meta:** `viz:order 6` · `viz:sections 3.5` · `viz:equations 3.22 3.23 3.24 3.25 3.26 3.27 3.28 3.29` · `viz:fluidpy
  ch03.vortex_profile ch03.rankine_vortex ch03.gaussian_vortex ch03.gaussian_vortex_max_radius ch03.polar_vorticity_z
  ch03.circulation_circle ch03.mean_vorticity_in_disc ch03.annular_sector_circulation` · `viz:derivations D17 D18 D19
  D20`.
- **Physics (JS ↔ Python):** `profile(kind, r, G, s)` → [u_θ, ω_z] ↔ `ch03.vortex_profile(kind, r, Gamma, sigma)` ·
  `rMax(s)` (`Viz.num.brentq` on 1 + 2x − eˣ in (0.5, 3)) ↔ `ch03.gaussian_vortex_max_radius(sigma)` · `polarVort(kind, r)`
  (central difference of r u_θ) ↔ `ch03.polar_vorticity_z` · `circ(kind, R, centre)` (512-point midpoint loop sum of
  u·t ds) ↔ `ch03.circulation_circle(kind, R, center=…, Gamma=, sigma=)` · `sector(kind, r, dr, dth)` (four legs) ↔
  `ch03.annular_sector_circulation` · `meanVort(kind, R)` ↔ `ch03.mean_vorticity_in_disc` · wheels and tracers: exact
  circular motion at angular rate u_θ(r)/r; wheel heading θ_w(t) = θ_w0 + ½ω_z(r) t.
- **Modes** (`kind` chips): **solid body** (ω₀ = Γ/2πσ²) · **line vortex** (B = Γ/2π) · **Rankine** · **Gaussian**.
- **Views** (rows [1.25, 0.85]):
  1. `plane` "The vortex" (`equal: true`; x, y ∈ [−3σ, 3σ]) — faint orange shading of ω_z, the core circle r = σ
     (dashed), 150 tracers, **paddle wheels** at r = 0.5σ, 1σ, 1.5σ, 2.5σ (orange crosses) carried round and turned at
     ½ω_z, the **loop** (purple circle or annular sector ABCD) with its Γ label; a marker × on the axis in line-vortex mode
     ("ω = ∞ here"). Pointer: drag the loop centre; drag its rim (sets radius); click a point (inspector).
  2. `profiles` "u_θ(r) and ω_z(r)" — u_θ (teal) and ω_z (orange) vs r/σ ∈ [0, 4] with the other kinds as grey ghosts,
     the core edge (dashed), the peak marker (at 1 or 1.1209), and the loop radius marker; title "peak u_θ = 0.638 m/s at
     r = 1.121 σ".
  3. `circ` "Circulation inside radius r" (`hidePortrait: true`) — Γ(r)/Γ vs r/σ (purple) with the 2πB line (line
     vortex) and, on a second log–log axis strip, the mean vorticity in a disc Γ(r)/πr² (slope −2 for the line vortex).
- **Controls:** `Gamma` "Circulation $\Gamma$" 0.5…20, step 0.1, default 6.283 m²/s · `sigma` "Core size $\sigma$"
  0.2…2, step 0.05, default 1 m · `R` "Loop radius" 0.05…3σ, default 1σ · `shape` chips circle / sector (*optional*) ·
  `realm` chips lab / bathtub / tornado / tropical cyclone (*optional*; sets σ and Γ in physical units from the table, the
  picture stays in r/σ). Transport `t` 0 → 4 × (2πσ²/Γ) (four core rotations), rate 1, `end: 'loop'`.
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (4 vortices) · **presets**
  "bathtub" {realm: 'bath'} · "tornado" · "tropical cyclone" · "loop off the axis" {loop centre (1.5σ, 0), R 0.5σ} ·
  "loop at the core edge" {R: σ} · **transport** · **status** (wheel or loop based: "🔄 inside the core: the wheel spins
  (½ω_z = …) as it orbits" / "🧭 irrotational here: the wheel keeps its heading" / "⚠️ the loop encloses the singular
  core: Γ = 2πB although ω = 0 on the loop") · **inspector** (click a point at r: "u_θ = (Γ/2πr)(1 − e^{−r²/σ²}) = … m/s;
  ω_z = (Γ/πσ²)e^{−r²/σ²} = … s⁻¹; spin ½ω_z = … vs orbit rate u_θ/r = … rad/s") · **notes** (a real-vortex table with
  the current row highlighted: bathtub σ 5 mm, u_max 0.3 m/s, Γ 0.015 m²/s · tornado 100 m, 60 m/s, 5.9 × 10⁴ ·
  tropical cyclone 36 km, 50 m/s, 1.8 × 10⁷ — "our order-of-magnitude numbers").
- **Readouts:** "u_θ at the loop" · "Γ of the loop" · "ω_z at the loop" · "wheel spin / orbit".
- **Explain:**
  0. *What the windows show* — "The **vortex**: tracers go round; <b class="c-orange">orange</b> paddle wheels turn at
     half the local vorticity (shaded orange); the <b class="c-accent">purple</b> loop measures the circulation. The
     **profiles** show the speed u_θ (<b class="c-teal">teal</b>) and vorticity ω_z (orange) against r/σ. **Circulation
     inside r** (hidden on phones) shows how Γ(r) builds up and how the mean vorticity in a disc behaves."
  1. *Speed at the loop radius* — Gaussian: "$u_\theta=\frac{\Gamma}{2\pi r}(1-e^{-r^2/\sigma^2})=\frac{6.283}{2\pi\times1}(1-e^{-1})=$
     **0.632** m/s at r = 1 m." (the formula of the chosen kind: (3.22), (3.25), (3.28) or (3.29) written out)
  2. *Circulation of the loop* — "Centred circle: $\Gamma=\oint\mathbf u\cdot d\mathbf s=2\pi ru_\theta=2\pi\times1\times0.632=$
     **3.972** m²/s = 63 % of the total Γ = 6.283 (by Stokes (3.18), the vorticity inside the loop)." (boxed)
  3. *Vorticity by (3.23)* — "$\omega_z=\frac1r\frac{d}{dr}(ru_\theta)=$ **0.736** s⁻¹ at r = 1 m ($=\frac{\Gamma}{\pi\sigma^2}e^{-1}$)."
     (boxed)
  4. *Spin vs orbit* — "A wheel at r = 1 m orbits at u_θ/r = **0.632** rad/s and spins at ½ω_z = **0.368** rad/s — it
     turns, but more slowly than it goes round. (Solid body: equal; line vortex: spin 0.)"
  5. *Where the wind peaks* — "Gaussian: $1+2r^2/\sigma^2=e^{r^2/\sigma^2}$ at r/σ = **1.1209**, u_max = **0.638** m/s;
     Rankine: at r = σ, u_max = Γ/2πσ = 1.000 m/s." (boxed)
  6. *Mean vorticity in the loop (3.27)* — "Γ(r)/πr² = **1.264** s⁻¹; for the line vortex it is 2B/r², which grows
     without bound as the loop shrinks — the δ core."
  7. *Off-centre loop / sector* (or hint) — "Loop centred at (1.5, 0) m, radius 0.5 m: Γ = **…** m²/s (line vortex: 0 —
     no vorticity inside; solid body: 2ω₀ × π(0.5)² = 1.571 m²/s for ω₀ = 1)."
  8. *Real vortices* (when `realm` ≠ lab) — "Tropical cyclone: σ = 36 km, Γ = 1.8 × 10⁷ m²/s: radius of maximum wind
     1.1209σ = **40** km, u_max = **50** m/s, core vorticity Γ/πσ² = **4.4 × 10⁻³** s⁻¹ ≈ 40 f at mid-latitudes."
  9. *At the current time* — "t = `live t` s; the inner wheel has orbited `live orb`° and turned `live spin`°."
  10. *Reading the current setting* — solid body: "Everything spins at the orbit rate: the tank turns as one body, no
      element deforms." · line vortex: "Elements race round near the axis yet never turn; all the vorticity sits on the
      axis. Any loop around it has Γ = 2πB, any loop beside it Γ = 0 — the region with a hole where ω = 0 does not give a
      single-valued potential (N29)." · Rankine/Gaussian with loop inside the core: "Inside the core the fluid turns
      almost like a solid body." · outside: "Outside the core the flow is almost irrotational and Γ has reached its
      full value."
- **Derivation tab:**
  - **D17** (10 steps) `view: 'plane'`; goal `set {kind: 'gaussian', shape: 'sector', R: 0.8}`; steps 2, 3, 5, 6 light each
    leg of the sector in turn (outer arc, inner arc, the two radial legs) with its contribution; step 4 **live** "[u_θ r]
    at r + dr minus at r = … × dθ"; step 9 **live** "ω_z = Γ_sector/(r dr dθ) = …, formula … s⁻¹"; step 10 `set {kind:
    'solid'}` `watch: "every wheel turns at the orbit rate: 2ω₀/2 = ω₀"`.
  - **D18** (5 steps) `view: 'profiles'`; goal `set {kind: 'rankine', R: 1}`; step 3 puts the loop at r = σ `watch: "the
    teal curve is continuous at the dashed line, the orange one jumps"`; step 5 **live** "u_max = Γ/2πσ = 1.000 m/s".
  - **D19** (7 steps) `view: 'circ'` (phones: `profiles`); goal `set {kind: 'gaussian'}`; step 1 `play` sweeps the loop
    radius outward; step 3 **live** "Γ(r) = Γ(1 − e^{−r²/σ²}) = …"; step 7 `watch: "far out the teal curve joins the grey
    line-vortex ghost"`.
  - **D20** (7 steps) `view: 'profiles'`; goal `set {kind: 'gaussian'}`; step 6 **live** "brentq: x* = 1.25643, r* =
    1.12091σ"; step 7 `watch: "the peak marker sits at r/σ = 1.1209 for every σ"`.
- **Code:**
  ```python
  kind, Gam, sig, R = "{{kind}}", {{Gamma}}, {{sigma}}, {{R}}   # vortex, Γ [m²/s], σ [m], loop radius [m]
  ut, wz = ch03.vortex_profile(kind, R, Gam, sig)        # u_θ = {{ut}} m/s, ω_z = {{wz}} 1/s
  G_loop = ch03.circulation_circle(kind, R, Gamma=Gam, sigma=sig)   # 2πR u_θ = {{GL}} m²/s
  w_num = ch03.polar_vorticity_z(lambda r, t: 0.0,
      lambda r, t: ch03.vortex_profile(kind, r, Gam, sig)[0], R, 0.0)  # (3.23): {{wn}} 1/s
  print(0.5*wz, ut/R)                                    # spin {{spin}} vs orbit {{orb}} rad/s
  rstar = ch03.gaussian_vortex_max_radius(sig)           # 1.1209σ = {{rstar}} m
  print(ch03.mean_vorticity_in_disc(kind, R, Gamma=Gam, sigma=sig))  # Γ(R)/πR² = {{mv}} 1/s
  ```
- **Walkthrough (7 steps):** 1. "Two tanks" — "One tank spun like a merry-go-round, one draining. Both go round. Put a
  paddle wheel in each: which turns?" `set {kind: 'solid'}`, `play: true` · 2. "Merry-go-round" — "Solid body: each wheel
  turns once per orbit. $\omega_z=\frac1r\frac{d}{dr}(r\cdot\omega_0r)=2\omega_0$ (3.23) everywhere." `eq: 'wz'`,
  `readouts: ['spin']` · 3. "Drain" — "Line vortex: wheels race round near the axis but keep pointing the same way —
  ω_z = 0 off the axis." `set {kind: 'line'}` · 4. "Where did the vorticity go?" — "Every centred loop has Γ = 2πB (3.26);
  shrink it and Γ/πr² = 2B/r² explodes — a δ core. Move the loop off the axis: Γ = 0." `set {R: 0.3}`, `controls:
  ['R']`, `play: false` · 5. "A real vortex" — "Rankine: solid-body core, irrotational outside, peak at r = σ." `set {kind:
  'rankine', R: 1}`, `derive: {id: 'D18', step: 3}` · 6. "Smooth core" — "Gaussian: the peak moves out to 1.1209σ,
  where $1+2r^2/\sigma^2=e^{r^2/\sigma^2}$." `set {kind: 'gaussian'}`, `derive: {id: 'D20', step: 6}` · 7. "Your turn" —
  "Pick 'tropical cyclone'. Predict the radius of maximum wind from σ = 36 km, then check in Explain." `set {realm:
  'cyclone'}`, `notes: true`.
- **Equations:** `sb` "Solid body and its circulation" ref 'Eqs. (3.22), (3.24)' `u_r=0,\ u_\theta=\omega_0r;\quad\Gamma=2\pi ru_\theta=2\pi r^2\omega_0`
  live · `wz` "Vorticity in polar coordinates" ref 'Eq. (3.23)'
  `\omega_z=\frac1r\frac{\partial}{\partial r}(ru_\theta)-\frac1r\frac{\partial u_r}{\partial\theta}` live "0.736 s⁻¹" · `lv` "Line vortex"
  ref 'Eqs. (3.25)–(3.27)' `u_\theta=B/r;\quad\Gamma=2\pi B;\quad[\omega_z]_{r\to0}=\lim_{r\to0}2B/r^2` (aligned) · `rk` "Rankine" ref
  'Eq. (3.28)' `u_\theta=\frac{\Gamma r}{2\pi\sigma^2}\ (r\le\sigma),\quad\frac{\Gamma}{2\pi r}\ (r>\sigma)` · `ga` "Gaussian" ref
  'Eq. (3.29)' `\omega_z=\frac{\Gamma}{\pi\sigma^2}e^{-r^2/\sigma^2},\quad u_\theta=\frac{\Gamma}{2\pi r}\big(1-e^{-r^2/\sigma^2}\big)`
  live "0.632 m/s at r = 1 m" note "peak where 1 + 2r²/σ² = e^{r²/σ²}".
- **Check yourself:** (1) "In the line vortex, a wheel at r = 0.5σ orbits fast. How fast does it spin?" — "Not at all:
  ω_z = 0 there, so ½ω_z = 0; it keeps its heading." `set {kind: 'line'}` · (2) "Solid body, loop of radius 1 m centred
  at (1.5, 0) m: Γ?" — "2ω₀ × π(1)² = 2π m²/s for ω₀ = 1 s⁻¹ — the same as a centred loop of the same area (N37)."
  `set {kind: 'solid', R: 1}` · (3) "Halve σ at fixed Γ in the Gaussian mode. What happens to the peak speed and its
  radius?" — "u_max doubles (∝ Γ/σ); its radius halves (1.1209σ)." `set {kind: 'gaussian', sigma: 0.5}` · (4) "Why does
  the Rankine speed have a kink at σ but the Gaussian does not?" — "Rankine's vorticity jumps from Γ/πσ² to 0 at σ;
  u_θ = Γ(r)/2πr then has a kink. The Gaussian's vorticity is smooth." `set {kind: 'rankine'}`.
- **Selftest parity rows:** `{name: 'r* / σ', js: rMax(1.0), py: 'ch03.gaussian_vortex_max_radius(1.0)', rtol: 1e-10}` ·
  `{name: 'Gaussian u_θ(1)', js: profile('gaussian', 1, 2*Math.PI, 1)[0], py: 'ch03.gaussian_vortex(1.0,
  6.283185307179586, 1.0)[0]', rtol: 1e-12}` · `{name: 'Rankine ω inside', js: profile('rankine', 0.5, 2*Math.PI, 1)[1],
  py: 'ch03.rankine_vortex(0.5, 6.283185307179586, 1.0)[1]', rtol: 1e-12}` · `{name: 'Γ off-centre (solid)', js:
  circ('solid', 1.0, [1.5, 0], 2*Math.PI, 1), py: 'ch03.circulation_circle("solid", 1.0, center=(1.5, 0.0),
  Gamma=6.283185307179586, sigma=1.0)', rtol: 1e-6}` · `{name: 'sector (line)', js: sector('line', 1.0, 0.1, 0.2),
  py: 'ch03.annular_sector_circulation("line", 1.0, 0.1, 0.2, Gamma=6.283185307179586, sigma=1.0)', atol: 1e-10}` ·
  `{name: 'mean vorticity (line, r=0.1)', js: meanVort('line', 0.1), py: 'ch03.mean_vorticity_in_disc("line", 0.1,
  Gamma=6.283185307179586, sigma=1.0)', rtol: 1e-8}` · `{name: 'polar ω (Gaussian, r=0.8)', js: polarVort('gaussian', 0.8),
  py: 'ch03.vortex_profile("gaussian", 0.8, 6.283185307179586, 1.0)[1]', rtol: 1e-6}` · invariant `{name: 'max equation',
  js: 1 + 2*rMax(1)**2 - Math.exp(rMax(1)**2), expect: 0, atol: 1e-12}`.
- **Fit plan:** portrait: `plane` + `profiles`; `circ` hidden (Γ of the loop in a readout, the peak in the `profiles`
  title); Explore: kind chips, Gamma, sigma, R (shape, realm optional); landscape phone one row; notebook three views.

### E7 · reynolds_transport_cv
- **Title:** "What changes inside a moving box?" · **Summary:** "A moving boundary sweeps a band coloured by the sign of
  b·n; the volume and surface terms of the Reynolds transport theorem add up to the measured rate of change of ∫F, and
  shrinking Δt shows the dropped second-order term vanish." · **CORE:** C15 (also N45–N55) · **Reference:**
  `fid_formula_lab.html` (terms that add up to a total, stage chips, code synced to the active line) with
  `forced_damped_vibrations.html`'s explanation panel.
- **meta:** `viz:order 7` · `viz:sections 3.6` · `viz:equations 3.30 3.31 3.32 3.33 3.34 3.35` · `viz:fluidpy
  ch03.leibniz_terms ch03.rtt_ellipse_2d ch03.rtt_field ch03.example_3_2 ch03.reynolds_transport
  ch03.volume_integral_rate_fd` · `viz:derivations D21 D22 D23`.
- **Physics (JS ↔ Python):** `F(name, x1, t)`, `dFdt(name, x1, t)` ↔ `ch03.rtt_field(name, dim)` (uniform / ramp /
  warming / carried — the fields depend on x₁ and t only) · `leibniz(p)` → [interior, upper, lower, total] (8-point
  Gauss–Legendre) ↔ `ch03.leibniz_terms` · `ellipseRTT(p)` → [volume, surface, total] (polar-grid volume sum, 256-point
  boundary sum of F b·n ds) ↔ `ch03.rtt_ellipse_2d` · `fdRate(p, dt)` (central difference of the area integral) ↔ the
  same idea as `ch03.volume_integral_rate_fd` · `fourTerms(p, dt)` → T1…T4 of (3.32) for the current shape (1-D exact, 2-D
  quadrature over the swept band) · `cone(h, r0, rdot)` → {direct, rtt} ↔ `ch03.example_3_2`.
- **Modes** (`mode` chips): **Leibniz 1-D** (interval [a(t), b(t)] over F(x, t)) · **RTT 2-D** (an ellipse with semi-axes
  a, b growing at ȧ, ḃ and centre moving at ċ, over a heat map of F) · **Ex. 3.2 cone** (side view of a cone of height
  h whose base radius grows at ṙ; F = 1).
- **Views** (rows [1.25, 0.8]):
  1. `cv` "The control volume" — Leibniz: F(x, t) (black) and F(x, t + Δt) (dashed), the area under F between a and b
     (light), the three strips of Fig. 3.17 (blue interior band ∂F/∂t·Δt, orange strip at b, rose strip at a), arrows ȧ,
     ḃ. RTT 2-D (`equal: true`): heat map of F, the ellipse at t (solid) and t + Δt (dashed), b arrows on 24 boundary
     points, the swept band shaded **blue where b·n > 0, rose where b·n < 0**, n drawn at one clickable boundary point.
     Cone: side section of the cone (apex at the origin), b arrows on the sloping side (length ∝ z/h), n arrows, the base
     disc, the swept sliver. Pointer: click a boundary point (inspector).
  2. `bars` "The budget" — waterfall: volume term ∫∂F/∂t (blue) + surface term ∮F b·n (orange; split into advancing
     (solid) and retreating (hatched, negative) parts) = total (purple), next to the finite-difference d/dt∫F (grey
     diamond); title "total 1.900 = measured 1.900".
  3. `limit` "Why the dropped term vanishes" (`hidePortrait: true`) — log–log |T4| = |∫_{ΔV}Δt ∂F/∂t| and the sliver
     error |∫_{ΔV}F − ∮FbΔt·n| vs Δt ∈ [10⁻⁴, 10⁻¹] s with slope-2 guides and the current Δt marker.
- **Controls:** `F` chips uniform / ramp / warming / carried (help: "F = 1; 1 + 0.5x; 1 + 0.5x + 0.3t; 1 + 0.5(x − 0.3t)")
  · `adot` "Lower end / x-axis rate $\dot a$" −0.5…0.5, default −0.2 m/s · `bdot` "Upper end / y-axis rate $\dot b$" −0.5…0.5,
  default 0.4 m/s · `cdot` "Translation $\dot c$" −0.5…0.5, default 0 m/s (*optional*) · `dt` "$\Delta t$" 10⁻⁴…0.3 s (log),
  default 0.05 · cone mode: `rdot` "$\dot r$" 0…0.3, default 0.1 m/s, `r0` 0.2…1 m, `h` 0.5…2 m (*optional*). Transport
  `t` 0 → 5 s, rate 0.5, `end: 'hold'`, card "Over 5 s: ∫F changed by …; ∫(volume term) dt = …, ∫(surface term) dt = …
  — they add up."
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **terms** (waterfall bars) · **transport** ·
  **modes** (1-D / 2-D / cone) · **presets** "fixed volume (b = 0)" {adot: 0, bdot: 0, cdot: 0, F: 'warming'} · "riding
  with a moving pattern" {mode: '2d', F: 'carried', adot: 0, bdot: 0, cdot: 0.3} (the two terms cancel) · "growing
  balloon" {mode: '2d', F: 'warming', adot: 0.2, bdot: 0.2, cdot: 0} · "one side retreats" {adot: 0.2, bdot: -0.2} · "Ex.
  3.2 cone" {mode: 'cone', h: 1, r0: 0.5, rdot: 0.1} · **status** ("📦 fixed volume: d/dt passes inside the integral" /
  "🎈 moving boundary: the surface term is …% of the total" / "↔ carried along: the two terms cancel") · **inspector**
  (click a boundary point: "b = (…, …), n = (…, …), b·n = … m/s → this patch sweeps (b·n)Δt ds = … m² per step, signed;
  F there = …; contribution F b·n ds = …").
- **Readouts:** "volume term" · "surface term" · "total" · "measured d/dt∫F" · "T4 (dropped)".
- **Explain:**
  0. *What the windows show* — "The **control volume** is the region whose content we track; its boundary moves with
     velocity b. The band between the solid (t) and dashed (t + Δt) boundaries is what the walls sweep:
     <b class="c-blue">blue</b> where they advance (b·n > 0), <b class="c-rose">rose</b> where they retreat. The **budget**
     adds the <b class="c-blue">volume term</b> and the <b class="c-orange">surface term</b> to the purple total, next to
     the measured rate (grey). **Why the dropped term vanishes** (hidden on phones) shows the second-order pieces."
  1. *The region now* (Leibniz) — "a = **1.000** m moving at ȧ = −0.2 m/s, b = **3.000** m moving at ḃ = +0.4 m/s; F =
     1 + 0.5x + 0.3t (warming)." (2-D: "semi-axes a = 2.000, b = 1.000 m, area πab = 6.283 m², ∇·b = ȧ/a + ḃ/b = …".)
  2. *The volume term* — "$\int_a^b\frac{\partial F}{\partial t}dx=0.3\times(3-1)=$ **0.600**: what changes in place."
     (boxed)
  3. *The surface term, piece by piece* — "Upper end: $\dot bF(b)=0.4\times2.5=$ **+1.000** (advancing, gains). Lower
     end: $-\dot aF(a)=-(-0.2)\times1.5=$ **+0.300** (the lower end moves *left*, which also gains). In 2-D: advancing arcs
     **+…**, retreating arcs **−…** (signed by b·n)." (boxed)
  4. *Total vs measured* — "(3.30): 0.600 + 1.000 + 0.300 = **1.900**; finite difference of ∫F over ±Δt: **1.900** —
     they agree to **2 × 10⁻⁹**." (boxed)
  5. *The error at this Δt* — "T4 = ∫_{ΔV}Δt ∂F/∂t = 0.05 × 0.3 × 0.03 = **4.5 × 10⁻⁴** at Δt = 0.05 s (∝ Δt²); dividing by Δt leaves
     **0.009**, which → 0 as Δt → 0 — why (3.33) drops it."
  6. *Ex. 3.2* (cone mode, or hint) — "$V=\tfrac13\pi hr^2$: $dV/dt=\tfrac23\pi hr_o\dot r=\tfrac23\pi\times1\times0.5\times0.1=$
     **0.1047** m³/s; RTT over the sloping side (b·n = 0 on the base): **0.1047** m³/s."
  7. *F = 1 and b = u (D23)* (2-D with F uniform, or hint) — "Then the total is dA/dt = ∮b·n ds = ∫∇·u dA with the
     linear flow u = (ȧ/a x, ḃ/b y): ∇·u = **0.2** s⁻¹ × area **6.283** = **1.257** m²/s — the (3.14) statement for a
     finite blob."
  8. *At the current time* — "t = `live t` s; ∫F = `live I`; region size `live L`."
  9. *Reading the current setting* — b = 0: "A fixed control volume: only the volume term is left, and d/dt may pass
     inside the integral." · carried: "The box rides along with a pattern: the volume term (F changes in place) and the
     surface term (the box moves into higher F) cancel — nothing inside changes." · growing: "The walls sweep new F in;
     the surface term dominates when F is large at the boundary." · one side retreating: "The rose band is negative
     volume: the same signed formula loses F there — no separate 'outflow' term needed." · cone: "Only the sloping side
     moves outward; the base slides within its own plane (b·n = 0)."
- **Derivation tab:**
  - **D21** (7 steps) `view: 'cv'`; goal `set {mode: '1d', F: 'warming', adot: -0.2, bdot: 0.4, t: 0}`; step 3 **live**
    "d/dt Φ(b(t), t) = F(b)ḃ + ∂Φ/∂t|_b = 2.5 × 0.4 + …"; step 6 shades the interior band; step 7 **live** "0.600 + 1.000 −
    (−0.300) = 1.900".
  - **D22** (12 steps) `view: 'cv'`; goal `set {mode: '2d', F: 'warming', adot: 0.2, bdot: -0.1, cdot: 0.3, dt: 0.1}`; step
    1 colours ΔV (blue/rose); step 3 shows the four (3.32) terms as four small bars (live values T1…T4); step 5 `set {dt:
    0.01}` `watch: "T4 falls 100× when Δt falls 10×"` (phones: the value in the `cv` title); step 8 draws (bΔt)·n on one
    boundary patch; step 9 **live** "∫_{ΔV}F = … vs ∮F(bΔt·n) = … (gap …)"; step 11 **live** "volume … + surface … = …";
    step 12 `set {mode: '1d'}` `watch: "two end points, n = ±1: Leibniz again"`.
  - **D23** (5 steps) `view: 'bars'`; goal `set {mode: '2d', F: 'uniform', adot: 0.2, bdot: 0.1, cdot: 0}`; step 2 **live**
    "∮u·n ds = … = ∫∇·u dA = (0.1 + 0.1) × 6.283"; step 5 `watch: "the bars reduce to the area growth rate"`.
- **Code:**
  ```python
  F, dF = ch03.rtt_field("{{F}}", dim={{dim}})          # the field and its time derivative
  # 1-D (Leibniz, (3.30)): interior, gain at b, loss at a, total
  print(ch03.leibniz_terms(F, dF, {{a}}, {{b}}, {{adot}}, {{bdot}}, {{t}}))   # ({{i}}, {{up}}, {{lo}}, {{tot}})
  # 2-D (RTT, (3.35)): an ellipse growing and translating
  vol, surf, tot = ch03.rtt_ellipse_2d(F, dF, {{ea}}, {{eb}}, {{adot}}, {{bdot}},
                                       (0.0, 0.0), ({{cdot}}, 0.0), {{t}})
  print(vol, surf, tot)                                   # {{vol}} + {{surf}} = {{tot2}}
  print(ch03.example_3_2({{h}}, {{r0}}, {{rdot}}))        # cone: direct = RTT = {{cone}} m³/s
  ```
- **Walkthrough (8 steps):** 1. "A box that moves" — "How fast does the amount of F inside change when the walls move?
  Two reasons: F changes in place, and the walls sweep F in or out." `set {mode: '1d', F: 'warming', t: 0}` · 2. "In
  place" — "The blue band: ∂F/∂t over the interval, 0.6 per second." `terms: true`, highlight `view:cv` · 3. "Swept in" —
  "The orange strip: the upper end advances at ḃ over F(b): +1.0. The lower end moves left: −ȧF(a) = +0.3." `eq: 'leib'`,
  `derive: {id: 'D21', step: 7}` · 4. "In 2-D and 3-D" — "Now a blob: blue where the wall advances, rose where it
  retreats. One signed surface integral replaces the two end terms." `set {mode: '2d', adot: 0.2, bdot: -0.1, cdot: 0.3}`,
  `eq: 'rtt'` · 5. "Check it" — "The grey diamond is d/dt∫F measured by finite differences: the bars land on it." `readouts:
  ['tot', 'fd']` · 6. "Second order" — "Shrink Δt: the dropped term T4 falls like Δt²." `set {dt: 0.01}`, `derive: {id:
  'D22', step: 5}` · 7. "Riding a pattern" — "Box and pattern move together: the two terms cancel exactly." `set {F:
  'carried', adot: 0, bdot: 0, cdot: 0.3}` · 8. "Your turn" — "Ex. 3.2: predict dV/dt for h = 1 m, r₀ = 0.5 m, ṙ = 0.1
  m/s from the formula, then open the cone." `set {mode: 'cone'}`, `controls: ['rdot', 'r0']`.
- **Equations:** `leib` "Leibniz's theorem" ref 'Eq. (3.30)'
  `\frac{d}{dt}\int_{a(t)}^{b(t)}F\,dx=\int_a^b\frac{\partial F}{\partial t}dx+\frac{db}{dt}F(b,t)-\frac{da}{dt}F(a,t)` live "0.600 + 1.000
  + 0.300" · `def` "Definition" ref 'Eq. (3.31)'
  `\frac{d}{dt}\int_{V^*}F\,dV=\lim_{\Delta t\to0}\frac1{\Delta t}\Big\{\int_{V^*(t+\Delta t)}F(t+\Delta t)dV-\int_{V^*(t)}F(t)dV\Big\}`
  (aligned over two rows) · `sliver` "Swept sliver" ref 'Eq. (3.34)' `\int_{\Delta V}F\,dV\cong\int_{A^*}F\,(\mathbf b\Delta t\cdot\mathbf n)\,dA`
  live "… vs …" · `rtt` "Reynolds transport theorem" ref 'Eq. (3.35)'
  `\frac{d}{dt}\int_{V^*}F\,dV=\int_{V^*}\frac{\partial F}{\partial t}dV+\int_{A^*}F\,\mathbf b\cdot\mathbf n\,dA` live "0.000 + 2.199 =
  2.199" · `cone` "Ex. 3.2" ref 'Ex. 3.2' `\frac{dV}{dt}=\int_{A^*}\mathbf b\cdot\mathbf n\,dA=\tfrac23\pi hr_o\dot r` live "0.1047 m³/s".
- **Check yourself:** (1) "With b = 0 (fixed volume) and F = warming, what are the two terms?" — "Volume term 0.3 ×
  length (or area); surface term 0: d/dt passes inside the integral." `set {adot: 0, bdot: 0, cdot: 0, F: 'warming'}` ·
  (2) "Make one wall retreat. What sign does its sliver carry, and do you need a separate formula?" — "Negative (b·n < 0,
  rose); no — the single signed surface integral handles it." `set {mode: '2d', adot: 0.2, bdot: -0.2}` · (3) "F uniform
  and b = u of a linear flow: what is the total?" — "The area growth rate ∮u·n = ∫∇·u dA — (3.14) for a finite blob
  (D23)." `set {mode: '2d', F: 'uniform'}` · (4) "In the cone, why does the base contribute nothing although it grows?"
  — "Its points slide within the base plane: b·n = 0 there (the book writes b = 0, which is not true)." `set {mode:
  'cone'}`.
- **Selftest parity rows:** `{name: 'Leibniz total', js: leibniz({F: 'warming', a: 1, b: 3, adot: -0.2, bdot: 0.4, t:
  0}).total, py: 'ch03.leibniz_terms(ch03.rtt_field("warming", 1)[0], ch03.rtt_field("warming", 1)[1], 1.0, 3.0, -0.2,
  0.4, 0.0)[3]', rtol: 1e-10}` · `{name: 'Leibniz lower', js: …lower, py: '…[2]', rtol: 1e-10}` · `{name: 'ellipse RTT
  total', js: ellipseRTT({F: 'ramp', a: 2, b: 1, adot: 0.2, bdot: 0.1, cdot: 0.3, t: 0}).total, py:
  'ch03.rtt_ellipse_2d(ch03.rtt_field("ramp", 2)[0], ch03.rtt_field("ramp", 2)[1], 2.0, 1.0, 0.2, 0.1, (0.0, 0.0), (0.3,
  0.0), 0.0)[2]', rtol: 1e-6}` · `{name: 'carried cancels', js: ellipseRTT({F: 'carried', a: 2, b: 1, adot: 0, bdot: 0,
  cdot: 0.3, t: 1}).total, py: 'ch03.rtt_ellipse_2d(ch03.rtt_field("carried", 2)[0], ch03.rtt_field("carried", 2)[1],
  2.0, 1.0, 0.0, 0.0, (0.0, 0.0), (0.3, 0.0), 1.0)[2]', atol: 1e-8}` · `{name: 'cone', js: cone(1, 0.5, 0.1).rtt, py:
  'ch03.example_3_2(1.0, 0.5, 0.1)["rtt_dblquad"]', rtol: 1e-8}` · `{name: 'F warming', js: F('warming', 2.0, 1.0), py:
  'ch03.rtt_field("warming", 1)[0](2.0, 1.0)', rtol: 1e-12}` · invariant `{name: 'FD = RTT (2-D)', js: fdRate(p, 1e-4) -
  ellipseRTT(p).total, expect: 0, atol: 1e-6}`.
- **Fit plan:** portrait: `cv` + `bars`; `limit` hidden (T4 in a readout and the `bars` title); Explore: mode chips, F
  chips, adot, bdot, dt (cdot, cone sizes optional); landscape phone one row; notebook three views.

### B1 · material_volume_divergence
- **Title:** "How fast does a blob of fluid swell?" · **Summary:** "A blob of tracers in a compressible flow grows or
  shrinks; its area follows e^{t tr G}, the surface integral ∮u·n equals ∫∇·u, and as the blob shrinks the fractional
  rate tends to the local divergence." · **CORE:** C09, C15 (also N53) · **Reference:** `overfitting_curves.html` (a
  minimal two-slider figure with a verdict line).
- **meta:** `viz:sections 3.4 3.6` · `viz:equations 3.14 3.35` · `viz:fluidpy ch03.material_volume_ratio
  ch03.volumetric_strain_rate ch03.material_volume_rate` · `viz:derivations D11 D23`.
- **Physics:** `areaRatio(G, t)` ↔ `ch03.material_volume_ratio` · nonlinear preset u = (x + 0.5x², y) (∇·u = 2 + x) with
  the blob tracked by RK4 · `fluxBlob(u, blob)` (∮u·n ds) and `divInt(u, blob)` ↔ `ch03.material_volume_rate` (2-D
  analogue, stated) .
- **Views:** `blob` (tracers, the blob outline at t and its ghost e^{t tr G} scaling, u·n colouring of the outline) ·
  `graph` (area ratio vs t: tracked black, ghost e^{t tr G} blue dashed, first order grey; for the nonlinear preset the
  fractional rate (1/A)dA/dt vs blob radius with ∇·u at the centre as the limit line).
- **Controls:** `tr` "tr G" −1…1 default 0.4 s⁻¹ (with a fixed traceless part), `radius` 0.05…1 m, `field` chips
  linear / nonlinear, transport t 0 → 3 s.
- **Depth features:** Explain, Code, Derivation + **linked views** (2) · **transport** · **status** ("🎈 swelling at
  tr G = 0.4 s⁻¹" / "🫧 shrinking" / "💧 area-preserving").
- **Explain:** 0 views; 1 tr G = ∇·u from the setting; 2 fractional rate from the tracked outline vs (3.14),
  $\frac1{\delta V}\frac{D(\delta V)}{Dt}=\frac{\partial u_i}{\partial x_i}$; 3 ∮u·n ds = ∫∇·u dA (Gauss) — the RTT with F = 1,
  b = u; 4 the small-blob limit (nonlinear field); 5 reading: "divergence is the growth rate of a blob per unit of its own
  size — the kinematic half of mass conservation (Ch. 4)".
- **Derivation tab:** D11 (7 steps, `view: 'graph'`, step 7 live "e^{t tr G} = …"), D23 (5 steps, `view: 'blob'`).
- **Code:** `ch03.volumetric_strain_rate(G)` · `ch03.material_volume_ratio(G, t)` · `ch03.material_volume_rate(u, cv, t)`
  with live values.
- **Walkthrough (5 steps):** a blob in still fluid → swelling flow → the area curve vs its ghost → the surface integral
  equals the divergence integral → your turn (predict the area after 2 s at tr G = −0.5).
- **Equations:** (3.14) $\frac1{\delta V}\frac{D}{Dt}(\delta V)=\frac{\partial u_i}{\partial x_i}=S_{ii}$ · (3.35) with F = 1,
  $\frac{dV^*}{dt}=\int_{A^*}\mathbf b\cdot\mathbf n\,dA$ · Jacobi $\det e^{\mathbf Gt}=e^{t\,\mathrm{tr}\,\mathbf G}$.
- **Check yourself:** (1) "tr G = 0.4 s⁻¹: area after 2 s?" — "e^{0.8} = 2.23." (2) "Why does the nonlinear blob's rate
  depend on its size?" — "∇·u = 2 + x varies across it; the average over a small blob → the centre value." (3) "tr G = 0
  but the blob changes shape — is that allowed?" — "Yes: deformation at constant area."
- **Selftest parity rows:** `{name: 'area ratio', js: areaRatio([[0.5, 0.2], [0, 0.3]], 2), py:
  'ch03.material_volume_ratio([[0.5, 0.2], [0.0, 0.3]], 2.0)', rtol: 1e-10}` (= e^{1.6} = 4.953) · `{name: 'tr G', js:
  trG([[0.5, 0.2], [0, 0.3]]), py: 'ch03.volumetric_strain_rate([[0.5, 0.2], [0.0, 0.3]])', rtol: 1e-12}`.
- **Fit plan:** two views side by side (landscape) or stacked (portrait); ≤ 4 controls.

---

## Part D — runtime budget (full run < 5 min on a laptop / Colab CPU)

This chapter is ODE- and quadrature-heavy but small: path lines (DOP853, rtol 1e-10) take ~5–20 ms each, a vectorised
streak line of 180 release times ~0.1 s, a `quad`/`dblquad` a few ms, a stencil evaluation µs. The heavy items are the
animations (≈ 0.12 s per frame at dpi 80, from the ch01/ch02 timings), the plotly slider figures that precompute
streamlines per step, and the sympy cells (D02, D06, D17, D22 checks, `lagrangian_to_eulerian`).

| Section | Heaviest cells | Full | FAST (`FLUIDPY_FAST=1`) |
|---|---|---|---|
| setup + imports | numpy/scipy/sympy/plotly/pint | 8 s | 8 s |
| §3.1 | `quad` average, profile figure, plotly 3-D unit vectors (9 cones) | 2 s | 2 s |
| §3.2 C01, C02 | sympy `lagrangian_to_eulerian` (≈ 1 s), chain-rule primer, D02 sympy line, front figure (2 × 200 points), time slider 25 steps × 3 traces × 97 points, float path line | 6 s | 5 s (13 steps) |
| §3.3 C03 | 3 streamlines + solid body, 4-panel streamplot (21²), t′ slider 24 steps (closed forms), 3-D stream tube (12 streamlines × 400 points) | 5 s | 4 s (12 steps, 8 tube lines) |
| §3.3 C04 | 2 path lines, streak line (181 releases, one vector ODE), RK4 loop (200 steps), Fig. 3.7, animation 60 frames (dye 60 dots, closed form per frame), N13 steady check (3 curves) | 12 s | 7 s (30 frames) |
| §3.3 C05 | chain-rule sympy primer, D06 sympy check (optional, ≈ 2 s), cylinder two-frame streamplot (121 × 81), term-bar figure (5 observers), observer slider 21 steps × 14 streamlines (≈ 0.02 s each → 6 s, **cached** in one list computed once) | 12 s | 6 s (11 steps, 8 seeds) |
| §3.4 C06–C09 | stencil G, Taylor remainder (6 points), ring figure, rose (72 tracked threads via expm, µs each), corner-angle figure (3 × 200 times), cube tracking, volume figure | 6 s | 5 s |
| §3.4 C10–C12 | `velocity_potential_2d` on 61² (≈ 0.3 s), γ slider 21 steps × 3 traces × 181 points, ring split figure, circle→ellipse video 48 frames (72 tracers, SVD per frame) | 10 s | 6 s (41² grid, 24 frames, 11 steps) |
| §3.5 head + C13, C14 | Fig. 3.14 frames player 16 frames, sympy (3.23) ×2, circulations (512-point loops), sector `quad` ×4, log–log disc figure (60 radii), paddle-wheel video 60 frames (8 wheels + 2 elements, closed-form circular motion), vortex profiles figure, σ slider 15 steps × 3 traces × 300 points, `brentq` + `lambertw` | 16 s | 9 s (10 + 30 frames, 8 steps) |
| §3.6 C15 | Leibniz frames player 12 frames, D21 sympy, D22 sympy ★★★ (series + integrate on a sphere ≈ 3 s), RTT on `GrowingSphere` (Gauss–Legendre 16 × 16 × 32 ≈ 10 ms), midpoint from-scratch (16 × 16 × 32 ≈ 20 ms), Δt slider 16 steps (`swept_terms_sphere`, ≈ 50 ms each → 0.8 s, cached), log–log figure, `example_3_2` (`dblquad` ≈ 10 ms), Fig. 3.18 blob figure | 14 s | 9 s (8 frames, 10 Δt steps, 12 × 12 × 24 grids) |
| explainers (7 × `show_viz`) | read the HTML files | 1 s | 1 s |
| **Total** | | **≈ 92 s** | **≈ 62 s** |

**FAST plan.** Every size-dependent choice is written `a if not FAST else b` in the cell: animations 60/48 frames → 30/24,
frames players 16/12 → 10/8; slider figures ≤ 25 steps → ≤ 13, ≤ 4 traces × ≤ 400 points (each < 250 kB); streak lines
181 → 91 release times; stream tube 12 → 8 lines; potential grid 61² → 41²; RTT quadrature 16 × 16 × 32 → 12 × 12 × 24.
**Cached arrays:** `u31 = ch03.preset_field("ex31")` and `d = ch03.example_3_1(…)` are built once in C03 and reused by
C04; the Ex. 3.1 streak ensemble of C04's code cell is reused by the animation (one ODE solve, positions interpolated
with `dense_output`); the cylinder streamlines for the 21 observer speeds are computed once in C05's figure cell and fed to
the slider; `Gsh`, `Gsb`, `Grand` are defined once in C10/C09 and reused in C11, C12; the swept-term table for the Δt
slider is computed once in C15. Outputs: 3 videos (< 2 MB each), 3 frames players (< 3 MB each), 7 slider figures, 2
plotly 3-D figures — page well under 15 MB. Sympy cells never simplify expressions with more than two unknown
functions (D06 check uses `sp.Function` components and `expand` + `simplify` on each component separately, ≈ 2 s).

---

## Part E — prerequisite ledger
Every concept, symbol, maths tool and Python function or idiom the notebook or its explainers use, with where it is
explained. "primer (in Cxx)" = a 📎 primer placed in that block before first use (the primer term is the Concept text
before any parenthesis — the builder uses it verbatim in `nb.primer`); "primer (in §3.1, before C01)" = a primer in the
§3.1 section, which has no A item; "gloss" = one sentence where used; "knowledge/primers.md: <term> (chNN)" = a one-line
reminder naming the earlier primer. IDs are this chapter's curation IDs (A items explain themselves; B/C items are
explained by their note, named in parentheses). Ch. 2 material that is neither a primer nor a ch03 RECAP is named by
section ("Ch. 2 §2.9"), never by a Ch. 2 C-number.

| Concept | First used in | Explained by |
|---|---|---|
| kinematics vs dynamics | §3.1 | C01 (N01 note in the §3.1 section) |
| steady vs unsteady flow, ∂(·)/∂t = 0 | §3.1 | C01 (N02 note); frame dependence shown in C05 |
| 1-D, 2-D, 3-D flow; cross-section average | §3.1 | C01 (N03 note) |
| scipy.integrate.quad and dblquad | §3.1 (N03 code) | primer (in §3.1, before C01) |
| ring area element 2πr dr in a pipe section | §3.1 (N03 code) | C01 (code comment in the §3.1 cell; the polar area element is primed in C13) |
| np.trapezoid (trapezoid rule) | §3.1 | knowledge/primers.md: trapezoid rule (ch01 P37) — reminder |
| cylindrical and spherical unit vectors | §3.1 (N05) | primer (in §3.1, before C01) |
| coordinate systems: plane polar, cylindrical (R, φ, z), spherical (r, θ, φ) | §3.1 | C01 (N05 note, with the ⚠️ on θ) |
| velocity components as projections on a local basis | §3.1 | knowledge/primers.md: orthonormal basis, projection and completeness (ch02 P65) — reminder |
| plotly 3-D arrows and lines | §3.1, C03 (N15) | knowledge/primers.md: plotly 3-D arrows, lines and meshes (ch02 P64) — reminder |
| Appendix B curvilinear operators | §3.1 | C01 (N06 note, pointer) |
| matplotlib figures | §3.1 | knowledge/primers.md: matplotlib figures (ch01 P01) — reminder |
| f-strings | §3.1 | knowledge/primers.md: f-strings (ch01 P04) — reminder |
| Lagrangian description; a label is not a variable | C01 | C01 (N07 note) |
| functions of time with parameters | C01 | primer (in C01) |
| particle velocity and acceleration (3.1) | C01 | C01 (N08 note) |
| ordinary and second derivative | C01 | knowledge/primers.md: ordinary derivative as a slope (ch01 P19) — reminder |
| Eulerian description F(x, t); compatibility (3.2) | C01 | C01 |
| inverse functions and sympy solve | C01 | primer (in C01) |
| sympy (symbols, diff, simplify, Function, series) | C01, C02, C05, C13, C15 | knowledge/primers.md: sympy (ch01 P40) — reminder |
| exponent rules | C01 (D01 step 5) | knowledge/primers.md: exponent rules (ch01 P43) — reminder |
| natural exponential e^{αt} and ∫e^{−s}ds | C01, C14 (D19) | knowledge/primers.md: natural logarithm and exponential (ch01 P36) — reminder |
| central finite differences | C01 check, C02 check | knowledge/primers.md: finite differences (ch01 P21) — reminder |
| functions as arguments and lambda | C01 | knowledge/primers.md: functions as arguments and lambda (ch01 P29) — reminder |
| assert np.allclose | C01 | knowledge/primers.md: assert np.allclose (ch01 P15) — reminder |
| tuple unpacking | C01 | knowledge/primers.md: tuple unpacking (ch01 P14) — reminder |
| continuum: every point hosts a particle | C02 (D02 step 6) | C02 (gloss in D02 step 6's *why*; continuum from Ch. 1 §1.4) |
| multivariable chain rule along a path | C02 (D02) | primer (in C02) |
| chain rule in two variables | C02 | knowledge/primers.md: chain rule (ch01 P49) — reminder inside the primer |
| partial derivative | C02 | knowledge/primers.md: partial derivative (ch01 P25) — reminder |
| gradient ∇F and the dot product u·∇F | C02 (D02 step 8) | C02 (D02 step 8's *why*; gradient from Ch. 2 §2.9, dot product from Ch. 2 §2.1) |
| material derivative D/Dt (3.3)–(3.5) | C02 | C02 (D02; N09, N11 notes) |
| local vs advective rate; advection vs convection | C02 | C02 (N10 note) |
| summation convention in (3.5), expand_indices_str | C02 | C02 (N11 note and code comment; convention from Ch. 2 §2.1) |
| directional derivative ∂F/∂s and the streamwise form (3.6) | C02 | knowledge/primers.md: level sets and the directional derivative (ch02 P75) — reminder in the N12 note |
| unit vector e_u = u/‖u‖ | C02 | C02 (N12 note) |
| warm and cold advection; K/s → K/h | C02 | C02 (worked example) |
| tanh front and sech² = 1 − tanh² | C02 figure, E2 | C02 (code comment in the figure cell) |
| scipy.integrate.solve_ivp (defaults) | C02 check (path line of the float) | knowledge/primers.md: scipy.integrate.solve_ivp (ch01 P31) — reminder |
| Python dictionaries | C02 (terms dict) | knowledge/primers.md: Python dictionaries (ch01 P23) — reminder |
| slider_figure | C02 | knowledge/primers.md: slider_figure (ch01 P17) — reminder |
| show_viz | C02 | knowledge/primers.md: show_viz (ch01 P18) — reminder |
| streamline; frozen time | C03 | C03 (D03) |
| parametric curves, tangent vector and arc length | C03 (D03) | primer (in C03) |
| parallel vectors and the cross-product test | C03 (D03) | primer (in C03) |
| cross product in components | C03 (D03 step 3's *why*), C11 | C03 (taught with parallel vectors before D03; cross product from Ch. 2 §2.7) |
| solve_ivp options: t_eval, dense_output, events, backward integration | C03 | primer (in C03) |
| stagnation point (‖u‖ = 0 ends a streamline) | C03 | C03 (the solve_ivp options text and D03's traps) |
| np.meshgrid and plt.streamplot | C03 figure | knowledge/primers.md: np.meshgrid and the project grid layout (ch02 P76); plt.contour, plt.quiver and plt.streamplot (ch02 P78) — reminders |
| stream tube; flux through a section ∫u·n dA | C03 | C03 (N15 note; flux from Ch. 2 §2.12) |
| path line (3.8) as an initial-value problem | C04 | C04 |
| streak line; release time as the parameter | C04 | C04 (N17 note) |
| Ex. 3.1 symbols t′, t_o, c_x, c_y, ξ_o | C04 | C04 (N18 note ⚠️) |
| integrating sin ωt and cos ωt | C04 (D04 step 2) | C04 (gloss in D04 step 2's *why*) |
| eliminating a parameter with sin² + cos² = 1 | C04 (D04, D05) | C04 (gloss in D04 steps 6–7) |
| equation of a circle (x − a)² + (y − b)² = r² | C04 | C04 (gloss in D04 step 7's *plain*) |
| implicit differentiation | C04 (D05 step 7) | C04 (gloss in D05 step 7's *why*) |
| RK4 by hand | C04 | primer (in C04) |
| explicit stepping (Euler) | C04 primer | knowledge/primers.md: explicit stepping (ch01 P30) — reminder |
| animate and show_animation | C04 | knowledge/primers.md: animate and show_animation (ch01 P16) — reminder |
| steady flow ⇒ the three lines coincide | C04 | C04 (N13 note) |
| PIV (particle image velocimetry) | C04 | C04 (gloss in the (3.8) markdown) |
| frames of reference and relative velocity | C05 | primer (in C05) |
| Galilean transformation u = U + u′, x = x′ + Ut + x′_o | C05 | C05 (N21 note) |
| chain rule with a moving frame | C05 (D06) | primer (in C05) |
| ideal flow past a cylinder (uniform stream + doublet) | C05 | C05 (N04 note; derived in Ch. 6) |
| acceleration of a particle Du/Dt = ∂u/∂t + (u·∇)u | C05 | C02 (D02 step 9) |
| Kronecker delta ∂x_j/∂x_k = δ_jk | C05 (D06 step 2), C08 (D10) | C05 (D06 step 2's *why*; δ from Ch. 2 §2.7) |
| linear vs quadratic terms (the nonlinearity) | C05 | C05 (N22 note) |
| frame-dependent split, invariant sum | C05 | C05 (N20, N23 notes) |
| traffic-light and roller-coaster analogies | C05 | C05 (N24 note) |
| Galilean vs rotating frames | C05 | C05 (⚠️ callout); rotating frame primed in C10 |
| velocity gradient G_ij = ∂u_i/∂x_j (row = component) | C06 | R01 (recap before C06) |
| split ∂u_i/∂x_j = S_ij + ½R_ij (3.11) | C06 | R01 |
| strain-rate tensor (3.12) | C06 | R02 |
| rotation tensor (3.13), book R = G − Gᵀ = 2A | C06 | R03 |
| multivariable first-order Taylor expansion | C06 (D07) | primer (in C06) |
| first-order Taylor expansion in one variable | C06 primer, C15 (D22 step 2) | knowledge/primers.md: first-order Taylor expansion (ch01 P26) — reminder |
| relative velocity (3.10) | C06 | C06 (D07) |
| matrix–vector product G·dx | C06 | knowledge/primers.md: matrix multiplication, transpose and identity (ch02 P63) — reminder |
| observed order of convergence, log–log slope | C06, C15 | knowledge/primers.md: power laws and log–log plots (ch01 P13) — reminder |
| limits and orders of smallness | C06 (D07), C07 (D08), C15 (D22) | knowledge/primers.md: limits and orders of smallness (ch02 P68) — reminder |
| deformation rates vs stress (why S matters) | C06 | C06 (N25, N26 notes) |
| material line element | C07 (D08) | primer (in C07) |
| scipy.linalg.expm and linear trajectories | C07 | knowledge/primers.md: scipy.linalg.expm (ch02 P79) — reminder |
| linear strain rate, n·S·n | C07 | C07 (D08) |
| n·A·n = 0 for antisymmetric A | C07 (D08 step 6) | C07 (gloss in D08 step 6's *why*); knowledge/primers.md: quadratic form and the Rayleigh quotient (ch02 P82) — reminder |
| Greek index = no summation | C07 | C07 (D08 step 5's *why*) |
| np.linalg.norm | C07 check | knowledge/primers.md: np.linalg.norm and np.linalg.qr (ch02 P67) — reminder |
| small-angle approximation | C08 (D09) | primer (in C08) |
| angles α (clockwise) and β (counterclockwise) of Fig. 3.11 | C08 | C08 (idea box and D09 step 1) |
| shear strain rate S₁₂ = ½D(α + β)/Dt; n₁·S·n₂ | C08 | C08 (D09) |
| rigid-body velocity Ω × x | C08 (D10) | primer (in C08) |
| ε_ijk, its antisymmetry and (a × b)_i = ε_ijk a_j b_k | C08 (D10), C11 (D15) | C08 (D10 steps 1 and 3's *why*; ε from Ch. 2 §2.7) |
| R_ij = −ε_ijk ω_k used before its recap | C08 (D10 step 4) | C08 (D10 step 4 writes (3.15) out with its Ch. 2 origin); R05 recaps it before C10 |
| rigid motion ⇒ S = 0; S independent of the observer | C08 | C08 (N27 note, D10) |
| volumetric strain rate (3.14) | C09 | C09 (D11) |
| product rule for three factors | C09 (D11 step 1) | C09 (gloss in D11 step 1); knowledge/primers.md: product rule for differentials (ch01 P38) — reminder |
| trace invariance under rotation | C09 (D11 step 6) | C09 (D11 step 6's *why*; invariants from Ch. 2 §2.5) |
| Jacobi's formula det e^{Gt} = e^{t tr G} | C09 | C09 (gloss in D11 step 7 + numeric check in the code cell) |
| np.linalg.det | C09 check | knowledge/primers.md: np.linalg.det and np.linalg.matrix_rank (ch01 P56) — reminder |
| np.random.default_rng | C09 code | knowledge/primers.md: np.random.default_rng (ch01 P10) — reminder |
| R antisymmetric ↔ a vector | C10 | R04 (recap before C10) |
| R_ij = −ε_ijk ω_k (3.15) | C10 | R05 |
| vorticity components (3.16) | C10 | R06 |
| circulation (3.18) and Stokes' theorem | C10, C13, C14 | R07 |
| angular velocity of a line | C10 (D12, D14) | primer (in C10) |
| element spin = ½ω; D12 | C10 | C10 |
| parallel shear flow u = (u₁(x₂), 0), γ = du₁/dx₂, ω₃ = −γ | C10 | C10 (markdown before D14); N34 restates it at the head of §3.5 |
| rotating frame of reference | C10 (D13) | primer (in C10) |
| curl of a rigid rotation = 2Ω | C10 (D13 step 3) | C08 (D10: ω = 2Ω) |
| vorticity depends on the frame, ω′ = ω − 2Ω | C10 | C10 (N28 note, D13) |
| planetary vorticity f = 2Ω sin φ; absolute vs relative vorticity | C10 | C10 (N28 note, gloss; developed in Ch. 13) |
| irrotational flow (3.17) and the velocity potential | C10 | C10 (N29 note) |
| simply connected region | C10 | C10 (gloss in the N29 note) |
| ∇×∇φ = 0 | C10 | C10 (N29 note; identity from Ch. 2 §2.13) |
| np.deg2rad | C10 code | knowledge/primers.md: cosines of angles between unit vectors (ch02 P66) — reminder |
| relative velocity = deformation + rigid rotation (3.19) | C11 | C11 (D15) |
| one index swap of ε flips the sign | C11 (D15 step 4) | C11 (D15 step 4's *why*; Ch. 2 §2.7) |
| numpy arrays with three axes (the ε array) | C11 check | knowledge/primers.md: numpy arrays with three axes and np.transpose (ch02 P73) — reminder |
| eigenvalues and eigenvectors of S; principal axes | C12 | knowledge/primers.md: eigenvalues and eigenvectors (ch02 P80) — reminder |
| principal frame (3.20)–(3.21) | C12 | C12 (N31, N32 notes) |
| linear map of a circle is an ellipse | C12 (D16) | primer (in C12) |
| singular values (np.linalg.svd) for finite time | C12 | C12 (gloss inside the linear-map-of-a-circle explanation) |
| equation of an ellipse/ellipsoid Σ x_α²/a_α² = 1 | C12 (D16 step 7) | C12 (gloss in D16 step 7's *why*) |
| tensor components in rotated axes | C12 (D16 step 2) | C12 (D16 step 2's *why*; rule (2.12) from Ch. 2 §2.4) |
| np.arctan2 | C12 code | knowledge/primers.md: np.arctan2 (ch02 P70) — reminder |
| shear-flow elements ABCD and PQRS (Fig. 3.14) | §3.5 head | C12 (N35 note) |
| polar coordinates as a moving basis | C13 (D17) | primer (in C13) |
| vorticity in polar coordinates (3.23) | C13 | C13 (D17) |
| solid-body rotation (3.22) and its circulation (3.24) | C13 | C13 (N36, N37 notes) |
| line vortex (3.25), Γ = 2πB (3.26) | C13 | C13 (N38, N39 notes) |
| delta-function core (3.27) | C13 | C13 (N40 note; delta gloss from Ch. 2 §2.12) |
| Γ_ABCD = 0 for a loop not enclosing the axis | C13 | C13 (N41 note) |
| line integral of a vector field around a loop | C13 (D17) | knowledge/primers.md: line integral of a vector field around a loop (ch02 P86) — reminder |
| Rankine vortex (3.28); piecewise functions with np.where | C14 | C14; knowledge/primers.md: np.where and np.select (ch01 P46) — reminder |
| continuity at a join; a kink as a maximum | C14 (D18) | C14 (gloss in D18 steps 3 and 5) |
| maximum of a function: derivative = 0, exclude trivial roots | C14 (D20) | C14 (gloss in D20 steps 3 and 5) |
| substitution in an integral | C14 (D19) | primer (in C14) |
| Gaussian vortex (3.29); Lamb–Oseen pointer | C14 | C14 (N43 note, D19) |
| np.expm1 and cancellation near zero | C14 | primer (in C14) |
| scipy.optimize.brentq | C14 (D20) | primer (in C14) |
| Lambert W function | C14 (D20 check) | C14 (gloss in D20's check) |
| real-vortex scales (bathtub, tornado, cyclone) | C14, E6 | C14 (N42 note and the code-cell table) |
| control volume V*, control surface A*, outward n, surface velocity b | C15 | C15 (N48 note) |
| differentiation under the integral sign | C15 (D21) | primer (in C15) |
| fundamental theorem of calculus; antiderivative | C15 (D21) | knowledge/primers.md: fundamental theorem of calculus (ch02 P84) — reminder |
| definite integral | C15 | knowledge/primers.md: definite integral (ch01 P27) — reminder |
| Leibniz's theorem (3.30) | C15 | C15 (N46 note, D21) |
| signed swept volume of a moving surface | C15 (D22) | primer (in C15) |
| Taylor expansion in time F(x, t + Δt) | C15 (D22 step 2) | C15 (gloss in D22 step 2's *why*) |
| mean-value theorem for integrals | C15 (D22, D23) | knowledge/primers.md: mean-value theorem for integrals (ch02 P85) — reminder |
| volume and surface integrals as midpoint sums | C15 check | knowledge/primers.md: volume and surface integrals as midpoint sums (ch02 P83) — reminder |
| Gauss' divergence theorem | C15 (D23, D24) | C15 (D23 step 2's *why*; Gauss (2.30) from Ch. 2 §2.12) |
| product rule ∇·(Fu) = u·∇F + F∇·u | C15 (D24 step 3) | C15 (gloss in D24 step 3's *why*) |
| Reynolds transport theorem (3.31)–(3.35) | C15 | C15 (D22; N49–N53 notes) |
| cone geometry: half-angle, slanted area element | C15 | C15 (N54 note) |
| sphere volume 4πR³/3 and area 4πR² | C15 | C15 (worked example, gloss) |
| continuity equation Dρ/Dt + ρ∇·u = 0 (preview) | C15 | C15 ("What would change if" cell; taught in Ch. 4) |
| Python objects with methods (GrowingSphere) | C15 | C15 (code comment: "a control volume object that knows its surface and its velocity b at any t") |
| RK4 by hand, in JavaScript (explainers) | E1, E3, E6 | primer (in C04) |
| heat maps with a colour scale | E2, E6, E7 | C02 (figure legends; `imshow` glossed in Ch. 2) |
| climate hooks: thermal advection, f + ζ, frontogenesis, cyclone profiles | C02, C10, C12, C14 | C02, C10, C12, C14 (the plain-words paragraphs; developed in Ch. 13) |

Ledger rows: 163. **Primers planned** (`nb.primer`, new in ch03, P87–P110): 24 — scipy.integrate.quad and dblquad ·
cylindrical and spherical unit vectors · functions of time with parameters · inverse functions and sympy solve ·
multivariable chain rule along a path · parametric curves, tangent vector and arc length · parallel vectors and the
cross-product test · solve_ivp options: t_eval, dense_output, events, backward integration · RK4 by hand · frames of
reference and relative velocity · chain rule with a moving frame · multivariable first-order Taylor expansion · material
line element · small-angle approximation · rigid-body velocity Ω × x · angular velocity of a line · rotating frame of
reference · linear map of a circle is an ellipse · polar coordinates as a moving basis · substitution in an integral ·
np.expm1 and cancellation near zero · scipy.optimize.brentq · differentiation under the integral sign · signed swept
volume of a moving surface. **Reminders** of ch01/ch02 primers: P01, P04, P10, P13, P14, P15, P16, P17, P18, P19, P21,
P23, P25, P26, P27, P29, P30, P31, P36, P37, P38, P40, P43, P46, P49, P56, P62–P70 (P63, P64, P65, P66, P67, P68, P70),
P73, P75, P76, P78, P79, P80, P82, P83, P84, P85, P86. **Glosses:** 17 (ring area, sech², integrating sin/cos, sin² + cos²,
circle equation, implicit differentiation, PIV, δ_jk, n·A·n = 0, product rule for three factors, Jacobi, planetary
vorticity, simply connected, ellipse equation, Lambert W, Taylor in time, ∇·(Fu)).

---

## Part F — derivation storyboards

Builders copy these word for word into `nb.derivation(key, title, goal=…, start=(tex, plain), plan=[…], uses=[…],
steps=[dict(did, tex, why, plain)], result=(tex, plain), interpret=…, check=…, check_src=…)` and into the explainer's
`derivations: [...]` (phones may shorten *why* to its first sentence; `live` and `set` are the explainer's). Every step is
one move; *why* names the rule and says why we make the move (≤ 35 words, ≥ 6); *did* ≤ 8 words. The book's own moves
were read on the rendered pages (p097 for D02, p098 for D03, p100 for D04/D05, p101 for D06, p103–p105 for D07–D09 and
D11, p106 for D12, p107 for D15, p108 for D16, p109 for D14, p110–p111 for D17, p112 for D18–D20, p113 for D21, p114–p115
for D22); the gaps listed in `analysis/ch03.md` §2b are filled and the notebook says so ("the book skips this move").
Colours: local blue, advective amber, total purple, strain teal, rotation orange, stream/path/streak teal/orange/rose,
volume term blue, surface term orange. Every equation named by number is written out. (No line in this part starts with
a table bar; absolute values are written with \lvert \rvert or in words.)

### D01 · The Eulerian velocity of a Lagrangian map: x = X e^{αt} ⇒ u = αx — ★, 6 steps, in C01 (notebook · `material_derivative_probe`)
- **Goal.** We know where every particle is at every time (a Lagrangian map). Find the velocity *field* a fixed probe
  would record, and the field of accelerations — the concrete version of the bridge (3.2) that the book never works out.
- **Start.** $x=r(t;X)=X\,e^{\alpha t}$ — *in words:* the particle labelled X (its position at t = 0) moves away from the
  origin at a rate proportional to its distance; α [1/s] is a constant stretching rate.
- **Plan.** (1) Differentiate the path with the label held fixed, (3.1). (2) Ask which particle is at x at time t
  (invert the map). (3) Substitute that label: the bridge (3.2). (4) Do the same for the acceleration.
- **Tools.** (3.1), $u=dr/dt$, $a=d^2r/dt^2$ (N08) · (3.2), $F[\mathbf r(t;\mathbf r_o,t_o),t]=F(\mathbf x,t)$ at
  $\mathbf x=\mathbf r$ (C01) · functions of time with parameters (primer in C01) · inverse functions and sympy solve
  (primer in C01) · exponent rules (ch01 P43).
- **Assumptions.** The map is invertible — two particles never occupy the same point (step 3); r is twice
  differentiable in t (steps 1–2).
- **Steps.**
  1. *did:* Differentiate the path, label fixed · *tex:* $u=\left.\dfrac{dx}{dt}\right|_X=\alpha X e^{\alpha t}$ · *why:*
     (3.1), $\mathbf u=d\mathbf r/dt$: a particle's velocity is the time derivative of its own path; X names the particle, so
     it is held constant. · *plain:* The particle labelled X speeds up exponentially. · *live:* "u = 0.5 × 2 × e^{0.5} =
     1.649 m/s".
  2. *did:* Differentiate once more · *tex:* $a=\left.\dfrac{d^2x}{dt^2}\right|_X=\alpha^2Xe^{\alpha t}$ · *why:* (3.1),
     $\mathbf a=d^2\mathbf r/dt^2$; differentiating $e^{\alpha t}$ brings down another factor α (chain rule for the
     exponential). · *plain:* Its acceleration also grows exponentially.
  3. *did:* Invert the map for the label · *tex:* $X=x\,e^{-\alpha t}$ · *why:* To use (3.2) we need the label of the
     particle that sits at x at time t; divide $x=Xe^{\alpha t}$ by $e^{\alpha t}\neq0$ — one answer, so the map is
     invertible. · *plain:* The particle now at x started at $xe^{-\alpha t}$. · *live:* "X = 3.297 × e^{−0.5} = 2.000 m".
  4. *did:* Substitute the label into u · *tex:* $u(x,t)=\alpha\,(xe^{-\alpha t})\,e^{\alpha t}$ · *why:* The bridge (3.2),
     $F[\mathbf r(t;\mathbf r_o,t_o),t]=F(\mathbf x,t)$ when $\mathbf x=\mathbf r$: the field at x is the velocity carried by
     the particle that is there now. · *plain:* The Eulerian velocity written with the particle's label replaced by its
     position.
  5. *did:* Simplify the exponentials · *tex:* $u(x,t)=\alpha x$ · *why:* $e^{-\alpha t}e^{\alpha t}=e^0=1$ (exponent rules,
     P43); t has disappeared — the field does not depend on time. · *plain:* A steady field: every point always has
     velocity α times its distance. · *live:* "u = 0.5 × 3.297 = 1.649 m/s".
  6. *did:* Same substitution for a · *tex:* $a(x,t)=\alpha^2x$ · *why:* Apply (3.2) to the acceleration of step 2 with
     the label of step 3; the exponentials cancel as in step 5. · *plain:* The acceleration field is also steady and
     never zero away from the origin.
- **Result.** $u(x,t)=\alpha x$, $a(x,t)=\alpha^2x$ — *in words:* the Eulerian velocity field is steady even though
  every particle accelerates.
- **Check.** Units: α [1/s] × x [m] = m/s ✓; α² x = m/s² ✓. Limit α = 0: nobody moves, u = a = 0 ✓. Numbers (X = 2 m,
  α = 0.5 s⁻¹, t = 1 s): x = 3.297 m, u = 1.649 m/s = αx ✓, a = 0.824 m/s² ✓. Forward check (C02): the material
  derivative $\frac{Du}{Dt}=\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x}=0+\alpha x\cdot\alpha=\alpha^2x$
  reproduces step 6 from the field alone (the notebook's sympy line).
- **What it means.** "Steady" is about points, "accelerating" about particles: a steady field can accelerate every
  particle. The missing link — how to get a from the Eulerian u without knowing the paths — is the material derivative
  (D02). The recipe fails if the map is not invertible (particles colliding, a shock).
- **Traps.** Differentiating $u=\alpha x$ with respect to t at fixed x gives 0 — that is ∂u/∂t, not the acceleration.
  Forgetting that X is constant along the path. Writing $u=\alpha Xe^{\alpha t}$ as "the field" — it is a function of the
  label, not of position.

### D02 · The material derivative: (3.3) → (3.4) → (3.5) — ★★, 9 steps, in C02 (notebook · `material_derivative_probe`)
- **Goal.** Find the rate at which a property F changes for an observer riding with a fluid particle, written only in
  terms of the Eulerian field F(x, t) and the velocity field u(x, t).
- **Start.** (3.2): $F[\mathbf r(t;\mathbf r_o,t_o),t]=F(\mathbf x,t)$ when $\mathbf x=\mathbf r(t;\mathbf r_o,t_o)$ — *in
  words:* the value a particle carries equals the field's value at the point it occupies.
- **Plan.** (1) Fix one particle, so F along its path is a function of t alone. (2) Differentiate with the chain rule —
  time enters through position and explicitly. (3) Recognise the particle's velocity and the field's slopes. (4) Argue
  the result holds at every point and name it.
- **Tools.** (3.2) (C01) · multivariable chain rule along a path (primer in C02) · (3.1), $\mathbf u=d\mathbf r/dt$ (N08)
  · partial derivative (ch01 P25) · summation convention and the gradient (Ch. 2 §2.1, §2.9).
- **Assumptions.** Continuum: at every instant every point is occupied by exactly one particle (step 6); F
  differentiable in x and t (step 3).
- **Steps.**
  1. *did:* Follow one particle · *tex:* $f(t)\equiv F[\mathbf r(t;\mathbf r_o,t_o),t]$ · *why:* The labels $\mathbf r_o$,
     $t_o$ are fixed for a chosen particle, so what its thermometer reads is an ordinary function of one variable, t. We
     want df/dt. · *plain:* The reading of a thermometer riding with the particle.
  2. *did:* List what depends on t · *tex:* $f(t)=F\big(r_1(t),r_2(t),r_3(t),t\big)$ · *why:* F has four arguments —
     three coordinates and time — and along the path all four change with t: the position because the particle moves,
     the last because the clock runs. · *plain:* Time enters twice: through where the particle is, and directly.
  3. *did:* Apply the chain rule · *tex:* $\dfrac{df}{dt}=\dfrac{\partial F}{\partial r_i}\dfrac{dr_i}{dt}+\dfrac{\partial F}{\partial t}$
     · *why:* Multivariable chain rule (primer): each changing argument contributes its sensitivity times its rate; the
     repeated i sums the three coordinates. This is (3.3) written out. · *plain:* The reading changes because the
     particle moves through the field and because the field itself changes.
  4. *did:* Recognise the particle's velocity · *tex:* $\dfrac{df}{dt}=\dfrac{\partial F}{\partial r_i}\,u_i+\dfrac{\partial F}{\partial t}$
     · *why:* (3.1), $\mathbf u=d\mathbf r/dt$: the rate of change of the particle's position components is its velocity
     components. · *plain:* The first term is 'how fast we move through the field'. · *live:* "u₁ = 0.0, u₂ = 10.0 m/s".
  5. *did:* Evaluate the slopes at x = r · *tex:* $\dfrac{df}{dt}=u_i\dfrac{\partial F}{\partial x_i}+\dfrac{\partial F}{\partial t}$
     · *why:* At time t the particle sits at $\mathbf x=\mathbf r$, so $\partial F/\partial r_i$ is the field's slope
     in direction i there; renaming the argument changes nothing. The book states this. · *plain:*
     Everything on the right is now an Eulerian quantity at (x, t).
  6. *did:* Extend to every point · *tex:* $\dfrac{df}{dt}\Big|_{\mathbf r=\mathbf x}=u_i\dfrac{\partial F}{\partial x_i}+\dfrac{\partial F}{\partial t}\quad\text{for all }\mathbf x,t$
     · *why:* The book skips this: at any instant every point is occupied by some particle (continuum, Ch. 1), so
     choosing that particle's labels makes the line hold at any x. · *plain:* The rule is a statement about the whole
     field, not one path. · *set:* moves the float to a new start.
  7. *did:* Give it a name · *tex:* $\dfrac{DF}{Dt}\equiv\dfrac{\partial F}{\partial t}+u_i\dfrac{\partial F}{\partial x_i}$ ·
     *why:* A symbol for 'd/dt following the particle, in Eulerian variables' — the book's (3.4) defines D/Dt exactly
     this way (material, substantial or particle derivative). · *plain:* D/Dt means: the rate seen by a moving parcel.
  8. *did:* Write the sum as a dot product · *tex:* $\dfrac{DF}{Dt}\equiv\dfrac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$
     · *why:* Summation convention (Ch. 2 §2.1): $u_i\,\partial F/\partial x_i$ is the dot product of u with the gradient
     $\nabla F$ (Ch. 2 §2.9). This is (3.5). · *plain:* Local rate plus the advective rate u·∇F. · *live:* "u·∇T = 10 ×
     (−1.0e-5) = −1.0e-4 K/s = −0.36 K/h".
  9. *did:* Apply it to each velocity component · *tex:* $\dfrac{D\mathbf u}{Dt}=\dfrac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u$
     · *why:* Steps 1–8 hold for any scalar field, so apply them with $F=u_j$ for j = 1, 2, 3 separately; the same
     operator acts on each component. We need this form in C05. · *plain:* A particle's acceleration = local + advective
     acceleration.
- **Result.** $\dfrac{DF}{Dt}\equiv\dfrac{\partial F}{\partial t}+\mathbf u\cdot\nabla F=\dfrac{\partial F}{\partial t}+u_i\dfrac{\partial F}{\partial x_i}$
  (3.5), and $\frac{D\mathbf u}{Dt}=\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u$ — *in words:* what a moving
  parcel feels = what a fixed probe sees + what it meets by moving through the gradient.
- **Check.** Units: both terms [F]/s ✓. Limits: F steady ⇒ DF/Dt = u·∇F (all advective); u = 0 ⇒ DF/Dt = ∂F/∂t; u ⟂ ∇F
  ⇒ no advective part ✓. Numbers (C02 worked example): ∂T/∂t = +1.0 × 10⁻⁴, u·∇T = −1.0 × 10⁻⁴, DT/Dt = 0 K/s ✓. D01's
  field: F = u = αx gives Du/Dt = 0 + αx·α = α²x — the acceleration found there ✓. Optional `check_src`: the sympy
  chain-rule check of the primer applied to $F=T_0-G(y-ct)$ along $y=vt$.
- **What it means.** Every conservation law of Ch. 4 is written with D/Dt because the laws are about particles while
  the equations are solved on fixed grids. The advective term is where fluid mechanics becomes nonlinear (C05, N22).
  The formula needs a continuum; at a point where particles do not exist (a free surface's edge, a shock) the
  derivative is taken one-sidedly.
- **Traps.** Reading ∂T/∂t as 'the air warms' (it is the fixed station). Holding the label fixed vs holding x fixed:
  df/dt in step 1 is at fixed label, ∂F/∂t at fixed x. Treating D/Dt of a vector as something new — it acts component by
  component in Cartesian coordinates (curvilinear components need Appendix B's extra terms).

### D03 · The streamline equations dx/u = dy/v = dz/w, Eq. (3.7) — ★, 5 steps, in C03 (notebook · `flow_lines_unsteady`)
- **Goal.** Turn "a curve that is tangent to the velocity everywhere at one instant" into equations we can integrate —
  the book leaves this to Exercise 3.3.
- **Start.** $d\mathbf s\parallel\mathbf u(\mathbf x,t)$ at a fixed t, with $d\mathbf s=(dx,dy,dz)$ a small step along the
  curve and $\mathbf u=(u,v,w)$ — *in words:* each little piece of the streamline points along the local velocity.
- **Plan.** (1) Write parallel as "one is a multiple of the other". (2) Take components. (3) Eliminate the multiple. (4)
  Fix the multiple by arc length to get an ODE we can integrate.
- **Tools.** Parallel vectors and the cross-product test (primer in C03) · parametric curves, tangent vector and arc
  length (primer in C03).
- **Assumptions.** The clock is frozen at t (every step); $\mathbf u\neq0$ at the points used (step 3 divides by
  components, step 4 by |u|).
- **Steps.**
  1. *did:* Write parallel as a multiple · *tex:* $d\mathbf s=\lambda\,\mathbf u$ · *why:* Two vectors are parallel when one
     is a scalar multiple of the other (primer); λ > 0 so the curve runs downstream, and λ has units of time (m ÷ m/s). ·
     *plain:* The small step is the velocity times a tiny 'time-like' number.
  2. *did:* Split into components · *tex:* $dx=\lambda u,\quad dy=\lambda v,\quad dz=\lambda w$ · *why:* Two vectors are
     equal when each component is equal. · *plain:* Each coordinate step is proportional to that velocity component.
  3. *did:* Eliminate λ · *tex:* $\dfrac{dx}{u}=\dfrac{dy}{v}=\dfrac{dz}{w}\;(=\lambda)$ · *why:* Solve each line for λ by
     dividing by u, v, w (allowed only where they are nonzero); all three equal the same λ. This is (3.7). Equivalently
     $\mathbf u\times d\mathbf s=0$, which never divides. · *plain:* The direction numbers of the step match those of
     the velocity. · *live:* "dy/dx = v/u = 0.500/0.866 = 0.577".
  4. *did:* Fix λ by the step's length · *tex:* $\lambda=\dfrac{ds}{\lvert\mathbf u\rvert}$ · *why:* Take the length of both
     sides of step 1: $ds=\lambda\lvert\mathbf u\rvert$ with ds the arc length (primer). Choosing arc length makes the next
     line well defined everywhere u ≠ 0. · *plain:* The step's length decides how big λ is.
  5. *did:* Write the arc-length ODE · *tex:* $\dfrac{d\mathbf x}{ds}=\dfrac{\mathbf u(\mathbf x,t)}{\lvert\mathbf u(\mathbf x,t)\rvert}$
     · *why:* Substitute λ into step 1 and divide by ds. The right side is a known unit vector field at frozen t, so this
     is an ODE `ch03.streamline` integrates both ways from a seed. · *plain:* Walk along the unit arrows at one instant.
- **Result.** $dx/u=dy/v=dz/w$ (3.7), equivalently $\mathbf u\times d\mathbf s=0$ and $d\mathbf x/ds=\mathbf u/\lvert\mathbf u\rvert$
  at frozen t — *in words:* a streamline follows the arrows of one instant.
- **Check.** Units: each ratio in s ✓. Plane flow: $dy/dx=v/u$; with u = (1, 2) the slope is 2 ✓. Solid-body rotation
  u = −y, v = x: $dy/dx=-x/y$ ⇒ $x^2+y^2=$ const, circles ✓. Ex. 3.1 at ωt′ = 30°: slope tan 30° = 0.577 ✓.
- **What it means.** At each instant the flow has a streamline pattern; in unsteady flow it changes, so streamlines are
  not where particles go (C04). The ODE stops at a stagnation point (u = 0, no direction) — `streamline` ends the curve
  there.
- **Traps.** Integrating (3.7) forward in *time* (that gives a path line, not a streamline — t is frozen). Using the
  ratio form where a component vanishes (use $\mathbf u\times d\mathbf s=0$ or the arc-length ODE instead).

### D04 · Ex. 3.1: the path line is a circle of radius ξ_o — ★★, 7 steps, in C04 (notebook · `flow_lines_unsteady`)
- **Goal.** Find the path of the particle that sits at the origin at the drawing instant t′ in the sloshing flow
  $u=\omega\xi_o\cos\omega t$, $v=\omega\xi_o\sin\omega t$ — and fill in the "little algebra" the book skips.
- **Start.** (3.8), $\dfrac{d\mathbf r}{dt}=\mathbf u(\mathbf r,t)$, $\mathbf r(t')=\mathbf 0$ — *in words:* the particle
  moves with the field's velocity at its current position and passes the origin at t′.
- **Plan.** (1) Write both components and integrate them (the field does not depend on position, so each equation
  integrates directly). (2) Fix the constants with the condition at t′. (3) Eliminate t with sin² + cos² = 1.
- **Tools.** (3.8) (C04) · integrating sin and cos (gloss in step 2) · eliminating a parameter with sin² + cos² = 1
  (gloss in steps 6–7) · equation of a circle (gloss in step 7).
- **Assumptions.** The field is spatially uniform (step 1: the right-hand sides depend on t only); ω ≠ 0.
- **Steps.**
  1. *did:* Write the two components · *tex:* $\dfrac{dx}{dt}=\omega\xi_o\cos\omega t,\quad\dfrac{dy}{dt}=\omega\xi_o\sin\omega t$
     · *why:* (3.8) component by component. The velocity is the same at every point, so the right sides are known
     functions of t alone — no coupling to x, y. · *plain:* The particle's velocity is just the field's velocity at that
     time.
  2. *did:* Integrate each once · *tex:* $x=\xi_o\sin\omega t+c_x,\quad y=-\xi_o\cos\omega t+c_y$ · *why:* An antiderivative
     of $\omega\cos\omega t$ is $\sin\omega t$ and of $\omega\sin\omega t$ is $-\cos\omega t$ (the ω cancels); $c_x,c_y$ are
     the constants (the book's x_o, y_o). · *plain:* Every particle traces the same shape; only the constants differ.
  3. *did:* Impose r(t′) = 0 · *tex:* $c_x=-\xi_o\sin\omega t',\quad c_y=\xi_o\cos\omega t'$ · *why:* The particle we want
     is at the origin at the drawing instant: set t = t′ and x = y = 0 in step 2 and solve for the constants. · *plain:*
     The condition picks one particle. · *live:* "c_x = −1 × sin 0 = 0, c_y = 1 × cos 0 = 1 m".
  4. *did:* Substitute the constants · *tex:* $x=\xi_o(\sin\omega t-\sin\omega t'),\quad y=\xi_o(\cos\omega t'-\cos\omega t)$ ·
     *why:* Put step 3 into step 2. These are the path's parametric equations, with t the parameter — the book prints
     them. · *plain:* Where that particle is at any time t.
  5. *did:* Isolate the parts that depend on t · *tex:* $x+\xi_o\sin\omega t'=\xi_o\sin\omega t,\quad y-\xi_o\cos\omega t'=-\xi_o\cos\omega t$
     · *why:* The book skips this: move the constants to the left so each right side holds only t — the form in which
     squaring and adding removes t. · *plain:* Offsets from a fixed point equal a sine and a cosine of the same angle.
  6. *did:* Square both and add · *tex:* $(x+\xi_o\sin\omega t')^2+(y-\xi_o\cos\omega t')^2=\xi_o^2(\sin^2\omega t+\cos^2\omega t)$
     · *why:* Squaring removes the minus sign; adding pairs sin² with cos² of the *same* angle ωt, ready for the identity.
     · *plain:* The squared distance from the fixed point.
  7. *did:* Use sin² + cos² = 1 · *tex:* $(x+\xi_o\sin\omega t')^2+(y-\xi_o\cos\omega t')^2=\xi_o^2$ · *why:* Pythagorean
     identity; t has gone, so this is the curve itself. $(x-a)^2+(y-b)^2=r^2$ is a circle of centre (a, b) and radius r.
     · *plain:* A circle of radius ξ_o centred at $(-\xi_o\sin\omega t',\ \xi_o\cos\omega t')$. · *set:* t = 2π (the trail
     closes).
- **Result.** $(x+\xi_o\sin\omega t')^2+(y-\xi_o\cos\omega t')^2=\xi_o^2$ — *in words:* the particle runs round a circle
  of radius ξ_o that touches the origin, once per period 2π/ω.
- **Check.** Units: every term m² ✓. The origin lies on it: $\xi_o^2\sin^2+\xi_o^2\cos^2=\xi_o^2$ ✓. Numbers (ξ_o = 1 m,
  ω = 1 s⁻¹, t′ = 0): centre (0, 1) m ✓. Direction: at t = 0 the particle is at the origin (the bottom of that circle)
  moving with u = (ωξ_o, 0) — to the right — so it runs counterclockwise ✓.
- **What it means.** In a uniform but rotating velocity field every particle circles: a caricature of the orbital
  motion under long surface waves (Ch. 7 has the real ones, which shrink with depth). The circle does not coincide with
  the straight streamline — the flow is unsteady.
- **Traps.** Confusing t′ (drawing instant), t_o (release time, D05), c_x, c_y (constants) and ξ_o (amplitude). Squaring
  before isolating the t-terms (the constants mix in). Forgetting the direction of travel.

### D05 · Ex. 3.1: the streak line is the mirror circle, and all three lines touch at the origin — ★★, 8 steps, in C04 (notebook · `flow_lines_unsteady`)
- **Goal.** Find where the dye that has been leaking from the origin sits at the instant t′, and prove the book's claim
  that the streamline, path line and streak line are tangent at the origin — both skipped in the book.
- **Start.** D04 steps 2–3 with the release time $t_o$ in place of t′: $x=\xi_o(\sin\omega t-\sin\omega t_o)$,
  $y=\xi_o(\cos\omega t_o-\cos\omega t)$ — *in words:* the path of the particle that was at the origin at time t_o.
- **Plan.** (1) Freeze the clock at t′ and let the release time vary. (2) Eliminate t_o as in D04. (3) Compare the three
  slopes at the origin.
- **Tools.** D04 · eliminating a parameter with sin² + cos² = 1 (gloss) · implicit differentiation (gloss in step 7) ·
  the Ex. 3.1 streamline $y=x\tan(\omega t')$ (N18).
- **Assumptions.** Dye has been released continuously for at least one period before t′ (step 6).
- **Steps.**
  1. *did:* Label particles by their release time · *tex:* $x=\xi_o(\sin\omega t-\sin\omega t_o),\ \ y=\xi_o(\cos\omega t_o-\cos\omega t)$
     · *why:* Same integration as D04 steps 1–2, but the constants are fixed by $\mathbf r(t_o)=\mathbf 0$: every release
     time labels a different dye particle. · *plain:* Where the dye released at t_o is at any later time t.
  2. *did:* Freeze the clock at t′ · *tex:* $x=\xi_o(\sin\omega t'-\sin\omega t_o),\ \ y=\xi_o(\cos\omega t_o-\cos\omega t')$ ·
     *why:* A streak line is a snapshot at t = t′ of all dye released at earlier times $t_o\le t'$ (N17): now t_o is the
     curve's parameter. · *plain:* The photograph of the dye at t′. · *set:* colours dye dots by release time.
  3. *did:* Isolate the t_o-terms · *tex:* $x-\xi_o\sin\omega t'=-\xi_o\sin\omega t_o,\quad y+\xi_o\cos\omega t'=\xi_o\cos\omega t_o$
     · *why:* Move the known terms to the left; the parameter to be eliminated — t_o, not t — stays on the right. ·
     *plain:* Offsets from a fixed point, now in terms of the release angle.
  4. *did:* Square both and add · *tex:* $(x-\xi_o\sin\omega t')^2+(y+\xi_o\cos\omega t')^2=\xi_o^2(\sin^2\omega t_o+\cos^2\omega t_o)$
     · *why:* As in D04 step 6, squaring removes signs and pairs sin² with cos² of the same angle $\omega t_o$. · *plain:*
     The squared distance of each dye particle from one fixed point.
  5. *did:* Use sin² + cos² = 1 · *tex:* $(x-\xi_o\sin\omega t')^2+(y+\xi_o\cos\omega t')^2=\xi_o^2$ · *why:* Pythagorean
     identity, as in D04 step 7; the release time t_o has disappeared, so this is the curve itself. · *plain:* A circle of radius ξ_o centred at $(\xi_o\sin\omega t',-\xi_o\cos\omega t')$ — the
     path circle reflected through the origin.
  6. *did:* Check which part is filled · *tex:* $t'-2\pi/\omega\le t_o\le t'\ \Rightarrow\ \text{whole circle}$ · *why:* The
     book is silent here: as t_o runs over one period, $(\sin\omega t_o,\cos\omega t_o)$ goes once round; older dye repeats
     the same points because the flow is periodic. · *plain:* After one period of leaking, the dye covers the full circle.
     · *live:* "releases over T = 6.283 s fill it".
  7. *did:* Slope of the path circle at 0 · *tex:* $\left.\dfrac{dy}{dx}\right|_{\rm path}=-\dfrac{x+\xi_o\sin\omega t'}{y-\xi_o\cos\omega t'}\Big|_{\mathbf 0}=\tan\omega t'$
     · *why:* Implicit differentiation of D04's circle: $2(x+\dots)dx+2(y-\dots)dy=0$; at the origin the ratio is
     $-\xi_o\sin\omega t'/(-\xi_o\cos\omega t')$. · *plain:* The path leaves the origin at angle ωt′. · *live:* "−0.500/(−0.866)
     = 0.577".
  8. *did:* Slope of the streak circle at 0 · *tex:* $\left.\dfrac{dy}{dx}\right|_{\rm streak}=-\dfrac{x-\xi_o\sin\omega t'}{y+\xi_o\cos\omega t'}\Big|_{\mathbf 0}=\tan\omega t'$
     · *why:* Implicit differentiation of step 5 at the origin gives $\xi_o\sin\omega t'/(\xi_o\cos\omega t')$; the
     streamline $y=x\tan(\omega t')$ has the same slope. · *plain:* All three curves touch at the origin, along the velocity
     there.
- **Result.** $(x-\xi_o\sin\omega t')^2+(y+\xi_o\cos\omega t')^2=\xi_o^2$, and the three lines share the slope tan ωt′ at
  the origin — *in words:* the dye lies on the mirror image of the particle's circle; streamline, path line and streak
  line only touch.
- **Check.** Units m² ✓. Numbers (ξ_o = 1, ωt′ = 0): centre (0, −1) — opposite the path centre (0, 1) ✓; the streak
  circle passes the origin ✓. `ch03.streakline` of 181 releases lies within 1e-7 of this circle ✓.
- **What it means.** A photograph of dye (streak line) and a long exposure of one speck (path line) are different
  curves in an unsteady flow, and neither is the instantaneous direction (streamline). They share only their tangent at
  the port — the velocity there at t′.
- **Traps.** Eliminating t instead of t_o. Thinking the streak line is only a partial arc (a full period of release
  fills it). Mixing up which circle is which: path centre (−ξ_o sin ωt′, ξ_o cos ωt′), streak centre (ξ_o sin ωt′,
  −ξ_o cos ωt′).

### D06 · Galilean invariance of the acceleration, Eq. (3.9) — ★★, 9 steps, in C05 (notebook · `galilean_frames_cylinder`)
- **Goal.** Show that local + advective acceleration is the same in two frames moving at constant relative velocity,
  although each term separately is not — the book leaves it to Exercise 3.12 ("it can be shown").
- **Start.** $\mathbf u(\mathbf x,t)=\mathbf U+\mathbf u'(\mathbf x',t')$, $t=t'$, $\mathbf x=\mathbf x'+\mathbf Ut+\mathbf x'_o$,
  with U and $\mathbf x'_o$ constant — *in words:* the primed frame moves at constant velocity U; velocities differ by U,
  positions by Ut plus a fixed offset, the clocks agree.
- **Plan.** (1) Express the primed coordinates through the unprimed ones and differentiate them. (2) Chain rule for the
  local term (the trap). (3) Chain rule for the space derivatives (they do not change). (4) Add: the two U-terms cancel.
- **Tools.** Galilean transformation (N21) · chain rule with a moving frame (primer in C05) · $\partial x_j/\partial x_k=\delta_{jk}$
  (Ch. 2 §2.7) · the acceleration $\frac{D\mathbf u}{Dt}=\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u$
  (D02 step 9) · summation convention.
- **Assumptions.** U constant in time (step 3: $\partial\mathbf U/\partial t=0$) and in space (step 5); axes parallel and
  not rotating (step 2: $\partial x'_j/\partial x_k=\delta_{jk}$); one common clock.
- **Steps.**
  1. *did:* Invert the transformation · *tex:* $x'_j=x_j-U_jt-x'_{o,j},\qquad t'=t$ · *why:* The field u′ is known in
     primed coordinates; to differentiate u with respect to x and t we need how x′, t′ depend on them — solve the start
     line for x′. · *plain:* A primed coordinate is the unprimed one minus how far O′ has moved.
  2. *did:* Differentiate the primed coordinates · *tex:* $\dfrac{\partial x'_j}{\partial t}=-U_j,\quad\dfrac{\partial x'_j}{\partial x_k}=\delta_{jk},\quad\dfrac{\partial t'}{\partial t}=1$
     · *why:* Differentiate step 1 term by term: U and $\mathbf x'_o$ are constants, $\partial x_j/\partial x_k=\delta_{jk}$
     (Ch. 2 §2.7), and t′ does not depend on x. These are the chain-rule factors. · *plain:* At a fixed lab point, the
     primed coordinate slides backwards at U.
  3. *did:* Chain rule for the local term · *tex:* $\dfrac{\partial u_i}{\partial t}=\dfrac{\partial u'_i}{\partial t'}-U_j\dfrac{\partial u'_i}{\partial x'_j}$
     · *why:* $u_i=U_i+u'_i(x',t')$ with U constant; t enters u′ through t′ (factor 1) and through each $x'_j$ (factor
     $-U_j$) — the moving-frame chain rule (primer). The book skips this. · *plain:* A fixed lab probe sees the primed
     pattern slide past. · *live:* "0 − (−1.0) × (0, −0.593) = (0, −0.593) m/s² (the body frame moves at U = (−1, 0) m/s relative to the lake)".
  4. *did:* Chain rule for space derivatives · *tex:* $\dfrac{\partial u_i}{\partial x_k}=\dfrac{\partial u'_i}{\partial x'_j}\delta_{jk}=\dfrac{\partial u'_i}{\partial x'_k}$
     · *why:* x enters only through x′ with factor $\delta_{jk}$ (step 2); the Kronecker delta substitutes its index. So
     $\nabla=\nabla'$ for a pure translation. · *plain:* Slopes in space are the same in both frames.
  5. *did:* Write the advective term · *tex:* $u_k\dfrac{\partial u_i}{\partial x_k}=(U_k+u'_k)\dfrac{\partial u'_i}{\partial x'_k}$ ·
     *why:* Substitute $u_k=U_k+u'_k$ (the transformation) and step 4 into $(\mathbf u\cdot\nabla)u_i$. · *plain:* The lab
     advects with the full velocity, U included.
  6. *did:* Expand the product · *tex:* $u_k\dfrac{\partial u_i}{\partial x_k}=U_k\dfrac{\partial u'_i}{\partial x'_k}+u'_k\dfrac{\partial u'_i}{\partial x'_k}$ ·
     *why:* Distribute the sum over k. The first piece is the extra advection by the frame speed; we expect it to meet
     step 3's extra term. · *plain:* Advective = advection by U + advection by u′. · *live:* "(−1.0) × (0, −0.593) + (0, −0.856) = (0, −0.263) m/s²".
  7. *did:* Add local and advective terms · *tex:* $\dfrac{\partial u_i}{\partial t}+u_k\dfrac{\partial u_i}{\partial x_k}=\dfrac{\partial u'_i}{\partial t'}-U_j\dfrac{\partial u'_i}{\partial x'_j}+U_k\dfrac{\partial u'_i}{\partial x'_k}+u'_k\dfrac{\partial u'_i}{\partial x'_k}$
     · *why:* The left side is $Du_i/Dt$ in Oxyz (D02 step 9); add step 3 and step 6. · *plain:* Four terms: two of them
     contain U.
  8. *did:* Cancel the two U-terms · *tex:* $\dfrac{\partial u_i}{\partial t}+u_k\dfrac{\partial u_i}{\partial x_k}=\dfrac{\partial u'_i}{\partial t'}+u'_k\dfrac{\partial u'_i}{\partial x'_k}$
     · *why:* j and k are dummy indices, so $U_j\,\partial u'_i/\partial x'_j$ and $U_k\,\partial u'_i/\partial x'_k$ are the same
     sum with opposite signs. · *plain:* What the lab calls 'local' and the frame's advection trade places and cancel. ·
     *set:* shows the two hatched half-bars.
  9. *did:* Return to vector notation · *tex:* $\dfrac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=\dfrac{\partial\mathbf u'}{\partial t'}+(\mathbf u'\cdot\nabla')\mathbf u'$
     · *why:* Summation convention read backwards: the index line for each i is the vector equation (3.9). · *plain:* The
     particle's acceleration is the same for both observers. · *set:* U_f = 0.5 (the purple bar does not move).
- **Result.** $\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=\left(\frac{D\mathbf u}{Dt}\right)_{Oxyz}=\left(\frac{D\mathbf u'}{Dt'}\right)_{O'x'y'z'}=\frac{\partial\mathbf u'}{\partial t'}+(\mathbf u'\cdot\nabla')\mathbf u'$
  (3.9) — *in words:* observers in uniform relative motion agree on every particle's acceleration.
- **Check.** Units m/s² ✓. U = 0: identity ✓. Numbers: lab sine wave local −0.707 + advective 0.832 = 0.125 m/s², wave
  frame 0 + 0.125 ✓; cylinder at P = (0, 1.5) m: fluid frame (0, −0.593) + (0, −0.263) = body frame (0, −0.856) ✓.
- **Sympy check intent (optional `check_src`, recommended):** build a general 2-D primed field from two
  `sp.Function`s, form u = U + u′(x − Ut, y − U₂t, t), compute ∂u/∂t + (u·∇)u and the primed acceleration substituted at
  the same point, simplify the difference → (0, 0):
  ```python
  import sympy as sp                                        # symbolic algebra
  x, y, t, U1, U2 = sp.symbols('x y t U_1 U_2')             # lab coordinates, time, frame velocity
  xp, yp, tp = sp.symbols('xp yp tp')                       # primed coordinates
  f, g = sp.Function('f'), sp.Function('g')                 # u' = (f, g): any smooth primed field
  up = [f(xp, yp, tp), g(xp, yp, tp)]                       # u'(x', t')
  sub = {xp: x - U1*t, yp: y - U2*t, tp: t}                 # x' = x − Ut (offset 0), t' = t
  u = [U1 + up[0].subs(sub), U2 + up[1].subs(sub)]          # u = U + u'  (step 0: the transformation)
  a_lab = [sp.diff(q, t) + u[0]*sp.diff(q, x) + u[1]*sp.diff(q, y) for q in u]   # ∂u/∂t + (u·∇)u
  a_pr = [(sp.diff(q, tp) + up[0]*sp.diff(q, xp) + up[1]*sp.diff(q, yp)).subs(sub) for q in up]  # primed, same point
  print([sp.simplify((a - b).doit()) for a, b in zip(a_lab, a_pr)])   # [0, 0]: (3.9)
  ```
- **What it means.** Newton's law can be written in any frame moving steadily — wind tunnels, wave frames and moving
  coordinate systems are legitimate. The split into 'unsteady' and 'advective' is a choice of observer. It fails for
  accelerating or rotating frames: then step 3 gains −dU/dt, or U = Ω × x varies in space and step 4 changes — Coriolis
  and centrifugal terms (Ch. 4 §4.7).
- **Traps.** Writing $\partial\mathbf u/\partial t=\partial\mathbf u'/\partial t'$ (forgets the $-\mathbf U\cdot\nabla'$
  term — the classic slip). Forgetting that U also advects in step 5. Applying the result to a rotating frame.

### D07 · Relative velocity near a point, Eq. (3.10) — ★, 4 steps, in C06 (notebook · `fluid_element_deformation`)
- **Goal.** Find how fast a neighbouring point moves relative to a given point, to first order in their separation —
  the relation everything in §3.4 is built on.
- **Start.** $u_i(\mathbf x+d\mathbf x,t)$ — *in words:* the velocity at the neighbour P, a small step dx away from O
  (Fig. 3.9), at the same instant.
- **Plan.** (1) Taylor-expand the velocity at P about O. (2) Subtract O's velocity. (3) Drop what vanishes faster than the
  separation itself. (4) Recognise a matrix–vector product.
- **Tools.** Multivariable first-order Taylor expansion (primer in C06) · orders of smallness (ch02 P68) · the velocity
  gradient G (R01) · matrix–vector product (ch02 P63).
- **Assumptions.** u is smooth (twice differentiable) near x (step 1); |dx| small compared with the distance over which
  G changes (step 3).
- **Steps.**
  1. *did:* Taylor-expand u about x · *tex:* $u_i(\mathbf x+d\mathbf x)=u_i(\mathbf x)+\dfrac{\partial u_i}{\partial x_j}dx_j+O(\lvert d\mathbf x\rvert^2)$
     · *why:* Multivariable Taylor expansion (primer): one first derivative per direction, summed over j; the remainder is
     of second order because u is smooth. · *plain:* Near O the velocity changes linearly with position. · *live:* "(0.0501,
     −0.0198) = (0, 0) + (0.05, −0.02) + (1, 2) × 10⁻⁴".
  2. *did:* Subtract the velocity at O · *tex:* $du_i\equiv u_i(\mathbf x+d\mathbf x)-u_i(\mathbf x)=\dfrac{\partial u_i}{\partial x_j}dx_j+O(\lvert d\mathbf x\rvert^2)$
     · *why:* The relative velocity of P seen from O is the difference of the two velocities; $u_i(\mathbf x)$ cancels.
     · *plain:* How fast the neighbour moves away from, toward or round O.
  3. *did:* Drop the second-order remainder · *tex:* $du_i=(\partial u_i/\partial x_j)\,dx_j$ · *why:* The remainder
     shrinks like $\lvert d\mathbf x\rvert^2$ while the kept term shrinks like $\lvert d\mathbf x\rvert$: their ratio → 0 as P
     approaches O (orders of smallness). This is (3.10). · *plain:* For close neighbours the linear term is all that
     matters. · *set:* ring radius 0.05 m.
  4. *did:* Write it as a matrix product · *tex:* $d\mathbf u=\mathbf G\cdot d\mathbf x,\qquad G_{ij}=\partial u_i/\partial x_j$ ·
     *why:* Row i of G times the column dx is the sum over j (Ch. 2 convention: row = velocity component). One matrix
     describes every neighbour. · *plain:* Nine numbers tell how the whole neighbourhood moves.
- **Result.** $du_i=(\partial u_i/\partial x_j)\,dx_j$ (3.10), $d\mathbf u=\mathbf G\cdot d\mathbf x$ — *in words:* near any
  point, relative velocity is a linear function of separation.
- **Check.** Units: (1/s)(m) = m/s ✓. Linear field u = G·x: exact for any dx (no remainder) ✓. Numbers: G =
  [[1, 2], [0, −1]], dx = (0.01, 0.02) m → du = (0.05, −0.02) m/s; the nonlinear field of the worked example differs by
  (1, 2) × 10⁻⁴ m/s and the difference falls with slope 2 ✓.
- **What it means.** A small blob of fluid moves as a linear map: all of its stretching, shearing, swelling and spinning
  is in G. The rest of §3.4 reads G's pieces. It fails for blobs as large as the scale on which G itself varies.
- **Traps.** Thinking the remainder is zero (it is only small). Transposing G ($G_{ij}=\partial u_i/\partial x_j$: row =
  component, column = derivative direction).

### D08 · Linear strain rate: stretching per length = ∂u₁/∂x₁, and n·S·n in any direction — ★, 6 steps, in C07 (notebook · `fluid_element_deformation`)
- **Goal.** Find how fast a short material thread stretches, per unit of its own length — first along x₁ (the book's
  Fig. 3.10), then along any direction n (which the book only mentions).
- **Start.** A material segment AB along x₁ of length δx₁; by D07 its ends move at $u_1$ (A) and
  $u_1+(\partial u_1/\partial x_1)\,\delta x_1$ (B) — *in words:* the far end moves faster by the slope times the length.
- **Plan.** (1) Move both ends for a short time dt. (2) Measure the new length. (3) Divide by the old length and dt and
  take the limit. (4) Generalise to any axis and any direction.
- **Tools.** Material line element (primer in C07) · D07, $du_i=(\partial u_i/\partial x_j)dx_j$ · limit of a difference
  quotient and orders of smallness (ch02 P68) · the strain-rate tensor (3.12) (R02).
- **Assumptions.** dt small (terms of order dt² dropped, step 2); δx₁ small (D07 applies).
- **Steps.**
  1. *did:* Move both ends for dt · *tex:* $AA'=u_1\,dt,\qquad BB'=\Big(u_1+\dfrac{\partial u_1}{\partial x_1}\delta x_1\Big)dt$ ·
     *why:* Each end moves with its own velocity; over a short dt the displacement is velocity × dt (first order in dt).
     · *plain:* B travels a little further than A.
  2. *did:* Book-keep the new length · *tex:* $A'B'=AB+BB'-AA'=\delta x_1+\dfrac{\partial u_1}{\partial x_1}\delta x_1\,dt$ ·
     *why:* Along x₁ the new length is the old one plus how much further B went than A (Fig. 3.10). Sideways velocity
     differences only tilt AB — a dt² change (the book skips this). · *plain:* The thread grows by
     the extra distance of its far end.
  3. *did:* Divide by length and time · *tex:* $\dfrac{A'B'-AB}{AB\,dt}=\dfrac{\partial u_1}{\partial x_1}$ · *why:* The
     fractional growth per unit time is what we want (a rate independent of the thread's length); δx₁ and dt cancel
     exactly. · *plain:* The thread stretches by ∂u₁/∂x₁ of its length per second. · *live:* "(1.0202 − 1)/(1 × 0.01) =
     2.02 → 2.00 as dt → 0".
  4. *did:* Take the limit dt → 0 · *tex:* $\dfrac{1}{\delta x_1}\dfrac{D}{Dt}(\delta x_1)=\dfrac{\partial u_1}{\partial x_1}$ ·
     *why:* The difference quotient becomes a derivative following the material thread — hence D/Dt; the neglected
     order-dt² pieces vanish after dividing by dt (P68). · *plain:* The book's linear strain rate.
  5. *did:* Same construction on each axis · *tex:* $\dfrac{1}{\delta x_\eta}\dfrac{D(\delta x_\eta)}{Dt}=\dfrac{\partial u_\eta}{\partial x_\eta}=S_{\eta\eta}$ ·
     *why:* Nothing in steps 1–4 used the choice of axis; a Greek index means no summation. The diagonal of G equals the
     diagonal of S, since $S_{\eta\eta}=\tfrac12(G_{\eta\eta}+G_{\eta\eta})$ by (3.12). · *plain:* The diagonal of S is the
     stretching rate along each axis.
  6. *did:* Any direction n · *tex:* $\dfrac{1}{\ell}\dfrac{D\ell}{Dt}=\mathbf n\cdot\mathbf G\cdot\mathbf n=\mathbf n\cdot\mathbf S\cdot\mathbf n$ ·
     *why:* For $\delta\mathbf x=\ell\mathbf n$: $\ell\,D\ell/Dt=\delta\mathbf x\cdot\mathbf G\cdot\delta\mathbf x$ (differentiate
     ℓ², D07); the antisymmetric part drops because $\mathbf n\cdot\mathbf A\cdot\mathbf n=0$. The book skips this. · *plain:*
     Stretching along any direction is a quadratic form of S. · *set:* θ = 45°.
- **Result.** $\frac{1}{\delta x_1}\frac{D}{Dt}(\delta x_1)=\frac{\partial u_1}{\partial x_1}$; in general
  $\frac1\ell\frac{D\ell}{Dt}=\mathbf n\cdot\mathbf S\cdot\mathbf n$ — *in words:* the diagonal entries of S (in any
  axes) are stretching rates per unit length.
- **Check.** Units 1/s ✓. Rigid rotation: n·S·n = 0 for every n ✓. Numbers: u = (2x, −2y): S₁₁ = 2 s⁻¹, a 1 cm thread
  → 1.0202 cm after 0.01 s (e^{0.02}) ✓; at 45°: ½(2 − 2) = 0 ✓.
- **What it means.** Stretching per length is a property of the flow at a point, not of the thread; its extreme values
  over all directions are the principal strain rates (C12). Vortex stretching (Ch. 5) and tracer filamentation (Ch. 12)
  are this rate at work.
- **Traps.** Using AB + BB′ + AA′ (A moves too, so subtract AA′). Believing the sideways shear stretches the thread at
  first order. Keeping R in n·G·n (it contributes nothing).

### D09 · Shear strain rate: S₁₂ = ½ D(α + β)/Dt — ★★, 8 steps, in C08 (notebook · `fluid_element_deformation`)
- **Goal.** Find how fast two initially perpendicular material threads close their right angle, and show that half of
  that rate is the off-diagonal entry S₁₂ (Fig. 3.11).
- **Start.** Two threads from the corner B: one along x₂ of length δx₂ (up to C), one along x₁ of length δx₁; dα is the
  clockwise tilt of the vertical thread, dβ the counterclockwise tilt of the horizontal one — *in words:* both tilts
  close the right angle.
- **Plan.** (1) Find each thread's tilt in dt from the relative velocity of its far end. (2) Use the small-angle
  approximation. (3) Add the tilts, halve, divide by dt. (4) Recognise S₁₂; generalise to any perpendicular pair.
- **Tools.** D07 · small-angle approximation (primer in C08) · (3.12), $S_{ij}=\tfrac12(\partial u_i/\partial x_j+\partial u_j/\partial x_i)$
  (R02) · tensor components in rotated axes (Ch. 2 §2.4).
- **Assumptions.** dt small (angles ∝ dt, step 3); threads short (D07).
- **Steps.**
  1. *did:* Relative velocity of C (top) · *tex:* $du_1=\dfrac{\partial u_1}{\partial x_2}\,\delta x_2$ · *why:* D07 with
     $d\mathbf x=(0,\delta x_2)$: the top of the vertical thread moves sideways (along x₁) faster than the bottom by this
     much. · *plain:* The top slides ahead of the bottom.
  2. *did:* Its sideways offset over its height · *tex:* $\tan d\alpha=\dfrac{(\partial u_1/\partial x_2)\,\delta x_2\,dt}{\delta x_2}$ ·
     *why:* In dt the top gains the offset du₁·dt; a right triangle with that offset over the height δx₂ gives the tilt
     angle. · *plain:* Tilt = sideways shift ÷ height.
  3. *did:* Small-angle approximation · *tex:* $d\alpha=\dfrac{\partial u_1}{\partial x_2}\,dt$ · *why:* tan dα ≈ dα with an
     error of order dt³ (primer); δx₂ cancels. Stretching of the thread changes its height only at order dt² — the book
     skips this. · *plain:* The vertical thread tilts clockwise at ∂u₁/∂x₂.
  4. *did:* Same for the horizontal thread · *tex:* $d\beta=\dfrac{\partial u_2}{\partial x_1}\,dt$ · *why:* D07 with
     $d\mathbf x=(\delta x_1,0)$: its right end rises by $(\partial u_2/\partial x_1)\delta x_1 dt$ over the run δx₁; small
     angle again. · *plain:* The horizontal thread tilts counterclockwise at ∂u₂/∂x₁.
  5. *did:* Add the two tilts · *tex:* $d\alpha+d\beta=\Big(\dfrac{\partial u_1}{\partial x_2}+\dfrac{\partial u_2}{\partial x_1}\Big)dt$ ·
     *why:* α is clockwise from the vertical and β counterclockwise from the horizontal: both move the threads toward
     each other, so their sum is how much the right angle closes. · *plain:* The corner shrinks by α + β.
  6. *did:* Halve and divide by dt · *tex:* $\tfrac12\dfrac{D(\alpha+\beta)}{Dt}=\tfrac12\Big(\dfrac{\partial u_1}{\partial x_2}+\dfrac{\partial u_2}{\partial x_1}\Big)$ ·
     *why:* The book defines the shear rate as the *average* closing rate of the two threads (hence ½); D/Dt because the
     threads are material; the limit is exact. · *plain:* The average rate at which each thread turns toward the other.
     · *live:* "½(1.00 + 0.00) = 0.50 s⁻¹".
  7. *did:* Recognise the strain-rate entry · *tex:* $\tfrac12\dfrac{D(\alpha+\beta)}{Dt}=S_{12}=S_{21}$ · *why:* (3.12) with
     i = 1, j = 2 is exactly the bracket of step 6; S is symmetric. · *plain:* S₁₂ is half the closing rate of a right
     angle.
  8. *did:* Any perpendicular pair · *tex:* $\tfrac12\dfrac{D(\alpha+\beta)}{Dt}=\mathbf n_1\cdot\mathbf S\cdot\mathbf n_2$ ·
     *why:* Repeat steps 1–7 in axes along $\mathbf n_1$, $\mathbf n_2$: the off-diagonal component in rotated axes is
     $\mathbf n_1\cdot\mathbf S\cdot\mathbf n_2$ (tensor rule, Ch. 2 §2.4). The book skips this. · *plain:* Every pair of
     directions has its own closing rate. · *set:* θ = 30°.
- **Result.** $\tfrac12\frac{D(\alpha+\beta)}{Dt}=\tfrac12\Big(\frac{\partial u_1}{\partial x_2}+\frac{\partial u_2}{\partial x_1}\Big)=S_{12}=S_{21}$ —
  *in words:* the off-diagonal strain rate is half the rate at which a right angle closes.
- **Check.** Units rad/s = 1/s ✓. Solid-body rotation: dα = −dβ (both threads turn the same way) ⇒ S₁₂ = 0 ✓. Numbers:
  shear u = (γy, 0), γ = 1 s⁻¹, dt = 0.01 s: dα = 0.01 rad, dβ = 0 ⇒ S₁₂ = 0.5 s⁻¹ = γ/2 ✓.
- **What it means.** S₁₂ measures shearing — what a Newtonian fluid resists with its viscosity (Ch. 4: τ₁₂ = 2μS₁₂,
  which reduces to Ch. 1's τ = μ du/dy for u = (u(y), 0)). Its value depends on the pair of directions; along the
  principal axes it vanishes (C12).
- **Traps.** Signs: α is clockwise, β counterclockwise — both *close* the angle. Forgetting the ½ (γ = 2S₁₂, not S₁₂; Ch.
  2's Ex. 2.4 used Γ = S₁₂). Keeping the stretching of the sides (second order).

### D10 · Rigid motion does not deform: U + Ω × x ⇒ S = 0, ω = 2Ω — ★, 5 steps, in C08 (notebook)
- **Goal.** Show that a rigid motion — any translation plus any rotation — has zero strain rate, so that S is the same
  for every translating or rotating observer; the book leaves it to Exercise 3.17.
- **Start.** $\mathbf u=\mathbf U+\boldsymbol\Omega\times\mathbf x$, U and Ω uniform in space — *in words:* every point moves
  like a point of a rigid body.
- **Plan.** (1) Write the velocity in index form. (2) Differentiate to get G. (3) Show its symmetric part vanishes. (4)
  Read off R and the vorticity. (5) Conclude what a change of frame does to S.
- **Tools.** Rigid-body velocity Ω × x (primer in C08) · $(\mathbf a\times\mathbf b)_i=\varepsilon_{ijk}a_jb_k$ and the
  antisymmetry of ε (Ch. 2 §2.7) · $\partial x_k/\partial x_m=\delta_{km}$ · (3.12), (3.13) (R02, R03) · (3.15)
  $R_{ij}=-\varepsilon_{ijk}\omega_k$ (from Ch. 2; recapped as R05 before C10).
- **Assumptions.** U and Ω do not depend on position (step 2); they may depend on time — only spatial derivatives
  enter.
- **Steps.**
  1. *did:* Write the velocity in index form · *tex:* $u_i=U_i+\varepsilon_{ijk}\Omega_jx_k$ · *why:* The cross product in
     components (Ch. 2 §2.7): $(\boldsymbol\Omega\times\mathbf x)_i=\varepsilon_{ijk}\Omega_jx_k$. Index form lets us
     differentiate term by term. · *plain:* A translation plus a rotation, written with ε.
  2. *did:* Differentiate with respect to x_m · *tex:* $\dfrac{\partial u_i}{\partial x_m}=\varepsilon_{ijm}\Omega_j$ · *why:*
     U and Ω are uniform, so only $x_k$ varies and $\partial x_k/\partial x_m=\delta_{km}$ picks k = m. · *plain:* The
     velocity gradient of a rigid motion.
  3. *did:* Symmetrise · *tex:* $S_{im}=\tfrac12(\varepsilon_{ijm}+\varepsilon_{mji})\Omega_j=0$ · *why:* (3.12) adds G and its
     transpose; swapping the first and last index of ε flips its sign, $\varepsilon_{mji}=-\varepsilon_{ijm}$, so the bracket
     vanishes. · *plain:* A rigid motion stretches and shears nothing.
  4. *did:* Read off R and ω · *tex:* $R_{im}=2\varepsilon_{ijm}\Omega_j=-2\varepsilon_{imj}\Omega_j\ \Rightarrow\ \boldsymbol\omega=2\boldsymbol\Omega$ ·
     *why:* (3.13): R = G − Gᵀ = 2G here; one index swap, compared with (3.15), $R_{ij}=-\varepsilon_{ijk}\omega_k$, gives the
     vorticity. · *plain:* The vorticity of a rigid rotation is twice its angular velocity.
  5. *did:* Add a rigid motion to any flow · *tex:* $\mathbf G_{\rm new}=\mathbf G+\mathbf G_{\rm rigid}\ \Rightarrow\ \mathbf S_{\rm new}=\mathbf S$ ·
     *why:* Changing to a translating or rotating observer adds (or subtracts) a rigid motion; G is linear in u and the
     rigid part has S = 0 (step 3), so only R changes. · *plain:* S is observer-independent; ω is not.
- **Result.** $\mathbf u=\mathbf U+\boldsymbol\Omega\times\mathbf x\ \Rightarrow\ S_{ij}=0,\ \boldsymbol\omega=2\boldsymbol\Omega$ —
  *in words:* rigid motion does not deform, and its vorticity is twice its rotation rate.
- **Check.** Units 1/s ✓. Numbers: Ω = (0, 0, 1) s⁻¹: G = [[0, −1], [1, 0]], S = 0, ω₃ = 1 − (−1) = 2 ✓; the random
  (U, Ω) of the notebook: S = 0 to 1e-12, ω = 2Ω = (1, −2, 4) ✓.
- **What it means.** A fluid in rigid rotation (a stirred cup after it settles, a spun-up tank) has no internal friction:
  the viscous stress can depend only on S (Ch. 4 §4.5). And "S is frame-independent" means: every observer who moves
  rigidly agrees on how fluid elements deform — while they disagree on how they spin (N28).
- **Traps.** Letting Ω depend on position (then it is no longer rigid). Concluding "ω is frame-independent too" — step 5
  changes R. Losing the factor 2 (ω = 2Ω, spin = Ω).

### D11 · Volumetric strain rate: (1/δV) D(δV)/Dt = S_ii, Eq. (3.14) — ★★, 7 steps, in C09 (notebook · `fluid_element_deformation`)
- **Goal.** Show that the fractional rate at which a small volume of fluid grows is the divergence of the velocity —
  the book states it and leaves the proof to Exercise 3.18.
- **Start.** $\delta V=\delta x_1\,\delta x_2\,\delta x_3$ — *in words:* a small box with material edges along the axes.
- **Plan.** (1) Differentiate the product following the fluid. (2) Use the linear strain rate of each edge (D08). (3)
  Divide by δV. (4) Explain why shear does not count and why the answer does not depend on orientation. (5) The exact
  finite-time factor.
- **Tools.** Product rule for three factors (gloss in step 1; ch01 P38) · D08,
  $\frac{1}{\delta x_\eta}\frac{D(\delta x_\eta)}{Dt}=\frac{\partial u_\eta}{\partial x_\eta}$ · trace invariance (Ch. 2 §2.5) ·
  Jacobi's formula (gloss in step 7).
- **Assumptions.** The box is small (D08 applies to each edge); first order in dt (step 5).
- **Steps.**
  1. *did:* Product rule for three factors · *tex:* $\dfrac{D(\delta V)}{Dt}=\dfrac{D\delta x_1}{Dt}\delta x_2\delta x_3+\delta x_1\dfrac{D\delta x_2}{Dt}\delta x_3+\delta x_1\delta x_2\dfrac{D\delta x_3}{Dt}$
     · *why:* Differentiate a product of three factors: one term per factor, the others untouched (the two-factor rule
     of Ch. 1 applied twice). · *plain:* The volume grows because each edge grows.
  2. *did:* Insert each edge's stretching rate · *tex:* $\dfrac{D\delta x_\eta}{Dt}=\dfrac{\partial u_\eta}{\partial x_\eta}\,\delta x_\eta\quad(\text{no sum})$ ·
     *why:* D08 step 5: each material edge stretches at its own linear strain rate; η is a Greek index, so no summation.
     · *plain:* Edge η grows at ∂u_η/∂x_η of its length.
  3. *did:* Substitute into the product rule · *tex:* $\dfrac{D(\delta V)}{Dt}=\Big(\dfrac{\partial u_1}{\partial x_1}+\dfrac{\partial u_2}{\partial x_2}+\dfrac{\partial u_3}{\partial x_3}\Big)\delta x_1\delta x_2\delta x_3$
     · *why:* Each term of step 1 becomes one rate times the same product δx₁δx₂δx₃; factor it out. · *plain:* Three
     stretching rates add up.
  4. *did:* Divide by δV · *tex:* $\dfrac{1}{\delta V}\dfrac{D(\delta V)}{Dt}=\dfrac{\partial u_i}{\partial x_i}=S_{ii}$ · *why:*
     δV ≠ 0; the sum is the summation convention's $\partial u_i/\partial x_i=\nabla\cdot\mathbf u$, equal to the trace of S
     since $S_{ii}=G_{ii}$. This is (3.14). · *plain:* The volume grows at ∇·u of itself per second. · *live:* "tr G = 0.5 +
     0.5 = 1.0 s⁻¹".
  5. *did:* Shear changes volume only at second order · *tex:* $\delta V(t+dt)=\delta V\,(1+S_{ii}\,dt)+O(dt^2)$ · *why:* The
     book skips this: shearing tilts the faces by angles ∝ dt; a parallelogram's area changes as cos of the tilt, 1 −
     O(dt²). Only edge stretching acts at first order. · *plain:* Tilting a box does not change its volume to first order.
  6. *did:* Independence of the box's orientation · *tex:* $S'_{ii}=S_{ii}$ · *why:* The trace is an invariant under a
     rotation of axes (Ch. 2 §2.5); a box turned any way gives the same rate. · *plain:* Any small blob, any orientation,
     same fractional growth rate.
  7. *did:* The exact factor for a linear flow · *tex:* $\dfrac{\delta V(t)}{\delta V(0)}=\det e^{\mathbf Gt}=e^{t\,\mathrm{tr}\,\mathbf G}$ ·
     *why:* A linear flow maps positions by $e^{\mathbf Gt}$ (ch02 P79) and volumes scale by its determinant; Jacobi's
     formula det e^M = e^{tr M} (gloss, checked numerically). Step 5 is its first-order tangent. · *plain:* Over finite
     times the growth compounds exponentially. · *live:* "e^{1.0 × 1} = 2.718".
- **Result.** $\frac{1}{\delta V}\frac{D}{Dt}(\delta V)=\frac{\partial u_1}{\partial x_1}+\frac{\partial u_2}{\partial x_2}+\frac{\partial u_3}{\partial x_3}=\frac{\partial u_i}{\partial x_i}=S_{ii}$
  (3.14) — *in words:* the divergence is the fractional growth rate of a small blob of fluid.
- **Check.** Units 1/s ✓. Simple shear: ∇·u = 0, volume kept ✓. Numbers: u = (x, y, z): 3 s⁻¹; a 1 cm³ box after 0.01 s
  → 1.03045 cm³ = e^{0.03} (first order 1.03) ✓; the random G of the notebook: det e^{Gt} = e^{t tr G} to 1e-12 ✓.
- **What it means.** Incompressible flow (∇·u = 0) keeps every blob's volume; compressible flow changes it. With mass
  conservation this becomes the continuity equation $D\rho/Dt=-\rho\nabla\cdot\mathbf u$ (Ch. 4). D23 derives the same
  result from the Reynolds transport theorem.
- **Traps.** Summing over the Greek index in step 2. Believing shear changes volume at first order. Using 1 + t tr G for
  finite times (the exact factor is e^{t tr G}).

### D12 · The spin of a fluid element is half the vorticity: ½ D(−α + β)/Dt = ½ω₃ — ★★, 8 steps, in C10 (notebook · `spin_and_principal_axes`)
- **Goal.** Define how fast a small fluid element rotates and show that the rate is half the vorticity component
  ω₃ = R₂₁ (Fig. 3.11).
- **Start.** The two perpendicular threads of D09: in dt the vertical thread tilts **clockwise** by
  $d\alpha=(\partial u_1/\partial x_2)dt$ and the horizontal thread **counterclockwise** by $d\beta=(\partial u_2/\partial x_1)dt$
  — *in words:* in general the two threads turn by different amounts.
- **Plan.** (1) Measure both turns counterclockwise-positive. (2) Define the element's rotation as their average (and
  say why). (3) Substitute the tilts and take the limit. (4) Recognise R₂₁ and ω₃.
- **Tools.** D09 steps 3–4 · angular velocity of a line (primer in C10) · (3.13) $R_{ij}=\partial u_i/\partial x_j-\partial u_j/\partial x_i$
  (R03) · (3.15) $R_{ij}=-\varepsilon_{ijk}\omega_k$ (R05) · (3.16) $\omega_3=\partial u_2/\partial x_1-\partial u_1/\partial x_2$ (R06).
- **Assumptions.** dt small; threads short (as in D09).
- **Steps.**
  1. *did:* Signed turns, counterclockwise positive · *tex:* $d\theta_{\rm horiz}=+d\beta,\qquad d\theta_{\rm vert}=-d\alpha$ · *why:*
     Angles of lines are measured counterclockwise (primer); α was defined clockwise, so the vertical thread's turn
     enters with a minus sign. · *plain:* One thread turns one way, the other the other way (in shear).
  2. *did:* Define the spin as their average · *tex:* $\Omega_{\rm el}\equiv\tfrac12\dfrac{D(-\alpha+\beta)}{Dt}$ · *why:* The book
     skips the reason: averaging a perpendicular pair cancels the closing (shear) part, which turns the two threads
     oppositely, and keeps what turns the element as a whole. · *plain:* The element's rotation = the average turning of
     two perpendicular threads.
  3. *did:* Substitute the two tilts · *tex:* $\Omega_{\rm el}=\lim_{dt\to0}\dfrac{1}{2dt}\Big(-\dfrac{\partial u_1}{\partial x_2}dt+\dfrac{\partial u_2}{\partial x_1}dt\Big)$ ·
     *why:* D09 steps 3–4 give dα and dβ to first order in dt. · *plain:* Turn of each thread in dt, averaged.
  4. *did:* Cancel dt · *tex:* $\Omega_{\rm el}=\tfrac12\Big(-\dfrac{\partial u_1}{\partial x_2}+\dfrac{\partial u_2}{\partial x_1}\Big)$ ·
     *why:* dt cancels exactly, so the limit is immediate; the neglected pieces are order dt². · *plain:* The spin rate in
     terms of two velocity slopes. · *live:* "½(−1.0 + (−0.3)) = −0.65 rad/s".
  5. *did:* Recognise the rotation tensor · *tex:* $\Omega_{\rm el}=-\dfrac{R_{12}}{2}=\dfrac{R_{21}}{2}$ · *why:* (3.13):
     $R_{12}=\partial u_1/\partial x_2-\partial u_2/\partial x_1$, which is minus the bracket; R is antisymmetric,
     $R_{21}=-R_{12}$. · *plain:* The spin is half an entry of the book's R.
  6. *did:* Recognise the vorticity · *tex:* $\Omega_{\rm el}=\tfrac12\omega_3$ · *why:* (3.15) with i = 2, j = 1:
     $R_{21}=-\varepsilon_{213}\omega_3=\omega_3$ (ε₂₁₃ = −1). · *plain:* The element spins at half the vorticity. · *set:*
     the paddle wheel turns at this rate.
  7. *did:* Cross-check with the curl · *tex:* $\omega_3=\dfrac{\partial u_2}{\partial x_1}-\dfrac{\partial u_1}{\partial x_2}$ ·
     *why:* (3.16) gives the same combination as step 4 without the ½ — the two routes agree. · *plain:* Consistent with
     ω = ∇ × u.
  8. *did:* Repeat in the other two planes · *tex:* $\boldsymbol\Omega_{\rm el}=\tfrac12\boldsymbol\omega$ · *why:* The same
     construction with threads in the (x₂, x₃) and (x₃, x₁) planes gives ½ω₁ and ½ω₂; the three rates form a vector. ·
     *plain:* A small element spins with angular velocity ω/2.
- **Result.** $\tfrac12\frac{D(-\alpha+\beta)}{Dt}=\tfrac12\Big(-\frac{\partial u_1}{\partial x_2}+\frac{\partial u_2}{\partial x_1}\Big)=-\frac{R_{12}}2=\frac{R_{21}}2=\tfrac12\omega_3$
  — *in words:* vorticity is twice the spin of a small fluid element.
- **Check.** Units rad/s ✓. Solid-body rotation at ω₀: G = [[0, −ω₀], [ω₀, 0]]: ½(ω₀ + ω₀) = ω₀ — the element spins as
  fast as it revolves ✓ (Fig. 3.15). Shear γ: ½(−γ + 0) = −γ/2 ✓. Pure strain: 0 ✓.
- **What it means.** A paddle wheel carried by the flow turns at ½ω; the vorticity field tells where fluid spins. Because
  single threads may turn at different rates, "the" rotation needs this average — D14 shows it is the same for *every*
  perpendicular pair. The spin depends on the observer's own rotation (D13).
- **Traps.** Averaging α + β (that is the strain rate, D09) instead of −α + β. Forgetting the ½ (ω is twice the spin).
  Sign: R₁₂ vs R₂₁.

### D13 · Vorticity depends on the observer's rotation: ω′_z = ω_z − 2Ω — ★★, 6 steps, in C10 (notebook · `spin_and_principal_axes`)
- **Goal.** Show how the vorticity changes for an observer who rotates — the book says only that it does (Exercise
  3.19); the answer is the relation between relative and absolute vorticity that geophysical fluid dynamics runs on.
- **Start.** $\mathbf u=\boldsymbol\Omega\times\mathbf x+\mathbf u'$ with $\boldsymbol\Omega=\Omega\,\mathbf e_z$ constant, at the
  instant the two sets of axes coincide — *in words:* the lab velocity is the velocity the rotating observer measures
  plus the velocity of the observer's own frame at that point.
- **Plan.** (1) Solve for the observer's velocity. (2) Take the curl of both sides. (3) Use the curl of a rigid rotation.
  (4) Read off the z-component and the special frame where the element looks still.
- **Tools.** Rotating frame of reference (primer in C10) · the curl is linear and a vector (Ch. 2 §2.9) · D10,
  ∇ × (Ω × x) = 2Ω (ω of a rigid rotation) · (3.16) (R06).
- **Assumptions.** Ω constant (step 3); the comparison is made at the instant the axes coincide, so components agree
  (step 2).
- **Steps.**
  1. *did:* Solve for the rotating observer's velocity · *tex:* $\mathbf u'=\mathbf u-\boldsymbol\Omega\times\mathbf x$ · *why:* A
     point fixed on the rotating frame moves with Ω × x in the lab (rigid-body velocity, primer); the observer measures
     motion relative to that. · *plain:* Subtract the frame's own motion.
  2. *did:* Take the curl of both sides · *tex:* $\nabla\times\mathbf u'=\nabla\times\mathbf u-\nabla\times(\boldsymbol\Omega\times\mathbf x)$ ·
     *why:* The curl is linear. At the instant the axes coincide both observers use the same coordinates, and the curl is a
     vector (same components), so ∇′× = ∇× here. · *plain:* The observer's vorticity is the lab vorticity minus that of the
     frame's rotation.
  3. *did:* Curl of a rigid rotation · *tex:* $\nabla\times(\boldsymbol\Omega\times\mathbf x)=2\boldsymbol\Omega$ · *why:* D10 step 4:
     a rigid rotation at Ω has vorticity 2Ω (also Ch. 2's Ex. 2.3, ∇ × (b × x) = 2b). · *plain:* The frame itself carries
     vorticity 2Ω.
  4. *did:* Combine · *tex:* $\boldsymbol\omega'=\boldsymbol\omega-2\boldsymbol\Omega$ · *why:* Substitute step 3 into step 2 with
     ω = ∇ × u and ω′ = ∇ × u′. · *plain:* A rotating observer sees the vorticity reduced by twice its rotation rate.
  5. *did:* The z-component · *tex:* $\omega'_z=\omega_z-2\Omega$ · *why:* Ω points along z, so only the z-component
     changes. · *plain:* The book's statement, with the factor 2. · *live:* "ω′₃ = 2.0 − 2 × 1.0 = 0".
  6. *did:* The frame where the element is still · *tex:* $\omega'_z=0\iff\Omega=\tfrac12\omega_z$ · *why:* Set step 5 to zero.
     Rotating with the element's own spin (½ω_z, D12) makes it look non-rotating — the book's 'co-rotating frame'. · *plain:*
     Spin with the paddle wheel and it stops. · *set:* Ω = ω₃/2.
- **Result.** $\omega'_z=\omega_z-2\Omega$, i.e. $\boldsymbol\omega=\boldsymbol\omega'+2\boldsymbol\Omega$ — *in words:* vorticity
  measured by a rotating observer is the absolute vorticity minus twice the frame's rotation rate.
- **Check.** Units 1/s ✓. Ω = 0: no change ✓. Solid body ω₀ = 1 s⁻¹ (ω_z = 2) seen from a turntable at Ω = 1 rad/s:
  ω′ = 0 — the tank looks at rest ✓. The notebook's numeric curl of u − Ω × x equals ω − 2Ω to 1e-8 ✓.
- **What it means.** On the Earth (Ω = 7.29 × 10⁻⁵ s⁻¹) the vorticity measured relative to the ground, ζ, and the
  absolute vorticity differ by the planetary vorticity (its local vertical part is f = 2Ω sin φ ≈ 10⁻⁴ s⁻¹) — often larger
  than ζ itself. Conservation laws for vorticity (Ch. 5, potential vorticity in Ch. 13) hold for the absolute one. Unlike
  S (D10), ω is not observer-independent.
- **Traps.** Forgetting the factor 2 (the frame's *vorticity* is 2Ω). Confusing this with Galilean invariance (C05):
  a rotating frame is not Galilean. Applying ω′ = ω − 2Ω to components not along Ω.

### D14 · In a parallel shear flow every perpendicular pair of threads spins at −γ/2 — ★★, 6 steps, in C10 (notebook · `spin_and_principal_axes`)
- **Goal.** Show the book's claim (§3.5, stated without proof) that in the shear flow u = (γx₂, 0) the average turning
  rate of two perpendicular threads is −γ/2 = ω₃/2 whatever pair you pick — while single threads turn at different
  rates.
- **Start.** Locally $\mathbf u=(\gamma x_2,0)$, so $\mathbf G=\begin{bmatrix}0&\gamma\\0&0\end{bmatrix}$; a material thread at
  angle θ has unit direction $\mathbf e(\theta)=(\cos\theta,\sin\theta)$ — *in words:* layers slide over each other at shear
  rate γ.
- **Plan.** (1) Turning rate of one thread from the relative velocity of its tip. (2) Evaluate for this G. (3) Do the
  same for its perpendicular partner. (4) Average.
- **Tools.** Angular velocity of a line (primer in C10) · D07, $d\mathbf u=\mathbf G\cdot d\mathbf x$ · sin² + cos² = 1 ·
  (3.16) (R06).
- **Assumptions.** Threads short enough for the linear field (D07).
- **Steps.**
  1. *did:* Turning rate from the tip's velocity · *tex:* $\dot\theta=\mathbf e_\theta\cdot\mathbf G\cdot\mathbf e(\theta),\qquad\mathbf e_\theta=(-\sin\theta,\cos\theta)$ ·
     *why:* The tip of a unit thread moves relative to its tail by G·e (D07); only the part perpendicular to the thread
     turns it (primer: angular velocity of a line). · *plain:* How fast a thread at angle θ turns.
  2. *did:* Multiply G by e(θ) · *tex:* $\mathbf G\cdot\mathbf e(\theta)=(\gamma\sin\theta,\ 0)$ · *why:* Matrix–vector product: the
     only nonzero entry, $G_{12}=\gamma$, multiplies the second component sin θ. · *plain:* The tip is pushed along the flow
     in proportion to its height.
  3. *did:* Take the perpendicular part · *tex:* $\dot\theta=-\gamma\sin^2\theta$ · *why:* Take the dot product of step 2 with $\mathbf e_\theta=(-\sin\theta,\cos\theta)$:
     $-\sin\theta\cdot\gamma\sin\theta+\cos\theta\cdot0$; only the perpendicular part turns the thread. · *plain:* Every thread turns clockwise, fastest when vertical. · *live:*
     "−1 × sin²30° = −0.25 rad/s".
  4. *did:* Check the book's two threads · *tex:* $\dot\theta(\pi/2)=-\gamma,\qquad\dot\theta(0)=0$ · *why:* Substitute θ = 90°
     (the vertical thread AB) and θ = 0 (the horizontal thread BC) — the values the book reads off Fig. 3.14. · *plain:* AB
     turns at −γ, BC not at all.
  5. *did:* The perpendicular partner · *tex:* $\dot\theta(\theta+\pi/2)=-\gamma\cos^2\theta$ · *why:* Step 3 at θ + π/2, with
     sin(θ + π/2) = cos θ. · *plain:* The partner turns at a different rate.
  6. *did:* Average the pair · *tex:* $\tfrac12\big(-\gamma\sin^2\theta-\gamma\cos^2\theta\big)=-\tfrac{\gamma}{2}=\tfrac12\omega_3$ ·
     *why:* sin²θ + cos²θ = 1 removes θ; (3.16) gives $\omega_3=\partial u_2/\partial x_1-\partial u_1/\partial x_2=-\gamma$. · *plain:*
     Every perpendicular pair spins at −γ/2, half the vorticity. · *set:* sweep θ 0 → 180°.
- **Result.** $\dot\theta=-\gamma\sin^2\theta$ for one thread; $\tfrac12[\dot\theta(\theta)+\dot\theta(\theta+\tfrac\pi2)]=-\gamma/2=\omega_3/2$
  for every θ — *in words:* single threads disagree, every perpendicular pair agrees.
- **Check.** Units rad/s ✓. γ = 1 s⁻¹, θ = 30°: −0.25 and −0.75, average −0.5 ✓. Solid-body G = [[0, −ω₀], [ω₀, 0]] in the
  same formula gives θ̇ = ω₀ for every θ ✓.
- **What it means.** "Spin" is well defined even when threads turn at different rates, and it equals ½ω. A flow in
  straight lines (a river near its bank, wind near the ground) spins every element clockwise at γ/2 — the half of the
  shear that is rotation (C11). This is why D12's average of two particular threads is the right definition.
- **Traps.** Sign: ω₃ = −γ is clockwise for γ > 0. Averaging threads that are not perpendicular. Confusing γ with ch02's
  Γ = S₁₂ (γ = 2S₁₂).

### D15 · Relative velocity = deformation + rigid rotation, Eq. (3.19) — ★★, 6 steps, in C11 (notebook · `spin_and_principal_axes`)
- **Goal.** Split the relative velocity of a neighbour into a pure deformation and a rigid rotation, and find the
  rotation's angular velocity.
- **Start.** (3.10): $du_i=(\partial u_i/\partial x_j)\,dx_j$ — *in words:* the neighbour's relative velocity is the velocity
  gradient times the separation.
- **Plan.** (1) Split the gradient into S and ½R. (2) Replace R by the vorticity. (3) Recognise a cross product —
  carefully, one index swap flips a sign. (4) Read the rotation's angular velocity.
- **Tools.** (3.11) $\partial u_i/\partial x_j=S_{ij}+\tfrac12R_{ij}$ (R01) · (3.15) $R_{ij}=-\varepsilon_{ijk}\omega_k$ (R05) ·
  antisymmetry of ε and $(\mathbf a\times\mathbf b)_i=\varepsilon_{ijk}a_jb_k$, Eq. (2.21) (Ch. 2 §2.7) · rigid-body velocity
  Ω × x (primer in C08).
- **Assumptions.** As in D07 (small dx).
- **Steps.**
  1. *did:* Split the gradient · *tex:* $du_i=\big(S_{ij}+\tfrac12R_{ij}\big)dx_j$ · *why:* (3.11),
     $\frac{\partial u_i}{\partial x_j}=S_{ij}+\tfrac12R_{ij}$: every velocity gradient is its symmetric part plus half the
     book's rotation tensor (R01). · *plain:* Two kinds of relative motion.
  2. *did:* Replace R by the vorticity · *tex:* $du_i=\big(S_{ij}-\tfrac12\varepsilon_{ijk}\omega_k\big)dx_j$ · *why:* (3.15),
     $R_{ij}=-\varepsilon_{ijk}\omega_k$ (the book's text cites (3.14), a slip — it is (3.15)). · *plain:* The antisymmetric
     part is carried by the three numbers of ω. · *live:* "ω = (0, 0, −1)".
  3. *did:* Distribute over dx_j · *tex:* $du_i=S_{ij}dx_j-\tfrac12\varepsilon_{ijk}\omega_kdx_j$ · *why:* Multiply out the
     bracket; the two terms can now be read separately. · *plain:* Deformation part and rotation part.
  4. *did:* Swap two indices of ε · *tex:* $-\tfrac12\varepsilon_{ijk}\omega_kdx_j=\tfrac12\varepsilon_{ikj}\omega_kdx_j$ · *why:* ε changes
     sign when two indices are swapped (Ch. 2 §2.7); swapping j and k absorbs the minus sign. The book skips this — the
     step that fixes the sign. · *plain:* Rewrite the rotation term with its indices in cross-product order. · *set:*
     flips the orange arrow on screen, then back.
  5. *did:* Recognise the cross product · *tex:* $\tfrac12\varepsilon_{ikj}\omega_kdx_j=\tfrac12(\boldsymbol\omega\times d\mathbf x)_i$ · *why:*
     (2.21), $(\mathbf a\times\mathbf b)_i=\varepsilon_{ijk}a_jb_k$, with a = ω carrying the second index (k) and b = dx the
     third (j). · *plain:* The rotation part is ½ω × dx.
  6. *did:* Assemble (3.19) · *tex:* $du_i=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i$ · *why:* Put steps 3–5
     together. The second term has the form of the rigid-body velocity Ω × x with $\boldsymbol\Omega=\boldsymbol\omega/2$
     (primer). · *plain:* Near any point the fluid deforms by S and turns rigidly at ω/2. · *live:* "(0.05, 0) + (0.05, 0) =
     (0.10, 0) m/s".
- **Result.** $du_i=\big(S_{ij}-\tfrac12\varepsilon_{ijk}\omega_k\big)dx_j=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i$ (3.19)
  — *in words:* relative velocity = pure deformation + rigid rotation at half the vorticity.
- **Check.** Units m/s ✓. The rotation part is perpendicular to dx (½ω × dx · dx = 0): it changes no distance ✓. Numbers
  (shear γ = 1, dx = (0, 1)): S·dx = (0.5, 0), ½ω × dx = ½(0, 0, −1) × (0, 1, 0) = (0.5, 0, 0), sum (1, 0) = G·dx ✓.
- **What it means.** Only the S part changes distances between particles, so only S can produce internal friction —
  the reason the Newtonian stress law (Ch. 4 §4.5) is built on S. The rotation part is invisible to a co-rotating
  observer (D13).
- **Traps.** Missing the sign flip in step 4 (getting −½ω × dx). Writing ω instead of ω/2 as the angular velocity.
  Citing (3.14) instead of (3.15).

### D16 · A small sphere becomes an ellipsoid on the principal axes of S — ★★, 8 steps, in C12 (notebook · `spin_and_principal_axes`)
- **Goal.** Show the book's closing claim of §3.4 (stated without proof): in a short time a small sphere of fluid
  becomes an ellipsoid whose axes are the principal axes of the strain-rate tensor.
- **Start.** (3.19): $du_i=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i$ for points dx on a sphere
  $\lvert d\mathbf x\rvert=\varepsilon$ — *in words:* each surface point moves by a strain part and a rigid rotation.
- **Plan.** (1) Set the rotation aside (it keeps shapes). (2) Go to the principal frame where S is diagonal: (3.20),
  (3.21). (3) Move each point for dt. (4) Substitute into the sphere's equation and recognise an ellipsoid. (5) Say what
  happens for finite times.
- **Tools.** D15 · rigid-body velocity (primer in C08) · eigenvalues and eigenvectors of a symmetric tensor (ch02 P80) ·
  tensor components in rotated axes (Ch. 2 §2.4) · linear map of a circle is an ellipse (primer in C12).
- **Assumptions.** First order in dt (step 4); the sphere is small (D07).
- **Steps.**
  1. *did:* Set the rigid rotation aside · *tex:* $d\mathbf u_{\rm rot}=\tfrac12\boldsymbol\omega\times d\mathbf x\ \Rightarrow\ \text{shape unchanged}$ ·
     *why:* A rigid rotation (primer) keeps all distances — it turns the finished shape by ½ω dt but cannot change a
     sphere into anything else. So only S shapes it. · *plain:* The spin turns the blob; the strain reshapes it.
  2. *did:* Rotate to the principal axes · *tex:* $d\bar{\mathbf u}=\bar{\mathbf S}\cdot d\bar{\mathbf x},\qquad\bar{\mathbf S}=\mathrm{diag}(\bar S_{11},\bar S_{22},\bar S_{33})$ ·
     *why:* S is symmetric, so it has three perpendicular eigenvectors (ch02 P80); in axes along them its components are
     diagonal (tensor rule, Ch. 2 §2.4). This is (3.20). · *plain:* In the right axes the strain is three plain
     stretchings.
  3. *did:* Write the three components · *tex:* $d\bar u_\alpha=\bar S_{\alpha\alpha}\,d\bar x_\alpha\quad(\text{no sum})$ · *why:*
     A diagonal matrix times a vector multiplies each component by its own diagonal entry. This is (3.21). · *plain:*
     Each principal coordinate grows in proportion to itself.
  4. *did:* Move each point for dt · *tex:* $d\bar x_\alpha(t+dt)=d\bar x_\alpha\,(1+\bar S_{\alpha\alpha}\,dt)$ · *why:* New position =
     old + velocity × dt (first order in dt), with step 3's velocity. · *plain:* Along the stretching axis points move out,
     along the compressing axis in.
  5. *did:* The sphere's equation before moving · *tex:* $\textstyle\sum_\alpha(d\bar x_\alpha)^2=\varepsilon^2$ · *why:* The points
     start on a sphere of radius ε; in any orthonormal axes the squared distance is the sum of squares. · *plain:* The
     starting shape.
  6. *did:* Substitute the moved positions · *tex:* $\textstyle\sum_\alpha\Big[\dfrac{d\bar x_\alpha(t+dt)}{1+\bar S_{\alpha\alpha}dt}\Big]^2=\varepsilon^2$ ·
     *why:* Invert step 4 for each coordinate (1 + S̄dt ≠ 0 for small dt) and put the old coordinates into step 5. · *plain:*
     The equation obeyed by the moved points.
  7. *did:* Recognise an ellipsoid · *tex:* $a_\alpha=\varepsilon\,(1+\bar S_{\alpha\alpha}\,dt)$ · *why:* $\sum_\alpha X_\alpha^2/a_\alpha^2=1$
     is an ellipsoid with semi-axes $a_\alpha$ along the coordinate axes — here the principal axes of S (gloss). · *plain:*
     Longest along the largest eigenvalue, shortest along the smallest. · *live:* "a = 1 + 0.5 × 0.1 = 1.05, b = 0.95 mm".
  8. *did:* Say how long it holds · *tex:* $\text{finite }t:\ \text{axes}=\text{singular vectors of }e^{\mathbf Gt}$ · *why:* The book
     skips this: over finite times rotation and strain act together; the exact shape comes from $e^{\mathbf Gt}$ (primer),
     whose axes drift from S's when ω ≠ 0 (shear: 43.6° at t = 0.1 s). · *plain:* The principal-axis picture is exact only
     at the first instant. · *set:* t = 2 s.
- **Result.** In a short time dt a sphere of radius ε becomes an ellipsoid with semi-axes $\varepsilon(1+\bar S_{\alpha\alpha}dt)$
  along the principal axes of S, turned by ½ω dt — *in words:* S decides the shape, ω/2 its orientation.
- **Check.** Units: S̄dt dimensionless ✓. Volume ratio $\prod(1+\bar S_{\alpha\alpha}dt)\approx1+S_{ii}dt$ = (3.14) ✓. Numbers
  (shear γ = 1, ε = 1 mm, dt = 0.1 s): 1.05 and 0.95 mm at ±45° (first order); strain alone e^{±0.05} = 1.0513, 0.9512;
  exact 1.0513, 0.9513 at 43.6° ✓.
- **What it means.** Blobs of dye, plankton patches and temperature anomalies are pulled out along the stretching axis
  — filaments in the ocean, fronts in the atmosphere (frontogenesis, Ch. 13). In pure strain the axes stay put; in shear
  the rotation keeps turning the ellipse toward the flow direction.
- **Traps.** Expecting the ellipse to stay at 45° in a shear flow for long times. Forgetting the rotation part turns the
  shape. Reading S̄ as "the S of the original axes".

### D17 · Vorticity in polar coordinates, Eq. (3.23), from the circulation round a small sector — ★★, 10 steps, in C13 (notebook · `vortex_paddle_wheels`)
- **Goal.** Derive the vorticity of a plane flow in polar coordinates, which the book takes from Appendix B (and leaves
  to Exercise 3.20), and apply it to solid-body rotation.
- **Start.** Vorticity is circulation per unit area (R07, from (3.18) $\Gamma=\oint_C\mathbf u\cdot d\mathbf s=\int_A\boldsymbol\omega\cdot\mathbf n\,dA$):
  $\omega_z=\lim_{A\to0}\frac{1}{A}\oint_C\mathbf u\cdot d\mathbf s$, taken round the small polar sector r…r + dr, θ…θ + dθ,
  counterclockwise — *in words:* add u·ds round a tiny sector, divide by its area.
- **Plan.** (1) The sector's area. (2) The two arcs (their lengths differ). (3) The two radial legs. (4) Add, divide,
  take the limit. (5) Apply to (3.22).
- **Tools.** Circulation and Stokes, (3.18) (R07) · polar coordinates as a moving basis (primer in C13) · line integral
  round a loop (ch02 P86) · first-order Taylor expansion (ch01 P26).
- **Assumptions.** u smooth inside the sector (step 9: the sector must not contain a singular point such as the axis of
  a line vortex).
- **Steps.**
  1. *did:* The sector's area · *tex:* $A=r\,dr\,d\theta$ · *why:* Polar area element (primer): an arc of length r dθ
     times a width dr, to leading order. · *plain:* The small area we divide by.
  2. *did:* The outer arc, counterclockwise · *tex:* $\int_{\rm outer}\mathbf u\cdot d\mathbf s=\big[u_\theta\,r\big]_{r+dr}\,d\theta$ ·
     *why:* On the arc at radius r + dr, $d\mathbf s=(r+dr)\,d\theta\,\mathbf e_\theta$ (primer), so u·ds = u_θ(r + dr)dθ, taken
     at the arc's midpoint. · *plain:* The outer arc collects the tangential velocity times its length.
  3. *did:* The inner arc, run backwards · *tex:* $\int_{\rm inner}\mathbf u\cdot d\mathbf s=-\big[u_\theta\,r\big]_{r}\,d\theta$ ·
     *why:* Going counterclockwise round the sector the inner arc is traversed clockwise, $d\mathbf s=-r\,d\theta\,\mathbf e_\theta$.
     · *plain:* The inner arc subtracts. · *set:* lights the inner arc.
  4. *did:* Combine the arcs by Taylor · *tex:* $\big[u_\theta r\big]_{r+dr}d\theta-\big[u_\theta r\big]_rd\theta=\dfrac{\partial(ru_\theta)}{\partial r}\,dr\,d\theta$ ·
     *why:* First-order Taylor in r of the product $ru_\theta$. The book skips this: the arcs have different lengths, which
     is why r stays inside the derivative. · *plain:* The arcs nearly cancel; what is left is how r u_θ changes outward. ·
     *live:* "[u_θ r] rises by … over dr".
  5. *did:* The radial leg at θ, outward · *tex:* $\int_{\theta}\mathbf u\cdot d\mathbf s=\big[u_r\big]_{\theta}\,dr$ · *why:* Along a ray
     $d\mathbf s=dr\,\mathbf e_r$, so only the radial component contributes. · *plain:* The outward leg collects u_r.
  6. *did:* The radial leg at θ + dθ, inward · *tex:* $\int_{\theta+d\theta}\mathbf u\cdot d\mathbf s=-\big[u_r\big]_{\theta+d\theta}\,dr$ ·
     *why:* Counterclockwise round the sector this leg runs inward, $d\mathbf s=-dr\,\mathbf e_r$. · *plain:* The inward leg
     subtracts.
  7. *did:* Combine the radial legs by Taylor · *tex:* $\big[u_r\big]_\theta dr-\big[u_r\big]_{\theta+d\theta}dr=-\dfrac{\partial u_r}{\partial\theta}\,d\theta\,dr$ ·
     *why:* First-order Taylor expansion of u_r in θ (P26), exactly as the two arcs were combined in r in step 4. ·
     *plain:* A radial velocity that changes with angle adds circulation.
  8. *did:* Add all four legs · *tex:* $d\Gamma=\Big[\dfrac{\partial(ru_\theta)}{\partial r}-\dfrac{\partial u_r}{\partial\theta}\Big]dr\,d\theta$ ·
     *why:* The circulation is the sum over the closed loop: steps 4 and 7. · *plain:* The circulation of the small sector.
  9. *did:* Divide by the area, shrink · *tex:* $\omega_z=\dfrac1r\dfrac{\partial}{\partial r}(ru_\theta)-\dfrac1r\dfrac{\partial u_r}{\partial\theta}$ ·
     *why:* Vorticity = circulation per unit area (R07): divide step 8 by $r\,dr\,d\theta$ and let dr, dθ → 0; the neglected
     Taylor terms vanish. This is (3.23). · *plain:* The polar formula. · *live:* "Γ_sector/(r dr dθ) = … vs formula …".
  10. *did:* Apply it to solid-body rotation · *tex:* $\omega_z=\dfrac1r\dfrac{d}{dr}(r\cdot\omega_0r)=2\omega_0$ · *why:* (3.22),
      $u_r=0$ and $u_\theta=\omega_0r$: the second term vanishes and $d(\omega_0r^2)/dr=2\omega_0r$. · *plain:* Every element
      of a solid-body rotation spins at ω₀ (half of 2ω₀) — as fast as it revolves. · *set:* kind = solid.
- **Result.** $\omega_z=\frac1r\frac{\partial}{\partial r}(ru_\theta)-\frac1r\frac{\partial u_r}{\partial\theta}$ (3.23), $=2\omega_0$ for
  (3.22) — *in words:* in polar coordinates the vorticity is the outward change of r u_θ minus the angular change of u_r,
  per unit r.
- **Check.** Units 1/s ✓. Line vortex (3.25) $u_\theta=B/r$: $\frac1r\frac{d}{dr}(B)=0$ ✓. Cartesian curl of
  (−ω₀y, ω₀x) = 2ω₀ ✓. Numbers: the notebook's four-leg sum for a Rankine core (Γ = 2π, σ = 2 m) gives 0.500 s⁻¹ = Γ/πσ² ✓.
  Optional `check_src`: `ch03.polar_vorticity_z_sym(u_r, u_th, r, th)` for a non-axisymmetric pair (u_r = r² sin θ,
  u_θ = r cos θ) vs the Cartesian curl of the same field written with $\mathbf e_r,\mathbf e_\theta$ → difference 0.
- **What it means.** "Going round" (u_θ ≠ 0) and "spinning" (ω_z ≠ 0) are different: what matters is how r u_θ changes
  with r. If r u_θ is constant (the line vortex) the fluid circles without spinning. Every vortex formula in Ch. 5, 6 and
  13 uses (3.23).
- **Traps.** Writing $\partial u_\theta/\partial r$ instead of $\frac1r\partial(ru_\theta)/\partial r$ (forgetting the arcs differ).
  Orientation of the legs. Applying it to a sector that contains the axis of a line vortex (u is singular there).

### D18 · The Rankine vortex is self-consistent: core vorticity, continuity at σ, peak at σ — ★, 5 steps, in C14 (notebook · `vortex_paddle_wheels`)
- **Goal.** Check that the two halves of (3.28) fit together — uniform vorticity inside, none outside, speed continuous,
  total circulation Γ — and find the peak speed; the book states all of this without working.
- **Start.** (3.28): $u_\theta(r)=\frac{\Gamma}{2\pi\sigma^2}r$ for $r\le\sigma$, $u_\theta(r)=\frac{\Gamma}{2\pi r}$ for $r>\sigma$,
  $u_r=0$ — *in words:* solid-body rotation inside a core of radius σ, a line vortex outside.
- **Plan.** (1) Vorticity in each region by (3.23). (2) Continuity of the speed at σ. (3) Total circulation. (4) Where the
  speed peaks.
- **Tools.** (3.23) (D17) · derivative of a product · (3.24)/(3.26) circulation of a centred circle (N37, N39) ·
  continuity at a join and a kink as a maximum (gloss).
- **Assumptions.** Axisymmetric (u_r = 0, no θ-dependence).
- **Steps.**
  1. *did:* Vorticity inside by (3.23) · *tex:* $\omega_z=\dfrac1r\dfrac{d}{dr}\Big(\dfrac{\Gamma r^2}{2\pi\sigma^2}\Big)=\dfrac{\Gamma}{\pi\sigma^2}$ ·
     *why:* u_r = 0 leaves the first term of (3.23); $d(r^2)/dr=2r$ and the r cancels. · *plain:* The core has uniform
     vorticity — a solid-body rotation with $\omega_0=\Gamma/2\pi\sigma^2$.
  2. *did:* Vorticity outside by (3.23) · *tex:* $\omega_z=\dfrac1r\dfrac{d}{dr}\Big(\dfrac{\Gamma}{2\pi}\Big)=0$ · *why:* Outside,
     r u_θ = Γ/2π is constant, so its derivative is zero. · *plain:* No vorticity outside the core.
  3. *did:* Compare the two speeds at σ · *tex:* $\dfrac{\Gamma}{2\pi\sigma^2}\sigma=\dfrac{\Gamma}{2\pi\sigma}$ · *why:* Evaluate each
     branch at r = σ (gloss: a join is continuous when the two sides agree there). The speed is continuous; the
     vorticity jumps from Γ/πσ² to 0. · *plain:* No jump in speed, a jump in spin. · *set:* loop at r = σ.
  4. *did:* Total circulation · *tex:* $\Gamma(r\ge\sigma)=2\pi r\,u_\theta=\Gamma=\dfrac{\Gamma}{\pi\sigma^2}\cdot\pi\sigma^2$ · *why:* The
     centred-circle circulation 2πr u_θ ((3.24)/(3.26)) outside the core; by Stokes it equals core vorticity × core area.
     · *plain:* All the circulation lives in the core.
  5. *did:* Locate the peak speed · *tex:* $u_{\max}=\dfrac{\Gamma}{2\pi\sigma}\ \text{at}\ r=\sigma$ · *why:* $du_\theta/dr>0$ inside
     (linear rise) and $<0$ outside (1/r fall): the derivative changes sign at σ — a maximum at a kink, not at a zero of
     the derivative (gloss). · *plain:* The fastest wind is at the edge of the core. · *live:* "Γ/2πσ = 6.283/(2π × 1) =
     1.000 m/s".
- **Result.** $\omega_z=\Gamma/\pi\sigma^2$ ($r\le\sigma$), 0 ($r>\sigma$); u_θ continuous at σ with $u_{\max}=\Gamma/2\pi\sigma$ —
  *in words:* a spinning core wrapped in an irrotational skirt, fastest at the core's edge.
- **Check.** Units: Γ/πσ² [1/s], Γ/2πσ [m/s] ✓. σ → 0 at fixed Γ: the core shrinks to the line vortex ✓. Numbers (Γ = 2π
  m²/s, σ = 1 m): inside u = r, ω = 2 s⁻¹; peak 1 m/s at 1 m ✓.
- **What it means.** The simplest model of a real vortex: bounded speed, finite core. The jump in vorticity is
  unphysical for a viscous fluid (viscosity smooths it — the Gaussian of D19).
- **Traps.** Setting du_θ/dr = 0 to find the maximum (there is no such point — the maximum is a kink). Thinking the
  vorticity is continuous because the speed is.

### D19 · The Gaussian vortex's speed from its vorticity, Eq. (3.29) — ★★, 7 steps, in C14 (notebook · `vortex_paddle_wheels`)
- **Goal.** Show that the speed profile of the Gaussian vortex (3.29) follows from its vorticity by Stokes' theorem —
  the book states both halves without connecting them — and find its two limits.
- **Start.** $\omega_z(r)=\frac{\Gamma}{\pi\sigma^2}\exp(-r^2/\sigma^2)$ — *in words:* the vorticity is largest on the axis and
  falls off smoothly over a core radius σ.
- **Plan.** (1) Circulation inside radius r as an area integral of vorticity. (2) Substitute to do the integral. (3)
  Circulation as u_θ times the circle's length. (4) Solve for u_θ. (5) Limits near and far.
- **Tools.** Stokes/circulation (3.18) (R07) · substitution in an integral (primer in C14) · ∫e^{−s}ds (ch01 P36) ·
  (3.24)-type circle integral (N37) · first-order Taylor of e^{−x} (ch01 P26).
- **Assumptions.** Axisymmetric flow, u_r = 0 (step 4: u_θ is the same all round the circle).
- **Steps.**
  1. *did:* Circulation inside radius r · *tex:* $\Gamma(r)=\displaystyle\int_0^r\omega_z(r')\,2\pi r'\,dr'$ · *why:* (3.18): circulation
     round the circle = vorticity flux through the disc; for axisymmetric ω split the disc into rings of area 2πr′dr′. ·
     *plain:* Add up the vorticity inside the circle. · *set:* the loop sweeps outward.
  2. *did:* Substitute s = r′²/σ² · *tex:* $\Gamma(r)=\dfrac{\Gamma}{\pi\sigma^2}\,\pi\sigma^2\displaystyle\int_0^{r^2/\sigma^2}e^{-s}\,ds$ · *why:*
     With $ds=2r'dr'/\sigma^2$, $2\pi r'dr'=\pi\sigma^2ds$ and the limits become 0 … r²/σ² (substitution primer). · *plain:* A plain
     exponential integral.
  3. *did:* Integrate the exponential · *tex:* $\Gamma(r)=\Gamma\big(1-e^{-r^2/\sigma^2}\big)$ · *why:* $\int_0^Xe^{-s}ds=1-e^{-X}$ (P36);
     the π σ² cancel. · *plain:* The fraction of the total circulation inside r. · *live:* "Γ(1) = Γ(1 − e⁻¹) = 0.632 Γ".
  4. *did:* Circulation of the circle directly · *tex:* $\Gamma(r)=\oint\mathbf u\cdot d\mathbf s=2\pi r\,u_\theta(r)$ · *why:* On a
     centred circle $d\mathbf s=r\,d\theta\,\mathbf e_\theta$ and u_θ is constant round it (axisymmetric) — as in (3.24),
     $\Gamma=\int_0^{2\pi}u_\theta r\,d\theta=2\pi ru_\theta$. · *plain:* The same circulation measured on the loop itself.
  5. *did:* Equate and solve for u_θ · *tex:* $u_\theta(r)=\dfrac{\Gamma}{2\pi r}\Big(1-e^{-r^2/\sigma^2}\Big)$ · *why:* Steps 3 and 4
     describe the same number; divide by 2πr (r > 0). This is (3.29). · *plain:* The speed profile of the Gaussian vortex.
  6. *did:* Near the axis, r ≪ σ · *tex:* $u_\theta\approx\dfrac{\Gamma}{2\pi\sigma^2}\,r$ · *why:* First-order Taylor
     $e^{-x}\approx1-x$ for small $x=r^2/\sigma^2$, so $1-e^{-x}\approx r^2/\sigma^2$. · *plain:* A solid-body core with
     $\omega_0=\Gamma/2\pi\sigma^2$ — like the Rankine core.
  7. *did:* Far away, r ≫ σ · *tex:* $u_\theta\to\dfrac{\Gamma}{2\pi r}$ · *why:* $e^{-r^2/\sigma^2}\to0$ much faster than any power
     of r. · *plain:* A line vortex with B = Γ/2π — irrotational outside the core. · *set:* the teal curve meets the grey
     ghost.
- **Result.** $u_\theta(r)=\frac{\Gamma}{2\pi r}\big(1-\exp(-r^2/\sigma^2)\big)$ (3.29), with $\Gamma(r)=\Gamma(1-e^{-r^2/\sigma^2})$
  — *in words:* the speed is the circulation enclosed divided by the circle's length.
- **Check.** Units m/s ✓. Γ(∞) = Γ ✓. Numbers (Γ = 2π, σ = 1): u_θ(1) = 0.632 m/s; u_θ(5) = 0.19999998 ≈ 1/5 ✓. (3.23)
  applied to (3.29) gives back ω_z (the notebook's sympy line) ✓.
- **What it means.** Any axisymmetric vorticity profile gives its speed the same way: $u_\theta=\Gamma(r)/2\pi r$. The
  Gaussian is the Lamb–Oseen vortex of viscous flow at one instant (σ² = 4νt, Ch. 5/8). Near r = 0 the formula subtracts
  nearly equal numbers — compute it with `-np.expm1(-x)` (primer).
- **Traps.** Forgetting the ring weight 2πr′ in step 1. Thinking the core has zero vorticity because the speed vanishes
  at the axis (the vorticity is largest there). Evaluating 1 − exp(−x) naively for tiny x.

### D20 · Where the Gaussian vortex's speed peaks: 1 + 2r²/σ² = e^{r²/σ²}, r ≈ 1.1209σ — ★★, 7 steps, in C14 (notebook · `vortex_paddle_wheels`)
- **Goal.** Find the radius of maximum speed of the Gaussian vortex — the book quotes r ≈ 1.12091σ and leaves the work
  to Exercise 3.26.
- **Start.** (3.29): $u_\theta(r)=\frac{\Gamma}{2\pi r}\big(1-e^{-r^2/\sigma^2}\big)$ — *in words:* zero on the axis and far
  away, so a maximum lies between.
- **Plan.** (1) Write u_θ with one dimensionless variable. (2) Differentiate and set to zero. (3) Simplify to the book's
  equation. (4) Solve numerically, excluding the trivial root.
- **Tools.** Product rule (ch01 P38) and chain rule · maximum of a function (gloss) · scipy.optimize.brentq (primer in
  C14) · Lambert W (gloss in the check).
- **Assumptions.** Γ, σ > 0 (steps 1, 3).
- **Steps.**
  1. *did:* Write u_θ with x = r²/σ² · *tex:* $u_\theta=\dfrac{\Gamma}{2\pi\sigma}\,f(x),\qquad f(x)=\dfrac{1-e^{-x}}{\sqrt x}$ · *why:*
     Substitute r = σ√x; the constant Γ/2πσ does not affect where the maximum is, and r ↦ x is increasing for r > 0. ·
     *plain:* One curve f(x) for every vortex; σ only rescales it.
  2. *did:* Differentiate f · *tex:* $f'(x)=\dfrac{e^{-x}}{\sqrt x}-\dfrac{1-e^{-x}}{2x^{3/2}}$ · *why:* Product rule (P38) on
     $(1-e^{-x})\,x^{-1/2}$: $(e^{-x})x^{-1/2}+(1-e^{-x})(-\tfrac12x^{-3/2})$; the peak is where this slope
     vanishes. · *plain:* The slope of the speed curve.
  3. *did:* Set f′ = 0, multiply by 2x^{3/2} · *tex:* $2x\,e^{-x}-(1-e^{-x})=0$ · *why:* A smooth maximum has zero slope
     (gloss); x > 0, so multiplying by $2x^{3/2}$ keeps the roots. · *plain:* The condition for the peak.
  4. *did:* Multiply by e^{x} · *tex:* $1+2x=e^{x}$ · *why:* $e^x\neq0$; rearrange $2x-e^x+1=0$. With x = r²/σ² this is the
     book's $1+2r^2/\sigma^2=\exp(r^2/\sigma^2)$. · *plain:* A straight line meets an exponential.
  5. *did:* Exclude the trivial root · *tex:* $x=0\ \text{solves it but is the axis}$ · *why:* 1 + 0 = e⁰, but r = 0 is where
     u_θ = 0 — a minimum. The line and the exponential cross once more for x > 0 (the exponential starts slower, then
     overtakes). · *plain:* We want the second crossing.
  6. *did:* Solve numerically in a bracket · *tex:* $x^*=1.25643\ \ (\text{brentq on }(0.5,3))$ · *why:* The equation has no
     elementary solution; $g(x)=1+2x-e^x$ is +0.351 at 0.5 and −13.1 at 3, so brentq (primer) finds the one root
     between. · *plain:* The peak is at x ≈ 1.2564. · *live:* "brentq: x* = 1.25643".
  7. *did:* Return to r and the peak speed · *tex:* $r^*=\sigma\sqrt{x^*}=1.12091\,\sigma,\qquad u_{\max}=0.6382\,\dfrac{\Gamma}{2\pi\sigma}$ · *why:*
     Undo x = r²/σ²; evaluate f(x*) = (1 − e^{−1.2564})/1.1209 = 0.6382. · *plain:* The strongest wind is about 12 % outside σ,
     36 % weaker than the Rankine peak. · *set:* peak marker.
- **Result.** $1+2\frac{r^2}{\sigma^2}=\exp\!\big(r^2/\sigma^2\big)$, $r^*\approx1.12091\,\sigma$, $u_{\max}\approx0.6382\,\Gamma/(2\pi\sigma)$ — *in
  words:* the Gaussian vortex's peak sits just outside its core radius.
- **Check.** Units: x dimensionless ✓. Residual $1+2x^*-e^{x^*}$ < 1e-14 ✓. Lambert W cross-check (gloss: W is the
  inverse of $we^w$): $x^*=-W_{-1}(-e^{-1/2}/2)-\tfrac12=1.2564312086$ ✓ (`scipy.special.lambertw(z, -1)`). Numbers (Γ = 2π,
  σ = 1): u_max = 0.638 m/s at 1.121 m ✓. A published fit of the Lamb–Oseen peak gives r²/σ² ≈ 1.256 (analysis §8) ✓.
- **What it means.** Measuring the radius of maximum wind of a vortex that looks Gaussian gives its core size σ = r*/1.121
  — how tornado and wake-vortex radars are read. For the Rankine model the same radius is σ itself.
- **Traps.** Accepting x = 0. Differentiating with respect to r in one term and x in another. Taking r* = 1.2564σ (that is
  x*, not √x*).

### D21 · Leibniz's theorem, Eq. (3.30) — ★★, 7 steps, in C15 (notebook · `reynolds_transport_cv`)
- **Goal.** Prove the rule for differentiating an integral whose limits move — the book cites it without proof and
  builds the Reynolds transport theorem on it.
- **Start.** $I(t)=\int_{a(t)}^{b(t)}F(x,t)\,dx$ — *in words:* the amount of F between two moving end points.
- **Plan.** (1) Write I with an antiderivative. (2) Differentiate each end with the chain rule. (3) Turn the leftover
  time derivatives into an integral of ∂F/∂t.
- **Tools.** Fundamental theorem of calculus (ch02 P84) · chain rule (ch01 P49) · differentiation under the integral
  sign (primer in C15).
- **Assumptions.** F and ∂F/∂t continuous (steps 1, 5); a(t), b(t) differentiable (step 3).
- **Steps.**
  1. *did:* Introduce an antiderivative in x · *tex:* $\Phi(x,t)\equiv\displaystyle\int_c^xF(x',t)\,dx',\qquad\dfrac{\partial\Phi}{\partial x}=F$ ·
     *why:* With a fixed lower limit c, Φ exists for continuous F and its x-derivative is F (fundamental theorem of
     calculus, P84). · *plain:* The running total of F from a fixed point.
  2. *did:* Write I with Φ · *tex:* $I(t)=\Phi(b(t),t)-\Phi(a(t),t)$ · *why:* Fundamental theorem of calculus: an integral is
     the antiderivative at the top minus at the bottom. · *plain:* The amount between the ends = total up to b minus total
     up to a.
  3. *did:* Chain rule at the upper end · *tex:* $\dfrac{d}{dt}\Phi(b(t),t)=\dfrac{\partial\Phi}{\partial x}\Big|_b\dfrac{db}{dt}+\dfrac{\partial\Phi}{\partial t}\Big|_b$ ·
     *why:* Φ depends on t through its first argument b(t) and directly — the chain rule (P49), as in D02 step 3. · *plain:*
     The total up to b changes because b moves and because F changes. · *live:* "2.5 × 0.4 + …".
  4. *did:* Both ends, with ∂Φ/∂x = F · *tex:* $\dfrac{dI}{dt}=F(b,t)\dfrac{db}{dt}-F(a,t)\dfrac{da}{dt}+\dfrac{\partial\Phi}{\partial t}\Big|_b-\dfrac{\partial\Phi}{\partial t}\Big|_a$ ·
     *why:* Step 3 for b minus the same for a, using step 1 for ∂Φ/∂x. · *plain:* The end terms have appeared.
  5. *did:* Differentiate Φ under the integral sign · *tex:* $\dfrac{\partial\Phi}{\partial t}(x,t)=\displaystyle\int_c^x\dfrac{\partial F}{\partial t}(x',t)\,dx'$ ·
     *why:* In Φ the limits c and x are held fixed while t varies, so the derivative may pass inside (primer; F and ∂F/∂t
     continuous). The book skips this. · *plain:* Φ's time change is the sum of F's local changes.
  6. *did:* Subtract at the two ends · *tex:* $\dfrac{\partial\Phi}{\partial t}\Big|_b-\dfrac{\partial\Phi}{\partial t}\Big|_a=\displaystyle\int_a^b\dfrac{\partial F}{\partial t}\,dx$ ·
     *why:* The part from c to a is common to both and cancels. · *plain:* The change of F in place, summed between the
     ends. · *set:* shades the interior band.
  7. *did:* Assemble (3.30) · *tex:* $\dfrac{d}{dt}\displaystyle\int_{a(t)}^{b(t)}F\,dx=\int_a^b\dfrac{\partial F}{\partial t}dx+\dfrac{db}{dt}F(b,t)-\dfrac{da}{dt}F(a,t)$ ·
     *why:* Put step 6 into step 4 and reorder. · *plain:* Change in place + what the upper end sweeps in − what the lower
     end sweeps out. · *live:* "0.600 + 1.000 − (−0.300) = 1.900".
- **Result.** $\frac{d}{dt}\int_{x=a(t)}^{x=b(t)}F(x,t)\,dx=\int_a^b\frac{\partial F}{\partial t}dx+\frac{db}{dt}F(b,t)-\frac{da}{dt}F(a,t)$ (3.30)
  — *in words:* interior change plus the gain at the moving upper limit minus the loss at the moving lower limit (Fig.
  3.17's three strips).
- **Check.** Units: [F]·m/s on every term ✓. Fixed limits: only the integral of ∂F/∂t ✓. F = 1: dI/dt = ḃ − ȧ, the rate the
  length changes ✓. Numbers: F = x²t, a = t, b = t² at t = 2: 18.667 + 128 − 8 = 138.667 = d/dt[(t⁷ − t⁴)/3] = (7·64 −
  4·8)/3 ✓ (`ch03.leibniz_example(2.0)`); E7's default 0.6 + 1.0 + 0.3 = 1.9 ✓. Optional `check_src`: sympy
  `sp.diff(sp.integrate(F, (x, a, b)), t)` vs the right side for F = x²t, a = t, b = t² → 0.
- **What it means.** The 1-D Reynolds transport theorem: (3.35) with the "surface" reduced to two end points with normals
  ±1 (D22 step 12). Layer budgets (Ch. 13) and the boundary-layer momentum integral (Ch. 9) are Leibniz in action.
- **Traps.** The sign of the lower-limit term (−ȧF(a): an end moving left *adds* F). Treating ∂Φ/∂t as ∂F/∂t (it is an
  integral of it). Passing d/dt inside an integral whose limits move.

### D22 · The Reynolds transport theorem: (3.31) → (3.32) → (3.33) → (3.34) → (3.35) — ★★★, 12 steps + sympy, in C15 (notebook · `reynolds_transport_cv`)
- **Goal.** Find the rate of change of the amount of F inside a volume V*(t) whose surface moves with any velocity b —
  the tool that turns every conservation law of Ch. 4 into equations. The book gives the steps; we fill the gaps (the
  sign of the sliver, the orders of smallness, the value of F in the sliver, the 1-D reduction).
- **Start.** (3.31): $\dfrac{d}{dt}\displaystyle\int_{V^*(t)}F(\mathbf x,t)dV=\lim_{\Delta t\to0}\dfrac{1}{\Delta t}\Big\{\int_{V^*(t+\Delta t)}F(\mathbf x,t+\Delta t)dV-\int_{V^*(t)}F(\mathbf x,t)dV\Big\}$
  — *in words:* the definition of a time derivative, applied to an integral whose region moves.
- **Plan.** (1) Split the new region into the old one plus a thin sliver ΔV, and Taylor-expand F in time → (3.32). (2)
  Cancel and drop what is second order → (3.33). (3) Turn the sliver into a surface integral → (3.34). (4) Take the
  limit → (3.35); check the 1-D case.
- **Tools.** (3.31) (N49) · signed swept volume of a moving surface (primer in C15) · first-order Taylor expansion in
  time (ch01 P26; gloss in step 2) · limits and orders of smallness (ch02 P68) · mean-value theorem for integrals (ch02
  P85) · Leibniz (3.30) (D21).
- **Assumptions.** F and ∂F/∂t continuous (steps 2, 9); b continuous and A* piecewise smooth, so the sliver has thickness
  O(Δt) everywhere (steps 5, 8); Δt → 0 at the end.
- **Steps.**
  1. *did:* Split the new region · *tex:* $\displaystyle\int_{V^*(t+\Delta t)}=\int_{V^*(t)}+\int_{\Delta V},\qquad\Delta V\equiv V^*(t+\Delta t)-V^*(t)$ ·
     *why:* Integrals add over regions. ΔV is **signed**: where the surface retreats (b·n < 0) the sliver is subtracted
     (primer) — one formula covers growth and shrinking; the book leaves this implicit. · *plain:* New region = old region +
     the band the walls swept. · *set:* colours ΔV blue/rose.
  2. *did:* Taylor-expand F in time · *tex:* $F(\mathbf x,t+\Delta t)=F(\mathbf x,t)+\Delta t\,\dfrac{\partial F}{\partial t}(\mathbf x,t)+O(\Delta t^2)$ ·
     *why:* First-order Taylor expansion in t at fixed x (gloss; P26); F is continuously differentiable. We want everything
     at time t. · *plain:* F a moment later ≈ F now + its local rate × Δt.
  3. *did:* Insert both into the first integral · *tex:* $\displaystyle\int_{V^*(t+\Delta t)}\!F(t+\Delta t)\,dV\cong\int_{V^*}\!F\,dV+\int_{V^*}\!\Delta t\,\dfrac{\partial F}{\partial t}dV+\int_{\Delta V}\!F\,dV+\int_{\Delta V}\!\Delta t\,\dfrac{\partial F}{\partial t}dV$ ·
     *why:* Apply step 2 inside both pieces of step 1 and distribute; this is (3.32). · *plain:* Four terms: old amount,
     local change, amount in the sliver, local change in the sliver. · *live:* "T1 … T4 = …".
  4. *did:* Cancel against the old amount · *tex:* $\text{bracket of (3.31)}=\displaystyle\int_{V^*}\Delta t\,\dfrac{\partial F}{\partial t}dV+\int_{\Delta V}F\,dV+\int_{\Delta V}\Delta t\,\dfrac{\partial F}{\partial t}dV$ ·
     *why:* The first term of (3.32) is exactly the $\int_{V^*(t)}F\,dV$ subtracted in (3.31). · *plain:* Only the changes are
     left.
  5. *did:* Size the sliver · *tex:* $\Delta V=O(\Delta t)\ \Rightarrow\ \displaystyle\int_{\Delta V}\Delta t\,\dfrac{\partial F}{\partial t}dV=O(\Delta t^2)$ ·
     *why:* The sliver's thickness is (b·n)Δt, so its volume is O(Δt); a bounded integrand times O(Δt) volume times Δt is
     O(Δt²). The book says "second order"; this is why. · *plain:* The last term is much smaller than the others. · *set:*
     Δt = 0.01 (T4 falls 100×).
  6. *did:* Divide by Δt and drop it · *tex:* $\dfrac{d}{dt}\displaystyle\int_{V^*}F\,dV=\lim_{\Delta t\to0}\dfrac{1}{\Delta t}\Big\{\int_{V^*}\Delta t\,\dfrac{\partial F}{\partial t}dV+\int_{\Delta V}F\,dV\Big\}$ ·
     *why:* After dividing by Δt the dropped term is O(Δt) → 0 (orders of smallness, P68). This is (3.33). · *plain:* Two
     contributions survive.
  7. *did:* Take the volume term's limit · *tex:* $\dfrac{1}{\Delta t}\displaystyle\int_{V^*}\Delta t\,\dfrac{\partial F}{\partial t}dV=\int_{V^*}\dfrac{\partial F}{\partial t}dV$ ·
     *why:* Δt is a constant for the integral, so it comes out and cancels exactly — no limit needed. · *plain:* The change
     of F in place, summed over the volume.
  8. *did:* The volume swept by one patch · *tex:* $dV_{\rm swept}=(\mathbf b\,\Delta t)\cdot\mathbf n\,dA$ · *why:* A patch dA
     moving at b for Δt sweeps a thin prism whose height is the normal part of its displacement; sliding along the
     surface sweeps nothing (primer). Positive where it advances. · *plain:* Each piece of wall adds (or removes) a thin
     slab. · *set:* draws (bΔt)·n on one patch.
  9. *did:* F in the sliver ≈ its surface value · *tex:* $\displaystyle\int_{\Delta V}F\,dV\cong\int_{A^*(t)}F\,(\mathbf b\,\Delta t\cdot\mathbf n)\,dA$ ·
     *why:* Add the slabs of step 8. Within a slab of thickness O(Δt), F differs from its surface value by O(Δt)
     (mean-value theorem, P85), an O(Δt²) error. This is (3.34). · *plain:* The sliver's content is F at the wall times
     the swept volume. · *live:* "∫_{ΔV}F = … vs ∮F(bΔt·n) = … (gap …)".
  10. *did:* Divide by Δt · *tex:* $\dfrac{1}{\Delta t}\displaystyle\int_{\Delta V}F\,dV\to\int_{A^*}F\,\mathbf b\cdot\mathbf n\,dA$ · *why:* Δt
      comes out of the surface integral and cancels; the O(Δt²) error of step 9 becomes O(Δt) and vanishes in the limit. ·
      *plain:* The rate at which the moving walls sweep F in or out.
  11. *did:* Assemble the theorem · *tex:* $\dfrac{d}{dt}\displaystyle\int_{V^*(t)}F\,dV=\int_{V^*(t)}\dfrac{\partial F}{\partial t}dV+\int_{A^*(t)}F\,\mathbf b\cdot\mathbf n\,dA$ ·
      *why:* Put steps 7 and 10 into (3.33). This is (3.35). · *plain:* Change inside = change in place + what the walls
      sweep. · *live:* "volume … + surface … = …".
  12. *did:* Check the 1-D case · *tex:* $A^*=\{a,b\},\ \mathbf n=\mp1\ \Rightarrow\ \dfrac{d}{dt}\displaystyle\int_a^bF\,dx=\int_a^b\dfrac{\partial F}{\partial t}dx+\dot bF(b)-\dot aF(a)$ ·
      *why:* In one dimension the 'surface' is the two end points with outward normals −1 at a and +1 at b; the surface
      integral becomes a two-term sum, and (3.30) returns. · *plain:* Leibniz is the Reynolds transport theorem in 1-D. ·
      *set:* mode = 1-D.
- **Result.** $\frac{d}{dt}\int_{V^*(t)}F(\mathbf x,t)dV=\int_{V^*(t)}\frac{\partial F(\mathbf x,t)}{\partial t}dV+\int_{A^*(t)}F(\mathbf x,t)\,\mathbf b\cdot\mathbf n\,dA$
  (3.35) — *in words:* the rate of change of what is inside a moving volume = what changes in place + what the moving
  surface sweeps in (b·n > 0) or out (b·n < 0).
- **Check.** Units: [F]·m³/s on every term ✓. b = 0: d/dt passes inside the integral ✓. F = 1: dV*/dt = ∮b·n dA ✓ (D23).
  Numbers: growing sphere R = 1 m, Ṙ = 0.1 m/s, F = t at t = 1 s: 4.189 + 1.257 = 5.445 = d/dt(tV) ✓; Ex. 3.2: 0.1047
  m³/s both ways ✓; the (3.32) terms for the growing sphere: T4 ∝ Δt² (slope 2.00) ✓.
- **Sympy check (required, `check_src`)** — re-runs the construction on a growing sphere with F = t x₁²:
  ```python
  import sympy as sp                                             # symbolic algebra
  t, dt, R0, Rd = sp.symbols('t Delta_t R_0 Rdot', positive=True) # time, step, radius at t = 0, growth rate
  r, th, ph = sp.symbols('r theta phi', nonnegative=True)         # spherical coordinates of a point
  R = R0 + Rd*t                                                  # radius of V*(t); on its surface b·n = Rdot
  F = t*(r*sp.sin(th)*sp.cos(ph))**2                             # F = t x1^2 in spherical coordinates
  dV = r**2*sp.sin(th)                                           # volume element r^2 sin(theta) dr dtheta dphi
  def vol(expr, r0, r1):                                         # integral of expr over the shell r0 < r < r1
      return sp.integrate(expr*dV, (r, r0, r1), (th, 0, sp.pi), (ph, 0, 2*sp.pi))
  lhs = sp.diff(vol(F, 0, R), t)                                 # d/dt of the integral over V*(t): left of (3.31)
  vterm = vol(sp.diff(F, t), 0, R)                               # step 7: integral of dF/dt over V*
  sterm = sp.integrate((F*Rd*dV).subs(r, R), (th, 0, sp.pi), (ph, 0, 2*sp.pi))  # step 10: surface F b.n dA, dA = R^2 sin
  print(sp.simplify(lhs - (vterm + sterm)))                      # 0: (3.35) holds (step 11)
  T3 = vol(F, R, R.subs(t, t + dt))                              # step 1: F over the swept shell (third term of (3.32))
  T4 = vol(dt*sp.diff(F, t), R, R.subs(t, t + dt))               # fourth term of (3.32)
  print(sp.series(T4, dt, 0, 3).removeO())                       # starts at dt^2: second order (step 5)
  print(sp.simplify(sp.limit(T3/dt, dt, 0) - sterm))             # 0: (3.34), the sliver becomes the surface term (steps 9-10)
  ```
- **What it means.** The theorem is pure kinematics — no physics law is used. With b = u (a material volume) it gives
  the rate of change of a lump of fluid's mass, momentum or energy, which Ch. 4 sets equal to zero, the net force, or
  the heat and work; with b = 0 it gives the fixed-control-volume forms engineers use. It fails only if F or the surface
  is not smooth enough (a shock inside V*, a surface that tears).
- **Traps.** Using an unsigned swept volume (and then adding a separate "outflow" term). Keeping the ΔV·Δt term.
  Evaluating F in the sliver at t + Δt. Moving d/dt inside ∫_{V*(t)} when b ≠ 0.

### D23 · With F = 1 and b = u the theorem gives the volumetric strain rate, Eq. (3.14) — ★★, 5 steps, in C15 (notebook · `reynolds_transport_cv`)
- **Goal.** Show that the Reynolds transport theorem, applied to a small material volume with F = 1, reproduces
  (3.14) — the book's first "physical interpretation", left to Exercise 3.28.
- **Start.** (3.35) with F = 1 and b = u: $\dfrac{d}{dt}\displaystyle\int_{\delta V}dV=\int_{\delta V}\dfrac{\partial 1}{\partial t}dV+\int_{A^*}\mathbf u\cdot\mathbf n\,dA$
  — *in words:* the volume of a lump of fluid changes by what its surface sweeps.
- **Plan.** (1) Simplify the left side and the volume term. (2) Turn the surface integral into a volume integral (Gauss).
  (3) Shrink the volume (mean value).
- **Tools.** (3.35) (D22) · Gauss' divergence theorem (2.30) (Ch. 2 §2.12) · mean-value theorem for integrals (ch02 P85).
- **Assumptions.** Material volume: its surface moves with the fluid, b = u (start); u smooth (step 2).
- **Steps.**
  1. *did:* Simplify both sides · *tex:* $\dfrac{D(\delta V)}{Dt}=\displaystyle\oint_{A^*}\mathbf u\cdot\mathbf n\,dA$ · *why:* ∫dV is the
     volume, and its rate following the fluid is D/Dt; ∂1/∂t = 0 removes the volume term. · *plain:* A lump's volume
     changes by the net outward flow through its skin.
  2. *did:* Apply Gauss' theorem · *tex:* $\displaystyle\oint_{A^*}\mathbf u\cdot\mathbf n\,dA=\int_{\delta V}\nabla\cdot\mathbf u\,dV$ ·
     *why:* Divergence theorem (2.30), $\int_V\nabla\cdot\mathbf Q\,dV=\oint_A\mathbf Q\cdot\mathbf n\,dA$ (Ch. 2 §2.12) with Q = u. ·
     *plain:* Outflow through the skin = divergence summed inside. · *live:* "∮u·n ds = 1.257 = ∫∇·u dA = 0.2 × 6.283".
  3. *did:* Mean-value theorem · *tex:* $\displaystyle\int_{\delta V}\nabla\cdot\mathbf u\,dV=(\nabla\cdot\mathbf u)(\mathbf x^*)\,\delta V$ · *why:*
     For a continuous integrand the integral equals its value at some point x* inside times the volume (P85). · *plain:*
     A small volume sees an average divergence.
  4. *did:* Divide by δV · *tex:* $\dfrac{1}{\delta V}\dfrac{D(\delta V)}{Dt}=(\nabla\cdot\mathbf u)(\mathbf x^*)$ · *why:* δV ≠ 0;
     combine steps 1–3. · *plain:* The fractional growth rate is the divergence somewhere inside.
  5. *did:* Shrink the volume to a point · *tex:* $\dfrac{1}{\delta V}\dfrac{D(\delta V)}{Dt}=\dfrac{\partial u_i}{\partial x_i}$ · *why:* As
     δV → 0 around x, x* → x and continuity gives ∇·u(x); in index form this is (3.14),
     $\frac{1}{\delta V}\frac{D}{Dt}(\delta V)=\frac{\partial u_i}{\partial x_i}=S_{ii}$. · *plain:* The same result as D11, now from the
     transport theorem. · *set:* the bars reduce to the area growth rate.
- **Result.** $\frac{1}{\delta V}\frac{D}{Dt}(\delta V)=\frac{\partial u_i}{\partial x_i}$ (3.14) — *in words:* divergence is the
  fractional growth rate of a small material volume, whichever route you take.
- **Check.** Units 1/s ✓. u = (x, y, z): 3 s⁻¹ ✓. E7's ellipse (u = (0.1x, 0.1y)): dA/dt = 0.2 × 6.283 = 1.257 m²/s ✓.
  `ch03.material_volume_rate` returns equal surface flux and ∫∇·u (1e-10) ✓.
- **What it means.** D11 (edge by edge) and D23 (whole surface) agree — the kinematic half of the continuity equation of
  Ch. 4 holds for blobs of any shape.
- **Traps.** Using b = 0 (then the surface term vanishes and the volume "does not change"). Dividing by δV after the
  limit instead of before.

### D24 · For a small material volume the theorem gives back the material derivative, Eq. (3.5) — ★★, 7 steps, in C15 (notebook)
- **Goal.** Show the book's second interpretation (Exercise 3.30) — that (3.35) "extends (3.5) to finite volumes" — by
  deriving the material derivative *from* the transport theorem, and expose the extra term F∇·u that the book's wording
  hides (the step Ch. 4 uses to derive the continuity and momentum equations).
- **Start.** (3.35) with b = u on a small material volume δV:
  $\dfrac{d}{dt}\displaystyle\int_{\delta V}F\,dV=\int_{\delta V}\dfrac{\partial F}{\partial t}dV+\oint_{A^*}F\,\mathbf u\cdot\mathbf n\,dA$ — *in words:* the
  amount of F carried by a lump of fluid changes in place and by what its moving skin sweeps.
- **Plan.** (1) Turn the surface term into a volume integral and expand the divergence of a product. (2) Shrink the lump
  so the integrals become values times δV. (3) Split the left side with the product rule and remove the lump's swelling
  with (3.14). (4) What is left is the rate following the lump.
- **Tools.** (3.35) (D22) · Gauss' divergence theorem (Ch. 2 §2.12) · product rule for ∇·(Fu) (gloss in step 3) ·
  mean-value theorem (ch02 P85) · product rule (ch01 P38) · (3.14), $\frac{1}{\delta V}\frac{D(\delta V)}{Dt}=\frac{\partial u_i}{\partial x_i}$
  (C09, D23).
- **Assumptions.** Material volume, b = u (start); F and u smooth (steps 1, 4); δV → 0 around a point x (step 4).
- **Steps.**
  1. *did:* Gauss on the surface term · *tex:* $\displaystyle\oint_{A^*}F\,\mathbf u\cdot\mathbf n\,dA=\int_{\delta V}\nabla\cdot(F\mathbf u)\,dV$ ·
     *why:* Divergence theorem with the vector field Fu (Ch. 2 §2.12); we want the whole budget as one volume integral. ·
     *plain:* What the skin sweeps = the divergence of the flux of F inside.
  2. *did:* Combine the two volume integrals · *tex:* $\dfrac{d}{dt}\displaystyle\int_{\delta V}F\,dV=\int_{\delta V}\Big[\dfrac{\partial F}{\partial t}+\nabla\cdot(F\mathbf u)\Big]dV$ ·
     *why:* Both terms of (3.35) are now integrals over the same volume; add the integrands. · *plain:* One integrand for
     the whole budget.
  3. *did:* Expand the divergence of a product · *tex:* $\nabla\cdot(F\mathbf u)=\mathbf u\cdot\nabla F+F\,\nabla\cdot\mathbf u$ · *why:*
     Product rule in index form: $\partial(Fu_i)/\partial x_i=u_i\,\partial F/\partial x_i+F\,\partial u_i/\partial x_i$ (gloss). ·
     *plain:* F carried around + F spread out by the swelling flow.
  4. *did:* Shrink the lump to a point · *tex:* $\dfrac{D}{Dt}(F\,\delta V)=\Big(\dfrac{\partial F}{\partial t}+\mathbf u\cdot\nabla F+F\,\nabla\cdot\mathbf u\Big)\delta V$ ·
     *why:* Mean-value theorem (P85): for small δV each integral is its integrand at x times δV; the time derivative of an
     integral over a material lump is the rate following the lump, written D/Dt. · *plain:* The budget of one tiny lump.
  5. *did:* Product rule on the left · *tex:* $\dfrac{D}{Dt}(F\,\delta V)=\delta V\,\dfrac{DF}{Dt}+F\,\dfrac{D(\delta V)}{Dt}$ · *why:*
     Product rule (P38) for the two factors F (the lump's value) and δV (its volume), both followed with the lump. ·
     *plain:* The content changes because F changes and because the lump swells.
  6. *did:* Insert the swelling rate (3.14) · *tex:* $\dfrac{D}{Dt}(F\,\delta V)=\delta V\,\dfrac{DF}{Dt}+F\,(\nabla\cdot\mathbf u)\,\delta V$ ·
     *why:* (3.14), $\frac{1}{\delta V}\frac{D(\delta V)}{Dt}=\frac{\partial u_i}{\partial x_i}$ (D23 derived it from this same
     theorem with F = 1). The book skips this. · *plain:* The swelling part is F times the divergence.
  7. *did:* Equate steps 4 and 6, cancel · *tex:* $\dfrac{DF}{Dt}=\dfrac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$ · *why:* Both lines
     describe the same quantity; the term F∇·u δV appears in each and cancels; divide by δV ≠ 0. This is (3.5). · *plain:*
     The rate following a lump is the local rate plus the advective rate — the material derivative, recovered.
- **Result.** $\frac{d}{dt}\int_{\delta V}F\,dV=\int_{\delta V}\big[\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F+F\nabla\cdot\mathbf u\big]dV$ and,
  for δV → 0, $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$ (3.5) — *in words:* the transport theorem for a
  small lump = the material derivative + the effect of the lump's volume changing.
- **Check.** Units [F]/s ✓. F = 1: steps 4 and 6 both give D(δV)/Dt = (∇·u)δV, i.e. D23 ✓. F = ρ with the lump's mass
  conserved: the left side of step 4 is 0, so $\frac{D\rho}{Dt}+\rho\nabla\cdot\mathbf u=0$ — Ch. 4's continuity equation
  (preview) ✓. Same result as D02, which reached (3.5) from the chain rule instead ✓.
- **What it means.** "(3.35) extends (3.5) to finite volumes" hides a term: per unit volume the content of a lump changes
  at DF/Dt + F∇·u. Setting F = ρ, ρu, ρe gives the three conservation laws of Ch. 4; with F = ρ × (something per unit
  mass) the F∇·u term is exactly what the continuity equation absorbs.
- **Traps.** Dropping F∇·u and concluding d/dt∫F dV = ∫DF/Dt dV (true only for incompressible flow). Using b = 0 instead
  of b = u. Treating F in D(FδV)/Dt as the field at a fixed point rather than the value carried by the lump.
