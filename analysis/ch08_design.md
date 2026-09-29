# Chapter 8 — Laminar Flow: lesson design
(from `analysis/ch08_curation.md` (A 15 · B 113 · C 15 = 143 rows; CORE 15 · NOTE 107 · RECAP 20 · SKIP 1; 33
derivations D01–D33 written out (★ 8 · ★★ 18 · ★★★ 7), 296 steps; E1–E9 + backup B1) and `analysis/ch08.md` (§2b
derivations a-D01…a-D49, §4 implementation rows I1–I41, §5 core modules `core/laminar.py`, `core/lubrication.py`,
`core/creeping.py`, §9 conventions and slips); every equation re-read on the rendered pages `chapters/pages/ch08/`
(printed = pdf − 27): (4.89), (8.1), (8.2, 8.3) p338 (printed 311); fully developed flow and Fig. 8.2 p339; (8.4a,b),
the twice-integrated line, (8.5), Q p340; V, Couette, Poiseuille, τ p341; tube equations, (8.6), τ_zR p342; (8.7),
(8.8), Q, V p343; circular Couette equations, (8.9), A, B, (8.10), (8.11) p344; σ_Rφ, (8.12) p345; (6.2), (8.13a, b),
(8.14), (8.15), (8.16a, b) p346; (8.17a, b), (8.18), (8.19) p347; Example 8.1 p348–p349; Example 8.2 p349–p351;
Example 8.3 p351–p353; (8.20)–(8.25) p354; chain-rule lines, (8.26)–(8.29) p355; (8.30) p356; (8.31), (8.32a, b)
p357; Example 8.4 p358; Example 8.5 p358–p361; Example 8.6 p361–p362; Example 8.7 p363–p364; (8.33)–(8.38) p364;
Fig. 8.16 and the 0.06U line p365; (8.39)–(8.42) p366; (8.43), Fig. 8.17 p367; ∇²ω = 0, ω_φ, (6.83), (8.44)–(8.47) p368;
f-ODE, (8.48)–(8.51), terminal balance p369; Millikan balance, (8.52), fluid-frame ψ p370; far-field orders p371;
Oseen substitution and equation p372; (8.53), Oseen C_D p373; §8.7 p374. 2026-09-28, lesson-designer. The implementer
works in parallel from analysis §4 + curation §9; **Part C is written first and is the contract both sides keep.**)

**Binding conventions for every builder (analysis §9, curation decisions 4–6, §8).**
1. **Imports and aliases.** `from fluidpy import ch08_laminar_flow as ch08`; `from fluidpy.core import laminar as LAM,
   lubrication as LUB, creeping as CRP`. **`ch08` re-exports every public name of `core/laminar.py`,
   `core/lubrication.py`, `core/creeping.py`, `core.diffusion.crank_nicolson_1d` and
   `core.diffusion.couette_startup_profile`** (the ch05–ch07 pattern), so explainer parity rows write `ch08.<name>`
   only. Units SI: y, h, a, R, r, L, δ [m]; t [s]; U, u, v, V [m/s]; ω [rad/s] (§8.5) or [1/s] (vorticity); ν [m²/s];
   μ [Pa s]; p [Pa]; dp/dx [Pa/m]; Q [m²/s] (per unit width) or [m³/s] (pipe); q [m²/s]; ψ [m²/s] (2-D) or [m³/s]
   (Stokes, axisymmetric); Γ [m²/s]; τ, σ [Pa]; W [N/m]; D [N]; torque per length [N m/m]; power per length [W/m].
2. **Defaults.** Water μ = 1.0 × 10⁻³ Pa s, ρ = 1000 kg/m³, ν = 1.0 × 10⁻⁶ m²/s unless a function is about air or oil;
   `core.creeping` and `core.lubrication` (thin films) use **g = `G0` = 9.80665** (settling physics, SI standard) —
   stated in every docstring (ch04, ch05, ch07 use `G_BOOK` = 9.81; the book prints no g-dependent number here).
   `circular_couette(R2=np.inf)` and `R1=0.0` are exact limit branches (never 1e300); sphere fields return NaN for
   r < a; erfc (never 1 − erf) for Stokes' first problem; −expm1 in the line vortex and in Oseen's 1 − e^{−s}.
3. **Scalar-callable and parity-friendly.** Every public function accepts Python floats and returns a float, a tuple of
   floats or a `dict` of floats (arrays only when the input is an array). Parity `py:` expressions use `ch08.…`, `np.pi`,
   numbers, strings, lists, dicts and keyword arguments only, indexed down to one float (`shot.py` evaluates them with no
   builtins). **`Viz.num.erf/erfc` are accurate to ~1.2 × 10⁻⁷ only** (CUMULATIVE): E5, E6 carry a local
   double-precision `erfcHP` (ch04 E4's continued fraction) and their erf parity rows use rtol 1e-10 with it.
4. **Signs and coordinates** (⚠️ callouts where first used, with numbers): the book passes **dp/dx** (favourable < 0);
   ch04 code used **G = −dp/dx** — new functions take `dpdx` and accept `G=` as an alias (a parity row pins the sign).
   Channel walls at **y = 0 (fixed) and y = h (moving at U)**. Capital **R** = cylindrical radius (§8.2), lower-case **r**
   = spherical (§8.6) or plane-polar (Example 8.6). **θ in §8.6 is measured from the downstream +x axis** (rear
   stagnation point θ = 0, front θ = π; `frame="body"` = sphere at rest in U e_x, `"fluid"` = sphere moving to −x);
   Example 8.6's θ is a plane-polar azimuth. **φ** is an azimuth except in Hele-Shaw (velocity potential). **ω** is
   vorticity (§8.1, Example 8.5, §8.6 ω_φ) but the oscillation frequency in §8.5. **η = y/√(νt)** in (8.25), but the book's
   Figs. 8.13–8.14 plot y/(2√(νt)): every figure and explainer axis says which (`similarity_variable(..., half=True)`).
   Stokes-layer depth: book **δ ~ 4√(ν/ω)** vs the literature's e-folding depth **δ_e = √(2ν/ω)** (ratio 2√2 = 2.83).
   Four Reynolds numbers: pipe **Ud/ν** (diameter, mean velocity); lubrication **Re_L = ρUL/μ** (passage length; the
   group that matters is ε²Re_L); generic **ρUL/μ** (8.40); sphere **Re = 2aU/ν** (diameter) in (8.52), (8.53) and the
   Oseen C_D — the radius-based Re_a = Re/2 turns Oseen's 3/16 into 3/8. Lubrication **p* = p/P_a** (atmospheric) makes
   the bearing number Λ = μUL/(P_a h²) large in a real bearing (Λ = 49 for our engine film; ≈ 10³ for thinner, longer
   films); the natural pressure scale is **μUL/h²** (the scale that makes Λ = 1).
5. **Book slips taught in corrected form** (each gets "the book prints X; the correct form is Y" where it is used, a
   discriminating wrong-variant option in code, never asserted as physics; analysis §9 labels in brackets):
   **(8.13b) prints ∂p/∂x** in the v-equation → ∂p/∂y [R6] (N24, D10 step 12) · **(8.17a) and Example 8.2 lack the ν**
   → $0\cong-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{\partial^2u}{\partial y^2}$ [R7] (C05, D11 step 6, N36) ·
   **Example 8.1**: the intermediate integrals print $(1-\alpha x/L)$ → $(1+\alpha x/L)$, and the final pressure prints the
   denominator to the first power → **squared** [R8] (C07, D14 steps 6 and 13; `slider_bearing(model="book")` is the
   ghost) · **(8.19) adds U₀** instead of U₀(1 − y/h) [R8b] (N34, D12 step 8; `lubrication_velocity(form="book")`) ·
   **channel V ≡ Q/h** prints the middle integral without 1/h [R9] (N10) · **∫₀^∞ω dy = +U**, printed −U [R10] (N53) ·
   **Example 8.5's ±2.76** is ±2.772 [R11] (N59, D22 check) · **rear pressure minimum −3μU/2a**, printed without the minus
   [R12] (N88, D30 step 11) · **Oseen's equation needs −∂p/∂x_i** [R13] (N95) · **cross-reference slips** [R14]: "(9.63)"
   → (8.43) (N93), "(9.68)" → (8.48) (N98), "Substitution of (8.33) into (8.20)" → (8.35) (N68), "the final equation of
   Example 8.2" → Example 8.3 (N63), Example 8.2's "y = 0 … y = h" → z (N36), §8.1 "kinematic viscosity" for μ → dynamic
   (C01), "(8.16a) can be integrated twice" → (8.17a) (N33, found in this design) · **power into the fluid is
   −2πR₁σ_Rφu_φ** [R15] (R11) · the book's lubrication example quotes the *atmospheric* Λ "near unity" although it is
   ≈ 10³ [R24] (N31).
6. **Colours (text, figures, derivation terms, explainers — one meaning each; curation §5):** pressure and the
   pressure-driven (Poiseuille) part `orange` · viscous and the wall-driven (Couette) part `rose` · inertia `teal` ·
   velocity profile `blue` · vorticity `accent` purple · backflow and warnings `amber` · similarity / rescaled curves purple
   dashed · analytic ghost muted dashed · numerical dots teal · printed-slip ghosts muted dotted, labelled "as printed".
7. **Every book equation is shown in full next to its number** (CLAUDE.md rule 3) — in markdown, derivation steps
   ("substitute (8.25), $u/U=F(y/\sqrt{\nu t})$"), traps, recaps, notes and every explainer text. Builders reuse
   `show_eqs(text, EQ)` from `notebooks/build_ch07.py` with an `EQ` dict for all 58 labels of ch08 plus the earlier
   labels cited here ((4.5), (4.10), (4.39b), (4.85), (4.89), (4.100), (4.107), (5.1), (5.2), (5.13), (3.29), (6.2),
   (6.10), (6.12), (6.77), (6.83), (6.86)); JS explainers use a local `showEqs`. `tools/eq_refs.py` must list 0
   offenders. Exercises are cited as "Exercise 8.19", never as bare equations.
8. **nbkit behaviour** (ch04–ch07 lesson): `nb.recap(...)` and `nb.section(...)` end the current CORE block, so every
   RECAP sits **before** the `nb.core` call of the block that uses it (R01, R02 before C01; R07, R08 before C02; R09
   before C03; R10–R12 before C04; R13 before C05; R14 before C09; R15, R16 before C12; R17–R19 before C13; R20 before
   C14), except R03–R06, which follow the completed C01 block as "the equations every solution below must satisfy" (C01
   has its code and visuals before them). Every `nb.derivation` sits inside its curation CORE block with `ref=` a bare
   equation label ("8.5").
9. **Parsing.** Part E is the only part with table rows after its heading; Part F has no line starting with a table
   bar (absolute values are written \lvert … \rvert or in words). Explainer headings are exactly `### E1 ·
   couette_poiseuille_backflow` … `### E9 · stokes_drag_settling` and `### B1 · rotating_cylinders_couette`. Primer
   terms in Part A are the exact Concept text of the Part E rows marked "primer" (coverage_check matches the first 18
   characters). "Explained by" never names an earlier chapter's C-number or an analysis slip label (it cites "Ch. 4
   §4.10" or a primers.md entry).
10. **Public repo (rule 9).** No exercise text; the book's lubrication example (0.1 mm gap, 25 cm, 10 m/s, 30-weight
    oil, ε²Re_L = 0.001), the δ₉₉ factor as a book value, the printed ±2.76, Example 8.1's printed constants, Fig. 8.15's
    νt values and every exercise answer stay in `tests/book_values_ch08.json`. The notebook shows **our** numbers from
    our inputs (engine film h = 50 µm, L = 5 cm, U = 5 m/s, ν = 10⁻⁴ m²/s, μ = 0.05 Pa s); numbers that coincide with a
    printed value (3.64, 0.06, 5.54) are built from the function and printed at 4 significant figures (3.643, 0.05911,
    5.544). Figures are drawn by our code (Reynolds's apparatus is our synthetic dye-streak sketch, labelled schematic).

Order of parts: C (contract) · A (notebook storyboard) · B (explainers) · D (runtime) · E (prerequisite ledger) · F
(derivations).

---

## Part C — functions the builders will call (the implementer's contract)

**Status column.** **I#** = planned in `analysis/ch08.md` §4 row I# (signature kept or refined here). **§9** = added by
the curation's §9 table (keys fixed here). **NEW** = added by this design (in neither) — flagged as the phase asks.
**reuse** = exists (ch01–ch07) and is only called. Docstrings cite § and Eq. (from the rendered page), list symbols with
units and assumptions, and carry the validation label. Return types: `float`/`ndarray` for single quantities, tuples
where the analysis fixed them, `dict` (keys listed) for several named results.

### C.0 Reused (existing; called, not changed)
| # | Callable (signature) | Returns | Used by | Status |
|---|---|---|---|---|
| 0.1 | `style.setup_notebook() -> bool` (FAST) · `style.COLORS` · `style.savefig(fig, "ch08", name)` | FAST · palette · path | setup, every figure | reuse |
| 0.2 | `anim.animate(update, frames, fig, interval)` · `anim.show_animation(anim, player="video"\|"frames")` | HTML | C08 (A3), C09 (A1), C10 (A5), C11 (A2), C13 (A4) | reuse |
| 0.3 | `interact.slider_figure(fn, name, values, *, unit, xlabel, ylabel, title, xrange, yrange, height)` · `interact.animate_figure(frame_fn, times, …)` · `interact.live(fn, **widgets)` | plotly Figure / widget | C02 (IF1), C03 (IF2), C04 (IF3), C07 (IF4, live), C09 (IF5), C11 (IF6), C12 (IF7), C15 (IF8) | reuse |
| 0.4 | `embed.show_viz("ch08", slug)` | display | 9 explainer cells | reuse |
| 0.5 | `tools.convergence.observed_order(h, err) -> float` | log–log slope | C09 (CN order 2), C10 (fitted exponents), C15 (ratio slope 1) | reuse |
| 0.6 | `core.navier_stokes.ns_incompressible_terms(u, p, x, t=0.0, rho=1000.0, mu=1e-3, …)` · `exact_solution(name, x, t=0.0, **p)` ("couette" with `G`, "pipe_poiseuille", "stokes_first", "lamb_oseen", "solid_body") · `viscous_force_forms(u, x, t, mu, h)` | term dict / (u, p) | R03 (C01 residual demo), C02, C03, C04 (N22), C09, C10 parity | reuse |
| 0.7 | `core.similarity.reynolds_number(U, l, nu=None, rho=None, mu=None)` · `nondimensional_ns_coefficients(pressure_scale="dynamic"\|"viscous", time_scale="advective")` · `sphere_drag_coefficient(Re, model="morrison"\|"stokes")` | [–] · sympy dict · C_D | R01, R16, C12, C14 (N92 parity), C15 (C_D figure) | reuse |
| 0.8 | `core.diffusion.ftcs_diffusion_1d(f0, D, dy, dt, nsteps, …)` · `stable_time_step(D, dy, safety=0.9)` · `couette_startup_profile(y, t, U, h, nu, nterms=200)` | arrays · dt · u [m/s] | R14 (C09 contrast), C09 (N39 FTCS vs CN), E5 preset | reuse (re-exported by `ch08`) |
| 0.9 | `core.vortices.gaussian_vortex(r, Gamma, sigma)` · `solid_body_rotation(r, omega0)` · `line_vortex(r, B)` | u_θ [m/s] | R12 (C04 parity), N61 (C10 parity σ = 2√(νt)) | reuse |
| 0.10 | `core.curvilinear.curl(A, system, coords)` · `vector_laplacian(A, system)` · `advective_acceleration(u, system)` · `strain_rate(u, system)` (sympy; `"cylindrical"`, `"spherical"`) | sympy lists / Matrix | R09 (C03), C04 (N22), R17 (C13), C14 (D31 check) | reuse |
| 0.11 | `core.potential.sphere(U, a)` (ideal-flow sphere, `AxisymFlow`) · `axisym_velocity_spherical(fn, r, theta, kind="psi")` · `axisym_velocity_spherical_sym(expr, r, theta)` · `Flow` (uniform + doublet cylinder for Hele-Shaw streamlines) | objects / (u_r, u_θ) | R18 (C13), N93 (C13 figure), N36 (C06 figure) | reuse |
| 0.12 | `core.laplace_solvers.solve_poisson(mask, f, bc, method="direct")` | p on a masked grid | N36 (C06 V3 check, FAST 64²) | reuse |
| 0.13 | `ch01.pi_groups(variables)` · `ch01.water_viscosity(T)` · `ch01.sutherland_viscosity(T)` · `ch01.kinematic_viscosity(mu, rho)` · `ch01.water_density(T)` | groups · μ · ν · ρ | C01 (D01), N43 (C09), N73 (C11), R20 (C14) | reuse |
| 0.14 | `ch03.pipe_profile(r, z, U_mean, R, n0, L_e)` (z → ∞ gives the parabola) | u [m/s] | C03 parity | reuse |
| 0.15 | `ch04.stokes_first_problem(y, t, U, nu)` (becomes a re-export of 1.16) · `ch04.plane_poiseuille(y, G, h, mu)` (**G = −dp/dx**) | u [m/s] | R08, C09 parity | reuse |
| 0.16 | `ch05.rotating_cylinder_flow(r, a, omega)` (**ω there = 2Ω₁**) · `ch05.diffusing_vortex_sheet(y, t, gamma, nu)` (**γ = u_below − u_above = −2U**) · `ch05.dissipation_outside_cylinder(a, R_out, Gamma, mu, rho)` | u_θ · u · dict | R10, R11 (C04), N59 (C10) parity | reuse |
| 0.17 | `core.thermo.G_BOOK` · `G0` · `P_ATM` | constants | throughout | reuse |

### C.1 `fluidpy/core/laminar.py` (NEW module, `core.LAM`; re-exported by `ch08`) — exact parallel and unsteady viscous flows for Ch. 8, 9, 11, 12, 13, 16
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 1.1 | `channel_flow(y, h, U=0.0, dpdx=0.0, mu=1e-3, G=None)` (`G` = −dp/dx alias; both given → `ValueError`) | (8.5) $u=\frac Uhy-\frac1{2\mu}\frac{dp}{dx}y(h-y)$ | u [m/s] | C02, E1, IF1, R03 residual, N03 | I5 |
| 1.2 | `channel_flow_rate(h, U=0.0, dpdx=0.0, mu=1e-3)` | $Q=\frac{Uh}{2}[1-\frac{h^2}{6\mu U}\frac{dp}{dx}]$ (written as $Uh/2-\frac{h^3}{12\mu}\frac{dp}{dx}$ so U = 0 works), V = Q/h | (Q [m²/s], V [m/s]) | C02 (N10), E1 | I6 |
| 1.3 | `channel_shear_stress(y, h, U=0.0, dpdx=0.0, mu=1e-3)` | τ = μ du/dy = μU/h − (h/2 − y) dp/dx (signed) | τ [Pa] | C02 (D05, N11), E1 | I6 |
| 1.4 | `channel_backflow_threshold(U, h, mu=1e-3)` | dp/dx* = 2μU/h² (ours, D05) | [Pa/m] | C02 (N09), E1 | I6 |
| 1.5 | `couette_poiseuille_state(h, U, dpdx, mu=1e-3)` | (8.5), Q, τ, D05 | dict(Q, V, Q_couette (Uh/2), Q_poiseuille (−h³dp/dx/12μ), tau_bottom, tau_top, backflow (bool, dp/dx > 2μU/h² for U > 0; mirrored for U < 0), threshold, ratio (dp/dx ÷ threshold; NaN if U = 0), u_max, y_umax, y_reversal (h − 2μU/(h dp/dx) inside (0, h), else NaN), zero_flow_dpdx (6μU/h²)) | E1, IF1, C02 figure | §9 (keys fixed here) |
| 1.6 | `pipe_poiseuille(R, a, dpdz, mu=1e-3, G=None)` | (8.6) $u_z=\frac{R^2-a^2}{4\mu}\frac{dp}{dz}$ | u_z [m/s] | C03, IF2 | I7 |
| 1.7 | `pipe_shear_stress(R, dpdz)` · `pipe_wall_stress(a, dpdz)` | (8.7) $\tau=\frac R2\frac{dp}{dz}$, (8.8) $\tau_0=\frac a2\frac{dp}{dz}$ | τ [Pa] | C03 | I8 |
| 1.8 | `pipe_flow_rate(a, dpdz, mu=1e-3)` | $Q=-\frac{\pi a^4}{8\mu}\frac{dp}{dz}$, $V=Q/\pi a^2$, u_max = 2V | (Q [m³/s], V, u_max [m/s]) | C03 (N16), IF2 | I8 |
| 1.9 | `pipe_friction_factor(Re)` | f = 64/Re (Darcy, Re = Vd/ν) | [–] | C03 (D07) | I8 |
| 1.10 | `circular_couette(R, R1, R2, Omega1, Omega2, return_coeffs=False)` (`R2=np.inf` → Ω₁R₁²/R, requires Ω₂ = 0; `R1=0.0` → Ω₂R) | (8.9) $u_\varphi=AR+B/R$, (8.10), (8.11), (8.12) | u_φ [m/s] (+ (A [1/s], B [m²/s])) | C04, IF3, B1 | I9 |
| 1.11 | `circular_couette_pressure(R, R1, R2, Omega1, Omega2, rho=1000.0, p1=0.0)` | N20 $p=p_1+\rho[\frac{A^2}2(R^2-R_1^2)+2AB\ln\frac R{R_1}-\frac{B^2}2(\frac1{R^2}-\frac1{R_1^2})]$ (ours) | p [Pa] | C04 (N20), B1 | I10 |
| 1.12 | `circular_couette_shear_stress(R, R1, R2, Omega1, Omega2, mu=1e-3)` | σ_Rφ = μR d(u_φ/R)/dR = −2μB/R² | [Pa] | R10 (C04) | I10 |
| 1.13 | `circular_couette_power(R1, R2, Omega1, Omega2, mu=1e-3)` (R2 = inf allowed) | R11: power in = −2πR₁σ_Rφu_φ\|_{R₁} (positive), ∫ε dA by `quad` | dict(torque_inner, torque_outer [N m/m], power_in, power_out, dissipation [W/m]) | R11 (C04) | I10 |
| 1.14 | `circular_couette_state(R1, R2, Omega1, Omega2, mu=1e-3, rho=1000.0)` | (8.10), N20, R10, R11; Rayleigh: (R u_φ)² increasing outward | dict(A, B, torque_inner, torque_outer, power_in, dissipation, dp_gap (p(R₂) − p(R₁)), rayleigh_stable (bool)) | B1, IF3 title | §9 |
| 1.15 | `similarity_variable(y, t, nu, half=False)` | (8.25) η = y/√(νt); `half=True` → y/(2√(νt)) (the figures' axis) | [–] | C09, E5 | I22 |
| 1.16 | `stokes_first_problem(y, t, U=1.0, nu=1e-6)` (moved from ch04, re-exported there) | (8.30) $\frac uU=1-\mathrm{erf}\frac{y}{2\sqrt{\nu t}}$ (via erfc) | u [m/s] | C09, E5, A1, IF5 | I22 |
| 1.17 | `stokes_first_vorticity(y, t, U=1.0, nu=1e-6)` | N53 ω = (U/√(πνt))e^{−y²/4νt} | ω [1/s] | C09 (N53), E5 | I22 |
| 1.18 | `stokes_first_stopped(y, t, T, U=1.0, nu=1e-6)` | N56 $u=U[\mathrm{erfc}\frac{y}{2\sqrt{\nu t}}-\mathrm{erfc}\frac{y}{2\sqrt{\nu(t-T)}}]$ for t > T (ours) | u [m/s] | C09 (N56), E5 preset | I24 |
| 1.19 | `stokes_first_state(t, U=1.0, nu=1e-6, level=0.01, rho=1000.0)` | (8.30), (8.31), N53 | dict(sqrt_nut, eta_edge (2 erfcinv(level)), delta, tau_w (ρνU/√(πνt)), vorticity_content (= U)) | E5 | §9 (`rho` added) |
| 1.20 | `diffusion_thickness(t, nu, level=0.01)` | (8.31) δ = 2 erfcinv(level)√(νt) | δ [m] | C09 (D20), N03 edge estimate, E5 | I23 |
| 1.21 | `transition_width(t, nu, level=0.95)` | N59 width = 4 erfinv(level)√(νt) | [m] | C10, E6 | I23 |
| 1.22 | `vortex_sheet_diffusion(y, t, U=1.0, nu=1e-6)` | N59 $u=U\,\mathrm{erf}\frac{y}{2\sqrt{\nu t}}$, $\omega_z=-\frac{U}{\sqrt{\pi\nu t}}e^{-y^2/4\nu t}$ | (u [m/s], ω_z [1/s]) | C10, E6, A5 | I26 |
| 1.23 | `temporal_bl_wall_stress(t, U, nu=1e-6, rho=1000.0)` | N60 τ_w = μU/√(πνt), $C_f=\frac2{\sqrt\pi}\sqrt{\frac{\nu}{U^2t}}$, Re_x = U²t/ν | dict(tau_w [Pa], Cf, Rex) | C10 (N60) | I27 (keys fixed) |
| 1.24 | `line_vortex_decay(r, t, Gamma, nu=1e-6)` · `line_vortex_spinup(r, t, Gamma, nu=1e-6)` | N61 $u_\theta=\frac{\Gamma}{2\pi r}[1-e^{-r^2/4\nu t}]$ (−expm1; r = 0 → 0), N62 $\frac{\Gamma}{2\pi r}e^{-r^2/4\nu t}$ | u_θ [m/s] | C10, E6, A5 | I28 |
| 1.25 | `stokes_second_problem(y, t, U=1.0, omega=2*np.pi, nu=1e-6)` | (8.38) $u=Ue^{-y\sqrt{\omega/2\nu}}\cos(\omega t-y\sqrt{\omega/2\nu})$ | u [m/s] | C11, E7, A2, IF6 | I29 |
| 1.26 | `stokes_layer(nu, omega)` | D25: δ_e = √(2ν/ω), δ_book = 4√(ν/ω), crest speed √(2νω), wavelength 2πδ_e | dict(delta_e, delta_book, phase_speed, wavelength) [m, m/s] | C11 | I29 (the callable `envelope` dropped: not parity-friendly) |
| 1.27 | `stokes_layer_state(nu, omega, y)` | (8.38), D25 | dict(delta_e, delta_book, amplitude (e^{−y/δ_e}), amp_at_book_depth (e^{−2√2}), phase_lag (y/δ_e rad), time_lag (s), crest_speed, wavelength, period) | E7 | §9 |
| 1.28 | `core.diffusion.crank_nicolson_1d(u0, y, dt, nsteps, D, bc_left, bc_right=0.0, startup_be=2, return_all=False)` (bc: float or callable of t; `startup_be` backward-Euler steps damp the impulsive-start oscillation) | (8.20) discretised, second order | u (ny,) or (t (nsteps+1,), U (nsteps+1, ny)) | C09 (N39, from scratch), A1, C11 check | I25 (`startup_be`, `return_all` added) |

### C.2 `fluidpy/core/lubrication.py` (NEW module, `core.LUB`; re-exported by `ch08`) — thin gaps and thin films for Ch. 8, 9, 13, 16
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 2.1 | `lubrication_scales(L, h, U, rho, mu, p_a=P_ATM)` | (8.14), (8.16a) groups | dict(eps, Re_L, eps2_Re_L, Lambda (μUL/(p_a h²)), p_visc (μUL/h² [Pa])) | C05 (N31), E2 | I12 |
| 2.2 | `lubrication_term_magnitudes(L, h, U, rho, mu, p_scale="viscous", p_a=P_ATM)` (`"viscous"`: Λ = 1; `"atm"`: Λ from p_a) | coefficients of (8.16a) and (8.16b) | dict(x_inertia (ε²Re_L), x_pressure (1/Λ), x_diff_along (ε²), x_diff_across (1), y_inertia (ε⁴Re_L), y_pressure (1/Λ), y_diff_along (ε⁴), y_diff_across (ε²)) | C05 figure, E2 | §9 (keys fixed here) |
| 2.3 | `lubrication_velocity(y, h, dpdx, U_h=0.0, U_0=0.0, mu=1e-3, form="consistent")` (`"book"` = printed U₀ placement, wrong-variant) | (8.19) consistent form $U_h\frac yh+U_0(1-\frac yh)$ | u [m/s] | C06, E3 | I14 |
| 2.4 | `lubrication_flux(h, dpdx, U_h=0.0, U_0=0.0, mu=1e-3)` | N32 $q=-\frac{h^3}{12\mu}\frac{\partial p}{\partial x}+\frac{(U_0+U_h)h}{2}$ | q [m²/s] | C06 (D13), E3 | I14 |
| 2.5 | `reynolds_pressure_1d(x, h, U_0=0.0, U_h=0.0, mu=1e-3, p_left=0.0, p_right=0.0)` (h: array on x or callable) | steady Reynolds: q constant, dp/dx = 12μ(q̄ − q)/h³, q̄ = (U₀ + U_h)h/2 | (p [Pa] on x, q [m²/s]) | C06 code, C07 check | I15 |
| 2.6 | `slider_bearing(x, h0, alpha, L, U, mu=1e-3, p_e=0.0, model="exact")` (`"linear"` O(α); `"book"` printed first-power denominator, wrong-variant) | Example 8.1 $p-p_e=\frac{6\mu LU}{h_o^2}\frac{\alpha(x/L)(1-x/L)}{(2+\alpha)(1+\alpha x/L)^2}$ | p [Pa] | C07, E3, IF4 | I16 |
| 2.7 | `slider_bearing_load(h0, alpha, L, U, mu=1e-3, model="exact")` (series for \|α\| < 1e-3) | $W=\frac{\alpha\mu L^2U}{2h_o^2}$ (linear); N35 $W=\frac{6\mu UL^2}{h_o^2\alpha^2}[\ln(1+\alpha)-\frac{2\alpha}{2+\alpha}]$ (exact, ours) | W [N/m] | C07, E3 | I16 |
| 2.8 | `slider_optimum_taper()` | max of the exact W(α) (`minimize_scalar`) | dict(alpha_opt (= 1.1889), K_opt (1 + α), W_star (W h₀²/(6μUL²) = 0.02671)) | C07 (N35), E3 preset | I16 (dict) |
| 2.9 | `slider_bearing_state(h0, alpha, L, U, mu=1e-3, p_e=0.0)` | Example 8.1, N35 | dict(C1 (the pad-frame flux ∫₀^h(u − U)dy [m²/s]; the ground-frame flux C₁ + Uh(x) varies with x), p_max, x_pmax, W_exact, W_linear, err_linear (%), p_max_atm, p_visc, inlet_backflow (bool: α > 1 for U > 0, pad-frame recirculation at the wide end, D14 check)) | E3, IF4 title | §9 (keys fixed here) |
| 2.10 | `slider_gap_velocity(x, y, h0, alpha, L, U, mu=1e-3, frame="ground")` (`"pad"` subtracts U) | (8.19) with the exact slider dp/dx | u [m/s] | C06 figure, E3 | §9 (`frame` added) |
| 2.11 | `hele_shaw_velocity(z, h, grad_p, mu=1e-3)` · `hele_shaw_potential(p, z, h, mu=1e-3)` | N36 $u=-\frac1{2\mu}\frac{\partial p}{\partial x}z(h-z)$, $\phi=-\frac{z(h-z)}{2\mu}p$ | (u, v) [m/s] · φ [m²/s] | C06 (N36) | I18 |
| 2.12 | `thin_film_flux(h, h_x, rho=1000.0, g=G0, mu=1.0)` | Example 8.3 flux $-\frac{\rho g}{3\mu}h^3\frac{\partial h}{\partial x}$ | q [m²/s] | C08, E4 | I19 |
| 2.13 | `thin_film_spread(h0, x, t_out, rho=1000.0, g=G0, mu=1.0, h_min=1e-6, cache=None)` (backward Euler + Picard, precursor film; `cache` = .npz path) | $\frac{\partial h}{\partial t}=\frac{\rho g}{3\mu}\frac{\partial}{\partial x}(h^3\frac{\partial h}{\partial x})$ | dict(t, h (nt, nx), volume (nt,), x_front (nt,)) | C08, A3, C10 (fit) | I19 (dict) |
| 2.14 | `viscous_current_similarity(x, t, area, rho=1000.0, g=G0, mu=1.0, return_front=False)` (`area` = **half** area ∫₀^{x_N}h dx) | N63 Huppert: x_N = η_N(βA³t)^{1/5}, η_N = 1.41124, $h=At^{-1/5}F(x/Dt^{1/5})$ | h [m] (or (h, x_N)) | C08, C10, E4, E6 | I20 (`return_front`) |
| 2.15 | `thin_film_state(t, area, rho=1000.0, g=G0, mu=1.0)` | N63 | dict(x_N, h_centre, beta (ρg/3μ), eta_N, effective_diffusivity (βh_c³ [m²/s])) | E4 | §9 |

### C.3 `fluidpy/core/creeping.py` (NEW module, `core.CRP`; re-exported by `ch08`) — creeping flow for Ch. 8, 13, 16
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 3.1 | `stokes_residual(u_fn, p_fn, x, mu=1e-3, h=1e-4)` | (8.43) ∇p − μ∇²u (central differences) | [Pa/m] (3,) per point | C12, C13 check | I32 |
| 3.2 | `E2(psi, r, theta)` · `E4_residual(psi, r, theta)` (sympy) | N82 E² of (6.77); (8.44) | sympy expr | C13 (D28, D29) | I33 |
| 3.3 | `stokes_sphere_streamfunction(r, theta, U=1.0, a=1.0, frame="body")` | (8.48); fluid frame N93 | ψ [m³/s] (NaN for r < a) | C13, E8 | I34 |
| 3.4 | `stokes_sphere_velocity(r, theta, U=1.0, a=1.0, frame="body")` · `stokes_sphere_velocity_xyz(x, y, z, U=1.0, a=1.0, frame="body")` | (8.49) | (u_r, u_θ) · (u, v, w) [m/s] | C13, C15, E8, A4 | I34 |
| 3.5 | `stokes_sphere_pressure(r, theta, U=1.0, a=1.0, mu=1e-3, p_inf=0.0)` | (8.50) $p-p_\infty=-\frac{3\mu aU\cos\theta}{2r^2}$ | p [Pa] | C14, E9 | I35 |
| 3.6 | `stokes_sphere_surface_stresses(theta, U=1.0, a=1.0, mu=1e-3, p_inf=0.0)` | N89 σ_rr = −p, σ_rθ = −(3μU/2a) sin θ, t_x = σ_rr cos θ − σ_rθ sin θ | (σ_rr, σ_rθ, t_x) [Pa] | C14, E9 | I35 |
| 3.7 | `stokes_drag(mu, a, U, parts=False)` | (8.51) D = 6πμaU | D [N] or dict(pressure, friction, total) | C14, E9 | I36 |
| 3.8 | `stokes_drag_running(theta, mu, a, U)` — drag collected from the rear stagnation point θ = 0 up to θ | D31 in closed form: pressure πμaU(1 − cos³θ), friction πμaU(2 − 3 cos θ + cos³θ) | dict(pressure, friction, total) [N] | C14 figure, E9 sweep | **NEW** |
| 3.9 | `sphere_drag_quadrature(stress_fn, a, n=64)` (Gauss–Legendre in cos θ) | ∮t_x dA | D [N] | C14 (from-scratch comparison) | I36 |
| 3.10 | `stokes_drag_coefficient(Re)` · `oseen_drag_coefficient(Re)` · `proudman_pearson_drag_coefficient(Re)` (Re = 2aU/ν; PP converts to Re_a = Re/2) | (8.52) 24/Re; N99 $\frac{24}{\mathrm{Re}}(1+\frac3{16}\mathrm{Re})$; N101 | C_D | C14, C15, E9 | I37 |
| 3.11 | `terminal_velocity(a, rho_p, rho, mu, g=G0)` (warns if Re > 0.1) · `radius_from_terminal_velocity(U, rho_p, rho, mu, g=G0)` · `millikan_charge(U_fall, U_rise, rho_p, rho, mu, E, g=G0)` | N90 $U_t=\frac{2(\rho'-\rho)ga^2}{9\mu}$; N91 | [m/s] · [m] · q [C] | C14, E9 | I38 |
| 3.12 | `settling_state(a, rho_p, rho, mu, g=G0)` | N90, (8.51), (8.52) | dict(U_t, Re, D, C_D, valid (Re < 0.1), D_pressure, D_friction) | E9 | §9 |
| 3.13 | `inertia_viscous_ratio(r, theta, U=1.0, a=1.0, nu=1e-6)` (central differences on the Cartesian Stokes field, h = 1e-4 r) | D32 \|u·∇u\|/\|ν∇²u\| | [–] | C15, E8 | I39 |
| 3.14 | `oseen_streamfunction(r, theta, U=1.0, a=1.0, Re=0.1, frame="body")` · `oseen_velocity(r, theta, U=1.0, a=1.0, Re=0.1, frame="body")` (analytic derivatives, −expm1 form) | (8.53) | ψ [m³/s] · (u_r, u_θ) [m/s] | C15, E8, IF8 | I40 · §9 |
| 3.15 | `side_line_speed(r, U=1.0, a=1.0, model="stokes", Re=0.1, frame="body")` (`"stokes"` (8.49), `"ideal"` (ch06 sphere: u_θ = −U(1 + a³/2r³) at θ = π/2), `"oseen"` (3.14)) | speed at θ = π/2 | \|u\| [m/s] | C13 worked number, E8 view 2 | **NEW** |
| 3.16 | `stokes_sphere_sympy()` (cached) | D27–D31 | dict(omega_phi, E2psi, E4psi (0), biharmonic_psi (≠ 0), f_ode, roots ([4, 2, 1, −1]), constants (A, B, C, D), psi, u_r, u_theta, curl_omega (r, θ), p, p_check_r (0), p_check_theta (0), sigma_rr_viscous_a (0), sigma_rtheta_a, D_pressure, D_friction, curlcurl_identity (0)) — sympy | C13 (D28, D29 checks), C14 (D30, D31 checks) | I33 (keys fixed here) |

### C.4 `fluidpy/ch08_laminar_flow.py` (chapter module; re-exports C.1–C.3)
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 4.1 | `pipe_flow_regime(U, d, nu)` | C01 Re = Ud/ν; labels "laminar" (< 2000), "transitional" (2000–3000), "turbulent" (> 3000) | (Re, label) | C01, figure | I1 |
| 4.2 | `inertia_viscous_scales(U, L, nu)` · `momentum_diffusivity(fluid, T=293.15)` (`"air"` Sutherland + p/RT, `"water"` IAPWS) · `diffusion_time(L, nu)` | R01 ratio; N01 ν; L²/ν | dict(inertia, viscous, ratio) · ν [m²/s] · t [s] | C01 (D01), R01 | I2 |
| 4.3 | `wall_bc_residuals(u_wall, U_s, n)` | (8.2), (8.3) | (normal, tangential) [m/s] | R04, R05 | I3 |
| 4.4 | `parallel_flow_sympy(case)` (`"channel"`, `"pipe"`, `"circular_couette"`) | D02–D04, D06, D08 | dict with at least: continuity, momentum (reduced equations), separation (d²p/dx² = 0 and u‴ = 0 lines), general (integrated profile with constants), constants, profile, residual (0), bc_residuals (list of 0); pipe adds `unbounded_term` (A ln R); circular Couette adds `A`, `B`, `pressure` | C02, C03, C04 checks | I4 (keys fixed) |
| 4.5 | `advective_acceleration_check(case)` | N22 | sympy vector (0 for channel, pipe; −u_φ²/R e_R) | C04 (N22) | I11 |
| 4.6 | `lubrication_nondim_sympy()` (cached) | D10, D11 | dict(continuity_star, x_coeffs {inertia, pressure, diff_along, diff_across}, y_coeffs {…}, printed_13b_y_coeffs (differ), limit_x (with ν), limit_y, eq_8_18) | C05 (D10 check), E2 | I13 (keys fixed) |
| 4.7 | `reynolds_equation_sympy()` (cached) | D13: Leibniz, kinematic walls, q, Reynolds equation | dict(q, leibniz_residual (0), kinematic_cancellation (0), reynolds_residual (0), steady_q_constant) — sympy | C06 (D13 check) | **NEW** (★★★ check needs a tested engine) |
| 4.8 | `slider_bearing_sympy()` (cached) | D14 | dict(C1_eq, dpdx, antiderivative_3, antiderivative_2, C1, C2, p_exact, ode_residual_exact (0), ode_residual_book (≠ 0), bc_residuals (0, 0), p_linear, W_linear, W_exact, W_exact_series) | C07 (D14 check) | I17 (keys fixed) |
| 4.9 | `hele_shaw_cylinder(x, y, z, U_mean, a, h, mu=1e-3)` (`z=None` → depth average) | N36 with the ideal cylinder | dict(u, v, p, phi) | C06 figure | I18 |
| 4.10 | `similarity_reduce_sympy(case)` (`"stokes1"`, `"stokes1_delta"`, `"vortex_sheet"`, `"line_vortex"`, `"spreading"`) · `similarity_ode_solve(case, eta_max=12.0)` (`solve_bvp`) | D17–D19, D21–D23, N61 | dict(ode, brackets, delta, n, m, bcs, F, residual (0)) · (eta, F) arrays | C09, C10 checks | I21 |
| 4.11 | `vorticity_content(t, U=1.0, nu=1e-6)` | N53 ∫₀^∞ω dy by `quad` | [m/s] (= U) | C09 | I22 |
| 4.12 | `similarity_collapse_error(case, n, m, times)` (cases `"stokes1"`, `"vortex_sheet"`, `"line_vortex"`, `"spreading"`; exact solutions sampled on 201 points of the rescaled variable) | E6 | spread [–] (≈ 0 at the right n, m) | C10 figure, E6 | §9 |
| 4.13 | `stokes_second_sympy()` | D24 | dict(ode, k_roots, bounded_root, f, u, residual (0)) | C11 (D24 check) | I30 |
| 4.14 | `low_re_scaling_sympy()` | D26 | dict(dynamic {inertia, pressure, viscous}, dynamic_times_Re {…}, viscous {…}) — sympy in Re | C12, IF7 | I31 |
| 4.15 | `oseen_linearisation_sympy()` · `oseen_limit_sympy()` | N95; D33 | dict(advection_split, dropped, oseen_x) · dict(s, series, limit_psi, stokes_psi, difference (0)) | C15 | I40 |
| 4.16 | `synthetic_millikan(n_drops=40, seed=0, noise=0.01)` | N91 demo | dict(q [C], e_est [C], n (ints)) | C14 (N91) | I38 |
| 4.17 | `G0`, `G_BOOK`, `P_ATM` (re-exported) | — | constants | throughout | reuse |

### C.5 `scripts/ch08_*.py` (runnable demos, each `--no-show` for headless runs; drawing helpers carry no physics)
`ch08_parallel_flows.py` (Figs. 8.2–8.4) · `ch08_pipe_and_couette.py` (Figs. 8.5–8.7) · `ch08_slider_bearing.py` (Figs. 8.8,
8.9) · `ch08_hele_shaw.py` (Fig. 8.10) · `ch08_spreading.py` (Fig. 8.11) · `ch08_stokes_first.py` (Figs. 8.12–8.13) ·
`ch08_vortex_diffusion.py` (Figs. 8.14–8.15) · `ch08_oscillating_plate.py` (Fig. 8.16) · `ch08_stokes_sphere.py` (Figs.
8.17, 8.19) · `ch08_millikan.py` (Fig. 8.18) · `ch08_oseen.py` (Fig. 8.20) · `ch08_drag_curve.py` — I41. One drawing
helper module **`scripts/ch08_drawings.py`** (NEW, no physics): `channel(ax, h, U)` (hatched fixed wall, moving wall
with arrow), `pipe_section(ax, a)`, `annulus(ax, R1, R2)`, `pad(ax, h0, alpha, L)`, `sphere(ax, a)`, `profile_arrows(ax,
y, u, x0, scale)` — the notebook imports them with `sys.path.insert(0, "scripts")`.

**Flag (functions not in analysis §4):** 1.5, 1.14, 1.19, 1.27, 2.2, 2.9, 2.10, 2.15, 3.12, 3.14 (`oseen_velocity`),
4.12 come from curation §9 (keys fixed above); **3.8 `stokes_drag_running`, 3.15 `side_line_speed`, 4.7
`reynolds_equation_sympy` and `scripts/ch08_drawings.py` are new in this design.** Refinements of §4 signatures: 1.19
gains `rho`; 1.26 drops the callable `envelope`; 1.28 gains `startup_be`, `return_all`; 2.8 returns a dict; 2.9 drops the
ambiguous `q` (C₁ is the pad-frame flux); 2.10 gains `frame`; 2.13 returns a dict and takes `cache`;
2.14 gains `return_front`; 4.4, 4.6, 4.8, 3.16 have their keys fixed; `ch08` re-exports `couette_startup_profile` (E5
preset parity).

---

## Part A — notebook storyboard (`notebooks/build_ch08.py` → `notebooks/ch08_laminar_flow.ipynb`)

**One line per book section** (cell numbers are estimates for the builder's budget; ≈ 560 cells in all):
- §8.1 → R01 R02, C01 (N01, N02, N103 · D01 · primer P185 · figure: dye streak + diffusion times), R03 R04 R05 R06 — cells ≈ 8–45
- §8.2 → R07 R08, C02 (N03–N11 · D02 D03 D04 D05 · IF1 · **E1**), R09, C03 (N12–N16 · D06 D07 · P186 · IF2), R10 R11 R12,
  C04 (N17–N22, N104 · D08 D09 · P187 · IF3) — cells ≈ 46–175
- §8.3 → R13, C05 (N23, N24, N26–N31 · D10 D11 · P188 · **E2**), C06 (N25, N32–N34, N36 · D12 D13 · P189 P190), C07 (N35,
  N105 · D14 · IF4 · live · **E3**), C08 (N37 · D15 · P191 P192 · **A3** · **E4**) — cells ≈ 176–305
- §8.4 → R14, C09 (N38–N56 · D16–D20 · P193 P194 P195 P196 · **A1** · IF5 · **E5**), C10 (N57–N64, N106 · D21 D22 D23 · P197
  · **A5** · **E6**) — cells ≈ 306–420
- §8.5 → C11 (N65–N74 · D24 D25 · **A2** · IF6 · **E7**) — cells ≈ 421–455
- §8.6 → R15 R16, C12 (N75–N80 · D26 · P198 · IF7), R17 R18 R19, C13 (N81–N87, N93, N107 · D27 D28 D29 · P199 · **A4** ·
  **E8**), R20, C14 (N88–N92 · D30 D31 · **E9**), C15 (N94–N101 · D32 D33 · IF8) — cells ≈ 456–550
- §8.7 → N102, S01 pointer, summary — cells ≈ 551–560

Every CORE block below follows the order: problem in plain words → idea → primers → maths (notes and derivations, Part F)
→ tiny example → code (fluidpy) + "What does the code above do?" → from-scratch check → visual(s) → notes and "What would
change if…". Code drafts are intent + exact calls; the builder writes every line commented (novice grade, units, the
equation written out next to its number). *expect* gives the numbers the executed cell must print (sanity values computed
for this design with water μ = 1.0 × 10⁻³ Pa s, ρ = 1000 kg/m³, ν = 10⁻⁶ m²/s unless stated). *see / read / change* are the
three figure notes. Note ids open every note text in bold (**N09 [B]**) so the lesson-reviewer can find them.

### A.0 Front matter
1. `nb.title(big_idea=…, roadmap=[…15…], prerequisites=[…])`. **Big idea (draft):** "Pour honey, drag a spoon through
   syrup, watch a raindrop's tiny cousin — a cloud droplet — drift down at a centimetre a second: when viscosity matters
   everywhere, the Navier–Stokes equations can often be solved exactly. The trick is always the same: a symmetry or a
   thin geometry kills the nonlinear term, and what is left is a balance between pressure (or wall motion) and viscous
   friction. This chapter solves that balance between plates, in pipes and between rotating cylinders; in the thin film
   under a bearing and under a spreading drop; above a plate that starts or shakes, where momentum *diffuses* a distance
   √(νt) — the scale behind every boundary layer, every Ekman layer and every spin-up time of the ocean and atmosphere;
   and round a sphere so small that inertia vanishes, which gives Stokes' drag 6πμaU and the fall speed of cloud
   droplets, aerosols and silt." **Roadmap (one line per CORE):** C01 laminar vs turbulent, ν as a diffusivity · C02
   Couette–Poiseuille flow and backflow · C03 Poiseuille pipe flow, Q ∝ a⁴ · C04 circular Couette flow and its two limits
   · C05 the lubrication approximation · C06 the lubrication profile and the Reynolds equation · C07 the slider bearing ·
   C08 the thin-film (viscous gravity current) equation · C09 Stokes' first problem and similarity · C10 the similarity
   ansatz and its exponents · C11 Stokes' second problem: the Stokes layer · C12 creeping flow: the Stokes equations ·
   C13 Stokes' flow round a sphere · C14 Stokes drag and settling · C15 where Stokes fails, and Oseen's fix.
   **Prerequisites:** Navier–Stokes and the no-slip condition (Ch. 4 §4.5–4.10) · Reynolds number and scaling (Ch. 4
   §4.11) · viscosity, ν = μ/ρ, plane Couette start-up (Ch. 1 §1.5) · dimensional analysis (Ch. 1 §1.11) ·
   vorticity diffusion, the ideal vortex, solid-body rotation (Ch. 5 §5.1–5.2) · the Gaussian vortex (Ch. 3 §3.5) · the
   Stokes stream function and the ideal-flow sphere (Ch. 6 §6.8) · complex amplitudes (Ch. 7 §7.7).
2. `nb.explainer_index([...])` — 9 rows: ("couette_poiseuille_backflow", "When does fluid flow backwards in a channel?",
   "a line plus a parabola: wall drag and pressure add, and backflow opens once dp/dx > 2μU/h²") · ("lubrication_scaling",
   "Why can a thin film ignore inertia?", "term bars of the scaled momentum equations reorder with ε and Re_L; only ε²Re_L
   must be small") · ("slider_bearing", "How does a film thinner than a hair carry a load?", "the same flux through a
   narrowing gap needs a pressure hump — and it turns into suction if the pad slides backwards") ·
   ("viscous_gravity_current", "Why does a honey drop spread like t^(1/5)?", "a nonlinear diffusion that chokes itself;
   every initial shape ends on one profile") · ("stokes_first_problem", "How can profiles at all times be one curve?",
   "momentum diffuses √(νt): plot against y/√(νt) and every profile collapses") · ("similarity_exponents", "Where do
   similarity exponents come from?", "guess At^(−n)F(ξ/δ): the equation and one conserved quantity fix n and δ(t)") ·
   ("oscillating_plate", "How deep does a shaking plate reach?", "a 'wave' made only of diffusion: each layer lags by
   y/δ and shrinks as e^(−y/δ), δ = √(2ν/ω)") · ("stokes_sphere_flow", "What does creeping flow round a sphere look
   like?", "fore–aft symmetric, reaching tens of radii, until inertia returns beyond r ~ a/Re: Oseen's wake") ·
   ("stokes_drag_settling", "Where does 6πμaU come from?", "one third pressure, two thirds friction; drag ∝ U makes
   a droplet's fall speed grow as a²").
3. `nb.setup()`.
4. `nb.code` — **chapter imports** (outside any block): `import numpy as np` · `import matplotlib.pyplot as plt` ·
   `import sympy as sp` · `from scipy import integrate, special, optimize, linalg` · `from fluidpy import
   ch08_laminar_flow as ch08` · `from fluidpy.core import laminar as LAM, lubrication as LUB, creeping as CRP, diffusion`
   · `from fluidpy import ch01_introduction as ch01, ch03_kinematics as ch03, ch04_conservation_laws as ch04,
   ch05_vorticity_dynamics as ch05` · `from fluidpy.core import navier_stokes as NS, similarity as SIM, vortices as VX,
   potential as PF, curvilinear as CL, laplace_solvers as LS` · `from fluidpy.core.interact import slider_figure,
   animate_figure, live` · `from fluidpy.core.anim import animate` · `from fluidpy.core.style import savefig` · `import
   sys; sys.path.insert(0, "scripts"); from ch08_drawings import channel, pipe_section, annulus, pad, sphere,
   profile_arrows` · `MU_W, RHO_W, NU_W = 1.0e-3, 1000.0, 1.0e-6`. *explain:* one line per import ("`ch08` re-exports
   the three new toolkits `core.laminar`, `core.lubrication`, `core.creeping`, so `ch08.channel_flow` and
   `LAM.channel_flow` are the same function").
5. `nb.md` — **⚠️ Conventions in this chapter** (a two-column table "symbol | meaning here (and before)", then numbers):
   dp/dx is the book's pressure gradient (favourable < 0) — Ch. 4 code used G = −dp/dx; channel walls y = 0 (fixed) and
   y = h (moving); R cylindrical radius, r spherical radius; θ in §8.6 from the **downstream** axis (rear stagnation
   point θ = 0); ω = vorticity (§8.1, §8.4, §8.6) but the oscillation frequency in §8.5 (as in Ch. 7); η = y/√(νt)
   (the book's figures plot y/(2√(νt)); here η is the similarity variable, the surface elevation in Ch. 7); ε = h/L the
   fineness ratio (a density ratio in Ch. 7, dissipation in Ch. 4); Λ bearing number; δ₉₉, δ_e, δ(t) three different
   thicknesses; four Reynolds numbers (pipe diameter, passage length, generic L, sphere **diameter** 2a); D is both an
   integration constant and the drag. "We compute with the book's conventions and name the other one wherever it
   differs." A second table lists the book slips of convention 5 with where each is taught.
6. `nb.md` — **🔁 Tools from earlier chapters used in this one** (one line each, primer number): partial derivative
   (P25), first-order Taylor (P26), definite integral (P27), linear second-order ODE (P44), separation of variables (P42),
   chain rule (P49/P91), exponent rules (P43), power laws and log–log plots (P13), limits and orders of smallness (P68),
   order-of-magnitude scaling (P130), scaled variables and the chain rule (P133), product rule (P38), fundamental theorem
   of calculus (P84), differentiation under the integral sign (P109), substitution in an integral (P106), complementary
   error function erfc (P123), square root of a negative number (P45), complex square roots (P159), complex amplitudes
   (P176), curl of a curl (P122), Schwarz's theorem (P121), surface integrals on a sphere (P163), line integral along a
   path (P35), cylindrical and spherical unit vectors (P88), polar coordinates as a moving basis (P105), boundary
   conditions (P20), np.linalg.solve (P57), sympy (P40), sympy expand/series (P117), quad (P87), brentq (P108), expm1
   (P107), minimize_scalar (P170), np.interp (P182), meshgrid (P76), broadcasting (P77), contour/streamplot (P78), finite
   differences (P21), explicit stepping (P30), trapezoid rule (P37), np.random.default_rng (P10), animate (P16),
   slider_figure (P17), show_viz (P18), live widgets (P47), assert np.allclose (P15). Each block repeats the ones it uses
   in a one-line reminder at first use.

### A.1 §8.1 Introduction — R01, R02, C01 (+N01, N02, N103 · D01 · P185), R03–R06
1. `nb.section("8.1", "Introduction", intro="**What is this section about?** Where this chapter lives: flows in which
   viscosity acts everywhere, not just in thin layers. We meet Reynolds's experiment (smooth versus chaotic flow in a
   pipe), learn to read the kinematic viscosity ν as a *diffusivity* of momentum and vorticity, and restate the equations
   and wall conditions every solution in the chapter must satisfy.")`
2. `nb.recap("R01", "The Reynolds number as inertia over viscosity", "Scaling the Navier–Stokes equation (Ch. 4) gives
   inertia ~ U²/L and viscous acceleration ~ μU/(ρL²); their ratio $\\frac{U^2/L}{\\mu U/\\rho L^2}=\\frac{\\rho UL}{\\mu}=
   \\mathrm{Re}$ is the Reynolds number. Earlier chapters dropped viscosity where Re ≫ 1 (ideal flow, Ch. 6); in this
   chapter viscosity is kept everywhere.", where="Ch. 4 §4.11")` + `nb.code`: `print(ch08.inertia_viscous_scales(0.1,
   0.01, NU_W))` · `print(SIM.reynolds_number(0.1, 0.01, nu=NU_W))`. *expect:* inertia 1.0 m/s², viscous 1.0e-3 m/s²,
   ratio 1000.0; Re = 1000.0.
3. `nb.recap("R02", "Vorticity diffuses with ν", "In a plane flow the vorticity equation (5.13) loses its stretching term,
   because ω = ω_z e_z and nothing depends on z, so (ω·∇)u = ω_z ∂u/∂z = 0; what is left is $D\\omega_z/Dt=\\nu\\nabla^2
   \\omega_z$ — the same form as the heat equation. A fluid with larger ν smooths out a vorticity pattern faster.",
   where="Ch. 5 §5.4")`
4. `nb.core("C01", "Laminar vs turbulent flow, and ν as a momentum diffusivity: $\\mathrm{Re}=Ud/\\nu\\sim2000\\text{–}3000$",
   question="When does a flow stay smooth — and what does the kinematic viscosity ν actually measure?")`
5. `nb.md` — **The problem in plain words:** "Open a tap a little and the stream is glassy; open it fully and it turns
   frothy. Reynolds put a thread of dye into water flowing through a glass tube and saw the same thing: at low flow the
   dye stayed a straight thread, at high flow it broke up and filled the tube. Every exact solution in this chapter is of
   the first, *laminar* kind — so we need to know when it applies, and what sets how fast viscosity spreads the motion."
6. `nb.md` — **The idea** (ASCII + table):
   ```
   laminar  (Re < 2000):   ─────────────── dye stays a thread; layers slide past each other
   turbulent (Re > 3000):  ~~~≈≈≈∿∿∿≈≈~~~ dye fills the tube; eddies mix it
   Re = U d / ν   (U = mean speed, d = diameter, ν = μ/ρ)
   ```
   | quantity | water (20 °C) | air (20 °C) |
   |---|---|---|
   | μ [Pa s] | 1.0 × 10⁻³ | 1.8 × 10⁻⁵ |
   | ρ [kg/m³] | 1000 | 1.2 |
   | ν = μ/ρ [m²/s] | 1.0 × 10⁻⁶ | 1.5 × 10⁻⁵ |
   | time to diffuse 1 cm, L²/ν | 100 s | 6.7 s |
   "**ν, not μ, says how fast motion spreads** — air is 55 times less viscous than water but spreads momentum 15 times
   faster."
7. `nb.md` — "> ⚠️ **Common confusion:** μ (dynamic viscosity, Pa s) is the friction coefficient in τ = μ du/dy; ν = μ/ρ
   (kinematic viscosity, m²/s) is the diffusivity. The book's §8.1 once calls μ 'the kinematic viscosity' — it means the
   dynamic one."
8. `nb.primer("diffusivity and the diffusion time L²/ν", "A diffusivity D (units m²/s) says how fast something spreads by
   random molecular exchange: in a time t it spreads a distance of order √(Dt), so crossing a distance L takes about
   L²/D. Heat has κ = k/ρC_p, dye its molecular diffusivity, and momentum has ν. Double the distance and the time
   quadruples.", code="L = [0.001, 0.01, 0.1]            # distances [m]: 1 mm, 1 cm, 10 cm\nnu = 1.0e-6        # water's
   kinematic viscosity [m^2/s]\nprint([l**2/nu for l in L])      # diffusion times [s]: [1.0, 100.0, 10000.0]")` (**P185**)
9. `nb.note` — **N01 [B]** "**Heat and momentum spread the same way.** The Boussinesq heat equation $DT/Dt=\\kappa
   \\nabla^2T$ (4.89), with the thermal diffusivity κ ≡ k/ρC_p, has the same form as $D\\mathbf u/Dt=-(1/\\rho)\\nabla p+
   \\nu\\nabla^2\\mathbf u$ (8.1); both κ and ν are in m²/s, so ν is the *momentum diffusivity*. The analogy is not
   complete: the pressure gradient in (8.1) also moves momentum around, so velocity is not simply diffused." equation
   `DT/Dt=\kappa\nabla^2T`, ref "4.89".
10. `nb.derivation("D01", …)` — Part F D01 (5 steps), ref "".
11. `nb.worked_example("water in a 1 cm tube at 10 cm/s", "1. $\\mathrm{Re}=Ud/\\nu=0.1\\times0.01/10^{-6}=1000$ —
    below 2000, laminar. 2. Double the speed: Re = 2000, the bottom of the transition band. 3. Diffusion across the
    tube: $t=d^2/\\nu=10^{-4}/10^{-6}=100$ s — much longer than the 0.1 s the water needs to travel 1 cm. 4. Same tube
    with air at 10 cm/s: $\\mathrm{Re}=0.1\\times0.01/1.5\\times10^{-5}=67$, $t=6.7$ s.")`
12. `nb.code` — `for U in (0.1, 0.25, 0.5): print(U, ch08.pipe_flow_regime(U, 0.01, NU_W))` · `nu_air =
    ch08.momentum_diffusivity("air", 293.15)` · `nu_water = ch08.momentum_diffusivity("water", 293.15)` · `print(nu_air,
    nu_water, nu_air/nu_water)` · `print(ch08.diffusion_time(0.01, nu_water), ch08.diffusion_time(0.01, nu_air))`.
    *expect:* (1000, 'laminar'), (2500, 'transitional'), (5000, 'turbulent'); ν_air ≈ 1.5 × 10⁻⁵ m²/s, ν_water ≈ 1.0 ×
    10⁻⁶ m²/s, ratio ≈ 15 (the property functions give ≈ 15.1); times ≈ 99.6 s and ≈ 6.6 s. *explain:* 1. `pipe_flow_regime`
    computes Re = Ud/ν and applies the 2000–3000 band; 2. `momentum_diffusivity` takes μ and ρ from Ch. 1's property
    functions (Sutherland air, IAPWS water) and returns ν = μ/ρ; 3. `diffusion_time` is L²/ν.
13. `nb.check_agree` — **from scratch (curation §7):** `Re = U*d/nu` and `label = "laminar" if Re < 2000 else
    ("turbulent" if Re > 3000 else "transitional")` for the three speeds; `assert np.allclose([Re…],
    [ch08.pipe_flow_regime(U, 0.01, NU_W)[0] for U in …], rtol=1e-12)`; `assert np.allclose(0.01**2/nu_air,
    ch08.diffusion_time(0.01, nu_air), rtol=1e-12)`.
14. `nb.figure` — **N103 [B] (Fig. 8.1 replaced by our synthetic sketch)** (7 × 3.2 in, two panels): (a) a tube drawn
    from x = 0 to 1 m, a dye thread injected at the axis; for Re = 1000 a straight thread (blue), for Re = 5000 the thread
    breaks into a seeded random walk that fills the tube (labelled "schematic — not a simulation", the regime from
    `ch08.pipe_flow_regime`); (b) diffusion time L²/ν vs L on log–log axes (L from 10 µm to 1 m) for water (blue) and air
    (teal), slope-2 guide, the 1 cm point marked on each line. Title "Smooth below Re ≈ 2000; momentum spreads as √(νt)".
    *see:* "a straight dye thread at Re = 1000 and a mixed tube at Re = 5000; two parallel lines of slope 2, air's 15 times
    lower." *read:* "on (b) read the time to diffuse a distance: 1 cm takes 100 s in water, 6.7 s in air; the slope 2 says
    ten times farther takes a hundred times longer." *change:* "…the fluid were glycerine (ν ≈ 10⁻³ m²/s): the line drops
    by 1000 — 1 cm in 0.1 s — and Re at 10 cm/s is 1, deeply laminar."
15. `nb.note` — **N02 [C]** "**Where this leads.** When viscosity matters only in thin layers next to walls, we get
    boundary layers (Ch. 9); whether these laminar flows survive small disturbances is Ch. 11's question; fully turbulent
    flow is Ch. 12; viscous layers in a rotating fluid — the Ekman layers of the ocean and atmosphere — are Ch. 13. The
    √(νt) spreading of this chapter reappears in all four."
16. `nb.recap("R03", "The Navier–Stokes equation for constant ρ and μ", "Every solution below satisfies continuity
    ∇·u = 0 (4.10) and $D\\mathbf u/Dt=-(1/\\rho)\\nabla p+\\nu\\nabla^2\\mathbf u$ (8.1), which is (4.85). We check each one
    with Ch. 4's term calculator: plug the velocity and pressure in and the residual must vanish.", where="Ch. 4 §4.6")` +
    `nb.code`: `u_fn = lambda x, t=0.0: np.stack([ch08.channel_flow(x[1], 0.01, U=0.1, dpdx=-4.0), 0*x[0], 0*x[0]])` (the channel velocity as a function of position, component-first layout as in Ch. 4: x[0] = x, x[1] = y) · `p_fn = lambda x, t=0.0: -4.0 * x[0]` (linear pressure, dp/dx = −4 Pa/m) · `terms = NS.ns_incompressible_terms(u_fn, p_fn, np.array([0.3, 0.004, 0.0]), mu=MU_W, g=(0.0, 0.0, 0.0))` · `print(terms)`. *expect:* the x-residual ≈ 0 (< 1e-9 relative to the pressure term 4 × 10⁻³ m/s²); advective term
    exactly 0. *explain:* the pressure term −(1/ρ)∂p/∂x = +4 × 10⁻³ m/s² and the viscous term ν d²u/dy² = −4 × 10⁻³ m/s²
    cancel; the advective term is zero (C02 explains why).
17. `nb.recap("R04", "No flow through a wall", "At a solid surface moving with velocity U_s the normal velocities agree:
    $\\mathbf n\\cdot\\mathbf U_s=(\\mathbf n\\cdot\\mathbf u)_{\\text{on the surface}}$ (8.2).", where="Ch. 4 §4.10")`
18. `nb.recap("R05", "No slip", "The tangential velocities agree too: $\\mathbf t\\cdot\\mathbf U_s=(\\mathbf t\\cdot
    \\mathbf u)_{\\text{on the surface}}$ (8.3) — viscous fluid sticks to walls. Every constant of integration in this
    chapter is fixed by (8.2) and (8.3).", where="Ch. 4 §4.10")` + `nb.code`: `print(ch08.wall_bc_residuals(np.array([0.1,
    0, 0]), np.array([0.1, 0, 0]), np.array([0, 1.0, 0])))` · `print(ch08.wall_bc_residuals(np.array([0.08, 0.01, 0]),
    np.array([0.1, 0, 0]), np.array([0, 1.0, 0])))`. *expect:* (0.0, 0.0); (0.01, −0.02) m/s — a leaky, slipping wall.
19. `nb.recap("R06", "Standing assumptions", "Constant density, an inertial frame, and gravity absorbed into the pressure
    (p → p + ρgz) whenever there is no free surface — only Example 8.3 (C08) and the settling sphere (C14) bring g
    back.", where="Ch. 4 §4.9")`

### A.2 §8.2 Exact Solutions for Steady Incompressible Viscous Flow — R07, R08, C02 (+N03–N11 · D02–D05 · IF1 · E1), R09, C03 (+N12–N16 · D06 D07 · P186 · IF2), R10, R11, R12, C04 (+N17–N22, N104 · D08 D09 · P187 · IF3)
1. `nb.section("8.2", "Exact Solutions for Steady Incompressible Viscous Flow", intro="**What is this section about?**
   Three flows where symmetry removes the nonlinear term u·∇u exactly, so the Navier–Stokes equation becomes a linear
   ODE that we can integrate by hand: flow between parallel plates driven by a moving wall and a pressure gradient, flow
   in a round pipe, and flow between two rotating cylinders. The same short argument — continuity kills one velocity,
   'a function of x equals a function of y' makes the pressure gradient constant — solves all three.")`
2. `nb.recap("R07", "Plane Couette flow", "A fluid between a fixed plate and a plate sliding at U, with no pressure
   gradient, moves with the straight profile u = Uy/h and a uniform shear stress τ = μ du/dy = μU/h (Ch. 1). C02 shows it
   is one half of a two-part formula.", where="Ch. 1 §1.5")`
3. `nb.recap("R08", "Plane Poiseuille flow", "With both plates at rest and a pressure gradient, the profile is the
   parabola $u(y)=-\\frac1{2\\mu}\\frac{dp}{dx}y(h-y)$, the stress $\\tau=\\mu\\frac{du}{dy}=-(\\frac h2-y)\\frac{dp}{dx}$
   is linear and its magnitude at both walls is (h/2)∣dp/dx∣. Ch. 4's `plane_poiseuille` takes G = −dp/dx.",
   where="Ch. 4 §4.6")` + `nb.code`: `y = np.linspace(0, 0.01, 5)` · `print(ch04.plane_poiseuille(y, G=2.0, h=0.01))` ·
   `print(ch08.channel_flow(y, 0.01, U=0.0, dpdx=-2.0))`. *expect:* identical arrays `[0, 0.01875, 0.025, 0.01875, 0]` m/s.
   Markdown: "> ⚠️ **Common confusion:** the book's dp/dx = −2 Pa/m (pressure falling downstream) is Ch. 4's G = +2 Pa/m.
   Same flow, opposite sign of the argument — the parity above pins it."
4. `nb.core("C02", "Plane Couette–Poiseuille flow $u(y)=\\frac Uh y-\\frac1{2\\mu}\\frac{dp}{dx}y(h-y)$ (8.5)",
   question="What profile results from a moving wall plus a pressure gradient — and when does fluid near the fixed wall
   flow backwards?")`
5. `nb.md` — **The problem in plain words:** "A belt drags oil along a channel while a pump pushes back. Near the belt the
   oil is carried forward; near the fixed floor the pump may win and push it backwards. Journal bearings, extrusion dies,
   the wind-driven surface layer of a lake with a return current underneath — all are this flow. We want the exact
   profile and the condition for the reversed current."
6. `nb.md` — **The idea** (ASCII): "Because the equation turns out linear, the two causes simply **add**:"
   ```
   y=h ═══════════► U        y=h ═════════        y=h ═══════════► U
        ─────────►                 ──►                    ─────────►
        ──────►        +           ────►        =         ──────►
        ───►                       ──►                    ─►      (backflow if the parabola
   y=0 ═════════             y=0 ═════════        y=0 ═══◄══       is pushed the other way)
     Couette: line            Poiseuille: parabola          sum (8.5)
   ```
7. `nb.note` — **N03 [B]** "**Fully developed flow.** Near the channel inlet, boundary layers grow on both walls and the
   profile changes with x — inside this *entrance length* ∂u/∂x ≠ 0, so continuity ∂u/∂x + ∂v/∂y = 0 forces v ≠ 0 and the
   flow is not parallel. Downstream the layers merge and the profile stops changing: u = u(y) alone. Our figure below draws
   the wall-layer edge with the diffusion estimate $\\delta_{99}\\sim3.64\\sqrt{\\nu t}$ (8.31) and t = x/U — our estimate
   (C09 derives it), not the book's; Ch. 9 treats entry flow properly." + `nb.figure`: a channel h = 1 cm, U = 1 cm/s
   uniform inflow, the two wall-layer edges δ(x) = `ch08.diffusion_thickness(x/0.01, NU_W)` (dashed, labelled "our
   estimate") meeting at x_e where δ = h/2 (printed: x_e ≈ 0.019 m, about twice the gap), profiles sketched at three stations (uniform-ish,
   developing, parabola from `channel_flow`). *see:* "two dashed edges growing like √x and meeting." *read:* "left of the
   meeting point the core still moves uniformly; right of it the profile is the parabola of (8.5) with U = 0." *change:*
   "…the inflow were 10 times faster: the edge δ ∝ √(νx/U) grows √10 times more slowly, so the entrance length
   grows tenfold (x_e ∝ Uh²/ν)."
   (Builder: compute x_e from δ(x_e) = h/2 by `brentq` and print it — x_e = Uh²/(4 × 3.643² ν) = 0.0188 m; never type it.)
8. `nb.note` — **N04 [B]** "**Continuity kills v.** With ∂u/∂x = 0, continuity gives ∂v/∂y = 0; v is then the same at
   every height and equals its wall value 0: u = (u(y), 0, 0). This is steps 1–3 of the derivation below." **N05 [B]**,
   **N06 [B]** "What is left of the momentum equations is $0=-\\frac1\\rho\\frac{\\partial p}{\\partial x}+\\nu
   \\frac{d^2u}{dy^2}$ (8.4a) and $0=-\\frac1\\rho\\frac{\\partial p}{\\partial y}$ (8.4b): no acceleration at all — the
   advective term vanishes *exactly*, which is why the problem is linear." equations (8.4a), (8.4b).
9. `nb.derivation("D02", …)` — Part F D02 (7 steps), ref "8.4".
10. `nb.note` — **N07 [B]** "**Separation of variables in one line.** From (8.4b) p depends on x only; then
    $\\frac1\\mu\\frac{dp}{dx}=\\frac{d^2u}{dy^2}$ says 'a function of x = a function of y', so both equal one constant:
    the pressure falls **linearly** along the channel. The same argument returns in the pipe (C03), the rotating
    cylinders (C04) and the lubrication gap (C06)." equation `\frac{1}{\mu}\frac{dp}{dx}=\frac{d^2u}{dy^2}=\text{const}`.
11. `nb.derivation("D03", …)` — Part F D03 (5 steps), ref "8.4".
12. `nb.md` — gloss **integrating an ODE twice** (not a new primer): "Integrating u″ = c twice brings two constants,
    u = cy²/2 + c₁y + c₂; two wall conditions fix them (a 2 × 2 linear system, P57). sympy's `dsolve` does it in one
    line:" + `nb.code` (3 lines): `y, c = sp.symbols("y c")` · `u = sp.Function("u")` · `print(sp.dsolve(sp.Eq(u(y).diff(y,
    2), c)))` *expect:* `Eq(u(y), C1 + C2*y + c*y**2/2)`.
13. `nb.note` — **N08 [B]** "**The book's constants.** The book integrates to $0=-\\frac{y^2}{2}\\frac{dp}{dx}+\\mu u+Ay+B$
    with B = 0 and A = (h/2)(dp/dx) − μU/h; its A and B are minus our c₁ and c₂ (D04 step 3)." equation
    `0=-\frac{y^2}{2}\frac{dp}{dx}+\mu u+Ay+B`.
14. `nb.derivation("D04", …)` — Part F D04 (9 steps), ref "8.5".
15. `nb.worked_example("a 1 cm channel, belt at 10 cm/s, adverse gradient 4 Pa/m", "Take h = 0.01 m, U = 0.1 m/s,
    μ = 10⁻³ Pa s, dp/dx = +4 Pa/m (pressure rising downstream), at the quarter height y = h/4 = 0.0025 m. 1. Couette
    part: $Uy/h=0.1\\times0.25=0.025$ m/s. 2. Poiseuille part: $-\\frac{1}{2\\mu}\\frac{dp}{dx}y(h-y)=-\\frac{4}{2\\times
    10^{-3}}\\times0.0025\\times0.0075=-0.0375$ m/s. 3. Sum: u = −0.0125 m/s — the fluid moves **backwards** at the quarter
    height. 4. Threshold: $2\\mu U/h^2=2\\times10^{-3}\\times0.1/10^{-4}=2$ Pa/m; we are at twice it. 5. Flow rate
    $Q=\\frac{Uh}2-\\frac{h^3}{12\\mu}\\frac{dp}{dx}=5\\times10^{-4}-3.33\\times10^{-4}=1.67\\times10^{-4}$ m²/s: still
    forward overall.")`
16. `nb.code` — `h, U = 0.01, 0.1` · `y = np.linspace(0, h, 9)` · `u = ch08.channel_flow(y, h, U=U, dpdx=4.0)` ·
    `print(np.round(u, 4))` · `print(ch08.channel_flow_rate(h, U=U, dpdx=4.0))` · `print(ch08.channel_backflow_threshold(U,
    h))` · `s = ch08.couette_poiseuille_state(h, U, 4.0)` · `print({k: s[k] for k in ("Q_couette", "Q_poiseuille",
    "tau_bottom", "tau_top", "backflow", "y_reversal")})` · `print(ch08.parallel_flow_sympy("channel")["profile"])`.
    *expect:* u at y = h/4 is −0.0125; (Q, V) = (1.667e-4 m²/s, 0.01667 m/s); threshold 2.0 Pa/m; Q_couette 5.0e-4,
    Q_poiseuille −3.33e-4, τ_bottom −0.01 Pa, τ_top +0.03 Pa, backflow True, y_reversal 0.005 m (the reversed layer fills
    the lower half); the sympy profile equals (8.5). *explain:* 1. `channel_flow` evaluates (8.5); 2. `channel_flow_rate`
    integrates it exactly; 3. `couette_poiseuille_state` bundles Q split into its two parts, the wall stresses, the
    threshold and where u changes sign; 4. `parallel_flow_sympy` repeats D02–D04 symbolically.
17. `nb.check_agree` — **from scratch (curation §7):** finite-difference solve of μu″ = dp/dx on 101 points: `N = 101; yy
    = np.linspace(0, h, N); dy = yy[1] - yy[0]`; build the tridiagonal matrix `Amat = (np.diag(-2*np.ones(N-2)) +
    np.diag(np.ones(N-3), 1) + np.diag(np.ones(N-3), -1)) / dy**2`; right side `rhs = np.full(N-2, 4.0/MU_W)` with the wall
    values moved over (`rhs[-1] -= U/dy**2`); `u_in = np.linalg.solve(Amat, rhs)`; `u_fd = np.r_[0.0, u_in, U]`; `assert
    np.allclose(u_fd, ch08.channel_flow(yy, h, U=U, dpdx=4.0), rtol=1e-10, atol=1e-12)`. Markdown: "Central differences
    are exact for a parabola, so the hand-built solver reproduces (8.5) to round-off."
18. `nb.note` — **N09 [B]** "**Favourable, adverse, backflow (Fig. 8.4).** dp/dx < 0 pushes with the wall: a fuller
    profile. dp/dx > 0 pushes against it; the fluid next to the fixed wall reverses as soon as the slope du/dy at y = 0
    turns negative, i.e. $\\frac{dp}{dx}>\\frac{2\\mu U}{h^2}$ (ours — the book only draws it). Ch. 9 meets the same
    fight between wall shear and an adverse gradient: that is how boundary layers separate." equation
    `\frac{dp}{dx}>\frac{2\mu U}{h^{2}}`.
19. `nb.derivation("D05", …)` — Part F D05 (7 steps), ref "8.5".
20. `nb.note` — **N10 [B]** "**Flow rate and mean velocity.** Two polynomial integrals of (8.5) give
    $Q=\\int_0^h u\\,dy=\\frac{Uh}{2}\\Big[1-\\frac{h^2}{6\\mu U}\\frac{dp}{dx}\\Big]$ and $V\\equiv Q/h=\\frac U2\\Big[1-
    \\frac{h^2}{6\\mu U}\\frac{dp}{dx}\\Big]$ — a favourable gradient raises Q, an adverse one lowers it, and Q = 0 at
    dp/dx = 6μU/h² (6 Pa/m above). ⚠️ The book writes the middle expression as $V=\\int_0^hu\\,dy$ without the 1/h; the
    units give it away (m²/s is not m/s)." equation (Q, V) + `nb.code`: `Q_num = integrate.quad(lambda yy:
    ch08.channel_flow(yy, h, U=U, dpdx=4.0), 0, h)[0]` · `print(Q_num, ch08.channel_flow_rate(h, U=U, dpdx=4.0)[0])`
    *expect:* both 1.6667e-4 m²/s.
21. `nb.figure` — **Fig. 8.4 remake** (9 × 3.4 in, two panels): (a) u(y)/U vs y/h for dp/dx = −4, 0, 2, 4, 8 Pa/m (h = 1 cm,
    U = 0.1 m/s; colours from blue to amber; the backflow part of the dp/dx = 4 and 8 curves shaded amber; the pure
    Poiseuille profile U = 0, dp/dx = −4 dashed orange for comparison); (b) the signed stress τ(y) from
    `channel_shear_stress` for the same cases (linear lines; Couette vertical at τ = μU/h = 0.01 Pa). Title "Adverse
    pressure beats the wall near the floor once dp/dx > 2μU/h²". *see:* "five profiles through (0, 0) and (1, 1); two of
    them dip below zero near the floor; the stresses are straight lines." *read:* "the slope of a profile at y = 0 is
    proportional to the floor stress; it crosses zero exactly at the 2 Pa/m curve (vertical at the floor), matching
    panel (b), where τ(0) = 0 for that case." *change:* "…the gap were halved to 5 mm: the threshold 2μU/h² quadruples to
    8 Pa/m, so the 4 Pa/m case no longer reverses."
22. `nb.plotly` — **IF1** `slider_figure(lambda G: {"sum u(y)": (u_sum(G), y/h), "Couette part": (U*y/h, y/h), "Poiseuille
    part": (u_pois(G), y/h)}, "dp/dx", np.linspace(-8, 12, 21), unit="Pa/m", xlabel="u [m/s]", ylabel="y/h [–]",
    title="A line plus a parabola")` with u from `ch08.channel_flow`; the trace name of the sum includes the verdict
    ("backflow" when `couette_poiseuille_state(...)["backflow"]`). *explain:* the slider re-evaluates (8.5); watch the blue
    sum cross u = 0 at the floor at dp/dx = 2 Pa/m.
23. `nb.explainer("couette_poiseuille_backflow", heading="When does fluid flow backwards in a channel?", why="Drag the
    pressure gradient through zero and past 2μU/h²: the straight Couette part and the parabolic Poiseuille part add up in
    front of you, the reversed layer opens at the floor and the flow rate splits into its two contributions — a family
    of curves that a static figure can only sample.", tries=["Start at pure Couette (dp/dx = 0) and slowly raise dp/dx —
    watch the floor slope readout reach 0 at 2 Pa/m.", "Press the 'zero net flow' preset: the Couette and Poiseuille bars
    cancel exactly at dp/dx = 6μU/h².", "Click at mid-height in the profile view and read the arithmetic of u from its two
    parts.", "Halve the gap h and predict the new threshold before you look."])`
24. `nb.note` — **N11 [C]** "A constant pressure gradient and a *linear* stress τ(y) are general for any fully developed
    channel flow — they survive in turbulent flow for the time-averaged stress (Ch. 12 uses exactly this)."
25. `nb.md` — **What would change if…** "…the channel were a round pipe? The same three moves (continuity, 'function of x
    = function of r', integrate twice) work, but the Laplacian brings a 1/R — and a new rule: the solution must stay
    finite on the axis (C03)."
26. `nb.recap("R09", "Shear stress in cylindrical coordinates", "In cylindrical coordinates the shear stress on a surface
    of constant R in the z-direction is $\\tau_{zR}=\\mu(\\frac{\\partial u_R}{\\partial z}+\\frac{\\partial u_z}{\\partial
    R})$ (Ch. 4's `strain_rate(..., 'cylindrical')`).", where="Ch. 4 §4.5")`
27. `nb.core("C03", "Poiseuille pipe flow $u_z(R)=\\frac{R^2-a^2}{4\\mu}\\frac{dp}{dz}$ (8.6) and Hagen–Poiseuille",
    question="How much flows through a pipe for a given pressure drop — and why does halving the radius cut the flow
    sixteen-fold?")`
28. `nb.md` — **The problem in plain words:** "Blood squeezing through a capillary, water in a garden hose, oil in a
    pipeline, air in the lungs' smallest airways: a pressure difference pushes fluid through a tube and viscosity at the
    wall holds it back. Doctors care because a slightly narrowed artery needs a much larger pressure; engineers because
    pumping costs follow the same law."
29. `nb.md` — **The idea** (ASCII): "a paraboloid of velocity, fastest on the axis, zero at the wall; the stress grows
    linearly from the axis to the wall:"
   ```
   R=a  ───────────────  u=0,  |τ| max = τ₀
          ─────────►
            ──────────────►    u_max = 2V on the axis
          ─────────►
   R=a  ───────────────
   ```
30. `nb.note` — **N12 [B]** "**Set-up in (R, φ, z).** Take u = (0, 0, u_z(R)): continuity is automatic, the R- and
    φ-equations give $0=\\partial p/\\partial\\varphi$ and $0=\\partial p/\\partial R$, so p = p(z). The z-equation is
    $0=-\\frac{dp}{dz}+\\frac\\mu R\\frac{d}{dR}\\big(R\\frac{du_z}{dR}\\big)$." equation (z-momentum).
31. `nb.primer("Laplacian in cylindrical coordinates", "For a function of R only, ∇²u = (1/R) d/dR(R du/dR) — the R inside
    the bracket comes from the circle's circumference growing with R. For an azimuthal velocity u_φ(R) the φ-component of
    the vector Laplacian is d/dR[(1/R) d(R u_φ)/dR] (it contains an extra −u_φ/R² because the unit vector e_φ turns). Both
    come from `core.curvilinear`.", code="import sympy as sp\nR = sp.symbols('R', positive=True)\nu = R**2              #
    a test profile\nprint(sp.simplify(sp.diff(R*sp.diff(u, R), R)/R))   # (1/R)(R u')' = 4 for u = R^2")` (**P186**)
32. `nb.md` — gloss **regularity at the axis** and **area element**: "ln R → −∞ as R → 0 (P36), so a solution containing
    A ln R is infinite on the axis; a physical velocity is finite, so A = 0. A thin ring of radius R and width dR has area
    2πR dR (P83), which turns a flow rate into ∫u 2πR dR."
33. `nb.derivation("D06", …)` — Part F D06 (10 steps), ref "8.6". **N13 [B]** is its step 7–8 ("$u_z=\\frac{R^2}{4\\mu}
    \\frac{dp}{dz}+A\\ln R+B$; A = 0 by boundedness, B from no slip") — tag it in a one-line note after the derivation:
    `nb.note("**N13 [B]** The general solution $u_z(R)=\\frac{R^2}{4\\mu}\\frac{dp}{dz}+A\\ln R+B$ keeps two constants; the
    axis removes A, the wall fixes B = −(a²/4μ)dp/dz (steps 7–9 above).")`.
34. `nb.note` — **N14 [B]** "**Linear stress.** Differentiating (8.6): $\\tau=\\mu\\frac{\\partial u_z}{\\partial R}=
    \\frac R2\\frac{dp}{dz}$ (8.7) — zero on the axis, largest at the wall." **N15 [B]** "**Wall stress** $\\tau_0=\\frac
    a2\\frac{dp}{dz}$ (8.8). It is negative for forward flow (dp/dz < 0): the wall pulls the fluid back. A force balance on
    a slug of fluid gives the same line without the profile, which is why (8.8) also holds for averaged turbulent pipe flow
    (D07 steps 4–5; Ch. 12 builds the friction velocity on it)." **N16 [B]** "**Hagen–Poiseuille.** $Q=\\int_0^au\\,2\\pi
    R\\,dR=-\\frac{\\pi a^4}{8\\mu}\\frac{dp}{dz}$, $V=\\frac{Q}{\\pi a^2}=-\\frac{a^2}{8\\mu}\\frac{dp}{dz}$, u_max = 2V, and
    the Darcy friction factor f = 64/Re (ours)." equations (8.7), (8.8), Q, V.
35. `nb.derivation("D07", …)` — Part F D07 (10 steps), ref "8.8".
36. `nb.worked_example("a 2 mm tube at dp/dz = −1000 Pa/m", "a = 1 mm, μ = 10⁻³ Pa s, water. 1. $V=-\\frac{a^2}{8\\mu}
    \\frac{dp}{dz}=\\frac{10^{-6}}{8\\times10^{-3}}\\times1000=0.125$ m/s. 2. u_max = 2V = 0.25 m/s. 3.
    $Q=\\pi a^2V=\\pi\\times10^{-6}\\times0.125=3.93\\times10^{-7}$ m³/s (0.39 mL/s). 4. $\\tau_0=\\frac a2\\frac{dp}{dz}=
    -0.5$ Pa. 5. Re = V·2a/ν = 0.125 × 0.002/10⁻⁶ = 250: laminar ✓. 6. f = 64/250 = 0.256, and the definition
    8∣τ₀∣/(ρV²) = 8 × 0.5/(1000 × 0.0156) = 0.256 ✓. 7. Halve a: V falls 4×, the area 4×, so Q falls 16×.")`
37. `nb.code` — `a, G = 1e-3, -1000.0` · `print(ch08.pipe_flow_rate(a, G))` · `print(ch08.pipe_wall_stress(a, G),
    ch08.pipe_shear_stress(0.5e-3, G))` · `Q, V, umax = ch08.pipe_flow_rate(a, G)` · `Re = V*2*a/NU_W` ·
    `print(Re, ch08.pipe_friction_factor(Re))` · `R = np.linspace(0, a, 5)` · `print(ch08.pipe_poiseuille(R, a, G))` ·
    `print(ch03.pipe_profile(R, z=1e3, U_mean=V, R=a))` (the developed profile of Ch. 3) · `print(NS.exact_solution(
    "pipe_poiseuille", np.array([[0.0, 0.0, 0.0]]), G=1000.0, R=a)[0])`. *expect:* (3.927e-7 m³/s, 0.125, 0.25); −0.5,
    −0.25 Pa; 250, 0.256; profile [0.25, 0.234, 0.1875, 0.109, 0] m/s, the same from ch03 (z far downstream) and 0.25 on
    the axis from ch04. *explain:* three independent earlier implementations agree with (8.6).
38. `nb.check_agree` — **from scratch:** trapezoid ∫₀^a u_z 2πR dR on 2001 points, `assert np.allclose(Q_trap,
    ch08.pipe_flow_rate(a, G)[0], rtol=1e-6)`; also the CV balance `assert np.isclose(-np.pi*a**2*G*1.0,
    -2*np.pi*a*1.0*ch08.pipe_wall_stress(a, G))` (pressure force on a 1 m slug = wall force, D07 step 5).
39. `nb.figure` — **Fig. 8.5 remake (N104 part)** (9 × 3.4 in): (a) the paraboloid u_z(R) across the tube (blue arrows at 9
    radii, the pipe walls grey) with the stress τ(R) as a second axis (rose, linear, τ₀ marked); (b) Q vs a on log–log
    (a from 0.1 mm to 1 cm at dp/dz = −1000 Pa/m), slope-4 guide, the points where Re crosses 2000 marked "turbulent
    beyond" (dashed continuation). Title "Flow rate grows like the fourth power of the radius". *see:* "a parabola, a
    straight stress line, and a straight log–log line of slope 4." *read:* "each decade of radius gives four decades of
    flow; past Re = 2000 (a ≈ 2.5 mm here) the laminar line is no longer trustworthy." *change:* "…the artery narrowed by
    20 %: 0.8⁴ = 0.41 — the same pressure pushes 59 % less blood, so the heart must raise the pressure 2.4× to keep Q."
40. `nb.plotly` — **IF2** `slider_figure(lambda a: {"u_z(R)": (R_over_a, ch08.pipe_poiseuille(R_over_a*a, a, -1000.0)),
    …}, "a", np.geomspace(2e-4, 2e-3, 13), unit="m", …)`, title carries Q and Re (the Q ∝ a⁴ surprise).
41. `nb.md` — **What would change if…** "…the walls moved instead of a pressure pushing — sideways, round the axis? Then
    the velocity is azimuthal, u_φ(R), and the same reduction gives circular Couette flow (C04)."
42. `nb.recap("R10", "A viscous flow with no net viscous force", "Outside a spinning cylinder the ideal-vortex velocity
    u_φ = Γ/2πR has shear stress $\\sigma_{R\\varphi}=\\mu R\\frac{d}{dR}(\\frac{u_\\varphi}R)=-\\frac{2\\mu\\Omega_1R_1^2}
    {R^2}$ but no net viscous force on any fluid element, because μ∇²u = −μ∇×ω = 0 when ω = 0 (Ch. 5).", where="Ch. 5
    §5.1")`
43. `nb.recap("R11", "Power in = dissipation", "The torque that keeps the cylinder spinning does work at the rate
    −2πR₁σ_Rφu_φ per unit length (positive: σ_Rφ < 0), and all of it is dissipated in the fluid: 4πμΩ₁²R₁² (Ch. 5's
    `dissipation_outside_cylinder`; Exercise 8.12). The book writes the power as (2πR₁)τ_Rφu_φ, which is negative with its
    own sign of σ_Rφ — the power delivered *to* the fluid carries a minus sign.", where="Ch. 5 §5.1")`
44. `nb.recap("R12", "Solid-body rotation", "A tank of fluid spun long enough rotates like a solid: $u_\\varphi=\\Omega R$
    (5.1), vorticity 2Ω everywhere, no shear at all.", where="Ch. 5 §5.1")`
45. `nb.core("C04", "Circular Couette flow $u_\\varphi=\\frac{1}{R_2^2-R_1^2}\\{[\\Omega_2R_2^2-\\Omega_1R_1^2]R-
    [\\Omega_2-\\Omega_1]\\frac{R_1^2R_2^2}{R}\\}$ (8.10) and its two limits", question="What swirl is set up between two
    rotating cylinders, and what happens when one of them is taken away?")`
46. `nb.md` — **The problem in plain words:** "Stir tea between a spoon handle and the cup wall; spin the inner cylinder of
    a lab rig (a Taylor–Couette cell, the classic apparatus of Ch. 11) or a viscometer. The fluid between two concentric
    rotating walls settles into a steady swirl. We want that swirl, the pressure that holds it in its circle, and what it
    becomes when the outer wall is at infinity or the inner one vanishes — the two vortices of Ch. 5."
47. `nb.md` — **The idea** (table): "Only two swirls satisfy the viscous balance, and the walls choose the mix:"
    | piece | u_φ | vorticity | shear | appears alone when |
    |---|---|---|---|---|
    | solid-body rotation | AR | 2A | none | no inner cylinder (R₁ → 0) |
    | free (line) vortex | B/R | 0 | yes, but no net force | outer wall at infinity (R₂ → ∞, Ω₂ = 0) |
48. `nb.note` — **N17 [B]** "**Two equations left.** With u = (0, u_φ(R), 0): the R-equation $-\\frac{u_\\varphi^2}{R}=-\\frac
    1\\rho\\frac{dp}{dR}$ says pressure rises outward to supply the centripetal acceleration; the φ-equation
    $0=\\mu\\frac{d}{dR}\\big[\\frac1R\\frac{d}{dR}(Ru_\\varphi)\\big]$ is a pure viscous balance (no pressure, no
    advection)." equations (R and φ momentum).
49. `nb.primer("Euler–Cauchy (equidimensional) ODE", "An ODE in which every term has the form R^k d^k u/dR^k (each
    derivative comes with the same power of R) is solved by trying u = R^λ: every term becomes a constant times R^λ, so
    λ must solve a polynomial. For R²u″ + Ru′ − u = 0 the polynomial is λ² − 1 = 0, so u = AR + B/R.", code="import sympy
    as sp\nR, lam = sp.symbols('R lambda')\nu = R**lam                          # trial power\nprint(sp.factor(sp.simplify((R**2*sp.diff(u, R, 2) + R*sp.diff(u, R) - u)/u)))   # (lambda - 1)*(lambda + 1)")`
    (**P187**)
50. `nb.derivation("D08", …)` — Part F D08 (10 steps), ref "8.10". Followed by `nb.note("**N18 [B]** The general solution
    $u_\\varphi(R)=AR+B/R$ (8.9) is solid-body rotation plus a line vortex. **N19 [B]** The walls give $A=\\frac{\\Omega_2
    R_2^2-\\Omega_1R_1^2}{R_2^2-R_1^2}$, $B=-\\frac{(\\Omega_2-\\Omega_1)R_1^2R_2^2}{R_2^2-R_1^2}$ (D08 steps 8–10).")`.
51. `nb.derivation("D09", …)` — Part F D09 (6 steps), ref "8.11".
52. `nb.note` — **N21 [B]** "**A viscous flow that is irrotational.** With R₂ → ∞ and Ω₂ = 0, $u_\\varphi=\\frac{\\Omega_1
    R_1^2}{R}$ (8.11) is the ideal vortex (5.2) with Γ = 2πΩ₁R₁² — shear stress everywhere, yet zero net viscous force
    (R10): the only completely irrotational viscous solution in the book. If gravity acts along the axis, the free surface
    dips toward the cylinder (Fig. 8.7). ⚠️ Ch. 5's `rotating_cylinder_flow(r, a, omega)` takes the cylinder's *vorticity*
    ω = 2Ω₁." **R12**'s solid-body rotation is the other limit: $u_\\varphi=\\Omega_2R$ (8.12).
53. `nb.worked_example("a 1 cm cylinder spinning at 1 rad/s inside a fixed 2 cm cup", "R₁ = 0.01 m, R₂ = 0.02 m, Ω₁ = 1
    rad/s, Ω₂ = 0. 1. $A=\\frac{0-1\\times10^{-4}}{4\\times10^{-4}-10^{-4}}=-\\frac13$ s⁻¹. 2. $B=-\\frac{(0-1)\\times10^{-4}
    \\times4\\times10^{-4}}{3\\times10^{-4}}=1.333\\times10^{-4}$ m²/s. 3. At mid-gap R = 0.015 m: $u_\\varphi=-\\frac13
    \\times0.015+\\frac{1.333\\times10^{-4}}{0.015}=-0.0050+0.0089=0.0039$ m/s. 4. Checks: at R₁, −0.00333 + 0.01333 =
    0.0100 = Ω₁R₁ ✓; at R₂, −0.00667 + 0.00667 = 0 ✓.")`
54. `nb.code` — `R1, R2 = 0.01, 0.02` · `u, (A, B) = ch08.circular_couette(0.015, R1, R2, 1.0, 0.0, return_coeffs=True)` ·
    `print(u, A, B)` · `print(ch08.circular_couette(0.05, R1, np.inf, 1.0, 0.0), ch05.rotating_cylinder_flow(0.05, R1,
    2*1.0))` · `print(ch08.circular_couette(0.01, 0.0, R2, 0.0, 3.0), VX.solid_body_rotation(0.01, 3.0))` ·
    `print(ch08.circular_couette_pressure(R2, R1, R2, 1.0, 0.0))` · `print(ch08.circular_couette_power(R1, R2, 1.0, 0.0,
    mu=MU_W))` · `print(ch08.parallel_flow_sympy("circular_couette")["profile"])`. *expect:* 0.003889 m/s, −0.3333 s⁻¹,
    1.333e-4 m²/s; 0.002 = 0.002 m/s (free vortex, parity with ch05 at ω = 2Ω₁); 0.03 = 0.03 m/s (solid body); pressure
    rise across the gap 0.0217 Pa; torque_inner = torque_outer = 1.676e-6 N m/m in magnitude, power_in = dissipation =
    1.676e-6 W/m (outer wall at rest, so power_out = 0); the sympy profile equals (8.10). *explain:* 1. the exact limit
    branches (R₂ = inf, R₁ = 0) agree with Ch. 5's vortex and solid-body functions; 2. the torque is the same at both
    walls — steady angular-momentum balance; 3. every watt the inner cylinder puts in is dissipated.
55. `nb.check_agree` — **from scratch:** `M = np.array([[R1, 1/R1], [R2, 1/R2]])`, `rhs = np.array([1.0*R1, 0.0*R2])`, `A_,
    B_ = np.linalg.solve(M, rhs)`; `assert np.allclose([A_, B_], [A, B], rtol=1e-12)`.
56. `nb.note` — **N20 [B]** "**The pressure that holds the swirl** (ours; the book says only that it 'can be determined').
    Integrating the R-balance $dp/dR=\\rho u_\\varphi^2/R$ with u_φ = AR + B/R: $p(R)=p_1+\\rho\\Big[\\frac{A^2}{2}(R^2-R_1^2)
    +2AB\\ln\\frac{R}{R_1}-\\frac{B^2}{2}\\Big(\\frac{1}{R^2}-\\frac{1}{R_1^2}\\Big)\\Big]$. In the atmosphere this
    centripetal balance is the *cyclostrophic* balance of tornadoes and dust devils (Ch. 13)." equation (p(R)) +
    `nb.code`: `Rg = np.linspace(R1, R2, 2001); pg = ch08.circular_couette_pressure(Rg, R1, R2, 1.0, 0.0); ug =
    ch08.circular_couette(Rg, R1, R2, 1.0, 0.0); assert np.allclose(np.gradient(pg, Rg)[5:-5], (1000*ug**2/Rg)[5:-5],
    rtol=1e-5)` — finite differences of p match ρu²/R.
57. `nb.note` — **N22 [B]** "**Why these three flows are solvable.** All three are confined between walls, and symmetry
    kills the advective acceleration: u·∇u = 0 in the channel and the pipe, while in circular Couette flow only the
    centripetal part −(u_φ²/R)e_R survives, balanced by the pressure. Other exact solutions exist (§8.4 and the exercises);
    Ch. 10 uses such exact solutions to test flow solvers." + `nb.code`: `for case in ("channel", "pipe",
    "circular_couette"): print(case, ch08.advective_acceleration_check(case))` *expect:* [0, 0, 0], [0, 0, 0], [−u_φ²/R,
    0, 0].
58. `nb.figure` — **Figs. 8.6–8.7 remakes (N104 [B])** (10 × 3.4 in, three panels): (a) u_φ(R) for R₁/R₂ = 0.5 and four
    (Ω₁, Ω₂) pairs — inner only (1, 0), outer only (0, 1), co-rotating (1, 1: solid body), counter-rotating (1, −1) — with
    the AR (orange dashed) and B/R (rose dashed) parts of the inner-only case; (b) the annulus seen from above (our
    remake of Fig. 8.6) with arrows ∝ u_φ at 12 angles × 5 radii; (c) R₂ → ∞: u_φ(R) = Ω₁R₁²/R (blue) with the
    free-surface dip $z_s=-\\Gamma^2/(8\\pi^2gR^2)$ (ours, from Ch. 5's funnel) drawn below as a profile (our Fig. 8.7).
    Title "Every swirl between cylinders is AR + B/R". *see:* "four curves joining the two wall speeds; the co-rotating one
    is a straight line; the counter-rotating one crosses zero inside the gap." *read:* "a straight line is pure AR (B = 0,
    rigid rotation); a curve bending like 1/R carries a vortex part; where u_φ changes sign the fluid turns the other way."
    *change:* "…the outer cylinder turned faster than the inner at the same sense: (Ru_φ)² grows outward and the flow is
    stable (Rayleigh, Ch. 11); spin only the inner one and it is unstable above a critical speed — the Taylor vortices."
59. `nb.plotly` — **IF3** `slider_figure(lambda r: {"u_φ(R)": …, "solid body Ω₂R": …, "free vortex Ω₁R₁²/R": …}, "Ω₂/Ω₁",
    np.linspace(-2, 2, 21), …)` at R₁/R₂ = 0.5; title carries A, B and the Rayleigh verdict from
    `ch08.circular_couette_state`.
60. `nb.md` — **What would change if…** "…the walls were not exactly parallel — a thin gap that narrows slowly, as under a
    bearing pad? The flow is then *nearly* parallel, the nonlinear term is small instead of zero, and we need a way to
    decide what may be dropped: the lubrication approximation (C05)." Hint: "(The backup explainer
    `rotating_cylinders_couette` animates this block if one of the nine explainers fails review.)"

### A.3 §8.3 Elementary Lubrication Theory — R13, C05 (+N23, N24, N26–N31 · D10 D11 · P188 · E2), C06 (+N25, N32–N34, N36 · D12 D13 · P189 P190), C07 (+N35, N105 · D14 · IF4 · live · E3), C08 (+N37 · D15 · P191 P192 · A3 · E4)
1. `nb.section("8.3", "Elementary Lubrication Theory", intro="**What is this section about?** Real walls are rarely
   exactly parallel, but in a thin gap they are *nearly* parallel. We learn how to decide, by scaling with two different
   lengths, which terms of the Navier–Stokes equations can be dropped; the result is a local Couette-plus-Poiseuille
   profile and one equation for the pressure (the Reynolds equation). It explains why an oil film under a bearing pad
   carries tonnes, why viscous flow between glass plates draws the streamlines of ideal flow, and how a honey drop
   spreads under its own weight. The same scaling argument becomes the boundary-layer approximation in Ch. 9.")`
2. `nb.recap("R13", "The field equations in the gap", "In two dimensions: continuity $\\frac{\\partial u}{\\partial x}+
   \\frac{\\partial v}{\\partial y}=0$ (6.2) and the x-component of (8.1), $\\frac{\\partial u}{\\partial t}+u\\frac{\\partial u}
   {\\partial x}+v\\frac{\\partial u}{\\partial y}=-\\frac1\\rho\\frac{\\partial p}{\\partial x}+\\frac\\mu\\rho\\big(\\frac{
   \\partial^2u}{\\partial x^2}+\\frac{\\partial^2u}{\\partial y^2}\\big)$ (8.13a) — nothing new, just written out.",
   where="Ch. 4 §4.6")`
3. `nb.core("C05", "The lubrication approximation $0\\cong-\\frac1\\rho\\frac{\\partial p}{\\partial x}+\\nu\\frac{\\partial^2u}
   {\\partial y^2}$ (8.17a)", question="Why can we throw away inertia in a thin film even when UL/ν is in the thousands?")`
4. `nb.md` — **The problem in plain words:** "Under a sliding bearing pad the oil film is 50 µm thick and 5 cm long —
   a thousand times longer than thick. The oil moves at metres per second and UL/ν ≈ 2500, which in a pipe would be
   almost turbulent. Yet lubrication engineers treat the film as purely viscous. We want the argument that makes this
   legitimate — and the number that must be small instead of 1/Re."
5. `nb.note` — **N23 [B]** "**The lubrication idea.** In a narrow passage of gap h and length L ≫ h (Fig. 8.8) the flow is
   nearly parallel; pressure and viscous forces dominate, and a gentle curvature of the passage does not matter as long as
   its radius is much larger than h. Exactly the same reasoning, with the gap replaced by a thin layer next to a wall,
   gives the boundary-layer approximation (Ch. 9 §9.1)."
6. `nb.md` — **The idea** (ASCII): "two length scales, and continuity fixes the size of v:"
   ```
   ↑ y ~ h = εL (small)          ∂/∂x ~ 1/L   (slow along the gap)
   │  ═══════════════════        ∂/∂y ~ 1/h   (fast across it)
   │  ──► ──► ──► u ~ U          continuity: U/L ~ v/h  ⇒  v ~ εU
   └──────────────────────► x ~ L
   ```
7. `nb.note` — **N24 [B]** "**The y-momentum equation** $\\frac{\\partial v}{\\partial t}+u\\frac{\\partial v}{\\partial x}+
   v\\frac{\\partial v}{\\partial y}=-\\frac1\\rho\\frac{\\partial p}{\\partial y}+\\frac\\mu\\rho\\big(\\frac{\\partial^2v}{
   \\partial x^2}+\\frac{\\partial^2v}{\\partial y^2}\\big)$ (8.13b). ⚠️ The book prints −(1/ρ)∂p/∂**x** here; the
   y-equation needs ∂p/∂**y** — the scaled (8.16b) confirms it, and `ch08.lubrication_nondim_sympy()` shows the printed
   form gives a coefficient set that does not match (8.16b)." equation (8.13b).
8. `nb.primer("anisotropic scaling with two length scales", "When a flow is long and thin, measure x in units of the length
   L and y in units of the thickness h = εL; derivatives then carry different sizes, ∂/∂x = (1/L)∂/∂x* but ∂/∂y =
   (1/εL)∂/∂y*. The velocity across the layer gets its own scale from continuity: U/L must balance v/h, so v ~ εU.
   Written this way every starred quantity is of order one and the size of each term sits in its coefficient.",
   code="L, eps, U = 0.05, 1e-3, 5.0            # length [m], fineness h/L, speed [m/s]\nh = eps*L\nprint(U/L, (eps*U)/h)                 # both 100 1/s: the two continuity terms balance")` (**P188**)
9. `nb.md` — gloss **bookkeeping of a small parameter** (P68 reminder): "Multiply the whole equation by the factor that makes
   the term you expect to dominate have coefficient 1; then every other coefficient is a power of ε times a group, and
   you read off which terms are small."
10. `nb.note` — **N26 [B]** "**The scalings (8.14):** $x^*=\\frac xL,\\ y^*=\\frac yh=\\frac{y}{\\varepsilon L},\\ t^*=\\frac{Ut}{L},
    \\ u^*=\\frac uU,\\ v^*=\\frac{v}{\\varepsilon U},\\ p^*=\\frac{p}{P_a}$ with the fineness ratio ε = h/L and the atmospheric
    pressure P_a." equation (8.14).
11. `nb.derivation("D10", …)` — Part F D10 (14 steps, ★★★, `check_src` below), ref "8.16". **N27 [B]** (scaled continuity
    (8.15) with no coefficient — steps 3–4), **N28 [B]** ((8.16a) — steps 5–10), **N29 [B]** ((8.16b) — steps 11–14) are
    tagged in one note right after it: `nb.note("**N27 [B]** $\\frac{\\partial u^*}{\\partial x^*}+\\frac{\\partial v^*}{
    \\partial y^*}=0$ (8.15) has no coefficient: mass is conserved exactly — v is small *because* u varies slowly.
    **N28 [B]** $\\varepsilon^2\\mathrm{Re}_L\\big(\\frac{\\partial u^*}{\\partial t^*}+u^*\\frac{\\partial u^*}{\\partial x^*}+
    v^*\\frac{\\partial u^*}{\\partial y^*}\\big)=-\\frac1\\Lambda\\frac{\\partial p^*}{\\partial x^*}+\\varepsilon^2\\frac{
    \\partial^2u^*}{\\partial x^{*2}}+\\frac{\\partial^2u^*}{\\partial y^{*2}}$ (8.16a), with Re_L = ρUL/μ and the bearing
    number Λ = μUL/(P_a h²). **N29 [B]** $\\varepsilon^4\\mathrm{Re}_L(\\ldots)=-\\frac1\\Lambda\\frac{\\partial p^*}{\\partial
    y^*}+\\varepsilon^4\\frac{\\partial^2v^*}{\\partial x^{*2}}+\\varepsilon^2\\frac{\\partial^2v^*}{\\partial y^{*2}}$ (8.16b):
    only the pressure term lacks a power of ε — the seed of 'pressure is constant across a thin layer' (Ch. 9, Ch. 13).")`.
    **check_src (D10):** `sol = ch08.lubrication_nondim_sympy()  # repeats steps 1–14 symbolically` · `eps, ReL, Lam =
    sp.symbols("epsilon Re_L Lambda", positive=True)` · `print(sol["continuity_star"])  # (8.15): no coefficient` ·
    `assert sp.simplify(sol["x_coeffs"]["inertia"] - eps**2*ReL) == 0  # (8.16a) inertia coefficient` · `assert
    sp.simplify(sol["x_coeffs"]["pressure"] - 1/Lam) == 0` · `assert sp.simplify(sol["y_coeffs"]["inertia"] - eps**4*ReL)
    == 0  # (8.16b)` · `assert sp.simplify(sol["y_coeffs"]["diff_across"] - eps**2) == 0` · `print(sol["printed_13b_y_coeffs"])
    # the printed ∂p/∂x gives a set that does not match (8.16b)` — every line commented.
12. `nb.derivation("D11", …)` — Part F D11 (7 steps), ref "8.17".
13. `nb.note` — **N30 [B]** "**Pressure is uniform across the gap:** $0\\cong-\\frac1\\rho\\frac{\\partial p}{\\partial y}$
    (8.17b), with an error of relative size ε². The same statement across a thin atmospheric or oceanic layer is the
    hydrostatic approximation of Ch. 13." equation (8.17b). Markdown: "> ⚠️ **Common confusion:** the book prints (8.17a) as
    $0\\cong-\\frac1\\rho\\frac{\\partial p}{\\partial x}+\\frac{\\partial^2u}{\\partial y^2}$ — without the ν. The units do not
    balance (m/s² against 1/(m s)); re-dimensionalising (D11 step 5) restores $\\nu\\frac{\\partial^2u}{\\partial y^2}$. Its
    own next line, (8.18), carries the 1/μ correctly."
14. `nb.worked_example("an engine-bearing oil film", "h = 50 µm, L = 5 cm, U = 5 m/s, a light oil with μ = 0.05 Pa s
    and ρ = 870 kg/m³, so ν = μ/ρ = 5.75 × 10⁻⁵ m²/s (numbers ours). 1. $\\varepsilon=h/L=5\\times10^{-5}/0.05=10^{-3}$.
    2. $\\mathrm{Re}_L=UL/\\nu=5\\times0.05/5.75\\times10^{-5}=4350$ — in a pipe this would be turbulent. 3.
    $\\varepsilon^2\\mathrm{Re}_L=10^{-6}\\times4350=4.35\\times10^{-3}$ — small: inertia is under half a percent of friction.
    4. The natural pressure scale $\\mu UL/h^2=0.05\\times5\\times0.05/(2.5\\times10^{-9})=5\\times10^6$ Pa ≈ 49 atm. 5. With
    the book's atmospheric scale, Λ = μUL/(P_a h²) = 5 × 10⁶/101 325 = 49 — not 'near unity'.")`
15. `nb.code` — `s = ch08.lubrication_scales(0.05, 50e-6, 5.0, 870.0, 0.05)` · `print(s)` ·
    `print(ch08.lubrication_term_magnitudes(0.05, 50e-6, 5.0, 870.0, 0.05, p_scale="viscous"))` · `print(ch08.
    lubrication_term_magnitudes(0.05, 50e-6, 5.0, 870.0, 0.05, p_scale="atm"))`. *expect:* eps 1e-3, Re_L 4350,
    eps2_Re_L 4.35e-3, Lambda 49.3, p_visc 5.0e6 Pa; viscous scale: x_inertia 4.35e-3, x_pressure 1, x_diff_along 1e-6,
    x_diff_across 1, y_inertia 4.35e-9, y_pressure 1, y_diff_along 1e-12, y_diff_across 1e-6; atmospheric scale: both
    pressure coefficients 1/49.3 = 0.0203. *explain:* the coefficients are exactly those of (8.16a) and (8.16b); with the
    viscous pressure scale μUL/h² every kept term is O(1).
16. `nb.check_agree` — **from scratch:** `eps = 50e-6/0.05; ReL = 870*5.0*0.05/0.05; Lam = 0.05*5.0*0.05/(101325*(50e-6)**2)`
    · `mine = [eps**2*ReL, 1/Lam, eps**2, 1.0, eps**4*ReL, 1/Lam, eps**4, eps**2]` · `lib = ch08.lubrication_term_magnitudes(
    0.05, 50e-6, 5.0, 870.0, 0.05, p_scale="atm")` · `assert np.allclose(mine, list(lib.values()), rtol=1e-12)` (keys in
    the order of Part C 2.2).
17. `nb.figure` — (8 × 3.4 in, two panels) term magnitudes of (8.16a) (left) and (8.16b) (right) against ε on log–log
    axes (ε from 10⁻⁴ to 0.3) at fixed Re_L = 4350 with the viscous pressure scale: inertia (teal, slopes 2 and 4),
    pressure (orange, flat at 1), streamwise diffusion (rose light, slopes 2 and 4), cross-gap diffusion (rose, flat at 1
    and slope 2); a vertical dashed line where ε²Re_L = 1 (ε = 0.015) labelled "lubrication fails to the right". Title "In
    a thin gap inertia is weighed with ε²Re_L, not Re_L". *see:* "two flat lines (pressure and cross-gap friction) and
    falling lines for everything else." *read:* "read each term's size at your ε: at ε = 10⁻³ inertia is 4.35 × 10⁻³ of
    the kept terms; only near ε ≈ 0.015 does it catch up." *change:* "…Re_L were 100 times larger: the teal line moves up
    two decades and lubrication fails already at ε ≈ 0.0015."
18. `nb.note` — **N31 [B]** "**How small is small?** For our engine film ε²Re_L = 4.35 × 10⁻³, so the lubrication balance
    is excellent. (The book's own example uses a different oil film; its number lives in our private test file.) ⚠️ With
    p* = p/P_a the bearing number Λ = μUL/(P_a h²) is 49 here and ≈ 10³ for thinner, longer films — not 'near unity' as
    the text assumes. Nothing in the ordering changes (Λ only rescales p*), but the physical pressure scale is μUL/h²,
    tens to hundreds of atmospheres."
19. `nb.explainer("lubrication_scaling", heading="Why can a thin film ignore inertia?", why="Drag the gap, the length and
    the speed and watch every term of the scaled momentum equations (8.16a) and (8.16b) as a bar on a log scale: inertia
    falls like ε²Re_L, the pressure and the cross-gap friction stay, and the y-equation leaves only the pressure — an
    ordering argument you can see move.", tries=["Start from the engine-film preset and double U three times — the teal
    inertia bar climbs but stays far below the pressure bar.", "Press 'thick gap': ε = 0.2 and Re_L = 10⁴ put inertia on
    top and the status turns amber.", "Switch the pressure scale to P_a (the book's) and read Λ — then back to μUL/h².",
    "Open the Derivation tab and step through D10: each substitution lights the bar it creates."])`
20. `nb.md` — **What would change if…** "…we now use (8.17a) to find the velocity everywhere in a gap of any slowly varying
    shape? Integrating twice across the gap gives a local Couette-plus-Poiseuille profile (C06)."
21. `nb.core("C06", "The lubrication profile $u\\cong-\\frac{h^2}{2\\mu}\\frac{\\partial p}{\\partial x}\\frac yh(1-\\frac yh)+
    U_h\\frac yh+U_0$ (8.19), the gap flux and the Reynolds equation", question="Given the shape of the gap, what is the
    velocity at every point — and what single equation decides the pressure?")`
22. `nb.md` — **The problem in plain words:** "Under a tilted pad, in a knee joint, in the gap of an injection mould, the
    gap height changes slowly along the flow. At each station the flow 'thinks' it is between parallel plates — so its
    profile is the Couette–Poiseuille profile of C02 with the local gap h(x) and the local pressure gradient. But the
    pressure gradient is unknown; mass conservation across the whole gap has to decide it."
23. `nb.md` — **The idea** (ASCII): "the channel formula, applied station by station:"
    ```
    station 1 (wide)       station 2            station 3 (narrow)
    ══════════ U_h         ══════════ U_h       ══════════ U_h
     ─►                     ──►                   ───►
     ►  + Poiseuille        ─► + smaller          ──►      same flux q must pass every station
    ═══════════ U₀        ═══════════ U₀         ═══════════ U₀
    ```
24. `nb.note` — **N25 [B]** "**Gap boundary conditions.** u = U₀(t) on the lower wall y = 0 and u = U_h(t) on the upper
    wall y = h(x, t); the pressure may depend on t. Time only enters as a parameter: the lubrication balance (8.17) has no
    ∂/∂t term, yet h, U₀, U_h and p may all change slowly."
25. `nb.derivation("D12", …)` — Part F D12 (8 steps), ref "8.19". Then `nb.note("**N33 [B]** The intermediate form
    $u(x,y,t)\\cong\\frac1\\mu\\frac{\\partial p(x,t)}{\\partial x}\\frac{y^2}{2}+Ay+B$ (8.18) has 'constants' A, B that may
    depend on x and t (step 3). ⚠️ The text says (8.16a) 'can be integrated twice' — it is the simplified (8.17a) that is
    integrated.", equation=r"u(x,y,t)\cong\frac{1}{\mu}\frac{\partial p(x,t)}{\partial x}\frac{y^{2}}{2}+Ay+B",
    ref="8.18")`.
26. `nb.note` — **N34 [B]** "**The U₀ term.** The printed (8.19) adds U₀ to the Couette part, so at the upper wall it gives
    u(h) = U_h + U₀ instead of U_h. The consistent profile is $U_h\\frac yh+U_0\\big(1-\\frac yh\\big)$ (D12 step 7). Every
    example of the chapter has U₀ = 0, so no result changes; our code uses the consistent form and keeps the printed one
    as `form='book'` for a test that must fail." + `nb.code`: `print(ch08.lubrication_velocity(50e-6, 50e-6, 0.0,
    U_h=2.0, U_0=1.0), ch08.lubrication_velocity(50e-6, 50e-6, 0.0, U_h=2.0, U_0=1.0, form="book"))` *expect:* 2.0 (the
    upper-wall speed ✓) and 3.0 (the printed form fails).
27. `nb.primer("Leibniz rule with a moving upper limit", "Differentiating an integral whose upper limit moves adds a
    boundary term: ∂/∂x ∫₀^{h(x)} u dy = ∫₀^{h} ∂u/∂x dy + u(x, h)·∂h/∂x. The second term counts what enters or leaves
    because the limit itself moved (Ch. 3 met the fixed-limit version, P109).", code="import sympy as sp\nx, y = sp.symbols('x y')\nh = 1 + x**2                 # a moving upper limit\nu = x*y                       # a test integrand\nlhs = sp.diff(sp.integrate(u, (y, 0, h)), x)\nrhs = sp.integrate(sp.diff(u, x), (y, 0, h)) + u.subs(y, h)*sp.diff(h, x)\nprint(sp.simplify(lhs - rhs))          # 0")` (**P189**)
28. `nb.md` — gloss **kinematic condition at a moving wall** (P132 reminder): "A fluid particle on the upper wall y = h(x, t)
    stays on it, so D(y − h)/Dt = 0 there: v = ∂h/∂t + u ∂h/∂x at y = h."
29. `nb.derivation("D13", …)` — Part F D13 (12 steps, ★★★), ref "". **check_src:** `r = ch08.reynolds_equation_sympy()
    # D13 built symbolically: profile → flux → integrated continuity` · `print(r["q"])  # −h³p_x/(12μ) + (U0 + Uh)h/2` ·
    `assert sp.simplify(r["leibniz_residual"]) == 0  # step 2: the Leibniz rule holds for the profile` · `assert
    sp.simplify(r["kinematic_cancellation"]) == 0  # steps 5–7: the U_h h_x terms cancel` · `assert
    sp.simplify(r["reynolds_residual"]) == 0  # step 12: h_t + q_x = 0 equals ∫(u_x + v_y)dy = 0` — every line commented.
30. `nb.note` — **N32 [B]** "**Completing a lubrication problem.** The profile alone is not a solution: we still need p(x).
    Integrating continuity across the gap gives $\\frac{\\partial h}{\\partial t}+\\frac{\\partial q}{\\partial x}=0$ with the
    gap flux $q=\\int_0^hu\\,dy=-\\frac{h^3}{12\\mu}\\frac{\\partial p}{\\partial x}+\\frac{(U_0+U_h)h}{2}$ — the 1-D **Reynolds
    equation** (the book leaves it to Exercises 8.19–8.20). For a steady gap q is the same at every station, which is one
    ODE for p(x); two end pressures close it (C07)." equations (q, Reynolds).
31. `nb.primer("scipy.integrate.cumulative_trapezoid", "cumulative_trapezoid(f, x, initial=0) returns the running integral
    ∫_{x₀}^{x_i} f dx at every sample — the numerical antiderivative. We use it to turn a pressure gradient into a
    pressure profile.", code="from scipy.integrate import cumulative_trapezoid\nimport numpy as np\nx = np.linspace(0, 1, 5)\nprint(cumulative_trapezoid(2*x, x, initial=0))   # ≈ x**2: [0, 0.0625, 0.25, 0.5625, 1.0]")` (**P190**)
32. `nb.worked_example("the flux through a 50 µm gap", "h = 50 µm, upper wall U_h = 5 m/s, lower wall fixed, μ = 0.05 Pa s.
    1. With no pressure gradient: $q=U_hh/2=5\\times5\\times10^{-5}/2=1.25\\times10^{-4}$ m²/s. 2. The gradient that halves
    it: $\\frac{h^3}{12\\mu}\\frac{dp}{dx}=\\frac{U_hh}{4}$ ⇒ $\\frac{dp}{dx}=\\frac{3\\mu U_h}{h^2}=\\frac{3\\times0.05\\times5}
    {2.5\\times10^{-9}}=3\\times10^8$ Pa/m — 3 bar per millimetre. 3. At dp/dx = 6μU_h/h² = 6 × 10⁸ Pa/m the flux is zero:
    the pressure-driven backflow exactly cancels the drag.")`
33. `nb.code` — `h, Uh, mu = 50e-6, 5.0, 0.05` · `print(ch08.lubrication_flux(h, 0.0, U_h=Uh, mu=mu),
    ch08.lubrication_flux(h, 3*mu*Uh/h**2, U_h=Uh, mu=mu))` · `yy = np.linspace(0, h, 5)` ·
    `print(ch08.lubrication_velocity(yy, h, 3*mu*Uh/h**2, U_h=Uh, mu=mu))` · `x = np.linspace(0, 0.05, 401); hx = h*(1 + 0.5*x/0.05)` ·
    `p, q = ch08.reynolds_pressure_1d(x, hx, U_0=-Uh, mu=mu)` (pad frame: the floor slides backwards at −U under the pad, which is at rest) · `print(p.max(), q)`.
    *expect:* 1.25e-4 and 6.25e-5 m²/s; the velocity profile at the halving gradient [0, −0.156, 0.625, 2.344, 5.0] m/s (a thin reversed layer at the fixed wall: 3μU/h² exceeds the backflow threshold 2μU/h² of C02); p.max = 1.00 MPa at x/L = 0.40 (the α = 0.5 pad of C07 — cross-checked there), q = −1.5e-4 m²/s (= C₁ = −(1 + α)Uh₀/(2 + α), the pad-frame flux).
    *explain:* 1. `lubrication_flux` is N32's q; 2. `reynolds_pressure_1d` finds the constant flux that makes the end
    pressures equal and integrates dp/dx = 12μ(q̄ − q)/h³ with `cumulative_trapezoid`.
34. `nb.check_agree` — **from scratch:** `q_num = integrate.quad(lambda yy: ch08.lubrication_velocity(yy, h, 1e8, U_h=Uh,
    mu=mu), 0, h)[0]`; `assert np.allclose(q_num, ch08.lubrication_flux(h, 1e8, U_h=Uh, mu=mu), rtol=1e-10)`.
35. `nb.figure` — **Fig. 8.8 remake with the slider gap** (9 × 3.2 in): the gap h(x) = h₀(1 + αx/L) for α = 1.5, drawn
    with the height exaggerated 200× (stated in the title), velocity profiles (blue arrows) at 7 stations from
    `ch08.slider_gap_velocity(x, y, 50e-6, 1.5, 0.05, 5.0, mu=0.05)` in the ground frame; the Couette part (rose dashed)
    and the Poiseuille part (orange dashed) at the widest and narrowest stations. Title "Couette plus Poiseuille, station by
    station". *see:* "profiles that are nearly straight in the middle of the pad, bulging forward at the narrow end and
    backward near the stationary floor at the wide end." *read:* "the area under every profile (the flux in the pad frame)
    is the same; where the gap is wide the pressure gradient must push fluid back to keep it so." *change:* "…the pad slid
    the other way (U < 0): every profile flips, the pressure hump becomes a suction dip (C07)."
36. `nb.note` — **N36 [B]** "**Hele-Shaw flow (Example 8.2): viscous flow that draws ideal streamlines.** Between two plates
    z = 0 and z = h the lubrication balances $0\\cong-\\frac1\\rho\\frac{\\partial p}{\\partial x}+\\nu\\frac{\\partial^2u}{
    \\partial z^2}$, $0\\cong-\\frac1\\rho\\frac{\\partial p}{\\partial y}+\\nu\\frac{\\partial^2v}{\\partial z^2}$,
    $0\\cong-\\frac1\\rho\\frac{\\partial p}{\\partial z}$ give $u\\cong-\\frac{1}{2\\mu}\\frac{\\partial p}{\\partial x}z(h-z)=
    \\frac{\\partial\\phi}{\\partial x}$ with $\\phi=-\\frac{z(h-z)}{2\\mu}p$. Integrating continuity over the gap (the 2-D
    Reynolds equation with constant h, depth-averaged velocity $\\bar{\\mathbf u}=-\\frac{h^2}{12\\mu}\\nabla p$) leaves
    $\\frac{\\partial^2p}{\\partial x^2}+\\frac{\\partial^2p}{\\partial y^2}=0$ — so φ obeys $u=\\partial\\phi/\\partial x$,
    $v=\\partial\\phi/\\partial y$, $\\nabla^2\\phi=0$ (6.10, 6.12), the equations of 2-D ideal flow (Ch. 6). Dye injected
    between glass plates therefore traces ideal-flow streamlines, except in layers of thickness ~h at obstacles, where
    no slip holds (Exercise 8.34). ⚠️ The book prints these lubrication equations without the ν, and writes the wall
    conditions at 'y = 0, h' — the gap coordinate here is z." equations (u, φ, Laplace) + `nb.figure`: (6 × 4 in) Hele-Shaw
    cell with a disc of radius a = 1 cm, depth-averaged dye lines from `ch08.hele_shaw_cylinder(X, Y, None, 1e-3, 0.01,
    1e-3)` (blue) overlaid on the ideal-flow cylinder streamlines of `PF.Flow` (uniform + doublet, muted dashed), and an
    inset of the parabola u(z) across the gap. *see:* "the blue dye lines lie on the dashed ideal streamlines." *read:*
    "the pressure is harmonic, so the depth-averaged velocity is a potential flow; only a thin layer at the disc (too thin
    to see here) differs." *change:* "…the gap were as wide as the disc: the ε ≪ 1 assumption fails, inertia and the
    no-slip layer spread and the pattern loses its fore–aft symmetry." V3 check cell: `LS.solve_poisson` on a masked
    128² grid (64² FAST) for p around the disc, compared with the ideal-flow pressure potential (max error printed, falls
    by ≈ 4 when the grid doubles).
37. `nb.md` — **What would change if…** "…the gap were the space under a real, tilted pad with the same pressure p_e at both
    ends? The flux must be constant along the pad, which forces a pressure hump inside — and the hump carries a load
    (C07)."
38. `nb.core("C07", "The slider bearing (Example 8.1): $p-p_e=\\frac{6\\mu LU}{h_o^2}\\frac{\\alpha(x/L)(1-x/L)}{(2+\\alpha)
    (1+\\alpha x/L)^2}$ and $W=\\frac{\\alpha\\mu L^2U}{2h_o^2}$", question="How does a thin viscous film carry a load — and
    why only in one direction?")`
39. `nb.md` — **The problem in plain words:** "A tilted pad skates on a film of oil over a flat surface, as in the thrust
    bearing of a ship's propeller shaft or a hydroelectric turbine carrying hundreds of tonnes. The pad never touches the
    surface. Where does the upward force come from, and why does the bearing fail — the pad is sucked down — if it slides
    the wrong way?"
40. `nb.md` — **The idea** (ASCII): "same flux through a narrowing gap ⇒ pressure must rise inside:"
    ```
            W ↓                       pad moves → at U   (pad frame: floor moves ← at U)
    p_e  ┌───────────────────┐  p_e
         │ h₀(1+α)  →  h₀     │     wide end ─── drag brings in more than the narrow end lets out
    ═════════════════════════════   ⇒ p rises until Poiseuille backflow evens the flux
    p(x): ____/‾‾‾\____             hump between the ends, ∫(p − p_e)dx = W
    ```
41. `nb.md` — gloss **antiderivatives and a small-α expansion** (P106, P26 reminders): "$\\int(1+\\alpha x/L)^{-n}dx=\\frac{L}
    {\\alpha}\\frac{(1+\\alpha x/L)^{1-n}}{1-n}$ (substitute s = 1 + αx/L); for small α, (1 + αx/L)^{-2} ≈ 1 − 2αx/L and
    2 + α ≈ 2."
42. `nb.derivation("D14", …)` — Part F D14 (15 steps, ★★★), ref "". **check_src:** `sb = ch08.slider_bearing_sympy()  # D14 in
    sympy` · `print(sb["C1"], sb["C2"])  # step 10–11: C1 = −(1+α)U h0/(2+α)` · `assert sp.simplify(sb["ode_residual_exact"])
    == 0  # the squared-denominator p satisfies dp/dx = −12μC1/h³ − 6μU/h²` · `assert sp.simplify(sb["ode_residual_book"])
    != 0  # the printed first-power form does not` · `assert all(sp.simplify(r) == 0 for r in sb["bc_residuals"])  # p(0)
    = p(L) = p_e` · `print(sb["W_linear"], sp.series(sb["W_exact"], sp.Symbol("alpha"), 0, 3))  # exact W → linear W as
    α → 0` — every line commented.
43. `nb.md` — "> ⚠️ **Common confusion (two printing slips in Example 8.1):** the intermediate integrals are printed with
    (1 − αx/L), although the gap is h₀(1 + αx/L); and the final exact pressure is printed with (1 + αx/L) to the first
    power in the denominator. Substituting back into dp/dx = −12μC₁/h³ − 6μU/h² shows the denominator must be
    **squared** (D14 step 13, checked above). The printed C₁, C₂, the linear-α pressure and the load W are right."
44. `nb.worked_example("a 5 cm pad on a 50 µm film", "L = 0.05 m, h₀ = 50 µm, U = 5 m/s, μ = 0.05 Pa s, α = 0.1. 1. Scale:
    $6\\mu LU/h_o^2=6\\times0.05\\times0.05\\times5/(2.5\\times10^{-9})=3\\times10^7$ Pa. 2. At the middle x = L/2:
    $\\frac{\\alpha(x/L)(1-x/L)}{(2+\\alpha)(1+\\alpha/2)^2}=\\frac{0.1\\times0.25}{2.1\\times1.1025}=0.0108$ ⇒ p − p_e =
    3.24 × 10⁵ Pa ≈ 3.2 atm. 3. Linear formula: $\\frac{3\\alpha\\mu LU}{h_o^2}\\times0.25=3.75\\times10^5$ Pa. 4. Load
    $W=\\alpha\\mu L^2U/(2h_o^2)=0.1\\times0.05\\times0.0025\\times5/(5\\times10^{-9})=1.25\\times10^4$ N/m — 12.5 kN per
    metre of pad width, about 1.3 tonnes. 5. Reverse the motion (U = −5 m/s): W = −12.5 kN/m, the pad is pulled down.")`
45. `nb.code` — `x = np.linspace(0, 0.05, 401)` · `for model in ("exact", "linear", "book"): print(model,
    ch08.slider_bearing(0.025, 50e-6, 0.1, 0.05, 5.0, mu=0.05, model=model))` · `st = ch08.slider_bearing_state(50e-6, 0.1,
    0.05, 5.0, mu=0.05)` · `print({k: st[k] for k in ("C1", "p_max", "x_pmax", "W_exact", "W_linear", "err_linear",
    "p_max_atm", "inlet_backflow")})` · `print(ch08.slider_optimum_taper())` · `print(ch08.slider_bearing_load(50e-6, 0.1,
    0.05, -5.0, mu=0.05, model="linear"))`. *expect:* exact 3.239e5 Pa, linear 3.75e5 Pa, as printed (first power) 3.401e5 Pa;
    C1 = −1.310e-4 m²/s, p_max = 3.247e5 Pa at x/L = 0.476, W_exact = 1.081e4, W_linear = 1.25e4 N/m, err ≈ 15.6 %,
    p_max_atm 3.20, inlet_backflow False; optimum α = 1.189 (1 + α = 2.189), W* = 0.02671; reversed U: −1.25e4 N/m.
    *explain:* 1. the three models: corrected exact, O(α), and the printed slip (a ghost); 2. the state bundles the pad-frame
    flux C₁, the peak, the two loads and their difference; 3. the optimum taper maximises the exact load.
46. `nb.check_agree` — **from scratch (curation §7):** for α = 0.1: define `dpdx(x, q) = 12*mu*(q - U*h(x)/2)/h(x)**3` in
    the pad frame… **use the ground-frame form of D14 step 5**: `dpdx = lambda x, C1: -12*mu*C1/h(x)**3 - 6*mu*U/h(x)**2`;
    find C₁ with `optimize.brentq(lambda C1: integrate.quad(lambda x: dpdx(x, C1), 0, L)[0], -U*h0, 0.0)` (p(L) = p(0));
    then `p_mine = cumulative_trapezoid(dpdx(x, C1), x, initial=0)`; `assert np.allclose(C1, st["C1"], rtol=1e-8)`;
    `assert np.allclose(p_mine, ch08.slider_bearing(x, 50e-6, 0.1, 0.05, 5.0, mu=0.05), rtol=1e-6, atol=1.0)`.
47. `nb.figure` — **Fig. 8.9 remake (N105 [B], with the pad of Fig. 8.9; Figs. 8.10–8.11 are remade in C06 and C08)**
    (9 × 3.4 in, two panels): (a) p(x) − p_e for α = 0.1, 0.5, 1.189 (orange, exact), their linear-α forms (dashed) and the
    printed first-power form for α = 0.5 (muted dotted, "as printed"), p in MPa, a pad sketch above; (b) W(α) exact (blue)
    and linear (dashed) for α from −0.9 to 3 with the optimum marked and W < 0 shaded amber for α < 0. Title "A pressure
    hump carries the load — best at inlet/outlet gap ratio 2.19". *see:* "humps peaking just before mid-pad, growing with α;
    the dotted printed curve misses the exact one; the exact load curve bends over and peaks near α = 1.19." *read:* "the
    load is the area under a hump; the linear formula is fine only for α ≲ 0.1 (15 % high at 0.1, 90 % high at 0.5)."
    *change:* "…the film were halved to 25 µm: every pressure and the load grow 4× (∝ 1/h₀²) — the bearing is 'stiff':
    extra load squeezes the film and raises the capacity."
48. `nb.note` — **N35 [B]** "**The exact load for any taper** (ours; the book stops at linear α). Integrating the exact
    pressure: $W=\\frac{6\\mu UL^2}{h_o^2\\alpha^2}\\Big[\\ln(1+\\alpha)-\\frac{2\\alpha}{2+\\alpha}\\Big]$, which tends to
    αμL²U/(2h₀²) as α → 0 (the series of ln(1 + α)). It peaks at an inlet/outlet gap ratio 1 + α = 2.189 — the classic
    optimum of lubrication texts (San Andrés, benchmark V5). With our numbers: α = 0.1 → 12.5 kN/m (linear) vs 10.8 kN/m
    (exact); α = 0.5 → 62.5 vs 32.8 kN/m." equation (W exact) + `nb.code`: `print(ch08.slider_bearing_load(50e-6, 0.5,
    0.05, 5.0, mu=0.05), integrate.quad(lambda xx: ch08.slider_bearing(xx, 50e-6, 0.5, 0.05, 5.0, mu=0.05), 0, 0.05)[0])`
    *expect:* both 3.279e4 N/m.
49. `nb.plotly` — **IF4** `slider_figure(lambda al: {"exact": (x/L, p_exact/1e6), "linear in α": (…), "as printed": (…)},
    "α", np.linspace(-0.8, 3.0, 20), …)`; the title shows W_exact and W_linear from `slider_bearing_state`.
50. `nb.live` — `live(lambda h0_um=50, alpha=0.5, U=5.0, mu=0.05: print(ch08.slider_bearing_state(h0_um*1e-6, alpha, 0.05,
    U, mu=mu)), h0_um=(10, 200, 5), alpha=(-0.8, 3.0, 0.05), U=(-10, 10, 0.5), mu=(0.005, 0.5, 0.005))` — prints p_max, W
    exact and linear and the error (paired with IF4 and E3 for the page).
51. `nb.explainer("slider_bearing", heading="How does a film thinner than a hair carry a load?", why="Drag the taper and
    flip the sliding direction: the pressure hump, the station-by-station Couette-plus-Poiseuille profiles and the load
    change together, the printed-slip ghost visibly departs from the exact hump, and past α = 1 a small recirculation appears
    at the wide inlet — the whole of Example 8.1 in one picture.", tries=["Set α = 0.1 and compare the exact and linear
    humps; then α = 1.19 (the optimum) and read the load.", "Flip U to negative: the hump turns into suction and the status
    says the pad is sucked down.", "Choose 'compare': the printed curve meets both end pressures but sits above the exact one — it fails the pressure equation.", "Push α past 1
    and look at the inlet profiles in the pad frame: a backward-moving layer appears under the pad."])`
52. `nb.md` — **What would change if…** "…there were no upper wall at all — a free surface — and gravity, not a pressure
    applied at the ends, pushed the fluid? The thin-layer balance still holds, with the hydrostatic pressure of the layer
    itself (C08)."
53. `nb.core("C08", "The thin-film equation $\\frac{\\partial h}{\\partial t}=\\frac{\\rho g}{3\\mu}\\frac{\\partial}{\\partial x}
    \\big(h^3\\frac{\\partial h}{\\partial x}\\big)$ (Example 8.3)", question="What equation governs a thin viscous layer
    spreading under its own weight?")`
54. `nb.md` — **The problem in plain words:** "Honey poured on a plate, paint on a wall, a lava lobe, a mud flow, and — with
    a stiffer, nonlinear rheology — an ice sheet: a thin layer of very viscous fluid spreads because its own weight
    presses harder where it is thicker. We want one equation for the thickness h(x, t) that predicts the spreading."
55. `nb.md` — **The idea** (ASCII): "thick places have higher hydrostatic pressure, so fluid is pushed from thick to thin; the
    flux is choked as the layer thins:"
    ```
         h ↑   ____
              /    \      p = ρg(h − y): higher under the crest
             /      \     flux q = −(ρg/3μ) h³ ∂h/∂x   ← h³: a thin film barely moves
    ════════╱════════╲════════
    ```
56. `nb.md` — gloss **stress-free surface and hydrostatic thin layer** (P20 reminder): "At a free surface with negligible air
    drag the shear stress vanishes, μ∂u/∂y = 0; in a thin, slowly varying layer the vertical balance is hydrostatic, p =
    p_a + ρg(h − y) — valid because the slope is small (the lubrication ordering of C05), not merely because accelerations
    are ignored."
57. `nb.derivation("D15", …)` — Part F D15 (10 steps), ref "".
58. `nb.primer("nonlinear diffusion in flux form", "∂h/∂t = ∂/∂x(D(h) ∂h/∂x) is a diffusion equation whose diffusivity
    depends on the unknown — here D = ρgh³/(3μ), large where the layer is thick, tiny where it is thin. Written as
    ∂h/∂t + ∂q/∂x = 0 with the flux q = −D ∂h/∂x, integrating over x shows the total volume ∫h dx never changes when no
    fluid leaves the ends.", code="import numpy as np\nrho, g, mu = 1260.0, 9.81, 1.0             # glycerol\nfor h in (0.01, 0.005, 0.001):                 # layer thickness [m]\n    print(h, rho*g*h**3/(3*mu))                # effective diffusivity [m^2/s] falls as h^3")` (**P191**)
59. `nb.primer("implicit time stepping with Picard iteration", "Explicit steps (P30) of a stiff diffusion equation need a
    tiny time step. Backward Euler evaluates the right-hand side at the *new* time, which is stable for any step; when the
    diffusivity depends on the unknown we guess it from the last iterate, solve the linear system, update the guess and
    repeat until it stops changing (Picard iteration). A thin 'precursor film' h_min keeps the diffusivity positive ahead
    of the front.", code="import numpy as np\nD, dt, dx = 1.0, 0.1, 0.1       # an explicit step would need dt ≤ dx²/(2D) = 0.005\nlam = D*dt/dx**2               # = 10: explicit would blow up\nprint(1/(1 + 4*lam))            # backward-Euler damping of the shortest mode: 0.024 < 1, stable")` (**P192**)
60. `nb.worked_example("a glycerol bead 1 cm high", "ρ = 1260 kg/m³, μ = 1 Pa s, g = 9.81 m/s². 1. Coefficient
    $\\rho g/3\\mu=1260\\times9.81/3=4.12\\times10^3$ m⁻¹s⁻¹. 2. Flux where h = 1 cm and the slope is ∂h/∂x = −0.5:
    $q=-4.12\\times10^3\\times10^{-6}\\times(-0.5)=2.06\\times10^{-3}$ m²/s, flowing outward. 3. Where h = 1 mm with the same
    slope the flux is 1000 times smaller — the edge barely moves. 4. Effective diffusivity at the crest ρgh³/3μ =
    4.1 × 10⁻³ m²/s, falling as h³ while the bead thins.")`
61. `nb.code` — `x = np.linspace(-0.1, 0.1, 400 if not FAST else 200)` · `h0 = np.where(np.abs(x) < 0.01, 0.01, 0.0)` (a 1 cm
    high, 2 cm wide box: half-area 10⁻⁴ m²) · `run = ch08.thin_film_spread(h0, x, [0, 1, 10, 100, 1000], rho=1260.0,
    mu=1.0, cache="outputs/ch08/thin_film_run.npz")` · `print(run["volume"])` · `print(run["x_front"])` ·
    `print(ch08.thin_film_flux(0.01, -0.5, rho=1260.0, mu=1.0))`. *expect:* volume constant to round-off (2.000e-4 m² at
    every time, relative change < 1e-12); fronts ≈ 0.01, …, the late ones close to the similarity law of C10 (x_N ≈ 4.7,
    7.5, 11.8 cm at 10, 100, 1000 s with g = G0; the builder prints both); flux 2.06e-3 m²/s. *explain:* 1. the solver uses
    a conservative finite-volume form (flux differences), so volume is conserved exactly; 2. implicit steps with Picard
    iteration handle the stiffness; 3. the run is cached so the animation and C10 reuse it.
62. `nb.check_agree` — **from scratch (curation §7):** 200 explicit finite-volume steps from the same box on 200 cells: face
    heights by averaging, face flux `-(rho*g/(3*mu))*hf**3*np.diff(h)/dx`, `dt = 0.2*dx**2/(rho*g/(3*mu)*h.max()**3)`,
    `h[1:-1] -= dt/dx*np.diff(F)`; then `assert abs(h.sum()*dx - h0.sum()*dx) < 1e-12*h0.sum()*dx` (volume) and
    `assert np.allclose(h, ch08.thin_film_spread(h0, x, [t_end], rho=1260.0, mu=1.0)["h"][-1], rtol=1e-2, atol=1e-5)`.
63. `nb.animation` — **A3** (frames player, 30 frames FAST 16): h(x, t) (blue fill) from the cached run on a log clock (t from
    0.1 s to 1000 s), the front marker (amber), the volume printed in the title each frame, an inset log–log x_N(t) with the
    slope-1/5 line (purple dashed) appearing once t > 10 s; `player="frames"` so the reader can stop at any time. Figure
    notes: *see:* "a box that slumps into a dome and keeps spreading ever more slowly; the volume never changes." *read:*
    "on the inset the front's points fall on a slope-1/5 line: ten times longer only buys 1.58 times more spread." *change:*
    "…the fluid were ten times more viscous: every time is ten times longer (t enters only as ρgt/μ), the curves are the
    same."
64. `nb.note` — **N37 [C]** "**A similarity solution** — a change of variables that turns this PDE into an ODE — exists for the
    spreading bead; §8.4 develops the idea (C09) and delivers this one (C10: h = At^{−1/5}F(x/Dt^{1/5}))."
65. `nb.explainer("viscous_gravity_current", heading="Why does a honey drop spread like t^(1/5)?", why="Play the spreading
    from any initial shape and watch it forget that shape and lock onto one profile, while the front marches along a
    slope-1/5 line on log–log axes and the volume readout never moves — convergence to self-similarity is a process in
    time, and C10 derives the exponent you measure here.", tries=["Start from 'two humps': they merge, then the profile
    turns into the same dome as the box.", "Watch the log–log front view: compare the slope-1/5 and slope-1/2 guides.",
    "Make the fluid 1000× more viscous (lava lobe) and predict how the curves change before you press play.", "Switch the
    rescaled view to the wrong exponent m = 1/2 and see the collapse fail."])`
66. `nb.md` — **What would change if…** "…a flow had no length or time scale at all — a plate started suddenly in a fluid
    that fills all space? Then the only way to make y dimensionless is to combine it with time, y/√(νt), and the whole
    PDE collapses to an ODE (C09)."

### A.4 §8.4 Similarity Solutions for Unsteady Incompressible Viscous Flow — R14, C09 (+N38–N56 · D16–D20 · P193–P196 · A1 · IF5 · E5), C10 (+N57–N64, N106 · D21 D22 D23 · P197 · A5 · E6)
1. `nb.section("8.4", "Similarity Solutions for Unsteady Incompressible Viscous Flow", intro="**What is this section
   about?** A plate starts moving under still fluid: how far does the motion reach after a time t? Because the problem
   has no length or time scale of its own, y and t can only appear together as y/√(νt); the PDE collapses to an ODE whose
   solution is the error function, and every profile at every time is one curve. We then turn this into a recipe — guess
   a power-law form, demand that every term scales the same way with t — and use it for a thickening shear layer, a
   decaying vortex and the spreading bead of C08. The √(νt) scale is the thickness of every laminar boundary layer
   (Ch. 9) and of the Ekman layers of Ch. 13.")`
2. `nb.recap("R14", "Couette start-up: a flow with a length scale", "Ch. 1 solved the start-up of plane Couette flow — a
   plate suddenly moving at U with a second, fixed plate at distance h — as a Fourier series (`couette_startup_profile`).
   Impulsively started parallel flows keep u ∂u/∂x = 0, so they are exact solutions (Exercise 8.31). The gap h is an
   imposed length: the profiles at different times are *not* one curve. Remove the second plate and that changes (C09).",
   where="Ch. 1 §1.5")`
3. `nb.core("C09", "Stokes' first problem $\\frac uU=1-\\mathrm{erf}\\big(\\frac{y}{2\\sqrt{\\nu t}}\\big)$ (8.30)",
   question="What is u(y, t) above a plate that suddenly starts moving — and why do all the profiles look the same when
   rescaled?")`
4. `nb.md` — **The problem in plain words:** "Yank a plate sideways under a still bath (or start the wind over a calm lake):
   at first only a thin layer moves; the 'news' that the wall moves spreads upward by viscous diffusion. How thick is the
   moving layer after one second, one minute, one hour? This is the prototype of every boundary layer."
5. `nb.md` — **The idea** (ASCII): "no ruler in the problem except √(νt), so the profile can only stretch:"
   ```
   t = 1 s        t = 4 s          t = 16 s        rescaled: y/√(νt)
   │▏             │▎                │▍               │ one curve F(η)
   │▎ thin        │▌ 2× thicker     │█ 4× thicker    │ for every t, U, ν
   ══►U           ══►U              ══►U
   ```
6. `nb.note` — **N38 [B]** "**Set-up (Fig. 8.12; often called Rayleigh's problem).** An infinite plate at y = 0 starts
   moving at U at t = 0 under fluid at rest in y > 0; nothing depends on x, so continuity gives v = 0, and the pressure is
   the same everywhere because far from the plate the fluid is still at rest." **N39 [B]** "What remains is the 1-D
   diffusion equation $\\frac{\\partial u}{\\partial t}=\\nu\\frac{\\partial^2u}{\\partial y^2}$ (8.20) — the heat equation of
   Ch. 1 with ν in place of κ." equation (8.20).
7. `nb.derivation("D16", …)` — Part F D16 (6 steps), ref "8.20".
8. `nb.note` — **N40 [B], N41 [B], N42 [B]** "**Conditions.** $u(y,t=0)=0$ (8.21), $u(y=0,t)=0$ for t < 0 and U for t ≥ 0
   (8.22), $u(y\\to\\infty,t)=0$ (8.23). The problem is well posed: (8.20) has one time derivative and two space
   derivatives, so it needs one condition in t and two in y." gloss **order of a PDE vs number of conditions**.
9. `nb.note` — **N43 [B]** "**Dimensional analysis.** u depends on U, y, t and ν; with two dimensions (L, T) there are three
   groups: $u/U=f(y/\\sqrt{\\nu t},\\,y/Ut)$ (8.24)." **N44 [B]** "**Linearity** removes the second group: $u/U=F(y/\\sqrt{\\nu
   t})\\equiv F(\\eta)$ (8.25) — one variable instead of two." equations (8.24), (8.25).
10. `nb.derivation("D17", …)` — Part F D17 (8 steps), ref "8.25". + `nb.code`: `print(ch01.pi_groups({"u": "m/s", "U":
    "m/s", "y": "m", "t": "s", "nu": "m^2/s"}))` *expect:* three groups equivalent to u/U, y/√(νt), y/(Ut) (the builder
    prints the function's form and checks each is dimensionless).
11. `nb.note` — **N45 [B]** "**Chain rule into the equation:** $\\frac{\\partial u}{\\partial t}=-\\frac{U\\eta}{2t}\\frac{dF}{d\\eta}$
    and $\\frac{\\partial^2u}{\\partial y^2}=\\frac{U}{\\nu t}\\frac{d^2F}{d\\eta^2}$." **N46 [B], N47 [B], N48 [B]** "They give the
    ODE $-\\frac\\eta2\\frac{dF}{d\\eta}=\\frac{d}{d\\eta}\\big(\\frac{dF}{d\\eta}\\big)$ (8.26) with $F(\\eta=0)=1$ (8.27) and
    $F(\\eta\\to\\infty)=0$ (8.28): the initial condition (8.21) and the far condition (8.23) both become η → ∞, so three
    conditions collapse into two — the test that the similarity form is right." equations (8.26)–(8.28).
12. `nb.derivation("D18", …)` — Part F D18 (9 steps), ref "8.26".
13. `nb.primer("Gaussian integral", "The bell curve e^(−ζ²) has total area √π: ∫_{−∞}^{∞} e^(−ζ²) dζ = √π, and by symmetry
    half of it, √π/2, lies on each side. No elementary antiderivative exists, which is why its running integral gets a
    name (next primer).", code="import numpy as np\nfrom scipy.integrate import quad\nprint(quad(lambda z: np.exp(-z**2), 0, np.inf)[0], np.sqrt(np.pi)/2)   # 0.886227 0.886227")` (**P194**)
14. `nb.primer("error function erf and its inverses", "erf(ζ) = (2/√π)∫₀^ζ e^(−ξ²)dξ rises from 0 to 1; erfc = 1 − erf
    (Ch. 4, P123) falls from 1 to 0 and is computed directly without cancellation. Their inverses answer 'where does the
    profile reach this level?': `scipy.special.erfinv`, `erfcinv`.", code="from scipy import special\nprint(special.erf(1.0), special.erfc(1.0))      # 0.8427 0.1573\nprint(2*special.erfcinv(0.01))                 # 3.643: where erfc(η/2) = 1 %")` (**P195**)
15. `nb.note` — **N49 [B], N50 [B], N51 [B]** "**Solving the ODE.** Separating gives $\\frac{dF}{d\\eta}=A\\exp(-\\eta^2/4)$,
    a second integration $F(\\eta)=A\\int_0^\\eta\\exp(-\\xi^2/4)\\,d\\xi+B$ (8.29), and the conditions give B = 1 and, with
    ξ = 2ζ and the Gaussian integral, A = −1/√π." equation (8.29).
16. `nb.derivation("D19", …)` — Part F D19 (11 steps), ref "8.30".
17. `nb.primer("scipy.integrate.solve_bvp", "solve_bvp solves an ODE with conditions at two ends (here F(0) = 1 and
    F(η_max) = 0 on a truncated domain) by collocation on a mesh; we give it the system as first-order equations and a
    guess.", code="import numpy as np\nfrom scipy.integrate import solve_bvp\neta = np.linspace(0, 12, 50)\nsol = solve_bvp(lambda e, Y: np.vstack([Y[1], -e/2*Y[1]]), lambda a, b: np.array([a[0]-1, b[0]]), eta, np.vstack([np.exp(-eta), -np.exp(-eta)]))\nprint(sol.sol(2.0)[0])          # 0.1573 = erfc(1)")` (**P196**)
18. `nb.note` — **N52 [B]** "**Figs. 8.12–8.13 remade:** profiles at several times (raw) and the same against the similarity
    variable. ⚠️ The book's figure axis is y/(2√(νt)) = η/2, while (8.25) defines η = y/√(νt); our axes say which." **N55
    [B]** "**Self-similarity:** profiles at any t, for any U and ν, fall on one curve when u/U is plotted against η — we
    measure the spread below (< 10⁻¹²)."
19. `nb.derivation("D20", …)` — Part F D20 (5 steps), ref "8.31". Then `nb.note("**N54 [B]** The 99 % thickness
    $\\delta_{99}\\sim3.64\\sqrt{\\nu t}$ (8.31) grows like t^{1/2}: in water 3.64 cm after 100 s, 21.9 cm after an hour —
    diffusion is slow over large distances.", equation=r"\delta_{99}\sim3.64\sqrt{\nu t}", ref="8.31")` (the builder prints
    2 erfcinv(0.01) = 3.643 from the code).
20. `nb.worked_example("a plate in water after 100 s", "ν = 10⁻⁶ m²/s, t = 100 s, U = 0.1 m/s. 1. √(νt) = √(10⁻⁴) = 0.01 m.
    2. At y = 1 cm: η = 1, η/2 = 0.5, u/U = erfc(0.5) = 0.4795 — half the plate speed. 3. δ₉₉ = 3.643 × 0.01 = 3.64 cm.
    4. After 400 s: √(νt) doubles to 2 cm, δ₉₉ = 7.29 cm, and u(2 cm) = 0.4795U again (same η). 5. Wall stress
    τ_w = μU/√(πνt) = 10⁻³ × 0.1/√(π × 10⁻⁴) = 5.64 × 10⁻³ Pa, falling like t^{−1/2}.")`
21. `nb.code` — `t = 100.0; y = np.array([0.0, 0.005, 0.01, 0.02, 0.0364])` · `print(ch08.stokes_first_problem(y, t, U=0.1,
    nu=NU_W))` · `print(ch08.similarity_variable(y, t, NU_W), ch08.similarity_variable(y, t, NU_W, half=True))` ·
    `print(ch08.diffusion_thickness(t, NU_W), ch08.diffusion_thickness(3600.0, NU_W))` · `print(ch08.stokes_first_state(t,
    U=0.1, nu=NU_W))` · `print(ch08.vorticity_content(t, U=0.1, nu=NU_W))` · `print(ch04.stokes_first_problem(y, t, 0.1,
    NU_W))` · `eta, F = ch08.similarity_ode_solve("stokes1"); print(np.max(np.abs(F - special.erfc(eta/2))))`. *expect:* u =
    [0.1, 0.0724, 0.0480, 0.0157, 0.0010] m/s; η = [0, 0.5, 1, 2, 3.64]; δ₉₉ = 0.03643 m and 0.2186 m; state: sqrt_nut
    0.01, eta_edge 3.643, delta 0.03643, tau_w 5.64e-3 Pa, vorticity_content 0.1; ∫ω dy = 0.1 m/s = U; ch04 gives the same
    u; solve_bvp error < 1e-7. *explain:* 1. `stokes_first_problem` is (8.30) via erfc; 2. `similarity_variable` has the
    `half` switch for the book's figure axis; 3. `stokes_first_state` bundles the numbers E5 shows; 4. the BVP solver,
    knowing nothing about erf, lands on the same curve.
22. `nb.primer("Crank–Nicolson with scipy.linalg.solve_banded", "Crank–Nicolson averages the diffusion term between the old
    and new time level: second-order accurate in time and stable for any step. Each step solves a tridiagonal system,
    which `scipy.linalg.solve_banded` does in O(N) from the three diagonals stored as rows.", code="import numpy as np\nfrom scipy.linalg import solve_banded\nab = np.array([[0, -1, -1], [4, 4, 4], [-1, -1, 0]], float)   # upper, main, lower diagonals\nprint(solve_banded((1, 1), ab, np.array([3.0, 2.0, 3.0])))  # [1. 1. 1.]")` (**P193**)
23. `nb.check_agree` — **from scratch (curation §7):** a hand-written Crank–Nicolson march on 400 points (200 FAST) up to
    y = 8√(νt_end), dt = t_end/200, the wall held at U after two backward-Euler start-up steps, tridiagonal
    `solve_banded`; `assert np.allclose(u_cn, ch08.stokes_first_problem(yy, t_end, 0.1, NU_W), rtol=1e-3, atol=2e-4)`; also
    `ch08.crank_nicolson_1d(...)` with the same inputs agrees with the hand version to 1e-12, and the FTCS scheme of Ch. 1
    (`diffusion.ftcs_diffusion_1d` at its stable step) agrees to 1e-3. A convergence table (N = 100, 200, 400) prints the
    observed order ≈ 2 (`observed_order`).
24. `nb.figure` — **Figs. 8.12–8.13 remade (N52, N55)** (9 × 3.4 in, three panels): (a) u/U vs y [cm] at t = 10, 30, 100,
    300, 1000 s (blue, darker with time) with δ₉₉ marked (amber ticks); (b) the same five profiles vs η = y/√(νt) — one
    curve, with CN dots (teal) at t = 1000 s, and a secondary axis in η/2 (the book's axis); (c) δ₉₉(t) on log–log axes
    (slope ½) for water and air. Title "One curve for all times: the profile only stretches like √(νt)". *see:* "five
    profiles of growing thickness; in (b) they coincide exactly; in (c) two straight lines of slope ½." *read:* "(b): the
    1 % point sits at η = 3.64 (η/2 = 1.82 on the book's axis); (c): 100 times longer gives 10 times thicker." *change:*
    "…the fluid were air: every thickness grows √15 ≈ 3.9 times faster, but (b) is unchanged."
25. `nb.code` — N55 collapse measured: profiles for (U, ν, t) ∈ {(0.1, 10⁻⁶, 10), (2, 1.5 × 10⁻⁵, 3), (0.01, 10⁻³, 500)}
    evaluated at the same η grid, `spread = max ∣u/U − erfc(η/2)∣`; `assert spread < 1e-12`.
26. `nb.note` — **N53 [B]** "**Reading it as vorticity.** At t = 0 the wall creates a vortex sheet (all the vorticity in an
    infinitely thin layer); afterwards it diffuses out: ω = −∂u/∂y = (U/√(πνt))e^{−y²/4νt} > 0. The total $\\int_0^\\infty
    \\omega\\,dy=u(0)-u(\\infty)=U$ never changes, so no new vorticity is made after t = 0. ⚠️ The book prints this integral as
    −U; since u decreases upward, ω is positive and the integral is +U (checked above by `vorticity_content`). The heat
    analogue: a cold solid whose face is suddenly heated." equation `\int_0^\infty\omega\,dy=U`.
27. `nb.note` — **N56 [B]** "**Similarity needs the absence of imposed scales.** Stop the plate at t = T (Exercise 8.30) and a
    time scale appears; by superposition (ours) $u=U\\big[\\mathrm{erfc}\\frac{y}{2\\sqrt{\\nu t}}-\\mathrm{erfc}\\frac{y}{2
    \\sqrt{\\nu(t-T)}}\\big]$ for t > T. Put a second plate at y = h (R14) and a length scale appears. Either way the
    profiles stop collapsing — spin-down of a stirred tank (Ch. 13) is of this kind." + `nb.figure`: (8 × 3 in, two panels)
    (a) profiles vs η for the stopped plate (T = 100 s) at t = 50, 150, 400 s — no longer one curve; (b) Couette start-up
    from `diffusion.couette_startup_profile(y, t, 0.1, 0.02, NU_W)` vs η at t = 10, 100, 400 s: the curves agree while δ₉₉ ≪
    h and part once the layer feels the far wall. *see / read / change* as usual ("…T → ∞: the stopped-plate curves
    rejoin the erfc curve").
28. `nb.animation` — **A1** (video, 60 frames FAST 30): the plate and fluid on the left with tracers at 8 heights and a dyed
    vertical line bending into the profile; u(y) on the right growing, CN dots riding the erfc curve and the δ₉₉ marker
    rising; log clock t = 1 … 1000 s. Figure notes: *see* "the dye line tilts and bends; the moving layer thickens ever more
    slowly" · *read* "the marker height at t is 3.64√(νt): it needs 4× the time to double" · *change* "…ν ten times
    larger: the same movie, ten times faster".
29. `nb.plotly` — **IF5** `animate_figure(frame_fn, times=np.geomspace(1, 1000, 25), …)` with two subplots: raw profiles
    (left) and rescaled (right) — the collapse without a kernel.
30. `nb.explainer("stokes_first_problem", heading="How can profiles at all times be one curve?", why="Start the plate, watch
    the profile spread with a moving δ₉₉ marker, then flip to the rescaled mode: every earlier profile falls on one curve,
    and changing ν or U does not break it — while stopping the plate or adding a second wall does.", tries=["Play with
    water, then switch to air: the layer grows faster but the rescaled curve does not move.", "Click a point in the
    profile view: read η and the erfc arithmetic of u/U.", "Choose 'stop the plate at T' and watch the collapse fail after
    T.", "Open the Derivation tab at D19 step 8 (ξ = 2ζ): the factor 2 in y/(2√(νt)) appears."])`
31. `nb.md` — **What would change if…** "…the initial state were a velocity jump inside the fluid (two streams), a line
    vortex, or the spreading bead of C08? The same trick works, but the similarity form may need a power of t in front —
    and finding the exponents becomes the whole game (C10)."
32. `nb.core("C10", "The similarity ansatz $\\gamma=At^{-n}F(\\xi/\\delta(t))\\equiv At^{-n}F(\\eta)$ (8.32a) and exponent
    matching", question="How do you find the exponents that make a problem self-similar?")`
33. `nb.md` — **The problem in plain words:** "Two streams slide past each other and the shear layer between them thickens;
    a vortex left alone decays; a honey bead spreads. Each looks self-similar, but with different growth laws (t^{1/2},
    t^{1/5}) and amplitudes that fall with time. We want a procedure that finds those laws instead of guessing them."
34. `nb.md` — **The idea** (table):
    | step | move |
    |---|---|
    | 1 | guess γ = At^{−n}F(ξ/δ(t)) (or Aξ^{−n}F(ξ/δ) if ξ appears in the initial condition) |
    | 2 | substitute into the PDE, divide by one coefficient |
    | 3 | demand every bracketed coefficient is the same power of t → one relation between n and δ(t) |
    | 4 | a conserved quantity (momentum jump, circulation, volume) → the second relation |
    "**Two unknown exponents, two conditions: nothing is left to choose.**"
35. `nb.primer("exponent matching for similarity forms", "If c₁ t^a F(η) + c₂ t^b G(η) = 0 must hold for every t and every η,
    the powers of t must agree (a = b) — otherwise dividing by t^a leaves a t that no function of η can cancel. Power laws
    δ = Dt^m turn every bracket into a power of t, and matching powers gives linear equations for the exponents.",
    code="import sympy as sp\nn, m = sp.symbols('n m')\nprint(sp.solve([sp.Eq(-n - 1, -4*n - 2*m), sp.Eq(-n + m, 0)], [n, m]))   # {n: 1/5, m: 1/5}")`
    (**P197**)
36. `nb.note` — **N57 [B]** "**The second form** $\\gamma=A\\xi^{-n}F(\\xi/\\delta(t))\\equiv A\\xi^{-n}F(\\eta)$ (8.32b) puts a power
    of the space coordinate in front; use it when ξ appears in the initial or boundary condition (the line vortex below,
    Γ/2πr)." equation (8.32b).
37. `nb.derivation("D21", …)` — Part F D21 (8 steps), ref "8.32". Then `nb.note("**N58 [B]** Example 8.4 recovers Stokes'
    first problem from the ansatz: A = 1, n = 0, δδ′ = C₁ν, δ = √(2C₁νt); C₁ = ½ gives η of (8.25).")`.
38. `nb.md` — gloss **spotting an exact derivative** (P38 reminder): "F + η dF/dη = d(ηF)/dη — the product rule read
    backwards."
39. `nb.derivation("D22", …)` — Part F D22 (14 steps, ★★★), ref "". **check_src:** `vs = ch08.similarity_reduce_sympy(
    "vortex_sheet")  # D22 in sympy` · `print(vs["brackets"])  # n/t, δ'/δ, ν/δ² with δ = √(νt)` · `print(vs["n"])  # 1/2
    from the conserved jump` · `assert sp.simplify(vs["residual"]) == 0  # ω = −U/√(πνt) e^{−y²/4νt} solves ∂ω/∂t = ν∂²ω/∂y²`
    · `y, t, U, nu = sp.symbols("y t U nu", positive=True)` · `omega = -U/sp.sqrt(sp.pi*nu*t)*sp.exp(-y**2/(4*nu*t))` ·
    `print(sp.simplify(-sp.integrate(omega, (y, -sp.oo, sp.oo))))  # 2U at every t: the jump is conserved` — every line
    commented.
40. `nb.note` — **N59 [B]** "**Example 8.5, the viscous vortex sheet.** $\\omega_z(y,t)=-\\frac{U}{\\sqrt{\\pi\\nu t}}\\exp\\{-
    \\frac{y^2}{4\\nu t}\\}$ and $u(y,t)=U\\,\\mathrm{erf}\\{\\frac{y}{2\\sqrt{\\nu t}}\\}$. Define the layer's width by u = ±0.95U:
    η = ±2 erfinv(0.95) = ±2.772 and the width is 5.544√(νt). ⚠️ The book prints ±2.76 but also the width 5.54 — the 2.76 is
    a rounding slip. ⚠️ Ch. 5's `diffusing_vortex_sheet(y, t, gamma, nu)` takes γ = u_below − u_above = −2U." equations
    (ω_z, u) + `nb.code`: `u, w = ch08.vortex_sheet_diffusion(0.001, 1.0, U=0.01, nu=NU_W)` · `print(u, w,
    ch05.diffusing_vortex_sheet(0.001, 1.0, -0.02, NU_W))` · `print(ch08.transition_width(1.0, NU_W), 2*special.erfinv(0.95))`.
    *expect:* u = 0.005205 m/s (erf(0.5) × 0.01), ω = −4.394 s⁻¹, ch05 equal; width 5.544e-3 m, 2.7718.
41. `nb.note` — **N60 [B]** "**A temporally developing boundary layer.** Shift the upper half of Example 8.5 by −U and flip
    the sign: it is Stokes' first problem. Its wall stress and skin friction $\\tau_w=\\mu\\big(\\frac{\\partial u}{\\partial y}
    \\big)_{y=0}=\\frac{\\mu U}{\\sqrt{\\pi\\nu t}}$, $C_f=\\frac{\\tau_w}{\\frac12\\rho U^2}=\\frac{2}{\\sqrt\\pi}\\sqrt{\\frac{\\nu}
    {U^2t}}$ become, with Ut read as a distance x, $C_f=1.128\\,\\mathrm{Re}_x^{-1/2}$ — the Re_x^{−1/2} law of a laminar
    boundary layer on a plate (Ch. 9 finds 0.664 for the real, spatially growing layer: same exponent, different
    constant)." equations + `nb.code`: `print(ch08.temporal_bl_wall_stress(1e4*NU_W/1.0**2, 1.0, NU_W))` *expect:* Rex =
    1e4, Cf = 0.01128 (Blasius would give 0.00664).
42. `nb.note` — **N61 [B]** "**Example 8.6, a decaying line vortex** (stated; the same moves as D22 with the form (8.32b),
    Ar^{−n} = Γ/2πr). The similarity ODE $(\\frac1\\eta-\\frac\\eta2)F'=F''$ gives F = 1 − e^{−η²/4} and $u_\\theta(r,t)=\\frac{
    \\Gamma}{2\\pi r}\\big[1-\\exp\\{-\\frac{r^2}{4\\nu t}\\}\\big]$ — the Gaussian (Lamb–Oseen) vortex (3.29) with σ² = 4νt (Ch. 3,
    Ch. 5): rigid rotation inside r ≈ 2.24√(νt), an ideal vortex outside." **N62 [B]** "**Exercise 8.26, a vortex switched
    on:** $u_\\theta=\\frac{\\Gamma}{2\\pi r}\\exp\\{-\\frac{r^2}{4\\nu t}\\}$ — the difference between the steady ideal vortex and
    the decaying one." equations + `nb.code` (sympy residual cell, every line commented): `r, t, G, nu = sp.symbols("r t
    Gamma nu", positive=True)` · for each `u` in the two forms: `res = sp.diff(u, t) - nu*sp.diff(sp.diff(r*u, r)/r, r)`;
    `print(sp.simplify(res))` → 0, 0 · `print(np.allclose(ch08.line_vortex_decay(0.01, 100.0, 1e-3, NU_W),
    VX.gaussian_vortex(0.01, 1e-3, 2*np.sqrt(NU_W*100.0))))` → True.
43. `nb.derivation("D23", …)` — Part F D23 (9 steps), ref "". Then `nb.note("**N63 [B]** Example 8.7: the spreading bead of
    C08 is self-similar with n = m = 1/5: $h(x,t)=At^{-1/5}F(x/Dt^{1/5})$. Solving the ODE for F (not done in the book) gives
    Huppert's (1982) dome F ∝ (1 − η²/η_N²)^{1/3} and the front x_N = η_N(βA³t)^{1/5}, η_N = 1.411 (benchmark V5). ⚠️ The
    book's last line says 'the final equation of Example 8.2'; it means Example 8.3.")`.
44. `nb.worked_example("a 1 cm/s shear layer after one second", "U = 0.01 m/s on each side, ν = 10⁻⁶ m²/s, t = 1 s. 1.
    √(νt) = 10⁻³ m. 2. Width 5.544 × 1 mm = 5.5 mm. 3. At y = 1 mm: η = 1, u = U erf(0.5) = 0.0052 m/s. 4. Peak vorticity
    −U/√(πνt) = −0.01/1.772 × 10⁻³ = −5.64 s⁻¹. 5. After 100 s: width 5.5 cm, peak −0.564 s⁻¹, and −∫ω dy = 2U still.")`
45. `nb.code` — `for case in ("stokes1", "vortex_sheet", "line_vortex", "spreading"): s = ch08.similarity_reduce_sympy(case);
    print(case, s["n"], s["m"] if "m" in s else s["delta"])` · `print(ch08.similarity_collapse_error("spreading", 0.2, 0.2,
    [10, 30, 100, 300, 1000]), ch08.similarity_collapse_error("spreading", 0.25, 0.5, [10, 30, 100, 300, 1000]))`.
    *expect:* stokes1 n = 0, δ ∝ √t; vortex_sheet n = 1/2; line_vortex n = 1 (space power) with δ = √(νt); spreading
    n = m = 1/5; collapse spread ≈ 1e-15 at (0.2, 0.2) and O(0.1) at the wrong pair.
46. `nb.check_agree` — **from scratch (curation §7):** log–log fits: the sheet width at t = 1, 10, 100, 1000 s from
    `vortex_sheet_diffusion` located by `brentq` (u = 0.95U) → slope 0.5; the bead front `run["x_front"]` of C08's cached
    run for t ≥ 10 s → slope 0.2; `assert abs(slope_sheet - 0.5) < 5e-3 and abs(slope_bead - 0.2) < 5e-3`.
47. `nb.figure` — **Figs. 8.14–8.15 remade (N106 [B])** (10 × 3.4 in, three panels): (a) ω_z(y) of the thickening sheet at
    t₁ = 1 s and t₂ = 4 s (purple) — the peak halves, the width doubles; (b) u/U vs y/(2√(νt)) (the book's axis, labelled)
    — one curve for both times; (c) the decaying line vortex u_θ(r) at νt = 0 (ideal vortex, dashed), and three later
    times (our values, not the book's), with the solid-body cores marked. Title "Three self-similar diffusions". *see / read /
    change* ("…Γ doubled: every curve in (c) doubles; its shape and the core radius 2.24√(νt) do not change").
48. `nb.animation` — **A5** (video, 60 frames FAST 30): left, the vortex sheet thickening (u and ω); right, the line vortex
    decaying; in the last 15 frames both switch to rescaled axes and collapse. Notes *see / read / change* as usual.
49. `nb.explainer("similarity_exponents", heading="Where do similarity exponents come from?", why="Set trial exponents n and m
    with two sliders: the rescaled profiles collapse onto one curve only at the right values, and the powers of t in the
    reduced equation turn green when they match — the exponent algebra of Examples 8.4–8.7 as a game you can win or lose.",
    tries=["In 'vortex sheet' mode, start at n = 0.4: the curves fan out; slide to 0.5 and watch them merge.", "In
    'spreading bead' mode try the diffusion guess m = 1/2 — the bracket chip turns red.", "Switch to 'line vortex': why is
    the prefactor now a power of r, not of t?", "Check the conserved-integral readout: it is flat only at the right n."])`
50. `nb.note` — **N64 [C]** "Diffusive lengths grow like (νt)^{1/2}; the bead's t^{1/5} is not a diffusion length but an
    advective one — how far the fluid itself has travelled."
51. `nb.md` — **What would change if…** "…the wall did not start once but oscillated forever? Then the problem has an imposed
    time 1/ω; no similarity variable exists, but the diffusion length √(ν/ω) sets the depth (C11)."

### A.5 §8.5 Flow Due to an Oscillating Plate — C11 (+N65–N74 · D24 D25 · A2 · IF6 · E7)
1. `nb.section("8.5", "Flow Due to an Oscillating Plate", intro="**What is this section about?** A plate shaking back and
   forth in its own plane drags a thin layer of fluid with it. The motion decays and lags with height, looking like a wave
   travelling away from the wall — but it is diffusion. The depth √(2ν/ω) and the phase structure (1 + i)y/δ are exactly
   those of tidal bottom boundary layers and of the Ekman layer of Ch. 13.")`
2. `nb.core("C11", "Stokes' second problem $u=U\\exp\\{-y\\sqrt{\\frac{\\omega}{2\\nu}}\\}\\cos(\\omega t-y\\sqrt{\\frac{\\omega}{2\\nu}})$
   (8.38)", question="How deep does an oscillation reach into a viscous fluid — and why do the deeper layers lag behind?")`
3. `nb.md` — **The problem in plain words:** "Shake a plate under water at one cycle per second; the tide rubs the sea
   floor back and forth twice a day; a sound wave slides air along a wall a thousand times a second. How far above the
   wall does the fluid feel the shaking, and how late does each layer follow?"
4. `nb.md` — **The idea** (ASCII + table): "each layer repeats the wall's motion, delayed by y/δ radians and shrunk by
   e^{−y/δ}:"
   ```
   y = 3δ   ~~ 5 %,  3 rad late
   y = 2δ   ~~~ 14 %, 2 rad late
   y = δ    ~~~~ 37 %, 1 rad late
   y = 0    ◄══►  U cos ωt           δ = δ_e = √(2ν/ω)
   ```
   | shaking | ν [m²/s] | δ_e |
   |---|---|---|
   | plate at 1 Hz in water | 10⁻⁶ | 0.56 mm |
   | M₂ tide, molecular ν | 10⁻⁶ | 0.12 m |
   | M₂ tide, eddy ν (turbulent) | 10⁻² | 12 m |
   | sound at 1 kHz in air | 1.5 × 10⁻⁵ | 0.069 mm |
5. `nb.note` — **N65 [B]** "**Set-up.** The plate at y = 0 oscillates in its own plane; we want only the periodic state after
   start-up transients have died, so there is no initial condition. The equation is (8.20) again, $\\frac{\\partial u}{\\partial t}
   =\\nu\\frac{\\partial^2u}{\\partial y^2}$, but now the wall imposes a time scale 1/ω." **N66 [B], N67 [B]** "$u(y=0,t)=U
   \\cos(\\omega t)$ (8.33) and $u(y\\to\\infty,t)$ bounded (8.34)." ⚠️ "ω is now a frequency [rad/s] (as in Ch. 7), not a
   vorticity."
6. `nb.md` — gloss **complex exponential solution** (P176 complex amplitudes, P159 complex square roots): "Write the wall
   motion as Re{Ue^{iωt}}; for a linear equation with real coefficients we can solve with the complex form and take the
   real part at the end. √i = (1 + i)/√2, since ((1 + i)/√2)² = i." + `nb.code`: `print(np.sqrt(1j), (1+1j)/np.sqrt(2))`
   *expect:* (0.7071+0.7071j) twice.
7. `nb.note` — **N68 [B]** "**The complex trial** $u(y,t)=\\mathrm{Re}\\{e^{i\\omega t}f(y)\\}$ (8.35). ⚠️ The text then says
   'substitution of (8.33) into (8.20)'; it is (8.35) that is substituted." **N69 [B]** "It turns (8.20) into
   $i\\omega f=\\nu(d^2f/dy^2)$ (8.36)." **N70 [B]** "$f=e^{ky}$ with $k=(i\\omega/\\nu)^{1/2}=\\pm(i+1)(\\omega/2\\nu)^{1/2}$."
   **N71 [B]** "$f(y)=A\\exp\\{-(i+1)y\\sqrt{\\omega/2\\nu}\\}+B\\exp\\{+(i+1)y\\sqrt{\\omega/2\\nu}\\}$ (8.37); bounded ⇒ B = 0,
   wall ⇒ A = U." equations (8.35)–(8.37).
8. `nb.derivation("D24", …)` — Part F D24 (10 steps), ref "8.38". + `nb.code`: `print(ch08.stokes_second_sympy()["residual"])`
   *expect:* 0.
9. `nb.note` — **N72 [B]** "**Reading (8.38).** The cosine moves toward +y — it looks like a damped transverse wave — but there
   is no restoring force: it is diffusion driven by a periodic wall. At y = 4(ν/ω)^{1/2} the amplitude is
   Ue^{−4/√2} = 0.059U, so the book calls δ ~ 4(ν/ω)^{1/2} the layer thickness; the literature's e-folding depth is
   δ_e = √(2ν/ω), 2√2 = 2.83 times smaller (ours: crests move at √(2νω), 'wavelength' 2πδ_e)." (Builder: print
   e^{−2√2} = 0.05911 from the code.)
10. `nb.derivation("D25", …)` — Part F D25 (6 steps), ref "".
11. `nb.worked_example("a plate shaken once a second in water", "ω = 2π rad/s, ν = 10⁻⁶ m²/s, U = 0.1 m/s. 1.
    δ_e = √(2ν/ω) = √(2 × 10⁻⁶/6.283) = 5.64 × 10⁻⁴ m = 0.56 mm. 2. Book depth 4√(ν/ω) = 1.60 mm. 3. At y = 1 mm:
    amplitude e^{−1/0.564} = e^{−1.77} = 0.17 → 1.7 cm/s; phase lag 1.77 rad, i.e. 1.77/(2π) × 1 s = 0.28 s late. 4.
    Crest speed ωδ_e = 3.5 mm/s. 5. Shake 100× faster: every length shrinks 10×.")`
12. `nb.code` — `y = np.array([0.0, 0.5e-3, 1e-3, 2e-3])` · `print(ch08.stokes_second_problem(y, 0.0, U=0.1,
    omega=2*np.pi, nu=NU_W))` · `print(ch08.stokes_layer(NU_W, 2*np.pi))` · `print(ch08.stokes_layer_state(NU_W, 2*np.pi,
    1e-3))` · `wM2 = 2*np.pi/(12.42*3600); print(ch08.stokes_layer(1e-6, wM2)["delta_e"], ch08.stokes_layer(1e-2,
    wM2)["delta_e"])` · `print(ch08.pi_groups({"u": "m/s", "U": "m/s", "y": "m", "t": "s", "nu": "m^2/s", "omega": "1/s"}))`
    (via `ch01.pi_groups`). *expect:* u(t = 0) = [0.1, 0.0260, −0.0034, −0.0027] m/s; delta_e 5.642e-4, delta_book 1.596e-3, phase_speed 3.545e-3 m/s, wavelength 3.545e-3 m; at 1 mm amplitude
    0.170, phase_lag 1.772 rad, time_lag 0.282 s; tide 0.119 m and 11.9 m; four groups, one of which contains U (e.g. yω/U) and is removed by linearity — **N73 [B]** "three groups u/U, ωt, y(ω/ν)^{1/2}: y and t never combine, so no similarity; the extent
    (ν/ω)^{1/2} is viscosity times the imposed period."
13. `nb.check_agree` — **from scratch (curation §7):** `de = np.sqrt(2*NU_W/(2*np.pi)); t = 0.3; u_mine = np.real(0.1 *
    np.exp(1j*2*np.pi*t) * np.exp(-(1+1j)*y/de))`; `assert np.allclose(u_mine, ch08.stokes_second_problem(y, t, U=0.1,
    omega=2*np.pi, nu=NU_W), rtol=1e-12, atol=1e-15)`. Plus a numerical twin: `ch08.crank_nicolson_1d` with the wall
    U cos ωt, run for 10 periods (5 FAST), matches (8.38) to 2 × 10⁻³ U over 0 ≤ y ≤ 6δ_e.
14. `nb.figure` — **Fig. 8.16 remade (N107 part)** (6 × 4 in): u/U vs y/δ_e (the book's axis y√(ω/2ν)) at ωt = 0, π/2, π, 3π/2
    (blue shades) inside the envelope ±e^{−y/δ_e} (muted dashed), the book depth 4√(ν/ω) = 2√2 δ_e marked (amber) with its
    0.059 amplitude. Title "Each layer is the wall's motion, delayed and shrunk". *see:* "four S-shaped curves inside a
    narrowing funnel." *read:* "at any height the curves' spread is the local amplitude; the height where a curve crosses
    zero moves up with ωt — that is the apparent 'wave'." *change:* "…the fluid were ten times more viscous: all lengths grow
    √10 ≈ 3.2 times; in these rescaled axes nothing moves."
15. `nb.animation` — **A2** (video, 60 frames FAST 30, two periods): left, tracers at 6 heights oscillating with growing lag and
    a dyed vertical line bending; right, the profile swinging inside its envelope with a probe dot at y = δ_e. Notes *see /
    read / change* ("…the probe moved to 2δ_e: its dot moves e⁻¹ as far and peaks 1 rad later").
16. `nb.plotly` — **IF6** `animate_figure(lambda wt: {"u/U": (…), "envelope": (…)}, np.linspace(0, 2*np.pi, 25), …)` — the
    page-surviving A2.
17. `nb.explainer("oscillating_plate", heading="How deep does a shaking plate reach?", why="The animation shows the phase lag
    and the envelope together; drag ω and the layer shrinks while the 0.059U marker at 4√(ν/ω) follows, and a probe clock
    shows a deep layer peaking later — none of which a four-phase static figure makes vivid.", tries=["Play the 1 Hz preset
    and click at y = 1 mm: read the amplitude and the time lag.", "Switch to the M₂ tide with molecular and then eddy
    viscosity: 12 cm versus 12 m.", "Double ω and predict δ_e before you look (it falls by √2).", "Open the 'Right now'
    notes: which real Stokes layer is closest to your setting?"])`
18. `nb.note` — **N74 [C]** "The same layer explains why sound is weakly absorbed at flat walls: the viscous acoustic boundary
    layer (Ch. 15). For Shammunul's thread: replace the plate's oscillation by the Coriolis force and the (1 + i)/δ structure
    of (8.37) becomes the Ekman spiral (Ch. 13)."
19. `nb.md` — **What would change if…** "…the plate stood still and a sphere moved slowly through the fluid instead? Then
    there is no single wall-normal direction and the flow is 3-D; for tiny Reynolds numbers inertia vanishes altogether and
    the equations become linear again (C12)."

### A.6 §8.6 Low Reynolds Number Viscous Flow Past a Sphere — R15, R16, C12 (+N75–N80 · D26 · P198 · IF7), R17, R18, R19, C13 (+N81–N87, N93, N107 · D27 D28 D29 · P199 · A4 · E8), R20, C14 (+N88–N92 · D30 D31 · E9), C15 (+N94–N101 · D32 D33 · IF8)
1. `nb.section("8.6", "Low Reynolds Number Viscous Flow Past a Sphere", intro="**What is this section about?** Tiny,
   slow things — cloud droplets, silt, bacteria, a bead sinking in syrup — live at Reynolds numbers far below one. We find
   the right way to drop inertia (pressure must be rescaled, not dropped), solve the resulting linear Stokes equations for
   a sphere, integrate the surface stresses to Stokes' drag D = 6πμaU, and use it for the fall speed of droplets and for
   Millikan's measurement of the electron's charge. Finally we see why the solution fails far from the sphere and how
   Oseen repaired it.")`
2. `nb.recap("R15", "Steady Navier–Stokes", "For steady flow round a body of size L at speed U: $\\rho\\mathbf u\\cdot\\nabla
   \\mathbf u+\\nabla p=\\mu\\nabla^2\\mathbf u$ (8.39) — (4.39b) without the time derivative and with gravity absorbed into
   p.", where="Ch. 4 §4.6")`
3. `nb.recap("R16", "The high-Re scaling", "With u* = u/U, x* = x/L and p* = (p − p∞)/ρU² (4.100), (8.39) becomes
   $\\mathbf u^*\\cdot\\nabla^*\\mathbf u^*+\\nabla^*p^*=\\frac1{\\mathrm{Re}}\\nabla^{*2}\\mathbf u^*$ (8.40), Re = ρUL/μ;
   for Re ≫ 1 the viscous term looks negligible (Euler, Ch. 6).", where="Ch. 4 §4.11")` + `nb.code`:
   `print(SIM.nondimensional_ns_coefficients("dynamic", "advective"))` · `print(SIM.nondimensional_ns_coefficients("viscous",
   "advective"))` *expect:* viscous coefficient μ/(ρUl) = 1/Re with the dynamic scale; pressure coefficient μ/(ρUl) with the
   viscous scale.
4. `nb.core("C12", "Creeping flow: the Stokes equations $\\nabla p=\\mu\\nabla^2\\mathbf u$ (8.43)", question="Which terms survive
   when Re → 0 — and why does pressure not disappear together with inertia?")`
5. `nb.md` — **The problem in plain words:** "A bacterium swims at a few body lengths per second; a grain of silt sinks
   through a river; a mist droplet drifts down. For them Re = UL/ν ≪ 1. Inertia is negligible — but if we simply
   multiply (8.40) by Re and let it go to zero, the pressure disappears too, and the resulting equation cannot hold a
   sphere in a stream. What went wrong, and what is the right limit?"
6. `nb.md` — **The idea** (table): "the size of pressure differences is set by whatever they balance:"
   | regime | pressure balances | pressure scale | dimensionless form |
   |---|---|---|---|
   | Re ≫ 1 | inertia | ρU² | (8.40) |
   | Re ≪ 1 | viscous stress | μU/L | (8.42) |
7. `nb.note` — **N75 [C]** "Many problems are solved as expansions in a small or large parameter (here Re); Ch. 9 (1/Re) and
   Ch. 11 (small disturbances) do the same." **N76 [C]** "At high Re the expansion in 1/Re fails near walls, where the
   inviscid solution cannot satisfy no slip — the boundary layer of Ch. 9; the two-length scaling (8.14) of C05 was the
   hint."
8. `nb.primer("dominant balance", "Choose the scale of a quantity from the term it must balance. If pressure differences are
   pushed by viscous stresses, then p − p∞ ~ μU/L, not ρU²; picking the wrong scale makes a term look negligible when it
   is not.", code="mu, U, L, rho = 1e-3, 1e-3, 1e-5, 1000.0   # a 10 µm particle at 1 mm/s in water\nprint(rho*U**2, mu*U/L)                    # dynamic 1e-3 Pa vs viscous 0.1 Pa: pressure is viscous here")`
   (**P198**)
9. `nb.note` — **N77 [B]** "Multiplying (8.40) by Re gives $\\mathrm{Re}(\\mathbf u^*\\cdot\\nabla^*\\mathbf u^*+\\nabla^*p^*)=
   \\nabla^{*2}\\mathbf u^*$ (8.41): as Re → 0 it leaves 0 = μ∇²u — pressure lost too." **N78 [B]** "The low-Re pressure scale
   p* = (p − p∞)L/(μU)." **N79 [B]** "It gives $\\mathrm{Re}(\\mathbf u^*\\cdot\\nabla^*\\mathbf u^*)=-\\nabla^*p^*+\\nabla^{*2}
   \\mathbf u^*$ (8.42), whose Re → 0 limit is (8.43)." equations (8.41), (8.42).
10. `nb.derivation("D26", …)` — Part F D26 (9 steps), ref "8.43".
11. `nb.worked_example("a 10 µm cloud droplet", "a = 10 µm, falling at U ≈ 1.2 cm/s in air (ν = 1.5 × 10⁻⁵ m²/s, μ = 1.81 ×
    10⁻⁵ Pa s; the speed is derived in C14). 1. Re = 2aU/ν = 2 × 10⁻⁵ × 0.012/1.5 × 10⁻⁵ = 0.016. 2. Dynamic pressure
    ρU² = 1.2 × 1.44 × 10⁻⁴ = 1.7 × 10⁻⁴ Pa. 3. Viscous pressure μU/a = 1.81 × 10⁻⁵ × 0.012/10⁻⁵ = 0.022 Pa — 125 times
    larger (= 2/Re). The viscous scale is the right one.")`
12. `nb.code` — `lr = ch08.low_re_scaling_sympy()` · `print(lr["dynamic"], lr["dynamic_times_Re"], lr["viscous"])` · `u_fn =
    lambda x, t=0.0: np.stack(ch08.stokes_sphere_velocity_xyz(x[0], x[1], x[2], U=1e-3, a=1e-5))` · `p_fn = lambda x, t=0.0:
    ch08.stokes_sphere_pressure(np.sqrt((x**2).sum(0)), np.arccos(x[0]/np.sqrt((x**2).sum(0))), U=1e-3, a=1e-5, mu=MU_W)` ·
    `print(ch08.stokes_residual(u_fn, p_fn, np.array([3e-5, 1e-5, 0.5e-5]), mu=MU_W))`. *expect:* dynamic {inertia 1,
    pressure 1, viscous 1/Re}; ×Re {Re, Re, 1}; viscous {Re, 1, 1}; Stokes residual ≈ 0 (≪ the pressure gradient, relative
    < 1e-6) — the sphere field of C13 satisfies (8.43). *explain:* 1. the sympy engine repeats D26's coefficients; 2.
    `stokes_residual` evaluates ∇p − μ∇²u by central differences.
13. `nb.check_agree` — **from scratch:** a coefficient table by hand, `Re = 0.016; dyn = [1, 1, 1/Re]; visc = [Re, 1, 1]`, and
    `assert np.allclose(dyn, [float(c.subs(sp.Symbol("Re"), Re)) for c in lr["dynamic"].values()])` (same for `visc`).
14. `nb.figure` — (7 × 3.4 in) the three coefficients (inertia teal, pressure orange, viscous rose) vs Re on log–log axes
    (10⁻³ … 10³): left panel the dynamic scaling ×Re (8.41) — pressure falls with inertia; right panel the viscous scaling
    (8.42) — pressure and viscous flat at 1, inertia ∝ Re. Title "Rescale pressure, not just inertia". *see / read / change*
    ("…Re = 1: the two scalings coincide — neither term is negligible").
15. `nb.plotly` — **IF7** `slider_figure(…, "log10 Re", np.linspace(-3, 3, 25), …)` bars of the three coefficients for both
    scalings, the page-surviving version of the figure.
16. `nb.note` — **N80 [C]** "The book's italic rule: the right length and time scales depend on the region of the flow and are
    found by balancing the terms that matter there. Ch. 9 (boundary layers) and Ch. 13 (rotating, stratified scales) live by
    it."
17. `nb.md` — **What would change if…** "…we solve (8.43) round a sphere? The equation is linear, so reversing the stream
    reverses the whole flow: no wake (C13)."
18. `nb.recap("R17", "Vorticity in spherical coordinates", "For an axisymmetric flow without swirl only the azimuthal vorticity
    survives, $\\omega_\\varphi=\\frac1r\\big[\\frac{\\partial(ru_\\theta)}{\\partial r}-\\frac{\\partial u_r}{\\partial\\theta}\\big]$
    (the spherical curl of `core.curvilinear`).", where="Ch. 3 §3.4")`
19. `nb.recap("R18", "Stokes' stream function", "Axisymmetric velocities come from one function: $u_r=\\frac{1}{r^2\\sin\\theta}
    \\frac{\\partial\\psi}{\\partial\\theta}$, $u_\\theta=-\\frac{1}{r\\sin\\theta}\\frac{\\partial\\psi}{\\partial r}$ (6.83).",
    where="Ch. 6 §6.8")`
20. `nb.recap("R19", "The uniform stream far away", "A uniform stream U along the axis has $\\psi=\\frac12Ur^2\\sin^2\\theta$
    (6.86); the far-field condition is $\\psi(r\\to\\infty,\\theta)=\\frac12Ur^2\\sin^2\\theta$ (8.47).", where="Ch. 6 §6.8")`
21. `nb.core("C13", "Stokes' solution for the sphere $\\psi=Ur^2\\sin^2\\theta\\big(\\frac12-\\frac{3a}{4r}+\\frac{a^3}{4r^3}\\big)$
    (8.48)", question="What is the flow round a sphere when viscosity dominates everywhere?")`
22. `nb.md` — **The problem in plain words:** "A tiny bead sinks through syrup. Drawn in the bead's frame, the fluid streams past
    it; drawn in the syrup's frame, the fluid is pushed aside in front and closes in behind. We want the whole velocity
    field — it gives the pressure, the stresses and the drag (C14), and it is the reference for every settling particle and
    swimming micro-organism."
23. `nb.md` — **The idea** (ASCII): "a uniform stream minus a slowly decaying disturbance, fore–aft symmetric:"
    ```
    ─────────►        ψ = ½Ur² sin²θ  −  (3a/4)Ur sin²θ  +  (a³/4)(U/r) sin²θ
    ────╮ ● ╭─►          stream          'Stokeslet' ~ a/r      'doublet' ~ a³/r³
    ─────────►        the a/r part decays so slowly that the sphere is felt tens of radii away
    ```
24. `nb.note` — **N81 [B]** "**Take the curl.** The curl of (8.43) kills the pressure (curl of a gradient) and leaves ∇²ω = 0
    — in Cartesian components; in spherical coordinates the operator means −∇×∇×ω (the book's footnote)." gloss **curl of a
    gradient is zero; curl and the Cartesian Laplacian commute** (P121, P122 reminders).
25. `nb.derivation("D27", …)` — Part F D27 (6 steps), ref "".
26. `nb.primer("the Stokes operator E² applied twice", "For axisymmetric flow the operator E² = ∂²/∂r² + (sin θ/r²)∂/∂θ(1/sin θ
    ∂/∂θ) (Ch. 6 used it once, (6.77)) turns ψ into the vorticity: ω_φ = −E²ψ/(r sin θ). Applying it twice, E²(E²ψ), is not
    the biharmonic ∇⁴ψ — a sympy test shows the difference.", code="import sympy as sp\nr, th = sp.symbols('r theta', positive=True)\nE2 = lambda f: sp.diff(f, r, 2) + sp.sin(th)/r**2*sp.diff(sp.diff(f, th)/sp.sin(th), th)\nprint(sp.simplify(E2(r**2*sp.sin(th)**2)))          # 0: the uniform stream has no vorticity")`
    (**P199**)
27. `nb.note` — **N82 [B]** "**Vorticity from ψ:** $\\omega_\\varphi=-\\frac1r\\big[\\frac{1}{\\sin\\theta}\\frac{\\partial^2\\psi}{\\partial
    r^2}+\\frac1{r^2}\\frac{\\partial}{\\partial\\theta}\\big(\\frac1{\\sin\\theta}\\frac{\\partial\\psi}{\\partial\\theta}\\big)\\big]=-
    \\frac{E^2\\psi}{r\\sin\\theta}$." **N83 [B]** "Combining with ∇²ω_φ = 0 in the footnote's sense gives $\\big[\\frac{\\partial^2}
    {\\partial r^2}+\\frac{\\sin\\theta}{r^2}\\frac{\\partial}{\\partial\\theta}\\big(\\frac1{\\sin\\theta}\\frac{\\partial}{\\partial\\theta}
    \\big)\\big]^2\\psi=0$ (8.44). ⚠️ It is the square of E², not the biharmonic (footnote 2)." equation (8.44).
28. `nb.derivation("D28", …)` — Part F D28 (12 steps, ★★★), ref "8.44". **check_src:** `S = ch08.stokes_sphere_sympy()  # D27–D31
    in one cached object` · `print(S["omega_phi"])  # −(3Ua/2r²) sin θ` · `assert sp.simplify(S["curlcurl_identity"]) == 0  #
    steps 7–9: curl-curl of A_φ e_φ equals −E²(r sinθ A_φ)/(r sinθ) e_φ` · `assert sp.simplify(S["E4psi"]) == 0  # (8.44)
    holds for (8.48)` · `assert sp.simplify(S["biharmonic_psi"]) != 0  # the biharmonic of (8.48) is not zero` — every line
    commented.
29. `nb.note` — **N84 [B], N85 [B]** "**Conditions on the sphere:** $\\psi(r=a,\\theta)=0$ (8.45) (no flow through it) and
    $\\partial\\psi(r=a,\\theta)/\\partial r=0$ (8.46) (no slip); far away (8.47)." **N86 [B]** "(8.47) suggests ψ = f(r) sin²θ;
    (8.44) becomes $f^{iv}-\\frac{4f''}{r^2}+\\frac{8f'}{r^3}-\\frac{8f}{r^4}=0$, solved by f = Ar⁴ + Br² + Cr + D/r with
    A = 0, B = U/2, C = −3Ua/4, D = Ua³/4." equations (8.45)–(8.47), f-ODE.
30. `nb.derivation("D29", …)` — Part F D29 (10 steps), ref "8.48".
31. `nb.note` — **N87 [B]** "**Velocities** (one derivative each through (6.83)): $u_r=U\\cos\\theta\\big(1-\\frac{3a}{2r}+\\frac{a^3}
    {2r^3}\\big)$, $u_\\theta=-U\\sin\\theta\\big(1-\\frac{3a}{4r}-\\frac{a^3}{4r^3}\\big)$ (8.49)." equation (8.49) + `nb.code`
    (sympy, every line commented): `print(sp.simplify(S["u_r"]), sp.simplify(S["u_theta"]))` · `print(S["u_r"].subs(r, a),
    S["u_theta"].subs(r, a))  # 0 0: no slip`.
32. `nb.worked_example("a bead in syrup, two radii out on the side", "a = 1 mm, U = 1 mm/s, ν = 10⁻³ m²/s (syrup: Re =
    2aU/ν = 0.002). At r = 2a, θ = π/2: u_r = 0, $u_\\theta=-U(1-\\frac38-\\frac1{32})=-0.594U$. Ideal flow (Ch. 6) would give
    −U(1 + 1/16) = −1.063U: the viscous fluid is *slowed* beside the sphere instead of sped up. At r = 10a the Stokes speed is
    0.925U — still a 7.5 % deficit; the ideal-flow disturbance is already 0.05 %.")`
33. `nb.code` — `for r in (2.0, 10.0): print(r, ch08.side_line_speed(r*1e-3, U=1e-3, a=1e-3, model="stokes"),
    ch08.side_line_speed(r*1e-3, U=1e-3, a=1e-3, model="ideal"))` · `print(ch08.stokes_sphere_velocity(2e-3, np.pi/2,
    U=1e-3, a=1e-3))` · `print(ch08.stokes_sphere_streamfunction(2e-3, np.pi/4, U=1e-3, a=1e-3, frame="fluid"))`.
    *expect:* 0.594e-3 vs 1.0625e-3 m/s at 2a; 0.925e-3 vs 1.0005e-3 at 10a; (0, −5.94e-4); the fluid-frame ψ value printed.
34. `nb.check_agree` — **from scratch (curation §7):** central differences of ψ through (6.83) at 20 random points outside the
    sphere (h = 1e-7 m), `assert np.allclose((ur_fd, uth_fd), ch08.stokes_sphere_velocity(rr, tt, U=1e-3, a=1e-3),
    rtol=1e-6)`.
35. `nb.note` — **N93 [B]** "**The sphere moving through still fluid** (fluid frame): subtract the stream, $\\psi=Ur^2\\sin^2\\theta
    \\big(-\\frac{3a}{4r}+\\frac{a^3}{4r^3}\\big)$ (Fig. 8.19). The pattern is fore–aft symmetric — no wake — because (8.43) is
    linear: reversing U maps u → −u and p − p∞ → −(p − p∞). ⚠️ The text calls the equation '(9.63)'; it means (8.43). Ch. 16
    returns to this reversibility (why a swimming bacterium cannot use a reciprocal stroke)."
36. `nb.figure` — **Figs. 8.17 (sphere) and 8.19 remade (N107 [B], the §8.6 figures; Fig. 8.16 was remade in C11)**
    (10 × 4 in, three panels): (a) body-frame streamlines of (8.48) with the speed as a heatmap; (b) fluid-frame streamlines
    (closed loops, symmetric front/back); (c) the ideal-flow sphere of Ch. 6 (`PF.sphere`) in the fluid frame for contrast.
    Title "Creeping flow: symmetric, and felt far away". *see:* "(b) loops that look the same in front and behind; (c) a much
    more compact disturbance." *read:* "count the radii where the loops still bend: in (b) the disturbance is visible beyond
    10a (decays like a/r), in (c) it is gone by 3a (a³/r³)." *change:* "…Re were 1: inertia breaks the symmetry and a wake
    forms behind (C15)."
37. `nb.animation` — **A4** (video, 60 frames FAST 30): tracers released on a vertical line ahead of a sphere moving left
    (fluid frame), next to the same for ideal flow; Stokes tracers are pushed aside and return almost to where they started,
    ideal-flow tracers barely move. Notes *see / read / change*.
38. `nb.explainer("stokes_sphere_flow", heading="What does creeping flow round a sphere look like?", why="Toggle body and fluid
    frames and Stokes, Oseen and ideal flow, and drag Re: symmetric loops, a disturbance that decays only like a/r, and then
    a wake that grows from far downstream once r reaches a/Re — three comparisons that need motion and switching to be
    believed. (Its Oseen mode previews C15.)", tries=["Fluid frame, Stokes: follow one tracer — it goes round a loop and
    comes back.", "Read the side-line speed view: at 10a Stokes still shows a 7.5 % deficit, ideal flow none.", "Switch to
    Oseen and raise Re to 1: where does the circle r = a/Re_a sit, and where does the wake begin?", "Open the Derivation tab
    at D29 step 8: zoom out to see why the r⁴ term must go."])`
39. `nb.md` — **What would change if…** "…we asked for the force? The pressure and the shear stress on the surface follow from
    (8.49); integrating them gives the drag (C14)."
40. `nb.recap("R20", "Dimensional analysis without ρ", "If inertia does not matter, ρ cannot appear: D = f(μ, U, a); four
    variables and three dimensions leave one group, D/(μUa) = constant, so D ∝ μUa and C_D ∝ 1/Re before any flow
    calculation (Ch. 1's Π method; Exercise 4.60).", where="Ch. 1 §1.11")` + `nb.code`: `print(ch01.pi_groups({"D": "N",
    "mu": "Pa*s", "U": "m/s", "a": "m"}))` *expect:* one group ∝ D/(μUa).
41. `nb.core("C14", "Stokes drag $D=6\\pi\\mu aU$ (8.51)", question="What force does the fluid exert on the sphere — and where
    does it come from, pressure or friction?")`
42. `nb.md` — **The problem in plain words:** "A cloud droplet falls at its terminal speed when the drag equals its weight minus
    buoyancy; Millikan used exactly this to weigh single oil drops and count electrons; sediment settles in rivers and the
    ocean; aerosols stay aloft for days. All need the drag of a slowly moving sphere — and the surprise is how it splits."
43. `nb.md` — **The idea** (ASCII): "the traction on the surface has a pressure part and a friction part; their x-components
    add up to a *uniform* push 3μU/2a over the whole sphere:"
    ```
    front (θ = π): high pressure +3μU/2a  ─►  ● ─►  rear (θ = 0): low pressure −3μU/2a
    sides: friction σ_rθ = −(3μU/2a) sin θ drags the surface along
    drag = ⅓ pressure (2πμaU) + ⅔ friction (4πμaU) = 6πμaU
    ```
44. `nb.md` — gloss **line integral of a gradient** (P35 reminder) and **traction projection** (Ch. 2 §2.6 Cauchy; P163 surface
    integrals on a sphere): "knowing ∂p/∂r and (1/r)∂p/∂θ, integrate one of them and check the other; the force per area on
    the sphere is t = σ·e_r, and its x-component is σ_rr cos θ − σ_rθ sin θ with θ measured from +x; dA = 2πa² sin θ dθ."
45. `nb.derivation("D30", …)` — Part F D30 (11 steps, ★★★), ref "8.50". **check_src:** `print(S["p"])  # −3μaU cos θ/(2r²)` ·
    `assert sp.simplify(S["p_check_r"]) == 0 and sp.simplify(S["p_check_theta"]) == 0  # both components of ∇p = μ∇²u` ·
    `print(S["p"].subs({r: a, th: 0}), S["p"].subs({r: a, th: sp.pi}))  # −3μU/2a at the rear, +3μU/2a at the front` —
    every line commented.
46. `nb.note` — **N88 [B]** "**The pressure** $p-p_\\infty=-\\frac{3\\mu aU\\cos\\theta}{2r^2}$ (8.50): maximum +3μU/2a at the front
    stagnation point (θ = π), minimum **−**3μU/2a at the rear (θ = 0). ⚠️ The book prints the minimum without the minus
    sign (and Fig. 8.17 labels both extremes 1.5)." **N89 [B]** "**Surface stresses** (ours): on r = a the viscous normal stress
    vanishes, σ_rr = −p, and $\\sigma_{r\\theta}=-\\frac{3\\mu U}{2a}\\sin\\theta$." equation (8.50).
47. `nb.derivation("D31", …)` — Part F D31 (13 steps, ★★★), ref "8.51". **check_src:** `print(S["sigma_rr_viscous_a"],
    S["sigma_rtheta_a"])  # 0 and −(3μU/2a) sin θ` · `print(sp.simplify(S["D_pressure"]), sp.simplify(S["D_friction"]))  #
    2πμaU and 4πμaU` · `assert sp.simplify(S["D_pressure"] + S["D_friction"] - 6*sp.pi*mu*a*U) == 0  # (8.51)` — every line
    commented.
48. `nb.worked_example("the drag and fall speed of a 10 µm cloud droplet", "a = 10 µm, water droplet ρ′ = 1000 kg/m³ in air
    ρ = 1.2 kg/m³, μ = 1.81 × 10⁻⁵ Pa s, g = 9.81 m/s². 1. Terminal balance $\\frac43\\pi a^3g(\\rho'-\\rho)=6\\pi\\mu aU$ ⇒
    $U_t=\\frac{2(\\rho'-\\rho)ga^2}{9\\mu}=\\frac{2\\times998.8\\times9.81\\times10^{-10}}{9\\times1.81\\times10^{-5}}=1.20$ cm/s.
    2. Re = 2aU_tρ/μ = 2 × 10⁻⁵ × 0.012 × 1.2/1.81 × 10⁻⁵ = 0.016 ✓ ≪ 1. 3. D = 6πμaU_t = 4.1 × 10⁻¹¹ N: 1.4 × 10⁻¹¹ N from
    pressure, 2.7 × 10⁻¹¹ N from friction. 4. C_D = 24/Re ≈ 1.5 × 10³. 5. A 100 µm drizzle drop: 100× faster by a² —
    1.2 m/s — but then Re ≈ 16 and Stokes no longer applies.")`
49. `nb.code` — `st = ch08.settling_state(10e-6, 1000.0, 1.2, 1.81e-5)` · `print(st)` · `print(ch08.stokes_drag(1.81e-5, 10e-6,
    st["U_t"], parts=True))` · `print(ch08.stokes_sphere_surface_stresses(np.array([0, np.pi/2, np.pi]), U=st["U_t"],
    a=10e-6, mu=1.81e-5))` · `print(ch08.stokes_drag_running(np.pi/2, 1.81e-5, 10e-6, st["U_t"]))` ·
    `print(ch08.stokes_drag_coefficient(st["Re"]), SIM.sphere_drag_coefficient(st["Re"], "stokes"))`. *expect:* U_t =
    1.203e-2 m/s (g = G0), Re = 0.0159, D = 4.10e-11 N, C_D = 1505, valid True, D_pressure 1.37e-11, D_friction 2.74e-11;
    tractions t_x = 3μU/2a = 0.0326 Pa at every θ; running drag to θ = π/2: pressure 6.84e-12, friction 1.37e-11 N;
    the two C_D equal. *explain:* 1. `settling_state` solves the terminal balance and checks Re < 0.1; 2. the surface
    tractions' x-component is the same everywhere — the pressure and friction parts trade off; 3. `stokes_drag_running`
    shows how the total is collected from the rear to the front.
50. `nb.check_agree` — **from scratch (curation §7):** midpoint sum over 2000 θ-cells of the pressure and friction
    x-tractions × 2πa² sin θ Δθ; `assert np.allclose([Dp, Df], [D["pressure"], D["friction"]], rtol=1e-6)`; also
    `ch08.sphere_drag_quadrature(...)` with n = 16 agrees to 1e-12 (spectral convergence).
51. `nb.note` — **N90 [B]** "**Terminal velocity:** $(4/3)\\pi a^3g(\\rho'-\\rho)=6\\pi\\mu aU$ ⇒ $U_t=\\frac{2(\\rho'-\\rho)ga^2}{9\\mu}$
    — ∝ a²: halve the droplet, quarter the speed. This sets how long aerosols and cloud droplets stay aloft and how fast silt
    settles (Ch. 13)." **N91 [B]** "**Millikan's oil drops (Fig. 8.18).** With the plates uncharged the fall speed gives the
    radius; with the upper plate negative the drop rises at U_u: $6\\pi\\mu U_ua+(4/3)\\pi a^3g(\\rho'-\\rho)=neE$, E = −V_b/L;
    charges from many drops are whole multiples of e." + `nb.code`: `m = ch08.synthetic_millikan(40, seed=0)` ·
    `print(m["e_est"], m["e_est"]/1.602176634e-19 - 1)` *expect:* e within 1 % of the CODATA value (benchmark V5).
52. `nb.note` — **N92 [B]** "**Drag coefficient:** with A = πa² (4.107), $C_D=\\frac{D}{\\frac12\\rho U^2\\pi a^2}=\\frac{24}{\\mathrm{Re}}$
    (8.52), Re = 2aU/ν based on the **diameter** — the low-Re branch of the drag curve of Ch. 4." equation (8.52).
53. `nb.figure` — **Fig. 8.17 remade (N107 part)** (10 × 3.4 in, three panels): (a) the sphere with traction arrows around it
    (pressure orange, friction rose, their sum black — all the same length in x); (b) p − p∞ and σ_rθ (in units of μU/a) vs θ
    from 0 (rear) to π (front), with the running drag P(θ), F(θ) of `stokes_drag_running` on a second axis reaching 2π and 4π
    (in μaU); (c) U_t vs radius on log–log (1 µm … 1 mm) for water droplets in air, quartz sand in water, bacteria in water,
    the Re = 0.1 limit marked on each line (Stokes invalid beyond, dashed). Title "A third pressure, two thirds friction;
    fall speed ∝ a²". *see / read / change* ("…the droplet were oil (ρ′ = 900): U_t falls by 10 %").
54. `nb.explainer("stokes_drag_settling", heading="Where does 6πμaU come from?", why="Sweep around the sphere and watch the local
    pressure and shear tractions as arrows and curves while their running integrals build up to 2πμaU and 4πμaU in the term
    bars; then drag the particle radius and read U_t and Re as the Stokes law's validity limit is crossed.", tries=["Play the
    sweep from the rear (θ = 0) to the front: which bar grows first?", "Pick 'drizzle 100 µm': the status turns amber — why?",
    "Switch the third view to C_D(Re) and find where Oseen and the correlation leave 24/Re.", "Click a point on the sphere
    and read the traction arithmetic."])`
55. `nb.md` — **What would change if…** "…we looked far from the sphere? The viscous force there decays like a/r³ but inertia
    only like a/r²; at some distance inertia wins even for tiny Re (C15)."
56. `nb.core("C15", "Where Stokes fails, and Oseen's fix: inertia/viscous $\\sim\\mathrm{Re}\\,\\frac ra$", question="If Re is tiny,
    how can inertia ever matter?")`
57. `nb.md` — **The problem in plain words:** "Far behind a slowly falling bead the fluid is barely disturbed — so inertia should
    be even less important there. It is the opposite: the viscous force falls off faster with distance than inertia does,
    so beyond some radius Stokes' solution is wrong. For a cylinder it is so wrong that no Stokes solution exists at all."
58. `nb.md` — **The idea** (table):
    | at distance r | size |
    |---|---|
    | disturbance velocity | ~ Ua/r |
    | inertia ρu·∇u | ~ ρU²a/r² |
    | viscous μ∇²u | ~ μUa/r³ |
    | ratio | ~ (ρUa/μ)(r/a) = Re_a r/a → equal at r ~ a/Re_a |
59. `nb.md` — gloss **asymptotic size of a term** (P130 reminder) and **series 1 − e^{−s} = s − s²/2 + …** (P26, P117).
60. `nb.derivation("D32", …)` — Part F D32 (7 steps), ref "".
61. `nb.note` — **N94 [C]** "**A singular perturbation at infinity.** Treated as the first term of an expansion in Re, Stokes'
    solution is not uniformly valid: the O(Re) correction grows without bound relative to it as r → ∞ (Whitehead's
    paradox); for a cylinder the Stokes equations cannot even meet the uniform stream (Stokes' paradox, Exercise 8.37).
    Unlike the boundary layer (Ch. 9), where 1/Re multiplies the highest derivative near the wall, here the trouble is far
    away."
62. `nb.note` — **N95 [B]** "**Oseen's linearisation (1910).** Write u = U + u′ with u′ small far away. The advective term is
    $u\\frac{\\partial u}{\\partial x}+v\\frac{\\partial u}{\\partial y}+w\\frac{\\partial u}{\\partial z}=U\\frac{\\partial u'}{\\partial x}+
    \\big[u'\\frac{\\partial u'}{\\partial x}+v'\\frac{\\partial u'}{\\partial y}+w'\\frac{\\partial u'}{\\partial z}\\big]$; dropping the
    bracket gives $\\rho U\\frac{\\partial u_i'}{\\partial x}=-\\frac{\\partial p}{\\partial x_i}+\\mu\\nabla^2u_i'$. ⚠️ The book prints
    +∂p/∂x_i; the minus of (8.43) is needed." **N96 [B]** "Conditions: u′, v′, w′ → 0 far away; u′ = −U, v′ = w′ = 0 on the
    sphere." + `nb.code`: `print(ch08.oseen_linearisation_sympy()["oseen_x"])`.
63. `nb.note` — **N97 [B]** "**Oseen's stream function** $\\frac{\\psi}{Ua^2}=\\big[\\frac{r^2}{2a^2}+\\frac{a}{4r}\\big]\\sin^2\\theta-
    \\frac{3}{\\mathrm{Re}}(1+\\cos\\theta)\\big\\{1-\\exp\\big[-\\frac{\\mathrm{Re}}{4}\\frac ra(1-\\cos\\theta)\\big]\\big\\}$ (8.53),
    Re = 2aU/ν (stated; no slip holds only to O(Re))." **N98 [B]** "Near the sphere the exponential's series returns (8.48).
    ⚠️ The text calls it '(9.68)'; it means (8.48)." equation (8.53).
64. `nb.derivation("D33", …)` — Part F D33 (7 steps), ref "8.53". + `nb.code`: `ol = ch08.oseen_limit_sympy();
    print(sp.simplify(ol["difference"]))` *expect:* 0.
65. `nb.worked_example("where inertia catches up", "Re = 2aU/ν = 0.02, so Re_a = 0.01. 1. Ratio ~ Re_a r/a = 1 at r = 100a. 2.
    For the 10 µm droplet (Re = 0.016) at r = 125a = 1.25 mm. 3. Oseen C_D at Re = 0.5: 48(1 + 3/32) = 52.5 vs Stokes 48 — a
    9 % correction.")`
66. `nb.code` — `for r in (10.0, 100.0, 1000.0): print(r, ch08.inertia_viscous_ratio(r*1e-3, np.pi/2, U=1e-3, a=1e-3,
    nu=1e-4))` (Re_a = Ua/ν = 0.01) · `print(ch08.oseen_drag_coefficient(0.5), ch08.stokes_drag_coefficient(0.5),
    ch08.proudman_pearson_drag_coefficient(0.5), SIM.sphere_drag_coefficient(0.5))` · `print(ch08.oseen_streamfunction(2e-3,
    np.pi/3, U=1e-3, a=1e-3, Re=1e-6), ch08.stokes_sphere_streamfunction(2e-3, np.pi/3, U=1e-3, a=1e-3))`. *expect:* ratios
    growing ∝ r with order Re_a r/a (≈ 0.1, 1, 10 up to an O(1) factor printed by the function); C_D 52.5, 48.0, the PP value
    and the Morrison correlation between them; the two ψ equal to ~1e-6 relative.
67. `nb.check_agree` — **from scratch (curation §7):** finite-difference ∣u·∇u∣/∣ν∇²u∣ at three points on the side line from
    `stokes_sphere_velocity_xyz`, `assert np.allclose(mine, ch08.inertia_viscous_ratio(...), rtol=1e-4)`; the log–log slope
    of the ratio for r/a ∈ [50, 500] is 1 ± 0.05 (`observed_order`).
68. `nb.figure` — (10 × 3.4 in, three panels): (a) inertia/viscous ratio vs r/a on log–log for Re_a = 0.001, 0.01, 0.1 (slope
    1), the line ratio = 1 marked; (b) Stokes vs Oseen streamlines in the fluid frame at Re = 1 (**N100 [B]**, Fig. 8.20 remade:
    the Oseen pattern is asymmetric with a wake behind); (c) C_D·Re/24 vs Re (0.01 … 10, log x): Stokes (1), Oseen (1 +
    3Re/16), Proudman–Pearson (**N101 [C]**, named: $D=6\\pi\\mu aU(1+\\frac38\\mathrm{Re}_a+\\frac9{40}\\mathrm{Re}_a^2\\ln\\mathrm{Re}_a+
    \\ldots)$, from matched asymptotic expansions — Kaplun, Proudman & Pearson 1957), and the Morrison correlation (the
    'experiments', benchmark band) lying between Stokes and Oseen for Re < 5 (**N99 [B]**: $C_D=\\frac{24}{\\mathrm{Re}}\\big(1+
    \\frac3{16}\\mathrm{Re}\\big)$; ⚠️ with the radius-based Re_a the same coefficient reads 3/8). Title "Stokes is right near
    the sphere and wrong far away". *see / read / change* ("…Re doubled: the crossover distance halves and the wake widens").
69. `nb.plotly` — **IF8** `slider_figure(…, "Re", [0.1, 0.2, 0.5, 1.0, 2.0], …)` — precomputed Oseen streamline polylines
    (contours of `oseen_streamfunction`, 200² grid, 100² FAST) and the C_D values in the title.
70. `nb.md` — "🎮 The explainer of C13 (`stokes_sphere_flow`) has an Oseen mode: open it again, choose Oseen and drag Re to see
    the wake grow from far downstream." (a pointer, not a second embed).

### A.7 §8.7 Final Remarks — N102, S01, summary
1. `nb.section("8.7", "Final Remarks", intro="**What is this section about?** Where the exact solutions stop.")`
2. `nb.note` — **N102 [C]** "Most laminar problems that can be solved with pencil and paper have been solved; the field moves on
   with perturbation methods (flows close to a known one — Ch. 9 boundary layers, Ch. 11 stability) and with numerical
   solution of the Navier–Stokes equations (Ch. 10, which uses this chapter's exact solutions as test cases)."
3. `nb.pointer` — **S01** "Exercises 8.1–8.38 and the literature are not reproduced. Ideas from them used here: the Reynolds
   equation (Exercises 8.19–8.20 → C06), the power–dissipation identity (8.12 → R11), the switched-on vortex (8.26 → C10),
   the stopped plate (8.30 → C09), Couette start-up (8.31 → R14), Hele-Shaw near an obstacle (8.34 → C06), Stokes drag from
   the stresses (8.35 → C14), Stokes' paradox for a cylinder (8.37 → C15), f = 64/Re (8.7 → C03)."
4. `nb.summary(clicked=[…15…], feeds_forward=[…], left_out=[…])` — **clicked** (one line per CORE): C01 "ν = μ/ρ is a
   diffusivity: motion spreads √(νt), and a flow stays laminar below Re ≈ 2000." · C02 "Between plates the profile is a
   Couette line plus a Poiseuille parabola; the floor flow reverses when dp/dx > 2μU/h²." · C03 "In a pipe Q = −(πa⁴/8μ)dp/dz:
   halve the radius, lose 94 % of the flow; f = 64/Re." · C04 "Between rotating cylinders the swirl is AR + B/R: rigid rotation
   plus a free vortex, chosen by the walls." · C05 "In a thin gap inertia is weighed with ε²Re_L, not Re_L; pressure is
   uniform across the gap." · C06 "Station by station the gap flow is Couette plus Poiseuille; integrating continuity gives
   the Reynolds equation h_t + q_x = 0." · C07 "A constant flux through a narrowing gap needs a pressure hump — that hump
   carries the load, but only if the pad slides the right way." · C08 "A spreading layer obeys h_t = (ρg/3μ)(h³h_x)_x, a
   diffusion that chokes itself as h shrinks." · C09 "With no imposed scale, y and t combine as y/√(νt): one erfc curve for
   every time, δ₉₉ = 3.64√(νt)." · C10 "Guess At^{−n}F(ξ/δ(t)); matching powers of t and one conserved quantity fix n and δ."
   · C11 "An oscillating wall drives a decaying, lagging layer of depth √(2ν/ω) — diffusion, not a wave." · C12 "At low Re
   rescale pressure by μU/L; the Stokes equations ∇p = μ∇²u are linear." · C13 "Stokes' sphere flow is fore–aft symmetric
   and decays only like a/r." · C14 "Drag = 6πμaU, ⅓ pressure and ⅔ friction; fall speed ∝ a²." · C15 "Inertia returns at
   r ~ a/Re_a; Oseen's linearisation adds the wake and C_D = (24/Re)(1 + 3Re/16)." **feeds_forward:** Ch. 9 (√(νt) ↔
   √(νx/U), Blasius vs the temporal layer, the lubrication scaling as the boundary-layer approximation, similarity exponents),
   Ch. 10 (exact solutions as code tests, Crank–Nicolson), Ch. 11 (Couette, Poiseuille, Taylor–Couette base states, the
   Rayleigh criterion in A and B), Ch. 12 (linear total stress, τ₀ and the friction velocity, f = 64/Re on the Moody chart),
   Ch. 13 (Ekman layers with the (1 + i)/δ structure of (8.37), spin-up times from √(νt), cyclostrophic balance, viscous
   gravity currents, settling of droplets and sediment), Ch. 16 (arteries, synovial lubrication, micro-swimmers).
   **left_out:** the exercises (S01); higher-order matched expansions for the sphere (named in N101).

### A.8 Placement check (every curation id has exactly one home)
| IDs | Block | Section |
|---|---|---|
| R01, R02 | before C01 | §8.1 |
| C01, N01, N02, N103, D01 | C01 | §8.1 |
| R03, R04, R05, R06 | after C01 | §8.1 |
| R07, R08 | before C02 | §8.2 |
| C02, N03–N11, D02–D05, IF1, E1 | C02 | §8.2 |
| R09 | before C03 | §8.2 |
| C03, N12–N16, D06, D07, IF2 | C03 | §8.2 |
| R10, R11, R12 | before C04 | §8.2 |
| C04, N17–N22, N104, D08, D09, IF3 | C04 | §8.2 |
| R13 | before C05 | §8.3 |
| C05, N23, N24, N26–N31, D10, D11, E2 | C05 | §8.3 |
| C06, N25, N32, N33, N34, N36, D12, D13 | C06 | §8.3 |
| C07, N35, N105, D14, IF4, live, E3 | C07 | §8.3 |
| C08, N37, D15, A3, E4 | C08 | §8.3 |
| R14 | before C09 | §8.4 |
| C09, N38–N56, D16–D20, A1, IF5, E5 | C09 | §8.4 |
| C10, N57–N64, N106, D21–D23, A5, E6 | C10 | §8.4 |
| C11, N65–N74, D24, D25, A2, IF6, E7 | C11 | §8.5 |
| R15, R16 | before C12 | §8.6 |
| C12, N75–N80, D26, IF7 | C12 | §8.6 |
| R17, R18, R19 | before C13 | §8.6 |
| C13, N81–N87, N93, N107, D27–D29, A4, E8 | C13 | §8.6 |
| R20 | before C14 | §8.6 |
| C14, N88–N92, D30, D31, E9 | C14 | §8.6 |
| C15, N94–N101, D32, D33, IF8 | C15 | §8.6 |
| N102, S01 | §8.7 section (after C15) | §8.7 |
Totals: 15 CORE blocks, 107 NOTE ids, 20 RECAP calls, 1 SKIP pointer, 33 derivations (296 steps), 9 explainers (backup B1
not embedded), 5 animations, 8 plotly figures, 1 live cell, 15 primers (P185–P199).

---

## Part B — explainer storyboards

Common to all ten: created with `tools/new_viz.py`; `<meta name="viz:chapter" content="ch08">`; tabs Walkthrough ·
Explore · Explain · Derivation · Equations · Code · Check; every displayed number is computed by a JS function that
mirrors a `ch08` callable and is proved by `selftest()` parity rows (`py:` expressions use only `ch08.…`, `np.pi`,
numbers, strings, lists, dicts and keywords, indexed down to one float — Part C convention 3). Explain is "Explanation &
interpretation" in numbered sections built with `Viz.work.step / line / box / table / hint / interpret`, modelled on
`forced_damped_vibrations.html` and `amplitude_phase_second_order_II_3.html`: **0** what the views show and what each colour
means (two sentences on phones) · **1…n** every displayed quantity from the controls ("formula = substituted = result —
why", results boxed) · a section per view hidden on phones, or a hint · the values at the current time (live) · **Reading
the current setting** (regime-dependent). Derivation steps are copied from Part F (same `did` titles, same step count;
phones shorten *why* to its first sentence; plain-text *why* and *watch* never contain raw TeX — write e^(−y/δ), not
`e^{…}`). Every tour, Explain, Derivation, notes, status, equation and quiz text that names a book equation **writes it
out** next to its number (`tools/eq_refs.py` → 0; no bare numbers in `<meta>` strings). Drafts below write equations in
Unicode for readability; builders set each in TeX (backslashes doubled in JS strings). Colours as in convention 6
(pressure/Poiseuille orange, viscous/Couette rose, inertia teal, velocity blue, vorticity purple, backflow amber, rescaled
purple dashed, analytic ghost muted dashed, numerical dots teal, printed slip muted dotted). Walkthrough texts ≤ 45 words,
step 1 ≤ 24 words, ≤ 2 extras per step, `play: false` on steps that quote numbers. A view hidden on portrait phones never
carries a step's key number (repeat it in a visible title or readout). **erf/erfc:** E5 and E6 define a local
double-precision `erfcHP` (continued fraction for x > 3, series below; ch04 E4 pattern) — `Viz.num.erfc` (1.2e-7) is not
used for displayed or parity numbers. Lengths and speeds print with a unit switch (µm/mm/cm/m) at 3–4 s.f.; numbers that
coincide with printed book values (3.64, 0.06, 5.54, 2.76) are shown at 4 s.f. from the function (3.643, 0.05911, 5.544,
2.772). Each explainer fits 360×640 … 1920×1080 and the 1000×700 notebook frame with no scrolling (fit plan per
explainer). Parallel builders use private scratch subfolders (`<scratchpad>/<slug>/`).

### E1 · couette_poiseuille_backflow
- **Title:** "When does fluid flow backwards in a channel?" · **Summary:** "A moving wall drags, a pressure gradient pushes:
  the profile is a straight line plus a parabola, and the floor layer reverses once dp/dx exceeds 2μU/h²." · **CORE:** C02
  (also N03–N11, R07, R08) · **Reference:** `forced_damped_vibrations.html` (system + graph linked on one state, live
  explanation computing every number) with the term bars of `fid_formula_lab.html`.
- **meta:** `viz:order 1` · `viz:sections 8.2` · `viz:equations 8.4 8.5` · `viz:fluidpy ch08.channel_flow
  ch08.channel_flow_rate ch08.channel_shear_stress ch08.channel_backflow_threshold ch08.couette_poiseuille_state` ·
  `viz:derivations D02 D03 D04 D05`.
- **Physics (JS ↔ Python):** `u(y, s)` = (U/h)y − (1/2μ)(dp/dx)y(h − y) ↔ `ch08.channel_flow(y, h, U, dpdx, mu)`; `flow(s)`
  → {Q, V} ↔ `ch08.channel_flow_rate`; `tau(y, s)` = μU/h − (h/2 − y)dp/dx ↔ `ch08.channel_shear_stress`; `thr(s)` = 2μU/h²
  ↔ `ch08.channel_backflow_threshold`; `state(s)` (Q_couette, Q_poiseuille, tau_bottom, tau_top, backflow, ratio, y_reversal,
  u_max, y_umax, zero_flow_dpdx) ↔ `ch08.couette_poiseuille_state`. μ = 10⁻³ Pa s (water) unless the optional slider moves.
- **Views** (rows [1.1, 1]): 1. `channel` "The channel" (row 0, flex 1.4) — x ∈ [0, 6h], y ∈ [0, h] (drawn in mm), hatched
  fixed floor, moving top wall with an arrow U (rose), 40 tracer dots advected at u(y) (blue, wrap around), a dye line that
  starts vertical at t = 0 and bends into the profile (x = x₀ + u(y)t, wrapped), the reversed layer 0 < y < y_rev shaded
  amber; the pressure along the channel as a thin orange bar at the top (dark at high p) with "p falls/rises downstream"
  label. Animates with the transport. 2. `profile` "u(y)" (row 0, flex 1) — u [cm/s] (x) vs y/h (y): Couette part rose
  dashed, Poiseuille part orange dashed, sum blue bold, the u = 0 line, the floor slope as a short black tangent at y = 0,
  y_rev marked amber, u_max dot; pointer: click → probe height. 3. `stress` "τ(y) and Q" (row 1, `hidePortrait: true`) —
  left half: τ(y) [mPa] line (black) with the two wall values labelled; right half: stacked bar Q = Q_couette (rose) +
  Q_poiseuille (orange, negative below zero) with the total (blue tick). Portrait: `channel` + `profile`; Q and τ(0) repeated
  in the `profile` title ("Q = 1.67 cm²/s · τ₀ = −10 mPa").
- **Controls (≤ 5 visible):** `dpdx` "Pressure gradient $dp/dx$" −30 … 30 Pa/m, step 0.1, default −4, help "negative =
  favourable (pushes with the wall)" · `U` "Wall speed $U$" −0.2 … 0.2 m/s, step 0.005, default 0.1 · `h` "Gap $h$"
  0.5 … 2 cm, step 0.05, default 1 cm · `mu` "Viscosity $\mu$" (optional) 0.5 … 5 mPa s, default 1 · `probe` (hidden, set by
  clicks) default h/2.
- **Transport:** `t` 0 → 20 s, rate 2 s/s, `end: 'loop'` (dye and tracers only; the steady profile does not change).
- **Presets:** "pure Couette" {dpdx: 0} · "pure Poiseuille" {U: 0, dpdx: −4} · "favourable" {dpdx: −4} · "onset of backflow"
  {dpdx: 2} (= 2μU/h² at the defaults; the preset computes it from the current U, h, μ) · "strong adverse" {dpdx: 8} ·
  "zero net flow" {dpdx: 6} (= 6μU/h², computed).
- **Status:** "➡️ all forward — floor slope U/h − (h/2μ)dp/dx = … 1/s > 0" · "⚠️ backflow near the fixed wall: dp/dx = … ×
  2μU/h², reversed layer 0 < y < … mm" · "⏸ zero net flow: Q = 0 at dp/dx = 6μU/h²" (|Q| < 1 % of Uh/2) · "pure Poiseuille
  (U = 0)".
- **Readouts:** "Flow rate $Q$" (cm²/s) · "Mean speed $V$" (cm/s) · "Floor stress $\tau(0)$" (mPa) · "Threshold $2\mu U/h^2$"
  (Pa/m) · "dp/dx ÷ threshold".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** · **status** · **terms**
  (Q = Couette part + Poiseuille part; click a bar to isolate that part in the profile) · **inspector** (click a height:
  "u = (U/h)y − (1/2μ)(dp/dx)y(h − y) = 0.1 × 0.25 − (4/0.002) × 0.0025 × 0.0075 = 0.025 − 0.0375 = **−0.0125** m/s").
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**The channel**: blue dots move at the local speed, the dye line bends into the profile, amber
     marks fluid moving backwards. **u(y)**: rose dashed = the wall-driven (Couette) part, orange dashed = the
     pressure-driven (Poiseuille) part, blue = their sum (8.5). **τ(y) and Q** (hidden on phones): the stress line and the
     flow rate split into its two parts."
  1. *The two parts at your probe height* — "Couette: U y/h = **live**; Poiseuille: −(1/2μ)(dp/dx)y(h − y) = **live**; sum
     u = **live** m/s" boxed; why: (8.5) $u=\frac Uhy-\frac1{2\mu}\frac{dp}{dx}y(h-y)$ is linear, so causes add.
  2. *Wall stresses* — "τ = μ du/dy = μU/h − (h/2 − y)dp/dx: floor **τ(0)**, top **τ(h)**" boxed; "for pure Poiseuille both
     are (h/2)∣dp/dx∣ in size".
  3. *Flow rate and mean speed* — "Q = Uh/2 − h³(dp/dx)/(12μ) = **Qc** + **Qp** = **Q** m²/s; V = Q/h = **V**" boxed; hint
     "the book writes V ≡ Q/h with the middle integral missing its 1/h — units tell (m²/s ≠ m/s)".
  4. *Backflow threshold* — "the floor fluid reverses when du/dy(0) = U/h − (h/2μ)dp/dx < 0, i.e. dp/dx > 2μU/h² =
     **thr** Pa/m; you are at **ratio** × threshold; reversed layer up to y = h − 2μU/(h dp/dx) = **yrev**" boxed.
  5. *At the current time* — "t = **live t** s: the dye at mid-height has moved u(h/2)·t = **live** cm; a floor particle at
     y = 0.1h moved **live** cm (negative = backwards)."
  6. *Reading the current setting* — favourable: "The pressure pushes with the wall: the profile is fuller than the straight
     line and every layer moves forward. Q exceeds the pure-drag value Uh/2." · adverse, below threshold: "The pressure pushes
     back but the wall's pull still wins at the floor: the profile sags but stays forward." · backflow: "Adverse gradient
     beyond 2μU/h²: next to the fixed floor the pressure wins and the fluid flows backwards — the same competition makes
     boundary layers separate (Ch. 9)." · zero net flow: "Forward drag near the wall and backflow near the floor cancel:
     a closed recirculation, as in a sealed extrusion die or a wind-driven lake with a return current."
- **Derivation tab:** **D02** (7 steps) `view: 'channel'`; goal `set` {dpdx: −4}; step 3 (v ≡ 0) `watch` "the tracers move
  only horizontally". **D03** (5 steps) `view: 'channel'`; step 5 `watch` "the orange pressure bar changes linearly along x".
  **D04** (9 steps) `view: 'profile'`; step 5 (apply u(h) = U) highlights the top wall; step 9 **live** "u(h/2) = U/2 −
  (h²/8μ)dp/dx = …". **D05** (7 steps) `view: 'profile'`; step 7 `set` preset "onset of backflow", **live** "2μU/h² = 2 ×
  0.001 × 0.1/0.0001 = 2 Pa/m", `watch` "the black floor tangent is vertical". Interpret: `s => "With your numbers the
  floor reverses above dp/dx = ${thr} Pa/m; you are at ${ratio}× that."`.
- **Code:**
  ```python
  h, U, dpdx, mu = {{h}}, {{U}}, {{dpdx}}, {{mu}}        # gap [m], wall speed [m/s], dp/dx [Pa/m], μ [Pa s]
  y = np.linspace(0, h, 101)                        # heights [m]
  u = ch08.channel_flow(y, h, U=U, dpdx=dpdx, mu=mu)  # (8.5): line + parabola; u(h/2) = {{umid}} m/s
  Q, V = ch08.channel_flow_rate(h, U=U, dpdx=dpdx, mu=mu)   # Q = {{Q}} m²/s, V = {{V}} m/s
  tau0 = ch08.channel_shear_stress(0.0, h, U=U, dpdx=dpdx, mu=mu)   # floor stress = {{tau0}} Pa
  thr = ch08.channel_backflow_threshold(U, h, mu)   # 2μU/h² = {{thr}} Pa/m
  print(dpdx > thr)                                 # backflow? {{back}}
  ```
- **Walkthrough (6 steps):** 1. "Belt versus pump" — "A belt drags fluid forward; a pump pushes back. Who wins near the
  floor?" `set` favourable, `play: true` · 2. "Two simple flows" — "The rose line is Couette flow u = Uy/h, the orange curve
  Poiseuille flow. Nothing else is possible here." `terms: true` · 3. "They add" — "Because the equation (8.4a)
  0 = −(1/ρ)∂p/∂x + ν d²u/dy² is linear, the profile is exactly their sum (8.5)." `eq: 'cp'`, `code: {lines: [3, 3]}`,
  `derive: {id: 'D04', step: 9}` · 4. "Push back" — "Drag dp/dx upward. The floor slope shrinks and reaches zero at
  2μU/h² = 2 Pa/m." `set` {dpdx: 2}, `controls: ['dpdx']`, `readouts: ['tau0']`, `derive: {id: 'D05', step: 7}` · 5.
  "Backflow" — "Past the threshold the amber layer opens: fluid next to the floor flows backwards although Q is still
  positive." `set` strong adverse, `inspect: true` · 6. "Your turn" — "Predict the dp/dx that stops the net flow, then press
  'zero net flow' to check." `controls: ['dpdx', 'h']`.
- **Equations:** `red` "Reduced momentum" ref 'Eq. (8.4a)' `0=-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{d^2u}{dy^2}` ·
  `cp` "Couette–Poiseuille profile" ref 'Eq. (8.5)' `u(y)=\frac Uhy-\frac1{2\mu}\frac{dp}{dx}y(h-y)` live "u(h/2) = …" ·
  `Q` "Flow rate" (unnumbered) `Q=\frac{Uh}{2}\big[1-\frac{h^2}{6\mu U}\frac{dp}{dx}\big]` · `tau` "Shear stress"
  `\tau=\mu\frac{du}{dy}=\mu\frac Uh-\big(\frac h2-y\big)\frac{dp}{dx}` · `bf` "Backflow (ours)" `\frac{dp}{dx}>
  \frac{2\mu U}{h^2}` · symbols U, h, μ, dp/dx, u, Q, V, τ with units.
- **Check yourself:** (1) "Halve the gap at fixed U. What happens to the backflow threshold?" — "It quadruples: 2μU/h²
  ∝ 1/h²." `set {h: 0.005}` · (2) "At dp/dx = 2μU/h², where is the fluid fastest?" — "At the top wall (u = U): the
  profile is vertical at the floor and bends toward the wall speed." · (3) "With U = 0, is there a threshold?" — "No: the
  pure Poiseuille profile is forward everywhere for dp/dx < 0 and backward everywhere for dp/dx > 0." `set {U: 0}` · (4)
  "Why does Q stay positive although the floor layer reverses?" — "Q = Uh/2 − h³(dp/dx)/12μ changes sign only at 6μU/h²,
  three times the backflow threshold."
- **Selftest parity rows:** `{name: 'u mid', js: u(0.005, S0), py: 'ch08.channel_flow(0.005, 0.01, U=0.1, dpdx=4.0)', rtol:
  1e-12}` · `{name: 'Q', js: flow(S0).Q, py: 'ch08.channel_flow_rate(0.01, U=0.1, dpdx=4.0)[0]', rtol: 1e-12}` · `{name:
  'tau floor', js: tau(0, S0), py: 'ch08.channel_shear_stress(0.0, 0.01, U=0.1, dpdx=4.0)', rtol: 1e-12}` · `{name: 'thr',
  js: thr(S0), py: 'ch08.channel_backflow_threshold(0.1, 0.01, 0.001)', rtol: 1e-12}` · `{name: 'y_rev', js:
  state(S0).y_reversal, py: 'ch08.couette_poiseuille_state(0.01, 0.1, 4.0)["y_reversal"]', rtol: 1e-12}` · `{name: 'G
  alias', js: u(0.003, {…S0, dpdx: -2}), py: 'ch08.channel_flow(0.003, 0.01, U=0.1, G=2.0)', rtol: 1e-12}` (S0 = {h: 0.01,
  U: 0.1, dpdx: 4, mu: 0.001}).
- **Fit plan:** 360×640: status one line, `channel` (45 %) above `profile` (55 %), `stress` hidden (its numbers in the
  profile title), presets as a wrapping chip row, walkthrough card paged. 844×390: `channel` | `profile` side by side.
  1000×700 and desktop: rows [1.1, 1] with `stress` under `profile`.

### E2 · lubrication_scaling
- **Title:** "Why can a thin film ignore inertia?" · **Summary:** "Scale the gap equations with two lengths and every term
  becomes a bar: inertia is weighed by ε²Re_L, the pressure and the cross-gap friction stay, and across the gap only the
  pressure is left." · **CORE:** C05 (also N23, N24, N26–N31, R13) · **Reference:** `fid_formula_lab.html` (formula terms
  as bars, stage chips) with the numbered live explanation of `amplitude_phase_second_order_II_3.html`.
- **meta:** `viz:order 2` · `viz:sections 8.3` · `viz:equations 8.13 8.14 8.15 8.16 8.17` · `viz:fluidpy
  ch08.lubrication_scales ch08.lubrication_term_magnitudes` · `viz:derivations D10 D11`.
- **Physics:** `scales(s)` → {eps, ReL, eps2ReL, Lambda, pvisc} ↔ `ch08.lubrication_scales(L, h, U, rho, mu)`; `terms(s,
  pscale)` → the eight coefficients ↔ `ch08.lubrication_term_magnitudes(L, h, U, rho, mu, p_scale)`.
- **Views** (rows [1, 1.2]): 1. `gap` "The gap (height × k)" (row 0) — the passage L long, h high drawn with a stated
  exaggeration factor k = 0.15 L/h (title "height × 150"), a Couette-plus-Poiseuille profile at three stations (blue), small
  v arrows ~ εU (teal, their length in the same exaggeration), dimension lines for h and L. Static except when the
  controls move. 2. `xbars` "x-momentum (8.16a)" (row 1, flex 1) — horizontal log bars from 10⁻¹² to 10³: inertia
  ε²Re_L (teal), pressure 1/Λ (orange), streamwise diffusion ε² (light rose), cross-gap diffusion 1 (rose); a vertical
  line at 1 and a band "dropped (≤ 1 %)" below 0.01. 3. `ybars` "y-momentum (8.16b)" (row 1, flex 1, `hidePortrait`) — the
  same for ε⁴Re_L, 1/Λ, ε⁴, ε². Portrait: `gap` + `xbars`; the y-equation verdict ("∂p/∂y = O(ε²) = …") in the `xbars`
  title.
- **Controls:** `h` "Gap $h$" log 1 µm … 5 mm, default 50 µm · `L` "Length $L$" log 1 mm … 1 m, default 5 cm · `U` "Speed
  $U$" log 1 mm/s … 30 m/s, default 5 m/s · `nu` "Kinematic viscosity $\nu$" log 10⁻⁶ … 10⁻² m²/s, default 5.75 × 10⁻⁵ ·
  `pscale` chips "μUL/h² (natural) · P_a (book)" (optional `rho` 800–1300 kg/m³, default 870).
- **Presets:** "engine film" {h: 50e-6, L: 0.05, U: 5, nu: 5.75e-5} (ε²Re_L = 4.35e-3) · "Hele-Shaw cell" {h: 1e-3, L: 0.1,
  U: 0.01, nu: 1e-3} (10⁻⁴) · "knee joint (rough)" {h: 1e-6, L: 0.01, U: 0.1, nu: 1e-5} (10⁻⁶; labelled order-of-magnitude)
  · "thick gap — fails" {h: 2e-3, L: 0.01, U: 1, nu: 1e-6} (ε = 0.2, Re_L = 10⁴, ε²Re_L = 400).
- **Status:** "✅ lubrication valid: ε²Re_L = … ≪ 1 and ε = … ≪ 1" (ε²Re_L < 0.01 and ε < 0.05) · "⚠️ marginal: ε²Re_L =
  …" (0.01–1) · "❌ inertia matters: ε²Re_L = … > 1" or "❌ not thin: ε = …".
- **Readouts:** "ε = h/L" · "Re_L" · "ε²Re_L" · "Λ (with P_a)" · "Pressure scale μUL/h²" (atm).
- **Depth features:** Explain + Code + Derivation · **terms** (the bars — click a bar to highlight its term in the equation
  card) · **presets** · **status** · **linked views** (3) · **modes** (`pscale` switches the pressure normalisation; bars and
  Λ recompute).
- **Explain:** 0. *What the views show* — "**The gap** is drawn with its height stretched; teal arrows are the small
  cross-gap velocity v ~ εU. **The bars** are the coefficients of each term of (8.16a) and (8.16b) on a log scale — the
  size of each term, since every starred derivative is of order one." 1. *Two length scales* — "ε = h/L = **live**; Re_L =
  UL/ν = **live**; v ~ εU = **live** m/s (continuity: U/L = v/h)." 2. *Inertia's real weight* — "ε²Re_L = **live**" boxed;
  why: the x-equation was multiplied by ε²L²/(νU) so that the cross-gap friction has coefficient 1 (D10 step 8). 3. *The
  pressure scale* — "with P_a: Λ = μUL/(P_a h²) = **live** (the book assumes ~1); with μUL/h² = **live** Pa (**live** atm)
  the pressure coefficient is exactly 1." 4. *The x-balance* — table of the four coefficients, the kept ones bold. 5. *The
  y-balance* (hint on phones: "turn sideways for the y-bars") — "only 1/Λ has no ε: ∂p/∂y = O(ε²) = **live** — the pressure
  is uniform across the gap (8.17b) $0\cong-\frac1\rho\frac{\partial p}{\partial y}$." 6. *Reading the current setting* —
  valid: "Inertia is ε²Re_L = … of friction: the flow is a local balance of pressure gradient and cross-gap friction,
  (8.17a) $0\cong-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{\partial^2u}{\partial y^2}$ — even though Re_L alone looks
  large." · marginal: "Inertia is a few percent: the lubrication profile is a first guess; corrections of order ε²Re_L
  matter for precise bearings." · fails: "The gap is not thin enough or the flow is too fast: inertia competes with friction
  and the full Navier–Stokes equations are needed (a separating, recirculating gap flow)."
- **Derivation tab:** **D10** (14 steps, ★★★) `view: 'xbars'` (phones), goal `set` preset "engine film"; steps 5–10 each
  highlight the bar they create (`highlight: ['term:x_inertia']` at step 9, `['term:x_pressure']` at step 10); step 13
  highlights `ybars` pressure; step 12 `watch` "the printed ∂p/∂x would not give this bar set". **D11** (7 steps) `view:
  'xbars'`; step 2 `set` {…engine film} with **live** "ε²Re_L = …"; step 6 `watch` "units: ν∂²u/∂y² is m/s², like
  (1/ρ)∂p/∂x".
- **Code:**
  ```python
  L, h, U, rho, mu = {{L}}, {{h}}, {{U}}, {{rho}}, {{mu}}   # SI
  s = ch08.lubrication_scales(L, h, U, rho, mu)     # eps = {{eps}}, Re_L = {{ReL}}
  print(s["eps2_Re_L"])                             # inertia weight ε²Re_L = {{e2}}
  print(s["Lambda"], s["p_visc"])                   # Λ = {{Lam}}, μUL/h² = {{pv}} Pa
  t = ch08.lubrication_term_magnitudes(L, h, U, rho, mu, p_scale="{{ps}}")
  print(t["x_inertia"], t["y_pressure"])            # {{xi}}, {{yp}}
  ```
- **Walkthrough (6 steps):** 1. "A paradox" — "Re_L is in the thousands, yet bearing films are treated as purely viscous. Why?"
  `set` engine film · 2. "Two rulers" — "Along the gap things change over L, across it over h = εL; continuity then makes
  v ~ εU." `derive: {id: 'D10', step: 4}` · 3. "Weigh the terms" — "Scaled so friction across the gap is 1, inertia gets the
  coefficient ε²Re_L = 0.0044." `terms: true`, `derive: {id: 'D10', step: 9}` · 4. "Across the gap" — "In the y-equation
  every term but pressure carries ε² or less, so the pressure is the same across the film (8.17b)." `eq: 'l2'` · 5. "Break it"
  — "Press 'thick gap': ε = 0.2, Re_L = 10⁴ — inertia jumps above friction." `set` thick gap · 6. "Your turn" — "Predict how
  much faster the engine oil could move before ε²Re_L reaches 1, then test with U." `controls: ['U', 'h']`.
- **Equations:** (8.14) scalings · (8.15) $\frac{\partial u^*}{\partial x^*}+\frac{\partial v^*}{\partial y^*}=0$ · (8.16a)
  with live coefficients · (8.16b) with live coefficients · `l1` ref 'Eq. (8.17a)' $0\cong-\frac1\rho\frac{\partial p}{
  \partial x}+\nu\frac{\partial^2u}{\partial y^2}$ (note: "printed without ν") · `l2` ref 'Eq. (8.17b)' · symbols.
- **Check yourself:** (1) "Double U. By what factor does the inertia bar grow?" — "2: ε²Re_L ∝ U." · (2) "Halve h at fixed
  L, U, ν. What happens to ε²Re_L and to the pressure scale?" — "ε²Re_L falls 4×; μUL/h² grows 4×: thinner films are more
  viscous-dominated and build higher pressures." `set {h: 25e-6}` · (3) "Why is Λ ≈ 49 not a problem?" — "Λ only rescales
  p*; choosing the natural scale μUL/h² makes it exactly 1 without changing which terms are small." · (4) "Which is the
  relevant smallness condition: Re_L ≪ 1 or ε²Re_L ≪ 1?" — "ε²Re_L ≪ 1 (with ε ≪ 1)."
- **Selftest parity rows:** `{name: 'eps2ReL', js: scales(S0).eps2ReL, py: 'ch08.lubrication_scales(0.05, 5e-5, 5.0, 870.0,
  0.05)["eps2_Re_L"]', rtol: 1e-12}` · `{name: 'Lambda', js: scales(S0).Lambda, py: 'ch08.lubrication_scales(0.05, 5e-5,
  5.0, 870.0, 0.05)["Lambda"]', rtol: 1e-12}` · `{name: 'y diff', js: terms(S0, 'atm').y_diff_across, py:
  'ch08.lubrication_term_magnitudes(0.05, 5e-5, 5.0, 870.0, 0.05, p_scale="atm")["y_diff_across"]', rtol: 1e-12}` ·
  invariant `{name: 'natural Lambda = 1', js: terms(S0, 'viscous').x_pressure, expect: 1, rtol: 1e-12}`.
- **Fit plan:** 360×640: `gap` (35 %) + `xbars` (65 %), `ybars` hidden; presets wrap. Landscape phones: `gap` | `xbars`.
  Desktop: `gap` over `xbars` | `ybars`.

### E3 · slider_bearing
- **Title:** "How does a film thinner than a hair carry a load?" · **Summary:** "The same flux must pass a narrowing gap,
  so pressure rises inside: drag the taper and the sliding direction and watch the hump, the gap profiles and the load
  respond — and the printed formula fail the pressure equation." · **CORE:** C07, C06 (also N25, N32–N35, N105) ·
  **Reference:** `amplitude_phase_second_order_II_3.html` (system + response curves + a numbered derivation with live
  numbers).
- **meta:** `viz:order 3` · `viz:sections 8.3` · `viz:equations 8.18 8.19` · `viz:fluidpy ch08.slider_bearing
  ch08.slider_bearing_load ch08.slider_bearing_state ch08.slider_gap_velocity ch08.lubrication_velocity
  ch08.lubrication_flux ch08.slider_optimum_taper` · `viz:derivations D12 D13 D14`.
- **Physics:** `pExact(x, s)`, `pLin`, `pBook` ↔ `ch08.slider_bearing(x, h0, alpha, L, U, mu, model)`; `W(s, model)` ↔
  `ch08.slider_bearing_load` (series for ∣α∣ < 1e-3); `state(s)` ↔ `ch08.slider_bearing_state`; `ugap(x, y, s, frame)` ↔
  `ch08.slider_gap_velocity`; `q(h, dpdx, Uh, U0)` ↔ `ch08.lubrication_flux`; α_opt = 1.1889 from a local golden-section
  search ↔ `ch08.slider_optimum_taper()["alpha_opt"]`.
- **Views** (rows [1, 1]): 1. `pad` "Pad and film (height × k), pad frame" (row 0, flex 1.4) — the pad (grey block) from
  x = 0 (narrow end) to L (wide end, h₀(1 + α)), the floor (hatched) sliding at −U (arrows), profiles u − U at 7 stations
  (blue arrows), tracers (teal) moving through the gap, recirculation near the pad at the inlet shaded amber when α > 1,
  the load arrow W above the pad (length ∝ W, amber and pointing up if W < 0); pointer: click → station for the inspector.
  2. `press` "p(x) − p_e" (row 1, flex 1.2) — MPa vs x/L: exact (orange), linear in α (orange dashed), as printed (muted
  dotted) when the model chip "compare" is on, p_max dot, the end condition p(L) = p_e marked with a ring on each curve.
  3. `load` "W(α)" (row 1, flex 1, `hidePortrait`) — exact (blue) and linear (dashed) for α ∈ [−0.9, 3], optimum marked,
  current dot, W < 0 region amber. Portrait: `pad` + `press`; W exact/linear in the `press` title.
- **Controls:** `alpha` "Taper $\alpha$" −0.9 … 3, step 0.01, default 0.5 · `U` "Pad speed $U$" −10 … 10 m/s, default 5 ·
  `h0` "Minimum gap $h_o$" 10 … 200 µm, default 50 · `mu` "Viscosity $\mu$" (optional) 0.005 … 0.5 Pa s, default 0.05 ·
  `model` chips "exact · linear · compare (with printed)". L = 5 cm fixed (stated).
- **Presets:** "small taper" {alpha: 0.1} · "optimum" {alpha: 1.1889} · "reversed motion" {U: −5} · "parallel plates"
  {alpha: 0} · "steep taper" {alpha: 3}.
- **Status:** "🏋️ load carried: W = … kN/m (… tonnes per metre)" · "⚠️ W < 0: the pad is sucked down (αU < 0)" · "— no
  taper, no load (α = 0)" · suffix "↺ inlet recirculation (α > 1)".
- **Readouts:** "Load $W$ (exact)" (kN/m) · "Load, linear in α" · "Peak pressure" (MPa, atm) · "Pad-frame flux $C_1$" (cm²/s).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **presets** · **status** · **inspector** (click a
  station: "h = h₀(1 + αx/L) = …; dp/dx = −12μC₁/h³ − 6μU/h² = …; u(y) = −(h²/2μ)(dp/dx)(y/h)(1 − y/h) + U y/h; q_pad = C₁
  = …") · **modes** (exact / linear / printed ghost).
- **Explain:** 0. *Views* — "**Pad and film**: in the pad's frame the floor slides backwards and drags oil into the gap from
  the wide end; arrows are u − U, the speed relative to the pad. **p(x)**: orange = exact pressure (with the corrected,
  squared denominator), dashed = the small-α formula, dotted = as printed. **W(α)** (hidden on phones): load for every
  taper." 1. *The flux is the same at every station* — "C₁ = ∫₀^h(u − U)dy = −(1 + α)Uh₀/(2 + α) = **live** m²/s" boxed. 2.
  *Pressure* — "at your station x = **live**: p − p_e = (6μLU/h₀²)·α(x/L)(1 − x/L)/[(2 + α)(1 + αx/L)²] = **live** MPa; peak
  **pmax** at x/L = **xmax**" boxed. 3. *Load* — "exact W = (6μUL²/h₀²α²)[ln(1 + α) − 2α/(2 + α)] = **Wex**; linear
  αμL²U/(2h₀²) = **Wlin**; the linear formula is **err** % off" boxed. 4. *How big* — "μUL/h₀² = **pv** MPa: the natural
  pressure scale; the peak is **atm** atmospheres." 5. *Inlet recirculation* (ours) — "at the wide end the wall shear under
  the pad changes sign when 3α/(2 + α) > 1, i.e. α > 1 (for U > 0): you are **above/below**." 6. *The printed slip* — "the
  printed (1 + αx/L) to the first power still meets p(0) = p(L) = p_e, but its slope misses dp/dx = −12μC₁/h³ − 6μU/h²:
  its peak is **live** MPa against the exact **live** (the dotted curve)". 7. *Reading the current setting* — load: "Oil dragged into a narrowing gap cannot all get out at the
  narrow end, so pressure builds until Poiseuille backflow evens the flux: the hump carries W." · optimum: "Near 1 + α =
  2.19 the load is largest for given h₀, U, μ, L — the classic design ratio." · reversed: "Sliding the other way makes the
  gap widen along the drag: pressure falls below p_e and pulls the pad down; the bearing works in one direction only." · no
  taper: "Parallel walls: pure Couette flow, uniform pressure, no load."
- **Derivation tab:** **D12** (8 steps) `view: 'pad'`, step 5 highlights the pad wall, step 8 `watch` "with U₀ ≠ 0 the printed
  profile would not meet the upper wall". **D13** (12 steps, ★★★) `view: 'pad'`; step 2 (Leibniz) highlights the sloping
  upper wall; step 7 **live** "h_t + q_x = 0; steady → q = C₁ = …". **D14** (15 steps, ★★★) `view: 'press'`; step 10 **live**
  "C₁ = −(1 + α)Uh₀/(2 + α) = …"; step 13 `set` {model: 'compare'}, `watch` "the dotted printed curve fails dp/dx = …";
  step 15 **live** "W = αμL²U/(2h₀²) = …".
- **Code:**
  ```python
  h0, alpha, L, U, mu = {{h0}}, {{alpha}}, 0.05, {{U}}, {{mu}}
  x = np.linspace(0, L, 201)
  p = ch08.slider_bearing(x, h0, alpha, L, U, mu=mu)       # exact p − p_e, peak {{pmax}} Pa
  W = ch08.slider_bearing_load(h0, alpha, L, U, mu=mu)     # exact load {{Wex}} N/m
  Wl = ch08.slider_bearing_load(h0, alpha, L, U, mu=mu, model="linear")  # {{Wlin}} N/m
  st = ch08.slider_bearing_state(h0, alpha, L, U, mu=mu)   # C1 = {{C1}} m²/s
  print(st["inlet_backflow"])                               # {{ib}}
  ```
- **Walkthrough (7 steps):** 1. "Floating on oil" — "A pad slides on a film 50 µm thick. What holds it up?" `set` {alpha: 0.5}
  · 2. "Pad frame" — "Ride with the pad: the floor drags oil in from the wide end toward the narrow end." `play: true` · 3.
  "Same flux everywhere" — "Continuity across the gap (the Reynolds equation) forces one flux C₁ at every station."
  `derive: {id: 'D13', step: 7}` · 4. "A pressure hump" — "Pure drag would carry more oil than the narrow end lets out, so
  pressure rises inside and pushes back." `eq: 'pex'` · 5. "The load" — "The area under the hump is W: 32.8 kN/m here,
  3.3 tonnes per metre." `readouts: ['W']` · 6. "Wrong way" — "Flip U: the hump becomes suction and the pad is pulled down."
  `set` reversed motion · 7. "Your turn" — "Find the taper that carries the most load, then press 'optimum'."
  `controls: ['alpha']`.
- **Equations:** `prof` ref 'Eq. (8.19)' (consistent form, note on the printed U₀) · `flux` "Gap flux" (ours)
  $q=-\frac{h^3}{12\mu}\frac{\partial p}{\partial x}+\frac{(U_0+U_h)h}{2}$ · `C1` "Pad-frame flux" $C_1=-\frac{h^3}{12\mu}
  \frac{dp}{dx}-\frac{Uh}2$ · `pex` "Pressure (Example 8.1)" $p-p_e=\frac{6\mu LU}{h_o^2}\frac{\alpha(x/L)(1-x/L)}{(2+\alpha)
  (1+\alpha x/L)^2}$ · `W` "Load" $W=\frac{\alpha\mu L^2U}{2h_o^2}$ and the exact form.
- **Check yourself:** (1) "Halve h₀. What happens to W?" — "It grows 4× (∝ 1/h₀²): a heavier load squeezes the film and the
  bearing pushes back harder — it is stable." · (2) "Why does the linear formula overestimate W at α = 0.5?" — "It drops the
  (1 + αx/L)² in the denominator, which reduces the pressure over most of the pad." · (3) "At α = 0 what is the flow?" —
  "Couette only, uniform pressure, zero load." `set {alpha: 0}` · (4) "Where is the peak pressure for α = 1?" — "At x/L = 1/3,
  toward the narrow end." `set {alpha: 1}`.
- **Selftest parity rows:** `{name: 'p exact mid', js: pExact(0.025, S0), py: 'ch08.slider_bearing(0.025, 5e-5, 0.5, 0.05,
  5.0, mu=0.05)', rtol: 1e-12}` · `{name: 'W exact', js: W(S0, 'exact'), py: 'ch08.slider_bearing_load(5e-5, 0.5, 0.05, 5.0,
  mu=0.05)', rtol: 1e-10}` · `{name: 'W small alpha', js: W({…S0, alpha: 1e-4}, 'exact'), py: 'ch08.slider_bearing_load(5e-5,
  1e-4, 0.05, 5.0, mu=0.05)', rtol: 1e-8}` · `{name: 'C1', js: state(S0).C1, py: 'ch08.slider_bearing_state(5e-5, 0.5, 0.05,
  5.0, mu=0.05)["C1"]', rtol: 1e-12}` · `{name: 'u gap', js: ugap(0.03, 3e-5, S0, 'ground'), py:
  'ch08.slider_gap_velocity(0.03, 3e-5, 5e-5, 0.5, 0.05, 5.0, mu=0.05)', rtol: 1e-10}` · `{name: 'alpha opt', js: aopt(),
  py: 'ch08.slider_optimum_taper()["alpha_opt"]', rtol: 1e-6}` · invariant `{name: 'p(L)=pe', js: pExact(0.05, S0), expect:
  0, atol: 1e-6}`.
- **Fit plan:** 360×640: `pad` (45 %) + `press` (55 %); `load` hidden. Landscape: `pad` | `press`. Desktop: `pad` full width,
  `press` | `load` below.

### E4 · viscous_gravity_current
- **Title:** "Why does a honey drop spread like t^(1/5)?" · **Summary:** "A thin viscous layer spreads under its own weight
  by a diffusion that chokes itself as it thins; every initial shape ends on one dome whose front moves as t^(1/5)." ·
  **CORE:** C08, C10 (also N37, N63, N64, N105) · **Reference:** `angular_frequency_explorer_1.html` (system animation +
  graphs on one clock, presets, end-of-run summary card).
- **meta:** `viz:order 4` · `viz:sections 8.3 8.4` · `viz:equations` (none numbered: Example 8.3 and 8.7 forms, written
  out) · `viz:fluidpy ch08.thin_film_flux ch08.viscous_current_similarity ch08.thin_film_state` · `viz:derivations D15 D23`.
- **Physics:** `flux(h, hx, s)` ↔ `ch08.thin_film_flux`; `huppert(x, t, s)` ↔ `ch08.viscous_current_similarity(x, t, area,
  rho, g, mu)`; `state(t, s)` ↔ `ch08.thin_film_state`; the JS march: explicit conservative finite volume on 240 cells of a
  symmetric half domain with sub-steps at dt = 0.2 dx²/(β h_max³), precursor film 10⁻⁴ × initial height, run in the
  dimensionless time τ = βh₀³t/x₀² and shown in real seconds (`Viz.fmtTime`); volume printed each frame.
- **Views** (rows [1.1, 1]): 1. `layer` "The layer h(x, t)" (row 0, flex 1.4) — x ∈ [−X, X] with X = 1.3 × the current
  similarity front, h as a blue fill, the similarity dome at the same t as a purple dashed ghost, flux arrows at 9 points
  (length ∝ q), front markers (amber), volume in the title; pointer: click → probe x. 2. `front` "Front x_N(t), log–log"
  (row 1, flex 1) — the simulated front (teal dots, "so far" bold), the law x_N = η_N(βA³t)^{1/5} (purple dashed), a
  slope-1/2 guide (muted) — diffusion's slope — labelled "if it were ordinary diffusion". 3. `rescaled` "h t^(1/5) vs
  x/t^(1/5)" (row 1, flex 1, `hidePortrait`) — the last 6 recorded profiles rescaled with the trial exponent m (chip 1/5 or
  1/2) collapsing (or not). Portrait: `layer` + `front`; the collapse spread in the front title.
- **Controls:** `shape` chips "box · Gaussian · two humps" · `A` "Half-area $A$" log 10⁻⁵ … 10 m², default 10⁻⁴ (1 cm × 1 cm)
  · `mu` "Viscosity $\mu$" log 0.1 … 10⁴ Pa s, default 1 (glycerol) · `rho` "Density $\rho$" (optional) 900 … 2700, default
  1260 · `m` chips "rescale with 1/5 · 1/2" (optional) · transport.
- **Transport:** `tau` (log clock) from 10⁻² to 10⁴ in τ, rate 0.5 decades/s, `end: 'hold'` with a `Viz.card` summary: "front
  slope = 0.20 (fitted over the last decade); volume change 0.0000 %; profile within … % of the similarity dome".
- **Presets:** "honey on a plate" {mu: 10, rho: 1400, A: 2e-4} · "glycerol bead" {mu: 1, rho: 1260, A: 1e-4} · "lava lobe"
  {mu: 1e3, rho: 2700, A: 10} · "ice-like slab (analogy)" {mu: 1e4, rho: 900, A: 10} with a note "real ice is
  non-Newtonian; this only shows the time scale" · "wrong exponent" {m: '1/2'}.
- **Status:** "⏳ still remembers its initial shape (profile … % from the dome)" (> 5 %) · "✅ self-similar: front ∝ t^0.20,
  volume conserved to 10⁻¹²" · always the volume.
- **Readouts:** "Front $x_N$ (simulated)" · "Front from the law" · "Centre height" · "Effective diffusivity ρgh³/3μ" (m²/s) ·
  "Volume 2A".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** · **status** ·
  **inspector** (click x: "h = …, ∂h/∂x = …, q = −(ρg/3μ)h³∂h/∂x = −4118 × … = … m²/s").
- **Explain:** 0. *Views* — "**The layer**: blue = the simulated thickness, purple dashed = the similarity dome at the same
  time, arrows = the flux. **Front**: teal dots are the simulated front on log–log axes; the purple line has slope 1/5,
  the grey one slope 1/2. **Rescaled** (hidden on phones): profiles multiplied by t^(1/5) against x/t^(1/5)." 1. *The flux*
  — "q = −(ρg/3μ)h³∂h/∂x: at your probe **live**" boxed; why: the half-parabola profile of D15. 2. *Volume* — "∂h/∂t +
  ∂q/∂x = 0 in flux form, so 2A = **live** m² stays fixed (change **live**)". 3. *The front law* — "x_N = η_N(βA³t)^{1/5} with
  β = ρg/3μ = **live**, η_N = 1.411: **law** vs simulated **sim**" boxed. 4. *Why so slow* — "effective diffusivity at the
  centre ρgh_c³/3μ = **live** m²/s — it falls as h_c ∝ t^(−1/5) cubed, so spreading chokes itself." 5. *At the current
  time* — "t = **live**; the centre is **live** mm thick". 6. *Reading the current setting* — early: "The layer still
  remembers its initial shape; the front has not yet joined the t^(1/5) line." · self-similar: "Every initial shape has
  forgotten itself: the profile is Huppert's dome and the front advances as t^(1/5) — the exponents forced by the equation
  and the conserved volume (Example 8.7)." · wrong exponent chosen: "With m = 1/2 (diffusion's guess) the rescaled curves
  drift apart: this is not a diffusion length."
- **Derivation tab:** **D15** (10 steps) `view: 'layer'`; step 3 `watch` "under the crest the hydrostatic pressure is
  highest"; step 9 **live** "q = −(ρg/3μ)h³h_x = …". **D23** (9 steps) `view: 'front'` (phones), step 9 `watch` "the teal dots
  follow the slope-1/5 line" with `set` {m: '1/5'}.
- **Code:**
  ```python
  rho, g, mu, A = {{rho}}, 9.80665, {{mu}}, {{A}}          # SI; A = half-area [m²]
  st = ch08.thin_film_state({{t}}, A, rho=rho, g=g, mu=mu)  # front {{xN}} m, centre {{hc}} m
  x = np.linspace(-st["x_N"], st["x_N"], 101)
  h = ch08.viscous_current_similarity(x, {{t}}, A, rho=rho, g=g, mu=mu)   # Huppert dome
  q = ch08.thin_film_flux(h[60], np.gradient(h, x)[60], rho=rho, g=g, mu=mu)  # flux {{q}} m²/s
  ```
- **Walkthrough (6 steps):** 1. "Honey on a plate" — "Why does a spreading drop slow down so much?" `set` glycerol bead ·
  2. "Weight pushes" — "Thicker places have higher hydrostatic pressure, so fluid flows from thick to thin." `play: true`,
  `derive: {id: 'D15', step: 9}` · 3. "Flux chokes" — "The flux carries h³: as the layer thins it slows down drastically."
  `readouts: ['Deff']` · 4. "One fifth" — "On log–log axes the front follows slope 1/5, not diffusion's 1/2." highlight
  `view:front` · 5. "Shape forgotten" — "Start from two humps: they merge and end on the same dome." `set {shape: 'two'}`,
  `play: true` · 6. "Your turn" — "A lava lobe is 1000× more viscous. Predict the time to reach the same front, then play."
  `set` lava lobe.
- **Equations:** "Thin-film equation (Example 8.3)" $\frac{\partial h}{\partial t}=\frac{\rho g}{3\mu}\frac{\partial}{\partial x}
  \big(h^3\frac{\partial h}{\partial x}\big)$ · "Flux" $\int_0^hu\,dy\cong-\frac{\rho g}{3\mu}h^3\frac{\partial h}{\partial x}$ ·
  "Similarity form (Example 8.7)" $h(x,t)=At^{-1/5}F(x/Dt^{1/5})$ · "Front (Huppert, ours)" $x_N=\eta_N(\beta A^3t)^{1/5}$.
- **Check yourself:** (1) "μ × 1000: how much longer to reach the same front?" — "1000×: t enters only as ρgt/μ." · (2)
  "Volume × 32: how much farther is the front at a given time?" — "32^{3/5} = 8×: x_N ∝ A^{3/5}." · (3) "Why is the front's
  slope not 1/2?" — "The diffusivity ρgh³/3μ falls as the layer thins; with constant volume the similarity exponents are
  n = m = 1/5." · (4) "Is volume exactly conserved by the simulation?" — "Yes, to round-off: the scheme updates by flux
  differences."
- **Selftest parity rows:** `{name: 'flux', js: flux(0.01, -0.5, S0), py: 'ch08.thin_film_flux(0.01, -0.5, rho=1260.0,
  g=9.80665, mu=1.0)', rtol: 1e-12}` · `{name: 'huppert h', js: huppert(0.02, 100, S0), py:
  'ch08.viscous_current_similarity(0.02, 100.0, 1e-4, rho=1260.0, g=9.80665, mu=1.0)', rtol: 1e-10}` · `{name: 'front',
  js: state(100, S0).x_N, py: 'ch08.thin_film_state(100.0, 1e-4, rho=1260.0, g=9.80665, mu=1.0)["x_N"]', rtol: 1e-10}` ·
  invariants `{name: 'volume conserved', js: marchVolumeDrift(), expect: 0, atol: 1e-12}` · `{name: 'late front', js:
  marchFront(1e3)/state(1e3, S0).x_N, expect: 1, rtol: 3e-2}`.
- **Fit plan:** 360×640: `layer` (55 %) + `front` (45 %); `rescaled` hidden. Landscape: `layer` | `front`. Desktop: `layer`
  above `front` | `rescaled`.

### E5 · stokes_first_problem
- **Title:** "How can profiles at all times be one curve?" · **Summary:** "A plate starts suddenly; momentum diffuses a
  distance √(νt); plot against y/√(νt) and every profile, for any speed and viscosity, collapses onto 1 − erf(η/2)." ·
  **CORE:** C09 (also N38–N56, R14) · **Reference:** `forced_damped_vibrations.html` (explanation panel, linked views on one
  time slider) with the raw/rescaled collapse pattern of `knowledge/viz_patterns.md`.
- **meta:** `viz:order 5` · `viz:sections 8.4` · `viz:equations 8.20 8.21 8.22 8.23 8.24 8.25 8.26 8.27 8.28 8.29 8.30 8.31`
  · `viz:fluidpy ch08.stokes_first_problem ch08.similarity_variable ch08.diffusion_thickness ch08.stokes_first_vorticity
  ch08.stokes_first_stopped ch08.stokes_first_state ch08.couette_startup_profile` · `viz:derivations D16 D17 D18 D19 D20`.
- **Physics:** `u(y, t, s)` = U erfcHP(y/2√(νt)) ↔ `ch08.stokes_first_problem`; `eta(y, t, s, half)` ↔
  `ch08.similarity_variable`; `d99(t, s, level)` = 2 erfcinvHP(level)√(νt) (Newton on erfcHP) ↔ `ch08.diffusion_thickness`;
  `w(y, t, s)` ↔ `ch08.stokes_first_vorticity`; `stopped(y, t, T, s)` ↔ `ch08.stokes_first_stopped`; `couette(y, t, h, s)`
  (200-term series) ↔ `ch08.couette_startup_profile(y, t, U, h, nu)`; `state(t, s)` ↔ `ch08.stokes_first_state`; a JS
  Crank–Nicolson march (200 points, Thomas algorithm) drawn as teal dots.
- **Views** (rows [1, 1]): 1. `fluid` "Above the plate" (row 0, flex 1.3) — y from 0 to 1.3 δ₉₉(t_end), the plate sliding at
  U (arrow), tracers at 8 heights moving at u(y, t), a dyed vertical line bending, vorticity shading (purple alpha ∝ ω);
  δ₉₉ as an amber line. 2. `profile` "u(y)" (row 0, flex 1) — raw mode: u/U vs y with faint earlier profiles (every 0.5
  decade) and the current one bold, CN dots; rescaled mode: all vs η (axis switch "η = y/√(νt) · η/2 (the book's figure
  axis)"), one purple curve 1 − erf(η/2); pointer: click → probe. 3. `thick` "δ₉₉(t), log–log" (row 1, `hidePortrait`) —
  slope-½ line for the current ν plus ghosts for water/air/honey, the current point. Portrait: `fluid` + `profile`; δ₉₉ in
  the profile title.
- **Controls:** `U` "Plate speed $U$" 0.01 … 1 m/s, default 0.1 · `nu` "Viscosity $\nu$" log 10⁻⁶ … 10⁻³ m²/s (chips water ·
  air · honey), default 10⁻⁶ · `mode` chips "raw · rescaled" · `level` chips "99 % · 95 %" (optional) · transport.
- **Transport:** `t` log clock 0.1 → 1000 s (in units scaled so the layer fills the view), `end: 'hold'`, summary card "δ₉₉
  grew ×100 in 10⁴× the time — slope ½".
- **Presets:** "water" · "air (15× faster)" {nu: 1.5e-5} · "honey" {nu: 1e-3} · "stop the plate at T" {stopT: t/2 of the run}
  (uses `stopped`; status: "⚠️ an imposed time T breaks similarity") · "Couette start-up (second wall)" {wall: 2 cm} (uses
  `couette`; collapse holds until δ₉₉ ~ h, then fails).
- **Status:** "✅ self-similar: profiles collapse in η (spread < 10⁻¹²)" · "⚠️ imposed scale (T or h): collapse breaks".
- **Readouts:** "√(νt)" · "δ₉₉ = 3.643√(νt)" · "Wall stress" (Pa) · "∫ω dy" (= U).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **modes** (raw / rescaled) ·
  **presets** · **inspector** (click a point: "η = y/√(νt) = 0.01/0.01 = 1.000; u/U = erfc(η/2) = erfc(0.5) = **0.4795**") ·
  **status**.
- **Explain:** 0. *Views* — "**Above the plate**: tracers and a dye line show the moving layer; purple shading is vorticity,
  amber the 99 % edge. **u(y)**: raw profiles, or all of them against η, where they coincide. **δ₉₉(t)** (hidden on
  phones): slope ½." 1. *The diffusion length* — "√(νt) = √(**nu** × **t**) = **live** m" boxed. 2. *Your probe* — "η =
  y/√(νt) = **live**; u/U = 1 − erf(η/2) = erfc(**live**) = **live**" boxed (8.30). 3. *Thickness* — "δ = 2 erfc⁻¹(level)√(νt)
  = **3.643** × **live** = **live** m" (8.31). 4. *Wall stress and vorticity* — "τ_w = μU/√(πνt) = **live** Pa; ∫₀^∞ω dy = U
  = **live** m/s (the book prints −U; it is +U)". 5. *At the current time* — live t, the dye displacement at the probe.
  6. *Reading the current setting* — self-similar: "Nothing in the problem sets a length except √(νt), so the profile can
  only stretch: one curve in η for every time, speed and fluid." · stopped: "Stopping the plate at T adds a time scale:
  after T the layer detaches from the wall and the curves no longer coincide." · Couette: "A second wall adds a length h:
  while δ₉₉ ≪ h the profiles still collapse; once the layer reaches the wall they approach the straight Couette line."
- **Derivation tab:** **D16** (6) `view: 'fluid'`; step 5 `watch` "far above the plate the fluid is still at rest". **D17** (8)
  `view: 'profile'`, step 7 `set` {mode: 'rescaled'}. **D18** (9) `view: 'profile'`, step 9 `watch` "t → 0 and y → ∞ both send
  η to ∞". **D19** (11) `view: 'profile'`; step 8 (ξ = 2ζ) `set` {mode: 'rescaled', axis: 'half'}, `watch` "the book's axis
  η/2"; step 11 **live** "u/U at your probe = …". **D20** (5) `view: 'profile'`; step 4 **live** "δ₉₉ = 3.643 × √(νt) = …".
- **Code:**
  ```python
  U, nu, t = {{U}}, {{nu}}, {{t}}                  # m/s, m²/s, s
  y = {{y}}                                        # probe height [m]
  print(ch08.similarity_variable(y, t, nu))        # η = {{eta}}
  u = ch08.stokes_first_problem(y, t, U=U, nu=nu)  # (8.30): {{u}} m/s
  d = ch08.diffusion_thickness(t, nu)              # (8.31): {{d99}} m
  st = ch08.stokes_first_state(t, U=U, nu=nu)      # τ_w = {{tw}} Pa
  ```
- **Walkthrough (7 steps):** 1. "Yank a plate" — "How far does the motion reach after one second, one minute?" · 2. "It
  diffuses" — "Only (8.20) ∂u/∂t = ν∂²u/∂y² is left: momentum spreads like heat." `play: true`, `derive: {id: 'D16', step: 6}`
  · 3. "No ruler" — "The only length is √(νt), so u/U = F(y/√(νt)) (8.25)." `derive: {id: 'D17', step: 7}` · 4. "One curve"
  — "Switch to rescaled: every profile lands on 1 − erf(η/2) (8.30)." `set {mode: 'rescaled'}` · 5. "How thick" — "δ₉₉ =
  3.643√(νt) (8.31): 3.6 cm after 100 s in water." `readouts: ['d99']` · 6. "Break it" — "Stop the plate: after T the curves
  part." `set` stop preset · 7. "Your turn" — "Predict δ₉₉ in air after 100 s, then choose air and check." `set` water.
- **Equations:** (8.20), (8.21)–(8.23), (8.24), (8.25), (8.26), (8.27)–(8.28), (8.29), (8.30) with erf's definition, (8.31),
  ∫ω dy = U — each written out, live values on (8.30) and (8.31).
- **Check yourself:** (1) "t × 4: how does δ₉₉ change?" — "×2 (∝ √t)." · (2) "Air instead of water at the same t?" — "√15 ≈
  3.9 × thicker." · (3) "Does U change the rescaled curve?" — "No: u/U is linear in U (the reason y/Ut dropped out)." · (4)
  "Why does stopping the plate break the collapse?" — "T is an imposed time scale: t/T is a second variable."
- **Selftest parity rows:** `{name: 'u', js: u(0.01, 100, S0), py: 'ch08.stokes_first_problem(0.01, 100.0, U=0.1,
  nu=1e-6)', rtol: 1e-10}` · `{name: 'd99', js: d99(100, S0, 0.01), py: 'ch08.diffusion_thickness(100.0, 1e-6)', rtol: 1e-9}`
  · `{name: 'omega', js: w(0.005, 100, S0), py: 'ch08.stokes_first_vorticity(0.005, 100.0, U=0.1, nu=1e-6)', rtol: 1e-10}` ·
  `{name: 'stopped', js: stopped(0.01, 150, 100, S0), py: 'ch08.stokes_first_stopped(0.01, 150.0, 100.0, U=0.1, nu=1e-6)',
  rtol: 1e-9}` · `{name: 'couette', js: couette(0.005, 50, 0.02, S0), py: 'ch08.couette_startup_profile(0.005, 50.0, 0.1,
  0.02, 1e-6)', rtol: 1e-8}` · `{name: 'tau_w', js: state(100, S0).tau_w, py: 'ch08.stokes_first_state(100.0, U=0.1,
  nu=1e-6)["tau_w"]', rtol: 1e-12}` · invariant `{name: 'CN vs erfc', js: cnMaxErr(), expect: 0, atol: 2e-3}`.
- **Fit plan:** 360×640: `fluid` (45 %) + `profile` (55 %); `thick` hidden. Landscape: side by side. Desktop: `fluid` |
  `profile` on top, `thick` full width below at 0.8 height.

### E6 · similarity_exponents
- **Title:** "Where do similarity exponents come from?" · **Summary:** "Guess γ = At^(−n)F(ξ/δ(t)) with δ ∝ t^m: the rescaled
  profiles collapse only at the right n and m, fixed by matching powers of t and one conserved quantity." · **CORE:** C10
  (also N57–N64, N106) · **Reference:** `stride_padding_playground.html` (a formula with numbers plugged in, badges that
  explain the result, presets at the classic settings) with the modes of `angular_frequency_explorer_1.html`.
- **meta:** `viz:order 6` · `viz:sections 8.4` · `viz:equations 8.20 8.25 8.32` · `viz:fluidpy ch08.similarity_collapse_error
  ch08.vortex_sheet_diffusion ch08.line_vortex_decay ch08.stokes_first_problem ch08.viscous_current_similarity
  ch08.transition_width ch08.temporal_bl_wall_stress` · `viz:derivations D21 D22 D23`.
- **Physics:** exact solutions per mode (JS mirrors with `erfcHP`): `stokes1(y, t)` ↔ `ch08.stokes_first_problem`;
  `sheet(y, t)` → {u, ω} ↔ `ch08.vortex_sheet_diffusion`; `vortex(r, t)` ↔ `ch08.line_vortex_decay`; `bead(x, t)` ↔
  `ch08.viscous_current_similarity`; `spread(mode, n, m, times)` (max ∣difference∣ of rescaled profiles on 201 points of
  ξ/t^m, normalised by the profile maximum) ↔ `ch08.similarity_collapse_error(case, n, m, times)`; `brackets(mode, n, m)`
  → the powers of t of each bracket in the reduced equation and of the conserved integral.
- **Modes** (`mode`, chips): "impulsive plate" (γ = u/U, ξ = y; n fixed by u(0) = U; δ ∝ t^m) · "vortex sheet" (γ = ω_z) ·
  "line vortex" ((8.32b): γ = u_θ, prefactor r^(−n)) · "spreading bead" (γ = h).
- **Views** (rows [1, 1]): 1. `raw` "Profiles at five times" (row 0) — γ vs ξ at t = 1, 3, 10, 30, 100 (units per mode; blue
  shades). 2. `scaled` "Rescaled with your n, m" (row 1, flex 1.2) — γ t^n (or γ ξ^n) vs ξ/t^m; the five curves, the
  spread printed in the title; purple dashed = the true similarity curve when `reveal` is on. 3. `powers` "Powers of t"
  (row 1, flex 1, `hidePortrait`) — chips for each bracket of the reduced equation showing its power of t with the trial
  exponents (green when all equal, red otherwise), and a small plot of the conserved integral (∫ω dy, ∫h dx, circulation)
  vs t (flat only at the right n). Portrait: `raw` + `scaled`; a one-line bracket verdict in the `scaled` title.
- **Controls:** `mode` chips · `n` "Trial $n$" −0.5 … 1.5, step 0.01 · `m` "Trial $m$ (δ ∝ t^m)" 0 … 1, step 0.01 · `reveal`
  toggle "show the true curve".
- **Presets:** "correct exponents" (per mode: n = 0, m = ½ · n = ½, m = ½ · n = 1 (space power), m = ½ · n = m = 1/5) · "near
  miss" {n: 0.4} · "diffusion guess" {m: 0.5} (bead mode).
- **Status:** "✅ collapsed: n = …, δ ∝ t^… (spread 10⁻¹⁵)" · "❌ spread = … — the brackets scale as t^…, t^…".
- **Readouts:** "Spread" · "Bracket powers" · "Conserved integral change" (%) · mode-specific: "Width 5.544√(νt)", "C_f
  temporal", "front t^(1/5)".
- **Depth features:** Explain + Code + Derivation · **modes** (4 systems) · **linked views** (3) · **presets** · **status** ·
  **terms** (bracket powers as chips).
- **Explain:** 0. *Views* — "**Profiles**: the exact solution at five times. **Rescaled**: the same curves with your trial
  exponents — one curve means you found the similarity form. **Powers of t** (hidden on phones): each bracket of the reduced
  equation, and the conserved quantity." 1. *The ansatz* — "(8.32a) γ = At^(−n)F(ξ/δ(t)), δ = Dt^m (or (8.32b) with ξ^(−n))."
  2. *Powers of t with your numbers* — table "bracket | power" (e.g. vortex sheet: n/t → t^(−1), δ′/δ → t^(−1), ν/δ² →
  t^(−2m)); "equal iff 2m = 1". 3. *The conserved quantity* — "vortex sheet: −∫ω dy = −At^(−n)δ∫F dη ∝ t^(m − n): constant
  iff n = m" (bead: ∫h dx ∝ t^(m − n); line vortex: circulation at infinity fixes n = 1 in (8.32b)). 4. *Your spread* —
  "**live**" boxed. 5. *Mode results* — vortex sheet: "u = U erf(y/2√(νt)), width 5.544√(νt) = **live**, C_f = (2/√π)√(ν/U²t)
  = **live**"; line vortex: "u_θ = (Γ/2πr)[1 − e^(−r²/4νt)] — the Gaussian vortex with σ² = 4νt"; bead: "n = m = 1/5".
  6. *Reading the current setting* — collapsed: "Your exponents make every bracket scale alike and the conserved quantity
  constant: nothing is left to choose — this is the similarity solution." · wrong m: "The brackets scale differently, so t
  cannot be divided out; the rescaled curves fan out." · wrong n: "The brackets may match but the conserved integral drifts
  in time: the amplitude law is wrong."
- **Derivation tab:** **D21** (8) `view: 'scaled'`, goal `set` {mode: 'impulsive plate'}; step 6 **live** "δ′/δ = C₁ν/δ² ⇒
  powers −1 and −2m". **D22** (14, ★★★) `set` {mode: 'vortex sheet'}; step 9 (n = ½ from the jump) `set` {n: 0.5, m: 0.5},
  `watch` "the conserved-integral line goes flat"; step 13 **live** "AD = −U/√(πν) = …". **D23** (9) `set` {mode: 'spreading
  bead'}; step 9 `set` {n: 0.2, m: 0.2}.
- **Code:**
  ```python
  times = [1, 3, 10, 30, 100]                     # s
  err = ch08.similarity_collapse_error("{{case}}", {{n}}, {{m}}, times)
  print(err)                                      # spread = {{err}}
  u, w = ch08.vortex_sheet_diffusion(1e-3, 1.0, U=0.01, nu=1e-6)   # {{u}}, {{w}}
  print(ch08.transition_width(1.0, 1e-6))        # 5.544 mm at t = 1 s
  ```
- **Walkthrough (6 steps):** 1. "Three problems, one trick" — "A shear layer, a vortex, a bead: all self-similar. How do we
  find their laws?" · 2. "Guess the form" — "Try γ = At^(−n)F(ξ/δ), δ ∝ t^m (8.32a). Two unknown exponents." · 3. "Match the
  powers" — "Every bracket must carry the same power of t: for the plate, 2m = 1." `set {mode: 'impulsive plate', n: 0, m:
  0.5}`, `derive: {id: 'D21', step: 6}` · 4. "A conserved quantity" — "For the vortex sheet the jump 2U is constant; that
  forces n = m = ½." `set {mode: 'vortex sheet', n: 0.4}`, `derive: {id: 'D22', step: 9}` · 5. "The bead" — "Volume is
  conserved and the bracket balance gives 3n + 2m = 1: n = m = 1/5." `set {mode: 'spreading bead'}` · 6. "Your turn" — "In
  line-vortex mode, find the exponents before pressing 'correct'." `controls: ['n', 'm']`.
- **Equations:** (8.32a) $\gamma=At^{-n}F(\xi/\delta(t))\equiv At^{-n}F(\eta)$ · (8.32b) $\gamma=A\xi^{-n}F(\xi/\delta(t))$ ·
  Example 8.4 bracket form $-\big[\frac1\delta\frac{d\delta}{dt}\big]\eta\frac{dF}{d\eta}=\big[\frac\nu{\delta^2}\big]\frac{d^2F}
  {d\eta^2}$ · Example 8.5 result $\omega_z=-\frac{U}{\sqrt{\pi\nu t}}e^{-y^2/4\nu t}$, $u=U\,\mathrm{erf}\frac{y}{2\sqrt{\nu t}}$ ·
  Example 8.6 $u_\theta=\frac{\Gamma}{2\pi r}[1-e^{-r^2/4\nu t}]$ · Example 8.7 $h=At^{-1/5}F(x/Dt^{1/5})$.
- **Check yourself:** (1) "Why must n = 0 for the impulsive plate?" — "u = U at η = 0 for every t; a factor t^(−n) would change
  the wall speed." · (2) "In vortex-sheet mode set n = 0.5, m = 0.4. Which condition fails?" — "The bracket ν/δ² scales as
  t^(−0.8), the others as t^(−1)." · (3) "Why is the bead's m not ½?" — "Its diffusivity ∝ h³ falls with time; with constant
  volume the exponent equations give 1/5." · (4) "Which quantity fixes n for the line vortex?" — "The circulation Γ far away
  (the initial Γ/2πr), hence the space prefactor r^(−1) of (8.32b)."
- **Selftest parity rows:** `{name: 'spread bead right', js: spread('spreading', 0.2, 0.2, T5), py:
  'ch08.similarity_collapse_error("spreading", 0.2, 0.2, [1, 3, 10, 30, 100])', atol: 1e-12}` · `{name: 'spread sheet
  wrong', js: spread('vortex_sheet', 0.4, 0.5, T5), py: 'ch08.similarity_collapse_error("vortex_sheet", 0.4, 0.5, [1, 3, 10,
  30, 100])', rtol: 1e-6}` · `{name: 'sheet u', js: sheet(1e-3, 1).u, py: 'ch08.vortex_sheet_diffusion(1e-3, 1.0, U=0.01,
  nu=1e-6)[0]', rtol: 1e-10}` · `{name: 'vortex', js: vortex(0.01, 100), py: 'ch08.line_vortex_decay(0.01, 100.0, 1e-3,
  nu=1e-6)', rtol: 1e-12}` · `{name: 'width', js: width(1), py: 'ch08.transition_width(1.0, 1e-6)', rtol: 1e-9}`.
- **Fit plan:** 360×640: `raw` (40 %) + `scaled` (60 %); `powers` hidden (verdict in the title); mode chips wrap. Desktop:
  `raw` | `scaled` over `powers`.

### E7 · oscillating_plate
- **Title:** "How deep does a shaking plate reach?" · **Summary:** "Each layer repeats the wall's motion delayed by y/δ radians
  and shrunk by e^(−y/δ), δ = √(2ν/ω): a 'wave' made only of diffusion." · **CORE:** C11 (also N65–N74) · **Reference:**
  `amplitude_phase_second_order_II_3.html` (amplitude and phase lag as the story; "evaluate → amplitude → phase → time lag →
  at time t").
- **meta:** `viz:order 7` · `viz:sections 8.5` · `viz:equations 8.33 8.34 8.35 8.36 8.37 8.38` · `viz:fluidpy
  ch08.stokes_second_problem ch08.stokes_layer ch08.stokes_layer_state` · `viz:derivations D24 D25`.
- **Physics:** `u(y, t, s)` ↔ `ch08.stokes_second_problem(y, t, U, omega, nu)`; `layer(s)` ↔ `ch08.stokes_layer(nu, omega)`;
  `state(y, s)` ↔ `ch08.stokes_layer_state(nu, omega, y)`.
- **Views** (rows [1, 1]): 1. `plate` "The plate and the fluid" (row 0, flex 1.3) — y from 0 to 6δ_e (unit switch), the
  plate moving ±U (rose arrow), tracers at 7 heights oscillating horizontally, a dyed vertical line bending, the probe
  height as a teal line. 2. `profile` "u(y, t)" (row 0, flex 1) — u/U vs y/δ_e, current profile blue, the envelope ±e^(−y/δ_e)
  muted dashed, the four book phases ωt = 0, π/2, π, 3π/2 as faint ghosts, the book depth 4√(ν/ω) = 2√2 δ_e with 0.059
  (amber), the growing-root ghost (muted dotted, clipped) on D24 step 7. 3. `clock` "u(t) at the wall and the probe" (row 1,
  `hidePortrait`) — two time series over two periods with the lag marked by an arrow. Portrait: `plate` + `profile`; the lag
  in the profile title ("probe: 17 % · 0.28 s late").
- **Controls:** `omega` "Frequency $\omega$" log 10⁻⁵ … 10⁴ rad/s, default 2π · `nu` "Viscosity $\nu$" log 10⁻⁶ … 10⁻¹ m²/s,
  default 10⁻⁶ · `yp` "Probe height $y$ / δ_e" 0 … 5, default 1.77 (1 mm at the defaults) · transport.
- **Transport:** `tau` in periods 0 → 4, rate 0.5 period/s, `end: 'loop'`; real time t = τT (`Viz.fmtTime`).
- **Presets:** "plate 1 Hz in water" · "M₂ tide, molecular ν" {omega: 1.405e-4, nu: 1e-6} (δ_e 0.12 m) · "M₂ tide, eddy ν"
  {nu: 1e-2} (12 m) · "sound 1 kHz in air" {omega: 6283, nu: 1.5e-5} (0.069 mm) · "Ekman teaser" (labelled "same (1 + i)/δ
  structure, Ch. 13 — not this physics").
- **Status:** "🌊 δ_e = √(2ν/ω) = …; the wall is felt to ≈ 4√(ν/ω) = …" · probe "at y = … the motion is …% of the wall's and …
  rad (… s) late".
- **Readouts:** "δ_e" · "Book depth 4√(ν/ω)" · "Amplitude at probe" · "Time lag at probe" · "Crest speed √(2νω)".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** · **inspector** (click a
  height: "amplitude e^(−y/δ_e) = e^(−1.772) = 0.170; phase lag y/δ_e = 1.772 rad; time lag = 1.772/ω = 0.282 s") · **notes**
  ("Right now": a table of real Stokes layers — lab plate, M₂ tide molecular/eddy, acoustic — with ω, ν, δ_e, the nearest
  row highlighted).
- **Explain:** 0. *Views* — "**Plate**: rose = wall motion, tracers follow at their height. **u(y, t)**: blue = now, dashed =
  the envelope ±e^(−y/δ_e), faint = the four phases of Fig. 8.16. **u(t)** (hidden on phones): wall vs probe." 1. *Evaluate
  the depth* — "δ_e = √(2ν/ω) = √(2 × **nu**/**omega**) = **live**" boxed; "book: 4√(ν/ω) = 2√2 δ_e = **live**, where the
  amplitude is e^(−2√2) = 0.05911". 2. *Amplitude at your probe* — "e^(−y/δ_e) = **live**". 3. *Phase and time lag* —
  "y/δ_e = **live** rad → **live** s late (period **T**)". 4. *Crest speed and 'wavelength'* — "√(2νω) = **live** m/s; 2πδ_e =
  **live** — but no restoring force: diffusion only." 5. *At the current time* — "ωt = **live**; wall u = U cos ωt =
  **live**; probe u = **live**". 6. *Reading the current setting* — "Faster shaking or a less viscous fluid confines the
  layer (δ_e ∝ √(ν/ω)); the lag grows linearly with height, so the profile looks like a wave travelling up at √(2νω). With
  eddy viscosity the tidal layer is metres thick — the bottom boundary layer of shelf seas; the Ekman layer of Ch. 13 has the
  same (1 + i)/δ form."
- **Derivation tab:** **D24** (10) `view: 'profile'`; step 7 (B = 0) `set` shows the growing root ghost with `watch` "the dotted
  curve blows up — not bounded"; step 10 **live** "u(y_p, t) = … ". **D25** (6) `view: 'profile'`; step 2 **live** "e^(−2√2)
  = 0.05911"; step 4 `watch` "the zero crossing climbs at √(2νω)".
- **Code:**
  ```python
  nu, omega, U = {{nu}}, {{omega}}, 1.0
  L = ch08.stokes_layer(nu, omega)          # δ_e = {{de}} m, book depth {{db}} m
  y = {{yp}}                                # probe [m]
  st = ch08.stokes_layer_state(nu, omega, y)   # amplitude {{amp}}, lag {{lag}} s
  u = ch08.stokes_second_problem(y, {{t}}, U=U, omega=omega, nu=nu)   # (8.38): {{u}}
  ```
- **Walkthrough (6 steps):** 1. "Shake a plate" — "How far above a shaking plate does the fluid move?" · 2. "It lags" —
  "Watch the tracers: higher ones move less and later." `play: true` · 3. "Complex trick" — "Write u = Re{e^(iωt)f(y)} (8.35):
  the PDE becomes iωf = νf″ (8.36)." `derive: {id: 'D24', step: 3}` · 4. "The layer" — "Only the decaying root survives:
  u = Ue^(−y/δ)cos(ωt − y/δ) (8.38), δ = 0.56 mm here." `eq: 'sol'` · 5. "Not a wave" — "The crest climbs at √(2νω) but no
  restoring force acts: it is diffusion with a delay." `notes: true` · 6. "Your turn" — "Predict the tidal layer depth with
  eddy viscosity 10⁻² m²/s, then press the preset." `controls: ['nu', 'omega']`.
- **Equations:** (8.33) $u(y=0,t)=U\cos(\omega t)$ · (8.34) bounded · (8.35) $u=\mathrm{Re}\{e^{i\omega t}f(y)\}$ · (8.36)
  $i\omega f=\nu(d^2f/dy^2)$ · $k=\pm(i+1)(\omega/2\nu)^{1/2}$ · (8.37) · `sol` (8.38) with live values.
- **Check yourself:** (1) "ω × 4: δ_e?" — "halves (∝ ω^(−1/2))." · (2) "At what height is the motion exactly half a period late?"
  — "y = πδ_e, where the amplitude is e^(−π) = 4 %." · (3) "Why is B = 0 in (8.37)?" — "e^(+(1+i)y/δ) grows without bound."
  · (4) "Why is this not self-similar?" — "ω imposes a time scale; y and t never combine."
- **Selftest parity rows:** `{name: 'u', js: u(1e-3, 0.3, S0), py: 'ch08.stokes_second_problem(1e-3, 0.3, U=1.0,
  omega=2*np.pi, nu=1e-6)', rtol: 1e-12}` · `{name: 'de', js: layer(S0).delta_e, py: 'ch08.stokes_layer(1e-6,
  2*np.pi)["delta_e"]', rtol: 1e-12}` · `{name: 'time lag', js: state(1e-3, S0).time_lag, py: 'ch08.stokes_layer_state(1e-6,
  2*np.pi, 1e-3)["time_lag"]', rtol: 1e-12}` · invariant `{name: 'book depth amplitude', js: Math.exp(-4*Math.sqrt(1/2)) -
  state(layer(S0).delta_book, S0).amplitude, expect: 0, atol: 1e-14}`.
- **Fit plan:** 360×640: `plate` (45 %) + `profile` (55 %); `clock` hidden. Landscape: side by side. Desktop: `plate` |
  `profile`, `clock` below.

### E8 · stokes_sphere_flow
- **Title:** "What does creeping flow round a sphere look like?" · **Summary:** "Stokes' flow is fore–aft symmetric and decays
  only like a/r; beyond r ~ a/Re_a inertia returns and Oseen's wake appears — compare with ideal flow." · **CORE:** C13, C15
  (also N81–N87, N93–N100, R17–R19) · **Reference:** `angular_frequency_explorer_1.html` (modes = one object in different
  models, linked views, presets, "Right now" notes).
- **meta:** `viz:order 8` · `viz:sections 8.6` · `viz:equations 8.43 8.44 8.45 8.46 8.47 8.48 8.49 8.53` · `viz:fluidpy
  ch08.stokes_sphere_streamfunction ch08.stokes_sphere_velocity ch08.oseen_streamfunction ch08.oseen_velocity
  ch08.side_line_speed ch08.inertia_viscous_ratio` · `viz:derivations D28 D29 D32 D33`.
- **Physics:** `psiS(r, th, frame)` ↔ `ch08.stokes_sphere_streamfunction`; `velS(r, th, frame)` ↔ `ch08.stokes_sphere_velocity`;
  `psiO(r, th, Re, frame)` (−expm1 form) ↔ `ch08.oseen_streamfunction`; `velO` ↔ `ch08.oseen_velocity`; ideal sphere
  ψ = ½U sin²θ(r² − a³/r) and its velocity (Ch. 6); `side(r, model, Re)` ↔ `ch08.side_line_speed`; `ratio(r, th, Re)` (central
  differences on the Cartesian Stokes field, h = 10⁻⁴r) ↔ `ch08.inertia_viscous_ratio`. U = 1, a = 1 (dimensionless; ν = 2/Re).
- **Views** (rows [1.2, 1]): 1. `flow` "Streamlines" (row 0, flex 1.4, `equal: true`) — the meridian plane x ∈ [−Z, Z],
  ϖ ∈ [−Z, Z] (Z = zoom r/a, 3 … 100, log slider), ψ contours (`Viz.field.contour`, 24 levels, mirrored), speed heatmap
  (optional), tracers (fluid or body frame) advected by RK4, the sphere, the circle r = a/Re_a (amber dashed, drawn when
  inside the window); pointer: click → probe. 2. `side` "Speed along the side line θ = π/2" (row 1, flex 1) — ∣u∣/U vs r/a on
  log x (1 … 100): Stokes (blue), ideal (muted dashed), Oseen (purple) at the current Re; the probe radius marked. 3.
  `ratio` "Inertia / viscous" (row 1, flex 1, `hidePortrait`) — log–log vs r/a: the computed ratio (teal dots) and Re_a r/a
  (dashed), the line ratio = 1. Portrait: `flow` + `side`; "inertia = viscosity at r ≈ …a" in the side title.
- **Controls:** `model` chips "Stokes · Oseen · ideal" · `frame` chips "body · fluid" · `Re` "Reynolds number $2aU/\nu$" log
  10⁻³ … 5, default 0.1 · `zoom` "View radius r/a" log 3 … 100, default 8 · transport (tracers).
- **Presets:** "Stokes, Re → 0" {model: 'Stokes', frame: 'body', zoom: 8} · "fluid frame loops" {frame: 'fluid'} · "Oseen,
  Re = 1 (wake)" {model: 'Oseen', Re: 1, zoom: 30} · "ideal flow for contrast" {model: 'ideal'} · "far field" {zoom: 100}.
- **Status:** "↔️ fore–aft symmetric (Stokes, linear)" · "🌊 Oseen wake: inertia wins beyond r ≈ … a" · "ideal: no viscosity,
  symmetric, compact".
- **Readouts:** "Side speed at probe" · "Deficit at 10a" (%) · "Crossover r/a = 1/Re_a" · "ψ at probe".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (Stokes / Oseen / ideal; body / fluid
  frame) · **transport** · **presets** · **inspector** (click: "u_r = U cos θ(1 − 3a/2r + a³/2r³) = … ; u_θ = −U sin θ(1 −
  3a/4r − a³/4r³) = …") · **status**.
- **Explain:** 0. *Views* — "**Streamlines** in a plane through the axis (the stream flows left to right in the body frame);
  amber circle = where inertia catches up. **Side-line speed**: blue Stokes, dashed ideal, purple Oseen. **Inertia/viscous**
  (hidden on phones)." 1. *ψ and velocity at your probe* — (8.48) and (8.49) with numbers, boxed. 2. *How far the sphere is
  felt* — "at 10a Stokes 0.925U (7.5 % deficit), ideal 1.0005U" live. 3. *Where Stokes fails* — "inertia/viscous ~ Re_a r/a
  = 1 at r/a = **live**" boxed. 4. *Oseen's correction* (Oseen mode) — "(8.53) with Re = **live**: the wake width …; near
  the sphere ψ_Oseen/ψ_Stokes = **live**". 5. *Reading the current setting* — Stokes: "Linear equations: reversing the stream
  reverses every arrow, so front and back look alike and there is no wake; the a/r term makes the disturbance reach tens of
  radii." · Oseen: "Far away inertia (U∂u′/∂x) sweeps the disturbance downstream: streamlines crowd behind the sphere — a
  wake — while near the sphere it is Stokes' flow again." · ideal: "Without viscosity the disturbance decays like a³/r³ and
  the fluid speeds up beside the sphere (1.5U at the equator), instead of slowing down."
- **Derivation tab:** **D28** (12, ★★★) `view: 'flow'`; step 6 **live** "ω_φ at the probe = −(3Ua/2r²) sin θ = …". **D29** (10)
  `view: 'flow'`; step 8 (A = 0) `set` {zoom: 100}, `watch` "an r⁴ term would outgrow the uniform stream"; step 10 **live**
  "C = −3Ua/4, D = Ua³/4". **D32** (7) `view: 'side'` (phones) / `ratio`; step 7 **live** "r/a = 1/Re_a = …". **D33** (7) `set`
  {model: 'Oseen', Re: 0.01}, step 7 `watch` "Oseen and Stokes streamlines coincide near the sphere".
- **Code:**
  ```python
  U, a, Re = 1.0, 1.0, {{Re}}
  r, th = {{r}}, {{th}}                               # probe (θ from downstream)
  print(ch08.stokes_sphere_streamfunction(r, th, U, a))   # (8.48): {{psi}}
  print(ch08.stokes_sphere_velocity(r, th, U, a))         # (8.49): {{ur}}, {{ut}}
  print(ch08.side_line_speed(10.0, U, a, model="stokes"))  # 0.9248 at 10a
  print(ch08.oseen_streamfunction(r, th, U, a, Re))       # (8.53): {{psio}}
  ```
- **Walkthrough (7 steps):** 1. "A bead in syrup" — "What does the flow round a sinking bead look like?" · 2. "Symmetric" —
  "Front and back look the same: the Stokes equations (8.43) ∇p = μ∇²u are linear." `set` Stokes · 3. "Loops" — "In the fluid
  frame particles are pushed aside and return." `set` fluid frame loops, `play: true` · 4. "Felt far away" — "At 10a the speed
  is still 7.5 % low; ideal flow has recovered." highlight `view:side` · 5. "Where it fails" — "Inertia ~ Re_a r/a equals
  viscosity at r ≈ a/Re_a." `derive: {id: 'D32', step: 7}` · 6. "Oseen's wake" — "Oseen keeps U∂u′/∂x: at Re = 1 the
  streamlines crowd behind." `set` Oseen Re = 1 · 7. "Your turn" — "Predict where the amber circle sits at Re = 0.02, then
  zoom out to check." `controls: ['Re', 'zoom']`.
- **Equations:** (8.43), (8.44) with the E² note, (8.45)–(8.47), (8.48), (8.49), inertia/viscous ~ Re r/a, (8.53) with
  Re = 2aU/ν.
- **Check yourself:** (1) "Reverse the stream: what happens to the Stokes streamlines?" — "Same lines, arrows reversed
  (linearity)." · (2) "Why is the side speed below U in Stokes flow but above U in ideal flow?" — "Viscosity drags the fluid
  near the sphere along with it; ideal fluid must speed up to get round." · (3) "Re = 0.01: where does inertia matter?" —
  "Beyond r ≈ 200a (Re_a = 0.005)." · (4) "Does Oseen's solution change the drag?" — "Yes, C_D = (24/Re)(1 + 3Re/16) — see E9."
- **Selftest parity rows:** `{name: 'psi', js: psiS(2, 1.0, 'body'), py: 'ch08.stokes_sphere_streamfunction(2.0, 1.0, 1.0,
  1.0)', rtol: 1e-12}` · `{name: 'u_theta', js: velS(2, 1.0, 'body')[1], py: 'ch08.stokes_sphere_velocity(2.0, 1.0, 1.0,
  1.0)[1]', rtol: 1e-12}` · `{name: 'oseen psi', js: psiO(3, 2.0, 0.5, 'body'), py: 'ch08.oseen_streamfunction(3.0, 2.0, 1.0,
  1.0, 0.5)', rtol: 1e-12}` · `{name: 'oseen u', js: velO(3, 2.0, 0.5, 'body')[0], py: 'ch08.oseen_velocity(3.0, 2.0, 1.0,
  1.0, 0.5)[0]', rtol: 1e-10}` · `{name: 'side', js: side(10, 'stokes'), py: 'ch08.side_line_speed(10.0, 1.0, 1.0,
  model="stokes")', rtol: 1e-12}` · `{name: 'ratio', js: ratio(50, 1.5708, 0.2), py: 'ch08.inertia_viscous_ratio(50.0,
  1.5708, 1.0, 1.0, 10.0)', rtol: 1e-4}` (ν = 2aU/Re = 10).
- **Fit plan:** 360×640: `flow` (60 %, square) + `side` (40 %); `ratio` hidden; mode chips wrap. Landscape: `flow` | `side`.
  Desktop: `flow` left full height, `side` over `ratio` right.

### E9 · stokes_drag_settling
- **Title:** "Where does 6πμaU come from?" · **Summary:** "Sweep round the sphere: the pressure and friction tractions add to a
  uniform push, collecting 2πμaU and 4πμaU; then drag the particle size and watch the fall speed grow as a²." · **CORE:**
  C14 (also N88–N92, R20, N99) · **Reference:** `fid_formula_lab.html` (term bars adding to a total, click to isolate) with the
  live explanation of `forced_damped_vibrations.html`.
- **meta:** `viz:order 9` · `viz:sections 8.6` · `viz:equations 8.50 8.51 8.52` · `viz:fluidpy ch08.stokes_sphere_pressure
  ch08.stokes_sphere_surface_stresses ch08.stokes_drag ch08.stokes_drag_running ch08.settling_state ch08.terminal_velocity
  ch08.stokes_drag_coefficient ch08.oseen_drag_coefficient` · `viz:derivations D30 D31`.
- **Physics:** `p(r, th, s)` ↔ `ch08.stokes_sphere_pressure`; `stresses(th, s)` → [σ_rr, σ_rθ, t_x] ↔
  `ch08.stokes_sphere_surface_stresses`; `running(th, s)` → {pressure, friction, total} ↔ `ch08.stokes_drag_running`;
  `settle(a, particle)` ↔ `ch08.settling_state(a, rho_p, rho, mu)`; `cdS(Re)`, `cdO(Re)` ↔ `ch08.stokes_drag_coefficient`,
  `oseen_drag_coefficient`; the Morrison correlation mirrored from `core.similarity.sphere_drag_coefficient` (ch04 E-code).
- **Views** (rows [1, 1]): 1. `sphere` "Tractions on the surface" (row 0, flex 1) — the sphere (meridian circle) with arrows at
  24 angles: pressure part (orange, normal), friction part (rose, tangential), their sum (black, all the same x-length), the
  sweep marker at θ; the swept arc highlighted; stream U from the left. 2. `curves` "Along the surface" (row 0, flex 1.2) —
  p − p∞ and σ_rθ in units of μU/a vs θ from 0 (rear) to π (front) (orange, rose), and on a second axis the running drag
  P(θ), F(θ) in μaU reaching 2π and 4π (dashed), the sweep line. 3. `settle` "Fall speed vs radius" / "C_D(Re)" (row 1,
  `hidePortrait`, mode chips) — U_t vs a on log–log for the four particle kinds with the Re = 0.1 validity limit, current
  particle dot; or C_D·Re/24 vs Re with Stokes, Oseen, Morrison. Portrait: `sphere` + `curves`; U_t and Re in the curves title.
- **Controls:** `theta` (the transport's sweep angle) · `particle` chips "cloud droplet · quartz sand · bacterium · oil drop" ·
  `a` "Radius $a$" log 1 µm … 1 mm, default 10 µm · `mu` "Viscosity $\mu$" (optional; set by the particle's fluid) · `view3`
  chips "settling · C_D".
- **Transport:** `theta` 0 → π, rate π/6 per s, `end: 'hold'` with a card "pressure 2πμaU (⅓) + friction 4πμaU (⅔) = 6πμaU".
- **Presets:** "cloud droplet 10 µm" · "drizzle 100 µm" (Re ≈ 16: invalid) · "silt 10 µm in a river" {particle: sand, a: 1e-5}
  (U_t 0.36 mm/s) · "Millikan drop" {particle: oil, a: 1e-6}.
- **Status:** "✅ Stokes regime: Re = … < 0.1" · "⚠️ Re = … > 0.1: Stokes underestimates the drag, so U_t is too high" ·
  sweep "collected so far: pressure …, friction …".
- **Readouts:** "Terminal speed $U_t$" · "Re = 2aU_tρ/μ" · "Drag D" · "Pressure part" · "Friction part".
- **Depth features:** Explain + Code + Derivation · **terms** (pressure ⅓ and friction ⅔ bars adding to 6πμaU; click to isolate
  in the sphere view) · **linked views** (3) · **transport** (sweep) · **presets** · **status** · **inspector** (click the
  sphere: "t_x = σ_rr cos θ − σ_rθ sin θ = (3μU/2a)cos²θ + (3μU/2a)sin²θ = 3μU/2a = …") · **notes** (table of real settling
  speeds: fog droplet, cloud droplet, silt, fine sand, bacterium, Millikan drop; current row highlighted).
- **Explain:** 0. *Views* — "**Tractions**: orange = pressure push, rose = friction, black = their sum (same x-part everywhere).
  **Along the surface**: the two stresses vs θ (θ = 0 is the rear) and the drag collected so far. **Fall speed** (hidden on
  phones)." 1. *Pressure at the swept point* — "(8.50) p − p∞ = −3μaU cos θ/(2a²) = **live** Pa; front +3μU/2a, rear
  −3μU/2a (the book drops this minus)". 2. *Friction* — "σ_rθ = −(3μU/2a) sin θ = **live**". 3. *The x-push* — "t_x = 3μU/2a
  = **live** Pa at every θ" boxed. 4. *Collected drag* — "P(θ) = πμaU(1 − cos³θ) = **live**; F(θ) = πμaU(2 − 3cos θ + cos³θ) =
  **live**; at θ = π: 2πμaU + 4πμaU = 6πμaU = **D**" boxed. 5. *Fall speed* — "U_t = 2(ρ′ − ρ)ga²/(9μ) = **live** m/s; Re =
  **live**; C_D = 24/Re = **live**" (8.52). 6. *Reading the current setting* — valid: "Drag grows in proportion to speed, so
  the particle settles when 6πμaU equals its weight in the fluid: U_t ∝ a² — halve the droplet, quarter the speed." ·
  invalid: "At Re above ~0.1 inertia adds drag (Oseen: (24/Re)(1 + 3Re/16)); the true fall speed is lower than Stokes'."
- **Derivation tab:** **D30** (11, ★★★) `view: 'curves'`; step 6 **live** "∂p/∂r = 3μUa cos θ/r³ = …"; step 8 `watch` "both
  components agree — g′(θ) = 0"; step 11 `set` {theta: 0} with **live** "p − p∞ = −3μU/2a = …". **D31** (13, ★★★) `view:
  'sphere'`; steps 11–13 `play` the sweep, `highlight: ['term:pressure']` at 12 and `['term:friction']` at 13.
- **Code:**
  ```python
  a, rho_p, rho, mu = {{a}}, {{rhop}}, {{rho}}, {{mu}}
  st = ch08.settling_state(a, rho_p, rho, mu)       # U_t = {{Ut}} m/s, Re = {{Re}}
  D = ch08.stokes_drag(mu, a, st["U_t"], parts=True)  # {{Dp}} + {{Df}} = {{D}} N
  th = {{th}}
  print(ch08.stokes_sphere_surface_stresses(th, st["U_t"], a, mu))   # σ_rr, σ_rθ, t_x
  print(ch08.stokes_drag_running(th, mu, a, st["U_t"]))              # collected so far
  ```
- **Walkthrough (6 steps):** 1. "Falling droplets" — "Why does a cloud droplet fall at only about a centimetre per second?" ·
  2. "Front and back" — "Pressure is high in front, low behind: +3μU/2a and −3μU/2a (8.50)." `derive: {id: 'D30', step: 11}` ·
  3. "Friction on the sides" — "Shear σ_rθ = −(3μU/2a) sin θ drags the surface along." `terms: true` · 4. "Add it up" — "Play
  the sweep: pressure gives 2πμaU, friction 4πμaU — total 6πμaU (8.51)." `play: true`, `derive: {id: 'D31', step: 13}` · 5.
  "Terminal speed" — "Drag ∝ U balances weight: U_t = 2(ρ′ − ρ)ga²/9μ = 1.2 cm/s." `readouts: ['Ut', 'Re']` · 6. "Your
  turn" — "Predict U_t for a 100 µm drizzle drop, then check the status." `set` drizzle.
- **Equations:** (8.50) $p-p_\infty=-\frac{3\mu aU\cos\theta}{2r^2}$ · surface stresses (ours) · (8.51) $D=6\pi\mu aU$ · terminal
  balance $(4/3)\pi a^3g(\rho'-\rho)=6\pi\mu aU$ · (8.52) $C_D=\frac{D}{\frac12\rho U^2\pi a^2}=\frac{24}{\mathrm{Re}}$ · Oseen
  $C_D=\frac{24}{\mathrm{Re}}(1+\frac3{16}\mathrm{Re})$.
- **Check yourself:** (1) "Double the radius: U_t?" — "×4 (∝ a²), and Re ×8." · (2) "Which contributes more drag, pressure or
  friction?" — "Friction, two thirds." · (3) "Is the x-traction larger at the front than at the sides?" — "No: it is 3μU/2a
  everywhere; pressure dominates the front, friction the sides." · (4) "For which particle is Stokes' law invalid at 100 µm?"
  — "The water droplet in air (Re ≈ 16); sand in water at 100 µm has Re ≈ 7 — also beyond 0.1."
- **Selftest parity rows:** `{name: 'p front', js: p(1e-5, Math.PI, S0), py: 'ch08.stokes_sphere_pressure(1e-5, np.pi, U=0.012,
  a=1e-5, mu=1.81e-5)', rtol: 1e-12}` · `{name: 'sigma rt', js: stresses(1.0, S0)[1], py:
  'ch08.stokes_sphere_surface_stresses(1.0, U=0.012, a=1e-5, mu=1.81e-5)[1]', rtol: 1e-12}` · `{name: 'running', js:
  running(Math.PI/2, S0).pressure, py: 'ch08.stokes_drag_running(np.pi/2, 1.81e-5, 1e-5, 0.012)["pressure"]', rtol: 1e-12}` ·
  `{name: 'Ut', js: settle(1e-5, 'droplet').U_t, py: 'ch08.settling_state(1e-5, 1000.0, 1.2, 1.81e-5)["U_t"]', rtol: 1e-12}` ·
  `{name: 'CD oseen', js: cdO(0.5), py: 'ch08.oseen_drag_coefficient(0.5)', rtol: 1e-12}` · invariant `{name: 'full sweep',
  js: running(Math.PI, S0).total / (6*Math.PI*S0.mu*S0.a*S0.U), expect: 1, rtol: 1e-12}`.
- **Fit plan:** 360×640: `sphere` (45 %) + `curves` (55 %); `settle` hidden (U_t and Re in the curves title). Landscape: side
  by side. Desktop: `sphere` | `curves` over `settle`.

### B1 · rotating_cylinders_couette
(Backup — built only if one of E1–E9 fails review.)
- **Title:** "What swirl fills the gap between two rotating cylinders?" · **Summary:** "Always AR + B/R: the walls choose A
  and B; send the outer wall to infinity and only a free vortex is left, remove the inner one and the fluid turns rigidly." ·
  **CORE:** C04 (also N17–N22, N104, R10–R12) · **Reference:** `angular_frequency_explorer_1.html` (rotating system + curve on
  one clock).
- **meta:** `viz:order 10` · `viz:sections 8.2` · `viz:equations 8.9 8.10 8.11 8.12` · `viz:fluidpy ch08.circular_couette
  ch08.circular_couette_pressure ch08.circular_couette_power ch08.circular_couette_state` · `viz:derivations D08 D09`.
- **Physics:** `uphi(R, s)` (exact limit branches for R₂ = ∞ and R₁ = 0) ↔ `ch08.circular_couette`; `pres(R, s)` ↔
  `ch08.circular_couette_pressure`; `state(s)` ↔ `ch08.circular_couette_state`.
- **Views:** (1) `annulus` top view with tracers moving at u_φ/R, the walls turning; (2) `profile` u_φ(R) with the AR (orange
  dashed) and B/R (rose dashed) parts and the two limit ghosts; (3) `rayleigh` (hidePortrait) (Ru_φ)² vs R — increasing =
  stable (Ch. 11 teaser).
- **Controls:** `Omega1` −2 … 2 rad/s, `Omega2` −2 … 2 rad/s, `ratio` R₁/R₂ 0.05 … 0.95, transport. **Presets:** inner only ·
  outer only · co-rotating (solid body) · counter-rotating · "R₂ → ∞" (free vortex) · "R₁ → 0" (solid body). **Status:**
  "rigid rotation (B = 0)" / "free vortex (A = 0)" / "mixed; Rayleigh stable/unstable". **Depth features:** linked views,
  transport, presets, status, inspector (click R: AR + B/R arithmetic).
- **Explain:** 0 views · 1 A and B from the walls (numbers) · 2 u_φ at the probe · 3 the pressure rise across the gap (N20) · 4
  torque and power = dissipation (R11) · 5 reading (which part dominates; the limits (8.11) $u_\varphi=\Omega_1R_1^2/R$ and
  (8.12) $u_\varphi=\Omega_2R$).
- **Derivation tab:** D08 (10), D09 (6) from Part F, `view: 'profile'`. **Code:** `circular_couette(R, R1, R2, O1, O2,
  return_coeffs=True)`, `circular_couette_state(...)`. **Walkthrough (5):** the question → two swirls → the walls choose →
  the two limits → your turn (predict the sign change of u_φ for counter-rotation). **Check yourself (3):** "Ω₁ = Ω₂: which
  flow?" — "solid body, B = 0" · "Outer at rest, R₂ → ∞: the vorticity?" — "zero everywhere (free vortex)" · "Why is the torque
  the same on both cylinders?" — "steady angular-momentum balance of the fluid ring".
- **Selftest parity rows:** `{name: 'u mid', js: uphi(0.015, S0), py: 'ch08.circular_couette(0.015, 0.01, 0.02, 1.0, 0.0)',
  rtol: 1e-12}` · `{name: 'free vortex', js: uphi(0.05, {…S0, R2: Infinity}), py: 'ch08.circular_couette(0.05, 0.01, np.inf,
  1.0, 0.0)', rtol: 1e-12}` · `{name: 'dp gap', js: state(S0).dp_gap, py: 'ch08.circular_couette_state(0.01, 0.02, 1.0,
  0.0)["dp_gap"]', rtol: 1e-10}`.
- **Fit plan:** 360×640: `annulus` + `profile`; `rayleigh` hidden.

---

## Part D — runtime budget (full run < 5 min on a laptop / Colab CPU)

Chapter 8 is mostly closed forms (µs per point) plus nine sympy engines and three numerical solvers. The costs: the
thin-film run (`thin_film_spread`, backward Euler + Picard on 400 cells to t = 1000 s ≈ 6 s, **cached** to
`outputs/ch08/thin_film_run.npz`), Crank–Nicolson marches (400 × 2000 steps ≈ 0.5 s; the convergence table 3 grids ≈ 1 s;
the Stokes-second check 10 periods ≈ 1 s), the Hele-Shaw Poisson solve (128² direct ≈ 2 s), Oseen contour grids (200² ×
5 Re ≈ 1 s), sympy (`parallel_flow_sympy` × 3 ≈ 2 s, `lubrication_nondim_sympy` ≈ 2 s, `reynolds_equation_sympy` ≈ 1 s,
`slider_bearing_sympy` ≈ 3 s, `similarity_reduce_sympy` × 5 ≈ 3 s, `stokes_second_sympy` < 1 s, `low_re_scaling_sympy` ≈ 1 s,
`stokes_sphere_sympy` ≈ 6 s, `oseen_*_sympy` ≈ 1 s — all `lru_cache`d), five animations and eight plotly figures. ch07 ran
in ≈ 190 s with 16 CORE blocks.

| Section | Heaviest cells | Full | FAST (`FLUIDPY_FAST=1`) |
|---|---|---|---|
| setup + imports | numpy/scipy/sympy/plotly/pint, three new core modules | 9 s | 9 s |
| §8.1 C01, R01–R06 | property functions, dye sketch (seeded random walk 2000 points), NS residual demo | 3 s | 2 s |
| §8.2 C02 | `parallel_flow_sympy("channel")`, FD solve 101², Fig. 8.4 remake, entrance sketch (brentq), IF1 (21 steps × 3 traces) | 5 s | 4 s |
| §8.2 C03 | `parallel_flow_sympy("pipe")`, trapezoid 2001, two-panel figure, IF2 (13 steps) | 4 s | 3 s |
| §8.2 C04 | `parallel_flow_sympy("circular_couette")`, pressure FD check 2001, power `quad`, three-panel figure, IF3 (21 steps) | 5 s | 4 s |
| §8.3 C05 | `lubrication_nondim_sympy` (★★★ check), coefficient figure | 4 s | 3 s |
| §8.3 C06 | `reynolds_equation_sympy` (★★★), `reynolds_pressure_1d` 401, `quad`, Fig. 8.8 remake, Hele-Shaw figure + Poisson 128² | 7 s | 4 s (64²) |
| §8.3 C07 | `slider_bearing_sympy` (★★★), brentq + quad from scratch, two-panel figure, IF4 (20 steps × 3 traces), live cell (not run on the page), `quad` of W | 7 s | 5 s (10 steps) |
| §8.3 C08 | `thin_film_spread` 400 cells (cached; 6 s first run), from-scratch 200 explicit steps, **A3 frames 30** | 12 s | 6 s (200 cells, 16 frames) |
| §8.4 C09 | CN by hand 400 × 200 + `crank_nicolson_1d` + FTCS, convergence table, `solve_bvp`, collapse check, three-panel figure, stopped/Couette figure (200-term series), **A1 video 60**, IF5 (25 frames × 2 subplots) | 18 s | 10 s (200 points, 30 frames, 12 IF5 frames) |
| §8.4 C10 | `similarity_reduce_sympy` × 5, collapse errors, log–log fits (reuses the cached film run), sympy residual cells, three-panel figure, **A5 video 60** | 12 s | 7 s (30 frames) |
| §8.5 C11 | `stokes_second_sympy`, CN 10 periods, Fig. 8.16 remake, **A2 video 60**, IF6 (25 frames) | 10 s | 6 s (5 periods, 30 frames) |
| §8.6 C12 | `low_re_scaling_sympy`, `stokes_residual` (3 × 19 stencils), coefficient figure, IF7 (25 steps) | 3 s | 2 s |
| §8.6 C13 | `stokes_sphere_sympy` (★★★, shared with C14), FD check 20 points, three streamline panels (contours 301², FAST 151²), **A4 video 60** (tracers precomputed once) | 14 s | 8 s |
| §8.6 C14 | settling, drag parts, midpoint 2000 + Gauss–Legendre, synthetic Millikan (40 drops), three-panel figure | 4 s | 3 s |
| §8.6 C15 | `inertia_viscous_ratio` 60 radii × 3 Re, Oseen contours 200² × 3 (FAST 100²), C_D figure, `oseen_limit_sympy`, IF8 (5 steps of precomputed polylines) | 8 s | 5 s |
| explainers (9 × `show_viz`) | read the HTML files | 1 s | 1 s |
| **Total** | | **≈ 126 s** | **≈ 82 s** |

**FAST plan.** Every size-dependent choice is written `a if not FAST else b` in the cell: thin-film cells 400 → 200 (the
cache key includes the grid, so FAST keeps its own cache); CN grids 400 → 200 points with 200 → 100 steps (still ≈ 20
points across δ₉₉); Stokes-second numerical twin 10 → 5 periods; Hele-Shaw grid 128² → 64²; Oseen contour grids 200² →
100²; streamline contours 301² → 151²; videos 60 → 30 frames, frame players 30 → 16; plotly sliders ≤ 25 → ≤ 12 steps (≤ 4
traces × ≤ 400 points, < 250 kB each). **Cached arrays / results:** every `*_sympy` engine is `functools.lru_cache`d (D10's
check and the E2 coefficient table share `lubrication_nondim_sympy`; D28–D31 share `stokes_sphere_sympy`); the thin-film run
is written once to `outputs/ch08/thin_film_run.npz` and reused by A3, C10's front fit and the C08 from-scratch comparison
(its final time); the C09 CN solution is reused by A1; the A4 tracer paths are integrated once for all frames (P137
pattern). **Outputs:** 4 videos (A1, A2, A4, A5) + 1 frames player (A3) at dpi 80 (each < 2.5 MB), 8 plotly figures, ≈ 28
static figures — the page stays under 15 MB. Sympy cells never call `simplify` on expressions with more than one generic
function (`reynolds_equation_sympy` builds u from a generic p_x(x, t) and h(x, t) and simplifies only the final residual).

---

## Part E — prerequisite ledger
Every concept, symbol, maths tool and Python function or idiom the notebook or its explainers use, with where it is
explained. "primer (in Cxx)" = a 📎 primer placed in that block before first use (the Concept text is the primer term,
used verbatim in `nb.primer`); "knowledge/primers.md: <term> (chNN Pnn) — reminder" = a one-line reminder naming the
earlier primer; a CORE/RECAP id alone = taught there; "Cxx (Nnn)" = the NOTE placed in that block; "Cxx (gloss …)" = one
sentence where it is used. Earlier chapters' material that is neither a primer nor a ch08 RECAP is named by section
("Ch. 4 §4.4"). New primers P185–P199 in first-use order: diffusivity and the diffusion time (C01), Laplacian in
cylindrical coordinates (C03), Euler–Cauchy ODE (C04), anisotropic scaling (C05), Leibniz rule with a moving limit and
cumulative_trapezoid (C06), nonlinear diffusion in flux form and implicit Picard stepping (C08), Crank–Nicolson with
solve_banded, the Gaussian integral, erf and its inverses, solve_bvp (C09), exponent matching (C10), dominant balance (C12),
the Stokes operator E² applied twice (C13) — 15 primers.

| Concept | First used in | Explained by |
|---|---|---|
| laminar vs turbulent flow, Reynolds's dye experiment | C01 | C01 (N103) |
| Reynolds number Re = Ud/ν (pipe, diameter, mean velocity) | C01 | C01 |
| inertia/viscous ratio = Re | C01 | R01 |
| core.similarity.reynolds_number, ch08.inertia_viscous_scales | C01 | R01 (code comment) |
| vorticity diffusion Dω_z/Dt = ν∇²ω_z (5.13) | C01 | R02 |
| dynamic vs kinematic viscosity μ, ν = μ/ρ | C01 | C01 (⚠️ callout); knowledge/primers.md: stress (ch01 P05) — reminder |
| diffusivity and the diffusion time L²/ν | C01 | primer (in C01) |
| heat equation (4.89), thermal diffusivity κ | C01 | C01 (N01) |
| momentum diffusivity ν of air and water | C01 | C01 (D01) |
| ratios of property values | C01 | C01 (D01 steps 3–4) |
| ch08.pipe_flow_regime, momentum_diffusivity, diffusion_time | C01 | C01 (code explain) |
| ch01.water_viscosity, sutherland_viscosity, water_density | C01 | C01 (code comment; Ch. 1 §1.5) |
| power laws and log–log plots | C01 | knowledge/primers.md: power laws and log–log plots (ch01 P13) — reminder |
| np.random.default_rng (dye sketch) | C01 | knowledge/primers.md: np.random.default_rng (ch01 P10) — reminder |
| matplotlib figures | C01 | knowledge/primers.md: matplotlib figures (ch01 P01) — reminder |
| Python dictionaries (fluidpy results) | C01 | knowledge/primers.md: Python dictionaries (ch01 P23) — reminder |
| f-strings | C01 | knowledge/primers.md: f-strings (ch01 P04) — reminder |
| tuple unpacking | C01 | knowledge/primers.md: tuple unpacking (ch01 P14) — reminder |
| assert np.allclose | C01 | knowledge/primers.md: assert np.allclose (ch01 P15) — reminder |
| road map to Ch. 9, 11, 12, 13 | C01 | C01 (N02) |
| Navier–Stokes (8.1) = (4.85), continuity (4.10) | C02 | R03 |
| NS residual, core.navier_stokes.ns_incompressible_terms | C02 | R03 (code explain) |
| functions as arguments and lambda | C02 | knowledge/primers.md: functions as arguments and lambda (ch01 P29) — reminder |
| component-first array layout (x[0], x[1]) | C02 | R03 (code comment); knowledge/primers.md: np.meshgrid and the project grid layout (ch02 P76) — reminder |
| no-through-flow (8.2) | C02 | R04 |
| no-slip (8.3) | C02 | R05 |
| unit normal and tangent vectors | C02 | R04 (gloss); knowledge/primers.md: orthonormal basis, projection and completeness (ch02 P65) — reminder |
| ch08.wall_bc_residuals | C02 | R05 (code explain) |
| gravity absorbed into pressure, constant ρ, inertial frame | C02 | R06 |
| plane Couette flow u = Uy/h | C02 | R07 |
| plane Poiseuille flow, linear τ | C02 | R08 |
| G = −dp/dx vs the book's dp/dx | C02 | R08 (⚠️ callout), C02 |
| ch04.plane_poiseuille | C02 | R08 (code comment) |
| fully developed flow, entrance length | C02 | C02 (N03) |
| continuity ⇒ v ≡ 0 | C02 | C02 (N04, D02) |
| partial derivative | C02 | knowledge/primers.md: partial derivative (ch01 P25) — reminder |
| reduced equations (8.4a), (8.4b) | C02 | C02 (N05, N06, D02) |
| 'function of x = function of y ⇒ constant' | C02 | C02 (N07, D03) |
| separation of variables | C02 | knowledge/primers.md: separation of variables (ch01 P42) — reminder |
| integrating an ODE twice and fixing two constants | C02 | C02 (gloss with sympy dsolve); knowledge/primers.md: linear second-order ODE (ch01 P44) — reminder |
| fundamental theorem of calculus | C02 | knowledge/primers.md: fundamental theorem of calculus (ch02 P84) — reminder |
| two linear equations in two unknowns, np.linalg.solve | C02 | knowledge/primers.md: np.linalg.solve (ch02 P57) — reminder |
| sympy and sympy.dsolve | C02 | knowledge/primers.md: sympy (ch01 P40) — reminder |
| Couette–Poiseuille profile (8.5) | C02 | C02 (D04) |
| book's constants A, B | C02 | C02 (N08) |
| superposition of solutions of a linear equation | C02 | C02 (idea cell) |
| shear stress τ = μ du/dy, signed | C02 | C02 (D05); knowledge/primers.md: stress (ch01 P05) — reminder |
| derivative as a slope | C02 | knowledge/primers.md: ordinary derivative as a slope (ch01 P19) — reminder |
| backflow threshold dp/dx > 2μU/h² | C02 | C02 (N09, D05) |
| inequalities under division by a positive number | C02 | knowledge/primers.md: inequalities under a sign change (ch01 P48) — reminder |
| flow rate Q and mean velocity V, the printed V slip | C02 | C02 (N10) |
| definite integral, scipy.integrate.quad | C02 | knowledge/primers.md: scipy.integrate.quad and dblquad (ch03 P87) — reminder |
| finite differences (tridiagonal solve from scratch) | C02 | knowledge/primers.md: finite differences (ch01 P21) — reminder |
| np.diag, np.r_ | C02 | C02 (gloss in the from-scratch cell) |
| linear stress in fully developed channel flow | C02 | C02 (N11) |
| slider_figure | C02 | knowledge/primers.md: slider_figure (ch01 P17) — reminder |
| show_viz and the explainer tabs | C02 | knowledge/primers.md: show_viz (ch01 P18) — reminder |
| ch08.channel_flow, channel_flow_rate, channel_shear_stress, channel_backflow_threshold, couette_poiseuille_state, parallel_flow_sympy | C02 | C02 (code explain) |
| cylindrical coordinates (R, φ, z) and unit vectors | C03 | knowledge/primers.md: cylindrical and spherical unit vectors (ch03 P88) — reminder |
| cylindrical shear stress τ_zR | C03 | R09 |
| core.curvilinear.strain_rate | C03 | R09 (code comment) |
| Laplacian in cylindrical coordinates | C03 | primer (in C03) |
| pipe set-up, p = p(z) | C03 | C03 (N12, D06) |
| natural logarithm and ln R → −∞ | C03 | knowledge/primers.md: natural logarithm and exponential (ch01 P36) — reminder |
| regularity (boundedness) at the axis | C03 | C03 (gloss, N13, D06 step 8) |
| Poiseuille pipe profile (8.6) | C03 | C03 (D06) |
| stress (8.7), wall stress (8.8) | C03 | C03 (N14, N15, D07) |
| control-volume force balance on a slug | C03 | C03 (D07 steps 4–5; Ch. 4 §4.4) |
| area element 2πR dR of a disc | C03 | C03 (gloss); knowledge/primers.md: volume and surface integrals as midpoint sums (ch02 P83) — reminder |
| Hagen–Poiseuille Q, V, u_max = 2V | C03 | C03 (N16, D07) |
| Darcy friction factor f = 64/Re | C03 | C03 (D07 steps 9–10) |
| trapezoid rule (from scratch) | C03 | knowledge/primers.md: trapezoid rule (ch01 P37) — reminder |
| ch03.pipe_profile, core.navier_stokes.exact_solution | C03 | C03 (code comment; Ch. 3 §3.2, Ch. 4 §4.6) |
| ch08.pipe_poiseuille, pipe_shear_stress, pipe_wall_stress, pipe_flow_rate, pipe_friction_factor | C03 | C03 (code explain) |
| ideal vortex (5.2), zero net viscous force | C04 | R10 |
| curl of a curl identity (∇²u = −∇×ω) | C04 | R10; knowledge/primers.md: curl of a curl identity (ch04 P122) — reminder |
| power in = dissipation, its sign | C04 | R11 |
| ch05.dissipation_outside_cylinder | C04 | R11 (code comment) |
| solid-body rotation (5.1) | C04 | R12 |
| core.vortices.solid_body_rotation, ch05.rotating_cylinder_flow (ω = 2Ω) | C04 | R12, C04 (N21 ⚠️) |
| centripetal acceleration in cylindrical coordinates | C04 | C04 (N17, D08 step 2); knowledge/primers.md: polar coordinates as a moving basis (ch03 P105) — reminder |
| R- and φ-momentum of circular Couette flow | C04 | C04 (N17, D08) |
| Euler–Cauchy (equidimensional) ODE | C04 | primer (in C04) |
| general solution AR + B/R (8.9) | C04 | C04 (N18, D08) |
| constants A, B from the walls | C04 | C04 (N19, D08) |
| circular Couette profile (8.10) | C04 | C04 (D08) |
| limits (8.11), (8.12) and Γ = 2πΩ₁R₁² | C04 | C04 (D09, N21) |
| limits and orders of smallness | C04 | knowledge/primers.md: limits and orders of smallness (ch02 P68) — reminder |
| circulation of a vortex | C04 | C04 (D09 step 4; Ch. 5 §5.1) |
| pressure from the centripetal balance (ours) | C04 | C04 (N20) |
| np.gradient | C04 | knowledge/primers.md: np.gradient (ch01 P22) — reminder |
| advective acceleration of the three confined flows | C04 | C04 (N22) |
| free-surface dip of a vortex | C04 | C04 (N104; Ch. 5 §5.1) |
| Rayleigh criterion (named) | C04 | C04 (figure note; Ch. 11) |
| ch08.circular_couette, circular_couette_pressure, circular_couette_shear_stress, circular_couette_power, circular_couette_state, advective_acceleration_check | C04 | C04 (code explain) |
| gap field equations (6.2), (8.13a) | C05 | R13 |
| lubrication idea h ≪ L | C05 | C05 (N23) |
| y-momentum (8.13b) and the printed ∂p/∂x slip | C05 | C05 (N24) |
| anisotropic scaling with two length scales | C05 | primer (in C05) |
| scaled variables and the chain rule | C05 | knowledge/primers.md: scaled variables and the chain rule (ch04 P133) — reminder |
| bookkeeping of a small parameter ε | C05 | C05 (gloss); knowledge/primers.md: limits and orders of smallness (ch02 P68) — reminder |
| scalings (8.14), fineness ratio ε | C05 | C05 (N26) |
| scaled continuity (8.15) | C05 | C05 (N27, D10) |
| scaled momentum (8.16a), (8.16b); Re_L, Λ | C05 | C05 (N28, N29, D10) |
| lubrication balance (8.17a), (8.17b), the printed missing ν | C05 | C05 (D11, N30) |
| order-of-magnitude scaling | C05 | knowledge/primers.md: order-of-magnitude scaling (ch04 P130) — reminder |
| bearing number Λ and the natural pressure scale μUL/h² | C05 | C05 (N31) |
| sympy symbols and simplify in check cells | C05 | knowledge/primers.md: sympy expand, series, removeO, collect and subs (ch04 P117) — reminder |
| ch08.lubrication_scales, lubrication_term_magnitudes, lubrication_nondim_sympy | C05 | C05 (code explain) |
| gap boundary conditions, time as a parameter | C06 | C06 (N25) |
| integration with x-dependent 'constants' A(x, t), B(x, t) | C06 | C06 (D12 step 3) |
| generic profile (8.18) | C06 | C06 (N33, D12) |
| lubrication profile (8.19) and the U₀ slip | C06 | C06 (D12, N34) |
| Leibniz rule with a moving upper limit | C06 | primer (in C06) |
| differentiation under the integral sign | C06 | knowledge/primers.md: differentiation under the integral sign (ch03 P109) — reminder |
| kinematic condition at a moving wall | C06 | C06 (gloss); knowledge/primers.md: moving level set and its normal speed (ch04 P132) — reminder |
| gap flux q and the 1-D Reynolds equation | C06 | C06 (N32, D13) |
| scipy.integrate.cumulative_trapezoid | C06 | primer (in C06) |
| Hele-Shaw flow, depth-averaged velocity, ∇²p = 0 | C06 | C06 (N36) |
| 2-D potential flow (6.10), (6.12) | C06 | C06 (N36; Ch. 6 §6.2) |
| core.potential.Flow cylinder, core.laplace_solvers.solve_poisson | C06 | C06 (N36 code comment; Ch. 6 §6.3, §6.7) |
| ch08.lubrication_velocity, lubrication_flux, reynolds_pressure_1d, slider_gap_velocity, hele_shaw_velocity, hele_shaw_potential, hele_shaw_cylinder, reynolds_equation_sympy | C06 | C06 (code explain) |
| moving control volume and (4.5) | C07 | C07 (D14 steps 1–2; Ch. 4 §4.2) |
| antiderivatives of (1 + αx/L)⁻ⁿ | C07 | C07 (gloss); knowledge/primers.md: substitution in an integral (ch03 P106) — reminder |
| Taylor expansion in a small parameter α | C07 | C07 (gloss); knowledge/primers.md: first-order Taylor expansion (ch01 P26) — reminder |
| slider-bearing pressure and load, the two printed slips | C07 | C07 (D14) |
| exact load for any taper, optimum 1 + α = 2.19 | C07 | C07 (N35) |
| scipy.optimize.brentq | C07 | knowledge/primers.md: scipy.optimize.brentq (ch03 P108) — reminder |
| scipy.optimize.minimize_scalar | C07 | knowledge/primers.md: `scipy.optimize.minimize_scalar` (ch07 P170) — reminder |
| np.log1p and cancellation for small α | C07 | knowledge/primers.md: np.expm1 and cancellation near zero (ch03 P107) — reminder |
| live widgets | C07 | knowledge/primers.md: live widgets (ch01 P47) — reminder |
| inlet recirculation for α > 1 (ours) | C07 | C07 (E3 Explain §5; D14 check) |
| ch08.slider_bearing, slider_bearing_load, slider_optimum_taper, slider_bearing_state, slider_bearing_sympy | C07 | C07 (code explain) |
| hydrostatic pressure in a thin layer, stress-free surface | C08 | C08 (gloss); knowledge/primers.md: boundary conditions (ch01 P20) — reminder |
| half-parabola film profile and its flux | C08 | C08 (D15) |
| thin-film equation (Example 8.3) | C08 | C08 (D15) |
| nonlinear diffusion in flux form | C08 | primer (in C08) |
| implicit time stepping with Picard iteration | C08 | primer (in C08) |
| explicit stepping and its stability limit | C08 | knowledge/primers.md: explicit stepping (ch01 P30) — reminder |
| np.where | C08 | knowledge/primers.md: np.where and np.select (ch01 P46) — reminder |
| animate and show_animation | C08 | knowledge/primers.md: animate and show_animation (ch01 P16) — reminder |
| similarity solution (named) | C08 | C08 (N37), C09 |
| ch08.thin_film_flux, thin_film_spread, viscous_current_similarity, thin_film_state | C08 | C08 (code explain) |
| Couette start-up series, imposed length | C09 | R14 |
| core.diffusion.couette_startup_profile | C09 | R14 (code comment) |
| Stokes' first problem set-up (Rayleigh) | C09 | C09 (N38, D16) |
| 1-D diffusion equation (8.20) | C09 | C09 (N39, D16) |
| conditions (8.21)–(8.23), order vs number of conditions | C09 | C09 (N40–N42, gloss) |
| Buckingham Π and ch01.pi_groups | C09 | C09 (N43, D17); Ch. 1 §1.11 |
| linearity of a solution (u ∝ U) | C09 | C09 (D17 steps 5–6) |
| similarity variable η = y/√(νt) (8.25) | C09 | C09 (N44, D17) |
| chain rule | C09 | knowledge/primers.md: chain rule (ch01 P49) — reminder |
| similarity ODE (8.26) and conditions (8.27), (8.28) | C09 | C09 (N45–N48, D18) |
| exponent rules | C09 | knowledge/primers.md: exponent rules (ch01 P43) — reminder |
| natural logarithm, exponentiating | C09 | knowledge/primers.md: natural logarithm and exponential (ch01 P36) — reminder |
| Gaussian integral | C09 | primer (in C09) |
| substitution ξ = 2ζ in an integral | C09 | knowledge/primers.md: substitution in an integral (ch03 P106) — reminder |
| error function erf and its inverses | C09 | primer (in C09) |
| complementary error function erfc | C09 | knowledge/primers.md: complementary error function erfc (ch04 P123) — reminder |
| solution (8.29), (8.30) | C09 | C09 (N49–N51, D19) |
| scipy.integrate.solve_bvp | C09 | primer (in C09) |
| 99 % thickness (8.31) | C09 | C09 (D20, N54) |
| self-similarity and the collapse | C09 | C09 (N52, N55) |
| vorticity content ∫ω dy = U, the printed −U | C09 | C09 (N53) |
| similarity destroyed by an imposed scale; stopped plate | C09 | C09 (N56) |
| Crank–Nicolson with scipy.linalg.solve_banded | C09 | primer (in C09) |
| observed order of convergence | C09 | C09 (code explain; tools.convergence) ; knowledge/primers.md: power laws and log–log plots (ch01 P13) — reminder |
| animate_figure (plotly time slider) | C09 | C09 (code explain); knowledge/primers.md: slider_figure (ch01 P17) — reminder |
| ch08.stokes_first_problem, similarity_variable, diffusion_thickness, stokes_first_vorticity, stokes_first_stopped, stokes_first_state, vorticity_content, crank_nicolson_1d, similarity_ode_solve | C09 | C09 (code explain) |
| similarity ansatz (8.32a), (8.32b) | C10 | C10 (N57) |
| exponent matching for similarity forms | C10 | primer (in C10) |
| Example 8.4, δ = √(2C₁νt) | C10 | C10 (N58, D21) |
| separable first-order ODE for δ(t) | C10 | C10 (D21 step 7); knowledge/primers.md: separation of variables (ch01 P42) — reminder |
| product rule read backwards, (ηF)′ | C10 | C10 (gloss); knowledge/primers.md: product rule for differentials (ch01 P38) — reminder |
| integral constraint (conserved jump) | C10 | C10 (D22 steps 7–9) |
| viscous vortex sheet (Example 8.5), width 5.544√(νt), printed 2.76 | C10 | C10 (N59, D22) |
| ch05.diffusing_vortex_sheet (γ = −2U) | C10 | C10 (N59 ⚠️) |
| Galilean shift, temporally developing boundary layer, C_f | C10 | C10 (N60); knowledge/primers.md: frames of reference and relative velocity (ch03 P96) — reminder |
| line-vortex decay (Example 8.6), Gaussian vortex (3.29) | C10 | C10 (N61; Ch. 3 §3.5) |
| line vortex switched on (Exercise 8.26) | C10 | C10 (N62) |
| core.vortices.gaussian_vortex | C10 | C10 (N61 code comment) |
| spreading-bead exponents (Example 8.7), Huppert dome | C10 | C10 (N63, D23) |
| diffusive vs advective length scales | C10 | C10 (N64) |
| ch08.similarity_reduce_sympy, similarity_collapse_error, vortex_sheet_diffusion, transition_width, temporal_bl_wall_stress, line_vortex_decay, line_vortex_spinup | C10 | C10 (code explain) |
| oscillating plate set-up, imposed time scale | C11 | C11 (N65) |
| conditions (8.33), (8.34) | C11 | C11 (N66, N67) |
| ω as a frequency again (not vorticity) | C11 | C11 (⚠️ callout) |
| complex exponential solution, taking the real part | C11 | C11 (gloss); knowledge/primers.md: complex amplitudes (ch07 P176) — reminder |
| square root of i | C11 | C11 (gloss); knowledge/primers.md: complex square roots and the quadratic formula (ch06 P159) — reminder |
| Euler's formula, i² = −1 | C11 | knowledge/primers.md: square root of a negative number (ch01 P45) — reminder |
| complex trial (8.35), ODE (8.36), roots k, (8.37) | C11 | C11 (N68–N71, D24) |
| Stokes layer (8.38) | C11 | C11 (D24) |
| e-folding depth vs the book's 4√(ν/ω), crest speed | C11 | C11 (N72, D25) |
| phase of a cosine and a time lag | C11 | knowledge/primers.md: phase of a wave (ch07 P165) — reminder |
| three groups, no similarity | C11 | C11 (N73) |
| acoustic boundary layer, Ekman link (named) | C11 | C11 (N74) |
| ch08.stokes_second_problem, stokes_layer, stokes_layer_state, stokes_second_sympy | C11 | C11 (code explain) |
| steady NS (8.39) | C12 | R15 |
| high-Re scaling (4.100), (8.40) | C12 | R16 |
| core.similarity.nondimensional_ns_coefficients | C12 | R16 (code comment) |
| perturbation expansions (named) | C12 | C12 (N75) |
| high-Re non-uniformity at walls (named) | C12 | C12 (N76) |
| dominant balance | C12 | primer (in C12) |
| (8.41), the viscous pressure scale, (8.42) | C12 | C12 (N77–N79, D26) |
| Stokes equations (8.43), creeping flow | C12 | C12 (D26) |
| scaling rule by region (named) | C12 | C12 (N80) |
| ch08.low_re_scaling_sympy, stokes_residual | C12 | C12 (code explain) |
| azimuthal vorticity ω_φ in spherical coordinates | C13 | R17 |
| Stokes stream function (6.83) | C13 | R18 |
| uniform stream (6.86), condition (8.47) | C13 | R19 |
| spherical coordinates with θ from the downstream axis | C13 | C13 (⚠️ in the front-matter conventions table); knowledge/primers.md: cylindrical and spherical unit vectors (ch03 P88) — reminder |
| curl of a gradient is zero, curl commutes with the Cartesian Laplacian | C13 | C13 (gloss); knowledge/primers.md: Schwarz's theorem (ch04 P121) — reminder |
| index notation ε_ijk | C13 | C13 (D27; Ch. 2 §2.7); knowledge/primers.md: permutations, cyclic order and parity (ch02 P72) — reminder |
| ∇²ω = 0 and its spherical meaning | C13 | C13 (N81, D27) |
| the Stokes operator E² applied twice | C13 | primer (in C13) |
| vorticity from ψ, ω_φ = −E²ψ/(r sin θ) | C13 | C13 (N82, D28) |
| stream-function equation (8.44), not the biharmonic | C13 | C13 (N83, D28) |
| sphere conditions (8.45), (8.46) | C13 | C13 (N84, N85) |
| separable ψ = f(r) sin²θ and the f-ODE | C13 | C13 (N86, D29) |
| Stokes' stream function (8.48) | C13 | C13 (D29) |
| velocities (8.49) | C13 | C13 (N87) |
| fluid-frame ψ, reversibility, no wake | C13 | C13 (N93) |
| ideal-flow sphere (Ch. 6) | C13 | C13 (figure; Ch. 6 §6.8) |
| meshgrid, contour and streamplot | C13 | knowledge/primers.md: plt.contour, plt.quiver and plt.streamplot (ch02 P78) — reminder |
| advecting many tracers in one solve_ivp call | C13 | knowledge/primers.md: advecting many points in one `solve_ivp` call (ch05 P137) — reminder |
| ch08.stokes_sphere_streamfunction, stokes_sphere_velocity, stokes_sphere_velocity_xyz, side_line_speed, E2, E4_residual, stokes_sphere_sympy | C13 | C13 (code explain) |
| ρ-free dimensional argument D = f(μ, U, a) | C14 | R20 |
| line integral of a gradient | C14 | C14 (gloss); knowledge/primers.md: line integral along a path (ch01 P35) — reminder |
| pressure (8.50) and the printed rear-minimum sign | C14 | C14 (N88, D30) |
| Newtonian stress in spherical coordinates | C14 | C14 (D31 steps 3–7; Ch. 4 §4.5) |
| traction on a surface and its x-projection | C14 | C14 (gloss, D31 step 2; Ch. 2 §2.6) |
| surface integrals on a sphere, dA = 2πa² sin θ dθ | C14 | knowledge/primers.md: surface integrals on a sphere (ch06 P163) — reminder |
| surface stresses σ_rr, σ_rθ (ours) | C14 | C14 (N89) |
| Stokes drag (8.51), ⅓ pressure + ⅔ friction | C14 | C14 (D31) |
| Gauss–Legendre quadrature | C14 | knowledge/primers.md: Gauss–Legendre quadrature in 3-D and a smoothed kernel (ch05 P143) — reminder |
| terminal velocity, effective weight, buoyancy | C14 | C14 (N90); knowledge/primers.md: weight and gravitational acceleration (ch01 P24) — reminder |
| Millikan's experiment, charge quantum | C14 | C14 (N91) |
| drag coefficient (8.52), Re with the diameter | C14 | C14 (N92) |
| core.similarity.sphere_drag_coefficient | C14 | C14 (code comment; Ch. 4 §4.11) |
| ch08.stokes_sphere_pressure, stokes_sphere_surface_stresses, stokes_drag, stokes_drag_running, sphere_drag_quadrature, settling_state, terminal_velocity, radius_from_terminal_velocity, millikan_charge, synthetic_millikan, stokes_drag_coefficient | C14 | C14 (code explain) |
| asymptotic size of a term far away | C15 | C15 (gloss, D32); knowledge/primers.md: order-of-magnitude scaling (ch04 P130) — reminder |
| inertia/viscous ~ Re r/a, far-field breakdown | C15 | C15 (D32) |
| singular perturbation at infinity, Stokes' and Whitehead's paradoxes (named) | C15 | C15 (N94) |
| Oseen's linearisation and the printed sign slip | C15 | C15 (N95) |
| Oseen boundary conditions | C15 | C15 (N96) |
| Oseen stream function (8.53) | C15 | C15 (N97) |
| series of 1 − e^(−s) | C15 | C15 (gloss); knowledge/primers.md: sympy expand, series, removeO, collect and subs (ch04 P117) — reminder |
| Oseen → Stokes near the sphere | C15 | C15 (N98, D33) |
| Oseen drag coefficient, radius vs diameter Re | C15 | C15 (N99) |
| Oseen wake (Fig. 8.20) | C15 | C15 (N100) |
| matched asymptotics, Proudman–Pearson (named) | C15 | C15 (N101) |
| Morrison drag correlation (benchmark band) | C15 | C15 (figure note; Ch. 4 §4.11) |
| ch08.inertia_viscous_ratio, oseen_streamfunction, oseen_velocity, oseen_drag_coefficient, proudman_pearson_drag_coefficient, oseen_linearisation_sympy, oseen_limit_sympy | C15 | C15 (code explain) |
| perturbation methods and CFD (final remarks) | C15 | C15 (N102, in §8.7) |

---

## Part F — derivation storyboards

Builders copy these word for word into `nb.derivation(key, title, goal=…, start=(tex, plain), plan=[…], uses=[…],
steps=[dict(did, tex, why, plain)], result=(tex, plain), interpret=…, check=…, check_src=…)` and into the explainer's
`derivations: [...]` (phones may shorten *why* to its first sentence; `live`, `set` and `watch` are the explainer's and are
listed in Part B). Every step is one move; *why* names the rule and says why we make it (≤ 35 words, ≥ 6); *did* ≤ 8 words.
The book's own moves were read on the rendered pages (p338 for D01; p340 for D02–D04; p341 for D05; p342–p343 for D06–D07;
p343–p345 for D08–D09; p346–p347 for D10–D12; D13 is not in the text (Exercises 8.19–8.20); p348–p349 for D14; p351–p353 for
D15; p353–p354 for D16–D17; p355–p356 for D18–D19; p357 for D20; p357–p358 for D21; p358–p360 for D22; p363 for D23; p364 for
D24; p365 for D25; p366–p367 for D26; p368 for D27–D28; p369 for D29–D31 (D30, D31 not written out by the book); p371 for
D32; p373 for D33); the gaps listed in `analysis/ch08.md` §2b are filled and the notebook says so ("the book skips this
move"). Every equation named by number is written out. Colours: pressure orange, viscous rose, inertia teal, velocity
blue, vorticity purple. No line of this part starts with a table bar; absolute values are written with \lvert \rvert or in
words.

### D01 · Air diffuses momentum ≈ 15× faster than water: ν = μ/ρ and t ~ L²/ν — ★, 5 steps, in C01 (notebook)
- **Goal.** Show that the rate at which motion spreads is set by ν = μ/ρ, not μ, and that air — 55 times less viscous than
  water — spreads momentum about 15 times faster.
- **Start.** $DT/Dt=\kappa\nabla^2T$ (4.89) and $D\mathbf u/Dt=-(1/\rho)\nabla p+\nu\nabla^2\mathbf u$ (8.1) — *in words:* in
  both, a diffusivity multiplies the Laplacian.
- **Plan.** (1) Identify the diffusivity in (8.1). (2) Turn a diffusion balance into a time scale. (3) Compare air and water
  with numbers.
- **Tools.** diffusivity and the diffusion time (primer P185) · order-of-magnitude scaling (P130) · ratios.
- **Assumptions.** Room temperature and pressure property values (step 4).
- **Steps.**
  1. *did:* Read the Laplacian's coefficient · *tex:* $\nu\equiv\frac{\mu}{\rho}\ \ [\mathrm{m^2/s}]$ · *why:* (8.1) is the
     momentum equation divided by ρ; the viscous term μ∇²u becomes (μ/ρ)∇²u, with the units m²/s of a diffusivity like κ.
     · *plain:* the diffusivity of momentum is ν, not μ.
  2. *did:* Balance change against diffusion · *tex:* $\frac{U}{t_d}\sim\nu\frac{U}{L^2}\ \Rightarrow\ t_d\sim\frac{L^2}{\nu}$ ·
     *why:* Order-of-magnitude scaling: ∂u/∂t ~ U/t and ∂²u/∂y² ~ U/L²; equating them gives the time diffusion needs to
     cross a distance L. · *plain:* spreading over L takes about L²/ν.
  3. *did:* Write the ratio of the two ν · *tex:* $\frac{\nu_{air}}{\nu_{w}}=\frac{\mu_{air}}{\mu_{w}}\cdot\frac{\rho_{w}}{\rho_{air}}$
     · *why:* Divide ν = μ/ρ for air by the same for water; the ratio splits into a viscosity ratio and an inverse density
     ratio. · *plain:* two effects compete: viscosity and density.
  4. *did:* Insert room-temperature values · *tex:* $\frac{\nu_{air}}{\nu_w}=\frac{1.8\times10^{-5}}{10^{-3}}\cdot\frac{1000}{1.2}=
     0.018\times833=15$ · *why:* μ and ρ at 20 °C (Ch. 1 property functions); air's density is 833 times smaller, which beats
     its 55 times smaller viscosity. · *plain:* air's ν is 15 times water's.
  5. *did:* Compare diffusion times across 1 cm · *tex:* $t_w=\frac{10^{-4}}{10^{-6}}=100\ \mathrm s,\ \ t_{air}=\frac{10^{-4}}{1.5\times
     10^{-5}}=6.7\ \mathrm s$ · *why:* Use step 2 with L = 1 cm; the ratio of times is the inverse ratio of ν. · *plain:*
     momentum crosses a centimetre of air 15 times sooner.
- **Result.** $\nu=\mu/\rho$, $t_d\sim L^2/\nu$; $\nu_{air}/\nu_w\approx15$ — *in words:* ν sets how fast motion spreads.
- **Check.** Units: (Pa s)/(kg/m³) = (kg/(m s))/(kg/m³) = m²/s ✓; t = m²/(m²/s) = s ✓. The property functions give 15.1.
- **What it means.** "More viscous" and "more diffusive" are different: honey is both, air is little viscous yet quickly
  diffusive. Every thickness in this chapter grows like √(νt).
- **Traps.** Comparing μ instead of ν; forgetting that the density ratio (≈ 830) is larger than the viscosity ratio (≈ 55).

### D02 · v ≡ 0 and the reduced equations $0=-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{d^2u}{dy^2}$ (8.4a), $0=-\frac1\rho\frac{\partial p}{\partial y}$ (8.4b) — ★, 7 steps, in C02 (notebook · `couette_poiseuille_backflow`)
- **Goal.** Show that in steady, fully developed flow between plates the vertical velocity vanishes and Navier–Stokes
  shrinks to two short linear equations.
- **Start.** $\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0$ (4.10 in 2-D) and $D\mathbf u/Dt=-(1/\rho)\nabla p+
  \nu\nabla^2\mathbf u$ (8.1) — *in words:* mass and momentum for the fluid between the plates.
- **Plan.** (1) Use "fully developed" in continuity. (2) Fix v with the walls. (3) Drop every term that vanishes.
- **Tools.** partial derivative (P25) · continuity (R03) · no through-flow (8.2) (R04).
- **Assumptions.** Steady (step 5); fully developed ∂u/∂x = 0 (step 1); 2-D, w = 0, ∂/∂z = 0 (step 1); impermeable walls
  (step 3).
- **Steps.**
  1. *did:* Use the fully developed condition · *tex:* $\frac{\partial u}{\partial x}=0\ \Rightarrow\ u=u(y)$ · *why:* Far past
     the entrance the profile no longer changes along the channel; with w = 0 and ∂/∂z = 0 the flow is 2-D. · *plain:* u
     depends on height only.
  2. *did:* Put it into continuity · *tex:* $\frac{\partial v}{\partial y}=0$ · *why:* (4.10) with ∂u/∂x = 0 leaves only the
     v-term; we want to learn what this says about v. · *plain:* v does not change with height.
  3. *did:* Use the impermeable walls · *tex:* $v(x,0)=0\ \Rightarrow\ v\equiv0$ · *why:* A function of y with zero y-derivative is
     the same at every height, and (8.2) makes it zero at the wall — so it is zero everywhere. · *plain:* the flow is
     parallel to the walls: u = (u(y), 0, 0).
  4. *did:* Evaluate the advective acceleration · *tex:* $u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=u\cdot0+0\cdot
     \frac{du}{dy}=0$ · *why:* Both factors vanish by steps 1 and 3 — exactly, not approximately; this is why the problem
     becomes linear. · *plain:* no fluid particle accelerates.
  5. *did:* Write the steady x-momentum · *tex:* $0=-\frac1\rho\frac{\partial p}{\partial x}+\nu\Big(\frac{\partial^2u}{\partial x^2}+
     \frac{\partial^2u}{\partial y^2}\Big)$ · *why:* The x-component of (8.1) with ∂u/∂t = 0 (steady) and step 4. · *plain:*
     pressure and friction balance.
  6. *did:* Drop the x-curvature · *tex:* $0=-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{d^2u}{dy^2}$ (8.4a) · *why:* u depends
     on y only, so ∂²u/∂x² = 0 and the partial y-derivative is an ordinary one. · *plain:* (8.4a): the pressure gradient is
     balanced by cross-channel friction.
  7. *did:* Write the y-momentum · *tex:* $0=-\frac1\rho\frac{\partial p}{\partial y}$ (8.4b) · *why:* Every v-term in the
     y-component of (8.1) is zero because v ≡ 0; only the pressure term survives. · *plain:* (8.4b): pressure does not
     change across the channel.
- **Result.** (8.4a) and (8.4b) — *in words:* a steady parallel flow is a pure pressure–friction balance.
- **Check.** Units of (8.4a): (1/ρ)∂p/∂x in (m³/kg)(Pa/m) = m/s², ν d²u/dy² in (m²/s)(1/(m s)) = m/s² ✓. The same steps
  applied in `ch08.parallel_flow_sympy("channel")` return these two equations.
- **What it means.** Because the nonlinear term vanishes identically, the solution will be exact, not approximate. In the
  entrance region (∂u/∂x ≠ 0) this fails: v ≠ 0 (N03).
- **Traps.** ∂v/∂y = 0 alone only says v is constant across the gap — the wall value fixes it to zero; the advective term is
  zero exactly, not "small".

### D03 · The pressure gradient is constant: $\frac1\mu\frac{dp}{dx}=\frac{d^2u}{dy^2}=\text{const}$ — ★, 5 steps, in C02 (notebook · `couette_poiseuille_backflow`)
- **Goal.** Show that in fully developed channel flow the pressure falls linearly along the channel.
- **Start.** $0=-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{d^2u}{dy^2}$ (8.4a) and $0=-\frac1\rho\frac{\partial p}{\partial y}$
  (8.4b) — *in words:* the reduced equations of D02.
- **Plan.** (1) Learn from (8.4b) what p depends on. (2) Separate (8.4a) into an x-side and a y-side. (3) Differentiate
  each side by the other variable.
- **Tools.** partial derivative (P25) · separation of variables (P42).
- **Assumptions.** Steady (p does not depend on t).
- **Steps.**
  1. *did:* Read (8.4b) · *tex:* $p=p(x)\ \Rightarrow\ \frac{\partial p}{\partial x}=\frac{dp}{dx}(x)$ · *why:* ∂p/∂y = 0 means p
     is the same at every height, so its x-derivative depends on x only. · *plain:* the pressure gradient can vary only
     along the channel.
  2. *did:* Multiply (8.4a) by ρ and rearrange · *tex:* $\frac{dp}{dx}(x)=\mu\frac{d^2u}{dy^2}(y)$ · *why:* ρν = μ; we isolate
     the pressure side, a function of x alone, and the velocity side, a function of y alone. · *plain:* an x-function equals
     a y-function at every point.
  3. *did:* Differentiate both sides by x · *tex:* $\frac{d^2p}{dx^2}=0$ · *why:* The right side does not contain x, so its
     x-derivative is zero; hence so is the left side's. · *plain:* the pressure gradient does not change along x.
  4. *did:* Differentiate both sides by y · *tex:* $\mu\frac{d^3u}{dy^3}=0$ · *why:* The left side does not contain y; so the
     y-derivative of the right side vanishes too — it is also a constant. · *plain:* the curvature of the profile is the same
     at every height.
  5. *did:* Name the common constant · *tex:* $\frac1\mu\frac{dp}{dx}=\frac{d^2u}{dy^2}=\text{const},\ \ p=p_0+\frac{dp}{dx}x$ ·
     *why:* Both sides are constants and equal each other; integrating d²p/dx² = 0 once more gives a linear pressure. ·
     *plain:* pressure falls (or rises) linearly, the profile has constant curvature.
- **Result.** $\frac1\mu\frac{dp}{dx}=\frac{d^2u}{dy^2}=\text{const}$ — *in words:* one constant sets both the pressure slope and
  the profile's curvature.
- **Check.** A pressure p = x² would make the left side 2x/μ, which cannot equal a function of y at every x —
  `parallel_flow_sympy` confirms no solution exists for non-constant dp/dx.
- **What it means.** The same "x-function = y-function" argument will fix the pressure gradient in the pipe (D06), between
  rotating cylinders (no pressure in the φ-equation, D08) and across the lubrication gap (D12).
- **Traps.** Differentiate each side by the *other* variable; the constant can have either sign (favourable or adverse).

### D04 · Couette–Poiseuille profile $u(y)=\frac Uhy-\frac1{2\mu}\frac{dp}{dx}y(h-y)$ (8.5) — ★★, 9 steps, in C02 (notebook · `couette_poiseuille_backflow`)
- **Goal.** Find the velocity at every height between a fixed floor and a top plate moving at U, with a constant pressure
  gradient.
- **Start.** $\mu\frac{d^2u}{dy^2}=\frac{dp}{dx}=\text{const}$ (D03), with u(0) = 0 and u(h) = U — *in words:* a constant
  curvature and two wall speeds.
- **Plan.** (1) Integrate twice. (2) Apply the two no-slip conditions. (3) Regroup into the book's form.
- **Tools.** integrating an ODE twice (gloss in C02) · 2 × 2 linear system (P57) · no slip (8.3) (R05).
- **Assumptions.** dp/dx constant (D03); no slip at both plates (steps 4–5).
- **Steps.**
  1. *did:* Integrate once in y · *tex:* $\mu\frac{du}{dy}=\frac{dp}{dx}y+c_1$ · *why:* Fundamental theorem of calculus; dp/dx
     is a constant (D03), so its antiderivative is (dp/dx)y plus a constant c₁. · *plain:* the shear stress changes linearly
     with height.
  2. *did:* Integrate again · *tex:* $\mu u=\frac{dp}{dx}\frac{y^2}{2}+c_1y+c_2$ · *why:* Integrate each term once more; a
     second constant c₂ appears. · *plain:* the profile is a parabola plus a straight line.
  3. *did:* Compare with the book's constants · *tex:* $0=-\frac{y^2}{2}\frac{dp}{dx}+\mu u+Ay+B,\ \ A=-c_1,\ B=-c_2$ · *why:*
     Moving everything to one side gives the book's line; its A and B are minus our constants — worth saying to avoid sign
     confusion. · *plain:* same equation, constants named differently.
  4. *did:* Apply no slip at the floor · *tex:* $u(0)=0\ \Rightarrow\ c_2=0$ · *why:* (8.3) on the fixed plate; at y = 0 the first
     two terms vanish, leaving c₂ = 0. · *plain:* the fluid at the floor is at rest.
  5. *did:* Apply no slip at the top · *tex:* $\mu U=\frac{dp}{dx}\frac{h^2}{2}+c_1h$ · *why:* (8.3) on the moving plate: at
     y = h the fluid moves at U. · *plain:* one equation for c₁.
  6. *did:* Solve for c₁ · *tex:* $c_1=\frac{\mu U}{h}-\frac h2\frac{dp}{dx}$ · *why:* Divide step 5 by h and subtract; the book's
     A = −c₁ = (h/2)dp/dx − μU/h matches the page. · *plain:* the wall stress constant mixes drag and pressure.
  7. *did:* Substitute the constants · *tex:* $\mu u=\frac{dp}{dx}\frac{y^2}{2}+\Big(\frac{\mu U}{h}-\frac h2\frac{dp}{dx}\Big)y$ ·
     *why:* Insert c₁ and c₂ = 0 into step 2 to get the profile. · *plain:* the profile with both constants fixed.
  8. *did:* Divide by μ and group · *tex:* $u=\frac Uhy+\frac{1}{2\mu}\frac{dp}{dx}(y^2-hy)$ · *why:* Collect the two dp/dx
     terms: y²/2 − hy/2 = (y² − hy)/2. · *plain:* a line plus a parabola.
  9. *did:* Rewrite the parabola · *tex:* $u(y)=\frac Uhy-\frac1{2\mu}\frac{dp}{dx}y(h-y)$ (8.5) · *why:* y² − hy = −y(h − y);
     this form makes the parabola positive inside the gap when dp/dx < 0. · *plain:* (8.5).
- **Result.** $u(y)=\frac Uhy-\frac1{2\mu}\frac{dp}{dx}y(h-y)$ (8.5) — *in words:* wall drag (a line) plus pressure push (a
  parabola).
- **Check.** u(0) = 0 ✓, u(h) = U ✓; units of (1/μ)(dp/dx)y² : (1/(Pa s))(Pa/m)m² = m/s ✓; dp/dx = 0 gives Couette,
  U = 0 gives Poiseuille (D05). Numbers: h = 1 cm, U = 0.1 m/s, dp/dx = 4 Pa/m at y = h/4: 0.025 − 0.0375 = −0.0125 m/s.
- **What it means.** Because the equation is linear, causes simply add; the book's Fig. 8.4 is four cases of one formula.
- **Traps.** The book's A, B have the opposite sign to c₁, c₂; dp/dx < 0 (favourable) gives a *forward* parabola because of
  the minus sign; check both walls at the end.

### D05 · Limits, stress, and the backflow threshold $\frac{dp}{dx}>\frac{2\mu U}{h^2}$ — ★, 7 steps, in C02 (notebook · `couette_poiseuille_backflow`)
- **Goal.** Read the two classic special cases and the shear stress off (8.5), and find exactly when the fluid near the
  fixed floor flows backwards (the book only draws it).
- **Start.** $u(y)=\frac Uhy-\frac1{2\mu}\frac{dp}{dx}y(h-y)$ (8.5).
- **Plan.** (1) Switch off one cause at a time. (2) Differentiate for the stress. (3) Ask when the floor slope turns negative.
- **Tools.** derivative as a slope (P19) · first-order Taylor (P26) · inequalities (P48).
- **Assumptions.** U > 0 in step 7 (for U < 0 the inequality flips).
- **Steps.**
  1. *did:* Switch off the pressure gradient · *tex:* $\frac{dp}{dx}=0\ \Rightarrow\ u=\frac{Uy}{h}$ · *why:* The parabola term is
     proportional to dp/dx; what is left is plane Couette flow (R07). · *plain:* pure wall drag gives a straight line.
  2. *did:* Switch off the wall · *tex:* $U=0\ \Rightarrow\ u=-\frac1{2\mu}\frac{dp}{dx}y(h-y),\ u_{max}=-\frac{h^2}{8\mu}\frac{dp}{dx}$ ·
     *why:* Plane Poiseuille flow (R08); the parabola peaks at mid-gap y = h/2 where y(h − y) = h²/4. · *plain:* pure
     pressure gives a symmetric parabola.
  3. *did:* Differentiate for the shear stress · *tex:* $\tau=\mu\frac{du}{dy}=\frac{\mu U}{h}-\Big(\frac h2-y\Big)\frac{dp}{dx}$ ·
     *why:* τ = μ du/dy (Newton's law); d/dy[y(h − y)] = h − 2y, times −(1/2μ)dp/dx·μ. · *plain:* the stress is linear in y.
  4. *did:* Evaluate at the walls (U = 0) · *tex:* $\tau(0)=-\frac h2\frac{dp}{dx},\ \ \tau(h)=+\frac h2\frac{dp}{dx}$ · *why:* Put
     y = 0 and y = h in step 3; for Poiseuille flow the two wall stresses are equal in size (h/2)∣dp/dx∣. · *plain:* both
     walls feel the same drag.
  5. *did:* State the backflow criterion · *tex:* $u(y)\approx\frac{du}{dy}\Big\rvert_0\,y\ \ (y\to0)$ · *why:* First-order Taylor
     with u(0) = 0: near the floor the sign of u is the sign of the floor slope. We want when it turns negative. · *plain:*
     backflow starts when the floor slope reverses.
  6. *did:* Evaluate the floor slope · *tex:* $\frac{du}{dy}\Big\rvert_0=\frac Uh-\frac{h}{2\mu}\frac{dp}{dx}$ · *why:* Step 3 at
     y = 0, divided by μ. · *plain:* wall drag pushes forward, an adverse gradient pushes back.
  7. *did:* Solve slope < 0 for dp/dx · *tex:* $\frac Uh-\frac{h}{2\mu}\frac{dp}{dx}<0\iff\frac{dp}{dx}>\frac{2\mu U}{h^2}$ · *why:*
     Add (h/2μ)dp/dx and multiply by 2μ/h > 0, which keeps the inequality's direction. · *plain:* an adverse gradient
     beyond 2μU/h² reverses the floor layer.
- **Result.** Couette u = Uy/h, Poiseuille u_max = −(h²/8μ)dp/dx, $\tau=\frac{\mu U}{h}-(\frac h2-y)\frac{dp}{dx}$, backflow iff
  $\frac{dp}{dx}>\frac{2\mu U}{h^2}$ — *in words:* the floor reverses when the pressure's pull there beats the wall's.
- **Check.** Units 2μU/h²: (Pa s)(m/s)/m² = Pa/m ✓. Numbers: h = 1 cm, U = 0.1 m/s → 2 Pa/m; at 4 Pa/m the reversed layer
  reaches y = h − 2μU/(h dp/dx) = 5 mm. Poiseuille u_max = 1.5V (V = −(h²/12μ)dp/dx).
- **What it means.** This is the channel version of boundary-layer separation: an adverse pressure gradient wins near a wall
  where the fluid is slow (Ch. 9).
- **Traps.** Backflow starts when u′(0) < 0, not when the velocity maximum moves; τ is signed (the book draws ∣τ∣).

### D06 · Poiseuille pipe profile $u_z(R)=\frac{R^2-a^2}{4\mu}\frac{dp}{dz}$ (8.6) — ★★, 10 steps, in C03 (notebook)
- **Goal.** Find the velocity across a round pipe driven by a constant pressure gradient.
- **Start.** u = (0, 0, u_z(R)) in cylindrical coordinates (R, φ, z); steady, fully developed, axisymmetric; (8.1) in
  cylindrical form (Appendix B).
- **Plan.** (1) Show continuity and advection are trivial. (2) Separate the z-equation. (3) Integrate twice; the axis and
  the wall fix the constants.
- **Tools.** Laplacian in cylindrical coordinates (primer P186) · regularity at the axis (gloss) · integrating twice (gloss).
- **Assumptions.** Steady, fully developed (∂/∂z of u = 0), axisymmetric (∂/∂φ = 0) — step 1; no slip at R = a (step 9).
- **Steps.**
  1. *did:* Check continuity · *tex:* $\nabla\cdot\mathbf u=\frac{\partial u_z}{\partial z}=0$ · *why:* With u_R = u_φ = 0 only the
     z-term of the cylindrical divergence is left, and u_z does not depend on z (fully developed). · *plain:* mass is
     conserved automatically.
  2. *did:* Evaluate the advective term · *tex:* $(\mathbf u\cdot\nabla)u_z=u_z\frac{\partial u_z}{\partial z}=0$ · *why:* Only
     the z-velocity is non-zero and nothing changes along z. · *plain:* no acceleration: the equation is linear.
  3. *did:* Read the R- and φ-equations · *tex:* $0=\frac{\partial p}{\partial R},\ \ 0=\frac{\partial p}{\partial\varphi}\ \Rightarrow\
     p=p(z)$ · *why:* All velocity terms in those components vanish, leaving only pressure. · *plain:* pressure is uniform
     over each cross-section.
  4. *did:* Write the z-equation · *tex:* $0=-\frac{dp}{dz}+\frac\mu R\frac{d}{dR}\Big(R\frac{du_z}{dR}\Big)$ · *why:* The viscous
     term is μ times the axisymmetric Laplacian (1/R)d/dR(R du/dR) (P186). · *plain:* pressure push balances friction
     across the pipe.
  5. *did:* Separate the variables · *tex:* $\frac{dp}{dz}=\frac\mu R\frac{d}{dR}\Big(R\frac{du_z}{dR}\Big)=\text{const}$ · *why:*
     The left side depends on z only, the right on R only — the D03 argument makes both one constant. · *plain:* the pressure
     falls linearly along the pipe.
  6. *did:* Multiply by R/μ and integrate · *tex:* $R\frac{du_z}{dR}=\frac{R^2}{2\mu}\frac{dp}{dz}+A$ · *why:* d/dR(R u′) =
     (R/μ)dp/dz; its antiderivative is (R²/2μ)dp/dz plus a constant A. Keep the R in front of u′. · *plain:* first
     integral.
  7. *did:* Divide by R and integrate again · *tex:* $u_z(R)=\frac{R^2}{4\mu}\frac{dp}{dz}+A\ln R+B$ · *why:* u′ =
     (R/2μ)dp/dz + A/R; the antiderivative of A/R is A ln R. · *plain:* the general solution has a parabola, a logarithm and
     a constant.
  8. *did:* Keep the velocity finite on the axis · *tex:* $A=0$ · *why:* ln R → −∞ as R → 0, but the velocity at the centre
     is finite; equivalently from step 6, R u′ → A must be 0 by symmetry. · *plain:* the logarithm is removed.
  9. *did:* Apply no slip at the wall · *tex:* $0=\frac{a^2}{4\mu}\frac{dp}{dz}+B\ \Rightarrow\ B=-\frac{a^2}{4\mu}\frac{dp}{dz}$ ·
     *why:* (8.3): the fluid at R = a is at rest. · *plain:* one condition fixes B.
  10. *did:* Combine · *tex:* $u_z(R)=\frac{R^2-a^2}{4\mu}\frac{dp}{dz}$ (8.6) · *why:* Insert A = 0 and B into step 7 and factor
      1/(4μ). · *plain:* (8.6), a paraboloid; positive when dp/dz < 0.
- **Result.** (8.6) — *in words:* the axial velocity is a paraboloid, fastest on the axis, zero at the wall.
- **Check.** u_z(a) = 0 ✓; units (m²)/(Pa s)(Pa/m) = m/s ✓; a = 1 mm, dp/dz = −1000 Pa/m → u(0) = 0.25 m/s; parity with
  Ch. 3's `pipe_profile` far downstream and Ch. 4's `exact_solution("pipe_poiseuille")`.
- **What it means.** The pipe profile is the channel's parabola rotated about the axis, but the factor is 1/(4μ), not
  1/(2μ), because a disc's circumference grows with R.
- **Traps.** Forgetting the R in R du/dR before the second integration; keeping A ln R; writing 1/(2μ).

### D07 · Stress (8.7), wall stress (8.8) as a force balance, Hagen–Poiseuille and f = 64/Re — ★★, 10 steps, in C03 (notebook)
- **Goal.** Get the shear stress, the wall stress (also from a force balance that holds for turbulent flow), the flow rate,
  the mean and maximum speeds and the friction factor.
- **Start.** $u_z(R)=\frac{R^2-a^2}{4\mu}\frac{dp}{dz}$ (8.6).
- **Plan.** (1) Differentiate for τ. (2) Confirm τ₀ by a force balance on a slug. (3) Integrate over the cross-section. (4)
  Write the friction factor.
- **Tools.** cylindrical shear stress (R09) · control-volume force balance (Ch. 4 §4.4) · area element 2πR dR (gloss).
- **Assumptions.** Fully developed (momentum flux in = out, step 4); Darcy definition of f (step 9).
- **Steps.**
  1. *did:* Write the shear stress · *tex:* $\tau=\tau_{zR}=\mu\frac{du_z}{dR}$ · *why:* $\tau_{zR}=\mu(\frac{\partial u_R}{\partial z}
     +\frac{\partial u_z}{\partial R})$ (R09) with u_R = 0. · *plain:* only the radial slope of u_z matters.
  2. *did:* Differentiate (8.6) · *tex:* $\tau=\frac{R}{2}\frac{dp}{dz}$ (8.7) · *why:* d/dR(R² − a²) = 2R; 2R/(4μ) × μ = R/2. ·
     *plain:* (8.7): stress grows linearly from zero on the axis.
  3. *did:* Evaluate at the wall · *tex:* $\tau_0=\frac a2\frac{dp}{dz}$ (8.8) · *why:* Put R = a; it is negative for forward flow
     (dp/dz < 0): the wall holds the fluid back. · *plain:* (8.8), the largest stress in magnitude.
  4. *did:* Balance forces on a fluid slug · *tex:* $\pi a^2[p(0)-p(L)]+2\pi aL\,\tau_0=0$ · *why:* Steady and fully developed, so
     the momentum flux in equals out (Ch. 4 CV momentum): pressure force on the ends plus wall shear on the side sum to zero.
     · *plain:* the pressure drop is spent on wall friction.
  5. *did:* Solve for the wall stress · *tex:* $\tau_0=\frac{a}{2}\frac{p(L)-p(0)}{L}=\frac a2\frac{dp}{dz}$ · *why:* Divide by
     2πaL; the profile was never used, so (8.8) holds for averaged turbulent pipe flow too. · *plain:* same result, more
     general.
  6. *did:* Set up the flow rate · *tex:* $Q=\int_0^au_z\,2\pi R\,dR$ · *why:* A thin ring of radius R and width dR has area
     2πR dR; summing velocity × area gives the volume per second. · *plain:* add up rings.
  7. *did:* Evaluate the integral · *tex:* $\int_0^a(R^2-a^2)R\,dR=\frac{a^4}4-\frac{a^4}2=-\frac{a^4}{4}\ \Rightarrow\ Q=-\frac{\pi
     a^4}{8\mu}\frac{dp}{dz}$ · *why:* Polynomial antiderivatives; then multiply by 2π/(4μ) and dp/dz. · *plain:* Q ∝ a⁴
     (Hagen–Poiseuille).
  8. *did:* Mean and maximum speeds · *tex:* $V=\frac{Q}{\pi a^2}=-\frac{a^2}{8\mu}\frac{dp}{dz},\ \ u_{max}=u_z(0)=-\frac{a^2}{4\mu}
     \frac{dp}{dz}=2V$ · *why:* Divide by the area; evaluate (8.6) on the axis. · *plain:* the centre moves twice the average.
  9. *did:* Express the wall stress through V · *tex:* $\lvert\tau_0\rvert=-\frac a2\frac{dp}{dz}=\frac{4\mu V}{a}$ · *why:* From step
     8, dp/dz = −8μV/a²; insert in (8.8). We want the Darcy factor f ≡ 8∣τ₀∣/(ρV²). · *plain:* wall friction ∝ V.
  10. *did:* Form the friction factor · *tex:* $f=\frac{8\lvert\tau_0\rvert}{\rho V^2}=\frac{32\mu}{\rho Va}=\frac{64}{\mathrm{Re}},\ \
      \mathrm{Re}=\frac{V\,2a}{\nu}$ · *why:* Substitute step 9 and write 2a = d, the diameter the pipe Re uses. · *plain:*
      the laminar line of the Moody chart.
- **Result.** (8.7), (8.8), $Q=-\frac{\pi a^4}{8\mu}\frac{dp}{dz}$, V, u_max = 2V, f = 64/Re — *in words:* the pipe's resistance
  grows as 1/a⁴.
- **Check.** Units: Q in m⁴/(Pa s) × Pa/m = m³/s ✓. Numbers a = 1 mm, dp/dz = −1000 Pa/m: V = 0.125 m/s, Q = 3.93 × 10⁻⁷ m³/s,
  τ₀ = −0.5 Pa, Re = 250, f = 0.256.
- **What it means.** Halving a pipe's radius at fixed pressure drop cuts the flow 16-fold — why narrowed arteries matter.
  f = 64/Re is where Ch. 12's turbulent friction curves start.
- **Traps.** τ₀ is negative for forward flow; ∫(R² − a²)R dR = −a⁴/4 (sign); Re uses the diameter and the mean velocity.

### D08 · Circular Couette flow $u_\varphi=AR+B/R$ (8.9) and (8.10) — ★★, 10 steps, in C04 (notebook)
- **Goal.** Find the swirl between an inner cylinder (R₁, Ω₁) and an outer cylinder (R₂, Ω₂).
- **Start.** u = (0, u_φ(R), 0) in (R, φ, z); steady; infinitely long cylinders.
- **Plan.** (1) Reduce the momentum equations. (2) Integrate the φ-equation twice. (3) Fix A and B at the walls.
- **Tools.** φ-component of the cylindrical vector Laplacian (primer P186) · Euler–Cauchy ODE (primer P187) · 2 × 2 system
  (P57).
- **Assumptions.** Steady, axisymmetric (∂/∂φ = 0), no axial flow, infinitely long (steps 1–4); no slip (step 8).
- **Steps.**
  1. *did:* Check continuity · *tex:* $\nabla\cdot\mathbf u=\frac1R\frac{\partial u_\varphi}{\partial\varphi}=0$ · *why:* Only the
     φ-term of the divergence could survive and nothing depends on φ. · *plain:* mass is conserved automatically.
  2. *did:* Evaluate the advective acceleration · *tex:* $(\mathbf u\cdot\nabla)\mathbf u=-\frac{u_\varphi^2}{R}\mathbf e_R$ · *why:*
     For a purely azimuthal flow the moving basis leaves only the centripetal part (P105, Ch. 4 §4.7). · *plain:* the fluid
     accelerates toward the axis.
  3. *did:* Write the R-momentum · *tex:* $-\frac{u_\varphi^2}{R}=-\frac1\rho\frac{dp}{dR}$ · *why:* The R-component of the viscous
     term, −(2/R²)∂u_φ/∂φ, vanishes; pressure supplies the centripetal force. · *plain:* pressure rises outward.
  4. *did:* Write the φ-momentum · *tex:* $0=\mu\frac{d}{dR}\Big[\frac1R\frac{d}{dR}(Ru_\varphi)\Big]$ · *why:* ∂p/∂φ = 0 by symmetry
     and the advective term has no φ-part; the viscous term is the φ-component of the vector Laplacian (P186). · *plain:* a
     pure friction balance.
  5. *did:* Integrate once · *tex:* $\frac1R\frac{d}{dR}(Ru_\varphi)=2A$ · *why:* Fundamental theorem of calculus; naming the
     constant 2A keeps the final form clean. · *plain:* first integral.
  6. *did:* Multiply by R and integrate · *tex:* $Ru_\varphi=AR^2+B$ · *why:* d(Ru_φ)/dR = 2AR has antiderivative AR² plus a new
     constant B. · *plain:* the angular momentum per mass is a quadratic in R.
  7. *did:* Divide by R · *tex:* $u_\varphi(R)=AR+\frac BR$ (8.9) · *why:* Solves the Euler–Cauchy ODE R²u″ + Ru′ − u = 0 (trial R^λ
     gives λ = ±1, P187). · *plain:* (8.9): solid-body rotation plus a free vortex.
  8. *did:* Apply no slip on both cylinders · *tex:* $\Omega_1R_1^2=AR_1^2+B,\ \ \Omega_2R_2^2=AR_2^2+B$ · *why:* The wall speed is
     ΩR (not Ω); we multiplied each condition by its radius so B stands alone. · *plain:* two equations for A, B.
  9. *did:* Subtract to find A · *tex:* $A=\frac{\Omega_2R_2^2-\Omega_1R_1^2}{R_2^2-R_1^2}$ · *why:* Subtracting eliminates B; divide by
     R₂² − R₁². · *plain:* A is set by the difference of the walls' angular momenta.
  10. *did:* Back-substitute for B and insert · *tex:* $B=-\frac{(\Omega_2-\Omega_1)R_1^2R_2^2}{R_2^2-R_1^2},\ \ u_\varphi=\frac{[\Omega_2
      R_2^2-\Omega_1R_1^2]R-[\Omega_2-\Omega_1]R_1^2R_2^2/R}{R_2^2-R_1^2}$ (8.10) · *why:* B = Ω₁R₁² − AR₁²; over a common
      denominator the Ω₂R₂²R₁² terms combine; then put A, B into (8.9). · *plain:* (8.10).
- **Result.** (8.10) with A, B above — *in words:* the walls choose how much rigid rotation and how much free vortex.
- **Check.** u_φ(R₁) = Ω₁R₁ and u_φ(R₂) = Ω₂R₂ ✓; Ω₁ = Ω₂ = Ω gives A = Ω, B = 0 (rigid rotation) ✓; numbers R₁ = 1 cm,
  R₂ = 2 cm, Ω₁ = 1 rad/s, Ω₂ = 0: A = −1/3 s⁻¹, B = 1.33 × 10⁻⁴ m²/s, u_φ(1.5 cm) = 3.9 mm/s.
- **What it means.** Vorticity is 2A everywhere (the B/R part is irrotational); Ch. 11's Rayleigh criterion for Taylor
  vortices is written in A and B.
- **Traps.** The φ-equation contains no pressure and no advection; the first integral is (1/R)(Ru_φ)′ = 2A; the wall speed
  is ΩR.

### D09 · The limits $u_\varphi=\Omega_1R_1^2/R$ (8.11) and $u_\varphi=\Omega_2R$ (8.12) — ★, 6 steps, in C04 (notebook)
- **Goal.** Show that circular Couette flow contains the ideal vortex and solid-body rotation as limits.
- **Start.** A and B of D08 in (8.9) $u_\varphi=AR+B/R$.
- **Plan.** (1) Send the outer wall to infinity at rest. (2) Remove the inner cylinder.
- **Tools.** limits (P68) · circulation (Ch. 5 §5.1).
- **Assumptions.** Ω₂ = 0 in the first limit (step 2); Ω₁ = 0 in the second (step 5).
- **Steps.**
  1. *did:* Divide through by R₂² · *tex:* $A=\frac{\Omega_2-\Omega_1R_1^2/R_2^2}{1-R_1^2/R_2^2},\ \ B=-\frac{(\Omega_2-\Omega_1)R_1^2}
     {1-R_1^2/R_2^2}$ · *why:* Dividing numerator and denominator by R₂² makes the limit R₂ → ∞ visible: every R₁²/R₂² → 0. ·
     *plain:* the same constants, ready for the limit.
  2. *did:* Let R₂ → ∞ with Ω₂ = 0 · *tex:* $A\to0,\ \ B\to\Omega_1R_1^2$ · *why:* R₁²/R₂² → 0; with Ω₂ = 0 the numerator of A
     vanishes (with Ω₂ ≠ 0, A → Ω₂ instead). · *plain:* only the vortex part survives.
  3. *did:* Insert into (8.9) · *tex:* $u_\varphi(R)=\frac{\Omega_1R_1^2}{R}$ (8.11) · *why:* u_φ = AR + B/R with A = 0. · *plain:*
     (8.11): speed falls as 1/R.
  4. *did:* Compute the circulation · *tex:* $\Gamma=\oint\mathbf u\cdot d\mathbf l=2\pi Ru_\varphi=2\pi\Omega_1R_1^2$ · *why:* Around
     a circle u is tangent and constant in size; the result does not depend on R — the signature of the ideal vortex (5.2)
     u = Γ/(2πR). · *plain:* a viscous flow identical to the ideal vortex.
  5. *did:* Let R₁ → 0 with Ω₁ = 0 · *tex:* $A=\frac{\Omega_2R_2^2}{R_2^2}=\Omega_2,\ \ B=-\frac{\Omega_2R_1^2R_2^2}{R_2^2-R_1^2}\to0$ ·
     *why:* Put Ω₁ = 0 and R₁ = 0 in D08's A and B. · *plain:* only rigid rotation survives.
  6. *did:* Insert into (8.9) · *tex:* $u_\varphi(R)=\Omega_2R$ (8.12) · *why:* u_φ = AR with B = 0. · *plain:* (8.12): the tank's
     fluid turns like a solid (5.1).
- **Result.** (8.11) and (8.12) — *in words:* remove one wall and you get one of the two basic vortices of Ch. 5.
- **Check.** (8.11) at R = R₁ gives Ω₁R₁ ✓; the R₂ = 10⁶R₁ value of (8.10) agrees with (8.11) to 1e-10; Ch. 5's
  `rotating_cylinder_flow(r, a, omega=2Ω₁)` gives the same.
- **What it means.** The flow outside a spinning rod is viscous but irrotational: shear exists (R10), no net viscous force
  acts, and all the power put in is dissipated (R11).
- **Traps.** Divide by R₂² before the limit; Ω₂ = 0 is needed or A → Ω₂; Ch. 5's function takes ω = 2Ω₁.

### D10 · The scaled gap equations (8.15), (8.16a), (8.16b) — ★★★, 14 steps, in C05 (notebook · `lubrication_scaling`)
- **Goal.** Rewrite the gap equations in variables of order one, so that the size of every term sits in its coefficient.
- **Start.** $\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0$ (6.2), $\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x}+
  v\frac{\partial u}{\partial y}=-\frac1\rho\frac{\partial p}{\partial x}+\frac\mu\rho\big(\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}
  {\partial y^2}\big)$ (8.13a), its v-twin (8.13b) with ∂p/∂y, and the scalings $x^*=\frac xL,\ y^*=\frac y{\varepsilon L},\ t^*=\frac
  {Ut}L,\ u^*=\frac uU,\ v^*=\frac v{\varepsilon U},\ p^*=\frac p{P_a}$ (8.14).
- **Plan.** (1) Turn the scalings into substitutions and derivative rules. (2) Scale continuity. (3) Scale x-momentum and
  divide by the cross-gap friction's coefficient. (4) Scale y-momentum and choose the factor that keeps the same 1/Λ.
- **Tools.** anisotropic scaling (primer P188) · scaled variables and the chain rule (P133) · bookkeeping of ε (gloss).
- **Assumptions.** One length along (L), one across (h = εL); pressure scaled by the atmospheric P_a (as the book does).
- **Steps.**
  1. *did:* Invert the scalings · *tex:* $x=Lx^*,\ y=\varepsilon Ly^*,\ t=\tfrac LU t^*,\ u=Uu^*,\ v=\varepsilon Uv^*,\ p=P_ap^*$ ·
     *why:* We will substitute these into the equations; each starred variable is of order one by design. · *plain:*
     dimensional quantities = scale × order-one variable.
  2. *did:* Write the derivative rules · *tex:* $\frac{\partial}{\partial x}=\frac1L\frac{\partial}{\partial x^*},\ \frac{\partial}
     {\partial y}=\frac{1}{\varepsilon L}\frac{\partial}{\partial y^*},\ \frac{\partial}{\partial t}=\frac UL\frac{\partial}{\partial t^*}$ ·
     *why:* Chain rule for scaled variables (P133): ∂x*/∂x = 1/L, and so on. · *plain:* derivatives across the gap are 1/ε
     times larger.
  3. *did:* Scale continuity · *tex:* $\frac UL\frac{\partial u^*}{\partial x^*}+\frac{\varepsilon U}{\varepsilon L}\frac{\partial v^*}
     {\partial y^*}=0$ · *why:* Substitute u = Uu*, v = εUv* and the rules of step 2 into (6.2). · *plain:* both terms carry
     U/L.
  4. *did:* Divide by U/L · *tex:* $\frac{\partial u^*}{\partial x^*}+\frac{\partial v^*}{\partial y^*}=0$ (8.15) · *why:* The ε
     cancels because v was scaled with εU — chosen exactly so that continuity balances (if v ~ U the u-term would drop). ·
     *plain:* (8.15): mass is conserved with no approximation.
  5. *did:* Scale the x-inertia terms · *tex:* $\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=
     \frac{U^2}{L}\Big(\frac{\partial u^*}{\partial t^*}+u^*\frac{\partial u^*}{\partial x^*}+v^*\frac{\partial u^*}{\partial y^*}\Big)$ ·
     *why:* Each term: (U/L)U, U(U/L), (εU)(U/εL) — all U²/L. · *plain:* inertia is of size U²/L.
  6. *did:* Scale the x-pressure term · *tex:* $-\frac1\rho\frac{\partial p}{\partial x}=-\frac{P_a}{\rho L}\frac{\partial p^*}{\partial
     x^*}$ · *why:* p = P_a p* and ∂/∂x = (1/L)∂/∂x*. · *plain:* pressure term of size P_a/(ρL).
  7. *did:* Scale the two x-viscous terms · *tex:* $\frac\mu\rho\frac{\partial^2u}{\partial x^2}=\frac{\mu U}{\rho L^2}\frac{\partial^2u^*}
     {\partial x^{*2}},\ \ \frac\mu\rho\frac{\partial^2u}{\partial y^2}=\frac{\mu U}{\rho\varepsilon^2L^2}\frac{\partial^2u^*}{\partial y^{*2}}$ ·
     *why:* Two derivatives along the gap bring 1/L², two across bring 1/(εL)². · *plain:* friction across the gap is 1/ε²
     larger.
  8. *did:* Multiply by ρε²L²/(μU) · *tex:* $\times\frac{\rho\varepsilon^2L^2}{\mu U}$ · *why:* This makes the cross-gap friction —
     the term we expect to dominate — have coefficient 1, so every other coefficient reads directly as a relative size. ·
     *plain:* measure everything against cross-gap friction.
  9. *did:* Read the inertia coefficient · *tex:* $\frac{U^2}{L}\cdot\frac{\rho\varepsilon^2L^2}{\mu U}=\varepsilon^2\frac{\rho UL}{\mu}=
     \varepsilon^2\mathrm{Re}_L$ · *why:* Simplify the product; Re_L = ρUL/μ uses the passage length. · *plain:* inertia is
     weighed by ε²Re_L.
  10. *did:* Read the other x-coefficients · *tex:* $\varepsilon^2\mathrm{Re}_L(\ldots)=-\frac1\Lambda\frac{\partial p^*}{\partial x^*}+
      \varepsilon^2\frac{\partial^2u^*}{\partial x^{*2}}+\frac{\partial^2u^*}{\partial y^{*2}}$ (8.16a) · *why:* Pressure:
      (P_a/ρL)(ρε²L²/μU) = P_a h²/(μUL) = 1/Λ; along-gap friction: ε²; across: 1. · *plain:* (8.16a).
  11. *did:* Scale the y-inertia terms · *tex:* $\frac{\partial v}{\partial t}+u\frac{\partial v}{\partial x}+v\frac{\partial v}{\partial y}=
      \frac{\varepsilon U^2}{L}(\ldots)^*$ · *why:* Each term now carries one factor εU from v: (U/L)εU, U(εU/L),
      (εU)(εU/εL). · *plain:* vertical inertia is ε times smaller.
  12. *did:* Scale the y-pressure and friction · *tex:* $-\frac{P_a}{\rho\varepsilon L}\frac{\partial p^*}{\partial y^*},\ \ \frac{\mu\varepsilon
      U}{\rho L^2}\frac{\partial^2v^*}{\partial x^{*2}},\ \ \frac{\mu U}{\rho\varepsilon L^2}\frac{\partial^2v^*}{\partial y^{*2}}$ · *why:*
      Use ∂/∂y = (1/εL)∂/∂y* on the pressure — the corrected (8.13b) has ∂p/∂y; the book prints ∂p/∂x, which would give a
      different set. · *plain:* pressure across the gap is large.
  13. *did:* Multiply by ρε³L²/(μU) · *tex:* $\times\frac{\rho\varepsilon^3L^2}{\mu U}$ · *why:* We choose the factor that gives the
      pressure term the same coefficient 1/Λ as in (8.16a): (P_a/ρεL)(ρε³L²/μU) = ε²LP_a/(μU) = 1/Λ. · *plain:* compare
      against the same pressure yardstick.
  14. *did:* Read the y-coefficients · *tex:* $\varepsilon^4\mathrm{Re}_L(\ldots)=-\frac1\Lambda\frac{\partial p^*}{\partial y^*}+
      \varepsilon^4\frac{\partial^2v^*}{\partial x^{*2}}+\varepsilon^2\frac{\partial^2v^*}{\partial y^{*2}}$ (8.16b) · *why:* Inertia
      (εU²/L)(ρε³L²/μU) = ε⁴Re_L; friction ε⁴ and ε² likewise. · *plain:* (8.16b): only pressure lacks a power of ε.
- **Result.** (8.15), (8.16a), (8.16b) — *in words:* in a thin gap inertia carries ε²Re_L and all y-terms except pressure
  carry ε² or less.
- **Check.** Every coefficient is dimensionless ✓. sympy (`lubrication_nondim_sympy`, check_src in C05) returns the sets
  {ε²Re_L, 1/Λ, ε², 1} and {ε⁴Re_L, 1/Λ, ε⁴, ε²}; the printed ∂p/∂x in (8.13b) gives a set that does not match (8.16b).
  Numbers (engine film): ε = 10⁻³, Re_L = 4350 → ε²Re_L = 4.35 × 10⁻³.
- **What it means.** "Large Re" is not the right test in a thin gap: inertia competes with ε²Re_L. The same bookkeeping
  becomes the boundary-layer approximation in Ch. 9.
- **Traps.** v scales with εU (from continuity), not U; the two momentum equations are multiplied by different factors
  (ε²L²ρ/μU and ε³L²ρ/μU) so that the pressure keeps 1/Λ; the printed ∂p/∂x in (8.13b) must be ∂p/∂y.
- **sympy check intent (★★★).** Build the scaled derivatives as in steps 1–2, substitute into (6.2), (8.13a), (8.13b),
  multiply by the two factors and extract each term's coefficient; assert the two sets; repeat with the printed (8.13b) and
  show the mismatch (check_src written in A.3 item 11).

### D11 · The lubrication balance $0\cong-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{\partial^2u}{\partial y^2}$ (8.17a), $0\cong-\frac1\rho\frac{\partial p}{\partial y}$ (8.17b) — ★★, 7 steps, in C05 (notebook · `lubrication_scaling`)
- **Goal.** Drop the small terms of (8.16) and return to dimensional form — restoring the ν the printed (8.17a) omits.
- **Start.** (8.16a) and (8.16b) from D10.
- **Plan.** (1) State the orders. (2) Drop the small terms in each equation. (3) Undo the scaling.
- **Tools.** orders of smallness (P68).
- **Assumptions.** ε ≪ 1, ε²Re_L ≪ 1, Λ = O(1) with a suitable pressure scale (step 1).
- **Steps.**
  1. *did:* State the orders · *tex:* $\varepsilon\ll1,\ \ \varepsilon^2\mathrm{Re}_L\ll1,\ \ \Lambda=O(1)$ · *why:* All starred
     derivatives are O(1) by construction; only the coefficients decide sizes. · *plain:* thin gap, small reduced Reynolds
     number.
  2. *did:* Drop the x-inertia · *tex:* $0=-\frac1\Lambda\frac{\partial p^*}{\partial x^*}+\varepsilon^2\frac{\partial^2u^*}{\partial x^{*2}}+
     \frac{\partial^2u^*}{\partial y^{*2}}+O(\varepsilon^2\mathrm{Re}_L)$ · *why:* The left side of (8.16a) is ε²Re_L × O(1). ·
     *plain:* acceleration is negligible — even the unsteady part.
  3. *did:* Drop the along-gap friction · *tex:* $0\cong-\frac1\Lambda\frac{\partial p^*}{\partial x^*}+\frac{\partial^2u^*}{\partial y^{*2}}$ ·
     *why:* ε²∂²u*/∂x*² is ε² times an order-one derivative, so it is as small as the inertia we just dropped. · *plain:* friction acts across the gap only.
  4. *did:* Reduce the y-equation · *tex:* $\frac1\Lambda\frac{\partial p^*}{\partial y^*}=O(\varepsilon^2)\ \Rightarrow\ 0\cong-\frac{\partial
     p^*}{\partial y^*}$ · *why:* In (8.16b) every term except the pressure carries ε² or smaller. · *plain:* pressure is the
     same across the gap to relative order ε².
  5. *did:* Undo the x-scaling · *tex:* $0\cong-\frac1\rho\frac{\partial p}{\partial x}+\frac\mu\rho\frac{\partial^2u}{\partial y^2}$ ·
     *why:* Multiply by μU/(ρε²L²) (the inverse of D10 step 8) and use p = P_a p*, ∂/∂y = (1/εL)∂/∂y*. · *plain:* back in
     dimensions.
  6. *did:* Write μ/ρ as ν · *tex:* $0\cong-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{\partial^2u}{\partial y^2}$ (8.17a) · *why:*
     ν = μ/ρ; the book prints this without the ν, which fails the units (m/s² vs 1/(m s)). · *plain:* (8.17a): pressure
     gradient = cross-gap friction.
  7. *did:* Undo the y-scaling · *tex:* $0\cong-\frac1\rho\frac{\partial p}{\partial y}$ (8.17b) · *why:* Multiply the reduced
     y-equation by P_a/(ρεL). · *plain:* (8.17b).
- **Result.** (8.17a), (8.17b) — *in words:* the lubrication approximation: a local Poiseuille balance with pressure
  uniform across the gap.
- **Check.** Units of (8.17a): both terms m/s² ✓. Time still appears as a parameter (no ∂/∂t term, yet h, p may vary
  slowly). Numbers: engine film ε²Re_L = 4.35 × 10⁻³.
- **What it means.** With P_a the bearing number is Λ = 49 for our film, not O(1); rescaling p by μUL/h² sets Λ = 1 and
  changes nothing in the ordering — the physical pressures are tens to hundreds of atmospheres.
- **Traps.** Dropping inertia needs ε²Re_L ≪ 1, not Re_L ≪ 1; ∂p*/∂y* is O(ε²), not exactly 0; re-dimensionalising restores ν.

### D12 · The lubrication profile: (8.18) → (8.19) — ★★, 8 steps, in C06 (notebook · `slider_bearing`)
- **Goal.** Find the velocity at every point of a slowly varying gap with moving walls.
- **Start.** $0\cong-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{\partial^2u}{\partial y^2}$ (8.17a), $0\cong-\frac1\rho\frac{\partial p}
  {\partial y}$ (8.17b), with u = U₀(t) on y = 0 and u = U_h(t) on y = h(x, t).
- **Plan.** (1) Use (8.17b) to treat ∂p/∂x as a constant across the gap. (2) Integrate twice. (3) Fix the "constants" with
  the two walls and regroup.
- **Tools.** integrating twice (gloss in C02) · 2 × 2 system (P57).
- **Assumptions.** Lubrication balance (D11); walls at y = 0 and y = h(x, t).
- **Steps.**
  1. *did:* Use (8.17b) · *tex:* $p=p(x,t)\ \Rightarrow\ \frac{\partial p}{\partial x}\ \text{is independent of }y$ · *why:* Pressure
     does not change across the gap, so at fixed x and t its x-slope is one number. · *plain:* each station has one pressure
     gradient.
  2. *did:* Rewrite (8.17a) · *tex:* $\frac{\partial^2u}{\partial y^2}=\frac1\mu\frac{\partial p}{\partial x}$ · *why:* Multiply
     (8.17a) by 1/ν and use ρν = μ. · *plain:* constant curvature across the gap.
  3. *did:* Integrate twice in y · *tex:* $u(x,y,t)\cong\frac1\mu\frac{\partial p(x,t)}{\partial x}\frac{y^2}{2}+Ay+B$ (8.18) · *why:*
     The integration is at fixed x and t, so the "constants" A and B may depend on x and t. · *plain:* (8.18).
  4. *did:* Apply the lower wall · *tex:* $u(x,0,t)=U_0\ \Rightarrow\ B=U_0$ · *why:* No slip at y = 0. · *plain:* the floor's speed.
  5. *did:* Apply the upper wall · *tex:* $A=\frac{U_h-U_0}{h}-\frac{h}{2\mu}\frac{\partial p}{\partial x}$ · *why:* No slip at y = h:
     U_h = (1/μ)p_x h²/2 + Ah + U₀, solved for A. · *plain:* the slope mixes the wall-speed difference and pressure.
  6. *did:* Substitute A and B · *tex:* $u=\frac{1}{2\mu}\frac{\partial p}{\partial x}(y^2-hy)+(U_h-U_0)\frac yh+U_0$ · *why:* Insert
     into (8.18) and collect the p_x terms. · *plain:* parabola plus a line.
  7. *did:* Regroup into gap fractions · *tex:* $u\cong-\frac{h^2}{2\mu}\frac{\partial p}{\partial x}\frac yh\Big(1-\frac yh\Big)+U_h\frac yh+
     U_0\Big(1-\frac yh\Big)$ · *why:* y² − hy = −h²(y/h)(1 − y/h) and (U_h − U₀)y/h + U₀ = U_h y/h + U₀(1 − y/h). · *plain:*
     local Poiseuille + local Couette.
  8. *did:* Compare with the printed (8.19) · *tex:* $U_0=0:\ \ u\cong-\frac{h^2}{2\mu}\frac{\partial p}{\partial x}\frac yh\Big(1-\frac yh\Big)+
     U_h\frac yh$ · *why:* The printed (8.19) adds U₀ instead of U₀(1 − y/h), which gives u(h) = U_h + U₀; with U₀ = 0 (every
     example) the two agree. · *plain:* the book's form is right whenever the floor is still.
- **Result.** $u\cong-\frac{h^2}{2\mu}\frac{\partial p}{\partial x}\frac yh(1-\frac yh)+U_h\frac yh+U_0(1-\frac yh)$ (the consistent
  (8.19)) — *in words:* at every station, the channel profile of (8.5) with the local gap.
- **Check.** u(0) = U₀ ✓, u(h) = U_h ✓ (the printed form gives U_h + U₀ ✗ when U₀ ≠ 0); U₀ = 0, h constant reproduces (8.5).
- **What it means.** The velocity field is known once p(x, t) is; D13 finds the equation for p.
- **Traps.** ∂p/∂x is constant only in the y-integration; A and B depend on x and t.

### D13 · The gap flux and the 1-D Reynolds equation $\frac{\partial h}{\partial t}+\frac{\partial q}{\partial x}=0$ — ★★★, 12 steps, in C06 (notebook · `slider_bearing`)
- **Goal.** Turn mass conservation into one equation for the pressure in a gap of any slowly varying shape (the book leaves
  it to Exercises 8.19–8.20).
- **Start.** $\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0$ (6.2) in 0 < y < h(x, t), and the consistent profile of
  D12.
- **Plan.** (1) Integrate continuity across the gap. (2) Move ∂/∂x outside with the Leibniz rule. (3) Use the kinematic
  conditions at the walls. (4) Compute the flux from the profile.
- **Tools.** Leibniz rule with a moving upper limit (primer P189) · kinematic condition at a moving wall (gloss) · polynomial
  integrals.
- **Assumptions.** Flat, impermeable, lower wall (step 4); the upper wall is a material surface y = h(x, t) (step 5); the
  lubrication profile (step 8).
- **Steps.**
  1. *did:* Integrate continuity from 0 to h · *tex:* $\int_0^h\frac{\partial u}{\partial x}dy+v(x,h,t)-v(x,0,t)=0$ · *why:*
     Integrate (6.2) over the gap; the v-term integrates exactly by the fundamental theorem of calculus. · *plain:* the net
     outflow of a thin slice is zero.
  2. *did:* Apply the Leibniz rule · *tex:* $\int_0^h\frac{\partial u}{\partial x}dy=\frac{\partial}{\partial x}\int_0^hu\,dy-u(x,h,t)
     \frac{\partial h}{\partial x}$ · *why:* The upper limit h depends on x, so differentiating the integral adds a boundary
     term (P189); we want the flux itself inside the derivative. · *plain:* the x-change of the flux, corrected for the
     moving limit.
  3. *did:* Name the flux · *tex:* $\frac{\partial q}{\partial x}-u(h)\frac{\partial h}{\partial x}+v(h)-v(0)=0,\ \ q\equiv\int_0^hu\,dy$ ·
     *why:* Insert step 2 into step 1. · *plain:* flux divergence plus wall terms.
  4. *did:* Use the flat, impermeable floor · *tex:* $v(x,0,t)=0$ · *why:* No flow through the lower wall (8.2), which does not
     move vertically. · *plain:* nothing enters through the floor.
  5. *did:* Use the kinematic upper wall · *tex:* $v(x,h,t)=\frac{\partial h}{\partial t}+u(x,h,t)\frac{\partial h}{\partial x}$ ·
     *why:* A particle on y = h(x, t) stays on it: D(y − h)/Dt = 0 there (P132). · *plain:* the fluid at the top moves with
     the wall.
  6. *did:* Substitute both wall conditions · *tex:* $\frac{\partial q}{\partial x}-u(h)\frac{\partial h}{\partial x}+\frac{\partial h}
     {\partial t}+u(h)\frac{\partial h}{\partial x}=0$ · *why:* Steps 4 and 5 into step 3. · *plain:* the boundary terms appear
     twice with opposite signs.
  7. *did:* Cancel the boundary terms · *tex:* $\frac{\partial h}{\partial t}+\frac{\partial q}{\partial x}=0$ · *why:* The Leibniz
     term and the kinematic term are equal and opposite — whatever u(h) is. · *plain:* the gap fills where more flux enters
     than leaves.
  8. *did:* Insert the profile in the flux · *tex:* $q=\int_0^h\Big[-\frac{h^2}{2\mu}\frac{\partial p}{\partial x}\frac yh\Big(1-\frac yh
     \Big)+U_h\frac yh+U_0\Big(1-\frac yh\Big)\Big]dy$ · *why:* Insert the consistent profile (8.19) of D12, the only place where the momentum balance enters. · *plain:* three pieces to
     integrate.
  9. *did:* Integrate the Poiseuille piece · *tex:* $\int_0^h\frac yh\Big(1-\frac yh\Big)dy=\frac h6\ \Rightarrow\ -\frac{h^3}{12\mu}
     \frac{\partial p}{\partial x}$ · *why:* ∫₀^h (y/h − y²/h²)dy = h/2 − h/3 = h/6, times −h²p_x/(2μ). · *plain:* pressure-
     driven flux ∝ h³.
  10. *did:* Integrate the Couette pieces · *tex:* $\int_0^h\frac yh\,dy=\int_0^h\Big(1-\frac yh\Big)dy=\frac h2\ \Rightarrow\
      \frac{(U_0+U_h)h}{2}$ · *why:* Each linear profile has the mean of its end values. · *plain:* drag flux = mean wall
      speed × gap.
  11. *did:* Add the pieces · *tex:* $q=-\frac{h^3}{12\mu}\frac{\partial p}{\partial x}+\frac{(U_0+U_h)h}{2}$ · *why:* Sum of steps 9
      and 10. · *plain:* the gap flux.
  12. *did:* Combine into one pressure equation · *tex:* $\frac{\partial h}{\partial t}+\frac{\partial}{\partial x}\Big[-\frac{h^3}{12\mu}
      \frac{\partial p}{\partial x}+\frac{(U_0+U_h)h}{2}\Big]=0$ · *why:* Put step 11 into step 7; for a steady gap ∂h/∂t = 0,
      so q is the same at every station. · *plain:* the 1-D Reynolds equation.
- **Result.** $\frac{\partial h}{\partial t}+\frac{\partial q}{\partial x}=0$, $q=-\frac{h^3}{12\mu}\frac{\partial p}{\partial x}+
  \frac{(U_0+U_h)h}{2}$ — *in words:* mass conservation of the gap; steady ⇒ constant flux.
- **Check.** Units of q: m³/(Pa s) × Pa/m = m²/s ✓ and (m/s)m ✓. h constant and ∂p/∂x = 0 → q = (U₀ + U_h)h/2 (Couette).
  Example 8.1's C₁ is this flux measured in the pad frame. sympy (check_src in C06 item 29): Leibniz identity, cancellation
  and final residual all 0.
- **What it means.** One ODE for p(x) in a steady gap: given the flux and two end pressures, the pressure is fixed — the
  lubrication method in one line (C07 applies it).
- **Traps.** The Leibniz boundary term −u(h)h_x cancels against the kinematic v(h) — do not drop only one of them; v(0) = 0
  only for a flat impermeable lower wall.
- **sympy check intent (★★★).** Build u from D12 with generic h(x, t), p_x(x, t); verify the Leibniz identity for this u,
  compute v(h) from continuity (v = −∫₀^y u_x dy), check v(h) − h_t − U_h h_x equals −(h_t + q_x), i.e. the residual of
  step 12 vanishes identically (`ch08.reynolds_equation_sympy`).

### D14 · The slider bearing: $p-p_e=\frac{6\mu LU}{h_o^2}\frac{\alpha(x/L)(1-x/L)}{(2+\alpha)(1+\alpha x/L)^2}$ and $W=\frac{\alpha\mu L^2U}{2h_o^2}$ — ★★★, 15 steps, in C07 (notebook · `slider_bearing`)
- **Goal.** Find the pressure under a sloped pad sliding at U over a flat surface, and the load it carries.
- **Start.** Gap $h(x)=h_o(1+\alpha x/L)$, pad moving at U, floor at rest, p = p_e at both ends; the profile of D12 with
  U₀ = 0, U_h = U; the CV mass balance (4.5).
- **Plan.** (1) Conserve mass in a CV riding with the pad. (2) Turn it into an ODE for p. (3) Integrate with the two end
  pressures. (4) Simplify, then expand for small α and integrate for W.
- **Tools.** moving control volume (Ch. 4 §4.2) · gap flux (D13 steps 9–10) · antiderivatives of (1 + αx/L)⁻ⁿ (gloss) ·
  2 × 2 system · Taylor in α (P26).
- **Assumptions.** Steady in the pad frame (step 1); lubrication profile (step 3); α ≪ 1 used only at step 14.
- **Steps.**
  1. *did:* Balance mass in a pad-fixed CV · *tex:* $\rho B\Big[-\int_0^{h(x_1)}(u-U)dy+\int_0^{h(x_2)}(u-U)dy\Big]=0$ · *why:* (4.5)
     for a CV of fixed shape moving at b = U e_x: its mass is constant, and the flux through a face uses the relative speed
     u − U. · *plain:* what enters the slice leaves it.
  2. *did:* Shrink the slice · *tex:* $\frac{d}{dx}\int_0^{h}(u-U)dy=0\ \Rightarrow\ \int_0^h(u-U)dy=C_1$ · *why:* Divide by ρB(x₂ −
     x₁) and let x₂ → x₁: a derivative of zero means a constant. · *plain:* the pad-frame flux C₁ is the same everywhere.
  3. *did:* Insert the profile's flux · *tex:* $\int_0^hu\,dy=-\frac{h^3}{12\mu}\frac{dp}{dx}+\frac{Uh}{2}$ · *why:* D13 steps 9–11 with
     U₀ = 0, U_h = U. · *plain:* ground-frame flux.
  4. *did:* Subtract Uh · *tex:* $C_1=-\frac{h^3}{12\mu}\frac{dp}{dx}-\frac{Uh}{2}$ · *why:* The pad speed U is the same at every height, so ∫₀^h U dy = Uh; subtracting it gives the pad-frame flux. · *plain:* the book's C₁
     equation.
  5. *did:* Solve for the pressure gradient · *tex:* $\frac{dp}{dx}=-\frac{12\mu C_1}{h^3}-\frac{6\mu U}{h^2}$ · *why:* Rearrange step
     4; this is a first-order ODE for p once h(x) is known. · *plain:* the pressure slope depends on the local gap.
  6. *did:* Insert the gap shape · *tex:* $\frac{dp}{dx}=-\frac{12\mu C_1}{h_o^3}\Big(1+\frac{\alpha x}{L}\Big)^{-3}-\frac{6\mu U}{h_o^2}
     \Big(1+\frac{\alpha x}{L}\Big)^{-2}$ · *why:* h = h₀(1 + αx/L). The book's intermediate line prints (1 − αx/L); the gap
     of the problem needs 1 + αx/L. · *plain:* two powers of the gap.
  7. *did:* Integrate in x · *tex:* $p=\frac{6\mu LC_1}{\alpha h_o^3}\Big(1+\frac{\alpha x}L\Big)^{-2}+\frac{6\mu UL}{\alpha h_o^2}\Big(1+
     \frac{\alpha x}L\Big)^{-1}+C_2$ · *why:* ∫(1 + αx/L)⁻ⁿdx = (L/α)(1 + αx/L)^{1−n}/(1 − n) for n = 3, 2 (substitute s = 1 +
     αx/L). · *plain:* pressure with two unknown constants.
  8. *did:* Apply p(0) = p_e · *tex:* $p_e=\frac{6\mu L}{\alpha h_o^2}\Big[\frac{C_1}{h_o}+U\Big]+C_2$ · *why:* At x = 0 both powers
     equal 1. · *plain:* first end condition.
  9. *did:* Apply p(L) = p_e · *tex:* $p_e=\frac{6\mu L}{\alpha h_o^2}\Big[\frac{C_1}{h_o(1+\alpha)^2}+\frac{U}{1+\alpha}\Big]+C_2$ · *why:*
     At x = L the gap is h₀(1 + α). · *plain:* second end condition.
  10. *did:* Subtract to find C₁ · *tex:* $C_1=-\frac{1+\alpha}{2+\alpha}Uh_o$ · *why:* C₂ cancels: (C₁/h₀)[1 − (1+α)⁻²] +
      U[1 − (1+α)⁻¹] = 0, i.e. (C₁/h₀)α(2+α)/(1+α)² = −Uα/(1+α). · *plain:* the pad-frame flux (negative: toward the
      narrow end).
  11. *did:* Back-substitute for C₂ · *tex:* $C_2=p_e-\frac{6\mu LU}{\alpha h_o^2}\cdot\frac{1}{2+\alpha}$ · *why:* Step 8 with C₁:
      U − (1 + α)U/(2 + α) = U/(2 + α). · *plain:* the book's C₂.
  12. *did:* Put p − p_e over one denominator · *tex:* $p-p_e=\frac{6\mu LU}{\alpha h_o^2}\cdot\frac{s(\alpha-s)}{(2+\alpha)(1+s)^2},\ \
      s=\frac{\alpha x}{L}$ · *why:* The numerator −(1 + α) + (2 + α)(1 + s) − (1 + s)² simplifies to s(α − s). · *plain:*
      the "some algebra" of the book, done.
  13. *did:* Write s back in terms of x · *tex:* $p-p_e=\frac{6\mu LU}{h_o^2}\frac{\alpha(x/L)(1-x/L)}{(2+\alpha)(1+\alpha x/L)^2}$ ·
      *why:* s(α − s) = α²(x/L)(1 − x/L); one α cancels. The printed result has the denominator to the first power; only the
      square satisfies step 5. · *plain:* the exact pressure hump.
  14. *did:* Keep the linear term in α · *tex:* $p-p_e\cong\frac{3\alpha\mu LU}{h_o^2}\frac xL\Big(1-\frac xL\Big)$ · *why:* For α ≪ 1,
      (2 + α)(1 + αx/L)² ≈ 2; (8.19) is only valid for small slopes anyway. · *plain:* a symmetric parabola.
  15. *did:* Integrate for the load · *tex:* $W=\int_0^L(p-p_e)dx=\frac{3\alpha\mu LU}{h_o^2}\cdot\frac L6=\frac{\alpha\mu L^2U}{2h_o^2}$ ·
      *why:* ∫₀^L (x/L)(1 − x/L)dx = L/2 − L/3 = L/6. · *plain:* the load per unit width.
- **Result.** Exact p − p_e above, $p-p_e\cong\frac{3\alpha\mu LU}{h_o^2}\frac xL(1-\frac xL)$, $W=\frac{\alpha\mu L^2U}{2h_o^2}$ —
  *in words:* a constant flux through a narrowing gap needs a pressure hump, and the hump carries the pad.
- **Check.** p(0) = p(L) = p_e ✓ (the printed first-power form also meets them, but it fails step 5's ODE); units μLU/h₀²:
  (Pa s)(m)(m/s)/m² = Pa ✓; W in Pa·m = N/m ✓; αU < 0 ⇒ W < 0. Numbers: L = 5 cm, h₀ = 50 µm, U = 5 m/s, μ = 0.05 Pa s,
  α = 0.1 → W = 12.5 kN/m (linear), 10.8 kN/m (exact, N35).
- **What it means.** Load ∝ 1/h₀²: a heavier load squeezes the film and raises its capacity — the bearing is stable. It works
  only if the pad slides toward the narrow end (αU > 0).
- **Traps.** The pad frame makes the flow steady; the integrands carry (1 + αx/L), not the printed (1 − αx/L); the final
  denominator is squared; α ≪ 1 enters only at step 14.
- **sympy check intent (★★★).** Rebuild steps 5–13 symbolically (`ch08.slider_bearing_sympy`): the squared form satisfies
  dp/dx = −12μC₁/h³ − 6μU/h² and both end conditions; the printed form's ODE residual is non-zero; the exact W(α) of N35
  expands to αμL²U/(2h₀²) + O(α²) (check_src in C07 item 42).

### D15 · The thin-film equation $\frac{\partial h}{\partial t}=\frac{\rho g}{3\mu}\frac{\partial}{\partial x}\big(h^3\frac{\partial h}{\partial x}\big)$ (Example 8.3) — ★★, 10 steps, in C08 (notebook · `viscous_gravity_current`)
- **Goal.** Find one equation for the thickness of a thin viscous layer spreading under its own weight.
- **Start.** The CV mass balance (4.5) on a slice of width dx, and the generic profile $u\cong\frac1\mu\frac{\partial p}{\partial x}
  \frac{y^2}{2}+Ay+B$ (8.18).
- **Plan.** (1) Conserve mass in a slice. (2) Get the pressure gradient from hydrostatics. (3) Fix A and B with no slip below and
  no stress on top. (4) Integrate the profile and insert.
- **Tools.** CV mass (Ch. 4 §4.2) · hydrostatic thin layer and stress-free surface (gloss in C08) · polynomial integral.
- **Assumptions.** Small slope and no inertia (lubrication, steps 3, 5); no surface tension; 2-D; air drag negligible (step 7).
- **Steps.**
  1. *did:* Balance mass in a slice · *tex:* $\rho\frac{\partial h}{\partial t}dx-\int_0^{h(x)}\rho u\,dy+\int_0^{h(x+dx)}\rho u\,dy=0$ ·
     *why:* (4.5) for a fixed slice of width dx and unit depth: storage change plus net outflow is zero. · *plain:* the layer
     thickens where more enters than leaves.
  2. *did:* Divide by ρ dx and shrink · *tex:* $\frac{\partial h}{\partial t}+\frac{\partial}{\partial x}\int_0^hu\,dy=0$ · *why:* The
     difference of the two integrals over dx becomes an x-derivative as dx → 0. · *plain:* thickness changes by flux
     divergence.
  3. *did:* Write the hydrostatic pressure · *tex:* $p=p_a+\rho g(h-y)$ · *why:* In a thin layer with small slope the vertical
     balance is hydrostatic (the lubrication ordering of C05), with the air pressure p_a on top. · *plain:* pressure grows
     with depth below the surface.
  4. *did:* Take its x-derivative · *tex:* $\frac{\partial p}{\partial x}=\rho g\frac{\partial h}{\partial x}$ · *why:* p_a is uniform
     and y does not depend on x. · *plain:* the surface slope drives the flow — same at every height.
  5. *did:* Put it into (8.18) · *tex:* $u\cong\frac{\rho g}{2\mu}\frac{\partial h}{\partial x}y^2+Ay+B$ · *why:* (8.18) with the
     pressure gradient of step 4. · *plain:* a parabola in y.
  6. *did:* Apply no slip on the plate · *tex:* $u(0)=0\ \Rightarrow\ B=0$ · *why:* No slip (8.3) on the stationary surface makes the velocity vanish at y = 0. · *plain:* the
     bottom layer sticks.
  7. *did:* Apply no stress on the surface · *tex:* $\mu\frac{\partial u}{\partial y}\Big\rvert_h=\rho g\frac{\partial h}{\partial x}h+\mu A=0\
     \Rightarrow\ A=-\frac{\rho g}{\mu}h\frac{\partial h}{\partial x}$ · *why:* The air exerts no shear on the free surface. ·
     *plain:* the profile is flat at the top.
  8. *did:* Assemble the profile · *tex:* $u\cong-\frac{\rho g}{2\mu}\frac{\partial h}{\partial x}y(2h-y)$ · *why:* Insert A, B:
     (ρg/2μ)h_x(y² − 2hy). · *plain:* a half-parabola, fastest at the surface, down the slope.
  9. *did:* Integrate across the layer · *tex:* $\int_0^hu\,dy=-\frac{\rho g}{2\mu}\frac{\partial h}{\partial x}\cdot\frac{2h^3}{3}=-\frac
     {\rho g}{3\mu}h^3\frac{\partial h}{\partial x}$ · *why:* ∫₀^h y(2h − y)dy = h³ − h³/3 = 2h³/3. · *plain:* flux ∝ h³ ×
     slope.
  10. *did:* Insert into the mass balance · *tex:* $\frac{\partial h}{\partial t}=\frac{\rho g}{3\mu}\frac{\partial}{\partial x}\Big(h^3
      \frac{\partial h}{\partial x}\Big)$ · *why:* Step 9 into step 2; the two minus signs cancel. · *plain:* a nonlinear
      diffusion equation for h.
- **Result.** $\frac{\partial h}{\partial t}=\frac{\rho g}{3\mu}\frac{\partial}{\partial x}(h^3\frac{\partial h}{\partial x})$ — *in
  words:* a diffusion whose diffusivity ρgh³/3μ collapses as the layer thins.
- **Check.** Units: ρg/3μ is (kg/m³)(m/s²)/(Pa s) = 1/(m s); ∂(h³∂h/∂x)/∂x is m³/m = m²; product m/s = ∂h/∂t ✓. Volume
  ∫h dx is conserved (flux form). Glycerol: ρg/3μ = 4.12 × 10³ m⁻¹s⁻¹.
- **What it means.** The layer spreads ever more slowly; the exponent (t^{1/5}) comes in C10. Lava, mud flows and, with a
  nonlinear rheology, ice sheets obey the same structure.
- **Traps.** The hydrostatic pressure needs a small slope, not only "no acceleration"; the top condition is ∂u/∂y = 0, not
  u = 0; ∫₀^h y(2h − y)dy = 2h³/3; the flux points down the slope.

### D16 · The diffusion equation $\frac{\partial u}{\partial t}=\nu\frac{\partial^2u}{\partial y^2}$ (8.20) for the impulsively started plate — ★, 6 steps, in C09 (notebook · `stokes_first_problem`)
- **Goal.** Reduce Navier–Stokes to a diffusion equation for the fluid above a plate that suddenly starts moving.
- **Start.** Infinite plate y = 0 moving at U for t ≥ 0, fluid at rest initially; (4.10) and (8.1).
- **Plan.** (1) Use invariance in x. (2) Show the pressure gradient vanishes. (3) Divide by ρ.
- **Tools.** continuity (D02's move) · partial derivative (P25).
- **Assumptions.** Infinite plate (nothing depends on x); fluid at rest at infinity at every finite t (step 5).
- **Steps.**
  1. *did:* Use invariance along the plate · *tex:* $\frac{\partial}{\partial x}=0\ \Rightarrow\ \frac{\partial v}{\partial y}=0\ \Rightarrow\
     v\equiv0$ · *why:* Continuity with ∂u/∂x = 0, and v = 0 at the impermeable plate (the D02 move). · *plain:* the flow is
     parallel to the plate.
  2. *did:* Drop the advective terms · *tex:* $u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=0$ · *why:* Both factors
     vanish; the equation becomes linear although u is not small. · *plain:* no advection.
  3. *did:* Write the two momentum equations · *tex:* $\rho\frac{\partial u}{\partial t}=-\frac{\partial p}{\partial x}+\mu\frac{\partial^2u}
     {\partial y^2},\ \ 0=-\frac{\partial p}{\partial y}$ · *why:* The x- and y-components of (8.1) times ρ with steps 1–2. ·
     *plain:* unsteady friction balance.
  4. *did:* Read the y-equation · *tex:* $p=p(x,t)$ · *why:* ∂p/∂y = 0: the x-gradient of pressure is the same at every
     height. · *plain:* what holds far away holds at the plate.
  5. *did:* Evaluate far from the plate · *tex:* $\frac{\partial p}{\partial x}=0$ · *why:* At any finite t there is a height where
     the fluid is still at rest, so ∂u/∂t = ∂²u/∂y² = 0 there, forcing ∂p/∂x = 0 — and by step 4 everywhere. · *plain:* no
     pressure gradient drives this flow.
  6. *did:* Divide by ρ · *tex:* $\frac{\partial u}{\partial t}=\nu\frac{\partial^2u}{\partial y^2}$ (8.20) · *why:* Dividing by the constant ρ turns μ into the kinematic viscosity ν = μ/ρ, the diffusivity of C01. · *plain:*
     (8.20): velocity diffuses like heat.
- **Result.** (8.20) with (8.21)–(8.23) — *in words:* the plate's motion spreads upward by pure diffusion.
- **Check.** Units: m/s² both sides ✓; the same equation as Ch. 1's Couette start-up (R14) without the upper wall.
- **What it means.** The problem is linear and has no length scale — the door to the similarity solution (D17).
- **Traps.** ∂p/∂x = 0 because a still region exists far away at every finite t; linearity comes from u∂u/∂x = 0, not from
  small u.

### D17 · From dimensional analysis to one variable: $u/U=f(y/\sqrt{\nu t},y/Ut)$ (8.24) → $u/U=F(y/\sqrt{\nu t})$ (8.25) — ★★, 8 steps, in C09 (notebook · `stokes_first_problem`)
- **Goal.** Show that the velocity can depend on y and t only through the single combination y/√(νt).
- **Start.** u = u(U, y, t, ν) solving (8.20)–(8.23).
- **Plan.** (1) Count variables and dimensions. (2) Form the groups. (3) Use linearity to drop one.
- **Tools.** Buckingham Π (Ch. 1 §1.11, `ch01.pi_groups`) · linearity (gloss).
- **Assumptions.** Linear equation and conditions (step 5); no other parameters (no imposed length or time).
- **Steps.**
  1. *did:* List the variables · *tex:* $u,\ U,\ y,\ t,\ \nu$ · *why:* The answer can depend only on what appears in (8.20)–(8.23):
     the plate speed, height, time and viscosity. · *plain:* five quantities.
  2. *did:* Count the dimensions · *tex:* $[u]=[U]=\mathrm{L/T},\ [y]=\mathrm L,\ [t]=\mathrm T,\ [\nu]=\mathrm{L^2/T}$ · *why:* Only
     length and time appear, so there are two base dimensions. · *plain:* two dimensions.
  3. *did:* Apply the Π theorem · *tex:* $5-2=3\ \text{groups}:\ \ \frac uU,\ \ \frac{y}{\sqrt{\nu t}},\ \ \frac{y}{Ut}$ · *why:* Buckingham:
     variables minus independent dimensions; each group is checked dimensionless. · *plain:* three dimensionless groups.
  4. *did:* Write the functional relation · *tex:* $u/U=f(y/\sqrt{\nu t},\ y/Ut)$ (8.24) · *why:* One group as a function of the
     other two. · *plain:* (8.24).
  5. *did:* Use linearity · *tex:* $u\ \text{solves with }U\ \Rightarrow\ \lambda u\ \text{solves with }\lambda U$ · *why:* (8.20) and its
     conditions are linear; only the wall value contains U, and it scales with it. · *plain:* doubling U doubles u.
  6. *did:* Conclude u/U cannot depend on U · *tex:* $\frac{\partial}{\partial U}\Big(\frac uU\Big)=0$ · *why:* By step 5, u/U is the
     same for every U; but y/Ut contains U, so f cannot depend on its second argument. · *plain:* the second group must go.
  7. *did:* Drop the second group · *tex:* $u/U=F(y/\sqrt{\nu t})\equiv F(\eta)$ (8.25) · *why:* What remains is a function of the
     first group only. · *plain:* (8.25): one similarity variable.
  8. *did:* Name the similarity variable · *tex:* $\eta=\frac{y}{\sqrt{\nu t}}$ · *why:* At fixed t, η ∝ y, so F reads as a profile;
     another choice (νt/y²) would only change F, not u. The book's figures use η/2. · *plain:* two variables became one.
- **Result.** (8.25) — *in words:* the PDE in y and t becomes an ODE in η.
- **Check.** η dimensionless: m/√((m²/s)s) ✓. Numbers: y = 1 cm, t = 100 s, ν = 10⁻⁶ → η = 1.
- **What it means.** No imposed scale exists, so y can only be measured in units of √(νt).
- **Traps.** Five variables and two dimensions give three groups; it is linearity, not dimensional analysis, that removes
  y/Ut.

### D18 · The similarity ODE $-\frac\eta2\frac{dF}{d\eta}=\frac{d}{d\eta}\big(\frac{dF}{d\eta}\big)$ (8.26) with $F(0)=1$ (8.27), $F(\infty)=0$ (8.28) — ★★, 9 steps, in C09 (notebook · `stokes_first_problem`)
- **Goal.** Turn (8.20) into an ODE for F(η) and its conditions.
- **Start.** u = U F(η), η = y/√(νt), in $\frac{\partial u}{\partial t}=\nu\frac{\partial^2u}{\partial y^2}$ (8.20).
- **Plan.** (1) Differentiate η. (2) Chain rule for ∂u/∂t and ∂²u/∂y². (3) Substitute, divide, check that t drops out.
  (4) Translate the conditions.
- **Tools.** chain rule (P49/P91) · exponent rules (P43).
- **Assumptions.** t > 0.
- **Steps.**
  1. *did:* Differentiate η in time · *tex:* $\frac{\partial\eta}{\partial t}=-\frac12\frac{y}{\sqrt\nu}t^{-3/2}=-\frac{\eta}{2t}$ · *why:*
     η = yν^{−1/2}t^{−1/2}; power rule at fixed y. · *plain:* η shrinks in time at a fixed height.
  2. *did:* Chain rule for ∂u/∂t · *tex:* $\frac{\partial u}{\partial t}=U\frac{dF}{d\eta}\frac{\partial\eta}{\partial t}=-\frac{U\eta}{2t}
     \frac{dF}{d\eta}$ · *why:* u depends on t only through η. · *plain:* the local rate of change.
  3. *did:* Differentiate η in y · *tex:* $\frac{\partial\eta}{\partial y}=\frac{1}{\sqrt{\nu t}}$ · *why:* η is linear in y with the slope 1/√(νt), which does not depend on y. ·
     *plain:* one unit of η is √(νt) of height.
  4. *did:* First y-derivative · *tex:* $\frac{\partial u}{\partial y}=\frac{U}{\sqrt{\nu t}}\frac{dF}{d\eta}$ · *why:* Chain rule: u depends on y only through η, so ∂u/∂y = U F′ ∂η/∂y. ·
     *plain:* the velocity slope.
  5. *did:* Second y-derivative · *tex:* $\frac{\partial^2u}{\partial y^2}=\frac{U}{\nu t}\frac{d^2F}{d\eta^2}$ · *why:* 1/√(νt) does
     not depend on y; differentiate F′ again with the chain rule. · *plain:* the curvature.
  6. *did:* Substitute into (8.20) · *tex:* $-\frac{U\eta}{2t}F'=\nu\frac{U}{\nu t}F''=\frac Ut F''$ · *why:* Insert steps 2 and 5 into (8.20) to see whether the time dependence can cancel. ·
     *plain:* every term carries U/t.
  7. *did:* Divide by U/t · *tex:* $-\frac\eta2\frac{dF}{d\eta}=\frac{d}{d\eta}\Big(\frac{dF}{d\eta}\Big)$ (8.26) · *why:* U/t ≠ 0; t
     disappears — the test that the similarity form is right. · *plain:* (8.26), an ODE.
  8. *did:* Translate the wall condition · *tex:* $F(\eta=0)=1$ (8.27) · *why:* y = 0 gives η = 0 for every t > 0, and u = U
     there (8.22). · *plain:* fluid at the plate moves with it.
  9. *did:* Merge the other two conditions · *tex:* $F(\eta\to\infty)=0$ (8.28) · *why:* y → ∞ at fixed t (8.23) and t → 0⁺ at
     fixed y > 0 (8.21) both send η → ∞; three conditions become two, as a second-order ODE needs. · *plain:* far away and
     initially, the fluid is at rest.
- **Result.** (8.26), (8.27), (8.28) — *in words:* a two-point boundary-value problem in η alone.
- **Check.** `ch08.similarity_ode_solve("stokes1")` (solve_bvp) finds F with F(0) = 1, F(12) = 0 matching erfc(η/2) to
  1e-7.
- **What it means.** Had t not cancelled at step 7, the guess u = UF(y/√(νt)) would have been wrong.
- **Traps.** ∂η/∂t = −η/(2t), not −η/t; t → 0 and y → ∞ both give η → ∞.

### D19 · Solving the ODE: (8.29) and $\frac uU=1-\mathrm{erf}\big(\frac{y}{2\sqrt{\nu t}}\big)$ (8.30) — ★★, 11 steps, in C09 (notebook · `stokes_first_problem`)
- **Goal.** Solve (8.26) with its conditions and express the answer through the error function.
- **Start.** $-\frac\eta2F'=F''$ (8.26), F(0) = 1 (8.27), F(∞) = 0 (8.28).
- **Plan.** (1) Solve for F′ by separation. (2) Integrate again. (3) Fix the constants with the Gaussian integral.
  (4) Recognise erf.
- **Tools.** separation of variables (P42) · natural log and exponential (P36) · Gaussian integral (primer P194) ·
  substitution (P106) · erf (primer P195).
- **Assumptions.** None beyond (8.26)–(8.28).
- **Steps.**
  1. *did:* Call the slope G · *tex:* $G\equiv F'\ \Rightarrow\ -\frac\eta2G=\frac{dG}{d\eta}$ · *why:* F itself does not appear in
     (8.26), so it is a first-order ODE for F′. · *plain:* first solve for the slope.
  2. *did:* Separate the variables · *tex:* $-\frac\eta2d\eta=\frac{dG}{G}$ · *why:* Divide by G and multiply by dη (P42). ·
     *plain:* η on one side, G on the other.
  3. *did:* Integrate both sides · *tex:* $-\frac{\eta^2}{4}=\ln G+\text{const}$ · *why:* ∫η dη = η²/2; ∫dG/G = ln G. · *plain:* a
     logarithm of the slope.
  4. *did:* Exponentiate · *tex:* $\frac{dF}{d\eta}=A\exp(-\eta^2/4)$ · *why:* e^{ln G} = G; the constant becomes a factor A. ·
     *plain:* the slope is a Gaussian.
  5. *did:* Integrate from 0 to η · *tex:* $F(\eta)=A\int_0^\eta\exp(-\xi^2/4)\,d\xi+B$ (8.29) · *why:* Fundamental theorem of calculus;
     ξ is a dummy variable, B = F(0). · *plain:* (8.29).
  6. *did:* Apply F(0) = 1 · *tex:* $B=1$ · *why:* The integral from 0 to 0 vanishes. · *plain:* the plate condition.
  7. *did:* Apply F(∞) = 0 · *tex:* $0=A\int_0^\infty e^{-\xi^2/4}d\xi+1$ · *why:* Apply (8.28), F → 0 as η → ∞, with B = 1 already fixed. · *plain:* one equation for A.
  8. *did:* Substitute ξ = 2ζ · *tex:* $0=2A\int_0^\infty e^{-\zeta^2}d\zeta+1$ · *why:* dξ = 2dζ and ξ²/4 = ζ²; this turns the
     integral into the standard Gaussian. · *plain:* the factor 2 that ends up in y/(2√(νt)).
  9. *did:* Use the Gaussian integral · *tex:* $0=2A\frac{\sqrt\pi}{2}+1\ \Rightarrow\ A=-\frac{1}{\sqrt\pi}$ · *why:* ∫₀^∞e^{−ζ²}dζ =
     √π/2, half of the full √π by symmetry (P194). · *plain:* A is negative: F falls.
  10. *did:* Substitute again inside F · *tex:* $F(\eta)=1-\frac{2}{\sqrt\pi}\int_0^{\eta/2}e^{-\zeta^2}d\zeta$ · *why:* Same ξ = 2ζ;
      the upper limit becomes η/2. · *plain:* ready to name.
  11. *did:* Recognise the error function · *tex:* $\frac{u}{U}=1-\mathrm{erf}\Big(\frac{y}{2\sqrt{\nu t}}\Big),\ \ \mathrm{erf}(\zeta)=
      \frac{2}{\sqrt\pi}\int_0^\zeta e^{-\xi^2}d\xi$ (8.30) · *why:* The definition of erf (P195); η/2 = y/(2√(νt)). ·
      *plain:* (8.30) = erfc.
- **Result.** (8.30) — *in words:* the velocity falls from U at the plate as the complementary error function of
  y/(2√(νt)).
- **Check.** F(0) = 1 − 0 = 1 ✓, F(∞) = 1 − 1 = 0 ✓; residual of (8.20) zero (sympy); y = 1 cm, t = 100 s: erfc(0.5) =
  0.4795; the CN solution converges to it at order 2.
- **What it means.** One curve for all t, U, ν — the collapse of Fig. 8.13; the plate's influence reaches ~ 3.6√(νt) (D20).
- **Traps.** The substitution carries dξ = 2dζ; ∫₀^∞ is √π/2 (half); A < 0; compute with erfc, not 1 − erf, for large η.

### D20 · The 99 % thickness $\delta_{99}\sim3.64\sqrt{\nu t}$ (8.31) — ★, 5 steps, in C09 (notebook · `stokes_first_problem`)
- **Goal.** Compute (instead of reading off a figure) the height where u falls to 1 % of U.
- **Start.** $\frac uU=1-\mathrm{erf}\big(\frac{y}{2\sqrt{\nu t}}\big)=\mathrm{erfc}\big(\frac{\eta}{2}\big)$ (8.30).
- **Plan.** (1) Set the level. (2) Invert erfc. (3) Convert to height.
- **Tools.** erf and its inverses (primer P195) · brentq (P108) as a cross-check.
- **Assumptions.** The 1 % level is a convention.
- **Steps.**
  1. *did:* Choose the edge level · *tex:* $\frac uU=0.01$ · *why:* The profile never reaches zero; the book (and Ch. 9) define
     the edge at 99 % of the change. · *plain:* the 1 % point.
  2. *did:* Write it with (8.30) · *tex:* $\mathrm{erfc}\Big(\frac{\eta_{99}}{2}\Big)=0.01$ · *why:* Substitute the 1 % level into (8.30), written with η/2 as its argument.
     · *plain:* one equation for η₉₉.
  3. *did:* Invert erfc · *tex:* $\frac{\eta_{99}}{2}=\mathrm{erfc}^{-1}(0.01)=1.8214$ · *why:* erfc is monotonic, so it has a
     unique inverse (`special.erfcinv`; brentq agrees). · *plain:* on the book's axis η/2 = 1.82.
  4. *did:* Double it · *tex:* $\eta_{99}=3.643\ \Rightarrow\ \delta_{99}=3.643\sqrt{\nu t}$ (8.31) · *why:* η = y/√(νt), so the height
     is η₉₉√(νt). · *plain:* (8.31).
  5. *did:* Put in numbers · *tex:* $\nu=10^{-6},\ t=100\ \mathrm s:\ \ \delta_{99}=3.64\ \mathrm{cm};\ \ t=3600\ \mathrm s:\ 21.9\ \mathrm{cm}$ ·
     *why:* √(νt) = 1 cm and 6 cm; the thickness grows as t^{1/2}. · *plain:* an hour of diffusion reaches only a hand's
     width in water.
- **Result.** δ₉₉ = 2 erfc⁻¹(0.01)√(νt) = 3.643√(νt) — *in words:* the moving layer thickens like √(νt).
- **Check.** Units m ✓; a 95 % level gives 2 erfc⁻¹(0.05) = 2.772 — the factor is a convention.
- **What it means.** The √(νt) law is the thickness of every laminar boundary layer (Ch. 9, with t = x/U) and of spin-up
  layers in the ocean (Ch. 13).
- **Traps.** The book's figure axis is η/2 = 1.82; (8.31) quotes η = 3.64.

### D21 · Example 8.4: the ansatz with A = 1, n = 0 gives δ = √(2C₁νt) — ★★, 8 steps, in C10 (notebook · `similarity_exponents`)
- **Goal.** Recover Stokes' first problem from the general similarity form, finding the growth law δ(t) instead of
  assuming it.
- **Start.** $\gamma=At^{-n}F(\xi/\delta(t))$ (8.32a) with γ = u/U, ξ = y, and (8.20).
- **Plan.** (1) Fix A and n from the wall. (2) Differentiate. (3) Demand that the brackets scale alike. (4) Solve for δ.
- **Tools.** exponent matching (primer P197) · chain rule (P49) · separable ODE (P42).
- **Assumptions.** A similarity solution exists (no imposed scales).
- **Steps.**
  1. *did:* Fix the amplitude law · *tex:* $A=1,\ n=0:\ \ \frac uU=F\Big(\frac{y}{\delta(t)}\Big)$ · *why:* u/U = 1 at η = 0 for every
     t > 0; any factor t^{−n} with n ≠ 0 would change the wall speed. · *plain:* only the width can change.
  2. *did:* Differentiate in time · *tex:* $\frac{\partial u}{\partial t}=U\frac{dF}{d\eta}\Big(-\frac{y}{\delta^2}\Big)\frac{d\delta}{dt}=-U\eta
     \frac{\delta'}{\delta}F'$ · *why:* Chain rule through η = y/δ(t). · *plain:* the profile stretches.
  3. *did:* Differentiate twice in y · *tex:* $\frac{\partial^2u}{\partial y^2}=\frac{U}{\delta^2}F''$ · *why:* Chain rule with ∂η/∂y = 1/δ, applied twice; δ does not depend on y. · *plain:*
     curvature in η units.
  4. *did:* Substitute into (8.20) · *tex:* $-U\eta\frac{\delta'}{\delta}F'=\nu\frac{U}{\delta^2}F''$ · *why:* Insert steps 2 and 3 into (8.20) to see which factors depend on time. · *plain:*
     one equation with t-dependent coefficients.
  5. *did:* Divide by U and bracket · *tex:* $-\Big[\frac{1}{\delta}\frac{d\delta}{dt}\Big]\eta\frac{dF}{d\eta}=\Big[\frac{\nu}{\delta^2}\Big]
     \frac{d^2F}{d\eta^2}$ · *why:* Put everything that depends on t into brackets. · *plain:* the book's bracketed form.
  6. *did:* Make the brackets proportional · *tex:* $\frac1\delta\frac{d\delta}{dt}=C_1\frac{\nu}{\delta^2}$ · *why:* Exponent matching
     (P197): dividing by one bracket must leave an ODE in η alone, so their ratio is a constant C₁. · *plain:* both sides
     scale alike in time.
  7. *did:* Solve for δ · *tex:* $\delta\frac{d\delta}{dt}=C_1\nu\ \Rightarrow\ \frac{\delta^2}{2}=C_1\nu t+C_2$ · *why:* Multiply by δ²;
     δδ′ = d(δ²/2)/dt, integrate in t. · *plain:* δ² grows linearly in time.
  8. *did:* Start from zero thickness · *tex:* $\delta(0)=0\ \Rightarrow\ \delta=\sqrt{2C_1\nu t};\ \ C_1=\tfrac12:\ \delta=\sqrt{\nu t}$ ·
     *why:* The layer starts with zero thickness; C₁ = ½ reproduces η of (8.25) and turns step 5 into (8.26). · *plain:* the
     √(νt) law found, not assumed.
- **Result.** δ = √(2C₁νt); with C₁ = ½, $-\frac\eta2F'=F''$ (8.26) — *in words:* the equation itself forces δ ∝ √t.
- **Check.** Units: δδ′ = C₁ν in m²/s ✓. `similarity_reduce_sympy("stokes1_delta")` returns δ = √(2C₁νt).
- **What it means.** The recipe: guess, bracket, match powers. D22 and D23 add a conserved quantity for the second exponent.
- **Traps.** ∂F(y/δ)/∂t = −ηF′δ′/δ; δ(0) = 0 fixes C₂ = 0; A = 1, n = 0 come from the wall condition, not the ODE.

### D22 · Example 8.5: the viscous vortex sheet, n = ½, $\omega_z=-\frac{U}{\sqrt{\pi\nu t}}e^{-y^2/4\nu t}$, $u=U\,\mathrm{erf}\frac{y}{2\sqrt{\nu t}}$ — ★★★, 14 steps, in C10 (notebook · `similarity_exponents`)
- **Goal.** Find how a velocity jump (+U above, −U below) thickens, using the ansatz with an amplitude that decays in time.
- **Start.** ω_z = −∂u/∂y, u = ±U for y ≷ 0 at t = 0, and (8.20).
- **Plan.** (1) Get the vorticity equation. (2) Insert the ansatz, bracket, take δ = √(νt). (3) Fix n by the conserved jump.
  (4) Integrate the ODE with an exact derivative. (5) Fix the amplitude; integrate for u.
- **Tools.** exponent matching (primer P197) · product rule read backwards (gloss) · Gaussian integral (primer P194) · erf
  (primer P195).
- **Assumptions.** Similarity (no imposed scales); ω and ηF vanish far away (step 11).
- **Steps.**
  1. *did:* Differentiate (8.20) with respect to y · *tex:* $\frac{\partial\omega_z}{\partial t}=\nu\frac{\partial^2\omega_z}{\partial y^2}$ ·
     *why:* ω_z = −∂u/∂y; y- and t-derivatives commute (P121). The initial jump is a delta of vorticity. · *plain:* vorticity
     diffuses too.
  2. *did:* Insert the ansatz · *tex:* $\omega_z=At^{-n}F(\eta),\ \ \eta=\frac{y}{\delta(t)}$ · *why:* (8.32a): y does not appear in the
     initial condition, so a time power is the right prefactor; the peak must fall as the sheet spreads. · *plain:* an
     amplitude and a width that change in time.
  3. *did:* Differentiate in time · *tex:* $\frac{\partial\omega_z}{\partial t}=-nAt^{-n-1}F-At^{-n}\eta\frac{\delta'}{\delta}F'$ · *why:*
     Product rule on t^{−n}F, chain rule through η. · *plain:* amplitude decay plus stretching.
  4. *did:* Differentiate twice in y · *tex:* $\frac{\partial^2\omega_z}{\partial y^2}=At^{-n}\frac{F''}{\delta^2}$ · *why:* Chain rule with ∂η/∂y = 1/δ, applied twice; the prefactor At^(−n) does not depend on y. · *plain:* curvature.
  5. *did:* Substitute and divide by At^{−n} · *tex:* $-\Big[\frac nt\Big]F-\Big[\frac{\delta'}{\delta}\Big]\eta F'=\Big[\frac{\nu}{\delta^2}\Big]
     F''$ · *why:* Insert steps 3 and 4 into the vorticity equation of step 1 and divide by At^(−n). · *plain:* three brackets.
  6. *did:* Take δ = √(νt) and multiply by t · *tex:* $-nF-\frac12\eta F'=F''$ · *why:* From D21 (C₁ = ½): δ′/δ = 1/(2t),
     ν/δ² = 1/t — every bracket ∝ 1/t. · *plain:* an ODE with n still free.
  7. *did:* State the conserved jump · *tex:* $-\int_{-\infty}^{\infty}\omega_z\,dy=u(\infty)-u(-\infty)=2U$ · *why:* ω_z = −∂u/∂y
     integrates to the velocity difference, which diffusion cannot change. · *plain:* the total vorticity is fixed.
  8. *did:* Insert the ansatz into it · *tex:* $-At^{-n}\delta\int_{-\infty}^{\infty}F(\eta)\,d\eta=2U$ · *why:* Change the integration variable to η: y = δη, so dy = δ dη at fixed t. · *plain:*
     amplitude × width × a number.
  9. *did:* Match the powers of t · *tex:* $t^{-n}\delta\propto t^{-n+1/2}=\text{const}\ \Rightarrow\ n=\tfrac12$ · *why:* ∫F dη is a
     number and the right side is constant for all t (P197). · *plain:* the peak falls as t^{−1/2}.
  10. *did:* Spot an exact derivative · *tex:* $-\frac12\Big(F+\eta\frac{dF}{d\eta}\Big)=-\frac12\frac{d}{d\eta}(\eta F)=\frac{d}{d\eta}\Big(\frac
      {dF}{d\eta}\Big)$ · *why:* Step 6 with n = ½; F + ηF′ = (ηF)′ by the product rule read backwards. · *plain:* both
      sides are derivatives.
  11. *did:* Integrate once · *tex:* $\frac{dF}{d\eta}+\frac12\eta F=C=0$ · *why:* F, F′ and ηF vanish as η → ∞ (checked once F
      is known), so the constant is zero. · *plain:* a first-order ODE.
  12. *did:* Separate and integrate · *tex:* $\frac{F'}{F}=-\frac\eta2\ \Rightarrow\ F=De^{-\eta^2/4}$ · *why:* ln F = −η²/4 + const;
      ηF → 0 indeed. · *plain:* a Gaussian.
  13. *did:* Fix the amplitude product AD · *tex:* $-At^{-1/2}\sqrt{\nu t}\,D\cdot2\sqrt\pi=2U\ \Rightarrow\ AD=-\frac{U}{\sqrt{\pi\nu}}$ ·
      *why:* ∫e^{−η²/4}dη = 2√π (Gaussian with η = 2ζ) in step 8 with δ = √(νt). · *plain:* the sign follows the jump.
  14. *did:* Assemble ω and integrate for u · *tex:* $\omega_z=-\frac{U}{\sqrt{\pi\nu t}}e^{-y^2/4\nu t},\ \ u=-\int_0^y\omega_z\,dy'=
      U\,\mathrm{erf}\frac{y}{2\sqrt{\nu t}}$ · *why:* u(0) = 0 by odd symmetry; the integral is the erf of P195. · *plain:*
      the sheet thickens like √(νt).
- **Result.** $\omega_z=-\frac{U}{\sqrt{\pi\nu t}}e^{-y^2/4\nu t}$, $u=U\,\mathrm{erf}\frac{y}{2\sqrt{\nu t}}$ — *in words:* a Gaussian
  vorticity layer of fixed total strength, widening as √(νt).
- **Check.** −∫ω dy = 2U at every t ✓ (sympy); u → ±U far away ✓; ±0.95U at η = ±2 erfinv(0.95) = ±2.772 (the book prints
  2.76) and width 5.544√(νt); U = 1 cm/s, t = 1 s: width 5.5 mm, peak ω −5.64 s⁻¹.
- **What it means.** The upper half is Stokes' first problem seen from a moving frame (N60): a temporally developing
  boundary layer with C_f ∝ Re_x^{−1/2}.
- **Traps.** It is the conserved jump, not the ODE, that fixes n = ½; C = 0 because ηF → 0; AD is negative; 2.772, not
  2.76.
- **sympy check intent (★★★).** `similarity_reduce_sympy("vortex_sheet")` reproduces the brackets and n = ½ and checks the
  final ω solves the diffusion equation; the check cell also integrates ω symbolically over all y to show the jump 2U
  (check_src in C10 item 39).

### D23 · Example 8.7: the spreading bead, n = m = 1/5, $h=At^{-1/5}F(x/Dt^{1/5})$ — ★★, 9 steps, in C10 (notebook · `viscous_gravity_current` · `similarity_exponents`)
- **Goal.** Find the similarity form of the spreading bead of C08 from the equation and its conserved volume.
- **Start.** $\frac{\partial h}{\partial t}=\frac{\rho g}{3\mu}\frac{\partial}{\partial x}\big(h^3\frac{\partial h}{\partial x}\big)$ (Example 8.3)
  and the ansatz $h=At^{-n}F(x/\delta(t))$ (8.32a).
- **Plan.** (1) Compute every derivative. (2) Bracket the t-dependence. (3) Two proportionality conditions and the volume.
  (4) Solve the exponent equations.
- **Tools.** exponent matching (primer P197) · product and chain rules (P38, P49).
- **Assumptions.** Constant volume (no source); similarity (step 1).
- **Steps.**
  1. *did:* Differentiate in time · *tex:* $\frac{\partial h}{\partial t}=-nAt^{-n-1}F-At^{-n}\eta\frac{\delta'}{\delta}F'$ · *why:*
     Product rule on t^{−n}F and chain rule through η = x/δ(t), as in D22 step 3. · *plain:* sinking plus stretching.
  2. *did:* Differentiate in space · *tex:* $\frac{\partial h}{\partial x}=At^{-n}\frac{F'}{\delta}$ · *why:* Chain rule with ∂η/∂x = 1/δ; the prefactor At^(−n) does not depend on x. · *plain:* the
     slope.
  3. *did:* Form the flux factor · *tex:* $h^3\frac{\partial h}{\partial x}=A^4t^{-4n}\frac{F^3F'}{\delta}$ · *why:* Multiply the cube of
     the ansatz by step 2. · *plain:* h³ brings A³t^{−3n}.
  4. *did:* Differentiate the flux factor · *tex:* $\frac{\partial}{\partial x}\Big(h^3\frac{\partial h}{\partial x}\Big)=\frac{A^4t^{-4n}}
     {\delta^2}(3F^2F'^2+F^3F'')$ · *why:* (F³F′)′ = 3F²F′² + F³F″ by the product rule; one more 1/δ. · *plain:* the diffusion
     side.
  5. *did:* Assemble with brackets · *tex:* $-[nAt^{-n-1}]F-\Big[At^{-n}\frac{\delta'}{\delta}\Big]\eta F'=\Big[\frac{\rho g}{3\mu}\frac{A^4
     t^{-4n}}{\delta^2}\Big](3F^2F'^2+F^3F'')$ · *why:* Steps 1 and 4 in the equation. · *plain:* three bracketed coefficients.
  6. *did:* Match the first two brackets · *tex:* $\frac{\delta'}{\delta}=\frac{Cn}{t}\ \Rightarrow\ \delta=Dt^m,\ \ m=Cn$ · *why:* Their
     ratio must be constant (P197); δ′/δ ∝ 1/t integrates to a power law. · *plain:* the width is a power of time.
  7. *did:* Match the second and third · *tex:* $t^{-n-1}\propto t^{-4n}\delta^{-2}=t^{-4n-2m}\ \Rightarrow\ -1=-3n-2m$ · *why:* Equal
     powers of t on both sides (the constant factors are absorbed into F). · *plain:* first exponent relation.
  8. *did:* Conserve the volume · *tex:* $\int_{-\infty}^{\infty}h\,dx=At^{-n}Dt^m\int F\,d\eta=\text{const}\ \Rightarrow\ -n+m=0$ · *why:*
     dx = δ dη; ∫F dη is a number, so the power of t must vanish. · *plain:* second relation — from a conserved quantity,
     not a boundary condition.
  9. *did:* Solve the two relations · *tex:* $n=m,\ \ 3n+2n=1\ \Rightarrow\ n=m=\tfrac15:\ \ h=At^{-1/5}F\big(x/Dt^{1/5}\big)$ · *why:* Two
     linear equations in n, m. · *plain:* height falls and width grows like t^{1/5}.
- **Result.** $h(x,t)=At^{-1/5}F(x/Dt^{1/5})$ — *in words:* the bead spreads as t^{1/5} and thins as t^{−1/5}.
- **Check.** Volume ∝ t^{−1/5}t^{1/5} = const ✓. The C08 run's fitted front slope is 0.200 ± 0.005; Huppert's closed form
  (N63) has the same exponents; `similarity_reduce_sympy("spreading")` returns n = m = 1/5.
- **What it means.** Not a diffusion length (t^{1/2}): the diffusivity ρgh³/3μ collapses as the layer thins. The book's last
  line says 'the final equation of Example 8.2'; it means Example 8.3.
- **Traps.** (h³h_x)_x brings δ^{−2} and A⁴t^{−4n}; the second condition is the volume, not a boundary condition.

### D24 · Stokes' second problem: (8.35) → (8.36) → (8.37) → $u=Ue^{-y\sqrt{\omega/2\nu}}\cos(\omega t-y\sqrt{\omega/2\nu})$ (8.38) — ★★, 10 steps, in C11 (notebook · `oscillating_plate`)
- **Goal.** Find the periodic flow above a plate oscillating as U cos ωt.
- **Start.** $\frac{\partial u}{\partial t}=\nu\frac{\partial^2u}{\partial y^2}$ (8.20), $u(y=0,t)=U\cos(\omega t)$ (8.33), u bounded as y →
  ∞ (8.34).
- **Plan.** (1) Complex form. (2) ODE in y. (3) Exponential roots with √i. (4) Pick the bounded root, fix the amplitude,
  take the real part.
- **Tools.** complex exponential solution (gloss, P176) · complex square roots (P159) · linear ODE and e^{λy} (P44).
- **Assumptions.** Periodic steady state (transients gone) — step 1.
- **Steps.**
  1. *did:* Try the complex form · *tex:* $u(y,t)=\mathrm{Re}\{e^{i\omega t}f(y)\}$ (8.35) · *why:* The wall motion is Re{Ue^{iωt}};
     (8.20) is linear with real coefficients, so we solve for ũ = e^{iωt}f and take the real part at the end. · *plain:* every
     layer oscillates at the wall's frequency.
  2. *did:* Differentiate ũ · *tex:* $\frac{\partial\tilde u}{\partial t}=i\omega e^{i\omega t}f,\ \ \frac{\partial^2\tilde u}{\partial y^2}=e^{i\omega t}
     f''$ · *why:* d/dt e^{iωt} = iωe^{iωt}; y enters only through f. · *plain:* time derivative → multiplication by iω.
  3. *did:* Substitute and cancel e^{iωt} · *tex:* $i\omega f=\nu\frac{d^2f}{dy^2}$ (8.36) · *why:* e^{iωt} ≠ 0. The text says "substitution
     of (8.33)" — it is (8.35) that is substituted. · *plain:* (8.36), an ODE with constant coefficients.
  4. *did:* Try f = e^{ky} · *tex:* $i\omega=\nu k^2\ \Rightarrow\ k^2=\frac{i\omega}{\nu}$ · *why:* The standard trial for constant
     coefficients (P44). · *plain:* k is a complex square root.
  5. *did:* Take the square root of i · *tex:* $k=\pm\frac{1+i}{\sqrt2}\sqrt{\frac\omega\nu}=\pm(1+i)\sqrt{\frac{\omega}{2\nu}}$ · *why:*
     ((1 + i)/√2)² = (1 + 2i − 1)/2 = i (P159). · *plain:* decay and oscillation at the same rate.
  6. *did:* Write the general solution · *tex:* $f=A\exp\{-(i+1)y\sqrt{\omega/2\nu}\}+B\exp\{+(i+1)y\sqrt{\omega/2\nu}\}$ (8.37) · *why:* The ODE is linear, so any combination of the two exponential solutions is a solution. · *plain:* (8.37).
  7. *did:* Discard the growing root · *tex:* $B=0$ · *why:* The second exponential grows like e^{+y√(ω/2ν)}, violating (8.34). ·
     *plain:* only the decaying layer survives.
  8. *did:* Match the wall motion · *tex:* $\tilde u(0,t)=Ae^{i\omega t}\ \Rightarrow\ A=U$ · *why:* Re{Ae^{iωt}} must be U cos ωt for all
     t. · *plain:* amplitude fixed.
  9. *did:* Name the depth · *tex:* $\tilde u=Ue^{-y/\delta_e}e^{i(\omega t-y/\delta_e)},\ \ \delta_e=\sqrt{2\nu/\omega}$ · *why:* Split
     e^{−(1+i)y/δ_e} into a real decay and a phase. · *plain:* amplitude and phase both change with y.
  10. *did:* Take the real part · *tex:* $u(y,t)=Ue^{-y\sqrt{\omega/2\nu}}\cos\big(\omega t-y\sqrt{\omega/2\nu}\big)$ (8.38) · *why:* Take the real part as promised in step 1: Re{e^(iθ)} = cos θ with θ = ωt − y/δ_e. · *plain:* (8.38).
- **Result.** (8.38) — *in words:* each layer repeats the wall's motion, shrunk by e^{−y/δ_e} and delayed by y/δ_e radians.
- **Check.** u(0, t) = U cos ωt ✓; bounded ✓; residual of (8.20) zero (`stokes_second_sympy`); the CN twin after 10 periods
  agrees to 2 × 10⁻³U.
- **What it means.** The same (1 + i)/δ structure gives the Ekman spiral (Ch. 13) and tidal bottom layers.
- **Traps.** ∂/∂t of e^{iωt}f is iωe^{iωt}f; keep only the decaying root; the phase lag carries a minus sign.

### D25 · Reading (8.38): 0.059U at 4√(ν/ω), crest speed √(2νω), e-folding depth √(2ν/ω) — ★, 6 steps, in C11 (notebook · `oscillating_plate`)
- **Goal.** Turn (8.38) into numbers a reader can picture: how deep, how late, how fast the apparent wave moves.
- **Start.** $u=Ue^{-y/\delta_e}\cos(\omega t-y/\delta_e)$, δ_e = √(2ν/ω) (8.38).
- **Plan.** (1) Envelope. (2) The book's depth. (3) Phase lag and crest speed.
- **Tools.** phase of a wave (P165) · exponent rules (P43).
- **Assumptions.** None.
- **Steps.**
  1. *did:* Read the envelope · *tex:* $\lvert u\rvert\le Ue^{-y/\delta_e}$ · *why:* The cosine is at most 1; δ_e is the e-folding
     depth. · *plain:* the motion shrinks by e every δ_e.
  2. *did:* Evaluate at the book's depth · *tex:* $y=4\sqrt{\nu/\omega}:\ \ \frac{y}{\delta_e}=\frac{4}{\sqrt2}=2\sqrt2,\ \ e^{-2\sqrt2}=
     0.0591$ · *why:* 4√(ν/ω)/√(2ν/ω) = 4/√2. The book's δ ~ 4√(ν/ω) is 2√2 = 2.83 δ_e. · *plain:* 6 % of the wall's motion
     left.
  3. *did:* Read the phase lag · *tex:* $\phi(y)=\frac{y}{\delta_e}\ \text{rad},\ \ \Delta t=\frac{y}{\omega\delta_e}$ · *why:* The cosine's
     argument is ωt − y/δ_e: a layer at y reaches the same phase y/δ_e radians later. · *plain:* deeper layers lag.
  4. *did:* Follow a crest · *tex:* $\omega t-\frac{y}{\delta_e}=\text{const}\ \Rightarrow\ \frac{dy}{dt}=\omega\delta_e=\sqrt{2\nu\omega}$ ·
     *why:* A point of fixed phase moves; differentiate in t. · *plain:* crests climb at √(2νω).
  5. *did:* Read the 'wavelength' · *tex:* $\lambda=2\pi\delta_e$ · *why:* The phase changes by 2π over Δy = 2πδ_e — but the amplitude
     has fallen by e^{−2π} ≈ 0.2 % by then. · *plain:* barely one 'wave' is ever visible.
  6. *did:* Put in numbers · *tex:* $\omega=2\pi,\ \nu=10^{-6}:\ \ \delta_e=0.56\ \mathrm{mm},\ \ 4\sqrt{\nu/\omega}=1.6\ \mathrm{mm},\ \ \sqrt{2\nu
     \omega}=3.5\ \mathrm{mm/s}$ · *why:* Direct evaluation with ω = 2π rad/s (1 Hz) and water's ν = 10⁻⁶ m²/s. · *plain:* a plate shaken once a second in water moves only a
     millimetre of fluid.
- **Result.** δ_e = √(2ν/ω), amplitude 0.0591U at 4√(ν/ω), lag y/δ_e, crest speed √(2νω) — *in words:* a diffusion layer
  with a phase delay.
- **Check.** Units: √(m²/s ÷ 1/s) = m ✓; √(m²/s × 1/s) = m/s ✓.
- **What it means.** It looks like a wave, but there is no restoring force: the pattern is diffusion from a periodic
  boundary.
- **Traps.** Book depth (4√(ν/ω)) and e-folding depth (√(2ν/ω)) differ by 2√2; calling it a real wave.

### D26 · From (8.39) to the Stokes equations $\nabla p=\mu\nabla^2\mathbf u$ (8.43) — ★★, 9 steps, in C12 (notebook)
- **Goal.** Show how to drop inertia correctly when Re → 0, and why pressure must be rescaled.
- **Start.** $\rho\mathbf u\cdot\nabla\mathbf u+\nabla p=\mu\nabla^2\mathbf u$ (8.39).
- **Plan.** (1) Scale with the high-Re pressure (4.100). (2) Multiply by Re and see the failure. (3) Rescale pressure by
  dominant balance. (4) Take the limit.
- **Tools.** scaled variables (P133) · dominant balance (primer P198).
- **Assumptions.** Steady, no body force; Re → 0 at the end.
- **Steps.**
  1. *did:* Scale with the dynamic pressure · *tex:* $\mathbf u=U\mathbf u^*,\ \mathbf x=L\mathbf x^*,\ p-p_\infty=\rho U^2p^*$ · *why:*
     The scaling (4.100), natural when inertia balances pressure. · *plain:* the high-Re yardsticks.
  2. *did:* Size each term · *tex:* $\frac{\rho U^2}{L}\mathbf u^*\cdot\nabla^*\mathbf u^*+\frac{\rho U^2}{L}\nabla^*p^*=\frac{\mu U}{L^2}\nabla^{*2}
     \mathbf u^*$ · *why:* ∇ = ∇*/L; two derivatives on the viscous term. · *plain:* inertia and pressure ~ ρU²/L, friction ~
     μU/L².
  3. *did:* Divide by ρU²/L · *tex:* $\mathbf u^*\cdot\nabla^*\mathbf u^*+\nabla^*p^*=\frac{1}{\mathrm{Re}}\nabla^{*2}\mathbf u^*$ (8.40) · *why:* Divide every term by ρU²/L; the viscous coefficient (μU/L²)/(ρU²/L) = μ/(ρUL) = 1/Re. · *plain:* (8.40).
  4. *did:* Multiply by Re · *tex:* $\mathrm{Re}(\mathbf u^*\cdot\nabla^*\mathbf u^*+\nabla^*p^*)=\nabla^{*2}\mathbf u^*$ (8.41) · *why:* Look
     at the low-Re end: the viscous term now has coefficient 1. · *plain:* (8.41).
  5. *did:* Let Re → 0 naively · *tex:* $0=\nabla^{*2}\mathbf u^*\ \Leftrightarrow\ 0=\mu\nabla^2\mathbf u$ · *why:* The pressure term
     carries Re too and disappears — leaving an equation with no force balance. · *plain:* wrong: pressure cannot vanish.
  6. *did:* Rescale pressure by the viscous stress · *tex:* $p-p_\infty=\frac{\mu U}{L}p^*$ · *why:* Dominant balance (P198): at low
     Re pressure differences push against viscous stresses ~ μU/L, not ρU². · *plain:* the right yardstick.
  7. *did:* Divide (8.39) by μU/L² · *tex:* $\mathrm{Re}\,\mathbf u^*\cdot\nabla^*\mathbf u^*+\nabla^*p^*=\nabla^{*2}\mathbf u^*$ · *why:*
     Inertia (ρU²/L)/(μU/L²) = Re; pressure (μU/L²)/(μU/L²) = 1. · *plain:* only inertia carries Re now.
  8. *did:* Rearrange · *tex:* $\mathrm{Re}(\mathbf u^*\cdot\nabla^*\mathbf u^*)=-\nabla^*p^*+\nabla^{*2}\mathbf u^*$ (8.42) · *why:* Move the
     pressure term to the right. · *plain:* (8.42).
  9. *did:* Let Re → 0 and restore dimensions · *tex:* $\nabla p=\mu\nabla^2\mathbf u$ (8.43) · *why:* The left side vanishes; multiply
     by μU/L² and undo the scalings. · *plain:* (8.43): pressure gradient balances friction everywhere.
- **Result.** (8.43), with (4.10) — *in words:* creeping flow is linear: pressure and viscosity balance.
- **Check.** `low_re_scaling_sympy`: dynamic {1, 1, 1/Re}, ×Re {Re, Re, 1}, viscous {Re, 1, 1}. Droplet: μU/a = 0.022 Pa vs
  ρU² = 1.7 × 10⁻⁴ Pa.
- **What it means.** Linearity: solutions add and reversing U reverses the flow (no wake). The italic rule — scales come
  from the terms that balance — is how Ch. 9 and Ch. 13 build their approximations.
- **Traps.** Multiplying (8.40) by Re and letting Re → 0 kills the pressure too; re-dimensionalise at the end.

### D27 · Curl of the Stokes equations: ∇²ω = 0 — ★★, 6 steps, in C13 (notebook)
- **Goal.** Remove the pressure from (8.43) and get an equation for the vorticity alone.
- **Start.** $\nabla p=\mu\nabla^2\mathbf u$ (8.43), in Cartesian components ∂_i p = μ∂_j∂_j u_i.
- **Plan.** (1) Take the curl. (2) The pressure drops out. (3) Swap the curl and the Laplacian. (4) Interpret in spherical
  coordinates.
- **Tools.** index notation ε_ijk (Ch. 2 §2.7) · Schwarz's theorem (P121) · curl of a curl (P122).
- **Assumptions.** Constant μ; Cartesian components for steps 1–5.
- **Steps.**
  1. *did:* Apply the curl · *tex:* $\varepsilon_{kli}\partial_l\partial_ip=\mu\,\varepsilon_{kli}\partial_l\partial_j\partial_ju_i$ · *why:* Multiply
     both sides by ε_kli∂_l and sum; this is the curl in index notation. We do it to eliminate p. · *plain:* take the curl
     of both sides.
  2. *did:* Drop the pressure · *tex:* $\varepsilon_{kli}\partial_l\partial_ip=0$ · *why:* ∂_l∂_i is symmetric in l, i (Schwarz, P121)
     while ε_kli is antisymmetric: the sum vanishes — the curl of a gradient is zero. · *plain:* pressure has no curl.
  3. *did:* Move the Laplacian outside · *tex:* $0=\mu\,\partial_j\partial_j(\varepsilon_{kli}\partial_lu_i)$ · *why:* In Cartesian
     components ε is constant and partial derivatives commute. · *plain:* curl and Laplacian can swap order.
  4. *did:* Recognise the vorticity · *tex:* $\omega_k=\varepsilon_{kli}\partial_lu_i\ \Rightarrow\ 0=\mu\nabla^2\omega_k$ · *why:* ω = ∇ × u
     in index form. · *plain:* each Cartesian component of ω is harmonic.
  5. *did:* Write it as a vector statement · *tex:* $\nabla^2\boldsymbol\omega=0$ · *why:* Divide by the constant viscosity μ, which is not zero. · *plain:* vorticity
     neither diffuses nor is advected here — it is simply harmonic.
  6. *did:* Say what it means in spherical form · *tex:* $\nabla^2\boldsymbol\omega=\nabla(\nabla\cdot\boldsymbol\omega)-\nabla\times\nabla\times
     \boldsymbol\omega=-\nabla\times\nabla\times\boldsymbol\omega=0$ · *why:* The curl-of-curl identity (P122) with ∇·ω = 0; in
     spherical components this, not the scalar Laplacian of each component, is meant (the book's footnote). · *plain:* the
     form used in D28.
- **Result.** $-\nabla\times\nabla\times\boldsymbol\omega=0$ (∇²ω = 0 in Cartesian components) — *in words:* in creeping flow the
  vorticity satisfies a pressure-free equation.
- **Check.** For Stokes' solution ω_φ = −(3Ua/2r²) sin θ (D30 step 3), `stokes_sphere_sympy()["curlcurl_identity"]` gives 0.
- **What it means.** The pressure is recovered afterwards from ∇p = −μ∇×ω (D30).
- **Traps.** Curl and Laplacian commute only in Cartesian components; applying the scalar Laplacian to ω_φ alone is wrong.

### D28 · The stream-function equation $\big[\frac{\partial^2}{\partial r^2}+\frac{\sin\theta}{r^2}\frac{\partial}{\partial\theta}\big(\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\big)\big]^2\psi=0$ (8.44) — ★★★, 12 steps, in C13 (notebook · `stokes_sphere_flow`)
- **Goal.** Write the vorticity equation of D27 as one scalar equation for Stokes' stream function ψ.
- **Start.** $u_r=\frac{1}{r^2\sin\theta}\frac{\partial\psi}{\partial\theta}$, $u_\theta=-\frac{1}{r\sin\theta}\frac{\partial\psi}{\partial r}$ (6.83);
  ω = ω_φ e_φ; $-\nabla\times\nabla\times\boldsymbol\omega=0$ (D27).
- **Plan.** (1) Express ω_φ through ψ with the operator E². (2) Show the curl–curl of any azimuthal field is an E²
  operation. (3) Combine.
- **Tools.** Stokes stream function (R18) · spherical vorticity (R17) · the Stokes operator E² (primer P199) · sympy.
- **Assumptions.** Axisymmetric (∂/∂φ = 0), no swirl (u_φ = 0); θ from the downstream axis.
- **Steps.**
  1. *did:* Keep only the azimuthal vorticity · *tex:* $\boldsymbol\omega=\omega_\varphi(r,\theta)\,\mathbf e_\varphi$ · *why:* With
     u_φ = 0 and ∂/∂φ = 0 the r- and θ-components of the spherical curl vanish. · *plain:* vorticity circles the axis.
  2. *did:* Write ω_φ from the velocity · *tex:* $\omega_\varphi=\frac1r\Big[\frac{\partial(ru_\theta)}{\partial r}-\frac{\partial u_r}{\partial
     \theta}\Big]$ · *why:* The φ-component of the spherical curl (R17). · *plain:* a rotation rate in the meridian plane.
  3. *did:* Insert u_θ from (6.83) · *tex:* $\frac{\partial(ru_\theta)}{\partial r}=-\frac{1}{\sin\theta}\frac{\partial^2\psi}{\partial r^2}$ ·
     *why:* ru_θ = −(1/sin θ)∂ψ/∂r, and sin θ does not depend on r. · *plain:* first half.
  4. *did:* Insert u_r from (6.83) · *tex:* $\frac{\partial u_r}{\partial\theta}=\frac{1}{r^2}\frac{\partial}{\partial\theta}\Big(\frac{1}{\sin\theta}
     \frac{\partial\psi}{\partial\theta}\Big)$ · *why:* u_r = (1/r²)(1/sin θ)∂ψ/∂θ; differentiate in θ at fixed r. · *plain:*
     second half.
  5. *did:* Combine · *tex:* $\omega_\varphi=-\frac1r\Big[\frac{1}{\sin\theta}\frac{\partial^2\psi}{\partial r^2}+\frac1{r^2}\frac{\partial}
     {\partial\theta}\Big(\frac1{\sin\theta}\frac{\partial\psi}{\partial\theta}\Big)\Big]$ · *why:* Insert steps 3 and 4 into the curl formula of step 2 and collect the bracket. · *plain:* the
     book's line.
  6. *did:* Factor out 1/sin θ · *tex:* $\omega_\varphi=-\frac{E^2\psi}{r\sin\theta},\ \ E^2=\frac{\partial^2}{\partial r^2}+\frac{\sin\theta}{r^2}
     \frac{\partial}{\partial\theta}\Big(\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\Big)$ · *why:* Multiply and divide the bracket by
     sin θ; the operator is Ch. 6's E² (P199). · *plain:* vorticity = −E²ψ/(r sin θ).
  7. *did:* Curl an azimuthal field · *tex:* $\nabla\times(A\,\mathbf e_\varphi)=\frac{\partial_\theta(\sin\theta A)}{r\sin\theta}\mathbf e_r-\frac
     {\partial_r(rA)}{r}\mathbf e_\theta$ · *why:* Spherical curl of A = A(r, θ)e_φ (`core.curvilinear`). We need the curl
     twice. · *plain:* an azimuthal field curls into a meridional one.
  8. *did:* Introduce Ψ = r sin θ A · *tex:* $\nabla\times(A\,\mathbf e_\varphi)=\frac{1}{r^2\sin\theta}\frac{\partial\Psi}{\partial\theta}
     \mathbf e_r-\frac{1}{r\sin\theta}\frac{\partial\Psi}{\partial r}\mathbf e_\theta$ · *why:* sin θ A = Ψ/r and rA = Ψ/sin θ
     turn step 7 into exactly the (6.83) pattern with Ψ in place of ψ. · *plain:* the curl looks like a velocity from a
     stream function Ψ.
  9. *did:* Curl again using step 6 · *tex:* $\nabla\times\nabla\times(A\,\mathbf e_\varphi)=-\frac{E^2\Psi}{r\sin\theta}\mathbf e_\varphi$ · *why:*
     Steps 2–6 showed: the curl of a (6.83)-type field from Ψ is −E²Ψ/(r sin θ) e_φ. · *plain:* curl–curl is an E²
     operation.
  10. *did:* Apply to the vorticity · *tex:* $\Psi=r\sin\theta\,\omega_\varphi=-E^2\psi$ · *why:* A = ω_φ, and step 6. · *plain:* the
      'stream function' of the vorticity is −E²ψ.
  11. *did:* Insert into D27's equation · *tex:* $-\nabla\times\nabla\times\boldsymbol\omega=\frac{E^2\Psi}{r\sin\theta}\mathbf e_\varphi=-\frac
      {E^2(E^2\psi)}{r\sin\theta}\mathbf e_\varphi=0$ · *why:* Insert step 10 into step 9, with A = ω_φ, into D27's equation. · *plain:* one scalar equation left.
  12. *did:* Divide by r sin θ · *tex:* $E^2(E^2\psi)=0$ (8.44) · *why:* r sin θ ≠ 0 off the axis. It is the square of E², not
      the biharmonic ∇⁴ (the book's footnote). · *plain:* (8.44).
- **Result.** $(E^2)^2\psi=0$ (8.44) — *in words:* creeping axisymmetric flow obeys the squared Stokes operator.
- **Check.** Units: E² is 1/m², so E⁴ψ ~ ψ/m⁴ ✓; sympy: E⁴ψ = 0 for (8.48), ∇⁴ψ ≠ 0 for (8.48), and the identity of step 9
  holds for a generic A(r, θ) (check_src in C13 item 28).
- **What it means.** A fourth-order equation needs four conditions: two on the sphere, the stream far away (and regularity)
  — D29 uses them.
- **Traps.** (E²)² is not ∇⁴; θ from the downstream axis; the identity of steps 7–9 must be checked, not assumed.
- **sympy check intent (★★★).** Build the spherical curl of A e_φ with `core.curvilinear.curl`, apply it twice, compare with
  −E²(r sin θ A)/(r sin θ) for a generic function A (difference 0); then evaluate E²E²ψ for (8.48) (0) and the scalar
  biharmonic of ψ (non-zero).

### D29 · Separable solution, f = Ar⁴ + Br² + Cr + D/r and Stokes' ψ (8.48) — ★★, 10 steps, in C13 (notebook · `stokes_sphere_flow`)
- **Goal.** Solve (8.44) with the sphere's conditions.
- **Start.** $(E^2)^2\psi=0$ (8.44), $\psi(a,\theta)=0$ (8.45), $\partial\psi(a,\theta)/\partial r=0$ (8.46), $\psi\to\frac12Ur^2\sin^2\theta$
  (8.47).
- **Plan.** (1) Separate ψ = f(r) sin²θ. (2) Get the ODE for f. (3) Euler–Cauchy roots. (4) Four constants from the
  conditions.
- **Tools.** E² (primer P199) · Euler–Cauchy ODE (primer P187) · 2 × 2 system (P57).
- **Assumptions.** The angular dependence of (8.47) persists everywhere (step 1).
- **Steps.**
  1. *did:* Try a separable form · *tex:* $\psi=f(r)\sin^2\theta$ · *why:* The far field (8.47) has sin²θ; if E² keeps sin²θ, the
     whole problem reduces to r. · *plain:* guess the angular shape from infinity.
  2. *did:* Apply the angular part · *tex:* $\sin\theta\frac{\partial}{\partial\theta}\Big(\frac{1}{\sin\theta}\frac{\partial\sin^2\theta}{\partial\theta}
     \Big)=\sin\theta\frac{\partial}{\partial\theta}(2\cos\theta)=-2\sin^2\theta$ · *why:* ∂_θ sin²θ = 2 sin θ cos θ; divided by
     sin θ it is 2 cos θ. · *plain:* sin²θ is reproduced.
  3. *did:* E² of the trial · *tex:* $E^2(f\sin^2\theta)=\Big(f''-\frac{2f}{r^2}\Big)\sin^2\theta$ · *why:* The radial part gives f″ sin²θ;
     the angular part gives (1/r²)(−2 sin²θ) f. · *plain:* E² acts on f alone.
  4. *did:* Apply E² again · *tex:* $f^{iv}-\frac{4f''}{r^2}+\frac{8f'}{r^3}-\frac{8f}{r^4}=0$ · *why:* With g = f″ − 2f/r², (8.44) is
     g″ − 2g/r² = 0; expanding (f/r²)″ = f″/r² − 4f′/r³ + 6f/r⁴ gives the book's ODE. · *plain:* a fourth-order Euler–Cauchy
     ODE.
  5. *did:* Try f = r^λ · *tex:* $\lambda(\lambda-1)(\lambda-2)(\lambda-3)-4\lambda(\lambda-1)+8\lambda-8=0$ · *why:* Every term is r^{λ−4}
     times a polynomial in λ (P187). · *plain:* an algebraic equation.
  6. *did:* Factor the polynomial · *tex:* $(\lambda-4)(\lambda-2)(\lambda-1)(\lambda+1)=0$ · *why:* 8λ − 8 = 8(λ − 1) and the rest also
     carries (λ − 1); the cubic λ³ − 5λ² + 2λ + 8 has the roots 4, 2, −1. · *plain:* four powers.
  7. *did:* Write the general solution · *tex:* $f=Ar^4+Br^2+Cr+\frac Dr$ · *why:* The ODE is linear, so any combination of the four power solutions solves it. · *plain:* four
     constants.
  8. *did:* Apply the far field (8.47) · *tex:* $A=0,\ \ B=\frac U2$ · *why:* f must tend to ½Ur²; an r⁴ term would outgrow the
     stream, and the r and 1/r terms are negligible far away. · *plain:* the uniform stream survives.
  9. *did:* Apply (8.45) and (8.46) · *tex:* $\frac{Ua^2}{2}+Ca+\frac Da=0,\ \ Ua+C-\frac{D}{a^2}=0$ · *why:* ψ = 0 and ∂ψ/∂r = 0 at
     r = a mean f(a) = f′(a) = 0. · *plain:* no flow through, no slip.
  10. *did:* Solve and assemble · *tex:* $C=-\frac{3Ua}{4},\ D=\frac{Ua^3}{4}:\ \ \psi=Ur^2\sin^2\theta\Big(\frac12-\frac{3a}{4r}+\frac{a^3}{4r^3}
      \Big)$ (8.48) · *why:* Divide the first by a and add the second: 3Ua/2 + 2C = 0; then D = a²(Ua + C). · *plain:* (8.48).
- **Result.** (8.48) — *in words:* stream, minus a 1/r 'Stokeslet', plus a 1/r³ doublet.
- **Check.** ψ(a) = Ua²(½ − ¾ + ¼) = 0 ✓; ∂ψ/∂r at a: Ua sin²θ(1 − ¾ − ¼) = 0 ✓; E⁴ψ = 0 (sympy).
- **What it means.** The −3a/4r term decays slowly: the sphere is felt tens of radii away, which is exactly why inertia
  returns far away (D32).
- **Traps.** E²(f sin²θ) = (f″ − 2f/r²) sin²θ; the roots are {4, 2, 1, −1}; A = 0 from the far field, not from the
  surface.

### D30 · The pressure $p-p_\infty=-\frac{3\mu aU\cos\theta}{2r^2}$ (8.50) by integrating ∇p = μ∇²u — ★★★, 11 steps, in C14 (notebook · `stokes_drag_settling`)
- **Goal.** Find the pressure round the sphere (the book gives only the result).
- **Start.** $\nabla p=\mu\nabla^2\mathbf u$ (8.43), ψ from (8.48), ω_φ = −E²ψ/(r sin θ) (D28).
- **Plan.** (1) Replace ∇²u by −∇×ω. (2) Compute ω and its curl. (3) Integrate ∂p/∂r; check ∂p/∂θ. (4) Evaluate extremes.
- **Tools.** curl of a curl (P122) · the azimuthal curl of D28 step 8 · line integral of a gradient (gloss).
- **Assumptions.** p → p∞ far away (step 9).
- **Steps.**
  1. *did:* Use the vector identity · *tex:* $\nabla^2\mathbf u=\nabla(\nabla\cdot\mathbf u)-\nabla\times\boldsymbol\omega=-\nabla\times\boldsymbol\omega\
     \Rightarrow\ \nabla p=-\mu\nabla\times\boldsymbol\omega$ · *why:* Curl-of-curl (P122) with ∇·u = 0. The vorticity is easier
     to curl than the vector Laplacian in spherical coordinates. · *plain:* pressure gradient = −μ curl of vorticity.
  2. *did:* Compute E²ψ · *tex:* $f''-\frac{2f}{r^2}=\Big(U+\frac{Ua^3}{2r^3}\Big)-\Big(U-\frac{3Ua}{2r}+\frac{Ua^3}{2r^3}\Big)=\frac{3Ua}{2r}$ ·
     *why:* With f = ½Ur² − (3Ua/4)r + (Ua³/4)/r and D29 step 3. · *plain:* E²ψ = (3Ua/2r) sin²θ.
  3. *did:* Get the vorticity · *tex:* $\omega_\varphi=-\frac{E^2\psi}{r\sin\theta}=-\frac{3Ua}{2r^2}\sin\theta$ · *why:* Use ω_φ = −E²ψ/(r sin θ) from D28 step 6 with step 2. ·
     *plain:* vorticity decays as 1/r².
  4. *did:* Form Ψ = r sin θ ω_φ · *tex:* $\Psi=-\frac{3Ua}{2r}\sin^2\theta$ · *why:* The curl of an azimuthal field is read off Ψ
     (D28 step 8). · *plain:* ready to curl.
  5. *did:* Radial component of ∇×ω · *tex:* $(\nabla\times\boldsymbol\omega)_r=\frac{1}{r^2\sin\theta}\frac{\partial\Psi}{\partial\theta}=-\frac{3Ua
     \cos\theta}{r^3}$ · *why:* ∂_θ sin²θ = 2 sin θ cos θ. · *plain:* radial part.
  6. *did:* Polar component of ∇×ω · *tex:* $(\nabla\times\boldsymbol\omega)_\theta=-\frac{1}{r\sin\theta}\frac{\partial\Psi}{\partial r}=-\frac{3Ua
     \sin\theta}{2r^3}$ · *why:* Use the D28 step 8 formula with ∂_r(−3Ua/2r) = +3Ua/(2r²). · *plain:* polar part.
  7. *did:* Write both pressure derivatives · *tex:* $\frac{\partial p}{\partial r}=\frac{3\mu Ua\cos\theta}{r^3},\ \ \frac1r\frac{\partial p}
     {\partial\theta}=\frac{3\mu Ua\sin\theta}{2r^3}$ · *why:* ∇p = −μ∇×ω, component by component. · *plain:* two equations for
     one function.
  8. *did:* Integrate the radial one · *tex:* $p=-\frac{3\mu Ua\cos\theta}{2r^2}+g(\theta)$ · *why:* Antiderivative of r⁻³ at fixed θ;
     the constant may depend on θ (gloss: line integral of a gradient). · *plain:* pressure up to a function of θ.
  9. *did:* Check the polar one · *tex:* $\frac1r\frac{\partial p}{\partial\theta}=\frac{3\mu Ua\sin\theta}{2r^3}+\frac{g'(\theta)}{r}\ \Rightarrow\
     g'=0$ · *why:* Compare with step 7: both components agree only if g is constant — a real check that ∇p is a gradient. ·
     *plain:* consistent.
  10. *did:* Fix the constant far away · *tex:* $p-p_\infty=-\frac{3\mu aU\cos\theta}{2r^2}$ (8.50) · *why:* The first term → 0 as
      r → ∞, so g = p∞. · *plain:* (8.50).
  11. *did:* Evaluate front and rear · *tex:* $\theta=\pi:\ +\frac{3\mu U}{2a},\ \ \theta=0:\ -\frac{3\mu U}{2a}$ · *why:* cos π = −1 at the
      front stagnation point, cos 0 = 1 at the rear. The book prints the minimum without its minus sign. · *plain:* high in
      front, low behind.
- **Result.** (8.50) — *in words:* the pressure disturbance is a dipole decaying as 1/r², +3μU/2a in front and −3μU/2a
  behind.
- **Check.** Units μaU/r²: (Pa s)(m)(m/s)/m² = Pa ✓; sympy: both components of ∇p = μ∇²u hold (check_src in C14 item 45).
  Droplet: 3μU/2a = 0.033 Pa.
- **What it means.** Pressure alone pushes the sphere downstream: its integral gives one third of the drag (D31).
- **Traps.** Integrate one component and *check* the other; p → p∞ fixes the constant; the rear value is negative.
- **sympy check intent (★★★).** From ψ (8.48): u_r, u_θ by (6.83), the spherical vector Laplacian of u from
  `core.curvilinear.vector_laplacian`, and ∇p of (8.50): both components of ∇p − μ∇²u simplify to 0.

### D31 · Stokes drag $D=6\pi\mu aU$ (8.51): 2πμaU from pressure + 4πμaU from friction — ★★★, 13 steps, in C14 (notebook · `stokes_drag_settling`)
- **Goal.** Integrate the stresses over the sphere and show where the drag comes from (Exercise 8.35).
- **Start.** (8.49) $u_r=U\cos\theta(1-\frac{3a}{2r}+\frac{a^3}{2r^3})$, $u_\theta=-U\sin\theta(1-\frac{3a}{4r}-\frac{a^3}{4r^3})$, (8.50), and
  the Newtonian stress σ = −pI + 2μS.
- **Plan.** (1) The traction and its x-projection. (2) σ_rr and σ_rθ on r = a. (3) Integrate pressure and friction parts.
- **Tools.** Newtonian stress in spherical coordinates (Ch. 4 §4.5) · traction projection (gloss) · surface integrals on a
  sphere (P163).
- **Assumptions.** Rigid sphere with no slip (u = 0 on r = a).
- **Steps.**
  1. *did:* Write the force on the sphere · *tex:* $\mathbf F=\oint\boldsymbol\sigma\cdot\mathbf e_r\,dA,\ \ D=F_x$ · *why:* Cauchy (Ch. 2
     §2.6): the fluid outside exerts the traction σ·n, with n = e_r pointing into the fluid. · *plain:* drag = x-force.
  2. *did:* Project the traction on x · *tex:* $t_x=\sigma_{rr}\cos\theta-\sigma_{r\theta}\sin\theta$ · *why:* e_r·e_x = cos θ and
     e_θ·e_x = −sin θ with θ measured from +x. · *plain:* normal and shear parts both push in x.
  3. *did:* Write σ_rr · *tex:* $\sigma_{rr}=-p+2\mu\frac{\partial u_r}{\partial r}$ · *why:* Normal stress of a Newtonian fluid in
     spherical coordinates. · *plain:* pressure plus viscous normal stress.
  4. *did:* Evaluate ∂u_r/∂r on the sphere · *tex:* $\frac{\partial u_r}{\partial r}\Big\rvert_a=U\cos\theta\Big(\frac{3}{2a}-\frac{3}{2a}\Big)=0$ ·
     *why:* ∂_r(1 − 3a/2r + a³/2r³) = 3a/2r² − 3a³/2r⁴, which vanishes at r = a. · *plain:* on the sphere σ_rr = −p.
  5. *did:* Write σ_rθ · *tex:* $\sigma_{r\theta}=\mu\Big[r\frac{\partial}{\partial r}\Big(\frac{u_\theta}{r}\Big)+\frac1r\frac{\partial u_r}{\partial
     \theta}\Big]$ · *why:* Shear stress in spherical coordinates (Ch. 4 §4.5). · *plain:* the friction.
  6. *did:* Drop the u_r term on the sphere · *tex:* $\frac1r\frac{\partial u_r}{\partial\theta}\Big\rvert_a=0$ · *why:* u_r = 0 at every θ on
     r = a, so its θ-derivative is zero. · *plain:* only the radial shear remains.
  7. *did:* Evaluate the radial shear · *tex:* $\sigma_{r\theta}\rvert_a=-\frac{3\mu U}{2a}\sin\theta$ · *why:* u_θ/r = −U sin θ(1/r −
     3a/4r² − a³/4r⁴); its r-derivative at a is −U sin θ(3/2)/a², times μr = μa. · *plain:* friction ∝ sin θ, largest at the
     sides.
  8. *did:* Insert the surface pressure · *tex:* $\sigma_{rr}\rvert_a=-p_\infty+\frac{3\mu U}{2a}\cos\theta$ · *why:* Evaluate the pressure (8.50) on the surface r = a, where σ_rr = −p (step 4). ·
     *plain:* normal stress on the surface.
  9. *did:* Assemble t_x · *tex:* $t_x=-p_\infty\cos\theta+\frac{3\mu U}{2a}\cos^2\theta+\frac{3\mu U}{2a}\sin^2\theta$ · *why:* Steps 7–8
     into step 2 (−σ_rθ sin θ = +(3μU/2a) sin²θ). · *plain:* pressure part ∝ cos²θ, friction part ∝ sin²θ.
  10. *did:* Write the area element · *tex:* $dA=2\pi a^2\sin\theta\,d\theta$ · *why:* A ring on the sphere at angle θ (P163). ·
      *plain:* integrate over θ from 0 to π.
  11. *did:* Integrate the ambient pressure · *tex:* $-p_\infty\int_0^\pi\cos\theta\,2\pi a^2\sin\theta\,d\theta=0$ · *why:* cos θ sin θ is
      odd about θ = π/2. · *plain:* a uniform pressure exerts no net force.
  12. *did:* Integrate the pressure part · *tex:* $\frac{3\mu U}{2a}2\pi a^2\int_0^\pi\cos^2\theta\sin\theta\,d\theta=3\pi\mu aU\cdot\frac23=2\pi
      \mu aU$ · *why:* ∫₀^π cos²θ sin θ dθ = 2/3 (substitute c = cos θ). · *plain:* one third of the drag.
  13. *did:* Integrate the friction part · *tex:* $3\pi\mu aU\int_0^\pi\sin^3\theta\,d\theta=3\pi\mu aU\cdot\frac43=4\pi\mu aU,\ \ D=6\pi\mu aU$
      (8.51) · *why:* ∫₀^π sin³θ dθ = 4/3; add step 12. · *plain:* two thirds friction; (8.51).
- **Result.** $D=6\pi\mu aU$ (8.51), ⅓ pressure + ⅔ friction — *in words:* drag is proportional to speed (Stokes' law of
  resistance).
- **Check.** Units (Pa s)(m)(m/s) = N ✓. Notice t_x = 3μU/2a − p∞ cos θ: the x-traction is uniform over the sphere. The
  running integrals πμaU(1 − cos³θ) and πμaU(2 − 3cos θ + cos³θ) (`stokes_drag_running`) reach 2π and 4π at θ = π; midpoint
  and Gauss–Legendre sums agree; C_D = 24/Re (8.52) follows.
- **What it means.** Drag ∝ U makes the fall speed ∝ a² (N90): a 10 µm droplet falls 1.2 cm/s.
- **Traps.** The viscous normal stress vanishes on the sphere; project with θ from the downstream axis (signs!); a bubble
  with a slipping surface gives 4πμaU — not this case.
- **sympy check intent (★★★).** From (8.49), (8.50) compute σ_rr and σ_rθ on r = a with `core.curvilinear.strain_rate`
  (spherical), integrate each x-projection over the sphere, and assert 2πμaU, 4πμaU and their sum 6πμaU (check_src in C14
  item 47).

### D32 · Far-field breakdown: inertia/viscous ~ $\frac{\rho Ua}{\mu}\frac ra=\mathrm{Re}_a\frac ra$ — ★★, 7 steps, in C15 (notebook · `stokes_sphere_flow`)
- **Goal.** Show that however small Re is, inertia matches viscosity at some distance from the sphere.
- **Start.** (8.49) for r ≫ a.
- **Plan.** (1) Size the disturbance far away. (2) Size inertia and viscous terms. (3) Take the ratio.
- **Tools.** asymptotic size of a term (gloss, P130).
- **Assumptions.** r ≫ a; orders of magnitude only.
- **Steps.**
  1. *did:* Keep the slowest-decaying disturbance · *tex:* $\lvert\mathbf u-U\mathbf e_x\rvert\sim\frac{Ua}{r}$ · *why:* In (8.49) the 3a/2r
     and 3a/4r terms decay like 1/r; the a³/r³ terms are much smaller far away. · *plain:* the sphere is felt as Ua/r.
  2. *did:* Size the velocity gradient · *tex:* $\frac{\partial u_\theta}{\partial r}\sim\frac{Ua}{r^2}$ · *why:* Differentiating a 1/r
     term brings another 1/r. · *plain:* gradients fall as 1/r².
  3. *did:* Size the inertia · *tex:* $\rho u_r\frac{\partial u_\theta}{\partial r}\sim\rho U\cdot\frac{Ua}{r^2}=\frac{\rho U^2a}{r^2}$ · *why:* Far
     away u_r ~ U (the stream carries the disturbance). · *plain:* inertia ~ ρU²a/r².
  4. *did:* Size the viscous force · *tex:* $\mu\nabla^2\mathbf u\sim\mu\frac{Ua/r}{r^2}=\frac{\mu Ua}{r^3}$ · *why:* Two derivatives of the
     Ua/r disturbance (the uniform part has none). · *plain:* friction ~ μUa/r³.
  5. *did:* Take the ratio · *tex:* $\frac{\text{inertia}}{\text{viscous}}\sim\frac{\rho U^2a/r^2}{\mu Ua/r^3}=\frac{\rho Ur}{\mu}$ · *why:*
     Divide step 3 by step 4. · *plain:* the ratio grows with r.
  6. *did:* Write it with a Reynolds number · *tex:* $\frac{\rho Ur}{\mu}=\frac{\rho Ua}{\mu}\frac ra=\mathrm{Re}_a\frac ra$ · *why:* Re_a
     uses the radius (Re = 2aU/ν = 2Re_a). · *plain:* small Re times a large distance.
  7. *did:* Find the crossover · *tex:* $\mathrm{Re}_a\frac ra\sim1\ \Rightarrow\ \frac ra\sim\frac{1}{\mathrm{Re}_a}$ · *why:* Set the ratio to one: that is where inertia stops being negligible. · *plain:* beyond r ~ a/Re_a Stokes' solution is wrong.
- **Result.** inertia/viscous ~ Re_a r/a — *in words:* the non-uniformity of Stokes flow is at infinity.
- **Check.** `inertia_viscous_ratio` has log–log slope 1 ± 0.05 for r/a ∈ [50, 500]; Re_a = 0.01 → r ≈ 100a; the 10 µm
  droplet (Re_a = 0.008) → r ≈ 125a = 1.25 mm.
- **What it means.** For a cylinder the disturbance decays even more slowly and no Stokes solution exists (Stokes'
  paradox); Oseen keeps U∂u′/∂x to fix the far field (N95).
- **Traps.** The disturbance far away is ~ Ua/r (the Stokeslet), not a³/r³; Re here uses the radius.

### D33 · Oseen's (8.53) reduces to Stokes' (8.48) near the sphere — ★★, 7 steps, in C15 (notebook · `stokes_sphere_flow`)
- **Goal.** Check that Oseen's improved solution agrees with Stokes' where Stokes is valid.
- **Start.** $\frac{\psi}{Ua^2}=\big[\frac{r^2}{2a^2}+\frac a{4r}\big]\sin^2\theta-\frac3{\mathrm{Re}}(1+\cos\theta)\big\{1-\exp\big[-\frac{\mathrm{Re}}4
  \frac ra(1-\cos\theta)\big]\big\}$ (8.53), Re = 2aU/ν.
- **Plan.** (1) Name the exponent. (2) Expand the exponential for small arguments. (3) Simplify with sin²θ = (1 − cos θ)(1 +
  cos θ).
- **Tools.** series of 1 − e^{−s} (gloss, P117).
- **Assumptions.** Re r/a ≪ 1 (near the sphere, small Re) — step 2.
- **Steps.**
  1. *did:* Name the exponent · *tex:* $s=\frac{\mathrm{Re}}{4}\frac ra(1-\cos\theta)$ · *why:* The whole Oseen correction depends on
     this one combination. · *plain:* s measures inertia's reach.
  2. *did:* State where s is small · *tex:* $s\ll1\iff\mathrm{Re}\,\frac ra\ll1$ · *why:* 0 ≤ 1 − cos θ ≤ 2, so s is small exactly
     where D32 says Stokes holds. · *plain:* near the sphere at small Re.
  3. *did:* Expand the exponential · *tex:* $1-e^{-s}=s-\frac{s^2}{2}+\ldots\approx s$ · *why:* Taylor series of e^{−s} (P117); keep
     the leading term. · *plain:* the bracket is linear in s.
  4. *did:* Insert into the Oseen term · *tex:* $\frac3{\mathrm{Re}}(1+\cos\theta)s=\frac3{\mathrm{Re}}(1+\cos\theta)\frac{\mathrm{Re}}4\frac ra
     (1-\cos\theta)$ · *why:* Replace 1 − e^(−s) by its leading term s from step 3. · *plain:* Re cancels.
  5. *did:* Use sin²θ = (1 + cos θ)(1 − cos θ) · *tex:* $=\frac34\frac ra\sin^2\theta$ · *why:* The product of the two brackets is
     1 − cos²θ. · *plain:* the Stokeslet term appears.
  6. *did:* Collect the terms · *tex:* $\frac{\psi}{Ua^2}\approx\Big[\frac{r^2}{2a^2}-\frac{3r}{4a}+\frac{a}{4r}\Big]\sin^2\theta$ · *why:*
     Subtract step 5 from the first bracket of (8.53). · *plain:* three terms, like Stokes.
  7. *did:* Multiply by Ua² · *tex:* $\psi\approx Ur^2\sin^2\theta\Big(\frac12-\frac{3a}{4r}+\frac{a^3}{4r^3}\Big)$ (8.48) · *why:* Ua² ×
     [r²/2a², −3r/4a, a/4r] = Ur²[½, −3a/4r, a³/4r³]. The text calls this "(9.68)"; it is (8.48). · *plain:* Stokes'
     solution recovered, with an O(Re) remainder.
- **Result.** ψ_Oseen = ψ_Stokes + O(Re r/a) — *in words:* Oseen changes only the far field.
- **Check.** `oseen_limit_sympy()["difference"]` = 0 in the limit Re → 0; numerically ψ agree to ~1e-6 at Re = 10⁻⁶.
  Oseen meets no slip on r = a only to O(Re).
- **What it means.** Near the sphere the two agree; far away Oseen's exponential sweeps the disturbance downstream into a
  wake, and the drag gains the factor (1 + 3Re/16).
- **Traps.** s is small only where Re r/a ≪ 1; the product (1 + cos θ)(1 − cos θ) is sin²θ; do not assert exact no slip for
  (8.53).
