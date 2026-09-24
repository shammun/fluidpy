# Chapter 6 — Ideal Flow: lesson design
(from `analysis/ch06_curation.md` (A 15 · B 116 · C 29 = 160 rows; CORE 15 · NOTE 106 · RECAP 32 · SKIP 7; 31 derivations
D01–D31 written out (★ 8 · ★★ 18 · ★★★ 5); E1–E9 + backup B1) and `analysis/ch06.md` (§2b derivations a-D01…a-D48, §4
implementation rows 1–30, §5 core modules, §9 slips and conventions); every equation placed below was re-read on the
rendered pages `chapters/pages/ch06/` (printed = pdf − 27): (6.1) and Fig. 6.1 p225, (6.2)–(6.6) p227, (6.7)–(6.23)
p228–p230, (6.24)–(6.28) p231–p232, the log expansion and (6.29) p233, (6.30)–(6.32), the half-body velocity, a, ψ = m/2,
h and h_max p234, (6.33)–(6.35) p235–p236, (6.36)–(6.38) and the off-body root p237–p238, (6.39) p238, the lift integral
and (6.40) p239, image rules p240, (6.41) and its streamline equation p241, Example 6.1 p242–p243, (6.42)–(6.44) p243,
(6.45)–(6.46) and dw/dz = (Aπ/α)z^{(π−α)/α} p244–p245, (6.47)–(6.52) p245, (6.53)–(6.56) p246, (6.57)–(6.60) p247,
the far field, (6.61), (6.62) p248, (6.63)–(6.64) p249–p250, the grid maps and cot z p250, (6.65)–(6.67) p251,
(6.68)–(6.69) p252, half-point differences and (6.70)–(6.71) p253, (6.72) p254, (6.73) p255, (6.74)–(6.79) p258,
(6.80)–(6.88) p259–p260, (6.89) p260, (6.90)–(6.92) p261, (6.93)–(6.95) p262, the axial system p263, (6.96) and Fig. 6.30
p264, (6.97)–(6.102) p265, (6.103)–(6.108) p266, (6.109) p267; 2026-09-23, lesson-designer. The implementer works in
parallel from analysis §4 + curation §8; **Part C is written first and is the contract both sides keep.**)

**Binding conventions for every builder (analysis §9, curation decisions 4–6 and §8, ch05 conventions carried over).**
1. **Imports and aliases.** `from fluidpy import ch06_ideal_flow as ch06`; `from fluidpy.core import potential as pf,
   conformal as cm, laplace_solvers as ls, panels as pn`. **`ch06` re-exports every public name of `core/potential.py`,
   `core/conformal.py`, `core/laplace_solvers.py` and `core/panels.py`** (the ch05 pattern), so explainer parity rows
   write `ch06.<name>` only. Units SI everywhere: x, y, z, r, R, a, b, h [m]; u, v, U [m/s]; plane ψ, φ, Γ, m, k [m²/s];
   **Stokes ψ and Q [m³/s]**; 3-D dipole d [m⁴/s], 2-D dipole d [m³/s]; p [Pa]; ρ [kg/m³]; D, L [N/m] (per unit depth);
   F_s [N]; M [kg]. Complex API: `z = x + 1j*y`, `dwdz(z) = u - 1j*v`.
2. **Every public function is scalar-callable** (scalar in → float out) and **parity-friendly**: a `py:` expression may
   use only `ch06.…`, `np.pi`, numbers, strings, lists and dicts (no builtins, no lambdas, no `abs`, no `.real`), so the
   "state" helpers of Part C.5 return dicts of floats and accept **element specs as dicts** (`{"kind": "source", "m":
   6.283185307179586, "x": 0.0, "y": 0.0}`) or **preset names** (`"cylinder"`, `"rankine"`).
3. **Γ has two signs in this chapter** (⚠️ callout in the notation cell, again at C07, C10, C11, E2, E5, E6). The
   project and (6.6) $\nabla^2\psi=-\Gamma\delta(x-x')\delta(y-y')$, (6.8) $\psi=-\frac{\Gamma}{2\pi}\ln r$ and (6.47)
   $w=-\frac{i\Gamma}{2\pi}\ln(z-z')$ take **Γ counterclockwise positive** (u_θ = +Γ/2πr). The book's (6.36)–(6.40),
   (6.52), (6.61)–(6.62), (6.68) and Example 6.1 use a **clockwise Γ** (the flow's circulation is −Γ) so that
   $L=\rho U\Gamma$ (6.40) comes out positive. Code: `pf.Vortex(Gamma)` is always counterclockwise; body helpers take
   keyword-only `Gamma_cw=` (the book's Γ) **or** `Gamma_ccw=` (project), with `Gamma_ccw = −Gamma_cw`. Numbers for the
   callout: Γ_cw = 2 m²/s, U = 10 m/s, a = 0.1 m, air → stagnation points **below** the axis at −9.16° and −170.84°,
   L = +24 N/m; the same flow has Γ_ccw = −2 m²/s.
4. **Doublet signs.** 2-D dipole vector $\mathbf d=\sum_i\mathbf x_im_i$ points **from the sink to the source**; (6.29)
   $\phi=-\mathbf d\cdot\mathbf x/2\pi r^2$; the cylinder needs $\mathbf d=-2\pi Ua^2\mathbf e_x$ (upstream); (6.49)'s
   scalar d means a dipole $-d\,\mathbf e_x$ (`pf.Doublet.from_book_scalar(d, z0)`); in complex form
   $w=-(d_x+id_y)/2\pi(z-z_0)$. 3-D: (6.88) is a dipole $-d\,\mathbf e_z$; the sphere has d = 2πa³U; (6.92) uses
   $\mathbf d=-d\,\mathbf e_z$; the moving sphere (6.97) has $\mathbf d(t)=+2\pi a^3\mathbf u_s$. Code stores vectors.
5. **Angles.** Polar θ from +x (downstream) in (6.31)–(6.39) and (6.89)–(6.91): the upstream stagnation point is θ = π.
   The half-body ψ of (6.31) uses **θ ∈ [0, 2π)** (`np.mod(np.arctan2(y, x), 2*np.pi)`), whose cut runs along the
   downstream +x axis **inside the body**, so ψ = m/2 on both halves of the dividing streamline (the principal
   `arctan2` range (−π, π] puts the cut on the *upstream* axis — the dividing streamline itself — and gives ψ = −m/2 on
   the lower half: a wrong-variant test). Fig. 6.10's angle is measured from the upstream stagnation point (= π − θ);
   (6.106)'s θ_s from the sphere's velocity (= π − θ of (6.91)). `_deg` in names at the interface only.
6. **Axes and letters.** §6.8: z is the **horizontal** symmetry axis along the stream (not up); §6.9: w is the
   **z-velocity** (not the complex potential), ξ = x − x_s (in §6.8 ξ is the axial coordinate of the line sink), ζ the
   Zhukhovsky variable. Symbols reused inside this chapter (each gets a one-line reminder where it changes meaning): a
   (radius; half-body stagnation distance m/2πU; image distance; line-sink length; Zhukhovsky circle), α (corner angle
   π/n; conformal angle; line-sink angle), b (Zhukhovsky constant), **B (span, not ch04's Bernoulli function)**, d, k
   (line-sink density), **M (added mass, not Mach)**, m (2-D source strength vs body mass in (6.109)), n (corner exponent
   vs unit normal), q (source density), Q (3-D source strength), R (cylindrical radius vs contour radius), z (complex
   variable vs axial coordinate).
7. **Forces.** D, L are forces **on the body** per unit depth; (6.54)'s F is on the fluid (F = −B(De_x + Le_y));
   (6.98)'s F_s is on the sphere. Contours for (6.56) are **counterclockwise** (asserted by a positive signed area).
   Pressure is measured from its local hydrostatic value in (6.1); gravity reappears only in the bubble/ball examples.
8. **Book slips taught in corrected form** (each gets "the book prints X; the correct form is Y" where it is used, and
   is never asserted in code): (6.61)'s 1/z² coefficient $(Ud/\pi-\Gamma^2/4\pi^2)$ and its extra outer square → the
   coefficient is $-(Ud/\pi+\Gamma^2/4\pi^2)$, and only the 1/z term matters (N62, D18 step 6, E5); (6.104)'s middle
   bracket $-\frac1{a^3}\mathbf u_s$ → $+\frac1{a^3}\mathbf u_s$ (N100, D29 step 8, E9); (6.108)'s stray dφ after the
   φ-integration (N104, D30 step 7); §6.3 "differentiation of (6.8)" for the source velocities → (6.15) (N24); §6.6
   "(see Figure 6.5)" → Fig. 6.3 (N66); §6.4 "(6.5) and (6.12), respectively" for φ and ψ → reversed, and "(6.43)
   ensures the equality" → (6.44) (N43); §6.5 "Section 3" → §6.3 (N53); §6.7 "first-order central differences" are
   first-*derivative* differences, second-order accurate (R19, D22 step 3); Example 6.2's FORTRAN loop over I should run
   over J and the inlet Δψ unit m² → m²/s (N74); (6.82)'s $\frac1r\frac{\partial}{\partial r}(r^2u_r)$ →
   $\frac1{r^2}$ (R28; the "= 0" is unaffected). **Two more found by this design** (flag for the implementer): the
   analysis/curation write "∇φ = e_z × ∇ψ" (D04, S01) — the correct relation is $\nabla\psi=\mathbf e_z\times\nabla\phi$,
   i.e. $\nabla\phi=\nabla\psi\times\mathbf e_z=-\mathbf e_z\times\nabla\psi$ (D04 step 5); and (6.54) calls F "the force
   applied to the fluid" while its A* hugs the body — we obtain (6.55) directly from Newton's third law with the
   body's outward normal (D16 steps 2–3), which is what the book's (6.56) uses.
9. **Colours (text, figures, explainers — one meaning each; curation §5):** uniform stream `blue` · sources `teal` ·
   sinks `rose` · vortices `amber` · doublets `accent` purple · body / dividing streamline black · pressure above p∞
   `orange`, suction `blue` · steady pressure part `teal`, acceleration part `orange` · lift `amber`, drag `muted` ·
   ψ contours solid, φ contours dashed · residual / error black · ghosts and references `muted` grey dashed · the
   "wrong variant" (principal root, printed coefficient) `rose` dashed.
10. **Every book equation is shown in full next to its number** (CLAUDE.md rule 3) — in markdown, derivation steps
    ("substitute (6.37), $u_\theta(a,\theta)=-2U\sin\theta-\Gamma/2\pi a$"), traps, recaps, and every explainer text.
    Builders reuse `show_eqs(text, EQ)` from `notebooks/build_ch05.py` with an `EQ` dict for all 110 labels of ch06 plus
    the earlier labels cited here ((3.7), (3.17), (3.18), (3.23), (4.7), (4.10), (4.12), (4.17), (4.38), (4.39b), (4.40),
    (4.72), (4.75), (4.83), (4.103), (4.106), (4.111), (5.2), (5.8), (5.11)); JS explainers use a local `showEqs`.
    `tools/eq_refs.py` must list 0 offenders. Exercise numbers are cited as "Exercise 6.49", never as bare equations.
11. **nbkit behaviour** (ch04/ch05 lesson): `nb.recap(...)` and `nb.section(...)` end the current CORE block, so every
    RECAP sits **before** the `nb.core` call of the block that uses it (R01 before C01; R02 R03 R04 R06 R07 R08 R09 before
    C02; R05 before C03; R10 R11 before C04; R12 R13 R14 before C06; R15 R16 R17 before C08; R18 before C10; R19 R20 R21
    before C12; R22–R30 before C13; R31 R32 before C15). Every `nb.derivation` sits inside its curation CORE block with
    `ref=` a bare equation label ("6.40"). Items taught in another block get a one-line `nb.pointer` in their own section
    (N90 is listed under §6.7 but taught in C14).
12. **Parsing.** Part E is the only part with table rows after its heading; Part F has no line starting with a table bar
    (absolute values are written \lvert … \rvert). Explainer headings are exactly `### E1 · superposition_sandbox` …
    `### E9 · added_mass_sphere` and `### B1 · flow_net_sources_vortices`. Primer terms in Part A are the exact Concept
    text of the Part E rows marked "primer" (coverage_check matches the first 18 characters). "Explained by" never names
    an earlier chapter's C-number (it cites "Ch. 4 §4.9" or a primers.md entry).
13. **Public repo.** No exercise text; Example 6.2's flow rate and grid values, Fig. 6.8's zero angle, Fig. 6.10's
    Reynolds number and the §6.1 Re threshold stay in `tests/book_values_ch06.json`. The notebook shows **our** numbers
    (Example 6.2 is run with Q = 1 m²/s and reported as ψ/Q), our derivations (D03, D08, D31 are ours by rule (b)) and a
    **qualitative** band (labelled so) for the measured cylinder pressure of Fig. 6.10.

Order of parts: C (contract) · A (notebook storyboard) · B (explainers) · D (runtime) · E (prerequisite ledger) · F
(derivations).

---

## Part C — functions the builders will call (the implementer's contract)

**Status column.** **§4 #n** = planned in `analysis/ch06.md` §4 row n (signature kept or refined here). **§8** = added
by the curation's notes for the implementer. **NEW** = added by this design (in neither) — flagged as the phase asks.
**reuse** = exists (ch01–ch05) and is only called. Docstrings cite § and Eq. (from the rendered page), list symbols with
units and assumptions, and carry the validation label. Return types: frozen dataclasses for elements and flows, `dict`
for state/scenario results (keys listed), `tuple` for (D, L) and (u, v).

### C.0 Reused (existing; called, not changed)
| # | Callable (signature) | Returns | Used by | Status |
|---|---|---|---|---|
| 0.1 | `style.setup_notebook() -> bool` (FAST) · `style.COLORS` · `style.savefig(fig, "ch06", name)` | FAST · palette · path | setup, every figure | reuse |
| 0.2 | `anim.animate(update, frames, fig, interval)` · `anim.show_animation(anim, player="video"\|"frames")` | HTML | C04, C08, C12, C15 | reuse |
| 0.3 | `interact.slider_figure(fn, name, values, *, unit, xlabel, ylabel, title, xrange, yrange, height)` · `interact.live(fn, **widgets)` | plotly Figure / widget | C05, C06, C07, C09, C11, C13, C14, C15 | reuse |
| 0.4 | `embed.show_viz("ch06", slug)` | display | 9 explainer cells | reuse |
| 0.5 | `tools.convergence.observed_order(h, err) -> float` | log–log slope | C03 (far field 1/R²), C04 (ε²), C12 (grid order), C14 (panels) | reuse |
| 0.6 | `core.streamfunction.velocity_from_streamfunction_2d(psi, x, y, rho=1.0, h=1e-5, **p)` · `velocity_from_streamfunction_axisym(psi, R, z, …)` · `flux_between_streamlines(psi, p1, p2, n)` | (u, v) · flux | R03, R24, N77 | reuse |
| 0.7 | `ch03.cylinder_flow(x, y, U, a, frame="body"\|"fluid", t=0.0)` · `ch03.cylinder_streamfunction(x, y, U, a)` | (u, v) · ψ | R12, R14 (parity with `pf.cylinder`) | reuse |
| 0.8 | `core.bernoulli`: `bernoulli_function(speed, p, z, rho, g, kind, p_o)`, `unsteady_bernoulli_pressure(dphi_dt, speed, z=0, rho, g, C)`, `viscous_irrotational_residual(u, x, t, mu, h)`, `dynamic_pressure(U, rho)` | floats | R05, R32, C01, C06 | reuse |
| 0.9 | `core.biot_savart`: `point_vortex_velocity(x, xv, Gamma)`, `wall_image_system(xv, Gamma, wall_y)`, `circle_image_system(xv, Gamma, a, …)`, `point_vortex_evolve(xv0, Gamma, t_eval, boundary="wall", **bp)` | arrays | R11, R15, R16, R17 (C08) | reuse |
| 0.10 | `core.vorticity.circle_loop_points(center, radius, n)` · `loop_circulation(u, pts, t)` | (2, n) · Γ | N38, B1 | reuse |
| 0.11 | `core.vortices.rankine_vortex(r, Gamma, sigma)` · `hill_stream_function(R, z, A, a, outside=True)` | (u_θ, ω_z) · ψ | C02 (Rankine), N83 (parity) | reuse |
| 0.12 | `ch04.accelerating_sphere_pressure(theta, a, dUdt, rho, U)` · `ch04.is_incompressible_regime(U, T, threshold)` · `ch05.kelvin_hypotheses_text(inviscid, barotropic, conservative, inertial)` | floats · bool · str | N101 (parity), R01 | reuse |
| 0.13 | `core.curvilinear.coordinates(system)`, `laplacian(f, system)`, `divergence(A, system)`, `curl(A, system)` (sympy) · `core.operators.laplacian(phi, h)` · `core.integral_theorems.gauss_legendre_nodes(lo, hi, n)` | sympy / arrays | R06–R09, R27–R30, D24, D25, D30 | reuse |
| 0.14 | `core.similarity` Reynolds helpers (`Scales`, `reynolds`) · `ch01.fluid_properties(name, T)` | floats | C01 (N01 numbers) | reuse |

### C.1 `fluidpy/core/potential.py` (NEW module, `core.PF` in the analysis) — superposable ideal-flow elements
Elements are `@dataclass(frozen=True)` objects with methods `w(z)`, `dwdz(z)` (= u − iv), `phi(x, y)`, `psi(x, y)`,
`velocity(x, y) -> (u, v)` and the attribute `singular_points` (complex tuple); `z0` is complex. Multivalued pieces
(ln) use the principal branch rotated by `cut_angle` so the cut can be put outside the plotted fluid.
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 1.1 | `Uniform(U: float, V: float = 0.0)` — w = (U − iV)z | (6.7) $\psi=-Vx+Uy$, (6.14) $\phi=Ux+Vy$ | m²/s, m/s | C02, C03, every body | §4 #2 |
| 1.2 | `Source(m: float, z0: complex = 0j, cut_angle: float = np.pi)` — w = (m/2π) ln(z − z0) | (6.15), (6.48) | m²/s | C02–C05, C08 | §4 #2 |
| 1.3 | `Vortex(Gamma: float, z0: complex = 0j, cut_angle: float = np.pi)` — **counterclockwise**, w = −(iΓ/2π) ln(z − z0) | (6.8), (6.47) | m²/s | C02, C04, C07 | §4 #2 |
| 1.4 | `Doublet(d_vec: tuple[float, float], z0: complex = 0j)` — w = −(d_x + i d_y)/(2π(z − z0)); `Doublet.from_book_scalar(d, z0=0j)` = dipole −d e_x | (6.29), (6.49) | m³/s | C04, C06 | §4 #2 |
| 1.5 | `Corner(A: float, n: float, cut_angle: float \| None = None, rotate: float = 0.0)` — w = A(z e^{−i·rotate})^n with the cut outside the wedge 0 ≤ θ ≤ π/n (default cut along θ = π/n + (2π − π/n)/2); ½ ≤ n ≤ 4 | (6.24)–(6.27), (6.46) | A in m^{2−n}/s | C04 (N21–N23), C09, E4 | §4 #2, #17; §8 |
| 1.6 | `Flow(elements: Sequence)` — sum; methods `w`, `dwdz`, `velocity(x, y) -> (u, v)`, `phi`, `psi`, `speed`, `cp(x, y, U)`, `pressure(x, y, rho=1.0, p_inf=0.0, U=None)`, `stagnation_points(guesses=None, box=None, tol=1e-12) -> np.ndarray[complex]` (complex Newton on dw/dz with deflation, residual-verified), `__add__` (Flow + element) | superposition (§6.2), (6.18) Bernoulli, (6.32) C_p | arrays | C03–C11 | §4 #2, #7 |
| 1.7 | `half_body(U: float, m: float) -> Flow` (ψ uses θ ∈ [0, 2π), convention 5) | (6.30), (6.31), (6.50) | Flow | C05, E1 | §4 #10 |
| 1.8 | `cylinder(U: float, a: float, *, Gamma_cw: float = 0.0, Gamma_ccw: float \| None = None) -> Flow` — w = U(z + a²/z) + (iΓ_cw/2π) ln(z/a); `velocity_polar(r, theta) -> (u_r, u_theta)` | (6.33), (6.34), (6.36), (6.37), (6.51), (6.52) | Flow | C06, C07, C10, E2, E5 | §4 #11 |
| 1.9 | `laplacian_residual(fn, x, y, h: float = 1e-3) -> ndarray` (4th-order 9-point stencil on a callable ψ or φ; NaN within 5h of a singular point) | (6.5), (6.12) | 1/s (for ψ in m²/s) | C02 (N05, N11, N26 table) | §4 #3 |
| 1.10 | `polar_velocity(fn, r, theta, kind: str = "psi", h: float = 1e-6) -> (u_r, u_theta)` · `polar_velocity_sym(expr, r, theta, kind="psi")` | (6.21), (6.22) | m/s | C02 (N18, N19), D09 step 4 | §4 #5 |
| 1.11 | `normal_velocity_on(flow, curve_pts: ndarray[complex], normals: ndarray[complex] \| None = None) -> float` (max \|u·n\| on a sampled closed curve; normals from the curve if None) | (6.16) | m/s | C03, C05, C06, C07, C11 | §4 #6 |
| 1.12 | `far_field_check(flow, R: float, n: int = 256, U: complex \| None = None) -> float` (max \|(u − iv) − U\| on \|z\| = R) | (6.17) | m/s | C03 (N17 slope) | §4 #6 |
| 1.13 | `pressure_coefficient(speed, U) -> ndarray` | (6.32) $C_p=1-\lvert\mathbf u\rvert^2/U^2$ | – | C05, C06, C13 | §4 #7 |
| 1.14 | `mirror(elements, wall: str = "y=0") -> list` (vortex → −Γ, source → +m, doublet reflected, uniform kept parallel) | image rules (Figs. 6.14–6.15), (6.41) | elements | C08 | §4 #14 |
| 1.15 | `circle_theorem(w_fn, a: float) -> Callable` (Milne-Thomson w(z) + conj(w(a²/conj z)); **our addition**, labelled) | R15 (named) | callable | C08 (R15 one line) | §4 #14 |
| 1.16 | `blasius_force(dwdz_fn, *, R: float \| None = None, center: complex = 0j, contour: ndarray[complex] \| None = None, rho: float = 1.0, n: int = 256) -> (D, L)` (periodic trapezoid on a circle or a closed counterclockwise polyline; raises `ValueError` if the polyline is clockwise) | (6.60) | N/m | C10, E5 | §4 #19; §8 |
| 1.17 | `laurent_coefficients(dwdz_fn, R: float, n: int = 64, center: complex = 0j, kmax: int = 4) -> dict[int, complex]` (FFT of samples on \|z\| = R: c₀, c₋₁ …) | N61, far field of (6.61) | 1/s·m^{k+1} | C10, E5 | §4 #19 |
| 1.18 | `AxisymUniform(U)`, `PointSource3D(Q, z0=0.0)`, `Doublet3D(d, z0=0.0)` (book dipole −d e_z) and `AxisymFlow(elements)` with `psi(R, z)`, `phi(R, z)`, `velocity_cyl(R, z) -> (u_R, u_z)`, `velocity_spherical(r, theta) -> (u_r, u_theta)` | (6.86), (6.87), (6.88), (6.75), (6.83) | m³/s, m/s | C13, C14 | §4 #24 |
| 1.19 | `sphere(U: float, a: float) -> AxisymFlow` (stream + doublet d = 2πa³U) | (6.89), (6.90) | AxisymFlow | C13, E8 | §4 #24 |
| 1.20 | `axisym_velocity_spherical(psi_fn, r, theta, h=1e-6) -> (u_r, u_theta)` · `axisym_velocity_spherical_sym(expr, r, theta)` | (6.83) | m/s | C13 (N79) | §4 #24 |
| 1.21 | `sphere_potential_vector(x, U_vec, d_vec) -> ndarray` (x shape (3,) or (3, N)) | (6.92) $\phi=(\mathbf U-\frac{\mathbf d}{4\pi\lvert\mathbf x\rvert^3})\cdot\mathbf x$ | m²/s | C13 (N85), D28 | §4 #24 |
| 1.22 | `moving_sphere_potential(x, xs, us, a) -> ndarray` | (6.97) | m²/s | C15, E9 | §4 #27 |
| 1.23 | `moving_sphere_velocity(x, xs, us, a) -> ndarray (3, N)` | (6.103) | m/s | C15 | §4 #27 |
| 1.24 | `moving_sphere_surface_pressure(e_xi, us, dus_dt, a, rho=1000.0, p_inf=0.0, split=False) -> ndarray \| dict(steady, acceleration, total)` (e_xi unit vectors (3,) or (3, N)) | (6.105) | Pa | C15, E9 | §4 #27; §8 (split) |
| 1.25 | `sphere_force_quadrature(p_fn, a, n_theta=32, n_phi=64) -> ndarray (3,)` (Gauss–Legendre in cos θ × periodic trapezoid in φ of −∮(p − p∞)n dA; `p_fn(e_xi)` returns p − p∞) | (6.98) | N | C15 | §4 #27 |
| 1.26 | `added_mass_sphere(a: float, rho: float) -> float` | (6.108) $M=\frac{2\pi a^3\rho}{3}$ | kg | C15, E9 | §4 #27 |

### C.2 `fluidpy/core/conformal.py` (NEW module, `core.CM`) — maps and flows carried through them
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 2.1 | `joukowski(zeta, b: float) -> complex ndarray` · `joukowski_derivative(zeta, b) -> complex` (dz/dζ = 1 − b²/ζ²) | (6.65) | m | C11, E6 | §4 #21 |
| 2.2 | `joukowski_inverse(z, b: float, branch: str = "outside") -> complex ndarray` — `"outside"`: ζ = ½[z + √(z − 2b)√(z + 2b)] (cut on the slit, \|ζ\| ≥ b asserted); `"principal"`: ½[z + √(z² − 4b²)] with numpy's principal root (**the wrong variant, kept callable for E6 and the wrong-variant test**) | (6.69) | m | C11, E6 | §4 #21; §8 |
| 2.3 | `mapped_flow(w_zeta, dw_dzeta, z_to_zeta, dzeta_dz) -> MappedFlow` with `w(z)`, `dwdz(z)` (chain rule), `velocity(x, y)`, `psi`, `phi`, `cp(x, y, U)` | u − iv = (dw/dζ)(dζ/dz) (after (6.69)) | object | C11, C10 (ellipse), E5, E6 | §4 #21 |
| 2.4 | `map_elements(f, dfdz, z0: complex, dz_list) -> ndarray[complex]` (images f′(z0)·δz) · `angle_preservation(f, dfdz, z0, dz1, dz2) -> (alpha, beta)` (angles in rad, before/after) | (6.63), (6.64) | rad | C11 (D19), E6 | §4 #20 |
| 2.5 | `grid_image(w_fn, xlim, ylim, n=200, levels=24) -> dict(X, Y, phi, psi, phi_levels, psi_levels)` (contour-ready Re/Im of w on a grid, cuts masked) | N66 (flow net as the image of a (φ, ψ) grid) | arrays | C11 figure | §4 #20 |

### C.3 `fluidpy/core/laplace_solvers.py` (NEW module, `core.LS`) — masked-grid Laplace and relaxation
Grids use the project layout `psi[j, i]` = (y, x); `mask[j, i]` is True at unknown (interior fluid) nodes; `bc` is an
array of the same shape holding the fixed values on non-mask nodes (NaN outside the domain).
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 3.1 | `laplacian_5pt(psi, dx, dy, mask=None) -> ndarray` | (6.70), (6.71) | 1/s | C12 (D22 order check) | §4 #22 |
| 3.2 | `node_update(psi, i, j, bc=None) -> float` (the average of the four neighbours) | (6.72) | m²/s | C12, E7 | §8 |
| 3.3 | `jacobi_sweep(psi, mask, bc) -> (psi_new, max_change)` · `gauss_seidel_sweep(psi, mask, bc, order="lex") -> (psi_new, max_change)` · `sor_sweep(psi, mask, bc, omega) -> (psi_new, max_change)` (one sweep each; lexicographic i fastest, j slowest — the book's order) | N73, (6.73) iterated | m²/s | C12, E7 | §8 |
| 3.4 | `residual_norm(psi, mask, dx=1.0, dy=1.0) -> float` (max \|ψ_{i−1,j} + ψ_{i+1,j} + ψ_{i,j−1} + ψ_{i,j+1} − 4ψ_{i,j}\|/4 over the mask) | (6.72) residual | m²/s | C12 | §4 #22 |
| 3.5 | `solve_laplace(mask, bc, method="gauss_seidel", tol=1e-10, max_iter=10**5, omega=None) -> (psi, history)`; `history` = dict(residual=[…], change=[…], sweeps=int); `method` ∈ {"jacobi", "gauss_seidel", "sor", "direct"} ("direct" = `scipy.sparse.linalg.spsolve`, history of length 0); `omega=None` for SOR → ω_opt = 2/(1 + sin(π/N)) | (6.72) solved | m²/s | C12, E7 | §4 #22 |
| 3.6 | `solve_poisson(mask, f, bc, dx, method="direct", …)` (for later chapters; one notebook line in C02's "what would change if") | ∇²ψ = f | m²/s | C02 pointer | §4 #22 |

### C.4 `fluidpy/core/panels.py` (NEW module, `core.PN`) — singularity-distribution methods
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 4.1 | `axial_singularity_solve(z_body, R_body, U: float, N: int \| None = None, z_nodes=None) -> dict(k, z_nodes, residual, cond, net_strength, A, rhs)` (segments from the body's nose to tail by default; collocation at the N body points) | the Fig. 6.29 system $\psi_m=-\sum_n\frac{k_n}{4\pi}(r^m_{n-1}-r^m_n)+\tfrac12UR_m^2=0$ | k [m²/s] | C14, E8 | §4 #25 |
| 4.2 | `axial_singularity_velocity(sol, R, z) -> (u_R, u_z)` · `axial_singularity_psi(sol, R, z) -> psi` | superposed (6.94) segments | m/s, m³/s | C14, E8 | §4 #25; NEW (`_psi`) |
| 4.3 | `source_panels(xb, yb, U=1.0, alpha=0.0) -> dict(lam, cp, xm, ym, normals, lengths)` (Hess & Smith constant-strength source panels, analytic influence, λ_i/2 self term; body points counterclockwise) | N90 (our extension) | λ [m/s] | C14 (N90), E8 panel mode | §4 #26 |
| 4.4 | `panel_velocity(sol, x, y) -> (u, v)` | N90 | m/s | C14 figure | §4 #26 |

### C.5 `fluidpy/ch06_ideal_flow.py` (chapter module; re-exports C.1–C.4)
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 5.1 | `ideal_flow_residuals(u_fn, p_fn=None, x=(0.2, 0.15), t=0.0, rho=1000.0, mu=1e-3, h=1e-4, **p) -> dict(continuity, euler, viscous_force, vorticity)`; `u_fn` may be a callable `u(x, t)` or a preset `"cylinder"` (U = 1, a = 0.1, p from Bernoulli) / `"poiseuille"` (plane Poiseuille, G = −dp/dx = 1 Pa/m, gap 2h₀ = 0.02 m) / `"corner"` (u = (2Ax, −2Ay)) | (6.1), (4.40) | 1/s, N/m³ | C01 | §4 #1 |
| 5.2 | `ideal_flow_applicability(Re: float, M: float = 0.0, baroclinic: bool = False, region: str = "outer") -> dict(ok: bool, reasons: list[str], verdict: str)` (regions "outer", "boundary_layer", "wake", "pipe"; reuses `ch04.is_incompressible_regime` logic via M and `ch05.kelvin_hypotheses_text(barotropic=not baroclinic)`) | §6.1, R01 | – | C01 | §4 #1 |
| 5.3 | `vorticity_from_psi(psi_fn, x, y, h=1e-3, **p) -> ndarray` (−∇²ψ, 9-point 4th order; `psi_fn` may be `"rankine"` with `Gamma`, `a`) · `vorticity_from_psi_sym(expr, x, y)` | (6.4) | 1/s | C02, B1 | §4 #3 |
| 5.4 | `rankine_vortex_psi(x, y, Gamma=1.0, a=0.1) -> ndarray` (ψ = −Γr²/(4πa²) inside, −(Γ/2π)[ln(r/a) + ½] outside; continuous with continuous slope) | (6.4) test field (Ch. 3 §3.5 Rankine) | m²/s | C02 figure, B1 | NEW |
| 5.5 | `delta_flux_check(fn, center=(0.0, 0.0), radii=(0.01, 0.1, 1.0, 10.0, 100.0), n=256, kind="psi", **p) -> ndarray` (∮∇f·n ds on each circle by the periodic trapezoid; `fn` a callable or `"vortex"` (keyword `Gamma`) / `"source"` (keyword `m`)) | (6.6), (6.13), D03 | m²/s | C02, B1 | §4 #4 |
| 5.6 | `orthogonality_check(flow, pts) -> float` (max \|∇φ·∇ψ\|/\|u\|² at the points; `flow` a Flow or a spec list) | (6.10), D04 | – | C02 | §4 #5 |
| 5.7 | `harmonic_polynomials(degree: int) -> list[tuple[sympy.Expr, sympy.Expr]]` (Re and Im of z^n for n ≤ degree, each asserted harmonic) | N20, (6.24)–(6.27) | sympy | C04 | §4 #8 |
| 5.8 | `source_sink_pair(x, y, m, eps) -> (phi, u, v)` · `doublet_limit_error(eps_list, d=2.0, pts=((1.0, 0.0), (0.0, 1.0), (0.7, 0.7))) -> ndarray` (max relative φ error vs (6.29) at 2mε = d) · `doublet_limit_frames(eps_list, d=2.0, n=161) -> list[dict(X, Y, psi, eps, err)]` | (6.28), (6.29) | m²/s | C04 (+ animation) | §4 #9; §8 (frames) |
| 5.9 | `half_body_shape(U, m, theta) -> (x, y)` (θ ∈ (0, 2π)) · `half_body_numbers(U, m) -> dict(a, h_max, psi_body)` · `half_body_surface_cp(theta) -> ndarray` · `half_body_cp_zero_angle() -> float` (rad, `brentq` on tan θ + 2(π − θ) in (π/2, π)) | (6.31), D07, D08 | m, – | C05, E1 | §4 #10 |
| 5.10 | `half_body_net_force(U, m, x_end, rho=1.0, n=4001) -> dict(D, L)` (pressure force on the body from the nose to x = x_end; → 0 as x_end grows) | N29 (Exercise 6.13 idea) | N/m | C05 | NEW |
| 5.11 | `rankine_oval(U, m, a) -> dict(half_length, half_width, psi_fn, stagnation)` (source +m at −a, sink −m at +a) | Exercise 6.19 preset (our derivation) | m | E1 preset, C05 "what would change" | §8 |
| 5.12 | `flow_from_spec(spec: list[dict]) -> Flow` (kinds "uniform" (U, V), "source" (m, x, y), "sink" (m, x, y), "vortex" (Gamma, x, y, ccw), "doublet" (dx, dy, x, y), "corner" (A, n)) | superposition | Flow | E1–E4 parity | NEW |
| 5.13 | `superposition_state(spec: list[dict], x: float, y: float, U_ref: float \| None = None) -> dict(u, v, psi, phi, cp, u_parts, v_parts, stagnation, net_source, closed, psi_dividing)` (`stagnation` a list of [x, y]; `closed` = abs(Σm) < 1e-12 and a body exists) | E1 terms/status | m/s, m²/s | E1 | §8 |
| 5.14 | `cylinder_surface_cp(theta, U=1.0, a=1.0, *, Gamma_cw=0.0) -> ndarray` · `cylinder_stagnation_points(U, a, *, Gamma_cw) -> ndarray[complex]` (closed form: two surface roots, one merged, or the off-body r₊) · `cylinder_surface_pressure(theta, U, a, *, Gamma_cw=0.0, rho=1.2, p_inf=0.0) -> ndarray` · `lift_per_span(rho, U, *, Gamma_cw) -> float` | (6.35), (6.37)–(6.40) | –, m, Pa, N/m | C06, C07, E2 | §4 #11 |
| 5.15 | `surface_pressure_force(flow, contour_pts: ndarray[complex], rho=1.2, p_inf=0.0, U=None) -> (D, L)` · `contour_force(p_fn, pts) -> (D, L)` (p_fn(z) → Pa; counterclockwise asserted) · `complex_force_from_pressure(p_fn, pts) -> complex` (D − iL) | (6.55)–(6.57) | N/m | C06, C07, C10 | §4 #12 |
| 5.16 | `circulation_family_check(U, a, Gammas_cw, R_far=None) -> dict(max_normal, far_field, circulation)` (arrays, one per Γ; far field measured on r = R_far, default 1000a) | N38 | m/s, m²/s | C07 | §4 #13 |
| 5.17 | `cylinder_circulation_state(U=10.0, a=0.1, *, Gamma_cw=2.0, rho=1.2, p_inf=0.0, n=128) -> dict(theta1_deg, theta2_deg, r_free, regime, D, L, L_KJ, speed_top, speed_bottom, cp_min)` (`regime` ∈ {"two surface points", "merged at the bottom", "free stagnation point"}; the θ's NaN when off the body) | E2 status/bars, C07 slider | °, m, N/m, m/s | E2, C07 | §8 |
| 5.18 | `two_sources(m, a) -> Flow` (via `mirror`) · `two_source_streamline(psi, m, a, x) -> ndarray` (y ≥ 0 of $x^2-y^2-2xy\cot(2\pi\psi/m)=a^2$) | (6.41), (6.53) | m | C08, E3 | §4 #14 |
| 5.19 | `example_6_1(t, Gamma=1.0, h=1.0, rho=1000.0, p_inf=0.0, route="closed") -> dict(xi_x, xi_y, p_origin, v_origin, dphidt, unsteady, speed_part)` (`route="numeric"`: moving `Vortex` + image, ∂φ/∂t by central difference in t, path from `point_vortex_evolve`) | Example 6.1 | m, Pa, m/s, m²/s² | C08, E3 | §4 #15 |
| 5.20 | `example_6_1_wall_pressure(y_wall, t, Gamma=1.0, h=1.0, rho=1000.0, p_inf=0.0, split=False) -> ndarray \| dict(unsteady, speed, total)` (p − p∞ along the wall x = 0; our extension of the origin value) | Example 6.1 extended | Pa | C08 animation, E3 | §8 |
| 5.21 | `cauchy_riemann_residual(w_fn, z, h=1e-6) -> (r1, r2)` (φ_x − ψ_y, φ_y + ψ_x by differences along x and iy; `w_fn` may be `"conj"` for the non-analytic control z*) · `cauchy_riemann_residual_sym(expr, x, y)` | (6.44) | 1/s | C09, E4 | §4 #16 |
| 5.22 | `corner_speed_exponent(n) -> float` (= n − 1 = (π − α)/α) | (6.46) | – | C09, E4 | §8 |
| 5.23 | `complex_potential_probe(kind, x, y, A=1.0, n=2.0, h=1e-4, **p) -> dict(phi, psi, dwdz_re, dwdz_im, u, v, qx_re, qx_im, qy_re, qy_im, cr1, cr2, speed)` (kinds "corner", "source", "vortex", "doublet", "cylinder", "conj"; `qx`, `qy` = difference quotients along x and along iy with step h) | (6.42)–(6.46), D14 | m/s | E4, C09 | NEW |
| 5.24 | `cv_force_on_body(flow, R_outer, rho=1.2, p_inf=0.0, U=None, n=512) -> (D, L)` (momentum flux + pressure on a large circle; Exercise 6.27 idea) | (6.54) route | N/m | C10 | §4 #18 |
| 5.25 | `kutta_zhukhovsky_sym() -> dict(square, coeff_z1, coeff_z2, coeff_z2_printed, residue, D, L)` (sympy: expands (U + iΓ/2πz − d/2πz²)², returns the true 1/z² coefficient beside the printed one) | (6.61), (6.62) | sympy | C10 (D18 check) | §4 #19 |
| 5.26 | `laurent_contributions(body="cylinder", R=0.2, *, Gamma_cw=2.0, U=10.0, a=0.1, b=None, kmax=4, rho=1.2) -> dict(powers, contrib_re, contrib_im)` (the share of each power z^k of (dw/dz)² in (iρ/2)∮(dw/dz)² dz; only k = −1 is nonzero) | (6.61) | N/m | E5, C10 figure | §8 |
| 5.27 | `blasius_state(body="cylinder", *, Gamma_cw=2.0, U=10.0, a=0.1, R=0.2, rho=1.2, n=256, b=None, alpha=0.0, offset_x=0.0, offset_y=0.0) -> dict(D, L, L_KJ, crosses_body, c0_re, cm1_re, cm1_im, D_pressure, L_pressure)` (bodies "cylinder", "ellipse" (Zhukhovsky, needs b), "tilted_ellipse", "rankine_oval_vortex") | (6.56), (6.60), (6.62) | N/m | E5 | NEW |
| 5.28 | `cot_flow(z) -> complex` (u − iv = cot z for w = ln sin z) | N66 | m/s | C11 | §4 #20 |
| 5.29 | `joukowski_ellipse(a, b) -> dict(A, B, foci)` · `elliptic_cylinder_flow(U, a, b, *, Gamma_cw=0.0, alpha=0.0) -> MappedFlow` | (6.66)–(6.69) | m | C10, C11, E5, E6 | §4 #21 |
| 5.30 | `joukowski_state(a=1.2, b=1.0, U=1.0, *, Gamma_cw=0.0, x=-3.0, y=0.5, branch="outside") -> dict(zeta_re, zeta_im, zeta_abs, u, v, A, B, foci, inside)` | (6.65)–(6.69) | m, m/s | E6 | NEW |
| 5.31 | `four_point_system(boundary="xy", values=None, as_grid=False) -> dict(A, b, psi)` or, with `as_grid=True`, `(mask, bc)` for `ls.solve_laplace` (ordering ψ₂₂, ψ₃₂, ψ₂₃, ψ₃₃; `boundary="xy"` fills ψ^B = (i − 1)(j − 1); `values` a dict {(i, j): ψ} overrides) | (6.73) | m²/s | C12, E7 | §4 #22 |
| 5.32 | `example_6_2(Q=1.0, n_iter=None, tol=1e-10, method="gauss_seidel", refine=1, omega=None) -> dict(psi, X, Y, mask, history, grid_shape, probe)` (`probe` = the (j, i) index of a fixed interior point, the same physical point at every `refine`; the contraction geometry of Fig. 6.24 on the book grid ×refine; our Q; book Q and grid values only in the private JSON) | Example 6.2 | m²/s | C12, E7 | §4 #22 |
| 5.33 | `relaxation_state(problem="four_point", method="gauss_seidel", sweeps=1, omega=None, refine=1) -> dict(residual, max_change, psi22, psi32, psi23, psi33, psi_probe, sweeps_to_tol)` (problems "four_point", "contraction", "xy_square") | N73, (6.73) | m²/s | E7 | NEW |
| 5.34 | `stokes_operator_residual(psi_fn, R, z, h=1e-4) -> ndarray` · `stokes_operator_residual_sym(expr, R, z)` · `axisym_flux_between(psi_fn, P1, P2, n=64) -> (quad, two_pi_dpsi)` · `axisym_laplacian_residual(phi_fn, R, z, h=1e-4)` | (6.77), (6.78), (6.80) | 1/s, m³/s | C13 | §4 #23 |
| 5.35 | `sphere_surface_cp(theta) -> ndarray` · `cylinder_vs_sphere(r_over_a) -> dict(cyl_pert, sph_pert, cyl_cp_min, sph_cp_min, cyl_speed_max, sph_speed_max)` | (6.91), (6.35) | – | C13 | §4 #24; NEW (`cylinder_vs_sphere`) |
| 5.36 | `line_sink_stream_function(R, z, k, a, method="closed") -> ndarray` (`"quad"` integrates (6.93)) · `airship(U, Q, a) -> dict(psi_fn, stagnation_front, stagnation_rear, length, R_max, contour, closure)` | (6.93)–(6.95) | m³/s, m | C14, E8 | §4 #25 |
| 5.37 | `axisym_body_target(kind="rankine_oval", N=40, **p) -> (z_body, R_body)` (kinds "airship", "rankine_oval", "sphere", "ellipsoid" with `fineness`) · `panel_cp_error(N, body="circle") -> float` | E8 targets, N90 convergence | m, – | C14, E8 | §8 |
| 5.38 | `axial_state(target="rankine_oval", N=20, U=1.0, fineness=3.0) -> dict(cond, net_strength, body_error, k_max, k_min)` (`cond` is the **1-norm** condition number `np.linalg.cond(A, 1)`, also in `axial_singularity_solve`, so JS and Python agree) | Fig. 6.29 system | –, m²/s | E8 | NEW |
| 5.39 | `moving_sphere_dphidt_sym() -> dict(direct, chain, printed_bracket, correct_bracket)` (sympy re-run of (6.101)–(6.104)) | (6.101)–(6.105) | sympy | C15 (D29 check) | §4 #27 |
| 5.40 | `added_mass_by_energy(a, rho, U=1.0, method="surface") -> float` (`"volume"` integrates ½ρ\|∇φ\|² over r > a with `quad`) | D31 (Exercise 6.49) | kg | C15, E9 | §4 #27 |
| 5.41 | `cylinder_added_mass(a, rho) -> float` (ρπa² per depth; force and energy routes in the tests) | N104 (Exercise 6.45) | kg/m | C15 | §4 #29 |
| 5.42 | `sphere_motion(m, a, rho=1000.0, F_E=None, t_eval=None, g=None, mode=None, amplitude=0.05, omega=2.0, added_mass=True) -> dict(t, x, u, du_dt, F_s)` (modes "ball" / "bubble" (buoyancy − weight with g), "oscillating" (prescribed x_s = amplitude sin ωt, returns F_s), None (F_E callable)) | (6.109) | s, m, m/s, N | C15, E9 | §4 #28; §8 |
| 5.43 | `added_mass_state(a=0.1, rho=1000.0, us=1.0, dus=1.0, angle_deg=0.0, theta_s_deg=0.0) -> dict(p_steady, p_accel, p_total, F_steady, F_accel, M, M_energy, cp_steady)` (angle between u_s and du_s/dt; θ_s = the surface point's angle from u_s in the plane of both; F_accel = the component along du_s/dt) | (6.105), (6.108) | Pa, N, kg | E9 | NEW |

### C.6 `scripts/ch06_*.py` (runnable demos; drawing helpers carry no physics)
| # | Callable | Purpose | Status |
|---|---|---|---|
| 6.1 | `scripts/ch06_drawings.py`: `flow_net(ax, flow, xlim, ylim, n=241, n_levels=21, mask_body=None)` (ψ solid, φ dashed, body black), `draw_body(ax, pts, **kw)`, `pressure_arrows(ax, pts, p, normals, scale)`, `separated_cp_band(theta_front_deg, sep_deg=80.0) -> (lo, hi)` (**qualitative sketch band**, labelled so in every figure) | notebook and scripts figures | NEW |
| 6.2 | `scripts/ch06_half_body.py`, `ch06_cylinder_cp.py`, `ch06_conformal_grid.py`, `ch06_contraction.py`, `ch06_airship.py`, `ch06_superposition_gallery.py`, `ch06_cylinder_circulation.py`, `ch06_added_mass.py`, `ch06_panels_convergence.py` — each `main()` writes `outputs/ch06/<name>.png` and asserts its headline number | DEMO rows | §4 #30 |

**Flags for the implementer** (things in neither the analysis nor the curation, or corrected here): NEW functions 5.4,
5.10, 5.12, 5.23, 5.27, 5.30, 5.33, 5.35 (`cylinder_vs_sphere`), 5.38, 5.43, 4.2 (`axial_singularity_psi`), 6.1; the
half-body θ range [0, 2π) (convention 5); ∇ψ = e_z × ∇φ (convention 8; do not assert the curation's form); `Doublet`
complex form w = −(d_x + i d_y)/2π(z − z0); `blasius_force` rejects clockwise polylines; `example_6_2` default Q = 1 m²/s
(ours); `solve_laplace(method="sor", omega=None)` uses ω_opt = 2/(1 + sin(π/N)); parity keys listed in Part B must exist
exactly as named.
---

## Part A — notebook storyboard (`notebooks/build_ch06.py` → `notebooks/ch06_ideal_flow.ipynb`)

**One line per book section** (cell numbers are estimates for the builder's budget; ≈ 560 cells in all):
- §6.1 → R01, C01 (N01 N02 N03 N04 · D01) — cells ≈ 8–30
- §6.2 → R02 R03 R04 R06 R07 R08 R09, C02 (N05–N15 N18 N19 · D02 D03 D04 · primer P149), R05, C03 (N16 N17 · D05) — cells ≈ 31–105
- §6.3 → R10 R11, C04 (N20–N26 · D06 · P150 · animation), C05 (N27–N29 · D07 D08 · E1), R12 R13 R14, C06 (N30 N31 · D09 · P151), C07 (N32–N38 · D10 D11 · P152 · E2), R15 R16 R17, C08 (N39 N40 · D12 D13 · animation · E3) — cells ≈ 106–280
- §6.4 → C09 (N41–N52 · D14 D15 · P153 P154 P155 · E4) — cells ≈ 281–325
- §6.5 → R18, C10 (N53–N63 · D16 D17 D18 · P156 P157 P158 · E5) — cells ≈ 326–375
- §6.6 → C11 (N64–N70 · D19 D20 D21 · P159 · E6) — cells ≈ 376–415
- §6.7 → R19 R20 R21, C12 (N71–N74 · D22 D23 · P160 P161 · animation · E7), pointer to N90 — cells ≈ 416–460
- §6.8 → R22–R30, C13 (N75–N85 · D24 D25), C14 (N86–N90 · D26 D27 · P162 · E8) — cells ≈ 461–520
- §6.9 → R31 R32, C15 (N91–N105 · D28 D29 D30 D31 · P163 P164 · animation · live · E9) — cells ≈ 521–555
- §6.10 → N106, S01–S07, summary — cells ≈ 556–562

Every CORE block below follows the order: problem in plain words → idea → primers → maths (notes and derivations,
Part F) → tiny example → code (fluidpy) + "What does the code above do?" → from-scratch check → visual(s) → notes and
"What would change if…". Code drafts are intent + exact calls; the builder writes every line commented (novice grade,
units, equation numbers written out). *expect* gives the numbers the executed cell must print (sanity values computed
for this design).

### A.0 Front matter
1. `nb.title(big_idea=…, roadmap=[…15…], prerequisites=[…])`. **Big idea (draft):** "Air over a wing, water past a
   bridge pier, the ocean round an island: far from walls and wakes these fluids move almost without friction and
   without spin. Then two scalar functions — the stream function ψ and the velocity potential φ — obey the *linear*
   Laplace equation, so flows can be *added*: a stream plus a source makes a blunt nose, a stream plus a doublet makes a
   cylinder, a vortex added to it makes lift. Complex numbers pack ψ and φ into one function w(z), turn force
   calculations into a residue and let a map carry the circle's flow onto other shapes. A grid of averages solves the
   same Laplace equation numerically, sources on an axis draw airships, and a moving sphere shows the one force an
   ideal fluid *does* exert: added mass." **Roadmap (one line per CORE):** C01 ideal-flow equations and their limits ·
   C02 ψ, φ and Laplace with point singularities · C03 superposition: a streamline can be a wall · C04 the element kit
   and the doublet · C05 the half-body · C06 the cylinder and d'Alembert's paradox · C07 circulation and lift ρUΓ · C08
   images and the wall pressure under a passing vortex · C09 the complex potential · C10 Blasius and Kutta–Zhukhovsky ·
   C11 conformal maps and the Zhukhovsky ellipse · C12 finite-difference Laplace and Gauss–Seidel · C13 axisymmetric flow
   and the sphere · C14 bodies from axial sources · C15 the accelerating sphere and added mass. **Prerequisites:**
   continuity and the stream function (Ch. 4 §4.2–4.3) · Bernoulli, steady and unsteady (Ch. 4 §4.9) · vorticity,
   circulation and Kelvin's theorem (Ch. 3 §3.4, Ch. 5 §5.2) · the point vortex and its image (Ch. 5 §5.7) · the
   Laplacian, Gauss' theorem, the Dirac delta (Ch. 2) · the cylinder flow in two frames (Ch. 3 §3.3).
2. `nb.explainer_index([...])` — 9 rows: ("superposition_sandbox", "Where does the body come from?", "add elements; the
   stagnation streamline becomes a wall; closed only when sources cancel") · ("cylinder_circulation_lift", "How does spin
   turn into lift?", "stagnation points slide, L = ρUΓ, drag stays zero") · ("vortex_wall_images", "What does a wall feel
   as an eddy passes?", "image vortex, drift Γ/4πh, suction then over-pressure") · ("complex_potential_corners", "Why is
   the velocity u − iv?", "Cauchy–Riemann as two difference quotients; calm vs violent corners") ·
   ("blasius_kutta_contour", "Why doesn't lift depend on the shape?", "the contour can move; only U × Γ/z survives") ·
   ("conformal_joukowski", "How does a map carry a flow?", "angles kept, circle → ellipse, the right square root") ·
   ("laplace_relaxation", "How does a grid solve Laplace's equation?", "average of neighbours, Jacobi vs Gauss–Seidel vs
   SOR, stop on the residual") · ("axial_singularity_bodies", "Given a shape, which sources draw it?", "N unknown
   strengths, ψ = 0 at N points, one linear solve") · ("added_mass_sphere", "Why does an ideal fluid resist
   acceleration?", "speed pressure cancels, acceleration pressure gives half the displaced mass").
3. `nb.setup()`.
4. `nb.code` — **chapter imports** (outside any block): `import numpy as np` · `import matplotlib.pyplot as plt` ·
   `import sympy as sp` · `from fluidpy import ch06_ideal_flow as ch06` · `from fluidpy.core import potential as pf,
   conformal as cm, laplace_solvers as ls, panels as pn` · `from fluidpy import ch03_kinematics as ch03,
   ch04_conservation_laws as ch04, ch05_vorticity_dynamics as ch05` · `from fluidpy.core.interact import
   slider_figure, live` · `from fluidpy.core.anim import animate` · `from fluidpy.core.style import savefig` ·
   `import sys; sys.path.insert(0, "scripts"); from ch06_drawings import flow_net, draw_body, pressure_arrows,
   separated_cp_band`. *explain:* one line per import ("the chapter module re-exports the four new core modules").
5. `nb.md` — **⚠️ Conventions in this chapter** (a two-column table, then numbers): Γ counterclockwise (project, (6.8)
   $\psi=-\frac{\Gamma}{2\pi}\ln r$) vs the book's clockwise Γ of (6.36)–(6.40) (cylinder: Γ_cw = 2 m²/s ⇔ Γ_ccw = −2 m²/s,
   stagnation points at −9.16° and −170.84°, L = +24 N/m for U = 10 m/s, a = 0.1 m, air); the 2-D dipole vector points
   from sink to source, cylinder $\mathbf d=-2\pi Ua^2\mathbf e_x$; θ from +x (the nose of a body in a stream from the
   left is at θ = π); in §6.8 z is horizontal along the stream; in §6.9 w is a velocity component; D, L are forces on the
   body; ψ in m²/s (plane) vs m³/s (Stokes); B is the span, M the added mass (not Mach). "We compute in the project
   convention and show the book's alongside wherever they differ."
6. `nb.md` — **🔁 Tools from earlier chapters used in this one** (one line each, primer number): partial derivatives
   (P25), Taylor series (P26/P98), definite integral (P27), net pressure force −∮p n dA (P28), `lambda` (P29), trapezoid
   (P37), product rule (P38), sympy (P40), i² = −1 and Euler's formula (P45), `np.where` (P46), chain rule (P49, P91),
   `np.linalg.solve` (P57), orders of smallness (P68), Vieta's formulas (P71), `np.arctan2` (P70), directional
   derivative (P75), `np.meshgrid` (P76), broadcasting (P77), contour/streamplot (P78), eigenvalues (P80), complex
   conjugate (P81), `quad` (P87), polar and spherical unit vectors (P88, P105), parametric curves (P92), substitution in
   an integral (P106), `brentq` (P108), dataclasses (P111), momentum flux (P114), Schwarz (P121), curl of a curl (P122),
   periodic trapezoid (P136), Poisson and Green's function (P139), gradient of 1/distance (P140), FFT (P142),
   Gauss–Legendre (P143), improper integral (P144), animate (P16), slider_figure (P17), show_viz (P18), `solve_ivp`
   (P31/P94), log–log slopes (P13), `assert np.allclose` (P15). Each block repeats the ones it uses in a one-line
   reminder at first use.

### A.1 §6.1 Relevance of Irrotational Constant-Density Flow Theory — R01, C01 (+N01 N02 N03 N04 · D01)
1. `nb.section("6.1", "Relevance of Irrotational Constant-Density Flow Theory", intro="**What is this section about?**
   Real fluids are viscous, so why study a theory with no friction? Because at high Reynolds number viscosity only
   matters in thin layers next to walls and in wakes; everywhere else the flow keeps no spin and obeys a much simpler
   pair of equations. This section says what those equations are, what they drop, and where they may be trusted.")`
2. `nb.recap("R01", "Why the fluid stays irrotational — Kelvin", "In a fluid of constant density pushed only by
   conservative forces, and away from viscous regions, the circulation of every material loop is constant, $D\Gamma/Dt=0$
   *(Eq. 5.8)* — so fluid that starts without vorticity (a uniform stream far upstream) stays without it until it enters
   a boundary layer, a wake or a separated region. Constant density needs a low Mach number $M=U/c\ll1$ *(Eq. 4.111)* and
   no baroclinic torque from a density field (Ch. 5 §5.2 lists the four restrictions).", where="Ch. 5 §5.2")` + `nb.code`:
   `print(ch05.kelvin_hypotheses_text(barotropic=False))` (Kelvin fails in a stratified fluid) · `print(
   ch04.is_incompressible_regime(10.0))` (a car at 10 m/s: M ≈ 0.03 → True). *expect:* a sentence naming the baroclinic
   term; `True`.
3. `nb.core("C01", "The ideal-flow equations $\\nabla\\cdot\\mathbf u=0$ and $\\rho\\,D\\mathbf u/Dt=-\\nabla p$ (6.1) — and
   where to trust them", question="Water and air are viscous. How can a theory with no friction describe them at all —
   and where does it fail?")`
4. `nb.md` — **The problem in plain words:** "Put your hand out of a car window at 10 m/s: the air feels the car only
   through a skin of slowed air about a millimetre thick on its surface; a few centimetres away the air moves as if it
   had no viscosity at all. Engineers design wings, hulls and cars from that outer flow; oceanographers use it for flow
   round islands and seamounts. We want the equations of that outer flow and an honest list of where they fail."
5. `nb.md` — **The idea** (ASCII):
   ```
   far from walls:   no spin (ω = 0)  ─►  viscous force μ∇²u = −μ∇×ω = 0   ─►  Euler + ∇·u = 0 = (6.1)
   near the wall:    thin boundary layer, ω ≠ 0, friction matters             ─►  Ch. 9
   behind blunt body: separated wake, ω ≠ 0 everywhere                        ─►  ideal flow fails
   ```
   "**Irrotational does not mean inviscid — it means the viscous forces cancel.**"
6. `nb.md` — 🔁 reminders: curl of a curl (P122), product rule for a divergence (P113 — Ch. 4 §4.2), hydrostatic split
   (Ch. 4 §4.9), boundary conditions (P20).
7. `nb.derivation("D01", …)` — Part F D01 (8 steps), ref "6.1".
8. `nb.note` — **N02 [B]** "**Only one boundary condition survives.** (6.1) contains only first derivatives of u, so a
   wall can impose one condition, no flow through it:" equation `\mathbf n\cdot\mathbf u=\mathbf n\cdot\mathbf U_s\quad(\text{kept}),\qquad \mathbf t\cdot\mathbf u=\mathbf t\cdot\mathbf U_s\quad(\text{dropped})`.
   "The fluid slides along the wall. D01 step 8 showed why: the term that needed a second condition, μ∇²u, is zero.
   The real fluid obeys no-slip, which is why a boundary layer must exist (Ch. 9)."
9. `nb.worked_example("a corner flow that is viscous but feels no viscous force", "Take $\\phi=A(x^2-y^2)$ with A = 1
   s⁻¹ (a stagnation flow, met again in C04). 1. $u=\\partial\\phi/\\partial x=2x$, $v=\\partial\\phi/\\partial y=-2y$ [m/s].
   2. $\\nabla\\cdot\\mathbf u=2-2=0$ ✓ incompressible. 3. $\\omega_z=\\partial v/\\partial x-\\partial u/\\partial y=0-0=0$ ✓
   irrotational. 4. $\\nabla^2u=\\partial^2(2x)/\\partial x^2+\\partial^2(2x)/\\partial y^2=0$, likewise $\\nabla^2v=0$: no net viscous
   force. 5. But the element is being stretched: the viscous normal stress $2\\mu\\,\\partial u/\\partial x=2\\times10^{-3}\\times2=
   0.004$ Pa in water is **not** zero — equal and opposite on the two faces of each element, so it cancels.")`
10. `nb.code` — *code:* `cases = [("corner", (0.2, 0.15)), ("cylinder", (0.2, 0.15)), ("poiseuille", (0.2, 0.005))]`
    (Poiseuille probed inside its 2 cm gap) · `for name, pt in cases: r = ch06.ideal_flow_residuals(name, x=pt);
    print(name, {k: float(np.max(np.abs(v))) for k, v in r.items()})`. *expect:* corner and cylinder: continuity,
    viscous_force, vorticity ≲ 1e-8 (stencil round-off), euler ≲ 1e-6; poiseuille: continuity 0, vorticity 5.0 s⁻¹
    (ω_z = −∂u/∂y = Gy/μ = 1 × 0.005/10⁻³), viscous_force of size 1.0 N/m³ (μ∇²u = −G e_x). *explain:* 1. three velocity
    fields as presets; 2. the residuals of continuity, of Euler with p from Bernoulli, the viscous force μ∇²u and the
    vorticity, by central differences; 3. the potential flows pass all four, the Poiseuille flow (rotational) has a net
    viscous force that balances the driving pressure gradient.
11. `nb.check_agree` — **from scratch (curation §7):** second central differences of u at one point: `h = 1e-3`;
    `u = lambda X, Y: ch06.pf.cylinder(1.0, 0.1).velocity(X, Y)[0]` (x-velocity of the cylinder flow); `lap = (u(x+h, y)
    + u(x-h, y) + u(x, y+h) + u(x, y-h) - 4*u(x, y))/h**2`; `mu = 1e-3`; `assert np.allclose(mu*lap,
    ch06.ideal_flow_residuals("cylinder", x=(x, y))["viscous_force"][0], atol=1e-8)`; same for Poiseuille (≈ −1 N/m³).
    Markdown: "Our five-point sum gives the library's viscous force: zero for the potential flow, −G for Poiseuille."
12. `nb.note` — **N01 [B]** "**Where (6.1) works** — with $\mathrm{Re}=\rho UL/\mu$ *(Eq. 4.103)* large, viscosity and spin
    stay in thin layers." Table | ideal flow predicts | it does not predict |: velocity away from walls · skin friction;
    pressure forces normal to a thin attached boundary layer · dissipation; streamline patterns that minimise form drag ·
    pipe and duct flow; unsteady inertia (added mass, C15) · wakes, turbulence, low-Re flow. "Number: a car, U = 10 m/s,
    L = 1 m, air ν = 1.5×10⁻⁵ m²/s: Re = UL/ν ≈ 6.7×10⁵; a laminar boundary layer is about L/√Re ≈ 1.2 mm thick." +
    `nb.code`: `for case in [dict(Re=6.7e5), dict(Re=50), dict(Re=6.7e5, M=0.5), dict(Re=6.7e5, baroclinic=True),
    dict(Re=6.7e5, region="wake")]: print(case, ch06.ideal_flow_applicability(**case)["verdict"])`. *expect:* ✅ outer
    flow · ⚠️ boundary layer not thin · ⚠️ compressible · ⚠️ Kelvin fails (baroclinic) · ⛔ wake (rotational).
13. `nb.figure` — two panels (7 × 3.2 in): left, the cylinder flow (U = 1 m/s, a = 0.1 m) with `flow_net`, coloured by
    log₁₀\|μ∇²u\| from `ideal_flow_residuals` on a 61 × 41 grid (round-off level ~1e-10, masked inside the body); right,
    plane Poiseuille in a 2 cm gap with u(y) arrows and a flat colour at log₁₀ 1 = 0 (|μ∇²u| = G = 1 N/m³). Title
    "No spin, no net viscous force". *see:* "left, a uniform round-off colour; right, the viscous force the pressure
    gradient must balance." *read:* "a potential flow has μ∇²u = −μ∇×ω = 0 everywhere, although the fluid is viscous;
    Poiseuille flow is all spin." *change:* "…μ multiplied by 100: the right panel's force stays G (the profile
    flattens by 100), the left panel stays at round-off — no μ can create a net viscous force where ω = 0."
14. `nb.note` — **N03 [C]** "**Separation breaks the limit.** As Re grows, an attached boundary layer thins toward zero
    and real flow approaches ideal flow. But where the layer separates (behind a cylinder, off a sharp corner) a wake of
    spinning fluid forms whose size does not shrink with Re: the limit μ → 0 of real flow is then *not* the flow with
    μ = 0. Ch. 9 explains separation; C06 shows its pressure signature on a cylinder."
15. `nb.note` — **N04 [C]** "**In and out.** Excluded: fluids of varying density, high subsonic and supersonic speeds,
    boundary layers, wakes, flows inside pipes, any region where elements spin. Included: flight (Ch. 14), water waves
    (Ch. 7), flow round vehicles and structures."
16. `nb.md` — "> ⚠️ **Common confusion:** 'ideal flow = inviscid fluid'. The fluid can be as viscous as honey; what
    matters is ω = 0. D01 shows the viscous *stress* survives (step 6 of the tiny example); only its net force vanishes.
    What ideal flow gives up is the *no-slip* condition."
17. `nb.md` — **What would change if…** "…the fluid were stratified (density varying with height)? Then the baroclinic
    torque ∇ρ × ∇p/ρ² (Ch. 5 §5.2) creates vorticity in the interior and Kelvin's theorem fails: the outer flow is no
    longer irrotational. That is the everyday case in the atmosphere and ocean — which is why Ch. 13 needs potential
    *vorticity*, not a potential."

### A.2 §6.2 Two-Dimensional Stream Function and Velocity Potential — R02 R03 R04 R06 R07 R08 R09, C02 (+N05–N15 N18 N19 · D02 D03 D04), R05, C03 (+N16 N17 · D05)
1. `nb.section("6.2", "Two-Dimensional Stream Function and Velocity Potential", intro="**What is this section about?**
   In a plane flow two scalar functions carry the whole velocity field: the stream function ψ (its contours are
   streamlines) and the velocity potential φ (its gradient is the velocity). Irrotational flow turns both into solutions
   of Laplace's equation, which is *linear* — and that is the key to the whole chapter.")`
2. `nb.recap("R02", "Plane continuity", "For constant density in two dimensions, $\\partial u/\\partial x+\\partial
   v/\\partial y=0$ *(Eq. 6.2)* — the plane form of $\\nabla\\cdot\\mathbf u=0$ *(Eq. 4.10)*.", where="Ch. 4 §4.2")`
3. `nb.recap("R03", "The stream function", "Define $u\\equiv\\partial\\psi/\\partial y$, $v\\equiv-\\partial\\psi/\\partial x$
   *(Eq. 6.3)*; then (6.2) holds for any smooth ψ (mixed partials commute). Along a curve ψ = const, $0=d\\psi=-v\\,dx+u\\,dy$, so
   $(dy/dx)_\\psi=v/u$ — the slope of a streamline. The difference of ψ between two streamlines is the volume flow between
   them per unit depth [m²/s].", where="Ch. 4 §4.3")` + `nb.code`: `from fluidpy.core import streamfunction as sf` ·
   `print(sf.velocity_from_streamfunction_2d(lambda x, y: 2.0*y, 0.3, 0.7))` (ψ = Uy with U = 2 m/s). *expect:* `(2.0,
   0.0)` to 1e-9.
4. `nb.recap("R04", "Plane irrotationality", "A 2-D flow has no spin when $\\partial v/\\partial x-\\partial u/\\partial y=0$
   *(Eq. 6.9)*, the z-component of ∇×u = 0 *(Eq. 3.23)*; then (and, in a simply connected region, only then) a velocity
   potential with $\\mathbf u=\\nabla\\phi$ *(Eq. 3.17)* exists.", where="Ch. 3 §3.4")`
5. `nb.recap("R06", "Polar continuity", "$\\frac1r\\frac{\\partial}{\\partial r}(ru_r)+\\frac1r\\frac{\\partial u_\\theta}{\\partial
   \\theta}=0$ *(Eq. 6.19)* — Appendix B's divergence in polar coordinates.", where="Ch. 4 (App. B operators)")`
6. `nb.recap("R07", "Polar irrotationality", "$\\frac1r\\frac{\\partial}{\\partial r}(ru_\\theta)-\\frac1r\\frac{\\partial
   u_r}{\\partial\\theta}=0$ *(Eq. 6.20)*.", where="Ch. 3 §3.4")`
7. `nb.recap("R08", "Polar Laplacian of ψ", "$\\nabla^2\\psi=\\frac1r\\frac{\\partial}{\\partial r}\\big(r\\frac{\\partial
   \\psi}{\\partial r}\\big)+\\frac1{r^2}\\frac{\\partial^2\\psi}{\\partial\\theta^2}=0$ *(Eq. 6.23a)* — used in D03.",
   where="Ch. 4 (App. B operators)")`
8. `nb.recap("R09", "Polar Laplacian of φ", "The same operator on φ, $\\frac1r\\frac{\\partial}{\\partial
   r}\\big(r\\frac{\\partial\\phi}{\\partial r}\\big)+\\frac1{r^2}\\frac{\\partial^2\\phi}{\\partial\\theta^2}=0$ *(Eq. 6.23b)*.",
   where="Ch. 4 (App. B operators)")` + `nb.code` (checks R06–R09 at once): `r, th = sp.symbols("r theta", positive=True)`
   · `lap = lambda f: sp.diff(r*sp.diff(f, r), r)/r + sp.diff(f, th, 2)/r**2` (the polar Laplacian of (6.23a, b)) ·
   `print([sp.simplify(lap(f)) for f in (r*sp.sin(th), sp.log(r), sp.cos(th)/r, r**2*sp.sin(2*th))])`. *expect:* `[0,
   0, 0, 0]` — "the stream, the vortex/source, the doublet and the corner are all harmonic (for r > 0)".

**C02 — ψ and φ: ω_z = −∇²ψ, two Laplace problems, vortices and sources as δ sources (6.4)–(6.15)**
9. `nb.core("C02", "Vorticity is minus the Laplacian of ψ: $\\omega_z=-\\nabla^2\\psi$ (6.4) — and vortices and sources are
   point sources of Laplace's equation", question="Every irrotational plane flow obeys one equation. Which — and what makes
   a whirlpool or a spring different from the flow around them?")`
10. `nb.md` — **The problem in plain words:** "Look down at a bath draining (a whirlpool) and at a garden sprinkler lying
    flat (a spring). Around both the water moves without spinning; yet one stirs the water round and the other pushes
    it out. We want one equation for all such flows, and a way to say exactly how strong the whirl or the spring is."
11. `nb.md` — **The idea** (table): | | stream function ψ | potential φ | — | defined by | u = ∂ψ/∂y, v = −∂ψ/∂x |
    u = ∂φ/∂x, v = ∂φ/∂y | — | automatic | continuity | irrotationality | — | the other law gives | ∇²ψ = −ω_z | ∇²φ = q
    (sources) | — | point singularity | vortex: −Γδ | source: mδ |. "**Away from a few special points both ψ and φ are
    harmonic: ∇² = 0. The special points carry all the information.**"
12. `nb.md` — 🔁 reminders: partial derivative and Laplacian (P25; Ch. 2 §2.9), Dirac delta as a concentrated total
    (Ch. 2 gloss, Ch. 5 P139), `np.meshgrid` (P76), contour plots (P78).
13. `nb.derivation("D02", …)` — Part F D02 (5 steps), ref "6.4".
14. `nb.note` — **N05 [B]** "In irrotational flow (6.4) becomes Laplace's equation for ψ:" equation `\nabla^2\psi=0`,
    ref "6.5". "Solutions of Laplace's equation are called *harmonic*."
15. `nb.note` — **N06 [B]** "A point (ideal) vortex of strength Γ at (x′, y′) is all its vorticity squeezed into one
    point — a delta function:" equation `\nabla^2\psi=-\Gamma\,\delta(x-x')\,\delta(y-y')`, ref "6.6". "D03 below shows
    that $\psi=-\frac{\Gamma}{2\pi}\ln r$ does exactly this, although the book never writes it out."
16. `nb.primer("2-D divergence theorem and the Dirac delta in the plane", "In the plane, Gauss' theorem says the total
    of ∇²f over a region equals the outward flux of ∇f through its edge: $\\int_A\\nabla^2 f\\,dA=\\oint_C\\nabla f\\cdot\\mathbf
    n\\,ds$. A 2-D delta $\\delta(x)\\delta(y)$ is zero everywhere except at the origin, with total 1 over any region that
    contains it. So a function whose Laplacian is zero away from the origin but whose flux through every circle round it
    is the same number F has $\\nabla^2 f=F\\,\\delta(x)\\delta(y)$.", code="th = np.linspace(0, 2*np.pi, 400,
    endpoint=False)   # 400 points on a circle\nR = 3.0                                   # radius [m]\nflux = np.sum(1/R *
    R) * (2*np.pi/400)     # ∇ln r·n = 1/R, ds = R dθ\nprint(flux, 2*np.pi)                      # 6.283185… both: the flux
    of ∇ ln r is 2π for any R")` — *expect:* `6.283185307179586 6.283185307179586`.
17. `nb.derivation("D03", …)` — Part F D03 (9 steps), no book number (ref "6.6"; the title says "the book never writes
    it out").
18. `nb.note` — **N07 [B]** "Uniform flow in any direction:" equation `\psi=-Vx+Uy`, ref "6.7". "U = 2, V = 1 m/s: ψ(1, 1)
    = −1 + 2 = 1 m²/s. Its twin:" + **N13 [B]** equation `\phi=Ux+Vy`, ref "6.14".
19. `nb.note` — **N08 [B]** "The ideal vortex at (x′, y′):" equation
    `\psi=-\frac{\Gamma}{2\pi}\ln\sqrt{(x-x')^2+(y-y')^2}`, ref "6.8". "It spins the fluid round at $u_\theta=\Gamma/2\pi r$
    *(Eq. 5.2)* (Ch. 3 §3.5, Ch. 5 §5.1). > ⚠️ Γ is **counterclockwise positive** here and in all our code; the cylinder
    with circulation (C07) will use the book's clockwise Γ."
20. `nb.note` — **N14 [B]** "The point source (m > 0) or sink (m < 0):" equation
    `\phi=\frac{m}{2\pi}\ln\sqrt{(x-x')^2+(y-y')^2}`, ref "6.15". "m is the volume flow out per unit depth [m²/s]: with
    m = 2π m²/s, $u_r=m/2\pi r=1/r$ — 1 m/s at 1 m, 0.5 m/s at 2 m."
21. `nb.note` — **N10 [B]** "Continuity in terms of φ is Poisson's equation with a source density q [1/s]:" equation
    `\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=\nabla^2\phi=q(x,y)`, ref "6.11". "It is the φ-twin of
    (6.4) (ch05 P139 did Poisson in 3-D)." + **N11 [B]** equation `\nabla^2\phi=0`, ref "6.12" + **N12 [B]** equation
    `\nabla^2\phi=m\,\delta(x-x')\,\delta(y-y')`, ref "6.13". "D03 step 8 derived this for the source."
22. `nb.worked_example("a whirlpool and a spring with round numbers", "1. Vortex, Γ = 2π m²/s at the origin: $\\psi=-\\ln r$
    [m²/s]; $u_\\theta=-\\partial\\psi/\\partial r=1/r$ → 1 m/s at r = 1 m. 2. Flux of ∇ψ out of the circle r = 2 m:
    $\\partial\\psi/\\partial r=-1/2$ m/s, times the length 2π × 2 = 4π m → −2π m²/s = −Γ ✓ (6.6), for **any** radius (try
    r = 5: −1/5 × 10π = −2π). 3. Source, m = 2π m²/s: $u_r=1/r$; flow out of the circle r = 2: 0.5 × 4π = 2π m²/s = m ✓
    (6.13). 4. Away from the centres, $\\nabla^2\\ln r=\\frac1r\\frac{d}{dr}(r\\cdot\\frac1r)=0$ (6.23a).")`
23. `nb.code` — *code:* `els = {"uniform": pf.Uniform(2.0, 1.0), "source": pf.Source(2*np.pi), "vortex":
    pf.Vortex(2*np.pi)}` · `pts = (np.array([0.7, -1.2, 2.0]), np.array([0.4, 0.9, -1.5]))` · `for k, e in els.items():
    print(k, np.max(np.abs(pf.laplacian_residual(e.psi, *pts))), np.max(np.abs(pf.laplacian_residual(e.phi, *pts))))` ·
    `print(ch06.delta_flux_check("vortex", Gamma=2*np.pi))` · `print(ch06.delta_flux_check("source", kind="phi",
    m=2*np.pi))` · `print(ch06.delta_flux_check("vortex", center=(5.0, 0.0), radii=(1.0,), Gamma=2*np.pi))` (circle not
    enclosing the vortex). *expect:* residuals ≲ 1e-8 (4th-order stencil); `[-6.2832 -6.2832 -6.2832 -6.2832 -6.2832]`;
    `[6.2832 … 6.2832]`; `[~0]` (≤ 1e-12). *explain:* 1. three elements as objects with ψ and φ; 2. the Laplacian of
    each at three ordinary points is zero (harmonic); 3. the flux of ∇ψ of the vortex through circles of radius 0.01 …
    100 m is −Γ every time — the delta of (6.6); 4. the source's flux is +m (6.13); 5. a circle that misses the
    centre gets 0.
24. `nb.code` — **Rankine vortex, ω = −∇²ψ:** `x = np.linspace(-0.3, 0.3, 121)`; `X, Y = np.meshgrid(x, x)` ·
    `om = ch06.vorticity_from_psi("rankine", X, Y, Gamma=1.0, a=0.1)` · `print(om[60, 60], om[60, 100])` (centre, and
    r = 0.2 m). *expect:* `31.831 0.0` (Γ/πa² = 1/(π·0.01) = 31.83 s⁻¹ inside, 0 outside; ≤ 1e-6). *explain:* the ψ of
    a vortex with a finite core (Ch. 3 §3.5); its −∇²ψ is the uniform core vorticity inside and zero outside — Laplace
    outside, Poisson inside.
25. `nb.check_agree` — **from scratch (curation §7):** five-point −∇²ψ by hand on the Rankine ψ: `psi =
    ch06.rankine_vortex_psi(X, Y, Gamma=1.0, a=0.1)`; `h = x[1] - x[0]`; `lap = (psi[1:-1, 2:] + psi[1:-1, :-2] +
    psi[2:, 1:-1] + psi[:-2, 1:-1] - 4*psi[1:-1, 1:-1])/h**2`; `mask = np.hypot(X, Y)[1:-1, 1:-1]`; `ok =
    np.abs(mask - 0.1) > 2*h` (away from the core edge, where ψ'' jumps); `assert np.allclose(-lap[ok], om[1:-1,
    1:-1][ok], atol=1e-6)` · loop sum of ∇ψ·n over a circle: `th = np.linspace(0, 2*np.pi, 512, endpoint=False)`; `R =
    0.7`; `xs, ys = R*np.cos(th), R*np.sin(th)`; `u, v = pf.Vortex(1.0).velocity(xs, ys)`; `dpsidn = -(u*(-np.sin(th)) +
    v*np.cos(th))` (∂ψ/∂n = −u_θ, from (6.22)); `assert np.isclose(np.sum(dpsidn)*R*2*np.pi/512,
    ch06.delta_flux_check("vortex", radii=(0.7,), Gamma=1.0)[0])`. Markdown: "Our own stencil and our own loop sum
    agree with the library."
26. `nb.md` — 🔁 reminder before D04: total differential dψ = ψ_x dx + ψ_y dy (Ch. 4 §4.3), dot product and
    perpendicular vectors (Ch. 2 §2.2), level sets (P75).
27. `nb.derivation("D04", …)` — Part F D04 (5 steps), ref "6.10".
28. `nb.note` — **N09 [B]** "The potential (in any irrotational flow; single-valued if the region is simply connected,
    Ch. 3 §3.4):" equation `u\equiv\partial\phi/\partial x,\ v\equiv\partial\phi/\partial y;\qquad (dy/dx)_{\phi=\text{const}}=-u/v`,
    ref "6.10". "D04 showed the lines of constant φ cross the streamlines at right angles, and ∇ψ is ∇φ turned by +90°."
    + `nb.code`: `print(ch06.orthogonality_check(pf.cylinder(1.0, 1.0, Gamma_cw=2.0), (np.array([1.5, -2, 0.3]),
    np.array([0.4, 1.1, -1.8]))))` → *expect:* ≲ 1e-10.
29. `nb.note` — **N18 [B], N19 [B]** "In polar coordinates (a chain-rule exercise the book gives 'for quick reference'):"
    equation `u_r=\frac{\partial\phi}{\partial r}=\frac1r\frac{\partial\psi}{\partial\theta},\qquad u_\theta=\frac1r\frac{\partial\phi}{\partial\theta}=-\frac{\partial\psi}{\partial r}`,
    ref "6.21, 6.22" (the builder writes the two numbers after the two formulas). "Number: the vortex with Γ = 2π has
    ψ = −ln r, so $u_\theta=-\partial\psi/\partial r=1/r$." + `nb.code`: `print(pf.polar_velocity_sym(-sp.log(r), r, th,
    kind="psi"))` → `(0, 1/r)`.
30. `nb.figure` — **flow nets** (three panels 9 × 3 in, `flow_net`, ψ solid, φ dashed): uniform (U = 1, V = 0.5),
    source (m = 2π), vortex (Γ = 2π); the source's φ circles and the vortex's ψ circles coincide in shape — "a source
    is a vortex with ψ and φ swapped". *see:* "in each panel two families of curves crossing at right angles." *read:*
    "solid lines are streamlines (equal Δψ, so crowded lines mean fast flow); dashed are equipotentials (equal Δφ);
    D04's right angles hold everywhere except at the centres." *change:* "…the vortex strength doubled: the same circles,
    twice as crowded — the velocity doubles (spacing = speed, the Ch. 4 `stream_function_spacing` explainer)."
31. `nb.figure` — **flux vs radius** (log-x, radii 0.01–100 m): ∮∇ψ·n ds for the vortex (amber, flat at −Γ), ∮∇φ·n ds
    for the source (teal, flat at m), and for a circle centred 5 m away (grey, 0 until the radius reaches 5 m, then −Γ).
    *see:* "flat lines." *read:* "the flux does not care about the radius: the whole Laplacian sits at one point — the
    delta of (6.6) and (6.13)." *change:* "…a Rankine core of radius 0.1 m: the amber line would rise from 0 at r = 0
    to −Γ at r = 0.1 m and stay flat — the delta is smeared into a disc."
32. `nb.figure` — **ω = −∇²ψ heatmap** of the Rankine vortex (Γ = 1 m²/s, a = 0.1 m) with streamlines. *see / read /
    change* — "uniform 31.8 s⁻¹ inside, zero outside; outside, ψ is harmonic and the flow is irrotational although it
    goes round; shrinking a with Γ fixed turns the core into the point vortex of (6.6)."
33. `nb.note` — **N15 [C]** "**Two descriptions, one flow.** ψ and φ each describe the same plane ideal flow and will fuse
    into one complex potential $w=\phi+i\psi$ in C09. ψ also works for rotational plane flow (C02's Rankine core); φ works
    for unsteady and 3-D flow (C13, C15)."
34. `nb.md` — **What would change if…** "…the vorticity were spread out (a smooth vortex)? Then $\nabla^2\psi=-\omega_z$ is a
    Poisson equation with a smooth right-hand side, solved by the Green's-function sum of Ch. 5 §5.5 or on a grid
    (`ls.solve_poisson`, C12's method with a source term). That is exactly how ocean and atmosphere models recover the
    stream function from the vorticity (Ch. 13's potential-vorticity inversion)."

35. `nb.recap("R05", "Bernoulli everywhere in irrotational flow", "Steady, constant density, irrotational: the Bernoulli
    constant is the same on every streamline, $p+\\tfrac12\\rho\\lvert\\nabla\\phi\\rvert^2=\\text{const}$ *(Eq. 6.18, from
    4.72)*; since $\\lvert\\nabla\\psi\\rvert=\\lvert\\nabla\\phi\\rvert$ (D04) the same holds with ψ. Unsteady flow adds ρ∂φ/∂t:
    $\\rho\\,\\partial\\phi/\\partial t+p+\\tfrac12\\rho\\lvert\\nabla\\phi\\rvert^2=\\text{const}$ *(Eq. 4.75)*. So once ψ or φ is known,
    the pressure is free — the nonlinear momentum equation is solved by one algebraic line.", where="Ch. 4 §4.9")` +
    `nb.code`: `from fluidpy.core import bernoulli as be` · `f = pf.half_body(1.0, 2*np.pi)` · `X = np.array([-3.0, 0.5, 2.0,
    -1.0])`; `Y = np.array([2.0, 3.5, -4.0, 1.5])` (points on different streamlines) · `p = f.pressure(X, Y, rho=1000.0,
    p_inf=0.0, U=1.0)` · `print(be.bernoulli_function(f.speed(X, Y), p, 0.0, rho=1000.0, g=0.0))`. *expect:* four equal
    values 0.5 (J/kg) — "the same Bernoulli constant across streamlines".

**C03 — Superposition and the no-through-flow condition: any streamline can be a wall**
36. `nb.core("C03", "Superposition: $\\nabla^2(\\phi_1+\\phi_2)=\\nabla^2\\phi_1+\\nabla^2\\phi_2=0$ — and any streamline can be a
    wall (6.16)", question="If every element solves Laplace's equation, where does the body come from — and who imposes
    its boundary condition?")`
37. `nb.md` — **The problem in plain words:** "Dye streaks in a river bend round a stone you cannot see. Could we produce
    that pattern without ever drawing the stone — by adding flows we already know and then *declaring* one of the
    streamlines to be a wall?"
38. `nb.md` — **The idea** (ASCII):
    ```
    uniform stream  +  source  =  a flow with a stagnation point S
                                    │
         the streamline through S ──┴──► call it a wall: no flow crosses a streamline,
                                         so the no-through-flow condition (6.16) already holds there
    ```
    "**Laplace's equation is linear, so solutions add; the boundary condition is met by choosing which streamline is
    the body.**"
39. `nb.md` — 🔁 reminders: directional derivative ∂φ/∂n = ∇φ·n (P75), unit normal and tangent of a curve (P92) +
    gloss "**points as complex numbers in code:** fluidpy stores a point (x, y) as the complex number x + iy, written
    `-1+0j` for (−1, 0), and returns stagnation points the same way; C09 explains complex numbers properly (P153) — here
    read x + iy simply as the pair of coordinates (i² = −1 was met in Ch. 1, P45)".
40. `nb.derivation("D05", …)` — Part F D05 (5 steps), ref "6.16".
41. `nb.note` — **N16 [B]** "No flow through a solid surface, in its general and its stationary form:" equation
    `\mathbf n\cdot\mathbf U_s=(\mathbf n\cdot\mathbf u)_{\text{surface}};\qquad \partial\phi/\partial n=0\ \text{or}\ \partial\psi/\partial s=0\ \text{on a stationary surface}`,
    ref "6.16". "A stationary wall is a streamline, and any streamline may be replaced by a wall."
42. `nb.note` — **N17 [B]** "Far away the flow must return to the given stream (U along x; U = 0 is still fluid):"
    equation `\partial\phi/\partial x=U,\quad\text{or}\quad\partial\psi/\partial y=U`, ref "6.17". "How fast? For a closed
    body the disturbance dies like 1/R² (a doublet, C04); for a half-body with a net source like 1/R." + `nb.code`: `Rs =
    np.array([10.0, 20.0, 40.0, 80.0])` · `e_cyl = [pf.far_field_check(pf.cylinder(1.0, 1.0), R) for R in Rs]` · `e_half =
    [pf.far_field_check(pf.half_body(1.0, 2*np.pi), R) for R in Rs]` · `from tools.convergence import observed_order`
    (repo `tools/` on the path) · `print(observed_order(Rs, e_cyl), observed_order(Rs, e_half))`. *expect:* `-2.00 -1.00`
    (± 0.05): "the log–log slopes of the far-field error (P13)".
43. `nb.worked_example("a stagnation point from two elements", "Stream U = 1 m/s plus a source m = 2π m²/s at (−1, 0) m.
    At the point (−2, 0): the source is 1 m away, straight downstream of the point, so it pushes the fluid there
    **upstream** at $u=\\frac{m}{2\\pi}\\frac{x-x'}{(x-x')^2+(y-y')^2}=\\frac{2\\pi}{2\\pi}\\cdot\\frac{-1}{1}=-1$ m/s. 1. Sum:
    $u=1+(-1)=0$. 2. On the axis y = 0 the source adds no v, so $v=0$. 3. The fluid stops: a stagnation point, 1 m in
    front of the source — the nose of a body that nobody drew (C05 finds its whole shape).")`
44. `nb.code` — *code:* `f = pf.Flow([pf.Uniform(1.0), pf.Source(2*np.pi, -1+0j)])` · `print(f.velocity(-2.0, 0.0))` ·
    `print(f.stagnation_points())` · `psi_s = f.psi(-2.0, 0.0)`; `print(psi_s)` · `th = np.linspace(0.05, 2*np.pi-0.05,
    400)` · `xb, yb = ch06.half_body_shape(1.0, 2*np.pi, th)` (the dividing streamline of a source at the origin) ·
    `body = (xb - 1.0) + 1j*yb` (shifted to the source at x = −1) · `print(pf.normal_velocity_on(f, body))`. *expect:* `(0.0, 0.0)` ·
    `[-2.+0.j]` · `3.14159` (= m/2) · ≲ 1e-10. *explain:* 1. a flow is a list of elements; 2. its velocity is the sum of
    theirs; 3. Newton's method on dw/dz finds the stagnation point; 4. the streamline through it has ψ = m/2; 5. no flow
    crosses that curve, so it can be a wall.
45. `nb.check_agree` — **from scratch (curation §7):** the two velocities added by hand at a point: `x, y = 0.5, 1.5`;
    `u_hand = 1.0 + (2*np.pi/(2*np.pi))*(x+1)/((x+1)**2 + y**2)`; `v_hand = (2*np.pi/(2*np.pi))*y/((x+1)**2 + y**2)`;
    `assert np.allclose((u_hand, v_hand), f.velocity(x, y))`. Markdown: "The library does nothing more than add."
46. `nb.figure` — **stream + source = sum** (three panels, 9 × 3 in): uniform stream (blue lines), source (teal rays),
    their sum with the stagnation point S (black dot) and the dividing streamline ψ = m/2 drawn bold black (the body),
    the region inside shaded light grey "fluid from the source". *see:* "straight lines, rays, and their sum: a
    rounded nose appears." *read:* "no element knows about a body; the bold line is simply the streamline through S —
    by D05 it is a wall." *change:* "…a sink of equal strength added 3 m downstream: the bold curve closes into an oval
    (the Rankine oval — try it in E1 after C05)."
47. `nb.md` — **What would change if…** "…we added a vortex instead of a source? A vortex creates no stagnation point
    in a uniform stream by itself far from its centre, but adds circulation: combined with a closed body it will
    produce lift (C07). The rule stays the same: add, then find the streamline that can be a wall."

### A.3 §6.3 Construction of Elementary Flows in Two Dimensions — R10 R11, C04 (+N20–N26 · D06), C05 (+N27–N29 · D07 D08 · E1), R12 R13 R14, C06 (+N30 N31 · D09), C07 (+N32–N38 · D10 D11 · E2), R15 R16 R17, C08 (+N39 N40 · D12 D13 · E3)
1. `nb.section("6.3", "Construction of Elementary Flows in Two Dimensions", intro="**What is this section about?** A
   kit of elementary flows — stream, source, vortex, corner, doublet — and what you get by adding them: a blunt
   half-body, a cylinder with zero drag (d'Alembert's paradox), a spinning cylinder with lift, and walls made by mirror
   images. Every body here is a streamline of a sum.")`
2. `nb.recap("R10", "The stagnation (corner) flow", "The simplest quadratic stream function, $\\psi=2Axy$ *(Eq. 6.24)*,
   gives u = 2Ax, v = −2Ay: fluid comes down the y-axis, stops at the origin and leaves along the x-axis; streamlines are
   the hyperbolae xy = const. It is Ch. 3's irrotational strain and the Ch. 4 'stagnation' preset. ⚠️ In Ch. 13 the same
   field is the deformation flow that sharpens weather fronts.", where="Ch. 3 §3.4")` + `nb.code`: `c = pf.Corner(1.0,
   n=2)` · `print(c.velocity(1.0, 0.5), c.psi(1.0, 0.5))` → *expect:* `(2.0, -1.0) 1.0` (u = 2Ax, v = −2Ay, ψ = 2Axy).
3. `nb.recap("R11", "The vortex velocity", "Differentiating (6.8) $\\psi=-\\frac{\\Gamma}{2\\pi}\\ln r$ gives
   $u=-\\frac{\\Gamma}{2\\pi}\\frac{y}{x^2+y^2}$, $v=\\frac{\\Gamma}{2\\pi}\\frac{x}{x^2+y^2}$, i.e. $u_\\theta=\\Gamma/2\\pi r$ *(Eq.
   5.2)*: circles round the centre, irrotational everywhere else.", where="Ch. 5 §5.1")` + `nb.code`: `from fluidpy.core
   import biot_savart as bs` · `print(pf.Vortex(1.0).velocity(0.3, 0.4), bs.point_vortex_velocity(np.array([0.3, 0.4]),
   np.array([[0.0], [0.0]]), np.array([1.0])))` → *expect:* the same (−0.25465, 0.19099) twice.

**C04 — The element kit and the doublet as a source–sink limit (6.29)**
4. `nb.core("C04", "The doublet: a source and a sink pushed together, $\\phi=-\\frac{\\mathbf d\\cdot\\mathbf x}{2\\pi
   r^2}=\\frac{\\lvert\\mathbf d\\rvert}{2\\pi}\\frac{\\cos\\theta}{r}$ (6.29)", question="What is left when a source and a sink
   merge — and why is that the piece that makes closed bodies?")`
5. `nb.md` — **The problem in plain words:** "Hold two garden hoses nose to nose, one blowing water out, one sucking the
   same amount in. From far away the water just loops from one to the other. Push them closer while turning both up so
   the loops do not vanish: the limit is a new element, the doublet — the only one that disturbs the flow without adding
   or removing water or circulation. It is what turns a stream into flow round a cylinder, a sphere, an airship."
6. `nb.md` — **The idea** (ASCII):
   ```
   source +m at (−ε,0)   sink −m at (+ε,0)         ε → 0 with 2mε = |d| fixed
        ●───────────────────○             ──►      ⊙  doublet: streamlines are circles tangent to the x-axis
   dipole vector d = Σ x_i m_i = −2mε e_x  (points from the sink to the source)
   ```
7. `nb.md` — 🔁 reminders: Taylor series ln(1 + s) ≈ s (P26), exponent and log rules (P43, P36; gloss: ln√s = ½ ln s,
   ln(ab) = ln a + ln b, d(ln r)/dx = x/r²), orders of smallness (P68).
8. `nb.note` — **N20 [B]** "**Harmonic polynomials.** A constant (no flow), the linear stream (6.7)/(6.14), and two
   quadratic families solve Laplace; in general Re zⁿ and Im zⁿ (C09). Higher powers make sharper corners, fractional
   powers wider ones." + `nb.code`: `for d in (2, 3): print(d, ch06.harmonic_polynomials(d)[-1])` → *expect:* `2 (x**2 -
   y**2, 2*x*y)` · `3 (x**3 - 3*x*y**2, 3*x**2*y - y**3)`. *explain:* "sympy splits zⁿ = (x + iy)ⁿ into real and imaginary
   parts; each is harmonic (a D14 preview)."
9. `nb.note` — **N21 [B], N22 [B], N23 [B]** "**The corner family.** Swapping which function is which, or turning the axes
   by 45°, gives the relatives of (6.24):" equation `\phi=2Axy\ \ (6.25),\qquad \psi=A(x^2-y^2)\ \ (6.26),\qquad \phi=A(x^2-y^2)\ \ (6.27)`.
   "(6.25) is the (6.24) flow turned by 45° (x² − y² in axes turned by 45° becomes 2x′y′); (6.27) *is* (6.24): both are
   the real and imaginary parts of Az²." + `nb.code`: `x, y = 0.7, 0.3` · `print(pf.Corner(1.0, n=2).phi(x, y), x**2 -
   y**2)` ((6.27) is φ of (6.24)) · `c45 = pf.Corner(1.0, n=2, rotate=np.pi/4)` (w = A(ze^{−iπ/4})² = −iAz²) · `print(
   c45.phi(x, y), 2*x*y, c45.psi(x, y), -(x**2 - y**2))` (turned by 45°: φ = 2Axy is (6.25), ψ = −A(x² − y²) is (6.26)
   with A → −A). *expect:* `0.4 0.4` · `0.42 0.42 -0.4 -0.4`.
10. `nb.note` — **N24 [B]** "**Source velocities.** Differentiating the source potential (6.15)
    $\phi=\frac{m}{2\pi}\ln\sqrt{x^2+y^2}$ gives $u=\frac{m}{2\pi}\frac{x}{x^2+y^2}$, $v=\frac{m}{2\pi}\frac{y}{x^2+y^2}$, so
    $u_r=m/2\pi r$, $u_\theta=0$. ⚠️ The book says 'differentiation of (6.8)' here; it means (6.15). ∇·u = 0 except at
    r = 0, where all m m²/s come out (D03)."
11. `nb.note` — **N25 [B]** "**The pair.** A source +m at (−ε, 0) and a sink −m at (+ε, 0):" equation
    `\phi=\frac{m}{2\pi}\ln\sqrt{(x+\varepsilon)^2+y^2}-\frac{m}{2\pi}\ln\sqrt{(x-\varepsilon)^2+y^2},\qquad\mathbf d=\sum_i\mathbf x_im_i=-2m\varepsilon\,\mathbf e_x`,
    ref "6.28". "The dipole vector points from the sink to the source."
12. `nb.primer("a limit with a product held fixed", "Sometimes two quantities go to extremes together: here the gap ε → 0
    while the strength m → ∞, with the product 2mε held at a fixed value |d|. Either limit alone gives nothing useful
    (ε → 0 at fixed m cancels the pair; m → ∞ at fixed ε blows up). Held together, the leading term survives and the
    rest shrinks like ε².", code="d = 2.0                                    # fixed dipole strength |d| = 2mε
    [m³/s]\nfor eps in (0.2, 0.02, 0.002):               # shrink the gap …\n    m = d/(2*eps)                          # … and
    raise the strength to keep 2mε = d\n    phi = m/(2*np.pi)*np.log((1+eps)/(1-eps))  # pair potential at (1, 0) m
    [m²/s]\n    print(eps, phi, d/(2*np.pi))           # → 0.31831 (the doublet value), error ∝ ε²")` — *expect:* 0.322659,
    0.318352, 0.3183103 against 0.3183099.
13. `nb.derivation("D06", …)` — Part F D06 (9 steps), ref "6.29".
14. `nb.worked_example("a strong pair, 4 cm apart", "m = 50 m²/s, ε = 0.02 m, so $\\lvert\\mathbf d\\rvert=2m\\varepsilon=2$
    m³/s. At (1, 0) m: 1. pair, (6.28): $\\phi=\\frac{50}{2\\pi}[\\ln1.02-\\ln0.98]=7.9577\\times0.0400053=0.318352$ m²/s. 2.
    doublet, (6.29): $\\phi=\\frac{\\lvert\\mathbf d\\rvert}{2\\pi}\\frac{\\cos0}{1}=0.318310$ m²/s. 3. relative difference
    1.33×10⁻⁴ ≈ ε²/3 — the dropped terms of D06 step 9.")`
15. `nb.code` — *code:* `eps = np.array([0.2, 0.1, 0.05, 0.025, 0.0125])` · `err = ch06.doublet_limit_error(eps,
    d=2.0)` · `print(err, observed_order(eps, err))` · `D = pf.Doublet((-2.0, 0.0))` (dipole −2 e_x, as the pair) ·
    `print(D.phi(1.0, 0.0), D.psi(0.0, 1.0))` · `print(pf.Doublet.from_book_scalar(2.0).phi(1.0, 0.0))`. *expect:* errors
    falling ×4 per halving, order 2.00 ± 0.05 · `0.31831 -0.31831` (ψ = −|d|y/2πr² for this d) · `0.31831` (the book's
    scalar d = 2 is the same dipole). *explain:* 1. the pair's φ against the doublet's at three points, for five gaps at
    fixed 2mε; 2. the log–log slope is 2: the error ∝ ε²; 3. the doublet object with its vector d; 4. (6.49)'s scalar d
    means the dipole −d e_x (convention 4).
16. `nb.check_agree` — **from scratch (curation §7):** the pair formula at shrinking ε by hand vs `pf.Doublet`: `for e in
    eps: m = 2.0/(2*e); pair = m/(2*np.pi)*(np.log(np.hypot(0.6+e, 0.8)) - np.log(np.hypot(0.6-e, 0.8)))`;
    `assert np.isclose(pair, D.phi(0.6, 0.8), rtol=e**2)` (relative error below ε²). Markdown: "The pair tends to the
    doublet with an error that falls like ε², as D06 predicted."
17. `nb.note` — **N26 [C]** "**Where elements are singular.** (6.7), (6.14), (6.24)–(6.27) are harmonic everywhere; (6.8),
    (6.15), (6.29) are harmonic for r > 0 only — their centres are where the flow's 'information' sits (C02). Our
    residual table prints NaN there on purpose (`laplacian_residual` refuses points within 5h of a centre)."
18. `nb.figure` — **the element kit** (5 panels, 11 × 2.6 in, `flow_net`): uniform (blue), source (teal), vortex
    (amber), corner n = 2 (grey walls), doublet with d = −2 e_x (purple; streamline circles tangent to the x-axis at the
    origin). Title "Five elements, all harmonic". *see:* "five flow nets; the doublet's streamlines are circles through
    the origin." *read:* "every panel solves Laplace; only the singular centres differ — the doublet's is a source and a
    sink merged." *change:* "…the doublet turned to point along +y: its circles turn by 90° and touch the y-axis."
19. `nb.animation` — **the pair becomes a doublet** (`player="frames"`, 12 frames, FAST 8): ε from 0.5 to 0.02 m with
    2mε = 2 m³/s fixed; ψ contours from `ch06.doublet_limit_frames(eps_list)` on a 161 × 161 grid (cached), the source
    teal and sink rose dots moving together, the relative error printed in the title. + see/read/change as a markdown
    cell: "the loops between the hoses become the doublet's tangent circles; the title's error falls ×4 each time ε
    halves."
20. `nb.md` — **What would change if…** "…the pair were a vortex pair (+Γ at (0, ε), −Γ at (0, −ε)) squeezed with 2Γε
    fixed? You get a doublet again, turned by 90° — the same element (Exercise-style, see S02)."

**C05 — The half-body: a stream plus a source (6.31)**
21. `nb.core("C05", "The half-body $\\psi=Uy+\\frac{m}{2\\pi}\\tan^{-1}\\frac yx=Ur\\sin\\theta+\\frac{m}{2\\pi}\\theta$ (6.31)",
    question="A stream meets a source. Where does the fluid stop, what shape is the body, and how wide does it get?")`
22. `nb.md` — **The problem in plain words:** "River water parts round the rounded nose of a bridge pier; wind parts round
    the leading edge of a wing or a cliff top. The simplest model of such a nose is a uniform stream plus one source:
    the source's fluid fills a body, the stream's fluid flows round it. We want the stagnation point, the shape, the
    width — and the pressure along the surface."
23. `nb.md` — **The idea** (ASCII):
    ```
    stream U →   ──────►          ψ = m/2 (dividing streamline, bold)
                 S ●─ a ─● source m       ────────────────────  h → h_max = m/2U
    stagnation where the source's push m/2πa equals U      inside: all the source's fluid, carried off at speed U
    ```
24. `nb.md` — 🔁 reminders: `np.arctan2` (P70) + gloss "**choosing the angle's range** — ψ contains θ itself, which jumps
    by 2π across a cut; we take θ ∈ [0, 2π), so the jump lies along +x behind the source, *inside* the body, and ψ is
    smooth everywhere in the fluid"; `brentq` (P108); polar coordinates (P105).
25. `nb.note` — **N27 [B]** "The same flow in φ:" equation
    `\phi=Ux+\frac{m}{2\pi}\ln\sqrt{x^2+y^2}=Ur\cos\theta+\frac{m}{2\pi}\ln r`, ref "6.30".
26. `nb.derivation("D07", …)` — Part F D07 (10 steps), ref "6.31".
27. `nb.worked_example("a pier nose with round numbers", "U = 1 m/s, m = 2π m²/s. 1. $a=m/2\\pi U=1$ m: the stagnation
    point is at (−1, 0) m. 2. The dividing streamline has ψ = m/2 = π m²/s. 3. Right above the source (θ = 90°):
    $h=m(\\pi-\\theta)/2\\pi U=2\\pi\\cdot(\\pi/2)/(2\\pi)=\\pi/2\\approx1.571$ m. 4. Far downstream $h_{\\max}=m/2U=\\pi\\approx3.142$ m. 5.
    Mass check: 2 × 3.142 m × 1 m/s = 6.283 m²/s = m ✓.")`
28. `nb.code` — *code:* `f = pf.half_body(1.0, 2*np.pi)` · `print(ch06.half_body_numbers(1.0, 2*np.pi))` · `print(
    f.stagnation_points())` · `th = np.array([np.pi/2, np.pi/4, 0.05, 3*np.pi/2])` · `xb, yb = ch06.half_body_shape(1.0,
    2*np.pi, th)` · `print(yb, f.psi(xb, yb))`. *expect:* `{'a': 1.0, 'h_max': 3.14159, 'psi_body': 3.14159}` · `[-1.+0.j]`
    · y = 1.5708, 2.3562, 3.0916, −1.5708 and ψ = 3.14159 at all four (both halves — the θ ∈ [0, 2π) choice).
    *explain:* 1. the half-body as a stream + a source; 2. its three numbers; 3. Newton finds the stagnation point; 4.
    points on the body from D07 step 6 lie on ψ = m/2, above and below.
29. `nb.note` — **N28 [B]** "**Pressure coefficient.** Bernoulli (6.18) with the constant set far upstream, p∞ + ½ρU²:"
    equation `C_p=\frac{p-p_\infty}{\frac12\rho U^2}=1-\frac{\lvert\mathbf u\rvert^2}{U^2}`, ref "6.32". "It is the Euler
    number of (4.106) (Ch. 4 §4.11). C_p = 1 at every stagnation point, 0 where the speed equals U, negative (suction)
    where the flow is faster than the stream."
30. `nb.derivation("D08", …)` — Part F D08 (8 steps), no book equation (ref "6.32"; "the book only plots this").
31. `nb.note` — **N29 [B]** "**Along the body** C_p starts at +1 at the nose, crosses zero where tan θ = −2(π − θ) — our
    root 113.2° from +x (the plotted curve in the book agrees; its value stays in the private file) — has its lowest value
    −0.587 near θ = 63°, and returns to 0 far downstream. The net pressure force on the whole (infinitely long) body is
    zero: the tail contributes nothing because C_p → 0 there." + `nb.code`: `print(np.degrees(
    ch06.half_body_cp_zero_angle()))` · `print(ch06.half_body_surface_cp(np.radians([179.0, 120.0, 90.0, 63.0, 10.0])))` ·
    `for xe in (10.0, 100.0, 1000.0): print(xe, ch06.half_body_net_force(1.0, 2*np.pi, xe, rho=1.0))`. *expect:*
    `113.218` · `[0.9997 0.1431 -0.4053 -0.5865 -0.1187]` · D shrinking toward 0 (roughly like 1/x_end), L = 0.
32. `nb.check_agree` — **from scratch (curation §7):** `from scipy.optimize import brentq` · `f_root = lambda t: np.tan(t)
    + 2*(np.pi - t)` (C_p = 0 on the body, D08 step 7) · `t0 = brentq(f_root, np.pi/2 + 1e-9, np.pi - 1e-9)` (tan
    changes sign between 90° and 180°) · `assert np.isclose(t0, ch06.half_body_cp_zero_angle())`. Markdown: "Our
    bracket and our function give the library's angle."
33. `nb.figure` — two panels (9 × 3.4 in): left, ψ contours of the half-body (equal Δψ, blue) with the dividing
    streamline ψ = m/2 bold black, S marked, a and h_max annotated, the source's fluid shaded; right, C_p on the body vs
    θ from 180° (nose) to 5°, the zero at 113.2° marked, the minimum marked, C_p = 1 at the nose. Title "Blunt nose:
    pressure rises at the front, suction on the shoulders". *see:* "a rounded nose opening to a width 2h_max; a C_p curve
    falling from +1 through zero to a suction dip and back to zero." *read:* "high pressure where the stream is stopped
    (nose), low where it is fastest (shoulder), recovery downstream; spacing of streamlines = speed." *change:* "…m
    doubled: every length (a, h, h_max) doubles, the C_p curve against θ is unchanged — only m/U sets the size."
34. `nb.plotly` — `slider_figure` over m/U from 0.5 to 5 m (20 steps, FAST 10): traces "body ψ = m/2" (black),
    "streamlines" (4 blue ψ contours as polylines from `contourpy`), "stagnation point" (black marker); axes fixed x ∈ [−3,
    8], y ∈ [−6, 6] m; the title shows a = m/2πU and h_max = m/2U. Title "The half-body scales with m/U". + a markdown
    note "What you see / How to read it / What would change if" in the same pattern.
35. `nb.explainer("superposition_sandbox", heading="Where does the body come from?", why="A static figure shows a finished
    sum. Here you build it: toggle elements, slide their strengths and watch the stagnation points and the dividing
    streamline appear, move and close; click a point to see each element's velocity add up.", tries=["Preset
    'half-body': read a = 1 m and h_max = 3.14 m in the Explain tab, then halve U — both double.", "Add the sink
    (preset 'Rankine oval'): the body closes and the status turns green: Σm = 0.", "Preset 'source–sink → doublet':
    press ▶ and watch the pair collapse into the doublet's circles.", "Click a point in the flow: the term bars show the
    stream's and the source's velocities adding to the total."])`
36. `nb.md` — **What would change if…** "…a sink of the same strength sat 2 m downstream of the source? The source's
    fluid would be swallowed again, the dividing streamline would close behind the body, and the drag of this closed
    body would be exactly zero — the Rankine oval (Exercise 6.19 idea; `ch06.rankine_oval(1.0, 2*np.pi, 1.0)` gives a
    half-length √3 ≈ 1.732 m and half-width 1.307 m). Closing the body is the step from C05 to C06."

37. `nb.recap("R12", "The cylinder = stream + doublet", "A uniform stream U plus a doublet $\\mathbf d=-2\\pi Ua^2\\mathbf e_x$
    at the origin gives $\\phi=U(r+a^2/r)\\cos\\theta$, $\\psi=U(r-a^2/r)\\sin\\theta$ *(Eq. 6.33)*; ψ = 0 on r = a for every θ,
    so the circle is a streamline, and the body closes because the doublet adds no net source. You met this flow in
    Ch. 3 (Fig. 3.2); D09 builds it.", where="Ch. 3 §3.3")` + `nb.code`: `print(pf.cylinder(2.0, 0.5).velocity(0.3, 0.8),
    ch03.cylinder_flow(0.3, 0.8, U=2.0, a=0.5))` → *expect:* the same pair twice (parity).
38. `nb.recap("R13", "The cylinder's velocity", "$u_r=U(1-a^2/r^2)\\cos\\theta$, $u_\\theta=-U(1+a^2/r^2)\\sin\\theta$ *(Eq.
    6.34)*: on r = a the radial part vanishes and the fluid slides round at $\\lvert u_\\theta\\rvert=2U\\lvert\\sin\\theta\\rvert$.",
    where="Ch. 3 §3.3")`
39. `nb.recap("R14", "A cylinder moving through still fluid", "Seen from the fluid far away, the moving cylinder's flow is
    the cylinder flow minus the stream: $w=U(z+a^2/z)-Uz=Ua^2/z$ — the instantaneous streamlines of a doublet (Ch. 3's
    `galilean_frames_cylinder` explainer shows both frames; open it from the Ch. 3 page).", where="Ch. 3 §3.3")` +
    `nb.code`: `print(ch03.cylinder_flow(0.3, 0.8, U=2.0, a=0.5, frame="fluid"), pf.Doublet((-2*np.pi*2.0*0.25,
    0.0)).velocity(0.3, 0.8))` → *expect:* equal pairs (the doublet d = −2πUa² e_x).

**C06 — The circular cylinder and d'Alembert's paradox**
40. `nb.core("C06", "d'Alembert's paradox: the cylinder's surface pressure $C_p=1-4\\sin^2\\theta$ (6.35) is symmetric, so
    the drag is zero", question="Wind presses hard on the front of a chimney. Why does ideal flow predict no net push at
    all?")`
41. `nb.md` — **The problem in plain words:** "Stand behind a chimney on a windy day and you are sheltered; the chimney
    itself feels a strong push downwind. Yet the ideal-flow solution for a circular cylinder — a flow you have already
    drawn in Ch. 3 — predicts that the push is exactly zero. This is d'Alembert's paradox (1752), and understanding why
    the prediction fails is the reason Ch. 9 exists."
42. `nb.md` — **The idea** (ASCII):
    ```
          suction (−3q)                 pressure on the surface, q = ½ρU²:
     +q  ( ○ )  +q        ─►            front +q, shoulders −3q, back +q
          suction (−3q)                 the back pushes forward exactly as hard as the front pushes back
    ```
    "**The pattern is the same front and back, top and bottom — so every push is cancelled.**"
43. `nb.md` — 🔁 reminders: polar velocities (6.21)–(6.22) (C02), net pressure force −∮p n dA (P28).
44. `nb.derivation("D09", …)` — Part F D09 (6 steps), ref "6.35".
45. `nb.note` — **N30 [B]** "**Surface pressure:**" equation `C_p(r=a,\theta)=1-4\sin^2\theta`, ref "6.35". "Stagnation
    points (C_p = 1) at θ = 0 and π; the fluid is twice as fast as the stream at the shoulders θ = ±90°, where
    C_p = −3. Air (ρ = 1.2 kg/m³), U = 10 m/s: ½ρU² = 60 Pa, so +60 Pa front and back, −180 Pa at the shoulders."
46. `nb.primer("integrals of sines and cosines over a full period", "Over one full turn, 0 ≤ θ ≤ 2π, the positive and
    negative halves of sin θ cancel, and so do those of sin³θ, cos θ, sin θ cos θ and sin²θ cos θ: all integrate to zero.
    Squares do not cancel: $\\int_0^{2\\pi}\\sin^2\\theta\\,d\\theta=\\int_0^{2\\pi}\\cos^2\\theta\\,d\\theta=\\pi$ (their average is ½).
    Force integrals round a circle pick out exactly the terms that survive.", code="th = np.linspace(0, 2*np.pi, 2001)
    # one full turn\nfor f in (np.sin(th), np.sin(th)**3, np.sin(th)*np.cos(th), np.sin(th)**2):\n    print(round(
    np.trapezoid(f, th), 6))   # 0, 0, 0, then π = 3.141593")` — *expect:* `0.0 0.0 0.0 3.141593`.
47. `nb.worked_example("the chimney in numbers", "Air, U = 10 m/s, a = 0.1 m, q = ½ρU² = 60 Pa. 1. Front (θ = 180°):
    C_p = 1 − 4·0 = 1 → p − p∞ = +60 Pa. 2. Shoulder (θ = 90°): C_p = 1 − 4 = −3 → −180 Pa. 3. Back (θ = 0): +60 Pa. 4.
    Drag = −∮(p − p∞) cos θ a dθ: the pair of points θ and π − θ have the same p but opposite cos θ, so they cancel in
    pairs → D = 0. 5. Lift: θ and −θ have the same p and opposite sin θ → L = 0.")`
48. `nb.code` — *code:* `f = pf.cylinder(10.0, 0.1)` · `th = np.radians([180.0, 135.0, 90.0, 45.0, 0.0])` ·
    `print(ch06.cylinder_surface_cp(th))` · `circle = 0.1*np.exp(1j*np.linspace(0, 2*np.pi, 64, endpoint=False))` ·
    `print(ch06.surface_pressure_force(f, circle, rho=1.2))` · `print(pf.normal_velocity_on(f, circle))`. *expect:* `[1.
    -1. -3. -1. 1.]` · `(≈0, ≈0)` (≤ 1e-12 × ½ρU²a = 6 N/m) · ≤ 1e-13. *explain:* 1. the cylinder object; 2. (6.35) at
    five angles; 3. the pressure force by the periodic trapezoid on 64 surface points: zero drag and lift; 4. the surface
    is a streamline.
49. `nb.check_agree` — **from scratch (curation §7):** `th = np.linspace(0, 2*np.pi, 64, endpoint=False)`; `p = 0.5*1.2*
    10.0**2*(1 - 4*np.sin(th)**2)` (p − p∞ from (6.35)); `D_hand = -np.sum(p*np.cos(th)*0.1)*(2*np.pi/64)` (−∮p n_x dl);
    `assert np.isclose(D_hand, ch06.surface_pressure_force(f, circle, rho=1.2)[0], atol=1e-10)`. Markdown: "Adding the
    64 pressure pushes by hand gives zero too."
50. `nb.figure` — two panels (9 × 3.6 in): left, the cylinder's streamlines (`flow_net`) with surface arrows −(p − p∞)n
    (orange pushing in where p > p∞, blue pulling out where suction) from `pressure_arrows`; right, C_p against the
    angle **from the front stagnation point** (0–180°, the Fig. 6.10 convention, labelled "angle from the front = 180° −
    θ") — the ideal curve 1 − 4 sin²θ (black) and the **qualitative** separated band from `separated_cp_band` (grey,
    labelled "qualitative sketch of a measured high-Re curve — not data"). Title "Symmetric pressure, no drag — until
    the flow separates". *see:* "arrows mirror-symmetric front/back; the grey band follows the ideal curve on the front
    and stays low on the back." *read:* "the ideal back half pushes forward as hard as the front pushes back; a real wake
    keeps the back at suction, so the front wins: drag." *change:* "…the flow separated earlier (drag the band's
    separation angle in the plotly below): the low-pressure plateau widens and the imbalance — the drag — grows."
51. `nb.plotly` — `slider_figure` over the separation angle 60°–120° from the front (13 steps): traces "ideal 1 − 4
    sin²θ" (black), "qualitative separated band" (grey fill), with the net "drag coefficient of the sketch" printed in
    the title computed by the trapezoid. Title "Where ideal and real flow part company (qualitative)".
52. `nb.note` — **N31 [B]** "**Ideal vs real.** A measured pressure distribution on a cylinder at high Reynolds number
    follows 1 − 4 sin²θ on the front, but behind the widest point the boundary layer separates and the rear pressure
    stays low (the book compares with a measured curve; no dataset is cited here, so we draw only a labelled band). ⚠️
    That comparison measures the angle from the upstream stagnation point, 180° − θ in our convention. Ch. 9 explains
    separation."
53. `nb.md` — "> ⚠️ **Common confusion:** 'no drag means no pressure'. The pressures are large (−180 Pa at the shoulders
    in the example); they just cancel in total. Each half of the cylinder feels a big force; the halves fight each other
    to a draw."
54. `nb.md` — **What would change if…** "…we added circulation? The top would speed up and the bottom slow down; the
    front–back symmetry survives, the top–bottom one does not — lift without drag (C07)."

**C07 — The cylinder with circulation: stagnation points move, lift L = ρUΓ (6.40)**
55. `nb.core("C07", "Circulation gives lift: $L=\\rho U\\Gamma$ (6.40)", question="A spinning ball curves; a rotor ship
    sails on spinning towers. How does circulation turn into a sideways force, and how big is it?")`
56. `nb.md` — **The problem in plain words:** "A football struck with spin bends in flight; Flettner's rotor ships used
    tall spinning cylinders instead of sails. In both, the air is carried round the body on one side and held back on
    the other. The ideal-flow model adds a vortex to the cylinder: the stagnation points slide round, the pressure drops
    on the fast side — and the net force comes out as a strikingly simple product."
57. `nb.md` — **The idea** (ASCII):
    ```
    stream →        + clockwise vortex        =  faster on top, slower below
       ( ○ )             ↻ (−Γ circulation)       stagnation points slide DOWN: sin θ = −Γ/4πaU
                                                   pressure lower on top ⇒ force UP:  L = ρUΓ, D = 0
    ```
58. `nb.md` — **⚠️ Convention callout (numbers).** "The book now names a **clockwise** circulation Γ (the flow's
    circulation is −Γ; its footnote says the sign gives the usual (6.40)). Our code keeps Γ counterclockwise in
    `pf.Vortex` and asks body functions for `Gamma_cw=` (book) or `Gamma_ccw=` (project). Γ_cw = 2 m²/s ⇔ Γ_ccw = −2 m²/s:
    stagnation points at −9.16° and −170.84° (below the axis), L = +24 N/m for U = 10 m/s, a = 0.1 m, air."
59. `nb.note` — **N32 [B]** "Stream + doublet + a clockwise vortex at the centre:" equation
    `\psi=U\Big(r-\frac{a^2}{r}\Big)\sin\theta+\frac{\Gamma}{2\pi}\ln\Big(\frac ra\Big)`, ref "6.36". "The a inside the
    logarithm only adds a constant to ψ, so that ψ = 0 on r = a still."
60. `nb.md` — 🔁 reminders: arcsin and the two solutions of sin θ = s (gloss: "θ and π − θ"), quadratic formula and
    Vieta's product of roots (P71), `brentq` (P108).
61. `nb.derivation("D10", …)` — Part F D10 (9 steps), ref "6.38".
62. `nb.note` — **N33 [B]** "Surface speed:" equation `u_\theta(r=a,\theta)=-2U\sin\theta-\Gamma/2\pi a`, ref "6.37". "On
    top (θ = 90°) the fluid moves at 2U + Γ/2πa, below at 2U − Γ/2πa: 23.18 and 16.82 m/s in the example."
63. `nb.note` — **N34 [B]** "**Stagnation points:**" equation `\sin\theta=-\Gamma/4\pi aU`, ref "6.38". "Two on the surface
    for Γ < 4πaU, one (at the bottom) when Γ = 4πaU, one off the body on the negative y-axis at
    $r=\frac1{4\pi U}[\Gamma+\sqrt{\Gamma^2-(4\pi aU)^2}]$ for Γ > 4πaU. U = 10 m/s, a = 0.1 m, Γ = 2 m²/s: sin θ = −0.159,
    θ = −9.16° and −170.84°; the critical circulation 4πaU = 12.57 m²/s; at Γ = 6πaU the free point sits at r = 0.2618 m
    (its partner inside at 0.0382 m, product a² = 0.01 m²)."
64. `nb.primer("Newton's method for complex zeros", "To find where a complex function f(z) vanishes, start at a guess
    z₀ and repeat $z\\leftarrow z-f(z)/f'(z)$: each step follows the tangent line to zero, and near a simple root the
    number of correct digits doubles each time. Here f = dw/dz (a stagnation point is where u − iv = 0) and f′ = d²w/dz².
    To find a second root, divide f by (z − z₁) ('deflation') so the first is not found again.", code="f = lambda z:
    z**2 + 1                      # zeros at ±i\nfp = lambda z: 2*z                         # its derivative\nz = 0.5 + 0.5j
    # a guess in the upper half\nfor _ in range(6): z = z - f(z)/fp(z)   # Newton steps\nprint(z)
    # 1j (to machine precision)")` — *expect:* `1j` (≈ 1e-16 real part).
65. `nb.note` — **N35 [B]** "**Surface pressure** from (6.18) with the constant p∞ + ½ρU² and (6.37):" equation
    `p(r=a,\theta)=p_\infty+\tfrac12\rho\Big[U^2-\Big(-2U\sin\theta-\frac{\Gamma}{2\pi a}\Big)^2\Big]`, ref "6.39".
66. `nb.note` — **N36 [B]** "**The lift integral.** With the outward normal n = e_r and the arc element dl = a dθ, the
    upward force is" equation `L=-\int_0^{2\pi}p(r=a,\theta)\,\mathbf n\,dl\cdot\mathbf e_y=-\int_0^{2\pi}p(r=a,\theta)\sin\theta\,a\,d\theta`
    (no number; it sits between (6.39) and (6.40)). "The minus sign: pressure on the upper surface, where n has an upward
    component, pushes the cylinder **down**."
67. `nb.derivation("D11", …)` — Part F D11 (9 steps), ref "6.40".
68. `nb.worked_example("a spinning cylinder in wind", "U = 10 m/s, a = 0.1 m, Γ = 2 m²/s (clockwise), air ρ = 1.2 kg/m³.
    1. Critical Γ = 4πaU = 12.57 m²/s — we are below it. 2. sin θ = −2/12.566 = −0.159 → θ = −9.16° and −170.84°. 3. Top
    speed 2U + Γ/2πa = 20 + 3.18 = 23.18 m/s; bottom 20 − 3.18 = 16.82 m/s. 4. Lift L = ρUΓ = 1.2 × 10 × 2 = 24 N per
    metre of cylinder, upward. 5. Drag 0.")`
69. `nb.code` — *code:* `f = pf.cylinder(10.0, 0.1, Gamma_cw=2.0)` · `print(ch06.cylinder_stagnation_points(10.0, 0.1,
    Gamma_cw=2.0))` · `print(np.degrees(np.angle(f.stagnation_points())))` · `th = np.linspace(0, 2*np.pi, 64,
    endpoint=False)`; `circle = 0.1*np.exp(1j*th)` · `print(ch06.surface_pressure_force(f, circle, rho=1.2),
    ch06.lift_per_span(1.2, 10.0, Gamma_cw=2.0))` · `print(f.velocity(0.0, 0.1), pf.cylinder(10.0, 0.1,
    Gamma_ccw=2.0).velocity(0.0, 0.1))` (the top of the cylinder, both signs). *expect:* two complex points at angles
    −9.158°, −170.842° (on \|z\| = 0.1) · the same angles from Newton · `(≈0, 24.0) 24.0` · `(23.183, ≈0)` and `(16.817,
    ≈0)`: the clockwise Γ speeds the top up, the counterclockwise one slows it (and gives L = −24 N/m). *explain:* 1. the
    cylinder with the book's clockwise Γ; 2. the closed-form stagnation points (6.38) and the Newton cross-check agree;
    3. the pressure integral gives L = ρUΓ = 24 N/m and D = 0; 4. a counterclockwise Γ of the same size slows the top
    and pushes the cylinder down.
70. `nb.check_agree` — **from scratch (curation §7):** `U, a, G, rho = 10.0, 0.1, 2.0, 1.2`; `th = np.linspace(0,
    2*np.pi, 64, endpoint=False)`; `p = 0.5*rho*(U**2 - (-2*U*np.sin(th) - G/(2*np.pi*a))**2)` ((6.39) with p∞ = 0);
    `L_hand = -np.sum(p*np.sin(th)*a)*(2*np.pi/64)` (the lift integral of N36); `assert np.isclose(L_hand, rho*U*G)` and
    `np.isclose(L_hand, ch06.lift_per_span(rho, U, Gamma_cw=G))`. Markdown: "64 pressure pushes added by hand give 24
    N/m to 13 digits: ρUΓ."
71. `nb.note` — **N38 [B]** "**Which Γ? Uniqueness needs topology.** In a region with no holes an ideal flow is fixed by
    its boundary conditions. Round a cylinder the region has a hole: every member of (6.36) has u·n = 0 on r = a and
    U far away, whatever Γ — the boundary conditions cannot choose. Real flow chooses, through viscosity at a sharp
    trailing edge: the Kutta condition (Ch. 14). (Simply connected regions: Ch. 3 §3.4.)" + `nb.code`:
    `chk = ch06.circulation_family_check(10.0, 0.1, [0.0, 2*np.pi*1.0, 4*np.pi*1.0, 8*np.pi*1.0])` (Γ = 0, 2πaU, 4πaU,
    8πaU with aU = 1) · `print(chk)`. *expect:* max_normal ≤ 1e-13 for all four; far_field 0, 0.01, 0.02, 0.04 m/s on r = 1000a = 100 m
    (the vortex part Γ/2πR decays only like 1/R — every Γ still returns to the stream); circulation = −Γ_cw on every
    circle (0, −6.28, −12.57, −25.13 m²/s).
72. `nb.figure` — **four regimes** (2 × 2 panels, 8 × 7 in): Γ = 0, 2πaU, 4πaU, 6πaU (U = 1, a = 1); streamlines, the
    cylinder, stagnation points (black dots), the free stagnation point's crossing streamlines in the last. Title "The
    stagnation points slide down, meet, and leave the body". *see:* "two points at the sides, lower, merged at the
    bottom, and a free point below." *read:* "sin θ = −Γ/4πaU: at Γ = 2πaU they sit at −30° and −150°; at 4πaU at −90°;
    above that the root r₊ = [Γ + √(Γ² − (4πaU)²)]/4πU = 2.618a." *change:* "…U doubled at fixed Γ: the points move back
    up (Γ/4πaU halves) and the lift doubles (ρUΓ)."
73. `nb.plotly` — `slider_figure` over Γ/(4πaU) from 0 to 2 (21 steps, FAST 11): traces "surface C_p" (black, from
    `cylinder_surface_cp`) against θ, "Γ = 0 ghost" (grey dashed), stagnation angles as vertical markers; the title
    prints L = ρUΓ and the regime word from `cylinder_circulation_state`. Title "Circulation shifts the suction to the
    top".
74. `nb.explainer("cylinder_circulation_lift", heading="How does spin turn into lift?", why="Figure 6.12's four
    snapshots sample a continuous story: dragging Γ moves the stagnation points round the body, makes them meet at
    4πaU and leave it, while the pressure arrows and the force bars move with them — the lift bar lands on ρUΓ and the
    drag bar never leaves zero.", tries=["Drag Γ/4πaU from 0 to 1: the two stagnation points slide down and meet.",
    "Keep going to 1.5: a free stagnation point leaves the body — find its radius in Explain.", "Switch the convention
    chips to 'project Γ (ccw)': the same picture, the numbers change sign.", "Open Derivation → D11 step 7: the only
    term of the pressure that survives the lift integral lights up on the C_p curve."])`
75. `nb.note` — **N37 [C]** "**History and spin.** Kutta and Zhukhovsky found (6.40) independently just after 1900.
    A rotating cylinder develops circulation through its boundary layer; at high spin the real flow resembles the
    Γ > 4πaU pattern, at low spin it separates. The Magnus effect on spinning balls is mostly delayed separation (Ch. 9
    §9.9); how a sharp edge fixes Γ is Ch. 14."
76. `nb.md` — **What would change if…** "…the body were an ellipse or an airfoil instead of a circle? The lift would
    still be ρUΓ — the shape only changes *which* Γ nature picks. C10 proves it for any body; C11 builds the ellipse."

77. `nb.recap("R15", "Images for circles", "A circle can be made a streamline by more than one image (a vortex outside a
    cylinder needs an image vortex at the inverse point and one at the centre), and images move when the flow is
    unsteady; Ch. 5's `circle_image_system` does it. Milne-Thomson's circle theorem, $w(z)+\\overline{w(a^2/\\bar z)}$, is the
    general rule (our addition, `pf.circle_theorem`).", where="Ch. 5 §5.7")`
78. `nb.recap("R16", "A vortex beside a wall", "A vortex of strength −Γ (clockwise) at (h, 0) beside the wall x = 0 is
    modelled by the vortex plus an image of the opposite sign at (−h, 0); then u = 0 on the wall.", where="Ch. 5 §5.7")`
79. `nb.recap("R17", "Its drift", "The vortex moves with the velocity its image induces at its centre (its own is taken
    as zero): $d\\xi_x/dt=0$, $d\\xi_y/dt=\\Gamma/4\\pi\\xi_x$, so $\\boldsymbol\\xi(t)=(h,\\ \\Gamma t/4\\pi h)$ — along the wall at
    Γ/4πh (0.0796 m/s for Γ = 1 m²/s, h = 1 m).", where="Ch. 5 §5.7")` + `nb.code`: `from fluidpy.core import
    biot_savart as bs` · `sol = bs.point_vortex_evolve(np.array([[0.0], [1.0]]), np.array([-1.0]), np.array([0.0, 10.0]),
    boundary="wall")` (Ch. 5's walls are horizontal, y = 0, so we turn Example 6.1 by 90° counterclockwise: the wall
    x = 0 becomes y = 0, the vortex at (h, 0) moves to (0, h), and the drift along +y becomes a drift along −x) ·
    `print(sol[-1])` → *expect:* x = −0.7958 m (= −Γt/4πh), y = 1.0 m — the drift Γ/4πh = 0.0796 m/s along the wall.

**C08 — The method of images and Example 6.1's wall-pressure signal**
80. `nb.core("C08", "Images, and what the wall feels as a vortex passes: $\\frac{p(0,0,t)-p_\\infty}{\\rho}=\\frac{\\Gamma^2}{4\\pi^2}
    \\frac{(\\Gamma t/4\\pi h)^2-h^2}{((\\Gamma t/4\\pi h)^2+h^2)^2}$ (Example 6.1)", question="An eddy sweeps along a wall. What does
    a pressure sensor on the wall record as it passes?")`
81. `nb.md` — **The problem in plain words:** "An eddy rolls along a seabed, a harbour wall or the side of a building;
    a pressure sensor mounted in the wall records a blip. The wall is solid, so the flow must slide along it — and a
    mirror trick, the image vortex, makes that automatic. With the image in place, unsteady Bernoulli tells us the whole
    pressure signal: suction as the eddy arrives, a small over-pressure after it has passed."
82. `nb.md` — **The idea** (ASCII):
    ```
     image +Γ ↺  │  ↻ vortex −Γ        mirror images: vortex → opposite sign, source → same sign
       (−h, s)   │    (h, s)            the pair drifts along the wall at Γ/4πh (each pushes the other)
                wall x = 0               p at the wall = −ρ ∂φ/∂t (unsteady, +) − ½ρv² (speed, −)
    ```
83. `nb.md` — 🔁 reminders + glosses: odd and even functions under a mirror ("ψ(x, −y) = −ψ(x, y) makes ψ = 0 on y = 0;
    an even φ has ∂φ/∂y = 0 there"), chain rule (P49), derivative of an arctangent with a moving argument (gloss:
    $\frac{d}{dt}\tan^{-1}(Y/X)=\frac{X\dot Y}{X^2+Y^2}$ at fixed X), unsteady Bernoulli (R32 later; Ch. 4 §4.9 (4.83)).
84. `nb.derivation("D12", …)` — Part F D12 (5 steps), ref "" (the image rules of Figs. 6.14–6.15 carry no number).
85. `nb.note` — **N39 [B]** "**The rules.** Wall y = 0: vortex images flip sign, $\psi_2=\psi_1(x,y)-\psi_1(x,-y)$; source
    images keep it, $\phi_2=\phi_1(x,y)+\phi_1(x,-y)$ — two different conditions (ψ₂ = 0, ∂φ₂/∂y = 0) that both make the wall
    impermeable." + `nb.code`: `els = pf.mirror([pf.Vortex(1.0, 0.5+1.0j), pf.Source(2.0, -1.0+0.7j)], wall="y=0")` ·
    `f = pf.Flow(els)` · `xw = np.linspace(-5, 5, 1001)` · `print(np.max(np.abs(f.velocity(xw, 0*xw)[1])))` → *expect:*
    ≤ 1e-14 (no v on the wall).
86. `nb.note` — **N40 [B]** "**A source by a wall = two sources.** Equal sources at x = ±a:" equation
    `\psi=\frac m{2\pi}\Big[\tan^{-1}\Big(\frac y{x+a}\Big)+\tan^{-1}\Big(\frac y{x-a}\Big)\Big],\quad \phi=\frac m{2\pi}\Big[\ln\sqrt{(x+a)^2+y^2}+\ln\sqrt{(x-a)^2+y^2}\Big]`,
    ref "6.41". "With the tangent addition formula the streamlines are $x^2-y^2-2xy\cot(2\pi\psi/m)=a^2$. The axes are
    streamlines and the origin a stagnation point, so the same picture is two sources, one source beside a wall (the
    y-axis), or a slit flow into a right-angled corner." + `nb.code`: `f2 = ch06.two_sources(2*np.pi, 1.0)` · `x = np.linspace(1.2, 4, 5)` ·
    `y = ch06.two_source_streamline(np.pi/2, 2*np.pi, 1.0, x)` · `print(f2.psi(x, y))` → *expect:* all π/2 (≤ 1e-12,
    modulo m where the branch is crossed — the builder prints `np.mod`).
87. `nb.derivation("D13", …)` — Part F D13 (11 steps), ref "" (the title says "Example 6.1").
88. `nb.worked_example("an eddy by a harbour wall", "Water ρ = 1000 kg/m³, Γ = 1 m²/s, h = 1 m. 1. Drift Γ/4πh =
    1/12.566 = 0.0796 m/s. 2. At t = 0 (eddy level with the sensor): $\\frac{p-p_\\infty}{\\rho}=\\frac{1}{39.48}\\cdot\\frac{0-1}{1}$ →
    p − p∞ = −25.33 Pa: suction. Split: the unsteady part −ρ∂φ/∂t = +25.33 Pa, the speed part −½ρv² = −50.66 Pa. 3. The
    signal crosses zero when Γt/4πh = h, i.e. t = 4πh²/Γ = 12.57 s (the eddy is 1 m past the sensor, at 45°). 4. The
    largest over-pressure comes at Γt/4πh = √3 h, t = 21.77 s: p − p∞ = ρΓ²/(32π²h²) = +3.17 Pa. 5. Then it fades as
    t⁻².")`
89. `nb.code` — *code:* `for t in (0.0, 4*np.pi, 4*np.sqrt(3)*np.pi, 60.0): r = ch06.example_6_1(t); print(t,
    r["xi_y"], r["p_origin"], r["unsteady"], r["speed_part"])` · `print(ch06.example_6_1(5.0, route="numeric")["p_origin"],
    ch06.example_6_1(5.0)["p_origin"])`. *expect:* t = 0: 0, −25.330, +25.330, −50.661 · t = 12.566: 1.0, 0.0 (≤ 1e-12),
    +12.665, −12.665 · t = 21.766: 1.732, +3.1663 · t = 60: 4.775, +0.975 (fading) · numeric vs closed agree to 1e-7
    (central difference in t, step 1e-3 s). *explain:* 1. the closed form of Example 6.1 and its split into the unsteady
    and the speed parts; 2. the numeric route moves a `Vortex` and its image, takes ∂φ/∂t by a central difference in time
    and uses Bernoulli — two independent routes, one number.
90. `nb.check_agree` — **from scratch (curation §7):** `Gm, h, rho = 1.0, 1.0, 1000.0`; `phi = lambda t:
    (Gm/(2*np.pi))*(np.arctan2(0 - Gm*t/(4*np.pi*h), 0 + h) - np.arctan2(0 - Gm*t/(4*np.pi*h), 0 - h))` (φ at the origin,
    D13 step 4, arctan2 for a branch that does not jump here); `dt = 1e-4`; `dphidt = (phi(5 + dt) - phi(5 - dt))/(2*dt)`;
    `v = Gm*h/(np.pi*(h**2 + (Gm*5/(4*np.pi*h))**2))`; `p = -rho*(dphidt + 0.5*v**2)`; `assert np.isclose(p,
    ch06.example_6_1(5.0)["p_origin"], rtol=1e-6)`. Markdown: "Bernoulli by hand at one instant gives the library's
    pressure."
91. `nb.figure` — two panels (9 × 3.4 in): left, the wall x = 0 (hatched), the vortex path (ξ = (1, Γt/4πh), amber),
    the image path (faint amber), and streamlines at t = 10 s; right, p(0, 0, t) − p∞ for t ∈ [−60, 60] s (the formula is
    even in t: the eddy approaching from below and leaving above), with markers at t = 0 (−25.3 Pa), ±12.57 s (zero),
    ±21.77 s (+3.17 Pa), and dashed curves for the unsteady (orange, +) and speed (blue, −) parts. Title "Suction as it
    passes, a small push before and after". *see:* "a deep negative dip flanked by two small positive humps." *read:* "the
    speed part is always suction; the unsteady part is always positive here and wins when the eddy is far (it decays
    like t⁻², the speed part like t⁻⁴)." *change:* "…h halved: the dip deepens ×4 (Γ²/h² scaling) and the whole signal
    runs 4× faster (time scale h²/Γ)."
92. `nb.figure` — **two sources** (Fig. 6.16 idea): ψ contours of `ch06.two_sources(2π, 1)` with the (6.41) curves for
    four ψ values overlaid as dashed lines, the y-axis drawn as a wall. *see / read / change* — "the dashed algebraic
    curves lie on the contours; the right half is a source beside a wall; moving the sources apart pushes the origin's
    stagnation point… it stays at the origin by symmetry."
93. `nb.animation` — **the eddy passes the sensor** (`player="video"`, 60 frames, FAST 30): the vortex and its image
    moving up, the wall points coloured by p − p∞ from `ch06.example_6_1_wall_pressure(y_wall, t)` (orange > 0, blue <
    0), and the p(0, 0, t) trace drawing itself in a second panel with a moving dot. + see/read/change markdown: "the
    blue patch travels with the eddy; orange patches lead and trail it."
94. `nb.explainer("vortex_wall_images", heading="What does a wall feel as an eddy passes?", why="The result is a time
    series tied to a moving vortex: scrubbing time links the vortex position to the sensor trace and to the pressure along
    the whole wall, and the source mode shows the other image rule on the same stage.", tries=["Press ▶ and watch the
    blue suction patch ride along the wall with the eddy.", "Preset 't = 4πh²/Γ': the signal crosses zero — the terms
    panel shows the unsteady and speed parts cancelling.", "Halve h: the dip is four times deeper and four times
    shorter.", "Switch to 'source by a wall': the image keeps its sign, and the wall still has no flow through it."])`
95. `nb.md` — **What would change if…** "…the wall were the sea floor and the eddy a storm-driven vortex in the ocean?
    The same image model gives the bottom-pressure signal that ocean-bottom pressure recorders see as eddies pass —
    a small signal (tens of pascals here) compared with the tides, which is why it is averaged out or modelled (Ch. 13
    adds rotation and stratification)."

### A.4 §6.4 Complex Potential — C09 (+N41–N52 · D14 D15 · E4)
1. `nb.section("6.4", "Complex Potential", intro="**What is this section about?** φ and ψ are two halves of one complex
   function w(z) of the complex position z = x + iy. If w has a derivative that does not depend on the direction you
   take it in, its real and imaginary parts are automatically a potential and a stream function of an ideal flow, and
   the derivative is the velocity — with v flipped: dw/dz = u − iv. Every element of §6.3 becomes one line.")`

**C09 — The complex potential w = φ + iψ, Cauchy–Riemann and the complex velocity (6.42)–(6.46)**
2. `nb.core("C09", "The complex potential $w\\equiv\\phi+i\\psi$ (6.42) and the complex velocity $dw/dz=u-iv$ (6.45)",
   question="Why does a single complex function carry a whole flow — and why is the velocity u − iv, not u + iv?")`
3. `nb.md` — **The problem in plain words:** "So far every flow needed two functions, one for streamlines and one for
   equipotentials, and every new body meant new algebra. Engineers designing wings in the 1900s found a shortcut: write
   the position as one complex number and the flow as one complex function. Then any 'nice' function of z is a flow,
   forces become a single contour integral (C10), and maps between planes carry flows onto new shapes (C11)."
4. `nb.md` — **The idea** (ASCII):
   ```
   z = x + iy ──► w(z) = φ + iψ          Re w: equipotentials (dashed)   Im w: streamlines (solid)
   derivative along x  =  derivative along iy   ⇔   Cauchy–Riemann   ⇔   φ ⟂ ψ and both harmonic
   dw/dz = u − iv   (the velocity, mirrored in the x-axis)
   ```
5. `nb.primer("the complex plane in numpy", "A complex number z = x + iy is a point (x, y); its modulus \\|z\\| = √(x² + y²) is
   the distance from 0 and its argument arg z = atan2(y, x) the angle, so $z=re^{i\\theta}$ (Euler's formula, P45).
   Multiplying two complex numbers multiplies their moduli and **adds their angles** — a rotation plus a stretch. The
   conjugate z* = x − iy mirrors in the x-axis (P81). numpy: `1j`, `np.abs`, `np.angle`, `np.conj`; arrays of complex
   numbers work like any other array.", code="z = 1 + 1j                              # the point (1, 1)\nprint(np.abs(z),
   np.degrees(np.angle(z)))   # 1.414, 45°\nq = z * 2j                               # times 2i: stretch ×2, turn +90°\nprint(
   q, np.abs(q), np.degrees(np.angle(q)), np.conj(z))   # (-2+2j) 2.828 135.0 (1-1j)")` — *expect:* `1.4142135623730951
   45.0` · `(-2+2j) 2.8284271247461903 135.0 (1-1j)`.
6. `nb.note` — **N41 [B]** "The complex position in Cartesian and polar form:" equation `z\equiv x+iy=re^{i\theta}`, ref
   "6.43". "r = (x² + y²)^½, θ = tan⁻¹(y/x) (numpy: `np.abs`, `np.angle`)."
7. `nb.primer("complex derivative and analytic functions", "The derivative of w(z) is the limit of $[w(z+\\delta z)-w(z)]/\\delta z$
   as the small complex step δz shrinks — but δz can point in any direction in the plane. If the limit is the same for
   every direction, w is **analytic** there. Polynomials, exp, log (away from its cut) and 1/z (away from 0) are
   analytic; z* is not: along x its quotient is 1, along iy it is −1.", code="w = lambda z: z**2
   # an analytic function\nz0, h = 1 + 1j, 1e-6                     # a point and a small step\nprint((w(z0 + h) - w(z0))/h,
   (w(z0 + 1j*h) - w(z0))/(1j*h))   # both ≈ 2+2j\ns = np.conj                                # not analytic\nprint((s(z0 +
   h) - s(z0))/h, (s(z0 + 1j*h) - s(z0))/(1j*h))   # 1 and -1: no derivative")` — *expect:* `(2+2j)` (to 1e-6) twice ·
   `(1+0j) (-1+0j)`.
8. `nb.derivation("D14", …)` — Part F D14 (8 steps), ref "6.44".
9. `nb.note` — **N42 [B]** "**Cauchy–Riemann:**" equation `\frac{\partial\phi}{\partial x}=\frac{\partial\psi}{\partial y},\qquad \frac{\partial\phi}{\partial y}=-\frac{\partial\psi}{\partial x}`,
   ref "6.44". "They make the φ- and ψ-lines cross at right angles — except where w or dw/dz is zero or infinite
   (stagnation points and singular points)."
10. `nb.note` — **N43 [B]** "**The complex velocity:**" equation `\frac{dw}{dz}=u-iv`, ref "6.45". "Applying (6.44) to it
    gives back continuity (6.2) and irrotationality (6.9), and Laplace's equation for both φ and ψ: *any analytic function
    of z is a plane ideal flow.* ⚠️ The book writes 'Laplace equations for φ and ψ, (6.5) and (6.12), respectively' — it
    is the other way round ((6.5) is ψ's, (6.12) is φ's); and '(6.43) ensures the equality of the u, v components' means
    the Cauchy–Riemann conditions (6.44)."
11. `nb.worked_example("w = z² at z = 1 + i", "1. $dw/dz=2z=2+2i$. 2. So $u-iv=2+2i$: u = 2 m/s and **v = −2 m/s** (the
    sign flips). 3. Check with ψ = Im z² = 2xy (the corner flow (6.24) with A = 1): $u=\\partial\\psi/\\partial y=2x=2$ ✓,
    $v=-\\partial\\psi/\\partial x=-2y=-2$ ✓. 4. Cauchy–Riemann: φ = x² − y², $\\partial\\phi/\\partial x=2x=2=\\partial\\psi/\\partial y$ ✓,
    $\\partial\\phi/\\partial y=-2y=-2=-\\partial\\psi/\\partial x$ ✓.")`
12. `nb.primer("complex logarithm, powers and branch cuts", "Because θ is only defined up to multiples of 2π, $\\ln
    z=\\ln r+i\\theta$ has many values; a program picks θ in a range (numpy: (−π, π]), so ln z jumps by 2πi across a **branch
    cut** (numpy's runs along the negative real axis). Powers inherit it: $z^n=e^{n\\ln z}=r^ne^{in\\theta}$. In a flow we
    rotate the cut so it lies inside a body or outside the fluid we draw — then every quantity we look at is smooth.",
    code="z = np.array([-1 + 1e-9j, -1 - 1e-9j])   # just above and below the negative real axis\nprint(np.log(z).imag)
    # +π and −π: a jump of 2π across the cut\nprint(np.sqrt(z))                         # ±i: the square root jumps
    too\nprint(np.log(z * np.exp(-1j*np.pi/2)).imag + np.pi/2)   # cut turned to the −y axis: no jump at z = −1")` —
    *expect:* `[ 3.14159265 -3.14159265]` · `[0.+1.j 0.-1.j]` (to 1e-9) · `[3.14159265 3.14159265]`.
13. `nb.derivation("D15", …)` — Part F D15 (7 steps), ref "6.46".
14. `nb.note` — **N44 [B]** "**Flow in a corner of angle α = π/n:**" equation `w(z)=Az^n=A(re^{i\theta})^n=Ar^n(\cos n\theta+i\sin n\theta),\quad n\ge\tfrac12`,
    ref "6.46". "Its velocity $dw/dz=nAz^{n-1}=(A\pi/\alpha)z^{(\pi-\alpha)/\alpha}$ is zero at the corner when α < π (a
    stagnation point) and infinite when α > π. n = 2 is (6.24)/(6.27); n = ½ is the flow round the edge of a flat plate.
    ⚠️ Ch. 9's wedge flows U ∝ x^m start from exactly this outer flow." + `nb.code`: `for n in (4, 2, 1, 2/3, 0.5):
    print(n, np.degrees(np.pi/n), ch06.corner_speed_exponent(n))` → *expect:* 45° 3 · 90° 1 · 180° 0 · 270° −1/3 · 360°
    −½.
15. `nb.note` — **N45 [B] – N51 [B]** "**Every element in one line.**" (a table; each row equation + number): vortex
    $w=-\frac{i\Gamma}{2\pi}\ln(z-z')=\frac{\Gamma}{2\pi}\theta'-i\frac{\Gamma}{2\pi}\ln r'$ (6.47) (⚠️ counterclockwise Γ; its φ
    = Γθ′/2π jumps by Γ across the cut — the circulation) · source $w=\frac{m}{2\pi}\ln(z-z')=\frac{m}{2\pi}\ln r'+i\frac{m\theta'}{2\pi}$
    (6.48) · doublet $w=\frac{d}{2\pi(z-z')}$ (6.49) (⚠️ scalar d = a dipole −d e_x; Re w is the φ of (6.29)) · half-body
    $w=Uz+\frac{m}{2\pi}\ln z$ (6.50) (Im w = (6.31)) · cylinder $w=U(z+a^2/z)$ (6.51) (Im w = (6.33)) · cylinder with a
    clockwise circulation $w=U(z+\frac{a^2}{z})+\frac{i\Gamma}{2\pi}\ln(z/a)$ (6.52) (Im w = (6.36)) · two sources
    $w=\frac m{2\pi}\ln\big(\frac{z^2-a^2}{a^2}\big)$ (6.53) (ln A + ln B = ln AB up to 2πi; Im w = (6.41) ψ mod m). "Each
    row is one line of algebra: $\ln(z-z')=\ln r'+i\theta'$ and $\frac1{z-z'}=\frac{(x-x')-i(y-y')}{r'^2}$."
16. `nb.code` — *code:* `z = np.array([1+1j, -0.5+2j, 2-0.3j])` · `flows = {"z^2": pf.Corner(1.0, n=2), "vortex":
    pf.Vortex(2*np.pi), "source": pf.Source(2*np.pi), "doublet": pf.Doublet.from_book_scalar(2*np.pi), "cylinder":
    pf.cylinder(1.0, 1.0, Gamma_cw=2.0)}` · `for k, e in flows.items(): r1, r2 = ch06.cauchy_riemann_residual(e.w, z);
    u, v = e.velocity(z.real, z.imag); print(k, np.max(np.abs(r1)), np.max(np.abs(e.dwdz(z) - (u - 1j*v))))` · `print(
    ch06.cauchy_riemann_residual("conj", z))`. *expect:* residuals ≲ 1e-8, dw/dz − (u − iv) ≲ 1e-12 for all five; `conj`:
    r1 = 2 everywhere (φ_x − ψ_y = 1 − (−1)) — "not a flow". *explain:* 1. five flows as complex functions; 2.
    Cauchy–Riemann holds for each (to difference-quotient accuracy); 3. dw/dz equals u − iv computed from the element's
    own velocity; 4. z* fails: it is not analytic.
17. `nb.check_agree` — **from scratch (curation §7):** complex difference quotients by hand: `c = pf.Corner(1.0, n=2)`;
    `z0, h = 1 + 1j, 1e-6`; `qx = (c.w(z0 + h) - c.w(z0))/h`; `qy = (c.w(z0 + 1j*h) - c.w(z0))/(1j*h)`; `assert
    np.allclose([qx, qy], c.dwdz(z0), atol=1e-5)`; and for z*: `qx_c = (np.conj(z0 + h) - np.conj(z0))/h`; `qy_c =
    (np.conj(z0 + 1j*h) - np.conj(z0))/(1j*h)`; `assert not np.isclose(qx_c, qy_c)`. Markdown: "Two directions, one
    derivative — that is what 'analytic' means; z* gives two different answers."
18. `nb.figure` — **corner flow nets** (4 panels, 10 × 2.8 in): w = Az^n for α = π/2, π, 3π/2 and 2π (n = 2, 1, ⅔, ½),
    walls drawn black, the cut placed outside the fluid (`pf.Corner(1, n)`), ψ solid, φ dashed, speed as a faint colour.
    Title "The corner angle decides: calm or violent". *see:* "a right-angle corner, a flat wall, a re-entrant corner, a
    plate edge." *read:* "crowded streamlines near the corner mean fast flow: none in the 90° corner (stagnation),
    everywhere near the 270° and 360° corners (infinite speed at the tip)." *change:* "…α = 60° (n = 3): the fluid near
    the corner is even calmer, speed ∝ r²."
19. `nb.figure` — log–log \|dw/dz\| against r along the corner's bisector for the four n (slopes n − 1: 1, 0, −⅓, −½
    fitted with `np.polyfit`, labelled). *see / read / change* in the same pattern: "the slopes are (π − α)/α; the
    270° corner of Example 6.2 (C12) is where the grid solution will struggle."
20. `nb.plotly` — `slider_figure` over n from ½ to 4 (22 steps, FAST 11): the corner's walls and five streamlines
    (polylines from ψ = Ar^n sin nθ = c solved for r(θ)), the title printing α = π/n and the speed exponent. Title "One
    exponent, every corner".
21. `nb.note` — **N52 [C]** "**Next: forces.** Complex variables also turn force calculations into one contour integral
    (§6.5, C10)."
22. `nb.explainer("complex_potential_corners", heading="Why is the velocity u − iv?", why="Click a point and watch the two
    difference quotients — along x and along iy — converge to the same dw/dz, with the velocity drawn as its mirror
    image; drag n and the walls bend from a 45° wedge to a flat plate while the speed at the corner flips from zero to
    infinite.", tries=["Click near the corner with n = 2: dw/dz → 0 (stagnation).", "Preset 'n = ⅔ (270°)': the speed
    readout explodes as you click closer to the tip.", "Preset 'z* — not a flow': the two quotients disagree and the
    Cauchy–Riemann bars do not match.", "Watch the purple dw/dz arrow and the teal velocity arrow: mirror images in the
    horizontal."])`
23. `nb.md` — "> ⚠️ **Common confusion:** plotting dw/dz as the velocity. It is the velocity *mirrored*: u − iv. Use
    `np.conj(dwdz)` (or `(dwdz.real, -dwdz.imag)`) for arrows. A wrong sign sends every arrow the wrong way in y."
24. `nb.md` — **What would change if…** "…w were not analytic somewhere inside the region — a source or a vortex? Then
    w has a singular point there, the flow has a source of mass or of circulation, and contour integrals round it pick
    up exactly that strength (C10's residues)."

### A.5 §6.5 Forces on a Two-Dimensional Body — R18, C10 (+N53–N63 · D16 D17 D18 · E5)
1. `nb.section("6.5", "Forces on a Two-Dimensional Body", intro="**What is this section about?** The pressure round any
   2-D body, written with complex numbers, becomes one contour integral of (dw/dz)² — Blasius's theorem — and that
   integral can be moved away from the body to a large circle, where every body looks alike: a uniform stream, its
   circulation and a doublet. Only the product of the stream and the circulation survives: no drag, and lift ρUΓ for
   every shape.")`
2. `nb.recap("R18", "Momentum balance on a control volume", "For a stationary control volume in steady flow, the
   momentum carried out through the surface equals the pressure force on the surface plus the other forces F:
   $\\int_{A^*}\\rho\\mathbf u(\\mathbf u\\cdot\\mathbf n)dA=-\\int_{A^*}p\\,\\mathbf n\\,dA+\\mathbf F$ *(Eq. 6.54, from 4.17)*.",
   where="Ch. 4 §4.4")` + `nb.code`: `print(ch06.cv_force_on_body(pf.cylinder(10.0, 0.1, Gamma_cw=2.0), R_outer=1.0,
   rho=1.2))` → *expect:* `(≈0, 24.0)` — "the momentum-flux route on a far circle already gives L = ρUΓ (Exercise 6.27's
   idea); we now derive the short route".

**C10 — Blasius's theorem and Kutta–Zhukhovsky lift for any body (6.60)–(6.62)**
3. `nb.core("C10", "Blasius $D-iL=\\frac{i\\rho}{2}\\oint_C(dw/dz)^2dz$ (6.60) and Kutta–Zhukhovsky $D=0$, $L=\\rho U\\Gamma$
   (6.62) for any body", question="A round tube and a flat ellipse with the same circulation feel exactly the same
   lift. Why does the shape not matter?")`
4. `nb.md` — **The problem in plain words:** "C07 found L = ρUΓ for a circle by integrating its surface pressure. Is that
   a coincidence of the circle? Aircraft wings, turbine blades and sails are not circles, and integrating their surface
   pressure by hand is hopeless. We want a proof that works for *any* cross-section — and it comes from moving the
   integral off the body."
5. `nb.md` — **The idea** (ASCII):
   ```
   pressure on the body ──complex form──► (iρ/2)∮(dw/dz)² dz on the body     (Blasius, (6.60))
          analytic between body and a big circle ──► same integral on the big circle (Cauchy)
   far away every body = U + iΓ/2πz − d/2πz² + …  ──square──► only U·(iΓ/πz) has a 1/z ──► L = ρUΓ, D = 0
   ```
6. `nb.primer("Cauchy's integral theorem", "If f(z) is analytic everywhere inside and on a closed curve C, then ∮_C f dz
   = 0. Consequence: for two curves round the same singular region, ∮ over one equals ∮ over the other when f is analytic
   in the ring between them — a contour can be squeezed or stretched freely, as long as it does not cross a
   singularity.", code="f = lambda z: (1 + 0.3j/z - 0.2/z**2)**2      # analytic except at z = 0\nfor R in (0.5, 1.0,
   4.0):                               # three different circles round 0\n    th = np.linspace(0, 2*np.pi, 256,
   endpoint=False); z = R*np.exp(1j*th)\n    print(R, np.sum(f(z)*1j*z)*(2*np.pi/256))   # dz = iz dθ: the same number
   each time")` — *expect:* −3.769911 (+ round-off × i) for all three radii; the markdown states "2πi × (the
   coefficient of 1/z) = 2πi × 0.6i = −1.2π = −3.770".
7. `nb.note` — **N53 [B]** "**The setting.** A body of span B at rest in a steady stream; D (along x) and L (along y) are
   the forces per unit depth **on the body**; the body pushes the fluid with F = −B(De_x + Le_y) (Newton's third law).
   ⚠️ B here is a length, not Ch. 4's Bernoulli function; 'Section 3' in the book means §6.3."
8. `nb.md` — 🔁 reminders: momentum flux (P114), parametric curves (P92) + gloss "**orientation and the outward normal:**
   walking a closed curve counterclockwise, the body is on your left; the tangent (dx, dy)/ds turned clockwise by 90°,
   (dy, −dx)/ds, points out of the body (check at the rightmost point, where dy > 0 and the normal is +x)".
9. `nb.derivation("D16", …)` — Part F D16 (8 steps), ref "6.56".
10. `nb.note` — **N54 [B], N55 [B]** "On the body u·n = 0 kills the momentum flux, leaving the pressure force:" equation
    `D\mathbf e_x+L\mathbf e_y=-\frac1B\int_{A^*}p\,\mathbf n\,dA`, ref "6.55" "and with $\mathbf n=(\mathbf e_xdy-\mathbf e_ydx)/ds$
    on a counterclockwise contour" equation `D=-\oint_Cp\,dy,\qquad L=\oint_Cp\,dx`, ref "6.56".
11. `nb.derivation("D17", …, check_src=…)` — Part F D17 (12 steps, ★★★), ref "6.60". The check cell (described in
    Part F) runs sympy on the cylinder with circulation.
12. `nb.note` — **N56 [B], N57 [C], N58 [B], N59 [B], N60 [B]** (one note, four equations): "the complex force"
    `D-iL=-i\oint_Cp\,dz^*` (6.57); "Bernoulli in complex form $p_\infty+\tfrac12\rho U^2=p+\tfrac12\rho(u-iv)(u+iv)$ gives"
    `D-iL=-i\oint_C\big[p_\infty+\tfrac12\rho U^2-\tfrac12\rho(u-iv)(u+iv)\big]dz^*` (6.58); "on the body the velocity is
    tangent" `(u+iv)dz^*=(u-iv)dz=(dw/dz)dz` (6.59); "so" `D-iL=\frac{i\rho}{2}\oint_C\Big(\frac{dw}{dz}\Big)^2dz` (6.60).
    "Blasius's theorem. The contour may be any curve round the body with no singularity of (dw/dz)² in between." +
    `nb.code`: `th = np.linspace(0, 2*np.pi, 16, endpoint=False)`; `zb = 0.1*np.exp(1j*th)`; `f = pf.cylinder(10.0, 0.1,
    Gamma_cw=2.0)`; `u, v = f.velocity(zb.real, zb.imag)`; `dz = 1j*zb` (the direction of dz along the circle);
    `print(np.max(np.abs((u + 1j*v)*np.conj(dz) - (u - 1j*v)*dz)))` → *expect:* ≤ 1e-13 (6.59) pointwise on the body.
13. `nb.primer("Laurent series and residues", "Outside a disc that contains all the singular points, an analytic function
    is a sum of powers of z including negative ones: $f(z)=\\sum_k c_kz^k$ (a Laurent series). Integrating term by term
    round a circle, every power gives zero **except** $z^{-1}$, which gives 2πi — so ∮f dz = 2πi c₋₁, and c₋₁ is called the
    **residue**. Sampling f on a circle and taking an FFT reads off all the c_k at once (P142).", code="R, n = 2.0, 64
    # a circle and n samples\nz = R*np.exp(2j*np.pi*np.arange(n)/n)\nf = 3 + 0.5/z - 2/z**2                     # c0 = 3,
    c-1 = 0.5, c-2 = -2\nck = np.fft.fft(f)/n                        # FFT coefficient of e^{ikθ} …\nprint(ck[0], ck[-1]*R,
    ck[-2]*R**2)          # … rescaled by R^k: 3, 0.5, -2")` — *expect:* `(3+0j) (0.5+0j) (-2+0j)` (to 1e-15).
14. `nb.note` — **N63 [B]** "**The residue theorem on a circle.** With z = Re^{iθ}, dz = iRe^{iθ}dθ:" equation
    `\oint z^{-n}dz=iR^{1-n}\int_0^{2\pi}e^{i(1-n)\theta}d\theta=2\pi i\,\delta_{n1}` "(no number). Only n = 1 survives: the
    integral round any closed curve is 2πi times the sum of the residues inside." + `nb.code`: `for n in (0, 1, 2, 3):
    th = np.linspace(0, 2*np.pi, 64, endpoint=False); z = 1.5*np.exp(1j*th); print(n, np.round(np.sum(z**(-n)*1j*z)*
    (2*np.pi/64), 12))` → *expect:* 0, 6.283185307180j, 0, 0.
15. `nb.note` — **N61 [B]** "**Every body from far away.** Outside a circle round the body, $w=Uz+\frac{m}{2\pi}\ln
    z+\frac{i\Gamma}{2\pi}\ln z+\frac{d}{2\pi z}+\dots$ (stream, net source, clockwise vortex, doublet); a closed body has
    m = 0." + `nb.code`: `c = pf.laurent_coefficients(pf.cylinder(10.0, 0.1, Gamma_cw=2.0).dwdz, R=0.3)` · `print(c[0],
    c[-1], 1j*2.0/(2*np.pi), c[-2])` → *expect:* `(10+0j)` · `0.3183j` twice · `-0.1+0j` (c₋₂ = −Ua² = −d/2π with
    d = 2πUa² = 0.628: the cylinder's doublet).
16. `nb.primer("sympy for complex series and residues", "sympy handles the algebra of D18: `sp.I` is i, `sp.expand`
    multiplies out a squared series, `.coeff(z, -1)` reads a coefficient, and `sp.residue(f, z, 0)` returns the
    coefficient of 1/z directly.", code="z, U, G, d = sp.symbols('z U Gamma d')\nf = sp.expand((U + sp.I*G/(2*sp.pi*z) -
    d/(2*sp.pi*z**2))**2)   # the squared far field\nprint(f.coeff(z, -1), f.coeff(z, -2))    # I*Gamma*U/pi,
    -Gamma**2/(4*pi**2) - U*d/pi\nprint(sp.residue(f, z, 0))                  # I*Gamma*U/pi")` — *expect:* the three
    printed expressions as in the comments.
17. `nb.derivation("D18", …, check_src=…)` — Part F D18 (11 steps, ★★★), ref "6.62".
18. `nb.note` — **N62 [B]** "**Blasius with the far field:**" equation
    `D-iL=\frac{i\rho}{2}\oint_C\Big(U+\frac{i\Gamma}{2\pi z}-\frac{d}{2\pi z^2}+\dots\Big)^2dz`, ref "6.61". "⚠️ The book prints
    the squared series with the 1/z² coefficient $(Ud/\pi-\Gamma^2/4\pi^2)$ and an extra outer square; expanding gives
    $-(Ud/\pi+\Gamma^2/4\pi^2)$ (D18 step 6 and its sympy cell). Nothing depends on it: only the 1/z term, $iU\Gamma/\pi$,
    has a residue."
19. `nb.worked_example("the same spinning cylinder, by residues", "U = 10 m/s, Γ = 2 m²/s (clockwise), ρ = 1.2 kg/m³. 1.
    Residue of (dw/dz)² at 0: $iU\\Gamma/\\pi=i\\times20/\\pi=6.366i$ m³/s². 2. $\\oint(dw/dz)^2dz=2\\pi i\\times6.366i=-40$ m⁴/s² (=
    −2UΓ). 3. $D-iL=\\frac{i\\rho}{2}\\times(-40)=-24i$ N/m. 4. D = 0, L = 24 N/m — the same as C07's pressure
    integral.")`
20. `nb.code` — *code:* `f = pf.cylinder(10.0, 0.1, Gamma_cw=2.0)` · `for R in (0.1, 0.2, 1.0): print(R,
    pf.blasius_force(f.dwdz, R=R, rho=1.2))` · `e = ch06.elliptic_cylinder_flow(10.0, 0.12, 0.1, Gamma_cw=2.0)` (a
    Zhukhovsky ellipse, C11) · `for R in (0.5, 1.0): print(R, pf.blasius_force(e.dwdz, R=R, rho=1.2))` ·
    `print(ch06.kutta_zhukhovsky_sym())` · `print(ch06.contour_force(lambda z: f.pressure(z.real, z.imag, rho=1.2,
    U=10.0), 0.1*np.exp(1j*np.linspace(0, 2*np.pi, 128, endpoint=False))))`. *expect:* (≈0, 24.0) for the three radii ·
    (≈0, 24.0) for the ellipse on two contours (the ellipse's semi-major axis 0.2033 m; R = 0.5 and 1.0 enclose it) ·
    the sympy dict with residue `I*Gamma*U/pi`, D 0, L `Gamma*U*rho` and the printed vs true 1/z² coefficients · (≈0,
    24.0) by the pressure route. *explain:* 1. Blasius on three circles: the same force (Cauchy); 2. a different body,
    the same Γ: the same lift; 3. sympy: the residue and the misprinted coefficient; 4. the surface-pressure route (6.56)
    agrees.
21. `nb.check_agree` — **from scratch (curation §7):** `R, n = 0.3, 128`; `th = np.linspace(0, 2*np.pi, n,
    endpoint=False)`; `z = R*np.exp(1j*th)`; `integral = np.sum(f.dwdz(z)**2 * 1j*z) * (2*np.pi/n)` (∮(dw/dz)² dz with
    dz = iz dθ); `F = 1j*1.2/2*integral` (D − iL); `assert np.allclose([F.real, -F.imag], pf.blasius_force(f.dwdz, R=R,
    rho=1.2))`. Markdown: "The periodic trapezoid round a circle is Blasius's integral; for an analytic integrand it
    converges exponentially (P136)."
22. `nb.figure` — two panels (9 × 3.4 in): left, D and L from `blasius_force` against the contour radius R/a from 1 to
    20 for the cylinder (amber L, grey D) and the Zhukhovsky ellipse (dashed), flat at ρUΓ = 24 N/m and 0, with the
    pressure-route ◇ at R = a; right, the Laurent term bars from `ch06.laurent_contributions("cylinder", R=0.2)`: the
    share of each power z^k (k = 0, −1, −2, −3, −4) in (iρ/2)∮(dw/dz)² dz — only k = −1 is nonzero. Title "Any contour,
    any body: only U × Γ/z survives". *see:* "flat lines; one bar." *read:* "Cauchy lets the contour move; the residue
    picks the 1/z term; its coefficient is iUΓ/π whatever the body's doublet d." *change:* "…a source hidden in the body
    (an open body, m ≠ 0): the 1/z coefficient gains U m/π and the drag becomes D = −ρUm — a thrust, pushing upstream
    (held-source force, Exercise-style)."
23. `nb.explainer("blasius_kutta_contour", heading="Why doesn't lift depend on the shape?", why="A proof cannot show
    invariance under deformation; here you drag the contour — radius, offset, shape — and watch the integral stay pinned
    while the integrand changes completely; the Laurent bars show which single term survives.", tries=["Drag the contour
    outward: the lift readout stays at ρUΓ.", "Switch the body to the ellipse with the same Γ: the same lift.", "Preset
    'contour cuts the body': the status turns ⚠️ and the number jumps — Cauchy needs analyticity in between.", "Preset
    'few quadrature points': see the trapezoid converge exponentially."])`
24. `nb.md` — **What would change if…** "…the flow were three-dimensional (a finite wing)? The 2-D proof no longer
    applies: the trailing vortices of a finite wing induce a downwash and a drag even in ideal flow — induced drag,
    Ch. 14. In 2-D there is no such escape: D = 0 always."

### A.6 §6.6 Conformal Mapping — C11 (+N64–N70 · D19 D20 D21 · E6)
1. `nb.section("6.6", "Conformal Mapping", intro="**What is this section about?** An analytic function maps one plane to
   another, turning and stretching every small figure but keeping its angles. A flow net — streamlines crossing
   equipotentials at right angles — therefore maps to another flow net. Solve the easy problem, flow round a circle, and
   let a map carry it: the Zhukhovsky map turns the circle into an ellipse (and, in Ch. 14, an airfoil).")`

**C11 — Conformal mapping and the Zhukhovsky transformation (6.63)–(6.69)**
2. `nb.core("C11", "Conformal maps keep angles, $\\delta w=\\frac{dw}{dz}\\delta z$ (6.63), and the Zhukhovsky map
   $z=\\zeta+b^2/\\zeta$ (6.65) carries the circle's flow onto an ellipse", question="We can solve the flow round a circle.
   How can that one solution give the flow round other shapes?")`
3. `nb.md` — **The problem in plain words:** "Draw a circle on a rubber sheet with a fine square grid, and stretch the
   sheet smoothly without tearing: the circle becomes an oval, but where you look closely the little squares are still
   little squares, just turned and resized. A conformal map is such a stretch done by an analytic function. Because
   streamlines and equipotentials form a grid of little squares, the map carries a flow into a flow."
4. `nb.md` — **The idea** (ASCII):
   ```
   ζ-plane: circle |ζ| = a, flow (6.68) known ──z = ζ + b²/ζ──► z-plane: ellipse (6.67), flow = the same w at ζ(z)
   small cross at ζ₀ ──× f′(ζ₀) = |f′|e^{iχ}──► turned by χ, stretched by |f′|, right angle kept
   inverse: two roots ζ, b²/ζ — keep the one OUTSIDE the circle
   ```
5. `nb.md` — 🔁 reminder: multiplying complex numbers adds their angles (primer "the complex plane in numpy", C09).
6. `nb.note` — **N64 [B]** "A small step δz at z₀ becomes" equation `\delta w=\frac{dw}{dz}\delta z`, ref "6.63".
7. `nb.derivation("D19", …)` — Part F D19 (6 steps), ref "6.64".
8. `nb.note` — **N65 [B]** "A second step δ′z at the same point becomes $\delta'w=\frac{dw}{dz}\delta'z$ *(Eq. 6.64)*, turned by
   the same angle arg(dw/dz), so the angle between the two is kept (α = β). It fails where dw/dz = 0 or ∞, and only
   *small* figures keep their shape." + `nb.code`: `f, fp = (lambda z: z**2), (lambda z: 2*z)` · `print(np.degrees(
   cm.angle_preservation(f, fp, 1+0.5j, 1e-3, 1e-3j)))` · `print(np.degrees(cm.angle_preservation(f, fp, 0j, 1e-3,
   1e-3*np.exp(1j*np.pi/4))))`. *expect:* `(90.0, 90.0)` · `(45.0, 90.0)` (doubled at the critical point z = 0).
9. `nb.note` — **N66 [B]** "**Flow nets as images of a grid.** A rectangular grid of φ = const and ψ = const lines in the
   w-plane is uniform flow; its image in the z-plane is the flow net of w = f(z). Chains of maps work too: w = ln ζ with
   ζ = sin z gives $u-iv=\frac{dw}{d\zeta}\frac{d\zeta}{dz}=\frac1\zeta\cos z=\cot z$; w = z² gives the 90° corner (the
   hyperbolae of (6.24) — ⚠️ the book's '(see Figure 6.5)' means Figure 6.3)." + `nb.code`: `print(ch06.cot_flow(0.7+0.2j),
   1/np.tan(0.7+0.2j))` → *expect:* equal.
10. `nb.derivation("D20", …)` — Part F D20 (8 steps), ref "6.67".
11. `nb.note` — **N67 [B], N68 [B]** "The circle of radius a > b maps to" equation `z=ae^{i\theta}+\frac{b^2}{a}e^{-i\theta}`,
    ref "6.66", "an ellipse" equation `\frac{x^2}{(a+b^2/a)^2}+\frac{y^2}{(a-b^2/a)^2}=1`, ref "6.67", "with foci at ±2b.
    b = 1 m, a = 1.2 m: semi-axes 2.033 m and 0.367 m, foci ±2 m." + `nb.code`: `print(ch06.joukowski_ellipse(1.2, 1.0))`
    → *expect:* `{'A': 2.0333, 'B': 0.3667, 'foci': 2.0}`.
12. `nb.note` — **N69 [B]** "In the ζ-plane, the circle of radius a in a stream U with a clockwise circulation Γ — (6.52)
    with z replaced by ζ:" equation `w=U\Big(\zeta+\frac{a^2}{\zeta}\Big)+\frac{i\Gamma}{2\pi}\ln(\zeta/a)`, ref "6.68".
13. `nb.primer("complex square roots and the quadratic formula", "A quadratic $\\zeta^2-z\\zeta+b^2=0$ has two roots
    $\\zeta=\\tfrac12[z\\pm(z^2-4b^2)^{1/2}]$, and in the complex plane 'the' square root is a choice: numpy's `np.sqrt`
    returns the principal root (real part ≥ 0), whose cut lies where the argument is a negative real number. The two
    roots here multiply to b² (Vieta, P71), so one is outside the circle \\|ζ\\| = b and one inside — and the principal
    root does **not** always pick the outside one.", code="b, z = 1.0, -3 + 0.5j              # a point to the left of the
    slit\nr = np.sqrt(z**2 - 4*b**2)                # numpy's principal root\nprint(abs(0.5*(z + r)), abs(0.5*(z - r)))
    # 0.370 and 2.701: '+' gives the INSIDE root here\nprint(abs(0.5*(z + np.sqrt(z - 2*b)*np.sqrt(z + 2*b))))   # 2.701:
    the product form picks the outside root")` — *expect:* `0.3702 2.7013` · `2.7013`.
14. `nb.derivation("D21", …, check_src=…)` — Part F D21 (11 steps, ★★★), ref "6.69".
15. `nb.note` — **N70 [B]** "**The inverse and the velocity:**" equation `\zeta=\tfrac12z+\tfrac12(z^2-4b^2)^{1/2}`, ref
    "6.69", "with the root that falls outside the circle, and $u-iv=\frac{dw}{dz}=\frac{dw}{d\zeta}\frac{d\zeta}{dz}$. ⚠️ Coded as
    $\zeta=\tfrac12[z+\sqrt{z-2b}\sqrt{z+2b}]$: at z = −3 + 0.5i (b = 1) it gives \|ζ\| = 2.70, while numpy's principal
    √(z² − 4b²) gives 0.37 — a point inside the cylinder, and a flow that makes no sense."
16. `nb.worked_example("an ellipse from a circle", "b = 1 m, a = 1.2 m. 1. The circle ζ = 1.2e^{iθ} maps to x = (1.2 +
    1/1.2) cos θ = 2.033 cos θ, y = (1.2 − 0.833) sin θ = 0.367 sin θ. 2. Foci: c² = 2.033² − 0.367² = 4.134 − 0.134 = 4
    → c = 2 = 2b. 3. At θ = 0 the map's derivative dz/dζ = 1 − b²/ζ² = 1 − 1/1.44 = 0.306: lengths along the circle
    near its right end shrink to 31 % — the ellipse's sharp ends. 4. Inverse at z = −3 + 0.5i: the product form gives
    ζ = −2.638 + 0.579i (\|ζ\| = 2.70 > 1.2, outside ✓).")`
17. `nb.code` — *code:* `b, a = 1.0, 1.2` · `th = np.linspace(0, 2*np.pi, 400, endpoint=False)` · `zc =
    cm.joukowski(a*np.exp(1j*th), b)` · `print(np.max(np.abs(zc.real**2/2.0333333**2 + zc.imag**2/0.3666667**2 - 1)))`
    · `zz = np.array([-3+0.5j, -3-0.5j, 3+0.5j, 0.1+2j])` · `print(np.abs(cm.joukowski_inverse(zz, b)),
    np.abs(cm.joukowski_inverse(zz, b, branch="principal")))` · `e = ch06.elliptic_cylinder_flow(1.0, a, b,
    Gamma_cw=1.0)` · `print(pf.normal_velocity_on(e, zc), pf.far_field_check(e, 200.0))`. *expect:* ≤ 1e-6 (the ellipse
    equation, 7-digit semi-axes) · outside: all \|ζ\| > 1 (2.70, 2.70, 2.70, …); principal: 0.37 and 0.37 for the two
    left points · ≤ 1e-12 · ~Γ/2πR ≈ 8e-4. *explain:* 1. the circle's image satisfies (6.67); 2. the two inverses
    compared in all four quadrants — the principal root fails on the left; 3. the mapped flow has no flow through the
    ellipse and the right stream far away.
18. `nb.check_agree` — **from scratch (curation §7):** `r1 = 0.5*(zz + np.sqrt(zz**2 - 4*b**2))`; `r2 = b**2/r1` (the
    other root: the product is b²); `mine = np.where(np.abs(r1) >= b, r1, r2)` (keep the root outside the circle); `assert
    np.allclose(mine, cm.joukowski_inverse(zz, b))`. Markdown: "Picking the root with \|ζ\| ≥ b by hand gives the
    library's branch everywhere."
19. `nb.figure` — **two planes** (2 panels, 10 × 3.6 in): left, the ζ-plane: the circle a = 1.2 (black), the critical
    points ±b (×), streamlines of (6.68) with Γ = 1 m²/s (U = 1 m/s); right, the z-plane: the ellipse, the slit
    [−2b, 2b] dashed inside, the mapped streamlines (the same ψ levels). Title "Map the circle, carry the flow". *see:*
    "the same streamlines, bent from a circle onto an ellipse." *read:* "ψ is the same number at a ζ and at its image z,
    so a streamline maps to a streamline; angles between streamlines and the body are kept." *change:* "…a → b: the
    ellipse flattens to the plate [−2, 2]; the speed at its ends blows up unless Γ is chosen to cancel it (the Kutta
    condition, Ch. 14)."
20. `nb.figure` — **grid images** (3 panels): the image of a (φ, ψ) grid under w = z² (hyperbolae), w = e^z (a polar net)
    and w = ln sin z (the cot flow), via `cm.grid_image`. *see / read / change* in the same pattern: "right angles
    everywhere except at the critical points (z = 0 for z²)."
21. `nb.plotly` — `slider_figure` over a/b from 1.0 to 3.0 (21 steps, FAST 11): the z-plane body and five mapped
    streamlines (Γ = 0), the title printing the semi-axes and the foci. Title "One map, a family of bodies: plate →
    ellipse → almost a circle".
22. `nb.explainer("conformal_joukowski", heading="How does a map carry a flow?", why="Two linked planes on one state is
    the natural interactive form: drag a small cross in the ζ-plane and watch its image stay a right angle (and fail at
    ζ = ±b), slide a/b from a plate to a circle, and toggle numpy's square root to see half the streamlines jump inside
    the body.", tries=["Drag the cross toward ζ = +b: the image angle opens toward 180°.", "Slide a/b to 1: the ellipse
    becomes a plate.", "Switch the branch chips to 'numpy principal': streamlines on the left jump inside — the status
    turns ⚠️.", "Mode 'z²' at 0: angles doubled."])`
23. `nb.md` — **What would change if…** "…the circle were shifted off-centre (up and a little left) before mapping? The
    image would be a cambered, round-nosed, sharp-tailed shape — a Zhukhovsky airfoil. Ch. 14 does exactly this, and
    fixes Γ so the flow leaves the sharp tail smoothly."

### A.7 §6.7 Numerical Solution Techniques in Two Dimensions — R19 R20 R21, C12 (+N71–N74 · D22 D23 · E7), pointer N90
1. `nb.section("6.7", "Numerical Solution Techniques in Two Dimensions", intro="**What is this section about?** Exact
   solutions exist only for simple shapes. On a grid, Laplace's equation becomes a rule — every value is the average of
   its four neighbours — and a sweep that keeps enforcing the rule converges to the solution. This is the book's first
   numerical PDE solver, and the ancestor of every elliptic solver in Ch. 10 and of bounded-domain streamfunction
   inversion in ocean and atmosphere models.")`
2. `nb.recap("R19", "Grids and half-point differences", "Store ψ at grid points, $\\psi_{i,j}=\\psi(i\\Delta x,j\\Delta y)$, and
   approximate a first derivative by the difference across a cell centred on the point, $(\\partial\\psi/\\partial x)_{i,j}\\simeq
   (\\psi_{i+\\frac12,j}-\\psi_{i-\\frac12,j})/\\Delta x$ (Ch. 2 §2.9 did central differences). ⚠️ The book calls these
   'first-order' differences: they approximate a *first derivative* and are **second-order accurate** (error ∝ Δx², D22
   step 3).", where="Ch. 2 §2.9")`
3. `nb.recap("R20", "Second difference in x", "Differencing twice: $(\\partial^2\\psi/\\partial x^2)_{i,j}\\simeq(\\psi_{i+1,j}-2\\psi_{i,j}+
   \\psi_{i-1,j})/\\Delta x^2$ *(Eq. 6.70)*.", where="Ch. 2 §2.9")`
4. `nb.recap("R21", "Second difference in y", "$(\\partial^2\\psi/\\partial y^2)_{i,j}\\simeq(\\psi_{i,j+1}-2\\psi_{i,j}+\\psi_{i,j-1})/
   \\Delta y^2$ *(Eq. 6.71)*.", where="Ch. 2 §2.9")` + `nb.code`: `hs = np.array([0.1, 0.05, 0.025])`; `err = []` · `for h in hs: x = np.arange(0.0, 1.0 + h/2, h);
   X, Y = np.meshgrid(x, x); psi = np.sin(np.pi*X)*np.sinh(np.pi*Y); lap = ls.laplacian_5pt(psi, h, h);
   err.append(np.max(np.abs(lap[1:-1, 1:-1])))` (ψ is harmonic, so the exact Laplacian is 0 and lap is pure error) ·
   `print(err, observed_order(hs, err))`. *expect:* errors falling ×4 per halving, order 2.00 ± 0.05.

**C12 — Finite-difference Laplace: the average rule and Gauss–Seidel relaxation (6.72)**
5. `nb.core("C12", "On a grid, Laplace says 'be the average of your neighbours': $\\psi_{i,j}=\\tfrac14[\\psi_{i-1,j}+\\psi_{i+1,j}+
   \\psi_{i,j-1}+\\psi_{i,j+1}]$ (6.72)", question="How does a grid of numbers solve Laplace's equation — and how do we know
   when it has?")`
6. `nb.md` — **The problem in plain words:** "Water flows through a channel that narrows abruptly. No formula gives
   its streamlines — but ψ must be constant on each wall and harmonic in between. Stretch a rubber membrane over a wire
   frame bent to the boundary values: every point settles at the average height of its neighbours. A computer can do
   the settling one point at a time."
7. `nb.md` — **The idea** (ASCII):
   ```
             ψ(i, j+1)
   ψ(i−1, j)   ψ(i, j)   ψ(i+1, j)       ψ(i, j) = ¼ (sum of the four)          ← (6.72)
             ψ(i, j−1)                   sweep: replace every value by the average of its neighbours, repeat
   Jacobi: use last sweep's values · Gauss–Seidel: use the newest · SOR: overshoot by ω
   ```
8. `nb.md` — 🔁 reminders: Taylor series (P26), `np.linalg.solve` (P57), eigenvalues (P80), log axes (P13).
9. `nb.note` — **N71 [C]** "**Other routes for complicated shapes:** distributions of sources and sinks (C14, and N90's
   panels), thin-body perturbation theory (Ch. 14), and numerical solution of Laplace's equation (here; Ch. 10)."
10. `nb.derivation("D22", …)` — Part F D22 (6 steps), ref "6.72".
11. `nb.worked_example("four unknowns with ψ = xy on the boundary", "Put the 4 × 4 grid of Fig. 6.23 on points x, y ∈ {0,
    1, 2, 3} (node (i, j) at x = i − 1, y = j − 1) and give the 12 boundary nodes the values of ψ = xy (a harmonic
    polynomial). 1. ψ₂₂ = ¼(ψ₁₂ + ψ₃₂ + ψ₂₁ + ψ₂₃) = ¼(0 + 2 + 0 + 2) = 1 = 1·1 ✓. 2. ψ₃₂ = ¼(ψ₂₂ + ψ₄₂ + ψ₃₁ + ψ₃₃) =
    ¼(1 + 3 + 0 + 4) = 2 ✓. 3. Likewise ψ₂₃ = 2, ψ₃₃ = 4. xy satisfies the average rule **exactly** (its fourth derivatives
    vanish, D22 step 3).")`
12. `nb.note` — **N72 [B]** "**The 16-point grid:** 12 known boundary values ψ^B and four unknowns give four linear
    equations —" equation `\begin{aligned}\psi_{2,2}&=\tfrac14[\psi^B_{1,2}+\psi_{3,2}+\psi^B_{2,1}+\psi_{2,3}],\\ \psi_{3,2}&=\tfrac14[\psi_{2,2}+\psi^B_{4,2}+\psi^B_{3,1}+\psi_{3,3}],\\ \psi_{2,3}&=\tfrac14[\psi^B_{1,3}+\psi_{3,3}+\psi_{2,2}+\psi^B_{2,4}],\\ \psi_{3,3}&=\tfrac14[\psi_{2,3}+\psi^B_{4,3}+\psi_{3,2}+\psi^B_{3,4}]\end{aligned}`,
    ref "6.73". + `nb.code`: `sysm = ch06.four_point_system("xy")` · `print(sysm["A"], sysm["b"])` · `print(
    np.linalg.solve(sysm["A"], sysm["b"]))`. *expect:* A = [[4, −1, −1, 0], [−1, 4, 0, −1], [−1, 0, 4, −1], [0, −1, −1,
    4]], b = [0, 3, 3, 12] · `[1. 2. 2. 4.]`.
13. `nb.primer("iterative solvers: Jacobi, Gauss–Seidel, SOR", "A big linear system Aψ = b from a stencil is solved by
    repeated sweeps instead of elimination. **Jacobi** updates every unknown from the *previous* sweep's values;
    **Gauss–Seidel** uses each new value as soon as it exists; **SOR** takes the Gauss–Seidel change and multiplies it by
    ω (1 < ω < 2). The error shrinks each sweep by a factor ρ < 1 (the spectral radius of the sweep's matrix): the closer
    ρ is to 1, the slower. Stop when the **residual** b − Aψ is small, not just when the change is.", code="A =
    np.array([[4., -1], [-1, 4]]); b = np.array([3., 3])   # a 2 × 2 toy system, answer (1, 1)\nx = np.zeros(2)\nfor k in
    range(4):\n    x[0] = (b[0] + x[1])/4; x[1] = (b[1] + x[0])/4   # Gauss–Seidel: x[1] uses the new x[0]\n    print(k,
    x, np.max(np.abs(b - A @ x)))            # residual falls ×16 per sweep")` — *expect:* residual 0.9375 → 0.0586 →
    0.00366 → 2.3e-4 (÷16 per sweep).
14. `nb.derivation("D23", …)` — Part F D23 (9 steps), ref "6.73".
15. `nb.note` — **N73 [B]** "**Gauss–Seidel:** sweep through the nodes, always using the latest value available:"
    equation `\psi^{(k+1)}_{i,j}=\tfrac14\big[\psi^{(k+1)}_{i-1,j}+\psi^{(k)}_{i+1,j}+\psi^{(k+1)}_{i,j-1}+\psi^{(k)}_{i,j+1}\big]`
    (no number). "The book stops when values stop changing; we stop on the residual of (6.72), which can be large while
    the change is small (D23 step 8)." + `nb.code`: `for method in ("jacobi", "gauss_seidel", "sor"): psi, hist =
    ls.solve_laplace(*ch06.four_point_system("xy", as_grid=True), method=method, tol=1e-10); print(method,
    hist["sweeps"], psi[1:3, 1:3].ravel())` (the builder uses the grid form of the 4-point problem: `four_point_system(...,
    as_grid=True)` → `(mask, bc)`; Part C 5.31 gets this keyword). *expect:* Jacobi 34, Gauss–Seidel 19, SOR 12 sweeps to a
    residual of 1e-10 (ω_opt = 1.072 for this grid); the four values 1, 2, 2, 4.
16. `nb.primer("boolean masks and scipy.sparse", "An L-shaped or stepped domain lives in a rectangular array: a boolean
    array `mask[j, i]` is True at the unknown nodes, and the boundary values sit in another array. For big grids we skip
    iteration and hand the whole system to a sparse direct solver: `scipy.sparse` stores only the non-zero entries (five
    per row here) and `scipy.sparse.linalg.spsolve` solves it.", code="import scipy.sparse as sps, scipy.sparse.linalg as
    spl\nmask = np.ones((3, 4), bool); mask[0, :2] = False       # a small stepped domain: True = unknown\nprint(mask.sum(),
    'unknowns')                              # 10\nA = sps.diags([4, -1, -1], [0, 1, -1], shape=(4, 4), format='csr')   #
    a 1-D toy stencil\nprint(spl.spsolve(A, np.ones(4)))                  # [0.3636 0.4545 0.4545 0.3636]")` —
    *expect:* `10 unknowns` · `[0.36363636 0.45454545 0.45454545 0.36363636]`.
17. `nb.note` — **N74 [B]** "**Example 6.2: a sharp contraction.** A channel narrows in one step. ψ = 0 on the lower wall
    and on the step, ψ = Q on the upper wall; across the inlet and the outlet the velocity is uniform, so ψ rises
    linearly there (u = ∂ψ/∂y constant). The region has no holes, so the solution is unique. At the step's corner the
    fluid turns through 270°: by (6.46) with α = 3π/2, the speed there is infinite (∝ r^{−1/3}), which slows grid
    convergence near it. We run it with Q = 1 m²/s and show ψ/Q (the book's flow rate and printed grid values stay in
    the private test file; our Gauss–Seidel reproduces them). ⚠️ The book's FORTRAN listing loops over I where it should
    loop over J (harmless there), and prints the inlet Δψ in m² — it is m²/s."
18. `nb.code` — *code:* `ex = ch06.example_6_2(Q=1.0)` · `print(ex["grid_shape"], ex["history"]["sweeps"],
    ex["history"]["residual"][-1])` · `for r in (1, 2, 4): e = ch06.example_6_2(Q=1.0, refine=r, method="direct" if r > 1
    else "gauss_seidel"); print(r, e["psi"][e["probe"]])` (`probe` = the index of a fixed interior point, the same physical
    point at every refinement — Part C 5.32). *expect:* the grid shape of the book's figure, the sweep count to
    1e-10, ψ/Q at the probe converging with refinement (differences falling ≈ ×4 away from the corner). *explain:* 1.
    the contraction on the book's grid, relaxed by Gauss–Seidel; 2. the same problem on grids 2× and 4× finer, solved
    directly; 3. the probe value settles: the solution converges.
19. `nb.check_agree` — **from scratch (curation §7):** a double-loop Gauss–Seidel on the 4-point problem: `P =
    np.array([[(i)*(j) for i in range(4)] for j in range(4)], float)` (ψ = xy on the 4 × 4 grid, row j = y, column i = x);
    `P[1:3, 1:3] = 0.0` (unknowns start at zero); `for sweep in range(40): for j in (1, 2): for i in (1, 2): P[j, i] =
    0.25*(P[j, i-1] + P[j, i+1] + P[j-1, i] + P[j+1, i])` (newest values used at once); `assert np.allclose([P[1, 1], P[1, 2], P[2, 1], P[2, 2]],
    np.linalg.solve(sysm["A"], sysm["b"]))` (ψ₂₂, ψ₃₂, ψ₂₃, ψ₃₃ are P[1,1], P[1,2], P[2,1], P[2,2] with rows = y,
    columns = x). Markdown: "Four lines of loops do what `solve_laplace` does."
20. `nb.figure` — two panels (9 × 3.4 in): left, the contraction with ψ/Q contours (0.1 … 0.9) from `example_6_2(refine=4,
    method="direct")`, the walls black, the re-entrant corner circled; right, residual (solid) and max change (dashed)
    against the sweep number on a log axis for Jacobi, Gauss–Seidel and SOR (ω_opt) on the book grid. Title "Streamlines
    by averaging; the smarter sweep wins". *see:* "streamlines squeezing into the narrow part; three straight lines on
    the log plot with different slopes." *read:* "a straight line on a log axis means a constant factor per sweep, the
    spectral radius; Gauss–Seidel's slope is twice Jacobi's (ρ_GS = ρ_J²); SOR is steeper still; the dashed change sits
    below the residual — stopping on it would stop too early." *change:* "…the grid refined ×2: Jacobi and Gauss–Seidel
    need ≈ 4× the sweeps (∝ N²), SOR ≈ 2× (∝ N)."
21. `nb.animation` — **Gauss–Seidel relaxing** (`player="frames"`, 24 frames, FAST 12): ψ on the book grid from zero,
    one frame per sweep for the first 12 sweeps then every 4th, colour + numbers on the nodes, the residual curve
    growing beside it. + see/read/change markdown: "the boundary values leak inward sweep by sweep; the residual line
    falls steadily."
22. `nb.pointer("§6.7 also names source panels — a boundary-only numerical method; our constant-strength source panels
    (N90) are taught with the axial singularity method in C14 (§6.8).")`
23. `nb.explainer("laplace_relaxation", heading="How does a grid solve Laplace's equation?", why="The method is a
    process: stepping one node at a time, then whole sweeps, with the residual curve growing beside the grid, makes the
    average rule and the difference between Jacobi, Gauss–Seidel and SOR visible.", tries=["Step node by node with
    Gauss–Seidel on the 4-point problem: the first values are 0, 0.75, 0.75, 3.375.", "Switch to Jacobi and compare
    the residual slopes.", "Preset 'SOR at ω_opt' on the refined contraction: an order of magnitude fewer sweeps.",
    "Preset 'ψ = xy': one solve is exact — the harmonic polynomial satisfies (6.72) exactly."])`
24. `nb.md` — **What would change if…** "…the right-hand side were not zero — the vorticity of a rotating flow, ∇²ψ = −ω
    (6.4)? The same sweeps solve Poisson's equation with ω h²/4 added to each average (`ls.solve_poisson`). That is how
    Ch. 10's vorticity–stream-function method, and ocean models that invert potential vorticity in a closed basin,
    recover ψ every time step."

### A.8 §6.8 Axisymmetric Ideal Flow — R22–R30, C13 (+N75–N85 · D24 D25), C14 (+N86–N90 · D26 D27 · E8)
1. `nb.section("6.8", "Axisymmetric Ideal Flow", intro="**What is this section about?** Flow round a sphere, a bullet or
   an airship is the same in every plane through the axis. One stream function — Stokes' — still describes it, but its
   equation is *not* Laplace's, so complex variables no longer help; elements are built with φ and ψ in spherical
   coordinates instead. ⚠️ In this section the symmetry axis z is horizontal, along the stream.")`
2. `nb.recap("R22", "Two stream functions in 3-D", "A general 3-D incompressible flow needs two stream functions (u =
   ∇χ × ∇ψ); an axisymmetric one needs only one.", where="Ch. 4 §4.3")`
3. `nb.recap("R23", "Axisymmetric continuity", "In cylindrical coordinates (R, φ, z) with no swirl and no φ-dependence:
   $\\frac1R\\frac{\\partial}{\\partial R}(Ru_R)+\\frac{\\partial u_z}{\\partial z}=0$ *(Eq. 6.74)*. ⚠️ z is now horizontal.",
   where="Ch. 4 §4.3")`
4. `nb.recap("R24", "The Stokes stream function", "Choosing χ = −φ (the azimuth angle) in u = ∇χ × ∇ψ gives
   $u_R=-\\frac1R\\frac{\\partial\\psi}{\\partial z}$, $u_z=\\frac1R\\frac{\\partial\\psi}{\\partial R}$ *(Eq. 6.75)*, which satisfies (6.74)
   identically. Same signs as Ch. 4's `velocity_from_streamfunction_axisym`.", where="Ch. 4 §4.3")` + `nb.code`: `from
   fluidpy.core import streamfunction as sf` · `print(sf.velocity_from_streamfunction_axisym(lambda R, z: 0.5*2.0*R**2,
   0.3, 0.1))` → *expect:* `(0.0, 2.0)` (ψ = ½UR² with U = 2 m/s is the uniform stream along z).
5. `nb.recap("R25", "The azimuthal vorticity", "$\\omega_\\varphi=\\frac{\\partial u_R}{\\partial z}-\\frac{\\partial u_z}{\\partial
   R}$ *(Eq. 6.76)* — Appendix B's curl for an axisymmetric field.", where="Ch. 4 (App. B operators)")`
6. `nb.recap("R26", "The axisymmetric potential", "$u_R=\\partial\\phi/\\partial R$, $u_z=\\partial\\phi/\\partial z$ *(Eq. 6.79)* — u =
   ∇φ (3.17) in cylindrical coordinates.", where="Ch. 3 §3.4")`
7. `nb.recap("R27", "Cylindrical and spherical coordinates", "(R, φ, z) and (r, θ, φ) with z along the axis: R = r sin θ,
   z = r cos θ, θ measured from +z (downstream) — the table the book gives in (6.81).", where="Ch. 3 §3.1")`
8. `nb.recap("R28", "Spherical continuity", "$\\frac1{r^2}\\frac{\\partial}{\\partial r}(r^2u_r)+\\frac1{r\\sin\\theta}\\frac{\\partial}{\\partial
   \\theta}(u_\\theta\\sin\\theta)=0$. ⚠️ The book's (6.82) prints $\\frac1r\\frac{\\partial}{\\partial r}(r^2u_r)+\\frac1{\\sin\\theta}\\frac{\\partial}
   {\\partial\\theta}(u_\\theta\\sin\\theta)=0$ — the same equation multiplied by r except for the first factor, which should be
   1/r² (Appendix B); '= 0' is unaffected, a residual would not be. Our code uses Appendix B.", where="Ch. 4 (App. B)")`
9. `nb.recap("R29", "Spherical vorticity", "$\\omega_\\varphi=\\frac1r\\Big[\\frac{\\partial}{\\partial r}(ru_\\theta)-\\frac{\\partial
   u_r}{\\partial\\theta}\\Big]$ *(Eq. 6.84)*.", where="Ch. 4 (App. B)")`
10. `nb.recap("R30", "Spherical axisymmetric Laplace", "$\\frac1{r^2}\\frac{\\partial}{\\partial r}\\Big(r^2\\frac{\\partial\\phi}{\\partial
    r}\\Big)+\\frac1{r^2\\sin\\theta}\\frac{\\partial}{\\partial\\theta}\\Big(\\sin\\theta\\frac{\\partial\\phi}{\\partial\\theta}\\Big)=0$ *(Eq.
    6.85)* — used in D25's check.", where="Ch. 4 (App. B)")` + `nb.code` (checks R28–R30 on the sphere later): `from
    fluidpy.core import curvilinear as cu` · `r_, t_, ph_ = cu.coordinates("spherical")` (the module's own symbols) · `print(
    cu.laplacian(r_*sp.cos(t_), "spherical"), cu.laplacian(sp.cos(t_)/r_**2, "spherical"))` → *expect:* `0 0` (the uniform stream and the 3-D doublet potentials are
    harmonic).

**C13 — Axisymmetric potential flow: the Stokes stream function, 3-D elements and the sphere (6.90)**
11. `nb.core("C13", "The sphere in a stream, $u_\\theta=-U\\big[1+\\tfrac12(a/r)^3\\big]\\sin\\theta$ (6.90), from a stream function
    whose equation (6.77) is not Laplace's", question="Air flows round a football and round a lamp-post. What changes when
    the flow can also go *over* the body?")`
12. `nb.md` — **The problem in plain words:** "A lamp-post is effectively 2-D: the air can only go round the sides. A
    football is 3-D: the air can go round, over and under, every way at once. That extra freedom should make the
    disturbance smaller and the fastest speed lower. We want the numbers — and we find along the way that the tidy
    complex-variable tools of §6.4–§6.6 do not carry over."
13. `nb.md` — **The idea** (table): | | cylinder (2-D) | sphere (3-D) | — | built from | stream + 2-D doublet | stream +
    3-D doublet | — | fastest surface speed | 2U | 1.5U | — | minimum C_p | −3 | −1.25 | — | disturbance decays like |
    (a/r)² | (a/r)³ | — | stream-function equation | Laplace | (6.77), not Laplace |. "**The fluid escapes sideways, so it
    needs to speed up less.**"
14. `nb.md` — 🔁 reminders: cylindrical and spherical unit vectors (P88), cross product (Ch. 2 §2.4) + gloss
    "**cylindrical cross products:** e_R × e_φ = e_z, e_φ × e_z = e_R, e_z × e_R = e_φ (right-handed, cyclic), so
    e_φ × e_R = −e_z", Taylor series (P98), 3-D point source flux (gloss "a source of Q m³/s sends Q through every sphere
    round it: 4πr²u_r = Q"), Poisson/Green in 3-D (P139).
15. `nb.derivation("D24", …)` — Part F D24 (8 steps), ref "6.77".
16. `nb.note` — **N75 [B]** "**The field equation of the Stokes stream function** in irrotational flow:" equation
    `\frac{\partial}{\partial R}\Big(\frac1R\frac{\partial\psi}{\partial R}\Big)+\frac1R\frac{\partial^2\psi}{\partial z^2}=-\omega_\varphi=0`,
    ref "6.77". "It is not the Laplacian of ψ, so ψ is not the imaginary part of an analytic function: no complex
    potential in axisymmetric flow." + `nb.code`: `sph = pf.sphere(1.0, 1.0)` · `R = np.array([0.5, 2.0, 1.5]); Z =
    np.array([2.0, -1.0, 0.8])` · `print(ch06.stokes_operator_residual(sph.psi, R, Z))` · `print(
    ch06.stokes_operator_residual(lambda R, z: R*z, R, Z))` (a control). *expect:* ≲ 1e-8 · clearly non-zero (the
    control fails the equation).
17. `nb.note` — **N76 [C]** "**Units.** The Stokes ψ is in m³/s (a volume flow), the plane ψ in m²/s (per unit depth);
    ψ = const is a surface of revolution."
18. `nb.note` — **N77 [B]** "**Flow between two stream surfaces:**" equation
    `dQ=2\pi R(\mathbf u\cdot\mathbf n)ds=2\pi R(-u_Rdz+u_zdR)=2\pi\,d\psi`, ref "6.78". "⚠️ 2π, not 1: dψ is the flow per
    radian round the axis." + `nb.code`: `print(ch06.axisym_flux_between(sph.psi, (1.5, 0.0), (3.0, 0.5)))` → *expect:* two
    equal numbers (the quadrature of 2πR u·n ds and 2π(ψ₂ − ψ₁)).
19. `nb.note` — **N78 [B]** "The potential obeys the true (axisymmetric) Laplace equation:" equation
    `\frac1R\frac{\partial}{\partial R}\Big(R\frac{\partial\phi}{\partial R}\Big)+\frac{\partial^2\phi}{\partial z^2}=0`, ref "6.80".
20. `nb.note` — **N79 [B]** "In spherical coordinates:" equation `u_r=\frac1{r^2\sin\theta}\frac{\partial\psi}{\partial\theta}=\frac{\partial\phi}{\partial r},\qquad u_\theta=-\frac1{r\sin\theta}\frac{\partial\psi}{\partial r}=\frac1r\frac{\partial\phi}{\partial\theta}`,
    ref "6.83". + `nb.code`: `print(pf.axisym_velocity_spherical_sym(sp.Rational(1, 2)*r_**2*sp.sin(t_)**2, r_, t_))` →
    *expect:* `(cos(theta), -sin(theta))` (the unit stream U = 1).
21. `nb.derivation("D25", …)` — Part F D25 (10 steps), ref "6.91".
22. `nb.note` — **N80 [B], N81 [B], N82 [B]** "**The 3-D elements** (a table; each equation with its number):" uniform
    flow along z: `\phi=Uz,\ \psi=\tfrac12UR^2;\qquad \phi=Ur\cos\theta,\ \psi=\tfrac12Ur^2\sin^2\theta` (6.86) · point source Q
    [m³/s]: `\phi=-\frac{Q}{4\pi\sqrt{R^2+z^2}},\ \psi=-\frac{Qz}{4\pi\sqrt{R^2+z^2}};\qquad \phi=-\frac{Q}{4\pi r},\ \psi=-\frac{Q}{4\pi}\cos\theta` (6.87)
    · doublet with dipole −d e_z: `\phi=\frac{d}{4\pi}\frac{z}{(R^2+z^2)^{3/2}},\ \psi=-\frac{d}{4\pi}\frac{R^2}{(R^2+z^2)^{3/2}};\qquad \phi=\frac{d}{4\pi r^2}\cos\theta,\ \psi=-\frac{d}{4\pi r}\sin^2\theta` (6.88).
    "In a plane through the axis their streamlines look like the 2-D ones."
23. `nb.note` — **N83 [B]** "**The sphere** = stream + an opposing doublet of strength d = 2πa³U:" equation
    `\psi=\tfrac12Ur^2\sin^2\theta-\frac{d}{4\pi r}\sin^2\theta=\tfrac12Ur^2\Big(1-\frac{a^3}{r^3}\Big)\sin^2\theta;\quad \phi=Ur\Big(1+\frac{a^3}{2r^3}\Big)\cos\theta`,
    ref "6.89". "ψ = 0 on the whole axis and on r = a. Outside, it is the flow round Hill's vortex seen from the vortex
    (Ch. 5 §5.4)." + `nb.code`: `print(sph.psi(1.3, 0.4), 0.5*1.3**2*(1 - 1/np.hypot(1.3, 0.4)**3))` (U = a = 1: ψ =
    ½R²(1 − a³/r³) since r² sin²θ = R²) → *expect:* two equal numbers, 0.4935.
24. `nb.note` — **N84 [B]** "**Surface speed and pressure:** with (6.90) $u_r=U[1-(a/r)^3]\cos\theta$,
    $u_\theta=-U[1+\frac12(a/r)^3]\sin\theta$ the surface speed is (3/2)U sin θ, and" equation
    `C_p=\frac{p-p_\infty}{\frac12\rho U^2}=1-\Big(\frac{u_\theta}{U}\Big)^2=1-\frac94\sin^2\theta`, ref "6.91". "Fore–aft
    symmetric: no drag again (3-D d'Alembert, Exercise 6.39)." + `nb.code`: `print(ch06.cylinder_vs_sphere(2.0))` →
    *expect:* `{'cyl_pert': 0.25, 'sph_pert': 0.125, 'cyl_cp_min': -3.0, 'sph_cp_min': -1.25, 'cyl_speed_max': 2.0,
    'sph_speed_max': 1.5}`.
25. `nb.note` — **N85 [B]** "**Without coordinates.** With x = r e_r, cos θ = e_z·e_r, U = U e_z and d = −d e_z:" equation
    `\phi=Ur\cos\theta+\frac{d}{4\pi r^2}\cos\theta=\mathbf U\cdot\mathbf x-\frac{\mathbf d\cdot\mathbf x}{4\pi\lvert\mathbf x\rvert^3}=\Big(\mathbf U-\frac{\mathbf d}{4\pi\lvert\mathbf x\rvert^3}\Big)\cdot\mathbf x`,
    ref "6.92". "Valid for a stream in any direction — the starting point of C15." + `nb.code`: `x = np.array([0.3,
    -1.2, 0.9])` · `print(pf.sphere_potential_vector(x, [0, 0, 1.0], [0, 0, -2*np.pi]), sph.phi(np.hypot(0.3, -1.2),
    0.9))` → *expect:* equal (U = 1, a = 1, d = −2π e_z).
26. `nb.worked_example("a football and a lamp-post, both a = 0.11 m, in wind U = 10 m/s", "1. Fastest surface speed:
    lamp-post 2U = 20 m/s; football 1.5U = 15 m/s. 2. Lowest pressure, air ½ρU² = 60 Pa: lamp-post C_p = −3 → −180 Pa;
    football −1.25 → −75 Pa. 3. At twice the radius (r = 0.22 m) the stream is disturbed by (a/r)² = 25 % past the
    lamp-post but only (a/r)³ = 12.5 % past the ball. 4. Both have zero drag in ideal flow.")`
27. `nb.code` — *code:* `sph = pf.sphere(10.0, 0.11)` · `th = np.radians([0, 45, 90, 135, 180])` · `print(
    sph.velocity_spherical(0.11, th))` · `print(ch06.sphere_surface_cp(th))` · `print(np.max(np.abs(sph.psi(0.11*np.sin(th),
    0.11*np.cos(th)))))`. *expect:* u_r ≈ 0, u_θ = −15 sin θ → (0, −10.61, −15, −10.61, 0) · C_p (1, −0.125, −1.25, −0.125,
    1) · ≤ 1e-14 (the sphere is the surface ψ = 0). *explain:* 1. the sphere as an axisymmetric flow; 2. the surface
    velocity: tangential, 1.5U at the equator; 3. (6.91); 4. ψ = 0 on the surface.
28. `nb.check_agree` — **from scratch (curation §7):** the sphere's ψ assembled by hand from (6.86) and (6.88):
    `U, a = 10.0, 0.11`; `d = 2*np.pi*a**3*U` (D25 step 7); `r, t = 0.3, 0.7`; `psi_hand = 0.5*U*r**2*np.sin(t)**2 -
    d/(4*np.pi*r)*np.sin(t)**2`; `assert np.isclose(psi_hand, sph.psi(r*np.sin(t), r*np.cos(t)))`. Markdown: "Stream plus
    doublet, added by hand, is the library's sphere."
29. `nb.figure` — two panels (9 × 3.4 in): left, meridian-plane streamlines of the sphere (ψ contours at equal Δψ,
    mirrored below the axis), the sphere black; right, surface C_p against θ for the cylinder (1 − 4 sin²θ, dashed) and
    the sphere (1 − (9/4) sin²θ, solid), and an inset of the velocity perturbation along the side (θ = 90°) against r/a
    on log–log axes (slopes −2 and −3). Title "In 3-D the flow escapes sideways". *see:* "streamlines closer together
    only near the equator; a shallower C_p dip for the sphere; two straight lines in the inset." *read:* "the sphere's
    disturbance is weaker and dies faster: 1.5U vs 2U, (a/r)³ vs (a/r)²." *change:* "…the equal-Δψ streamlines of a
    Stokes ψ are not equally spaced in flux per unit area: the flux between them is 2πΔψ spread over a ring of
    circumference 2πR — they crowd near the axis; read speed from spacing × R, not spacing alone."
30. `nb.plotly` — `slider_figure` over r/a from 1 to 5 (17 steps): traces "cylinder u_θ/U at θ = 90°" and "sphere"
    against θ (the angular profile at radius r), with the far-stream reference 1; title prints the two perturbations.
    Title "3-D relief: the sphere's disturbance dies as (a/r)³".
31. `nb.md` — **What would change if…** "…the stream were not along the axis but across it? The flow would no longer be
    axisymmetric about the stream — but the sphere is round, so (6.92) handles any direction of U at once. For a
    non-spherical body we would need the full 3-D problem (Ch. 10 methods, or the singularity distributions of C14)."

**C14 — Bodies of revolution from axial singularities: the airship and the axial singularity method**
32. `nb.core("C14", "Given a body, find its sources: the axial singularity method, $\\psi_m=-\\sum_n\\frac{k_n}{4\\pi}(r^m_{n-1}-
    r^m_n)+\\tfrac12UR_m^2=0$", question="An airship hull is drawn on paper. Which sources and sinks, hidden on its axis,
    produce exactly that shape?")`
33. `nb.md` — **The problem in plain words:** "So far we picked sources and found the body. A designer works the other
    way: the shape of a hull or a fuselage is given, and the flow round it is wanted. Hide unknown sources and sinks along
    the axis, demand that the body be a stream surface at as many points as there are unknowns, and solve — a linear
    system. The same 'unknown strengths + one condition per point' idea runs every panel code in aerodynamics."
34. `nb.md` — **The idea** (ASCII):
    ```
    point source Q at the nose + line sink k over length a + stream U   ─►  closed airship when Q = ak (6.95)
    inverse: N axial segments of unknown k_n ──ψ = 0 at N body points──► N × N linear system ──► k_n
    closure check: Σ k_n Δξ ≈ 0 (sources and sinks cancel, or the body does not close)
    ```
35. `nb.note` — **N86 [C]** "**Closed bodies need zero net source:** everything the sources emit must be swallowed
    again, or fluid would stream out of the body forever (the half-body of C05 is open for that reason). A sink spread
    over a length gives a gently tapering tail."
36. `nb.md` — 🔁 reminders: substitution in an integral (P106) + gloss "**z − ξ = R cot α** with z, R fixed: dξ = R dα/sin²α
    because d(cot α)/dα = −1/sin²α, and the limits turn around"; point-source ψ (6.87) (N81); `brentq` (P108); `quad`
    (P87).
37. `nb.derivation("D26", …)` — Part F D26 (10 steps), ref "6.95".
38. `nb.note` — **N87 [B], N88 [B], N89 [B]** "**The line sink and the airship.** An element kdξ of a line sink of
    density k [m²/s] on the axis contributes dψ = (k dξ/4π) cos α (the opposite sign of a source (6.87)), so" equation
    `\psi_{\text{sink}}=\frac{k}{4\pi}\int_0^a\cos\alpha\,d\xi`, ref "6.93", "which the substitution turns into" equation
    `\psi_{\text{sink}}=\frac{kR}{4\pi}\Big[\frac1{\sin\theta}-\frac1{\sin\alpha_1}\Big]=\frac{k}{4\pi}(r-r_1)`, ref "6.94"; "with
    a point source Q at O, k = Q/a and the stream," equation `\psi=-\frac{Q}{4\pi}\cos\theta+\frac{Q}{4\pi a}(r-r_1)+\tfrac12Ur^2\sin^2\theta`,
    ref "6.95". "Q = 1 m³/s, a = 1 m → k = 1 m²/s."
39. `nb.worked_example("a small airship", "Q = 1 m³/s, a = 1 m, U = 1 m/s. 1. Closure: k = Q/a = 1 m²/s, so the sink
    swallows 1 × 1 = 1 m³/s = Q ✓. 2. Nose: on the axis ahead of O the source's push Q/4πz² nearly balances U; with the
    line sink's pull the stagnation point is at z = −0.252 m (a pure source would give −√(Q/4πU) = −0.282 m). 3. Tail:
    just behind the sink's end, z = 1.070 m. 4. Length 1.322 m — `brentq` on u_z(z) = 0 along the axis.")`
40. `nb.code` — *code:* `R = np.array([0.1, 0.3]); Z = np.array([0.4, -0.5])` · `print(ch06.line_sink_stream_function(R, Z,
    1.0, 1.0), ch06.line_sink_stream_function(R, Z, 1.0, 1.0, method="quad"))` · `A = ch06.airship(1.0, 1.0, 1.0)` ·
    `print(A["stagnation_front"], A["stagnation_rear"], A["length"], A["R_max"], A["closure"])`. *expect:* closed form =
    quadrature to 1e-10 · −0.2521, 1.0696, 1.3217, the maximum radius (printed), closure 0. *explain:* 1. (6.94) against
    a numerical integral of (6.93); 2. the airship's nose, tail and thickness from its axis stagnation points and the
    ψ = 0 surface.
41. `nb.primer("collocation and the condition number", "To pin down N unknown strengths, demand the condition (here ψ = 0)
    at N chosen points: N equations, N unknowns — **collocation**. How trustworthy the answer is depends on the matrix:
    `np.linalg.cond(A)` is the factor by which relative errors in the data can grow in the solution. Around 10³ is
    harmless; above 10¹⁰ the digits are noise.", code="for N in (4, 8, 16):\n    x = np.linspace(0, 1, N)
    # N collocation points\n    V = np.vander(x, N)                    # a notoriously ill-conditioned matrix\n    print(N,
    f'{np.linalg.cond(V):.1e}')          # grows explosively with N")` — *expect:* 9.9e+01, 2.7e+05, 3.1e+12.
42. `nb.derivation("D27", …)` — Part F D27 (6 steps), no number (ref "" — the Fig. 6.29 system).
43. `nb.code` — *code:* `zb, Rb = ch06.axisym_body_target("rankine_oval", N=40)` · `sol = pn.axial_singularity_solve(zb, Rb,
    U=1.0, N=40)` · `print(sol["cond"], sol["net_strength"], np.max(np.abs(sol["residual"])))` · `for N in (10, 20, 40, 80):
    s = ch06.axial_state("rankine_oval", N=N); print(N, s["body_error"], s["net_strength"], s["cond"])` · `print(
    ch06.axial_state("sphere", N=40))`. *expect:* residual ≲ 1e-10 at the collocation points; body error falling with N;
    net strength → 0; cond growing; the sphere: small residual but a large condition number (fenced: "a sphere is a point
    doublet — a line of sources can only approximate it"). *explain:* 1. a target body; 2. the collocation system and its
    solve; 3. convergence with N; 4. the blunt body's warning sign.
44. `nb.check_agree` — **from scratch (curation §7):** the influence matrix by a double loop: `N = 20`; `zb, Rb =
    ch06.axisym_body_target("rankine_oval", N=N)`; `xi = np.linspace(zb.min(), zb.max(), N + 1)` (segment ends — the
    builder uses `sol["z_nodes"]` so both use the same nodes); `A = np.zeros((N, N))`; `for m in range(N): for n in range(N):
    r0 = np.hypot(zb[m] - xi[n], Rb[m]); r1 = np.hypot(zb[m] - xi[n+1], Rb[m]); A[m, n] = (r0 - r1)/(4*np.pi)` (a unit
    source segment's contribution with the minus sign moved to the right side); `k = np.linalg.solve(A, 0.5*1.0*Rb**2)`;
    `assert np.allclose(k, pn.axial_singularity_solve(zb, Rb, U=1.0, z_nodes=xi)["k"])`. Markdown: "Our own N × N matrix
    gives the library's strengths."
45. `nb.note` — **N90 [B]** "**The same idea on the body surface (our extension, not in the book):** cover a 2-D body with
    N straight panels carrying unknown uniform source strengths λ_j, demand u·n = 0 at each panel's midpoint, and solve
    (Hess & Smith 1967). On a circle the midpoint C_p converges to 1 − 4 sin²θ (6.35) as N grows, with error ∝ 1/N²,
    and Σλ_j s_j = 0 (closed body)." + `nb.code`: `for N in (8, 16, 32, 64): print(N, ch06.panel_cp_error(N))` +
    `print(observed_order(1/np.array([8, 16, 32, 64]), [ch06.panel_cp_error(N) for N in (8, 16, 32, 64)]))` → *expect:*
    errors falling ×4, order 2.0 ± 0.25.
46. `nb.figure` — three panels (11 × 3.2 in): left, the airship (6.95) in the body frame (streamlines, the ψ = 0 body)
    and, mirrored below the axis, the same body in the fluid frame (the stream subtracted: instantaneous streamlines
    looping from nose to tail); centre, the Rankine-oval target (dashed) and the computed ψ = 0 surface (solid) with k_n
    bars along the axis (teal sources, rose sinks); right, body error and condition number against N (log–log). Title
    "Hidden sources draw the body". *see:* "a streamlined hull; bars positive at the front, negative at the back; error
    falling while the condition number climbs." *read:* "sources push the flow out at the nose, sinks let it close at the
    tail; the bars sum to zero; more segments fit the shape better but the system gets touchier." *change:* "…a blunter
    target (the sphere): the error plateaus and the condition number explodes — its exact source is a point doublet,
    which a smooth line cannot copy."
47. `nb.plotly` — `slider_figure` over N from 4 to 64 (16 steps): the target, the computed body and the k_n bars (as a
    step line), the title printing the body error and cond. Title "The inverse method converging".
48. `nb.explainer("axial_singularity_bodies", heading="Given a shape, which sources draw it?", why="Convergence is a
    sequence, not a picture: raise N and watch the computed ψ = 0 surface snap onto the target while the k_n bars settle
    and the condition number climbs; switch to panels to see the same idea on a 2-D body's surface.", tries=["Preset
    'airship (6.95)': the bars show a point source and a flat line sink.", "Target 'Rankine oval', N from 4 to 64: the
    error falls, cond rises.", "Target 'sphere': the ⚠️ status shows why a blunt body is hard.", "Mode 'panels': C_p on a
    circle converges to 1 − 4 sin²θ."])`
49. `nb.md` — **What would change if…** "…we used vortices on the axis or on the surface instead of sources? Sources give
    thickness but no lift; vortex panels add circulation — the lifting version is how Ch. 14's panel codes compute an
    airfoil's lift (with the Kutta condition fixing the total Γ)."

### A.9 §6.9 Three-Dimensional Potential Flow and Apparent Mass — R31 R32, C15 (+N91–N105 · D28 D29 D30 D31 · E9)
1. `nb.section("6.9", "Three-Dimensional Potential Flow and Apparent Mass", intro="**What is this section about?** Steady
   ideal flow gives no drag in 3-D either. But when a body *accelerates*, the fluid round it must be accelerated too,
   and the pressure field that does this pushes back on the body. For a sphere the push is exactly that of an extra
   half of the displaced fluid's mass: the added mass.")`
2. `nb.recap("R31", "The 3-D potential", "u ≡ ∇φ with all three components, so the z-component is w ≡ ∂φ/∂z *(Eq. 3.17)*.
   ⚠️ In this section w is a velocity component, not the complex potential of §6.4.", where="Ch. 3 §3.4")`
3. `nb.recap("R32", "Unsteady Bernoulli between the surface and far away", "For unsteady irrotational flow of constant
   density with the fluid at rest and at pressure p∞ far away: $\\Big[\\frac{\\partial\\phi}{\\partial t}+\\frac12\\lvert\\nabla\\phi\\rvert^2+
   \\frac p\\rho\\Big]_{\\text{sphere's surface}}=\\frac{p_\\infty}{\\rho}$ *(Eq. 6.99, from 4.75)*. ⚠️ ∂φ/∂t is taken at a *fixed* point in
   space.", where="Ch. 4 §4.9")` + `nb.code`: `from fluidpy.core import bernoulli as be` · `print(
   be.unsteady_bernoulli_pressure(dphi_dt=-0.05, speed=1.0, rho=1000.0, g=0.0, C=0.0))` → *expect:* −450.0 Pa
   (p/ρ = −∂φ/∂t − ½|u|² = 0.05 − 0.5).

**C15 — The accelerating sphere: surface pressure and added mass (6.109)**
4. `nb.core("C15", "Newton's law for a submerged sphere, $\\mathbf F_E=\\big(m+\\frac{2\\pi}{3}\\rho a^3\\big)\\frac{d\\mathbf u_s}{dt}$
   (6.109)", question="An ideal fluid exerts no drag on a steadily moving body. So why is it hard to shake a ball under
   water — and why does a released bubble not shoot up infinitely fast?")`
5. `nb.md` — **The problem in plain words:** "Push a beach ball under water and let go: it accelerates up, but far less
   violently than its tiny mass would suggest. A bubble has almost no mass, yet it does not leave with infinite
   acceleration. A submarine, a fish, a kite and a ship in waves all feel the same effect: to speed up, a body must also
   speed up the fluid it pushes aside."
6. `nb.md` — **The idea** (ASCII):
   ```
   pressure on a moving sphere = steady part (∝ speed², symmetric front/back)  → no net force (d'Alembert)
                               + acceleration part (∝ du_s/dt, + in front, − behind) → force −M du_s/dt
   M = ½ × (mass of the displaced fluid) = (2π/3)ρa³       →   F_E = (m + M) du_s/dt
   ```
7. `nb.note` — **N91 [C]** "**3-D d'Alembert and beyond.** A closed body moving steadily through ideal fluid feels no
   drag in 3-D too (Exercise 6.39); forces appear when the motion is unsteady or vorticity is present (Ch. 14). The
   inertia of the fluid is what vehicles, fish and bubbles feel: apparent or added mass."
8. `nb.note` — **N92 [B]** "**The set-up** (Fig. 6.30 idea):" ASCII sketch of a sphere of radius a at x_s(t) with velocity
   $\mathbf u_s=d\mathbf x_s/dt$ and a known acceleration; the fluid is at rest far away at p∞; an external force F_E acts on
   the sphere; we want the fluid's force F_s. "Think of a submarine or a fish manoeuvring."
9. `nb.md` — 🔁 reminders: Galilean frames (Ch. 3 §3.3), the coordinate-free sphere potential (6.92) (N85), chain rule
   along a path (P91), gradient of 1/distance (P140) + glosses "**a function of x − x_s:** ∂/∂x_s = −∇" and "**∇ of
   u·ξ/\|ξ\|³:** product rule with ∇(u·ξ) = u and ∇\|ξ\|⁻³ = −3ξ/\|ξ\|⁵".
10. `nb.derivation("D28", …)` — Part F D28 (5 steps), ref "6.97".
11. `nb.note` — **N93 [B], N94 [B]** "With x → x − x_s and no stream:" equation
    `\phi=-\frac{1}{4\pi\lvert\mathbf x-\mathbf x_s(t)\rvert^3}\mathbf d\cdot(\mathbf x-\mathbf x_s(t))`, ref "6.96", "and the dipole
    following the motion, d(t) = 2πa³u_s:" equation `\phi(\mathbf x,\mathbf x_s,\mathbf u_s)=-\frac{a^3}{2\lvert\mathbf x-\mathbf x_s\rvert^3}\mathbf u_s\cdot(\mathbf x-\mathbf x_s)=-\frac{a^3}{2\lvert\boldsymbol\xi\rvert^3}\mathbf u_s\cdot\boldsymbol\xi`,
    ref "6.97". + `nb.code`: `xs, us, a = np.zeros(3), np.array([0.0, 0.0, 2.0]), 0.1` · `e = np.array([[0.6, 0.0, 0.8],
    [0.0, 1.0, 0.0]]).T` (two surface directions) · `u_a = pf.moving_sphere_velocity(xs[:, None] + a*e, xs, us, a)` ·
    `print(np.sum(u_a*e, axis=0), us @ e)` → *expect:* equal pairs (1.6, 0.0): "the fluid's normal velocity equals the
    sphere's, (6.16) with U_s = u_s".
12. `nb.note` — **N95 [B]** "**The force** is the pressure integrated over the surface (the 3-D form of (6.55); a uniform p∞
    gives nothing on a closed surface):" equation `\mathbf F_s=-\int_{\text{sphere's surface}}(p-p_\infty)\mathbf n\,dA`, ref
    "6.98".
13. `nb.derivation("D29", …, check_src=…)` — Part F D29 (14 steps, ★★★), ref "6.105".
14. `nb.note` — **N96 [B] – N101 [B]** (one note; the chain of D29 with each equation and number): `\frac{p_a-p_\infty}{\rho}=-\Big(\frac{\partial\phi}{\partial t}\Big)_a-\frac12\lvert\nabla\phi\rvert_a^2`
    (6.100) · `\frac{\partial\phi}{\partial t}=-\mathbf u\cdot\mathbf u_s-\frac{a^3}{2\lvert\boldsymbol\xi\rvert^3}\boldsymbol\xi\cdot\frac{d\mathbf u_s}{dt}`
    (6.101) · `\Big(\frac{\partial\phi}{\partial t}\Big)_a=-\mathbf u_a\cdot\mathbf u_s-\frac a2\mathbf e_\xi\cdot\frac{d\mathbf u_s}{dt}`
    (6.102) · `\nabla\phi=-\frac{a^3}{2}\Big[-\frac{3(\mathbf x-\mathbf x_s)}{\lvert\mathbf x-\mathbf x_s\rvert^5}\mathbf u_s\cdot(\mathbf x-\mathbf x_s)+\frac{\mathbf u_s}{\lvert\mathbf x-\mathbf x_s\rvert^3}\Big]`
    (6.103) · `\mathbf u_a=\tfrac32(\mathbf u_s\cdot\mathbf e_\xi)\mathbf e_\xi-\tfrac12\mathbf u_s` (6.104) ·
    `\frac{p_a-p_\infty}{\rho}=\frac12\lvert\mathbf u_s\rvert^2\Big(\frac94\frac{(\mathbf u_s\cdot\mathbf e_\xi)^2}{\lvert\mathbf u_s\rvert^2}-\frac54\Big)+\frac a2\mathbf e_\xi\cdot\frac{d\mathbf u_s}{dt}`
    (6.105). "⚠️ The book's middle expression of (6.104) prints −u_s/a³ inside the bracket; +u_s/a³ is what (6.103) gives
    and what the final form needs (sympy in D29's check: the printed bracket would give +½u_s). Check with your hands:
    at the sphere's side (e ⟂ u_s) the fluid moves at −½u_s — backwards, filling in." + `nb.code`: `print(
    ch06.moving_sphere_dphidt_sym())` → *expect:* the dict with `correct_bracket` equal to (3/2)(u·e)e − ½u and
    `printed_bracket` equal to (3/2)(u·e)e + ½u.
15. `nb.note` — **N102 [B]** "**Steady motion:**" equation `\Big(\frac{p_a-p_\infty}{\frac12\rho\lvert\mathbf u_s\rvert^2}\Big)_{\text{steady}}=\frac94\cos^2\theta_s-\frac54=1-\frac94\sin^2\theta_s`,
    ref "6.106" "— (6.91) again: a sphere moving steadily through still fluid feels the same pressures as a still sphere in
    a stream (Galilean invariance), and no drag. ⚠️ θ_s is measured from the sphere's velocity, i.e. from its front; in
    (6.91) θ runs from the downstream side — sin² is the same either way."
16. `nb.primer("surface integrals on a sphere", "On a sphere of radius a use the polar angle θ (from the chosen axis) and
    the azimuth φ: the area element is $dA=a^2\\sin\\theta\\,d\\theta\\,d\\varphi$ and the outward normal is $\\mathbf e_r=(\\sin\\theta
    \\cos\\varphi,\\ \\sin\\theta\\sin\\varphi,\\ \\cos\\theta)$. Integrate φ first: anything with cos φ or sin φ vanishes. A handy result:
    $\\int_0^\\pi\\cos^2\\theta\\sin\\theta\\,d\\theta=\\big[-\\tfrac13\\cos^3\\theta\\big]_0^\\pi=\\tfrac23$.", code="from scipy.integrate import
    dblquad\narea = dblquad(lambda ph, th: np.sin(th), 0, np.pi, 0, 2*np.pi)[0]    # ∫∫ sin θ dφ dθ on the unit
    sphere\nI = dblquad(lambda ph, th: np.cos(th)**2*np.sin(th), 0, np.pi, 0, 2*np.pi)[0]\nprint(area, 4*np.pi, I,
    4*np.pi/3)                   # 12.566 = 4π, 4.189 = 2π × 2/3")` — *expect:* `12.566370614359172 12.566370614359172
    4.1887902047863905 4.1887902047863905`.
17. `nb.derivation("D30", …)` — Part F D30 (9 steps), ref "6.108".
18. `nb.note` — **N103 [B], N104 [B]** "With the acceleration along e_z," equation
    `\mathbf F_s=-\rho\frac a2\Big\lvert\frac{d\mathbf u_s}{dt}\Big\rvert\int_{\theta=0}^{\pi}\int_{\varphi=0}^{2\pi}\cos\theta(\mathbf e_x\sin\theta\cos\varphi+\mathbf e_y\sin\theta\sin\varphi+\mathbf e_z\cos\theta)a^2\sin\theta\,d\varphi\,d\theta`,
    ref "6.107", "and after the φ-integration" equation `\mathbf F_s=-\frac23\pi\rho a^3\frac{d\mathbf u_s}{dt}=-M\frac{d\mathbf u_s}{dt},\qquad M=\frac{2\pi a^3\rho}{3}`,
    ref "6.108". "⚠️ The book's (6.108) still shows dφ in the integrand after the φ-integration has produced its 2π. M is
    **half the mass of the displaced fluid**: a = 0.1 m in water, M = 2.094 kg vs 4.189 kg displaced. In 2-D a circular
    cylinder's added mass is ρπa² per unit depth — the *whole* displaced mass (Exercise 6.45; `ch06.cylinder_added_mass(0.1,
    1000.0)` = 31.42 kg/m). In general M is a tensor, $(F_s)_i=-M_{ij}\,d(u_s)_j/dt$."
19. `nb.note` — **N105 [C]** "**Why the fluid pushes back:** in front of an accelerating sphere the fluid must be pushed
    aside faster and faster, behind it the hole must be filled faster and faster; the pressure that does this is high in
    front and low behind."
20. `nb.primer("Green's first identity", "For any smooth φ, $\\nabla\\cdot(\\phi\\nabla\\phi)=\\lvert\\nabla\\phi\\rvert^2+\\phi\\nabla^2\\phi$ (the
    product rule for a divergence, P113). If ∇²φ = 0 the last term drops, so the kinetic-energy density ½ρ\\|∇φ\\|² is a
    divergence — and Gauss turns its volume integral into a surface integral. For the fluid outside a body the surface
    is the body (with the normal pointing **into** the body, out of the fluid) plus a far sphere.", code="x, y, z =
    sp.symbols('x y z')\nphi = x*y - z**2 + x**2/2 + y**2/2         # a harmonic polynomial: ∇²φ = 1 + 1 − 2 = 0\ngrad =
    [sp.diff(phi, v) for v in (x, y, z)]\nlhs = sum(sp.diff(phi*g, v) for g, v in zip(grad, (x, y, z)))   # ∇·(φ∇φ)\nprint(
    sp.simplify(lhs - sum(g**2 for g in grad)))    # 0: it equals |∇φ|²")` — *expect:* `0`.
21. `nb.derivation("D31", …, check_src=…)` — Part F D31 (10 steps, ★★★), no book number (ref "6.108" — "Exercise 6.49's
    route, which the book leaves as an exercise").
22. `nb.md` — gloss: "**Kinetic energy and added mass.** A body of mass m moving at U has kinetic energy ½mU²; the fluid it
    sets moving carries ½MU² more. Pushing both up to speed costs ½(m + M)U² — the same M as the force route."
23. `nb.worked_example("a ball, a sphere of water and a bubble", "a = 0.1 m, water ρ = 1000 kg/m³. 1. Displaced mass ρ(4/3)πa³
    = 4.189 kg; added mass M = 2.094 kg. 2. A steel ball (7800 kg/m³): m = 32.67 kg, behaves as m + M = 34.77 kg; released
    from rest it accelerates at (m − ρV)g/(m + M) = (32.67 − 4.19) × 9.81/34.77 = 8.04 m/s² (8.55 m/s² without added mass).
    3. A bubble (m ≈ 0): buoyancy ρVg = 41.1 N accelerates only M: 41.1/2.094 = 19.6 m/s² = **2g**, not infinity.")`
24. `nb.code` — *code:* `a, rho = 0.1, 1000.0` · `print(pf.added_mass_sphere(a, rho), ch06.added_mass_by_energy(a, rho),
    ch06.added_mass_by_energy(a, rho, method="volume"))` · `us = np.array([0.0, 0.0, 1.0]); dus = np.array([0.0, 0.0,
    2.0])` · `F = pf.sphere_force_quadrature(lambda e: pf.moving_sphere_surface_pressure(e, us, dus, a, rho), a)` ·
    `print(F, -pf.added_mass_sphere(a, rho)*dus)` · `sp_ = pf.moving_sphere_surface_pressure(np.array([0, 0, 1.0]), us,
    dus, a, rho, split=True)` · `print(sp_)` · `print(ch06.sphere_motion(0.0, a, rho, g=9.81, mode="bubble",
    t_eval=np.array([0.0, 0.01]))["du_dt"][0]/9.81)`. *expect:* `2.0944 2.0944 2.0944` (pressure route, energy by the
    surface form, energy by the volume integral) · F = (0, 0, −4.1888) N twice · at the front point: steady
    ½ρ\|u\|²(9/4 − 5/4) = 500 Pa, acceleration ρ(a/2)·2 = 100 Pa, total 600 Pa · `2.0`. *explain:* 1. M three ways; 2. the
    force by quadrature over the sphere = −M du_s/dt; 3. the two pressure parts at the nose; 4. a bubble starts at 2g.
25. `nb.check_agree` — **from scratch (curation §7):** a Gauss–Legendre (cos θ) × trapezoid (φ) sum of −(p − p∞)n dA:
    `mu, w = np.polynomial.legendre.leggauss(24)` (nodes in μ = cos θ); `ph = np.linspace(0, 2*np.pi, 48,
    endpoint=False)`; `Mu, Ph = np.meshgrid(mu, ph, indexing="ij")`; `st = np.sqrt(1 - Mu**2)`; `E = np.array([st*np.cos(Ph),
    st*np.sin(Ph), Mu])` (unit normals); `p = pf.moving_sphere_surface_pressure(E.reshape(3, -1), us, dus, a, rho).reshape(
    Mu.shape)`; `Fz = -np.sum(p*E[2]*w[:, None])*a**2*(2*np.pi/48)` (dA = a² dμ dφ); `assert np.isclose(Fz,
    -pf.added_mass_sphere(a, rho)*2.0)`. Markdown: "Our own surface sum is −M du_s/dt: the steady part cancelled, the
    acceleration part did not."
26. `nb.figure` — two panels (9 × 3.4 in): left, the surface pressure p − p∞ against θ_s (angle from the direction of
    motion) for u_s = 1 m/s and du_s/dt = 2 m/s² (a = 0.1 m, water): steady part (teal, symmetric about 90°), acceleration
    part (orange, antisymmetric: + in front, − behind) and total (black); right, velocity against time for a bubble and a
    steel ball released from rest, with added mass (solid) and without (dashed, "vacuum-like"). Title "Speed pushes
    symmetrically; acceleration pushes back". *see:* "a symmetric teal curve, an antisymmetric orange one; solid curves
    rising more slowly than dashed ones." *read:* "only the orange part survives the force integral: F = −M du_s/dt; the
    lighter the body, the more its motion is ruled by M (a bubble's initial acceleration 2g instead of ∞)." *change:*
    "…the body a long cylinder moving sideways: M doubles relative to its displaced mass (the whole displaced mass per
    unit length) — a bubble column would start at g, not 2g."
27. `nb.animation` — **an oscillating sphere** (`player="video"`, 60 frames, FAST 30): x_s = A sin ωt (A = 0.05 m, ω = 2
    rad/s), the sphere section with surface pressure colour (steady teal + acceleration orange shading), an arrow for
    F_s, and a panel with F_s(t) and −M du_s/dt (dashed) on top of each other; the steady part flickers twice per period,
    the force follows the acceleration. + see/read/change markdown: "the force is in phase with the acceleration, never
    with the speed: no energy is lost on average (no drag), but the inertia is larger."
28. `nb.plotly` — `slider_figure` over the density ratio ρ_body/ρ from 0 (bubble) to 8 (steel) (17 steps): velocity after
    release (u(t) from `sphere_motion`) with and without added mass, title printing the initial acceleration in units of
    g. Title "When added mass matters: light bodies".
29. `nb.live` — `live(lambda m, a, F_E: …, m=(0.0, 40.0, 0.5), a=(0.02, 0.3, 0.01), F_E=(0.0, 50.0, 1.0))` that plots
    `ch06.sphere_motion(m, a, rho=1000.0, F_E=lambda t: F_E, t_eval=np.linspace(0, 2, 101))["u"]` with and without added
    mass and prints (m + M). (Paired with the plotly figure above — the page shows that one.)
30. `nb.explainer("added_mass_sphere", heading="Why does an ideal fluid resist acceleration?", why="Set speed and
    acceleration independently, even at an angle, and see the two pressure parts on the sphere with their separate
    forces — one always zero, the other always −M du_s/dt; the ball and bubble modes integrate (6.109) in time.",
    tries=["Preset 'steady motion': the teal pressure is strong, the force bar zero.", "Preset 'starting from rest':
    only the orange part, and the force equals −M du_s/dt.", "Turn the acceleration 90° away from the velocity: the
    force follows the acceleration, not the motion.", "Preset 'air bubble': the transport shows it starting at 2g."])`
31. `nb.md` — **What would change if…** "…the body were not a sphere? M becomes a tensor: a flat plate accelerated
    broadside drags along far more fluid than edgewise — why a kite or a sail reacts differently to gusts from different
    directions, and why fish shape their bodies as they do (Ch. 16)."

### A.10 §6.10 Concluding Remarks — N106, S01–S07, summary
1. `nb.section("6.10", "Concluding Remarks", intro="**What is this section about?** Where ideal-flow theory came from,
   where it fails, and where it lives on.")`
2. `nb.note` — **N106 [C]** "**Two and a half centuries.** Euler, Bernoulli, d'Alembert, Lagrange, Stokes, Helmholtz,
   Kirchhoff and Kelvin built potential-flow theory; the same mathematics describes heat conduction, elasticity and
   electrostatics. Its zero-drag prediction contradicted every observation until Prandtl's boundary layer (Ch. 9)
   explained where viscosity hides. It remains the working tool for pressure on streamlined bodies (Ch. 9) and lift
   (Ch. 14, with conformal maps and the Kutta condition)." (Curation parent C15: this closes the force story d'Alembert →
   lift → added mass.)
3. `nb.pointer("Exercises 6.1–6.5, 6.8 and 6.33 check the §6.2 statements (delta-source integrals, orthogonality of ∇φ and
   ∇ψ, Bernoulli from (6.1)); their ideas are D03 and D04 above (S01).")`
4. `nb.pointer("Exercises 6.6–6.7, 6.9, 6.14–6.15, 6.17 and 6.22: harmonic polynomials (generated by
   `ch06.harmonic_polynomials`, N20), the doublet's stream function (stated in N25 and N47), sketched flows and a
   vortex-pair doublet — left to the reader (S02).")`
5. `nb.pointer("Exercises 6.10–6.13, 6.18–6.20 and 6.23–6.24: forces on held singularities, the half-body's width by a
   mass balance (D07 step 9) and its zero net force (N29), the Rankine oval (E1's preset) — the rest left to the reader
   (S03).")`
6. `nb.pointer("Exercises 6.21, 6.25–6.27 and 6.29–6.32: Blasius for the cylinder (E5's cylinder preset) and a
   control-volume proof of Kutta–Zhukhovsky (`ch06.cv_force_on_body`, R18); corners, images and vortex pairs left to
   the reader (S04).")`
7. `nb.pointer("Exercises 6.16 and 6.41 are kitchen experiments (a paper cylinder vs an airfoil in a draught; a vacuum
   nozzle over sugar grains — a sink and its image). Try them (S05).")`
8. `nb.pointer("Exercises 6.34–6.40 and 6.42–6.43: the 3-D source and doublet (D25), 3-D d'Alembert (N91), the airship's
   length (D26 and the worked example of C14) — the rest left to the reader (S06).")`
9. `nb.pointer("Exercises 6.44–6.49: cavity collapse, the cylinder's added mass (N104), a bubble's 2g start (E9's bubble
   preset), oscillation frequencies, and the kinetic-energy route to added mass (D31) (S07).")`
10. `nb.summary(clicked=[…15…], feeds_forward=[…], left_out=[…])` — **clicked** (one sentence per CORE): C01 ω = 0 makes
    the viscous force vanish, so outside boundary layers and wakes the flow obeys ∇·u = 0 and ρDu/Dt = −∇p (6.1) · C02
    ω_z = −∇²ψ (6.4): irrotational means Laplace for ψ and φ, with vortices and sources as delta sources · C03 Laplace is
    linear, so flows add, and any streamline can be a wall (6.16) · C04 a source and a sink merged at fixed 2mε give the
    doublet φ = |d| cos θ/2πr (6.29) · C05 stream + source = half-body: stagnation at m/2πU, width m/2U far downstream
    (6.31) · C06 the cylinder's C_p = 1 − 4 sin²θ (6.35) is symmetric: zero drag, d'Alembert's paradox · C07 circulation
    moves the stagnation points (sin θ = −Γ/4πaU (6.38)) and gives L = ρUΓ (6.40) · C08 images turn walls into symmetry
    lines; a passing eddy is felt as suction then over-pressure · C09 w = φ + iψ (6.42) is analytic; dw/dz = u − iv (6.45)
    · C10 Blasius (6.60) + residues: D = 0 and L = ρUΓ for any body (6.62) · C11 analytic maps keep angles (6.64);
    z = ζ + b²/ζ (6.65) carries the circle onto an ellipse — take the outside root · C12 on a grid Laplace is "average of
    your neighbours" (6.72); Gauss–Seidel relaxes to it; stop on the residual · C13 the Stokes ψ obeys (6.77), not
    Laplace; the sphere has max speed 1.5U and C_p = 1 − (9/4) sin²θ (6.91) · C14 unknown axial sources + ψ = 0 at N
    points = one linear solve for a given body · C15 acceleration pressure gives F_s = −M du_s/dt with M = 2πρa³/3
    (6.108). **feeds forward:** velocity potentials of water waves and images at the bottom (Ch. 7); ideal outer flow
    U_e(x) that drives boundary layers, wedge flows (6.46) and separation (Ch. 9); elliptic solvers and panel methods
    (Ch. 10); ψ–vorticity inversion, deformation fields, coastal images and bounded-basin Poisson problems (Ch. 13); the
    Kutta condition, Zhukhovsky airfoils and induced drag (Ch. 14); added mass of bubbles, fish and structures (Ch. 16).
    **left out:** exercise details (S01–S07); the measured cylinder pressure curve (qualitative band only); the moment
    form of Blasius's theorem; non-circular added-mass tensors.

### A.11 Placement table (ID → notebook section → host block → call)
Every curation ID appears once. `note` = `nb.note`; `md` = a markdown cell inside the block; `recap` = `nb.recap` before
the block's `nb.core`; `deriv` = `nb.derivation`; `pointer` = `nb.pointer`.
- §6.1: R01 recap (before C01) · C01 core · N02 note (item 8), N01 note (12), N03 note (14), N04 note (15) · D01 deriv.
- §6.2: R02, R03, R04, R06, R07, R08, R09 recap (before C02) · C02 core · N05, N06, N07, N13, N08, N14, N10, N11, N12, N09,
  N18, N19, N15 note · D02, D03, D04 deriv · R05 recap (before C03) · C03 core · N16, N17 note · D05 deriv.
- §6.3: R10, R11 recap (before C04) · C04 core · N20, N21, N22, N23, N24, N25, N26 note · D06 deriv · C05 core · N27, N28,
  N29 note · D07, D08 deriv · **E1** explainer (C05) · R12, R13, R14 recap (before C06) · C06 core · N30, N31 note · D09
  deriv · C07 core · N32, N33, N34, N35, N36, N38, N37 note · D10, D11 deriv · **E2** (C07) · R15, R16, R17 recap (before
  C08) · C08 core · N39, N40 note · D12, D13 deriv · **E3** (C08).
- §6.4: C09 core · N41, N42, N43, N44, N45, N46, N47, N48, N49, N50, N51, N52 note · D14, D15 deriv · **E4** (C09).
- §6.5: R18 recap (before C10) · C10 core · N53, N54, N55, N56, N57, N58, N59, N60, N63, N61, N62 note · D16, D17, D18
  deriv · **E5** (C10).
- §6.6: C11 core · N64, N65, N66, N67, N68, N69, N70 note · D19, D20, D21 deriv · **E6** (C11).
- §6.7: R19, R20, R21 recap (before C12) · C12 core · N71, N72, N73, N74 note · D22, D23 deriv · **E7** (C12) · pointer to
  N90.
- §6.8: R22, R23, R24, R25, R26, R27, R28, R29, R30 recap (before C13) · C13 core · N75, N76, N77, N78, N79, N80, N81, N82,
  N83, N84, N85 note · D24, D25 deriv · C14 core · N86, N87, N88, N89, N90 note · D26, D27 deriv · **E8** (C14).
- §6.9: R31, R32 recap (before C15) · C15 core · N91, N92, N93, N94, N95, N96, N97, N98, N99, N100, N101, N102, N103, N104,
  N105 note · D28, D29, D30, D31 deriv · **E9** (C15).
- §6.10: N106 note · S01, S02, S03, S04, S05, S06, S07 pointer · summary. B1 `flow_net_sources_vortices` not embedded
  (backup).
Counts: 15 `nb.core`, 32 `nb.recap`, 106 NOTE ids placed (as 88 `nb.note` calls — some notes carry several ids), 7
`nb.pointer` for SKIP + 1 for N90's section, 31 `nb.derivation`, 16 `nb.primer`, 9 `nb.explainer`, 4 animations, 9
plotly figures, 1 live cell, 27 static figures.

## Part B — explainer storyboards

Common to all ten: created with `tools/new_viz.py`; `<meta name="viz:chapter" content="ch06">`; tabs Walkthrough ·
Explore · Explain · Derivation · Equations · Code · Check; `autoplay: false` and `play: false` on every walkthrough step
that quotes numbers; every displayed number is computed by a JS function that mirrors a `ch06` callable and is proved by
`selftest()` parity rows (`py:` expressions use only `ch06.…`, `np.pi`, numbers, strings, lists and dicts — no builtins,
no lambdas; results indexed down to one number; element lists are dict specs, Part C convention 2). Explain is
"Explanation & interpretation" in numbered sections built with `Viz.work.step / line / box / table / hint / interpret`,
modelled on `forced_damped_vibrations.html` and `amplitude_phase_second_order_II_3.html` (0 what the views show and what
each colour means · 1…n every displayed quantity from the controls, "formula = substituted = result — why", results
boxed · a section per view hidden on phones or a hint · the values at the current time (live) · **Reading the current
setting**). Derivation steps are copied from Part F (same `did` titles, same step count; phones shorten *why* to its
first sentence; plain-text *why* and *watch* never contain raw TeX). Every tour, Explain, Derivation, notes, status,
equation and quiz text that names a book equation **writes it out** next to its number (convention 10; `tools/eq_refs.py`
→ 0). Drafts below write equations in Unicode for readability; builders set each in TeX (`$…$`, backslashes doubled in JS
strings). **Walkthrough texts ≤ 30 words.** Colours as in convention 9. A view hidden on portrait phones never carries a
step's key number (it is repeated in a visible title or readout). Mode-dependent controls are hidden per mode with a
class toggled by the chapter script (the ch05 `.xxx-off` pattern). Γ-convention chips appear wherever a clockwise Γ is
used, and the status line shows both conventions (the ch01 two-convention pattern). Complex arithmetic uses a small local
helper (`cadd, cmul, cdiv, clog, csqrt`, or `Viz.num.C` where it suffices). Parallel builders use private scratch
subfolders. Each explainer fits 360×640 … 1920×1080 and the 1000×700 notebook frame with no scrolling (fit plan per
explainer).

### E1 · superposition_sandbox
- **Title:** "Where does the body come from?" · **Summary:** "Add a stream, sources, sinks, a doublet and a vortex; the
  streamline through the stagnation point becomes a wall — open when the sources win, closed when they cancel." ·
  **CORE:** C03, C04, C05 (also N16 (6.16), N17 (6.17), R05 (6.18), N24, N25 (6.28), N27 (6.30), N28 (6.32), N29;
  closed bodies preview C06) · **Reference:** `fid_formula_lab.html` (a sum whose terms are clickable bars with a
  total) with `angular_frequency_explorer_1.html`'s presets and "Right now" notes.
- **meta:** `viz:order 1` · `viz:sections 6.2 6.3` · `viz:equations 6.7 6.8 6.15 6.16 6.28 6.29 6.30 6.31 6.32` ·
  `viz:fluidpy ch06.superposition_state ch06.half_body_numbers ch06.rankine_oval ch06.doublet_limit_error
  ch06.flow_from_spec` · `viz:derivations D05 D06 D07`.
- **Physics (JS ↔ Python):** `elements(s)` builds the element list of the current kit ↔ `ch06.flow_from_spec(spec)`;
  `wdw(els, z)` returns (w, dw/dz) by complex sums (uniform (U)z; source (m/2π)ln(z − z₀); vortex −(iΓ/2π)ln(z − z₀);
  doublet −D/(2π(z − z₀)) with D = d_x + i d_y) ↔ `pf.Flow.w`, `.dwdz`; `state(els, x, y)` ↔
  `ch06.superposition_state(spec, x, y)` (per-element u, v, total, ψ, φ, C_p); `stagnation(els)` (Newton on dw/dz from
  seeds on the axis and the body, deflation, residual < 1e-12) ↔ `pf.Flow.stagnation_points`; `halfBody(U, m)` ↔
  `ch06.half_body_numbers`; `oval(U, m, a)` ↔ `ch06.rankine_oval` (half-length a√(1 + m/πUa), half-width by `brentq`
  on Uh = (m/π)tan⁻¹(a/h)); `pairError(eps, d)` ↔ `ch06.doublet_limit_error`. ψ uses θ ∈ [0, 2π) for a single source
  (convention 5).
- **Modes (kit chips, `display: 'chips'`):** "stream + source" (half-body) · "+ sink" (Rankine oval, source at −ε, sink at
  +ε) · "source–sink pair" (no stream; the ε-squeeze) · "stream + doublet" (cylinder) · "+ vortex" (lifting cylinder, Γ
  counterclockwise in this explainer, stated) · "custom" (each element toggled).
- **Views** (rows [1.35, 1]):
  1. `flow` "The flow" (equal aspect, x ∈ [−4, 6] m, y ∈ [−3, 3] m) — faint speed map (`Viz.field.heatmap`, or C_p when
     the `map` chip says so), ψ contours at equal Δψ (blue), the dividing streamline(s) bold black, stagnation points as
     black dots labelled S, element markers (teal source ●, rose sink ○, purple doublet ⊙, amber vortex ↺), tracer
     particles advected by RK4 (`Viz.field.particles`), the probe ✚ with its velocity arrow (black) and the element
     arrows (their colours) tip to tail. Pointer: click → probe; drag a source/sink marker along the axis (sets `eps`).
  2. `terms` "Velocity at the probe" — term bars u (and v below) per element, coloured, plus the total (black); click a
     bar to isolate that element's arrow in `flow`.
  3. `cp` "C_p along the body" (`hidePortrait: true`) — C_p = 1 − \|u\|²/U² along the dividing streamline against arc
     length from S (upper half), the zero line, the C_p = 1 nose, and for the half-body the 113.2° zero marked.
  Portrait: `flow` + `terms`; the key numbers (a, h_max or half-length) sit in the `flow` title.
- **Controls (≤ 5 visible):** `kit` (chips above) · `U` "Stream $U$" 0…3 m/s, step 0.05, default 1, help "the uniform flow
  from the left" · `m` "Source strength $m$" 0…20 m²/s, step 0.1, default 6.28, help "volume per second per metre of
  depth" · `eps` "Half-gap $\varepsilon$" 0.02…3 m, log, default 1, help "source at −ε, sink at +ε" (oval, pair) · `G`
  "Vortex $\Gamma$ (ccw)" −20…20 m²/s, default 6.28 (+ vortex, custom; optional) · `d` "Doublet $\lvert d\rvert$" 0…40 m³/s,
  default 2πUa² with a = 1 (doublet kits; optional, shown in place of `m`) · `map` chips "speed · C_p" (optional).
- **Transport:** `t` 0 → 20 s, rate 1, `end: 'loop'` — advects 150 tracers seeded on the left edge; in the pair kit the
  walkthrough's squeeze step instead animates `eps` from 1 to 0.05 with 2mε fixed (its `enter()` runs a local ticker;
  said in the step).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **presets** "stream only" {kit: custom, all
  off}, "half-body" {kit: 'stream + source', U: 1, m: 6.283}, "Rankine oval" {kit: '+ sink', U: 1, m: 6.283, eps: 1},
  "source–sink → doublet" {kit: 'source–sink pair', m: 1, eps: 1} (then ▶ squeeze), "cylinder" {kit: 'stream + doublet',
  U: 1, d: 6.283}, "lifting cylinder" {kit: '+ vortex', U: 1, d: 6.283, G: −6.283} · **status** ("🟢 closed body: Σm = 0 ·
  half-length 1.732 m" / "🔵 open body: net source 6.28 m²/s ⇒ it widens to 2h_max = 6.28 m" / "⚪ no body: nothing stops
  the stream") · **terms** (per-element u, v) · **inspector** (click: "u = U + (m/2π)(x − x_s)/r_s² − … = 1 + 0.9549·(…)
  = **0.382** m/s" line by line) · **transport** (tracers).
- **Readouts:** "Stagnation $x_S$" (m) · "Width $2h_{\max}$" (m) or "Half-length" (m) · "$C_p$ at probe" (–) · "Σ sources"
  (m²/s).
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**The flow**: blue lines are streamlines at equal steps of ψ (closer = faster); the bold
     black line is the one through the stagnation point S; colours mark the elements (teal source, rose sink, purple
     doublet, amber vortex, blue stream). **Velocity at the probe** splits the velocity at ✚ into one bar per element.
     **C_p along the body** (hidden on phones) follows the pressure over the surface."
  1. *Each element at the probe* — one `Viz.work.line` per active element: stream "u = U = 1.00 m/s"; source
     "(m/2π)(x − x_s)/r² = (6.28/6.28)(0.5 − 0)/(0.5² + 1.5²) = 0.200 m/s", v likewise; … — "Laplace is linear, so each
     element's velocity is computed as if it were alone" (C03).
  2. *The sum* — "u = Σu_i = **live u**, v = **live v**; speed **live q** m/s; ψ = **live psi** m²/s" (boxed).
  3. *Where the fluid stops* — half-body: "on the axis v = 0, so solve u = U + m/2πx = 0: x = −m/2πU = −6.28/(6.28 × 1) =
     **−1.00 m**"; oval: "U + (m/2π)[1/(x + ε) − 1/(x − ε)] = 0 ⇒ x = ±ε√(1 + m/πUε) = **±1.732 m**"; cylinder: "±a"; with
     a vortex: the two roots of (6.38) in project sign. Boxed.
  4. *The body* — "ψ at S = **live psiS** m²/s; every point with that ψ is on the black line, and by (6.16) ∂ψ/∂s = 0 there
     is a wall." Half-body: "far downstream the width is 2h_max = m/U = **6.28 m** — all the source's fluid, carried at
     U". Oval: "half-width **1.307 m** (solve Uh = (m/π)tan⁻¹(ε/h))". "Σm = **live sum** m²/s: the body closes only if
     it is 0."
  5. *Pressure at the probe* — "(6.32) C_p = 1 − \|u\|²/U² = 1 − (**live q**/1.00)² = **live cp**"; "R05: the same
     Bernoulli constant on every streamline, because the flow is irrotational."
  6. *C_p along the body* — hint on phones "turn the phone sideways for the surface-pressure view"; else "C_p = 1 at S,
     falls to **live cpmin** on the shoulder, returns to 0 far downstream (half-body: crosses 0 at 113.2° from +x)."
  7. *Right now* — "tracers released **live t** s ago; the fastest are the ones hugging the shoulders."
  8. *Reading the current setting* — half-body: "An open body: the source keeps adding m m²/s, so the body widens to
     m/U and never closes. Halve U and every length doubles." · oval: "The sink swallows what the source emits: a
     closed body, symmetric front and back — its drag is zero (C06)." · pair (no stream): "No stream, no body: the flow
     loops from source to sink. Shrink ε with 2mε fixed and the loops become the doublet's circles (error ∝ ε² = **live
     err**)." · cylinder: "A doublet in a stream: the circle r = a is the ψ = 0 streamline." · lifting cylinder: "The
     vortex breaks the top–bottom symmetry: the stagnation points slide (to the bottom for a clockwise vortex) — lift
     without drag (C07)." · no body: "Nothing stops the stream; every streamline is straight."
- **Derivation tab:** **D05** (5 steps) `view: 'flow'`; goal page `set {kit: 'stream + source', U: 1, m: 6.283}`; step 3
  `set` the probe on the body with `watch` "the black tangent arrow and the teal normal: u·n is 0"; step 5 **live**
  "∂ψ/∂s along the black line = **0.000** m/s at your probe". **D06** (9 steps) `view: 'flow'`, goal page `set {kit:
  'source–sink pair', m: 1, eps: 1}`; step 5 `set {eps: 0.1, m: 10}` with `watch` "the loops tighten toward circles"; step
  9 **live** "at your probe the pair differs from the doublet by **live err** (relative) — about ε²/3". **D07** (10 steps)
  `view: 'flow'`, goal page `set {kit: 'stream + source'}`; step 3 **live** "x = −m/2πU = −(live m)/(2π × live U) = **live
  xS** m"; step 7 `watch` "the width readout approaches m/U"; step 9 `watch` "the teal region's width × U equals m";
  interpret `s => "With your numbers the nose sits ${a} m ahead of the source and the body is ${2hmax} m wide far
  downstream."`.
- **Code:**
  ```python
  spec = [{"kind": "uniform", "U": {{U}}},
          {"kind": "source", "m": {{m}}, "x": 0.0, "y": 0.0}]     # the kit you chose
  s = ch06.superposition_state(spec, {{px}}, {{py}})            # probe ✚
  print(s["u_parts"], s["u"])        # {{uparts}} → u = {{u}} m/s
  print(s["stagnation"])             # {{stag}}  (u = v = 0)
  print(s["psi_dividing"], s["closed"])   # ψ on the body {{psiS}}; closed? {{closed}}
  print(ch06.half_body_numbers({{U}}, {{m}}))   # a = m/2πU, h_max = m/2U (6.31)
  ```
- **Walkthrough (6 steps):** 1. "A body from nothing" — "Dye bends round a stone nobody drew. Can adding flows produce a
  body?" `set {kit: 'custom', all off}` · 2. "Add a source" — "The source pushes fluid out; the stream carries it off. A
  stagnation point S appears." `set` half-body, `readouts: ['xS']`, highlight `readout:xS` · 3. "The wall" — "The bold
  streamline through S never lets fluid cross: call it a wall (6.16), ∂ψ/∂s = 0." `derive: {id: 'D05', step: 5}` · 4.
  "Why the sizes" — "Stagnation at x = −m/2πU; width m/U far away (6.31)." `terms: true`, `derive: {id: 'D07', step: 3}`,
  `code: {lines: [5, 6]}` · 5. "Close it" — "Add an equal sink: the body closes (Σm = 0). Squeeze the pair: a doublet."
  `set` Rankine oval, then `play` in pair kit · 6. "Your turn" — "Predict the half-body's width for m = 10 m²/s, U = 2 m/s,
  then check." `controls: ['U', 'm']`.
- **Equations:** `uni` "Uniform flow" ref 'Eq. (6.7)' `\psi=-Vx+Uy` · `src` "Source" ref 'Eq. (6.15)'
  `\phi=\frac{m}{2\pi}\ln r` live "u_r = m/2πr = … m/s at your probe" · `dbl` "Doublet" ref 'Eq. (6.29)'
  `\phi=-\frac{\mathbf d\cdot\mathbf x}{2\pi r^2}=\frac{\lvert\mathbf d\rvert}{2\pi}\frac{\cos\theta}{r}` · `hb` "Half-body" ref 'Eq.
  (6.31)' `\psi=Ur\sin\theta+\frac{m}{2\pi}\theta` live "a = … m, h_max = … m" · `cp` "Pressure coefficient" ref 'Eq.
  (6.32)' `C_p=1-\lvert\mathbf u\rvert^2/U^2` · `wall` "No flow through the body" ref 'Eq. (6.16)'
  `\partial\psi/\partial s=0` · symbols U, m, d, Γ, ψ, φ, C_p with units.
- **Check yourself:** (1) "Double m at fixed U. How far does the nose move, and how wide does the body get?" — "Both
  double: a = m/2πU and 2h_max = m/U are proportional to m." `set {m: 12.57}` · (2) "In the oval preset, make the sink
  slightly weaker than the source (custom kit). What happens to the body?" — "It opens again: the leftover source
  strength must leave, so the body stretches downstream without closing." · (3) "Click a probe on the bold line. Why is
  C_p there not 1?" — "Only S is a stagnation point; elsewhere on the body the fluid slides along it (ideal flow keeps
  no no-slip)." · (4) "Squeeze the pair to ε = 0.05 m. By how much does it differ from the doublet?" — "About ε²/3 ≈
  8×10⁻⁴ relative: the error falls like ε² (D06 step 9)." `set {kit: 'source–sink pair', eps: 0.05}`.
- **Selftest parity rows:** `{name: 'half-body a', js: halfBody(1, 2*Math.PI).a, py:
  'ch06.half_body_numbers(1.0, 6.283185307179586)["a"]', rtol: 1e-12}` · `{name: 'half-body h_max', js: halfBody(2,
  10).h_max, py: 'ch06.half_body_numbers(2.0, 10.0)["h_max"]', rtol: 1e-12}` · `{name: 'probe u', js: state(HB, 0.5,
  1.5).u, py: 'ch06.superposition_state([{"kind": "uniform", "U": 1.0}, {"kind": "source", "m": 6.283185307179586, "x":
  0.0, "y": 0.0}], 0.5, 1.5)["u"]', rtol: 1e-12}` · `{name: 'probe psi', js: state(HB, 0.5, 1.5).psi, py: 'ch06.
  superposition_state([{"kind": "uniform", "U": 1.0}, {"kind": "source", "m": 6.283185307179586, "x": 0.0, "y": 0.0}],
  0.5, 1.5)["psi"]', rtol: 1e-12}` · `{name: 'oval half-length', js: oval(1, 2*Math.PI, 1).half_length, py:
  'ch06.rankine_oval(1.0, 6.283185307179586, 1.0)["half_length"]', rtol: 1e-12}` · `{name: 'oval half-width', js:
  oval(1, 2*Math.PI, 1).half_width, py: 'ch06.rankine_oval(1.0, 6.283185307179586, 1.0)["half_width"]', rtol: 1e-9}` ·
  `{name: 'pair error eps=0.1', js: pairError(0.1, 2), py: 'ch06.doublet_limit_error([0.1], 2.0)[0]', rtol: 1e-9}` ·
  invariant `{name: 'cp at S', js: state(HB, -1, 0).cp, expect: 1, atol: 1e-12}`.
- **Fit plan:** 360×640: `flow` (60 %) + `terms` (bars in two rows u, v); kit chips wrap to two rows; presets hidden on
  short portrait screens (scoped CSS, Derivation tab exempt); status one line ("🟢 closed · 1.732 m"). Landscape phone:
  `flow` + `terms` side by side. Notebook 1000×700 and desktop: rows [1.35, 1], `terms` and `cp` share the second row.

### E2 · cylinder_circulation_lift
- **Title:** "How does spin turn into lift?" · **Summary:** "Drag the circulation round a cylinder: the stagnation points
  slide down, meet at Γ = 4πaU and leave the body; the pressure difference integrates to exactly ρUΓ while the drag stays
  zero." · **CORE:** C06, C07 (also R12 (6.33), R13 (6.34), N30 (6.35), N31, N32–N36 (6.36)–(6.39), N38; links ch03
  `galilean_frames_cylinder`) · **Reference:** `amplitude_phase_second_order_II_3.html` (system + graphs linked by one
  state; a numbered derivation with live numbers in the explanation).
- **meta:** `viz:order 2` · `viz:sections 6.3` · `viz:equations 6.33 6.34 6.35 6.36 6.37 6.38 6.39 6.40` · `viz:fluidpy
  ch06.cylinder_circulation_state ch06.cylinder_surface_cp ch06.cylinder_surface_pressure ch06.lift_per_span
  ch06.surface_pressure_force` · `viz:derivations D09 D10 D11`.
- **Physics:** `cyl(U, a, Gcw)` gives w, dw/dz (w = U(z + a²/z) + (iΓ_cw/2π)ln(z/a)) ↔ `pf.cylinder(U, a, Gamma_cw=…)`;
  `stag(U, a, Gcw)` closed form (6.38) and the off-body root ↔ `ch06.cylinder_stagnation_points`; `surfP(th, s)` (6.39)
  ↔ `ch06.cylinder_surface_pressure`; `forces(s, n)` periodic trapezoid of −∮p n a dθ ↔ `ch06.surface_pressure_force`;
  `termsL(s)` the lift from each of the four terms of the expanded (6.39) (constant, −2ρU² sin²θ, −ρUΓ sin θ/πa,
  −ρΓ²/8π²a²) — each by the same quadrature; `regime(s)` ↔ `ch06.cylinder_circulation_state(...)["regime"]`.
- **Views** (rows [1.3, 1]):
  1. `flow` "Round the cylinder" (equal aspect, x, y ∈ [−4a, 4a]) — streamlines (blue), the cylinder, stagnation points
     (black dots; the free point with its crossing streamlines), surface arrows −(p − p∞)n (orange push-in where
     C_p > 0, blue pull-out where C_p < 0), a curved arrow showing the circulation's sense with its label in the chosen
     convention. Title carries "L = … N/m · D = 0". Pointer: click the surface → inspector at that angle.
  2. `cp` "Surface C_p(θ)" — C_p against θ ∈ [0°, 360°]: current (black), Γ = 0 ghost (grey dashed), the stagnation
     angles as vertical markers, and in D11 step 7 the cross term's contribution shaded amber; optional qualitative
     real band (grey, "qualitative sketch").
  3. `forces` "Forces per metre" (`hidePortrait: true`) — term bars: lift contributions of the four pressure terms
     (constant, sin², cross, Γ²; only "cross" non-zero, amber) and the total lift with a ◇ at ρUΓ; the drag bar (muted)
     at 0.
- **Controls:** `gam` "Circulation $\Gamma/4\pi aU$" −2…2, step 0.01, default 0.159, help "the book's clockwise Γ in units of
  the critical value" · `U` "Stream $U$" 1…30 m/s, default 10 · `a` "Radius $a$" 0.02…0.5 m, default 0.1 (optional) ·
  `conv` chips "book Γ (clockwise) · project Γ (ccw)" (display only: flips the sign shown, not the flow) · `band` toggle
  "real flow (qualitative)" (optional) · `frame` chips "body · fluid" (optional; fluid frame subtracts U: a doublet + vortex).
- **Transport:** `gam` 0 → 2, rate 0.25 /s, `end: 'hold'`, hold 1.5 s; end-of-run `Viz.card`: "Γ went 0 → 8πaU: two
  stagnation points slid down, merged at Γ = 4πaU and left; L grew linearly to 8πρaU²; D stayed 0."
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** "d'Alembert: Γ = 0"
  {gam: 0}, "Γ = 2πaU" {gam: 0.5}, "merging: Γ = 4πaU" {gam: 1}, "free point: Γ = 6πaU" {gam: 1.5}, "lift down" {gam:
  −0.5}, "moving cylinder" {gam: 0, frame: 'fluid'} · **status** ("two surface stagnation points · sin θ = −0.159 ·
  Γ_cw = 2.00 ⇔ Γ_ccw = −2.00 m²/s" / "merged at the bottom" / "free stagnation point at r = 0.262 m") · **terms** (lift
  by pressure term) · **inspector** (click at θ: "u_θ = −2U sin θ − Γ/2πa = −20 × 0.707 − 3.18 = −17.32 m/s; C_p = 1 −
  (17.32/10)² = −2.00; p − p∞ = C_p × ½ρU² = −120 Pa") · **notes** (regime text).
- **Readouts:** "Stagnation θ₁, θ₂" (°) · "Lift $L$" (N/m) · "Drag $D$" (N/m) · "Top / bottom speed" (m/s).
- **Explain:**
  0. *What the views show* — "**Round the cylinder**: streamlines, the stagnation points (black), and the push of the
     pressure on the surface — <b class="c-orange">orange</b> arrows push in where the pressure is above p∞, <b
     class="c-blue">blue</b> pull out where it is below. **Surface C_p(θ)**: black now, grey dashed without circulation.
     **Forces** (hidden on phones): which part of the pressure makes the lift."
  1. *The circulation, both ways* — "Γ/4πaU = **live gam** ⇒ Γ_cw = **live G** m²/s (book, clockwise) ⇔ Γ_ccw = **−live G**
     (project). The flow's circulation measured counterclockwise is −Γ_cw."
  2. *Surface speed (6.37)* — "u_θ(a, θ) = −2U sin θ − Γ/2πa: top (θ = 90°) −20 − **live dG** = **live top** m/s; bottom
     **live bot** m/s — the vortex adds on top and subtracts below."
  3. *Stagnation points (6.38)* — "sin θ = −Γ/4πaU = **−live gam** ⇒ θ = **live th1**° and **live th2**°" (two roots,
     both below the axis) or "Γ = 4πaU: one point at −90°" or "Γ > 4πaU: on θ = −90°, U(1 + a²/r²) − Γ/2πr = 0 ⇒ r =
     [Γ + √(Γ² − (4πaU)²)]/4πU = **live rf** m (its partner inside at a²/r)". Boxed.
  4. *Pressure at the clicked point (6.39)* — "p − p∞ = ½ρ[U² − (−2U sin θ − Γ/2πa)²] = ½ × 1.2 × [100 − (**live
     uth**)²] = **live dp** Pa".
  5. *The lift integral* — "L = −∫₀^{2π} p sin θ a dθ. Expanding (6.39): only −ρUΓ sin θ/πa survives against sin θ;
     ∫ sin²θ dθ = π gives L = ρUΓ = 1.2 × **live U** × **live G** = **live L** N/m" (boxed); "our quadrature on 128
     points: **live Lq** N/m. Drag: every term against cos θ integrates to 0 ⇒ D = **live D**."
  6. *Forces view* — hint on phones; else "the amber bar is the cross term; the other three are zero — as in D11."
  7. *Right now* — "transport at Γ/4πaU = **live gam**."
  8. *Reading the current setting* — Γ = 0: "d'Alembert: symmetric front/back and top/bottom — no force at all, despite
     ±180 Pa on the surface." · 0 < Γ < 4πaU: "Two stagnation points below the axis; faster on top ⇒ lower pressure on
     top ⇒ lift **live L** N/m upward; still no drag." · Γ = 4πaU: "The points merge at the bottom; the cylinder carries
     a closed ring of fluid round it." · Γ > 4πaU: "A free stagnation point below the body; a ring of fluid circulates
     round the cylinder for ever (ideal flow)." · Γ < 0: "Counterclockwise circulation (in the book's sign, negative):
     faster below, lift downward." · fluid frame: "Seen from the still fluid: a doublet (+ vortex) moving left — the same
     forces (Galilean invariance, Ch. 3 §3.3)."
- **Derivation tab:** **D09** (6 steps) `view: 'flow'`, goal `set {gam: 0}`; step 2 **live** "D = 2πUa² = 2π × **live U** × **live a**² =
  **live d** m³/s"; step 6 `watch` "C_p = −3 at θ = ±90°". **D10** (9 steps) `view: 'cp'` (phones keep `cp`), goal `set
  {gam: 0.159}`; step 5 **live** "sin θ = −**live gam**"; step 7 `set {gam: 1}`, `watch` "the two markers merge at −90°";
  step 9 `set {gam: 1.5}`, **live** "r₊ = **live rf** m, r₋ = **live rm** m, r₊r₋ = a²". **D11** (9 steps) `view: 'forces'`
  (phones: `cp`), goal `set {gam: 0.5}`; step 7 highlights `term:cross`, **live** "ρUΓ/π × π = **live L** N/m"; step 9
  `watch` "the drag bar stays at 0". Interpret `s => "Your cylinder: L = ρUΓ = ${L} N/m; stagnation points at
  ${th1}° and ${th2}°."`.
- **Code:**
  ```python
  st = ch06.cylinder_circulation_state({{U}}, {{a}}, Gamma_cw={{G}}, rho=1.2)
  print(st["theta1_deg"], st["theta2_deg"], st["r_free"])   # {{th1}} {{th2}} {{rf}}
  print(st["regime"])                    # "{{regime}}"
  print(st["L"], st["L_KJ"], st["D"])    # {{Lq}} = ρUΓ = {{L}} N/m, D = {{D}}
  th = np.radians({{thdeg}})             # your clicked angle
  print(ch06.cylinder_surface_pressure(th, {{U}}, {{a}}, Gamma_cw={{G}}, rho=1.2))   # {{dp}} Pa, (6.39)
  ```
- **Walkthrough (6 steps):** 1. "A spinning body" — "A spinning ball curves. Where does a sideways force come from?"
  `set {gam: 0}` · 2. "No spin, no force" — "Front and back pressures mirror each other: D = 0, L = 0 — d'Alembert."
  `derive: {id: 'D09', step: 6}`, `readouts: ['L', 'D']` · 3. "Add circulation" — "Press ▶: the stagnation points slide
  down; sin θ = −Γ/4πaU (6.38)." `play: true`, `controls: ['gam']` · 4. "Why lift" — "Faster on top, lower pressure on top.
  Only the cross term of (6.39) survives: L = ρUΓ." `terms: true`, `derive: {id: 'D11', step: 7}`, `code: {lines: [4, 4]}`
  · 5. "Which Γ?" — "Every Γ gives a valid flow. Nature picks it at a sharp edge (Ch. 14)." `set {gam: 1.5}`, `notes:
  true` · 6. "Your turn" — "Predict L for U = 20 m/s at the same Γ, then check." `controls: ['U']`.
- **Equations:** `cyl` ref 'Eq. (6.33)' `\psi=U(r-a^2/r)\sin\theta` · `cp0` ref 'Eq. (6.35)' `C_p=1-4\sin^2\theta` · `cylG` ref
  'Eq. (6.36)' `\psi=U\big(r-\frac{a^2}{r}\big)\sin\theta+\frac{\Gamma}{2\pi}\ln\frac ra` · `us` ref 'Eq. (6.37)'
  `u_\theta(a,\theta)=-2U\sin\theta-\Gamma/2\pi a` · `stag` ref 'Eq. (6.38)' `\sin\theta=-\Gamma/4\pi aU` live · `p` ref 'Eq.
  (6.39)' `p(a,\theta)=p_\infty+\tfrac12\rho\big[U^2-\big(-2U\sin\theta-\frac{\Gamma}{2\pi a}\big)^2\big]` · `L` ref 'Eq. (6.40)'
  `L=\rho U\Gamma` live "= 1.2 × 10 × 2 = 24 N/m".
- **Check yourself:** (1) "Double U at fixed Γ. What happens to the stagnation points and to L?" — "They move back up
  (Γ/4πaU halves) while L doubles (ρUΓ)." `set {U: 20}` · (2) "Set Γ = 4πaU and switch the convention chips. Which
  numbers change?" — "Only the sign shown for Γ; the flow, the stagnation point and L are the same object." · (3)
  "With Γ ≠ 0, why is the drag still zero?" — "The circulation changes top vs bottom, not front vs back: the pressure is
  still fore–aft symmetric, so ∫p cos θ = 0." · (4) "Turn on the real-flow band at Γ = 0. Where do ideal and real differ,
  and what would that do to the drag?" — "On the back half: the real pressure stays low behind the separation point,
  so the front push wins — drag (Ch. 9)." `set {gam: 0, band: true}`.
- **Selftest parity rows:** `{name: 'theta1', js: stag(10, 0.1, 2).th1deg, py:
  'ch06.cylinder_circulation_state(10.0, 0.1, Gamma_cw=2.0)["theta1_deg"]', rtol: 1e-10}` · `{name: 'lift quadrature',
  js: forces(S0, 128).L, py: 'ch06.cylinder_circulation_state(10.0, 0.1, Gamma_cw=2.0)["L"]', rtol: 1e-10}` · `{name:
  'L_KJ', js: 1.2*10*2, py: 'ch06.lift_per_span(1.2, 10.0, Gamma_cw=2.0)', rtol: 1e-14}` · `{name: 'free point', js:
  stag(10, 0.1, 1.5*4*Math.PI*0.1*10).r, py: 'ch06.cylinder_circulation_state(10.0, 0.1,
  Gamma_cw=18.84955592153876)["r_free"]', rtol: 1e-10}` · `{name: 'cp at 90', js: cpSurf(Math.PI/2, 10, 0.1, 0), py:
  'ch06.cylinder_surface_cp(1.5707963267948966, 10.0, 0.1, Gamma_cw=0.0)', rtol: 1e-12}` · `{name: 'p at 45deg', js:
  surfP(Math.PI/4, S0), py: 'ch06.cylinder_surface_pressure(0.7853981633974483, 10.0, 0.1, Gamma_cw=2.0, rho=1.2)', rtol:
  1e-12}` · exact-text `{name: 'regime word', js: regime(S15), expect: 'free stagnation point'}` · invariant `{name: 'drag
  zero', js: forces(S0, 128).D, expect: 0, atol: 1e-9}`.
- **Fit plan:** 360×640: `flow` (55 %) + `cp`; status one line ("L 24.0 N/m · θ −9.2°, −170.8°"); presets hidden on short
  portrait. Landscape phone: `flow` + `cp` side by side, `forces` hidden. Desktop/notebook: `forces` beside `cp`.

### E3 · vortex_wall_images
- **Title:** "What does a wall feel as an eddy passes?" · **Summary:** "A clockwise vortex drifts along a wall beside its
  image; a sensor on the wall records suction as it passes and a small over-pressure before and after — from unsteady
  Bernoulli." · **CORE:** C08 (also N39, N40 (6.41), R16, R17, R32 (6.99); links ch05 `point_vortex_lab`, ch04
  `which_bernoulli`) · **Reference:** `forced_damped_vibrations.html` (one time slider driving the system and its x(t)
  graph, with the explanation panel computing the value at the current time).
- **meta:** `viz:order 3` · `viz:sections 6.3` · `viz:equations 6.41 6.99` · `viz:fluidpy ch06.example_6_1
  ch06.example_6_1_wall_pressure ch06.two_sources ch06.two_source_streamline` · `viz:derivations D12 D13`.
- **Physics:** `xi(t, s)` = (h, Γt/4πh) ↔ `ch06.example_6_1(t)["xi_y"]`; `pOrigin(t, s)` closed form + its parts
  ↔ `ch06.example_6_1(t)["p_origin"]`, `["unsteady"]`, `["speed_part"]`; `pWall(y, t, s)` ↔
  `ch06.example_6_1_wall_pressure(y, t, …, split=True)` (φ(0, y, t) differentiated in t analytically, v(0, y, t) from
  the image pair); `psiVortex(x, y, t)` for the streamlines; source mode: `twoSources(m, a)` ↔ `ch06.two_sources`, the
  (6.41) curve ↔ `ch06.two_source_streamline`.
- **Modes:** "vortex by a wall" · "source by a wall (6.41)".
- **Views** (rows [1.2, 1]):
  1. `scene` "The eddy and its image" (equal aspect, x ∈ [−3h, 3h], y ∈ [−4h, 4h]) — the wall x = 0 as a thick line
     whose points are coloured by p − p∞ (orange > 0, blue < 0), the vortex (amber ↻) and its faded image (↺), their
     paths (dotted), instantaneous streamlines, the sensor at the origin (⬤) and a movable wall probe; source mode: two
     sources, ψ contours and the dashed (6.41) curves. Title carries "p(0, 0, t) − p∞ = … Pa".
  2. `trace` "Sensor pressure p(0, 0, t) − p∞" — t ∈ [−60, 60] s (scaled by h²/Γ): the faint whole curve, the bold part
     so far, the moving dot, markers at t = 0, ±4πh²/Γ (zero) and ±4√3πh²/Γ (maximum); dashed orange (unsteady part)
     and blue (speed part).
  3. `wall` "Along the wall now" (`hidePortrait: true`) — p − p∞ against y on the wall at the current t, with the
     vortex's height marked.
- **Controls:** `G` "Vortex strength $\Gamma$" 0.1…5 m²/s, default 1 (vortex mode) · `h` "Distance $h$" 0.2…3 m, default 1 ·
  `yw` "Wall probe $y$" −5…5 m, default 0 (optional; also by clicking the wall) · `m` "Source $m$" 0.5…10 m²/s, default
  6.28 (source mode) · `mode` chips.
- **Transport:** `t` −60 → 60 s (the builder scales the range by h²/Γ and says so), rate 5 s/s, `end: 'hold'`; end-of-run
  card: "The eddy passed: suction −25.3 Pa at t = 0, zero at ±12.6 s, +3.17 Pa at ±21.8 s."
- **Depth features:** Explain + Code + Derivation · **transport** · **linked views** (3) · **modes** · **presets**
  "strongest suction: t = 0" {t: 0}, "crossing: t = 4πh²/Γ" {t: 12.566}, "largest push: t = 4√3πh²/Γ" {t: 21.766},
  "closer eddy (h = 0.5)" {h: 0.5}, "source by a wall" {mode: 'source'}, "two sources" {mode: 'source', showBoth: true}
  · **terms** (unsteady vs speed part at the probe, and their sum) · **inspector** (click the wall: the arithmetic of
  ∂φ/∂t and v there) · **status** ("🌀 approaching: over-pressure" / "🔵 passing: suction" / "🌀 leaving: over-pressure").
- **Readouts:** "Drift speed" (m/s) · "Vortex height $\xi_y$" (m) · "$p-p_\infty$ at sensor" (Pa) · "Time" (s).
- **Explain:**
  0. *What the views show* — "The vortex (<b class="c-amber">amber</b>, clockwise) and its image (faded, opposite sign)
     sit at mirror points; the wall is coloured by p − p∞ (<b class="c-orange">orange</b> above p∞, <b
     class="c-blue">blue</b> below). The trace is the sensor's record; dashed curves are its two parts."
  1. *Why an image* — "An opposite image at (−h, ξ_y) makes u = 0 on the wall: the two vortices' x-velocities cancel
     there (D12 and D13 step 1)."
  2. *The drift* — "The image pushes the vortex along the wall at Γ/4πh = **live G**/(4π × **live h**) = **live drift**
     m/s; its own velocity is taken as zero (a finite core spins symmetrically)."
  3. *Position now* — "ξ = (h, Γt/4πh) = (**live h**, **live xiy**) m."
  4. *The two pressure parts at the probe* — "unsteady: −ρ∂φ/∂t = ρΓ²/[4π²(h² + s²)] = **live pu** Pa (with s = ξ_y −
     y_probe); speed: −½ρv², v = Γh/[π(h² + s²)] = **live v** m/s ⇒ **live ps** Pa" (for the origin; the wall-probe
     version uses the general formula and says so).
  5. *The sum* — "p − p∞ = **live pu** + (**live ps**) = **live p** Pa" (boxed); "(Example 6.1) =
     (ρΓ²/4π²)(s² − h²)/(s² + h²)²".
  6. *Special times* — "zero when s = h: t = 4πh²/Γ = **live t0** s; largest push when s = √3h: t = **live t1** s,
     p − p∞ = ρΓ²/32π²h² = **live pmax** Pa."
  7. *The wall view* — hint on phones; else "the suction patch travels with the eddy, the pushes lead and trail it."
  8. *Reading the current setting* — \|s\| < h: "Passing: the speed part wins — suction, deepest when the eddy is level
     with the sensor." · \|s\| > h: "Far away: the unsteady part wins (it decays like 1/s², the speed part like 1/s⁴) —
     a small over-pressure." · source mode: "A source's image keeps its sign: the wall still has no flow through it
     (∂φ/∂y = 0 by symmetry); the origin is a stagnation point (6.41)."
- **Derivation tab:** **D12** (5 steps) `view: 'scene'`, goal `set {mode: 'vortex', t: 0}`; step 1 `watch` "flip the
  image's sign in your head: with the same sign, the wall would leak" (a toggle `wrongImage` in the step's `set` shows
  the leaking wall in rose); step 4 `set {mode: 'source'}`. **D13** (11 steps) `view: 'trace'` (phones keep `trace`),
  goal `set {mode: 'vortex', t: 0}`; step 3 **live** "ξ_y = Γt/4πh = **live xiy** m"; step 8 **live** "∂φ/∂t at the
  origin = −Γ²/4π²(h² + s²) = **live dphidt** m²/s²"; step 11 **live** "p − p∞ = **live p** Pa"; `watch` on step 11 "the
  dot on the trace". Interpret `s => "At t = ${t} s the sensor reads ${p} Pa: ${regimeText}."`.
- **Code:**
  ```python
  r = ch06.example_6_1({{t}}, Gamma={{G}}, h={{h}}, rho=1000.0)
  print(r["xi_y"])                        # vortex height {{xiy}} m (drift Γ/4πh)
  print(r["unsteady"], r["speed_part"])   # {{pu}} + {{ps}} Pa
  print(r["p_origin"])                    # = {{p}} Pa at the sensor
  pw = ch06.example_6_1_wall_pressure({{yw}}, {{t}}, Gamma={{G}}, h={{h}})
  print(pw)                               # {{pwall}} Pa at your wall probe
  ```
- **Walkthrough (6 steps):** 1. "An eddy by a wall" — "An eddy drifts along a harbour wall. What does a pressure sensor
  in the wall record?" `set {t: -40}` · 2. "The mirror trick" — "An opposite vortex behind the wall cancels the flow
  through it (D12)." `derive: {id: 'D12', step: 2}` · 3. "It drifts" — "The image pushes the eddy along at Γ/4πh."
  `play: true`, `readouts: ['drift']` · 4. "Two pressure parts" — "Speed part always suction; the changing flow adds
  −ρ∂φ/∂t > 0 (6.99)." `terms: true`, `derive: {id: 'D13', step: 10}`, `code: {lines: [3, 4]}` · 5. "The signal" — "Suction
  at t = 0, zero at 4πh²/Γ, a small push at 4√3πh²/Γ." `set {t: 21.766}` · 6. "Your turn" — "Predict the dip if h halves,
  then drag h." `controls: ['h']`.
- **Equations:** `img` "Image stream function" (no ref; Example 6.1) `\psi=\frac{\Gamma}{2\pi}\big[-\ln r_{\text{image}}+\ln r_{\text{vortex}}\big]`
  · `path` "Drift" `\boldsymbol\xi(t)=(h,\ \Gamma t/4\pi h)` · `bern` "Unsteady Bernoulli" ref 'Eq. (6.99)'
  `\frac{\partial\phi}{\partial t}+\frac12\lvert\nabla\phi\rvert^2+\frac p\rho=\frac{p_\infty}{\rho}` · `ex61` "Sensor pressure" (Example
  6.1) `\frac{p(0,0,t)-p_\infty}{\rho}=\frac{\Gamma^2}{4\pi^2}\frac{(\Gamma t/4\pi h)^2-h^2}{((\Gamma t/4\pi h)^2+h^2)^2}` live · `two` "Two
  sources" ref 'Eq. (6.41)' `x^2-y^2-2xy\cot(2\pi\psi/m)=a^2`.
- **Check yourself:** (1) "Halve h. How do the depth and the duration of the dip change?" — "The dip is 4× deeper
  (Γ²/h²) and 4× shorter (time scale h²/Γ)." `set {h: 0.5}` · (2) "Why is the pressure positive long before and long
  after the eddy passes?" — "Then the unsteady part −ρ∂φ/∂t (∝ 1/s²) beats the speed part (∝ 1/s⁴)." · (3) "Give the
  image the same sign as the vortex (Derivation D12 step 1). What goes wrong?" — "The wall leaks: u ≠ 0 on x = 0; the
  vortex rule must flip the sign." · (4) "In source mode, is the origin a stagnation point?" — "Yes: the two sources' pushes
  cancel there by symmetry."
- **Selftest parity rows:** `{name: 'p origin t=0', js: pOrigin(0, P0).p, py: 'ch06.example_6_1(0.0)["p_origin"]', rtol:
  1e-12}` · `{name: 'p origin t=21.77', js: pOrigin(4*Math.sqrt(3)*Math.PI, P0).p, py:
  'ch06.example_6_1(21.765592370810612)["p_origin"]', rtol: 1e-10}` · `{name: 'unsteady part t=5', js: pOrigin(5,
  P0).pu, py: 'ch06.example_6_1(5.0)["unsteady"]', rtol: 1e-12}` · `{name: 'xi_y t=10', js: xi(10, P0).y, py:
  'ch06.example_6_1(10.0)["xi_y"]', rtol: 1e-12}` · `{name: 'wall p y=0.5 t=3', js: pWall(0.5, 3, P0).total, py:
  'ch06.example_6_1_wall_pressure(0.5, 3.0, split=True)["total"]', rtol: 1e-10}` · `{name: '(6.41) curve', js:
  twoCurve(1.570796, 6.283185, 1, 2.0), py: 'ch06.two_source_streamline(1.570796, 6.283185, 1.0, 2.0)', rtol: 1e-10}`.
- **Fit plan:** 360×640: `scene` (portrait, tall and narrow — the wall vertical) + `trace`; the wall view hidden; status
  one line. Landscape: the three side by side is too tight — `scene` + `trace` in row 1, `wall` below only ≥ 700 px tall.

### E4 · complex_potential_corners
- **Title:** "Why is the velocity u − iv?" · **Summary:** "Click a point: the derivative of w taken along x and along iy
  converge to the same number only for analytic w — that is Cauchy–Riemann — and its mirror image is the velocity; drag
  the exponent and a corner turns from calm to violent." · **CORE:** C09 (also N41–N51 (6.43)–(6.53)) · **Reference:**
  `amplitude_phase_second_order_II_3.html` (a complex number evaluated live, its conjugate and the geometric reading in
  linked views).
- **meta:** `viz:order 4` · `viz:sections 6.4` · `viz:equations 6.42 6.43 6.44 6.45 6.46 6.47 6.48 6.49` · `viz:fluidpy
  ch06.complex_potential_probe ch06.cauchy_riemann_residual ch06.corner_speed_exponent ch06.harmonic_polynomials` ·
  `viz:derivations D14 D15`.
- **Physics:** `wOf(kind, z, p)` (corner Az^n with the cut outside the wedge; source, vortex, doublet, cylinder; `conj`:
  f = z*) ↔ `pf.Corner`, `pf.Source`, …; `probe(kind, x, y, p)` returns φ, ψ, dw/dz, u, v, the two difference quotients
  q_x, q_iy at step h and the Cauchy–Riemann residuals ↔ `ch06.complex_potential_probe(kind, x, y, A, n, h)`;
  `speedExp(n)` ↔ `ch06.corner_speed_exponent`.
- **Views** (rows [1.3, 1]):
  1. `plane` "The z-plane" (equal aspect) — flow net (ψ solid blue, φ dashed grey) from the analytic function (or, for
     z*, the "net" of Re and Im of z* with the label "not a flow"), the walls (black) for the corner, the probe ✚, the
     arrow of dw/dz (purple) and of the velocity u + iv (teal) — mirror images in the horizontal. Pointer: click → probe.
  2. `quot` "Two difference quotients" — against log₁₀ h: Re and Im of q_x (solid) and q_iy (dashed) converging (or not);
     beside them two Cauchy–Riemann bar pairs: ∂φ/∂x vs ∂ψ/∂y and ∂φ/∂y vs −∂ψ/∂x.
  3. `loglog` "Speed near the corner" (`hidePortrait: true`) — \|dw/dz\| against r along the bisector, log–log, slope
     n − 1 marked.
- **Controls:** `kind` chips "corner Azⁿ · source · vortex · doublet · cylinder · z* (not analytic)" · `n` "Exponent $n$
  (α = π/n)" 0.5…4, step 0.01, default 2 (corner) · `A` "Strength $A$" 0.2…3, default 1 (optional) · `lh` "Step
  $\log_{10}h$" −6…−0.3, default −2, help "the size of the complex step" · probe by click.
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **inspector** (the dw/dz arithmetic at the
  probe: "(w(z + h) − w(z))/h = ((2.01 + 2.02i) − …)/0.01 = 2.01 + 2.00i") · **terms** (Cauchy–Riemann pairs as bars) ·
  **presets** "n = 2: stagnation corner (6.24)" {kind: corner, n: 2}, "n = 1: uniform" {n: 1}, "n = ⅔: 270° corner" {n:
  0.6667}, "n = ½: flat plate" {n: 0.5}, "z*: not a flow" {kind: 'conj'} · **status** ("✅ analytic: a flow — dw/dz = u −
  iv" / "⛔ z* is not analytic: the two quotients disagree") · **notes** (corner regime).
- **Readouts:** "$dw/dz$" (1/s) · "$u$, $v$" (m/s) · "Corner α" (°) · "Speed exponent" (–).
- **Explain:**
  0. *What the views show* — "Solid blue: ψ (streamlines); dashed grey: φ; black: the corner's walls. <b
     class="c-accent">Purple</b> arrow: dw/dz; <b class="c-teal">teal</b> arrow: the velocity. The second view compares the
     derivative taken along x and along iy."
  1. *w at the probe* — "w = Azⁿ = 1 × (**live z**)^**live n** = **live w**: φ = **live phi**, ψ = **live psi**."
  2. *Two quotients* — "along x: (w(z + h) − w(z))/h = **live qx**; along iy: (w(z + ih) − w(z))/(ih) = **live qy**;
     exact dw/dz = nAz^{n−1} = **live dw**" (boxed) — "equal ⇒ w is analytic here".
  3. *Velocity* — "(6.45) dw/dz = u − iv ⇒ u = **live u**, v = −Im(dw/dz) = **live v** m/s — the teal arrow is the purple
     one mirrored."
  4. *Cauchy–Riemann (6.44)* — "∂φ/∂x − ∂ψ/∂y = **live cr1**, ∂φ/∂y + ∂ψ/∂x = **live cr2**" (≈ 0; for z*: 2 and 0).
  5. *The corner* — "α = π/n = **live alpha**°; speed ∝ r^{n−1} = r^**live ex** near the tip: **live cornerword**."
  6. *Speed near the corner* — hint on phones; else "the slope of the log–log line is n − 1 = (π − α)/α."
  7. *Reading the current setting* — n > 1: "A corner narrower than a straight wall: the fluid creeps into it and stops
     at the tip — a stagnation point (n = 2 is the 90° corner of (6.24))." · n = 1: "No corner: uniform flow." · ½ ≤ n <
     1: "A corner wider than 180°: the fluid must whip round the tip — infinite speed there (the 270° step of Example
     6.2; the plate edge at n = ½)." · z*: "Not analytic: the two quotients differ, Cauchy–Riemann fails; Re and Im of
     z* are not a potential and a stream function of any flow."
- **Derivation tab:** **D14** (8 steps) `view: 'quot'` (phones keep `quot`), goal `set {kind: corner, n: 2, lh: -2}`;
  step 1 `set {lh: -1}` and `watch` "the solid curve (along x)"; step 2 `watch` "the dashed curve (along iy)"; step 4
  `set {lh: -5}` **live** "φ_x = **live px**, ψ_y = **live qy**"; step 5 `watch` "teal = purple mirrored". **D15** (7
  steps) `view: 'plane'`, goal `set {n: 2}`; step 3 `set {n: 3}` `watch` "walls at 0° and 60°"; step 6 `set {n: 0.6667}`
  **live** "\|dw/dz\| ∝ r^(−1/3)"; step 7 `set {n: 0.5}`. Interpret `s => "At your probe dw/dz = ${dw}: u = ${u}, v =
  ${v} m/s."`.
- **Code:**
  ```python
  p = ch06.complex_potential_probe("{{kind}}", {{x}}, {{y}}, A={{A}}, n={{n}}, h={{h}})
  print(p["qx_re"], p["qx_im"])   # along x:  {{qx}}
  print(p["qy_re"], p["qy_im"])   # along iy: {{qy}}
  print(p["u"], p["v"])           # dw/dz = u − iv (6.45): u = {{u}}, v = {{v}}
  print(p["cr1"], p["cr2"])       # Cauchy–Riemann (6.44): {{cr1}}, {{cr2}}
  print(ch06.corner_speed_exponent({{n}}))   # speed ∝ r^{{ex}} at the corner
  ```
- **Walkthrough (6 steps):** 1. "One function, two families" — "Re w draws equipotentials, Im w streamlines. Why must w be
  special?" · 2. "Two ways to differentiate" — "Along x and along iy: for z² both give 2z." `derive: {id: 'D14', step: 3}`
  · 3. "The velocity is mirrored" — "dw/dz = u − iv (6.45): teal is purple flipped." `readouts: ['u']`, `code: {lines: [4,
  4]}` · 4. "A non-flow" — "z*: the two derivatives differ. No flow." `set {kind: 'conj'}` · 5. "Corners" — "w = Azⁿ:
  walls at 0 and π/n (6.46). Speed ∝ r^{n−1}." `set {kind: corner, n: 0.667}`, `derive: {id: 'D15', step: 6}` · 6. "Your
  turn" — "Predict the speed exponent of a 60° corner, then set n." `controls: ['n']`.
- **Equations:** `w` ref 'Eq. (6.42)' `w\equiv\phi+i\psi` · `z` ref 'Eq. (6.43)' `z=x+iy=re^{i\theta}` · `cr` ref 'Eq. (6.44)'
  `\phi_x=\psi_y,\ \phi_y=-\psi_x` live · `cv` ref 'Eq. (6.45)' `dw/dz=u-iv` live · `corner` ref 'Eq. (6.46)'
  `w=Az^n,\ dw/dz=(A\pi/\alpha)z^{(\pi-\alpha)/\alpha}` · `vort` ref 'Eq. (6.47)' `w=-\frac{i\Gamma}{2\pi}\ln(z-z')` · `src` ref
  'Eq. (6.48)' · `dbl` ref 'Eq. (6.49)' `w=\frac{d}{2\pi(z-z')}`.
- **Check yourself:** (1) "At z = 1 + i with w = z², what are u and v?" — "dw/dz = 2 + 2i ⇒ u = 2, v = −2 m/s." `set
  {kind: corner, n: 2, x: 1, y: 1}` · (2) "Shrink h for z*. Do the two quotients ever agree?" — "No: 1 along x, −1 along
  iy at any h — z* has no derivative." · (3) "Which n makes the corner's tip a point of infinite speed?" — "n < 1
  (α > 180°): the exponent n − 1 is negative." · (4) "Choose the vortex. Is its φ single-valued round the centre?" —
  "No: φ = Γθ/2π jumps by Γ across the cut — the circulation (6.47)."
- **Selftest parity rows:** `{name: 'dwdz re z^2 at 1+i', js: probe('corner', 1, 1, {A: 1, n: 2, h: 1e-4}).dw.re, py:
  'ch06.complex_potential_probe("corner", 1.0, 1.0, A=1.0, n=2.0)["dwdz_re"]', rtol: 1e-12}` · `{name: 'v z^2 at 1+i',
  js: probe('corner', 1, 1, P2).v, py: 'ch06.complex_potential_probe("corner", 1.0, 1.0)["v"]', rtol: 1e-12}` · `{name:
  'CR for conj', js: probe('conj', 0.7, 0.2, P2).cr1, py: 'ch06.complex_potential_probe("conj", 0.7, 0.2)["cr1"]', rtol:
  1e-9}` · `{name: 'speed exponent n=2/3', js: speedExp(2/3), py: 'ch06.corner_speed_exponent(0.6666666666666666)',
  rtol: 1e-14}` · `{name: 'qx along x n=0.5', js: probe('corner', 0.3, 0.8, {A: 1, n: 0.5, h: 1e-4}).qx.re, py:
  'ch06.complex_potential_probe("corner", 0.3, 0.8, n=0.5, h=0.0001)["qx_re"]', rtol: 1e-10}`.
- **Fit plan:** 360×640: `plane` + `quot` (the CR bars as two rows under the curves); kind chips wrap; the loglog view
  hidden. Desktop: three views, loglog beside quot.

### E5 · blasius_kutta_contour
- **Title:** "Why doesn't lift depend on the shape?" · **Summary:** "Drag a contour round a body with circulation: Blasius's
  integral does not move; square the far field and only U × Γ/z has a residue — every body gets L = ρUΓ and no drag." ·
  **CORE:** C10 (also R18 (6.54), N53–N63 (6.55)–(6.61); the ellipse from C11) · **Reference:** `fid_formula_lab.html`
  (term-by-term bars and a total, clickable terms) with `amplitude_phase_second_order_II_3.html`'s numbered live
  explanation.
- **meta:** `viz:order 5` · `viz:sections 6.5` · `viz:equations 6.54 6.56 6.57 6.59 6.60 6.61 6.62` · `viz:fluidpy
  ch06.blasius_state ch06.laurent_contributions ch06.blasius_force ch06.laurent_coefficients` · `viz:derivations D16 D17
  D18`.
- **Physics:** `dwdzBody(body, z, p)` — cylinder (6.52); Zhukhovsky ellipse (chain rule through the outside-root inverse);
  tilted ellipse (the same with the stream at angle α); Rankine oval + vortex (four elements) ↔ `ch06.blasius_state`'s
  bodies; `blasius(contour, p, n)` periodic trapezoid (circle) or polygon rule ↔ `pf.blasius_force`; `laurent(p, R, K)`
  FFT (a small radix-2 JS FFT, n = 64) ↔ `pf.laurent_coefficients`; `contrib(p, R)` ↔ `ch06.laurent_contributions`;
  `crosses(contour, body)` (min distance of the contour to the body's surface < 0 ⇒ ⚠️) ↔ `blasius_state(...)["crosses_body"]`;
  the pressure route on the body ↔ `ch06.contour_force`.
- **Modes (body chips):** "cylinder" · "ellipse (Zhukhovsky)" · "tilted ellipse" · "Rankine oval + vortex".
- **Views** (rows [1.3, 1]):
  1. `body` "The body and your contour" (equal aspect) — streamlines, the body (black), the contour (amber, dashed,
     draggable: drag its rim to change R, drag its centre to offset it), the integrand direction arrows along it; the
     title shows "D = … · L = … N/m (Blasius on your contour)".
  2. `bars` "Which term survives?" — the contribution of each power z^k (k = 0, −1, −2, −3, −4) of (dw/dz)² to
     (iρ/2)∮(dw/dz)² dz, split into drag (muted) and lift (amber) parts; only k = −1 non-zero; clicking a bar lights its
     term in the Equations tab.
  3. `radius` "Force vs contour size" (`hidePortrait: true`) — D and L against R/a (log) with a ◇ for the on-body
     pressure route (6.56); flat lines; a rose segment where the contour cuts the body.
- **Controls:** `body` chips · `G` "Circulation $\Gamma$ (clockwise)" −4…4 m²/s, default 2 · `U` "Stream $U$" 1…20 m/s,
  default 10 (optional) · `R` "Contour radius $R$" 0.05…1 m (log), default 0.2 (also by dragging) · `n` "Quadrature
  points" 8…512 (log2 steps), default 256 (optional).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **terms** (Laurent contributions) · **modes**
  (bodies) · **presets** "cylinder, Γ = 0 (d'Alembert)" {body: cylinder, G: 0}, "cylinder with Γ" {G: 2}, "ellipse, same
  Γ" {body: ellipse, G: 2}, "contour cuts the body ⚠️" {R: 0.07}, "few quadrature points" {n: 8} · **status** ("✅ contour
  encloses the body and no other singularity" / "⚠️ contour crosses the body — Cauchy does not apply") · **inspector**
  (click a contour point: "(dw/dz)² dz there = … ; its share of the sum").
- **Readouts:** "Lift $L$" (N/m) · "Drag $D$" (N/m) · "$\rho U\Gamma$" (N/m) · "Residue" (m³/s²).
- **Explain:**
  0. *What the views show* — "Your amber contour round the body; the bars split Blasius's integral by powers of z; the
     third view (hidden on phones) plots the force against the contour's size."
  1. *Far away every body looks alike* — "Laurent coefficients of dw/dz on your contour: c₀ = **live c0** (= U), c₋₁ =
     **live cm1** (= iΓ/2π with Γ clockwise), c₋₂ = **live cm2** (the body's doublet)" (N61).
  2. *Square it* — "(dw/dz)² = U² + (iUΓ/π)/z + [−(Ud/π + Γ²/4π²)]/z² + …; the 1/z coefficient is iUΓ/π = **live res**.
     ⚠️ The book prints the 1/z² coefficient as (Ud/π − Γ²/4π²) inside an extra square — it does not matter: only 1/z
     counts."
  3. *The residue theorem* — "∮z^{−k}dz = 2πi only for k = 1, so ∮(dw/dz)² dz = 2πi × iUΓ/π = −2UΓ = **live I**; (6.60)
     D − iL = (iρ/2)(−2UΓ) = −iρUΓ ⇒ D = 0, L = ρUΓ = **live LKJ** N/m" (boxed).
  4. *Your contour, numerically* — "trapezoid with n = **live n** points: D = **live D**, L = **live L** N/m (error
     **live err**)."
  5. *The body route* — "(6.56) D = −∮p dy, L = ∮p dx on the body itself: L = **live Lp** N/m — the same."
  6. *Force vs size* — hint on phones; else "flat lines: Cauchy lets the contour move anywhere the integrand is
     analytic."
  7. *Reading the current setting* — "✅ The contour encloses the body and nothing else: its size and shape do not
     matter; any body with circulation Γ gets L = **live LKJ** N/m and D = 0 — a more general d'Alembert." · ⚠️ "Your
     contour cuts the body: between it and the body lies no fluid, the integrand is not analytic there, and the number
     is meaningless." · few points: "A smooth analytic integrand on a circle: the trapezoid error falls exponentially
     with n (P136)."
- **Derivation tab:** **D16** (8 steps) `view: 'body'`, goal `set {body: cylinder, G: 2, R: 0.1}` (contour on the body);
  step 5 `watch` "the outward normal arrows"; step 7 **live** "D = **live Dp**, L = **live Lp** N/m". **D17** (12 steps,
  ★★★) `view: 'body'`; step 8 `watch` "on the body, u + iv and dz point the same way"; step 11 `set {R: 0.3}` `watch`
  "the readout does not move"; step 12 `set {R: 0.07}` `watch` "cutting the body breaks it". **D18** (11 steps, ★★★)
  `view: 'bars'` (phones keep `bars`); step 6 highlights `term:k-1`, **live** "iUΓ/π = **live res**"; step 10 **live** "D −
  iL = −i × **live LKJ**". Interpret `s => "Your body (${body}) with Γ = ${G} m²/s: L = ${L} N/m on a contour of radius ${R}
  m."`.
- **Code:**
  ```python
  s = ch06.blasius_state("{{body}}", Gamma_cw={{G}}, U={{U}}, R={{R}}, rho=1.2, n={{n}})
  print(s["D"], s["L"])            # Blasius on your contour: {{D}}, {{L}} N/m
  print(s["L_KJ"])                 # ρUΓ = {{LKJ}} N/m (6.62)
  print(s["cm1_im"])               # c₋₁ = iΓ/2π → {{cm1}}
  print(s["crosses_body"])         # {{crosses}}
  c = ch06.laurent_contributions("{{body}}", R={{R}}, Gamma_cw={{G}}, U={{U}})
  print(c["contrib_im"])           # only k = −1 non-zero: {{bars}}
  ```
- **Walkthrough (6 steps):** 1. "Two different wings" — "A tube and a flat ellipse, the same circulation. Same lift?"
  `set {body: cylinder, G: 2}` · 2. "Pressure as one integral" — "Complex numbers turn the pressure integral into
  (iρ/2)∮(dw/dz)² dz (6.60)." `derive: {id: 'D17', step: 10}` · 3. "Move the contour" — "Drag the amber circle: L stays
  put — Cauchy's theorem." `controls: ['R']`, `readouts: ['L']` · 4. "Only one term" — "Square the far field: only U ×
  iΓ/2πz makes a 1/z term." `terms: true`, `derive: {id: 'D18', step: 6}` · 5. "Any body" — "Switch to the ellipse: same
  L = ρUΓ (6.62)." `set {body: ellipse}`, `code: {lines: [3, 3]}` · 6. "Your turn" — "Predict what happens if the contour
  cuts the body, then try." `set {R: 0.07}`.
- **Equations:** `cv` "Force on a body" ref 'Eq. (6.56)' `D=-\oint p\,dy,\ L=\oint p\,dx` · `cf` ref 'Eq. (6.57)'
  `D-iL=-i\oint p\,dz^*` · `tan` ref 'Eq. (6.59)' `(u+iv)dz^*=(u-iv)dz` · `bl` ref 'Eq. (6.60)'
  `D-iL=\frac{i\rho}{2}\oint(dw/dz)^2dz` live · `far` ref 'Eq. (6.61)' (corrected coefficient, the printed one noted) ·
  `kj` ref 'Eq. (6.62)' `D=0,\ L=\rho U\Gamma` live.
- **Check yourself:** (1) "Set Γ = 0 on the ellipse. What are D and L?" — "Both zero: with no circulation there is no
  1/z term at all — d'Alembert for any shape." `set {body: ellipse, G: 0}` · (2) "Double U. What happens to the
  residue?" — "It doubles (iUΓ/π), and so does L." · (3) "Why does the contour cutting the body give nonsense?" — "Inside
  the body the formula continues into a region that is not fluid: the integrand has singularities (the doublet, the
  vortex) that the contour now passes, so Cauchy's theorem no longer applies." · (4) "With n = 8 points is the lift
  already right?" — "Nearly (error ~1e-4 relative): for an analytic integrand the trapezoid converges exponentially."
- **Selftest parity rows:** `{name: 'L cylinder R=0.2', js: blasius(circ(0.2), P0, 256).L, py:
  'ch06.blasius_state("cylinder", Gamma_cw=2.0, U=10.0, R=0.2)["L"]', rtol: 1e-10}` · `{name: 'D cylinder R=1', js:
  blasius(circ(1), P0, 256).D, py: 'ch06.blasius_state("cylinder", Gamma_cw=2.0, U=10.0, R=1.0)["D"]', rtol: 0, atol:
  1e-9}` · `{name: 'L ellipse', js: blasius(circ(0.5), PE, 256).L, py: 'ch06.blasius_state("ellipse", Gamma_cw=2.0,
  U=10.0, a=0.12, b=0.1, R=0.5)["L"]', rtol: 1e-9}` · `{name: 'c-1 imag', js: laurent(P0, 0.3, 4)[-1].im, py:
  'ch06.blasius_state("cylinder", Gamma_cw=2.0, U=10.0, R=0.3)["cm1_im"]', rtol: 1e-10}` · `{name: 'k=-1 bar', js:
  contrib(P0, 0.2).lift[1], py: 'ch06.laurent_contributions("cylinder", R=0.2, Gamma_cw=2.0, U=10.0)["contrib_im"][1]',
  rtol: 1e-9}` · exact `{name: 'crossing flag', js: crosses(circ(0.07), 'cylinder'), expect: true}`.
- **Fit plan:** 360×640: `body` + `bars` (bars horizontal, 5 rows of 22 px); status one line; body chips wrap. Desktop:
  third view beside the bars.

### E6 · conformal_joukowski
- **Title:** "How does a map carry a flow?" · **Summary:** "Two linked planes: a circle's flow in the ζ-plane and its image
  under z = ζ + b²/ζ; small crosses keep their right angles, the circle becomes an ellipse — and only the outside square
  root keeps the flow outside the body." · **CORE:** C11 (also N64–N70 (6.63)–(6.69), C09's w, N66 grid maps) ·
  **Reference:** `angular_frequency_explorer_1.html` (linked views on one state, modes, a highlighted table of special
  cases).
- **meta:** `viz:order 6` · `viz:sections 6.6` · `viz:equations 6.63 6.64 6.65 6.66 6.67 6.68 6.69` · `viz:fluidpy
  ch06.joukowski_state ch06.joukowski_ellipse ch06.joukowski_inverse ch06.angle_preservation ch06.elliptic_cylinder_flow`
  · `viz:derivations D19 D20 D21`.
- **Physics:** `jk(zeta, b)` ↔ `cm.joukowski`; `jkInv(z, b, branch)` (outside: ½[z + √(z − 2b)√(z + 2b)]; principal:
  ½[z + √(z² − 4b²)] with the principal complex root) ↔ `cm.joukowski_inverse`; `ellipse(a, b)` ↔ `ch06.joukowski_ellipse`;
  `Wzeta(zeta, p)` = U(ζe^{−iα} + a²e^{iα}/ζ) + (iΓ/2π)ln(ζ/a) and its derivative; velocity by the chain rule ↔
  `ch06.elliptic_cylinder_flow(...).dwdz`; `angles(f, fp, z0, d1, d2)` ↔ `cm.angle_preservation`; the probe state ↔
  `ch06.joukowski_state`.
- **Modes (map chips):** "Zhukhovsky z = ζ + b²/ζ" · "z²" · "e^z" (the last two for angle preservation and grid images;
  flows drawn as images of a uniform grid).
- **Views** (rows [1.25, 1] on wide; two columns):
  1. `zeta` "ζ-plane (the circle)" (equal aspect) — circle \|ζ\| = a (black), critical points ±b (×), streamlines of
     (6.68), a small cross (two short orange and teal arms) at the draggable element point, the probe.
  2. `z` "z-plane (the body)" — the image body (ellipse or slit), the dashed slit [−2b, 2b], the mapped streamlines (same
     ψ levels), the image cross with its angle readout; with the principal branch, the streamlines computed through the
     wrong root drawn rose where \|ζ\| < a.
  3. `jac` "Stretch round the circle" (`hidePortrait: true`) — \|dz/dζ\| against the angle on the circle (zero at ±b when
     a = b).
- **Controls:** `ab` "Circle $a/b$" 1…3, step 0.01, default 1.2 · `G` "Circulation $\Gamma$ (clockwise)" 0…4 m²/s,
  default 0 · `al` "Stream angle $\alpha$" −15°…15°, default 0 (optional) · `branch` chips "outside root · numpy
  principal" · `map` chips (modes).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (three maps) · **presets** "plate
  a = b" {ab: 1}, "ellipse a = 1.2b" {ab: 1.2}, "almost a circle a = 3b" {ab: 3}, "no circulation" {G: 0}, "wrong
  branch" {branch: 'principal'}, "z² at 0" {map: 'z2', elementAt: 0} · **inspector** (the element's mapping arithmetic:
  "f′(ζ₀) = 1 − b²/ζ₀² = 0.31 − 0.12i: stretch ×0.33, turn −21°") · **status** ("✅ outside root: \|ζ\| = 2.70 > a" / "⚠️
  principal root lands inside the circle (\|ζ\| = 0.37)").
- **Readouts:** "Semi-axes $A$, $B$" (m) · "Foci" (m) · "\|ζ\| at probe" (m) · "Angle kept?" (°).
- **Explain:**
  0. *What the views show* — "Left: the circle's flow in the ζ-plane; right: its image under z = ζ + b²/ζ. The <b
     class="c-orange">orange</b> and <b class="c-teal">teal</b> arms of the small cross and their images show the angle
     between two small steps."
  1. *The local stretch and turn (6.63)* — "f′(ζ₀) = 1 − b²/ζ₀² = **live fp**: every small step is stretched by \|f′\| =
     **live mag** and turned by arg f′ = **live arg**°."
  2. *Angles kept (6.64)* — "before α = **live a1**°, after β = **live a2**° — equal because both arms turn by the same
     arg f′ (fails where f′ = 0: ζ = ±b)."
  3. *The body (6.66)–(6.67)* — "a = **live a**, b = 1: semi-axes A = a + b²/a = **live A**, B = a − b²/a = **live B** m;
     foci ±√(A² − B²) = ±2b = **±2** m" (boxed).
  4. *Back from z to ζ (6.69)* — "at the probe z = **live z**: roots ζ₁ = **live r1** (\|ζ\| = **live m1**), ζ₂ = b²/ζ₁ =
     **live r2** (\|ζ\| = **live m2**); their product is b²; the outside root is **live chosen**" (boxed).
  5. *Velocity by the chain rule* — "u − iv = (dW/dζ)(dζ/dz) = **live dW** × **live dz** = **live uv** ⇒ u = **live
     u**, v = **live v** m/s."
  6. *Stretch round the circle* — hint on phones; else "\|dz/dζ\| is smallest near ζ = ±b: the ellipse's sharp ends."
  7. *Reading the current setting* — a = b: "The circle passes through the critical points: the image is a flat plate
     and the speed at its edges is infinite unless Γ cancels it (the Kutta condition, Ch. 14)." · a slightly > b:
     "A thin ellipse with sharp but rounded ends." · a ≫ b: "The map is nearly the identity: almost a circle." · wrong
     branch: "⚠️ numpy's principal root picks the root inside the circle for points left of the slit: those points are
     sent into the cylinder, and the 'flow' there is not the flow round the ellipse."
- **Derivation tab:** **D19** (6 steps) `view: 'z'`, goal `set {map: 'jk'}`; step 2 **live** "\|f′\| = **live mag**, arg
  f′ = **live arg**°"; step 5 `set {map: 'z2', elementAt: 0}` `watch` "the image angle doubles". **D20** (8 steps) `view:
  'z'`; step 2 `set {ab: 1}` `watch` "the plate"; step 5 `set {ab: 1.2}` **live** "A = **live A**, B = **live B**"; step 7
  `watch` "the × marks are the plate's ends". **D21** (11 steps, ★★★) `view: 'z'`; step 5 `set {branch: 'principal'}`
  `watch` "rose streamlines on the left"; step 6 `set {branch: 'outside'}`; step 11 **live** "u − iv = **live uv**".
  Interpret `s => "Your ellipse: ${A} × ${B} m; at the probe the outside root has |ζ| = ${m}."`.
- **Code:**
  ```python
  st = ch06.joukowski_state(a={{a}}, b=1.0, Gamma_cw={{G}}, x={{x}}, y={{y}},
                            branch="{{branch}}")
  print(st["A"], st["B"], st["foci"])    # {{A}} {{B}} ±{{foci}} m  (6.67)
  print(st["zeta_abs"], st["inside"])    # |ζ| = {{m}}; inside the circle? {{inside}}
  print(st["u"], st["v"])                # velocity at the probe {{u}}, {{v}} m/s
  z = complex({{x}}, {{y}})               # ζ = ½[z + √(z − 2b)√(z + 2b)]  (6.69)
  ```
- **Walkthrough (6 steps):** 1. "Borrow a solution" — "We know the flow round a circle. Can a map give other bodies?" ·
  2. "Small crosses stay crosses" — "Every small step turns by the same arg f′: right angles survive (6.64)." `derive:
  {id: 'D19', step: 4}` · 3. "Circle → ellipse" — "z = ζ + b²/ζ sends \|ζ\| = a to an ellipse with foci ±2b (6.67)."
  `readouts: ['A']`, `derive: {id: 'D20', step: 5}` · 4. "Carry the flow" — "ψ at ζ = ψ at its image: streamlines map to
  streamlines." `set {G: 1}` · 5. "The right root" — "Two roots, product b²: take the outside one. numpy's may not."
  `set {branch: 'principal'}`, `code: {lines: [6, 6]}` · 6. "Your turn" — "Predict the ellipse for a = 2b, then slide."
  `controls: ['ab']`.
- **Equations:** `dm` ref 'Eq. (6.63)' `\delta w=\frac{dw}{dz}\delta z` · `ang` ref 'Eq. (6.64)' `\alpha=\beta` · `jk` ref 'Eq.
  (6.65)' `z=\zeta+b^2/\zeta` · `circ` ref 'Eq. (6.66)' `z=ae^{i\theta}+\frac{b^2}{a}e^{-i\theta}` · `ell` ref 'Eq. (6.67)' live ·
  `wz` ref 'Eq. (6.68)' `w=U(\zeta+a^2/\zeta)+\frac{i\Gamma}{2\pi}\ln(\zeta/a)` · `inv` ref 'Eq. (6.69)'
  `\zeta=\tfrac12z+\tfrac12(z^2-4b^2)^{1/2}` (with "outside root" in the note).
- **Check yourself:** (1) "Where on the circle do angles fail to be kept?" — "At ζ = ±b, where dz/dζ = 1 − b²/ζ² = 0 — the
  plate's ends when a = b." · (2) "For a = 1.2 b, what are the semi-axes?" — "2.033 b and 0.367 b; foci ±2b." `set {ab:
  1.2}` · (3) "Switch to numpy's principal root. Which points go wrong?" — "Those left of the slit (Re z < 0): the
  principal √(z² − 4b²) returns the root that lands inside the circle." `set {branch: 'principal'}` · (4) "Set a = b and
  Γ = 0. Where is the flow fastest?" — "At the plate's edges: dζ/dz → ∞ there."
- **Selftest parity rows:** `{name: 'semi-axis A', js: ellipse(1.2, 1).A, py: 'ch06.joukowski_ellipse(1.2, 1.0)["A"]',
  rtol: 1e-14}` · `{name: 'outside |zeta| at -3+0.5i', js: jkState(P0, -3, 0.5, 'outside').m, py:
  'ch06.joukowski_state(a=1.2, b=1.0, x=-3.0, y=0.5, branch="outside")["zeta_abs"]', rtol: 1e-12}` · `{name: 'principal
  |zeta| at -3+0.5i', js: jkState(P0, -3, 0.5, 'principal').m, py: 'ch06.joukowski_state(a=1.2, b=1.0, x=-3.0, y=0.5,
  branch="principal")["zeta_abs"]', rtol: 1e-12}` · `{name: 'u at probe', js: jkState({...P0, G: 1}, 0.5, 1.5,
  'outside').u, py: 'ch06.joukowski_state(a=1.2, b=1.0, Gamma_cw=1.0, x=0.5, y=1.5)["u"]', rtol: 1e-10}` · `{name:
  'v at probe', js: jkState({...P0, G: 1}, 0.5, 1.5, 'outside').v, py: 'ch06.joukowski_state(a=1.2, b=1.0,
  Gamma_cw=1.0, x=0.5, y=1.5)["v"]', rtol: 1e-10}`.
- **Fit plan:** 360×640: `zeta` above `z` (both square-ish, 45 % each); the jac view hidden; status one line; map chips
  optional on phones (Zhukhovsky default). Landscape/desktop: `zeta` and `z` side by side, `jac` below.

### E7 · laplace_relaxation
- **Title:** "How does a grid solve Laplace's equation?" · **Summary:** "Every value becomes the average of its four
  neighbours, again and again: step one node at a time or whole sweeps, compare Jacobi, Gauss–Seidel and SOR, and stop on
  the residual." · **CORE:** C12 (also R19–R21 (6.70)–(6.71), N72 (6.73), N73, N74; N71 named) · **Reference:**
  `pixels_as_parameters.html` (an algorithm in slow motion, step back/forward, a clickable loss curve) with
  `forced_damped_vibrations.html`'s explanation panel.
- **meta:** `viz:order 7` · `viz:sections 6.7` · `viz:equations 6.70 6.71 6.72 6.73` · `viz:fluidpy ch06.relaxation_state
  ch06.four_point_system ch06.node_update ch06.gauss_seidel_sweep ch06.jacobi_sweep ch06.sor_sweep` · `viz:derivations D22
  D23`.
- **Physics:** `problem(name, refine)` builds the mask and boundary values — "four-point" (ψ = xy boundary, 4 × 4),
  "contraction" (the Example 6.2 grid, Q = 1, the same geometry array as `ch06.example_6_2`; the builder exports it once
  from Python into the HTML as a small integer mask), "ψ = xy square" (8 × 8) ↔ `ch06.relaxation_state`; `nodeUpdate`,
  `sweep(method, ω)` ↔ `ls.node_update`, `ls.gauss_seidel_sweep`, `ls.jacobi_sweep`, `ls.sor_sweep`; `residual(psi)` ↔
  `ls.residual_norm`; `rhoJacobi(N)` = cos(π/N) and ω_opt = 2/(1 + sin(π/N)) for the predicted rates.
- **Modes (method chips):** "Jacobi" · "Gauss–Seidel" · "SOR".
- **Views** (rows [1.3, 1]):
  1. `grid` "The grid" — nodes as squares coloured by ψ/Q (numbers printed on grids ≤ 8 × 8), boundary nodes outlined,
     the node being updated highlighted with its four neighbours joined by lines, ψ contours drawn once converged. Pointer:
     click a node → inspector.
  2. `hist` "Residual and change per sweep" — log₁₀ residual (solid) and max change (dashed) against sweep for the
     current method, faint ghosts of the other two methods (precomputed), a moving dot at the current sweep; click the
     curve to jump to a sweep.
  3. `matrix` "The four equations" (`hidePortrait: true`; four-point problem only) — the 4 × 4 system Aψ = b with the
     current iterate beside the exact solution; otherwise the spectral-radius numbers.
- **Controls:** `method` chips · `omega` "SOR factor $\omega$" 1…1.95, step 0.01, default ω_opt (SOR only) · `prob`
  chips "4-point · contraction · ψ = xy" · `refine` chips "×1 · ×2 · ×4" (optional; contraction) · `init` chips "zero ·
  linear" (optional).
- **Transport:** `k` (node updates) 0 → N_nodes × 60, rate = one sweep per 0.5 s, `end: 'hold'`; ◀ ▶ step one node;
  "sweep" buttons as a custom control (two buttons in the transport strip, hidden on phones); end card: "Converged in 19
  sweeps (Gauss–Seidel) — residual 7e-11; Jacobi needed 34."
- **Depth features:** Explain + Code + Derivation · **transport** (node / sweep stepping) · **linked views** (3) ·
  **modes** (methods) · **presets** "4-point system" {prob: '4-point', method: 'Gauss–Seidel'}, "book grid" {prob:
  'contraction', refine: 1}, "refined ×4" {refine: 4}, "SOR at ω_opt" {method: 'SOR'}, "harmonic test ψ = xy" {prob: 'xy'}
  · **inspector** (click a node: "¼(ψ_W + ψ_E + ψ_S + ψ_N) = ¼(0 + 2 + 0 + 2) = 1.000") · **status** ("⏳ iterating ·
  sweep 7 · residual 3.1e-5" / "✅ converged: residual < 10⁻¹⁰ after 19 sweeps").
- **Readouts:** "Sweep" · "Residual" (m²/s) · "Max change" (m²/s) · "Predicted ρ" (–).
- **Explain:**
  0. *What the views show* — "Coloured squares: ψ/Q at each node (outlined = fixed boundary values). The highlighted
     node is being replaced by the average of its four linked neighbours. The second view: how far the grid is from
     satisfying (6.72) after each sweep (solid) and how much the last sweep changed it (dashed)."
  1. *The rule (6.72)* — "ψ_{i,j} = ¼(ψ_{i−1,j} + ψ_{i+1,j} + ψ_{i,j−1} + ψ_{i,j+1}) = ¼(**live nW** + **live nE** + **live
     nS** + **live nN**) = **live avg**" (boxed) — "from (6.70) + (6.71) with Δx = Δy (D22)".
  2. *The change* — "old **live old** → new **live avg**: change **live dlt**; Gauss–Seidel already used the new
     values of the nodes before this one in the sweep (Jacobi would not)."
  3. *The residual* — "max over the grid of \|average − value\| = **live res**; this, not the change, says how far from
     the solution we are."
  4. *How fast* — "error × ρ per sweep: Jacobi ρ = **live rJ**, Gauss–Seidel ρ² = **live rG**, SOR ω_opt − 1 = **live
     rS**; sweeps to 10⁻¹⁰ ≈ ln(10⁻¹⁰)/ln ρ = **live pred**."
  5. *The matrix view* — hint on phones; else "four equations, four unknowns: A = [[4, −1, −1, 0], …], b from the
     boundary; the iterate beside the exact [1, 2, 2, 4]."
  6. *Right now* — "sweep **live k** of this method; residual **live res**."
  7. *Reading the current setting* — iterating: "The boundary information is still spreading inward — each sweep moves
     it about one node." · converged: "Every value is the average of its neighbours to 10⁻¹⁰: the grid solves Laplace's
     equation." · Jacobi vs GS: "Gauss–Seidel needs about half Jacobi's sweeps (ρ_GS = ρ_J²); SOR at ω_opt needs far fewer
     — about N instead of N²." · ψ = xy: "A quadratic harmonic polynomial satisfies the five-point rule exactly (its
     fourth derivatives vanish, D22): the converged grid is exact."
- **Derivation tab:** **D22** (6 steps) `view: 'grid'`, goal `set {prob: '4-point'}`; step 2 highlights the node (2,2)
  and its neighbours; step 5 **live** "¼(**live nW** + … ) = **live avg**". **D23** (9 steps) `view: 'hist'` (phones keep
  `hist`); step 4 `set {method: 'Jacobi'}` `play`; step 5 `set {method: 'Gauss–Seidel'}` **live** "first sweep: 0, 0.75,
  0.75, 3.375"; step 6 **live** "ρ_J = ½, ρ_GS = ¼"; step 8 `watch` "the dashed change below the solid residual".
  Interpret `s => "${method}: residual ${res} after ${k} sweeps; predicted ${pred} sweeps to 1e-10."`.
- **Code:**
  ```python
  s = ch06.relaxation_state("{{prob}}", method="{{method}}", sweeps={{k}},
                            omega={{omega}})
  print(s["residual"], s["max_change"])   # {{res}}, {{dlt}}
  print(s["psi22"], s["psi32"], s["psi23"], s["psi33"])   # 4-point: {{psi4}}
  print(s["sweeps_to_tol"])               # {{pred}} sweeps to 1e-10
  # one node: psi[j, i] = 0.25*(psi[j, i-1] + psi[j, i+1] + psi[j-1, i] + psi[j+1, i])
  ```
- **Walkthrough (7 steps):** 1. "A membrane on a frame" — "Every point of a stretched membrane sits at its neighbours'
  average. So does a harmonic ψ." `set {prob: '4-point', k: 0}` · 2. "The rule" — "Two second differences, Δx = Δy:
  ψ_{i,j} = ¼ × (four neighbours) (6.72)." `derive: {id: 'D22', step: 5}` · 3. "Four equations" — "The 16-point grid has
  four unknowns: (6.73) is Aψ = b." `derive: {id: 'D23', step: 2}` · 4. "Step by step" — "Press ▶ one node: Gauss–Seidel
  uses the newest values." `controls: ['method']`, `code: {lines: [5, 5]}` · 5. "Sweeps and the residual" — "The residual
  falls by a constant factor per sweep: ¼ here." `play: true` · 6. "A real problem" — "The contraction: SOR at ω_opt beats
  Gauss–Seidel by far." `set {prob: 'contraction', method: 'SOR'}` · 7. "Your turn" — "Predict: refine ×2. Sweeps ×4 for
  Gauss–Seidel? Try it." `controls: ['refine']`.
- **Equations:** `d2x` ref 'Eq. (6.70)' `(\psi_{i+1,j}-2\psi_{i,j}+\psi_{i-1,j})/\Delta x^2` · `d2y` ref 'Eq. (6.71)' · `avg`
  ref 'Eq. (6.72)' live · `four` ref 'Eq. (6.73)' (the four equations, `aligned`) · `gs` "Gauss–Seidel update" (no ref).
- **Check yourself:** (1) "After one Gauss–Seidel sweep from zero on the 4-point problem, what is ψ₃₃?" — "3.375: it
  already uses the new ψ₃₂ = 0.75 and ψ₂₃ = 0.75 (¼(0.75 + 6 + 0.75 + 6))." `set {prob: '4-point', method:
  'Gauss–Seidel', k: 4}` · (2) "Why can the max change be tiny while the answer is still wrong?" — "When ρ ≈ 1 each sweep
  changes little although the error is large: change ≈ (1 − ρ) × error." · (3) "Which method needs fewest sweeps on the
  refined contraction, and roughly how many fewer?" — "SOR at ω_opt: ≈ N sweeps instead of ≈ N² for Gauss–Seidel." · (4)
  "Is the converged ψ = xy grid exact or approximate?" — "Exact: xy is a harmonic polynomial the five-point rule
  reproduces without error."
- **Selftest parity rows:** `{name: 'GS 1 sweep psi32', js: relax('4-point', 'GS', 1).psi[1], py:
  'ch06.relaxation_state("four_point", method="gauss_seidel", sweeps=1)["psi32"]', rtol: 1e-14}` · `{name: 'GS 1 sweep
  psi33', js: relax('4-point', 'GS', 1).psi[3], py: 'ch06.relaxation_state("four_point", method="gauss_seidel",
  sweeps=1)["psi33"]', rtol: 1e-14}` · `{name: 'Jacobi 1 sweep psi33', js: relax('4-point', 'J', 1).psi[3], py:
  'ch06.relaxation_state("four_point", method="jacobi", sweeps=1)["psi33"]', rtol: 1e-14}` · `{name: 'GS sweeps to tol',
  js: relax('4-point', 'GS', 200).sweepsToTol, py: 'ch06.relaxation_state("four_point",
  method="gauss_seidel", sweeps=200)["sweeps_to_tol"]', rtol: 0}` · `{name: 'contraction residual after 10 GS', js:
  relax('contraction', 'GS', 10).residual, py: 'ch06.relaxation_state("contraction", method="gauss_seidel",
  sweeps=10)["residual"]', rtol: 1e-10}`.
- **Fit plan:** 360×640: `grid` + `hist`; numbers on nodes only when a node is ≥ 26 px; the matrix view hidden (its key
  number, the exact solution, in the grid title); transport step/sweep buttons hidden on phones (play only). Desktop:
  three views.

### E8 · axial_singularity_bodies
- **Title:** "Given a shape, which sources draw it?" · **Summary:** "Hide unknown sources on the axis, demand ψ = 0 at N
  points on the body, and one linear solve draws the hull; raise N and watch the fit improve while the system grows
  touchy — or put the sources on the surface (panels)." · **CORE:** C14 (also N86–N89 (6.93)–(6.95), N90, C13's sphere as a
  target, N84) · **Reference:** `stride_padding_playground.html` (presets that are the classic settings, the formula with
  numbers plugged in, badges explaining the result) with `angular_frequency_explorer_1.html`'s modes.
- **meta:** `viz:order 8` · `viz:sections 6.7 6.8` · `viz:equations 6.87 6.93 6.94 6.95` · `viz:fluidpy ch06.axial_state
  ch06.airship ch06.line_sink_stream_function ch06.axisym_body_target ch06.panel_cp_error` · `viz:derivations D26 D27`.
- **Physics:** `lineSinkPsi(R, z, k, a)` (6.94) ↔ `ch06.line_sink_stream_function`; `airship(U, Q, a)` (axis stagnation
  points by `Viz.num.brentq`, the ψ = 0 contour) ↔ `ch06.airship`; `target(kind, N, fineness)` ↔
  `ch06.axisym_body_target`; `axialSolve(zb, Rb, U, N)` (dense Gaussian elimination with partial pivoting, local helper
  `lu`) and `cond1(A)` (1-norm condition number from the explicit inverse, N ≤ 80) ↔ `pn.axial_singularity_solve` /
  `ch06.axial_state` (**Part C: `cond` is the 1-norm condition number, `np.linalg.cond(A, 1)`, on both sides**); panel
  mode `panels(xb, yb, U)` (constant-strength source panels, analytic influence) ↔ `pn.source_panels`, error ↔
  `ch06.panel_cp_error`.
- **Modes:** "axial sources (book §6.8)" · "source panels (2-D, our extension)".
- **Views** (rows [1.2, 1]):
  1. `body` "The body" — axial mode: meridian half-plane (z horizontal, R up, mirrored faintly below), the target body
     (dashed grey), the computed ψ = 0 surface (black), the N axial segments coloured by k_n (teal > 0 sources, rose < 0
     sinks, thickness ∝ \|k_n\|), the N collocation points (dots), a few streamlines; panel mode: the 2-D body with its
     panels, midpoints, and C_p arrows. Pointer: click a collocation point → inspector.
  2. `bars` "Strengths along the axis" — k_n bars against z with the exact distribution (Rankine oval: two spikes drawn
     as delta arrows with their areas; airship: point source + flat line sink −k) as a grey reference; panel mode: C_p at
     midpoints vs 1 − 4 sin²θ.
  3. `conv` "Fit vs N" (`hidePortrait: true`) — body error (black) and condition number (amber) against N (log–log),
     the current N marked.
- **Controls:** `mode` chips · `tgt` chips "airship (6.95) · Rankine oval · sphere · ellipsoid" (axial) · `N` "Segments
  $N$" 4…80 (axial) / "Panels $N$" 8…128 (panels), integer, default 20 · `fin` "Fineness L/D" 1.5…6, default 3 (ellipsoid,
  Rankine; optional) · `U` "Stream $U$" 0.5…3 m/s, default 1 (optional).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** · **presets** "airship (6.95)" {tgt:
  airship}, "Rankine oval (exact recovery)" {tgt: rankine, N: 40}, "sphere (hard)" {tgt: sphere, N: 40}, "panels N = 8"
  {mode: panels, N: 8}, "panels N = 64" {mode: panels, N: 64} · **inspector** (ψ_m as a sum of segment terms: "ψ_m =
  −Σ(k_n/4π)(r^m_{n−1} − r^m_n) + ½UR_m² = −(0.52 + 0.31 − …) + 0.72 = 3e-15") · **status** ("✅ closed: Σk_nΔξ =
  −2e-14 · cond 4.1e3" / "⚠️ cond 3.8e11: the digits are noise — a blunt body" / panels: "C_p error 1.7e-3 (∝ 1/N²)").
- **Readouts:** "Body error" (m) · "Σ k_n Δξ" (m³/s) · "cond" (–) · "Length" (m, airship).
- **Explain:**
  0. *What the views show* — "Dashed grey: the body we want; black: the stream surface ψ = 0 our sources produce; bars:
     the strength of each axial segment (<b class="c-teal">teal</b> sources, <b class="c-rose">rose</b> sinks)."
  1. *One segment's contribution (6.94)* — "a source segment from ξ_{n−1} to ξ_n gives ψ_mn = −(k_n/4π)(r^m_{n−1} − r^m_n) at
     body point m (a sink (6.94) with the sign flipped): at your clicked point, segment 3 gives −(**live k3**/4π)(**live
     r2** − **live r3**) = **live psi3** m³/s."
  2. *The system* — "ψ_m = −Σ_n A_mn k_n + ½UR_m² = 0 for m = 1…N ⇒ A k = ½UR² with A_mn = (r^m_{n−1} − r^m_n)/4π: an
     N × N linear system (**live N** unknowns)."
  3. *Closure* — "Σ k_n Δξ = **live net** m³/s: sources and sinks cancel, so the body closes (N86)."
  4. *How good* — "max distance between the computed and target surfaces = **live err** m."
  5. *How touchy* — "cond₁(A) = **live cond**: relative errors in the data can grow by this factor."
  6. *The airship (6.95)* — "k = Q/a = **live k** m²/s; nose at z = **live zf**, tail at **live zr** m, length **live L**
     m (u_z = 0 on the axis, `brentq`)."
  7. *Fit vs N* — hint on phones; else "error falls with N while cond climbs."
  8. *Reading the current setting* — Rankine: "Its exact sources are a point source and a point sink on the axis; the
     segments converge to two narrow spikes." · airship: "The fitted bars reproduce a point source and a flat line sink:
     the inverse method recovers (6.95)." · sphere: "A sphere's exact singularity is a point doublet; a line of sources can
     only approximate it, and the system becomes ill-conditioned — trust the shape, not the individual bars." · panels:
     "On the surface itself the method works for any 2-D shape: C_p converges to the exact 1 − 4 sin²θ like 1/N²."
- **Derivation tab:** **D26** (10 steps) `view: 'body'`, goal `set {tgt: airship}`; step 1 highlights a sink element; step
  5 `watch` "α runs from θ at O to α₁ at A"; step 9 **live** "Q − ak = **live closure**"; step 10 **live** "ψ = 0 surface:
  length **live L** m". **D27** (6 steps) `view: 'bars'` (phones keep `bars`); step 4 `set {tgt: rankine, N: 20}` **live**
  "one row of A at your point"; step 6 `set {tgt: sphere}` `watch` "cond explodes". Interpret `s => "N = ${N}: body error
  ${err} m, Σk_nΔξ = ${net}, cond ${cond}."`.
- **Code:**
  ```python
  s = ch06.axial_state("{{tgt}}", N={{N}}, U={{U}})
  print(s["body_error"], s["net_strength"], s["cond"])   # {{err}}, {{net}}, {{cond}}
  zb, Rb = ch06.axisym_body_target("{{tgt}}", N={{N}})
  sol = ch06.axial_singularity_solve(zb, Rb, U={{U}}, N={{N}})
  print(sol["k"][:3])          # first strengths {{k3}} m²/s
  A = ch06.airship({{U}}, 1.0, 1.0)            # (6.95) with Q = a = 1
  print(A["length"], A["closure"])             # {{L}} m, {{closure}}
  ```
- **Walkthrough (6 steps):** 1. "An airship on paper" — "The hull is given. Which hidden sources make the stream flow round
  it?" `set {tgt: airship}` · 2. "A line sink" — "A sink spread over the axis tapers the tail: ψ = (k/4π)(r − r₁) (6.94)."
  `derive: {id: 'D26', step: 8}` · 3. "Close the body" — "Source Q, sink k over length a: closed only if Q = ak."
  `readouts: ['net']` · 4. "Turn it round" — "Unknown k_n, ψ = 0 at N points: an N × N system." `set {tgt: rankine}`,
  `derive: {id: 'D27', step: 4}`, `code: {lines: [4, 5]}` · 5. "Converge" — "Raise N: the fit improves, cond climbs." `play`
  over N (the step's `enter` animates N from 4 to 60) · 6. "Your turn" — "Predict: the sphere. Easy or hard? Then try it."
  `set {tgt: sphere}`.
- **Equations:** `src3` ref 'Eq. (6.87)' `\psi=-\frac{Q}{4\pi}\cos\theta` · `ls` ref 'Eq. (6.93)'
  `\psi_{\text{sink}}=\frac{k}{4\pi}\int_0^a\cos\alpha\,d\xi` · `lsc` ref 'Eq. (6.94)' `\psi_{\text{sink}}=\frac{k}{4\pi}(r-r_1)` · `air` ref
  'Eq. (6.95)' live · `sys` "The axial system" `\psi_m=-\sum_n\frac{k_n}{4\pi}(r^m_{n-1}-r^m_n)+\tfrac12UR_m^2=0` · `close` "Closure"
  `\sum_nk_n\Delta\xi=0`.
- **Check yourself:** (1) "In the airship, make the line sink weaker than Q/a. What happens?" — "The body no longer
  closes: net outflow Q − ak must leave downstream, like the half-body." · (2) "Why do the bars add up to zero?" — "A closed
  body cannot emit fluid: Σk_nΔξ = 0 is forced by the shape, not imposed." · (3) "Why is the sphere hard?" — "Its exact
  singularity is a point doublet; smooth axial sources can only approximate it, and the matrix becomes ill-conditioned."
  `set {tgt: sphere}` · (4) "In panel mode, how does the error change from N = 16 to N = 32?" — "It falls about 4×:
  second-order convergence."
- **Selftest parity rows:** `{name: 'airship length', js: airship(1, 1, 1).length, py: 'ch06.airship(1.0, 1.0,
  1.0)["length"]', rtol: 1e-9}` · `{name: 'line sink psi', js: lineSinkPsi(0.3, -0.5, 1, 1), py:
  'ch06.line_sink_stream_function(0.3, -0.5, 1.0, 1.0)', rtol: 1e-12}` · `{name: 'axial rankine N=20 error', js:
  axialState('rankine', 20).err, py: 'ch06.axial_state("rankine_oval", N=20)["body_error"]', rtol: 1e-6}` · `{name:
  'axial rankine N=20 cond', js: axialState('rankine', 20).cond, py: 'ch06.axial_state("rankine_oval", N=20)["cond"]',
  rtol: 1e-8}` · `{name: 'net strength', js: axialState('rankine', 20).net, py: 'ch06.axial_state("rankine_oval",
  N=20)["net_strength"]', rtol: 0, atol: 1e-10}` · `{name: 'panel error N=16', js: panelErr(16), py:
  'ch06.panel_cp_error(16, body="ellipse")', rtol: 1e-8}`.
- **Fit plan:** 360×640: `body` + `bars`; the convergence view hidden (its two numbers in the status line); target chips
  wrap. Desktop: three views.

### E9 · added_mass_sphere
- **Title:** "Why does an ideal fluid resist acceleration?" · **Summary:** "Set a sphere's speed and acceleration (even at
  an angle): the speed pressure is symmetric and cancels, the acceleration pressure pushes back with half the displaced
  mass — so a bubble starts upward at 2g." · **CORE:** C15 (also N92–N104 (6.96)–(6.108), R32 (6.99), N84 (6.91), N105) ·
  **Reference:** `forced_damped_vibrations.html` (the oscillating sphere as a system + response graph with the full
  explanation panel).
- **meta:** `viz:order 9` · `viz:sections 6.9` · `viz:equations 6.97 6.98 6.99 6.104 6.105 6.106 6.108 6.109` ·
  `viz:fluidpy ch06.added_mass_state ch06.moving_sphere_surface_pressure ch06.added_mass_sphere ch06.added_mass_by_energy
  ch06.sphere_motion` · `viz:derivations D28 D29 D30 D31`.
- **Physics:** `surfVel(e, us)` (6.104) ↔ `pf.moving_sphere_velocity` on \|ξ\| = a; `surfP(e, us, dus, a, rho)` split
  steady/acceleration (6.105) ↔ `pf.moving_sphere_surface_pressure(..., split=True)`; `forceQuad(p, a)` (Gauss–Legendre 16 ×
  trapezoid 32) ↔ `pf.sphere_force_quadrature`; `M(a, rho)` ↔ `pf.added_mass_sphere`; `energyM(a, rho)` (the surface form
  of D31) ↔ `ch06.added_mass_by_energy`; `motion(mode, …)` RK4 on (m + M)du/dt = F_E (+ (ρV − m)g for ball/bubble) ↔
  `ch06.sphere_motion`; the probe state ↔ `ch06.added_mass_state`.
- **Modes:** "prescribed motion" (sliders u_s, du_s/dt, angle) · "oscillating" (x_s = A sin ωt) · "ball / bubble" (released
  from rest; density ratio chips).
- **Views** (rows [1.25, 1]):
  1. `sphere` "The sphere" — a section through the sphere in the plane of u_s and du_s/dt: the surface coloured by
     p − p∞ (steady teal shading + acceleration orange shading toggles, or the total), pressure arrows, fluid-frame
     streamlines of the dipole (6.97), arrows u_s (black) and du_s/dt (orange), the clicked surface point with its
     velocity u_a (teal); the title carries "F = … N".
  2. `ptheta` "Surface pressure vs θ_s" — steady part (teal, symmetric about 90°), acceleration part (orange, + front, −
     back), total (black) against θ_s ∈ [0°, 180°] measured from u_s (in the plane of motion); markers at the clicked
     angle.
  3. `time` "Force / velocity vs time" (`hidePortrait: true`) — oscillating: F_s(t) (black) and −M du_s/dt (dashed
     orange) on top of each other, x_s(t) faint; ball/bubble: u(t) with added mass (solid) and without (dashed), the
     initial slope labelled "2g" for the bubble.
- **Controls:** `us` "Speed $u_s$" 0…3 m/s, default 1 · `dus` "Acceleration $du_s/dt$" −5…5 m/s², default 2 · `ang`
  "Angle between them" 0…180°, default 0 · `a` "Radius $a$" 0.02…0.5 m, default 0.1 (optional) · `mode` chips · `dens`
  chips "bubble 0 · wood 0.6 · water 1 · steel 7.8" (ball/bubble mode).
- **Transport:** `t` 0 → 10 s (oscillating: 2 periods; ball/bubble: 1 s), rate 1, `end: 'loop'` (oscillating) /
  `'hold'` (release); end card (release): "Bubble: u = 2g·t at first; without added mass it would start at ∞."
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **terms** (steady force ≡ 0,
  acceleration force −M du_s/dt, and the energy route's M as a check bar) · **presets** "steady motion" {us: 1, dus: 0},
  "starting from rest" {us: 0, dus: 2}, "oscillating sphere" {mode: 'oscillating'}, "air bubble" {mode: 'release', dens:
  0}, "steel ball" {mode: 'release', dens: 7.8}, "acceleration sideways" {ang: 90} · **modes** · **status** ("⚖️ steady:
  no force (d'Alembert)" / "🚀 accelerating: the fluid pushes back with M = 2.09 kg") · **inspector** (click the surface:
  "u_a = (3/2)(u_s·e)e − ½u_s = …; p_steady = ½ρ\|u_s\|²((9/4)cos²θ_s − 5/4) = … Pa; p_acc = ρ(a/2)e·du_s/dt = … Pa").
- **Readouts:** "Added mass $M$" (kg) · "Displaced mass" (kg) · "Force $F_s$" (N) · "Accel. (release)" (m/s²).
- **Explain:**
  0. *What the views show* — "The sphere's surface is coloured by pressure: <b class="c-teal">teal</b> the part set by
     the speed, <b class="c-orange">orange</b> the part set by the acceleration. The middle view plots both round the
     surface; the third (hidden on phones) follows the force or the velocity in time."
  1. *Surface velocity (6.104)* — "at your point (θ_s = **live ths**° from u_s): u_a = (3/2)(u_s·e)e − ½u_s ⇒ \|u_a\| =
     **live ua** m/s — at the front (e ∥ u_s) the fluid moves with the sphere, u_a = u_s; at the side (e ⟂ u_s) it moves
     backwards at −½u_s."
  2. *Speed part of the pressure* — "½ρ\|u_s\|²((9/4)cos²θ_s − 5/4) = ½ × 1000 × **live us**² × (**live bracket**) =
     **live ps** Pa — symmetric front/back (6.106)."
  3. *Acceleration part* — "ρ(a/2) e·du_s/dt = 1000 × **live a**/2 × **live edu** = **live pa** Pa — positive in front
     of the acceleration, negative behind."
  4. *The forces (6.98)* — "∮ of the speed part × n = **0** (symmetric); ∮ of the acceleration part × n = −(2π/3)ρa³
     du_s/dt = **live F** N along −du_s/dt" (boxed; quadrature **live Fq** N).
  5. *Added mass (6.108)* — "M = (2π/3)ρa³ = **live M** kg = ½ × displaced **live md** kg; by kinetic energy (D31)
     ½MU² = ½ρ∫\|∇φ\|²dV ⇒ M = **live Me** kg — the same."
  6. *Newton (6.109)* — ball/bubble: "(m + M) du/dt = (ρV − m)g ⇒ du/dt = **live acc** m/s² (**live accg** g); without M:
     **live acc0** m/s²."
  7. *At the current time* — oscillating: "x_s = A sin ωt = **live xs** m, du_s/dt = −Aω² sin ωt = **live dusd**, F_s =
     **live Ft** N."
  8. *Reading the current setting* — steady: "No acceleration: only the teal part, symmetric — zero force, the same
     pressures as a still sphere in a stream (6.106)." · accelerating: "The orange part is high in front, low behind:
     the fluid being pushed aside and filled in resists — M = **live M** kg." · angle ≠ 0: "The force follows the
     acceleration, not the velocity: the speed part never contributes." · bubble: "Almost no mass: its motion is ruled by
     M, so it starts at 2g." · steel: "Heavy: M adds only 6 % to its inertia."
- **Derivation tab:** **D28** (5 steps) `view: 'sphere'`, goal `set {us: 1, dus: 0}`; step 5 `watch` "the normal velocity
  of the fluid equals the sphere's all round". **D29** (14 steps, ★★★) `view: 'ptheta'` (phones keep `ptheta`); step 8
  `watch` "at θ_s = 90° the fluid moves backwards at ½u_s"; step 12 **live** "(9/8)(u·e)² − (5/8)\|u\|² = **live ps**/ρ";
  step 13 `set {dus: 2}` `watch` "the orange curve appears". **D30** (9 steps) `view: 'sphere'`; step 2 `set {dus: 0}`
  `watch` "the steady force bar stays 0"; step 8 **live** "F = −(2/3)π × 1000 × **live a**³ × **live dus** = **live F** N".
  **D31** (10 steps, ★★★) `view: 'time'` (phones: `sphere`); step 7 **live** "(2π/3)a³U² = **live T2**"; step 10 highlights
  `term:energy` "the energy bar equals M". Interpret `s => "Your sphere: M = ${M} kg (half of ${md} kg); force ${F} N
  against the acceleration."`.
- **Code:**
  ```python
  st = ch06.added_mass_state(a={{a}}, rho=1000.0, us={{us}}, dus={{dus}},
                             angle_deg={{ang}}, theta_s_deg={{ths}})
  print(st["p_steady"], st["p_accel"])   # {{ps}} + {{pa}} Pa at your point
  print(st["F_steady"], st["F_accel"])   # 0 and −M du/dt = {{F}} N
  print(st["M"], st["M_energy"])         # {{M}} kg by pressure and by energy
  m = ch06.sphere_motion({{mb}}, {{a}}, mode="{{pymode}}", g=9.81,
                         t_eval=[0.0, 0.01])
  print(m["du_dt"][0])                   # {{acc}} m/s² at release
  ```
- **Walkthrough (7 steps):** 1. "Shaking a ball in water" — "No drag in ideal flow — so why is a ball hard to shake
  under water?" `set {mode: 'prescribed', us: 1, dus: 0}` · 2. "The moving dipole" — "A moving sphere carries a dipole
  d = 2πa³u_s (6.97)." `derive: {id: 'D28', step: 4}` · 3. "Two pressure parts" — "Speed part symmetric; acceleration part
  + in front, − behind (6.105)." `set {dus: 2}`, `derive: {id: 'D29', step: 13}` · 4. "The force" — "The teal part
  integrates to zero, the orange to −M du_s/dt (6.108)." `terms: true`, `code: {lines: [4, 4]}` · 5. "Half the displaced
  mass" — "M = (2π/3)ρa³: the fluid adds half its displaced mass. Energy says the same (D31)." `readouts: ['M']` · 6. "A
  bubble" — "Released from rest, a bubble starts at 2g — not infinity." `set {mode: 'release', dens: 0}`, `play: true` ·
  7. "Your turn" — "Predict the steel ball's start, then press ▶." `set {dens: 7.8}`.
- **Equations:** `phi` ref 'Eq. (6.97)' `\phi=-\frac{a^3}{2\lvert\boldsymbol\xi\rvert^3}\mathbf u_s\cdot\boldsymbol\xi` · `F` ref 'Eq. (6.98)' ·
  `ua` ref 'Eq. (6.104)' `\mathbf u_a=\tfrac32(\mathbf u_s\cdot\mathbf e_\xi)\mathbf e_\xi-\tfrac12\mathbf u_s` · `p` ref 'Eq. (6.105)'
  live · `st` ref 'Eq. (6.106)' · `M` ref 'Eq. (6.108)' `M=\frac{2\pi a^3\rho}{3}` live · `N` ref 'Eq. (6.109)'
  `\mathbf F_E=(m+\frac{2\pi}{3}\rho a^3)\frac{d\mathbf u_s}{dt}`.
- **Check yourself:** (1) "Set du_s/dt = 0. What is the force, whatever u_s?" — "Zero: the speed part of the pressure is
  fore–aft symmetric (6.106)." `set {dus: 0}` · (2) "Turn the acceleration 90° from the velocity. Which way does the force
  point?" — "Against the acceleration: only the acceleration part contributes." `set {ang: 90}` · (3) "Why does a bubble
  start at 2g?" — "Buoyancy ρVg accelerates only the added mass ½ρV: a = ρVg/(½ρV) = 2g." `set {mode: 'release', dens:
  0}` · (4) "Double the radius. By how much does M grow?" — "8× (a³), like the displaced mass: M/displaced stays ½."
- **Selftest parity rows:** `{name: 'M a=0.1', js: M(0.1, 1000), py: 'ch06.added_mass_sphere(0.1, 1000.0)', rtol: 1e-14}` ·
  `{name: 'M by energy', js: energyM(0.1, 1000), py: 'ch06.added_mass_by_energy(0.1, 1000.0)', rtol: 1e-10}` · `{name:
  'p front total', js: surfP(E0, P0).total, py: 'ch06.added_mass_state(a=0.1, rho=1000.0, us=1.0, dus=2.0, angle_deg=0.0,
  theta_s_deg=0.0)["p_total"]', rtol: 1e-12}` · `{name: 'p at 60 deg steady', js: surfP(E60, P0).steady, py:
  'ch06.added_mass_state(theta_s_deg=60.0)["p_steady"]', rtol: 1e-12}` · `{name: 'F accel by quadrature', js:
  forceQuad(P0).acc, py: 'ch06.added_mass_state(a=0.1, rho=1000.0, us=1.0, dus=2.0)["F_accel"]', rtol: 1e-10}` · `{name:
  'bubble start', js: motion('release', 0, 0.1).a0, py: 'ch06.sphere_motion(0.0, 0.1, mode="bubble", g=9.81,
  t_eval=[0.0, 0.01])["du_dt"][0]', rtol: 1e-8}` · invariant `{name: 'bubble = 2g', js: motion('release', 0, 0.1).a0/9.81,
  expect: 2, rtol: 1e-9}`.
- **Fit plan:** 360×640: `sphere` + `ptheta`; the time view hidden (the release acceleration in the status line: "🚀
  bubble: 19.6 m/s² = 2g"); mode chips wrap. Desktop: three views.

### B1 · flow_net_sources_vortices
- **Title:** "Where does a vortex hide its strength?" · **Summary:** "A draggable loop measures the flux of ∇ψ and ∇φ round a
  vortex and a source: the same number for every radius, zero when the loop misses the centre; ω = −∇²ψ is zero everywhere
  except in the core." · **CORE:** C02 (also N05–N14 (6.5)–(6.15), N09, R08) · **Reference:** ch02 `gauss_flux_box`
  pattern with `forced_damped_vibrations.html`'s explanation panel. **Backup only — built only if one of E1–E9 fails
  review.**
- **meta:** `viz:order 10` · `viz:sections 6.2` · `viz:equations 6.4 6.6 6.10 6.13 6.15` · `viz:fluidpy
  ch06.delta_flux_check ch06.vorticity_from_psi ch06.orthogonality_check ch06.rankine_vortex_psi` · `viz:derivations
  D02 D03 D04`.
- **Physics:** `flux(kind, strength, cx, cy, R)` (periodic trapezoid of ∇f·n) ↔ `ch06.delta_flux_check`; `omegaRankine`
  ↔ `ch06.vorticity_from_psi("rankine", …)`.
- **Views:** (1) `net` flow net (ψ solid, φ dashed) with the draggable loop; (2) `flux` flux vs loop radius (flat) with
  the current loop's dot; (3) `omega` (`hidePortrait`) ω = −∇²ψ heatmap, point vortex vs Rankine core.
- **Controls:** kind chips (source · vortex · both · Rankine) · strength · loop radius · loop centre (drag) · core a
  (Rankine).
- **Depth features:** Explain, Code, Derivation, linked views (3), presets (source, vortex, both, Rankine), inspector
  (the flux arithmetic: "∂ψ/∂r = −Γ/2πR = −0.159 m/s × 2πR = −1.000 m²/s"), status ("🎯 loop encloses the centre: flux = −Γ"
  / "⭕ loop misses it: flux 0").
- **Explain** (outline): 0 views · 1 ∂ψ/∂n on the loop with numbers · 2 flux = −Γ for any R · 3 the source's flux m · 4 the
  loop that misses the centre · 5 ω = −∇²ψ in the Rankine core (Γ/πa²) · 6 reading: "Laplace everywhere but the centre;
  the loop reads the strength."
- **Walkthrough (5 steps):** "A whirlpool and a spring" · "ψ = −(Γ/2π) ln r (6.8) is harmonic for r > 0" · "the loop's flux
  is −Γ (6.6)" · "move the loop off the centre: 0" · "your turn: the Rankine core".
- **Check yourself (3):** flux independent of R (why) · zero when the loop misses the centre · inside a Rankine core the
  flux grows like r² (vorticity spread out).
- **Selftest parity rows:** `{name: 'vortex flux R=2', js: flux('vortex', 1, 0, 0, 2), py: 'ch06.delta_flux_check("vortex",
  radii=[2.0], Gamma=1.0)[0]', rtol: 1e-12}` · `{name: 'source flux R=0.3', js: flux('source', 2, 0, 0, 0.3), py:
  'ch06.delta_flux_check("source", kind="phi", radii=[0.3], m=2.0)[0]', rtol: 1e-12}` · `{name: 'Rankine omega centre',
  js: omegaRankine(0, 0, 1, 0.1), py: 'ch06.vorticity_from_psi("rankine", 0.0, 0.0, Gamma=1.0, a=0.1)', rtol: 1e-6}`.
- **Fit plan:** as E1 (two views on phones, the ω view hidden).

## Part D — runtime budget (full run < 5 min on a laptop / Colab CPU)

Chapter 6 is mostly closed forms in complex arithmetic (µs per point); the costs are the flow-net figures (contours of ψ
and φ on grids ≤ 241 × 241, several panels each), Newton stagnation searches, the plotly sliders that precompute every
step (contour polylines via `contourpy`), four animations, the Laplace solves (pure-Python Gauss–Seidel only on small
grids), the axial method (dense solves N ≤ 80), the panel convergence (N ≤ 128), `quad`-based energy integrals and five
★★★ sympy checks (D17, D18, D21, D29, D31). The ch05 notebook ran in 67 s alone; this one has 15 CORE blocks, 31
derivations and 9 explainers.

| Section | Heaviest cells | Full | FAST (`FLUIDPY_FAST=1`) |
|---|---|---|---|
| setup + imports | numpy/scipy/sympy/plotly/pint, the chapter modules | 9 s | 9 s |
| §6.1 C01 | residuals at points (stencils), applicability table, two-panel residual map (61 × 41 stencil evaluations × 4 terms) | 4 s | 3 s (41 × 27) |
| §6.2 C02 | residual table, flux checks (5 radii × 256 points), Rankine ω (121² 9-point stencil), from-scratch 5-point + loop, three flow nets (3 × 241² w evaluations), flux-vs-radius, ω heatmap, sympy polar checks | 6 s | 4 s (161² grids) |
| §6.2 C03 | Flow sum, Newton stagnation, far-field slopes (4 radii × 256), three-panel figure | 3 s | 2 s |
| §6.3 C04 | harmonic polynomials (sympy degree 3), pair error (5 ε × 3 points), element gallery (5 × 201²), **frames animation 12 frames** (161² ψ per frame, cached by `doublet_limit_frames`) | 8 s | 5 s (8 frames, 121²) |
| §6.3 C05 | half-body numbers, `brentq`, net force (3 lengths × 4001 points), two-panel figure, **slider 20 steps** (contours of 161² per step → polylines) | 7 s | 4 s (10 steps, 121²) |
| §6.3 C06 | surface forces (64-point trapezoid), arrows + C_p figure, qualitative-band slider 13 steps | 3 s | 3 s |
| §6.3 C07 | stagnation points, Newton, family check (4 Γ), 2 × 2 regime figure (4 × 201²), **slider 21 steps** (C_p curves: cheap) | 5 s | 4 s (11 steps) |
| §6.3 C08 | Example 6.1 closed + numeric routes, `point_vortex_evolve`, two figures, **video animation 60 frames** (201 wall points + 121² streamlines per frame) | 14 s | 8 s (30 frames, 81²) |
| §6.4 C09 | CR residuals (5 flows), corner nets (4 × 201²), log–log slopes, **slider 22 steps** (analytic polylines r(θ): cheap) | 5 s | 4 s |
| §6.5 C10 | Blasius on 3 + 2 contours (n = 256), Laurent FFT, **D17 sympy ★★★** (cylinder surface identity, ∮ const dz* and the Blasius integral by residues: ≈ 2 s), **D18 sympy ★★★** (series square, coefficients, residue: ≈ 1 s), contour-force figure (R/a 1–20 × 2 bodies × n = 256), Laurent bars | 8 s | 7 s |
| §6.6 C11 | ellipse check, inverse branches (4 points), mapped flow checks (normal velocity on 400 points, far field), **D21 sympy/wrong-variant ★★★** (roots, product b², principal-root failure on a 41 × 41 grid in the left half-plane: < 1 s), two-plane figure (2 × 241² with the inverse), grid images (3 × 201²), **slider 21 steps** (mapped streamlines: contours of ψ(ζ(z)) on 161² per step) | 10 s | 6 s (11 steps, 121²) |
| §6.7 C12 | 5-point order (3 grids), 4-point system, three solvers to 1e-10 (tiny), `example_6_2` book grid Gauss–Seidel (≈ 0.2 s) + refined ×2, ×4 by `spsolve`, from-scratch double loop, two-panel figure (3 residual histories), **frames animation 24 frames** | 9 s | 6 s (12 frames; refine ×1, ×2) |
| §6.8 C13 | Stokes residuals, flux quadrature, sympy spherical checks (`core.curvilinear`, cached), sphere figure (meridian ψ on 201²), cylinder-vs-sphere figure, slider 17 steps (analytic curves) | 6 s | 5 s |
| §6.8 C14 | line sink closed vs `quad` (2 points), airship (`brentq` × 2 + ψ = 0 contour on 201²), collocation primer, axial solves N = 10 … 80 (dense, ms each), sphere target, from-scratch double loop (N = 20), panels N = 8 … 64 (+128 in the figure), three-panel figure, **slider 16 steps** (a solve per step) | 8 s | 5 s (panels to 64, slider 8 steps) |
| §6.9 C15 | added mass three routes (`added_mass_by_energy(method="volume")` with `dblquad`: ≈ 1 s), force quadrature 24 × 48, from-scratch sum, **D29 sympy ★★★** (symbolic motion x_s(t), direct ∂φ/∂t vs the chain, gradient on the surface, the printed bracket: ≈ 3 s), **D31 sympy ★★★** (volume integral over r > a in spherical coordinates vs the surface form: ≈ 2 s), two-panel figure (`sphere_motion` × 4), **video animation 60 frames**, density slider 17 steps (17 × 2 `solve_ivp` runs), live cell (not executed on the page) | 16 s | 11 s (30 frames, 9 steps) |
| §6.10 + summary | markdown | 0 s | 0 s |
| explainers (9 × `show_viz`) | read the HTML files | 1 s | 1 s |
| **Total** | | **≈ 132 s** | **≈ 97 s** |

**FAST plan.** Every size-dependent choice is written `a if not FAST else b` in the cell: flow-net grids 241² → 161²
(figures) and 161² → 121² (sliders, animations); animations 60 → 30 frames (videos), 24/12 → 12/8 (frames players);
slider figures ≤ 22 steps → ≤ 11 (≤ 4 traces × ≤ 400 points each, < 250 kB); Example 6.2 refinement ×4 → ×2; panels to
N = 64; density slider 17 → 9 steps. **Cached arrays / results:** `ch06.doublet_limit_frames` caches its ψ grids
(`functools.lru_cache` on the ε tuple) for the C04 animation and figure; the C05 half-body grid is built once and reused
by its figure and slider; `ch06.example_6_1_wall_pressure` is vectorised over wall points and times (one call for the
whole C08 animation); Example 6.2's Gauss–Seidel history is computed once and reused by the C12 figure, animation and
the E7 ghost export; the sympy results of D17, D18, D21, D29, D31 and the `core.curvilinear` spherical simplifications
are cached with `lru_cache` in `ch06.kutta_zhukhovsky_sym` and `ch06.moving_sphere_dphidt_sym`. **Complex arithmetic
budget:** every flow-net cell evaluates `Flow.w` once on the grid (complex, vectorised) and takes Re/Im for φ/ψ —
never separate φ and ψ calls. **Outputs:** 2 videos (C08, C15) + 2 frames players (C04, C12), 9 plotly sliders, ≈ 27
static figures — page well under 15 MB (each video < 2 MB at dpi 80). Sympy cells never `simplify` expressions with more
than two generic functions at once; D29's check uses a concrete polynomial motion x_s(t) and substitutes before
simplifying.

---

## Part E — prerequisite ledger
Every concept, symbol, maths tool and Python function or idiom the notebook or its explainers use, with where it is
explained. "primer (in Cxx)" = a 📎 primer placed in that block before first use (the primer term is the Concept text —
the builder uses it verbatim in `nb.primer`); "knowledge/primers.md: <term> (chNN Pnn) — reminder" = a one-line reminder
naming the earlier primer; a CORE/RECAP id alone = taught there; "Cxx (Nnn)" = the NOTE of this chapter placed in that
block; "Cxx (gloss …)" = one sentence where it is used; "Cxx (P1nn …)" = explained by this chapter's primer of that
number, named without repeating its title. Earlier chapters' material that is neither a primer nor a ch06 RECAP is named
by section ("Ch. 4 §4.9"). New primers P149–P164 in first-use order: 2-D divergence theorem and the Dirac delta in the
plane (C02), a limit with a product held fixed (C04), integrals of sines and cosines over a full period (C06), Newton's
method for complex zeros (C07), the complex plane in numpy, complex derivative and analytic functions, complex
logarithm, powers and branch cuts (C09), Cauchy's integral theorem, Laurent series and residues, sympy for complex
series and residues (C10), complex square roots and the quadratic formula (C11), iterative solvers: Jacobi,
Gauss–Seidel, SOR, boolean masks and scipy.sparse (C12), collocation and the condition number (C14), surface integrals
on a sphere, Green's first identity (C15) — 16 primers.

| Concept | First used in | Explained by |
|---|---|---|
| partial derivative ∂/∂x | C01 | knowledge/primers.md: partial derivative (ch01 P25) — reminder |
| Laplacian ∇² | C01 | C01 (reminder; Ch. 2 §2.9) |
| vorticity ω = ∇×u | C01 | C01 (D01; Ch. 3 §3.4) |
| continuity (4.7) and (4.10) | C01 | C01 (D01 step 1; Ch. 4 §4.2) |
| product rule for a divergence | C01 | knowledge/primers.md: product rule for a divergence (ch04 P113) — reminder |
| Navier–Stokes (4.38), (4.39b) | C01 | C01 (D01 step 2; Ch. 4 §4.6) |
| hydrostatic split p = p_h + p′ | C01 | C01 (D01 step 3; Ch. 4 §4.9) |
| curl of a curl identity | C01 | knowledge/primers.md: curl of a curl identity (ch04 P122) — reminder |
| viscous force −μ∇×ω (4.40) | C01 | C01 (D01 step 5; Ch. 4 §4.6) |
| ideal-flow equations (6.1) | C01 | C01 |
| Kelvin's circulation theorem (5.8) | C01 | R01 |
| Mach number M = U/c (4.111) | C01 | R01 |
| baroclinic torque ∇ρ×∇p/ρ² | C01 | R01 (Ch. 5 §5.2) |
| Reynolds number Re = ρUL/μ (4.103) | C01 | C01 (N01) |
| boundary layer, thickness L/√Re | C01 | C01 (idea, N01) |
| separation and wake | C01 | C01 (N03) |
| no-slip vs no-through-flow | C01 | C01 (N02, D01 step 8) |
| boundary conditions | C01 | knowledge/primers.md: boundary conditions (ch01 P20) — reminder |
| viscous stress vs net viscous force | C01 | C01 (tiny example, ⚠️ callout) |
| plane Poiseuille flow | C01 | C01 (code explain; Ch. 4 §4.6) |
| central differences at a point | C01 | knowledge/primers.md: finite differences (ch01 P21) — reminder |
| Python dictionaries (fluidpy results) | C01 | knowledge/primers.md: Python dictionaries (ch01 P23) — reminder |
| functions as arguments and lambda | C01 | knowledge/primers.md: functions as arguments and lambda (ch01 P29) — reminder |
| assert np.allclose | C01 | knowledge/primers.md: assert np.allclose (ch01 P15) — reminder |
| np.meshgrid and the grid layout | C01 | knowledge/primers.md: np.meshgrid and the project grid layout (ch02 P76) — reminder |
| contour and streamline plots | C01 | knowledge/primers.md: plt.contour, plt.quiver and plt.streamplot (ch02 P78) — reminder |
| logarithmic colour and axes | C01 | knowledge/primers.md: power laws and log–log plots (ch01 P13) — reminder |
| f-strings and printing | C01 | knowledge/primers.md: f-strings (ch01 P04) — reminder |
| ch06.ideal_flow_residuals, ideal_flow_applicability | C01 | C01 (code explain) |
| ch05.kelvin_hypotheses_text, ch04.is_incompressible_regime | C01 | R01 (code comment) |
| 2-D continuity (6.2) | C02 | R02 |
| stream function ψ (6.3) and streamline slope | C02 | R03 |
| Δψ = volume flow per unit depth | C02 | R03 |
| 2-D irrotationality (6.9) | C02 | R04 |
| velocity potential φ (6.10), simply connected region | C02 | R04, C02 (N09) |
| polar continuity (6.19) | C02 | R06 |
| polar irrotationality (6.20) | C02 | R07 |
| polar Laplacian (6.23a), (6.23b) | C02 | R08, R09 |
| sympy diff, simplify, symbols | C02 | knowledge/primers.md: sympy (ch01 P40) — reminder |
| ω_z = −∇²ψ (6.4) | C02 | C02 (D02) |
| Laplace's equation and harmonic functions (6.5), (6.12) | C02 | C02 (N05, N11) |
| 2-D divergence theorem and the Dirac delta in the plane | C02 | primer (in C02) |
| flux of ∇ ln r through a circle = 2π | C02 | C02 (D03 steps 3–4) |
| point vortex as a delta source (6.6) | C02 | C02 (N06, D03) |
| point source as a delta source (6.13) | C02 | C02 (N12, D03) |
| Poisson equation (6.11) | C02 | C02 (N10); knowledge/primers.md: Poisson equation and Green's function (ch05 P139) — reminder |
| natural logarithm, d(ln r)/dr = 1/r | C02 | knowledge/primers.md: natural logarithm and exponential (ch01 P36) — reminder |
| uniform flow (6.7), (6.14) | C02 | C02 (N07, N13) |
| vortex ψ (6.8), Γ counterclockwise | C02 | C02 (N08), front-matter conventions cell |
| source φ (6.15), m in m²/s | C02 | C02 (N14) |
| total differential dψ | C02 | R03 (Ch. 4 §4.3) |
| perpendicular slopes, dot product | C02 | C02 (D04; Ch. 2 §2.2) |
| level sets | C02 | knowledge/primers.md: level sets and the directional derivative (ch02 P75) — reminder |
| equipotentials ⟂ streamlines, ∇ψ = e_z × ∇φ | C02 | C02 (D04, N09) |
| polar velocities (6.21), (6.22) | C02 | C02 (N18, N19) |
| chain rule | C02 | knowledge/primers.md: chain rule (ch01 P49) — reminder |
| Rankine vortex | C02 | C02 (code explain; Ch. 3 §3.5) |
| pf.Uniform, pf.Source, pf.Vortex objects | C02 | C02 (code explain) |
| dataclasses and named results | C02 | knowledge/primers.md: dataclasses and named results (ch04 P111) — reminder |
| pf.laplacian_residual, ch06.delta_flux_check, vorticity_from_psi | C02 | C02 (code explain) |
| numpy slicing for stencils | C02 | knowledge/primers.md: numpy broadcasting (ch02 P77) — reminder |
| np.hypot | C02 | C02 (gloss in a code comment: √(x² + y²)) |
| periodic trapezoid rule on a closed loop | C02 | knowledge/primers.md: periodic trapezoid rule on a closed loop (ch05 P136) — reminder |
| stagnation point | C03 | C03 (tiny example; Ch. 3 §3.3) |
| Bernoulli everywhere (6.18) | C03 | R05 |
| unsteady Bernoulli (4.75) | C03 | R05 |
| core.bernoulli.bernoulli_function | C03 | R05 (code comment) |
| superposition (linearity of Laplace) | C03 | C03 |
| directional derivative ∂φ/∂n | C03 | knowledge/primers.md: level sets and the directional derivative (ch02 P75) — reminder |
| unit normal and tangent of a curve | C03 | knowledge/primers.md: parametric curves, tangent vector and arc length (ch03 P92) — reminder |
| no-through-flow condition (6.16) | C03 | C03 (D05, N16) |
| far-field condition (6.17) | C03 | C03 (N17) |
| dividing streamline and body | C03 | C03 (figure), C05 |
| complex numbers as points x + iy in code (`-1+0j`) | C03 | C03 (gloss; full treatment C09 P153); knowledge/primers.md: square root of a negative number (ch01 P45) — reminder |
| Flow.stagnation_points (u = v = 0 by iteration) | C03 | C03 (code explain; method in C07 P152) |
| observed order from a log–log slope | C03 | knowledge/primers.md: power laws and log–log plots (ch01 P13) — reminder |
| corner flow ψ = 2Axy (6.24) | C04 | R10 |
| vortex velocities, u_θ = Γ/2πr (5.2) | C04 | R11 |
| harmonic polynomials | C04 | C04 (N20) |
| rotation of axes by 45° | C04 | C04 (N21–N23; Ch. 2 §2.2) |
| relatives (6.25)–(6.27) | C04 | C04 (N21–N23) |
| source velocities | C04 | C04 (N24) |
| source–sink pair (6.28), dipole vector | C04 | C04 (N25) |
| a limit with a product held fixed | C04 | primer (in C04) |
| first-order Taylor expansion, ln(1 + s) ≈ s | C04 | knowledge/primers.md: first-order Taylor expansion (ch01 P26) — reminder |
| log rules ln√s = ½ ln s, ln(ab) = ln a + ln b | C04 | C04 (gloss) |
| limits and orders of smallness O(ε²) | C04 | knowledge/primers.md: limits and orders of smallness (ch02 P68) — reminder |
| doublet (6.29) | C04 | C04 (D06) |
| singular points of elements | C04 | C04 (N26) |
| animate and show_animation | C04 | knowledge/primers.md: animate and show_animation (ch01 P16) — reminder |
| half-body (6.30), (6.31) | C05 | C05 (D07, N27) |
| np.arctan2 | C05 | knowledge/primers.md: np.arctan2 (ch02 P70) — reminder |
| angle range [0, 2π) and the cut inside the body | C05 | C05 (gloss) |
| mass balance 2h_max U = m | C05 | C05 (D07 step 9; Ch. 4 §4.2) |
| pressure coefficient C_p (6.32) | C05 | C05 (N28; Ch. 4 §4.11) |
| half-body surface C_p and its zero | C05 | C05 (D08, N29) |
| scipy.optimize.brentq | C05 | knowledge/primers.md: scipy.optimize.brentq (ch03 P108) — reminder |
| slider_figure | C05 | knowledge/primers.md: slider_figure (ch01 P17) — reminder |
| show_viz | C05 | knowledge/primers.md: show_viz (ch01 P18) — reminder |
| contourpy polylines for plotly | C05 | C05 (code comment) |
| Rankine oval | C05 | C05 (what would change if), E1 preset |
| cylinder = stream + doublet (6.33) | C06 | R12 |
| cylinder velocity (6.34) | C06 | R13 |
| moving cylinder, Galilean frames | C06 | R14 |
| cylinder surface C_p (6.35) | C06 | C06 (D09, N30) |
| integrals of sines and cosines over a full period | C06 | primer (in C06) |
| net pressure force −∮p n dA | C06 | knowledge/primers.md: net force from pressure (ch01 P28) — reminder |
| trapezoid rule, np.trapezoid | C06 | knowledge/primers.md: trapezoid rule (ch01 P37) — reminder |
| d'Alembert's paradox | C06 | C06 |
| drag and lift per unit depth D, L | C06 | C06, C10 (N53) |
| measured vs ideal C_p, angle from the front | C06 | C06 (N31) |
| clockwise vs counterclockwise Γ | C07 | C07 (⚠️ callout), front-matter conventions cell |
| cylinder with circulation (6.36) | C07 | C07 (N32) |
| surface speed (6.37) | C07 | C07 (N33, D10) |
| stagnation points (6.38) and the free point | C07 | C07 (D10, N34) |
| two solutions of sin θ = s | C07 | C07 (gloss) |
| Vieta's formulas (product of roots) | C07 | knowledge/primers.md: Vieta's formulas (ch02 P71) — reminder |
| Newton's method for complex zeros | C07 | primer (in C07) |
| surface pressure (6.39) | C07 | C07 (N35) |
| lift integral with n = e_r, dl = a dθ | C07 | C07 (N36, D11) |
| lift L = ρUΓ (6.40) | C07 | C07 (D11) |
| circulation Γ = ∮u·dx (3.18) | C07 | C07 (N38; Ch. 3 §3.4) |
| uniqueness in multiply connected regions, Kutta condition (named) | C07 | C07 (N38) |
| Magnus effect | C07 | C07 (N37) |
| images for circles, circle theorem | C08 | R15 |
| vortex beside a wall | C08 | R16 |
| drift Γ/4πh | C08 | R17 |
| core.biot_savart.point_vortex_evolve | C08 | R17 (code comment; Ch. 5 §5.7) |
| odd and even functions under a mirror | C08 | C08 (gloss) |
| image rules | C08 | C08 (D12, N39) |
| two sources (6.41) and its streamlines | C08 | C08 (N40) |
| tangent addition formula | C08 | C08 (N40, gloss) |
| derivative of tan⁻¹(Y/X) with moving Y | C08 | C08 (gloss) |
| unsteady Bernoulli (4.83), ∂φ/∂t at a fixed point | C08 | C08 (D13), R05 |
| Example 6.1 wall pressure | C08 | C08 (D13) |
| the complex plane in numpy | C09 | primer (in C09) |
| np.abs, np.angle, np.conj, 1j | C09 | C09 (P153) |
| Euler's formula e^{iθ} = cos θ + i sin θ | C09 | knowledge/primers.md: square root of a negative number (ch01 P45) — reminder |
| complex conjugate | C09 | knowledge/primers.md: complex conjugate (ch02 P81) — reminder |
| complex derivative and analytic functions | C09 | primer (in C09) |
| complex potential w (6.42) | C09 | C09 |
| z = x + iy = re^{iθ} (6.43) | C09 | C09 (N41) |
| Cauchy–Riemann (6.44) | C09 | C09 (D14, N42) |
| complex velocity dw/dz = u − iv (6.45) | C09 | C09 (D14, N43) |
| Schwarz's theorem | C09 | knowledge/primers.md: Schwarz's theorem (ch04 P121) — reminder |
| complex logarithm, powers and branch cuts | C09 | primer (in C09) |
| corner flow w = Azⁿ (6.46) | C09 | C09 (D15, N44) |
| element complex potentials (6.47)–(6.53) | C09 | C09 (N45–N51) |
| np.polyfit for a slope | C09 | knowledge/primers.md: power laws and log–log plots (ch01 P13) — reminder |
| control-volume momentum (6.54) | C10 | R18 |
| momentum flux through a surface | C10 | knowledge/primers.md: momentum flux through a surface (ch04 P114) — reminder |
| force on body vs on fluid, span B | C10 | C10 (N53) |
| outward normal of a counterclockwise contour | C10 | C10 (gloss) |
| force components (6.55), (6.56) | C10 | C10 (D16, N54, N55) |
| Cauchy's integral theorem | C10 | primer (in C10) |
| complex force (6.57)–(6.59) | C10 | C10 (D17, N56–N59) |
| Blasius's theorem (6.60) | C10 | C10 (D17, N60) |
| Laurent series and residues | C10 | primer (in C10) |
| residue theorem ∮z⁻ⁿdz = 2πiδ_{n1} | C10 | C10 (N63) |
| FFT of samples on a circle | C10 | C10 (P157); knowledge/primers.md: Fourier modes and the FFT Poisson solver (ch05 P142) — reminder |
| far field of any body | C10 | C10 (N61) |
| sympy for complex series and residues | C10 | primer (in C10) |
| Kutta–Zhukhovsky theorem (6.61), (6.62) | C10 | C10 (D18, N62) |
| induced drag (named) | C10 | C10 (what would change if; Ch. 14) |
| conformal map δw = (dw/dz)δz (6.63) | C11 | C11 (N64) |
| angle preservation (6.64) and critical points | C11 | C11 (D19, N65) |
| flow nets as images of a grid | C11 | C11 (N66) |
| Zhukhovsky transformation (6.65) | C11 | C11 (D20) |
| ellipse, semi-axes and foci | C11 | C11 (D20, gloss) |
| cos² + sin² = 1, e^{iθ} + e^{−iθ} = 2 cos θ | C11 | C11 (D20 steps 2 and 5) |
| complex square roots and the quadratic formula | C11 | primer (in C11) |
| inverse Zhukhovsky map (6.69) and its branch | C11 | C11 (D21, N70) |
| derivative of an inverse function | C11 | C11 (D21 step 10) |
| circle flow with circulation (6.68) | C11 | C11 (N69) |
| grid points ψ_{i,j}, half-point differences | C12 | R19 |
| second differences (6.70), (6.71) | C12 | R20, R21 |
| five-point average rule (6.72) | C12 | C12 (D22) |
| truncation error (Δx²/12)ψ_xxxx | C12 | C12 (D22 step 3) |
| discrete maximum principle | C12 | C12 (D22 step 6) |
| four-point system (6.73) as Aψ = b | C12 | C12 (N72, D23) |
| diagonal dominance | C12 | C12 (D23 step 3) |
| iterative solvers: Jacobi, Gauss–Seidel, SOR | C12 | primer (in C12) |
| spectral radius | C12 | C12 (D23 step 6); knowledge/primers.md: eigenvalues and eigenvectors (ch02 P80) — reminder |
| residual vs change as a stopping rule | C12 | C12 (D23 step 8, N73) |
| np.linalg.solve | C12 | knowledge/primers.md: np.linalg.solve (ch01 P57) — reminder |
| boolean masks and scipy.sparse | C12 | primer (in C12) |
| Example 6.2, re-entrant corner | C12 | C12 (N74) |
| two stream functions in 3-D | C13 | R22 |
| axisymmetric continuity (6.74) | C13 | R23 |
| Stokes stream function (6.75) | C13 | R24 |
| azimuthal vorticity (6.76) | C13 | R25 |
| axisymmetric potential (6.79) | C13 | R26 |
| cylindrical and spherical coordinates (6.81) | C13 | R27 |
| cylindrical and spherical unit vectors | C13 | knowledge/primers.md: cylindrical and spherical unit vectors (ch03 P88) — reminder |
| spherical continuity (6.82) with the 1/r² slip | C13 | R28 |
| spherical vorticity (6.84) | C13 | R29 |
| spherical Laplace (6.85) | C13 | R30 |
| core.curvilinear (sympy operators) | C13 | R30 (code comment; Ch. 4 App. B) |
| cross products of cylindrical unit vectors | C13 | C13 (gloss) |
| Stokes field equation (6.77), E² operator | C13 | C13 (D24, N75) |
| flux 2πdψ between stream surfaces (6.78) | C13 | C13 (N77) |
| axisymmetric Laplace (6.80) | C13 | C13 (N78) |
| spherical velocities (6.83) | C13 | C13 (N79) |
| flux of a 3-D point source 4πr²u_r = Q | C13 | C13 (gloss) |
| 3-D elements (6.86)–(6.88) | C13 | C13 (D25, N80–N82) |
| Taylor series in 3-D (1/\|x ∓ εe_z\|) | C13 | knowledge/primers.md: multivariable first-order Taylor expansion (ch03 P98) — reminder |
| sphere flow (6.89), (6.90), C_p (6.91) | C13 | C13 (D25, N83, N84) |
| coordinate-free potential (6.92) | C13 | C13 (N85) |
| Hill's spherical vortex (exterior) | C13 | C13 (N83; Ch. 5 §5.4) |
| closed body needs zero net source | C14 | C14 (N86) |
| line sink (6.93), (6.94) | C14 | C14 (D26, N87, N88) |
| substitution in an integral | C14 | knowledge/primers.md: substitution in an integral (ch03 P106) — reminder |
| substitution z − ξ = R cot α | C14 | C14 (gloss) |
| airship (6.95) | C14 | C14 (D26, N89) |
| scipy.integrate.quad and dblquad | C14 | knowledge/primers.md: scipy.integrate.quad and dblquad (ch03 P87) — reminder |
| collocation and the condition number | C14 | primer (in C14) |
| np.vander, np.linalg.cond | C14 | C14 (P162) |
| axial singularity method | C14 | C14 (D27) |
| source panels (our extension) | C14 | C14 (N90) |
| 3-D potential, w as the z-velocity | C15 | R31 |
| unsteady Bernoulli between surface and far field (6.99) | C15 | R32 |
| moving-sphere set-up, x_s, u_s | C15 | C15 (N92) |
| Galilean frames | C15 | C15 (reminder; Ch. 3 §3.3), R14 |
| moving-sphere potential (6.96), (6.97) | C15 | C15 (D28, N93, N94) |
| kinematic condition on a moving body | C15 | C15 (D28 step 5), C03 (N16) |
| force on the sphere (6.98) | C15 | C15 (N95) |
| multivariable chain rule along a path | C15 | knowledge/primers.md: multivariable chain rule along a path (ch03 P91) — reminder |
| a function of x − x_s: ∂/∂x_s = −∇ | C15 | C15 (gloss) |
| gradient of u·ξ/\|ξ\|³ | C15 | C15 (gloss); knowledge/primers.md: gradient of 1/distance with respect to the source point (ch05 P140) — reminder |
| surface pressure (6.100)–(6.105) | C15 | C15 (D29, N96–N101) |
| steady limit (6.106) and θ_s | C15 | C15 (N102) |
| surface integrals on a sphere | C15 | primer (in C15) |
| added mass (6.107), (6.108) | C15 | C15 (D30, N103, N104) |
| Newton's law with added mass (6.109) | C15 | C15 |
| apparent-mass tensor M_ij | C15 | C15 (N104) |
| Green's first identity | C15 | primer (in C15) |
| divergence theorem in 3-D | C15 | C15 (D31 step 3; Ch. 2 §2.12) |
| kinetic energy of the fluid ½MU² | C15 | C15 (gloss) |
| buoyancy ρVg | C15 | C15 (tiny example; Ch. 1 §1.7) |
| Gauss–Legendre quadrature | C15 | knowledge/primers.md: Gauss–Legendre quadrature in 3-D and a smoothed kernel (ch05 P143) — reminder |
| scipy.integrate.solve_ivp | C15 | knowledge/primers.md: scipy.integrate.solve_ivp (ch01 P31) — reminder |
| live widgets | C15 | knowledge/primers.md: live widgets (ch01 P47) — reminder |
| symbol Γ circulation [m²/s], both signs | C02 | C02 (N08), C07 (⚠️ callout) |
| symbol m (2-D source) [m²/s] vs body mass m | C02 | C02 (N14), C15 |
| symbol d (2-D dipole [m³/s], 3-D [m⁴/s]) | C04 | C04 (N25), C13 (N82) |
| symbols U, V (stream components) | C02 | C02 (N07) |
| symbol A (corner strength) | C04 | R10 |
| symbol n (corner exponent) vs unit normal n | C09 | C09 (N44), front-matter conventions cell |
| symbol α (corner angle, conformal angle, line-sink angle) | C09 | C09 (N44), C11 (D19), C14 (D26) |
| symbol a (radius, stagnation distance, image distance, line length, circle) | C05 | C05, front-matter conventions cell |
| symbols h, h_max | C05 | C05 (D07) |
| symbols ρ, p∞ | C03 | R05 |
| symbols z, ζ, w (complex) | C09 | C09 (N41), C11 (D20) |
| symbol w as z-velocity in §6.9 | C15 | R31 |
| symbol b (Zhukhovsky constant) | C11 | C11 (D20) |
| symbol B (span) | C10 | C10 (N53) |
| symbols ψ_{i,j}, Δx, Δy | C12 | R19 |
| symbol ω (over-relaxation factor) | C12 | C12 (N73) |
| symbol ρ (spectral radius) vs density ρ | C12 | C12 (D23 step 6) |
| symbols R, r, θ, φ (azimuth) in §6.8 | C13 | R27 |
| symbols Q (3-D source) and k (line density) | C13 | C13 (N81), C14 (N87) |
| symbols M (added mass), not Mach | C15 | C15 (N104), front-matter conventions cell |
| symbols u_s, x_s, ξ, e_ξ, θ_s | C15 | C15 (N92, D29) |
| symbols F_E, F_s | C15 | C15 (N92, N95) |
| explainer controls, presets, walkthrough and tabs | C05 | C05 (the explainer cell's "What to try"; Ch. 1 §1.4 explainer guide) |

---

## Part F — derivation storyboards

Builders copy these word for word into `nb.derivation(key, title, goal=…, start=(tex, plain), plan=[…], uses=[…],
steps=[dict(did, tex, why, plain)], result=(tex, plain), interpret=…, check=…, check_src=…)` and into the explainer's
`derivations: [...]` (phones may shorten *why* to its first sentence; `live`, `set` and `watch` are the explainer's).
Every step is one move; *why* names the rule and says why we make the move (≤ 35 words, ≥ 6); *did* ≤ 8 words. The
book's own moves were read on the rendered pages (p225 for D01, p227 for D02, none for D03 (Exercises 6.1b, 6.3b), p228
for D04, p229 for D05, p233 for D06, p234 for D07, none for D08 (the book only plots Fig. 6.8), p235–p236 for D09,
p237–p238 for D10, p238–p239 for D11, p240 for D12, p242–p243 for D13, p243–p244 for D14, p244–p245 for D15, p246 for
D16, p247 for D17, p248 for D18, p249–p250 for D19, p251 for D20, p252 for D21, p253–p254 for D22, p254–p255 for D23,
p258 for D24, p259–p260 for D25, p262 for D26, p263 for D27, p264–p265 for D28, p265–p266 for D29, p266–p267 for D30,
none for D31 (Exercise 6.49)); the gaps listed in `analysis/ch06.md` §2b are filled and the notebook says so ("the book
skips this move"). Every equation named by number is written out. Colours: stream blue, source teal, sink rose, vortex
amber, doublet purple (D05–D11); steady pressure teal, acceleration pressure orange (D29–D31); lift amber, drag muted
(D11, D16–D18). (No line in this part starts with a table bar; absolute values are written with \lvert \rvert or in
words.)

### D01 · The ideal-flow equations (6.1) from continuity and Navier–Stokes — ★★, 8 steps, in C01 (notebook)
- **Goal.** Show that a constant-density flow without vorticity obeys ∇·u = 0 and ρDu/Dt = −∇p *(6.1)* — even though the
  fluid is viscous — and see which boundary condition this costs.
- **Start.** $\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0$ *(4.7)* and $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf
  g+\mu\nabla^2\mathbf u$ *(4.39b, the constant-μ, ∇·u = 0 form of 4.38)* — *in words:* mass is conserved, and mass ×
  acceleration equals pressure, weight and viscous forces.
- **Plan.** (1) Use constant density in continuity. (2) Absorb gravity into a hydrostatic pressure. (3) Rewrite the
  viscous force with the curl-of-a-curl identity. (4) Set ω = 0 and read off (6.1) and its boundary condition.
- **Tools.** product rule for a divergence ∇·(ρu) = u·∇ρ + ρ∇·u (ch04 P113, reminder) · hydrostatic balance ∇p_h = ρg
  (Ch. 4 §4.9) · curl of a curl ∇×(∇×u) = ∇(∇·u) − ∇²u (ch04 P122, reminder) · Kelvin's theorem (R01).
- **Assumptions.** ρ constant (low Mach number, no stratification: step 1); gravity conservative (step 3); ω = 0 in the
  region (step 6); constant μ (step 2).
- **Steps.**
  1. *did:* Use constant density in continuity · *tex:* $\frac{D\rho}{Dt}+\rho\nabla\cdot\mathbf u=0\ \Rightarrow\ \nabla\cdot\mathbf u=0$ ·
     *why:* The product rule splits ∇·(ρu) into u·∇ρ + ρ∇·u, and ∂ρ/∂t + u·∇ρ = Dρ/Dt = 0 for constant ρ; divide by ρ ≠ 0.
     · *plain:* Fluid of fixed density cannot be squeezed: what flows into a point flows out.
  2. *did:* Write the constant-μ momentum equation · *tex:* $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ ·
     *why:* With ∇·u = 0 and constant μ the general Navier–Stokes (4.38) reduces to (4.39b) (Ch. 4 §4.6); we start from the
     simplest correct form. · *plain:* Mass times acceleration equals pressure push, weight and viscous force.
  3. *did:* Absorb gravity into hydrostatic pressure · *tex:* $p=p_h+p',\ \nabla p_h=\rho\mathbf g\ \Rightarrow\ \rho\frac{D\mathbf u}{Dt}=-\nabla p'+\mu\nabla^2\mathbf u$ ·
     *why:* With constant ρ and conservative g a hydrostatic pressure p_h exists with ∇p_h = ρg; subtracting it removes
     gravity from the equation (the book's "measured relative to its hydrostatic value"). · *plain:* Weight is balanced by
     a background pressure; only the leftover pressure p′ moves the fluid.
  4. *did:* Apply the curl-of-a-curl identity · *tex:* $\nabla^2\mathbf u=\nabla(\nabla\cdot\mathbf u)-\nabla\times(\nabla\times\mathbf u)$ ·
     *why:* A vector identity valid for any smooth field (ch04 P122); we want the viscous term expressed through
     divergence and vorticity, the two quantities we control. · *plain:* The Laplacian of velocity splits into a
     compression part and a spin part.
  5. *did:* Insert ∇·u = 0 and ω = ∇×u · *tex:* $\mu\nabla^2\mathbf u=-\mu\nabla\times\boldsymbol\omega$ · *why:* Step 1 kills the
     gradient term; ∇×u is the vorticity by definition. This is (4.40): the net viscous force comes from spatial
     changes of vorticity only. · *plain:* Viscosity pushes on an element only where the spin varies.
  6. *did:* Set the vorticity to zero · *tex:* $\boldsymbol\omega=0\ \Rightarrow\ \mu\nabla^2\mathbf u=0\quad(\mu\neq0)$ · *why:* The flow
     is irrotational by assumption (kept so by Kelvin, R01); the curl of a zero field is zero, whatever the value of μ.
     · *plain:* A viscous fluid that does not spin feels no net viscous force — the stresses on each element cancel.
  7. *did:* Rename p′ as p · *tex:* $\nabla\cdot\mathbf u=0\quad\text{and}\quad\rho\frac{D\mathbf u}{Dt}=-\nabla p$ · *why:* Steps 1, 3 and
     6 leave exactly these two equations; calling the leftover pressure p is the book's convention. This is (6.1). ·
     *plain:* These are the equations of ideal flow.
  8. *did:* Count the boundary conditions left · *tex:* $\mathbf n\cdot\mathbf u=\mathbf n\cdot\mathbf U_s\ \text{kept},\quad \mathbf t\cdot\mathbf u=\mathbf t\cdot\mathbf U_s\ \text{dropped}$ ·
     *why:* (6.1) has only first derivatives of u, so a wall can impose one condition; the second-derivative term that
     needed no-slip has vanished (the book says this only in words). · *plain:* The fluid may slide along walls; it just
     cannot pass through them.
- **Result.** $\nabla\cdot\mathbf u=0$ and $\rho\,D\mathbf u/Dt=-\nabla p$ *(6.1)* with only $\mathbf n\cdot\mathbf u=\mathbf n\cdot\mathbf U_s$ on
  walls — *in words:* irrotational constant-density flow is governed by Euler's equation even in a viscous fluid, at the
  price of letting it slip along walls.
- **Check.** Units: each term of step 2 in N/m³ ✓. Limit: μ = 0 gives the same (6.1) — so an ideal flow is a solution
  of both the viscous and the inviscid equations; only the boundary conditions tell them apart. Numbers: the cylinder
  flow (U = 1 m/s, a = 0.1 m) has μ∇²u ≈ 1e-10 N/m³ by stencils; plane Poiseuille has 1 N/m³ (C01 code).
- **What it means.** Ideal flow is exact for the outer flow; the dropped no-slip condition is repaired by a thin boundary
  layer (Ch. 9). It fails where vorticity exists: wakes, separated regions, rotating or stratified fluids (C01's table).
- **Traps.** Thinking ∇·u = 0 needs steady flow (it needs constant ρ). Forgetting that the p in (6.1) excludes the
  hydrostatic part. Believing "irrotational ⇒ no viscous stress" (the stress survives; its net force vanishes).

### D02 · Vorticity is minus the Laplacian of the stream function: (6.4) → (6.5), (6.6) — ★, 5 steps, in C02 (notebook)
- **Goal.** Express the vorticity of a plane flow through its stream function, and conclude which equation ψ satisfies
  in irrotational flow.
- **Start.** $\omega_z=\frac{\partial v}{\partial x}-\frac{\partial u}{\partial y}$ with $u\equiv\partial\psi/\partial y$, $v\equiv-\partial\psi/\partial x$
  *(6.3)* — *in words:* vorticity is the local spin; ψ gives the velocity.
- **Plan.** (1) Substitute (6.3) into each derivative. (2) Collect the second derivatives. (3) Set ω_z = 0, then allow a
  point of concentrated vorticity.
- **Tools.** ψ (R03) · the plane vorticity (R04, Ch. 3 §3.4) · the Laplacian (Ch. 2 §2.9) · the 2-D delta (C02 primer
  "2-D divergence theorem and the Dirac delta in the plane", placed before D03 — this derivation only names it).
- **Assumptions.** ψ twice continuously differentiable (steps 1–2); irrotational except at isolated points (steps 4–5).
- **Steps.**
  1. *did:* Substitute v = −∂ψ/∂x · *tex:* $\frac{\partial v}{\partial x}=\frac{\partial}{\partial x}\Big(-\frac{\partial\psi}{\partial x}\Big)=-\frac{\partial^2\psi}{\partial x^2}$ ·
     *why:* (6.3) defines v through ψ; differentiating it once more in x is allowed for a smooth ψ. The minus sign comes
     from the definition of v. · *plain:* How fast v changes across x is minus the curvature of ψ in x.
  2. *did:* Substitute u = ∂ψ/∂y · *tex:* $\frac{\partial u}{\partial y}=\frac{\partial}{\partial y}\Big(\frac{\partial\psi}{\partial y}\Big)=\frac{\partial^2\psi}{\partial y^2}$ ·
     *why:* Same move for the other term of ω_z, using u = ∂ψ/∂y from (6.3). · *plain:* How fast u changes up the page is
     the curvature of ψ in y.
  3. *did:* Subtract: ω_z = step 1 − step 2 · *tex:* $\omega_z=-\frac{\partial^2\psi}{\partial x^2}-\frac{\partial^2\psi}{\partial y^2}=-\nabla^2\psi$ ·
     *why:* The definition ω_z = ∂v/∂x − ∂u/∂y; both terms carry a minus, so they combine into minus the Laplacian. This
     is (6.4). · *plain:* The spin of a plane flow is minus the Laplacian of its stream function.
  4. *did:* Set the vorticity to zero · *tex:* $\nabla^2\psi=0$ · *why:* In irrotational flow ω_z = 0 everywhere in the region;
     (6.4) then gives Laplace's equation, (6.5). · *plain:* An irrotational stream function is harmonic.
  5. *did:* Concentrate vorticity Γ at one point · *tex:* $\omega_z=\Gamma\,\delta(x-x')\delta(y-y')\ \Rightarrow\ \nabla^2\psi=-\Gamma\,\delta(x-x')\delta(y-y')$ ·
     *why:* A point vortex has zero vorticity everywhere except at x′, with total ∫ω_z dA = Γ (Stokes: the circulation);
     a delta of weight Γ is exactly that. Substitute into (6.4). This is (6.6). · *plain:* A point vortex is a point
     source of ψ's Laplacian, with strength −Γ.
- **Result.** $\omega_z=-\nabla^2\psi$ *(6.4)*; irrotational flow: $\nabla^2\psi=0$ *(6.5)*; a point vortex: $\nabla^2\psi=-\Gamma\,\delta(x-x')\delta(y-y')$
  *(6.6)* — *in words:* away from vortices ψ is harmonic; each vortex is a point singularity of strength −Γ.
- **Check.** Units: ω_z [1/s], ψ [m²/s], ∇² [1/m²] ✓; Γδδ [m²/s × 1/m²] = 1/s ✓. Example: ψ = −(r²)Γ/(4πa²) inside a
  Rankine core gives −∇²ψ = Γ/πa² (the uniform core vorticity, 31.8 s⁻¹ for Γ = 1, a = 0.1) ✓.
- **What it means.** Ch. 10's vorticity–stream-function method and Ch. 13's potential-vorticity inversion start here:
  knowing ω, solve a Poisson equation for ψ.
- **Traps.** The sign: the minus appears because v = −∂ψ/∂x. "Irrotational" in a region with a vortex means "except at
  the vortex".

### D03 · ln r is harmonic, yet its flux through every circle is 2π: the vortex and source as deltas — ★★, 9 steps, in C02 (notebook)
- **Goal.** Show that ψ = −(Γ/2π) ln r and φ = (m/2π) ln r really solve (6.6) and (6.13) — Laplace everywhere except at the
  origin, with the right strength there. The book leaves this to Exercises 6.1b and 6.3b; we write it out.
- **Start.** $\psi=-\frac{\Gamma}{2\pi}\ln r$ *(6.8)* and $\phi=\frac{m}{2\pi}\ln r$ *(6.15)*, $r=\sqrt{x^2+y^2}$ — *in words:* both
  elements are multiples of ln r.
- **Plan.** (1) Compute ∇² ln r away from the origin (zero). (2) Compute the flux of ∇ ln r through a circle (2π, any
  radius). (3) Use the 2-D divergence theorem to identify a delta. (4) Multiply by the elements' constants.
- **Tools.** polar Laplacian (6.23a) (R08) · 2-D divergence theorem and the Dirac delta in the plane (primer, in C02) ·
  d(ln r)/dr = 1/r (ch01 P36, reminder).
- **Assumptions.** r > 0 for steps 1–2; distributions (deltas) for steps 5–6.
- **Steps.**
  1. *did:* Use the polar Laplacian for f(r) · *tex:* $\nabla^2\ln r=\frac1r\frac{d}{dr}\Big(r\frac{d\ln r}{dr}\Big)$ · *why:* (6.23a) with no
     θ-dependence drops the θ-term; ln r depends on r only. · *plain:* For a function of distance alone, only the radial
     part of the Laplacian remains.
  2. *did:* Differentiate: d(ln r)/dr = 1/r · *tex:* $\frac1r\frac{d}{dr}\Big(r\cdot\frac1r\Big)=\frac1r\frac{d}{dr}(1)=0\quad(r>0)$ · *why:* r × (1/r) = 1
     is constant, so its derivative vanishes; valid wherever r ≠ 0. · *plain:* Away from the origin ln r is harmonic.
  3. *did:* Take the gradient of ln r · *tex:* $\nabla\ln r=\frac{1}{r}\mathbf e_r$ · *why:* The gradient of a function of r alone
     points along e_r with size d/dr (polar gradient); we need it to compute a flux. · *plain:* ln r increases outward,
     by 1/r per metre.
  4. *did:* Integrate its flux over a circle R · *tex:* $\oint_{r=R}\nabla\ln r\cdot\mathbf n\,ds=\int_0^{2\pi}\frac1R\,R\,d\theta=2\pi$ · *why:* On the
     circle n = e_r and ds = R dθ; the 1/R and the R cancel. The book skips this move. · *plain:* The outward flux is 2π
     for every circle — it does not depend on R.
  5. *did:* Apply the 2-D divergence theorem · *tex:* $\int_{r<R}\nabla^2\ln r\,dA=\oint_{r=R}\nabla\ln r\cdot\mathbf n\,ds=2\pi$ · *why:* Gauss in
     the plane turns the area total of a Laplacian into the boundary flux of the gradient (primer); step 4 gives the
     value. · *plain:* The total Laplacian inside any disc round the origin is 2π.
  6. *did:* Recognise a delta of weight 2π · *tex:* $\nabla^2\ln r=2\pi\,\delta(x)\delta(y)$ · *why:* Zero everywhere except at the
     origin (step 2), yet 2π in total over any disc containing it (step 5): that is the definition of 2π times the
     delta. · *plain:* ln r is harmonic except at one point, which carries all of its Laplacian.
  7. *did:* Multiply by −Γ/2π and shift to x′ · *tex:* $\nabla^2\Big(-\frac{\Gamma}{2\pi}\ln\lvert\mathbf x-\mathbf x'\rvert\Big)=-\Gamma\,\delta(x-x')\delta(y-y')$ ·
     *why:* ∇² is linear and does not care where the origin is; the constant −Γ/2π multiplies the delta's weight 2π.
     This is (6.6) for (6.8). · *plain:* The vortex potential has exactly the point vorticity Γ.
  8. *did:* Multiply by m/2π for the source · *tex:* $\nabla^2\Big(\frac{m}{2\pi}\ln\lvert\mathbf x-\mathbf x'\rvert\Big)=m\,\delta(x-x')\delta(y-y')$ · *why:* Same
     move with the source's constant from (6.15); this is (6.13). · *plain:* The source is a point where volume m is
     created per second per metre of depth.
  9. *did:* Read m as a volume flux · *tex:* $\oint\mathbf u\cdot\mathbf n\,ds=\oint\nabla\phi\cdot\mathbf n\,ds=\frac{m}{2\pi}\cdot2\pi=m$ · *why:* u = ∇φ, and
     step 4 gives the flux of ∇ ln r; so the flow out of any circle round the source is m [m²/s]. · *plain:* m is the
     volume leaving the source per second per metre of depth.
- **Result.** $\nabla^2\ln r=2\pi\delta(x)\delta(y)$, so $\psi=-\frac{\Gamma}{2\pi}\ln r$ solves $\nabla^2\psi=-\Gamma\delta\delta$ *(6.6)* and
  $\phi=\frac m{2\pi}\ln r$ solves $\nabla^2\phi=m\delta\delta$ *(6.13)* — *in words:* the elements are harmonic everywhere but at their
  centres, and a loop anywhere round a centre measures its strength.
- **Check.** Units: flux of ∇ψ [m²/s] = Γ ✓. Numbers: Γ = 2π, ∂ψ/∂r = −1/r; at r = 2: −0.5 × 4π = −2π = −Γ ✓ (C02
  tiny example). `delta_flux_check` gives −Γ on radii 0.01…100 m and 0 on a circle that misses the centre.
- **What it means.** The whole information of a vortex or a source sits in one point; that is why superposing a few
  points can build bodies (C03–C05) and why a numerical flux loop can measure circulation (Ch. 14).
- **Traps.** "∇² ln r = 0" is true only for r > 0 — the total is not zero. m is flux per unit depth, not per radian. The
  sign of the vortex term follows from ψ = −(Γ/2π) ln r with Γ counterclockwise.

### D04 · Equipotentials cross streamlines at right angles — ★, 5 steps, in C02 (notebook)
- **Goal.** Show that the lines φ = const and ψ = const are perpendicular wherever the flow moves, and find how ∇φ and ∇ψ
  are related.
- **Start.** $u=\frac{\partial\phi}{\partial x}=\frac{\partial\psi}{\partial y},\ v=\frac{\partial\phi}{\partial y}=-\frac{\partial\psi}{\partial x}$ *(6.10)* and *(6.3)* —
  *in words:* one velocity, two ways to write it.
- **Plan.** (1) Slopes of both families. (2) Their product. (3) The vector form.
- **Tools.** total differential (R03) · level sets (ch02 P75, reminder) · dot product and perpendicular vectors (Ch. 2 §2.2).
- **Assumptions.** u and v not both zero at the point (step 3): stagnation points excluded.
- **Steps.**
  1. *did:* Slope of a streamline · *tex:* $d\psi=-v\,dx+u\,dy=0\ \Rightarrow\ \Big(\frac{dy}{dx}\Big)_\psi=\frac vu$ · *why:* Along a
     streamline ψ does not change; the total differential with (6.3) gives the slope (R03). · *plain:* Streamlines point
     along the velocity.
  2. *did:* Slope of an equipotential · *tex:* $d\phi=u\,dx+v\,dy=0\ \Rightarrow\ \Big(\frac{dy}{dx}\Big)_\phi=-\frac uv$ · *why:* Same move for φ
     with (6.10): along an equipotential φ does not change. · *plain:* Equipotentials run across the velocity.
  3. *did:* Multiply the two slopes · *tex:* $\frac vu\cdot\Big(-\frac uv\Big)=-1$ · *why:* Two lines are perpendicular exactly when
     their slopes multiply to −1; valid when u, v ≠ 0 (otherwise one slope is 0 and the other infinite — still
     perpendicular). · *plain:* The two families meet at right angles.
  4. *did:* Write both gradients as vectors · *tex:* $\nabla\phi=(u,v),\quad\nabla\psi=(-v,u)\ \Rightarrow\ \nabla\phi\cdot\nabla\psi=-uv+vu=0$ · *why:* The
     components come from (6.10) and (6.3); a zero dot product means perpendicular gradients, the vector version of step
     3. · *plain:* The gradients of φ and ψ are perpendicular and equally long, \lvert∇φ\rvert = \lvert∇ψ\rvert = speed.
  5. *did:* Recognise a quarter turn · *tex:* $\nabla\psi=\mathbf e_z\times\nabla\phi,\qquad\nabla\phi=\nabla\psi\times\mathbf e_z$ · *why:* e_z × (u, v, 0) =
     (−v, u, 0) = ∇ψ: turning ∇φ by +90° gives ∇ψ. (The analysis writes ∇φ = e_z × ∇ψ; that has the wrong sign.) · *plain:*
     ∇ψ is the velocity turned a quarter turn counterclockwise.
- **Result.** $\nabla\phi\cdot\nabla\psi=0$, $\lvert\nabla\phi\rvert=\lvert\nabla\psi\rvert=\lvert\mathbf u\rvert$, $\nabla\psi=\mathbf e_z\times\nabla\phi$ — *in words:* the flow net is a
  grid of small squares, streamlines along the flow, equipotentials across it.
- **Check.** Corner flow ψ = 2xy, φ = x² − y² at (1, 1): ∇φ = (2, −2), ∇ψ = (2, 2), dot 4 − 4 = 0 ✓, both lengths
  2.83 ✓, e_z × (2, −2) = (2, 2) ✓. `orthogonality_check` ≲ 1e-10 on the cylinder with circulation.
- **What it means.** Crowded streamlines (fast flow) come with crowded equipotentials; conformal maps (C11) keep this
  square grid square.
- **Traps.** At stagnation points both gradients vanish and "perpendicular" loses meaning. ∇ψ is u turned by +90°, not
  −90°.

### D05 · No flow through a wall ⇔ the wall is a streamline: (6.16) — ★, 5 steps, in C03 (notebook · `superposition_sandbox`)
- **Goal.** Show that the no-through-flow condition on a stationary wall is the same as saying ψ is constant along the
  wall — so any streamline of any flow can be turned into a solid wall.
- **Start.** $\mathbf n\cdot\mathbf U_s=(\mathbf n\cdot\mathbf u)_{\text{on the surface}}$ *(6.16)* — *in words:* the fluid's normal velocity
  matches the wall's.
- **Plan.** (1) Stationary wall. (2) Write n·u with φ. (3) Write it with ψ via the tangent. (4) Conclude.
- **Tools.** directional derivative ∂f/∂n = n·∇f (ch02 P75, reminder) · unit tangent and normal (ch03 P92, reminder) ·
  D04's ∇ψ = e_z × ∇φ.
- **Assumptions.** The wall does not move (step 1); a plane flow (step 3).
- **Steps.**
  1. *did:* Hold the wall still · *tex:* $\mathbf U_s=0\ \Rightarrow\ \mathbf n\cdot\mathbf u=0\ \text{on the wall}$ · *why:* (6.16) with a wall at
     rest; this is the only condition ideal flow keeps (D01 step 8). · *plain:* No fluid crosses the wall.
  2. *did:* Express n·u through φ · *tex:* $\mathbf n\cdot\mathbf u=\mathbf n\cdot\nabla\phi=\frac{\partial\phi}{\partial n}=0$ · *why:* u = ∇φ, and n·∇φ is
     the rate of change of φ along n (directional derivative). · *plain:* φ does not change across the wall.
  3. *did:* Take the tangent t = e_z × n · *tex:* $\frac{\partial\psi}{\partial s}=\mathbf t\cdot\nabla\psi=(\mathbf e_z\times\mathbf n)\cdot(\mathbf e_z\times\nabla\phi)$ ·
     *why:* ∂ψ/∂s is the directional derivative along the wall's tangent; D04 step 5 gives ∇ψ = e_z × ∇φ. The sign of t
     depends on the walking direction. · *plain:* The change of ψ along the wall is a quarter-turned version of the flow
     through it.
  4. *did:* Remove the double quarter turn · *tex:* $(\mathbf e_z\times\mathbf n)\cdot(\mathbf e_z\times\nabla\phi)=\mathbf n\cdot\nabla\phi=\mathbf n\cdot\mathbf u$ · *why:*
     Turning two in-plane vectors by the same angle keeps their dot product (a rotation preserves angles and lengths). ·
     *plain:* How fast ψ changes along the wall equals the flow through the wall.
  5. *did:* Combine steps 2 and 4 · *tex:* $\frac{\partial\phi}{\partial n}=0\iff\frac{\partial\psi}{\partial s}=0\iff\psi=\text{const on the wall}$ · *why:*
     Both derivatives equal n·u; a function whose derivative along a curve vanishes is constant on it. This is the second
     half of (6.16). · *plain:* A stationary wall is a streamline, and any streamline may be replaced by a wall. · *live:*
     "∂ψ/∂s along the black line = 0.000 m/s at your probe".
- **Result.** $\partial\phi/\partial n=0$ or $\partial\psi/\partial s=0$ on a stationary surface *(6.16)* — *in words:* walls and streamlines
  are interchangeable in ideal flow.
- **Check.** Half-body: ψ = m/2 on the dividing streamline and `normal_velocity_on` ≲ 1e-10 there (C03 code); cylinder:
  ψ = 0 on r = a, u_r(a, θ) = 0 ✓.
- **What it means.** The construction method of the chapter: add elements, pick the streamline you like, call it a body.
  A moving body needs n·u = n·U_s instead (D28).
- **Traps.** ∂ψ/∂s = ±n·u — the sign flips with the walking direction, not the physics. Only a *stationary* wall is a
  streamline.

### D06 · The doublet as the limit of a source–sink pair: (6.28) → (6.29) — ★★, 9 steps, in C04 (notebook · `superposition_sandbox`)
- **Goal.** Find the potential left when a source and a sink of equal strength are pushed together while their strength
  grows so that 2mε stays fixed — the doublet — and see how fast the pair approaches it.
- **Start.** $\phi=\frac{m}{2\pi}\ln\sqrt{(x+\varepsilon)^2+y^2}-\frac{m}{2\pi}\ln\sqrt{(x-\varepsilon)^2+y^2}$ *(6.28)* — *in words:* a source +m at
  (−ε, 0) and a sink −m at (+ε, 0).
- **Plan.** (1) Factor r² out of each square. (2) Expand each logarithm for small ε/r. (3) Subtract. (4) Take the limit
  with 2mε = \lvert d\rvert fixed and write it with the dipole vector.
- **Tools.** first-order Taylor ln(1 + s) ≈ s (ch01 P26, reminder) · log rules ln√s = ½ ln s, ln(ab) = ln a + ln b
  (gloss) · a limit with a product held fixed (primer, in C04) · the dipole vector d = Σx_i m_i (N25).
- **Assumptions.** ε ≪ r, i.e. we look at points far from the pair compared with its size (steps 3–4); the order of
  limits: ε → 0 and m → ∞ together (step 7).
- **Steps.**
  1. *did:* Factor r² out of the square · *tex:* $(x\pm\varepsilon)^2+y^2=r^2\Big(1\pm\frac{2\varepsilon x}{r^2}+\frac{\varepsilon^2}{r^2}\Big)$ · *why:* Expand the
     square and use r² = x² + y²; we want a "1 + small" form that a Taylor series can handle. · *plain:* Seen from far
     away, moving a point by ε changes its distance only slightly.
  2. *did:* Split the logarithm · *tex:* $\ln\sqrt{(x\pm\varepsilon)^2+y^2}=\ln r+\tfrac12\ln\Big(1\pm\frac{2\varepsilon x}{r^2}+\frac{\varepsilon^2}{r^2}\Big)$ · *why:*
     ln√(ab) = ½ln a + ½ln b and ½ ln r² = ln r (log rules). The book writes this with the square root kept; we halve
     explicitly. · *plain:* Each element's potential is the centred one plus a small correction.
  3. *did:* Expand ln(1 + s) to first order · *tex:* $\tfrac12\ln(1+s)\approx\tfrac12s,\quad s=\pm\frac{2\varepsilon x}{r^2}+\frac{\varepsilon^2}{r^2}$ · *why:* Taylor
     ln(1 + s) = s − s²/2 + … for small s (P26); ε ≪ r makes s small. · *plain:* A small correction to a log is just the
     small quantity itself.
  4. *did:* Keep the order-ε part · *tex:* $\ln\sqrt{(x\pm\varepsilon)^2+y^2}\cong\ln r\pm\frac{\varepsilon x}{r^2}$ · *why:* The ½ cancels the 2 of
     2εx/r²; terms of order ε² are the same for both signs (step 9 shows they cancel in the difference). This is the
     book's first displayed line. · *plain:* Moving the source by ε changes ln r by ±εx/r².
  5. *did:* Subtract sink from source · *tex:* $\phi\cong\frac{m}{2\pi}\Big[\Big(\ln r+\frac{\varepsilon x}{r^2}\Big)-\Big(\ln r-\frac{\varepsilon x}{r^2}\Big)\Big]=\frac{m\varepsilon}{\pi}\frac{x}{r^2}$ ·
     *why:* Substitute step 4 into (6.28); the ln r terms cancel exactly because the strengths are equal and opposite. ·
     *plain:* What is left is a small term proportional to mε. · *set:* {eps: 0.1, m: 10} · *watch:* "the loops tighten
     toward circles".
  6. *did:* Name the dipole vector · *tex:* $\mathbf d=\sum_i\mathbf x_im_i=(-\varepsilon)(m)\mathbf e_x+(\varepsilon)(-m)\mathbf e_x=-2m\varepsilon\,\mathbf e_x$ · *why:*
     The dipole moment weights each position by its strength (N25); it points from the sink to the source (−x here). ·
     *plain:* The pair's strength and orientation fit into one vector of size 2mε.
  7. *did:* Let ε → 0 with 2mε fixed · *tex:* $\lim_{\varepsilon\to0,\ 2m\varepsilon=\lvert\mathbf d\rvert}\frac{m\varepsilon}{\pi}\frac{x}{r^2}=\frac{\lvert\mathbf d\rvert}{2\pi}\frac{x}{r^2}$ · *why:*
     mε = \lvert d\rvert/2 stays constant by construction (primer: a limit with a product held fixed); ε → 0 alone would
     give 0, m → ∞ alone infinity. · *plain:* Squeezing the pair while strengthening it leaves a finite flow.
  8. *did:* Write it with d and in polar form · *tex:* $\phi=-\frac{\mathbf d\cdot\mathbf x}{2\pi r^2}=\frac{\lvert\mathbf d\rvert}{2\pi}\frac{\cos\theta}{r}$ · *why:* With
     d = −\lvert d\rvert e_x, −d·x = \lvert d\rvert x; and x/r² = cos θ/r. This is (6.29). · *plain:* The doublet's potential
     falls off like 1/r and depends on the direction through cos θ.
  9. *did:* Check what was dropped · *tex:* $\phi_{\text{pair}}-\phi_{\text{doublet}}=O\Big(\frac{\lvert\mathbf d\rvert\,\varepsilon^2}{r^3}\Big)$ · *why:* The
     order-ε² terms of step 3 are even in ±ε and cancel in the difference; the first surviving correction is order mε³,
     i.e. \lvert d\rvert ε². Our addition. · *plain:* The pair differs from the doublet by a relative error of about ε²/r² —
     it vanishes quickly. · *live:* "at your probe the pair differs from the doublet by ε²/3 relative".
- **Result.** $\phi=\frac{m\varepsilon}{\pi}\frac{x}{r^2}=-\frac{\mathbf d\cdot\mathbf x}{2\pi r^2}=\frac{\lvert\mathbf d\rvert}{2\pi}\frac{\cos\theta}{r}$ *(6.29)* — *in words:* a
  merged source–sink pair is a doublet, whose streamlines are circles tangent to the dipole's axis at the centre.
- **Check.** Units: d [m³/s], φ [m³/s × 1/m] = m²/s ✓. Numbers (C04 tiny example): m = 50, ε = 0.02, at (1, 0): pair
  0.318352, doublet 0.318310, relative gap 1.33×10⁻⁴ ≈ ε²/3 ✓. `doublet_limit_error` has observed order 2.00.
- **What it means.** The doublet is the element that disturbs a stream without adding fluid or circulation: stream +
  doublet is a closed body (C06, C13). Its 1/r potential is the far field of every closed 2-D body (C10).
- **Traps.** ln√(1 + s) = ½ ln(1 + s): forgetting the ½ doubles the answer. d points from sink to source, so φ = −d·x/2πr²
  has a minus sign. The limit needs m → ∞ as ε → 0.

### D07 · The half-body: stagnation point, dividing streamline, width — (6.30) → (6.31) — ★★, 10 steps, in C05 (notebook · `superposition_sandbox`)
- **Goal.** For a uniform stream plus a source, find where the fluid stops, which streamline forms the body, how wide the
  body is at each angle, and how wide it becomes far downstream — and check the last result by a mass balance.
- **Start.** $\phi=Ux+\frac{m}{2\pi}\ln\sqrt{x^2+y^2}$ *(6.30)*, $\psi=Uy+\frac{m}{2\pi}\tan^{-1}\frac yx=Ur\sin\theta+\frac{m}{2\pi}\theta$ *(6.31)*,
  θ ∈ [0, 2π) — *in words:* a stream from the left and a source at the origin.
- **Plan.** (1) Velocity. (2) Stagnation point on the axis. (3) ψ there. (4) The body ψ = m/2 and its half-width. (5) The
  far-downstream width, then the mass-balance check.
- **Tools.** superposition (C03) · source velocities (N24) · `np.arctan2` and the angle range (ch02 P70, reminder; gloss in
  C05) · mass conservation (Ch. 4 §4.2).
- **Assumptions.** Steady, plane, ideal (all steps); θ ∈ [0, 2π) so the cut lies inside the body (step 4).
- **Steps.**
  1. *did:* Differentiate φ for the velocity · *tex:* $u=U+\frac{m}{2\pi}\frac{x}{x^2+y^2},\qquad v=\frac{m}{2\pi}\frac{y}{x^2+y^2}$ · *why:* u = ∂φ/∂x,
     v = ∂φ/∂y (6.10), with ∂(ln r)/∂x = x/r² (N24). · *plain:* The stream's velocity plus the source's outward push.
  2. *did:* Look on the axis y = 0 · *tex:* $v(x,0)=0$ · *why:* The source pushes straight out along the axis and the stream
     is horizontal, so only u can vanish there; we search the axis first because symmetry puts S on it. · *plain:* On
     the axis the flow is purely horizontal.
  3. *did:* Solve u = 0 on the axis · *tex:* $U+\frac{m}{2\pi x}=0\ \Rightarrow\ x=-\frac{m}{2\pi U}\equiv-a$ · *why:* At y = 0, x/(x² + y²) = 1/x; a
     negative x (upstream) lets the source's push cancel the stream. · *plain:* The fluid stops a distance a = m/2πU in
     front of the source. · *live:* "x = −m/2πU = −(m)/(2π × U) = xS m".
  4. *did:* Evaluate ψ at the stagnation point · *tex:* $\psi_S=U\cdot0+\frac{m}{2\pi}\cdot\pi=\frac m2$ · *why:* At S, y = 0 and θ = π
     (upstream on the axis); the streamline through S is the dividing streamline, so its ψ value names the body. ·
     *plain:* The body is the streamline ψ = m/2.
  5. *did:* Set ψ = m/2 in (6.31) · *tex:* $Ur\sin\theta+\frac{m}{2\pi}\theta=\frac m2$ · *why:* Every point on the body has the same ψ as
     S (D05). · *plain:* This one equation describes the whole body.
  6. *did:* Solve for the half-width y · *tex:* $h=y=\frac{m(\pi-\theta)}{2\pi U}$ · *why:* Move the θ-term to the right and
     divide by U: m/2 − mθ/2π = m(π − θ)/2π. · *plain:* The body's half-width at polar angle θ grows as θ decreases from
     π.
  7. *did:* Go far downstream, θ → 0 · *tex:* $h_{\max}=\lim_{\theta\to0}\frac{m(\pi-\theta)}{2\pi U}=\frac{m}{2U}$ · *why:* Far downstream on the
     body the angle seen from the source tends to 0 while y stays finite. · *plain:* The body widens to a total width
     m/U. · *watch:* "the width readout approaches m/U".
  8. *did:* Check the nose, θ → π · *tex:* $r=\frac{h}{\sin\theta}=\frac{m}{2\pi U}\frac{\pi-\theta}{\sin\theta}\ \to\ \frac{m}{2\pi U}=a$ · *why:* sin θ ≈
     π − θ near θ = π (small-angle, ch03 P100), so the ratio tends to 1 and the body passes through S. · *plain:* The
     body's nose is the stagnation point, as it must be.
  9. *did:* Balance mass far downstream · *tex:* $2h_{\max}\,U=m$ · *why:* All the source's fluid (m per second per metre)
     stays inside the body (it cannot cross ψ = m/2) and far away moves at U across the width 2h_max. The book leaves
     this to Exercise 6.12. · *plain:* The body is exactly wide enough to carry the source's output away. · *watch:* "the
     teal region's width × U equals m".
  10. *did:* Compare the two routes · *tex:* $h_{\max}=\frac{m}{2U}\quad\text{(steps 7 and 9)}$ · *why:* Two independent arguments —
      geometry of ψ and conservation of mass — must agree; they do. · *plain:* The width far downstream is m/2U either
      way.
- **Result.** $a=\frac{m}{2\pi U}$, body $\psi=\frac m2$, $h(\theta)=\frac{m(\pi-\theta)}{2\pi U}$, $h_{\max}=\frac{m}{2U}$ — *in words:* the source's
  fluid fills a blunt, open body whose nose sits m/2πU upstream and whose width tends to m/U.
- **Check.** Units: m/U [m²/s ÷ m/s] = m ✓. Numbers: U = 1, m = 2π → a = 1 m, h(90°) = 1.571 m, h_max = 3.142 m, mass
  2 × 3.142 × 1 = 6.283 = m ✓. Limit: m → 0 gives a vanishing body.
- **What it means.** The blunt noses of piers, wing leading edges and cliff tops look like this; the body is open because
  the source never stops adding fluid (compare the closed Rankine oval).
- **Traps.** ψ at S is m/2, not m (θ = π, not 2π). Only u = 0 needs solving on the axis. With the principal
  `arctan2` range the lower half of the body shows ψ = −m/2 — choose θ ∈ [0, 2π) so the cut lies inside the body.

### D08 · The half-body's surface pressure and where it crosses zero — ★★, 8 steps, in C05 (notebook)
- **Goal.** Find the pressure coefficient along the half-body's surface as a function of the polar angle, and the angle
  where the surface pressure equals the free-stream pressure. The book only plots this (Fig. 6.8); we derive it.
- **Start.** $C_p=1-\frac{\lvert\mathbf u\rvert^2}{U^2}$ *(6.32)* and the body $r=\frac{m(\pi-\theta)}{2\pi U\sin\theta}$ (D07 step 6 with
  y = r sin θ) — *in words:* the pressure follows the speed, and we know where the body is.
- **Plan.** (1) Write the velocity in polar pieces. (2) Put the point on the body. (3) Compute \lvert u\rvert²/U². (4)
  Read C_p, its limits, and its zero.
- **Tools.** C_p (N28) · D07's body shape · `brentq` (ch03 P108, reminder).
- **Assumptions.** Points on the body only (step 3 onward); θ in radians inside formulas.
- **Steps.**
  1. *did:* Define the shorthand k(θ) · *tex:* $k\equiv\frac{\sin\theta}{\pi-\theta}$ · *why:* The body shape contains this ratio; naming it
     keeps the algebra short. · *plain:* k measures how close a body point is to the source, relative to the stream.
  2. *did:* Evaluate the source's push on the body · *tex:* $\frac{m}{2\pi r}=\frac{m}{2\pi}\cdot\frac{2\pi U\sin\theta}{m(\pi-\theta)}=Uk$ · *why:* Substitute
     the body's r from D07; the m and 2π cancel. · *plain:* On the body, the source's speed m/2πr is U times k.
  3. *did:* Write u and v in polar pieces · *tex:* $u=U+\frac{m}{2\pi r}\cos\theta,\qquad v=\frac{m}{2\pi r}\sin\theta$ · *why:* D07 step 1 with
     x/r² = cos θ/r and y/r² = sin θ/r. · *plain:* The source adds a radial push of size m/2πr.
  4. *did:* Insert step 2 · *tex:* $u=U(1+k\cos\theta),\qquad v=Uk\sin\theta$ · *why:* Replace m/2πr by Uk — allowed only on the
     body. · *plain:* On the body both velocity components are U times simple functions of θ.
  5. *did:* Square and add · *tex:* $\frac{\lvert\mathbf u\rvert^2}{U^2}=(1+k\cos\theta)^2+k^2\sin^2\theta=1+2k\cos\theta+k^2$ · *why:* Expand the square and use
     cos² + sin² = 1. · *plain:* The surface speed squared, relative to the stream.
  6. *did:* Insert into (6.32) · *tex:* $C_p=-(2k\cos\theta+k^2)$ · *why:* C_p = 1 − \lvert u\rvert²/U²; the 1s cancel. · *plain:* The
     surface pressure coefficient as a function of θ alone.
  7. *did:* Check the nose, θ → π · *tex:* $k\to1,\ \cos\theta\to-1\ \Rightarrow\ C_p\to-(-2+1)=1$ · *why:* sin θ ≈ π − θ near π, so k → 1;
     C_p = 1 is the stagnation value — a consistency check. · *plain:* Full stagnation pressure at the nose, as expected.
  8. *did:* Find the zero · *tex:* $C_p=0\iff2\cos\theta+k=0\iff\tan\theta=-2(\pi-\theta)$ · *why:* Factor C_p = −k(2cos θ + k), k ≠ 0
     on the body; divide by cos θ(π − θ). Solved by `brentq` on (π/2, π) where tan changes sign. · *plain:* The surface
     pressure equals p∞ at θ ≈ 113.2° from +x.
- **Result.** $C_p=-(2k\cos\theta+k^2)$, $k=\frac{\sin\theta}{\pi-\theta}$, with $C_p=0$ where $\tan\theta=-2(\pi-\theta)$ (θ ≈ 113.2°) — *in
  words:* over-pressure at the nose, zero at 113°, suction on the shoulder (minimum −0.587 near 63°), back to zero far
  downstream.
- **Check.** Units: C_p dimensionless ✓. Numbers: θ = 90° → k = 0.6366, C_p = −0.405 ✓ (C05 code prints −0.4053). Limit
  θ → 0: k → 0, C_p → 0 — far downstream the body is parallel to the stream at speed U.
- **What it means.** Because C_p → 0 downstream, the infinitely long body still feels no net pressure force: the suction
  on the shoulders balances the nose's push (Exercise 6.13; `half_body_net_force` → 0).
- **Traps.** Using the field formula off the body (step 4 holds only on it). Degrees vs radians in tan θ = −2(π − θ).

### D09 · The circular cylinder: why d = −2πUa² e_x, and C_p = 1 − 4 sin²θ (6.33)–(6.35) — ★, 6 steps, in C06 (notebook · `cylinder_circulation_lift`)
- **Goal.** Choose the doublet that turns a uniform stream into flow round a circle of radius a, and find the velocity and
  pressure on its surface.
- **Start.** $\psi=Ur\sin\theta-\frac{D}{2\pi}\frac{\sin\theta}{r}$ — stream (6.7) plus a doublet of strength D pointing upstream,
  $\mathbf d=-D\mathbf e_x$ (from Im of (6.49)) — *in words:* a stream and an unknown doublet.
- **Plan.** (1) Demand ψ = 0 on r = a. (2) Read off D. (3) Velocities by (6.21)–(6.22). (4) Surface speed and C_p.
- **Tools.** superposition (C03) · doublet (C04) · polar velocities (N18, N19) · C_p (N28).
- **Assumptions.** Steady, plane, no circulation (added in D10).
- **Steps.**
  1. *did:* Require ψ = 0 on the circle · *tex:* $\psi(a,\theta)=\Big(Ua-\frac{D}{2\pi a}\Big)\sin\theta=0\ \ \text{for all }\theta$ · *why:* D05: the
     circle is a body if it is a streamline; the value 0 is the one the x-axis already has, so the axis and circle join.
     · *plain:* The circle must be one streamline, for every angle.
  2. *did:* Solve for the doublet strength · *tex:* $D=2\pi Ua^2\ \Rightarrow\ \mathbf d=-2\pi Ua^2\,\mathbf e_x$ · *why:* The bracket must vanish,
     since sin θ does not; the doublet points upstream, opposing the stream. · *plain:* The right doublet has strength
     2πUa². · *live:* "D = 2πUa² = 2π × U × a² = d m³/s".
  3. *did:* Substitute back · *tex:* $\psi=U\Big(r-\frac{a^2}{r}\Big)\sin\theta,\qquad\phi=U\Big(r+\frac{a^2}{r}\Big)\cos\theta$ · *why:* Insert D into the
     start; φ follows the same way from (6.29) with the same d. This is (6.33). · *plain:* The cylinder's stream function
     and potential.
  4. *did:* Differentiate for the velocity · *tex:* $u_r=\frac1r\frac{\partial\psi}{\partial\theta}=U\Big(1-\frac{a^2}{r^2}\Big)\cos\theta,\quad u_\theta=-\frac{\partial\psi}{\partial r}=-U\Big(1+\frac{a^2}{r^2}\Big)\sin\theta$ ·
     *why:* The polar forms (6.21)–(6.22). This is (6.34). · *plain:* Far away the stream; near the body the flow bends
     round it.
  5. *did:* Evaluate on the surface r = a · *tex:* $u_r(a,\theta)=0,\qquad u_\theta(a,\theta)=-2U\sin\theta$ · *why:* a²/r² = 1 on the
     surface: u_r = 0 confirms no flow through the body (6.16). · *plain:* On the surface the fluid slides round at
     2U\lvert sin θ\rvert — twice the stream at the shoulders.
  6. *did:* Insert the surface speed in (6.32) · *tex:* $C_p(a,\theta)=1-\frac{4U^2\sin^2\theta}{U^2}=1-4\sin^2\theta$ · *why:* C_p = 1 − \lvert u\rvert²/U²
     with \lvert u\rvert = 2U\lvert sin θ\rvert. This is (6.35). · *plain:* +1 at the front and back, −3 at the shoulders. ·
     *watch:* "C_p = −3 at θ = ±90°".
- **Result.** $\psi=U(r-a^2/r)\sin\theta$ *(6.33)*, $u_\theta(a,\theta)=-2U\sin\theta$, $C_p=1-4\sin^2\theta$ *(6.35)* — *in words:* a doublet
  2πUa² pointing upstream turns the stream into flow round a cylinder; the surface pressure is symmetric front–back and
  top–bottom.
- **Check.** Units ✓ (Ua² m³/s). Numbers: air, U = 10 m/s → ½ρU² = 60 Pa, +60 Pa at θ = 0, π, −180 Pa at ±90° ✓. C_p is
  unchanged by θ → −θ and θ → π − θ: the symmetry behind d'Alembert.
- **What it means.** Both symmetries mean every push has a mirror push: zero drag and zero lift (C06), until
  circulation breaks the top–bottom symmetry (C07).
- **Traps.** d points upstream ("opposing the stream"). The surface speed is 2U at the shoulders, not U.

### D10 · The cylinder with circulation: from (6.36) to the stagnation points (6.38) and the free point — ★★, 9 steps, in C07 (notebook · `cylinder_circulation_lift`)
- **Goal.** Add a clockwise point vortex to the cylinder flow, keep the circle a streamline, and find where the fluid stops
  — on the surface for weak circulation, off it for strong circulation.
- **Start.** $\psi=U\big(r-\frac{a^2}{r}\big)\sin\theta$ *(6.33)* and a vortex at the origin whose flow circulation is −Γ (the book's
  clockwise Γ) — *in words:* the cylinder plus a clockwise whirl.
- **Plan.** (1) Write the clockwise vortex's ψ. (2) Shift it so ψ(a) = 0 → (6.36). (3) Surface speed (6.37). (4) Solve
  u_θ = 0 on the surface (6.38). (5) When that fails, solve on the negative y-axis.
- **Tools.** vortex ψ (6.8) with the sign convention (N08, ⚠️ callout) · polar velocity (6.22) (N19) · two solutions of
  sin θ = s (gloss) · quadratic formula and Vieta's product (ch02 P71, reminder).
- **Assumptions.** Γ ≥ 0 in the book's clockwise sense (the negative case mirrors it); steady flow.
- **Steps.**
  1. *did:* Write the clockwise vortex's ψ · *tex:* $\psi_v=-\frac{(-\Gamma)}{2\pi}\ln r=\frac{\Gamma}{2\pi}\ln r$ · *why:* (6.8) is written for a
     counterclockwise strength; a clockwise vortex of size Γ has counterclockwise strength −Γ. The book's footnote chooses
     this sign so (6.40) comes out +ρUΓ. · *plain:* A clockwise whirl has ψ = +(Γ/2π) ln r.
  2. *did:* Shift ψ by a constant · *tex:* $\psi=U\Big(r-\frac{a^2}{r}\Big)\sin\theta+\frac{\Gamma}{2\pi}\ln\Big(\frac ra\Big)$ · *why:* Subtracting
     (Γ/2π) ln a changes no velocity (only derivatives of ψ matter) and makes ψ = 0 on r = a again. This is (6.36). ·
     *plain:* The circle is still the streamline ψ = 0.
  3. *did:* Differentiate for u_θ · *tex:* $u_\theta=-\frac{\partial\psi}{\partial r}=-U\Big(1+\frac{a^2}{r^2}\Big)\sin\theta-\frac{\Gamma}{2\pi r}$ · *why:*
     (6.22); d/dr of ln(r/a) is 1/r. · *plain:* The swirl adds a clockwise speed Γ/2πr everywhere.
  4. *did:* Evaluate on the surface · *tex:* $u_\theta(a,\theta)=-2U\sin\theta-\frac{\Gamma}{2\pi a}$ · *why:* r = a in step 3; u_r = 0 there
     still (the vortex adds no radial velocity). This is (6.37). · *plain:* On top (θ = 90°) the speeds add; below they
     subtract.
  5. *did:* Set the surface speed to zero · *tex:* $\sin\theta=-\frac{\Gamma}{4\pi aU}$ · *why:* A stagnation point on the surface
     needs u_θ = 0; solve step 4 for sin θ. This is (6.38). · *plain:* The stagnation points sit where sin θ equals
     −Γ/4πaU — below the axis. · *live:* "sin θ = −gam".
  6. *did:* Take both roots for Γ < 4πaU · *tex:* $\theta_1=-\arcsin\frac{\Gamma}{4\pi aU},\qquad\theta_2=-\pi+\arcsin\frac{\Gamma}{4\pi aU}$ · *why:*
     sin θ = s has two solutions in one turn, θ and π − θ (gloss); both lie in the lower half since s < 0. · *plain:*
     Two stagnation points, mirror images about the vertical axis, sliding down as Γ grows.
  7. *did:* Let them meet at Γ = 4πaU · *tex:* $\sin\theta=-1\ \Rightarrow\ \theta=-\frac\pi2$ · *why:* The right-hand side of (6.38)
     reaches −1; the two roots coincide at the bottom. · *plain:* One stagnation point, at the bottom of the cylinder.
     · *set:* {gam: 1} · *watch:* "the two markers merge at −90°".
  8. *did:* For Γ > 4πaU search the −y axis · *tex:* $u_\theta\Big(r,-\frac\pi2\Big)=U\Big(1+\frac{a^2}{r^2}\Big)-\frac{\Gamma}{2\pi r}=0\ \Rightarrow\ Ur^2-\frac{\Gamma}{2\pi}r+Ua^2=0$ ·
     *why:* (6.38) has no solution when \lvert sin θ\rvert would exceed 1; on θ = −π/2, cos θ = 0 makes u_r = 0 for free, so
     only u_θ = 0 is needed; multiply by r²/1. · *plain:* The stagnation point leaves the body along the downward axis.
  9. *did:* Solve the quadratic, keep r > a · *tex:* $r=\frac{1}{4\pi U}\Big[\Gamma\pm\sqrt{\Gamma^2-(4\pi aU)^2}\Big],\qquad r_+r_-=a^2$ · *why:*
     Quadratic formula; Vieta's product of roots is Ua²/U = a², so one root is outside the circle and one inside — only
     r₊ > a is in the fluid. · *plain:* A free stagnation point below the body at r₊. · *set:* {gam: 1.5} · *live:* "r₊ =
     rf m, r₋ = rm m, r₊r₋ = a²".
- **Result.** $u_\theta(a,\theta)=-2U\sin\theta-\Gamma/2\pi a$ *(6.37)*, stagnation at $\sin\theta=-\Gamma/4\pi aU$ *(6.38)* for Γ ≤ 4πaU, and at
  $r=\frac{\Gamma+\sqrt{\Gamma^2-(4\pi aU)^2}}{4\pi U}$ on the negative y-axis for Γ > 4πaU — *in words:* circulation slides the
  stagnation points down, merges them, then lifts one off the body.
- **Check.** Units: Γ/(aU) dimensionless ✓. Numbers: U = 10, a = 0.1, Γ = 2 → sin θ = −0.159, θ = −9.16° and −170.84° ✓;
  Γ = 6πaU = 18.85 → r₊ = 0.2618 m, r₋ = 0.0382 m, product 0.0100 = a² ✓. Limit Γ = 0: θ = 0 and −π (the front and back
  points of C06).
- **What it means.** Every Γ gives a valid flow round the same cylinder (N38); something else — viscosity at a sharp
  edge, the Kutta condition — must choose Γ (Ch. 14).
- **Traps.** A clockwise vortex has ψ = +(Γ/2π) ln r (the sign flip of (6.8)). Both surface roots are below the axis for
  Γ > 0. Of the two r-roots keep the one with r > a.

### D11 · The lift on the cylinder: from the surface pressure (6.39) to L = ρUΓ (6.40), and D = 0 — ★★, 9 steps, in C07 (notebook · `cylinder_circulation_lift`)
- **Goal.** Integrate the surface pressure of the cylinder with circulation to find the sideways force per unit length,
  and show that the drag is zero.
- **Start.** $p(r=a,\theta)=p_\infty+\tfrac12\rho\Big[U^2-\Big(-2U\sin\theta-\frac{\Gamma}{2\pi a}\Big)^2\Big]$ *(6.39)* — *in words:* Bernoulli
  (6.18) with the surface speed (6.37).
- **Plan.** (1) Expand the square. (2) Write the force as −∮p n dl. (3) Integrate term by term over a full turn. (4) Only
  the cross term survives.
- **Tools.** Bernoulli (R05) · net pressure force −∮p n dA (ch01 P28, reminder) · integrals of sines and cosines over a
  full period (primer, in C06).
- **Assumptions.** Steady ideal flow; forces per unit length of cylinder; p∞ uniform.
- **Steps.**
  1. *did:* Expand the square · *tex:* $\Big(2U\sin\theta+\frac{\Gamma}{2\pi a}\Big)^2=4U^2\sin^2\theta+\frac{2U\Gamma}{\pi a}\sin\theta+\frac{\Gamma^2}{4\pi^2a^2}$ ·
     *why:* (A + B)² = A² + 2AB + B²; the sign inside does not matter once squared. The book skips this. · *plain:* Three
     pieces: the stream's, a mixed piece, and the vortex's.
  2. *did:* Write the pressure term by term · *tex:* $p=p_\infty+\tfrac12\rho U^2-2\rho U^2\sin^2\theta-\frac{\rho U\Gamma}{\pi a}\sin\theta-\frac{\rho\Gamma^2}{8\pi^2a^2}$ ·
     *why:* Multiply step 1 by −½ρ and add to p∞ + ½ρU². · *plain:* Four kinds of term: constants, a sin² term, a sin θ
     term and a constant from the vortex alone.
  3. *did:* Write the force on the body · *tex:* $\mathbf F=-\oint p\,\mathbf n\,dl,\qquad\mathbf n=\mathbf e_r=(\cos\theta,\sin\theta),\ dl=a\,d\theta$ · *why:*
     Pressure pushes along the inward normal −n (P28); on the circle the outward normal is e_r and an arc element is
     a dθ. · *plain:* Add the pushes all round the surface.
  4. *did:* Take the upward component · *tex:* $L=-\int_0^{2\pi}p(a,\theta)\sin\theta\,a\,d\theta$ · *why:* n·e_y = sin θ; the minus sign
     because pressure on the top surface (sin θ > 0) pushes down (Fig. 6.13). · *plain:* Lift is the upward sum of the
     pushes.
  5. *did:* Drop the constant terms · *tex:* $\int_0^{2\pi}\sin\theta\,d\theta=0$ · *why:* Over a full period sin θ integrates to zero
     (primer); p∞, ½ρU² and −ρΓ²/8π²a² are constants, so they give no force. · *plain:* A uniform pressure pushes equally
     from all sides.
  6. *did:* Drop the sin² term · *tex:* $\int_0^{2\pi}\sin^2\theta\sin\theta\,d\theta=\int_0^{2\pi}\sin^3\theta\,d\theta=0$ · *why:* sin³θ is odd about θ = π
     (it changes sign between the upper and lower halves), so it integrates to zero (primer). · *plain:* The stream's own
     pressure is the same on top and bottom — no lift from it.
  7. *did:* Keep the cross term · *tex:* $L=-\int_0^{2\pi}\Big(-\frac{\rho U\Gamma}{\pi a}\sin\theta\Big)\sin\theta\,a\,d\theta=\frac{\rho U\Gamma}{\pi}\int_0^{2\pi}\sin^2\theta\,d\theta$ ·
     *why:* The only term with a single sin θ survives against the sin θ of the lift; the a's cancel. · *plain:* The lift
     comes entirely from the mixing of stream and swirl. · *live:* "ρUΓ/π × π = L N/m".
  8. *did:* Integrate sin² over a turn · *tex:* $\int_0^{2\pi}\sin^2\theta\,d\theta=\pi\ \Rightarrow\ L=\rho U\Gamma$ · *why:* The average of sin²
     over a period is ½ (primer); ½ × 2π = π. This is (6.40). · *plain:* Lift per metre = density × stream speed ×
     circulation.
  9. *did:* Repeat for the drag · *tex:* $D=-\int_0^{2\pi}p\cos\theta\,a\,d\theta=0$ · *why:* Every term of step 2 times cos θ — cos θ,
     sin²θ cos θ, sin θ cos θ — integrates to zero over a full period (primer): the pressure is symmetric front to back. ·
     *plain:* Circulation gives lift but still no drag. · *watch:* "the drag bar stays at 0".
- **Result.** $L=\rho U\Gamma$ *(6.40)*, $D=0$ — *in words:* the lift per unit length is density × speed × circulation, whatever
  the radius; the drag stays zero.
- **Check.** Units: ρUΓ [kg/m³ × m/s × m²/s] = N/m ✓. Numbers: ρ = 1.2, U = 10, Γ = 2 → 24 N/m; the 64-point trapezoid
  gives 24.000… (C07 code) ✓. Limit Γ = 0: L = 0 (C06).
- **What it means.** The radius dropped out — the first hint that the lift does not depend on the shape (C10 proves it
  for any body). The direction: clockwise circulation in a stream from the left lifts up.
- **Traps.** Forgetting the minus sign of step 4 (n·e_y = sin θ but the push is −p n). Not expanding the square fully:
  only the cross term 2UΓ sin θ/πa survives. Using a counterclockwise Γ in (6.40) flips the sign of L.

### D12 · Images: vortex images flip sign, source images keep it — ★, 5 steps, in C08 (notebook · `vortex_wall_images`)
- **Goal.** Turn a flow in the whole plane into a flow above a wall y = 0 by adding mirror images, and see why a vortex's
  image has the opposite sign while a source's has the same sign.
- **Start.** An unbounded flow with stream function ψ₁(x, y) (vorticity ω₁) or potential φ₁(x, y) (sources q₁), all its
  singularities in y > 0 — *in words:* a flow that ignores the wall.
- **Plan.** (1) Build an odd ψ₂. (2) Check the wall. (3) Check the fluid's field equation. (4) Build an even φ₂. (5) Check
  both.
- **Tools.** odd and even functions under a mirror (gloss) · superposition (C03) · a wall is a streamline (D05) · (6.4),
  (6.11).
- **Assumptions.** A straight wall y = 0; the fluid is y > 0; the images lie in y < 0 (outside the fluid).
- **Steps.**
  1. *did:* Subtract the mirror image of ψ₁ · *tex:* $\psi_2(x,y)=\psi_1(x,y)-\psi_1(x,-y)$ · *why:* We want a ψ that vanishes on
     y = 0; a function minus its mirror is odd in y, and an odd function is zero on the mirror line. · *plain:* The flow
     plus its upside-down, reversed copy.
  2. *did:* Evaluate on the wall · *tex:* $\psi_2(x,0)=\psi_1(x,0)-\psi_1(x,0)=0$ · *why:* Direct substitution y = 0; ψ₂ is constant
     along the wall, so the wall is a streamline (D05). · *plain:* No fluid crosses the wall.
  3. *did:* Check the field equation in the fluid · *tex:* $\nabla^2\psi_2=-\omega_1(x,y)+\omega_1(x,-y)$ · *why:* ∇² does not change under
     y → −y, so each part obeys (6.4); the image's vorticity lies in y < 0, outside the fluid. A vortex Γ at (x₀, y₀) gets
     −Γ at (x₀, −y₀). · *plain:* In the fluid nothing changed; behind the wall sits a vortex of opposite spin.
  4. *did:* Add the mirror image of φ₁ · *tex:* $\phi_2(x,y)=\phi_1(x,y)+\phi_1(x,-y)$ · *why:* For sources we control φ; an even
     function has zero slope across its mirror line, so v = ∂φ₂/∂y = 0 on y = 0. · *plain:* The flow plus its upside-down
     copy, same sign.
  5. *did:* Check both conditions · *tex:* $\frac{\partial\phi_2}{\partial y}\Big\rvert_{y=0}=0,\qquad\nabla^2\phi_2=q_1(x,y)+q_1(x,-y)$ · *why:* Differentiating
     an even function gives an odd one, zero on y = 0; the image source (same m) lies outside the fluid. · *plain:* A
     source's image has the same sign; the wall again has no flow through it.
- **Result.** Vortex images: $\psi_2=\psi_1(x,y)-\psi_1(x,-y)$; source images: $\phi_2=\phi_1(x,y)+\phi_1(x,-y)$ — *in words:* mirror the
  flow in the wall; vortices come back with opposite spin, sources with the same strength.
- **Check.** A vortex at (0, 1) and its image: v on the wall = 0 at 1001 points (C08 code, ≤ 1e-14) ✓; with the same-sign
  image (the wrong variant) the wall leaks.
- **What it means.** A wall is a symmetry line. The same trick gives the bottom condition of water waves (Ch. 7), coasts
  in ocean models (Ch. 13) and ground effect (Ch. 14).
- **Traps.** Two different conditions (ψ₂ = 0 vs ∂φ₂/∂y = 0) give the same wall. Mirroring a counterclockwise vortex
  makes it clockwise — hence the sign.

### D13 · Example 6.1: a vortex beside a wall and the pressure it makes on the wall — ★★, 11 steps, in C08 (notebook · `vortex_wall_images`)
- **Goal.** For a clockwise vortex a distance h from a wall, find how it moves and the pressure at the wall point
  opposite its starting position, as a function of time.
- **Start.** $\psi=\frac{\Gamma}{2\pi}\Big[-\ln\sqrt{(x+\xi_x)^2+(y-\xi_y)^2}+\ln\sqrt{(x-\xi_x)^2+(y-\xi_y)^2}\Big]$ — the vortex (strength −Γ,
  i.e. clockwise Γ) at ξ = (ξ_x, ξ_y) and its image (+Γ) at (−ξ_x, ξ_y), wall x = 0 — *in words:* D12's rule with the wall
  vertical.
- **Plan.** (1) Check the wall. (2) Move the vortex with its image's velocity. (3) Build φ. (4) Unsteady Bernoulli at the
  origin: ∂φ/∂t and v there. (5) Combine into one fraction.
- **Tools.** wall image (R16), drift (R17) · unsteady Bernoulli (4.83) with g = 0 (Ch. 4 §4.9, R05) · derivative of
  tan⁻¹(Y/X) with a moving Y (gloss) · chain rule (ch01 P49, reminder).
- **Assumptions.** The vortex's own (self-induced) velocity is zero — a real core is symmetric (step 2); fluid at rest at
  p∞ far away; no gravity; ∂φ/∂t at a fixed point (step 7).
- **Steps.**
  1. *did:* Check u = 0 on the wall · *tex:* $u(0,y,t)=\frac{\partial\psi}{\partial y}\Big\rvert_{x=0}=0$ · *why:* At x = 0 the two logarithms have the
     same argument and opposite signs, so ψ and its y-derivative vanish on the wall (D12). · *plain:* No flow through the
     wall at any time.
  2. *did:* Move the vortex with its image · *tex:* $\frac{d\xi_x}{dt}=0,\qquad\frac{d\xi_y}{dt}=\frac{\Gamma}{2\pi}\frac{1}{2\xi_x}$ · *why:* A free vortex
     moves with the local fluid; its own velocity is taken as zero, the image (+Γ at distance 2ξ_x) induces Γ/2π(2ξ_x)
     upward. · *plain:* The eddy slides along the wall, never toward it.
  3. *did:* Integrate from (h, 0) · *tex:* $\boldsymbol\xi(t)=\Big(h,\ \frac{\Gamma t}{4\pi h}\Big)$ · *why:* ξ_x stays h; then dξ_y/dt = Γ/4πh is
     constant, integrated with ξ_y(0) = 0. · *plain:* It climbs the wall at Γ/4πh (0.0796 m/s for Γ = 1 m²/s, h = 1 m). ·
     *live:* "ξ_y = Γt/4πh = xiy m".
  4. *did:* Write φ from the two vortices · *tex:* $\phi=\frac{\Gamma}{2\pi}\tan^{-1}\Big(\frac{y-\Gamma t/4\pi h}{x+h}\Big)-\frac{\Gamma}{2\pi}\tan^{-1}\Big(\frac{y-\Gamma t/4\pi h}{x-h}\Big)$ ·
     *why:* A counterclockwise vortex has φ = (Γ/2π) × angle (6.47); the image (+Γ) and the vortex (−Γ) each contribute
     one. Near the origin neither arctan jumps. · *plain:* The potential of the pair at time t.
  5. *did:* Write unsteady Bernoulli · *tex:* $\frac{\partial\phi}{\partial t}+\frac12\lvert\nabla\phi\rvert^2+\frac p\rho=\frac{p_\infty}{\rho}$ · *why:* (4.83) with g = 0;
     the constant is fixed far away, where the fluid is still at p∞. · *plain:* Pressure is lowered by speed and by a
     growing potential.
  6. *did:* Evaluate at the wall origin · *tex:* $\frac{p(0,0,t)-p_\infty}{\rho}=-\Big(\frac{\partial\phi}{\partial t}+\tfrac12(u^2+v^2)\Big)_{x=y=0}$ · *why:*
     Rearrange step 5 at the sensor point; u = 0 there by step 1. · *plain:* The sensor reads two contributions.
  7. *did:* Differentiate φ at a fixed point · *tex:* $\frac{\partial}{\partial t}\tan^{-1}\frac{Y}{X}=\frac{X\dot Y}{X^2+Y^2},\quad Y=y-\frac{\Gamma t}{4\pi h},\ \dot Y=-\frac{\Gamma}{4\pi h}$ ·
     *why:* Chain rule with d(tan⁻¹s)/ds = 1/(1 + s²) (gloss); X = x ± h does not change, only Y moves with the vortex.
     The book skips this move. · *plain:* The potential changes because the vortices move past.
  8. *did:* Put x = y = 0 · *tex:* $\Big(\frac{\partial\phi}{\partial t}\Big)_0=\frac{\Gamma}{2\pi}\Big(-\frac{\Gamma}{4\pi h}\Big)\Big[\frac{h}{h^2+s^2}+\frac{h}{h^2+s^2}\Big]=-\frac{\Gamma^2}{4\pi^2}\frac{1}{h^2+s^2}$ ·
     *why:* With s = Γt/4πh: X = h for the first term, X = −h for the second (whose minus sign turns the second bracket
     positive), Y = −s in both. · *plain:* The potential at the sensor keeps falling while the eddy passes. · *live:*
     "∂φ/∂t at the origin = dphidt m²/s²".
  9. *did:* Find v at the origin · *tex:* $v(0,0,t)=-\frac{\partial\psi}{\partial x}\Big\rvert_0=\frac{\Gamma h}{\pi}\frac{1}{h^2+s^2}$ · *why:* Differentiate ψ in x
     (both terms give h/(h² + s²) with the same sign at the origin) and use v = −∂ψ/∂x (6.3). · *plain:* The fluid at the
     sensor slides along the wall, fastest when the eddy is level with it.
  10. *did:* Substitute into step 6 · *tex:* $\frac{p-p_\infty}{\rho}=\frac{\Gamma^2}{4\pi^2}\frac{1}{h^2+s^2}-\frac{\Gamma^2h^2}{2\pi^2}\frac{1}{(h^2+s^2)^2}$ · *why:* −∂φ/∂t
      from step 8 and −½v² from step 9 (v² = Γ²h²/π²(h² + s²)², halved). · *plain:* A positive unsteady part minus a
      speed part.
  11. *did:* Put over one denominator · *tex:* $\frac{p(0,0,t)-p_\infty}{\rho}=\frac{\Gamma^2}{4\pi^2}\frac{s^2-h^2}{(s^2+h^2)^2},\quad s=\frac{\Gamma t}{4\pi h}$ · *why:*
      Multiply the first term by (h² + s²)/(h² + s²) and combine: (h² + s²) − 2h² = s² − h². This is the book's result.
      · *plain:* Suction while the eddy is closer than h along the wall, over-pressure when farther. · *live:* "p − p∞ = p
      Pa" · *watch:* "the dot on the trace".
- **Result.** $\frac{p(0,0,t)-p_\infty}{\rho}=\frac{\Gamma^2}{4\pi^2}\frac{(\Gamma t/4\pi h)^2-h^2}{((\Gamma t/4\pi h)^2+h^2)^2}$ — *in words:* the sensor feels
  −ρΓ²/4π²h² when the eddy passes, zero at t = 4πh²/Γ and a maximum push ρΓ²/32π²h² at t = 4√3πh²/Γ (our addition: set
  d/ds of (s² − h²)/(s² + h²)² to zero, s = √3h).
- **Check.** Units: Γ²/h² [m²/s²] ✓ (p/ρ). Numbers: water, Γ = 1, h = 1: −25.33 Pa at t = 0 (= +25.33 unsteady − 50.66
  speed), 0 at 12.57 s, +3.17 Pa at 21.77 s ✓ (C08 code); the numerical route (moving vortices, ∂φ/∂t by central
  differences) agrees to 1e-7. Limit t → ∞: p → p∞.
- **What it means.** A wall-mounted sensor (or an ocean-bottom recorder) sees a passing eddy as a dip flanked by small
  humps; the humps come from the flow changing in time, not from its speed.
- **Traps.** ∂φ/∂t is at a fixed point while the vortex moves. The self-induced velocity is set to zero by a modelling
  argument, not by the equations. The arctan branch must not jump between the vortex and its image (it does not at the
  origin).

### D14 · Cauchy–Riemann (6.44) and the complex velocity dw/dz = u − iv (6.45) — ★★, 8 steps, in C09 (notebook · `complex_potential_corners`)
- **Goal.** Show that if w(z) = φ + iψ has a derivative independent of direction, then φ and ψ satisfy the Cauchy–Riemann
  equations, the derivative is u − iv, and both φ and ψ are harmonic — so every analytic function is a flow.
- **Start.** $w(z)=\phi(x,y)+i\psi(x,y)$ *(6.42)*, $\frac{dw}{dz}=\lim_{\delta z\to0}\frac{w(z+\delta z)-w(z)}{\delta z}$ the same for every direction of
  δz — *in words:* w is analytic.
- **Plan.** (1) Take the limit along x. (2) Take it along iy. (3) Equate. (4) Read off u − iv. (5) Differentiate once more
  for Laplace and orthogonality.
- **Tools.** complex derivative and analytic functions (primer, in C09) · 1/i = −i · Schwarz's theorem (ch04 P121,
  reminder) · ψ and φ (R03, N09).
- **Assumptions.** φ, ψ twice continuously differentiable (steps 6–8).
- **Steps.**
  1. *did:* Step along x: δz = δx · *tex:* $\frac{dw}{dz}=\frac{\partial\phi}{\partial x}+i\frac{\partial\psi}{\partial x}$ · *why:* With δz real, the quotient is
     the ordinary partial derivative of φ + iψ in x. · *plain:* The derivative seen walking east.
  2. *did:* Step along iy: δz = iδy · *tex:* $\frac{dw}{dz}=\frac{1}{i}\Big(\frac{\partial\phi}{\partial y}+i\frac{\partial\psi}{\partial y}\Big)$ · *why:* Now the step is
     i times δy, so we divide the change in w by iδy. · *plain:* The derivative seen walking north.
  3. *did:* Use 1/i = −i · *tex:* $\frac{dw}{dz}=\frac{\partial\psi}{\partial y}-i\frac{\partial\phi}{\partial y}$ · *why:* 1/i = −i since i × (−i) = 1; multiply
     out. The book skips this move. · *plain:* Walking north, the roles of φ and ψ swap.
  4. *did:* Equate real and imaginary parts · *tex:* $\frac{\partial\phi}{\partial x}=\frac{\partial\psi}{\partial y},\qquad\frac{\partial\phi}{\partial y}=-\frac{\partial\psi}{\partial x}$ ·
     *why:* Analytic means steps 1 and 3 are the same complex number; two complex numbers are equal when both parts are.
     This is (6.44). · *plain:* The Cauchy–Riemann equations: φ and ψ are locked together.
  5. *did:* Read the velocity from step 1 · *tex:* $\frac{dw}{dz}=u-iv$ · *why:* ∂φ/∂x = u by (6.10) and ∂ψ/∂x = −v by (6.3). This is
     (6.45). · *plain:* The complex derivative is the velocity with v flipped.
  6. *did:* Differentiate the first CR equation in x · *tex:* $\frac{\partial^2\phi}{\partial x^2}=\frac{\partial^2\psi}{\partial x\partial y}$ · *why:* Differentiate
     both sides of the first equation of (6.44) once more; allowed for smooth φ, ψ. · *plain:* φ's curvature in x
     equals a mixed derivative of ψ.
  7. *did:* Add the second, differentiated in y · *tex:* $\frac{\partial^2\phi}{\partial x^2}+\frac{\partial^2\phi}{\partial y^2}=\frac{\partial^2\psi}{\partial x\partial y}-\frac{\partial^2\psi}{\partial y\partial x}=0$ ·
     *why:* The second equation of (6.44) differentiated in y gives φ_yy = −ψ_xy; mixed partials commute (Schwarz), so
     they cancel. The same move with the roles swapped gives ∇²ψ = 0. · *plain:* Both φ and ψ are harmonic —
     irrotational and incompressible automatically.
  8. *did:* Check orthogonality with CR · *tex:* $\nabla\phi\cdot\nabla\psi=\phi_x\psi_x+\phi_y\psi_y=\phi_x(-\phi_y)+\phi_y\phi_x=0$ · *why:* Replace ψ_x by
     −φ_y and ψ_y by φ_x from (6.44). · *plain:* Equipotentials and streamlines of any analytic w cross at right angles
     (D04 again).
- **Result.** $\frac{\partial\phi}{\partial x}=\frac{\partial\psi}{\partial y},\ \frac{\partial\phi}{\partial y}=-\frac{\partial\psi}{\partial x}$ *(6.44)*, $\frac{dw}{dz}=u-iv$ *(6.45)*,
  $\nabla^2\phi=\nabla^2\psi=0$ — *in words:* every analytic function of z is a plane ideal flow, and its derivative is the
  mirrored velocity.
- **Check.** w = z² at z = 1 + i: dw/dz = 2 + 2i ⇒ u = 2, v = −2 m/s; ψ = 2xy gives u = 2x = 2, v = −2y = −2 ✓. Control:
  f = z* gives 1 along x and −1 along iy — no derivative, not a flow ✓.
- **What it means.** One complex function replaces two real ones; complex analysis (residues, maps) now works for flows
  (C10, C11). ⚠️ The book writes "(6.5) and (6.12), respectively" for φ and ψ — the other way round.
- **Traps.** Along iy the quotient divides by i. The velocity arrow is conj(dw/dz), not dw/dz.

### D15 · Flow in a corner, w = Azⁿ (6.46): walls and the speed at the tip — ★★, 7 steps, in C09 (notebook · `complex_potential_corners`)
- **Goal.** Show that w = Azⁿ is the flow in a corner of angle α = π/n, and find how the speed behaves at the corner's tip.
- **Start.** $w(z)=Az^n$, A real, $z=re^{i\theta}$ *(6.43)* — *in words:* a power of z.
- **Plan.** (1) Polar form of zⁿ. (2) Read ψ. (3) Find the walls. (4) Recover n = 2. (5) Speed near the tip. (6) The
  plate.
- **Tools.** complex logarithm, powers and branch cuts (primer, in C09) · polar form (N41) · D14's dw/dz = u − iv.
- **Assumptions.** 0 ≤ θ ≤ π/n inside the corner; the branch cut of zⁿ placed outside it (step 7).
- **Steps.**
  1. *did:* Raise the polar form to the n-th power · *tex:* $z^n=r^ne^{in\theta}$ · *why:* zⁿ = e^{n ln z} and ln z = ln r + iθ (primer),
     valid on the branch where θ runs across the corner. · *plain:* The distance is raised to the power n, the angle
     multiplied by n.
  2. *did:* Split into φ and ψ · *tex:* $w=Ar^n\cos n\theta+i\,Ar^n\sin n\theta\ \Rightarrow\ \psi=Ar^n\sin n\theta$ · *why:* Euler's formula (ch01
     P45) and w = φ + iψ (6.42). This is (6.46). · *plain:* The stream function of the corner flow.
  3. *did:* Find where ψ = 0 · *tex:* $\sin n\theta=0\ \Rightarrow\ \theta=0\ \text{and}\ \theta=\frac{\pi}{n}\equiv\alpha$ · *why:* Walls are streamlines
     (D05); ψ = 0 on these two rays, which bound the region between them. · *plain:* A corner of angle α = π/n.
  4. *did:* Check n = 2 · *tex:* $w=Az^2=A(x^2-y^2)+i\,2Axy$ · *why:* Expand (x + iy)²; Re is (6.27), Im is (6.24), walls at 0 and
     90°. · *plain:* The right-angle stagnation flow is the case n = 2.
  5. *did:* Differentiate · *tex:* $\frac{dw}{dz}=nAz^{n-1}=\frac{A\pi}{\alpha}z^{(\pi-\alpha)/\alpha}$ · *why:* Power rule for an analytic
     function; n = π/α rewrites the exponent n − 1 as (π − α)/α. · *plain:* The complex velocity of the corner flow.
  6. *did:* Read the speed at the tip · *tex:* $\lvert\mathbf u\rvert=nAr^{n-1}\ \to\ \begin{cases}0,&\alpha<\pi\\ \infty,&\alpha>\pi\end{cases}\quad(r\to0)$ · *why:* \lvert
     dw/dz\rvert = \lvert u − iv\rvert = speed; r^{n−1} vanishes for n > 1 and blows up for n < 1. · *plain:* Narrow corners
     are calm (a stagnation point); corners wider than a straight wall are violent. · *live:* "\lvert dw/dz\rvert ∝
     r^(−1/3)" (n = 2/3).
  7. *did:* Take n = ½: the flat plate · *tex:* $w=Az^{1/2},\quad\alpha=2\pi$ · *why:* The walls θ = 0 and θ = 2π are both the
     positive x-axis: a plate; z^{1/2} needs its cut along the plate itself so the fluid round it is smooth. · *plain:*
     Flow round the edge of a thin plate, with infinite speed at the edge.
- **Result.** $w=Az^n$ *(6.46)* is flow in a corner of angle α = π/n with $dw/dz=(A\pi/\alpha)z^{(\pi-\alpha)/\alpha}$ — *in words:* the
  corner's angle alone decides whether the tip is a stagnation point (α < π) or a point of infinite speed (α > π).
- **Check.** n = 2 gives (6.24)/(6.27) ✓; log–log slopes of \lvert dw/dz\rvert vs r equal n − 1 for n = 2, 1, ⅔, ½ (C09
  figure) ✓. Units: A [m^{2−n}/s] so that w is m²/s ✓.
- **What it means.** Example 6.2's re-entrant 270° step (C12) has α = 3π/2, speed ∝ r^{−1/3} — why grids converge slowly
  there; Ch. 9's wedge flows U ∝ x^m grow from this outer flow; a sharp trailing edge (α = 2π) needs the Kutta condition
  (Ch. 14).
- **Traps.** ψ = Arⁿ sin nθ vanishes on θ = π/n, not on θ = π. n < 1 means infinite speed at the tip. The branch cut must
  lie outside the fluid.

### D16 · The force on a body from its surface pressure: (6.54) → (6.55) → (6.56) — ★★, 8 steps, in C10 (notebook · `blasius_kutta_contour`)
- **Goal.** Write the drag and lift per unit depth on a stationary 2-D body as contour integrals of the pressure round its
  cross-section — the starting point of Blasius's theorem.
- **Start.** $\int_{A^*}\rho\mathbf u(\mathbf u\cdot\mathbf n)dA=-\int_{A^*}p\mathbf n\,dA+\mathbf F$ *(6.54)* — *in words:* steady momentum balance for a
  stationary control surface A* (R18).
- **Plan.** (1) Put A* on the body. (2) The flux term vanishes. (3) Newton's third law gives the force on the body.
  (4) Parametrise the contour and its outward normal. (5) Split into components.
- **Tools.** CV momentum (R18) · momentum flux (ch04 P114, reminder) · outward normal of a counterclockwise contour
  (gloss in C10) · net pressure force −∮p n dA (ch01 P28, reminder).
- **Assumptions.** Steady flow; stationary body (step 2); span B large, forces per unit depth (step 4).
- **Steps.**
  1. *did:* Let A* hug the body · *tex:* $\mathbf u\cdot\mathbf n=0\ \text{on the body}$ · *why:* A stationary body is a streamline
     (D05), so no fluid crosses a control surface lying on it. · *plain:* Nothing flows through the body's skin.
  2. *did:* Drop the momentum-flux term · *tex:* $\int_{A^*}\rho\mathbf u(\mathbf u\cdot\mathbf n)dA=0$ · *why:* The integrand contains u·n = 0; with
     no flux, only pressure transmits force between fluid and body. · *plain:* Across the body's skin the fluid pushes
     only by pressure.
  3. *did:* Write the force on the body · *tex:* $B(D\mathbf e_x+L\mathbf e_y)=-\int_{A^*}p\,\mathbf n\,dA$ · *why:* Pressure pushes on the
     body along −n, n the body's outward normal (P28); by Newton's third law this is minus the force the body applies to
     the fluid, F = −B(De_x + Le_y). · *plain:* Drag and lift are the sum of the pressure pushes on the skin.
  4. *did:* Divide by the span B · *tex:* $D\mathbf e_x+L\mathbf e_y=-\frac1B\int_{A^*}p\,\mathbf n\,dA,\qquad dA=B\,ds$ · *why:* The body does not
     vary along its span, so each strip of length B contributes equally; this is (6.55). · *plain:* Forces per metre of
     span come from a line integral round the cross-section.
  5. *did:* Walk the contour counterclockwise · *tex:* $d\mathbf s=\mathbf e_xdx+\mathbf e_ydy,\qquad\mathbf n=\frac{\mathbf e_xdy-\mathbf e_ydx}{ds}$ · *why:*
     The tangent (dx, dy)/ds turned clockwise by 90° is (dy, −dx)/ds, which points out of the body when the body is on
     the left (at its rightmost point dy > 0, so n = +e_x ✓). · *plain:* The outward normal from the contour's steps. ·
     *watch:* "the outward normal arrows".
  6. *did:* Substitute n and dA · *tex:* $D\mathbf e_x+L\mathbf e_y=-\oint_Cp\,(\mathbf e_xdy-\mathbf e_ydx)$ · *why:* B ds × (e_x dy − e_y dx)/ds / B:
     the B's and ds's cancel. · *plain:* Only the steps dx and dy along the contour remain.
  7. *did:* Split into components · *tex:* $D=-\oint_Cp\,dy,\qquad L=\oint_Cp\,dx$ · *why:* Match the e_x and e_y parts; the minus of
     −e_y dx turns positive. This is (6.56). · *plain:* Drag collects pressure against vertical steps, lift against
     horizontal ones. · *live:* "D = Dp, L = Lp N/m".
  8. *did:* Check on a circle · *tex:* $x=a\cos\theta\ \Rightarrow\ dx=-a\sin\theta\,d\theta,\quad L=-\int_0^{2\pi}p\sin\theta\,a\,d\theta$ · *why:* The
     counterclockwise circle; the lift formula becomes D11 step 4 exactly — a consistency check. · *plain:* The general
     formula contains the cylinder calculation.
- **Result.** $D\mathbf e_x+L\mathbf e_y=-\frac1B\int_{A^*}p\mathbf n\,dA$ *(6.55)*, $D=-\oint_Cp\,dy$, $L=\oint_Cp\,dx$ *(6.56)* — *in words:* on a
  counterclockwise contour round the body, drag and lift per unit depth are the pressure integrated against dy and dx.
- **Check.** Uniform pressure: ∮dy = ∮dx = 0 on a closed curve ⇒ no force ✓. A clockwise contour flips both signs (the
  code rejects it). Cylinder: 24 N/m by (6.56) (C10 code) ✓.
- **What it means.** (6.56) needs the pressure on the body; the next step (D17) rewrites it so the body surface is no
  longer needed.
- **Traps.** n is outward only for counterclockwise traversal. F in (6.54) is the force on the fluid; D and L are on the
  body. The flux term vanishes only because u·n = 0 on a stationary body.

### D17 · Blasius's theorem: $D-iL=\frac{i\rho}{2}\oint_C(dw/dz)^2dz$ (6.60) — ★★★, 12 steps, in C10 (notebook · `blasius_kutta_contour`)
- **Goal.** Turn the pressure integrals (6.56) into one complex contour integral of (dw/dz)², and show the contour can then
  be moved off the body.
- **Start.** $D=-\oint_Cp\,dy,\quad L=\oint_Cp\,dx$ *(6.56)* — *in words:* drag and lift from the surface pressure.
- **Plan.** (1) Combine D and L into one complex number and use dz* = dx − i dy. (2) Insert Bernoulli with |q|² = qq*.
  (3) Drop the constant. (4) Use that the velocity is tangent to the body. (5) Use Cauchy's theorem to move the contour.
- **Tools.** complex conjugate (ch02 P81, reminder) · Bernoulli (R05) · complex velocity (6.45) (C09) · Cauchy's integral
  theorem (primer, in C10) · the complex plane in numpy (primer, in C09).
- **Assumptions.** Steady plane ideal flow; the body is a streamline (step 7); (dw/dz)² analytic between the body and the
  new contour (step 11).
- **Steps.**
  1. *did:* Combine into D − iL · *tex:* $D-iL=-\oint_Cp\,dy-i\oint_Cp\,dx$ · *why:* One complex number carries both components; the
     minus on L is chosen so the result is compact. · *plain:* Drag in the real part, minus lift in the imaginary part.
  2. *did:* Factor out −i · *tex:* $D-iL=-i\oint_Cp\,(dx-i\,dy)$ · *why:* −i(dx − i dy) = −i dx + i² dy = −i dx − dy, which
     reproduces step 1 term by term. · *plain:* The two integrals become one.
  3. *did:* Recognise dz* · *tex:* $D-iL=-i\oint_Cp\,dz^*$ · *why:* dz = dx + i dy, so its conjugate is dx − i dy. This is (6.57).
     · *plain:* The force is the pressure integrated against the conjugate step.
  4. *did:* Insert Bernoulli with the far-field constant · *tex:* $p=p_\infty+\tfrac12\rho U^2-\tfrac12\rho(u^2+v^2)$ · *why:* (6.18) with the
     constant evaluated far upstream (R05). · *plain:* Pressure is high where the fluid is slow.
  5. *did:* Write the speed squared as a product · *tex:* $u^2+v^2=(u-iv)(u+iv)$ · *why:* \lvert q\rvert² = q q* for any complex q = u −
     iv (P81); this lets dw/dz appear. · *plain:* Speed squared is the complex velocity times its conjugate.
  6. *did:* Integrate the constant to zero · *tex:* $\oint_C\big(p_\infty+\tfrac12\rho U^2\big)dz^*=\big(p_\infty+\tfrac12\rho U^2\big)\Big(\oint_Cdz\Big)^*=0$ ·
     *why:* The integral of dz round a closed curve is the end point minus the start, zero; the constant comes out. The
     book says this without the reason. · *plain:* A uniform pressure gives no force.
  7. *did:* Keep the velocity term, fix the sign · *tex:* $D-iL=-i\oint_C\Big(-\tfrac12\rho\Big)(u-iv)(u+iv)\,dz^*=\frac{i\rho}{2}\oint_C(u-iv)(u+iv)\,dz^*$ · *why:*
     Steps 3–6 combined (this is (6.58) with the constant removed); −i × (−½ρ) = +iρ/2. · *plain:* Only the velocity term
     is left.
  8. *did:* Use that the velocity is tangent · *tex:* $u+iv=\lvert q\rvert e^{i\vartheta},\quad dz=\lvert dz\rvert e^{i\vartheta}$ · *why:* On the body the
     velocity points along the surface (D05), so the velocity and the step dz have the same direction angle ϑ. ·
     *plain:* On the body, flow and surface point the same way. · *watch:* "on the body, u + iv and dz point the same
     way".
  9. *did:* Swap the conjugate · *tex:* $(u+iv)dz^*=\lvert q\rvert\lvert dz\rvert=(u-iv)dz$ · *why:* The phases e^{iϑ} and e^{−iϑ} cancel either
     way round; this is (6.59). · *plain:* On the body the conjugate can move from the step to the velocity.
  10. *did:* Use dw/dz = u − iv · *tex:* $D-iL=\frac{i\rho}{2}\oint_C\Big(\frac{dw}{dz}\Big)^2dz$ · *why:* Substitute step 9 into step 7 and
      (6.45). This is (6.60), Blasius's theorem. · *plain:* The force is a contour integral of the complex velocity
      squared.
  11. *did:* Apply Cauchy between two contours · *tex:* $\oint_C\Big(\frac{dw}{dz}\Big)^2dz-\oint_{C'}\Big(\frac{dw}{dz}\Big)^2dz=0$ · *why:* (dw/dz)² is
      analytic in the fluid between the body C and a larger curve C′ (no singularity there); Cauchy's theorem (primer)
      makes the integral round the ring zero. · *plain:* The same number on any curve round the body. · *set:* {R: 0.3}
      · *watch:* "the readout does not move".
  12. *did:* Note what was needed · *tex:* $D-iL=\frac{i\rho}{2}\oint_{C'}\Big(\frac{dw}{dz}\Big)^2dz\quad\text{for any }C'\text{ enclosing the body}$ · *why:*
      Step 9 held only on the body, but the final integrand is analytic, so step 11 carries it anywhere outside —
      provided C′ encloses the body and no other singularity. · *plain:* We may compute the force far away. · *set:* {R:
      0.07} · *watch:* "cutting the body breaks it".
- **Result.** $D-iL=\frac{i\rho}{2}\oint_C\Big(\frac{dw}{dz}\Big)^2dz$ *(6.60)*, on any contour enclosing the body — *in words:* the force on
  a 2-D body is a single contour integral of the squared complex velocity, and the contour may be moved freely through
  the fluid.
- **Check.** Units: ρ (dw/dz)² dz [kg/m³ × m²/s² × m] = N/m ✓. Cylinder with Γ_cw = 2, U = 10, ρ = 1.2: (0, 24) N/m on R =
  a, 2a, 10a (C10 code) ✓; the trapezoid on a circle converges exponentially.
- **sympy check intent (★★★, `check_src`):** re-run steps 7–10 on the cylinder with circulation. `import sympy as sp` ·
  `th, U, a, G, rho = sp.symbols('theta U a Gamma rho', positive=True)` (symbols) · `z = a*sp.exp(sp.I*th)` (a point on
  the body) · `dwdz = U*(1 - a**2/z**2) + sp.I*G/(2*sp.pi*z)` ((6.52) differentiated) · `u_minus_iv = dwdz` and
  `u_plus_iv = sp.conjugate(dwdz)` (the velocity and its conjugate) · `dz = sp.diff(z, th)` (the step along the circle,
  times dθ) · `print(sp.simplify(sp.expand_complex(u_plus_iv*sp.conjugate(dz) - u_minus_iv*dz)))` (step 9 on the body →
  0) · `print(sp.integrate(sp.conjugate(dz), (th, 0, 2*sp.pi)))` (step 6: ∮dz* → 0) · `I = sp.integrate(sp.expand(dwdz**2
  * dz), (th, 0, 2*sp.pi))` (Blasius's integral on the body) · `print(sp.simplify(sp.I*rho/2*I))` (→ −I*Gamma*U*rho:
  D = 0, L = ρUΓ) · and the same integral with a → 2a in z (a bigger contour, same dw/dz) gives the same value (step 11).
  Every line commented in the builder.
- **What it means.** The body's shape enters only through dw/dz; on a large circle every body looks alike, which D18
  exploits. Blasius's theorem is the complex-variable form of a momentum balance (Exercise 6.27's control-volume proof
  gives the same number, `cv_force_on_body`).
- **Traps.** −i(dx − i dy) = −i dx − dy (a sign slip here flips D). −i × (−½ρ) = +iρ/2. Step 9 holds only on the body; the
  contour may leave the body only because the final integrand is analytic.

### D18 · Kutta–Zhukhovsky for any body: D = 0, L = ρUΓ (6.61)–(6.62) — ★★★, 11 steps, in C10 (notebook · `blasius_kutta_contour`)
- **Goal.** Evaluate Blasius's integral for an arbitrary 2-D body with circulation, using only what the flow looks like far
  away.
- **Start.** $D-iL=\frac{i\rho}{2}\oint_C\Big(\frac{dw}{dz}\Big)^2dz$ *(6.60)* on a large circle — *in words:* D17's result with the contour
  moved far out.
- **Plan.** (1) Write dw/dz outside the body as a Laurent series. (2) Identify its first coefficients physically. (3)
  Square it. (4) Integrate term by term with the residue rule.
- **Tools.** Laurent series and residues (primer, in C10) · residue theorem on a circle (N63) · sympy for complex series
  and residues (primer, in C10) · the element potentials (6.47)–(6.49) (N45–N47).
- **Assumptions.** Closed body (step 4); all singularities inside the body; clockwise Γ (book sign).
- **Steps.**
  1. *did:* Expand dw/dz outside the body · *tex:* $\frac{dw}{dz}=c_0+\frac{c_{-1}}{z}+\frac{c_{-2}}{z^2}+\dots$ · *why:* dw/dz is analytic outside a
     circle containing the body, so it has a Laurent series there (primer); no positive powers, because the velocity
     stays finite far away. The book asserts this form. · *plain:* Far away, the velocity is a sum of decaying powers.
  2. *did:* Match the uniform stream · *tex:* $c_0=U$ · *why:* As \lvert z\rvert → ∞, dw/dz → u − iv = U (a stream along x). ·
     *plain:* The leading term is the free stream.
  3. *did:* Integrate to w and name the terms · *tex:* $w=Uz+\frac{m}{2\pi}\ln z+\frac{i\Gamma}{2\pi}\ln z+\frac{d}{2\pi z}+\dots,\quad c_{-1}=\frac{m+i\Gamma}{2\pi},\ c_{-2}=-\frac{d}{2\pi}$ ·
     *why:* The antiderivative of c₋₁/z is c₋₁ ln z; its real part is a source (6.48), its imaginary part a clockwise
     vortex ((6.47) with −Γ); c₋₂/z² integrates to a doublet (6.49). · *plain:* From far away any body is a stream, a
     source, a vortex and a doublet.
  4. *did:* Close the body: m = 0 · *tex:* $c_{-1}=\frac{i\Gamma}{2\pi}$ · *why:* A closed body emits no net fluid: the flux of u out of
     a large circle is m, and it must vanish. · *plain:* Only the circulation survives in the 1/z term.
  5. *did:* Write the far-field velocity · *tex:* $\frac{dw}{dz}=U+\frac{i\Gamma}{2\pi z}-\frac{d}{2\pi z^2}+\dots$ · *why:* Substitute steps 2–4 into the series of step 1; the doublet term −d/2πz² is c₋₂/z². ·
     *plain:* Stream + circulation + doublet, in powers of 1/z.
  6. *did:* Square the series · *tex:* $\Big(\frac{dw}{dz}\Big)^2=U^2+\frac{iU\Gamma}{\pi}\frac1z-\Big(\frac{Ud}{\pi}+\frac{\Gamma^2}{4\pi^2}\Big)\frac1{z^2}+\dots$ · *why:*
     Multiply out: 2 × U × iΓ/2πz, and (iΓ/2πz)² + 2U(−d/2πz²). ⚠️ The book prints the 1/z² coefficient as (Ud/π −
     Γ²/4π²) inside an extra square; this is the correct one (the sympy cell prints both). · *plain:* The square has a
     1/z term made of U times Γ. · *live:* "iUΓ/π = res".
  7. *did:* Integrate powers round a circle · *tex:* $\oint z^{-n}dz=iR^{1-n}\int_0^{2\pi}e^{i(1-n)\theta}d\theta=2\pi i\,\delta_{n1}$ · *why:* z = Re^{iθ}, dz =
     iz dθ (N63); the exponential integrates to zero unless its exponent is zero. · *plain:* Only the 1/z term gives
     anything.
  8. *did:* Keep the residue · *tex:* $\oint_C\Big(\frac{dw}{dz}\Big)^2dz=2\pi i\cdot\frac{iU\Gamma}{\pi}$ · *why:* By step 7 the constant, 1/z², … terms
     vanish; the coefficient of 1/z is the residue. · *plain:* The integral is 2πi times the residue.
  9. *did:* Multiply out · *tex:* $2\pi i\cdot\frac{iU\Gamma}{\pi}=-2U\Gamma$ · *why:* i × i = −1 and the π cancels. · *plain:* A real
     number: −2UΓ.
  10. *did:* Insert into Blasius · *tex:* $D-iL=\frac{i\rho}{2}(-2U\Gamma)=-i\rho U\Gamma$ · *why:* (6.60) with step 9. This is the first form
      of (6.62). · *plain:* The force is purely imaginary. · *live:* "D − iL = −i × LKJ".
  11. *did:* Read off drag and lift · *tex:* $D=0,\qquad L=\rho U\Gamma$ · *why:* The real part is D, the imaginary part is −L;
      compare both sides. This is (6.62). · *plain:* No drag, and lift ρUΓ, for any shape.
- **Result.** $D-iL=\frac{i\rho}{2}2\pi i\Big(\frac{iU\Gamma}{\pi}\Big)=-i\rho U\Gamma$, so $D=0$ and $L=\rho U\Gamma$ *(6.62)* — *in words:* the
  Kutta–Zhukhovsky theorem: every 2-D body with circulation Γ in a stream U feels lift ρUΓ and no drag.
- **Check.** Units ✓. The body's doublet d dropped out — the shape does not matter. Numbers: Γ = 2, U = 10, ρ = 1.2 → 24
  N/m for the cylinder and the Zhukhovsky ellipse (C10 code) ✓. Open body (m ≠ 0): the 1/z coefficient becomes U(m +
  iΓ)/π and D = −ρUm (a thrust) — the held-source force.
- **sympy check intent (★★★, `check_src`):** `z, U, G, d = sp.symbols('z U Gamma d')` · `f = U + sp.I*G/(2*sp.pi*z) -
  d/(2*sp.pi*z**2)` (step 5) · `sq = sp.expand(f**2)` (step 6) · `print(sq.coeff(z, -1))` (→ I*Gamma*U/pi) ·
  `true_c2 = sq.coeff(z, -2)`; `printed_c2 = U*d/sp.pi - G**2/(4*sp.pi**2)` (the book's coefficient) ·
  `print(sp.simplify(true_c2 - (-(U*d/sp.pi + G**2/(4*sp.pi**2)))), sp.simplify(true_c2 - printed_c2))` (→ 0 and a
  non-zero difference: the printed coefficient is wrong) · `res = sp.residue(sq, z, 0)` (step 8's residue) · `rho =
  sp.symbols('rho', positive=True)`; `F = sp.I*rho/2*2*sp.pi*sp.I*res` (steps 9–10) · `print(sp.re(F), -sp.im(F))` (→ 0 and
  Gamma*U*rho, with U, Γ, ρ declared real) · and a numeric cross-check: `ch06.kutta_zhukhovsky_sym()` returns the same
  dict. Every line commented.
- **What it means.** Lift depends only on U and Γ; the shape decides which Γ nature chooses (Ch. 14's Kutta condition).
  The 2-D d'Alembert paradox holds for every body.
- **Traps.** The ln z term's real coefficient is a source, the imaginary one a vortex; m = 0 only for a closed body. The
  book's 1/z² coefficient and outer square are misprints — only 1/z matters. The clockwise Γ sign gives L = +ρUΓ.

### D19 · Conformal maps keep angles: (6.63) → (6.64) — ★★, 6 steps, in C11 (notebook · `conformal_joukowski`)
- **Goal.** Show that an analytic map turns and stretches every small figure at a point by the same amount, so angles
  between small line elements are preserved — except where the derivative is zero or infinite.
- **Start.** $\delta w=\frac{dw}{dz}\delta z$ *(6.63)* — *in words:* near a point, a small step in z maps to a small step in w by
  multiplication with the derivative.
- **Plan.** (1) Write the derivative in polar form. (2) See what it does to one element. (3) And to a second. (4)
  Subtract the angles. (5) Find where it fails.
- **Tools.** the complex plane in numpy: multiplying adds angles (primer, in C09) · polar form (N41).
- **Assumptions.** w = f(z) analytic at z₀ with f′(z₀) ≠ 0, ∞ (steps 1–4); small elements (first order in δz).
- **Steps.**
  1. *did:* Write f′ in polar form · *tex:* $\frac{dw}{dz}=\Big\lvert\frac{dw}{dz}\Big\rvert e^{i\chi}$ · *why:* Any non-zero complex number has a modulus
     and an argument (P153); χ = arg f′(z₀). · *plain:* The derivative is a stretch and a turn.
  2. *did:* Apply it to one element · *tex:* $\delta w=\Big\lvert\frac{dw}{dz}\Big\rvert\lvert\delta z\rvert e^{i(\chi+\arg\delta z)}$ · *why:* Multiplying complex
     numbers multiplies moduli and adds arguments. · *plain:* The small step is stretched by \lvert f′\rvert and turned by χ.
     · *live:* "\lvert f′\rvert = mag, arg f′ = arg°".
  3. *did:* Apply it to a second element · *tex:* $\delta'w=\frac{dw}{dz}\delta'z$ · *why:* The same point, the same derivative — (6.64);
     a second direction δ′z from z₀. · *plain:* Every small step at z₀ gets the same stretch and turn.
  4. *did:* Subtract the two angles · *tex:* $\arg\delta'w-\arg\delta w=(\chi+\arg\delta'z)-(\chi+\arg\delta z)=\arg\delta'z-\arg\delta z$ · *why:* χ appears in
     both and cancels. So α = β: the angle between the images equals the angle between the originals. · *plain:* Small
     crosses stay crosses.
  5. *did:* Look where f′ = 0 · *tex:* $w=z^2:\ \arg w=2\arg z\ \text{near }0$ · *why:* At a critical point χ is undefined and the
     first-order argument fails; for z² the next term doubles every angle at the origin. · *plain:* At a critical point
     angles are multiplied, not kept. · *set:* {map: 'z2', elementAt: 0} · *watch:* "the image angle doubles".
  6. *did:* Apply to a flow net · *tex:* $\phi=\text{const}\ \perp\ \psi=\text{const}\ \text{in }w\ \Rightarrow\ \text{also in }z$ · *why:* A uniform grid of φ-
     and ψ-lines crosses at right angles in the w-plane; by step 4 their images still do (N66). · *plain:* A map carries
     a flow net to a flow net.
- **Result.** $\delta'w/\delta w=\delta'z/\delta z$ in angle, i.e. α = β *(6.64)*, wherever $0<\lvert dw/dz\rvert<\infty$ — *in words:* an analytic
  map is locally a rotation plus a scaling, so it keeps angles between small elements.
- **Check.** f = z² at z₀ = 1 + 0.5i: two perpendicular steps stay perpendicular (90° → 90°); at z₀ = 0: 45° → 90° (C11
  code) ✓. Large figures distort (the ellipse is not a scaled circle).
- **What it means.** Streamlines and equipotentials, perpendicular in any flow, stay perpendicular under an analytic map:
  a solved flow becomes a new solved flow (D20, D21).
- **Traps.** The statement is local: only small figures keep their shape. At f′ = 0 or ∞ (ζ = ±b for Zhukhovsky) angles
  are not kept — those are the plate's ends.

### D20 · The Zhukhovsky map sends circles to a slit and to ellipses (6.65)–(6.67) — ★★, 8 steps, in C11 (notebook · `conformal_joukowski`)
- **Goal.** Find the images of circles centred at the origin under z = ζ + b²/ζ: the circle of radius b and a larger one of
  radius a.
- **Start.** $z=\zeta+\frac{b^2}{\zeta}$ *(6.65)* — *in words:* the Zhukhovsky transformation.
- **Plan.** (1) Far away. (2) The circle of radius b. (3) A circle of radius a. (4) Separate x and y, eliminate θ. (5)
  Foci, critical points, two-to-one.
- **Tools.** Euler's formula (ch01 P45) · e^{iθ} + e^{−iθ} = 2 cos θ, cos² + sin² = 1 (gloss) · ellipse foci c² = A² − B²
  (gloss).
- **Assumptions.** b > 0 real; a > b for the ellipse (steps 3–6).
- **Steps.**
  1. *did:* Look far from the origin · *tex:* $\lvert\zeta\rvert\gg b\ \Rightarrow\ z\approx\zeta$ · *why:* b²/ζ → 0 as \lvert ζ\rvert grows; the map
     leaves the far field — the free stream — unchanged. · *plain:* Far away the two planes coincide.
  2. *did:* Map the circle ζ = b e^{iθ} · *tex:* $z=be^{i\theta}+be^{-i\theta}=2b\cos\theta$ · *why:* b²/(b e^{iθ}) = b e^{−iθ}; the sum of the
     two exponentials is twice the cosine. · *plain:* The circle of radius b becomes the segment from −2b to 2b, traced
     twice. · *set:* {ab: 1} · *watch:* "the plate".
  3. *did:* Map a circle of radius a > b · *tex:* $z=ae^{i\theta}+\frac{b^2}{a}e^{-i\theta}$ · *why:* Substitute ζ = a e^{iθ} in (6.65). This is
     (6.66). · *plain:* A larger circle maps to a closed curve.
  4. *did:* Separate real and imaginary parts · *tex:* $x=\Big(a+\frac{b^2}{a}\Big)\cos\theta,\qquad y=\Big(a-\frac{b^2}{a}\Big)\sin\theta$ · *why:* Euler's
     formula on each exponential; e^{−iθ} has sine of the opposite sign. · *plain:* x and y oscillate with different
     amplitudes.
  5. *did:* Eliminate θ · *tex:* $\frac{x^2}{(a+b^2/a)^2}+\frac{y^2}{(a-b^2/a)^2}=\cos^2\theta+\sin^2\theta=1$ · *why:* Divide each by its amplitude, square
     and add (gloss). This is (6.67). · *plain:* The image is an ellipse with semi-axes a + b²/a and a − b²/a. · *live:*
     "A = A, B = B".
  6. *did:* Locate the foci · *tex:* $c^2=\Big(a+\frac{b^2}{a}\Big)^2-\Big(a-\frac{b^2}{a}\Big)^2=4b^2\ \Rightarrow\ c=2b$ · *why:* For an ellipse c² = A² −
     B² (gloss); the cross terms ±2b² differ by 4b². · *plain:* Every ellipse of the family has its foci at ±2b — the
     slit's ends.
  7. *did:* Find the critical points · *tex:* $\frac{dz}{d\zeta}=1-\frac{b^2}{\zeta^2}=0\ \Rightarrow\ \zeta=\pm b$ · *why:* Angles are not kept
     where the derivative vanishes (D19 step 5); these points lie on the circle of radius b, inside or on our circle
     a ≥ b. · *plain:* The slit's ends are where the map pinches. · *watch:* "the × marks are the plate's ends".
  8. *did:* Notice the map is two-to-one · *tex:* $z(\zeta)=z\Big(\frac{b^2}{\zeta}\Big)$ · *why:* Replacing ζ by b²/ζ swaps the two terms of
     (6.65); a point outside the circle \lvert ζ\rvert = b and its inverse point inside give the same z. · *plain:* Inside
     and outside the circle each cover the whole z-plane — we use only the outside.
- **Result.** $z=2b\cos\theta$ for \lvert ζ\rvert = b (a slit of length 4b); $\frac{x^2}{(a+b^2/a)^2}+\frac{y^2}{(a-b^2/a)^2}=1$ *(6.67)* for \lvert ζ\rvert =
  a, foci ±2b — *in words:* the Zhukhovsky map flattens circles into ellipses sharing the same foci, and the critical
  circle into a flat plate.
- **Check.** b = 1, a = 1.2: A = 2.033, B = 0.367, c² = 4.134 − 0.134 = 4 ✓. a → ∞: A/a, B/a → 1 (a circle again).
- **What it means.** Any flow round a circle becomes flow round an ellipse (D21); an off-centre circle becomes an airfoil
  (Ch. 14).
- **Traps.** The map is two-to-one; we keep the exterior of the circle. ζ = ±b are critical points: the plate's ends,
  where the speed will blow up unless the Kutta condition holds.

### D21 · Back from z to ζ: the right square root, and the flow round an elliptic cylinder (6.68)–(6.69) — ★★★, 11 steps, in C11 (notebook · `conformal_joukowski`)
- **Goal.** Invert the Zhukhovsky map so that each point z outside the ellipse is sent to the point ζ outside the circle,
  and compute the velocity round the elliptic cylinder with circulation by the chain rule.
- **Start.** $z=\zeta+\frac{b^2}{\zeta}$ *(6.65)* and the circle flow $w=U\Big(\zeta+\frac{a^2}{\zeta}\Big)+\frac{i\Gamma}{2\pi}\ln(\zeta/a)$ *(6.68)* — *in
  words:* the map and the flow we know in the ζ-plane.
- **Plan.** (1) Solve the quadratic for ζ. (2) Pick the outside root. (3) See numpy fail and fix it. (4) Chain rule for
  the velocity.
- **Tools.** complex square roots and the quadratic formula (primer, in C11) · Vieta's formulas (ch02 P71, reminder) · chain
  rule (ch01 P49, reminder) · derivative of an inverse function (gloss in step 10).
- **Assumptions.** a > b (the ellipse encloses the slit); points z outside the ellipse; a stream along x (the explainer
  adds an angle by rotating U).
- **Steps.**
  1. *did:* Multiply the map by ζ · *tex:* $\zeta^2-z\zeta+b^2=0$ · *why:* ζ ≠ 0 outside the circle; multiplying (6.65) by ζ gives a
     quadratic in ζ for each z. · *plain:* Each z has two candidate ζ's.
  2. *did:* Apply the quadratic formula · *tex:* $\zeta=\tfrac12\big[z\pm(z^2-4b^2)^{1/2}\big]$ · *why:* ζ = [z ± √(z² − 4b²)]/2 for a monic
     quadratic (primer). · *plain:* Two roots, differing in the sign of a square root.
  3. *did:* Multiply the two roots · *tex:* $\zeta_+\zeta_-=b^2\ \Rightarrow\ \lvert\zeta_+\rvert\,\lvert\zeta_-\rvert=b^2$ · *why:* Vieta: the product of the roots
     is the constant term (P71). So unless both lie on \lvert ζ\rvert = b, one is outside the circle and one inside. ·
     *plain:* One root is the physical point, the other its inverse inside the circle.
  4. *did:* Choose the root outside · *tex:* $\zeta(z)=\text{the root with }\lvert\zeta\rvert\ge b$ · *why:* The flow (6.68) lives outside the
     circle; the inside root is excluded (the book's "the negative root … has been excluded"). This gives (6.69). ·
     *plain:* Keep the root that lands in the fluid.
  5. *did:* Test numpy's principal root · *tex:* $z=-3+0.5i,\ b=1:\ \tfrac12\big[z+\sqrt{z^2-4b^2}\big]_{\text{principal}}=-0.362-0.079i,\ \lvert\zeta\rvert=0.37$ · *why:*
     `np.sqrt` returns the root with non-negative real part; for Re z < 0 that makes '+' give the inside root. Our check,
     not the book's. · *plain:* The textbook formula with numpy's square root sends left-half points inside the cylinder. ·
     *set:* {branch: 'principal'} · *watch:* "rose streamlines on the left".
  6. *did:* Split the root into two factors · *tex:* $(z^2-4b^2)^{1/2}=\sqrt{z-2b}\,\sqrt{z+2b}$ · *why:* Each principal factor has its cut
     along the real axis to the left of ±2b; in the product the two cuts cancel beyond −2b, leaving a cut only on the slit
     [−2b, 2b]. · *plain:* A square root whose only jump is hidden inside the body. · *set:* {branch: 'outside'}
  7. *did:* Check that it behaves like z far away · *tex:* $\sqrt{z-2b}\sqrt{z+2b}\approx z\ \Rightarrow\ \zeta\approx z\quad(\lvert z\rvert\to\infty)$ · *why:* Both
     factors tend to √z, whose product is z on every side, so '+' gives ζ ≈ z, the outside root, in all four quadrants
     (z = −3 + 0.5i: \lvert ζ\rvert = 2.70 ✓). · *plain:* The product form picks the outside root everywhere.
  8. *did:* Compose the flow · *tex:* $w(z)=W(\zeta(z)),\quad W(\zeta)=U\Big(\zeta+\frac{a^2}{\zeta}\Big)+\frac{i\Gamma}{2\pi}\ln\frac\zeta a$ · *why:* A map carries
     a flow net to a flow net (D19); the same w value at a point and at its image means the same φ and ψ. · *plain:* The
     flow round the ellipse is the circle flow read at ζ(z).
  9. *did:* Apply the chain rule · *tex:* $u-iv=\frac{dw}{dz}=\frac{dW}{d\zeta}\frac{d\zeta}{dz}$ · *why:* Derivative of a composition of analytic
     functions (the book's last line). · *plain:* Velocity = circle velocity × the map's local stretch.
  10. *did:* Invert the map's derivative · *tex:* $\frac{d\zeta}{dz}=\frac{1}{dz/d\zeta}=\frac{1}{1-b^2/\zeta^2}$ · *why:* The derivative of an inverse
      function is the reciprocal of the derivative (gloss), from (6.65). · *plain:* Where the map squeezes, the flow
      speeds up.
  11. *did:* Write the velocity round the ellipse · *tex:* $u-iv=\frac{U(1-a^2/\zeta^2)+i\Gamma/2\pi\zeta}{1-b^2/\zeta^2}$ · *why:* dW/dζ from (6.68)
      differentiated, times step 10; for a > b the denominator vanishes only at ζ = ±b, inside the circle. · *plain:* A
      finite velocity everywhere round the ellipse. · *live:* "u − iv = uv".
- **Result.** $\zeta=\tfrac12z+\tfrac12(z^2-4b^2)^{1/2}$ with the root outside the circle *(6.69)*, coded as
  $\zeta=\tfrac12\big[z+\sqrt{z-2b}\sqrt{z+2b}\big]$, and $u-iv=\frac{dW}{d\zeta}\frac{d\zeta}{dz}$ — *in words:* invert the map with the square root
  cut along the slit, and the circle's flow gives the flow round the ellipse.
- **Check.** `joukowski_inverse(joukowski(ζ)) = ζ` for 10⁴ random ζ outside the circle in all four quadrants (tests);
  normal velocity on the ellipse ≲ 1e-12 and Blasius lift ρUΓ (C10, C11 code) ✓. a = b: the denominator vanishes at the
  plate's ends — infinite speed unless the numerator also vanishes there.
- **sympy check intent (★★★, `check_src`, a wrong-variant check):** `zeta, z, b = sp.symbols('zeta z b')` · `roots =
  sp.solve(zeta**2 - z*zeta + b**2, zeta)` (step 2) · `print(sp.simplify(roots[0]*roots[1]))` (→ b**2, step 3) · numeric
  part: `zz = (np.linspace(-4, -0.5, 41)[:, None] + 1j*np.linspace(-2, 2, 41)[None, :]).ravel()` (left half-plane grid,
  points outside the slit) · `zz = zz[np.abs(cm.joukowski_inverse(zz, 1.0)) > 1.0]` · `bad = np.abs(0.5*(zz +
  np.sqrt(zz**2 - 4))) < 1` (the principal-root variant) · `print(bad.mean())` (→ 1.0: it fails at every left-half point)
  · `good = cm.joukowski_inverse(zz, 1.0)`; `assert np.allclose(cm.joukowski(good, 1.0), zz)` and `assert
  np.all(np.abs(good) >= 1)`. Every line commented.
- **What it means.** Conformal maps turn one solved flow into a family; the one numerical trap is the branch of the
  inverse. With a = b and Γ chosen to cancel the numerator at the trailing edge, this becomes the flat-plate airfoil of
  Ch. 14.
- **Traps.** The two roots multiply to b² — one inside, one outside. numpy's principal √(z² − 4b²) returns the inside root
  for Re z < 0. dζ/dz blows up at the slit ends when a = b.

### D22 · Laplace on a grid: the five-point average rule (6.70)–(6.72) — ★, 6 steps, in C12 (notebook · `laplace_relaxation`)
- **Goal.** Replace Laplace's equation by an algebraic rule on a square grid, and find how accurate the rule is.
- **Start.** $\Big(\frac{\partial\psi}{\partial x}\Big)_{i,j}\simeq\frac{\psi_{i+\frac12,j}-\psi_{i-\frac12,j}}{\Delta x}$ (the half-point difference of R19) — *in
  words:* a slope from values half a cell to each side.
- **Plan.** (1) Difference the difference in x. (2) Estimate the error by Taylor series. (3) Same in y and add. (4) Take
  Δx = Δy.
- **Tools.** grid values (R19) · Taylor series (ch01 P26, reminder) · Laplace (6.5).
- **Assumptions.** ψ smooth (four derivatives, step 3); uniform square grid (step 5).
- **Steps.**
  1. *did:* Difference the first differences · *tex:* $\Big(\frac{\partial^2\psi}{\partial x^2}\Big)_{i,j}\simeq\frac1{\Delta x}\Big[\Big(\frac{\partial\psi}{\partial x}\Big)_{i+\frac12,j}-\Big(\frac{\partial\psi}{\partial x}\Big)_{i-\frac12,j}\Big]$ ·
     *why:* A second derivative is the rate of change of the first; apply the half-point rule to it. · *plain:* Curvature
     = change of slope across a cell.
  2. *did:* Use whole-point values for each slope · *tex:* $\Big(\frac{\partial^2\psi}{\partial x^2}\Big)_{i,j}\simeq\frac{\psi_{i+1,j}-2\psi_{i,j}+\psi_{i-1,j}}{\Delta x^2}$ · *why:*
     The slope at i + ½ uses ψ_{i+1} and ψ_i; the slope at i − ½ uses ψ_i and ψ_{i−1}; the half points never need values.
     This is (6.70). · *plain:* Left, centre twice, right.
  3. *did:* Measure the error with Taylor · *tex:* $\psi_{i+1,j}+\psi_{i-1,j}-2\psi_{i,j}=\Delta x^2\psi_{xx}+\frac{\Delta x^4}{12}\psi_{xxxx}+\dots$ · *why:* Expand ψ(x
     ± Δx) to fourth order; the odd powers cancel in the sum. So the error of (6.70) is (Δx²/12)ψ_xxxx: second-order
     accurate (the book's "first-order" means first-derivative). · *plain:* Halving the grid spacing cuts the error by
     four.
  4. *did:* Do the same in y and add · *tex:* $\frac{\psi_{i+1,j}-2\psi_{i,j}+\psi_{i-1,j}}{\Delta x^2}+\frac{\psi_{i,j+1}-2\psi_{i,j}+\psi_{i,j-1}}{\Delta y^2}=0$ · *why:*
     (6.71) is (6.70) in y; their sum approximates ∇²ψ, which is zero by (6.5). · *plain:* The grid version of Laplace's
     equation.
  5. *did:* Take Δx = Δy and solve for ψ_{i,j} · *tex:* $\psi_{i,j}=\tfrac14\big[\psi_{i-1,j}+\psi_{i+1,j}+\psi_{i,j-1}+\psi_{i,j+1}\big]$ · *why:* Multiply
     by Δx², move the −4ψ_{i,j} to the other side, divide by 4. This is (6.72). · *plain:* Every value is the average of
     its four neighbours. · *live:* "¼(nW + nE + nS + nN) = avg".
  6. *did:* Read off a consequence · *tex:* $\min(\text{neighbours})\le\psi_{i,j}\le\max(\text{neighbours})$ · *why:* An average lies between
     the smallest and largest of the values averaged; so no interior node can be a strict maximum or minimum (discrete
     maximum principle). · *plain:* A harmonic grid function takes its extremes on the boundary, like a stretched
     membrane.
- **Result.** $\psi_{i,j}=\tfrac14[\psi_{i-1,j}+\psi_{i+1,j}+\psi_{i,j-1}+\psi_{i,j+1}]$ *(6.72)*, with truncation error $(\Delta x^2/12)(\psi_{xxxx}+\psi_{yyyy})$ —
  *in words:* on a square grid Laplace's equation says "be the average of your neighbours", to second-order accuracy.
- **Check.** Observed order 2.00 ± 0.05 on sin(πx) sinh(πy) (R21 code) ✓. ψ = xy has zero fourth derivatives, so it
  satisfies (6.72) exactly: ψ₂₂ = ¼(0 + 2 + 0 + 2) = 1 ✓.
- **What it means.** A PDE has become one linear equation per node; solving them is D23. The same stencil with a source
  term solves Poisson problems (Ch. 10, Ch. 13).
- **Traps.** The average rule needs Δx = Δy. "First-order central differences" in the book means first-derivative, not
  first-order accurate.

### D23 · The four-point system (6.73) as Aψ = b, and why Gauss–Seidel converges — ★★, 9 steps, in C12 (notebook · `laplace_relaxation`)
- **Goal.** Write the four interior equations of the 16-point grid as a matrix system, and understand why the book's
  sweep ("always use the latest value") converges and how fast.
- **Start.** (6.72) at the nodes (2,2), (3,2), (2,3), (3,3), i.e. $\psi_{2,2}=\tfrac14[\psi^B_{1,2}+\psi_{3,2}+\psi^B_{2,1}+\psi_{2,3}]$, … *(6.73)* —
  *in words:* four averages, twelve known boundary values.
- **Plan.** (1) Matrix form. (2) Diagonal dominance. (3) Jacobi vs Gauss–Seidel. (4) Error per sweep: the spectral radius.
  (5) When to stop.
- **Tools.** iterative solvers: Jacobi, Gauss–Seidel, SOR (primer, in C12) · eigenvalues (ch02 P80, reminder) ·
  `np.linalg.solve` (ch01 P57, reminder).
- **Assumptions.** Dirichlet data (all boundary values given); the ordering ψ₂₂, ψ₃₂, ψ₂₃, ψ₃₃.
- **Steps.**
  1. *did:* Move the unknowns to the left · *tex:* $4\psi_{2,2}-\psi_{3,2}-\psi_{2,3}=\psi^B_{1,2}+\psi^B_{2,1}$ (and three more) · *why:* Multiply (6.73)
     by 4 and collect unknowns on one side, knowns on the other. · *plain:* Each node: four times itself minus its
     unknown neighbours equals its known neighbours.
  2. *did:* Write it as Aψ = b · *tex:* $A=\begin{pmatrix}4&-1&-1&0\\-1&4&0&-1\\-1&0&4&-1\\0&-1&-1&4\end{pmatrix},\quad\boldsymbol\psi=(\psi_{22},\psi_{32},\psi_{23},\psi_{33})$ · *why:*
     Each row is one node's equation; the −1s mark neighbours that are unknowns. · *plain:* A small linear system, solvable
     directly (np.linalg.solve). · *live:* "b = (0, 3, 3, 12) for ψ = xy".
  3. *did:* Note the diagonal dominance · *tex:* $\lvert a_{ii}\rvert=4>\sum_{j\ne i}\lvert a_{ij}\rvert=2$ · *why:* Each row's diagonal beats the sum of its
     other entries; such a matrix is invertible and simple sweeps converge on it. · *plain:* Each node listens mostly to
     itself: the system is well behaved.
  4. *did:* Jacobi: update from the old sweep · *tex:* $\boldsymbol\psi^{(k+1)}=\tfrac14\big(\mathbf b+N\boldsymbol\psi^{(k)}\big),\quad N=4I-A$ · *why:* Solve each row
     for its diagonal unknown using last sweep's values; N holds the neighbour links. · *plain:* Everybody updates at once
     from the old picture. · *set:* {method: 'Jacobi'}
  5. *did:* Gauss–Seidel: use new values at once · *tex:* $\psi^{(k+1)}_{3,2}=\tfrac14\big[\psi^{(k+1)}_{2,2}+\psi^B_{4,2}+\psi^B_{3,1}+\psi^{(k)}_{3,3}\big]$ · *why:*
     The book's rule "always use the latest available value": ψ₂₂ was just updated, so it is used immediately. ·
     *plain:* Information travels further in one sweep. · *live:* "first sweep: 0, 0.75, 0.75, 3.375".
  6. *did:* Follow the error per sweep · *tex:* $\mathbf e^{(k+1)}=G\,\mathbf e^{(k)},\quad\rho_J=\tfrac14\max\lvert\lambda(N)\rvert=\tfrac24=\tfrac12,\quad\rho_{GS}=\rho_J^2=\tfrac14$ ·
     *why:* The error obeys the sweep's linear map G, so its largest eigenvalue size (spectral radius) is the factor
     per sweep; N has eigenvalues ±2, 0, 0, and for this matrix Gauss–Seidel squares Jacobi's factor. · *plain:* Jacobi halves the error each sweep, Gauss–Seidel quarters it. · *live:* "ρ_J = ½, ρ_GS = ¼".
  7. *did:* Count sweeps to a tolerance · *tex:* $k\approx\frac{\ln(10^{-8})}{\ln\rho}:\quad27\ \text{(Jacobi)},\ 14\ \text{(Gauss–Seidel)}$ · *why:* Solve ρᵏ = 10⁻⁸
     for k. On an N × N grid ρ_J ≈ 1 − π²/2N², so sweeps grow like N²; SOR needs only about N. · *plain:* Finer grids converge more slowly; a smarter sweep helps a lot.
  8. *did:* Stop on the residual · *tex:* $r=\max_{i,j}\Big\lvert\tfrac14(\text{neighbours})-\psi_{i,j}\Big\rvert$ · *why:* The change per sweep is about (1
     − ρ) × error: when ρ ≈ 1 a small change hides a large error. The residual measures how badly (6.72) is violated. ·
     *plain:* Ask "is it solved?", not "did it stop moving?". · *watch:* "the dashed change below the solid residual".
  9. *did:* Check with harmonic boundary values · *tex:* $\psi^B=xy\ \Rightarrow\ \boldsymbol\psi=(1,2,2,4)$ · *why:* xy satisfies (6.72) exactly (D22
     step 3), so the converged grid must reproduce it; np.linalg.solve and every sweep method agree. · *plain:* A test
     problem whose answer we know.
- **Result.** $A\boldsymbol\psi=\mathbf b$ with A diagonally dominant; Jacobi's error shrinks by ρ_J = ½ and Gauss–Seidel's by ρ_GS = ¼ per
  sweep on the 16-point grid — *in words:* relaxation converges because each node is dominated by itself; using the newest
  values roughly halves the work, and the residual says when to stop.
- **Check.** Sweeps to a residual of 10⁻¹⁰ from zero (C12 code): Jacobi 34, Gauss–Seidel 19, SOR (ω_opt = 1.072) 12 ✓;
  solution 1, 2, 2, 4 ✓.
- **What it means.** Every elliptic solver of Ch. 10 (pressure Poisson, ω–ψ) and every bounded-basin streamfunction
  inversion in ocean models is this idea made faster (SOR, multigrid).
- **Traps.** Jacobi uses only old values; Gauss–Seidel uses new ones in the same sweep. Stopping on the change instead of
  the residual stops too early on fine grids.

### D24 · The Stokes stream function (6.75) and its field equation (6.77), which is not Laplace's — ★★, 8 steps, in C13 (notebook)
- **Goal.** Derive the velocity from the axisymmetric (Stokes) stream function and the equation it obeys in irrotational
  flow — and see why it is not the Laplace equation.
- **Start.** $\mathbf u=\nabla\chi\times\nabla\psi,\ \chi=-\varphi$ (φ here the azimuth angle) and $\omega_\varphi=\frac{\partial u_R}{\partial z}-\frac{\partial u_z}{\partial R}$
  *(6.76)* — *in words:* Ch. 4's two-stream-function form with the first chosen as the angle.
- **Plan.** (1) ∇χ. (2) The cross product. (3) Components (6.75). (4) Substitute into ω_φ. (5) Compare with ∇²ψ.
- **Tools.** Stokes stream function (R24) · azimuthal vorticity (R25) · cross products of cylindrical unit vectors (gloss) ·
  cylindrical unit vectors (ch03 P88, reminder).
- **Assumptions.** Axisymmetric, no swirl (u_φ = 0, nothing depends on φ); irrotational in step 7.
- **Steps.**
  1. *did:* Take the gradient of χ = −φ · *tex:* $\nabla\chi=-\frac{1}{R}\mathbf e_\varphi$ · *why:* In cylindrical coordinates the φ-part of
     the gradient is (1/R)∂/∂φ e_φ, and ∂φ/∂φ = 1. · *plain:* The first stream function points round the axis.
  2. *did:* Write ∇ψ for an axisymmetric ψ · *tex:* $\nabla\psi=\frac{\partial\psi}{\partial R}\mathbf e_R+\frac{\partial\psi}{\partial z}\mathbf e_z$ · *why:* ψ does not depend
     on φ. · *plain:* ψ changes only with distance from the axis and along it.
  3. *did:* Cross the two · *tex:* $\mathbf u=-\frac1R\,\mathbf e_\varphi\times\Big(\frac{\partial\psi}{\partial R}\mathbf e_R+\frac{\partial\psi}{\partial z}\mathbf e_z\Big)$ · *why:* u = ∇χ ×
     ∇ψ with steps 1–2. · *plain:* The velocity is perpendicular to both gradients.
  4. *did:* Evaluate the cylindrical cross products · *tex:* $u_R=-\frac1R\frac{\partial\psi}{\partial z},\qquad u_z=\frac1R\frac{\partial\psi}{\partial R}$ · *why:*
     e_φ × e_R = −e_z and e_φ × e_z = e_R because (R, φ, z) is right-handed and cyclic (gloss); the two minus
     signs combine. This is (6.75). ·
     *plain:* Axial speed from how ψ grows outward, radial speed from how it changes along the axis.
  5. *did:* Differentiate u_R in z · *tex:* $\frac{\partial u_R}{\partial z}=-\frac1R\frac{\partial^2\psi}{\partial z^2}$ · *why:* R is held fixed in ∂/∂z, so 1/R
     comes out. · *plain:* The first half of the vorticity.
  6. *did:* Differentiate u_z in R · *tex:* $\frac{\partial u_z}{\partial R}=\frac{\partial}{\partial R}\Big(\frac1R\frac{\partial\psi}{\partial R}\Big)$ · *why:* Here 1/R depends
     on R, so it must stay inside the derivative (product rule, P38). · *plain:* The second half.
  7. *did:* Substitute into ω_φ and set it zero · *tex:* $\frac{\partial}{\partial R}\Big(\frac1R\frac{\partial\psi}{\partial R}\Big)+\frac1R\frac{\partial^2\psi}{\partial z^2}=-\omega_\varphi=0$ ·
     *why:* (6.76) with steps 5–6, the sign flipped; irrotational flow has ω_φ = 0. This is (6.77). · *plain:* The field
     equation of the Stokes stream function.
  8. *did:* Compare with the Laplacian · *tex:* $R\times(6.77)=\psi_{RR}-\frac{\psi_R}{R}+\psi_{zz}\ \ne\ \nabla^2\psi=\psi_{RR}+\frac{\psi_R}{R}+\psi_{zz}$ · *why:*
     Expand ∂/∂R(ψ_R/R) = ψ_RR/R − ψ_R/R² and multiply by R; the first-derivative term has the opposite sign. · *plain:*
     Not Laplace's equation — so ψ is not the imaginary part of an analytic function, and complex variables do not apply.
- **Result.** $u_R=-\frac1R\frac{\partial\psi}{\partial z},\ u_z=\frac1R\frac{\partial\psi}{\partial R}$ *(6.75)* and $\frac{\partial}{\partial R}\Big(\frac1R\frac{\partial\psi}{\partial R}\Big)+\frac1R\frac{\partial^2\psi}{\partial z^2}=0$
  *(6.77)* — *in words:* one stream function still describes axisymmetric flow, but its equation is not Laplace's.
- **Check.** ψ = ½UR² (uniform stream): (6.77) gives ∂/∂R(U) = 0 ✓ while ∇²ψ = 2U ≠ 0. The sphere's ψ (6.89) passes
  `stokes_operator_residual`; the control Rz fails it (residual −z/R²; note R²z is itself a solution — axisymmetric stagnation flow) (C13 code).
- **What it means.** Axisymmetric flows are built from elements in real coordinates (D25), not from complex functions;
  φ, which does obey Laplace (6.80), stays useful. The zonal-mean overturning streamfunction in climate science has this
  Stokes geometry.
- **Traps.** u_R = −(1/R)∂ψ/∂z carries the minus. The 1/R inside ∂/∂R must not be pulled out.

### D25 · 3-D elements and the sphere: why d = 2πa³U, and C_p = 1 − (9/4) sin²θ (6.86)–(6.91) — ★★, 10 steps, in C13 (notebook)
- **Goal.** Build the 3-D point source and doublet from mass conservation and a limit, add a doublet to a stream to make a
  sphere, and find its surface speed and pressure.
- **Start.** $u_r=\frac{1}{r^2\sin\theta}\frac{\partial\psi}{\partial\theta}=\frac{\partial\phi}{\partial r},\ u_\theta=-\frac{1}{r\sin\theta}\frac{\partial\psi}{\partial r}=\frac1r\frac{\partial\phi}{\partial\theta}$
  *(6.83)* and the uniform stream $\phi=Ur\cos\theta,\ \psi=\tfrac12Ur^2\sin^2\theta$ *(6.86)* — *in words:* spherical velocities and the
  simplest flow.
- **Plan.** (1) Check the stream. (2) The point source from its flux. (3) The doublet as a merged pair. (4) Stream +
  doublet with ψ = 0 on r = a. (5) Surface speed and C_p.
- **Tools.** spherical velocities (N79) · flux of a 3-D point source (gloss) · multivariable Taylor (ch03 P98, reminder) ·
  a limit with a product held fixed (primer, in C04) · C_p (N28).
- **Assumptions.** Axisymmetric about z (horizontal, along the stream); r > 0.
- **Steps.**
  1. *did:* Check the uniform stream · *tex:* $u_z=\frac1R\frac{\partial}{\partial R}\Big(\tfrac12UR^2\Big)=U,\qquad u_R=0$ · *why:* (6.75) with ψ = ½UR²;
     R = r sin θ turns it into ½Ur² sin²θ. This is (6.86). · *plain:* A stream U along the axis.
  2. *did:* Get the source speed from its flux · *tex:* $4\pi r^2u_r=Q\ \Rightarrow\ u_r=\frac{Q}{4\pi r^2}$ · *why:* All Q m³/s crosses every
     sphere round the source (gloss). · *plain:* The outflow thins out over ever larger spheres.
  3. *did:* Integrate for φ and ψ · *tex:* $\phi=-\frac{Q}{4\pi r},\qquad\frac{\partial\psi}{\partial\theta}=r^2\sin\theta\,u_r=\frac{Q}{4\pi}\sin\theta\ \Rightarrow\ \psi=-\frac{Q}{4\pi}\cos\theta$ ·
     *why:* ∂φ/∂r = u_r and (6.83) for ψ; constants of integration set to zero. This is (6.87). · *plain:* The 3-D source,
     in both functions.
  4. *did:* Merge a source and a sink along z · *tex:* $\phi=-\frac{Q}{4\pi\lvert\mathbf x+\varepsilon\mathbf e_z\rvert}+\frac{Q}{4\pi\lvert\mathbf x-\varepsilon\mathbf e_z\rvert}\approx\frac{Q}{4\pi}\frac{2\varepsilon z}{r^3}$ ·
     *why:* Source at −ε e_z, sink at +ε e_z; Taylor 1/\lvert x ∓ εe_z\rvert ≈ 1/r ± εz/r³ (P98); the 1/r terms cancel. ·
     *plain:* The pair's potential, to first order in ε.
  5. *did:* Hold d = 2εQ fixed · *tex:* $\phi=\frac{d}{4\pi}\frac{z}{r^3}=\frac{d}{4\pi r^2}\cos\theta,\qquad\psi=-\frac{d}{4\pi r}\sin^2\theta$ · *why:* The limit with
     2εQ fixed (primer); ψ from (6.83): ∂ψ/∂θ = r² sin θ ∂φ/∂r = −(d/2πr) sin θ cos θ. The dipole is −d e_z. This is
     (6.88). · *plain:* The 3-D doublet.
  6. *did:* Add a stream and an opposing doublet · *tex:* $\psi=\tfrac12Ur^2\sin^2\theta-\frac{d}{4\pi r}\sin^2\theta$ · *why:* Superposition (C03);
     the doublet opposes the stream as for the cylinder (D09). · *plain:* A stream disturbed by a doublet.
  7. *did:* Make r = a a stream surface · *tex:* $\tfrac12Ua^2=\frac{d}{4\pi a}\ \Rightarrow\ d=2\pi a^3U$ · *why:* ψ(a, θ) = 0 for every θ
     requires the bracket multiplying sin²θ to vanish (D05). · *plain:* The doublet strength that makes a sphere.
  8. *did:* Substitute back · *tex:* $\psi=\tfrac12Ur^2\Big(1-\frac{a^3}{r^3}\Big)\sin^2\theta,\qquad\phi=Ur\Big(1+\frac{a^3}{2r^3}\Big)\cos\theta$ · *why:* d/4πr =
     ½Ua³/r. This is (6.89). · *plain:* The flow round a sphere.
  9. *did:* Differentiate for the velocity · *tex:* $u_r=U\Big[1-\Big(\frac ar\Big)^3\Big]\cos\theta,\qquad u_\theta=-U\Big[1+\frac12\Big(\frac ar\Big)^3\Big]\sin\theta$ · *why:*
     (6.83) applied to step 8: ∂ψ/∂r = ½U sin²θ(2r + a³/r²). This is (6.90). · *plain:* The disturbance dies like
     (a/r)³ — faster than the cylinder's (a/r)².
  10. *did:* Evaluate on the surface · *tex:* $u_\theta(a,\theta)=-\tfrac32U\sin\theta\ \Rightarrow\ C_p=1-\frac94\sin^2\theta$ · *why:* r = a gives u_r = 0 and
      \lvert u\rvert = (3/2)U\lvert sin θ\rvert; C_p = 1 − \lvert u\rvert²/U² (6.32). This is (6.91). · *plain:* The fastest
      surface speed is 1.5U, the lowest C_p −1.25.
- **Result.** Source $\phi=-Q/4\pi r$ *(6.87)*, doublet $\phi=\frac{d}{4\pi r^2}\cos\theta$ *(6.88)*, sphere with $d=2\pi a^3U$ *(6.89)*,
  $u_\theta=-U[1+\frac12(a/r)^3]\sin\theta$ *(6.90)*, $C_p=1-\frac94\sin^2\theta$ *(6.91)* — *in words:* the 3-D flow escapes sideways, so the
  sphere disturbs the stream less than the cylinder does.
- **Check.** Units: Stokes ψ [m/s × m²] = m³/s ✓, d [m⁴/s] ✓. Numbers: sphere max speed 1.5U, C_p min −1.25 vs cylinder 2U,
  −3; at r = 2a: 0.125 vs 0.25 ✓ (C13 code). Spherical Laplace (6.85) of both potentials = 0 (R30 code).
- **What it means.** The sphere is the 3-D reference body: the Stokes-flow comparison of Ch. 8, the separation of Ch. 9,
  and the added mass of C15 all start here. Still no drag (fore–aft symmetric C_p).
- **Traps.** φ = −Q/4πr (minus for a source). ψ of the source is −(Q/4π) cos θ plus any constant. The 3-D doublet's ψ has
  sin²θ, not sin θ.

### D26 · A line sink and the airship: (6.93) → (6.94) → (6.95) — ★★, 10 steps, in C14 (notebook · `axial_singularity_bodies`)
- **Goal.** Add up a uniform line of sinks on the axis in closed form, and combine it with a point source and a stream so
  that the body closes — the airship.
- **Start.** A point source Q at the origin has $\psi=-\frac{Q}{4\pi}\cos\theta$ *(6.87)*, θ measured at the source from +z — *in
  words:* the stream function of one 3-D source.
- **Plan.** (1) One sink element. (2) Sum along the line (6.93). (3) Substitute z − ξ = R cot α and integrate (6.94). (4)
  Close the body and add the stream (6.95).
- **Tools.** point-source ψ (N81) · substitution in an integral (ch03 P106, reminder) and the cot substitution (gloss) ·
  superposition (C03).
- **Assumptions.** Uniform sink density k [m²/s] on 0 ≤ ξ ≤ a; P = (z, R) off the axis; θ at O, α at the element, α₁ at A
  are the angles between +z and the lines to P (Fig. 6.28 geometry, our notation as the book's).
- **Steps.**
  1. *did:* Take one sink element kdξ · *tex:* $d\psi_{\text{sink}}=+\frac{k\,d\xi}{4\pi}\cos\alpha$ · *why:* A sink is a source of strength −kdξ;
     (6.87) about the element's own position, with α its angle to P, gives the plus sign. · *plain:* Each short piece of
     sink contributes like a point sink.
  2. *did:* Add the elements from O to A · *tex:* $\psi_{\text{sink}}=\frac{k}{4\pi}\int_0^a\cos\alpha\,d\xi$ · *why:* Superposition: the stream
     function of a line of sinks is the integral of the elements. This is (6.93). · *plain:* The whole line's stream
     function, still as an integral.
  3. *did:* Relate ξ to α at fixed P · *tex:* $z-\xi=R\cot\alpha$ · *why:* The line from the element (ξ, 0) to P = (z, R) has
     horizontal run z − ξ and height R; cot α = run/height. · *plain:* Walking along the sink changes the angle at which
     P is seen.
  4. *did:* Differentiate at fixed P · *tex:* $-d\xi=R\,d(\cot\alpha)=-\frac{R}{\sin^2\alpha}d\alpha\ \Rightarrow\ d\xi=\frac{R\,d\alpha}{\sin^2\alpha}$ · *why:* z and R
     do not change along the sink; d(cot α)/dα = −1/sin²α (gloss). · *plain:* A step along the sink is a step in angle.
  5. *did:* Change the limits · *tex:* $\xi=0\leftrightarrow\alpha=\theta,\qquad\xi=a\leftrightarrow\alpha=\alpha_1$ · *why:* At O the angle is θ, at A it is α₁
     (Fig. 6.28). · *plain:* The integral now runs over the angles seen from P. · *watch:* "α runs from θ at O to α₁ at A".
  6. *did:* Substitute · *tex:* $\psi_{\text{sink}}=\frac{kR}{4\pi}\int_\theta^{\alpha_1}\frac{\cos\alpha}{\sin^2\alpha}d\alpha=\frac{kR}{4\pi}\int_\theta^{\alpha_1}\frac{d(\sin\alpha)}{\sin^2\alpha}$ · *why:* Steps
     4–5 in (6.93); cos α dα = d(sin α). · *plain:* An integral we can do at sight.
  7. *did:* Integrate · *tex:* $\psi_{\text{sink}}=\frac{kR}{4\pi}\Big[-\frac{1}{\sin\alpha}\Big]_\theta^{\alpha_1}=\frac{kR}{4\pi}\Big[\frac1{\sin\theta}-\frac1{\sin\alpha_1}\Big]$ · *why:*
     The antiderivative of 1/s² is −1/s (s = sin α); evaluate it between the two limits. · *plain:* The line sink in
     closed form.
  8. *did:* Use the geometry R/sin = distance · *tex:* $\frac{R}{\sin\theta}=r,\quad\frac{R}{\sin\alpha_1}=r_1\ \Rightarrow\ \psi_{\text{sink}}=\frac{k}{4\pi}(r-r_1)$ · *why:*
     r and r₁ are the distances from O and A to P; R is the height of P above the axis. This is (6.94). · *plain:* The
     line sink's ψ is k/4π times the difference of the distances to its ends.
  9. *did:* Close the body: Q = ak · *tex:* $Q-ak=0$ · *why:* The line swallows ka m³/s; a closed body cannot emit fluid, so the
     source must supply exactly that (N86). · *plain:* Source and sink balance. · *live:* "Q − ak = closure".
  10. *did:* Add the source and the stream · *tex:* $\psi=-\frac{Q}{4\pi}\cos\theta+\frac{Q}{4\pi a}(r-r_1)+\tfrac12Ur^2\sin^2\theta$ · *why:* Superposition of
      (6.87), (6.94) with k = Q/a, and (6.86). This is (6.95). · *plain:* The airship's stream function; its body is the
      stream surface through the axis stagnation points. · *live:* "ψ = 0 surface: length L m".
- **Result.** $\psi_{\text{sink}}=\frac{k}{4\pi}(r-r_1)$ *(6.94)* and $\psi=-\frac{Q}{4\pi}\cos\theta+\frac{Q}{4\pi a}(r-r_1)+\frac12Ur^2\sin^2\theta$ *(6.95)* — *in
  words:* a point source at the nose and a uniform line sink behind it, in a stream, draw a closed streamlined hull.
- **Check.** Closed form = `quad` of (6.93) to 1e-10 (C14 code) ✓. Numbers: Q = 1 m³/s, a = 1 m, U = 1 m/s → k = 1 m²/s,
  nose at z = −0.252 m, tail at 1.070 m, length 1.322 m (`brentq` on u_z = 0 along the axis) ✓. Far from the line (r ≫ a)
  r − r₁ → a cos θ and ψ_sink → (ka/4π) cos θ = −ψ of a point source Q: the source is cancelled far away ✓.
- **What it means.** Distributing the sink tapers the tail; a non-uniform k(ξ) would draw other shapes — which D27 turns
  into a design method.
- **Traps.** A sink's dψ has the opposite sign of a source's. d(cot α) = −dα/sin²α reverses the direction of the limits.
  R/sin θ = r and R/sin α₁ = r₁ (distances, not heights).

### D27 · The axial singularity method: N unknown strengths from ψ = 0 at N body points — ★, 6 steps, in C14 (notebook · `axial_singularity_bodies`)
- **Goal.** Turn the airship construction round: for a given body of revolution, find the strengths of sources and sinks
  on its axis that make the body a stream surface.
- **Start.** One segment from ξ_{n−1} to ξ_n with uniform source density k_n [m²/s]; (6.94) for a sink — *in words:* D26's
  line, cut into N pieces with unknown strengths.
- **Plan.** (1) One segment's ψ at a body point. (2) Sum and add the stream. (3) Demand ψ = 0 at N points. (4) Solve and
  check closure.
- **Tools.** (6.94) (D26) · collocation and the condition number (primer, in C14) · `np.linalg.solve` (ch01 P57, reminder).
- **Assumptions.** Body of revolution; singularities on the axis inside it; N body points ↔ N segments.
- **Steps.**
  1. *did:* Write one source segment's ψ · *tex:* $\psi_{mn}=-\frac{k_n}{4\pi}\big(r^m_{n-1}-r^m_n\big)$ · *why:* (6.94) is for a sink from its
     start to its end; a source has the opposite sign; r^m_{n−1}, r^m_n are the distances from body point m to the
     segment's ends. · *plain:* Each segment's contribution at each body point.
  2. *did:* Superpose N segments and the stream · *tex:* $\psi_m=-\sum_{n=1}^N\frac{k_n}{4\pi}\big(r^m_{n-1}-r^m_n\big)+\tfrac12UR_m^2$ · *why:*
     Superposition (C03) with the uniform stream ψ = ½UR² (6.86). · *plain:* ψ at body point m as a sum of N unknown
     pieces.
  3. *did:* Demand the body be ψ = 0 · *tex:* $\psi_m=0,\qquad m=1,\dots,N$ · *why:* The axis ahead of and behind the body has ψ = 0;
     the body surface joins it, so it is the same stream surface. Collocation: one condition per unknown (primer). ·
     *plain:* N equations for N strengths.
  4. *did:* Write it as a matrix system · *tex:* $\sum_nA_{mn}k_n=\tfrac12UR_m^2,\qquad A_{mn}=\frac{r^m_{n-1}-r^m_n}{4\pi}$ · *why:* Move the stream term
     to the right; the geometry-only coefficients form the influence matrix A. · *plain:* A linear system: shape in A,
     stream on the right. · *set:* {tgt: rankine, N: 20}
  5. *did:* Solve for the strengths · *tex:* $\mathbf k=A^{-1}\big(\tfrac12U\mathbf R^2\big)$ · *why:* A is square and (for smooth slender bodies)
     invertible; `np.linalg.solve` or, as the book suggests, the iteration of §6.7. · *plain:* One solve gives every
     segment's strength.
  6. *did:* Check closure and conditioning · *tex:* $\sum_nk_n\Delta\xi\approx0,\qquad\text{cond}(A)$ · *why:* A closed body must have zero net
     source (N86) — a free check, since it was not imposed; a large condition number warns that blunt bodies make the
     system nearly singular. · *plain:* The strengths should cancel; if cond is huge, trust the shape, not the bars. ·
     *set:* {tgt: sphere} · *watch:* "cond explodes".
- **Result.** $\psi_m=-\sum_{n=1}^N\frac{k_n}{4\pi}(r^m_{n-1}-r^m_n)+\frac12UR_m^2=0$ at N body points — an N × N linear system for the
  axial strengths — *in words:* hide unknown sources on the axis, ask for a streamline at N points, solve.
- **Check.** Rankine-oval target: the recovered k_n converge to a point source and sink as N grows; Σk_nΔξ → 0 (C14 code) ✓.
  Airship target: k_n reproduce the flat line sink of (6.95). Sphere: small residual but cond ≫ 10⁸ — fenced as
  qualitative.
- **What it means.** The inverse ("design") method of the chapter and the ancestor of panel methods (N90, Ch. 14):
  unknown strengths + one condition per point + a linear solve.
- **Traps.** Source vs sink sign. ψ = 0 on the axis and on the body is one surface. Blunt bodies (a sphere is a point
  doublet) make A ill-conditioned as N grows.

### D28 · The potential of a moving sphere: d(t) = 2πa³u_s (6.96)–(6.97) — ★, 5 steps, in C15 (notebook · `added_mass_sphere`)
- **Goal.** Write the velocity potential of a sphere moving through still fluid with an arbitrary velocity u_s(t), and
  check that the fluid moves with the sphere's surface.
- **Start.** $\phi=\Big(\mathbf U-\frac{\mathbf d}{4\pi\lvert\mathbf x\rvert^3}\Big)\cdot\mathbf x$ *(6.92)* — *in words:* the sphere in a stream U, with its
  doublet d = −2πa³U (for U = U e_z, d = −d e_z with d = 2πa³U).
- **Plan.** (1) Move the centre to x_s(t). (2) Remove the stream. (3) Choose the dipole in the sphere's frame. (4)
  Substitute. (5) Check the surface condition.
- **Tools.** coordinate-free potential (N85) · Galilean frames (Ch. 3 §3.3, R14) · no through-flow on a moving body (6.16)
  (N16).
- **Assumptions.** Rigid sphere; fluid at rest far away; at each instant the flow is the steady one for the current
  velocity (ideal flow has no memory).
- **Steps.**
  1. *did:* Centre the sphere at x_s(t) · *tex:* $\mathbf x\to\boldsymbol\xi=\mathbf x-\mathbf x_s(t)$ · *why:* (6.92) is for a sphere at the origin;
     shifting the argument moves it. · *plain:* Measure positions from the sphere's centre.
  2. *did:* Drop the stream: still fluid far away · *tex:* $\phi=-\frac{\mathbf d\cdot(\mathbf x-\mathbf x_s)}{4\pi\lvert\mathbf x-\mathbf x_s\rvert^3}$ · *why:* With U = 0
     the fluid far away is at rest, as required; only the doublet remains. This is (6.96). · *plain:* A moving sphere
     looks like a moving doublet.
  3. *did:* Pick the dipole from the sphere's frame · *tex:* $\mathbf d(t)=-2\pi a^3\big(-\mathbf u_s(t)\big)=2\pi a^3\mathbf u_s(t)$ · *why:* Seen from the
     sphere the fluid streams past at −u_s; the steady sphere in a stream V needs d = −2πa³V (6.89, 6.92). · *plain:* The
     dipole points along the sphere's velocity and follows it.
  4. *did:* Substitute d(t) · *tex:* $\phi=-\frac{a^3}{2\lvert\boldsymbol\xi\rvert^3}\mathbf u_s\cdot\boldsymbol\xi$ · *why:* 2πa³/4π = a³/2. This is (6.97). · *plain:* The
     potential depends on time only through x_s(t) and u_s(t).
  5. *did:* Check the surface condition · *tex:* $\frac{\partial\phi}{\partial\lvert\boldsymbol\xi\rvert}\Big\rvert_{\lvert\boldsymbol\xi\rvert=a}=\frac{a^3}{\lvert\boldsymbol\xi\rvert^3}\mathbf u_s\cdot\mathbf e_\xi\Big\rvert_a=\mathbf u_s\cdot\mathbf e_\xi$ ·
     *why:* Write u_s·ξ = \lvert ξ\rvert u_s·e_ξ, so φ ∝ \lvert ξ\rvert⁻²; its radial derivative on \lvert ξ\rvert = a is the fluid's
     normal velocity, which must equal the sphere's (6.16 with U_s = u_s). · *plain:* The fluid moves out of the way
     exactly as fast as the surface arrives. · *watch:* "the normal velocity of the fluid equals the sphere's all
     round".
- **Result.** $\phi(\mathbf x,\mathbf x_s,\mathbf u_s)=-\frac{a^3}{2\lvert\mathbf x-\mathbf x_s\rvert^3}\mathbf u_s\cdot(\mathbf x-\mathbf x_s)$ *(6.97)* — *in words:* a sphere
  moving at u_s carries a dipole 2πa³u_s with it.
- **Check.** Units: a³u_s/ξ² [m²/s] ✓. Numerically n·∇φ = n·u_s at surface points (C15 code: 1.6 and 0.0 for u_s = 2
  e_z) ✓. u_s constant: the pressure is steady in the sphere's frame (Galilean, N102).
- **What it means.** The fluid's motion is set instantly by the sphere's velocity; its time dependence (through u_s and
  x_s) is what makes the pressure of D29 depend on the acceleration.
- **Traps.** In the sphere's frame the oncoming stream is −u_s, so the dipole is +2πa³u_s. w in this section is the
  z-velocity, not the complex potential.

### D29 · The pressure on an arbitrarily moving sphere: (6.99) → (6.105) — ★★★, 14 steps, in C15 (notebook · `added_mass_sphere`)
- **Goal.** Find the pressure at every point of the surface of a sphere moving with any velocity and acceleration through
  still ideal fluid.
- **Start.** $\phi=-\frac{a^3}{2\lvert\boldsymbol\xi\rvert^3}\mathbf u_s\cdot\boldsymbol\xi$ *(6.97)* and $\Big[\frac{\partial\phi}{\partial t}+\frac12\lvert\nabla\phi\rvert^2+\frac p\rho\Big]_{\text{surface}}=\frac{p_\infty}{\rho}$
  *(6.99)* — *in words:* the moving dipole and unsteady Bernoulli between the surface and far away.
- **Plan.** (1) Bernoulli on the surface (6.100). (2) ∂φ/∂t by the chain rule through x_s and u_s (6.101)–(6.102). (3) ∇φ
  and its surface value (6.103)–(6.104). (4) Combine and simplify (6.105). (5) The steady limit (6.106).
- **Tools.** unsteady Bernoulli (R32) · multivariable chain rule along a path (ch03 P91, reminder) · a function of
  x − x_s: ∂/∂x_s = −∇ (gloss) · gradient of u·ξ/\lvert ξ\rvert³ (gloss; ch05 P140 reminder) · \lvert A\rvert² = A·A.
- **Assumptions.** Ideal, constant density; fluid at rest at p∞ far away; rigid sphere; ∂/∂t at a fixed point x.
- **Steps.**
  1. *did:* Evaluate Bernoulli on the surface · *tex:* $\frac{p_a-p_\infty}{\rho}=-\Big(\frac{\partial\phi}{\partial t}\Big)_a-\frac12\lvert\nabla\phi\rvert_a^2$ · *why:* Rearrange
     (6.99); subscript a means on \lvert ξ\rvert = a. This is (6.100). · *plain:* Surface pressure is lowered by the changing
     potential and by the speed.
  2. *did:* Differentiate φ at a fixed point · *tex:* $\frac{\partial\phi}{\partial t}=\frac{\partial\phi}{\partial(x_s)_i}\frac{d(x_s)_i}{dt}+\frac{\partial\phi}{\partial(u_s)_i}\frac{d(u_s)_i}{dt}$ ·
     *why:* At fixed x, φ changes only because x_s and u_s change (chain rule along a path, P91). · *plain:* The potential
     at a point changes because the sphere moves and because it speeds up.
  3. *did:* Use ∂φ/∂x_s = −∇φ · *tex:* $\frac{\partial\phi}{\partial(x_s)_i}=-\frac{\partial\phi}{\partial x_i}=-u_i$ · *why:* φ depends on x and x_s only through ξ = x − x_s,
     so moving x_s forward is like moving x backward (gloss). · *plain:* Moving the sphere is like moving the observer the
     other way.
  4. *did:* Insert dx_s/dt = u_s · *tex:* $\frac{\partial\phi}{\partial(x_s)_i}\frac{d(x_s)_i}{dt}=-\mathbf u\cdot\mathbf u_s$ · *why:* The sphere's velocity is the
     rate of change of its position. · *plain:* The first part of ∂φ/∂t.
  5. *did:* Differentiate with respect to u_s · *tex:* $\frac{\partial\phi}{\partial\mathbf u_s}\cdot\frac{d\mathbf u_s}{dt}=-\frac{a^3}{2\lvert\boldsymbol\xi\rvert^3}\boldsymbol\xi\cdot\frac{d\mathbf u_s}{dt}$ · *why:*
     φ is linear in u_s, so its u_s-gradient is −a³ξ/2\lvert ξ\rvert³. With step 4 this is (6.101). · *plain:* The second
     part: how the potential responds to acceleration.
  6. *did:* Put the point on the surface · *tex:* $\Big(\frac{\partial\phi}{\partial t}\Big)_a=-\mathbf u_a\cdot\mathbf u_s-\frac a2\mathbf e_\xi\cdot\frac{d\mathbf u_s}{dt}$ · *why:* ξ = a e_ξ
     on the surface, so a³/(2a³) × a = a/2. This is (6.102). · *plain:* On the surface ∂φ/∂t has a speed part and an
     acceleration part.
  7. *did:* Take the gradient of φ · *tex:* $\nabla\phi=-\frac{a^3}{2}\Big[-\frac{3(\mathbf x-\mathbf x_s)}{\lvert\mathbf x-\mathbf x_s\rvert^5}\mathbf u_s\cdot(\mathbf x-\mathbf x_s)+\frac{\mathbf u_s}{\lvert\mathbf x-\mathbf x_s\rvert^3}\Big]$ ·
     *why:* Product rule on (u_s·ξ)\lvert ξ\rvert⁻³ with ∇(u_s·ξ) = u_s and ∇\lvert ξ\rvert⁻³ = −3ξ/\lvert ξ\rvert⁵ (gloss). This is (6.103). ·
     *plain:* The fluid velocity round the moving sphere.
  8. *did:* Evaluate on the surface ξ = a e_ξ · *tex:* $\mathbf u_a=-\frac{a^3}{2}\Big[-\frac{3a\mathbf e_\xi}{a^5}\mathbf u_s\cdot a\mathbf e_\xi+\frac{1}{a^3}\mathbf u_s\Big]=\tfrac32(\mathbf u_s\cdot\mathbf e_\xi)\mathbf e_\xi-\tfrac12\mathbf u_s$ ·
     *why:* Substitute and cancel powers of a. ⚠️ The book's middle expression prints −u_s/a³ in the bracket; (6.103) gives
     +u_s/a³, which the final form of (6.104) needs. · *plain:* In front the fluid moves with the sphere; at the side it
     moves backwards at ½u_s. · *watch:* "at θ_s = 90° the fluid moves backwards at ½u_s".
  9. *did:* Substitute into step 1 · *tex:* $\frac{p_a-p_\infty}{\rho}=\mathbf u_a\cdot\mathbf u_s+\frac a2\mathbf e_\xi\cdot\frac{d\mathbf u_s}{dt}-\frac12\lvert\mathbf u_a\rvert^2$ · *why:* The minus of
     (6.100) flips both terms of (6.102). · *plain:* Three terms: two from the speed, one from the acceleration.
  10. *did:* Expand u_a·u_s · *tex:* $\mathbf u_a\cdot\mathbf u_s=\tfrac32(\mathbf u_s\cdot\mathbf e_\xi)^2-\tfrac12\lvert\mathbf u_s\rvert^2$ · *why:* Dot step 8 with u_s. The book skips
      this move. · *plain:* The first speed term.
  11. *did:* Expand \lvert u_a\rvert² carefully · *tex:* $\lvert\mathbf u_a\rvert^2=\tfrac94(\mathbf u_s\cdot\mathbf e_\xi)^2-\tfrac32(\mathbf u_s\cdot\mathbf e_\xi)^2+\tfrac14\lvert\mathbf u_s\rvert^2=\tfrac34(\mathbf u_s\cdot\mathbf e_\xi)^2+\tfrac14\lvert\mathbf u_s\rvert^2$ ·
      *why:* \lvert A − B\rvert² = A·A − 2A·B + B·B with A = (3/2)(u_s·e)e (e a unit vector) and B = ½u_s; the cross term is
      2 × (3/2)(½)(u_s·e)². · *plain:* The surface speed squared.
  12. *did:* Collect the speed terms · *tex:* $\tfrac32c^2-\tfrac12\lvert\mathbf u_s\rvert^2-\tfrac38c^2-\tfrac18\lvert\mathbf u_s\rvert^2=\tfrac98c^2-\tfrac58\lvert\mathbf u_s\rvert^2,\quad c=\mathbf u_s\cdot\mathbf e_\xi$ ·
      *why:* Step 10 minus half of step 11. · *plain:* The speed part of the surface pressure. · *live:* "(9/8)c² −
      (5/8)\lvert u_s\rvert² = p_steady/ρ".
  13. *did:* Factor out ½\lvert u_s\rvert² · *tex:* $\frac{p_a-p_\infty}{\rho}=\frac12\lvert\mathbf u_s\rvert^2\Big(\frac94\frac{(\mathbf u_s\cdot\mathbf e_\xi)^2}{\lvert\mathbf u_s\rvert^2}-\frac54\Big)+\frac a2\mathbf e_\xi\cdot\frac{d\mathbf u_s}{dt}$ ·
      *why:* 9/8 = ½ × 9/4 and 5/8 = ½ × 5/4. This is (6.105). · *plain:* Surface pressure = a speed part + an
      acceleration part. · *set:* {dus: 2} · *watch:* "the orange curve appears".
  14. *did:* Set the acceleration to zero · *tex:* $\Big(\frac{p_a-p_\infty}{\frac12\rho\lvert\mathbf u_s\rvert^2}\Big)_{\text{steady}}=\frac94\cos^2\theta_s-\frac54=1-\frac94\sin^2\theta_s$ ·
      *why:* u_s·e_ξ = \lvert u_s\rvert cos θ_s, θ_s the angle from the velocity; cos² = 1 − sin². This is (6.106) = (6.91). ·
      *plain:* A steadily moving sphere feels exactly the pressures of a still sphere in a stream.
- **Result.** $\frac{p_a-p_\infty}{\rho}=\frac12\lvert\mathbf u_s\rvert^2\Big(\frac94\frac{(\mathbf u_s\cdot\mathbf e_\xi)^2}{\lvert\mathbf u_s\rvert^2}-\frac54\Big)+\frac a2\mathbf e_\xi\cdot\frac{d\mathbf u_s}{dt}$ *(6.105)* — *in
  words:* the pressure on a moving sphere is a fore–aft-symmetric part set by its speed plus a part set by its
  acceleration, high on the side it accelerates toward.
- **Check.** Units: \lvert u_s\rvert² and a du_s/dt both m²/s² ✓. Numbers (C15 code): u_s = 1 m/s, du_s/dt = 2 m/s², a = 0.1
  m, water, at the front point: steady 500 Pa + acceleration 100 Pa = 600 Pa ✓. Parity with `ch04.accelerating_sphere_pressure`
  for u_s ∥ du_s/dt; steady limit = (6.91) with θ_s = π − θ ✓.
- **sympy check intent (★★★, `check_src`):** re-run steps 2–13 for a concrete motion. `t, a, x, y, z = sp.symbols('t a x y
  z', real=True)` · `xs = sp.Matrix([sp.Rational(1, 2)*t**2, 0, t])` (a motion with a turning acceleration) · `us =
  xs.diff(t)` · `X = sp.Matrix([x, y, z])`; `xi = X - xs`; `r = sp.sqrt(xi.dot(xi))` · `phi = -a**3/(2*r**3)*us.dot(xi)`
  ((6.97)) · `dphidt_direct = sp.diff(phi, t)` (∂φ/∂t at fixed x, no chain rule) · `grad = sp.Matrix([sp.diff(phi, v) for v
  in (x, y, z)])` ((6.103)) · `chain = -grad.dot(us) - a**3/(2*r**3)*xi.dot(us.diff(t))` ((6.101)) · `print(sp.simplify(
  (dphidt_direct - chain).subs({x: 0.3, y: 0.4, z: 1.2, t: 0.7, a: 0.5})))` (→ 0 at a sample point) · on the surface: pick
  e = (0.6, 0, 0.8), substitute X = xs + a e, compare `grad` with `sp.Rational(3, 2)*us.dot(e)*e - us/2` (→ 0: (6.104)) and
  with the printed bracket `-(a**3/2)*(-3*e*us.dot(e)/a**3 - us/a**3)` (→ (3/2)(u·e)e + ½u ≠: the printed middle term is
  wrong) · finally `p = -(dphidt_direct) - grad.dot(grad)/2` on the surface vs (6.105) (→ 0). `ch06.moving_sphere_dphidt_sym()`
  wraps the same. Every line commented.
- **What it means.** The speed part cannot push the sphere (it is symmetric, D30 step 2); the acceleration part does —
  the origin of added mass. For a sphere accelerated from rest the pressure is purely the acceleration part.
- **Traps.** ∂/∂t at fixed x acts through x_s(t) and u_s(t). ∂φ/∂x_s = −∇φ, not +. The book's middle bracket of (6.104)
  has the wrong sign on u_s/a³. Expand \lvert(3/2)(u·e)e − ½u\rvert² with the cross term.

### D30 · The force on an accelerating sphere and its added mass: (6.98) → (6.108) → (6.109) — ★★, 9 steps, in C15 (notebook · `added_mass_sphere`)
- **Goal.** Integrate the surface pressure (6.105) over the sphere to find the fluid's force, and write Newton's second law
  for a submerged sphere.
- **Start.** $\mathbf F_s=-\int_{\text{sphere's surface}}(p-p_\infty)\mathbf n\,dA$ *(6.98)* with *(6.105)*, n = e_ξ — *in words:* the force is the
  pressure excess pushed along the inward normal.
- **Plan.** (1) Split the pressure into its two parts. (2) The speed part gives nothing. (3) Orient the axes and integrate
  the acceleration part (6.107)–(6.108). (4) Newton (6.109).
- **Tools.** surface integrals on a sphere (primer, in C15) · spherical unit vectors (ch03 P88, reminder) · symmetry
  (even integrands).
- **Assumptions.** As D29; the acceleration direction fixed at the instant considered.
- **Steps.**
  1. *did:* Split the pressure excess · *tex:* $p-p_\infty=S(\mathbf e_\xi)+\frac{\rho a}{2}\mathbf e_\xi\cdot\frac{d\mathbf u_s}{dt},\quad S=\frac{\rho}{2}\lvert\mathbf u_s\rvert^2\Big(\frac94\frac{(\mathbf u_s\cdot\mathbf e_\xi)^2}{\lvert\mathbf u_s\rvert^2}-\frac54\Big)$ ·
     *why:* (6.105) times ρ; the integral of a sum is the sum of integrals. · *plain:* A speed part and an acceleration
     part, integrated separately.
  2. *did:* Show the speed part gives no force · *tex:* $S(-\mathbf e_\xi)=S(\mathbf e_\xi)\ \Rightarrow\ \int S(\mathbf e_\xi)\,\mathbf e_\xi\,dA=0$ · *why:* S depends on
     (u_s·e)², unchanged when e → −e; opposite surface points carry equal pressures on opposite normals, which cancel. The
     book says "no drag" without integrating. · *plain:* Steady motion: no force (d'Alembert in 3-D). · *set:* {dus: 0}
     · *watch:* "the steady force bar stays 0".
  3. *did:* Point the acceleration along e_z · *tex:* $\frac{d\mathbf u_s}{dt}=\Big\lvert\frac{d\mathbf u_s}{dt}\Big\rvert\mathbf e_z\ \Rightarrow\ \mathbf e_\xi\cdot\frac{d\mathbf u_s}{dt}=\Big\lvert\frac{d\mathbf u_s}{dt}\Big\rvert\cos\theta$ ·
     *why:* We may choose axes at the instant of interest; θ is then the polar angle from the acceleration. · *plain:* The
     acceleration part is largest at the front, zero at the equator, negative behind.
  4. *did:* Write e_ξ and dA in spherical coordinates · *tex:* $\mathbf e_\xi=(\sin\theta\cos\varphi,\ \sin\theta\sin\varphi,\ \cos\theta),\quad dA=a^2\sin\theta\,d\theta\,d\varphi$ ·
     *why:* The outward unit normal and area element of a sphere of radius a (primer). · *plain:* Everything is now a
     function of the two angles.
  5. *did:* Substitute into (6.98) · *tex:* $\mathbf F_s=-\rho\frac a2\Big\lvert\frac{d\mathbf u_s}{dt}\Big\rvert\int_0^\pi\int_0^{2\pi}\cos\theta\,\mathbf e_\xi\,a^2\sin\theta\,d\varphi\,d\theta$ · *why:* Only the
     acceleration part remains (step 2). This is (6.107). · *plain:* A double integral over the sphere.
  6. *did:* Integrate over φ first · *tex:* $\int_0^{2\pi}\cos\varphi\,d\varphi=\int_0^{2\pi}\sin\varphi\,d\varphi=0,\qquad\int_0^{2\pi}d\varphi=2\pi$ · *why:* Over a full turn cos φ and
     sin φ cancel (primer, C06): the x and y components vanish; the z component gains 2π. · *plain:* The force points
     along the acceleration.
  7. *did:* Keep the z-component · *tex:* $\mathbf F_s=-\pi\rho a^3\Big\lvert\frac{d\mathbf u_s}{dt}\Big\rvert\mathbf e_z\int_0^\pi\cos^2\theta\sin\theta\,d\theta$ · *why:* ρ(a/2) × a² × 2π
     = πρa³. ⚠️ The book's (6.108) still writes dφ inside after this integration. · *plain:* One integral left.
  8. *did:* Integrate over θ · *tex:* $\int_0^\pi\cos^2\theta\sin\theta\,d\theta=\Big[-\tfrac13\cos^3\theta\Big]_0^\pi=\tfrac23\ \Rightarrow\ \mathbf F_s=-\frac23\pi\rho a^3\frac{d\mathbf u_s}{dt}\equiv-M\frac{d\mathbf u_s}{dt}$ ·
     *why:* Substitute s = cos θ (P106). This is (6.108), with M = 2πρa³/3. · *plain:* The fluid pushes back with half the
     displaced fluid's mass times the acceleration. · *live:* "F = −(2/3)π × 1000 × a³ × dus = F N".
  9. *did:* Write Newton's law for the sphere · *tex:* $\mathbf F_E+\mathbf F_s=m\frac{d\mathbf u_s}{dt}\ \Rightarrow\ \mathbf F_E=\Big(m+\frac{2\pi}{3}\rho a^3\Big)\frac{d\mathbf u_s}{dt}$ ·
     *why:* The sphere of mass m feels the external force F_E and the fluid force F_s; move F_s to the right. This is
     (6.109). · *plain:* Under water a sphere behaves as if it were heavier by half the displaced mass.
- **Result.** $\mathbf F_s=-M\frac{d\mathbf u_s}{dt}$, $M=\frac{2\pi a^3\rho}{3}$ *(6.108)* and $\mathbf F_E=\big(m+\frac{2\pi}{3}\rho a^3\big)\frac{d\mathbf u_s}{dt}$ *(6.109)* — *in
  words:* an ideal fluid exerts no force on a sphere moving steadily, but opposes its acceleration as an extra mass of
  half the displaced fluid.
- **Check.** Units: ρa³ du/dt [kg × m/s²] = N ✓. Numbers: a = 0.1 m, water → M = 2.094 kg (displaced 4.189 kg); du_s/dt = 2
  m/s² → F_s = −4.189 N; the Gauss–Legendre × trapezoid surface sum agrees (C15 code) ✓. Bubble (m → 0) with buoyancy
  ρVg: du/dt = ρVg/M = 2g ✓.
- **What it means.** Light bodies (bubbles, balloons, fish) are dominated by their added mass; heavy ones (a steel ball:
  +6 %) hardly notice. For other shapes M becomes a tensor M_ij.
- **Traps.** The steady part is not zero pointwise — it integrates to zero. The stray dφ in the printed (6.108). M is half
  the displaced mass for a sphere; for a cylinder (per unit length) it is the whole displaced mass.

### D31 · Added mass by kinetic energy: ½ρ∫\lvert∇φ\rvert² dV = ½MU² — ★★★, 10 steps, in C15 (notebook · `added_mass_sphere`)
- **Goal.** Find the kinetic energy of the fluid set moving by a sphere at speed U, and show it equals ½MU² with the same
  M = 2πρa³/3 as the force route — an independent derivation the book leaves to Exercise 6.49.
- **Start.** $T=\frac12\rho\int_{r>a}\lvert\nabla\phi\rvert^2dV$ with $\phi=-\frac{a^3U\cos\theta}{2r^2}$ (6.97 with u_s = U e_z, the sphere at the origin) — *in
  words:* the fluid's kinetic energy outside the sphere.
- **Plan.** (1) Turn \lvert∇φ\rvert² into a divergence (Green's first identity). (2) Gauss: volume → surfaces. (3) The far
  surface vanishes. (4) Evaluate on the sphere. (5) Read off M.
- **Tools.** Green's first identity (primer, in C15) · divergence theorem (Ch. 2 §2.12) · surface integrals on a sphere
  (primer, in C15) · kinetic energy and added mass (gloss).
- **Assumptions.** ∇²φ = 0 in the fluid (step 2); fluid at rest far away; the instant when the sphere is at the origin.
- **Steps.**
  1. *did:* Apply Green's first identity · *tex:* $\nabla\cdot(\phi\nabla\phi)=\lvert\nabla\phi\rvert^2+\phi\nabla^2\phi$ · *why:* The product rule for a
     divergence (P113) with the vector φ∇φ (primer). · *plain:* The energy density is almost a divergence.
  2. *did:* Use Laplace's equation · *tex:* $\lvert\nabla\phi\rvert^2=\nabla\cdot(\phi\nabla\phi)$ · *why:* ∇²φ = 0 for a potential flow of constant
     density (C13, (6.80)). · *plain:* Exactly a divergence.
  3. *did:* Apply the divergence theorem on the fluid · *tex:* $\int_{a<r<R}\lvert\nabla\phi\rvert^2dV=\oint_{r=R}\phi\frac{\partial\phi}{\partial r}dA-\oint_{r=a}\phi\frac{\partial\phi}{\partial r}dA$ ·
     *why:* Gauss on the shell between the sphere and a big sphere R; on the inner surface the normal pointing out of the
     fluid is −e_r, hence the minus. · *plain:* The fluid's energy is set by what happens on its two boundaries.
  4. *did:* Let the outer sphere grow · *tex:* $\phi\frac{\partial\phi}{\partial r}=O(R^{-5}),\ \ dA=O(R^2)\ \Rightarrow\ \oint_{r=R}\to0$ · *why:* φ ~ r⁻², ∂φ/∂r ~
     r⁻³; the area grows only as R², so the integral falls as R⁻³. · *plain:* The far fluid barely moves and carries no
     energy.
  5. *did:* Evaluate φ and ∂φ/∂r on r = a · *tex:* $\phi(a,\theta)=-\tfrac12Ua\cos\theta,\qquad\frac{\partial\phi}{\partial r}\Big\rvert_a=\frac{a^3U\cos\theta}{r^3}\Big\rvert_a=U\cos\theta$ · *why:*
     Substitute r = a; the second is the kinematic condition u_r = U cos θ (D28 step 5). · *plain:* On the surface the
     fluid moves out at the sphere's normal speed.
  6. *did:* Insert into the inner surface term · *tex:* $-\oint_{r=a}\phi\frac{\partial\phi}{\partial r}dA=\frac12U^2a\int_0^{2\pi}\!\!\int_0^\pi\cos^2\theta\,a^2\sin\theta\,d\theta\,d\varphi$ ·
     *why:* The product is −½U²a cos²θ; the two minus signs cancel; dA = a² sin θ dθ dφ (primer). · *plain:* One surface
     integral left.
  7. *did:* Integrate · *tex:* $\frac12U^2a^3\cdot2\pi\cdot\frac23=\frac{2\pi}{3}a^3U^2$ · *why:* ∫dφ = 2π and ∫cos²θ sin θ dθ = 2/3 (primer). ·
     *plain:* The volume integral of \lvert∇φ\rvert² equals (2π/3)a³U². · *live:* "(2π/3)a³U² = T2".
  8. *did:* Multiply by ½ρ · *tex:* $T=\frac12\rho\cdot\frac{2\pi}{3}a^3U^2$ · *why:* The kinetic energy per unit volume is ½ρ\lvert∇φ\rvert², and step 7 gave the volume integral. · *plain:* The fluid's kinetic
     energy.
  9. *did:* Define M by T = ½MU² · *tex:* $\tfrac12MU^2=\tfrac12\cdot\frac{2\pi}{3}\rho a^3U^2\ \Rightarrow\ M=\frac{2\pi}{3}\rho a^3$ · *why:* The fluid's energy
     looks like that of an extra mass M moving with the sphere (gloss). · *plain:* The added mass, from energy alone.
  10. *did:* Compare with the force route · *tex:* $M_{\text{energy}}=M_{\text{force}}=\frac{2\pi a^3\rho}{3}\quad\text{(6.108)}$ · *why:* Pushing the sphere
      at dU/dt, the power −F_s·U must feed dT/dt = MU dU/dt; so F_s = −M dU/dt: the two routes must agree, and they do. ·
      *plain:* Two independent routes, one number.
- **Result.** $T=\frac12\rho\int\lvert\nabla\phi\rvert^2dV=\frac12\Big(\frac{2\pi}{3}\rho a^3\Big)U^2$, so $M=\frac{2\pi a^3\rho}{3}$ — *in words:* the energy the sphere
  must give the fluid to move at U is that of half the displaced mass moving at U.
- **Check.** Units: ρa³U² [J] ✓. Numbers: a = 0.1 m, water, U = 1 m/s → T = 1.047 J, M = 2.094 kg = the D30 value;
  `added_mass_by_energy` by the surface form and by a `dblquad` volume integral agree to 1e-10 (C15 code) ✓.
- **sympy check intent (★★★, `check_src`):** `r, th, ph, a, U, rho = sp.symbols('r theta phi a U rho', positive=True)` ·
  `phi = -a**3*U*sp.cos(th)/(2*r**2)` (the potential) · `gr = sp.diff(phi, r)`; `gt = sp.diff(phi, th)/r` (spherical
  gradient components) · `vol = sp.integrate(sp.integrate((gr**2 + gt**2)*r**2*sp.sin(th), (th, 0, sp.pi)), (r, a,
  sp.oo))*2*sp.pi` (step 3's left side, the volume integral to infinity) · `surf = -sp.integrate((phi*gr).subs(r,
  a)*a**2*sp.sin(th), (th, 0, sp.pi))*2*sp.pi` (the inner surface form, steps 5–7) · `print(sp.simplify(vol - surf),
  sp.simplify(vol))` (→ 0 and 2*pi*U**2*a**3/3) · `M = sp.simplify(rho*vol/U**2)`; `print(M)` (→ 2*pi*a**3*rho/3) ·
  `print(sp.simplify(sp.diff(r**2*sp.diff(phi, r), r)/r**2 + sp.diff(sp.sin(th)*sp.diff(phi, th), th)/(r**2*sp.sin(th))))`
  (→ 0: the ∇²φ = 0 used in step 2). Every line commented.
- **What it means.** Added mass is the fluid's kinetic energy in disguise: whatever accelerates a body must also fill the
  fluid round it with energy ½MU². The energy route generalises to any shape (M becomes a tensor, T = ½M_ij U_i U_j).
- **Traps.** ∇·(φ∇φ) = \lvert∇φ\rvert² needs ∇²φ = 0. The normal out of the fluid points into the sphere, which fixes the sign.
  The far surface vanishes only because φ∂φ/∂r ~ r⁻⁵ beats the r² growth of the area.

## Part G — errata from verification (2026-09-23; these override Parts A, B, F where they conflict)
- **G1 · D21 sympy check (O1).** The principal-root failure fraction is 0.998, not 1.0: the three grid points on the slit (y = 0, |x| ≤ 2b) have |ζ| = b for both roots. Say "fails at every left-half point off the slit" or exclude the slit in the cell.
- **G2 · C12 / N74 refinement order (O2).** Example 6.2 converges at order **4/3 everywhere**, not "≈ 2 away from the corner": the 270° corner flow behaves as ψ ~ r^{2/3} sin(2θ/3) (n = ⅔ in $w = A z^n$ (6.46)), and this singularity pollutes the whole grid. The notebook's refinement cell must state 4/3 and this reason; second order holds for the smooth test problem (sin πx sinh πy: 1.995).
- **G3 · E8 text and ranges (O3).** (a) Rankine oval: the axial strengths k_n *alternate in sign* (±9.6 at N = 20, ±18.6 at N = 40); only their moments converge (dipole error 1.4e-3 → 2.9e-6). Do not say "the segments converge to two narrow spikes". (b) Airship: the fitted bars do *not* reproduce a point source plus a flat line sink (k from −4.9 to 14.1 at N = 20, Σk_nΔξ = −0.094). (c) Panels: "C_p converges like 1/N²" holds on the **ellipse**; on the circle constant-source panels are exact at every N. (d) Cap N at 40 in axial mode (the Rankine oval breaks down by N = 80, cond ≈ 3.5e17) or show a conditioning warning; the ellipsoid target is the smooth, convergent one (k_n order ≈ 1).
- **G4 · Panel velocity off the body** converges at about first order (0.93 → 1), while control-point C_p is exact (circle) or second order (ellipse).
- **G5 · Fig. 6.10's measured curve** has no fetched dataset: the separated band (`separated_cp_band`) is a labelled *qualitative* sketch.
