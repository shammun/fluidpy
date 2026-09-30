# Chapter 10 — Computational Fluid Dynamics: lesson design
(from `analysis/ch10_curation.md` (A 15 · B 105 · C 23 = 143 rows; CORE 15 · NOTE 113 · RECAP 13 · SKIP 2; 23 derivations
D01–D23 (★ 4 · ★★ 18 · ★★★ 1); E1–E8 + backup B1; primers from P221) and `analysis/ch10.md` (§2b derivations a-d1…a-d48, §4
implementation rows 1–47, §5 core modules `core/fd.py`, `core/fem1d.py`, `core/mac.py`, `core/maccormack.py`, `core/fem2d.py`, §9
slips R1–R12 and symbol collisions). Every equation placed here was re-read on the rendered pages `chapters/pages/ch10/`
(printed = pdf − 27): (10.4)–(10.7) p451, (10.8)–(10.12) p452, (10.13)–(10.17) p453, (10.18)–(10.25) p454, (10.26)–(10.31)
p455, (10.32)–(10.36) p456, (10.37)–(10.41) p457, (10.42)–(10.49) p458, (10.50)–(10.58) p459, (10.59)–(10.63) + Fig. 10.2 p460,
(10.64)–(10.66) + Fig. 10.3 p461, (10.67)–(10.77) p462, (10.82)–(10.86) p464, (10.87)–(10.92) p465, (10.93)–(10.94) p466,
(10.95)–(10.102) p467, (10.110)–(10.114) p469, (10.115)–(10.120) p470, Fig. 10.4 + (10.121)–(10.125) p471, (10.126) + the
no-pressure-BC paragraph p472, (10.134)–(10.136) + Fig. 10.5 p474, (10.137) p475; the rest from the analyst's page-image
transcription (analysis §2). 2026-09-30, lesson-designer. The implementer works in parallel from analysis §4 + curation §8–§9;
**Part C is written first and is the contract both sides keep.** Numbers marked *expect* were computed for this design with a
scratch script (numpy); the builder compares an executed cell against them.)

**Binding conventions for every builder (analysis §9, curation decisions 1–7 and §8).**
1. **Imports and aliases.** `from fluidpy import ch10_computational_fluid_dynamics as ch10`; `from fluidpy.core import fd as
   FD, fem1d as FEM1, mac as MAC, maccormack as MCK, fem2d as FEM2`. **`ch10` re-exports every public name of `FD`, `FEM1`,
   `MAC`, `MCK`, `FEM2`**, so explainer parity rows write `ch10.<name>` only. Dimensional functions (§10.2–10.3) take SI:
   x, L, Δx [m]; t, Δt [s]; u [m/s]; D [m²/s]; T in any unit (a temperature, a concentration). Everything in §10.4–10.5 is
   **dimensionless** (lengths in the cavity side / block side / cylinder diameter, speeds in the lid / inflow speed, Re, Ma).
2. **Scalar-callable and parity-friendly** (the shot.py evaluator has no builtins, only `np`, `math`, `ch10`): every public
   function accepts Python floats and returns a float, a tuple or a `dict` of floats (arrays only when an array goes in or
   when the result is a field). Complex results have a real-valued twin (`amplification_modulus`). Parity `py:` expressions
   use `ch10.…`, `np.…`, numbers, strings, lists and keywords, indexed down to one number.
3. **Symbols (overloading resolved; a ⚠️ table cell in the front matter, and a one-line reminder where each is first used).**

| Book symbol | Meanings in the book | Notebook / explainer symbol | Code name |
|---|---|---|---|
| i | grid index (10.4)–(10.30) and √−1 (10.20)–(10.24) | index $i$ (italic) in quoted book lines; **our own lines in D04–D08 rename the index $j$**; imaginary unit always upright $\mathrm i$ | `i`, `j`; `1j` |
| α, β | FTCS numbers (10.11); Θ-scheme weights (10.129)–(10.133); time weights (10.163) | FTCS: $\alpha=u\Delta t/(2\Delta x)$, $\beta=D\Delta t/\Delta x^2$ (the book's); Courant number $C=2\alpha=u\Delta t/\Delta x$; Θ-scheme $\alpha_\Theta,\beta_\Theta$; time weights $\alpha_t,\beta_t$ | `alpha, beta`, `C`; `alpha_split, beta_split`; `alpha_t, beta_t` |
| θ | Fourier angle θ = kπΔx (10.26); Θ-scheme fraction | Fourier angle $\theta$ (radians per cell); Θ-scheme fraction $\Theta$ (capital) | `theta`; `theta_split` |
| g, gⁿ(k) | Dirichlet value (10.2); Fourier amplitude (10.20); body force **g** (10.79); its y-component (10.119) | Dirichlet value $g$; Fourier amplitude written $\hat\xi^n$ in our lines; body force $\mathbf g=(g_x,g_y)$; gravity in the climate CFL note $g_0$ | `g`; `xi_hat`; `body` |
| D | diffusivity (10.1); strain-rate tensor D[u] (10.136); cavity side and block side (§10.5) | $D$ diffusivity; $\mathbf D[\mathbf u]$ bold with brackets; cavity side $L$; block side $d$ | `D`; `strain_rate`; `L`, `d` |
| M | mass matrix (10.55); Mach number (§10.5) | $\mathbf M$ (bold); $Ma=U/c$ | `M`; `Ma` |
| K | error constant (10.15); stiffness matrix (10.56) | error constant $K_e$; $\mathbf K$ (bold) | `K_err`; `K` |
| S | trial space (10.32); Strouhal number (p. 468) | trial space $\mathcal S$; Strouhal $St$ | —; `St` |
| n | time level; number of elements/cells; shedding frequency n in S = nd/U | superscript $n$ = time level; $n_{el}$ elements, $N$ cells; shedding frequency $f_s$ [cycles per unit time] | `n_el`, `n`; `f_s` |
| c | sound speed (10.95), (10.99); MacCormack coefficients c₁–c₅ (10.109); FE test coefficients c_A (10.45) | $c$ sound speed; $c_1\dots c_5$; test coefficients $c_A$ (capital-letter index) | `c`; `c1…c5`; `c_A` |
| E, e | truncation error Eᵢⁿ (10.16); flux vector **E** (10.100); solution error eᵢⁿ (10.14) | $E^n_i$ italic; $\mathbf E$ bold; $e^n_i$ | `E_trunc`; `E_flux`; `err` |
| G | amplification factor (our name for gⁿ⁺¹/gⁿ); Galerkin bracket G_A (10.52) | $G(\theta)$; $G_A$ | `G`; `G_A` |
| δ | layer thickness (10.89); Kronecker δ_AB (10.61) | $\delta$; $\delta_{AB}$ | `delta`; — |
| h | mesh size; superscript "discrete" (Sʰ, Tʰ) | $h$; superscript $h$ | `h` |
| R | global Péclet (10.87); R_cell (10.31) | $R=uL/D$; $R_{cell}=u\Delta x/D$ | `R`, `R_cell` |
| p′, u′ | Newton corrections (10.164) (ch04/ch07: perturbations) | $p'$, $\mathbf u'$ with the word "correction" | `dp`, `du` |
| a(·,·), a | bilinear form (10.42); (none) | $a(w,v)$; the linear-advection speed in D18 is $u$ (never $a$) | `bilinear_form` |

4. **Book slips taught in corrected form** (analysis §9 R1–R12 = curation §8): each is a printed-vs-correct box
   (`> ⚠️ **The book prints** … **; the correct form is** …`) where it is used, and, where computable, a planted wrong variant
   in code that a test must fail. **Two values, the book's first** (CUMULATIVE rule) where the book gives a number.

| Slip | The book prints | Correct | Taught in | Code option that must fail |
|---|---|---|---|---|
| R1 (10.67)–(10.68) | on the element $[x_{A-1},x_A]$: $\frac{dN_A}{dx}=\frac{2}{x_A-x_{A-1}}\frac{dN_1}{d\xi}=\frac{-1}{x_A-x_{A-1}}$ and $\frac{dN_{A+1}}{dx}=\frac{1}{x_A-x_{A-1}}$ | the local shapes there are $N_{A-1}\leftrightarrow N_1$ and $N_A\leftrightarrow N_2$: $\frac{dN_{A-1}}{dx}=-\frac1{h^e}$, $\frac{dN_A}{dx}=+\frac1{h^e}$ | C08 (N35, D13 step 7) | `FEM1.shape_slopes(xa, xb, printed=True)` fails the chain-rule test |
| R2 text after (10.73) | "the nonzero ones require that A = e or e + 1 and B = e or e + 1" | $A,B\in\{e-1,e\}$ (node 0 is the Dirichlet node, hence $f^1_a=-g\,k^1_{a1}$) | C08 (N37, D14 step 9) | `FEM1.connectivity(n_el, printed=True)` gives an out-of-range node |
| R3 text after (10.93) | "a forward-difference scheme is used to discretize the convective term" | $T_j-T_{j-1}$ is a **backward** (upwind for u > 0) difference | C09 (N48, D17 step 1) | ⚠️ box; `FD.steady_cd_fd(scheme="forward")` (downwind) wiggles for every R_cell |
| R4 no-pressure-BC paragraph | "$p_{0,2}$ will not appear in equation (10.120)" | "…in equation (10.124)" — (10.120) is the v-predictor | C12 (N70, D20 step 9) | ⚠️ box |
| R5 cavity algorithm Step 5 | $2\rho^{n+1}=(\rho^n+\rho^*)-a_1+[(\rho u)^*_{i,j}-(\rho u)^*_{i-1,j}]-\dots$ | $-a_1[(\rho u)^*_{i,j}-(\rho u)^*_{i-1,j}]$ (compare (10.106)) | C14 (N85) | `MCK.cavity_coefficients(..., printed=True)` path in `weakly_compressible_step` breaks mass conservation |
| R6 (10.172) | second sum $\sum_{A'}u_{A'}\int N^u_{A',y}N^p_B\,d\Omega$ | $\sum_{A'}v_{A'}\int N^u_{A',y}N^p_B\,d\Omega$ (from $\partial v'/\partial y$; pairs with $\mathbf B^T_{vp}$) | C13 (N101) | `FEM2.assemble_newton_system(..., printed_10_172=True)` fails the Stokes/Poiseuille test |
| R7 (10.186) | third expansion "$v'=\sum p^e_b\psi_b$" | $p'=\sum_{b=1}^3p^e_b\psi_b$ | C13 (N106) | ⚠️ box |
| R8 p. 458, Fig. 10.14 caption | "the fourth example later in this section" | §10.5 has three examples; the cylinder is the third | C15 (N92) | ⚠️ box |
| R9 p. 470 | computed confined St compared with an unbounded value | the computed cylinder is confined (a channel a few diameters wide, sliding walls): trend only; the unbounded literature is 0.164–0.165 (Williamson 1996 exp.; Kravchenko; Posdziech & Grundmann) and the book's quoted value is a little higher — **never compare the confined C_D or St numerically with unbounded data** | C13 (R13) | ⚠️ box; no numeric comparison in any cell |
| R10 (10.166) | $\beta\,\partial v^*/\partial t(t_n)$ | $\beta\,\partial v/\partial t(t_n)$ — known data at $t_n$, not a Newton iterate | C13 (N99) | ⚠️ box |
| R11 (10.29)–(10.30) | $T^{n+1}_i=T^n_i-2\alpha(T^n_i-T^n_{i-1})$, $u\Delta t/\Delta x\le1$ with no sign of u | valid for u > 0; for u < 0 the upwind side is $i+1$; the CFL condition is $\lvert u\rvert\Delta t/\Delta x\le1$ | C05 (N16, D08 step 9) | `FD.transport_1d_step(..., scheme="upwind_printed")` blows up for u < 0 |
| R12 (10.155) | $\Delta t\le\frac{\sigma}{\sqrt2}Ma\,\Delta x$ follows "with large grid Reynolds numbers" | it also needs $Ma\ll1$ ($\lvert u\rvert/\Delta x+\lvert v\rvert/\Delta y\ll\sqrt2\,c/\Delta x$) | C15 (N91) | `MCK.maccormack_dt_asymptotic` vs `MCK.maccormack_dt` differ by > 10 % at Ma = 0.5 |

5. **Colours (text, figures, derivation terms, explainers — one meaning each; curation §5):** exact / analytic solution
   `muted` dashed ghost · FTCS and centred schemes `orange` · upwind `teal` · implicit (BTCS, Crank–Nicolson, backward Euler)
   `blue` · MacCormack / Lax–Wendroff `accent` purple · unstable growth, wiggles, divergence `rose` · physical diffusivity
   `blue`, numerical diffusivity `rose` · pressure and the pressure correction −Δt∇p `orange` · velocity arrows `blue` ·
   operator A₁ (convection–diffusion) `teal`, A₂ (pressure/continuity) `orange` · Taylor terms: the kept derivative `blue`,
   cancelled terms `muted`, the leading error term `rose` · FE: hat functions `muted`, weighted hats `teal`, their sum Tʰ
   `blue`, the active element `amber` · benchmark points (Ghia, Hou) black dots with a white halo · printed-slip ghosts
   `muted` dotted, labelled "as printed".
6. **Every book equation is shown in full next to its number** (CLAUDE.md rule 3) — markdown, derivation steps
   ("substitute (10.4), $T^n_{i+1}=T^n_i+\Delta x\,T_x+\dots$"), traps, recaps, notes, figure notes and every explainer text.
   Builders reuse the ch07 `show_eqs(text, EQ)` helper with an `EQ` dict holding all 199 labels of ch10 (the TeX of analysis
   §2 column 5, checked against the pages listed above) plus the earlier labels cited here ((1.x) FTCS of ch01 is quoted by
   formula, (3.5), (4.7), (4.39b), (4.102), (4.107), (4.108), (7.28)); JS explainers use a local `showEqs`. `tools/eq_refs.py`
   must list 0 offenders. Exercises are cited as "Exercise 10.3", never as bare numbers; **(10.199)** sits in Exercise 10.3:
   the equation itself is shown (N113) with **our** wall temperature and grid — no exercise text, inputs or answers.
7. **nbkit behaviour** (ch04–ch09 lesson): `nb.recap(...)` and `nb.section(...)` end the current CORE block, so every RECAP
   sits **before** the `nb.core` call of the block that uses it (R01–R03 before C01, R04–R05 before C04, R06–R07 before C09,
   R09 before C10, R08 before C11, R10 before C13, R11 before C15), except R12 and R13, which open the "C13 (continued): the
   cylinder by mixed finite elements" part of §10.5 (no block may be restarted; notes after a section break carry the
   parent id in bold, "**C13 (continued)**", as ch09 did for §9.9). Every `nb.derivation` sits inside its curation CORE
   block with `ref=` a bare label ("10.6").
8. **Parsing.** Part E is the only part with table rows after its heading; Part F has no line starting with a table bar
   (absolute values are written \lvert … \rvert or in words). Explainer headings are exactly `### E1 · fd_stencil_order` …
   `### E8 · mixed_fe_lbb` and `### B1 · operator_splitting_theta`. Primer terms in Part A are the exact Concept text of the
   Part E rows marked "primer" (coverage_check matches the first 18 characters). "Explained by" never names an earlier
   chapter's C-number (it cites "Ch. 4 §4.x" or a primers.md entry).
9. **Public repo (rule 9).** No exercise text or answers; the Moore's-law time, the book's rounded percentages for e⁻¹ and e⁻², the book's safety factor σ,
   the printed decimal value of Glowinski's θ, the §10.5 geometries, Mach numbers, grids, C_D/C_L/torque values, periods, St and mesh counts stay in
   `tests/book_values_ch10.json`. The notebook shows **our** numbers from our inputs: e^{−1} and e^{−2} are computed; σ is an argument with
   **our** default 0.8 (labelled "safety factor, ours"); the cavity runs at Ma = 0.08 (ours)
   on our grids (16² … 128²); the block (our geometry: side d = 1, channel height 4, 8 sides upstream, 20 downstream) runs at
   Ma = 0.06 on Δx = 1/8 (coarse) and cached 1/16, 1/32; the cylinder figures come from our
   mesh (d = 1, channel width W = 5, 8d upstream, 16d downstream, the confined sliding-wall set-up) with our force history in `reference/ch10/cylinder_fe_re100_forces.csv`.
   **Ghia (1982) Table I and the Hou et al. (1995) centres are public published data**, loaded from `reference/ch10/`
   (never typed into prose), cited to their authors. Figures are drawn by our code (every "Fig. 10.n" is remade, none copied).
10. **Honesty labels.** (10.110) is semi-empirical ("Tannehill, Anderson & Pletcher 1997, as cited"); (10.27) is Noye's result
    — derived here in D06; (10.127)–(10.128) are Peyret & Taylor's — stated and checked numerically, not derived; the
    Θ-scheme order is Glowinski's — checked by a sympy series cell; the block and cylinder wall-density closures and the
    corner averaging are **heuristic**; "fully implicit and unconditionally stable" (Newton–FE) holds for linear model
    problems, not as a theorem for Navier–Stokes; the confined C_D and St are *qualitative* in every comparison.
11. **Verifier notes found while designing (must reach the verifier).** (i) The upwind steady solution does **not** equal the
    exact solution of the modified equation (10.94) at finite R_cell: with n = 4, R_cell = 1 the upwind nodes are 0.0667,
    0.200, 0.467 and the modified-equation solution 0.0708, 0.209, 0.477 (agreement is O(Δx), improving with n at fixed R) —
    the notebook shows the gap shrinking, and `test_ch10` must assert convergence of the gap, not equality. (ii) The
    Lax–Wendroff/MacCormack local error is $-\frac{u\Delta t\Delta x^2}{6}(1-C^2)T_{xxx}$ (D18 step 12): it vanishes at C = 1
    (exact shift) — a V1 test. (iii) The FTCS zigzag growth factor at β = 0.51 is ∣1 − 4β∣ = 1.04: a 10⁻¹⁶ round-off seed needs
    ≈ 940 steps to reach O(1), a 10⁻¹⁰ kick ≈ 470 steps to reach 10⁻²; the Lax demo runs 2000 steps. (iv) For P1–P1 on the
    structured diagonal mesh with all-Dirichlet velocity, N_u = 2(n − 1)² and N_p = (n + 1)² − 1: for n ≤ 4 there are more
    pressure than velocity unknowns, so β_h = 0 exactly (count-forced spurious modes) — `FEM2.infsup_constant` reports both
    β_h and the number of spurious modes.

Order of parts: C (contract) · A (notebook storyboard) · B (explainers) · D (runtime) · E (prerequisite ledger) · F
(derivations).

---
## Part C — functions the builders will call (the implementer's contract)

**Status column.** **A#** = planned in `analysis/ch10.md` §4 row # (signature kept or refined here). **§8** = added by the
curation's §8. **NEW** = added by this design (in neither) — flagged as the phase asks. **CHG** = analysis signature changed
here (reason given). **reuse** = exists (ch01–ch09) and is only called. Docstrings cite § and Eq. (from the rendered page),
list symbols with units and assumptions, and carry the validation label. Return types: `float`/`ndarray` for single
quantities, `dict` (keys listed) for several named results. Every function is scalar-callable unless "arrays" is stated.
Grid layout `[j, i]` (y then x) as in `core/grids.py`.

### C.0 Reused (existing; called, not changed)
| # | Callable (signature) | Returns | Used by | Status |
|---|---|---|---|---|
| 0.1 | `style.setup_notebook() -> bool` (FAST) · `style.savefig(fig, "ch10", name)` · `style.COLORS` | FAST · path · palette | setup, every figure | reuse |
| 0.2 | `anim.animate(update, frames, fig, interval)` · `anim.show_animation(anim, player="video" or "frames")` | HTML | A1–A4 | reuse |
| 0.3 | `interact.slider_figure(fn, name, values, *, unit, xlabel, ylabel, title, xrange, yrange, height)` · `interact.animate_figure(frame_fn, times, …)` · `interact.live(fn, **widgets)` | plotly figure / widget | F1–F7, two live cells | reuse |
| 0.4 | `embed.show_viz("ch10", slug)` | display | 8 explainer cells | reuse |
| 0.5 | `tools.convergence.observed_order(h, err)` · `pairwise_orders(h, err)` · `richardson(coarse, fine, r=2, p=2)` · `refinement_study(solve, ns, h_of_n)` | float / list / float / Study | C01, C02, C03, C09, C14, C15 | reuse |
| 0.6 | `core.diffusion.stable_time_step(D, dy, safety=0.9)` · `ftcs_diffusion_1d(f0, D, dy, dt, nsteps, …)` · `crank_nicolson_1d(u0, y, dt, nsteps, D, bc_left, bc_right=0.0, startup_be=2, …)` · `couette_startup_profile(y, t, U, h, nu, nterms=200, tol=1e-12)` | float / arrays | C02 (u = 0 parity), C04 (R05), C05 (N113 pattern) | reuse |
| 0.7 | `core.laplace_solvers.solve_poisson(mask, f, bc, method="direct" or "sor", **kw)` · `sor_sweep(…)` · `assemble_laplace(mask, bc, dx, dy, f, order)` | arrays / csr | C12 (SOR option), C14 (ψ from ω) | reuse |
| 0.8 | `core.navier_stokes.exact_solution(name, x, t=0.0, **p)` (`"taylor_green"`, `"poiseuille"`, `"couette"`) · `navier_stokes_sym(rho, u_exprs, p_expr, coords, t, mu, mu_v=0, g=None, form="4.38")` · `conservative_to_advective_sym(rho_expr, u_exprs, coords, t)` | dict / sympy list | C09 (R06), C10 (R09), C11 (R08), C12 (Taylor–Green) | reuse |
| 0.9 | `core.similarity.reynolds_number(U, l, nu=None, rho=None, mu=None)` · `mach_number(U, c=None, T=None)` · `strouhal_number(Omega, l, U)` (Ω form; pass 2π f_s) · `drag_coefficient(F, rho, U, A)` | float | C10, C13 (R13), C15 (R11) | reuse |
| 0.10 | `ch04.is_incompressible_regime(U, T=288.15, threshold=0.3)` | bool | C10 (N53) | reuse |
| 0.11 | `core.waves.phase_speed(k, H=np.inf, g=G_BOOK)` (shallow limit k → 0 gives √(gH)) | float [m/s] | C05 (climate CFL note) | reuse |
| 0.12 | `core.bluff_body.cylinder_flow_regime(Re)` · `shedding_frequency(U, d, St=0.2)` | dict | C13 (R12, R13) | reuse |
| 0.13 | `ch07.kdv_solve(eta0, x, t_out, H, g, dt=None, dealias=True)` | dict | C15 (N112 pointer, one line, optional) | reuse |

### C.1 `fluidpy/core/fd.py` (NEW module, `FD`; re-exported by `ch10`) — stencils, schemes, stability, accuracy (Ch. 10–13, 15)
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 1.1 | `fd_weights(offsets, m, symbolic=True)` | Taylor matching (10.4)–(10.7): solve $\sum_j w_js_j^k/k!=\delta_{km}$, k = 0…len−1 (Vandermonde, exact `sympy.Rational`) | list of Rational (multiply by $1/h^m$); `symbolic=False` → ndarray — *expect* [−1, 0, 1], m = 1 → (−1/2, 0, 1/2); [−1, 0, 1], m = 2 → (1, −2, 1); [0, 1, 2], m = 1 → (−3/2, 2, −1/2); [0, −1, −2, −3], m = 2 → (2, −5, 4, −1) | C01, C14 (N80), C15 (N88), E1 | A#1 |
| 1.2 | `taylor_table_sympy(offsets, order)` | the table of $s_j^k/k!$ | sympy Matrix (rows = offsets, columns = k) | C01 (D01 check) | A#2 |
| 1.3 | `stencil_taylor_coefficients(kind, m=1, n_terms=6)` — kinds `"forward"`, `"backward"`, `"central"`, `"central2"` (second derivative), `"onesided2"` ([0, 1, 2]) | $\sum_jw_js_j^k/k!$ for k = 0…n_terms−1 (coefficient of $h^{k-m}f^{(k)}$) | dict(offsets, weights (floats), coeffs (list of floats), order p (first k > m with a nonzero coefficient, minus m), leading (that coefficient)) — *expect* forward: coeffs (0, 1, 0.5, 1/6, 1/24, 1/120), order 1, leading 0.5; central: (0, 1, 0, 1/6, 0, 1/120), order 2, leading 1/6; central2: (0, 0, 1, 0, 1/12, 0), order 2, leading 1/12 | C01, E1 (term bars) | **NEW** |
| 1.4 | `test_function(name, x, deriv=0)` — names `"sin"` (sin x), `"exp"` (eˣ), `"gauss"` ($e^{-x^2}$) | exact derivatives 0–5 | float / ndarray | C01, E1 | **NEW** |
| 1.5 | `stencil_derivative(func, x0, h, kind="central", m=1)` (`func` a callable or a `test_function` name) · `stencil_error(name, x0, h, kind="central", m=1)` (signed: stencil − exact) | (10.6), (10.7), one-sided second order (p. 450) | float — *expect* sin at x₀ = 1, h = 0.1: forward 0.497364 (error −4.294e-2), backward error +4.114e-2, central 0.539402 (−9.001e-4), central2 error +7.010e-4, onesided2 error +1.585e-3; h = 0.05: −2.126e-2, −2.251e-4, +1.753e-4 | C01, E1 | **NEW** |
| 1.6 | `fd_derivative(f, dx, kind="central", m=1, axis=-1)` (arrays) | explicit stencils, never `np.gradient` | ndarray, same shape; NaN where the stencil does not fit | C01, C02 | A#3 |
| 1.7 | `mixed_derivative_onesided(f, dx, dy, side)` (arrays) | (10.150) pattern | ndarray | C15 (N88) | A#4 |
| 1.8 | `transport_rhs(T, u, D, dx, periodic=True)` (arrays) · `advected_gaussian(x, t, u, D, x0=0.3, s0=0.05, L=1.0, n_images=3)` | semi-discrete RHS of (10.1); exact $T=\frac{s_0}{s(t)}\sum_m\exp\!\big(-\frac{(x-x_0-ut-mL)^2}{2s(t)^2}\big)$, $s^2=s_0^2+2Dt$ | ndarray | C02, C03, E1 (mode FTCS), E3 | A#5 |
| 1.9 | `ftcs_coefficients(u, D, dx, dt)` | (10.11) | tuple (alpha, beta) — *expect* (0.1, 0.01, 0.1, 0.2) → (0.1, 0.2) | C02, E2 | A#6 |
| 1.10 | `cfl_number(u, dt, dx)` (= ∣u∣Δt/Δx) · `diffusion_number(D, dt, dx)` (= β) · `cell_peclet(u, dx, D)` (= ∣u∣Δx/D) | (10.30), (10.11), (10.31) | float — *expect* cell_peclet(1, 0.01, 0.005) = 2 | C02, C05, C09 | A#14 |
| 1.11 | `transport_1d_step(T, alpha, beta, scheme="ftcs", periodic=True, g=None, q=None, dx=None)` — schemes `"ftcs"` (10.10), `"btcs"` (10.13), `"upwind"` (10.29) with the side from sign(α), `"upwind_printed"` (always $i-1$; R11 wrong variant), `"cn"`, `"lax_wendroff"`; non-periodic: left Dirichlet `g`, right Neumann `q` by a ghost node $T_{N+1}=T_{N-1}+2\Delta x\,q$ | (10.10), (10.13), (10.29) | ndarray (new level) — *expect* FTCS on [0, 0, 1, 0, 0] (periodic) with (α, β) = (0.1, 0.2) → [0, 0.1, 0.6, 0.3, 0] (sum 1); upwind with α = 0.25 (C = 0.5) → [0, 0, 0.5, 0.5, 0] | C02, C04, C05, E2, E3 | A#7, **CHG** (periodic default; `"upwind_printed"` added for R11) |
| 1.12 | `solve_transport_1d(T0, x, u, D, dt, nsteps, scheme="ftcs", g=None, q=None, periodic=False, save_every=None, check_stability=True, growth_cap=1e6)` (arrays) | loop of 1.11; raises `ValueError` before running when (10.27)/(10.30) is violated unless `check_stability=False`; stops when max∣T∣ > growth_cap × max∣T0∣ | dict(t, T (final), history (list), times, blew_up (bool), step_blown (int or None), max_abs (per saved step)) | C02, C05 (Lax demo), A1 | A#8, **CHG** (growth cap and dict return) |
| 1.13 | `error_norm(num, exact, kind="rms")` — kinds `"rms"`, `"max"`, `"l1"` | (10.14)–(10.15) | float | C03, C14, C15 | A#9 |
| 1.14 | `convergence_study(scheme="ftcs", n_list=(20, 40, 80, 160), u=0.5, D=0.01, t_end=0.2, dt_rule="diffusive", L=1.0, beta=0.25, C=0.5)` — `dt_rule` `"diffusive"` (Δt = βΔx²/D) or `"advective"` (Δt = CΔx/∣u∣) | (10.15) with the exact Gaussian | dict(dx (list), dt (list), err (list, rms), order (observed, in Δx), pairwise) — *expect* FTCS diffusive order 2.0 ± 0.15; upwind advective order 1.0 ± 0.15; BTCS diffusive 2.0 ± 0.15 | C03, C05 | A#9 |
| 1.15 | `truncation_error_sympy(scheme)` — `"ftcs"`, `"btcs"`, `"upwind"`, `"upwind_steady"`, `"central_steady"`, `"lax_wendroff"` | D03, D17, D18: substitute the two-variable Taylor series into the stencil, divide by the time-step factor, subtract the PDE | dict(E (sympy expression, leading terms), order_t, order_x, modified_equation (sympy Eq)) — *expect* FTCS: E = Δt/2 T_tt + uΔx²/6 T_xxx − DΔx²/12 T_xxxx; upwind_steady: extra diffusion uΔx/2 | C03, C09, C10 | A#10 |
| 1.16 | `truncation_terms(u, D, dx, dt, x, t, x0=0.3, s0=0.05, L=1.0)` | (10.17) evaluated on `advected_gaussian` (exact derivatives by the closed form) | dict(time = Δt/2·T_tt, conv = uΔx²/6·T_xxx, diff = −DΔx²/12·T_xxxx, total, measured = one-step FTCS residual) — scalar per (x, t) | C03, E1 (mode FTCS), F1 | §8.1 |
| 1.17 | `propagate_error(xi0, alpha, beta, nsteps, scheme="ftcs", periodic=True)` (arrays) | (10.19) = 1.11 applied to the disturbance | dict(xi (final), max_abs (per step, ndarray)) | C04, E2 | A#11 |
| 1.18 | `fourier_mode(x, k, convention="book")` | (10.21): $e^{\mathrm i\pi kx}$ (book, θ = kπΔx) or $e^{\mathrm ikx}$ (`"standard"`, θ = kΔx) | complex ndarray | C04 (R04) | A#11 |
| 1.19 | `amplification_factor(theta, alpha=0.0, beta=0.0, scheme="ftcs")` — FTCS (10.24) $G=(\alpha+\beta)e^{-\mathrm i\theta}+(1-2\beta)+(\beta-\alpha)e^{\mathrm i\theta}$; BTCS $G=[1+4\beta\sin^2\frac\theta2+2\mathrm i\alpha\sin\theta]^{-1}$; upwind $G=1-C(1-e^{-\mathrm i\theta})-4\beta\sin^2\frac\theta2$ with C = 2α (sign-aware: for α < 0 the mirror form); CN $G=\frac{1-2\beta\sin^2(\theta/2)-\mathrm i\alpha\sin\theta}{1+2\beta\sin^2(\theta/2)+\mathrm i\alpha\sin\theta}$; Lax–Wendroff $G=1-\mathrm iC\sin\theta-C^2(1-\cos\theta)$ | (10.24), D07, D08, D18 | complex (ndarray for array θ) — *expect* FTCS (θ = π/2, α = 0.1, β = 0.2) → 0.6 − 0.2i; BTCS (π, 0, 100) → 1/401 = 0.0024938 | C04, C05, C10, E2, E3 | A#12 |
| 1.20 | `amplification_modulus(theta, alpha=0.0, beta=0.0, scheme="ftcs")` · `ftcs_amplification_modulus2(theta, alpha, beta)` | ∣G∣ (float twin of 1.19) · (10.26) $\big(1-4\beta\sin^2\frac\theta2\big)^2+(2\alpha\sin\theta)^2$ | float — *expect* modulus2(π/2, 0.1, 0.2) = 0.40 | C04, E2, E3 (parity) | **NEW** (twin), A#12 |
| 1.21 | `amplification_curve(scheme, alpha, beta, n_theta=361)` | G on θ ∈ [0, π] | dict(theta, G (complex ndarray), modulus) | E2 table, F2 | §8.2 |
| 1.22 | `max_amplification(scheme, alpha, beta, n_theta=2001)` · `worst_theta(scheme, alpha, beta, n_theta=2001)` · `is_von_neumann_stable(scheme, alpha, beta, n_theta=2001, tol=1e-12)` · `ftcs_stable(alpha, beta)` (closed form (10.27)) | (10.25), (10.27) | float / float [rad] / bool / bool — *expect* max_amplification("ftcs", 0, 0.51) = 1.04 at θ = π; ftcs_stable(0.1, 0.2) True, (0.4, 0.2) False (4α² = 0.64 > 2β = 0.4), (0, 0.51) False | C04, E2 | A#13, **NEW** (`worst_theta`) |
| 1.23 | `stability_verdict(scheme, alpha, beta)` | Noye / CFL reading | dict(stable (bool), reason (exact ASCII strings the explainers copy: "stable", "unstable: 4a^2 > 2b, long waves grow", "unstable: 2b > 1, the zigzag grows", "unstable: pure convection, every wave grows", "unstable: C > 1, the characteristic leaves the stencil", "stable for every dt (implicit)"), Gmax, theta_worst) | C04, C05, E2, E3 (status) | **NEW** |
| 1.24 | `phase_error(scheme, C, theta)` | relative phase speed $\frac{-\arg G}{C\theta}$ | float — *expect* LW (C = 0.5, θ = π/2) 0.7487; upwind 1.0000 (C = 0.5); LW (0.8, π/4) 0.9679 | C10, E3 | §8.2 |
| 1.25 | `numerical_diffusivity(u, dx, D=0.0, scheme="upwind_steady", C=0.0)` | (10.94): steady upwind 0.5 R_cell D = ∣u∣Δx/2; unsteady upwind ∣u∣Δx(1 − C)/2 (ours, D08 note) | float [m²/s] — *expect* (1, 0.01, C = 0.5, `"upwind"`) = 0.0025 | C05, C09, E3, E4 | A#14, **CHG** (C and unsteady form) |
| 1.26 | `ftcs2d_max_amplification(u, v, Re, dx, dt, n=721)` | 2-D von Neumann of the explicit convection–diffusion step (the MAC predictor) | float | C12 (N73) | A#15 |
| 1.27 | `steady_cd_exact(x, R, L=1.0, limit=None)` (overflow-safe: $e^{-R(1-x/L)}\frac{1-e^{-Rx/L}}{1-e^{-R}}$ with `expm1` for R > 0; x/L for R → 0; `limit="large_R"` → (10.88)) · `cd_layer_thickness(R, L=1.0, level=np.exp(-1))` | (10.86)–(10.89) | ndarray / float [m] — *expect* R = 4 at x = 0.25, 0.5, 0.75: 0.03206, 0.1192, 0.3561; finite at R = 10⁴ | C09, E4 | A#16 |
| 1.28 | `steady_cd_fd(n, R, scheme="central", grid="uniform", L=1.0, beta_s=2.0)` — schemes `"central"` (10.91), `"upwind"` (10.93), `"forward"` (downwind; R3 illustration); grids `"uniform"`, `"stretched"` (ours: nodes clustered at x = L by a tanh map with strength `beta_s`, nonuniform 3-point stencils) | (10.90), (10.91), (10.93) | dict(x, T) (tridiagonal `solve_banded`) — *expect* n = 4, R = 4 central → [0, 0.025, 0.1, 0.325, 1]; R = 16 central → [0, −0.05, 0.1, −0.35, 1]; upwind → [0, 0.00641, 0.03846, 0.1987, 1] | C09, E4 | A#17 |
| 1.29 | `steady_cd_discrete_exact(n, R, scheme="central")` · `discrete_root(R_cell, scheme="central")` · `wiggle_indicator(T)` | D16: $T_j=\frac{r^j-1}{r^n-1}$, $r=\frac{1+R_{cell}/2}{1-R_{cell}/2}$ (central), $r=1+R_{cell}$ (upwind); R_cell = 2 degenerate (returns the limit $T_j=0$ for j < n) | dict(x, T, r) / float / dict(min, n_sign_changes, wiggles (bool)) — *expect* discrete_root(1) = 3, (4) = −3, upwind (4) = 5 | C09, E4 | A#17 |
| 1.30 | `advect_periodic(profile="square", C=0.8, n_cells=100, n_rev=1.0, scheme="upwind", growth_cap=1e6)` — profiles `"square"` (width 0.2), `"gauss"` (s₀ = 0.05); schemes `"upwind"`, `"maccormack"` (FB), `"maccormack_bf"`, `"lax_wendroff"`, `"ftcs"`, `"exact"` | (10.29), (10.101)–(10.102) on T_t + uT_x = 0 | dict(x, T, T0, exact, amplitude_ratio = max T / max T0, rms_error, blew_up, steps) | C05, C10, E3, A2, F3 | §8.3 |
| 1.31 | `ode_scheme_order(scheme, lam=-1.0, t_end=1.0, n_list=(20, 40, 80, 160))` — `"forward"`, `"backward"`, `"leapfrog"`, `"trapezoidal"` | (10.8) on y′ = λy | float (observed order) — *expect* 1, 1, 2, 2 (± 0.1) | C01 (N07), C13 (N97) | A#9 |
| 1.32 | `lax_demo(betas=(0.45, 0.50, 0.51), n_cells=20, nsteps=2000, T_wall=1.0, D=1.0)` | FTCS on the heated rod (N113) with Δt = βΔx²/D | dict per β: (max_error_history, blew_up, step_blown, final_error) — *expect* 0.45 and 0.50 bounded; 0.51 blown (∣G∣ = 1.04 at the zigzag, ≈ 900–1000 steps from round-off) | C05 (N18), A1 | **NEW** (was only a script, analysis row 28) |

### C.2 `fluidpy/core/fem1d.py` (NEW module, `FEM1`; re-exported by `ch10`) — 1-D Galerkin FE (Ch. 10, 11, 13, 16)
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 2.1 | `hat(x, x_nodes, A)` · `hat_basis(x_nodes, x)` (arrays) | (10.59)–(10.61); N₀ included as column 0 | ndarray / ndarray shape (len(x), n+1) — partition of unity, $N_A(x_B)=\delta_{AB}$ | C07, E5 | A#18 |
| 2.2 | `interpolate(d, g, x_nodes, x)` | (10.49), (10.62) | ndarray | C07, E5 | A#18 |
| 2.3 | `parent_shapes(xi)` · `map_to_element(xi, xa, xb)` · `map_to_parent(x, xa, xb)` · `shape_slopes(xa, xb, printed=False)` | (10.64)–(10.68); `printed=True` returns the book's labels (A, A+1) — slip R1 | tuple (N1, N2) / float [m] / float / dict(labels, slopes) — *expect* slopes (−1/h, +1/h) with labels (A−1, A) | C08, D13 | A#19 |
| 2.4 | `element_matrices_linear(h, u, D)` | D14: $\mathbf m^e=\frac h6\begin{pmatrix}2&1\\1&2\end{pmatrix}$, $\mathbf k^e=\frac u2\begin{pmatrix}-1&1\\-1&1\end{pmatrix}+\frac Dh\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$ | tuple (m, k) of 2 × 2 ndarrays — *expect* (0.25, 1, 0.25): m = [[0.08333, 0.04167], [0.04167, 0.08333]], k = [[0.5, −0.5], [−1.5, 1.5]] | C08, E5 | A#20 |
| 2.5 | `element_force(e, n_el, h, u, D, g=0.0, q=0.0)` | (10.77) | ndarray(2): e = 1 → −g k¹_{a1}; interior 0; e = n_el → Dq δ_{a2} | C08, E5 | **NEW** (split out of 2.7) |
| 2.6 | `element_integrals(xa, xb, u, D, quad=2)` | (10.75)–(10.76) by Gauss–Legendre | tuple (m, k) — equal to 2.4 to 1e-15 | C08 (from scratch) | A#20 |
| 2.7 | `connectivity(n_el, printed=False)` | (10.78): local a = 1 → A = e − 1, a = 2 → A = e (node 0 = Dirichlet); `printed=True` → A = e, e + 1 (R2) | ndarray (n_el, 2) int | C08, E5 | A#21, **CHG** (`printed`) |
| 2.8 | `assemble_1d(x_nodes, u, D, g=0.0, q=0.0, lumped=False, dirichlet_right=None)` (arrays) | (10.69)–(10.78): element loop, COO scatter-add → CSR; node 0 eliminated; `dirichlet_right=T_L` also eliminates node n (for (10.85)) | tuple (M, K, F) (scipy.sparse CSR, ndarray) | C07, C08, E5 | A#21, **CHG** (`dirichlet_right`) |
| 2.9 | `assembly_trace(n_el, u, D, L=1.0, g=0.0, q=0.0)` | the element loop step by step | list of dicts (e, nodes (global), m, k, f, M_after (dense), K_after (dense)) | E5 transport, C08 figure | §8.4 |
| 2.10 | `solve_steady(x_nodes, u, D, g=0.0, q=None, T_L=None)` | K d = F (∂/∂t = 0) | dict(x, T) — *expect* nodes np.linspace(0, 1, 5), u = 1, D = 0.25, T_L = 1 → [0, 0.025, 0.1, 0.325, 1] (= centred FD exactly) | C07, C08, C09, E5 | A#22 |
| 2.11 | `solve_transport(x_nodes, T0, u, D, g, q, dt, nsteps, theta=0.5, method="theta", lumped=False)` — methods `"theta"` (θ-scheme on M ḋ + K d = F), `"solve_ivp"` (method of lines, N29) | (10.58) | dict(t, T (final), history) | C07 (N29) | A#22 |
| 2.12 | `interior_stencil(h, u, D)` | (10.63) × h | dict(mass = (h/6, 2h/3, h/6), stiff = (−u/2 − D/h, 2D/h, u/2 − D/h)) — *expect* (0.25, 1, 0.25): stiff (−1.5, 2.0, −0.5) | C07, E5 | A#23 |
| 2.13 | `bilinear_form(w, v, x, u, D)` (callables or samples) · `weak_residual(T_fn, w_fn, u, D, q, L=1.0, g=0.0)` | (10.42); the steady (10.36) residual | float | C06, C07 | A#23 |
| 2.14 | `galerkin_equations_sympy(n=3)` | D11 steps (10.50)–(10.58) for n hats on a uniform mesh | dict(G_A (list), M, K, F (sympy Matrices), interior_row) | C07 (check cell) | A#24 |

### C.3 `fluidpy/core/mac.py` (NEW module, `MAC`; re-exported by `ch10`) — staggered-grid projection (Ch. 10–13)
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 3.1 | `MacGrid(nx, ny, Lx=1.0, Ly=1.0, periodic=(False, False))` (dataclass: dx, dy, xc, yc, xu, yv, shapes `p[ny, nx]`, `u[ny, nx+1]`, `v[ny+1, nx]`; `.zeros()`) | Fig. 10.4 | grid object | C12, E6, E7 | A#33 |
| 3.2 | `divergence(u, v, g)` · `gradient(p, g)` · `curl(u, v, g)` (corners) · `vorticity(u, v, g)` | (10.123), (10.126) | ndarray / tuple (gx on u-faces, gy on v-faces) / ndarray | C11, C12, C14, E6 | A#33 |
| 3.3 | `convective_terms(u, v, g, form="advective")` · `predictor(u, v, g, Re, dt, lid=1.0, body=None)` | (10.116), (10.119)–(10.120); face averages (ours, documented); tangential no-slip by ghost values | tuple (us, vs) | C11, C12, E7 | A#34 |
| 3.4 | `pressure_poisson_matrix(g, pin="corner")` · `solve_pressure(rhs, g, method="splu")` (methods `"splu"` (factor cached per grid), `"fft"` (periodic), `"sor"` (ch06 sweeps, `n_sweeps`)) · `pin_pressure(p, where="mean")` | (10.124) with boundary neighbours dropped (Neumann built in), one value pinned | csr / ndarray / ndarray | C12, E6 | A#35 |
| 3.5 | `correct(us, vs, p, g, dt)` · `project(us, vs, g, dt, return_parts=False)` | (10.117)–(10.118), (10.121)–(10.122) | tuple (u, v) / tuple (u, v, p) or dict (u, v, p, div_before, div_after, gx, gy) | C11, C12, E6 | A#36 |
| 3.6 | `projection_stages(us, vs, g, dt)` | the four stages of one projection | dict(div_before, rhs, rhs_sum, p, gx, gy, u, v, div_after, curl_correction) | C12, E6, A3 | §8.5 |
| 3.7 | `step(state, g, Re, dt, lid=1.0)` · `run(state, g, Re, dt, nsteps, callback=None)` · `dt_limit(umax, vmax, Re, dx, safety=0.8)` | predictor → Poisson → correction; (10.127) $\frac12(u^2+v^2)\Delta t\,Re\le1$ and (10.128) $\frac{4\Delta t}{Re\,\Delta x^2}\le1$ | dict state / float — *expect* dt_limit(1, 0, 100, 1/32, safety=1) = 0.02; (1, 0, 100, 1/64, 1) = 0.0061035 | C12, C14 | A#37 |
| 3.8 | `cavity(Re=100.0, n=32, t_end=20.0, tol_steady=1e-6, dt=None, save_every=None, cache=True)` | lid-driven cavity, lid u = 1 at y = 1 | dict(u, v, p, psi, omega, t, steps, converged, g, history (t, max∣Δu∣/Δt), snapshots) — cached to `outputs/ch10/mac_cavity_Re{Re}_n{n}.npz` | C14, C15, E7, A4 | A#38 |
| 3.9 | `cavity_centreline(state)` (u(y) at x = ½, interpolated between the two u-columns when nx is odd) · `cavity_centreline_v(state)` · `streamfunction(u, v, g)` (ψ = 0 on the walls, corners) · `primary_vortex_centre(psi, g)` | Fig. 10.8 analogue | dict(y, u) / dict(x, v) / ndarray / dict(x, y, psi_min) | C14, C15, E7 | §8.5 |
| 3.10 | `taylor_green(n, Re, t_end, periodic=True)` · `channel_poiseuille(nx, ny, Re, dpdx)` | V1 fields | dict(err_u, err_p, order) / dict(u, exact, max_err) | C12 (verification cell) | A#38 |
| 3.11 | `mac_projection_summary(n=8, kind="divergent", dt=1.0)` — deterministic fields (no RNG): `"divergent"` $u^*=\sin\pi x\cos\pi y+\tfrac12\sin2\pi x$, $v^*=\cos\pi x\sin\pi y$ sampled on the faces with wall-normal faces set to 0; `"shear"` $u^*=\sin^2\pi x\,\sin2\pi y$, $v^*=0.3\sin\pi x\,\sin\pi y$ | one projection | dict(div_before_max, div_after_max, rhs_sum, p_max, p_min (mean-zero p), curl_correction_max) | E6 (parity), C12 | **NEW** |

### C.4 `fluidpy/core/maccormack.py` (NEW module, `MCK`; re-exported by `ch10`) — explicit predictor–corrector (Ch. 10, 15)
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 4.1 | `maccormack_step(U, flux_E, flux_F, dt, dx, dy, arrangement="FF/BB")` (arrays; `flux_F=None` in 1-D) | (10.101)–(10.102) | ndarray | C10 | A#27 |
| 4.2 | `maccormack_advection_1d(T0, C, nsteps, arrangement="FB")` (periodic) | (10.101)–(10.102) with E = uT | ndarray — *expect* [0, 0, 1, 0, 0], C = 0.5, one step → [0, −0.125, 0.75, 0.375, 0] (= Lax–Wendroff) | C10, E3 | A#27 |
| 4.3 | `ns_fluxes(rho, rhou, rhov, c)` · `pressure_isothermal(rho, c)` | (10.96)–(10.99) | tuple of flux arrays / ndarray | C10 | A#26 |
| 4.4 | `ns_coefficients(dt, dx, dy, mu)` · `ns_predictor(state, coeffs, arrangement)` · `ns_corrector(state, star, coeffs, arrangement)` | (10.103)–(10.109) | dict c1…c5 / dict state | C10 | A#28 |
| 4.5 | `maccormack_dt(u, v, c, dx, dy, rho, mu, sigma=0.8)` · `maccormack_dt_asymptotic(Ma, dx, sigma=0.8)` (σ default ours) | (10.110), (10.155) | float — *expect* asymptotic(0.08, 1/64) = 7.071e-4; (0.08, 1/32) = 1.414e-3; maccormack_dt(1, 0, 12.5, 1/64, 1/64, 1, 0.01) = 2.935e-4 | C10, C15 (N91) | A#29 |
| 4.6 | `cavity_density_bc(state, side, stage, dt, dx, dy, U=1.0)` · `cavity_coefficients(dt, dx, dy, Ma, Re)` (→ a1…a11) · `weakly_compressible_step(state, coeffs, bc_fn, arrangement, dtype=np.float64, perturbation=True, printed_step5=False)` | (10.138)–(10.146), six substeps; `printed_step5=True` = slip R5 | dict / dict / dict state | C14 (N79–N85) | A#30, **CHG** (`printed_step5`) |
| 4.7 | `cavity_maccormack(Re=100.0, Ma=0.08, n=32, t_end=20.0, arrangement="FF/BB", cycle=False, tol_steady=1e-7, cache=True)` | explicit march to steady state | dict(u, v, rho, psi, t, steps, converged) — cached `outputs/ch10/mck_cavity_Re{Re}_n{n}.npz` | C10, C14, E7 (table) | A#31 |
| 4.8 | `block_wall_density(state, side, dx, dy, Ma, Re)` · `block_channel(Re=20.0, Ma=0.06, dx=0.125, t_end=30.0, H=4.0, ahead=8.0, behind=20.0, cycle=True, cache=True)` (our geometry — the book's is private) · `body_forces(state, mask, Re, Ma)` | (10.147)–(10.154) (heuristic); C_D (4.107), C_L (4.108) | dict / dict(t, CD, CL, state) / dict(CD, CL) | C15 (R11, N86–N90) | A#32 |

### C.5 `fluidpy/core/fem2d.py` (NEW module, `FEM2`; re-exported by `ch10`) — triangles, mixed elements (Ch. 10, 11, 14, 16)
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 5.1 | `p2_shape(xi, eta)` · `p2_shape_grad(xi, eta)` · `p1_shape(xi, eta)` | (10.185), (10.187) | ndarray(6) / ndarray(6, 2) / ndarray(3) — $\phi_a(\text{node}_b)=\delta_{ab}$, Σφ = 1 | C13, E8 | A#42 |
| 5.2 | `iso_map(xe, ye, xi, eta)` · `jacobian(xe, ye, xi, eta)` | (10.184), J = x_ξy_η − x_ηy_ξ | tuple (x, y) / float — *expect* affine element: J = 2 × area | C13 (N104) | A#42 |
| 5.3 | `tri_quad_7pt()` · `integrate_element(f, xe, ye)` | (10.198), Radon's degree-5 rule, ΣW = 1 | tuple (points (7, 2), W (7)) / float | C13 (N108), from scratch | A#42 |
| 5.4 | `structured_square_mesh(n, pair="P2P1")` | n × n squares, each cut along one diagonal | dict(vertices, tris, edges, vel_nodes, p_nodes, boundary_vel, counts = dict(V, E, T, N_u, N_p)) | C13, E8 | §8.7 |
| 5.5 | `p2_node_counts(V, E, T)` · `euler_check(mesh)` | N_u = V + E, N_p = V; V − E + T = 1 − holes | dict / int — *expect* n = 4 square: V 25, E 56, T 32, N_u 81 (= (2n + 1)²), N_p 25 | C13 (N94) | A#43 |
| 5.6 | `stokes_cavity(n, pair="P2P1", lid=1.0)` (dense ≤ 400 unknowns allowed; sparse above) | Stokes (Re → 0) lid-driven cavity, velocity Dirichlet everywhere, mean-zero pressure | dict(p (at pressure nodes), u, v, p_range, p_std_checker (the alternating component), n_u, n_p, residual) | C13, E8 | §8.7 |
| 5.7 | `infsup_constant(n, pair="P2P1")` (mesh = `structured_square_mesh(n, pair)`) | $\beta_h^2=\lambda_{\min}(\mathbf B^T\mathbf A^{-1}\mathbf B,\mathbf M_p)$ on mean-zero pressures (`scipy.linalg.eigh`) | dict(beta, n_spurious (eigenvalues < 1e-10 beyond the constant), n_u, n_p) — P2P1: β bounded (≈ constant in n), n_spurious 0; P1P1: n_spurious > 0 for n ≤ 4, β falling on refinement | C13, E8, F7 | A#45, **CHG** (takes n, returns counts) |
| 5.8 | `infsup_table(pairs=("P1P1", "P2P1"), n_list=(2, 3, 4, 6, 8, 12, 16))` | 5.7 in a loop; written to `reference/ch10/infsup_table.csv` (ours) | dict(pair → list of dicts) | E8 view 3, F7 | §8.7 |
| 5.9 | `stokes_p2p1(mesh, bc)` · `kovasznay_test(n, Re=40.0)` · `assemble_saddle(mesh, Re, …)` | Taylor–Hood Stokes/NS; Kovasznay exact field | dict / dict(err_u, err_p, orders) | C13 (verification cell) | A#45 |
| 5.10 | `cylinder_channel_mesh(d=1.0, W=5.0, xmin=-8.0, xmax=16.0, n_theta=48, n_r=12, h_far=0.4)` (our geometry — the book's is private) · `channel_bcs(mesh, U=1.0)` · `element_newton_matrices(…)` · `assemble_newton_system(mesh, state_star, state_n, Re, dt, alpha_t, beta_t, printed_10_172=False)` · `apply_dirichlet(A, b, nodes, values, method="replace")` · `newton_solve(mesh, state_n, Re, dt, alpha_t, beta_t, tol=1e-10, max_iter=12)` · `time_derivative(u_new, u_old, dudt_old, dt, alpha_t, beta_t)` · `cylinder_steady(Re, mesh_level="coarse", cache=True)` | (10.156)–(10.198) | dicts (mesh, residual history, fields, C_D) | C13 (N93–N110, R12) | A#43–44, **NEW** (`cylinder_steady` convenience + cache) |

### C.6 `fluidpy/ch10_computational_fluid_dynamics.py` (chapter module; re-exports C.1–C.5)
| # | Callable (signature) | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 6.1 | `weak_form_sympy()` · `weak_to_strong_sympy()` | D09, D10 with polynomial T, w on [0, L] | dict(steps, residual (0)) | C06 (check cells) | A#24 |
| 6.2 | `compressible_ns_sympy()` | (10.96)–(10.98) vs `navier_stokes_sym(mu_v=0)` | dict(difference (0)) | C10 (R09) | A#41 |
| 6.3 | `maccormack_linear_sympy()` | D18: two stages → Lax–Wendroff; local error; ∣G∣² | dict(update, lw_difference (0), local_error, G2, G2_factorised) | C10 (D18 check cell) | **NEW** |
| 6.4 | `projection_sympy()` | D19: ∇·(u* − Δt∇p) = 0 ⇒ ∇²p = ∇·u*/Δt; curl of ∇p = 0 | dict(poisson, curl (0)) | C11 | **NEW** |
| 6.5 | `weak_ns_identity_sympy()` · `weak_ns_components_sympy()` | D22; (10.156)–(10.159) | dict | C13 | A#41 |
| 6.6 | `split_linear_system(A1=None, A2=None, f1=None, f2=None)` (default: the non-commuting pair $A_1=\begin{pmatrix}1&1\\0&1\end{pmatrix}$, $A_2=\begin{pmatrix}1&0\\-1&2\end{pmatrix}$ — ours) · `marchuk_yanenko(A1, A2, f1, f2, phi0, dt, nsteps)` · `theta_scheme_linear(L_stokes, L_conv, y0, dt, nsteps, theta_split=1-1/np.sqrt(2), alpha_split=None)` · `theta_scheme_amplification_sympy()` · `splitting_order(method, dt_list=(0.1, 0.05, 0.025, 0.0125), t_end=1.0, theta_split=None)` | (10.111)–(10.114), (10.129)–(10.133) | dict(exact via `expm`, commutator) / ndarray / ndarray / dict(R_minus_exp series, z2_coeff, theta_root) / float (observed order) — *expect* MY order 1.0 ± 0.1; Θ-scheme 2.0 at θ = 1 − 1/√2, 1.0 at θ = 0.25 | C11 (N60, N61, N74), B1 | A#40, **NEW** (`splitting_order`) |
| 6.7 | `collocated_gradient(p, dx, dy)` · `checkerboard(nx, ny)` · `gradient_null_space(nx, ny, kind="collocated")` | (10.125), (10.126), SVD | tuple / ndarray / int (null-space dimension) — *expect* collocated ≥ 4 (periodic even grids), staggered 1 | C12, E6 | A#39 |
| 6.8 | `ghia_centreline(Re=100)` · `ghia_vortex_centre(Re=100)` · `hou_centres()` | read `reference/ch10/ghia1982_table1.csv` (Ghia, Ghia & Shin 1982, Table I; u on x = ½ and v on y = ½) and `hou1995_centres.csv` (Hou et al. 1995) | dict(y, u) / dict(x, y) / dict(Re → (x, y)) | C14, E7 | §8.6 |
| 6.9 | `cavity_error_vs_ghia(state, Re=100)` | interpolate our centreline at Ghia's y (`np.interp`) | dict(max_dev, rms_dev, rel_max (÷ max∣u∣)) | C14, C15, E7 | **NEW** |
| 6.10 | `strouhal_from_period(tau_bar)` · `dominant_frequency(t, signal, method="zero_crossings")` | S = 1/τ̄ (cyclic) | float | C13 (R13) | A#46 |
| 6.11 | `rod_heating_exact(x, t, D, T_wall, T_init=0.0, L=1.0, tol=1e-12)` | (10.199) with **our** wall temperature (default in examples T_wall = 1) and a truncation test | ndarray | C05 (N113), A1 | A#46, **CHG** (`L`) |
| 6.12 | `artificial_compressibility_channel(ny, Re, c=10.0, dtau=None, n_iter=20000, tol=1e-10)` | (10.95) pseudo-time march → plane Poiseuille | dict(u, exact, max_err, div_history) | C10 (N52) | A#25 |
| 6.13 | `book_slips()` | the slips table of convention 4 (id, printed, correct, evaluator) | list of dicts | callouts, tests | §8.9 |
| 6.14 | `cfl_time_step(c, dx, courant=1.0)` | Δt = C Δx / c (the climate CFL note) | float [s] — *expect* c = √(9.81 × 4000) = 198.1 m/s, Δx = 100 km → 504.8 s; Δx = 25 km → 126.2 s; c = 320 m/s, Δx = 25 km → 78.1 s | C05 | **NEW** |
| 6.15 | `derive_all()` | runs every sympy engine | dict of residual summaries | tests | A#41 |

### C.7 `tools/convergence.py` (machinery; orchestrator-approved addition)
| # | Callable | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 7.1 | `grid_convergence_index(f1, f2, f3, r=2.0, p=None, Fs=1.25)` (f1 finest) | observed $p=\ln\frac{f_3-f_2}{f_2-f_1}/\ln r$, Richardson $f_{ext}=f_1+\frac{f_1-f_2}{r^p-1}$, Roache GCI | dict(p, f_ext, gci_fine, gci_coarse, asymptotic_ratio, monotone (bool)) — *expect* f = 1 + h² at h = 0.025, 0.05, 0.1 → p = 2.000, f_ext = 1.000000 | C15, E7 | analysis §5 (NEW in tools) |

### C.8 `scripts/ch10_*.py` (runnable demos, each `--no-show`; heavy runs cache to `outputs/ch10/` with a parameter hash)
`ch10_drawings.py` (helpers, no physics: `spacetime_grid(ax)`, `stencil_diagram(ax, kind)`, `characteristic_diagram(ax, C)`,
`hat_functions(ax, n_el)`, `element_map(ax)`, `staggered_grid(ax, n=3, highlight=(2, 2))`, `p2p1_triangle(ax, iso=False)`,
`cavity_sketch(ax)`, `block_channel_sketch(ax)`, `cylinder_channel_sketch(ax)`) · `ch10_lax_demo.py` (A1) · `ch10_cavity.py`
(MAC 16–128 and MacCormack 32/64 at Re = 100, 400; writes the cache and the JS tables of E7) · `ch10_block.py` (Re = 20 at
Δx = 1/8, 1/16, 1/32 and Re = 100, our geometry; `--convergence`) · `ch10_cylinder_fem.py` (`--steady` Re = 1, 10, 40; `--unsteady` Re = 100
→ `reference/ch10/cylinder_fe_re100_forces.csv`, ours) · `ch10_tables.py` (writes `reference/ch10/infsup_table.csv` and the
JSON tables embedded in E7/E8, each value ≤ 4 significant figures, labelled "ours").

### C.9 Not in `analysis/ch10.md` §4 (flag list for the implementer)
`FD.stencil_taylor_coefficients`, `test_function`, `stencil_derivative`, `stencil_error`, `amplification_modulus`,
`worst_theta`, `stability_verdict`, `lax_demo`, the `"upwind_printed"` scheme (R11), `solve_transport_1d` growth cap and dict
return (CHG), `numerical_diffusivity(C=…)` (CHG) · `FEM1.element_force`, `connectivity(printed=)` (CHG), `assemble_1d(dirichlet_right=)`
(CHG) · `MAC.mac_projection_summary`, `vorticity` · `MCK.weakly_compressible_step(printed_step5=)` (CHG) · `FEM2.infsup_constant(n, pair)`
(CHG), `cylinder_steady` · `ch10.maccormack_linear_sympy`, `projection_sympy`, `splitting_order`, `cavity_error_vs_ghia`,
`cfl_time_step`, `richardson_three`, `rod_heating_exact(L=)` (CHG) · curation §8 items 1–10 are all included above (`truncation_terms`,
`amplification_curve`, `phase_error`, `advect_periodic`, `assembly_trace`, `projection_stages`, `cavity_centreline`,
`primary_vortex_centre`, `streamfunction`, `ghia_centreline`, `hou_centres`, `structured_square_mesh`, `stokes_cavity`,
`infsup_table`, `book_slips`, `grid_convergence_index`). Nothing in Parts A, B or F calls an unlisted name.

### C.10 Parity conventions (the explainers reproduce these exactly; the implementer fixes them in the docstrings)
| # | Callable / convention | Exact definition | Used by | Status |
|---|---|---|---|---|
| 10.1 | `advect_periodic` profiles and grid | nodes $x_i=i/N$, i = 0…N − 1 on [0, 1) periodic, u = 1; `"square"`: T = 1 where 0.2 ≤ xᵢ < 0.4, else 0; `"gauss"`: $T=e^{-(x-0.3)^2/(2\cdot0.05^2)}$ (nearest periodic image); one revolution = round(N/C) steps (C adjusted to N/steps so the pulse returns exactly); `amplitude_ratio` = max T(end)/max T0; `rms_error` against the exact shifted profile | E3, F3, A2, C05, C10 | **NEW** |
| 10.2 | `advected_gaussian` defaults | x₀ = 0.3, s₀ = 0.05, L = 1, n_images = 3, amplitude s₀/s(t) | E1 (FTCS mode), C02–C03 | **NEW** |
| 10.3 | `steady_cd_fd` grid | uniform nodes $x_j=jL/n$, j = 0…n, T₀ = 0, T_n = 1; stretched: $x_j=L\big[1-\frac{\tanh(\beta_s(1-j/n))}{\tanh\beta_s}\big]$ (clustered at x = L) with the nonuniform 3-point stencils for $T_x$ and $T_{xx}$ | E4, C09 | **NEW** |
| 10.4 | `MAC.cavity` time step and stopping | Δt = 0.8 × `dt_limit(1, 0, Re, 1/n, safety=1)` unless given; steady when max∣Δu∣/Δt < `tol_steady` (0 = run to t_end); lid u = 1 applied through the ghost value $u_{ghost}=2U-u_{top}$; predictor with central convective differences and the face averages of `MAC.predictor` | E7 (live 16²–32² runs), C14, C15 | **NEW** |
| 10.5 | `mac_projection_summary` fields | as in 3.11; grid `MacGrid(n, n)`, Δx = 1/n, wall-normal faces zero, p pinned by zero mean | E6 | **NEW** |
| 10.6 | `FEM2.structured_square_mesh(n, pair)` | vertices $(i/n, j/n)$; each square split by the diagonal from (i, j) to (i + 1, j + 1); P2 mid-edge nodes numbered after the vertices in edge order (horizontal, vertical, diagonal); all boundary velocities Dirichlet (lid u = 1 on y = 1 for `stokes_cavity`, corners u = 0) | E8 | **NEW** |
| 10.7 | `ch10.richardson_three(f1, f2, f3, r=2.0)` | the dict(p, f_ext) part of `grid_convergence_index`, available from the chapter module for parity rows (tools are not importable in shot.py's namespace) | E7, C15 | **NEW** |
| 10.8 | `ch10.ghia_centreline(Re)` column order | y ascending (0 … 1), 17 rows; `ch10.ghia_vortex_centre(Re)` → dict(x, y) | E7 | §8.6 |

---

## Part A — notebook storyboard (`notebooks/build_ch10.py` → `notebooks/ch10_computational_fluid_dynamics.ipynb`)

**One line per book section** (cell numbers are estimates for the builder's budget; ≈ 560 cells in all):
- §10.1 → N01 N02 N03 N04 (they open the C01 story; no A item of its own) — cells ≈ 8–16
- §10.2 → R01 R02 R03, C01 (N05 N07 · D01 · P221 P222 · F1 · **E1**), C02 (N06 N08 N09 · D02 · P223), C03 (N10 · D03 · P224),
  R04 R05, C04 (N11–N15, N17 · D04 D05 D06 D07 · P225 P226 · F2 · live · **E2**), C05 (N16 N18 N19 N113 · D08 · P227 P228 P229 ·
  **A1** · F3 · **E3**) — cells ≈ 17–215
- §10.3 → C06 (N20 N21 N22 · D09 D10 · P230 P231), C07 (N23–N33 · D11 D12 · P232 P233 · F4), C08 (N34–N39 · D13 D14 · P234
  P235 · **E5**) — cells ≈ 216–320
- §10.4 → R06 R07, C09 (N40 N42–N50 · D15 D16 D17 · P236 · F5 · live · **E4**), R09, C10 (N51–N59 · D18 · P237 P238 P239 ·
  **A2**), R08, C11 (N41 N60–N63 N72 N74 · D19 · P240 P241 · **A3**), C12 (N64–N71 N73 · D20 D21 · P242 P243 P244 P245 · **E6**),
  R10, C13 (N75 N76 N96 N100–N103 N107 in outline · D22 · P246 P247 P248 P249 P250 · F7 · **E8**) — cells ≈ 321–470
- §10.5 → C14 (N77–N85 · P251 P252 · **A4** · **E7**), R11, C15 (N86–N91 · D23 · F6), N92, **C13 (continued)**: R12 R13 N93–N95 N97–N99
  N104–N106 N108–N110 — cells ≈ 471–540
- §10.6 → **C15 (continued)**: N111 (P253), N112; S01 S02; summary — cells ≈ 541–560

Every CORE block below follows the order: problem in plain words → idea → primers → maths (notes and derivations, Part F) →
tiny example → code (fluidpy) + "What does the code above do?" → from-scratch check → visual(s) → notes and "What would
change if…". Code drafts are intent + exact calls; the builder writes every line commented (novice grade, units, the
equation written out next to its number). *expect* gives the numbers the executed cell must print. *see / read / change*
are the three figure notes. Note ids open every note text in bold (**N09 [B]**) so the lesson-reviewer can find them.

### A.0 Front matter
1. `nb.title(big_idea=…, roadmap=[…15…], prerequisites=[…])`. **Big idea (draft):** "Every flow in this book so far was
   solved with pencil and paper — and only because we could simplify it: a straight wall, a thin layer, a small wave. Real
   flows (a whole ocean basin, the air over a city, blood in a branching artery) do not simplify, so we hand the conservation
   laws of Ch. 4 to a computer. A computer cannot store a function; it stores numbers at points (finite differences) or the
   weights of simple building blocks (finite elements), and it replaces every derivative by arithmetic on those numbers.
   This chapter is about what that replacement costs and how to control it: how wrong a stencil is (order of accuracy), when
   an update rule blows up (von Neumann stability, the CFL condition), why a perfectly consistent scheme can still produce
   negative temperatures (the cell Péclet number), how the pressure is found in incompressible flow (the projection on a
   staggered grid — the same C-grid ocean and atmosphere models use), and finally how to earn trust in a CFD answer: a
   published benchmark and a grid-convergence study." **Roadmap (one line per CORE):** C01 stencils and their order ·
   C02 the explicit FTCS scheme · C03 consistency and the truncation error · C04 von Neumann stability · C05 upwinding and
   the CFL condition · C06 the weak form · C07 Galerkin with hat functions · C08 element matrices and assembly · C09
   convection-dominated flow: wiggles and numerical diffusion · C10 weak compressibility and MacCormack · C11 operator
   splitting and projection · C12 the staggered (MAC) grid · C13 mixed finite elements and the LBB condition · C14 the
   lid-driven cavity benchmark · C15 grid convergence and Richardson extrapolation. **Prerequisites:** Navier–Stokes and
   continuity (Ch. 4 §4.2, §4.5), dimensionless Navier–Stokes and Re (Ch. 4 §4.11), the diffusion equation and FTCS (Ch. 1
   §1.5, Ch. 8 §8.4 with Crank–Nicolson), Laplace/Poisson on a grid and relaxation (Ch. 6 §6.7), Fourier modes and plane
   waves (Ch. 5, Ch. 7), the lid-driven and cylinder flows' physics (Ch. 8, Ch. 9 §9.8).
2. `nb.explainer_index([...])` — 8 rows: ("fd_stencil_order", "What does 'second-order accurate' really mean?", "a stencil
   is a weighted sum of Taylor series; the first term it fails to cancel is its error, visible as a slope on log–log axes") ·
   ("von_neumann_amplification", "Why does β = 0.51 blow up when 0.50 does not?", "every Fourier mode of the round-off is
   multiplied by the same complex G(θ) each step; stability is ∣G∣ ≤ 1 for all θ") · ("upwind_cfl_advection", "Why must the time
   step shrink with the grid?", "the characteristic through the new point must start inside the stencil (CFL); first order
   pays in amplitude, second order in phase") · ("cell_peclet_wiggles", "Why does a centred scheme give negative
   temperatures?", "the discrete solution is rʲ with r < 0 once R_cell > 2; upwinding adds a diffusivity 0.5 R_cell D") ·
   ("fem_hat_assembly", "Where do finite-element matrices come from?", "each element adds a 2 × 2 block; the assembled rows
   are FD stencils with a ⅙–⅔–⅙ mass") · ("mac_projection_staggered", "What is the pressure doing in incompressible flow?",
   "one Poisson solve removes the divergence; on a staggered grid a zigzag pressure cannot hide") · ("lid_driven_cavity",
   "How do I know my CFD answer is right?", "reproduce a benchmark and watch the error fall at the expected order") ·
   ("mixed_fe_lbb", "Why can't velocity and pressure use the same elements?", "equal order leaves pressure modes the divergence
   cannot see; P2–P1 keeps the inf–sup constant bounded").
3. `nb.setup()`.
4. `nb.code` — **chapter imports** (outside any block): `import numpy as np` · `import sympy as sp` · `import matplotlib.pyplot
   as plt` · `import scipy.sparse as sps` · `from scipy import linalg, sparse` · `from fluidpy import
   ch10_computational_fluid_dynamics as ch10` · `from fluidpy.core import fd as FD, fem1d as FEM1, mac as MAC, maccormack as
   MCK, fem2d as FEM2` · `from fluidpy.core import diffusion as DIF, navier_stokes as NS, similarity as SIM, waves as WV,
   bluff_body as BB` · `from fluidpy import ch04_conservation_laws as ch04` · `from fluidpy.core.interact import slider_figure,
   animate_figure, live` · `from fluidpy.core.anim import animate` · `from fluidpy.core.style import COLORS, savefig` · `from
   tools.convergence import observed_order, pairwise_orders, richardson, grid_convergence_index` · `import sys;
   sys.path.insert(0, "scripts"); from ch10_drawings import spacetime_grid, stencil_diagram, characteristic_diagram,
   hat_functions, element_map, staggered_grid, p2p1_triangle, cavity_sketch, block_channel_sketch, cylinder_channel_sketch`.
   *explain:* one line per import ("`ch10` re-exports the five new toolkits, so `ch10.fd_weights` and `FD.fd_weights` are the
   same function"; "`scripts/ch10_drawings.py` only draws — no physics lives there").
5. `nb.md` — **⚠️ Conventions in this chapter**: the symbol table of convention 3 (book symbol · meanings · what we write),
   then "We compute with the book's symbols wherever the book's equations are quoted, and rename only where one letter would
   mean two things in the same line." A second table lists the twelve printed slips R1–R12 of convention 4 with where each
   is taught ("the book prints X, correct is Y"). A third short table: **three grids** (node-based FD §10.2 and MacCormack,
   element meshes §10.3 and §10.5, cell-centred staggered MAC §10.4) — which function uses which.
6. `nb.md` — **🔁 Tools from earlier chapters used in this one** (one line each, primer number): finite differences and FTCS
   (P21), `np.gradient` (P22), partial derivative (P25), Taylor expansion (P26, P98), power laws and log–log slopes (P13),
   `assert np.allclose` (P15), animate (P16), `slider_figure` (P17), `show_viz` (P18), sympy (P40) and its series (P117),
   Euler's formula and complex numbers (P45, P153), linear second-order ODE by an exponential trial (P44), chain rule (P49),
   product rule (P38), integration by parts (P218a), substitution in an integral (P106), `np.linalg.solve` (P57), null space
   and rank (P58), eigenvalues (P80), `expm` (P79), `np.expm1` (P107), Fourier modes and FFT Poisson (P142), Gauss–Legendre
   (P143), Poisson equation (P139), Newton's method (P152), Jacobi/Gauss–Seidel/SOR (P160), masks and `scipy.sparse` (P161),
   implicit stepping (P192), Crank–Nicolson and `solve_banded` (P193), `solve_ivp` (P31/P94), linear stability by eigenvalues
   (P214), characteristics of a first-order wave equation (P174), `np.interp` (P182), dicts (P23), f-strings (P04), lambda
   (P29), meshgrid/contour/streamplot (P76/P78), live widgets (P47). Each block repeats the ones it uses in a one-line
   reminder at first use.

---

### A.1 §10.1 Introduction — N01 N02 N03 N04 (open the C01 story)
1. `nb.section("10.1", "Introduction", intro="**What is this section about?** Computational fluid dynamics (CFD) turns the
   conservation laws of Ch. 4 into arithmetic a computer can do, and so predicts velocity, pressure and temperature fields,
   flow rates and forces for flows no formula can reach. This short section says what CFD is for, what can make its answers
   wrong, and which families of methods exist. The rest of the chapter builds two of them — finite differences and finite
   elements — from scratch on one small model problem, then applies them to the Navier–Stokes equations.")`
2. `nb.note` — **N01 [C]** "**Computational fluid dynamics** = quantitative flow predictions by computer from the
   conservation laws (mass $\\frac{\\partial\\rho}{\\partial t}+\\nabla\\cdot(\\rho\\mathbf u)=0$ (4.7), momentum $\\rho\\frac{D\\mathbf u}{Dt}=-\\nabla p+\\rho\\mathbf g+\\mu\\nabla^2\\mathbf u$ (4.39b), energy — Ch. 4). Pointer: the laws are Ch. 4; the methods are C01–C13 of
   this chapter."
3. `nb.note` — **N02 [B]** "**Four ways a CFD answer can be wrong** (a four-row table in our words): | source | example | how
   we control it |: *discretisation* (derivatives replaced by differences on a finite grid — the one we can **measure**, C03
   and C15) · *input data* (a viscosity or roughness known only to 10 %) · *initial and boundary conditions* (an outflow
   condition placed too close to a body — tested in N111) · *modelling* (a turbulence, cloud or sea-ice parameterisation in a
   climate model; Ch. 12). Only the first shrinks when you buy a bigger computer." + `nb.code`: a 3-line demo — the derivative
   of sin at 1 by the forward difference with h = 0.1 and 0.01 (`ch10.stencil_error("sin", 1.0, h, "forward")`). *expect:*
   −4.29e-2 and −4.2e-3 — "discretisation error shrinks with h; the other three would not".
4. `nb.note` — **N03 [C]** "**Why CFD is attractive:** cheap and fast compared with experiments, gives every quantity at
   every point, easy to change a parameter, can run conditions no laboratory can (a planet, an explosion). Computing power has
   grown exponentially for decades (Moore's law), which makes C15's question — *is the answer right?* — more, not less,
   important."
5. `nb.note` — **N04 [C]** "**Families of methods:** finite difference (C01–C05, C09–C12), finite element (C06–C08, C13),
   finite volume (a cell-balance cousin of both; the MAC grid of C12 is one), spectral (Fourier series; used for the KdV
   equation in Ch. 7 — `kdv_solve` — and for global weather models)."
6. `nb.md` — **ASCII pipeline** "from PDE to numbers":
   ```
   PDE  ∂T/∂t + u ∂T/∂x = D ∂²T/∂x²      (10.1)
     │  choose a grid  x_i = iΔx, t_n = nΔt
     ▼
   stencils  ∂T/∂x ≈ (T_{i+1} − T_{i−1})/(2Δx) + O(Δx²)     ── C01: how wrong is each?
     │  substitute
     ▼
   update rule  T_i^{n+1} = T_i^n − α(T_{i+1} − T_{i−1}) + β(T_{i+1} − 2T_i + T_{i−1})   ── C02
     │  ask three questions
     ▼
   consistent? (C03)   stable? (C04, C05)   ⇒ convergent (Lax, C05)
   ```

---

### A.2 §10.2 Finite-Difference Method — R01 R02 R03, C01 · C02 · C03 · R04 R05, C04 · C05
1. `nb.section("10.2", "Finite-Difference Method", intro="**What is this section about?** One model problem — a scalar T
   (a temperature, a dye, a pollutant) carried by a steady current u and spread by diffusion D — and everything a
   finite-difference method needs to solve it: stencils and their errors (C01), an explicit update rule (C02), the leftover
   when the exact solution is put into it (C03), the test that decides whether round-off errors grow (C04) and the time-step
   limit set by the flow speed (C05). The same five ideas return in every CFD code, including the dynamical core of every
   weather and climate model.")`
2. `nb.recap("R01", "The space–time grid", "Ch. 1 §1.5 stepped the diffusion equation on a grid: positions
   $x_i=i\\Delta x$, times $t_n=n\\Delta t$, and $T^n_i\\approx T(x_i,t_n)$ — superscript time level, subscript position. Our
   Fig. 10.1 analogue below draws three time levels.", where="Ch. 1 §1.5")` + `nb.figure` — `spacetime_grid(ax)` (three
   horizontal lines t_{n−1}, t_n, t_{n+1}; dots at x_0 = 0 … x_N = L; the FTCS stencil — three dots at t_n feeding one at
   t_{n+1} — highlighted orange). *see:* "a lattice of points in (x, t)"; *read:* "a scheme is a rule that computes the top dot
   from dots below it"; *change:* "…Δx halves: twice as many dots per row and (C04) four times as many rows for FTCS".
3. `nb.recap("R02", "Taylor series at the neighbours", "Ch. 1 P26 and Ch. 3 P98 expanded a smooth function about a point.
   At the grid neighbours $x_{i\\pm1}=x_i\\pm\\Delta x$: $T^n_{i\\pm1}=T^n_i\\pm\\Delta x\\big[\\tfrac{\\partial T}{\\partial x}\\big]^n_i
   +\\tfrac{\\Delta x^2}{2}\\big[\\tfrac{\\partial^2T}{\\partial x^2}\\big]^n_i\\pm\\tfrac{\\Delta x^3}{6}\\big[\\tfrac{\\partial^3T}
   {\\partial x^3}\\big]^n_i+\\tfrac{\\Delta x^4}{24}\\big[\\tfrac{\\partial^4T}{\\partial x^4}\\big]^n_i+O(\\Delta x^5)$ (10.4)–(10.5):
   the signs alternate only on the odd powers. These two lines are the starting line of D01.", where="Ch. 1 P26, Ch. 3 P98")`
4. `nb.recap("R03", "The centred second derivative", "Ch. 1's FTCS and Ch. 6's 5-point Laplacian used
   $\\big[\\tfrac{\\partial^2T}{\\partial x^2}\\big]^n_i=\\tfrac{T^n_{i+1}-2T^n_i+T^n_{i-1}}{\\Delta x^2}+O(\\Delta x^2)$ (10.7). D01
   step 8 re-derives it, with its error term.", where="Ch. 1 §1.5, Ch. 6 §6.7")`

#### C01 — Finite-difference stencils and their order of accuracy
5. `nb.core("C01", "Finite-difference stencils and their order of accuracy: $\\big[\\frac{\\partial T}{\\partial x}\\big]_i=
   \\frac{T_{i+1}-T_{i-1}}{2\\Delta x}+O(\\Delta x^2)$ (10.6)", question="A computer only knows T at grid points. How do we get a
   derivative from those numbers — and how wrong is the answer?")`
6. `nb.md` — **The problem in plain words:** "A string of thermometers hangs in a lake every 10 cm. You want the temperature
   gradient at one of them — the heat flux depends on it. You have three readings: the thermometer and its two neighbours.
   The gradient could be 'right minus me', 'me minus left' or 'right minus left over twice the spacing'. All three look
   reasonable; which is best, and how does the error shrink if you hang the thermometers twice as close? A climate model
   asks the same question every time step for every one of its millions of grid boxes (25 km apart for a modern global
   atmosphere model)."
7. `nb.md` — **The idea** (ASCII):
   ```
   T(x_i+Δx) = T + Δx T′ + Δx²/2 T″ + Δx³/6 T‴ + …     ← Taylor at the right neighbour
   T(x_i−Δx) = T − Δx T′ + Δx²/2 T″ − Δx³/6 T‴ + …     ← Taylor at the left neighbour
   ───────────────────────────────────────────────
   subtract:  T_{i+1} − T_{i−1} = 2Δx T′ + 0 + Δx³/3 T‴ + …   the T″ terms CANCEL
   ÷ 2Δx:     (T_{i+1} − T_{i−1})/(2Δx) = T′ + Δx²/6 T‴ + …    first survivor ∝ Δx²  ⇒ second order
   ```
   "A stencil is a weighted sum of Taylor series. Choose the weights so the low terms cancel; **the first term that survives is
   the error**, and its power of Δx is the order."
8. `nb.primer("big-O notation and the order of accuracy", "O(Δx²) means 'a quantity no bigger than a constant times Δx² when
   Δx is small' — it shrinks by 4 when Δx halves. A stencil is *p-th order accurate* when its error is O(Δxᵖ). Unlike the
   'orders of smallness' of Ch. 2 (P68), where small terms were *dropped*, here the leftover is *kept* and measured: it is the
   price of the approximation.", code="import numpy as np\nfor h in (0.1, 0.05, 0.025):                      # halve the spacing twice\n    err = (np.sin(1 + h) - np.sin(1 - h))/(2*h) - np.cos(1)   # centred difference minus the exact slope\n    print(h, err, err/h**2)                        # err/h^2 stays ≈ -0.090: the error is O(h^2)")`
   (**P221**)
9. `nb.note` — **N05 [B]** "**The model problem** used from here to C05: $\\frac{\\partial T}{\\partial t}+u\\frac{\\partial T}
   {\\partial x}=D\\frac{\\partial^2T}{\\partial x^2}$, $0\\le x\\le L$ (10.1) — T carried by a constant current u [m/s] and
   spread with diffusivity D [m²/s]. Its exact solution for a Gaussian start is a Gaussian that travels at u and widens like
   $s(t)^2=s_0^2+2Dt$ (Ch. 8's diffusion spreading), our test field: `FD.advected_gaussian`." equation (10.1) + `nb.code`
   `x = np.linspace(0, 1, 201); T0 = FD.advected_gaussian(x, 0.0, 0.5, 0.01); T1 = FD.advected_gaussian(x, 0.4, 0.5, 0.01)`;
   print the peak position and height. *expect:* peak moves 0.3 → 0.5 m, height 1 → s₀/s = 0.05/√(0.0025 + 0.008) = 0.488.
10. `nb.derivation("D01", …)` — Part F D01 (9 steps), ref "10.6". (Uses R02 Taylor, P221 big-O.)
11. `nb.worked_example("the slope of sin at x = 1 with h = 0.1", "1. Exact slope cos 1 = 0.540302. 2. Forward: (sin 1.1 − sin
    1)/0.1 = (0.891207 − 0.841471)/0.1 = 0.497364; error −0.04294, and the leading term of D01 predicts
    $\\frac h2T''=\\frac{0.1}{2}(-0.841471)=-0.04207$ (the small rest is the next term, $\\frac{h^2}6T'''$). 3. Centred: (sin 1.1 − sin 0.9)/0.2 = (0.891207 − 0.783327)/0.2 =
    0.539402; error −0.000900, predicted $\\frac{h^2}6T'''=\\frac{0.01}6(-\\cos1)=-0.000900$ ✓. 4. Halve h to 0.05: forward
    error −0.0213 (÷ 2.02), centred −0.000225 (÷ 4.00) — first and second order, as the Taylor bookkeeping said.")`
12. `nb.code` — `for kind, offs, m in (("forward", [0, 1], 1), ("backward", [-1, 0], 1), ("central", [-1, 0, 1], 1),
    ("central2", [-1, 0, 1], 2), ("onesided2", [0, 1, 2], 1)): print(kind, FD.fd_weights(offs, m),
    FD.stencil_taylor_coefficients(kind, m)["order"])`; then a table over h ∈ {0.1, 0.05, 0.025, 0.0125} of
    `FD.stencil_error("sin", 1.0, h, kind, m)` and `observed_order(h_list, abs(errors))`. *expect:* weights (−1, 1), (−1, 1),
    (−1/2, 0, 1/2), (1, −2, 1), (−3/2, 2, −1/2); orders 1, 1, 2, 2, 2; observed 1.00, 1.00, 2.00, 2.00, 1.99 (± 0.02).
    *explain:* 1. `fd_weights` solves the Taylor-matching equations exactly (fractions); 2. `stencil_taylor_coefficients`
    lists which Taylor terms survive; 3. `stencil_error` applies the stencil to sin and subtracts the exact derivative;
    4. `observed_order` is the log–log slope (P13).
13. `nb.check_agree` — **from scratch (curation §7):** build the 3 × 3 Taylor-matching matrix for offsets [−1, 0, 1] (rows
    k = 0, 1, 2: $s^k/k!$) and for [0, 1, 2], solve with `np.linalg.solve` (P57) for the first-derivative weights, compare:
    `assert np.allclose(w_mine, [float(c) for c in FD.fd_weights([-1, 0, 1], 1)])` and the same for [0, 1, 2]. Then the three
    stencils by array slicing on `f = np.sin(x)` vs `FD.fd_derivative(f, dx, kind)`: `assert np.allclose(mine[1:-1],
    lib[1:-1])`. *explain:* "the library does exactly this Vandermonde solve".
14. `nb.note` — **N07 [B]** "**Time differences (10.8):** the same three moves in t give $\\big[\\frac{\\partial T}{\\partial
    t}\\big]^n_i=\\frac{T^{n+1}_i-T^n_i}{\\Delta t}+O(\\Delta t)$ (forward), $\\frac{T^n_i-T^{n-1}_i}{\\Delta t}+O(\\Delta t)$
    (backward) and $\\frac{T^{n+1}_i-T^{n-1}_i}{2\\Delta t}+O(\\Delta t^2)$ (centred, 'leapfrog' — the time step of many older
    weather models). Measured on y′ = −y below." equation (10.8) + `nb.code`: `[FD.ode_scheme_order(s) for s in ("forward",
    "backward", "leapfrog")]`. *expect:* ≈ 1.0, 1.0, 2.0.
15. `nb.primer("floating-point round-off and machine epsilon", "A computer stores about 16 significant digits (double
    precision): `np.finfo(float).eps` = 2.2e-16 is the gap after 1.0. A difference quotient subtracts two nearly equal numbers
    and divides by h, so its round-off error is about ε·∣f∣/h — it *grows* as h shrinks. Truncation falls like hᵖ, round-off
    rises like 1/h: the total error has a floor, near h ≈ ε^{1/2} for a first-order and ε^{1/3} for a second-order stencil.
    Single precision (float32, ε ≈ 1.2e-7) hits the floor much sooner.", code="import numpy as np\nprint(np.finfo(float).eps, np.finfo(np.float32).eps)   # 2.2e-16 and 1.2e-07\nfor h in (1e-4, 1e-8, 1e-12):                       # shrink h far too much\n    print(h, (np.sin(1 + h) - np.sin(1))/h - np.cos(1))   # error falls, then RISES: about 4e-5, 3e-9, 4e-5")`
    (**P222**)
16. `nb.figure` — **"The order is a slope"** (7 × 4 in, log–log): ∣error∣ of the five stencils on sin at x₀ = 1 for
    h = logspace(−12, −0.5, 120) (forward teal, backward teal dashed, central orange, central2 orange dashed, onesided2 purple);
    slope-1 and slope-2 guide lines (muted, labelled "∝ h", "∝ h²"); a shaded round-off band ε/h (rose, labelled "round-off
    ≈ ε/h"); dots at the best h of forward (≈ 1.6e-9, error 4.5e-11) and central (≈ 2.2e-6, error 6.6e-13). Title "Halve h:
    first order halves the error, second order quarters it — until round-off wins". *see:* "straight lines of slope 1 and 2 on
    the right, a V-shaped floor on the left"; *read:* "the slope is the order; the lowest point is the best h — about 10⁻⁸ for
    forward, 10⁻⁵ for central"; *change:* "…we computed in float32: the floor rises to ~10⁻³ and moves to h ≈ 10⁻³ (P222) —
    the reason §10.5 insists on double precision (N90)".
17. `nb.plotly` — **F1** `slider_figure` over h (25 log-spaced values 0.2 … 1e-3): left trace set = the five stencil errors
    vs h (static) with a moving marker at the current h; a second panel as bars = the Taylor terms of the central stencil at
    that h (T′ kept blue, T‴ term rose): uses `FD.stencil_taylor_coefficients` and `FD.test_function`. *explain:* "the page
    version of E1's term bars".
18. `nb.explainer("fd_stencil_order", heading="What does 'second-order accurate' really mean?", why="Dragging h on a log
    slider moves one dot along the error curve while the stencil nodes and their secant close in on x₀, and the Taylor-term
    bars show which terms cancel — the order becomes something you see happen, not a label.", tries=["Pick 'forward' and halve
    h twice: watch the error halve each time; switch to 'central' and watch it quarter.", "Drag h below 10⁻⁶: the status turns
    to 'round-off dominated' — the dot climbs back up.", "Open the FTCS mode: the three truncation terms of (10.17) appear as
    bars; set Δt so that DΔt/Δx² = 1/6 and watch the time and diffusion terms cancel."])`
19. `nb.md` — **What would change if…** "…we used the stencils on a PDE instead of a known function? The error no longer
    shows up in one derivative but in the whole update rule — C02 builds that rule for (10.1), and C03 finds its leftover."

#### C02 — The explicit FTCS scheme (with BTCS as its implicit twin)
20. `nb.core("C02", "The explicit FTCS scheme $T^{n+1}_i=T^n_i-\\alpha(T^n_{i+1}-T^n_{i-1})+\\beta(T^n_{i+1}-2T^n_i+T^n_{i-1})$
    (10.10)", question="How do we march a temperature carried and spread by a current forward in time, one grid point at a
    time?")`
21. `nb.md` — **The problem in plain words:** "A factory releases a puff of dye into a river flowing at 0.1 m/s. Downstream the
    puff drifts with the current and slowly spreads. You know the concentration along the river now, at points 10 cm apart,
    and want it 0.2 s later — then 0.2 s after that, and so on. The simplest rule computes each new value from the three old
    values around it: forward in time, centred in space (FTCS)."
22. `nb.md` — **The idea** (ASCII stencil):
    ```
    t_{n+1}              ● T_i^{n+1}            = T_i − α(T_{i+1} − T_{i−1}) + β(T_{i+1} − 2T_i + T_{i−1})
                       ╱ │ ╲
    t_n          ●─────●─────●                  α = uΔt/(2Δx)  (how far the current moves in Δt, per 2 cells)
              T_{i−1}  T_i  T_{i+1}              β = DΔt/Δx²    (how far diffusion reaches in Δt, per cell²)
    ```
    "Every new value is a weighted average of three old ones. With u = 0 this is Ch. 1's FTCS with r = β."
23. `nb.derivation("D02", …)` — Part F D02 (6 steps), ref "10.10". (Uses N07 forward time difference, D01 stencils.)
24. `nb.worked_example("one FTCS step on five points", "u = 0.1 m/s, D = 0.01 m²/s, Δx = 0.1 m, Δt = 0.2 s. 1. α = uΔt/(2Δx) =
    0.1 × 0.2/0.2 = 0.1; β = DΔt/Δx² = 0.01 × 0.2/0.01 = 0.2. 2. Start T = [0, 0, 1, 0, 0] (a spike in the middle). 3. Middle:
    T₂ = 1 − 0.1(0 − 0) + 0.2(0 − 2 + 0) = 0.6. 4. Left neighbour: T₁ = 0 − 0.1(1 − 0) + 0.2(1 − 0 + 0) = 0.1. 5. Right
    neighbour: T₃ = 0 − 0.1(0 − 1) + 0.2(0 − 0 + 1) = 0.3. 6. New T = [0, 0.1, 0.6, 0.3, 0]: the sum is still 1 (nothing is
    lost), the peak dropped (diffusion) and the right side got more than the left (the current carries it right).")`
25. `nb.code` — `alpha, beta = FD.ftcs_coefficients(0.1, 0.01, 0.1, 0.2)`; `T = np.array([0, 0, 1, 0, 0.0])`;
    `print(alpha, beta, FD.transport_1d_step(T, alpha, beta, "ftcs", periodic=True))`; then the Gaussian run: `x =
    np.linspace(0, 1, 101)[:-1]` (periodic, Δx = 0.01), u = 0.5, D = 0.01, β = 0.25 ⇒ Δt = 0.0025 s, `run =
    FD.solve_transport_1d(FD.advected_gaussian(x, 0, 0.5, 0.01), x, 0.5, 0.01, 0.0025, 80, "ftcs", periodic=True)`; print
    `FD.error_norm(run["T"], FD.advected_gaussian(x, 0.2, 0.5, 0.01))`. *expect:* α = 0.1, β = 0.2; [0, 0.1, 0.6, 0.3, 0];
    rms error ≈ 1e-3 (the cell prints it; the assert is < 5e-3). *explain:* 1. `ftcs_coefficients` = (10.11); 2. one step of
    (10.10); 3. 80 steps to t = 0.2 s and the rms error (10.14).
26. `nb.check_agree` — **from scratch:** the FTCS update in one line with `np.roll` (periodic neighbours, P77):
    `T_new = T - alpha*(np.roll(T, -1) - np.roll(T, 1)) + beta*(np.roll(T, -1) - 2*T + np.roll(T, 1))`; `assert
    np.allclose(T_new, FD.transport_1d_step(T, alpha, beta, "ftcs", periodic=True))`. Then u = 0 parity with Ch. 1: `assert
    np.allclose(FD.solve_transport_1d(T0, x, 0.0, D, dt, 40, "ftcs", g=0.0, q=None)["T"], DIF.ftcs_diffusion_1d(T0, D, dx, dt,
    40)[...])` (same Dirichlet ends; the builder matches the ch01 call's boundary arguments).
27. `nb.primer("ghost node for a Neumann boundary", "At x = L the book prescribes a slope, ∂T/∂x = q (10.2), not a value.
    Invent one extra point beyond the end, $T_{N+1}$, chosen so the centred difference at the end gives the slope:
    $(T_{N+1}-T_{N-1})/(2\\Delta x)=q$, i.e. $T_{N+1}=T_{N-1}+2\\Delta x\\,q$. The ordinary FTCS formula then runs at the last node
    too and stays second order; q = 0 is an insulated end.", code="import numpy as np\nT, dx, q = np.array([0.0, 0.4, 0.7, 0.9]), 0.1, 0.5   # last value 0.9 at x = L; required slope 0.5\nghost = T[-2] + 2*dx*q                            # T_{N+1} = T_{N-1} + 2Δx q = 0.7 + 0.1 = 0.8\nprint(ghost, (ghost - T[-2])/(2*dx))                # 0.8 and the centred slope 0.5 ✓")` (**P223**)
28. `nb.note` — **N06 [B]** "**Boundary and initial conditions (10.2)–(10.3):** $T(0,t)=g$ (Dirichlet: the value is
    prescribed — a heater held at fixed temperature) and $\\frac{\\partial T}{\\partial x}(L,t)=q$ (Neumann: the gradient, i.e.
    the diffusive flux, is prescribed — q = 0 is an insulated end); $T(x,0)=T_0(x)$ must agree with them. In code: `g=` sets
    node 0, `q=` uses the ghost node (P223)." equation $T(0,t)=g,\\ \\partial_xT(L,t)=q,\\ T(x,0)=T_0(x)$ (10.2)–(10.3) +
    `nb.code`: an insulated rod (q = 0, g = 1, u = 0) marched 400 steps; print T at x = L. *expect:* rises toward 1 (the rod
    fills up from the heated end); the ghost node keeps the end slope at 0 to round-off.
29. `nb.note` — **N08 [B]** "**The implicit twin — BTCS (10.12)–(10.13):** evaluate the space differences at the *new* level:
    $T^n_i+\\alpha(T^n_{i+1}-T^n_{i-1})-\\beta(T^n_{i+1}-2T^n_i+T^n_{i-1})\\approx T^{n-1}_i$ (10.13). Now each equation contains
    three unknowns, so a whole tridiagonal system is solved per step (Ch. 8 P193 `solve_banded` reminder). One step on our five
    points (α = 0.1, β = 0.2) gives the matrix row (−0.3, 1.4, −0.1) and the solve is printed." equation (10.13) + `nb.code`:
    `FD.transport_1d_step(T, 0.1, 0.2, "btcs", periodic=True)` and the 5 × 5 matrix built by hand, `np.linalg.solve`, `assert
    np.allclose`. *expect:* row coefficients (−(α + β), 1 + 2β, α − β) = (−0.3, 1.4, −0.1); the result sums to 1.
30. `nb.note` — **N09 [B]** "**Explicit or implicit?**" a two-column table | | explicit (FTCS) | implicit (BTCS) |: one step
    = a formula per node | a linear solve; cost per step tiny | larger (tridiagonal: still O(N)); time step limited by
    stability (C04: β ≤ ½) | any Δt is stable (D07) but accuracy still limits it; parallelises trivially | needs a solver.
    "Climate models use both: explicit for fast advection, implicit (or semi-implicit) for the fastest gravity waves and
    vertical diffusion."
31. `nb.figure` — **"FTCS and BTCS follow the drifting, spreading puff"** (8 × 3.4 in, two panels): (a) profiles at t = 0,
    0.1, 0.2 s — exact (muted dashed), FTCS dots (orange), BTCS dots (blue) — for u = 0.5, D = 0.01, Δx = 0.01, β = 0.25;
    (b) the same at β = 2.5 (Δt ten times larger): BTCS still smooth (but smeared, first order in time), FTCS not run
    ("unstable — see C04", a rose label). *see:* "the peak moves right and flattens; both schemes sit on the exact curve at
    β = 0.25"; *read:* "at β = 2.5 BTCS survives but lags and smears: stability is not accuracy"; *change:* "…u = 0: the peak
    stays put and only spreads — Ch. 1's diffusion".
32. `nb.md` — **What would change if…** "…we asked *how close* FTCS is to (10.1) itself, not to one solution? Put the exact
    solution into the rule and see what is left over — C03."

#### C03 — Consistency and the truncation error
33. `nb.core("C03", "Consistency and the truncation error $E^n_i=\\frac{\\Delta t}{2}T_{tt}+u\\frac{\\Delta x^2}{6}T_{xxx}-
    D\\frac{\\Delta x^2}{12}T_{xxxx}+O(\\Delta t^2,\\Delta x^4)$ (10.17)", question="Does the update rule really approximate the
    PDE — and which equation does it solve exactly?")`
34. `nb.md` — **The problem in plain words:** "You wrote a rule that *looks* like (10.1). A typo, a wrong sign or a missing
    factor 2 would also look plausible. Consistency is the test: put the true solution into the rule. If what is left over
    shrinks to zero as Δx and Δt shrink, the rule approximates (10.1); the size of the leftover tells how fast. The leftover
    also tells the *character* of the error — whether it smears (like extra diffusion) or ripples (like dispersion)."
35. `nb.md` — **The idea:** "scheme(exact T) = PDE(exact T) + E = 0 + E" as a two-line box, and the three words: **convergent**
    (the numbers approach the exact solution), **consistent** (the equations approach the PDE), **stable** (errors do not
    grow) — a three-row table with how each is tested (measure the error; Taylor-expand the scheme; von Neumann, C04).
36. `nb.primer("norms of an error array: rms, max, L1", "An error is an array — one number per grid point. To judge it by one
    number: the root mean square √(mean e²) (the book's choice in (10.15)), the maximum ∣e∣ (the worst point) or the mean ∣e∣.
    For a smooth error they shrink at the same rate; near a shock or a wiggle the max norm is the honest one.",
    code="import numpy as np\ne = np.array([0.0, 0.01, -0.02, 0.01])      # an error at four points\nprint(np.sqrt(np.mean(e**2)), np.max(np.abs(e)), np.mean(np.abs(e)))   # 0.0122, 0.02, 0.01")` (**P224**)
37. `nb.note` — **N10 [B]** "**Convergence and its rates (10.14)–(10.15):** the solution error $e^n_i=T^n_i-T(x_i,t_n)$ (10.14)
    and, for a convergent scheme, $\\lVert e^n_i\\rVert\\le K_e\\,\\Delta x^a\\,\\Delta t^b$ (10.15) with rates a, b (we write $K_e$
    so it is not confused with the stiffness matrix $\\mathbf K$ of C07). To measure a alone, tie Δt to Δx (Δt ∝ Δx² makes the
    Δt error as small as the Δx one); D23 turns (10.15) into the slope of a log–log plot." equations (10.14), (10.15) +
    `nb.code`: `st = FD.convergence_study("ftcs", (20, 40, 80, 160), dt_rule="diffusive")`; print `st["dx"], st["err"],
    st["order"], st["pairwise"]`. *expect:* order 2.0 ± 0.15 (pairwise ≈ 1.9–2.05).
38. `nb.derivation("D03", …)` — Part F D03 (11 steps), ref "10.17". (Uses P98 two-variable Taylor, D01 results, P117 sympy
    series; `check_src` optional — ★★, but the executed sympy cell of item 40 repeats it.)
39. `nb.worked_example("the leftover for pure diffusion", "Take u = 0, D = 1, and the exact solution T = e^{−π²t} sin πx (a
    sine that decays). Its derivatives: $T_{tt}=\\pi^4T$ (differentiate $e^{-\\pi^2t}$ twice) and $T_{xxxx}=\\pi^4T$. 1. Insert in
    (10.17): $E=\\frac{\\Delta t}2\\pi^4T-\\frac{\\Delta x^2}{12}\\pi^4T=\\pi^4T\\Delta x^2\\big(\\frac\\beta2-\\frac1{12}\\big)$ using
    Δt = βΔx² (D = 1). 2. At x = ½, t = 0 (T = 1), Δx = 0.1, β = 0.25: E = 97.41 × 0.01 × (0.125 − 0.0833) = 0.0406. 3. Surprise:
    at β = 1/6 the bracket is zero — the two leading errors cancel and FTCS becomes fourth order in space for pure
    diffusion. (Try it in E1's FTCS mode.)")`
40. `nb.code` — `res = FD.truncation_error_sympy("ftcs")`; `print(res["E"], res["order_t"], res["order_x"])`; then
    `tt = FD.truncation_terms(0.5, 0.01, 0.01, 0.0025, 0.35, 0.1)`; print the four entries. *expect:* E printed with the three
    terms of (10.17); orders (1, 2); `tt["total"]` agrees with `tt["measured"]` to within 5 % (the next Taylor terms).
    *explain:* 1. sympy substitutes two-variable Taylor series (P98) into (10.10), divides by Δt and subtracts (10.1); 2. the
    numeric terms on the travelling Gaussian; 3. the measured one-step residual.
41. `nb.check_agree` — **from scratch:** the one-step residual of FTCS on the exact Gaussian,
    `(T_ex(t+dt) - FD.transport_1d_step(T_ex(t), alpha, beta, "ftcs", periodic=True))/dt`, at the node nearest x = 0.35,
    against `tt["total"]`: `assert abs(resid - tt["total"]) < 0.05*abs(tt["total"])`; and the sign convention stated: E is on
    the *left* of (10.16), so "exact minus scheme" = +E Δt.
42. `nb.figure` — **"Consistent: the error falls at the rate the leftover predicts"** (8 × 3.4 in, two panels, log–log):
    (a) rms error vs Δx with Δt ∝ Δx² (β = 0.25 fixed: slope 2, orange) and with Δt ∝ Δx (C = 0.5 fixed: slope 1 once the
    Δt/2 T_tt term dominates, orange dashed) — `FD.convergence_study` both rules; (b) the three terms of (10.17) at the pulse
    centre vs Δx (time term amber, dispersive T_xxx term purple, diffusive T_xxxx term teal) — a bar trio at three Δx. *see:*
    "two straight lines of different slope"; *read:* "tie Δt to Δx² and the scheme looks second order; tie it to Δx and the
    first-order time error shows"; *change:* "…β = 1/6: in pure diffusion the bars of the time and T_xxxx terms cancel (see the
    tiny example)".
43. `nb.md` — **What would change if…** "…β = 0.51? The leftover is still tiny and still shrinks with Δx — the scheme is
    consistent — yet the computation explodes. Consistency is necessary, not sufficient: C04 finds the missing ingredient."
44. `nb.recap("R04", "Fourier modes", "Ch. 5 P142 and Ch. 7 wrote any periodic grid function as a sum of waves $e^{\\mathrm
    ikx}$. The book writes the error as $\\xi^n_i=\\sum_kg^n(k)e^{\\mathrm i\\pi kx_i}$ (10.20) — its k counts *half*-waves per unit
    length, so the phase change per cell is $\\theta=k\\pi\\Delta x$ (10.26). ⚠️ The book's i is both the grid index and
    $\\sqrt{-1}$; we write the imaginary unit upright ($\\mathrm i$) and, in our own lines, call the grid index j.
    `FD.fourier_mode(x, k, convention='book' or 'standard')` shows both conventions.", where="Ch. 5 (P142), Ch. 7 §7.1")`
    + `nb.code`: `x = np.arange(8)*0.125; print(np.round(FD.fourier_mode(x, 8, "book").real, 3))` *expect:* the zigzag
    [1, −1, 1, −1, …] — θ = π, the shortest wave a grid can carry.
45. `nb.recap("R05", "The diffusion limit β ≤ ½", "Ch. 1's FTCS for pure diffusion was stable only for
    $r=D\\Delta t/\\Delta y^2\\le\\frac12$ (`core.diffusion.stable_time_step`). Here it reappears as $0\\le\\beta\\le\\frac12$, i.e.
    $\\Delta t\\le\\frac12\\frac{\\Delta x^2}D$ (10.28); D06 derives it, and adds the new half: with D = 0 FTCS is never stable.",
    where="Ch. 1 §1.5")`

#### C04 — Von Neumann stability: the amplification factor and Noye's region
46. `nb.core("C04", "Von Neumann stability: $G=(\\alpha+\\beta)e^{-\\mathrm i\\theta}+(1-2\\beta)+(\\beta-\\alpha)e^{\\mathrm
    i\\theta}$ (10.24), stable iff $0\\le4\\alpha^2\\le2\\beta\\le1$ (10.27)", question="Why does a round-off error of 10⁻¹⁶ explode
    into garbage with β = 0.51 but die away with β = 0.50?")`
47. `nb.md` — **The problem in plain words:** "Run FTCS for pure diffusion on a rod, twice: with β = 0.50 and with β = 0.51.
    The two runs agree for hundreds of steps, then the second one grows a zigzag — neighbouring nodes alternating up and down —
    that doubles every ~18 steps until the numbers overflow. Nothing in the physics changed by 2 %. What happened is that the
    computer's own rounding errors (always present, about 10⁻¹⁶) contain every wavelength, and for β = 0.51 the shortest one
    is amplified a little every step." (A 2-line code cell runs `FD.lax_demo(betas=(0.50, 0.51))` and prints `step_blown`.)
    *expect:* None for 0.50, ≈ 900–1000 for 0.51.
48. `nb.md` — **The idea:** "Errors obey the same linear rule as the solution (N12). Break the error into waves; the rule
    treats each wave separately and multiplies its amplitude by one complex number G(θ) per step. After n steps the amplitude
    is ∣G∣ⁿ times the start: **stable ⇔ ∣G(θ)∣ ≤ 1 for every wave**." ASCII:
    ```
    error now = Σ  amplitude_θ · wave_θ          one step later: Σ G(θ)·amplitude_θ · wave_θ
    |G| = 0.98 : 10⁻¹⁶ → 10⁻¹⁶·0.98ⁿ  dies           |G| = 1.04 : 10⁻¹⁶ → 10⁻¹⁶·1.04ⁿ ≈ 1 after ~940 steps
    ```
49. `nb.note` — **N11 [B]** "**Stable** (§10.2): a scheme is stable when small disturbances — round-off made at any step —
    decay, or at least stay bounded, instead of taking over the solution. It is a property of the *rule*, not of the PDE:
    (10.1) itself is perfectly well behaved at both β values." + `nb.figure` — **"a kick that stays small vs one that grows"**:
    `FD.propagate_error` of a 10⁻¹⁰ random kick (seeded, P10) for 600 steps at β = 0.50 and 0.51, log₁₀ max∣ξ∣ against step
    (blue, rose); *see* "one flat line, one straight rising line"; *read* "the slope of the rising line is log₁₀ 1.04 = 0.017 per
    step"; *change* "…β = 0.45: the line falls".
50. `nb.note` — **N12 [B]** "**The disturbance obeys the same rule (10.18)–(10.19):** $\\xi^n_i=T^n_i-\\overline T^n_i$ (10.18)
    (exact solution of the discrete system minus the computed one) satisfies $\\xi^{n+1}_i=(\\alpha+\\beta)\\xi^n_{i-1}+(1-2\\beta)
    \\xi^n_i+(\\beta-\\alpha)\\xi^n_{i+1}$ (10.19) — because the rule is linear (D04 steps 1–3)." equations (10.18), (10.19).
51. `nb.primer("half-angle identities", "Two trigonometric facts do all the work in this block: $1-\\cos\\theta=2\\sin^2\\frac\\theta2$
    and $\\sin\\theta=2\\sin\\frac\\theta2\\cos\\frac\\theta2$ (so $\\sin^2\\theta=4s(1-s)$ with $s=\\sin^2\\frac\\theta2$). They turn a
    cosine that runs from 1 to −1 into a square that runs from 0 to 1.", code="import numpy as np\nth = np.linspace(0, np.pi, 5)                     # a few angles\nprint(np.allclose(1 - np.cos(th), 2*np.sin(th/2)**2))   # True\ns = np.sin(th/2)**2                                  # the new variable, between 0 and 1\nprint(np.allclose(np.sin(th)**2, 4*s*(1 - s)))       # True")` (**P225**)
52. `nb.md` — **🔁 Reminder:** Euler's formula $e^{\\mathrm i\\theta}=\\cos\\theta+\\mathrm i\\sin\\theta$ (Ch. 1 P45) and complex
    numbers in numpy (`1j`, `np.abs`, `np.angle` — Ch. 6 P153); the modulus $\\lvert a+\\mathrm ib\\rvert^2=a^2+b^2$.
53. `nb.derivation("D04", …)` — Part F D04 (11 steps), ref "10.24". (Uses R04, P45, P153.)
54. `nb.note` — **N13 [B]** "**The stability condition (10.25):** $\\big\\lvert g^{n+1}/g^n\\big\\rvert\\le1$ for **every**
    wavenumber — one bad wave is enough. `FD.max_amplification(scheme, α, β)` scans θ ∈ [0, π]; `FD.is_von_neumann_stable`
    compares with 1." equation (10.25).
55. `nb.derivation("D05", …)` — Part F D05 (8 steps), ref "10.26". (Uses P45, P225.)
56. `nb.note` — **N14 [B]** "**∣G∣² of FTCS (10.26):** $\\big(1-4\\beta\\sin^2\\frac\\theta2\\big)^2+(2\\alpha\\sin\\theta)^2\\le1$,
    $\\theta=k\\pi\\Delta x$ — the real part carries the diffusion, the imaginary part the convection." equation (10.26).
57. `nb.primer("sign of a linear function on an interval", "A straight line f(s) = a + bs is ≤ 0 everywhere on an interval
    if and only if it is ≤ 0 at both ends — a line cannot poke up in the middle. So one inequality 'for all s' becomes two
    inequalities you can read off. (If the interval is open at an end, test the limit there.)", code="import numpy as np\na, b = -0.4, 0.3                         # f(s) = -0.4 + 0.3 s on [0, 1]\ns = np.linspace(0, 1, 11)\nprint((a + b*s <= 0).all(), a <= 0 and a + b <= 0)   # True True: checking the two ends was enough")` (**P226**)
58. `nb.derivation("D06", …)` — Part F D06 (12 steps), ref "10.27". (Uses P225, P226.)
59. `nb.note` — **N15 [B]** "**Noye's region (10.27):** $0\\le4\\alpha^2\\le2\\beta\\le1$ — the book cites Noye (1983) without
    proof; D06 derives it. Two edges: $4\\alpha^2\\le2\\beta$ (enough diffusion to damp the longest waves that convection kicks)
    and $2\\beta\\le1$ (diffusion must not overshoot the shortest wave)." equation (10.27).
60. `nb.derivation("D07", …)` — Part F D07 (7 steps), ref "10.13".
61. `nb.note` — **N17 [B]** "**BTCS is unconditionally stable** (the book: 'it can easily be shown'; D07 shows it):
    $G_{BTCS}=\\big[1+4\\beta\\sin^2\\frac\\theta2+2\\mathrm i\\alpha\\sin\\theta\\big]^{-1}$ has ∣G∣ ≤ 1 for every Δt. At β = 100 the zigzag
    is damped to 1/401 per step." + `nb.code`: `print(abs(FD.amplification_factor(np.pi, 0.0, 100.0, "btcs")))`. *expect:*
    0.0024938.
62. `nb.worked_example("G by hand", "α = 0.1, β = 0.2 (the C02 example), θ = π/2 (a wave four cells long). 1. $e^{-\\mathrm i\\pi/2}
    =-\\mathrm i$, $e^{\\mathrm i\\pi/2}=\\mathrm i$. 2. (10.24): G = 0.3(−i) + 0.6 + 0.1(i) = 0.6 − 0.2i. 3. ∣G∣² = 0.36 + 0.04 = 0.40.
    4. Check with (10.26): (1 − 4 × 0.2 × ½)² + (2 × 0.1 × 1)² = 0.6² + 0.2² = 0.40 ✓ — this wave loses 37 % of its amplitude per
    step. 5. Now β = 0.51, α = 0, θ = π (the zigzag): G = 1 − 4β = −1.04: the zigzag flips sign and grows 4 % per step; from
    10⁻¹⁶ it reaches 1 after ln(10¹⁶)/ln 1.04 ≈ 940 steps. 6. Noye for (0.1, 0.2): 4α² = 0.04 ≤ 2β = 0.4 ≤ 1 ✓ stable.")`
63. `nb.code` — `G = FD.amplification_factor(np.pi/2, 0.1, 0.2, "ftcs")`; print G, `FD.ftcs_amplification_modulus2(np.pi/2,
    0.1, 0.2)`; then a table over (α, β) ∈ {(0.1, 0.2), (0.4, 0.2), (0, 0.5), (0, 0.51), (0.2, 0)} of
    `FD.stability_verdict("ftcs", a, b)` (reason, Gmax, θ_worst in units of π). *expect:* 0.6−0.2j, 0.40; verdicts "stable"
    (Gmax 1.000 at θ = 0), "unstable: 4a^2 > 2b, long waves grow" (θ_worst small), "stable" (Gmax 1.000), "unstable: 2b > 1, the
    zigzag grows" (1.04, θ = π), "unstable: pure convection, every wave grows". *explain:* 1. G at one θ; 2. its modulus two ways;
    3. the verdict function scans θ and names the edge that fails.
64. `nb.check_agree` — **from scratch (curation §7):** G(θ) on 2001 angles by complex arithmetic from (10.24) and ∣G∣² by
    `np.abs(G)**2` vs the closed form (10.26): `assert np.allclose(np.abs(G)**2, FD.ftcs_amplification_modulus2(th, a, b))`
    for 200 random (α, β) (seeded); then Noye by brute force: `stable_scan = np.abs(G).max() <= 1 + 1e-12` vs
    `FD.ftcs_stable(a, b)` on a 60 × 60 grid of (α, β) ∈ [0, 0.6] × [0, 0.6] (points within 10⁻³ of an edge skipped): `assert
    (scan == closed).all()`.
65. `nb.figure` — **"One bad wave is enough"** (10 × 3.6 in, three panels): (a) complex plane: the unit circle (muted), G(θ)
    for θ ∈ [0, π] at (α, β) = (0.1, 0.2) (orange, inside), (0.4, 0.2) (rose, leaves the circle near θ = 0), (0, 0.51) (rose,
    leaves at G = −1.04) and BTCS at (0.1, 0.2) (blue); (b) ∣G(θ)∣ vs θ/π for the same cases with the line ∣G∣ = 1; (c) the (α, β)
    plane with Noye's region 4α² ≤ 2β ≤ 1 shaded (teal), its two edges labelled, the five table points as dots. *see:* "curves
    inside or poking out of the unit circle; a lens-shaped region"; *read:* "a curve leaving the circle near θ = 0 means long
    waves grow (the 4α² ≤ 2β edge); leaving at θ = π means the zigzag grows (the 2β ≤ 1 edge)"; *change:* "…D = 0 (β = 0): the
    region shrinks to the single point α = 0 — pure convection is never stable".
66. `nb.plotly` — **F2** `slider_figure` over β (30 values 0.05 … 0.6) at fixed α = 0.15: traces ∣G(θ)∣ (orange) and the
    line 1 (muted); a second trace set = the region point. *explain:* "the stability edge on the page, no kernel needed".
67. `nb.live` — `live(lambda alpha, beta, scheme: <plot ∣G(θ)∣ and the verdict>, alpha=(0.0, 0.6, 0.01), beta=(0.0, 0.8,
    0.01), scheme=["ftcs", "btcs", "upwind", "cn"])` — paired with F2 (the page shows F2).
68. `nb.explainer("von_neumann_amplification", heading="Why does β = 0.51 blow up when 0.50 does not?", why="Dragging the
    point (α, β) across the edge of Noye's region, you watch the G(θ) curve cross the unit circle and, on the same clock, a
    10⁻¹⁰ kick grow into a zigzag — the causal chain from a number in a formula to garbage on the screen.", tries=["Press the
    'β = 0.51' preset and play: which wave grows first? (the status names it)", "Set β = 0 and drag α: every curve leaves the
    circle — pure convection is hopeless for FTCS.", "Switch the scheme to BTCS and set β = 100: the curve shrinks towards the
    origin."])`
69. `nb.md` — **What would change if…** "…there is no diffusion at all (a pollutant in a fast river, a tracer in the upper
    ocean)? FTCS is useless — but taking the convective difference from the upstream side fixes it, at a price. C05."

#### C05 — Upwind differencing and the CFL condition (closed by the Lax equivalence theorem)
70. `nb.core("C05", "Upwinding and the CFL condition $u\\frac{\\Delta t}{\\Delta x}\\le1$ (10.30)", question="Why must the time
    step shrink whenever the grid is refined — and what does the flow speed have to do with it?")`
71. `nb.md` — **The problem in plain words:** "A weather model with 25 km grid boxes must never let the fastest signal — a
    sound or gravity wave at 200–300 m/s — cross more than one box per time step, so its time step is at most about a
    minute and a half. Halve the box size and the step must halve too: that is why doubling a climate model's resolution
    costs about eight times the computing (twice the columns in each direction and twice the steps). Where does such a rule
    come from?"
72. `nb.md` — **The idea** (ASCII x–t diagram):
    ```
    t_{n+1}            ● new value at x_i
                      ╱                    the characteristic x − ut = const through the new point
                     ╱  slope 1/u           lands at x_i − uΔt at the old level
    t_n     ●───────◆───────●             ◆ = the foot of the characteristic
          x_{i−1}         x_i
    CFL:  the foot must lie between the stencil's points  ⇔  0 ≤ uΔt/Δx ≤ 1
    ```
    "Pure convection copies T along lines x − ut = const. The update at x_i may only use numbers from the side the flow comes
    from (upwind), and those numbers must bracket the foot of the line."
73. `nb.primer("domain of dependence and characteristics", "For $T_t+uT_x=0$ the value at (x, t) is the value that was at
    (x − uΔt, t − Δt): information travels along the lines x − ut = const (the characteristics of Ch. 7 P174). The *domain of
    dependence* of a point is the set of earlier data that can influence it — for this PDE a single point upstream; for a
    stencil it is the stencil's points. A scheme can only be right if its domain contains the PDE's.", code="u, dt, x = 2.0, 0.1, 1.0      # speed [m/s], step [s], where we want the new value [m]\nfoot = x - u*dt                   # the value at x now came from here one step ago\nprint(foot)                       # 0.8 m: with dx = 0.1 m that is 2 cells upstream — outside a 1-cell stencil (C = 2)")`
    (**P227**)
74. `nb.note` — **N16 [B]** "**First-order upwind (10.29):** $T^{n+1}_i=T^n_i-2\\alpha(T^n_i-T^n_{i-1})$ — the convective
    difference is taken backwards, from the upstream neighbour (u > 0), and 2α = uΔt/Δx = C is the **Courant number**. > ⚠️ **The
    book prints (10.29)–(10.30) without the sign of u; for u < 0 the upstream neighbour is $i+1$ and the condition is
    $\\lvert u\\rvert\\Delta t/\\Delta x\\le1$.** `scheme='upwind'` takes the side from sign(u); `scheme='upwind_printed'` always
    uses i − 1 and blows up for u < 0 (a test checks it)." equation (10.29).
75. `nb.derivation("D08", …)` — Part F D08 (9 steps), ref "10.30". (Uses D04 moves, P45, P225, P227.)
76. `nb.worked_example("one upwind step, and the CFL edge", "u = 1 m/s, Δx = 0.1 m, Δt = 0.05 s: C = 0.5. 1. Start [0, 0, 1,
    0, 0]. 2. Node 2: 1 − 0.5(1 − 0) = 0.5; node 3: 0 − 0.5(0 − 1) = 0.5; others 0. 3. New [0, 0, 0.5, 0.5, 0]: the spike moved half a
    cell (as it should, uΔt = 0.05 m) but is now spread over two cells — upwind smears. 4. The zigzag θ = π: ∣G∣² = 1 −
    2C(1 − C)(1 − cos π) = 1 − 4 × 0.25 = 0 — killed in one step. 5. C = 1: G = e^{−iθ}, an exact shift by one cell. 6. C = 1.1:
    ∣G(π)∣² = 1 − 2(1.1)(−0.1)(2) = 1.44, ∣G∣ = 1.2 — the zigzag grows 20 % per step.")`
77. `nb.code` — `for C in (0.5, 1.0, 1.05): r = FD.advect_periodic("square", C=C, n_cells=100, n_rev=1.0, scheme="upwind");
    print(C, r["amplitude_ratio"], r["rms_error"], r["blew_up"])`; `print(FD.stability_verdict("upwind", 0.525, 0.0)["reason"])`;
    `print(FD.numerical_diffusivity(1.0, 0.01, scheme="upwind", C=0.5))`. *expect:* C = 0.5: amplitude ≈ 0.5–0.6 of the square
    after one revolution (smeared), not blown; C = 1.0: amplitude 1.000, rms error ≈ 1e-15 (exact shift); C = 1.05: blew_up
    True; reason "unstable: C > 1, the characteristic leaves the stencil"; D_num = 0.0025 m²/s. *explain:* 1. one revolution =
    N/C steps; 2. the verdict; 3. the numerical diffusivity ∣u∣Δx(1 − C)/2 that upwind adds (ours, from the modified equation
    — C09 does the steady version).
78. `nb.check_agree` — **from scratch (curation §7):** an upwind loop over nodes for u > 0 (`for i in range(N): new[i] =
    T[i] - C*(T[i] - T[i-1])`, Python's `T[-1]` wraps round — periodic), then `assert np.allclose(new,
    FD.transport_1d_step(T, C/2, 0.0, "upwind", periodic=True))`; and u < 0: `assert not np.allclose(...)` against
    `"upwind_printed"` after 50 steps (it has grown by > 10³).
79. `nb.md` — **Climate hook (N-less, our note):** "**Every explicit weather and ocean model lives under this rule.** The
    fastest signal in a shallow-water or hydrostatic model is the gravity wave, speed $\\sqrt{g_0H}$ (the shallow limit kH ≪ 1 of Ch. 7's $\\omega^2=gk\\tanh kH$ (7.28)),
    `ch10.cfl_time_step(c, dx)`: ocean, H = 4000 m, c = 198 m/s: Δx = 100 km → Δt ≤ 505 s; Δx = 25 km → 126 s.
    Atmosphere, Lamb-wave speed ≈ 320 m/s, Δx = 25 km → 78 s. This is why models treat the fastest waves implicitly (C02's
    BTCS idea) or split them off (C11)." + `nb.code` printing the three numbers with `WV.phase_speed(1e-7, 4000.0)`. *expect:*
    198.1 m/s; 504.8 s, 126.2 s, 78.1 s.
80. `nb.primer("well-posed problem", "A problem is well posed (Hadamard) when a solution exists, is unique, and depends
    continuously on the data — a small change in the start or boundary values gives a small change in the answer. The heat
    and advection equations with sensible boundary conditions are; running the heat equation backwards in time is not (tiny
    wiggles explode). The Lax theorem below assumes it.", code="import numpy as np\nk = np.array([1, 10, 100])          # wavenumbers of a small wiggle\nprint(np.exp(-k**2*0.01))           # forward heat equation, t = 0.01: wiggles shrink (well posed)\nprint(np.exp(+k**2*0.01))           # backward: the k = 100 wiggle grows by e^100 (ill posed)")` (**P228**)
81. `nb.note` — **N18 [B]** "**The Lax equivalence theorem** (Richtmyer & Morton 1967, our words): *for a well-posed linear
    initial-value problem and a consistent scheme, stability is necessary and sufficient for convergence.* So the three
    properties of C03 collapse into two checks we can do on paper: consistency (Taylor, C03) and stability (von Neumann, C04–
    C05). The book gives no proof; we show it at work." + `nb.animation` — **A1** (frames player, 40 frames; FAST 24): the
    heated rod (N113, our T_w = 1) by FTCS at β = 0.50 and 0.51 side by side (profiles over the exact series, and log₁₀ max
    error vs step underneath); the right panel grows a zigzag after ≈ 900 steps. Built with `FD.lax_demo` +
    `ch10.rod_heating_exact`. *see/read/change:* "same rule, same consistency; the 2 % larger step crosses ∣G∣ = 1 and the error
    stops converging"; "…upwind at C = 0.9 vs 1.1: the same story at the CFL edge (a second code line prints both final
    errors)".
82. `nb.note` — **N113 [B]** "**An exact solution to test against — the heated rod (Exercise 10.3 uses this form):** a rod
    0 ≤ x ≤ 1 at temperature 0, whose ends are suddenly held at $T_w$, has
    $T=T_w-\\sum_{m=1}^{\\infty}\\frac{4T_w}{(2m-1)\\pi}\\sin[(2m-1)\\pi x]\\,e^{-D(2m-1)^2\\pi^2t}$ — equation (10.199) is this
    series (a Fourier sine series, as for Couette start-up in Ch. 8). We use our own $T_w=1$ and grid; the exercise's inputs
    and answers are not reproduced." equation (the series above) + `nb.code`: `ch10.rod_heating_exact(np.array([0.0, 0.25,
    0.5]), 0.05, 1.0, 1.0)`. *expect:* [1 (the wall), 0.4468, 0.2277] (the cell prints; the assert checks the PDE residual by
    finite differences < 1e-6 and the walls = T_w).
83. `nb.primer("Péclet number (global R and cell R_cell)", "Advection against diffusion, as one ratio: over a length ℓ,
    advection moves heat at rate ~ uT/ℓ, diffusion at ~ DT/ℓ²; their ratio is uℓ/D. With ℓ = L (the domain) it is the global
    Péclet number R = uL/D (10.87) — the Reynolds number's twin for a scalar. With ℓ = Δx it is the cell Péclet number
    $R_{cell}=u\\Delta x/D$ (10.31): how advective one grid cell is.", code="u, L, dx, D = 1.0, 1.0, 0.01, 0.005   # speed [m/s], length [m], cell [m], diffusivity [m^2/s]\nprint(u*L/D, u*dx/D)                     # R = 200 (strongly advective), R_cell = 2 (the edge of C09)")` (**P229**)
84. `nb.note` — **N19 [B]** "**The cell Péclet condition (10.31):** $R_{cell}=u\\frac{\\Delta x}{D}\\le2$ to avoid wiggles next to a
    thin layer; with u = 1 m/s, Δx = 0.01 m, D = 0.005 m²/s, R_cell = 2 exactly. The reason — the discrete solution changes
    sign — is worked out in C09 (D16)." equation (10.31) + `nb.code`: `FD.cell_peclet(1.0, 0.01, 0.005)`. *expect:* 2.0.
85. `nb.figure` — **"CFL is about information"** (10 × 3.4 in, three panels): (a)–(c) `characteristic_diagram(ax, C)` for
    C = 0.5, 1.0, 1.2: the two-point upwind stencil, the characteristic through the new point, its foot (◆) inside, on, or
    outside the stencil (teal, teal, rose); underneath each, the square pulse after one revolution (exact muted, upwind teal;
    C = 1.2 a rose "blew up at step n" label from `advect_periodic`). *see:* "the foot slides left as C grows and leaves the
    stencil at C = 1"; *read:* "inside: stable (and smeared); on the edge: exact; outside: unstable"; *change:* "…u < 0 with the
    printed (10.29): the foot is on the *other* side — unstable for every C (R11)".
86. `nb.plotly` — **F3** `slider_figure` over C (24 values 0.1 … 1.0): upwind (teal) and MacCormack (purple) after one
    revolution of the square pulse vs exact (muted) — `FD.advect_periodic`. *explain:* "smearing shrinks as C → 1 for upwind;
    MacCormack keeps the height but ripples (C10)".
87. `nb.explainer("upwind_cfl_advection", heading="Why must the time step shrink with the grid?", why="The Courant slider
    moves the characteristic across the stencil's edge at exactly C = 1 while, on the same clock, the pulse stays clean,
    smears, ripples or explodes; toggling upwind, MacCormack and FTCS separates numerical diffusion from dispersion.",
    tries=["Press 'C = 1 exact shift' and play: nothing smears.", "Drag C to 1.05 and watch the characteristic's foot leave the
    stencil — then the blow-up.", "Compare upwind and MacCormack at C = 0.8 on the square pulse: which loses height, which
    ripples?"])`
88. `nb.md` — **What would change if…** "…the problem were steady and both effects mattered (heat carried to a hot wall and
    conducted away)? The time step disappears, but the grid still has to resolve the thin layer where diffusion wins: C09.
    First, §10.3 builds the second family of methods — finite elements — on the same model problem."

---

### A.3 §10.3 Finite-Element Method — C06 · C07 · C08
1. `nb.section("10.3", "Finite-Element Method", intro="**What is this section about?** Finite differences ask the equation
   to hold at grid points. Finite elements ask something weaker and more flexible: that the equation hold *on average*,
   weighted by every member of a family of test functions. Written that way (the weak form, C06), the problem can be solved
   with simple building blocks — tent-shaped 'hat' functions — whose weights become the unknowns (Galerkin, C07); the
   integrals are computed one small element at a time and added up (assembly, C08). On a uniform 1-D mesh the result is almost
   the FD scheme you already know; the payoff comes in 2-D and 3-D, where elements can follow a curved body (the cylinder of
   §10.5, a coastline, an artery).")`

#### C06 — The weak form; essential and natural boundary conditions
2. `nb.core("C06", "The weak form $\\int_0^LT_tw\\,dx+u\\int_0^LT_xw\\,dx+D\\int_0^LT_xw_x\\,dx=Dq\\,w(L)$ (10.36)",
   question="How can an equation 'hold on average' — and why does that make one boundary condition automatic?")`
3. `nb.md` — **The problem in plain words:** "A finite-element approximation is made of straight pieces joined at kinks. Its
   second derivative is zero on every piece and infinite at the kinks — so it can never satisfy $D\\,\\partial^2T/\\partial
   x^2$ point by point. We need a version of (10.1) that only asks for *first* derivatives. The trick: multiply the equation
   by a smooth-enough weight w, integrate over the rod, and move one derivative from T onto w. Nothing is lost — if the
   averaged equation holds for every weight, the pointwise one holds too — and the insulated/flux end condition drops into the
   equation by itself."
4. `nb.md` — **The idea** (ASCII, three moves):
   ```
   strong:  T_t + u T_x − D T_xx = 0  at every x          (10.1)
      × w(x), ∫₀ᴸ dx          (w = 0 where T is prescribed: x = 0)
   ∫ T_t w + u ∫ T_x w − D ∫ T_xx w = 0                    (10.34)
      ∫ T_xx w = [T_x w]₀ᴸ − ∫ T_x w_x      (integration by parts)
   ∫ T_t w + u ∫ T_x w + D ∫ T_x w_x = D T_x(L) w(L) = D q w(L)      (10.35)–(10.36)
            only first derivatives left ──┘          └── the Neumann datum q enters here: "natural"
   ```
5. `nb.note` — **N20 [C]** "**Strong (classical) form** = the PDE (10.1) with its boundary conditions (10.2) and initial condition
   (10.3) — 'strong' because it asks for the equation at every point."
6. `nb.primer("test functions and the spaces H¹, S and V", "A *test function* w(x) is a weight we multiply the equation by; a
   weak statement says 'for every w in a family'. The family here is $H^1$: functions whose slope is square-integrable,
   $\\int_0^L(\\partial w/\\partial x)^2dx<\\infty$ — kinks allowed, jumps not ('finite slope energy'). The *trial* space $\\mathcal S$
   = the $H^1$ functions with the Dirichlet value built in, T(0) = g (10.32); the *test* space V = the same with w(0) = 0 (10.33).",
   code="import numpy as np\nx = np.linspace(0, 1, 100001)\nhat = np.maximum(0, 1 - np.abs(x - 0.5)/0.25)       # a tent: kinks at 0.25, 0.5, 0.75\nslope = np.gradient(hat, x)                            # ±4 on the sides, 0 elsewhere\nprint(np.trapezoid(slope**2, x))                       # ≈ 8 (16 × 0.5): finite, so the tent is in H^1")`
   (**P230**)
7. `nb.note` — **N21 [B]** "**The two spaces (10.32)–(10.33):** $\\mathcal S=\\{T\\,\\vert\\,T\\in H^1,\\ T(0)=g\\}$ and
   $V=\\{w\\,\\vert\\,w\\in H^1,\\ w(0)=0\\}$ (the book writes S; we write $\\mathcal S$ so it is not the Strouhal number). The test
   functions vanish exactly where T is prescribed — there we already know the answer and must not 'test' it." equations
   (10.32), (10.33).
8. `nb.md` — **🔁 Reminder:** integration by parts $\\int_a^b f\\,g'\\,dx=[fg]_a^b-\\int_a^b f'g\\,dx$ (Ch. 9 P218a), the product
   rule read backwards (Ch. 1 P38).
9. `nb.derivation("D09", …)` — Part F D09 (9 steps), ref "10.36". (Uses P218a, P38, P230.)
10. `nb.worked_example("integration by parts with numbers", "Take T = x² and w = x on [0, 1] (w(0) = 0 ✓). 1. Left side of the
    move: $\\int_0^1T_{xx}w\\,dx=\\int_0^12x\\,dx=1$. 2. Right side: $[T_xw]_0^1-\\int_0^1T_xw_x\\,dx=(2\\cdot1\\cdot1-0)-\\int_0^12x\\cdot
    1\\,dx=2-1=1$ ✓. 3. The boundary term came only from x = 1, because w(0) = 0 killed the other end — exactly how the Dirichlet
    end disappears from (10.35) and the Neumann end survives as $Dq\\,w(L)$.")`
11. `nb.code` — `ws = ch10.weak_form_sympy()`; print `ws["residual"]` (0) and the three displayed lines; then the steady
    exact solution of $uT_x=DT_{xx}$ with T(0) = 0, $T_x(1)=q$ (u = 1, D = 0.25, q = 1: $T=\\frac{qD}{u}e^{-u/D}(e^{ux/D}-1)$) tested
    against hat functions: `[FEM1.weak_residual(T_fn, lambda x, A=A: FEM1.hat(x, nodes, A), 1.0, 0.25, 1.0) for A in
    range(1, 9)]`. *expect:* sympy residual 0; the eight weak residuals ≈ 0 (< 1e-10): the exact solution passes every test.
    *explain:* 1. sympy repeats D09 on polynomials; 2. `weak_residual` evaluates (10.36) (steady) by Gauss quadrature; 3. every
    hat test gives zero.
12. `nb.primer("fundamental lemma of the calculus of variations", "If $\\int_0^Lf(x)\\,w(x)\\,dx=0$ for *every* smooth w that
    vanishes at the ends, then f = 0 everywhere (for continuous f). Reason: choose w = f times a bump that is zero at the ends;
    the integral becomes $\\int f^2\\times$bump, which is zero only if f is. A weak statement 'for all w' therefore pins down f
    point by point.", code="import numpy as np\nx = np.linspace(0, 1, 2001)\nf = np.sin(3*x) - 0.5                        # some nonzero f\nw = f*x*(1 - x)                               # the clever test function: f times a bump\nprint(np.trapezoid(f*w, x))                   # 0.0218 > 0: this w exposes f ≠ 0")` (**P231**)
13. `nb.derivation("D10", …)` — Part F D10 (9 steps), ref "10.38". (Uses P218a, P231.)
14. `nb.note` — **N22 [B]** "**Weak ⇒ strong (10.37)–(10.38):** reversing the integration by parts gives
    $\\int_0^L(T_t+uT_x-DT_{xx})w\\,dx+D[T_x(L)-q]w(L)=0$ (10.37); holding for every w in V it forces the PDE on (0, L) *and*
    $T_x(L)=q$ (10.38). **Essential** boundary condition (Dirichlet): built into the trial space $\\mathcal S$ — the weak form never
    sees it. **Natural** (Neumann): comes out of the weak form by itself." equations (10.37), (10.38) + `nb.code`:
    `ch10.weak_to_strong_sympy()["pde"], ["natural_bc"]`.
15. `nb.figure` — **"The natural condition is learned, not imposed"** (8 × 3.4 in, two panels): (a) the steady FE solution
    (`FEM1.solve_steady(np.linspace(0, 1, n + 1), 1.0, 0.25, g=0.0, q=1.0)`) for n = 4, 8, 16 (teal dots, lines) on the exact
    curve (muted dashed) with the prescribed slope q = 1 drawn as a short orange tangent at x = 1; (b) ∣last-element slope − q∣
    vs h on log–log (orange dots, slope-1 guide). *see:* "the value at x = 0 is exact on every mesh; the slope at x = 1 approaches
    q"; *read:* "essential: imposed exactly; natural: satisfied in the limit, error ∝ h"; *change:* "…we prescribed T(1) instead:
    it would have to be built into the trial space as well (the `T_L=` option), and both ends would be exact".
16. `nb.md` — **What would change if…** "…we only allowed a *finite* set of test and trial functions? The weak form becomes a
    finite system of equations — the Galerkin method, C07."

#### C07 — Galerkin with hat functions: M ḋ + K d = F
17. `nb.core("C07", "Galerkin's method with hat functions: $\\mathbf M\\dot{\\mathbf d}+\\mathbf K\\mathbf d=\\mathbf F$ (10.58)",
    question="If the unknown is a whole function, how does it become a handful of numbers the computer can solve for?")`
18. `nb.md` — **The problem in plain words:** "Build T out of Lego: tent-shaped pieces, one centred on each node, each as tall as
    the temperature there. Adding the tents gives a broken-line profile through the nodal values. The heights are the unknowns.
    To find them, demand the weak form not for *every* test function (impossible — infinitely many) but for every tent. One tent,
    one equation: as many equations as unknowns."
19. `nb.md` — **The idea** (ASCII):
    ```
    N₀    N₁    N₂    N₃    N₄        hats: 1 at their own node, 0 at every other (N_A(x_B) = δ_AB)
    /\    /\    /\    /\    /|
    T^h(x) = g·N₀ + d₁N₁ + d₂N₂ + d₃N₃ + d₄N₄     a broken line through (x_A, d_A)
    test with w = N_A for A = 1…n   ⇒   n equations   M ḋ + K d = F
    ```
20. `nb.note` — **N23 [B]** "**The discrete weak problem (10.39)–(10.42):** find $T^h\\in\\mathcal S^h$ with
    $\\int_0^LT^h_tw^h\\,dx+u\\int_0^LT^h_xw^h\\,dx+D\\int_0^LT^h_xw^h_x\\,dx=Dq\\,w^h(L)$ for all $w^h\\in V^h$ (10.39). Split off the
    boundary value: $v^h=T^h-g^h$ (10.40) with $g^h(0)=g$, so $v^h$ and $w^h$ live in the *same* space; then
    $\\int_0^Lv^h_tw^h\\,dx+a(w^h,v^h)=Dq\\,w^h(L)-a(w^h,g^h)$ (10.41) with the bilinear form
    $a(w,v)=u\\int_0^Lv_xw\\,dx+D\\int_0^Lv_xw_x\\,dx$ (10.42)." equations (10.39)–(10.42) + `nb.code`: `FEM1.bilinear_form` on two
    hats vs the assembled K entry (`assert np.isclose`).
21. `nb.note` — **N24 [C]** "**Galerkin vs Petrov–Galerkin:** Galerkin = test functions from the same space as the solution;
    Petrov–Galerkin = a different test space. Pointer: SUPG (N50) is Petrov–Galerkin — it tilts the test functions upstream."
22. `nb.primer("bilinear form a(w, v)", "A function of two functions that is linear in each slot separately: $a(c_1w_1+c_2w_2,
    v)=c_1a(w_1,v)+c_2a(w_2,v)$, and the same in v. So sums and constants can be pulled out of it, just like out of an integral —
    the move that turns (10.50) into (10.52). Here a(w, v) is *not* symmetric (the convective part differentiates only v).",
    code="import numpy as np\nx = np.linspace(0, 1, 2001); u, D = 1.0, 0.1\na = lambda w, v: np.trapezoid(u*np.gradient(v, x)*w + D*np.gradient(v, x)*np.gradient(w, x), x)\nw1, w2, v = x, x**2, np.sin(x)\nprint(np.isclose(a(2*w1 + 3*w2, v), 2*a(w1, v) + 3*a(w2, v)))   # True: linear in the first slot\nprint(np.isclose(a(w1, v), a(v, w1)))                             # False: not symmetric")`
    (**P232**)
23. `nb.note` — **N25 [B], N26 [B], N27 [B]** "**The basis (10.43)–(10.49):** $N_A(x)$, A = 1…n, with $N_A(0)=0$ (10.43)–(10.44), so
    every test function is $w^h=\\sum_{A=1}^nc_AN_A(x)$ (10.45); one extra $N_0$ with $N_0(0)=1$ (10.46) carries the boundary value,
    $g^h=gN_0$ (10.47); the solution is $v^h=\\sum_Ad_A(t)N_A$ (10.48) and $T^h=\\sum_{A=1}^nd_A(t)N_A(x)+gN_0(x)$ (10.49)."
    equations (10.45), (10.49) + `nb.code`: `B = FEM1.hat_basis(nodes, x)`; print `B.sum(axis=1).min(), .max()` and `B[nodes
    index]` (the identity). *expect:* 1.0, 1.0 (partition of unity), identity matrix.
24. `nb.primer("arbitrary coefficients: every bracket is zero", "If $\\sum_{A=1}^nc_AG_A=0$ must hold for *every* choice of the
    numbers c_A, then each $G_A=0$: choose c = (1, 0, …, 0) to get $G_1=0$, then (0, 1, 0, …) for $G_2$, and so on. It is how one
    weak statement 'for all w' becomes n separate equations.", code="import numpy as np\nG = np.array([0.3, -0.1, 0.0])           # suppose sum(c*G) were 0 for every c ...\nfor c in np.eye(3):                       # ... then try the unit vectors\n    print(c, c @ G)                       # 0.3 ≠ 0 exposes G_1: so all G_A must vanish")`
    (**P233**)
25. `nb.derivation("D11", …)` — Part F D11 (13 steps), ref "10.58". (Uses P232, P233, N23–N27.)
26. `nb.note` — **N28 [B]** "**The middle steps (10.50)–(10.53):** substituting (10.45), (10.48) into (10.41) gives (10.50);
    rearranged, $\\sum_Ac_AG_A=0$ (10.51) with $G_A=\\sum_B\\dot d_B\\int_0^LN_AN_B\\,dx+\\sum_Bd_B\\,a(N_A,N_B)-DqN_A(L)+g\\,a(N_A,N_0)$
    (10.52); arbitrary $c_A$ ⇒ $\\sum_B\\dot d_B\\int_0^LN_BN_A\\,dx+\\sum_Bd_B\\,a(N_A,N_B)=DqN_A(L)-g\\,a(N_A,N_0)$ (10.53): n ODEs."
    equations (10.51)–(10.53) + `nb.code`: `ch10.galerkin_equations_sympy(3)` printing M, K, F for three hats.
27. `nb.note` — **N29 [C]** "**Time next:** (10.58) is a system of ODEs — integrate it with an ODE solver (method of lines; `solve_ivp`,
    Ch. 1 P31) or with finite differences in time (the θ-scheme). Both are options of `FEM1.solve_transport(method=…)`."
28. `nb.note` — **N30 [B]** "**Hat functions (10.59)–(10.61):** $N_A=\\frac{x-x_{A-1}}{x_A-x_{A-1}}$ on $[x_{A-1},x_A)$,
    $\\frac{x_{A+1}-x}{x_{A+1}-x_A}$ on $[x_A,x_{A+1}]$, 0 elsewhere (10.59); half-hats at the ends, $N_n$ (10.60) and $N_0$ (10.61).
    They are *compact*: each is nonzero on two elements only, and $N_A(x_B)=\\delta_{AB}$." equation (10.59) + `nb.figure` —
    **Fig. 10.2 remade** (`hat_functions(ax, n_el=6)`: N₀ … N₆ muted, one interior hat N₃ highlighted amber; underneath, the
    weighted hats $d_AN_A$ (teal) of a sample profile and their sum $T^h$ (blue) through the nodal dots). *see:* "tents that
    overlap only with their neighbours; a broken line made of them"; *read:* "the height of each tent is the value at its node
    (N31)"; *change:* "…quadratic pieces (P2, C13): three nodes per element, curved tents".
29. `nb.note` — **N31 [B]** "**The weights are nodal values (10.62):** $d_A=T^h(x_A)=T_A$ — because only $N_A$ is nonzero at
    $x_A$. Read d₂ off the plot above." equation (10.62).
30. `nb.derivation("D12", …)` — Part F D12 (10 steps), ref "10.63". (Uses N30, D11.)
31. `nb.note` — **N32 [B]** "**The interior row (10.63):** $\\frac d{dt}\\Big(\\frac{T_{A-1}}6+\\frac{2T_A}3+\\frac{T_{A+1}}6\\Big)+
    \\frac u{2h}(T_{A+1}-T_{A-1})-\\frac D{h^2}(T_{A-1}-2T_A+T_{A+1})=0$ — the centred convection and diffusion stencils of C01,
    but the time derivative acts on a ⅙–⅔–⅙ average of three nodes (the **consistent mass**). Replacing the average by T_A
    alone ('lumping') gives exactly the semi-discrete FD scheme." equation (10.63).
32. `nb.note` — **N33 [B]** "**'Galerkin FE is equivalent to an FD method'** (§10.3) — true, with a caveat: for *uniform linear*
    elements the steady equations are identical (the mass matrix drops out when ∂/∂t = 0); in time the consistent mass differs;
    on nonuniform or higher-order meshes the two differ. FE's real advantage is geometry." + `nb.code`: steady FE vs centred FD,
    `np.max(np.abs(FEM1.solve_steady(nodes, 1.0, 0.25, T_L=1.0)["T"] - FD.steady_cd_fd(4, 4.0)["T"]))`. *expect:* < 1e-13.
33. `nb.worked_example("four elements, by hand", "L = 1, n = 4 elements, h = 0.25, u = 1 m/s, D = 0.25 m²/s (R_cell = uh/D = 1),
    steady, T(0) = 0, T(1) = 1. 1. Interior row of (10.63) × h: mass (h/6, 2h/3, h/6) = (0.0417, 0.1667, 0.0417) (unused when
    steady). 2. Convection: (u/2)(T_{A+1} − T_{A−1}) → (−0.5, 0, 0.5). 3. Diffusion: −(D/h)(T_{A−1} − 2T_A + T_{A+1}) → (−1, 2, −1).
    4. Sum: (−1.5, 2, −0.5). 5. Three equations for T₁, T₂, T₃: 2T₁ − 0.5T₂ = 0; −1.5T₁ + 2T₂ − 0.5T₃ = 0; −1.5T₂ + 2T₃ = 0.5 (the
    known T₄ = 1 moved right). 6. Solve: T = (0.025, 0.1, 0.325). 7. Exact (10.86) at x = 0.25, 0.5, 0.75: 0.0321, 0.1192, 0.3561 —
    close, with R_cell = 1 < 2 no wiggles.")`
34. `nb.code` — `nodes = np.linspace(0, 1, 5)`; `M, K, F = FEM1.assemble_1d(nodes, 1.0, 0.25, g=0.0, dirichlet_right=1.0)`;
    print `M.toarray()`, `K.toarray()`, `F`; `print(FEM1.interior_stencil(0.25, 1.0, 0.25))`; `print(FEM1.solve_steady(nodes, 1.0,
    0.25, T_L=1.0)["T"])`; `print(FD.steady_cd_exact(nodes, 4.0))`. *expect:* K rows (2, −0.5, 0), (−1.5, 2, −0.5), (0, −1.5, 2); F =
    (0, 0, 0.5); M rows (h/6)(4, 1, 0)…; stiff (−1.5, 2.0, −0.5); T = [0, 0.025, 0.1, 0.325, 1]; exact [0, 0.03206, 0.1192,
    0.3561, 1]. *explain:* 1. assembly (C08 opens the box); 2. the interior row is (10.63) × h; 3. the steady solve.
35. `nb.figure` — **"FE builds the answer from tents"** (8 × 3.4 in, two panels): (a) the steady FE solution for n = 4 as
    weighted hats (teal) summing to Tʰ (blue) on the exact (10.86) curve (muted dashed), R = 4; (b) a transient: an initial
    step advected-diffused on n = 20 elements with consistent (blue) and lumped (orange) mass vs a fine reference (muted) at
    t = 0.3 (`FEM1.solve_transport(..., lumped=…)`). *see:* "a broken line through the nodal dots; two transient curves that
    differ slightly"; *read:* "steady: FE = FD exactly; unsteady: the consistent mass is a little less diffusive (sharper
    front)"; *change:* "…nonuniform elements (smaller near x = 1): FE still assembles the same way; FD would need new stencils".
36. `nb.plotly` — **F4** `slider_figure` over n (2, 3, 4, 6, 8, 12, 16, 24, 32): FE steady solution with its hats (teal),
    centred FD dots (orange) and exact (muted), R = 10 (wiggles for n ≤ 4, gone by n = 8). *explain:* "FE = FD on uniform meshes
    for every n; both converge".
37. `nb.md` — **What would change if…** "…we computed each integral $M_{AB}$, $K_{AB}$ over the whole rod? Most are zero, and the
    nonzero ones are all alike. C08 computes one small block per element and adds them up."

#### C08 — Element matrices and assembly
38. `nb.core("C08", "Element matrices and assembly: $\\mathbf m^e=\\frac h6\\begin{pmatrix}2&1\\\\1&2\\end{pmatrix}$,
    $\\mathbf k^e=\\frac u2\\begin{pmatrix}-1&1\\\\-1&1\\end{pmatrix}+\\frac Dh\\begin{pmatrix}1&-1\\\\-1&1\\end{pmatrix}$ from
    (10.74)–(10.77)", question="How is a finite-element matrix actually built — and why does that make FE good at complicated
    shapes?")`
39. `nb.md` — **The problem in plain words:** "A hat function touches only two elements, so the integral $\\int N_AN_B\\,dx$ is
    zero unless A and B are neighbours, and on each element only two hats are alive. Instead of looping over pairs of hats
    (mostly wasted work), loop over elements: compute a small 2 × 2 block on each, then add it into the big matrix at the rows and
    columns of that element's two nodes. The same loop works for triangles on an unstructured mesh around a cylinder, which is
    why FE handles complex geometry."
40. `nb.md` — **The idea** (ASCII):
    ```
    element e = [x_{e−1}, x_e]      local nodes a = 1, 2  ↔  global nodes A = e − 1, e       (10.78)
             ┌ k₁₁ k₁₂ ┐                          column: e−1   e
    k^e  =   └ k₂₁ k₂₂ ┘     scatter-add  ───►   row e−1  [ +k₁₁  +k₁₂ ]
                                                  row e    [ +k₂₁  +k₂₂ ]
    every interior node belongs to two elements, so its row collects one piece from each
    ```
41. `nb.note` — **N34 [B]** "**Element point of view, Fig. 10.3 (10.64)–(10.66):** every element is a stretched copy of the
    *parent element* ξ ∈ [−1, 1], with shapes $N_1(\\xi)=\\frac12(1-\\xi)$, $N_2(\\xi)=\\frac12(1+\\xi)$ (10.64), the map
    $x(\\xi)=N_1(\\xi)x^e_1+N_2(\\xi)x^e_2=\\frac12[(x_A-x_{A-1})\\xi+x_A+x_{A-1}]$ (10.65) and its inverse
    $\\xi(x)=\\frac{2x-x_A-x_{A-1}}{x_A-x_{A-1}}$ (10.66)." equations (10.64)–(10.66) + `nb.figure` — **Fig. 10.3 remade**
    (`element_map(ax)`: the element $[x_{A-1},x_A]$ with $N_{A-1}$, $N_A$ and the parent element with $N_1$, $N_2$, arrows for
    x(ξ)). *see:* "two identical pictures at different scales"; *read:* "compute on the right, use on the left"; *change:* "…a
    curved triangle (C13): the map becomes quadratic (10.184)".
42. `nb.primer("affine map to a parent element", "An affine map is 'stretch and shift': $x(\\xi)=\\frac h2\\xi+x_{mid}$ sends
    ξ ∈ [−1, 1] to [x_a, x_b] with h = x_b − x_a. Its derivative is the constant $dx/d\\xi=h/2$, so $dx=\\frac h2d\\xi$ (substitution,
    Ch. 3 P106) and $d\\xi/dx=2/h$ (chain rule, P49). Every element integral becomes an integral over [−1, 1].",
    code="xa, xb = 0.5, 0.75                         # an element of length h = 0.25\nh = xb - xa\nx = lambda xi: h/2*xi + (xa + xb)/2          # the map\nprint(x(-1), x(1), h/2, 2/h)                 # 0.5 0.75 0.125 8.0 (dx/dξ and dξ/dx)")` (**P234**)
43. `nb.derivation("D13", …)` — Part F D13 (7 steps), ref "10.68". (Uses P234, P49.)
44. `nb.note` — **N35 [B]** "**Shape-function slopes (10.67)–(10.68)** by the chain rule. > ⚠️ **The book prints**, on the element
    $[x_{A-1},x_A]$, $\\frac{dN_A}{dx}=\\frac{dN_A}{d\\xi}\\frac{d\\xi}{dx}=\\frac2{x_A-x_{A-1}}\\frac{dN_1}{d\\xi}=\\frac{-1}{x_A-x_{A-1}}$
    and $\\frac{dN_{A+1}}{dx}=\\frac1{x_A-x_{A-1}}$; **the correct labels are** $N_{A-1}$ (↔ $N_1$, slope $-1/h^e$) and $N_A$ (↔ $N_2$,
    slope $+1/h^e$) — the text two lines earlier says so, and $N_{A+1}$ is zero on this element (slip R1)." + `nb.code`:
    `FEM1.shape_slopes(0.5, 0.75)` and `FEM1.shape_slopes(0.5, 0.75, printed=True)`. *expect:* labels (A−1, A) with slopes (−4, 4);
    the printed labels (A, A+1) flagged "fails the chain-rule check".
45. `nb.note` — **N36 [B]** "**Global = sum of element pieces (10.69)–(10.70):** $\\mathbf M=\\sum_{e=1}^{n_{el}}\\mathbf M^e$,
    $\\mathbf K=\\sum_e\\mathbf K^e$, $\\mathbf F=\\sum_e\\mathbf F^e$ — because an integral over [0, L] is the sum of integrals over
    the elements." equation (10.69).
46. `nb.primer("scatter-add assembly (np.add.at, COO duplicates)", "Adding a small block into a big matrix at chosen rows and
    columns, where several blocks hit the same entry: `np.add.at(K, (rows, cols), block)` adds *every* contribution (plain
    `K[rows, cols] += block` would keep only the last one when indices repeat). For sparse matrices, collect (row, col, value)
    triples and build `scipy.sparse.coo_matrix(...).tocsr()` — duplicate entries are summed.", code="import numpy as np\nK = np.zeros((3, 3))\nblock = np.array([[1.0, -1.0], [-1.0, 1.0]])    # one element's stiffness\nfor e in (0, 1):                                  # two elements: nodes (0,1) and (1,2)\n    idx = np.array([e, e + 1])\n    np.add.at(K, (idx[:, None], idx[None, :]), block)\nprint(K)                                          # middle diagonal 2: both elements added there")`
    (**P235**)
47. `nb.note` — **N37 [B]** "**Element integrals (10.71)–(10.73):** $M^e_{AB}=\\int_{\\Omega^e}N_AN_B\\,dx$ (10.71),
    $K^e_{AB}=u\\int_{\\Omega^e}N_{B,x}N_A\\,dx+D\\int_{\\Omega^e}N_{B,x}N_{A,x}\\,dx$ (10.72), and the force (10.73) whose $Dq$ term lives
    only on the last element. > ⚠️ **The book prints** 'the nonzero ones require that A = e or e + 1 and B = e or e + 1'; **with**
    $\\Omega^e=[x_{A-1},x_A]$ and the map (10.78) **the nonzero entries are** $A,B\\in\\{e-1,e\\}$ (slip R2) — node 0 is the Dirichlet
    node, which is why the first element's force is $-g\\,k^1_{a1}$." equations (10.71)–(10.72).
48. `nb.note` — **N38 [B]** "**Local → global node map (10.78):** A = e − 1 for a = 1, A = e for a = 2." + `nb.code`:
    `FEM1.connectivity(4)` and `FEM1.connectivity(4, printed=True)`. *expect:* [[0, 1], [1, 2], [2, 3], [3, 4]]; the printed version
    [[1, 2], …, [4, 5]] points at a node 5 that does not exist.
49. `nb.derivation("D14", …)` — Part F D14 (11 steps), ref "10.77". (Uses P234, P106, D13, P235.)
50. `nb.worked_example("the element blocks with h = 0.25, u = 1, D = 0.25", "1. Mass: (h/6)[[2, 1], [1, 2]] = [[0.0833, 0.0417],
    [0.0417, 0.0833]]. 2. Convection: (u/2)[[−1, 1], [−1, 1]] = [[−0.5, 0.5], [−0.5, 0.5]] — not symmetric. 3. Diffusion: (D/h)[[1, −1],
    [−1, 1]] = [[1, −1], [−1, 1]]. 4. k^e = [[0.5, −0.5], [−1.5, 1.5]]; each row sums to 0 (a constant T has no convection and no
    diffusion). 5. Node 2 sits in element 2 (as local node 2) and element 3 (as local node 1): its row collects (k₂₁, k₂₂) = (−1.5,
    1.5) at columns (1, 2) and (k₁₁, k₁₂) = (0.5, −0.5) at columns (2, 3): total (−1.5, 2.0, −0.5) — the row C07 found from (10.63).")`
51. `nb.code` — `m, k = FEM1.element_matrices_linear(0.25, 1.0, 0.25)`; `mg, kg = FEM1.element_integrals(0.5, 0.75, 1.0, 0.25,
    quad=2)`; `assert np.allclose(m, mg) and np.allclose(k, kg)`; `trace = FEM1.assembly_trace(4, 1.0, 0.25)`; print, for each
    element, its global nodes and the running K row of node 2. *expect:* m, k as the tiny example; the row of node 2 is (0, 0, 0)
    → (−1.5, 1.5, 0) after e = 2 → (−1.5, 2.0, −0.5) after e = 3. *explain:* 1. closed-form block (D14); 2. the same by 2-point
    Gauss–Legendre (Ch. 5 P143) — exact for these polynomial integrands; 3. the element loop step by step.
52. `nb.check_agree` — **from scratch (curation §7):** an element loop with `np.add.at` building dense M and K for 4 elements
    (5 nodes), then deleting the Dirichlet rows/columns and solving the steady problem with `np.linalg.solve`: `assert
    np.allclose(K_mine[1:-1, 1:-1], FEM1.assemble_1d(nodes, 1.0, 0.25, dirichlet_right=1.0)[1].toarray())` and `assert
    np.allclose(T_mine, FEM1.solve_steady(nodes, 1.0, 0.25, T_L=1.0)["T"])`.
53. `nb.figure` — **"Assembly, one element at a time"** (10 × 3.2 in, four small panels = the global 5 × 5 K after e = 1, 2, 3, 4):
    heat-coloured cells, the current element's 2 × 2 block outlined amber, numbers printed in the cells. *see:* "a band filling in
    from the top-left"; *read:* "interior diagonal entries receive two pieces; the matrix is tridiagonal because hats only
    overlap neighbours"; *change:* "…triangles in 2-D: the same loop with 3 × 3 (P1) or 6 × 6 (P2) blocks and an irregular pattern
    (N102's sparsity picture)".
54. `nb.explainer("fem_hat_assembly", heading="Where do finite-element matrices come from?", why="The transport steps element
    by element, lighting up the element on the mesh, its two hat halves and the four matrix cells it adds to — the sum sign of
    (10.69) happens in front of you, and the interior row becomes the FD stencil with a ⅙–⅔–⅙ mass.", tries=["Press the 'n = 4
    by hand' preset and step through the four elements; stop at e = 3 and read node 2's row.", "Click a matrix cell to see which
    elements contributed and the integral behind each number.", "Raise R until the FE solution wiggles — it wiggles exactly
    where centred FD does (C09)."])`
55. `nb.note` — **N39 [C]** "**2-D and 3-D** follow the same steps — parent element, map, element matrices, scatter-add — with
    triangles or tetrahedra. Pointer: C13 (mixed elements) and the cylinder mesh of §10.5 (N93–N109)."
56. `nb.md` — **What would change if…** "…the flow were fast and the diffusion small? Both FD and FE produce the same centred
    equations — and the same wiggles. §10.4 starts there."

---

### A.4 §10.4 Incompressible Viscous Fluid Flow — R06 R07, C09 · R09, C10 · R08, C11 · C12 · R10, C13
1. `nb.section("10.4", "Incompressible Viscous Fluid Flow", intro="**What is this section about?** Now the real equations:
   Navier–Stokes with ∇·u = 0. Two difficulties stand between them and a working code. First, when convection dominates, thin
   layers form that a grid cannot resolve, and centred schemes answer with wiggles (C09). Second, incompressibility is not an
   evolution equation but a *constraint*: nothing tells the pressure how to change in time, yet the pressure must be whatever
   keeps every cell's net outflow zero. Four answers follow — let the fluid be slightly compressible (C10), split each time step
   and project the velocity onto the divergence-free fields (C11) on a staggered grid (C12), or choose finite-element spaces
   that respect the constraint (C13). The staggered grid and the projection are exactly what ocean and atmosphere models use.")`
2. `nb.recap("R06", "Incompressible Navier–Stokes", "Ch. 4 gave $\\rho\\Big(\\frac{\\partial\\mathbf u}{\\partial t}+(\\mathbf u\\cdot
   \\nabla)\\mathbf u\\Big)=\\rho\\mathbf g-\\nabla p+\\mu\\nabla^2\\mathbf u$ (10.79) and $\\nabla\\cdot\\mathbf u=0$ (10.80) — the
   incompressible momentum and mass balances (Ch. 4 (4.39b), (4.10)). `core.navier_stokes` residuals are the code-verification
   tool for every solver below.", where="Ch. 4 §4.2, §4.6")`
3. `nb.recap("R07", "Dimensionless Navier–Stokes", "Scaling lengths by L, speeds by U, time by L/U and pressure by ρU² (Ch. 4
   §4.11) leaves one number: $\\frac{\\partial\\mathbf u}{\\partial t}+(\\mathbf u\\cdot\\nabla)\\mathbf u=\\mathbf g-\\nabla p+
   \\frac1{Re}\\nabla^2\\mathbf u$ (10.81), Re = UL/ν. Every solver of §10.4–10.5 integrates this form.", where="Ch. 4 §4.11")`

#### C09 — Convection-dominated flow: wiggles for R_cell > 2 and numerical diffusion
4. `nb.core("C09", "Convection-dominated flow: $\\delta=O\\big(\\frac{\\Delta x}{R_{cell}}\\big)$ (10.92), wiggles for
   $R_{cell}>2$, and upwinding's extra diffusivity $u\\frac{\\partial T}{\\partial x}=D(1+0.5R_{cell})\\frac{\\partial^2T}{\\partial
   x^2}$ (10.94)", question="Why does a consistent, centred scheme produce negative temperatures next to a hot wall — and why is the
   upwind cure 'stable but wrong'?")`
5. `nb.note` — **N40 [C]** "**Scope of §10.4:** primitive variables (velocity and pressure); streamfunction–vorticity forms exist
   (the ψ Poisson solve of Ch. 6 `solve_poisson`) but do not extend to 3-D; laminar flow only — turbulence models are Ch. 12."
6. `nb.note` — **N42 [B]** "**The two difficulties** (the map of this section): (1) convection-dominated problems oscillate on grids
   that do not resolve thin layers → C09; (2) continuity is a constraint that determines the pressure → C10 (weak compressibility),
   C11–C12 (projection), C13 (mixed elements)."
7. `nb.md` — **The problem in plain words:** "Warm water flows at speed u towards a wall held at temperature 1, while the inflow
   is at 0. Heat carried towards the wall must be conducted away through a thin layer next to it; far from the wall nothing
   happens. If the grid spacing is wider than that layer, the centred scheme cannot balance 'heat in' and 'heat out' at the last
   interior node — and it answers with a temperature *below* the inflow value, then above, then below: a zigzag. The same
   happens to a sharp tracer front in an ocean model. Upwinding removes the zigzag, but only by quietly making the fluid more
   diffusive."
8. `nb.md` — **The idea** (ASCII):
   ```
   T                                  exact: e^{−R(1−x/L)}, a layer of thickness L/R at the wall
   1 ┤                         ╭─●          centred:  T_j = (r^j − 1)/(r^n − 1),  r = (1 + R_cell/2)/(1 − R_cell/2)
     │                        ╱               R_cell < 2: r > 0  ⇒ smooth     R_cell > 2: r < 0  ⇒ signs alternate
   0 ┼──●───●───●───●───●───●╯              upwind:   r = 1 + R_cell > 0 always — but it solves D(1 + 0.5 R_cell)
     └─────────── x ───────────┘ L
   ```
9. `nb.note` — **N43 [B]** "**The steady test problem (10.84)–(10.85):** $u\\frac{\\partial T}{\\partial x}=D\\frac{\\partial^2T}
   {\\partial x^2}$, $0\\le x\\le L$, with T(0) = 0 and T(L) = 1 — steady (10.1): convection towards the wall at x = L balances
   conduction." equations (10.84), (10.85).
10. `nb.md` — **🔁 Reminder:** a linear ODE with constant coefficients is solved by trying $e^{mx}$ (Ch. 1 P44); the Péclet number
    R = uL/D (P229, C05).
11. `nb.derivation("D15", …)` — Part F D15 (8 steps), ref "10.89". (Uses P44, P229, P107.)
12. `nb.note` — **N44 [B], N45 [B]** "**Exact solution and its layer (10.86)–(10.89):** $T=\\frac{e^{Rx/L}-1}{e^R-1}$ (10.86) with the
    global Péclet number $R=uL/D$ (10.87); for large R, $T=e^{-R(1-x/L)}$ (10.88), a layer of thickness $\\frac\\delta L=
    O\\big(\\frac1{\\lvert R\\rvert}\\big)$ (10.89). One layer thickness from the wall T = e⁻¹ = 0.368 of the wall value, two thicknesses
    e⁻² = 0.135 (computed, not quoted). ⚠️ The printed form (10.86) overflows for R ≳ 700 — `FD.steady_cd_exact` uses the scaled
    form $e^{-R(1-x/L)}\\frac{1-e^{-Rx/L}}{1-e^{-R}}$ (`np.expm1`, Ch. 3 P107)." equations (10.86)–(10.89) + `nb.code`: `print(np.exp(-1),
    np.exp(-2), FD.steady_cd_exact(np.array([0.99, 0.999]), 1e4), FD.cd_layer_thickness(100.0))`. *expect:* 0.3679, 0.1353,
    [e⁻¹⁰⁰ ≈ 3.7e-44, e⁻¹⁰ = 4.54e-5], 0.01.
13. `nb.note` — **N46 [B]** "**Centred differences (10.90)–(10.91):** $\\frac{u\\Delta x}{2D}(T_{j+1}-T_{j-1})=T_{j+1}-2T_j+T_{j-1}$
    (10.90), i.e. $0.5R_{cell}(T_{j+1}-T_{j-1})=T_{j+1}-2T_j+T_{j-1}$ (10.91) with $\\Delta x=L/n$ and $R_{cell}=u\\Delta x/D=R/n$.
    (Grid index j here, as in the book.)" equations (10.90), (10.91).
14. `nb.primer("linear recurrence with constant coefficients (geometric trial)", "An equation linking neighbours with constant
    weights, $aT_{j+1}+bT_j+cT_{j-1}=0$, is the discrete twin of a constant-coefficient ODE (Ch. 1 P44). Try $T_j=r^j$: dividing by
    $r^{j-1}$ gives the quadratic $ar^2+br+c=0$. With two roots $r_1\\ne r_2$ the general solution is $A r_1^j+Br_2^j$; the two end
    values fix A and B. A **negative** root makes $r^j$ flip sign at every step.", code="import numpy as np\na, b, c = 1.0, -2.0, -3.0                 # T_{j+1} - 2T_j - 3T_{j-1} = 0\nr = np.roots([a, b, c]); print(r)           # [ 3. -1.]: one root is negative\nprint([(-1.0)**j for j in range(5)])        # 1, -1, 1, -1, 1: the zigzag it produces")`
    (**P236**)
15. `nb.derivation("D16", …)` — Part F D16 (12 steps), ref "10.92". (Uses P236.)
16. `nb.md` — **The A-item result, in one box:** "Centred: $r=\\frac{1+R_{cell}/2}{1-R_{cell}/2}$, negative ⇔ $R_{cell}>2$ ⇔ the layer
    $\\delta=O\\big(\\frac L{nR_{cell}}\\big)=O\\big(\\frac{\\Delta x}{R_{cell}}\\big)$ (10.92) is thinner than half a cell. The wiggle is the
    scheme's honest warning that the grid does not resolve the layer."
17. `nb.note` — **N47 [B]** "**Cures that keep second order:** refine the grid until R_cell < 2, or refine only where the layer is
    (a stretched grid). At the same n = 10 and R = 40 (R_cell = 4) a stretched grid removes the wiggles." + `nb.code`:
    `FD.wiggle_indicator(FD.steady_cd_fd(10, 40.0, "central", grid="stretched", beta_s=2.5)["T"])`. *expect:* wiggles False (the
    uniform run: True, min T ≈ −0.3).
18. `nb.note` — **N48 [B]** "**First-order upwind (10.93):** $R_{cell}(T_j-T_{j-1})=T_{j+1}-2T_j+T_{j-1}$. > ⚠️ **The book calls
    this 'a forward-difference scheme'; $T_j-T_{j-1}$ is a *backward* difference** — upwind for u > 0 (slip R3). A truly forward
    (downwind) convective difference wiggles at every R_cell (`scheme='forward'` in the code, for the curious). Its discrete root
    is r = 1 + R_cell > 0: no sign change, ever." equation (10.93).
19. `nb.derivation("D17", …)` — Part F D17 (8 steps), ref "10.94". (Uses R02 Taylor, D01.)
20. `nb.note` — **N49 [B]** "**Upwinding = extra diffusion (10.94):** the upwind scheme is consistent with
    $u\\frac{\\partial T}{\\partial x}=D(1+0.5R_{cell})\\frac{\\partial^2T}{\\partial x^2}$ — a numerical diffusivity $0.5R_{cell}D=u\\Delta x/2$
    on top of the physical D. Accurate only if $0.5R_{cell}\\ll1$ — which is exactly the grid on which centred differences did not
    wiggle anyway. **Climate hook:** first-order upwind tracer advection in an ocean model with u = 0.1 m/s and Δx = 10 km adds
    500 m²/s of horizontal diffusion — often more than the physical eddy value; that is why models use higher-order or
    flux-limited schemes." equation (10.94) + `nb.code`: `FD.numerical_diffusivity(0.1, 1e4, scheme="upwind_steady")`. *expect:*
    500.0 m²/s.
21. `nb.note` — **N50 [C]** "**Better cures, named:** higher-order upwind schemes (Fletcher 1988), streamline-upwind/Petrov–Galerkin
    finite elements (Brooks & Hughes 1982; the Petrov–Galerkin idea of N24) and stabilised least-squares FE (Franca et al. 1992;
    the same GLS returns in C13 against spurious pressures)."
22. `nb.worked_example("four cells and a thin layer", "L = 1, n = 4 (Δx = 0.25), R = 16 so R_cell = 4. 1. Centred root r =
    (1 + 2)/(1 − 2) = −3. 2. $T_j=\\frac{(-3)^j-1}{(-3)^4-1}=\\frac{(-3)^j-1}{80}$: T₁ = −4/80 = −0.05, T₂ = 8/80 = 0.1, T₃ = −28/80 = −0.35 —
    negative temperatures, alternating. 3. Upwind root r = 1 + 4 = 5: $T_j=\\frac{5^j-1}{624}$: 0.0064, 0.0385, 0.199 — smooth but far too
    big (the exact values are 0.00001, 0.00034, 0.0183: the layer is only 1/16 thick). 4. Upwind's numerical diffusivity: 0.5 × 4 × D
    = 2D — the scheme behaves as if D were three times larger, R_eff = 16/3. 5. With n = 16 (R_cell = 1) both schemes are fine.")`
23. `nb.code` — `for n in (4, 16): print(n, FD.steady_cd_fd(n, 16.0, "central")["T"].round(4),
    FD.steady_cd_fd(n, 16.0, "upwind")["T"].round(4))`; `print(FD.discrete_root(4.0), FD.discrete_root(4.0, "upwind"))`;
    `print(FD.wiggle_indicator(FD.steady_cd_fd(4, 16.0)["T"]))`; then a scan: `[FD.wiggle_indicator(FD.steady_cd_fd(20, 20*Rc)["T"])
    ["wiggles"] for Rc in (1.9, 1.99, 2.01, 2.5, 4)]`. *expect:* n = 4 central [0, −0.05, 0.1, −0.35, 1], upwind [0, 0.0064, 0.0385,
    0.1987, 1]; roots −3.0, 5.0; wiggles True; scan [False, False, True, True, True]. *explain:* 1. tridiagonal solve of (10.91) and
    (10.93); 2. the root decides the sign pattern; 3. the switch sits exactly at R_cell = 2.
24. `nb.check_agree` — **from scratch (curation §7):** the geometric closed form $T_j=\\frac{r^j-1}{r^n-1}$ with
    $r=\\frac{1+R_{cell}/2}{1-R_{cell}/2}$ written in two lines, against the tridiagonal solve for (n, R) ∈ {(10, 5), (10, 40), (40, 40)}:
    `assert np.allclose(T_closed, FD.steady_cd_fd(n, R)["T"], atol=1e-13)`.
25. `nb.figure` — **"Wiggles are the grid confessing"** (10 × 3.6 in, three panels): (a) R_cell = 1 (n = 16, R = 16): exact (muted),
    centred (orange dots), upwind (teal dots) — all close; (b) R_cell = 4 (n = 4): centred zigzags below zero (rose shading
    under T = 0), upwind smooth but too thick, the modified-equation solution with D(1 + 0.5 R_cell) as a dashed teal line close to
    the upwind dots; (c) r vs R_cell ∈ [0, 6] for centred (orange, crossing zero at 2 and diverging there) and upwind (teal, the
    line 1 + R_cell), the band r < 0 shaded rose. *see:* "one clean panel, one zigzag, one pair of root curves"; *read:* "the zigzag
    appears exactly where the orange root goes negative"; *change:* "…a stretched grid with the same n (N47): the orange dots
    follow the layer without a single negative value".
26. `nb.plotly` — **F5** `slider_figure` over R_cell (30 values 0.5 … 8, n = 10): exact, centred, upwind profiles on [0.5, 1].
    *explain:* "the page version of E4".
27. `nb.live` — `live(lambda R, n, grid: <the (b) panel>, R=(1, 500, 1), n=(4, 100, 1), grid=["uniform", "stretched"])` —
    paired with F5.
28. `nb.explainer("cell_peclet_wiggles", heading="Why does a centred scheme give negative temperatures?", why="Sliding R or n
    moves R_cell across 2: the centred dots flip into a zigzag at exactly that point while, in the second view, the root r crosses
    zero; the modified-equation ghost lies on the upwind dots — upwind solves a different problem well.", tries=["Start at
    R_cell = 1 and drag n down until the first negative value appears; read R_cell.", "Switch on the diffusivity bars: at
    R_cell = 4 the numerical part is twice the physical one.", "Try the stretched grid at R_cell = 4: same n, no wiggles."])`
29. `nb.md` — **What would change if…** "…the unknown were a velocity field that must also stay divergence-free? The wiggle
    problem remains (it returns in Ch. 12 as a resolution question and in Ch. 13 for tracer advection); the new difficulty is the
    pressure. C10 takes the gentlest route: let the fluid compress a little."
30. `nb.recap("R09", "Compressible Navier–Stokes in two dimensions", "With Stokes' hypothesis μ_v = 0 (Ch. 4 §4.5, the ★★★
    constitutive derivation) the 2-D equations are $\\frac{\\partial\\rho}{\\partial t}+\\frac{\\partial(\\rho u)}{\\partial x}+
    \\frac{\\partial(\\rho v)}{\\partial y}=0$ (10.96) and $\\frac{\\partial}{\\partial t}(\\rho u)+\\frac{\\partial}{\\partial x}(\\rho u^2)+
    \\frac{\\partial}{\\partial y}(\\rho vu)=\\rho g_x-\\frac{\\partial p}{\\partial x}+\\mu\\nabla^2u+\\frac\\mu3\\frac{\\partial}{\\partial x}
    \\Big(\\frac{\\partial u}{\\partial x}+\\frac{\\partial v}{\\partial y}\\Big)$ (10.97), with (10.98) the same for v. `ch10.compressible_ns_sympy()`
    checks them against Ch. 4's `navier_stokes_sym(mu_v=0)`.", where="Ch. 4 §4.5–§4.6")` + `nb.code`: print the difference (0).

#### C10 — Weakly compressible Navier–Stokes and the MacCormack predictor–corrector
31. `nb.core("C10", "MacCormack's predictor–corrector for $\\mathbf U_t+\\mathbf E(\\mathbf U)_x+\\mathbf F(\\mathbf U)_y=0$
    (10.100)–(10.102)", question="If pressure has no equation of its own in incompressible flow, can we give it one by letting the
    fluid be very slightly compressible — and how do we step such a system accurately?")`
32. `nb.md` — **The problem in plain words:** "Real water *is* compressible — sound travels through it at 1500 m/s. A pressure pulse
    tells the whole tank about a moving lid by sound waves. If we keep a small, finite sound speed, the pressure gets an evolution
    equation and an explicit code can march everything forward — at the price of time steps short enough to follow the sound
    (C05's CFL with the sound speed) and an O(Ma²) error in the density. MacCormack's two-stage scheme is a classic way to do the
    marching with second-order accuracy."
33. `nb.md` — **The idea** (ASCII): "predict with forward differences, correct with backward ones from the predicted state, average":
    ```
    U*      = Uⁿ − (Δt/Δx)(E_{i+1} − E_i)ⁿ               predictor: forward difference (a guess at t_{n+1})
    Uⁿ⁺¹    = ½[Uⁿ + U* − (Δt/Δx)(E*_i − E*_{i−1})]       corrector: backward difference of the guess, then average
    forward + backward, averaged  ⇒  centred and second order (D18)
    ```
34. `nb.note` — **N51 [B]** "**Incompressibility, three ways of looking** (§10.4): continuity is a *constraint* on u that determines
    p; equivalently the pressure is a *Lagrange multiplier* that enforces it; physically, an incompressible fluid has infinite sound
    speed, so pressure information arrives everywhere instantly — which is why a Poisson equation for p (continuous or discrete)
    appears, and why solving it is often the most expensive step. Three answers follow: weak compressibility (here), projection
    (C11–C12), mixed elements (C13)."
35. `nb.primer("Lagrange multiplier", "To minimise or balance something *subject to a constraint*, add a new unknown λ times the
    constraint; λ adjusts itself until the constraint holds, and its value measures how hard the constraint pushes back. In
    incompressible flow the pressure is that λ: it takes whatever values make ∇·u = 0, and −∇p is the force needed.",
    code="import sympy as sp\nx, y, lam = sp.symbols('x y lam')\nL = x**2 + y**2 + lam*(x + y - 1)              # minimise x²+y² subject to x + y = 1\nprint(sp.solve([sp.diff(L, v) for v in (x, y, lam)], [x, y, lam]))   # x = y = 1/2, lam = -1")` (**P239**)
36. `nb.note` — **N52 [B]** "**Artificial compressibility (Chorin 1967) (10.95):** replace continuity by
    $\\frac{\\partial p}{\\partial t}+c^2\\nabla\\cdot\\mathbf u=0$ with an arbitrary c — continuity of a fluid with $p=c^2\\rho$,
    linearised, in words: a pressure that grows wherever fluid piles up. It has no physical meaning before steady state (a
    *pseudo-time*); at steady state ∇·u = 0 whatever c is. **Climate hook:** some atmospheric models slow sound down on purpose
    ('reduced speed of sound') for the same reason." equation (10.95) + `nb.code`: `r = ch10.artificial_compressibility_channel(ny=16,
    Re=10.0)`; print `r["max_err"]` and the last divergence. *expect:* plane Poiseuille reproduced to < 1e-10 (the 3-point
    Laplacian is exact for a parabola); divergence → 0.
37. `nb.note` — **N53 [B]** "**Isothermal equation of state (10.99):** $p=c^2\\rho$; at low Mach number Ma = U/c and nearly constant
    temperature, the density varies only by O(Ma²) — Ma = 0.1 gives ≈ 1 % (Ch. 4's incompressibility criterion Ma < 0.3,
    `ch04.is_incompressible_regime`)." equation (10.99) + `nb.code`: `print(0.1**2, ch04.is_incompressible_regime(34.0))`
    (U = 34 m/s in air, Ma ≈ 0.1). *expect:* 0.01, True.
38. `nb.primer("conservation (flux) form U_t + E_x + F_y = 0", "Stack the conserved quantities in one vector U = (ρ, ρu, ρv) and
    write every conservation law as 'rate of change + divergence of a flux = 0': $\\mathbf U_t+\\mathbf E(\\mathbf U)_x+\\mathbf
    F(\\mathbf U)_y=0$ with $\\mathbf E=(\\rho u,\\ \\rho u^2+p,\\ \\rho uv)$ and $\\mathbf F=(\\rho v,\\ \\rho uv,\\ \\rho v^2+p)$ (viscous
    terms added separately). One scheme then updates all components at once, and summing a flux difference over the cells
    telescopes — what flows out of one cell flows into the next.", code="import numpy as np\nrho, u, v, c = 1.0, 0.2, 0.0, 10.0          # a state; c = 1/Ma (Ma = 0.1)\np = c**2*rho                                 # (10.99)\nE = np.array([rho*u, rho*u**2 + p, rho*u*v])   # x-flux of (ρ, ρu, ρv)\nprint(E)                                     # [0.2, 100.04, 0.0]: the pressure dominates the momentum flux")`
    (**P237**)
39. `nb.primer("predictor–corrector (Heun's second-order idea)", "Estimate the next state with a cheap first-order step (predict),
    evaluate the slope again at the prediction, and average the two slopes (correct). For y′ = f(y): $y^*=y^n+\\Delta tf(y^n)$,
    $y^{n+1}=y^n+\\frac{\\Delta t}2[f(y^n)+f(y^*)]$ — second order, like the RK2 methods of Ch. 3 (P95 did RK4).",
    code="import numpy as np\nf = lambda y: -y; dt = 0.1                 # y' = -y, exact e^{-t}\ny = 1.0; ys = y + dt*f(y)                   # predict: 0.9\ny1 = y + dt/2*(f(y) + f(ys))                # correct: 0.905\nprint(y1, np.exp(-0.1))                     # 0.905 vs 0.904837: error 1.6e-4 ~ dt^3 per step")` (**P238**)
40. `nb.derivation("D18", …)` — Part F D18 (14 steps, ★★★, with `check_src`), ref "10.102". (Uses P238, P98, P117, D04 moves, P225.)
41. `nb.worked_example("one MacCormack step on five points", "Linear advection T_t + uT_x = 0 (E = uT), C = uΔt/Δx = 0.5, periodic,
    start [0, 0, 1, 0, 0]. 1. Predictor $T^*_i=T_i-C(T_{i+1}-T_i)$: [0, −0.5, 1.5, 0, 0]. 2. Corrector $T^{n+1}_i=\\frac12[T_i+T^*_i-
    C(T^*_i-T^*_{i-1})]$: node 1 ½[0 − 0.5 − 0.5(−0.5 − 0)] = −0.125; node 2 ½[1 + 1.5 − 0.5(1.5 + 0.5)] = 0.75; node 3 ½[0 + 0 −
    0.5(0 − 1.5)] = 0.375. 3. Result [0, −0.125, 0.75, 0.375, 0], sum 1 (conserved). 4. The Lax–Wendroff formula of D18 gives the same
    numbers: node 1: 0 − 0.25(1 − 0) + 0.125(1) = −0.125 ✓. 5. Compare upwind at the same C: [0, 0, 0.5, 0.5, 0] — no undershoot but
    a flatter peak. MacCormack keeps the peak higher but dips below zero upstream: dispersion, not diffusion.")`
42. `nb.code` — `print(MCK.maccormack_advection_1d(np.array([0, 0, 1, 0, 0.0]), 0.5, 1))`; `lin = ch10.maccormack_linear_sympy()`;
    print `lin["lw_difference"], lin["local_error"], lin["G2_factorised"]`; `for sch in ("upwind", "maccormack"): print(sch,
    [FD.advect_periodic("gauss", C=0.5, n_cells=N, n_rev=1.0, scheme=sch)["rms_error"] for N in (50, 100, 200, 400)])` and the
    observed orders. *expect:* [0, −0.125, 0.75, 0.375, 0]; difference 0; local error $-\\frac{u\\Delta t\\Delta x^2}6(1-C^2)T_{xxx}$;
    ∣G∣² = 1 − 4C²(1 − C²)sin⁴(θ/2); orders ≈ 1.0 (upwind) and ≈ 2.0 (MacCormack). *explain:* 1. one step; 2. sympy repeats D18;
    3. the measured orders.
43. `nb.check_agree` — **from scratch (curation §7):** MacCormack for linear advection in two lines (`Ts = T - C*(np.roll(T, -1) -
    T)`; `Tn = 0.5*(T + Ts - C*(Ts - np.roll(Ts, 1)))`) vs the Lax–Wendroff formula (`T - C/2*(np.roll(T, -1) - np.roll(T, 1)) +
    C**2/2*(np.roll(T, -1) - 2*T + np.roll(T, 1))`) vs `MCK.maccormack_advection_1d(T, C, 1)` on a random profile (seeded): `assert
    np.allclose(Tn, T_lw, atol=1e-15)` and `assert np.allclose(Tn, lib)`.
44. `nb.note` — **N54 [B], N55 [B], N56 [B]** "**MacCormack for the 2-D Navier–Stokes equations (10.103)–(10.109).** The same
    predictor (forward) and corrector (backward) applied to (10.96)–(10.98); the pressure is folded into the flux as $\\rho u^2+c^2\\rho$;
    the viscous second derivatives are always centred (to keep second order), with the normal one weighted 4/3 and a cross term.
    A colour-coded table (flux differences purple · viscous normal $\\frac43c_3$ blue · viscous tangential $c_4$ blue · cross term
    $c_5$ amber) shows the x-momentum predictor (10.104) $(\\rho u)^*_{i,j}=(\\rho u)^n_{i,j}-c_1[(\\rho u^2+c^2\\rho)^n_{i+1,j}-
    (\\rho u^2+c^2\\rho)^n_{i,j}]-c_2[(\\rho uv)^n_{i,j+1}-(\\rho uv)^n_{i,j}]+\\frac43c_3(u^n_{i+1,j}-2u^n_{i,j}+u^n_{i-1,j})+
    c_4(u^n_{i,j+1}-2u^n_{i,j}+u^n_{i,j-1})+c_5(v^n_{i+1,j+1}+v^n_{i-1,j-1}-v^n_{i+1,j-1}-v^n_{i-1,j+1})$, the corrector (10.107)
    in the same colours, and the coefficients $c_1=\\frac{\\Delta t}{\\Delta x}$, $c_2=\\frac{\\Delta t}{\\Delta y}$, $c_3=
    \\frac{\\mu\\Delta t}{\\Delta x^2}$, $c_4=\\frac{\\mu\\Delta t}{\\Delta y^2}$, $c_5=\\frac{\\mu\\Delta t}{12\\Delta x\\Delta y}$
    (10.109). Why 1/12: the cross derivative's centred stencil carries $\\frac1{4\\Delta x\\Delta y}$ and the Stokes term carries $\\frac\\mu3$."
    equations (10.104), (10.109) + `nb.code`: `MCK.ns_coefficients(1e-3, 1/64, 1/64, 0.01)`.
45. `nb.note` — **N57 [B]** "**Arrangements FF/BB, BB/FF, FB/BF, BF/FB:** which direction is forward in the predictor (x then y);
    cycling them from step to step removes a small bias. On the advection test FB and BF differ by O(Δx²)." + `nb.code`:
    `FD.advect_periodic("gauss", C=0.8, n_cells=N, scheme="maccormack")` vs `"maccormack_bf"` for N = 50, 100, 200; print the
    difference norms. *expect:* the difference falls by ≈ 4 per halving.
46. `nb.note` — **N58 [B]** "**MacCormack's time-step limit (10.110)** (semi-empirical, Tannehill, Anderson & Pletcher 1997, as cited):
    $\\Delta t\\le\\frac{\\sigma}{1+2/Re_\\Delta}\\Big[\\frac{\\lvert u\\rvert}{\\Delta x}+\\frac{\\lvert v\\rvert}{\\Delta y}+c\\sqrt{\\frac1{\\Delta x^2}+
    \\frac1{\\Delta y^2}}\\Big]^{-1}$ with a safety factor σ (we pass our σ = 0.8) and the minimum mesh Reynolds
    number $Re_\\Delta$. The sound speed dominates the bracket at low Ma: a CFL condition built on c (C05)." equation (10.110) +
    `nb.code`: `MCK.maccormack_dt(1.0, 0.0, 12.5, 1/64, 1/64, 1.0, 0.01, sigma=0.8)` (lid speed, c = 1/Ma = 12.5 for our Ma = 0.08, μ = 1/Re =
    0.01). *expect:* ≈ 2.9e-4 (the bracket ≈ 1195, 1 + 2/Re_Δ ≈ 2.28); the asymptotic (10.155) gives 7.1e-4 (N91).
47. `nb.note` — **N59 [C]** "**Density (pressure) boundary conditions are the key difficulty** of explicit MacCormack for
    weakly compressible flow — worked out for the cavity (C14, N79–N84) and the block (C15, N87–N89)."
48. `nb.figure` — **"Diffusion versus dispersion"** (10 × 3.4 in, three panels): (a) square pulse after one revolution, C = 0.8,
    N = 100: exact (muted), upwind (teal, rounded and lower), MacCormack (purple, sharp with ripples trailing — waves lag);
    (b) ∣G(θ)∣ vs θ for upwind and Lax–Wendroff at C = 0.8 (`FD.amplification_modulus`); (c) relative phase speed vs θ
    (`FD.phase_error`): upwind ≈ 1, LW < 1 for short waves (0.91 at θ = π/2). *see:* "one curve too smooth, one too wiggly";
    *read:* "upwind damps short waves (∣G∣ < 1) but moves them at the right speed; MacCormack keeps their amplitude but lets them
    fall behind — the ripples are short waves arriving late"; *change:* "…C = 1: both curves in (b) and (c) are exactly 1 — the
    exact shift".
49. `nb.animation` — **A2** (video, 60 frames; FAST 30): the square pulse travelling one revolution at C = 0.8 with exact (muted),
    upwind (teal) and MacCormack (purple) on one clock (`FD.advect_periodic` snapshots). *explain:* "watch the upwind pulse sag
    while the MacCormack one grows a trailing ripple".
50. `nb.md` — **What would change if…** "…we refuse to let the fluid compress at all? Then the time step is no longer tied to a
    (fake) sound speed, but the pressure must be solved for every step: C11 splits the step so that one part is exactly that
    solve."
51. `nb.recap("R08", "The conservative convective term", "Because ∇·u = 0, $(\\mathbf u\\cdot\\nabla)\\mathbf u=\\nabla\\cdot(\\mathbf
    u\\mathbf u)$ (10.82) — the product rule $\\nabla\\cdot(\\mathbf u\\mathbf u)=(\\mathbf u\\cdot\\nabla)\\mathbf u+\\mathbf u(\\nabla\\cdot
    \\mathbf u)$ with the last term zero (Ch. 4, `conservative_to_advective_sym`). `MAC.convective_terms(form=…)` offers both.",
    where="Ch. 4 §4.4")`

#### C11 — Operator splitting and the projection method
52. `nb.core("C11", "Operator splitting and projection: $\\frac{\\mathbf u^{n+1}-\\mathbf u^{n+1/2}}{\\Delta t}+\\nabla p^{n+1}=
    \\mathbf 0$ (10.117), $\\nabla\\cdot\\mathbf u^{n+1}=0$ (10.118)", question="How can one time step move the flow forward *and*
    keep it exactly divergence-free?")`
53. `nb.md` — **The problem in plain words:** "Step the momentum equation forward ignoring the pressure: the fluid is pushed around
    by its own inertia and viscosity, and in some cells more fluid now arrives than leaves — a divergence that an incompressible
    fluid cannot have. Now ask: what is the smallest correction that removes it? Push with a pressure gradient: fluid is pushed out
    of the cells that were over-full. The pressure that does exactly this solves a Poisson equation, and the push it gives is a
    pure gradient — it removes the divergent part of the velocity and touches nothing else."
54. `nb.md` — **The idea** (ASCII):
    ```
    uⁿ ──[ A₁: convect + diffuse, explicit ]──► u^{n+½}  (∇·u^{n+½} ≠ 0)
                                                     │  solve  ∇²p^{n+1} = ∇·u^{n+½}/Δt        (Poisson)
                                                     ▼
    u^{n+1} = u^{n+½} − Δt ∇p^{n+1}   ──►   ∇·u^{n+1} = 0     (A₂: pressure/continuity, implicit)
    any field = divergence-free part + gradient part — projection keeps the first and deletes the second
    ```
55. `nb.primer("operator splitting and the commutator [A₁, A₂]", "To solve dφ/dt + (A₁ + A₂)φ = 0, take the step in two parts:
    first with A₁ alone, then with A₂ alone. Exact would be $e^{-\\Delta t(A_1+A_2)}$; splitting gives $e^{-\\Delta tA_2}e^{-\\Delta
    tA_1}$. For numbers these agree; for matrices they differ by $\\frac{\\Delta t^2}2[A_1,A_2]+O(\\Delta t^3)$, where the commutator
    $[A_1,A_2]=A_1A_2-A_2A_1$ measures how much the order matters. An O(Δt²) error per step is O(Δt) overall: first order.",
    code="import numpy as np\nfrom scipy.linalg import expm                      # matrix exponential (Ch. 2 P79)\nA1 = np.array([[1., 1.], [0., 1.]]); A2 = np.array([[1., 0.], [-1., 2.]])\nprint(A1 @ A2 - A2 @ A1)                           # [[-1, 1], [0, 1]]: they do not commute\nfor dt in (0.1, 0.05):\n    print(dt, np.abs(expm(-dt*A2) @ expm(-dt*A1) - expm(-dt*(A1 + A2))).max())   # error ÷ ≈ 4 (4.1e-3 → 1.1e-3) when dt halves: O(dt²) per step")`
    (**P241**)
56. `nb.note` — **N60 [B]** "**Operator splitting (10.111)–(10.112):** $\\frac{d\\phi}{dt}+A(\\phi)=f$, $\\phi(0)=\\phi_0$, with
    $A(\\phi)=A_1(\\phi)+A_2(\\phi)$ — each substep handles one part with the method best for it. **Climate hook:** every GCM splits
    'dynamics' (advection, pressure, Coriolis) from 'physics' (radiation, clouds, turbulence) in exactly this way." equations
    (10.111), (10.112) + `nb.code`: `sys2 = ch10.split_linear_system()`; print its commutator.
57. `nb.note` — **N61 [B]** "**Marchuk–Yanenko fractional steps (10.113)–(10.114):** $\\frac{\\phi^{n+1/2}-\\phi^n}{\\Delta t}+
    A_1(\\phi^{n+1/2})=f_1^{n+1}$ then $\\frac{\\phi^{n+1}-\\phi^{n+1/2}}{\\Delta t}+A_2(\\phi^{n+1})=f_2^{n+1}$ (both implicit, with
    $f_1+f_2=f$) — first order in time (each substep is backward Euler, and the split adds $\\frac{\\Delta t^2}2[A_1,A_2]$ per step)."
    equations (10.113), (10.114) + `nb.code`: `ch10.splitting_order("marchuk_yanenko")`. *expect:* 1.0 ± 0.1.
58. `nb.note` — **N62 [B]** "**The MAC split (10.115):** $\\mathbf A_1(\\mathbf u,p)=\\begin{pmatrix}(\\mathbf u\\cdot\\nabla)\\mathbf u
    -\\frac1{Re}\\nabla^2\\mathbf u\\\\\\mathbf 0\\end{pmatrix}$ (teal: convection–diffusion) and $\\mathbf A_2(\\mathbf u,p)=
    \\begin{pmatrix}\\nabla p\\\\\\nabla\\cdot\\mathbf u\\end{pmatrix}$ (orange: pressure and continuity)." equation (10.115).
59. `nb.note` — **N63 [B]** "**Substep 1, explicit (10.116):** $\\frac{\\mathbf u^{n+1/2}-\\mathbf u^n}{\\Delta t}+(\\mathbf u^n\\cdot\\nabla)
    \\mathbf u^n-\\frac1{Re}\\nabla^2\\mathbf u^n=\\mathbf g^{n+1}$ — an FTCS-like step (C02) in 2-D; `MAC.predictor`." equation (10.116).
60. `nb.primer("Helmholtz–Hodge decomposition", "Any smooth vector field can be split into a divergence-free part and a gradient:
    $\\mathbf w=\\mathbf u+\\nabla\\phi$ with ∇·u = 0. Taking the divergence gives $\\nabla^2\\phi=\\nabla\\cdot\\mathbf w$ (a Poisson problem,
    Ch. 5 P139), so φ — and then u = w − ∇φ — is found by one solve. The gradient part is curl-free (curl grad = 0, Ch. 2); the
    boundary condition on the normal velocity makes the split unique.", code="import sympy as sp\nx, y = sp.symbols('x y')\nw = sp.Matrix([x + y, 0])                  # a field with divergence 1\nphi = x**2/2                                 # solves phi_xx + phi_yy = div w = 1\nu = w - sp.Matrix([sp.diff(phi, x), sp.diff(phi, y)])\nprint(u.T, sp.diff(u[0], x) + sp.diff(u[1], y))   # (y, 0) and divergence 0")` (**P240**)
61. `nb.derivation("D19", …)` — Part F D19 (9 steps), ref "10.118". (Uses P240, P139, curl of a gradient (Ch. 2).)
62. `nb.note` — **N72 [B]** "**The projection method** (Chorin 1968, Temam 1969): applied explicitly on the staggered grid it is the
    MAC scheme (apart from boundary details). Its meaning (D19's result): the first step leaves a field with divergence; the
    second adds an **irrotational** correction — the gradient of a potential proportional to p — that removes exactly that
    divergence. The discrete curl of the correction is zero to round-off (checked in C12)."
63. `nb.note` — **N41 [B]** "**Initial and boundary conditions for Navier–Stokes (10.83):** $\\mathbf u(\\mathbf x,t=0)=\\mathbf u_0(\\mathbf
    x)$ with $\\nabla\\cdot\\mathbf u_0=0$; a table | boundary | velocity | pressure |: solid wall | no slip (u = wall velocity) | **none**;
    inflow | u given | (p may be given); outflow | zero tangential velocity and zero normal stress, or zero derivatives —
    artificial, so test that the answer does not depend on where it is put (N111); and p fixed at one point, because only ∇p
    appears: p + const gives the same flow." equation (10.83) + `nb.code`: one MAC step with p and with p + 7 → `assert
    np.allclose(u1, u2)` (`MAC.pin_pressure`).
64. `nb.worked_example("a projection by hand", "Δt = 1, and a predicted field $\\mathbf u^{n+1/2}=(x+y,\\ 0)$ (unbounded plane, no
    walls). 1. Its divergence: ∂(x + y)/∂x + 0 = 1 — not allowed. 2. Poisson: ∇²p = ∇·u^{n+1/2}/Δt = 1; one solution p = x²/2. 3.
    Correction −Δt∇p = −(x, 0). 4. $\\mathbf u^{n+1}=(x+y-x,\\ 0)=(y,\\ 0)$: a simple shear flow, divergence 0 ✓. 5. The removed part
    (x, 0) = ∇(x²/2) has zero curl ✓ — it is pure 'squeeze', no rotation; the vorticity of the flow (−1) is untouched. 6. With walls,
    the boundary condition would pick which solution of ∇²p = 1 is the right one (C12, D20).")`
65. `nb.code` — `ps = ch10.projection_sympy()`; print `ps["poisson"]`, `ps["curl"]`; then on a 16 × 16 MAC grid:
    `g = MAC.MacGrid(16, 16)`; the deterministic divergent field (`MAC.mac_projection_summary` builds it; here `parts =
    MAC.project(us, vs, g, dt=1.0, return_parts=True)`); print `np.abs(parts["div_before"]).max(), np.abs(parts["div_after"]).max()`.
    *expect:* ∇²p = ∇·u*/Δt; curl 0; div before O(1), after < 1e-12. *explain:* 1. sympy repeats D19; 2. one discrete projection.
66. `nb.figure` — **"Projection: keep the swirl, delete the squeeze"** (10 × 3.4 in, three quiver panels on [−1, 1]²): (a) $\\mathbf
    u^{n+1/2}=(x+y,0)$ (blue arrows) over its divergence (rose, uniform 1); (b) the correction −Δt∇p = (−x, 0) (orange arrows);
    (c) $\\mathbf u^{n+1}=(y,0)$ (blue) with divergence 0 (white). *see:* "arrows that spread out, arrows that squeeze back,
    arrows that only shear"; *read:* "the correction has no curl: it cannot change the vorticity, only the divergence"; *change:*
    "…the predicted field had no divergence: p = const, the projection does nothing".
67. `nb.animation` — **A3** (frames player, 4 stages × 8 frames; FAST 4 × 4): one MAC step in slow motion on a 16² cavity grid
    started from rest with the lid moving: uⁿ (quiver) → u^{n+½} with its divergence heatmap (rose) → p (orange contours) →
    u^{n+1} with divergence < 1e-12 (`MAC.projection_stages`). *explain:* "the rose blotches appear after the predictor and vanish
    after the projection".
68. `nb.note` — **N74 [B]** "**Glowinski's Θ-scheme (10.129)–(10.133)** — a better split: three substeps (Stokes to n + Θ, a nonlinear
    convection–diffusion step to n + 1 − Θ, Stokes to n + 1) with weights $\\alpha_\\Theta+\\beta_\\Theta=1$. It is second order only for
    $\\Theta=1-\\frac1{\\sqrt2}$ with $\\beta_\\Theta=\\frac\\Theta{1-\\Theta}$ (Glowinski; checked here by a sympy series: the z² term of
    R(z) − e^z vanishes there) — ⚠️ this Θ is not the Fourier angle θ of C04. The first Stokes substep, as printed:
    $\\frac{\\mathbf u^{n+\\Theta}-\\mathbf u^n}{\\Theta\\Delta t}-\\frac{\\alpha_\\Theta}{Re}\\nabla^2\\mathbf u^{n+\\Theta}+\\nabla p^{n+\\Theta}=
    \\mathbf g^{n+\\Theta}+\\frac{\\beta_\\Theta}{Re}\\nabla^2\\mathbf u^n-(\\mathbf u^n\\cdot\\nabla)\\mathbf u^n$, $\\nabla\\cdot\\mathbf
    u^{n+\\Theta}=0$ (10.129)–(10.130)." equation (10.129) + `nb.code` (a `check`-style sympy cell, every line commented):
    `th = ch10.theta_scheme_amplification_sympy()`; print `th["z2_coeff"]`, `th["theta_root"]`; then `ch10.splitting_order("theta",
    theta_split=1 - 1/np.sqrt(2))`, `ch10.splitting_order("theta", theta_split=0.25)`. *expect:* θ_root = 1 − 1/√2 printed by the cell to 4 significant figures (computed,
    not typed); orders ≈ 2.0 and ≈ 1.0.
69. `nb.md` — **What would change if…** "…we wrote the projection on a grid? Where the numbers live — p at cell centres, velocities
    on faces — decides whether the discrete projection is exact and whether a zigzag pressure can hide. C12."

#### C12 — The staggered (MAC) grid: discrete continuity, the pressure Poisson equation and the checkerboard
70. `nb.core("C12", "The staggered (MAC) grid: $\\nabla^2_dp^{n+1}_{i,j}=\\frac1{\\Delta t}\\Big(\\frac{u^{n+1/2}_{i+1/2,j}-u^{n+1/2}_{i-1/2,j}}
    {\\Delta x}+\\frac{v^{n+1/2}_{i,j+1/2}-v^{n+1/2}_{i,j-1/2}}{\\Delta y}\\Big)$ (10.124)", question="Why store pressure at cell
    centres and velocities on cell faces — and what goes wrong if all three sit at the same points?")`
71. `nb.md` — **The problem in plain words:** "An ocean model divides the sea into boxes. What matters for mass is the water
    crossing each face of each box; what drives that water is the pressure difference between the two boxes the face separates.
    So put the pressure in the middle of each box and the velocity on each face — then 'net outflow of a box' and 'push across a
    face' each use two neighbouring numbers, the most compact and accurate choice. Put everything at the same points instead and
    the pressure difference across a point must skip a neighbour: a pressure that zigzags (+, −, +, −) looks perfectly flat to the
    momentum equation and grows unchecked. The staggered layout is Arakawa's **C-grid**, used by most ocean and many atmosphere
    models."
72. `nb.md` — **The idea** (ASCII of one cell, Fig. 10.4 in text):
    ```
                 v_{i,j+½}  ↑
            ┌───────────────┼───────────────┐
            │                               │
     u_{i−½,j} →          p_{i,j}          → u_{i+½,j}      net outflow = (u_{i+½} − u_{i−½})/Δx + (v_{j+½} − v_{j−½})/Δy
            │                               │                push on a face = (p_{i+1} − p_i)/Δx   (adjacent values)
            └───────────────┼───────────────┘
                 v_{i,j−½}  ↑
    ```
73. `nb.primer("half-index notation and staggered array shapes", "On a staggered grid a velocity sits halfway between two pressure
    points: $u_{i+1/2,j}$ lives on the vertical face between cells (i, j) and (i + 1, j). In code there are no half indices, so each
    field gets its own array: `p[ny, nx]` (centres), `u[ny, nx+1]` (vertical faces, including both side walls), `v[ny+1, nx]`
    (horizontal faces) — the project's `[j, i]` layout (Ch. 2 P76). `u[j, i]` means $u_{i-1/2,j}$.", code="import numpy as np\nnx, ny = 3, 2\np = np.zeros((ny, nx)); u = np.zeros((ny, nx + 1)); v = np.zeros((ny + 1, nx))\nprint(p.shape, u.shape, v.shape)          # (2, 3) (2, 4) (3, 2): one more face than cells in each direction\ndiv = (u[:, 1:] - u[:, :-1]) + (v[1:, :] - v[:-1, :])   # net outflow per cell (dx = dy = 1): shape (2, 3)\nprint(div.shape)")`
    (**P242**)
74. `nb.figure` — **Fig. 10.4 remade** (`staggered_grid(ax, n=3, highlight=(2, 2))`: pressure dots at centres, u arrows on vertical
    faces, v arrows on horizontal faces, the control cell around $p_{2,2}$ shaded, the boundary Γ with the wall-normal faces drawn
    thick). *see:* "three kinds of points"; *read:* "each face velocity sits between exactly the two pressures that push it"; *change:*
    "…a collocated grid: every symbol at the same dot — see the checkerboard below".
75. `nb.note` — **N64 [B]** "**Staggered predictor (10.119)–(10.120):** $u^{n+1/2}_{i+1/2,j}=u^n_{i+1/2,j}-\\Delta t\\big(uu_x+vu_y-\\frac1{Re}
    \\nabla^2u\\big)^n_{i+1/2,j}+\\Delta t\\,g_x{}^{n+1}_{i+1/2,j}$ (10.119) and the same for v on its faces (10.120) (the book writes the body
    force as **g** = (f, g); we write $(g_x,g_y)$). Quantities not stored on a face are averaged to it (two- and four-point averages —
    our documented choice in `MAC.predictor`)." equations (10.119), (10.120).
76. `nb.note` — **N65 [B]** "**Velocity correction (10.121)–(10.122):** $u^{n+1}_{i+1/2,j}=u^{n+1/2}_{i+1/2,j}-\\frac{\\Delta t}{\\Delta x}
    (p^{n+1}_{i+1,j}-p^{n+1}_{i,j})$ and $v^{n+1}_{i,j+1/2}=v^{n+1/2}_{i,j+1/2}-\\frac{\\Delta t}{\\Delta y}(p^{n+1}_{i,j+1}-p^{n+1}_{i,j})$ —
    (10.117) with the face gradient (10.126)." equations (10.121), (10.122).
77. `nb.note` — **N66 [B]** "**Discrete continuity per cell (10.123):** $\\frac{u^{n+1}_{i+1/2,j}-u^{n+1}_{i-1/2,j}}{\\Delta x}+
    \\frac{v^{n+1}_{i,j+1/2}-v^{n+1}_{i,j-1/2}}{\\Delta y}=0$ — the flux balance of one box (Gauss' theorem, Ch. 2, applied to a single
    cell): in minus out equals zero. Second order, O(Δx², Δy²), with only four values." equation (10.123).
78. `nb.primer("singular linear systems and the compatibility condition", "The pure-Neumann Poisson matrix (every boundary a wall)
    sends a constant vector to zero — adding a constant to p changes nothing. So A is singular: Ax = b has a solution only if b is
    orthogonal to that null vector, i.e. **Σ b = 0** (net inflow through all walls is zero), and then infinitely many (x + const).
    Fix it by pinning one value or asking for zero mean.", code="import numpy as np\nA = np.array([[1., -1, 0], [-1, 2, -1], [0, -1, 1]])   # 1-D Neumann Laplacian (3 cells)\nprint(np.linalg.matrix_rank(A), A @ np.ones(3))          # rank 2; A·(1,1,1) = 0\nb = np.array([1., 0, -1])                                   # sums to 0: compatible\nx = np.linalg.lstsq(A, b, rcond=None)[0]; print(x, A @ x)  # one solution (mean zero), reproduces b")`
    (**P243**)
79. `nb.primer("scipy.sparse.diags, kron and a cached splu factorisation", "A 1-D second-difference matrix is three diagonals:
    `sparse.diags([1, -2, 1], [-1, 0, 1], shape=(n, n))`. The 2-D 5-point Laplacian is built from two 1-D ones by Kronecker products:
    `sparse.kron(I_y, Lx) + sparse.kron(Ly, I_x)` (Ch. 6 P161 introduced sparse storage). The matrix never changes during a run, so
    factorise it once, `lu = splu(A.tocsc())`, and each time step only calls `lu.solve(b)` — the reason MAC codes are fast.",
    code="import numpy as np, scipy.sparse as sp\nfrom scipy.sparse.linalg import splu\nn = 4; L1 = sp.diags([1, -2, 1], [-1, 0, 1], shape=(n, n))\nA = sp.kron(sp.identity(n), L1) + sp.kron(L1, sp.identity(n))   # 16 × 16 Dirichlet Laplacian\nlu = splu(A.tocsc())                          # factorise once ...\nprint(A.shape, lu.solve(np.ones(16))[:4])    # ... solve many times")` (**P244**)
80. `nb.derivation("D20", …)` — Part F D20 (12 steps), ref "10.124". (Uses P242, P243, N65, N66.)
81. `nb.note` — **N67 [B], N70 [B]** "**The discrete pressure Poisson equation (10.124)** is the 5-point Laplacian of Ch. 6 acting on
    p, equal to the discrete divergence of the predicted velocity over Δt — and **it needs no pressure boundary condition**: at a
    wall the normal velocity is known and is not corrected, so the pressure outside the wall never appears. > ⚠️ **The book says**
    '$p_{0,2}$ will not appear in equation (10.120)'; **it means (10.124)** — (10.120) is the v-predictor (slip R4). Equivalent
    view: a zero normal pressure gradient at the wall; the result does not depend on the boundary values of $u^{n+1/2}$ (Peyret &
    Taylor), which a test checks by perturbing them." equation (10.124).
82. `nb.worked_example("three cells in a closed pipe", "One row of three cells, Δx = 1, Δt = 1, walls at both ends ($u_{1/2}=u_{7/2}=0$,
    fixed). Predicted interior faces $u_{3/2}=u_{5/2}=1$. 1. Divergences: cell 1 (1 − 0) = 1, cell 2 (1 − 1) = 0, cell 3 (0 − 1) = −1;
    sum 0 ✓ (compatible). 2. Poisson rows: cell 1: $p_2-p_1=1$ (the wall face is not corrected, so there is no $p_0$); cell 2:
    $p_1-2p_2+p_3=0$; cell 3: $p_2-p_3=-1$. 3. Pin $p_1=0$: $p_2=1$, $p_3=2$. 4. Correct: $u_{3/2}=1-(p_2-p_1)=0$, $u_{5/2}=1-(p_3-p_2)
    =0$. 5. All divergences 0: in a closed pipe of incompressible fluid nothing can move, and the pressure rises in the direction
    the fluid was trying to go.")`
83. `nb.code` — `g = MAC.MacGrid(8, 8)`; `A = MAC.pressure_poisson_matrix(g)`; print its shape, rank deficiency before pinning (via
    `np.linalg.matrix_rank(A.toarray())` on the unpinned version) and the row of a corner cell; `st = MAC.projection_stages(us, vs, g,
    dt=0.01)` on the deterministic divergent field; print `st["rhs_sum"], np.abs(st["div_before"]).max(),
    np.abs(st["div_after"]).max(), np.abs(st["curl_correction"]).max()`; then `MAC.mac_projection_summary(8, "divergent")`.
    *expect:* 64 × 64, rank 63 unpinned; the corner row has 3 nonzero entries (two neighbours dropped); rhs_sum ≈ 1e-15; div
    after < 1e-12; curl of the correction < 1e-12. *explain:* 1. the matrix of (10.124); 2. the four stages; 3. the parity summary
    E6 uses.
84. `nb.check_agree` — **from scratch (curation §7):** the 2-D Neumann Laplacian from `sparse.diags` + `sparse.kron` (1-D Neumann
    blocks with −1 in the corners, P244), one pressure pinned, the divergence of the staggered field by array slicing (P242), the
    solve with `spsolve`, the face correction by slicing, and `assert np.abs(div_after).max() < 1e-12`; `assert np.allclose(p_mine -
    p_mine.mean(), st["p"] - st["p"].mean())`.
85. `nb.primer("null space by SVD (count the tiny singular values)", "Every matrix factors as $A=U\\Sigma V^T$ (singular value
    decomposition); the columns of V whose singular values are (numerically) zero span the null space — the inputs A cannot see.
    `np.linalg.svd(A)` returns the singular values; count those below a tolerance (Ch. 1 P58 did null spaces symbolically).",
    code="import numpy as np\nA = np.array([[1., -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1]])   # differences of neighbours\ns = np.linalg.svd(A, compute_uv=False)\nprint(s, (s < 1e-12).sum() + (A.shape[1] - len(s)))   # one invisible input: the constant")` (**P245**)
86. `nb.derivation("D21", …)` — Part F D21 (7 steps), ref "10.126". (Uses D01 stencils, P245.)
87. `nb.note` — **N68 [B], N69 [B]** "**Collocated vs staggered gradient (10.125)–(10.126):** on a normal grid the centred pressure
    gradient $\\big(\\frac{\\partial p}{\\partial x}\\big)_{i,j}=\\frac{p_{i+1,j}-p_{i-1,j}}{2\\Delta x}$ (10.125) skips the neighbour, so a zigzag
    $p=(-1)^{i+j}$ has zero gradient everywhere — the momentum equation feels it as a uniform pressure. On the staggered grid the face
    gradient $\\big(\\frac{\\partial p}{\\partial x}\\big)_{i+1/2,j}=\\frac{p_{i+1,j}-p_{i,j}}{\\Delta x}$ (10.126) sees ±2/Δx and is second
    order at the face. Likewise the collocated divergence cannot see a zigzag velocity." equations (10.125), (10.126) + `nb.code`:
    `p = ch10.checkerboard(8, 8)`; `print(np.abs(ch10.collocated_gradient(p, 1, 1)[0][1:-1, 1:-1]).max(),
    np.abs(MAC.gradient(p, MAC.MacGrid(8, 8))[0][:, 1:-1]).max())`; `print(ch10.gradient_null_space(8, 8, "collocated"),
    ch10.gradient_null_space(8, 8, "staggered"))`. *expect:* 0.0 and 2.0; null-space dimensions ≥ 4 and 1 (only the constant).
88. `nb.figure` — **"The checkerboard hides on a collocated grid"** (8 × 3.4 in, two panels): (a) the zigzag pressure as a
    heatmap (orange/white) with the collocated gradient arrows (none — a caption "all zero"); (b) the same pressure on the staggered
    grid with face-gradient arrows ±2/Δx (orange) alternating. *see:* "an empty panel and a panel full of arrows"; *read:* "on the
    staggered grid the momentum equation pushes back against a zigzag, so it cannot grow"; *change:* "…mixed finite elements with
    equal order (C13): the same blindness appears, cured the same way — give the velocity more freedom".
89. `nb.note` — **N71 [B]** "**One MAC time step:** predictor (10.119)–(10.120) → Poisson (10.124) → correction (10.121)–(10.122).
    The Poisson solve is the costliest step." ASCII pipeline + `nb.code`: time each stage of `MAC.step` on a 64² grid (a
    `time.perf_counter` loop over 50 steps). *expect:* Poisson solve ≳ half the time with `splu` (the cell prints the shares;
    the builder states whatever it measures).
90. `nb.note` — **N73 [B]** "**MAC stability limits (Peyret & Taylor, Δx = Δy) (10.127)–(10.128):** $\\frac12(u^2+v^2)\\Delta t\\,Re\\le1$
    (10.127) and $\\frac{4\\Delta t}{Re\\,\\Delta x^2}\\le1$ (10.128) — the first is convection's, the second diffusion's (the 2-D β ≤ ¼).
    Stated, not derived (★★★): checked by a 2-D von Neumann scan and by a Taylor–Green run at 0.95 and 1.05 times the limit. For the
    cavity at Re = 100: Δx = 1/32 → Δt ≤ 0.02; Δx = 1/64 → Δt ≤ 0.0061." equations (10.127), (10.128) + `nb.code`:
    `MAC.dt_limit(1.0, 0.0, 100.0, 1/32, safety=1.0), MAC.dt_limit(1.0, 0.0, 100.0, 1/64, safety=1.0)`; `FD.ftcs2d_max_amplification(1.0,
    0.0, 100.0, 1/32, dt)` for dt = 0.95 and 1.05 × the limit. *expect:* 0.02, 0.0061035; ≤ 1 and > 1.
91. `nb.code` — **verification cell:** `tg = MAC.taylor_green(32, 100.0, 0.5)` and at n = 16, 64 → `observed_order`; `cp =
    MAC.channel_poiseuille(16, 16, 10.0, -1.0)`. *expect:* spatial order ≈ 2 (with Δt ∝ Δx², the projection's first-order time error
    kept small); Poiseuille max error < 1e-12 (exact for the 3-point Laplacian).
92. `nb.explainer("mac_projection_staggered", heading="What is the pressure doing in incompressible flow?", why="Stepping the
    transport through predictor → Poisson → correction shows the divergence appear and vanish to 10⁻¹⁵ on the same small grid;
    the collocated/staggered switch makes the checkerboard's invisibility visible; clicking a cell adds up its four face fluxes.",
    tries=["Step through the four stages and watch the rose divergence blotches vanish at stage 4.", "Click the cell next to a
    wall before and after the projection: which face never changes?", "Switch to the checkerboard mode: count the arrows on the
    collocated side."])`
93. `nb.md` — **Climate hook:** "The MAC staggering is the **Arakawa C-grid**: u and v on the faces, the sea-surface height (or
    pressure, or layer thickness) at the centres — used by NEMO, MOM, MITgcm and many atmosphere cores; Ch. 13's Kelvin and Rossby
    waves are computed on exactly this layout. The projection is the pressure solve of every Boussinesq ocean model (Ch. 4 §4.9's
    Boussinesq approximation)."
94. `nb.md` — **What would change if…** "…we used finite elements for the same problem? The zigzag has an FE twin — spurious
    pressure modes — and the cure is again to give the velocity more freedom than the pressure. C13."
95. `nb.recap("R10", "The rate-of-strain tensor", "$\\mathbf D[\\mathbf u]=\\frac12[\\nabla\\mathbf u+(\\nabla\\mathbf u)^T]$ (10.136) is
    the symmetric part of the velocity gradient of Ch. 3 (S_ij, `core.kinematics`); D22 uses its symmetry.", where="Ch. 3 §3.4")`

#### C13 — Mixed finite elements and the LBB (inf–sup) condition
96. `nb.core("C13", "Mixed finite elements: the saddle-point system $\\begin{pmatrix}\\mathbf M\\dot{\\mathbf u}\\\\\\mathbf 0\\end{pmatrix}+
    \\begin{pmatrix}\\mathbf A&\\mathbf B\\\\\\mathbf B^T&\\mathbf 0\\end{pmatrix}\\begin{pmatrix}\\mathbf u\\\\\\mathbf p\\end{pmatrix}=
    \\begin{pmatrix}\\mathbf f_u\\\\\\mathbf f_p\\end{pmatrix}$ (10.137) and the LBB condition (Taylor–Hood P2–P1)", question="Why
    can't velocity and pressure be built from the same simple elements?")`
97. `nb.md` — **The problem in plain words:** "Build a finite-element Navier–Stokes solver with the same straight-sided tents for
    velocity and pressure — the obvious first try. The velocity looks fine, but the pressure comes out as noise: a mottled pattern
    that changes wildly from node to node. The reason is counting. Every pressure unknown adds one continuity constraint; if the
    velocity space is not rich enough to satisfy all of them, some pressure patterns are left that the discrete divergence simply
    cannot see — the FE version of the checkerboard. Give the velocity quadratic elements (six nodes per triangle) and the pressure
    linear ones (three), and the noise disappears."
98. `nb.md` — **The idea** (ASCII + counting table):
    ```
          ●             P2–P1 (Taylor–Hood): velocity at 6 nodes (● vertices, ○ mid-edges), pressure at 3 (●)
         ○ ○            iso-P2/P1: split the triangle into 4, linear velocity on the small ones, pressure on the big one
        ●──○──●
    ```
    | mesh (n × n squares, all walls) | velocity unknowns P1 | pressure (−1 pinned) | velocity unknowns P2 |
    | n = 2 | 2 | 8 | 18 |
    | n = 4 | 18 | 24 | 98 |
    "When pressure unknowns outnumber velocity unknowns, some pressure modes are certainly invisible."
99. `nb.note` — **N75 [B]** "**The weak Navier–Stokes equations (10.134)–(10.135):** $\\int_\\Omega\\Big(\\frac{\\partial\\mathbf u}{\\partial t}
    +\\mathbf u\\cdot\\nabla\\mathbf u-\\mathbf g\\Big)\\cdot\\tilde{\\mathbf u}\\,d\\Omega+\\frac2{Re}\\int_\\Omega\\mathbf D[\\mathbf u]:
    \\mathbf D[\\tilde{\\mathbf u}]\\,d\\Omega-\\int_\\Omega p\\,(\\nabla\\cdot\\tilde{\\mathbf u})\\,d\\Omega=0$ (10.134) and
    $\\int_\\Omega\\tilde p\\,\\nabla\\cdot\\mathbf u\\,d\\Omega=0$ (10.135), with velocity and pressure variations ũ, p̃ — C06's recipe in
    2-D." equations (10.134), (10.135).
100. `nb.md` — **🔁 Reminder:** the divergence theorem (Ch. 2) and the double dot product $\\mathbf A:\\mathbf B=A_{ij}B_{ij}$ (Ch. 2).
101. `nb.derivation("D22", …)` — Part F D22 (11 steps), ref "10.134". (Uses D09 moves, R10, the divergence theorem, P113 product rule
     for a divergence.)
102. `nb.primer("saddle-point (KKT) matrix", "A block matrix $\\begin{pmatrix}\\mathbf A&\\mathbf B\\\\\\mathbf B^T&\\mathbf 0\\end{pmatrix}$ —
     the shape of every 'minimise subject to a constraint' problem (the Lagrange multiplier of P239 is the second block). It has
     positive and negative eigenvalues (indefinite), a zero block on the diagonal, and is invertible only if $\\mathbf B$ has full
     column rank — no multiplier pattern may be invisible to the constraint.", code="import numpy as np\nA = np.eye(2); B = np.array([[1.], [1.]])            # 2 velocities, 1 pressure\nK = np.block([[A, B], [B.T, np.zeros((1, 1))]])\nprint(np.linalg.eigvalsh(K))                         # [-1, 1, 2]: one negative eigenvalue\nB2 = np.array([[1., 1.], [1., 1.]])                  # 2 pressures that the constraint cannot tell apart\nprint(np.linalg.matrix_rank(np.block([[A, B2], [B2.T, np.zeros((2, 2))]])))   # 3 < 4: singular")`
     (**P246**)
103. `nb.primer("GMRES in one line", "For large nonsymmetric systems (like the Newton steps of the cylinder problem) direct
     elimination is too costly; GMRES (Saad) builds the best approximation in the growing space spanned by b, Ab, A²b, … and stops
     when the residual is small. `scipy.sparse.linalg.gmres(A, b)` returns (x, info); info = 0 means converged.",
     code="import numpy as np, scipy.sparse as sp\nfrom scipy.sparse.linalg import gmres\nA = sp.diags([-1, 2.5, -1.2], [-1, 0, 1], shape=(50, 50)).tocsr()   # nonsymmetric, well conditioned\nx, info = gmres(A, np.ones(50))\nprint(info, np.abs(A @ x - 1).max())                 # 0 and a small residual")` (**P250**)
104. `nb.note` — **N76 [B]** "**The semi-discrete system (10.137):** velocity unknowns **u** and pressure unknowns **p**; **M** the mass
     matrix, **A** depends on u (convection), **B** comes from the pressure term of (10.134) and **Bᵀ** from (10.135) — the symmetric
     arrangement is why the book writes (10.135) with that sign. Discretised in time and linearised by Newton, each step is a sparse
     saddle-point solve — direct elimination for small systems, GMRES for large ones." equation (10.137) + `nb.figure` — the sparsity
     pattern (`plt.spy`) of the Stokes P2–P1 matrix for n = 4 with the A, B, Bᵀ, 0 blocks outlined. *see:* "a big A block, thin B
     strips, an empty corner"; *read:* "the empty corner is the missing 'pressure equation' — continuity has no p in it"; *change:*
     "…GLS stabilisation: a small negative-definite block appears in the corner".
105. `nb.primer("generalised symmetric eigenproblem scipy.linalg.eigh(A, B)", "Solve $\\mathbf Ax=\\lambda\\mathbf Bx$ for symmetric A and
     positive-definite B: the eigenvalues measure A 'in units of' B (Ch. 2 P80 was the case B = I). The discrete inf–sup constant is
     $\\beta_h^2=\\lambda_{\\min}(\\mathbf B^T\\mathbf A^{-1}\\mathbf B,\\ \\mathbf M_p)$ over pressures with zero mean: how strongly the
     worst pressure pattern is felt by the velocity.", code="import numpy as np\nfrom scipy.linalg import eigh\nA = np.array([[2., 0], [0, 1]]); B = np.array([[1., 0], [0, 4]])\nprint(eigh(A, B, eigvals_only=True))            # [0.25, 2.0]: A measured in units of B")` (**P247**)
106. `nb.md` — **The LBB (inf–sup) condition, in words and one formula:** "A pair of spaces is stable if every pressure pattern q
     can be 'felt' by some velocity v with a force not vanishing on refinement: $\\inf_q\\sup_{\\mathbf v}\\frac{\\int q\\,\\nabla\\cdot\\mathbf
     v\\,d\\Omega}{\\lVert q\\rVert\\,\\lVert\\mathbf v\\rVert_1}\\ge\\beta>0$ independent of h (Babuška–Brezzi; the book states it by name). Equal
     order (P1–P1) fails; P2–P1 (Taylor–Hood, Fig. 10.5a) and iso-P2/P1 (Fig. 10.5b) pass; GLS stabilisation (Tezduyar; Franca & Frey)
     keeps equal order by adding least-squares terms — the mixed method needs no tuning parameter, which the book prefers."
107. `nb.primer("barycentric (area) coordinates on a triangle", "Any point of a triangle is a weighted average of its three corners
     with weights (ζ, ξ, η) ≥ 0 that add to 1 — on the parent triangle ζ = 1 − ξ − η. Each weight is 1 at 'its' corner and 0 on the
     opposite side, so they are the linear shape functions (10.187), and products like 4ξζ are the quadratic mid-edge ones (10.185).",
     code="xi, eta = 0.2, 0.3                     # a point inside the parent triangle\nzeta = 1 - xi - eta                     # 0.5\nprint(zeta + xi + eta)                  # 1.0: a partition of unity\nprint(4*xi*zeta, 4*xi*eta, 4*eta*zeta)  # mid-edge P2 shapes at this point: 0.4, 0.24, 0.6")` (**P248**)
108. `nb.primer("Jacobian determinant of a 2-D map", "A map (ξ, η) → (x, y) stretches a tiny parent square of area dξ dη into a
     parallelogram of area J dξ dη with $J=x_\\xi y_\\eta-x_\\eta y_\\xi$ — the 2-D version of dx = (h/2) dξ (P234). For a straight-sided
     triangle J = 2 × its area (the parent triangle has area ½); for a curved (isoparametric) one J varies inside.",
     code="import numpy as np\nxe, ye = np.array([0., 2, 0]), np.array([0., 0, 1])     # corners of a triangle with area 1\nxxi, xeta = xe[1] - xe[0], xe[2] - xe[0]; yxi, yeta = ye[1] - ye[0], ye[2] - ye[0]\nprint(xxi*yeta - xeta*yxi)                               # J = 2 = 2 × area")` (**P249**)
109. `nb.worked_example("counting unknowns on a 2 × 2 mesh", "The unit square cut into 2 × 2 squares, each split into two triangles:
     V = 9 vertices, E = 16 edges, T = 8 triangles (V − E + T = 1 ✓). All walls have prescribed velocity. 1. P1–P1: velocity nodes =
     vertices; only the centre vertex is interior ⇒ 2 velocity unknowns. Pressure: 9 values, minus 1 pinned ⇒ 8. Eight constraints,
     two unknowns: at least 6 pressure patterns cannot be felt — spurious modes. 2. P2–P1: velocity nodes = vertices + edge midpoints
     = 25, interior (2n − 1)² = 9 ⇒ 18 unknowns; pressure still 8. 3. Counting is necessary, not sufficient: at n = 8 P1–P1 has 98
     velocity and 80 pressure unknowns, yet β_h still falls with h — the inf–sup test below decides.")`
110. `nb.code` — `for pair in ("P1P1", "P2P1"): for n in (2, 4, 8): r = FEM2.infsup_constant(n, pair); print(pair, n, r["n_u"],
     r["n_p"], r["n_spurious"], r["beta"])`; `for pair in ("P1P1", "P2P1"): s = FEM2.stokes_cavity(6, pair); print(pair,
     s["p_range"], s["p_std_checker"])`. *expect:* P1P1 n = 2: n_u 2, n_p 8, spurious ≥ 6, β = 0; P2P1: spurious 0 and β roughly
     constant across n (the cell prints the values; the assert is β(P2P1, 8) > 0.5 β(P2P1, 2) and β(P1P1, 8) < 0.5 β(P2P1, 8));
     the P1P1 cavity pressure has a large alternating component, P2P1's does not. *explain:* 1. counts from the mesh; 2. the
     generalised eigenvalue (P247); 3. the Stokes cavity.
111. `nb.check_agree` — **from scratch (curation §7):** the six P2 shape functions (10.185) typed from the formula and evaluated at the
     six nodes (vertices (0, 0), (1, 0), (0, 1) and mid-edges) → the identity matrix; partition of unity at random points; the 7-point
     rule (points and weights from `FEM2.tri_quad_7pt`) integrating ξ²η² over the parent triangle: `assert np.isclose(½ Σ W f,
     2!·2!/6!)` = 1/180; and ξ³η³ (degree 6) *not* exact. `assert np.allclose(phi_mine, FEM2.p2_shape(xi, eta))`.
112. `nb.figure` — **Fig. 10.5 remade and the spurious pressure** (10 × 3.4 in, three panels): (a) `p2p1_triangle(ax)` and
     `p2p1_triangle(ax, iso=True)` (the two element types); (b) the Stokes-cavity pressure for P1–P1 (n = 6) and P2–P1 (n = 6) as
     coloured triangulations side by side (orange scale); (c) β_h vs n on log axes for both pairs (P1P1 rose, P2P1 blue). *see:*
     "a mottled pressure and a smooth one; one falling curve, one flat"; *read:* "the flat curve is the LBB condition holding"; *change:*
     "…iso-P2/P1: its curve is flat too, at a slightly lower level".
113. `nb.plotly` — **F7** `slider_figure` over n (2 … 16 from `FEM2.infsup_table`, read from `reference/ch10/infsup_table.csv`):
     β_h for both pairs with the current n marked. *explain:* "the LBB trend without running anything".
114. `nb.explainer("mixed_fe_lbb", heading="Why can't velocity and pressure use the same elements?", why="Switching the element pair
     and refining the mesh shows the Stokes-cavity pressure change from spurious noise to a smooth field while the inf–sup constant
     falls or stays flat, and the node view counts unknowns as they grow — the trend with refinement is the whole point.",
     tries=["Press 'P1–P1, n = 4': count velocity and pressure unknowns in the status.", "Switch to P2–P1 at the same n: watch the
     pressure noise disappear.", "Open one element: click a triangle to see its six velocity nodes and three pressure nodes."])`
115. `nb.md` — **Where C13 goes next:** "§10.5's third example — flow past a cylinder in a channel — is computed with exactly these
     P2–P1 elements, curved to fit the cylinder, and Newton's method for the nonlinear term: 'C13 (continued)' below." **What would
     change if…** "…we now wanted evidence that any of these solvers is right? §10.5 answers with a benchmark (C14) and a
     grid-convergence study (C15)."

---

### A.5 §10.5 Three Examples — C14 · R11, C15 · N92 · C13 (continued)
1. `nb.section("10.5", "Three Examples", intro="**What is this section about?** Three computations put the methods to work: a
   square cavity whose lid slides (the most-used benchmark in CFD), a square block in a channel (drag, lift, shedding — and the
   question every simulation must answer: does the answer change if the grid is refined?), and a cylinder in a channel computed
   with the mixed finite elements of C13. We reproduce the first with our own solvers against published data, study grid
   convergence on the first two, and show the third from our own (smaller, cached) runs. All numbers here are ours; the book's
   own values are kept aside for private tests.")`

#### C14 — The lid-driven cavity benchmark
2. `nb.core("C14", "The lid-driven cavity benchmark: $u(0.5,y)$ against Ghia, Ghia & Shin (1982) and the primary-eddy centre",
   question="How can we tell whether a flow solver solves the equations we meant — before trusting it on a problem nobody has
   solved?")`
3. `nb.note` — **N77 [C]** "**The three examples** of §10.5: cavity and block computed with explicit MacCormack on weakly compressible
   Navier–Stokes (C10), the cylinder with mixed finite elements (C13). We add the MAC projection solver (C11–C12) for the cavity."
4. `nb.md` — **The problem in plain words:** "Fill a square box with a viscous liquid and drag the lid sideways at constant speed.
   The liquid starts to turn: one big eddy fills the box, with small counter-rotating eddies tucked into the bottom corners. The flow
   is simple to set up (no inflow, no outflow, one number Re), yet has no formula — so dozens of groups have computed it on very
   fine grids and published their numbers. If your code reproduces those numbers, you have **evidence** (not proof) that it solves
   the Navier–Stokes equations you meant."
5. `nb.md` — **The idea** (ASCII):
   ```
        lid u = 1 ────────────►
       ┌───────────────────────┐      the lid drags the top layer; continuity turns it down the right wall,
       │      ↺ primary eddy   │      back along the floor and up the left wall: one big eddy
       │  (centre up and to the │      corner eddies ↻ in the bottom corners (much weaker)
       │   right of the middle) │      benchmark: u(0.5, y) on the vertical centreline + the eddy centre
       │↻                     ↻│
       └───────────────────────┘
   ```
6. `nb.primer("verification versus validation", "*Verification*: are we solving the equations right? (code bugs, grid and time-step
   errors — checked against exact solutions, manufactured solutions and other accurate numerical solutions such as Ghia's).
   *Validation*: are we solving the right equations? (checked against experiments). A benchmark computation like the cavity is
   verification; Ch. 12's comparisons with turbulence measurements are validation.", code="checks = {'exact solution (Poiseuille, Taylor-Green)': 'verification', 'Ghia et al. 1982 cavity data': 'verification', 'wind-tunnel drag of a real cylinder': 'validation'}\nfor k, v in checks.items(): print(f'{v:12s} <- {k}')")` (**P253**)
7. `nb.note` — **N78 [B]** "**Set-up and scales (Fig. 10.6 analogue):** square of side L (the book calls it D), lid speed U; lengths
   in L, velocities in U, time in L/U, density in ρ₀, pressure in ρ₀U²; for the weakly compressible solver the dimensionless pressure
   is $p=\\rho/Ma^2$ with $Ma=U/c$, and $Re=\\rho_0UL/\\mu$. The two top corners are singular (the lid velocity jumps from 1 to 0) — the
   grid smooths the jump; results near the corners converge slowly." + `nb.figure` — `cavity_sketch(ax)` (our drawing: lid arrow,
   walls, the centreline x = ½ dashed, the eddy cartoon). *see/read/change:* "the lid is the only forcing"; "the centreline cuts the
   primary eddy through its core"; "…a deep cavity (height 2): a second eddy stacks below the first".
8. `nb.note` — **N79 [B], N80 [B], N81 [B], N82 [B], N83 [B], N84 [B]** "**Density on the walls from continuity (10.138)–(10.146)** —
   the weakly compressible solver needs ρ on the walls, and continuity supplies it. On the left wall v = 0 along the wall, so
   $\\frac{\\partial\\rho}{\\partial t}+\\frac{\\partial(\\rho u)}{\\partial x}=0$ (10.138); ∂(ρu)/∂x is taken one-sided and second order,
   $\\big(\\frac{\\partial f}{\\partial x}\\big)_i=\\frac1{2\\Delta x}(-f_{i+2}+4f_{i+1}-3f_i)+O(\\Delta x^2)$ (N80, weights (−3/2, 2, −1/2) by C01's
   `FD.fd_weights([0, 1, 2], 1)`), giving the predictor $\\rho^*_{i,j}=\\rho^n_{i,j}-\\frac{\\Delta t}{2\\Delta x}[-(\\rho u)^n_{i+2,j}+
   4(\\rho u)^n_{i+1,j}-3(\\rho u)^n_{i,j}]$ (10.139) and the corrector (10.140). The right wall (10.141)–(10.142) mirrors it (backward
   stencil, sign flips); the bottom (10.143)–(10.144) is the y-version; on the lid u = U, so $\\frac{\\partial(\\rho u)}{\\partial x}=
   U\\frac{\\partial\\rho}{\\partial x}$, taken centred (10.145)–(10.146)." equations (10.138), (10.139) + `nb.code`:
   `FD.fd_weights([0, 1, 2], 1)`, `FD.fd_weights([0, -1, -2], 1)`. *expect:* (−3/2, 2, −1/2) and (3/2, −2, 1/2).
9. `nb.note` — **N85 [B]** "**The six-substep MacCormack cavity algorithm** (as a pipeline): 1 velocities from momenta → 2 interior
   predictor with coefficients a₁ … a₁₁ ($a_1=\\frac{\\Delta t}{\\Delta x}$, …, $a_5=\\frac{4\\Delta t}{3Re\\Delta x^2}$, $a_9=\\frac{\\Delta
   t}{12Re\\Delta x\\Delta y}$, $a_{10}=2(a_5+a_6)$, $a_{11}=2(a_7+a_8)$) → 3 wall densities at $t_{n+1}$ → 4 starred velocities → 5
   corrector → 6 wall densities. > ⚠️ **The book prints Step 5 as** $2\\rho^{n+1}=(\\rho^n+\\rho^*)-a_1+[(\\rho u)^*_{i,j}-(\\rho u)^*_{i-1,j}]-
   \\dots$; **the correct form is** $-a_1[(\\rho u)^*_{i,j}-(\\rho u)^*_{i-1,j}]$ as in (10.106) (slip R5 — the printed version does not
   conserve mass; `printed_step5=True` shows it drifting)." + `nb.code`: `MCK.cavity_coefficients(1e-3, 1/32, 1/32, 0.08, 100.0)`.
10. `nb.primer("reading reference data from a file and interpolating (np.interp)", "Published benchmark numbers live in a small
    text table (`reference/ch10/ghia1982_table1.csv`, with the citation in its header). `np.loadtxt(path, delimiter=',',
    skiprows=…)` or `pandas.read_csv` reads it; to compare with our grid, interpolate our profile at *their* points with
    `np.interp(their_y, our_y, our_u)` (Ch. 7 P182) — never the other way round, so no benchmark value is invented.",
    code="import numpy as np\nour_y = np.linspace(0, 1, 11); our_u = our_y**2           # a stand-in profile on our grid\ntheir_y = np.array([0.0625, 0.5, 0.9531])                  # three points where a benchmark reports values\nprint(np.interp(their_y, our_y, our_u))                     # our values at their points")` (**P251**)
11. `nb.primer("caching expensive runs (np.savez and a parameter key)", "A converged 128² cavity takes about a minute; re-running it
    at every notebook execution wastes time. Save the result once — `np.savez(path, u=u, v=v, ...)` in `outputs/ch10/` with the
    parameters in the file name (Re, n, …) — and load it next time with `np.load(path)`. `functools.lru_cache` does the same in
    memory within one session. The solvers here do it for you (`cache=True`); the notebook says whether a result was loaded or
    computed.", code="import numpy as np, pathlib, tempfile\np = pathlib.Path(tempfile.gettempdir()) / 'demo_Re100_n16.npz'\nif not p.exists(): np.savez(p, u=np.arange(3.0))   # 'expensive' run happens once\nprint(np.load(p)['u'])                                # later runs just load it")` (**P252**)
12. `nb.worked_example("sizes and step counts for a Re = 100 cavity", "Water-like honey: L = 0.1 m, U = 0.1 m/s, ν = 10⁻⁴ m²/s ⇒
    Re = UL/ν = 100. 1. MAC on 64 × 64 cells (Δx = 1/64): the diffusion limit (10.128) gives Δt ≤ ReΔx²/4 = 100/(4 × 4096) = 0.0061,
    the convection limit (10.127) Δt ≤ 2/(Re·1²) = 0.02; take Δt = 0.8 × 0.0061 = 0.0049 ⇒ about 4100 steps to t = 20 (20 lid-transit
    times, enough to settle). 2. MacCormack with our Ma = 0.08 (sound speed c = 1/Ma = 12.5 lid speeds): (10.155) gives Δt ≤ (σ/√2) ×
    0.08/64 = 7.1 × 10⁻⁴ with our σ = 0.8 ⇒ 28 000 steps — seven times more, the price of following sound waves. 3. Density error ~ Ma² = 0.6 %.")`
13. `nb.code` — `st = MAC.cavity(Re=100.0, n=64 if not FAST else 32)` (cached); `cl = MAC.cavity_centreline(st)`; `gh =
    ch10.ghia_centreline(100)`; `err = ch10.cavity_error_vs_ghia(st, 100)`; `c = MAC.primary_vortex_centre(st["psi"], st["g"])`; print
    `st["steps"], st["converged"], err`, `c`, `ch10.ghia_vortex_centre(100)`, `ch10.hou_centres()[100]`; then `mc =
    MCK.cavity_maccormack(Re=100.0, Ma=0.08, n=32)` (cached; 64 if a cache from `scripts/ch10_cavity.py` exists) and its error.
    *expect:* converged True; MAC 64²: rel_max ≤ 0.02 (≤ 0.04 at 32²); centre within one cell of Ghia's (the cell prints both — Ghia's
    loaded from the file, never typed); MacCormack 32²: rel_max ≤ 0.05. *explain:* 1. the MAC run to steady state (predictor →
    Poisson → correction, 4100 steps); 2. the centreline and the benchmark from the file; 3. the eddy centre (ψ minimum); 4. the
    weakly compressible run for comparison.
14. `nb.check_agree` — **from scratch (curation §7):** interpolate our centreline at Ghia's y-points with `np.interp` (P251) and compute
    the maximum deviation by hand: `dev = np.abs(np.interp(gh["y"], cl["y"], cl["u"]) - gh["u"]).max()`; `assert np.isclose(dev,
    err["max_dev"])`.
15. `nb.figure` — **Figs. 10.7–10.8 remade** (10 × 4 in, two panels): (a) streamlines (ψ contours, blue, with the two corner eddies
    in teal) over a vorticity heatmap for our 64² MAC run, the lid arrow, our eddy centre (✕) and Ghia's (black dot); (b) u(0.5, y):
    MAC 32² (muted), MAC 64² (blue), MacCormack 32² (purple dashed), Ghia's points (black dots with white halo). *see:* "one big eddy
    and two tiny ones; curves that pass through the dots"; *read:* "the 64² curve lies on the benchmark within 2 %; the coarse curves
    miss the minimum near y ≈ 0.45 most"; *change:* "…Re = 400: the eddy centre moves toward the middle and the minimum of u deepens
    (a cached run, if present, is overlaid)".
16. `nb.animation` — **A4** (video, 60 frames; FAST 30): the cavity spinning up from rest at Re = 100 (64², snapshots every 0.33 time
    units from the cached run's `snapshots`): streamlines + vorticity (left) and the centreline profile approaching Ghia's dots
    (right). *explain:* "a steady benchmark is the end of a transient — watch the profile settle onto the dots".
17. `nb.explainer("lid_driven_cavity", heading="How do I know my CFD answer is right?", why="Running the MAC solver live on 16², 24² and
    32² grids (and loading cached 64²/128² and MacCormack results), you watch the eddy form from rest and the centreline approach Ghia's
    points, while the error-vs-h view fills in dot by dot and the Richardson estimate updates.", tries=["Run 16², then 32²: how much
    closer is the minimum to Ghia's?", "Load the cached 128² run and read the observed order in the status.", "Switch Re to 400: which
    part of the profile changes most?"])`
18. `nb.md` — **What would change if…** "…there were no benchmark at all — a new geometry, a new Re? Then the only evidence left is
    internal: refine the grid and show that the answer stops changing at the rate the scheme promises. C15."
19. `nb.recap("R11", "Force coefficients and a steady-to-shedding wake", "Ch. 4 defined $C_D=\\frac{F_D}{\\frac12\\rho U^2A}$ (4.107) and
    $C_L=\\frac{F_L}{\\frac12\\rho U^2A}$ (4.108) (here per unit span with A = the block side); Ch. 9 §9.8 showed a bluff body's wake is
    steady at low Re and sheds vortices above a threshold. For the block in a channel we compute C_D(t), C_L(t) at Re = 20 (steady after
    an acoustic transient) and Re = 100 (periodic shedding) — coarse runs here, finer ones cached by `scripts/ch10_block.py`.",
    where="Ch. 4 §4.11, Ch. 9 §9.8")`

#### C15 — Grid convergence, Richardson extrapolation and the verification checklist
20. `nb.core("C15", "Grid convergence and Richardson extrapolation: $f_{h\\to0}\\approx f_1+\\frac{f_1-f_2}{r^p-1}$ from the error model
    $\\lVert e\\rVert\\le K_e\\Delta x^a$ (10.15)", question="When no exact answer exists, how can three runs on three grids tell us both
    that the code works and what the grid-independent answer is?")`
21. `nb.note` — **N86 [B]** "**A square block in a channel (Fig. 10.9 analogue, our geometry):** block side d = 1, channel height
    H = 4 (blockage 1/4), 8 sides upstream, 20 downstream, walls sliding at the stream speed (the block moves through still fluid —
    the same flow seen from the block, Ch. 3's Galilean change of frame), uniform inflow, zero-gradient outflow, Ma = 0.06, Re = 20
    and 100 on the block side." + `nb.figure` — `block_channel_sketch(ax)`; *see/read/change:* "a small square in a long channel";
    "the walls move so no boundary layer grows on them"; "…fixed walls: the channel's own Poiseuille profile would meet the block".
22. `nb.note` — **N87 [B], N88 [B], N89 [B]** "**Density on the block's faces (10.147)–(10.154)** (heuristic, as the book says 'it gives
    better results'): on the front face u = v = 0, so the convective terms and the time derivative of the x-momentum equation vanish and
    $\\frac{\\partial\\rho}{\\partial x}=\\frac{Ma^2}{Re}\\Big(\\frac43u_{xx}+\\frac13v_{xy}\\Big)$ (10.147) — only viscous terms set the
    density gradient; one-sided stencils $\\big(\\frac{\\partial\\rho}{\\partial x}\\big)_{i,j}=\\frac{-1}{2\\Delta x}(-\\rho_{i-2,j}+4\\rho_{i-1,j}-
    3\\rho_{i,j})$ (10.148), $(u_{xx})_{i,j}=\\frac1{\\Delta x^2}(2u_{i,j}-5u_{i-1,j}+4u_{i-2,j}-u_{i-3,j})$ (10.149) (weights (2, −5, 4, −1)
    recomputed by `FD.fd_weights([0, -1, -2, -3], 2)`) and a one-sided mixed derivative (10.150); solving for the wall value gives the
    front-face density (10.151) $\\rho_{i,j}=\\frac13(4\\rho_{i-1,j}-\\rho_{i-2,j})+\\frac{8}{9\\Delta x}\\frac{Ma^2}{Re}(-5u_{i-1,j}+4u_{i-2,j}-
    u_{i-3,j})-\\frac1{18\\Delta y}\\frac{Ma^2}{Re}\\big[-(v_{i-2,j+1}-v_{i-2,j-1})+4(v_{i-1,j+1}-v_{i-1,j-1})-3(v_{i,j+1}-v_{i,j-1})\\big]$ —
    8/9 = (2/3)(4/3) and 1/18 = (2/3)(1/3)(1/4); the back, top and bottom faces (10.152)–(10.154) by symmetry; corners averaged."
    equations (10.147), (10.151) + `nb.code`: `FD.fd_weights([0, -1, -2, -3], 2)`. *expect:* (2, −5, 4, −1).
23. `nb.note` — **N90 [B]** "**Practice:** long explicit runs accumulate round-off, so use double precision (P222), compute the density
    perturbation ρ′ = ρ − 1 instead of ρ (fewer digits cancel), and cycle the MacCormack arrangements (FB/BF then BF/FB, N57)." +
    `nb.code`: 20 000 steps of `MCK.weakly_compressible_step` on a uniform state (nothing should happen) in float32 and float64; print
    max∣ρ − 1∣. *expect:* float64 ~ 1e-13, float32 ~ 1e-4 or larger — the drift is pure round-off.
24. `nb.note` — **N91 [B]** "**The asymptotic MacCormack step (10.155):** for Δx = Δy and large grid Reynolds number (10.110) becomes
    $\\Delta t\\le\\frac{\\sigma}{\\sqrt2}\\,Ma\\,\\Delta x$. > ⚠️ **The book states it follows 'with large grid Reynolds numbers'; it also needs
    $Ma\\ll1$** — so that $\\lvert u\\rvert/\\Delta x+\\lvert v\\rvert/\\Delta y\\ll\\sqrt2c/\\Delta x$ (slip R12). At Ma = 0.5 the two limits differ by
    ≈ 35 %." equation (10.155) + `nb.code`: `for Ma in (0.08, 0.5): print(MCK.maccormack_dt_asymptotic(Ma, 1/64),
    MCK.maccormack_dt(1.0, 0.0, 1/Ma, 1/64, 1/64, 1.0, 1e-6))` (μ tiny ⇒ large Re_Δ). *expect:* 7.07e-4 vs 6.69e-4 at 0.08 (within 6 %); 4.42e-3 vs 3.27e-3 at 0.5 (the asymptotic form ≈ 35 % too
    large).
25. `nb.md` — **The problem in plain words:** "You computed the drag on the block: C_D = 3.4 on one grid, 3.1 on a finer one. Which is
    right? Neither — both carry discretisation error. But if the scheme is second order, the error should fall by 4 each time the
    spacing halves. Three grids let you (1) *measure* the order (if it is not ≈ 2, something is wrong: a bug, or a grid too coarse to
    be in the 'asymptotic range'), and (2) *extrapolate* to zero spacing — Richardson's trick — giving a better number than any of the
    three runs and an honest error bar."
26. `nb.md` — **The idea** (ASCII):
    ```
    f(h) ≈ f₀ + K hᵖ          (the error model (10.15), one term)
    h₃ = 4h   f₃ ─┐
    h₂ = 2h   f₂ ─┤ differences shrink by rᵖ:   (f₃ − f₂)/(f₂ − f₁) = rᵖ  ⇒  p
    h₁ = h    f₁ ─┘ extrapolate:  f₀ ≈ f₁ + (f₁ − f₂)/(rᵖ − 1)
    ```
27. `nb.derivation("D23", …)` — Part F D23 (8 steps), ref "10.15". (Uses N10, P13 log–log slope.)
28. `nb.worked_example("Richardson on a made-up answer", "Suppose f(h) = 1 + h² exactly (so the true answer is 1 and p = 2). 1. Three
    grids h = 0.1, 0.05, 0.025: f₃ = 1.01, f₂ = 1.0025, f₁ = 1.000625. 2. Differences: f₃ − f₂ = 0.0075, f₂ − f₁ = 0.001875, ratio 4. 3.
    p = ln 4/ln 2 = 2 ✓. 4. Richardson: f₀ ≈ 1.000625 + (1.000625 − 1.0025)/(2² − 1) = 1.000625 − 0.000625 = 1.000000 ✓ — the leading
    error is removed exactly. 5. With a real code the next term (h³, h⁴) remains, so the extrapolated value is better but not exact.")`
29. `nb.code` — `gci = grid_convergence_index(1.000625, 1.0025, 1.01, r=2.0)` (the tiny example); then the cavity study: `umin =
    [MAC.cavity_centreline(MAC.cavity(Re=100.0, n=n))["u"].min() for n in ns]` with `ns = (16, 32, 64)` (FAST (12, 24, 48)); `gci =
    grid_convergence_index(umin[2], umin[1], umin[0], r=2.0)`; print p, f_ext, gci_fine, and Ghia's minimum from the file for comparison;
    then, if `outputs/ch10/block_cd_*.npz` exist, the same for the block's steady C_D at Re = 20 on Δx = 1/8, 1/16, 1/32 (else a
    one-line "run scripts/ch10_block.py --convergence" note). *expect:* tiny example p = 2.000, f_ext = 1.000000; cavity observed p
    between 1.5 and 2.5 (the cell prints it; the assert is monotone convergence and 1.2 < p < 2.8), extrapolated minimum within 1 % of
    Ghia's; block p ≈ 2 if cached (qualitative: the corners of the block limit the order — the book notes they do not spoil it much).
    *explain:* 1. the formula on known numbers; 2. three cavity grids; 3. the block if cached.
30. `nb.check_agree` — **from scratch (curation §7):** observed order and Richardson in four lines — `p = np.log((f3 - f2)/(f2 - f1))/
    np.log(r)`, `f0 = f1 + (f1 - f2)/(r**p - 1)` — against `grid_convergence_index` and `tools.convergence.richardson(f2, f1, r, p)`:
    `assert np.isclose(p, gci["p"]) and np.isclose(f0, gci["f_ext"]) and np.isclose(f0, richardson(f2, f1, 2.0, p))`.
31. `nb.figure` — **"Refine, measure the order, extrapolate"** (10 × 3.4 in, three panels): (a) the cavity's u_min vs h (dots, blue) with
    the Richardson estimate (purple star at h = 0) and Ghia's value (black line), plus the fitted f₀ + Khᵖ (muted dashed); (b) ∣u_min −
    f_ext∣ vs h on log–log with the slope printed; (c) the block's C_D(t) at Re = 20 on the coarse grid (orange; steady after the acoustic
    transient) and, if cached, the Fig. 10.13 analogue C_D vs Δx with its extrapolation. *see:* "three dots marching toward a line; a
    straight line on log–log"; *read:* "slope ≈ 2 = the scheme's order: the code is doing what it promises; the star is our best
    estimate, with the GCI as its error bar"; *change:* "…the three grids were too coarse (outside the asymptotic range): the
    differences would not shrink by a constant factor and p would be meaningless".
32. `nb.plotly` — **F6** `slider_figure` over the grid (16, 24, 32, 64, 128 — cached): the cavity centreline u(y) with Ghia's points.
    *explain:* "the page version of E7's convergence".
33. `nb.md` — **What would change if…** "…the quantity oscillated in time (the block at Re = 100)? Refine both Δx and Δt and compare
    time-averages and the Strouhal number — the checklist of §10.6."
34. `nb.pointer` — **N92 [C]** "A fourth computation — a cylinder at Re = 1000 by MacCormack, with a smoke line carried by an extra
    convection–diffusion equation (the 2-D version of (10.1)) — is shown in the book but not reproduced here (far beyond our run-time
    budget). > ⚠️ **The book calls the cylinder 'the fourth example'; §10.5 has three examples** — the cylinder by finite elements is the
    third (slip R8)."

#### C13 (continued) — the cylinder in a channel by mixed finite elements
35. `nb.md` — "**C13 (continued).** The third example uses C13's P2–P1 elements on a curved mesh around a cylinder, with Newton's method
    at every time step. The machinery is stated (not derived), coded and tested; our steady runs are small and cached, and the unsteady
    Re = 100 run was done once by `scripts/ch10_cylinder_fem.py` (our force history is in `reference/ch10/`)."
36. `nb.recap("R12", "Cylinder wake regimes", "Ch. 9 §9.8: at Re ≈ 1 the flow is nearly symmetric front-to-back (creeping flow, Ch. 8);
    by Re ≈ 10–40 two steady eddies sit behind the cylinder and grow with Re; above a threshold the wake sheds a Kármán street
    (`core.bluff_body.cylinder_flow_regime`). Our steady FE runs at Re = 1, 10, 40 show the same trend — **confined** (channel width 5d,
    sliding walls), so the numbers are not those of an unbounded cylinder.", where="Ch. 9 §9.8, Ch. 8 §8.6")` + `nb.figure` — streamlines and
    vorticity at Re = 1, 10, 40 from `FEM2.cylinder_steady(Re, mesh_level="coarse")` (cached; FAST: Re = 40 only). *see/read/change:*
    "a symmetric pattern at Re = 1, a growing attached eddy pair at 40"; "vorticity is made at the front and carried downstream"; "…a
    mesh that is not up–down symmetric: the steady solution is slightly asymmetric (as the book notes for its own mesh)".
37. `nb.recap("R13", "Strouhal number with the cyclic frequency", "St = f_s d/U = 1/τ̄ with the shedding frequency f_s in cycles per unit
    time and the period τ̄ in units of d/U (Ch. 4's $St=\\frac{\\Omega d}U$ (4.102) with Ω the cyclic frequency f_s; Ch. 9's `shedding_frequency`). > ⚠️ **The book compares its
    computed (confined) St with an unbounded value; the configurations differ** (a channel only a few diameters wide with sliding walls vs an
    open stream), so the comparison is qualitative only (slip R9). For an unbounded cylinder at Re = 100 the literature gives
    0.164–0.165 (Williamson 1996, experiments; Kravchenko et al.; Posdziech & Grundmann) — a different configuration.", where="Ch. 4
    §4.11, Ch. 9 §9.8")` + `nb.code`: read `reference/ch10/cylinder_fe_re100_forces.csv` (ours), `tau = ...;
    ch10.strouhal_from_period(tau)`, `ch10.dominant_frequency(t, CL)`. *expect:* the two estimates agree to 0.5 %; the value is printed
    with the label "confined, ours".
38. `nb.note` — **N93 [B], N94 [B]** "**Domain and mesh (Fig. 10.15–10.16 analogues, ours):** a channel of width W = 5d with the cylinder
    8d from the inflow Γ₁ and 16d from the outflow Γ₂, sliding walls Γ₃, Γ₄, the cylinder Γ₅; a triangular mesh refined near the body.
    Node counts: P2 velocity nodes = vertices + edges, P1 pressure nodes = vertices; Euler's check for a domain with one hole,
    V − E + T = 0." + `nb.code`: `mesh = FEM2.cylinder_channel_mesh()`; `print(FEM2.euler_check(mesh), FEM2.p2_node_counts(...))`.
    *expect:* 0 and the counts (printed). + `nb.figure` — the mesh (`plt.triplot`) with the boundaries coloured.
39. `nb.note` — **N95 [B], N96 [C]** "**Cartesian weak equations (10.156)–(10.159):** expanding $\\mathbf D:\\tilde{\\mathbf D}=u_x\\tilde u_x+
    \\frac12(u_y+v_x)(\\tilde u_y+\\tilde v_x)+v_y\\tilde v_y$ and separating ũ and ṽ gives $\\int_\\Omega(u_t+uu_x+vu_y)\\tilde u\\,d\\Omega-
    \\int_\\Omega p\\tilde u_x\\,d\\Omega+\\frac1{Re}\\int_\\Omega[2u_x\\tilde u_x+(u_y+v_x)\\tilde u_y]\\,d\\Omega=0$ (10.157), its v-twin (10.158)
    and $-\\int_\\Omega(u_x+v_y)\\tilde p\\,d\\Omega=0$ (10.159) — the minus sign is chosen so that the pressure blocks come out as B and Bᵀ.
    The Galerkin statements (10.160)–(10.162) are the same with every field replaced by its finite-element version on Ωʰ." equation
    (10.157) + `nb.code`: `ch10.weak_ns_components_sympy()["difference"]` (0).
40. `nb.note` — **N97 [B]** "**Time derivative at $t_{n+1}$ (10.163):** $\\frac{\\partial\\mathbf u}{\\partial t}(\\mathbf x,t_{n+1})\\approx
    \\alpha_t\\frac{\\mathbf u(\\mathbf x,t_{n+1})-\\mathbf u(\\mathbf x,t_n)}{\\Delta t}-\\beta_t\\frac{\\partial\\mathbf u}{\\partial t}(\\mathbf
    x,t_n)$ — $\\alpha_t=1,\\ \\beta_t=0$ is backward Euler (first order); $\\alpha_t=2,\\ \\beta_t=1$ is the trapezoidal rule (second order,
    'a variation of Crank–Nicolson'). (The book writes α, β; the subscript t keeps them apart from the FTCS numbers.)" equation
    (10.163) + `nb.code`: `FD.ode_scheme_order("backward"), FD.ode_scheme_order("trapezoidal")`. *expect:* ≈ 1, ≈ 2.
41. `nb.note` — **N98 [B], N99 [B], N100 [C], N101 [C], N102 [B], N103 [C]** "**Newton at every time step (10.164)–(10.183):** write
    each field as guess + correction, $\\mathbf u^h=\\mathbf u^*+\\mathbf u'$, $p^h=p^*+p'$ (10.164) (Newton, Ch. 6 P152; here u′, p′ are
    corrections, not perturbations); substitute, drop the product of two corrections $\\mathbf u'\\cdot\\nabla\\mathbf u'$ — the one idea of
    (10.165)–(10.167) — and the right-hand sides are the residuals of the current guess (zero at convergence), the left-hand sides the
    Jacobian. > ⚠️ **(10.166) prints** $\\beta\\,\\partial v^*/\\partial t(t_n)$; **it is** $\\beta\\,\\partial v/\\partial t(t_n)$ — known data, not
    an iterate (slip R10). Expanding in velocity shapes $N^u_A$ and pressure shapes $N^p_B$ (10.168)–(10.169) gives one equation per node
    (10.170)–(10.172); > ⚠️ **(10.172) prints the second sum with $u_{A'}$; it must be $v_{A'}$** (it comes from ∂v′/∂y; slip R6 — the
    printed version fails our Poiseuille test). Grouped by unknown they form the block system $\\begin{pmatrix}\\mathbf A_{uu}&\\mathbf
    A_{uv}&\\mathbf B_{up}\\\\\\mathbf A_{vu}&\\mathbf A_{vv}&\\mathbf B_{vp}\\\\\\mathbf B^T_{up}&\\mathbf B^T_{vp}&0\\end{pmatrix}\\begin{pmatrix}
    \\mathbf u\\\\\\mathbf v\\\\\\mathbf p\\end{pmatrix}=\\begin{pmatrix}\\mathbf f_u\\\\\\mathbf f_v\\\\\\mathbf f_p\\end{pmatrix}$ (10.173), with the
    entries (10.174)–(10.183) coded in `FEM2.assemble_newton_system`." equations (10.164), (10.173) + `nb.code`: the residual history of
    a steady Re = 40 solve (`FEM2.cylinder_steady(40.0)["residuals"]`) and the ratios r_{k+1}/r_k². *expect:* roughly constant ratios —
    quadratic convergence (5–7 iterations to 1e-10).
42. `nb.note` — **N104 [B], N105 [B], N106 [B], N107 [C], N108 [B], N109 [B]** "**One curved element (10.184)–(10.198):** the six-node
    isoparametric map $x(\\xi,\\eta)=\\sum_{a=1}^6x^e_a\\phi_a(\\xi,\\eta)$, $y=\\sum_ay^e_a\\phi_a$ (10.184) bends the element to the circle; the
    quadratic shapes $\\phi_1=\\zeta(2\\zeta-1),\\ \\phi_2=\\xi(2\\xi-1),\\ \\phi_3=\\eta(2\\eta-1),\\ \\phi_4=4\\xi\\zeta,\\ \\phi_5=4\\xi\\eta,\\
    \\phi_6=4\\eta\\zeta$ with ζ = 1 − ξ − η (10.185) (barycentric coordinates, P248); element expansions $u'=\\sum_{a=1}^6u^e_a\\phi_a$,
    $v'=\\sum_av^e_a\\phi_a$, $p'=\\sum_{b=1}^3p^e_b\\psi_b$ (10.186) with $\\psi_1=\\zeta,\\ \\psi_2=\\xi,\\ \\psi_3=\\eta$ (10.187) — > ⚠️ **the book
    prints the third expansion as '$v'=\\sum p^e_b\\psi_b$'; it is $p'$** (slip R7); the element matrices (10.188)–(10.197) are (10.175)–
    (10.183) with $N^u\\to\\phi$, $N^p\\to\\psi$; integrals by the 7-point rule on the parent triangle,
    $\\int_{\\Omega^e}f\\,d\\Omega=\\frac12\\sum_{l=1}^{7}f(\\xi_l,\\eta_l)J(\\xi_l,\\eta_l)W_l$ with $J=x_\\xi y_\\eta-x_\\eta y_\\xi$ (10.198)
    (P249; exact for polynomials of degree ≤ 5, no longer exact on curved elements); assembled by the node map, with essential
    conditions imposed by row replacement or a large penalty (both give the same Poiseuille solution)." equations (10.184), (10.185),
    (10.198) + `nb.figure` — one curved element: the parent triangle, its image with the mid-node on the circle (`FEM2.iso_map`), the 7
    quadrature points mapped. *see/read/change:* "a straight triangle becomes one with a curved side"; "J varies inside the curved
    element, so the rule is approximate there"; "…straight elements: the circle becomes a polygon — the error that curved elements remove".
43. `nb.note` — **N110 [B]** "**Re = 100, unsteady:** started from rest the wake is symmetric at first, then a small asymmetry grows and
    the wake settles into periodic shedding — a *supercritical Hopf bifurcation* (a steady state losing stability to an oscillation;
    Ch. 11) that produces a Kármán street. Our force history (drag, lift and torque against t̄ = tU/d, confined, ours) is plotted from the
    cached file; the torque's non-zero mean in the book comes from its mesh's asymmetry." + `nb.figure` — C_D(t̄), C_L(t̄), torque(t̄)
    from `reference/ch10/cylinder_fe_re100_forces.csv` (three stacked panels, blue/orange/muted) with the period τ̄ marked. *see:* "lift
    oscillating about zero, drag at twice the frequency"; *read:* "drag peaks twice per cycle because each shed vortex pulls back";
    *change:* "…Re = 40: all three curves flatten to constants — below the Hopf point".

---

### A.6 §10.6 Concluding Remarks — C15 (continued): N111 N112 · S01 S02 · summary
1. `nb.section("10.6", "Concluding Remarks", intro="**What is this section about?** CFD is a tool, not an oracle. The chapter ends
   with a checklist for trusting a result and a glance at methods not covered. We apply the checklist to our own cavity solver.")`
2. `nb.note` — **N111 [B]** "**C15 (continued) — the verification checklist**, applied to our MAC cavity (a five-row table with our
   results): | check | what we did | outcome |: benchmark with a known answer (Ghia 1982 centreline, 64²: within 2 %) · mesh refinement
   (16 → 32 → 64, observed order printed by C15) · time-step refinement (Δt and Δt/2 at 32²: steady result unchanged to 1e-6) ·
   insensitivity to boundary treatment (the ghost-value vs one-sided tangential no-slip options: centreline differs < 0.5 %) · parameter
   sensitivity (Re = 100 vs 400 from the cache: the eddy centre moves as Ghia reports)." + `nb.code`: the time-step refinement run
   (`MAC.cavity(Re=100, n=32, dt=…)` twice) and the difference printed. *expect:* < 1e-5 in u.
3. `nb.note` — **N112 [C]** "**Other methods, named:** spectral and spectral-element methods (global smooth basis functions; our Ch. 7
   `kdv_solve` is pseudo-spectral); lattice-gas and lattice-Boltzmann methods and dissipative particle dynamics — kinetic models that do
   not start from the continuum Navier–Stokes equations. (Hou et al.'s cavity centres used in C14 are lattice-Boltzmann results.)"
4. `nb.pointer` — **S01** "Exercises 10.1–10.6 are in the book; our own versions of their ideas are C04 (Noye's region), N17 (an implicit
   scheme's stability), N113 (the heated rod), C08 (element matrices), C09 (central vs upwind) and C14 (the cavity). No exercise text,
   input or answer is reproduced here."
5. `nb.pointer` — **S02** "Literature and supplemental reading: see the book's list; the data we use are cited where they are used —
   Ghia, Ghia & Shin (1982), J. Comput. Phys. 48, 387–411, and Hou, Zou, Chen, Doolen & Cogley (1995), J. Comput. Phys. 118, 329–347."
6. `nb.summary(clicked=[…15…], feeds_forward=[…], left_out=[…])`. **Clicked (one per CORE):** C01 "A stencil is a weighted sum of
   Taylor series; the first term it fails to cancel is its error, and the power of Δx there is its order — a slope on log–log axes." ·
   C02 "FTCS computes each new value as a weighted average of three old ones, with α = uΔt/(2Δx) and β = DΔt/Δx²." · C03 "Put the exact
   solution into the scheme: what is left over, E = Δt/2 T_tt + uΔx²/6 T_xxx − DΔx²/12 T_xxxx, must vanish as the grid shrinks —
   consistency." · C04 "Every Fourier mode of the round-off is multiplied by G(θ) each step; FTCS is stable only inside 0 ≤ 4α² ≤ 2β ≤ 1,
   BTCS always." · C05 "Upwinding is stable exactly when the characteristic's foot stays inside the stencil, C = uΔt/Δx ≤ 1; with
   consistency, that stability is what buys convergence (Lax)." · C06 "Multiply by a test function and integrate by parts: only first
   derivatives remain, the Neumann condition comes free (natural) and the Dirichlet one is built into the space (essential)." · C07
   "Galerkin with hats turns the weak form into M ḋ + K d = F, whose rows are FD stencils with a ⅙–⅔–⅙ mass." · C08 "Compute a 2 × 2
   block on each element and scatter-add it — the same loop assembles any mesh." · C09 "The centred discrete solution is rʲ with r < 0
   once R_cell > 2 — wiggles; upwinding removes them by adding a diffusivity 0.5 R_cell D." · C10 "A slightly compressible fluid gives
   the pressure its own equation; MacCormack's forward-then-backward stages average to a second-order (Lax–Wendroff) step, stable for
   C ≤ 1 on the sound speed." · C11 "Split the step: move the flow ignoring pressure, then add the one gradient that removes the
   divergence — a Poisson solve and a curl-free correction." · C12 "Pressure at centres and velocities on faces make divergence and
   gradient two-point differences: no pressure boundary condition, no hidden checkerboard — the C-grid of climate models." · C13
   "Velocity and pressure spaces must satisfy the inf–sup condition; P2–P1 does, P1–P1 does not." · C14 "Reproducing a published
   benchmark within a stated tolerance is evidence that the code solves the intended equations." · C15 "Three grids give the observed
   order (a code check) and a Richardson estimate of the grid-independent answer (with an error bar)." **Feeds forward:** Ch. 11 (FD
   operators, eigenvalue stability of discretised base flows, the Hopf bifurcation of N110), Ch. 12 (numerical vs eddy diffusivity,
   resolution), **Ch. 13** (the C-grid, CFL for gravity waves √(gH), projection = Boussinesq pressure solve, upwind tracer advection's
   implicit diffusion), Ch. 15 (MacCormack/Lax–Wendroff for shocks). **Left out:** the Re = 1000 smoke-line cylinder (N92), the full
   derivations of the Newton and quadrature machinery (§4c of the curation), turbulence models (Ch. 12).
7. `nb.save()`.

### A.7 Placement check (every curation id has exactly one home)
- **C:** C01 A.2 5 · C02 A.2 20 · C03 A.2 33 · C04 A.2 46 · C05 A.2 70 · C06 A.3 2 · C07 A.3 17 · C08 A.3 38 · C09 A.4 4 · C10 A.4 31 ·
  C11 A.4 52 · C12 A.4 70 · C13 A.4 96 (+ continued A.5 35–43) · C14 A.5 2 · C15 A.5 20 (+ continued A.6 2) — 15.
- **R:** R01 R02 R03 A.2 2–4 · R04 R05 A.2 44–45 · R06 R07 A.4 2–3 · R09 A.4 30 · R08 A.4 51 · R10 A.4 95 · R11 A.5 19 · R12 R13 A.5
  36–37 — 13.
- **S:** S01 S02 A.6 4–5.
- **N:** §10.1 N01–N04 (A.1) · C01 N05 N07 · C02 N06 N08 N09 · C03 N10 · C04 N11 N12 N13 N14 N15 N17 · C05 N16 N18 N19 N113 · C06 N20
  N21 N22 · C07 N23–N33 · C08 N34–N39 · C09 N40 N42–N50 · C10 N51–N59 · C11 N41 N60–N63 N72 N74 · C12 N64–N71 N73 · C13 N75 N76 (core)
  and N93–N110 (continued) · C14 N77–N85 · C15 N86–N91 + N111 N112 (continued) · N92 pointer (A.5 34) — 113.
- **D:** D01 C01 · D02 C02 · D03 C03 · D04–D07 C04 · D08 C05 · D09 D10 C06 · D11 D12 C07 · D13 D14 C08 · D15–D17 C09 · D18 C10 · D19 C11
  · D20 D21 C12 · D22 C13 · D23 C15 — 23, each inside its curation CORE block.
- **Explainers:** E1 C01 · E2 C04 · E3 C05 · E5 C08 · E4 C09 · E6 C12 · E8 C13 · E7 C14 — 8, each embedded once. B1 not embedded.
- **Animations:** A1 C05 · A2 C10 · A3 C11 · A4 C14. **Plotly:** F1 C01 · F2 C04 · F3 C05 · F4 C07 · F5 C09 · F6 C15 · F7 C13. **Live:**
  C04, C09.
- **Primers P221–P253 (33):** P221 P222 C01 · P223 C02 · P224 C03 · P225 P226 C04 · P227 P228 P229 C05 · P230 P231 C06 · P232 P233 C07 ·
  P234 P235 C08 · P236 C09 · P237 P238 P239 C10 · P240 P241 C11 · P242 P243 P244 P245 C12 · P246 P247 P248 P249 P250 C13 · P251 P252 P253
  C14. (The streamfunction of a computed field is a reminder of Ch. 4 ψ + Ch. 6 `solve_poisson`, not a new primer.)
- **Visual per CORE block:** C01 fig + F1 + E1 · C02 fig · C03 fig · C04 figs + F2 + live + E2 · C05 fig + F3 + A1 + E3 · C06 fig · C07
  figs + F4 · C08 figs + E5 · C09 fig + F5 + live + E4 · C10 fig + A2 · C11 fig + A3 · C12 figs + E6 · C13 figs + F7 + E8 · C14 figs + A4
  + E7 · C15 fig + F6.
- **Code per CORE block:** every block has ≥ 2 `nb.code` cells and one `check_agree` (from-scratch or parity).

---

## Part B — explainer storyboards

Common to all nine: created with `tools/new_viz.py`; `<meta name="viz:chapter" content="ch10">`; tabs Walkthrough · Explore ·
Explain · Derivation · Equations · Code · Check; every displayed number is computed by a JS function that mirrors a `ch10` callable
and is proved by `selftest()` parity rows (`py:` expressions use only `ch10.…`, `np.…`, numbers, strings, lists and keywords,
indexed down to one number — Part C convention 2 and C.10). Explain is "Explanation & interpretation" in numbered sections built
with `Viz.work.step / line / box / table / hint / interpret`, modelled on `forced_damped_vibrations.html` (the regime-dependent
reading) and `amplitude_phase_second_order_II_3.html` (a numbered derivation with live numbers, one section per optional view):
**0** what the views show and what each colour means (two sentences on phones) · **1…n** every displayed quantity from the controls
("formula = substituted = result — why", results boxed) · a section per view hidden on phones, or a hint · the values at the
current time (live) · **Reading the current setting** (regime-dependent). Derivation steps are copied from Part F (same `did`
titles, same step count, same order; phones shorten *why* to its first sentence; plain-text *why* and *watch* never contain raw
TeX). Every tour, Explain, Derivation, notes, status, equation and quiz text that names a book equation **writes it out** next to
its number (`tools/eq_refs.py` → 0; no bare numbers in `<meta>` strings). Drafts below write equations in Unicode for readability;
builders set each in TeX (backslashes doubled in JS strings). Colours as in convention 5. Walkthrough texts ≤ 45 words, step 1
≤ 24 words, ≤ 2 extras per step, `play: false` on steps that quote numbers, `autoplay: false`. A view hidden on portrait phones
never carries a step's key number (repeat it in a visible title or readout). Symbols as in convention 3 (upright $\mathrm i$; grid
index j in our own lines; $K_e$, $\mathcal S$, $Ma$). **Computation budget in the browser:** 1-D grids ≤ 400 nodes; MAC ≤ 32² live
with a banded LU of the pinned Poisson matrix factorised once per grid; dense FE ≤ 400 unknowns; anything larger comes from an
embedded table labelled "ours, computed by fluidpy" with a parity row at a small size. Data embedded in JS: `CAVITY` (E7:
centreline u(y) and ψ on a 33 × 33 sampling for MAC 64² and 128² at Re = 100 and 400, MacCormack 64² at Re = 100, Ghia's Table I
columns and centres, Hou's centres — ≤ 4 significant figures, ≈ 40 kB; from `scripts/ch10_cavity.py`), `INFSUP` (E8: β_h and
counts for P1P1/P2P1, n = 2…16, from `reference/ch10/infsup_table.csv`). Parallel builders use private scratch subfolders
(`<scratchpad>/<slug>/`). **Library note:** E5, E6, E8 need small dense/banded solvers: E5 a tridiagonal (Thomas) solve; E6 a banded
LU (bandwidth n) — local helpers, or `Viz.num.lu` if the orchestrator adds it (viz_patterns item 11); E8 a dense LU and a symmetric
generalised eigen-solver (Cholesky of M_p + cyclic Jacobi, ≤ 49 × 49; ~60 lines, local).

### E1 · fd_stencil_order
- **Title:** "What does 'second-order accurate' really mean?" · **Summary:** "Shrink h and watch the stencil's error fall along a
  slope of 1 or 2 — the first Taylor term it fails to cancel — until round-off takes over; in FTCS mode the three truncation terms of
  (10.17) appear as bars." · **CORE:** C01, C03 (also R02, R03, N07, N10, N80) · **Reference:** `fid_formula_lab.html` (terms as bars
  that add to a total, story + code with the active line) and `amplitude_phase_second_order_II_3.html` (numbered live Explain).
- **meta:** `viz:order 1` · `viz:sections 10.2` · `viz:equations 10.4 10.5 10.6 10.7 10.8 10.16 10.17` · `viz:fluidpy ch10.fd_weights
  ch10.stencil_taylor_coefficients ch10.stencil_error ch10.test_function ch10.truncation_terms` · `viz:derivations D01 D03`.
- **Physics (JS ↔ Python):** `weights(offsets, m)` — Gauss elimination of the Taylor-matching system in doubles ↔
  `ch10.fd_weights(offsets, m)`; `taylorCoeffs(kind, m)` → {coeffs, order, leading} ↔ `ch10.stencil_taylor_coefficients(kind, m)`;
  `fk(name, x, k)` the exact k-th derivative of sin, exp, e^{−x²} ↔ `ch10.test_function(name, x, k)`; `stencilErr(name, x0, h, kind, m)`
  ↔ `ch10.stencil_error`; `roundoff(name, x0, h, kind, m)` = ε·Σ∣w∣·∣f(x₀)∣/h^m (an estimate, labelled so); `truncTerms(u, D, dx, dt, x,
  t)` with the closed-form derivatives of the travelling Gaussian ↔ `ch10.truncation_terms`.
- **Views** (rows [1.15, 1]): 1. `stencil` "The stencil at x₀" (row 0, flex 1.3) — f(x) on x₀ ± 0.6 (muted), the exact tangent (or the
  exact parabola for m = 2) as a muted dashed ghost, the stencil's secant/parabola through its nodes (orange for centred, teal for
  one-sided), nodes as dots labelled with their weights ("−½", "0", "+½"); when h < 0.02 the nodes are drawn in a magnified inset
  ("×1/h zoom") so they stay visible; title "h = 0.1 · error −9.0e-4". 2. `error` "∣error∣ against h" (row 0, flex 1) — log–log, h ∈
  [1e-12, 0.5], ∣error∣ ∈ [1e-17, 1]; the chosen stencil's full curve (bold; faint beyond the current h — the 'so far' idiom), the other
  four stencils as faint ghosts, slope-1 and slope-2 guides (muted, labelled), the round-off band ε/h (rose, shaded), the current dot,
  the best-h marker (◆); pointer: click → inspector at that h. 3. `terms` "Taylor terms of the stencil" (row 1, `hidePortrait: true`)
  — bars for k = 0…5 of $c_kh^{k-m}f^{(k)}(x_0)$: the kept derivative (blue), cancelled terms (muted, height 0, labelled "cancels"),
  the leading error term (rose), later terms (light rose); in mode **FTCS** the bars become the three terms of (10.17) — time
  $\frac{\Delta t}2T_{tt}$ (amber), dispersive $u\frac{\Delta x^2}6T_{xxx}$ (purple), diffusive $-D\frac{\Delta x^2}{12}T_{xxxx}$ (teal) —
  and their total next to the measured one-step residual (muted outline).
- **Controls (≤ 5 visible):** `stencil` chips forward / backward / central / central2 (second derivative) / onesided2 (default
  central) · `lh` "Step $\log_{10}h$" −12 … −0.3, step 0.05, default −1 (h = 0.1; the transport parameter) · `func` chips sin / exp /
  Gaussian (optional, default sin) · `x0` "Point $x_0$" −1 … 2, step 0.05, default 1 (optional) · `mode` chips stencil / FTCS
  truncation; in FTCS mode `dx` "$\Delta x$" 0.005 … 0.05 m (default 0.01) and `beta` "$\beta=D\Delta t/\Delta x^2$" 0.05 … 0.5 (default
  0.25) replace `lh` and `stencil` (u = 0.5 m/s, D = 0.01 m²/s fixed and printed).
- **Transport:** `lh` from −0.3 to −12, rate −1.2 decades/s, `end: 'hold'`, hold 2 s; end card (`Viz.card`): "best h ≈ 2e-6 for
  central: error 7e-13 — below this, round-off wins".
- **Presets:** "central, h = 0.1" {stencil: central, lh: −1} · "forward, h = 0.1" {stencil: forward, lh: −1} · "round-off: forward,
  h = 1e-10" {stencil: forward, lh: −10} · "one-sided at a wall" {stencil: onesided2, lh: −1.3} · "FTCS β = 1/6" {mode: ftcs, beta:
  0.1667, func: sin}.
- **Status:** "📉 truncation-dominated: slope 2 — halving h divides the error by 4.00" (truncation estimate > 10 × round-off) ·
  "⚖️ near the best h: truncation ≈ round-off" · "🎲 round-off dominated: error ≈ ε/h = 2.2e-8 — a smaller h makes it worse" · FTCS mode:
  "✨ β = 1/6: the time and diffusion terms cancel (fourth order for pure diffusion)" (∣β − 1/6∣ < 0.005) else "E ≈ 3 terms, first order in
  Δt, second in Δx".
- **Readouts:** "Error" (signed) · "Leading term" · "Error ratio h→h/2" · "Order p" · "Best h".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **terms** (Taylor or truncation bars) · **presets** ·
  **status** · **transport** · **inspector** (click on the error curve: "h = 0.05: stencil = (sin 1.05 − sin 0.95)/0.1 = 0.540077;
  exact cos 1 = 0.540302; error = −2.251e-4; leading term h²/6·(−cos 1) = −2.251e-4") · **modes** (stencil / FTCS truncation).
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**The stencil**: grey is f, dashed grey its exact tangent at x₀, orange the straight line (or parabola)
     through the stencil's points — its slope is the stencil's answer. **∣error∣ against h**: the bold curve is your stencil's error for
     every h, the dot is your h; the rose band is where round-off lives. **Taylor terms** (hidden on phones): what the weighted sum of
     Taylor series contains — blue survives as the derivative, grey cancels, rose is the error."
  1. *The weights of your stencil* — "offsets **s** = (−1, 0, 1), derivative m = **1**: the Taylor-matching conditions
     $\sum_jw_js_j^k/k!=\delta_{km}$ for k = 0, 1, 2 give $w=(-\tfrac12,0,\tfrac12)$" boxed; why: the stencil must return f′ exactly for
     1, x, x² (and fail first at x³).
  2. *Which Taylor terms cancel* — "c₀ = Σw = **0** (f cancels), c₁ = Σw s = **1** (the derivative survives), c₂ = Σw s²/2 = **0**
     (cancels), c₃ = Σw s³/6 = **1/6** ← first survivor ⇒ error ≈ (1/6)h²f‴(x₀): order p = **2**" (table, current stencil's row lit);
     why: (10.4)–(10.5) subtracted.
  3. *The error with your h* — "stencil = **0.539402**; exact f′(1) = cos 1 = **0.540302**; error = **−9.001e-4**; the leading term
     predicts (1/6)(0.1)²(−cos 1) = **−9.005e-4** (they differ by the next term, O(h⁴)); halving h: error ratio = **4.00** = 2^p" boxed.
  4. *Round-off* — "each value of f carries a relative error ε = 2.2e-16; the stencil multiplies them by Σ∣w∣/h^m = **10**, so round-off
     ≈ **2e-15**; truncation ≈ **9e-4**: truncation wins by 10¹¹. The two are equal near h* ≈ (ε/c)^{1/(p+m)} = **2.2e-6**."
  5. *FTCS mode* (or the hint "switch the mode to FTCS truncation to see (10.17)") — "u = 0.5 m/s, D = 0.01 m²/s, Δx = **0.01** m, β =
     **0.25** ⇒ Δt = βΔx²/D = **0.0025** s. At the pulse centre: time term (Δt/2)T_tt = **…**, dispersive uΔx²/6 T_xxx = **…**, diffusive
     −DΔx²/12 T_xxxx = **…**; E = **…** vs the measured one-step residual **…**. For pure diffusion E ∝ Δx²(β/2 − 1/12): zero at β = 1/6."
  6. *Reading the current setting* — truncation: "You are on the straight part: every halving of h cuts the error by 2^p = **4**. This is
     what 'second order' means." · near-best: "You are at the bottom of the V: the best this stencil can do in double precision." ·
     round-off: "Round-off rules: the difference of two nearly equal numbers keeps only ~16 digits, divided by a tiny h. Use a larger h
     — or a higher-order stencil." · one-sided: "One-sided second order costs three points instead of two but keeps p = 2 at a wall —
     the cavity's wall densities use exactly this (N80)."
- **Derivation tab:** **D01** (9 steps) `view: 'stencil'`, goal `set` {stencil: 'central', lh: −1}; step 3 (forward difference) `set`
  {stencil: 'forward'} with `live` "error −0.0429, leading term −0.0421"; step 6 (centred first derivative) `set` {stencil: 'central'}
  `watch` "the orange line is parallel to the dashed tangent"; step 9 (second derivative) `set` {stencil: 'central2'}. **D03**
  (11 steps) `view: 'terms'`, goal `set` {mode: 'ftcs', beta: 0.25}; step 10 `live` "E = Δt/2 T_tt + uΔx²/6 T_xxx − DΔx²/12 T_xxxx =
  … with your Δx"; step 11 `watch` "set β = 1/6: two bars cancel". Interpret: `s => "With h = " + h + " your stencil's error is " +
  err + "; it is " + (trunc ? "truncation" : "round-off") + "-dominated."`.
- **Code:**
  ```python
  w = ch10.fd_weights({{offs}}, {{m}})           # Taylor matching: {{wtxt}}
  co = ch10.stencil_taylor_coefficients("{{kind}}", {{m}})
  print(co["order"], co["leading"])              # p = {{p}}, first survivor {{lead}}
  h = {{h}}                                       # step [same units as x]
  err = ch10.stencil_error("{{func}}", {{x0}}, h, "{{kind}}", {{m}})
  print(err)                                      # {{err}}  (≈ lead·h^p·f^(m+p))
  ```
- **Walkthrough (6 steps):** 1. "Slope from three numbers" — "A computer knows f only at points. How good is a slope from neighbours?
  Press ▶: h shrinks." `play: true` · 2. "Taylor, twice" — "Expand f at x₀ ± h (10.4)–(10.5). Subtracting kills the f″ terms; the first
  survivor is (h²/6)f‴." `terms: true`, `derive: {id: 'D01', step: 6}` · 3. "Order = slope" — "Halve h twice: the central error drops
  by 4 each time, the forward one by 2. The slope on log axes is the order." `set` {stencil: 'central', lh: −1}, `readouts: ['ratio',
  'order']` · 4. "Round-off floor" — "Below h ≈ 10⁻⁶ the error climbs again: differences of nearly equal numbers lose digits (ε/h)."
  `set` {lh: −9}, `highlight: ['view:error']` · 5. "The scheme's leftover" — "FTCS mode: put the exact solution into the scheme; three
  terms remain (10.17). They set the scheme's order." `set` {mode: 'ftcs'}, `derive: {id: 'D03', step: 10}` · 6. "Your turn" — "Predict:
  at β = 1/6, which two bars cancel? Then drag β and check." `controls: ['beta', 'dx']`.
- **Equations:** `taylor` "Taylor at the neighbours" ref 'Eq. (10.4)–(10.5)' $T^n_{i\pm1}=T^n_i\pm\Delta x\,T_x+\frac{\Delta x^2}2T_{xx}\pm
  \frac{\Delta x^3}6T_{xxx}+\frac{\Delta x^4}{24}T_{xxxx}+O(\Delta x^5)$ · `first` "First-derivative stencils" ref 'Eq. (10.6)'
  $\frac{T_{i+1}-T_i}{\Delta x}+O(\Delta x)$, $\frac{T_i-T_{i-1}}{\Delta x}+O(\Delta x)$, $\frac{T_{i+1}-T_{i-1}}{2\Delta x}+O(\Delta x^2)$ live
  "error = …" · `second` "Second derivative" ref 'Eq. (10.7)' $\frac{T_{i+1}-2T_i+T_{i-1}}{\Delta x^2}+O(\Delta x^2)$ · `time` "In time"
  ref 'Eq. (10.8)' $\frac{T^{n+1}-T^{n-1}}{2\Delta t}+O(\Delta t^2)$ · `trunc` "Truncation error of FTCS" ref 'Eq. (10.17)'
  $E^n_i=\frac{\Delta t}2T_{tt}+u\frac{\Delta x^2}6T_{xxx}-D\frac{\Delta x^2}{12}T_{xxxx}+O(\Delta t^2,\Delta x^4)$ live "E = …" · symbols h/Δx
  (m), Δt (s), u (m/s), D (m²/s), p (–), ε (–).
- **Check yourself:** (1) "Halve h for the forward stencil. By what factor does the error fall?" — "About 2: first order, the leading
  term is (h/2)f″." `set {stencil: 'forward', lh: −1}` · (2) "Why is the central stencil's error zero for f = x² but not for x³?" — "Its
  weights cancel the f″ term exactly; the first survivor involves f‴, which is 0 for x² and 6 for x³." · (3) "Where is the best h for
  the forward stencil, roughly?" — "Near 10⁻⁸ (√ε): truncation (h/2)∣f″∣ equals round-off ε∣f∣/h there." `set {stencil: 'forward', lh:
  −8}` · (4) "In FTCS mode, what β makes the scheme more accurate for pure diffusion?" — "β = 1/6: the Δt/2·T_tt term cancels
  −DΔx²/12·T_xxxx." `set {mode: 'ftcs'}`.
- **Selftest parity rows:** `{name: 'central error h=0.1', js: stencilErr('sin', 1, 0.1, 'central', 1), py: 'ch10.stencil_error("sin",
  1.0, 0.1, "central", 1)', rtol: 1e-10}` · `{name: 'forward error h=0.05', js: stencilErr('sin', 1, 0.05, 'forward', 1), py:
  'ch10.stencil_error("sin", 1.0, 0.05, "forward", 1)', rtol: 1e-10}` · `{name: 'onesided weight w0', js: weights([0, 1, 2], 1)[0], py:
  'ch10.fd_weights([0, 1, 2], 1)[0]', rtol: 1e-12}` · `{name: 'central leading', js: taylorCoeffs('central', 1).leading, py:
  'ch10.stencil_taylor_coefficients("central", 1)["leading"]', rtol: 1e-12}` · `{name: 'central2 weight', js: weights([-1, 0, 1], 2)[1],
  py: 'ch10.fd_weights([-1, 0, 1], 2)[1]', rtol: 1e-12}` · `{name: 'FTCS E total', js: truncTerms(0.5, 0.01, 0.01, 0.0025, 0.35,
  0.1).total, py: 'ch10.truncation_terms(0.5, 0.01, 0.01, 0.0025, 0.35, 0.1)["total"]', rtol: 1e-9}` · invariant `{name: 'weights sum
  0', js: weights([-1, 0, 1], 1).reduce((a, b) => a + b), expect: 0, atol: 1e-15}`.
- **Fit plan:** 360×640: status (1 line), `stencil` (45 %) over `error` (55 %), `terms` hidden (the leading term and order in the
  `error` title); chips wrap; walkthrough card paged. 844×390: `stencil` | `error`, row 2 hidden. 1000×700 and desktop: rows
  [1.15, 1] with `terms` under both.

### E2 · von_neumann_amplification
- **Title:** "Why does β = 0.51 blow up when 0.50 does not?" · **Summary:** "Every Fourier mode of the round-off is multiplied by the
  same complex G(θ) each step; drag (α, β) across Noye's region and watch G leave the unit circle and a 10⁻¹⁰ kick grow into a
  zigzag." · **CORE:** C04, C02 (also N08, N11–N15, N17, R04, R05) · **Reference:** `angular_frequency_explorer_1.html` (linked views
  on one clock, presets, "right now" table) and `forced_damped_vibrations.html` (regime status + interpretation panel).
- **meta:** `viz:order 2` · `viz:sections 10.2` · `viz:equations 10.10 10.11 10.13 10.19 10.24 10.25 10.26 10.27 10.28` ·
  `viz:fluidpy ch10.amplification_factor ch10.amplification_modulus ch10.max_amplification ch10.worst_theta ch10.ftcs_stable
  ch10.stability_verdict ch10.propagate_error` · `viz:derivations D04 D05 D06 D07`.
- **Physics:** `G(theta, a, b, scheme)` → {re, im} for "ftcs", "btcs", "upwind", "cn" (Part C 1.19 formulas; complex arithmetic by
  hand, `Viz.cx` if present) ↔ `ch10.amplification_factor` (`.real`/`.imag`) and `ch10.amplification_modulus`; `Gmax(scheme, a, b)` scan
  of 2001 angles ↔ `ch10.max_amplification`, `worstTheta` ↔ `ch10.worst_theta`; `noye(a, b)` ↔ `ch10.ftcs_stable`; `verdict(scheme, a, b)`
  exact reason strings ↔ `ch10.stability_verdict`; `march(xi, a, b, scheme)` one periodic step of (10.19) (BTCS/CN by a cyclic
  tridiagonal solve, N = 64) ↔ `ch10.propagate_error` (seeded initial kick from `Viz.rng(7)`, amplitude 1e-10; parity uses a
  deterministic kick $\xi_j=10^{-10}\cos(\pi j)+10^{-10}\sin(2\pi j/N)$).
- **Views** (rows [1, 1]): 1. `plane` "G in the complex plane" (row 0, flex 1, `equal: true`) — unit circle (muted), the curve G(θ) for
  θ ∈ [0, π] (orange FTCS, blue BTCS, teal upwind, purple CN), the part outside the circle drawn rose, the probe θ as a dot with its ∣G∣
  printed; the point G = 1 at θ = 0 marked. 2. `region` "Noye's region in the (α, β) plane" (row 0, flex 1) — α ∈ [0, 0.6], β ∈ [0,
  0.8]; the region 0 ≤ 4α² ≤ 2β ≤ 1 shaded teal with the two edges labelled "4α² = 2β (long waves)" and "2β = 1 (zigzag)"; the current
  point (draggable: pointerdown/move sets α, β); in BTCS mode the whole quadrant shaded blue. 3. `march` "The error, step by step"
  (row 1, `hidePortrait: false`, flex 1) — left half: ξ_j at the current step on 64 nodes (bars, symmetric log scale ±10^k with the
  current max labelled); right half: log₁₀ max∣ξ∣ against step n (0 … 800) with the 'so far' curve bold, the predicted slope log₁₀ Gmax
  as a dashed ghost; the dominant wavelength labelled ("2Δx zigzag" / "long wave"). On portrait phones the `march` right half only.
- **Controls (≤ 5 visible):** `alpha` "$\alpha=u\Delta t/(2\Delta x)$" 0 … 0.6, step 0.005, default 0.1 · `beta` "$\beta=D\Delta t/\Delta
  x^2$" 0 … 0.8, step 0.005, default 0.2 · `scheme` chips FTCS / BTCS / upwind / Crank–Nicolson (default FTCS) · `theta` "Probe
  $\theta=k\pi\Delta x$" 0 … π, step π/180, default π/2 (optional) · `n` step count (the transport).
- **Transport:** `n` 0 → 800 steps, rate 60 steps/s, `end: 'hold'`, hold 2 s; end card: "after 800 steps: max∣ξ∣ = 1.8e-4 (grew
  ×1.8e6) — unstable" or "… = 3e-40 — stable". Growth is capped at 10⁶ × the start (a rose "capped" label), never NaN.
- **Presets:** "β = 0.50 (edge, stable)" {alpha: 0, beta: 0.5} · "β = 0.51" {alpha: 0, beta: 0.51} · "pure convection" {alpha: 0.2,
  beta: 0} · "Noye edge 4α² = 2β" {alpha: 0.3, beta: 0.18} · "C02 example" {alpha: 0.1, beta: 0.2} · "BTCS, β = 100" {scheme: 'btcs', beta:
  100 (the β slider extends to 100 on a log scale in BTCS mode), alpha: 0.1}.
- **Status:** exact reason strings from `stability_verdict`: "✅ stable: max∣G∣ = 1.000 (θ = 0)" · "💥 unstable: 4a^2 > 2b, long waves
  grow — max∣G∣ = 1.02 at θ = 0.2π" · "💥 unstable: 2b > 1, the zigzag grows — max∣G∣ = 1.04 at θ = π" · "💥 unstable: pure convection,
  every wave grows" · "🛡️ stable for every dt (implicit)".
- **Readouts:** "∣G(θ)∣" · "max∣G∣" · "worst θ/π" · "4α² vs 2β" · "steps to ×10⁶".
- **Depth features:** Explain + Code + Derivation · **linked views** (3 on one clock) · **transport** · **presets** · **status** ·
  **inspector** (click on the plane curve or pick θ: "G = (α+β)e^{−iθ} + (1−2β) + (β−α)e^{iθ} = 0.3(−i) + 0.6 + 0.1(i) = 0.6 − 0.2i; ∣G∣² =
  0.40") · **modes** (FTCS / BTCS / upwind / CN).
- **Explain:**
  0. *What the views show* — "**Complex plane**: each point of the orange curve is the factor G one Fourier wave of the error is
     multiplied by per step; inside the grey circle it shrinks, outside (rose) it grows. **(α, β) plane**: the teal lens is where FTCS is
     stable; the dot is your setting. **The error**: a tiny random kick marched by your scheme — bars now, and its size over time."
  1. *Your numbers* — "α = uΔt/(2Δx) = **0.1**, β = DΔt/Δx² = **0.2** (10.11); Courant C = 2α = **0.2**."
  2. *G at your probe θ* — "(10.24): G = (α+β)e^{−iθ} + (1−2β) + (β−α)e^{iθ} = **0.3**(cos θ − i sin θ) + **0.6** + **0.1**(cos θ + i sin θ)
     = **0.6 − 0.2i**" boxed; "∣G∣² = (1 − 4β sin²(θ/2))² + (2α sin θ)² (10.26) = **0.36 + 0.04 = 0.40**, ∣G∣ = **0.632**: this wave keeps
     63 % per step."
  3. *The worst wave* — "scanning θ ∈ [0, π]: max∣G∣ = **1.000** at θ = **0**; Noye (10.27): 4α² = **0.04** ≤ 2β = **0.4** ≤ 1 ✓" (table of
     the two edges with the failing one lit).
  4. *What happens to the kick* — "after n = **n** steps the worst wave is multiplied by max∣G∣ⁿ = **…**: from 10⁻¹⁰ to **…**; to reach 1 it
     needs ln(10¹⁰)/ln(max∣G∣) = **…** steps" (live; "never" when max∣G∣ ≤ 1).
  5. *Implicit and other schemes* (mode-dependent) — BTCS: "G = 1/(1 + 4β sin²(θ/2) + 2iα sin θ): the denominator's modulus is at least
     its real part ≥ 1, so ∣G∣ ≤ 1 for every Δt (D07); at β = 100, the zigzag keeps 1/401 per step." Upwind: "G = 1 − C(1 − e^{−iθ}), stable
     iff C ≤ 1 (E3)." CN: "G = (1 − …)/(1 + …): ∣G∣ ≤ 1 always, and ∣G(π)∣ → 1 as β grows (slow-dying zigzags, Ch. 8)."
  6. *Reading the current setting* — stable: "Every wave shrinks or stays: round-off stays at round-off level." · long-wave unstable:
     "Convection kicks the long waves harder than diffusion can damp them: 4α² > 2β. Reduce Δt, or add diffusion." · zigzag unstable:
     "Diffusion overshoots on the shortest wave: 1 − 4β < −1 ⇒ the zigzag flips and grows by ∣1 − 4β∣ = **1.04** per step. Δt must satisfy
     Δt ≤ Δx²/(2D) (10.28)." · pure convection: "With D = 0, ∣G∣² = 1 + 4α² sin²θ > 1 for every wave: FTCS cannot do pure convection —
     upwind can (E3)."
- **Derivation tab:** **D04** (11 steps) `view: 'march'`; goal `set` {alpha: 0.1, beta: 0.2, scheme: 'ftcs'}; step 3 `watch` "the error
  obeys the same rule"; step 8 `live` "G = 0.6 − 0.2i at θ = π/2". **D05** (8 steps) `view: 'plane'`; step 8 `live` "∣G∣² = 0.36 + 0.04".
  **D06** (12 steps) `view: 'region'`; step 9 `set` {alpha: 0.3, beta: 0.18} `watch` "on the long-wave edge"; step 10 `set` {alpha: 0,
  beta: 0.5}; step 12 `set` {alpha: 0.2, beta: 0} `watch` "the whole curve is outside the circle". **D07** (7 steps) `view: 'plane'`,
  `set` {scheme: 'btcs'}. Interpret: `s => "With α = " + a + ", β = " + b + ": max∣G∣ = " + gmax + " — " + verdict.reason`.
- **Code:**
  ```python
  a, b = {{alpha}}, {{beta}}                       # alpha = u dt/(2 dx), beta = D dt/dx^2
  th = np.linspace(0, np.pi, 2001)                 # every wave the grid can carry
  G = ch10.amplification_factor(th, a, b, "{{scheme}}")   # (10.24)
  print(np.abs(G).max())                           # {{gmax}} at theta = {{worst}}
  print(ch10.ftcs_stable(a, b))                    # Noye (10.27): {{noye}}
  print(ch10.stability_verdict("{{scheme}}", a, b)["reason"])
  ```
- **Walkthrough (7 steps):** 1. "Same physics, two fates" — "β = 0.50 and 0.51 differ by 2 %. One run is fine, one explodes. Why?"
  `set` β = 0.51 preset, `play: true` · 2. "Errors are waves" — "Any error is a sum of waves (10.20). A linear scheme treats each wave
  alone." `derive: {id: 'D04', step: 5}` · 3. "One number per wave" — "Each step multiplies a wave by G(θ). The orange curve shows G for
  every θ." `readouts: ['G']`, `inspect: true` · 4. "Leaving the circle" — "Outside the circle ∣G∣ > 1: that wave grows every step. At
  β = 0.51 it happens at θ = π — the zigzag." `set` {theta: π}, `derive: {id: 'D05', step: 8}` · 5. "Noye's lens" — "All the (α, β) that
  keep the curve inside: 0 ≤ 4α² ≤ 2β ≤ 1 (10.27). Drag the dot out through each edge." `controls: ['alpha', 'beta']`, `derive: {id:
  'D06', step: 11}` · 6. "Implicit rescue" — "BTCS divides instead of multiplying: ∣G∣ ≤ 1 for any Δt." `set` BTCS preset · 7. "Your turn"
  — "Predict: with α = 0.2, how small can β be? Check with the dot." `controls: ['alpha', 'beta']`.
- **Equations:** `err` "Error equation" ref 'Eq. (10.19)' $\xi^{n+1}_i=(\alpha+\beta)\xi^n_{i-1}+(1-2\beta)\xi^n_i+(\beta-\alpha)\xi^n_{i+1}$ ·
  `G` "Amplification factor" ref 'Eq. (10.24)' $\frac{g^{n+1}}{g^n}=(\alpha+\beta)e^{-\mathrm i\pi k\Delta x}+(1-2\beta)+(\beta-\alpha)
  e^{\mathrm i\pi k\Delta x}$ live "= … at your θ" · `G2` "Its modulus" ref 'Eq. (10.26)' $\big(1-4\beta\sin^2\frac\theta2\big)^2+(2\alpha\sin
  \theta)^2\le1$ · `noye` "Noye's region" ref 'Eq. (10.27)' $0\le4\alpha^2\le2\beta\le1$ live "0.04 ≤ 0.4 ≤ 1 ✓" · `diff` "Pure diffusion"
  ref 'Eq. (10.28)' $\Delta t\le\frac12\frac{\Delta x^2}D$ · `btcs` "BTCS (ours, D07)" $G=[1+4\beta\sin^2\frac\theta2+2\mathrm i\alpha\sin\theta]^{-1}$ ·
  symbols α, β, θ, G, ξ.
- **Check yourself:** (1) "At α = 0, which wave is the first to grow as β passes 0.5?" — "θ = π, the 2Δx zigzag: G(π) = 1 − 4β." `set
  {alpha: 0, beta: 0.52}` · (2) "Set β = 0. Is any α > 0 stable?" — "No: ∣G∣² = 1 + 4α² sin²θ > 1 for every θ ≠ 0, π." `set {beta: 0}` ·
  (3) "With α = 0.3, what is the smallest stable β?" — "2β ≥ 4α² = 0.36 ⇒ β ≥ 0.18." `set {alpha: 0.3}` · (4) "Why can BTCS use β = 100?"
  — "Its G has a denominator whose modulus exceeds 1; stability is not accuracy, though." `set {scheme: 'btcs'}`.
- **Selftest parity rows:** `{name: '|G| ftcs pi/2', js: Gmod(Math.PI/2, 0.1, 0.2, 'ftcs'), py: 'ch10.amplification_modulus(np.pi/2, 0.1,
  0.2, "ftcs")', rtol: 1e-12}` · `{name: 'Im G ftcs', js: G(Math.PI/2, 0.1, 0.2, 'ftcs').im, py: 'ch10.amplification_factor(np.pi/2, 0.1,
  0.2, "ftcs").imag', rtol: 1e-12}` · `{name: 'Gmax beta .51', js: Gmax('ftcs', 0, 0.51), py: 'ch10.max_amplification("ftcs", 0.0, 0.51)',
  rtol: 1e-9}` · `{name: '|G| btcs', js: Gmod(Math.PI, 0, 100, 'btcs'), py: 'ch10.amplification_modulus(np.pi, 0.0, 100.0, "btcs")', rtol:
  1e-12}` · `{name: 'worst theta long wave', js: worstTheta('ftcs', 0.3, 0.1), py: 'ch10.worst_theta("ftcs", 0.3, 0.1)', rtol: 1e-6}` ·
  `{name: 'march 50 steps', js: marchMax(0, 0.51, 50), py: 'ch10.propagate_error(1e-10*np.cos(np.pi*np.arange(64)) +
  1e-10*np.sin(2*np.pi*np.arange(64)/64), 0.0, 0.51, 50)["max_abs"][-1]', rtol: 1e-9}` · invariant `{name: 'Noye edge', js: noye(0.3, 0.18)
  ? 1 : 0, expect: 1}`.
- **Fit plan:** 360×640: status, `plane` | `region` side by side (square, 48 % each) above the `march` log-plot; the ξ bars hidden
  (the dominant wavelength in the status); walkthrough paged. 844×390: three views in one row. Desktop: rows [1, 1].

### E3 · upwind_cfl_advection
- **Title:** "Why must the time step shrink with the grid?" · **Summary:** "Move the Courant number across 1 and watch the characteristic
  leave the stencil and the pulse explode; below 1, upwind smears and MacCormack ripples — diffusion versus dispersion on one clock." ·
  **CORE:** C05, C10 (also N16, N18, N57, C03's modified equation) · **Reference:** `angular_frequency_explorer_1.html` (one clock,
  linked views, end-of-run summary card) and `overfitting_curves.html` (a minimal verdict line that changes with the regime).
- **meta:** `viz:order 3` · `viz:sections 10.2 10.4` · `viz:equations 10.29 10.30 10.100 10.101 10.102` · `viz:fluidpy
  ch10.advect_periodic ch10.transport_1d_step ch10.maccormack_advection_1d ch10.amplification_modulus ch10.phase_error
  ch10.numerical_diffusivity ch10.stability_verdict` · `viz:derivations D08 D18`.
- **Physics:** `advect(profile, C, N, steps, scheme)` exactly as Part C 10.1 (upwind (10.29) with the side from sign(u); MacCormack FB
  (10.101)–(10.102) and BF; Lax–Wendroff; FTCS) ↔ `ch10.advect_periodic`; `Gup(θ, C)`, `Glw(θ, C)` ↔ `ch10.amplification_modulus(θ, C/2,
  0, "upwind" | "lax_wendroff")`; `phase(scheme, C, θ)` ↔ `ch10.phase_error`; `dnum(u, dx, C)` = ∣u∣Δx(1 − C)/2 ↔
  `ch10.numerical_diffusivity(u, dx, 0, "upwind", C)`; growth capped at 10⁶.
- **Views** (rows [1.2, 1]): 1. `pulse` "The pulse going round" (row 0, flex 1.6) — x ∈ [0, 1) periodic, T ∈ [−0.4, 1.4]; exact
  (muted dashed), upwind (teal), MacCormack (purple), FTCS (orange; rose once it exceeds the frame, with "blew up at step n"); a small
  clock "revolution 0.43"; toggles per scheme. 2. `xt` "Stencil and characteristic" (row 0, flex 1) — the x–t diagram around node i:
  stencil points of the chosen scheme at level n (dots), the new point at n + 1, the characteristic line through it with slope 1/u,
  its foot ◆ at x_i − CΔx, the domain of dependence of the stencil shaded (teal inside, rose when ◆ is outside); the label "C = 0.8".
  3. `spec` "∣G∣ and phase speed vs θ" (row 1, `hidePortrait: true`) — two small panels: ∣G(θ)∣ for upwind and LW; relative phase speed
  vs θ; the pulse's dominant θ (2πΔx/width) marked.
- **Controls (≤ 5 visible):** `C` "Courant number $C=u\Delta t/\Delta x$" 0.1 … 1.2, step 0.01, default 0.8 · `schemes` toggles upwind /
  MacCormack / FTCS (default upwind + MacCormack) · `profile` chips square / Gaussian · `N` "Cells $N$" 25 … 400, step 25, default 100
  (optional) · `arr` chips FB / BF (MacCormack arrangement, optional) · transport `rev`.
- **Transport:** `rev` (revolutions) 0 → 3, rate 0.5 rev/s, `end: 'hold'`; end card: "after 3 revolutions: upwind keeps 38 % of the
  height (rms error 0.19); MacCormack 96 % with ripples (rms 0.08)" (numbers live).
- **Presets:** "C = 1: exact shift" {C: 1} · "C = 0.5, upwind smear" {C: 0.5, schemes: upwind} · "C = 0.8, MacCormack ripples" {C: 0.8,
  profile: square} · "C = 1.05: blow-up" {C: 1.05} · "FTCS at any C" {schemes: ftcs, C: 0.5} · "u < 0 as printed (R11)" (a mode flag
  `printed: true` that uses the i − 1 side with u = −1: blows up; its status names slip R11).
- **Status:** "✅ C = 0.80 ≤ 1: the foot is inside the stencil — stable" · "🎯 C = 1: the foot hits a grid point — exact shift" · "💥 C = 1.05 >
  1: the characteristic leaves the stencil — unstable" · FTCS: "💥 FTCS with no diffusion is unstable at every C (C04)" · printed: "💥 the
  printed (10.29) takes the downstream side when u < 0".
- **Readouts:** "C" · "Foot at" (cells) · "Height left (upwind)" · "Height left (MacCormack)" · "D_num (upwind)".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** + end card · **presets** · **status** ·
  **terms** (bars: amplitude loss per revolution upwind vs MacCormack, phase lag) · **modes** (as printed / sign-aware).
- **Explain:**
  0. *What the views show* — "**The pulse**: grey dashed is the exact answer (the pulse just moves); teal is first-order upwind, purple
     MacCormack. **Stencil and characteristic**: the line is the path information travels along; its foot ◆ must fall between the dots the
     scheme uses. **∣G∣ and phase** (hidden on phones): how much each wave keeps and how fast it moves, per step."
  1. *Courant number* — "C = uΔt/Δx = **0.8** (10.30); the foot of the characteristic is C cells upstream: **0.8** cells ⇒ inside the
     two-point upwind stencil" boxed.
  2. *Upwind loses height* — "(10.29): G = 1 − C(1 − e^{−iθ}); for the pulse's main wave θ = 2πΔx/λ = **…** rad (λ = twice the pulse width): ∣G∣ = **…** per step; one
     revolution = N/C = **125** steps ⇒ ∣G∣^125 = **…**; the modified equation adds D_num = ∣u∣Δx(1 − C)/2 = **1e-3** (with u = 1, Δx = 0.01, C = 0.8)"
     boxed.
  3. *MacCormack keeps height but lags* — "(10.101)–(10.102) = Lax–Wendroff for this equation (D18); ∣G∣² = 1 − 4C²(1 − C²)sin⁴(θ/2) =
     **…** ≈ 1; relative phase speed at θ = π/4: **0.968** — short waves arrive late, as ripples behind the pulse."
  4. *At the current time* — "revolution **r**: height upwind **…**, MacCormack **…**; rms errors **…**, **…**" (live).
  5. *Beyond C = 1* — "∣G(π)∣² = 1 − 2C(1 − C)·2 = **1.44** at C = 1.1: the zigzag grows 20 % per step; the characteristic's foot is outside
     the stencil — the scheme cannot know where the information came from."
  6. *Reading the current setting* — C < 1: "Stable, but every scheme pays: first order in amplitude (smear), second order in phase
     (ripples)." · C = 1: "Exact shift: G = e^{−iθ} — a lucky special case that only exists for constant u." · C > 1: "CFL violated:
     reduce Δt — which is why a finer grid forces a smaller time step (and why doubling a weather model's resolution costs 8×)." ·
     printed: "The printed (10.29) assumes u > 0; with u < 0 it takes the downstream side — unstable for every C (slip R11)."
- **Derivation tab:** **D08** (9 steps) `view: 'xt'`; goal `set` {C: 0.5}; step 6 `live` "∣G(π)∣² = 1 − 2C(1−C)(1−cos π) = 0 at C = 0.5";
  step 8 `set` {C: 1.05} `watch` "the foot leaves the stencil". **D18** (14 steps, ★★★) `view: 'pulse'`; goal `set` {C: 0.8, schemes:
  ['maccormack']}; step 5 `live` "the one-step numbers of the tiny example"; step 12 `live` "local error −(uΔtΔx²/6)(1 − C²)T_xxx: zero at
  C = 1"; step 14 `set` {C: 1.05}. Interpret: `s => "At C = " + C + ": upwind keeps " + hu + " of the height per revolution, MacCormack "
  + hm + "."`.
- **Code:**
  ```python
  C, N = {{C}}, {{N}}                              # Courant number, cells
  up = ch10.advect_periodic("{{profile}}", C=C, n_cells=N, scheme="upwind")
  mc = ch10.advect_periodic("{{profile}}", C=C, n_cells=N, scheme="maccormack")
  print(up["amplitude_ratio"], mc["amplitude_ratio"])   # {{hu}}, {{hm}}
  print(ch10.numerical_diffusivity(1.0, 1/N, 0.0, "upwind", C))  # {{dnum}}
  print(ch10.stability_verdict("upwind", C/2, 0.0)["reason"])
  ```
- **Walkthrough (6 steps):** 1. "Why shrink Δt?" — "A pulse rides a current. Finer grid, smaller time step — why? Press ▶." `play:
  true` · 2. "Where information comes from" — "The value at the new point came from C cells upstream (the ◆). Upwind uses the upstream
  neighbour." `derive: {id: 'D08', step: 2}` · 3. "The edge" — "Drag C past 1: the ◆ leaves the stencil and the pulse explodes. That is
  CFL (10.30)." `controls: ['C']`, `set` {C: 1.05} · 4. "Stable is not exact" — "At C = 0.5 upwind smears: it adds a diffusivity
  uΔx(1 − C)/2." `set` {C: 0.5}, `readouts: ['dnum']` · 5. "Second order ripples" — "MacCormack predicts forward, corrects backward:
  Lax–Wendroff. Height kept, but short waves lag." `set` {C: 0.8, profile: 'square'}, `derive: {id: 'D18', step: 9}` · 6. "Your turn" —
  "Predict which scheme wins on the Gaussian at C = 0.9, then play three revolutions." `controls: ['C', 'profile']`.
- **Equations:** `up` "Upwind" ref 'Eq. (10.29)' $T^{n+1}_i=T^n_i-2\alpha(T^n_i-T^n_{i-1})$ (u > 0) · `cfl` "CFL condition" ref 'Eq. (10.30)'
  $u\frac{\Delta t}{\Delta x}\le1$ live "C = 0.8 ✓" (note: $\lvert u\rvert$ with the side from sign(u)) · `pred` "Predictor" ref 'Eq. (10.101)'
  $\mathbf U^*_{i}=\mathbf U^n_i-\frac{\Delta t}{\Delta x}(\mathbf E^n_{i+1}-\mathbf E^n_i)$ · `corr` "Corrector" ref 'Eq. (10.102)'
  $\mathbf U^{n+1}_i=\frac12\big[\mathbf U^n_i+\mathbf U^*_i-\frac{\Delta t}{\Delta x}(\mathbf E^*_i-\mathbf E^*_{i-1})\big]$ · `lw` "Their sum
  (ours, D18)" $T^{n+1}_i=T_i-\frac C2(T_{i+1}-T_{i-1})+\frac{C^2}2(T_{i+1}-2T_i+T_{i-1})$ · symbols C, u, Δt, Δx, θ.
- **Check yourself:** (1) "At what C does upwind become exact, and why?" — "C = 1: the foot of the characteristic lands on the upstream
  node, so copying it is exact." `set {C: 1}` · (2) "Which scheme undershoots below zero on the square pulse?" — "MacCormack: dispersion
  makes short waves lag and overshoot; upwind never goes below zero for C ≤ 1 (it is a weighted average)." `set {C: 0.8}` · (3) "Double
  N at fixed C. What happens to upwind's smearing per revolution?" — "It shrinks: D_num ∝ Δx — first order." · (4) "Why is FTCS useless
  here at any C?" — "With no diffusion ∣G∣² = 1 + C² sin²θ > 1 (C04)."
- **Selftest parity rows:** `{name: 'upwind |G| pi C.5', js: Gup(Math.PI, 0.5), py: 'ch10.amplification_modulus(np.pi, 0.25, 0.0,
  "upwind")', atol: 1e-12, rtol: 0}` · `{name: 'LW phase', js: phase('lax_wendroff', 0.8, Math.PI/4), py: 'ch10.phase_error("lax_wendroff",
  0.8, np.pi/4)', rtol: 1e-10}` · `{name: 'MacCormack one step', js: mcStep([0, 0, 1, 0, 0], 0.5)[1], py:
  'ch10.maccormack_advection_1d(np.array([0.0, 0, 1, 0, 0]), 0.5, 1)[1]', rtol: 1e-12}` · `{name: 'upwind height 1 rev', js:
  advect('gauss', 0.5, 100, 'upwind').amplitude_ratio, py: 'ch10.advect_periodic("gauss", C=0.5, n_cells=100, n_rev=1.0,
  scheme="upwind")["amplitude_ratio"]', rtol: 1e-9}` · `{name: 'MacCormack height', js: advect('square', 0.8, 100,
  'maccormack').amplitude_ratio, py: 'ch10.advect_periodic("square", C=0.8, n_cells=100, n_rev=1.0, scheme="maccormack")["amplitude_ratio"]',
  rtol: 1e-9}` · `{name: 'D_num', js: dnum(1, 0.01, 0.5), py: 'ch10.numerical_diffusivity(1.0, 0.01, 0.0, "upwind", 0.5)', rtol: 1e-12}`.
- **Fit plan:** 360×640: status, `pulse` (58 %) over `xt` (42 %), `spec` hidden (the phase number in the `pulse` title); scheme
  toggles as chips. 844×390: `pulse` | `xt`. Desktop rows [1.2, 1].

### E4 · cell_peclet_wiggles
- **Title:** "Why does a centred scheme give negative temperatures?" · **Summary:** "Slide R_cell across 2: the centred discrete solution
  rʲ flips sign as the root r turns negative, while upwind stays smooth by solving a more diffusive problem." · **CORE:** C09 (also N19,
  N42–N49, N47, N50) · **Reference:** `forced_damped_vibrations.html` (regime status + interpretation) and `overfitting_curves.html`
  (a minimal two-slider figure with a verdict line).
- **meta:** `viz:order 4` · `viz:sections 10.2 10.4` · `viz:equations 10.31 10.84 10.85 10.86 10.87 10.88 10.89 10.90 10.91 10.92 10.93
  10.94` · `viz:fluidpy ch10.steady_cd_exact ch10.steady_cd_fd ch10.steady_cd_discrete_exact ch10.discrete_root ch10.wiggle_indicator
  ch10.numerical_diffusivity ch10.cd_layer_thickness` · `viz:derivations D15 D16 D17`.
- **Physics:** `exact(x, R)` scaled form ↔ `ch10.steady_cd_exact`; `fd(n, R, scheme, grid)` Thomas solve of (10.91)/(10.93) (and the
  stretched grid of C.10 10.3) ↔ `ch10.steady_cd_fd`; `root(Rc, scheme)` ↔ `ch10.discrete_root`; `closed(n, R, scheme)` rʲ form ↔
  `ch10.steady_cd_discrete_exact`; `wig(T)` ↔ `ch10.wiggle_indicator`; `dnum(u, dx, D)` ↔ `ch10.numerical_diffusivity(u, dx, D,
  "upwind_steady")`; `modexact(x, R, Rc)` = exact with R/(1 + 0.5 R_cell) (the modified-equation ghost).
- **Views** (rows [1.2, 1]): 1. `profile` "T near the hot wall" (row 0, flex 1.5) — x/L ∈ [0, 1] or zoomed to [1 − 6/R, 1] (toggle
  `zoom`), T ∈ [−0.6, 1.1]; exact (10.86) (muted dashed), centred dots joined (orange; negative values in rose with a rose band under
  T = 0), upwind dots (teal), modified-equation solution (teal dashed); the layer thickness L/R marked by a bracket; the e⁻¹ level
  dotted. 2. `roots` "The discrete root r" (row 0, flex 1) — r vs R_cell ∈ [0, 8]: centred (orange: r = (1 + R_cell/2)/(1 − R_cell/2),
  → +∞ at 2⁻ and from −∞ at 2⁺), upwind (teal: 1 + R_cell), the band r < 0 shaded rose, current dots; a mini strip below the axis showing
  the signs of r⁰ … r⁶ as ± squares. 3. `diff` "Physical vs numerical diffusivity" (row 1, `hidePortrait: true`) — two bars D (blue) and
  0.5 R_cell D (rose) stacked to D(1 + 0.5 R_cell) — and a small log–log error-vs-n panel for both schemes (slopes 2 and 1).
- **Controls (≤ 5 visible):** `logR` "Péclet $R=uL/D$ (log)" 0 … 2.7 (R = 1 … 500), step 0.01, default log10 40 · `n` "Cells $n$" 4 … 100,
  step 1, default 10 · `schemes` toggles centred / upwind (both on) · `grid` chips uniform / stretched (optional) · `zoom` toggle
  "zoom to the layer" (optional).
- **Transport:** none (a steady problem); instead a **"sweep R_cell"** play button that animates n from 40 down to 4 at fixed R (stops
  at the first wiggle and holds with a card "first negative T at R_cell = 2.1").
- **Presets:** "R_cell = 1" {R: 10, n: 10} · "R_cell = 2 (edge)" {R: 20, n: 10} · "R_cell = 4" {R: 40, n: 10} · "R_cell = 20" {R: 200,
  n: 10} · "stretched, R_cell = 4" {R: 40, n: 10, grid: 'stretched'} · "tiny example" {R: 16, n: 4}.
- **Status:** "✅ R_cell = 1.0 < 2: centred r = 3.0 > 0, no wiggles" · "⚠️ R_cell = 4.0 > 2: centred r = −3.0 — wiggles (min T = −0.35)" ·
  "🟦 upwind: smooth, but D_eff = D(1 + 0.5 R_cell) = 3.0 D — the layer is 3× too thick" · stretched: "✅ stretched grid: smallest cell h_min (live) ⇒
  local R_cell = u h_min/D (live) < 2 — no wiggles at the same n".
- **Readouts:** "R_cell = R/n" · "r centred" · "r upwind" · "min T (centred)" · "D_num/D".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **presets** · **status** · **terms** (D and numerical D bars) ·
  **inspector** (click node j: "T_j = (r^j − 1)/(r^n − 1) = ((−3)³ − 1)/((−3)⁴ − 1) = −28/80 = −0.35") · **modes** (uniform / stretched).
- **Explain:**
  0. *What the views show* — "**T near the hot wall**: dashed grey is the exact temperature — flat, then a thin rise at x = L; orange dots
     are centred differences, teal upwind. **The discrete root**: each scheme's solution is rʲ; the rose band is where r < 0 (signs
     alternate). **Diffusivity** (hidden on phones): the physical D and what upwinding adds."
  1. *Two Péclet numbers* — "R = uL/D = **40** (10.87); R_cell = R/n = **4.0** (10.31): the layer thickness L/R = **0.025** L is
     **0.25** of a cell (10.92)" boxed; why: (10.89) δ/L = O(1/R).
  2. *The exact layer* — "T(x) = (e^{Rx/L} − 1)/(e^R − 1) (10.86) ≈ e^{−R(1−x/L)} (10.88): one layer thickness from the wall T = e⁻¹ =
     0.368, two thicknesses e⁻² = 0.135; at the last interior node x = L − Δx: T = e^{−R_cell} = **0.018**."
  3. *Centred root* — "(10.91) with T_j = r^j: (1 − R_cell/2)r² − 2r + (1 + R_cell/2) = 0 ⇒ r = 1 or r = (1 + R_cell/2)/(1 − R_cell/2) =
     **−3.0**; T_j = (r^j − 1)/(r^n − 1): T at the last interior node = **−0.35**" boxed (table of j, r^j, T_j with the clicked row lit).
  4. *Upwind root and its price* — "(10.93): r = 1 + R_cell = **5.0** > 0 always; but (10.94) says it solves uT_x = D(1 + 0.5R_cell)T_xx:
     D_num = 0.5 R_cell D = **2.0** D; effective R = R/(1 + 0.5 R_cell) = **13.3**."
  5. *Refinement* (the hidden view's numbers) — "error ∝ n^{−2} (centred, once R_cell < 2) and ∝ n^{−1} (upwind): at n = **…**: **…**, **…**."
  6. *Reading the current setting* — R_cell < 2: "The grid resolves the layer: centred differences are second order and accurate." ·
     R_cell ≈ 2: "On the edge: the centred root is infinite — the scheme copies the boundary value into one node." · R_cell > 2: "The layer is
     thinner than half a cell: the discrete balance at the last node cannot hold with smooth values, so it alternates. The wiggle is the
     grid confessing — refine, stretch, or accept upwind's extra diffusion." · stretched: "Same number of cells, placed where the layer is:
     local R_cell < 2 near the wall."
- **Derivation tab:** **D15** (8 steps) `view: 'profile'`; step 7 `set` {zoom: true} `live` "T(L − δ) = e⁻¹ = 0.368". **D16** (12 steps)
  `view: 'roots'`; step 7 `live` "r = (1 + 2)/(1 − 2) = −3 at R_cell = 4"; step 11 `set` {R: 40, n: 10}. **D17** (8 steps) `view:
  'profile'`; step 8 `live` "D(1 + 0.5 × 4) = 3D". Interpret: `s => "R_cell = " + Rc + ": centred r = " + rc + (rc < 0 ? " — wiggles" :
  " — smooth") + "; upwind adds " + dn + " D."`.
- **Code:**
  ```python
  n, R = {{n}}, {{R}}                              # cells, global Peclet uL/D
  xc = ch10.steady_cd_fd(n, R, "central")          # (10.91): T = {{tc}}
  xu = ch10.steady_cd_fd(n, R, "upwind")           # (10.93): T = {{tu}}
  print(ch10.discrete_root(R/n), ch10.discrete_root(R/n, "upwind"))  # {{rc}}, {{ru}}
  print(ch10.wiggle_indicator(xc["T"])["wiggles"]) # {{wig}}
  print(ch10.numerical_diffusivity(1.0, 1/n, 1/R, "upwind_steady"))  # {{dn}} (u = 1, L = 1)
  ```
- **Walkthrough (6 steps):** 1. "Negative temperatures?" — "Heat flows to a hot wall. Can a correct scheme give T < 0? Look at the
  orange dots." `set` R_cell = 4 preset · 2. "The exact layer" — "The temperature jumps up within L/R of the wall (10.88). Here that is a
  quarter of a cell." `derive: {id: 'D15', step: 7}` · 3. "rʲ" — "The centred equation (10.91) is solved exactly by rʲ. Negative r
  alternates signs." `derive: {id: 'D16', step: 8}`, `inspect: true` · 4. "The switch at 2" — "Press 'sweep': the first negative value
  appears exactly when R_cell passes 2." `play: true` · 5. "Upwind's price" — "Upwind is smooth but solves D(1 + 0.5R_cell) (10.94): the
  rose bar." `terms: true`, `derive: {id: 'D17', step: 8}` · 6. "Your turn" — "Keep n = 10, R = 40. Predict: does the stretched grid
  wiggle? Try it." `controls: ['grid']`.
- **Equations:** `cd` "Steady problem" ref 'Eq. (10.84)–(10.85)' $uT_x=DT_{xx},\ T(0)=0,\ T(L)=1$ · `ex` "Exact solution" ref 'Eq. (10.86)'
  $T=\frac{e^{Rx/L}-1}{e^R-1}$ live · `layer` "Layer thickness" ref 'Eq. (10.89)' $\frac\delta L=O\big(\frac1{\lvert R\rvert}\big)$ · `cen`
  "Centred" ref 'Eq. (10.91)' $0.5R_{cell}(T_{j+1}-T_{j-1})=T_{j+1}-2T_j+T_{j-1}$ · `upw` "Upwind" ref 'Eq. (10.93)' $R_{cell}(T_j-T_{j-1})=
  T_{j+1}-2T_j+T_{j-1}$ · `mod` "What upwind solves" ref 'Eq. (10.94)' $uT_x=D(1+0.5R_{cell})T_{xx}$ · `pe` "Cell Péclet" ref 'Eq. (10.31)'
  $R_{cell}=u\frac{\Delta x}D\le2$ · symbols R, R_cell, r, δ, D_num.
- **Check yourself:** (1) "At R_cell = 2 exactly, what is the centred root?" — "Infinite (1 − R_cell/2 = 0): the equation degenerates;
  just above 2, r is large and negative." `set {R: 20, n: 10}` · (2) "Double n at fixed R = 40. Do the wiggles go?" — "Yes: R_cell = 2 at
  n = 20 — on the edge; n = 21 is clean." · (3) "How thick does upwind make the layer at R_cell = 20?" — "D_eff = 11D: eleven times the
  physical layer." `set {R: 200, n: 10}` · (4) "Why are the book's words 'forward difference' for (10.93) wrong?" — "T_j − T_{j−1} uses the
  node behind: a backward (upwind) difference (slip R3)."
- **Selftest parity rows:** `{name: 'exact R=4 x=.75', js: exact(0.75, 4), py: 'ch10.steady_cd_exact(np.array([0.75]), 4.0)[0]', rtol:
  1e-12}` · `{name: 'exact R=1e4 finite', js: exact(0.999, 1e4), py: 'ch10.steady_cd_exact(np.array([0.999]), 1e4)[0]', rtol: 1e-9}` ·
  `{name: 'root centred 4', js: root(4, 'central'), py: 'ch10.discrete_root(4.0, "central")', rtol: 1e-12}` · `{name: 'centred T3', js:
  fd(4, 16, 'central', 'uniform').T[3], py: 'ch10.steady_cd_fd(4, 16.0, "central")["T"][3]', rtol: 1e-12}` · `{name: 'upwind T3', js: fd(4,
  16, 'upwind', 'uniform').T[3], py: 'ch10.steady_cd_fd(4, 16.0, "upwind")["T"][3]', rtol: 1e-12}` · `{name: 'stretched min', js: fd(10, 40,
  'central', 'stretched').T.reduce((a, b) => Math.min(a, b)), py: 'ch10.steady_cd_fd(10, 40.0, "central", grid="stretched")["T"].min()',
  rtol: 1e-9, atol: 1e-12}`.
- **Fit plan:** 360×640: status, `profile` (60 %) over `roots` (40 %), `diff` hidden (D_num/D in a readout); 844×390 side by side;
  desktop rows [1.2, 1].

### E5 · fem_hat_assembly
- **Title:** "Where do finite-element matrices come from?" · **Summary:** "Step through the elements: each adds a 2 × 2 block into the
  global mass and stiffness matrices; the assembled rows are centred FD stencils with a ⅙–⅔–⅙ mass, and the weighted hats sum to the
  FE solution." · **CORE:** C06, C07, C08 (also N21–N38, N32, N33) · **Reference:** `stride_padding_playground.html` (a size formula
  with numbers plugged in; click an output cell) and `amplitude_phase_second_order_II_3.html` (numbered live Explain).
- **meta:** `viz:order 5` · `viz:sections 10.3` · `viz:equations 10.36 10.49 10.55 10.56 10.58 10.59 10.63 10.64 10.65 10.74 10.75 10.76
  10.77 10.78` · `viz:fluidpy ch10.hat ch10.element_matrices_linear ch10.element_force ch10.connectivity ch10.assemble_1d
  ch10.assembly_trace ch10.solve_steady ch10.interior_stencil ch10.steady_cd_fd ch10.steady_cd_exact` · `viz:derivations D09 D11 D12 D14`.
- **Physics:** `hat(x, nodes, A)` ↔ `ch10.hat`; `elem(h, u, D)` → {m, k} ↔ `ch10.element_matrices_linear`; `force(e, nel, h, u, D, g, q)`
  ↔ `ch10.element_force`; `conn(nel)` ↔ `ch10.connectivity`; `assemble(nel, u, D, upTo)` scatter-add of elements 1…upTo (dense, n ≤ 20)
  ↔ `ch10.assembly_trace(...)[e]["K_after"]`; `solveSteady(nel, u, D, bc)` Thomas solve after eliminating the Dirichlet node(s) ↔
  `ch10.solve_steady`; `row(h, u, D)` ↔ `ch10.interior_stencil`; `fdCentral(nel, R)` ↔ `ch10.steady_cd_fd`; `exact(x, R)` ↔
  `ch10.steady_cd_exact`. L = 1, u = 1, D = 1/R.
- **Views** (rows [1.2, 1]): 1. `mesh` "Hats and the FE solution" (row 0, flex 1.4) — x ∈ [0, 1]; all hats N_A (muted), the probe test
  function N_A (amber fill), the weighted hats d_A N_A (teal, drawn only for the elements assembled so far), their sum Tʰ (blue, bold)
  and the exact (10.86) (muted dashed) once assembly is complete; the element being assembled shaded amber with its two local shapes
  N₁, N₂ drawn over it; node dots with d_A. 2. `matrix` "Global K (and M)" (row 0, flex 1) — the (n + 1) × (n + 1) matrix as coloured
  cells (diverging colour, 0 white), numbers printed for n ≤ 8; the current element's 2 × 2 block outlined amber; a chip switches K / M;
  eliminated Dirichlet rows greyed; pointer: click a cell → inspector. 3. `row` "Nodal values and the interior row" (row 1,
  `hidePortrait: true`) — FE dots (blue), centred FD crosses (orange), exact line (muted); beneath, the interior row of node A written
  with numbers: "(1/6, 2/3, 1/6)·h d/dt(T) + (−1.5, 2, −0.5)·T = 0".
- **Controls (≤ 5 visible):** `nel` "Elements $n_{el}$" 2 … 20, step 1, default 4 · `R` "Péclet $R=uL/D$" 0.1 … 40, step 0.1, default 4 ·
  `e` "Assembled up to element" 0 … n_el (the transport) · `mass` chips consistent / lumped (optional) · `bc` chips "T(1) = 1" /
  "∂T/∂x(1) = 1 (natural)" (optional) · `A` "Test function $N_A$" 1 … n_el − 1 (optional).
- **Transport:** `e` 0 → n_el, rate 1 element/s, `end: 'hold'`; end card: "assembled: K is tridiagonal; FE = centred FD to 1e-15".
- **Presets:** "n = 4 by hand" {nel: 4, R: 4} · "pure diffusion" {nel: 6, R: 0.1} · "R = 10: FE wiggles too" {nel: 4, R: 10} · "fine mesh"
  {nel: 20, R: 10} · "natural BC" {bc: 'neumann', nel: 8, R: 4} · "lumped mass" {mass: 'lumped', nel: 4}.
- **Status:** assembling: "🧱 element 3 of 4 adds its block to rows/columns 2 and 3" · complete: "✅ assembled — FE equals centred FD
  (max difference 1.1e-16); R_cell = uh/D = 1.0" · "⚠️ R_cell = 2.5 > 2: the FE solution wiggles exactly like centred FD (C09)" · natural:
  "🔁 natural BC: the slope at x = 1 is learned, … vs q = 1 at n = 8 (live)".
- **Readouts:** "h" · "R_cell" · "K row of node A" · "d at probe node" · "max ∣FE − FD∣".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** (assembly) · **presets** · **status** ·
  **inspector** (click a matrix cell: "K₂₂ = k⁽²⁾₂₂ + k⁽³⁾₁₁ = (u/2 + D/h) + (−u/2 + D/h) = 1.5 + 0.5 = 2.0"; off-diagonal "K₂₃ = k⁽³⁾₁₂ =
  u/2 − D/h = −0.5") · **modes** (consistent / lumped mass; Dirichlet / natural end).
- **Explain:**
  0. *What the views show* — "**Hats**: every grey tent is one basis function N_A (1 at its node, 0 at the others); teal are the tents
     scaled by the nodal values; blue is their sum Tʰ — the FE solution. The amber element is the one being assembled. **Global K**: the
     matrix fills in one 2 × 2 block per element. **Nodal values** (hidden on phones): FE, FD and exact side by side."
  1. *The element and its map* — "h = 1/n_el = **0.25**; on the parent element ξ ∈ [−1, 1]: N₁ = ½(1 − ξ), N₂ = ½(1 + ξ) (10.64);
     x(ξ) = (h/2)ξ + x_mid (10.65) ⇒ dx = (h/2)dξ, dN₁/dx = −1/h = **−4**, dN₂/dx = **+4**."
  2. *Element mass* — "m_ab = ∫N_aN_b dx (10.75) = (h/6)[[2, 1], [1, 2]] = [[**0.0833**, **0.0417**], [**0.0417**, **0.0833**]]" boxed (lumped:
     (h/2)I).
  3. *Element stiffness* — "k_ab = u∫N_{b,x}N_a dx + D∫N_{b,x}N_{a,x} dx (10.76) = (u/2)[[−1, 1], [−1, 1]] + (D/h)[[1, −1], [−1, 1]] =
     [[**0.5**, **−0.5**], [**−1.5**, **1.5**]]" boxed; why: the convective block is not symmetric (only v is differentiated in a(w, v)).
  4. *Scatter-add* — "element e maps local nodes (1, 2) to global (e − 1, e) (10.78); node A gets k⁽ᴬ⁾₂₁, k⁽ᴬ⁾₂₂ + k⁽ᴬ⁺¹⁾₁₁, k⁽ᴬ⁺¹⁾₁₂ =
     (**−1.5**, **2.0**, **−0.5**) — (10.63) × h: (−u/2 − D/h, 2D/h, u/2 − D/h) ✓" (table of the elements that touched the row, current lit).
  5. *The solution* — "K d = F with T(0) = 0, T(1) = 1 moved to F: d = (**0.025**, **0.1**, **0.325**); exact (10.86): (**0.0321**, **0.1192**,
     **0.3561**); centred FD (10.91): identical to FE (difference **1e-16**)."
  6. *Reading the current setting* — R_cell < 2: "FE and FD are the same equations here — FE's advantage is geometry, not accuracy in
     1-D." · R_cell > 2: "The FE solution wiggles because its equations are the centred ones (C09); the Petrov–Galerkin fix (SUPG, N50)
     tilts the test functions upstream." · natural BC: "The Dirichlet end is exact on every mesh (built into the space); the flux end is
     only approached as h → 0 (natural)." · lumped: "Lumping replaces the ⅙–⅔–⅙ average by the node value: the time-dependent scheme
     becomes FD exactly; steady solutions do not change."
- **Derivation tab:** **D09** (9 steps) `view: 'mesh'`, goal `set` {bc: 'neumann'}; step 6 `watch` "the test function vanishes at x = 0".
  **D11** (13 steps) `view: 'matrix'`; step 11 `set` {e: nel} `watch` "the rows are the n equations". **D12** (10 steps) `view: 'row'`;
  step 4 `live` "∫N_A² = 2h/3 = 0.1667"; step 10 `live` "row (−1.5, 2, −0.5)". **D14** (11 steps) `view: 'matrix'`; step 8 `set` {e: 1};
  step 10 `set` {e: 3} `live` "node 2 collects (−1.5, 2.0, −0.5)". Interpret: `s => "n = " + nel + ", R = " + R + ": the assembled row is ("
  + rowtxt + "), and FE = FD to " + dmax + "."`.
- **Code:**
  ```python
  nodes = np.linspace(0, 1, {{nel}} + 1)          # the mesh, h = {{h}}
  m, k = ch10.element_matrices_linear({{h}}, 1.0, {{D}})   # one element's blocks
  M, K, F = ch10.assemble_1d(nodes, 1.0, {{D}}, dirichlet_right=1.0)
  print(K.toarray()[1])                            # row of node 2: {{row}}
  T = ch10.solve_steady(nodes, 1.0, {{D}}, T_L=1.0)["T"]
  print(T)                                         # {{T}}  (= centred FD)
  ```
- **Walkthrough (7 steps):** 1. "Tents, not points" — "FE builds T from tents. Where do its matrices come from? Press ▶ to assemble."
  `play: true` · 2. "Weak form" — "Multiply by a test tent, integrate by parts: only first derivatives remain (10.36)." `derive: {id: 'D09',
  step: 6}` · 3. "One element" — "On one element only two tents live. Their overlaps give a 2 × 2 block (10.75)–(10.76)." `set` {e: 1},
  `derive: {id: 'D14', step: 6}` · 4. "Scatter-add" — "Each block adds into rows e − 1 and e. Click a cell to see who contributed."
  `set` {e: 3}, `inspect: true` · 5. "A familiar row" — "Node 2's row is the centred FD stencil — but the time term averages ⅙–⅔–⅙
  (10.63)." `derive: {id: 'D12', step: 10}`, `readouts: ['row']` · 6. "FE = FD, here" — "Steady and uniform: identical solutions. Push R
  to 10: both wiggle." `set` {R: 10} · 7. "Your turn" — "Switch to the natural end: predict whether the slope at x = 1 is exact on a
  coarse mesh." `controls: ['bc', 'nel']`.
- **Equations:** `weak` "Weak form" ref 'Eq. (10.36)' $\int_0^LT_tw\,dx+u\int_0^LT_xw\,dx+D\int_0^LT_xw_x\,dx=Dq\,w(L)$ · `mat` "Matrix form"
  ref 'Eq. (10.58)' $\mathbf M\dot{\mathbf d}+\mathbf K\mathbf d=\mathbf F$ with $M_{AB}=\int N_AN_B\,dx$ (10.55) · `hat` "Hat function" ref
  'Eq. (10.59)' $N_A=\frac{x-x_{A-1}}{x_A-x_{A-1}}$ on the left, $\frac{x_{A+1}-x}{x_{A+1}-x_A}$ on the right · `row` "Interior row" ref
  'Eq. (10.63)' $\frac d{dt}\big(\frac{T_{A-1}}6+\frac{2T_A}3+\frac{T_{A+1}}6\big)+\frac u{2h}(T_{A+1}-T_{A-1})-\frac D{h^2}(T_{A-1}-2T_A+T_{A+1})=0$ ·
  `elem` "Element blocks" ref 'Eq. (10.75)–(10.77)' $m^e_{ab}=\int_{\Omega^e}N_aN_b\,dx$, $k^e_{ab}=u\int N_{b,x}N_a+D\int N_{b,x}N_{a,x}$ live
  "k = [[0.5, −0.5], [−1.5, 1.5]]" · symbols h, N_A, d_A, M, K, F.
- **Check yourself:** (1) "Why is K tridiagonal?" — "Hats overlap only with their neighbours, so ∫ of a product is zero unless ∣A − B∣ ≤ 1."
  · (2) "Which entry of node 2's row gets contributions from two elements?" — "The diagonal: node 2 is local node 2 of element 2 and local
  node 1 of element 3." `set {nel: 4, e: 3}` · (3) "Is K symmetric? Why not?" — "No: the convective part differentiates only the trial
  function, (u/2)[[−1, 1], [−1, 1]]." · (4) "Does lumping the mass change the steady solution?" — "No — M multiplies ḋ, which is 0 at steady
  state." `set {mass: 'lumped'}`.
- **Selftest parity rows:** `{name: 'k21', js: elem(0.25, 1, 0.25).k[1][0], py: 'ch10.element_matrices_linear(0.25, 1.0, 0.25)[1][1][0]',
  rtol: 1e-12}` · `{name: 'm11', js: elem(0.25, 1, 0.25).m[0][0], py: 'ch10.element_matrices_linear(0.25, 1.0, 0.25)[0][0][0]', rtol: 1e-12}` ·
  `{name: 'steady T2', js: solveSteady(4, 1, 0.25, 'dirichlet')[2], py: 'ch10.solve_steady(np.linspace(0.0, 1.0, 5), 1.0, 0.25, T_L=1.0)["T"][2]',
  rtol: 1e-12}` · `{name: 'row stiff', js: row(0.25, 1, 0.25).stiff[0], py: 'ch10.interior_stencil(0.25, 1.0, 0.25)["stiff"][0]', rtol: 1e-12}`
  · `{name: 'hat value', js: hat(0.3, [0, 0.25, 0.5, 0.75, 1], 1), py: 'ch10.hat(np.array([0.3]), np.linspace(0.0, 1.0, 5), 1)[0]', rtol:
  1e-12}` · invariant `{name: 'FE = FD', js: maxDiff(solveSteady(10, 1, 0.1), fdCentral(10, 10).T), expect: 0, atol: 1e-13}`.
- **Fit plan:** 360×640: status, `mesh` (55 %) over `matrix` (45 %; numbers hidden above n = 6), `row` hidden (the row in a readout);
  844×390 side by side; desktop rows [1.2, 1].

### E6 · mac_projection_staggered
- **Title:** "What is the pressure doing in incompressible flow?" · **Summary:** "Step through predictor → Poisson → correction on a small
  staggered grid: the divergence appears and vanishes to round-off; switch to the checkerboard to see why pressure lives at centres and
  velocity on faces." · **CORE:** C11, C12 (also N41, N60–N74, N68 checkerboard) · **Reference:** `pixels_as_parameters.html` (an
  algorithm in slow motion, stages you can step back and forth) and `angular_frequency_explorer_1.html` (modes).
- **meta:** `viz:order 6` · `viz:sections 10.4` · `viz:equations 10.115 10.116 10.117 10.118 10.121 10.122 10.123 10.124 10.125 10.126` ·
  `viz:fluidpy ch10.MacGrid ch10.divergence ch10.gradient ch10.pressure_poisson_matrix ch10.project ch10.projection_stages
  ch10.mac_projection_summary ch10.collocated_gradient ch10.checkerboard ch10.gradient_null_space` · `viz:derivations D19 D20 D21`.
- **Physics:** staggered arrays `p[ny][nx]`, `u[ny][nx+1]`, `v[ny+1][nx]` (C.10 10.5 fields); `div` ↔ `ch10.divergence`; `grad` ↔
  `ch10.gradient`; `poissonLU(n)` banded LU of the corner-pinned matrix (10.124), cached per n ↔ `ch10.pressure_poisson_matrix`;
  `project(us, vs, dt)` (rhs = div/Δt, solve, subtract the mean, correct (10.121)–(10.122)) ↔ `ch10.project` / `ch10.projection_stages`;
  `sor(k)` k SOR sweeps (ω = 1.8) from p = 0 with the residual history ↔ `ch10.solve_pressure(method="sor")`; `summary(n, kind)` ↔
  `ch10.mac_projection_summary`; `checker(n)` ↔ `ch10.checkerboard`; `colGrad(p)` ↔ `ch10.collocated_gradient`; null-space dimensions
  from the embedded values of `ch10.gradient_null_space` (n = 4, 8).
- **Views** (rows [1.25, 1]): 1. `grid` "The staggered grid" (row 0, flex 1.3, `equal: true`) — n × n cells (4 … 16); p as an orange
  heatmap at centres (stage ≥ 2), u and v as blue arrows on the faces, each cell's border coloured by its divergence (rose ±, white 0),
  wall faces thick black (never corrected); stage banner "① uⁿ ② u* ③ p ④ uⁿ⁺¹"; pointer: click a cell → inspector. 2. `press` "Pressure
  and correction" (row 0, flex 1, `equal: true`) — p contours with the correction arrows −Δt∇p (orange) on the faces; in **checkerboard**
  mode: left half the collocated grid with the zigzag p and its (zero) centred gradient, right half the staggered grid with ±2/Δx face
  gradients. 3. `hist` "How much divergence is left" (row 1, `hidePortrait: true`) — bars of log₁₀ max∣∇·u∣ for the four stages (rose,
  the last ≈ −15); with the SOR solver, the residual vs sweep on log axes and the divergence left after k sweeps.
- **Controls (≤ 5 visible):** `n` "Cells per side $n$" 4 … 16, step 2, default 8 · `field` chips divergent / shear / cavity first step ·
  `stage` 0 … 3 (the transport) · `solver` chips direct / SOR (with `sweeps` 1 … 200, optional) · `mode` chips projection / checkerboard.
- **Transport:** `stage` 0 → 3, rate 0.6 stage/s, `end: 'hold'`, hold 2 s; end card "max∣∇·u∣: … → ~1e-16; the correction's curl:
  ~1e-16" (live).
- **Presets:** "8² divergent field" {n: 8, field: 'divergent'} · "16² shear" {n: 16, field: 'shear'} · "cavity, first step" {n: 16,
  field: 'cavity'} · "SOR, only 20 sweeps" {solver: 'sor', sweeps: 20} · "checkerboard" {mode: 'checkerboard', n: 8} · "3-cell pipe" (a
  1 × 3 special case that reproduces C12's tiny example: p = 0, 1, 2).
- **Status:** stage-dependent: "① uⁿ is divergence-free" · "② predictor: max∣∇·u*∣ = … — not allowed" · "③ Poisson: Σ rhs = 3e-16 ✓
  compatible; p pinned to zero mean" · "④ corrected: max∣∇·u∣ = 4e-16 ✅" · SOR: "⚠️ SOR 20 sweeps: residual … — divergence left …" ·
  checkerboard: "🏁 collocated gradient of the zigzag: 0 everywhere (invisible); staggered: ±2/Δx".
- **Readouts:** "max∣∇·u∣" · "Σ rhs" · "p range" · "curl of correction" · "clicked cell".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** (stages) · **modes** (projection / checkerboard) ·
  **inspector** (click a cell: "(u_{i+½} − u_{i−½})/Δx + (v_{j+½} − v_{j−½})/Δy = (a − b)/0.125 + (c − d)/0.125 = …" with the four face values;
  at a wall cell the wall term shows "(fixed, not corrected)") · **presets** · **status**.
- **Explain:**
  0. *What the views show* — "**Staggered grid**: orange colours are pressures at cell centres, blue arrows velocities on faces, a rose
     border means fluid is piling up in (or draining from) that cell. **Pressure and correction**: the orange arrows are the push −Δt∇p
     that the projection adds. **Divergence left** (hidden on phones): how much net outflow remains after each stage."
  1. *The grid* — "n = **8** cells per side, Δx = Δy = **0.125**; arrays p **8 × 8**, u **8 × 9**, v **9 × 8**; the wall-normal faces (left,
     right columns of u; top, bottom rows of v) are fixed."
  2. *Divergence of the predicted field* — "(10.123) per cell: max∣∇·u*∣ = **…**; in your clicked cell **…**" boxed.
  3. *The Poisson equation* — "(10.124): ∇²_d p = ∇·u*/Δt with Δt = **1**; at a wall cell the wall face's velocity is not corrected, so
     the outside pressure drops out — no pressure boundary condition; Σ rhs = **3e-16** (compatible: net flow through the walls is zero);
     one value pinned, then the mean removed."
  4. *The correction* — "(10.121): u_{i+½} ← u*_{i+½} − (Δt/Δx)(p_{i+1} − p_i): on the face you clicked **…** − (1/0.125)(**…**) = **…**."
  5. *After the projection* — "max∣∇·u∣ = **4e-16** (round-off); the correction −Δt∇p has discrete curl **3e-16**: it changed only the
     divergent part (D19)."
  6. *Checkerboard mode* (or the hint "switch the mode to see why the grid is staggered") — "p = (−1)^{i+j}: collocated (10.125):
     (p_{i+1} − p_{i−1})/(2Δx) = (±1 ∓ ... ) = **0** everywhere; staggered (10.126): (p_{i+1} − p_i)/Δx = **±16** (= 2/Δx). Null spaces:
     collocated **≥ 4** patterns invisible, staggered **1** (the constant)."
  7. *Reading the current setting* — direct: "The projection is exact to round-off in one solve: incompressibility is enforced, not
     approximated." · SOR not converged: "A partly converged Poisson solve leaves divergence behind — real codes iterate to a tolerance
     (or use multigrid); the cost of that solve is why it dominates a MAC step." · checkerboard: "On a collocated grid a zigzag pressure
     is invisible to the momentum equation and can grow unchecked; the staggered grid (Arakawa's C-grid in ocean models) cannot miss it."
- **Derivation tab:** **D19** (9 steps) `view: 'press'`; goal `set` {field: 'divergent', stage: 1}; step 5 `set` {stage: 2} `watch` "the
  Poisson equation appears"; step 8 `set` {stage: 3} `live` "curl of the correction = 3e-16". **D20** (12 steps) `view: 'grid'`; step 3
  `set` {stage: 1}; step 9 `watch` "click a wall cell: the wall face is not corrected"; step 11 `live` "Σ rhs = 3e-16". **D21** (7 steps)
  `view: 'press'`, `set` {mode: 'checkerboard'}; step 4 `live` "collocated: (1 − 1)/(2 × 0.125) = 0". Interpret: `s => "After the
  projection the largest net outflow of any cell is " + divAfter + " — incompressible to round-off."`.
- **Code:**
  ```python
  g = ch10.MacGrid({{n}}, {{n}})                   # p {{n}}x{{n}}, u {{n}}x{{n1}}, v {{n1}}x{{n}}
  st = ch10.projection_stages(us, vs, g, dt=1.0)   # the four stages
  print(abs(st["div_before"]).max())               # {{db}}  (predictor)
  print(st["rhs_sum"])                             # {{rs}}  compatible
  print(abs(st["div_after"]).max())                # {{da}}  after (10.121)-(10.122)
  print(abs(st["curl_correction"]).max())          # {{cc}}  curl-free push
  ```
- **Walkthrough (7 steps):** 1. "Pressure's job" — "Incompressible fluid cannot pile up. What does the pressure do about it? Step
  through one time step." `set` 8² preset · 2. "Predict" — "Move the flow ignoring pressure (10.116). Rose borders: cells that now gain or
  lose fluid." `set` {stage: 1} · 3. "Poisson" — "Find p with ∇²p = ∇·u*/Δt (10.124). Walls need no pressure condition." `set` {stage: 2},
  `derive: {id: 'D20', step: 9}` · 4. "Correct" — "Push with −Δt∇p (10.121): every border turns white." `set` {stage: 3}, `inspect: true` ·
  5. "Only a gradient" — "The push has no curl: vorticity is untouched, only divergence removed." `derive: {id: 'D19', step: 8}` · 6. "Why
  stagger?" — "Checkerboard mode: on a collocated grid the zigzag's gradient is zero — invisible." `set` checkerboard preset, `derive: {id:
  'D21', step: 4}` · 7. "Your turn" — "Predict how many SOR sweeps make the divergence < 10⁻⁶ on 16², then try." `controls: ['solver',
  'sweeps', 'n']`.
- **Equations:** `split` "MAC split" ref 'Eq. (10.115)' $\mathbf A_1=\begin{pmatrix}(\mathbf u\cdot\nabla)\mathbf u-\frac1{Re}\nabla^2\mathbf
  u\\\mathbf 0\end{pmatrix},\ \mathbf A_2=\begin{pmatrix}\nabla p\\\nabla\cdot\mathbf u\end{pmatrix}$ · `proj` "Projection" ref 'Eq.
  (10.117)–(10.118)' $\frac{\mathbf u^{n+1}-\mathbf u^{n+1/2}}{\Delta t}+\nabla p^{n+1}=\mathbf 0,\ \nabla\cdot\mathbf u^{n+1}=0$ · `corr`
  "Face correction" ref 'Eq. (10.121)' $u^{n+1}_{i+1/2,j}=u^{n+1/2}_{i+1/2,j}-\frac{\Delta t}{\Delta x}(p^{n+1}_{i+1,j}-p^{n+1}_{i,j})$ · `cont`
  "Cell continuity" ref 'Eq. (10.123)' $\frac{u_{i+1/2,j}-u_{i-1/2,j}}{\Delta x}+\frac{v_{i,j+1/2}-v_{i,j-1/2}}{\Delta y}=0$ live · `pois`
  "Discrete Poisson" ref 'Eq. (10.124)' $\nabla^2_dp^{n+1}_{i,j}=\frac1{\Delta t}\big(\frac{u_{i+1/2,j}-u_{i-1/2,j}}{\Delta x}+\frac{v_{i,j+1/2}-
  v_{i,j-1/2}}{\Delta y}\big)^{n+1/2}$ · `grads` "Collocated vs staggered" ref 'Eq. (10.125)–(10.126)' $\frac{p_{i+1,j}-p_{i-1,j}}{2\Delta x}$
  vs $\frac{p_{i+1,j}-p_{i,j}}{\Delta x}$ · symbols p, u, v, Δt, Δx.
- **Check yourself:** (1) "Click a cell touching the left wall before and after stage 4. Which face never changes?" — "The wall face:
  its velocity is prescribed and is not corrected — which is why no p outside the wall is needed." · (2) "Why must Σ rhs be zero?" — "The
  Neumann Poisson matrix annihilates constants; a solution exists only if the net divergence (net flow through the walls) is zero." ·
  (3) "Does the correction change the vorticity?" — "No: it is a gradient; its discrete curl is 10⁻¹⁶." · (4) "In checkerboard mode, what
  does the collocated momentum equation feel?" — "Nothing: the centred gradient of (−1)^{i+j} is zero at every node." `set {mode:
  'checkerboard'}`.
- **Selftest parity rows:** `{name: 'div before 8', js: summary(8, 'divergent').div_before_max, py:
  'ch10.mac_projection_summary(8, "divergent")["div_before_max"]', rtol: 1e-12}` · `{name: 'p max 8', js: summary(8, 'divergent').p_max, py:
  'ch10.mac_projection_summary(8, "divergent")["p_max"]', rtol: 1e-9}` · `{name: 'p max shear 12', js: summary(12, 'shear').p_max, py:
  'ch10.mac_projection_summary(12, "shear")["p_max"]', rtol: 1e-9}` · `{name: 'div after', js: summary(8, 'divergent').div_after_max, py:
  'ch10.mac_projection_summary(8, "divergent")["div_after_max"]', atol: 1e-12, rtol: 0}` · `{name: 'staggered null', js: NULL.staggered[8],
  py: 'ch10.gradient_null_space(8, 8, "staggered")', rtol: 0}` · `{name: 'collocated checker grad', js: colGrad(checker(4))[0][1][1], py:
  'ch10.collocated_gradient(ch10.checkerboard(4, 4), 1.0, 1.0)[0][1][1]', atol: 1e-15, rtol: 0}`.
- **Fit plan:** 360×640: status, `grid` (58 %) over `press` (42 %), `hist` hidden (max∣∇·u∣ in a readout and in the status); n capped at
  12 on phones (arrows stay ≥ 6 px); 844×390 side by side; desktop rows [1.25, 1].

### E7 · lid_driven_cavity
- **Title:** "How do I know my CFD answer is right?" · **Summary:** "Run a MAC solver on 16²–32² grids (and load cached finer and
  MacCormack runs): the eddy forms, the centreline approaches Ghia's benchmark, and three grids give the observed order and a Richardson
  estimate." · **CORE:** C14, C15 (also C11–C12 as the engine, N71, N78, N111, D23) · **Reference:** `angular_frequency_explorer_1.html`
  (one clock, a "right now" table with the current row lit) and `amplitude_phase_second_order_II_3.html` (numbered Explain).
- **meta:** `viz:order 7` · `viz:sections 10.4 10.5 10.6` · `viz:equations 10.15 10.81 10.116 10.124 10.127 10.128` · `viz:fluidpy
  ch10.cavity ch10.cavity_centreline ch10.primary_vortex_centre ch10.streamfunction ch10.dt_limit ch10.ghia_centreline
  ch10.ghia_vortex_centre ch10.hou_centres ch10.cavity_error_vs_ghia ch10.richardson_three` · `viz:derivations D23`.
- **Physics:** a JS MAC solver identical to `ch10.cavity` for n ≤ 32 (C.10 10.4: Δt = 0.8 × `dt_limit`, lid ghost value, central
  convective differences with the face averages of `MAC.predictor`, corner-pinned Poisson banded LU factorised once, correction);
  `psi` by the cumulative sum of u along y at the cell corners (ψ = 0 on the walls) ↔ `ch10.streamfunction`; `centre(psi)` ↔
  `ch10.primary_vortex_centre`; `centreline(state)` ↔ `ch10.cavity_centreline`; `errGhia(state, Re)` ↔ `ch10.cavity_error_vs_ghia`;
  `rich(f1, f2, f3, r)` ↔ `ch10.richardson_three`; `dtLim(umax, vmax, Re, dx)` ↔ `ch10.dt_limit`. Cached tables in `CAVITY` (above) for
  64², 128² (Re = 100, 400) and MacCormack 64² (Re = 100) — each labelled "ours, cached, computed by fluidpy".
- **Views** (rows [1.25, 1]): 1. `cavity` "The cavity" (row 0, flex 1.1, `equal: true`) — vorticity heatmap (diverging, clipped at ±5),
  ψ contours (blue; the corner eddies teal), the lid arrow, our eddy centre ✕ and Ghia's (black dot), a clock "t = 7.4"; pointer: click
  → the local (u, v, ω). 2. `centre` "u on the vertical centreline" (row 0, flex 1) — u ∈ [−0.4, 1] against y ∈ [0, 1] (y vertical, like
  the cavity): the current solution (bold blue, the 'so far' curve during spin-up), other grids (faint, one colour each), MacCormack
  (purple dashed, cached), Ghia's 17 points (black with white halo). 3. `conv` "Error against grid size" (row 1, `hidePortrait: true`) —
  log–log max deviation from Ghia (at their points) vs h for every grid run so far (dots), the fitted slope printed, and the Richardson
  estimate of u_min vs h as a second small panel (star at h = 0, Ghia's value as a line).
- **Controls (≤ 5 visible):** `Re` chips 100 / 400 · `grid` chips 16 / 24 / 32 (live) / 64 / 128 (cached) / MacCormack 64 (cached) · `t`
  (the transport) · `show` toggles for other grids (optional) · `speed` (optional).
- **Transport:** `t` 0 → 20 (lid transit times), rate 2 per second of wall time for live grids (the solver runs as many steps per frame
  as the budget allows; ≤ 30 ms/frame), `end: 'hold'`; cached grids replay their stored snapshots (every 0.5 time units). End card:
  "steady at t = … (max∣Δu∣/Δt < 1e-6); max deviation from Ghia … (… % of max∣u∣) at 32²" (all live; Ghia's numbers read from the table).
- **Presets:** "Re 100, 16² (live)" · "Re 100, 32² (live)" · "Re 100, 128² (cached)" · "Re 400, 32² (live)" · "MacCormack 64² (cached)" ·
  "three grids → Richardson" (runs 16, 32 live and loads 64 cached, then fills the `conv` view).
- **Status:** "⏳ spinning up: max∣Δu∣/Δt = … at t = …" · "✅ steady: max∣Δu∣/Δt < 1e-6 at t = …" · "📏 within … % of Ghia (max
  deviation … at y = …)" · "📈 observed order p = … from 16/32/64; Richardson u_min = … (Ghia …)" — every number live, every Ghia
  number read from the embedded table (never typed, not even in this storyboard).
- **Readouts:** "t" · "steps" · "max dev. from Ghia" · "eddy centre" · "observed p".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** + end card · **presets** · **status** · **notes**
  (a benchmark table: Ghia centre, Hou centre, ours per grid — the current grid's row highlighted; `Viz.work.table`).
- **Explain:**
  0. *What the views show* — "**The cavity**: colours are vorticity (red = counter-clockwise), blue lines are streamlines — one big eddy
     and two small corner ones; ✕ is our eddy centre, the black dot Ghia's. **Centreline**: the horizontal velocity along the vertical
     line through the middle; black dots are Ghia, Ghia & Shin (1982). **Error vs h** (hidden on phones): how the deviation shrinks as the
     grid is refined."
  1. *Grid and time step* — "n = **32**, Δx = **0.03125**; (10.127): Δt ≤ 2/(Re·1²) = **0.02**; (10.128): Δt ≤ ReΔx²/4 = **0.0244**; we use
     0.8 × the smaller = **0.016** ⇒ **1250** steps to t = 20" boxed.
  2. *One step* — "predictor (10.116) → Poisson (10.124) (one back-substitution of the factorised matrix) → correction (10.121)–(10.122);
     after each step max∣∇·u∣ ≈ **1e-15**."
  3. *Is it steady?* — "max∣Δu∣/Δt = **…** (live); steady below 10⁻⁶."
  4. *Against the benchmark* — "at Ghia's 17 points: max deviation **…** at y = **…**, i.e. **…** % of max∣u∣; rms **…**" boxed;
     "eddy centre ours (**…**, **…**) vs Ghia (**…**) vs Hou (**…**) — within one cell" (table, current row lit).
  5. *Three grids* (D23) — "u_min on 16/32/64: **…**, **…**, **…**; ratio of differences = **…** ⇒ p = ln(ratio)/ln 2 = **…**; Richardson f₀ ≈ f₁ +
     (f₁ − f₂)/(2^p − 1) = **…** vs Ghia **…**."
  6. *Reading the current setting* — spinning up: "A benchmark is a steady state: wait for the transient to die (about 10–15 lid transits
     at Re = 100)." · coarse steady: "A 16² grid already has the right pattern but the minimum of u is too shallow: numerical diffusion
     of the coarse grid." · fine steady: "Within 1–2 % of an independent, published solution — evidence that the code solves the intended
     equations (verification, not proof)." · Re 400: "Faster flow, thinner layers: the same grid is less converged; the eddy centre moves
     toward the middle." · MacCormack: "The weakly compressible solver agrees with the projection solver to O(Ma²) + O(Δx²)."
- **Derivation tab:** **D23** (8 steps) `view: 'conv'` (on phones `centre`); goal `set` the "three grids" preset; step 5 `live` "ratio =
  (f₃ − f₂)/(f₂ − f₁) = … ⇒ p = …"; step 8 `live` "f₀ ≈ …". Interpret: `s => "On " + n + "² the centreline is within " + pct + " % of
  Ghia; three grids say p = " + p + "."`.
- **Code:**
  ```python
  st = ch10.cavity(Re={{Re}}, n={{n}})              # MAC projection to steady state
  cl = ch10.cavity_centreline(st)                   # u(0.5, y)
  err = ch10.cavity_error_vs_ghia(st, {{Re}})       # at Ghia's 17 points
  print(st["steps"], err["max_dev"], err["rel_max"])   # {{steps}}, {{dev}}, {{rel}}
  f = [ch10.cavity_centreline(ch10.cavity(Re={{Re}}, n=m))["u"].min() for m in (16, 32, 64)]
  print(ch10.richardson_three(f[2], f[1], f[0]))    # p = {{p}}, f_ext = {{fext}}
  ```
- **Walkthrough (7 steps):** 1. "Is it right?" — "A computed flow has no answer key. How do we trust it? Start the lid." `set` Re 100
  16², `play: true` · 2. "The engine" — "Each step: predict, solve for p, correct (10.116)–(10.124). Divergence stays ~10⁻¹⁵."
  `readouts: ['steps']` · 3. "Spin-up" — "The eddy forms in ~10 lid transits; the benchmark is the steady state." `play: true` · 4. "The
  benchmark" — "Black dots: Ghia et al. (1982). 16² misses the minimum; try 32²." `set` {grid: 32} · 5. "Refine" — "Three grids: the error
  falls like h^p. Read p." `set` three-grids preset, `derive: {id: 'D23', step: 5}` · 6. "Extrapolate" — "Richardson: f₀ ≈ f₁ + (f₁ − f₂)/
  (2^p − 1). Compare with Ghia." `derive: {id: 'D23', step: 8}`, `notes: true` · 7. "Your turn" — "Predict: at Re = 400 on 32², larger or
  smaller error? Run it." `controls: ['Re', 'grid']`.
- **Equations:** `ns` "Dimensionless NS" ref 'Eq. (10.81)' $\partial_t\mathbf u+(\mathbf u\cdot\nabla)\mathbf u=-\nabla p+\frac1{Re}\nabla^2\mathbf u$
  · `pred` "Predictor" ref 'Eq. (10.116)' $\frac{\mathbf u^{n+1/2}-\mathbf u^n}{\Delta t}+(\mathbf u^n\cdot\nabla)\mathbf u^n-\frac1{Re}\nabla^2\mathbf
  u^n=\mathbf g^{n+1}$ · `pois` "Pressure" ref 'Eq. (10.124)' (as in E6) · `dt` "Stability" ref 'Eq. (10.127)–(10.128)' $\frac12(u^2+v^2)\Delta
  t\,Re\le1,\ \frac{4\Delta t}{Re\,\Delta x^2}\le1$ live · `err` "Error model" ref 'Eq. (10.15)' $\lVert e\rVert\le K_e\Delta x^a$ · `rich`
  "Richardson (ours, D23)" $f_0\approx f_1+\frac{f_1-f_2}{r^p-1}$ live · symbols Re, n, Δx, Δt, ψ, ω, p (order).
- **Check yourself:** (1) "Why does the 16² solution miss the depth of the u-minimum?" — "Coarse-grid numerical diffusion smooths the
  profile; the minimum sits in a thin region." · (2) "Halving Δx at fixed Re: by how much does the MAC time step shrink?" — "By 4 once the
  diffusive limit (10.128) dominates (Δt ∝ Δx²)." · (3) "What does an observed order far from 2 tell you?" — "Either a bug or grids outside
  the asymptotic range (differences not shrinking by a constant factor)." · (4) "Is matching Ghia verification or validation?" —
  "Verification: another accurate solution of the same equations, not an experiment." 
- **Selftest parity rows:** `{name: 'dt limit 32', js: dtLim(1, 0, 100, 1/32), py: 'ch10.dt_limit(1.0, 0.0, 100.0, 1/32, 1.0)', rtol: 1e-12}` ·
  `{name: 'cavity 16 u_min t=5', js: runCavity(100, 16, 5).umin, py: 'ch10.cavity_centreline(ch10.cavity(Re=100.0, n=16, t_end=5.0,
  tol_steady=0.0, cache=False))["u"].min()', rtol: 1e-8}` · `{name: 'cavity 16 psi_min t=5', js: runCavity(100, 16, 5).psimin, py:
  'ch10.primary_vortex_centre(ch10.cavity(Re=100.0, n=16, t_end=5.0, tol_steady=0.0, cache=False)["psi"], ch10.MacGrid(16, 16))["psi_min"]',
  rtol: 1e-8}` · `{name: 'richardson', js: rich(1.000625, 1.0025, 1.01, 2).f_ext, py: 'ch10.richardson_three(1.000625, 1.0025, 1.01, 2.0)["f_ext"]',
  rtol: 1e-12}` · `{name: 'ghia table', js: CAVITY.ghia[100].u[8], py: 'ch10.ghia_centreline(100)["u"][8]', rtol: 1e-6}` · `{name: 'cached 64
  u_min', js: CAVITY.mac[100][64].umin, py: 'ch10.cavity_centreline(ch10.cavity(Re=100.0, n=64))["u"].min()', rtol: 1e-4}`.
- **Fit plan:** 360×640: status, `cavity` (square, 50 %) over `centre` (50 %), `conv` hidden (p and the Richardson value in the status);
  grid chips wrap; 844×390: `cavity` | `centre`; desktop rows [1.25, 1]. Live runs are throttled on phones (n ≤ 24).

### E8 · mixed_fe_lbb
- **Title:** "Why can't velocity and pressure use the same elements?" · **Summary:** "Switch between P1–P1 and P2–P1 and refine the mesh:
  count velocity and pressure unknowns, watch the Stokes-cavity pressure turn from noise to smooth, and follow the inf–sup constant." ·
  **CORE:** C13 (also N75, N76, N94, N100, N105, N108, R12) · **Reference:** `angular_frequency_explorer_1.html` (modes + presets) and
  `amplitude_phase_second_order_II_3.html` (numbered Explain).
- **meta:** `viz:order 8` · `viz:sections 10.4 10.5` · `viz:equations 10.134 10.135 10.137 10.185 10.187 10.198` · `viz:fluidpy
  ch10.structured_square_mesh ch10.p2_shape ch10.p1_shape ch10.stokes_cavity ch10.infsup_constant ch10.p2_node_counts ch10.tri_quad_7pt` ·
  `viz:derivations D22`.
- **Physics:** `mesh(n, pair)` exactly as C.10 10.6 ↔ `ch10.structured_square_mesh`; `phi2(xi, eta)`, `phi1` ↔ `ch10.p2_shape`,
  `ch10.p1_shape`; `quad7()` ↔ `ch10.tri_quad_7pt`; `stokes(n, pair)` dense assembly (Re → 0: A = vector Laplacian, B from ∫q∇·v) with
  the 7-point rule, all-Dirichlet velocity (lid u = 1), pressure pinned by zero mean (a Lagrange row), dense LU (≤ 400 unknowns: P2P1
  n ≤ 6 → 291) ↔ `ch10.stokes_cavity`; `infsup(n, pair)` BᵀA⁻¹B and M_p on mean-zero pressures, Cholesky + cyclic Jacobi (≤ 49 × 49) ↔
  `ch10.infsup_constant`; larger n from `INFSUP`.
- **Views** (rows [1.15, 1]): 1. `mesh` "Mesh and unknowns" (row 0, flex 1, `equal: true`) — the triangulated unit square; velocity nodes
  (blue dots; boundary ones hollow = prescribed), pressure nodes (orange squares); the counting box "N_u = 98, N_p = 24"; pointer: click a
  triangle → inspector (its six/three nodes and φ values at its centroid); in *element* mode the parent triangle with the six φ_a as
  small contour maps. 2. `press` "Stokes cavity pressure" (row 0, flex 1, `equal: true`) — coloured triangulation of p (orange diverging,
  symmetric range), the lid arrow; a caption with p range and its alternating ('checker') part. 3. `beta` "Inf–sup constant β_h" (row 1,
  `hidePortrait: true`) — log–log β_h vs n (2 … 16) for P1P1 (rose; zeros drawn at the axis floor with "0 (spurious modes)") and P2P1
  (blue), the current n marked.
- **Controls (≤ 5 visible):** `pair` chips P1–P1 / P2–P1 · `n` "Mesh $n\times n$" 2 … 16 (live pressure for n ≤ 6; table beyond) · `count`
  toggle "show counting" · `elem` toggle "one element" (optional).
- **Transport:** `n` 2 → 16, rate 1 per second, `end: 'hold'` (the β view fills from the table; the pressure view updates for n ≤ 6 and
  shows "table only" beyond).
- **Presets:** "P1–P1, n = 2 (counting)" · "P1–P1, n = 4" · "P2–P1, n = 4" · "P2–P1, n = 6" · "refinement sweep" (play the transport) · "one
  element".
- **Status:** "💥 P1–P1, n = 2: 2 velocity unknowns for 8 pressure constraints — 6 spurious modes, β_h = 0: LBB fails" · "⚠️ P1–P1, n = 8:
  β_h = … and falling" · "✅ P2–P1: β_h ≈ … (bounded as n grows) — a stable pair".
- **Readouts:** "N_u" · "N_p" · "spurious modes" · "β_h" · "p range".
- **Depth features:** Explain + Code + Derivation · **modes** (P1–P1 / P2–P1; mesh / element) · **linked views** (3) · **presets** ·
  **status** · **inspector** (click a triangle: "nodes 12, 13, 17 (vertices), 40, 41, 52 (mid-edges); φ at the centroid = (−1/9, −1/9,
  −1/9, 4/9, 4/9, 4/9); J = 2 × area = 0.0625") · **transport**.
- **Explain:**
  0. *What the views show* — "**Mesh**: blue dots carry velocity, orange squares pressure; hollow dots sit on walls where the velocity is
     given. **Pressure**: the computed pressure of a slow (Stokes) flow in a lid-driven box — smooth is right, speckled is spurious.
     **β_h** (hidden on phones): how strongly the worst pressure pattern is felt by the velocity; it must not fall to zero."
  1. *Counting* — "n = **4**: V = (n + 1)² = **25**, E = 3n² + 2n = **56**, T = 2n² = **32** (V − E + T = 1); P1 velocity nodes = V, interior
     (n − 1)² = **9** ⇒ N_u = **18**; P2 nodes V + E = **81**, interior (2n − 1)² = **49** ⇒ N_u = **98**; pressure N_p = V − 1 = **24**" boxed.
  2. *Constraints vs freedom* — "each pressure unknown is one discrete continuity constraint (10.135); if N_p > N_u (P1–P1 at n ≤ 4) at
     least N_p − N_u patterns cannot be felt: **6** at n = 2."
  3. *The inf–sup constant* — "β_h² = λ_min(BᵀA⁻¹B, M_p) (P247) = **…** ⇒ β_h = **…**; spurious modes (λ < 10⁻¹⁰): **…**."
  4. *The pressure you see* — "range **…**; its alternating component **…** — for P1–P1 the noise is as large as the physical pressure; for
     P2–P1 it is **…**."
  5. *One element* (or the hint "switch on 'one element'") — "φ₁ = ζ(2ζ − 1), …, φ₄ = 4ξζ (10.185) with ζ = 1 − ξ − η; ψ₁ = ζ, ψ₂ = ξ, ψ₃ = η
     (10.187); at the centroid the vertex shapes are −1/9 and the mid-edge ones 4/9; the 7-point rule (10.198) integrates degree 5
     exactly."
  6. *Reading the current setting* — P1P1: "Equal order: the velocity space is too poor to feel every pressure pattern — the FE twin of
     the checkerboard (C12). Stabilise (GLS) or enrich the velocity." · P2P1: "Taylor–Hood: quadratic velocities feel every linear
     pressure; β_h stays bounded on refinement — the LBB condition holds, as for the staggered grid."
- **Derivation tab:** **D22** (11 steps) `view: 'mesh'`; goal `set` {pair: 'P2P1', n: 4}; step 6 `watch` "ũ vanishes on the hollow
  (Dirichlet) nodes"; step 11 `watch` "the last line is the continuity row Bᵀ". Interpret: `s => pair + " at n = " + n + ": N_u = " + nu + ",
  N_p = " + np + ", β_h = " + beta + "."`.
- **Code:**
  ```python
  mesh = ch10.structured_square_mesh({{n}}, "{{pair}}")
  print(mesh["counts"])                            # N_u = {{nu}}, N_p = {{np}}
  r = ch10.infsup_constant({{n}}, "{{pair}}")      # generalised eigenvalue
  print(r["beta"], r["n_spurious"])                # {{beta}}, {{nsp}}
  s = ch10.stokes_cavity({{n}}, "{{pair}}")        # slow lid-driven box
  print(s["p_range"], s["p_std_checker"])          # {{prange}}, {{checker}}
  ```
- **Walkthrough (6 steps):** 1. "Same elements?" — "Use the same simple triangles for velocity and pressure. What goes wrong? Look at the
  pressure." `set` P1–P1 n = 4 · 2. "Weak Navier–Stokes" — "The pressure enters as a multiplier of ∇·ũ (10.134) and tests continuity
  (10.135)." `derive: {id: 'D22', step: 9}` · 3. "Count" — "P1–P1 at n = 2: 2 velocity unknowns, 8 pressure constraints. Some patterns
  cannot be felt." `set` counting preset, `readouts: ['nu', 'np']` · 4. "Taylor–Hood" — "Give velocity six nodes per triangle (P2), keep
  pressure on three. The noise vanishes." `set` {pair: 'P2P1'} · 5. "The trend" — "Refine: P1–P1's β_h falls to zero, P2–P1's stays flat. That
  is the LBB condition." `play: true` · 6. "Your turn" — "Click a triangle: predict its φ values at the centroid, then check." `inspect:
  true`.
- **Equations:** `wns` "Weak momentum" ref 'Eq. (10.134)' $\int_\Omega(\partial_t\mathbf u+\mathbf u\cdot\nabla\mathbf u-\mathbf g)\cdot\tilde{\mathbf u}
  +\frac2{Re}\int_\Omega\mathbf D[\mathbf u]:\mathbf D[\tilde{\mathbf u}]-\int_\Omega p\,\nabla\cdot\tilde{\mathbf u}=0$ · `wc` "Weak continuity" ref
  'Eq. (10.135)' $\int_\Omega\tilde p\,\nabla\cdot\mathbf u\,d\Omega=0$ · `kkt` "Saddle point" ref 'Eq. (10.137)' $\begin{pmatrix}\mathbf A&\mathbf
  B\\\mathbf B^T&\mathbf 0\end{pmatrix}\begin{pmatrix}\mathbf u\\\mathbf p\end{pmatrix}=\begin{pmatrix}\mathbf f_u\\\mathbf f_p\end{pmatrix}$ (steady) ·
  `p2` "Quadratic shapes" ref 'Eq. (10.185)' $\phi_1=\zeta(2\zeta-1),\dots,\phi_4=4\xi\zeta$ · `p1` "Linear pressure shapes" ref 'Eq. (10.187)'
  $\psi_1=\zeta,\ \psi_2=\xi,\ \psi_3=\eta$ · `infsup` "Inf–sup (ours)" $\beta_h^2=\lambda_{\min}(\mathbf B^T\mathbf A^{-1}\mathbf B,\mathbf M_p)$ live
  · symbols N_u, N_p, β_h, φ_a, ψ_b.
- **Check yourself:** (1) "For P1–P1 at n = 2, how many spurious pressure modes must there be at least?" — "8 − 2 = 6: more constraints
  than velocity unknowns." `set {pair: 'P1P1', n: 2}` · (2) "Does P1–P1 become stable once N_u > N_p (n ≥ 5)?" — "No: counting is necessary,
  not sufficient; β_h keeps falling with h." `set {pair: 'P1P1', n: 8}` · (3) "What is the FD cousin of this problem?" — "The checkerboard on a
  collocated grid (C12); staggering is the FD cure, P2–P1 the FE cure." · (4) "At the centroid, what are the P2 vertex shape values?" —
  "−1/9 each; the mid-edge ones are 4/9."
- **Selftest parity rows:** `{name: 'phi4', js: phi2(0.2, 0.3)[3], py: 'ch10.p2_shape(0.2, 0.3)[3]', rtol: 1e-12}` · `{name: 'N_u P2 n4', js:
  mesh(4, 'P2P1').counts.N_u, py: 'ch10.structured_square_mesh(4, "P2P1")["counts"]["N_u"]', rtol: 0}` · `{name: 'beta P2P1 n2', js:
  infsup(2, 'P2P1').beta, py: 'ch10.infsup_constant(2, "P2P1")["beta"]', rtol: 1e-6}` · `{name: 'spurious P1P1 n2', js: infsup(2,
  'P1P1').n_spurious, py: 'ch10.infsup_constant(2, "P1P1")["n_spurious"]', rtol: 0}` · `{name: 'Stokes p range', js: stokes(4, 'P2P1').p_range,
  py: 'ch10.stokes_cavity(4, "P2P1")["p_range"]', rtol: 1e-6}` · `{name: 'quad weight sum', js: quad7().W.reduce((a, b) => a + b), py:
  'ch10.tri_quad_7pt()[1].sum()', rtol: 1e-14}`.
- **Fit plan:** 360×640: status, `mesh` | `press` side by side (squares, 48 % each) over the counting readouts; `beta` hidden (β_h in the
  status); 844×390: three views in one row; desktop rows [1.15, 1]. Live pressure solves only for n ≤ 6 (≤ 291 unknowns, < 50 ms).

### B1 · operator_splitting_theta
- **Title:** "Why is a split time step only first order?" · **Summary:** "On a two-variable system whose operators do not commute, compare
  the exact solution with Marchuk–Yanenko splitting and Glowinski's Θ-scheme: the error falls like Δt or Δt², and only Θ = 1 − 1/√2 gives
  the second." · **CORE:** C11 (also N60, N61, N74) · **Reference:** `amplitude_phase_second_order_II_3.html` (linked phase plane + graph).
  Built only if an explainer above fails review (the notebook does not embed it).
- **Physics:** `exact(t)` by a 2 × 2 matrix exponential (closed form via eigen-decomposition) ↔ `ch10.split_linear_system()`; `my(dt, n)`
  two backward-Euler substeps (10.113)–(10.114) ↔ `ch10.marchuk_yanenko`; `theta(dt, n, Θ)` three substeps ↔ `ch10.theta_scheme_linear`;
  `order(method)` ↔ `ch10.splitting_order`; `z2(Θ)` the z² coefficient of R(z) − e^z ↔ `ch10.theta_scheme_amplification_sympy()`.
- **Views:** `phase` (trajectories in the (φ₁, φ₂) plane: exact muted, Marchuk–Yanenko teal, Θ-scheme purple) · `err` (log–log error vs Δt,
  slopes 1 and 2) · `z2` (hidePortrait: the z² coefficient vs Θ with its zero at 1 − 1/√2 marked, the value computed in JS). **Controls:** Θ (0.05–0.45), Δt (log), the
  commutator strength ε (A₂ → A₂ + εC). **Presets:** "Θ = 1 − 1/√2", "Θ = 0.25", "commuting operators (ε = 0)". **Status:** "order 2 (Θ at the
  magic value)" / "order 1". **Depth features:** linked views, presets, status, inspector (click a Δt: the two errors and their ratio).
  **Derivations:** none. **Walkthrough:** 5 steps (problem → split → error ∝ [A₁, A₂] → Θ-scheme → your turn). **Check:** (1) "With ε = 0 is
  Marchuk–Yanenko exact?" — "No: each substep is backward Euler, first order anyway." (2) "Which Θ cancels the Δt² error?" — "1 − 1/√2."
  (3) "Why do splitting errors matter for climate models?" — "Dynamics/physics splitting is first order unless symmetrised." **Selftest:**
  `{name: 'MY order', js: order('marchuk_yanenko'), py: 'ch10.splitting_order("marchuk_yanenko")', rtol: 1e-6}` · `{name: 'theta root', js:
  1 - 1/Math.SQRT2, py: 'ch10.theta_scheme_amplification_sympy()["theta_root"]', rtol: 1e-12}`. **Fit plan:** as E4.

---

## Part D — runtime budget (full run < 5 min on a laptop / Colab CPU)

Chapter 10 is numerics, so the cost is in the solvers, not in closed forms. The heavy items and their budgets: the MAC cavity
(64²: Δt = 0.8 × 0.0061 ⇒ ≈ 4100 steps; each step = a vectorised predictor + one back-substitution of the `splu`-factorised
64²-unknown Poisson matrix (factorised once, ≈ 0.05 s) ≈ 3 ms ⇒ ≈ 12 s; 32² ≈ 2 s; 16² < 0.5 s), the MacCormack cavity (32² at our Ma =
0.08, σ = 0.8: Δt = 1.4e-3 ⇒ 14 100 steps × ≈ 0.4 ms ≈ 6 s; 64² ≈ 28 000 steps × 1.5 ms ≈ 40 s — script and cache only), the block in a channel
(Δx = 1/8: 233 × 33 nodes, Δt = 4.8e-3 ⇒ 6300 steps ≈ 5 s, only when not FAST; finer grids script-only), the FE cylinder (steady Re = 40 on
the coarse mesh ≈ 6 Newton iterations × ≈ 1 s; Re = 1, 10 similar; unsteady Re = 100 script-only, 30–60 min), the mixed-FE inf–sup
eigenproblems (n ≤ 8 in the notebook; the table to n = 16 is precomputed in `reference/ch10/infsup_table.csv`), nine sympy engines (each
`lru_cache`d, 1–3 s), 4 animations and 7 plotly figures. ch09 ran in ≈ 219 s with 14 CORE blocks; ch10 has 15.

| Section | Heaviest cells | Full (cold cache) | Full (warm cache) | FAST (`FLUIDPY_FAST=1`) |
|---|---|---|---|---|
| setup + imports | numpy/scipy/sympy/plotly, five new core modules | 10 s | 10 s | 10 s |
| §10.1 | two stencil numbers | 1 s | 1 s | 1 s |
| §10.2 C01 | 5 stencils × 4 h, 120-point log–log sweep, Vandermonde from scratch, F1 (25 steps) | 4 s | 4 s | 3 s (12 steps) |
| §10.2 C02 | FTCS 100 nodes × 80 steps, BTCS at β = 2.5, u = 0 parity with ch01, one BTCS solve | 2 s | 2 s | 2 s |
| §10.2 C03 | `truncation_error_sympy("ftcs")`, two convergence studies (n ≤ 160), figure | 5 s | 5 s | 3 s (n ≤ 80) |
| §10.2 C04 | 200 random (α, β) × 2001 θ, 60 × 60 Noye brute force, 600-step kick, `lax_demo`, F2 (30), live | 5 s | 5 s | 3 s (30 × 30, 15 steps) |
| §10.2 C05 | `advect_periodic` × 6 (N = 100), rod series, **A1 frames 40** (lax_demo 2000 steps × 3), F3 (24) | 9 s | 9 s | 5 s (24 frames, 12 steps) |
| §10.3 C06–C08 | `weak_form_sympy`, `weak_to_strong_sympy`, `galerkin_equations_sympy(3)`, FE runs n ≤ 32, transient n = 20, F4 (9 steps), assembly figure | 9 s | 9 s | 7 s |
| §10.4 C09 | steady FD/FE tiny systems, wiggle scan, stretched grid, F5 (30), live | 3 s | 3 s | 2 s |
| §10.4 C10 | `maccormack_linear_sympy`, advection orders (N ≤ 400), artificial compressibility (ny = 16), `compressible_ns_sympy`, **A2 video 60** | 13 s | 13 s | 8 s (N ≤ 200, 30 frames) |
| §10.4 C11 | `projection_sympy`, one 16² projection, splitting orders, `theta_scheme_amplification_sympy`, **A3 frames 32** | 9 s | 9 s | 6 s (16 frames) |
| §10.4 C12 | Poisson matrix 8² (rank), stage timing 50 steps at 64², `ftcs2d_max_amplification` × 2, Taylor–Green 16/32/64 (t = 0.5), checkerboard SVD | 7 s | 7 s | 4 s (TG 16/32, timing at 32²) |
| §10.4 C13 | `infsup_constant` n = 2, 4, 8 × 2 pairs, `stokes_cavity(6)` × 2, spy, P2/quadrature from scratch, F7 (from CSV) | 6 s | 6 s | 4 s (n ≤ 4) |
| §10.5 C14 | `MAC.cavity` 64² (32² FAST), `MCK.cavity_maccormack` 32², Ghia comparison, **A4 video 60** (cached snapshots) | 35 s | 12 s | 12 s (32², 30 frames) |
| §10.5 C15 | cavity 16/32/64 (64 reused), `grid_convergence_index`, block Re = 20 at Δx = 1/8 (not FAST), float32/64 drift (5000 steps on 16²), F6 (cached grids or live 16/24/32) | 18 s | 10 s | 6 s (12/24/48 grids, cached; drift 2000 steps; no block run) |
| §10.5 C13 (continued) | `cylinder_channel_mesh`, `cylinder_steady` Re = 1, 10, 40 (cached; cold coarse ≈ 6 s each), force CSV, curved-element figure | 22 s | 3 s | 8 s (Re = 40 only) |
| §10.6 | time-step refinement at 32² (2 runs) | 4 s | 4 s | 4 s |
| explainers (8 × `show_viz`) | read the HTML files | 1 s | 1 s | 1 s |
| **Total** | | **≈ 163 s** | **≈ 118 s** | **≈ 89 s** |

**FAST plan.** Every size-dependent choice is written `a if not FAST else b` in the cell: convergence studies n ≤ 160 → 80; brute-force
Noye grid 60² → 30²; MacCormack advection N ≤ 400 → 200; MAC cavity 64² → 32² (C14) and the C15 three-grid study 16/32/64 → 12/24/48
(r = 2 kept; 48² ≈ 5 s, cached); Taylor–Green 16/32/64 → 16/32 (+ a stated two-grid order); block run skipped (the cached C_D figure is used
if present, else a sentence); cylinder Re = 1, 10 skipped (Re = 40 only, coarse); videos 60 → 30 frames, frame players 40/32 → 24/16; plotly
sliders ≤ 30 → ≤ 15 steps (≤ 4 traces × ≤ 400 points, < 250 kB each); sympy engines unchanged (cached). **Cached arrays / results**
(`outputs/ch10/*.npz`, keyed by the parameters, P252): MAC cavity at every (Re, n) (C14, C15, §10.6, A4 and F6 share them), MacCormack
cavity 32²/64², block C_D histories, cylinder steady fields and residual histories, the A1 lax-demo histories; `reference/ch10/`
(committed, public, ours or cited): Ghia Table I, Hou centres, `infsup_table.csv`, `cylinder_fe_re100_forces.csv`. On a fresh Colab
runtime nothing is cached: the notebook computes MAC 64² (≈ 12 s), MacCormack 32² (≈ 5 s) and the coarse cylinder Re = 40 (≈ 6 s) and
says "computed now (no cache)"; the 128² cavity, MacCormack 64², block convergence and the unsteady cylinder appear only from
`reference/ch10/` or with the note "run `scripts/ch10_cavity.py` / `ch10_block.py` to add this". Sympy cells never `simplify` expressions
with generic functions beyond the final residual. **Outputs:** 2 videos (A2, A4) + 2 frame players (A1, A3) at dpi 80 (each < 3 MB), 7
plotly figures, ≈ 40 static figures — the page stays under 15 MB.

---

## Part E — prerequisite ledger
Every concept, symbol, maths tool and Python function or idiom the notebook or its explainers use, with where it is explained.
"primer (in Cxx)" = a 📎 primer placed in that block before first use (the Concept text is the primer term, used verbatim in
`nb.primer`); "knowledge/primers.md: <term> (chNN Pnn) — reminder" = a one-line reminder naming the earlier primer; a CORE/RECAP id alone
= taught there; "Cxx (Nnn)" = the NOTE placed in that block; "Cxx (D0n)" = the derivation where it is used; "Cxx (gloss …)" = one
sentence where it is used. Earlier chapters' material that is neither a primer nor a ch10 RECAP is named by section ("Ch. 2 §2.9"). New
primers P221–P253 (33) in first-use order: big-O and order (C01), round-off (C01), ghost node (C02), error norms (C03), half-angle
identities, sign of a line (C04), domain of dependence, well-posed, Péclet (C05), test functions and spaces, fundamental lemma (C06),
bilinear form, arbitrary coefficients (C07), affine map, scatter-add (C08), linear recurrence (C09), conservation form, predictor–corrector,
Lagrange multiplier (C10), Helmholtz–Hodge, operator splitting (C11), half-index arrays, singular systems, sparse diags/kron/splu, SVD null
space (C12), saddle-point matrix, GMRES, generalised eigenproblem, barycentric coordinates, Jacobian determinant (C13), verification vs
validation, reference data + np.interp, caching (C14).

| Concept | First used in | Explained by |
|---|---|---|
| computational fluid dynamics (what it predicts) | C01 | C01 (N01) |
| conservation laws of mass and momentum | C01 | Ch. 4 §4.2, §4.4 (named in C01, N01) |
| four error sources (discretisation, data, conditions, models) | C01 | C01 (N02) |
| Moore's law and why verification still matters | C01 | C01 (N03) |
| method families (FD, FE, FV, spectral) | C01 | C01 (N04) |
| uniform space–time grid, T_i^n notation | C01 | R01 |
| Taylor series at the grid neighbours (10.4)–(10.5) | C01 | R02 |
| first-order Taylor expansion | C01 | knowledge/primers.md: first-order Taylor expansion (ch01 P26) — reminder |
| multivariable Taylor expansion | C03 | knowledge/primers.md: multivariable first-order Taylor expansion (ch03 P98) — reminder |
| centred second derivative (10.7) | C01 | R03 |
| finite differences and FTCS for diffusion | C01 | knowledge/primers.md: finite differences (ch01 P21) — reminder |
| big-O notation and the order of accuracy | C01 | primer (in C01) |
| limits and orders of smallness | C01 | knowledge/primers.md: limits and orders of smallness (ch02 P68) — reminder |
| forward, backward and centred stencils (10.6) | C01 | C01 (D01) |
| Taylor-matching weights (Vandermonde system) | C01 | C01 (D01, from-scratch cell) |
| np.linalg.solve | C01 | knowledge/primers.md: np.linalg.solve (ch01 P57) — reminder |
| fractions and exact rational weights (sympy Rational) | C01 | knowledge/primers.md: fractions.Fraction (ch01 P60) — reminder |
| power laws and log–log slopes, observed order | C01 | knowledge/primers.md: power laws and log–log plots (ch01 P13) — reminder |
| tools.convergence.observed_order, pairwise_orders | C01 | C01 (gloss: slope of log error vs log h, P13) |
| floating-point round-off and machine epsilon | C01 | primer (in C01) |
| float32 vs float64 | C01 | C01 (round-off note and figure) |
| time differences, leapfrog (10.8) | C01 | C01 (N07) |
| model problem (10.1): convection–diffusion of a scalar | C01 | C01 (N05) |
| travelling, spreading Gaussian as an exact solution | C01 | C01 (N05) |
| diffusion spreading s² = s₀² + 2Dt | C01 | knowledge/primers.md: diffusivity and the diffusion time L²/ν (ch08 P185) — reminder |
| numpy arrays and slicing of neighbours | C01 | knowledge/primers.md: numpy broadcasting (ch02 P77) — reminder |
| np.linspace and np.logspace | C01 | knowledge/primers.md: np.linspace and np.logspace (ch01 P06) — reminder |
| matplotlib figures, log axes | C01 | knowledge/primers.md: matplotlib figures (ch01 P01) — reminder |
| f-strings | C01 | knowledge/primers.md: f-strings (ch01 P04) — reminder |
| Python dictionaries (returned results) | C01 | knowledge/primers.md: Python dictionaries (ch01 P23) — reminder |
| functions as arguments and lambda | C01 | knowledge/primers.md: functions as arguments and lambda (ch01 P29) — reminder |
| assert np.allclose | C01 | knowledge/primers.md: assert np.allclose (ch01 P15) — reminder |
| slider_figure | C01 | knowledge/primers.md: slider_figure (ch01 P17) — reminder |
| show_viz | C01 | knowledge/primers.md: show_viz (ch01 P18) — reminder |
| partial derivative | C01 | knowledge/primers.md: partial derivative (ch01 P25) — reminder |
| explicit FTCS scheme (10.9)–(10.11), α and β | C02 | C02 (D02) |
| explicit stepping and its limit | C02 | knowledge/primers.md: explicit stepping (ch01 P30) — reminder |
| Dirichlet and Neumann conditions (10.2), initial condition (10.3) | C02 | C02 (N06) |
| boundary conditions | C02 | knowledge/primers.md: boundary conditions (ch01 P20) — reminder |
| ghost node for a Neumann boundary | C02 | primer (in C02) |
| periodic neighbours with np.roll | C02 | C02 (gloss: np.roll shifts an array cyclically) |
| implicit BTCS scheme (10.12)–(10.13), tridiagonal system | C02 | C02 (N08) |
| Crank–Nicolson and solve_banded | C02 | knowledge/primers.md: Crank–Nicolson with `scipy.linalg.solve_banded` (ch08 P193) — reminder |
| implicit time stepping | C02 | knowledge/primers.md: implicit time stepping with Picard iteration (ch08 P192) — reminder |
| explicit vs implicit algorithms | C02 | C02 (N09) |
| conservation of the sum under FTCS (weights add to 1) | C02 | C02 (tiny example) |
| consistency and the truncation error (10.16)–(10.17) | C03 | C03 (D03) |
| modified equation (the scheme solves the PDE + E) | C03 | C03 (D03) |
| convergence, error (10.14), rates a, b (10.15) | C03 | C03 (N10) |
| error constant K_e (not the stiffness K) | C03 | C03 (N10) |
| norms of an error array: rms, max, L1 | C03 | primer (in C03) |
| sympy symbols, subs, series, removeO | C03 | knowledge/primers.md: sympy expand, series, removeO, collect and subs (ch04 P117) — reminder |
| sympy basics | C03 | knowledge/primers.md: sympy (ch01 P40) — reminder |
| β = 1/6 special FTCS (fourth order for diffusion) | C03 | C03 (D03 step 11) |
| stability (round-off must not grow) | C04 | C04 (N11) |
| disturbance and error equation (10.18)–(10.19) | C04 | C04 (N12, D04) |
| linearity and superposition | C04 | C04 (D04) |
| Fourier modes of a grid error, the book's e^{iπkx} convention | C04 | R04 |
| Fourier modes and the FFT Poisson solver | C04 | knowledge/primers.md: Fourier modes and the FFT Poisson solver (ch05 P142) — reminder |
| Euler's formula, imaginary unit | C04 | knowledge/primers.md: square root of a negative number (ch01 P45) — reminder |
| complex numbers in numpy (1j, np.abs, np.angle) | C04 | knowledge/primers.md: the complex plane in numpy (ch06 P153) — reminder |
| complex conjugate (G(−θ)) | C04 | knowledge/primers.md: complex conjugate (ch02 P81) — reminder |
| amplification factor G (10.24) | C04 | C04 (D04) |
| stability condition ∣G∣ ≤ 1 (10.25) | C04 | C04 (N13, D04) |
| half-angle identities | C04 | primer (in C04) |
| ∣G∣² of FTCS (10.26) | C04 | C04 (N14, D05) |
| sign of a linear function on an interval | C04 | primer (in C04) |
| Noye's region (10.27) | C04 | C04 (N15, D06) |
| pure diffusion limit β ≤ ½ (10.28) | C04 | R05 |
| pure convection never stable for FTCS | C04 | C04 (D06) |
| BTCS unconditionally stable | C04 | C04 (N17, D07) |
| random generator with a seed | C04 | knowledge/primers.md: np.random.default_rng (ch01 P10) — reminder |
| live widgets | C04 | knowledge/primers.md: live widgets (ch01 P47) — reminder |
| upwind scheme (10.29), Courant number C = 2α | C05 | C05 (N16, D08) |
| CFL condition (10.30) | C05 | C05 (D08) |
| sign of u and the upwind side | C05 | C05 (N16, D08 step 9) |
| domain of dependence and characteristics | C05 | primer (in C05) |
| first-order wave equation and characteristics | C05 | knowledge/primers.md: first-order wave equation and characteristics (ch07 P174) — reminder |
| numerical diffusivity of upwind advection ∣u∣Δx(1 − C)/2 | C05 | C05 (code note, D17 twin) |
| climate CFL: gravity-wave speed √(gH), model time steps | C05 | C05 (climate note) |
| shallow-water wave speed √(gH) | C05 | Ch. 7 §7.2 (recalled in C05) |
| well-posed problem | C05 | primer (in C05) |
| Lax equivalence theorem | C05 | C05 (N18) |
| heated rod: Fourier sine series solution (10.199) | C05 | C05 (N113) |
| separation of variables for a PDE | C05 | knowledge/primers.md: separation of variables for a PDE (ch07 P167) — reminder |
| Péclet number (global R and cell R_cell) | C05 | primer (in C05) |
| cell Péclet condition (10.31) | C05 | C05 (N19) |
| animate and show_animation | C05 | knowledge/primers.md: animate and show_animation (ch01 P16) — reminder |
| strong (classical) form | C06 | C06 (N20) |
| test functions and the spaces H¹, S and V | C06 | primer (in C06) |
| trial space 𝒮 and test space V (10.32)–(10.33) | C06 | C06 (N21) |
| product rule | C06 | knowledge/primers.md: product rule for differentials (ch01 P38) — reminder |
| integration by parts | C06 | knowledge/primers.md: integration by parts (ch09 P218a) — reminder |
| fundamental theorem of calculus | C06 | knowledge/primers.md: fundamental theorem of calculus (ch02 P84) — reminder |
| weak form (10.34)–(10.36) | C06 | C06 (D09) |
| fundamental lemma of the calculus of variations | C06 | primer (in C06) |
| weak ⇒ strong, essential vs natural conditions (10.37)–(10.38) | C06 | C06 (N22, D10) |
| trapezoid rule np.trapezoid | C06 | knowledge/primers.md: `scipy.integrate.simpson` and `np.trapezoid` (ch09 P203) — reminder |
| discrete weak problem, lifting, Galerkin form (10.39)–(10.42) | C07 | C07 (N23, D11) |
| Galerkin vs Petrov–Galerkin | C07 | C07 (N24) |
| bilinear form a(w, v) | C07 | primer (in C07) |
| basis, shape functions N_A, N₀ (10.43)–(10.49) | C07 | C07 (N25, N26, N27) |
| linear combinations of basis functions | C07 | Ch. 2 §2.1 (recalled in C07) |
| arbitrary coefficients: every bracket is zero | C07 | primer (in C07) |
| n ODEs (10.50)–(10.53) | C07 | C07 (N28, D11) |
| mass matrix, stiffness matrix, force vector, M ḋ + K d = F (10.54)–(10.58) | C07 | C07 (D11) |
| method of lines vs FD in time | C07 | C07 (N29) |
| solve_ivp | C07 | knowledge/primers.md: scipy.integrate.solve_ivp (ch01 P31) — reminder |
| hat functions (10.59)–(10.61), compact support, δ_AB | C07 | C07 (N30) |
| Kronecker delta | C07 | Ch. 2 §2.7 (recalled in C07, N30) |
| nodal values d_A = T_A (10.62) | C07 | C07 (N31) |
| interior FE row (10.63), consistent mass ⅙–⅔–⅙ | C07 | C07 (N32, D12) |
| lumped mass | C07 | C07 (N32) |
| FE = FD on uniform linear meshes (and its caveat) | C07 | C07 (N33) |
| partition of unity | C07 | C07 (N25–N27 code) |
| scipy.sparse matrices and spsolve | C07 | knowledge/primers.md: boolean masks and scipy.sparse (ch06 P161) — reminder |
| element viewpoint, parent element (10.64)–(10.66) | C08 | C08 (N34) |
| affine map to a parent element | C08 | primer (in C08) |
| substitution in an integral | C08 | knowledge/primers.md: substitution in an integral (ch03 P106) — reminder |
| chain rule | C08 | knowledge/primers.md: chain rule (ch01 P49) — reminder |
| shape-function slopes (10.67)–(10.68), the printed index shift | C08 | C08 (N35, D13) |
| global = sum of element matrices (10.69)–(10.70) | C08 | C08 (N36) |
| scatter-add assembly (np.add.at, COO duplicates) | C08 | primer (in C08) |
| element integrals (10.71)–(10.73) and the nonzero entries | C08 | C08 (N37) |
| element matrices and force (10.74)–(10.77) | C08 | C08 (D14) |
| local → global node map (10.78) | C08 | C08 (N38) |
| Gauss–Legendre quadrature | C08 | knowledge/primers.md: Gauss–Legendre quadrature in 3-D and a smoothed kernel (ch05 P143) — reminder |
| 2-D and 3-D finite elements follow the same steps | C08 | C08 (N39) |
| incompressible Navier–Stokes and continuity (10.79)–(10.80) | C09 | R06 |
| dimensionless Navier–Stokes (10.81), Re | C09 | R07 |
| scope of §10.4 (primitive variables, ψ–ω, laminar) | C09 | C09 (N40) |
| the two difficulties of incompressible CFD | C09 | C09 (N42) |
| steady convection–diffusion problem (10.84)–(10.85) | C09 | C09 (N43) |
| linear second-order ODE by an exponential trial | C09 | knowledge/primers.md: linear second-order ODE (ch01 P44) — reminder |
| exact steady solution (10.86), global Péclet R (10.87) | C09 | C09 (N44, D15) |
| large-R form (10.88), layer thickness (10.89) | C09 | C09 (N45, D15) |
| np.expm1 and overflow-safe exponentials | C09 | knowledge/primers.md: np.expm1 and cancellation near zero (ch03 P107) — reminder |
| natural logarithm and exponential (e⁻¹, e⁻²) | C09 | knowledge/primers.md: natural logarithm and exponential (ch01 P36) — reminder |
| centred FD of the steady problem (10.90)–(10.91) | C09 | C09 (N46) |
| linear recurrence with constant coefficients (geometric trial) | C09 | primer (in C09) |
| quadratic formula and its roots | C09 | knowledge/primers.md: complex square roots and the quadratic formula (ch06 P159) — reminder |
| discrete root r, wiggles for R_cell > 2, δ = O(Δx/R_cell) (10.92) | C09 | C09 (D16) |
| nonuniform (stretched) grids | C09 | C09 (N47) |
| first-order upwind for the steady problem (10.93) | C09 | C09 (N48) |
| upwind modified equation, numerical diffusivity (10.94) | C09 | C09 (N49, D17) |
| higher-order upwind, SUPG, stabilised FE (named) | C09 | C09 (N50) |
| tridiagonal solve (Thomas algorithm) | C09 | knowledge/primers.md: Crank–Nicolson with `scipy.linalg.solve_banded` (ch08 P193) — reminder |
| incompressibility as a constraint, pressure as multiplier, infinite sound speed | C10 | C10 (N51) |
| Lagrange multiplier | C10 | primer (in C10) |
| artificial compressibility (10.95), pseudo-time | C10 | C10 (N52) |
| compressible Navier–Stokes in 2-D with μ_v = 0 (10.96)–(10.98) | C10 | R09 |
| Stokes' hypothesis (bulk viscosity zero) | C10 | Ch. 4 §4.5 (recalled in R09) |
| isothermal equation of state p = c²ρ (10.99), low-Mach error O(Ma²) | C10 | C10 (N53) |
| Mach number and the incompressibility criterion Ma < 0.3 | C10 | Ch. 4 §4.2 (recalled in C10, N53) |
| sound speed | C10 | Ch. 1 §1.8 (recalled in C10) |
| conservation (flux) form U_t + E_x + F_y = 0 | C10 | primer (in C10) |
| predictor–corrector (Heun's second-order idea) | C10 | primer (in C10) |
| RK4 by hand | C10 | knowledge/primers.md: RK4 by hand (ch03 P95) — reminder |
| MacCormack predictor–corrector (10.100)–(10.102) | C10 | C10 (D18) |
| Lax–Wendroff scheme | C10 | C10 (D18) |
| Schwarz's theorem (swap ∂t and ∂x) | C10 | knowledge/primers.md: Schwarz's theorem (ch04 P121) — reminder |
| dispersion (phase error) vs numerical diffusion | C10 | C10 (figure notes, D18 step 12) |
| NS predictor, corrector and coefficients c₁–c₅ (10.103)–(10.109) | C10 | C10 (N54, N55, N56) |
| cross-derivative stencil | C10 | C10 (N56) |
| FF/BB … arrangements and cycling | C10 | C10 (N57) |
| MacCormack time-step limit (10.110), mesh Reynolds number Re_Δ | C10 | C10 (N58) |
| density boundary conditions as the key issue | C10 | C10 (N59) |
| conservative convective term (10.82) | C11 | R08 |
| product rule for a divergence | C11 | knowledge/primers.md: product rule for a divergence (ch04 P113) — reminder |
| operator splitting (10.111)–(10.112) | C11 | C11 (N60) |
| operator splitting and the commutator [A₁, A₂] | C11 | primer (in C11) |
| scipy.linalg.expm | C11 | knowledge/primers.md: scipy.linalg.expm (ch02 P79) — reminder |
| matrix multiplication | C11 | knowledge/primers.md: matrix multiplication, transpose and identity (ch02 P63) — reminder |
| Marchuk–Yanenko fractional steps (10.113)–(10.114) | C11 | C11 (N61) |
| MAC split A₁, A₂ (10.115) | C11 | C11 (N62) |
| explicit convection–diffusion predictor (10.116) | C11 | C11 (N63) |
| projection step (10.117)–(10.118) | C11 | C11 (D19) |
| Helmholtz–Hodge decomposition | C11 | primer (in C11) |
| Poisson equation | C11 | knowledge/primers.md: Poisson equation and Green's function (ch05 P139) — reminder |
| divergence of a gradient = Laplacian, curl of a gradient = 0 | C11 | Ch. 2 §2.9 (recalled in C11, D19) |
| Gauss' divergence theorem | C11 | Ch. 2 §2.12 (recalled in C11, D19 and C12, N66) |
| projection method (Chorin, Temam), irrotational correction | C11 | C11 (N72) |
| initial and boundary conditions for NS (10.83), pressure up to a constant | C11 | C11 (N41) |
| Glowinski Θ-scheme (10.129)–(10.133), Θ = 1 − 1/√2 | C11 | C11 (N74) |
| dynamics/physics splitting in climate models | C11 | C11 (N60) |
| staggered grid (Fig. 10.4), Arakawa C-grid | C12 | C12 |
| half-index notation and staggered array shapes | C12 | primer (in C12) |
| np.meshgrid and the [j, i] layout | C12 | knowledge/primers.md: np.meshgrid and the project grid layout (ch02 P76) — reminder |
| staggered predictor with face averages (10.119)–(10.120) | C12 | C12 (N64) |
| velocity correction (10.121)–(10.122) | C12 | C12 (N65) |
| discrete continuity per cell (10.123) | C12 | C12 (N66) |
| singular linear systems and the compatibility condition | C12 | primer (in C12) |
| null space and rank | C12 | knowledge/primers.md: null space and rank–nullity theorem (ch01 P58) — reminder |
| scipy.sparse.diags, kron and a cached splu factorisation | C12 | primer (in C12) |
| discrete pressure Poisson equation (10.124) | C12 | C12 (N67, D20) |
| no pressure boundary condition on the staggered grid | C12 | C12 (N70, D20) |
| 5-point Laplacian and relaxation solvers (Jacobi, Gauss–Seidel, SOR) | C12 | knowledge/primers.md: iterative solvers: Jacobi, Gauss–Seidel, SOR (ch06 P160) — reminder |
| checkerboard pressure, collocated gradient (10.125) | C12 | C12 (N68, D21) |
| staggered face gradient (10.126) | C12 | C12 (N69, D21) |
| null space by SVD (count the tiny singular values) | C12 | primer (in C12) |
| MAC algorithm: predictor → Poisson → correction | C12 | C12 (N71) |
| timing code with time.perf_counter | C12 | C12 (gloss: a wall-clock timer in seconds) |
| MAC stability limits (10.127)–(10.128) | C12 | C12 (N73) |
| 2-D von Neumann scan | C12 | C12 (N73; the 1-D idea is C04) |
| Taylor–Green vortex and plane Poiseuille as exact tests | C12 | Ch. 4 §4.6 and Ch. 8 §8.2 (recalled in C12, verification cell) |
| rate-of-strain tensor D[u] (10.136) | C13 | R10 |
| weak Navier–Stokes (10.134)–(10.135) | C13 | C13 (N75, D22) |
| double dot product A:B | C13 | Ch. 2 §2.5 (recalled in C13, D22) |
| symmetric and antisymmetric parts of a tensor | C13 | Ch. 2 §2.10 (recalled in C13, D22 step 8) |
| saddle-point (KKT) matrix | C13 | primer (in C13) |
| semi-discrete saddle-point system (10.137), B and Bᵀ | C13 | C13 (N76) |
| GMRES in one line | C13 | primer (in C13) |
| sparsity pattern (plt.spy) | C13 | C13 (gloss: a dot for every nonzero entry) |
| generalised symmetric eigenproblem scipy.linalg.eigh(A, B) | C13 | primer (in C13) |
| eigenvalues and eigenvectors | C13 | knowledge/primers.md: eigenvalues and eigenvectors (ch02 P80) — reminder |
| LBB (inf–sup) condition, spurious pressure modes | C13 | C13 |
| Taylor–Hood P2–P1, iso-P2/P1, GLS stabilisation | C13 | C13 |
| counting velocity and pressure unknowns | C13 | C13 (tiny example) |
| barycentric (area) coordinates on a triangle | C13 | primer (in C13) |
| Jacobian determinant of a 2-D map | C13 | primer (in C13) |
| quadratic triangle shape functions (10.185), linear pressure shapes (10.187) | C13 | C13 (N105, N106) |
| isoparametric map (10.184), curved elements | C13 | C13 (N104) |
| seven-point triangle quadrature (10.198), exactness degree | C13 | C13 (N108, from-scratch cell) |
| cylinder-in-channel domain and mesh (Figs. 10.15–10.16 analogues) | C13 | C13 (N93, N94) |
| Euler's formula for meshes V − E + T | C13 | C13 (N94 gloss: counts of a triangulation with one hole) |
| plt.triplot and triangulations | C13 | C13 (gloss: draws the triangles of a mesh) |
| Cartesian weak equations (10.156)–(10.159), the sign of (10.159) | C13 | C13 (N95) |
| Galerkin statements on Ωʰ (10.160)–(10.162) | C13 | C13 (N96) |
| time derivative at t_{n+1} (10.163), backward Euler vs trapezoidal | C13 | C13 (N97) |
| Newton's method | C13 | knowledge/primers.md: Newton's method for complex zeros (ch06 P152) — reminder |
| Newton linearisation (10.164)–(10.167), residual as right side | C13 | C13 (N98, N99) |
| expansions and algebraic equations per node (10.168)–(10.172) | C13 | C13 (N100, N101) |
| block matrix form (10.173)–(10.183) | C13 | C13 (N102, N103) |
| element matrices (10.188)–(10.197) | C13 | C13 (N107) |
| assembly by node map, row replacement and penalty for essential BCs | C13 | C13 (N109) |
| cylinder wake regimes, steady eddies | C13 | R12 |
| creeping (Stokes) flow | C13 | Ch. 8 §8.6 (recalled in R12) |
| supercritical Hopf bifurcation (named) | C13 | C13 (N110) |
| Kármán vortex street | C13 | Ch. 9 §9.8 (recalled in C13, N110) |
| drag, lift, torque histories; Strouhal with the cyclic frequency | C13 | R13 |
| confined vs unbounded cylinder comparison | C13 | R13 |
| np.fft.rfft for a dominant frequency | C13 | knowledge/primers.md: `np.fft.rfft` and `np.fft.rfftfreq` (ch07 P181) — reminder |
| zero crossings with np.sign | C13 | knowledge/primers.md: `np.sign` and `np.nonzero`: where a curve crosses zero (ch07 P180) — reminder |
| lid-driven cavity set-up, scales, p = ρ/Ma² (Fig. 10.6 analogue) | C14 | C14 (N78) |
| corner singularities of the cavity | C14 | C14 (N78) |
| density on the walls from continuity (10.138)–(10.146) | C14 | C14 (N79, N81, N82, N83, N84) |
| one-sided second-order first derivative | C14 | C14 (N80; weights by the C01 method) |
| six-substep MacCormack cavity algorithm, coefficients a₁–a₁₁ | C14 | C14 (N85) |
| verification versus validation | C14 | primer (in C14) |
| benchmark (Ghia 1982, Hou 1995) | C14 | C14 |
| reading reference data from a file and interpolating (np.interp) | C14 | primer (in C14) |
| np.interp | C14 | knowledge/primers.md: `np.interp`: reading a curve between samples (ch07 P182) — reminder |
| caching expensive runs (np.savez and a parameter key) | C14 | primer (in C14) |
| streamfunction of a computed field (ψ contours) | C14 | Ch. 4 §4.3 and Ch. 6 §6.7 (recalled in C14: ψ from u, `solve_poisson`) |
| vorticity of a computed field | C14 | Ch. 5 §5.1 (recalled in C14) |
| plt.contour, plt.quiver, plt.streamplot | C14 | knowledge/primers.md: plt.contour, plt.quiver and plt.streamplot (ch02 P78) — reminder |
| steady state test max∣Δu∣/Δt | C14 | C14 (gloss in the code cell) |
| overview of the three examples | C14 | C14 (N77) |
| drag and lift coefficients (4.107)–(4.108) | C15 | R11 |
| block in a channel, sliding walls, Galilean frame change | C15 | C15 (N86) |
| frames of reference and relative velocity | C15 | knowledge/primers.md: frames of reference and relative velocity (ch03 P96) — reminder |
| front-face density from viscous terms (10.147) | C15 | C15 (N87) |
| one-sided stencils (10.148)–(10.150) | C15 | C15 (N88) |
| wall densities on the block (10.151)–(10.154), corners averaged | C15 | C15 (N89) |
| double precision, ρ′ = ρ − 1, cycling | C15 | C15 (N90) |
| asymptotic MacCormack step (10.155) and its hidden Ma ≪ 1 | C15 | C15 (N91) |
| grid convergence, observed order from three grids | C15 | C15 (D23) |
| Richardson extrapolation | C15 | C15 (D23) |
| grid convergence index (GCI) | C15 | C15 (D23 step 8) |
| tools.convergence.richardson, grid_convergence_index | C15 | C15 (code cell gloss) |
| cylinder at Re = 1000 with a smoke line (named, not reproduced) | C15 | C15 (N92) |
| verification checklist (benchmark, mesh, time step, boundaries, parameters) | C15 | C15 (N111) |
| other methods: spectral, spectral element, lattice Boltzmann, DPD | C15 | C15 (N112) |
| exercises of the book (not reproduced) | C15 | C15 (pointer, S01) |
| literature and supplemental reading | C15 | C15 (pointer, S02) |

## Part F — derivation storyboards

Builders copy these word for word into `nb.derivation(key, title, goal=…, start=(tex, plain), plan=[…], uses=[…],
steps=[dict(did, tex, why, plain)], result=(tex, plain), interpret=…, check=…, check_src=…)` and into the explainer's
`derivations: [...]` (phones may shorten *why* to its first sentence; `live`, `set` and `watch` are the explainer's and are listed
in Part B). Every step is one move; *why* names the rule and says why we make it (≤ 35 words); *did* ≤ 8 words. The book's own moves
were read on the rendered pages (p451 for D01; p452 for D02; p453 for D03; p454–p455 for D04–D06; D07 and D08 are ours (the book
asserts BTCS stability and states CFL); p456 for D09; p457 for D10; p457–p459 for D11; p460 for D12; p461–p462 for D13–D14; p464–p465
for D15–D16; p466 for D17; p467 for D18 (the book only says 'second order'); D19 is ours (the book writes only the discrete (10.124));
p471–p472 for D20–D21; p474 for D22; D23 is ours from (10.15)); the gaps listed in `analysis/ch10.md` §2b are filled and the notebook
says so ("the book skips this move; we add it"). Every equation named by number is written out. Colours: kept derivative blue,
cancelled terms muted, leading error rose; convection orange, diffusion teal, time amber; pressure orange, velocity blue. **Index
convention:** in D04–D08 our own lines use the grid index j and the imaginary unit $\mathrm i$ (upright), so i means only
$\sqrt{-1}$; the book's quoted equations keep their i. No line of this part starts with a table bar; absolute values are written
with \lvert \rvert or in words. The ★★★ derivation (D18) carries a `check_src` cell, every line commented; D03, D06 and D16 carry an
optional one.

### D01 · Stencils and their errors: forward, backward and centred $\big[\frac{\partial T}{\partial x}\big]_i$ (10.6) and the centred $\big[\frac{\partial^2T}{\partial x^2}\big]_i$ (10.7) from the Taylor series (10.4)–(10.5) — ★, 9 steps, in C01 (notebook · `fd_stencil_order`)
- **Goal.** Find how to compute a slope and a curvature from grid values, and exactly how wrong each formula is.
- **Start.** The Taylor series at the two neighbours, $T^n_{i+1}=T^n_i+\Delta x\,T_x+\frac{\Delta x^2}2T_{xx}+\frac{\Delta x^3}6T_{xxx}+\frac{\Delta
  x^4}{24}T_{xxxx}+O(\Delta x^5)$ (10.4) and $T^n_{i-1}=T^n_i-\Delta x\,T_x+\frac{\Delta x^2}2T_{xx}-\frac{\Delta x^3}6T_{xxx}+\frac{\Delta x^4}{24}
  T_{xxxx}+O(\Delta x^5)$ (10.5), all derivatives at $(x_i,t_n)$ — *in words:* what the neighbours know about the point.
- **Plan.** (1) Solve one series for $T_x$ (forward and backward). (2) Subtract the two series: the even terms cancel (centred).
  (3) Add them: the odd terms cancel (second derivative). (4) Read the power of Δx of the first term left over.
- **Tools.** Taylor series (R02; Ch. 1 P26) · big-O notation and the order of accuracy (primer P221).
- **Assumptions.** T has four bounded derivatives near $x_i$ (used in step 1); uniform spacing Δx (steps 5, 7).
- **Steps.**
  1. *did:* Start from the right neighbour · *tex:* $T_{i+1}=T_i+\Delta x\,T_x+\frac{\Delta x^2}2T_{xx}+\frac{\Delta x^3}6T_{xxx}+O(\Delta x^4)$ ·
     *why:* Taylor's theorem about $x_i$ (R02), valid because T is smooth; we keep enough terms to see the first error of every stencil. ·
     *plain:* the value next door is the value here plus slope, curvature and more, each times a power of Δx.
  2. *did:* Solve for the slope · *tex:* $T_x=\frac{T_{i+1}-T_i}{\Delta x}-\frac{\Delta x}2T_{xx}-\frac{\Delta x^2}6T_{xxx}+O(\Delta x^3)$ · *why:*
     Subtract $T_i$ from both sides and divide by Δx ≠ 0 — plain algebra; we want $T_x$ alone. · *plain:* the slope is 'right minus
     here over Δx' minus a correction.
  3. *did:* Name the forward difference and its error · *tex:* $\frac{T_{i+1}-T_i}{\Delta x}-T_x=+\frac{\Delta x}2T_{xx}+O(\Delta x^2)$ · *why:*
     Rearrange step 2 as stencil minus exact; the first surviving term is proportional to Δx¹, so the forward difference is first order
     (P221) — the O(Δx) of (10.6). · *plain:* halve Δx and the forward error halves.
  4. *did:* Repeat with the left neighbour (10.5) · *tex:* $\frac{T_i-T_{i-1}}{\Delta x}-T_x=-\frac{\Delta x}2T_{xx}+O(\Delta x^2)$ · *why:* The
     same two moves on (10.5), where the odd powers carry minus signs; the error has the opposite sign but the same size. · *plain:* the
     backward difference is also first order.
  5. *did:* Subtract (10.5) from (10.4) · *tex:* $T_{i+1}-T_{i-1}=2\Delta x\,T_x+\frac{\Delta x^3}3T_{xxx}+O(\Delta x^5)$ · *why:* Terms with even
     powers (T, Δx²T_xx, Δx⁴T_xxxx) have the same sign in both series and cancel; odd ones double. The curvature term is gone. · *plain:*
     'right minus left' contains twice the slope and nothing of the curvature.
  6. *did:* Divide by 2Δx · *tex:* $\frac{T_{i+1}-T_{i-1}}{2\Delta x}-T_x=\frac{\Delta x^2}6T_{xxx}+O(\Delta x^4)$ · *why:* Isolates the stencil
     minus the slope; the first survivor is proportional to Δx², so the centred difference is second order — the third line of (10.6). ·
     *plain:* halve Δx and the centred error drops by four.
  7. *did:* Add (10.4) and (10.5) · *tex:* $T_{i+1}+T_{i-1}=2T_i+\Delta x^2T_{xx}+\frac{\Delta x^4}{12}T_{xxxx}+O(\Delta x^6)$ · *why:* Now the
     odd powers cancel (the Δx⁵ terms too, being odd), and even ones double: 2 × Δx⁴/24 = Δx⁴/12. · *plain:* the sum of the neighbours
     sees the curvature, not the slope.
  8. *did:* Solve for the curvature · *tex:* $\frac{T_{i+1}-2T_i+T_{i-1}}{\Delta x^2}-T_{xx}=\frac{\Delta x^2}{12}T_{xxxx}+O(\Delta x^4)$ · *why:*
     Subtract 2T_i and divide by Δx²; the leftover starts at Δx² (dividing O(Δx⁶) by Δx² leaves O(Δx⁴)) — second order, (10.7). · *plain:*
     the three-point curvature formula is second order.
  9. *did:* Collect the orders · *tex:* $\big[T_x\big]_i=\frac{T_{i+1}-T_i}{\Delta x}+O(\Delta x)=\frac{T_i-T_{i-1}}{\Delta x}+O(\Delta x)=
     \frac{T_{i+1}-T_{i-1}}{2\Delta x}+O(\Delta x^2)$ · *why:* Steps 3, 4 and 6 written as the book writes (10.6), with the error terms folded
     into O(·); (10.7) is step 8. · *plain:* one-sided slopes are first order, centred ones second.
- **Result.** (10.6) and (10.7) with their leading errors: forward $+\frac{\Delta x}2T_{xx}$, backward $-\frac{\Delta x}2T_{xx}$, centred
  $\frac{\Delta x^2}6T_{xxx}$, second derivative $\frac{\Delta x^2}{12}T_{xxxx}$ — *in words:* the order is the power of Δx of the first Taylor
  term a stencil fails to cancel.
- **Check.** Units: each stencil has the units of the derivative it approximates (T/m or T/m²) ✓. Numbers (sin at 1, Δx = 0.1): forward
  error −0.0429 vs predicted (0.05)(−0.841) = −0.0421; centred −9.001e-4 vs (0.01/6)(−0.5403) = −9.005e-4 ✓. `FD.stencil_taylor_coefficients`
  lists the same coefficients (1/2, 1/6, 1/12).
- **What it means.** Every scheme of the chapter is assembled from these stencils; their orders set the scheme's order (C03). Fails
  where T is not smooth (a shock, a kink) — there no stencil beats first order.
- **Traps.** Getting the sign of the forward error backwards (it depends on whether you write stencil − exact or exact − stencil);
  claiming O(Δx³) for (10.7) because the series was cut at O(Δx⁵) — the Δx⁴ term is kept, so the error is O(Δx²); forgetting the 2 in
  2Δx.

### D02 · The explicit FTCS update $T^{n+1}_i=T^n_i-\alpha(T^n_{i+1}-T^n_{i-1})+\beta(T^n_{i+1}-2T^n_i+T^n_{i-1})$ (10.10) with $\alpha=u\frac{\Delta t}{2\Delta x}$, $\beta=D\frac{\Delta t}{\Delta x^2}$ (10.11) — ★, 6 steps, in C02 (notebook)
- **Goal.** Turn the PDE (10.1) into a rule that computes the new temperature at each point from the old ones.
- **Start.** $\frac{\partial T}{\partial t}+u\frac{\partial T}{\partial x}=D\frac{\partial^2T}{\partial x^2}$ (10.1) — *in words:* change in time =
  carried in by the current + spread by diffusion.
- **Plan.** (1) Ask the PDE at one grid point. (2) Replace each derivative by a stencil. (3) Solve for the unknown new value.
- **Tools.** The time difference (10.8) (N07) · the stencils (10.6)–(10.7) (D01).
- **Assumptions.** u and D constant (step 5 names constant coefficients); a uniform grid.
- **Steps.**
  1. *did:* Evaluate (10.1) at a grid point · *tex:* $\big[T_t\big]^n_i+u\big[T_x\big]^n_i=D\big[T_{xx}\big]^n_i$ · *why:* A discretisation
     demands the equation only at $(x_i,t_n)$; we choose the old time level so every space derivative uses known values. · *plain:* the
     PDE, asked at one point of the grid.
  2. *did:* Forward difference in time · *tex:* $\big[T_t\big]^n_i=\frac{T^{n+1}_i-T^n_i}{\Delta t}+O(\Delta t)$ · *why:* The first line of (10.8),
     derived like D01 step 3 but in t; it is the only choice that puts a single unknown, $T^{n+1}_i$, in each equation. · *plain:* the rate of
     change is 'next minus now over Δt'.
  3. *did:* Centred differences in space · *tex:* $\frac{T^{n+1}_i-T^n_i}{\Delta t}+u\frac{T^n_{i+1}-T^n_{i-1}}{2\Delta x}=D\frac{T^n_{i+1}-2T^n_i+
     T^n_{i-1}}{\Delta x^2}+O(\Delta t,\Delta x^2)$ (10.9) · *why:* Substitute (10.6) (centred) and (10.7); their errors, O(Δx²), join the
     O(Δt) of step 2. · *plain:* forward in time, centred in space — FTCS.
  4. *did:* Multiply through by Δt · *tex:* $T^{n+1}_i-T^n_i=-\frac{u\Delta t}{2\Delta x}(T^n_{i+1}-T^n_{i-1})+\frac{D\Delta t}{\Delta x^2}(T^n_{i+1}-
     2T^n_i+T^n_{i-1})+O(\Delta t^2,\Delta t\Delta x^2)$ · *why:* Clears the fraction on the left and moves the convective term right; the
     remainder is multiplied by Δt too. · *plain:* the change over one step, written out.
  5. *did:* Name the two numbers α, β · *tex:* $\alpha=u\frac{\Delta t}{2\Delta x},\quad\beta=D\frac{\Delta t}{\Delta x^2}$ (10.11) · *why:* They
     are dimensionless and constant (u, D, Δx, Δt fixed); the 2 in α comes from the centred difference's 2Δx. · *plain:* α measures how far
     the current moves in a step, β how far diffusion reaches.
  6. *did:* Drop the remainder, solve for the new value · *tex:* $T^{n+1}_i\approx T^n_i-\alpha(T^n_{i+1}-T^n_{i-1})+\beta(T^n_{i+1}-2T^n_i+T^n_{i-1})$
     (10.10) · *why:* Discarding O(Δt², ΔtΔx²) turns '=' into '≈'; the new value appears only on the left, so each node is updated by a
     formula — an explicit scheme. · *plain:* new value = old value, pushed by the current, smoothed by diffusion.
- **Result.** (10.10) with (10.11) — *in words:* each new value is a fixed weighted combination of three old ones.
- **Check.** α, β dimensionless: (m/s)(s)/m and (m²/s)(s)/m² ✓. The five-point example: (0.1, 0.2) gives [0, 0.1, 0.6, 0.3, 0], the sum
  stays 1 (the weights (α + β) + (1 − 2β) + (β − α) = 1). u = 0 recovers Ch. 1's FTCS with r = β.
- **What it means.** The simplest scheme there is; C03 asks how close it is to (10.1) and C04 when it blows up.
- **Traps.** Writing α = uΔt/Δx (forgetting the 2 of the centred difference); keeping '=' after dropping the error; mixing time levels
  (a $T^{n+1}$ on the right would make the scheme implicit — that is BTCS, N08).

### D03 · The truncation error of FTCS $E^n_i=\frac{\Delta t}2\big[\frac{\partial^2T}{\partial t^2}\big]^n_i+u\frac{\Delta x^2}6\big[\frac{\partial^3T}{\partial x^3}\big]^n_i-D\frac{\Delta x^2}{12}\big[\frac{\partial^4T}{\partial x^4}\big]^n_i+O(\Delta t^2,\Delta x^4)$ (10.16)–(10.17) — ★★, 11 steps, in C03 (notebook · `fd_stencil_order`)
- **Goal.** Find what is left over when the exact solution is put into the FTCS rule — the truncation error — and show that it vanishes
  as the grid shrinks (consistency).
- **Start.** The FTCS scheme (10.9) without its error term, applied to a smooth function T: $\mathcal R=\frac{T^{n+1}_i-T^n_i}{\Delta t}+u
  \frac{T^n_{i+1}-T^n_{i-1}}{2\Delta x}-D\frac{T^n_{i+1}-2T^n_i+T^n_{i-1}}{\Delta x^2}$ — *in words:* the scheme's residual.
- **Plan.** (1) Expand the time difference in Taylor series (the book skips this). (2) Reuse D01 for the space differences. (3) Collect:
  the PDE plus a leftover. (4) Read the orders, then look at a special case.
- **Tools.** Taylor series in t and x (Ch. 3 P98) · D01's results · sympy `series` (Ch. 4 P117) for the optional check.
- **Assumptions.** T smooth in x and t (steps 2–5); constant u, D.
- **Steps.**
  1. *did:* Put a smooth T into the scheme · *tex:* $\mathcal R=\frac{T^{n+1}_i-T^n_i}{\Delta t}+u\frac{T^n_{i+1}-T^n_{i-1}}{2\Delta x}-D\frac{T^n_{i+1}
     -2T^n_i+T^n_{i-1}}{\Delta x^2}$ · *why:* Consistency asks what the scheme's left side becomes for a smooth function; we will expand each
     difference about $(x_i,t_n)$. · *plain:* measure how badly a smooth function fails the rule.
  2. *did:* Taylor-expand the new time level · *tex:* $T^{n+1}_i=T+\Delta t\,T_t+\frac{\Delta t^2}2T_{tt}+\frac{\Delta t^3}6T_{ttt}+O(\Delta t^4)$ ·
     *why:* Taylor in t at fixed $x_i$ (P98). The book never writes this line — we add it; it is the source of the Δt/2 term. · *plain:*
     the value one step later, as a series in Δt.
  3. *did:* Form the time difference · *tex:* $\frac{T^{n+1}_i-T^n_i}{\Delta t}=T_t+\frac{\Delta t}2T_{tt}+O(\Delta t^2)$ · *why:* Subtract T and
     divide by Δt — the same moves as D01 step 2 in time. · *plain:* the forward difference is the rate plus a Δt/2 correction.
  4. *did:* Reuse the centred first difference · *tex:* $\frac{T^n_{i+1}-T^n_{i-1}}{2\Delta x}=T_x+\frac{\Delta x^2}6T_{xxx}+O(\Delta x^4)$ · *why:*
     D01 step 6, read the other way round. · *plain:* the convective stencil is the slope plus a Δx² correction.
  5. *did:* Reuse the centred second difference · *tex:* $\frac{T^n_{i+1}-2T^n_i+T^n_{i-1}}{\Delta x^2}=T_{xx}+\frac{\Delta x^2}{12}T_{xxxx}+O(\Delta x^4)$ ·
     *why:* D01 step 8, read the other way round. · *plain:* the diffusive stencil is the curvature plus a Δx² correction.
  6. *did:* Substitute steps 3–5 into $\mathcal R$ · *tex:* $\mathcal R=T_t+uT_x-DT_{xx}+\frac{\Delta t}2T_{tt}+u\frac{\Delta x^2}6T_{xxx}-D\frac{\Delta
     x^2}{12}T_{xxxx}+O(\Delta t^2,\Delta x^4)$ · *why:* Each difference replaced by its series; the diffusive correction enters with −D because
     of the minus sign in $\mathcal R$. · *plain:* the residual is the PDE plus three small terms.
  7. *did:* Name the leftover $E^n_i$ · *tex:* $E^n_i=\frac{\Delta t}2T_{tt}+u\frac{\Delta x^2}6T_{xxx}-D\frac{\Delta x^2}{12}T_{xxxx}+O(\Delta t^2,
     \Delta x^4)$ (10.17) · *why:* Collect everything that is not the PDE; this is the truncation error. · *plain:* E is what the scheme adds to
     the equation.
  8. *did:* Read the scheme as a modified equation · *tex:* $\big[T_t\big]^n_i+u\big[T_x\big]^n_i-D\big[T_{xx}\big]^n_i+E^n_i=0$ (10.16) · *why:* The
     numerical solution satisfies the scheme exactly, $\mathcal R=0$; if it is smooth, step 6 says it obeys the PDE plus E — the book's
     sign convention (E on the left). · *plain:* FTCS solves (10.1) with E added.
  9. *did:* Let Δx, Δt shrink · *tex:* $E^n_i\to0\quad(\Delta t,\Delta x\to0)$ · *why:* Every term of E carries a positive power of Δt or Δx and
     bounded derivatives; hence (10.16) tends to (10.1): the scheme is consistent. · *plain:* on a fine enough grid the scheme is the PDE.
  10. *did:* Read the orders · *tex:* $E=O(\Delta t)+O(\Delta x^2)$ · *why:* The lowest powers in (10.17) are Δt¹ and Δx²: first order in time,
      second in space, as the book states after (10.17). · *plain:* halving Δt halves the time error; halving Δx quarters the space error.
  11. *did:* Special case: pure diffusion · *tex:* $E=D\,\Delta x^2\Big(\frac\beta2-\frac1{12}\Big)T_{xxxx}$ · *why:* With u = 0 the PDE gives
      $T_{tt}=D\,\partial_{xx}T_t=D^2T_{xxxx}$ and $\Delta t=\beta\Delta x^2/D$; substituting into (10.17) combines the two terms. Ours, not the
      book's. · *plain:* at β = 1/6 the two leading errors cancel exactly.
- **Result.** (10.16)–(10.17) — *in words:* FTCS is consistent, first order in time and second order in space; its leftover contains a
  time term, a dispersive term (odd derivative) and a diffusive term (even derivative).
- **Check.** Units: every term of E has units of T/s (e.g. $\frac{\Delta t}2T_{tt}$: s · T/s²) ✓. Tiny example (u = 0, D = 1, Δx = 0.1,
  β = 0.25, T = e^{−π²t} sin πx at x = ½, t = 0): E = π⁴ × 0.01 × (0.125 − 0.0833) = 0.0406 ✓; `FD.truncation_error_sympy("ftcs")` returns
  (10.17); the measured one-step residual on the Gaussian agrees with `FD.truncation_terms` within 5 %.
- **Optional sympy check (`check_src`, every line commented):** build T(x, t) as a symbolic function, form $\mathcal R$ with `subs`, call
  `series` in Δt and Δx to second order, subtract $T_t+uT_x-DT_{xx}$ and print the remainder — the three terms of (10.17).
- **What it means.** Consistency is the first of the three properties; the modified-equation reading (step 8) returns in D17 (upwind's
  numerical diffusion) and D18 (MacCormack's dispersion). Fails if the solution is not smooth.
- **Traps.** Leaving out the time Taylor series (then the Δt/2 T_tt term is missing and FTCS looks second order in time); the sign
  convention (E on the left of (10.16)); thinking consistency guarantees a good answer (C04).

### D04 · The error equation $\xi^{n+1}_i=(\alpha+\beta)\xi^n_{i-1}+(1-2\beta)\xi^n_i+(\beta-\alpha)\xi^n_{i+1}$ (10.19) and the amplification factor $\frac{g^{n+1}}{g^n}=(\alpha+\beta)e^{-\mathrm i\pi k\Delta x}+(1-2\beta)+(\beta-\alpha)e^{\mathrm i\pi k\Delta x}$ (10.24) — ★★, 11 steps, in C04 (notebook · `von_neumann_amplification`)
- **Goal.** Show that a round-off error obeys the scheme itself, and that each of its Fourier waves is multiplied by one complex number
  per step.
- **Start.** The FTCS rule (10.10), satisfied by the exact solution T of the discrete system and — up to the round-off made once — by
  the computed $\overline T$; the disturbance $\xi^n_i=T^n_i-\overline T^n_i$ (10.18) — *in words:* two runs of the same rule that differ by a
  tiny error.
- **Plan.** (1) Subtract the two runs (the book: 'it can be shown'). (2) Write the error as a sum of waves and follow one. (3) Divide
  out the wave shape to get G.
- **Tools.** Linearity (superposition) · Fourier modes (R04) · Euler's formula and complex exponentials (Ch. 1 P45, Ch. 6 P153).
- **Assumptions.** The scheme is linear with constant α, β (step 3); a periodic or infinite grid, so every wave fits (step 5); uniform
  spacing (step 8). **Index renamed:** grid index j from step 5 on, so that i is only $\sqrt{-1}$.
- **Steps.**
  1. *did:* Write FTCS for the exact discrete solution · *tex:* $T^{n+1}_i=T^n_i-\alpha(T^n_{i+1}-T^n_{i-1})+\beta(T^n_{i+1}-2T^n_i+T^n_{i-1})$ · *why:*
     (10.10) with '=' because T is by definition the exact solution of the discrete system. · *plain:* the rule, obeyed perfectly.
  2. *did:* Write it for the computed values · *tex:* $\overline T^{\,n+1}_i=\overline T^n_i-\alpha(\overline T^n_{i+1}-\overline T^n_{i-1})+\beta(\overline
     T^n_{i+1}-2\overline T^n_i+\overline T^n_{i-1})$ · *why:* The computer applies the same rule to its own (slightly wrong) numbers; we follow
     how an error present at step n travels on. · *plain:* the rule, applied to the numbers we actually hold.
  3. *did:* Subtract step 2 from step 1 · *tex:* $\xi^{n+1}_i=\xi^n_i-\alpha(\xi^n_{i+1}-\xi^n_{i-1})+\beta(\xi^n_{i+1}-2\xi^n_i+\xi^n_{i-1})$ · *why:*
     The rule is linear, so the difference of two solutions is again a solution, with ξ = T − T̄ (10.18). This is the move the book
     skips. · *plain:* the error obeys the same rule as the temperature.
  4. *did:* Collect the coefficients · *tex:* $\xi^{n+1}_i=(\alpha+\beta)\xi^n_{i-1}+(1-2\beta)\xi^n_i+(\beta-\alpha)\xi^n_{i+1}$ (10.19) · *why:* Gather
     the terms of $\xi_{i-1}$ (+α + β), $\xi_i$ (1 − 2β) and $\xi_{i+1}$ (−α + β). · *plain:* each new error is a weighted mix of three old ones.
  5. *did:* Follow one Fourier wave · *tex:* $\xi^n_j=\hat\xi^n\,e^{\mathrm i\pi kx_j}$ (10.21) · *why:* Any grid error is a sum of such waves
     (10.20) (R04); the rule is linear and a wave stays a wave, so each can be followed alone. We rename the index j and the amplitude
     $\hat\xi^n$ (the book's $g^n(k)$). · *plain:* one wave, with an amplitude that may change from step to step.
  6. *did:* The same wave one step later · *tex:* $\xi^{n+1}_j=\hat\xi^{n+1}e^{\mathrm i\pi kx_j}$ (10.22) · *why:* The shape in space is fixed by k;
     only the amplitude can change — that is what we want to find. · *plain:* same wave, new amplitude.
  7. *did:* Substitute into (10.19) · *tex:* $\hat\xi^{n+1}e^{\mathrm i\pi kx_j}=\hat\xi^n\big[(\alpha+\beta)e^{\mathrm i\pi kx_{j-1}}+(1-2\beta)e^{\mathrm i\pi kx_j}+
     (\beta-\alpha)e^{\mathrm i\pi kx_{j+1}}\big]$ (10.23) · *why:* Steps 5 and 6 inserted in step 4, point by point. · *plain:* the rule applied to
     the wave.
  8. *did:* Shift the neighbours · *tex:* $e^{\mathrm i\pi kx_{j\pm1}}=e^{\mathrm i\pi kx_j}\,e^{\pm\mathrm i\theta},\quad\theta=k\pi\Delta x$ · *why:* $x_{j\pm1}=x_j\pm
     \Delta x$ and $e^{a+b}=e^ae^b$; θ is the phase change per cell (the book's convention: k counts half-waves per unit length). · *plain:* a
     neighbour's wave value is ours, turned by θ.
  9. *did:* Divide by $\hat\xi^ne^{\mathrm i\pi kx_j}$ · *tex:* $G\equiv\frac{\hat\xi^{n+1}}{\hat\xi^n}=(\alpha+\beta)e^{-\mathrm i\theta}+(1-2\beta)+(\beta-\alpha)
     e^{\mathrm i\theta}$ (10.24) · *why:* $e^{\mathrm i\pi kx_j}$ is never zero and appears in every term, so it cancels; we assume $\hat\xi^n\ne0$ (a
     wave that is present). The result no longer depends on j. · *plain:* one complex number G per wave.
  10. *did:* Repeat n times · *tex:* $\hat\xi^n=G^n\hat\xi^0,\qquad\lvert\hat\xi^n\rvert=\lvert G\rvert^n\lvert\hat\xi^0\rvert$ · *why:* G is the same at every
      step (constant α, β), so the amplitude is multiplied by G again and again; the modulus of a product is the product of moduli. ·
      *plain:* after n steps the wave has grown or shrunk by $\lvert G\rvert^n$.
  11. *did:* Demand no growth for every wave · *tex:* $\Big\lvert\frac{g^{n+1}}{g^n}\Big\rvert=\lvert G(\theta)\rvert\le1\quad\text{for all }\theta\in[0,\pi]$
      (10.25) · *why:* One wave with $\lvert G\rvert>1$ grows without bound from round-off; $G(-\theta)$ is the complex conjugate of $G(\theta)$, so
      θ ∈ [0, π] covers every wave the grid can carry. · *plain:* stable means no wave grows.
- **Result.** (10.19), (10.24) and the condition (10.25) — *in words:* the error of a linear scheme obeys the scheme, and every wave of it
  is multiplied by G(θ) each step.
- **Check.** G is dimensionless ✓. θ = 0 (a constant error): G = (α + β) + (1 − 2β) + (β − α) = 1 — a constant offset neither grows nor
  decays ✓. θ = π/2, α = 0.1, β = 0.2: G = 0.6 − 0.2i ✓ (C04 tiny example).
- **What it means.** Stability is now a property of one function G(θ); D05–D08 read it for FTCS, BTCS and upwind. Fails for nonlinear
  schemes and variable coefficients (then it is a local, 'frozen-coefficient' guide only) and near boundaries.
- **Traps.** Confusing the index i with $\sqrt{-1}$ (we renamed it j); using θ = kΔx (the standard convention) with the book's k;
  dividing by zero when the wave is absent; forgetting that one bad θ is enough.

### D05 · The modulus $\lvert G\rvert^2=\big(1-4\beta\sin^2\frac\theta2\big)^2+(2\alpha\sin\theta)^2$ (10.26) — ★★, 8 steps, in C04 (notebook · `von_neumann_amplification`)
- **Goal.** Turn the complex G of (10.24) into a real formula whose size can be compared with 1.
- **Start.** $G=(\alpha+\beta)e^{-\mathrm i\theta}+(1-2\beta)+(\beta-\alpha)e^{\mathrm i\theta}$ (10.24) — *in words:* the factor one wave is multiplied by.
- **Plan.** (1) Split G into real and imaginary parts with Euler's formula. (2) Rewrite 1 − cos θ with a half angle. (3) Square and add.
- **Tools.** Euler's formula (Ch. 1 P45) · half-angle identities (primer P225) · the modulus of a complex number (Ch. 6 P153).
- **Assumptions.** α, β real (they are).
- **Steps.**
  1. *did:* Insert Euler's formula · *tex:* $G=(\alpha+\beta)(\cos\theta-\mathrm i\sin\theta)+(1-2\beta)+(\beta-\alpha)(\cos\theta+\mathrm i\sin\theta)$ ·
     *why:* $e^{\pm\mathrm i\theta}=\cos\theta\pm\mathrm i\sin\theta$ (P45); we want G as 'real + i·real'. · *plain:* the two turned neighbours written with
     cos and sin.
  2. *did:* Collect the real part · *tex:* $\mathrm{Re}\,G=(\alpha+\beta+\beta-\alpha)\cos\theta+1-2\beta=1-2\beta(1-\cos\theta)$ · *why:* The α terms
     cancel in the real part; 2β cos θ + 1 − 2β regroups as 1 − 2β(1 − cos θ). · *plain:* diffusion alone sets the real part.
  3. *did:* Collect the imaginary part · *tex:* $\mathrm{Im}\,G=\big[-(\alpha+\beta)+(\beta-\alpha)\big]\sin\theta=-2\alpha\sin\theta$ · *why:* Now the β
     terms cancel; only convection survives. · *plain:* convection alone sets the imaginary part.
  4. *did:* Use the half-angle identity · *tex:* $1-\cos\theta=2\sin^2\tfrac\theta2$ · *why:* P225; it turns a quantity running from 0 to 2 into
     2 × (a square between 0 and 1), which makes step 7 and D06 clean. · *plain:* a cosine rewritten as a square.
  5. *did:* Rewrite the real part · *tex:* $\mathrm{Re}\,G=1-4\beta\sin^2\tfrac\theta2$ · *why:* Substitute step 4 into step 2: 2β × 2 sin² = 4β sin². ·
     *plain:* the real part drops from 1 as the wave gets shorter.
  6. *did:* Assemble G · *tex:* $G=\big(1-4\beta\sin^2\tfrac\theta2\big)-2\mathrm i\alpha\sin\theta$ · *why:* Steps 5 and 3 together; a compact form for
     plotting (the orange curve of E2). · *plain:* G = (damping part) − i (turning part).
  7. *did:* Square and add · *tex:* $\lvert G\rvert^2=\big(1-4\beta\sin^2\tfrac\theta2\big)^2+(2\alpha\sin\theta)^2$ · *why:* $\lvert a+\mathrm ib\rvert^2=a^2+b^2$
     (P153); the sign of the imaginary part disappears in the square. · *plain:* the squared size of the factor.
  8. *did:* Apply the stability condition · *tex:* $\big(1-4\beta\sin^2\tfrac\theta2\big)^2+(2\alpha\sin\theta)^2\le1$ (10.26) · *why:* (10.25) squared (both
     sides are non-negative, so squaring keeps the inequality). · *plain:* the book's stability condition for FTCS, for every θ.
- **Result.** (10.26) — *in words:* the real part carries diffusion, the imaginary part convection, and the squared size must not exceed 1.
- **Check.** θ = π/2, α = 0.1, β = 0.2: (1 − 0.4)² + 0.2² = 0.40 ✓ (= ∣0.6 − 0.2i∣²). α = 0, θ = π: ∣G∣ = ∣1 − 4β∣, the zigzag factor ✓.
- **What it means.** A formula in two parameters and one angle — D06 turns "for every θ" into two inequalities.
- **Traps.** Writing the real part as 1 − 2β cos θ (losing the 1 − cos θ grouping); carrying the minus sign of Im G into the square.

### D06 · Noye's region $0\le4\alpha^2\le2\beta\le1$ (10.27), the diffusion limit $\Delta t\le\frac12\frac{\Delta x^2}D$ (10.28) and pure convection never stable — ★★, 13 steps, in C04 (notebook · `von_neumann_amplification`)
- **Goal.** Replace "∣G∣ ≤ 1 for every wave" by conditions on α and β alone — the book cites Noye (1983) without proof; we derive it.
- **Start.** $\big(1-4\beta\sin^2\frac\theta2\big)^2+(2\alpha\sin\theta)^2\le1$ for all θ ∈ [0, π] (10.26) — *in words:* no wave may grow.
- **Plan.** (1) Change variable to s = sin²(θ/2) ∈ [0, 1]. (2) Factor ∣G∣² − 1 = s × (a straight line in s). (3) A line is ≤ 0 on an interval
  iff it is ≤ 0 at both ends. (4) Read the special cases.
- **Tools.** Half-angle identities (P225) · sign of a linear function on an interval (primer P226).
- **Assumptions.** β ≥ 0 (D ≥ 0, a physical diffusivity); α any real number (its sign does not matter — only α² appears).
- **Steps.**
  1. *did:* Name $s=\sin^2(\theta/2)$ · *tex:* $s=\sin^2\tfrac\theta2\in[0,1]\quad\text{for }\theta\in[0,\pi]$ · *why:* As θ runs over [0, π], θ/2 runs over
     [0, π/2] and its sine squared over [0, 1]: one variable instead of two trigonometric functions. · *plain:* s = 0 is the longest
     wave, s = 1 the zigzag.
  2. *did:* Express sin²θ through s · *tex:* $\sin^2\theta=4\sin^2\tfrac\theta2\cos^2\tfrac\theta2=4s(1-s)$ · *why:* sin θ = 2 sin(θ/2)cos(θ/2) and
     cos² = 1 − sin² (P225). · *plain:* the convective term, in the same variable.
  3. *did:* Rewrite ∣G∣² in s · *tex:* $\lvert G\rvert^2=(1-4\beta s)^2+16\alpha^2s(1-s)$ · *why:* Substitute steps 1–2 into (10.26): $(2\alpha)^2\cdot4s(1-s)=
     16\alpha^2s(1-s)$. · *plain:* ∣G∣² is a quadratic in s.
  4. *did:* Expand the square · *tex:* $\lvert G\rvert^2=1-8\beta s+16\beta^2s^2+16\alpha^2s-16\alpha^2s^2$ · *why:* $(1-4\beta s)^2=1-8\beta s+16\beta^2s^2$;
     multiply out the second term. · *plain:* every term written separately.
  5. *did:* Subtract 1 and factor out s · *tex:* $\lvert G\rvert^2-1=s\big[16\alpha^2-8\beta+s(16\beta^2-16\alpha^2)\big]$ · *why:* The constant 1 cancels
     and every remaining term contains s. · *plain:* growth minus one equals s times a straight line in s.
  6. *did:* Handle the longest wave, s = 0 · *tex:* $s=0:\ \lvert G\rvert^2-1=0$ · *why:* θ = 0 is a constant error, G = 1 (D04 check): it never grows,
     so it imposes no condition. · *plain:* a constant offset is harmless.
  7. *did:* Divide by s > 0 · *tex:* $f(s)\equiv16\alpha^2-8\beta+s\,(16\beta^2-16\alpha^2)\le0\quad\text{for all }s\in(0,1]$ · *why:* Dividing an
     inequality by a positive number keeps its direction (s > 0 here). · *plain:* stability ⇔ a straight line stays below zero.
  8. *did:* Use the endpoint rule for a line · *tex:* $f\le0\text{ on }(0,1]\iff f(0^+)\le0\ \text{and}\ f(1)\le0$ · *why:* P226: a linear function takes its
     largest value at an end of the interval; at the open end we need the limit s → 0⁺, i.e. f(0) ≤ 0 by continuity. · *plain:* check the
     two ends only.
  9. *did:* The long-wave end, s → 0⁺ · *tex:* $16\alpha^2-8\beta\le0\iff4\alpha^2\le2\beta$ · *why:* Put s = 0 in f and divide by 4. This is the edge
     where convection kicks long waves harder than diffusion damps them. · *plain:* enough diffusion for the long waves.
  10. *did:* The zigzag end, s = 1 · *tex:* $16\beta^2-8\beta\le0\iff8\beta(2\beta-1)\le0\iff0\le\beta\le\tfrac12$ · *why:* At s = 1 the α² terms cancel
      (16α² − 16α²); a product of two factors is ≤ 0 when they have opposite signs, which with β ≥ 0 gives 2β ≤ 1. · *plain:* diffusion must
      not overshoot the shortest wave.
  11. *did:* Combine the two ends · *tex:* $0\le4\alpha^2\le2\beta\le1$ (10.27) · *why:* 4α² ≥ 0 always; steps 9 and 10 chained; both are needed and
      together sufficient (step 8). · *plain:* Noye's region: a lens in the (α, β) plane.
  12. *did:* Pure diffusion, u = 0 · *tex:* $0\le\beta\le\tfrac12\iff\Delta t\le\tfrac12\frac{\Delta x^2}D$ (10.28) · *why:* α = 0 makes step 9 automatic;
      step 10 remains; β = DΔt/Δx² (10.11) solved for Δt. Ch. 1's r ≤ ½ (R05). · *plain:* the diffusion time-step limit.
  13. *did:* Pure convection, D = 0 · *tex:* $\beta=0:\ \lvert G\rvert^2=1+4\alpha^2\sin^2\theta>1\quad(0<\theta<\pi,\ \alpha\ne0)$ · *why:* With β = 0, step 9
      demands 4α² ≤ 0, i.e. u = 0; directly, (10.26) becomes 1 + (2α sin θ)², larger than 1 for every wave except θ = 0, π. The book's 'never
      satisfied'. · *plain:* FTCS cannot do pure convection at any time step.
- **Result.** (10.27), (10.28), and pure convection unstable — *in words:* FTCS is stable exactly when diffusion is strong enough for the
  long waves (4α² ≤ 2β) and gentle enough for the zigzag (2β ≤ 1).
- **Check.** (α, β) = (0.1, 0.2): 0.04 ≤ 0.4 ≤ 1 ✓ stable; (0.4, 0.2): 0.64 > 0.4 ✗ (the scan finds ∣G∣ > 1 at small θ); (0, 0.51): 1.02 > 1 ✗
  (∣G(π)∣ = 1.04). The brute-force scan over a 60 × 60 (α, β) grid agrees with the closed form (C04 from-scratch cell). In terms of the
  physics: 4α² ≤ 2β ⇔ $u^2\Delta t\le2D$ — the time step must not exceed 2D/u².
- **Optional sympy check (`check_src`):** `expand((1-4*b*s)**2 + 16*a**2*s*(1-s) - 1)` and `factor` show s·(…); substitute s = 1 and
  `solve` the two end conditions.
- **What it means.** Two separate dangers: long waves (convection-driven, 4α² ≤ 2β ⇔ Δt ≤ 2D/u²) and the zigzag (diffusion-driven,
  Δt ≤ Δx²/2D). Refining the grid tightens the second quadratically — the cost of explicit diffusion.
- **Traps.** Dividing by s without excluding s = 0; testing only θ = π (misses the long-wave edge); treating the condition as α ≤ …
  (only α² appears — u of either sign behaves the same for FTCS).

### D07 · BTCS is unconditionally stable: $G=\big[1+4\beta\sin^2\frac\theta2+2\mathrm i\alpha\sin\theta\big]^{-1}$, $\lvert G\rvert\le1$ for every Δt — ★★, 7 steps, in C04 (notebook · `von_neumann_amplification`)
- **Goal.** Show what the book asserts ('it can easily be shown'): the implicit scheme (10.13) never lets an error grow, whatever Δt.
- **Start.** $T^n_i+\alpha(T^n_{i+1}-T^n_{i-1})-\beta(T^n_{i+1}-2T^n_i+T^n_{i-1})\approx T^{n-1}_i$ (10.13) — *in words:* the space differences taken at the
  new level.
- **Plan.** (1) The error obeys the same rule (as D04). (2) Insert one wave. (3) Divide instead of multiply. (4) Bound the denominator.
- **Tools.** D04 steps 1–9 · Euler's formula (P45) · half-angle identity (P225) · modulus of a quotient (P153).
- **Assumptions.** Linear, constant α, β; β ≥ 0 (step 6).
- **Steps.**
  1. *did:* Error equation for BTCS · *tex:* $\xi^n_j+\alpha(\xi^n_{j+1}-\xi^n_{j-1})-\beta(\xi^n_{j+1}-2\xi^n_j+\xi^n_{j-1})=\xi^{n-1}_j$ · *why:* (10.13) is
     linear, so the difference of two solutions obeys it (D04 step 3); index renamed j. · *plain:* the error follows the implicit rule.
  2. *did:* Insert one wave · *tex:* $\hat\xi^n\big[1+\alpha(e^{\mathrm i\theta}-e^{-\mathrm i\theta})-\beta(e^{\mathrm i\theta}-2+e^{-\mathrm i\theta})\big]=\hat\xi^{n-1}$ · *why:* Put
     $\xi^n_j=\hat\xi^ne^{\mathrm i\pi kx_j}$, shift the neighbours by $e^{\pm\mathrm i\theta}$ and divide by $e^{\mathrm i\pi kx_j}$ (D04 steps 5–9). · *plain:* the new
     amplitude times a bracket equals the old amplitude.
  3. *did:* Simplify the convective bracket · *tex:* $e^{\mathrm i\theta}-e^{-\mathrm i\theta}=2\mathrm i\sin\theta$ · *why:* Euler's formula: the cosines cancel,
     the sines add. · *plain:* convection contributes an imaginary part.
  4. *did:* Simplify the diffusive bracket · *tex:* $e^{\mathrm i\theta}-2+e^{-\mathrm i\theta}=2\cos\theta-2=-4\sin^2\tfrac\theta2$ · *why:* Euler's formula,
     then 1 − cos θ = 2 sin²(θ/2) (P225). · *plain:* diffusion contributes a real, positive part (after the minus sign).
  5. *did:* Divide to get G · *tex:* $G=\frac{\hat\xi^n}{\hat\xi^{n-1}}=\frac1{1+4\beta\sin^2\frac\theta2+2\mathrm i\alpha\sin\theta}$ · *why:* Steps 3–4 in step 2; the
     unknowns sit at the new level, so the bracket multiplies $\hat\xi^n$ and we divide by it. · *plain:* the implicit scheme divides instead of
     multiplying.
  6. *did:* Bound the denominator · *tex:* $\big\lvert1+4\beta s+2\mathrm i\alpha\sin\theta\big\rvert\ge1+4\beta s\ge1,\quad s=\sin^2\tfrac\theta2$ · *why:* The modulus
     of a complex number is at least its real part; with β ≥ 0 and s ≥ 0 the real part is at least 1. · *plain:* the denominator is never
     smaller than 1.
  7. *did:* Conclude · *tex:* $\lvert G\rvert=\frac1{\lvert1+4\beta s+2\mathrm i\alpha\sin\theta\rvert}\le1\quad\text{for every }\Delta t$ · *why:* The modulus of a
     quotient is the quotient of moduli; a number ≥ 1 in the denominator gives ≤ 1. No condition on α or β remains. · *plain:* no wave can
     grow — unconditional stability.
- **Result.** $\lvert G_{BTCS}\rvert\le1$ for all θ, α and β ≥ 0 — *in words:* BTCS is stable for any time step (it is consistent too, by the
  moves of D03 with a backward time difference).
- **Check.** β = 100, θ = π, α = 0: G = 1/401 = 0.00249 ✓ (N17 code). α = β = 0: G = 1 (nothing happens) ✓.
- **What it means.** Implicit schemes buy freedom in Δt with a linear solve per step (N09); accuracy still limits Δt (the smeared puff
  of C02's figure at β = 2.5). Climate models treat their fastest waves this way.
- **Traps.** Multiplying by the bracket instead of dividing (the unknown level is on the left); forgetting that ∣z∣ ≥ Re z needs the real
  part positive (β ≥ 0).

### D08 · Upwind and the CFL condition: $\lvert G\rvert^2=1-2C(1-C)(1-\cos\theta)$, stable iff $u\frac{\Delta t}{\Delta x}\le1$ (10.30) — ★★, 9 steps, in C05 (notebook · `upwind_cfl_advection`)
- **Goal.** Find when the upwind scheme (10.29) is stable, and see why the answer is 'the characteristic must stay inside the stencil'.
- **Start.** $T^{n+1}_i=T^n_i-2\alpha(T^n_i-T^n_{i-1})$ (10.29), u > 0 — *in words:* the convective difference taken from the upstream side.
- **Plan.** (1) Write the scheme with the Courant number. (2) Find G (D04 moves). (3) Compute ∣G∣² and read its sign. (4) Special cases
  and the sign of u.
- **Tools.** D04 moves · Euler's formula (P45) · cos² + sin² = 1 · domain of dependence (primer P227).
- **Assumptions.** Pure convection with constant u > 0 (step 1); periodic grid.
- **Steps.**
  1. *did:* Name the Courant number · *tex:* $C=2\alpha=u\frac{\Delta t}{\Delta x},\qquad T^{n+1}_j=T^n_j-C(T^n_j-T^n_{j-1})$ · *why:* The book's α carries a ½
     from the centred difference; the upwind difference has none, so 2α = C appears. Index renamed j. · *plain:* C is the number of cells the
     flow moves per step.
  2. *did:* Insert one wave · *tex:* $G=1-C\big(1-e^{-\mathrm i\theta}\big)$ · *why:* The error obeys the linear rule (D04 step 3); with
     $\xi_j=\hat\xi e^{\mathrm i\pi kx_j}$ the upstream neighbour is $e^{-\mathrm i\theta}$ times it; divide out the wave (D04 step 9). · *plain:* one complex factor
     per wave.
  3. *did:* Split real and imaginary parts · *tex:* $G=(1-C+C\cos\theta)-\mathrm iC\sin\theta$ · *why:* $e^{-\mathrm i\theta}=\cos\theta-\mathrm i\sin\theta$ (P45). · *plain:*
     real part: what is kept; imaginary part: the turn.
  4. *did:* Square and add · *tex:* $\lvert G\rvert^2=(1-C+C\cos\theta)^2+C^2\sin^2\theta$ · *why:* ∣a + ib∣² = a² + b² (P153). · *plain:* the squared factor.
  5. *did:* Expand and use cos² + sin² = 1 · *tex:* $\lvert G\rvert^2=(1-C)^2+2C(1-C)\cos\theta+C^2$ · *why:* $(1-C+C\cos\theta)^2=(1-C)^2+2C(1-C)
     \cos\theta+C^2\cos^2\theta$ and $C^2\cos^2\theta+C^2\sin^2\theta=C^2$. · *plain:* three terms left.
  6. *did:* Regroup around 1 · *tex:* $\lvert G\rvert^2=1-2C(1-C)(1-\cos\theta)$ · *why:* $(1-C)^2+C^2=1-2C(1-C)$; the remaining terms share the
     factor 2C(1 − C). · *plain:* one minus a product.
  7. *did:* Read the sign · *tex:* $\lvert G\rvert^2\le1\ \forall\theta\iff C(1-C)\ge0\iff0\le C\le1$ (10.30) · *why:* 1 − cos θ lies in [0, 2] and is positive for
     θ ≠ 0, so ∣G∣² ≤ 1 for all waves exactly when C(1 − C) ≥ 0. · *plain:* stable iff the flow moves at most one cell per step: the CFL
     condition.
  8. *did:* Look at the edges · *tex:* $C=1:\ G=e^{-\mathrm i\theta};\qquad C>1:\ \lvert G(\pi)\rvert^2=(2C-1)^2>1$ · *why:* C = 1 turns the rule into
     $T^{n+1}_j=T^n_{j-1}$, an exact shift; for C > 1 the zigzag grows (C = 1.1: ∣G∣ = 1.2). Physically (P227): the foot of the characteristic,
     $x_j-u\Delta t$, lies C cells upstream and must stay within the stencil $[x_{j-1},x_j]$. · *plain:* on the edge exact, beyond it explosive.
  9. *did:* Let u be negative · *tex:* $u<0:\quad T^{n+1}_j=T^n_j-C(T^n_{j+1}-T^n_j),\ \ \lvert C\rvert=\frac{\lvert u\rvert\Delta t}{\Delta x}\le1$ · *why:* With the printed
     stencil and u < 0, C < 0 and C(1 − C) < 0: unstable for every step. The upstream side is now j + 1 (slip R11); the same algebra gives
     ∣C∣ ≤ 1. · *plain:* take the difference from where the flow comes from, whichever way that is.
- **Result.** $\lvert G\rvert^2=1-2C(1-C)(1-\cos\theta)$, stable iff 0 ≤ C = uΔt/Δx ≤ 1 (10.30), with the side chosen by the sign of u —
  *in words:* a fluid particle must not cross more than one cell per step.
- **Check.** C dimensionless ✓. C = ½, θ = π: ∣G∣² = 1 − 2·½·½·2 = 0 (the zigzag killed in one step) ✓ (tiny example). The numerical diffusivity
  of the stable scheme, ∣u∣Δx(1 − C)/2, vanishes at C = 1 — consistent with the exact shift.
- **What it means.** The CFL condition sets the time step of every explicit weather and ocean model (with the gravity- or sound-wave
  speed, C05's climate note) and of MacCormack (D18, with the sound speed).
- **Traps.** Writing C = α (the ½ belongs to the centred difference); thinking C = 1 is 'marginal' (it is exact); using the upwind side
  of (10.29) when u < 0 (slip R11).

### D09 · The weak form $\int_0^LT_tw\,dx+u\int_0^LT_xw\,dx+D\int_0^LT_xw_x\,dx=Dq\,w(L)$ (10.34)–(10.36) — ★★, 9 steps, in C06 (notebook · `fem_hat_assembly`)
- **Goal.** Rewrite the transport problem so that only first derivatives of T appear and the flux condition at x = L is built in.
- **Start.** The strong form: $T_t+uT_x=DT_{xx}$ on (0, L) (10.1) with $T(0,t)=g$ and $T_x(L,t)=q$ (10.2) — *in words:* the equation at every point
  plus its two end conditions.
- **Plan.** (1) Multiply by a test function that vanishes where T is prescribed. (2) Integrate over the rod. (3) Move one derivative from
  T onto w (integration by parts). (4) Insert the end conditions.
- **Tools.** Test functions and the spaces H¹, 𝒮, V (primer P230) · product rule (Ch. 1 P38) · integration by parts (Ch. 9 P218a) · the
  fundamental theorem of calculus (Ch. 2 P84).
- **Assumptions.** T and w in H¹ (square-integrable slopes); w(0) = 0 (w ∈ V, used in step 7); T satisfies the Neumann datum (step 8).
- **Steps.**
  1. *did:* Multiply the PDE by w · *tex:* $\big(T_t+uT_x-DT_{xx}\big)\,w=0$ · *why:* A true equation stays true when both sides are multiplied by
     any function; w ∈ V is our 'test'. · *plain:* weight the equation's residual by w.
  2. *did:* Integrate over the rod · *tex:* $\int_0^LT_tw\,dx+u\int_0^LT_xw\,dx=D\int_0^LT_{xx}w\,dx$ (10.34) · *why:* An integral of zero is zero;
     constants u, D come out of the integrals; the diffusion term is moved to the right as the book writes it. · *plain:* the weighted
     residual vanishes on average.
  3. *did:* Product rule for $(T_xw)_x$ · *tex:* $\frac{d}{dx}\big(T_xw\big)=T_{xx}w+T_xw_x$ · *why:* P38; it links the unwanted $T_{xx}w$ to a total
     derivative and a term with only first derivatives. · *plain:* the second derivative can be traded.
  4. *did:* Solve for $T_{xx}w$ · *tex:* $T_{xx}w=\frac{d}{dx}\big(T_xw\big)-T_xw_x$ · *why:* Subtract $T_xw_x$ from both sides of step 3. · *plain:* curvature
     times w = a derivative of something minus slope times slope.
  5. *did:* Integrate: integration by parts · *tex:* $\int_0^LT_{xx}w\,dx=\big[T_xw\big]_0^L-\int_0^LT_xw_x\,dx$ · *why:* The integral of a derivative is
     the difference of end values (P84) — together this is integration by parts (P218a). · *plain:* one derivative moved from T to w, at the
     price of an end term.
  6. *did:* Substitute into (10.34) · *tex:* $\int_0^LT_tw\,dx+u\int_0^LT_xw\,dx+D\int_0^LT_xw_x\,dx=D\big[T_xw\big]_0^L$ · *why:* Step 5 times D, the
     integral term moved to the left (its sign flips from − to +). · *plain:* only first derivatives of T remain.
  7. *did:* Use w(0) = 0 at the left end · *tex:* $D\big[T_xw\big]_0^L=D\,T_x(L)\,w(L)-D\,T_x(0)\cdot0$ · *why:* Test functions vanish where T is
     prescribed (w ∈ V); so the unknown flux at x = 0 never enters. · *plain:* the Dirichlet end drops out.
  8. *did:* Use $T_x(L)=q$ at the right end · *tex:* $\int_0^LT_tw\,dx+u\int_0^LT_xw\,dx+D\int_0^LT_xw_x\,dx=Dq\,w(L)$ (10.35) · *why:* The Neumann condition
     (10.2) supplies the slope at L; the flux datum now sits on the right-hand side. · *plain:* the flux condition becomes a source term.
  9. *did:* State the weak problem · *tex:* $\text{find }T\in\mathcal S:\ \text{(10.35) holds for all }w\in V$ (10.36) · *why:* 'For all w' replaces 'at every x';
     $\mathcal S$ carries T(0) = g (P230). · *plain:* the weak form: find the T that passes every test.
- **Result.** (10.36) — *in words:* the equation holds on average against every admissible test function; only first derivatives appear,
  the Neumann datum appears as Dq w(L) and the Dirichlet datum lives in the trial space.
- **Check.** Units: every term is (T/s)·(w)·m ✓ (Dq w(L): m²/s · T/m · w). Tiny example T = x², w = x on [0, 1]: $\int T_{xx}w=1$ and
  $[T_xw]-\int T_xw_x=2-1=1$ ✓. `ch10.weak_form_sympy()["residual"]` = 0.
- **What it means.** The starting line of every finite-element method (C07, C13). It fails only for trial functions with jumps (not in H¹).
- **Traps.** The sign after moving the derivative (−∫T_x w_x on the right becomes +D∫T_x w_x on the left); forgetting w(0) = 0 (then the unknown
  flux at x = 0 appears); putting the Dirichlet value into the weak form instead of the space.

### D10 · Weak ⇒ strong: $\int_0^L(T_t+uT_x-DT_{xx})w\,dx+D\big[T_x(L)-q\big]w(L)=0$ (10.37), hence the PDE and $T_x(L)=q$ (10.38) — ★★, 9 steps, in C06 (notebook)
- **Goal.** Show that a (smooth) solution of the weak problem solves the PDE and satisfies the flux condition — nothing was lost.
- **Start.** (10.36) for T ∈ 𝒮 and every w ∈ V — *in words:* T passes every test.
- **Plan.** (1) Undo the integration by parts. (2) Test with functions that vanish at both ends → the PDE (fundamental lemma). (3) Test
  with one function that does not vanish at L → the boundary condition.
- **Tools.** Integration by parts (P218a) · the fundamental lemma of the calculus of variations (primer P231).
- **Assumptions.** T smooth enough for $T_{xx}$ to exist (T ∈ H², step 2); D > 0 (step 8).
- **Steps.**
  1. *did:* Start from the weak form · *tex:* $\int_0^LT_tw\,dx+u\int_0^LT_xw\,dx+D\int_0^LT_xw_x\,dx=Dq\,w(L)$ · *why:* (10.36), assumed to hold for every
     w ∈ V. · *plain:* the tests T passes.
  2. *did:* Integrate the diffusion term back by parts · *tex:* $D\int_0^LT_xw_x\,dx=D\big[T_xw\big]_0^L-D\int_0^LT_{xx}w\,dx$ · *why:* D09 step 5 read the
     other way (P218a); allowed because T is now assumed twice differentiable. · *plain:* move the derivative back onto T.
  3. *did:* Substitute and use w(0) = 0 · *tex:* $\int_0^L(T_t+uT_x-DT_{xx})w\,dx+D\,T_x(L)\,w(L)=Dq\,w(L)$ · *why:* The lower end term vanishes because
     w(0) = 0; collect the three integrals into one. · *plain:* the residual of the PDE appears inside the integral.
  4. *did:* Move Dq w(L) left · *tex:* $\int_0^L(T_t+uT_x-DT_{xx})w\,dx+D\big[T_x(L)-q\big]w(L)=0$ (10.37) · *why:* Subtract the right side; this is the
     book's (10.37). · *plain:* two pieces must add to zero: the PDE residual and the flux mismatch.
  5. *did:* Test with w vanishing at L too · *tex:* $\int_0^L(T_t+uT_x-DT_{xx})\,w\,dx=0\quad\text{for all }w\text{ with }w(0)=w(L)=0$ · *why:* Such w belong to
     V, so (10.37) applies to them, and the boundary term vanishes. · *plain:* the residual is orthogonal to every interior bump.
  6. *did:* Apply the fundamental lemma · *tex:* $T_t+uT_x-DT_{xx}=0\quad\text{on }(0,L)$ · *why:* P231: if ∫f w dx = 0 for every such w and f is
     continuous, f = 0. · *plain:* the PDE holds at every interior point.
  7. *did:* Return to (10.37) with any w · *tex:* $D\big[T_x(L)-q\big]w(L)=0\quad\text{for all }w\in V$ · *why:* By step 6 the integral in (10.37) is now zero
     for every w, leaving only the boundary term. · *plain:* only the flux mismatch is left.
  8. *did:* Choose w with w(L) = 1 · *tex:* $T_x(L)-q=0$ (10.38) · *why:* w = x/L is in V (w(0) = 0) with w(L) = 1; D > 0 can be divided out. ·
     *plain:* the flux condition holds — it came out by itself.
  9. *did:* Name the two kinds of condition · *tex:* $T(0)=g:\ \text{essential (in }\mathcal S);\qquad T_x(L)=q:\ \text{natural (from the weak form)}$ · *why:* The
     Dirichlet value was never derived — it was imposed by choosing T ∈ 𝒮; the Neumann one followed from the equations. · *plain:*
     essential conditions are built in, natural ones come free.
- **Result.** (10.37)–(10.38) — *in words:* weak and strong forms are equivalent for smooth T; the Neumann condition is natural, the
  Dirichlet one essential.
- **Check.** `ch10.weak_to_strong_sympy()` returns the PDE and $T_x(L)-q$; the FE figure of C06 shows the natural slope approached at order 1
  while T(0) is exact on every mesh.
- **What it means.** In FE codes, flux (Neumann) conditions cost nothing — they are the default at any boundary where nothing is
  imposed; value (Dirichlet) conditions must be enforced (row replacement or penalty, N109).
- **Traps.** Choosing a w that is nonzero at x = 0 (not in V); expecting the Dirichlet condition to come out; needing $T_{xx}$ without saying
  T is smooth enough.

### D11 · Galerkin to matrix form $\mathbf M\dot{\mathbf d}+\mathbf K\mathbf d=\mathbf F$ (10.58) with $M_{AB}=\int_0^LN_AN_B\,dx$ (10.55), $K_{AB}=u\int_0^LN_{B,x}N_A\,dx+D\int_0^LN_{B,x}N_{A,x}\,dx$ (10.56), $F_A=DqN_A(L)-gu\int_0^LN_{0,x}N_A\,dx-gD\int_0^LN_{0,x}N_{A,x}\,dx$ (10.57) — ★★, 13 steps, in C07 (notebook · `fem_hat_assembly`)
- **Goal.** Turn the weak form, restricted to finitely many basis functions, into a system of ODEs for the nodal weights.
- **Start.** The weak form in finite spaces: $\int_0^LT^h_tw^h\,dx+u\int_0^LT^h_xw^h\,dx+D\int_0^LT^h_xw^h_x\,dx=Dq\,w^h(L)$ for all $w^h\in V^h$ (10.39) —
  *in words:* the weak problem with a finite menu of trial and test functions.
- **Plan.** (1) Lift off the boundary value so trial and test functions share one space. (2) Write the bilinear form. (3) Insert the
  expansions and pull the sums out. (4) Arbitrary test coefficients ⇒ one equation per basis function. (5) Name the matrices.
- **Tools.** Bilinear form (primer P232) · linear combinations (Ch. 2) · arbitrary coefficients ⇒ each bracket zero (primer P233).
- **Assumptions.** $\mathcal S^h\subset\mathcal S$, $V^h\subset V$ (step 1); the boundary value g constant in time (step 3); the basis $N_A$ with
  $N_A(0)=0$ and one extra $N_0$ (N25–N27).
- **Steps.**
  1. *did:* Split off the boundary value · *tex:* $v^h=T^h-g^h,\qquad g^h(0)=g$ (10.40) · *why:* Then $v^h(0)=0$, so $v^h$ lives in $V^h$ like the test
     functions — the Galerkin requirement (trial = test space). · *plain:* the unknown part vanishes where T is prescribed.
  2. *did:* Name the bilinear form · *tex:* $a(w,v)=u\int_0^Lv_xw\,dx+D\int_0^Lv_xw_x\,dx$ (10.42) · *why:* The two space terms of (10.39) together,
     linear in each slot (P232); a shorthand for everything that is not the time derivative. · *plain:* one symbol for convection plus
     diffusion.
  3. *did:* Use $g^h$ constant in time · *tex:* $T^h_t=v^h_t$ · *why:* g is a fixed number, so $g^h(x)=gN_0(x)$ does not depend on t. · *plain:* only the
     unknown part changes in time.
  4. *did:* Split a(·,·) over the sum · *tex:* $a(w^h,T^h)=a(w^h,v^h)+a(w^h,g^h)$ · *why:* Linearity in the second slot (P232). · *plain:* convection and
     diffusion act on each part separately.
  5. *did:* Write the Galerkin form · *tex:* $\int_0^Lv^h_tw^h\,dx+a(w^h,v^h)=Dq\,w^h(L)-a(w^h,g^h)$ (10.41) · *why:* Steps 2–4 in (10.39), the known
     $a(w^h,g^h)$ moved to the right. · *plain:* unknowns left, known data right.
  6. *did:* Insert the expansions · *tex:* $w^h=\sum_Ac_AN_A,\ \ v^h=\sum_Bd_B(t)N_B,\ \ g^h=gN_0$ (10.45), (10.48), (10.47) · *why:* Every function in
     the finite spaces is a combination of the basis (N25–N27); the time dependence sits in the weights $d_B(t)$. · *plain:* functions
     become lists of numbers.
  7. *did:* Substitute into (10.41) · *tex:* $\int_0^L\Big(\sum_B\dot d_BN_B\sum_Ac_AN_A\Big)dx+a\Big(\sum_Ac_AN_A,\sum_Bd_BN_B\Big)=Dq\sum_Ac_AN_A(L)-a\Big(\sum_Ac_AN_A,gN_0\Big)$ (10.50) ·
     *why:* Straight substitution; $\dot d_B=dd_B/dt$ because $v^h_t=\sum\dot d_BN_B$. · *plain:* the Galerkin form in terms of the weights.
  8. *did:* Pull the sums out of the integral · *tex:* $\int_0^L\sum_A\sum_Bc_A\dot d_BN_AN_B\,dx=\sum_Ac_A\sum_B\dot d_B\int_0^LN_AN_B\,dx$ · *why:* The integral of a
     finite sum is the sum of integrals, and $c_A$, $\dot d_B$ do not depend on x. · *plain:* the mass integrals separate from the numbers.
  9. *did:* Pull the sums out of a(·,·) · *tex:* $a\Big(\sum_Ac_AN_A,\sum_Bd_BN_B\Big)=\sum_Ac_A\sum_Bd_B\,a(N_A,N_B)$ · *why:* Bilinearity (P232), in each slot in
     turn; likewise $a(\sum c_AN_A,gN_0)=\sum_Ac_A\,g\,a(N_A,N_0)$. · *plain:* the stiffness integrals separate too.
  10. *did:* Collect everything under one sum · *tex:* $\sum_{A=1}^nc_AG_A=0$ (10.51), $\ G_A=\sum_B\dot d_B\int_0^LN_AN_B\,dx+\sum_Bd_B\,a(N_A,N_B)-DqN_A(L)+g\,a(N_A,N_0)$ (10.52) ·
      *why:* Move all terms to the left and factor $c_A$ out of each. · *plain:* a weighted sum of brackets must vanish.
  11. *did:* Arbitrary $c_A$ ⇒ each bracket zero · *tex:* $\sum_B\dot d_B\int_0^LN_BN_A\,dx+\sum_Bd_B\,a(N_A,N_B)=DqN_A(L)-g\,a(N_A,N_0),\ \ A=1,\dots,n$ (10.53) · *why:* (10.51) must
      hold for every test function, i.e. every choice of the $c_A$; taking unit vectors gives $G_A=0$ one by one (P233). · *plain:* one ODE
      per basis function — as many equations as unknowns.
  12. *did:* Name the matrices · *tex:* $M_{AB}=\int_0^LN_AN_B\,dx,\ \ K_{AB}=a(N_A,N_B),\ \ F_A=DqN_A(L)-g\,a(N_A,N_0)$ (10.55)–(10.57) · *why:* Read the
      coefficients of $\dot d_B$, of $d_B$ and the known right side off (10.53); expanding a(·,·) gives the book's (10.56) and (10.57). M is
      symmetric, K is not (convection differentiates only $N_B$). · *plain:* mass, stiffness and force.
  13. *did:* Write it as one matrix equation · *tex:* $\mathbf M\dot{\mathbf d}+\mathbf K\mathbf d=\mathbf F$ (10.58) · *why:* The n equations of (10.53) are the
      rows of a matrix–vector product (10.54). · *plain:* a system of ODEs, ready for any time stepper.
- **Result.** (10.58) with (10.55)–(10.57) — *in words:* Galerkin turns the PDE into mass × rate + stiffness × weights = force.
- **Check.** `ch10.galerkin_equations_sympy(3)` returns M, K, F for three hats; units: $M_{AB}$ [m], $K_{AB}$ [m/s], $F_A$ [T·m/s] — each row of (10.58)
  in T·m/s ✓. With u = 0, K is symmetric (pure diffusion) ✓.
- **What it means.** Any basis gives this form; hat functions make M and K tridiagonal (D12). The same recipe with P2/P1 triangles gives
  (10.137) in C13.
- **Traps.** Forgetting that K differentiates the *trial* function in the convective part (swapping A and B makes K the transpose);
  leaving the g N₀ terms on the left (they are known data); treating $\dot d_B$ as time derivatives of the basis.

### D12 · The interior row $\frac d{dt}\Big(\frac{T_{A-1}}6+\frac{2T_A}3+\frac{T_{A+1}}6\Big)+\frac u{2h}(T_{A+1}-T_{A-1})-\frac D{h^2}(T_{A-1}-2T_A+T_{A+1})=0$ (10.63) — ★★, 10 steps, in C07 (notebook · `fem_hat_assembly`)
- **Goal.** Evaluate the integrals of row A of (10.53) for hat functions on a uniform mesh and compare with finite differences.
- **Start.** Row A of (10.53) for an interior node, $\sum_B\dot d_B\int_0^LN_BN_A\,dx+\sum_Bd_B\,a(N_A,N_B)=0$ ($N_A(L)=0$ and $N_A$ does not overlap $N_0$ for
  2 ≤ A ≤ n − 2) — *in words:* one row of the Galerkin system.
- **Plan.** (1) Only three neighbours overlap. (2) Compute the three mass integrals, the three convective and the three diffusive ones.
  (3) Assemble and divide by h.
- **Tools.** Hat functions (10.59) (N30) · integrals of linear functions · D11.
- **Assumptions.** Uniform mesh $x_A=Ah$ (step 2); interior node away from both ends (the right side is zero).
- **Steps.**
  1. *did:* Keep only overlapping neighbours · *tex:* $B\in\{A-1,A,A+1\}$ · *why:* $N_A$ is nonzero only on $[x_{A-1},x_{A+1}]$ (compact support, N30);
     every other product $N_AN_B$ and $N_{A,x}N_{B,x}$ is zero everywhere. · *plain:* three terms per row.
  2. *did:* Mass, diagonal · *tex:* $\int N_A^2\,dx=2\int_0^h\big(\tfrac sh\big)^2ds=\tfrac{2h}3$ · *why:* On each of its two elements $N_A$ is a straight ramp from
     0 to 1; ∫₀ʰ(s/h)²ds = h/3, twice. · *plain:* a hat overlaps itself on two elements.
  3. *did:* Mass, off-diagonal · *tex:* $\int N_AN_{A\pm1}\,dx=\int_0^h\big(1-\tfrac sh\big)\tfrac sh\,ds=\tfrac h6$ · *why:* The neighbours overlap on one element,
     where one rises and the other falls; ∫₀¹(1 − σ)σ dσ = 1/6. · *plain:* neighbours share one element.
  4. *did:* Slopes of the hats · *tex:* $N_{A,x}=+\tfrac1h$ on $[x_{A-1},x_A]$, $-\tfrac1h$ on $[x_A,x_{A+1}]$ · *why:* Differentiate (10.59); likewise
     $N_{A+1,x}=+1/h$ on $[x_A,x_{A+1}]$ and $N_{A-1,x}=-1/h$ on $[x_{A-1},x_A]$. · *plain:* up-slope then down-slope.
  5. *did:* Convective integrals · *tex:* $\int N_{A+1,x}N_A\,dx=\tfrac1h\cdot\tfrac h2=\tfrac12,\quad\int N_{A-1,x}N_A\,dx=-\tfrac12,\quad\int N_{A,x}N_A\,dx=0$ · *why:*
     Constant slope times the area under half a hat (h/2); the last is $\frac12[N_A^2]$ between two zeros. · *plain:* convection couples only the
     two neighbours, with opposite signs.
  6. *did:* Diffusive integrals · *tex:* $\int N_{A,x}^2\,dx=\tfrac2h,\qquad\int N_{A\pm1,x}N_{A,x}\,dx=-\tfrac1h$ · *why:* (1/h)² over length 2h; for neighbours
     (±1/h)(∓1/h) over one element h. · *plain:* the diffusion weights 2/h and −1/h.
  7. *did:* Assemble the mass part · *tex:* $h\Big(\tfrac16\dot d_{A-1}+\tfrac23\dot d_A+\tfrac16\dot d_{A+1}\Big)$ · *why:* Steps 2–3 times $\dot d_B$. · *plain:* a
     weighted average of three rates.
  8. *did:* Assemble the stiffness part · *tex:* $\tfrac u2(d_{A+1}-d_{A-1})+\tfrac Dh(-d_{A-1}+2d_A-d_{A+1})$ · *why:* $K_{AB}=u\int N_{B,x}N_A+D\int N_{B,x}N_{A,x}$
     (10.56) with steps 5–6. · *plain:* centred convection and diffusion, times h.
  9. *did:* Set the row to zero · *tex:* $h\Big(\tfrac16\dot d_{A-1}+\tfrac23\dot d_A+\tfrac16\dot d_{A+1}\Big)+\tfrac u2(d_{A+1}-d_{A-1})-\tfrac Dh(d_{A-1}-2d_A+d_{A+1})=0$ ·
     *why:* $F_A=0$ for an interior node ($N_A(L)=0$, no overlap with $N_0$). · *plain:* the row of (10.53).
  10. *did:* Divide by h, use $d_A=T_A$ · *tex:* $\frac d{dt}\Big(\frac{T_{A-1}}6+\frac{2T_A}3+\frac{T_{A+1}}6\Big)+\frac u{2h}(T_{A+1}-T_{A-1})-\frac D{h^2}(T_{A-1}-2T_A+T_{A+1})=0$ (10.63) ·
      *why:* The weights are nodal values (10.62) (N31); dividing by h makes the stencils look like (10.6)–(10.7). · *plain:* centred FD in
      space, with a ⅙–⅔–⅙ averaged time derivative.
- **Result.** (10.63) — *in words:* uniform linear FE = centred FD for convection and diffusion, with a consistent-mass average in time.
- **Check.** The mass weights sum to 1 (⅙ + ⅔ + ⅙) ✓; the stiffness row sums to 0 (a constant T has no convection or diffusion) ✓. Numbers
  (h = 0.25, u = 1, D = 0.25): stiffness row × h = (−1.5, 2, −0.5) ✓ (`FEM1.interior_stencil`); sympy confirms ⅙, ⅔ (analysis note).
- **What it means.** For steady problems FE and FD give identical equations (N33); in time the consistent mass is slightly less
  diffusive. On nonuniform or quadratic meshes the rows differ.
- **Traps.** ∫N_A² = h/3 (forgetting the second element); dividing by h too early (the mass row is then h/6, 2h/3 … mixed up); a sign in
  the convective pair.

### D13 · Parent element, the map (10.65) and the slopes (10.67)–(10.68), with the printed index shift corrected — ★, 7 steps, in C08 (notebook)
- **Goal.** Compute shape-function slopes on any element from one standard 'parent' element.
- **Start.** $N_1(\xi)=\frac12(1-\xi)$, $N_2(\xi)=\frac12(1+\xi)$ on ξ ∈ [−1, 1] (10.64) — *in words:* two straight ramps on a standard interval.
- **Plan.** (1) Map the parent element onto $[x_{A-1},x_A]$. (2) Invert the map. (3) Chain rule for the slopes. (4) Say which global hat
  each local ramp is.
- **Tools.** Affine map to a parent element (primer P234) · chain rule (Ch. 1 P49).
- **Assumptions.** A straight (affine) element of length $h^e=x_A-x_{A-1}>0$.
- **Steps.**
  1. *did:* Build the map from the shapes · *tex:* $x(\xi)=N_1(\xi)x^e_1+N_2(\xi)x^e_2,\quad x^e_1=x_{A-1},\ x^e_2=x_A$ · *why:* The same ramps that
     interpolate T interpolate the coordinate (isoparametric idea); ξ = −1 ↦ $x_{A-1}$, ξ = 1 ↦ $x_A$. · *plain:* the parent interval
     stretched onto the element.
  2. *did:* Multiply out · *tex:* $x(\xi)=\tfrac12\big[(x_A-x_{A-1})\xi+x_A+x_{A-1}\big]$ (10.65) · *why:* Insert (10.64) and collect the ξ terms. · *plain:*
     stretch by $h^e/2$, shift to the midpoint.
  3. *did:* Invert the map · *tex:* $\xi(x)=\frac{2x-x_A-x_{A-1}}{x_A-x_{A-1}}$ (10.66) · *why:* Solve the linear equation (10.65) for ξ. · *plain:* where in
     the parent a point x sits.
  4. *did:* Differentiate the inverse · *tex:* $\frac{d\xi}{dx}=\frac2{h^e}$ · *why:* (10.66) is linear in x with slope $2/(x_A-x_{A-1})$ (P234: dx = (h/2)dξ). ·
     *plain:* one unit of x is 2/h units of ξ.
  5. *did:* Parent slopes · *tex:* $\frac{dN_1}{d\xi}=-\tfrac12,\qquad\frac{dN_2}{d\xi}=+\tfrac12$ · *why:* Differentiate (10.64). · *plain:* one ramp falls, one
     rises.
  6. *did:* Chain rule · *tex:* $\frac{dN_1}{dx}=-\tfrac12\cdot\tfrac2{h^e}=-\frac1{h^e},\qquad\frac{dN_2}{dx}=+\frac1{h^e}$ · *why:* $\frac{dN}{dx}=\frac{dN}{d\xi}
     \frac{d\xi}{dx}$ (P49). · *plain:* slopes ∓1/h on the element.
  7. *did:* Name the global hats · *tex:* $\frac{dN_{A-1}}{dx}=-\frac1{x_A-x_{A-1}},\qquad\frac{dN_A}{dx}=+\frac1{x_A-x_{A-1}}\quad\text{on }[x_{A-1},x_A]$ · *why:* On this
     element $N_1\leftrightarrow N_{A-1}$ and $N_2\leftrightarrow N_A$ (the book says so above (10.64)). The book prints (10.67)–(10.68) with the labels
     $N_A$ and $N_{A+1}$ — an index shift (slip R1); $N_{A+1}$ is zero here. · *plain:* the falling ramp is the left node's hat, the rising
     one the right node's.
- **Result.** On $[x_{A-1},x_A]$: $\frac{dN_{A-1}}{dx}=-\frac1{h^e}$, $\frac{dN_A}{dx}=+\frac1{h^e}$ — *in words:* every element's slopes come from the
  parent's ±½ times 2/h.
- **Check.** Units 1/m ✓; the two slopes sum to zero (the hats sum to 1 on the element) ✓; `FEM1.shape_slopes(0.5, 0.75)` → (−4, +4) with labels
  (A − 1, A); the `printed=True` labels fail the chain-rule test.
- **What it means.** Every element integral can now be done on [−1, 1] (D14); in 2-D the same idea gives the isoparametric triangles of
  (10.184).
- **Traps.** dξ/dx = h/2 (it is 2/h — that is dx/dξ); copying the printed labels A, A + 1.

### D14 · Element matrices $\mathbf m^e=\frac h6\begin{pmatrix}2&1\\1&2\end{pmatrix}$, $\mathbf k^e=\frac u2\begin{pmatrix}-1&1\\-1&1\end{pmatrix}+\frac Dh\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$, the force (10.77) and assembly into (10.63) — ★★, 11 steps, in C08 (notebook · `fem_hat_assembly`)
- **Goal.** Write out the numbers the book never shows — the 2 × 2 element blocks of (10.75)–(10.77) — and check that adding them up gives
  the interior row (10.63).
- **Start.** $m^e_{ab}=\int_{\Omega^e}N_aN_b\,dx$ (10.75), $k^e_{ab}=u\int_{\Omega^e}N_{b,x}N_a\,dx+D\int_{\Omega^e}N_{b,x}N_{a,x}\,dx$ (10.76), $f^e_a$ (10.77) — *in
  words:* the element versions of M, K, F.
- **Plan.** (1) Change variables to the parent element. (2) Integrate the mass, convective and diffusive products. (3) Place the boundary
  data. (4) Scatter-add and compare with D12.
- **Tools.** Affine map dx = (h/2)dξ (P234) and substitution (Ch. 3 P106) · D13 slopes · scatter-add (primer P235).
- **Assumptions.** Linear elements of length h; node 0 Dirichlet (value g), node n Neumann (datum q).
- **Steps.**
  1. *did:* Change variables · *tex:* $\int_{\Omega^e}f\,dx=\frac h2\int_{-1}^1f\,d\xi$ · *why:* dx = (h/2)dξ on the affine element (P234, P106). · *plain:* every
     element integral becomes one on [−1, 1].
  2. *did:* Mass, diagonal · *tex:* $m^e_{11}=\frac h2\int_{-1}^1\tfrac14(1-\xi)^2d\xi=\frac h2\cdot\frac14\cdot\frac83=\frac h3$ · *why:* ∫₋₁¹(1 − ξ)²dξ = 8/3; the
     same for $m^e_{22}$ by symmetry. · *plain:* each ramp overlaps itself: h/3.
  3. *did:* Mass, off-diagonal · *tex:* $m^e_{12}=\frac h2\int_{-1}^1\tfrac14(1-\xi^2)d\xi=\frac h2\cdot\frac14\cdot\frac43=\frac h6$ · *why:* $(1-\xi)(1+\xi)=1-\xi^2$;
     ∫₋₁¹(1 − ξ²)dξ = 4/3. · *plain:* the two ramps overlap: h/6.
  4. *did:* Collect the mass block · *tex:* $\mathbf m^e=\frac h6\begin{pmatrix}2&1\\1&2\end{pmatrix}$ · *why:* Steps 2–3 (h/3 = 2h/6); symmetric because $N_aN_b=N_bN_a$. ·
     *plain:* the element mass matrix.
  5. *did:* Area under one ramp · *tex:* $\int_{\Omega^e}N_a\,dx=\frac h2\int_{-1}^1\tfrac12(1\mp\xi)\,d\xi=\frac h2$ · *why:* Each ramp is a triangle of height 1 and
     base h. · *plain:* half the element length.
  6. *did:* Convective block · *tex:* $u\int N_{b,x}N_a\,dx=u\cdot\frac{dN_b}{dx}\cdot\frac h2\ \Rightarrow\ \frac u2\begin{pmatrix}-1&1\\-1&1\end{pmatrix}$ · *why:* The slope
     of $N_b$ is constant (∓1/h, D13) and comes out; row a, column b: every row gets (−u/2, +u/2). Not symmetric. · *plain:* convection's block
     depends only on the column.
  7. *did:* Diffusive block · *tex:* $D\int N_{b,x}N_{a,x}\,dx=D\,(\mp\tfrac1h)(\mp\tfrac1h)\,h\ \Rightarrow\ \frac Dh\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$ · *why:* Product of two
     constant slopes times the length h; equal slopes give +, opposite −. · *plain:* diffusion's block is symmetric.
  8. *did:* Add them · *tex:* $\mathbf k^e=\frac u2\begin{pmatrix}-1&1\\-1&1\end{pmatrix}+\frac Dh\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$ · *why:* (10.76) is the sum of the
     two integrals. · *plain:* the element stiffness matrix; each row sums to zero.
  9. *did:* Place the boundary data · *tex:* $f^1_a=-g\,k^1_{a1},\qquad f^{n_{el}}_a=Dq\,\delta_{a2},\qquad f^e_a=0\ \text{otherwise}$ (10.77) · *why:* On element 1
     local node 1 is the Dirichlet node 0: its known value g times its column moves to the right side; the flux datum lives at the last
     node. The book's text says 'A = e or e + 1'; with (10.78) it is A ∈ {e − 1, e} (slip R2). · *plain:* known values become forces.
  10. *did:* Scatter-add into node A's row · *tex:* $\text{row }A:\ \big(k^{(A)}_{21},\ k^{(A)}_{22}+k^{(A+1)}_{11},\ k^{(A+1)}_{12}\big)$ at columns (A − 1, A, A + 1) · *why:*
      Node A is local node 2 of element A and local node 1 of element A + 1 (10.78); `np.add.at` adds both contributions (P235). · *plain:* an
      interior row collects pieces from its two elements.
  11. *did:* Compare with (10.63) · *tex:* $\big(-\tfrac u2-\tfrac Dh,\ \tfrac{2D}h,\ \tfrac u2-\tfrac Dh\big),\qquad h\big(\tfrac16,\ \tfrac23,\ \tfrac16\big)$ · *why:* Insert steps 8
      and 4 into step 10 (for M: h/6, h/3 + h/3, h/6). This is (10.63) multiplied by h — D12 recovered element by element. · *plain:* assembly
      rebuilds the same row.
- **Result.** The element blocks above and $\mathbf M=\sum_e\mathbf M^e$, $\mathbf K=\sum_e\mathbf K^e$, $\mathbf F=\sum_e\mathbf F^e$ (10.69) — *in words:* one 2 × 2 block per
  element, added at the element's two global nodes, rebuilds the global matrices.
- **Check.** (h, u, D) = (0.25, 1, 0.25): m = [[0.0833, 0.0417], [0.0417, 0.0833]], k = [[0.5, −0.5], [−1.5, 1.5]], node 2's row (−1.5, 2, −0.5) ✓; the
  2-point Gauss–Legendre integration (`FEM1.element_integrals`) gives the same to 1e-15.
- **What it means.** The loop 'compute on the parent, scatter-add' is the whole of FE assembly — the same code shape handles the curved
  P2 triangles of the cylinder (N104–N109).
- **Traps.** Writing the convective block symmetric (row a, column b: the slope belongs to b); the force on element 1 (the Dirichlet
  column, not a flux); the book's A = e, e + 1 labels.

### D15 · The exact steady solution $T=\frac{e^{Rx/L}-1}{e^R-1}$ (10.86), $R=uL/D$ (10.87), its large-R form $T=e^{-R(1-x/L)}$ (10.88) and $\frac\delta L=O\big(\frac1{\lvert R\rvert}\big)$ (10.89) — ★, 8 steps, in C09 (notebook · `cell_peclet_wiggles`)
- **Goal.** Solve the steady convection–diffusion problem exactly and see the thin layer at the wall.
- **Start.** $u\frac{dT}{dx}=D\frac{d^2T}{dx^2}$ on 0 ≤ x ≤ L (10.84), T(0) = 0, T(L) = 1 (10.85) — *in words:* heat carried to a hot wall balances heat
  conducted.
- **Plan.** (1) Exponential trial. (2) Two roots, general solution. (3) Fit the end values. (4) Take R large.
- **Tools.** Linear ODE with constant coefficients by an exponential trial (Ch. 1 P44) · Péclet number (primer P229) · `np.expm1` for the
  numbers (Ch. 3 P107).
- **Assumptions.** u, D constant, D > 0; R > 0 (flow toward x = L) in steps 7–8.
- **Steps.**
  1. *did:* Rewrite as a homogeneous ODE · *tex:* $D\,T''-u\,T'=0$ · *why:* Move uT′ to the right-hand side's side; a constant-coefficient linear
     ODE. · *plain:* diffusion minus convection is zero.
  2. *did:* Try $T=e^{mx}$ · *tex:* $Dm^2-um=0$ · *why:* P44: the exponential's derivatives are multiples of itself; divide by $e^{mx}\ne0$. ·
     *plain:* an algebraic equation for m.
  3. *did:* Solve for m · *tex:* $m=0\quad\text{or}\quad m=\frac uD$ · *why:* Factor m(Dm − u) = 0. · *plain:* a constant and an exponential.
  4. *did:* Write the general solution · *tex:* $T=A+B\,e^{ux/D}$ · *why:* Two distinct roots, superposition of the two solutions. · *plain:* two
     constants to fix.
  5. *did:* Use T(0) = 0 · *tex:* $A=-B$ · *why:* At x = 0, A + B = 0. · *plain:* the inflow value.
  6. *did:* Use T(L) = 1 · *tex:* $T=\frac{e^{Rx/L}-1}{e^R-1},\qquad R=\frac{uL}D$ (10.86)–(10.87) · *why:* $B(e^{uL/D}-1)=1$, so $B=1/(e^R-1)$; write
     $ux/D=Rx/L$. · *plain:* the exact profile; one number R controls it.
  7. *did:* Let R be large · *tex:* $T\approx\frac{e^{Rx/L}}{e^R}=e^{-R(1-x/L)}$ (10.88) · *why:* For R ≫ 1, $e^R\gg1$; and wherever T is not tiny, $e^{Rx/L}\gg1$ too,
     so the 1s are negligible. · *plain:* an exponential rise at the wall.
  8. *did:* Read the layer thickness · *tex:* $1-\frac xL=\frac1R:\ T=e^{-1}=0.368;\quad\frac2R:\ T=e^{-2}=0.135\ \Rightarrow\ \frac\delta L=O\Big(\frac1{\lvert R\rvert}\Big)$ (10.89) ·
     *why:* T falls by a factor e every L/R away from the wall; e⁻¹ and e⁻² are computed, not quoted. · *plain:* the layer
     is L/R thick.
- **Result.** (10.86)–(10.89) — *in words:* at high Péclet number all the change happens in a layer of thickness L/R at the downstream wall.
- **Check.** R → 0: T → x/L (pure conduction, via e^{ε} ≈ 1 + ε) ✓; T(0) = 0, T(L) = 1 ✓; R = 4, x = 0.75: 0.3561 ✓ (the tiny example of C07). The
  printed form overflows at R ≈ 710 (e^R > 1.8 × 10³⁰⁸) — the code uses $e^{-R(1-x/L)}\frac{1-e^{-Rx/L}}{1-e^{-R}}$.
- **What it means.** A grid must put points inside this layer — C09's cell Péclet condition (D16). The same boundary-layer structure is
  Ch. 9's, with R in the role of √Re's square.
- **Traps.** Coding (10.86) as printed for large R (overflow → NaN); calling the layer thickness L·R.

### D16 · The centred discrete solution $T_j=\frac{r^j-1}{r^n-1}$, $r=\frac{1+R_{cell}/2}{1-R_{cell}/2}$: wiggles iff $R_{cell}>2$, i.e. $\delta=O\big(\frac\Delta x{R_{cell}}\big)$ (10.92) below half a cell — ★★, 12 steps, in C09 (notebook · `cell_peclet_wiggles`)
- **Goal.** Solve the centred difference equations exactly (the book only argues with a heat balance) and find precisely when they
  produce sign-alternating values.
- **Start.** Centred differences of (10.84): $u\frac{T_{j+1}-T_{j-1}}{2\Delta x}=D\frac{T_{j+1}-2T_j+T_{j-1}}{\Delta x^2}$, $T_0=0$, $T_n=1$ — *in words:* the
  discrete balance at every interior node.
- **Plan.** (1) Scale to the book's (10.90)–(10.91). (2) Try $T_j=r^j$. (3) Solve the quadratic. (4) Fit the ends. (5) Read the sign of r;
  compare with upwind.
- **Tools.** Linear recurrence with constant coefficients (primer P236) · quadratic formula.
- **Assumptions.** Uniform grid Δx = L/n; u > 0, so R_cell > 0.
- **Steps.**
  1. *did:* Multiply by Δx²/D · *tex:* $\frac{u\Delta x}{2D}(T_{j+1}-T_{j-1})=T_{j+1}-2T_j+T_{j-1}$ (10.90) · *why:* Clears both denominators. · *plain:* the
     balance in grid units.
  2. *did:* Name the cell Péclet number · *tex:* $0.5R_{cell}(T_{j+1}-T_{j-1})=T_{j+1}-2T_j+T_{j-1},\quad R_{cell}=\frac{u\Delta x}D=\frac Rn$ (10.91) · *why:* $\frac{u\Delta
     x}{2D}=\frac12R_{cell}$ and Δx = L/n. We write P = R_cell below. · *plain:* one number per cell.
  3. *did:* Try a geometric sequence · *tex:* $T_j=r^j$ · *why:* P236: constant coefficients linking neighbours ⇒ powers solve it, as $e^{mx}$
     solved the ODE (D15). · *plain:* each node is r times the previous.
  4. *did:* Substitute · *tex:* $0.5P\,(r^{j+1}-r^{j-1})=r^{j+1}-2r^j+r^{j-1}$ · *why:* Step 3 in step 2. · *plain:* the recurrence for the trial.
  5. *did:* Divide by $r^{j-1}$ · *tex:* $0.5P\,(r^2-1)=r^2-2r+1$ · *why:* $r\ne0$ (a zero root would give T ≡ 0). · *plain:* an equation for r alone.
  6. *did:* Collect as a quadratic · *tex:* $(1-0.5P)\,r^2-2r+(1+0.5P)=0$ · *why:* Move everything to one side and group the powers of r. · *plain:*
     a quadratic in r.
  7. *did:* Solve it · *tex:* $r=\frac{1\pm0.5P}{1-0.5P}:\qquad r_1=1,\quad r_2=\frac{1+P/2}{1-P/2}$ · *why:* Quadratic formula: discriminant 4 − 4(1 − P²/4) = P², so
     $r=\frac{2\pm P}{2(1-P/2)}$. · *plain:* one root is the constant, the other decides everything.
  8. *did:* Fit the two end values · *tex:* $T_j=A+B\,r_2^j,\ T_0=0,\ T_n=1\ \Rightarrow\ T_j=\frac{r^j-1}{r^n-1}\ (r=r_2)$ · *why:* Superposition of the two
     solutions; A = −B from j = 0, B = 1/(rⁿ − 1) from j = n — the discrete twin of D15 step 6. · *plain:* the exact discrete solution.
  9. *did:* Find when $r_2<0$ · *tex:* $r_2<0\iff1-\tfrac P2<0\iff P=R_{cell}>2$ · *why:* The numerator 1 + P/2 is positive for P > 0, so the sign of $r_2$ is the
     sign of the denominator. · *plain:* the root turns negative exactly at R_cell = 2.
  10. *did:* Read the consequence · *tex:* $r<0\ \Rightarrow\ r^j=\lvert r\rvert^j(-1)^j$ · *why:* Powers of a negative number alternate in sign, so $T_j$ oscillates
      node to node — the wiggles. At P = 2 exactly the r² coefficient vanishes (a degenerate case: $T_j=0$ for j < n). · *plain:* negative r means
      a zigzag.
  11. *did:* Compare with the layer thickness · *tex:* $\delta=O\Big(\frac L{nR_{cell}}\Big)=O\Big(\frac{\Delta x}{R_{cell}}\Big)$ (10.92) · *why:* (10.89) with R = nR_cell
      and L = nΔx; R_cell > 2 ⇔ δ < Δx/2: the layer is thinner than half a cell. · *plain:* wiggles mean the grid cannot see the layer.
  12. *did:* Same trial for upwind (10.93) · *tex:* $P\,(r-1)=(r-1)^2\ \Rightarrow\ r=1\ \text{or}\ r=1+P>0$ · *why:* Insert $r^j$ in $R_{cell}(T_j-T_{j-1})=T_{j+1}-2T_j+T_{j-1}$
      and divide by $r^{j-1}$; the right side is a perfect square. · *plain:* upwind's root is always positive — it never wiggles.
- **Result.** Centred: $T_j=\frac{r^j-1}{r^n-1}$ with $r=\frac{1+R_{cell}/2}{1-R_{cell}/2}$, wiggles iff R_cell > 2; upwind: r = 1 + R_cell, never — *in words:* the
  sign of one root decides whether the scheme oscillates.
- **Check.** n = 4, R = 16 (P = 4): r = −3, T = (−0.05, 0.1, −0.35) at j = 1, 2, 3 ✓ (tiny example; `FD.steady_cd_fd` agrees to 1e-13). P → 0: $r_2\approx1+P$,
  $T_j\approx j/n$ (conduction) ✓.
- **Optional sympy check (`check_src`):** `solve((1 - P/2)*r**2 - 2*r + (1 + P/2), r)` → [1, (P + 2)/(2 − P)]; substitute $T_j$ back into (10.91) and
  `simplify` to 0.
- **What it means.** R_cell ≤ 2 (10.31) is not a stability condition (the steady problem has no time step) but a resolution condition;
  the oscillation is a warning, not a bug.
- **Traps.** Missing the degenerate case P = 2; concluding wiggles need 'instability'; forgetting that $r_1=1$ is always a root (the
  constant solution).

### D17 · Upwind's modified equation $u\frac{\partial T}{\partial x}=D(1+0.5R_{cell})\frac{\partial^2T}{\partial x^2}$ (10.94) — ★★, 8 steps, in C09 (notebook · `cell_peclet_wiggles`)
- **Goal.** Show what the book says 'can be easily shown': the upwind scheme is consistent not with (10.84) but with a more diffusive
  equation.
- **Start.** $R_{cell}(T_j-T_{j-1})=T_{j+1}-2T_j+T_{j-1}$ (10.93) — *in words:* upwind convection balances centred diffusion.
- **Plan.** (1) Taylor-expand both sides. (2) Divide by Δx² keeping R_cell fixed. (3) Read the extra diffusion.
- **Tools.** Taylor series (R02) · D01 results.
- **Assumptions.** T smooth; R_cell held fixed while Δx → 0 (step 7 — the book's parenthesis).
- **Steps.**
  1. *did:* Note which difference this is · *tex:* $T_j-T_{j-1}\ \text{is a backward (upwind) difference}$ · *why:* It uses the node behind j, where the
     flow (u > 0) comes from. The book calls it 'forward-difference' (slip R3). · *plain:* information taken from upstream.
  2. *did:* Expand the convective difference · *tex:* $T_j-T_{j-1}=\Delta x\,T_x-\frac{\Delta x^2}2T_{xx}+\frac{\Delta x^3}6T_{xxx}+O(\Delta x^4)$ · *why:* (10.5)
     rearranged (R02), all derivatives at $x_j$. · *plain:* slope, minus half a curvature, …
  3. *did:* Expand the diffusive difference · *tex:* $T_{j+1}-2T_j+T_{j-1}=\Delta x^2T_{xx}+\frac{\Delta x^4}{12}T_{xxxx}+O(\Delta x^6)$ · *why:* D01 step 7. · *plain:* the
     curvature term.
  4. *did:* Substitute into (10.93) · *tex:* $R_{cell}\Big(\Delta x\,T_x-\frac{\Delta x^2}2T_{xx}+\dots\Big)=\Delta x^2T_{xx}+O(\Delta x^4)$ · *why:* Steps 2 and 3 in (10.93). ·
     *plain:* both sides as series.
  5. *did:* Divide by Δx² · *tex:* $\frac{R_{cell}}{\Delta x}T_x-\frac{R_{cell}}2T_{xx}+O(R_{cell}\Delta x)=T_{xx}+O(\Delta x^2)$ · *why:* To compare with a differential
     equation per unit length². · *plain:* the balance per unit Δx².
  6. *did:* Use $R_{cell}/\Delta x=u/D$ · *tex:* $\frac uDT_x=\Big(1+\frac{R_{cell}}2\Big)T_{xx}+O(R_{cell}\Delta x,\Delta x^2)$ · *why:* R_cell = uΔx/D by definition; move
     the $T_{xx}$ term to the right. · *plain:* the convection now balances a larger diffusion.
  7. *did:* Multiply by D, drop the O(Δx) terms · *tex:* $u\frac{\partial T}{\partial x}=D\,(1+0.5R_{cell})\frac{\partial^2T}{\partial x^2}$ (10.94) · *why:* With R_cell held
     fixed, the remaining terms vanish as Δx → 0 — the book's 'when the cell Péclet number is kept finite'. · *plain:* upwind solves a more
     diffusive problem.
  8. *did:* Name the numerical diffusivity · *tex:* $D_{num}=0.5R_{cell}D=\frac{u\Delta x}2$ · *why:* The difference between the two diffusivities in (10.94). Small
     compared with D only if 0.5 R_cell ≪ 1. · *plain:* upwinding adds a fake diffusivity u Δx/2.
- **Result.** (10.94) — *in words:* first-order upwinding is centred differencing of a fluid whose diffusivity has been raised by uΔx/2.
- **Check.** Units m²/s ✓. R_cell = 4: D_num = 2D, so the layer is three times too thick (tiny example: R_eff = 16/3) ✓. The upwind nodes approach
  the modified-equation solution as n grows at fixed R (verifier note (i)).
- **What it means.** 'Stable but wrong': upwinding removes wiggles by smearing; ocean tracer schemes pay exactly this (C09 climate
  hook). In time-dependent advection the same move gives ∣u∣Δx(1 − C)/2 (C05).
- **Traps.** Calling T_j − T_{j−1} a forward difference; letting R_cell → 0 with Δx (then upwind looks consistent with (10.84) and the
  message is lost); dropping the ½.

### D18 · MacCormack = Lax–Wendroff: second order with local error $-\frac{u\Delta t\Delta x^2}6(1-C^2)T_{xxx}$, and $\lvert G\rvert^2=1-4C^2(1-C^2)\sin^4\frac\theta2$ ⇒ stable iff C ≤ 1, from (10.100)–(10.102) — ★★★, 14 steps, in C10 (notebook · `upwind_cfl_advection`)
- **Goal.** Show what the book only states — that the predictor–corrector (10.101)–(10.102) is second order in space and time — and find
  its stability limit, on the simplest system it applies to.
- **Start.** The scheme for $\mathbf U_t+\mathbf E(\mathbf U)_x+\mathbf F(\mathbf U)_y=0$ (10.100): predictor $\mathbf U^*_{i,j}=\mathbf U^n_{i,j}-\frac{\Delta t}{\Delta x}(\mathbf E^n_{i+1,j}-
  \mathbf E^n_{i,j})-\frac{\Delta t}{\Delta y}(\mathbf F^n_{i,j+1}-\mathbf F^n_{i,j})$ (10.101), corrector $\mathbf U^{n+1}_{i,j}=\frac12\big[\mathbf U^n_{i,j}+\mathbf U^*_{i,j}-\frac{\Delta t}{\Delta x}
  (\mathbf E^*_{i,j}-\mathbf E^*_{i-1,j})-\frac{\Delta t}{\Delta y}(\mathbf F^*_{i,j}-\mathbf F^*_{i,j-1})\big]$ (10.102) — *in words:* predict with forward differences,
  correct with backward ones, average.
- **Plan.** (1) Apply it to one scalar in 1-D, T_t + uT_x = 0 (E = uT, F = 0). (2) Eliminate the predicted values: Lax–Wendroff. (3) Compare
  with the exact solution's Taylor series (using the PDE twice). (4) Find G and its modulus.
- **Tools.** Predictor–corrector (primer P238) · Taylor series in t and x (Ch. 3 P98) · Schwarz's theorem, to swap ∂t and ∂x (Ch. 4 P121) ·
  D04 moves · half-angle identities (P225) · sympy (Ch. 1 P40, Ch. 4 P117).
- **Assumptions.** Constant u > 0, $C=u\Delta t/\Delta x$ (step 1); smooth solution (steps 6–11); periodic grid (steps 13–14).
- **Steps.**
  1. *did:* Specialise the predictor · *tex:* $T^*_i=T^n_i-C\,(T^n_{i+1}-T^n_i)$ · *why:* (10.101) with $\mathbf E=uT$, no y-direction; $\frac{\Delta t}{\Delta x}u=C$. ·
     *plain:* a forward-difference guess at the new level.
  2. *did:* Specialise the corrector · *tex:* $T^{n+1}_i=\tfrac12\big[T^n_i+T^*_i-C\,(T^*_i-T^*_{i-1})\big]$ · *why:* (10.102) with $\mathbf E^*=uT^*$. · *plain:* a backward
     difference of the guess, averaged with the old value.
  3. *did:* Difference the predicted values · *tex:* $T^*_i-T^*_{i-1}=(T_i-T_{i-1})-C\,(T_{i+1}-2T_i+T_{i-1})$ · *why:* Write step 1 at i and at i − 1 and
     subtract; the C terms combine into a second difference. (Superscript n dropped.) · *plain:* the guess's backward difference.
  4. *did:* Substitute into the corrector · *tex:* $T^{n+1}_i=\tfrac12\big[2T_i-C(T_{i+1}-T_i)-C(T_i-T_{i-1})+C^2(T_{i+1}-2T_i+T_{i-1})\big]$ · *why:* Insert step 1
     for $T^*_i$ and step 3; this is where C² enters — the corrector differences the predictor. · *plain:* no starred values left.
  5. *did:* Combine the first differences · *tex:* $T^{n+1}_i=T_i-\frac C2(T_{i+1}-T_{i-1})+\frac{C^2}2(T_{i+1}-2T_i+T_{i-1})$ · *why:* $(T_{i+1}-T_i)+(T_i-T_{i-1})=T_{i+1}-T_{i-1}$;
     halve everything. This is the Lax–Wendroff scheme. · *plain:* a centred convective step plus a C²/2 diffusion-like correction.
  6. *did:* Taylor-expand the exact solution in time · *tex:* $T(t+\Delta t)=T+\Delta t\,T_t+\frac{\Delta t^2}2T_{tt}+\frac{\Delta t^3}6T_{ttt}+O(\Delta t^4)$ · *why:* P98 at
     fixed x; we will compare the scheme with this. · *plain:* what the exact solution does in one step.
  7. *did:* Trade time derivatives for space ones · *tex:* $T_t=-uT_x,\quad T_{tt}=u^2T_{xx},\quad T_{ttt}=-u^3T_{xxx}$ · *why:* Differentiate the PDE in t and use
     $\partial_t\partial_x=\partial_x\partial_t$ (P121) with $T_t=-uT_x$ again; each time derivative becomes −u ∂x. · *plain:* in pure advection, time
     change is space change carried along.
  8. *did:* Exact step in space derivatives · *tex:* $T(t+\Delta t)=T-u\Delta t\,T_x+\frac{u^2\Delta t^2}2T_{xx}-\frac{u^3\Delta t^3}6T_{xxx}+\dots$ · *why:* Step 7 in
     step 6. · *plain:* the target, term by term.
  9. *did:* Expand the scheme's differences · *tex:* $\tfrac12(T_{i+1}-T_{i-1})=\Delta x\,T_x+\frac{\Delta x^3}6T_{xxx}+\dots,\quad T_{i+1}-2T_i+T_{i-1}=\Delta x^2T_{xx}+O(\Delta x^4)$ ·
     *why:* D01 steps 5 and 7. · *plain:* the stencils as series.
  10. *did:* Scheme step in space derivatives · *tex:* $T^{n+1}_i=T-u\Delta t\,T_x+\frac{u^2\Delta t^2}2T_{xx}-\frac{u\Delta t\,\Delta x^2}6T_{xxx}+\dots$ · *why:* Step 9 in
      step 5 with CΔx = uΔt and C²Δx² = u²Δt². · *plain:* the scheme's step, term by term.
  11. *did:* Subtract exact from scheme · *tex:* $T^{n+1}_i-T(t+\Delta t)=-\frac{u\Delta t\,\Delta x^2}6T_{xxx}+\frac{u^3\Delta t^3}6T_{xxx}+\dots$ · *why:* Steps 10 − 8: the
      T, T_x and T_xx terms are identical and cancel — the whole point of the C²/2 correction. · *plain:* only a third-derivative term is left.
  12. *did:* Factor and read the order · *tex:* $T^{n+1}_i-T(t+\Delta t)=-\frac{u\Delta t\,\Delta x^2}6\,(1-C^2)\,T_{xxx}+O(\Delta t\,\Delta x^4)$ · *why:* $u^3\Delta t^3=u\Delta t\cdot
      C^2\Delta x^2$. Per unit time (÷ Δt) the error is O(Δx², Δt²): second order in both, as the book states; an odd derivative, so
      dispersive; zero at C = 1. · *plain:* second order — waves keep their height but short ones lag.
  13. *did:* Insert one Fourier wave · *tex:* $G=1-\mathrm iC\sin\theta-C^2(1-\cos\theta)=1-2C^2s-\mathrm iC\sin\theta,\quad s=\sin^2\tfrac\theta2$ · *why:* D04 moves in step
      5: $\frac12(e^{\mathrm i\theta}-e^{-\mathrm i\theta})=\mathrm i\sin\theta$, $e^{\mathrm i\theta}-2+e^{-\mathrm i\theta}=-2(1-\cos\theta)$, then the half-angle identity (P225). ·
      *plain:* one complex factor per wave.
  14. *did:* Square, add, read the limit · *tex:* $\lvert G\rvert^2=(1-2C^2s)^2+4C^2s(1-s)=1-4C^2(1-C^2)\sin^4\tfrac\theta2\le1\iff C\le1$ · *why:* $\sin^2\theta=4s(1-s)$;
      expanding, the 4C²s terms cancel, leaving $1-4C^2s^2+4C^4s^2$. It is ≤ 1 for every θ exactly when C² ≤ 1: a CFL condition. ·
      *plain:* stable exactly when the Courant number is at most 1.
- **Result.** For linear advection MacCormack is Lax–Wendroff, $T^{n+1}_i=T_i-\frac C2(T_{i+1}-T_{i-1})+\frac{C^2}2(T_{i+1}-2T_i+T_{i-1})$, second order with a
  dispersive leading error $-\frac{u\Delta x^2}6(1-C^2)T_{xxx}$ per unit time, and stable iff C ≤ 1 — *in words:* forward-then-backward averages
  to a centred, second-order step whose price is phase error, not amplitude loss.
- **Check.** Units: the local error is in T ✓ (u Δt Δx² T_xxx: m/s · s · m² · T/m³). Tiny example (C = 0.5, [0, 0, 1, 0, 0]) gives [0, −0.125,
  0.75, 0.375, 0] by both routes ✓. C = 1: G = e^{−iθ}, the exact shift, and the local error vanishes ✓. C = 0: G = 1 ✓. Measured order ≈ 2 on
  the Gaussian (C10 code cell).
- **sympy check intent (`check_src`, ★★★, every line commented):** re-run the construction, not just the answer —
  ```python
  import sympy as sp                                      # symbolic algebra (Ch. 1 P40)
  C, th, a = sp.symbols('C theta a', real=True)           # Courant number, Fourier angle, half-angle a = theta/2
  Tm, T0, Tp = sp.symbols('T_m T_0 T_p')                  # T at i-1, i, i+1 (a stencil of three values)
  pred = lambda left, right: left - C*(right - left)      # predictor (10.101) with E = uT: T*_i from T_i and T_{i+1}
  Ts0, Tsm = pred(T0, Tp), pred(Tm, T0)                   # T*_i and T*_{i-1}
  corr = sp.Rational(1, 2)*(T0 + Ts0 - C*(Ts0 - Tsm))     # corrector (10.102)
  lw = T0 - C/2*(Tp - Tm) + C**2/2*(Tp - 2*T0 + Tm)       # Lax-Wendroff, step 5
  print(sp.simplify(corr - lw))                           # 0: steps 1-5 are right
  x, dx = sp.symbols('x dx', positive=True)               # position and grid step (C is held fixed)
  f = sp.Function('f')                                    # any smooth profile; the exact solution is T = f(x - u t)
  scheme = lw.subs({Tm: f(x - dx), T0: f(x), Tp: f(x + dx)})   # the scheme applied to exact values at t = 0
  exact = f(x - C*dx)                                     # exact value one step later: shifted by u dt = C dx
  err = sp.series(scheme - exact, dx, 0, 4).removeO().doit()   # Taylor in dx (steps 6-11 done by machine)
  print(sp.factor(sp.simplify(err)))                      # C (C-1)(C+1) dx^3 f_xxx/6 = -(u dt dx^2/6)(1 - C^2) T_xxx: step 12
  G = 1 - sp.I*C*sp.sin(th) - C**2*(1 - sp.cos(th))       # step 13
  G2 = sp.expand(sp.re(G)**2 + sp.im(G)**2)               # |G|^2 with real C and theta
  d = G2 - (1 - 4*C**2*(1 - C**2)*sp.sin(th/2)**4)        # difference from the claimed form of step 14
  print(sp.simplify(sp.expand_trig(d.subs(th, 2*a))))     # 0 once theta = 2a is expanded by double-angle rules
  ```
  (`ch10.maccormack_linear_sympy()` packages the same checks for the tests.)
- **What it means.** MacCormack's second order is what lets the cavity and block computations of §10.5 run on modest grids; its
  dispersive error is the ripple behind the square pulse (C10 figure, E3). With the Navier–Stokes fluxes the same structure holds, with
  the sound speed setting C (10.110), (10.155). It fails at discontinuities (Gibbs-like overshoots — Ch. 15's shocks need limiters).
- **Traps.** Using $\mathbf E^n$ instead of $\mathbf E^*$ in the corrector (then C² never appears and the scheme is first order — FTCS-like and
  unstable); forgetting $T_{tt}=u^2T_{xx}$ (a spurious first-order term remains); reading the local error as the order (divide by Δt).

### D19 · The continuous projection: $\nabla^2p^{n+1}=\frac1{\Delta t}\nabla\cdot\mathbf u^{n+1/2}$ from (10.117)–(10.118), and the correction is curl-free — ★★, 9 steps, in C11 (notebook · `mac_projection_staggered`)
- **Goal.** Find the equation that fixes the pressure in the projection step and show that the correction removes divergence without
  touching the vorticity. (The book writes only the discrete (10.124); this continuous version is ours.)
- **Start.** $\frac{\mathbf u^{n+1}-\mathbf u^{n+1/2}}{\Delta t}+\nabla p^{n+1}=\mathbf 0$ (10.117) and $\nabla\cdot\mathbf u^{n+1}=0$ (10.118) — *in words:* add a pressure push
  so that the new velocity is divergence-free.
- **Plan.** (1) Take the divergence of (10.117). (2) Impose (10.118) → a Poisson equation. (3) Boundary and solvability conditions.
  (4) Take the curl → vorticity untouched.
- **Tools.** Divergence of a gradient = Laplacian (Ch. 2 §2.9) · curl of a gradient = 0 (Ch. 2 §2.9) · Gauss' theorem (Ch. 2 §2.12) ·
  Helmholtz–Hodge decomposition (primer P240) · Poisson equation (Ch. 5 P139).
- **Assumptions.** Δt > 0 constant; walls with prescribed normal velocity (step 6).
- **Steps.**
  1. *did:* Solve (10.117) for the new velocity · *tex:* $\mathbf u^{n+1}=\mathbf u^{n+1/2}-\Delta t\,\nabla p^{n+1}$ · *why:* Multiply by Δt and rearrange; the new
     velocity is the predicted one minus a gradient. · *plain:* correct the prediction with a pressure push.
  2. *did:* Take the divergence of both sides · *tex:* $\nabla\cdot\mathbf u^{n+1}=\nabla\cdot\mathbf u^{n+1/2}-\Delta t\,\nabla\cdot\nabla p^{n+1}$ · *why:* The divergence is
     linear; Δt is a constant. · *plain:* how much the push changes the net outflow.
  3. *did:* Recognise the Laplacian · *tex:* $\nabla\cdot\nabla p=\nabla^2p$ · *why:* Ch. 2: ∂/∂x_i(∂p/∂x_i) is the Laplacian. · *plain:* the push changes the
     outflow by Δt × the curvature of p.
  4. *did:* Impose incompressibility (10.118) · *tex:* $0=\nabla\cdot\mathbf u^{n+1/2}-\Delta t\,\nabla^2p^{n+1}$ · *why:* ∇·u^{n+1} = 0 is *imposed* — it is what we want
     the new field to satisfy, not something derived. · *plain:* the push must cancel the predicted divergence.
  5. *did:* Rearrange: the pressure Poisson equation · *tex:* $\nabla^2p^{n+1}=\frac1{\Delta t}\nabla\cdot\mathbf u^{n+1/2}$ · *why:* Divide by Δt. A Poisson equation
     (P139): p is found by one elliptic solve, everywhere at once — the 'instantaneous' pressure of N51. · *plain:* where fluid piles up,
     the pressure has a maximum.
  6. *did:* Dot (10.117) with the wall normal · *tex:* $\frac{\partial p^{n+1}}{\partial n}=\frac1{\Delta t}\,\mathbf n\cdot(\mathbf u^{n+1/2}-\mathbf u^{n+1})$ · *why:* At a wall
     $\mathbf n\cdot\mathbf u^{n+1}$ is prescribed; if the predictor already carries the wall value, $\partial p/\partial n=0$ — the Neumann condition that the
     staggered grid builds in (D20). · *plain:* the boundary condition for p comes from the normal velocity.
  7. *did:* Check solvability with Gauss · *tex:* $\int_V\nabla^2p\,dV=\oint\frac{\partial p}{\partial n}dA\ \Rightarrow\ \frac1{\Delta t}\int_V\nabla\cdot\mathbf u^{n+1/2}dV=\oint\frac{\partial p}
     {\partial n}dA$ · *why:* Gauss' theorem applied to ∇p; with ∂p/∂n = 0 the predicted field must have zero net outflow — true when the walls'
     normal velocities balance. · *plain:* a closed box can only redistribute fluid, not create it.
  8. *did:* Take the curl of step 1 · *tex:* $\nabla\times\mathbf u^{n+1}=\nabla\times\mathbf u^{n+1/2}-\Delta t\,\nabla\times\nabla p^{n+1}=\nabla\times\mathbf u^{n+1/2}$ · *why:* The curl
     of a gradient is zero (Ch. 2). · *plain:* the projection does not change the vorticity at all.
  9. *did:* Read it as a decomposition · *tex:* $\mathbf u^{n+1/2}=\underbrace{\mathbf u^{n+1}}_{\nabla\cdot=0}+\underbrace{\nabla(\Delta t\,p^{n+1})}_{\nabla\times=0}$ · *why:*
     Step 1 rearranged: the Helmholtz–Hodge split (P240); the projection keeps the divergence-free part and deletes the gradient part —
     the book's 'irrotational correction field … proportional to the pressure'. · *plain:* keep the swirl, remove the squeeze.
- **Result.** $\nabla^2p^{n+1}=\frac1{\Delta t}\nabla\cdot\mathbf u^{n+1/2}$ with $\partial p/\partial n$ from the normal velocity, and $\nabla\times(\mathbf u^{n+1}-\mathbf u^{n+1/2})=0$ —
  *in words:* the pressure is whatever removes the divergence; its push is a pure gradient.
- **Check.** Units: with p the kinematic (or dimensionless) pressure, ∇²p is (m²/s²)/m² = 1/s² and ∇·u/Δt is (1/s)/s = 1/s² ✓. Tiny example:
  u* = (x + y, 0), Δt = 1 ⇒ p = x²/2, u = (y, 0) ✓ (divergence 0; vorticity −1 before and after). `ch10.projection_sympy()`.
- **What it means.** Every incompressible (and Boussinesq/anelastic) code solves this Poisson problem each step — the pressure solve of
  ocean models. It fails to be exact in time: splitting makes the scheme first order in Δt (N61).
- **Traps.** Thinking ∇·u^{n+1} = 0 is derived (it is imposed); forgetting that p needs a boundary condition in the continuous
  problem; taking the divergence of (10.116) instead of (10.117).

### D20 · The discrete pressure Poisson equation (10.124) from (10.121)–(10.123), no pressure boundary condition, and solvability Σ rhs = 0 — ★★, 12 steps, in C12 (notebook · `mac_projection_staggered`)
- **Goal.** Derive (10.124) cell by cell on the staggered grid, including the boundary cells the book discusses in words, and find when it
  can be solved.
- **Start.** The corrections $u^{n+1}_{i+1/2,j}=u^{n+1/2}_{i+1/2,j}-\frac{\Delta t}{\Delta x}(p^{n+1}_{i+1,j}-p^{n+1}_{i,j})$ (10.121), $v^{n+1}_{i,j+1/2}=v^{n+1/2}_{i,j+1/2}-\frac{\Delta t}
  {\Delta y}(p^{n+1}_{i,j+1}-p^{n+1}_{i,j})$ (10.122) and the cell continuity $\frac{u^{n+1}_{i+1/2,j}-u^{n+1}_{i-1/2,j}}{\Delta x}+\frac{v^{n+1}_{i,j+1/2}-v^{n+1}_{i,j-1/2}}{\Delta y}=0$
  (10.123) — *in words:* face pushes and a cell budget.
- **Plan.** (1) Insert the face corrections into the budget of an interior cell. (2) Repeat for a wall cell. (3) Look at the matrix: rows sum
  to zero ⇒ solvability and a free constant.
- **Tools.** Half-index notation (primer P242) · singular systems and compatibility (primer P243) · sparse matrices (Ch. 6 P161, primer
  P244).
- **Assumptions.** Uniform staggered grid; the normal velocity on every boundary face prescribed and not corrected (step 7).
- **Steps.**
  1. *did:* Start from one cell's budget · *tex:* $\frac{u^{n+1}_{i+1/2,j}-u^{n+1}_{i-1/2,j}}{\Delta x}+\frac{v^{n+1}_{i,j+1/2}-v^{n+1}_{i,j-1/2}}{\Delta y}=0$ (10.123) · *why:* The
     new velocity must be divergence-free in every cell. · *plain:* what flows in flows out.
  2. *did:* Insert (10.121) on the two x-faces · *tex:* $u^{n+1}_{i\pm1/2,j}=u^{n+1/2}_{i\pm1/2,j}-\frac{\Delta t}{\Delta x}(p_{i\pm1/2+1/2,j}-p_{i\pm1/2-1/2,j})$ · *why:* Each interior
     face is corrected by the difference of the two pressures beside it (superscript n + 1 on p dropped). · *plain:* both side faces
     pushed.
  3. *did:* Collect the x-part · *tex:* $\frac{u^{n+1/2}_{i+1/2,j}-u^{n+1/2}_{i-1/2,j}}{\Delta x}-\frac{\Delta t}{\Delta x^2}\big[(p_{i+1,j}-p_{i,j})-(p_{i,j}-p_{i-1,j})\big]$ · *why:*
     Subtract the two face expressions and divide by Δx. · *plain:* predicted x-outflow minus the pressure's second difference.
  4. *did:* The same in y · *tex:* $\frac{v^{n+1/2}_{i,j+1/2}-v^{n+1/2}_{i,j-1/2}}{\Delta y}-\frac{\Delta t}{\Delta y^2}\big[p_{i,j+1}-2p_{i,j}+p_{i,j-1}\big]$ · *why:* (10.122) on the two
     y-faces. · *plain:* the other direction.
  5. *did:* Add and set to zero · *tex:* $\Delta t\Big[\frac{p_{i+1,j}-2p_{i,j}+p_{i-1,j}}{\Delta x^2}+\frac{p_{i,j+1}-2p_{i,j}+p_{i,j-1}}{\Delta y^2}\Big]=\nabla_d\cdot\mathbf u^{n+1/2}_{i,j}$ · *why:*
     Steps 3 + 4 = 0 by step 1; move the pressure terms to one side. · *plain:* Δt × Laplacian of p = predicted divergence.
  6. *did:* Divide by Δt · *tex:* $\nabla^2_dp^{n+1}_{i,j}=\frac1{\Delta t}\Big(\frac{u^{n+1/2}_{i+1/2,j}-u^{n+1/2}_{i-1/2,j}}{\Delta x}+\frac{v^{n+1/2}_{i,j+1/2}-v^{n+1/2}_{i,j-1/2}}{\Delta y}\Big)$ (10.124) ·
     *why:* The 5-point Laplacian of Ch. 6 on p, equal to the divergence over Δt — the discrete twin of D19 step 5. · *plain:* the book's
     Poisson equation.
  7. *did:* Look at a cell next to the left wall · *tex:* $u^{n+1}_{1/2,j}=u_{wall}\quad(\text{not corrected})$ · *why:* The wall face's normal velocity is
     prescribed, so (10.121) is not applied there (Fig. 10.4's $u_{1/2,2}$). · *plain:* the wall face stays as it is.
  8. *did:* Redo step 3 for that cell · *tex:* $\frac{u^{n+1/2}_{3/2,j}-u_{wall}}{\Delta x}-\frac{\Delta t}{\Delta x^2}(p_{2,j}-p_{1,j})$ · *why:* Only the interior face carries a
     pressure difference; a pressure $p_{0,j}$ outside the wall never enters. · *plain:* one pressure difference instead of two.
  9. *did:* Read the boundary row · *tex:* $\frac{p_{2,j}-p_{1,j}}{\Delta x^2}+(\text{y-part})=\text{rhs}_{1,j}$ · *why:* No pressure boundary condition was needed:
     the row simply lacks its outside neighbour — equivalent to $\partial p/\partial n=0$ with a ghost $p_{0,j}=p_{1,j}$ (D19 step 6). The book
     says '(10.120)' here; it means (10.124) (slip R4). · *plain:* the wall is built into the matrix.
  10. *did:* Notice every row sums to zero · *tex:* $\mathbf A\,\mathbf 1=\mathbf 0$ · *why:* Interior rows have −2, 1, 1 per direction; boundary rows have one
      fewer neighbour *and* one fewer −1 on the diagonal. So a constant p is invisible. · *plain:* adding a constant to p changes nothing.
  11. *did:* Sum all equations: compatibility · *tex:* $\sum_{i,j}\text{rhs}_{i,j}=\frac1{\Delta t}\sum_{\text{cells}}\nabla_d\cdot\mathbf u^{n+1/2}=\frac1{\Delta t\,\Delta x\Delta y}(\text{net
      boundary flux})=0$ · *why:* A singular system is solvable only if the right side is orthogonal to the null vector (P243); the divergence
      telescopes to the wall fluxes, whose sum is zero in a closed box. · *plain:* the box cannot gain or lose fluid.
  12. *did:* Fix the free constant · *tex:* $p_{1,1}=0\ \ \text{or}\ \ \textstyle\sum p_{i,j}=0$ · *why:* One extra condition removes the null space; velocities use
      only differences of p, so they do not depend on the choice (N41). · *plain:* pin one pressure.
- **Result.** (10.124) with boundary rows lacking their outside neighbour, a constant null vector, solvable iff Σ rhs = 0, one value pinned —
  *in words:* the staggered grid turns incompressibility into a Neumann Poisson problem with no pressure boundary condition to invent.
- **Check.** 3-cell pipe (C12 tiny example): rows (−1, 1, 0), (1, −2, 1), (0, 1, −1) sum to zero; rhs (1, 0, −1) sums to 0; p = (0, 1, 2); corrected
  velocities 0 ✓. The 8 × 8 matrix has rank 63 before pinning ✓ (C12 code); perturbing the boundary u^{n+1/2} leaves the result unchanged (the
  Peyret–Taylor remark, a test).
- **What it means.** The MAC grid's great convenience: the pressure equation writes itself, walls included. The same Neumann Poisson
  problem is solved by every Boussinesq ocean model each step.
- **Traps.** Correcting the wall face (then $p_{0,j}$ appears and seems to need a condition); forgetting the compatibility condition (the
  solver fails or drifts); pinning a pressure and then comparing pressures between codes without removing the constant.

### D21 · The checkerboard $p_{i,j}=(-1)^{i+j}$: zero collocated gradient (10.125), ±2/Δx staggered gradient (10.126) — ★★, 7 steps, in C12 (notebook · `mac_projection_staggered`)
- **Goal.** Show why a zigzag pressure is invisible on a collocated grid and visible on the staggered one — the book says it in words.
- **Start.** The two gradients $\big(\frac{\partial p}{\partial x}\big)_{i,j}=\frac{p_{i+1,j}-p_{i-1,j}}{2\Delta x}$ (10.125) and $\big(\frac{\partial p}{\partial x}\big)_{i+1/2,j}=
  \frac{p_{i+1,j}-p_{i,j}}{\Delta x}$ (10.126) — *in words:* skip the neighbour, or use it.
- **Plan.** (1) Define the zigzag. (2) Evaluate both gradients. (3) Conclude about null spaces; repeat for velocity.
- **Tools.** Stencil arithmetic (D01) · null space (Ch. 1 P58; SVD primer P245).
- **Assumptions.** Uniform grid (both directions alike).
- **Steps.**
  1. *did:* Define the checkerboard · *tex:* $p_{i,j}=(-1)^{i+j}$ · *why:* The shortest pattern a grid can hold: +1 and −1 alternating in both directions,
     like a chessboard. · *plain:* a pressure that zigzags.
  2. *did:* Look at the collocated neighbours · *tex:* $p_{i\pm1,j}=(-1)^{i\pm1+j}=-(-1)^{i+j}$ · *why:* Moving one cell flips the sign; moving one cell
     left or right flips it the same way. · *plain:* both neighbours have the same value.
  3. *did:* Evaluate (10.125) · *tex:* $\frac{p_{i+1,j}-p_{i-1,j}}{2\Delta x}=\frac{-(-1)^{i+j}+(-1)^{i+j}}{2\Delta x}=0$ · *why:* Step 2 in (10.125); the same in y. ·
     *plain:* the collocated gradient of the zigzag is zero everywhere.
  4. *did:* Conclude: an invisible mode · *tex:* $\nabla_{coll}\,p_{checker}=\mathbf 0\ \Rightarrow\ p_{checker}\in\ker\nabla_{coll}$ · *why:* The momentum equation
     sees p only through its gradient, so this p behaves like a constant: nothing resists its growth, and the Poisson problem has extra
     null vectors (the SVD count ≥ 4 on an even periodic grid, P245). · *plain:* the zigzag hides.
  5. *did:* Evaluate (10.126) · *tex:* $\frac{p_{i+1,j}-p_{i,j}}{\Delta x}=\frac{-(-1)^{i+j}-(-1)^{i+j}}{\Delta x}=-\frac{2(-1)^{i+j}}{\Delta x}$ · *why:* Adjacent values differ in
     sign, so they do not cancel. · *plain:* the staggered gradient is ±2/Δx.
  6. *did:* Conclude: only the constant is invisible · *tex:* $\lvert\nabla_{stag}\,p_{checker}\rvert=\frac2{\Delta x}\ne0,\qquad\ker\nabla_{stag}=\{\text{const}\}$ · *why:* A
     zero face difference at every face forces all neighbours equal, i.e. a constant. · *plain:* on the staggered grid the zigzag is
     pushed back.
  7. *did:* The same for velocity · *tex:* $\frac{u_{i+1}-u_{i-1}}{2\Delta x}=0\ \text{for}\ u_i=(-1)^i,\qquad\frac{u_{i+1/2}-u_{i-1/2}}{\Delta x}=\pm\frac2{\Delta x}$ · *why:* A
     collocated divergence skips the neighbour too, so a zigzag velocity passes continuity; the staggered (10.123) uses adjacent faces. ·
     *plain:* staggering protects both equations.
- **Result.** Collocated centred gradient of the checkerboard = 0; staggered = ±2/Δx — *in words:* compact two-point differences see every
  pattern except a constant.
- **Check.** `ch10.collocated_gradient(ch10.checkerboard(8, 8), 1, 1)` is 0 in the interior; `MAC.gradient` gives ±2 ✓; null-space dimensions ≥ 4
  and 1 (`ch10.gradient_null_space`).
- **What it means.** The FD root of the FE problem of C13 (spurious pressure modes); it is why ocean and atmosphere models use the C-grid
  and why collocated codes need special interpolation (Rhie–Chow, named only).
- **Traps.** Testing a smooth pressure (both gradients agree); thinking the zigzag is an instability of the time stepping (it is a
  null space of the operator).

### D22 · The weak Navier–Stokes equations (10.134)–(10.135) — ★★, 11 steps, in C13 (notebook · `mixed_fe_lbb`)
- **Goal.** Apply C06's recipe to the Navier–Stokes equations and obtain the form in which velocity and pressure appear in mixed finite
  elements.
- **Start.** $\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=\mathbf g-\nabla p+\frac1{Re}\nabla^2\mathbf u$ (10.81) and $\nabla\cdot\mathbf u=0$ (10.80) — *in words:*
  dimensionless momentum and mass.
- **Plan.** (1) Dot with a velocity test function and integrate. (2) Move one derivative off the pressure and off the viscous term
  (divergence theorem). (3) Use ∇·u = 0 and the symmetry of D. (4) Test continuity with a pressure test function.
- **Tools.** Product rule for a divergence (Ch. 4 P113) · the divergence theorem (Ch. 2 §2.12) · the double dot product (Ch. 2) · the
  strain-rate tensor (R10) · D09's moves.
- **Assumptions.** ũ = 0 on the boundary where u is prescribed (steps 3, 7); ∇·u = 0 (step 5); Newtonian, constant Re.
- **Steps.**
  1. *did:* Dot with ũ and integrate · *tex:* $\int_\Omega\Big(\frac{\partial\mathbf u}{\partial t}+\mathbf u\cdot\nabla\mathbf u-\mathbf g\Big)\cdot\tilde{\mathbf u}\,d\Omega=\int_\Omega
     \Big(-\nabla p+\frac1{Re}\nabla^2\mathbf u\Big)\cdot\tilde{\mathbf u}\,d\Omega$ · *why:* As in D09 steps 1–2, with a vector test function (the velocity
     variation ũ). · *plain:* the momentum equation, weighted and averaged.
  2. *did:* Product rule on the pressure term · *tex:* $\nabla p\cdot\tilde{\mathbf u}=\nabla\cdot(p\tilde{\mathbf u})-p\,\nabla\cdot\tilde{\mathbf u}$ · *why:* P113: ∇·(pũ) = ∇p·ũ +
     p∇·ũ. · *plain:* trade the pressure gradient for a divergence.
  3. *did:* Integrate; the surface term vanishes · *tex:* $-\int_\Omega\nabla p\cdot\tilde{\mathbf u}\,d\Omega=-\oint p\,\tilde{\mathbf u}\cdot\mathbf n\,dA+\int_\Omega p\,\nabla\cdot
     \tilde{\mathbf u}\,d\Omega=\int_\Omega p\,\nabla\cdot\tilde{\mathbf u}\,d\Omega$ · *why:* Divergence theorem; ũ = 0 where velocity is prescribed (elsewhere the term is
     part of a natural traction condition). · *plain:* the pressure now multiplies the divergence of the test function.
  4. *did:* Rewrite the Laplacian · *tex:* $\nabla\cdot\big(\nabla\mathbf u+(\nabla\mathbf u)^T\big)=\nabla^2\mathbf u+\nabla(\nabla\cdot\mathbf u)$ · *why:* In components
     $\partial_j(\partial_ju_i+\partial_iu_j)=\partial_j\partial_ju_i+\partial_i\partial_ju_j$, swapping derivatives (Schwarz). · *plain:* the Laplacian plus a
     gradient of the divergence.
  5. *did:* Use ∇·u = 0 · *tex:* $\frac1{Re}\nabla^2\mathbf u=\frac2{Re}\nabla\cdot\mathbf D[\mathbf u],\qquad\mathbf D[\mathbf u]=\tfrac12\big[\nabla\mathbf u+(\nabla\mathbf u)^T\big]$ (10.136) ·
     *why:* The second term of step 4 vanishes for incompressible flow; the bracket is 2D (R10). · *plain:* viscous force = divergence of twice
     the strain rate.
  6. *did:* Integrate the viscous term by parts · *tex:* $\int_\Omega\frac2{Re}(\nabla\cdot\mathbf D)\cdot\tilde{\mathbf u}\,d\Omega=\frac2{Re}\oint\tilde{\mathbf u}\cdot\mathbf D\cdot\mathbf n\,dA
     -\frac2{Re}\int_\Omega\mathbf D:\nabla\tilde{\mathbf u}\,d\Omega$ · *why:* The tensor version of D09 steps 3–5 (product rule + divergence theorem). · *plain:* one
     derivative moved from u to ũ.
  7. *did:* Drop the surface term · *tex:* $\frac2{Re}\oint\tilde{\mathbf u}\cdot\mathbf D\cdot\mathbf n\,dA=0$ · *why:* ũ = 0 on the Dirichlet parts of the boundary;
     elsewhere it would be a traction condition (natural, as in D10). · *plain:* no boundary contribution where velocity is given.
  8. *did:* Use the symmetry of D · *tex:* $\mathbf D[\mathbf u]:\nabla\tilde{\mathbf u}=\mathbf D[\mathbf u]:\mathbf D[\tilde{\mathbf u}]$ · *why:* A symmetric tensor contracted
     with any tensor sees only that tensor's symmetric part; the antisymmetric (rotation) part of ∇ũ gives zero. · *plain:* only the
     strain of the test function matters.
  9. *did:* Collect the momentum statement · *tex:* $\int_\Omega\Big(\frac{\partial\mathbf u}{\partial t}+\mathbf u\cdot\nabla\mathbf u-\mathbf g\Big)\cdot\tilde{\mathbf u}\,d\Omega+
     \frac2{Re}\int_\Omega\mathbf D[\mathbf u]:\mathbf D[\tilde{\mathbf u}]\,d\Omega-\int_\Omega p\,(\nabla\cdot\tilde{\mathbf u})\,d\Omega=0$ (10.134) · *why:* Steps 3, 6–8 in step 1, all
     terms moved left (signs flip). · *plain:* the weak momentum equation.
  10. *did:* Test continuity with p̃ · *tex:* $\int_\Omega\tilde p\,\nabla\cdot\mathbf u\,d\Omega=0$ (10.135) · *why:* Multiply ∇·u = 0 by any pressure test function
      and integrate — no derivative to move. · *plain:* the weak continuity equation.
  11. *did:* Notice the shared pairing · *tex:* $b(\mathbf v,q)=-\int_\Omega q\,\nabla\cdot\mathbf v\,d\Omega\ \ \text{appears as}\ b(\tilde{\mathbf u},p)\ \text{and}\ b(\mathbf u,\tilde p)$ ·
      *why:* The same bilinear form couples pressure to velocity in (10.134) and velocity to pressure in (10.135) (with the sign chosen as in
      (10.159)); discretised, it gives B and Bᵀ in the saddle-point system (10.137). · *plain:* pressure is the multiplier of continuity.
- **Result.** (10.134)–(10.135) — *in words:* in weak form, viscosity acts through strain rates, the pressure multiplies the divergence of
  the test velocity, and continuity is tested by pressure functions — the structure of (10.137).
- **Check.** Units (dimensionless): every integrand ~ U²/L × U over the domain ✓. `ch10.weak_ns_identity_sympy()` verifies the identity of step
  4–8 for polynomial fields vanishing on the boundary; the Cartesian expansion gives (10.156)–(10.159) (`weak_ns_components_sympy`).
- **What it means.** The pressure's role as a Lagrange multiplier (P239) is now explicit; whether the discrete B can 'see' every pressure is
  the LBB question of C13.
- **Traps.** The sign of the pressure term; dropping ∇(∇·u) without saying ∇·u = 0 (the form changes for compressible flow); writing
  D : ∇ũ = D : ∇ũᵀ without the symmetry argument.

### D23 · Observed order and Richardson extrapolation $f_0\approx f_1+\frac{f_1-f_2}{r^p-1}$ from the error model $\lVert e\rVert\le K_e\Delta x^a$ (10.15) — ★★, 8 steps, in C15 (notebook · `lid_driven_cavity`)
- **Goal.** From three runs on grids refined by a factor r, measure the order p of the code and estimate the grid-independent answer (the
  book only says grid-independence gives confidence).
- **Start.** The error model (10.15) for one output f (a drag coefficient, a velocity minimum): $f(h)=f_0+K_eh^p+\dots$ — *in words:* the
  computed value is the true one plus an error shrinking like hᵖ.
- **Plan.** (1) Three grids h, rh, r²h. (2) Differences cancel the unknown f₀; their ratio cancels K_e. (3) Logs give p. (4) Two equations
  give f₀.
- **Tools.** Log–log slopes (Ch. 1 P13) · the convergence model (N10).
- **Assumptions.** All three grids in the *asymptotic range* (the leading term dominates, step 1); the same refinement ratio r; Δt refined
  consistently or its error negligible.
- **Steps.**
  1. *did:* Write the model on three grids · *tex:* $f_1=f_0+K_eh^p,\quad f_2=f_0+K_e(rh)^p,\quad f_3=f_0+K_e(r^2h)^p$ · *why:* (10.15) with one dominant term; h₁ = h
     is the finest grid. · *plain:* three runs, one unknown truth.
  2. *did:* Difference neighbouring grids · *tex:* $f_2-f_1=K_eh^p(r^p-1),\qquad f_3-f_2=K_eh^pr^p(r^p-1)$ · *why:* Subtracting removes f₀; factor out
     $K_eh^p$. · *plain:* the changes between grids.
  3. *did:* Take their ratio · *tex:* $\frac{f_3-f_2}{f_2-f_1}=r^p$ · *why:* $K_eh^p(r^p-1)$ cancels (nonzero if the grids differ and the scheme is not
     exact). · *plain:* each refinement shrinks the change by rᵖ.
  4. *did:* Solve for p · *tex:* $p=\frac{\ln\big[(f_3-f_2)/(f_2-f_1)\big]}{\ln r}$ · *why:* Take natural logs of both sides (P13 — a slope on log axes); the
     ratio must be positive. · *plain:* the observed order of accuracy.
  5. *did:* Compare p with the scheme's order · *tex:* $p_{observed}\approx p_{formal}\ \ (=2\ \text{for MAC in space})$ · *why:* Agreement is the evidence the
     code solves its equations as designed; disagreement means a bug or grids outside the asymptotic range. · *plain:* the code check.
  6. *did:* Use the two finest grids for $K_eh^p$ · *tex:* $K_eh^p=\frac{f_2-f_1}{r^p-1}$ · *why:* The first equation of step 2 divided by rᵖ − 1. · *plain:* the
     size of the error on the finest grid.
  7. *did:* Subtract it from f₁ · *tex:* $f_0\approx f_1-K_eh^p=f_1+\frac{f_1-f_2}{r^p-1}$ · *why:* From step 1, f₀ = f₁ − K_e hᵖ; this is Richardson's
     extrapolation (`tools.convergence.richardson`). · *plain:* the grid-independent estimate.
  8. *did:* Attach an error bar · *tex:* $\mathrm{GCI}_{fine}=F_s\,\frac{\lvert f_1-f_2\rvert}{(r^p-1)\lvert f_1\rvert},\quad F_s=1.25$ · *why:* Roache's grid convergence
     index turns step 6 into a relative uncertainty with a safety factor; meaningful only if the ratio of step 3 is positive and steady
     (monotone convergence). · *plain:* how far the fine-grid answer may be from the truth.
- **Result.** $p=\ln\frac{f_3-f_2}{f_2-f_1}/\ln r$ and $f_0\approx f_1+\frac{f_1-f_2}{r^p-1}$ — *in words:* three grids measure the order and remove the
  leading error.
- **Check.** f = 1 + h² at h = 0.1, 0.05, 0.025: ratio 4, p = 2, f₀ = 1.000000 ✓ (C15 tiny example); units: p dimensionless, f₀ in f's units ✓.
  `grid_convergence_index` and `ch10.richardson_three` return the same.
- **What it means.** The verification habit the chapter ends on (N111): refine, measure p, extrapolate. It fails for oscillatory
  convergence (negative ratio), for grids too coarse (ratio drifting), and when a singularity (the cavity's lid corners, the block's
  corners) lowers the local order.
- **Traps.** Using two grids only (p must then be assumed); mixing refinement ratios; reading a p far from the formal order as 'good
  enough'; forgetting that Δt must be refined too for time-dependent outputs.

---
