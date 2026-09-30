# Chapter 10 — Computational Fluid Dynamics: curation
(from `analysis/ch10.md` — 143 inventory rows (#1–#143), 199 equation numbers (10.1)–(10.199), 48 derivations in §2b (a-d1…a-d48),
10 explainer seeds; `book.yaml` ch10 (6 sections, 5 viz seeds) and `policy.tier_a_full_treatment: [12, 18]`, `coverage: exhaustive`.
Curated 2026-09-30, concept-curator. Read: `knowledge/CUMULATIVE.md`, `concept_map.md`, `primers.md` (P01–P220 + P218a; **ch10 primers
start at P221**), `notation.md`, `viz_patterns.md`, `analysis/ch09_curation.md` (format), the `interactive-viz` skill §4–§5.)

Counts: A 15 · B 105 · C 23 · RECAP 13 · SKIP 2 · derivations written out 23 (★ 4 ★★ 18 ★★★ 1) · demoted to statements 21

Tier words (parsed by `tools/nbkit.py`): **CORE 15 · NOTE 113 · RECAP 13 · SKIP 2** = 143 rows; **DERIVATION 23**.
Depth: A 15 (all CORE) · B 105 (95 NOTE + 10 RECAP) · C 23 (18 NOTE + 3 RECAP: R01, R08, R10 + 2 SKIP: S01, S02).
**Reconciliation with analysis §2:** the 143 rows `#1…#143` each appear exactly once in §2 below (last column). A = 15 of 143 rows (10.5 %).

**Decisions that shape this chapter.**
1. **Fifteen A items from the analyst's sixteen.** The analyst's set is kept with these merges: (a) Noye's region (#23), the
   error equation (#18), ∣G∣² (#22) and BTCS stability (#27) ride inside one von Neumann block **C04** whose A row is the
   amplification factor (#20) — they are one chain of algebra (D04–D07) and one explainer; (b) the Lax equivalence theorem (#28) is a
   B item inside **C05** (the book states it without proof; its measured demo — FTCS at β = 0.50 vs 0.51, upwind at C = 0.9 vs 1.1 —
   closes §10.2); (c) hat functions (#42) and the FE = FD interior stencil (#44) are B items of the Galerkin block **C07** (their
   derivation D12 is written out there); (d) the steady convection–diffusion layer, the centred-scheme wiggles and upwinding's
   numerical diffusion are **one** A block **C09** (#63) — one picture (exact layer vs two schemes), one explainer, three derivations;
   (e) operator splitting (#79–#80) and the Θ-scheme (#95) are B items of the projection block **C11**; (f) the staggered grid, the
   discrete Poisson equation and the checkerboard are **C12** (#84). Not promoted (next candidates, all B inside C13): the P2
   triangle and quadrature (#131, #134) and the Newton linearisation (#125) — the FE cylinder is script-only and cached, so its
   machinery is stated, coded and tested but not derived.
2. **SEEN rows are RECAP** (thirteen, R01–R13, inventory order): the space–time grid, Taylor series, centred second derivative,
   Fourier modes, the diffusion limit β ≤ ½, incompressible and dimensionless Navier–Stokes, the conservative convective term,
   compressible NS with μ_v = 0, the strain-rate tensor, C_D/C_L histories, cylinder regimes, and the Strouhal number. Each reuses the
   earlier chapter's function; where the row adds a new instance (pure convection unstable, the block's force histories, the
   confined St) the new part is stated.
3. **SKIP is minimal** — S01 exercises (10.1–10.6), S02 literature. The analyst's prose SKIP rows hold new ideas and are C items with a
   pointer; the heated-rod series (10.199) is a B item (our own wall temperature, used as the exact field of the C05 Lax demo; exercise
   parameters and answers never printed).
4. **Derivations (D01–D23).** 27 of the 48 analysis derivations are written out, merged into 23 D rows (a-d6+a-d7, a-d9+a-d10,
   a-d15+a-d16, a-d21+a-d22, a-d30+a-d32, a-d47 into D23) plus one of ours (D19, the continuous projection Poisson equation that the
   book never writes). 21 are demoted to statements (§4c): the whole 2-D FE cylinder machinery, the MacCormack coefficient
   bookkeeping, the wall-density closures, the MAC stability limit and the Θ-scheme order (both ★★★, checked numerically / by sympy
   in a code cell instead). The one ★★★ written out (D18, MacCormack = Lax–Wendroff, second order, stable for C ≤ 1) is shown in
   `upwind_cfl_advection`.
5. **Printed slips (analysis §9 R1–R12) become named callouts** "book prints X, correct is Y" (table in §8), each with a planted wrong
   variant a test must fail where it is computable (R1 slopes, R5 corrector sign, R6 (10.172) coefficient, R11 upwind side for u < 0).
6. **Climate hooks are named where they are real:** C05 (the CFL condition sets the time step of every explicit weather model;
   gravity-wave CFL √(gH)Δt/Δx ≤ 1 with ch07 speeds), C09 (upwind tracer advection = implicit diffusion 0.5 R_cell D), C11 (operator
   splitting dynamics/physics in GCMs; projection = the Boussinesq/anelastic pressure solve of ch04 C13), C12 (the MAC staggering is
   the Arakawa C-grid of ocean and atmosphere models), C10 (a sound speed kept artificially low ≈ "reduced speed of sound" models).
7. **Book values stay private** (`tests/book_values_ch10.json`): Moore's-law time, 37 %/13.5 % as quoted, σ, θ = 0.2929, all §10.5
   geometry, M, grids, C_D/C_L, periods and St, exercise inputs. Ghia (1982) Table I and Hou et al. (1995) centres are public
   published data, cited to their authors and stored in `reference/ch10/`.

## 1. Teaching order (A IDs grouped by book section, B/C IDs under each; one sentence each: "once you see X, Y follows")
Order = book order. Two prerequisite notes: the cell Péclet number is stated at the end of §10.2 (N19 inside C05) and worked out in
C09; the continuous projection (C11) is taught before the staggered grid (C12), as in the book. RECAP rows are reminders inside the
block named as their A parent.

**§10.1 Introduction** — no A item of its own: N01–N04 open the C01 block (what CFD predicts, the four error sources, the method families).

**§10.2 Finite-Difference Method**
- **C01** Finite-difference stencils from Taylor series and their order of accuracy (10.4)–(10.8): once you add and subtract the Taylor expansions of the neighbours, every stencil comes with its leading error term, and "first order" or "second order" is simply the power of Δx that survives — the ruler every later scheme is measured with.
  - B: N02, N05, N07
  - C: N01, N03, N04   ·   RECAP: R01, R02, R03
- **C02** Explicit FTCS scheme for the convection–diffusion equation (10.9)–(10.11), with BTCS as its implicit twin (10.12)–(10.13): once you replace each derivative in (10.1) by a stencil, the PDE becomes an update rule with two numbers α = uΔt/(2Δx) and β = DΔt/Δx² — explicit (one formula per node) or implicit (one tridiagonal solve per step).
  - B: N06, N08, N09
  - C: –   ·   RECAP: –
- **C03** Consistency and the truncation error of FTCS (10.16)–(10.17), with convergence and its rates (10.14)–(10.15): once you put the exact solution into the scheme, what is left over is the truncation error E ~ Δt/2 T_tt + uΔx²/6 T_xxx − DΔx²/12 T_xxxx; it vanishes as Δx, Δt → 0 (consistency), but that alone does not make the computed numbers converge.
  - B: N10
  - C: –   ·   RECAP: –
- **C04** Von Neumann stability: the amplification factor G (10.19)–(10.24) and the FTCS region 0 ≤ 4α² ≤ 2β ≤ 1 (10.25)–(10.28): once you see that a linear scheme multiplies each Fourier mode of the round-off by the same complex number G every step, stability is ∣G∣ ≤ 1 for every wavenumber — which gives Noye's region for FTCS, β ≤ ½ for pure diffusion, "never" for pure convection, and "always" for BTCS.
  - B: N11, N12, N13, N14, N15, N17
  - C: –   ·   RECAP: R04, R05
- **C05** Upwind differencing and the CFL condition uΔt/Δx ≤ 1 (10.29)–(10.30), closed by the Lax equivalence theorem: once you take the convective difference from the side the flow comes from, the scheme is stable exactly when the characteristic through the new point starts inside the stencil (CFL); and consistency + stability is precisely what buys convergence (Lax).
  - B: N16, N18, N19, N113
  - C: –   ·   RECAP: –

**§10.3 Finite-Element Method**
- **C06** Weak (variational) form of the transport problem and natural vs essential boundary conditions (10.34)–(10.38): once you multiply by a test function and integrate the diffusion term by parts, the Neumann condition enters by itself (natural) while the Dirichlet one must be built into the trial space (essential) — the starting line of every FE method.
  - B: N21, N22
  - C: N20   ·   RECAP: –
- **C07** Galerkin method with hat functions → M ḋ + K d = F (10.39)–(10.63): once you expand T in hat functions and demand the weak form for every hat, you get a small ODE system whose rows, on a uniform mesh, are the centred FD stencils with a ⅙–⅔–⅙ averaged time derivative.
  - B: N23, N25, N26, N27, N28, N30, N31, N32, N33
  - C: N24, N29   ·   RECAP: –
- **C08** Element matrices and assembly (10.64)–(10.78): once you compute one 2 × 2 mass and stiffness block on a parent element and scatter-add it into the global matrix, any mesh (and in 2-D any geometry) is assembled the same way — this is why FE handles complex shapes.
  - B: N34, N35, N36, N37, N38
  - C: N39   ·   RECAP: –

**§10.4 Incompressible Viscous Fluid Flow**
- **C09** Convection-dominated flow: the exact boundary layer, centred-scheme wiggles for R_cell > 2 and upwinding's numerical diffusion 0.5 R_cell D (10.84)–(10.94): once you solve the centred scheme exactly (T_j = (r^j − 1)/(r^n − 1)), wiggles appear precisely when r = (1 + R_cell/2)/(1 − R_cell/2) turns negative, and upwinding cures them only by adding a fake diffusivity 0.5 R_cell D.
  - B: N42, N43, N44, N45, N46, N47, N48, N49
  - C: N40, N50   ·   RECAP: R06, R07
- **C10** Weakly compressible Navier–Stokes and the MacCormack predictor–corrector (10.95)–(10.110): once you let the density change a little (p = c²ρ, errors O(M²)), the pressure has its own time derivative and a forward-then-backward two-stage explicit step advances everything — second order, stable for a CFL number ≤ 1 built on the sound speed.
  - B: N51, N52, N53, N54, N55, N56, N57, N58
  - C: N59   ·   RECAP: R09
- **C11** Operator splitting and the projection method (10.111)–(10.118): once you split one time step into "convect and diffuse" then "add the pressure gradient that removes the divergence", the pressure is the solution of a Poisson equation and the correction is a pure gradient (curl-free) — the MAC/projection idea used by every incompressible and Boussinesq code.
  - B: N41, N60, N61, N62, N63, N72, N74
  - C: –   ·   RECAP: R08
- **C12** The staggered (MAC) grid: discrete continuity, the pressure Poisson equation and the checkerboard (10.119)–(10.128): once you put p at cell centres and velocities on faces, the discrete divergence and gradient are compact two-point differences, the Poisson equation needs no pressure boundary condition, and a zigzag pressure can no longer hide from the momentum equation.
  - B: N64, N65, N66, N67, N68, N69, N70, N71, N73
  - C: –   ·   RECAP: –
- **C13** Mixed finite elements and the LBB (inf–sup) condition: Taylor–Hood P2–P1 (10.134)–(10.137): once you see the pressure as a Lagrange multiplier for continuity, the velocity space must be rich enough to "feel" every pressure mode — equal-order elements fail (the FE version of the checkerboard) and one-order-higher velocities pass.
  - B: N75, N76, N93, N94, N95, N97, N98, N99, N102, N104, N105, N106, N108, N109, N110
  - C: N96, N100, N101, N103, N107   ·   RECAP: R10, R12, R13

**§10.5 Three Examples**
- **C14** The lid-driven cavity benchmark: MAC projection (and MacCormack) against Ghia's centreline data (Figs. 10.6–10.8): once your solver reproduces the primary eddy and the centreline velocity of a published benchmark within a stated tolerance, you have evidence — not proof — that the code solves the equations you meant.
  - B: N78, N79, N80, N81, N82, N83, N84, N85
  - C: N77   ·   RECAP: –
- **C15** Grid convergence, Richardson extrapolation and the verification checklist (Fig. 10.13, §10.6): once you run three grids and the differences shrink like h^p, the observed order checks the code and Richardson extrapolation estimates the grid-independent answer — the habit the chapter ends on and every later simulation needs.
  - B: N86, N87, N88, N89, N90, N91, N111
  - C: N92, N112   ·   RECAP: R11

**§10.6 Concluding Remarks** — no A item of its own: N111 (the verification checklist) and N112 (other methods named) sit inside C15.

**End matter** — S01 (Exercises 10.1–10.6), S02 (literature and supplemental reading).

## 2. Chapter map (depth) — every inventory row exactly once
Rows are in analysis §2 order (inventory row number in the last column). IDs: `C` = A/CORE (teaching order), `N` = B or C NOTE (inventory
order), `R` = RECAP, `S` = SKIP. No ASCII pipe inside cells (∣G∣ is written with U+2223).

| ID | Item | § | Depth | Tier | A parent | Reason (A) / treatment (B) / pointer (C, SKIP) / source chapter (RECAP) | Row |
|---|---|---|---|---|---|---|---|
| N01 | Def: computational fluid dynamics — quantitative flow predictions by computer from the conservation laws | 10.1 | C | NOTE | C01 | named in the chapter opener; pointer: the laws are ch04, the methods C01–C13 | #1 |
| N02 | Four error sources: discretisation, input data, initial/boundary conditions, modelling | 10.1 | B | NOTE | C01 | stated as a four-row table; discretisation error is the one we measure (C03, C15); boundary sensitivity returns in N111 | #2 |
| N03 | Advantages of CFD; computing power doubling (Moore's law) | 10.1 | C | NOTE | C01 | named in one sentence (number private); pointer: C15 says why cheap runs still need verification | #3 |
| N04 | Method families named: finite difference, finite element, finite volume, spectral | 10.1 | C | NOTE | C01 | named; FD in C01–C05, FE in C06–C08/C13, spectral = ch07 `kdv_solve` (N112) | #4 |
| N05 | Eq. (10.1) 1-D convection–diffusion of T(x, t) with constant u and D, the model problem | 10.2 | B | NOTE | C01 | stated with the exact travelling-diffusing Gaussian `FD.advected_gaussian` as the test field of C02–C05 | #5 |
| N06 | Eq. (10.2) Dirichlet T(0,t) = g and Neumann ∂T/∂x(L,t) = q; Eq. (10.3) initial condition | 10.2 | B | NOTE | C02 | stated; Neumann by a ghost node (primer) in `FD.solve_transport_1d`; one number (q = 0: insulated end) | #6 |
| R01 | Fig. 10.1 uniform space–time grid x_i = iΔx, t_n = nΔt, T_i^n | 10.2 | C | RECAP | C01 | ch01 C12 (P21) FTCS grid; our schematic `scripts/ch10_grids.py` | #7 |
| R02 | Eq. (10.4), Eq. (10.5) Taylor expansions of T at the neighbours i ± 1 to O(Δx⁵) | 10.2 | B | RECAP | C01 | ch01 P26, ch03 P98 Taylor; the starting line of D01 (`FD.taylor_table_sympy`) | #8 |
| C01 | Eq. (10.6) forward, backward O(Δx) and centred O(Δx²) first-derivative stencils; order of accuracy | 10.2 | A | CORE | – | load-bearing: every scheme of the chapter is built from these stencils and judged by their order; observed-order plot is the chapter's first verification | #9 |
| R03 | Eq. (10.7) centred second derivative O(Δx²) | 10.2 | B | RECAP | C01 | ch01 C12 / ch06 5-point Laplacian; re-derived as the last move of D01 | #10 |
| N07 | Eq. (10.8) time differences: forward, backward O(Δt), centred (leapfrog) O(Δt²) | 10.2 | B | NOTE | C01 | stated (same moves as D01 in t, §4c); `FD.ode_scheme_order` on y′ = λy with one number | #11 |
| C02 | Eq. (10.9), Eq. (10.10), Eq. (10.11) explicit FTCS scheme with α = uΔt/(2Δx), β = DΔt/Δx² | 10.2 | A | CORE | – | load-bearing: the object of the consistency (C03) and stability (C04) analyses; u = 0 reduces to ch01 FTCS | #12 |
| N08 | Eq. (10.12), Eq. (10.13) implicit BTCS scheme; a tridiagonal system per step | 10.2 | B | NOTE | C02 | stated with `solve_banded` (P193 reminder) and one step traced on 5 nodes; stability in D07 | #13 |
| N09 | Def: explicit vs implicit algorithm | 10.2 | B | NOTE | C02 | two-column table (update formula vs system solve; cost vs time-step freedom) | #14 |
| N10 | Def: convergence; Eq. (10.14) error; Eq. (10.15) error ≤ K Δx^a Δt^b with rates a, b | 10.2 | B | NOTE | C03 | stated with `FD.error_norm` (rms) and a measured (a, b) = (2, 1) for FTCS; the log move is in D23 | #15 |
| C03 | Def: consistency; Eq. (10.16), Eq. (10.17) FTCS truncation error E_i^n | 10.2 | A | CORE | – | load-bearing: the modified-equation idea behind numerical diffusion (C09), MacCormack's order (C10) and Lax (N18) | #16 |
| N11 | Def: stability (round-off disturbances decay) | 10.2 | B | NOTE | C04 | stated with a picture: a 1e-10 kick that stays small vs grows | #17 |
| N12 | Eq. (10.18) disturbance ξ; Eq. (10.19) ξ obeys the same homogeneous scheme | 10.2 | B | NOTE | C04 | stated; its subtraction is the first move of D04; `FD.propagate_error` | #18 |
| R04 | Eq. (10.20), Eq. (10.21), Eq. (10.22) Fourier decomposition of the error, one mode e^{iπkx_i} at a time | 10.2 | B | RECAP | C04 | ch05 P142, ch07 Fourier modes; ⚠️ book convention θ = kπΔx vs e^{ikx} (`FD.fourier_mode`) | #19 |
| C04 | Eq. (10.23), Eq. (10.24) amplification factor G = g^{n+1}/g^n of FTCS (von Neumann analysis) | 10.2 | A | CORE | – | load-bearing: every stability statement of the chapter (FTCS, BTCS, upwind, MacCormack, MAC) is a statement about G | #20 |
| N13 | Eq. (10.25) stability condition ∣G∣ ≤ 1 for every wavenumber | 10.2 | B | NOTE | C04 | stated; `FD.max_amplification`, `FD.is_von_neumann_stable` | #21 |
| N14 | Eq. (10.26) ∣G∣² of FTCS with θ = kπΔx | 10.2 | B | NOTE | C04 | derived in D05 (A chain); plotted as ∣G(θ)∣ curves | #22 |
| N15 | Eq. (10.27) FTCS stable ⇔ 0 ≤ 4α² ≤ 2β ≤ 1 (Noye) | 10.2 | B | NOTE | C04 | the book cites Noye without proof: derived in D06 (A chain); the (α, β) stability map | #23 |
| R05 | Eq. (10.28) pure diffusion β ≤ ½ ⇔ Δt ≤ Δx²/(2D); pure convection never stable | 10.2 | B | RECAP | C04 | ch01 C12 `stable_time_step` (r ≤ ½); the new half (D = 0 unstable) is the last step of D06 | #24 |
| N16 | Eq. (10.29) first-order upwind scheme for pure convection | 10.2 | B | NOTE | C05 | stated with the one-line update; ⚠️ slip R11: for u < 0 the upwind side flips (`scheme="upwind"` uses sign(u)) | #25 |
| C05 | Eq. (10.30) CFL condition uΔt/Δx ≤ 1 | 10.2 | A | CORE | – | load-bearing: sets the time step of every explicit solver here (MacCormack, MAC) and of weather/ocean models (gravity-wave CFL) | #26 |
| N17 | Claim: implicit (10.13) is consistent and unconditionally stable | 10.2 | B | NOTE | C04 | the book asserts it: shown in D07 (G_BTCS) and by a run at β = 100 | #27 |
| N18 | Thm: Lax equivalence — for a well-posed linear problem, consistent + stable ⇔ convergent | 10.2 | B | NOTE | C05 | stated (the book gives no proof) with the measured demo: FTCS β = 0.50 converges, 0.51 blows up; upwind C = 0.9 vs 1.1 (`scripts/ch10_lax_demo.py`) | #28 |
| N19 | Eq. (10.31) cell Péclet number R_cell = uΔx/D ≤ 2 | 10.2 | B | NOTE | C05 | stated with one number (u = 1 m/s, Δx = 0.01 m, D = 0.005 m²/s → 2); the reason is C09 (D16) | #29 |
| N20 | Def: strong (classical) form = (10.1) + (10.2) + (10.3) | 10.3 | C | NOTE | C06 | named in one sentence at the top of C06 | #30 |
| N21 | Def: H¹; Eq. (10.32) trial space S (Dirichlet built in), Eq. (10.33) test space V | 10.3 | B | NOTE | C06 | stated in plain words ("finite slope energy"); primer on test functions and spaces | #31 |
| C06 | Eq. (10.34), Eq. (10.35), Eq. (10.36) weak form: multiply by w, integrate, integrate by parts; boundary term Dq w(L) | 10.3 | A | CORE | – | load-bearing: both FE parts (1-D C07–C08 and the mixed NS form C13) start here; essential vs natural BCs | #32 |
| N22 | Eq. (10.37), Eq. (10.38) weak ⇒ strong: Neumann natural, Dirichlet essential | 10.3 | B | NOTE | C06 | derived in D10 (A chain of C06); sympy `weak_to_strong_sympy` | #33 |
| N23 | Eq. (10.39) discrete weak problem; Eq. (10.40) lifting v^h = T^h − g^h; Eq. (10.41) Galerkin form; Eq. (10.42) bilinear form a(w, v) | 10.3 | B | NOTE | C07 | stated as the first lines of D11; `FEM1.bilinear_form` | #34 |
| N24 | Def: Galerkin (test = trial space) vs Petrov–Galerkin | 10.3 | C | NOTE | C07 | named; Petrov–Galerkin returns as SUPG (N50) | #35 |
| N25 | Eq. (10.43), Eq. (10.44), Eq. (10.45) basis N_A with N_A(0) = 0; w^h = Σ c_A N_A | 10.3 | B | NOTE | C07 | stated with `FEM1.hat_basis` | #36 |
| N26 | Eq. (10.46), Eq. (10.47) extra function N₀ with N₀(0) = 1; g^h = g N₀ | 10.3 | B | NOTE | C07 | stated; drawn on the hat figure | #37 |
| N27 | Eq. (10.48), Eq. (10.49) v^h = Σ d_A N_A; T^h = Σ d_A N_A + g N₀ | 10.3 | B | NOTE | C07 | stated with `FEM1.interpolate` | #38 |
| N28 | Eq. (10.50), Eq. (10.51), Eq. (10.52), Eq. (10.53) substitution, Σ c_A G_A = 0 ⇒ G_A = 0: n ODEs | 10.3 | B | NOTE | C07 | the middle steps of D11; sympy `galerkin_equations_sympy(3)` | #39 |
| C07 | Eq. (10.54), Eq. (10.55), Eq. (10.56), Eq. (10.57), Eq. (10.58) matrix form M ḋ + K d = F | 10.3 | A | CORE | – | load-bearing: the Galerkin recipe (expand, test, integrate) that produces every FE system, 1-D here and 2-D in C13 | #40 |
| N29 | Remark: time by an ODE solver (method of lines) or by FD in time first | 10.3 | C | NOTE | C07 | named; both routes are options of `FEM1.solve_transport` (`method="theta"` or `"solve_ivp"`) | #41 |
| N30 | Eq. (10.59), Eq. (10.60), Eq. (10.61) piecewise-linear hat functions; Fig. 10.2; N_A(x_B) = δ_AB | 10.3 | B | NOTE | C07 | stated with our Fig. 10.2 analogue and the partition-of-unity check | #42 |
| N31 | Eq. (10.62) coefficients are nodal values d_A = T^h(x_A) | 10.3 | B | NOTE | C07 | stated with one number (d_2 read off the plot) | #43 |
| N32 | Eq. (10.63) interior FE equation: centred FD stencils with a ⅙–⅔–⅙ consistent mass | 10.3 | B | NOTE | C07 | derived in D12 (A chain of C07); `FEM1.interior_stencil`; lumped-mass variant for comparison | #44 |
| N33 | Remark: Galerkin FE is equivalent to an FD method; FE's advantage is geometry | 10.3 | B | NOTE | C07 | stated with the caveat (uniform linear elements, apart from the mass matrix); steady FE = centred FD to 1e-13 | #45 |
| N34 | Element viewpoint, Fig. 10.3; Eq. (10.64) parent shapes; Eq. (10.65) map x(ξ); Eq. (10.66) inverse | 10.3 | B | NOTE | C08 | stated; the map is D13's first steps; our Fig. 10.3 analogue | #46 |
| N35 | Eq. (10.67), Eq. (10.68) shape-function slopes by the chain rule (printed index shift) | 10.3 | B | NOTE | C08 | derived in D13; ⚠️ slip R1 callout, `shape_slopes(printed=True)` must fail the chain-rule test | #47 |
| N36 | Eq. (10.69), Eq. (10.70) global matrices = sums of element contributions | 10.3 | B | NOTE | C08 | stated; scatter-add primer | #48 |
| N37 | Eq. (10.71), Eq. (10.72), Eq. (10.73) element integrals; Dq only on the last element | 10.3 | B | NOTE | C08 | stated; ⚠️ slip R2 (A ∈ {e − 1, e}); `FEM1.element_integrals` (2-point Gauss) | #49 |
| C08 | Eq. (10.74), Eq. (10.75), Eq. (10.76), Eq. (10.77) local 2 × 2 element mass, stiffness and force | 10.3 | A | CORE | – | load-bearing: element-by-element assembly is how FE meshes any geometry (the cylinder of §10.5); the book never writes the numbers | #50 |
| N38 | Eq. (10.78) local → global node map | 10.3 | B | NOTE | C08 | stated with `FEM1.connectivity(n_el)` printed for 4 elements | #51 |
| N39 | Remark: 2-D/3-D FE follow the same steps | 10.3 | C | NOTE | C08 | named; pointer: C13 and the cylinder (N93–N110) | #52 |
| N40 | Scope of §10.4: primitive variables (u, p); ψ–ω formulations; laminar only; turbulence models | 10.4 | C | NOTE | C09 | named; ψ–ω Poisson = ch06 `solve_poisson`; turbulence models → Ch. 12 | #53 |
| R06 | Eq. (10.79), Eq. (10.80) incompressible Navier–Stokes and continuity | 10.4 | B | RECAP | C09 | ch04 C08 (4.39b); `core.navier_stokes` residuals verify every solver of C11–C14 | #54 |
| R07 | Eq. (10.81) dimensionless Navier–Stokes with Re | 10.4 | B | RECAP | C09 | ch04 D30, `core.similarity.Scales`; the form every §10.4–10.5 solver integrates | #55 |
| R08 | Eq. (10.82) (u·∇)u = ∇·(uu) because ∇·u = 0 | 10.4 | C | RECAP | C11 | ch04 `conservative_to_advective_sym`; `MAC.convective_terms(form=…)` | #56 |
| N41 | Eq. (10.83) divergence-free initial velocity; wall, inflow and outflow BCs; pressure defined up to a constant | 10.4 | B | NOTE | C11 | stated as a BC table; V7 demo: p + const gives the same velocity (`MAC.pin_pressure`) | #57 |
| N42 | The two difficulties: convection-dominated oscillations; continuity as a constraint that fixes p | 10.4 | B | NOTE | C09 | stated as the §10.4 opener; difficulty 1 → C09, difficulty 2 → C10–C13 | #58 |
| N43 | Eq. (10.84), Eq. (10.85) steady 1-D convection–diffusion with T(0) = 0, T(L) = 1 | 10.4 | B | NOTE | C09 | stated; the test problem of C09 | #59 |
| N44 | Eq. (10.86), Eq. (10.87) exact solution; global Péclet number R = uL/D | 10.4 | B | NOTE | C09 | derived in D15 (A chain); overflow-safe `FD.steady_cd_exact` | #60 |
| N45 | Eq. (10.88) large-R form; Eq. (10.89) layer thickness δ/L = O(1/R); T ≈ 37 % and 13.5 % | 10.4 | B | NOTE | C09 | last steps of D15, with e⁻¹, e⁻² computed (not quoted) | #61 |
| N46 | Eq. (10.90), Eq. (10.91) centred FD of (10.84) with R_cell = R/n | 10.4 | B | NOTE | C09 | first steps of D16; `FD.steady_cd_fd(scheme="central")` | #62 |
| C09 | Eq. (10.92) δ = O(Δx/R_cell); wiggles for R_cell > 2 | 10.4 | A | CORE | – | load-bearing: the convection-dominated difficulty of all flow codes; wiggles vs numerical diffusion recurs in Ch. 12 (resolution) and Ch. 13 (tracer advection) | #63 |
| N47 | Remark: refine the grid; nonuniform grids fine in the layer | 10.4 | B | NOTE | C09 | stated with one stretched-grid run (`grid="stretched"`) that removes the wiggles at the same n | #64 |
| N48 | Eq. (10.93) first-order upwind for the steady problem | 10.4 | B | NOTE | C09 | stated; ⚠️ slip R3 (the book says "forward", it is backward/upwind); discrete r = 1 + R_cell > 0 | #65 |
| N49 | Eq. (10.94) upwind modified equation: numerical diffusivity 0.5 R_cell D | 10.4 | B | NOTE | C09 | derived in D17 (A chain); `FD.numerical_diffusivity`; the upwind solution equals the exact solution of the modified equation | #66 |
| N50 | Higher-order upwind, SUPG (Petrov–Galerkin), stabilised least-squares FE | 10.4 | C | NOTE | C09 | named; pointer: the Petrov–Galerkin idea N24; GLS also cures LBB (N76) | #67 |
| N51 | Incompressibility: continuity as a constraint, pressure as a Lagrange multiplier, instantaneous pressure ⇒ Poisson equation | 10.4 | B | NOTE | C10 | stated as the "second difficulty, three answers" overview (C10 weak compressibility, C11–C12 projection, C13 mixed FE) | #68 |
| N52 | Eq. (10.95) artificial compressibility (Chorin) ∂p/∂t + c²∇·u = 0 | 10.4 | B | NOTE | C10 | stated with the two-line motivation in words (§4c); `ch10.artificial_compressibility_channel` → plane Poiseuille, ∇·u → 0 | #69 |
| R09 | Eq. (10.96), Eq. (10.97), Eq. (10.98) 2-D compressible continuity and NS in conservation form, μ_v = 0 | 10.4 | B | RECAP | C10 | ch04 D08–D13 (Stokes' hypothesis), `navier_stokes_sym(mu_v=0)`; sympy parity in a test | #70 |
| N53 | Eq. (10.99) isothermal state p = c²ρ; low Mach ≈ incompressible, density error O(M²) | 10.4 | B | NOTE | C10 | stated with one number (M = 0.1 → ~1 %); Mach reminder of ch04 C02 | #71 |
| C10 | Eq. (10.100), Eq. (10.101), Eq. (10.102) MacCormack predictor–corrector for U_t + E_x + F_y = 0 | 10.4 | A | CORE | – | key method: the explicit engine of the cavity and block examples; the classic shock-capturing scheme of Ch. 15 | #72 |
| N54 | Eq. (10.103), Eq. (10.104), Eq. (10.105) NS predictor with centred viscous terms | 10.4 | B | NOTE | C10 | stated in a colour-coded table (flux / viscous / cross term); `MCK.ns_predictor` | #73 |
| N55 | Eq. (10.106), Eq. (10.107), Eq. (10.108) NS corrector | 10.4 | B | NOTE | C10 | stated beside N54; `MCK.ns_corrector` | #74 |
| N56 | Eq. (10.109) coefficients c₁ … c₅ | 10.4 | B | NOTE | C10 | stated; c₅ = μΔt/(12ΔxΔy) explained in one line (the cross stencil has 1/(4ΔxΔy) and the μ/3 factor) | #75 |
| N57 | FF/BB, BB/FF, FB/BF, BF/FB arrangements; cycling | 10.4 | B | NOTE | C10 | stated; overlay of FF/BB vs BB/FF on the advection test (difference O(Δx²)) | #76 |
| N58 | Eq. (10.110) semi-empirical MacCormack time-step limit with safety factor σ (book value private; ours 0.8) and Re_Δ | 10.4 | B | NOTE | C10 | stated as semi-empirical; `MCK.maccormack_dt` with one number for the 64² cavity | #77 |
| N59 | Remark: density (pressure) boundary conditions are the key issue | 10.4 | C | NOTE | C10 | named; pointer: C14 (N79–N84) and C15 (N87–N89) | #78 |
| N60 | Eq. (10.111), Eq. (10.112) operator splitting dφ/dt + A(φ) = f with A = A₁ + A₂ | 10.4 | B | NOTE | C11 | stated with a 2 × 2 non-commuting test system (`ch10.split_linear_system`, exact by `expm`) | #79 |
| N61 | Eq. (10.113), Eq. (10.114) Marchuk–Yanenko fractional steps; first order | 10.4 | B | NOTE | C11 | stated (§4c); measured order 1 on the 2 × 2 system | #80 |
| N62 | Eq. (10.115) MAC split: A₁ = convection–diffusion, A₂ = (∇p, ∇·u) | 10.4 | B | NOTE | C11 | stated; colour code A₁ teal, A₂ orange (reused in E6) | #81 |
| N63 | Eq. (10.116) explicit convection–diffusion step to u^{n+1/2} | 10.4 | B | NOTE | C11 | stated; `MAC.predictor` | #82 |
| C11 | Eq. (10.117), Eq. (10.118) projection step: implicit pressure gradient enforcing ∇·u^{n+1} = 0 | 10.4 | A | CORE | – | load-bearing: the pressure–Poisson projection behind MAC, the cavity (C14) and the Boussinesq/ocean model pressure solve (Ch. 11–13) | #83 |
| C12 | Staggered grid, Fig. 10.4: p at centres, u and v on faces | 10.4 | A | CORE | – | load-bearing: the discrete Poisson equation and the checkerboard cure; the Arakawa C-grid of climate models (Ch. 13) | #84 |
| N64 | Eq. (10.119), Eq. (10.120) staggered predictor with terms interpolated to faces | 10.4 | B | NOTE | C12 | stated; face averages documented as our choice in `MAC.predictor` | #85 |
| N65 | Eq. (10.121), Eq. (10.122) velocity correction with adjacent-pressure differences | 10.4 | B | NOTE | C12 | stated; the first move of D20; `MAC.correct` | #86 |
| N66 | Eq. (10.123) discrete continuity per cell | 10.4 | B | NOTE | C12 | stated as a flux balance of one cell (ch02 Gauss); `MAC.divergence` | #87 |
| N67 | Eq. (10.124) discrete pressure Poisson equation | 10.4 | B | NOTE | C12 | derived in D20 (A chain); `MAC.pressure_poisson_matrix`, `MAC.solve_pressure` | #88 |
| N68 | Eq. (10.125) collocated centred pressure gradient: a checkerboard p is invisible | 10.4 | B | NOTE | C12 | derived in D21 (A chain); SVD null-space count collocated vs staggered | #89 |
| N69 | Eq. (10.126) staggered face gradient | 10.4 | B | NOTE | C12 | stated; `MAC.gradient` | #90 |
| N70 | No pressure boundary condition needed on the staggered grid | 10.4 | B | NOTE | C12 | the boundary steps of D20; ⚠️ slip R4 ((10.120) should read (10.124)); V1 demo: perturbing the boundary u^{n+1/2} changes nothing | #91 |
| N71 | Algorithm: predictor → Poisson → correction; the Poisson solve is the costliest step | 10.4 | B | NOTE | C12 | stated as a pipeline diagram; `MAC.step` timed per stage | #92 |
| N72 | Projection method (Chorin 1968, Temam 1969): an irrotational correction removes the divergence | 10.4 | B | NOTE | C11 | stated as D19's result; discrete curl of the correction = 0 to round-off | #93 |
| N73 | Eq. (10.127), Eq. (10.128) MAC stability limits (Δx = Δy) | 10.4 | B | NOTE | C12 | stated (§4c, ★★★ not written out); checked by `FD.ftcs2d_max_amplification` and a Taylor–Green run at 0.95× vs 1.05× the limit | #94 |
| N74 | Eq. (10.129), Eq. (10.130), Eq. (10.131), Eq. (10.132), Eq. (10.133) Glowinski Θ-scheme; second order for θ = 1 − 1/√2 | 10.4 | B | NOTE | C11 | stated (§4c, ★★★); sympy series cell and measured orders 2 vs 1 on the 2 × 2 system | #95 |
| N75 | Eq. (10.134), Eq. (10.135) weak Navier–Stokes with velocity and pressure variations | 10.4 | B | NOTE | C13 | derived in D22 (A chain); sympy `weak_ns_identity_sympy` | #96 |
| R10 | Eq. (10.136) rate-of-strain tensor D[u] | 10.4 | C | RECAP | C13 | ch02/ch03 S_ij (`core.kinematics`); used inside D22 | #97 |
| C13 | Spurious pressure with equal order; LBB (inf–sup) condition; Taylor–Hood P2–P1 (Fig. 10.5); iso-P2/P1; GLS | 10.4 | A | CORE | – | key method: the FE answer to the incompressibility constraint; parent of the whole cylinder computation; same lesson as the checkerboard | #98 |
| N76 | Eq. (10.137) semi-discrete saddle-point system; Newton; direct or GMRES | 10.4 | B | NOTE | C13 | stated with the block sparsity pattern (`spy`) of a small P2–P1 system; GMRES named with one line | #99 |
| N77 | §10.5 overview: cavity and block by MacCormack, cylinder by mixed FE | 10.5 | C | NOTE | C14 | named; pointer: C14, C15, C13 (N93) | #100 |
| N78 | Lid-driven cavity setup, Fig. 10.6: scales, p = ρ/M², M = U/c, Re; corner singularities | 10.5 | B | NOTE | C14 | stated with our geometry figure and the scales table (M, grid private) | #101 |
| N79 | Eq. (10.138) continuity on the left wall (v = 0) | 10.5 | B | NOTE | C14 | stated; `MCK.cavity_density_bc` | #102 |
| N80 | One-sided second-order first derivative (unnumbered) | 10.5 | B | NOTE | C14 | stated; computed with C01's Taylor-matching `FD.fd_weights([0, 1, 2], 1)` | #103 |
| N81 | Eq. (10.139), Eq. (10.140) left-wall density predictor/corrector | 10.5 | B | NOTE | C14 | stated (one representative wall, §4c) | #104 |
| N82 | Eq. (10.141), Eq. (10.142) right wall | 10.5 | B | NOTE | C14 | stated as the mirror of N81 (sign flip from the backward stencil) | #105 |
| N83 | Eq. (10.143), Eq. (10.144) bottom wall | 10.5 | B | NOTE | C14 | stated as the y-version of N81 | #106 |
| N84 | Eq. (10.145), Eq. (10.146) moving lid: ∂(ρu)/∂x = U∂ρ/∂x, centred | 10.5 | B | NOTE | C14 | stated with the one new move (u = U on the lid) | #107 |
| N85 | Six-substep MacCormack cavity algorithm with coefficients a₁–a₁₁ (unnumbered) | 10.5 | B | NOTE | C14 | stated as a pipeline; ⚠️ slip R5 (stray "+" in Step 5); `MCK.cavity_coefficients`, `weakly_compressible_step` | #108 |
| C14 | Cavity results: primary-eddy centres, Fig. 10.7 streamlines, Fig. 10.8 centreline u vs published data | 10.5 | A | CORE | – | load-bearing: the benchmark that turns a solver into evidence (V5 Ghia 1982, Hou 1995); reused as the test case of every later solver | #109 |
| N86 | Square block in a channel, Fig. 10.9: geometry, sliding walls, inflow, outflow, M | 10.5 | B | NOTE | C15 | stated with our domain figure (geometry numbers private); Galilean frame change (ch03) | #110 |
| N87 | Eq. (10.147) x-momentum at the block's front face ⇒ density gradient from viscous terms only | 10.5 | B | NOTE | C15 | stated with the one-line reason (u = v = 0 kills convective terms) (§4c) | #111 |
| N88 | Eq. (10.148), Eq. (10.149), Eq. (10.150) one-sided ∂ρ/∂x, u_xx (2, −5, 4, −1), mixed derivative | 10.5 | B | NOTE | C15 | stated; weights recomputed by `FD.fd_weights` (C01 method) | #112 |
| N89 | Eq. (10.151), Eq. (10.152), Eq. (10.153), Eq. (10.154) wall density on the block's four faces; corners averaged | 10.5 | B | NOTE | C15 | stated (front face in full, the rest by symmetry); labelled heuristic; `MCK.block_wall_density` | #113 |
| N90 | Practice: double precision, ρ′ = ρ − 1, FB/BF–BF/FB cycling | 10.5 | B | NOTE | C15 | stated with a float32 vs float64 drift demo (floating-point primer) | #114 |
| R11 | Block results: C_D (4.107), C_L (4.108) histories; Re = 20 steady, Re = 100 shedding | 10.5 | B | RECAP | C15 | ch04 `drag_coefficient`, ch09 shedding; coarse run in the notebook, fine run cached (values private) | #115 |
| C15 | Fig. 10.13 grid convergence of C_D; grid-independent results give confidence | 10.5 | A | CORE | – | load-bearing: grid convergence + Richardson is the verification habit of the chapter and of every later simulation | #116 |
| N91 | Eq. (10.155) asymptotic MacCormack step Δt ≤ (σ/√2) M Δx | 10.5 | B | NOTE | C15 | stated (§4c); ⚠️ slip R12: also needs M ≪ 1; `MCK.maccormack_dt_asymptotic` vs (10.110) | #117 |
| N92 | Cylinder at Re = 1000 by MacCormack with a smoke line (Fig. 10.14) | 10.5 | C | NOTE | C15 | named; not reproduced (cost); ⚠️ slip R8 ("fourth example"); passive scalar = (10.1) in 2-D | #118 |
| N93 | FE cylinder in a channel, Fig. 10.15: domain, boundaries Γ₁–Γ₅ | 10.5 | B | NOTE | C13 | stated with our mesh figure (`FEM2.cylinder_channel_mesh`) | #119 |
| N94 | Fig. 10.16 mesh counts; P2 nodes = vertices + edges, P1 = vertices | 10.5 | B | NOTE | C13 | stated with Euler's check V − E + T = 0 on our meshes (§4c); counts private | #120 |
| N95 | Eq. (10.156), Eq. (10.157), Eq. (10.158), Eq. (10.159) Cartesian weak momentum and continuity | 10.5 | B | NOTE | C13 | stated (§4c); sympy `weak_ns_components_sympy`; the minus sign of (10.159) explained (symmetric B, Bᵀ) | #121 |
| N96 | Eq. (10.160), Eq. (10.161), Eq. (10.162) Galerkin statements on Ω^h | 10.5 | C | NOTE | C13 | named: (10.157)–(10.159) with every field replaced by its FE version | #122 |
| N97 | Eq. (10.163) time derivative at t_{n+1}: backward Euler (α = 1, β = 0) or trapezoidal (α = 2, β = 1) | 10.5 | B | NOTE | C13 | stated (§4c); orders 1 and 2 measured on y′ = λy | #123 |
| N98 | Eq. (10.164) Newton: fields = guess + correction | 10.5 | B | NOTE | C13 | stated with the Newton reminder (ch06 P152) and residual ratios from a steady Re = 40 solve | #124 |
| N99 | Eq. (10.165), Eq. (10.166), Eq. (10.167) linearised correction equations; right sides = residuals | 10.5 | B | NOTE | C13 | stated (§4c, ★★★): the one idea "drop u′·∇u′" in words; ⚠️ slip R10 (spurious star in (10.166)) | #125 |
| N100 | Eq. (10.168), Eq. (10.169) expansions in velocity and pressure shape functions | 10.5 | C | NOTE | C13 | named; LBB ⇒ velocity one order higher | #126 |
| N101 | Eq. (10.170), Eq. (10.171), Eq. (10.172) algebraic equations per node | 10.5 | C | NOTE | C13 | named; ⚠️ slip R6 (second sum in (10.172) must be v_{A′}), planted-variant test | #127 |
| N102 | Eq. (10.173), Eq. (10.174) block matrix form | 10.5 | B | NOTE | C13 | stated with the 3 × 3 block picture (`scipy.sparse.bmat`) | #128 |
| N103 | Eq. (10.175)–(10.183) global block entries (9 equations: (10.175), (10.176), (10.177), (10.178), (10.179), (10.180), (10.181), (10.182), (10.183)) | 10.5 | C | NOTE | C13 | named; coded inside `FEM2.assemble_newton_system` | #129 |
| N104 | Eq. (10.184) isoparametric map with six nodes; Fig. 10.17 | 10.5 | B | NOTE | C13 | stated with a curved element drawn by `FEM2.iso_map` | #130 |
| N105 | Eq. (10.185) quadratic triangle shape functions with ζ = 1 − ξ − η | 10.5 | B | NOTE | C13 | stated; φ_a(node_b) = δ_ab and Σφ = 1 checked (barycentric primer) | #131 |
| N106 | Eq. (10.186) element expansions; Eq. (10.187) linear pressure shapes | 10.5 | B | NOTE | C13 | stated; ⚠️ slip R7 (printed v′ for p′) | #132 |
| N107 | Eq. (10.188)–(10.197) element matrices and vectors (10 equations: (10.188), (10.189), (10.190), (10.191), (10.192), (10.193), (10.194), (10.195), (10.196), (10.197)) | 10.5 | C | NOTE | C13 | named; `FEM2.element_newton_matrices` | #133 |
| N108 | Eq. (10.198) triangle quadrature, J = x_ξ y_η − x_η y_ξ, seven-point degree-5 rule | 10.5 | B | NOTE | C13 | stated; exactness test on monomials (degree 5 exact, 6 not); Jacobian primer | #134 |
| N109 | Assembly by node map; essential BCs by row replacement or penalty | 10.5 | B | NOTE | C13 | stated; both BC methods give the same Poiseuille solution | #135 |
| R12 | Steady cylinder Re = 1, 10, 40: streamlines and vorticity (Figs. 10.18–10.19) | 10.5 | B | RECAP | C13 | ch09 C10 regimes, ch08 Stokes flow; our steady FE run (confined, trend only) | #136 |
| N110 | Re = 100 unsteady: symmetric wake → Hopf bifurcation → Kármán street (Fig. 10.20) | 10.5 | B | NOTE | C13 | stated with our cached force history (script-only run); Hopf named → Ch. 11 | #137 |
| R13 | Fig. 10.21 drag, lift, torque; Strouhal S = nd/U = 1/τ̄ | 10.5 | B | RECAP | C13 | ch04 (4.102), ch09 `shedding_frequency` (cyclic n); ⚠️ slip R9: confined vs unbounded St | #138 |
| N111 | Concluding remarks: benchmarks, mesh refinement, time-step refinement, boundary insensitivity | 10.6 | B | NOTE | C15 | stated as a five-row checklist applied to our MAC cavity (verification vs validation primer) | #139 |
| N112 | Other methods: spectral, spectral element, lattice gas/Boltzmann, DPD | 10.6 | C | NOTE | C15 | named; spectral = ch07 `kdv_solve`; Hou's cavity data (C14) are lattice Boltzmann | #140 |
| S01 | Exercises 10.1–10.6 | Ex. | C | SKIP | – | pointer: our own versions are C04 (N15), N17, N113, C08, C09, C14; no exercise text or answers in public files | #141 |
| N113 | Eq. (10.199) Fourier-series solution of a heated rod | Ex. 10.3 | B | NOTE | C05 | stated with our own wall temperature (sine series, ch08) as the exact field of the Lax demo (§4c) | #142 |
| S02 | Literature cited; supplemental reading | – | C | SKIP | – | pointer: one bibliography line; Ghia (1982), Hou et al. (1995) cited where their data are used | #143 |

## 3. Section coverage

| § | Title | A | B | C | RECAP | SKIP |
|---|---|---|---|---|---|---|
| 10.1 | Introduction | – (covered by B/C items opening the C01 block) | N02 | N01, N03, N04 | – | – |
| 10.2 | Finite-Difference Method | C01, C02, C03, C04, C05 | N05, N06, N07, N08, N09, N10, N11, N12, N13, N14, N15, N16, N17, N18, N19, N113 | – | R01, R02, R03, R04, R05 | – |
| 10.3 | Finite-Element Method | C06, C07, C08 | N21, N22, N23, N25, N26, N27, N28, N30, N31, N32, N33, N34, N35, N36, N37, N38 | N20, N24, N29, N39 | – | – |
| 10.4 | Incompressible Viscous Fluid Flow | C09, C10, C11, C12, C13 | N41, N42, N43, N44, N45, N46, N47, N48, N49, N51, N52, N53, N54, N55, N56, N57, N58, N60, N61, N62, N63, N64, N65, N66, N67, N68, N69, N70, N71, N72, N73, N74, N75, N76 | N40, N50, N59 | R06, R07, R08, R09, R10 | – |
| 10.5 | Three Examples | C14, C15 | N78, N79, N80, N81, N82, N83, N84, N85, N86, N87, N88, N89, N90, N91, N93, N94, N95, N97, N98, N99, N102, N104, N105, N106, N108, N109, N110 | N77, N92, N96, N100, N101, N103, N107 | R11, R12, R13 | – |
| 10.6 | Concluding Remarks | – (covered by B/C items inside C15) | N111 | N112 | – | S01, S02 |

Notes. N113 (Eq. (10.199), in Exercise 10.3) is taught in §10.2 inside C05 and listed there. The FE cylinder rows of §10.5 (N93–N110,
R12, R13) sit in the C13 block (mixed FE) and are shown from our steady runs and a cached unsteady run. S01/S02 close the chapter and
are listed with §10.6.

## 4. Prerequisites needing primers (concept or tool | needed by | why it is not A/B/RECAP)
Numbering continues at **P221** (the designer assigns numbers). Already primed and only *reminded* (one sentence, no new primer): P13
log–log slope and `np.polyfit`, P15 `assert np.allclose`, P16 animate, P17 slider_figure, P18 show_viz, P21/P22 finite differences
and FTCS, P25 partial derivative, P26/P98 Taylor (one and several variables), P29 lambda, P31/P94 `solve_ivp`, P38 product rule,
P40 sympy, P44 linear ODE by e^{mx} trial, P45 Euler's formula, P49 chain rule, P57 `np.linalg.solve`, P58 null space and rank,
P62 `np.einsum`, P67 `np.linalg.norm`, P68 orders of smallness, P76 `meshgrid` and the `[j, i]` layout, P77 broadcasting and
neighbour slices, P78 contour/quiver/streamplot, P79 `expm`, P80 eigenvalues, P84 fundamental theorem of calculus, P106 substitution,
P107 `np.expm1`, P117 sympy series, P142 Fourier modes and FFT Poisson, P143 Gauss–Legendre, P152 Newton's method, P153 the complex
plane in numpy (`1j`, `np.abs`), P160 Jacobi/Gauss–Seidel/SOR, P161 masks and `scipy.sparse`/`spsolve`, P192 implicit stepping,
P193 Crank–Nicolson and `solve_banded`, P214 linear stability by eigenvalues, P218a integration by parts; ch04 C02 Mach regime
(`is_incompressible_regime`), ch04 C13 Boussinesq, ch06 C12 relaxation solvers.

| Concept or tool | Needed by | Why it is not A/B/RECAP |
|---|---|---|
| Big-O notation for truncation errors (O(Δx²) means "shrinks like Δx² as Δx → 0") and "order of accuracy p" | C01, C03, D01 | new vocabulary; extends P68 orders of smallness to a remainder that is kept, not dropped |
| Floating-point round-off: machine epsilon, why the error of a difference quotient rises again for tiny h (∼ ε/h), float32 vs float64 | C01 (log–log floor), C15 (N90) | new tool; the round-off floor appears in C01's first figure |
| Ghost node for a Neumann boundary (T_{n+1} = T_{n−1} + 2Δx q keeps second order) | C02 (N06) | new numerical trick |
| Norms of an error array: rms, max, L1 (`np.sqrt(np.mean(e**2))`, `np.max(np.abs(e))`) | C03 (N10), C14, C15 | new; P67 was the length of one vector |
| Half-angle identities 1 − cos θ = 2 sin²(θ/2), sin θ = 2 sin(θ/2) cos(θ/2) | C04 (D05, D06), C10 (D18) | new trig tool inside the ∣G∣² algebra |
| Sign of a linear function on an interval: it is ≤ 0 everywhere iff it is ≤ 0 at both ends | C04 (D06) | new inequality move behind Noye's condition |
| Well-posed problem (a solution exists, is unique and depends continuously on the data) | C05 (N18 Lax) | new vocabulary, a hypothesis of the Lax theorem |
| Domain of dependence and characteristics of T_t + uT_x = 0 (the value at (x, t) came from x − uΔt) | C05 (D08), C10 | new physical picture behind CFL; ch03 material derivative is the bridge |
| Péclet number (advection vs diffusion, global R = uL/D and per cell R_cell = uΔx/D) | C05 (N19), C09 | new dimensionless group (Π groups ch01 recalled) |
| Test functions, "for every w" statements and the spaces H¹, S, V in plain words | C06, C07 | new functional-analysis vocabulary |
| Fundamental lemma of the calculus of variations (∫f w dx = 0 for every w ⇒ f = 0) | C06 (D10), C13 | new |
| Bilinear form a(w, v): linear in each slot, so sums and constants come out | C07 (D11), C13 | new |
| "Arbitrary coefficients ⇒ every bracket is zero" (choose c = unit vectors) | C07 (D11), C13 (N101) | new move in every Galerkin derivation |
| Affine map to a parent element ξ ∈ [−1, 1] and dx = (h/2)dξ | C08 (D13, D14) | new; P106 substitution recalled |
| Scatter-add assembly: `np.add.at` and COO duplicates summed (`scipy.sparse.coo_matrix(...).tocsr()`) | C08, C13 | new Python idiom |
| Linear recurrence with constant coefficients: geometric trial T_j = r^j, two roots, general solution | C09 (D16) | new; the discrete twin of P44 |
| Predictor–corrector (Heun / second-order Runge–Kutta idea: predict, re-evaluate, average) | C10 (D18) | new; P95 RK4 recalled |
| Conservation (flux) form U_t + E(U)_x + F(U)_y = 0 with state vector U and flux vectors E, F | C10 | new notation; ch04 conservation laws recalled |
| Lagrange multiplier: an extra unknown that enforces a constraint (pressure ↔ ∇·u = 0) | C10 (N51), C11, C13 | new |
| Helmholtz–Hodge decomposition: any field = divergence-free part + gradient part | C11 (D19) | new; ch02 curl of a gradient = 0 recalled |
| Operator splitting and the commutator [A₁, A₂] = A₁A₂ − A₂A₁ (why splitting errors appear) | C11 (N60, N61, N74) | new; P79 `expm` recalled for the exact reference |
| Half-index notation u_{i+1/2, j} and staggered array shapes `p[ny, nx]`, `u[ny, nx+1]`, `v[ny+1, nx]` | C12 | new; extends P76 layout |
| Singular linear systems and the compatibility condition (pure-Neumann Poisson: Σ rhs = 0, pin one value) | C12 (D20), C14 | new; extends P58 null space |
| `scipy.sparse.diags`, `scipy.sparse.kron` to build a 2-D Laplacian; `splu` factorise once, solve many | C12, C14 | extends P161 with factor caching |
| Null space by `np.linalg.svd` (count singular values below a tolerance) | C12 (N68), C13 | SVD was only glossed (ch03); first computational use |
| Saddle-point (KKT) matrix [[A, B], [Bᵀ, 0]] and why it is indefinite | C13 (N76, N102) | new |
| Generalised symmetric eigenproblem `scipy.linalg.eigh(A, B)` (the inf–sup constant) | C13 | extends P80 |
| Barycentric (area) coordinates on a triangle (ζ = 1 − ξ − η) | C13 (N105) | new |
| Jacobian determinant of a 2-D map and dΩ = J dξ dη | C13 (N104, N108) | new (2-D change of variables) |
| GMRES in one line (Krylov iteration for non-symmetric systems) | C13 (N76) | name only; `scipy.sparse.linalg.gmres` |
| Streamfunction of a computed 2-D velocity field for plotting (Poisson ∇²ψ = −ω with ψ = 0 on the walls) | C14 | reminder of ch04 ψ and ch06 `solve_poisson` with a new use |
| Reading reference data from a file (`np.loadtxt` / `pandas.read_csv` of `reference/ch10/*.csv`) and interpolating to our grid (`np.interp`) | C14 | new tool for benchmark comparison |
| Caching expensive runs (`np.savez` to `outputs/ch10/` keyed by a parameter hash; `functools.lru_cache`) | C14, C15, C13 | new Python idiom; budget rule |
| Verification vs validation (solving the equations right vs solving the right equations) | C15 (N111) | new vocabulary |

## 4b. Derivations written out (parsed by tools: ID first, CORE id in a column, ★★★ for hard, explainer slugs backticked in the LAST column)
Analysis §2b numbers appear as (a-dNN) in the Result column. Every D row stands for one `nb.derivation` block.
| ID | Result (Eq.) | CORE | Difficulty | Steps | Tools used | Traps | Shown in |
|---|---|---|---|---|---|---|---|
| D01 | forward, backward, centred first derivative and centred second derivative with their leading errors, (10.6)–(10.7) from (10.4)–(10.5) (a-d1) | C01 | ★ | 9 | Taylor series (R02, P26), big-O (primer), adding and subtracting expansions | the forward difference's leading error is −(Δx/2)T_xx (sign); in the sum the odd Δx⁵ term cancels so the second-derivative remainder is O(Δx²) after dividing by Δx², not O(Δx³) | notebook · `fd_stencil_order` |
| D02 | FTCS update (10.9)–(10.11) with α = uΔt/(2Δx), β = DΔt/Δx² (a-d3) | C02 | ★ | 6 | forward time difference (N07), D01 stencils | the 2 in α comes from the centred difference; "=" becomes "≈" once O(Δt, Δx²) is dropped; T^{n+1} only on the left (explicit) | notebook |
| D03 | truncation error of FTCS E = (Δt/2)T_tt + u(Δx²/6)T_xxx − D(Δx²/12)T_xxxx + O(Δt², Δx⁴), (10.16)–(10.17) (a-d5; the book skips the time series) | C03 | ★★ | 11 | Taylor in t and x (P98), D01 results, sympy `series` (P117) | the time Taylor series T^{n+1} = T + ΔtT_t + ½Δt²T_tt is never written by the book; sign convention (E on the left); the Δt/2 T_tt term makes FTCS first order in time | notebook · `fd_stencil_order` |
| D04 | error equation (10.19) and the FTCS amplification factor (10.20)–(10.24) (a-d6 + a-d7) | C04 | ★★ | 11 | linearity, complex exponential (P45, P153), Fourier modes (R04) | the error obeys the same scheme only because the scheme is linear; e^{iπkx_{i±1}} = e^{iπkx_i}e^{±iθ}; book convention θ = kπΔx; i is both the index and √−1 | notebook · `von_neumann_amplification` |
| D05 | ∣G∣² = (1 − 4β sin²(θ/2))² + (2α sin θ)², (10.26) (a-d8) | C04 | ★★ | 8 | Euler's formula (P45), half-angle identities (primer), modulus of a complex number (P153) | real part is 1 − 2β(1 − cos θ), not 1 − 2β cos θ; the α terms give −2iα sin θ (sign irrelevant after squaring) | notebook · `von_neumann_amplification` |
| D06 | Noye's condition 0 ≤ 4α² ≤ 2β ≤ 1 (10.27), pure diffusion β ≤ ½ (10.28) and pure convection unstable (a-d9 + a-d10; the book cites Noye without proof) | C04 | ★★ | 12 | substitution s = sin²(θ/2), sin²θ = 4s(1 − s), endpoint test for a linear function (primer) | dividing by s is allowed only for s > 0; the s → 0⁺ end is the long-wave condition 4α² ≤ 2β; D = 0 gives ∣G∣² = 1 + 4α² sin²θ > 1 for every Δt | notebook · `von_neumann_amplification` |
| D07 | BTCS unconditionally stable: G = 1/(1 + 4β sin²(θ/2) + 2iα sin θ), ∣G∣ ≤ 1 for every Δt (claim after (10.30), a-d12; the book asserts it) | C04 | ★★ | 7 | D04 moves, reciprocal of a complex number (P153) | the unknowns sit at the new level so G multiplies the stencil on the left and we divide; ∣denominator∣ ≥ its real part ≥ 1 needs β ≥ 0 | notebook · `von_neumann_amplification` |
| D08 | upwind G = 1 − C(1 − e^{−iθ}), ∣G∣² = 1 − 2C(1 − C)(1 − cos θ) ⇒ stable iff 0 ≤ C = uΔt/Δx ≤ 1, (10.29)–(10.30) (a-d11) | C05 | ★★ | 9 | D04 moves, Euler's formula (P45), domain of dependence (primer) | C = 2α (the book's α carries a ½); u < 0 flips the upwind side (slip R11); C = 1 is an exact shift, not "marginal" | notebook · `upwind_cfl_advection` |
| D09 | weak form (10.34)–(10.36): multiply by w, integrate, integrate D T_xx w by parts, insert w(0) = 0 and T_x(L) = q (a-d13) | C06 | ★★ | 9 | integration by parts (P218a), product rule (P38), test functions (primer) | sign after moving the derivative onto w; w(0) = 0 kills the Dirichlet end; the Neumann datum enters as D q w(L) on the right | notebook · `fem_hat_assembly` |
| D10 | weak ⇒ strong (10.37)–(10.38): the PDE on (0, L) and the natural condition T_x(L) = q (a-d14) | C06 | ★★ | 9 | integration by parts reversed (P218a), fundamental lemma (primer) | needs T smooth enough (H²); first choose w with w(L) = 0 to get the PDE, then any w for the boundary term; the Dirichlet condition is not recovered — it was built into S | notebook |
| D11 | Galerkin with lifting to M ḋ + K d = F, (10.39)–(10.58) (a-d15 + a-d16) | C07 | ★★ | 13 | bilinear form (primer), linear combinations (ch02), arbitrary coefficients ⇒ each bracket zero (primer) | g^h is time independent so v^h_t = T^h_t; K_AB = a(N_A, N_B) differentiates N_B in the convective part, so K is not symmetric; the g N₀ terms move to F | notebook · `fem_hat_assembly` |
| D12 | interior row on a uniform mesh (10.63): ⅙–⅔–⅙ mass, centred convection, centred diffusion (a-d17) | C07 | ★★ | 10 | piecewise integrals of hat functions (N30), D11 | ∫N_A² = 2h/3 (two halves), ∫N_A N_{A±1} = h/6; ∫N_{A±1,x}N_A = ±½; divide by h only at the end | notebook · `fem_hat_assembly` |
| D13 | parent-element map and shape slopes (10.64)–(10.68), with the printed index shift corrected (slip R1) (a-d18) | C08 | ★ | 7 | affine map (primer), chain rule (P49) | on [x_{A−1}, x_A] the local shapes are N_{A−1}, N_A (the book writes N_A, N_{A+1}); dξ/dx = 2/h, not h/2 | notebook |
| D14 | element matrices m^e = (h/6)[[2, 1], [1, 2]], k^e = (u/2)[[−1, 1], [−1, 1]] + (D/h)[[1, −1], [−1, 1]], f^e; scatter-add reproduces (10.63) (a-d19; the book never writes the numbers) | C08 | ★★ | 11 | dx = (h/2)dξ (P106, primer), D13 slopes, scatter-add (primer) | the convective block is not symmetric (row a, column b); f¹ = −g k¹_{a1} carries the Dirichlet column; A ∈ {e − 1, e} (slip R2) | notebook · `fem_hat_assembly` |
| D15 | exact steady solution T = (e^{Rx/L} − 1)/(e^R − 1), R = uL/D, large-R form and δ/L = O(1/R), e⁻¹ and e⁻² levels, (10.86)–(10.89) (a-d20) | C09 | ★ | 8 | constant-coefficient linear ODE by e^{mx} trial (P44), exponentials | the printed form overflows for R ≳ 700 (use the scaled form, P107); R → 0 gives x/L; the 37 % point is one layer thickness from the wall | notebook · `cell_peclet_wiggles` |
| D16 | centred discrete solution T_j = (r^j − 1)/(r^n − 1), r = (1 + R_cell/2)/(1 − R_cell/2), wiggles iff R_cell > 2, δ = O(Δx/R_cell), (10.90)–(10.92) (a-d21 + a-d22; the book gives only a heat-balance argument) | C09 | ★★ | 12 | linear recurrence with a geometric trial (primer), quadratic roots | r = 1 is the constant solution; at R_cell = 2 the r² coefficient vanishes (degenerate case); r < 0 alternates sign node to node | notebook · `cell_peclet_wiggles` |
| D17 | upwind modified equation u T_x = D(1 + 0.5 R_cell) T_xx, (10.93)–(10.94) (a-d23; "it can be easily shown") | C09 | ★★ | 8 | Taylor (R02), D01 | the book calls T_j − T_{j−1} "forward" (slip R3); R_cell stays fixed while reading the extra term; uΔx/2 = 0.5 R_cell D | notebook · `cell_peclet_wiggles` |
| D18 | MacCormack on linear advection = Lax–Wendroff; truncation O(Δt², Δx²); ∣G∣² = 1 − 4C²(1 − C²) sin⁴(θ/2) ⇒ stable iff C ≤ 1, (10.100)–(10.102) (a-d26; the book only states "second order") | C10 | ★★★ | 14 | predictor–corrector (primer), Taylor in t and x (P98), U_tt = a²U_xx from the PDE twice, D04 moves, half-angle identities (primer), sympy (P40, P117) | forward predictor + backward corrector average to a centred scheme; E* uses U*, which is where C² enters; forgetting U_tt = a²U_xx leaves a spurious first-order term | notebook · `upwind_cfl_advection` |
| D19 | continuous projection: divergence of (10.117) with (10.118) gives ∇²p^{n+1} = ∇·u^{n+1/2}/Δt; the correction −Δt∇p is curl-free (ours; the book writes only the discrete (10.124)) | C11 | ★★ | 9 | divergence of a gradient = Laplacian (ch02), curl of a gradient = 0 (ch02), Helmholtz–Hodge (primer) | take the divergence of both sides; ∇·u^{n+1} = 0 is imposed, not derived; the pressure's boundary condition comes from the normal velocity | notebook · `mac_projection_staggered` |
| D20 | discrete Poisson equation (10.124) from (10.121)–(10.123), boundary cells need no pressure BC, solvability Σ rhs = 0 (a-d30 + a-d32) | C12 | ★★ | 12 | half-index notation (primer), sparse matrices (P161), singular systems and compatibility (primer) | the wall-face velocity is known and not corrected, so the outside neighbour drops out (book's slip R4 names (10.120)); the matrix has a constant null vector — pin one p | notebook · `mac_projection_staggered` |
| D21 | checkerboard p = (−1)^{i+j}: zero collocated gradient (10.125), ±2/Δx staggered gradient (10.126) (a-d31) | C12 | ★★ | 7 | null space (P58), stencil arithmetic (D01) | (−1)^{i+1} − (−1)^{i−1} = 0 for every i; the same zigzag passes the collocated divergence of u; the staggered face difference sees it | notebook · `mac_projection_staggered` |
| D22 | weak Navier–Stokes (10.134)–(10.135): multiply by ũ, integrate, −∇p·ũ = −∇·(pũ) + p∇·ũ, (1/Re)∇²u = (2/Re)∇·D[u] when ∇·u = 0, D:∇ũ = D:D[ũ] (a-d35) | C13 | ★★ | 11 | divergence theorem (ch02), double dot product (ch02), strain rate (R10), D09 moves | ũ = 0 on Dirichlet boundaries kills the surface terms; D:∇ũ = D:D[ũ] only because D is symmetric; sign of the pressure term | notebook · `mixed_fe_lbb` |
| D23 | observed order from three grids and Richardson extrapolation f_h→0 ≈ f₁ + (f₁ − f₂)/(r^p − 1), from the error model (10.15) (ours + a-d47; the book only says grid independence gives confidence) | C15 | ★★ | 8 | log–log slope (P13), error model e = K h^p (N10), ratio of differences | grids must be in the asymptotic range (differences of one sign, shrinking); r is the refinement ratio; oscillatory convergence gives a negative ratio | notebook · `lid_driven_cavity` |

## 4c. Derivations demoted to statements (first column is the A parent in bold, e.g. **C20** — never a bare ID or a D id: the parser reads a bare first-cell ID as an item and blanks its tier)
| A parent | Analysis §2b item | Result stated (Eq.) | Stated in (B item) | Why not written out |
|---|---|---|---|---|
| **C01** | a-d2 time differences | forward/backward O(Δt), leapfrog O(Δt²) (10.8) | N07 | the same three moves as D01 with t for x; stated and measured on y′ = λy |
| **C02** | a-d4 BTCS rearrangement | (10.12)–(10.13), tridiagonal system | N08 | the rearrangement mirrors D02; the tridiagonal structure is shown in code (`solve_banded`) |
| **C05** | a-d45 heated-rod series | T = T_w − Σ (4T_w/((2m−1)π)) sin((2m−1)πx) e^{−D(2m−1)²π²t} with our T_w (10.199) | N113 | a Fourier sine series already derived in ch08 (Couette start-up); used only as an exact test field |
| **C10** | a-d24 artificial compressibility | ∂p/∂t + c²∇·u = 0 (10.95) | N52 | two lines in words (continuity with ρ = p/c², linearised); meaningful only at steady state |
| **C10** | a-d25 2-D compressible NS | (10.96)–(10.98) with μ_v = 0 | R09 | recap of ch04 D08–D13; sympy parity test instead |
| **C10** | a-d27 predictor/corrector coefficients | c₁–c₅ (10.109), a₁–a₁₁ (cavity algorithm) | N54, N55, N56, N85 | index bookkeeping with no new idea; each coefficient checked by a sympy stencil expansion in the tests |
| **C11** | a-d29 Marchuk–Yanenko is first order | local error O(Δt²) incl. ½Δt²[A₁, A₂], global O(Δt) (10.113)–(10.114) | N61 | a B item; measured order 1 on the 2 × 2 system |
| **C11** | a-d34 Θ-scheme second order (★★★) | R(z) − e^z = O(z³) iff θ = 1 − 1/√2, β = θ/(1 − θ) (10.129)–(10.133) | N74 | a B item; a commented sympy series cell shows the z² coefficient vanishing, plus measured orders 2 vs 1 |
| **C11** | a-d46 conservative convective term | ∇·(uu) = (u·∇)u + u(∇·u) (10.82) | R08 | recap of ch04 (product rule) |
| **C12** | a-d33 MAC stability (★★★) | ½(u² + v²)ΔtRe ≤ 1 and 4Δt/(ReΔx²) ≤ 1 (10.127)–(10.128) | N73 | cited result (Peyret & Taylor); checked numerically by a 2-D von Neumann scan and runs at 0.95× / 1.05× the limit |
| **C13** | a-d36 Cartesian weak components | (10.156)–(10.159) | N95 | expansion of D:D̃ after D22; sympy `weak_ns_components_sympy` |
| **C13** | a-d37 time-derivative variant | α = 1, β = 0 backward Euler; α = 2, β = 1 trapezoidal (10.163) | N97 | two-line reading of the trapezoidal rule; orders measured |
| **C13** | a-d38 Newton linearisation (★★★) | (10.164)–(10.167) | N98, N99 | the cylinder FE is script-only; the one idea (drop u′·∇u′, right side = residual) is stated and quadratic convergence is shown |
| **C13** | a-d39 algebraic block system | (10.168)–(10.183) | N100, N101, N102, N103 | index bookkeeping; block structure shown as a sparsity picture |
| **C13** | a-d40 isoparametric map, Jacobian, quadrature (★★★) | (10.184)–(10.198) | N104, N105, N108 | stated with exactness tests (degree 5 exact, 6 not; J = 2 × area for affine elements) |
| **C13** | a-d44 Strouhal from the period | S = nd/U = 1/τ̄ | R13 | one line (cyclic frequency); recap of ch09 |
| **C13** | a-d48 P2/P1 node counts | N_u = V + E, N_p = V, V − E + T = 0 for one hole | N94 | counting check on our meshes, stated |
| **C14** | a-d42 one-sided stencils | p. 450 formula, (10.148)–(10.150) weights | N80, N88 | computed by C01's Taylor-matching `FD.fd_weights` (the D01 method in general form) |
| **C14** | a-d43 cavity wall density updates | (10.138)–(10.146) | N79, N81, N82, N83, N84 | one wall written with its single new move (a velocity component vanishes); the others by symmetry |
| **C15** | a-d28 asymptotic MacCormack step | Δt ≤ (σ/√2) M Δx (10.155) | N91 | a limit of (10.110) in two lines; slip R12 (also needs M ≪ 1) stated |
| **C15** | a-d41 block wall density | (10.147)–(10.154) | N87, N89 | heuristic closures; sympy confirms (10.151) in the tests; front face stated in full |

## 5. Interactive explainers (5–10 + backup)
Eight explainers; every A item except C02 and C06 (both carried as "also shows" inside E2 and E5) is attached to one. Every window
shows its equations with numbers, has an Explain tab ("Explanation & interpretation", numbered sections with the reader's numbers and
a "Reading the current setting" paragraph), a synced Code tab, a Derivation tab for its D rows, a 4–8-step walkthrough and ≥ 3 check
questions. Colour code across the chapter: exact solution = muted ghost, FTCS/centred = orange, upwind = teal, implicit/BTCS = blue,
MacCormack/Lax–Wendroff = purple, unstable/wiggle = rose. Reference explainers (`interactive-viz` §5): FDV =
`forced_damped_vibrations.html`, AFE = `angular_frequency_explorer_1.html`, APS = `amplitude_phase_second_order_II_3.html`.
All JS computations are small (≤ 400 grid points in 1-D, ≤ 32² in 2-D, dense FE systems ≤ 400 unknowns).

### E1 · fd_stencil_order
- A: C01, C03 (also shows R02, R03, N07, N10, N80 one-sided stencil) · **Confusion removed:** "what does 'second-order accurate' actually mean, and why does the error stop falling for very small h?" — the order is the power of h of the first Taylor term the stencil fails to cancel; halving h divides the error by 2^p until round-off (∼ ε/h) takes over.
- **Why interactive:** dragging h on a log slider moves one dot along the log–log error curve while the stencil points and the secant/parabola through them shrink around x₀; the Taylor-term bars show live which terms cancel. A static log–log plot shows the slopes but not *why* they are 1 or 2.
- **Stage:** (1) f(x) with the stencil nodes at x₀ ± h and the difference line (secant for first derivative, parabola for the second), true tangent as a ghost; (2) log–log error vs h for the chosen stencil with slope-1 and slope-2 guides, a round-off band and the current dot; (3) (hidePortrait) term bars: coefficients of f, f′, f″, f‴, f⁗ in the stencil combination (cancelled terms greyed, the leading error term highlighted). Mode "scheme": the same bars show the three FTCS truncation terms of (10.17) for the advected Gaussian with sliders Δx, Δt.
- **Controls:** stencil (chips: forward, backward, centred, second-derivative centred, one-sided 2nd order) · h (log slider 1e-8 … 0.5) · test function (sin, exp, Gaussian) · x₀ · mode (stencil / FTCS truncation).
- **Equations:** (10.4)–(10.5) Taylor, (10.6) $\frac{T_{i+1}-T_{i-1}}{2\Delta x}+O(\Delta x^2)$, (10.7), (10.8), (10.16)–(10.17) $E=\tfrac{\Delta t}{2}T_{tt}+u\tfrac{\Delta x^2}{6}T_{xxx}-D\tfrac{\Delta x^2}{12}T_{xxxx}$.
- **Mirrors** `FD.fd_weights`, `FD.fd_derivative`, `FD.truncation_terms` (numeric, §8).
- derivations: D01, D03 · depth features: explain (the Taylor sum with your h, the error, the ratio when h halves, the round-off estimate ε/h, the reading), code (stencil + observed order), **linked views**, **terms**, **presets** (centred at h = 0.1, forward at h = 0.1, round-off h = 1e-8, one-sided at a wall), **inspector** (click a point on the error curve: stencil arithmetic with numbers), **status** ("truncation-dominated, slope 2" / "round-off-dominated") · follows: `fid_formula_lab` (term bars) + APS (numbered live explanation) · **aha:** a stencil is a weighted sum of Taylor series; the weights are chosen to kill the low terms, and the first survivor is the error — you can see its slope and where round-off beats it.

### E2 · von_neumann_amplification
- A: C04, C02 (also shows N08 BTCS, N11–N15, N17, R04, R05) · **Confusion removed:** "why does a tiny round-off error explode in one scheme and die in another, and why does β = 0.51 fail when 0.50 works?" — every Fourier mode of the error is multiplied by the same complex G(θ) each step; one mode with ∣G∣ > 1 grows like ∣G∣ⁿ however small it starts.
- **Why interactive:** the reader drags (α, β) across the edge of Noye's region and watches the G(θ) curve cross the unit circle, the max ∣G∣ pass 1 and the live march blow up from a 1e-10 kick; a static region plot cannot show the causal link from a point in the (α, β) plane to the growth of the zigzag mode.
- **Stage:** (1) complex plane: G(θ) traced for θ ∈ [0, π] with the unit circle, the current θ dot; (2) the (α, β) plane with Noye's region 4α² ≤ 2β ≤ 1 shaded and the current point; (3) a live FTCS (or BTCS/upwind) march of a random initial error on a periodic grid: ξ_i^n and log₁₀ max∣ξ∣ vs n, with the dominant wavelength labelled.
- **Controls:** α (or u) · β (or Δt) · scheme (FTCS / BTCS / upwind / Crank–Nicolson) · θ probe · transport (time steps).
- **Equations:** (10.19), (10.24) $G=(\alpha+\beta)e^{-\mathrm i\theta}+(1-2\beta)+(\beta-\alpha)e^{\mathrm i\theta}$, (10.25), (10.26), (10.27) $0\le4\alpha^2\le2\beta\le1$, (10.28), G_BTCS.
- **Mirrors** `FD.amplification_factor`, `FD.ftcs_stable`, `FD.max_amplification`, `FD.propagate_error`.
- derivations: D04, D05, D06, D07 · depth features: explain (G at your θ with α, β substituted; ∣G∣²; the worst θ; growth after n steps; the regime reading), code, **linked views** (plane + region + march on one clock), **transport**, **presets** (β = 0.50, β = 0.51, pure convection β = 0, Noye edge 4α² = 2β, BTCS β = 100), **status** ("stable: max ∣G∣ = 0.98" / "unstable: 4α² > 2β, long waves grow" / "unstable: 2β > 1, zigzag grows"), **inspector** (click a θ: the arithmetic of G) · follows: AFE (linked views on one clock, presets, "right now") + FDV (Explain tab) · **aha:** stability is a property of G alone: pure diffusion fails at the shortest wave (θ = π) when β > ½, pure convection fails at every wave because ∣G∣² = 1 + 4α² sin²θ, and the implicit scheme divides instead of multiplying, so ∣G∣ ≤ 1 always.

### E3 · upwind_cfl_advection
- A: C05, C10 (also shows N16, N18 Lax, C03's modified-equation idea, N57 FF/BB vs BB/FF) · **Confusion removed:** "why must the time step shrink with the grid, and why does one scheme smear a pulse while another makes ripples?" — the stencil must contain the point the characteristic came from (C ≤ 1); upwind loses amplitude (numerical diffusion uΔx(1 − C)/2), MacCormack/Lax–Wendroff keeps amplitude but lets short waves lag (dispersion).
- **Why interactive:** the CFL slider moves the characteristic line in the x–t stencil view across the stencil's edge at exactly C = 1 while the pulse in the main view stays clean, smears, ripples or explodes; toggling schemes on the same clock separates diffusion from dispersion — a static comparison cannot show the threshold.
- **Stage:** (1) a square pulse and a Gaussian advected round a periodic domain: exact (ghost), upwind (teal), MacCormack (purple), FTCS (orange, blows up); (2) x–t stencil diagram at the current node: the three stencil points, the characteristic through (x_i, t_{n+1}) and its foot, shaded domain of dependence; (3) (hidePortrait) ∣G(θ)∣ and relative phase speed vs θ for the selected schemes.
- **Controls:** Courant number C (0.1–1.2) · scheme toggles (upwind, MacCormack FF/BB or BB/FF, FTCS) · pulse (square / Gaussian) · cells N (25–200) · transport (revolutions).
- **Equations:** (10.29) $T^{n+1}_i=T^n_i-2\alpha(T^n_i-T^n_{i-1})$, (10.30) $u\Delta t/\Delta x\le1$, (10.100)–(10.102) predictor–corrector, Lax–Wendroff form, $\lvert G\rvert^2=1-4C^2(1-C^2)\sin^4(\theta/2)$.
- **Mirrors** `FD.transport_1d_step(scheme="upwind")`, `MCK.maccormack_advection_1d`, `FD.amplification_factor(scheme="upwind"|"lax_wendroff")`, `FD.phase_error` (§8).
- derivations: D08, D18 (★★★, all 14 steps) · depth features: explain (C, the foot of the characteristic, ∣G∣ at the pulse's dominant θ, amplitude left after one revolution, numerical diffusivity uΔx(1 − C)/2 with your numbers), code, **linked views**, **transport**, **presets** (C = 1 exact shift, C = 0.5 upwind smear, C = 0.8 MacCormack ripples on the square pulse, C = 1.05 blow-up, FTCS any C), **status** ("C = 1.05 > 1: characteristic leaves the stencil — unstable"), end-of-run summary card (amplitude and phase error after N revolutions) · follows: AFE (one clock, linked views, end-of-run card) · **aha:** CFL is about information, not accuracy — and inside the stable range every scheme still pays: first order pays in amplitude, second order in phase.

### E4 · cell_peclet_wiggles
- A: C09 (also shows N19, N42–N49, N47 stretched grid, N50 named) · **Confusion removed:** "why does a perfectly consistent centred scheme produce negative temperatures, and why is upwinding 'stable but wrong'?" — the centred discrete solution is r^j with r = (1 + R_cell/2)/(1 − R_cell/2), which is negative for R_cell > 2; upwinding keeps r = 1 + R_cell > 0 by adding a diffusivity 0.5 R_cell D.
- **Why interactive:** sliding R (or n) moves R_cell across 2 and the centred dots flip into a zigzag at exactly that point while the root r crosses zero in the second view; the "modified equation" ghost lies on the upwind dots, showing that upwind solves a different problem exactly. Static pictures at two R_cell values hide the threshold.
- **Stage:** (1) T(x) near the outflow wall: exact (10.86) ghost, centred dots (orange), upwind dots (teal), exact solution of the modified equation (dashed teal); (2) r vs R_cell for both schemes with the r < 0 band and the current dots; (3) (hidePortrait) error vs n on log–log (slopes 2 and 1) or D_num/D bars (physical vs numerical diffusivity).
- **Controls:** global Péclet R (1–500, log) · cells n (4–100) · scheme toggles · grid (uniform / stretched) · zoom to the layer.
- **Equations:** (10.84)–(10.87) $T=\frac{e^{Rx/L}-1}{e^R-1}$, (10.89), (10.91) $0.5R_{cell}(T_{j+1}-T_{j-1})=T_{j+1}-2T_j+T_{j-1}$, (10.92), (10.93), (10.94) $uT_x=D(1+0.5R_{cell})T_{xx}$, (10.31).
- **Mirrors** `FD.steady_cd_exact`, `FD.steady_cd_fd`, `FD.steady_cd_discrete_exact`, `FD.numerical_diffusivity`.
- derivations: D15, D16, D17 · depth features: explain (R_cell = R/n, r for each scheme, T at the last interior node, layer thickness vs Δx, D_num with your numbers), code, **linked views**, **presets** (R_cell = 1, 2, 4, 20; stretched grid at R_cell = 4), **status** ("R_cell = 4 > 2: centred r = −3, wiggles" / "upwind: smeared layer, D_num = 2D"), **terms** (bars: physical D and numerical 0.5 R_cell D), **inspector** (click node j: r^j arithmetic) · follows: FDV (regime status + interpretation) + `overfitting_curves` (minimal two-slider figure with a verdict line) · **aha:** wiggles are the discrete equation telling the truth — the layer is thinner than the grid; upwinding hides the message by thickening the layer.

### E5 · fem_hat_assembly
- A: C06, C07, C08 (also shows N21–N38, N32 the ⅙–⅔–⅙ row, N33 FE = FD) · **Confusion removed:** "where do finite-element matrices come from, and how is FE different from FD?" — each element contributes a 2 × 2 block computed once on the parent element; scatter-adding the blocks gives the global matrices, whose interior rows are the centred FD stencils except that the time term is averaged ⅙–⅔–⅙.
- **Why interactive:** the transport steps element by element, lighting up the element on the mesh, its two hat halves and the four matrix entries it adds; sliders change n, u, D and the FE solution rebuilds from weighted hats. The reader sees the assembly happen instead of reading a sum sign.
- **Stage:** (1) the mesh with hat functions N_A (faint), the weighted hats d_A N_A and their sum T^h (bold) against the exact steady solution; the test function w currently used is highlighted; (2) the global M and K (heat-coloured cells) with the current element's 2 × 2 block outlined; (3) (hidePortrait) nodal values: FE vs centred FD vs exact (they coincide in the steady case) and the interior row written out.
- **Controls:** elements n (2–20) · u · D (or R) · mass (consistent / lumped) · transport (element being assembled) · test-function probe A.
- **Equations:** (10.36) weak form, (10.54)–(10.58) $\mathbf M\dot{\mathbf d}+\mathbf K\mathbf d=\mathbf F$, (10.59) hats, (10.63) interior row, (10.74)–(10.77) element blocks with the numbers of D14.
- **Mirrors** `FEM1.hat`, `FEM1.element_matrices_linear`, `FEM1.assemble_1d`, `FEM1.solve_steady`, `FEM1.assembly_trace` (§8).
- derivations: D09, D11, D12, D14 · depth features: explain (element integrals with your h, u, D; the assembled row; the solution at node 2; FE vs FD difference), code (element loop with `np.add.at`), **linked views**, **transport**, **inspector** (click a matrix cell: which elements contributed and the integral arithmetic), **presets** (n = 4 hand example, pure diffusion, R = 10, lumped vs consistent) · follows: `stride_padding_playground` (formula with numbers + click a cell) + APS (numbered Explain) · **aha:** FE builds the same tridiagonal equations as FD, one small block at a time — which is exactly why it can do curved, unstructured meshes.

### E6 · mac_projection_staggered
- A: C11, C12 (also shows N41, N60–N74, N68 checkerboard) · **Confusion removed:** "what is the pressure *doing* in incompressible flow, and why put u, v and p at different places?" — the pressure is whatever it takes to make every cell's net outflow zero: one Poisson solve removes the divergence left by the predictor, the correction is a gradient (curl-free), and on the staggered grid a zigzag pressure cannot hide.
- **Why interactive:** stepping the transport through predictor → Poisson → correction shows the divergence heatmap appear and vanish to 1e-15 on the same small grid; the mode switch collocated/staggered makes the checkerboard's null mode visible (collocated gradient = 0 everywhere); clicking a cell shows its four face fluxes adding up.
- **Stage:** (1) a small staggered grid (8² to 16²): p at centres (heatmap), u and v arrows on faces, cell divergence as coloured cell borders; (2) the pressure solution and the correction arrows −Δt∇p, or in checkerboard mode the zigzag p with collocated (zero) vs staggered (±2/Δx) gradients; (3) (hidePortrait) log₁₀ max∣∇·u∣ over the stages (and the SOR residual if the iterative solver is chosen).
- **Controls:** grid n · initial field (random divergent / cavity first step / shear layer) · stage (transport: u^n, predictor u*, Poisson p, corrected u^{n+1}) · Poisson solver (direct / SOR sweeps) · mode (staggered / collocated checkerboard).
- **Equations:** (10.115)–(10.118), (10.121)–(10.124) $\nabla^2_dp^{n+1}=\frac{1}{\Delta t}\nabla_d\cdot\mathbf u^{n+1/2}$, (10.125) vs (10.126), continuous projection (D19).
- **Mirrors** `MAC.divergence`, `MAC.gradient`, `MAC.pressure_poisson_matrix`, `MAC.project`, `MAC.projection_stages` (§8), `ch10.collocated_gradient`, `ch10.checkerboard`.
- derivations: D19, D20, D21 · depth features: explain (divergence of your chosen cell before and after, Σ rhs = 0 check, the pinned value, curl of the correction, the reading), code, **linked views**, **transport** (stages), **modes** (staggered / collocated), **inspector** (click a cell: the four face fluxes and their sum), **status** ("max ∣∇·u∣ = 3e-16 after projection"), **presets** · follows: `pixels_as_parameters` (algorithm in slow motion, stages) + AFE (modes) · **aha:** the projection is one linear solve that deletes exactly the divergent part of the velocity and nothing else; the staggered grid is what lets the discrete gradient and divergence see every pressure pattern — the same C-grid ocean and atmosphere models use.

### E7 · lid_driven_cavity
- A: C14, C15 (also shows C11/C12 as the engine, N71, N78, N111 checklist, D23) · **Confusion removed:** "how do I know my CFD answer is right?" — reproduce a published benchmark within a stated tolerance, and show the error falling at the expected rate as the grid is refined; Richardson's estimate from three grids lands near the benchmark.
- **Why interactive:** the reader runs the MAC solver live on 16², 24² and 32² grids (and loads cached 64²/128² and MacCormack results), watching the primary eddy form from rest and the centreline profile approach Ghia's points; the error-vs-h view fills in dot by dot and the Richardson estimate updates. Static figures hide the convergence process.
- **Stage:** (1) the cavity: streamlines (ψ contours) over a vorticity heatmap, lid arrow, eddy-centre marker vs Ghia's centre; (2) centreline u(y) at x = ½: current solution (bold), other grids (faint), Ghia (1982) points; (3) (hidePortrait) error vs h on log–log with the observed slope and the Richardson extrapolate.
- **Controls:** Re (chips 100 / 400) · grid (16 / 24 / 32 live; 64 / 128 cached; MacCormack 64 cached) · transport (time to steady state) · show grids (toggles).
- **Equations:** (10.81), (10.116)–(10.124) (the algorithm in one line each), error model (10.15) $\lVert e\rVert\le K\Delta x^a$, Richardson (D23), the steady-state test max∣Δu∣/Δt < tol.
- **Mirrors** `MAC.cavity` (small n), `MAC.cavity_centreline`, `MAC.primary_vortex_centre`, `ch10.ghia_centreline`, `tools.convergence.richardson`.
- derivations: D23 · depth features: explain (grid, Δt from (10.127)–(10.128), steps to steady state, max deviation from Ghia, observed order, Richardson estimate with your three grids, reading), code, **linked views**, **transport**, **presets** (Re 100 on 16², 32², cached 128²; Re 400 on 32²), **status** ("steady: max∣Δu∣ < 1e-6" / "within 2 % of Ghia"), end-of-run card, highlighted benchmark table (Ghia centre vs ours, current grid row) · follows: AFE (one clock, "right now" table) + APS (numbered Explain) · **aha:** a 16² grid already has the right flow pattern but the wrong numbers; convergence at the expected order toward an independent benchmark is what earns trust.

### E8 · mixed_fe_lbb
- A: C13 (also shows N75, N76, N94, N100, N105, N108 on a single element, R12 via a cached steady cylinder picture) · **Confusion removed:** "why can't I use the same simple elements for velocity and pressure?" — with equal order, some pressure patterns are invisible to the discrete divergence (like the checkerboard of C12), the saddle-point system has spurious pressure modes and the inf–sup constant falls to zero; P2 velocities with P1 pressures keep it bounded.
- **Why interactive:** switching the pair and refining the mesh shows the Stokes cavity pressure change from a noisy spurious pattern (P1–P1) to a smooth field (P2–P1), while the inf–sup constant β_h falls or stays flat; the node view counts velocity and pressure unknowns as the mesh grows. A static pair of pictures cannot show the trend with refinement.
- **Stage:** (1) the triangulated unit square with velocity nodes (dots) and pressure nodes (squares) for the chosen pair; counts; (2) the computed pressure of a Stokes lid-driven cavity (dense JS solve on ≤ 6 × 6 squares) as a coloured triangulation; (3) (hidePortrait) β_h vs n on log axes for both pairs (JS generalised eigenvalue for small n + fluidpy table for larger n).
- **Controls:** pair (chips: P1–P1, P2–P1) · mesh n (2–6 live; up to 16 from the table) · show counting · show one element (shape functions ϕ_a on the parent triangle).
- **Equations:** (10.134)–(10.135) weak NS, (10.137) $\begin{pmatrix}\mathbf A&\mathbf B\\ \mathbf B^T&\mathbf 0\end{pmatrix}$, (10.185) ϕ_a, (10.187) ψ_b, the inf–sup condition in words and as the generalised eigenvalue $\beta_h^2=\lambda_{\min}(\mathbf B^T\mathbf A^{-1}\mathbf B,\ \mathbf M_p)$.
- **Mirrors** `FEM2.structured_square_mesh`, `FEM2.stokes_cavity` (§8), `FEM2.infsup_constant`, `FEM2.p2_node_counts`, `FEM2.p2_shape`.
- derivations: D22 · depth features: explain (unknown counts, constraint count, β_h with your n, pressure range, reading), code, **modes** (P1–P1 / P2–P1), **linked views**, **presets** (P1–P1 n = 4, P2–P1 n = 4, refinement sweep), **status** ("β_h = 0.03 and falling: LBB fails" / "β_h ≈ 0.2 bounded: stable pair"), **inspector** (click a triangle: its nodes and shape values) · follows: AFE (modes + presets) + APS (numbered Explain) · **aha:** the pressure is a Lagrange multiplier, and a multiplier is only meaningful if the velocity space can respond to it — P2–P1 is the FE twin of the staggered grid.

### B1 · operator_splitting_theta (backup)
- A: C11 (also shows N60 splitting, N61 Marchuk–Yanenko, N74 Θ-scheme) · **Confusion removed:** "why is a split scheme only first order, and what is special about θ = 1 − 1/√2?" — splitting errors come from the commutator [A₁, A₂]; Glowinski's three-substep scheme cancels the Δt² term only at θ = 1 − 1/√2.
- Stage: a 2 × 2 non-commuting linear system: phase-plane trajectories (exact by expm vs Marchuk–Yanenko vs Θ-scheme) · error vs Δt on log–log with slopes 1 and 2 · the z² coefficient of R(z) − e^z vs θ with its zero marked. Controls: θ, Δt, commutator strength. Mirrors `ch10.marchuk_yanenko`, `ch10.theta_scheme_linear`, `ch10.split_linear_system`. Derivations none · depth: linked views, presets, status, inspector. Built only if an explainer above fails review.

## 6. Python animations and interactive figures (A ID → what, why, player/figure kind)
Animations (4, `animate` + `show_animation`; ≤ 90 frames, FAST → 40; dpi 80):
| # | A | What moves | Why | Player |
|---|---|---|---|---|
| A1 | C04, C05 (N18) | FTCS on our heated rod at β = 0.50 and 0.51 side by side: profiles and log max-error vs step; a zigzag appears and grows in the right panel | "consistent but unstable does not converge" made visible (Lax) | frames (stop at the first visible zigzag) |
| A2 | C05, C10 | square pulse advected one revolution at C = 0.8: exact, upwind, MacCormack | numerical diffusion vs dispersion on one clock | video |
| A3 | C11, C12 | one MAC step in slow motion on a 16² grid: divergence heatmap of u^n, u*, then after projection; pressure appearing | the projection as a sequence of stages | frames |
| A4 | C14 | the cavity spinning up from rest at Re = 100 (64², computed once and cached): streamlines + vorticity, centreline profile approaching Ghia's points | a steady benchmark is the end of a transient | video |

Interactive figures (plotly `slider_figure` / `animate_figure`, precomputed; ≤ 40 steps):
| # | A | What the slider controls | Why |
|---|---|---|---|
| F1 | C01, C03 | h: error of five stencils vs h (log–log) with the current h marked, plus the FTCS truncation-term magnitudes at that h | published-page version of E1 |
| F2 | C04 | β at fixed α: ∣G(θ)∣ curves and the region point | published-page version of E2 (the stability edge) |
| F3 | C05 | Courant number C: upwind and MacCormack after one revolution vs exact | the smear/ripple trade at every C |
| F4 | C07, C08 | elements n: FE steady solution with its hats vs exact and centred FD | FE = FD on uniform meshes, convergence with n |
| F5 | C09 | R_cell: exact vs centred vs upwind steady profiles | published-page version of E4 |
| F6 | C14, C15 | grid (16, 24, 32, 64, 128): cavity centreline u vs Ghia points | published-page version of E7's convergence |
| F7 | C13 | mesh n: inf–sup constant for P1–P1 and P2–P1 (precomputed table) | the LBB trend with refinement |
A live `ipywidgets` cell (kernel only) pairs with F2 (free α, β, scheme) and one with F5 (free R, n, grid). Static figures accompany every
A block (our analogues of Figs. 10.1–10.8, 10.13, 10.15–10.19, all generated by our code); C15's static figure is the block C_D vs Δx
from the cached script run plus our cavity convergence plot.

## 7. From-scratch moments
Each shows a hand-written version beside the tested function, followed by `assert np.allclose(...)`.
| A | § | Hand-written | Compared with |
|---|---|---|---|
| **C01** | 10.2 | Taylor-matching weights by `np.linalg.solve` on the Vandermonde system for offsets [−1, 0, 1] and [0, 1, 2]; three stencils by slicing | `FD.fd_weights`, `FD.fd_derivative` |
| **C02** | 10.2 | FTCS update in one line of neighbour slices with periodic `np.roll` | `FD.transport_1d_step(scheme="ftcs")` |
| **C04** | 10.2 | G(θ) by complex arithmetic from (10.24) and ∣G∣² by `np.abs(G)**2` vs the closed form (10.26) | `FD.amplification_factor`, `FD.ftcs_stable` |
| **C05** | 10.2 | upwind loop over nodes for u > 0 | `FD.transport_1d_step(scheme="upwind")` |
| **C07/C08** | 10.3 | element loop with `np.add.at` building M and K for 4 elements, then the steady solve | `FEM1.assemble_1d`, `FEM1.solve_steady` |
| **C09** | 10.4 | the geometric closed form T_j = (r^j − 1)/(r^n − 1) | `FD.steady_cd_fd` (tridiagonal solve) |
| **C10** | 10.4 | MacCormack for linear advection in two lines (predictor forward, corrector backward) vs the Lax–Wendroff formula | `MCK.maccormack_step`, `MCK.maccormack_advection_1d` |
| **C12** | 10.4 | 2-D Neumann Laplacian from `scipy.sparse.diags` + `kron`, divergence of a random staggered field before/after the projection | `MAC.pressure_poisson_matrix`, `MAC.project`, `MAC.divergence` |
| **C13** | 10.4 | P2 shape functions and the 7-point rule by hand: partition of unity and exact ∫ξ²η² on the parent triangle | `FEM2.p2_shape`, `FEM2.tri_quad_7pt` |
| **C14** | 10.5 | interpolate our centreline to Ghia's y-points with `np.interp` and compute the max deviation | `MAC.cavity_centreline`, `ch10.ghia_centreline` |
| **C15** | 10.5 | observed order and Richardson from three grids in four lines | `tools.convergence.observed_order`, `richardson`, `grid_convergence_index` |
§10.1 and §10.6 have no computable A item; C03 and C06 have sympy checks (`FD.truncation_error_sympy`, `ch10.weak_form_sympy`) instead of a
second numeric version; C11 is checked through C12's projection (discrete curl of the correction = 0).

## 8. Notes for the implementer
Functions and helpers that the A items' figures and explainers need beyond analysis §4 (B/C functions optional):
1. `FD.truncation_terms(scheme, u, D, dx, dt, x, t)` → dict of the numeric truncation-error terms of (10.17) on the advected Gaussian (E1 term bars, F1); must agree with `FD.truncation_error_sympy` evaluated numerically.
2. `FD.amplification_curve(scheme, alpha, beta, n_theta)` → complex G on θ ∈ [0, π] (E2 plane, F2) and `FD.phase_error(scheme, C, theta)` → relative phase speed (E3 view 3); scalar-callable for parity rows.
3. `FD.advect_periodic(profile, C, n_cells, n_rev, scheme)` convenience wrapper (E3, A2, F3) with growth capped at 1e6 (no NaN plots).
4. `FEM1.assembly_trace(n_el, u, D)` → list of (element, local m^e, k^e, global indices) for E5's transport and inspector.
5. `MAC.projection_stages(u, v, grid, dt)` → (div_before, p, correction, div_after) for E6, A3, the C12 from-scratch cell; `MAC.cavity_centreline(state, grid)`, `MAC.primary_vortex_centre(psi, grid)`, `MAC.streamfunction(u, v, grid)`; `MAC.cavity(..., save_every=)` history for A4.
6. `ch10.ghia_centreline(Re)` reading `reference/ch10/ghia1982_table1.csv` (public, cited: Ghia, Ghia & Shin 1982, Table I; cross-check 3 values against a second transcription before use) and `ch10.hou_centres()` (Hou et al. 1995, public).
7. `FEM2.structured_square_mesh(n, pair)` and `FEM2.stokes_cavity(n, pair)` (small dense-capable Stokes solve) for E8 parity; `FEM2.infsup_table(pairs, n_list)` precomputed and written to `reference/ch10/infsup_table.csv` (our data) for E8 view 3 and F7.
8. JS tables for the explainers (our data, labelled as ours): cavity centreline and ψ at 64² and 128² for Re = 100, 400 (MAC) and 64² MacCormack; the inf–sup table; each embedded with ≤ 4 significant figures and a parity row against the fluidpy function at a small n.
9. `ch10.book_slips()` → table (slip id, printed form, correct form, evaluator) for the callouts and the planted-variant tests.
10. `tools/convergence.grid_convergence_index(f1, f2, f3, r, p=None)` (analysis §5) used by C15.

Printed-slip callouts (shown as "book prints X, correct is Y" or ⚠️ Common confusion):
| slip | Where shown | Callout |
|---|---|---|
| R1 (10.67)–(10.68) slopes | C08 / N35, D13 | labels N_A, N_{A+1} on [x_{A−1}, x_A]; correct N_{A−1}, N_A (planted `printed=True` fails) |
| R2 nonzero element entries | C08 / N37, D14 | "A = e or e + 1"; correct A ∈ {e − 1, e} |
| R3 "forward difference" after (10.93) | C09 / N48, D17 | it is a backward (upwind for u > 0) difference |
| R4 "(10.120)" in the no-pressure-BC argument | C12 / N70, D20 | should be (10.124) |
| R5 cavity Step 5 stray "+" | C14 / N85 | −a₁[(ρu)*_{i,j} − (ρu)*_{i−1,j}] (planted variant fails mass conservation) |
| R6 (10.172) second sum | C13 / N101 | u_{A′} must be v_{A′} (planted variant fails the Stokes test) |
| R7 (10.186) | C13 / N106 | printed v′ for the pressure expansion; correct p′ |
| R8 "fourth example" | C15 / N92 | §10.5 has three examples |
| R9 unbounded St comparison | C13 / R13 | two values (the book's first, then the literature 0.164–0.165); the computed case is confined, so qualitative only |
| R10 (10.166) spurious star | C13 / N99 | β ∂v/∂t(t_n) is known data, not an iterate |
| R11 (10.29)–(10.30) sign of u | C05 / N16, D08 | upwind side and CFL use ∣u∣ with the side from sign(u) (planted variant fails for u < 0) |
| R12 (10.155) | C15 / N91 | also needs M ≪ 1, not only large grid Reynolds number |

## 9. Implementation guidance (modules, runtime budget, caching)
Consistent with analysis §4 (47 rows) and §5: `fluidpy/core/fd.py` (FD), `fluidpy/core/fem1d.py` (FEM1), `fluidpy/core/mac.py` (MAC),
`fluidpy/core/maccormack.py` (MCK), `fluidpy/core/fem2d.py` (FEM2), chapter wrappers and sympy engines in
`fluidpy/ch10_computational_fluid_dynamics.py` (re-exporting every core name for `nb.core` imports), scripts `scripts/ch10_*.py`
(all `--no-show`). Additions from this curation: items 1–10 of §8.

**Computational-cost rules (carried forward from analysis §4 rows 31, 32, 38, 44, 47 and §9; notebook < 5 min on Colab CPU with FAST):**
- **MacCormack cavity:** 256² explicit to steady state is *hours* in numpy — never in the notebook or an explainer. The notebook runs
  **64²** (≈ 20 s; FAST: 32² ≈ 3 s) at Re = 100; 128² (≈ 3 min) and Re = 400 run in `scripts/ch10_cavity.py` and are cached to
  `outputs/ch10/` (npz keyed by a parameter hash); the notebook loads the cache if present, else falls back to the 64² run and says so.
  State the resolution beside every comparison with Ghia; never type Ghia numbers into prose (load them from `reference/ch10/`).
- **MAC cavity:** 64² ≈ 10 s (FAST: 32² ≈ 2 s), 128² ≈ 60 s cached by script; the Poisson matrix is factorised once per grid (`splu`,
  cached) and each step is one back-substitution. Taylor–Green (periodic, FFT pressure) 64² < 2 s is the V1 test.
- **Block in a channel (C15):** Re = 20 at Δx = 0.1 (≈ 15 s) in the notebook only when not FAST; Δx = 0.05 (≈ 2 min) and Re = 100 fine
  (≈ 20 min) are script-only and cached; C15's convergence figure uses the cached C_D values when present and the MAC cavity
  three-grid study (16/32/64, FAST 16/24/32) otherwise, so the Richardson demonstration always runs.
- **FE cylinder (C13):** steady Re = 1, 10, 40 on a coarse mesh (≈ 5 Newton iterations × 1 s) may run in the notebook when not FAST
  (FAST: Re = 40 only, coarser mesh). The **unsteady Re = 100 run is script-only** (`scripts/ch10_cylinder_fem.py --unsteady`,
  30–60 min) and writes our own small force-history file `reference/ch10/cylinder_fe_re100_forces.csv` (our data, public) so the
  notebook plots C_D, C_L, torque and St without running it; no vorticity animation of the unsteady FE run in the notebook.
- **Re = 1000 smoke-line cylinder (N92):** not reproduced anywhere.
- **FTCS blow-up demos:** cap growth at 1e6 and stop (report the step), never plot NaN; `check_stability=False` only in demos.
- **Explainers:** JS stays within the sizes in §5 (1-D ≤ 400 nodes; MAC ≤ 32² live with a banded/direct Poisson solve factorised once;
  FE dense ≤ 400 unknowns); anything larger comes from precomputed tables with parity rows.
- **Overflow:** `steady_cd_exact` uses the scaled form (R up to 10⁴); the printed form overflows near R ≈ 700.

**Evidence for A items** follows analysis §6 (minimum: stencils V1 V2 V3; FTCS/BTCS/upwind V1 V3; truncation error V2 V3; von Neumann
V1 V3 (empirical blow-up thresholds); weak form/Galerkin/assembly V1 V2 V3; cell Péclet V1 V3; MacCormack V1 V3 (+ V5 cavity);
MAC projection V4 V1 V3 V5; mixed FE V1 V2 V3 (Poiseuille exact, Kovasznay orders, inf–sup trend); cavity V5 (Ghia, Hou) V3; grid
convergence V3 with Richardson), with planted wrong variants for slips R1, R5, R6, R11.
