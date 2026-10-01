# Chapter 11 — Instability: lesson design
(from `analysis/ch11_curation.md` (A 15 · B 131 · C 24 = 170 rows; CORE 15 · NOTE 129 · RECAP 24 · SKIP 2; 25 derivations D01–D25
(★ 7 · ★★ 13 · ★★★ 5: D10, D14, D18, D19, D24); E1–E9 + backup B1; primers from **P255**) and `analysis/ch11.md` (§2b derivations
a-D1…a-D37, §4 implementation rows 1–42, §5 the new core module `core/stability.py`, §8 benchmarks, §9 slips S1–S12 and symbol
collisions). Every equation placed here was re-read on the rendered pages `chapters/pages/ch11/` (printed = pdf − 27): (11.1) +
stability classes p503; (11.6)–(11.12) p505; (11.13)–(11.18) p506; the instability inequality, (11.19), (11.20) p507; (11.21) +
the Boussinesq set p511; (11.22)–(11.27) + Fig. 11.8 p512; (11.28)–(11.33) p513; (11.34)–(11.38) p514; (11.39)–(11.42) + Fig.
11.9 p515; the even solution, the 3 × 3 determinant and the stress-free conditions p516 (slip #1 seen on the derivative line);
(11.43), (11.44), the printed dRa/dK² p517 (slips #2, #3 seen); 27π⁴/4 p518; the EOS p519; (11.45), (11.46) p521; the ring
interchange p523; (11.47)–(11.49) p524; (11.50) p525; (11.51)–(11.53) p526; the stratified set-up p529; (11.55)–(11.57) p530
(slip #4 seen); (11.58)–(11.63) p531; (11.64)–(11.67) p532; (11.68)–(11.70) p533 (slip #11 seen); (11.71), (11.72) and the
semicircle p534 (the semicircle itself is **unnumbered**; (11.72) is $\int[U^2-c_r^2-c_i^2]Q\,dz>0$); the growth bound p535;
(11.73) p535; (11.74)–(11.78) p536; (11.79), (11.80) p537; (11.81), (11.82) p538; (11.83)–(11.86) p539; Fig. 11.21 p540;
(11.87) and the cat's eye p541; (11.88) p546; the 2-D energy equation p547; (11.89), (11.90) p553; Fig. 11.28 p554; (11.91) and
the fixed points p555. 2026-09-30, lesson-designer. **Part C is written first and is the contract both the implementer and the
builders keep.** Numbers marked *expect* were computed for this design with scratch scripts (numpy/scipy, Chebyshev collocation
with N = 40–180, the book's own determinant, `solve_ivp` DOP853); the builder compares an executed cell against them.)

**Binding conventions for every builder (analysis §9, curation decisions 1–8, §8 and §9).**
1. **Imports and aliases.** `from fluidpy import ch11_instability as ch11`; `from fluidpy.core import stability as ST`.
   **`ch11` re-exports every public name of `ST`**, so explainer parity rows write `ch11.<name>` only. Reused modules:
   `from fluidpy.core import waves as WV, laminar as LAM, boundary_layer as BL, jets as JET, stratification as STRAT,
   kinematics as KIN`; `from fluidpy import ch05_vorticity_dynamics as ch05, ch07_gravity_waves as ch07`.
2. **Scalar-callable and parity-friendly** (the shot.py evaluator has no builtins, only `np`, `math`, `ch11`): every public
   function accepts Python floats and returns a float, a complex, a tuple or a `dict` of floats; arrays only when an array goes
   in or a field/spectrum comes out. Parity `py:` rows extract real numbers with `np.real(…)`, `np.imag(…)`, `np.abs(…)`,
   dict keys and integer indices (never ndarray methods, never builtins such as `abs`, `max`, `float`).
3. **Units.** Dimensional SI only in the thin wrappers that say so (KH §11.3, `rayleigh_number`, the double-diffusive
   gradients, `taylor_number`, `ring_interchange_energy`, `rayleigh_taylor_cutoff`, `gradient_richardson` with SI inputs);
   every eigen-solver is **non-dimensional** in the scales of its section (convention 6). Default g: **G0 = 9.80665 m/s²**
   (`core.thermo.G0`); the ch07 parity rows pass `g=` explicitly because `core.waves` defaults to G_BOOK = 9.81. Densities in
   kg/m³, surface tension `surface_tension` σ_s in N/m.
4. **Symbols (overloading resolved; a ⚠️ table cell in the front matter, a one-line reminder where each is first used).**

| Book symbol | Meanings in the book | Notebook / explainer symbol | Code name |
|---|---|---|---|
| K, k, κ | Bénard horizontal wavenumber ∣K∣ (§11.4); streamwise / axial wavenumber (§11.3, §11.6–11.10); Lorenz roll wavenumber (§11.14); thermal diffusivity | $K$ (Bénard, non-dim); $k$ (all others); $k_L$ in our own Lorenz lines when it could meet $K$; $\kappa$ heat, $\kappa_s$ salt | `K`, `k`, `kappa`, `kappa_s` |
| α | thermal expansion (§11.4–11.5); $\Omega_2/\Omega_1-1$ (§11.6) | $\alpha$ thermal; $\alpha_\Omega\equiv\mu-1$ with $\mu\equiv\Omega_2/\Omega_1$ | `alpha`; `mu` (never `alpha` for rotation) |
| β | haline contraction (§11.5) (ch10: FTCS β) | $\beta$ (haline) | `beta_S` |
| σ | growth rate (11.1); surface tension (Ex. 11.1, 11.2) | $\sigma$ growth rate; $\sigma_s$ surface tension | `sigma`; `surface_tension` |
| Γ | $-d\bar T/dz$ (11.21) and $\Delta T/d$ (11.24); circulation $2\pi rU_\theta$ (§11.6); ch01 lapse rate $dT/dz$ | $\Gamma_T\equiv-d\bar T/dz$ in our lines (the book's Γ quoted in its equations); $\Gamma=2\pi rU_\theta$ circulation | code never takes Γ_T: `dT` = T_bottom − T_top; `Gamma1, Gamma2` circulation |
| U₁ | upper-stream speed (§11.3); U at the inflection point (Fjørtoft, p. 512) | $U_1,U_2$ KH streams; $U_I$ Fjørtoft | `U1, U2`; `U_I` |
| r, R | Lorenz $r=\mathrm{Ra}/\mathrm{Ra}_{cr}$; radius $R$ (§11.6), $r$ (ring interchange); bifurcation parameter $R$ (§11.14) | $r$ Lorenz; $R,R_1,R_2$ radii; $r_1,r_2$ ring radii; $R_n$ bifurcation values (written $A_n$ for the logistic map) | `r`; `R1, R2`; `r1, r2`; `A` |
| φ | perturbation velocity potential (§11.3); $\hat\psi/(U-c)^{1/2}$ (§11.7); OS stream-function amplitude (§11.8–11.9); cylindrical angle (§11.6) | $\phi_1,\phi_2$ potentials (§11.3); $\phi$ TG-transformed (C09, D18 — a callout says "not a potential"); $\phi(y)$ OS/Rayleigh amplitude; $\varphi$ angle | `phi1, phi2`; `phi` |
| ψ | §11.7 $u=\partial\psi/\partial z$, $w=-\partial\psi/\partial x$; §11.8 $u=\partial\psi/\partial y$, $v=-\partial\psi/\partial x$; §11.14 $u=-\partial\psi/\partial z$, $w=\partial\psi/\partial x$ (slip #9) | each block states its convention in a ⚠️ line | keyword `convention="11.7"/"11.8"/"11.14"` where a function returns velocities |
| c, σ | complex phase speed (§11.3, §11.7–11.11); growth rate (§11.4, §11.6, §11.14) | $\sigma=-\mathrm i\lvert\mathbf K\rvert c$ | `sigma_from_c`, `c_from_sigma` |
| Ra | $g\alpha\Gamma d^4/(\kappa\nu)$ > 0 heated below (11.21); $g\alpha d^4(d\bar T/dz)/(\nu\kappa)$ < 0 heated below (§11.5) | Ra (§11.4 sign); $\mathrm{Ra}_{\S11.5}$ in C06 and D13 | `rayleigh_number` (§11.4 sign); `thermal_rayleigh_signed` (§11.5 sign) |
| Pr | Prandtl number; in (11.91) it is the "σ" of the Lorenz literature | Pr | `Pr` |
| i | √−1 | upright $\mathrm i$ in our lines | `1j` |

5. **Book slips taught in corrected form** (analysis §9 S1–S12 = curation §8): each is a printed-vs-correct box
   (`> ⚠️ **slip #k — the book prints** … **; the correct form is** …`) where it is used, and, where computable, a planted wrong
   variant that a test must fail. In code the keys stay `S1…S12` (`ch11.book_slips()`); in every reader-facing text they are
   "slip #1 … slip #12" (no collision with recap ids R01–R24). **Two values, the book's first** only where the book's number is
   itself public-standard (Ra_c ≈ 1708 vs 1707.76; 657.5 vs 27π⁴/4 = 657.51; (11.46)'s "657"); the book's other rounded values
   (rigid–free pair, Poiseuille, Blasius, Ri of the jet, Taylor's radius ratio, wavelengths, Tollmien coefficients, exercise
   inputs) never appear — we show the published benchmark with its citation or our own computed value.

| Slip | The book prints | Correct | Taught in | Code option that must fail |
|---|---|---|---|---|
| #1 p. 489 | $(d^2/dz^2-K^2)^2W=\dots+B(q^{*2}-K^2)^2\cosh q^*z$ | $C(q^{*2}-K^2)^2\cosh q^*z$ (the matrix below it is right) | C04 (N46, D10 step 9) | `benard_determinant(…, printed=True)` uses the wrong root in column 3 (det no longer purely imaginary; root moves > 1 %) |
| #2 p. 490 | $W=A\sin(n\pi z)$ | $W=A\sin n\pi(z+\tfrac12)$ (n = 1: $\cos\pi z$) | C05 (N48, D11 step 7) | `benard_free_free_mode(z, 1, printed=True)` ≠ 0 at z = ±½ |
| #3 p. 490 | $\frac{d\mathrm{Ra}}{dK^2}=\frac{3(\pi^2+K^2)^2}{K^2}-\frac{3(\pi^2+K^2)^3}{K^4}$ | no 3 on the second term | C05 (N49, D11 step 9) | `benard_free_free_sympy()["printed_root"] == []` |
| #4 (11.55) | w-equation with $-\frac1{\rho_0}\frac{\partial p}{\partial x}$ | $-\frac1{\rho_0}\frac{\partial p}{\partial z}$ | C08 (N71, D16 step 5) | `stratified_shear_sympy(printed=True)` does not reach (11.61) |
| #5 p. 480 | "see (7.96)" | Ch. 7 (7.95), $\omega^2=gk\frac{\rho_2-\rho_1}{\rho_2+\rho_1}$ | C02 (R09) | ⚠️ box |
| #6 p. 520 | Tollmien middle branch with the square inside the bracket | $1-b(1-\eta)^2$ (continuous, matching slope) | C13 (N105) | `tollmien_profile(eta, printed=True)` fails continuity |
| #7 (11.93) | $(d^2/dR^2-k^2)^2\hat u_\varphi=-\mathrm{Ta}\,k^2\hat u_R$ | first power, as (11.51) at σ = 0 | C07 (N124) | `taylor_galerkin_Ta(…, printed=True)` differs from `taylor_marginal_Ta` (the verifier measures and reports the gap) |
| #8 p. 507 | unweighted $\int[U_{\min}-U][U_{\max}-U]dz\le0$, "(U_min − U) < 0" | weight $Q\ge0$ inside from the start; ≤ | C10 (D19 step 12) | ⚠️ box |
| #9 §11.7/11.8/11.14 | three stream-function sign conventions | not an error — each block states its own | C08, C11, C15 | `convention=` keyword |
| #10 (11.21) | $\Gamma=-d\bar T/dz$ (opposite to ch01's Kundu Γ ≡ dT/dz); §11.5 redefines Ra with $d\bar T/dz$ | both conventions shown side by side (two-row table) | C03 (N25), C06 (N55) | `rayleigh_number(dT=…)` takes T_bottom − T_top; a planted Γ = +dT/dz gives Ra < 0 |
| #11 p. 506 | $\frac d{dz}[(U-c)^2F']-k^2(U-c)F+N^2F=0$ | $-k^2(U-c)^2F$ | C10 (N80, D19 step 5) | `stratified_shear_sympy()["printed_slip11_residual"] != 0` |
| #12 (11.47), (11.50) | $\frac{\partial}{\partial R}(R\tilde u_R)+\frac{\partial\tilde u_z}{\partial z}=0$ | $\frac1R\frac{\partial}{\partial R}(R\tilde u_R)+\frac{\partial\tilde u_z}{\partial z}=0$ | C07 (R15, D14 step 4) | `taylor_perturbation_sympy()["printed_continuity_units_ok"] is False` |

6. **Non-dimensionalisations (a four-row table in C03, a one-line reminder in C07, C11, C15).** §11.4: lengths d, time
   $d^2/\kappa$, $w$ kept dimensional until $W\equiv(\Gamma_Td^2/\kappa)\hat w$; z ∈ [−½, ½]. §11.6: lengths by the gap d,
   $x=(R-R_1)/d\in[0,1]$, k means $kd$, σ means $\sigma d^2/\nu$, "d/dR" in (11.51) means d/dx. §11.8–11.11: L and U₀ per flow
   (Poiseuille: half-width and centreline speed, U = 1 − y² on [−1, 1]; Blasius: δ* and U∞; tanh/sech²: L and U₀), time L/U₀,
   pressure ρU₀². §11.14: d, $d^2/\kappa$, then time × $(\pi^2+k^2)$ and rescaled X, Y, Z (D24).
7. **Colours (text, figures, derivation terms, explainers — one meaning each; curation §5):** growing / unstable `rose`;
   decaying / stable `teal`; neutral / marginal curve `accent` purple; base state and ghosts `muted` (dashed); buoyancy and
   density `blue`; shear and kinetic energy `orange`; salt `amber`; surface tension (E2 term bar) `teal`; Reynolds-stress
   production `orange`, viscous dissipation `blue`, dE/dt `rose`/`teal` by sign; printed-slip ghosts `muted` dotted, labelled
   "as printed"; benchmark points black dots with a white halo.
8. **Every book equation is shown in full next to its number** (CLAUDE.md rule 3) — markdown, derivation steps
   ("substitute (11.14), $\phi_j=A_j(z)e^{\mathrm ik(x-ct)}$"), traps, recaps, notes, figure notes and every explainer text. The
   builder keeps an `EQ` dict with all 96 labels of ch11 (TeX from analysis §2 column 5, checked against the pages listed above;
   (11.72) as on p534) plus the earlier labels cited here ((4.10), (4.86), (4.89), (4.75), (7.95), (7.128), (8.9), (8.10),
   (9.71), (1.29)) and uses the ch10 `self_check_near` + `self_check_ctrl`. Exercises are cited as "Exercise 11.6", never as bare
   numbers; (11.92)–(11.96) sit in Exercises 11.9, 11.12, 11.13 and are shown as equations (N124, N127, N128) with **our** inputs.
9. **nbkit behaviour** (ch04–ch10 lesson): `nb.recap(...)` and `nb.section(...)` end the current CORE block, so every RECAP
   sits **before** the `nb.core` call of the block that uses it (R01–R11 before C02, R12 R13 before C03, R14–R16 before C07,
   R17–R19 before C08, R20 before C09, R21 R22 before C11, R23 R24 before C13). Notes after a section break that belong to an
   earlier block carry the parent id in bold ("**C13 (continued)**" for §11.11, "**C14 (continued)**" for §11.12–11.13). Every
   `nb.derivation` sits inside its curation CORE block with `ref=` a bare label ("11.18").
10. **Parsing.** Part E is the only part with table rows after its heading; Part F has no line starting with a table bar
    (absolute values are written \lvert … \rvert or in words). Explainer headings are exactly `### E1 · normal_mode_growth` …
    `### E9 · lorenz_attractor` and `### B1 · period_doubling_route`. Primer terms in Part A are the exact Concept text of the Part
    E rows marked "primer" (coverage_check matches the first 18 characters). "Explained by" never names an earlier chapter's
    C-number (it cites "Ch. N §N.M" or a primers.md entry).
11. **Public repo (rule 9).** The strings in `tests/book_values_ch11.json → _forbidden_public` never appear (no book run
    parameters, no rounded book values that are not public standards, no exercise inputs). Our own inputs: water layer
    d = 5 mm, ΔT = 2 K, α = 2.1 × 10⁻⁴ K⁻¹, κ = 1.4 × 10⁻⁷ m²/s, ν = 1.0 × 10⁻⁶ m²/s; air layer d = 1 cm (α = 1/293 K⁻¹,
    ν = 1.5 × 10⁻⁵, κ = 2.1 × 10⁻⁵ m²/s); air over water ρ₁ = 1.2, ρ₂ = 1000 kg/m³, σ_s = 0.074 N/m; seawater α = 2 × 10⁻⁴ K⁻¹,
    β = 7.6 × 10⁻⁴ (g/kg)⁻¹, κ_s = 1.5 × 10⁻⁹ m²/s (τ = κ_s/κ = 0.0107); Taylor–Couette radius ratio R₂/R₁ = 1.05 (ours),
    water ν = 10⁻⁶ m²/s; Tollmien breakpoint η₁ = 0.2 (ours; a, b from continuity); Lorenz Pr = 10, b = 8/3, r = 28 (Lorenz 1963,
    public). Published benchmarks (Chandrasekhar 1961; Orszag 1971; Thomas via Gallagher, Griffiths & Stephen 2016; Jordinson 1970;
    Tatsumi & Kakutani 1958; Michalke 1964; Drazin & Reid 1981; Lorenz 1963; Feigenbaum 1978) live in `reference/ch11/` and are
    cited where used. Figures are generated by our code (every "Fig. 11.n" is remade, none copied; photographs are named only).
12. **Honesty labels.** (11.54) is approximate (0.77 % high at μ = 0; 6.5 % high at μ = −0.5, where σ = 0 is itself only
    assumed); exchange of stabilities is proved for Bénard (Ra > 0) and co-rotating Taylor only; Squire's theorem holds for
    parallel non-rotating, unstratified flows; the parallel-flow Orr–Sommerfeld analysis of a boundary layer is an approximation
    (N109); the double-diffusive margin σ = 0 is for fingers only (the diffusive onset is oscillatory — computed with our cubic,
    labelled "ours"); the Lorenz truncation is valid only slightly above onset; the Lorenz largest Lyapunov exponent ≈ 0.906 is
    a literature value used qualitatively (our single-run fit is printed with its window); TG/Rayleigh growth rates near a neutral
    curve converge slowly (critical layer) and are labelled "approximate".
13. **Verifier notes found while designing (must reach the verifier).** (i) The book's (11.45) amplitudes $\hat T,\hat s$ are
    *scaled by the gradients* (our D13 step 7: $\hat T\to-\hat T/(d\,d\bar T/dz)$); the identity $\hat T=\kappa_s\hat s/\kappa$
    holds for the scaled amplitudes, so `double_diffusive_sigma` must be built from the dimensional equations (D13) — its
    τ = 1, Rs = 0 limit reproduces `benard_free_free_sigma` with Ra_§11.4 = −Ra_§11.5 (checked: 11.0155 at K = π/√2, Ra = 2000,
    Pr = 1). (ii) The KH quadratic printed on p. 506 contains $(-\mathrm ikc+\mathrm ikU_1)^2=-k^2(U_1-c)^2$; a test that drops the
    $\mathrm i$ gets the sign of the shear term wrong — plant it. (iii) The Taylor narrow-gap Ta is $-4A\Omega_1d^4/\nu^2$ with A of
    (11.49) (D14 step 13): it is ≤ 0 on and beyond the Rayleigh line, so `taylor_number` must not be fed μ ≥ (R₁/R₂)² silently —
    return the value and a `rayleigh_stable` flag. (iv) The book's semicircle is unnumbered; (11.72) is the integral inequality —
    explainers' `ref:` badges say "unnumbered, p. 507" for the semicircle.

Order of parts: C (contract) · A (notebook storyboard) · B (explainers) · D (runtime) · E (prerequisite ledger) · F (derivations).

---
## Part C — functions the builders will call (the implementer's contract)

**Status column.** **A#** = planned in `analysis/ch11.md` §4 row # (signature kept or refined here). **§8.n** = added by the
curation's §8 item n. **NEW** = added by this design (in neither) — flagged as the phase asks. **CHG** = the analysis signature is
changed here (reason given). **reuse** = exists (ch01–ch10) and is only called. Docstrings cite § and Eq. (from the rendered page),
list symbols with units and assumptions, carry the validation label and **state which convention they use** (Γ sign, ψ sign,
σ or c, scales). Return types: `float`/`complex`/`ndarray` for single quantities, `dict` (keys listed) for several named results.

### C.0 Reused (existing; called, not changed)
| Callable | Used for | Where |
|---|---|---|
| `WV.interface_omega(k, rho1, rho2, g)` (ρ₁ upper) | R09 parity: KH with U₁ = U₂ = 0 is ch07 (7.95) | C02 check cell, E2 parity row |
| `WV.real_field(amp, phase)` | real part of a complex mode | C01 from-scratch |
| `ch07.interface_fields(x, z, t, a, k, rho1, rho2, …)` | the potentials $Ae^{\mp kz}$ shared with (11.15) | R08 |
| `ch05.sheet_rollup(N, gamma, amplitude, delta, t_eval, L)` | nonlinear roll-up (R11, A2) | C02 |
| `core.bernoulli.unsteady_bernoulli_pressure(dphi_dt, speed, z, rho, g, C)` | R06 reminder | C02 |
| `KIN.galilean_transform(u, U)` | N19 moving frame | C02 |
| `STRAT.brunt_vaisala_sq(rho0, drho_dz, drho_a_dz=0, g)` | R18, N² of a thermocline | C08 |
| `LAM.circular_couette(R, R1, R2, Omega1, Omega2, return_coeffs=True)` | R16 base state, A and B of (11.49) | C07, E5 |
| `LAM.channel_flow(y, h, U, dpdx, mu, G)` | Poiseuille / Couette base profiles (sketches) | C13 |
| `BL.blasius_profile(eta)` → (f, f′, f″); `BL.blasius_delta_star`; `BL.falkner_skan(m, …)`; `BL.profile_inflection(y, u)` | Blasius and FS base flows (U″ from the ODE, never a numerical second derivative) | C13 (N107, N108, R24) |
| `JET.free_jet_profile(eta)` → dict(f, fp, fpp) | Bickley sech² (R23) | C13 |
| `core.similarity.prandtl_number(nu, kappa)` | Pr | C03 |
| `ch01.seawater_density_linear(T, S, rho0, alpha_T, beta_S, T0, S0)` | linear EOS recalled (N53) | C06 |
| `ch10.dominant_frequency(t, signal)` | N112 spectrum peak | C15 |
| `FD.fd_weights`, `FD.fd_derivative` (ch10) | the FD cross-check of `cheb` (C03 primer) | C03 |
| `tools.convergence.observed_order` | spectral convergence shown as digits gained, not an order | C03 |

### C.1 `fluidpy/core/stability.py` (NEW module, `ST`; re-exported by `ch11`) — the eigen-solver toolkit (Ch. 11, 12, 13, 14)
Chebyshev–Gauss–Lobatto collocation (Trefethen ordering: $x_j=\cos(j\pi/N)$ mapped to the domain, so **x[0] is the right/top
end**), boundary rows replaced in A and zeroed in B, `scipy.linalg.eig(A, B)`, infinite eigenvalues dropped, spurious ones
removed by **N-convergence** (keep eigenvalues that move < `tol` when N → round(1.5 N)), never by a tight ∣c∣ cap alone (the TS
mode is weak). Never N > 150 in the notebook.

| # | Callable (signature) | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 1.1 | `cheb(N, domain=(-1.0, 1.0))` | Trefethen's matrix (off-diagonal $c_i(-1)^{i+j}/(c_j(x_i-x_j))$, negative-sum diagonal), scaled by 2/(b − a) | (D (N+1)×(N+1), x (N+1,)) — *expect* N = 4: x = [1, 0.7071, 0, −0.7071, −1]; max∣D x³ − 3x²∣ = 4.4e-16 | C03 primer, from-scratch, every solver | A#1 |
| 1.2 | `cheb_matrices(N, orders=(1, 2, 4), domain=(-1.0, 1.0))` | powers of D | dict(x, D1, D2, D4) | solvers | A#2 |
| 1.3 | `clenshaw_curtis_weights(N, domain=(-1.0, 1.0))` | $\int f\,dx\approx\sum w_jf(x_j)$ (exact for degree ≤ N) | w (N+1,) — *expect* Σw = b − a; ∫x² on [−1, 1] = 2/3 to 1e-15 | C09, C10, C14, P271 | A#26 (name kept from §5) |
| 1.4 | `apply_bc_rows(A, B, rows)` (rows = list of (row_index, row_vector)) | replace row r of A by the BC vector, zero row r of B | (A, B) copies | C03 from-scratch, solvers | A#3 |
| 1.5 | `generalized_eigs(A, B, sort="imag", c_max=None, return_vectors=False)` | `scipy.linalg.eig`, drop non-finite, optional ∣λ∣ < c_max, sort descending by imag ("imag") or real ("real") | λ (and V) | solvers | A#3 |
| 1.6 | `converged_eigs(eig_fn, N, factor=1.5, tol=1e-6)` | keep λ of `eig_fn(N)` within tol·max(1, ∣λ∣) of some λ of `eig_fn(round(factor·N))` | λ sorted | C03 figure (spurious modes), every spectrum | **NEW** |
| 1.7 | `normal_mode(x, y, t, uhat, k, m=0.0, sigma=None, c=None)` | (11.1): first form $\mathrm{Re}\{\hat u e^{\mathrm i(kx+my)+\sigma t}\}$ if `sigma` given, second $\mathrm{Re}\{\hat u e^{\mathrm i\lvert\mathbf K\rvert(\mathbf e_K\cdot\mathbf x-ct)}\}$ if `c` given (exactly one) | real ndarray | C01 | A#4 |
| 1.8 | `sigma_from_c(K, c)` · `c_from_sigma(K, sigma)` | σ = −i∣K∣c; c = iσ/∣K∣ | complex — *expect* K = 2, c = 1 + 0.5i → σ = 1 − 2i | C01, E1 | A#4 |
| 1.9 | `stability_class(sigma=None, c=None, tol=1e-10)` | sign of σ_r (or c_i) of one mode | "stable" / "neutral" / "unstable" | C01, E1 | A#4 |
| 1.10 | `stability_verdict(growth, tol=1e-10)` (growth = array over k) | "for every k": unstable if any > tol, neutral if max within tol, else stable (P255) | dict(verdict, k_index_max, max_growth) | C01, E1, E2 | **NEW** |
| 1.11 | `marginal_type(sigma, tol=1e-8)` | σ_r = 0 with σ_i = 0 → "stationary", else "oscillatory" | str | C01 (N06), C06 | A#4 |
| 1.12 | `rayleigh_eigs(k, U, Upp, domain=(-1.0, 1.0), N=120, bc="wall", map_scale=None, parity=None, return_vectors=False, filter=True)` | (11.81) with (11.82); `bc="decay"` uses the algebraic map $y=a\xi/\sqrt{1-\xi^2}$ (a = `map_scale`, default 3) with φ → 0; `parity="even"/"odd"` restricts to symmetric/antisymmetric φ (Ex. 11.12, (11.95)) | c sorted by c_i (vectors: dict(y, phi)) | C12, C13, E7 tables | A#31 |
| 1.13 | `taylor_goldstein_eigs(k, U, Upp, N2, domain=(0.0, 1.0), N=80, bc="wall", map_scale=3.0, return_vectors=False, filter=True)` | (11.61) multiplied by (U − c): $c^2M_2+cM_1+M_0=0$, $M_2=D^2-k^2$, $M_1=-2U(D^2-k^2)+U''$, $M_0=U^2(D^2-k^2)-UU''+N^2$; companion linearisation (P268) | c sorted by c_i | C08, C09, E6 tables | A#24 |
| 1.14 | `orr_sommerfeld_eigs(k, Re, U, Upp, domain=(-1.0, 1.0), N=100, bc="wall", y_max=None, map_scale=None, return_vectors=False, filter=True)` | (11.79) with (11.80); `bc="semi_infinite"` (Blasius, FS) maps $[0,y_{max}]$ algebraically, φ = φ′ = 0 at both ends | c sorted by c_i — *expect* Poiseuille (U = 1 − y², U″ = −2), Re = 10⁴, k = 1, N = 100: **c = 0.23752649 + 0.00373967i** (N = 60/80/100/120 agree to 1e-8) | C11, C13, C14, E8 tables | A#28 |
| 1.15 | `os_mode(k, Re, U, Upp, domain=(-1.0, 1.0), N=100, **kw)` | leading mode, normalised so max∣û∣ = 1 with û real and positive at that point | dict(c, y, phi, u_hat (= φ′), v_hat (= −ikφ)) | C14, A4, E8 | §8.8 |
| 1.16 | `max_growth(eig_fn, k_bounds, n_scan=24)` (eig_fn(k) → leading c) | scan then bounded `minimize_scalar` of −k c_i | dict(k, c, growth) — *expect* tanh Rayleigh: k = 0.4449, kc_i = 0.1897 (Michalke 1964: 0.4446, 0.1897) | C12, C13 | A#34 |
| 1.17 | `neutral_curve(eig_fn2, Re_values, k_bounds, n_k=40)` (eig_fn2(k, Re) → leading c) | for each Re: the k-interval where c_i > 0 (brentq on both edges) | dict(Re, k_lower, k_upper) | C13, F7 | A#34 |
| 1.18 | `critical_point(eig_fn2, k_bounds, Re_bounds, k_guess=None)` | brentq on c_i(Re) = 0 at fixed k inside bounded minimisation over k (P263) | dict(Re_c, k_c, c_r) — *expect* Poiseuille 5772.22, 1.02056, 0.26400 | C13 | A#34 |
| 1.19 | `howard_semicircle(Umin, Umax, n=200)` · `in_howard_semicircle(c, Umin, Umax, tol=1e-9)` | the unnumbered semicircle on p. 507 and (11.71) | (c_r, c_i) arrays · bool | C10, E7 | A#26 |
| 1.20 | `inflection_points(y, U=None, Upp=None)` (callables → brentq on sign changes; arrays → linear interpolation) | zeros of U″ where it changes sign | ndarray of y_I | C12, E7 | A#32 |
| 1.21 | `squire_transform(k, m, Re, uhat=None, what=None, phat=None)` | (11.78): $\bar k=\sqrt{k^2+m^2}$, $\bar c=c$, $\bar k\bar u=k\hat u+m\hat w$, $\bar p/\bar k=\hat p/k$, $\bar k\overline{\mathrm{Re}}=k\mathrm{Re}$ | dict(kbar, Rebar, ubar, pbar, angle_deg) — *expect* m/k = tan 30°: Rebar = 0.8660 Re | C11, E8 | A#27 |
| 1.22 | `disturbance_energy_budget(k, c, phi, y, Up, Re, weights=None)` (§11.8 convention u = φ′, v = −ikφ; wavelength-averaged; `weights` default Clenshaw–Curtis on the Chebyshev y) | (11.88) 2-D form: $E=\tfrac14\int(\lvert\hat u\rvert^2+\lvert\hat v\rvert^2)dy$, production $P=-\tfrac12\int\mathrm{Re}(\hat u\hat v^*)U'dy$, dissipation $\Lambda=\tfrac1{2\mathrm{Re}}\int(k^2\lvert\hat u\rvert^2+\lvert\hat u'\rvert^2+k^2\lvert\hat v\rvert^2+\lvert\hat v'\rvert^2)dy$ | dict(E, dEdt (= 2kc_iE), production, dissipation, residual, ratio P/Λ) — *expect* Poiseuille Re = 10⁴, k = 1 (mode normalised as returned by eig): P/Λ = **1.616**, residual dE/dt − (P − Λ) ≈ 1e-12; Re = 5000, k = 1: P/Λ = **0.765** (decays) | C14, A4, E8 | A#36 |

### C.2 `fluidpy/ch11_instability.py` (chapter module; re-exports C.1)
**§11.1–11.3 — normal modes and Kelvin–Helmholtz** (dimensional SI; layer 1 = upper, 2 = lower, as the book and ch07)
| # | Callable (signature) | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 2.1 | `potential_well_demo(shape, x0, v0=0.0, damping=0.3, t_end=20.0, n=400)` (shape ∈ "bowl", "cap", "plane", "dimple") | ball in V(x) (bowl x²/2, cap −x²/2, plane 0, dimple x²/2 − x⁴/4 with rim at ∣x∣ = 1), $\ddot x=-V'(x)-\gamma\dot x$ | dict(t, x, v, V (callable), escaped) — *expect* dimple: x0 = 0.3 stays (∣x∣ < 1), x0 = 1.2 escapes | C01 A1 | A#39 |
| 2.2 | `normal_mode_growth(system, k, **params)` (system "interface": rho1, rho2, g, surface_tension; "kh": U1, U2, rho1, rho2, g, surface_tension, h; "benard_free": Ra, Pr) | σ = −ikc₊ for interface/KH (c₊ = the root with c_i ≥ 0), larger root of D12 for Bénard (K = k) | complex σ — *expect* ("kh", k = 1, U1 = 6, U2 = 0, rho1 = 1, rho2 = 3, g = 10): σ = 1.3229 − 1.5i | C01, E1, F1 | §8.3 |
| 2.3 | `kh_phase_speed(k, U1, U2, rho1, rho2, g=G0, surface_tension=0.0, h=None)` | (11.18); with σ_s and/or lower depth h the Exercise 11.1 form (N120, D04) | (c_plus, c_minus) complex, c_plus.imag ≥ 0 — *expect* k = 1, U1 = 6, U2 = 0, ρ = 1, 3, g = 10: **1.5 ± 1.3229i**; U1 = 4: **2.4142, −0.4142** (neutral) | C02, E2 | A#5 |
| 2.4 | `kh_discriminant_terms(k, U1, U2, rho1, rho2, g=G0, surface_tension=0.0, h=None)` | the three terms under the root of (11.18)/(Ex. 11.1): gravity $\frac{\rho_2-\rho_1}{\rho_2+\rho_1}\frac gk$, tension $\frac{\sigma_sk}{\rho_1+\rho_2}$, shear $-\frac{\rho_1\rho_2(U_2-U_1)^2}{(\rho_1+\rho_2)^2}$ (coth kh factors when h is given) | dict(gravity, tension, shear, total, mean) — *expect* tiny example: 5, 0, −6.75, −1.75, 1.5 | C02, E2 terms | **NEW** |
| 2.5 | `kh_growth_rate(k, U1, U2, rho1, rho2, g=G0, surface_tension=0.0, h=None)` | k·max(c_i, 0) [1/s] | float / ndarray | C02, F1 | A#5 |
| 2.6 | `kh_critical_k(U1, U2, rho1, rho2, g=G0)` | $k_c=g(\rho_2^2-\rho_1^2)/(\rho_1\rho_2(U_2-U_1)^2)$ (no surface tension) | float [1/m] — *expect* air/water ΔU = 5 m/s: 327 m⁻¹ (λ_c = 1.92 cm) | C02 | A#5 |
| 2.7 | `kh_stability_boundary(k, rho1, rho2, g=G0, surface_tension=0.0, h=None)` | ΔU_min(k) where the discriminant is zero | float / ndarray [m/s] | C02, E2 | §8.2 |
| 2.8 | `kh_min_shear(rho1, rho2, g=G0, surface_tension=0.074)` | $\Delta U_{\min}^2=2\sqrt{g\Delta\rho\,\sigma_s}(\rho_1+\rho_2)/(\rho_1\rho_2)$ at $k_*=\sqrt{g\Delta\rho/\sigma_s}$ (requires σ_s > 0; σ_s = 0 raises with a message pointing to `kh_critical_k`) | dict(dU_min, k_star, wavelength) — *expect* air over water: **6.70 m/s at λ = 1.73 cm** | C02, E2 | §8.2, CHG (dict return; σ_s = 0 raises) |
| 2.9 | `kh_unstable_band(dU, rho1, rho2, g=G0, surface_tension=0.0)` | the k-interval with negative discriminant (quadratic in k with σ_s) | (k1, k2) (k2 = inf without σ_s) — *expect* ΔU = 8 m/s air/water: k ∈ [149.2, 887.4] m⁻¹, λ from 0.708 to 4.21 cm; ΔU = 5: empty (stable) | C02, E2 status | **NEW** |
| 2.10 | `vortex_sheet_c(U1, U2)` | (11.20) | (c_plus, c_minus) | C02 (R10) | A#5 |
| 2.11 | `rayleigh_taylor_cutoff(surface_tension, rho_heavy, rho_light, g=G0)` | Exercise 11.2 (ours): $\lambda_c=2\pi\sqrt{\sigma_s/((\rho_h-\rho_l)g)}$ | float [m] — *expect* water over air, σ_s = 0.074: **1.73 cm** | C02 (N121) | A#5 |
| 2.12 | `kh_amplitudes(k, c, U1, U2, zeta0)` | from (11.16): $A_-=-\mathrm i(U_1-c)\zeta_o$, $A_+=\mathrm i(U_2-c)\zeta_o$ | (A_minus, A_plus) complex | C02 (N14) | A#6 |
| 2.13 | `kh_fields(x, z, t, k, U1, U2, rho1, rho2, g=G0, zeta0=0.01, branch="+")` | (11.15) with (11.16) amplitudes, real parts | dict(zeta, phi1, phi2, u1, w1, u2, w2, c) | C02 figures, A2 | A#6 |
| 2.14 | `kh_residuals(k, U1, U2, rho1, rho2, g=G0)` | residuals of (11.3), (11.9), (11.13) for the computed mode at random points | dict(laplace, kinematic, dynamic) (all ≤ 1e-12) | C02 check | A#6 |
| 2.15 | `kh_sympy()` | D02 + D03 from (11.6)–(11.17) symbolically: linearised conditions, A±, the quadratic, (11.18), the discriminant identity $(\rho_1U_1+\rho_2U_2)^2-(\rho_1+\rho_2)(\rho_1U_1^2+\rho_2U_2^2)=-\rho_1\rho_2(U_1-U_2)^2$ | dict(kinematic_lin, dynamic_lin, quadratic, roots, identity (0), criterion) | C02 D03 check_src | A#6 |
| 2.16 | `kh_depth_tension_sympy()` | D04: lower mode cosh k(z + h), pressure jump, coth factors, h → ∞ limit | dict(dispersion, limit_residual (0), dU_min_sq) | C02 D04 | **NEW** |
| 2.17 | `kh_mixing_energy(U1, h, rho, profile="linear")` | N22: $E_i=\tfrac\rho2U_1^2h$, $E_f=\tfrac\rho2\int_{-h}^hU^2dz=\tfrac\rho3U_1^2h$, momentum ∫ρU dz both = ρU₁h | dict(E_i, E_f, ratio, M_i, M_f) — *expect* ratio 2/3 | C02 | A#7 |

**§11.4 — Bénard** (§11.4 scales; Ra > 0 heated from below; z ∈ [−½, ½]; bc tuple = (bottom, top))
| # | Callable (signature) | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 3.1 | `rayleigh_number(alpha, dT, d, kappa, nu, g=G0)` (**dT = T_bottom − T_top [K]**; Γ_T = dT/d) | (11.21) $\mathrm{Ra}=g\alpha\Gamma d^4/(\kappa\nu)$ = gα dT d³/(κν) | float — *expect* water d = 5 mm, ΔT = 2 K: **3679** (with g = 9.81; 3677 with G0) | C03, E3 | A#8, **CHG** (argument order: g last with a default; every caller uses keywords) |
| 3.2 | `gamma_conventions(dT, d)` | the three signs of one gradient: (11.21) $\Gamma_T=-d\bar T/dz=+dT/d$; ch01 Kundu $dT/dz=-dT/d$; meteorology $-dT/dz=+dT/d$ | dict(Gamma_11_21, dTdz_kundu, Gamma_met) [K/m] — *expect* 400, −400, 400 for d = 5 mm, ΔT = 2 K | C03 two-convention table, E3 badge | **NEW** |
| 3.3 | `benard_base_state(z, T0, dT, d, rho0=1000.0, alpha=2.1e-4, g=G0)` (z dimensional, centred) | (11.23)–(11.24) $\bar T=T_0-\tfrac12\Delta T-\Gamma z$ (`T0` = bottom temperature; see the note below the table) | dict(T, P) | C03 (R13, N28) | A#8 |
| 3.4 | `benard_scales(alpha, dT, d, kappa, nu, g=G0)` | N31 scaling | dict(Ra, Pr, Gamma, w_scale (κ/d), t_scale (d²/κ)) — *expect* water: w ~ 2.8e-5 m/s, t ~ 179 s | C03 | A#8 |
| 3.5 | `benard_perturbation_sympy()` | D05, D06, D07, D09: (11.25)–(11.29), (11.31)–(11.37), (11.39)–(11.41); planted ∇² for ∇_H² fails | dict(eq_11_27, eq_11_29_residual (0), planted_residual (≠ 0), eq_11_36, eq_11_37, eq_11_40_residual (0)) | C03, C04 check cells | A#9 |
| 3.6 | `exchange_of_stabilities_sympy()` | D08: the two energy relations and the imaginary-part identity $\sigma_i(J_1/\Pr+\mathrm{Ra}K^2I_1)=0$, built on polynomial trial functions that meet (11.38) | dict(relation_T, relation_W, imag_identity (0)) | C03 D08 check | A#15 |
| 3.7 | `benard_growth_rate(K, Ra, Pr, bc=("rigid", "rigid"), N=40, all=False)` | (11.36)–(11.37) with (11.38) as a 2(N+1) generalised eigenproblem in (W, T̂) | leading real σ (or all, sorted) — *expect* free–free check vs 3.10 to 1e-8; all ∣σ_i∣ < 1e-10 | C03, A3 | A#10 |
| 3.8 | `benard_marginal_Ra(K, bc=("rigid", "rigid"), mode="even", N=40)` (mode "even"/"odd"/"any") | (11.39) as a generalised eigenproblem for Ra | float — *expect* rigid–rigid K = 2: **2177.41**, K = 3.1163: **1707.76**, K = 5: **2439.32** | C04, F2 | A#11 |
| 3.9 | `benard_critical(bc=("rigid", "rigid"), mode="even", N=40)` · `benard_neutral_curve(Ks, bc=("rigid", "rigid"), mode="even", N=40)` | bounded minimisation over K | dict(Ra_c, K_c) — *expect* rigid–rigid **1707.762, 3.1163** (Chandrasekhar 1961: 1707.76, 3.117); free–free 657.511, 2.2214; rigid–free **1100.650, 2.6823**; odd **17610.39, 5.3647** · ndarray | C04, E3 tables | A#11 |
| 3.10 | `benard_free_free_sigma(K, Ra, Pr, n=1)` | D12: $(\sigma+a^2)(\sigma/\Pr+a^2)a^2=\mathrm{Ra}K^2$, $a^2=n^2\pi^2+K^2$ | (σ_plus, σ_minus) — *expect* K = π/√2, Ra = 2000, Pr = 1: **11.0155, −40.6243**; Ra = 657.511: σ₊ = 0 | C05, E1 | §8.1 |
| 3.11 | `benard_char_roots(Ra, K)` | (11.42) and $q_0$ | (q0 float, q complex, q_conj) — *expect* K = 2, Ra = 1000: q₀ = 3.4459, q = 3.8822 + 1.7705i | C04 (D10), E3 | A#12 |
| 3.12 | `benard_determinant(Ra, K, mode="even", printed=False)` | the 3 × 3 determinant of p. 489 (odd mode: sin/sinh analogue) | complex (purely imaginary when printed=False) — *expect* K = 3.1163, Ra = 1707.762: ≈ 0.0318i | C04, E3 | A#12 |
| 3.13 | `benard_marginal_Ra_det(K, mode="even", Ra_max=2e4, n_scan=400)` | first sign change of Im det above Ra = K⁴ (bracket starts at 1.0001 K⁴ + 1), then brentq | float — *expect* equals 3.8 to 1e-8 at K = 2, 3.1163, 5 | C04 from-scratch, E3 | A#12 |
| 3.14 | `benard_free_free_Ra(K, n=1)` · `benard_free_free_critical()` | (11.44); dRa/dK² = 0 | float · dict(Ra_c = 27π⁴/4 = 657.511, K_c = π/√2 = 2.2214, K_c2 = π²/2) — *expect* Ra(1) = 1284.23, Ra(2) = 667.01, Ra(3) = 746.53, Ra(4) = 1082.06 | C05, E1, E3 | A#13 |
| 3.15 | `benard_free_free_mode(z, n=1, printed=False)` | $\sin n\pi(z+\tfrac12)$ (printed: $\sin n\pi z$) | ndarray | C05 slip #2 | A#13 |
| 3.16 | `benard_free_free_sympy()` | D11: (11.43) from stress-free walls, residual of (11.40), dRa/dK², the root; the printed (slip #3) form and its (empty) root set | dict(bc_W4 (0), residual (0), dRa_dK2, root (π²/2), Ra_c, printed_dRa_dK2, printed_root ([])) | C05 D11 | A#13 |
| 3.17 | `benard_eigenfunction(K, bc=("rigid", "rigid"), mode="even", x=None, z=None, N=40)` | eigenvector of 3.8 (W, T̂), roll fields ψ (w = ∂ψ/∂x), w, T′ on a grid, normalised max∣w∣ = 1 | dict(x, z, psi, w, T, W_profile, T_profile, Ra) | C04 figures (Fig. 11.9 analogue), A3 | A#14 |
| 3.18 | `planform(x, y, K, kind="rolls")` (rolls, squares, hexagons) | f with $\nabla_H^2f=-K^2f$ (P264) | ndarray | C04 (N51) | A#14 |
| 3.19 | `benard_neutral_table(Ks=None)` | rigid–rigid, free–free, rigid–free and odd neutral curves on one K grid (0.5 … 10, 96 points) for E3/F2; writes `reference/ch11/benard_neutral_curves.csv` | dict(K, rigid, free, rigid_free, odd) | E3 table, F2 | **NEW** |

(3.3 note — the book's Fig. 11.8: bottom temperature $T_0$, top $T_0-\Delta T$, $\bar T=T_0-\Gamma(z+d/2)$, which is (11.24)
$\bar T=T_0-\tfrac12\Delta T-\Gamma z$. The docstring states it exactly so; `T0` = bottom temperature.)

**§11.5 — double diffusion** (SI gradients, z up; §11.5 sign of Ra; S in g/kg with β per (g/kg))
| # | Callable (signature) | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 4.1 | `linear_eos(T, S, rho0=1027.0, alpha=2e-4, beta_S=7.6e-4, T0=10.0, S0=35.0)` | $\tilde\rho=\rho_0[1-\alpha(\tilde T-T_0)+\beta(\tilde s-s_0)]$ (p. 492) | ρ [kg/m³] | C06 (N53) | A#16 |
| 4.2 | `thermal_rayleigh_signed(dTdz, d, alpha, nu, kappa, g=G0)` · `salinity_rayleigh(dSdz, d, beta_S, nu, kappa_s, g=G0)` · `salinity_rayleigh_prime(dSdz, d, beta_S, nu, kappa, g=G0)` | §11.5: $\mathrm{Ra}=g\alpha d^4(d\bar T/dz)/(\nu\kappa)$ (< 0 heated from below), $\mathrm{Rs}=g\beta d^4(dS/dz)/(\nu\kappa_s)$, $\mathrm{Rs}'=g\beta d^4(dS/dz)/(\nu\kappa)$ | float — *expect* our thermocline (dT/dz = 0.01 K/m, dS/dz = 0.002 (g/kg)/m), d = 5 cm: Ra = 875.9, Rs = 62 130 | C06, E4 | A#16 |
| 4.3 | `double_diffusive_margin(Ra, Rs)` | Rs − Ra − 27π⁴/4 (> 0: finger-unstable) | float | C06 | A#16 |
| 4.4 | `salt_finger_unstable(dTdz, dSdz, d, alpha=2e-4, beta_S=7.6e-4, nu=1e-6, kappa=1.4e-7, kappa_s=1.5e-9, g=G0)` | (11.46) | dict(unstable, lhs, margin, density_stable, R_rho, Ra, Rs) — *expect* d = 5 cm: lhs = 61 254, margin = 60 597, density_stable True, R_ρ = αT_z/(βS_z) = 1.316; thinnest finger-unstable layer **d = 1.61 cm** | C06, E4 | A#16 |
| 4.5 | `double_diffusive_sigma(K2, Ra, Rs, Pr, tau, n=1)` (Ra §11.5 sign, Rs with κ_s, τ = κ_s/κ) | **ours** (D13 note): free–free cubic $(\sigma/\Pr+a^2)a^2(\sigma+a^2)(\sigma+\tau a^2)=K^2[-\mathrm{Ra}(\sigma+\tau a^2)+\tau\mathrm{Rs}(\sigma+a^2)]$ | 3 complex roots sorted by real part (desc) — *expect* K² = π²/2, Pr = 7, τ = 0.0107: (Ra, Rs) = (1000, 2000): real root **+0.0308** (fingers); (−2×10⁴, −1.9×10⁶): **9.65 ± 70.39i** (diffusive, oscillatory, statically stable) | C06, E4 | A#17, §8.4 |
| 4.6 | `salt_finger_regime(dTdz, dSdz, d, alpha=2e-4, beta_S=7.6e-4, nu=1e-6, kappa=1.4e-7, kappa_s=1.5e-9, g=G0, Pr=7.0, K2=None)` | regime label from the static stability and the cubic (K2 default π²/2) | dict(regime ∈ "stable", "fingers", "diffusive", "overturning"; R_rho, margin, sigma_max (complex), density_stable, text) | C06, E4 status | §8.4, **CHG** (argument order = 4.4) |

**§11.6 — Taylor–Couette** (narrow gap: lengths by d, x ∈ [0, 1], k ≡ kd, σ ≡ σd²/ν; μ ≡ Ω₂/Ω₁)
| # | Callable (signature) | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 5.1 | `ring_interchange_energy(Gamma1, Gamma2, r1, r2)` | p. 496: $E=\Gamma^2/(8\pi^2r^2)$, ΔE | dict(E_i, E_f, dE) — *expect* Γ₁ = 4, Γ₂ = 2, r₁ = 1, r₂ = 2: E_i = 0.21531, E_f = 0.10132, **ΔE = −0.11399** (released) | C07, E5 inspector | A#18 |
| 5.2 | `rayleigh_circulation_criterion(r, U_theta)` | sign of dΓ²/dr on samples (Γ = 2πrU_θ) | dict(unstable, where (r values with dΓ²/dr < 0), dGamma2_dr) | C07 (R14) | A#18 |
| 5.3 | `couette_rayleigh_line(R1, R2)` | μ_R = (R₁/R₂)² | float — *expect* R₂/R₁ = 1.05: 0.9070 | C07, E5 | A#18 |
| 5.4 | `taylor_number(Omega1, Omega2, R1, R2, nu)` | (11.52) $\mathrm{Ta}=4\frac{\Omega_1R_1^2-\Omega_2R_2^2}{R_2^2-R_1^2}\frac{\Omega_1d^4}{\nu^2}$ | dict(Ta, rayleigh_stable (μ ≥ (R₁/R₂)²), mu) — *expect* R₁ = 0.1 m, R₂ = 0.105 m, Ω₁ = 0.5 rad/s, Ω₂ = 0, ν = 1e-6: **Ta = 6098** | C07, E5 | A#19, **CHG** (dict with the Rayleigh flag, verifier note iii) |
| 5.5 | `taylor_number_narrow_inner(Omega1, R1, d, nu)` | $2(\Omega_1R_1d/\nu)^2(d/R_1)$ | float — *expect* 6250 for the case above (ratio 0.976 = 1 − O(d/R₁)) | C07 (D15) | A#19 |
| 5.6 | `taylor_perturbation_sympy()` | D14 steps 1–13 (linearisation, the factor 2, $dU_\varphi/dR+U_\varphi/R=2A$, elimination of û_z, p̂, narrow-gap limit, Ta) + slip #12 units check | dict(term_2A (0 residual), operator_identity (0), narrow_gap_eqs, Ta_coefficient (−4AΩ₁d⁴/ν²), printed_continuity_units_ok (False)) | C07 D14 check_src | A#20 |
| 5.7 | `taylor_marginal_Ta(k, mu, N=40)` | (11.51) at σ = 0 with (11.53) | float — *expect* μ = 0: Ta(3.1266) = 3389.90 | C07 | A#21 |
| 5.8 | `taylor_critical(mu, N=40)` · `taylor_critical_approx(mu)` · `taylor_critical_table(mus=None)` | minimum over k; (11.54); a cached table (`reference/ch11/taylor_critical.csv`, μ from −0.5 to 1 in 0.05, ours) | dict(Ta_c, k_c) · float · dict(mu, Ta_c, k_c, approx) — *expect* μ = 1: **1707.76 at 3.1163**; 0.9: 1797.61 (approx 1797.9); 0.5: **2275.09 at 3.1175** (2277.3, +0.10 %); 0: **3389.90 at 3.1266** (3416.0, +0.77 %); −0.25: 4461.3 (+2.1 %); −0.5: 6413.7 at 3.1985 (6832, +6.5 %) | C07, F4, E5 | A#21, §8.5 |
| 5.9 | `taylor_growth_rate(k, Ta, mu, N=40)` | leading σ of (11.51) | complex σ (real for μ ≥ 0 — test) | C07 (N67) | A#21 |
| 5.10 | `taylor_eigenfunction(k, mu, x=None, z=None, N=40)` | marginal eigenvector → meridional stream function of the vortices (u_R = ∂ψ/∂z convention stated) | dict(x, z, u_R, u_phi, psi) | C07 figure, E5 table | §8.5 |
| 5.11 | `taylor_galerkin_Ta(k, mu, n_modes=4, printed=False)` | Exercise 11.9 route: $\hat u_\varphi=\sum C_m\sin m\pi x$ (11.94), û_R from the clamped 4th-order problem, Galerkin projection | float — agrees with 5.7 to 0.1 % at 4 modes | C07 (N124) | A#22 |
| 5.12 | `taylor_stability_boundary(R2_over_R1=1.05, mus=None)` | narrow-gap Ta_c(μ) converted to (Ω₂R₂²/ν, Ω₁R₂²/ν) for our radius ratio, with the Rayleigh line | dict(x_outer, y_inner, rayleigh_x, rayleigh_y) | C07 figure (Fig. 11.17 analogue) | A#21, **CHG** (our ratio, not the book's) |

**§11.7 — stratified parallel flows** (§11.7: z up, u = ∂ψ/∂z, w = −∂ψ/∂x)
| # | Callable (signature) | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 6.1 | `stratified_shear_sympy(printed=False)` | D16, D17, D18 (steps 3–8), D19 (steps 2–5): (11.55)–(11.61); (11.61) ⇔ (11.64) under (11.63); (11.61) ⇔ the divergence form under (11.68) with (U − c)²; printed ∂p/∂x fails | dict(tg_residual (0), self_adjoint_residual (0), divergence_residual (0), printed_slip11_residual (≠ 0), printed_slip4_reaches_tg (False)) | C08–C10 check cells | A#23 |
| 6.2 | `gradient_richardson(z, U=None, N2=None, dUdz=None)` (callables or arrays) | (11.66) $\mathrm{Ri}=N^2/(dU/dz)^2$ | ndarray | C09, E6 | A#25 |
| 6.3 | `miles_howard_stable(z, U=None, N2=None, dUdz=None)` | (11.67) | dict(guaranteed_stable, Ri_min, z_min, text) | C09, E6 status | A#25 |
| 6.4 | `richardson_profiles(kind="tanh", J=0.1, R=1.0)` | U = tanh z, N² = J sech²(Rz) (R = shear/density thickness ratio; R = 1: Ri(z) = J cosh²z, Ri_min = J) | dict(U, Up, Upp, N2, Ri) callables | C09, E6 | §8.6 |
| 6.5 | `tg_growth(k, J, R=1.0, N=100, map_scale=3.0)` · `tg_growth_map(ks, Js, R=1.0, N=100, cache=True)` | leading kc_i of 1.13 on the infinite domain (decay map) for 6.4; map cached in `reference/ch11/tg_growth_map.csv` (ours) | float · dict(k, J, kci (2-D)) — *expect* J = 0: max kc_i = **0.1897 at k = 0.445**; J = 0.1, k = 0.4: ≈ 0.124 (approximate); J ≥ ¼: 0 everywhere (Miles–Howard) | C09, F5, E6 table | §8.6 |
| 6.6 | `richardson_identity_check(k, c, psi, z, U, Up, Upp, N2, w=None)` | (11.65) and its imaginary part for a computed TG mode (φ = ψ̂/(U − c)^{1/2}) | dict(lhs, rhs, residual, imag_lhs, imag_rhs) (residual ≤ 1e-8) | C09 | A#26 |
| 6.7 | `howard_identity_check(k, c, psi, z, U, N2, w=None)` | (11.69), (11.70) with F = ψ̂/(U − c), Q = ∣F′∣² + k²∣F∣² | dict(res_69, res_70, cr_mean (∫UQ/∫Q)) | C10 | A#26 |
| 6.8 | `stratified_energy_budget(k, c, psi, z, U, N2)` | Exercise 11.10 (N125, named): kinetic + available potential energy rates | dict(dKdt, dPdt, production, residual) | C14 (N125) | A#36 |

**§11.8–11.10 — Squire, Orr–Sommerfeld, Rayleigh, viscous results** (§11.8: u = ∂ψ/∂y, v = −∂ψ/∂x; û = φ′, v̂ = −ikφ)
| # | Callable (signature) | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 7.1 | `os_derivation_sympy(printed_v_sign=False)` | D20, D21: (11.74)–(11.77) and (11.79) from (11.77) with m = ŵ = 0; planted v̂ = +ikφ fails | dict(residual (0), planted_residual (≠ 0), steps) | C11 D21 check | A#30 |
| 7.2 | `os_3d_eigs(k, m, Re, U, Upp, N=60)` | (11.77) in primitive variables (only for the Squire check) | c sorted | C11 (N88) | A#29 |
| 7.3 | `inviscid_profile(name, **p)` (names "blasius_like", "couette", "poiseuille", "wall_vorticity_max", "jet", "shear_layer" (tanh), "bickley", "sin" (p: b)) | analytic stand-ins of Fig. 11.21 and the §11.9 examples on their domains | dict(U, Up, Upp (callables), domain, bc, label) | C12, E7 | §8.7 |
| 7.4 | `rayleigh_criterion(y, U=None, Upp=None)` · `fjortoft_criterion(y, U=None, Upp=None)` | (11.84) necessary condition; (11.86) $(U-U_I)U''<0$ somewhere | dict(has_inflection, y_I) · dict(satisfied, U_I, y_I, min_product) | C12, E7 badges | A#32 |
| 7.5 | `rayleigh_identity_check(k, c, phi, y, U, Upp, w=None)` | (11.83) real and imaginary parts for a computed mode | dict(res_real, res_imag) | C12 | A#32 |
| 7.6 | `rayleigh_spectrum_table(names=None, ks=None, N=120, cache=True)` | leading c for each profile and k (`reference/ch11/rayleigh_spectra.json`, ours) | dict(name → dict(k, c_r, c_i)) | E7, F6 | §8.7 |
| 7.7 | `sin_profile_max_growth(b, ks=None, N=80)` | N95: U = sin y on ∣y∣ ≤ b, walls; max kc_i over k | float (→ 0 as 2b → π⁺) | C12 | **NEW** |
| 7.8 | `critical_layer(y, U, c)` · `cats_eye_streamfunction(x, y, y_c=0.0, A=0.1, phi_c=1.0, k=1.0, Uy_c=1.0, U=None, c=None, exact=False)` · `cats_eye_width(A, phi_c, Uy_c)` | N96, N97, (11.87) and its expansion; eye width $2\sqrt{A\phi_c/U'_c}$ (ours) | float · ndarray · float — *expect* A = 0.1, φ_c = 1, U′_c = 1: width 0.632 | C12, E7 | A#33, **NEW** (`cats_eye_width`) |
| 7.9 | `piecewise_shear_layer_c(kh, dU=1.0)` · `piecewise_neutral_kh()` | Exercise 11.11 formula (N126); neutral root of $(kh-1)^2=e^{-2kh}$ | complex · float — *expect* **1.2785** | C12 | A#38 |
| 7.10 | `tanh_shear_layer_neutral_curve(Re_values=None, N=100, cache=True)` · `tanh_max_growth()` | N99: OS with decay map for U = tanh y | dict(Re, k_upper) · dict(k, kci, c) | C13, E8 table | A#35 |
| 7.11 | `bickley_critical(parity="sinuous", cache=True)` | R23 (Tatsumi & Kakutani 1958: Re_c ≈ 4 at k ≈ 0.2, approximate) | dict(Re_c, k_c) | C13 | A#35 |
| 7.12 | `poiseuille_critical(N=100, cache=True)` · `poiseuille_neutral_curve(Re_values=None, N=80, cache=True)` · `poiseuille_spectrum(k, Re, N=100)` | C13 | dict(Re_c, k_c, c_r) — *expect* **5772.22, 1.02056, 0.26400** · dict(Re, k_lower, k_upper) — *expect* Re = 10⁴: unstable k from ≈ 0.80 to ≈ 1.07 · c array | C13, F7, E8 | A#35 |
| 7.13 | `couette_max_growth(k, Re, N=80)` | N100 | float (< 0) — *expect* k = 1, Re = 10³: c_i = −0.1192; Re = 10⁴: −0.0521 | C13 | A#35 |
| 7.14 | `blasius_base(y_max=20.0)` · `blasius_critical(N=100, y_max=20.0, cache=True)` · `blasius_neutral_curve(Re_values=None, in_frequency=False, cache=True)` | δ*-scaled U = f′(1.7208 y), U″ = (1.7208)²(−½ff″); N107, N108 | dict(U, Upp) · dict(Re_c, k_c, omega_c) (Thomas via Gallagher et al. 2016: 519.2, 0.303, 0.120; Jordinson 1970: 520) · dict(Re, k or F lower/upper) | C13 | A#35 |
| 7.15 | `falkner_skan_neutral_curve(m, Re_values=None, cache=True)` | R24 (FS base flows from `BL.falkner_skan`) | dict(Re, k_lower, k_upper) | C13 | A#35 |
| 7.16 | `table_11_1(cache=True)` | N102 recomputed: jet (Bickley), shear layer (tanh), Blasius, plane Poiseuille, pipe (stated, not computed), plane Couette | list of dict(flow, U, Re_c_ours, k_c_ours, benchmark, source, length_scale, remark) | C13 | A#35 |
| 7.17 | `ts_wave_fields(x, y, t, k, c, phi, phi_y, amp=0.05, U=None)` (§11.8 convention) | real perturbation ψ, u, v (+ U(y) + u if `U` is given) | dict(psi, u, v, u_total) | C14, A4, E8 table | §8.8 |
| 7.18 | `energy_equation_sympy()` | D23: (11.88) and its 2-D form from (11.96), with periodic averaging | dict(residual (0), terms) | C14 D23 check | A#36 |
| 7.19 | `tollmien_profile(eta, eta1=0.2, printed=False)` | N105: $a\eta$, $1-b(1-\eta)^2$, 1 with $a=2/(1+\eta_1)$, $b=1/(1-\eta_1^2)$ from value + slope continuity (**our breakpoint**; the book's coefficients stay private) | ndarray — *expect* η₁ = 0.2: a = 1.6667, b = 1.0417; printed fails continuity | C13 (N105) | A#37, **CHG** (coefficients from continuity) |
| 7.20 | `neutral_curve_tables(out_dir="reference/ch11")` (script-level) | writes the Poiseuille, Blasius, tanh, Bickley neutral curves, the Poiseuille spectrum at k = 1 for 24 Re, mode samples and budgets on a (Re, k) grid for E8 (≤ 4 s.f.) | dict of written paths | scripts, E8 | §8.8 |

**§11.14 — chaos**
| # | Callable (signature) | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 8.1 | `pendulum_rhs(t, s, g_over_l=1.0, damping=0.0)` · `pendulum_energy(s, g_over_l=1.0)` · `phase_portrait(rhs, starts, t_end=20.0, n=400, **kw)` | (11.89) | list · float · list of dict(t, X, Y) — energy conserved to 1e-9 (damping 0) | C15 (N113) | A#39 |
| 8.2 | `hopf_normal_form(t, s, mu, omega=1.0)` · `limit_cycle_amplitude(mu)` | N114 (ours): $\dot z=(\mu+\mathrm i\omega)z-\lvert z\rvert^2z$ | list · √μ (0 for μ ≤ 0) | C15 | A#39 |
| 8.3 | `lorenz_rhs(t, s, Pr=10.0, r=28.0, b=8/3)` · `lorenz_integrate(s0=(1.0, 1.0, 1.0), t_end=40.0, Pr=10.0, r=28.0, b=8/3, n=4001, rtol=1e-10, atol=1e-12, method="DOP853")` · `lorenz_rk4(s0, dt, n_steps, Pr=10.0, r=28.0, b=8/3)` | (11.91); fixed-step RK4 twin of the JS | list · dict(t, X, Y, Z) · ndarray (n+1, 3) — *expect* rhs(1,1,1) = (0, 26, −1.667); one RK4 step dt = 0.01 from (1,1,1): (1.012567, 1.259918, 0.984891) | C15, E9 parity | A#40, **NEW** (`lorenz_rk4`) |
| 8.4 | `lorenz_b(k)` · `lorenz_r(Ra, k)` | $b=4\pi^2/(\pi^2+k^2)$; $r=\mathrm{Ra}\,k^2/(\pi^2+k^2)^3$ | float — *expect* b(π/√2) = 8/3 | C15 (D24) | A#40, **NEW** (`lorenz_r`) |
| 8.5 | `lorenz_fixed_points(r, b=8/3)` · `lorenz_jacobian(s, Pr=10.0, r=28.0, b=8/3)` · `lorenz_eigs(r, Pr=10.0, b=8/3, which="C")` · `lorenz_hopf_r(Pr=10.0, b=8/3)` | N116, D25 | list of 3-tuples · 3×3 · eigenvalues sorted by real part · float — *expect* r = 28: C± = (±8.4853, ±8.4853, 27), eig(C) = −13.8546, 0.0940 ± 10.1945i; origin 11.8277, −2.6667, −22.8277; **r_H = 24.7368** | C15, E9 | A#40, **NEW** (`lorenz_eigs`) |
| 8.6 | `lorenz_separation(s0=(1.0, 1.0, 1.0), delta0=1e-8, t_end=40.0, **kw)` · `lorenz_predictability_time(delta0=1e-8, threshold=1.0, s0=(1.0, 1.0, 1.0))` · `lorenz_largest_lyapunov(t_end=200.0, renorm_dt=1.0, cache=True)` (optional) | N112 sensitivity | dict(t, sep, slope, window) · float — *expect* δ₀ = 1e-8 from (1,1,1) + 1e-8 eₓ: separation reaches 1 at **t ≈ 29.1**, 10 at t ≈ 33.5 · float (≈ 0.9, literature 0.906, qualitative) | C15, E9 | A#40, §8.9 |
| 8.7 | `lorenz_fields(x, z, X, Y, Z, k=np.pi/np.sqrt(2))` | (11.90) with the §11.14 convention u = −∂ψ/∂z, w = ∂ψ/∂x | dict(psi, T, u, w) | C15, A5, E9 | A#40 |
| 8.8 | `lorenz_r_sweep(r_values, t_end=30.0, n=1500, cache=True)` | X–Z samples for F8 | dict(r, X (n_r, n), Z) | F8 | §8.9 |
| 8.9 | `lorenz_sympy()` | D24: projection of the truncated equations, the scalings and (11.91) | dict(dA, dB, dC, scaled_residuals (0), b_at_kc (8/3)) | C15 D24 check_src | **NEW** |
| 8.10 | `logistic_map(A, x0, n)` · `cobweb(A, x0, n)` · `bifurcation_diagram(A_values, n_transient=500, n_keep=100, x0=0.5)` · `period_doubling_points(n_max=6)` · `feigenbaum_estimate(n_max=6)` | N117, N129 | ndarray · dict(xs, ys) · dict(A, x) · dict(A_n (3, 1 + √6 = 3.449490, 3.544090, 3.564407 …), S_n superstable (2, 3.236068, 3.498562, 3.554641, 3.566667, 3.569244, 3.569795)) · list — *expect* δ estimates 4.709, 4.681, 4.663, 4.668, 4.669 → 4.6692 | C15, B1 | A#41, **NEW** (`cobweb`) |

**Book slips, engines**
| # | Callable | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 9.1 | `book_slips()` | the slips table of convention 5 (key S1…S12, printed, correct, where, evaluator) | list of dicts | callouts, tests | A#42, §8.10 |
| 9.2 | `derive_all()` | runs every sympy engine (`kh_sympy`, `kh_depth_tension_sympy`, `benard_perturbation_sympy`, `exchange_of_stabilities_sympy`, `benard_free_free_sympy`, `taylor_perturbation_sympy`, `stratified_shear_sympy`, `os_derivation_sympy`, `energy_equation_sympy`, `lorenz_sympy`) | dict of residual summaries | tests | **NEW** |

### C.3 `scripts/ch11_*.py` (runnable demos, each `--no-show`; heavy runs cache to `outputs/ch11/` with a parameter hash)
`ch11_drawings.py` (helpers, no physics: `kh_sketch(ax)`, `benard_sketch(ax)` (Fig. 11.8 analogue), `taylor_sketch(ax)` (Fig.
11.16 analogue), `cv_sketch(ax)` (Fig. 11.25 analogue), `semicircle_sketch(ax)`, `six_profiles(ax_grid)` (Fig. 11.21 analogue from
`inviscid_profile`)) · `ch11_neutral_curves.py` (Poiseuille, Blasius, tanh, Bickley, FS neutral curves and critical points; writes
`reference/ch11/os_neutral_*.csv` and the E8 tables) · `ch11_tables.py` (Bénard neutral table, Taylor critical table, TG growth map,
Rayleigh spectra, Lorenz sweep; writes `reference/ch11/explainer_tables.json` with ≤ 4 significant figures, labelled "ours") ·
`ch11_lorenz.py` (attractor, separation, optional Lyapunov) · `ch11_kh_rollup.py` (A2 frames from ch05 `sheet_rollup`).

### C.4 Not in `analysis/ch11.md` §4 (flag list for the implementer)
`ST.converged_eigs`, `ST.stability_verdict`, `ch11.kh_discriminant_terms`, `kh_unstable_band`, `kh_depth_tension_sympy`,
`gamma_conventions`, `benard_neutral_table`, `salinity_rayleigh_prime`, `sin_profile_max_growth`, `cats_eye_width`, `lorenz_rk4`,
`lorenz_r`, `lorenz_eigs`, `lorenz_sympy`, `cobweb`, `derive_all` (NEW); `rayleigh_number` argument order, `kh_min_shear` dict
return, `salt_finger_regime` argument order, `taylor_number` dict return, `taylor_stability_boundary` with our ratio,
`tollmien_profile` coefficients from continuity (CHG); curation §8 items 1–11 are all included above (`benard_free_free_sigma`,
`kh_min_shear`, `kh_stability_boundary`, `normal_mode_growth`, `salt_finger_regime`, `taylor_critical_table`,
`taylor_eigenfunction`, `richardson_profiles`, `tg_growth_map`, `inviscid_profile`, `rayleigh_spectrum_table`, `os_mode`,
`ts_wave_fields`, `neutral_curve_tables`, `lorenz_r_sweep`, `lorenz_predictability_time`, `book_slips`, the JS tables). Nothing in
Parts A, B or F calls an unlisted name.

### C.5 Parity conventions (the explainers reproduce these exactly; the implementer fixes them in the docstrings)
| # | Callable / convention | Exact definition | Used by | Status |
|---|---|---|---|---|
| 5.1 | complex roots in JS | `[re, im]` pairs; the principal square root of a negative real is $+\mathrm i\sqrt{\lvert x\rvert}$ (numpy `np.lib.scimath.sqrt`); c₊ is the root with c_i ≥ 0 | E1, E2, E3 | **NEW** |
| 5.2 | `normal_mode_growth` branches | σ = −ik c₊; for a neutral pair (two real c) σ = −ik c₊ with c₊ = mean + √disc | E1 | **NEW** |
| 5.3 | `benard_determinant` | z ∈ [−½, ½]; rows = BCs at z = +½; q = √(K²[1 + ½s(1 + i√3)]) principal branch (Re q > 0); columns (cos q₀z, cosh qz, cosh q*z); returns the complex det (JS uses its imaginary part) | E3 live | **NEW** |
| 5.4 | Bénard neutral table | K grid `np.linspace(0.5, 10, 96)`; columns rigid, free, rigid_free, odd (Chebyshev N = 40) | E3, F2 | **NEW** |
| 5.5 | Taylor table | μ grid −0.5 … 1.0 step 0.05; Ta_c, k_c (N = 40); eigenfunction ψ on 41 × 41 (x ∈ [0, 1], z over one wavelength 2π/k_c) for μ = 0, 0.5, 1 | E5, F4 | **NEW** |
| 5.6 | TG growth map | U = tanh z, N² = J sech²z (R = 1); k = 0.05 … 1.0 step 0.05, J = 0 … 0.3 step 0.01; kc_i (0 where no converged unstable mode); N = 100, map scale 3; values ≤ 4 s.f. | E6, F5 | **NEW** |
| 5.7 | Rayleigh spectra | profiles of 7.3 with k = 0.1 … 2.0 step 0.1: leading c (c_r, c_i), N = 120 | E7, F6 | **NEW** |
| 5.8 | OS tables | Poiseuille (U = 1 − y², [−1, 1]), Blasius (δ*), tanh, Bickley: neutral curves on 40 log-spaced Re; leading c and budget (P, Λ, E) on a 24 Re × 25 k grid per flow; mode φ, û, v̂ on 41 y points for the presets | E8, F7 | **NEW** |
| 5.9 | Lorenz | RK4 with dt = 0.005 in the explainer; parity rows compare 200 steps with `lorenz_rk4` (same dt) to 1e-10, never the adaptive solver beyond t ≈ 5 | E9 | **NEW** |
| 5.10 | `reference/ch11/` | `SOURCES.md` (benchmarks of analysis §8 with URLs), `benchmarks.json` (Ra_c rigid 1707.76/3.117, Poiseuille 5772.22/1.02056/0.264 and c(10⁴, 1), Blasius 519.2/0.303/0.120, Bickley 4.0/0.2, tanh 0.4446/0.1897, Lorenz r_H, Feigenbaum δ, logistic A_n), our tables (above) — public, no book table | all | A (analysis §8) |

---

## Part A — notebook storyboard (`notebooks/build_ch11.py` → `notebooks/ch11_instability.ipynb`)

**One line per book section** (cell numbers are estimates for the builder's budget; ≈ 620 cells in all):
- §11.1 → N01 N02 N03 (open the C01 story) — cells ≈ 8–15
- §11.2 → C01 (N04 N05 N06 · D01 · P255 · A1 · fig · **E1**) — cells ≈ 16–45
- §11.3 → R01–R11, C02 (N07–N19 N22 N120 N121 named N20 N21 N23 · D02 D03 D04 · P256 · figs · F1 · A2 · **E2**), S01 — cells ≈ 46–120
- §11.4 → R12 R13, C03 (N24–N40 N122 · D05 D06 D07 D08 · P257 P258 P259 P260 · figs), C04 (N41–N46 N50 N51 N52 N123 · D09 D10 ·
  P261 P262 P263 P264 · figs · F2 · live · A3 · **E3**), C05 (N47 N48 N49 · D11 D12 · P265 · fig) — cells ≈ 121–260
- §11.5 → C06 (N53–N57 · D13 · P266 · fig · F3 · **E4**) — cells ≈ 261–290
- §11.6 → R14 R15 R16, C07 (N58–N69 N124 · D14 D15 · P267 · figs · F4 · **E5**) — cells ≈ 291–340
- §11.7 → R17 R18 R19, C08 (N70–N74 · D16 D17 · P268 P269 P270 · fig), R20, C09 (N75–N78 · D18 · P271 · fig · F5 · **E6**), C10
  (N79–N82 · D19 · P272 · fig · F6) — cells ≈ 341–430
- §11.8 → R21 R22, C11 (N83–N90 · D20 D21 · fig) — cells ≈ 431–460
- §11.9 → C12 (N91–N97 N126 N127 · D22 · figs · **E7**) — cells ≈ 461–495
- §11.10 → R23 R24, C13 (N98–N103 N107 N108 · figs · F7 · table), C14 (N104 N125 N128 · D23 · P273 · fig · A4 · **E8**) — cells ≈ 496–560
- §11.11 → **C13 (continued)**: N105 N106 N109 (+ N107 N108 pointer back) — cells ≈ 561–568
- §11.12 → **C14 (continued)**: N110 — cells ≈ 569–571
- §11.13 → **C14 (continued)**: N111 — cells ≈ 572–575
- §11.14 → C15 (N112–N119 N129 · D24 D25 · P274–P279 · figs · A5 · F8 · live · **E9**), S02; summary — cells ≈ 576–620

Every CORE block below follows the order: problem in plain words → idea → primers → maths (notes and derivations, Part F) →
tiny example → code (fluidpy) + "What does the code above do?" → from-scratch check → visual(s) → notes and "What would
change if…". Code drafts are intent + exact calls; the builder writes every line commented (novice grade, units, the
equation written out next to its number). *expect* gives the numbers the executed cell must print. *see / read / change*
are the three figure notes. Note ids open every note text in bold (**N09 [B]**) so the lesson-reviewer can find them.
Slip callouts are `> ⚠️ **slip #k — the book prints** … **; the correct form is** …` boxes.

### A.0 Front matter
1. `nb.title(big_idea=…, roadmap=[…15…], prerequisites=[…])`. **Big idea (draft):** "Every flow we solved so far was a
   *possible* flow — it satisfied the equations. Nature only shows us the possible flows that *survive a nudge*. A pencil can
   balance on its tip in a textbook, never on a desk. This chapter asks of each steady flow the same question: put a tiny
   wave-shaped disturbance on it, does the wave die or grow? Because the disturbance is tiny, the equations become linear, a
   wave of one wavelength evolves on its own, and the question becomes an eigenvalue problem: a growth rate σ (or a complex
   wave speed c) for every wavenumber. We meet the classic mechanisms one at a time — shear between two streams
   (Kelvin–Helmholtz billows, the start of wind waves), heating from below (Bénard convection, the atmosphere's shallow
   convection), two diffusivities (salt fingers in the subtropical ocean), rotation (Taylor vortices), shear against
   stratification (the Richardson number ¼ behind every ocean and atmosphere mixing scheme), and plain viscous shear
   (Tollmien–Schlichting waves) — learn the theorems that decide stability without solving anything (Rayleigh, Fjørtoft,
   Howard, Miles–Howard, Squire), compute every number with one numerical tool (Chebyshev collocation), and end with three
   equations of convection that never settle down: Lorenz's chaos, the reason weather has a predictability limit."
   **Roadmap (one line per CORE):** C01 normal modes and the stability vocabulary · C02 Kelvin–Helmholtz: shear vs gravity ·
   C03 the Bénard eigenproblem and its growth rate · C04 the rigid–rigid neutral curve and Ra_c ≈ 1708 · C05 free walls: the
   neutral curve by hand · C06 salt fingers: diffusion that destabilises · C07 Taylor vortices between rotating cylinders ·
   C08 the Taylor–Goldstein equation · C09 the Richardson number ¼ (Miles–Howard) · C10 Howard's semicircle · C11 Squire and
   the Orr–Sommerfeld equation · C12 Rayleigh's inflection-point theorem (and Fjørtoft) · C13 plane Poiseuille: viscosity
   destabilises (Re_c = 5772) · C14 the disturbance-energy budget (Reynolds-stress production) · C15 Lorenz: deterministic
   chaos. **Prerequisites:** potential flow and Laplace's equation (Ch. 6), interface waves (Ch. 7 §7.7, (7.95)), vortex
   sheets and their roll-up (Ch. 5 §5.8), unsteady Bernoulli (Ch. 4 §4.9), the Boussinesq equations (Ch. 4 §4.9), N² and
   stratification (Ch. 1 §1.10, Ch. 7 §7.8), circular Couette flow (Ch. 8 §8.2), Blasius, Falkner–Skan and jets (Ch. 9),
   perturb–linearise–eigenvalues (Ch. 9, P214), Fourier modes and amplification factors (Ch. 10 §10.2), the generalised
   eigenproblem (P247).
2. `nb.explainer_index([...])` — 9 rows: ("normal_mode_growth", "When is a flow 'unstable'?", "a disturbance is a sum of
   e^{ikx+σt} waves; the flow is unstable if σ_r > 0 for any k, and onset is where the σ_r(k) curve first touches zero") ·
   ("kelvin_helmholtz_boundary", "Why does wind raise waves only above a threshold?", "gravity holds long waves, surface tension
   short ones, shear pushes all; two real wave speeds collide and become a growing/decaying pair") · ("benard_neutral_curve",
   "Why does convection start near Ra ≈ 1708?", "each cell width has its own marginal Ra; onset is the bottom of that valley")
   · ("salt_fingers", "How can a stably stratified column overturn?", "heat diffuses 100× faster than salt, so a displaced
   parcel keeps its salt and keeps sinking") · ("taylor_couette_onset", "Why do stacked vortices appear between spinning
   cylinders?", "rings swap and release energy when Γ² falls outward; viscosity adds the same 1708-type threshold") ·
   ("richardson_shear_instability", "Why does Ri = ¼ decide whether a shear layer billows?", "Ri > ¼ everywhere guarantees
   stability; below ¼ instability becomes possible, not certain") · ("inviscid_shear_criteria", "Which profiles can be
   unstable without viscosity?", "an inflection point with a vorticity maximum, and every unstable c inside Howard's
   semicircle") · ("orr_sommerfeld_neutral_curve", "How can viscosity make a channel flow unstable?", "the viscous layer
   shifts v against u so the Reynolds stress draws energy from the shear: a thumb-shaped neutral curve above Re = 5772") ·
   ("lorenz_attractor", "How can three exact equations be unpredictable?", "past r_H ≈ 24.74 every steady state is unstable;
   nearby starts separate exponentially on a strange attractor").
3. `nb.setup()`.
4. `nb.code` — **chapter imports** (outside any block): `import numpy as np` · `import sympy as sp` · `import matplotlib.pyplot as
   plt` · `from scipy import linalg` · `from scipy.optimize import brentq, minimize_scalar` · `from scipy.integrate import
   solve_ivp` · `from fluidpy import ch11_instability as ch11` · `from fluidpy.core import stability as ST` · `from fluidpy.core
   import waves as WV, laminar as LAM, boundary_layer as BL, jets as JET, stratification as STRAT, kinematics as KIN` · `from
   fluidpy import ch05_vorticity_dynamics as ch05, ch07_gravity_waves as ch07` · `from fluidpy.core.interact import slider_figure,
   animate_figure, live` · `from fluidpy.core.anim import animate` · `from fluidpy.core.style import COLORS, savefig` · `from
   fluidpy.core.thermo import G0` · `import sys; sys.path.insert(0, "scripts"); from ch11_drawings import kh_sketch,
   benard_sketch, taylor_sketch, cv_sketch, semicircle_sketch, six_profiles`. *explain:* one line per import ("`ch11`
   re-exports the new eigen-solver toolkit `core/stability.py`, so `ch11.orr_sommerfeld_eigs` and `ST.orr_sommerfeld_eigs` are
   the same function"; "`scripts/ch11_drawings.py` only draws — no physics lives there").
5. `nb.md` — **⚠️ Conventions in this chapter**: the symbol table of convention 4 (book symbol · meanings · what we write); the
   twelve slips of convention 5 as a table "slip #k · the book prints · correct · where we fix it"; the four
   non-dimensionalisations of convention 6 as a four-row table; and the **two Γ conventions** as a two-row table (water layer
   d = 5 mm, ΔT = 2 K): (11.21) $\Gamma_T=-d\bar T/dz=+400$ K/m (heated from below ⇒ positive); Ch. 1's Kundu lapse rate
   $\Gamma\equiv dT/dz=-400$ K/m; meteorology's $\Gamma_{met}=-dT/dz=+400$ K/m — "negate the number, flip the inequality; the
   code never takes a Γ, it takes `dT` = T_bottom − T_top" (project memory, P48 recalled).
6. `nb.md` — **🔁 Tools from earlier chapters used in this one** (one line each, primer number): complex numbers and Euler's
   formula (P45), complex conjugate (P81), the complex plane in numpy (P153), complex powers and branch cuts (P155), complex
   square roots and the quadratic formula (P159), complex amplitudes (P176), linear second-order ODE by an exponential trial
   (P44), hyperbolic functions (P168), separation of variables for a PDE (P167), Taylor transfer of a boundary condition (P166),
   level sets (P75), operator elimination (P177), mean of a product of real parts (P178), eigenvalues (P80), the generalised
   symmetric eigenproblem (P247), collocation (P162), linear stability by eigenvalues (P214), inflection point (P212), sech
   (P215), integration by parts (P218a), integrals of sines over a period (P151), product rule (P38), chain rule (P49), Taylor
   expansion (P26/P98), orders of smallness (P68), order-of-magnitude scaling (P130), scaled variables (P133), anisotropic
   scaling (P188), cylindrical Laplacian (P186), determinants (P53, P56), null space (P58), sympy (P40) and its series (P117),
   `brentq` (P108), `minimize_scalar` (P170), `solve_ivp` (P31/P94), RK4 by hand (P95), `np.sign`/`np.nonzero` (P180), Fourier
   modes (P142), Simpson/trapezoid (P203), power laws and log–log axes (P13), `assert np.allclose` (P15), animate (P16),
   `slider_figure` (P17), `show_viz` (P18), plotly 3-D (P41/P64), live widgets (P47), caching runs (P252), reading reference data
   (P251), verification vs validation (P253), Python dictionaries (P23), lambda (P29). Each block repeats the ones it uses in a
   one-line reminder at first use.

---

### A.1 §11.1 Introduction — N01 N02 N03 (open the C01 story)
1. `nb.section("11.1", "Introduction", intro="**What is this section about?** A steady flow that satisfies the equations of
   motion is not necessarily a flow you will ever see. If a tiny disturbance — a ripple, a puff of warmer air, the vibration of a
   pump — grows instead of dying, the steady flow is *unstable* and something else takes its place: billows, convection cells,
   vortices, turbulence. This section sets up the question; §11.2 gives the method used for the rest of the chapter.")`
2. `nb.note` — **N01 [B]** "**Basic state and disturbance.** The *basic* (background) state is the steady flow we already know —
   two streams, a heated layer at rest, Couette flow, Poiseuille flow. We add a *disturbance* so small that products of
   disturbances can be dropped (the same move as Ch. 1's parcel argument for N² and Ch. 9's perturb–linearise–eigenvalues,
   P214). *Linear stability* asks whether infinitesimal disturbances grow; *finite-amplitude* (nonlinear) instability is when
   only large enough kicks grow. Both can happen to the same flow (the dimple below)."
3. `nb.note` — **N02 [B]** "**Four mechanical pictures** (our Fig. 11.1 analogue, animated below): a ball in a bowl (every kick
   dies: stable), on an upturned bowl (every kick grows: unstable), on a flat table (it rolls to a new rest: neutral), and in a
   small dimple on a hilltop (small kicks die, a large kick throws it out: stable to small, unstable to large disturbances)."
4. `nb.animation` — **A1** `potential_well_demo` for the four shapes, each with a small kick (x₀ = 0.3) and a large kick (x₀ =
   1.2), damping 0.3, 60 frames (FAST 30), `player="frames"` (curation §6). Layout 4 panels (2 × 2): the potential V(x) (muted)
   with the ball (rose if it is escaping, teal if returning) for both kicks. *explain:* "1. `potential_well_demo` integrates
   $\ddot x=-V'(x)-\gamma\dot x$ with `solve_ivp`; 2. `escaped` is True once ∣x∣ passes the rim." Notes: *see:* "two balls per
   panel"; *read:* "in the dimple the small kick rings down, the large one rolls away — stability can depend on the size of the
   disturbance"; *change:* "…we remove the damping: the bowl's ball oscillates forever (neutral in the energy sense) — the
   difference between 'neutral' and 'stable' of §11.2". *expect:* dimple escaped = (False, True); cap escaped for both.
5. `nb.note` — **N03 [C]** "**What this chapter is for, and what it leaves out.** Linear theory predicts *when* a flow first
   becomes unstable and *which wavelength* appears; it cannot say what the flow becomes afterwards (§11.12–11.14). We follow
   disturbances growing in *time* everywhere at once (temporal instability); disturbances that grow as they travel downstream
   (spatial instability, Huerre & Monkewitz 1990) are only named. No rotation here: the Coriolis force and baroclinic
   instability of the atmosphere and ocean come in Ch. 13 §13.17."

---

### A.2 §11.2 Method of Normal Modes — C01
#### C01 — The normal mode (11.1) and the stability vocabulary
1. `nb.section("11.2", "Method of Normal Modes", intro="**What is this section about?** One recipe used in every section of the
   chapter: write the disturbance as a sum of waves, follow each wave alone, and read off its growth rate. This section gives
   the recipe and the words — stable, neutral, unstable, marginal, stationary, oscillatory — that the rest of the chapter uses.")`
2. `nb.core("C01", "The normal mode $u=\\hat u(z)\\exp\\{\\mathrm ikx+\\mathrm imy+\\sigma t\\}=\\hat u(z)\\exp\\{\\mathrm
   i\\lvert\\mathbf K\\rvert(\\mathbf e_K\\cdot\\mathbf x-ct)\\}$ (11.1) and the stability vocabulary", question="How can one
   number per wavelength decide whether a whole flow is stable?")`
3. `nb.md` — **The problem in plain words:** "A light breeze blows over a pond. The surface is never perfectly flat: there are
   ripples of every length, from millimetres to metres, all at once. Some grow into waves, others die. To predict which, we would
   like to follow each ripple length separately — and for small ripples we can, because small disturbances obey *linear*
   equations, where a sum of solutions is a solution. A weather forecaster asks the same of a jet stream: which wavelength of
   meander grows fastest?"
4. `nb.md` — **The idea** (ASCII):
   ```
   any small disturbance  =  Σ over k   û(z) · e^{ikx} · e^{σ(k) t}          (Fourier, Ch. 10)
                                         ─────   ──────   ─────────
                                         shape   wave in x  grows or decays in time
   linear equations ⇒ each k evolves alone ⇒ one number σ(k) = σ_r + iσ_i per wavelength
     σ_r < 0 for every k   → stable        (every ripple dies)
     σ_r > 0 for some k    → unstable      (that wavelength grows like e^{σ_r t})
     σ_i ≠ 0               → the pattern travels/oscillates;  σ_i = 0 → it grows in place (cells)
   ```
   "Stability is not a property of one wave but of the whole curve σ_r(k): the flow is stable only if the curve stays below
   zero everywhere, and as a control parameter rises, the wavelength where the curve first touches zero is the pattern you see."
5. `nb.primer("necessary vs sufficient conditions, \"for every k\", and proof by contradiction", "Three pieces of logic the
   chapter's theorems use. *Sufficient*: 'if A then B' — A guarantees B (Ri > ¼ everywhere guarantees stability). *Necessary*:
   'B only if A' — without A no B, but A alone does not force B (an inflection point is needed for inviscid instability, but
   not enough). 'Stable' means σ_r ≤ 0 **for every** k; 'unstable' needs **one** k with σ_r > 0. *Proof by contradiction*:
   assume the opposite (a growing mode exists), derive something impossible (a positive number equal to a negative one), so
   the assumption was false.", code="import numpy as np\nsig_r = np.array([-0.3, -0.1, 0.02, -0.5])   # growth rates at four wavenumbers\nprint(np.all(sig_r < 0))                    # False: 'stable' needs EVERY k decaying\nprint(np.any(sig_r > 0))                    # True: one growing k is enough for 'unstable'")`
   (**P255**)
6. `nb.note` — **N05 [B]** "**Why one wave at a time is allowed.** The coefficients of the linearised equations do not depend
   on x, y or t (the basic state is uniform in those directions), so $e^{\mathrm i(kx+my)+\sigma t}$ passes through every
   derivative unchanged: $\partial/\partial x\to\mathrm ik$, $\partial/\partial t\to\sigma$ — the same trick as Ch. 10's
   Fourier error modes, where one mode was multiplied by G(θ) each step and here by $e^{\sigma\Delta t}$. Linearity lets a
   sum of modes evolve as the sum of their evolutions, so any disturbance is covered. What is left is an ODE in z with σ as an
   unknown — an **eigenvalue problem**." + `nb.code`: superposition demo on an interface: three modes k = 1, 2, 3 with σ = −0.2,
   0.1, −0.5 summed with `ST.normal_mode`, and the k = 2 term alone; print the amplitude at t = 10 of each. *expect:* e^{−2} =
   0.135, e^{1} = 2.718, e^{−5} = 0.0067 — "after a while only the growing mode is visible".
7. `nb.derivation("D01", …)` — Part F D01 (5 steps), ref "11.1". (Uses Euler's formula P45, complex amplitudes P176.)
8. `nb.note` — **N04 [B]** "**Stable, neutral, unstable** (one mode): $\sigma_r<0$ or $c_i<0$ stable; $\sigma_r=0$ or $c_i=0$
   neutrally stable; $\sigma_r>0$ or $c_i>0$ unstable. For the *flow*: stable if every k is stable, unstable if any k is (P255)."
   + `nb.code`: `[ST.stability_class(sigma=s) for s in (-0.2, 0.0, 0.3+1j)]` and `ST.stability_verdict(np.array([-0.3, -0.1,
   0.02, -0.5]))`. *expect:* ['stable', 'neutral', 'unstable']; verdict 'unstable' at index 2.
9. `nb.note` — **N06 [B]** "**Marginal state; stationary or oscillatory onset.** The *marginal* state sits on the border
   between stable and unstable: $\sigma_r=c_i=0$ **and** a small change of the control parameter (Ra, Re, Ta) makes σ_r > 0.
   (A *neutral* mode has σ_r = 0 but need not sit on such a border — a stable interface wave is neutral for every parameter
   value.) If also $\sigma_i=0$ at the margin, the instability appears as a **stationary** pattern (cells: Bénard, Taylor);
   if $\sigma_i\neq0$ it appears as growing **oscillations** (the book's 'oscillatory mode', older name 'overstability': the
   diffusive regime of double diffusion, C06)." A two-row table | onset | σ at the margin | examples |. + `nb.code`:
   `ST.marginal_type(0.0)`, `ST.marginal_type(0.0 + 2.5j)`. *expect:* 'stationary', 'oscillatory'.
10. `nb.worked_example("both forms of (11.1) for one mode", "Take K = (k, m) = (2, 0) m⁻¹ and c = 1 + 0.5i m/s. 1. Second form:
    $e^{\\mathrm iK(x-ct)}=e^{\\mathrm i2x}e^{-\\mathrm i2(1+0.5\\mathrm i)t}=e^{\\mathrm i2x}e^{(1-2\\mathrm i)t}$. 2. Compare
    with $e^{\\mathrm ikx+\\sigma t}$: σ = 1 − 2i s⁻¹, i.e. $\\sigma=-\\mathrm iKc$ (D01). 3. Growth rate σ_r = 1 s⁻¹ = K c_i =
    2 × 0.5 ✓: the amplitude multiplies by e ≈ 2.72 every second. 4. σ_i = −2 = −K c_r: the crests move at c_r = 1 m/s toward +x
    (a crest at kx − 2t = 0 moves right). 5. Verdict: unstable, oscillatory (travelling).")`
11. `nb.code` — `sig = ST.sigma_from_c(2.0, 1 + 0.5j)`; `print(sig, ST.c_from_sigma(2.0, sig), ST.stability_class(sigma=sig),
    ST.marginal_type(sig))`; then the field: `x = np.linspace(0, 2*np.pi, 200)`; `u_a = ST.normal_mode(x, 0.0, 1.0, 1.0, 2.0,
    sigma=sig)`; `u_b = ST.normal_mode(x, 0.0, 1.0, 1.0, 2.0, c=1+0.5j)`. *expect:* σ = (1−2j), c back = (1+0.5j), 'unstable',
    'oscillatory'; max∣u_a − u_b∣ ≈ 1e-16. *explain:* 1. (11.1) second form → first form; 2. classes of N04, N06; 3. the two
    field forms on a grid agree.
12. `nb.check_agree` — **from scratch (curation §7):** `u_mine = np.real(np.exp(1j*(2.0*x) + sig*1.0))` and `u_mine2 =
    np.real(np.exp(1j*2.0*(x - (1+0.5j)*1.0)))`; `assert np.allclose(u_mine, u_a) and np.allclose(u_mine2, u_b)`; parity with Ch. 7:
    `assert np.allclose(u_mine, WV.real_field(np.exp(sig.real*1.0), 2.0*x + sig.imag*1.0))`.
13. `nb.figure` — **"Growth is a property of the whole σ_r(k) curve"** (two panels, 9 × 3.4 in): (a) σ_r(k) for the static
    two-layer interface turned upside down (Rayleigh–Taylor, water over air with surface tension, `ch11.normal_mode_growth
    ("interface", k, rho1=1000, rho2=1.2, surface_tension=0.074)`) — rose above zero, teal below, the cut-off k marked; (b) a
    space–time diagram (x horizontal, t vertical, colour = η) of three modes summed: the fastest-growing k takes over. *see:* "a
    hump of growth rates, positive only for long waves"; *read:* "the top of the hump is the wavelength you will see; the zero
    crossing is the shortest wave that can grow (≈ 1.73 cm, N121)"; *change:* "…we remove surface tension: the hump never comes
    back down — every short wave grows faster (the KH and RT 'short waves always lose' of C02)".
14. `nb.explainer("normal_mode_growth", heading="When is a flow 'unstable'?", why="Dragging k along the σ_r(k) curve makes the
    interface grow, decay or travel at once, and raising the control parameter lifts the whole curve through zero — the wavelength
    that goes first appears in front of you; a static σ(k) plot hides the link between one point on the curve and the motion it
    means.", tries=["Pick 'upside-down interface' and drag k from long to short waves: watch the growth rate peak, then die where
    surface tension wins.", "Switch to 'Bénard (free walls)' and raise Ra slowly from 500: note the K where σ first crosses zero
    (≈ 2.22 at Ra ≈ 657.5).", "Open the complex σ-plane: at the onset the root crosses the imaginary axis — on the real axis for
    Bénard (stationary), off it for a travelling interface wave."])`
15. `nb.md` — **What would change if…** "…the basic state moved — two streams sliding past each other? Then c is no longer a pure
    imaginary or real number: it has a travelling part and a growing part at once, and the growth comes from the shear itself.
    That is the Kelvin–Helmholtz problem (C02), the first complete normal-mode calculation."

---

### A.3 §11.3 Kelvin–Helmholtz Instability — R01–R11, C02, S01
1. `nb.section("11.3", "Kelvin-Helmholtz Instability", intro="**What is this section about?** Two layers of fluid slide past
   each other: wind over the sea, a fast current over a slow one in the ocean, the warm air of a front over cold air. Gravity
   wants the interface flat (if the heavy fluid is below); the velocity jump wants to roll it up. We put a small wave on the
   interface, apply the conditions you met in Ch. 7, and get a formula for the wave speed whose square root turns imaginary when
   shear wins. The result explains billow clouds, wind-driven ripples and the start of mixing in the thermocline.")`
2. `nb.recap("R01", "Total potentials and Kelvin's theorem", "Each layer is irrotational: far upstream the flow is uniform, and
   Kelvin's theorem (Ch. 5 §5.2, (5.8): circulation of a material loop is constant in inviscid, barotropic flow with conservative
   forces) keeps it irrotational, so each layer has a potential, $\\tilde\\phi_1=U_1x+\\phi_1,\\ \\tilde\\phi_2=U_2x+\\phi_2$
   (11.2), with $\\phi_j$ the small disturbance.", where="Ch. 5 §5.2, Ch. 6")`
3. `nb.recap("R02", "Laplace's equation for the disturbances", "Incompressible + irrotational ⇒ $\\nabla^2\\phi_1=0,\\
   \\nabla^2\\phi_2=0$ (11.3) — exactly Ch. 6's potential flow and Ch. 7 §7.1.", where="Ch. 6, Ch. 7 §7.1")`
4. `nb.recap("R03", "Decay far away", "$\\phi_1\\to0$ as $z\\to+\\infty$ (11.4) and $\\phi_2\\to0$ as $z\\to-\\infty$ (11.5): the
   disturbance lives near the interface, as for Ch. 7's interface waves; this decides the signs $e^{\\mp kz}$ in D03.",
   where="Ch. 7 §7.7")`
5. `nb.recap("R04", "Kinematic and dynamic interface conditions", "Fluid on either side moves with the interface, $\\mathbf
   n\\cdot\\nabla\\tilde\\phi_1=\\mathbf n\\cdot\\mathbf U_s=\\mathbf n\\cdot\\nabla\\tilde\\phi_2$ on $z=\\zeta$ (11.6), and
   (no surface tension) the pressure is continuous, $p_1=p_2$ on $z=\\zeta$ (11.7) — the analogues of Ch. 7's (7.14) and
   (7.20). They are the starting line of D02.", where="Ch. 4 §4.10, Ch. 7 §7.1")`
6. `nb.recap("R05", "The normal of a level set", "The interface is the level set $f=z-\\zeta(x,t)=0$, so its unit normal is
   $\\mathbf n=\\nabla f/\\lvert\\nabla f\\rvert=[-(\\partial\\zeta/\\partial x)\\mathbf e_x+\\mathbf e_z]/\\sqrt{1+(\\partial
   \\zeta/\\partial x)^2}$, which turns (11.6) into (11.8) (Ch. 2 P75, Ch. 7 (7.16)).", where="Ch. 2 P75, Ch. 7 §7.1")`
7. `nb.recap("R06", "Unsteady Bernoulli in each layer", "Irrotational, unsteady: $\\frac{\\partial\\tilde\\phi_j}{\\partial t}+
   \\tfrac12\\lvert\\nabla\\tilde\\phi_j\\rvert^2+\\frac{p_j}{\\rho_j}+gz=C_j$ (11.10) — Ch. 4 (4.75); `core.bernoulli` computes
   it.", where="Ch. 4 §4.9")`
8. `nb.recap("R07", "A″ = k²A and its exponentials", "A constant-coefficient second-order ODE is solved by the trial $e^{mz}$
   (Ch. 1 P44): $m^2=k^2$, so $A=A_\\pm e^{\\pm kz}$ (separation of variables, Ch. 7 P167).", where="Ch. 1 P44, Ch. 7 P167")`
9. `nb.recap("R08", "The decaying solutions", "With (11.4)–(11.5): $\\phi_1=A_-\\exp\\{\\mathrm ik(x-ct)-kz\\}$, $\\phi_2=A_+
   \\exp\\{\\mathrm ik(x-ct)+kz\\}$ (11.15) — the same potentials as Ch. 7's `interface_fields`.", where="Ch. 7 §7.7")`
10. `nb.recap("R09", "The static interface (U₁ = U₂ = 0)", "Then $c=\\pm\\big[\\frac{\\rho_2-\\rho_1}{\\rho_2+\\rho_1}\\frac
    gk\\big]^{1/2}$ (11.19): interface waves of Ch. 7, $\\omega^2=gk\\frac{\\rho_2-\\rho_1}{\\rho_2+\\rho_1}$ (7.95); neutral if the
    heavy fluid is below, Rayleigh–Taylor unstable if it is above. A parity cell checks `ch11.kh_phase_speed` against
    `WV.interface_omega`.", where="Ch. 7 §7.7")` + **slip #5 box**: "the book writes 'see (7.96)'; the interface dispersion
    relation is Ch. 7 (7.95), $\\omega^2=gk(\\rho_2-\\rho_1)/(\\rho_2+\\rho_1)$ ((7.96) is the energy)."
11. `nb.recap("R10", "The vortex sheet", "Equal densities ρ₁ = ρ₂: the interface is a vortex sheet of strength γ = U₂ − U₁ (Ch. 5
    §5.8, γ = u_below − u_above). **New here:** (11.18) becomes $c=\\frac{U_2+U_1}2\\pm\\mathrm i\\frac{U_2-U_1}2$ (11.20), so
    the growth rate $kc_i=k\\lvert U_2-U_1\\rvert/2$ is positive at **every** k — a vortex sheet is unstable to all wavelengths,
    fastest for the shortest.", where="Ch. 5 §5.8")`
12. `nb.recap("R11", "Roll-up of a perturbed vortex sheet", "Ch. 5's `sheet_rollup` (point vortices with Krasny's δ-smoothing)
    follows the sheet past the linear stage into cat's-eye rolls (our Fig. 11.6 analogue). Animation A2 below puts the linear
    growth $e^{k\\Delta U t/2}$ next to it.", where="Ch. 5 §5.8")`

#### C02 — The Kelvin–Helmholtz dispersion relation (11.18) and its stability boundary
13. `nb.core("C02", "The Kelvin–Helmholtz dispersion relation $c=\\frac{\\rho_2U_2+\\rho_1U_1}{\\rho_2+\\rho_1}\\pm\\Big[\\frac{\\rho_2
    -\\rho_1}{\\rho_2+\\rho_1}\\frac gk-\\frac{\\rho_2\\rho_1}{(\\rho_2+\\rho_1)^2}(U_2-U_1)^2\\Big]^{1/2}$ (11.18)",
    question="When does a velocity jump across an interface make waves grow, and which ones?")`
14. `nb.md` — **The problem in plain words:** "Blow gently across a glass of water: nothing. Blow harder: ripples appear and grow.
    In the sky, a fast layer of air over a slow, denser one sometimes paints a row of breaking-wave 'billow clouds'. In the ocean
    thermocline, dye released by divers rolls up into the same shapes. All are one instability: shear across an interface. We
    want the threshold — how fast must the wind be? — and the wavelength that grows."
15. `nb.md` — **The idea** (ASCII):
    ```
          U₁, ρ₁ (light)   ───────────►            over a crest the upper stream is squeezed → faster → lower pressure
       ~~~~~~~~~~~~~~~~~ interface ζ = ζ₀ e^{ik(x−ct)} ~~~~   ⇒ the crest is SUCKED up (Bernoulli): shear destabilises
          U₂, ρ₂ (heavy)   ──►                     gravity pulls the crest back down: stratification stabilises
       c = mean speed ± √( gravity term − shear term )   →  √(negative) = imaginary ⇒ one root grows
    ```
    "Everything is decided by the sign under the square root: restoring (gravity, later surface tension) minus driving (the
    kinetic energy of the shear)."
16. `nb.note` — **N07 [B]** "**The set-up** (our Fig. 11.2 analogue below): upper fluid speed U₁ [m/s], density ρ₁ [kg/m³];
    lower U₂, ρ₂; interface z = ζ(x, t); inviscid, irrotational in each layer, no surface tension (added in N120)." + `nb.figure`
    `kh_sketch(ax)` driven by `ch11.kh_fields` (interface and velocity arrows; layer 1 orange, layer 2 blue). *see:* "two streams,
    a wavy interface"; *read:* "arrows longer above crests in the upper layer (squeezed streamlines)"; *change:* "…U₁ = U₂: no
    squeezing difference — Ch. 7's interface wave".
17. `nb.derivation("D02", …)` — Part F D02 (12 steps), ref "11.13" — contains **N08** (exact kinematic condition, step 5),
    **N09** (11.9) (step 7), **N10** (11.11) (step 8), **N11** (11.12) (step 9), **N12** (11.13) (step 12). (Uses R04, R05, R06,
    Taylor transfer P166, orders of smallness P68.)
18. `nb.note` — **N13 [B]** "**The normal-mode trial** (C01's mode with K = (k, 0, 0)): $\\phi_j=A_j(z)\\exp\\{\\mathrm ik(x-ct)
    \\}$ (11.14) and $\\zeta=\\zeta_o\\exp\\{\\mathrm ik(x-ct)\\}$; inserted into (11.3) it gives R07's $A''=k^2A$."
19. `nb.primer("np.roots and np.lib.scimath.sqrt", "`np.roots(coeffs)` returns all roots of a polynomial given its
    coefficients, highest power first (complex if needed). `np.sqrt(-1.0)` gives `nan` with a warning, but
    `np.lib.scimath.sqrt(-1.0)` returns the principal complex root `1j` — what we want when a discriminant can turn negative.",
    code="import numpy as np\nprint(np.roots([1, -3, 4]))                 # c^2 - 3c + 4 = 0 -> 1.5 ± 1.3229j\nprint(np.lib.scimath.sqrt(-1.75))           # 1.3229j: the principal square root of a negative number")`
    (**P256**)
20. `nb.derivation("D03", …)` — Part F D03 (14 steps), ref "11.18" — contains **N14** (11.16)–(11.17) (steps 5–7), **N15** the
    quadratic (step 9), **N16** the instability inequality (step 13), **R09**/(11.19), **R10**/(11.20) and **N19** (step 14,
    c_r = mean); ★★ with an optional `check_src` (`ch11.kh_sympy()["identity"] == 0`).
21. `nb.note` — **N16 [B]** "**The stability boundary.** Unstable exactly when the square root's argument is negative:
    $\\frac{\\rho_2-\\rho_1}{\\rho_2+\\rho_1}\\frac gk<\\frac{\\rho_2\\rho_1}{(\\rho_2+\\rho_1)^2}(U_2-U_1)^2$, i.e. $g(\\rho_2^2-
    \\rho_1^2)<k\\rho_1\\rho_2(U_2-U_1)^2$ — unstable if the velocity jump is large, the density jump small, or the wave short:
    $k>k_c=g(\\rho_2^2-\\rho_1^2)/(\\rho_1\\rho_2\\Delta U^2)$." + `nb.code`: `ch11.kh_critical_k(5.0, 0.0, 1.2, 1000.0)`.
    *expect:* 327 m⁻¹, λ_c = 2π/k_c = 1.92 cm.
22. `nb.note` — **N17 [B]** "**Growing and decaying twins.** (11.18) has real coefficients, so its two roots are complex
    conjugates when they are not real: c₊ = c₋*. For every growing wave there is a decaying one; a nonzero c_i therefore always
    means instability (stated, a-D20; shown on E2's c-plane)."
23. `nb.note` — **N18 [B]** "**Short waves always lose.** Without surface tension k_c is finite for any ΔU ≠ 0, so all waves
    shorter than 2π/k_c grow — 'the flow is always unstable to short waves when U₁ ≠ U₂' — and the growth rate kc_i rises without
    bound as k grows (figure F1). Surface tension (N120) cures it."
24. `nb.worked_example("a two-layer flow with easy numbers", "ρ₁ = 1, ρ₂ = 3 kg/m³ (light over heavy), g = 10 m/s², k = 1 m⁻¹,
    U₂ = 0. **Case U₁ = 4 m/s.** 1. Mean speed $(\\rho_2U_2+\\rho_1U_1)/(\\rho_1+\\rho_2)=4/4=1$ m/s. 2. Gravity term
    $\\frac{3-1}{4}\\cdot\\frac{10}{1}=5$ m²/s². 3. Shear term $\\frac{3\\cdot1}{16}\\cdot16=3$ m²/s². 4. Under the root 5 − 3 = 2
    > 0 ⇒ c = 1 ± √2 = 2.414 or −0.414 m/s: two real speeds, neutral waves. **Case U₁ = 6 m/s.** 5. Mean 6/4 = 1.5 m/s. 6.
    Shear term $\\frac3{16}\\cdot36=6.75$. 7. Under the root 5 − 6.75 = −1.75 < 0 ⇒ c = 1.5 ± 1.3229i m/s. 8. Growth rate kc_i =
    1.3229 s⁻¹: the amplitude grows by e every 0.76 s. 9. Check with the inequality: g(ρ₂² − ρ₁²) = 10 × 8 = 80; kρ₁ρ₂ΔU² = 3 × 16
    = 48 (stable) and 3 × 36 = 108 (unstable) ✓; k_c = 80/108 = 0.741 m⁻¹ for ΔU = 6.")`
25. `nb.code` — `cp, cm = ch11.kh_phase_speed(1.0, 6.0, 0.0, 1.0, 3.0, g=10.0)`; `print(cp, cm)`;
    `print(ch11.kh_discriminant_terms(1.0, 6.0, 0.0, 1.0, 3.0, g=10.0))`; `print(ch11.kh_phase_speed(1.0, 4.0, 0.0, 1.0, 3.0,
    g=10.0))`; then air over water: `k = np.logspace(1, 4, 400)` and `ch11.kh_growth_rate(k, 5.0, 0.0, 1.2, 1000.0)`. *expect:*
    (1.5+1.3229j), (1.5−1.3229j); dict(gravity 5.0, tension 0.0, shear −6.75, total −1.75, mean 1.5); (2.4142, −0.4142); growth
    zero below k_c = 327 m⁻¹. *explain:* 1. (11.18) with our tiny numbers; 2. the three pieces under the root (E2's bars);
    3. the neutral case; 4. a whole growth curve for air over water.
26. `nb.check_agree` — **from scratch (curation §7):** expand the quadratic ρ₁(U₁ − c)² + ρ₂(U₂ − c)² = (g/k)(ρ₂ − ρ₁) into
    $(\\rho_1+\\rho_2)c^2-2(\\rho_1U_1+\\rho_2U_2)c+(\\rho_1U_1^2+\\rho_2U_2^2)-(g/k)(\\rho_2-\\rho_1)=0$ and solve with `np.roots`
    (P256): `mine = np.roots([4.0, -2*6.0, 36.0 - 10*2.0])`; `assert np.allclose(sorted(mine, key=np.imag), sorted([cp, cm],
    key=np.imag))`. Parity with Ch. 7 (R09): `assert np.isclose(np.real(ch11.kh_phase_speed(2.0, 0.0, 0.0, 1.0, 3.0, g=9.81)[0]),
    WV.interface_omega(2.0, 1.0, 3.0, g=9.81)/2.0)`; residual check `ch11.kh_residuals(1.0, 6.0, 0.0, 1.0, 3.0, g=10.0)` all <
    1e-12.
27. `nb.figure` — **"Which waves grow: the (k, ΔU) plane"** (two panels): (a) the stability boundary ΔU_min(k) for air over water
    without (orange dashed) and with surface tension (purple, its minimum ◆ at 6.70 m/s, λ = 1.73 cm), the unstable region rose,
    log k axis, our two winds (5 and 8 m/s) as horizontal lines; (b) growth rate kc_i(k) for ΔU = 8 m/s with and without surface
    tension; the band 0.71–4.21 cm shaded. *see:* "a V-shaped boundary with surface tension, a falling line without"; *read:* "a
    horizontal wind line inside the rose region = waves of those lengths grow; at 5 m/s with surface tension there are none";
    *change:* "…the densities differ by only 0.2 % (an ocean thermocline, no surface tension): with ΔU = 0.1 m/s every wave shorter than ≈ 1.6 m grows (k_c = 3.93 m⁻¹)".
28. `nb.plotly` — **F1** `slider_figure` over ΔU (24 values 1 … 12 m/s, FAST 12): traces kc_i(k) without and with surface tension
    (air over water), the k_c marker. *explain:* "the page version of E2's growth panel".
29. `nb.derivation("D04", …)` — Part F D04 (9 steps), ref "11.18" (Exercise 11.1; the book gives only the result) — states **N120**.
30. `nb.note` — **N120 [B]** "**Depth and surface tension (Exercise 11.1).** With a lower layer of depth h and surface tension
    σ_s [N/m]: $c=\\frac{\\rho_1U_1+\\rho_2U_2\\coth kh}{\\rho_1+\\rho_2\\coth kh}\\pm\\Big[\\frac{(g/k)(\\rho_2-\\rho_1)+\\sigma_s k}
    {\\rho_1+\\rho_2\\coth kh}-\\frac{\\rho_1\\rho_2(U_1-U_2)^2\\coth kh}{(\\rho_1+\\rho_2\\coth kh)^2}\\Big]^{1/2}$; h → ∞ gives
    (11.18) plus the σ_s k term. Minimum wind for deep water: $\\Delta U_{\\min}^2=2\\sqrt{g\\Delta\\rho\\,\\sigma_s}(\\rho_1+\\rho_2)/
    (\\rho_1\\rho_2)$ at $k_*=\\sqrt{g\\Delta\\rho/\\sigma_s}$." + `nb.code`: `ch11.kh_min_shear(1.2, 1000.0)` and
    `ch11.kh_unstable_band(8.0, 1.2, 1000.0, surface_tension=0.074)`. *expect:* dU_min 6.70 m/s, k* 364 m⁻¹, wavelength 1.73
    cm; band 149–887 m⁻¹ (0.71–4.21 cm); "real wind waves start at lower winds — the turbulent wind's pressure fluctuations
    (Miles, Phillips) do the job KH cannot; named only".
31. `nb.note` — **N121 [B]** "**Rayleigh–Taylor with surface tension (Exercise 11.2, our answer).** Heavy fluid above (paint on
    a ceiling): from (11.19) with ρ₁ > ρ₂ every wave grows unless surface tension holds the short ones; the longest wave that
    stays neutral is $\\lambda_c=2\\pi\\sqrt{\\sigma_s/((\\rho_h-\\rho_l)g)}$ — the capillary length times 2π." + `nb.code`:
    `ch11.rayleigh_taylor_cutoff(0.074, 1000.0, 1.2)`. *expect:* 0.0173 m: "a water film on a ceiling holds only for patches
    smaller than ≈ 1.7 cm".
32. `nb.note` — **N19 [B]** "**The moving frame** (our Fig. 11.3 analogue): seen from a frame moving at (U₁ + U₂)/2 (Ch. 3's
    Galilean shift, `KIN.galilean_transform`), the two streams are ±ΔU/2 — a symmetric picture with no preferred direction, so
    the unstable wave cannot travel: c_r equals the mean speed when ρ₁ = ρ₂." + small figure: the two velocity profiles in the lab
    and moving frames side by side (*see/read/change* short: "symmetric → c_r = 0 in that frame; with ρ₁ ≠ ρ₂ the
    density-weighted mean of (11.18) replaces the plain mean").
33. `nb.note` — **N20 [C]** "**Shear against stratification in nature.** Laboratory tilting-tube experiments (Thorpe 1971),
    billow clouds and dye in ocean thermoclines (Woods 1969) all show this instability; it is a main source of internal waves and
    of mixing across density interfaces. With *continuous* stratification the question becomes the Richardson number — C09.
    Climate hook: KH billows in the thermocline and in the stable night-time atmosphere are why every mixing parameterisation
    switches on below Ri ≈ ¼ (Ch. 12, Ch. 13)."
34. `nb.animation` — **A2** (curation §6): left, the linear KH interface $\\zeta_0e^{k\\Delta Ut/2}\\cos(kx)$ (equal densities,
    frame moving with the mean); right, `ch05.sheet_rollup(N=200 (FAST 100), gamma=1.0, amplitude=0.01, delta=0.05)` frames on
    the same clock with the linear envelope overlaid (muted dashed); 60 frames (FAST 30), `player="video"`. *explain:* "1. the
    linear mode grows exponentially forever; 2. the point-vortex sheet follows it early, then rolls up into cat's-eyes when the
    slope is O(1)". Notes: *see:* "two growing waves that agree, then part ways"; *read:* "linear theory is right until the
    amplitude is about a tenth of the wavelength"; *change:* "…we halve δ (less smoothing): the roll-up core tightens, the
    linear stage is the same".
35. `nb.note` — **N21 [C]** + **N22 [B]** "**Where the energy comes from.** The billows feed on the kinetic energy of the two
    streams (our Fig. 11.7 analogue: before and after profiles). Mixing example (p. 482): a stream U₁ over still fluid, both of
    thickness h, mixed into a linear profile $U(z)=U_1(\\tfrac12+\\tfrac z{2h})$, $-h\\le z\\le h$, keeps its momentum
    ∫ρU dz = ρU₁h but its kinetic energy drops from $E_i=\\tfrac\\rho2U_1^2h$ to $E_f=\\tfrac\\rho2\\int_{-h}^hU^2dz=\\tfrac\\rho3U_1^2
    h$ — a third is released, available to lift heavy fluid (mixing raises the potential energy)." + `nb.code`:
    `ch11.kh_mixing_energy(1.0, 1.0, 1000.0)` with a `quad` cross-check of ∫U². *expect:* ratio 0.6667; M_i = M_f = 1000 kg/s per
    m. + bar figure (E_i, E_f, released, orange/muted).
36. `nb.note` — **N23 [C]** "**A general rule** (named): for a fixed momentum ∫U dz, smoothing the velocity gradients always lowers
    ∫U² dz (Cauchy–Schwarz: $\\int U^2dz\\ge(\\int U\\,dz)^2/(2h)$), so mixing a shear layer always releases kinetic energy. C14
    shows how a growing wave taps that energy (the Reynolds stress)."
37. `nb.explainer("kelvin_helmholtz_boundary", heading="Why does wind raise waves only above a threshold?", why="Two sliders (ΔU,
    density ratio) and two toggles (surface tension, finite depth) move both the stability boundary and the current point; on the
    c-plane the two real wave speeds slide together, collide and split into a growing/decaying pair exactly when the point crosses
    the boundary — the collision is the instability, and only motion shows it.", tries=["Preset 'air over water 5 m/s': stable.
    Drag ΔU up and stop when the status turns rose — read ΔU_min ≈ 6.7 m/s and λ ≈ 1.7 cm.", "Watch the c-plane while you cross
    the boundary: two dots on the real axis meet and leave it vertically.", "Turn surface tension off: the boundary becomes a
    falling line — every short wave grows (N18)."])`
38. `nb.pointer("**S01** — Exercises 11.3–11.5 (a porous surface, a compliant surface, membrane flutter) extend the interface
    conditions of C02 to other boundaries; they are not solved here.")`
39. `nb.md` — **What would change if…** "…instead of a velocity jump the fluid were at rest and heated from below? Then there is
    no shear to tap; the energy source is buoyancy, and viscosity and heat diffusion fight it. The answer is a threshold number,
    the Rayleigh number — C03."

---

### A.4 §11.4 Thermal Instability: The Bénard Problem — R12 R13, C03 · C04 · C05
1. `nb.section("11.4", "Thermal Instability: The Bénard Problem", intro="**What is this section about?** Heat a layer of fluid
   from below: the bottom fluid expands, becomes lighter, and wants to rise — yet a thin layer of honey on a warm plate stays
   still. Viscosity and heat diffusion resist; only when the heating beats them, measured by one number, the Rayleigh number,
   does the layer overturn into convection cells. We set up the linear problem (C03), find the threshold for rigid walls,
   Ra ≈ 1708 (C04), and solve the stress-free case by hand (C05). The same physics drives shallow atmospheric convection on a
   sunny day, cloud streets, and convection in the ocean's mixed layer.")`
2. `nb.recap("R12", "The Boussinesq equations", "Density changes are kept only where they meet gravity (Ch. 4 §4.9): $\\nabla
   \\cdot\\tilde{\\mathbf u}=0$, $\\frac{\\partial\\tilde{\\mathbf u}}{\\partial t}+(\\tilde{\\mathbf u}\\cdot\\nabla)\\tilde{\\mathbf
   u}=-\\frac1{\\rho_0}\\nabla\\tilde p-g[1-\\alpha(\\tilde T-T_0)]\\mathbf e_z+\\nu\\nabla^2\\tilde{\\mathbf u}$, $\\frac{\\partial
   \\tilde T}{\\partial t}+(\\tilde{\\mathbf u}\\cdot\\nabla)\\tilde T=\\kappa\\nabla^2\\tilde T$ (4.10, 4.86, 4.89), with
   $\\rho=\\rho_0[1-\\alpha(\\tilde T-T_0)]$. They are the starting line of D05.", where="Ch. 4 §4.9")`
3. `nb.recap("R13", "The conduction state", "With no motion the momentum equation is hydrostatic and the heat equation is
   $0=\\kappa\\,\\partial^2\\bar T/\\partial z^2$: $0=-\\frac1{\\rho_0}\\nabla P-g[1-\\alpha(\\bar T-T_0)]\\mathbf e_z$ and
   $0=\\kappa\\frac{\\partial^2\\bar T}{\\partial z^2}$ (11.23) — Ch. 1's hydrostatics and Ch. 8's steady conduction: a straight
   temperature line between the plates.", where="Ch. 1 §1.7, Ch. 8")`

#### C03 — The Bénard amplitude problem (11.36)–(11.37) and its growth rate σ(K, Ra)
4. `nb.core("C03", "The Bénard amplitude problem $\\big(\\sigma+K^2-\\frac{d^2}{dz^2}\\big)\\hat T=W$ and $\\big(\\frac\\sigma\\Pr
   +K^2-\\frac{d^2}{dz^2}\\big)\\big(\\frac{d^2}{dz^2}-K^2\\big)W=-\\mathrm{Ra}K^2\\hat T$ (11.36, 11.37)", question="A layer
   heated from below: which disturbances grow, and how fast?")`
5. `nb.md` — **The problem in plain words:** "A pan of water on a stove, a 5 mm layer of oil on a hot plate, the lowest kilometre
   of the atmosphere on a sunny morning: warm fluid underneath, cool on top. A blob nudged upward finds itself warmer and lighter
   than its new surroundings and keeps rising — unless viscosity slows it and heat leaks out of it first. We want to know, for a
   given layer, whether that tug-of-war is won, and if so how fast the convection starts."
6. `nb.md` — **The idea** (ASCII):
   ```
   T₀ − ΔT  ─────────────── cold plate
             ↑ warm blob: buoyancy  gαT′           pushes it up      (∝ ΔT, d)
             │ viscosity ν∇²w                       holds it back
             │ heat leaks out: κ∇²T′                erases T′
   T₀       ─────────────── hot plate
   one number compares them:  Ra = buoyancy / (viscous × diffusive) = gαΓd⁴/(κν)
   small disturbance + normal modes in x, y  ⇒  two ODEs in z  ⇒  σ is an EIGENVALUE that depends on (K, Ra, Pr)
   ```
7. `nb.note` — **N24 [C]** "**History.** Bénard's 1900 hexagons in millimetre-thin layers with a free surface were mostly driven
   by the variation of surface tension with temperature (Marangoni convection), not by buoyancy (Drazin & Reid 1981); Rayleigh
   (1916) solved the buoyancy problem with free surfaces, Jeffreys (1928) with rigid ones. The patterns themselves are N51."
8. `nb.note` — **N25 [B]** "**The Rayleigh number** $\\mathrm{Ra}=g\\alpha\\Gamma d^4/\\kappa\\nu$ (11.21), with g [m/s²], α
   thermal expansion [1/K], $\\Gamma=-d\\bar T/dz$ [K/m] (positive when heated from below), d depth [m], κ thermal diffusivity
   [m²/s], ν kinematic viscosity [m²/s] — dimensionless. Since Γ = ΔT/d (N28), Ra = gαΔT d³/(κν)." + slip #10 box: **⚠️ Common
   confusion — which Γ?** the two-row table of front-matter item 5 for this layer, plus "(11.21) uses the *opposite* sign to
   Ch. 1's Kundu lapse rate Γ ≡ dT/dz; here Γ happens to agree with the meteorological −dT/dz. The code avoids the trap: it
   takes `dT` = T_bottom − T_top." + `nb.code`: `ch11.rayleigh_number(alpha=2.1e-4, dT=2.0, d=0.005, kappa=1.4e-7, nu=1.0e-6,
   g=9.81)`; `ch11.gamma_conventions(2.0, 0.005)`; a planted wrong sign `ch11.rayleigh_number(..., dT=-2.0, ...)`. *expect:* 3679;
   dict(400, −400, 400); −3679 "(negative Ra = heated from above: stable for every disturbance)".
9. `nb.note` — **N26 [B]** "**Geometry** (our Fig. 11.8 analogue): a layer of depth d between two plates, z measured from the
   middle (−d/2 ≤ z ≤ d/2), bottom at T₀, top at T₀ − ΔT." + `nb.figure` `benard_sketch(ax)` (the layer with $\\bar T(z)$ beside it,
   blue). *see / read / change:* "a straight temperature line, hot at the bottom" / "Γ is its (positive) slope downward" / "…heat
   from above: the line tilts the other way, Ra < 0, nothing can grow (D08's identity then says nothing; energy arguments show
   stability)".
10. `nb.note` — **N27 [B]** + **N28 [B]** "**Base state plus disturbance:** $\\tilde{\\mathbf u}=0+\\mathbf u(\\mathbf x,t)$,
    $\\tilde T=\\bar T(z)+T'(\\mathbf x,t)$, $\\tilde p=P(z)+p(\\mathbf x,t)$ (11.22); integrating (11.23) gives the conduction
    profile $\\bar T(z)=T_0-\\tfrac12\\Delta T-\\Gamma z$, $\\Gamma\\equiv\\Delta T/d$ (11.24)." + `nb.code`:
    `ch11.benard_base_state(np.array([-0.0025, 0.0, 0.0025]), T0=300.0, dT=2.0, d=0.005)`. *expect:* T = 300, 299, 298 K.
11. `nb.derivation("D05", …)` — Part F D05 (8 steps), ref "11.27" — contains **N29** (the nonlinear perturbation equations with
    −wΓ, step 6) and **N30** (11.25)–(11.27) (step 8). (Uses R12, R13, orders of smallness P68.)
12. `nb.note` — **N31 [B]** "**Where Ra comes from (scaling).** In (11.27) $\\frac{\\partial T'}{\\partial t}-w\\Gamma=\\kappa\\nabla^2
    T'$ balance $w\\Gamma\\sim\\kappa\\Delta T/d^2$ with $T'\\sim\\Delta T$, $\\nabla\\sim1/d$: $w\\sim\\kappa/d$. Then in (11.26)
    buoyancy over viscous force $\\sim\\frac{g\\alpha T'}{\\nu w/d^2}\\sim\\frac{g\\alpha\\Gamma d^4}{\\nu\\kappa}=\\mathrm{Ra}$
    (P130 recalled)." + `nb.code`: `ch11.benard_scales(alpha=2.1e-4, dT=2.0, d=0.005, kappa=1.4e-7, nu=1e-6, g=9.81)` and the
    air layer `ch11.rayleigh_number(alpha=1/293, dT=1.0, d=0.01, kappa=2.1e-5, nu=1.5e-5, g=9.81)`. *expect:* w ~ 2.8e-5 m/s,
    t ~ 179 s, Pr = 7.14; air 106 per kelvin ⇒ ≈ 16 K needed for Ra_c ≈ 1708 in a 1 cm air gap, < 1 K for 5 mm of water.
13. `nb.derivation("D06", …)` — Part F D06 (8 steps), ref "11.29" — contains **N32** (11.28), **N33** (the pressure Poisson
    equation and its z-derivative), **N34** (11.29). (Uses ∇·, commuting operators P121, Ch. 10's pressure Poisson D19 recalled.)
14. `nb.note` — **N35 [B]** "**Rigid isothermal walls:** $w=\\partial w/\\partial z=T'=0$ on $z=\\pm d/2$ (11.30) — no slip gives u =
    v = w = 0 on the whole plate, so ∂u/∂x = ∂v/∂y = 0 there, and continuity forces ∂w/∂z = 0."
15. `nb.md` — **⚠️ Four ways this chapter makes things dimensionless** — the four-row table of convention 6 (§11.4 · §11.6 ·
    §11.8–11.11 · §11.14), with "we are here: §11.4 — lengths by d, time by d²/κ, w left dimensional until the very last step".
16. `nb.derivation("D07", …)` — Part F D07 (9 steps), ref "11.37" — contains **N36** (11.31)–(11.33) with Pr ≡ ν/κ, **N37**
    (operator substitutions) and **N38** (11.34)–(11.35). (Uses scaled variables P133, operator elimination P177.)
17. `nb.primer("eigenvalue problem for a differential equation", "In Ch. 2 an eigenvalue problem was a matrix equation Av = λv
    with a nonzero v. Here it is an ODE with boundary conditions, like −u″ = λu with u(0) = u(π) = 0: for most λ the only
    solution is u = 0; only for special λ (here λ = 1, 4, 9, …, u = sin nz) does a nonzero solution meet every boundary
    condition. Those λ are the eigenvalues. In (11.36)–(11.37) the growth rate σ plays λ: only special σ allow a nonzero (W, T̂).",
    code="import numpy as np\nz = np.linspace(0, np.pi, 5)\nfor n in (1, 2):                                  # u = sin(n z) meets u(0) = u(pi) = 0\n    print(n**2, np.allclose(np.sin(n*z)[[0, -1]], 0))   # eigenvalue n^2 = 1, 4 and True")` (**P257**)
18. `nb.note` — **N39 [B]** "**Boundary conditions for the amplitudes:** $W=\\partial W/\\partial z=\\hat T=0$ on $z=\\pm1/2$
    (11.38). With (11.36)–(11.37) this is a sixth-order problem in z whose eigenvalue is σ."
19. `nb.primer("Chebyshev–Gauss–Lobatto points and the differentiation matrix", "To turn an ODE into a matrix problem we sample
    the unknown at N + 1 points $x_j=\\cos(j\\pi/N)$ — crowded near the ends — and replace d/dx by a matrix D that differentiates
    the polynomial through those values exactly. Because the points are Chebyshev, the error falls *exponentially* with N for
    smooth solutions (spectral accuracy) — far faster than Ch. 10's O(h²) stencils. D² is D @ D; z ∈ [−½, ½] just rescales D by
    2.", code="import numpy as np\nfrom fluidpy.core import stability as ST\nD, x = ST.cheb(4)                         # 5 points: 1, 0.707, 0, -0.707, -1\nprint(np.max(np.abs(D @ x**3 - 3*x**2)))   # ~4e-16: exact for polynomials up to degree N")` (**P258**)
20. `nb.primer("boundary-row replacement and the generalised non-symmetric eigenproblem scipy.linalg.eig(A, B)", "Collocation
    turns the ODEs into $A\\mathbf v=\\sigma B\\mathbf v$. Boundary conditions replace the rows at the end points: in A the row
    becomes the condition (e.g. [1, 0, …, 0] for W = 0), in B the row becomes zero — so those rows say '0 = condition' for every
    σ. B is then singular and `scipy.linalg.eig(A, B)` returns some infinite eigenvalues (drop them) and, at larger N, a few huge
    spurious ones caused by round-off (keep only eigenvalues that do not move when N grows). Unlike Ch. 10's `eigh` (P247), A and
    B are not symmetric, so eigenvalues can be complex.", code="import numpy as np\nfrom scipy.linalg import eig\nA = np.array([[2.0, 1.0], [1.0, 1.0]]); B = np.array([[1.0, 0.0], [0.0, 0.0]])   # second row: a 'boundary row' (B row zero)\nprint(eig(A, B, right=False))                # [1, inf]: the boundary row gives an infinite eigenvalue")` (**P259**)
21. `nb.worked_example("a 5 mm water layer on a hot plate", "d = 5 mm, ΔT = 2 K, α = 2.1 × 10⁻⁴ K⁻¹, κ = 1.4 × 10⁻⁷ m²/s, ν = 1.0 ×
    10⁻⁶ m²/s, g = 9.81 m/s². 1. Γ = ΔT/d = 400 K/m. 2. Ra = gαΔT d³/(κν) = 9.81 × 2.1e-4 × 2 × 1.25e-7/(1.4e-7 × 1e-6) = 5.15e-10/
    1.4e-13 ≈ 3679. 3. Pr = ν/κ = 7.14. 4. Time unit d²/κ = 25e-6/1.4e-7 = 179 s. 5. The code below finds the largest growth
    rate σ ≈ 19.9 (in units of κ/d²) at K = 3.12: dimensional σ = 19.9/179 s ≈ 0.111 s⁻¹, so a disturbance grows by e every ≈ 9 s.
    6. Halve ΔT to 1 K: Ra = 1839, just above 1708 — the growth almost stops (C04).")`
22. `nb.code` — `Ra = ch11.rayleigh_number(alpha=2.1e-4, dT=2.0, d=0.005, kappa=1.4e-7, nu=1e-6, g=9.81)`; `Pr = 1e-6/1.4e-7`;
    `sig = ch11.benard_growth_rate(3.1163, Ra, Pr, all=True)`; print the top three; a convergence table for N = 16, 24, 32, 40
    (`ch11.benard_growth_rate(3.1163, Ra, Pr, N=N)`); `print(np.max(np.abs(np.imag(sig[:5]))))`. *expect:* σ₁ = 19.8666 (agrees to
    7 digits for N = 16–40); all leading σ real (imag 0) — "exchange of stabilities in action (D08)"; e-folding 179/19.87 = 8.99 s.
    *explain:* 1. Ra from (11.21); 2. the generalised eigenproblem of (11.36)–(11.38) built and solved; 3. digits gained with N.
23. `nb.check_agree` — **from scratch (curation §7):** Trefethen's `cheb` in ten lines (points, $c_i$, off-diagonal formula,
    negative-sum diagonal), `assert np.allclose(D_mine, ST.cheb(8)[0])` and the derivative of x³ on 9 points; then the (W, T̂)
    matrices for N = 24 assembled by hand with `apply_bc_rows` logic written out and `scipy.linalg.eig`, keeping finite σ with
    ∣σ∣ < 10⁴: `assert np.isclose(np.max(np.real(lam)), ch11.benard_growth_rate(3.1163, Ra, Pr, N=24), rtol=1e-8)`.
24. `nb.primer("real and imaginary parts of a complex identity", "One complex equation is two real equations: its real parts
    agree and its imaginary parts agree. Three facts do most of the work in this chapter: (i) $\\int\\lvert f\\rvert^2dz>0$ unless
    f ≡ 0 (a sum of non-negative numbers); (ii) $\\frac1{U-c}=\\frac{U-c^*}{\\lvert U-c\\rvert^2}$, so
    $\\mathrm{Im}\\frac1{U-c}=\\frac{c_i}{\\lvert U-c\\rvert^2}$ (multiply top and bottom by the conjugate, Ch. 2 P81); (iii) a
    real quantity times $c_i$ has imaginary part zero only if $c_i=0$ or the quantity vanishes.", code="import numpy as np\nU, c = 0.3, 0.1 + 0.2j                       # a real velocity and a complex wave speed\nprint(np.imag(1/(U - c)), c.imag/abs(U - c)**2)   # both 2.5: Im{1/(U-c)} = c_i/|U-c|^2")` (**P260**)
25. `nb.derivation("D08", …)` — Part F D08 (12 steps), ref "11.37" — contains **N40** (σ real for Ra > 0; exchange of
    stabilities ⇒ marginal state σ = 0 with stationary cells) and **N122** (Exercise 11.6's integrals I₁, I₂, J₁, J₂); ★★ with
    an optional `check_src` (`ch11.exchange_of_stabilities_sympy()["imag_identity"] == 0`). (Uses P218a, P260.)
26. `nb.figure` — **"The spectrum, and which eigenvalues to trust"** (two panels): (a) all finite σ for N = 24 (dots) and N = 36
    (crosses) at K = 3.1163, Ra = 3679, Pr = 7.14 on a symlog axis: the leading ones coincide, a few huge spurious ones move;
    (b) the leading eigenfunctions W(z) (blue) and T̂(z) (rose) normalised to 1. *see:* "a column of real eigenvalues; the top few
    sit on top of each other"; *read:* "trust what does not move with N; the top σ ≈ 19.9 > 0 grows, everything else decays";
    *change:* "…Ra = 1000: the top eigenvalue drops to ≈ −7.8 — every disturbance decays".
27. `nb.md` — **What would change if…** "…we only want to know *when* convection starts, not how fast? Since σ is real (D08), it
    passes through zero, not around it: setting σ = 0 removes Pr and leaves one equation for Ra(K) — C04."

#### C04 — The rigid–rigid neutral curve and the critical Rayleigh number
28. `nb.core("C04", "The rigid–rigid neutral curve Ra(K) and the critical point $\\mathrm{Ra}_c=1707.76$, $K_c=3.117$ (Fig. 11.10;
    Chandrasekhar 1961)", question="At what heating does convection start, and how wide are the first cells?")`
29. `nb.md` — **The problem in plain words:** "Turn up the stove slowly. For a while nothing moves; then, at a sharp temperature
    difference, rolls appear — and they are all about as wide as the layer is deep. Why that width, and why that threshold? Each
    possible cell width is a separate normal mode, each with its own threshold; convection starts at the cheapest one."
30. `nb.md` — **The idea** (ASCII valley):
    ```
     Ra
      ▲  \                                /      narrow cells (large K): viscosity and diffusion act over a short distance
      │   \     marginal Ra(K)          /        wide cells (small K): buoyancy drives them only weakly
      │    \___                    ___/
      │        \___  ●  ______/        ← the bottom of the valley: Ra_c ≈ 1708 at K_c ≈ 3.12  (λ_c = 2π/K_c ≈ 2d)
      └────────────────────────────────► K
     heat until Ra reaches the valley bottom: the first K to go unstable is K_c
    ```
31. `nb.derivation("D09", …)` — Part F D09 (6 steps), ref "11.40" — contains **N41** (11.39), **N42** (11.40), **N43** (11.41).
32. `nb.primer("even and odd functions; parity under z → −z", "f is *even* if f(−z) = f(z) (cos, cosh, z²) and *odd* if f(−z) =
    −f(z) (sin, sinh, z). A problem that looks the same after flipping z → −z (both walls alike) has eigenfunctions that are
    either even or odd, so we can look for them separately — even W gives one row of cells, odd W two.", code="import numpy as np\nz = np.linspace(-0.5, 0.5, 5)\nprint(np.allclose(np.cosh(2*z), np.cosh(-2*z)), np.allclose(np.sinh(2*z), -np.sinh(-2*z)))   # True True")` (**P261**)
33. `nb.note` — **N44 [B]** "**Pr drops out; even and odd modes** (our Fig. 11.9 analogue): at σ = 0 the Prandtl number has
    disappeared from (11.39)–(11.41), so the onset does not depend on the fluid's Pr — only Ra matters. Both walls are alike, so
    the eigenfunctions are even (W symmetric about z = 0: one row of cells) or odd (two rows); the even one goes unstable first."
    + `nb.figure` `ch11.benard_eigenfunction(3.1163, mode="even")` and `ch11.benard_eigenfunction(5.3647, mode="odd")`:
    streamlines (ψ, muted) over T′ (blue–rose), W(z) beside each. *see / read / change:* "one row vs two rows of counter-rotating
    rolls" / "the odd mode's W changes sign at mid-depth" / "…free walls (C05): the even mode becomes cos πz, exactly".
34. `nb.primer("cube roots of a negative number", "Every nonzero number has three cube roots. Write −1 = e^{iπ}: its cube roots
    are $e^{\\mathrm i\\pi/3},e^{\\mathrm i\\pi}=-1,e^{\\mathrm i5\\pi/3}$, i.e. −1 and ½(1 ± i√3). So $s^3=-a$ (a > 0) has
    $s=-a^{1/3}$ and $s=\\tfrac12a^{1/3}(1\\pm\\mathrm i\\sqrt3)$ — one real, two complex conjugates (Euler P45, complex powers
    P155).", code="import numpy as np\nprint(np.roots([1, 0, 0, 1]))     # s^3 + 1 = 0 -> -1 and 0.5 ± 0.866j\nprint(0.5*(1 + 1j*np.sqrt(3))**3)  # (-1+0j)")` (**P262**)
35. `nb.derivation("D10", …)` — Part F D10 (**★★★, 14 steps**), ref "11.42" — contains **N45** ((q² − K²)³ = −RaK², (11.42), q₀)
    and **N46** (the even solution and the 3 × 3 determinant) with the **slip #1 box** at step 9: "the book prints
    $B(q^{*2}-K^2)^2\\cosh q^*z$ as the third term of $(d^2/dz^2-K^2)^2W$; the correct coefficient is C (the matrix under it is
    right)"; `check_src` = Part F D10 sympy cell.
36. `nb.primer("neutral curve as the zero contour of the growth rate; critical point as its minimum", "For each wavenumber K find
    the parameter value where the leading growth rate crosses zero (`brentq`, P108): that is the neutral (marginal) curve
    Ra(K). Then minimise Ra(K) over K (`minimize_scalar`, P170): the minimum is the critical point (Ra_c, K_c). Inside the
    curve σ > 0, outside σ < 0. The same two nested searches give Ta_c (C07) and Re_c (C13).", code="from scipy.optimize import minimize_scalar\nRa = lambda K: (3.1416**2 + K**2)**3/K**2       # a neutral curve given in closed form (C05)\nr = minimize_scalar(Ra, bounds=(0.5, 6), method='bounded')\nprint(round(r.x, 3), round(r.fun, 1))            # 2.221 657.5")` (**P263**)
37. `nb.worked_example("the characteristic roots at K = 2, Ra = 1000", "1. Ra/K⁴ = 1000/16 = 62.5; s = (62.5)^{1/3} = 3.9685. 2.
    First root of (11.42): q² = −K²(s − 1) = −4 × 2.9685 = −11.874, so q = ±iq₀ with q₀ = √11.874 = 3.4459. 3. The other two:
    q² = K²[1 + ½s(1 ± i√3)] = 4[1 + 1.9843 ± 3.4369i] = 11.937 ± 13.748i, so q = 3.8822 + 1.7705i and its conjugate. 4. Check:
    (q² − K²)³ for the first root = (−11.874 − 4)³ = (−15.874)³ = −4000 = −RaK² ✓. 5. Below K⁴ = 16 the cube root s < 1 and q₀
    turns imaginary — why the determinant search starts above Ra = K⁴.")`
38. `nb.code` — `print(ch11.benard_char_roots(1000.0, 2.0))`; `for K in (2.0, 3.1163, 5.0): print(K,
    ch11.benard_marginal_Ra(K), ch11.benard_marginal_Ra_det(K))`; `crit = ch11.benard_critical()`; print with the published
    benchmark `(1707.76, 3.117)` (`reference/ch11/benchmarks.json`, Chandrasekhar 1961) and the relative difference. *expect:*
    (3.4459, 3.8822+1.7705j, 3.8822−1.7705j); Ra(2) = 2177.41 (both routes agree to 1e-8), Ra(3.1163) = 1707.76, Ra(5) = 2439.32;
    Ra_c = 1707.762 at K_c = 3.1163; difference from Chandrasekhar < 1e-5; "the book rounds to 1708 and 3.12 — same numbers";
    λ_c = 2π/K_c = 2.016 d. *explain:* 1. (11.42); 2. two independent routes to the marginal Ra (Chebyshev eigenproblem and the
    book's determinant); 3. P263's recipe.
39. `nb.check_agree` — **from scratch (curation §7):** the 3 × 3 determinant with `np.cos`/`np.cosh` of complex arguments and
    `brentq` on its imaginary part at K = 3.1163 between Ra = 1.0001K⁴ + 1 and the first sign change of a scan: `assert
    np.isclose(Ra_det, ch11.benard_marginal_Ra(3.1163), rtol=1e-8)`; and the slip: the same with the printed column (q instead of
    q* in row 3, column 3) gives a determinant that is no longer purely imaginary (print its real part).
40. `nb.figure` — **"Onset is the bottom of a valley"** (our Fig. 11.10 analogue, drawn as the book does:
    K vertical, Ra horizontal on a log axis): neutral curves rigid–rigid (purple, bold), rigid–free (purple dashed),
    free–free (muted), the odd mode (purple dotted); minima ◆ with labels (1707.76, 3.117; 1100.65, 2.682; 657.51, 2.221); the
    unstable region of the rigid curve shaded rose, labelled σ > 0 / σ < 0. *see:* "three nested valleys, rigid walls the
    highest"; *read:* "for a given Ra, the K inside the curve grow; at Ra = 2500 rigid walls let K from 1.76 to 5.08 grow";
    *change:* "…we double the depth d: Ra grows 8× (d³), so a layer twice as deep convects at 1/8 of the temperature difference".
41. `nb.plotly` — **F2** `slider_figure` over Ra (30 values log-spaced 500 … 2 × 10⁴, FAST 15): the three neutral curves (static,
    from `ch11.benard_neutral_table()`) with the unstable band of K at that Ra highlighted for each boundary pair. *explain:*
    "the band opens at the valley bottom".
42. `nb.live` — free K, Ra (log), boundaries, Pr: prints σ (`ch11.benard_growth_rate`) and the marginal Ra(K), draws W(z) —
    paired with F2 (python-viz lesson).
43. `nb.animation` — **A3** (curation §6): rolls (streamlines over T′) for K = K_c at Ra = 1.3 Ra_c (amplitude $e^{\\sigma t}$, rose
    title "grows") and Ra = 0.7 Ra_c (teal, "decays"), 60 frames (FAST 30), `player="video"`, Pr = 7.14, time in units of d²/κ.
    Notes: *see:* "the same roll pattern brightening on the left, fading on the right"; *read:* "onset is a sign change of σ at
    fixed K"; *change:* "…Pr = 0.7 (air): the growth rates change, the threshold does not (N44)".
44. `nb.note` — **N50 [B]** "**Rigid bottom, free top** (a liquid layer open to air): our computation gives Ra_c = 1100.65 at
    K_c = 2.682 (Chandrasekhar 1961 tabulates this case); one free wall lowers the threshold because a stress-free surface
    resists less." + `nb.code` `ch11.benard_critical(bc=("rigid", "free"))`. *expect:* 1100.650, 2.6823.
45. `nb.note` — **N123 [B]** "**The gravest odd mode** (Exercise 11.7): two rows of cells need Ra = 17 610.39 at K = 5.365
    (Chandrasekhar 1961, as cited by the exercise) — ten times the even threshold, so it is never the first to appear." +
    `nb.code` `ch11.benard_critical(mode="odd")`. *expect:* 17610.39 at 5.3647.
46. `nb.primer("Helmholtz equation in the plane; planforms as sums of cosines", "Linear theory fixes only ∣K∣. Any horizontal
    pattern f(x, y) with $\\nabla_H^2f=-K^2f$ has the same growth rate: rolls f = cos Kx, squares cos Kx + cos Ky, hexagons
    $\\sum_{j=1}^3\\cos(\\mathbf k_j\\cdot\\mathbf x)$ with three wave vectors of length K at 120°.", code="import numpy as np\nfrom fluidpy import ch11_instability as ch11\nx = y = np.linspace(0, 2, 41); X, Y = np.meshgrid(x, y)\nf = ch11.planform(X, Y, 3.0, 'hexagons')        # a hexagonal pattern with |K| = 3\nprint(f.shape, round(f.max(), 3))              # the maximum 3 sits at the cell centres")` (**P264**)
47. `nb.note` — **N51 [B]** "**Which pattern?** Linear theory cannot choose between rolls, squares and hexagons (P264); experiments
    near onset often show hexagons first, rolls as Ra rises, then time-dependent and finally turbulent convection at much larger
    Ra. In liquids the cell centres usually rise (viscosity falls with temperature), in gases they sink." + `nb.figure` three
    small panels of `ch11.planform` (rolls, squares, hexagons) at K = K_c. *see / read / change:* "three patterns, one ∣K∣" /
    "same growth rate for each" / "…the nonlinear terms decide (named, N52)". + climate hook: "cloud streets are rolls aligned with
    the wind; open and closed hexagonal cells cover the subtropical oceans in satellite pictures".
48. `nb.note` — **N52 [C]** "**After onset** (named): the cells grow until the nonlinear terms balance them — kinetic energy
    produced by buoyancy (warm fluid rising) equals viscous dissipation, fed by the heating. The energy bookkeeping is C14's; the
    simplest nonlinear model of it is Lorenz's (C15)."
49. `nb.explainer("benard_neutral_curve", heading="Why does convection start near Ra ≈ 1708?", why="Dragging a point in the (K, Ra)
    plane shows that every point is a whole eigenproblem — inside the curve the rolls grow, outside they decay, and the third view
    shows the determinant's imaginary part crossing zero as Ra passes the curve; switching boundaries moves the whole valley (657.5
    → 1100.65 → 1707.76).", tries=["Preset 'rigid onset': read Ra_c = 1707.76 at K = 3.117; nudge Ra up by 5 % and watch the rolls
    start growing.", "Keep Ra = 2500 and drag K from 1 to 6: the status flips twice — the unstable band.", "Switch the walls to
    free–free: the valley drops to 657.5 and moves to K = π/√2 ≈ 2.22."])`
50. `nb.md` — **What would change if…** "…both walls were stress-free? Every even derivative of W vanishes at the walls, sine
    modes solve the problem exactly, and the neutral curve becomes a formula you can minimise by hand — C05."

#### C05 — Free–free boundaries: the neutral curve (11.44) by hand
51. `nb.core("C05", "Free–free boundaries: $\\mathrm{Ra}=(n^2\\pi^2+K^2)^3/K^2$ (11.44), $K_c^2=\\pi^2/2$, $\\mathrm{Ra}_c=\\tfrac{27}
    {4}\\pi^4$", question="Can we find the onset of convection with pencil and paper?")`
52. `nb.md` — **The problem in plain words:** "Rigid plates are easy to build but hard to solve. A layer whose boundaries cannot
    hold any shear stress — oil floating on mercury, or a slice of the atmosphere with no lid — is the case Rayleigh solved in
    1916. It gives an exact formula, and its minimum, 27π⁴/4 ≈ 657.5, is the number double diffusion (C06) and Lorenz's model
    (C15) are built on."
53. `nb.md` — **The idea**: "With W = W″ = W⁗ = 0 at the walls, a sine wave sin nπ(z + ½) vanishes at both walls together with all
    its even derivatives, and $(d^2/dz^2-K^2)$ acting on it just multiplies by $-(n^2\\pi^2+K^2)$. Every operator becomes a
    number, and (11.40) becomes algebra."
54. `nb.primer("quotient rule, and differentiating with respect to K² as the variable", "$\\frac{d}{dx}\\frac{f}{g}=\\frac{f'g-fg'}{g^2}$.
    When a formula depends on K only through K², call x = K² and differentiate in x: the minimum over K and over x is the same
    point (for K > 0).", code="import sympy as sp\nx = sp.symbols('x', positive=True)            # x stands for K^2\nRa = (sp.pi**2 + x)**3/x\nprint(sp.solve(sp.diff(Ra, x), x))           # [pi**2/2]")` (**P265**)
55. `nb.derivation("D11", …)` — Part F D11 (11 steps), ref "11.44" — contains **N47** (11.43) (steps 1–5), **N48** with the **slip #2
    box** (step 7: "the book prints W = A sin(nπz); on z ∈ [−½, ½] that vanishes at both walls only for even n; the family is
    sin nπ(z + ½), whose n = 1 member is cos πz") and **N49** with the **slip #3 box** (step 9: "the book prints a factor 3 on the
    second term of dRa/dK²; with it the root equation reads K² = π² + K², which has no solution").
56. `nb.derivation("D12", …)` — Part F D12 (8 steps), ref "11.44" (ours: the free–free growth rate; the book never writes it).
57. `nb.worked_example("the free–free neutral curve with easy K", "1. K = 2: Ra = (π² + 4)³/4 = (13.870)³/4 = 667.0. 2. K = 1:
    (10.870)³/1 = 1284.2. 3. K = 3: (18.870)³/9 = 746.5. 4. The minimum is between 1 and 3: dRa/dK² = 0 at K² = π²/2 = 4.935, K =
    2.2214. 5. There Ra = (3π²/2)³/(π²/2) = (27π⁶/8)·(2/π²) = 27π⁴/4 = 657.51. 6. Growth rate at K = π/√2, Ra = 2000, Pr = 1
    (D12): a² = 3π²/2 = 14.804; (σ + a²)² = RaK²/a² = 2000 × 4.935/14.804 = 666.7; σ + a² = 25.82; σ = 11.02 (in units κ/d²).")`
58. `nb.code` — `print([ch11.benard_free_free_Ra(K) for K in (1.0, 2.0, 3.0, 4.0)])`; `print(ch11.benard_free_free_critical())`;
    `print(ch11.benard_free_free_sigma(np.pi/np.sqrt(2), 2000.0, 1.0))`; `print(ch11.benard_growth_rate(np.pi/np.sqrt(2), 2000.0,
    1.0, bc=("free", "free")))`; `print(ch11.benard_critical(bc=("free", "free")))`; `print(ch11.benard_free_free_sympy()["root"],
    ch11.benard_free_free_sympy()["printed_root"])`. *expect:* [1284.23, 667.01, 746.53, 1082.06]; Ra_c = 657.511, K_c = 2.2214;
    (11.0155, −40.6243); 11.0155 (Chebyshev, free walls — the closed form and the eigen-solver agree); 657.511, 2.2214; π²/2 and
    [] (slip #3 has no root). *explain:* 1. (11.44); 2. its minimum; 3. D12's quadratic vs the general solver.
59. `nb.check_agree` — **from scratch (curation §7):** `r = minimize_scalar(lambda K: (np.pi**2 + K**2)**3/K**2, bounds=(0.5, 6),
    method="bounded")`; `assert np.isclose(r.fun, 27*np.pi**4/4) and np.isclose(r.x, np.pi/np.sqrt(2), rtol=1e-5)`; and slip #2:
    `assert not np.allclose(ch11.benard_free_free_mode(np.array([-0.5, 0.5]), 1, printed=True), 0)`.
60. `nb.figure` — **"Free walls: a valley you can draw by hand"** (two panels): (a) (11.44) for n = 1, 2 (purple, n = 2 dashed),
    the minimum ◆ 657.51 at 2.221 and the rigid curve as a muted ghost; (b) σ(K) from D12 for Ra = 500 (teal, all negative), 657.5
    (purple, touching zero at K_c), 2000 (rose, positive for 0.754 < K < 5.358), Pr = 1. *see:* "the growth curve lifts through zero
    at one K"; *read:* "this is C01's picture with a real formula: the band of growing K opens at K_c"; *change:* "…Pr = 7: the
    band edges stay (σ = 0 does not involve Pr), the growth rates double (22.26 instead of 11.02 at Ra = 2000)".
61. `nb.md` — **What would change if…** "…the density depended on a second, slowly diffusing ingredient — salt? The same marginal
    problem reappears with Ra replaced by a difference of two Rayleigh numbers, and diffusion can now *destabilise* a column that
    is lighter on top — C06."

---

### A.5 §11.5 Double-Diffusive Instability — C06
1. `nb.section("11.5", "Double-Diffusive Instability", intro="**What is this section about?** Seawater's density depends on
   temperature and salt, and heat diffuses about a hundred times faster than salt. That difference alone can make a column that is
   lighter on top overturn — in thin 'salt fingers' or in oscillating layers. We reuse the free–free Bénard solution to find the
   criterion, and meet the subtropical ocean's thermohaline staircases.")`

#### C06 — The salt-finger criterion (11.46)
2. `nb.core("C06", "The salt-finger criterion $\\frac{gd^4}{\\nu}\\Big[\\frac\\beta{\\kappa_s}\\frac{dS}{dz}-\\frac\\alpha\\kappa\\frac{d
   \\bar T}{dz}\\Big]=657$ (11.46)", question="How can a column that is lighter on top still overturn?")`
3. `nb.md` — **The problem in plain words:** "In the subtropical Atlantic, warm salty water from the surface sits on top of cooler,
   fresher water. The warmth makes the top lighter, the salt makes it heavier, and the warmth wins: the column is stable by its
   density. Yet profilers find the water arranged in 'staircases' — uniform layers metres thick separated by sharp steps — the
   fingerprint of salt fingers. How can diffusion, which usually smooths everything, start an overturn?"
4. `nb.md` — **The idea** (ASCII, the parcel argument):
   ```
   warm, salty  ↑   a parcel pushed DOWN:   heat leaks in/out fast (κ)  → it takes the cold temperature of its new level
   ─────────────┼                           salt stays (κ_s ≈ κ/100)     → it keeps its extra salt
   cold, fresh  ↓                           ⇒ colder AND saltier than its neighbours ⇒ heavier ⇒ keeps sinking: a FINGER
   condition: the salt's destabilising effect (÷ its slow diffusivity κ_s) beats heat's stabilising one (÷ κ) by 27π⁴/4
   ```
5. `nb.note` — **N53 [B]** "**Two ingredients of density.** $\\tilde\\rho=\\rho_0[1-\\alpha(\\tilde T-T_0)+\\beta(\\tilde s-s_0)]$
   (p. 492), α thermal expansion [1/K], β haline contraction [per g/kg], both positive (Ch. 1's linear seawater EOS); κ_s salt
   diffusivity ≈ 1.5 × 10⁻⁹ m²/s vs κ ≈ 1.4 × 10⁻⁷ m²/s (τ = κ_s/κ ≈ 0.011)." + `nb.code`: `ch11.linear_eos(np.array([20.0,
   10.0]), np.array([36.5, 35.0]))` (warm salty vs cold fresh) and the parity `ch01.seawater_density_linear`-style check with the
   same coefficients. *expect:* the warm salty water is lighter by 0.88 kg/m³.
6. `nb.note` — **N54 [B]** "**Two regimes** (our Fig. 11.13 analogue, a two-row table): *fingers* — hot salty over cold fresh
   (dT̄/dz > 0, dS/dz > 0): a displaced parcel keeps its salt, onset is stationary (σ real), long thin fingers; *diffusive* —
   cold fresh over hot salty (dT̄/dz < 0, dS/dz < 0): a displaced parcel keeps its salt but loses heat, overshoots and
   oscillates with growing amplitude (σ complex, oscillatory onset, N06), ending in convecting layers separated by sharp
   interfaces — Arctic water under sea ice." + `nb.code` our free–free cubic (D13's equations with σ kept):
   `ch11.double_diffusive_sigma(np.pi**2/2, 1000.0, 2000.0, 7.0, 0.0107)` and `(…, -2e4, -1.9e6, …)`. *expect:* finger case
   real root +0.0308 (others −59.3 ± 18.0i); diffusive case 9.65 ± 70.39i with the column statically stable — "oscillatory
   growth, labelled ours (the book gives no dispersion relation)".
7. `nb.primer("uniqueness of a linear boundary-value problem", "If two unknowns obey the same linear ODE with the same right-hand
   side and the same boundary conditions, and the homogeneous problem has only the zero solution, they are equal: their
   difference solves the homogeneous problem, so it is zero. Here $(d^2/dz^2-K^2)f=0$ with f = 0 at both walls has only f = 0
   (its solutions are e^{±Kz}, which cannot vanish twice).", code="import numpy as np\nK = 2.0\nM = np.array([[np.exp(-K/2), np.exp(K/2)], [np.exp(K/2), np.exp(-K/2)]])   # f = A e^{Kz} + B e^{-Kz} at z = -1/2, +1/2\nprint(np.linalg.det(M))     # nonzero -> only A = B = 0: the solution is unique")` (**P266**)
8. `nb.derivation("D13", …)` — Part F D13 (12 steps), ref "11.46" (the book writes "repeat the derivation") — contains **N55**
   (11.45) with the **second Ra sign box** (slip #10, step 5: "§11.5 redefines $\\mathrm{Ra}\\equiv g\\alpha d^4(d\\bar T/dz)/(\\nu
   \\kappa)$, negative when heated from below — the opposite of (11.21); that is why (11.45) carries −Ra") and **N56** (T̂ =
   κ_sŝ/κ ⇒ Rs − Ra = 27π⁴/4).
9. `nb.worked_example("a finger layer with easy Rayleigh numbers", "§11.5 signs. 1. Take Ra = 1000 (warm on top: stabilising) and
   Rs = 2000 (salty on top: destabilising). 2. Margin Rs − Ra − 27π⁴/4 = 1000 − 657.5 = 342.5 > 0 ⇒ finger-unstable. 3. Is the
   column statically stable? Density gradient ∝ −αT_z + βS_z, and αT_z/(βS_z) = (Ra κ)/(Rs κ_s) = 0.5/0.0107 = 46.7 > 1: the heat
   term wins by far — lighter on top. 4. Still it overturns: the salt gradient counts κ/κ_s ≈ 93 times more in (11.46) than in the
   density.")`
10. `nb.code` — our thermocline: `r = ch11.salt_finger_unstable(dTdz=0.01, dSdz=0.002, d=0.05)`; print every entry; the thinnest
    unstable layer `d_min = (27*np.pi**4/4*1e-6/(9.81*(7.6e-4*0.002/1.5e-9 - 2e-4*0.01/1.4e-7)))**0.25`; `ch11.salt_finger_regime(...)["regime"]` for
    (0.01, 0.002, d = 0.05) fingers, (−0.01, −0.0027, d = 0.11) diffusive (Ra ≈ −20 500, Rs ≈ −1.96 × 10⁶, statically stable,
    σ = 9.04 ± 72.0i), (0.01, 0.0) stable, (−0.01, 0.0) overturning. *expect:* lhs 61 254, margin 60 597, density_stable True, R_ρ = 1.316, Ra = 875.9, Rs = 62 130;
    d_min = 1.61 cm. *explain:* 1. (11.46) with SI gradients; 2. static stability and R_ρ; 3. the four regimes of E4's presets.
11. `nb.check_agree` — **from scratch (curation §7):** (11.46) in one line, `lhs = 9.81*0.05**4/1e-6*(7.6e-4*0.002/1.5e-9 -
    2e-4*0.01/1.4e-7)`; `assert np.isclose(lhs, r["lhs"])`; and the single-component limit: `assert
    np.isclose(np.max(np.real(ch11.double_diffusive_sigma(np.pi**2/2, -2000.0, 0.0, 1.0, 1.0))), ch11.benard_free_free_sigma(np.pi/
    np.sqrt(2), 2000.0, 1.0)[0])` (verifier note i).
12. `nb.figure` — **"Stable by density, unstable by diffusion"** (regime map, our Fig. 11.13 companion): axes αdT̄/dz (x) and
    βdS/dz (y) [1/m], d = 5 cm; the static-stability line αT_z = βS_z (blue), the finger line (11.46) (amber), the diffusive
    sector; regions coloured (stable teal, fingers amber/rose, diffusive rose hatched, overturning rose); our thermocline ●.
    *see:* "a thin wedge between the two lines: lighter on top, yet fingering"; *read:* "the wedge is wide because the finger line
    is nearly the horizontal axis — salt counts κ/κ_s ≈ 93 times more"; *change:* "…κ_s = κ: the finger line rotates onto the
    static line and the wedge closes — no double diffusion without two diffusivities (F3)".
13. `nb.plotly` — **F3** `slider_figure` over κ_s/κ (20 values log-spaced 0.005 … 1): the regime map's two lines and the wedge.
14. `nb.note` — **N57 [C]** "**Fingers, staircases and layers** (named): at onset the cells are as wide as the layer, but well
    above it the fingers are long and thin; a thick finger layer breaks into a staircase of convecting layers with fingers only at
    the steps (observed in the subtropical North Atlantic and the Tyrrhenian Sea); heating a salt gradient from below builds a
    stack of diffusive layers (the Arctic under sea ice). The two requirements: two diffusivities, and opposite contributions to
    the density gradient. Climate hook: salt fingers mix heat and salt at different rates — a term ocean models parameterise
    separately (Ch. 12)."
15. `nb.explainer("salt_fingers", heading="How can a stably stratified column overturn?", why="Dragging (dT/dz, dS/dz) across the
    regime map crosses the density-stable line and the (11.46) line separately, so the reader sees a region that is stable by
    density yet finger-unstable; the parcel animation shows the parcel's temperature relaxing while its salt stays; the κ_s/κ slider
    opens and closes the wedge.", tries=["Preset 'subtropical thermocline': read 'statically stable but finger-unstable'.", "Drag
    κ_s/κ to 1: the wedge closes and the same gradients are stable.", "Preset 'Arctic': the σ-roots leave the real axis as a
    complex pair with positive real part — growing oscillations."])`
16. `nb.md` — **What would change if…** "…the 'stratification' were not of heat or salt but of angular momentum — fluid spinning
    between two cylinders? The centrifugal force plays gravity, and the same viscous threshold 1708 comes back — C07."

---

### A.6 §11.6 Centrifugal Instability: Taylor Problem — R14 R15 R16, C07
1. `nb.section("11.6", "Centrifugal Instability: Taylor Problem", intro="**What is this section about?** Fluid between two
   coaxial cylinders, the inner one spinning: above a certain speed the smooth circular flow of Ch. 8 breaks into a stack of
   doughnut-shaped vortices. Rotation acts like gravity — the centrifugal force pulls fast-spinning fluid outward — and the
   threshold has the same 1708 as Bénard convection. In the atmosphere and ocean the same mechanism is 'inertial instability'
   (Ch. 13).")`
2. `nb.recap("R14", "Rayleigh's circulation criterion", "Ch. 8 flagged circular Couette flow as inviscidly stable when
   $\\Omega_2/\\Omega_1>(R_1/R_2)^2$. In general (Rayleigh 1916): an inviscid swirling flow is unstable if the square of the
   circulation, $\\Gamma^2=(2\\pi rU_\\theta)^2$, decreases outward somewhere — $d\\Gamma^2/dr<0$ — just as a fluid is
   statically unstable when $d\\bar\\rho/dz>0$. C07 adds the reason (the ring interchange, N59).", where="Ch. 8 §8.2")`
3. `nb.recap("R15", "Axisymmetric Navier–Stokes in cylindrical coordinates", "With ∂/∂φ = 0, $\\frac{D\\tilde u_R}{Dt}-\\frac{\\tilde
   u_\\varphi^2}R=-\\frac1\\rho\\frac{\\partial\\tilde p}{\\partial R}+\\nu\\Big(\\nabla^2\\tilde u_R-\\frac{\\tilde u_R}{R^2}\\Big)$,
   $\\frac{D\\tilde u_\\varphi}{Dt}+\\frac{\\tilde u_R\\tilde u_\\varphi}R=\\nu\\Big(\\nabla^2\\tilde u_\\varphi-\\frac{\\tilde u_\\varphi}
   {R^2}\\Big)$, $\\frac{D\\tilde u_z}{Dt}=-\\frac1\\rho\\frac{\\partial\\tilde p}{\\partial z}+\\nu\\nabla^2\\tilde u_z$ and continuity
   (11.47), with $\\nabla^2=\\partial_R^2+\\frac1R\\partial_R+\\partial_z^2$ (Ch. 4's Appendix-B operators in `core.curvilinear`, P186).", where="Ch. 4
   Appendix B, Ch. 8 P186")` + **slip #12 box**: "the book prints continuity as $\\frac{\\partial}{\\partial R}(R\\tilde u_R)+
   \\frac{\\partial\\tilde u_z}{\\partial z}=0$ (in (11.47) and (11.50)); the units do not match (m/s against 1/s). The correct
   axisymmetric form is $\\frac1R\\frac{\\partial}{\\partial R}(R\\tilde u_R)+\\frac{\\partial\\tilde u_z}{\\partial z}=0$."
4. `nb.recap("R16", "Circular Couette flow", "Between cylinders R₁ < R < R₂ turning at Ω₁, Ω₂: $U_R=U_z=0$, $U_\\varphi=AR+B/R$,
   $\\frac1\\rho\\frac{dP}{dR}=\\frac{U_\\varphi^2}R$ (11.49) with $A\\equiv\\frac{\\Omega_2R_2^2-\\Omega_1R_1^2}{R_2^2-R_1^2}$,
   $B\\equiv\\frac{(\\Omega_1-\\Omega_2)R_1^2R_2^2}{R_2^2-R_1^2}$ — Ch. 8 (8.9)–(8.10), `LAM.circular_couette(…,
   return_coeffs=True)`.", where="Ch. 8 §8.2")`

#### C07 — Taylor–Couette onset: the narrow-gap equations (11.51), the Taylor number (11.52) and Ta_c (11.54)
5. `nb.core("C07", "Taylor–Couette onset: $\\mathrm{Ta}_{cr}=\\frac{1708}{\\tfrac12(1+\\Omega_2/\\Omega_1)}$ (11.54) from the
   narrow-gap equations (11.51) and the Taylor number (11.52)", question="Why does spinning the inner cylinder stack up
   vortices, while spinning the outer one does not?")`
6. `nb.md` — **The problem in plain words:** "G. I. Taylor (1923) filled the gap between two glass cylinders with water and turned
   the inner one. At low speed the water simply swirls in circles. Past a sharp speed, the flow organises into a stack of
   counter-rotating ring vortices, each about as tall as the gap is wide. Turn the outer cylinder instead and nothing happens. The
   same physics governs a journal bearing, a Couette viscometer — and, with Earth's rotation, inertial instability of jets in the
   atmosphere and ocean (Ch. 13)."
7. `nb.md` — **The idea** (ASCII):
   ```
   swap two thin fluid rings (equal mass) at r₁ < r₂, each keeping its circulation Γ = 2πrU_θ (Kelvin):
       E = Γ²/(8π²r²) per unit mass   ⇒   ΔE = (Γ₂² − Γ₁²)(1/r₁² − 1/r₂²)/(8π²)
       Γ² falls outward (Γ₂² < Γ₁²)  ⇒ ΔE < 0: energy RELEASED ⇒ unstable   (inner cylinder spinning)
       Γ² rises outward              ⇒ ΔE > 0: you must PAY         ⇒ stable     (outer cylinder spinning)
   centrifugal force ↔ gravity,  Γ² ↔ −ρ̄  — and viscosity adds a threshold, Ta_c ≈ 1708 / ((1 + Ω₂/Ω₁)/2)
   ```
8. `nb.note` — **N58 [C]** "**The set-up**: coaxial cylinders, radii R₁ < R₂, angular speeds Ω₁, Ω₂ (counter-clockwise +), gap
   d = R₂ − R₁; experiments show the first instability is axisymmetric (∂/∂φ = 0), so we only look at disturbances that are rings.
   Wavy vortices (∂/∂φ ≠ 0) come later (N69)."
9. `nb.note` — **N59 [B]** "**The ring interchange** (p. 496): with $E=U_\\theta^2/2=\\Gamma^2/(8\\pi^2r^2)$, swapping rings 1 and 2
   gives $E_{\\rm final}=\\frac1{8\\pi^2}\\big[\\frac{\\Gamma_2^2}{r_1^2}+\\frac{\\Gamma_1^2}{r_2^2}\\big]$, $E_{\\rm initial}=\\frac1{8\\pi^2}
   \\big[\\frac{\\Gamma_1^2}{r_1^2}+\\frac{\\Gamma_2^2}{r_2^2}\\big]$, so $\\Delta E=\\frac1{8\\pi^2}(\\Gamma_2^2-\\Gamma_1^2)\\Big(\\frac1{r_1^2}
   -\\frac1{r_2^2}\\Big)$ — since r₂ > r₁ the second bracket is positive and the sign of ΔE is the sign of Γ₂² − Γ₁²." +
   `nb.code`: `ch11.ring_interchange_energy(4.0, 2.0, 1.0, 2.0)`. *expect:* E_i = 0.21531, E_f = 0.10132, ΔE = −0.11399 m²/s²
   (per unit mass) — released: unstable; + `ch11.rayleigh_circulation_criterion(R, LAM.circular_couette(R, 0.1, 0.105, 0.5, 0.0))`
   (inner only: unstable everywhere) and with Ω₂ = 0.5, Ω₁ = 0.4 (outer faster: stable).
10. `nb.note` — **N61 [B]** "**Geometry** (our Fig. 11.16 analogue): `taylor_sketch(ax)` — the gap, the axis, a stack of
    counter-rotating rolls of height ≈ d." + figure notes (*see/read/change*: "a meridional cut through the gap" / "neighbouring
    rolls turn opposite ways, a pair spans one wavelength 2π/k ≈ 2d" / "…the outer cylinder turns faster than the Rayleigh line
    allows: no rolls at any speed").
11. `nb.primer("narrow-gap (small-curvature) approximation: expand in d/R and keep the leading order", "When the gap d is much
    smaller than the radius R₁, curvature terms such as $\\frac1R\\frac{d}{dR}$ or $\\frac1{R^2}$ are smaller than
    $\\frac{d^2}{dR^2}$ by factors d/R₁, $(d/R_1)^2$; we keep the leading order and drop the rest — like Ch. 8's lubrication
    (P188). With x = (R − R₁)/d ∈ [0, 1], d/dR = (1/d) d/dx.", code="R1, d = 0.10, 0.005                 # 10 cm inner radius, 5 mm gap\nprint(d/R1, (d/R1)**2)            # 0.05 and 0.0025: the sizes of the dropped terms relative to the kept one")`
    (**P267**)
12. `nb.derivation("D14", …)` — Part F D14 (**★★★, 15 steps**), ref "11.51" — contains **N60** (11.48) (step 2), **N62** (11.50)
    (steps 3–5), **N63** (normal modes and the narrow gap, steps 6, 10), **N64** (11.51) (steps 11–15); `check_src` = Part F D14
    sympy cell (`ch11.taylor_perturbation_sympy()`).
13. `nb.derivation("D15", …)` — Part F D15 (6 steps), ref "11.52" — contains **N65** (11.52) and its narrow-gap inner-only form.
14. `nb.note` — **N66 [B]** + **N67 [B]** "**Conditions and the marginal state.** No slip on both cylinders:
    $\\hat u_R=d\\hat u_R/dR=\\hat u_\\varphi=0$ at $x=0,1$ (11.53). Taylor assumed the onset is stationary (σ = 0); this is proved
    for cylinders turning the same way (μ ≥ 0) and checked numerically here, but not in general — for counter-rotation (μ < 0) the
    σ = 0 curve is only what (11.51) gives if σ = 0 is assumed." + `nb.code`: `ch11.taylor_growth_rate(3.12, 4000.0, 0.0)` and
    at μ = 0.5. *expect:* real σ > 0 (imag 0) — "unstable, stationary".
15. `nb.worked_example("water between two cylinders", "R₁ = 10 cm, R₂ = 10.5 cm (d = 5 mm, R₂/R₁ = 1.05, our choice), water ν
    = 10⁻⁶ m²/s, outer cylinder at rest, inner at Ω₁ = 0.5 rad/s. 1. (11.52): Ta = 4(Ω₁R₁² − 0)/(R₂² − R₁²) · Ω₁d⁴/ν² = 4 ×
    0.5 × 0.01/0.001025 × 0.5 × 6.25e-10/1e-12 = 6098. 2. Narrow-gap shortcut 2(Ω₁R₁d/ν)²(d/R₁) = 2 × (250)² × 0.05 = 6250 —
    2.4 % higher (the O(d/R₁) we dropped). 3. Critical value for μ = 0: (11.54) gives 1708/0.5 = 3416; the exact narrow-gap value
    is 3389.9 (0.77 % lower). 4. 6098 > 3390 ⇒ vortices. 5. Onset speed: Ta ∝ Ω₁², so Ω₁,c = 0.5 √(3389.9/6098) = 0.373 rad/s
    (3.6 revolutions per minute).")`
16. `nb.code` — `print(ch11.taylor_number(0.5, 0.0, 0.1, 0.105, 1e-6))`; `print(ch11.taylor_number_narrow_inner(0.5, 0.1, 0.005,
    1e-6))`; `for mu in (1.0, 0.5, 0.0, -0.5): print(mu, ch11.taylor_critical(mu), ch11.taylor_critical_approx(mu))`.
    *expect:* Ta = 6098 (rayleigh_stable False); 6250; μ = 1: 1707.76 at 3.1163 (= Bénard, D15) vs 1708; μ = 0.5: 2275.09 at
    3.1175 vs 2277.3; μ = 0: 3389.90 at 3.1266 vs 3416.0; μ = −0.5: 6413.7 at 3.1985 vs 6832.0 (+6.5 %: (11.54) is a co-rotation
    formula). *explain:* 1. (11.52); 2. its narrow-gap shortcut; 3. the exact narrow-gap minimum (P263's recipe) next to (11.54).
17. `nb.check_agree` — **from scratch (curation §7):** Ta from (11.52) typed out, Ta_c from (11.54) typed out, ΔE of N59 typed
    out: `assert np.isclose(Ta_mine, ch11.taylor_number(0.5, 0.0, 0.1, 0.105, 1e-6)["Ta"])`, `assert
    np.isclose(1708/(0.5*(1 + 0.5)), ch11.taylor_critical_approx(0.5))`, `assert np.isclose(dE_mine,
    ch11.ring_interchange_energy(4.0, 2.0, 1.0, 2.0)["dE"])`; and μ → 1 parity with Bénard (D15): `assert
    np.isclose(ch11.taylor_critical(1.0)["Ta_c"], ch11.benard_critical()["Ra_c"], rtol=1e-6)`.
18. `nb.figure` — **"Rotation's threshold, exact and approximate"** (two panels): (a) Ta_c(μ) exact narrow gap (purple dots from
    `ch11.taylor_critical_table()`) and (11.54) (muted dashed), μ from −0.5 to 1, the Bénard value 1707.76 marked at μ = 1, the
    error (11.54)/exact − 1 on a twin axis (rose); (b) our Fig. 11.17 analogue for R₂/R₁ = 1.05: `ch11.taylor_stability_boundary()`
    in the (Ω₂R₂²/ν, Ω₁R₂²/ν) plane with the Rayleigh line Ω₁/Ω₂ = R₂²/R₁² (blue dashed); region above the curve unstable (rose),
    Rayleigh-unstable-but-viscously-stable band hatched. *see:* "a curve that rises to the right and much more steeply to the
    left"; *read:* "co-rotation (right) needs more inner speed but becomes unstable just before the Rayleigh line; counter-rotation
    (left) is Rayleigh-unstable yet viscously held"; *change:* "…a wider gap (R₂/R₁ = 1.5): the narrow-gap theory no longer
    applies; the true curve moves (named)".
19. `nb.figure` — **"Taylor vortices"**: the marginal eigenfunction `ch11.taylor_eigenfunction(3.1266, 0.0)` as streamlines of the
    meridional flow (x across the gap, z up two wavelengths), û_φ in colour (orange). *see / read / change:* "pairs of
    counter-rotating cells" / "each cell is about as tall as the gap: λ = 2π/k_c ≈ 2.0 d" / "…μ → 1: the cells become exactly
    Bénard's rolls (D15)".
20. `nb.plotly` — **F4** `slider_figure` over μ (21 values −0.5 … 1.0): the marginal curves Ta(k) exact (purple) with their minimum
    and the constant (11.54) value (muted dashed), with the Rayleigh verdict in the trace name.
21. `nb.note` — **N68 [B]** "**Theory against Taylor's experiment** (Fig. 11.17): Taylor's measured onset speeds lay on his
    narrow-gap curve — one of the first quantitative triumphs of linear stability theory. Our panel (b) above redraws the theory
    for our radius ratio (the book's ratio is not used); the agreement shown in the book is with Taylor's own data."
22. `nb.note` — **N124 [B]** "**A second route (Exercises 11.8–11.9).** At σ = 0, (11.51) reads $(d^2/dR^2-k^2)^2\\hat u_R=(1+\\alpha
    x)\\hat u_\\varphi$ (11.92) and $(d^2/dR^2-k^2)\\hat u_\\varphi=-\\mathrm{Ta}\\,k^2\\hat u_R$; expanding $\\hat u_\\varphi=\\sum_{m=1}^
    \\infty C_m\\sin(m\\pi x)$ (11.94) and projecting (Galerkin, ch10 C07's idea with sines) gives Ta_c with four modes to 0.1 %."
    + **slip #7 box**: "the book prints (11.93) as $(d^2/dR^2-k^2)^2\\hat u_\\varphi=-\\mathrm{Ta}k^2\\hat u_R$; (11.51) at σ = 0 has
    the first power — with the square the problem is of eighth order and the threshold changes." + `nb.code`:
    `ch11.taylor_galerkin_Ta(3.1266, 0.0, n_modes=4)` vs `ch11.taylor_marginal_Ta(3.1266, 0.0)` and the `printed=True` variant.
    *expect:* agreement < 0.1 %; the printed variant differs (the cell prints by how much).
23. `nb.note` — **N69 [C]** "**After onset** (named; photographs in the book): as Ω₁ rises the vortices become wavy (∂/∂φ ≠ 0),
    then modulated, then turbulent (Coles 1965); the steady equations have many solutions for the same speeds (non-uniqueness).
    Relatives: Dean vortices in curved channels, Görtler vortices on concave walls (our sketch, Fig. 11.19 named). Pointer: the
    route to turbulence, N111 and Ch. 12."
24. `nb.explainer("taylor_couette_onset", heading="Why do stacked vortices appear between spinning cylinders?", why="Dragging (Ω₂/Ω₁,
    Ta) moves the point across the Rayleigh line and the viscous boundary separately; the profile view shows where dΓ²/dr < 0 and
    the ring-pair inspector gives ΔE for any two radii; vortices appear only past the viscous curve.", tries=["Preset 'inner only':
    raise Ta until the vortices appear — compare the onset with (11.54)'s 3416 and the exact 3390.", "Preset 'counter-rotation':
    the status says Rayleigh-unstable but viscously stable — explain why with the ring inspector near each wall.", "Drag μ to 1:
    the threshold becomes Bénard's 1708."])`
25. `nb.md` — **What would change if…** "…the 'restoring' stratification and the 'driving' shear lived in the same fluid at the
    same place — a current over a stably stratified thermocline? Then gravity and shear compete continuously across the layer, and
    one ODE (Taylor–Goldstein) decides — C08."

---

### A.7 §11.7 Instability of Continuously Stratified Parallel Flows — R17 R18 R19, C08 · R20, C09 · C10
1. `nb.section("11.7", "Instability of Continuously Stratified Parallel Flows", intro="**What is this section about?** C02 had a
   sharp interface; real oceans and atmospheres have smooth profiles of velocity U(z) and density ρ̄(z). We derive one equation
   for small waves on such a flow (Taylor–Goldstein, C08), then prove two theorems that need no solving: if the gradient
   Richardson number exceeds ¼ everywhere the flow is stable (C09), and any unstable wave speed lies inside a semicircle set by the
   slowest and fastest fluid (C10).")`
2. `nb.recap("R17", "The stratified parallel flow", "Basic state $\\tilde{\\mathbf u}=U(z)\\mathbf e_x+\\mathbf u$,
   $\\tilde p=P+p$, $\\tilde\\rho=\\bar\\rho(z)+\\rho$, inviscid Boussinesq momentum $\\frac{\\partial\\tilde{\\mathbf u}}{\\partial t}+
   (\\tilde{\\mathbf u}\\cdot\\nabla)\\tilde{\\mathbf u}=-\\frac1{\\rho_0}\\nabla\\tilde p-g\\frac{\\bar\\rho+\\rho}{\\rho_0}\\mathbf e_z$ and
   the hydrostatic balance $0=-\\frac1{\\rho_0}\\frac{\\partial P}{\\partial z}-g\\frac{\\bar\\rho}{\\rho_0}$ (p. 502–503) — Ch. 4's
   Boussinesq set and Ch. 7 §7.8's base state. Only 2-D disturbances are considered (Squire's theorem, proved in C11 for
   unstratified flow, is *assumed* here).", where="Ch. 4 §4.9, Ch. 7 §7.8")`
3. `nb.recap("R18", "The buoyancy frequency", "$N^2\\equiv-\\frac g{\\rho_0}\\frac{d\\bar\\rho}{dz}$ (7.128) — Ch. 1's (1.29) for an
   incompressible fluid; the linearised density equation becomes $\\frac{\\partial\\rho}{\\partial t}+U\\frac{\\partial\\rho}{\\partial
   x}-\\frac{\\rho_0N^2w}g=0$ (11.56).", where="Ch. 1 §1.10, Ch. 7 §7.8")` + `nb.code`: a thermocline `STRAT.brunt_vaisala_sq(1025.0,
   -0.002*1025/10.0, 0.0)` — the builder sets the arguments as `core.stratification` defines them; *expect:* N² ≈ 2 × 10⁻³ s⁻²
   for a 0.2 % density change over 10 m (the builder prints the computed number and N in rad/s).
4. `nb.recap("R19", "The stream function, §11.7 version", "$u=\\partial\\psi/\\partial z$, $w=-\\partial\\psi/\\partial x$ (11.57) —
   Ch. 7's sign (u = ∂ψ/∂z). ⚠️ slip #9: §11.8 uses u = ∂ψ/∂y, v = −∂ψ/∂x (y is the cross-stream coordinate there) and §11.14
   uses u = −∂ψ/∂z, w = ∂ψ/∂x; each block repeats its own.", where="Ch. 4 §4.3, Ch. 7 §7.8")`

#### C08 — The Taylor–Goldstein equation (11.61)
5. `nb.core("C08", "The Taylor–Goldstein equation $(U-c)\\big(\\frac{d^2}{dz^2}-k^2\\big)\\hat\\psi-\\frac{d^2U}{dz^2}\\hat\\psi+
   \\frac{N^2}{U-c}\\hat\\psi=0$ (11.61)", question="What single equation decides whether a stratified current is stable?")`
6. `nb.md` — **The problem in plain words:** "A tidal current flows over a stably stratified thermocline; a jet stream blows over a
   cold inversion. The velocity changes smoothly with height, and so does the density. Shear wants to roll the layers up; the
   density stratification resists lifting heavy water. We want one equation that contains both, for one wave at a time."
7. `nb.md` — **The idea**: "Exactly C02's recipe with smooth profiles: linearise, use a stream function so continuity is automatic,
   insert $e^{\\mathrm ik(x-ct)}$, eliminate pressure and density. One second-order ODE in z remains, with c as the eigenvalue: U″
   is the shear profile's curvature (C12 will make it the hero), N²/(U − c) the stratification. With N² = 0 it is Rayleigh's
   equation (11.81)."
8. `nb.note` — **N70 [C]** "**History.** Taylor (1915) conjectured from examples that a gradient Richardson number below ¼ is
   needed for instability; Prandtl, Goldstein, Richardson, Synge and Chandrasekhar suggested other values; Miles (1961) proved
   Taylor's conjecture and Howard (1961) gave the short proof we follow (C09)."
9. `nb.derivation("D16", …)` — Part F D16 (7 steps), ref "11.57" — contains **N71** (11.55) with the **slip #4 box** at step 5
   ("the book's w-equation in (11.55) has $-\\frac1{\\rho_0}\\frac{\\partial p}{\\partial x}$; the vertical momentum equation needs
   $-\\frac1{\\rho_0}\\frac{\\partial p}{\\partial z}$, as (11.57) then uses").
10. `nb.derivation("D17", …)` — Part F D17 (9 steps), ref "11.61" — contains **N72** (11.58)–(11.60) and (11.61).
11. `nb.note` — **N73 [B]** "**Pairs again.** (11.61) contains no i, so if (ψ̂, c) solves it so does (ψ̂*, c*): every growing
    mode has a decaying twin, and any c_i ≠ 0 means instability (as N17)." + **N74 [B]** "**Rigid lids:** w = 0 at z = 0 and d
    means $\\partial\\psi/\\partial x=\\mathrm ik\\hat\\psi e^{\\mathrm ik(x-ct)}=0$, so $\\hat\\psi(0)=\\hat\\psi(d)=0$ (11.62); for an
    unbounded layer, ψ̂ → 0 far away."
12. `nb.primer("quadratic eigenvalue problem c²M₂ + cM₁ + M₀ and its companion linearisation", "Multiplying (11.61) by (U − c)
    leaves c² in it, so collocation gives $(c^2M_2+cM_1+M_0)\\mathbf v=0$. Introduce $\\mathbf w=c\\mathbf v$; then
    $\\begin{pmatrix}0&I\\\\-M_0&-M_1\\end{pmatrix}\\begin{pmatrix}\\mathbf v\\\\\\mathbf w\\end{pmatrix}=c\\begin{pmatrix}I&0\\\\0&M_2
    \\end{pmatrix}\\begin{pmatrix}\\mathbf v\\\\\\mathbf w\\end{pmatrix}$ — an ordinary generalised eigenproblem twice as large
    (P259). For 1 × 1 'matrices' it is the quadratic formula in disguise.", code="import numpy as np\nfrom scipy.linalg import eig\nM2, M1, M0 = 1.0, -3.0, 4.0                     # c^2 - 3c + 4 = 0 (the KH tiny example's quadratic)\nA = np.array([[0, 1], [-M0, -M1]]); B = np.array([[1, 0], [0, M2]])\nprint(eig(A, B, right=False))                   # 1.5 ± 1.3229j")` (**P268**)
13. `nb.primer("singular point of an ODE", "Where the coefficient of the highest derivative vanishes, an ODE is *singular*: in
    (11.61) and (11.81) that happens where U(z) = c. For a growing mode (c_i ≠ 0) it never happens; for a neutral mode (real c
    inside the velocity range) the solution may have a kink or a log there — the critical layer of C12. Viscosity (C11) removes
    the singularity by restoring the fourth derivative.", code="import numpy as np\nz = np.linspace(-1, 1, 5); U = z; c = 0.5          # Couette U = z, a real c inside the range\nprint(U - c)                                      # changes sign: the coefficient (U - c) vanishes at z = 0.5")` (**P269**)
14. `nb.primer("mapping an infinite domain to [−1, 1]", "A tanh shear layer or a jet extends to z = ±∞. The algebraic map
    $z=a\\xi/\\sqrt{1-\\xi^2}$ sends ξ ∈ (−1, 1) to all z; Chebyshev points in ξ crowd near the layer (|z| ≲ a) and still reach
    far away. Derivatives follow from the chain rule (P49): $\\frac{d}{dz}=\\frac{d\\xi}{dz}\\frac{d}{d\\xi}$. Always check that the
    answer does not change when a or N changes.", code="import numpy as np\nxi = np.cos(np.pi*np.arange(9)/8)[1:-1]          # interior Chebyshev points\nprint(np.round(3*xi/np.sqrt(1 - xi**2), 2))       # z values: dense near 0, out to about +-15")` (**P270**)
15. `nb.worked_example("two small cases by hand", "**(a) Couette, no stratification:** U = z on [0, 1], U″ = 0, N² = 0. Then
    (11.61) is (U − c)(ψ̂″ − k²ψ̂) = 0. If c is not between 0 and 1, U − c ≠ 0, so ψ̂″ = k²ψ̂ with ψ̂(0) = ψ̂(1) = 0, whose only
    solution is ψ̂ = 0 (P266): no eigenvalue outside the velocity range, no instability. **(b) the companion trick:** for the KH
    quadratic c² − 3c + 4 = 0 (C02's tiny example with U₁ = 6), the 2 × 2 companion matrix [[0, 1], [−4, 3]] has eigenvalues
    1.5 ± 1.3229i — the same roots (P268).")`
16. `nb.code` — `prof = ch11.richardson_profiles("tanh", J=0.1)`; `c = ch11.taylor_goldstein_eigs(0.4, prof["U"], prof["Upp"],
    prof["N2"], bc="decay", map_scale=3.0, N=100)`; print the leading c and check that its conjugate is in the spectrum; then N² = 0:
    `cR = ch11.rayleigh_eigs(0.4, np.tanh, lambda z: -2*np.tanh(z)/np.cosh(z)**2, bc="decay", N=100)` vs
    `ch11.taylor_goldstein_eigs(0.4, …, lambda z: 0*z, …)`. *expect:* leading c ≈ 0 + 0.311i (kc_i ≈ 0.124, approximate); its
    conjugate present to 1e-15; at J = 0 the two solvers agree (c ≈ 0 + 0.470i at k = 0.4, kc_i = 0.188). *explain:* 1. the tanh
    shear layer with N² = J sech²z (C09's family); 2. (11.61) on the mapped infinite domain (P270, P268); 3. N² = 0 ⇒ Rayleigh.
17. `nb.check_agree` — **sympy instead of a second numeric version (curation §7):** `out = ch11.stratified_shear_sympy()`;
    `assert out["tg_residual"] == 0`; and slip #4: `assert ch11.stratified_shear_sympy(printed=True)["printed_slip4_reaches_tg"] is
    False` — "with ∂p/∂x in the w-equation, eliminating p does not give (11.61)". The cell shows `out` line by line.
18. `nb.figure` — **"A Taylor–Goldstein spectrum"** (c-plane): all converged eigenvalues for k = 0.4, J = 0.1 (N = 100 and 140
    overlaid): the unstable mode (rose) and its conjugate (teal), the dense string of real c in [−1, 1] (muted — the continuous
    spectrum: critical-layer modes, not instabilities), the semicircle of C10 as a muted ghost. *see:* "two mirror dots off the axis
    and a line of dots on it"; *read:* "growing and decaying twins; the line is where U(z) = c somewhere (P269)"; *change:* "…J =
    0.3: the pair falls onto the axis — no growth (C09 proves it must)".
19. `nb.md` — **What would change if…** "…we wanted a *guarantee* of stability without computing any spectrum? Multiply the
    equation by the conjugate, integrate, and look at the imaginary part — C09 turns (11.61) into the Richardson-number theorem."

#### C09 — The Miles–Howard criterion: Ri > ¼ everywhere ⇒ stable (11.67)
20. `nb.recap("R20", "The gradient Richardson number", "$\\mathrm{Ri}(z)\\equiv N^2/(dU/dz)^2$ (11.66): stratification's
    restoring strength over the shear's — Ch. 4's Ri_g and Ch. 1's N² together. Large Ri: buoyancy wins; small Ri: shear wins.",
    where="Ch. 4 §4.11 (notation Ri_g), Ch. 1 §1.10")`
21. `nb.core("C09", "The Miles–Howard criterion $\\mathrm{Ri}>\\tfrac14$ everywhere ⇒ stable (11.67)", question="How much
    stratification guarantees that a shear layer cannot break into billows?")`
22. `nb.md` — **The problem in plain words:** "Every ocean and climate model has to decide, grid box by grid box, whether the shear
    is strong enough to mix the water or air. The rule they use is a number: if the gradient Richardson number is above about ¼,
    no mixing from shear instability. Where does ¼ come from — and is it a guarantee, a trigger, or both?"
23. `nb.md` — **The idea** (ASCII):
    ```
    TG equation  ──substitute φ = ψ̂/(U − c)^{1/2}──►  self-adjoint form (11.64)
                 ──× φ*, integrate, by parts──────►  (11.65): ∫ (N² − ¼U′²)/(U − c) |φ|² = ∫(U − c)(|φ′|² + k²|φ|²) + ∫ ½U″|φ|²
                 ──imaginary part──────────────────►  c_i · ∫ (N² − ¼U′²)/|U − c|² |φ|²  =  − c_i · ∫ (|φ′|² + k²|φ|²)
    if N² > ¼U′² everywhere: c_i × (positive) = c_i × (negative)  ⇒  c_i = 0   ⇒  Ri > ¼ everywhere guarantees stability
    ```
24. `nb.primer("integrating on a Chebyshev grid (Clenshaw–Curtis weights)", "On Chebyshev points, the integral of the polynomial
    through the samples is a weighted sum $\\int f\\,dx\\approx\\sum_jw_jf(x_j)$; the weights are exact for polynomials up to
    degree N, so smooth integrands converge as fast as the collocation itself. We use them to check the integral identities
    (11.65), (11.69)–(11.70) and the energy budget (11.88) on computed modes (Simpson, P203, would be far less accurate on these
    clustered points).", code="import numpy as np\nfrom fluidpy.core import stability as ST\nD, x = ST.cheb(8); w = ST.clenshaw_curtis_weights(8)\nprint(w @ x**2, w @ np.cos(x))                   # 0.666667 (= 2/3) and 1.682942 (= 2 sin 1)")` (**P271**)
25. `nb.derivation("D18", …)` — Part F D18 (**★★★, 14 steps**), ref "11.67" — contains **N75** (11.63) and the two derivative
    formulas (steps 2–4), **N76** (11.64) (steps 5–8, the book's 'after some rearrangement' written out), **N77** (11.65) and its
    imaginary part (steps 9–12); `check_src` = Part F D18 sympy cell.
26. `nb.worked_example("Ri for a thermocline and for the tanh layer", "**Thermocline:** a current changing by ΔU = 0.1 m/s over
    h = 2 m, so dU/dz ≈ 0.05 s⁻¹, with N² = 2 × 10⁻⁴ s⁻². Ri = 2e-4/0.05² = 2e-4/0.0025 = 0.08 < ¼: instability is *allowed* (not
    guaranteed). **tanh layer:** U = tanh z, N² = J sech²z. U′ = sech²z, so Ri(z) = J sech²z/sech⁴z = J cosh²z: smallest at the
    centre, Ri_min = J. J = 0.1: Ri_min = 0.1 < ¼ (the computed growth is positive, C08); J = 0.3: Ri > ¼ everywhere ⇒ stable by
    (11.67), whatever k.")`
27. `nb.code` — `z = np.linspace(-4, 4, 801)`; `for J in (0.1, 0.24, 0.3): p = ch11.richardson_profiles("tanh", J=J); print(J,
    ch11.miles_howard_stable(z, U=p["U"], N2=p["N2"]))`; `print([ch11.tg_growth(0.4, J) for J in (0.0, 0.1, 0.3)])`; identity:
    `c, vec = ch11.taylor_goldstein_eigs(0.4, …J = 0.1…, return_vectors=True)` and `ch11.richardson_identity_check(0.4, c, vec["psi"],
    vec["z"], …)`. *expect:* J = 0.1: guaranteed_stable False, Ri_min = 0.1 at z = 0; 0.24: False; 0.3: True; growth 0.188, ≈ 0.124,
    0.0; identity residual < 1e-8, imaginary parts equal. *explain:* 1. (11.66) on the profile; 2. the theorem's verdict; 3. the
    computed growth agrees; 4. (11.65) holds for the computed mode (P271).
28. `nb.check_agree` — **from scratch (curation §7):** `Ri_mine = J/np.cosh(z)**2/(1/np.cosh(z)**2)**2`; `i = np.argmin(Ri_mine)`;
    `assert np.isclose(Ri_mine[i], ch11.miles_howard_stable(z, U=p["U"], N2=p["N2"])["Ri_min"]) and z[i] == 0.0`.
29. `nb.figure` — **"Ri < ¼ allows billows; Ri > ¼ forbids them"** (two panels): (a) U(z) (orange), N²(z) (blue), Ri(z) (purple) with
    the ¼ line for J = 0.1 and 0.3; (b) the growth map kc_i(k, J) (cached `ch11.tg_growth_map`, contour, rose) with the J = ¼ line
    (purple) and Drazin's analytic neutral curve for this family overlaid only if the builder's computed map supports it (else the
    zero contour of the map). *see:* "a growth tongue that closes at J = ¼"; *read:* "above the line nothing grows — the guarantee;
    below it growth is allowed and here actually happens"; *change:* "…make the density layer thinner than the shear layer (R > 1):
    Ri_min sits off-centre and the tongue changes shape (E6's chips)".
30. `nb.note` — **N78 [B]** "**Necessary, not sufficient.** Ri < ¼ somewhere is needed for instability but does not force it; there
    is no universal critical Ri (walls and profile shapes lower it). For common shear layers (linear, tanh, erf) instability does
    start when Ri_min drops below ¼, and as Ri_min → 0 the fastest wave has kL ≈ 0.445 for U = U₀ tanh(z/L), i.e. λ ≈ 14.1 L
    (Michalke 1964; our solver 0.4449, 0.1897). Lab (Scotti & Corcos 1972) and ocean (Eriksen 1978) data support Ri_min < ¼ as a
    useful guide."
31. `nb.plotly` — **F5** `slider_figure` over J (26 values 0 … 0.3, FAST 13): Ri(z) with the ¼ line, and kc_i(k) from the cached map.
32. `nb.explainer("richardson_shear_instability", heading="Why does Ri = ¼ decide whether a shear layer billows?", why="Sliding the
    stratification J shows the Ri(z) profile dip below the ¼ line at the same moment the growth tongue of the (k, J) map is reached
    and the billow animation starts growing — the necessary-not-sufficient logic is visible.", tries=["Preset 'J = 0.24' vs 'J =
    0.26': read the status before and after ¼.", "Drag k at J = 0.1: growth only in a band of wavenumbers — Ri < ¼ allows,
    the wavelength decides.", "Click the Ri(z) profile: the inspector computes N²/U′² at that height."])`
33. `nb.md` — **What would change if…** "…we asked not *whether* a wave can grow but *how fast and how fast it travels*? A second
    substitution, F = ψ̂/(U − c), bounds both at once — C10."

#### C10 — Howard's semicircle
34. `nb.core("C10", "Howard's semicircle $\\big[c_r-\\tfrac12(U_{\\max}+U_{\\min})\\big]^2+c_i^2\\le\\big[\\tfrac12(U_{\\max}-U_{\\min})
    \\big]^2$ (p. 507, from $\\int[U^2-c_r^2-c_i^2]Q\\,dz>0$ (11.72)) and $kc_i<\\tfrac k2(U_{\\max}-U_{\\min})$", question="Where in
    the complex plane can the wave speed of an unstable mode be?")`
35. `nb.md` — **The problem in plain words:** "Before computing a spectrum it helps to know where to look — and how fast anything
    could possibly grow. In a current whose speed varies between 0.2 and 1 m/s, can an unstable wave travel at 3 m/s? Can it grow
    in a microsecond? Howard's answer: no — every unstable wave speed lies in a half-disc spanned by the slowest and fastest fluid."
36. `nb.md` — **The idea** (ASCII):
    ```
     c_i ▲            ___
         │        .-'     '-.          every unstable c lies under this arc
         │      /             \        radius = ½(U_max − U_min),  centre = ½(U_max + U_min)
         │     |       ●       |       ⇒ U_min < c_r < U_max   and   k c_i ≤ (k/2)(U_max − U_min)
         └─────┴───────────────┴────► c_r
             U_min            U_max
    ```
37. `nb.primer("completing the square into a circle (x − a)² + y² ≤ R²", "An inequality x² + y² − 2ax + b ≤ 0 is a disc: add and
    subtract a², $(x-a)^2+y^2\\le a^2-b$, centre (a, 0), radius $\\sqrt{a^2-b}$ (P129 did this for tensors). Howard's last step
    is exactly this with x = c_r, y = c_i.", code="import numpy as np\na, b = 1.0, 0.0                                  # x^2 + y^2 - 2x <= 0\nprint(np.sqrt(a**2 - b))                          # radius 1 around (1, 0): the half-disc on [0, 2]")` (**P272**)
38. `nb.derivation("D19", …)` — Part F D19 (**★★★, 14 steps**), ref "11.72" — contains **N79** (11.68) and its derivatives,
    **N80** the divergence form with the **slip #11 box** (step 5: "the book prints $-k^2(U-c)F$; expanding the line above it gives
    $-k^2(U-c)^2F$, and the next integral (with $(U-c)^2\\lvert F\\rvert^2$) confirms it"), **N81** (11.69)–(11.70), **N82** (11.71),
    and the **slip #8 box** at step 12 ("the weight Q ≥ 0 must be inside the integral from the start, and the inequality is ≤");
    `check_src` = Part F D19 sympy cell.
39. `nb.worked_example("a semicircle with easy numbers", "U between U_min = 0 and U_max = 2 m/s. 1. Centre ½(0 + 2) = 1, radius
    ½(2 − 0) = 1. 2. Is c = 1.2 + 0.5i possible? (1.2 − 1)² + 0.5² = 0.04 + 0.25 = 0.29 ≤ 1 ✓. 3. c = 2.5 + 0.1i? (1.5)² + 0.01 =
    2.26 > 1 ✗ — no unstable wave can outrun the fastest fluid. 4. Growth bound for k = 2 m⁻¹: kc_i ≤ (2/2)(2 − 0) = 2 s⁻¹, so the
    e-folding time is at least 0.5 s. 5. tanh layer (U from −1 to 1): its fastest mode c = 0 + 0.4267i is well inside the unit
    half-disc.")`
40. `nb.code` — `cr, ci = ST.howard_semicircle(-1.0, 1.0)`; collect the unstable eigenvalues of the tanh layer for k = 0.1 … 0.9 and
    J = 0, 0.1 (`ch11.taylor_goldstein_eigs`); `print(all(ST.in_howard_semicircle(c, -1, 1) for c in unstable))`;
    `ch11.howard_identity_check(…)` for one mode. *expect:* True; res_69, res_70 < 1e-8; cr_mean = c_r (≈ 0 for the symmetric
    layer). *explain:* 1. the arc; 2. every computed unstable c tested; 3. (11.69)–(11.70) hold for a computed mode (P271).
41. `nb.check_agree` — **from scratch (curation §7):** `inside = (np.real(cs) - 0.0)**2 + np.imag(cs)**2 <= 1.0 + 1e-9`; `assert
    np.all(inside) and np.all(inside == np.array([ST.in_howard_semicircle(c, -1.0, 1.0) for c in cs]))`.
42. `nb.figure` — **"Every unstable wave speed sits inside the half-disc"** (our Fig. 11.20 analogue): the semicircle (purple), the
    velocity range on the c_r axis (orange), the computed unstable eigenvalues of the tanh layer for several k and J (rose dots),
    the growth bound kc_i ≤ k(U_max − U_min)/2 as a dashed line on a second panel of kc_i(k). *see:* "dots on the imaginary axis,
    under the arc"; *read:* "the symmetric layer's waves do not travel (c_r = 0, N19); the arc is far from tight — a sanity
    window, not a prediction"; *change:* "…an asymmetric layer U = 1 + tanh z: the arc and the dots shift right together by 1
    (Galilean shift)".
43. `nb.plotly` — **F6** `slider_figure` over the sin-profile half-width b (16 values 1.5 … 3.0): Rayleigh eigenvalues of
    U = sin y on ∣y∣ ≤ b on the c-plane with the semicircle on [−sin b, sin b] (cached `ch11.rayleigh_spectrum_table`); the
    unstable dots appear only when 2b > π (N95). *explain:* "the semicircle confines them; growth vanishes as 2b → π".
44. `nb.md` — **What would change if…** "…viscosity came back? The equation gains a fourth derivative, the critical-layer
    singularity disappears, conjugate pairs break — and, surprisingly, a profile with no inflection point can become unstable. That
    is the Orr–Sommerfeld story, C11–C14."

---

### A.8 §11.8 Squire's Theorem and the Orr–Sommerfeld Equation — R21 R22, C11
1. `nb.section("11.8", "Squire's Theorem and the Orr-Sommerfeld Equation", intro="**What is this section about?** Viscous
   parallel flows — water in a channel, air in a boundary layer. Two simplifications make them tractable: Squire's theorem says
   two-dimensional disturbances go unstable first, and a stream function then turns the linearised Navier–Stokes equations into one
   fourth-order equation, the Orr–Sommerfeld equation, whose eigenvalue c(k, Re) decides stability.")`
2. `nb.recap("R21", "Dimensionless Navier–Stokes", "Lengths by L, velocities by U₀, time by L/U₀, pressure by ρU₀², Re = U₀L/ν
   (Ch. 4 §4.11, D30): the perturbed x-momentum equation is $\\frac{\\partial u}{\\partial t}+(U+u)\\frac{\\partial}{\\partial x}(U+u)+
   v\\frac{\\partial}{\\partial y}(U+u)=-\\frac{\\partial}{\\partial x}(P+p)+\\frac1{\\mathrm{Re}}\\nabla^2(U+u)$ (11.73) — the starting
   line of D20. ⚠️ fourth scaling of the chapter (convention 6).", where="Ch. 4 §4.11")`
3. `nb.recap("R22", "The stream function, §11.8 version", "Here y is across the flow: $u=\\partial\\psi/\\partial y$,
   $v=-\\partial\\psi/\\partial x$ (Ch. 4 (4.12)). With $\\psi=\\phi(y)e^{\\mathrm ik(x-ct)}$: $\\hat u=\\phi'$, $\\hat v=-\\mathrm
   ik\\phi$ — and φ is now the stream-function amplitude, **not** a potential (slip #9).", where="Ch. 4 §4.3")`

#### C11 — The Orr–Sommerfeld equation (11.79)
4. `nb.core("C11", "The Orr–Sommerfeld equation $(U-c)\\big(\\frac{d^2\\phi}{dy^2}-k^2\\phi\\big)-\\frac{d^2U}{dy^2}\\phi=\\frac1{\\mathrm
   ik\\mathrm{Re}}\\big(\\frac{d^4\\phi}{dy^4}-2k^2\\frac{d^2\\phi}{dy^2}+k^4\\phi\\big)$ (11.79)", question="What equation decides
   whether a viscous channel or boundary-layer flow is stable?")`
5. `nb.md` — **The problem in plain words:** "Reynolds showed in 1883 that water in a pipe stays smooth at low speed and turns
   turbulent at high speed. The wing of an airliner keeps a laminar boundary layer only near its nose. To predict where smooth flow
   breaks down we need the viscous version of C08's equation: small waves on U(y), now with friction. Friction should only damp
   — but C13 will show that it can do the opposite."
6. `nb.md` — **The idea** (ASCII):
   ```
   3-D wave e^{i(kx + mz − kct)} at Re   ──Squire──►  2-D wave with k̄ = √(k² + m²) at the LOWER Re̅ = k Re / k̄
        ⇒ the first instability (lowest Re) is two-dimensional
   2-D: u = ∂ψ/∂y, v = −∂ψ/∂x, ψ = φ(y)e^{ik(x−ct)}   ──eliminate p──►   one 4th-order ODE for φ:  Orr–Sommerfeld
        left side: Rayleigh's inviscid operator      right side: viscosity, ∝ 1/(ikRe)
   ```
7. `nb.note` — **N83 [C]** "**Viscosity can destabilise.** In Bénard and Taylor flows viscosity only raised thresholds; in channels
   and boundary layers it can be the very cause of instability (C13), for a reason C14 explains."
8. `nb.derivation("D20", …)` — Part F D20 (7 steps), ref "11.77" — contains **N84** (11.74), **N85** (11.75), **N86** (11.76),
   **N87** (11.77).
9. `nb.note` — **N88 [B]** "**Squire's theorem.** The transformation $\\bar k=\\sqrt{k^2+m^2}$, $\\bar c=c$, $\\bar k\\bar u=k\\hat
   u+m\\hat w$, $\\bar v=\\hat v$, $\\bar p/\\bar k=\\hat p/k$, $\\bar k\\overline{\\mathrm{Re}}=k\\mathrm{Re}$ (11.78) turns the
   three-dimensional equations (11.77) into the two-dimensional ones with m = ŵ = 0 (multiply the x-equation by k, the z-equation by
   m, add, divide by k̄ — a-D24, stated). Because k̄ ≥ k, the equivalent 2-D problem has the lower Reynolds number $\\overline{
   \\mathrm{Re}}=k\\mathrm{Re}/\\bar k\\le\\mathrm{Re}$: the critical Re is found with 2-D waves. ⚠️ Not valid with rotation or
   stratification (Ch. 13); §11.7 assumed it." + `nb.code`: `ST.squire_transform(1.0, np.tan(np.pi/6), 1e4)` and the numerical
   check `ch11.os_3d_eigs(1.0, 0.5774, 1e4, …Poiseuille…, N=60)` vs `ch11.orr_sommerfeld_eigs(1.1547, 8660.25, …)`. *expect:*
   kbar = 1.1547, Rebar = 8660.25, angle 30°; leading c equal to 1e-8.
10. `nb.note` — **N89 [B]** "**What Squire means physically.** An oblique wave sees only the component of the basic flow along its
    crests' normal — a slower flow, a lower effective Reynolds number — and its growth rate $\\bar k\\bar c_i$ exceeds $kc_i$. So
    two-dimensional disturbances are the most dangerous (for parallel, non-rotating, unstratified flows)."
11. `nb.derivation("D21", …)` — Part F D21 (10 steps), ref "11.79" — the book writes only 'this effort yields'; ★★ with an
    optional `check_src` (`ch11.os_derivation_sympy()["residual"] == 0` and the planted v̂ = +ikφ ≠ 0).
12. `nb.note` — **N90 [B]** "**No slip:** $\\phi=d\\phi/dy=0$ at $y=y_1$ and $y_2$ (11.80) — v̂ = −ikφ = 0 and û = φ′ = 0 on each
    wall. Four conditions for a fourth-order equation; the eigenvalue is c."
13. `nb.worked_example("a channel and an oblique wave", "1. Water (ν = 10⁻⁶ m²/s) in a channel of half-width L = 1 cm with
    centreline speed U₀ = 0.58 m/s: Re = U₀L/ν = 0.58 × 0.01/1e-6 = 5800 — right at the edge of instability (C13: 5772). 2. An
    oblique wave at 30° to the flow, k = 1, m = tan 30° = 0.577 (in units of 1/L): k̄ = √(1 + 0.333) = 1.1547. 3. Its 2-D twin
    sees Re̅ = Re · k/k̄ = 5800 × 0.866 = 5023 < 5772: the oblique wave is stable even though the 2-D wave at the same Re is
    (just) not.")`
14. `nb.code` — `U = lambda y: 1 - y**2; Upp = lambda y: -2 + 0*y`; `c = ch11.orr_sommerfeld_eigs(1.0, 1e4, U, Upp, N=100)`;
    `print(c[0])`; convergence `for N in (60, 80, 100, 120)`; the benchmark from `reference/ch11/benchmarks.json` (Orszag 1971:
    0.23752649 + 0.00373967i). *expect:* c₁ = 0.23752649 + 0.00373967i (all four N agree to 1e-8); c_i > 0: this wave grows,
    kc_i = 0.00374 per L/U₀. *explain:* 1. the plane-Poiseuille base flow in channel units (half-width, centreline speed);
    2. (11.79)–(11.80) by Chebyshev collocation with clamped rows; 3. convergence; 4. the published benchmark.
15. `nb.check_agree` — **from scratch (curation §7):** assemble the OS matrices for Poiseuille from `ST.cheb(100)` by hand
    ($A=U(D^2-k^2)-U''-\\frac1{\\mathrm ik\\mathrm{Re}}(D^4-2k^2D^2+k^4)$, $B=D^2-k^2$, rows 0, 1, N − 1, N replaced by φ = 0, φ′ = 0),
    `scipy.linalg.eig`, keep ∣c∣ < 2, take the largest c_i: `assert np.isclose(c_mine, c[0], atol=1e-8) and np.isclose(c_mine,
    0.23752649 + 0.00373967j, atol=1e-7)`.
16. `nb.figure` — **"The Orr–Sommerfeld spectrum is a Y"** (c-plane, Re = 10⁴, k = 1, N = 100 dots and N = 140 crosses): the three
    branches (wall modes A, centre modes P, the S branch descending), the one unstable mode (rose, c_i = 0.0037) at the top of the A
    branch, spurious modes that move with N (muted). *see:* "a Y of eigenvalues, one just above the axis"; *read:* "that weak mode
    is the Tollmien–Schlichting wave; everything else decays"; *change:* "…Re = 5000: the TS mode drops below the axis (c_i =
    −0.0018) — C13 finds where it crosses".
17. `nb.md` — **What would change if…** "…Re → ∞? The right side of (11.79) disappears, the order drops from four to two, two
    boundary conditions must be abandoned, and the equation is Rayleigh's — whose theorems need no computer at all (C12)."

---

### A.9 §11.9 Inviscid Stability of Parallel Flows — C12
1. `nb.section("11.9", "Inviscid Stability of Parallel Flows", intro="**What is this section about?** Drop viscosity from the
   disturbances (the basic profile may still be a viscous one): the Orr–Sommerfeld equation becomes Rayleigh's equation. Two
   theorems then classify profiles by their shape alone — a growing wave needs an inflection point (Rayleigh) at which the vorticity
   peaks (Fjørtoft) — and neutral waves come with a 'critical layer' and Kelvin's cat's-eye streamlines.")`

#### C12 — Rayleigh's inflection-point theorem (11.83)–(11.84) and Fjørtoft's sharpening (11.85)–(11.86)
2. `nb.core("C12", "Rayleigh's inflection-point theorem $c_i\\int\\frac{U''}{\\lvert U-c\\rvert^2}\\lvert\\phi\\rvert^2dy=0$ (11.84):
   an unstable inviscid parallel flow needs $U''$ to change sign", question="Can we tell from the shape of a velocity profile
   whether it can be unstable without viscosity?")`
3. `nb.md` — **The problem in plain words:** "Jets, wakes and mixing layers break into eddies within a few widths downstream;
   channel and boundary-layer flows hold out much longer. A forecaster sees jet streams meander into weather systems; an engineer
   sees a wake shed vortices. Is there something about the *shape* of a profile that makes it fragile?"
4. `nb.md` — **The idea**: "U″ is the slope of the vorticity −U′ across the flow. An inflection point (U″ = 0 with a sign change)
   is where the vorticity has an extremum — the smooth version of C02's vortex sheet. Rayleigh's identity says: no such point, no
   inviscid growth. Fjørtoft adds that the extremum must be a *maximum* of the vorticity magnitude."
5. `nb.note` — **N91 [B]** "**Rayleigh's equation** (Re → ∞ in (11.79)): $(U-c)\\big(\\frac{d^2\\phi}{dy^2}-k^2\\phi\\big)-\\frac{d^2U}
   {dy^2}\\phi=0$ (11.81) — C08's (11.61) with N² = 0. The limit is *singular*: the order drops from 4 to 2 and only two boundary
   conditions can be kept." + **N92 [B]** "$\\phi=0$ at $y=y_1,y_2$ (11.82) (no flow through the walls; slip is now allowed). No i
   appears, so eigenvalues come in conjugate pairs — a property the viscous term (with its i) destroys." + `nb.figure`: the
   Rayleigh and Orr–Sommerfeld spectra of the tanh layer side by side (k = 0.4; OS at Re = 100 and 1000 on the mapped domain):
   mirror pairs vs no pairs. *see / read / change:* "symmetric dots vs a lopsided cloud" / "viscosity breaks the c ↔ c* symmetry"
   / "…Re → ∞: the OS unstable mode approaches the Rayleigh one".
6. `nb.derivation("D22", …)` — Part F D22 (9 steps), ref "11.84" (uses P218a, P260, P255).
7. `nb.note` — **N93 [B]** "**Fjørtoft's theorem.** The real part of (11.83) is $\\int\\frac{U-c_r}{\\lvert U-c\\rvert^2}\\frac{d^2U}
   {dy^2}\\lvert\\phi\\rvert^2dy=-\\int(\\lvert\\phi'\\rvert^2+k^2\\lvert\\phi\\rvert^2)dy<0$ (11.85); because the integral in (11.84) is
   zero for a growing mode, any constant times it can be added: $(c_r-U_I)\\int\\frac1{\\lvert U-c\\rvert^2}\\frac{d^2U}{dy^2}
   \\lvert\\phi\\rvert^2dy=0$ (11.86); the sum gives $\\int\\frac{U-U_I}{\\lvert U-c\\rvert^2}\\frac{d^2U}{dy^2}\\lvert\\phi\\rvert^2dy<0$, so
   $(U-U_I)U''<0$ somewhere — the vorticity magnitude has a **maximum** at the inflection point. ⚠️ The book writes U₁ for U at the
   inflection point; we write U_I (U₁ was the upper stream in C02)."
8. `nb.worked_example("three profiles by hand", "1. **Poiseuille** U = 1 − y²: U″ = −2 everywhere — no sign change, no inflection
   ⇒ inviscidly stable (Rayleigh). 2. **tanh** U = tanh y: U″ = −2 tanh y sech²y — negative for y > 0, positive for y < 0 ⇒
   inflection at y = 0, U_I = 0; (U − U_I)U″ = −2 tanh²y sech²y ≤ 0 ⇒ Fjørtoft satisfied: instability possible (and real: kc_i =
   0.19). 3. **sin y on |y| ≤ b**: U″ = −sin y changes sign at 0, and (U − 0)U″ = −sin²y < 0: both criteria pass — yet for 2b < π
   the flow is stable (N95). The criteria are necessary, not sufficient.")`
9. `nb.code` — `names = ("blasius_like", "couette", "poiseuille", "wall_vorticity_max", "jet", "shear_layer")`; for each
   `p = ch11.inviscid_profile(name)`, `y = np.linspace(*p["domain"], 801)`, print `ch11.rayleigh_criterion(y, p["U"], p["Upp"])`
   and `ch11.fjortoft_criterion(y, p["U"], p["Upp"])`; then `ch11.tanh_max_growth()`. *expect:* the verdict table of Fig. 11.21 —
   (a)–(c): no inflection; (d) wall-vorticity-max: Rayleigh yes, Fjørtoft no; (e) jet and (f) shear layer: both yes; tanh: k =
   0.4449, kc_i = 0.1897. *explain:* 1. our analytic stand-ins of the six profiles; 2. the two necessary conditions; 3. the
   fastest tanh wave (Michalke 1964).
10. `nb.check_agree` — **from scratch (curation §7):** `s = np.sign(Upp(y)); idx = np.nonzero(np.diff(s))[0]` (P180) and the
    Fjørtoft product `np.min((U(y) - U(y[idx[0]]))*Upp(y))`; `assert np.allclose(y[idx], ST.inflection_points(y, Upp=Upp(y)),
    atol=y[1] - y[0])`.
11. `nb.figure` — **"Six profiles, two verdicts"** (our Fig. 11.21 analogue via `six_profiles`): each panel U(y) (orange) with U″
    (muted) and the inflection point ● (purple), badges "Rayleigh ✓/✗", "Fjørtoft ✓/✗" (rose/teal). *see:* "only the bottom two
    pass both"; *read:* "(d) has an inflection point but its vorticity is largest at the wall, not inside"; *change:* "…an adverse
    pressure gradient bends a boundary-layer profile into (d)/(f)-like shapes (Ch. 9 (9.51)–(9.52)): it becomes inviscidly
    unstable (R24, C13)".
12. `nb.note` — **N94 [B]** "**Which flows are suspect?** Jets, wakes, mixing layers and adverse-pressure boundary layers have an
    inflection point with a vorticity maximum: potentially unstable at any Re. Couette, Poiseuille and favourable or zero-gradient
    boundary layers have none: inviscidly stable — any instability they have must come from viscosity (C13)."
13. `nb.note` — **N95 [B]** "**Not sufficient.** U = sin y between walls at y = ±b has an inflection point that satisfies both
    criteria, yet it is stable for 2b < π (Tollmien): " + `nb.code` `[ch11.sin_profile_max_growth(b) for b in (1.4, 1.6, 2.0, 3.0)]`.
    *expect:* 0 for b = 1.4 (2b < π), small and growing beyond (≈ 0.046 at b = 2.0, approximate) — the builder prints the
    computed values.
14. `nb.note` — **N96 [B]** "**Critical layers.** A neutral mode has a real c, and Howard's (11.71) puts it inside the velocity range,
    so U(y_c) = c at some y_c: the *critical layer*, a singular point of (11.81) (P269) where φ may be discontinuous. The full
    equation (11.79) has no such singularity; a thin viscous layer forms there instead, thinning as Re grows." + `nb.code`:
    `ch11.critical_layer(y, np.tanh(y), 0.3)` → y_c = artanh 0.3 = 0.3095.
15. `nb.note` — **N97 [B]** "**Kelvin's cat's eye.** Seen by an observer moving at c, the basic flow is U − c, with stream function
    ∫(U − c)dy; adding the neutral wave, $\\hat\\psi=\\int(U-c)dy+A\\phi(y)\\exp\\{\\mathrm ikx\\}$ (11.87). Near y_c a Taylor
    expansion (U(y_c) = c kills the linear term) gives $\\hat\\psi\\cong\\frac{(y-y_c)^2}2\\big[\\frac{dU}{dy}\\big]_{y=y_c}+A\\phi(y_c)
    \\cos(kx)$: closed eyes of half-width $2\\sqrt{A\\phi_c/U'_c}$ (ours) around the critical layer (our Fig. 11.22 analogue)." +
    `nb.figure` `ch11.cats_eye_streamfunction(X, Y, y_c=0, A=0.1, phi_c=1, k=1, Uy_c=1)` contours, the separatrix bold, the
    computed width marked. *expect:* width = 2√0.1 = 0.632. *see / read / change:* "a chain of eyes along y = y_c" / "fluid inside
    an eye circulates with the wave, outside it passes by — the billows of C02 in their linear infancy" / "…double A: the eyes
    widen by √2".
16. `nb.note` — **N126 [B]** "**A shear layer you can solve by hand (Exercise 11.11).** A piecewise-linear layer (speed U₁ above,
    U₃ below, linear in between, thickness h) gives $c_o^2=\\big(\\frac{U_1-U_3}{2kh}\\big)^2\\{(kh-1)^2-e^{-2kh}\\}$: unstable when
    the brace is negative, neutral at kh = 1.2785 (our root of $(kh-1)^2=e^{-2kh}$) — a V1 target for the Rayleigh solver on a
    smoothed profile." + `nb.code`: `ch11.piecewise_neutral_kh()`, `ch11.piecewise_shear_layer_c(0.797)`. *expect:* 1.2785; c
    purely imaginary with kc_i h/ΔU ≈ 0.2012 at kh = 0.797 (analyst, the builder prints it).
17. `nb.note` — **N127 [B]** "**Symmetric and antisymmetric profiles (Exercise 11.12).** Written for v̂ = −ikφ, Rayleigh's equation
    is $(U-c)\\big(\\frac{d^2\\hat v}{dy^2}-k^2\\hat v\\big)-\\frac{d^2U}{dy^2}\\hat v=0$ (11.95). If U is odd (a shear layer) and c is an
    eigenvalue so is −c* (waves come in mirror pairs); if U is even (a jet) the modes split into sinuous (v̂ even: the jet wiggles)
    and varicose (v̂ odd: it pulses) — `ch11.rayleigh_eigs(…, parity="even"/"odd")` on the Bickley jet (R23)."
18. `nb.explainer("inviscid_shear_criteria", heading="Which profiles can be unstable without viscosity?", why="Switching profiles
    moves the inflection marker, flips the verdict badges and moves the eigenvalues on the c-plane — always inside the semicircle —
    and clicking a neutral mode opens the cat's eye in the frame moving with c: each theorem is checked against an example instead
    of being read.", tries=["Preset 'wall-vorticity-max': Rayleigh ✓ but Fjørtoft ✗ — find where (U − U_I)U″ is positive.",
    "Preset 'sin y' and drag b across π/2: the eigenvalue leaves the real axis only when 2b > π.", "Open the cat's eye and double
    A: the eye width grows by √2."])`
19. `nb.md` — **What would change if…** "…the profile has no inflection point at all, like plane Poiseuille flow? Rayleigh says no
    inviscid instability — but C11's spectrum already showed a growing mode at Re = 10⁴. Viscosity itself must be doing it: C13."

---

### A.10 §11.10 Results for Parallel and Nearly Parallel Viscous Flows — R23 R24, C13 · C14
1. `nb.section("11.10", "Results for Parallel and Nearly Parallel Viscous Flows", intro="**What is this section about?** With the
   Orr–Sommerfeld solver of C11 we can trace, for any parallel flow, the curve in the (Re, k) plane where waves neither grow nor
   decay. Inflectional flows (jets, mixing layers) go unstable at tiny Re; plane Poiseuille flow, which Rayleigh calls stable,
   goes unstable at Re = 5772 because of viscosity; Couette and pipe flow never do, linearly. C14 explains how viscosity can
   destabilise, through the energy budget of the disturbance.")`
2. `nb.recap("R23", "The Bickley jet", "The plane laminar jet of Ch. 9, $U=U_0\\,\\mathrm{sech}^2(y/L)$ ((9.71), `JET.free_jet_profile`).
   **New here:** it is inflectional and very unstable: its sinuous mode goes unstable at Re = U₀L/ν ≈ 4 at kL ≈ 0.2 (Tatsumi &
   Kakutani 1958, approximate); inviscid neutral modes at k = 2 (sinuous) and k = 1 (varicose) with c = 2/3 (Drazin & Reid 1981;
   `bickley_critical`).", where="Ch. 9 (free jets, (9.71))")`
3. `nb.recap("R24", "Boundary layers with a pressure gradient", "Ch. 9's rule (9.51)–(9.52): a favourable pressure gradient keeps
   the profile free of inflection points; an adverse one creates one. **New here** (our Fig. 11.24 analogue): the neutral curves
   of Falkner–Skan profiles — favourable and zero gradient close into a loop as Re → ∞ (inviscidly stable), adverse ones keep a
   flat upper branch (inflectional, inviscidly unstable) and go unstable at lower Re (`falkner_skan_neutral_curve`).", where="Ch. 9
   §9.7")`

#### C13 — Plane Poiseuille flow: inviscidly stable, viscously unstable at Re_c = 5772.22
4. `nb.core("C13", "Plane Poiseuille flow: viscously unstable above $\\mathrm{Re}_c=5772.22$ at $k_c=1.02056$ (Orszag 1971), and the
   Table 11.1 family", question="How can a flow that Rayleigh's theorem calls stable become unstable when viscosity is added?")`
5. `nb.md` — **The problem in plain words:** "Water flowing between two plates has the parabolic profile of Ch. 8 — no inflection
   point, so by C12 no inviscid instability. Yet the Orr–Sommerfeld spectrum of C11 had a growing wave at Re = 10⁴. Where exactly
   does it start, which wavelength comes first, and how do other flows of the same family compare?"
6. `nb.md` — **The idea** (ASCII):
   ```
    k ▲      ____                    neutral curve c_i(k, Re) = 0: a 'thumb'
      │   .-'    '----.___           inside: Tollmien–Schlichting (TS) waves grow
      │  (   unstable      '-----    outside: everything decays
      │   '-.___ ◆_______________    ◆ = the critical point (Re_c, k_c): first instability as Re rises
      └──────────────────────────► log Re
    inviscid profiles (tanh, jets): the thumb reaches Re → 0;  Poiseuille: starts at 5772;  Couette, pipe: no thumb at all
   ```
7. `nb.note` — **N98 [C]** "**Viscosity cuts both ways** (named): at very large Re it damps the shortest waves; near the walls it
   shifts the phase of the motion so that the wave can extract energy (C14). The asymptotic theory (Heisenberg, Lin, Shen, Yih —
   wall and critical layers as singular perturbations) is notoriously hard; we compute."
8. `nb.worked_example("what Re_c means for water", "Channel half-width L = 1 cm, water ν = 10⁻⁶ m²/s. 1. Re = U₀L/ν with U₀ the
   centreline speed: Re_c = 5772 ⇒ U₀ = 5772 × 1e-6/0.01 = 0.577 m/s. 2. The mean speed of a parabola is ⅔U₀ = 0.385 m/s; a
   Reynolds number built on the mean speed and the full width 2L would read (⅔)(2)(5772) = 7696 — always say which length
   and speed (Table 11.1 mixes them, P200). 3. At Re = 10⁴, k = 1 (C11): kc_i = 0.00374 in units U₀/L; with L = 1 cm and U₀ =
   1 m/s that is 0.374 s⁻¹, an e-folding time of 2.7 s, for a wave 2πL = 6.3 cm long travelling at c_r U₀ = 0.24 m/s.")`
9. `nb.code` — `crit = ch11.poiseuille_critical()` (cached; recomputed coarse if absent, with a note); `nc =
   ch11.poiseuille_neutral_curve()`; print crit and the band at Re = 10⁴; `print(ch11.couette_max_growth(1.0, 1e3),
   ch11.couette_max_growth(1.0, 1e4))`. *expect:* Re_c = 5772.22, k_c = 1.02056, c_r = 0.26400 (Orszag 1971: all digits); unstable
   k from ≈ 0.80 to ≈ 1.07 at Re = 10⁴; Couette −0.119, −0.052 (decays). *explain:* 1. P263's two nested searches on the OS
   eigenvalue; 2. the thumb; 3. Couette for contrast.
10. `nb.check_agree` — **from scratch (curation §7):** `Re_mine = brentq(lambda Re: np.imag(ch11.orr_sommerfeld_eigs(1.02056, Re,
    U, Upp, N=80)[0]), 5000, 7000, xtol=1e-3)`; `assert np.isclose(Re_mine, 5772.22, rtol=1e-4)` — "the cached critical point is
    not a black box: one brentq reproduces it (a cached run is invisible to a mutant, ch10 lesson)".
11. `nb.figure` — **"A thumb-shaped neutral curve"** (two panels): (a) plane Poiseuille neutral curve (purple) on log Re with the
    unstable region rose and ◆ at (5772.22, 1.02056); the Re = 10⁴ line with its band; (b) neutral curves of the tanh layer (Re_c
    = 0, k → 1 as Re → ∞), the Bickley jet (Re_c ≈ 4) and Blasius (Re_δ* ≈ 519) on one log axis (cached
    `reference/ch11/os_neutral_*.csv`). *see:* "Poiseuille's thumb starts far right; the inflectional flows' curves reach down to
    tiny Re"; *read:* "inflectional = inviscid instability, viscosity merely trims it; no inflection = instability only through
    viscosity, at large Re"; *change:* "…plane Couette: no curve at all (N100)".
12. `nb.plotly` — **F7** `slider_figure` over Re (24 values log-spaced 10³ … 10⁵, cached `os_spectrum_poiseuille.json`): the OS
    spectrum at k = 1 on the c-plane with the TS eigenvalue (rose) crossing c_i = 0 between Re = 5000 and 6000. *explain:* "the
    Y-shaped spectrum and the one mode that crosses".
13. `nb.note` — **N99 [B]** "**The two-stream shear layer** $U=U_0\\tanh(y/L)$ (our Fig. 11.23 analogue in panel (b)): unstable for
    0 < kL < k_u(Re), with k_u → 1 as Re → ∞ and Re_c = 0 — viscosity is only a small correction to an inviscid instability. The
    inviscid neutral mode at kL = 1 has c = 0 and φ = sech(y/L) (check: residual of (11.81) is zero — the cell does it with sympy)."
    + `nb.code`: `sp.simplify(...)` residual for φ = sech y, U = tanh y, k = 1, c = 0; `ch11.tanh_shear_layer_neutral_curve(
    [10, 100, 1000])`. *expect:* residual 0; k_u rising toward 1.
14. `nb.note` — **N100 [B]** "**Plane Couette flow** is linearly stable at every Re (every OS eigenvalue has c_i < 0 on our grid k ∈
    [0.1, 3], Re ∈ [10², 10⁶]); in experiments it still becomes turbulent near Re ≈ a few hundred (based on half the wall speed
    difference and half the gap) through finite disturbances — named." + `nb.code`: a 6 × 5 grid of `ch11.couette_max_growth`;
    *expect:* all negative.
15. `nb.note` — **N101 [C]** "**Pipe flow** (named): linearly stable too; it becomes turbulent through finite disturbances whose
    threshold depends on how quiet the inlet is; recent work describes the transition as a 'chaotic saddle' (Eckhardt et al. 2007).
    Pointer: Ch. 12."
16. `nb.note` — **N102 [B]** "**Table 11.1, recomputed.**" + `nb.code`: `import pandas as pd; pd.DataFrame(ch11.table_11_1())` —
    columns flow, U(y)/U₀, length scale, our Re_c, our k_c, published benchmark and source: jet (sech², half-width, ≈ 4 at ≈ 0.2,
    Tatsumi & Kakutani 1958), shear layer (tanh, 0), Blasius (δ*, 519.2 at αδ* = 0.303, Thomas; Jordinson 1970: 520), plane
    Poiseuille (half-width, 5772.22, Orszag 1971), pipe (∞, linearly stable, stated), plane Couette (∞, our grid). "The book's
    rounded column is not reproduced (its Poiseuille value is a rounding of Orszag's)." *expect:* the computed rows within the
    stated benchmark tolerances.
17. `nb.note` — **N103 [C]** "**A caveat on 'Re_c = 0'** (named): the mixing layer's zero critical Re assumes a parallel flow; a
    spreading mixing layer has a finite one (Bhattacharya et al. 2006)."
18. `nb.md` — **What would change if…** "…we want to know *why* viscosity can feed the wave? Follow the wave's kinetic energy: what
    produces it and what dissipates it — C14."

#### C14 — The disturbance kinetic-energy equation (11.88)
19. `nb.core("C14", "The disturbance kinetic-energy equation $\\frac{d}{dt}\\int\\tfrac12u_i^2dV=-\\int u_iu_j\\frac{\\partial U_i}{\\partial
    x_j}dV-\\Lambda$ (11.88)", question="Where does a growing wave get its energy, and how can viscosity help it?")`
20. `nb.md` — **The problem in plain words:** "A growing wave gains kinetic energy; something must supply it. The only reservoir is
    the mean flow's shear. And something always takes energy away: viscous dissipation. Bookkeeping the two tells us when a wave
    grows — and reveals the trick by which viscosity, a pure damper, can make the supply larger than the loss."
21. `nb.md` — **The idea** (ASCII):
    ```
     d/dt ∫ ½(u² + v²) dV   =   −∫ uv dU/dy dV      −   Λ
          wave's energy          PRODUCTION            DISSIPATION ν∫(∂u_i/∂x_j)² ≥ 0
                                 Reynolds stress −uv working on the shear U′
     u and v exactly out of phase (u = sin, v = cos): mean of uv = 0 → NO production (inviscid, no inflection)
     viscosity near the wall tilts the phase of v against u → mean uv < 0 where U′ > 0 → production > 0
    ```
22. `nb.primer("the integral of an x-derivative of a periodic function over one period is zero", "If f(x) repeats every λ, then
    $\\int_0^\\lambda\\frac{\\partial f}{\\partial x}dx=f(\\lambda)-f(0)=0$ (fundamental theorem of calculus, P84). Averaged over a
    whole number of wavelengths, every term of the form ∂(…)/∂x disappears — the two ends of the control volume cancel.",
    code="import numpy as np\nx = np.linspace(0, 2*np.pi, 2001)\nf = np.sin(3*x)**2 + np.cos(x)                 # periodic with period 2*pi\nprint(np.trapezoid(np.gradient(f, x), x))        # ~0 (to discretisation error)")` (**P273**)
23. `nb.note` — **N128 [B]** "**The starting line (Exercise 11.13):** the disturbed Navier–Stokes equation in index form,
    $\\frac{\\partial}{\\partial t}(U_i+u_i)+(U_j+u_j)\\frac{\\partial}{\\partial x_j}(U_i+u_i)=-\\frac1\\rho\\frac{\\partial}{\\partial
    x_i}(P+p)+\\nu\\frac{\\partial^2}{\\partial x_j\\partial x_j}(U_i+u_i)$ (11.96)."
24. `nb.derivation("D23", …)` — Part F D23 (12 steps), ref "11.88" (the book leaves it to Exercise 11.13); ★★ with an optional
    `check_src` (`ch11.energy_equation_sympy()["residual"] == 0`).
25. `nb.note` — **N104 [B]** "**The 2-D form and what it says:** $\\frac d{dt}\\int\\tfrac12(u^2+v^2)dV=-\\int uv\\frac{\\partial U}{\\partial
    y}dV-\\Lambda$ (p. 520). The product uv averaged over a wavelength (the *Reynolds shear stress*, Ch. 12) vanishes when u and v
    are 90° out of phase (P178). For an inviscid, inflection-free profile the phases are such that the wave cannot take energy from
    the shear; viscosity changes the phase relation near the walls and the critical layer so that −⟨uv⟩U′ integrates to a positive
    number larger than the dissipation — the destabilising mechanism."
26. `nb.worked_example("why the phase matters", "Average over one period (P178): 1. u = sin t, v = cos t (90° apart): mean uv =
    ½ sin 2t averaged = 0 — no production. 2. u = sin t, v = −sin(t − π/6) (tilted by 30° from anti-phase): mean uv = −½ cos(π/6)
    = −0.433; with U′ = +1 the production −⟨uv⟩U′ = +0.433 > 0. 3. For the TS wave at Re = 10⁴, k = 1 the code below finds P/Λ =
    1.62: production beats dissipation by 62 %, so dE/dt = 2kc_iE > 0. 4. At Re = 5000 the same mode has P/Λ = 0.77: it decays.")`
27. `nb.code` — `m = ch11.os_mode(1.0, 1e4, U, Upp)`; `b = ch11.disturbance_energy_budget(1.0, m["c"], m["phi"], m["y"],
    -2*m["y"], 1e4)`; print every entry; the same at Re = 5000. *expect:* Re = 10⁴: dEdt = 2kc_iE, residual ≈ 1e-12, P/Λ = 1.616;
    Re = 5000: c = 0.26813 − 0.00175i, P/Λ = 0.765. *explain:* 1. the leading OS mode with û = φ′, v̂ = −ikφ (§11.8 convention);
    2. (11.88) averaged over a wavelength with Clenshaw–Curtis weights (P271); 3. the closing check dE/dt = P − Λ, independent of the
    eigen-solve.
28. `nb.check_agree` — **from scratch (curation §7):** interpolate û, v̂ to 4001 uniform y with `np.polynomial.chebyshev.Chebyshev.fit`
    (degree N) on the collocation values (real and imaginary parts separately), then `P = -0.5*np.trapezoid(np.real(uh*np.conj(vh))*(-2*y), y)` and the dissipation with
    `np.gradient`; `assert np.isclose(P, b["production"], rtol=1e-4)` (P203 reminder: trapezoid on a fine uniform grid).
29. `nb.figure` — **"Viscosity tilts the wave so that it can feed"** (three panels): (a) the Reynolds stress −⟨uv⟩(y) (orange) and
    its product with U′ (the local production) across the channel, with the wall layers and the critical layers (where U = c_r, ×)
    marked; (b) the phase of v̂ relative to û across y (purple), the 90° line muted; (c) bars P, Λ, dE/dt = P − Λ (orange, blue,
    rose/teal) for Re = 10⁴ and 5000. *see:* "production concentrated near the walls and the critical layers"; *read:* "where the
    phase departs from 90° the stress is nonzero; at Re = 10⁴ P > Λ"; *change:* "…Re = 5000: P shrinks below Λ (bars flip)".
30. `nb.animation` — **A4** (curation §6): the TS wave at Re = 10⁴, k = 1 travelling at c_r over the parabolic profile:
    perturbation streamlines (`ch11.ts_wave_fields`, amplitude 0.05), the Reynolds-stress profile beside it and P, Λ numbers in the
    title, 60 frames (FAST 30), `player="video"`. Notes: *see:* "tilted cells sliding downstream near each wall"; *read:* "the tilt
    against the shear is the visible sign of −⟨uv⟩ > 0"; *change:* "…k = 0.6 (outside the thumb): the tilt weakens and the wave
    decays".
31. `nb.note` — **N125 [C]** "**With stratification (Exercise 11.10, named):** the budget gains the available potential energy,
    $\\tfrac12\\frac d{dt}\\int\\big(u^2+w^2+\\frac{g^2\\rho^2}{\\rho_0^2N^2}\\big)dV=-\\int uw\\frac{\\partial U}{\\partial z}dV$ — the
    shear's production now has to lift fluid too (C09's energy view). Pointer: available potential energy in Ch. 13."
32. `nb.explainer("orr_sommerfeld_neutral_curve", heading="How can viscosity make a channel flow unstable?", why="Clicking a point
    (Re, k) on the neutral-curve map shows that Tollmien–Schlichting mode travelling over the profile, its production/dissipation
    bars and the u–v phase; crossing the curve flips the balance. Switching flows recomputes Table 11.1 row by row, and the oblique
    slider applies Squire's map live.", tries=["Preset 'Orszag (Re = 10⁴, k = 1)': read c = 0.2375 + 0.0037i and P/Λ ≈ 1.6.", "Drag
    Re down through 5772 at k = 1.02: the bars flip when the dot leaves the thumb.", "Set the oblique angle to 30°: the effective Re
    drops by cos 30° and the 3-D wave falls outside the curve (Squire)."])`
33. `nb.md` — **What would change if…** "…the wave keeps growing? Its own Reynolds stress starts to change the mean flow it feeds
    on — the linear theory's end, §11.12–11.13."

---

### A.11 §11.11 Experimental Verification of Boundary-Layer Instability — **C13 (continued)**: N105 N106 N107 N108 N109
1. `nb.section("11.11", "Experimental Verification of Boundary-Layer Instability", intro="**What is this section about?** The
   Blasius boundary layer on a flat plate (Ch. 9) is nearly parallel, so the Orr–Sommerfeld machinery applies to it as an
   approximation. Its neutral curve was computed long before anyone saw a Tollmien–Schlichting wave; a famous wind-tunnel
   experiment then confirmed it. (These notes continue C13.)")`
2. `nb.note` — **N105 [B]** "**C13 (continued).** **Tollmien's stand-in for Blasius**: a piecewise profile — linear near the wall,
    a parabola reaching 1 with zero slope at y = δ, then 1 — with zero curvature at the wall like Blasius. We fix the breakpoint
    η₁ = 0.2 ourselves and get the coefficients from continuity of value and slope: a = 2/(1 + η₁), b = 1/(1 − η₁²)." + **slip #6
    box**: "the book prints the middle branch as 1 − b[1 − (y/δ)²] — square inside the bracket; that branch jumps at the breakpoint
    (it is near 0 there while the linear branch is near 0.3). With the square outside, 1 − b(1 − y/δ)², value and slope match." +
    `nb.code` `ch11.tollmien_profile(eta)` vs `BL.blasius_profile` (δ = the 99 % thickness, η ≈ 4.91 in
    Blasius units) and the printed variant's jump. *expect:* a = 1.6667, b = 1.0417; jump of the printed form ≈ 0.3.
3. `nb.note` — **N106 [C]** "**C13 (continued).** **Seen at last** (named): TS waves stayed undetected until Schubauer & Skramstad
    (1947) excited them with a vibrating ribbon in a very quiet wind tunnel and measured them with hot wires; their growth and decay
    matched the computed neutral curve."
4. `nb.note` — **N107 [B]** "**C13 (continued).** **The Blasius neutral curve in frequency** (our Fig. 11.26 analogue): experiments
    fix the frequency ω of the ribbon, so the curve is drawn as the dimensionless frequency $F=\\omega\\nu/U_\\infty^2=kc_r/
    \\mathrm{Re}_{\\delta^*}$ against $\\mathrm{Re}_{\\delta^*}=U_\\infty\\delta^*/\\nu$; a wave of fixed F travelling downstream
    (Re_δ* rising like √x) enters the curve, grows, and leaves it." + `nb.figure` from
    `ch11.blasius_neutral_curve(in_frequency=True)` (cached) with a constant-F ray. *see / read / change:* "a tongue in (Re_δ*, F)" /
    "a horizontal line = one ribbon frequency followed downstream" / "…a higher F: the line misses the tongue — those waves never
    grow".
5. `nb.note` — **N108 [B]** "**C13 (continued).** **The Blasius critical point** with the parallel-flow assumption: Re_δ*,c = 519.2
    at αδ* = 0.303 (Thomas; via Gallagher, Griffiths & Stephen 2016; Jordinson 1970: 520). Our base flow uses U = f′ and
    U″ = −½ff″ from the Blasius ODE (Ch. 9), not a numerical second derivative." + `nb.code` `ch11.blasius_critical()` (cached).
    *expect:* within 1 % of 519.2 and 0.303.
6. `nb.note` — **N109 [C]** "**C13 (continued).** **What the parallel model misses** (named): the boundary layer grows downstream
    (non-parallel corrections, Nayfeh & Saric 1975); disturbances can grow for a while even where every eigenmode decays (transient
    growth from the non-normal OS operator, Reshotko 2001); suction thins the layer and delays transition. Pointer: Ch. 12."

### A.12 §11.12 Comments on Nonlinear Effects — **C14 (continued)**: N110
1. `nb.section("11.12", "Comments on Nonlinear Effects", intro="**What is this section about?** Linear theory says a disturbance
   grows like e^{σt} forever; real disturbances stop growing. The reason is that a finite wave changes the flow it lives on.")`
2. `nb.note` — **N110 [B]** "**C14 (continued).** **The wave changes its own supply.** The Reynolds stress −⟨uv⟩ of C14 (and the heat
   flux ⟨uT′⟩ in convection) is a *rectified* (non-zero average) flux: it transports momentum and heat across the layer and so
   modifies U(y) or T̄(z) — usually in the direction that weakens the instability, and the growth stops at a finite amplitude
   (N52's equilibrium cells). Some flows go the other way: finite disturbances destabilise flows that are linearly stable (pipe,
   N101). A rotating annulus heated at the rim shows the whole sequence — steady waves, then irregular ones — the laboratory
   model of the atmosphere's baroclinic waves (pointer: Ch. 13 §13.17)."

### A.13 §11.13 Transition — **C14 (continued)**: N111
1. `nb.section("11.13", "Transition", intro="**What is this section about?** How a laminar flow becomes turbulent: a sequence of
   instabilities, each feeding on the state the previous one left behind.")`
2. `nb.note` — **N111 [B]** "**C14 (continued).** **A sequence of instabilities** (Landau's idea), drawn as two pipelines:
   ```
   boundary layer / channel:  TS waves (2-D, C13) ──► saturated, then a secondary 3-D instability (Klebanoff's
                              peak–valley pattern) ──► Λ-shaped vortices ──► turbulent spots ──► turbulence
   free shear layer:          KH roll-up (C02, A2) ──► vortex pairing ──► 3-D 'braid' vortices ──► turbulence
                              (shear layers transition at Re of order 10; boundary layers need order 10³)
   ```
   Pointer: Ch. 12 (turbulence) and the routes to chaos of C15."

---

### A.14 §11.14 Deterministic Chaos — C15, S02, summary
1. `nb.section("11.14", "Deterministic Chaos", intro="**What is this section about?** Far above onset the steady cells and waves of
   this chapter give way to irregular motion. Remarkably, irregularity does not need many degrees of freedom: three ordinary
   differential equations from a three-mode model of convection (Lorenz 1963) never settle down and amplify any tiny difference in
   the start. This is why weather has a predictability limit of about two weeks, and why forecasters run ensembles.")`

#### C15 — The Lorenz system (11.91)
2. `nb.core("C15", "The Lorenz system $\\dot X=\\Pr(Y-X)$, $\\dot Y=-XZ+rX-Y$, $\\dot Z=XY-bZ$ (11.91), $r=\\mathrm{Ra}/\\mathrm{Ra}_{cr}$,
   $b=4\\pi^2/(\\pi^2+k^2)$", question="How can three exact equations be unpredictable?")`
3. `nb.md` — **The problem in plain words:** "In 1961 Edward Lorenz restarted a small weather model from numbers he had typed with
   three decimals instead of six. Within a couple of simulated months the two runs had nothing in common. He then boiled the
   model down to convection in a layer — the Bénard problem of C05 — keeping only three modes. The three equations still behaved
   the same way: completely determined, and yet unpredictable beyond a horizon."
4. `nb.md` — **The idea** (ASCII):
   ```
   free–free convection (C05)  ──keep 1 roll + 2 temperature modes (11.90)──►  3 ODEs (11.91) for X (flow), Y, Z (temperature)
      r < 1            : conduction (origin attracts)
      1 < r < r_H      : steady rolls C± (two senses of rotation)
      r > r_H ≈ 24.74  : every fixed point unstable ⇒ the state wanders forever between the two rolls: a STRANGE ATTRACTOR
   two starts 10⁻⁸ apart: separation ∝ e^{λt} (λ ≈ 0.9) ⇒ after ~20–30 time units they are unrelated
   ```
5. `nb.note` — **N112 [B]** "**Deterministic chaos** = extreme sensitivity to the initial state + aperiodic motion + a broadband
   spectrum (not a few sharp frequencies). In a *linear* dissipative system steady forcing gives a steady response; with
   nonlinearity it can give a periodic or an aperiodic one (the heated rotating annulus; the cylinder wake that starts to oscillate
   near Re ≈ 40, Ch. 9)." + spectrum figure of X(t) at r = 28 vs r = 10 via `np.fft.rfft` (P181 reminder) with
   `ch10.dominant_frequency` on the r = 10 transient. *see / read / change:* "a broad hump vs a decaying single line" / "chaos =
   energy at all frequencies" / "…r = 24 (just below r_H): long chaotic transients, then a steady state".
6. `nb.note` — **N113 [B]** "**Phase space** with the pendulum: $\\ddot X+(g/l)\\sin X=0$ becomes $\\dot X=Y$, $\\dot Y=-(g/l)\\sin X$
   (11.89); the plane (X, Y) is its phase space, a solution is a trajectory there, and the number of independent initial values (2)
   is its number of degrees of freedom." + `nb.figure` `ch11.phase_portrait(ch11.pendulum_rhs, starts)` with energy contours,
   energy conserved to 1e-9 (printed). *see / read / change:* "closed loops (swinging) and wavy lines (rotating over the top)" /
   "each curve is one constant energy" / "…add damping: every trajectory spirals into the point (0, 0) — an attractor".
7. `nb.primer("Jacobian matrix of a nonlinear ODE system and the stability of its fixed points", "For $\\dot{\\mathbf s}=\\mathbf
   f(\\mathbf s)$ a fixed point has f(s*) = 0. Near it write s = s* + δ: to first order $\\dot\\delta=J\\delta$ with the Jacobian
   $J_{ij}=\\partial f_i/\\partial s_j$ at s* (multivariable Taylor, P98). The fixed point is stable if every eigenvalue of J has
   negative real part (P214, P80); a complex pair crossing the imaginary axis starts an oscillation (a Hopf bifurcation).",
   code="import numpy as np\nPr, r, b = 10.0, 28.0, 8/3\nJ0 = np.array([[-Pr, Pr, 0], [r, -1, 0], [0, 0, -b]])   # Lorenz Jacobian at the origin\nprint(np.sort(np.linalg.eigvals(J0)))     # -22.83, -2.67, 11.83: one positive -> the origin is unstable")` (**P274**)
8. `nb.note` — **N114 [B]** "**Attractors and bifurcations** (our Fig. 11.28 analogue): in a dissipative system trajectories crowd
   onto *attractors* — a fixed point (a steady flow) or a limit cycle (a steady oscillation); as a parameter R passes R_cr a fixed
   point can turn into a repeller and a limit cycle grow around it whose size increases with R − R_cr (a *bifurcation*). Our model
   of that picture is the Hopf normal form $\\dot z=(\\mu+\\mathrm i\\omega)z-\\lvert z\\rvert^2z$ (ours): stable fixed point for μ < 0,
   limit cycle of radius √μ for μ > 0." + `nb.figure` three panels (stable point μ = −0.5; limit cycle μ = 0.5; bifurcation
   diagram amplitude vs μ with `ch11.limit_cycle_amplitude`). *see / read / change:* "a spiral in, a spiral to a circle, a
   parabola-shaped branch" / "linear theory gives the growth at μ > 0; the nonlinear term fixes the final amplitude √μ" / "…the
   sign of the cubic term reversed: no small stable cycle (a subcritical bifurcation — Lorenz's case, D25)".
9. `nb.note` — **N115 [B]** "**Lorenz's truncation**: stress-free walls, rolls along y, $u=-\\partial\\psi/\\partial z$, $w=\\partial
   \\psi/\\partial x$ (⚠️ slip #9: the opposite sign to §11.7), and $\\psi\\propto X(t)\\cos(\\pi z)\\sin(kx)$, $T'\\propto Y(t)\\cos(\\pi z)\\cos(kx)
   +Z(t)\\sin(2\\pi z)$ (11.90), z ∈ [−½, ½]: X measures the roll speed, Y the temperature difference between rising and sinking
   fluid, Z the distortion of the mean temperature profile. cos πz is the free–free mode sin π(z + ½) of C05 (slip #2's correct
   family)."
10. `nb.primer("Galerkin truncation: keep a few modes and project with orthogonality of sines", "Write the unknown as a sum of a few
    fixed shapes with time-dependent amplitudes, insert it, and multiply the equation by each shape and integrate: orthogonality
    ($\\int\\sin m\\pi z\\sin n\\pi z\\,dz=0$ for m ≠ n, P151) keeps one equation per amplitude, and everything that does not fit the
    chosen shapes is thrown away. Ch. 10's finite elements used hats as shapes; here we use sines and cosines.",
    code="import numpy as np\nz = np.linspace(-0.5, 0.5, 20001)\nprint(np.trapezoid(np.cos(np.pi*z)*np.cos(3*np.pi*z), z), np.trapezoid(np.cos(np.pi*z)**2, z))   # 0 and 0.5")` (**P275**)
11. `nb.derivation("D24", …)` — Part F D24 (**★★★, 15 steps**), ref "11.91" (the book writes only "Lorenz finally obtained");
    `check_src` = Part F D24 sympy cell (`ch11.lorenz_sympy()`).
12. `nb.primer("purely imaginary roots of a cubic λ³ + a₂λ² + a₁λ + a₀: exactly when a₂a₁ = a₀", "Put λ = iω: the real part
    gives −a₂ω² + a₀ = 0 and the imaginary part −ω³ + a₁ω = 0, so ω² = a₁ and a₂a₁ = a₀. For positive coefficients, all three roots
    have negative real part when a₂a₁ > a₀ (Routh–Hurwitz); at a₂a₁ = a₀ a pair sits on the imaginary axis — the Hopf point.",
    code="import numpy as np\na2, a1 = 13.6667, 101.3333                       # Lorenz convection-state cubic at r = 28\nprint(a2*a1, 2*(8/3)*10*(28 - 1))                 # 1384.9 < 1440 = a0: a pair has crossed (unstable)\nprint(np.roots([1, a2, a1, 1440.0]))              # -13.85 and 0.094 +- 10.19j")` (**P276**)
13. `nb.derivation("D25", …)` — Part F D25 (11 steps), ref "11.91" (the book says only 'if r is large') — contains **N116**
    (fixed points, the pitchfork at r = 1, the Hopf point r_H).
14. `nb.worked_example("the Lorenz numbers at r = 28", "Pr = 10, b = 8/3, r = 28 (Lorenz 1963). 1. Convection states: X̄ = Ȳ =
    ±√(b(r − 1)) = ±√(8/3 × 27) = ±√72 = ±8.485, Z̄ = r − 1 = 27. 2. Hopf point: r_H = Pr(Pr + b + 3)/(Pr − b − 1) = 10 ×
    15.667/6.333 = 24.74 < 28 ⇒ the convection states are unstable (eigenvalues 0.094 ± 10.19i). 3. Origin: eigenvalues −8/3 and
    the roots of λ² + 11λ − 270 = 0: 11.83 and −22.83 ⇒ unstable. 4. Volume: every blob of states shrinks at the rate
    −(Pr + 1 + b) = −13.67 per unit time, so the attractor has zero volume — yet it is not a point or a loop: it is strange.
    5. b = 8/3 comes from k² = π²/2, the free–free K_c of C05: 4π²/(π² + π²/2) = 8/3.")`
15. `nb.code` — `print(ch11.lorenz_fixed_points(28.0)); print(ch11.lorenz_eigs(28.0, which="C")); print(ch11.lorenz_eigs(28.0,
    which="origin")); print(ch11.lorenz_hopf_r()); print(ch11.lorenz_b(np.pi/np.sqrt(2)))`; `sol = ch11.lorenz_integrate((1.0,
    1.0, 1.0), 40.0)`; `sep = ch11.lorenz_separation(delta0=1e-8)`; `print(ch11.lorenz_predictability_time(1e-8, 1.0))`.
    *expect:* (0, 0, 0), (±8.4853, ±8.4853, 27); −13.8546, 0.0940 ± 10.1945i; 11.8277, −2.6667, −22.8277; 24.7368; 2.6667;
    separation reaches 1 at t ≈ 29.1 from (1, 1, 1) with δ₀ = 1e-8 along X. *explain:* 1. N116 and D25; 2. (11.91) with DOP853,
    rtol 1e-10 (P31/P94); 3. two runs 1e-8 apart.
16. `nb.check_agree` — **from scratch (curation §7):** an RK4 loop for (11.91) (P95), dt = 0.01, 4000 steps from (1, 1, 1);
    `assert np.allclose(X_rk4[:501], np.interp(t[:501], sol["t"], sol["X"]), atol=1e-3)` (agreement to t = 5) and print the first
    time the two differ by more than 1 (≈ 18.5) — "both are correct numerical solutions; their tiny truncation differences grow
    like any other perturbation: the disagreement *is* the lesson"; plus `assert np.allclose(ch11.lorenz_rk4((1, 1, 1), 0.01,
    1)[1], [1.012567, 1.259918, 0.984891], atol=1e-6)`.
17. `nb.primer("Lyapunov exponent: the slope of log(separation) against time", "If two nearby trajectories separate like
    $\\lvert\\delta(t)\\rvert\\approx\\lvert\\delta_0\\rvert e^{\\lambda t}$, then ln∣δ∣ grows linearly with slope λ, the (largest)
    Lyapunov exponent; λ > 0 means chaos. The time to grow from δ₀ to an error of size Δ is about ln(Δ/δ₀)/λ — doubling the
    precision of the start only adds ln 2/λ to the forecast horizon. Fit the slope only while the separation is small (it
    saturates at the attractor's size).", code="import numpy as np\nt = np.linspace(0, 10, 11); sep = 1e-8*np.exp(0.9*t)   # an ideal exponential separation\nprint(np.polyfit(t, np.log(sep), 1)[0])        # 0.9: the slope of ln(sep) is lambda")` (**P277**)
18. `nb.primer("matplotlib 3-D line plots", "`ax = fig.add_subplot(projection='3d')` makes 3-D axes; `ax.plot(X, Y, Z)` draws a
    curve through space, `ax.view_init(elev, azim)` sets the viewpoint. (Ch. 1–2's 3-D figures were plotly, P41/P64; for a static
    page figure matplotlib is lighter.)", code="import numpy as np, matplotlib.pyplot as plt\nfig = plt.figure(figsize=(4, 3)); ax = fig.add_subplot(projection='3d')\nt = np.linspace(0, 6*np.pi, 300); ax.plot(np.cos(t), np.sin(t), t/10)   # a helix\nax.view_init(20, 35)")` (**P278**)
19. `nb.figure` — **"The Lorenz attractor and two runs that part"** (three panels, our Figs. 11.29–11.30 analogues): (a) 3-D attractor
    X–Y–Z (orbit faint muted, the last 5 time units bold orange) with C± and the origin marked; (b) X(t) for the two runs 10⁻⁸ apart
    (orange and blue); (c) log₁₀ of the separation with the fitted slope (rose) and its window printed. *see:* "two lobes; two X(t)
    curves on top of each other until t ≈ 25, then unrelated; a straight rising line in (c) that flattens"; *read:* "X switches lobes
    irregularly (the roll reverses its sense); the slope in (c) is λ (our fit ≈ 0.6–0.9 depending on the window; literature 0.906)";
    *change:* "…δ₀ = 10⁻¹⁶ (double precision): the two runs agree only ≈ 20 time units longer — ln(10⁸)/λ".
20. `nb.animation` — **A5** (curation §6): two trajectories 10⁻⁸ apart moving on the X–Z projection, the convection roll (ψ, T′
    from `ch11.lorenz_fields`) driven by run 1 beside it, and log₁₀∣δ∣(t) underneath, 90 frames (FAST 40), `player="video"`.
    Notes: *see:* "two dots glued together, then flying apart; the roll reversing"; *read:* "the reversal of the roll is a lobe
    switch"; *change:* "…r = 20: both dots spiral into the same convection state — no chaos (transient only)".
21. `nb.plotly` — **F8** `slider_figure` over r (30 values 0.5 … 30, cached `ch11.lorenz_r_sweep`): the X–Z trajectory after a
    transient, with the verdict (conduction / steady convection / transient chaos / chaos) in the title.
22. `nb.live` — free r, Pr, b, δ₀: the attractor and the separation — paired with F8.
23. `nb.primer("iterated maps: fixed point, stability ∣f′(x*)∣ < 1, cobweb diagram", "A map $x_{n+1}=f(x_n)$ is a rule applied
    again and again. A fixed point has x* = f(x*); a small error δ becomes f′(x*)δ next step, so x* attracts if ∣f′(x*)∣ < 1
    (compare ch10's ∣G∣ ≤ 1). A cobweb draws the iteration: up to the curve y = f(x), across to the diagonal y = x, repeat.",
    code="A = 2.8; x = 0.2\nfor n in range(30):\n    x = A*x*(1 - x)                  # logistic map\nprint(round(x, 6), 1 - 1/A, abs(2 - A))  # 0.642857 = fixed point 1 - 1/A; |f'| = 0.8 < 1: stable")` (**P279**)
24. `nb.note` — **N117 [B]** + **N129 [B]** "**One route to chaos: period doubling.** The logistic map $x_{n+1}=Ax_n(1-x_n)$
    (Exercise 11.14, a transition toy; background state x* = 1 − 1/A) has a stable fixed point for 1 < A < 3 (∣f′(x*)∣ = ∣2 − A∣ <
    1), a stable 2-cycle from A₁ = 3, a 4-cycle from A₂ = 1 + √6 = 3.4495, and so on; the gaps between doublings shrink by a
    universal ratio, $\\frac{R_n-R_{n-1}}{R_{n+1}-R_n}\\to4.6692$ (p. 530; Feigenbaum 1978), the same in any one-hump map and in
    convection experiments (our Fig. 11.31 analogue). Our estimate from superstable parameters (where x = ½ is on the cycle):
    4.709, 4.681, 4.663, 4.668, 4.669." + `nb.code` `ch11.period_doubling_points(6)`, `ch11.feigenbaum_estimate(6)` + `nb.figure`
    (two panels): cobweb at A = 2.8, 3.2, 3.5 (`ch11.cobweb`), and the bifurcation tree `ch11.bifurcation_diagram(np.linspace(2.8,
    4, 1200))` with A₁, A₂, A₃ marked. *expect:* A_n = 3, 3.449490, 3.544090, 3.564407; δ estimates as above. *see / read /
    change:* "a single line splitting into 2, 4, 8 … then a cloud" / "each split comes sooner by the factor 4.669" / "…a different
    one-hump map (sine map in B1's chips): the same ratio".
25. `nb.note` — **N118 [C]** "**Other routes** (named): quasi-periodicity — after two incommensurate frequencies the motion becomes
    chaotic (Ruelle & Takens 1971), unlike Landau's picture of an endless sequence of new frequencies; the Lorenz system follows
    neither — its chaos appears after a subcritical Hopf point (D25). Pointer: the turbulence-onset debate in Ch. 12."
26. `nb.note` — **N119 [C]** "**What chaos means for prediction** (named): a deterministic system can still be unpredictable in
    practice — not because of quantum uncertainty but because no measurement is exact (Poincaré saw this in 1908). Climate hook:
    weather forecasts lose skill after about two weeks; centres therefore run *ensembles* of forecasts from slightly different
    starts and forecast probabilities; climate projections ask about statistics of the attractor, not one trajectory (Ch. 13). Whether
    turbulence itself is 'chaos in a few modes' remains open."
27. `nb.explainer("lorenz_attractor", heading="How can three exact equations be unpredictable?", why="The 3-D orbit, the X(t) traces
    and the log-separation panel run on one clock from two starts 10⁻⁸ apart; the r slider moves the fixed points and changes the
    verdict while you watch; a static attractor picture hides both the divergence and the r-dependence.", tries=["Preset 'r = 10':
    both runs spiral into a convection state; raise r past 24.74 and watch the status turn rose.", "At r = 28 read the time at which
    the runs differ by 50 % (end card); make δ₀ 100 times smaller and see how few time units you gain.", "Click the orbit: the
    inspector shows the Jacobian's eigenvalues there."])`
28. `nb.pointer("**S02** — Literature: the benchmarks used here are cited where they appear (Chandrasekhar 1961; Orszag 1971;
    Thomas, via Gallagher, Griffiths & Stephen 2016; Jordinson 1970; Tatsumi & Kakutani 1958; Michalke 1964; Drazin & Reid 1981;
    Lorenz 1963; Feigenbaum 1978); the book's bibliography is not reproduced.")`
29. `nb.summary(clicked=[…15…], feeds_forward=[…], left_out=[…])` — **clicked** (one line per CORE): C01 "A small disturbance is a
    sum of waves e^{ikx+σt}; a flow is stable only if σ_r ≤ 0 for every k, and onset happens where σ_r(k) first touches zero." ·
    C02 "Across a velocity jump, gravity's restoring term fights shear's kinetic term under a square root; when shear wins the two
    wave speeds become a complex pair and one grows — short waves first, unless surface tension holds them (6.7 m/s for wind over
    water)." · C03 "Heating from below gives two ODEs in z with one parameter Ra; σ is a real eigenvalue, computed by Chebyshev
    collocation." · C04 "Each cell width has its own marginal Ra; the valley's bottom, Ra_c = 1707.76 at K_c = 3.117, is the onset,
    with cells about as wide as the layer is deep." · C05 "With stress-free walls sine modes make everything algebra:
    Ra = (π² + K²)³/K², minimum 27π⁴/4 at K² = π²/2." · C06 "Because heat diffuses ~100× faster than salt, a column lighter on top
    can still overturn into salt fingers when Rs − Ra > 27π⁴/4." · C07 "Swapping rings releases energy when Γ² falls outward; the
    centrifugal force plays gravity and viscosity adds the 1708-type threshold Ta_c ≈ 1708/(½(1 + Ω₂/Ω₁))." · C08 "Stratified
    shear flow obeys one ODE, Taylor–Goldstein, with c as the eigenvalue; its unstable modes come in conjugate pairs." · C09
    "Ri > ¼ everywhere guarantees stability (Miles–Howard); Ri < ¼ only allows instability." · C10 "Every unstable wave speed lies in
    the half-disc on [U_min, U_max]; growth is at most k(U_max − U_min)/2." · C11 "Squire makes 2-D waves the first to go unstable;
    eliminating pressure leaves the fourth-order Orr–Sommerfeld equation." · C12 "Inviscid instability needs an inflection point with
    a vorticity maximum; neutral waves have a critical layer and cat's eyes." · C13 "Plane Poiseuille flow, inviscidly stable, is
    viscously unstable above Re_c = 5772.22 (k_c = 1.02056) — Tollmien–Schlichting waves." · C14 "A disturbance grows when its
    Reynolds-stress production −∫uvU′ beats viscous dissipation; viscosity can supply the phase shift that makes production
    positive." · C15 "Three modes of convection give the Lorenz system; past r_H ≈ 24.74 it is chaotic — deterministic yet
    unpredictable beyond ~20 time units." **feeds_forward:** "Ch. 12: Reynolds stresses −⟨uv⟩ and the production term of the
    turbulent kinetic energy (C14), transition routes (N111), Ri-based mixing (C09)" · "Ch. 13: baroclinic instability (C01's recipe
    with rotation; the semicircle and inflection theorems reappear with U″ − β, Rayleigh–Kuo), inertial instability (C07), stratified
    shear (C08–C09), predictability (C15)" · "Ch. 14: boundary-layer transition on airfoils (C13)". **left_out:** "exercises 11.3–11.5
    (S01)", "pipe-flow stability analysis (named, N101)", "spatial instability and absolute/convective instability (named, N03)",
    "weakly nonlinear theory (named, N110)".

### A.15 Placement check (every curation id has exactly one home)
- **C:** C01 A.2 2 · C02 A.3 13 · C03 A.4 4 · C04 A.4 28 · C05 A.4 51 · C06 A.5 2 · C07 A.6 5 · C08 A.7 5 · C09 A.7 21 · C10 A.7 34 ·
  C11 A.8 4 · C12 A.9 2 · C13 A.10 4 (+ continued A.11) · C14 A.10 19 (+ continued A.12, A.13) · C15 A.14 2 — **15**.
- **R:** R01–R11 A.3 2–12 · R12 R13 A.4 2–3 · R14 R15 R16 A.6 2–4 · R17 R18 R19 A.7 2–4 · R20 A.7 20 · R21 R22 A.8 2–3 · R23 R24 A.10
  2–3 — **24**, each before the `nb.core` of the block that uses it.
- **S:** S01 A.3 38 · S02 A.14 28.
- **N:** §11.1 N01 N02 N03 (A.1) · C01 N04 N05 N06 · C02 N07–N23 (N08–N12 inside D02, N13–N17 N19 inside D03 and stated) N120 N121 ·
  C03 N24–N40 N122 (N29 N30 in D05, N32–N34 in D06, N36–N38 in D07, N40 N122 in D08) · C04 N41–N46 N50 N51 N52 N123 · C05 N47 N48
  N49 · C06 N53–N57 · C07 N58–N69 N124 · C08 N70–N74 · C09 N75–N78 · C10 N79–N82 · C11 N83–N90 · C12 N91–N97 N126 N127 · C13 N98–N103
  (+ N105–N109 continued) · C14 N104 N125 N128 (+ N110, N111 continued) · C15 N112–N119 N129 — **129**.
- **D:** D01 C01 · D02 D03 D04 C02 · D05 D06 D07 D08 C03 · D09 D10 C04 · D11 D12 C05 · D13 C06 · D14 D15 C07 · D16 D17 C08 · D18 C09 ·
  D19 C10 · D20 D21 C11 · D22 C12 · D23 C14 · D24 D25 C15 — **25**, each inside its curation CORE block (C13 has none, as in the
  curation).
- **Explainers:** E1 C01 · E2 C02 · E3 C04 · E4 C06 · E5 C07 · E6 C09 · E7 C12 · E8 C14 · E9 C15 — **9**, each embedded once. B1 not
  embedded.
- **Animations:** A1 C01 (§11.1 cell, before the C01 block — it is a §11.1 note visual; the builder may move it into C01 if nbkit needs
  every animation inside a block) · A2 C02 · A3 C04 · A4 C14 · A5 C15. **Plotly:** F1 C02 · F2 C04 · F3 C06 · F4 C07 · F5 C09 · F6 C10 ·
  F7 C13 · F8 C15. **Live:** C04, C15.
- **Primers P255–P279 (25):** P255 C01 · P256 C02 · P257 P258 P259 P260 C03 · P261 P262 P263 P264 C04 · P265 C05 · P266 C06 · P267 C07 ·
  P268 P269 P270 C08 · P271 C09 · P272 C10 · P273 C14 · P274 P275 P276 P277 P278 P279 C15 (first-use order).
- **Visual per CORE block:** C01 fig + E1 (+ A1) · C02 figs + F1 + A2 + E2 · C03 figs · C04 figs + F2 + live + A3 + E3 · C05 fig · C06
  fig + F3 + E4 · C07 figs + F4 + E5 · C08 fig · C09 fig + F5 + E6 · C10 fig + F6 · C11 fig · C12 figs + E7 · C13 figs + F7 · C14 fig
  + A4 + E8 · C15 figs + A5 + F8 + live + E9.
- **Code per CORE block:** every block has ≥ 2 `nb.code` cells and one `check_agree` (from scratch, parity or — C08 — sympy).

---

## Part B — explainer storyboards

Common to all nine: created with `tools/new_viz.py`; `<meta name="viz:chapter" content="ch11">`; tabs Walkthrough · Explore ·
Explain · Derivation · Equations · Code · Check; every displayed number is computed by a JS function that mirrors a `ch11`
callable and is proved by `selftest()` parity rows (`py:` expressions use only `ch11.…`, `np.…`, numbers, strings, lists, dict keys
and integer indices — Part C convention 2 and C.5). Explain is "Explanation & interpretation" in numbered sections built with
`Viz.work.step / line / box / table / hint / interpret`, modelled on `forced_damped_vibrations.html` (the regime-dependent reading)
and `amplitude_phase_second_order_II_3.html` (a numbered derivation with live numbers, one section per optional view): **0** what the
views show and what each colour means · **1…n** every displayed quantity from the controls ("formula = substituted = result — why",
results boxed) · a section per view hidden on phones, or a hint · the values at the current time (live) · **Reading the current
setting** (regime-dependent). Derivation steps are copied from Part F (same `did` titles, same step count, same order; phones shorten
*why* to its first sentence; plain-text *why* and *watch* never contain raw TeX). Every tour, Explain, Derivation, notes, status,
equation and quiz text that names a book equation **writes it out** next to its number (`tools/eq_refs.py` → 0; no bare numbers in
`<meta>` strings). Drafts below write equations in Unicode for readability; builders set each in TeX (backslashes doubled in JS
strings; never through a shell heredoc). Colours as in convention 7. Walkthrough texts ≤ 45 words, step 1 ≤ 24 words, ≤ 2 extras per
step, `play: false` on steps that quote numbers, `autoplay: false`. A view hidden on portrait phones never carries a step's key number
(repeat it in a visible title or readout). Symbols as in convention 4 (upright $\mathrm i$, $\Gamma_T$, $U_I$, $\mu$, $\sigma_s$).
**Library note:** there is no complex or linear-algebra helper in `assets/viz_lib.js` beyond `Viz.num.C` (add, sub, mul, div, exp,
log, pow, conj, abs, arg); builders write local `csqrt`, `ccos`, `ccosh`, `csinh`, `csin`, `cdet3` (E1–E3), and no explainer solves
an eigenproblem live: spectra and neutral curves come from tables written by `scripts/ch11_tables.py` / `ch11_neutral_curves.py`
(`reference/ch11/explainer_tables.json`, ≤ 4 s.f., labelled "ours, computed by fluidpy") with ≥ 2 parity rows each. Parallel builders
use private scratch subfolders (`<scratchpad>/<slug>/`).

### E1 · normal_mode_growth
- **Title:** "When is a flow 'unstable'?" · **Summary:** "A disturbance is a sum of waves e^{ikx+σt}: drag k along the σ_r(k) curve
  and watch the wave grow, decay or travel; raise the control parameter and see which wavelength crosses zero first." · **CORE:**
  C01, C05 (also N04, N05, N06, R09 (11.19), N121, R10 (11.20), N40) · **Reference:** `angular_frequency_explorer_1.html` (modes +
  linked views on one clock, presets) + `forced_damped_vibrations.html` (Explain with regime reading).
- **meta:** `viz:order 1` · `viz:sections 11.2 11.3 11.4` · `viz:equations 11.1 11.19 11.20 11.44` · `viz:fluidpy
  ch11.normal_mode_growth ch11.sigma_from_c ch11.stability_class ch11.marginal_type ch11.benard_free_free_sigma ch11.kh_phase_speed`
  · `viz:derivations D01 D12`.
- **Physics (JS ↔ Python):** `sigmaInterface(k, rho1, rho2, g, st)` = −ik c₊ with c₊ from the static two-layer relation with
  surface tension, c² = (ρ₂ − ρ₁)g/((ρ₁ + ρ₂)k) + σ_s k/(ρ₁ + ρ₂) (principal complex root, C.5 5.1) ↔ `ch11.normal_mode_growth
  ("interface", k, rho1=…, rho2=…, g=…, surface_tension=…)`; `sigmaKH(k, U1, U2, rho1, rho2, g, st)` ↔ `normal_mode_growth("kh", …)`
  (via `ch11.kh_phase_speed`); `sigmaBenard(K, Ra, Pr)` = larger root of (σ + a²)(σ/Pr + a²)a² = RaK², a² = π² + K² ↔
  `ch11.benard_free_free_sigma(K, Ra, Pr)[0]`; `stabilityClass(sigma)` ↔ `ST.stability_class`; `marginalType(sigma)` ↔
  `ST.marginal_type`; `sigmaFromC(K, c)` ↔ `ST.sigma_from_c`.
- **Views** (rows [1.1, 1]): 1. `scene` "The disturbance" (row 0, flex 1.3): interface/KH systems — the two layers (upper orange-tint,
  lower blue-tint), the interface $\eta=\zeta_0e^{\sigma_rt}\cos(kx+\sigma_it)$ over two wavelengths, a ghost of the t = 0 shape
  (muted dashed), crest arrows showing c_r; Bénard — the layer with roll streamlines and T′ colour, amplitude $e^{\sigma t}$ (rose
  frame if growing, teal if decaying); amplitude capped at 0.3 of the view height with "×e^{σt} = …" printed. 2. `curve` "σ_r(k) and
  σ_i(k)" (row 0, flex 1): log k axis (interface: 10 … 2000 m⁻¹; Bénard: K 0.3 … 10), σ_r bold (rose above 0, teal below), σ_i
  dashed purple, the zero line, the unstable band shaded rose, the current k dot, the fastest k ◆ and the cut-off ▲; pointer: click
  → set k (inspector). 3. `plane` "The σ-plane" (row 1, `hidePortrait`): Re σ horizontal, Im σ vertical, both roots as dots moving
  with k (faint track of the whole k range), the imaginary axis drawn purple as the stability border.
- **Controls (≤ 5 visible):** `system` chips interface / shear (KH) / Bénard (free walls) (mode) · `lk` "Wavenumber $\log_{10}k$"
  (interface/KH: 1 … 3.3, step 0.01, default 2.32 (k = 210 m⁻¹); Bénard: `K` 0.3 … 10, default 2.22) · control parameter:
  interface `ratio` "$\rho_1/\rho_2$" log 0.001 … 1000, default 0.0012 (air over water; the heavier fluid is
  1000 kg/m³: ratio ≤ 1 sets ρ₂ = 1000, ρ₁ = 1000·ratio, ratio > 1 sets ρ₁ = 1000, ρ₂ = 1000/ratio — 833 is water over air); KH `dU` "$\Delta U=U_1-U_2$" 0 … 12 m/s,
  default 5; Bénard `Ra` log 100 … 10⁴, default 657.5 · `Pr` (Bénard, optional) 0.1 … 100 log, default 1 · `tension` toggle
  "surface tension σ_s = 0.074 N/m" (optional, default on) · transport `t`.
- **Transport:** `t` 0 … 4/∣σ_r∣ (capped to 0.2 s for interface, 0.5 in units d²/κ for Bénard), rate auto, `end: 'hold'`; end card
  (`Viz.card`): "amplitude × e^{σ_r t} = 54.6 (four e-foldings)" or "decayed to 2 %".
- **Presets:** "stable interface (air over water)" {system: interface, ratio: 0.0012, lk: 2} · "upside down: Rayleigh–Taylor"
  {ratio: 833, lk: 2.32} · "vortex sheet" {system: kh, dU: 1, ratio: 1, tension: false} · "Bénard marginal Ra = 657.5"
  {system: benard, Ra: 657.51, K: 2.2214} · "Bénard Ra = 2000" {Ra: 2000, K: 2.2214}.
- **Status:** "🌊 neutral: σ = ±… i — the wave travels at c = … m/s, nothing grows" (all σ_r ≤ 1e-12, σ_i ≠ 0) · "📈 unstable:
  σ_r = 37.0 s⁻¹ at λ = 2.99 cm; band λ > 1.73 cm" (RT preset) · "⚖️ marginal: σ = 0 at K = π/√2 — stationary onset" (Bénard at Ra_c)
  · "📉 stable: every k decays (max σ_r = …)". Exact strings from `ST.stability_verdict` / `ST.marginal_type` pinned by parity rows.
- **Readouts:** "Growth σ_r" [1/s or κ/d²] · "Frequency σ_i" · "Wave speed c_r" · "Wavelength λ" · "e-folding time".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **modes** (three systems) · **presets** ·
  **status** · **inspector** (click the σ curve: "k = 210 m⁻¹: σ² = [gk(ρ₁ − ρ₂) − σ_s k³]/(ρ₁ + ρ₂) = [9.81 × 210 × 998.8 −
  0.074 × 210³]/1001.2 = 1370 ⇒ σ_r = 37.0 s⁻¹").
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**The disturbance**: one normal mode, its amplitude multiplied by e^{σ_r t} (rose frame: growing,
     teal: decaying), its crests moving at c_r = −σ_i/k. **σ_r(k)**: one growth rate per wavelength — bold rose above the zero line
     means that wavelength grows. **The σ-plane** (hidden on phones): the roots; the purple vertical axis is the border between
     decay (left) and growth (right)."
  1. *The mode you picked* — "u = û e^{ikx + σt} (11.1) with k = **210** m⁻¹ (λ = 2π/k = **2.99** cm)."
  2. *Its σ with your numbers* — interface: "c² = (ρ₂ − ρ₁)g/((ρ₁ + ρ₂)k) + σ_s k/(ρ₁ + ρ₂) = … ⇒ c = …; σ = −ikc (D01) = **…**";
     Bénard: "a² = π² + K² = **14.80**; (σ + a²)(σ/Pr + a²)a² = RaK² (D12) ⇒ σ = **…**, **…** (both real — exchange of
     stabilities)". Boxed.
  3. *The whole curve* — "fastest k = **…** with σ_r = **…**; growing band **…**; cut-off where σ_r = 0: k_c = √(gΔρ/σ_s) = **364**
     m⁻¹ (λ_c = 2π√(σ_s/(Δρ g)) = **1.73** cm, Exercise 11.2) / Bénard: Ra(K) = (π² + K²)³/K² (11.44) = **…** at your K".
  4. *Stable, neutral, unstable* — the N04 table with your σ's row lit (stable: σ_r < 0 for every k; neutral: σ_r = 0; unstable:
     σ_r > 0 for some k) and the N06 line "σ_i at the margin = **0** ⇒ stationary (cells) / ≠ 0 ⇒ oscillatory".
  5. *The σ-plane* (or the hint "turn the phone sideways to see the roots move") — "the two roots are complex conjugates for the
     interface (N17) / both real for Bénard (D08)".
  6. *At the current time* — "t = `Viz.live('t')`; amplitude factor e^{σ_r t} = `Viz.live('amp')`; crest moved c_r t = …".
  7. *Reading the current setting* — interface stable: "Heavy below, light above: every wave just travels. This is Ch. 7's interface
     wave, ω² = gk(ρ₂ − ρ₁)/(ρ₂ + ρ₁) (7.95)." · RT: "Heavy on top: long waves fall fastest-growing near λ ≈ 3 cm; surface tension
     saves only the waves shorter than 1.73 cm — the reason water drips from a ceiling in drops of about that size." · vortex sheet:
     "Equal densities: c = (U₁ + U₂)/2 ± i(U₂ − U₁)/2 (11.20) — every wavelength grows, the shortest fastest." · Bénard below Ra_c:
     "Ra < 657.5: no K grows." · at Ra_c: "The curve touches zero at K = π/√2 only: that is the cell width you would see first." ·
     above: "A band of K grows; the fastest is near the middle of the band."
- **Derivation tab:** **D01** (5 steps) `view: 'curve'`, goal `set` {system: interface, ratio: 833}; step 3 (σ = −i∣K∣c) `live`
  "σ = −i × 210 × (… i) = 37.0"; step 5 `watch` "the rose dot sits where c_i > 0". **D12** (8 steps) `view: 'curve'`, goal `set`
  {system: benard, Ra: 2000, K: 2.2214, Pr: 1}; step 2 `live` "a² = π² + K² = 14.80"; step 8 `live` "σ = 11.02, −40.62 at Ra = 2000" and
  `watch` "set Ra = 657.5: the bold curve touches zero at K = 2.22". Interpret: `s => "At your k the mode " + verdict + "."`.
- **Code:**
  ```python
  sig = ch11.normal_mode_growth("{{sys}}", {{k}}, {{params}})   # σ = {{sig}}
  print(ST.stability_class(sigma=sig), ST.marginal_type(sig))    # {{cls}}, {{mt}}
  c = ST.c_from_sigma({{k}}, sig)                                 # c = iσ/k = {{c}}
  ks = np.logspace({{k0}}, {{k1}}, 200)
  growth = np.real([ch11.normal_mode_growth("{{sys}}", kk, {{params}}) for kk in ks])
  print(ST.stability_verdict(growth))                             # {{verdict}}
  ```
- **Walkthrough (6 steps):** 1. "A ripple on an interface" — "Will a small ripple grow or die? Press ▶ and watch one wave." `play: true`
  · 2. "One number per wavelength" — "Every wavelength evolves alone with its own σ. Real part: growth; imaginary part: travel. Drag
  k." `controls: ['lk']`, `derive: {id: 'D01', step: 3}` · 3. "Turn it upside down" — "Heavy water over air: σ_r > 0 for long waves
  (37 s⁻¹ at λ ≈ 3 cm); surface tension stops waves below 1.73 cm." `set` {ratio: 833, lk: 2.32}, `readouts: ['sr', 'lam']` · 4.
  "Stable means every k" — "One growing wavelength is enough for 'unstable' (P255). The rose band is where the flow fails." `highlight:
  ['view:curve']`, `inspect: true` · 5. "Bénard: the curve rises through zero" — "Heat a layer: raise Ra. The σ_r(K) curve first
  touches zero at K = π/√2 when Ra = 27π⁴/4 ≈ 657.5." `set` {system: benard, Ra: 657.51, K: 2.2214}, `derive: {id: 'D12', step: 8}` ·
  6. "Your turn" — "Predict: at Ra = 2000, which K grows fastest? Then find the band edges (≈ 0.75 and 5.36)." `controls: ['Ra',
  'K']`.
- **Equations:** `mode` "The normal mode" ref 'Eq. (11.1)' $u=\hat u(z)e^{\mathrm ikx+\mathrm imy+\sigma t}=\hat u(z)e^{\mathrm
  i\lvert\mathbf K\rvert(\mathbf e_K\cdot\mathbf x-ct)}$, live "σ = −i∣K∣c = …" · `static` "Interface waves" ref 'Eq. (11.19)'
  $c=\pm\big[\frac{\rho_2-\rho_1}{\rho_2+\rho_1}\frac gk\big]^{1/2}$ · `sheet` "Vortex sheet" ref 'Eq. (11.20)' $c=\frac{U_2+U_1}2\pm
  \mathrm i\frac{U_2-U_1}2$ · `ff` "Free–free marginal curve" ref 'Eq. (11.44)' $\mathrm{Ra}=(n^2\pi^2+K^2)^3/K^2$ live · `growth`
  "Free–free growth rate (ours, D12)" ref 'D12 (no book number)' $(\sigma+a^2)(\sigma/\Pr+a^2)a^2=\mathrm{Ra}K^2$ · symbols k
  (1/m), σ (1/s), c (m/s), ρ (kg/m³), σ_s (N/m), Ra, Pr (–).
- **Check yourself:** (1) "With the heavy fluid on top and no surface tension, which wavelength grows fastest?" — "The shortest:
  σ_r = √(gkΔρ/(ρ₁ + ρ₂)) grows with k; only surface tension caps it." `set {ratio: 833, tension: false}` · (2) "Why is the Bénard
  onset 'stationary' and the interface wave 'oscillatory'?" — "At Ra_c the Bénard σ is exactly 0 (real); the interface roots are
  ±iω (pure travel)." · (3) "At Ra = 2000 (free walls, Pr = 1), is K = 6 unstable?" — "No: Ra(6) = (π² + 36)³/36 = 2672 > 2000,
  outside the band 0.75–5.36." `set {system: benard, Ra: 2000, K: 6}` · (4) "Does changing Pr move the Bénard onset?" — "No: σ = 0
  makes Pr drop out of D12's equation; it only changes the growth rate above onset." `set {Pr: 7}`.
- **Selftest parity rows:** `{name: 'KH sigma real', js: sigmaKH(1, 6, 0, 1, 3, 10, 0)[0], py: 'np.real(ch11.normal_mode_growth("kh",
  1.0, U1=6.0, U2=0.0, rho1=1.0, rho2=3.0, g=10.0))', rtol: 1e-10}` · `{name: 'RT sigma at k=210', js: sigmaInterface(210, 1000, 1.2,
  9.80665, 0.074)[0], py: 'np.real(ch11.normal_mode_growth("interface", 210.0, rho1=1000.0, rho2=1.2, surface_tension=0.074))', rtol:
  1e-9}` · `{name: 'Benard sigma', js: sigmaBenard(2.2214415, 2000, 1), py: 'ch11.benard_free_free_sigma(2.2214415, 2000.0, 1.0)[0]',
  rtol: 1e-10}` · `{name: 'sigma from c', js: sigmaFromC(2, [1, 0.5])[1], py: 'np.imag(ch11.sigma_from_c(2.0, 1+0.5j))', rtol: 1e-12}`
  · invariant `{name: 'marginal sigma ~ 0', js: sigmaBenard(Math.PI/Math.SQRT2, 27*Math.PI**4/4, 1), expect: 0, atol: 1e-9}`.
- **Fit plan:** 360×640: status (1 line), `scene` (55 %) over `curve` (45 %), `plane` hidden (σ in the readouts); chips wrap to two
  lines; the walkthrough card paged. 844×390: `scene` | `curve`. 1000×700 and desktop: rows [1.1, 1] with `plane` under `curve`.

### E2 · kelvin_helmholtz_boundary
- **Title:** "Why does wind raise waves only above a threshold?" · **Summary:** "Gravity holds long waves, surface tension short
  ones, shear pushes them all: watch the two wave speeds collide and turn into a growing/decaying pair as the shear crosses the
  boundary." · **CORE:** C02 (also R09, R10, N16, N17, N18, N19, N120, N121; R11 named) · **Reference:** `forced_damped_vibrations.html`
  (Explain with numbers, regime reading) + `amplitude_phase_second_order_II_3.html` (a complex-plane window like its Nyquist view).
- **meta:** `viz:order 2` · `viz:sections 11.3` · `viz:equations 11.9 11.13 11.15 11.16 11.17 11.18 11.19 11.20` · `viz:fluidpy
  ch11.kh_phase_speed ch11.kh_discriminant_terms ch11.kh_critical_k ch11.kh_stability_boundary ch11.kh_min_shear ch11.kh_unstable_band
  ch11.vortex_sheet_c` · `viz:derivations D02 D03 D04`.
- **Physics:** `khRoots(k, U1, U2, rho1, rho2, g, st, h)` (with coth kh when h finite) ↔ `ch11.kh_phase_speed`; `khTerms(...)` ↔
  `ch11.kh_discriminant_terms`; `boundary(k, rho1, rho2, g, st, h)` ↔ `ch11.kh_stability_boundary`; `minShear(rho1, rho2, g, st)`
  ↔ `ch11.kh_min_shear`; `band(dU, …)` ↔ `ch11.kh_unstable_band`.
- **Views** (rows [1.15, 1]): 1. `flow` "Two streams" (row 0, flex 1.4): upper layer (orange tint) with arrows U₁, lower (blue
  tint) with U₂, the interface growing (linear mode, amplitude capped) over two wavelengths at the chosen λ, a frame toggle (lab /
  moving with (U₁ + U₂)/2, N19) — in the moving frame the arrows become ±ΔU/2 and the crests stand still when ρ₁ = ρ₂. 2. `plane`
  "Stability boundary" (row 0, flex 1): log k (10 … 3000 m⁻¹) vs ΔU (0 … 15 m/s); boundary ΔU_min(k) (purple), its minimum ◆
  (6.70 m/s at λ = 1.73 cm with surface tension), unstable region rose, the current point (k, ΔU); pointer: click/drag → set k
  and ΔU (inspector). 3. `cplane` "Wave speeds c₊, c₋" (row 1, `hidePortrait`): c_r horizontal (around the mean), c_i vertical; the
  two roots (rose and teal) with their paths as ΔU changes (faint), the collision point marked. A 'growth' sub-mode of view 2 shows
  kc_i(k) instead (chips "boundary / growth").
- **Controls (≤ 5):** `dU` "Shear $\Delta U=U_1-U_2$" 0 … 15 m/s, step 0.05, default 5 · `fluids` chips air/water (1.2, 1000),
  thermocline (ρ₁/ρ₂ = 0.998, ρ₂ = 1025), equal densities, upside down (1000 over 1.2) · `lam` "Wavelength λ" log 0.3 cm … 10 m,
  default 1.73 cm (k = 2π/λ) · `tension` toggle (default on, σ_s = 0.074 N/m; off for the thermocline) · `h` "Lower depth h"
  (optional) log 1 mm … ∞ (∞ = deep).
- **Transport:** `t` 0 … 5 e-folding times (or 2 wave periods when neutral), `end: 'hold'`; end card "amplitude ×148 after 5
  e-folding times (… s)" / "a neutral wave: travelled … wavelengths".
- **Presets:** "air over water 5 m/s" {fluids: air, dU: 5, lam: 0.0173} · "air over water 8 m/s" {dU: 8, lam: 0.0101} · "equal
  densities (vortex sheet)" {fluids: equal, dU: 1, tension: false} · "ocean thermocline" {fluids: thermocline, dU: 0.1, lam: 1.26,
  tension: false} · "Rayleigh–Taylor" {fluids: upside, dU: 0}.
- **Status:** "🌊 stable: gravity + surface tension beat shear (discriminant +0.0239 m²/s² at λ = 1.73 cm)" (5 m/s) · "📈 unstable
  for 0.71 < λ < 4.21 cm; fastest λ ≈ 1.01 cm, growth 76 s⁻¹" (8 m/s; numbers computed by `kh_unstable_band` and a scan) ·
  "🌀 vortex sheet: unstable at every λ, growth kΔU/2" · thermocline: "📈 unstable for λ < 1.60 m (k_c = 3.93 m⁻¹)".
- **Readouts:** "c₊" · "c₋" · "Growth kc_i" · "Critical k_c" · "ΔU_min".
- **Terms** (bar chart, the discriminant under the root of (11.18)/(N120)): gravity $\frac{\rho_2-\rho_1}{\rho_2+\rho_1}\frac gk$
  (blue), surface tension $\frac{\sigma_sk}{\rho_1+\rho_2}$ (teal), shear $-\frac{\rho_1\rho_2\Delta U^2}{(\rho_1+\rho_2)^2}$ (orange);
  total (purple outline) — e.g. air/water, λ = 1.73 cm, 5 m/s: 0.02689 + 0.02689 − 0.02993 = +0.02385 (stable); 8 m/s: shear
  −0.07662, total −0.02284 (unstable). Click a bar to isolate it.
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** · **status** · **terms** ·
  **inspector** (click a (k, ΔU) point: "discriminant = 0.02689 + 0.02689 − 0.07662 = −0.02284 < 0 ⇒ c = 0.00959 ± 0.1511i m/s,
  kc_i = 55.0 s⁻¹" — the builder computes the exact digits) · **modes** (boundary / growth).
- **Explain:**
  0. *What the views show* — "**Two streams**: the interface wave with its amplitude growing (rose) or steady (teal); arrows are the
     streams. **Stability boundary**: below the purple curve every wave is neutral; above it (rose) waves of that k grow. **Wave
     speeds** (hidden on phones): the two roots c₊, c₋ of (11.18) — they collide on the real axis exactly on the boundary."
  1. *The three terms under the square root* — "gravity (ρ₂ − ρ₁)g/((ρ₂ + ρ₁)k) = **…**; surface tension σ_s k/(ρ₁ + ρ₂) = **…**;
     shear −ρ₁ρ₂ΔU²/(ρ₁ + ρ₂)² = **…**; total **…** m²/s²" boxed; why: (11.18) with Exercise 11.1's σ_s term.
  2. *The wave speeds* — "mean (ρ₂U₂ + ρ₁U₁)/(ρ₂ + ρ₁) = **…**; c = mean ± √(total) = **…**, **…**" — "a negative total gives a
     complex pair c₊ = c₋* (N17)".
  3. *Growth* — "kc_i = **…** s⁻¹; e-folding time 1/(kc_i) = **…**".
  4. *The boundary for these fluids* — "k_c = g(ρ₂² − ρ₁²)/(ρ₁ρ₂ΔU²) = **…** (no tension); with tension the minimum wind is
     ΔU_min = √(2√(gΔρσ_s)(ρ₁ + ρ₂)/(ρ₁ρ₂)) = **6.70** m/s at λ = 2π√(σ_s/(gΔρ)) = **1.73** cm" boxed.
  5. *Depth* (when h is finite) — "coth kh = **…**: the lower layer feels the bottom only for λ ≳ 2πh".
  6. *The c-plane* (or a hint) — "the two dots: **…**, **…**".
  7. *At the current time* — "amplitude factor e^{kc_i t} = `Viz.live('amp')`".
  8. *Reading the current setting* — stable: "Restoring terms win. Real wind-wave generation starts at lower winds through turbulent
     pressure fluctuations (named in C02) — KH alone needs 6.7 m/s." · unstable: "Shear wins in the band …; the fastest wave will
     dominate the pattern you see." · vortex sheet: "Nothing restores the sheet: every wavelength grows like kΔU/2 (11.20) — the
     roll-up of Ch. 5." · thermocline: "A 0.2 % density step is held only by gravity: 10 cm/s of shear makes every wave shorter than
     1.6 m grow — billows in the thermocline." · RT: "Heavy on top with no shear: Ch. 7's (7.95) with ρ₁ > ρ₂ — growth without
     shear, cut off by surface tension at 1.73 cm."
- **Derivation tab:** **D02** (12 steps) `view: 'flow'` (the linearisation, quoted on the moving-frame picture); step 7 `watch` "the
  interface is moved to z = 0: the arrows stay, the wave stays small". **D03** (14 steps) `view: 'cplane'` on desktop, `'plane'` on
  phones; goal `set` {fluids: air, dU: 8}; step 9 (the quadratic) `live` "ρ₁(U₁ − c)² + ρ₂(U₂ − c)² = (g/k)(ρ₂ − ρ₁) with your
  numbers"; step 12 `live` "the discriminant = …"; step 13 `watch` "drag ΔU through 6.7 m/s: the dots collide". **D04** (9 steps)
  `view: 'plane'`; step 9 `live` "ΔU_min = 6.70 m/s at k* = 364 m⁻¹". Interpret: `s => "With ΔU = " + dU + " m/s the interface is "
  + verdict + "."`.
- **Code:**
  ```python
  k = 2*np.pi/{{lam}}                                    # wavenumber [1/m]
  terms = ch11.kh_discriminant_terms(k, {{U1}}, 0.0, {{r1}}, {{r2}}, surface_tension={{st}})
  print(terms)                                           # {{terms}}
  cp, cm = ch11.kh_phase_speed(k, {{U1}}, 0.0, {{r1}}, {{r2}}, surface_tension={{st}})
  print(cp, k*np.imag(cp))                               # c+ = {{cp}}, growth {{gr}} 1/s
  print(ch11.kh_min_shear({{r1}}, {{r2}}, surface_tension=0.074))   # 6.70 m/s at 1.73 cm
  ```
- **Walkthrough (7 steps):** 1. "Wind over water" — "Does a 5 m/s wind make ripples grow? Watch." `play: true` · 2. "Three terms" —
  "Under the square root of (11.18): gravity and surface tension (+) against shear (−). Their sum decides." `terms: true`, `derive:
  {id: 'D03', step: 12}` · 3. "Stable at 5 m/s" — "Sum +0.024 m²/s² at λ = 1.73 cm: two real speeds, a travelling wave." `set` {dU:
  5}, `readouts: ['cp', 'cm']` · 4. "Cross the boundary" — "At 8 m/s the shear bar outgrows the others: the sum turns negative,
  c₊ = c₋* — one wave grows." `set` {dU: 8}, `highlight: ['view:cplane']` · 5. "The lowest wind" — "The boundary's minimum: 6.70 m/s
  at λ = 1.73 cm — where gravity and surface tension are equal." `derive: {id: 'D04', step: 8}` · 6. "Without surface tension" —
  "Turn it off: the boundary falls for short waves — every shear is unstable to some wavelength (N18)." `set` {tension: false} · 7.
  "Your turn" — "Predict the thermocline's critical wavelength for 10 cm/s, then pick the preset and check (≈ 1.6 m)." `controls:
  ['dU', 'lam']`.
- **Equations:** `kin` "Kinematic condition" ref 'Eq. (11.9)' $-U_1\frac{\partial\zeta}{\partial x}+\frac{\partial\phi_1}{\partial z}=
  \frac{\partial\zeta}{\partial t}=-U_2\frac{\partial\zeta}{\partial x}+\frac{\partial\phi_2}{\partial z}$ on z = 0 · `dyn` "Dynamic
  condition" ref 'Eq. (11.13)' $\rho_1(\partial_t\phi_1+U_1\partial_x\phi_1+g\zeta)=\rho_2(\partial_t\phi_2+U_2\partial_x\phi_2+g\zeta)$ ·
  `disp` "Kelvin–Helmholtz" ref 'Eq. (11.18)' (full) live "c = … ± …" · `crit` "Instability" ref 'p. 480 (unnumbered)'
  $g(\rho_2^2-\rho_1^2)<k\rho_1\rho_2(U_2-U_1)^2$ live · `static` ref 'Eq. (11.19)' · `sheet` ref 'Eq. (11.20)' · `ex` "Depth and
  surface tension (Exercise 11.1)" the N120 formula · symbols.
- **Check yourself:** (1) "Why does the boundary have a minimum when surface tension is on?" — "Gravity's term falls like 1/k,
  tension's rises like k; their sum is smallest at k* = √(gΔρ/σ_s) where the two are equal." · (2) "With equal densities, can any
  wind be too weak to make waves grow?" — "No: with ρ₁ = ρ₂ the gravity term vanishes; any ΔU gives c_i = ΔU/2 (11.20)."
  `set {fluids: equal}` · (3) "How fast does the unstable wave travel in the frame moving with the mean speed (equal densities)?" —
  "Not at all: c_r = (U₁ + U₂)/2, so it is at rest in that frame (N19)." `set {frame: 'moving'}` · (4) "Why can a thermocline
  billow with only 10 cm/s of shear?" — "Its density step is 0.2 %, so the gravity term is tiny; k_c = g(ρ₂² − ρ₁²)/(ρ₁ρ₂ΔU²) ≈ 3.9
  m⁻¹."
- **Selftest parity rows:** `{name: 'c+ imag tiny', js: khRoots(1, 6, 0, 1, 3, 10, 0)[0][1], py: 'np.imag(ch11.kh_phase_speed(1.0, 6.0,
  0.0, 1.0, 3.0, g=10.0)[0])', rtol: 1e-12}` · `{name: 'shear term 8 m/s', js: khTerms(363.88, 8, 0, 1.2, 1000, 9.80665, 0.074).shear,
  py: 'ch11.kh_discriminant_terms(363.88, 8.0, 0.0, 1.2, 1000.0, surface_tension=0.074)["shear"]', rtol: 1e-12}` · `{name: 'dU min',
  js: minShear(1.2, 1000, 9.80665, 0.074).dU, py: 'ch11.kh_min_shear(1.2, 1000.0, surface_tension=0.074)["dU_min"]', rtol: 1e-12}` ·
  `{name: 'band upper k', js: band(8, 1.2, 1000, 9.80665, 0.074)[1], py: 'ch11.kh_unstable_band(8.0, 1.2, 1000.0, surface_tension=
  0.074)[1]', rtol: 1e-10}` · `{name: 'boundary at k=100', js: boundary(100, 1.2, 1000, 9.80665, 0.074, Infinity), py:
  'ch11.kh_stability_boundary(100.0, 1.2, 1000.0, surface_tension=0.074)', rtol: 1e-10}`.
- **Fit plan:** 360×640: status, `flow` (50 %) over `plane` (50 %), `cplane` hidden (c₊, c₋ in readouts and in the `plane` title);
  terms open as an overlay card from the status. 844×390: `flow` | `plane`. Desktop: rows [1.15, 1] with `cplane` under `plane`.

### E3 · benard_neutral_curve
- **Title:** "Why does convection start near Ra ≈ 1708?" · **Summary:** "Each cell width K has its own marginal Rayleigh number;
  drag a point in the (K, Ra) plane, watch the rolls grow or decay and the determinant's imaginary part cross zero, and switch the
  walls to move the whole valley." · **CORE:** C04, C03, C05 (also N25 with the Γ sign, N41–N46, N44, N50, N51, N123) ·
  **Reference:** `amplitude_phase_second_order_II_3.html` (numbered derivation with live numbers) + `angular_frequency_explorer_1.html`
  (linked views, presets).
- **meta:** `viz:order 3` · `viz:sections 11.4` · `viz:equations 11.21 11.36 11.37 11.39 11.40 11.41 11.42 11.43 11.44` ·
  `viz:fluidpy ch11.benard_char_roots ch11.benard_determinant ch11.benard_marginal_Ra_det ch11.benard_free_free_Ra
  ch11.benard_free_free_critical ch11.benard_neutral_table ch11.rayleigh_number ch11.gamma_conventions` · `viz:derivations D09 D10
  D11`.
- **Physics:** `charRoots(Ra, K)` ↔ `ch11.benard_char_roots`; `detEven(Ra, K)` / `detOdd(Ra, K)` with local complex cos/cosh (C.5
  5.3) ↔ `ch11.benard_determinant`; `margRaDet(K, mode)` (scan above K⁴ + brentq on Im det, `Viz.num.brentq`) ↔
  `ch11.benard_marginal_Ra_det`; `ffRa(K, n)` ↔ `ch11.benard_free_free_Ra`; rigid–free and odd curves from the table `NEUTRAL`
  (C.5 5.4) ↔ `ch11.benard_neutral_table`; `rayleigh(alpha, dT, d, kappa, nu, g)` ↔ `ch11.rayleigh_number`.
- **Views** (rows [1.1, 1]): 1. `layer` "The layer" (row 0, flex 1.3): roll streamlines over T′ (blue–rose) for the current K (even:
  one row; odd: two rows), cells of width π/K, amplitude ∝ e^{σt} with σ from the sign of Ra − Ra(K) (exact σ for free walls via
  D12, a qualitative rate elsewhere — labelled), plates hatched (rigid) or plain (free); title "cell width 2π/K = … d". 2. `curve`
  "Neutral curves Ra(K)" (row 0, flex 1): K horizontal (0.5 … 10), Ra vertical (log 400 … 2 × 10⁴); rigid–rigid (purple bold),
  rigid–free (purple dashed), free–free (muted), odd (dotted, toggle); minima ◆ labelled; the unstable region of the chosen boundary
  rose; the current point; pointer: drag → set (K, Ra), click on a curve → inspector. 3. `det` "Im det(Ra) at your K" (row 1,
  `hidePortrait`): the imaginary part of the 3 × 3 determinant against Ra (log), its first zero (purple ●) = the marginal Ra, the
  current Ra line; for free walls the panel shows the bars (π² + K²)³ and K² whose ratio is Ra (11.44).
- **Controls (≤ 5):** `walls` chips rigid–rigid / free–free / rigid–free (mode) · `K` "Wavenumber $K$" 0.5 … 10, step 0.01, default
  3.117 · `lRa` "$\log_{10}\mathrm{Ra}$" 2.6 … 4.5, step 0.005, default log10(1707.76) · `parity` chips even / odd (optional) · `fluid`
  (optional, Explain only) chips water 5 mm / air 1 cm with ΔT slider giving Ra through (11.21).
- **Transport:** `t` 0 … 3 (units d²/κ), `end: 'hold'`; end card "rolls grew ×… / decayed to …".
- **Presets:** "rigid onset" {walls: rigid, K: 3.1163, lRa: log10(1707.76)} · "free onset" {walls: free, K: 2.2214, lRa: log10(657.51)}
  · "rigid–free" {walls: rigidfree, K: 2.682, lRa: log10(1100.65)} · "Ra = 2500 band" {walls: rigid, lRa: log10(2500), K: 3.1} ·
  "odd mode" {parity: odd, K: 5.365, lRa: log10(17610.4)}.
- **Status:** "🔥 Ra = 2500 > Ra(K = 4) = 1879: rolls grow" / "❄️ Ra = 1500 < Ra(K) = 1708: conduction holds" / "⚖️ on the curve: marginal
  (σ = 0, stationary)"; exact Ra(K) computed live (rigid by the determinant; others from the table, flagged "table").
- **Readouts:** "Marginal Ra(K)" · "Ra / Ra(K)" · "Cell width λ/d" · "q₀" · "Im det".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (walls) · **presets** · **status** · **inspector**
  (click the rigid curve at K: "q₀ = K[(Ra/K⁴)^{1/3} − 1]^{1/2} = …; q = …; det = … i; Ra(K) = …") · **transport**.
- **Explain:**
  0. *What the views show* — "**The layer**: rolls drawn for your K (blue cold, rose warm); they brighten if your Ra is above the
     curve. **Neutral curves**: for each K the Ra where σ = 0 — rigid walls purple, free walls grey; inside a curve (rose) that K
     grows. **Im det** (hidden on phones): the book's 3 × 3 determinant along Ra — its first zero is the marginal Ra."
  1. *Your Ra from a real layer* (fluid chips) — "Ra = gαΔT d³/(κν) (11.21) = 9.81 × 2.1e-4 × ΔT × (0.005)³/(1.4e-7 × 1e-6) =
     **…**" with the two-convention line "Γ = −dT̄/dz = ΔT/d = **400** K/m in (11.21); Kundu Ch. 1: dT/dz = **−400** K/m;
     meteorology: −dT/dz = **400** K/m" boxed.
  2. *The roots at your (K, Ra)* — "s = (Ra/K⁴)^{1/3} = **…**; q₀ = K√(s − 1) = **…**; q = K√(1 + ½s(1 + i√3)) = **…**" ((11.42)).
  3. *The determinant* — the matrix with numbers (rows: W, W′, (D² − K²)²W at z = ½), "det = **…** i (purely imaginary — the
     second and third columns are conjugates)".
  4. *Marginal Ra and the distance from onset* — "Ra(K) = **…**; your Ra/Ra(K) = **…** ⇒ grows/decays; cell width 2π/K = **…** d".
  5. *Free walls* — "Ra = (π² + K²)³/K² (11.44) = **…**; minimum at K² = π²/2: 27π⁴/4 = **657.51**".
  6. *Im det panel* (or the hint) — "first zero at Ra = **…**".
  7. *Reading the current setting* — below: "Heating too weak for this cell width: buoyancy loses to viscosity and diffusion." ·
     above: "This K grows. The band of growing K at your Ra is **…**–**…** (rigid: 1.76–5.08 at Ra = 2500)." · at the minimum: "The
     cheapest cell: K_c = 3.117, cells ≈ 2.0 d wide — what you see first." · odd: "Two rows of cells need ten times more heating
     (17 610)." · free: "Stress-free walls resist less: 657.5."
- **Derivation tab:** **D09** (6 steps) `view: 'curve'`. **D10** (★★★, 14 steps) `view: 'det'` on desktop, `'curve'` on phones;
  goal `set` {walls: rigid, K: 2}; step 1 `live` "(q² − K²)³ = −RaK² = −… "; step 5 `live` "q₀ = 3.4459 at K = 2, Ra = 1000"; step 9
  quotes slip #1 (`watch`: "the book's B in the third term is a C"); step 13 `live` "Im det changes sign at Ra = 2177.41"; step 14
  `set` {K: 3.1163, lRa: log10(1707.76)} `watch` "the ◆ is the minimum". **D11** (11 steps) `view: 'curve'`, `set` {walls: free};
  step 7 (slip #2) and step 9 (slip #3) quote their boxes. Interpret: `s => "At K = " + K + " convection needs Ra = " + RaK + "."`.
- **Code:**
  ```python
  print(ch11.benard_char_roots({{Ra}}, {{K}}))               # q0, q, q* = {{roots}}
  print(ch11.benard_determinant({{Ra}}, {{K}}))              # {{det}} (purely imaginary)
  RaK = ch11.benard_marginal_Ra_det({{K}})                   # marginal Ra at your K = {{RaK}}
  print(ch11.benard_critical(bc=("{{b0}}", "{{b1}}")))       # {{crit}}
  print(ch11.benard_free_free_Ra({{K}}))                     # (11.44): {{ffRa}}
  ```
- **Walkthrough (7 steps):** 1. "A heated layer" — "Warm below, cool above: when does it start to turn over?" `play: true` · 2. "One
  curve, many cell widths" — "Each K has its own marginal Ra. Drag K along the purple curve." `controls: ['K']` · 3. "Six roots" —
  "(q² − K²)³ = −RaK²: one imaginary pair ±iq₀ and two complex pairs (11.42)." `derive: {id: 'D10', step: 5}` · 4. "The
  determinant" — "Three wall conditions on A cos q₀z + B cosh qz + C cosh q*z: nonzero only if the 3 × 3 determinant vanishes."
  `derive: {id: 'D10', step: 11}`, `highlight: ['view:det']` · 5. "The bottom of the valley" — "Minimise over K: Ra_c = 1707.76 at
  K_c = 3.117 — cells about 2d wide." `set` {K: 3.1163, lRa: 3.2324}, `readouts: ['RaK', 'width']` · 6. "Change the walls" — "Free
  walls: 657.5 at π/√2. One free wall: 1100.65 at 2.682." `set` {walls: free} · 7. "Your turn" — "Predict the band of growing K at
  Ra = 2500 with rigid walls; then drag K to find its edges (1.76 and 5.08)." `controls: ['K', 'lRa']`.
- **Equations:** `Ra` ref 'Eq. (11.21)' $\mathrm{Ra}=g\alpha\Gamma d^4/\kappa\nu$ live · `amp` ref 'Eqs. (11.36, 11.37)' (both) ·
  `marg` ref 'Eq. (11.39)' · `six` ref 'Eq. (11.40)' $(\frac{d^2}{dz^2}-K^2)^3W=-\mathrm{Ra}K^2W$ · `bc` ref 'Eq. (11.41)' · `roots`
  ref 'Eq. (11.42)' live · `ff` ref 'Eq. (11.44)' live · symbols K, Ra, q₀, q (–), d (m), ΔT (K).
- **Check yourself:** (1) "Why does the curve rise for very narrow cells (large K)?" — "Narrow cells have steep gradients: viscosity
  and heat diffusion act over short distances and win." `set {K: 8}` · (2) "Doubling the layer depth: what happens to the ΔT
  needed?" — "Ra ∝ ΔT d³, so ΔT_c falls by 8." · (3) "Why is the determinant purely imaginary?" — "Columns 2 and 3 are complex
  conjugates (q and q*): swapping them conjugates the matrix, so det* = −det." · (4) "Which onset is lower, rigid–free or free–free?
  Why?" — "Free–free (657.5 < 1100.65): each stress-free wall removes a no-slip constraint." `set {walls: rigidfree}`.
- **Selftest parity rows:** `{name: 'q0 K=2 Ra=1000', js: charRoots(1000, 2)[0], py: 'ch11.benard_char_roots(1000.0, 2.0)[0]', rtol:
  1e-12}` · `{name: 'Im det at Ra_c', js: detEven(1707.762, 3.1163)[1], py: 'np.imag(ch11.benard_determinant(1707.762, 3.1163))', rtol:
  1e-8}` · `{name: 'marginal Ra K=2 (det)', js: margRaDet(2, 'even'), py: 'ch11.benard_marginal_Ra(2.0)', rtol: 1e-7}` · `{name: 'ff Ra
  K=2', js: ffRa(2, 1), py: 'ch11.benard_free_free_Ra(2.0)', rtol: 1e-12}` · `{name: 'Ra water', js: rayleigh(2.1e-4, 2, 0.005,
  1.4e-7, 1e-6, 9.81), py: 'ch11.rayleigh_number(alpha=2.1e-4, dT=2.0, d=0.005, kappa=1.4e-7, nu=1e-6, g=9.81)', rtol: 1e-12}` ·
  `{name: 'table rigid-free min', js: Math.min(...NEUTRAL.rigid_free), py: 'np.min(ch11.benard_neutral_table()["rigid_free"])', rtol:
  1e-3}`.
- **Fit plan:** 360×640: status, `layer` (45 %) over `curve` (55 %), `det` hidden (Im det and marginal Ra in readouts); 844×390:
  `layer` | `curve`; desktop: rows [1.1, 1] with `det` under `curve`.

### E4 · salt_fingers
- **Title:** "How can a stably stratified column overturn?" · **Summary:** "Heat diffuses about a hundred times faster than salt: a
  displaced parcel keeps its salt and keeps sinking. Drag the two gradients across a regime map and see a column that is stable by
  density yet finger-unstable." · **CORE:** C06 (also N53–N57, N54 the diffusive regime, N55 the second Ra sign) · **Reference:**
  `forced_damped_vibrations.html` (regime status + interpretation) + `overfitting_curves.html` (a minimal two-slider map with a
  verdict line).
- **meta:** `viz:order 4` · `viz:sections 11.5` · `viz:equations 11.45 11.46` · `viz:fluidpy ch11.salt_finger_unstable
  ch11.double_diffusive_margin ch11.double_diffusive_sigma ch11.salt_finger_regime ch11.thermal_rayleigh_signed
  ch11.salinity_rayleigh ch11.linear_eos` · `viz:derivations D13`.
- **Physics:** `RaS(dTdz, d)` ↔ `ch11.thermal_rayleigh_signed`; `RsS(dSdz, d, ks)` ↔ `ch11.salinity_rayleigh`; `finger(dTdz, dSdz,
  d, …)` ↔ `ch11.salt_finger_unstable`; `cubicRoots(K2, Ra, Rs, Pr, tau)` (closed-form cubic via a local Cardano/`Viz.num`-free
  routine or Durand–Kerner) ↔ `ch11.double_diffusive_sigma`; `regime(…)` ↔ `ch11.salt_finger_regime` (exact `text` pinned).
  Fixed: α = 2 × 10⁻⁴ K⁻¹, β = 7.6 × 10⁻⁴ (g/kg)⁻¹, ν = 10⁻⁶, κ = 1.4 × 10⁻⁷ m²/s, g = 9.80665 m/s².
- **Views** (rows [1.1, 1]): 1. `column` "A displaced parcel" (row 0, flex 1): a vertical column with T (rose–blue) and S (amber)
  colour strips, a parcel pushed down (or up) by the transport: its inside T relaxes toward the surroundings with rate ∝ κ, its S
  with rate ∝ κ_s (two small bars inside the parcel: ΔT and ΔS against the surroundings, and its buoyancy arrow, blue up / rose
  down). 2. `map` "Regime map" (row 0, flex 1.2): x = α dT̄/dz, y = β dS/dz [10⁻⁶ m⁻¹]; the static-stability line αT_z = βS_z
  (blue), the finger line (11.46) (amber), regions stable (teal), fingers (amber–rose), diffusive (rose hatched), overturning (rose);
  the current point; pointer: drag → set the gradients; click → inspector. 3. `roots` "σ of the free–free cubic" (row 1,
  `hidePortrait`): the three roots in the complex plane, the imaginary axis purple.
- **Controls (≤ 5):** `dTdz` "$d\bar T/dz$" −0.02 … 0.02 K/m, default 0.01 · `dSdz` "$dS/dz$" −0.005 … 0.005 (g/kg)/m, default 0.002 ·
  `tau` "$\kappa_s/\kappa$" log 0.003 … 1, default 0.0107 · `d` "Layer depth d" log 0.005 … 0.3 m, default 0.05 · `Pr` (optional) 1 …
  20, default 7.
- **Transport:** parcel displacement and relaxation over 3 diffusion times of the parcel, `end: 'hold'`; end card "after the heat
  has gone: ΔT ≈ 0, ΔS = 98 % of the start — the parcel is heavier and keeps sinking" (finger case).
- **Presets:** "subtropical thermocline (fingers)" {dTdz: 0.01, dSdz: 0.002, d: 0.05} · "Arctic: cold fresh over warm salty
  (diffusive)" {dTdz: −0.01, dSdz: −0.0027, d: 0.11} · "κ_s = κ: no fingers" {tau: 1} · "statically unstable" {dTdz: −0.01, dSdz: 0}
  · "thinnest finger layer" {dTdz: 0.01, dSdz: 0.002, d: 0.0161}.
- **Status:** "🧂 statically stable (R_ρ = 1.32) but finger-unstable: Rs − Ra = 61 254 > 657.5" (thermocline) · "〰️ diffusive regime:
  statically stable, oscillatory growth σ = 9.04 ± 72.0 i" (Arctic preset: Ra ≈ −20 500, Rs ≈ −1.96 × 10⁶) · "✅ stable: Rs − Ra =
  … < 657.5" · "⚠️ statically unstable: heavy on top — ordinary convection". Exact strings from `salt_finger_regime(...)["text"]`.
- **Readouts:** "Ra (§11.5 sign)" · "Rs" · "Rs − Ra − 27π⁴/4" · "R_ρ = αT_z/(βS_z)" · "max Re σ".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **presets** · **status** · **transport** · **inspector**
  (click the map: "Rs − Ra = (gd⁴/ν)[βS_z/κ_s − αT_z/κ] = (9.81 × 6.25e-6/1e-6)(1013.3 − 14.3) = 61 254 > 657.5").
- **Explain:**
  0. *What the views show* — "**Parcel**: its two bars show how much warmer/saltier it is than its new neighbours; heat leaks fast,
     salt slowly. **Regime map**: blue line = neutral density, amber line = (11.46); the dot is your column. **σ roots** (hidden on
     phones): growth rates of the free–free cubic (ours)."
  1. *Density* — "dρ̄/dz ∝ −αT_z + βS_z = −**2.0** + **1.52** (×10⁻⁶ m⁻¹) < 0: lighter on top — statically stable; density ratio
     R_ρ = αT_z/(βS_z) = **1.32**."
  2. *The two Rayleigh numbers (§11.5 signs)* — "Ra = gαd⁴(dT̄/dz)/(νκ) = **875.9** (positive here: warm on top stabilises); Rs =
     gβd⁴(dS/dz)/(νκ_s) = **62 130**" with the ⚠️ line "§11.5's Ra has the opposite sign to (11.21)'s; heated from below it would be
     negative".
  3. *The finger criterion (11.46)* — "Rs − Ra = **61 254** vs 27π⁴/4 = **657.5** ⇒ margin **60 597** > 0: fingers" boxed; "the salt
     term is divided by κ_s, the heat term by κ: salt counts κ/κ_s = **93** times more than in the density".
  4. *Thinnest unstable layer* — "d_min = [657.5ν/(g(βS_z/κ_s − αT_z/κ))]^{1/4} = **1.61** cm".
  5. *The σ roots* (or the hint) — "**…**, **…**, **…** (one real positive: stationary fingers / a complex pair: oscillatory)".
  6. *Reading the current setting* — fingers: "Warm salty water over cold fresh water, stable by density, overturns in thin fingers
     — the subtropical North Atlantic's thermohaline staircases." · diffusive: "Cold fresh over warm salty (under Arctic sea ice): a
     displaced parcel loses heat and overshoots — growing oscillations that build convecting layers." · κ_s = κ: "With equal
     diffusivities the finger line lies on the density line: only a top-heavy column can overturn — diffusion stabilises a
     single-component fluid." · overturning: "Top-heavy: ordinary convection (C03–C05)."
- **Derivation tab:** **D13** (12 steps) `view: 'map'`; goal `set` thermocline preset; step 5 (the second Ra sign) `watch` "§11.5's Ra is
  negative when heated from below"; step 11 `live` "T̂ = κ_sŝ/κ ⇒ (Rs − Ra)K²T̂ with your numbers"; step 12 `live` "Rs − Ra = 61 254 >
  657.5". Interpret: `s => "Your column is " + regime + "."`.
- **Code:**
  ```python
  r = ch11.salt_finger_unstable(dTdz={{Tz}}, dSdz={{Sz}}, d={{d}}, kappa_s={{ks}})
  print(r["density_stable"], r["R_rho"])        # {{ds}}, {{Rrho}}
  print(r["Ra"], r["Rs"], r["margin"])            # {{Ra}}, {{Rs}}, margin {{m}}
  sig = ch11.double_diffusive_sigma(np.pi**2/2, r["Ra"], r["Rs"], {{Pr}}, {{tau}})
  print(sig)                                      # {{sig}}
  ```
- **Walkthrough (6 steps):** 1. "Lighter on top" — "Warm salty water over cold fresh water: stable by density. Can it overturn?"
  `play: true` · 2. "A parcel's two clocks" — "Push a parcel down: heat equalises in seconds, salt stays. It is now heavier — and
  sinks on." `highlight: ['view:column']` · 3. "Two Rayleigh numbers" — "Rs counts salt with its slow diffusivity; (11.46): fingers
  when Rs − Ra > 27π⁴/4 ≈ 657.5." `derive: {id: 'D13', step: 12}`, `readouts: ['margin']` · 4. "The wedge" — "Between the blue density
  line and the amber finger line: stable by density, unstable by diffusion." `highlight: ['view:map']` · 5. "Close the wedge" — "Set
  κ_s = κ: the lines merge; no double diffusion." `set` {tau: 1} · 6. "Your turn" — "Predict the thinnest finger-unstable layer for
  the thermocline gradients (≈ 1.6 cm); then shrink d and check." `controls: ['d']`.
- **Equations:** `eos` "Linear EOS" ref 'p. 492 (unnumbered)' $\tilde\rho=\rho_0[1-\alpha(\tilde T-T_0)+\beta(\tilde s-s_0)]$ · `marg`
  "Double-diffusive marginal equations" ref 'Eq. (11.45)' (all three) · `crit` "Salt fingers" ref 'Eq. (11.46)' live · `rs` "Rs − Ra"
  ref 'p. 494 (unnumbered)' $\mathrm{Rs}-\mathrm{Ra}=\tfrac{27}4\pi^4=657$ · `cubic` "Free–free cubic (ours)" ref 'no book number'
  · symbols.
- **Check yourself:** (1) "Why does the salt gradient count ~93 times more in (11.46) than in the density?" — "It is divided by κ_s,
  the heat gradient by κ, and κ/κ_s ≈ 93." · (2) "Does a layer get more or less finger-unstable if it is thicker?" — "More: the
  bracket is multiplied by d⁴." `set {d: 0.2}` · (3) "In the Arctic preset, why does the parcel oscillate?" — "It loses its heat
  quickly but keeps its (low) salt; stratification pushes it back and it overshoots." `set` Arctic · (4) "What does κ_s = κ do?" —
  "Removes the difference in clocks: the finger line becomes the density line." `set {tau: 1}`.
- **Selftest parity rows:** `{name: 'lhs 11.46', js: finger(0.01, 0.002, 0.05).lhs, py: 'ch11.salt_finger_unstable(dTdz=0.01, dSdz=0.002,
  d=0.05)["lhs"]', rtol: 1e-12}` · `{name: 'Rs', js: RsS(0.002, 0.05, 1.5e-9), py: 'ch11.salinity_rayleigh(0.002, 0.05, 7.6e-4, 1e-6,
  1.5e-9)', rtol: 1e-12}` · `{name: 'finger root', js: maxRe(cubicRoots(Math.PI**2/2, 1000, 2000, 7, 0.0107)), py:
  'np.max(np.real(ch11.double_diffusive_sigma(np.pi**2/2, 1000.0, 2000.0, 7.0, 0.0107)))', rtol: 1e-8}` · `{name: 'diffusive
  imag', js: maxIm(cubicRoots(Math.PI**2/2, -2.0518e4, -1.9648e6, 7, 1.5e-9/1.4e-7)), py: 'np.max(np.imag(ch11.double_diffusive_sigma(
  np.pi**2/2, -2.0518e4, -1.9648e6, 7.0, 1.5e-9/1.4e-7)))', rtol: 1e-6}` · exact-text row for `regime(...).text` vs
  `ch11.salt_finger_regime(0.01, 0.002, 0.05)["text"]` (string equality as an `expect` row).
- **Fit plan:** 360×640: status, `column` (40 %) over `map` (60 %), `roots` hidden (max Re σ in readouts); 844×390: `column` | `map`;
  desktop rows [1.1, 1] with `roots` under `map`.

### E5 · taylor_couette_onset
- **Title:** "Why do vortices stack up between spinning cylinders?" · **Summary:** "Swapping two fluid rings releases energy when Γ²
  falls outward; drag (Ω₂/Ω₁, Ta) across the Rayleigh line and the viscous boundary, and watch Taylor vortices appear only past the
  second." · **CORE:** C07 (also R14, N59, R16, N64–N68, N124) · **Reference:** `angular_frequency_explorer_1.html` (system animation
  + response curve + presets).
- **meta:** `viz:order 5` · `viz:sections 11.6` · `viz:equations 11.49 11.51 11.52 11.53 11.54` · `viz:fluidpy ch11.taylor_number
  ch11.taylor_critical_approx ch11.taylor_critical_table ch11.rayleigh_circulation_criterion ch11.ring_interchange_energy
  ch11.couette_rayleigh_line ch11.taylor_eigenfunction` · `viz:derivations D14 D15`.
- **Physics:** `couette(R, R1, R2, O1, O2)` ↔ `LAM.circular_couette`; `taylorNumber(O1, O2, R1, R2, nu)` ↔ `ch11.taylor_number`;
  `TaApprox(mu)` ↔ `ch11.taylor_critical_approx`; `TaExact(mu)` = interpolation in the table `TAYLOR` (C.5 5.5) ↔
  `ch11.taylor_critical`; `ringDE(G1, G2, r1, r2)` ↔ `ch11.ring_interchange_energy`; vortex picture from `TAYLOR.psi` (μ = 0, 0.5,
  1; nearest shown, labelled). Dimensional frame: R₁ = 0.1 m, R₂ = 0.105 m (R₂/R₁ = 1.05, ours), ν = 10⁻⁶ m²/s.
- **Views** (rows [1.15, 1]): 1. `gap` "The gap" (row 0, flex 1.2): left, a top view of the two cylinders turning (arrows ∝ Ω₁,
  Ω₂); right, a meridional cut of the gap with Taylor vortices (streamlines of the eigenfunction, û_φ colour orange) whose
  amplitude grows if Ta > Ta_c, else a plain Couette shear (muted arrows). 2. `plane` "(μ, Ta) stability plane" (row 0, flex 1):
  μ = Ω₂/Ω₁ from −1 to 1.2, Ta log 10² … 10⁵; the exact narrow-gap Ta_c(μ) (purple dots/line), (11.54) (muted dashed), the Rayleigh
  line μ = (R₁/R₂)² = 0.907 (blue; to its right Rayleigh-stable), regions shaded; the current point; drag → set (μ, Ta). 3. `profile`
  "U_φ(R) and Γ²(R)" (row 1, `hidePortrait`): U_φ (orange) and Γ² = (2πRU_φ)² (purple) across the gap; two draggable ring markers r₁,
  r₂ and the ΔE bars (E_i, E_f, ΔE).
- **Controls (≤ 5):** `mu` "$\mu=\Omega_2/\Omega_1$" −1 … 1.2, step 0.01, default 0 · `lTa` "$\log_{10}\mathrm{Ta}$" 2 … 5, default
  log10(6098) · `O1` "$\Omega_1$" (optional; sets Ta through (11.52) for our cylinders) 0.05 … 2 rad/s · `rings` (inspector, two
  markers) · `fluid` chips water / glycerol-water (ν = 10⁻⁶, 10⁻⁵) (optional).
- **Transport:** `t` vortex growth over 4 growth times (σ from the sign of Ta − Ta_c, qualitative rate labelled), `end: 'hold'`; end
  card "Ta/Ta_c = 1.80: vortices" / "below onset: Couette flow survives".
- **Presets:** "inner only" {mu: 0, lTa: log10(6098)} · "co-rotation μ = 0.5" {mu: 0.5, lTa: log10(2736)} · "counter-rotation,
  viscously held" {mu: −0.5, lTa: log10(3405)} · "on the Rayleigh line" {mu: 0.907} · "μ → 1 = Bénard" {mu: 0.999, lTa:
  log10(1708.6)}.
- **Status:** "🌀 Taylor vortices: Ta = 6098 > Ta_c = 3390 (exact; (11.54) 3416)" · "🧲 Rayleigh-unstable but viscously stable: Ta =
  3405 < Ta_c = 6414" (counter-rotation preset) · "✅ Rayleigh-stable (μ > 0.907): no axisymmetric instability" · "⚖️ marginal".
- **Readouts:** "Ta" · "Ta_c exact" · "Ta_c (11.54)" · "error of (11.54)" · "ΔE rings".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **presets** · **status** · **inspector** (two rings: "Γ₁ =
  2πr₁U_φ(r₁) = 0.02526, Γ₂ = 0.00641 m²/s; ΔE = (Γ₂² − Γ₁²)(1/r₁² − 1/r₂²)/(8π²) = −4.21 × 10⁻⁵ m²/s²: released" for r₁ = 0.101,
  r₂ = 0.104 m, inner-only) · **transport** · **modes** (narrow-gap theory / Galerkin route of N124 shown as a second dot set,
  optional).
- **Explain:**
  0. *What the views show* — "**The gap**: top view of the cylinders and a cut through the gap; vortices appear when your point is
     above the purple curve. **(μ, Ta) plane**: purple = exact narrow-gap onset, dashed = (11.54), blue = Rayleigh's line. **Profiles**
     (hidden on phones): the base flow U_φ and Γ² across the gap, with the two rings you can swap."
  1. *The base flow* — "A = (Ω₂R₂² − Ω₁R₁²)/(R₂² − R₁²) = **…** s⁻¹, B = (Ω₁ − Ω₂)R₁²R₂²/(R₂² − R₁²) = **…** m²/s; U_φ = AR + B/R
     (11.49)".
  2. *Rayleigh's verdict* — "dΓ²/dr < 0 at **…** of the gap; Rayleigh line μ = (R₁/R₂)² = **0.907**".
  3. *The Taylor number (11.52)* — "Ta = 4(Ω₁R₁² − Ω₂R₂²)/(R₂² − R₁²) · Ω₁d⁴/ν² = **…**" boxed (with the narrow inner-only shortcut
     when μ = 0).
  4. *The threshold* — "(11.54): Ta_c = 1708/(½(1 + μ)) = **…**; exact narrow gap (our table): **…** at k = **…**; error **…** %".
  5. *Ring swap* (or the hint) — the inspector numbers.
  6. *Reading the current setting* — inner only: "Γ² falls outward everywhere: rings want to swap; viscosity resists until Ta ≈
     3390." · co-rotation: "Less shear, a higher threshold; (11.54) is within 0.1 % here." · counter-rotation: "Γ² falls outward near
     the inner wall only; viscosity holds the flow until Ta ≈ 6414 — and (11.54) overestimates by 6.5 % (σ = 0 is only assumed
     here)." · Rayleigh-stable: "Γ² grows outward: swapping costs energy. No axisymmetric instability at any speed." · μ → 1: "The
     equations become Bénard's (D15): Ta_c → 1707.76, cells as wide as the gap." Climate line: "Earth's rotation turns this into
     inertial instability of jets (Ch. 13)."
- **Derivation tab:** **D14** (★★★, 15 steps) `view: 'profile'` on desktop, `'plane'` on phones; step 4 quotes slip #12; step 9
  `live` "dU_φ/dR + U_φ/R = 2A = **…** s⁻¹ for your cylinders"; step 12 `watch` "the (1 + αx) factor with α = μ − 1 = **…**";
  step 15 `live` "Ta = −4AΩ₁d⁴/ν² = **…**". **D15** (6 steps) `view: 'plane'`; step 6 `set` {mu: 0.999} `watch` "the purple curve meets
  1708 at μ = 1". Interpret: `s => verdict`.
- **Code:**
  ```python
  r = ch11.taylor_number({{O1}}, {{O2}}, 0.1, 0.105, {{nu}})    # Ta = {{Ta}}, Rayleigh-stable: {{rs}}
  print(ch11.taylor_critical({{mu}}))                          # exact narrow gap: {{Tac}}
  print(ch11.taylor_critical_approx({{mu}}))                   # (11.54): {{Taa}}
  print(ch11.couette_rayleigh_line(0.1, 0.105))                # mu on Rayleigh's line: 0.907
  print(ch11.ring_interchange_energy({{G1}}, {{G2}}, {{r1}}, {{r2}}))   # dE = {{dE}}
  ```
- **Walkthrough (6 steps):** 1. "Spin the inner cylinder" — "Water between two cylinders; the inner one turns. When do vortices
  appear?" `play: true` · 2. "Swap two rings" — "Each ring keeps its circulation; the energy change has the sign of Γ₂² − Γ₁². Inner
  spinning: energy released." `inspect: true` · 3. "The Taylor number" — "Ta (11.52) compares centrifugal driving with viscosity: 6098
  for Ω₁ = 0.5 rad/s." `readouts: ['Ta']`, `derive: {id: 'D14', step: 15}` · 4. "The threshold" — "Exact narrow gap 3390 at k = 3.13;
  (11.54) gives 3416. Above it, vortices about 2d tall." `set` {mu: 0} · 5. "Counter-rotation" — "Rayleigh says unstable near the
  inner wall, yet viscosity holds it up to Ta ≈ 6414." `set` counter preset · 6. "Your turn" — "Predict Ta_c at μ = 1, then drag
  there: it is Bénard's 1708 (D15)." `controls: ['mu', 'lTa']`.
- **Equations:** `base` ref 'Eq. (11.49)' · `ray` "Rayleigh" ref 'p. 497 (unnumbered)' $d\Gamma^2/dr<0$ · `narrow` ref 'Eq. (11.51)'
  (both lines) · `Ta` ref 'Eq. (11.52)' live · `bc` ref 'Eq. (11.53)' · `crit` ref 'Eq. (11.54)' live · symbols Ω (rad/s), R, d (m), ν
  (m²/s), Γ (m²/s), Ta, μ, k (–).
- **Check yourself:** (1) "Why is the outer-only case stable at any speed?" — "Γ² = (2πRU_φ)² grows outward; swapping rings costs
  energy." `set {mu: 1.1}` · (2) "Where does (11.54) fail most?" — "Counter-rotation: 6.5 % high at
  μ = −0.5, where σ = 0 is only assumed." · (3) "What is the vortex height at onset?" — "Half a wavelength π/k_c ≈ 1.0 gap: square
  cells." · (4) "Why does μ → 1 give 1708?" — "With α = 0 the narrow-gap equations are exactly the Bénard pair with Ta for Ra (D15)."
- **Selftest parity rows:** `{name: 'Ta inner', js: taylorNumber(0.5, 0, 0.1, 0.105, 1e-6), py: 'ch11.taylor_number(0.5, 0.0, 0.1,
  0.105, 1e-6)["Ta"]', rtol: 1e-12}` · `{name: 'Ta_c approx mu=0.5', js: TaApprox(0.5), py: 'ch11.taylor_critical_approx(0.5)', rtol:
  1e-14}` · `{name: 'Ta_c exact mu=0 (table)', js: TaExact(0), py: 'ch11.taylor_critical(0.0)["Ta_c"]', rtol: 1e-3}` · `{name: 'ring dE',
  js: ringDE(4, 2, 1, 2), py: 'ch11.ring_interchange_energy(4.0, 2.0, 1.0, 2.0)["dE"]', rtol: 1e-12}` · `{name: 'Rayleigh line', js:
  Math.pow(0.1/0.105, 2), py: 'ch11.couette_rayleigh_line(0.1, 0.105)', rtol: 1e-14}`.
- **Fit plan:** 360×640: status, `gap` (50 %) over `plane` (50 %), `profile` hidden (ΔE in readouts); 844×390: `gap` | `plane`;
  desktop rows [1.15, 1].

### E6 · richardson_shear_instability
- **Title:** "Why does Ri = ¼ decide whether shear makes billows?" · **Summary:** "Slide the stratification J of a tanh shear layer:
  the Ri(z) profile dips below ¼ exactly when the growth tongue of the (k, J) map is reached and the billows start to grow." ·
  **CORE:** C08, C09 (also R18, R20, N73, N78, N20) · **Reference:** `forced_damped_vibrations.html` (live Explain) +
  `angular_frequency_explorer_1.html` (linked views on one clock).
- **meta:** `viz:order 6` · `viz:sections 11.7` · `viz:equations 11.56 11.61 11.63 11.64 11.65 11.66 11.67` · `viz:fluidpy
  ch11.gradient_richardson ch11.miles_howard_stable ch11.richardson_profiles ch11.tg_growth_map ch11.tg_growth` · `viz:derivations
  D17 D18`.
- **Physics:** `profiles(J, R)` (U = tanh z, U′ = sech²z, N² = J sech²(Rz)) ↔ `ch11.richardson_profiles`; `Ri(z, J, R)` ↔
  `ch11.gradient_richardson`; `mh(J, R)` (Ri_min and where, by a fine scan) ↔ `ch11.miles_howard_stable`; growth kc_i(k, J) bilinear
  in the table `TGMAP` (C.5 5.6; R = 1 only; for R ≠ 1 the map is hidden and the status uses Ri only) ↔ `ch11.tg_growth`. Billow
  picture: the linear mode's streamlines are not available live; the scene draws density layers displaced by a sinusoidal
  displacement with amplitude e^{kc_i t} (label "linear growth, schematic shape").
- **Views** (rows [1.1, 1]): 1. `billow` "The shear layer" (row 0, flex 1.3): coloured density layers (blue shades), velocity arrows
  (orange) ±1 far away, the layers waving with amplitude ∝ e^{kc_i t} at the chosen k (capped), the billow-cloud preset drawn as a
  cloud band. 2. `prof` "U, N², Ri" (row 0, flex 1): z vertical (−4 … 4); U(z) orange, N²(z) blue, Ri(z) purple on a log x-axis with
  the ¼ line (rose dashed), Ri_min ●; click → inspector. 3. `map` "Growth map kc_i(k, J)" (row 1, `hidePortrait`): k 0.05 … 1, J 0 …
  0.3; colour = kc_i (white = 0), the J = ¼ line, the current (k, J) point; click → set.
- **Controls (≤ 5):** `J` "Bulk Richardson $J$ (= Ri at the centre)" 0 … 0.4, step 0.005, default 0.1 · `k` "Wavenumber $k$" 0.05 …
  1, default 0.44 · `R` chips "density layer as thick / half / a third" (R = 1, 2, 3) · transport `t` · `dims` (optional) chips
  "non-dimensional / thermocline (ΔU = 0.2 m/s, L = 1 m) / atmosphere (ΔU = 10 m/s, L = 100 m)" for the e-folding time in seconds.
- **Transport:** `t` over 4 e-folding times of the current mode (or a fixed 50 time units when stable), `end: 'hold'`; end card
  "amplitude ×55" / "no growth: Ri > ¼ everywhere".
- **Presets:** "J = 0: pure shear layer" {J: 0, k: 0.445} · "J = 0.1" {J: 0.1, k: 0.44} · "J = 0.24" {J: 0.24} · "J = 0.26"
  {J: 0.26} · "thermocline" {J: 0.08, dims: thermo} · "billow-cloud layer" {J: 0.05, dims: atmo}.
- **Status:** "🌀 Ri_min = 0.10 < ¼: instability allowed; here kc_i = 0.124 at k = 0.40 (table)" · "🛡️ Ri > ¼ everywhere (Ri_min =
  0.26): stable by Miles–Howard (11.67) — for every k" · "⚖️ Ri_min < ¼ but this k is outside the tongue: allowed, not happening".
- **Readouts:** "Ri_min" · "where" · "kc_i" · "e-folding time" · "verdict".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **presets** · **status** · **transport** · **inspector**
  (click Ri(z): "N² = J sech²(Rz) = …, U′ = sech²z = …, Ri = … at z = …").
- **Explain:**
  0. *What the views show* — "**Shear layer**: density layers (blue) carried by the tanh current (orange); they wave and grow if the
     current (k, J) is unstable. **U, N², Ri**: the profiles and the Richardson number, with ¼ dashed. **Growth map** (hidden on
     phones): computed kc_i for every (k, J) — ours, from `ch11.tg_growth_map`."
  1. *At the centre* — "N² = J = **0.10**, U′ = 1 ⇒ Ri = N²/U′² (11.66) = **0.10**."
  2. *The minimum* — "Ri(z) = N²/U′² = J sech²(Rz)/sech⁴z = J sech²(Rz) cosh⁴z (= J cosh²z for R = 1) = **…**; Ri_min = **…** at
     z = **…**" boxed.
  3. *The guarantee* — "Ri_min > ¼? **no/yes** ⇒ (11.67) **does not apply / guarantees stability for every k**".
  4. *What the solver finds* — "kc_i = **…** at k = **…** (table, approximate near the edge); e-folding 1/(kc_i) = **…** L/U₀ =
     **…** s for the thermocline / atmosphere".
  5. *The growth map* (or the hint) — "the tongue closes at J = ¼ for R = 1".
  6. *Reading the current setting* — "Ri_min well below ¼: shear wins, billows (KH in a continuous profile) — the start of mixing in
     the thermocline and of billow clouds." · "Just below ¼: allowed, but only a narrow band of k grows, slowly." · "Above ¼: the
     theorem forbids growth whatever k — mixing schemes switch shear mixing off here." · R > 1: "A thin density layer inside a thicker
     shear layer: Ri_min sits off-centre."
- **Derivation tab:** **D17** (9 steps) `view: 'prof'`. **D18** (★★★, 14 steps) `view: 'prof'`; goal `set` {J: 0.1}; step 12 `live` "Im
  part with your profile: c_i ∫(N² − ¼U′²)/∣U − c∣² ∣φ∣² …"; step 14 `set` {J: 0.26} `watch` "N² − ¼U′² > 0 everywhere: the left side is
  positive, the right negative". Interpret: `s => verdict`.
- **Code:**
  ```python
  p = ch11.richardson_profiles("tanh", J={{J}}, R={{R}})
  z = np.linspace(-4, 4, 801)
  print(ch11.miles_howard_stable(z, U=p["U"], N2=p["N2"]))   # Ri_min = {{Rimin}} at z = {{zmin}}
  print(ch11.tg_growth({{k}}, {{J}}))                          # k c_i = {{kci}}
  ```
- **Walkthrough (6 steps):** 1. "A stratified shear layer" — "A current over denser water. Will it roll up into billows?" `play: true`
  · 2. "The Richardson number" — "Ri = N²/U′² (11.66): buoyancy's restoring against shear's driving; smallest at the centre, = J."
  `highlight: ['view:prof']` · 3. "The theorem" — "If Ri > ¼ everywhere, c_i = 0 for every k: no growth (Miles–Howard)." `derive:
  {id: 'D18', step: 13}` · 4. "Cross ¼" — "J = 0.24: allowed (a sliver of growth); J = 0.26: guaranteed stable." `set` {J: 0.26} ·
  5. "Necessary, not sufficient" — "Below ¼ only some k grow: drag k at J = 0.1." `controls: ['k']` · 6. "Your turn" — "Make the
  density layer thinner (R = 2) and predict where Ri_min moves; then check the profile." `controls: ['R', 'J']`.
- **Equations:** `N2` ref 'Eq. (11.56)' with (7.128) · `tg` ref 'Eq. (11.61)' · `phi` ref 'Eq. (11.63)' · `sa` ref 'Eq. (11.64)' ·
  `id` ref 'Eq. (11.65)' · `Ri` ref 'Eq. (11.66)' live · `mh` ref 'Eq. (11.67)' · symbols.
- **Check yourself:** (1) "Is Ri_min < ¼ enough for growth?" — "No: necessary only. At J = 0.1 waves with k ≳ 0.7 do not grow." ·
  (2) "Where is Ri smallest for R = 1, and why?" — "At z = 0: Ri = J cosh²z grows away from the centre." · (3) "Why are the computed
  eigenvalues symmetric about the real axis?" — "(11.61) has real coefficients: c and c* both solve it (N73)." · (4) "Thermocline
  preset: how long until a billow grows by e?" — "1/(kc_i) in units L/U₀ — the readout gives seconds."
- **Selftest parity rows:** `{name: 'Ri at z=1 J=0.1', js: Ri(1, 0.1, 1), py: 'ch11.gradient_richardson(1.0,
  N2=ch11.richardson_profiles("tanh", J=0.1)["N2"], dUdz=ch11.richardson_profiles("tanh", J=0.1)["Up"])', rtol: 1e-12}` · `{name: 'Ri_min
  R=2', js: mh(0.1, 2).Ri_min, py: 'ch11.miles_howard_stable(np.linspace(-4, 4, 8001), N2=ch11.richardson_profiles("tanh",
  J=0.1, R=2.0)["N2"], dUdz=ch11.richardson_profiles("tanh", J=0.1, R=2.0)["Up"])["Ri_min"]', rtol: 1e-4}` · `{name: 'map kci k=0.45 J=0', js: tgTable(0.45, 0), py: 'ch11.tg_growth(0.45, 0.0)',
  rtol: 2e-3}` · `{name: 'map kci k=0.4 J=0.1', js: tgTable(0.4, 0.1), py: 'ch11.tg_growth(0.4, 0.1)', rtol: 5e-3}`.
- **Fit plan:** 360×640: status, `billow` (50 %) over `prof` (50 %), `map` hidden (kc_i in readouts); 844×390: `billow` | `prof`;
  desktop rows [1.1, 1].

### E7 · inviscid_shear_criteria
- **Title:** "Which profiles can be unstable without viscosity?" · **Summary:** "Switch profiles: the inflection point, Rayleigh's and
  Fjørtoft's verdicts and the eigenvalues on the c-plane move together, always inside Howard's semicircle; open a neutral mode's cat's
  eye." · **CORE:** C12, C10 (also N91, N93, N94, N95, N96, N97, N82, N126, N127) · **Reference:** `fid_formula_lab.html` (verdict
  badges + term colours) + `amplitude_phase_second_order_II_3.html` (numbered live explanation).
- **meta:** `viz:order 7` · `viz:sections 11.7 11.9` · `viz:equations 11.71 11.72 11.81 11.83 11.84 11.85 11.86 11.87` · `viz:fluidpy
  ch11.inflection_points ch11.rayleigh_criterion ch11.fjortoft_criterion ch11.howard_semicircle ch11.in_howard_semicircle
  ch11.cats_eye_streamfunction ch11.cats_eye_width ch11.rayleigh_spectrum_table ch11.inviscid_profile` · `viz:derivations D19 D22`.
- **Physics:** `profile(name, b)` (U, U″ closed forms) ↔ `ch11.inviscid_profile`; `inflections(name)` ↔ `ST.inflection_points`;
  `rayleighOK(name)`, `fjortoft(name)` ↔ `ch11.rayleigh_criterion`, `ch11.fjortoft_criterion`; `inSemicircle(c, Umin, Umax)` ↔
  `ST.in_howard_semicircle`; `catsEye(x, y, A, phic, k, Uyc)` ↔ `ch11.cats_eye_streamfunction`; leading c(k) per profile from the table
  `RAYLEIGH` (C.5 5.7).
- **Views** (rows [1.15, 1]): 1. `prof` "U(y), U″ and the Fjørtoft product" (row 0, flex 1): U (orange), U″ (muted), inflection ●
  (purple) with U_I, the product (U − U_I)U″ shaded rose where negative, teal where positive; badges "Rayleigh ✓/✗", "Fjørtoft ✓/✗";
  click → inspector. 2. `cplane` "Wave speeds and Howard's semicircle" (row 0, flex 1): the semicircle on [U_min, U_max] (purple),
  the leading eigenvalue for the current k (rose) and its track over k (faint), the velocity range on the axis (orange). 3. `eye`
  "Cat's eye at the critical layer" (row 1, `hidePortrait`): streamlines of (11.87)'s expansion in the frame moving with c, the
  separatrix bold, width marked.
- **Controls (≤ 5):** `profile` chips Blasius-like / Poiseuille / wall-vorticity-max / tanh shear layer / Bickley jet / sin y · `k`
  0.05 … 2, default 0.445 · `b` (sin y half-width, optional) 1.2 … 3.0, default 1.4 · `A` (cat's-eye amplitude, optional) 0.01 … 0.5,
  default 0.1 · `parity` chips sinuous / varicose (Bickley only).
- **Presets:** "tanh: both pass" {profile: tanh, k: 0.445} · "Poiseuille: no inflection" {profile: poiseuille} · "wall-vorticity-max:
  Rayleigh ✓, Fjørtoft ✗" · "Bickley sinuous" {profile: bickley, parity: sinuous, k: 1} · "sin y, 2b < π" {profile: sin, b: 1.4}.
- **Status:** "🌀 inflection at y = 0 with a vorticity maximum: instability possible; c = 0 + 0.4267i (inside the semicircle)" (tanh,
  k = 0.445) · "🛡️ no inflection point: inviscidly stable (Rayleigh)" · "⚠️ inflection, but the vorticity peaks at the wall:
  Fjørtoft rules instability out" · "🧩 both criteria pass, yet stable (2b < π): necessary ≠ sufficient".
- **Readouts:** "y_I" · "U_I" · "c" · "distance to the arc" · "eye width".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **presets** · **status** · **inspector** (click U(y): "U″ =
  …, (U − U_I)U″ = … × … = …") · **modes** (profile family) · **terms** (the two integrals of (11.85) + (11.86) as bars: (U − c_r)
  part and (c_r − U_I) part summing to a negative total).
- **Explain:**
  0. *What the views show* — "**Profile**: U orange, U″ grey, the inflection point purple; rose shading = where (U − U_I)U″ < 0 (what
     Fjørtoft needs). **Wave speeds**: the computed leading eigenvalue (rose) inside Howard's half-disc (purple). **Cat's eye** (hidden
     on phones)."
  1. *Inflection points* — "U″ changes sign at y = **…** ⇒ Rayleigh's necessary condition (11.84) is **met / not met**".
  2. *Fjørtoft* — "U_I = **…**; min of (U − U_I)U″ = **…** ⇒ (11.86) **met / not met**".
  3. *Howard* — "centre ½(U_max + U_min) = **…**, radius ½(U_max − U_min) = **…**; your c = **…**: (c_r − m)² + c_i² = **…** ≤ R² ✓;
     growth bound k(U_max − U_min)/2 = **…** vs actual kc_i = **…**".
  4. *Cat's eye* (or the hint) — "near y_c: ψ ≈ ½(y − y_c)²U′(y_c) + Aφ(y_c)cos kx; width 2√(Aφ_c/U′_c) = **…**".
  5. *Reading the current setting* — per profile (tanh: "the classic mixing layer: unstable for 0 < k < 1, fastest k = 0.445");
     (Poiseuille: "stable without viscosity — its instability in C13 is viscous"); (wall-vorticity-max: "an inflection point is not
     enough: the vorticity must peak inside"); (Bickley: "jets are very unstable; sinuous mode first"); (sin y: "a counterexample:
     both conditions hold and still nothing grows for 2b < π").
- **Derivation tab:** **D22** (9 steps) `view: 'prof'`; step 7 `live` "Im{U″/(U − c)} = c_iU″/∣U − c∣²"; step 9 `watch` "the rose and
  teal shading must both appear". **D19** (★★★, 14 steps) `view: 'cplane'`; step 5 quotes slip #11; step 11 quotes slip #8; step 14
  `live` "your semicircle: centre …, radius …". Interpret.
- **Code:**
  ```python
  p = ch11.inviscid_profile("{{name}}")
  y = np.linspace(*p["domain"], 801)
  print(ch11.rayleigh_criterion(y, p["U"], p["Upp"]))   # {{ray}}
  print(ch11.fjortoft_criterion(y, p["U"], p["Upp"]))   # {{fj}}
  print(ST.in_howard_semicircle({{c}}, {{Umin}}, {{Umax}}))   # {{inside}}
  ```
- **Walkthrough (6 steps):** 1. "Shape decides?" — "Jets break up fast, channel flows slowly. Can the profile's shape tell?" ·
  2. "Rayleigh: an inflection point" — "Multiply Rayleigh's equation (11.81) by φ*, integrate, take the imaginary part: c_i∫U″/∣U −
  c∣²∣φ∣² = 0 (11.84) — U″ must change sign." `derive: {id: 'D22', step: 8}` · 3. "Fjørtoft: a maximum" — "Wall-vorticity-max has an
  inflection point, yet (U − U_I)U″ > 0 everywhere: stable." `set` wall preset, `terms: true` · 4. "Howard's half-disc" — "Every
  unstable c lies under the arc on [U_min, U_max]: tanh's c = 0.4267i at k = 0.445." `set` tanh, `derive: {id: 'D19', step: 14}` ·
  5. "Necessary only" — "sin y between walls: both criteria pass, yet stable for 2b < π." `set` sin preset · 6. "Your turn" —
  "Open the cat's eye and predict how its width changes when A is doubled (×√2)." `controls: ['A']`.
- **Equations:** `ray` ref 'Eq. (11.81)' · `id` ref 'Eq. (11.83)' · `im` ref 'Eq. (11.84)' · `fj` ref 'Eqs. (11.85), (11.86)' · `range`
  ref 'Eq. (11.71)' $U_{\min}<c_r<U_{\max}$ · `ineq` ref 'Eq. (11.72)' $\int[U^2-c_r^2-c_i^2]Q\,dz>0$ · `semi` "Howard's semicircle"
  ref 'p. 507 (unnumbered)' · `eye` ref 'Eq. (11.87)' · symbols.
- **Check yourself:** (1) "Why does Poiseuille flow get no inviscid instability?" — "U″ = −2 never changes sign; (11.84) then forces
  c_i = 0." · (2) "Could an unstable wave on the tanh layer travel at c_r = 1.5?" — "No: c_r must lie in (−1, 1) (11.71)." · (3) "What
  limits the growth rate of the tanh layer at k = 0.5?" — "kc_i ≤ (k/2)(U_max − U_min) = 0.5; actual ≈ 0.19." · (4) "Why does the
  cat's eye need a critical layer?" — "Only where U = c does the fluid move with the wave; streamlines close around that level."
- **Selftest parity rows:** `{name: 'tanh inflection', js: inflections('shear_layer')[0], py: 'ch11.rayleigh_criterion(np.linspace(-5,
  5, 1001), np.tanh, ch11.inviscid_profile("shear_layer")["Upp"])["y_I"][0]', rtol: 0, atol: 1e-6}` · `{name: 'semicircle point', js:
  inSemicircle([0.2, 0.9], -1, 1) ? 1 : 0, py: '1.0 if ch11.in_howard_semicircle(0.2+0.9j, -1.0, 1.0) else 0.0', rtol: 0, atol: 0}`
  (a conditional expression is syntax, not a builtin) · `{name: 'cats eye width', js: eyeWidth(0.1, 1, 1), py:
  'ch11.cats_eye_width(0.1, 1.0, 1.0)', rtol: 1e-14}` · `{name: 'tanh c_i table', js: rayTable('shear_layer', 0.45).ci, py:
  'np.imag(ch11.rayleigh_eigs(0.45, np.tanh, ch11.inviscid_profile("shear_layer")["Upp"], bc="decay", N=120)[0])', rtol: 2e-3}`.
- **Fit plan:** 360×640: status, `prof` over `cplane`, `eye` hidden; 844×390: `prof` | `cplane`; desktop rows [1.15, 1].

### E8 · orr_sommerfeld_neutral_curve
- **Title:** "How can viscosity make a channel flow unstable?" · **Summary:** "Click a point (Re, k) on the neutral-curve map: the
  Tollmien–Schlichting wave travels over the profile, and its production and dissipation bars show why it grows inside the thumb
  and decays outside." · **CORE:** C13, C11, C14 (also N88, N99, R23, N100, N102, N104, N105, N107, N108, R24) · **Reference:**
  `forced_damped_vibrations.html` (live Explain) + `fid_formula_lab.html` (term bars that add up to a total).
- **meta:** `viz:order 8` · `viz:sections 11.8 11.10 11.11` · `viz:equations 11.77 11.78 11.79 11.80 11.88` · `viz:fluidpy
  ch11.orr_sommerfeld_eigs ch11.poiseuille_neutral_curve ch11.disturbance_energy_budget ch11.squire_transform ch11.table_11_1
  ch11.os_mode ch11.ts_wave_fields` · `viz:derivations D21 D23`.
- **Physics:** tables `OS` (C.5 5.8): per flow (Poiseuille, Blasius, tanh, Bickley, Couette) the neutral curve (Re, k_lower,
  k_upper), and on a 24 Re × 25 k grid the leading c and the budget (E, P, Λ); mode samples (φ, û, v̂ on 41 y) at the presets ↔
  `ch11.orr_sommerfeld_eigs`, `ch11.poiseuille_neutral_curve`, `ch11.disturbance_energy_budget`, `ch11.os_mode`; `squire(k, m, Re)`
  live ↔ `ST.squire_transform`; `tsFields(x, y, t, k, c, mode)` ↔ `ch11.ts_wave_fields`. Bilinear interpolation in (log Re, k) between
  grid nodes, labelled "table"; presets sit on grid nodes or sampled modes.
- **Views** (rows [1.1, 1]): 1. `map` "Neutral curve (Re, k)" (row 0, flex 1): log Re (10 … 10⁵; Blasius Re_δ*), k; the flow's neutral
  curve (purple), unstable region rose, Re_c ◆, the current point; for an oblique wave a ghost point at (Re̅, k̄) (Squire) with an
  arrow; click → set (Re, k). 2. `wave` "The wave over the profile" (row 0, flex 1.3): the base profile U(y) (orange, muted fill),
  perturbation streamlines (ψ contours, rose/teal by sign) travelling at c_r, the critical levels y_c (U = c_r, ×) and the wall layers
  marked; amplitude e^{kc_i t}. 3. `budget` "Energy budget" (row 1, `hidePortrait`): bars production P (orange), dissipation −Λ
  (blue), dE/dt = 2kc_iE (rose if > 0, teal if < 0), each ÷ E; a small inset: the Reynolds stress −⟨uv⟩(y) profile.
- **Controls (≤ 5):** `flow` chips Poiseuille / Blasius / tanh / Bickley / Couette (mode) · `lRe` "$\log_{10}\mathrm{Re}$" 1 … 5,
  default 4 · `k` 0.1 … 2, default 1 · `theta` "oblique angle" 0 … 60° (optional, Squire) · transport.
- **Transport:** `t` over 3 wave periods 2π/(kc_r), `end: 'loop'`.
- **Presets:** "Poiseuille Re_c" {flow: poiseuille, lRe: log10(5772.22), k: 1.02056} · "Orszag: Re = 10⁴, k = 1" {lRe: 4, k: 1} ·
  "Re = 5000 (decays)" {lRe: log10(5000), k: 1} · "Blasius Re_δ* = 519" {flow: blasius, lRe: log10(519.2), k: 0.303} · "Couette:
  stable" {flow: couette, lRe: 4, k: 1} · "tanh: Re_c = 0" {flow: tanh, lRe: 1, k: 0.5}.
- **Status:** "📈 inside the neutral curve: P/Λ = 1.62 — the TS wave grows (kc_i = 0.00374)" (Orszag) · "📉 outside: P/Λ = 0.77 — it
  decays" (Re = 5000) · "🛡️ plane Couette: every mode decays at every Re (N100)" · "🌀 inflectional profile: unstable down to tiny Re".
- **Readouts:** "c (table)" · "kc_i" · "P/Λ" · "Re̅ (Squire)" · "Re_c (this flow)".
- **Terms:** P, −Λ, dE/dt (sum check: P − Λ = dE/dt to the table's precision; the residual printed).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **terms** · **inspector** (click the map: "Re =
  …, k = …: c = …, E = …, P = …, Λ = …, P − Λ = … = 2kc_iE ✓") · **presets** · **status** · **modes** (flows).
- **Explain:**
  0. *What the views show* — "**Map**: the thumb is where waves grow (rose). **Wave**: the disturbance streamlines over the base
     profile; × marks where U = c_r. **Budget** (hidden on phones): what feeds the wave (orange, Reynolds-stress production) and what
     drains it (blue, viscous dissipation); their difference is the growth."
  1. *Your mode (table)* — "c = **…** at Re = **…**, k = **…** (Chebyshev, N = 100, ours); kc_i = **…** per L/U₀".
  2. *Is it inside?* — "neutral curve at this Re: k from **…** to **…** ⇒ **inside / outside**; Re_c of this flow = **…**".
  3. *The energy budget (11.88)* — "dE/dt = −∫uvU′ dV − Λ: P = **…**, Λ = **…**, P − Λ = **…** = 2kc_iE = **…**" boxed.
  4. *Squire (11.78)* — "oblique angle θ = **…**: k̄ = k/cos θ … Re̅ = kRe/k̄ = Re cos θ = **…**: the 3-D wave behaves like the 2-D wave
     at the ghost point".
  5. *Table 11.1* — the six-row table with the current flow lit (our Re_c beside the published benchmark).
  6. *Reading the current setting* — inside: "Viscosity near the wall shifts v against u (the phase in the inset departs from 90°):
     −⟨uv⟩U′ > 0 beats dissipation — a Tollmien–Schlichting wave." · outside: "The tilt is too weak (low Re) or the wave too short
     (dissipation wins)." · Couette: "No thumb: linear theory never finds an unstable wave; transition needs finite disturbances." ·
     tanh/Bickley: "Inflectional: viscosity only trims an inviscid instability (C12)." · Blasius: "Parallel-flow approximation: the
     real layer grows downstream (N109)."
- **Derivation tab:** **D21** (10 steps) `view: 'wave'`; step 1 `watch` "û = φ′, v̂ = −ikφ"; step 8 `live` "with your k and Re the
  viscous factor 1/(ikRe) = …". **D23** (12 steps) `view: 'budget'` on desktop, `'map'` on phones; step 11 `live` "P = …, Λ = …";
  step 12 `watch` "the bars add up to dE/dt". Interpret.
- **Code:**
  ```python
  U = lambda y: 1 - y**2; Upp = lambda y: -2 + 0*y           # plane Poiseuille, half-width units
  c = ch11.orr_sommerfeld_eigs({{k}}, {{Re}}, U, Upp, N=100)[0]   # {{c}}
  m = ch11.os_mode({{k}}, {{Re}}, U, Upp)
  b = ch11.disturbance_energy_budget({{k}}, m["c"], m["phi"], m["y"], -2*m["y"], {{Re}})
  print(b["production"]/b["dissipation"])                    # P/Λ = {{ratio}}
  print(ST.squire_transform({{k}}, {{m}}, {{Re}}))           # Re̅ = {{Rebar}}
  ```
- **Walkthrough (7 steps):** 1. "A channel flow" — "No inflection point, so no inviscid instability. Yet…" `play: true` · 2. "The
  thumb" — "The neutral curve c_i = 0: rose inside, waves grow. It starts at Re_c = 5772.22, k_c = 1.02 (Orszag 1971)." `set` Re_c
  preset · 3. "Orszag's mode" — "Re = 10⁴, k = 1: c = 0.2375 + 0.0037i — a slow, weakly growing TS wave." `set` Orszag,
  `readouts: ['c', 'kci']` · 4. "Who feeds it?" — "Production −∫uvU′ beats dissipation: P/Λ = 1.62." `terms: true`, `derive: {id: 'D23',
  step: 12}` · 5. "Viscosity tilts the wave" — "The u–v phase departs from 90° near the walls: that tilt is the energy pump."
  `highlight: ['view:wave']` · 6. "Squire" — "Tilt the wave 30°: it sees Re cos 30° — less unstable. 2-D waves go first." `set`
  {theta: 30} · 7. "Your turn" — "Predict whether Re = 8000, k = 1.2 grows; check the map (edge ≈ inside?)." `controls: ['lRe', 'k']`.
- **Equations:** `3d` ref 'Eq. (11.77)' · `squire` ref 'Eq. (11.78)' live · `os` ref 'Eq. (11.79)' · `bc` ref 'Eq. (11.80)' · `energy`
  ref 'Eq. (11.88)' and the 2-D form (p. 520) live · symbols.
- **Check yourself:** (1) "Why can't plane Couette flow have a thumb?" — "Our grid finds c_i < 0 everywhere; linear theory has no
  growing mode (N100)." `set Couette` · (2) "At Re = 10⁴, which k grow?" — "About 0.80 to 1.07 (the thumb's width there)." · (3) "What
  does a 30° oblique wave at Re = 6000 do?" — "It behaves like a 2-D wave at Re̅ = 5196 < 5772: it decays." `set {lRe: 3.778, theta:
  30}` · (4) "Why is Blasius's Re_c (519) not comparable with Poiseuille's (5772)?" — "Different lengths and speeds: δ* and U∞ vs the
  half-width and centreline speed (P200)."
- **Selftest parity rows:** `{name: 'Orszag c_r (table)', js: osTable('poiseuille', 1e4, 1).cr, py:
  'np.real(ch11.orr_sommerfeld_eigs(1.0, 1e4, lambda y: 1 - y**2, lambda y: -2 + 0*y, N=100)[0])', rtol: 1e-6}` · `{name: 'Orszag c_i
  (table)', js: osTable('poiseuille', 1e4, 1).ci, py: 'np.imag(ch11.orr_sommerfeld_eigs(1.0, 1e4, lambda y: 1 - y**2, lambda y: -2 +
  0*y, N=100)[0])', rtol: 1e-4}` — NB: `lambda` is a keyword, allowed in eval · `{name: 'Squire Rebar', js: squire(1, Math.tan(Math.PI/6),
  1e4).Rebar, py: 'ch11.squire_transform(1.0, np.tan(np.pi/6), 1e4)["Rebar"]', rtol: 1e-12}` · `{name: 'Re_c', js: OS.poiseuille.Re_c,
  py: 'ch11.poiseuille_critical()["Re_c"]', rtol: 1e-5}` · `{name: 'P/Lambda Orszag', js: osTable('poiseuille', 1e4, 1).ratio, py:
  'ch11.disturbance_energy_budget(1.0, ch11.os_mode(1.0, 1e4, lambda y: 1 - y**2, lambda y: -2 + 0*y)["c"], ch11.os_mode(1.0, 1e4,
  lambda y: 1 - y**2, lambda y: -2 + 0*y)["phi"], ch11.os_mode(1.0, 1e4, lambda y: 1 - y**2, lambda y: -2 + 0*y)["y"], -2*ch11.os_mode(
  1.0, 1e4, lambda y: 1 - y**2, lambda y: -2 + 0*y)["y"], 1e4)["ratio"]', rtol: 1e-4}`.
- **Fit plan:** 360×640: status, `map` (45 %) over `wave` (55 %), `budget` hidden (P/Λ in the status and readouts); 844×390: `map` |
  `wave`; desktop rows [1.1, 1] with `budget` under `wave`.

### E9 · lorenz_attractor
- **Title:** "How can three exact equations be unpredictable?" · **Summary:** "Two starts 10⁻⁸ apart on the Lorenz attractor, on one
  clock: the 3-D orbit, X(t) of both runs and the log of their separation; move r through 1 and 24.74 and watch the verdict change."
  · **CORE:** C15 (also N112–N117, N119; 3-D) · **Reference:** `ddpm3_unet_3d.html` (orbitable 3-D scene, inspector) +
  `angular_frequency_explorer_1.html` (linked views, end-of-run summary).
- **meta:** `viz:order 9` · `viz:sections 11.14` · `viz:equations 11.90 11.91` · `viz:fluidpy ch11.lorenz_rhs ch11.lorenz_rk4
  ch11.lorenz_fixed_points ch11.lorenz_hopf_r ch11.lorenz_b ch11.lorenz_eigs ch11.lorenz_fields` · `viz:derivations D24 D25`.
- **Physics:** `rhs(s, Pr, r, b)` + `Viz.num.rk4Step` with dt = 0.005 (C.5 5.9) ↔ `ch11.lorenz_rhs`, `ch11.lorenz_rk4`; `fixed(r, b)` ↔
  `ch11.lorenz_fixed_points`; `rH(Pr, b)` ↔ `ch11.lorenz_hopf_r`; `bOf(k)` ↔ `ch11.lorenz_b`; `eigC(r, Pr, b)` (closed-form cubic
  roots) ↔ `ch11.lorenz_eigs`; `fields(x, z, X, Y, Z, k)` ↔ `ch11.lorenz_fields`.
- **Views** (rows [1.2, 1]): 1. `scene` "The attractor (3-D)" (row 0, flex 1.4; `Viz.three`, orbit/zoom, click picks a point): run 1
  orange, run 2 blue (trails of the last 5 time units bold, the rest faint), C± (purple spheres) and the origin (muted); offline
  fallback: X–Z projection on the 2-D canvas. 2. `traces` "X(t) and log₁₀∣δ∣" (row 0, flex 1): top half X(t) for both runs, bottom half
  log₁₀ of their separation with the fitted slope (rose dashed) and the 'differ by 1' line. 3. `roll` "The convection roll" (row 1,
  `hidePortrait`): ψ and T′ of (11.90) driven by run 1's X, Y, Z (the roll reverses when X changes sign).
- **Controls (≤ 5):** `r` 0.5 … 35, step 0.1, default 28 · `Pr` 1 … 20 (optional), default 10 · `b` 1 … 4 (optional; chip "b from
  k_c: 8/3") · `ld0` "$\log_{10}\delta_0$" −12 … −2, default −8 · transport.
- **Transport:** `t` 0 … 40, rate 2 time units/s, `end: 'hold'`; end card: "the runs differ by 1 after t = 29.1 (δ₀ = 10⁻⁸); fitted
  growth rate ≈ … per unit time".
- **Presets:** "r = 0.5 conduction" · "r = 10 steady convection" · "r = 20 transient chaos" · "r = 28 Lorenz" · "r = 24: a long chaotic
  transient below r_H".
- **Status:** "🌪️ r = 28 > r_H = 24.74: both convection states unstable — chaos" · "🔁 1 < r < r_H: steady convection (C±) after a
  transient" (r = 20: "long chaotic transient, then C±") · "😴 r < 1: conduction, the origin attracts".
- **Readouts:** "r_H" · "C± = (±√(b(r − 1)), …, r − 1)" · "eig(C±)" · "δ(t)" · "volume rate −(Pr + 1 + b)".
- **Depth features:** Explain + Code + Derivation · **3-D** · **linked views** (3) · **transport** · **presets** · **status** ·
  **inspector** (click the orbit: the state and the Jacobian's eigenvalues there) · end-of-run card.
- **Explain:**
  0. *What the views show* — "**Attractor**: two runs, orange and blue, starting 10⁻⁸ apart; purple balls = the steady convection
     states. **Traces**: X(t) for both (on top of each other at first) and the log of their distance. **Roll** (hidden on phones): the
     convection cell that X, Y, Z describe — X > 0 turns it one way, X < 0 the other."
  1. *The model* — "(11.91) with Pr = **10**, r = **28**, b = **8/3** (= 4π²/(π² + k²) at k² = π²/2, the free–free K_c of C05)".
  2. *Steady states* — "origin; C± = (±√(b(r − 1)), ±√(b(r − 1)), r − 1) = (±**8.485**, ±**8.485**, **27**)" boxed.
  3. *Their stability* — "origin: eigenvalues **11.83**, **−2.667**, **−22.83** (unstable for r > 1); C±: the cubic λ³ + (Pr + b +
     1)λ² + b(r + Pr)λ + 2bPr(r − 1) = 0 has roots **−13.85**, **0.094 ± 10.19i** — Hopf point r_H = Pr(Pr + b + 3)/(Pr − b − 1) =
     **24.74**" boxed.
  4. *Volume* — "∇·ṡ = −(Pr + 1 + b) = **−13.67**: volumes shrink by e every 0.073 time units — the attractor has zero volume."
  5. *Sensitivity* — "δ₀ = 10⁻⁸ → ∣δ∣ = 1 at t ≈ **29.1** (this run); fitted slope ≈ **…** (literature λ₁ ≈ 0.906); each factor 10 of
     extra precision buys only ln 10/λ ≈ **2.5** time units".
  6. *At the current time* — "t = `Viz.live('t')`, X₁ = …, X₂ = …, ∣δ∣ = …".
  7. *Reading the current setting* — r < 1: "Heating below Ra_c: conduction." · 1 < r < 13.93 (the homoclinic value, Sparrow 1982 — cited, not computed here): "Steady
     rolls." · 13.93 < r < 24.74: "Chaotic transients, then steady rolls — preturbulence." · r > 24.74: "No stable steady state: the
     roll reverses irregularly forever. Predictability is limited, as in weather (ensembles, Ch. 13)."
- **Derivation tab:** **D24** (★★★, 15 steps) `view: 'roll'` on desktop, `'traces'` on phones; step 6 `watch` "the three modes of
  (11.90)"; step 14 `live` "r = Ra/Ra_c(k) with your Ra"; step 15 `live` "b = 4π²/(π² + k²) = 8/3 at k² = π²/2". **D25** (11 steps)
  `view: 'scene'`; step 6 `live` "λ² + (Pr + 1)λ + Pr(1 − r) = 0 at your r"; step 11 `live` "r_H = 24.74 with your Pr, b".
- **Code:**
  ```python
  s = ch11.lorenz_integrate((1.0, 1.0, 1.0), 40.0, Pr={{Pr}}, r={{r}}, b={{b}})
  print(ch11.lorenz_fixed_points({{r}}, {{b}}))       # C± = {{fp}}
  print(ch11.lorenz_eigs({{r}}, {{Pr}}, {{b}}))       # {{eig}}
  print(ch11.lorenz_hopf_r({{Pr}}, {{b}}))            # r_H = {{rH}}
  print(ch11.lorenz_predictability_time({{d0}}))      # t at |δ| = 1: {{tp}}
  ```
- **Walkthrough (7 steps):** 1. "Convection in three numbers" — "Lorenz kept one roll and two temperature modes. Press ▶." `play:
  true` · 2. "Where the equations come from" — "Insert (11.90) into the Boussinesq equations, project, rescale: (11.91)." `derive: {id:
  'D24', step: 15}` · 3. "Steady rolls" — "r = 10: both runs spiral into a convection state C±." `set` {r: 10} · 4. "The Hopf point" —
  "Past r_H = 24.74 the spirals turn outward: no steady state is stable." `derive: {id: 'D25', step: 11}`, `set` {r: 28} · 5. "Two runs
  part" — "δ₀ = 10⁻⁸: the traces agree until t ≈ 25, then are unrelated. ln∣δ∣ grows linearly." `highlight: ['view:traces']` · 6.
  "Volume shrinks, motion never settles" — "Volumes contract at 13.67 per unit time: a zero-volume, never-repeating set." · 7. "Your
  turn" — "Predict how much later the runs part with δ₀ = 10⁻¹². Then set it." `controls: ['ld0']`.
- **Equations:** `trunc` ref 'Eq. (11.90)' · `lorenz` ref 'Eq. (11.91)' live · `fp` "Steady states" ref 'p. 528 (unnumbered)' · `rH`
  "Hopf point (ours, D25)" ref 'no book number' · `div` "Volume contraction (ours)" · symbols X, Y, Z, r, b, Pr (–).
- **Check yourself:** (1) "Why does halving δ₀ hardly help?" — "The error grows like e^{λt}: halving δ₀ buys ln 2/λ ≈ 0.8 time units." ·
  (2) "At r = 20, is the motion chaotic forever?" — "No: a long chaotic transient, then C± (r < r_H)." `set {r: 20}` · (3) "What does
  X < 0 mean for the roll?" — "The roll turns the other way (ψ ∝ X)." · (4) "Why must the attractor have zero volume?" — "Every volume
  shrinks at the constant rate Pr + 1 + b."
- **Selftest parity rows:** `{name: 'rk4 200 steps X', js: rk4Run([1, 1, 1], 0.005, 200, 10, 28, 8/3)[0], py:
  'ch11.lorenz_rk4((1.0, 1.0, 1.0), 0.005, 200)[200][0]', rtol: 1e-10}` · `{name: 'r_H', js: rH(10, 8/3), py: 'ch11.lorenz_hopf_r(10.0,
  8/3)', rtol: 1e-14}` · `{name: 'C+ x', js: fixed(28, 8/3)[1][0], py: 'ch11.lorenz_fixed_points(28.0, 8/3)[1][0]', rtol: 1e-14}` ·
  `{name: 'b(kc)', js: bOf(Math.PI/Math.SQRT2), py: 'ch11.lorenz_b(np.pi/np.sqrt(2))', rtol: 1e-14}` · `{name: 'eig C real', js:
  eigC(28, 10, 8/3)[1][0], py: 'np.max(np.real(ch11.lorenz_eigs(28.0, 10.0, 8/3)))', rtol: 1e-8}`.
- **Fit plan:** 360×640: status, `scene` (55 %) over `traces` (45 %), `roll` hidden; 844×390: `scene` | `traces`; desktop rows [1.2, 1].

### B1 · period_doubling_route
- **Title:** "What does a route to chaos look like?" (backup; built only if an explainer above fails review) · **Summary:** "Cobweb
  the logistic map and watch its fixed point split into 2, 4, 8 … cycles, each split sooner by Feigenbaum's 4.669." · **CORE:** C15 (also
  N117, N129, N114) · **Reference:** `pixels_as_parameters.html` (an algorithm stepped in slow motion) + `overfitting_curves.html`.
- **meta:** `viz:order 10` · `viz:sections 11.14` · `viz:equations` none numbered (the Feigenbaum ratio is unnumbered, p. 530) ·
  `viz:fluidpy ch11.logistic_map ch11.cobweb ch11.bifurcation_diagram ch11.period_doubling_points ch11.feigenbaum_estimate` ·
  `viz:derivations none`.
- **Views:** `cobweb` (the map y = Ax(1 − x), the diagonal, the cobweb path; fixed point and ∣f′∣) · `tree` (bifurcation diagram 2.8 …
  4 with A₁, A₂, A₃ markers and the current A line) · `table` (hidePortrait: A_n, superstable S_n, δ estimates 4.709, 4.681, 4.663,
  4.668, 4.669). **Controls:** A 2.5 … 4 · x₀ · iterates · map chips logistic / sine. **Presets:** A = 2.8, 3.2, 3.5, 3.57, 3.83.
  **Status:** "fixed point x* = 1 − 1/A, ∣f′∣ = ∣2 − A∣ = 0.8: stable" · "period 2" · "period 4" · "chaos" · "period-3 window".
  **Depth:** linked views, presets, status, inspector (click the cobweb: x_n → x_{n+1} arithmetic). **Walkthrough** (5 steps), **Explain**
  (fixed point and its stability; the 2-cycle from f(f(x)) = x; the δ table; reading), **Check** (3), **parity rows**: `logistic
  x after 30 at A = 2.8` vs `ch11.logistic_map(2.8, 0.2, 30)[-1]`; `A₂` vs `ch11.period_doubling_points(3)["A_n"][1]`.

---

## Part D — runtime budget (full run < 5 min on a laptop / Colab CPU)

Chapter 11's cost sits in four places: the sympy engines (ten, each `lru_cache`d, 2–15 s), the eigen-solver sweeps (TG growth map,
Rayleigh spectra, Orr–Sommerfeld neutral curves and critical points — tens of seconds to minutes cold), the animations (5) and the
Lorenz integrations. Single eigen-solves are cheap (measured for this design: Bénard N = 40 ≈ 5 ms, Taylor N = 40 ≈ 30 ms, Rayleigh/TG
N = 100 ≈ 0.2–0.5 s, Orr–Sommerfeld N = 100 ≈ 0.02 s for eigenvalues only, a Poiseuille brentq at N = 80 ≈ 0.5 s in all).

| Section | Heaviest cells | Full (cold `outputs/` cache; `reference/ch11/` present) | Full (warm cache) | FAST (`FLUIDPY_FAST=1`) |
|---|---|---|---|---|
| setup + imports | numpy/scipy/sympy/plotly, `core.stability`, `ch11` | 12 s | 12 s | 12 s |
| §11.1 | A1 (4 shapes × 2 kicks, 60 frames) | 8 s | 8 s | 4 s (30 frames) |
| §11.2 C01 | three-mode sum, σ curve figure | 3 s | 3 s | 3 s |
| §11.3 C02 | `kh_sympy`, `kh_depth_tension_sympy`, residuals, 2 figures, F1 (24 × 2 traces), **A2** (`sheet_rollup` N = 200, 60 frames) | 30 s | 22 s | 14 s (N = 100, 30 frames, F1 12) |
| §11.4 C03 | `benard_perturbation_sympy`, `exchange_of_stabilities_sympy`, convergence N ≤ 40, spectrum figure | 18 s | 10 s | 10 s |
| §11.4 C04 | three `benard_critical`, determinant scans, `benard_neutral_table` (96 K × 4), eigenfunction figures, F2 (30), live, **A3** (60 frames) | 16 s | 12 s | 8 s |
| §11.4 C05 | `benard_free_free_sympy`, two figures | 4 s | 3 s | 3 s |
| §11.5 C06 | closed forms, cubic, regime map, F3 (20) | 4 s | 4 s | 3 s |
| §11.6 C07 | `taylor_perturbation_sympy`, `taylor_critical` ×4, `taylor_critical_table` (31 μ, cached), F4 (21 μ × 30 k), Galerkin, 2 figures | 45 s | 14 s | 10 s (F4 11 μ × 20 k) |
| §11.7 C08 | `stratified_shear_sympy`, TG N = 100/140 spectra | 22 s | 12 s | 10 s (N = 80/120) |
| §11.7 C09 | Miles–Howard scans, `tg_growth` ×3, identity check, growth map (from `reference/ch11/tg_growth_map.csv`), F5 (26) | 8 s | 6 s | 5 s |
| §11.7 C10 | 18 TG solves for the semicircle figure, identity check, F6 (from `rayleigh_spectra.json`) | 12 s | 8 s | 5 s (9 solves) |
| §11.8 C11 | OS N = 60…140, `os_3d_eigs` (N = 60), `os_derivation_sympy`, spectrum figure | 14 s | 10 s | 8 s |
| §11.9 C12 | six-profile verdicts, `tanh_max_growth`, `sin_profile_max_growth` (4 b × 30 k), Rayleigh/OS comparison, cat's eye | 22 s | 12 s | 8 s (3 b × 15 k) |
| §11.10 C13 | `poiseuille_critical` (cached; cold coarse N = 60 ≈ 20 s), neutral curves (from `reference/ch11/os_neutral_*.csv`), from-scratch brentq (N = 80), Couette grid (30), `table_11_1` (cached), F7 (from JSON) | 40 s | 10 s | 8 s (Couette 12) |
| §11.10 C14 | `energy_equation_sympy`, two budgets, fine-grid from-scratch, 3-panel figure, **A4** (60 frames) | 16 s | 12 s | 8 s |
| §11.11–11.13 | Tollmien figure, Blasius curves (cached) | 4 s | 3 s | 3 s |
| §11.14 C15 | `lorenz_sympy`, 4 DOP853 runs (t = 40), RK4 loop, spectra, pendulum/Hopf, logistic tree (1200 A × 600 iterations, vectorised), **A5** (90 frames), F8 (30 r, cached), live | 35 s | 22 s | 14 s (A5 40 frames, F8 15) |
| explainers (9 × `show_viz`) | read the HTML files | 2 s | 2 s | 2 s |
| **Total** | | **≈ 315 s** | **≈ 185 s** | **≈ 138 s** |

**FAST plan.** Every size-dependent choice is written `a if not FAST else b` in the cell: animations 60/90 → 30/40 frames (dpi 80),
`sheet_rollup` N = 200 → 100, TG/Rayleigh N = 100/140 → 80/120 (still N-converged for the leading modes), the sin-profile sweep 4 × 30
→ 3 × 15, the Couette grid 30 → 12 points, F1/F4/F8 slider counts halved (≤ 15 steps), Taylor F4 30 → 20 k per μ. **The cold full run
exceeds 300 s only because of three one-off computations** — the Taylor critical table (31 μ), the Poiseuille critical point and the
F8 Lorenz sweep; these are cached in `outputs/ch11/*.npz` (keyed by their parameters, P252) **and** shipped as small public tables in
`reference/ch11/` (`taylor_critical.csv`, `benchmarks.json` with our Poiseuille numbers, `os_neutral_*.csv`, `tg_growth_map.csv`,
`rayleigh_spectra.json`, `benard_neutral_curves.csv`), so a fresh Colab runtime reads them and stays ≈ 200 s; the notebook says
"read from reference/ch11 (computed by scripts/ch11_tables.py)" and keeps one live recomputation per table (one Taylor μ, one
Poiseuille brentq, one TG point) so a reader sees the method run. The Blasius and Falkner–Skan neutral curves and the optional Lyapunov
exponent are script-only (`scripts/ch11_neutral_curves.py`, `scripts/ch11_lorenz.py`). Sympy cells never `simplify` expressions with
generic functions beyond the final residual. **Outputs:** 4 videos (A2–A5) + 1 frame player (A1), each < 3 MB at dpi 80; 8 plotly
figures (≤ 30 steps × ≤ 4 traces × ≤ 400 points); ≈ 45 static figures — the page stays under 15 MB.

---

## Part E — prerequisite ledger
Every concept, symbol, maths tool and Python function or idiom the notebook or its explainers use, with where it is explained.
"primer (in Cxx)" = a 📎 primer placed in that block before first use (the Concept text is the primer term, used verbatim in
`nb.primer`); "knowledge/primers.md: <term> (chNN Pnn) — reminder" = a one-line reminder naming the earlier primer; a CORE/RECAP id
alone = taught there; "Cxx (Nnn)" = the NOTE placed in that block; "Cxx (D0n)" = the derivation where it is used; "Cxx (gloss …)" =
one sentence where it is used. Earlier chapters' material that is neither a primer nor a ch11 RECAP is named by section ("Ch. 4 §4.9").
New primers P255–P279 (25) in first-use order: necessary/sufficient (C01), np.roots (C02), ODE eigenvalue problem, Chebyshev, BC rows +
eig(A, B), real/imaginary parts (C03), parity, cube roots, neutral-curve recipe, Helmholtz planforms (C04), quotient rule (C05),
uniqueness (C06), narrow gap (C07), quadratic eigenproblem, singular point, infinite-domain map (C08), Clenshaw–Curtis (C09),
completing the square (C10), periodic x-derivative (C14), Jacobian and fixed points, Galerkin truncation, Hopf cubic condition,
Lyapunov exponent, matplotlib 3-D, iterated maps (C15).

| Concept | First used in | Explained by |
|---|---|---|
| basic state and infinitesimal disturbance; linear vs finite-amplitude instability | C01 | C01 (N01) |
| four potential wells (bowl, cap, plane, dimple) | C01 | C01 (N02, A1) |
| temporal vs spatial instability; no Coriolis here | C01 | C01 (N03) |
| normal mode (11.1), both forms | C01 | C01 |
| σ = −i∣K∣c, σ_r = ∣K∣c_i | C01 | C01 (D01) |
| stable, neutral, unstable (for every k) | C01 | C01 (N04) |
| superposition of modes; modes evolve independently | C01 | C01 (N05) |
| marginal state; stationary vs oscillatory onset (overstability) | C01 | C01 (N06) |
| necessary vs sufficient conditions, "for every k", and proof by contradiction | C01 | primer (in C01) |
| complex numbers and Euler's formula | C01 | knowledge/primers.md: square root of a negative number (ch01 P45) — reminder |
| complex amplitudes and Re{·} | C01 | knowledge/primers.md: complex amplitudes (ch07 P176) — reminder |
| complex conjugate | C01 | knowledge/primers.md: complex conjugate (ch02 P81) — reminder |
| the complex plane in numpy (1j, np.real, np.imag) | C01 | knowledge/primers.md: the complex plane in numpy (ch06 P153) — reminder |
| Fourier modes and amplification factors (Ch. 10 bridge) | C01 | knowledge/primers.md: Fourier modes and the FFT Poisson solver (ch05 P142) — reminder; Ch. 10 §10.2 (G(θ)) |
| linear stability by eigenvalues (perturb, linearise) | C01 | knowledge/primers.md: linear stability of a steady configuration (perturb, linearise, eigenvalues) (ch09 P214) — reminder |
| solve_ivp | C01 | knowledge/primers.md: scipy.integrate.solve_ivp (ch01 P31) — reminder |
| animate and show_animation (frames player) | C01 | knowledge/primers.md: animate and show_animation (ch01 P16) — reminder |
| show_viz | C01 | knowledge/primers.md: show_viz (ch01 P18) — reminder |
| np.all, np.any | C01 | C01 (gloss in the P255 demo) |
| ST.normal_mode, sigma_from_c, stability_class, stability_verdict, marginal_type | C01 | C01 (code cells, explain lists) |
| velocity potential and Kelvin's theorem, (11.2) | C02 | R01 |
| Laplace's equation for the disturbances (11.3) | C02 | R02 |
| decay far away (11.4)–(11.5) | C02 | R03 |
| kinematic and dynamic interface conditions (11.6)–(11.7) | C02 | R04 |
| level sets and their normals, (11.8) | C02 | R05; knowledge/primers.md: level sets and the directional derivative (ch02 P75) — reminder |
| unsteady Bernoulli (11.10) | C02 | R06 |
| A″ = k²A and exponential solutions | C02 | R07; knowledge/primers.md: linear second-order ODE (ch01 P44) — reminder |
| separation of variables | C02 | knowledge/primers.md: separation of variables for a PDE (ch07 P167) — reminder |
| decaying potentials (11.15) | C02 | R08 |
| static interface waves (11.19) = (7.95); Rayleigh–Taylor | C02 | R09 |
| vortex sheet (11.20); growth kΔU/2 at every k | C02 | R10 |
| nonlinear roll-up of a vortex sheet | C02 | R11 |
| Kelvin–Helmholtz set-up (two streams, interface) | C02 | C02 (N07) |
| exact kinematic condition, square-root factor cancels | C02 | C02 (N08, D02) |
| linearised kinematic condition (11.9) | C02 | C02 (N09, D02) |
| Taylor transfer of a boundary condition to z = 0 | C02 | knowledge/primers.md: Taylor transfer of a boundary condition (ch07 P166) — reminder |
| orders of smallness, dropping products | C02 | knowledge/primers.md: limits and orders of smallness (ch02 P68) — reminder |
| pressure matching (11.11), undisturbed balance (11.12), (11.13) | C02 | C02 (N10, N11, N12, D02) |
| normal-mode trial (11.14) | C02 | C02 (N13) |
| remnants (11.16)–(11.17), amplitudes A± | C02 | C02 (N14, D03) |
| the KH quadratic and its discriminant | C02 | C02 (N15, D03) |
| quadratic formula with complex roots | C02 | knowledge/primers.md: complex square roots and the quadratic formula (ch06 P159) — reminder |
| KH dispersion relation (11.18) | C02 | C02 (D03) |
| instability criterion g(ρ₂² − ρ₁²) < kρ₁ρ₂ΔU², k_c | C02 | C02 (N16, D03) |
| inequalities under a sign change | C02 | knowledge/primers.md: inequalities under a sign change (ch01 P48) — reminder |
| conjugate pairs of growing/decaying modes | C02 | C02 (N17) |
| short waves always unstable | C02 | C02 (N18) |
| moving frame, c_r = mean | C02 | C02 (N19); knowledge/primers.md: frames of reference and relative velocity (ch03 P96) — reminder |
| shear instability in nature (billows, thermocline) | C02 | C02 (N20) |
| energy source, mixing energy 2/3 | C02 | C02 (N21, N22) |
| Cauchy–Schwarz: smoothing lowers ∫U² | C02 | C02 (N23) |
| depth h and surface tension (Exercise 11.1), coth kh | C02 | C02 (N120, D04) |
| hyperbolic functions cosh, sinh, coth | C02 | knowledge/primers.md: hyperbolic functions cosh, sinh, tanh (ch07 P168) — reminder |
| Laplace pressure jump σ_s ∂²ζ/∂x² | C02 | C02 (D04 step 3); knowledge/primers.md: curvature of a plane curve (ch07 P169) — reminder |
| minimum of a/k + bk | C02 | C02 (D04 step 8) |
| Rayleigh–Taylor cut-off (Exercise 11.2) | C02 | C02 (N121) |
| np.roots and np.lib.scimath.sqrt | C02 | primer (in C02) |
| sympy | C02 | knowledge/primers.md: sympy (ch01 P40) — reminder |
| sympy expand, collect, subs | C02 | knowledge/primers.md: sympy expand, series, removeO, collect and subs (ch04 P117) — reminder |
| scipy.integrate.quad (mixing energy cross-check) | C02 | knowledge/primers.md: scipy.integrate.quad and dblquad (ch03 P87) — reminder |
| slider_figure | C02 | knowledge/primers.md: slider_figure (ch01 P17) — reminder |
| np.logspace | C02 | knowledge/primers.md: np.linspace and np.logspace (ch01 P06) — reminder |
| power laws and log axes | C02 | knowledge/primers.md: power laws and log–log plots (ch01 P13) — reminder |
| Python dictionaries (returned results) | C02 | knowledge/primers.md: Python dictionaries (ch01 P23) — reminder |
| assert np.allclose | C01 | knowledge/primers.md: assert np.allclose (ch01 P15) — reminder |
| Boussinesq equations (4.10, 4.86, 4.89) | C03 | R12 |
| conduction base state (11.23) | C03 | R13 |
| Rayleigh number (11.21) | C03 | C03 (N25) |
| sign of Γ: (11.21) vs Kundu ch01 vs meteorology (slip #10) | C03 | C03 (N25 box); knowledge/primers.md: inequalities under a sign change (ch01 P48) — reminder |
| geometry, centred z (Fig. 11.8 analogue) | C03 | C03 (N26) |
| decomposition (11.22), conduction profile (11.24) | C03 | C03 (N27, N28) |
| linearised perturbation equations (11.25)–(11.27), −wΓ term | C03 | C03 (N29, N30, D05) |
| scaling w ~ κ/d, buoyancy/viscous ~ Ra | C03 | C03 (N31); knowledge/primers.md: order-of-magnitude scaling (ch04 P130) — reminder |
| pressure elimination (11.28) → (11.29), horizontal Laplacian ∇_H² | C03 | C03 (N32–N34, D06) |
| divergence of a vector equation; commuting constant-coefficient operators | C03 | knowledge/primers.md: Schwarz's theorem (ch04 P121) — reminder; Ch. 10 §10.4 (pressure Poisson, D19) |
| rigid isothermal walls (11.30) | C03 | C03 (N35) |
| non-dimensional equations (11.31)–(11.33), Pr = ν/κ | C03 | C03 (N36, D07) |
| scaled variables and the chain rule | C03 | knowledge/primers.md: scaled variables and the chain rule (ch04 P133) — reminder |
| operator substitution ∂_t → σ, ∇_H² → −K² | C03 | C03 (N37, D07); knowledge/primers.md: operator elimination for linear PDEs (ch07 P177) — reminder |
| amplitude equations (11.34)–(11.38), W ≡ (Γd²/κ)ŵ | C03 | C03 (N38, N39, D07) |
| eigenvalue problem for a differential equation | C03 | primer (in C03) |
| Chebyshev–Gauss–Lobatto points and the differentiation matrix | C03 | primer (in C03) |
| boundary-row replacement and the generalised non-symmetric eigenproblem scipy.linalg.eig(A, B) | C03 | primer (in C03) |
| eigenvalues and eigenvectors | C03 | knowledge/primers.md: eigenvalues and eigenvectors (ch02 P80) — reminder |
| generalised symmetric eigenproblem (contrast) | C03 | knowledge/primers.md: generalised symmetric eigenproblem scipy.linalg.eigh(A, B) (ch10 P247) — reminder |
| collocation | C03 | knowledge/primers.md: collocation and the condition number (ch06 P162) — reminder |
| spurious eigenvalues, N-convergence filter | C03 | C03 (P259 text, spectrum figure) |
| real and imaginary parts of a complex identity | C03 | primer (in C03) |
| integration by parts with boundary terms | C03 | knowledge/primers.md: integration by parts (ch09 P218a) — reminder |
| σ real for Ra > 0 (exchange of stabilities), Exercise 11.6 | C03 | C03 (N40, N122, D08) |
| growth rate σ(K, Ra, Pr) | C03 | C03 |
| marginal pair (11.39), sixth-order (11.40), conditions (11.41) | C04 | C04 (N41–N43, D09) |
| Pr drops out at the margin | C04 | C04 (N44) |
| even and odd functions; parity under z → −z | C04 | primer (in C04) |
| cube roots of a negative number | C04 | primer (in C04) |
| characteristic roots (11.42), q₀ | C04 | C04 (N45, D10) |
| even solution, 3 × 3 determinant, slip #1 | C04 | C04 (N46, D10) |
| matrices and determinants; nonzero solution ⇔ det = 0 | C04 | knowledge/primers.md: matrices, determinants and minors (ch01 P53) — reminder; knowledge/primers.md: null space and rank–nullity theorem (ch01 P58) — reminder |
| np.linalg.det | C04 | knowledge/primers.md: np.linalg.det and np.linalg.matrix_rank (ch01 P56) — reminder |
| complex cos and cosh | C04 | C04 (D10 step 7 gloss; P45/P168 recalled) |
| scipy.optimize.brentq | C04 | knowledge/primers.md: scipy.optimize.brentq (ch03 P108) — reminder |
| scipy.optimize.minimize_scalar | C04 | knowledge/primers.md: `scipy.optimize.minimize_scalar` (ch07 P170) — reminder |
| neutral curve as the zero contour of the growth rate; critical point as its minimum | C04 | primer (in C04) |
| Ra_c = 1707.76, K_c = 3.117 (Chandrasekhar 1961) | C04 | C04 |
| rigid–free onset 1100.65 | C04 | C04 (N50) |
| gravest odd mode (Exercise 11.7) | C04 | C04 (N123) |
| Helmholtz equation in the plane; planforms as sums of cosines | C04 | primer (in C04) |
| np.meshgrid | C04 | knowledge/primers.md: np.meshgrid and the project grid layout (ch02 P76) — reminder |
| contour and streamline plots | C04 | knowledge/primers.md: plt.contour, plt.quiver and plt.streamplot (ch02 P78) — reminder |
| rolls, squares, hexagons; pattern selection | C04 | C04 (N51) |
| nonlinear equilibrium of cells | C04 | C04 (N52) |
| live widgets | C04 | knowledge/primers.md: live widgets (ch01 P47) — reminder |
| Bénard history, Marangoni convection | C03 | C03 (N24) |
| stress-free walls ⇒ W = W″ = W⁗ = 0 (11.43) | C05 | C05 (N47, D11) |
| sine modes sin nπ(z + ½), slip #2 | C05 | C05 (N48, D11) |
| integrals of sines over a period, orthogonality | C05 | knowledge/primers.md: integrals of sines and cosines over a full period (ch06 P151) — reminder |
| free–free relation (11.44) | C05 | C05 (D11) |
| quotient rule, and differentiating with respect to K² as the variable | C05 | primer (in C05) |
| dRa/dK² = 0, K_c² = π²/2, 27π⁴/4, slip #3 | C05 | C05 (N49, D11) |
| free–free growth rate (ours) | C05 | C05 (D12) |
| linear equation of state with salt | C06 | C06 (N53) |
| salt fingers and the diffusive regime | C06 | C06 (N54) |
| uniqueness of a linear boundary-value problem | C06 | primer (in C06) |
| double-diffusive marginal equations (11.45), §11.5 Ra sign | C06 | C06 (N55, D13) |
| Rs − Ra = 27π⁴/4, criterion (11.46) | C06 | C06 (N56, D13) |
| density ratio R_ρ | C06 | C06 (code cell gloss) |
| staircases, layering | C06 | C06 (N57) |
| free–free cubic dispersion relation (ours) | C06 | C06 (N54, D13 note) |
| Rayleigh's circulation criterion dΓ²/dr < 0 | C07 | R14 |
| axisymmetric NS in cylindrical coordinates (11.47), slip #12 | C07 | R15 |
| Laplacian in cylindrical coordinates | C07 | knowledge/primers.md: Laplacian in cylindrical coordinates (ch08 P186) — reminder |
| circular Couette flow (11.49), A and B | C07 | R16 |
| circulation Γ = 2πrU_θ and Kelvin's theorem for rings | C07 | C07 (N59); Ch. 5 §5.2 |
| ring interchange energy | C07 | C07 (N59) |
| Taylor problem set-up, axisymmetric disturbances | C07 | C07 (N58) |
| decomposition (11.48), linearised (11.50) | C07 | C07 (N60, N62, D14) |
| narrow-gap (small-curvature) approximation: expand in d/R and keep the leading order | C07 | primer (in C07) |
| anisotropic scaling with two length scales | C07 | knowledge/primers.md: anisotropic scaling with two length scales (ch08 P188) — reminder |
| narrow-gap equations (11.51), α = μ − 1 | C07 | C07 (N63, N64, D14) |
| Taylor number (11.52) and its inner-only form | C07 | C07 (N65, D15) |
| no-slip conditions (11.53); σ = 0 assumed | C07 | C07 (N66, N67) |
| Ta_c (11.54); μ → 1 gives Bénard | C07 | C07 (D15) |
| theory vs Taylor's experiment; Rayleigh line | C07 | C07 (N68) |
| Galerkin route, (11.92)–(11.94), slip #7 | C07 | C07 (N124) |
| wavy vortices, Görtler and Dean vortices | C07 | C07 (N69) |
| stratified parallel flow set-up | C08 | R17 |
| buoyancy frequency N² (11.56), (7.128) | C08 | R18 |
| stream function, §11.7 sign (11.57), slip #9 | C08 | R19 |
| history of the Richardson criterion | C08 | C08 (N70) |
| perturbation equations (11.55), slip #4 | C08 | C08 (N71, D16) |
| material derivative of density | C08 | Ch. 3 §3.4 ((3.5) recalled in D16 step 6) |
| normal-mode equations (11.58)–(11.60) | C08 | C08 (N72, D17) |
| Taylor–Goldstein equation (11.61) | C08 | C08 (D17) |
| conjugate symmetry of (11.61) | C08 | C08 (N73) |
| rigid lids (11.62) | C08 | C08 (N74) |
| quadratic eigenvalue problem c²M₂ + cM₁ + M₀ and its companion linearisation | C08 | primer (in C08) |
| singular point of an ODE | C08 | primer (in C08) |
| mapping an infinite domain to [−1, 1] | C08 | primer (in C08) |
| continuous spectrum (critical-layer modes) | C08 | C08 (spectrum figure notes) |
| gradient Richardson number (11.66) | C09 | R20 |
| transformation (11.63) and its derivatives | C09 | C09 (N75, D18) |
| complex powers (U − c)^{1/2} and branches | C09 | knowledge/primers.md: complex logarithm, powers and branch cuts (ch06 P155) — reminder |
| product rule | C09 | knowledge/primers.md: product rule for differentials (ch01 P38) — reminder |
| chain rule | C09 | knowledge/primers.md: chain rule (ch01 P49) — reminder |
| self-adjoint form (11.64) | C09 | C09 (N76, D18) |
| integral identity (11.65) and its imaginary part | C09 | C09 (N77, D18) |
| Miles–Howard criterion (11.67) | C09 | C09 (D18) |
| integrating on a Chebyshev grid (Clenshaw–Curtis weights) | C09 | primer (in C09) |
| Ri < ¼ necessary, not sufficient; tanh layer results | C09 | C09 (N78) |
| Simpson and trapezoid (contrast) | C09 | knowledge/primers.md: `scipy.integrate.simpson` and `np.trapezoid` (ch09 P203) — reminder |
| transformation F = ψ̂/(U − c) (11.68) | C10 | C10 (N79, D19) |
| divergence form, slip #11; weight Q | C10 | C10 (N80, D19) |
| real and imaginary parts (11.69)–(11.70) | C10 | C10 (N81, D19) |
| c_r inside the velocity range (11.71) | C10 | C10 (N82, D19) |
| completing the square into a circle (x − a)² + y² ≤ R² | C10 | primer (in C10) |
| completing the square (tensors, recalled) | C10 | knowledge/primers.md: completing the square for tensors (ch04 P129) — reminder |
| Howard's semicircle, (11.72), growth bound, slip #8 | C10 | C10 (D19) |
| dimensionless Navier–Stokes, Re = U₀L/ν (11.73) | C11 | R21 |
| six Reynolds numbers (which length?) | C11 | knowledge/primers.md: six Reynolds numbers (which length?) (ch09 P200) — reminder |
| stream function, §11.8 sign; û = φ′, v̂ = −ikφ | C11 | R22 |
| viscosity can destabilise (named) | C11 | C11 (N83) |
| perturbation equations (11.74)–(11.77) | C11 | C11 (N84–N87, D20) |
| Squire's theorem and transformation (11.78) | C11 | C11 (N88) |
| consequences of Squire (oblique waves) | C11 | C11 (N89) |
| Orr–Sommerfeld equation (11.79) | C11 | C11 (D21) |
| no-slip conditions (11.80) | C11 | C11 (N90) |
| Orszag's benchmark c = 0.23752649 + 0.00373967i | C11 | C11 (code cell, `reference/ch11/benchmarks.json`) |
| reading reference data from a file | C11 | knowledge/primers.md: reading reference data from a file and interpolating (np.interp) (ch10 P251) — reminder |
| verification versus validation | C11 | knowledge/primers.md: verification versus validation (ch10 P253) — reminder |
| Rayleigh equation (11.81), singular limit | C12 | C12 (N91) |
| inviscid conditions (11.82); pairs broken by viscosity | C12 | C12 (N92) |
| inflection point | C12 | knowledge/primers.md: inflection point (ch09 P212) — reminder |
| Rayleigh's theorem (11.83)–(11.84) | C12 | C12 (D22) |
| Fjørtoft's theorem (11.85)–(11.86), U_I | C12 | C12 (N93) |
| six profiles, verdict table | C12 | C12 (N94) |
| sin y counterexample | C12 | C12 (N95) |
| critical layer | C12 | C12 (N96) |
| Kelvin's cat's eye (11.87) | C12 | C12 (N97) |
| sech and (tanh)′ = sech² | C12 | knowledge/primers.md: sech, arccosh and (tanh)′ = sech² (ch09 P215) — reminder |
| np.sign and np.nonzero | C12 | knowledge/primers.md: `np.sign` and `np.nonzero`: where a curve crosses zero (ch07 P180) — reminder |
| piecewise-linear shear layer (Exercise 11.11) | C12 | C12 (N126) |
| sinuous and varicose modes (11.95) | C12 | C12 (N127) |
| Bickley jet | C13 | R23 |
| boundary layers with pressure gradients, neutral loops | C13 | R24 |
| viscosity stabilising and destabilising (named) | C13 | C13 (N98) |
| plane Poiseuille critical point Re_c = 5772.22 | C13 | C13 |
| Tollmien–Schlichting wave | C13 | C13 |
| tanh layer neutral curve, sech neutral mode | C13 | C13 (N99) |
| plane Couette linearly stable | C13 | C13 (N100) |
| pipe flow (named) | C13 | C13 (N101) |
| Table 11.1 recomputed | C13 | C13 (N102) |
| pandas DataFrame for a table | C13 | C13 (gloss: `pd.DataFrame(list_of_dicts)` shows rows and columns) |
| spreading mixing layer (named) | C13 | C13 (N103) |
| caching expensive runs | C13 | knowledge/primers.md: caching expensive runs (np.savez and a parameter key) (ch10 P252) — reminder |
| Tollmien profile, slip #6 | C13 | C13 (N105) |
| Blasius profile and δ* | C13 | Ch. 9 §9.3 (recalled in N105, N108) |
| Schubauer & Skramstad (named) | C13 | C13 (N106) |
| neutral curve in frequency F = ων/U∞² | C13 | C13 (N107) |
| Blasius critical Re_δ* (Thomas; Jordinson) | C13 | C13 (N108) |
| non-parallel corrections, transient growth (named) | C13 | C13 (N109) |
| disturbed NS in index form (11.96) | C14 | C14 (N128) |
| index notation and the divergence theorem | C14 | Ch. 2 §2.13 (Gauss; recalled in D23) |
| the integral of an x-derivative of a periodic function over one period is zero | C14 | primer (in C14) |
| fundamental theorem of calculus | C14 | knowledge/primers.md: fundamental theorem of calculus (ch02 P84) — reminder |
| disturbance kinetic-energy equation (11.88) | C14 | C14 (D23) |
| Reynolds stress −⟨uv⟩ and production | C14 | C14 (N104) |
| mean of a product of real parts | C14 | knowledge/primers.md: mean of a product of real parts (ch07 P178) — reminder |
| viscous dissipation Λ ≥ 0 | C14 | C14 (D23 step 9) |
| stratified energy equation (Exercise 11.10, named) | C14 | C14 (N125) |
| rectified fluxes, saturation, annulus (named) | C14 | C14 (N110) |
| sequence of instabilities in transition | C14 | C14 (N111) |
| np.polynomial.chebyshev fit (re-interpolation) | C14 | C14 (gloss in the from-scratch cell) |
| deterministic chaos (sensitivity, aperiodicity, broadband) | C15 | C15 (N112) |
| np.fft.rfft spectrum | C15 | knowledge/primers.md: `np.fft.rfft` and `np.fft.rfftfreq` (ch07 P181) — reminder |
| phase space, trajectory, degrees of freedom; pendulum (11.89) | C15 | C15 (N113) |
| Jacobian matrix of a nonlinear ODE system and the stability of its fixed points | C15 | primer (in C15) |
| multivariable first-order Taylor expansion | C15 | knowledge/primers.md: multivariable first-order Taylor expansion (ch03 P98) — reminder |
| attractors, limit cycles, bifurcations; Hopf normal form | C15 | C15 (N114) |
| Lorenz truncation (11.90), slip #9 | C15 | C15 (N115) |
| Galerkin truncation: keep a few modes and project with orthogonality of sines | C15 | primer (in C15) |
| vorticity–stream function form of 2-D Boussinesq flow | C15 | C15 (D24 steps 2–4); Ch. 5 §5.1 |
| Jacobian J(ψ, f) = ψ_x f_z − ψ_z f_x (advection) | C15 | C15 (D24 step 4) |
| Lorenz system (11.91) | C15 | C15 (D24) |
| purely imaginary roots of a cubic λ³ + a₂λ² + a₁λ + a₀: exactly when a₂a₁ = a₀ | C15 | primer (in C15) |
| fixed points, pitchfork at r = 1, Hopf point r_H | C15 | C15 (N116, D25) |
| RK4 by hand | C15 | knowledge/primers.md: RK4 by hand (ch03 P95) — reminder |
| solve_ivp options (DOP853, rtol, t_eval) | C15 | knowledge/primers.md: solve_ivp options: t_eval, dense_output, events, backward integration (ch03 P94) — reminder |
| Lyapunov exponent: the slope of log(separation) against time | C15 | primer (in C15) |
| np.polyfit | C15 | C15 (gloss in the P277 demo: least-squares straight line) |
| matplotlib 3-D line plots | C15 | primer (in C15) |
| plotly 3-D (explainers, alternatives) | C15 | knowledge/primers.md: plotly 3-D surface (ch01 P41) — reminder |
| iterated maps: fixed point, stability ∣f′(x*)∣ < 1, cobweb diagram | C15 | primer (in C15) |
| period doubling, Feigenbaum ratio | C15 | C15 (N117) |
| logistic map (Exercise 11.14) | C15 | C15 (N129) |
| quasi-periodic route, Landau, Ruelle–Takens (named) | C15 | C15 (N118) |
| predictability and ensembles (named) | C15 | C15 (N119) |
| exercises 11.3–11.5 (not solved) | C02 | C02 (pointer, S01) |
| literature cited | C15 | C15 (pointer, S02) |
| matplotlib figures | C01 | knowledge/primers.md: matplotlib figures (ch01 P01) — reminder |
| numpy arrays and broadcasting | C01 | knowledge/primers.md: numpy broadcasting (ch02 P77) — reminder |
| f-strings | C01 | knowledge/primers.md: f-strings (ch01 P04) — reminder |
| functions as arguments and lambda | C02 | knowledge/primers.md: functions as arguments and lambda (ch01 P29) — reminder |
| tuple unpacking | C02 | knowledge/primers.md: tuple unpacking (ch01 P14) — reminder |
| explainer tabs (Explain, Derivation, Code) | C01 | C01 (gloss in the explainer cell: "the Explain tab works every number out with your settings") |

---

## Part F — derivation storyboards

Builders copy these word for word into `nb.derivation(key, title, goal=…, start=(tex, plain), plan=[…], uses=[…],
steps=[dict(did, tex, why, plain)], result=(tex, plain), interpret=…, check=…, check_src=…)` and into the explainer's
`derivations: [...]` (phones may shorten *why* to its first sentence; `live`, `set` and `watch` are the explainer's and are listed
in Part B). Every step is one move; *why* names the rule and says why we make it (≤ 35 words); *did* ≤ 8 words. The book's own moves
were read on the rendered pages (p503 for D01; p505–p506 for D02–D03; D04 is ours (Exercise 11.1 gives only the result); p511–p514
for D05–D07; D08 is ours (Exercise 11.6 gives an outline); p515–p516 for D09–D10; p516–p518 for D11; D12 is ours; p521 for D13 (the
book says "repeat the derivation"); D14 is ours (the book defers to Chandrasekhar 1961); p526 for D15; p529–p531 for D16–D17;
p531–p532 for D18; p533–p535 for D19; p535–p536 for D20; p537 for D21 ("this effort yields"); p538–p539 for D22; p546–p547 for D23
(Exercise 11.13); p553–p555 for D24 ("Lorenz finally obtained") and D25 (the book says only "if r is large")). The gaps listed in
`analysis/ch11.md` §2b are filled and the notebook says so ("the book skips this move; we add it"). Every equation named by number
is written out. Colours: growing rose, decaying teal, marginal purple, buoyancy blue, shear/kinetic orange, salt amber. No line of
this part starts with a table bar; absolute values are written with \lvert \rvert or in words. The ★★★ derivations (D10, D14, D18,
D19, D24) carry a `check_src` cell, every line commented; D03, D08, D11, D21, D23 carry an optional one.

### D01 · The two forms of the normal mode (11.1): $\sigma=-\mathrm i\lvert\mathbf K\rvert c$, $\sigma_r=\lvert\mathbf K\rvert c_i$, $\sigma_i=-\lvert\mathbf K\rvert c_r$ — ★, 5 steps, in C01 (notebook · `normal_mode_growth`)
- **Goal.** Show that the two ways of writing a normal mode in (11.1) are the same mode, and read off how growth and travel
  appear in each.
- **Start.** $u=\hat u(z)\exp\{\mathrm ikx+\mathrm imy+\sigma t\}=\hat u(z)\exp\{\mathrm i\lvert\mathbf K\rvert(\mathbf e_K\cdot
  \mathbf x-ct)\}$ (11.1), $\mathbf K=(k,m,0)$, $\mathbf e_K=\mathbf K/\lvert\mathbf K\rvert$ — *in words:* a wave in x, y whose
  amplitude changes in time, written with a growth rate σ or with a complex wave speed c.
- **Plan.** (1) Expand the second exponent. (2) Equate it with the first. (3) Split σ and c into real and imaginary parts.
- **Tools.** Euler's formula and complex numbers (P45 reminder); complex amplitudes and Re{·} (P176 reminder).
- **Assumptions.** k, m real (the disturbance stays bounded in x and y) — used in step 4.
- **Steps.**
  1. *did:* Expand the dot product · *tex:* $\lvert\mathbf K\rvert\,\mathbf e_K\cdot\mathbf x=\mathbf K\cdot\mathbf x=kx+my$ ·
     *why:* e_K is K divided by its length, so the length times e_K is K itself; we want both exponents in the same variables. ·
     *plain:* the second form's space part is the same wave in x and y.
  2. *did:* Rewrite the second exponent · *tex:* $\mathrm i\lvert\mathbf K\rvert(\mathbf e_K\cdot\mathbf x-ct)=\mathrm i(kx+my)
     -\mathrm i\lvert\mathbf K\rvert c\,t$ · *why:* distribute i∣K∣ over the bracket and use step 1; this isolates the part that
     depends on t. · *plain:* the exponent is 'wave in space' plus 'something times t'.
  3. *did:* Equate the time parts · *tex:* $\sigma=-\mathrm i\lvert\mathbf K\rvert c$ · *why:* both forms must give the same
     function of t for every t, so the coefficients of t agree (the space parts already do). · *plain:* a growth rate is a complex
     wave speed times −i∣K∣.
  4. *did:* Insert c = c_r + i c_i · *tex:* $\sigma=\lvert\mathbf K\rvert c_i-\mathrm i\lvert\mathbf K\rvert c_r$ · *why:* multiply
     out, using i × i = −1; real and imaginary parts are then visible (∣K∣ is real by the assumption). · *plain:* the real part of σ
     comes from c_i, the imaginary part from c_r.
  5. *did:* Read growth and travel · *tex:* $\sigma_r=\lvert\mathbf K\rvert c_i,\quad\sigma_i=-\lvert\mathbf K\rvert c_r$ · *why:*
     $e^{\sigma t}=e^{\sigma_rt}e^{\mathrm i\sigma_it}$: the first factor is the amplitude, the second a phase; a crest kx + σ_i t = const
     moves at −σ_i/k = c_r. · *plain:* c_i > 0 means growth; c_r is the speed of the crests.
- **Result.** $\sigma=-\mathrm i\lvert\mathbf K\rvert c$; growth rate $\sigma_r=\lvert\mathbf K\rvert c_i$, crest speed $c_r=-\sigma_i/
  \lvert\mathbf K\rvert$ — *in words:* "unstable" means σ_r > 0, equivalently c_i > 0.
- **Check.** Units: σ [1/s] = ∣K∣ [1/m] × c [m/s] ✓. Numbers: K = 2, c = 1 + 0.5i ⇒ σ = 1 − 2i (`ST.sigma_from_c`) ✓.
- **What it means.** Stationary patterns (cells) have σ_i = 0; travelling waves have σ_i ≠ 0. Fails if k or m is complex (spatial
  instability, N03).
- **Traps.** The minus sign in σ_i: a wave moving toward +x has σ_i < 0; only the real part of (11.1) is physical.

### D02 · The linearised interface conditions (11.9) and (11.13) from (11.6)–(11.8) and (11.10)–(11.12) — ★★, 12 steps, in C02 (notebook · `kelvin_helmholtz_boundary`)
- **Goal.** Turn the exact conditions on the moving interface into simple linear conditions on the flat plane z = 0.
- **Start.** $\mathbf n\cdot\nabla\tilde\phi_1=\mathbf n\cdot\mathbf U_s=\mathbf n\cdot\nabla\tilde\phi_2$ on $z=\zeta$ (11.6), $p_1=p_2$
  on $z=\zeta$ (11.7), and the Bernoulli equations $\frac{\partial\tilde\phi_j}{\partial t}+\tfrac12\lvert\nabla\tilde\phi_j\rvert^2+
  \frac{p_j}{\rho_j}+gz=C_j$ (11.10) — *in words:* fluid on each side moves with the interface, and the pressure is continuous.
- **Plan.** (1) Write the normal and the interface velocity. (2) Do the dot products and cancel the common factor. (3) Drop
  products of small quantities and move the conditions to z = 0. (4) Same for the pressures via Bernoulli.
- **Tools.** Level sets and their normals (R05, P75); unsteady Bernoulli (R06); orders of smallness (P68); Taylor transfer of a
  boundary condition (P166).
- **Assumptions.** Small slope ∂ζ/∂x ≪ 1 and small disturbance potentials (steps 6, 11); no surface tension (step 8);
  irrotational layers (R01).
- **Steps.**
  1. *did:* Describe the interface as a level set · *tex:* $f=z-\zeta(x,t)=0,\quad\nabla f=-\frac{\partial\zeta}{\partial x}\mathbf e_x
     +\mathbf e_z$ · *why:* the interface is where f vanishes; the gradient of a level-set function is normal to the level set
     (P75). · *plain:* we have a vector that points straight out of the interface.
  2. *did:* Normalise it · *tex:* $\mathbf n=\frac{-(\partial\zeta/\partial x)\mathbf e_x+\mathbf e_z}{\sqrt{1+(\partial\zeta/\partial
     x)^2}}$ · *why:* dividing by the length gives the unit normal used in (11.6); the square root is the length of ∇f. ·
     *plain:* the unit normal tilts with the local slope.
  3. *did:* Choose the interface velocity · *tex:* $\mathbf U_s=\frac{\partial\zeta}{\partial t}\mathbf e_z$ · *why:* only the normal
     part of U_s enters (11.6), so we may let points of the interface move vertically; it keeps the algebra simple. · *plain:* we
     follow the interface up and down at fixed x.
  4. *did:* Take the three dot products · *tex:* $\frac{-\zeta_x\partial_x\tilde\phi_1+\partial_z\tilde\phi_1}{\sqrt{1+\zeta_x^2}}=
     \frac{\zeta_t}{\sqrt{1+\zeta_x^2}}=\frac{-\zeta_x\partial_x\tilde\phi_2+\partial_z\tilde\phi_2}{\sqrt{1+\zeta_x^2}}$ · *why:* the
     velocity is ∇φ̃ = (∂_xφ̃, ∂_zφ̃), and dot products add componentwise; this is (11.8) written out. · *plain:* the normal
     velocity of each fluid equals the normal speed of the interface.
  5. *did:* Cancel the common square root · *tex:* $-(U_1+\partial_x\phi_1)\zeta_x+\partial_z\phi_1=\zeta_t=-(U_2+\partial_x
     \phi_2)\zeta_x+\partial_z\phi_2$ · *why:* the same nonzero factor divides all three members, so multiplying through by it is
     allowed; then insert (11.2), ∂_xφ̃_j = U_j + ∂_xφ_j. The book skips why the root disappears. · *plain:* the exact kinematic
     condition, still on the wavy interface.
  6. *did:* Drop the quadratic term · *tex:* $-U_1\zeta_x+\partial_z\phi_1=\zeta_t=-U_2\zeta_x+\partial_z\phi_2\quad(z=\zeta)$ ·
     *why:* ∂_xφ_j ζ_x is a product of two small quantities, one order smaller than the rest (P68); linear theory keeps first order
     only. · *plain:* each stream carries the slope; the disturbance adds a vertical velocity.
  7. *did:* Move the condition to z = 0 · *tex:* $-U_1\frac{\partial\zeta}{\partial x}+\frac{\partial\phi_1}{\partial z}=\frac{\partial
     \zeta}{\partial t}=-U_2\frac{\partial\zeta}{\partial x}+\frac{\partial\phi_2}{\partial z}\ \text{on }z=0$ (11.9) · *why:* Taylor
     transfer (P166): $\partial_z\phi(x,\zeta)=\partial_z\phi(x,0)+\zeta\,\partial_z^2\phi(x,0)+\dots$; the correction is a product of
     small quantities. · *plain:* the linear kinematic condition on the flat plane.
  8. *did:* Equate the two pressures · *tex:* $\rho_1\big(C_1-\partial_t\tilde\phi_1-\tfrac12\lvert\nabla\tilde\phi_1\rvert^2-gz\big)=
     \rho_2\big(C_2-\partial_t\tilde\phi_2-\tfrac12\lvert\nabla\tilde\phi_2\rvert^2-gz\big)\ (z=\zeta)$ (11.11) · *why:* solve each
     Bernoulli equation (11.10) for p_j and use (11.7) (no surface tension, so no pressure jump). · *plain:* the pressure from
     above equals the pressure from below at the interface.
  9. *did:* Write the undisturbed balance · *tex:* $\rho_1(C_1-\tfrac12U_1^2)=\rho_2(C_2-\tfrac12U_2^2)$ (11.12) · *why:* with
     φ_j = ζ = 0 (the plain streams) (11.11) must hold too; it fixes how the two Bernoulli constants are related. · *plain:* the
     flat interface is in pressure balance.
  10. *did:* Subtract (11.11) from (11.12) · *tex:* $\rho_1\big(\partial_t\phi_1+U_1\partial_x\phi_1+\tfrac12\lvert\nabla\phi_1\rvert^2+
      g\zeta\big)=\rho_2\big(\partial_t\phi_2+U_2\partial_x\phi_2+\tfrac12\lvert\nabla\phi_2\rvert^2+g\zeta\big)$ · *why:* both are
      equalities between the same two sides, so their difference is one; expand $\lvert\nabla\tilde\phi_j\rvert^2=U_j^2+2U_j\partial_x
      \phi_j+\lvert\nabla\phi_j\rvert^2$ — the U_j² cancels. · *plain:* only the disturbance's part of the pressure is left.
  11. *did:* Drop the quadratic terms · *tex:* $\rho_1(\partial_t\phi_1+U_1\partial_x\phi_1+g\zeta)=\rho_2(\partial_t\phi_2+U_2\partial_x
      \phi_2+g\zeta)\ (z=\zeta)$ · *why:* ½∣∇φ_j∣² is quadratic in the disturbance (P68). Note: it is $U\partial_x\phi$ that
      survives, not ½(∂_xφ)². · *plain:* unsteadiness, advection by the stream and gravity balance.
  12. *did:* Evaluate on z = 0 · *tex:* $\rho_1\Big(\frac{\partial\phi_1}{\partial t}+U_1\frac{\partial\phi_1}{\partial x}+g\zeta\Big)
      =\rho_2\Big(\frac{\partial\phi_2}{\partial t}+U_2\frac{\partial\phi_2}{\partial x}+g\zeta\Big)\ \text{on }z=0$ (11.13) · *why:*
      the same Taylor transfer as step 7; moving from z = ζ to z = 0 changes each term only at second order. · *plain:* the linear
      dynamic condition.
- **Result.** (11.9) and (11.13) — *in words:* two linear conditions on z = 0 linking the interface shape ζ to the two potentials.
- **Check.** Units: (11.9) every term m/s ✓; (11.13) every term Pa ✓. Limit U₁ = U₂ = 0: Ch. 7's (7.18) and (7.21) ✓.
  `ch11.kh_residuals` gives residuals ≤ 1e-12 for the computed mode.
- **What it means.** The streams enter only through U_j∂/∂x — advection of the disturbance by each stream; that is where shear
  will act. Fails when the slope is not small (breaking billows, A2).
- **Traps.** Keeping ½(∂_xφ)² instead of U∂_xφ; forgetting that the square-root factor multiplies all three members; thinking the
  shift to z = 0 is an approximation of the same order as the linearisation (it is one order higher).

### D03 · The Kelvin–Helmholtz dispersion relation (11.18), the instability inequality and the limits (11.19), (11.20) — ★★, 14 steps, in C02 (notebook · `kelvin_helmholtz_boundary`)
- **Goal.** Find the complex wave speed c of a small interface wave between two streams, and when it has a growing part.
- **Start.** $\nabla^2\phi_j=0$ (11.3), decay (11.4)–(11.5), (11.9), (11.13) and the trial $\phi_j=A_j(z)\exp\{\mathrm ik(x-ct)\}$
  (11.14) — *in words:* the linear problem of D02 and a wave-shaped guess.
- **Plan.** (1) Solve Laplace's equation for the depth structure. (2) Use the kinematic condition to fix the amplitudes. (3) Put
  them into the dynamic condition: a quadratic for c. (4) Solve and read off when c is complex.
- **Tools.** Linear ODE by an exponential trial (P44, R07); complex algebra (P153); the quadratic formula with complex roots (P159);
  inequalities (P48); sympy (P40) for the check.
- **Assumptions.** k > 0 real (step 3); infinitely deep layers (step 3); ζ₀ ≠ 0 (step 9).
- **Steps.**
  1. *did:* Put the trial into Laplace · *tex:* $A_j''-k^2A_j=0$ · *why:* ∂²/∂x² of $e^{\mathrm ik(x-ct)}$ gives −k²; the common
     exponential is never zero, so divide it out. · *plain:* the depth structure obeys a simple ODE.
  2. *did:* Solve the ODE · *tex:* $A_j=a_je^{kz}+b_je^{-kz}$ · *why:* the trial $e^{mz}$ gives m² = k² (R07). · *plain:* each layer's
     disturbance grows or decays exponentially away from the interface.
  3. *did:* Keep the decaying parts · *tex:* $\phi_1=A_-e^{\mathrm ik(x-ct)-kz},\quad\phi_2=A_+e^{\mathrm ik(x-ct)+kz}$ (11.15) · *why:*
     (11.4) needs φ₁ → 0 as z → +∞ (only e^{−kz}), (11.5) needs φ₂ → 0 as z → −∞ (only e^{+kz}); k > 0. · *plain:* the disturbance
     lives near the interface.
  4. *did:* Give the interface the same form · *tex:* $\zeta=\zeta_o e^{\mathrm ik(x-ct)}$, $\partial_x\to\mathrm ik$, $\partial_t\to
     -\mathrm ikc$ · *why:* every term must carry the same exponential for the conditions to hold at every x and t. · *plain:* the
     interface is the same travelling wave.
  5. *did:* Insert into the kinematic condition (11.9) · *tex:* $-\mathrm iU_1k\zeta_o-kA_-=-\mathrm ikc\zeta_o=-\mathrm iU_2k\zeta_o
     +kA_+$ (11.16) · *why:* at z = 0, ∂_zφ₁ = −kφ₁ and ∂_zφ₂ = +kφ₂; then divide by the common exponential. · *plain:* three
     numbers that must be equal.
  6. *did:* Solve for the amplitudes · *tex:* $A_-=-\mathrm i(U_1-c)\zeta_o,\quad A_+=\mathrm i(U_2-c)\zeta_o$ · *why:* each equality of
     (11.16) has one unknown; divide by k. · *plain:* the potentials are fixed by the interface's height and the speeds.
  7. *did:* Insert into the dynamic condition (11.13) · *tex:* $\rho_1(-\mathrm ikcA_-+\mathrm ikU_1A_-+g\zeta_o)=\rho_2(-\mathrm
     ikcA_++\mathrm ikU_2A_++g\zeta_o)$ (11.17) · *why:* ∂_t → −ikc, ∂_x → ik on each potential at z = 0. · *plain:* the pressure
     balance for the wave.
  8. *did:* Substitute A± · *tex:* $\rho_1\big(k(U_1-c)^2+g\big)\zeta_o=\rho_2\big(-k(U_2-c)^2+g\big)\zeta_o$ · *why:*
     $\mathrm ik(U_1-c)\cdot(-\mathrm i)(U_1-c)=k(U_1-c)^2$ and $\mathrm ik(U_2-c)\cdot\mathrm i(U_2-c)=-k(U_2-c)^2$ (i² = −1). The
     book prints this as $(-\mathrm ikc+\mathrm ikU)^2=-k^2(U-c)^2$. · *plain:* the shear appears squared.
  9. *did:* Divide by kζ_o and rearrange · *tex:* $\rho_1(U_1-c)^2+\rho_2(U_2-c)^2=\frac gk(\rho_2-\rho_1)$ · *why:* ζ_o ≠ 0 (a wave
     exists) and k > 0; move the gravity terms to one side. · *plain:* kinetic-energy-like terms balance gravity's restoring term.
  10. *did:* Expand in powers of c · *tex:* $(\rho_1+\rho_2)c^2-2(\rho_1U_1+\rho_2U_2)c+\rho_1U_1^2+\rho_2U_2^2-\frac gk(\rho_2-\rho_1)=0$ ·
      *why:* multiply out the squares; a quadratic in standard form is ready for the formula. · *plain:* c solves a quadratic.
  11. *did:* Apply the quadratic formula · *tex:* $c=\frac{\rho_1U_1+\rho_2U_2\pm\sqrt{\Delta}}{\rho_1+\rho_2}$ · *why:* P159, with
      $\Delta=(\rho_1U_1+\rho_2U_2)^2-(\rho_1+\rho_2)\big(\rho_1U_1^2+\rho_2U_2^2-\frac gk(\rho_2-\rho_1)\big)$. · *plain:* two wave
      speeds around a density-weighted mean.
  12. *did:* Simplify the discriminant · *tex:* $\frac{\Delta}{(\rho_1+\rho_2)^2}=\frac{\rho_2-\rho_1}{\rho_2+\rho_1}\frac gk-
      \frac{\rho_1\rho_2(U_2-U_1)^2}{(\rho_1+\rho_2)^2}$ — this is (11.18) · *why:* the identity $(\rho_1U_1+\rho_2U_2)^2-(\rho_1+\rho_2)
      (\rho_1U_1^2+\rho_2U_2^2)=-\rho_1\rho_2(U_1-U_2)^2$ (expand both; the book skips it). · *plain:* gravity term minus shear term.
  13. *did:* Ask when the root is imaginary · *tex:* $g(\rho_2^2-\rho_1^2)<k\rho_1\rho_2(U_2-U_1)^2$ · *why:* c is complex (one root with
      c_i > 0) exactly when the bracket is negative; multiply by $k(\rho_1+\rho_2)^2>0$, which keeps the inequality's direction (P48).
      · *plain:* shear wins for strong shear, small density step or short waves.
  14. *did:* Read off the limits · *tex:* $U_1=U_2=0:\ c=\pm\big[\frac{\rho_2-\rho_1}{\rho_2+\rho_1}\frac gk\big]^{1/2}$ (11.19);
      $\rho_1=\rho_2:\ c=\frac{U_2+U_1}2\pm\mathrm i\frac{U_2-U_1}2$ (11.20) · *why:* set the speeds or the density difference to zero
      in (11.18); the mean becomes (U₁ + U₂)/2. · *plain:* Ch. 7's interface waves, and a vortex sheet unstable at every k.
- **Result.** $c=\frac{\rho_2U_2+\rho_1U_1}{\rho_2+\rho_1}\pm\Big[\frac{\rho_2-\rho_1}{\rho_2+\rho_1}\frac gk-\frac{\rho_2\rho_1}{(\rho_2+
  \rho_1)^2}(U_2-U_1)^2\Big]^{1/2}$ (11.18), unstable iff $g(\rho_2^2-\rho_1^2)<k\rho_1\rho_2(U_2-U_1)^2$ — *in words:* gravity holds the
  interface, shear tears it; short waves always lose.
- **Check.** Units: every term under the root m²/s² ✓. Tiny example (ρ = 1, 3; g = 10; k = 1; U₁ = 6): c = 1.5 ± 1.3229i ✓. Galilean:
  adding V to both U's shifts c by V ✓. `ch11.kh_sympy()["identity"] == 0` (optional `check_src`: build the quadratic from (11.16)–(11.17)
  with sympy, solve, and compare with (11.18)).
- **What it means.** The mean term says where the wave drifts; the root says whether it grows. For air over water, gravity alone
  cannot hold short waves (N18); surface tension must (N120). Fails for finite slopes (roll-up, A2) and for thick shear layers
  (continuous profiles: C08).
- **Traps.** Choosing e^{+kz} above (it grows upward); dropping the i in (ik)² and getting the shear term's sign wrong (verifier
  note ii); dividing by ζ_o without saying it is nonzero; multiplying an inequality by a quantity of unknown sign.

### D04 · Kelvin–Helmholtz with a lower layer of depth h and surface tension (Exercise 11.1), and the minimum wind $\Delta U_{\min}^2=2\sqrt{g\Delta\rho\,\sigma_s}(\rho_1+\rho_2)/(\rho_1\rho_2)$ — ★★, 9 steps, in C02 (notebook · `kelvin_helmholtz_boundary`)
- **Goal.** Add a bottom at depth h and surface tension to D03 and find the weakest shear that can make waves grow.
- **Start.** D03's (11.15)–(11.17) with the lower potential now meeting a wall: ∂φ₂/∂z = 0 at z = −h — *in words:* the same two
  streams, a finite lower layer, a stretched interface.
- **Plan.** (1) Replace the lower exponential by a cosh that meets the wall. (2) Add the surface-tension pressure jump. (3) Repeat
  D03's algebra. (4) Minimise the restoring terms over k.
- **Tools.** D03's moves; hyperbolic functions (P168); the curvature of a plane curve (P169) for the pressure jump; minimising
  a/k + bk (step 8).
- **Assumptions.** Small slopes (curvature ≈ ∂²ζ/∂x², step 4); constant surface tension σ_s [N/m].
- **Steps.**
  1. *did:* Choose the lower mode · *tex:* $\phi_2=A_+\cosh k(z+h)\,e^{\mathrm ik(x-ct)}$ · *why:* cosh k(z + h) solves A″ = k²A
     (R07) and its z-derivative k sinh k(z + h) vanishes at z = −h: no flow through the bottom. · *plain:* the lower disturbance
     fills the layer and stops at the floor.
  2. *did:* Apply the lower kinematic condition · *tex:* $-\mathrm iU_2k\zeta_o+kA_+\sinh kh=-\mathrm ikc\zeta_o$ · *why:* as D03 step
     5, now with ∂_zφ₂ = kA₊ sinh kh at z = 0. · *plain:* the lower fluid follows the interface.
  3. *did:* Solve for A₊ · *tex:* $A_+=\frac{\mathrm i(U_2-c)\zeta_o}{\sinh kh},\quad\phi_2(z=0)=\mathrm i(U_2-c)\zeta_o\coth kh\,
     e^{\mathrm ik(x-ct)}$ · *why:* divide by k sinh kh ≠ 0; cosh kh/sinh kh = coth kh. A₋ is unchanged from D03 step 6. · *plain:*
     depth enters only through coth kh.
  4. *did:* Add the surface-tension jump · *tex:* $p_1-p_2=\sigma_s\frac{\partial^2\zeta}{\partial x^2}=-\sigma_sk^2\zeta$ · *why:* a
     stretched surface pushes toward its concave side with pressure σ_s × curvature (P169); for a small slope the curvature is
     ∂²ζ/∂x². · *plain:* a crest is pushed down by the skin.
  5. *did:* Rewrite the dynamic condition · *tex:* $\rho_1(\partial_t\phi_1+U_1\partial_x\phi_1+g\zeta)-\rho_2(\partial_t\phi_2+U_2
     \partial_x\phi_2+g\zeta)=\sigma_sk^2\zeta$ · *why:* linearised Bernoulli gives $p_j=-\rho_j(\partial_t\phi_j+U_j\partial_x\phi_j+
     g\zeta)$ plus constants; insert into step 4 at z = 0. · *plain:* (11.13) plus the skin's push.
  6. *did:* Substitute the amplitudes · *tex:* $\rho_1(U_1-c)^2+\rho_2\coth kh\,(U_2-c)^2=\frac gk(\rho_2-\rho_1)+\sigma_sk$ · *why:* the
     same products as D03 step 8 (the lower one carries coth kh); divide by kζ_o. · *plain:* surface tension joins gravity on the
     restoring side.
  7. *did:* Solve the quadratic · *tex:* $c=\frac{\rho_1U_1+\rho_2U_2\coth kh}{\rho_1+\rho_2\coth kh}\pm\Big[\frac{\frac gk(\rho_2-\rho_1)+
     \sigma_sk}{\rho_1+\rho_2\coth kh}-\frac{\rho_1\rho_2(U_1-U_2)^2\coth kh}{(\rho_1+\rho_2\coth kh)^2}\Big]^{1/2}$ · *why:* D03 steps
     10–12 with ρ₂ replaced by ρ₂ coth kh (same identity). · *plain:* Exercise 11.1's formula.
  8. *did:* Minimise the restoring side (deep) · *tex:* $\min_k\Big[\frac gk\Delta\rho+\sigma_sk\Big]=2\sqrt{g\Delta\rho\,\sigma_s}\ \text{at}\
     k_*=\sqrt{g\Delta\rho/\sigma_s}$ · *why:* a/k + bk has derivative −a/k² + b = 0 at k = √(a/b), value 2√(ab); with h → ∞ (coth → 1)
     instability needs $\frac{\rho_1\rho_2\Delta U^2}{\rho_1+\rho_2}>\frac gk\Delta\rho+\sigma_sk$. · *plain:* the weakest spot of the
     restoring forces is at one wavelength.
  9. *did:* Read off the minimum wind · *tex:* $\Delta U_{\min}^2=\frac{2\sqrt{g\Delta\rho\,\sigma_s}(\rho_1+\rho_2)}{\rho_1\rho_2}$ ·
     *why:* step 8's minimum equals the shear side there; solve for ΔU. With h → ∞ step 7 reduces to (11.18) plus σ_s k. · *plain:*
     no shear below this value can make any wave grow.
- **Result.** Exercise 11.1's c (step 7) and $\Delta U_{\min}=6.70$ m/s at λ* = 2π/k* = 1.73 cm for air over water (σ_s = 0.074 N/m)
  — *in words:* gravity guards the long waves, surface tension the short ones; between them a minimum shear exists.
- **Check.** Units: σ_s k/(ρ) [N/m × 1/m ÷ kg/m³ = m²/s²] ✓. Limits: σ_s = 0, h → ∞ gives (11.18) ✓; ρ₁ = ρ₂, σ_s = 0 gives (11.20) ✓.
  Numbers: `ch11.kh_min_shear(1.2, 1000.0)` ✓.
- **What it means.** KH alone needs 6.7 m/s of wind; real ripples start at lower winds because turbulent pressure fluctuations
  (not in this model) do the work. Fails when the wind profile is not a jump (a boundary layer in the air).
- **Traps.** Using e^{kz} for the lower layer (it ignores the bottom); the sign of the pressure jump (the skin pushes a crest down:
  it stabilises); forgetting that coth kh → 1 only when kh ≳ 3.

### D05 · The linear perturbation equations (11.25)–(11.27) from the Boussinesq set and (11.22)–(11.24) — ★, 8 steps, in C03 (notebook)
- **Goal.** Get the equations obeyed by small disturbances of a fluid at rest heated from below.
- **Start.** $\nabla\cdot\tilde{\mathbf u}=0$, $\frac{\partial\tilde{\mathbf u}}{\partial t}+(\tilde{\mathbf u}\cdot\nabla)\tilde{\mathbf u}=
  -\frac1{\rho_0}\nabla\tilde p-g[1-\alpha(\tilde T-T_0)]\mathbf e_z+\nu\nabla^2\tilde{\mathbf u}$, $\frac{\partial\tilde T}{\partial t}+
  (\tilde{\mathbf u}\cdot\nabla)\tilde T=\kappa\nabla^2\tilde T$ (4.10, 4.86, 4.89) — *in words:* the Boussinesq equations of Ch. 4.
- **Plan.** (1) Insert base state + disturbance. (2) Subtract the base state's own equations. (3) Drop products of disturbances.
- **Tools.** The Boussinesq set (R12); the conduction state (R13); orders of smallness (P68).
- **Assumptions.** Boussinesq (ρ varies only in the buoyancy term); small disturbances (step 8).
- **Steps.**
  1. *did:* Split every field · *tex:* $\tilde{\mathbf u}=0+\mathbf u,\ \tilde T=\bar T(z)+T',\ \tilde p=P(z)+p$ (11.22) · *why:* the
     base state is at rest with conduction; u, T′, p are the disturbance. · *plain:* total = background + disturbance.
  2. *did:* Insert into momentum · *tex:* $\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=-\frac{\nabla P+\nabla p}
     {\rho_0}-g[1-\alpha(\bar T-T_0)]\mathbf e_z+g\alpha T'\mathbf e_z+\nu\nabla^2\mathbf u$ · *why:* substitute (11.22); the buoyancy
     bracket splits into its base part and $+g\alpha T'\mathbf e_z$. · *plain:* the disturbance feels the base forces plus its own
     buoyancy.
  3. *did:* Subtract the base balance · *tex:* $\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=-\frac1{\rho_0}\nabla p
     +g\alpha T'\mathbf e_z+\nu\nabla^2\mathbf u$ · *why:* (11.23) says $0=-\frac1{\rho_0}\nabla P-g[1-\alpha(\bar T-T_0)]\mathbf e_z$:
     the base pressure gradient cancels the base buoyancy exactly. · *plain:* only the extra buoyancy of a warm blob drives motion.
  4. *did:* Insert into the heat equation · *tex:* $\frac{\partial T'}{\partial t}+(\mathbf u\cdot\nabla)(\bar T+T')=\kappa\nabla^2\bar T+
     \kappa\nabla^2T'$ · *why:* T̄ does not depend on t, so ∂_tT̃ = ∂_tT′. · *plain:* the disturbance carries both temperatures.
  5. *did:* Use the conduction profile · *tex:* $\kappa\nabla^2\bar T=\kappa\frac{d^2\bar T}{dz^2}=0$ · *why:* (11.23)–(11.24): T̄ is
     linear in z, $\bar T=T_0-\tfrac12\Delta T-\Gamma z$. · *plain:* the background conducts heat steadily.
  6. *did:* Evaluate (u·∇)T̄ · *tex:* $(\mathbf u\cdot\nabla)\bar T=w\frac{d\bar T}{dz}=-w\Gamma$ · *why:* T̄ depends on z only, so only
     w∂_z survives; dT̄/dz = −Γ by (11.24). The sign comes from Γ = −dT̄/dz (slip #10's convention). · *plain:* rising fluid (w > 0)
     brings warmer fluid from below, which the −wΓ term counts as a local warming.
  7. *did:* Collect the nonlinear set · *tex:* $\frac{\partial T'}{\partial t}-w\Gamma+(\mathbf u\cdot\nabla)T'=\kappa\nabla^2T',\quad
     \nabla\cdot\mathbf u=0$ · *why:* steps 4–6 combined; continuity holds for u because the base flow is zero. · *plain:* the exact
     disturbance equations (p. 485).
  8. *did:* Drop products of disturbances · *tex:* $\nabla\cdot\mathbf u=0,\ \frac{\partial\mathbf u}{\partial t}=-\frac1{\rho_0}\nabla p+
     g\alpha T'\mathbf e_z+\nu\nabla^2\mathbf u,\ \frac{\partial T'}{\partial t}-w\Gamma=\kappa\nabla^2T'$ (11.25)–(11.27) · *why:*
     (u·∇)u and (u·∇)T′ are quadratic in the small disturbance (P68). · *plain:* linear equations — waves of different shapes
     will not interact.
- **Result.** (11.25)–(11.27) — *in words:* a warm blob is pushed up by gαT′ and its temperature is fed by −wΓ.
- **Check.** Units: each term of (11.26) m/s², of (11.27) K/s ✓. Sign: with w > 0 (rising) and Γ > 0, (11.27) gives ∂T′/∂t =
  +wΓ > 0: rising fluid is warmer than its surroundings — destabilising ✓ (`ch11.benard_perturbation_sympy`).
- **What it means.** The feedback loop of convection: w makes T′, T′ makes w. Fails for large ΔT (Boussinesq) and finite amplitude.
- **Traps.** The sign of the −wΓ term (Γ = −dT̄/dz, not dT̄/dz — slip #10); forgetting that the base pressure cancels the *whole*
  base buoyancy.


### D06 · Eliminating the pressure: (11.28) → $\frac{\partial}{\partial t}\nabla^2w=+g\alpha\nabla_H^2T'+\nu\nabla^4w$ (11.29) — ★★, 8 steps, in C03 (notebook)
- **Goal.** Get one equation for the vertical velocity w alone, without the unknown pressure.
- **Start.** The z-component of (11.26), $\frac{\partial w}{\partial t}=-\frac1{\rho_0}\frac{\partial p}{\partial z}+g\alpha T'+\nu\nabla^2w$,
  and $\nabla\cdot\mathbf u=0$ (11.25) — *in words:* vertical momentum and continuity of the disturbance.
- **Plan.** (1) Take the Laplacian of the w-equation. (2) Find what the pressure must satisfy from the divergence of the momentum
  equation. (3) Subtract to remove the pressure.
- **Tools.** Divergence and Laplacian (Ch. 2); constant-coefficient operators commute (Schwarz, P121); Ch. 10's pressure Poisson idea.
- **Assumptions.** Smooth fields (derivatives commute); incompressible disturbance (step 4).
- **Steps.**
  1. *did:* Take the Laplacian of the w-equation · *tex:* $\frac{\partial}{\partial t}\nabla^2w=-\frac1{\rho_0}\nabla^2\frac{\partial p}
     {\partial z}+g\alpha\nabla^2T'+\nu\nabla^4w$ (11.28) · *why:* ∇² commutes with ∂_t and ∂_z because the coefficients are
     constants (P121); we want a form in which the pressure appears as ∇²∂_zp. · *plain:* the same balance, one level of
     derivatives higher.
  2. *did:* Take the divergence of (11.26) · *tex:* $\frac{\partial}{\partial t}\nabla\cdot\mathbf u=-\frac1{\rho_0}\nabla^2p+g\alpha
     \frac{\partial T'}{\partial z}+\nu\nabla^2(\nabla\cdot\mathbf u)$ · *why:* ∇·∇p = ∇²p, ∇·(T′e_z) = ∂_zT′, and the divergence
     commutes with ∂_t and ∇². · *plain:* what the momentum equation says about the pressure.
  3. *did:* Use continuity · *tex:* $0=-\frac1{\rho_0}\nabla^2p+g\alpha\frac{\partial T'}{\partial z}$ · *why:* (11.25) ∇·u = 0 kills both
     the ∂_t and the ν terms — the book uses this without comment. · *plain:* the pressure is set by the buoyancy (a Poisson
     equation, as in Ch. 10's projection).
  4. *did:* Differentiate in z · *tex:* $0=-\frac1{\rho_0}\nabla^2\frac{\partial p}{\partial z}+g\alpha\frac{\partial^2T'}{\partial z^2}$
     · *why:* we need ∇²∂_zp to match step 1; ∂_z commutes with ∇² (P121). · *plain:* the pressure term of step 1, expressed by T′.
  5. *did:* Subtract step 4 from step 1 · *tex:* $\frac{\partial}{\partial t}\nabla^2w=g\alpha\Big(\nabla^2T'-\frac{\partial^2T'}{\partial
     z^2}\Big)+\nu\nabla^4w$ · *why:* the pressure terms are identical and cancel; that was the aim. · *plain:* no pressure left.
  6. *did:* Recognise the horizontal Laplacian · *tex:* $\nabla^2-\frac{\partial^2}{\partial z^2}=\frac{\partial^2}{\partial x^2}+
     \frac{\partial^2}{\partial y^2}\equiv\nabla_H^2$ · *why:* definition of ∇² in Cartesian coordinates. · *plain:* only horizontal
     variations of T′ drive w.
  7. *did:* Write the result · *tex:* $\frac{\partial}{\partial t}\nabla^2w=+g\alpha\nabla_H^2T'+\nu\nabla^4w$ (11.29) · *why:* steps 5–6.
     · *plain:* vertical motion is driven by horizontal temperature differences and damped by viscosity.
  8. *did:* Pair it with the heat equation · *tex:* $\frac{\partial T'}{\partial t}-w\Gamma=\kappa\nabla^2T'$ (11.27) · *why:* (11.27)
     and (11.29) are two equations for two unknowns, w and T′. · *plain:* the closed linear problem of convection.
- **Result.** (11.29) with (11.27) — *in words:* a uniformly warm layer does not move; side-to-side temperature differences do.
- **Check.** Units: every term of (11.29) 1/(m s²) ✓. A horizontally uniform T′ (∇_H²T′ = 0) gives no forcing of w — correct, it only
  changes the hydrostatic pressure. Planted variant ∇² instead of ∇_H² fails `benard_perturbation_sympy` ✓.
- **What it means.** Convection cells need horizontal structure: the K of every later step.
- **Traps.** Forgetting to differentiate the Poisson equation in z (leaves ∇²T′ instead of ∇_H²T′ — the planted mutant); dropping the
  ∂_t and ν terms of the divergence without citing ∇·u = 0.

### D07 · Non-dimensional equations (11.31)–(11.33), normal modes and the amplitude problem (11.34)–(11.38) — ★, 9 steps, in C03 (notebook)
- **Goal.** Turn (11.27) and (11.29) into two ODEs in z with a single parameter, the Rayleigh number.
- **Start.** (11.27), (11.29) and the walls $w=\partial w/\partial z=T'=0$ on $z=\pm d/2$ (11.30) — *in words:* the linear
  problem of D06 in dimensional form.
- **Plan.** (1) Scale lengths by d and time by d²/κ (keep w dimensional for now). (2) Insert normal modes in x, y, t. (3) Scale w to
  make Ra appear.
- **Tools.** Scaled variables and the chain rule (P133); operator substitution on $e^{\mathrm i(kx+ly)+\sigma t}$ (P177).
- **Assumptions.** k, l real (bounded in x, y) — step 5.
- **Steps.**
  1. *did:* Choose the scales · *tex:* $t\to\frac{d^2}{\kappa}t,\ (x,y,z)\to(xd,yd,zd)\ \Rightarrow\ \frac{\partial}{\partial t}\to\frac\kappa
     {d^2}\frac{\partial}{\partial t},\ \nabla\to\frac1d\nabla$ · *why:* d is the only length and d²/κ the time heat needs to diffuse
     across it (P133); the book reuses the same symbols for the new variables. · *plain:* measure in layer depths and diffusion times.
  2. *did:* Scale the heat equation · *tex:* $\Big(\frac{\partial}{\partial t}-\nabla^2\Big)T'=\frac{\Gamma d^2}{\kappa}w$ (11.31) · *why:*
     every term of (11.27) except wΓ gets κ/d²; divide by it. w stays dimensional on purpose. · *plain:* heating by vertical motion
     against diffusion.
  3. *did:* Scale the w-equation · *tex:* $\Big(\frac1\Pr\frac{\partial}{\partial t}-\nabla^2\Big)\nabla^2w=\frac{g\alpha d^2}\nu\nabla_H^2T'$
     (11.32) · *why:* in (11.29) the left side gets κ/d⁴, the viscous term ν/d⁴, the buoyancy term gα/d²; divide by ν/d⁴; Pr ≡ ν/κ. ·
     *plain:* buoyancy against viscosity, with the Prandtl number on the time derivative.
  4. *did:* Scale the walls · *tex:* $w=\partial w/\partial z=T'=0\ \text{on}\ z=\pm\tfrac12$ (11.33) · *why:* z = ±d/2 becomes ±½. ·
     *plain:* the same no-slip, fixed-temperature plates.
  5. *did:* Insert normal modes · *tex:* $w=\hat w(z)e^{\mathrm ikx+\mathrm ily+\sigma t},\ T'=\hat T(z)e^{\mathrm ikx+\mathrm ily+\sigma t}$ ·
     *why:* the coefficients do not depend on x, y, t (C01); k, l real keep the disturbance bounded. · *plain:* one horizontal wave at
     a time.
  6. *did:* Replace the operators · *tex:* $\frac{\partial}{\partial t}\to\sigma,\ \nabla_H^2\to-k^2-l^2\equiv-K^2,\ \nabla^2\to\frac{d^2}{dz^2}-K^2$
     · *why:* each derivative of the exponential multiplies it by the corresponding factor (P177). · *plain:* derivatives in x, y, t
     become numbers.
  7. *did:* Write the amplitude equations · *tex:* $\Big(\sigma+K^2-\frac{d^2}{dz^2}\Big)\hat T=\frac{\Gamma d^2}\kappa\hat w,\ \Big(\frac\sigma\Pr+
     K^2-\frac{d^2}{dz^2}\Big)\Big(\frac{d^2}{dz^2}-K^2\Big)\hat w=-\frac{g\alpha d^2K^2}\nu\hat T$ (11.34, 11.35) · *why:* steps 2–3 with
     step 6, divided by the common exponential. · *plain:* two linked ODEs in z.
  8. *did:* Rescale the velocity · *tex:* $W\equiv\frac{\Gamma d^2}\kappa\hat w\ \Rightarrow\ \Big(\frac\sigma\Pr+K^2-\frac{d^2}{dz^2}\Big)\Big(
     \frac{d^2}{dz^2}-K^2\Big)W=-\mathrm{Ra}K^2\hat T$ · *why:* multiply (11.35) by Γd²/κ: the right side becomes $-\frac{g\alpha\Gamma
     d^4}{\nu\kappa}K^2\hat T=-\mathrm{Ra}K^2\hat T$ by (11.21); the first equation becomes $(\sigma+K^2-\frac{d^2}{dz^2})\hat T=W$. ·
     *plain:* all the physics now sits in one number, Ra (11.36, 11.37).
  9. *did:* Write the conditions · *tex:* $W=\frac{dW}{dz}=\hat T=0\ \text{on}\ z=\pm\tfrac12$ (11.38) · *why:* (11.33) for the amplitudes.
     · *plain:* a sixth-order eigenvalue problem with eigenvalue σ.
- **Result.** $\big(\sigma+K^2-\frac{d^2}{dz^2}\big)\hat T=W$, $\big(\frac\sigma\Pr+K^2-\frac{d^2}{dz^2}\big)\big(\frac{d^2}{dz^2}-K^2\big)W=
  -\mathrm{Ra}K^2\hat T$ (11.36, 11.37) with (11.38) — *in words:* growth rate σ(K, Ra, Pr).
- **Check.** Dimensionless: every term of (11.36)–(11.37) is a pure number ✓. Pr appears only next to σ ✓ (so σ = 0 removes it).
- **What it means.** The whole of Bénard convection is one eigenvalue problem with three numbers; C03's code solves it.
- **Traps.** Expecting Ra before W is introduced (it appears only at step 8); forgetting that ∂_t brings κ/d², not 1/d.

### D08 · σ is real for Ra > 0 — exchange of stabilities (Exercise 11.6) — ★★, 12 steps, in C03 (notebook)
- **Goal.** Show that heated-from-below convection can only start as a stationary pattern: the growth rate σ has no imaginary
  part.
- **Start.** (11.36), (11.37) and (11.38), with complex σ allowed — *in words:* the amplitude problem of D07.
- **Plan.** (1) Multiply each equation by the conjugate of its unknown and integrate over the layer. (2) Integrate by parts to get
  positive integrals. (3) Combine the two relations so the mixed integral cancels. (4) Take the imaginary part.
- **Tools.** Integration by parts (P218a reminder); real and imaginary parts of a complex identity (P260); ∫∣f∣² > 0 (P260).
- **Assumptions.** Ra > 0 (heated from below), used in step 12; the rigid conditions (11.38) (free walls work too).
- **Steps.** (all integrals over −½ ≤ z ≤ ½; D ≡ d/dz)
  1. *did:* Multiply (11.36) by T̂* and integrate · *tex:* $\sigma\int\lvert\hat T\rvert^2dz+\int\hat T^*(K^2-D^2)\hat T\,dz=\int\hat T^*W\,dz$
     · *why:* multiplying by the conjugate and integrating turns functions into numbers whose signs we can read (P260). · *plain:*
     an 'energy' version of the heat equation.
  2. *did:* Integrate the D² term by parts · *tex:* $-\int\hat T^*\hat T''dz=\big[-\hat T^*\hat T'\big]_{-1/2}^{1/2}+\int\lvert\hat T'\rvert^2dz=
     \int\lvert\hat T'\rvert^2dz$ · *why:* P218a; the boundary term vanishes because T̂ = 0 at both walls (11.38). · *plain:* a
     curvature term becomes a positive gradient term.
  3. *did:* Name the positive integrals · *tex:* $\sigma I_1+I_2=\int\hat T^*W\,dz,\quad I_1=\int\lvert\hat T\rvert^2,\ I_2=\int(\lvert
     \hat T'\rvert^2+K^2\lvert\hat T\rvert^2)$ · *why:* steps 1–2; I₁, I₂ are real and positive for a nonzero T̂ (P260). · *plain:*
     Exercise 11.6's first relation.
  4. *did:* Rewrite (11.37) · *tex:* $\frac\sigma\Pr(D^2-K^2)W-(D^2-K^2)^2W=-\mathrm{Ra}K^2\hat T$ · *why:* expand
     $(\frac\sigma\Pr+K^2-D^2)=\frac\sigma\Pr-(D^2-K^2)$. · *plain:* the same equation, sorted by σ.
  5. *did:* Multiply by −W* and integrate · *tex:* $-\frac\sigma\Pr\int W^*(D^2-K^2)W+\int W^*(D^2-K^2)^2W=\mathrm{Ra}K^2\int W^*\hat T$ · *why:*
     the same idea as step 1, with a sign chosen so the integrals come out positive. · *plain:* an energy version of the w-equation.
  6. *did:* Integrate the first term by parts · *tex:* $-\int W^*(D^2-K^2)W\,dz=\int(\lvert W'\rvert^2+K^2\lvert W\rvert^2)dz\equiv J_1$ · *why:*
     as step 2, using W = 0 at the walls. · *plain:* a positive number.
  7. *did:* Move both D² across in the second · *tex:* $\int W^*(D^2-K^2)G\,dz=\int G\,(D^2-K^2)W^*dz,\quad G\equiv(D^2-K^2)W$ · *why:* two
     integrations by parts; the boundary terms contain W and W′, both zero at rigid walls (11.38). · *plain:* the operator can act on
     the other factor.
  8. *did:* Recognise a squared modulus · *tex:* $\frac\sigma\Pr J_1+J_2=\mathrm{Ra}K^2\int W^*\hat T\,dz,\quad J_2=\int\lvert(D^2-K^2)W\rvert^2dz$
     · *why:* $(D^2-K^2)W^*$ is the conjugate of G, so G times it is ∣G∣² ≥ 0. · *plain:* Exercise 11.6's second relation.
  9. *did:* Conjugate the first relation · *tex:* $\sigma^*I_1+I_2=\int\hat T\,W^*dz$ · *why:* I₁, I₂ are real; conjugating swaps the
     stars on the right — now it matches the right side of step 8. · *plain:* the same mixed integral appears in both.
  10. *did:* Eliminate the mixed integral · *tex:* $\frac\sigma\Pr J_1+J_2=\mathrm{Ra}K^2(\sigma^*I_1+I_2)$ · *why:* substitute step 9 into
      step 8. · *plain:* one equation with σ and σ*.
  11. *did:* Take the imaginary part · *tex:* $\sigma_i\Big(\frac{J_1}\Pr+\mathrm{Ra}K^2I_1\Big)=0$ · *why:* Im σ = σ_i, Im σ* = −σ_i;
      J₁, J₂, I₁, I₂ real (P260). · *plain:* σ_i times a bracket is zero.
  12. *did:* Use Ra > 0 · *tex:* $\sigma_i=0$ · *why:* for Ra > 0 the bracket is a sum of positive terms, never zero; so σ_i must be
      zero. · *plain:* σ is real: the onset (σ = 0) is a stationary pattern — exchange of stabilities.
- **Result.** σ real for every Bénard mode when Ra > 0 — *in words:* convection cells appear in place, they do not oscillate; the
  margin is σ = 0.
- **Check.** Every computed σ in C03 has ∣σ_i∣ < 1e-10 ✓. For Ra < 0 (heated from above) the argument fails, and indeed the free–free
  quadratic (D12) can give complex roots there. Optional `check_src`: `ch11.exchange_of_stabilities_sympy()["imag_identity"] == 0`.
- **What it means.** At onset we may set σ = 0 (C04). With a second diffusing component the bracket can vanish and oscillatory
  onset appears (C06's diffusive regime).
- **Traps.** Forgetting to conjugate one relation before combining (the mixed integrals then do not match); needing both W = 0
  and W′ = 0 for the double integration by parts (free walls use W = W″ = 0 instead).

### D09 · The marginal pair (11.39) → the sixth-order equation (11.40) with the rigid conditions (11.41) — ★, 6 steps, in C04 (notebook · `benard_neutral_curve`)
- **Goal.** At the onset (σ = 0) reduce the two amplitude equations to one equation for W alone.
- **Start.** (11.36), (11.37), (11.38) and σ = 0 (D08) — *in words:* the marginal state of convection.
- **Plan.** (1) Set σ = 0. (2) Apply one operator to eliminate T̂. (3) Translate the conditions.
- **Tools.** Operator elimination (P177).
- **Assumptions.** σ = 0 at the margin (D08).
- **Steps.**
  1. *did:* Set σ = 0 in (11.36) · *tex:* $\Big(\frac{d^2}{dz^2}-K^2\Big)\hat T=-W$ · *why:* the onset is σ = 0 by D08; multiply by −1 to
     write the operator as (D² − K²). · *plain:* the temperature is driven by the vertical velocity.
  2. *did:* Set σ = 0 in (11.37) · *tex:* $\Big(\frac{d^2}{dz^2}-K^2\Big)^2W=\mathrm{Ra}K^2\hat T$ (11.39) · *why:* $(K^2-D^2)(D^2-K^2)=
     -(D^2-K^2)^2$; multiply by −1. Pr has disappeared. · *plain:* the velocity is driven by the temperature.
  3. *did:* Apply (D² − K²) to the second · *tex:* $\Big(\frac{d^2}{dz^2}-K^2\Big)^3W=\mathrm{Ra}K^2\Big(\frac{d^2}{dz^2}-K^2\Big)\hat T$ · *why:*
     constant-coefficient operators can be applied to both sides (P177); this produces the combination step 1 knows. · *plain:* one
     more derivative level.
  4. *did:* Use step 1 · *tex:* $\Big(\frac{d^2}{dz^2}-K^2\Big)^3W=-\mathrm{Ra}K^2W$ (11.40) · *why:* replace (D² − K²)T̂ by −W. · *plain:*
     one sixth-order equation for W.
  5. *did:* Keep the velocity conditions · *tex:* $W=\frac{dW}{dz}=0\ \text{on}\ z=\pm\tfrac12$ · *why:* directly from (11.38). · *plain:*
     no slip, no flow through the plates.
  6. *did:* Translate T̂ = 0 · *tex:* $\Big(\frac{d^2}{dz^2}-K^2\Big)^2W=0\ \text{on}\ z=\pm\tfrac12$ (11.41) · *why:* by step 2,
     (D² − K²)²W = RaK²T̂, and T̂ = 0 at the walls. · *plain:* the third condition, written with W only.
- **Result.** $\big(\frac{d^2}{dz^2}-K^2\big)^3W=-\mathrm{Ra}K^2W$ (11.40) with $W=\frac{dW}{dz}=\big(\frac{d^2}{dz^2}-K^2\big)^2W=0$ on
  $z=\pm\tfrac12$ (11.41) — *in words:* six conditions for a sixth-order equation; only special Ra(K) allow a nonzero W.
- **Check.** No Pr ✓; the free–free sine mode satisfies (11.40) with Ra = (π² + K²)³/K² (D11) ✓.
- **What it means.** Onset depends only on Ra and K: a neutral curve Ra(K).
- **Traps.** Getting the third condition from the wrong equation; keeping Pr.

### D10 · Characteristic roots (11.42), the even solution, the 3 × 3 determinant and $\mathrm{Ra}_c=1707.76$ at $K_c=3.117$ — ★★★, 14 steps, in C04 (notebook · `benard_neutral_curve`)
- **Goal.** Solve (11.40)–(11.41) exactly for rigid walls and find the smallest Rayleigh number at which convection can start.
- **Start.** $\big(\frac{d^2}{dz^2}-K^2\big)^3W=-\mathrm{Ra}K^2W$ (11.40) with $W=\frac{dW}{dz}=\big(\frac{d^2}{dz^2}-K^2\big)^2W=0$ on
  $z=\pm\tfrac12$ (11.41) — *in words:* the marginal problem of D09.
- **Plan.** (1) Exponential trial → a cubic for q² − K². (2) Its three cube roots → six exponents. (3) Even solution and its three wall
  conditions → a 3 × 3 determinant. (4) Find its zero in Ra for each K, then minimise over K.
- **Tools.** Exponential trial (P44); cube roots of a negative number (P262); even and odd functions (P261); determinants and
  nonzero solutions (P53, P58); `brentq` (P108); `minimize_scalar` (P170); neutral-curve recipe (P263); sympy (P40).
- **Assumptions.** Ra > K⁴ (q₀ real, step 5); the even mode (step 7; the odd mode is N123).
- **Steps.**
  1. *did:* Try W = e^{qz} · *tex:* $(q^2-K^2)^3=-\mathrm{Ra}K^2$ · *why:* constant coefficients: each d²/dz² becomes q² (P44); divide by
     e^{qz} ≠ 0. · *plain:* an algebraic equation for the exponent q.
  2. *did:* Call S = q² − K² · *tex:* $S^3=-\mathrm{Ra}K^2=\mathrm{Ra}K^2e^{\mathrm i\pi}$ · *why:* a cube equation for S; writing −1 as
     e^{iπ} prepares the three cube roots (P262). · *plain:* we need the cube roots of a negative number.
  3. *did:* Take the three cube roots · *tex:* $S=-(\mathrm{Ra}K^2)^{1/3},\quad S=\tfrac12(\mathrm{Ra}K^2)^{1/3}(1\pm\mathrm i\sqrt3)$ · *why:* the
     cube roots of −1 are −1 and e^{±iπ/3} = ½(1 ± i√3) (P262). · *plain:* one real and two complex conjugate values.
  4. *did:* Return to q² · *tex:* $q^2=-K^2\Big[\Big(\frac{\mathrm{Ra}}{K^4}\Big)^{1/3}-1\Big],\quad q^2=K^2\Big[1+\frac12\Big(\frac{\mathrm{Ra}}
     {K^4}\Big)^{1/3}(1\pm\mathrm i\sqrt3)\Big]$ (11.42) · *why:* q² = K² + S and $(\mathrm{Ra}K^2)^{1/3}=K^2(\mathrm{Ra}/K^4)^{1/3}$. ·
     *plain:* the book's three roots.
  5. *did:* Name the first pair · *tex:* $q=\pm\mathrm iq_0,\quad q_0=K\Big[\Big(\frac{\mathrm{Ra}}{K^4}\Big)^{1/3}-1\Big]^{1/2}$ · *why:* for
     Ra > K⁴ the first q² is negative, so its square roots are imaginary (P45). · *plain:* an oscillating part, cos q₀z.
  6. *did:* Name the other two pairs · *tex:* $\pm q,\ \pm q^*$ · *why:* the last two q² are complex conjugates, so their square roots
     are q and q* (principal branch, Re q > 0) with their negatives. · *plain:* six exponents in all.
  7. *did:* Keep the even combination · *tex:* $W=A\cos q_0z+B\cosh qz+C\cosh q^*z$ · *why:* both walls are alike, so even and odd
     solutions separate (P261); e^{±iq₀z} pair into cos q₀z, e^{±qz} into cosh qz. W is real when A is real and C = B*. · *plain:*
     one row of cells.
  8. *did:* Differentiate once · *tex:* $\frac{dW}{dz}=-Aq_0\sin q_0z+Bq\sinh qz+Cq^*\sinh q^*z$ · *why:* derivatives of cos and cosh
     (complex arguments allowed). · *plain:* for the no-slip condition.
  9. *did:* Apply (D² − K²) twice · *tex:* $\Big(\frac{d^2}{dz^2}-K^2\Big)^2W=A(q_0^2+K^2)^2\cos q_0z+B(q^2-K^2)^2\cosh qz+C(q^{*2}-K^2)^2\cosh q^*z$
     · *why:* (D² − K²)cos q₀z = −(q₀² + K²)cos q₀z and (D² − K²)cosh qz = (q² − K²)cosh qz; squaring removes the sign. ⚠️ slip #1:
     the book prints B in the third term; it must be C. · *plain:* for the temperature condition.
  10. *did:* Use the symmetry · *tex:* $\text{conditions at }z=+\tfrac12\ \Leftrightarrow\ \text{conditions at }z=-\tfrac12$ · *why:* W,
      d²W/dz² … are even and dW/dz is odd, so each condition at −½ repeats the one at +½. · *plain:* three conditions remain.
  11. *did:* Write them as a matrix · *tex:* $\begin{bmatrix}\cos\frac{q_0}2&\cosh\frac q2&\cosh\frac{q^*}2\\-q_0\sin\frac{q_0}2&q\sinh\frac q2&
      q^*\sinh\frac{q^*}2\\(q_0^2+K^2)^2\cos\frac{q_0}2&(q^2-K^2)^2\cosh\frac q2&(q^{*2}-K^2)^2\cosh\frac{q^*}2\end{bmatrix}\begin{bmatrix}A\\B\\C
      \end{bmatrix}=0$ · *why:* steps 7–9 at z = ½. A nonzero (A, B, C) exists only if the determinant vanishes (P53, P58). ·
      *plain:* the eigenvalue condition.
  12. *did:* See that det is imaginary · *tex:* $\overline{\det M}=-\det M\ \Rightarrow\ \det M=\mathrm i\,\mathrm{Im}\det M$ · *why:* column 3
      is the conjugate of column 2 and column 1 is real; conjugating M equals swapping two columns, which flips the sign of a
      determinant. So solve Im det M = 0. · *plain:* one real equation in Ra for each K.
  13. *did:* Find the root in Ra at fixed K · *tex:* $\mathrm{Im}\det M(\mathrm{Ra};K)=0\ \Rightarrow\ \mathrm{Ra}(K)$, e.g. Ra(2) = 2177.41 ·
      *why:* scan Ra upward from just above K⁴ (below it q₀ is not real and spurious sign changes appear), then `brentq` on the
      first sign change (P108). · *plain:* one point of the neutral curve.
  14. *did:* Minimise over K · *tex:* $\mathrm{Ra}_c=\min_K\mathrm{Ra}(K)=1707.76\ \text{at}\ K_c=3.117$ · *why:* every K is allowed, so the
      layer first goes unstable at the lowest point of Ra(K) (P263, `minimize_scalar`). · *plain:* the critical Rayleigh number.
- **Result.** $\mathrm{Ra}_c=1707.76$, $K_c=3.117$ (Chandrasekhar 1961; the book rounds to 1708 and 3.12); cell width
  $\lambda_c=2\pi d/K_c\approx2.02d$ — *in words:* rigid plates start convecting at Ra ≈ 1708 in cells about twice as wide as the
  layer is deep (a pair of counter-rotating rolls).
- **Check.** The determinant root equals the Chebyshev eigenvalue at K = 2, 3.1163, 5 to 1e-8 (two independent routes) ✓; with the
  printed B (slip #1) the determinant is no longer imaginary and the root moves ✓. Dimensionless ✓.
- **sympy check intent (`check_src`).** Verify each move, not just the answer:
  ```python
  import sympy as sp                                   # symbolic algebra
  z, K, Ra = sp.symbols("z K Ra", positive=True)       # height and the two parameters
  s = (Ra/K**4)**sp.Rational(1, 3)                     # the cube-root factor of (11.42)
  for q2 in (-K**2*(s - 1), K**2*(1 + s*(1 + sp.sqrt(3)*sp.I)/2), K**2*(1 + s*(1 - sp.sqrt(3)*sp.I)/2)):
      print(sp.simplify(sp.expand((q2 - K**2)**3 + Ra*K**2)))   # step 1: each root of (11.42) satisfies (q^2-K^2)^3 = -Ra K^2 -> 0
  q = sp.symbols("q")                                  # a generic exponent
  L = lambda f: sp.diff(f, z, 2) - K**2*f              # the operator d^2/dz^2 - K^2
  print(sp.simplify(L(L(sp.cosh(q*z))) - (q**2 - K**2)**2*sp.cosh(q*z)))   # step 9: (D^2-K^2)^2 cosh qz -> 0
  import numpy as np                                   # numbers for steps 12-14
  from fluidpy import ch11_instability as ch11         # the tested functions
  d = ch11.benard_determinant(1707.762, 3.1163)        # step 12: the determinant at the critical point
  print(abs(d.real) < 1e-9*abs(d))                     # True: purely imaginary
  print(np.isclose(ch11.benard_marginal_Ra_det(3.1163), ch11.benard_marginal_Ra(3.1163), rtol=1e-8))   # step 13: two routes agree
  ```
- **What it means.** The heating must beat a threshold set jointly by viscosity and diffusion; the preferred cell is the one that
  minimises the cost. Fails above onset (nonlinear rolls, N52) and for non-Boussinesq or very deep layers.
- **Traps.** Using q twice instead of q and q* (W would not be real); the printed B (slip #1); trusting a sign change of Im det below
  Ra = K⁴; forgetting that the minimum over K, not a single K, gives Ra_c.

### D11 · Stress-free walls (11.43), sine modes, (11.44) and $\frac{d\mathrm{Ra}}{dK^2}=0\Rightarrow K_c^2=\pi^2/2$, $\mathrm{Ra}_c=\tfrac{27}{4}\pi^4$ — ★★, 11 steps, in C05 (notebook · `benard_neutral_curve`)
- **Goal.** Solve the marginal problem by hand when neither wall can hold a shear stress.
- **Start.** (11.40) and the stress-free, isothermal walls: $w=T'=\mu\big(\frac{\partial u}{\partial z}+\frac{\partial w}{\partial x}\big)=
  \mu\big(\frac{\partial v}{\partial z}+\frac{\partial w}{\partial y}\big)=0$ (p. 489) — *in words:* the walls stop vertical flow but
  let the fluid slide.
- **Plan.** (1) Turn zero stress into a condition on W. (2) Show all even derivatives vanish. (3) Use sine modes. (4) Minimise over K².
- **Tools.** Differentiating continuity; Fourier sine modes and orthogonality (P151); the quotient rule (P265); sympy (P40).
- **Assumptions.** Stress-free, isothermal walls; the gravest mode n = 1 (step 9).
- **Steps.**
  1. *did:* Use w = 0 along the wall · *tex:* $w=0\ \text{for all}\ x,y\ \Rightarrow\ \frac{\partial w}{\partial x}=\frac{\partial w}{\partial y}
     =0$ · *why:* a function that is zero everywhere on the plane has zero derivatives along it. · *plain:* the wall-normal velocity
     does not vary along the wall.
  2. *did:* Simplify zero stress · *tex:* $\frac{\partial u}{\partial z}=\frac{\partial v}{\partial z}=0$ · *why:* step 1 removes the ∂w terms
     from the two stress conditions. · *plain:* no shear at the wall.
  3. *did:* Differentiate continuity in z · *tex:* $\frac{\partial^2w}{\partial z^2}=-\frac{\partial}{\partial x}\frac{\partial u}{\partial z}-
     \frac{\partial}{\partial y}\frac{\partial v}{\partial z}=0$ · *why:* ∂_z of ∇·u = 0, then step 2 (derivatives along the wall of zero
     are zero). · *plain:* W″ = 0 at a free wall.
  4. *did:* Collect the conditions · *tex:* $W=\Big(\frac{d^2}{dz^2}-K^2\Big)^2W=\frac{d^2W}{dz^2}=0\ \text{on}\ z=\pm\tfrac12$ · *why:* W = 0,
     step 3, and T̂ = 0 translated as in D09 step 6. · *plain:* three conditions per wall.
  5. *did:* Expand the operator · *tex:* $\Big(\frac{d^2}{dz^2}-K^2\Big)^2W=W''''-2K^2W''+K^4W=W''''$ · *why:* W = W″ = 0 there. So
     $W=W''=W''''=0$ (11.43) — the book's "expanding the products of operators". · *plain:* the fourth derivative vanishes too.
  6. *did:* Show every even derivative vanishes · *tex:* $W^{(6)}=-\mathrm{Ra}K^2W+\dots=0,\ \text{and so on}$ · *why:* (11.40) expresses
     W⁽⁶⁾ through W, W″, W⁗, all zero at the wall; differentiating (11.40) twice repeats the argument for W⁽⁸⁾, … . · *plain:* only
     odd derivatives survive at the walls.
  7. *did:* Choose sine modes · *tex:* $W=A\sin n\pi\big(z+\tfrac12\big)$ · *why:* sin nπ(z + ½) and all its even derivatives vanish at
     z = ±½. ⚠️ slip #2: the book prints sin(nπz), which vanishes at z = ±½ only for even n; n = 1 here is cos πz. · *plain:* n
     half-waves across the layer.
  8. *did:* Insert into (11.40) · *tex:* $-(n^2\pi^2+K^2)^3=-\mathrm{Ra}K^2\ \Rightarrow\ \mathrm{Ra}=\frac{(n^2\pi^2+K^2)^3}{K^2}$ (11.44) · *why:*
     d²/dz² → −n²π², so (D² − K²) → −(n²π² + K²); cube it; divide by −K². · *plain:* the neutral curve as a formula.
  9. *did:* Differentiate in x = K² (n = 1) · *tex:* $\frac{d\mathrm{Ra}}{dK^2}=\frac{3(\pi^2+K^2)^2}{K^2}-\frac{(\pi^2+K^2)^3}{K^4}$ · *why:*
     quotient rule in the variable K² (P265); the minimum over K is the minimum over K². ⚠️ slip #3: the book prints a 3 on the second
     term; with it there is no root. · *plain:* the slope of the neutral curve.
  10. *did:* Set the slope to zero · *tex:* $3K^2-(\pi^2+K^2)=0\ \Rightarrow\ K_c^2=\frac{\pi^2}2$ · *why:* multiply by
      K⁴/(π² + K²)² > 0. · *plain:* the cheapest cell, K_c = π/√2 ≈ 2.22.
  11. *did:* Evaluate Ra there · *tex:* $\mathrm{Ra}_c=\frac{(3\pi^2/2)^3}{\pi^2/2}=\frac{27\pi^6/8}{\pi^2/2}=\frac{27}4\pi^4\approx657.5$ · *why:*
      substitute K² = π²/2 into (11.44) with n = 1; higher n only raise Ra. · *plain:* the free–free threshold.
- **Result.** $\mathrm{Ra}=(n^2\pi^2+K^2)^3/K^2$ (11.44), $K_c^2=\pi^2/2$, $\mathrm{Ra}_c=\tfrac{27}4\pi^4=657.51$ — *in words:* without
  wall friction the threshold drops from 1708 to 657.5 and the cells widen (2π/K_c = 2.83 d).
- **Check.** Units: dimensionless ✓. Numbers: Ra(2) = 667.01 > 657.51 ✓. The Chebyshev solver with free walls gives 657.511 at
  2.2214 ✓. Optional `check_src`: `ch11.benard_free_free_sympy()` (residual of (11.40) is 0; root π²/2; printed root []).
- **What it means.** The formula used by double diffusion (C06: Rs − Ra = 27π⁴/4) and Lorenz (C15: b = 8/3 from K² = π²/2).
- **Traps.** sin nπz (slip #2); the extra 3 (slip #3); forgetting that higher n give larger Ra; cubing −(π² + K²) and losing the sign.

### D12 · The free–free growth rate $(\sigma+a^2)(\sigma/\Pr+a^2)a^2=\mathrm{Ra}K^2$, $a^2=\pi^2+K^2$ (ours; the book never writes it) — ★★, 8 steps, in C05 (notebook · `normal_mode_growth`)
- **Goal.** Find how fast a free–free convection mode grows, not only where it is marginal.
- **Start.** (11.36), (11.37) with stress-free isothermal walls — *in words:* the full amplitude problem with σ kept.
- **Plan.** (1) Try the same sine for W and T̂. (2) Turn operators into numbers. (3) Eliminate the amplitudes. (4) Solve the quadratic.
- **Tools.** Sine modes from D11; operator substitution (P177); quadratic formula (P159).
- **Assumptions.** Stress-free isothermal walls; n = 1.
- **Steps.**
  1. *did:* Try one sine for both · *tex:* $W=W_0\sin\pi\big(z+\tfrac12\big),\ \hat T=T_0\sin\pi\big(z+\tfrac12\big)$ · *why:* this shape meets
     every free–free condition (D11 steps 4–7). · *plain:* the gravest mode.
  2. *did:* Turn the operators into numbers · *tex:* $\frac{d^2}{dz^2}\to-\pi^2,\ \ a^2\equiv\pi^2+K^2$ · *why:* d²/dz² of the sine is −π²
     times it; so (K² − d²/dz²) → a² and (d²/dz² − K²) → −a². · *plain:* derivatives become multiplications.
  3. *did:* Use (11.36) · *tex:* $(\sigma+a^2)T_0=W_0$ · *why:* $(\sigma+K^2-\frac{d^2}{dz^2})\to\sigma+a^2$. · *plain:* the temperature
     amplitude follows the velocity.
  4. *did:* Use (11.37) · *tex:* $\Big(\frac\sigma\Pr+a^2\Big)a^2W_0=\mathrm{Ra}K^2T_0$ · *why:* $(\frac\sigma\Pr+a^2)(-a^2)W_0=-\mathrm{Ra}K^2T_0$;
     multiply by −1. · *plain:* buoyancy drives the velocity.
  5. *did:* Eliminate W₀ · *tex:* $(\sigma+a^2)\Big(\frac\sigma\Pr+a^2\Big)a^2=\mathrm{Ra}K^2$ · *why:* insert step 3 into step 4 and divide by
     T₀ ≠ 0. · *plain:* the growth rate obeys a quadratic.
  6. *did:* Expand · *tex:* $\frac{\sigma^2}\Pr+a^2\Big(1+\frac1\Pr\Big)\sigma+a^4-\frac{\mathrm{Ra}K^2}{a^2}=0$ · *why:* multiply out and divide
     by a². · *plain:* a standard quadratic in σ.
  7. *did:* Check the discriminant · *tex:* $a^4\Big(1-\frac1\Pr\Big)^2+\frac{4\mathrm{Ra}K^2}{\Pr\,a^2}>0$ · *why:* $(1+\frac1\Pr)^2-\frac4\Pr=
     (1-\frac1\Pr)^2$; the rest is positive for Ra > 0. · *plain:* both roots are real — exchange of stabilities again (D08).
  8. *did:* Find where σ = 0 · *tex:* $\sigma=0\ \Leftrightarrow\ a^6=\mathrm{Ra}K^2\ \Leftrightarrow\ \mathrm{Ra}=\frac{(\pi^2+K^2)^3}{K^2}$ ·
     *why:* set σ = 0 in step 5; it is (11.44) with n = 1; the larger root is positive exactly when Ra exceeds it. · *plain:* the
     growth curve crosses zero on the neutral curve.
- **Result.** $\sigma_\pm=\frac{\Pr}2\Big[-a^2\big(1+\tfrac1\Pr\big)\pm\sqrt{a^4\big(1-\tfrac1\Pr\big)^2+\frac{4\mathrm{Ra}K^2}{\Pr a^2}}\Big]$ —
  *in words:* one mode grows when Ra beats (11.44), the other always decays.
- **Check.** K = π/√2, Ra = 2000, Pr = 1: σ = 11.0155, −40.6243 (agrees with the Chebyshev solver to 1e-8) ✓; Ra = 27π⁴/4: σ₊ = 0 ✓.
- **What it means.** The σ(K) curve of E1's Bénard mode; Pr changes how fast, not whether.
- **Traps.** (d²/dz² − K²) on the sine is −a², not +a²; forgetting T₀ ≠ 0 when dividing.

### D13 · The double-diffusive marginal equations (11.45) and the salt-finger criterion (11.46) — ★★, 12 steps, in C06 (notebook · `salt_fingers`)
- **Goal.** Repeat the Bénard derivation with salt added (the book only says "repeat the derivation") and find when a column that
  is lighter on top overturns.
- **Start.** Boussinesq with $\tilde\rho=\rho_0[1-\alpha(\tilde T-T_0)+\beta(\tilde s-s_0)]$ (p. 492) and the salt equation
  $\frac{\partial\tilde s}{\partial t}+(\tilde{\mathbf u}\cdot\nabla)\tilde s=\kappa_s\nabla^2\tilde s$ — *in words:* two diffusing
  ingredients of density.
- **Plan.** (1) Linearise about linear profiles of T̄ and S̄ (as D05). (2) Eliminate pressure (as D06), scale (as D07) and set σ = 0.
  (3) Scale the amplitudes to the book's (11.45). (4) Use uniqueness to tie T̂ to ŝ, then reuse C05.
- **Tools.** D05–D07 moves; uniqueness of a linear boundary-value problem (P266); D11's free–free result.
- **Assumptions.** Stress-free walls at fixed T and S (step 11); σ = 0 at the margin — valid for fingers, not for the oscillatory
  diffusive regime (step 6); dT̄/dz, dS/dz ≠ 0 (step 7).
- **Steps.**
  1. *did:* Linearise about linear profiles · *tex:* $\frac{\partial T'}{\partial t}+w\frac{d\bar T}{dz}=\kappa\nabla^2T',\ \frac{\partial s'}{\partial t}+
     w\frac{dS}{dz}=\kappa_s\nabla^2s'$ · *why:* as D05 steps 4–8, with conduction and diffusion profiles T̄(z), S̄(z) of constant
     gradients. · *plain:* vertical motion carries both ingredients.
  2. *did:* Write the buoyancy · *tex:* $-\frac{g\rho'}{\rho_0}=g(\alpha T'-\beta s')$ · *why:* the linear EOS: warm is light (α), salty
     is heavy (β). · *plain:* heat lifts, salt sinks.
  3. *did:* Eliminate the pressure · *tex:* $\frac{\partial}{\partial t}\nabla^2w=g\nabla_H^2(\alpha T'-\beta s')+\nu\nabla^4w$ · *why:* the
     D06 moves with the new buoyancy. · *plain:* (11.29) with two sources.
  4. *did:* Scale (d, d²/κ, w by κ/d) · *tex:* $(\partial_t-\nabla^2)T'=-d\frac{d\bar T}{dz}W,\ \ (\partial_t-\tfrac{\kappa_s}\kappa\nabla^2)s'=-d\frac{dS}{dz}W$
     · *why:* as D07 step 2, now also making w dimensionless, w = (κ/d)W; heat diffusion sets the time unit. · *plain:* salt
     diffuses κ_s/κ times slower in these units.
  5. *did:* Note the sign of §11.5's Ra · *tex:* $\mathrm{Ra}\equiv\frac{g\alpha d^4(d\bar T/dz)}{\nu\kappa},\ \ \mathrm{Rs}'\equiv\frac{g\beta d^4(dS/dz)}{\nu\kappa}$
     · *why:* §11.5 defines Ra with dT̄/dz itself — negative when heated from below, the opposite of (11.21) (slip #10's second
     convention). · *plain:* here Ra > 0 means warm on top (stabilising).
  6. *did:* Insert normal modes, σ = 0 · *tex:* $(D^2-K^2)\hat T=d\frac{d\bar T}{dz}W,\ \tfrac{\kappa_s}\kappa(D^2-K^2)\hat s=d\frac{dS}{dz}W,\ (D^2-K^2)^2W=
     \frac{gd^3K^2}{\nu\kappa}(\alpha\hat T-\beta\hat s)$ · *why:* D07 steps 5–6 at the margin. σ = 0 is right for fingers; the diffusive
     onset is oscillatory (N54). · *plain:* three linked marginal equations.
  7. *did:* Scale the amplitudes · *tex:* $\theta\equiv-\frac{\hat T}{d\,d\bar T/dz},\ \ \varsigma\equiv-\frac{\hat s}{d\,dS/dz}$ · *why:* this makes
     both right-hand sides −W, the form the book prints (it reuses the letters T̂, ŝ for θ, ς). · *plain:* temperature and salt
     measured in units of their background change over d.
  8. *did:* Rewrite the three equations · *tex:* $(D^2-K^2)\theta=-W,\ \tfrac{\kappa_s}\kappa(D^2-K^2)\varsigma=-W,\ (D^2-K^2)^2W=-\mathrm{Ra}K^2\theta+\mathrm{Rs}'K^2\varsigma$
     (11.45) · *why:* substitute step 7 into step 6; the factors combine into Ra and Rs′ of step 5. · *plain:* the book's (11.45).
  9. *did:* Compare θ and (κ_s/κ)ς · *tex:* $(D^2-K^2)\big[\theta-\tfrac{\kappa_s}\kappa\varsigma\big]=0,\ \ \theta-\tfrac{\kappa_s}\kappa\varsigma=0\
     \text{at the walls}$ · *why:* subtract the first two equations; both amplitudes vanish at the walls. · *plain:* their difference
     solves a problem with zero data.
  10. *did:* Use uniqueness · *tex:* $\theta=\frac{\kappa_s}\kappa\varsigma\ \text{everywhere}$ · *why:* (D² − K²)f = 0 with f = 0 at both walls
      has only f = 0 (P266). The book writes T̂ = κ_sŝ/κ. · *plain:* temperature and salt amplitudes are locked together.
  11. *did:* Substitute into the third equation · *tex:* $(D^2-K^2)^2W=(\mathrm{Rs}-\mathrm{Ra})K^2\theta,\ \ \mathrm{Rs}\equiv\frac\kappa{\kappa_s}\mathrm{Rs}'=
      \frac{g\beta d^4(dS/dz)}{\nu\kappa_s}$ · *why:* ς = (κ/κ_s)θ. With (D² − K²)θ = −W this is exactly (11.39) with Ra replaced by
      Rs − Ra. · *plain:* double diffusion is Bénard with an effective Rayleigh number.
  12. *did:* Use the free–free threshold · *tex:* $\mathrm{Rs}-\mathrm{Ra}=\tfrac{27}4\pi^4\ \Leftrightarrow\ \frac{gd^4}\nu\Big[\frac\beta{\kappa_s}\frac{dS}{dz}-
      \frac\alpha\kappa\frac{d\bar T}{dz}\Big]=657$ (11.46) · *why:* identical equations and conditions give identical eigenvalues: D11's
      27π⁴/4 ≈ 657.5 (the book rounds to 657). · *plain:* the salt-finger criterion.
- **Result.** (11.45) and (11.46) — *in words:* the salt gradient, divided by the slow salt diffusivity, must beat the heat gradient,
  divided by the fast one, by 27π⁴/4.
- **Check.** κ_s = κ and dS/dz = 0: Rs = 0 and the criterion becomes −Ra = 657.5, i.e. (11.21)'s Ra = 657.5 ✓ (single component).
  Thermocline numbers (C06): lhs = 61 254 ✓. Dimensionless ✓.
- **What it means.** Because κ_s ≪ κ, a column can be lighter on top (αT_z > βS_z) and still satisfy (11.46): diffusion destabilises
  (C06's wedge). Fails for the diffusive regime (oscillatory onset), and for rigid walls (then 1708 replaces 657.5).
- **Traps.** Mixing the two Ra signs ((11.21) vs §11.5); forgetting that the book's T̂, ŝ in (11.45) are scaled amplitudes (verifier
  note i); using Rs′ (with κ) instead of Rs (with κ_s) in the final criterion.

### D14 · The linearised Taylor–Couette equations (11.50) and the narrow-gap system (11.51) from (11.47)–(11.49) — ★★★, 15 steps, in C07 (notebook · `taylor_couette_onset`)
- **Goal.** Derive the equations Taylor solved — the book only says "see Chandrasekhar (1961) for details".
- **Start.** The axisymmetric equations (11.47) (with continuity corrected, slip #12), the base flow $U_R=U_z=0$, $U_\varphi=AR+B/R$,
  $\frac1\rho\frac{dP}{dR}=\frac{U_\varphi^2}R$ (11.49) — *in words:* circular Couette flow and the Navier–Stokes equations in (R, φ, z).
- **Plan.** (1) Linearise about Couette flow. (2) Insert normal modes in z and t and eliminate û_z and p̂. (3) Drop curvature in a
  narrow gap. (4) Scale and rescale so one number, Ta, is left.
- **Tools.** Cylindrical Laplacian (P186); linearisation (P68); operator elimination (P177); narrow-gap approximation (P267); sympy.
- **Assumptions.** Axisymmetric disturbances (∂/∂φ = 0); narrow gap d ≪ R₁ (steps 11–12).
- **Steps.**
  1. *did:* Split the fields · *tex:* $\tilde{\mathbf u}=\mathbf U+\mathbf u,\ \tilde p=P+p$ (11.48), $\mathbf U=(0,U_\varphi(R),0)$ · *why:*
     base flow plus small disturbance, as in every block. · *plain:* swirl plus a small ring disturbance.
  2. *did:* Linearise the centrifugal term · *tex:* $\frac{(U_\varphi+u_\varphi)^2}R\approx\frac{U_\varphi^2}R+\frac{2U_\varphi u_\varphi}R$ · *why:*
     expand the square and drop u_φ² (P68); the factor 2 comes from the cross term. · *plain:* a faster ring feels twice the extra
     centrifugal pull.
  3. *did:* Write the radial equation · *tex:* $\frac{\partial u_R}{\partial t}-\frac{2U_\varphi u_\varphi}R=-\frac1\rho\frac{\partial p}{\partial R}+
     \nu\Big(\nabla^2u_R-\frac{u_R}{R^2}\Big)$ · *why:* subtract the base balance $\frac1\rho\frac{dP}{dR}=\frac{U_\varphi^2}R$ (11.49);
     U_R = 0, so advection by the base flow is absent (∂/∂φ = 0). · *plain:* first line of (11.50).
  4. *did:* Write the other three · *tex:* $\frac{\partial u_\varphi}{\partial t}+\Big(\frac{dU_\varphi}{dR}+\frac{U_\varphi}R\Big)u_R=\nu\Big(\nabla^2u_\varphi-
     \frac{u_\varphi}{R^2}\Big),\ \frac{\partial u_z}{\partial t}=-\frac1\rho\frac{\partial p}{\partial z}+\nu\nabla^2u_z,\ \frac1R\frac{\partial(Ru_R)}{\partial R}+
     \frac{\partial u_z}{\partial z}=0$ · *why:* u_R∂_RU_φ from advection and u_RU_φ/R from the Coriolis-like term are the only linear
     products. ⚠️ slip #12: continuity needs 1/R. · *plain:* (11.50).
  5. *did:* Insert normal modes · *tex:* $(u_R,u_\varphi,u_z,p)=(\hat u_R,\hat u_\varphi,\hat u_z,\hat p)(R)\,e^{\mathrm ikz+\sigma t}$ · *why:*
     coefficients depend on R only; k real keeps the disturbance bounded in z. Write $D=\frac d{dR}$, $D_*=\frac d{dR}+\frac1R$. ·
     *plain:* rings stacked along the axis.
  6. *did:* Express the operators · *tex:* $\nabla^2-\frac1{R^2}\to DD_*-k^2,\qquad\nabla^2\to D_*D-k^2$ · *why:* $DD_*f=f''+\frac{f'}R-\frac f{R^2}$
     and $D_*Df=f''+\frac{f'}R$ (product rule); ∂²/∂z² → −k². · *plain:* compact names for the radial operators.
  7. *did:* Solve continuity for û_z · *tex:* $\hat u_z=\frac{\mathrm i}kD_*\hat u_R$ · *why:* $D_*\hat u_R+\mathrm ik\hat u_z=0$. · *plain:* the
     axial flow closes the rings.
  8. *did:* Solve the axial equation for p̂ · *tex:* $\frac{\hat p}\rho=\frac1{k^2}\big(\nu(D_*D-k^2)-\sigma\big)D_*\hat u_R$ · *why:* $\sigma\hat
     u_z=-\frac{\mathrm ik}\rho\hat p+\nu(D_*D-k^2)\hat u_z$, then insert step 7 ($\frac{\mathrm i}k\cdot\frac1{\mathrm ik}=\frac1{k^2}$). ·
     *plain:* the pressure is slaved to the radial flow.
  9. *did:* Use the base vorticity · *tex:* $\frac{dU_\varphi}{dR}+\frac{U_\varphi}R=2A\ \Rightarrow\ \big(\nu(DD_*-k^2)-\sigma\big)\hat u_\varphi=2A\,\hat u_R$ ·
     *why:* $U_\varphi=AR+B/R$ gives $(A-B/R^2)+(A+B/R^2)=2A$; this is the azimuthal equation with normal modes. · *plain:* the radial
     flow feeds the swirl disturbance at a rate set by 2A.
  10. *did:* Insert p̂ into the radial equation · *tex:* $\big(\nu(DD_*-k^2)-\sigma\big)(DD_*-k^2)\hat u_R=2k^2\frac{U_\varphi}R\hat u_\varphi$ · *why:*
      $D(D_*D-k^2)=(DD_*-k^2)D$, so $-D\hat p/\rho=-\frac1{k^2}(\nu(DD_*-k^2)-\sigma)DD_*\hat u_R$; multiply by k² and collect. · *plain:*
      the swirl disturbance drives the radial flow through the centrifugal force.
  11. *did:* Drop curvature (narrow gap) · *tex:* $DD_*\approx\frac{d^2}{dR^2}$ · *why:* the dropped terms are smaller by d/R₁ (P267). ·
      *plain:* the gap looks like a flat channel.
  12. *did:* Linearise the angular speed · *tex:* $\frac{U_\varphi}R=\Omega(R)\approx\Omega_1(1+\alpha x),\ \ x=\frac{R-R_1}d,\ \ \alpha=\frac{\Omega_2}{\Omega_1}-1$ ·
      *why:* Ω goes from Ω₁ to Ω₂ across the gap; over a narrow gap a straight line is the leading approximation. · *plain:* the
      fluid's spin varies linearly across the gap.
  13. *did:* Scale lengths and time · *tex:* $\Big(\frac{d^2}{dx^2}-k^2-\sigma\Big)\Big(\frac{d^2}{dx^2}-k^2\Big)\hat u_R=\frac{2k^2d^2\Omega_1}\nu(1+\alpha x)
      \hat u_\varphi,\ \ \Big(\frac{d^2}{dx^2}-k^2-\sigma\Big)\hat u_\varphi=\frac{2Ad^2}\nu\hat u_R$ · *why:* d/dR = (1/d)d/dx; from now on k
      means kd and σ means σd²/ν (the book does not say so); divide by ν. · *plain:* dimensionless in the gap width.
  14. *did:* Rescale the swirl · *tex:* $\hat u_\varphi\to\frac{\nu}{2k^2d^2\Omega_1}\hat u_\varphi$ · *why:* chosen so the first equation's
      right side becomes (1 + αx)û_φ; the constant moves to the second equation. · *plain:* only one number will remain.
  15. *did:* Read off Ta · *tex:* $\Big(\frac{d^2}{dR^2}-k^2-\sigma\Big)\hat u_\varphi=-\mathrm{Ta}\,k^2\hat u_R,\ \ \mathrm{Ta}=-\frac{4A\Omega_1d^4}{\nu^2}=
      4\Big(\frac{\Omega_1R_1^2-\Omega_2R_2^2}{R_2^2-R_1^2}\Big)\frac{\Omega_1d^4}{\nu^2}$ (11.51), (11.52) · *why:* $\frac{2Ad^2}\nu\cdot\frac{2k^2d^2
      \Omega_1}\nu=\frac{4A\Omega_1d^4}{\nu^2}k^2$ and A of (11.49). The book writes d/dR for d/dx. · *plain:* Taylor's narrow-gap pair
      with the Taylor number.
- **Result.** $\big(\frac{d^2}{dR^2}-k^2-\sigma\big)\big(\frac{d^2}{dR^2}-k^2\big)\hat u_R=(1+\alpha x)\hat u_\varphi$, $\big(\frac{d^2}{dR^2}-k^2-\sigma
  \big)\hat u_\varphi=-\mathrm{Ta}k^2\hat u_R$ (11.51) with Ta (11.52) and $\hat u_R=\frac{d\hat u_R}{dR}=\hat u_\varphi=0$ at x = 0, 1 (11.53) —
  *in words:* centrifugal forcing (Ta) against viscosity, in a flat gap.
- **Check.** Ta is dimensionless ✓; Ta ≤ 0 on and beyond the Rayleigh line (A ≥ 0 there) — no forcing ✓ (verifier note iii); μ → 1
  gives Bénard (D15) ✓; `ch11.taylor_marginal_Ta` reproduces 3389.90 at μ = 0 ✓.
- **sympy check intent (`check_src`).**
  ```python
  import sympy as sp                                         # symbolic algebra
  R, A, B, k = sp.symbols("R A B k", positive=True)          # radius, Couette constants, axial wavenumber
  U = A*R + B/R                                              # (11.49): circular Couette flow
  print(sp.simplify(sp.diff(U, R) + U/R - 2*A))              # step 9: dU/dR + U/R = 2A  -> 0
  u = sp.symbols("u")                                        # a small swirl disturbance
  print(sp.expand((U + u)**2/R) - sp.expand(U**2/R))         # step 2: exact difference = 2Uu/R + u^2/R (the u^2 is dropped)
  f = sp.Function("f")(R)                                    # any radial profile
  D = lambda g: sp.diff(g, R); Ds = lambda g: sp.diff(g, R) + g/R    # D and D_* of step 5
  lhs = D(Ds(D(f)) - k**2*f); rhs = D(Ds(D(f))) - k**2*D(f)  # step 10: D(D_*D - k^2) versus (DD_* - k^2)D
  print(sp.simplify(lhs - rhs))                              # -> 0: the operators can be swapped as claimed
  O1, O2, R1, R2, d, nu = sp.symbols("Omega1 Omega2 R1 R2 d nu", positive=True)
  Aco = (O2*R2**2 - O1*R1**2)/(R2**2 - R1**2)                # A of (11.49)
  Ta_book = 4*(O1*R1**2 - O2*R2**2)/(R2**2 - R1**2)*O1*d**4/nu**2   # (11.52)
  print(sp.simplify(-4*Aco*O1*d**4/nu**2 - Ta_book))         # step 15: Ta = -4 A Omega1 d^4/nu^2  -> 0
  ```
- **What it means.** The forcing is the product of the base flow's vorticity (2A) and its angular speed (Ω) — the centrifugal analogue
  of buoyancy times stratification. Fails for wide gaps (curvature) and for non-axisymmetric (wavy) modes.
- **Traps.** The factor 2 in 2U_φu_φ/R; continuity without 1/R (slip #12); treating k and σ as dimensional after step 13; forgetting
  that Ω(x) is only linear to leading order.

### D15 · The Taylor number's narrow-gap form $2(\Omega_1R_1d/\nu)^2(d/R_1)$ and the Bénard limit μ → 1 — ★, 6 steps, in C07 (notebook · `taylor_couette_onset`)
- **Goal.** Simplify (11.52) for an inner cylinder spinning alone, and show why co-rotating cylinders give Bénard's 1708.
- **Start.** $\mathrm{Ta}\equiv4\Big(\frac{\Omega_1R_1^2-\Omega_2R_2^2}{R_2^2-R_1^2}\Big)\frac{\Omega_1d^4}{\nu^2}$ (11.52) and (11.51) — *in words:* D14's
  results.
- **Plan.** (1) Set Ω₂ = 0 and expand R₂² − R₁². (2) Set α = 0 in (11.51) and compare with (11.39).
- **Tools.** Narrow-gap approximation (P267); D09's (11.39).
- **Assumptions.** d ≪ R₁ (step 3); σ = 0 (step 5).
- **Steps.**
  1. *did:* Set Ω₂ = 0 · *tex:* $\mathrm{Ta}=\frac{4\Omega_1^2R_1^2d^4}{\nu^2(R_2^2-R_1^2)}$ · *why:* only the inner cylinder turns. · *plain:* the
     inner-only Taylor number.
  2. *did:* Factor the denominator · *tex:* $R_2^2-R_1^2=(R_2-R_1)(R_2+R_1)=d(2R_1+d)$ · *why:* difference of squares, R₂ = R₁ + d. ·
     *plain:* gap width times the sum of radii.
  3. *did:* Use the narrow gap · *tex:* $\mathrm{Ta}\approx\frac{4\Omega_1^2R_1^2d^4}{2\nu^2R_1d}=2\Big(\frac{\Omega_1R_1d}\nu\Big)^2\frac d{R_1}$ · *why:*
     2R₁ + d ≈ 2R₁ to leading order in d/R₁ (P267). · *plain:* a Reynolds number squared times the gap ratio.
  4. *did:* Let the cylinders co-rotate (μ → 1) · *tex:* $\alpha=\mu-1\to0:\ \ (1+\alpha x)\to1$ · *why:* the angular speed becomes uniform
     across the gap; we look at the *shape* of (11.51), not at the flow at μ = 1 (where the leading Ta itself → 0). · *plain:* the
     gap's spin no longer varies.
  5. *did:* Compare with Bénard at σ = 0 · *tex:* $\Big(\frac{d^2}{dx^2}-k^2\Big)^2\hat u_R=\hat u_\varphi,\ \Big(\frac{d^2}{dx^2}-k^2\Big)\hat u_\varphi=-\mathrm{Ta}k^2\hat u_R$
     · *why:* put $\hat u_\varphi=\mathrm{Ta}\,k^2\hat T$, $\hat u_R=W$: the pair becomes $(D^2-k^2)^2W=\mathrm{Ta}k^2\hat T$ and $(D^2-k^2)\hat T=-W$ —
     (11.39) with Ta for Ra. · *plain:* the same equations as heated convection.
  6. *did:* Compare the conditions · *tex:* $\mathrm{Ta}_c(\mu\to1)=\mathrm{Ra}_c=1707.76\ \text{at}\ k=3.117$ · *why:* (11.53) is W = W′ = T̂ = 0
     at both walls, as (11.38); same equations and conditions give the same eigenvalues (x ∈ [0, 1] is z ∈ [−½, ½] shifted). Our
     heuristic for (11.54): the gap average of (1 + αx) is 1 + α/2 = ½(1 + μ). · *plain:* 1708 again.
- **Result.** $\mathrm{Ta}\approx2(\Omega_1R_1d/\nu)^2(d/R_1)$ (narrow gap, inner only); $\mathrm{Ta}_c\to1707.76$ as μ → 1, and
  $\mathrm{Ta}_{cr}\approx\frac{1708}{\frac12(1+\Omega_2/\Omega_1)}$ (11.54) — *in words:* rotation's threshold is Bénard's, divided by the
  average spin factor.
- **Check.** R₁ = 0.1 m, d = 5 mm, Ω₁ = 0.5 rad/s: exact 6098 vs shortcut 6250 (2.4 % = O(d/R₁)) ✓; `ch11.taylor_critical(1.0)` = 1707.76 ✓.
- **What it means.** Centrifugal instability is buoyancy in disguise. (11.54) is a good fit for co-rotation (0.77 % at μ = 0), poor
  for counter-rotation (6.5 % at μ = −0.5).
- **Traps.** Reading μ = 1 as "unstable at 1708" — at μ = 1 the flow is solid-body rotation and Ta ≈ 0: nothing drives it.

### D16 · The stratified perturbation equations (11.55)–(11.57) — ★, 7 steps, in C08 (notebook)
- **Goal.** Write the linear equations of small 2-D disturbances on a stratified current U(z), ρ̄(z).
- **Start.** $\frac{\partial\tilde{\mathbf u}}{\partial t}+(\tilde{\mathbf u}\cdot\nabla)\tilde{\mathbf u}=-\frac1{\rho_0}\nabla\tilde p-g\frac{\bar\rho+\rho}{\rho_0}\mathbf e_z$
  and $0=-\frac1{\rho_0}\frac{\partial P}{\partial z}-g\frac{\bar\rho}{\rho_0}$ (p. 502–503), $\frac{D\tilde\rho}{Dt}=0$ — *in words:* inviscid,
  non-diffusive Boussinesq flow and its hydrostatic base.
- **Plan.** (1) Subtract the base state and linearise. (2) Write components. (3) Linearise the density equation. (4) Introduce ψ.
- **Tools.** Linearisation (P68); the material derivative (Ch. 3 (3.5)); N² (R18); the stream function (R19).
- **Assumptions.** Inviscid, non-diffusive, Boussinesq, 2-D disturbances (Squire assumed).
- **Steps.**
  1. *did:* Subtract the base and linearise · *tex:* $\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)(U\mathbf e_x)+U(\mathbf e_x\cdot\nabla)\mathbf u=
     -\frac1{\rho_0}\nabla p-g\frac\rho{\rho_0}\mathbf e_z$ · *why:* insert ũ = Ue_x + u, p̃ = P + p, ρ̃ = ρ̄ + ρ; the base pressure
     cancels the base weight; drop products of disturbances. · *plain:* the disturbance is advected by U and shears U.
  2. *did:* Evaluate (u·∇)(Ue_x) · *tex:* $(\mathbf u\cdot\nabla)(U\mathbf e_x)=w\frac{dU}{dz}\mathbf e_x$ · *why:* U depends on z only. · *plain:*
     vertical motion moves fast fluid down or slow fluid up.
  3. *did:* Take the x-component · *tex:* $\frac{\partial u}{\partial t}+w\frac{dU}{dz}+U\frac{\partial u}{\partial x}=-\frac1{\rho_0}\frac{\partial p}{\partial x}$
     (11.55) first · *why:* e_x-parts of step 1 with step 2. · *plain:* horizontal momentum.
  4. *did:* Take the z-component · *tex:* $\frac{\partial w}{\partial t}+U\frac{\partial w}{\partial x}=-\frac1{\rho_0}\frac{\partial p}{\partial z}-g\frac\rho{\rho_0}$
     (11.55) second · *why:* e_z-parts of step 1. · *plain:* vertical momentum with buoyancy.
  5. *did:* Correct the printed form · *tex:* $-\frac1{\rho_0}\frac{\partial p}{\partial z}\ \text{(not}\ \frac{\partial p}{\partial x})$ · *why:* ⚠️ slip #4: the
     book's (11.55) prints ∂p/∂x in the w-equation; the z-component of ∇p is ∂p/∂z, as (11.57) then uses. · *plain:* a misprint
     the next equation silently corrects.
  6. *did:* Linearise the density equation · *tex:* $\frac{\partial\rho}{\partial t}+U\frac{\partial\rho}{\partial x}+w\frac{d\bar\rho}{dz}=0\ \Leftrightarrow\
     \frac{\partial\rho}{\partial t}+U\frac{\partial\rho}{\partial x}-\frac{\rho_0N^2w}g=0$ (11.56) · *why:* Dρ̃/Dt (Ch. 3 (3.5)) with ρ̄(z); drop
     uρ_x, wρ_z; N² = −(g/ρ₀)dρ̄/dz (7.128). · *plain:* rising fluid brings heavier fluid up.
  7. *did:* Introduce the stream function · *tex:* $u=\frac{\partial\psi}{\partial z},\ w=-\frac{\partial\psi}{\partial x}$ ⇒ (11.57) · *why:* this satisfies
     ∂u/∂x + ∂w/∂z = 0 automatically (R19, §11.7 sign). · *plain:* three equations for ψ, p, ρ.
- **Result.** (11.55)–(11.57) — *in words:* linear advection by U, production by shear (wU′), buoyancy and stratification (N²).
- **Check.** Units: each momentum term m/s², density terms kg/(m³ s) ✓. With U = 0: Ch. 7's internal-wave equations ✓.
- **What it means.** The input to D17. Fails with diffusion or viscosity (C11) and for 3-D stratified disturbances (Squire is not a
  theorem here).
- **Traps.** Slip #4; the sign of w with ψ (w = −∂ψ/∂x here).

### D17 · The Taylor–Goldstein equation (11.61) from (11.58)–(11.60) — ★★, 9 steps, in C08 (notebook · `richardson_shear_instability`)
- **Goal.** Reduce (11.57) to one ODE for the stream-function amplitude ψ̂(z).
- **Start.** (11.57) and the normal modes $[\rho,p,\psi]=[\hat\rho,\hat p,\hat\psi](z)e^{\mathrm ik(x-ct)}$ — *in words:* the stratified
  disturbance equations, one wave at a time.
- **Plan.** (1) Insert the modes. (2) Remove the pressure by differentiating and subtracting. (3) Remove the density.
- **Tools.** Normal modes (C01); differentiation and the product rule (P38); complex algebra (P153).
- **Assumptions.** U(z) ≠ c everywhere (step 7; true for c_i ≠ 0).
- **Steps.**
  1. *did:* Replace the derivatives · *tex:* $\frac{\partial}{\partial x}\to\mathrm ik,\ \ \frac{\partial}{\partial t}\to-\mathrm ikc$ · *why:* the
     mode's exponential (C01). · *plain:* every x, t derivative becomes a factor.
  2. *did:* Insert into the x-equation · *tex:* $(U-c)\hat\psi'-U'\hat\psi=-\frac{\hat p}{\rho_0}$ (11.58) · *why:* $-\mathrm ikc\hat\psi'-\mathrm
     ik\hat\psi U'+\mathrm ikU\hat\psi'=-\mathrm ik\hat p/\rho_0$; divide by ik. · *plain:* horizontal momentum of the mode.
  3. *did:* Insert into the z-equation · *tex:* $k^2(U-c)\hat\psi=-\frac{g\hat\rho}{\rho_0}-\frac{\hat p'}{\rho_0}$ (11.59) · *why:*
     $-\partial_t\partial_x\psi\to-(-\mathrm ikc)(\mathrm ik)\hat\psi=-k^2c\hat\psi$ and $-U\partial_x^2\psi\to k^2U\hat\psi$. · *plain:* vertical
     momentum of the mode.
  4. *did:* Insert into the density equation · *tex:* $(U-c)\hat\rho+\frac{\rho_0N^2}g\hat\psi=0$ (11.60) · *why:* $-\mathrm ikc\hat\rho+\mathrm ikU
     \hat\rho+\frac{\rho_0N^2}g\mathrm ik\hat\psi=0$; divide by ik. · *plain:* the density disturbance from vertical displacement.
  5. *did:* Differentiate (11.58) in z · *tex:* $(U-c)\hat\psi''+U'\hat\psi'-U''\hat\psi-U'\hat\psi'=-\frac{\hat p'}{\rho_0}$ · *why:* product rule
     (P38); we need p̂′ to match (11.59). · *plain:* the two U′ψ̂′ terms cancel.
  6. *did:* Subtract (11.59) · *tex:* $(U-c)(\hat\psi''-k^2\hat\psi)-U''\hat\psi=\frac{g\hat\rho}{\rho_0}$ · *why:* the p̂′ terms cancel; that was
     the aim (the book's single sentence). · *plain:* pressure gone, density left.
  7. *did:* Solve (11.60) for ρ̂ · *tex:* $\hat\rho=-\frac{\rho_0N^2\hat\psi}{g(U-c)}$ · *why:* divide by U − c, allowed only where U ≠ c
     (not at a critical layer, P269). · *plain:* density follows the displacement.
  8. *did:* Substitute · *tex:* $(U-c)(\hat\psi''-k^2\hat\psi)-U''\hat\psi=-\frac{N^2}{U-c}\hat\psi$ · *why:* $\frac g{\rho_0}\hat\rho=-\frac{N^2\hat\psi}{U-c}$. ·
     *plain:* one equation for ψ̂.
  9. *did:* Move everything left · *tex:* $(U-c)\Big(\frac{d^2}{dz^2}-k^2\Big)\hat\psi-\frac{d^2U}{dz^2}\hat\psi+\frac{N^2}{U-c}\hat\psi=0$ (11.61) · *why:*
     add the right side to both sides. · *plain:* the Taylor–Goldstein equation.
- **Result.** (11.61) with $\hat\psi(0)=\hat\psi(d)=0$ (11.62) — *in words:* shear curvature U″ and stratification N² decide the
  eigenvalue c.
- **Check.** N² = 0: Rayleigh's equation (11.81) ✓; units of every term 1/s × 1/m (×ψ̂) ✓; `ch11.stratified_shear_sympy()["tg_residual"]
  == 0` ✓.
- **What it means.** The master equation of inviscid stratified shear instability (C09, C10, Ch. 13). Fails where U = c (critical
  layers need viscosity or a careful limit).
- **Traps.** (ik)² = −k² sign bookkeeping; forgetting that the U′ψ̂′ terms cancel; dividing by U − c without saying it is nonzero.

### D18 · The Miles–Howard criterion: (11.63) → (11.64) → (11.65) → Ri > ¼ everywhere ⇒ stable (11.66)–(11.67) — ★★★, 14 steps, in C09 (notebook · `richardson_shear_instability`)
- **Goal.** Prove that a stratified shear flow with gradient Richardson number above ¼ everywhere has no growing wave.
- **Start.** $(U-c)\big(\frac{d^2}{dz^2}-k^2\big)\hat\psi-\frac{d^2U}{dz^2}\hat\psi+\frac{N^2}{U-c}\hat\psi=0$ (11.61) with $\hat\psi(0)=\hat\psi(d)=0$
  (11.62) — *in words:* the Taylor–Goldstein problem.
- **Plan.** (1) Assume a growing mode (c_i ≠ 0) and substitute φ = ψ̂/(U − c)^{1/2}. (2) Rearrange into a self-adjoint form. (3) Multiply
  by φ*, integrate, integrate by parts. (4) Take the imaginary part and reach a contradiction when Ri > ¼.
- **Tools.** Complex powers (P155); product and chain rules (P38, P49); integration by parts (P218a); real and imaginary parts
  (P260); proof by contradiction (P255); sympy (P40).
- **Assumptions.** c_i ≠ 0, so U − c never vanishes (steps 2, 12); U, N² smooth; rigid lids (step 10).
- **Steps.**
  1. *did:* Suppose a growing mode exists · *tex:* $c_i\neq0\ \Rightarrow\ U(z)-c\neq0\ \text{for all}\ z$ · *why:* proof by contradiction
     (P255): assume instability and look for something impossible; U is real, so U − c has imaginary part −c_i ≠ 0. · *plain:* we may
     divide by U − c and take its square root.
  2. *did:* Substitute the new unknown · *tex:* $\phi\equiv\frac{\hat\psi}{(U-c)^{1/2}},\ \ \hat\psi=(U-c)^{1/2}\phi$ (11.63) · *why:* a clever change of
     variable (Howard) that will make the equation self-adjoint; any fixed branch of the root works since U − c ≠ 0 (P155). ·
     *plain:* φ is not a potential — just a scaled ψ̂.
  3. *did:* Differentiate once · *tex:* $\hat\psi'=(U-c)^{1/2}\phi'+\frac{\phi\,U'}{2(U-c)^{1/2}}$ · *why:* product rule and chain rule:
     $\frac d{dz}(U-c)^{1/2}=\frac{U'}{2(U-c)^{1/2}}$ (P38, P49). · *plain:* the first derivative of ψ̂.
  4. *did:* Differentiate twice · *tex:* $\hat\psi''=(U-c)^{1/2}\phi''+\frac{\phi'U'+\tfrac12\phi U''}{(U-c)^{1/2}}-\frac{\phi U'^2}{4(U-c)^{3/2}}$ ·
     *why:* differentiate step 3 term by term; the last term comes from $\frac d{dz}(U-c)^{-1/2}=-\frac{U'}{2(U-c)^{3/2}}$. · *plain:* the book's
     two derivative formulas (p. 504).
  5. *did:* Multiply (11.61) through · *tex:* $(U-c)\hat\psi''=(U-c)^{3/2}\phi''+(U-c)^{1/2}\big(\phi'U'+\tfrac12\phi U''\big)-\frac{\phi U'^2}{4(U-c)^{1/2}}$ ·
     *why:* (11.61) needs (U − c)ψ̂″; multiply step 4 by U − c. · *plain:* the first term of (11.61) in the new unknown.
  6. *did:* Insert every term into (11.61) · *tex:* $(U-c)^{3/2}\phi''+(U-c)^{1/2}\phi'U'+\big(\tfrac12-1\big)(U-c)^{1/2}U''\phi-k^2(U-c)^{3/2}\phi-\frac{\big(\tfrac14U'^2-N^2\big)\phi}{(U-c)^{1/2}}=0$
     · *why:* $-U''\hat\psi=-U''(U-c)^{1/2}\phi$ and $\frac{N^2\hat\psi}{U-c}=\frac{N^2\phi}{(U-c)^{1/2}}$; collect the two U″ terms. · *plain:* the
     book's "after some rearrangement", written out.
  7. *did:* Divide by (U − c)^{1/2} · *tex:* $(U-c)\phi''+U'\phi'-\Big\{k^2(U-c)+\tfrac12U''+\frac{\tfrac14U'^2-N^2}{U-c}\Big\}\phi=0$ · *why:* nonzero by
     step 1; ½ − 1 = −½ gives the U″ coefficient. · *plain:* every term is now simple.
  8. *did:* Recognise a derivative · *tex:* $\frac d{dz}\Big[(U-c)\frac{d\phi}{dz}\Big]-\Big\{k^2(U-c)+\frac12\frac{d^2U}{dz^2}+\frac{\tfrac14(dU/dz)^2-N^2}{U-c}\Big\}\phi=0$ (11.64)
     · *why:* product rule backwards: $[(U-c)\phi']'=(U-c)\phi''+U'\phi'$. · *plain:* the self-adjoint form.
  9. *did:* Multiply by φ* and integrate · *tex:* $\int_0^d\big[(U-c)\phi'\big]'\phi^*dz-\int_0^d\Big\{k^2(U-c)+\tfrac12U''+\frac{\tfrac14U'^2-N^2}{U-c}\Big\}\lvert\phi\rvert^2dz=0$
     · *why:* φ*φ = ∣φ∣² (P81); integrating turns the ODE into one number equation we can read. · *plain:* an 'energy' identity.
  10. *did:* Integrate the first term by parts · *tex:* $\int_0^d\big[(U-c)\phi'\big]'\phi^*dz=-\int_0^d(U-c)\lvert\phi'\rvert^2dz$ · *why:* P218a; the
      boundary term $[(U-c)\phi'\phi^*]_0^d$ vanishes because ψ̂ = 0 at the lids (11.62) makes φ = 0 there. · *plain:* a curvature
      becomes a gradient squared.
  11. *did:* Rearrange · *tex:* $\int_0^d\frac{N^2-\tfrac14U'^2}{U-c}\lvert\phi\rvert^2dz=\int_0^d(U-c)\big\{\lvert\phi'\rvert^2+k^2\lvert\phi\rvert^2\big\}dz+\int_0^d\tfrac12U''\lvert\phi\rvert^2dz$ (11.65)
      · *why:* steps 9–10; move the N² term to the left (sign flip) and the rest to the right. · *plain:* the book's integral identity.
  12. *did:* Take the imaginary part · *tex:* $c_i\int_0^d\frac{N^2-\tfrac14U'^2}{\lvert U-c\rvert^2}\lvert\phi\rvert^2dz=-c_i\int_0^d\big\{\lvert\phi'\rvert^2+k^2\lvert\phi\rvert^2\big\}dz$
      · *why:* $\mathrm{Im}\frac1{U-c}=\frac{c_i}{\lvert U-c\rvert^2}$ (P260), Im(U − c) = −c_i, the U″ integral is real. The text extraction
      drops this minus; the page has it. · *plain:* c_i times one integral equals −c_i times a positive one.
  13. *did:* Divide by c_i · *tex:* $\int_0^d\frac{N^2-\tfrac14U'^2}{\lvert U-c\rvert^2}\lvert\phi\rvert^2dz=-\int_0^d\big\{\lvert\phi'\rvert^2+k^2\lvert\phi\rvert^2\big\}dz<0$ ·
      *why:* c_i ≠ 0 by step 1; the right side is minus a positive integral (φ ≢ 0). · *plain:* the left integral must be negative.
  14. *did:* Reach the contradiction · *tex:* $N^2>\tfrac14\Big(\frac{dU}{dz}\Big)^2\ \text{everywhere}\ \Leftrightarrow\ \mathrm{Ri}\equiv\frac{N^2}{(dU/dz)^2}>\frac14$ (11.66), (11.67)
      · *why:* then the left integrand is positive and the left side positive, contradicting step 13; so the assumption c_i ≠ 0 was
      false (P255). · *plain:* Ri > ¼ everywhere ⇒ no growing mode.
- **Result.** If $\mathrm{Ri}(z)=N^2/(dU/dz)^2>\tfrac14$ everywhere, the flow is linearly stable (11.67) — *in words:* enough stratification
  everywhere guarantees stability; Ri < ¼ somewhere is necessary (not sufficient) for instability.
- **Check.** Units: N² and U′² both 1/s² ✓. Limit N² = 0: the left side is negative definite — no conclusion, consistent with
  inflectional instability (C12). tanh layer with J = 0.3 (Ri_min = 0.3): no growth computed ✓; identity (11.65) holds to 1e-8 for a
  computed mode ✓.
- **sympy check intent (`check_src`).**
  ```python
  import sympy as sp                                            # symbolic algebra
  z, k, c = sp.symbols("z k c")                                 # height, wavenumber, complex wave speed
  U, N2, phi = sp.Function("U")(z), sp.Function("N2")(z), sp.Function("phi")(z)   # generic profiles and the new unknown
  s = sp.sqrt(U - c)                                            # (U - c)^(1/2) of (11.63)
  psi = s*phi                                                   # step 2: psi-hat in terms of phi
  TG = (U - c)*(sp.diff(psi, z, 2) - k**2*psi) - sp.diff(U, z, 2)*psi + N2*psi/(U - c)     # (11.61)
  SA = sp.diff((U - c)*sp.diff(phi, z), z) - (k**2*(U - c) + sp.diff(U, z, 2)/2
        + (sp.diff(U, z)**2/4 - N2)/(U - c))*phi                # (11.64)
  print(sp.simplify(TG/s - SA))                                 # steps 3-8: (11.61) divided by (U-c)^(1/2) is (11.64) -> 0
  cr, ci, Ur = sp.symbols("c_r c_i U_r", real=True)             # real parts for step 12
  print(sp.simplify(sp.im(1/(Ur - cr - sp.I*ci)) - ci/((Ur - cr)**2 + ci**2)))   # Im{1/(U-c)} = c_i/|U-c|^2 -> 0
  ```
- **What it means.** Ri = ¼ is a *guarantee*, not a trigger — it is why mixing parameterisations switch shear mixing off above about ¼.
  Fails for non-Boussinesq or viscous/diffusive flows, and says nothing about finite-amplitude instabilities.
- **Traps.** Dropping the ½U″ or the −¼U′²/(U − c) term in step 6; forgetting that φ = 0 at the lids follows from (11.62); the minus
  sign in step 12; claiming Ri < ¼ implies instability.

### D19 · Howard's semicircle (from (11.68)–(11.72)) and the growth bound $kc_i\le\frac k2(U_{\max}-U_{\min})$ — ★★★, 14 steps, in C10 (notebook · `inviscid_shear_criteria`)
- **Goal.** Show that the wave speed of any unstable mode lies inside the half-disc spanned by the slowest and fastest fluid.
- **Start.** (11.61) with (11.62), a growing mode c_i ≠ 0, and stable stratification N² ≥ 0 — *in words:* the Taylor–Goldstein
  problem again.
- **Plan.** (1) Substitute F = ψ̂/(U − c) and write a divergence form. (2) Multiply by F*, integrate: an identity with a positive weight
  Q. (3) Take real and imaginary parts. (4) Combine with the trivial inequality (U − U_min)(U − U_max) ≤ 0 and complete the square.
- **Tools.** Product rule (P38); integration by parts (P218a); real and imaginary parts (P260); completing the square into a circle
  (P272).
- **Assumptions.** c_i ≠ 0 (steps 1, 9, 13); N² ≥ 0 (step 10 — the book uses it without saying); rigid lids.
- **Steps.**
  1. *did:* Substitute the new unknown · *tex:* $F\equiv\frac{\hat\psi}{U-c},\ \ \hat\psi=(U-c)F$ (11.68) · *why:* allowed because U − c ≠ 0
     (c_i ≠ 0); this choice will cancel the U″ term. · *plain:* ψ̂ seen in the frame moving with the wave.
  2. *did:* Differentiate · *tex:* $\hat\psi'=(U-c)F'+U'F,\ \ \hat\psi''=(U-c)F''+2U'F'+U''F$ · *why:* product rule twice (P38). · *plain:* the book's
     derivative formulas.
  3. *did:* Insert into (11.61) · *tex:* $(U-c)\big[(U-c)F''+2U'F'+U''F-k^2(U-c)F\big]-U''(U-c)F+N^2F=0$ · *why:* every ψ̂ replaced; the
     N²/(U − c) term becomes N²F. · *plain:* (11.61) in F.
  4. *did:* Cancel the U″ terms · *tex:* $(U-c)\big[(U-c)F''+2U'F'-k^2(U-c)F\big]+N^2F=0$ · *why:* +U″(U − c)F and −U″(U − c)F cancel —
     the reason for the substitution. · *plain:* the profile's curvature has gone.
  5. *did:* Write a divergence form · *tex:* $\frac d{dz}\Big[(U-c)^2\frac{dF}{dz}\Big]-k^2(U-c)^2F+N^2F=0$ · *why:* $[(U-c)^2F']'=(U-c)^2F''+2(U-c)U'F'$,
     exactly the first two terms. ⚠️ slip #11: the book prints −k²(U − c)F; expanding step 4 gives (U − c)². · *plain:* a flux-like
     form ready for integration by parts.
  6. *did:* Multiply by F* and integrate by parts · *tex:* $-\int(U-c)^2\lvert F'\rvert^2dz-k^2\int(U-c)^2\lvert F\rvert^2dz+\int N^2\lvert F\rvert^2dz=0$ ·
     *why:* P218a; the boundary term $[(U-c)^2F'F^*]$ vanishes because F = 0 at the lids (ψ̂ = 0). · *plain:* one complex number
     equation.
  7. *did:* Name the weight · *tex:* $\int(U-c)^2Q\,dz=\int N^2\lvert F\rvert^2dz,\ \ Q\equiv\lvert F'\rvert^2+k^2\lvert F\rvert^2\ge0$ · *why:* collect step
     6; Q is a sum of squares (P260). · *plain:* a positive weight across the layer.
  8. *did:* Expand (U − c)² · *tex:* $(U-c)^2=(U-c_r)^2-c_i^2-2\mathrm ic_i(U-c_r)$ · *why:* c = c_r + ic_i, square the binomial. · *plain:* real
     and imaginary parts of the factor.
  9. *did:* Take the imaginary part · *tex:* $c_i\int(U-c_r)Q\,dz=0\ \Rightarrow\ \int UQ\,dz=c_r\int Q\,dz$ (11.70) · *why:* the right side of step 7
     is real; c_i ≠ 0. Since Q ≥ 0 and U − c_r integrates to zero against it, U − c_r changes sign: $U_{\min}<c_r<U_{\max}$ (11.71). ·
     *plain:* an unstable wave travels at a speed found somewhere in the flow.
  10. *did:* Take the real part · *tex:* $\int\big[(U-c_r)^2-c_i^2\big]Q\,dz=\int N^2\lvert F\rvert^2dz\ \ge0$ (11.69) · *why:* the real part of step 7
      with step 8; the right side is ≥ 0 because N² ≥ 0 — the hidden assumption. · *plain:* a second, real identity.
  11. *did:* Eliminate c_r's cross term · *tex:* $\int\big[U^2-c_r^2-c_i^2\big]Q\,dz\ \ge0$ (11.72) · *why:* expand (U − c_r)² and replace
      $\int UQ$ by $c_r\int Q$ (step 9): $-2c_r\int UQ+c_r^2\int Q=-c_r^2\int Q$. (The book writes > 0.) · *plain:* ∫U²Q is at least
      (c_r² + c_i²)∫Q.
  12. *did:* Use the velocity range, weighted · *tex:* $\int(U-U_{\min})(U-U_{\max})\,Q\,dz\le0$ · *why:* U − U_min ≥ 0 and U − U_max ≤ 0 everywhere,
      Q ≥ 0. ⚠️ slip #8: Q must be inside from the start and the inequality is ≤. · *plain:* the flow cannot be faster than its
      fastest part.
  13. *did:* Combine with steps 9 and 11 · *tex:* $\big[c_r^2+c_i^2-c_r(U_{\max}+U_{\min})+U_{\max}U_{\min}\big]\int Q\,dz\le0$ · *why:* expand step 12:
      $\int U^2Q\le(U_{\max}+U_{\min})\int UQ-U_{\max}U_{\min}\int Q$; replace ∫U²Q from below by step 11 and ∫UQ by step 9. · *plain:* a
      condition on c alone.
  14. *did:* Complete the square · *tex:* $\Big[c_r-\tfrac12(U_{\max}+U_{\min})\Big]^2+c_i^2\le\Big[\tfrac12(U_{\max}-U_{\min})\Big]^2$ · *why:* ∫Q > 0, so the
      bracket is ≤ 0; add and subtract ¼(U_max + U_min)² (P272), using $\tfrac14(U_{\max}+U_{\min})^2-U_{\max}U_{\min}=\tfrac14(U_{\max}-U_{\min})^2$.
      · *plain:* a half-disc on [U_min, U_max]; hence $kc_i\le\frac k2(U_{\max}-U_{\min})$.
- **Result.** Howard's semicircle (p. 507) and $kc_i\le\frac k2(U_{\max}-U_{\min})$ — *in words:* every unstable wave speed sits under the
  arc whose diameter is the velocity range; growth is bounded by half the velocity jump times k.
- **Check.** Units ✓. tanh layer: c = 0.4267i at k = 0.445, inside the unit half-disc ✓; all computed unstable TG eigenvalues pass
  `ST.in_howard_semicircle` ✓; (11.69)–(11.70) hold to 1e-8 ✓.
- **sympy check intent (`check_src`).**
  ```python
  import sympy as sp                                           # symbolic algebra
  z, k, c = sp.symbols("z k c")                                # height, wavenumber, complex speed
  U, N2, F = sp.Function("U")(z), sp.Function("N2")(z), sp.Function("F")(z)   # generic profiles and the unknown
  psi = (U - c)*F                                              # step 1: (11.68)
  TG = (U - c)*(sp.diff(psi, z, 2) - k**2*psi) - sp.diff(U, z, 2)*psi + N2*psi/(U - c)   # (11.61)
  DIV = sp.diff((U - c)**2*sp.diff(F, z), z) - k**2*(U - c)**2*F + N2*F                   # step 5, corrected
  DIV_printed = sp.diff((U - c)**2*sp.diff(F, z), z) - k**2*(U - c)*F + N2*F              # slip #11 as printed
  print(sp.simplify(TG - DIV))                                 # -> 0: the corrected divergence form is (11.61)
  print(sp.simplify(TG - DIV_printed) == 0)                    # False: the printed form is not
  cr, ci, Umax, Umin = sp.symbols("c_r c_i U_max U_min", real=True)
  circle = (cr - (Umax + Umin)/2)**2 + ci**2 - ((Umax - Umin)/2)**2   # step 14's left minus right side
  print(sp.expand(circle - (cr**2 + ci**2 - cr*(Umax + Umin) + Umax*Umin)))   # -> 0: completing the square
  ```
- **What it means.** A search window for every eigen-solver (C08–C13) and a hard growth bound; c_r inside the velocity range means
  a critical layer exists for neutral modes (N96). The same semicircle holds for the quasi-geostrophic problems of Ch. 13. Fails if
  N² < 0 somewhere (step 10).
- **Traps.** Slip #11 (U − c instead of (U − c)²); slip #8 (the weight); forgetting the N² ≥ 0 assumption; using step 9 only once.

### D20 · The perturbation equations of a viscous parallel flow (11.73)–(11.77) — ★, 7 steps, in C11 (notebook)
- **Goal.** Write the linear equations for small 3-D disturbances on a parallel flow U(y).
- **Start.** $\frac{\partial u}{\partial t}+(U+u)\frac{\partial}{\partial x}(U+u)+v\frac{\partial}{\partial y}(U+u)=-\frac{\partial}{\partial x}(P+p)+\frac1{\mathrm{Re}}
  \nabla^2(U+u)$ (11.73) and the base balance $0=-\frac{\partial P}{\partial x}+\frac1{\mathrm{Re}}\nabla^2U$ — *in words:* dimensionless
  Navier–Stokes for the total and for the basic flow.
- **Plan.** (1) Subtract the base. (2) Drop products. (3) Repeat for y, z and continuity. (4) Insert 3-D normal modes.
- **Tools.** Dimensionless NS (R21); linearisation (P68); normal modes (C01).
- **Assumptions.** Parallel basic flow — exact for Poiseuille and Couette, approximate for a growing boundary layer (N109).
- **Steps.**
  1. *did:* Use U = U(y) · *tex:* $(U+u)\frac{\partial(U+u)}{\partial x}=U\frac{\partial u}{\partial x}+u\frac{\partial u}{\partial x},\ \ v\frac{\partial(U+u)}{\partial y}=v\frac{dU}{dy}+v\frac{\partial u}{\partial y}$
     · *why:* ∂U/∂x = 0 for a parallel flow. · *plain:* the base flow only advects and is sheared.
  2. *did:* Subtract the base balance · *tex:* $\frac{\partial u}{\partial t}+U\frac{\partial u}{\partial x}+v\frac{dU}{dy}+u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=-\frac{\partial p}{\partial x}+\frac1{\mathrm{Re}}\nabla^2u$
     · *why:* the P and U terms satisfy the base equation exactly (Poiseuille: −dP/dx = 2/Re). · *plain:* only the disturbance is left.
  3. *did:* Drop products · *tex:* $\frac{\partial u}{\partial t}+U\frac{\partial u}{\partial x}+v\frac{\partial U}{\partial y}=-\frac{\partial p}{\partial x}+\frac1{\mathrm{Re}}\nabla^2u$ (11.74)
     · *why:* u∂_xu and v∂_yu are quadratic (P68). · *plain:* linear x-momentum.
  4. *did:* Do the same for y, z, continuity · *tex:* $\frac{\partial v}{\partial t}+U\frac{\partial v}{\partial x}=-\frac{\partial p}{\partial y}+\frac{\nabla^2v}{\mathrm{Re}},\ \frac{\partial w}{\partial t}+U\frac{\partial w}{\partial x}=
     -\frac{\partial p}{\partial z}+\frac{\nabla^2w}{\mathrm{Re}},\ \frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}+\frac{\partial w}{\partial z}=0$ (11.75) · *why:*
     the base has no v, w and its pressure varies only in x, so nothing is subtracted but the advection U∂_x remains. · *plain:* the
     other three equations.
  5. *did:* Insert 3-D normal modes · *tex:* $[\mathbf u,p]=[\hat{\mathbf u}(y),\hat p(y)]\exp\{\mathrm i(kx+mz-kct)\}$ (11.76) · *why:* coefficients depend on y
     only; k, m ≥ 0 without loss of generality. · *plain:* an oblique travelling wave.
  6. *did:* Replace the operators · *tex:* $\partial_t\to-\mathrm ikc,\ \partial_x\to\mathrm ik,\ \partial_z\to\mathrm im,\ \nabla^2\to\frac{d^2}{dy^2}-(k^2+m^2)$ · *why:*
     derivatives of the exponential (C01). · *plain:* numbers for x, z, t derivatives.
  7. *did:* Collect · *tex:* $\mathrm ik(U-c)\hat u+\hat vU'=-\mathrm ik\hat p+\frac{\hat u''-(k^2+m^2)\hat u}{\mathrm{Re}}$, … , $\mathrm ik\hat u+\hat v'+\mathrm im\hat w=0$ (11.77)
     · *why:* −ikc + ikU = ik(U − c) in every advective term. · *plain:* the 3-D normal-mode equations.
- **Result.** (11.77) — *in words:* four ODEs in y for û, v̂, ŵ, p̂ with eigenvalue c.
- **Check.** Units: dimensionless ✓; Re → ∞ removes all second derivatives (C12's limit) ✓.
- **What it means.** The input to Squire (N88) and Orr–Sommerfeld (D21). Fails for a growing boundary layer at O(1/Re) (N109).
- **Traps.** Forgetting the vU′ term; thinking the base pressure gradient must be zero (it balances viscosity).

### D21 · The Orr–Sommerfeld equation (11.79) from (11.77) — ★★, 10 steps, in C11 (notebook · `orr_sommerfeld_neutral_curve`)
- **Goal.** Eliminate the pressure and reduce the 2-D disturbance equations to one equation for the stream-function amplitude φ(y)
  — the book writes only "this effort yields".
- **Start.** (11.77) with m = ŵ = 0 (Squire): $\mathrm ik(U-c)\hat u+\hat vU'=-\mathrm ik\hat p+\frac1{\mathrm{Re}}(\hat u''-k^2\hat u)$,
  $\mathrm ik(U-c)\hat v=-\hat p'+\frac1{\mathrm{Re}}(\hat v''-k^2\hat v)$, $\mathrm ik\hat u+\hat v'=0$ — *in words:* 2-D x-, y-momentum and
  continuity of one wave.
- **Plan.** (1) Use the stream function. (2) Differentiate the x-equation in y and multiply the y-equation by ik. (3) Subtract to remove
  p̂. (4) Divide by ik.
- **Tools.** Stream function (R22); product rule (P38); complex algebra (P153); operator elimination (P177); sympy.
- **Assumptions.** Parallel flow (D20); 2-D waves (Squire, N88).
- **Steps.**
  1. *did:* Use the stream function · *tex:* $\hat u=\phi',\ \ \hat v=-\mathrm ik\phi$ · *why:* u = ∂ψ/∂y, v = −∂ψ/∂x with ψ = φe^{ik(x−ct)} (R22);
     continuity ikφ′ − ikφ′ = 0 holds automatically. · *plain:* one unknown instead of two.
  2. *did:* Rewrite the x-equation · *tex:* $\mathrm ik(U-c)\phi'-\mathrm ikU'\phi=-\mathrm ik\hat p+\frac1{\mathrm{Re}}(\phi'''-k^2\phi')$ · *why:* insert
     step 1; v̂U′ = −ikU′φ. · *plain:* x-momentum in φ.
  3. *did:* Rewrite the y-equation · *tex:* $k^2(U-c)\phi=-\hat p'-\frac{\mathrm ik}{\mathrm{Re}}(\phi''-k^2\phi)$ · *why:* ik(U − c)(−ikφ) = k²(U − c)φ. ·
     *plain:* y-momentum in φ.
  4. *did:* Differentiate step 2 in y · *tex:* $\mathrm ik(U-c)\phi''-\mathrm ikU''\phi=-\mathrm ik\hat p'+\frac1{\mathrm{Re}}(\phi''''-k^2\phi'')$ · *why:*
     product rule: ik[U′φ′ + (U − c)φ″] − ik[U″φ + U′φ′]; the U′φ′ terms cancel. We need ik p̂′, which the y-equation can supply. ·
     *plain:* now both equations contain p̂′.
  5. *did:* Multiply step 3 by ik · *tex:* $\mathrm ik^3(U-c)\phi=-\mathrm ik\hat p'+\frac{k^2}{\mathrm{Re}}(\phi''-k^2\phi)$ · *why:* (ik)(−ik) = k²; this makes
     the pressure terms identical. · *plain:* the same −ik p̂′ in both.
  6. *did:* Subtract step 5 from step 4 · *tex:* $\mathrm ik(U-c)(\phi''-k^2\phi)-\mathrm ikU''\phi=\frac1{\mathrm{Re}}\big(\phi''''-k^2\phi''-k^2\phi''+k^4\phi\big)$
     · *why:* the pressure cancels — the aim; group the (U − c) terms. · *plain:* pressure eliminated.
  7. *did:* Collect the viscous terms · *tex:* $\phi''''-2k^2\phi''+k^4\phi=\Big(\frac{d^2}{dy^2}-k^2\Big)^2\phi$ · *why:* add like terms; it is the
     biharmonic of the mode. · *plain:* viscosity acts on the vorticity's curvature.
  8. *did:* Divide by ik · *tex:* $(U-c)\Big(\frac{d^2\phi}{dy^2}-k^2\phi\Big)-\frac{d^2U}{dy^2}\phi=\frac1{\mathrm ik\mathrm{Re}}\Big(\frac{d^4\phi}{dy^4}-2k^2\frac{d^2\phi}
     {dy^2}+k^4\phi\Big)$ (11.79) · *why:* k > 0. · *plain:* the Orr–Sommerfeld equation.
  9. *did:* Read it as a vorticity balance · *tex:* $\hat\omega=-(\phi''-k^2\phi)$ · *why:* the disturbance vorticity of a stream function is
     −∇²ψ; (11.79) says: vorticity advected by U − c, plus the base vorticity gradient advected by v, equals its viscous diffusion. ·
     *plain:* a vorticity equation (the book's remark).
  10. *did:* State the conditions · *tex:* $\phi=\frac{d\phi}{dy}=0\ \text{at}\ y=y_1,y_2$ (11.80) · *why:* no slip: v̂ = −ikφ = 0 and û = φ′ = 0 on
      each wall. · *plain:* four conditions for a fourth-order equation; c is the eigenvalue.
- **Result.** (11.79) with (11.80) — *in words:* inviscid Rayleigh operator on the left, viscosity ∝ 1/(ikRe) on the right.
- **Check.** Re → ∞ gives (11.81) ✓; planted v̂ = +ikφ changes the sign of the viscous side and fails `ch11.os_derivation_sympy` ✓;
  Orszag's c at Re = 10⁴, k = 1 ✓. Optional `check_src`: the sympy elimination of `os_derivation_sympy()`.
- **What it means.** The viscous term carries an i: it breaks the c ↔ c* symmetry (N92) and lets viscosity both damp and destabilise.
- **Traps.** v̂ = +ikφ (wrong sign); forgetting to differentiate the x-equation (pressure would not cancel); dropping the factor 2 in
  2k²φ″.

### D22 · Rayleigh's inflection-point theorem (11.83)–(11.84) — ★★, 9 steps, in C12 (notebook · `inviscid_shear_criteria`)
- **Goal.** Show that an inviscid parallel flow can have a growing wave only if its velocity profile has an inflection point.
- **Start.** $(U-c)\big(\frac{d^2\phi}{dy^2}-k^2\phi\big)-\frac{d^2U}{dy^2}\phi=0$ (11.81), φ = 0 at y₁, y₂ (11.82) — *in words:* Rayleigh's equation.
- **Plan.** (1) Assume a growing mode and divide by U − c. (2) Multiply by φ*, integrate, by parts. (3) Imaginary part.
- **Tools.** Integration by parts (P218a); real and imaginary parts (P260); proof by contradiction (P255).
- **Assumptions.** c_i ≠ 0 (steps 1, 8).
- **Steps.**
  1. *did:* Suppose c_i ≠ 0 · *tex:* $U-c\neq0\ \text{for all}\ y$ · *why:* U is real, so Im(U − c) = −c_i ≠ 0 (P255). · *plain:* we may divide.
  2. *did:* Divide by U − c · *tex:* $\frac{d^2\phi}{dy^2}-k^2\phi-\frac1{U-c}\frac{d^2U}{dy^2}\phi=0$ · *why:* step 1. · *plain:* the book's
     rewritten form.
  3. *did:* Multiply by φ* and integrate · *tex:* $\int\phi^*\phi''dy-k^2\int\lvert\phi\rvert^2dy-\int\frac{U''}{U-c}\lvert\phi\rvert^2dy=0$ · *why:* conjugate
     times function gives a modulus (P81); integrate over [y₁, y₂]. · *plain:* an identity between numbers.
  4. *did:* Integrate by parts · *tex:* $\int\phi^*\phi''dy=\big[\phi^*\phi'\big]_{y_1}^{y_2}-\int\lvert\phi'\rvert^2dy=-\int\lvert\phi'\rvert^2dy$ · *why:* P218a;
     φ = 0 at the walls (11.82). · *plain:* curvature becomes gradient squared.
  5. *did:* Collect · *tex:* $\int(\lvert\phi'\rvert^2+k^2\lvert\phi\rvert^2)dy+\int\frac1{U-c}\frac{d^2U}{dy^2}\lvert\phi\rvert^2dy=0$ (11.83) · *why:* multiply
     step 3 by −1 after step 4. · *plain:* the book's identity.
  6. *did:* Use the conjugate trick · *tex:* $\frac1{U-c}=\frac{U-c^*}{\lvert U-c\rvert^2}=\frac{U-c_r+\mathrm ic_i}{\lvert U-c\rvert^2}$ · *why:* multiply top and
     bottom by the conjugate (P260). · *plain:* the imaginary part is now visible.
  7. *did:* Take the imaginary part · *tex:* $c_i\int\frac1{\lvert U-c\rvert^2}\frac{d^2U}{dy^2}\lvert\phi\rvert^2dy=0$ (11.84) · *why:* the first integral of
     (11.83) is real; only the c_i part of step 6 is imaginary. · *plain:* c_i times a weighted integral of U″ is zero.
  8. *did:* Divide by c_i · *tex:* $\int\frac{\lvert\phi\rvert^2}{\lvert U-c\rvert^2}\frac{d^2U}{dy^2}dy=0$ · *why:* c_i ≠ 0 (step 1). · *plain:* the weighted average
     of U″ vanishes.
  9. *did:* Conclude · *tex:* $U''\ \text{changes sign in}\ (y_1,y_2)$ · *why:* the weight ∣φ∣²/∣U − c∣² is positive inside; a function of one sign
     cannot integrate to zero against it. · *plain:* an inflection point is necessary for inviscid instability.
- **Result.** Rayleigh's theorem: a growing inviscid mode requires U″ = 0 with a sign change somewhere inside the flow — *in words:*
  inflection-free profiles (Couette, Poiseuille, favourable boundary layers) are inviscidly stable.
- **Check.** Poiseuille (U″ = −2): no unstable Rayleigh eigenvalue found ✓; tanh: inflection at 0 and growth ✓; identity (11.83)
  holds for a computed mode (`rayleigh_identity_check`) ✓.
- **What it means.** Classifies profiles without solving; Fjørtoft (N93) sharpens it. Ch. 13: with the β-effect U″ is replaced by
  U″ − β (Rayleigh–Kuo). Necessary, not sufficient (N95).
- **Traps.** Forgetting that the theorem needs c_i ≠ 0; reading "inflection point" as sufficient; dropping the boundary term
  without citing φ = 0.

### D23 · The disturbance kinetic-energy equation (11.88) and its 2-D form from (11.96) — ★★, 12 steps, in C14 (notebook · `orr_sommerfeld_neutral_curve`)
- **Goal.** Derive how a disturbance's kinetic energy changes (the book leaves it to Exercise 11.13): what produces it, what
  destroys it.
- **Start.** $\frac{\partial}{\partial t}(U_i+u_i)+(U_j+u_j)\frac{\partial}{\partial x_j}(U_i+u_i)=-\frac1\rho\frac{\partial}{\partial x_i}(P+p)+\nu\frac{\partial^2}
  {\partial x_j\partial x_j}(U_i+u_i)$ (11.96) — *in words:* Navier–Stokes for basic flow plus disturbance, in index form.
- **Plan.** (1) Subtract the basic flow's equation. (2) Multiply by u_i. (3) Turn every term that can be into a divergence. (4)
  Integrate over a box of whole wavelengths: the divergences vanish.
- **Tools.** Index notation and the divergence theorem (Ch. 2); product rule (P38); the integral of a periodic x-derivative (P273).
- **Assumptions.** Incompressible (∂_jU_j = ∂_ju_j = 0); steady basic flow; constant ν; u = 0 on walls or far away; a box of an
  integer number of wavelengths (Fig. 11.25).
- **Steps.** (∂_j ≡ ∂/∂x_j, summation over repeated indices)
  1. *did:* Subtract the basic-flow equation · *tex:* $\partial_tu_i+U_j\partial_ju_i+u_j\partial_jU_i+u_j\partial_ju_i=-\frac1\rho\partial_ip+\nu\partial_j\partial_ju_i$ ·
     *why:* U, P satisfy (11.96) without u, p (steady). We keep the u_j∂_ju_i term — it will vanish on its own. · *plain:* the exact
     disturbance equation.
  2. *did:* Multiply by u_i · *tex:* $\partial_t\big(\tfrac12u_iu_i\big)+u_iU_j\partial_ju_i+u_iu_j\partial_jU_i+u_iu_j\partial_ju_i=-\frac1\rho u_i\partial_ip+\nu u_i\partial_j\partial_ju_i$ ·
     *why:* u_i∂_tu_i = ∂_t(½u_iu_i) (chain rule). · *plain:* an equation for the disturbance's kinetic energy per unit mass.
  3. *did:* Rewrite the advection by U · *tex:* $u_iU_j\partial_ju_i=\partial_j\big(\tfrac12u_iu_iU_j\big)$ · *why:* product rule:
     ∂_j(½u_i²U_j) = U_j∂_j(½u_i²) + ½u_i²∂_jU_j, and ∂_jU_j = 0. · *plain:* the base flow only carries energy around.
  4. *did:* Rewrite self-advection · *tex:* $u_iu_j\partial_ju_i=\partial_j\big(\tfrac12u_iu_iu_j\big)$ · *why:* the same move with ∂_ju_j = 0. · *plain:* the
     disturbance moves its own energy around.
  5. *did:* Rewrite the pressure work · *tex:* $u_i\partial_ip=\partial_i(pu_i)$ · *why:* product rule and ∂_iu_i = 0. · *plain:* pressure only
     moves energy.
  6. *did:* Split the viscous term · *tex:* $\nu u_i\partial_j\partial_ju_i=\partial_j(\nu u_i\partial_ju_i)-\nu(\partial_ju_i)(\partial_ju_i)$ · *why:* product
     rule backwards. · *plain:* viscous transport plus a loss.
  7. *did:* Collect · *tex:* $\partial_t\big(\tfrac12u_i^2\big)=-u_iu_j\partial_jU_i-\nu(\partial_ju_i)^2-\partial_j\Big[\tfrac12u_i^2U_j+\tfrac12u_i^2u_j+\frac{pu_j}\rho-\nu u_i\partial_ju_i\Big]$
     · *why:* steps 2–6. · *plain:* production + loss + transport.
  8. *did:* Integrate over the fixed box · *tex:* $\frac d{dt}\int\tfrac12u_i^2dV=-\int u_iu_j\partial_jU_i\,dV-\nu\int(\partial_ju_i)^2dV-\oint[\dots]_jn_j\,dA$ ·
     *why:* the box does not move, so d/dt passes outside; Gauss's divergence theorem (Ch. 2) turns the bracket into a surface
     integral. · *plain:* only the surface can carry transport terms.
  9. *did:* Remove the walls · *tex:* $u_i=0\ \text{on the walls (or far away)}$ · *why:* no slip (or decay) makes every term of the bracket
     vanish there. · *plain:* nothing crosses the walls.
  10. *did:* Remove the ends · *tex:* $\oint_{\text{ends}}=0$ · *why:* the box spans an integer number of wavelengths: the bracket repeats,
      and the outward normals are opposite, so the two ends cancel (P273). · *plain:* what flows out at one end flows in at the other.
  11. *did:* Name the dissipation · *tex:* $\frac d{dt}\int\tfrac12u_i^2dV=-\int u_iu_j\frac{\partial U_i}{\partial x_j}dV-\Lambda,\ \ \Lambda=\nu\int\Big(\frac{\partial u_i}
      {\partial x_j}\Big)^2dV\ge0$ (11.88) · *why:* steps 8–10. · *plain:* production by the Reynolds stress against the shear, minus
      viscous dissipation.
  12. *did:* Specialise to 2-D shear · *tex:* $\frac d{dt}\int\tfrac12(u^2+v^2)dV=-\int uv\frac{\partial U}{\partial y}dV-\Lambda$ · *why:* with U = (U(y), 0, 0)
      only ∂U_1/∂x_2 = U′ is nonzero, so only u_1u_2∂_2U_1 = uvU′ survives. · *plain:* the 2-D form (p. 520).
- **Result.** (11.88) and its 2-D form — *in words:* a disturbance grows only if its Reynolds stress −⟨uv⟩ works against the mean shear
  faster than viscosity dissipates it.
- **Check.** Units: m⁵/s³ per unit density on every term ✓. Computed TS mode at Re = 10⁴, k = 1: dE/dt = P − Λ to 1e-12, P/Λ = 1.616 ✓.
  Optional `check_src`: `ch11.energy_equation_sympy()["residual"] == 0`.
- **What it means.** The production term −⟨uv⟩U′ is the turbulence production of Ch. 12; viscosity can change the u–v phase so that
  production becomes positive (C14). Fails if the box is not an integer number of wavelengths or the walls move.
- **Traps.** Forgetting ∂_jU_j = 0 in step 3; the sign of the dissipation (it is a loss); claiming the self-advection term contributes
  (it is a pure divergence).

### D24 · The Lorenz system (11.91) from the truncation (11.90): $r=\mathrm{Ra}/\mathrm{Ra}_c(k)$, $b=4\pi^2/(\pi^2+k^2)$ — ★★★, 15 steps, in C15 (notebook · `lorenz_attractor`)
- **Goal.** Derive Lorenz's three equations from free–free convection — the book writes only "Lorenz finally obtained".
- **Start.** 2-D Boussinesq convection between stress-free isothermal walls (C05), rolls invariant in y, $u=-\frac{\partial\psi}{\partial z}$,
  $w=\frac{\partial\psi}{\partial x}$, and the truncation $\psi\propto X(t)\cos(\pi z)\sin(kx)$, $T'\propto Y(t)\cos(\pi z)\cos(kx)+Z(t)\sin(2\pi z)$ (11.90),
  z ∈ [−½, ½] — *in words:* one roll and two temperature shapes.
- **Plan.** (1) Write the vorticity and temperature equations with ψ. (2) Make them dimensionless. (3) Insert (11.90) and project
  onto each shape (Galerkin). (4) Rescale time and amplitudes so the constants vanish.
- **Tools.** Vorticity form (Ch. 5); Galerkin truncation (P275); orthogonality of sines (P151); (11.44) from C05; sympy.
- **Assumptions.** Stress-free walls; only the three modes kept — valid slightly above onset (step 10); 2-D rolls.
- **Steps.**
  1. *did:* Take the curl of the momentum equation · *tex:* $\frac{D}{Dt}\nabla^2\psi=g\alpha\frac{\partial T'}{\partial x}+\nu\nabla^4\psi$ · *why:*
     ∂_z of the x-equation minus ∂_x of the z-equation removes the pressure; with this sign convention (slip #9) the y-vorticity is
     −∇²ψ. · *plain:* the vorticity of the roll is made by horizontal temperature differences.
  2. *did:* Write the temperature equation · *tex:* $\frac{DT'}{Dt}-\Gamma\frac{\partial\psi}{\partial x}=\kappa\nabla^2T'$ · *why:* (11.27) with its
     nonlinear advection kept and w = ∂ψ/∂x. · *plain:* rising fluid warms the layer, diffusion smooths it.
  3. *did:* Write advection as a Jacobian · *tex:* $\frac{D}{Dt}f=\frac{\partial f}{\partial t}+J(\psi,f),\ \ J(\psi,f)=\psi_xf_z-\psi_zf_x$ · *why:*
     u∂_x + w∂_z with u = −ψ_z, w = ψ_x. · *plain:* the nonlinearity is a Jacobian.
  4. *did:* Note J(ψ, ∇²ψ) vanishes here · *tex:* $\nabla^2\psi=-a^2\psi,\ a^2=\pi^2+k^2\ \Rightarrow\ J(\psi,\nabla^2\psi)=-a^2J(\psi,\psi)=0$ · *why:* the
     single roll is an eigenfunction of ∇²; J(ψ, ψ) = 0. · *plain:* the roll does not advect its own vorticity.
  5. *did:* Make it dimensionless · *tex:* $\frac1\Pr\partial_t\nabla^2\psi=\mathrm{Ra}\,\theta_x+\nabla^4\psi,\ \ \partial_t\theta+J(\psi,\theta)-\psi_x=\nabla^2\theta$ · *why:*
     lengths d, time d²/κ, ψ by κ, T′ by ΔT (θ); Ra of (11.21) appears as in D07. · *plain:* two equations, two numbers (Ra, Pr).
  6. *did:* Insert the truncation · *tex:* $\psi=A(t)\cos\pi z\,\sin kx,\ \ \theta=B(t)\cos\pi z\cos kx+C(t)\sin2\pi z$ · *why:* (11.90) with amplitudes
     A, B, C; each shape meets the free–free, isothermal conditions at z = ±½ (C05). · *plain:* three time-dependent amplitudes.
  7. *did:* Project the vorticity equation · *tex:* $\dot A=\Pr\Big(\frac{\mathrm{Ra}\,k}{a^2}B-a^2A\Big)$ · *why:* ∇²ψ = −a²ψ, ∇⁴ψ = a⁴ψ, θ_x =
     −kB cos πz sin kx (the C term has no x); match the coefficients of cos πz sin kx. · *plain:* buoyancy B spins up the roll,
     viscosity slows it.
  8. *did:* Compute the Jacobian J(ψ, θ) · *tex:* $J(\psi,\theta)=-\frac{\pi k}2AB\sin2\pi z+2\pi kAC\cos\pi z\cos2\pi z\cos kx$ · *why:* product rule
     on each factor; sin² + cos² = 1 collects the AB terms; 2 sin πz cos πz = sin 2πz. · *plain:* the roll moves heat between the
     shapes.
  9. *did:* Project on cos πz cos kx · *tex:* $\dot B=-a^2B+kA-\pi kAC$ · *why:* cos πz cos 2πz = ½(cos πz + cos 3πz) (P151); keep cos πz, drop
     the cos 3πz harmonic (the Galerkin truncation, P275); −ψ_x gives +kA. · *plain:* the roll makes the temperature difference B,
     and its interaction with C limits it.
  10. *did:* Project on sin 2πz · *tex:* $\dot C=\frac{\pi k}2AB-4\pi^2C$ · *why:* the x-averaged part of J feeds the mean-profile mode C;
      ∇²(sin 2πz) = −4π² sin 2πz. Every other harmonic is dropped. · *plain:* the roll flattens the temperature profile.
  11. *did:* Rescale time · *tex:* $\tau=a^2t:\ \ \frac{dA}{d\tau}=\Pr\Big(\frac{\mathrm{Ra}\,k}{a^4}B-A\Big),\ \frac{dB}{d\tau}=-B+\frac k{a^2}A-\frac{\pi k}{a^2}AC,\ \frac{dC}{d\tau}=
      \frac{\pi k}{2a^2}AB-\frac{4\pi^2}{a^2}C$ · *why:* divide each equation by a² so the decay of B has rate 1. · *plain:* time in units
      of d²/((π² + k²)κ).
  12. *did:* Rescale the amplitudes · *tex:* $X=\alpha_LA,\ Y=\beta_LB,\ Z=\gamma_LC$ · *why:* free scale factors chosen next so that every
      coefficient but r, b, Pr becomes 1. · *plain:* measure each mode in its natural unit.
  13. *did:* Fix the factors · *tex:* $\beta_L=\alpha_L\frac{\mathrm{Ra}\,k}{a^4},\ \ r\equiv\frac{\beta_Lk}{\alpha_La^2}=\frac{\mathrm{Ra}\,k^2}{a^6},\ \ \gamma_L=\pi r,\ \ \alpha_L=\frac{\pi k}{\sqrt2\,a^2}$
      · *why:* demand dX/dτ = Pr(Y − X), the rX term in dY/dτ, the −XZ term (α_Lγ_L = β_Lπk/a²) and the +XY term in dZ/dτ
      (α_Lβ_L = γ_Lπk/(2a²)). · *plain:* four conditions, four unknowns (r is the leftover number).
  14. *did:* Recognise r · *tex:* $r=\frac{\mathrm{Ra}\,k^2}{(\pi^2+k^2)^3}=\frac{\mathrm{Ra}}{\mathrm{Ra}_c(k)}$ · *why:* (11.44) with n = 1:
      Ra_c(k) = (π² + k²)³/k². · *plain:* r measures how far above its own onset the roll is heated.
  15. *did:* Write the system · *tex:* $\dot X=\Pr(Y-X),\ \ \dot Y=-XZ+rX-Y,\ \ \dot Z=XY-bZ,\ \ b=\frac{4\pi^2}{\pi^2+k^2}$ (11.91) · *why:* steps
      11–14; the last coefficient is 4π²/a². At the critical roll k² = π²/2, b = 8/3. · *plain:* Lorenz's equations.
- **Result.** (11.91) with r = Ra/Ra_cr and b = 4π²/(π² + k²) — *in words:* three numbers describe the roll speed (X), the temperature
  contrast (Y) and the flattening of the mean profile (Z).
- **Check.** r = 1 is the onset (C05) ✓; b(π/√2) = 8/3 ✓; conduction X = Y = Z = 0 is a solution ✓; `ch11.lorenz_sympy()` residuals 0.
- **sympy check intent (`check_src`).**
  ```python
  import sympy as sp                                              # symbolic algebra
  x, z, t, k, Pr, Ra = sp.symbols("x z t k Pr Ra", positive=True)   # coordinates, wavenumber, parameters
  A, B, C = [sp.Function(n)(t) for n in "ABC"]                    # the three amplitudes
  psi = A*sp.cos(sp.pi*z)*sp.sin(k*x)                             # (11.90), roll
  th = B*sp.cos(sp.pi*z)*sp.cos(k*x) + C*sp.sin(2*sp.pi*z)        # (11.90), temperature
  J = lambda f, g: sp.diff(f, x)*sp.diff(g, z) - sp.diff(f, z)*sp.diff(g, x)   # step 3
  lap = lambda f: sp.diff(f, x, 2) + sp.diff(f, z, 2)             # Laplacian
  heat = sp.diff(th, t) + J(psi, th) - sp.diff(psi, x) - lap(th)  # step 5, temperature residual
  proj = lambda e, s: sp.integrate(sp.integrate(sp.expand(e*s), (x, 0, 2*sp.pi/k)), (z, -sp.Rational(1, 2), sp.Rational(1, 2)))   # Galerkin projection
  eqB = sp.solve(proj(heat, sp.cos(sp.pi*z)*sp.cos(k*x)), sp.diff(B, t))[0]    # step 9
  eqC = sp.solve(proj(heat, sp.sin(2*sp.pi*z)), sp.diff(C, t))[0]               # step 10
  a2 = sp.pi**2 + k**2
  print(sp.simplify(eqB - (-a2*B + k*A - sp.pi*k*A*C)))           # -> 0
  print(sp.simplify(eqC - (sp.pi*k/2*A*B - 4*sp.pi**2*C)))        # -> 0
  print(sp.simplify(4*sp.pi**2/a2).subs(k, sp.pi/sp.sqrt(2)))     # b at k^2 = pi^2/2 -> 8/3
  ```
- **What it means.** The simplest possible model of convection that keeps the nonlinear heat transport (Z) — enough to produce
  chaos (D25). Fails for strong convection (many modes matter); the real Bénard layer is not chaotic at r = 28.
- **Traps.** The sign convention u = −∂ψ/∂z (slip #9); keeping cos 3πz (not in the truncation); forgetting that X, Y, Z are rescaled.

### D25 · The Lorenz fixed points, the pitchfork at r = 1 and the Hopf point $r_H=\Pr(\Pr+b+3)/(\Pr-b-1)$ — ★★, 11 steps, in C15 (notebook · `lorenz_attractor`)
- **Goal.** Find the steady states of (11.91) and the value of r at which the steady rolls lose stability — the book says only "if r
  is large".
- **Start.** $\dot X=\Pr(Y-X),\ \dot Y=-XZ+rX-Y,\ \dot Z=XY-bZ$ (11.91) — *in words:* Lorenz's system.
- **Plan.** (1) Set the right sides to zero. (2) Linearise with the Jacobian. (3) Origin: a 2 × 2 block. (4) Rolls: a cubic, and the
  condition for a pair of imaginary roots.
- **Tools.** Jacobian and fixed-point stability (P274); purely imaginary roots of a cubic (P276); P214.
- **Assumptions.** Pr > b + 1 (step 11).
- **Steps.**
  1. *did:* Set Ẋ = 0 · *tex:* $Y=X$ · *why:* Pr ≠ 0. · *plain:* in a steady roll the two amplitudes match.
  2. *did:* Set Ż = 0 · *tex:* $Z=\frac{XY}b=\frac{X^2}b$ · *why:* b > 0; use step 1. · *plain:* the profile flattens with the roll's strength.
  3. *did:* Set Ẏ = 0 · *tex:* $X\Big(r-1-\frac{X^2}b\Big)=0$ · *why:* −XZ + rX − Y with steps 1–2. · *plain:* either no roll, or a roll of a
     definite size.
  4. *did:* List the fixed points · *tex:* $(0,0,0);\ \ C_\pm=\big(\pm\sqrt{b(r-1)},\ \pm\sqrt{b(r-1)},\ r-1\big)\ (r>1)$ · *why:* step 3's two
     factors. · *plain:* conduction, and two rolls turning opposite ways.
  5. *did:* Write the Jacobian · *tex:* $J=\begin{pmatrix}-\Pr&\Pr&0\\r-Z&-1&-X\\Y&X&-b\end{pmatrix}$ · *why:* partial derivatives of the three right sides
     (P274). · *plain:* the linear response near any state.
  6. *did:* Linearise at the origin · *tex:* $\lambda=-b,\ \ \lambda^2+(\Pr+1)\lambda+\Pr(1-r)=0$ · *why:* at X = Y = Z = 0 the Z-row decouples; the 2 × 2
     block gives the quadratic. · *plain:* the conduction state's growth rates.
  7. *did:* Read the pitchfork · *tex:* $r>1\ \Rightarrow\ \lambda_1\lambda_2=\Pr(1-r)<0$ · *why:* a negative product means one positive root; the
     origin becomes unstable at r = 1, exactly where C± appear (Ra = Ra_c). · *plain:* convection starts at r = 1.
  8. *did:* Linearise at C± · *tex:* $J_C=\begin{pmatrix}-\Pr&\Pr&0\\1&-1&-\bar X\\\bar X&\bar X&-b\end{pmatrix}$ · *why:* r − Z̄ = 1 and X̄ = Ȳ. · *plain:*
     the rolls' response matrix.
  9. *did:* Compute the characteristic cubic · *tex:* $\lambda^3+(\Pr+b+1)\lambda^2+b(r+\Pr)\lambda+2b\Pr(r-1)=0$ · *why:* expand det(λI − J_C) along
     the first row; X̄² = b(r − 1). · *plain:* three growth rates of the rolls.
  10. *did:* Impose a pair of imaginary roots · *tex:* $(\Pr+b+1)\,b(r+\Pr)=2b\Pr(r-1)$ · *why:* λ = iω is a root iff a₂a₁ = a₀ (P276). · *plain:*
      the rolls' spirals neither shrink nor grow.
  11. *did:* Solve for r · *tex:* $r_H=\frac{\Pr(\Pr+b+3)}{\Pr-b-1}=24.74\ (\Pr=10,\ b=8/3)$ · *why:* collect r: r(b + 1 − Pr) = −Pr(Pr + b + 3);
      positive only if Pr > b + 1. Above r_H, a₂a₁ < a₀ and a pair has positive real part (P276). · *plain:* beyond 24.74 no steady
      state is stable.
- **Result.** Fixed points (0, 0, 0) and C±; the origin unstable for r > 1; C± unstable for r > r_H = Pr(Pr + b + 3)/(Pr − b − 1) =
  24.74 — *in words:* above r_H the motion has nowhere to settle.
- **Check.** r = 28: eig(C±) = −13.85, 0.094 ± 10.19i (positive real part) ✓; r = 20: −13.36, −0.155 ± 8.71i (stable) ✓; at r_H
  the real part crosses zero ✓ (`ch11.lorenz_eigs`, `lorenz_hopf_r`).
- **What it means.** Together with volume contraction −(Pr + 1 + b), the motion is confined to a zero-volume set yet cannot settle:
  a strange attractor. The Hopf bifurcation is subcritical (no small stable cycle is born — the chaos is not a gentle oscillation).
- **Traps.** Expanding the 3 × 3 determinant with a sign error in the cofactors; forgetting r − Z̄ = 1; applying a₂a₁ = a₀ with the
  inequality reversed.
