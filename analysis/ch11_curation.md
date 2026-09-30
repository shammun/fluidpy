# Chapter 11 — Instability: curation
(from `analysis/ch11.md` — 170 inventory rows (#1–#170), 96 equation numbers (11.1)–(11.96), 37 derivations in §2b (a-D1…a-D37),
the analyst's proposed A spine of 16 (§3 notes); `book.yaml` ch11 (14 sections, 5 viz seeds) and `policy.tier_a_full_treatment: [12, 18]`,
`coverage: exhaustive`. Curated 2026-09-30, concept-curator. Read: `knowledge/CUMULATIVE.md`, `concept_map.md`, `primers.md` (P01–P254 +
P218a; **ch11 primers start at P255**), `notation.md`, `viz_patterns.md`, `analysis/ch10_curation.md` (format), the `interactive-viz` skill §4–§5.)

Counts: A 15 · B 131 · C 24 · RECAP 24 · SKIP 2 · derivations written out 25 (★ 7 ★★ 13 ★★★ 5) · demoted to statements 11

Tier words (parsed by `tools/nbkit.py`): **CORE 15 · NOTE 129 · RECAP 24 · SKIP 2** = 170 rows; **DERIVATION 25**.
Depth: A 15 (all CORE) · B 131 (110 NOTE + 21 RECAP) · C 24 (19 NOTE + 3 RECAP: R01, R02, R07 + 2 SKIP: S01, S02).
**Reconciliation with analysis §2:** the 170 rows `#1…#170` each appear exactly once in §2 below (last column). A = 15 of 170 rows (8.8 %).

**Decisions that shape this chapter.**
1. **Fifteen A items from the analyst's sixteen.** Merges: (a) the vortex-sheet limit (#30, SEEN in ch05/ch07) is a RECAP inside the
   KH block **C02** — its new content (growth kΔU/2 at every k) is one line of (11.18), and its roll-up animation reuses ch05; (b) the
   Rayleigh number (#38) is a B item opening the Bénard block **C03**, whose A row is the amplitude problem (11.36)–(11.37) — Ra is
   *defined and explained* there (with the Γ sign callout) but its only derivation is the scaling argument; (c) Rayleigh's circulation
   criterion (#79, SEEN in ch08 as a flag) is a RECAP inside the Taylor block **C07** with the ring-interchange energy (#78) as its new
   reason; (d) Squire's theorem (#119) is a B item inside the Orr–Sommerfeld block **C11** (stated with the transformation and a numerical
   3-D ↔ 2-D check — it does not hold for the rotating/stratified flows of ch13, which is said there); (e) Fjørtoft (#127) is a B item
   of the inflection-point block **C12**; (f) Table 11.1, Blasius and the tanh/Bickley layers are B items of the plane-Poiseuille block
   **C13**; (g) the Lorenz fixed points, Hopf point and attractor (#155), attractors and bifurcations (#152) and period doubling (#156)
   are B items of the Lorenz block **C15**. Next candidates if a reviewer asks for one more A (all B now): Squire (#119), the Rayleigh
   equation (#124), Table 11.1 (#139).
2. **Bénard gets three A blocks** (C03 amplitude/growth problem, C04 rigid–rigid neutral curve and Ra_c, C05 free–free closed form)
   because each carries a different method that later sections reuse: the Chebyshev eigen-solver (C03 → C07, C08, C11), the
   "neutral curve + minimise over k" recipe (C04 → C07, C13) and the hand-solvable modes (C05 → C06's 27π⁴/4 and Lorenz's b = 8/3).
3. **SEEN rows are RECAP** (twenty-four, R01–R24, inventory order): KH's potential-flow machinery (ch04–ch07), the Boussinesq set and the
   conduction base state, circular Couette and Rayleigh's criterion (ch08), cylindrical and dimensionless Navier–Stokes (ch04), N², ψ and
   Ri (ch01/ch04/ch07), Bickley and Falkner–Skan (ch09). Each reuses the earlier chapter's function; where the row adds something new
   (vortex-sheet growth at every k, the inflection rule now as a stability statement) that part is stated.
4. **SKIP is minimal** — S01 (Exercises 11.3–11.5, porous/compliant/membrane), S02 (literature). The analyst's prose SKIP rows hold ideas
   and are B or C items with a pointer (history, photographs, the Marangoni remark, nonlinear effects, transition, closure on
   predictability). Exercises the lesson uses are B items: Ex. 11.1 (KH with depth and surface tension — the explainer seed), Ex. 11.2
   (Rayleigh–Taylor cut-off, our answer), Ex. 11.6 (σ real), Ex. 11.7 (odd mode), Ex. 11.9 (slip #7 and the Galerkin route), Ex. 11.11
   (piecewise shear layer as a V1 check), Ex. 11.12 ((11.95) parity), Ex. 11.13 ((11.96), the start of D23), Ex. 11.14 (logistic map).
   Exercise inputs and answers printed by the book stay private.
5. **Derivations (D01–D25).** 26 of the 37 analysis derivations are written out, merged into 24 D rows (a-D2 + a-D3 in D03, a-D8 + a-D9
   in D07), plus one of ours (D12, the free–free growth-rate quadratic that drives the "growth rate vs wavenumber" explainer). Written out
   because the book never writes them: D08 (σ real, Ex. 11.6), D12, D13 ((11.45)), D14 (narrow-gap Taylor, ★★★), D21 (Orr–Sommerfeld),
   D23 ((11.88)), D24 (Lorenz from (11.90), ★★★), D25 (Hopf point r_H), D04 (Ex. 11.1). Eleven are demoted to statements (§4c). The five
   ★★★ (D10 determinant, D14 Taylor narrow gap, D18 Miles–Howard, D19 Howard semicircle, D24 Lorenz) are each shown in an explainer.
6. **Printed slips (analysis §9 S1–S12) become named callouts "slip #k"** in the notebook ("book prints X, correct is Y"; table in §8);
   `S1…S12` stay the code keys of `ch11.book_slips()`. Planted wrong variants a test must fail: slip #1 (B for C in the determinant),
   #2 (sin nπz), #3 (factor 3 in dRa/dK²), #4 (∂p/∂x), #6 (Tollmien middle branch), #7 (squared operator), #11 (−k²(U − c)F).
7. **Climate hooks named where they are real:** C02/C09 (KH billows in the thermocline and billow clouds; Ri < ¼ in every mixing
   parameterisation), C06 (salt fingers and thermohaline staircases in the subtropical thermocline; the diffusive regime under polar
   ice), C03–C05 (Rayleigh–Bénard ↔ shallow atmospheric convection, cloud streets as rolls), C07 (centrifugal ↔ inertial instability in
   ch13), C12 (Rayleigh–Kuo barotropic instability of jets in ch13 is this theorem with U″ → U″ − β), C15 (Lorenz's predictability limit
   for weather; ensemble forecasting), N110 (rotating annulus → baroclinic waves, ch13 §13.17).
8. **Book values stay private** (`tests/book_values_ch11.json`): the book's rounded critical numbers, Taylor's radius ratio, the sech²-jet
   Ri, most-unstable wavelength ratio, pipe and suction numbers, Tollmien coefficients, low-Pr experiment numbers, exercise inputs.
   Published benchmarks are cited to their authors (Chandrasekhar 1961 Ra_c = 1707.76; Orszag 1971 Re_c = 5772.22; Thomas / Jordinson
   for Blasius; Tatsumi & Kakutani 1958; Michalke 1964; Lorenz 1963; Feigenbaum 1978) and stored in `reference/ch11/`.

## 1. Teaching order (A IDs grouped by book section, B/C IDs under each; one sentence each: "once you see X, Y follows")
Order = book order. Prerequisite notes: N40 (σ real, exchange of stabilities) closes C03 before the marginal problem (C04); C05's
growth-rate quadratic (D12) uses C05's sine modes, so it comes after them; the Howard semicircle (C10) is taught before the
Orr–Sommerfeld block (C11) as in the book and is reused as the eigenvalue window of C11–C13. RECAP rows are reminders inside the block
named as their A parent.

**§11.1 Introduction** — no A item of its own: N01–N03 open the C01 block (basic state and disturbance, the four potential wells of
Fig. 11.1 animated, temporal vs spatial instability, no Coriolis until ch13).

**§11.2 Method of Normal Modes**
- **C01** The normal mode (11.1) and the stability vocabulary (σ_r, c_i, stable/neutral/unstable, marginal, stationary vs oscillatory): once you write any small disturbance as a sum of e^{ikx + σt} pieces, linearity lets each piece evolve alone, so "is the flow stable?" becomes "is Re σ(k) ≤ 0 for every k?" — and the critical point is where the σ_r(k) curve first touches zero.
  - B: N01, N02, N04, N05, N06
  - C: N03   ·   RECAP: –

**§11.3 Kelvin–Helmholtz Instability**
- **C02** The Kelvin–Helmholtz dispersion relation (11.18) and the stability boundary g(ρ₂² − ρ₁²) < kρ₁ρ₂(U₂ − U₁)²: once you linearise the two interface conditions and insert decaying modes, the wave speed solves a real quadratic whose discriminant is "gravity's restoring term minus shear's kinetic term" — negative discriminant means a complex-conjugate pair, one of which grows; short waves always lose, surface tension and finite depth move the boundary.
  - B: N07, N08, N09, N10, N11, N12, N13, N14, N15, N16, N17, N18, N19, N22, N120, N121
  - C: N20, N21, N23   ·   RECAP: R01, R02, R03, R04, R05, R06, R07, R08, R09, R10, R11

**§11.4 Thermal Instability: The Bénard Problem**
- **C03** The Bénard amplitude problem (11.36)–(11.37) and its growth rate σ(K, Ra): once you split the Boussinesq set into conduction + perturbation, eliminate the pressure and insert normal modes, two coupled ODEs in z with one parameter Ra = gαΓd⁴/(κν) remain, σ is an eigenvalue (computed by Chebyshev collocation), and σ is real (exchange of stabilities), so onset is σ = 0.
  - B: N25, N26, N27, N28, N29, N30, N31, N32, N33, N34, N35, N36, N37, N38, N39, N40, N122
  - C: N24   ·   RECAP: R12, R13
- **C04** The rigid–rigid neutral curve and the critical Rayleigh number Ra_c ≈ 1707.76 at K_c ≈ 3.117 (Fig. 11.10): once σ = 0 the pair collapses into one sixth-order equation whose even solution is a cos + two cosh terms, three boundary conditions give a 3 × 3 determinant, each K has its own marginal Ra(K), and the bottom of that valley is the onset (cells about as wide as the layer is deep).
  - B: N41, N42, N43, N44, N45, N46, N50, N51, N123
  - C: N52   ·   RECAP: –
- **C05** Free–free boundaries: Ra = (π² + K²)³/K² (11.44), K_c² = π²/2, Ra_c = 27π⁴/4: once stress-free walls make every even derivative of W vanish, sine modes solve the problem exactly, the neutral curve is a formula you can minimise by hand, and its growth rate is a quadratic in σ — the closed form that double diffusion and Lorenz reuse.
  - B: N47, N48, N49
  - C: –   ·   RECAP: –

**§11.5 Double-Diffusive Instability**
- **C06** The salt-finger criterion (11.46): once salt is added, the marginal problem is identical to free–free Bénard with Ra replaced by Rs − Ra, so a column that is denser at the bottom can still overturn when the salt Rayleigh number (divided by the slow salt diffusivity) beats the heat one by 27π⁴/4.
  - B: N53, N54, N55, N56
  - C: N57   ·   RECAP: –

**§11.6 Centrifugal Instability: Taylor Problem**
- **C07** Taylor–Couette onset: the narrow-gap equations (11.51), the Taylor number (11.52) and Ta_c ≈ 1708/(½(1 + Ω₂/Ω₁)) (11.54): once you see that swapping two fluid rings releases energy when Γ² falls outward (Rayleigh), the centrifugal force plays gravity, the narrow-gap equations reduce to Bénard's when the cylinders co-rotate, and viscosity adds the same 1708-type threshold.
  - B: N59, N60, N61, N62, N63, N64, N65, N66, N67, N68, N124
  - C: N58, N69   ·   RECAP: R14, R15, R16

**§11.7 Instability of Continuously Stratified Parallel Flows**
- **C08** The Taylor–Goldstein equation (11.61): once you linearise the inviscid Boussinesq equations about U(z), ρ̄(z), use ψ and normal modes, one second-order ODE in z with the eigenvalue c remains; with N² = 0 it is Rayleigh's equation, and real coefficients make unstable modes come in growing/decaying pairs.
  - B: N71, N72, N73, N74
  - C: N70   ·   RECAP: R17, R18, R19
- **C09** The Miles–Howard criterion Ri > ¼ everywhere ⇒ stable (11.67): once you substitute φ = ψ̂/(U − c)^{1/2}, multiply by φ*, integrate and take the imaginary part, c_i times a sum of positive integrals equals c_i times ∫(N² − ¼U′²)/∣U − c∣²∣φ∣² — impossible for c_i ≠ 0 when Ri > ¼ everywhere; ¼ is a guarantee of stability, not a trigger of instability.
  - B: N75, N76, N77, N78
  - C: –   ·   RECAP: R20
- **C10** Howard's semicircle (11.72): once you substitute F = ψ̂/(U − c), the real and imaginary parts of one integral identity confine every unstable c to the half-disc on [U_min, U_max] — c_r inside the velocity range, growth rate kc_i ≤ k(U_max − U_min)/2.
  - B: N79, N80, N81, N82
  - C: –   ·   RECAP: –

**§11.8 Squire's Theorem and the Orr–Sommerfeld Equation**
- **C11** The Orr–Sommerfeld equation (11.79) with no-slip (11.80): once Squire's transformation shows that every unstable 3-D mode has a more unstable 2-D twin at lower Re, a stream function and normal modes turn the linearised Navier–Stokes equations into one fourth-order eigenvalue problem for c(k, Re) — the tool behind every viscous result of §11.10–11.11.
  - B: N84, N85, N86, N87, N88, N89, N90
  - C: N83   ·   RECAP: R21, R22

**§11.9 Inviscid Stability of Parallel Flows**
- **C12** Rayleigh's inflection-point theorem (11.83)–(11.84), with Fjørtoft's sharpening (11.85)–(11.86): once the viscous term is dropped (Rayleigh's equation (11.81)), the imaginary part of one integral identity says c_i∫U″∣φ∣²/∣U − c∣² = 0, so a growing mode needs U″ to change sign — an inflection point — and Fjørtoft adds that the vorticity must peak there; neutral modes carry a critical layer and Kelvin's cat's eyes.
  - B: N91, N92, N93, N94, N95, N96, N97, N126, N127
  - C: –   ·   RECAP: –

**§11.10 Results for Parallel and Nearly Parallel Viscous Flows**
- **C13** Plane Poiseuille flow: inviscidly stable, viscously unstable at Re_c = 5772.22, k_c = 1.02056 (Orszag 1971), and the Table 11.1 family: once you trace c_i(k, Re) = 0 with the Orr–Sommerfeld solver, a thumb-shaped neutral curve appears even though U = 1 − y² has no inflection point — the Tollmien–Schlichting wave; inflectional profiles (tanh, Bickley jet) go unstable at tiny Re, Couette and pipe never do linearly.
  - B: N99, N100, N102, N105, N107, N108
  - C: N98, N101, N103, N106, N109   ·   RECAP: R23, R24
- **C14** The disturbance kinetic-energy equation (11.88): once you multiply the perturbation momentum equation by u_i and average over a wavelength, dE/dt = production −∫uv U′ minus dissipation Λ — viscosity destabilises by shifting the phase of v against u so that the Reynolds stress −⟨uv⟩ can draw energy from the mean shear; the same −⟨uv⟩U′ is ch12's turbulence production, and nonlinear Reynolds stresses are what saturate the growth (§11.12) and start transition (§11.13).
  - B: N104, N110, N111, N128
  - C: N125   ·   RECAP: –

**§11.11 Experimental Verification of Boundary-Layer Instability** — no A item of its own: N105 (Tollmien profile, slip #6), N106
(Schubauer & Skramstad), N107 (Blasius neutral curve in frequency), N108 (Blasius Re_δ*,c ≈ 519) and N109 (non-parallel corrections,
non-normality) sit in the C13 block.

**§11.12 Comments on Nonlinear Effects** — no A item of its own: N110 (rectified fluxes modify the basic state and arrest growth;
annulus → ch13) closes the C14 block.

**§11.13 Transition** — no A item of its own: N111 (sequence of instabilities, secondary 3-D instability of TS waves, roll-up and pairing
in free shear layers) sits in the C14 block with a pointer to ch12.

**§11.14 Deterministic Chaos**
- **C15** The Lorenz system (11.91) from a three-mode truncation of free–free convection (11.90): once you keep one roll mode and two temperature modes, Bénard convection becomes three nonlinear ODEs whose fixed points (conduction and two convection states) all lose stability past r_H = Pr(Pr + b + 3)/(Pr − b − 1) ≈ 24.74, and the trajectory wanders forever on a strange attractor where nearby states separate exponentially — deterministic but unpredictable.
  - B: N112, N113, N114, N115, N116, N117, N129
  - C: N118, N119   ·   RECAP: –

**End matter** — S01 (Exercises 11.3–11.5), S02 (literature cited).

## 2. Chapter map (depth) — every inventory row exactly once
Rows are in analysis §2 order (inventory row number in the last column). IDs: `C` = A/CORE (teaching order), `N` = B or C NOTE (inventory
order), `R` = RECAP, `S` = SKIP. No ASCII pipe inside cells (∣K∣ is written with U+2223). "slip #k" = analysis §9 Sk.

| ID | Item | § | Depth | Tier | A parent | Reason (A) / treatment (B) / pointer (C, SKIP) / source chapter (RECAP) | Row |
|---|---|---|---|---|---|---|---|
| N01 | Def: basic (background) state, infinitesimal disturbance; linear stability vs finite-amplitude (nonlinear) instability | 11.1 | B | NOTE | C01 | stated as the opening paragraph of C01 with the parcel picture of ch01 (N²) and ch09's perturb–linearise–eigenvalues (P214) | #1 |
| N02 | Fig. 11.1 four mechanical analogies: bowl (stable), cap (unstable), plane (neutral), dimple (stable to small, unstable to large kicks) | 11.1 | B | NOTE | C01 | stated with our animation A1 (`potential_well_demo`): a damped ball in four potentials, small vs large kick | #2 |
| N03 | Aims (onset of transition); temporal vs spatial instability (Huerre & Monkewitz); successes of linear theory; no Coriolis here | 11.1 | C | NOTE | C01 | named; pointer: rotation and baroclinic instability in Ch. 13 §13.17; spatial growth named only | #3 |
| C01 | Eq. (11.1) normal mode û(z) exp{ikx + imy + σt} = û(z) exp{i∣K∣(e_K·x − ct)}, σ = −i∣K∣c | 11.2 | A | CORE | – | load-bearing: every section of the chapter is this one recipe (perturb with e^{ikx+σt}, solve for σ or c); ch13's baroclinic growth rates use it | #4 |
| N04 | Def: stable / neutrally stable / unstable as σ_r (c_i) < 0, = 0, > 0 for every wavenumber | 11.2 | B | NOTE | C01 | stated with `ST.stability_class` and the "for every k" quantifier (logic primer); one number per class | #5 |
| N05 | Normal-mode method: any disturbance = superposition of modes; linearity ⇒ modes evolve independently ⇒ eigenvalue problem | 11.2 | B | NOTE | C01 | stated with ch10's Fourier error modes (G ↔ e^{σΔt}) as the bridge; superposition shown on an interface with three modes | #6 |
| N06 | Def: marginal state σ_r = 0; stationary (σ_i = 0, cells) vs oscillatory onset (overstability); neutral ≠ marginal | 11.2 | B | NOTE | C01 | stated with `ST.marginal_type`; a two-row table (stationary: Bénard, Taylor; oscillatory: diffusive double diffusion) | #7 |
| N07 | Def: Kelvin–Helmholtz set-up: streams U₁, ρ₁ above and U₂, ρ₂ below, interface z = ζ(x, t); inertial ideal flow (Fig. 11.2) | 11.3 | B | NOTE | C02 | stated with our sketch from `kh_fields` (interface and velocity arrows) | #8 |
| R01 | Eq. (11.2) total potentials φ̃ = Ux + φ; perturbation irrotational by Kelvin's theorem | 11.3 | C | RECAP | C02 | ch06 velocity potential; ch05 Kelvin's theorem (5.8) | #9 |
| R02 | Eq. (11.3) Laplace's equation for the perturbation potentials | 11.3 | C | RECAP | C02 | ch06 / ch07 §7.1 (∇²φ = 0) | #10 |
| R03 | Eq. (11.4), Eq. (11.5) decay far above and below | 11.3 | B | RECAP | C02 | ch07 §7.7 interface waves (e^{∓kz}); the sign choice is a step of D03 | #11 |
| R04 | Eq. (11.6), Eq. (11.7) kinematic and dynamic interface conditions | 11.3 | B | RECAP | C02 | ch07 (7.14), (7.20); ch04 §4.10 kinematic condition; the starting line of D02 | #12 |
| R05 | Eq. (11.8) kinematic condition with n = ∇f/∣∇f∣, f = z − ζ | 11.3 | B | RECAP | C02 | ch07 (7.16), ch03 level sets (P75); first moves of D02 | #13 |
| N08 | Exact kinematic condition after the dot products (unnumbered, p. 478) | 11.3 | B | NOTE | C02 | stated as a middle line of D02 (why the square root cancels: it multiplies all three members) | #14 |
| N09 | Eq. (11.9) linearised kinematic condition on z = 0 | 11.3 | B | NOTE | C02 | derived in D02 (A chain); residual checked by `kh_residuals` | #15 |
| R06 | Eq. (11.10) unsteady Bernoulli in each layer | 11.3 | B | RECAP | C02 | ch04 (4.75) (D26), `core.bernoulli.unsteady_bernoulli` | #16 |
| N10 | Eq. (11.11) pressure matching on z = ζ | 11.3 | B | NOTE | C02 | stated as a step of D02 | #17 |
| N11 | Eq. (11.12) undisturbed pressure balance (fixes C₁, C₂) | 11.3 | B | NOTE | C02 | stated as a step of D02 | #18 |
| N12 | Eq. (11.13) linearised dynamic condition | 11.3 | B | NOTE | C02 | derived in D02 (½∣∇φ̃∣² − ½U² → U∂φ/∂x); `kh_residuals` | #19 |
| N13 | Eq. (11.14) normal-mode trial φ_j = A_j(z) e^{ik(x − ct)} | 11.3 | B | NOTE | C02 | stated as the first move of D03 (C01's mode with K = (k, 0, 0)) | #20 |
| R07 | Amplitude ODE A″ = k²A and its exponentials | 11.3 | C | RECAP | C02 | ch07 separation of variables (P167), linear ODE by e^{mz} trial (P44) | #21 |
| R08 | Eq. (11.15) decaying solutions φ₁ = A₋e^{ik(x−ct)−kz}, φ₂ = A₊e^{ik(x−ct)+kz} | 11.3 | B | RECAP | C02 | ch07 §7.7 `interface_fields`; reused in `kh_fields` | #22 |
| N14 | Eq. (11.16), Eq. (11.17) remnants of the interface conditions | 11.3 | B | NOTE | C02 | derived in D03; `kh_amplitudes` gives A₋ = −i(U₁ − c)ζₒ, A₊ = i(U₂ − c)ζₒ | #23 |
| N15 | Amplitudes and the quadratic ρ₁(U₁ − c)² + ρ₂(U₂ − c)² = (g/k)(ρ₂ − ρ₁) (p. 479) | 11.3 | B | NOTE | C02 | derived in D03; `kh_sympy` re-derives it | #24 |
| C02 | Eq. (11.18) KH dispersion relation c = (ρ₂U₂ + ρ₁U₁)/(ρ₂ + ρ₁) ± [((ρ₂ − ρ₁)/(ρ₂ + ρ₁))(g/k) − ρ₁ρ₂(U₂ − U₁)²/(ρ₂ + ρ₁)²]^{1/2} | 11.3 | A | CORE | – | load-bearing: the first complete normal-mode calculation and the prototype of shear instability (billows, wind waves); feeds the Richardson criterion (C09) and the vortex-sheet roll-up | #25 |
| N16 | Instability criterion g(ρ₂² − ρ₁²) < kρ₁ρ₂(U₂ − U₁)² (p. 480) | 11.3 | B | NOTE | C02 | derived in the last steps of D03; `kh_critical_k`; the stability-boundary curve ΔU_min(k) | #26 |
| N17 | Real coefficients ⇒ growing and decaying modes come in conjugate pairs | 11.3 | B | NOTE | C02 | stated (§4c a-D20) with the two roots c₊, c₋ = c₊* on the c-plane of E2 | #27 |
| R09 | Eq. (11.19) static-interface limit c = ±[((ρ₂ − ρ₁)/(ρ₂ + ρ₁))(g/k)]^{1/2}; Rayleigh–Taylor when ρ₁ > ρ₂ | 11.3 | B | RECAP | C02 | ch07 (7.95) `core.waves.interface_omega` (parity test); ⚠️ slip #5: the book points to (7.96) | #28 |
| N18 | Short waves always unstable when U₁ ≠ U₂ (k > k_c) | 11.3 | B | NOTE | C02 | stated with the growth curve kc_i(k) from `kh_growth_rate` (F1) and the surface-tension cure (N120) | #29 |
| R10 | Eq. (11.20) vortex-sheet limit ρ₁ = ρ₂: c = (U₁ + U₂)/2 ± i(U₂ − U₁)/2, growth kΔU/2 at every k | 11.3 | B | RECAP | C02 | ch05 §5.8 vortex sheet (γ = U₂ − U₁), ch07 interface vortex sheet; new: unstable at every wavelength, `vortex_sheet_c` | #30 |
| N19 | Fig. 11.3 frame moving with (U₁ + U₂)/2: symmetric, the wave is stationary (c_r = mean) | 11.3 | B | NOTE | C02 | stated with a frame toggle in E2 and ch03's Galilean shift (`core.kinematics.galilean_transform`) | #31 |
| N20 | Shear vs stratification; tilting-tube experiment (Thorpe), billow clouds, thermocline dye (Woods); shear instability as an internal-wave source | 11.3 | C | NOTE | C02 | named with the climate hook; pointer: Ri criterion in C09, internal waves ch07 §7.8, mixing ch12 | #32 |
| R11 | Fig. 11.6 nonlinear roll-up of a perturbed vortex sheet | 11.3 | B | RECAP | C02 | ch05 `sheet_rollup` (and explainer `vortex_sheet_rollup`); our animation A2 overlays the linear growth e^{kΔU t/2} | #33 |
| N21 | Energy source of KH = kinetic energy of the streams; Fig. 11.7 profiles (PE rises, KE falls) | 11.3 | C | NOTE | C02 | named; quantified in N22 | #34 |
| N22 | Mixing example: linear profile after mixing, E_final/E_initial = 2/3 with momentum unchanged (p. 482) | 11.3 | B | NOTE | C02 | stated (§4c a-D4) with `kh_mixing_energy` and bars before/after; one number per energy | #35 |
| N23 | General statement: fixed ∫U dz ⇒ ∫U² dz falls when gradients are smoothed (Cauchy–Schwarz) | 11.3 | C | NOTE | C02 | named; a random-smoothing property test backs it; pointer: energy budget C14 | #36 |
| N24 | Bénard history (Bénard 1900 hexagons were surface-tension/Marangoni driven; Rayleigh 1916; Jeffreys 1928) | 11.4 | C | NOTE | C03 | named, with the Marangoni remark; pointer: planforms N51 | #37 |
| N25 | Eq. (11.21) Rayleigh number Ra = gαΓd⁴/(κν), Γ = −dT̄/dz (positive when heated from below) | 11.4 | B | NOTE | C03 | stated with `rayleigh_number(g, alpha, dT, d, kappa, nu)` (dT = T_bottom − T_top), a water-layer number and the ⚠️ two-convention table (slip #10, §9) | #38 |
| R12 | Boussinesq set (4.10), (4.86), (4.89) with ρ = ρ₀[1 − α(T̃ − T₀)] | 11.4 | B | RECAP | C03 | ch04 C13 (D27, D28); the starting line of D05 | #39 |
| N26 | Fig. 11.8 geometry: layer of depth d, z centred, T̄ = T₀ − Γ(z + d/2) | 11.4 | B | NOTE | C03 | stated with our sketch (T̄(z) beside the layer) | #40 |
| N27 | Eq. (11.22) decomposition ũ = 0 + u, T̃ = T̄(z) + T′, p̃ = P(z) + p | 11.4 | B | NOTE | C03 | stated as the first move of D05 | #41 |
| R13 | Eq. (11.23) base state: hydrostatic P(z) and conduction ∂²T̄/∂z² = 0 | 11.4 | B | RECAP | C03 | ch01 hydrostatics, ch08 steady conduction (linear profile); `benard_base_state` | #42 |
| N28 | Eq. (11.24) conduction profile T̄ = T₀ − ½ΔT − Γz, Γ ≡ ΔT/d | 11.4 | B | NOTE | C03 | stated with `benard_base_state` and one number | #43 |
| N29 | Nonlinear perturbation equations with the −wΓ term (p. 485) | 11.4 | B | NOTE | C03 | a middle line of D05 ((u·∇)T̄ = w dT̄/dz = −wΓ) | #44 |
| N30 | Eq. (11.25), Eq. (11.26), Eq. (11.27) linearised perturbation equations | 11.4 | B | NOTE | C03 | derived in D05; sympy `benard_perturbation_sympy` | #45 |
| N31 | Scaling w ~ κ/d and buoyancy/viscous ratio ~ Ra (p. 485–486) | 11.4 | B | NOTE | C03 | stated (§4c a-D6) with numbers for a water and an air layer | #46 |
| N32 | Eq. (11.28) Laplacian of the z-momentum equation | 11.4 | B | NOTE | C03 | a step of D06 | #47 |
| N33 | Pressure Poisson equation and its z-derivative (p. 486) | 11.4 | B | NOTE | C03 | steps of D06 (recalls ch10's pressure Poisson, D19 of ch10) | #48 |
| N34 | Eq. (11.29) ∂_t∇²w = gα∇_H²T′ + ν∇⁴w | 11.4 | B | NOTE | C03 | derived in D06 (A chain); planted ∇² instead of ∇_H² fails the sympy test | #49 |
| N35 | Eq. (11.30) rigid isothermal walls w = ∂w/∂z = T′ = 0 | 11.4 | B | NOTE | C03 | stated (∂w/∂z = 0 from continuity at a no-slip wall) | #50 |
| N36 | Eq. (11.31), Eq. (11.32), Eq. (11.33) non-dimensional equations and Pr = ν/κ | 11.4 | B | NOTE | C03 | derived in D07; ⚠️ w stays dimensional until W (first of three scalings, §9) | #51 |
| N37 | Normal modes in x, y and operator substitutions ∂_t → σ, ∇_H² → −K², ∇² → d²/dz² − K² | 11.4 | B | NOTE | C03 | steps of D07 | #52 |
| N38 | Eq. (11.34), Eq. (11.35) amplitude equations for T̂ and ŵ | 11.4 | B | NOTE | C03 | derived in D07 | #53 |
| C03 | Eq. (11.36), Eq. (11.37) scaled amplitude problem in W ≡ (Γd²/κ)ŵ with the single parameter Ra; growth rate σ(K, Ra) | 11.4 | A | CORE | – | key method: the Bénard linear problem as an eigenvalue problem, solved by the Chebyshev collocation reused for Taylor, Taylor–Goldstein and Orr–Sommerfeld; Ra reappears in double diffusion and Lorenz | #54 |
| N39 | Eq. (11.38) W = ∂W/∂z = T̂ = 0 on z = ±½ | 11.4 | B | NOTE | C03 | stated with the boundary-row replacement (primer) in `benard_growth_rate` | #55 |
| N40 | σ real for Ra > 0 (exchange of stabilities) ⇒ marginal state σ = 0 with stationary cells | 11.4 | B | NOTE | C03 | derived in D08 (A chain; the book leaves it to Ex. 11.6); all computed σ real to 1e-10 | #56 |
| N41 | Eq. (11.39) marginal pair (d²/dz² − K²)T̂ = −W, (d²/dz² − K²)²W = RaK²T̂ | 11.4 | B | NOTE | C04 | derived in D09 (A chain); `benard_marginal_Ra` | #57 |
| N42 | Eq. (11.40) sixth-order equation (d²/dz² − K²)³W = −RaK²W | 11.4 | B | NOTE | C04 | derived in D09 | #58 |
| N43 | Eq. (11.41) rigid conditions in W only | 11.4 | B | NOTE | C04 | last step of D09 (third condition from T̂ = 0) | #59 |
| N44 | Pr drops out at the margin; even/odd modes about z = 0 (Fig. 11.9: one vs two rows of cells) | 11.4 | B | NOTE | C04 | stated with `benard_eigenfunction` streamlines for both parities (parity primer) | #60 |
| N45 | Characteristic equation (q² − K²)³ = −RaK² and roots (11.42), q₀ | 11.4 | B | NOTE | C04 | derived in D10 (A chain; cube roots of −1, primer); `benard_char_roots` | #61 |
| N46 | Even solution A cos q₀z + B cosh qz + C cosh q*z and the 3 × 3 determinant (p. 489) | 11.4 | B | NOTE | C04 | derived in D10; ⚠️ slip #1 (the derivative line prints B for C; planted variant moves the root by > 1 %) | #62 |
| C04 | Fig. 11.10 rigid–rigid neutral curve Ra(K) and critical point Ra_c = 1707.76, K_c = 3.117 (Chandrasekhar 1961) | 11.4 | A | CORE | – | load-bearing: the chapter's signature result and its recipe "neutral curve + minimise over the wavenumber", reused by Taylor (μ → 1 gives the same 1708) and by every viscous critical Reynolds number (C13) | #63 |
| N47 | Free (stress-free) surfaces ⇒ ∂²w/∂z² = 0; Eq. (11.43) W = W″ = W⁗ = 0 | 11.4 | B | NOTE | C05 | derived in D11 (A chain) | #64 |
| N48 | Eigenfunctions W = A sin nπ(z + ½) (slip #2: printed sin nπz) | 11.4 | B | NOTE | C05 | stated in D11; ⚠️ slip #2 callout, `benard_free_free_mode(printed=True)` fails W(±½) = 0 at n = 1 | #65 |
| C05 | Eq. (11.44) free–free eigenvalue relation Ra = (n²π² + K²)³/K² | 11.4 | A | CORE | – | load-bearing: the one Bénard case solvable by hand; its minimum 27π⁴/4 is the double-diffusive threshold (C06) and its K_c² = π²/2 gives Lorenz's b = 8/3 (C15) | #66 |
| N49 | dRa/dK² = 0 ⇒ K_c² = π²/2, Ra_c = 27π⁴/4 ≈ 657.5 (slip #3: spurious factor 3) | 11.4 | B | NOTE | C05 | derived in the last steps of D11; ⚠️ slip #3 callout (the printed form has no root); sympy `benard_free_free_sympy` | #67 |
| N50 | Rigid bottom, free top: Ra_c ≈ 1100.65 at K_c ≈ 2.682 (computed; Chandrasekhar 1961) | 11.4 | B | NOTE | C04 | stated with `benard_critical(bc=("rigid","free"))` and its curve on the C04 figure (book's rounding private) | #68 |
| N51 | Planform undetermined by linear theory: rolls, squares, hexagons (Figs. 11.11–11.12); rising centre in liquids | 11.4 | B | NOTE | C04 | stated with `planform` (∇_H²f = −K²f, Helmholtz primer) and three small panels; turbulence at large Ra named | #69 |
| N52 | Nonlinear equilibrium: KE generation balanced by dissipation, PE release | 11.4 | C | NOTE | C04 | named; pointer: energy budget C14, saturation N110, Lorenz C15 | #70 |
| N53 | Double diffusion: EOS ρ̃ = ρ₀[1 − α(T̃ − T₀) + β(s̃ − s₀)], κ_s ≪ κ; diffusion can destabilise a statically stable column | 11.5 | B | NOTE | C06 | stated with `linear_eos` (ch01 seawater coefficients recalled) and one density number | #71 |
| N54 | Two regimes (Fig. 11.13): salt fingers (hot salty over cold fresh, σ real) and the diffusive regime (cold fresh over hot salty, oscillatory, layering) | 11.5 | B | NOTE | C06 | stated with parcel arguments and our free–free cubic `double_diffusive_sigma` (three roots, regime label; the D12 moves with salt added) | #72 |
| N55 | Eq. (11.45) double-diffusive marginal equations, with Ra's sign flipped (dT̄/dz) and Rs′ | 11.5 | B | NOTE | C06 | derived in D13 (A chain; the book writes "repeat the derivation"); ⚠️ second Ra sign convention (§9) | #73 |
| N56 | Reduction T̂ = κ_sŝ/κ ⇒ (11.39) with Ra → Rs − Ra; free–free threshold Rs − Ra = 27π⁴/4 | 11.5 | B | NOTE | C06 | derived in D13 (uniqueness primer); `double_diffusive_margin` | #74 |
| C06 | Eq. (11.46) salt-finger criterion (gd⁴/ν)[(β/κ_s)dS/dz − (α/κ)dT̄/dz] = 27π⁴/4 | 11.5 | A | CORE | – | key physics: diffusion can destabilise a statically stable column — salt fingers and thermohaline staircases of the subtropical ocean; the cleanest use of C05 | #75 |
| N57 | Finger width ~ d at onset, thin fingers above it; staircases; layering by heating a salinity gradient (Fig. 11.15); the two requirements | 11.5 | C | NOTE | C06 | named with the climate hook (thermohaline staircases, Arctic diffusive layering); pointer: mixing in ch12 | #76 |
| N58 | Taylor problem set-up: rotating coaxial cylinders, axisymmetric disturbances ∂/∂θ = 0 | 11.6 | C | NOTE | C07 | named at the top of C07; pointer: wavy vortices N69 | #77 |
| N59 | Rayleigh ring interchange: Γ = 2πrU_θ conserved, E = Γ²/(8π²r²), ΔE = (Γ₂² − Γ₁²)(1/r₁² − 1/r₂²)/(8π²) (p. 496) | 11.6 | B | NOTE | C07 | stated (§4c a-D15) with `ring_interchange_energy` and bars for one ring pair (E5 inspector) | #78 |
| R14 | Rayleigh's circulation criterion: inviscidly unstable iff dΓ²/dr < 0 somewhere (analogy dρ̄/dz > 0) | 11.6 | B | RECAP | C07 | ch08 Rayleigh-stability flag Ω₂/Ω₁ > (R₁/R₂)² in `core.laminar`; new: the ring-interchange reason (N59) and `rayleigh_circulation_criterion` | #79 |
| R15 | Eq. (11.47) axisymmetric Navier–Stokes in (R, φ, z) | 11.6 | B | RECAP | C07 | ch04 cylindrical operators `core.curvilinear`, P186; ⚠️ slip #12 continuity needs 1/R | #80 |
| N60 | Eq. (11.48) decomposition ũ = U + u, p̃ = P + p | 11.6 | B | NOTE | C07 | the first move of D14 | #81 |
| R16 | Eq. (11.49) circular Couette base state U_φ = AR + B/R with A, B | 11.6 | B | RECAP | C07 | ch08 (8.9)–(8.10) `core.laminar.circular_couette(return_coeffs=True)` | #82 |
| N61 | Fig. 11.16 geometry: gap d, counter-rotating toroidal rolls | 11.6 | B | NOTE | C07 | stated with our cross-section sketch (E5 view 1) | #83 |
| N62 | Eq. (11.50) linearised perturbation equations, the (dU_φ/dR + U_φ/R)u_R term | 11.6 | B | NOTE | C07 | derived in D14 (book defers to Chandrasekhar); sympy `taylor_perturbation_sympy` | #84 |
| N63 | Normal modes exp{ikz + σt}; narrow gap d ≪ (R₁ + R₂)/2 | 11.6 | B | NOTE | C07 | steps of D14 (narrow-gap primer) | #85 |
| N64 | Eq. (11.51) narrow-gap Taylor equations with α ≡ Ω₂/Ω₁ − 1, x = (R − R₁)/d | 11.6 | B | NOTE | C07 | derived in D14 (A chain, ★★★); `taylor_marginal_Ta` (Chebyshev on x ∈ [0, 1]) | #86 |
| N65 | Eq. (11.52) Taylor number; narrow gap, inner only Ta = 2(Ω₁R₁d/ν)²(d/R₁) | 11.6 | B | NOTE | C07 | derived in D15; `taylor_number`, `taylor_number_narrow_inner` with one number | #87 |
| N66 | Eq. (11.53) û_R = dû_R/dR = û_φ = 0 at x = 0, 1 | 11.6 | B | NOTE | C07 | stated (no slip on both cylinders) | #88 |
| N67 | Marginal state σ = 0 assumed (proved for co-rotation only) | 11.6 | B | NOTE | C07 | stated; σ checked real for μ ≥ 0 by `taylor_growth_rate`; open for counter-rotation | #89 |
| N68 | Fig. 11.17 theory vs Taylor's experiment in the (Ω₂, Ω₁) plane; Rayleigh line | 11.6 | B | NOTE | C07 | stated with `taylor_stability_boundary` for our radius ratio (the book's ratio private) and the Rayleigh line | #90 |
| C07 | Eq. (11.54) critical Taylor number Ta_cr = 1708/(½(1 + Ω₂/Ω₁)); μ → 1 reduces to Bénard; counter-rotation Rayleigh-unstable but viscously stable | 11.6 | A | CORE | – | load-bearing: rotation's instability is buoyancy in disguise (Γ² plays −ρ), with the same 1708 threshold; the analogue of inertial instability in ch13 | #91 |
| N69 | Taylor vortices, wavy vortices (∂/∂φ ≠ 0), turbulence (Coles); Dean and Görtler secondary flows; non-uniqueness of NS solutions | 11.6 | C | NOTE | C07 | named; pointer: transition N111 and Ch. 12 | #92 |
| N70 | Stratified shear-flow history: Taylor 1915, Miles 1961, Howard 1961 | 11.7 | C | NOTE | C08 | named at the top of C08 | #93 |
| R17 | Set-up U(z)e_x, ρ̄(z); inviscid Boussinesq momentum and basic balance (p. 502–503) | 11.7 | B | RECAP | C08 | ch04 Boussinesq, ch07 §7.8 base state (`boussinesq_linear_sympy`); 2-D by Squire is an assumption here (N88) | #94 |
| N71 | Eq. (11.55) perturbation momentum (slip #4: the w-equation prints ∂p/∂x) | 11.7 | B | NOTE | C08 | derived in D16; ⚠️ slip #4 callout, `stratified_shear_sympy(printed=True)` fails | #95 |
| R18 | Linearised density equation and Eq. (11.56) N² ≡ −(g/ρ₀)dρ̄/dz | 11.7 | B | RECAP | C08 | ch01 (1.29), ch07 (7.128) `core.stratification.brunt_vaisala_sq` | #96 |
| R19 | Eq. (11.57) stream function u = ∂ψ/∂z, w = −∂ψ/∂x | 11.7 | B | RECAP | C08 | ch04 (4.12) ψ; ⚠️ slip #9 sign conventions change between §11.7, §11.8 and §11.14 | #97 |
| N72 | Normal modes [ρ, p, ψ] = [ρ̂, p̂, ψ̂]e^{ik(x−ct)} → Eq. (11.58), Eq. (11.59), Eq. (11.60) | 11.7 | B | NOTE | C08 | derived in D17 | #98 |
| C08 | Eq. (11.61) Taylor–Goldstein equation (U − c)(ψ̂″ − k²ψ̂) − U″ψ̂ + N²ψ̂/(U − c) = 0 | 11.7 | A | CORE | – | load-bearing: the master equation of inviscid stratified shear flow (contains Rayleigh's equation at N² = 0); source of C09, C10; stratified jets and thermocline billows in ch13 | #99 |
| N73 | Conjugate symmetry: ψ̂*, c* also solve ⇒ growing/decaying pairs | 11.7 | B | NOTE | C08 | stated (§4c a-D20); spectrum closed under conjugation (property test) | #100 |
| N74 | Eq. (11.62) rigid lids ψ̂(0) = ψ̂(d) = 0 | 11.7 | B | NOTE | C08 | stated with the BC rows of `taylor_goldstein_eigs` | #101 |
| N75 | Eq. (11.63) φ = ψ̂/(U − c)^{1/2} and the two derivative formulas (p. 504) | 11.7 | B | NOTE | C09 | derived as steps of D18 (complex powers P155) | #102 |
| N76 | Eq. (11.64) self-adjoint form | 11.7 | B | NOTE | C09 | derived in D18 ("after some rearrangement" written out); sympy check | #103 |
| N77 | Eq. (11.65) integral identity and its imaginary part | 11.7 | B | NOTE | C09 | derived in D18; residual of a computed TG mode by `richardson_identity_check` | #104 |
| R20 | Eq. (11.66) gradient Richardson number Ri(z) ≡ N²/(dU/dz)² | 11.7 | B | RECAP | C09 | ch04 notation (Ri_g), ch01 N²; `gradient_richardson` on the tanh/sech² family | #105 |
| C09 | Eq. (11.67) Miles–Howard criterion: Ri > ¼ everywhere ⇒ stable (Ri < ¼ somewhere is necessary for instability) | 11.7 | A | CORE | – | load-bearing: the ¼ threshold behind ocean and atmosphere mixing parameterisations (ch12, ch13); the first "multiply by the conjugate and integrate" proof | #106 |
| N78 | Not sufficient: no universal Ri_c; tanh/erf shear layers unstable iff Ri_min < ¼; most-amplified wavelength; Scotti & Corcos, Eriksen | 11.7 | B | NOTE | C09 | stated with the (k, J) growth map of `taylor_goldstein_eigs` for U = tanh z, N² = J sech²z (cached; tongue closes at J = ¼) | #107 |
| N79 | Eq. (11.68) F ≡ ψ̂/(U − c) and its derivatives | 11.7 | B | NOTE | C10 | steps of D19 | #108 |
| N80 | Transformed divergence form (p. 506; slip #11: printed −k²(U − c)F, must be (U − c)²) and the weight Q | 11.7 | B | NOTE | C10 | derived in D19; ⚠️ slip #11 callout (sympy decides) | #109 |
| N81 | Eq. (11.69), Eq. (11.70) real and imaginary parts | 11.7 | B | NOTE | C10 | derived in D19; residuals by `howard_identity_check` | #110 |
| N82 | Eq. (11.71) U_min < c_r < U_max for an unstable mode | 11.7 | B | NOTE | C10 | derived in D19; `in_howard_semicircle` | #111 |
| C10 | Eq. (11.72) Howard's semicircle [c_r − ½(U_max + U_min)]² + c_i² ≤ [½(U_max − U_min)]² and kc_i ≤ (k/2)(U_max − U_min) (Fig. 11.20) | 11.7 | A | CORE | – | load-bearing: bounds every unstable eigenvalue — the sanity window for every eigen-solver here and in ch13 (Eady/Charney–Stern semicircles) | #112 |
| N83 | Viscous parallel flows: viscosity can destabilise (intro to §11.8) | 11.8 | C | NOTE | C11 | named at the top of C11; pointer: C13 (Poiseuille) and C14 (why) | #113 |
| R21 | Eq. (11.73) non-dimensional perturbed x-momentum; Re = U₀L/ν | 11.8 | B | RECAP | C11 | ch04 D30 dimensionless NS, `core.similarity`; the starting line of D20 | #114 |
| N84 | Eq. (11.74) perturbation x-momentum | 11.8 | B | NOTE | C11 | derived in D20 | #115 |
| N85 | Eq. (11.75) y-, z-momentum and continuity | 11.8 | B | NOTE | C11 | derived in D20 | #116 |
| N86 | Eq. (11.76) 3-D normal modes exp{i(kx + mz − kct)} | 11.8 | B | NOTE | C11 | a step of D20 | #117 |
| N87 | Eq. (11.77) 3-D normal-mode equations | 11.8 | B | NOTE | C11 | derived in D20 (A chain); `os_3d_eigs` uses them only for the Squire check | #118 |
| N88 | Squire's theorem and transformation Eq. (11.78) k̄ = √(k² + m²), k̄Re̅ = kRe (2-D equations, p. 510) | 11.8 | B | NOTE | C11 | stated (§4c a-D24) with `ST.squire_transform` and the numerical check `os_3d_eigs(k, m, Re)` = `orr_sommerfeld_eigs(k̄, Re̅)`; ⚠️ not valid with rotation or stratification (ch13) | #119 |
| N89 | Consequences: equivalent 2-D problem at lower Re; larger 2-D growth; oblique wave sees only the flow along it | 11.8 | B | NOTE | C11 | stated with one number (oblique angle 30° → Re̅ = Re cos 30°) and the E8 oblique slider | #120 |
| R22 | 2-D stream function u = ∂ψ/∂y, v = −∂ψ/∂x; û = φ′, v̂ = −ikφ (φ is not a potential) | 11.8 | B | RECAP | C11 | ch04 (4.12) ψ; ⚠️ slip #9 sign (here y is the cross-stream coordinate) | #121 |
| C11 | Eq. (11.79) Orr–Sommerfeld equation (U − c)(φ″ − k²φ) − U″φ = (φ⁗ − 2k²φ″ + k⁴φ)/(ikRe) | 11.8 | A | CORE | – | key method: the viscous parallel-flow eigenproblem behind TS waves, Table 11.1 and boundary-layer transition (ch12, ch14); `ST.orr_sommerfeld_eigs` reused later | #122 |
| N90 | Eq. (11.80) no slip φ = φ′ = 0 at the walls | 11.8 | B | NOTE | C11 | stated (clamped BC rows) | #123 |
| N91 | Eq. (11.81) Rayleigh equation (Re → ∞) | 11.9 | B | NOTE | C12 | stated (§4c a-D26) as the singular limit (order 4 → 2, φ′ = 0 dropped); `ST.rayleigh_eigs` | #124 |
| N92 | Eq. (11.82) inviscid conditions φ = 0; c(k) eigenvalue; conjugate pairs, broken by the viscous term | 11.9 | B | NOTE | C12 | stated with a spectrum comparison (Rayleigh pairs vs OS no pairs) | #125 |
| C12 | Eq. (11.83), Eq. (11.84) Rayleigh's inflection-point theorem: c_i∫U″∣φ∣²/∣U − c∣² dy = 0 ⇒ U″ changes sign | 11.9 | A | CORE | – | load-bearing: classifies profiles without solving anything; ch13's barotropic (Rayleigh–Kuo) instability is this theorem with U″ − β | #126 |
| N93 | Fjørtoft's theorem: Eq. (11.85), Eq. (11.86) ⇒ (U − U_I)U″ < 0 somewhere (vorticity maximum at the inflection point) | 11.9 | B | NOTE | C12 | stated (§4c a-D28) as "add (c_r − U_I) times (11.84) to the real part"; `fjortoft_criterion`; ⚠️ the book calls U_I "U₁" | #127 |
| N94 | Fig. 11.21 six profiles; only the jet-like and shear-layer ones pass both criteria | 11.9 | B | NOTE | C12 | stated as a verdict table from `rayleigh_criterion`/`fjortoft_criterion` on our analytic stand-ins | #128 |
| N95 | Criteria are not sufficient: U = sin y on ∣y∣ ≤ b stable for 2b < π | 11.9 | B | NOTE | C12 | stated with a `rayleigh_eigs` sweep in b (max c_i → 0 as 2b → π⁺) | #129 |
| N96 | Critical layer U(y_c) = c_r; critical point of (11.81); viscous critical layer thins as Re → ∞ | 11.9 | B | NOTE | C12 | stated with `critical_layer` (singular-point primer) | #130 |
| N97 | Kelvin cat's eye: Eq. (11.87) and its expansion ψ̂ ≈ ½(y − y_c)²U′(y_c) + Aφ(y_c)cos kx (Fig. 11.22) | 11.9 | B | NOTE | C12 | stated (§4c a-D29) with `cats_eye_streamfunction`, streamlines and the eye width 2√(Aφ_c/U′_c) | #131 |
| N98 | Viscosity both stabilises and destabilises; OS asymptotics hard, numerics used | 11.10 | C | NOTE | C13 | named; pointer: C14 explains the destabilising mechanism | #132 |
| N99 | Two-stream shear layer U₀ tanh(y/L): Re_c = 0, neutral curve k_u(Re) → 1/L, sech neutral mode at k = 1/L (Fig. 11.23) | 11.10 | B | NOTE | C13 | stated with `tanh_shear_layer_neutral_curve` (cached) and the analytic neutral mode (sympy residual 0) | #133 |
| R23 | Bickley jet U₀ sech²(y/L): small Re_c ≈ 4 at kL ≈ 0.2 (Tatsumi & Kakutani 1958), sinuous mode | 11.10 | B | RECAP | C13 | ch09 (9.71) `core.jets.free_jet_profile`; new: `bickley_critical`, inviscid neutral modes k = 2, 1 with c = 2/3 | #134 |
| C13 | Plane Poiseuille flow: inviscidly stable, viscously unstable (TS waves); Re_c = 5772.22, k_c = 1.02056 (Orszag 1971); singular-perturbation wall and critical layers | 11.10 | A | CORE | – | load-bearing: "viscosity destabilises" — the central result of §11.10–11.11 and the benchmark that proves the OS solver (Re = 10⁴, k = 1: c = 0.23752649 + 0.00373967i) | #135 |
| N100 | Plane Couette flow linearly stable at every Re | 11.10 | B | NOTE | C13 | stated with `couette_max_growth` on a (k, Re) grid (always < 0); finite-amplitude transition named | #136 |
| N101 | Pipe flow linearly stable; transition by finite disturbances, inlet effects; chaotic saddle (Eckhardt et al. 2007) | 11.10 | C | NOTE | C13 | named (numbers private); pointer: Ch. 12 transition and turbulence | #137 |
| R24 | Boundary layers with pressure gradient (Fig. 11.24): favourable loops close, adverse has an inflection and a flat upper branch | 11.10 | B | RECAP | C13 | ch09 (9.51)–(9.52) inflection rule, `core.boundary_layer.falkner_skan`; new: `falkner_skan_neutral_curve` (cached) | #138 |
| N102 | Table 11.1 critical Reynolds numbers (jet, shear layer, Blasius, plane Poiseuille, pipe, plane Couette) | 11.10 | B | NOTE | C13 | stated as our recomputed table `table_11_1()` beside the published benchmarks (book's rounded column private) | #139 |
| N103 | Finite Re_c for a spreading mixing layer (Bhattacharya et al. 2006) | 11.10 | C | NOTE | C13 | named in one sentence (parallel-flow assumption vs spreading) | #140 |
| C14 | Eq. (11.88) disturbance kinetic-energy equation d/dt∫½u_i² dV = −∫u_iu_j ∂U_i/∂x_j dV − Λ (Fig. 11.25) | 11.10 | A | CORE | – | load-bearing: explains why viscosity can destabilise (Reynolds-stress production) and seeds ch12's production term −⟨uv⟩U′ | #141 |
| N104 | 2-D form and interpretation: production −∫uv U′ dV by the Reynolds stress; out-of-phase u, v produce nothing; viscosity shifts the phase | 11.10 | B | NOTE | C14 | stated with `disturbance_energy_budget` for a TS mode: Reynolds-stress profile and the u–v phase across the channel (P178 recalled) | #142 |
| N105 | Tollmien–Schlichting approximate Blasius profile (slip #6: the printed middle branch is discontinuous; 1 − b(1 − η)² is continuous) | 11.11 | B | NOTE | C13 | stated with `tollmien_profile` vs Blasius; ⚠️ slip #6 (`printed=True` fails continuity; coefficients private) | #143 |
| N106 | TS waves undetected until Schubauer & Skramstad (1947): vibrating ribbon, hot wires | 11.11 | C | NOTE | C13 | named (history) | #144 |
| N107 | Fig. 11.26 Blasius neutral curve in frequency ων/U∞² vs Re_δ* | 11.11 | B | NOTE | C13 | stated with `blasius_neutral_curve(in_frequency=True)` (cached) and F = kc_r/Re_δ* | #145 |
| N108 | Blasius parallel-flow OS: Re_δ*,c ≈ 519.2 at αδ* ≈ 0.303 (Thomas; Jordinson 520) | 11.11 | B | NOTE | C13 | stated with `blasius_critical` (cached; U″ = −½ff″ from ch09, no numerical second derivative) | #146 |
| N109 | Non-parallel corrections (Nayfeh & Saric), transient growth and non-normality (Reshotko), suction delays transition | 11.11 | C | NOTE | C13 | named (fractions private); pointer: non-normal growth in Ch. 12 | #147 |
| N110 | Nonlinear effects: rectified fluxes ⟨uv⟩, ⟨uT′⟩ change the basic state and arrest growth; finite-amplitude destabilisation; rotating annulus → chaos | 11.12 | B | NOTE | C14 | stated as the paragraph after the budget (the production term changes U itself); pointer: Ch. 13 §13.17 annulus, Ch. 12 | #148 |
| N111 | Transition: sequence of instabilities (Landau), secondary 3-D instability of TS waves (Klebanoff, Fig. 11.27 peak–valley), roll-up and pairing in free shear layers | 11.13 | B | NOTE | C14 | stated as a five-box ASCII pipeline (TS wave → 3-D → Λ-vortices → spots → turbulence; shear layer: roll-up → pairing); pointer: Ch. 12 | #149 |
| N112 | Deterministic chaos defined: sensitivity, aperiodicity, broadband spectrum; steady forcing can give aperiodic response | 11.14 | B | NOTE | C15 | stated with the C15 two-trajectory figure and a spectrum (ch10 `dominant_frequency` recalled) | #150 |
| N113 | Pendulum Ẍ + (g/l)sin X = 0 as Eq. (11.89); phase space, trajectory, degrees of freedom | 11.14 | B | NOTE | C15 | stated with `pendulum_rhs`, `phase_portrait` (energy conserved to 1e-9) | #151 |
| N114 | Attractors: fixed point, limit cycle, repeller; bifurcation at R_cr; Fig. 11.28 | 11.14 | B | NOTE | C15 | stated with the Hopf normal form (ours, `hopf_normal_form`, amplitude √μ) and three panels | #152 |
| N115 | Eq. (11.90) Lorenz truncation: ψ ∝ X cos πz sin kx, T′ ∝ Y cos πz cos kx + Z sin 2πz | 11.14 | B | NOTE | C15 | the starting line of D24; `lorenz_fields` drives the roll animation; ⚠️ slip #9 sign (u = −∂ψ/∂z here) | #153 |
| C15 | Eq. (11.91) Lorenz system Ẋ = Pr(Y − X), Ẏ = −XZ + rX − Y, Ż = XY − bZ with r = Ra/Ra_c, b = 4π²/(π² + k²) | 11.14 | A | CORE | – | load-bearing: deterministic chaos from three modes of convection — the predictability limit behind weather forecasting and climate ensembles | #154 |
| N116 | Steady states (origin; X = Y = ±√(b(r − 1)), Z = r − 1); instability at large r (Hopf r_H ≈ 24.74, ours); Lorenz's parameters; aperiodic switching (Fig. 11.29); strange attractor (Fig. 11.30) | 11.14 | B | NOTE | C15 | stated with `lorenz_fixed_points`, `lorenz_jacobian`, `lorenz_hopf_r`; derivation D25 (A chain; the book says only "if r is large") | #155 |
| N117 | Routes to chaos (1): period doubling, Feigenbaum ratio → 4.6692; Fig. 11.31 tree | 11.14 | B | NOTE | C15 | stated (§4c a-D34) with `bifurcation_diagram`, `period_doubling_points` and a δ table (low-Pr experiment numbers private) | #156 |
| N118 | Routes to chaos (2): quasi-periodic (Ruelle–Takens), Landau's conjecture; Lorenz follows neither | 11.14 | C | NOTE | C15 | named; pointer: turbulence onset debate in Ch. 12 | #157 |
| N119 | Closure: unpredictability (weather), deterministic ≠ quantum uncertainty; Poincaré; relevance to turbulence unclear | 11.14 | C | NOTE | C15 | named with the climate hook (predictability ~ two weeks, ensembles); pointer: Ch. 13 | #158 |
| N120 | Ex. 11.1 KH with a lower layer of depth h and surface tension (formula for c and the instability inequality) | Ex. 11.1 | B | NOTE | C02 | derived in D04 (A chain; the book gives only the result); `kh_phase_speed(surface_tension, h)`, `kh_min_shear` (air over water ΔU_min ≈ 6.7 m/s at λ ≈ 1.7 cm) | #159 |
| N121 | Ex. 11.2 Rayleigh–Taylor (heavy over light) with surface tension: longest neutral wavelength λ_c = 2π√(σ_s/((ρ_h − ρ_l)g)) | Ex. 11.2 | B | NOTE | C02 | stated (ours; the book prints no answer) with `rayleigh_taylor_cutoff` and one number (water over air ≈ 1.7 cm) | #160 |
| S01 | Ex. 11.3–11.5 porous surface, compliant surface, membrane flutter | Ex. | C | SKIP | – | pointer: one line ("three exercises extend C02's interface conditions"); not solved | #161 |
| N122 | Ex. 11.6 σ real for Bénard (I₁, I₂, J₁, J₂ positive definite) | Ex. 11.6 | B | NOTE | C03 | the exercise is D08's construction (`exchange_of_stabilities_sympy`) | #162 |
| N123 | Ex. 11.7 gravest odd mode Ra ≈ 17610.39 at K ≈ 5.365 (Chandrasekhar 1961) | Ex. 11.7 | B | NOTE | C04 | stated with `benard_critical(mode="odd")` and the odd determinant (computed) | #163 |
| N124 | Ex. 11.8/11.9 narrow-gap algebra; Eq. (11.92), Eq. (11.93) (slip #7: squared operator), Eq. (11.94) Galerkin sine series | Ex. 11.9 | B | NOTE | C07 | stated with `taylor_galerkin_Ta` as the second route (agrees to 0.1 % with 4 modes); ⚠️ slip #7 callout | #164 |
| N125 | Ex. 11.10 stratified KH energy equation (kinetic + available potential energy) | Ex. 11.10 | C | NOTE | C14 | named with the equation in one line; pointer: available potential energy in Ch. 13 | #165 |
| N126 | Ex. 11.11 piecewise-linear shear layer: c₀² = ((U₁ − U₃)/(2kh))²{(kh − 1)² − e^{−2kh}}, neutral kh ≈ 1.2785 | Ex. 11.11 | B | NOTE | C12 | stated (§4c a-D36) with `piecewise_shear_layer_c` as a V1 target for `rayleigh_eigs` | #166 |
| N127 | Ex. 11.12 Eq. (11.95) Rayleigh equation for v̂; antisymmetric U ⇒ ±c; symmetric U ⇒ sinuous/varicose parity | Ex. 11.12 | B | NOTE | C12 | stated (§4c a-D37) with `rayleigh_eigs(parity=...)` on the Bickley jet | #167 |
| N128 | Ex. 11.13 Eq. (11.96) disturbed Navier–Stokes in index form | Ex. 11.13 | B | NOTE | C14 | the starting line of D23 | #168 |
| N129 | Ex. 11.14 logistic map x_{n+1} = Ax_n(1 − x_n) as a transition toy; background state x = 1 − 1/A | Ex. 11.14 | B | NOTE | C15 | stated (§4c a-D34) with `logistic_map` and a cobweb (exercise inputs private) | #169 |
| S02 | Literature cited | – | C | SKIP | – | pointer: one bibliography line; benchmarks cited where used (Chandrasekhar, Orszag, Thomas, Jordinson, Tatsumi & Kakutani, Michalke, Lorenz, Feigenbaum) | #170 |

## 3. Section coverage

| § | Title | A | B | C | RECAP | SKIP |
|---|---|---|---|---|---|---|
| 11.1 | Introduction | – (covered by B/C items opening the C01 block) | N01, N02 | N03 | – | – |
| 11.2 | Method of Normal Modes | C01 | N04, N05, N06 | – | – | – |
| 11.3 | Kelvin-Helmholtz Instability | C02 | N07, N08, N09, N10, N11, N12, N13, N14, N15, N16, N17, N18, N19, N22, N120, N121 | N20, N21, N23 | R01, R02, R03, R04, R05, R06, R07, R08, R09, R10, R11 | S01 |
| 11.4 | Thermal Instability: The Bénard Problem | C03, C04, C05 | N25, N26, N27, N28, N29, N30, N31, N32, N33, N34, N35, N36, N37, N38, N39, N40, N41, N42, N43, N44, N45, N46, N47, N48, N49, N50, N51, N122, N123 | N24, N52 | R12, R13 | – |
| 11.5 | Double-Diffusive Instability | C06 | N53, N54, N55, N56 | N57 | – | – |
| 11.6 | Centrifugal Instability: Taylor Problem | C07 | N59, N60, N61, N62, N63, N64, N65, N66, N67, N68, N124 | N58, N69 | R14, R15, R16 | – |
| 11.7 | Instability of Continuously Stratified Parallel Flows | C08, C09, C10 | N71, N72, N73, N74, N75, N76, N77, N78, N79, N80, N81, N82 | N70 | R17, R18, R19, R20 | – |
| 11.8 | Squire's Theorem and the Orr-Sommerfeld Equation | C11 | N84, N85, N86, N87, N88, N89, N90 | N83 | R21, R22 | – |
| 11.9 | Inviscid Stability of Parallel Flows | C12 | N91, N92, N93, N94, N95, N96, N97, N126, N127 | – | – | – |
| 11.10 | Results for Parallel and Nearly Parallel Viscous Flows | C13, C14 | N99, N100, N102, N104, N128 | N98, N101, N103, N125 | R23, R24 | – |
| 11.11 | Experimental Verification of Boundary-Layer Instability | – (covered by B/C items inside C13) | N105, N107, N108 | N106, N109 | – | – |
| 11.12 | Comments on Nonlinear Effects | – (covered by a B item closing C14) | N110 | – | – | – |
| 11.13 | Transition | – (covered by a B item inside C14) | N111 | – | – | – |
| 11.14 | Deterministic Chaos | C15 | N112, N113, N114, N115, N116, N117, N129 | N118, N119 | – | S02 |

Notes. Exercise rows are listed with the section where they are taught: Ex. 11.1–11.5 with §11.3 (N120, N121, S01), Ex. 11.6–11.7 with
§11.4 (N122, N123), Ex. 11.8–11.9 with §11.6 (N124), Ex. 11.10 and 11.13 with §11.10 (N125, N128), Ex. 11.11–11.12 with §11.9
(N126, N127), Ex. 11.14 with §11.14 (N129). S02 closes the chapter and is listed with §11.14.

## 4. Prerequisites needing primers (concept or tool | needed by | why it is not A/B/RECAP)
Numbering continues at **P255** (the designer assigns numbers). Already primed and only *reminded* (one sentence, no new primer): P13
log–log, P15 `assert np.allclose`, P16 animate, P17 slider_figure, P18 show_viz, P25 partial derivative, P26/P98 Taylor, P31/P94
`solve_ivp`, P38 product rule, P40 sympy, P41/P64 plotly 3-D, P44 linear ODE by e^{mx} trial, P45 Euler's formula, P49 chain rule, P53/P56
determinants and `np.linalg.det`, P58 null space, P75 level sets, P78 contour/quiver/streamplot, P80 eigenvalues, P81 complex conjugate,
P95 RK4 by hand, P108 `brentq`, P117 sympy series/collect, P129 completing the square, P130 order-of-magnitude scaling, P133 scaled
variables, P142 Fourier modes, P151 integrals of sines over a period, P153 the complex plane in numpy, P155 complex powers and branch
cuts, P159 complex square roots and the quadratic formula, P162 collocation, P167 separation of variables, P168 cosh/sinh/tanh, P170
`minimize_scalar`, P176 complex amplitudes, P177 operator elimination, P178 mean of a product of real parts, P186 cylindrical Laplacian,
P188 anisotropic scaling, P203 `simpson`/`trapezoid`, P212 inflection point, P214 linear stability by eigenvalues, P215 sech, P218a
integration by parts, P233 arbitrary coefficients ⇒ every bracket zero, P247 generalised symmetric eigenproblem, P249 Jacobian
determinant, P252 caching expensive runs; ch10 C07 Galerkin (FE) and C04 amplification factor G ↔ e^{σΔt}.

| Concept or tool | Needed by | Why it is not A/B/RECAP |
|---|---|---|
| Necessary vs sufficient conditions, "for every k", and proof by contradiction | C01 (N04), C09, C10, C12 | new logic vocabulary; the theorems of §11.7–11.9 are all of the form "if unstable then …" |
| Real and imaginary parts of a complex identity (one complex equation = two real ones; Im{1/(U − c)} = c_i/∣U − c∣² via the conjugate; ∫∣f∣² > 0 unless f ≡ 0) | C09 (D18), C10 (D19), C12 (D22), C03 (D08) | the central move of every integral-identity proof; extends P81, P153 |
| Cube roots of a negative number (−1 = e^{iπ}: roots −1, e^{±iπ/3}) | C04 (D10) | new; P45, P155 recalled |
| Even and odd functions; parity under z → −z | C04 (N44, D10), C13 (R23 sinuous/varicose), C12 (N127) | new |
| Eigenvalue problem for a differential equation (a parameter value for which a nonzero solution meets every boundary condition) | C03, C04, C07, C08, C11 | new; extends P80 from matrices to ODEs |
| Chebyshev–Gauss–Lobatto points and the differentiation matrix (spectral collocation; error falls exponentially in N) | C03, C04 (second route), C07, C08, C11, C13 | new numerical tool; P162 recalled |
| Boundary-row replacement and the generalised non-symmetric eigenproblem `scipy.linalg.eig(A, B)` (infinite eigenvalues, spurious modes, keep eigenvalues that do not move when N grows) | C03, C07, C08, C11, C13 | extends P247 (symmetric) to non-symmetric with singular B |
| Quadratic eigenvalue problem c²M₂ + cM₁ + M₀ and its companion linearisation | C08 | new |
| Mapping an infinite domain to [−1, 1] (algebraic map) and checking truncation by y_max | C11, C13 (tanh, Bickley, Blasius) | new |
| Integrating on a Chebyshev grid (Clenshaw–Curtis weights) | C14, C09/C10 identity checks | new; P203 recalled |
| Neutral curve as the zero contour of the growth rate; critical point as its minimum (`brentq` in Re inside `minimize_scalar` in k) | C04, C07, C13 | combines P108 and P170 into one recipe; new |
| Quotient rule, and differentiating with respect to K² as the variable | C05 (D11, slip #3) | new (only glossed in ch04 D22) |
| Uniqueness of a linear boundary-value problem (same equations + same conditions ⇒ same solution) | C06 (D13) | new |
| Narrow-gap (small-curvature) approximation: expand in d/R and keep the leading order | C07 (D14, D15) | new; P188 recalled |
| Completing the square into a circle (x − a)² + y² ≤ R² | C10 (D19) | extends P129 (tensors) to the scalar circle |
| Singular point of an ODE (the highest-derivative coefficient vanishes: U = c) | C08 (division by U − c), C12 (N96) | new |
| The integral of an x-derivative of a periodic function over one period is zero | C14 (D23) | new move; P151 recalled |
| Jacobian matrix of a nonlinear ODE system and the stability of its fixed points | C15 (D25, N113, N114) | extends P214 (stability by eigenvalues) and P249 (Jacobian determinant) |
| Purely imaginary roots of a cubic λ³ + a₂λ² + a₁λ + a₀: exactly when a₂a₁ = a₀ (the Hopf condition) | C15 (D25) | new |
| Galerkin truncation: keep a few modes and project with orthogonality of sines | C15 (D24) | extends ch10's FE Galerkin (C07) and P151, P233 to a modal basis |
| Lyapunov exponent: the slope of log(separation) against time | C15 | new |
| Iterated maps: fixed point, stability ∣f′(x*)∣ < 1, cobweb diagram | C15 (N117, N129) | new |
| Helmholtz equation in the plane ∇_H²f = −K²f; planforms as sums of cosines | C04 (N51) | new |
| `np.roots` (polynomial roots) and `np.lib.scimath.sqrt` (principal complex root of a negative real) | C02, C06 (N54), C15 | new Python tools; P159 recalled |
| matplotlib 3-D line plots (`fig.add_subplot(projection="3d")`) | C15 | new Python; P41/P64 are plotly |

## 4b. Derivations written out (parsed by tools: ID first, CORE id in a column, ★★★ for hard, explainer slugs backticked in the LAST column)
Analysis §2b numbers appear as (a-DNN) in the Result column. Every D row stands for one `nb.derivation` block.
| ID | Result (Eq.) | CORE | Difficulty | Steps | Tools used | Traps | Shown in |
|---|---|---|---|---|---|---|---|
| D01 | the two forms of the normal mode (11.1): σ = −i∣K∣c, so σ_r = ∣K∣c_i and σ_i = −∣K∣c_r (a-D33) | C01 | ★ | 5 | Euler's formula (P45), complex exponentials, equating exponents | the minus sign in σ_i: a mode moving to +x has σ_i < 0; only the real part of the field is physical | notebook · `normal_mode_growth` |
| D02 | linearised kinematic (11.9) and dynamic (11.13) interface conditions from (11.6)–(11.8), (11.10)–(11.12) (a-D1) | C02 | ★★ | 12 | level sets and normals (P75), Taylor transfer to z = 0 (P166), unsteady Bernoulli (R06), orders of smallness (P68) | the 1/√(1 + ζ_x²) factor multiplies all three members and cancels (book: "remove the square root"); ½∣∇φ̃∣² − ½U² → U∂φ/∂x, not ½(∂φ/∂x)²; the O(ζ) shift to z = 0 is quadratic | notebook · `kelvin_helmholtz_boundary` |
| D03 | KH dispersion relation (11.18) from (11.14)–(11.17), its instability inequality g(ρ₂² − ρ₁²) < kρ₁ρ₂(U₂ − U₁)², the limits (11.19), (11.20) and c_r = mean (a-D2 + a-D3) | C02 | ★★ | 14 | linear ODE by e^{mz} trial (P44), complex algebra (P153), quadratic formula (P159), inequalities (P48), sympy (P40) | decaying signs: e^{−kz} above, e^{+kz} below; (ik)² = −k²; the discriminant simplifies via (ρ₁U₁ + ρ₂U₂)² − (ρ₁ + ρ₂)(ρ₁U₁² + ρ₂U₂²) = −ρ₁ρ₂(U₁ − U₂)²; divide by ζₒ ≠ 0 and keep k > 0 | notebook · `kelvin_helmholtz_boundary` |
| D04 | KH with lower depth h and surface tension (Ex. 11.1): coth kh factors and the σ_s k term; minimum shear ΔU_min² = 2√(gΔρσ_s)(ρ₁ + ρ₂)/(ρ₁ρ₂) at k = √(gΔρ/σ_s) (a-D35 + ours; the book gives only the result) | C02 | ★★ | 9 | D03 moves, cosh/sinh (P168), Laplace pressure jump (ch01/ch04), minimisation of a/k + bk | the lower mode is cosh k(z + h), not e^{kz}; the pressure jump is −σ_sζ_xx = +σ_sk²ζ with the sign of the curvature; h → ∞ must give (11.18) plus σ_s k | notebook · `kelvin_helmholtz_boundary` |
| D05 | linear perturbation equations (11.25)–(11.27) from the Boussinesq set and (11.22)–(11.24) (a-D5) | C03 | ★ | 8 | Boussinesq set (R12), linearisation (drop products, P68) | (u·∇)T̄ = w dT̄/dz = −wΓ (sign from Γ = −dT̄/dz); the base pressure gradient cancels the reference buoyancy exactly, leaving +gαT′ e_z | notebook |
| D06 | pressure elimination (11.28) → (11.29): ∂_t∇²w = gα∇_H²T′ + ν∇⁴w (a-D7) | C03 | ★★ | 8 | Laplacian and divergence (ch02), commuting constant-coefficient operators (P121), ∇² − ∂_z² = ∇_H² | the divergence of (11.26) kills ∂_t and ν terms only because ∇·u = 0; forgetting the z-derivative of the Poisson equation leaves ∇²T′ instead of ∇_H²T′ (planted variant) | notebook |
| D07 | non-dimensional equations (11.31)–(11.33), normal modes and the scaled amplitude problem (11.34)–(11.38) (a-D8 + a-D9) | C03 | ★ | 9 | scaled variables and the chain rule (P133), operator substitution on e^{i(kx+ly)} (P177) | w is left dimensional on purpose; each ∂ brings 1/d and ∂_t brings κ/d²; Ra appears only after W ≡ (Γd²/κ)ŵ | notebook |
| D08 | σ is real for Ra > 0 (exchange of stabilities, Ex. 11.6; the book gives only the outline) (a-D13) | C03 | ★★ | 12 | integration by parts with the wall conditions (P218a), real/imaginary parts (primer), positive-definite integrals (primer) | integrate by parts twice for I₂, four times for J₂ and use every boundary condition; conjugate one relation before combining so the mixed integral cancels; the conclusion needs Ra > 0 | notebook |
| D09 | marginal pair (11.39) → sixth-order equation (11.40) with rigid conditions (11.41) (a-D10) | C04 | ★ | 6 | operator algebra (P177) | apply (d²/dz² − K²) to the second equation, then use the first; the third boundary condition comes from T̂ = 0 through the second equation | notebook · `benard_neutral_curve` |
| D10 | characteristic roots (11.42), even solution, the 3 × 3 determinant and Ra_c by minimising over K (a-D11; slip #1 corrected) | C04 | ★★★ | 14 | e^{qz} trial (P44), cube roots of −1 (primer), even functions (primer), determinants and nontrivial solutions (P53, P58), `brentq` (P108), `minimize_scalar` (P170), sympy (P40) | the cosh terms use the conjugate pair q, q* so W is real; the derivative line prints B for C (slip #1); the determinant is purely imaginary — use its imaginary part; spurious sign changes near Ra = K⁴ | notebook · `benard_neutral_curve` |
| D11 | free–free conditions (11.43), sine modes, (11.44), and dRa/dK² = 0 ⇒ K_c² = π²/2, Ra_c = 27π⁴/4 (a-D12; slips #2, #3 corrected) | C05 | ★★ | 11 | zero tangential stress, differentiating continuity, Fourier sine modes, quotient rule (primer), sympy (P40) | expanding (d²/dz² − K²)²W with W = W″ = 0 gives W⁗ = 0; the family is sin nπ(z + ½), not sin nπz (slip #2); no factor 3 on the second term of the derivative (slip #3); (3π²/2)³/(π²/2) = 27π⁴/4 | notebook · `benard_neutral_curve` |
| D12 | free–free growth rate: (σ + a²)(σ/Pr + a²)a² = RaK² with a² = π² + K², its two real roots and σ = 0 on (11.44) (ours; the book never writes it) | C05 | ★★ | 8 | sine modes from D11, operator substitution (P177), quadratic formula (P159) | (d²/dz² − K²) on sin nπ(z + ½) gives −a²; the discriminant is positive (both roots real — exchange of stabilities again); σ = 0 reproduces (11.44) | notebook · `normal_mode_growth` |
| D13 | double-diffusive marginal equations (11.45) and the finger criterion (11.46) via T̂ = κ_sŝ/κ and Rs − Ra = 27π⁴/4 (a-D14; the book says "repeat the derivation") | C06 | ★★ | 12 | D05–D09 moves with a salt equation, linear EOS (N53), uniqueness of a linear BVP (primer) | Ra is redefined with dT̄/dz (negative when heated from below) so (11.45) carries −Ra; Rs uses κ_s, Rs′ uses κ; identical equations and conditions ⇒ T̂ ∝ ŝ; σ = 0 is valid for fingers only (the diffusive onset is oscillatory) | notebook · `salt_fingers` |
| D14 | linearised Taylor–Couette equations (11.50) and the narrow-gap system (11.51) from (11.47)–(11.49) (a-D16; the book defers to Chandrasekhar) | C07 | ★★★ | 15 | cylindrical operators (R15, P186), linearisation, normal modes in z, operator elimination of û_z and p̂ (P177), narrow-gap approximation (primer), sympy (P40) | linearising −ũ_φ²/R gives −2U_φu_φ/R (factor 2); (dU_φ/dR + U_φ/R)u_R is 2A, the base vorticity; continuity needs 1/R (slip #12); in the narrow gap U_φ/R → Ω(x) is linear in x, giving (1 + αx); û_φ is rescaled so Ta appears on one side only | notebook · `taylor_couette_onset` |
| D15 | the Taylor number (11.52), its narrow-gap inner-only form 2(Ω₁R₁d/ν)²(d/R₁), and μ → 1 reducing (11.51) to (11.39) with Ta for Ra (a-D17) | C07 | ★ | 6 | narrow-gap approximation (primer), algebra | R₂² − R₁² ≈ 2R₁d; α = 0 gives exactly the Bénard pair with the same rigid conditions, hence 1708; μ = 1 itself is solid-body rotation (Ta = 0) — it is the shape of the equations that reduces, not the flow | notebook · `taylor_couette_onset` |
| D16 | stratified perturbation equations (11.55)–(11.57) from the inviscid Boussinesq set (a-D18; slip #4 corrected) | C08 | ★ | 7 | linearisation, material derivative (ch03), N² (R18), stream function (R19) | (u·∇)(Ue_x) = wU′e_x; the w-equation has ∂p/∂z, not ∂p/∂x (slip #4); signs with w = −∂ψ/∂x | notebook |
| D17 | Taylor–Goldstein equation (11.61) from (11.58)–(11.60) (a-D19) | C08 | ★★ | 9 | normal modes (C01), differentiation, complex algebra (P153) | the U′ψ̂′ terms cancel after differentiating (11.58); (ik)² = −k²; dividing by U − c assumes c is not a value of U (critical layer) | notebook · `richardson_shear_instability` |
| D18 | Miles–Howard: (11.63) → (11.64) → (11.65) → Ri > ¼ everywhere ⇒ stable (11.66)–(11.67) (a-D21) | C09 | ★★★ | 14 | product and chain rule with fractional powers (P38, P49, P155), integration by parts (P218a), real/imaginary parts (primer), contradiction (primer), sympy (P40) | "after some rearrangement" is written out (the ½U″ and −¼U′²/(U − c) terms); φ(0) = φ(d) = 0 follows from (11.62); the imaginary part of (U − c) is −c_i; the branch of (U − c)^{1/2} does not matter when c_i ≠ 0; the extraction drops a minus the page has | notebook · `richardson_shear_instability` |
| D19 | Howard's semicircle (11.68)–(11.72) and the growth bound kc_i ≤ (k/2)(U_max − U_min) (a-D22; slip #11 corrected) | C10 | ★★★ | 14 | product rule (P38), integration by parts (P218a), real/imaginary parts (primer), completing the square into a circle (primer) | multiply the transformed equation by 1/(U − c) first to reach the divergence form; the k² term carries (U − c)², not (U − c) (slip #11); the weight Q ≥ 0 must sit inside the inequality from the start; the hidden assumption N² ≥ 0; use (11.70) ⇒ ∫UQ = c_r∫Q twice | notebook · `inviscid_shear_criteria` |
| D20 | perturbation equations (11.73)–(11.77) of a parallel viscous flow with 3-D normal modes (a-D23) | C11 | ★ | 7 | dimensionless NS (R21), linearisation, normal modes (C01) | the parallel base flow satisfies 0 = −∂P/∂x + U″/Re exactly (Poiseuille, Couette) but only approximately for boundary layers; ∇² → d²/dy² − (k² + m²) | notebook |
| D21 | Orr–Sommerfeld equation (11.79) from (11.77) with m = ŵ = 0, û = φ′, v̂ = −ikφ (a-D25; the book writes "this effort yields") | C11 | ★★ | 10 | differentiation, complex algebra (P153), operator elimination of p̂ (P177), sympy (P40) | differentiate the x-equation in y, multiply the y-equation by ik, subtract; U″φ appears from (U φ′)′; the viscous side is (φ⁗ − 2k²φ″ + k⁴φ)/(ikRe) — the sign and the i are the classic slip (planted v̂ = +ikφ fails) | notebook · `orr_sommerfeld_neutral_curve` |
| D22 | Rayleigh's inflection-point theorem (11.83)–(11.84): divide by U − c, multiply by φ*, integrate by parts, take the imaginary part (a-D27) | C12 | ★★ | 9 | integration by parts with φ = 0 at the walls (P218a), real/imaginary parts (primer), sign-definite integrals (primer) | Im{U″/(U − c)} = c_iU″/∣U − c∣² needs the conjugate written out; the theorem needs c_i ≠ 0; a sign-definite integrand forces U″ to change sign inside the open interval | notebook · `inviscid_shear_criteria` |
| D23 | disturbance kinetic-energy equation (11.88) and its 2-D form from (11.96) (a-D30; the book leaves it to Ex. 11.13) | C14 | ★★ | 12 | index notation and the divergence theorem (ch02), product rule (P38), periodic averaging (primer) | subtract the basic-state equation first; u_iU_j∂_ju_i = ∂_j(½u_i²U_j) needs ∂_jU_j = 0; νu_i∇²u_i = ∂_j(νu_i∂_ju_i) − ν(∂_ju_i)²; divergence terms vanish on walls and cancel between periodic ends; in 2-D only −∫uvU′ survives | notebook · `orr_sommerfeld_neutral_curve` |
| D24 | the Lorenz system (11.91) from the truncation (11.90): r = Ra/Ra_c(k), b = 4π²/(π² + k²), time scaled by (π² + k²)κ/d² (a-D31; the book writes "Lorenz finally obtained") | C15 | ★★★ | 15 | vorticity–stream function form (ch06), Galerkin truncation (primer), orthogonality of sines (P151), Ra_c(k) = (π² + k²)³/k² from C05, sympy (P40) | the Jacobian term ∂(ψ, T′)/∂(x, z) produces the sin 2πz mode (Z) — every other harmonic is dropped; the sign convention u = −∂ψ/∂z (slip #9); X, Y, Z are rescaled so the constants vanish; b = 8/3 exactly at k² = π²/2 | notebook · `lorenz_attractor` |
| D25 | Lorenz fixed points, the pitchfork at r = 1 and the Hopf point r_H = Pr(Pr + b + 3)/(Pr − b − 1) ≈ 24.74 (a-D32; the book says only "if r is large") | C15 | ★★ | 11 | Jacobian matrix of an ODE system (primer), characteristic polynomial, purely imaginary roots of a cubic (primer), P214 | the origin's Jacobian splits into a 2 × 2 block and −b; the convection states' cubic λ³ + (Pr + b + 1)λ² + b(r + Pr)λ + 2bPr(r − 1) = 0; the Hopf condition a₂a₁ = a₀ needs Pr > b + 1; it is subcritical (no small stable cycle) | notebook · `lorenz_attractor` |

## 4c. Derivations demoted to statements (first column is the A parent in bold, e.g. **C20** — never a bare ID or a D id: the parser reads a bare first-cell ID as an item and blanks its tier)
| A parent | Analysis §2b item | Result stated (Eq.) | Stated in (B item) | Why not written out |
|---|---|---|---|---|
| **C02** | a-D4 mixing energy | E_final = (2/3)E_initial with ∫U dz unchanged (p. 482) | N22 | the book writes both energies; one definite integral shown in the code cell (`quad` cross-check) |
| **C03** | a-D6 Ra from scaling | w ~ κ/d, buoyancy/viscous ~ Ra | N31 | an order-of-magnitude argument the book states; numbers for water and air instead |
| **C07** | a-D15 ring interchange | ΔE = (Γ₂² − Γ₁²)(1/r₁² − 1/r₂²)/(8π²) ⇒ dΓ²/dr < 0 unstable | N59, R14 | the book writes the energies; the factorisation is one line in the B text and checked by `ring_interchange_energy` vs E_f − E_i |
| **C08** | a-D20 conjugate pairs | c eigenvalue ⇒ c* eigenvalue for (11.18), (11.61), (11.81) | N17, N73, N92 | one move (conjugate an equation with real coefficients, P81); shown by the spectra |
| **C11** | a-D24 Squire's theorem | (11.78): k̄ = √(k² + m²), k̄Re̅ = kRe, Re̅ ≤ Re, k̄c_i ≥ kc_i | N88, N89 | a B item; the transformation is stated and checked numerically (3-D solver = 2-D solver at Re̅) |
| **C12** | a-D26 Rayleigh equation | (11.81) as Re → ∞ in (11.79) | N91 | one move (drop the 1/Re side); its singular nature (order 4 → 2) is said in words |
| **C12** | a-D28 Fjørtoft | (11.85) + (c_r − U_I) × (11.84) ⇒ (11.86) | N93 | the book writes the moves; a B item stated in three lines (why any constant may multiply (11.84)) |
| **C12** | a-D29 cat's-eye expansion | ψ̂ ≈ ½(y − y_c)²U′(y_c) + Aφ(y_c)cos kx; eye width 2√(Aφ_c/U′_c) | N97 | one Taylor move (U(y_c) = c kills the linear term); checked by `cats_eye_streamfunction(exact=True)` |
| **C12** | a-D36 piecewise-linear shear layer (★★★) | c₀² = ((U₁ − U₃)/(2kh))²{(kh − 1)² − e^{−2kh}} (Ex. 11.11) | N126 | a long exercise used only as a V1 check of `rayleigh_eigs`; the jump-condition idea is named |
| **C12** | a-D37 parity of symmetric profiles | (11.95): sinuous/varicose modes, ±c for antisymmetric U | N127 | exercise hint; shown by `rayleigh_eigs(parity=...)` |
| **C15** | a-D34 logistic map | x* = 1 − 1/A stable for 1 < A < 3; period 2 at A = 3; A₂ = 1 + √6 | N117, N129 | exercise material; ∣f′(x*)∣ = ∣2 − A∣ stated with the iterated-map primer and the cobweb |

## 5. Interactive explainers (5–10 + backup)
Nine explainers; every A item is attached to at least one. Every window shows its equations with numbers, has an Explain tab
("Explanation & interpretation", numbered sections with the reader's numbers and a "Reading the current setting" paragraph), a synced
Code tab, a Derivation tab for its D rows, a 4–8-step walkthrough and ≥ 3 check questions. **Colour code across the chapter:** growing /
unstable = rose, decaying / stable = teal, neutral / marginal curve = purple (accent), base state and ghosts = muted, buoyancy or
density = blue, shear or kinetic energy = orange, salt = amber. Reference explainers (`interactive-viz` §5): FDV =
`forced_damped_vibrations.html`, AFE = `angular_frequency_explorer_1.html`, APS = `amplitude_phase_second_order_II_3.html`. **JS physics
without linear algebra** for E1, E2, E4, E5 (approximate formula + table), E6 (Ri and verdict live; growth from a table), E7 (criteria live;
eigenvalues from a table), E9; E3 computes the rigid determinant live with complex cos/cosh; E8 reads cached Orr–Sommerfeld tables. Every
table is ours, labelled, parity-checked against the fluidpy function at ≥ 2 points (`selftest` rows).

### E1 · normal_mode_growth
- A: C01, C05 (also shows N04, N05, N06, R09 (11.19), N121 Rayleigh–Taylor, R10 vortex sheet, N40 exchange of stabilities) · **Confusion removed:** "what does it mean for a flow to be 'unstable', and why does the answer depend on the wavelength?" — a normal mode e^{ikx + σt} grows if σ_r > 0 and travels or oscillates through σ_i; a flow is stable only if σ_r ≤ 0 for every k, and the critical point is where the σ_r(k) curve first touches zero.
- **Why interactive:** dragging k along the σ_r(k) curve makes the picture grow, decay or travel at once; raising the control parameter (density contrast, ΔU, Ra) lifts the whole curve through zero and the reader sees *which* wavelength goes first — a static σ(k) plot hides the link between one point on the curve and the motion it means.
- **Stage:** (1) the disturbance itself: an interface η = ζ₀e^{σ_r t}cos(kx + σ_i t) or a Bénard roll (ψ, T′) whose amplitude follows e^{σt}; (2) σ_r(k) (bold) and σ_i(k) (dashed) with the current k dot, zero line, unstable band shaded rose; (3) (hidePortrait) the complex σ-plane with the roots moving as k changes and the imaginary axis as the stability border.
- **Controls:** system (modes: interface / Rayleigh–Taylor, vortex sheet with gravity, Bénard free–free) · k · control parameter (ρ₁/ρ₂, ΔU or Ra, log) · Pr (Bénard) · transport t.
- **Equations:** (11.1) $u=\hat u(z)e^{ikx+imy+\sigma t}$, $\sigma=-i\lvert\mathbf K\rvert c$; (11.19) $c=\pm[(\rho_2-\rho_1)g/((\rho_2+\rho_1)k)]^{1/2}$; (11.20); (11.44) $\mathrm{Ra}=(\pi^2+K^2)^3/K^2$ with the growth quadratic of D12 $(\sigma+a^2)(\sigma/\Pr+a^2)a^2=\mathrm{Ra}K^2$.
- **Mirrors** `ST.sigma_from_c`, `ST.stability_class`, `ST.marginal_type`, `ch11.kh_phase_speed`, `ch11.benard_free_free_sigma` (§8).
- derivations: D01, D12 · depth features: explain (σ at your k with every term substituted; the band of growing k; the fastest k and its e-folding time; the regime reading), code, **linked views**, **transport**, **modes**, **presets** (stable interface, upside-down Rayleigh–Taylor, vortex sheet, Bénard Ra = 657.5 marginal, Ra = 2000), **status** ("unstable: σ_r = 0.42 s⁻¹ at λ = 3.1 cm; band 2–9 cm" / "marginal: σ_r = 0 at K = π/√2"), **inspector** (click the σ curve: the arithmetic of σ at that k) · follows: AFE (modes + linked views on one clock) + FDV (Explain tab) · **aha:** instability is a property of the whole curve σ_r(k): as the control parameter rises the curve touches zero first at one wavenumber, and that wavelength is the pattern you see at onset.

### E2 · kelvin_helmholtz_boundary
- A: C02 (also shows R09, R10, N16, N17, N18, N19, N120, N121, R11 named) · **Confusion removed:** "why does wind make waves on water only above a threshold, while a shear between equal densities is always unstable?" — gravity stabilises long waves and surface tension short ones, shear destabilises in proportion to kΔU²; the boundary is where they balance and its lowest point is the minimum wind shear (≈ 6.7 m/s for air over water at λ ≈ 1.7 cm).
- **Why interactive:** two sliders (ΔU, density ratio) and two toggles (surface tension, finite depth) move both the stability-boundary curve ΔU_min(k) and the current point; on the c-plane the two real wave speeds slide together, collide on the real axis and split into a growing/decaying pair exactly when the point crosses the boundary — the collision is the instability, and only motion shows it.
- **Stage:** (1) two streams with the interface wave growing (animated), arrows for U₁, U₂, frame toggle lab / moving with (U₁ + U₂)/2 (Fig. 11.3); (2) the (k, ΔU) plane with the boundary curve and the unstable region shaded, or kc_i(k) with k_c marked; (3) (hidePortrait) the complex c-plane with c₊ and c₋.
- **Controls:** ΔU · ρ₁/ρ₂ (chips: air/water, thermocline 0.998, equal) · surface tension on/off · lower depth h (optional) · wavelength · transport.
- **Equations:** (11.18), $g(\rho_2^2-\rho_1^2)<k\rho_1\rho_2(U_2-U_1)^2$, (11.19), (11.20) $c=\tfrac12(U_1+U_2)\pm\tfrac i2(U_2-U_1)$, the Ex. 11.1 form with $\sigma_s k$ and $\coth kh$, $\Delta U_{\min}^2=2\sqrt{g\Delta\rho\,\sigma_s}(\rho_1+\rho_2)/(\rho_1\rho_2)$.
- **Mirrors** `ch11.kh_phase_speed`, `ch11.kh_critical_k`, `ch11.kh_min_shear` (§8), `ch11.vortex_sheet_c`.
- derivations: D02, D03, D04 · depth features: explain (the discriminant with your numbers term by term — gravity, surface tension, shear; c₊, c₋; k_c; growth rate and e-folding time; the reading), code, **linked views**, **transport**, **presets** (air over water 5 m/s and 8 m/s, equal densities, ocean thermocline, Rayleigh–Taylor upside down), **status** ("stable: gravity + surface tension beat shear by 30 %" / "unstable for 0.9 < λ < 3.4 cm"), **terms** (bars: gravity term, surface-tension term, shear term in the discriminant), **inspector** (click a (k, ΔU) point: the discriminant arithmetic) · follows: FDV (Explain with numbers) + APS (a complex-plane view like its Nyquist window) · **aha:** the discriminant decides everything — when shear's kinetic term beats the restoring terms the two wave speeds merge and become a complex pair, and one of them grows.

### E3 · benard_neutral_curve
- A: C04, C03, C05 (also shows N25 Ra with the Γ sign, N41–N46, N44 even/odd, N50 rigid–free, N51 planforms, N123 odd mode) · **Confusion removed:** "why does convection start near Ra ≈ 1708 with cells about as wide as the layer is deep, and why does the number change with the boundaries?" — every K has its own marginal Ra(K): narrow cells pay viscosity and diffusion, wide cells couple buoyancy weakly; onset is the bottom of that valley, and rigid walls cost more than free ones.
- **Why interactive:** the reader drags a point in the (K, Ra) plane: inside the curve the rolls grow, outside they decay, and the third view shows the determinant's imaginary part crossing zero as Ra passes the curve; switching boundaries moves the whole valley (657.5 → 1100.65 → 1707.76). A static Fig. 11.10 shows the curve but not that each point is a whole eigenproblem.
- **Stage:** (1) the layer with roll streamlines over the temperature anomaly (even mode one row, odd mode two), amplitude following e^{σt}; (2) neutral curves Ra(K) for the three boundary pairs, minima marked, current point, unstable region shaded; (3) (hidePortrait) Im det(Ra) at the current K with its root (rigid), or the free–free bars (π² + K²)³ vs K².
- **Controls:** boundaries (chips: rigid–rigid, free–free, rigid–free) · K · Ra (log) · parity (even/odd) · transport.
- **Equations:** (11.21), (11.36)–(11.37), (11.39)–(11.42), the 3 × 3 determinant, (11.44), dRa/dK² = 0, Ra_c = 27π⁴/4.
- **Mirrors** `ch11.benard_determinant`, `ch11.benard_marginal_Ra_det`, `ch11.benard_free_free_Ra`, `ch11.benard_free_free_critical`, `ch11.benard_critical` (rigid–free and odd from a table).
- derivations: D09, D10 (★★★, all 14 steps), D11 · depth features: explain (Ra from your ΔT, d and fluid with both Γ conventions; q₀, q, q* at your K; det value; marginal Ra(K); distance from onset; cell width 2π/K; the reading), code, **linked views**, **modes** (boundaries), **presets** (rigid onset 1707.76 at K = 3.117, free onset 657.5 at π/√2, Ra = 5000 band, odd mode), **status** ("Ra = 2500 > Ra(K) = 1784: rolls grow"), **inspector** (click the Ra(K) curve: the determinant arithmetic), **transport** · follows: APS (numbered derivation with live numbers) + AFE (linked views, presets) · **aha:** onset is the bottom of a valley — the first wavenumber whose marginal Ra the heating reaches.

### E4 · salt_fingers
- A: C06 (also shows N53–N57, N54 the diffusive regime, N55 Ra sign flip) · **Confusion removed:** "how can a column that is lighter on top (statically stable) still overturn?" — heat diffuses about a hundred times faster than salt: a sinking warm salty finger loses its heat to the surroundings but keeps its salt, so it stays heavy and keeps sinking; (11.46) says when this wins.
- **Why interactive:** dragging (dT/dz, dS/dz) across a regime map crosses the density-stable line and the (11.46) line separately, so the reader sees a region that is stable by density yet finger-unstable; the parcel animation shows its temperature relaxing while its salt stays; the κ_s/κ slider opens and closes the finger wedge.
- **Stage:** (1) a column with a displaced parcel (inside/outside T and S colour bars) moving under its own buoyancy; (2) the regime map (α dT/dz vs β dS/dz, or R_ρ) with the static-stability line, the finger line (11.46) and the diffusive region; (3) (hidePortrait) the three roots σ of the free–free cubic in the complex plane.
- **Controls:** dT/dz · dS/dz (or density ratio R_ρ) · κ_s/κ · Pr · K (optional).
- **Equations:** EOS $\tilde\rho=\rho_0[1-\alpha(\tilde T-T_0)+\beta(\tilde s-s_0)]$, (11.45), $\mathrm{Rs}-\mathrm{Ra}=\tfrac{27}{4}\pi^4$, (11.46).
- **Mirrors** `ch11.salt_finger_unstable`, `ch11.double_diffusive_margin`, `ch11.double_diffusive_sigma`, `ch11.salt_finger_regime` (§8).
- derivations: D13 · depth features: explain (Ra and Rs with your gradients, the margin Rs − Ra − 657.5, R_ρ, the density gradient, the fastest σ; the reading), code, **linked views**, **presets** (subtropical thermocline: warm salty over cold fresh — fingers; Arctic: cold fresh over warm salty — diffusive; κ_s = κ: no fingers; statically unstable), **status** ("statically stable (R_ρ = 1.8) but finger-unstable: Rs − Ra = 2.4 × 10⁴ > 657.5"), **transport**, **inspector** (click the map: the margin arithmetic) · follows: FDV (regime status + interpretation) + `overfitting_curves` (minimal map with a verdict line) · **aha:** static stability is not enough when the two ingredients of density diffuse at different speeds.

### E5 · taylor_couette_onset
- A: C07 (also shows R14 Rayleigh's criterion, N59 ring interchange, R16 circular Couette, N64–N68, N124 Galerkin route) · **Confusion removed:** "why does spinning the inner cylinder make stacked doughnut vortices while spinning the outer one does not, and why does 1708 come back?" — swapping two rings releases energy when Γ² falls outward (Rayleigh); the centrifugal force plays gravity, and viscosity sets a threshold Ta_c that tends to the Bénard value for co-rotation.
- **Why interactive:** dragging (Ω₂/Ω₁, Ta) moves the point across the Rayleigh line and the viscous boundary separately; the profile view shows where dΓ²/dr < 0 and the ring-pair inspector gives ΔE for any two radii; Taylor vortices appear in the gap only past the viscous curve.
- **Stage:** (1) meridional cross-section of the gap with Taylor vortices (narrow-gap eigenfunction from a table) and a top view of the rotating cylinders; (2) the (μ = Ω₂/Ω₁, Ta) plane with the Rayleigh line, (11.54) curve, exact narrow-gap dots and the current point; (3) (hidePortrait) U_φ(R) and Γ²(R) with the ring-interchange ΔE bars.
- **Controls:** Ω₁ · Ω₂ (or μ) · gap ratio d/R₁ · ν (chips: water, glycerol) · ring pair r₁, r₂ (inspector).
- **Equations:** (11.49) $U_\varphi=AR+B/R$, $d\Gamma^2/dr<0$, (11.51), (11.52), (11.54) $\mathrm{Ta}_{cr}=1708/(\tfrac12(1+\Omega_2/\Omega_1))$.
- **Mirrors** `ch11.taylor_number`, `ch11.taylor_critical_approx`, `ch11.taylor_critical` (table), `ch11.rayleigh_circulation_criterion`, `ch11.ring_interchange_energy`, `core.laminar.circular_couette`.
- derivations: D14 (★★★, all 15 steps), D15 · depth features: explain (A, B and Γ² at your speeds, Ta, Ta_c from (11.54) and from the exact table with the error, ΔE for the inspected rings, the reading), code, **linked views**, **presets** (inner only, co-rotation μ = 0.5, counter-rotation μ = −0.5, on the Rayleigh line), **status** ("Rayleigh-unstable but Ta = 1200 < Ta_c = 3416: viscosity holds it"), **inspector** (click two radii: ΔE arithmetic), **transport** (vortex animation) · follows: AFE (system animation + response curve + presets) · **aha:** centrifugal instability is buoyancy in disguise — Γ² plays −ρ, the Rayleigh line is neutral stratification, and viscosity adds the 1708-type threshold.

### E6 · richardson_shear_instability
- A: C08, C09 (also shows R18 N², R20 Ri, N73 conjugate pairs, N78 the (k, J) growth map, N20 billow clouds) · **Confusion removed:** "why does Ri = ¼ decide whether a stratified shear layer breaks into billows — and why is it not a sharp switch?" — Ri = N²/U′² compares the work against buoyancy with the shear's energy; if Ri > ¼ everywhere no mode can grow (a guarantee), but Ri < ¼ somewhere only allows growth (necessary, not sufficient).
- **Why interactive:** the reader slides the stratification J (= Ri at the centre) and watches the Ri(z) profile dip below the ¼ line at the same moment the unstable tongue of the (k, J) map is reached and the billow animation starts growing; moving the density layer's thickness shows Ri_min and the growth tongue changing together.
- **Stage:** (1) a shear layer with density layers deforming into KH billows (linear mode with amplitude e^{kc_i t}); (2) profiles U(z), N²(z) and Ri(z) with the ¼ line and Ri_min marked; (3) (hidePortrait) growth map kc_i(k, J) for U = tanh z, N² = J sech²z (cached table) with the current point and J = ¼.
- **Controls:** J · k · thickness ratio of density to shear layer (chips) · transport.
- **Equations:** (11.56) $N^2=-\frac{g}{\rho_0}\frac{d\bar\rho}{dz}$, (11.61), (11.64), (11.65), (11.66) $\mathrm{Ri}=N^2/U'^2$, (11.67) $\mathrm{Ri}>\tfrac14$.
- **Mirrors** `ch11.gradient_richardson`, `ch11.miles_howard_stable`, `ST.taylor_goldstein_eigs` (table parity), `ch11.tg_growth_map` (§8).
- derivations: D17, D18 (★★★, all 14 steps) · depth features: explain (N², U′ and Ri at the centre with your J; Ri_min and where; the growth rate from the map; the e-folding time for a thermocline and an atmosphere case; the reading), code, **linked views**, **presets** (J = 0 pure KH, J = 0.1, J = 0.24, J = 0.26, thermocline, billow-cloud layer), **status** ("Ri_min = 0.18 < ¼: instability allowed; kc_i = 0.06" / "Ri > ¼ everywhere: stable by Miles–Howard"), **transport**, **inspector** (click Ri(z): its arithmetic) · follows: FDV (live Explain) + AFE (linked views on one clock) · **aha:** ¼ is a guarantee of stability, not a trigger of instability.

### E7 · inviscid_shear_criteria
- A: C12, C10 (also shows N91 Rayleigh's equation, N93 Fjørtoft, N94 six profiles, N95 sin y, N96 critical layer, N97 cat's eye, N82 c_r in range, N126 piecewise layer, N127 parity) · **Confusion removed:** "which velocity profiles can be unstable without viscosity, and where can the speed of an unstable wave lie?" — a growing mode needs an inflection point (Rayleigh) where the vorticity peaks (Fjørtoft), and every unstable c lies inside Howard's semicircle on [U_min, U_max]; neutral modes carry a critical layer with cat's-eye streamlines.
- **Why interactive:** switching profiles moves the inflection marker, flips the verdict badges and moves the eigenvalues (from a table) on the c-plane, always inside the semicircle; clicking a neutral mode opens the cat's eye in the frame moving with c — the reader checks each theorem against an example instead of reading it.
- **Stage:** (1) U(y) with U″ and vorticity, the inflection point and U_I, the Fjørtoft integrand (U − U_I)U″ shaded by sign; (2) complex c-plane with the semicircle, the velocity range and the profile's eigenvalues; (3) (hidePortrait) cat's-eye streamlines around the critical layer.
- **Controls:** profile (chips: Blasius-like, Poiseuille, tanh shear layer, Bickley jet, wall-vorticity-max, sin y) · k · b (sin y half-width) · cat's-eye amplitude A.
- **Equations:** (11.81), (11.83)–(11.84), (11.85)–(11.86), (11.71), (11.72), (11.87).
- **Mirrors** `ST.inflection_points`, `ch11.rayleigh_criterion`, `ch11.fjortoft_criterion`, `ST.howard_semicircle`, `ST.in_howard_semicircle`, `ch11.cats_eye_streamfunction`, `ST.rayleigh_eigs` (table).
- derivations: D19 (★★★, all 14 steps), D22 · depth features: explain (U″ sign change, U_I, the Fjørtoft sign, the semicircle centre and radius with your profile, the leading eigenvalue and its distance from the circle, the eye width; the reading), code, **linked views**, **presets** (tanh: both pass, Poiseuille: no inflection, wall-vorticity-max: Rayleigh yes / Fjørtoft no, Bickley sinuous, sin y with 2b < π), **status** ("inflection at y = 0, vorticity maximum: instability possible; c = 0 + 0.43i inside the semicircle"), **inspector** (click U(y): U″ and (U − U_I)U″ arithmetic) · follows: `fid_formula_lab` (verdict badges + term colours) + APS (numbered live explanation) · **aha:** inviscid instability is written in the profile's curvature, and the unstable wave speeds cannot leave a half-disc set by the slowest and fastest fluid.

### E8 · orr_sommerfeld_neutral_curve
- A: C13, C11, C14 (also shows N88 Squire, N99 tanh, R23 Bickley, N100 Couette, N102 Table 11.1, N104 Reynolds stress, N105 Tollmien profile, N107 frequency form, N108 Blasius, R24 Falkner–Skan) · **Confusion removed:** "how can viscosity, which only damps, make plane Poiseuille flow unstable?" — the viscous wall layer shifts the phase of v against u so that the Reynolds stress −⟨uv⟩ draws energy from the mean shear; inside the thumb-shaped neutral curve, above Re_c = 5772, production beats dissipation.
- **Why interactive:** clicking a point (Re, k) on the cached neutral-curve map shows that TS mode travelling over the base profile, its production/dissipation bars with dE/dt = 2kc_iE, and the u–v phase; crossing the curve flips the bars' balance. Switching flows (Poiseuille, Blasius, tanh, Bickley, Couette) recomputes Table 11.1 row by row; the oblique-angle slider applies Squire's map live.
- **Stage:** (1) neutral curve(s) in (Re, k) on log Re with the unstable region shaded, Re_c marked and the current point; (2) the mode: perturbation streamfunction over the base profile, travelling at c_r; (3) (hidePortrait) energy bars (production, dissipation, dE/dt) + Reynolds-stress profile.
- **Controls:** flow (chips) · Re (log) · k · oblique angle (Squire) · transport.
- **Equations:** (11.77), (11.78) $\bar k\,\overline{\mathrm{Re}}=k\,\mathrm{Re}$, (11.79), (11.80), (11.88) and $\frac{d}{dt}\int\tfrac12(u^2+v^2)dV=-\int uv\,U'dV-\Lambda$, Table 11.1.
- **Mirrors** `ST.orr_sommerfeld_eigs` (table + one small-N live parity row), `ch11.poiseuille_neutral_curve` (table), `ST.disturbance_energy_budget` (table), `ST.squire_transform` (live).
- derivations: D21, D23 · depth features: explain (c at your (Re, k) from the table; growth 2kc_i; production, dissipation and their difference; Re̅ for the oblique wave; Table 11.1 row highlighted; the reading), code, **linked views**, **transport**, **terms** (production, dissipation, dE/dt), **inspector** (click the map: the mode's numbers), **presets** (Poiseuille Re_c, Re = 10⁴ k = 1 (Orszag), Blasius Re_δ* = 519, Couette stable, tanh Re_c = 0), **status** ("inside the neutral curve: production 1.4 × dissipation — TS wave grows") · follows: FDV (live Explain) + `fid_formula_lab` (term bars that add to a total) · **aha:** viscosity creates the phase lag that lets the wave tap the mean shear — damping and driving come from the same term.

### E9 · lorenz_attractor
- A: C15 (also shows N112–N117, N119; 3-D) · **Confusion removed:** "how can three deterministic equations be unpredictable?" — past r_H ≈ 24.74 every fixed point is unstable, so the state wanders aperiodically between two convection lobes, and nearby starts separate like e^{λt} until they are unrelated.
- **Why interactive:** the 3-D orbit, the X(t) traces and the log-separation panel run on one clock from two starts 10⁻⁸ apart; the r slider moves the fixed points and changes the verdict (conduction → steady convection → chaos) while the reader watches; a static attractor picture hides both the divergence and the r-dependence.
- **Stage:** (1) 3-D attractor (orbit, two trajectories, fixed points C± and the origin); (2) X(t) for both trajectories and log₁₀∣δ∣ vs t with a slope line; (3) (hidePortrait) the convection roll (ψ, T′) from (11.90) driven by X, Y, Z.
- **Controls:** r · Pr · b (or the roll wavenumber k) · δ₀ (log) · transport.
- **Equations:** (11.90), (11.91), fixed points $\bar X=\bar Y=\pm\sqrt{b(r-1)},\ \bar Z=r-1$, $r_H=\Pr(\Pr+b+3)/(\Pr-b-1)$, $\nabla\cdot\dot{\mathbf s}=-(\Pr+1+b)$.
- **Mirrors** `ch11.lorenz_rhs` (RK4 parity over a short time), `ch11.lorenz_fixed_points`, `ch11.lorenz_hopf_r`, `ch11.lorenz_b`.
- derivations: D24 (★★★, all 15 steps), D25 · depth features: explain (r, b and the fixed points with your settings; Jacobian eigenvalues at C±; r_H; volume contraction rate; the measured divergence slope and the predictability time; the reading), code, **3-D**, **linked views**, **transport**, **presets** (r = 0.5 conduction, r = 10 steady convection, r = 20 transient chaos, r = 28 Lorenz), **status** ("r = 28 > r_H = 24.74: both convection states unstable — chaos"), end-of-run card (time until the two runs differ by 50 %), **inspector** (click the orbit: the state and the Jacobian there) · follows: `ddpm3_unet_3d` (orbitable 3-D scene) + AFE (linked views, end-of-run summary) · **aha:** deterministic is not predictable — the rule is exact, but the error doubles every ≈ 0.8 time units.

### B1 · period_doubling_route (backup)
- A: C15 (also shows N117 Feigenbaum, N129 logistic map, N114 bifurcations) · **Confusion removed:** "what does a 'route to chaos' look like, and why is 4.669 universal?" — each period doubling happens when the fixed point of the twice-iterated map loses stability (∣f′∣ = 1); the gaps between doublings shrink by the same ratio in any one-hump map.
- Stage: cobweb diagram of x_{n+1} = Ax_n(1 − x_n) · bifurcation tree with A_n markers · the δ table from superstable points. Controls: A, x₀, number of iterates, map (logistic / sine). Mirrors `ch11.logistic_map`, `ch11.bifurcation_diagram`, `ch11.period_doubling_points`, `ch11.feigenbaum_estimate`. Derivations none · depth: linked views, presets (A = 2.8, 3.2, 3.5, 3.57, 3.83), status, inspector. Built only if an explainer above fails review.

## 6. Python animations and interactive figures (A ID → what, why, player/figure kind)
Animations (5, `animate` + `show_animation`; ≤ 90 frames, FAST → 40; dpi 80):
| # | A | What moves | Why | Player |
|---|---|---|---|---|
| A1 | C01 (N02) | a damped ball in the four potentials of Fig. 11.1 (bowl, cap, plane, dimple), small and large kicks | stability as a question about small vs large disturbances, before any equation | frames (stop at the escape from the dimple) |
| A2 | C02 (R10, R11) | a KH interface growing at e^{kc_i t} (linear) beside ch05's `sheet_rollup` of the same sheet (nonlinear) | where linear theory is right (early) and where it stops (roll-up) | video |
| A3 | C04, C05 | Bénard rolls (streamlines over T′) for one Ra above and one below the neutral curve, amplitudes e^{σt} | onset as a sign change of σ at fixed K | video |
| A4 | C13, C14 | a TS wave in plane Poiseuille flow at Re = 10⁴, k = 1: perturbation ψ travelling at c_r, with the Reynolds-stress profile and the production/dissipation numbers | the phase shift that feeds the wave, on one clock | video |
| A5 | C15 | two Lorenz trajectories 10⁻⁸ apart on the attractor (X–Z projection) with the roll field from (11.90) and log₁₀∣δ∣(t) | sensitive dependence, made visible | video |

Interactive figures (plotly `slider_figure` / `animate_figure`, precomputed; ≤ 40 steps):
| # | A | What the slider controls | Why |
|---|---|---|---|
| F1 | C01, C02 | ΔU: KH growth rate kc_i(k) for air over water with and without surface tension, k_c marked | "which wavelengths grow" (seed 1) on the published page |
| F2 | C04, C05 | Ra: neutral curves Ra(K) (rigid, free, rigid–free) with the unstable band of K at that Ra highlighted | the band opens at the minimum (seed 3) |
| F3 | C06 | κ_s/κ: regime map with the finger line (11.46) and the static-stability line | the finger wedge opening as salt diffuses more slowly |
| F4 | C07 | μ = Ω₂/Ω₁: Ta_c(k) curves exact vs (11.54), with the Rayleigh verdict | the approximation's quality and the co-rotation limit |
| F5 | C08, C09 | J: Ri(z) profile and the TG growth curve kc_i(k) for tanh/sech² (cached) | growth vanishing as Ri_min passes ¼ (seed 4) |
| F6 | C10, C12 | sin y half-width b (or profile family parameter): Rayleigh eigenvalues on the c-plane with the semicircle | eigenvalues stay inside; growth → 0 as 2b → π |
| F7 | C11, C13 | Re: plane Poiseuille OS spectrum at k = 1 (cached), the TS eigenvalue crossing c_i = 0 | the Y-shaped spectrum and the one mode that crosses |
| F8 | C15 | r: Lorenz trajectory (X–Z) for r from 0.5 to 30 (precomputed) | conduction → convection → chaos (seed 5) |
A live `ipywidgets` cell (kernel only) pairs with F2 (free K, Ra, boundaries) and one with F8 (free r, Pr, b, δ₀). Static figures accompany
every A block (our analogues of Figs. 11.1, 11.3, 11.9–11.10, 11.13, 11.17, 11.20–11.24, 11.26, 11.28–11.31, all generated by our code).

## 7. From-scratch moments
Each shows a hand-written version beside the tested function, followed by `assert np.allclose(...)`.
| A | § | Hand-written | Compared with |
|---|---|---|---|
| **C01** | 11.2 | both forms of (11.1) evaluated with complex arithmetic on a grid; σ from c by σ = −i∣K∣c | `ST.normal_mode`, `ST.sigma_from_c` |
| **C02** | 11.3 | c from the quadratic ρ₁(U₁ − c)² + ρ₂(U₂ − c)² = (g/k)(ρ₂ − ρ₁) by `np.roots` | `ch11.kh_phase_speed` |
| **C03** | 11.4 | Trefethen's Chebyshev matrix in ten lines (points, off-diagonal formula, negative-sum diagonal) and d/dx of x³ | `ST.cheb` |
| **C04** | 11.4 | the 3 × 3 determinant with complex cos/cosh and `brentq` on its imaginary part at K = 3.117 | `ch11.benard_marginal_Ra` (Chebyshev) |
| **C05** | 11.4 | (11.44) minimised by `minimize_scalar` vs 27π⁴/4 and π/√2 | `ch11.benard_free_free_critical` |
| **C06** | 11.5 | the criterion (11.46) in one line for a thermocline example | `ch11.salt_finger_unstable` |
| **C07** | 11.6 | Ta from (11.52) and Ta_c from (11.54) by hand; ΔE of one ring pair | `ch11.taylor_number`, `ch11.taylor_critical_approx`, `ch11.ring_interchange_energy` |
| **C09** | 11.7 | Ri(z) for U = tanh z, N² = J sech²z from the formulas, its minimum by `np.argmin` | `ch11.gradient_richardson`, `ch11.miles_howard_stable` |
| **C10** | 11.7 | the semicircle test (c_r − m)² + c_i² ≤ R² for every computed unstable eigenvalue | `ST.in_howard_semicircle` |
| **C11** | 11.8 | the Orr–Sommerfeld matrices for plane Poiseuille assembled from `cheb` (clamped rows), `scipy.linalg.eig`, leading c at Re = 10⁴, k = 1 | `ST.orr_sommerfeld_eigs` and Orszag's c = 0.23752649 + 0.00373967i |
| **C12** | 11.9 | inflection points by the sign change of U″ (`np.sign`, `np.nonzero`) and the Fjørtoft product | `ST.inflection_points`, `ch11.fjortoft_criterion` |
| **C13** | 11.10 | `brentq` on c_i(Re) at k = 1.02056 near Re = 5772 | `ch11.poiseuille_critical` (cached) |
| **C14** | 11.10 | production −∫uvU′ and dissipation with `np.trapezoid` on a fine re-interpolated grid | `ST.disturbance_energy_budget` (Clenshaw–Curtis) |
| **C15** | 11.14 | an RK4 loop for (11.91); agreement with the adaptive solver up to t ≈ 5, then divergence — itself the lesson | `ch11.lorenz_integrate` |
C08 has a sympy check (`stratified_shear_sympy`) instead of a second numeric version; §11.1 has no computable A item; §11.11–11.13 are
covered inside C13/C14.

## 8. Notes for the implementer
Functions and helpers that the A items' figures and explainers need beyond analysis §4 (B/C functions optional):
1. `ch11.benard_free_free_sigma(K, Ra, Pr, n=1)` → the two real roots of (σ + a²)(σ/Pr + a²)a² = RaK², a² = n²π² + K² (D12; E1 Bénard mode, V1 test of `benard_growth_rate`).
2. `ch11.kh_min_shear(rho1, rho2, g=G0, surface_tension=0.0)` → (ΔU_min, k*) with ΔU_min² = 2√(g(ρ₂ − ρ₁)σ_s)(ρ₁ + ρ₂)/(ρ₁ρ₂) at k* = √(g(ρ₂ − ρ₁)/σ_s) (∞ boundary without surface tension: returns k_c(ΔU) instead); `ch11.kh_stability_boundary(k, rho1, rho2, g, surface_tension, h)` → ΔU_min(k) for E2's boundary curve.
3. `ch11.normal_mode_growth(system, k, **params)` → complex σ for system ∈ {"interface", "kh", "benard_free"} — one entry point for E1 and F1.
4. `ch11.salt_finger_regime(dTdz, dSdz, alpha, beta, kappa, kappa_s, nu, d, g)` → {regime: "stable" / "fingers" / "diffusive" / "overturning", R_rho, margin} for E4's status and F3.
5. `ch11.taylor_critical_table(mus)` (exact narrow gap, cached to `reference/ch11/taylor_critical.csv`, ours) and `ch11.taylor_eigenfunction(k, mu, x)` for E5's vortex picture.
6. `ch11.richardson_profiles(kind, J, R)` → (U, U″, N², Ri) callables for the tanh/sech² family; `ch11.tg_growth_map(ks, Js, R)` cached (E6, F5).
7. `ch11.inviscid_profile(name, **p)` → (U, U″) for the six Fig. 11.21 stand-ins, tanh, Bickley, sin y; `ch11.rayleigh_spectrum_table(names, ks)` cached (E7, F6).
8. `ST.os_mode(k, Re, profile, N)` → (c, φ normalised so max∣û∣ = 1) and `ch11.ts_wave_fields(x, y, t, k, c, phi, amp)` → (ψ, u, v) for A4 and E8; `ch11.neutral_curve_tables()` writing Poiseuille, Blasius, tanh, Bickley neutral curves + mode samples + budget terms to `reference/ch11/` (ours, labelled).
9. `ch11.lorenz_r_sweep(r_values, t_end)` (precomputed X–Z samples for F8) and `ch11.lorenz_predictability_time(delta0, threshold)`.
10. `ch11.book_slips()` keyed "S1"…"S12" (the notebook prints them as "slip #1"…"slip #12"), each with printed form, correct form and an evaluator for the planted-variant tests.
11. JS tables for the explainers (ours): rigid–free and odd-mode neutral curves (E3), exact Taylor Ta_c(μ) and one eigenfunction (E5), TG (k, J) map (E6), Rayleigh eigenvalues per profile (E7), OS neutral curves, mode samples and budgets (E8); each with ≤ 4 significant figures and ≥ 2 parity rows against the fluidpy function.

Printed-slip callouts (shown as "book prints X, correct is Y" or ⚠️ Common confusion):
| slip | Where shown | Callout |
|---|---|---|
| slip #1 (p. 489, derivative line of the even solution) | C04 / N46, D10 | third term printed with B; correct C(q*² − K²)² cosh q*z (planted variant moves Ra_c by > 1 %) |
| slip #2 (p. 490) | C05 / N48, D11 | W = A sin nπz vanishes at z = ±½ only for even n; correct sin nπ(z + ½) (n = 1: cos πz) |
| slip #3 (p. 490) | C05 / N49, D11 | dRa/dK² printed with a factor 3 on the second term; the printed form has no root |
| slip #4 ((11.55)) | C08 / N71, D16 | w-equation prints ∂p/∂x; correct ∂p/∂z |
| slip #5 (p. 480) | C02 / R09 | "(7.96)" should be ch07 (7.95), the interface dispersion relation |
| slip #6 (p. 520) | C13 / N105 | Tollmien middle branch 1 − b[1 − η²] is discontinuous; 1 − b(1 − η)² is continuous with matching slope |
| slip #7 ((11.93)) | C07 / N124 | operator printed squared; (11.51) at σ = 0 gives the first power |
| slip #8 (p. 507) | C10 / D19 | the weight Q must sit inside the inequality from the start; "<" should be "≤" |
| slip #9 (stream-function signs) | C08 / R19, C11 / R22, C15 / N115, D24 | not an error: three conventions (§9); each function states its own |
| slip #10 (sign of Γ) | C03 / N25, C06 / N55 | (11.21) Γ = −dT̄/dz, opposite to ch01's Kundu lapse rate; §11.5 flips Ra again (§9) |
| slip #11 (p. 506) | C10 / N80, D19 | −k²(U − c)F should be −k²(U − c)²F (sympy decides) |
| slip #12 ((11.47), (11.50)) | C07 / R15, D14 | continuity needs (1/R)∂(Rũ_R)/∂R |

## 9. Implementation guidance (conventions, runtime budget, caching)
Consistent with analysis §4 (42 rows) and §5: `fluidpy/core/stability.py` (ST: `cheb`, `generalized_eigs`, normal-mode vocabulary,
`orr_sommerfeld_eigs`, `rayleigh_eigs`, `taylor_goldstein_eigs`, `max_growth`, `neutral_curve`, `critical_point`, semicircle,
inflection, Squire, energy budget) and `fluidpy/ch11_instability.py` (chapter physics, sympy engines, wrappers; re-exports ST names
for `nb.core` imports); scripts `scripts/ch11_*.py` (all `--no-show`). Additions from this curation: items 1–11 of §8.

**Conventions to settle (state each in a ⚠️ callout at first use and in every docstring):**
- **Slips** are "slip #k" in the notebook, page and explainers; `S1…S12` only in code and analysis (no collision with recap IDs R01…).
- **The sign of Γ (slip #10; project rule "show both conventions").** (11.21) uses Γ = −dT̄/dz (> 0 when heated from below); ch01's Kundu
  lapse rate is Γ ≡ dT/dz (Γ_a ≈ −9.8 K/km) and meteorology's is Γ_met ≡ −dT/dz. Code never takes a Γ: `rayleigh_number(..., dT, d, ...)`
  with dT = T_bottom − T_top. The notebook shows a two-row table for one layer (water, d = 5 mm, ΔT = 2 K): Kundu ch01 dT/dz = −400 K/m;
  (11.21) Γ = +400 K/m (which happens to equal the meteorological sign); a worked conversion (negate the number) and Ra ≈ 3.7 × 10³ from
  either. §11.5 changes again: there Ra ≡ gαd⁴(dT̄/dz)/(νκ) is negative when heated from below, so (11.45) carries −Ra — code
  `thermal_rayleigh_signed(dTdz)` and the callout in C06. Never pass `core.stratification`'s Γ into ch11 functions.
- **Four non-dimensionalisations** (one table in C03, repeated as a one-line reminder in C07, C11, C15): §11.4 lengths d, time d²/κ, w
  dimensional until W ≡ (Γd²/κ)ŵ; §11.6 gap d, x = (R − R₁)/d ∈ [0, 1], σ by d²/ν (inferred), d/dR in (11.51) means d/dx; §11.8–11.11 L
  and U₀ per flow (half-width and centreline speed for Poiseuille, δ* and U∞ for Blasius, L and U₀ for tanh/sech²), time L/U₀, pressure
  ρU₀² — Table 11.1's Re uses different lengths per row; §11.14 Lorenz time (π² + k²)κ/d² and rescaled X, Y, Z.
- **Two eigenvalue conventions**: σ (e^{σt}) in §11.4, §11.6, §11.14; c (e^{ik(x−ct)}) in §11.3, §11.7–11.11; σ = −ikc. Functions return
  what the section uses and `sigma_from_c` converts.
- **Stream-function signs (slip #9)**: §11.7 u = ∂ψ/∂z, w = −∂ψ/∂x; §11.8 u = ∂ψ/∂y, v = −∂ψ/∂x; §11.14 u = −∂ψ/∂z, w = ∂ψ/∂x.
- **Notebook symbols for the overloaded ones** (add to `knowledge/notation.md`): K Bénard horizontal wavenumber (non-dim), k streamwise
  or axial wavenumber, κ thermal and κ_s salt diffusivity; α thermal expansion, and the Taylor parameter written α_Ω ≡ Ω₂/Ω₁ − 1 with
  μ ≡ Ω₂/Ω₁ (code `mu`); β haline coefficient (code `beta_S`; not ch10's FTCS β); σ growth rate, σ_s surface tension (code
  `surface_tension`); Γ_T for (11.21)'s temperature gradient where it could meet Γ = 2πrU_θ circulation in §11.6; U₁, U₂ the KH streams,
  U_I the speed at the inflection point (the book writes U₁ in Fjørtoft); r Lorenz's Ra/Ra_c vs R radius; R_n the bifurcation values in
  §11.14; φ is the KH potential in §11.3, the transformed TG variable in §11.7 and the OS amplitude in §11.8 (callout at each); Pr in the
  Lorenz system is the "σ" of the wider literature (callout, since σ is a growth rate here).
- **Hidden assumptions to state in the blocks**: Howard needs N² ≥ 0; Miles–Howard divides by U − c (c_i ≠ 0); exchange of stabilities is
  proved for Bénard (Ra > 0) and co-rotating Taylor only; Squire holds for parallel non-rotating, unstratified flows (§11.7 assumes it);
  parallel-flow OS for boundary layers is an approximation (N109); the double-diffusive margin σ = 0 is for fingers only.

**Computational-cost rules (notebook < 5 min on Colab CPU with FAST):**
- **Critical-point searches** (plane Poiseuille ≈ 30 s, Blasius ≈ 30 s, Bickley, tanh neutral curves) run in `scripts/ch11_neutral_curves.py`
  and are cached to `outputs/ch11/` (npz keyed by a parameter hash) and, for the explainers, to `reference/ch11/*.csv` (ours, labelled);
  the notebook loads the cache and falls back to a coarse run (N = 60, 8 Re values) when absent, saying so. One **slow** test recomputes a
  coarse Poiseuille critical point (≤ 1e-3) so a mutant cannot hide behind the cache (ch10 lesson); fast tests use the Orszag point at
  Re = 10⁴, k = 1 (N = 100, ~0.3 s).
- **Eigen-solvers**: Bénard and Taylor N = 40 (FAST 24) ~10–30 ms each; TG/Rayleigh N = 80–120 (FAST 60); OS N = 100 (FAST 60); never
  N > 150 (ill-conditioning). Every spectrum is filtered by N-convergence, never by a tight ∣c∣ cap alone (the TS mode is weak).
- **Lorenz**: `solve_ivp` DOP853 rtol 1e-10, t_end = 50 (FAST 25); Lyapunov estimate optional and cached; assert invariants and
  statistics, never pointwise values beyond t ≈ 20.
- **Animations**: ≤ 90 frames (FAST 40), dpi 80, each < 6 MB; the ch05 roll-up for A2 at N = 200 points (FAST 100).
- **Explainers**: live JS only for closed forms, the 3 × 3 determinant and RK4; everything spectral from tables with parity rows.

**Evidence for A items** follows analysis §6 (minimum: normal modes V1 V7; KH V1 V2 (+ ch07 parity); Bénard growth V1 V3; rigid Ra_c V5 V1
(determinant vs Chebyshev) V3; free–free V1 V2; salt fingers V1 V7; Taylor V1 V3 (μ → 1 = Bénard, Galerkin route); TG V1 V4; Miles–Howard
V2 V1 (random profiles); semicircle V4 V1; OS V5 (Orszag) V3 V4 (energy budget closes); inflection V1 V7; Poiseuille V5 V3; energy budget V4
V2; Lorenz V1 V4 V5), with planted wrong variants for slips #1, #2, #3, #4, #6, #7, #11 and for v̂ = +ikφ, Γ = +dT/dz, Re̅ = k̄Re/k.
