# Chapter 4 — Conservation Laws: curation
(from `analysis/ch04.md` — 185 inventory rows, 124 equation labels (4.1)–(4.119), 60 derivations in §2b (59 + 1 recap
row), 17 candidate A ideas; `book.yaml` ch04, `policy.tier_a_full_treatment: [12, 18]`, `coverage: exhaustive`.
Curated 2026-09-23, concept-curator. Read: `knowledge/CUMULATIVE.md`, `concept_map.md`, `primers.md` (P01–P110),
`notation.md`, `viz_patterns.md`, `knowledge/ch03.md` §8, `analysis/ch03_curation.md` (format).)

Counts: A 15 · B 151 · C 19 · RECAP 12 · SKIP 2 · derivations written out 30 (★ 5 ★★ 23 ★★★ 2) · demoted to statements 27

Tier words (parsed by `tools/nbkit.py`): **CORE 15 · NOTE 156 · RECAP 12 · SKIP 2** = 185 rows. Depth: A 15 (all CORE) · B 151 (142 NOTE + R01, R02, R03, R05, R07, R08, R09, R10, R12) · C 19 (14 NOTE + R04, R06, R11 + S01 + S02). Every inventory row number `[#n]` appears exactly once in §2.

**Decisions that shape this chapter.**
1. **Fifteen A items from the analyst's seventeen.** Two merges keep the A tier a clear minority (15 of 185 rows):
   (a) the rotating-frame acceleration (4.43) and the Coriolis force are B items inside one A block, **C09 = Navier–Stokes
   in a noninertial frame (4.45)** — its derivation chain D14 → D15 ★★★ → D16 is written out in full, and the Coriolis
   projectile (N61, D17) is C09's worked number and animation, so nothing is lost for Ch. 13; (b) the internal-energy
   equation (4.57), the dissipation ε ≥ 0 (4.58) and the entropy production (4.63) form one A block, **C10**, whose
   derivations D19–D23 all stay written out (D20, D21, D23 under rule (b): the book defers them to Exercises
   4.45–4.47). Energy Bernoulli (4.78) is B inside C11 (as the analyst suggested); continuity's corollary ∇·u = 0 (4.10)
   is B inside C02; stress symmetry (4.25) is B inside C07 with its proof D08 kept under rule (b).
2. **Three Bernoulli A items, on purpose.** They answer three different questions a reader will meet again: C05 (4.19)
   *what is constant along a streamline* (pitot, orifice — the engineering workhorse, derived from the stream-tube
   element of Ex. 4.2); C11 (4.71) *when and why Bernoulli holds* (Lamb identity → Bernoulli function → constant on
   streamlines **and** vortex lines, everywhere if irrotational — the gateway to Ch. 5); C12 (4.75) *the unsteady
   potential form* that is Ch. 7's dynamic free-surface condition. One explainer (`which_bernoulli`) ties the four
   forms and their hypotheses together (N102 decision table).
3. **SEEN rows are never A.** The analyst marked 21 rows SEEN. Pure restatements become RECAPs (R01–R12, each reusing
   the earlier chapter's function): material volume (ch03), Gauss (ch02), Cauchy traction (ch02), surface tension (ch01,
   twice), Galilean invariance ⇒ σ(S) (ch03), non-Newtonian (ch01), σ:R = 0 (ch02), Gibbs (ch01), u = ∇φ (ch03),
   hydrostatics (ch01), sphere drag Π groups (ch01). SEEN rows that carry a *new numbered equation* ((4.2), (4.3),
   (4.14), (4.15), (4.26), (4.47)) or a new role (dH/dt = M, Laplace's jump now *derived*, two routes to similarity) are
   B/C NOTEs inside the nearest A block.
4. **Derivations written out: 30** (★ 5 · ★★ 23 · ★★★ 2), all A-block chains plus twelve results the book never writes
   out (rule (b): (4.25) Ex. 4.30, (4.37) add-and-subtract, (4.39a/b) coefficient bookkeeping, (4.40) Ex. 4.38, (4.43)
   Ex. 4.42, the Coriolis Ωut², the centrifugal potential Ex. 4.43, (4.55) Ex. 4.45, (4.56) Ex. 4.46, (4.58) Ex. 4.47,
   (4.68) Ex. 4.50, and ψ₂ − ψ₁ = flux Ex. 4.8). **27 demoted to statements** (§4c), including the two ★★★ surface-tension
   results (a-D50 Laplace from the cap, a-D51 meniscus) — both results are SEEN or examples, the book writes them out,
   and C14 states them with numeric checks (the analyst's "can stay B" note).
5. **Book slips handled in our own words** (analysis §9): (4.15)'s trailing "= 0" (N21); (4.51)'s dA for dV (N71);
   (4.74)'s gauge sign φ → φ − ∫B dt′ (N93, D26, a test that the printed sign doubles B); (4.63)'s "μ, κ, k" = μ, μ_v, k
   ≥ 0 (N81); the curve C's sign in §4.10 (N126); Ex. 4.7's dropped minus and γ = h/δ (N129); Ex. 4.2's ½ρU² statement
   (D06 trap); Cauchy prose ∂τ_ij/∂x_j vs the first index (D07 trap, non-symmetric τ test); "Section 3.6" → §3.4 (C07);
   (4.100) ↔ (4.101) and (4.106)/(4.107) → (4.109)/(4.112) references (N150); Ex. 4.8's rounded total (N156).
6. **Conventions stated where first used** (⚠️ callouts): tensor divergence contracts the **first** index (C06, D07);
   the prime has four meanings in one chapter — rotating frame (C09), perturbation (C13, N106), dummy variable (N87),
   ch03's translating frame (notation callout in C09); Coriolis *acceleration* +2Ω × u′ vs *force* −2Ω × u′, centripetal
   vs centrifugal (C09, D15); σ traceless only when μ_v = 0 or ∇·u = 0 (N41); mean vs thermodynamic pressure (N47);
   2-D ψ sign (u = ∂ψ/∂y; many GFD texts use u = −∂ψ/∂y) (C03); several meanings of σ, ε, γ, Φ, Ψ, Ω, h, M, R (notation
   table in the §4.1 notebook section); lapse rate/N² in Kundu's Γ ≡ dT/dz with the meteorological form alongside
   (N137 Ri, N114).
7. **Climate hooks, where real:** the Coriolis projectile at the pole and flow round highs and lows (C09); effective
   gravity and the geopotential (N64); dissipation ε and the kinetic-energy budget (C10 → Ch. 12, 13); the Boussinesq set
   as the ocean–atmosphere equations with the ocean's liquid-C_p argument (C13); internal Froude and Richardson numbers
   (N137); the stream function of geostrophic flow (C03 → Ch. 13); a Rossby-number mode in E9 as a forward pointer
   (defined there from the scale of the Coriolis term in (4.45), taught in Ch. 13).

## 1. Teaching order (A IDs grouped by book section, B/C IDs under each; one sentence each: "once you see X, Y follows")
A items in **bold**; B and C items listed where they are taught, with depth in brackets. Equations are written next to
their numbers so downstream agents have them.

**§4.1 Introduction** (own short notebook section; items tagged C01, C06)
- N01 [C → C01] integral (control-volume) vs differential (point) forms, equivalent via Gauss + localisation · N02
  [B → C06] the equations-vs-unknowns ledger, printed empty here and filled at C06 ($6<13$), C08 ($4<5$) and C10
  ($7=7$) · notation table (symbols with several meanings; primes).

**§4.2 Conservation of Mass**
- **C01 Mass for an arbitrarily moving control volume (4.5)** [#8]: once the mass of a material volume is constant,
  $\frac{d}{dt}\int_{V(t)}\rho\,dV=0$ *(Eq. 4.1)*, and a control volume V* is chosen to coincide with it at one instant
  (same region, same integrand ⇒ equal ∂ρ/∂t integrals and fluxes (4.4)), the RTT turns the material law into
  $\frac{d}{dt}\int_{V^*}\rho\,dV+\int_{A^*}\rho(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=0$ *(Eq. 4.5)* for any moving,
  deforming CV (D01) — the template reused for momentum, energy and the pillbox. Inside: R01 [B] material volume, N03 [B]
  (4.1) with the expanding-flow number, N04 [B] (4.2), N05 [B] (4.3), N06 [B] (4.4).
- **C02 The continuity equation (4.7)** [#11]: once Gauss turns the flux into $\int_V\nabla\cdot(\rho\mathbf u)dV$
  *(Eq. 4.6)* and the integral vanishes for *every* material volume, the localisation lemma forces
  $\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0$ *(Eq. 4.7)* (D02); the product rule gives
  $\frac1\rho\frac{D\rho}{Dt}+\nabla\cdot\mathbf u=0$ *(Eq. 4.8)* and incompressibility $D\rho/Dt=0$ *(Eq. 4.9)* ⇒
  $\nabla\cdot\mathbf u=0$ *(Eq. 4.10)* (D03). Inside: R02 [B] (4.6), N07 [B] localisation lemma, N08 [B] flux
  divergence, N09 [B] (4.8), N10 [B] (4.9) (stratified example), N11 [B] (4.10), N12 [B] incompressible flow vs fluid,
  M < 0.3.

**§4.3 Stream Functions**
- **C03 The 2-D stream function** [#21]: once steady continuity $\nabla\cdot(\rho\mathbf u)=0$ *(Eq. 4.11)* is satisfied
  identically by $\rho\mathbf u=\nabla\chi\times\nabla\psi$ *(Eq. 4.12)*, choosing χ = −z gives
  $\rho u=\partial\psi/\partial y,\ \rho v=-\partial\psi/\partial x$: ψ is constant along streamlines and the flux
  between two streamlines is $\psi_2-\psi_1$, so crowded contours mean fast flow (D04). Inside: N13 [B] (4.11), N14 [B]
  (4.12), N15 [B] two stream functions and stream surfaces (Fig. 4.1, plotly 3-D), N16 [B] stream-tube flux
  (b − a)(d − c), N17 [B] axisymmetric (Stokes) ψ, N18 [C] constant density.

**§4.4 Conservation of Momentum**
- **C04 Momentum for an arbitrarily moving control volume (4.17)** [#31]: once Newton's law for a material volume
  $\frac{d}{dt}\int_V\rho\mathbf u\,dV=\int_V\rho\mathbf g\,dV+\int_A\mathbf f\,dA$ *(Eq. 4.13)* goes through the same
  coincidence trick (D05), $\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV+\int_{A^*}\rho\mathbf u(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA$
  *(Eq. 4.17)* gives forces from fluxes — the wake behind a bar (Ex. 4.1), a bore (Ex. 4.3), a rocket (Ex. 4.4), a
  sprinkler's torque (§4.9). Inside: N19–N25 [B] (4.13)–(4.16d), N26 [B] body vs surface forces, N27 [B]
  $\mathbf g=-\nabla\Phi$ *(Eq. 4.18)*, R03 [B] $f_j=n_i\tau_{ij}$, R04 [C] surface tension via BCs, N28 [B] drag sign,
  N29 [B] Ex. 4.1 (worked number), N31 [B] Ex. 4.3, N32 [B] Ex. 4.4, N82 [C] (4.64), N83 [B] (4.65), N84 [B] Ex. 4.6
  (the §4.9 angular-momentum rows are taught here, as the rotational twin of (4.17)).
- **C05 Bernoulli along a streamline (4.19)** [#39]: once the momentum balance is applied to a thin stream-tube element
  (Ex. 4.2) and first-order terms are kept (D06), $\tfrac12U^2+gz+p/\rho$ is constant along a streamline *(Eq. 4.19)*
  for steady, inviscid, constant-density flow — a pitot tube reads speed from a pressure difference and a tank drains
  at $\sqrt{2gh}$. Inside: N30 [B] Ex. 4.2, N103 [B] pitot tube, N104 [B] stagnation and dynamic pressure, N105 [B]
  orifice (the §4.9 applications are taught with their equation).
- **C06 Cauchy's equation of motion (4.24)** [#47]: once Gauss converts both surface integrals of (4.14), contracting
  the stress's **first** index with n (D07), localisation gives the flux form
  $\frac{\partial}{\partial t}(\rho u_j)+\frac{\partial}{\partial x_i}(\rho u_iu_j)=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$
  *(Eq. 4.22)*, and subtracting $u_j\times$ continuity leaves
  $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ *(Eq. 4.24)* — Newton for every continuum, with
  6 equations for 13 unknowns. Inside: N02 [B] the ledger, N33–N37 [B] (4.20a)–(4.23).

**§4.5 Constitutive Equation for a Newtonian Fluid**
- **C07 The Newtonian stress law (4.31)** [#56]: once the stress is symmetric (N38, D08), splits as
  $\tau_{ij}=-p\delta_{ij}+\sigma_{ij}$ *(Eq. 4.27)* with σ depending only on S (R05), and σ is linear and isotropic,
  $K_{ijmn}=\lambda\delta_{ij}\delta_{mn}+\mu\delta_{im}\delta_{jn}+\gamma\delta_{in}\delta_{jm}$ *(Eq. 4.29)*, the δ
  substitutions collapse 81 coefficients to two:
  $\tau_{ij}=-p\,\delta_{ij}+2\mu S_{ij}+\lambda S_{mm}\delta_{ij}$ *(Eq. 4.31)* (D09 ★★★), rewritten with bulk viscosity
  as $\tau_{ij}=-p\delta_{ij}+2\mu(S_{ij}-\tfrac13S_{mm}\delta_{ij})+\mu_vS_{mm}\delta_{ij}$ *(Eq. 4.37)* (D10) — which
  is τ = μ du/dy (1.3) for a parallel flow. Inside: N38 [B] (4.25), N39 [C] constitutive equation, N40 [B] (4.26), N41
  [B] (4.27), R05 [B] σ(S), N42–N44 [B] (4.28)–(4.30), N45–N47 [B] (4.32)–(4.34), N48 [B] (4.35), N49 [B] μ_v, N50 [B]
  (4.36), N51 [B] (4.37), R06 [C] non-Newtonian.

**§4.6 Navier–Stokes Momentum Equation**
- **C08 Incompressible Navier–Stokes (4.39b)** [#67]: once the Newtonian stress (4.37) is put into Cauchy (4.24) (D11)
  we get (4.38); constant viscosity and ∇·u = 0 (D12) leave
  $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ *(Eq. 4.39b)*, whose viscous force is also
  $-\mu\nabla\times\boldsymbol\omega$ *(Eq. 4.40)* (D13) and whose μ = 0 limit is Euler (4.41). Inside: N52 [B] (4.38)
  + ledger 4 < 5, N53 [B] (4.39a), N54 [B] (4.40), N55 [B] (4.41).

**§4.7 Noninertial Frame of Reference**
- **C09 Navier–Stokes in a translating, rotating frame (4.45)** [#74]: once the rotating basis is seen to turn,
  $d\mathbf e'_i/dt=\boldsymbol\Omega\times\mathbf e'_i$, velocity picks up $\boldsymbol\Omega\times\mathbf x'$ (4.42)
  (D14) and differentiating once more gives five acceleration terms with the famous factor 2 (D15 ★★★); moved to the
  right of (4.39b) (D16) they become apparent body forces:
  $\rho\frac{D'\mathbf u'}{Dt}=-\nabla'p+\rho[\mathbf g-\frac{d\mathbf U}{dt}-2\boldsymbol\Omega\times\mathbf u'-\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')]+\mu\nabla'^2\mathbf u'$
  *(Eq. 4.45)* — Coriolis deflects a projectile by $\Omega ut^2$ (D17) and the centrifugal term folds into effective
  gravity (D18). Inside: N56 [C] inertial frame, N57 [B] (4.42), N58 [B] (4.43), N59 [B] (4.44), N60 [B] frame
  acceleration, N61 [B] Coriolis (worked number), N62 [B] highs and lows, N63 [C] angular acceleration, N64 [B]
  centrifugal and effective gravity, N65 [B] Ex. 4.5 pump.

**§4.8 Conservation of Energy**
- **C10 The internal-energy equation (4.57) and viscous dissipation** [#92]: once the first law is written for a
  material volume (4.46) and taken through RTT, Gauss and localisation to the total-energy equation (4.53) (D19, D20),
  subtracting the mechanical-energy equation — Cauchy dotted with u, (4.56) (D21) — leaves
  $\frac{De}{Dt}=-p\frac{Dv}{Dt}+\frac1\rho\sigma_{ij}S_{ij}-\frac1\rho\frac{\partial q_i}{\partial x_i}$ *(Eq. 4.57)*
  (D22); its middle term is the dissipation
  $\varepsilon=2\nu(S_{ij}-\tfrac13S_{mm}\delta_{ij})^2+\frac{\mu_v}\rho S_{mm}^2\ge0$ *(Eq. 4.58)* (D23), the one-way
  pipe from kinetic energy to heat — so entropy is produced, $\frac{k}{\rho T^2}\lvert\nabla T\rvert^2+\frac\varepsilon T\ge0$
  *(Eq. 4.63)*, and the system closes at 7 = 7. Inside: N66–N73 [B] (4.46)–(4.53), N74 [B] (4.54), N75 [B] (4.55), N76
  [B] (4.56), R07 [B] σ:R = 0, N77 [B] (4.58), N78 [B] (4.59), N79 [B] (4.60) + ledger, R08 [B] Gibbs (4.61), N80 [B]
  (4.62), N81 [B] (4.63).

**§4.9 Special Forms of the Equations**
- (angular momentum N82–N84 are taught inside C04; pitot, stagnation pressure and orifice N103–N105 inside C05.)
- **C11 The Bernoulli function: constant on streamlines and vortex lines (4.71)** [#109]: once the Lamb identity
  $u_i\frac{\partial u_j}{\partial x_i}=-(\mathbf u\times\boldsymbol\omega)_j+\frac{\partial}{\partial x_j}(\tfrac12u_i^2)$
  *(Eq. 4.68)* and a barotropic pressure function (4.67) turn Euler into $\partial\mathbf u/\partial t+\nabla B=\mathbf u\times\boldsymbol\omega$
  *(Eq. 4.69)* (D24), steady flow has $\nabla B\perp\mathbf u,\boldsymbol\omega$, so B is constant along streamlines and
  vortex lines *(Eq. 4.71)* and everywhere when ω = 0 *(Eq. 4.72)* (D25); the energy form $h+\tfrac12\lvert\mathbf u\rvert^2+gz$
  *(Eq. 4.78)* is a different statement with different hypotheses. Inside: N85 [C] not new laws, N86 [B] (4.66), N87 [B]
  (4.67), N88 [B] (4.68), N89 [B] (4.69), N90 [B] (4.70) Lamb surfaces, N91 [B] (4.72), N92 [C] irrotational stays
  irrotational (→ Ch. 5), N94–N96 [B] (4.76)–(4.78), N102 [B] the decision table of the four forms.
- **C12 Unsteady Bernoulli for potential flow (4.75)** [#114]: once u = ∇φ *(Eq. 4.73)* makes the Lamb term vanish,
  (4.69) is the gradient of one bracket, so the bracket depends on t only and can be absorbed into φ (with the sign the
  book gets wrong) (D26): $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int\frac{dp}{\rho}+gz$ =
  constant *(Eq. 4.75)* — the pressure of an accelerating flow, and Ch. 7's free-surface condition; with constant ρ and
  μ viscosity drops out of irrotational flow entirely (4.80) and a streamline form (4.82) follows. Inside: R09 [B]
  (4.73), N93 [B] (4.74), N97–N101 [B] (4.79)–(4.83).
- **C13 The Boussinesq approximation (4.86) + (4.89)** [#131]: once the hydrostatic state is subtracted (4.84) and the
  scalings show αδT ≪ 1 (N108), density matters only where it is multiplied by g (D27):
  $\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u$ *(Eq. 4.86)*, and the
  small pressure work of the expansion turns C_v into C_p (D28), leaving $\frac{DT}{Dt}=\kappa\nabla^2T$ *(Eq. 4.89)* —
  the ocean–atmosphere equation set. Inside: R10 [B] hydrostatic reference, N106 [B] (4.84), N107 [B] (4.85), N108 [B]
  validity conditions, N109 [B] (4.87), N110 [B] the C_p mechanism, N111 [B] (4.88), N112 [B] viscous heating
  negligible, N113 [B] (4.89) (Gaussian-blob animation), N114 [B] the Boussinesq set with the linear equation of state.

**§4.10 Boundary Conditions**
- **C14 The kinematic boundary condition (4.91), with the wall and interface conditions** [#142]: once a surface is
  written as η(x, t) = 0 and a point riding on it keeps η = 0 (4.90), "no fluid crosses it" means the fluid's normal
  velocity equals the surface's: $D\eta/Dt=0$ on η = 0 *(Eq. 4.91)* (D29); with the pillbox jump conditions (continuous
  mass flux, traction and heat flux), no-slip and no temperature jump, every later problem knows what happens at its
  edges; surface tension adds the Laplace jump, now derived from a force balance on a curved cap. Inside: N115 [C] BC
  specification, N116 [B] pillbox conditions (two worked cases), N117 [B] no-slip, N118 [B] (4.90), N119 [B] (4.92),
  N120 [B] (4.93), R11 [C] surface tension, N121–N123 [B] (4.94)–(4.96), N124 [C] Marangoni, N125 [B] (4.97), N126
  [B] (4.98), N127 [B] Laplace (1.5) derived, N128 [B] capillary length, N129 [B] Ex. 4.7 meniscus.

**§4.11 Dimensionless Forms of the Equations and Dynamic Similarity**
- **C15 Dimensionless Navier–Stokes (4.101) and dynamic similarity** [#160]: once every variable is scaled by the
  problem's own sizes (4.100) and (4.39b) is divided by ρU²/l (D30),
  $\big[\frac{\Omega l}{U}\big]\frac{\partial\mathbf u^*}{\partial t^*}+(\mathbf u^*\cdot\nabla^*)\mathbf u^*=-\nabla^*p^*+\big[\frac{gl}{U^2}\big]\mathbf g^*+\big[\frac{\mu}{\rho Ul}\big]\nabla^{*2}\mathbf u^*$
  *(Eq. 4.101)* shows that flows with equal St, Fr, Re and the same shape are the same flow — model tests work, sphere
  drag data collapse on one C_D(Re) curve, and a ship model cannot match Re and Fr at once. Inside: N130 [B] two routes,
  R12 [B] (4.99), N131 [B] Fig. 4.21, N132 [C] flow parameters, N133 [B] (4.100), N134–N136 [B] St, Re, Fr, N137 [B]
  Fr′, Ri, Ri_g, N138 [C], N139 [B] C_p, N140 [C], N141 [B] oscillatory body, N142–N143 [B] C_D, C_L, N144 [C],
  N145–N147 [B] (4.109)–(4.111) and M, N148–N152 [B] (4.112)–(4.116), N153–N155 [B] We, Bo, Ca, N156 [B] Ex. 4.8.

**Exercises and literature** — S01, S02 (one pointer line each at the end of the notebook).

## 2. Chapter map (depth) — every inventory row exactly once
One row per inventory row (185, `[#n]` = analysis §2 row). `A parent` names the A block a B or C item is written in.
IDs: CORE C01–C15, NOTE N01–N156 and RECAP R01–R12 in inventory order, SKIP S01–S02.

| ID | Item | § | Depth | Tier | A parent | Reason (A) / treatment (B) / pointer (C, SKIP) / source chapter (RECAP) |
|---|---|---|---|---|---|---|
| N01 | Integral form (a control volume and the fluxes through its surface) vs differential form (a point); the two are equivalent [#1] | 4.1 | C | NOTE | C01 | named in the chapter's opening paragraph and again at C02 (Gauss + localisation turns one into the other); the integral form is used in Ch. 9 (momentum integral), Ch. 13 (layer budgets), Ch. 15 (shock control volumes) |
| N02 | Equations vs unknowns: $6<13$ after Cauchy (4.24), $4<5$ after Navier–Stokes (4.38), $7=7$ after energy (4.60); thermal and caloric equations of state (1.12) supply two [#2] | 4.1 | B | NOTE | C06 | stated as a running ledger printed at C06, C08 and C10 (`ch04.closure_count(stage)` → a three-row table); the §4.1 notebook section introduces the ledger with an empty table that fills in |
| R01 | Material volume $V(t)$ and material surface $A(t)$ moving with the fluid, $\mathbf b=\mathbf u$ (sealed-balloon picture) [#3] | 4.2 | B | RECAP | C01 | ch03 N48/N53 (C15): reminded with `core.transport.ControlVolume` and a material interval tracked by `core.kinematics.pathline` in the expanding test flow |
| N03 | Eq. (4.1): $\frac{d}{dt}\int_{V(t)}\rho(\mathbf x,t)\,dV=0$ — the mass of a material volume is constant [#4] | 4.2 | B | NOTE | C01 | stated as D01's starting line; number: expanding flow $u=ax/(1+at)$, $\rho=\rho_0/(1+at)$ with a = 1 s⁻¹, ρ₀ = 1 kg/m³: the material interval [1, 2] m at t = 0 becomes [2, 4] m at t = 1 s with ρ = 0.5 kg/m³ — mass 1 kg/m² both times (`ch04.material_mass`) |
| N04 | Eq. (4.2): $\int_{V(t)}\frac{\partial\rho}{\partial t}dV+\int_{A(t)}\rho\,\mathbf u\cdot\mathbf n\,dA=0$ (RTT with F = ρ, b = u) [#5] | 4.2 | B | NOTE | C01 | stated as step 2 of D01 (the RTT (3.35) with F = ρ and b = u, recap of ch03 C15) |
| N05 | Eq. (4.3): $\frac{d}{dt}\int_{V^*(t)}\rho\,dV-\int_{V^*(t)}\frac{\partial\rho}{\partial t}dV-\int_{A^*(t)}\rho\,\mathbf b\cdot\mathbf n\,dA=0$ (RTT for an arbitrary CV) [#6] | 4.2 | B | NOTE | C01 | stated as step 3 of D01 |
| N06 | Eq. (4.4): at the coincidence instant $\int_{V^*}\frac{\partial\rho}{\partial t}dV=\int_{V}\frac{\partial\rho}{\partial t}dV=-\int_{A}\rho\mathbf u\cdot\mathbf n\,dA=-\int_{A^*}\rho\mathbf u\cdot\mathbf n\,dA$ — the instantaneously coincident control volume [#7] | 4.2 | B | NOTE | C01 | stated as steps 4–5 of D01 with the gap the book skips (same region, same integrand, same instant ⇒ equal integrals; the d/dt∫ρ terms differ because the two volumes part company a moment later); animation frame in E1 |
| C01 | Eq. (4.5): mass conservation for an arbitrarily moving control volume, $\frac{d}{dt}\int_{V^*(t)}\rho\,dV+\int_{A^*(t)}\rho(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=0$ [#8] | 4.2 | A | CORE | – | load-bearing: the template for every integral budget (momentum (4.17), energy (4.48), the pillbox of §4.10) and for Ch. 9 (momentum integral), Ch. 13 (layer mass budgets), Ch. 14–15 (CV lift, drag, shocks); b = u gives (4.1), b = 0 a fixed CV |
| R02 | Eq. (4.6): Gauss' theorem turns the mass flux into a volume integral, $\int_V\frac{\partial\rho}{\partial t}dV+\int_A\rho\mathbf u\cdot\mathbf n\,dA=\int_V\{\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)\}dV=0$ [#9] | 4.2 | B | RECAP | C02 | ch02 C14 (Gauss (2.30), D25): reminded with `core.integral_theorems.divergence_theorem_box` on ρu of the expanding flow; step 3 of D02 |
| N07 | Localisation lemma: $\int_Vf\,dV=0$ for every V and f continuous ⇒ $f\equiv0$ [#10] | 4.2 | B | NOTE | C02 | stated with its proof as steps 5–7 of D02 (the book gives one sentence): a positive f(x₀) stays > f(x₀)/2 on a small ball whose integral is then positive; demo: integrals of a bump over shrinking balls ∝ f(x₀)·volume (`scipy.integrate.tplquad`); reused in C06, C10 |
| C02 | Eq. (4.7): the continuity equation $\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0$, or $\frac{\partial\rho}{\partial t}+\frac{\partial}{\partial x_i}(\rho u_i)=0$ [#11] | 4.2 | A | CORE | – | load-bearing: one of the four field equations of every later chapter; every exact solution of Ch. 5–16 is checked against it, and ∇·u = 0 (4.10) is the constraint of incompressible flow (Ch. 6 potential flow, Ch. 8, Ch. 13) |
| N08 | Flux-divergence (transport) term $\nabla\cdot(\rho\mathbf u)$: net local loss; integrates to zero over a domain with no boundary flux [#12] | 4.2 | B | NOTE | C02 | stated with a periodic-box demo (`core.operators.divergence(bc='periodic')` sums to round-off; ρ rises where ∇·(ρu) < 0); one number printed |
| N09 | Eq. (4.8): $\frac{1}{\rho}\frac{D\rho}{Dt}+\nabla\cdot\mathbf u=0$ (product rule on ∇·(ρu)) [#13] | 4.2 | B | NOTE | C02 | stated as D03's result; number: expanding flow, (1/ρ)Dρ/Dt = −a/(1 + at) = −∇·u (= −1 s⁻¹ at t = 0); `core.navier_stokes.continuity_terms` |
| N10 | Eq. (4.9): incompressible flow $\frac{D\rho}{Dt}\equiv\frac{\partial\rho}{\partial t}+\mathbf u\cdot\nabla\rho=0$; re-displayed in §4.11 as $\nabla\cdot\mathbf u=-\frac1\rho\frac{D\rho}{Dt}=-\frac{1}{\rho c^2}\frac{Dp}{Dt}$ [#14] | 4.2 | B | NOTE | C02 | stated with the stratified shear flow u = (U(z), 0, 0), ρ = ρ(z): Dρ/Dt = 0 although ρ varies (`ch04.stratified_shear_flow`); ⚠️ the §4.11 re-display is a general identity, not the incompressible condition (taught at N146) |
| N11 | Eq. (4.10): $\nabla\cdot\mathbf u=0$ for incompressible flow [#15] | 4.2 | B | NOTE | C02 | stated as D03's last step with the ch03 C09 reading (zero volumetric strain rate); `core.operators.is_solenoidal`, `core.navier_stokes.divergence_free_check` on the cylinder flow (0 to round-off) |
| N12 | Constant-density flow ⊂ incompressible flow; liquids nearly incompressible; gases incompressible at low Mach number (M < 0.3) [#16] | 4.2 | B | NOTE | C02 | stated with a two-row table and one number: air at 288 K, U = 100 m/s → M ≈ 0.29 (`core.similarity.mach_number`, `ch04.is_incompressible_regime`); ⚠️ incompressible flow vs incompressible fluid; threshold from the book stays in the JSON |
| N13 | Eq. (4.11): steady continuity $\nabla\cdot(\rho\mathbf u)=0$ [#17] | 4.3 | B | NOTE | C03 | stated as C03's starting line (C02 with ∂/∂t = 0) |
| N14 | Eq. (4.12): $\rho\mathbf u=\nabla\times\boldsymbol\Psi$ with $\boldsymbol\Psi=\chi\nabla\psi$, so $\rho\mathbf u=\nabla\chi\times\nabla\psi$ [#18] | 4.3 | B | NOTE | C03 | stated (∇·∇× ≡ 0 and ∇×∇ ≡ 0 recalled from ch02 D26 corollary); step 1–2 of D04; sympy check `core.streamfunction.mass_flux_from_stream_functions` |
| N15 | Two stream functions χ, ψ: $\rho\mathbf u\cdot\nabla\chi=\rho\mathbf u\cdot\nabla\psi=0$; 3-D streamlines are intersections of stream surfaces (Fig. 4.1) [#19] | 4.3 | B | NOTE | C03 | stated with a plotly 3-D figure of the two surface families and traced streamlines on their intersections (`core.streamfunction.stream_surface_check`); result of a-D05 given (§4c) |
| N16 | Stream-tube mass flux by Stokes: $\dot m=\oint_C\chi\,d\psi=(b-a)(d-c)$ [#20] | 4.3 | B | NOTE | C03 | stated with the leg-by-leg picture (on ψ = const legs dψ = 0) and a numeric patch flux (`core.streamfunction.stream_tube_mass_flux`) equal to (b − a)(d − c); a-D06 given (§4c) |
| C03 | 2-D stream function (χ = −z): $\rho u=\partial\psi/\partial y$, $\rho v=-\partial\psi/\partial x$; ψ constant on streamlines and $\psi_2-\psi_1$ = flux between them [#21] | 4.3 | A | CORE | – | load-bearing: the tool of every 2-D flow picture that follows — potential flows (Ch. 6), Couette/Poiseuille (Ch. 8), Blasius (Ch. 9), geostrophic and quasi-geostrophic flow (Ch. 13); contour spacing = speed |
| N17 | Axisymmetric (Stokes) stream function (χ = −φ): $\rho u_R=-R^{-1}\partial\psi/\partial z$, $\rho u_z=R^{-1}\partial\psi/\partial R$ [#22] | 4.3 | B | NOTE | C03 | stated with one test field (uniform stream ψ = ½UR² → u_z = U) and `core.streamfunction.velocity_from_streamfunction_axisym`; used for the sphere in Ch. 6 and pipe flow in Ch. 8 |
| N18 | Constant density: the same construction holds for u, and ψ-differences are volume fluxes [#23] | 4.3 | C | NOTE | C03 | named in one sentence (all functions take `rho=1.0`); used throughout Ch. 6 |
| N19 | Eq. (4.13): Newton's second law for a material volume, $\frac{d}{dt}\int_{V(t)}\rho\mathbf u\,dV=\int_{V(t)}\rho\mathbf g\,dV+\int_{A(t)}\mathbf f(\mathbf n,\mathbf x,t)\,dA$ [#24] | 4.4 | B | NOTE | C04 | stated as D05's starting line (Newton II recap P09 in words: momentum per volume ρu, body force per mass g, surface force per area f) |
| N20 | Eq. (4.14): $\int_{V}\frac{\partial}{\partial t}(\rho\mathbf u)dV+\int_{A}\rho\mathbf u(\mathbf u\cdot\mathbf n)dA=\int_{V}\rho\mathbf g\,dV+\int_{A}\mathbf f\,dA$ (RTT, component by component) [#25] | 4.4 | B | NOTE | C04 | stated as step 2 of D05; the start of C06's D07 as well |
| N21 | Eq. (4.15): $\int_{V^*}\frac{\partial}{\partial t}(\rho\mathbf u)dV=\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV-\int_{A^*}\rho\mathbf u\,\mathbf b\cdot\mathbf n\,dA$ [#26] | 4.4 | B | NOTE | C04 | stated as step 3 of D05; ⚠️ book prints a spurious trailing “= 0” (analysis §9.1), taught as “book prints X; the identity needs Y” |
| N22 | Eq. (4.16a): $\int_{V(t)}\frac{\partial}{\partial t}(\rho\mathbf u)dV=\int_{V^*(t)}\frac{\partial}{\partial t}(\rho\mathbf u)dV$ [#27] | 4.4 | B | NOTE | C04 | stated as step 4 of D05 (the coincidence trick of D01, now for a vector) |
| N23 | Eq. (4.16b): $\int_{A(t)}\rho\mathbf u(\mathbf u\cdot\mathbf n)dA=\int_{A^*(t)}\rho\mathbf u(\mathbf u\cdot\mathbf n)dA$ [#28] | 4.4 | B | NOTE | C04 | stated as step 4 of D05 |
| N24 | Eq. (4.16c): $\int_{V(t)}\rho\mathbf g\,dV=\int_{V^*(t)}\rho\mathbf g\,dV$ [#29] | 4.4 | B | NOTE | C04 | stated as step 4 of D05 |
| N25 | Eq. (4.16d): $\int_{A(t)}\mathbf f\,dA=\int_{A^*(t)}\mathbf f\,dA$ [#30] | 4.4 | B | NOTE | C04 | stated as step 4 of D05 |
| C04 | Eq. (4.17): momentum conservation for an arbitrarily moving control volume, $\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV+\int_{A^*}\rho\mathbf u(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA$ [#31] | 4.4 | A | CORE | – | load-bearing: forces from fluxes without knowing the flow inside — wake drag, bores, thrust (this chapter), the von Kármán momentum integral (Ch. 9), lift and drag from a CV (Ch. 14), shock relations (Ch. 15) |
| N26 | Body forces (no contact, ∝ mass, per unit mass an acceleration; fictitious ones in noninertial frames) vs surface forces (contact, ∝ area, the stress) [#32] | 4.4 | B | NOTE | C04 | stated in a two-column table (per mass vs per area; gravity vs pressure/viscous) at the head of C04 |
| N27 | Eq. (4.18): conservative body force $\mathbf g=-\nabla\Phi$ or $g_j=-\partial\Phi/\partial x_j$; gravity Φ = gz with z up [#33] | 4.4 | B | NOTE | C04 | stated with the path-independence check (work of g round a closed loop = 0, `core.integral_theorems.circulation`) and `ch04.body_force_from_potential`; reused by C11 |
| R03 | Surface force from the stress tensor $f_j=n_i\tau_{ij}$ (2.15); normal part $n_if_i$, tangential $f_k-(n_if_i)n_k$ [#34] | 4.4 | B | RECAP | C04 | ch02 C05 (Cauchy's traction, D05, first index = face normal): reminded with `core.tensors.traction`, `normal_shear_stress` |
| R04 | Surface tension acts on lines within interfaces and enters through boundary conditions, not the field equations [#35] | 4.4 | C | RECAP | C04 | ch01 C15 (surface tension): one reminder sentence with a pointer to §4.10 (C14, the cap force balance) |
| N28 | Drag sign convention: drag is the force on the object; the CV law needs the force on the fluid, $-F_D\mathbf e_x$ (Newton III) [#36] | 4.4 | B | NOTE | C04 | stated as a ⚠️ Common confusion callout in C04's worked example; the docstring of `ch04.wake_drag_per_span` repeats it |
| N29 | Ex. 4.1: drag of a long bar from the wake deficit, $F_D/l=\rho\int_{-H/2}^{+H/2}U(y)\big(U_\infty-U(y)\big)dy$ (Fig. 4.2) [#37] | 4.4 | B | NOTE | C04 | C04's worked number (result given, a-D09 in §4c, face-by-face fluxes shown as bars): Gaussian deficit U∞ = 10 m/s, Δ = 2 m/s, b = 0.1 m, ρ = 1.2 kg/m³ → side leakage 0.354 m²/s and F_D/l ≈ 3.65 N/m (`ch04.wake_drag_per_span`, `gaussian_wake`); momentum thickness → Ch. 9; E1 mode |
| N30 | Ex. 4.2: Bernoulli from a stream-tube element (Fig. 4.3): first-order mass and streamwise momentum balance, $U\frac{\partial U}{\partial s}ds=-g\,dz-\frac1\rho\frac{\partial p}{\partial s}ds$ [#38] | 4.4 | B | NOTE | C05 | stated as the construction of D06 (C05's derivation, with the conical-side pressure force and the dropped (ds)² terms the book skips); sympy re-run `ch04.stream_tube_element_balance_sym` |
| C05 | Eq. (4.19): Bernoulli, $\tfrac12U^2+gz+p/\rho$ = constant along a streamline (steady, inviscid, constant density) [#39] | 4.4 | A | CORE | – | load-bearing: the most-used result of the book (pitot tubes, orifices, pressure on bodies in Ch. 6 and 14, the outer flow of Ch. 9, the free surface of Ch. 7); the constant-ρ case of (4.71) |
| N31 | Ex. 4.3: small wave (bore) with a moving CV, $U=\sqrt{\frac{gh_{out}}{2h_{in}}(h_{in}+h_{out})}\approx\sqrt{gh}$ (Fig. 4.4) [#40] | 4.4 | B | NOTE | C04 | stated with a-D11 given (§4c) and one number: h_in = 1 m, h_out = 1.1 m → U = 3.37 m/s vs √(g·1 m) = 3.13 m/s (`ch04.bore_speed`); E1 mode (moving CV b = U e_x); shallow-water speed → Ch. 7, 13 |
| N32 | Ex. 4.4: rocket with an accelerating CV, $\frac{dM}{dt}+\rho_eV_eA_e=0$, $M\frac{d^2z_R}{dt^2}=-V_e\frac{dM}{dt}-Mg+F_S$ (Fig. 4.5) [#41] | 4.4 | B | NOTE | C04 | stated with a-D12 given (§4c); number: Tsiolkovsky Δb = V_e ln(M₀/M₁) when g = F_S = 0 (`ch04.rocket_trajectory`, `rocket_delta_v`); E1 mode (accelerating CV) |
| N33 | Eq. (4.20a): $\int_A\rho\mathbf u(\mathbf u\cdot\mathbf n)dA=\int_V\frac{\partial}{\partial x_i}(\rho u_iu_j)dV$ [#42] | 4.4 | B | NOTE | C06 | stated as step 2 of D07 (Gauss per component j, the contracted index is the normal's) |
| N34 | Eq. (4.20b): $\int_A n_i\tau_{ij}\,dA=\int_V\frac{\partial\tau_{ij}}{\partial x_i}dV$ (first index) [#43] | 4.4 | B | NOTE | C06 | stated as step 3 of D07; ⚠️ first-index divergence (`core.operators.tensor_divergence(index=0)`), a non-symmetric test τ pins it |
| N35 | Eq. (4.21): $\int_V\{\frac{\partial}{\partial t}(\rho u_j)+\frac{\partial}{\partial x_i}(\rho u_iu_j)-\rho g_j-\frac{\partial\tau_{ij}}{\partial x_i}\}dV=0$ [#44] | 4.4 | B | NOTE | C06 | stated as step 4 of D07 |
| N36 | Eq. (4.22): conservative (flux) form $\frac{\partial}{\partial t}(\rho u_j)+\frac{\partial}{\partial x_i}(\rho u_iu_j)=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ [#45] | 4.4 | B | NOTE | C06 | stated as step 5 of D07 (localisation, N07); the finite-volume form of Ch. 10; `core.navier_stokes.momentum_conservative_residual` |
| N37 | Eq. (4.23): the flux form expanded, the bracket is continuity: $\ldots=\rho\frac{Du_j}{Dt}$ [#46] | 4.4 | B | NOTE | C06 | stated as steps 6–8 of D07; sympy identity `core.navier_stokes.conservative_to_advective_sym` (difference = u_j × continuity) |
| C06 | Eq. (4.24): Cauchy's equation of motion $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ [#47] | 4.4 | A | CORE | – | load-bearing: Newton's law for any continuum — Navier–Stokes (4.38) is Cauchy plus a stress law, and the same template (material volume → RTT → Gauss → localise) gives the energy equation; ledger 6 equations vs 13 unknowns |
| N38 | Eq. (4.25): the stress tensor is symmetric, $\tau_{ij}=\tau_{ji}$ (6 independent components; broken only by body couples) [#48] | 4.5 | B | NOTE | C07 | stated with D08 written out (book defers it to Exercise 4.30): torque on a shrinking cube, angular acceleration 6(τ₁₂ − τ₂₁)/(ρh²) → ∞ as h → 0; `ch04.cube_spin_acceleration` log–log slope −2 |
| N39 | Constitutive equation; the Newtonian fluid is the simplest linear stress–strain-rate law [#49] | 4.5 | C | NOTE | C07 | named in C07's opening scene; non-Newtonian fluids → Ch. 16 |
| N40 | Eq. (4.26): fluid at rest, isotropic stress $\tau_{ij}=-p\,\delta_{ij}$ (p the thermodynamic pressure; minus because tension is positive) [#50] | 4.5 | B | NOTE | C07 | stated (δ is the only isotropic 2nd-order tensor, ch02 N42; a-D16 given in §4c); traction on any n is −pn (`core.constitutive.static_stress` + `core.tensors.traction`) |
| N41 | Eq. (4.27): $\tau_{ij}=-p\,\delta_{ij}+\sigma_{ij}$ (σ the viscous or “deviatoric” stress) [#51] | 4.5 | B | NOTE | C07 | stated; ⚠️ σ is traceless only when μ_v = 0 or ∇·u = 0 (tr σ = 3μ_v∇·u); `core.constitutive.total_stress` |
| R05 | Galilean invariance ⇒ σ depends on ∇u not u; only S enters, R (rigid rotation) does not [#52] | 4.5 | B | RECAP | C07 | ch03 C05, C11, D10 (rigid motion ⇒ S = 0): reminded with `ch03.rigid_body_velocity` → σ = 0 in `core.constitutive.viscous_stress` |
| N42 | Eq. (4.28): most general linear law $\sigma_{ij}=K_{ijmn}S_{mn}$ (81 coefficients) [#53] | 4.5 | B | NOTE | C07 | stated as step 1 of D09; `core.constitutive.linear_stress` (`np.einsum('ijmn,mn->ij')`) |
| N43 | Eq. (4.29): isotropic 4th-order tensor $K_{ijmn}=\lambda\delta_{ij}\delta_{mn}+\mu\delta_{im}\delta_{jn}+\gamma\delta_{in}\delta_{jm}$ [#54] | 4.5 | B | NOTE | C07 | stated (cited, Aris) as step 3 of D09 with a numeric check instead of a proof: invariant under 50 random rotations (`core.constitutive.isotropic_fourth_order`, `core.tensors.transform_tensor`) |
| N44 | Eq. (4.30): symmetry requires $\gamma=\mu$ [#55] | 4.5 | B | NOTE | C07 | stated as steps 6–7 of D09 with the subtlety taught: on a symmetric S only μ + γ acts, so “γ = μ” names that sum 2μ (sympy) |
| C07 | Eq. (4.31): Newtonian stress $\tau_{ij}=-p\,\delta_{ij}+2\mu S_{ij}+\lambda S_{mm}\delta_{ij}$ (from linearity, isotropy, symmetry) [#56] | 4.5 | A | CORE | – | load-bearing: the stress law that turns Cauchy into Navier–Stokes; wall shear and heating in Ch. 8, τ_w in Ch. 9, CFD fluxes in Ch. 10, the Reynolds-stress analogy in Ch. 12; generalises τ = μ du/dy (1.3) |
| N45 | Eq. (4.32): $p=-\tfrac13\tau_{ii}+\big(\tfrac23\mu+\lambda\big)\nabla\cdot\mathbf u$ [#57] | 4.5 | B | NOTE | C07 | stated (trace of (4.31), δ_ii = 3; a-D18 given in §4c); `core.constitutive.thermodynamic_pressure_from_stress` |
| N46 | Eq. (4.33): mean (mechanical) pressure $\bar p\equiv-\tfrac13\tau_{ii}$ [#58] | 4.5 | B | NOTE | C07 | stated; `core.constitutive.mean_pressure` |
| N47 | Eq. (4.34): $p-\bar p=\big(\tfrac23\mu+\lambda\big)\nabla\cdot\mathbf u$; incompressible p is only mechanical, known up to a constant [#59] | 4.5 | B | NOTE | C07 | stated with a number (μ_v = 0 ⇒ p = p̄) and the demo that adding a constant to p leaves (4.39b) unchanged; ⚠️ mean vs thermodynamic pressure callout |
| N48 | Eq. (4.35): incompressible Newtonian stress $\tau_{ij}=-p\,\delta_{ij}+2\mu S_{ij}$ [#60] | 4.5 | B | NOTE | C07 | stated; `core.constitutive.newtonian_stress(..., incompressible=True)` raises if tr G ≠ 0 |
| N49 | Bulk viscosity $\mu_v=\lambda+\tfrac23\mu$ (sound absorption, shock structure; non-zero in polyatomic gases) [#61] | 4.5 | B | NOTE | C07 | stated; `core.constitutive.bulk_viscosity`; → Ch. 15 |
| N50 | Eq. (4.36): Stokes assumption $\lambda+\tfrac23\mu=0$ [#62] | 4.5 | B | NOTE | C07 | stated as the code default `mu_v=0.0`; → Ch. 15 |
| N51 | Eq. (4.37): $\tau_{ij}=-p\,\delta_{ij}+2\mu\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)+\mu_vS_{mm}\delta_{ij}$; parallel flow gives τ = μ du/dy (1.3) [#63] | 4.5 | B | NOTE | C07 | stated with D10 written out (the book states it; add and subtract (2μ/3)S_mm δ_ij); number: u = (γy, 0, 0), γ = 10 s⁻¹, water μ = 1.0×10⁻³ Pa s → τ₁₂ = 0.01 Pa = `ch01.newton_shear_stress` |
| R06 | Non-Newtonian fluids (shear thinning, memory, viscoelasticity) [#64] | 4.5 | C | RECAP | C07 | ch01 C02 (Bingham and Maxwell materials named): one sentence, pointer to Ch. 16 |
| N52 | Eq. (4.38): Navier–Stokes with variable μ, μ_v, $\rho\big(\frac{\partial u_j}{\partial t}+u_i\frac{\partial u_j}{\partial x_i}\big)=-\frac{\partial p}{\partial x_j}+\rho g_j+\frac{\partial}{\partial x_i}\big[\mu\big(\frac{\partial u_j}{\partial x_i}+\frac{\partial u_i}{\partial x_j}\big)+\big(\mu_v-\tfrac23\mu\big)\frac{\partial u_m}{\partial x_m}\delta_{ij}\big]$; 4 equations, 5 unknowns; barotropic closure [#65] | 4.6 | B | NOTE | C08 | stated as D11's result (written out: C08's derivation chain); ledger 4 vs 5; `core.navier_stokes.navier_stokes_residual`, `navier_stokes_sym` |
| N53 | Eq. (4.39a): constant μ, μ_v (compressible), $\rho\frac{Du_j}{Dt}=-\frac{\partial p}{\partial x_j}+\rho g_j+\mu\frac{\partial^2u_j}{\partial x_i^2}+\big(\mu_v+\tfrac13\mu\big)\frac{\partial}{\partial x_j}\frac{\partial u_m}{\partial x_m}$ [#66] | 4.6 | B | NOTE | C08 | stated as D12's middle result (Schwarz swap, coefficient bookkeeping μ + μ_v − ⅔μ); sympy identity |
| C08 | Eq. (4.39b): incompressible Navier–Stokes $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ [#67] | 4.6 | A | CORE | – | load-bearing: the equation of the book — vorticity (Ch. 5), potential flow (Ch. 6, as Euler), waves (Ch. 7), every exact solution (Ch. 8), boundary layers (Ch. 9), CFD (Ch. 10), instability (Ch. 11), turbulence (Ch. 12), GFD (Ch. 13) |
| N54 | Eq. (4.40): viscous force three ways, $\mu\nabla^2\mathbf u=2\mu\,\partial S_{ij}/\partial x_i=-\mu\nabla\times\boldsymbol\omega$ (incompressible); solid-body rotation feels none [#68] | 4.6 | B | NOTE | C08 | stated with D13 written out (book defers it to Exercise 4.38; ε–δ recap ch02 D09) and the resolved “paradox” in numbers: solid-body rotation has uniform ω ⇒ zero viscous force (`core.navier_stokes.viscous_force_forms`); → Ch. 5 vorticity diffusion |
| N55 | Eq. (4.41): Euler equation $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g$ (viscosity negligible away from boundaries) [#69] | 4.6 | B | NOTE | C08 | stated (μ = 0 in (4.39b)) with the check that cylinder potential flow + Bernoulli pressure satisfies it (`core.navier_stokes.ns_incompressible_terms(mu=0)`); → Ch. 5, 6, 7, 14; E4 preset |
| N56 | Inertial frame (fixed to the distant stars; the lab usually suffices; the Earth frame is noninertial on geophysical scales) [#70] | 4.7 | C | NOTE | C09 | named in C09's opening scene (a turntable and the rotating Earth); used throughout Ch. 13 |
| N57 | Eq. (4.42): $\mathbf u=\mathbf U+\mathbf u'+\boldsymbol\Omega\times\mathbf x'$, with $d\mathbf e'_i/dt=\boldsymbol\Omega\times\mathbf e'_i$ (Figs. 4.6–4.7) [#71] | 4.7 | B | NOTE | C09 | stated as D14's result (the cone construction and the product rule on x′_i e′_i); `core.rotating.rotating_basis`, `inertial_velocity`; numeric d e′/dt = Ω × e′ |
| N58 | Eq. (4.43): $\mathbf a=\frac{d\mathbf U}{dt}+\mathbf a'+2\boldsymbol\Omega\times\mathbf u'+\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'+\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')$ [#72] | 4.7 | B | NOTE | C09 | stated with D15 written out (★★★; the book states the result and defers it to Exercise 4.42): why the Coriolis term has a 2 (one Ω × u′ from the turning basis, one from d(Ω × x′)/dt); `core.rotating.frame_acceleration_terms` vs finite differences |
| N59 | Eq. (4.44): $\big(\frac{D\mathbf u}{Dt}\big)_{O123}=\big(\frac{D'\mathbf u'}{Dt}\big)_{O'1'2'3'}+\frac{d\mathbf U}{dt}+2\boldsymbol\Omega\times\mathbf u'+\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'+\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')$ [#73] | 4.7 | B | NOTE | C09 | stated as D15's last step (particle derivative → D/Dt in each frame) |
| C09 | Eq. (4.45): incompressible Navier–Stokes in a noninertial frame, $\rho\frac{D'\mathbf u'}{Dt}=-\nabla'p+\rho\big[\mathbf g-\frac{d\mathbf U}{dt}-2\boldsymbol\Omega\times\mathbf u'-\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')\big]+\mu\nabla'^2\mathbf u'$ [#74] | 4.7 | A | CORE | – | load-bearing: the frame of every geophysical flow — Coriolis and effective gravity are the starting point of Ch. 13 (f-plane, geostrophy, Ekman layers, Rossby waves) and of rotating flows in Ch. 8; the climate-dynamics entry point |
| N60 | Frame-acceleration term $-d\mathbf U/dt$ (seat push; weightlessness in parabolic flight when dU/dt = g) [#75] | 4.7 | B | NOTE | C09 | stated as one bracket term with one number (dU/dt = g cancels g exactly); E5 preset |
| N61 | Coriolis force per mass $-2\boldsymbol\Omega\times\mathbf u'$: velocity-dependent, deflects right in the NH, does no work; projectile from the pole: deflection $\Omega ut^2$, angle Ωt (Fig. 4.8) [#76] | 4.7 | B | NOTE | C09 | stated with D17 written out (the book asserts Ωut² without the ½(2Ωu)t² step) and C09's worked number: Ω = 7.292×10⁻⁵ rad/s, u = 10 m/s, t = 1 h → acceleration 1.46×10⁻³ m/s², forward 36 km, deflection 9.45 km, angle 15.0°; `ch04.coriolis_projectile`; ⚠️ acceleration +2Ω × u′ (4.43) vs force −2Ω × u′ (4.45) |
| N62 | Flow out of a high and into a low: $-2\boldsymbol\Omega\times\mathbf u=-2\Omega_zu_R\,\mathbf e_\varphi$ — clockwise round highs, counter-clockwise round lows (NH) [#77] | 4.7 | B | NOTE | C09 | stated with a quiver figure (radial outflow + its Coriolis arrows, `core.rotating.coriolis_acceleration`); E5 mode; cyclones and gradient wind → Ch. 13 |
| N63 | Angular-acceleration term $-(d\boldsymbol\Omega/dt)\times\mathbf x'$, unimportant for steady rotation [#78] | 4.7 | C | NOTE | C09 | named as one bar of E5's term list (zero for the Earth) |
| N64 | Centrifugal term $-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')=+\Omega^2R\,\mathbf e_R$ with potential $-\tfrac12\Omega^2R^2$; effective gravity $\mathbf g_e=\mathbf g+\Omega^2R\,\mathbf e_R$, equipotentials and the Earth's bulge (Fig. 4.9) [#79] | 4.7 | B | NOTE | C09 | stated with D18 written out (Exercise 4.43: triple product, then the gradient of −½Ω²R²); number: Ω²a = 0.0339 m/s² at the equator; `core.rotating.effective_gravity`; WGS-84 bulge 42.77 km (book's rounded value in the JSON); geopotential → Ch. 13 |
| N65 | Ex. 4.5: von Kármán viscous pump — the three cylindrical momentum equations with rotation terms $\rho[2\Omega_zu_\varphi+\Omega_z^2R]$ and $\rho[-2\Omega_zu_R]$ (Fig. 4.10) [#80] | 4.7 | B | NOTE | C09 | stated with the rotation terms highlighted in a sympy printout (`ch04.rotating_pump_equations`, `core.curvilinear`); a-D28 given in §4c; Appendix-B operators glossed; → Ch. 8 rotating flows |
| N66 | Eq. (4.46): energy of a material volume, $\frac{d}{dt}\int_V\rho\big(e+\tfrac12\lvert\mathbf u\rvert^2\big)dV=\int_V\rho\mathbf g\cdot\mathbf u\,dV+\int_A\mathbf f\cdot\mathbf u\,dA-\int_A\mathbf q\cdot\mathbf n\,dA$ [#81] | 4.8 | B | NOTE | C10 | stated as D19's starting line (first law (1.10) per material volume, recap ch01 C25: heat in + work done = rise of internal + kinetic energy) |
| N67 | Eq. (4.47): the RTT form of (4.46) [#82] | 4.8 | B | NOTE | C10 | stated as step 2 of D19 |
| N68 | Eq. (4.48): energy for an arbitrarily moving control volume (instantaneously coincident) [#83] | 4.8 | B | NOTE | C10 | stated (a-D29 given in §4c: the D01 coincidence algebra with F = ρ(e + ½u²)); `core.conservation.energy_budget`; → Ch. 15 stagnation enthalpy |
| N69 | Eq. (4.49): Gauss on the energy flux, $\int_A\big(\rho e+\tfrac\rho2\lvert\mathbf u\rvert^2\big)(\mathbf u\cdot\mathbf n)dA=\int_V\frac{\partial}{\partial x_i}\big(\rho(e+\tfrac12u_j^2)u_i\big)dV$ [#84] | 4.8 | B | NOTE | C10 | stated as step 3 of D19 |
| N70 | Eq. (4.50): $\int_A\mathbf f\cdot\mathbf u\,dA=\int_An_i\tau_{ij}u_j\,dA=\int_V\frac{\partial}{\partial x_i}(\tau_{ij}u_j)dV$ [#85] | 4.8 | B | NOTE | C10 | stated as step 4 of D19 |
| N71 | Eq. (4.51): $\int_A\mathbf q\cdot\mathbf n\,dA=\int_V\nabla\cdot\mathbf q\,dV$ [#86] | 4.8 | B | NOTE | C10 | stated as step 5 of D19; ⚠️ the book prints dA inside the volume integrals (must be dV) |
| N72 | Eq. (4.52): all terms under one volume integral [#87] | 4.8 | B | NOTE | C10 | stated as step 6 of D19 |
| N73 | Eq. (4.53): total-energy equation $\frac{\partial}{\partial t}\big(\rho[e+\tfrac12u_j^2]\big)+\frac{\partial}{\partial x_i}\big(\rho[e+\tfrac12u_j^2]u_i\big)=\rho g_iu_i+\frac{\partial}{\partial x_i}(\tau_{ij}u_j)-\frac{\partial q_i}{\partial x_i}$ [#88] | 4.8 | B | NOTE | C10 | stated as D19's result; → Ch. 10 (conservative form), Ch. 15 |
| N74 | Eq. (4.54): stress work = deformation work $(-p\,\partial u_j/\partial x_j+\sigma_{ij}\partial u_j/\partial x_i)$ + work of the net forces $(-u_j\partial p/\partial x_j+u_j\partial\sigma_{ij}/\partial x_i)$ [#89] | 4.8 | B | NOTE | C10 | stated as step 2 of D20 (product rule); `core.navier_stokes.stress_work_split` sums back to ∂(τ_ij u_j)/∂x_i on a random field; the split behind the turbulent kinetic-energy budget of Ch. 12 |
| N75 | Eq. (4.55): total energy in D/Dt form [#90] | 4.8 | B | NOTE | C10 | stated with D20 written out (book defers it to Exercise 4.45: the (4.23) move with E = e + ½u_j²) |
| N76 | Eq. (4.56): mechanical (kinetic) energy $\rho\frac{D}{Dt}\big(\tfrac12u_j^2\big)=\rho g_ju_j-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}$ [#91] | 4.8 | B | NOTE | C10 | stated with D21 written out (book defers it to Exercise 4.46; multiply Cauchy (4.24) by u_j); `core.navier_stokes.kinetic_energy_budget`; mean and turbulent KE → Ch. 12 |
| C10 | Eq. (4.57): internal-energy equation $\frac{De}{Dt}=-p\frac{Dv}{Dt}+\frac1\rho\sigma_{ij}S_{ij}-\frac1\rho\frac{\partial q_i}{\partial x_i}$ — the thermal half of the energy budget, with the viscous dissipation that feeds it [#92] | 4.8 | A | CORE | – | load-bearing: splits energy into a mechanical part that viscosity drains and a thermal part that it heats — viscous heating in Ch. 8, thermal boundary layers in Ch. 9, the dissipation ε of turbulence (Ch. 12) and of the ocean–atmosphere energy cycle (Ch. 13) |
| R07 | $\sigma_{ij}\,\partial u_j/\partial x_i=\sigma_{ij}S_{ij}$ for symmetric σ [#93] | 4.8 | B | RECAP | C10 | ch02 N61 (τ_ij A_ij = 0 for symmetric τ, antisymmetric A — the ch02 pointer to exactly this step): reminded with `core.tensors.symmetric_double_contraction`; step of D22 |
| N77 | Eq. (4.58): dissipation rate $\varepsilon\equiv\frac1\rho\sigma_{ij}S_{ij}=2\nu\big(S_{ij}-\tfrac13\frac{\partial u_m}{\partial x_m}\delta_{ij}\big)^2+\frac{\mu_v}{\rho}\big(\frac{\partial u_m}{\partial x_m}\big)^2\ge0$ [#94] | 4.8 | B | NOTE | C10 | stated with D23 written out (book defers it to Exercise 4.47: completing the square); C10's worked number: plane Couette U = 1 m/s, h = 1 mm, water → ρε = μ(U/h)² = 1000 W/m³; `core.constitutive.dissipation_rate(form='both')`; ε of turbulence → Ch. 12 |
| N78 | Eq. (4.59): Newtonian viscous stress $\sigma_{ij}=\mu\big(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\big)+\big(\mu_v-\tfrac23\mu\big)\frac{\partial u_m}{\partial x_m}\delta_{ij}$ [#95] | 4.8 | B | NOTE | C10 | stated (= the σ part of (4.37), C07); `core.constitutive.viscous_stress` |
| N79 | Eq. (4.60): internal energy with Fourier conduction; field equations (4.7), (4.38), (4.60) + two state equations: 7 equations, 7 unknowns [#96] | 4.8 | B | NOTE | C10 | stated (substitute (4.58) and q = −k∇T, recap ch01 C11); ledger 7 = 7; number: Couette with viscous heating, both walls at T₀, U = 1 m/s, h = 1 mm, water k = 0.6 W/(m K) → peak rise μU²/(8k) ≈ 2.1×10⁻⁴ K, and the heat conducted out through the two walls together equals the shear work μU²/h = 1 W/m² (`core.navier_stokes.internal_energy_residual`, `ch04.couette_heating`) |
| R08 | Eq. (4.61): Gibbs relation along a particle $\frac{De}{Dt}=T\frac{Ds}{Dt}-p\frac{D(1/\rho)}{Dt}$ [#97] | 4.8 | B | RECAP | C10 | ch01 C35 (Gibbs (1.18), D10): reminded in one line with D/Dt |
| N80 | Eq. (4.62): entropy equation $\frac{Ds}{Dt}=-\frac1\rho\frac{\partial}{\partial x_i}\big(\frac{q_i}{T}\big)-\frac{q_i}{\rho T^2}\frac{\partial T}{\partial x_i}+\frac\varepsilon T$ [#98] | 4.8 | B | NOTE | C10 | stated (a-D35 given in §4c: combine (4.57) with (4.61), quotient rule); `core.navier_stokes.entropy_terms` |
| N81 | Eq. (4.63): entropy production $\frac{k}{\rho T^2}\lvert\nabla T\rvert^2+\frac\varepsilon T\ge0$ ⇒ μ, μ_v, k ≥ 0; inviscid, non-conducting ⇒ isentropic particles [#99] | 4.8 | B | NOTE | C10 | stated with the second-law demo (a negative μ gives negative production, `core.navier_stokes.entropy_production`); ⚠️ the book's “μ, κ, k > 0” means μ, μ_v, k ≥ 0 (analysis §9.4); → Ch. 15 |
| N82 | Eq. (4.64): angular momentum of a body, $d\mathbf H/dt=\mathbf M$ [#100] | 4.9 | C | NOTE | C04 | named (pre-book mechanics) as the starting line of the angular-momentum paragraph in C04 |
| N83 | Eq. (4.65): angular-momentum principle for a stationary CV, $\frac{d}{dt}\int_{V_o}(\mathbf r\times\rho\mathbf u)dV+\int_{A_o}(\mathbf r\times\rho\mathbf u)(\mathbf u\cdot\mathbf n)dA=\int_{V_o}(\mathbf r\times\rho\mathbf g)dV+\int_{A_o}(\mathbf r\times\mathbf f)dA$ (Fig. 4.11) [#101] | 4.9 | B | NOTE | C04 | stated (a-D36 given in §4c: the RTT with F = r × ρu); `core.conservation.angular_momentum_flux`; absolute angular momentum → Ch. 13 |
| N84 | Ex. 4.6: lawn sprinkler held still, torque $M=2a\rho AU^2\cos\alpha$ (Fig. 4.12) [#102] | 4.9 | B | NOTE | C04 | stated with a-D37 given in §4c and one number: a = 0.2 m, water, A = 1 cm², U = 5 m/s, α = 30° → M ≈ 0.87 N m (`ch04.sprinkler_torque`); our extension: free-spinning rate U cos α/a |
| N85 | Bernoulli equations are not new laws: consequences of (4.38) and (4.60) under stated conditions [#103] | 4.9 | C | NOTE | C11 | named as C11's opening sentence and the header of the decision table (N102) |
| N86 | Eq. (4.66): Euler with a gravity potential, $\frac{\partial u_j}{\partial t}+u_i\frac{\partial u_j}{\partial x_i}=-\frac1\rho\frac{\partial p}{\partial x_j}-\frac{\partial\Phi}{\partial x_j}$ [#104] | 4.9 | B | NOTE | C11 | stated as D24's starting line ((4.41) with (4.18)) |
| N87 | Eq. (4.67): barotropic ρ = ρ(p): $\frac1\rho\frac{\partial p}{\partial x_j}=\frac{\partial}{\partial x_j}\int_{p_o}^p\frac{dp'}{\rho(p')}$ (constant-density, isothermal and isentropic flows) [#105] | 4.9 | B | NOTE | C11 | stated as step 2 of D24 (FTC with a variable limit, P84 reminder); closed forms: (p − p_o)/ρ, RT ln(p/p_o), h − h_o (`core.bernoulli.pressure_function`); ⚠️ the prime here is a dummy variable |
| N88 | Eq. (4.68): Lamb identity $u_i\frac{\partial u_j}{\partial x_i}=-(\mathbf u\times\boldsymbol\omega)_j+\frac{\partial}{\partial x_j}\big(\tfrac12u_i^2\big)$ [#106] | 4.9 | B | NOTE | C11 | stated with its proof as steps 3–6 of D24 (book defers it to Exercise 4.50; ε–δ recap ch02 D09); `core.navier_stokes.lamb_vector`; the vorticity equation of Ch. 5 starts here |
| N89 | Eq. (4.69): $\frac{\partial u_j}{\partial t}+\frac{\partial}{\partial x_j}\big[\tfrac12u_i^2+\int_{p_o}^p\frac{dp'}{\rho(p')}+gz\big]=(\mathbf u\times\boldsymbol\omega)_j$ — the Bernoulli function B [#107] | 4.9 | B | NOTE | C11 | stated as D24's result; `core.bernoulli.bernoulli_function` |
| N90 | Eq. (4.70): steady flow $\nabla B=\mathbf u\times\boldsymbol\omega$; B constant on Lamb surfaces containing streamlines and vortex lines (Fig. 4.13) [#108] | 4.9 | B | NOTE | C11 | stated as step 1 of D25 with a Lamb-surface picture (swirling Rankine vortex); `core.bernoulli.lamb_surface_check` |
| C11 | Eq. (4.71): Bernoulli 1 — steady, inviscid, barotropic, conservative body force: $\tfrac12u_i^2+\int_{p_o}^p\frac{dp'}{\rho(p')}+gz$ constant along streamlines and vortex lines [#109] | 4.9 | A | CORE | – | load-bearing: explains when Bernoulli holds and when it fails (across streamlines of a rotational flow); the Lamb identity behind it starts the vorticity equation (Ch. 5); pressure in potential flow (Ch. 6), Crocco's theorem (Ch. 14–15) |
| N91 | Eq. (4.72): irrotational flow — B constant everywhere [#110] | 4.9 | B | NOTE | C11 | stated as D25's last step (ω = 0 ⇒ ∇B = 0); number: B identical at 1000 random points of the cylinder flow vs varying across circles inside a Rankine core (`core.bernoulli.bernoulli_along_line`); E6 presets |
| N92 | Irrotational flow stays irrotational (barotropic, inviscid, non-rotating frame); flow outside a boundary layer (Fig. 4.14) [#111] | 4.9 | C | NOTE | C11 | named with pointers: proof by Kelvin's theorem in Ch. 5, the boundary layer in Ch. 9 |
| R09 | Eq. (4.73): velocity potential $\mathbf u\equiv\nabla\phi$ [#112] | 4.9 | B | RECAP | C12 | ch03 N29 ((3.17), irrotational flow, simply-connected caveat): reminded with `ch03` potential helpers; step 1 of D26 |
| N93 | Eq. (4.74): $\nabla\big[\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int_{p_o}^p\frac{dp'}{\rho}+gz\big]=0$, so the bracket is B(t); φ redefined to absorb B(t) [#113] | 4.9 | B | NOTE | C12 | stated as steps 3–6 of D26 with the sign fix: φ_new = φ − ∫B dt′ (the book prints +; a test shows the printed sign doubles B) |
| C12 | Eq. (4.75): Bernoulli 2 — unsteady, inviscid, irrotational, barotropic: $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int_{p_o}^p\frac{dp'}{\rho(p')}+gz$ = constant [#114] | 4.9 | A | CORE | – | load-bearing: the dynamic free-surface condition of every wave in Ch. 7, added mass of accelerating bodies (Ch. 6), acoustics and bubbles (Ch. 15); the only Bernoulli form valid in unsteady flow everywhere |
| N94 | Eq. (4.76): steady energy equation with σ = q = 0 [#115] | 4.9 | B | NOTE | C11 | stated as the first line of the (4.78) paragraph (a-D41 given in §4c) |
| N95 | Eq. (4.77): $\rho u_i\frac{\partial}{\partial x_i}\big(e+\frac p\rho+\tfrac12u_j^2+gz\big)=0$ [#116] | 4.9 | B | NOTE | C11 | stated as the second line of the (4.78) paragraph |
| N96 | Eq. (4.78): Bernoulli 3 (energy form) — $h+\tfrac12\lvert\mathbf u\rvert^2+gz$ constant on streamlines (steady, inviscid, non-conducting); isentropic ⇒ ∫dp/ρ = h links it to (4.71) [#117] | 4.9 | B | NOTE | C11 | stated with one number: air at 300 K, U = 100 m/s → stagnation temperature T + U²/2C_p ≈ 305 K (`core.bernoulli.stagnation_temperature`); → Ch. 15; E6 preset “high-speed gas” |
| N97 | Eq. (4.79): Lamb form of (4.39b), $\rho\frac{\partial\mathbf u}{\partial t}+\rho\nabla(\tfrac12\lvert\mathbf u\rvert^2)-\rho\mathbf u\times\boldsymbol\omega=-\nabla p+\rho\mathbf g-\mu\nabla\times\boldsymbol\omega$ [#118] | 4.9 | B | NOTE | C12 | stated (combines (4.68) and (4.40); a-D42 given in §4c) |
| N98 | Eq. (4.80): irrotational, constant ρ: the viscous term vanishes, $\rho\frac{\partial\mathbf u}{\partial t}+\nabla\big(\tfrac12\rho\lvert\mathbf u\rvert^2+\rho gz+p\big)=0$ [#119] | 4.9 | B | NOTE | C12 | stated with the demo that −μ∇×ω of any potential flow is exactly 0 (`core.bernoulli.viscous_irrotational_residual`); viscous wave damping → Ch. 7 |
| N99 | Eq. (4.81): (4.80) dotted with ds along a streamline and integrated from 1 to 2 [#120] | 4.9 | B | NOTE | C12 | stated as the middle line of the (4.82) paragraph (directional derivative P75 reminder) |
| N100 | Eq. (4.82): Bernoulli 4 — $\int_1^2\frac{\partial\mathbf u}{\partial t}\cdot d\mathbf s+\big(\tfrac12\lvert\mathbf u\rvert^2+gz+\frac p\rho\big)_2=\big(\ldots\big)_1$ (constant ρ, μ; irrotational; along a streamline at one instant) [#121] | 4.9 | B | NOTE | C12 | stated with C12's worked number: water column L = 1 m accelerating at 1 m/s² needs p₁ − p₂ = ρL dU/dt = 1000 Pa (`core.bernoulli.unsteady_streamline_bernoulli`); U-tube and seiches → Ch. 7 |
| N101 | Eq. (4.83): Bernoulli 4 in potential form $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+gz+\frac p\rho$ = constant [#122] | 4.9 | B | NOTE | C12 | stated as the constant-ρ case of (4.75), now valid with constant μ too |
| N102 | Summary of the four Bernoulli equations and their hypotheses ((4.19), (4.71), (4.75), (4.78), (4.82)/(4.83)) [#123] | 4.9 | B | NOTE | C11 | stated as a decision table (steady? viscous? irrotational? barotropic? isentropic? → which form, constant along what); `core.bernoulli.BERNOULLI_FORMS`, `which_bernoulli`; E6's status line and quiz |
| N103 | Pitot tube (Fig. 4.15): $\lvert\mathbf u\rvert_1=\sqrt{2(p_2-p_1)/\rho}=\sqrt{2g(h_2-h_1)}$ [#124] | 4.9 | B | NOTE | C05 | stated (a-D43 given in §4c) with one number: Δp = 500 Pa in air (ρ = 1.2 kg/m³) → 28.9 m/s (`core.bernoulli.pitot_speed`); E6 preset; airspeed indicators → Ch. 14 |
| N104 | Stagnation (total) pressure $p+\tfrac12\rho\lvert\mathbf u\rvert^2$ and dynamic pressure $\tfrac12\rho u^2$ [#125] | 4.9 | B | NOTE | C05 | stated with the same number (28.9 m/s ⇒ dynamic pressure 500 Pa); `core.bernoulli.stagnation_pressure`, `dynamic_pressure` |
| N105 | Orifice (Figs. 4.16–4.17): $u=\sqrt{2gh}$, $\dot m=\rho A_c\sqrt{2gh}$; vena contracta, sharp vs rounded [#126] | 4.9 | B | NOTE | C05 | stated with one number: h = 1 m → 4.43 m/s (`core.bernoulli.torricelli_speed`, `ch04.orifice_mass_flow`); tank-draining ODE figure; contraction coefficient 0.611 (cited; book value in the JSON) |
| R10 | Hydrostatic reference state $0=-\nabla p_s+\rho_s\mathbf g$ [#127] | 4.9 | B | RECAP | C13 | ch01 C20 (hydrostatic law (1.8), D05): reminded with `core.statics.integrate_hydrostatic`; step 1 of D27 |
| N106 | Eq. (4.84): $\rho\frac{D\mathbf u}{Dt}=-\nabla p'+\rho'\mathbf g+\mu\nabla^2\mathbf u$ with $p'=p-p_s$, $\rho'=\rho-\rho_s$ [#128] | 4.9 | B | NOTE | C13 | stated as steps 1–3 of D27; ⚠️ primes now mean perturbations (four meanings of ′ in this chapter, notation callout); `core.navier_stokes.perturbation_fields` |
| N107 | Eq. (4.85): constant density — gravity disappears, $\rho\frac{D\mathbf u}{Dt}=-\nabla p'+\mu\nabla^2\mathbf u$; it returns with a free surface, interface or density variation [#129] | 4.9 | B | NOTE | C13 | stated with the demo: the same Poiseuille velocity solves (4.39b) with gravity and (4.85) without (`core.navier_stokes.ns_incompressible_terms`) |
| N108 | Boussinesq conditions: low Mach, no sound, $L\ll c^2/g$, small δT so αδT ≪ 1; then $\frac{(1/\rho)(D\rho/Dt)}{\nabla\cdot\mathbf u}\sim\alpha\,\delta T\ll1$ ⇒ ∇·u = 0 [#130] | 4.9 | B | NOTE | C13 | stated (a-D45 given in §4c) with a validity table: water α ≈ 2×10⁻⁴ K⁻¹, δT = 10 K → 2×10⁻³; air α = 1/T, δT = 10 K at 300 K → 0.033; c²/g for air at 288 K ≈ 11.8 km (`ch04.boussinesq_validity`); E8 stage |
| C13 | Eq. (4.86): Boussinesq momentum $\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u$ with the heat equation (4.89) $\frac{DT}{Dt}=\kappa\nabla^2T$ [#131] | 4.9 | A | CORE | – | load-bearing: the equation set of the ocean and atmosphere in Ch. 13 and of internal waves (Ch. 7 §7.8), convection and stratified instability (Ch. 11), buoyant turbulence (Ch. 12) — Shammunul's daily workhorse |
| N109 | Eq. (4.87): internal energy in vector form $\rho\frac{De}{Dt}=-p\nabla\cdot\mathbf u+\rho\varepsilon-\nabla\cdot\mathbf q$ [#132] | 4.9 | B | NOTE | C13 | stated as D28's starting line ((4.60) restated) |
| N110 | Pressure work is not negligible although ∇·u ≈ 0: $-p\nabla\cdot\mathbf u\cong-p\alpha\frac{DT}{Dt}=-\rho(C_p-C_v)\frac{DT}{Dt}$ (perfect gas), turning C_v into C_p [#133] | 4.9 | B | NOTE | C13 | stated as steps 2–5 of D28 (sympy: e = C_vT, p = ρRT, α = 1/T ⇒ ρC_v DT/Dt + p∇·u = ρC_p DT/Dt); ⚠️ liquids: the term is small and C_p ≈ C_v anyway (the ocean uses the same (4.89)) |
| N111 | Eq. (4.88): $\rho C_p\frac{DT}{Dt}=\rho\varepsilon-\nabla\cdot\mathbf q$ [#134] | 4.9 | B | NOTE | C13 | stated as step 6 of D28; `core.navier_stokes.heat_equation_terms` |
| N112 | Viscous heating negligible under Boussinesq: $\frac{\rho\varepsilon}{\rho C_p(DT/Dt)}\sim\frac{\nu U}{C_p\delta T\,L}$ [#135] | 4.9 | B | NOTE | C13 | stated with one number (water, U = 0.1 m/s, L = 1 m, δT = 1 K: ~2×10⁻¹¹ — tiny); inside `ch04.boussinesq_validity` |
| N113 | Eq. (4.89): Boussinesq heat equation $\frac{DT}{Dt}=\kappa\nabla^2T$, $\kappa\equiv k/\rho C_p$ [#136] | 4.9 | B | NOTE | C13 | stated as D28's result (κ recap ch01); exact test solution: a Gaussian blob advected at U and spreading with σ² = σ₀² + 2κt (`core.navier_stokes.temperature_equation_residual`, `ch04.gaussian_blob_advection_diffusion`) — C13's animation |
| N114 | Boussinesq set: (4.10), (4.86) with g = −g e_z, (4.89), linear equation of state $\rho=\rho_0[1-\alpha(T-T_0)]$ [#137] | 4.9 | B | NOTE | C13 | stated as a boxed summary (`ch04.boussinesq_density`; link to `ch01.seawater_density_linear`); buoyancy b = −gρ′/ρ₀ and N² = ∂b/∂z (recap ch01 C51); → Ch. 7, 11, 13 |
| N115 | Boundary-condition specification: velocity on all bounding surfaces; external flows need the state on a distant closed surface [#138] | 4.10 | C | NOTE | C14 | named at the head of C14; realised by the solvers of Ch. 8–10 |
| N116 | Interface conditions from a pillbox (Fig. 4.18): $\rho_1\mathbf u_1\cdot\mathbf n=\rho_2\mathbf u_2\cdot\mathbf n$, traction $n_i\tau^{(1)}_{ij}=n_i\tau^{(2)}_{ij}$ continuous, heat flux $k_1\partial T_1/\partial n=k_2\partial T_2/\partial n$ [#139] | 4.10 | B | NOTE | C14 | stated (a-D48 given in §4c: volume and side terms vanish as l → 0) with two worked cases: composite wall (flux continuous, T kinked in the ratio k₂/k₁) and two-fluid Couette (shear stress continuous, u kinked) — `ch04.two_layer_conduction`, `two_fluid_couette`, `core.interfaces.pillbox_limit` |
| N117 | No-slip $\mathbf u_1\cdot\mathbf t=0$ and no temperature jump $T_1=T_2$ at a solid wall — not from conservation laws; exceptions (superfluid helium, super-hydrophobic textures, rarefied gases Kn ~ 1) [#140] | 4.10 | B | NOTE | C14 | stated with a Couette profile with no-slip vs a Navier slip length (our extension, labelled) and `ch01.knudsen_number` to show when slip matters; why boundary layers exist → Ch. 9 |
| N118 | Eq. (4.90): a point riding on a moving surface η(x, t) = 0 sees $d\eta/dt=\partial\eta/\partial t+(\mathbf u_s\cdot\nabla)\eta=0$ [#141] | 4.10 | B | NOTE | C14 | stated as steps 1–3 of D29; `core.interfaces.surface_normal_speed` = −(∂η/∂t)/∣∇η∣ |
| C14 | Eq. (4.91): kinematic boundary condition $\partial\eta/\partial t+(\mathbf u\cdot\nabla)\eta\equiv D\eta/Dt=0$ on η = 0 — with the wall and interface conditions it completes the problem [#142] | 4.10 | A | CORE | – | load-bearing: the free-surface condition of every wave problem (Ch. 7), layer interfaces in shallow-water and two-layer models (Ch. 13), moving shocks (Ch. 15); paired with no-slip and the pillbox jump conditions it tells every later solver what happens at an edge |
| N119 | Eq. (4.92): relative normal velocity $(u_{rel})_n=\big(\mathbf u\cdot\nabla\eta+\partial\eta/\partial t\big)/\lvert\nabla\eta\rvert=(1/\lvert\nabla\eta\rvert)D\eta/Dt$ [#143] | 4.10 | B | NOTE | C14 | stated as D29's last step; `core.interfaces.relative_normal_velocity` |
| N120 | Eq. (4.93): mass flux per area through the moving surface $(\rho/\lvert\nabla\eta\rvert)D\eta/Dt$ on η = 0 [#144] | 4.10 | B | NOTE | C14 | stated; `core.interfaces.interface_mass_flux`; shock mass flux → Ch. 15 |
| R11 | Surface tension revisited: unbalanced molecular attraction; the surface contracts as far as allowed [#145] | 4.10 | C | RECAP | C14 | ch01 C15 (surface tension, σ(T) of water): one reminder sentence |
| N121 | Eq. (4.94): Helmholtz free energy per unit mass $f=e-Ts$ [#146] | 4.10 | B | NOTE | C14 | stated (a Legendre device like ch01's Gibbs g = h − Ts, P50 reminder); `core.thermo.helmholtz_free_energy` |
| N122 | Eq. (4.95): $df=de-T\,ds-s\,dT$; isothermal reversible: dF = −p dv (work done on the system) [#147] | 4.10 | B | NOTE | C14 | stated with a one-line sympy check for a perfect gas at constant T |
| N123 | Eq. (4.96): $F=\rho_1V_1f_1+\rho_2V_2f_2+A\sigma$ — σ is free energy per area; σ > 0 immiscible (area minimised) [#148] | 4.10 | B | NOTE | C14 | stated (a-D52 given in §4c) with one figure: spheroid area at fixed volume is least at aspect 1 (`ch04.spheroid_area`) |
| N124 | Marangoni flows driven by surface-tension gradients [#149] | 4.10 | C | NOTE | C14 | named in one sentence; not treated in the book (pointer to Ch. 16 and interfacial flows) |
| N125 | Eq. (4.97): net pressure force on a small curved cap $\eta=z-x^2/2R_1-y^2/2R_2$, $(F_p)_z=-\pi\Delta p\sqrt{2R_1\zeta}\sqrt{2R_2\zeta}$ (Fig. 4.19) [#150] | 4.10 | B | NOTE | C14 | stated (a-D50 given in §4c: n_z dA = dx dy, so the force is Δp × ellipse area) with `core.interfaces.cap_pressure_force` (dblquad vs closed form) |
| N126 | Eq. (4.98): surface-tension force $\sigma\oint_C\mathbf t\times\mathbf n\,ds$, small-cap limit $(F_{st})_z=\pi\sigma\sqrt{2R_1\zeta}\sqrt{2R_2\zeta}\big(\frac1{R_1}+\frac1{R_2}\big)$ [#151] | 4.10 | B | NOTE | C14 | stated with the exact-ellipse line integral converging to the closed form, ratio − 1 = O(ζ) (`core.interfaces.cap_surface_tension_force`); ⚠️ the book's curve C has a sign slip (analysis §9.5) |
| N127 | Laplace's jump (1.5) derived from the cap force balance: $\Delta p=\sigma(1/R_1+1/R_2)$ [#152] | 4.10 | B | NOTE | C14 | stated as the balance F_p + F_st = 0 (ch01 N07 stated it; now it is derived); number: water drop R = 1 mm, σ = 0.0728 N/m → 146 Pa (`core.interfaces.laplace_jump_from_balance` = `ch01.laplace_pressure_jump`); capillary waves → Ch. 7 |
| N128 | Bubble below a free surface; capillary length $\ell_c=(\sigma/\rho g)^{1/2}$ [#153] | 4.10 | B | NOTE | C14 | stated with one number: water at 20 °C → 2.72 mm (`core.interfaces.capillary_length`; book value in the JSON); gravity–capillary crossover → Ch. 7 |
| N129 | Ex. 4.7: meniscus at a vertical wall, $\big(\frac{\rho g}{2\sigma}\big)\zeta^2+(1+\zeta'^2)^{-1/2}=1$, $h^2=\frac{2\sigma}{\rho g}(1-\sin\theta)$ (Fig. 4.20) [#154] | 4.10 | B | NOTE | C14 | stated (a-D51 ★★★ given in §4c) with the wall height and profiles for several contact angles in units of δ (closed form vs `solve_ivp`, `ch04.meniscus_profile_x`, `meniscus_profile_ode`); θ = 0 → h = √2 δ ≈ 3.85 mm for water |
| N130 | Two routes to dimensionless groups (from the equations vs dimensional analysis); dynamic similarity = equal groups + geometric similarity [#155] | 4.11 | B | NOTE | C15 | stated at the head of C15 with the “same equation, same solution” picture (two dimensional Stokes-layer solutions collapsing on one curve) |
| R12 | Eq. (4.99): sphere drag $\frac{F_D}{\rho U^2d^2}=\Psi\big(\frac{\rho Ud}{\mu}\big)$ or $\frac{F_D\rho}{\mu^2}=\Phi\big(\frac{\mu}{\rho Ud}\big)$ [#156] | 4.11 | B | RECAP | C15 | ch01 C69 (Π theorem, D28, `SPHERE_DRAG` preset): reminded with `core.dimensional.pi_groups` and two repeating sets (one group is the other times Re²; a-D53 recap in §4c) |
| N131 | Fig. 4.21: sphere $C_D=D/\tfrac12\rho U^2A$ vs $\mathrm{Re}=\rho Ud/\mu$ — 24/Re at low Re, nearly constant above ~10³, drag crisis [#157] | 4.11 | B | NOTE | C15 | stated with our figure: Morrison correlation (V5) + Stokes line on log–log axes and synthetic “experiments” at many (d, U, ρ, μ) collapsing onto one curve (`core.similarity.sphere_drag_coefficient`); number: d = 1 cm, U = 1 m/s in water → Re = 10⁴; drag crisis → Ch. 9; E9 view |
| N132 | Flow parameters l, U, Ω, ρ, μ of a general unsteady flow (pulsating tube flow, swimmer, turbomachine) [#158] | 4.11 | C | NOTE | C15 | named as the setting of (4.100) |
| N133 | Eq. (4.100): scaled variables $x_i^*=x_i/l$, $t^*=\Omega t$, $u_j^*=u_j/U$, $p^*=(p-p_\infty)/\rho U^2$, $g_j^*=g_j/g$ [#159] | 4.11 | B | NOTE | C15 | stated as step 1 of D30; `core.similarity.Scales`; ⚠️ time and pressure scales change inside the chapter (Ωt vs Ut/l; ρU², μU/l, ρgl) |
| C15 | Eq. (4.101): dimensionless Navier–Stokes $\big[\frac{\Omega l}{U}\big]\frac{\partial\mathbf u^*}{\partial t^*}+(\mathbf u^*\cdot\nabla^*)\mathbf u^*=-\nabla^*p^*+\big[\frac{gl}{U^2}\big]\mathbf g^*+\big[\frac{\mu}{\rho Ul}\big]\nabla^{*2}\mathbf u^*$ and dynamic similarity [#160] | 4.11 | A | CORE | – | load-bearing: where Re, Fr, St (and later Ro, Ek, Ri, Ma, Pr) come from — every later chapter classifies its flows by these groups (Stokes vs boundary-layer flow in Ch. 8–9, Fr in Ch. 7, Ri in Ch. 11–13, Ro in Ch. 13, M in Ch. 15), and model testing rests on it |
| N134 | Eq. (4.102): Strouhal number $\mathrm{St}\propto\frac{\partial u/\partial t}{u\,\partial u/\partial x}\propto\frac{\Omega l}U$ [#161] | 4.11 | B | NOTE | C15 | stated as a bracket of D30; `core.similarity.strouhal_number`; vortex shedding St ≈ 0.2 → Ch. 9 |
| N135 | Eq. (4.103): Reynolds number $\mathrm{Re}\propto\frac{\rho u\,\partial u/\partial x}{\mu\,\partial^2u/\partial x^2}\propto\frac{\rho Ul}{\mu}$ [#162] | 4.11 | B | NOTE | C15 | stated with a table of real values (bacterium ~10⁻⁵, sphere above 10⁴, ocean liner ~10⁹ — our computed numbers); `core.similarity.reynolds_number` |
| N136 | Eq. (4.104): Froude number $\mathrm{Fr}\propto\big[\frac{\rho U^2/l}{\rho g}\big]^{1/2}=\frac U{\sqrt{gl}}$; gravity matters only with a free surface or density differences [#163] | 4.11 | B | NOTE | C15 | stated (square root of a force ratio); `core.similarity.froude_number`; hydraulic jumps and ship waves → Ch. 7 |
| N137 | Eq. (4.105): internal Froude number $\mathrm{Fr}'=U/\sqrt{g'l}$, $g'=g(\rho_2-\rho_1)/\rho_1$, or $U/Nl$; Richardson $\mathrm{Ri}=1/\mathrm{Fr}'^2$, gradient Richardson $N^2/(dU/dz)^2$ [#164] | 4.11 | B | NOTE | C15 | stated with one number (ocean thermocline: Δρ/ρ = 10⁻³, l = 100 m, U = 0.1 m/s → g′ ≈ 9.8×10⁻³ m/s², Ri ≈ 98); N² recap ch01 C51 (Kundu Γ ≡ dT/dz convention); `core.similarity.richardson_number`, `gradient_richardson_number`; Ri > ¼ → Ch. 11; E9 mode |
| N138 | Under dynamic similarity every dimensionless quantity matches; local vs overall relations [#165] | 4.11 | C | NOTE | C15 | named as the sentence that introduces (4.106) |
| N139 | Eq. (4.106): pressure coefficient $C_p=\frac{p-p_\infty}{\tfrac12\rho U^2}=\Psi\big(\mathrm{St},\mathrm{Fr},\mathrm{Re};\frac{\mathbf x}l,\Omega t\big)$ [#166] | 4.11 | B | NOTE | C15 | stated with the ideal cylinder: C_p = 1 − 4 sin²θ identical for five (U, a) pairs (`core.similarity.pressure_coefficient`) |
| N140 | Steady boundary conditions: time scale l/U, $t^*=Ut/l$; spontaneous unsteadiness scales with l/U [#167] | 4.11 | C | NOTE | C15 | named (a `Scales` option); vortex shedding → Ch. 9 |
| N141 | Purely oscillatory body: U = lΩ ⇒ St = 1, Re = Ωl²/ν, Fr = Ω(l/g)^{1/2} [#168] | 4.11 | B | NOTE | C15 | stated (a-D56 given in §4c); `core.similarity.Scales.from_oscillation`; Stokes layer and Womersley → Ch. 8, 16 |
| N142 | Eq. (4.107): drag coefficient $C_D\equiv\frac{F_D}{\tfrac12\rho U^2A}$ [#169] | 4.11 | B | NOTE | C15 | stated; `core.similarity.drag_coefficient` |
| N143 | Eq. (4.108): lift coefficient $C_L\equiv\frac{F_L}{\tfrac12\rho U^2A}$; reference areas (frontal for blunt bodies, planform for plates and airfoils) [#170] | 4.11 | B | NOTE | C15 | stated; `core.similarity.lift_coefficient`, `reference_area`; → Ch. 14 |
| N144 | Ship: C_D = C_D(Fr, Re); far from a free surface and incompressible, C_D = C_D(Re) [#171] | 4.11 | C | NOTE | C15 | named as the bridge to Ex. 4.8 |
| N145 | Eq. (4.109): compressible scalings $t^*=Ut/l$, $p^*=(p-p_\infty)/\rho_oU^2$, $\rho^*=\rho/\rho_o$ [#172] | 4.11 | B | NOTE | C15 | stated; `core.similarity.Scales(time_scale='advective')` |
| N146 | Eq. (4.110): $\nabla^*\cdot\mathbf u^*=-\big[\frac{U^2}{c^2}\big]\frac1{\rho^*}\frac{Dp^*}{Dt^*}$ — departure from ∇·u = 0 scales with M² [#173] | 4.11 | B | NOTE | C15 | stated (a-D57 given in §4c) with the §4.11 re-display of (4.9); number: M = 0.3 → M² = 0.09; `core.similarity.compressibility_parameter` |
| N147 | Eq. (4.111): Mach number $M=U/c$ (square root of inertia over compressibility force); incompressible when M < 0.3 [#174] | 4.11 | B | NOTE | C15 | stated; `core.similarity.mach_number` (c from `core.thermo.perfect_gas_sound_speed`, recap ch01 C36); → Ch. 14, 15 |
| N148 | Eq. (4.112): enthalpy form $\rho\frac{Dh}{Dt}=\frac{Dp}{Dt}+\rho\varepsilon+\frac{\partial}{\partial x_i}\big(k\frac{\partial T}{\partial x_i}\big)$ [#175] | 4.11 | B | NOTE | C15 | stated (a-D58 given in §4c; sympy one-liner `core.navier_stokes.energy_forms_sym` proves (4.60) ⇔ (4.112)) |
| N149 | Eq. (4.113): thermal scalings $\varepsilon^*=\rho_ol^2\varepsilon/\mu_oU^2$, $\mu^*=\mu/\mu_o$, $k^*=k/k_o$, $T^*=(T-T_o)/(T_w-T_o)$ [#176] | 4.11 | B | NOTE | C15 | stated; `core.similarity.Scales(T_o, T_w, …)` |
| N150 | Eq. (4.114): dimensionless energy equation with coefficients Ec, Ec/Re and 1/(Pr Re) [#177] | 4.11 | B | NOTE | C15 | stated (a-D59 given in §4c); `core.similarity.nondimensional_energy_coefficients`; ⚠️ the book cites (4.106), (4.107) where (4.109), (4.112) are meant |
| N151 | Eq. (4.115): Eckert number $\mathrm{Ec}=U^2/C_p(T_w-T_o)$; low Ec ⇒ (4.112) → (4.89) [#178] | 4.11 | B | NOTE | C15 | stated; `core.similarity.eckert_number`; the Boussinesq heat equation's hidden assumption (C13) |
| N152 | Eq. (4.116): Prandtl number $\mathrm{Pr}=\nu/\kappa=\mu_oC_p/k_o$ (air ≈ 0.7, water ≈ 7) [#179] | 4.11 | B | NOTE | C15 | stated with our computed values from `ch01.FLUIDS` (`ch04.prandtl_of`) and the Eucken line 4γ/(9γ − 5); thermal boundary layers → Ch. 9, Rayleigh–Bénard → Ch. 11 |
| N153 | Eq. (4.117): Weber number $\mathrm{We}=\rho U^2l/\sigma$ (ratio of forces) [#180] | 4.11 | B | NOTE | C15 | stated; `core.similarity.weber_number`; → Ch. 7, 16 |
| N154 | Eq. (4.118): Bond number $\mathrm{Bo}=\rho l^2g/\sigma=(l/\ell_c)^2$ [#181] | 4.11 | B | NOTE | C15 | stated (link to the capillary length of C14); `core.similarity.bond_number` |
| N155 | Eq. (4.119): capillary number $\mathrm{Ca}=\mu U/\sigma$ (= We/Re) [#182] | 4.11 | B | NOTE | C15 | stated; `core.similarity.capillary_number`; → Ch. 16 |
| N156 | Ex. 4.8: ship model — Froude matching $U_m=U_p\sqrt{l_m/l_p}$, friction subtracted and re-added, wave drag scaled by $(\rho_p/\rho_m)(l_p/l_m)^2(U_p/U_m)^2$ [#183] | 4.11 | B | NOTE | C15 | stated (a-D60 given in §4c) with a bar chart model total → friction + wave → scaled wave + prototype friction (`ch04.ship_drag_extrapolation`, book numbers in the JSON); why Re cannot also match (ratio λ^{3/2}); E9 mode |
| S01 | Exercises 4.1–≈4.63 [#184] | Ex. | C | SKIP | – | pointer line: practise on the book's exercises; the derivations the text defers to Exercises 4.7, 4.8, 4.30, 4.38, 4.42, 4.43, 4.45–4.47, 4.50 are written out in our words (§4b); exercise text is not reproduced |
| S02 | Literature cited and supplemental reading [#185] | Lit. | C | SKIP | – | pointer line: bibliography (Aris for isotropic tensors, Batchelor, Lamb, Spiegel & Veronis for Boussinesq) |

### 2a. What each A block contains (for the lesson-designer)
- **C01** picture: a sealed balloon of air drifting and swelling vs a fixed box it passes through · question: "the
  balloon keeps its mass — what does the box see?" · D01 · number: expanding flow $u=ax/(1+at)$,
  $\rho=\rho_0/(1+at)$, a = 1 s⁻¹, ρ₀ = 1 kg/m³: material interval [1, 2] m keeps 1 kg/m²; the fixed box [1, 2] m at
  t = 0 has dM/dt = −(ρu)₂ + (ρu)₁ = −2 + 1 = −1 kg/(m² s) = ∫∂ρ/∂t dx · code: `ch04.expanding_flow`, `material_mass`,
  `core.conservation.mass_budget` (fixed, moving, material CV; residual → 0) + from-scratch midpoint sums · figure:
  one flow, three CVs (fixed, moving, material) with budget bars; animation of the material interval stretching ·
  explainer E1 · B/C: N01, R01, N03–N06.
- **C02** picture: a crowd leaving a square — the density drops where more people leave than arrive · question: "what
  does 'mass is conserved' say at a single point?" · D02, D03 · number: the expanding flow at t = 0: ∂ρ/∂t = −1,
  ∂(ρu)/∂x = +1 kg/(m³ s), sum 0; (1/ρ)Dρ/Dt = −1 s⁻¹ = −∇·u · code: `core.navier_stokes.continuity_residual(_sym)`,
  `continuity_terms`, `density_material_rate`, `ch04.stratified_shear_flow` + from-scratch central stencils · figure:
  `slider_figure` over time of the two terms of (4.7) along x; stratified shear flow (ρ varies, Dρ/Dt = 0) · explainer
  (none of its own; E2 and E4 show ∇·u = 0) · B/C: R02, N07–N12.
- **C03** picture: a weather map's streamlines — tight where the wind is strong · question: "can one scalar field carry
  both the direction and the speed of a 2-D flow?" · D04 · number: stagnation flow ψ = kxy, k = 1 s⁻¹: u = kx,
  v = −ky; between ψ = 1 and ψ = 2 m²/s flows 1 m²/s per metre of depth; where those contours are 0.1 m apart the speed
  is ≈ 10 m/s · code: `core.streamfunction.velocity_from_streamfunction_2d(_sym)`, `flux_between_streamlines`,
  `velocity_from_streamfunction_axisym` + from-scratch stencils and a trapezoid flux · figure: ψ contours with equal
  Δψ over a speed heatmap (cylinder flow, `ch03.cylinder_streamfunction`); plotly 3-D stream surfaces (N15) · explainer
  E2 · B/C: N13–N18.
- **C04** picture: a flat plate held in a wind tunnel — you measure only the wake, never the plate · question: "can you
  weigh a force by watching what flows out of a box?" · D05 · number: Ex. 4.1 Gaussian wake U∞ = 10 m/s, Δ = 2 m/s,
  b = 0.1 m, ρ = 1.2 kg/m³ → side leakage 0.354 m²/s, F_D/l ≈ 3.65 N/m; bore h_in = 1 m, h_out = 1.1 m → 3.37 m/s; the
  sprinkler 0.87 N m · code: `core.conservation.momentum_budget`, `ch04.wake_drag_per_span`, `gaussian_wake`,
  `bore_speed`, `rocket_trajectory`, `sprinkler_torque` + from-scratch face sums with `np.trapezoid` · figure: the CV
  with the four faces' mass and momentum fluxes as bars (our Fig. 4.2), F_D/l vs box height H converging; animation of
  the bore with the CV riding along · explainer E1 · B/C: N19–N29, N31, N32, R03, R04, N82–N84.
- **C05** picture: an airliner's pitot probe; a water tank with a hole · question: "why does fast flow have low
  pressure?" · D06 · number: pitot Δp = 500 Pa in air → 28.9 m/s; Torricelli h = 1 m → 4.43 m/s · code:
  `core.bernoulli.bernoulli_head`, `bernoulli_solve`, `pitot_speed`, `torricelli_speed`, `ch04.orifice_mass_flow`,
  `stream_tube_element_balance_sym` · figure: p along streamlines of the cylinder flow from (4.19) (C_p = 1 − 4 sin²θ on
  the surface); tank-draining curve with C_c = 0.611 vs 1 · explainer E6 · B/C: N30, N103–N105.
- **C06** picture: a small cube of fluid pushed by its neighbours on six faces and pulled by gravity · question: "what
  is F = ma for a fluid particle, before we know anything about the fluid?" · D07 · number: rigid rotation Ω = 1 rad/s
  with p = ρΩ²r²/2, water: at r = 0.1 m the stress divergence −∂p/∂r = −100 N/m³ equals ρ times the centripetal
  acceleration −Ω²r · code: `core.navier_stokes.cauchy_terms`, `momentum_conservative_residual`,
  `conservative_to_advective_sym`, `ch04.closure_count` + from-scratch first-index divergence (a non-symmetric τ shows
  why the index matters) · figure: term bars (inertia = body + stress divergence) at points of the rotating fluid; the
  ledger table · explainer E4 (Cauchy terms) · B/C: N02, N33–N37.
- **C07** picture: honey sheared between two plates, then squeezed · question: "how does the fluid decide its stress
  from how it deforms?" · D08, D09 ★★★, D10 · number: shear γ = 10 s⁻¹, water μ = 1.0×10⁻³ Pa s → τ₁₂ = 0.01 Pa;
  extension u = (εx, −εy, 0), ε = 1 s⁻¹ → τ₁₁ + p = 2 mPa, τ₂₂ + p = −2 mPa · code: `core.constitutive.newtonian_stress`,
  `viscous_stress`, `isotropic_fourth_order`, `linear_stress`, `mean_pressure`, `ch04.cube_spin_acceleration` +
  from-scratch loops over i, j (τ_ij = −pδ_ij + 2μS_ij + λS_mmδ_ij) · figure: normal and shear stress on a plane vs its
  angle (Mohr-type curve) for shear vs extension; cube spin rate vs h on log–log (slope −2) · explainer E3 · B/C:
  N38–N51, R05, R06.
- **C08** picture: honey in a pipe and air round a wing — same equation, different terms matter · question: "which of
  the five terms of Navier–Stokes balance where?" · D11, D12, D13 · number: plane Poiseuille, water, h = 1 mm,
  G = −dp/dx = 100 Pa/m → u_max = Gh²/(8μ) = 12.5 mm/s; pressure force +100 N/m³ balanced by viscous −100 N/m³ at every
  point · code: `core.navier_stokes.ns_incompressible_terms`, `navier_stokes_residual`, `navier_stokes_sym`,
  `viscous_force_forms` + from-scratch Poiseuille residual · figure: term bars across the channel; Stokes' first problem
  profiles (`slider_figure` in t); solid-body rotation vs Lamb–Oseen viscous force · explainer E4 · B/C: N52–N55.
- **C09** picture: a ball rolled across a spinning merry-go-round; a hurricane turning · question: "why does the ball
  curve for the rider but not for someone standing outside?" · D14, D15 ★★★, D16, D17, D18 · number: Ω = 7.292×10⁻⁵
  rad/s at the pole, u = 10 m/s, t = 1 h → Coriolis acceleration 1.46×10⁻³ m/s², forward 36 km, deflection 9.45 km,
  angle 15.0°; centrifugal Ω²a = 0.0339 m/s² at the equator · code: `core.rotating.frame_acceleration_terms`,
  `apparent_body_forces`, `coriolis_acceleration`, `effective_gravity`, `ch04.coriolis_projectile` + from-scratch
  `np.cross` terms vs finite differences of a path seen in both frames · figure: animation of the projectile in the
  inertial (straight) and rotating (curved) frames side by side; quiver of Coriolis arrows on flow out of a high;
  effective-gravity arrows vs latitude · explainer E5 · B/C: N56–N65.
- **C10** picture: a spoon stirring coffee — the swirl dies and the coffee warms (by a tiny amount) · question: "where
  does kinetic energy go when viscosity stops a flow?" · D19–D23 · number: plane Couette U = 1 m/s, h = 1 mm, water:
  ρε = μ(U/h)² = 1000 W/m³; heat out through both walls = shear work μU²/h = 1 W/m²; peak temperature rise μU²/(8k) ≈
  2.1×10⁻⁴ K · code: `core.constitutive.dissipation_rate(form='both')`, `core.navier_stokes.kinetic_energy_budget`,
  `internal_energy_terms`, `internal_energy_residual`, `entropy_production`, `ch04.couette_heating` + from-scratch
  double sum σ_ijS_ij/ρ · figure: Couette/Poiseuille profiles of u, ε(y) and T(y) with the energy budget bars; ε vs
  random G (always ≥ 0) · explainer E7 · B/C: N66–N81, R07, R08.
- **C11** picture: a whirlpool — the water surface dips at the centre, so "Bernoulli" cannot hold across the circles ·
  question: "Bernoulli is constant along what, exactly?" · D24, D25 · number: Rankine vortex (Γ = 2π m²/s, σ = 1 m,
  water): inside the core B changes from circle to circle, outside it is the same everywhere; cylinder potential flow:
  B equal at 1000 random points · code: `core.navier_stokes.lamb_vector`, `lamb_identity_terms`,
  `core.bernoulli.bernoulli_function`, `bernoulli_along_line`, `pressure_function`, `stagnation_temperature`,
  `which_bernoulli` + from-scratch u × ω by `np.cross` · figure: B along vs across streamlines of the Rankine vortex;
  the decision table · explainer E6 · B/C: N85–N92, N94–N96, N102.
- **C12** picture: water sloshing in a U-tube; a balloon suddenly pushed through water · question: "what does the
  pressure do when the flow speeds up?" · D26 · number: a water column L = 1 m accelerated at 1 m/s² needs p₁ − p₂ =
  1000 Pa; a sphere accelerating from rest feels an extra force of half its displaced mass times dU/dt · code:
  `core.bernoulli.unsteady_bernoulli_pressure`, `unsteady_bernoulli_B`, `unsteady_streamline_bernoulli`,
  `viscous_irrotational_residual` · figure: animation of a U-tube column with the pressure difference bar tracking
  ρL dU/dt; surface pressure on an accelerating sphere (added-mass term) · explainer E6 · B/C: R09, N93, N97–N101.
- **C13** picture: warm water rising in a cold lake; the sea breeze · question: "density hardly changes — so why does it
  drive the flow?" · D27, D28 · number: water α ≈ 2×10⁻⁴ K⁻¹, δT = 10 K → αδT = 2×10⁻³ and reduced gravity
  g′ = 0.0196 m/s²; air at 300 K → αδT = 0.033; c²/g ≈ 11.8 km for air at 288 K · code:
  `core.navier_stokes.boussinesq_momentum_terms`, `buoyancy`, `perturbation_fields`, `heat_equation_terms`,
  `temperature_equation_residual`, `ch04.boussinesq_validity`, `boussinesq_density`,
  `gaussian_blob_advection_diffusion` + from-scratch heat-equation residual · figure: animation of a warm Gaussian blob
  advected and spreading (σ² = σ₀² + 2κt); validity table for ocean, lake, boundary layer, deep atmosphere · explainer
  E8 · B/C: R10, N106–N114.
- **C14** picture: a wave on the sea — the water on the surface stays on the surface · question: "what does 'the
  surface is made of fluid' say as an equation?" · D29 · number: a rigid wall moving at V = 1 m/s (η = x − Vt): the
  surface's normal speed −(∂η/∂t)/∣∇η∣ = 1 m/s equals u·n, residual 0; a linear wave with ka = 0.1: the full
  condition's residual ~ (ka)² ≈ 0.01, a quarter of that at ka = 0.05; water drop R = 1 mm → Δp = 146 Pa; capillary
  length 2.72 mm · code: `core.interfaces.kinematic_bc_residual`, `surface_normal_speed`, `relative_normal_velocity`,
  `pillbox_limit`, `cap_pressure_force`, `cap_surface_tension_force`, `laplace_jump_from_balance`, `capillary_length`,
  `ch04.two_layer_conduction`, `two_fluid_couette`, `meniscus_profile_x` + from-scratch Dη/Dt stencil · figure:
  animation of a surface wave with particles riding on it; composite wall and two-fluid Couette (kinked profiles,
  continuous fluxes); meniscus profiles · explainer backup B1 · B/C: N115–N129, R11.
- **C15** picture: a ship model in a towing tank; a toy airplane in a wind tunnel · question: "when does a small model
  behave like the real thing?" · D30 · number: sphere d = 1 cm, U = 1 m/s in water → Re = 10⁴ (C_D ≈ 0.4 on the
  plateau); Froude scaling at λ = 1/25: U_m = U_p/5 and Re_m/Re_p = λ^{3/2} = 1/125 — Re cannot also match · code:
  `core.similarity.nondimensional_ns_coefficients`, `Scales`, `reynolds_number`, `froude_number`, `richardson_number`,
  `sphere_drag_coefficient`, `ch04.ship_drag_extrapolation`, `stokes_first_problem` + from-scratch sympy substitution
  of (4.100) into (4.39b) · figure: two dimensional Stokes-layer solutions collapsing onto one dimensionless curve;
  sphere C_D(Re) with synthetic experiments collapsing; ship-drag bar chart · explainer E9 · B/C: N130–N156, R12.

## 3. Section coverage
| § | Title | A | B | C | RECAP | SKIP |
|---|---|---|---|---|---|---|
| 4.1 | Introduction | — (N01 inside C01; N02 inside C06) | N02 | N01 | — | — |
| 4.2 | Conservation of Mass | C01, C02 | N03, N04, N05, N06, N07, N08, N09, N10, N11, N12 | — | R01 (B), R02 (B) | — |
| 4.3 | Stream Functions | C03 | N13, N14, N15, N16, N17 | N18 | — | — |
| 4.4 | Conservation of Momentum | C04, C05, C06 | N19, N20, N21, N22, N23, N24, N25, N26, N27, N28, N29, N30, N31, N32, N33, N34, N35, N36, N37 | — | R03 (B), R04 (C) | — |
| 4.5 | Constitutive Equation for a Newtonian Fluid | C07 | N38, N40, N41, N42, N43, N44, N45, N46, N47, N48, N49, N50, N51 | N39 | R05 (B), R06 (C) | — |
| 4.6 | Navier-Stokes Momentum Equation | C08 | N52, N53, N54, N55 | — | — | — |
| 4.7 | Noninertial Frame of Reference | C09 | N57, N58, N59, N60, N61, N62, N64, N65 | N56, N63 | — | — |
| 4.8 | Conservation of Energy | C10 | N66, N67, N68, N69, N70, N71, N72, N73, N74, N75, N76, N77, N78, N79, N80, N81 | — | R07 (B), R08 (B) | — |
| 4.9 | Special Forms of the Equations | C11, C12, C13 | N83, N84, N86, N87, N88, N89, N90, N91, N93, N94, N95, N96, N97, N98, N99, N100, N101, N102, N103, N104, N105, N106, N107, N108, N109, N110, N111, N112, N113, N114 | N82, N85, N92 | R09 (B), R10 (B) | — |
| 4.10 | Boundary Conditions | C14 | N116, N117, N118, N119, N120, N121, N122, N123, N125, N126, N127, N128, N129 | N115, N124 | R11 (C) | — |
| 4.11 | Dimensionless Forms of the Equations and Dynamic Similarity | C15 | N130, N131, N133, N134, N135, N136, N137, N139, N141, N142, N143, N145, N146, N147, N148, N149, N150, N151, N152, N153, N154, N155, N156 | N132, N138, N140, N144 | R12 (B) | S01, S02 (end of chapter) |

Totals: A 15 · B 151 (142 NOTE + 9 RECAP) · C 19 (14 NOTE + 3 RECAP + 2 SKIP) = 185. No section is empty; §4.1 has no A item and is its own short notebook section (N01 tagged C01, N02 tagged C06). Rows of §4.9 that belong to earlier A blocks are listed under §4.9 here and taught in C04 (N82–N84, angular momentum) and C05 (N103–N105, pitot, stagnation pressure, orifice), as §1 says.

## 4. Prerequisites needing primers (concept or tool | needed by | why it is not A/B/RECAP)
Rows marked **gloss** are one plain sentence where the item is used; all others are full 📎 primers (one or two plain
sentences, what it means here, a 2–4-line runnable demo with easy numbers). New primers continue at **P111**. Items
already in `knowledge/primers.md` get a one-line reminder naming the earlier primer (listed at the end).

| Concept or tool | Needed by | Why it is not A/B/RECAP |
|---|---|---|
| continuity of a function and the "small ball" argument (ε–δ reasoning in words) | C02 (D02, N07) | pre-book analysis; the localisation lemma is the B item, this is the tool it rests on |
| product rule for a divergence $\nabla\cdot(\rho\mathbf u)=\mathbf u\cdot\nabla\rho+\rho\nabla\cdot\mathbf u$ | C02 (D03), C06 (D07), C10 (D20) | ch03 used it as a **gloss** in D24; here it drives three derivations — promote to a primer with a sympy demo |
| tensor divergence over the **first** index, $\partial\tau_{ij}/\partial x_i$, and why it differs from $\partial\tau_{ij}/\partial x_j$ for a non-symmetric τ | C06 (D07), C08 (D11), C10 (D19) | ch02 `tensor_divergence` defaulted to the second index; the convention is new to this chapter |
| curl of a product $\nabla\times(\chi\nabla\psi)=\nabla\chi\times\nabla\psi+\chi\nabla\times\nabla\psi$ | C03 (D04), N14 | pre-book vector calculus — **gloss** with a sympy check |
| ∇·(∇×A) = 0 and ∇×∇φ = 0 | C03 (D04) | ch02 D26 corollary and Exercise 2.19 — **gloss** reminder |
| momentum flux: ρu carried through a surface at rate (u − b)·n, a vector per area per time | C04 (D05), N29 | new physics idea of CV momentum budgets |
| Newton's third law (force on the body = − force on the fluid) | C04 (N28) | pre-book mechanics — **gloss** inside the ⚠️ callout |
| torque, moment arm r × F, moment of inertia of a cube ρh⁵/6 | C07 (D08), C04 (N83, N84) | pre-book mechanics not primed in ch01–ch03 |
| side pressure force on a slowly widening tube (mean pressure × projected area) | C05 (D06) | geometric step of Ex. 4.2 — **gloss** with a sketch |
| sympy `expand`, `series`/`removeO`, `collect`, `subs` (dropping (ds)² terms; collecting coefficients) | C05 (D06), C15 (D30), C09 (D15 check) | P40 primed `diff`/`simplify` only |
| isotropic 4th-order tensor: the only rotation-invariant forms are products of two δ's (λ, μ, γ) | C07 (D09) | cited result (Aris); shown with a numeric rotation test, not proved |
| δ substitution in products: $\delta_{im}\delta_{jn}S_{mn}=S_{ij}$, $\delta_{mn}S_{mn}=S_{mm}$ | C07 (D09, D10), C08 (D11) | ch02 C01 taught δ; the two-δ products are new — **gloss** with `np.einsum` (P62 reminder) |
| deviatoric (traceless) part of a tensor $S_{ij}-\tfrac13S_{mm}\delta_{ij}$ | C07 (D10), C10 (D23) | new tensor idea |
| Schwarz's theorem: mixed partial derivatives commute (for smooth fields) | C08 (D12), C12 (D26) | pre-book calculus |
| curl of a curl $\nabla\times(\nabla\times\mathbf u)=\nabla(\nabla\cdot\mathbf u)-\nabla^2\mathbf u$ | C08 (D13), C12 (N97) | built from ε–δ (ch02 D09 reminder); the identity itself is new |
| vector Laplacian acting component by component (Cartesian only) | C08 (D12), C09 (D16) | **gloss**; the curvilinear forms are Appendix B (N65) |
| derivative of a rotating unit vector $d\mathbf e/dt=\boldsymbol\Omega\times\mathbf e$ (cone construction) | C09 (D14, D15) | new; the rotating frame of ch03 P103 was one instant only |
| product rule for a cross product $\frac{d}{dt}(\mathbf a\times\mathbf b)=\dot{\mathbf a}\times\mathbf b+\mathbf a\times\dot{\mathbf b}$ | C09 (D15) | pre-book vector calculus |
| latitude, Earth's rotation rate Ω = 7.292×10⁻⁵ rad/s, "up" and the local vertical component Ω sin φ | C09 (N61, N62, N64) | physics vocabulary for the climate hooks; f = 2Ω sin φ is named as a forward pointer only |
| constant-acceleration kinematics $s=\tfrac12at^2$ | C09 (D17) | pre-book — **gloss** |
| conservative force and its potential (work independent of path, $\mathbf g=-\nabla\Phi$) | C04 (N27), C09 (D18), C11 (D24) | pre-book mechanics; the book equation (4.18) is the B item |
| power of a force (f·u) and heat flux through a surface (q·n) | C10 (D19) | physics vocabulary for the energy budget |
| completing the square for tensors; Frobenius sum of squares $A_{ij}A_{ij}\ge0$ | C10 (D23) | new algebraic move |
| quotient rule $D(1/\rho)/Dt=-(1/\rho^2)D\rho/Dt$ | C10 (D22), N80 | P38 product rule reminder — **gloss** |
| chain rule for $D(\tfrac12u_j^2)/Dt=u_jDu_j/Dt$ | C10 (D21) | P49/P91 reminder — **gloss** |
| a·(a × b) = 0 (triple product with a repeated vector) | C11 (D25) | ch02 D09 — **gloss** |
| a function with zero gradient everywhere is constant in space (connected region), so it can depend on t only | C11 (D25), C12 (D26) | pre-book — **gloss** |
| barotropic fluid (ρ a function of p only) and the pressure function ∫dp/ρ(p) | C11 (D24) | the B item N87; its tool FTC with a variable upper limit is P84 (reminder) |
| order-of-magnitude scaling: derivatives ≈ (change)/(length over which it happens), U/l, U/l² | C13 (D27, N108, N112), C15 (D30, N134–N137) | new estimation tool used by two A blocks |
| reduced gravity g′ = gΔρ/ρ₀ and buoyancy of a light parcel | C13, N137 | physics vocabulary (N² is ch01 C51, reminder) |
| level set moving in time η(x, t) = 0 and its normal speed −(∂η/∂t)/∣∇η∣ | C14 (D29) | P75 primed static level sets; the moving case is new |
| scaled (dimensionless) variables and the chain rule $\partial/\partial t=\Omega\,\partial/\partial t^*$, $\nabla=\nabla^*/l$ | C15 (D30) | new; P59 primed scaling of base units, not of equations |
| `dataclasses` results with named fields (`MassBudget`, `MomentumBudget`, `Scales`) | C01, C04, C15 | Python tool (P23 primed dicts only) |
| `np.cross` with broadcasting over arrays of vectors | C09 (from scratch), C11 | Python tool — **gloss** (P77 broadcasting reminder) |
| `scipy.integrate.tplquad` | N07 demo | Python tool — **gloss** (P87 primed `quad`/`dblquad`) |

Reminders only (already primed): partial derivative P25 · chain rule P49 and along a path P91 · Taylor P26 and
multivariable P98 · definite integral P27 · trapezoid rule P37 · product rule P38 · sympy P40 · `solve_ivp` P31/P94 ·
lambda P29 · `assert np.allclose` P15 · `animate` P16 · `slider_figure` P17 · `show_viz` P18 · `np.einsum` P62 ·
orders of smallness P68 · FTC P84 · mean-value theorem P85 · line integral round a loop P86 · `quad`/`dblquad` P87 ·
cylindrical unit vectors P88 · rigid-body velocity P101 · rotating frame P103 · signed swept volume P110 ·
frames of reference P96 · level sets P75 · `meshgrid` P76 · broadcasting P77 · `contour/quiver/streamplot` P78 ·
eigenvalues P80 (principal stresses) · right-hand rule P74 · log–log slope P13 · plotly 3-D P41/P64 · Newton II P09 ·
net pressure force P28 · internal and kinetic energy P32 · Gibbs free energy P50 (for Helmholtz f) · perfect gas P33 ·
`np.where` P46 · `brentq` P108.

## 4b. Derivations written out (parsed by tools: ID first, CORE id in a column, ★★★ for hard, explainer slugs backticked in the LAST column)
Written out because (a) the result belongs to an A item, or (b) the book never writes it out and the lesson needs it
(marked "(b) book never writes it out"). One small move per step; the book's skipped moves (analysis §2b) are filled in
and the notebook says so. Analysis §2b numbers in brackets (`a-Dnn`). 30 rows.

| ID | Result (Eq.) | CORE | Difficulty | Steps | Tools used | Traps | Shown in |
|---|---|---|---|---|---|---|---|
| D01 | mass for an arbitrarily moving CV: (4.1) → (4.2) → (4.3) + (4.4) → (4.5) [a-D01, a-D02] | C01 | ★★ | 9 | RTT (3.35) (recap ch03 C15), material volume b = u (R01), coincident control volume (N06), signed b·n (P110 reminder) | why ∫∂ρ/∂t over V and V* are equal (same region, same integrand, same instant) while the d/dt∫ρ terms are not; u and b in the same frame; u − b is the velocity relative to the moving wall; checks b = u → (4.1), b = 0 → fixed CV | notebook · `control_volume_budgets` |
| D02 | continuity equation (4.2) → (4.6) → (4.7) with the localisation lemma [a-D03] | C02 | ★★ | 8 | Gauss (R02, ch02 C14), localisation lemma (N07) with continuity of a function (primer) | V must be *every* material volume (that licenses the lemma); ρu must be continuously differentiable for Gauss; the lemma needs a continuous integrand (the book's one-sentence proof is filled in) | notebook |
| D03 | (4.7) → (4.8) → (4.9) → (4.10): continuity with D/Dt; incompressible ⇒ ∇·u = 0 [a-D04] | C02 | ★ | 5 | product rule for ∇·(ρu) (primer), D/Dt (3.5) (recap ch03 C02) | dividing by ρ needs ρ > 0; Dρ/Dt = 0 is weaker than ρ = const (stratified example N10); the §4.11 re-display of (4.9) is not the incompressible condition | notebook |
| D04 | 2-D stream function from $\rho\mathbf u=\nabla\chi\times\nabla\psi$ with χ = −z: $\rho u=\partial\psi/\partial y$, $\rho v=-\partial\psi/\partial x$, and $\psi_2-\psi_1$ = flux per unit depth — (b) book never writes it out (Exercise 4.8) [a-D07] | C03 | ★★ | 8 | curl of a product (gloss), ∇·∇× = 0 (gloss), cross product with e_z (ch02 C08 reminder), chain rule along a curve dψ = ∇ψ·ds (P91 reminder), line integral (P86 reminder) | sign convention (many GFD texts use u = −∂ψ/∂y); orientation of the normal along the connecting curve; flux per unit depth has units m²/s (ρ = 1) | notebook · `stream_function_spacing` |
| D05 | momentum for an arbitrarily moving CV: (4.13) → (4.14) → (4.15) → (4.16a–d) → (4.17) [a-D08] | C04 | ★★ | 8 | RTT component by component (gloss), D01's coincidence move, momentum flux (primer) | (4.15)'s printed "= 0" is spurious; ρu(u − b)·n is a vector times a scalar; the forces are on the *fluid* (drag sign, N28); f is the traction on A*, including pressure | notebook · `control_volume_budgets` |
| D06 | Bernoulli (4.19) from the stream-tube element of Ex. 4.2 [a-D10] | C05 | ★★ | 11 | first-order Taylor (P26), orders of smallness (P68), side pressure force (gloss), sympy `expand`/`series` (primer) | the extra pressure force on the conical side (p + ½∂p/∂s ds)(∂A/∂s)ds; the mean area in the gravity term; U∂U/∂s ds = d(U²/2); constant along *one* streamline only; the book's statement "½ρU² + gz + p/ρ" is dimensionally inconsistent | notebook · `which_bernoulli` |
| D07 | Cauchy's equation: (4.14) → (4.20a, b) → (4.21) → (4.22) → (4.23) → (4.24) [a-D13, a-D14] | C06 | ★★ | 10 | Gauss per component (R02), first-index tensor divergence (primer), localisation (N07), product rule for a divergence (primer) | contract the **first** index of τ with n_i (f_j = n_iτ_ij); the book's prose writes ∂τ_ij/∂x_j — harmless only because τ is symmetric, which is proved *after* (4.24); the bracket in (4.23) is continuity | notebook |
| D08 | stress symmetry (4.25) from the torque on a shrinking cube — (b) book never writes it out (Exercise 4.30) [a-D15] | C07 | ★★ | 8 | torque and moment of inertia of a cube (primer), orders of smallness (P68) | normal stresses and body forces give higher-order (or zero) torque; the moment of inertia ρh⁵/6 shrinks faster than the torque (τ₁₂ − τ₂₁)h³, so the angular acceleration 6(τ₁₂ − τ₂₁)/(ρh²) blows up unless τ₁₂ = τ₂₁; body couples are the exception | notebook · `newtonian_stress_lab` |
| D09 | Newtonian constitutive law (4.28) → (4.29) → (4.30) → (4.31) [a-D17] | C07 | ★★★ | 13 | σ depends on S only (R05), isotropic 4th-order tensor (primer, cited form + rotation test), δ substitution in products (gloss), symmetric S (ch03 C06), sympy check cell | "γ = μ" really names μ + γ = 2μ (only the sum acts on a symmetric S); δ_inδ_jmS_mn = S_ji = S_ij; σ must vanish when S = 0 (no u, no R); λ is not the bulk viscosity | notebook · `newtonian_stress_lab` |
| D10 | (4.31) → (4.37): add and subtract (2μ/3)S_mmδ_ij; bulk viscosity μ_v = λ + ⅔μ, Stokes (4.36) — (b) book states the result [a-D19] | C07 | ★ | 4 | deviatoric part (primer), δ_ii = 3 (ch02 reminder) | the traceless (shear) part is S − ⅓S_mmδ, not S; μ_v = λ + ⅔μ (not λ); the Stokes assumption sets μ_v = 0, not λ = 0 | notebook · `newtonian_stress_lab` |
| D11 | Navier–Stokes (4.24) + (4.37) → (4.38) [a-D20] | C08 | ★★ | 7 | δ substitution (gloss), first-index tensor divergence (primer), product rule (P38) | ∂(pδ_ij)/∂x_i = ∂p/∂x_j; keep a variable μ(T) *inside* the derivative; 2μS_ij = μ(∂u_i/∂x_j + ∂u_j/∂x_i) | notebook · `navier_stokes_term_balance` |
| D12 | (4.38) → (4.39a) → (4.39b): constant viscosity, then incompressible — (b) book states the results [a-D21] | C08 | ★★ | 7 | Schwarz's theorem (primer), vector Laplacian (gloss) | ∂/∂x_i(∂u_i/∂x_j) = ∂/∂x_j(∇·u) needs smooth u; coefficient bookkeeping μ + μ_v − ⅔μ = μ_v + ⅓μ; needs μ, μ_v uniform; ∇·u = 0 removes the last term | notebook · `navier_stokes_term_balance` |
| D13 | viscous force three ways (4.40): $\mu\nabla^2\mathbf u=2\mu\partial S_{ij}/\partial x_i=-\mu\nabla\times\boldsymbol\omega$ — (b) book never writes it out (Exercise 4.38) [a-D22] | C08 | ★★ | 9 | ε–δ (2.19) (recap ch02 D09), curl of a curl (primer), vorticity ω = ∇×u (recap ch03 C10) | 2∂S_ij/∂x_i = ∇²u_j + ∂_j(∇·u) — equality needs ∇·u = 0; index order in ε_jik (a swap flips the sign); solid-body rotation has ω ≠ 0 but ∇×ω = 0, so no viscous force ("the paradox") | notebook · `navier_stokes_term_balance` |
| D14 | velocity in a translating–rotating frame (4.42), with $d\mathbf e'_i/dt=\boldsymbol\Omega\times\mathbf e'_i$ [a-D23] | C09 | ★★ | 7 | derivative of a rotating unit vector (primer), product rule (P38), cross product (ch02 C08 reminder) | the rotating frame's basis vectors are functions of time (the whole difference from ch03's Galilean frame); u′ is the rate of the *components* only; Ω × x′ uses the rotating-frame position | notebook · `rotating_frame_coriolis` |
| D15 | acceleration in a noninertial frame (4.42) → (4.43) → (4.44) — (b) book never writes it out (Exercise 4.42) [a-D24] | C09 | ★★★ | 12 | D14 applied twice, product rule for a cross product (primer), vector triple product (recap ch02 D09), sympy check cell | the **two** Ω × u′ (one from the turning basis, one from d(Ω × x′)/dt) make the Coriolis factor 2; d/dt(Ω × x′) = Ω̇ × x′ + Ω × (u′ + Ω × x′); acceleration +2Ω × u′ (4.43) vs force −2Ω × u′ (4.45); "centripetal" Ω × (Ω × x′) vs "centrifugal" −Ω × (Ω × x′) | notebook · `rotating_frame_coriolis` |
| D16 | Navier–Stokes in a noninertial frame (4.44) + (4.39b) → (4.45) [a-D25] | C09 | ★★ | 6 | vector Laplacian (gloss), gradient of a scalar is frame-independent at the coincidence instant (gloss) | ∇²(U + Ω × x′) = 0 because it is linear in x′ — needs U, Ω uniform in space; derived from the incompressible, constant-μ (4.39b) (the book says "(4.39)"); moving the frame terms to the right flips their signs | notebook · `rotating_frame_coriolis` |
| D17 | Coriolis projectile from the pole: deflection $\Omega ut^2$, angle Ωt — (b) the book asserts it [a-D26] | C09 | ★ | 5 | constant-acceleration kinematics (gloss), cross product direction (right-hand rule P74) | valid for small Ωt only (path nearly straight); the Coriolis force does no work (⟂ u); Ω up in the NH ⇒ deflection to the right | notebook · `rotating_frame_coriolis` |
| D18 | centrifugal term $-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')=\Omega^2R\,\mathbf e_R$ and its potential $-\tfrac12\Omega^2R^2$; effective gravity — (b) book never writes it out (Exercise 4.43) [a-D27] | C09 | ★ | 5 | vector triple product (recap ch02 D09), conservative force (primer), cylindrical gradient (gloss) | −Ω × (Ω × x) points *away* from the axis; R is the distance from the axis, not from the Earth's centre; the potential's sign | notebook · `rotating_frame_coriolis` |
| D19 | total-energy equation (4.46) → (4.47) → (4.49)–(4.52) → (4.53) [a-D30] | C10 | ★★ | 8 | RTT (recap ch03 C15), Gauss three times (R02), first-index contraction (primer), localisation (N07), power of a force (primer) | the book prints dA inside the volume integrals of (4.51) (must be dV); u_j² is a sum over j; heat flowing out is a loss (the minus sign on q·n) | notebook |
| D20 | (4.53) + (4.54) → (4.55): total energy in D/Dt form — (b) book never writes it out (Exercise 4.45) [a-D31] | C10 | ★★ | 7 | product rule for a divergence (primer), the (4.23) move with E = e + ½u_j², δ substitution (gloss) | subtract E × continuity (not E × ∇·u); −p∂u_j/∂x_j comes from the δ_ij in τ_ij | notebook |
| D21 | mechanical-energy equation (4.56) — (b) book never writes it out (Exercise 4.46) [a-D32] | C10 | ★★ | 6 | chain rule for D(½u_j²)/Dt (gloss), (4.27) split (N41) | start from Cauchy (4.24) (from the flux form (4.22) you must first subtract u_j × continuity); u_j∂τ_ij/∂x_i = −u_j∂p/∂x_j + u_j∂σ_ij/∂x_i | notebook · `viscous_dissipation_heating` |
| D22 | internal-energy equation (4.55) − (4.56) → (4.57) [a-D33] | C10 | ★★ | 8 | quotient rule (gloss), (4.8) (N09), σ:R = 0 for symmetric σ (R07) | the force-work terms cancel exactly; −(p/ρ)∂u_j/∂x_j = −p Dv/Dt via (4.8); σ_ij∂u_j/∂x_i = σ_ijS_ij *only* because σ is symmetric | notebook · `viscous_dissipation_heating` |
| D23 | dissipation as a sum of squares, ε ≥ 0 (4.58) — (b) book never writes it out (Exercise 4.47) [a-D34] | C10 | ★★ | 8 | completing the square for tensors (primer), deviatoric part (primer), δ_ijδ_ij = 3 | the cross term needs δ_ijδ_ij = 3 (not 1); ε ≥ 0 iff μ ≥ 0 and μ_v ≥ 0; incompressible flow gives ε = 2νS_ijS_ij | notebook · `viscous_dissipation_heating` |
| D24 | Bernoulli function: (4.66) + (4.67) + Lamb identity (4.68) → (4.69) — (b) book never writes (4.68) out (Exercise 4.50) [a-D38] | C11 | ★★ | 10 | ε–δ (recap ch02 D09), FTC with a variable upper limit (P84 reminder), chain rule (P49), barotropic pressure function (N87), conservative force (primer) | the sign of u × ω; the prime in (4.67) is a dummy integration variable; barotropy is what makes dp/ρ an exact differential (a baroclinic fluid has no single B) | notebook · `which_bernoulli` |
| D25 | steady flow ∇B = u × ω (4.70) ⇒ B constant on streamlines and vortex lines (4.71), everywhere if irrotational (4.72) [a-D39] | C11 | ★ | 5 | triple product with a repeated vector (gloss), zero gradient ⇒ constant (gloss) | constant along each line but different constants on different lines; (4.72) needs ω = 0 everywhere in a connected region | notebook · `which_bernoulli` |
| D26 | unsteady Bernoulli (4.73) → (4.74) → (4.75) [a-D40] | C12 | ★★ | 8 | Schwarz's theorem to swap ∂_t and ∇ (primer), zero gradient ⇒ function of t only (gloss), velocity potential (R09) | the gauge change must be φ_new = φ − ∫B dt′ (the book prints +, which doubles B); the "constant" may depend on t before it is absorbed; needs an irrotational, barotropic flow | notebook · `which_bernoulli` |
| D27 | Boussinesq momentum: (4.39b) − hydrostatics → (4.84) → (4.85) → (4.86) [a-D44, a-D46] | C13 | ★★ | 9 | hydrostatic base state (R10), order-of-magnitude scaling (primer), reduced gravity (primer) | why ρ′ is kept next to g but dropped in ρ Du/Dt (g ≫ ∣Du/Dt∣, so ρ′g/ρ₀ competes with the other forces); replacing ρ_s(z) by ρ₀ needs a small base-state variation; primes are perturbations here | notebook · `boussinesq_buoyancy` |
| D28 | Boussinesq heat equation (4.87) → (4.88) → (4.89) [a-D47] | C13 | ★★ | 8 | perfect gas p = ρRT, R = C_p − C_v, α = 1/T (recap ch01 C40, C42, C48), Fourier (recap ch01 C11), order-of-magnitude scaling (primer) | the sign bookkeeping that moves ρ(C_p − C_v)DT/Dt to the left and turns C_v into C_p; for liquids the argument is different (p∇·u small and C_p ≈ C_v); dropping ρε needs low Ec (N151) | notebook · `boussinesq_buoyancy` |
| D29 | kinematic boundary condition (4.90) → (4.91) → (4.92) → (4.93) [a-D49] | C14 | ★★ | 8 | moving level set and its normal speed (primer), chain rule along a path (P91 reminder), unit normal ∇η/∣∇η∣ (P75 reminder) | only the normal part of the surface velocity u_s is defined (and needed); n = ∇η/∣∇η∣ points toward increasing η; the linearised condition is exact only to O(ka) | notebook · `kinematic_free_surface` |
| D30 | dimensionless Navier–Stokes (4.100) → (4.101) [a-D54] | C15 | ★★ | 8 | scaled variables and the chain rule (primer), sympy `subs`/`collect` (primer) | ∂/∂t = Ω∂/∂t*, ∇ = ∇*/l; p∞ drops because only ∇p enters; divide by ρU²/l (the advective scale); other pressure scalings (μU/l, ρgl) give products of the same groups | notebook · `dynamic_similarity_models` |

Counts: ★ 5 (D03, D10, D17, D18, D25) · ★★ 23 · ★★★ 2 (D09 in `newtonian_stress_lab`, D15 in `rotating_frame_coriolis`;
each has a sympy check cell that re-runs the construction: D09 builds K_ijmn from δ's, contracts with a symbolic
symmetric S and shows only μ + γ survives, then checks isotropy under a symbolic rotation about z; D15 differentiates
x = X(t) + R(t)x′(t) with a sympy rotation matrix R(t) and collects the five terms). Written out under rule (b): D04,
D08, D10, D12, D13, D15, D17, D18, D20, D21, D23, D24 (12); under rule (a) only: the other 18. Analysis rows merged into
one written-out derivation: a-D01 + a-D02 → D01, a-D13 + a-D14 → D07, a-D44 + a-D46 → D27.

## 4c. Derivations demoted to statements (first column is the A parent in bold, e.g. **C20** — never a bare ID or a D id: the parser reads a bare first-cell ID as an item and blanks its tier)
Results **given, not derived**, inside the named A block (a paragraph, the equation, a number). 27 rows (26 demoted +
1 recap row of analysis §2b).

| A parent | Analysis §2b item | Result stated (Eq.) | Stated in (B item) | Why not written out |
|---|---|---|---|---|
| **C03** | a-D05 ρu = ∇χ × ∇ψ satisfies (4.11); ρu ⟂ ∇χ, ∇ψ | $\nabla\cdot(\nabla\chi\times\nabla\psi)=0$, $\rho\mathbf u\cdot\nabla\chi=\rho\mathbf u\cdot\nabla\psi=0$ | N14, N15 | two ch02 identities (∇·∇× = 0, ∇×∇ = 0) in one line each; the existence of such χ, ψ (Clebsch) is only asserted by the book — said so; the 2-D case is derived in D04 |
| **C03** | a-D06 stream-tube flux $\dot m=(b-a)(d-c)$ | $\dot m=\oint_C\chi\,d\psi=b(d-c)+a(c-d)$ | N16 | the book writes the chain; the leg-by-leg reading (dψ = 0 on ψ = const legs) is one sentence plus a numeric patch flux |
| **C04** | a-D09 Ex. 4.1 wake drag | $F_D/l=\rho\int U(U_\infty-U)dy$ | N29 | the book writes it; shown as C04's worked number with the four faces' fluxes as bars and the side leakage computed (0.354 m²/s), not as symbolic algebra |
| **C04** | a-D11 Ex. 4.3 bore speed | $U=\sqrt{gh_{out}(h_{in}+h_{out})/2h_{in}}\approx\sqrt{gh}$ | N31 | example of C04's law; stated with the moving CV and one number; the √(gh) limit is Ch. 7's |
| **C04** | a-D12 Ex. 4.4 rocket | $M\,d^2z_R/dt^2=-V_e\,dM/dt-Mg+F_S$, Tsiolkovsky $\Delta b=V_e\ln(M_0/M_1)$ | N32 | example; stated with a trajectory figure from `solve_ivp` |
| **C04** | a-D36 angular momentum for a stationary CV (4.65) | $\frac{d}{dt}\int(\mathbf r\times\rho\mathbf u)dV+\int(\mathbf r\times\rho\mathbf u)(\mathbf u\cdot\mathbf n)dA=\int\mathbf r\times\rho\mathbf g\,dV+\int\mathbf r\times\mathbf f\,dA$ | N83 | the RTT with F = r × ρu, the same move as D05 — one sentence |
| **C04** | a-D37 Ex. 4.6 sprinkler | $M=2a\rho AU^2\cos\alpha$ | N84 | example with one number (0.87 N m) |
| **C05** | a-D43 pitot and orifice speeds | $\lvert\mathbf u\rvert_1=\sqrt{2(p_2-p_1)/\rho}$, $u=\sqrt{2gh}$ | N103, N105 | two lines of (4.19) with u₂ = 0 or p = p_atm; stated with numbers |
| **C07** | a-D16 static stress −pδ (4.26) | $\tau_{ij}=-p\,\delta_{ij}$ | N40 | one-line reason (δ is the only isotropic 2nd-order tensor, ch02 N42) and the sign convention |
| **C07** | a-D18 trace of (4.31) → (4.32)–(4.34) | $p-\bar p=(\tfrac23\mu+\lambda)\nabla\cdot\mathbf u$ | N45, N46, N47 | a contraction with δ_ii = 3; stated with the ⚠️ "Section 3.6" → §3.4 note |
| **C09** | a-D28 Ex. 4.5 pump in cylindrical components | the three momentum equations with $\rho[2\Omega_zu_\varphi+\Omega_z^2R]$, $\rho[-2\Omega_zu_R]$ | N65 | Appendix-B operators; printed by sympy (`ch04.rotating_pump_equations`) with the rotation terms highlighted |
| **C10** | a-D29 energy for an arbitrarily moving CV (4.48) | (4.48) | N68 | the D01/D05 coincidence algebra a third time — one sentence |
| **C10** | a-D35 entropy equation (4.62)–(4.63) | $\frac{Ds}{Dt}=\frac1\rho\nabla\cdot(\frac kT\nabla T)+\frac{k}{\rho T^2}\lvert\nabla T\rvert^2+\frac\varepsilon T$ | N80, N81 | the book shows the moves (combine with Gibbs, quotient rule); stated with the second-law payoff μ, μ_v, k ≥ 0 and a numeric sign demo |
| **C11** | a-D41 energy Bernoulli (4.76) → (4.77) → (4.78) | $h+\tfrac12\lvert\mathbf u\rvert^2+gz$ = constant on streamlines | N94, N95, N96 | B inside C11; stated with the stagnation-temperature number; Ch. 15 derives its compressible use |
| **C12** | a-D42 viscous irrotational Bernoulli (4.79)–(4.83) | $\int_1^2\frac{\partial\mathbf u}{\partial t}\cdot d\mathbf s+(\tfrac12\lvert\mathbf u\rvert^2+gz+\frac p\rho)_2=(\ldots)_1$ | N97, N98, N99, N100, N101 | follows from D13 and D24 in two sentences (ω = 0 kills both the Lamb and the viscous term); stated with the U-tube number |
| **C13** | a-D45 Boussinesq continuity scaling | $\frac{(1/\rho)D\rho/Dt}{\nabla\cdot\mathbf u}\sim\alpha\,\delta T$, $L\ll c^2/g$ | N108 | the book writes the scaling chain; stated with the validity table (c²/g ≈ 11.8 km for air) |
| **C14** | a-D48 pillbox interface conditions | continuous $\rho\mathbf u\cdot\mathbf n$, $n_i\tau_{ij}$, $k\,\partial T/\partial n$ | N116 | stated with the l → 0 picture and two worked cases (composite wall, two-fluid Couette) |
| **C14** | a-D50 ★★★ Laplace jump from the cap force balance (4.97) + (4.98) → (1.5) | $(F_p)_z=-\pi\Delta p\sqrt{2R_1\zeta}\sqrt{2R_2\zeta}$, $(F_{st})_z=\pi\sigma\sqrt{2R_1\zeta}\sqrt{2R_2\zeta}(\frac1{R_1}+\frac1{R_2})$ ⇒ $\Delta p=\sigma(\frac1{R_1}+\frac1{R_2})$ | N125, N126, N127 | the result is SEEN (ch01 N07) and the book writes the construction; C14 states both forces with numeric checks (dblquad; exact-ellipse line integral converging at O(ζ)) and the one-line balance |
| **C14** | a-D51 ★★★ Ex. 4.7 meniscus | $h^2=\frac{2\sigma}{\rho g}(1-\sin\theta)$ and the closed-form profile x(ζ) | N129 | worked example; the closed form is verified numerically (ODE vs formula) instead of derived; the book's dropped minus sign is noted |
| **C14** | a-D52 free energy ⇒ interface contracts | $F=\rho_1V_1f_1+\rho_2V_2f_2+A\sigma$ | N121, N122, N123 | thermodynamic framing; stated with the spheroid least-area figure |
| **C15** | a-D53 (RECAP) sphere-drag groups (4.99) | $F_D\rho/\mu^2=(F_D/\rho U^2d^2)\,\mathrm{Re}^2$ | R12 | recap of ch01 D28 with two repeating sets (`pi_groups`) |
| **C15** | a-D55 groups as force ratios (4.102)–(4.105), (4.111), (4.117)–(4.119) | $\mathrm{St}=\Omega l/U$, $\mathrm{Re}=\rho Ul/\mu$, $\mathrm{Fr}=U/\sqrt{gl}$, $M=U/c$, We, Bo, Ca | N134, N135, N136, N137, N147, N153, N154, N155 | each is one ratio chain the book prints; the order-of-magnitude primer explains the move once |
| **C15** | a-D56 purely oscillatory body | U = lΩ ⇒ St = 1, Re = Ωl²/ν | N141 | a substitution |
| **C15** | a-D57 compressibility ∝ M² (4.110) | $\nabla^*\cdot\mathbf u^*=-[U^2/c^2]\frac1{\rho^*}\frac{Dp^*}{Dt^*}$ | N146 | the D30 scaling move once more with dp = c²dρ; stated with M² = 0.09 at M = 0.3 |
| **C15** | a-D58 enthalpy form (4.112) | $\rho\frac{Dh}{Dt}=\frac{Dp}{Dt}+\rho\varepsilon+\nabla\cdot(k\nabla T)$ | N148 | the book calls it "a mild revision" and never writes it out; the lesson needs only the result — stated with a sympy one-liner (`energy_forms_sym`) proving (4.60) ⇔ (4.112) |
| **C15** | a-D59 dimensionless energy → Ec, Pr (4.114)–(4.116) | coefficients Ec, Ec/Re, 1/(Pr Re) | N150, N151, N152 | the D30 move applied to (4.112); printed by `nondimensional_energy_coefficients` |
| **C15** | a-D60 Ex. 4.8 ship model | $U_m=U_p\sqrt{l_m/l_p}$, wave drag ×λ³ | N156 | example; bar chart with numbers; the Re mismatch λ^{3/2} is the lesson |

## 5. Interactive explainers (5–10 + backup)
Nine explainers on the A items where manipulation or motion teaches most; one backup. Each has the required live
**Explain** tab ("Explanation & interpretation", numbered sections computing every number on screen with the reader's
settings, ending in "Reading the current setting"), a synced **Code** tab, a **Derivation** tab for its D ids, a 4–8-step
walkthrough, ≥ 3 check questions and ≥ 2 more depth features. Physics mirrors scalar-callable `fluidpy` functions (parity
rows). Colours: inflow blue, outflow orange, storage purple, force rose; pressure orange, viscous rose, inertia teal,
gravity/buoyancy blue, Coriolis amber, centrifugal purple; mechanical energy teal, heat rose.
Not duplicated from earlier chapters: ch03 `reynolds_transport_cv` (the RTT itself) is *linked from* E1, which adds
forces and applications; ch02 `cauchy_traction_principal_axes` (traction on a plane) is linked from E3, which adds the
constitutive law; ch03 `galilean_frames_cylinder` (translating observer) is linked from E5, which adds rotation; ch01
`buckingham_pi_machine` (Π groups) is linked from E9, which adds the equations and model testing.

### E1 · control_volume_budgets
- **A:** C01, C04 (also shows N03 (4.1) material CV, N06 the coincident CV, N26 body vs surface forces, N28 drag sign,
  N29 Ex. 4.1 wake, N31 Ex. 4.3 bore, N32 Ex. 4.4 rocket, N84 Ex. 4.6 sprinkler as a mode)
- **Confusion removed:** "to find a force you must know the flow everywhere" — and the sign mess of fluxes through a
  moving wall: (u − b)·n is the flux *relative to the wall*, and the force in the budget is the force *on the fluid*.
- **Why interactive:** a budget is bookkeeping per face; dragging the box height H, moving the CV at b and switching
  scenarios while the face bars and the residual update shows the terms trade places yet always sum to zero — the
  moving-CV case (bore) becomes steady only when b matches the wave speed, which the reader finds by dragging.
- **Stage:** (1) the scene with the control volume drawn on it (wake behind a bar with the U(y) profile; a bore with
  the CV riding at b; a jet striking a plate; a rocket) and flux arrows per face coloured by sign; (2) waterfall bars:
  storage d/dt∫ + outflow − inflow per face (mass, then momentum), body and surface forces, total = 0; (3)
  (hidePortrait) the result curve (F_D/l vs H converging; U vs h_out/h_in with the √(gh) ghost; rocket speed vs t with
  the Tsiolkovsky ghost).
- **Controls:** mode: wake (Ex. 4.1) / bore (Ex. 4.3) / jet on a plate (ours) / rocket (Ex. 4.4) · box height H or CV
  speed b · deficit Δ or h_out/h_in · budget: mass / momentum · time (transport, rocket and bore).
- **Equations shown (live):** (4.1) $\frac{d}{dt}\int_{V(t)}\rho\,dV=0$; (4.5)
  $\frac{d}{dt}\int_{V^*}\rho\,dV+\int_{A^*}\rho(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=0$; (4.17)
  $\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV+\int_{A^*}\rho\mathbf u(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA$;
  Ex. 4.1 $F_D/l=\rho\int U(U_\infty-U)dy$; Ex. 4.3 $U=\sqrt{gh_{out}(h_{in}+h_{out})/2h_{in}}$; Ex. 4.4
  $M\,d^2z_R/dt^2=-V_e\,dM/dt-Mg+F_S$.
- **Mirrors:** `core.conservation.mass_budget`, `momentum_budget`, `ch04.wake_drag_per_span`, `gaussian_wake`,
  `bore_speed`, `rocket_trajectory`, `ch04.cv_scenario` (§8).
- **Derivations:** D01 (the coincidence step freezes the material CV on top of the fixed one; the b = u / b = 0 checks
  switch the CV type), D05 (the four coincidence equalities light the four bar groups).
- **Depth features:** Explain tab (face by face: area, (u − b)·n, mass flux, momentum flux → storage → forces → the
  residual → the drag or wave speed boxed → reading the current setting), synced Code tab, + linked views (3), term bars
  (waterfall), modes (four scenarios), presets (fixed box / box riding with the bore / material CV b = u; tall vs short
  box), status verdict ("✅ budget closes" / "⚠️ box too short: pressure and shear on the far faces matter"), transport.
- **Follows reference:** `fid_formula_lab.html` (terms that add up to a total, stage chips) with
  `forced_damped_vibrations.html`'s explanation panel; the E7 pattern of ch03 (`reynolds_transport_cv`) for the CV.
- **Aha:** you can weigh a force without touching the body — count what flows in and out of any box around it, relative
  to the box's own walls.

### E2 · stream_function_spacing
- **A:** C03 (also shows N11 ∇·u = 0 holds automatically, N13 (4.11), N16 flux between stream surfaces, N17
  axisymmetric mode, N18)
- **Confusion removed:** "streamlines only show direction" — with equal Δψ between contours, their spacing *is* the
  speed and the flux between two contours is Δψ, whatever path you cross them on.
- **Why interactive:** drag a gate (a line segment) across the flow and watch the flux through it stay equal to
  ψ₂ − ψ₁ however you tilt or bend it; click a point to read u = ∂ψ/∂y, v = −∂ψ/∂x from the contour slope — the
  "spacing = speed" rule is felt, not read.
- **Stage:** (1) ψ contours at equal Δψ over a speed heatmap, with tracer particles moving along them; a draggable
  gate with its flux readout; (2) (hidePortrait) the profile of ψ along the gate and its slope (= normal velocity);
  (3) the equal-flux ribbon between two chosen contours (area proportional to Δψ).
- **Controls:** preset flow: uniform stream / stagnation point ψ = kxy / source + stream (half-body) / cylinder /
  line vortex / shear flow / axisymmetric uniform stream · strength (k, U, m or Γ) · number of contours · gate
  endpoints (drag) · time (tracers).
- **Equations shown (live):** (4.11) $\nabla\cdot(\rho\mathbf u)=0$; (4.12) $\rho\mathbf u=\nabla\chi\times\nabla\psi$;
  $\rho u=\partial\psi/\partial y$, $\rho v=-\partial\psi/\partial x$; $\psi_2-\psi_1=\int_1^2\mathbf u\cdot\mathbf n\,ds$;
  axisymmetric $\rho u_R=-R^{-1}\partial\psi/\partial z$, $\rho u_z=R^{-1}\partial\psi/\partial R$.
- **Mirrors:** `core.streamfunction.velocity_from_streamfunction_2d`, `flux_between_streamlines`,
  `velocity_from_streamfunction_axisym`, `core.streamfunction.streamfunction_preset` (§8), `ch03.cylinder_streamfunction`.
- **Derivations:** D04 (the χ = −z step sets the uniform stream; the flux step sweeps the gate between two contours).
- **Depth features:** Explain tab (ψ at the gate ends → Δψ → the numeric flux along the gate → local spacing and speed
  u ≈ Δψ/Δn → reading the current setting), synced Code tab, + linked views (3), presets (seven flows), inspector (click
  a point: ∂ψ/∂y and −∂ψ/∂x by central differences with the numbers), transport (tracers), status ("↔ fast where
  contours crowd" / "⚠️ stagnation point: contours cross").
- **Follows reference:** `stride_padding_playground.html` (classic presets, the formula with numbers plugged in, click
  a cell to move the window) with `forced_damped_vibrations.html`'s explanation panel.
- **Aha:** a 2-D flow fits in one scalar: its contours are the streamlines, their crowding is the speed, and the flux
  between any two of them is just the difference of their labels.

### E3 · newtonian_stress_lab
- **A:** C07 (also shows N38 (4.25) with the spinning-cube view, N40 (4.26), N41 (4.27), R05 σ depends on S only,
  N42–N44 (4.28)–(4.30), N45–N47 p vs p̄, N49 μ_v, N50 Stokes, N51 (4.37) and (1.3); C06's traction; links to ch02
  `cauchy_traction_principal_axes`)
- **Confusion removed:** "viscous stress depends on velocity" (it depends on the *strain rate*; a rigid rotation feels
  none), and "pressure is the mean normal stress" (only when μ_v = 0 or ∇·u = 0).
- **Why interactive:** set a velocity gradient G with sliders or presets (shear, extension, rotation, expansion) and
  see the element deform, S and R split, and the traction on a draggable plane swing through its normal and shear
  values — toggling μ_v changes p̄ vs p only for expansion. No static figure shows the chain G → S → τ → traction
  responding together.
- **Stage:** (1) a fluid square deforming under G with the traction arrow on a rotatable plane; (2) normal and shear
  stress on the plane vs its angle (Mohr-type curve) with the current angle marked and principal directions; (3)
  (hidePortrait) the 3×3 matrices G, S, τ with the entries coloured by source (−pδ orange, 2μS rose, λS_mmδ purple).
- **Controls:** preset flow: shear / extension / pure rotation / isotropic expansion / custom · strain magnitude
  (s⁻¹) · μ (Pa s) · μ_v (Pa s, 0 = Stokes) · plane angle (drag) · pressure p.
- **Equations shown (live):** (4.25) $\tau_{ij}=\tau_{ji}$; (4.26) $\tau_{ij}=-p\delta_{ij}$; (4.27)
  $\tau_{ij}=-p\delta_{ij}+\sigma_{ij}$; (4.28) $\sigma_{ij}=K_{ijmn}S_{mn}$; (4.29)
  $K_{ijmn}=\lambda\delta_{ij}\delta_{mn}+\mu\delta_{im}\delta_{jn}+\gamma\delta_{in}\delta_{jm}$; (4.31)
  $\tau_{ij}=-p\delta_{ij}+2\mu S_{ij}+\lambda S_{mm}\delta_{ij}$; (4.33) $\bar p=-\tfrac13\tau_{ii}$; (4.37)
  $\tau_{ij}=-p\delta_{ij}+2\mu(S_{ij}-\tfrac13S_{mm}\delta_{ij})+\mu_vS_{mm}\delta_{ij}$.
- **Mirrors:** `core.constitutive.newtonian_stress`, `viscous_stress`, `mean_pressure`, `isotropic_fourth_order`,
  `ch04.cube_spin_acceleration`, `core.constitutive.stress_on_plane` (§8).
- **Derivations:** D08 (a spinning-cube mode: shrink h and watch the angular acceleration blow up unless τ₁₂ = τ₂₁),
  D09 ★★★ (13 steps; the δ-substitution steps light the matrix entries they produce; the "γ = μ" step shows μ + γ acting
  on a symmetric S; the rotation-test step spins the plane), D10 (the add-and-subtract step splits the diagonal colours).
- **Depth features:** Explain tab (G → S and R → each τ entry with numbers → traction on the plane → p̄ vs p → the
  regime reading), synced Code tab, + linked views (3), presets (five flows), inspector (click a matrix entry: its
  arithmetic), status ("🔄 rigid rotation: no viscous stress" / "↔ shear: τ₁₂ = μγ" / "🫧 expansion: p̄ ≠ p unless
  μ_v = 0"), modes (stress lab / spinning cube).
- **Follows reference:** `amplitude_phase_second_order_II_3.html` (several windows on one state, crosshair readouts,
  numbered live derivation in the explanation).
- **Aha:** a Newtonian fluid feels only how fast it is being deformed — spin it rigidly and the viscous stress is zero;
  two numbers, μ and μ_v, are all the constitutive law needs.

### E4 · navier_stokes_term_balance
- **A:** C08 (also shows C06 Cauchy's terms, N52 (4.38) and the 4-vs-5 ledger, N53 (4.39a), N54 (4.40) viscous force
  three ways, N55 Euler (4.41); N107 (4.85) gravity absorbed into p′)
- **Confusion removed:** "all five terms always matter" — in Poiseuille flow only pressure and viscosity balance, in
  potential flow viscosity does no net work though μ ≠ 0, in a decaying vortex the local term balances viscosity; the
  term that matters depends on the flow and the place.
- **Why interactive:** pick an exact solution, click any point and see the five terms of (4.39b) as bars summing to
  zero; scrub time for the unsteady solutions and watch the balance shift from local vs viscous to advective vs
  pressure — a static figure can show one point at one time only.
- **Stage:** (1) the flow field (speed heatmap + streamlines or the profile u(y, t)) with a probe; (2) term bars at the
  probe: local ∂u/∂t (blue), advective (u·∇)u (teal), pressure −∇p/ρ (orange), gravity (grey), viscous ν∇²u (rose),
  residual (black, ≈ 0); (3) (hidePortrait) the viscous force three ways (∇²u, 2∂S/∂x_i, −∇×ω) as three equal arrows,
  or the profile's time history.
- **Controls:** solution: Couette / plane Poiseuille / Stokes' first problem / Taylor–Green / Lamb–Oseen / cylinder
  potential flow (Euler) · ν (log) · probe position (click) · time (transport) · component (x / y).
- **Equations shown (live):** (4.24) $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$; (4.38)
  (variable μ); (4.39b) $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$; (4.40)
  $\mu\nabla^2\mathbf u=2\mu\partial S_{ij}/\partial x_i=-\mu\nabla\times\boldsymbol\omega$; (4.41)
  $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g$.
- **Mirrors:** `core.navier_stokes.ns_incompressible_terms`, `viscous_force_forms`, `core.navier_stokes.exact_solution`
  (§8), `core.vortices.gaussian_vortex`.
- **Derivations:** D11 (steps light the bar each term becomes), D12 (the Schwarz step sets a compressible toy field to
  show the extra term, then ∇·u = 0 removes it), D13 (the ε–δ step draws ∇×ω next to ∇²u; the paradox step sets solid
  body rotation).
- **Depth features:** Explain tab (the solution's formula at the probe → each term with numbers → residual → which
  terms dominate (ratios) → local Reynolds number → reading the current setting), synced Code tab, + linked views (3),
  term bars, transport, presets (six solutions), inspector (click: the finite-difference arithmetic of ∇²u), status
  ("⚖️ pressure ↔ viscous" / "⏳ local ↔ viscous (diffusion)" / "🌀 inviscid: viscous force exactly 0").
- **Follows reference:** `fid_formula_lab.html` (clickable terms, bars adding to a total) with
  `forced_damped_vibrations.html`'s explanation panel.
- **Aha:** Navier–Stokes is a budget of five accelerations; each exact solution is a flow where only two or three of
  them are awake — and viscosity can be present yet do nothing.

### E5 · rotating_frame_coriolis
- **A:** C09 (also shows N56 inertial frame, N57 (4.42), N58 (4.43), N59 (4.44), N60 frame acceleration, N61 Coriolis
  projectile, N62 highs and lows, N63 angular acceleration, N64 centrifugal and effective gravity; links to ch03
  `galilean_frames_cylinder`)
- **Confusion removed:** "Coriolis is a real force pushing the air" — it is the rotating observer's bookkeeping of a
  straight inertial path; and "the Coriolis acceleration deflects to the right" (the *force* −2Ω × u′ does; the
  acceleration term in (4.43) has the opposite sign).
- **Why interactive:** the same ball, one clock, two observers: an inertial view (straight line, turntable spinning
  under it) beside the rotating view (curved path) with the five apparent-force arrows on the ball; dragging Ω, the
  launch speed or the hemisphere and watching deflection = Ωut² appear is the derivation made visible.
- **Stage:** (1) inertial view: turntable (or polar cap) rotating, ball on a straight path; (2) rotating view: the
  curved path with arrows for Coriolis (amber), centrifugal (purple), frame acceleration and angular acceleration, and
  the deflection δ with its Ωut² ghost; (3) (hidePortrait) term bars per unit mass, or the flow-out-of-a-high mode
  (radial outflow + Coriolis arrows turning it clockwise in the NH).
- **Controls:** mode: turntable / Earth's pole / flow out of a high (or into a low) / effective gravity vs latitude · Ω
  (or hemisphere sign) · launch speed u · time (transport) · show terms (toggle).
- **Equations shown (live):** (4.42) $\mathbf u=\mathbf U+\mathbf u'+\boldsymbol\Omega\times\mathbf x'$; (4.43)
  $\mathbf a=\frac{d\mathbf U}{dt}+\mathbf a'+2\boldsymbol\Omega\times\mathbf u'+\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'+\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')$;
  (4.45) with the bracket $\mathbf g-\frac{d\mathbf U}{dt}-2\boldsymbol\Omega\times\mathbf u'-\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')$;
  $\delta=\Omega ut^2$; $-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')=\Omega^2R\,\mathbf e_R$;
  $-2\boldsymbol\Omega\times\mathbf u=-2\Omega_zu_R\mathbf e_\varphi$.
- **Mirrors:** `core.rotating.frame_acceleration_terms`, `apparent_body_forces`, `coriolis_acceleration`,
  `effective_gravity`, `ch04.coriolis_projectile`, `core.rotating.projectile_paths` (§8).
- **Derivations:** D14 (the cone step draws e′₁ sweeping its cone), D15 ★★★ (12 steps; the two Ω × u′ steps each add
  one half of the Coriolis arrow; the triple-product step draws the centrifugal arrow), D16 (moving the terms right
  flips the arrows from accelerations to forces), D17 (sets the pole mode; the ½at² step overlays the parabola), D18
  (effective-gravity mode; the potential step draws equipotentials).
- **Depth features:** Explain tab (Ω and u → Coriolis magnitude 2Ωu → deflection and angle at time t → centrifugal
  Ω²R → effective gravity → which way the NH/SH turns → reading the current setting), synced Code tab, + linked views
  (3), transport (play/step/scrub, end-of-run card "the ball went straight; the table turned by Ωt"), modes (four),
  presets (Earth pole 1 h, merry-go-round, weightless parabola dU/dt = g, Southern Hemisphere), term bars, status
  ("↪ NH: deflected right" / "↩ SH: deflected left" / "⭕ Ω = 0: straight line").
- **Follows reference:** `angular_frequency_explorer_1.html` (linked views on one clock, modes, ghost reference,
  end-of-run summary).
- **Aha:** the ball flies straight; the floor turns under it — Coriolis is what that looks like from the floor, and its
  factor 2 is half from the turning axes and half from the moving position.

### E6 · which_bernoulli
- **A:** C05, C11, C12 (also shows N30 Ex. 4.2, N88 Lamb identity, N90 Lamb surfaces, N91 (4.72), N96 energy form
  (4.78), N97–N101 viscous irrotational forms, N102 the decision table, N103 pitot, N105 orifice)
- **Confusion removed:** "Bernoulli holds everywhere" — (4.19) holds along a streamline in steady inviscid flow; across
  streamlines only if the flow is irrotational (4.72); in unsteady flow only in the potential form (4.75); and (4.78) is
  an energy statement with its own hypotheses.
- **Why interactive:** toggle the flow's properties (steady? irrotational? viscous? compressible?) and watch the status
  pick the valid equation while B is plotted along and across streamlines — in the Rankine vortex B is flat along each
  circle and changes across them, in the cylinder flow it is flat everywhere, in the U-tube the ∂φ/∂t term carries the
  difference. The applicability map is a decision the reader must make, not read.
- **Stage:** (1) the scene with streamlines and two probes (pitot tube and static port; tank and orifice; Rankine
  vortex; cylinder; U-tube column; nozzle with hot gas); (2) B along the chosen streamline (flat?) and across
  streamlines (flat?) with the equation's terms as stacked bars (½u² teal, p/ρ orange, gz blue, ∂φ/∂t amber); (3)
  (hidePortrait) the decision table with the current row highlighted.
- **Controls:** scenario (six presets) · flow speed / head / swirl Γ · probe positions (drag) · time (U-tube) ·
  hypothesis toggles (steady, irrotational, inviscid, constant ρ).
- **Equations shown (live):** (4.19) $\tfrac12U^2+gz+p/\rho$ = const along a streamline; (4.68)
  $u_i\partial u_j/\partial x_i=-(\mathbf u\times\boldsymbol\omega)_j+\partial(\tfrac12u_i^2)/\partial x_j$; (4.69)
  $\partial\mathbf u/\partial t+\nabla B=\mathbf u\times\boldsymbol\omega$; (4.71), (4.72); (4.75)
  $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int\frac{dp}{\rho}+gz$ = const; (4.78)
  $h+\tfrac12\lvert\mathbf u\rvert^2+gz$ = const on streamlines; (4.82); pitot $\lvert\mathbf u\rvert=\sqrt{2\Delta p/\rho}$;
  Torricelli $u=\sqrt{2gh}$.
- **Mirrors:** `core.bernoulli.which_bernoulli`, `bernoulli_function`, `bernoulli_along_line`, `pitot_speed`,
  `torricelli_speed`, `unsteady_streamline_bernoulli`, `stagnation_temperature`, `core.bernoulli.bernoulli_scenario`
  (§8).
- **Derivations:** D06 (stream-tube element: the side-pressure step highlights the conical wall), D24 (Lamb identity
  steps draw u × ω on the Rankine vortex), D25 (the triple-product step shows ∇B ⟂ streamline), D26 (U-tube; the gauge
  step shows the book's sign doubling B).
- **Depth features:** Explain tab (the hypotheses of the scenario → the valid equation → B at each probe with numbers →
  the speed or pressure it predicts → where it would fail → reading the current setting), synced Code tab, + linked
  views (3), presets (six scenarios), status verdict ("✅ (4.71): constant along this streamline" / "✅ (4.72):
  constant everywhere" / "⚠️ rotational: B differs across streamlines" / "⏱ unsteady: use (4.75)"), term bars
  (stacked B), inspector (click a point: B's arithmetic), check questions built on the decision table.
- **Follows reference:** `overfitting_curves.html` (verdict line that changes with the regime) with
  `angular_frequency_explorer_1.html`'s highlighted "Right now" table.
- **Aha:** "Bernoulli" is four statements, not one — each is constant along something specific, and the flow's
  properties decide which one you may use.

### E7 · viscous_dissipation_heating
- **A:** C10 (also shows N74 (4.54) stress-work split, N76 (4.56) mechanical energy, N77 (4.58) ε ≥ 0, N79 (4.60) and
  the 7 = 7 ledger, N81 (4.63) entropy production)
- **Confusion removed:** "viscous forces just move energy around" — viscosity's *deformation work* ρε always turns
  mechanical energy into heat (never back), while its *force work* only moves kinetic energy from place to place.
- **Why interactive:** in a sheared channel the reader drags the wall speed U or the pressure gradient and watches, on
  one clock, the velocity profile, ε(y) (largest where the shear is), the temperature profile rising to a new steady
  state, and energy-budget bars that show shear work in = heat out; flipping the sign of μ (a "forbidden" preset)
  makes the entropy production negative and the status flags the second law.
- **Stage:** (1) the channel with u(y) (teal) and the dissipation ε(y) as a heat strip (rose); (2) T(y, t) relaxing to
  its steady parabola with the ghost of the exact steady solution; (3) (hidePortrait) budget bars: work of the moving
  wall, kinetic-energy change, dissipation ρε, heat conducted out; entropy production ≥ 0 readout.
- **Controls:** flow: Couette / Poiseuille / Couette + pressure gradient · U (m/s) or dp/dx · μ (log) · k (log) ·
  time (transport).
- **Equations shown (live):** (4.53); (4.54) stress work = deformation + force work; (4.56)
  $\rho\frac{D}{Dt}(\tfrac12u_j^2)=\rho g_ju_j-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}$;
  (4.57) $\frac{De}{Dt}=-p\frac{Dv}{Dt}+\frac1\rho\sigma_{ij}S_{ij}-\frac1\rho\frac{\partial q_i}{\partial x_i}$; (4.58)
  $\varepsilon=2\nu(S_{ij}-\tfrac13S_{mm}\delta_{ij})^2+\frac{\mu_v}\rho S_{mm}^2$; (4.60); (4.63).
- **Mirrors:** `core.constitutive.dissipation_rate`, `core.navier_stokes.internal_energy_terms`,
  `entropy_production`, `ch04.couette_heating` (§8).
- **Derivations:** D21 (dotting Cauchy with u lights the kinetic-energy bars), D22 (the subtraction step moves the
  deformation-work bar into the heat side), D23 (completing the square draws ε(y) as a square of the shear).
- **Depth features:** Explain tab (shear rate → ε = ν(du/dy)² with numbers → total dissipation per wall area → heat
  flux at each wall → peak temperature rise → entropy production → reading the current setting), synced Code tab, +
  linked views (3), term bars (energy budget), transport (heating up to steady state, end-of-run card "shear work in =
  heat out"), presets (water Couette, oil bearing, air, "negative μ" second-law violation), status ("🔥 heating: ε > 0" /
  "⛔ second law violated (μ < 0)").
- **Follows reference:** `forced_damped_vibrations.html` (transient → steady state on one time slider, boxed numbers,
  regime reading).
- **Aha:** viscosity is a one-way valve from motion to heat: ε is a sum of squares, so it can only be ≥ 0 — tiny in
  water, decisive in a hot bearing and in turbulence.

### E8 · boussinesq_buoyancy
- **A:** C13 (also shows R10 hydrostatic base state, N106 (4.84), N107 (4.85), N108 validity conditions, N110 the C_p
  mechanism, N112 viscous heating negligible, N113 (4.89), N114 the Boussinesq set and linear equation of state; N137
  reduced gravity and N² link)
- **Confusion removed:** "if density hardly changes, drop it everywhere" — the Boussinesq approximation drops it in the
  inertia but keeps it where it multiplies g, because g is huge compared with the accelerations of the flow.
- **Why interactive:** the reader sets the fluid, temperature contrast and depth and sees three small numbers (αδT,
  L/(c²/g), νU/(C_pδT L)) against their thresholds, while a warm blob rises and spreads under exactly the kept terms;
  switching the buoyancy term off stops the blob although ρ changed by only 0.2 % — the "small but multiplied by g"
  point needs the side-by-side.
- **Stage:** (1) a tank or atmosphere column with a warm (or cold) blob moving under buoyancy and diffusing (T field
  heatmap); (2) term bars per unit mass in the vertical momentum equation: inertia Dw/Dt, buoyancy ρ′g/ρ₀, pressure,
  viscous, and the dropped term ρ′Du/Dt/ρ₀ (grey, tiny); (3) (hidePortrait) the validity panel: three bars on log
  axes with their thresholds for the chosen scenario, table of scenarios with the current row lit.
- **Controls:** scenario: lake / ocean thermocline / lab tank / atmospheric boundary layer / deep atmosphere (fails
  L ≪ c²/g) · δT (K) · depth scale L (m, log) · keep buoyancy (toggle) · time (transport).
- **Equations shown (live):** (4.84) $\rho\frac{D\mathbf u}{Dt}=-\nabla p'+\rho'\mathbf g+\mu\nabla^2\mathbf u$; validity
  $\alpha\delta T\ll1$, $L\ll c^2/g$; (4.86) $\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u$;
  (4.88) $\rho C_p\frac{DT}{Dt}=\rho\varepsilon-\nabla\cdot\mathbf q$; (4.89) $\frac{DT}{Dt}=\kappa\nabla^2T$; linear
  equation of state $\rho=\rho_0[1-\alpha(T-T_0)]$.
- **Mirrors:** `ch04.boussinesq_validity`, `boussinesq_density`, `core.navier_stokes.buoyancy`,
  `boussinesq_momentum_terms`, `ch04.gaussian_blob_advection_diffusion`, `ch04.boussinesq_scenario` (§8).
- **Derivations:** D27 (the subtraction step draws the hydrostatic base state; the "why keep ρ′g" step shows the two
  grey vs coloured bars), D28 (the C_v → C_p step shows the pressure-work bar folding into the heating).
- **Depth features:** Explain tab (α and δT → αδT → g′ → c and c²/g → viscous-heating ratio → the blob's rise speed
  scale √(g′L) → verdict → reading the current setting), synced Code tab, + linked views (3), presets (five scenarios),
  status verdict ("✅ Boussinesq valid" / "⚠️ deep layer: L ≈ c²/g — use the anelastic or compressible equations"),
  term bars, transport (blob), a small table of real values with the current row highlighted.
- **Follows reference:** `angular_frequency_explorer_1.html` ("Right now" notes with a highlighted table, modes).
- **Aha:** density differences of a fraction of a percent drive the ocean and the atmosphere, because they are
  multiplied by g — Boussinesq keeps exactly that product and nothing else.

### E9 · dynamic_similarity_models
- **A:** C15 (also shows N130 two routes, R12 (4.99), N131 sphere C_D(Re), N133 (4.100), N134–N137 St, Re, Fr, Fr′/Ri,
  N139 C_p, N142 C_D, N147 M, N152 Pr, N156 Ex. 4.8 ship model; links to ch01 `buckingham_pi_machine`)
- **Confusion removed:** "a scale model just needs the right shape" — it needs the right dimensionless groups, and a
  model cannot usually match all of them (Froude and Reynolds conflict for ships; hence the friction correction).
- **Why interactive:** choose a prototype and a model scale; sliders for the model's speed and fluid show which groups
  match (badges) and the consequences — the sphere data points of many sizes and fluids collapse onto one C_D(Re)
  curve only when plotted against Re; the ship's drag split updates as the Froude-matched speed forces a Reynolds
  mismatch.
- **Stage:** (1) prototype and model side by side (ship hulls, spheres, a building in wind, a rotating tank vs the
  atmosphere) with their l, U, ν, g; (2) the groups as paired bars (prototype vs model: Re, Fr, St, M, Ri or Ro) with
  "matched" badges; (3) (hidePortrait) the sphere C_D(Re) curve with synthetic experiments (dimensional axes ⇄
  collapsed axes toggle) or the ship-drag waterfall (model total → friction + wave → scaled wave + prototype friction).
- **Controls:** mode: sphere / ship (Ex. 4.8) / stratified flow (Fr′, Ri) / rotating tank (Rossby number, forward
  pointer) · scale λ = l_m/l_p (log) · model speed U_m · model fluid (water, air, glycerine) · axes: dimensional /
  dimensionless.
- **Equations shown (live):** (4.100) scalings; (4.101)
  $[\frac{\Omega l}U]\frac{\partial\mathbf u^*}{\partial t^*}+(\mathbf u^*\cdot\nabla^*)\mathbf u^*=-\nabla^*p^*+[\frac{gl}{U^2}]\mathbf g^*+[\frac{\mu}{\rho Ul}]\nabla^{*2}\mathbf u^*$;
  (4.102)–(4.105) St, Re, Fr, Fr′, Ri; (4.107) $C_D=F_D/\tfrac12\rho U^2A$; (4.111) M = U/c; Ex. 4.8
  $U_m=U_p\sqrt{l_m/l_p}$; Ro = U/(Ωl) defined in the explainer from the ratio of the advective to the Coriolis term of
  (4.45) (taught in Ch. 13).
- **Mirrors:** `core.similarity.reynolds_number`, `froude_number`, `richardson_number`, `sphere_drag_coefficient`,
  `froude_scaled_speed`, `ch04.ship_drag_extrapolation`, `core.similarity.model_prototype` (§8).
- **Derivations:** D30 (each scaling step turns one dimensional term into a bracket; the final step lights the three
  group bars).
- **Depth features:** Explain tab (prototype groups → model groups with the chosen λ and fluid → which match → the
  speed Froude matching demands → the Re mismatch λ^{3/2} → the drag extrapolation with numbers → reading the current
  setting), synced Code tab, + linked views (3), modes (four), presets (1:25 ship in water, 1:10 car in air, sphere in
  glycerine vs air, lab tank vs atmosphere), status ("✅ dynamically similar" / "⚠️ Fr matched, Re off by ×125: correct
  friction"), inspector (click a data point: its d, U, ν and the Re arithmetic).
- **Follows reference:** `stride_padding_playground.html` (presets that are the classic settings, the formula with
  numbers plugged in, badges explaining the result).
- **Aha:** the equations only know the groups — match St, Re, Fr (or whichever the physics keeps) and the model is the
  prototype in miniature; when two groups cannot both match, you correct for the one you had to give up.

### B1 · kinematic_free_surface (backup)
- **A:** C14 (also shows N116 pillbox jump conditions, N117 no-slip, N118 (4.90), N119 (4.92), N120 (4.93))
- **Confusion removed:** "the free surface is where the fluid happens to be" — it is a material surface: the fluid's
  normal velocity equals the surface's, Dη/Dt = 0, and only that normal part is defined.
- **Why interactive:** a wave surface with particles riding on it; the reader moves a probe along the surface and sees
  u·n = u_s·n; switching to a "leaky" surface (a moving interface with mass flux, like evaporation or a shock) makes the
  residual (4.92) non-zero and particles cross.
- **Stage:** surface + particles + normal-velocity bars; linearised vs full condition residual vs ka (log–log, slope 2).
  **Controls:** amplitude ka, wavelength, time, mode (free surface / moving wall / leaky interface / pillbox).
  **Equations:** (4.90) $d\eta/dt=\partial\eta/\partial t+(\mathbf u_s\cdot\nabla)\eta=0$; (4.91) $D\eta/Dt=0$; (4.92)
  $(u_{rel})_n=(1/\lvert\nabla\eta\rvert)D\eta/Dt$; (4.93). **Mirrors:** `core.interfaces.kinematic_bc_residual`,
  `surface_normal_speed`, `relative_normal_velocity`, `ch04.linear_wave_surface` (§8).
- **Derivations:** D29. **Depth features:** explain, code, + linked views, transport, status, inspector. **Follows
  reference:** `forward_noising_lab.html` (click a point to see its exact arithmetic). **Aha:** a free surface is made of
  fluid particles, so it moves exactly as fast as the fluid pushes it — no faster, no slower.

## 6. Python animations and interactive figures (A ID → what, why, player/figure kind)
Animations (`animate` + `show_animation`, ≤ 120 frames, dpi 80, FAST halves the frames):
- **C01** — the expanding flow: a material interval stretching while a fixed box beside it loses mass, with the running
  budget bars; why: "same mass, different box" must be watched; `player="frames"`.
- **C04** — Ex. 4.3 bore with the CV riding along (h_in, h_out, the hydrostatic side forces as arrows); why: the moving
  CV turns an unsteady wave into a steady budget; `player="video"`.
- **C09** — the Coriolis projectile at the pole in the inertial (straight) and rotating (curved) frames side by side;
  why: the whole idea is two observers of one motion; `player="video"`.
- **C12** — a U-tube column oscillating with the pressure-difference bar ρL dU/dt tracking the acceleration; why: the
  ∂φ/∂t term is invisible in a snapshot; `player="video"`.
- **C13** — a warm Gaussian blob advected by a uniform current and spreading (σ² = σ₀² + 2κt), the exact solution of
  (4.89); why: advection and diffusion together; `player="video"`.
- **C14** — a linear surface wave with particles on the surface riding up and down (never leaving it); why: the
  kinematic condition is a statement about motion; `player="video"`.

Interactive figures (`slider_figure`, plotly, precomputed; 3-D plotly where noted):
- **C02** — slider over time: ∂ρ/∂t and ∂(ρu)/∂x along x for the expanding flow (mirror images summing to 0); why: the
  balance holds at every point and time.
- **C03** — slider over the source strength of source + stream: ψ contours at equal Δψ with the speed heatmap; why:
  spacing = speed as the body grows.
- **C04** — slider over the box height H: F_D/l from the four faces converging to ρ∫U(U∞ − U)dy; why: the far faces
  must be far enough.
- **C07** — slider over μ_v: normal stress on a plane vs angle for an expanding flow (p̄ vs p gap); why: bulk viscosity
  acts only on expansion.
- **C08** — slider over time: Stokes' first problem profile u(y, t) with the local and viscous terms; why: diffusion of
  momentum in time.
- **C10** — slider over U: Couette ε(y) and T(y) with the budget numbers; why: heating grows as U².
- **C11** — slider over the Rankine core radius: B along vs across circles; why: rotational vs irrotational regions.
- **C15** — slider over Re: sphere C_D(Re) (Morrison) with the Stokes line and the current point; plus a two-parameter
  collapse figure (Stokes layer at two (U, ν, l) on one dimensionless curve); why: similarity at a glance.
- **N15 (→ C03)** — plotly 3-D: two stream-surface families and their intersection streamlines (Fig. 4.1); why: a
  3-D statement.
- **N64 (→ C09)** — plotly 3-D: exaggerated-Ω Earth with g, Ω²R and g_e arrows at several latitudes and the
  equipotential; why: the bulge geometry.

## 7. From-scratch moments
Every book section with a computable A item has at least one (§4.1 has no A item; its ledger is a printed table):
- §4.2 **C01** — midpoint sums for the mass in a fixed and a material interval vs `core.conservation.mass_budget`;
  **C02** — hand-written central stencils for ∂ρ/∂t + ∂(ρu)/∂x vs `core.navier_stokes.continuity_residual`
  (`assert np.allclose`, residual ≈ 0 for the expanding flow, ≠ 0 for a wrong ρ).
- §4.3 **C03** — u, v from ψ by central differences vs `velocity_from_streamfunction_2d`; flux along a gate by
  `np.trapezoid` vs ψ₂ − ψ₁.
- §4.4 **C04** — the four face integrals of Ex. 4.1 by `np.trapezoid` vs `ch04.wake_drag_per_span`; **C06** — a
  first-index stress divergence by hand (non-symmetric τ) vs `cauchy_terms` (and the second-index version failing).
- §4.5 **C07** — τ_ij = −pδ_ij + 2μS_ij + λS_mmδ_ij by explicit loops over i, j vs `core.constitutive.newtonian_stress`.
- §4.6 **C08** — the Poiseuille residual of (4.39b) by hand (pressure + viscous) vs `ns_incompressible_terms`.
- §4.7 **C09** — the five acceleration terms with `np.cross` vs `core.rotating.frame_acceleration_terms` and vs a
  finite-difference acceleration of a path seen in both frames.
- §4.8 **C10** — ε by the explicit double sum σ_ijS_ij/ρ vs `core.constitutive.dissipation_rate` (and the
  sum-of-squares route).
- §4.9 **C11** — u × ω by `np.cross` and the Lamb identity checked pointwise vs `core.navier_stokes.lamb_vector`;
  **C13** — the heat-equation residual of the Gaussian blob by hand vs `temperature_equation_residual`.
- §4.10 **C14** — Dη/Dt by central differences for a moving wall and a linear wave vs `kinematic_bc_residual`.
- §4.11 **C15** — Re and Fr by hand for the ship model and a sympy substitution of (4.100) into (4.39b) vs
  `core.similarity.nondimensional_ns_coefficients`.

## 8. Notes for the implementer
Functions the A items' figures and the explainers need that analysis §4 did not plan (all scalar-callable for parity
rows, SI units, docstrings citing § and Eq.):
- `ch04.cv_scenario(name, **p) -> dict(faces, mass_flux, momentum_flux, body, surface, residual)` — E1's four
  scenarios (wake, bore, jet on a plate, rocket) with per-face numbers; `ch04.jet_plate_force(rho, V, A, theta)` (our
  extension, labelled; normal plate F = ρV²A).
- `core.streamfunction.streamfunction_preset(name, x, y, **p) -> psi` (uniform, stagnation, source + stream, cylinder,
  line vortex, shear, axisymmetric uniform) — E2 and the C03 figures; parity: velocities vs `ch03.cylinder_flow`.
- `core.constitutive.stress_on_plane(G, p, mu, mu_v, n) -> (normal, shear)` — E3's rotatable plane (wraps `traction`).
- `core.navier_stokes.exact_solution(name, x, t, **p) -> (u, p)` — Couette, plane Poiseuille, Stokes' first problem,
  Taylor–Green, Lamb–Oseen, cylinder potential flow — for E4 and the C08 figures (each already in the validation plan;
  one dispatcher).
- `core.rotating.projectile_paths(u0, Omega, t) -> (inertial_xy, rotating_xy)` — E5 and the C09 animation;
  `core.rotating.coriolis_parameter(lat)` = 2Ω sin φ (named forward pointer only).
- `core.bernoulli.bernoulli_scenario(name, **p) -> dict(B_along, B_across, valid_forms, …)` — E6's six scenarios.
- `ch04.couette_heating(y, U, h, mu, k, T0, dpdx=0.0) -> dict(u, eps, T, q_walls)` and the transient T(y, t) for E7
  (series solution or FTCS via `core.diffusion`).
- `ch04.boussinesq_scenario(name) -> dict(alpha, dT, L, U, nu, cp, c)` — E8's five scenarios feeding
  `boussinesq_validity`.
- `core.similarity.model_prototype(l_p, U_p, scale, fluid_p, fluid_m, match="Fr") -> dict(groups_p, groups_m, U_m,
  mismatch)` — E9; `core.similarity.rossby_number(U, Omega, l)` for E9's forward-pointer mode (docstring: defined from
  (4.45) scaling, taught in Ch. 13).
- `ch04.stokes_first_problem(y, t, U, nu)` (erfc profile) — the C08/C15 collapse figure and an E4 solution.
- `ch04.linear_wave_surface(x, z, t, a, k, g, H) -> (eta, u, w)` — the C14 animation and B1 (test field only; Ch. 7
  derives it).
- `ch04.u_tube_column(t, L, h0, g)` and `ch04.accelerating_sphere_pressure(theta, a, U, dUdt, rho)` — C12's worked
  numbers and animation.
- Parity rows must include the wake drag closed form (E1), Δψ = flux (E2), τ₁₂ = μγ (E3), the Poiseuille balance (E4),
  deflection → Ωut² (E5), (4.72) B uniform for the cylinder (E6), ρε = μ(U/h)² (E7), αδT and c²/g (E8), U_m = U_p√λ
  (E9). Book numbers stay in `tests/book_values_ch04.json`.
- Budget: exact-solution fields evaluated on ≤ 200×200 grids; the sphere "experiments" and Stokes-layer collapse
  precomputed once; sympy derivation checks (D09, D15, D30) cached with `functools.lru_cache`; animations use FAST = 30
  frames; the whole notebook < 5 min on Colab CPU.
- Machinery note: the 3-digit NOTE ids (N100–N156) are accepted by `nbkit.curation_items` (`\d{2,3}`); keep them
  zero-padded to three digits only when ≥ 100.
