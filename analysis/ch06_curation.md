# Chapter 6 — Ideal Flow: curation
(from `analysis/ch06.md` — 160 inventory rows, 110 equation labels (6.1)–(6.109) with (6.23a, b), 48 derivations in
§2b, 17 candidate A ideas in §3; `book.yaml` ch06, `policy.tier_a_full_treatment: [12, 18]`, `coverage: exhaustive`.
Curated 2026-09-23, concept-curator. Read: `knowledge/CUMULATIVE.md`, `concept_map.md`, `primers.md` (P01–P148),
`notation.md`, `viz_patterns.md`, `knowledge/ch05.md` §8, `analysis/ch05_curation.md` (format).)

Counts: A 15 · B 116 · C 29 · RECAP 32 · SKIP 7 · derivations written out 31 (★ 8 ★★ 18 ★★★ 5) · demoted to statements 14

Tier words (parsed by `tools/nbkit.py`): **CORE 15 · NOTE 106 · RECAP 32 · SKIP 7** = 160 rows; **DERIVATION 31**.
Depth: A 15 (all CORE) · B 116 (93 NOTE + 23 RECAP) · C 29 (13 NOTE: N03, N04, N15, N26, N37, N52, N57, N71, N76,
N86, N91, N105, N106 · 9 RECAP: R02, R04, R06, R07, R08, R09, R15, R22, R31 · 7 SKIP S01–S07). Every inventory row
number `[#n]` appears exactly once in §2.

**Decisions that shape this chapter.**
1. **Fifteen A items from the analyst's seventeen (9 % of the rows).** Two merges: (a) the conformal-mapping idea
   (6.63)–(6.64) and the Zhukhovsky map (6.65)–(6.69) are one A block, **C11**, because the angle-preservation argument
   is only worth its derivation when it immediately carries a circle flow onto an ellipse; (b) the analyst's "Bernoulli
   and C_p" (A4) is not an A item — Bernoulli everywhere is ch04's (4.72) (RECAP R05 inside **C03**), and C_p (6.32) is
   a B item of **C05** (the half-body is where the book defines it). The analyst's A16 is split the other way: the
   sphere (a result Ch. 8, 9 and §6.9 build on) goes with the Stokes stream function into **C13**, and the airship plus
   the axial singularity method (the chapter's *inverse* method) become **C14**. Images and Example 6.1 stay an A item
   (**C08**): the method is load-bearing (Ch. 7 bottom conditions, Ch. 13 walls, Ch. 14 ground effect) and Example 6.1
   is the chapter's only unsteady-Bernoulli pressure prediction.
2. **SEEN rows are RECAP** (32 rows, R01–R32, in inventory order), reminded where needed with the earlier chapter's
   function. Rows the analyst marks "SEEN + NEW" whose NEW part is a construction or a derivation (the vortex ψ form,
   φ with orthogonality, the Stokes flux 2πdψ, the sphere construction, the moving-sphere pressure and added mass) are
   **B NOTE items** inside their A block. **No CORE row is a SEEN row**: where the headline flow was already met (the
   cylinder (6.33), ch03 C05), the CORE row is its new consequence (d'Alembert, #49) and the flow itself is a B RECAP.
3. **Derivations written out: 31** (★ 8 · ★★ 18 · ★★★ 5), covering 34 analysis rows (three pairs merged: a-D21 + a-D22
   → D14, a-D37 + a-D38 → D25, a-D40 + a-D41 → D26). Rule (b) ("the book never writes it out and the lesson needs
   it"): D03 (∇² ln r = 2πδ, Exercises 6.1b/6.3b), D08 (half-body surface C_p and its zero, our reduction), D31 (added
   mass by kinetic energy, Exercise 6.49 — the independent second route). **14 demoted to statements** (§4c): the SEEN
   ψ identity, Bernoulli everywhere, the polar forms, the quadratic harmonic polynomials, the element velocities, the
   doublet's circles, the (6.41) streamline equation, the element complex potentials, Example 6.2's boundary values,
   the flux 2πdψ, the spherical table, the coordinate-free sphere potential, the steady limit (6.106) and the cylinder's
   added mass (Exercise 6.45).
4. **Book slips taught in our own words** (analysis §9, carried into §8 below): (6.61)'s 1/z² coefficient and extra
   outer square (N62, D18 — the correct coefficient is −(Ud/π + Γ²/4π²), and only the 1/z term matters); (6.104)'s
   middle sign (N100, D29 — the bracket's last term is +u_s/a³); (6.108)'s stray dφ (N104, D30); "(6.8)" for "(6.15)" in
   the source velocities (N24); "Figure 6.5" for "Figure 6.3" (N66); "(6.5) and (6.12), respectively" reversed and
   "(6.43)" for "(6.44)" (N43); "Section 3" for "Section 6.3" (N53); "first-order" central differences that are
   second-order accurate (R19, D22); the Example 6.2 FORTRAN loop index and Δψ units m² → m²/s (N74). (The earlier "(6.82) slip" is retracted: the printed form is r × App. B, correct.)
5. **Conventions stated where first used** (⚠️ callouts with numbers, §8): **Γ counterclockwise** in (6.6), (6.8),
   (6.47) and in all project code vs the **clockwise Γ** of (6.36)–(6.40), (6.52), (6.61)–(6.62), (6.68) and Example
   6.1 (the flow's circulation is −Γ, so L = +ρUΓ); the **2-D doublet vector** d = Σx_i m_i points from sink to source
   (cylinder d = −2πUa² e_x), (6.49)'s scalar d is a dipole −d e_x; the 3-D dipole −d e_z (sphere d = 2πa³U), moving
   sphere d(t) = +2πa³u_s; θ from +x (downstream) in polar formulas, from the upstream stagnation point in the Fig. 6.10
   comparison, from the sphere's velocity in (6.106); **z horizontal along the stream in §6.8**, w = the z-velocity in
   §6.9 (not the complex potential); the Stokes ψ in m³/s vs plane ψ in m²/s; forces D, L on the body vs F on the fluid.
6. **Climate hooks, where real:** ω_z = −∇²ψ (6.4) is the ψ–vorticity inversion of QG models (Ch. 13 PV inversion,
   C02); Gauss–Seidel/SOR relaxation is how bounded-domain streamfunctions and pressure Poisson problems are solved
   (Ch. 10, Ch. 13, C12); the corner/stagnation flow ψ = 2Axy (6.24) is the deformation field of frontogenesis (Ch. 13,
   C04); flow past an obstacle and images off a wall are the building blocks of topographic flow and coastal images
   (C06, C08); the Stokes stream function (6.75) is the zonal-mean overturning (Hadley/Ferrel) streamfunction's
   geometry (C13); added mass is why a rising thermal or bubble accelerates at most at 2g (C15).

## 1. Teaching order (A IDs grouped by book section, B/C IDs under each; one sentence each: "once you see X, Y follows")
A items in **bold**; B and C items listed where they are taught, with depth in brackets. Equations written beside their
numbers so downstream agents have them.

**§6.1 Relevance of irrotational constant-density flow theory**
- **C01 The ideal-flow equations and where they can be trusted**, $\nabla\cdot\mathbf u=0$ and
  $\rho\,D\mathbf u/Dt=-\nabla p$ *(Eq. 6.1)* [#1]: once you see that ω = 0 makes the net viscous force
  $-\mu\nabla\times\boldsymbol\omega$ vanish even though μ ≠ 0, you see why only the no-through-flow condition survives
  (N02 [B]), why the theory holds outside thin attached boundary layers at large Re (N01 [B]) and fails in wakes and
  separated regions (N03 [C]), and why Kelvin keeps the outer flow irrotational (R01 [B]); exclusions named (N04 [C]).

**§6.2 Two-dimensional stream function and velocity potential**
- **C02 ψ and φ: ω_z = −∇²ψ, two Laplace problems, vortices and sources as δ sources** [#9]: once
  $\omega_z=-\nabla^2\psi$ *(Eq. 6.4)* is written, irrotational flow is $\nabla^2\psi=0$ *(Eq. 6.5)* (N05 [B]) except at
  point vortices $\nabla^2\psi=-\Gamma\delta(x-x')\delta(y-y')$ *(Eq. 6.6)* (N06 [B]), and the same story for φ
  (continuity $\nabla^2\phi=q$ *(Eq. 6.11)*, N10 [B]; Laplace (6.12), N11 [B]; sources (6.13), N12 [B]); the elementary
  solutions follow: uniform (6.7)/(6.14) (N07, N13 [B]), vortex $\psi=-\frac{\Gamma}{2\pi}\ln r$ *(Eq. 6.8)* (N08 [B]),
  source $\phi=\frac{m}{2\pi}\ln r$ *(Eq. 6.15)* (N14 [B]); equipotentials cross streamlines at right angles (6.10) (N09
  [B]); polar forms (6.21)–(6.22) (N18, N19 [B]) with the recaps of ψ (6.3) (R03 [B]), continuity (6.2) (R02 [C]),
  irrotationality (6.9), (6.20) (R04, R07 [C]), polar continuity (6.19) (R06 [C]) and polar Laplacians (6.23a, b) (R08,
  R09 [C]); ψ and φ fuse into w later (N15 [C] → C09).
- **C03 Superposition and the no-through-flow condition: any streamline can be a wall** [#24]: once
  $\nabla^2(\phi_1+\phi_2)=\nabla^2\phi_1+\nabla^2\phi_2=0$ is seen, flows are built by adding elements and the body is
  whichever streamline the sum makes impermeable — $\partial\phi/\partial n=0$ or $\partial\psi/\partial s=0$ on the
  surface *(Eq. 6.16)* (N16 [B]), the far field $\partial\psi/\partial y=U$ *(Eq. 6.17)* (N17 [B]) — and pressure comes
  free from Bernoulli everywhere, $p+\tfrac12\rho\lvert\nabla\phi\rvert^2=$ const *(Eq. 6.18)* (R05 [B]).

**§6.3 Construction of elementary flows in two dimensions**
- **C04 The element kit and the doublet as a source–sink limit** [#40]: once a source and a sink are pushed together
  with 2mε held fixed, $\phi=-\mathbf d\cdot\mathbf x/2\pi r^2=\frac{\lvert\mathbf d\rvert}{2\pi}\frac{\cos\theta}{r}$
  *(Eq. 6.29)* appears — the element that makes closed bodies; the kit around it: harmonic polynomials (N20 [B]),
  corner flow $\psi=2Axy$ *(Eq. 6.24)* (R10 [B]) and its three relatives (6.25)–(6.27) (N21–N23 [B]), vortex and source
  velocities (R11 [B], N24 [B]), the pair (6.28) (N25 [B]); where each is singular (N26 [C]).
- **C05 The half-body: a stream plus a source** [#43]: once $\psi=Ur\sin\theta+\frac{m}{2\pi}\theta$ *(Eq. 6.31)* is
  contoured, the stagnation point $a=m/2\pi U$, the dividing streamline ψ = m/2 and the width $h_{\max}=m/2U$ follow;
  (6.30) (N27 [B]); the pressure coefficient $C_p=1-\lvert\mathbf u\rvert^2/U^2$ *(Eq. 6.32)* (N28 [B]) and its curve
  along the body (N29 [B]).
- **C06 The circular cylinder and d'Alembert's paradox** [#49]: once a doublet $\mathbf d=-2\pi Ua^2\mathbf e_x$ is
  added to the stream, $\psi=U(r-a^2/r)\sin\theta$ *(Eq. 6.33)* (R12 [B]) makes r = a a streamline, the velocity (6.34)
  (R13 [B]) gives $C_p=1-4\sin^2\theta$ *(Eq. 6.35)* (N30 [B]), and the fore–aft symmetric pressure gives D = 0; real
  flows separate (N31 [B], → Ch. 9); the moving cylinder is a doublet (R14 [B]).
- **C07 The cylinder with circulation: stagnation points move, lift L = ρUΓ** [#57]: once a clockwise vortex is added,
  $\psi=U(r-a^2/r)\sin\theta+\frac{\Gamma}{2\pi}\ln(r/a)$ *(Eq. 6.36)* (N32 [B]), the surface speed (6.37) (N33 [B]) moves
  the stagnation points to $\sin\theta=-\Gamma/4\pi aU$ *(Eq. 6.38)* (N34 [B]), the surface pressure (6.39) (N35 [B])
  integrated round the body (N36 [B]) leaves only $L=\rho U\Gamma$ *(Eq. 6.40)*; every Γ satisfies the same boundary
  conditions — uniqueness needs the Kutta condition (N38 [B], → Ch. 14); Magnus and history (N37 [C]).
- **C08 The method of images and Example 6.1's wall-pressure signal** [#65]: once the rule "vortex images flip sign,
  source images keep it" (N39 [B]) turns a wall into a symmetry line, a vortex beside a wall (R16 [B]) drifts at
  Γ/4πh (R17 [B]) and the unsteady Bernoulli equation gives the pressure at the wall,
  $\frac{p(0,0,t)-p_\infty}{\rho}=\frac{\Gamma^2}{4\pi^2}\frac{(\Gamma t/4\pi h)^2-h^2}{((\Gamma t/4\pi h)^2+h^2)^2}$;
  two sources / source by a wall / corner (6.41) (N40 [B]); circle images (R15 [C]).

**§6.4 Complex potential**
- **C09 The complex potential w = φ + iψ, Cauchy–Riemann and the complex velocity** [#66]: once $w\equiv\phi+i\psi$
  *(Eq. 6.42)* of $z=x+iy=re^{i\theta}$ *(Eq. 6.43)* (N41 [B]) is analytic, Cauchy–Riemann (6.44) (N42 [B]) makes φ- and
  ψ-lines orthogonal and $dw/dz=u-iv$ *(Eq. 6.45)* (N43 [B]) — every analytic function is a flow: corners
  $w=Az^n$ *(Eq. 6.46)* (N44 [B]) and every element in one line (6.47)–(6.53) (N45–N51 [B]); bridge to forces (N52 [C]).

**§6.5 Forces on a two-dimensional body**
- **C10 Blasius's theorem and Kutta–Zhukhovsky lift for any body** [#91]: once the pressure integral
  (6.54)–(6.56) (R18 [B], N53–N55 [B]) is written as the complex force (6.57) (N56 [B]) and Bernoulli is inserted
  (N57 [C], (6.58) N58 [B]) with the tangency trick (6.59) (N59 [B]), Blasius
  $D-iL=\frac{i\rho}{2}\oint_C(dw/dz)^2dz$ *(Eq. 6.60)* (N60 [B]) meets the far field of any body (N61 [B]), and the
  residue theorem (N63 [B]) keeps only the 1/z term of (6.61) (N62 [B]): $D=0$, $L=\rho U\Gamma$ *(Eq. 6.62)*.

**§6.6 Conformal mapping**
- **C11 Conformal mapping and the Zhukhovsky transformation** [#95]: once $\delta w=\frac{dw}{dz}\delta z$ *(Eq. 6.63)*
  (N64 [B]) shows that every small element is turned and stretched alike, angles are kept (6.64) (N65 [B]), a
  rectangular (φ, ψ) grid maps to a flow net (N66 [B]), and $z=\zeta+b^2/\zeta$ *(Eq. 6.65)* carries a circle of radius a
  (6.66) (N67 [B]) onto an ellipse (6.67) (N68 [B]) together with its flow (6.68) (N69 [B]) — provided the inverse
  $\zeta=\tfrac12z+\tfrac12(z^2-4b^2)^{1/2}$ *(Eq. 6.69)* takes the root outside the circle (N70 [B]).

**§6.7 Numerical solution techniques in two dimensions**
- **C12 Finite-difference Laplace: the average rule and Gauss–Seidel relaxation** [#104]: once the half-point
  differences (R19 [B]) give (6.70)–(6.71) (R20, R21 [B]), Laplace on a square grid says
  $\psi_{i,j}=\tfrac14[\psi_{i-1,j}+\psi_{i+1,j}+\psi_{i,j-1}+\psi_{i,j+1}]$ *(Eq. 6.72)*; a 16-point grid gives four
  equations (6.73) (N72 [B]), sweeping with the latest values converges (N73 [B]), and Example 6.2's contraction is
  solved (N74 [B]); other numerical routes named (N71 [C]); the boundary-only alternative (source panels) is N90 in C14.

**§6.8 Axisymmetric ideal flow**
- **C13 Axisymmetric potential flow: the Stokes stream function, 3-D elements and the sphere** [#126]: once
  $u_R=-\frac1R\frac{\partial\psi}{\partial z}$, $u_z=\frac1R\frac{\partial\psi}{\partial R}$ *(Eq. 6.75)* (R24 [B]) is
  inserted in ω_φ (6.76) (R25 [B]), the field equation (6.77) (N75 [B]) is *not* the Laplacian (no complex variables);
  flux 2πdψ (6.78) (N77 [B]), units (N76 [C]); φ-side (6.79)–(6.80) (R26 [B], N78 [B]); coordinate tables (6.81)–(6.85)
  (R27–R30 [B], N79 [B]); elements (6.86)–(6.88) (N80–N82 [B]); the sphere (6.89) (N83 [B]) with velocity
  $u_\theta=-U[1+\tfrac12(a/r)^3]\sin\theta$ *(Eq. 6.90)* and $C_p=1-\tfrac94\sin^2\theta$ *(Eq. 6.91)* (N84 [B]);
  coordinate-free form (6.92) (N85 [B]) for §6.9; recaps (6.74) (R23 [B]) and two stream functions (R22 [C]).
- **C14 Bodies of revolution from axial singularities: airship and the axial singularity method** [#133]: once a
  closed body needs zero net source strength (N86 [C]), a uniform line sink (6.93)–(6.94) (N87, N88 [B]) plus a point
  source makes the airship (6.95) (N89 [B]); turned round, unknown segment strengths with ψ = 0 at N body points give
  $\psi_m=-\sum_n\frac{k_n}{4\pi}(r^m_{n-1}-r^m_n)+\tfrac12UR_m^2=0$, an N × N linear system; the 2-D analogue on the
  body surface, source panels, is our labelled extension (N90 [B]).

**§6.9 Three-dimensional potential flow and apparent mass**
- **C15 The accelerating sphere: surface pressure and added mass** [#152]: once the sphere's dipole follows the motion
  (setup N92 [B]; (6.96)–(6.97) N93, N94 [B]) and the force is the surface pressure integral (6.98) (N95 [B]), unsteady
  Bernoulli (6.99) (R32 [B]) with the chain rule for ∂φ/∂t (6.100)–(6.102) (N96–N98 [B]) and the surface velocity
  (6.103)–(6.104) (N99, N100 [B]) gives the pressure (6.105) (N101 [B]); its steady part is (6.91) again, no drag (6.106)
  (N102 [B]); its acceleration part integrates (6.107) (N103 [B]) to $M=\tfrac{2\pi}{3}\rho a^3$ (6.108) (N104 [B]) and
  $\mathbf F_E=(m+\tfrac{2\pi}{3}\rho a^3)\,d\mathbf u_s/dt$ *(Eq. 6.109)*; framing and why (N91, N105 [C]); w
  notation (R31 [C]).

**§6.10 Concluding remarks** — covered by N106 [C] (history, potential theory elsewhere, pointers to Ch. 9 and Ch. 14)
at the end of C15's block, and the exercise pointers S01–S07.

## 2. Chapter map (depth) — every inventory row exactly once
One row per inventory row (160, `[#n]` = analysis §2 row). `A parent` names the A block a B or C item is written in.
IDs: CORE C01–C15 (teaching order), NOTE N01–N106 and RECAP R01–R32 in inventory order, SKIP S01–S07 (exercises).

| ID | Item | § | Depth | Tier | A parent | Reason (A) / treatment (B) / pointer (C, SKIP) / source chapter (RECAP) |
|---|---|---|---|---|---|---|
| C01 | Eq. (6.1): the ideal-flow equations $\nabla\cdot\mathbf u=0$ and $\rho\,D\mathbf u/Dt=-\nabla p$ (p from its local hydrostatic value); they hold even if μ ≠ 0 because the net viscous force $-\mu\nabla\times\boldsymbol\omega$ vanishes when ω = 0 [#1] | 6.1 | A | CORE | – | load-bearing: the premise of the whole chapter, of linear water waves (Ch. 7), of the outer flow that drives every boundary layer (Ch. 9) and of lift theory (Ch. 14); D01 |
| N01 | Applicability map: $\mathrm{Re}=\rho UL/\mu\gg1$ confines vorticity to thin boundary layers; ideal flow predicts outer velocity, normal pressure forces, minimum-drag shapes, unsteady inertia — not skin friction, dissipation, duct flow, wakes, turbulence [#2] | 6.1 | B | NOTE | C01 | stated as a two-column table (predicts / does not) with one number: car at 10 m/s, L = 1 m in air → Re ≈ 6.7×10⁵, boundary layer ~ L/√Re ≈ 1.2 mm; `ch06.ideal_flow_applicability(Re, M, baroclinic, region)` verdicts; ch04 C15 recap of Re |
| N02 | Boundary conditions: (6.1) has only first derivatives, so only $\mathbf n\cdot\mathbf u=\mathbf n\cdot\mathbf U_s$ is kept and no-slip $\mathbf t\cdot\mathbf u=\mathbf t\cdot\mathbf U_s$ is dropped [#3] | 6.1 | B | NOTE | C01 | stated as the last step of D01 (why "μ ≠ 0 is allowed" silently needs no-slip dropped); the order-of-PDE ↔ number-of-conditions reminder (ch04 N115) |
| N03 | Real vs ideal at high Re: the boundary layer thins as Re grows, but when it separates the limit μ → 0 of real flow is not the flow with μ = 0 (Figs. 6.1–6.2) [#4] | 6.1 | C | NOTE | C01 | named with a pointer to Ch. 9 (separation) and to N31, where Fig. 6.10's rear pressure shows it |
| R01 | Why the hypotheses hold: $M=U/c\ll1$ for constant ρ, no baroclinic density field, and Kelvin $D\Gamma/Dt=0$ (5.8) keeps fluid irrotational unless it enters a boundary layer, wake or separated region [#5] | 6.1 | B | RECAP | C01 | ch05 C03 (Kelvin, its four restrictions) and ch04 C15 (M): reminded with `ch05.kelvin_hypotheses_text` inside `ch06.ideal_flow_applicability` (baroclinic → "Kelvin fails") |
| N04 | Summary of exclusions (inhomogeneous fluid, high Mach number, boundary layers, wakes, interior flows, rotational regions) and inclusions (flight, water waves, vehicles) [#6] | 6.1 | C | NOTE | C01 | named in one sentence closing C01; pointers Ch. 7 (waves), Ch. 14 (flight), Ch. 15 (compressible) |
| R02 | Eq. (6.2): 2-D continuity $\partial u/\partial x+\partial v/\partial y=0$ [#7] | 6.2 | C | RECAP | C02 | ch04 C02 (4.10): one sentence opening C02 |
| R03 | Eq. (6.3): stream function $u\equiv\partial\psi/\partial y$, $v\equiv-\partial\psi/\partial x$; ψ = const ⇒ $(dy/dx)_{\psi}=v/u$ [#8] | 6.2 | B | RECAP | C02 | ch04 C03 (D04) and ch03 C03: reminded with `core.streamfunction.velocity_from_streamfunction_2d` on ψ = Uy; a-D02 given (§4c) |
| C02 | Eq. (6.4): vorticity of a ψ-flow is minus its Laplacian, $\omega_z=\partial v/\partial x-\partial u/\partial y=-\nabla^2\psi$ — irrotational flow is then a Laplace problem for ψ and for φ, with point vortices and sources as δ sources [#9] | 6.2 | A | CORE | – | load-bearing: every §6.3–6.7 flow is a solution of these Laplace problems; the ψ–ω inversion is Ch. 10's vorticity–stream-function method and Ch. 13's PV inversion; D02, D03, D04 |
| N05 | Eq. (6.5): irrotational ⇒ $\nabla^2\psi=0$ [#10] | 6.2 | B | NOTE | C02 | stated as D02's result; `core.potential.laplacian_residual` = 0 on every element (printed table) |
| N06 | Eq. (6.6): a point vortex of strength Γ at x′ is a δ source of ψ, $\nabla^2\psi=-\Gamma\,\delta(x-x')\,\delta(y-y')$ [#11] | 6.2 | B | NOTE | C02 | stated with D03 (the book never shows why ln r carries a δ): `ch06.delta_flux_check` gives ∮∇ψ·n ds = −Γ on circles of radius 0.01…100 m |
| N07 | Eq. (6.7): uniform flow $\psi=-Vx+Uy$ [#12] | 6.2 | B | NOTE | C02 | stated with one number (U = 2, V = 1 m/s → ψ(1, 1) = 1 m²/s); `core.potential.Uniform` |
| N08 | Eq. (6.8): ψ of an ideal (counterclockwise) vortex, $\psi=-\frac{\Gamma}{2\pi}\ln\sqrt{(x-x')^2+(y-y')^2}$ [#13] | 6.2 | B | NOTE | C02 | stated with u_θ = Γ/2πr (ch03 C13, ch05 R03 reminder); ⚠️ Γ counterclockwise here and in all project code (C07 uses the book's clockwise Γ); `core.potential.Vortex` |
| R04 | Eq. (6.9): 2-D irrotationality $\partial v/\partial x-\partial u/\partial y=0$ [#14] | 6.2 | C | RECAP | C02 | ch03 C13 (3.23): one sentence before (6.10) |
| N09 | Eq. (6.10): velocity potential $u\equiv\partial\phi/\partial x$, $v\equiv\partial\phi/\partial y$; equipotentials $(dy/dx)_\phi=-u/v$ cross streamlines at right angles [#15] | 6.2 | B | NOTE | C02 | stated with D04 (slopes multiply to −1; ∇φ·∇ψ = 0, ∣∇φ∣ = ∣∇ψ∣) and ch03 N29 (φ exists only where ω = 0, simply connected); `ch06.orthogonality_check`; the flow-net figure |
| N10 | Eq. (6.11): continuity in φ is Poisson, $\nabla^2\phi=q(x,y)$ with a source density q [#16] | 6.2 | B | NOTE | C02 | stated as the φ-twin of (6.4); ch05 P139 (Poisson) reminder |
| N11 | Eq. (6.12): $\nabla^2\phi=0$ [#17] | 6.2 | B | NOTE | C02 | stated; same residual table as N05 |
| N12 | Eq. (6.13): a point source of strength m, $\nabla^2\phi=m\,\delta(x-x')\,\delta(y-y')$ [#18] | 6.2 | B | NOTE | C02 | stated with D03 (flux of ∇φ through any circle = m, `ch06.delta_flux_check`) |
| N13 | Eq. (6.14): uniform flow $\phi=Ux+Vy$ [#19] | 6.2 | B | NOTE | C02 | stated beside N07 (same flow, the other function) |
| N14 | Eq. (6.15): source potential $\phi=\frac{m}{2\pi}\ln\sqrt{(x-x')^2+(y-y')^2}$, m = volume flow per unit depth [m²/s], m < 0 a sink [#20] | 6.2 | B | NOTE | C02 | stated with one number (m = 2π m²/s → u_r = 1/r m/s, 1 m/s at r = 1 m); `core.potential.Source` |
| N15 | ψ and φ each describe a 2-D ideal flow and fuse into one complex potential; ψ extends to rotational flow, φ to unsteady and 3-D flow [#21] | 6.2 | C | NOTE | C02 | named with a pointer to C09 (w = φ + iψ) and C13/C15 (φ in 3-D) |
| N16 | Eq. (6.16): no flow through a solid surface $\mathbf n\cdot\mathbf U_s=(\mathbf n\cdot\mathbf u)$; stationary body ⇒ $\partial\phi/\partial n=0$ or $\partial\psi/\partial s=0$: a wall is a streamline and a streamline may be replaced by a wall [#22] | 6.2 | B | NOTE | C03 | stated as D05's result; `core.potential.normal_velocity_on` = 0 on r = a for the cylinder and on ψ = m/2 for the half-body |
| N17 | Eq. (6.17): far-field condition $\partial\phi/\partial x=U$ or $\partial\psi/\partial y=U$ (U = 0 quiescent) [#23] | 6.2 | B | NOTE | C03 | stated with `core.potential.far_field_check` (error ∝ 1/R² for a closed body, 1/R for the half-body — a log–log slope) |
| C03 | Superposition: Laplace is linear, $\nabla^2(\phi_1+\phi_2)=\nabla^2\phi_1+\nabla^2\phi_2=0$, so sums of elements are flows — the boundary conditions belong to the sum, and the body is the streamline the sum makes impermeable [#24] | 6.2 | A | CORE | – | load-bearing: the chapter's construction method (every body in §§6.3–6.9, the numerical methods of §§6.7–6.8), the wave superpositions of Ch. 7 and the panel methods of Ch. 14; D05 |
| R05 | Eq. (6.18): Bernoulli everywhere in steady irrotational flow, $p+\tfrac12\rho(u^2+v^2)=p+\tfrac12\rho\lvert\nabla\phi\rvert^2=p+\tfrac12\rho\lvert\nabla\psi\rvert^2=$ const; unsteady flow adds ρ∂φ/∂t [#25] | 6.2 | B | RECAP | C03 | ch04 C11 (4.72) and C12 (4.75): reminded with `core.bernoulli.bernoulli_function` (constant across streamlines for the half-body) and `core.potential.Flow.pressure`; a-D07 given (§4c) |
| R06 | Eq. (6.19): polar continuity $\frac1r\frac{\partial}{\partial r}(ru_r)+\frac1r\frac{\partial u_\theta}{\partial\theta}=0$ [#26] | 6.2 | C | RECAP | C02 | ch04 `core.curvilinear` (App. B): one line in the polar-forms note |
| R07 | Eq. (6.20): polar irrotationality $\frac1r\frac{\partial}{\partial r}(ru_\theta)-\frac1r\frac{\partial u_r}{\partial\theta}=0$ [#27] | 6.2 | C | RECAP | C02 | ch03 C13 (3.23): one line |
| N18 | Eq. (6.21): $u_r=\partial\phi/\partial r=\frac1r\partial\psi/\partial\theta$ [#28] | 6.2 | B | NOTE | C02 | stated (a-D08 given, §4c) with a sympy check through the chain rule; `core.potential.polar_velocity` |
| N19 | Eq. (6.22): $u_\theta=\frac1r\partial\phi/\partial\theta=-\partial\psi/\partial r$ [#29] | 6.2 | B | NOTE | C02 | stated beside N18; one number: vortex Γ = 2π → u_θ = −∂ψ/∂r = 1/r |
| R08 | Eq. (6.23a): $\nabla^2\psi=\frac1r\frac{\partial}{\partial r}(r\frac{\partial\psi}{\partial r})+\frac1{r^2}\frac{\partial^2\psi}{\partial\theta^2}=0$ [#30] | 6.2 | C | RECAP | C02 | ch04 `core.curvilinear` polar Laplacian: one line; used in D03 |
| R09 | Eq. (6.23b): $\nabla^2\phi=\frac1r\frac{\partial}{\partial r}(r\frac{\partial\phi}{\partial r})+\frac1{r^2}\frac{\partial^2\phi}{\partial\theta^2}=0$ [#31] | 6.2 | C | RECAP | C02 | as R08 |
| N20 | Harmonic polynomials: constants, linear (6.7)/(6.14), two quadratic families; higher powers ⇒ sharper corners, fractional ⇒ wider [#32] | 6.3 | B | NOTE | C04 | stated with `ch06.harmonic_polynomials(degree)` (sympy Re/Im of zⁿ, two per degree; a-D09 given, §4c); forward pointer to (6.46) in C09 |
| R10 | Eq. (6.24): stagnation / 90° corner flow $\psi=2Axy$, u = 2Ax, v = −2Ay, streamlines xy = const (Fig. 6.3) [#33] | 6.3 | B | RECAP | C04 | ch03 irrotational strain preset (`core.kinematics.velocity_gradient_preset`) and ch04 E2: reminded with `core.potential.Corner(A, n=2)`; ⚠️ Ch. 13's frontogenetic deformation field |
| N21 | Eq. (6.25): $\phi=2Axy$ — the same flow turned by 45° [#34] | 6.3 | B | NOTE | C04 | stated with the rotation x² − y² = 2x′y′ (ch02 C02 reminder) and an assertion on `Corner` with complex A |
| N22 | Eq. (6.26): $\psi=A(x^2-y^2)$ [#35] | 6.3 | B | NOTE | C04 | stated as the rotated family; `Corner` with complex coefficient |
| N23 | Eq. (6.27): $\phi=A(x^2-y^2)$ — equals (6.24) [#36] | 6.3 | B | NOTE | C04 | stated: both are Re/Im of Az² (preview of C09) |
| R11 | Vortex velocities from (6.8): $u=-\frac{\Gamma}{2\pi}\frac{y}{x^2+y^2}$, $v=\frac{\Gamma}{2\pi}\frac{x}{x^2+y^2}$, so u_θ = Γ/2πr (Fig. 6.4) [#37] | 6.3 | B | RECAP | C04 | ch03 C13 and ch05 R03 (5.2): reminded with `core.potential.Vortex` vs `core.biot_savart.point_vortex_velocity` parity |
| N24 | Source velocities $u=\frac{m}{2\pi}\frac{x}{x^2+y^2}$, $v=\frac{m}{2\pi}\frac{y}{x^2+y^2}$ ⇒ u_r = m/2πr; ∇·u = 0 except at r = 0 (Fig. 6.5) [#38] | 6.3 | B | NOTE | C04 | stated (a-D10 given, §4c); ⚠️ the book says "differentiation of (6.8)" — it is (6.15); total outflow m through any circle (links N12) |
| N25 | Eq. (6.28): source +m at (−ε, 0), sink −m at (+ε, 0), dipole vector $\mathbf d=\sum\mathbf x_im_i=-2m\varepsilon\mathbf e_x$ (from sink to source) [#39] | 6.3 | B | NOTE | C04 | stated as the start of D06; `ch06.source_sink_pair` |
| C04 | Eq. (6.29): the doublet as the limit ε → 0, m → ∞ with 2mε fixed, $\phi=\frac{m\varepsilon}{\pi}\frac{x}{r^2}=-\frac{\mathbf d\cdot\mathbf x}{2\pi r^2}=\frac{\lvert\mathbf d\rvert}{2\pi}\frac{\cos\theta}{r}$ (Fig. 6.6), the last member of the element kit (uniform, source, vortex, corner, doublet) [#40] | 6.3 | A | CORE | – | load-bearing: the doublet closes the cylinder (C06), the sphere (C13) and the moving sphere (C15); its far field 1/r is the signature of every closed body in the Laurent series (C10); the element kit is reused by Ch. 7 and Ch. 14; D06 |
| N26 | (6.7), (6.14), (6.24)–(6.27) are harmonic everywhere; (6.8), (6.15), (6.29) are singular at their centre [#41] | 6.3 | C | NOTE | C04 | named with the residual table of N05 (NaN at the centres by design) |
| N27 | Eq. (6.30): half-body potential $\phi=Ux+\frac{m}{2\pi}\ln\sqrt{x^2+y^2}=Ur\cos\theta+\frac{m}{2\pi}\ln r$ [#42] | 6.3 | B | NOTE | C05 | stated as the superposition; `core.potential.half_body(U, m)` |
| C05 | Eq. (6.31): the half-body $\psi=Uy+\frac{m}{2\pi}\tan^{-1}(y/x)=Ur\sin\theta+\frac{m}{2\pi}\theta$; stagnation point at x = −a = −m/2πU, dividing streamline ψ = m/2, half-width h = m(π − θ)/2πU → h_max = m/2U (Fig. 6.7) [#43] | 6.3 | A | CORE | – | load-bearing: the first body built by superposition — stagnation point, dividing streamline, C_p along a surface — the pattern every later body repeats; the blunt leading edge of Ch. 9; D07, D08 |
| N28 | Eq. (6.32): $C_p=\frac{p-p_\infty}{\frac12\rho U^2}=1-\frac{\lvert\mathbf u\rvert^2}{U^2}$ from Bernoulli with const = p∞ + ½ρU² [#44] | 6.3 | B | NOTE | C05 | stated with the ch04 C15 C_p recap; C_p = 1 at every stagnation point; `core.potential.pressure_coefficient` |
| N29 | Half-body surface C_p (Fig. 6.8): +1 at the nose, negative beyond, zero where the surface speed is U; $C_p=-(2k\cos\theta+k^2)$, $k=\sin\theta/(\pi-\theta)$; zero net force (Exercise 6.13) [#45] | 6.3 | B | NOTE | C05 | stated with D08 (our reduction, rule (b)); our root ≈ 113.2° from tan θ = −2(π − θ) (`ch06.half_body_cp_zero_angle`); the net-force integral → 0 as the body length grows (printed) |
| R12 | Eq. (6.33): cylinder = stream + doublet $\mathbf d=-2\pi Ua^2\mathbf e_x$, $\phi=U(r+a^2/r)\cos\theta$, $\psi=U(r-a^2/r)\sin\theta$; closed because the net source is zero (Figs. 6.9, 3.2a) [#46] | 6.3 | B | RECAP | C06 | ch03 C05 (`ch03.cylinder_flow`, Fig. 3.2): reminded, now *built* in D09 (why this d) with `core.potential.cylinder` ↔ `ch03.cylinder_flow` parity |
| R13 | Eq. (6.34): $u_r=U(1-a^2/r^2)\cos\theta$, $u_\theta=-U(1+a^2/r^2)\sin\theta$ [#47] | 6.3 | B | RECAP | C06 | ch03 C05: reminded as D09 step; `core.potential.cylinder(...).velocity_polar` |
| N30 | Eq. (6.35): surface pressure coefficient $C_p(r=a,\theta)=1-4\sin^2\theta$; stagnation points (a, 0), (a, π), minimum −3 at θ = ±π/2 [#48] | 6.3 | B | NOTE | C06 | stated as D09's result; number: air, U = 10 m/s → ½ρU² = 60 Pa, suction −180 Pa at the shoulders; `ch06.cylinder_surface_cp` |
| C06 | d'Alembert's paradox: the cylinder's pressure is fore–aft and top–bottom symmetric, so D = 0 (and L = 0); in general a steadily moving body feels no drag in 2-D ideal flow; real drag comes from skin friction and separation [#49] | 6.3 | A | CORE | – | load-bearing: the chapter's central surprise, proved for any body in C10 and in 3-D in C15; the reason Ch. 9 (boundary layers, separation) exists; the cylinder is the circle plane of Ch. 14's airfoils; D09 |
| N31 | Fig. 6.10: ideal C_p = 1 − 4 sin²θ vs a measured high-Re curve whose rear pressure stays low (separation) [#50] | 6.3 | B | NOTE | C06 | stated with our ideal curve and a labelled *qualitative* band (no cited dataset yet, analysis §8); ⚠️ Fig. 6.10's angle is measured from the upstream stagnation point (= π − θ); pointer Ch. 9 |
| R14 | A cylinder moving through still fluid = cylinder flow + uniform −U: the instantaneous streamlines of a doublet, $w=Ua^2/z$ (Figs. 6.11, 3.3b) [#51] | 6.3 | B | RECAP | C06 | ch03 C05 (Galilean frames) and E3 `galilean_frames_cylinder` (linked, not duplicated): reminded with `ch03.cylinder_flow(frame="fluid")` |
| N32 | Eq. (6.36): cylinder with a **clockwise** point vortex (the flow's circulation is −Γ), $\psi=U(r-a^2/r)\sin\theta+\frac{\Gamma}{2\pi}\ln(r/a)$, $u_\theta=-U(1+a^2/r^2)\sin\theta-\frac{\Gamma}{2\pi r}$ (Fig. 6.12) [#52] | 6.3 | B | NOTE | C07 | stated as D10's start; ⚠️ convention callout with numbers (book Γ clockwise vs project Γ_ccw = −Γ_cw); `core.potential.cylinder(U, a, Gamma_cw=…)` |
| N33 | Eq. (6.37): surface speed $u_\theta(r=a,\theta)=-2U\sin\theta-\Gamma/2\pi a$ [#53] | 6.3 | B | NOTE | C07 | stated as D10 step; faster on top (θ ≈ π/2 has ∣u_θ∣ = 2U + Γ/2πa) |
| N34 | Eq. (6.38): stagnation points $\sin\theta=-\Gamma/4\pi aU$ (two for Γ < 4πaU, one at Γ = 4πaU, off-body at $r=[\Gamma\pm\sqrt{\Gamma^2-(4\pi aU)^2}]/4\pi U$ for Γ > 4πaU) [#54] | 6.3 | B | NOTE | C07 | stated as D10's result; number: U = 10 m/s, a = 0.1 m, Γ = 2 m²/s → sin θ = −0.159, θ = −9.16° and −170.84°; critical Γ = 4πaU = 12.57 m²/s; `ch06.cylinder_stagnation_points` vs `Flow.stagnation_points` (Newton) |
| N35 | Eq. (6.39): surface pressure $p(r=a,\theta)=p_\infty+\tfrac12\rho[U^2-(-2U\sin\theta-\frac{\Gamma}{2\pi a})^2]$ [#55] | 6.3 | B | NOTE | C07 | stated as D11's start; `ch06.cylinder_surface_pressure` |
| N36 | Lift integral $L=-\int_0^{2\pi}p(r=a,\theta)\sin\theta\,a\,d\theta$ with n = e_r, dl = a dθ (Fig. 6.13) [#56] | 6.3 | B | NOTE | C07 | stated as D11 steps (why the minus sign); `ch06.surface_pressure_force` (periodic trapezoid) |
| C07 | Eq. (6.40): lift on the cylinder with circulation, $L=\rho U\Gamma$ — the circulation speeds the top and slows the bottom, the stagnation points slide down (6.38), and only the Γ-term of the pressure survives the integral [#57] | 6.3 | A | CORE | – | load-bearing: the first lift result, generalised to any body in C10 and to airfoils in Ch. 14 (Kutta condition); the non-uniqueness it exposes is why Ch. 14 needs the Kutta condition; D10, D11 |
| N37 | Kutta–Zhukhovsky history; circulation is set by viscosity at a sharp edge but its size by U, shape and angle; rotating cylinders, the Magnus effect, delayed separation on spinning balls [#58] | 6.3 | C | NOTE | C07 | named with pointers to Ch. 9 §9.9 (spinning balls) and Ch. 14 (Kutta condition) |
| N38 | Uniqueness and topology: unique in a singly connected domain; round a body every Γ in (6.36) meets the same boundary conditions — the Kutta condition (Ch. 14) picks one [#59] | 6.3 | B | NOTE | C07 | stated with `ch06.circulation_family_check` (u·n = 0 on r = a and far field U for Γ ∈ {0, 2πaU, 4πaU, 8πaU}; loop circulation −Γ_cw on every circle); ch03 N29 simply connected reminder |
| N39 | Method of images (Figs. 6.14–6.15): wall y = 0 — vortex image with the **opposite** sign, $\psi_2=\psi_1(x,y)-\psi_1(x,-y)$; source image with the **same** sign, $\phi_2=\phi_1(x,y)+\phi_1(x,-y)$ [#60] | 6.3 | B | NOTE | C08 | stated as D12's result (vortex half recalled from ch05 C13, D22; the source half new); `core.potential.mirror` |
| N40 | Eq. (6.41): source near a wall = two equal sources at x = ±a (Fig. 6.16); streamlines $x^2-y^2-2xy\cot(2\pi\psi/m)=a^2$; three readings (two sources, source by a wall, slit flow into a corner) [#61] | 6.3 | B | NOTE | C08 | stated (a-D19 given, §4c: the tangent addition formula in one line); `ch06.two_sources`, `ch06.two_source_streamline` curve drawn on the ψ contours |
| R15 | Images for circles (more than one image) and moving images in unsteady flow [#62] | 6.3 | C | RECAP | C08 | ch05 C13 (`circle_image_system`): one sentence; Milne-Thomson circle theorem named (our addition, `core.potential.circle_theorem`) |
| R16 | Example 6.1 setup: vortex of strength −Γ at (h, 0) beside the wall x = 0, ψ = original + opposite image [#63] | 6.3 | B | RECAP | C08 | ch05 C13 (D22, `wall_image_system`): reminded with `core.biot_savart.point_vortex_evolve(boundary="wall")` |
| R17 | Example 6.1 trajectory: self-induced velocity zero, $d\xi_y/dt=\Gamma/4\pi\xi_x$ ⇒ ξ(t) = (h, Γt/4πh) [#64] | 6.3 | B | RECAP | C08 | ch05 D22 (drift Γ/4πh): reminded as D13 steps 1–3; number Γ = 1 m²/s, h = 1 m → 0.0796 m/s |
| C08 | Example 6.1: the pressure at the wall under a passing vortex, from its image and unsteady Bernoulli, $\frac{p(0,0,t)-p_\infty}{\rho}=\frac{\Gamma^2}{4\pi^2}\frac{(\Gamma t/4\pi h)^2-h^2}{((\Gamma t/4\pi h)^2+h^2)^2}$ — suction first, over-pressure after [#65] | 6.3 | A | CORE | – | load-bearing: the method of images (walls as symmetry lines; Ch. 7 bottoms, Ch. 13 coasts, Ch. 14 ground effect) and the chapter's only unsteady-Bernoulli prediction (wall-pressure signals, Ch. 12); D12, D13 |
| C09 | Eq. (6.42): the complex potential $w\equiv\phi+i\psi$ — analytic in z, so Cauchy–Riemann (6.44) holds, φ- and ψ-lines are orthogonal and $dw/dz=u-iv$ (6.45): every analytic function of z is a 2-D ideal flow [#66] | 6.4 | A | CORE | – | load-bearing: the language of §§6.4–6.6 (forces by residues, conformal maps) and of Ch. 14's airfoil theory; the conjugate u − iv is the chapter's most common trap; D14, D15 |
| N41 | Eq. (6.43): $z\equiv x+iy=re^{i\theta}$ (Fig. 6.17) [#67] | 6.4 | B | NOTE | C09 | stated with the complex-plane primer (modulus, argument, `np.angle`, `np.abs`) |
| N42 | Eq. (6.44): Cauchy–Riemann $\partial\phi/\partial x=\partial\psi/\partial y$, $\partial\phi/\partial y=-\partial\psi/\partial x$; lines orthogonal except where w or dw/dz is 0 or ∞ [#68] | 6.4 | B | NOTE | C09 | stated as D14's middle result; `ch06.cauchy_riemann_residual` ≈ 0 for every element, ≠ 0 for z* (control) |
| N43 | Eq. (6.45): complex velocity $dw/dz=u-iv$; with (6.44) gives (6.2), (6.9) and Laplace for φ and ψ [#69] | 6.4 | B | NOTE | C09 | stated as D14's result; number: w = z² at z = 1 + i → dw/dz = 2 + 2i ⇒ u = 2, v = −2 m/s; ⚠️ the book writes "(6.5) and (6.12), respectively" for φ and ψ — reversed — and "(6.43)" where (6.44) is meant |
| N44 | Eq. (6.46): corner flow $w=Az^n=Ar^n(\cos n\theta+i\sin n\theta)$, α = π/n, $dw/dz=(A\pi/\alpha)z^{(\pi-\alpha)/\alpha}$: stagnation at the corner for α < π, infinite speed for α > π; n = ½ is the flat plate [#70] | 6.4 | B | NOTE | C09 | stated as D15's result (written out as part of C09, a-D23); `core.potential.Corner(A, n, cut)`; ⚠️ Ch. 9 wedge flows U ∝ x^m |
| N45 | Eq. (6.47): vortex $w=-\frac{i\Gamma}{2\pi}\ln(z-z')=\frac{\Gamma}{2\pi}\theta'-i\frac{\Gamma}{2\pi}\ln r'$ [#71] | 6.4 | B | NOTE | C09 | stated (a-D24 given, §4c) with the complex-log primer (branch cut; φ of a vortex jumps by Γ) |
| N46 | Eq. (6.48): source $w=\frac{m}{2\pi}\ln(z-z')=\frac{m}{2\pi}\ln r'+i\frac{m\theta'}{2\pi}$ [#72] | 6.4 | B | NOTE | C09 | stated beside N45 |
| N47 | Eq. (6.49): doublet $w=\frac{d}{2\pi(z-z')}$ (dipole strength −d e_x) [#73] | 6.4 | B | NOTE | C09 | stated with Re w = the doublet φ of (6.29); ⚠️ scalar d here is a dipole −d e_x (`Doublet.from_book_scalar`) |
| N48 | Eq. (6.50): half-body $w=Uz+\frac{m}{2\pi}\ln z$ [#74] | 6.4 | B | NOTE | C09 | stated: Im w reproduces (6.31) |
| N49 | Eq. (6.51): cylinder $w=U(z+a^2/z)$ [#75] | 6.4 | B | NOTE | C09 | stated: Im w reproduces (6.33); pointer to C11 (circle plane) |
| N50 | Eq. (6.52): cylinder with clockwise circulation $w=U(z+a^2/z)+\frac{i\Gamma}{2\pi}\ln(z/a)$ [#76] | 6.4 | B | NOTE | C09 | stated: Im w reproduces (6.36); the input of C10 and C11 |
| N51 | Eq. (6.53): images in complex form $w=\frac{m}{2\pi}\ln(\frac{z^2-a^2}{a^2})$ [#77] | 6.4 | B | NOTE | C09 | stated: ln A + ln B = ln AB modulo 2πi; Im w = (6.41) ψ mod m |
| N52 | Complex variables give very general force results [#78] | 6.4 | C | NOTE | C09 | named as the bridge sentence into C10 |
| N53 | Setting of §6.5: stationary body of span B; drag D (x) and lift L (y) per unit depth on the body; force on the fluid $\mathbf F=-B(D\mathbf e_x+L\mathbf e_y)$ [#79] | 6.5 | B | NOTE | C10 | stated with the force-direction callout (on body vs on fluid); ⚠️ "Section 3" in the book means §6.3; B = span, not ch04's Bernoulli function |
| R18 | Eq. (6.54): steady momentum for a stationary CV, $\int_{A^*}\rho\mathbf u(\mathbf u\cdot\mathbf n)dA=-\int_{A^*}p\mathbf n\,dA+\mathbf F$ [#80] | 6.5 | B | RECAP | C10 | ch04 C04 (4.17), `core.conservation.momentum_budget`: reminded as D16 step 1; `ch06.cv_force_on_body` (large-circle route, Exercise 6.27) |
| N54 | Eq. (6.55): on the body u·n = 0, so $D\mathbf e_x+L\mathbf e_y=-\frac1B\int_{A^*}p\mathbf n\,dA$ [#81] | 6.5 | B | NOTE | C10 | stated as D16 step |
| N55 | Eq. (6.56): with the outward normal $(\mathbf e_xdy-\mathbf e_ydx)/ds$ of a counterclockwise contour, $D=-\oint_Cp\,dy$, $L=\oint_Cp\,dx$ (Fig. 6.18) [#82] | 6.5 | B | NOTE | C10 | stated as D16's result; `ch06.contour_force` (orientation asserted by the signed area) |
| N56 | Eq. (6.57): the complex force $D-iL=-i\oint_Cp\,dz^*$ [#83] | 6.5 | B | NOTE | C10 | stated as D17 step; `ch06.complex_force_from_pressure` |
| N57 | Bernoulli in complex form $p_\infty+\tfrac12\rho U^2=p+\tfrac12\rho(u-iv)(u+iv)$ [#84] | 6.5 | C | NOTE | C10 | named as D17 step 3 (∣q∣² = q q*) |
| N58 | Eq. (6.58): $D-iL=-i\oint_C[p_\infty+\tfrac12\rho U^2-\tfrac12\rho(u-iv)(u+iv)]dz^*$ [#85] | 6.5 | B | NOTE | C10 | stated as D17 step (∮ const dz* = 0) |
| N59 | Eq. (6.59): on the body the velocity is tangent, so $(u+iv)dz^*=(u-iv)dz=(dw/dz)dz$ [#86] | 6.5 | B | NOTE | C10 | stated as D17 step; asserted pointwise on the cylinder in code |
| N60 | Eq. (6.60): Blasius's theorem $D-iL=\frac{i\rho}{2}\oint_C(dw/dz)^2dz$; the contour may be moved to any curve enclosing the body if (dw/dz)² is analytic in between [#87] | 6.5 | B | NOTE | C10 | stated as D17's result (written out in C10); `core.potential.blasius_force` flat in R for the cylinder |
| N61 | Far field of any body with clockwise circulation Γ: $w=Uz+\frac{m}{2\pi}\ln z+\frac{i\Gamma}{2\pi}\ln z+\frac{d}{2\pi z}+\dots$, m = 0 for a closed body [#88] | 6.5 | B | NOTE | C10 | stated as D18 step (Laurent primer); `core.potential.laurent_coefficients` (FFT on a circle: c₀ = U, c₋₁ = iΓ_cw/2π) |
| N62 | Eq. (6.61): Blasius with the far-field series; the book prints the 1/z² coefficient as $(Ud/\pi-\Gamma^2/4\pi^2)$ and an extra outer square — the correct coefficient is $-(Ud/\pi+\Gamma^2/4\pi^2)$, and only the 1/z term, $iU\Gamma/\pi$, matters [#89] | 6.5 | B | NOTE | C10 | stated as D18 steps with the slip taught ("book prints X; Y is right", sympy expansion); the printed coefficient is never asserted |
| N63 | Residue theorem $\oint_Cf\,dz=2\pi i\sum\mathrm{Res}$; here Res at z = 0 is iUΓ/π; $\oint z^{-n}dz=2\pi i\,\delta_{n1}$ on a circle [#90] | 6.5 | B | NOTE | C10 | stated with the one-line circle computation ∫ i e^{i(1−n)θ} dθ and a numeric check (periodic trapezoid) before D18 |
| C10 | Eq. (6.62): Kutta–Zhukhovsky lift theorem for any cross-section, $D-iL=\frac{i\rho}{2}2\pi i(\frac{iU\Gamma}{\pi})=-i\rho U\Gamma$, i.e. D = 0 and L = ρUΓ, via Blasius's theorem (6.60) [#91] | 6.5 | A | CORE | – | load-bearing: generalises C06 and C07 to every 2-D body (d'Alembert in general, lift per circulation) — the foundation of airfoil theory (Ch. 14); the contour-deformation and residue moves recur in Ch. 7 and Ch. 11; D16, D17, D18 |
| N64 | Eq. (6.63): conformal map $\delta w=\frac{dw}{dz}\delta z$ — each small element is scaled by ∣dw/dz∣ and turned by arg(dw/dz) [#92] | 6.6 | B | NOTE | C11 | stated as D19 step; `core.conformal.map_elements` |
| N65 | Eq. (6.64): a second element $\delta'w=\frac{dw}{dz}\delta'z$ is turned by the same angle, so α = β; fails where dw/dz = 0 or ∞; large figures distort (Fig. 6.19) [#93] | 6.6 | B | NOTE | C11 | stated as D19's result; `core.conformal.angle_preservation` (z² doubles angles at 0) |
| N66 | A rectangular (φ, ψ) grid in the w-plane maps to the flow net in z (Fig. 6.20); intermediate maps w = ln ζ, ζ = sin z ⇒ u − iv = cot z; w = z² is the 90° corner [#94] | 6.6 | B | NOTE | C11 | stated with `core.conformal.grid_image` and `ch06.cot_flow`; ⚠️ the book's "(see Figure 6.5)" for the hyperbolae means Fig. 6.3 |
| C11 | Eq. (6.65): the Zhukhovsky transformation $z=\zeta+b^2/\zeta$ — identity far away, circle ∣ζ∣ = b → slit 2b cos θ, circle a > b → ellipse — used with angle preservation (6.63)–(6.64) to carry the circle flow (6.68) onto an elliptic cylinder through the correct inverse (6.69) [#95] | 6.6 | A | CORE | – | load-bearing: the method that turns the solved circle into new bodies, and the direct ancestor of Zhukhovsky airfoils and the Kutta condition in Ch. 14; D19, D20, D21 |
| N67 | Eq. (6.66): a circle of radius a > b maps to $z=ae^{i\theta}+\frac{b^2}{a}e^{-i\theta}$ [#96] | 6.6 | B | NOTE | C11 | stated as D20 step |
| N68 | Eq. (6.67): the ellipse $\frac{x^2}{(a+b^2/a)^2}+\frac{y^2}{(a-b^2/a)^2}=1$, foci at ±2b [#97] | 6.6 | B | NOTE | C11 | stated as D20's result; number b = 1, a = 1.2 → semi-axes 2.033 and 0.367, foci ±2; `ch06.joukowski_ellipse` |
| N69 | Eq. (6.68): circle flow with clockwise circulation in the ζ-plane $w=U(\zeta+a^2/\zeta)+\frac{i\Gamma}{2\pi}\ln(\zeta/a)$ [#98] | 6.6 | B | NOTE | C11 | stated (it is (6.52) with z → ζ); `core.conformal.mapped_flow` with `core.potential.cylinder` |
| N70 | Eq. (6.69): inverse map $\zeta=\tfrac12z+\tfrac12(z^2-4b^2)^{1/2}$, root outside the circle; $u-iv=\frac{dw}{d\zeta}\frac{d\zeta}{dz}$ [#99] | 6.6 | B | NOTE | C11 | stated as D21's result; ⚠️ numpy's principal √(z² − 4b²) picks the inside root for Re z < 0 (∣ζ∣ = 0.37 instead of 2.70 at z = −3 + 0.5i, b = 1) — use ½[z + √(z − 2b)√(z + 2b)]; `core.conformal.joukowski_inverse`, `ch06.elliptic_cylinder_flow` |
| N71 | Numerical routes for complex geometry: singularity distributions, thin-body perturbation, numerical Laplace [#100] | 6.7 | C | NOTE | C12 | named with pointers to C14 (singularity distributions, N90 panels) and Ch. 10 (CFD) |
| R19 | Grid $\psi_{i,j}=\psi(i\Delta x,j\Delta y)$ and half-point central differences $(\partial\psi/\partial x)_{i,j}\simeq(\psi_{i+\frac12,j}-\psi_{i-\frac12,j})/\Delta x$ (Fig. 6.22) [#101] | 6.7 | B | RECAP | C12 | ch02 C09 (`core.grids`, `core.operators.partial`): reminded; ⚠️ the book's "first-order" means first-*derivative*: the error is O(Δx²) (D22) |
| R20 | Eq. (6.70): $(\partial^2\psi/\partial x^2)_{i,j}\simeq(\psi_{i+1,j}-2\psi_{i,j}+\psi_{i-1,j})/\Delta x^2$ [#102] | 6.7 | B | RECAP | C12 | ch02 C09/C10 (`core.operators.laplacian`, 5-point): reminded as D22 steps; `core.laplace_solvers.laplacian_5pt` (masked grid) observed order 2 |
| R21 | Eq. (6.71): $(\partial^2\psi/\partial y^2)_{i,j}\simeq(\psi_{i,j+1}-2\psi_{i,j}+\psi_{i,j-1})/\Delta y^2$ [#103] | 6.7 | B | RECAP | C12 | as R20 |
| C12 | Eq. (6.72): with Δx = Δy, the finite-difference Laplace equation says every value is the average of its four neighbours, $\psi_{i,j}=\tfrac14[\psi_{i-1,j}+\psi_{i+1,j}+\psi_{i,j-1}+\psi_{i,j+1}]$, solved by Gauss–Seidel relaxation (Example 6.2) [#104] | 6.7 | A | CORE | – | load-bearing: the book's first numerical PDE solver — the ancestor of every elliptic solver in Ch. 10 (pressure Poisson, ω–ψ) and of bounded-domain PV inversion in Ch. 13; D22, D23 |
| N72 | Eq. (6.73): the 16-point grid (Fig. 6.23): 12 boundary values ψ^B, four unknowns, four linear equations [#105] | 6.7 | B | NOTE | C12 | stated as D23's start; number: boundary ψ = xy on a unit grid → interior 1, 2, 2, 4 exactly (ψ(1,1) = ¼(0 + 2 + 0 + 2) = 1); `ch06.four_point_system` (A, b) and `np.linalg.solve` |
| N73 | Gauss–Seidel iteration: sweep using the latest value at each point, $\psi^{(k+1)}_{i,j}=\tfrac14[\psi^{(k+1)}_{i-1,j}+\psi^{(k)}_{i+1,j}+\psi^{(k+1)}_{i,j-1}+\psi^{(k)}_{i,j+1}]$; stop on the residual, not only the change [#106] | 6.7 | B | NOTE | C12 | stated as D23 steps; `core.laplace_solvers.solve_laplace(method=…)` with Jacobi and SOR ghosts and residual history |
| N74 | Example 6.2: flow through a sharp contraction — linear ψ at inlet and outlet, ψ = 0 on the lower wall, ψ = Q on the upper; singly connected ⇒ unique (Figs. 6.24–6.25) [#107] | 6.7 | B | NOTE | C12 | stated (a-D33 given, §4c); `ch06.example_6_2` book grid + refinement; the 270° re-entrant corner (infinite speed, α = 3π/2 in (6.46)) lowers the convergence order near it; ⚠️ the FORTRAN loop runs over I but should run over J (harmless), Δψ units m²/s; book grid values only in `tests/book_values_ch06.json` |
| R22 | Two stream functions describe 3-D flow; axisymmetric flow needs one [#108] | 6.8 | C | RECAP | C13 | ch04 C03 (§4.3): one sentence opening C13 |
| R23 | Eq. (6.74): axisymmetric continuity $\frac1R\frac{\partial}{\partial R}(Ru_R)+\frac{\partial u_z}{\partial z}=0$ [#109] | 6.8 | B | RECAP | C13 | ch04 C03 axisymmetric mode (`core.streamfunction`): reminded with its numbers; ⚠️ z is now horizontal, along the stream |
| R24 | Eq. (6.75): Stokes stream function $u_R=-\frac1R\frac{\partial\psi}{\partial z}$, $u_z=\frac1R\frac{\partial\psi}{\partial R}$ (from χ = −φ, $\mathbf u=\nabla\chi\times\nabla\psi$) [#110] | 6.8 | B | RECAP | C13 | ch04 C03 (`core.streamfunction.velocity_from_streamfunction_axisym`): reminded as D24 step 1 |
| R25 | Eq. (6.76): $\omega_\varphi=\frac{\partial u_R}{\partial z}-\frac{\partial u_z}{\partial R}$ [#111] | 6.8 | B | RECAP | C13 | App. B curl (`core.curvilinear`), ch05 Hill ω_φ: reminded as D24 step |
| N75 | Eq. (6.77): irrotational Stokes ψ satisfies $\frac{\partial}{\partial R}(\frac1R\frac{\partial\psi}{\partial R})+\frac1R\frac{\partial^2\psi}{\partial z^2}=-\omega_\varphi=0$ — **not** the Laplacian, so complex variables do not apply [#112] | 6.8 | B | NOTE | C13 | stated as D24's result; `ch06.stokes_operator_residual` = 0 for (6.86)–(6.89), ≠ 0 for R²z (control) |
| N76 | Stokes ψ has units m³/s (plane ψ m²/s); ψ = const are surfaces of revolution [#113] | 6.8 | C | NOTE | C13 | named in the units line of C13 |
| N77 | Eq. (6.78): flow rate between neighbouring stream surfaces $dQ=2\pi R(\mathbf u\cdot\mathbf n)ds=2\pi\,d\psi$ (Fig. 6.26) [#114] | 6.8 | B | NOTE | C13 | stated (a-D35 given, §4c) with `ch06.axisym_flux_between` (quadrature = 2πΔψ); ⚠️ 2π, not 1 |
| R26 | Eq. (6.79): $u_R=\partial\phi/\partial R$, $u_z=\partial\phi/\partial z$ [#115] | 6.8 | B | RECAP | C13 | ch03 N29 (3.17) in cylindrical form: one line |
| N78 | Eq. (6.80): axisymmetric Laplace $\frac1R\frac{\partial}{\partial R}(R\frac{\partial\phi}{\partial R})+\frac{\partial^2\phi}{\partial z^2}=0$ [#116] | 6.8 | B | NOTE | C13 | stated beside N75 (φ obeys the true Laplacian, ψ does not); `ch06.axisym_laplacian_residual` |
| R27 | Eq. (6.81): cylindrical (R, φ, z) and spherical (r, θ, φ) coordinates, z-axis horizontal in this section [#117] | 6.8 | B | RECAP | C13 | ch03 (`core.coords`), P88 unit vectors: reminded with a two-column table |
| R28 | Eq. (6.82): spherical axisymmetric continuity (book prints r × the App. B divergence — correct) [#118] | 6.8 | B | RECAP | C13 | App. B (`core.curvilinear`): reminded; the book's form is r × App. B (correct; not a slip); code uses the App. B normalisation |
| N79 | Eq. (6.83): $u_r=\frac1{r^2\sin\theta}\frac{\partial\psi}{\partial\theta}=\frac{\partial\phi}{\partial r}$, $u_\theta=-\frac1{r\sin\theta}\frac{\partial\psi}{\partial r}=\frac1r\frac{\partial\phi}{\partial\theta}$ [#119] | 6.8 | B | NOTE | C13 | stated (a-D36 given, §4c) with a sympy check; `core.potential.axisym_velocity_spherical` |
| R29 | Eq. (6.84): $\omega_\varphi=\frac1r[\frac{\partial}{\partial r}(ru_\theta)-\frac{\partial u_r}{\partial\theta}]$ [#120] | 6.8 | B | RECAP | C13 | App. B spherical curl: one line, used in the residual tests |
| R30 | Eq. (6.85): spherical axisymmetric Laplace $\frac1{r^2}\frac{\partial}{\partial r}(r^2\frac{\partial\phi}{\partial r})+\frac1{r^2\sin\theta}\frac{\partial}{\partial\theta}(\sin\theta\frac{\partial\phi}{\partial\theta})=0$ [#121] | 6.8 | B | RECAP | C13 | App. B (`core.curvilinear`): reminded, used in D25's check |
| N80 | Eq. (6.86): uniform flow along z, $\phi=Uz$, $\psi=\tfrac12UR^2$; $\phi=Ur\cos\theta$, $\psi=\tfrac12Ur^2\sin^2\theta$ [#122] | 6.8 | B | NOTE | C13 | stated as D25 step 1; `core.potential.AxisymUniform` |
| N81 | Eq. (6.87): 3-D point source Q [m³/s], $\phi=-\frac{Q}{4\pi r}$, $\psi=-\frac{Q}{4\pi}\cos\theta$ (and the (R, z) forms) [#123] | 6.8 | B | NOTE | C13 | stated as D25 steps (flux through a sphere); `core.potential.PointSource3D` |
| N82 | Eq. (6.88): 3-D doublet with dipole −d e_z, $\phi=\frac{d}{4\pi r^2}\cos\theta$, $\psi=-\frac{d}{4\pi r}\sin^2\theta$ [#124] | 6.8 | B | NOTE | C13 | stated as D25 step (3-D limit); `core.potential.Doublet3D` |
| N83 | Eq. (6.89): sphere = stream + opposing doublet d = 2πa³U, $\psi=\tfrac12Ur^2(1-a^3/r^3)\sin^2\theta$, $\phi=Ur(1+\frac{a^3}{2r^3})\cos\theta$ (Fig. 6.27) [#125] | 6.8 | B | NOTE | C13 | stated as D25's construction; parity with the Hill exterior (`core.vortices.hill_stream_function`, ch05 N23) and the ch04 accelerating-sphere φ |
| C13 | Eq. (6.90): the sphere in a uniform stream, $u_r=U[1-(a/r)^3]\cos\theta$, $u_\theta=-U[1+\tfrac12(a/r)^3]\sin\theta$, from the Stokes stream function (6.75) whose field equation (6.77) is not Laplace's — axisymmetric ideal flow and its 3-D elements [#126] | 6.8 | A | CORE | – | load-bearing: the sphere is the 3-D reference body (Stokes-flow comparison in Ch. 8, separation in Ch. 9, added mass in C15) and the Stokes ψ is the tool for every axisymmetric flow (Ch. 13 overturning streamfunction); D24, D25 |
| N84 | Eq. (6.91): sphere surface $C_p=1-(u_\theta/U)^2=1-\tfrac94\sin^2\theta$; max speed 1.5U; fore–aft symmetric ⇒ no drag [#127] | 6.8 | B | NOTE | C13 | stated as D25's result with the comparison table cylinder (2U, C_p min −3, decay 1/r²) vs sphere (1.5U, −1.25, 1/r³) — 3-D relief; `ch06.sphere_surface_cp` |
| N85 | Eq. (6.92): coordinate-free sphere potential $\phi=(\mathbf U-\frac{\mathbf d}{4\pi\lvert\mathbf x\rvert^3})\cdot\mathbf x$ [#128] | 6.8 | B | NOTE | C13 | stated (a-D39 given, §4c: cos θ/r² = e_z·x/∣x∣³); `core.potential.sphere_potential_vector` — the input of C15 |
| N86 | Closed bodies of revolution: stream + sources and sinks of zero net strength; a distributed sink gives a streamlined tail (Fig. 6.28) [#129] | 6.8 | C | NOTE | C14 | named as C14's opening rule (Σ Q_i = 0), checked by `airship` and `axial_singularity_solve` |
| N87 | Eq. (6.93): uniform line sink of density k on the axis from O to A: $d\psi_{\text{sink}}=\frac{k\,d\xi}{4\pi}\cos\alpha$, $\psi_{\text{sink}}=\frac{k}{4\pi}\int_0^a\cos\alpha\,d\xi$ [#130] | 6.8 | B | NOTE | C14 | stated as D26's start; `ch06.line_sink_stream_function(method="quad")` |
| N88 | Eq. (6.94): closed form by z − ξ = R cot α, $\psi_{\text{sink}}=\frac{kR}{4\pi}[\frac1{\sin\theta}-\frac1{\sin\alpha_1}]=\frac{k}{4\pi}(r-r_1)$ [#131] | 6.8 | B | NOTE | C14 | stated as D26 steps; closed form = `quad` to 1e-10 |
| N89 | Eq. (6.95): airship = point source Q at O + line sink k = Q/a + uniform stream, $\psi=-\frac{Q}{4\pi}\cos\theta+\frac{Q}{4\pi a}(r-r_1)+\tfrac12Ur^2\sin^2\theta$ (Fig. 6.28, both frames) [#132] | 6.8 | B | NOTE | C14 | stated as D26's result; number: Q = 1 m³/s, a = 1 m → k = 1 m²/s; body length from the axis stagnation points (`ch06.airship`, `brentq`) |
| C14 | The axial singularity (inverse) method: N axial segments of unknown density k_n, $\psi_m=-\sum_{n=1}^N\frac{k_n}{4\pi}(r^m_{n-1}-r^m_n)+\tfrac12UR_m^2=0$ at N body points — an N × N linear system that finds the sources for a given body of revolution (Fig. 6.29) [#133] | 6.8 | A | CORE | – | load-bearing: the chapter's design method (shape given, sources unknown) — the ancestor of panel and boundary-element methods in Ch. 10 and Ch. 14 and of slender-body theory; D26, D27 |
| N90 | Extension (not in the book): 2-D constant-strength source panels (Hess & Smith 1967) — u·n = 0 at panel midpoints, λ_i/2 self term, C_p on the body converging to 1 − 4 sin²θ for a circle [#134] | 6.7 | B | NOTE | C14 | stated inside C14 as "the same idea on the body surface in 2-D", labelled *our extension*; `core.panels.source_panels` (observed order 2 on the circle, Σλ_js_j = 0); the E8 panel mode; listed under §6.7 (the book's numerical section) |
| N91 | 3-D d'Alembert for closed bodies in steady flow (Exercise 6.39); unsteadiness gives forces; apparent (added) mass for vehicles, fish, bubbles [#135] | 6.9 | C | NOTE | C15 | named as C15's question; pointers Ch. 14 (vorticity-induced drag), Ch. 16 (swimming, bubbles) |
| R31 | 3-D potential u ≡ ∇φ with w ≡ ∂φ/∂z — w is the z-velocity here, not the complex potential [#136] | 6.9 | C | RECAP | C15 | ch03 N29 (3.17): one sentence with the ⚠️ symbol clash |
| N92 | Setup (Fig. 6.30): sphere of radius a at x_s(t), velocity $\mathbf u_s=d\mathbf x_s/dt$, known acceleration; external force F_E on the sphere; find the fluid force F_s [#137] | 6.9 | B | NOTE | C15 | stated with the sketch (ASCII) and the question |
| N93 | Eq. (6.96): moving-sphere potential from (6.92) with x → x − x_s(t), U = 0: $\phi=-\frac{1}{4\pi\lvert\mathbf x-\mathbf x_s\rvert^3}\mathbf d\cdot(\mathbf x-\mathbf x_s)$ [#138] | 6.9 | B | NOTE | C15 | stated as D28 step; `core.potential.moving_sphere_potential` |
| N94 | Eq. (6.97): the dipole follows the motion, $\mathbf d(t)=2\pi a^3\mathbf u_s(t)$, $\phi=-\frac{a^3}{2\lvert\boldsymbol\xi\rvert^3}\mathbf u_s\cdot\boldsymbol\xi$ [#139] | 6.9 | B | NOTE | C15 | stated as D28's result; check n·∇φ = n·u_s on ∣ξ∣ = a |
| N95 | Eq. (6.98): fluid force $\mathbf F_s=-\int(p-p_\infty)\mathbf n\,dA$ over the sphere [#140] | 6.9 | B | NOTE | C15 | stated (3-D analogue of (6.55)); `core.potential.sphere_force_quadrature` |
| R32 | Eq. (6.99): unsteady Bernoulli between the surface and far away, $[\frac{\partial\phi}{\partial t}+\frac12\lvert\nabla\phi\rvert^2+\frac p\rho]_{\text{surface}}=\frac{p_\infty}{\rho}$ [#141] | 6.9 | B | RECAP | C15 | ch04 C12 (4.75), D26: reminded as D29 step 1 with `core.bernoulli.unsteady_bernoulli_pressure` |
| N96 | Eq. (6.100): $\frac{p_a-p_\infty}{\rho}=-(\frac{\partial\phi}{\partial t})_a-\frac12\lvert\nabla\phi\rvert_a^2$ [#142] | 6.9 | B | NOTE | C15 | stated as D29 step |
| N97 | Eq. (6.101): ∂φ/∂t at a fixed point acts through x_s(t) and u_s(t): $\frac{\partial\phi}{\partial t}=-\mathbf u\cdot\mathbf u_s-\frac{a^3}{2\lvert\boldsymbol\xi\rvert^3}\boldsymbol\xi\cdot\frac{d\mathbf u_s}{dt}$ [#143] | 6.9 | B | NOTE | C15 | stated as D29 steps (∂φ/∂x_s = −∇φ); `ch06.moving_sphere_dphidt_sym` |
| N98 | Eq. (6.102): on the surface $(\frac{\partial\phi}{\partial t})_a=-\mathbf u_a\cdot\mathbf u_s-\frac a2\mathbf e_\xi\cdot\frac{d\mathbf u_s}{dt}$ [#144] | 6.9 | B | NOTE | C15 | stated as D29 step |
| N99 | Eq. (6.103): gradient of the dipole potential $\nabla\phi=-\frac{a^3}{2}[-\frac{3(\mathbf x-\mathbf x_s)}{\lvert\mathbf x-\mathbf x_s\rvert^5}\mathbf u_s\cdot(\mathbf x-\mathbf x_s)+\frac{\mathbf u_s}{\lvert\mathbf x-\mathbf x_s\rvert^3}]$ [#145] | 6.9 | B | NOTE | C15 | stated as D29 step; `core.potential.moving_sphere_velocity` |
| N100 | Eq. (6.104): surface velocity $\mathbf u_a=\tfrac32(\mathbf u_s\cdot\mathbf e_\xi)\mathbf e_\xi-\tfrac12\mathbf u_s$ (the book's middle bracket prints −u_s/a³; +u_s/a³ is right) [#146] | 6.9 | B | NOTE | C15 | stated as D29 step with the slip taught (sympy: the printed bracket gives +½u_s); physical check: at the sphere's side the fluid moves at −½u_s |
| N101 | Eq. (6.105): surface pressure on an arbitrarily moving sphere $\frac{p_a-p_\infty}{\rho}=\frac12\lvert\mathbf u_s\rvert^2(\frac94\frac{(\mathbf u_s\cdot\mathbf e_\xi)^2}{\lvert\mathbf u_s\rvert^2}-\frac54)+\frac a2\mathbf e_\xi\cdot\frac{d\mathbf u_s}{dt}$ [#147] | 6.9 | B | NOTE | C15 | stated as D29's result (★★★, written out in C15); parity with `ch04.accelerating_sphere_pressure` and oblique acceleration (new); `core.potential.moving_sphere_surface_pressure` |
| N102 | Eq. (6.106): no acceleration ⇒ $\frac94\cos^2\theta_s-\frac54=1-\frac94\sin^2\theta_s$ = (6.91): Galilean invariance, no drag [#148] | 6.9 | B | NOTE | C15 | stated (a-D45 given, §4c); ⚠️ θ_s is measured from the sphere's velocity (= π − θ of (6.91)) |
| N103 | Eq. (6.107): force with the acceleration along e_z, $\mathbf F_s=-\rho\frac a2\lvert\frac{d\mathbf u_s}{dt}\rvert\int\int\cos\theta(\mathbf e_x\sin\theta\cos\varphi+\mathbf e_y\sin\theta\sin\varphi+\mathbf e_z\cos\theta)a^2\sin\theta\,d\varphi\,d\theta$ [#149] | 6.9 | B | NOTE | C15 | stated as D30 step (φ-integral kills x, y) |
| N104 | Eq. (6.108): the added mass of a sphere, $\mathbf F_s=-\frac23\pi\rho a^3\frac{d\mathbf u_s}{dt}=-M\frac{d\mathbf u_s}{dt}$, $M=\frac{2\pi a^3\rho}{3}$ = half the displaced mass (the printed integrand keeps a stray dφ) [#150] | 6.9 | B | NOTE | C15 | stated as D30's result and confirmed by D31 (kinetic energy); number: a = 0.1 m in water → M = 2.094 kg, displaced 4.189 kg; cylinder: ρπa² per depth = the whole displaced mass (a-D48 given, §4c); `core.potential.added_mass_sphere`, `ch06.added_mass_by_energy` |
| N105 | Why: the fluid ahead must be pushed aside faster and the fluid behind filled in [#151] | 6.9 | C | NOTE | C15 | named in the interpretation after D30 |
| C15 | Eq. (6.109): Newton's law for a submerged sphere, $\mathbf F_E+\mathbf F_s=\mathbf F_E-M\frac{d\mathbf u_s}{dt}=m\frac{d\mathbf u_s}{dt}$, i.e. $\mathbf F_E=(m+\frac{2\pi}{3}\rho a^3)\frac{d\mathbf u_s}{dt}$ — the sphere behaves as if heavier by half the displaced fluid; apparent mass is in general a tensor [#152] | 6.9 | A | CORE | – | load-bearing: the one force ideal flow predicts on a body in 3-D (unsteady inertia) — bubbles, fish, ships, oscillating structures (Ch. 16), and the unsteady side of Ch. 14; closes the chapter's force story (d'Alembert → lift → added mass); D28, D29, D30, D31 |
| N106 | Concluding remarks: history (Euler … Kelvin), potential theory in other fields, zero drag vs observation, viscous theory (Prandtl), airfoils and conformal maps [#153] | 6.10 | C | NOTE | C15 | named in the closing cell with pointers Ch. 9 (Prandtl, boundary layers) and Ch. 14 (airfoils, Kutta condition) |
| S01 | Exercises 6.1–6.5, 6.8, 6.33 — checks of §6.2 (δ-source integrals, ∇ψ·∇φ = 0, e_z × ∇ψ = ∇φ, orthogonality, Bernoulli from (6.1)) [#154] | Ex. | C | SKIP | C02 | pointer: their ideas are D03 and D04 (tests of N06, N09); exercise text not reproduced |
| S02 | Exercises 6.6–6.7, 6.9, 6.14–6.15, 6.17, 6.22 — harmonic polynomials, doublet ψ, perturbed cylinder, sketched flows, vortex-pair doublet limit [#155] | Ex. | C | SKIP | C04 | pointer: 6.6/6.7 generated by `ch06.harmonic_polynomials` (N20), 6.9's result stated in N25/§4c; rest left to the reader |
| S03 | Exercises 6.10–6.13, 6.18–6.20, 6.23–6.24 — forces on held singularities, half-body width and zero drag, Rankine oval, vortex over a plate, source above a plane [#156] | Ex. | C | SKIP | C05 | pointer: 6.12 is D07's mass balance, 6.13 is N29's net force, the Rankine oval is the E1 preset (`ch06.rankine_oval`); others left to the reader |
| S04 | Exercises 6.21, 6.25–6.27, 6.29–6.32 — Blasius for the cylinder, cylinder in a corner, circle images, CV proof of Kutta–Zhukhovsky, vortex pairs and corners [#157] | Ex. | C | SKIP | C10 | pointer: 6.21 is the E5 cylinder preset, 6.27 is `ch06.cv_force_on_body` (R18); others left to the reader |
| S05 | Exercises 6.16, 6.41 — kitchen experiments (paper cylinder vs airfoil; vacuum nozzle over sugar grains) [#158] | Ex. | C | SKIP | C08 | pointer: 6.41's sink + image is the 3-D `mirror`; try them at home |
| S06 | Exercises 6.34–6.40, 6.42–6.43 — 3-D source and doublet limits, 3-D half-body, perturbation radii, 3-D d'Alembert, hemisphere lift-off, airship length, line-source half-body [#159] | Ex. | C | SKIP | C13 | pointer: 6.34/6.37 are D25, 6.42 is D26's length equation, 6.39 is N91; others left to the reader |
| S07 | Exercises 6.44–6.49 — cavity collapse, cylinder added mass, bubble acceleration, oscillation frequencies, unsteady stream past a sphere, kinetic-energy added mass [#160] | Ex. | C | SKIP | C15 | pointer: 6.45 stated in N104, 6.46 is the E9 bubble preset (2g), 6.49 is D31; the rest left to the reader |

### 2a. What each A block contains (for the lesson-designer)
- **C01** picture: air over a car — a thin skin of slowed air on the body, everything outside it behaving as if
  frictionless · question: "water is viscous — how can a theory with no friction describe it at all?" · D01 · number:
  car at 10 m/s, L = 1 m, air ν = 1.5×10⁻⁵ m²/s → Re ≈ 6.7×10⁵, boundary layer ~ 1.2 mm; cylinder potential flow:
  μ∇²u = 0 to round-off while the viscous *stress* is not zero · code: `ch06.ideal_flow_residuals` (cylinder vs plane
  Poiseuille), `ch06.ideal_flow_applicability` + from-scratch central differences of μ∇²u at a point · figure: two
  panels — μ∇²u heatmap for the cylinder (zero) vs Poiseuille (nonzero); applicability table · B/C: N01–N04, R01.
- **C02** picture: a bathtub vortex and a garden sprinkler seen from above — one stirs, one spills · question: "what
  single equation do all irrotational 2-D flows obey, and what makes a vortex or a source special?" · D02, D03, D04 ·
  number: vortex Γ = 2π m²/s → u_θ = 1 m/s at r = 1 m, ∮∇ψ·n ds = −2π on every circle; source m = 2π m²/s → u_r = 1/r ·
  code: `ch06.vorticity_from_psi`, `core.potential.laplacian_residual`, `ch06.delta_flux_check`, `ch06.orthogonality_check`,
  `core.potential.Uniform/Source/Vortex` + from-scratch 5-point −∇²ψ on the Rankine core vs `vorticity_from_psi` ·
  figure: flow nets (ψ solid, φ dashed) for stream, source, vortex; flux-vs-radius plot (flat); ω = −∇²ψ heatmap of a
  Rankine vortex (zero outside the core) · B/C: N05–N15, N18, N19, R02–R04, R06–R09.
- **C03** picture: a stream meeting a hidden obstacle — dye lines bend around a shape nobody drew · question: "if every
  element solves Laplace, where does the body come from?" · D05 · number: U = 1 m/s + source m = 2π m²/s at (−1, 0) m:
  u = 1 + (1/2π)(2π)(−1)/1 = 0 → a stagnation point (preview of C05) · code: `core.potential.Flow` (sum of elements),
  `normal_velocity_on`, `far_field_check`, `Flow.pressure` + from-scratch sum of element velocities at a point vs
  `Flow.velocity` · figure: three panels stream + source = sum, with the ψ contour through the stagnation point drawn as
  a wall · explainer E1 · B/C: N16, N17, R05.
- **C04** picture: two garden hoses, one blowing out and one sucking in, pushed nose to nose · question: "what is left
  when a source and a sink merge?" · D06 · number: m = 50 m²/s, ε = 0.02 m → ∣d∣ = 2mε = 2 m³/s; at (1, 0) m the pair's φ
  differs from the doublet's by ~ε² relative · code: `ch06.source_sink_pair`, `ch06.doublet_limit_error` (slope 2 on
  log–log), `core.potential.Doublet`, `Corner` + from-scratch pair formula vs `Doublet` at shrinking ε · figure: the
  element gallery (uniform, source, vortex, corner, doublet — streamlines circles tangent at the origin) · animation
  (§6) · explainer E1 · B/C: N20–N26, R10, R11.
- **C05** picture: river water parting round a bridge pier's rounded nose · question: "where does the stream stop,
  and how wide does the body get?" · D07, D08 · number: U = 1 m/s, m = 2π m²/s → a = 1 m, h(θ = 90°) = π/2 ≈ 1.571 m,
  h_max = π ≈ 3.142 m; C_p = 0 on the body at our root ≈ 113.2° · code: `core.potential.half_body`,
  `ch06.half_body_shape`, `half_body_numbers`, `half_body_surface_cp`, `half_body_cp_zero_angle` + from-scratch `brentq`
  on tan θ + 2(π − θ) · figure: ψ contours with the dividing streamline and stagnation point; C_p along the body vs θ with
  the zero marked; plotly slider over m/U · explainer E1 · B/C: N27–N29.
- **C06** picture: wind round a chimney · question: "the air pushes hard on the front — why is there no net push?" ·
  D09 · number: air, U = 10 m/s, a = 0.1 m → ½ρU² = 60 Pa, front stagnation +60 Pa, shoulders −180 Pa, back +60 Pa;
  D = 0 to 1e-13 of ½ρU²a from the spectral trapezoid · code: `core.potential.cylinder`, `ch06.cylinder_surface_cp`,
  `ch06.surface_pressure_force` + from-scratch trapezoid ∮ −p cos θ a dθ vs `surface_pressure_force` · figure: C_p(θ)
  ideal vs qualitative real band (Fig. 6.10 idea, angle from the front); pressure arrows round the body · explainer E2 ·
  B/C: N30, N31, R12–R14.
- **C07** picture: a spinning ball curving in flight, or a rotor ship's tall spinning cylinders · question: "how does
  circulation turn into sideways force, and how big is it?" · D10, D11 · number: U = 10 m/s, a = 0.1 m, Γ = 2 m²/s,
  air → stagnation points at −9.16° and −170.84°, L = 1.2 × 10 × 2 = 24 N/m; merge at Γ = 12.57 m²/s · code:
  `core.potential.cylinder(Gamma_cw=…)`, `ch06.cylinder_stagnation_points`, `cylinder_surface_pressure`, `lift_per_span`,
  `circulation_family_check` + from-scratch trapezoid lift integral vs ρUΓ · figure: four panels Γ = 0, 2πaU, 4πaU,
  6πaU (Fig. 6.12 idea); plotly slider over Γ (surface C_p, stagnation angles, L) · explainer E2 · B/C: N32–N38.
- **C08** picture: an eddy sweeping along a seabed or a wall — a pressure sensor on the wall feels it pass · question:
  "what does the wall feel as a vortex goes by?" · D12, D13 · number: water, Γ = 1 m²/s, h = 1 m → drift 0.0796 m/s,
  p − p∞ = −25.33 Pa at t = 0, zero at t = 4πh²/Γ = 12.57 s, maximum +3.17 Pa at t = 4√3πh²/Γ = 21.77 s · code:
  `core.potential.mirror`, `ch06.two_sources`, `ch06.example_6_1` (closed form + numeric route with ∂φ/∂t by central
  differences in t) · figure: vortex path beside the wall and p(0, 0, t) with the markers; two-source streamlines with
  the (6.41) curve overlaid · animation (§6) · explainer E3 · B/C: N39, N40, R15–R17.
- **C09** picture: one complex function whose real part draws the equipotentials and whose imaginary part draws the
  streamlines · question: "why does a single complex function carry a whole flow, and why is the velocity u − iv?" ·
  D14, D15 · number: w = z² at z = 1 + i → dw/dz = 2 + 2i ⇒ u = 2, v = −2 m/s (ψ = 2xy agrees) · code: `Flow.w`,
  `Flow.dwdz`, `ch06.cauchy_riemann_residual`, `core.potential.Corner(A, n)` + from-scratch complex difference quotients
  along x and along iy vs `Flow.dwdz` · figure: flow nets of z^n for α = π/2, π, 3π/2, 2π; log–log ∣dw/dz∣ vs r near
  the corner (slope n − 1) · explainer E4 · B/C: N41–N52.
- **C10** picture: two very different wings (a round tube, a flat ellipse) with the same circulation feeling the same
  lift · question: "why does the lift not depend on the shape?" · D16, D17, D18 · number: cylinder U = 10 m/s, Γ = 2
  m²/s, ρ = 1.2 kg/m³: residue iUΓ/π = 6.366i, D − iL = −iρUΓ → L = 24 N/m by Blasius on R = a, 2a, 10a (flat) · code:
  `core.potential.blasius_force`, `laurent_coefficients`, `ch06.contour_force`, `cv_force_on_body`,
  `kutta_zhukhovsky_sym` + from-scratch periodic trapezoid of (dw/dz)² on a circle vs `blasius_force` · figure: force vs
  contour radius (flat) for the cylinder and the Zhukhovsky ellipse; Laurent coefficient bars · explainer E5 · B/C:
  N53–N63, R18.
- **C11** picture: a circle drawn on a rubber sheet that is stretched without tearing — small crosses stay crosses ·
  question: "how can the solved circle give the flow round other shapes?" · D19, D20, D21 · number: b = 1 m, a = 1.2 m →
  ellipse semi-axes 2.033 m and 0.367 m, foci ±2 m; at z = −3 + 0.5i the right inverse has ∣ζ∣ = 2.70, numpy's principal
  root 0.37 (inside — wrong) · code: `core.conformal.joukowski`, `joukowski_inverse`, `mapped_flow`,
  `angle_preservation`, `grid_image`, `ch06.joukowski_ellipse`, `elliptic_cylinder_flow` + from-scratch quadratic formula
  with the outside root picked by ∣ζ∣ > b vs `joukowski_inverse` · figure: ζ-plane and z-plane side by side (circle →
  ellipse with streamlines); w-grid → z-plane for z², e^z · explainer E6 · B/C: N64–N70.
- **C12** picture: a stretched rubber membrane over a bent wire frame — each point sits at the average height of its
  neighbours · question: "how does a grid of numbers solve Laplace's equation?" · D22, D23 · number: ψ = xy on a unit
  grid: ψ(1, 1) = ¼(0 + 2 + 0 + 2) = 1 exactly (a harmonic polynomial is exact for the 5-point stencil); Gauss–Seidel
  needs roughly half Jacobi's sweeps, SOR a small fraction · code: `core.laplace_solvers.laplacian_5pt`, `solve_laplace`
  (Jacobi, Gauss–Seidel, SOR, direct), `ch06.four_point_system`, `ch06.example_6_2` + from-scratch double-loop
  Gauss–Seidel on the 4-point system vs `np.linalg.solve` · figure: Example 6.2 streamlines in the contraction; residual
  vs sweep (log) for three methods; refinement study with the corner singularity fenced · animation (§6) · explainer E7 ·
  B/C: N71–N74, R19–R21.
- **C13** picture: air round a football (3-D) vs round a lamp-post (2-D) — the 3-D flow can also go *over* the body ·
  question: "what changes when the flow can escape sideways?" · D24, D25 · number: sphere max speed 1.5U (C_p −1.25) vs
  cylinder 2U (−3); perturbation at r = 2a: sphere (a/r)³ = 0.125 vs cylinder (a/r)² = 0.25 · code:
  `ch06.stokes_operator_residual`, `ch06.axisym_flux_between`, `core.potential.AxisymUniform/PointSource3D/Doublet3D/sphere`,
  `ch06.sphere_surface_cp` + from-scratch ψ superposition (stream + doublet) vs `sphere` · figure: meridian streamlines
  of the sphere; cylinder vs sphere C_p and velocity decay; the flux ring (2πdψ) · B/C: N75–N85, R22–R30.
- **C14** picture: an airship hull shaped by hiding a source in its nose and a spread-out sink along its body ·
  question: "given a shape, which sources draw it?" · D26, D27 · number: Q = 1 m³/s, a = 1 m, k = 1 m²/s, closure Q = ak;
  the axial method on a Rankine-oval target recovers the point source and sink as N grows (condition number printed) ·
  code: `ch06.line_sink_stream_function`, `ch06.airship`, `core.panels.axial_singularity_solve`,
  `core.panels.source_panels` + from-scratch N × N influence matrix built by a double loop vs `axial_singularity_solve` ·
  figure: airship in body and fluid frames (Fig. 6.28 idea); target body vs ψ = 0 contour with k_n bars; panel C_p on a
  circle converging to 1 − 4 sin²θ (N90) · explainer E8 · B/C: N86–N90.
- **C15** picture: pushing a beach ball under water and letting it go — or a bubble rising · question: "an ideal fluid
  gives no drag, so why does it resist being accelerated?" · D28, D29, D30, D31 · number: a = 0.1 m sphere in water:
  M = (2/3)π(1000)(0.1)³ = 2.094 kg (half of 4.189 kg); a steel ball (7800 kg/m³, m = 32.67 kg) behaves as 34.77 kg; an
  air bubble starts upward at 2g, not ∞ · code: `core.potential.moving_sphere_potential`, `moving_sphere_velocity`,
  `moving_sphere_surface_pressure`, `sphere_force_quadrature`, `added_mass_sphere`, `ch06.added_mass_by_energy`,
  `ch06.sphere_motion` + from-scratch Gauss–Legendre × trapezoid force integral vs `added_mass_sphere` · figure: surface
  pressure split steady / acceleration on the sphere; force vs time for an oscillating sphere with −M du/dt; bubble vs
  ball velocity with and without added mass · animation (§6) · explainer E9 · B/C: N91–N106, R31, R32.

## 3. Section coverage
| § | Title | A | B | C | RECAP | SKIP |
|---|---|---|---|---|---|---|
| 6.1 | Relevance of Irrotational Constant-Density Flow Theory | C01 | N01, N02 | N03, N04 | R01 (B) | — |
| 6.2 | Two-Dimensional Stream Function and Velocity Potential | C02, C03 | N05, N06, N07, N08, N09, N10, N11, N12, N13, N14, N16, N17, N18, N19 | N15 | R02 (C), R03 (B), R04 (C), R05 (B), R06 (C), R07 (C), R08 (C), R09 (C) | — |
| 6.3 | Construction of Elementary Flows in Two Dimensions | C04, C05, C06, C07, C08 | N20, N21, N22, N23, N24, N25, N27, N28, N29, N30, N31, N32, N33, N34, N35, N36, N38, N39, N40 | N26, N37 | R10 (B), R11 (B), R12 (B), R13 (B), R14 (B), R15 (C), R16 (B), R17 (B) | — |
| 6.4 | Complex Potential | C09 | N41, N42, N43, N44, N45, N46, N47, N48, N49, N50, N51 | N52 | — | — |
| 6.5 | Forces on a Two-Dimensional Body | C10 | N53, N54, N55, N56, N58, N59, N60, N61, N62, N63 | N57 | R18 (B) | — |
| 6.6 | Conformal Mapping | C11 | N64, N65, N66, N67, N68, N69, N70 | — | — | — |
| 6.7 | Numerical Solution Techniques in Two Dimensions | C12 | N72, N73, N74, N90 (taught in C14) | N71 | R19 (B), R20 (B), R21 (B) | — |
| 6.8 | Axisymmetric Ideal Flow | C13, C14 | N75, N77, N78, N79, N80, N81, N82, N83, N84, N85, N87, N88, N89 | N76, N86 | R22 (C), R23 (B), R24 (B), R25 (B), R26 (B), R27 (B), R28 (B), R29 (B), R30 (B) | — |
| 6.9 | Three-Dimensional Potential Flow and Apparent Mass | C15 | N92, N93, N94, N95, N96, N97, N98, N99, N100, N101, N102, N103, N104 | N91, N105 | R31 (C), R32 (B) | — |
| 6.10 | Concluding Remarks | — | — | N106 | — | S01, S02, S03, S04, S05, S06, S07 (end of chapter) |

Totals: A 15 · B 116 (93 NOTE + 23 RECAP) · C 29 (13 NOTE + 9 RECAP + 7 SKIP) = 160. No section is empty; every
section except §6.10 (two pages of remarks) has at least one A item, and §6.10 is covered by N106 at the end of C15's
block plus the exercise pointers. N90 (source panels, our extension) is listed under §6.7 — the book's numerical section
— but taught inside C14, where the singularity-distribution idea lives; the §6.7 notebook section names it with a
forward pointer.

## 4. Prerequisites needing primers (concept or tool | needed by | why it is not A/B/RECAP)
Rows marked **gloss** are one plain sentence where the item is used; all others are full 📎 primers (one or two plain
sentences, what it means here, a 2–4-line runnable demo with easy numbers). New primers continue at **P149**. Items
already in `knowledge/primers.md` get a one-line reminder naming the earlier primer (listed at the end). The complex-
analysis primers must come **before C09** (the analyst's teaching risk): place the complex-plane and complex-log primers
at the top of §6.4, Cauchy's theorem and Laurent series at the top of §6.5.

| Concept or tool | Needed by | Why it is not A/B/RECAP |
|---|---|---|
| 2-D divergence theorem for a gradient, $\oint\nabla f\cdot\mathbf n\,ds=\int\nabla^2f\,dA$, and the 2-D Dirac delta as "concentrated total" | C02 (D03) | ch02 C14 was 3-D and ch05 P139 did the 3-D Green's function; the 2-D flux-of-ln r version is new |
| logarithm rules $\ln\sqrt{s}=\tfrac12\ln s$, $\ln(ab)=\ln a+\ln b$, $d(\ln r)/dx=x/r^2$ | C02, C04 (D06) | **gloss** |
| limit with a product held fixed (m → ∞, ε → 0, 2mε = const) and why the order matters | C04 (D06) | new maths move (P26 Taylor reminder for ln(1 + s) ≈ s) |
| `np.arctan2` branch for the half-body ψ (θ ∈ (−π, π] with the cut inside the body) | C05 (D07) | P70 primed arctan2; the choice of cut for a multivalued ψ is new — **gloss** |
| trigonometric integrals over a period: $\int_0^{2\pi}\sin\theta\,d\theta=\int_0^{2\pi}\sin^3\theta\,d\theta=\int_0^{2\pi}\sin\theta\cos\theta\,d\theta=0$, $\int_0^{2\pi}\sin^2\theta\,d\theta=\pi$ (odd/even symmetry) | C06, C07 (D11) | new maths tool (orthogonality of sines) |
| the two solutions of sin θ = s (θ and π − θ) and choosing the root of a quadratic by a physical condition (r > a) | C07 (D10) | **gloss** (P108 `brentq` reminder for the numerical cross-check) |
| complex Newton iteration for zeros of dw/dz (stagnation points) with deflation | C07 (code) | new numerics idea — **gloss** with a 3-line demo |
| odd and even functions under a mirror, ψ(x, −y) = −ψ(x, y) makes y = 0 a streamline | C08 (D12) | **gloss** |
| derivative of arctan(Y/X) with a moving Y(t), $\frac{d}{dt}\tan^{-1}(Y/X)=\frac{X\dot Y}{X^2+Y^2}$ | C08 (D13) | **gloss** on top of the chain rule (P49) |
| complex plane: modulus, argument, $z=re^{i\theta}$; multiplying complex numbers multiplies moduli and adds arguments; numpy `1j`, `np.abs`, `np.angle`, `np.conj`, complex arrays | C09, C11 (D19) | P45 primed i² = −1 and Euler's formula, P81 the conjugate; geometry of multiplication and numpy's complex API are new |
| complex derivative as a limit along any direction; analytic function (same limit from every direction) | C09 (D14) | new maths tool (the heart of D14) |
| complex logarithm $\ln z=\ln r+i\theta$, multivalued θ and branch cuts; complex powers $z^n=e^{n\ln z}$ and rotating the cut out of the fluid | C09 (D15, N45–N47), C10 | new maths tool |
| Cauchy's integral theorem: ∮ of an analytic function round a closed curve is 0, so a contour can be deformed through regions where the integrand is analytic | C10 (D17) | new maths tool |
| Laurent series (powers of z including negative ones outside a disc) and reading off coefficients; FFT of samples on a circle gives them | C10 (D18, N61) | new maths tool (P142 FFT reminder) |
| parametrised closed curve and its outward normal $(\mathbf e_xdy-\mathbf e_ydx)/ds$ for counterclockwise traversal (signed area > 0) | C10 (D16) | P92 primed parametric curves; the orientation–normal link is new — **gloss** |
| sympy for complex series: `sp.series`, `sp.residue`, `sp.I`, `sp.expand` | C10 (D18 check) | Python idiom not yet primed (P40 sympy reminder) |
| trig identities $e^{i\theta}+e^{-i\theta}=2\cos\theta$, $\cos^2+\sin^2=1$ to eliminate a parameter; foci of an ellipse $c^2=A^2-B^2$ | C11 (D20) | **gloss** |
| quadratic formula in ℂ and the branch of a complex square root; why `np.sqrt` of a complex array returns the principal root; the product of the two roots | C11 (D21) | new maths/Python tool — the chapter's main numerical trap |
| matrix form Aψ = b of a stencil, diagonal dominance; Jacobi vs Gauss–Seidel vs SOR; convergence rate set by the spectral radius; stopping on the residual | C12 (D23) | new numerics (P57 `np.linalg.solve` and P80 eigenvalues reminders) |
| boolean masks for an L-shaped domain (`mask[j, i]`), fixed boundary values, `scipy.sparse` matrices and `scipy.sparse.linalg.spsolve` | C12 (code) | new Python tools (P46 `np.where`, P76 grid layout reminders) |
| cross products of cylindrical unit vectors ($\mathbf e_\varphi\times\mathbf e_R=-\mathbf e_z$, $\mathbf e_\varphi\times\mathbf e_z=\mathbf e_R$) | C13 (D24) | **gloss** on P88 |
| flux of a 3-D point source through a sphere (4πr² u_r = Q) and 1/r potentials | C13 (D25) | **gloss** (ch05 P139 reminder) |
| substitution z − ξ = R cot α with $d(\cot\alpha)=-d\alpha/\sin^2\alpha$ and reversed limits | C14 (D26) | **gloss** on P106 substitution |
| collocation (demand the condition at N chosen points to get N equations) and the condition number `np.linalg.cond` | C14 (D27, code) | new numerics idea and Python tool |
| a function of x − x_s: $\partial/\partial\mathbf x_s=-\nabla$, and the multivariable chain rule through x_s(t), u_s(t) | C15 (D29) | P49/P91 chain rule reminders; the moving-argument trick is new — **gloss** |
| gradient of $\mathbf u\cdot\boldsymbol\xi/\lvert\boldsymbol\xi\rvert^3$ (product rule + ∇∣ξ∣⁻³ = −3ξ/∣ξ∣⁵) | C15 (D29) | **gloss** (P140 reminder) |
| surface integral on a sphere, $dA=a^2\sin\theta\,d\theta\,d\varphi$, e_r in Cartesian components, $\int_0^\pi\cos^2\theta\sin\theta\,d\theta=2/3$ | C15 (D30) | new maths tool (P88 spherical unit vectors reminder) |
| Green's first identity $\nabla\cdot(\phi\nabla\phi)=\lvert\nabla\phi\rvert^2+\phi\nabla^2\phi$ and turning a volume integral over the exterior of a sphere into a surface integral (normal into the sphere) | C15 (D31) | new vector identity (ch02 C14 Gauss reminder) |
| kinetic energy of the fluid as ½ × (added mass) × speed² | C15 (D31) | physics vocabulary — **gloss** |

Reminders only (already primed): partial derivative P25 · Taylor P26/P98 · definite integral P27 · net pressure force
P28 · lambda P29 · `solve_ivp` P31/P94 · trapezoid P37 · product rule P38 · sympy P40 · Euler's formula P45 · `np.where`
P46 · chain rule P49/P91 · `np.linalg.solve` P57 · orthonormal basis P65 · `arctan2` P70 · level sets and directional
derivative P75 · `meshgrid` and grid layout P76 · broadcasting P77 · contour/quiver/streamplot P78 · eigenvalues P80 ·
complex conjugate P81 · `quad` P87 · cylindrical/spherical unit vectors P88 · parametric curves P92 · substitution P106
· `brentq` P108 · dataclasses P111 · momentum flux P114 · Schwarz P121 · curl of a curl P122 · periodic trapezoid on a
loop P136 · Poisson and Green's function P139 · gradient of 1/distance P140 · FFT P142 · Gauss–Legendre P143 · improper
integral P144 · `animate` P16 · `slider_figure` P17 · `show_viz` P18.

## 4b. Derivations written out (parsed by tools: ID first, CORE id in a column, ★★★ for hard, explainer slugs backticked in the LAST column)
Written out because (a) the result belongs to an A item, or (b) the book never writes it out and the lesson needs it
(marked "(b) book never writes it out"). One small move per step; the book's skipped moves (analysis §2b) are filled in
and the notebook says so. Analysis §2b numbers in brackets (`a-Dnn`). 31 rows.

| ID | Result (Eq.) | CORE | Difficulty | Steps | Tools used | Traps | Shown in |
|---|---|---|---|---|---|---|---|
| D01 | the ideal-flow equations (6.1) from continuity (4.7) and Navier–Stokes (4.38) [a-D01] | C01 | ★★ | 8 | continuity (recap ch04 C02), hydrostatic split p = p_h + p′ (recap ch04), curl of a curl ∇²u = ∇(∇·u) − ∇×ω (P122), Kelvin (R01) | ∇·u = 0 needs constant ρ (M ≪ 1, no baroclinic field), not steadiness; p in (6.1) is measured from hydrostatic; μ∇²u = −μ∇×ω = 0 although μ ≠ 0 — this only works because no-slip is dropped (N02) | notebook |
| D02 | ω_z = −∇²ψ (6.4) ⇒ Laplace (6.5), point vortices as δ sources (6.6) [a-D03] | C02 | ★ | 5 | ψ (R03), vorticity in 2-D (R04), Laplacian (recap ch02 C10) | the minus sign: ∂/∂x(−∂ψ/∂x) − ∂/∂y(∂ψ/∂y); "irrotational except at points" is what turns Laplace into Laplace-plus-deltas | notebook |
| D03 | ln r is harmonic for r > 0 but the flux of ∇ ln r through any circle is 2π, so ψ = −(Γ/2π) ln r ↔ −Γδ and φ = (m/2π) ln r ↔ mδ — (b) book never writes it out (Exercises 6.1b, 6.3b) [a-D04] | C02 | ★★ | 9 | polar Laplacian (R08), 2-D divergence theorem and the 2-D delta (primer), logarithm rules (gloss) | ∇² ln r = 0 away from the origin does not mean zero total — the flux is independent of r; the sign of the vortex term follows from ψ = −(Γ/2π) ln r (counterclockwise Γ); m is the volume flux per depth, not per radian | notebook |
| D04 | equipotentials ⟂ streamlines: slopes multiply to −1, ∇φ·∇ψ = 0, ∣∇φ∣ = ∣∇ψ∣, ∇φ = e_z × ∇ψ [a-D05] | C02 | ★ | 5 | total differential (recap ch04 D04), dot product (recap ch02), P75 level sets | (v/u)(−u/v) = −1 fails at stagnation points (u = v = 0); ∇ψ = (−v, u) is u turned by +90°, not −90° | notebook |
| D05 | no through-flow (6.16): ∂φ/∂n = 0 ⇔ ∂ψ/∂s = 0, so a wall is a streamline and a streamline may be a wall [a-D06] | C03 | ★ | 5 | directional derivative (P75), rotating a vector by 90°, unit normal and tangent (P92) | ∂ψ/∂s = t·∇ψ = ±n·u — the sign depends on the orientation of t; a moving body needs n·u = n·U_s, not 0 | notebook · `superposition_sandbox` |
| D06 | the doublet as the limit of a source–sink pair (6.28) → (6.29) [a-D11] | C04 | ★★ | 9 | Taylor ln(1 + s) ≈ s (P26), logarithm rules (gloss), limit with a product held fixed (primer), dipole vector (N25) | ln√(1 + s) = ½ ln(1 + s), so the ½ cancels the 2 of 2εx/r²; take m → ∞ as ε → 0 with 2mε = ∣d∣ fixed (either limit alone gives 0 or ∞); d points from sink to source (−x here), so φ = −d·x/2πr²; the neglected term is O(ε²) | notebook · `superposition_sandbox` |
| D07 | the half-body (6.30)–(6.31): velocity, stagnation point x = −m/2πU, dividing streamline ψ = m/2, h(θ), h_max = m/2U by the limit and by mass balance [a-D13] | C05 | ★★ | 10 | superposition (C03), `arctan2` branch (P70, gloss), polar ψ (N18), mass conservation (recap ch04 C01) | the stagnation point is on the upstream axis where θ = π, so ψ there is m/2 (not m); v = 0 on y = 0 so only u = 0 must be solved; the θ-branch must put the cut inside the body; mass balance 2h_maxU = m checks the θ → 0 limit | notebook · `superposition_sandbox` |
| D08 | half-body surface speed and C_p: $C_p=-(2k\cos\theta+k^2)$, $k=\sin\theta/(\pi-\theta)$, zero where tan θ = −2(π − θ) (our root ≈ 113.2°) — (b) book never writes it out [a-D14] | C05 | ★★ | 8 | C_p (N28), body shape r = m(π − θ)/(2πU sin θ) (D07), `brentq` (P108) | use the speed *on the body*, not anywhere in the field; degrees vs radians in the root; C_p → 0 far downstream, so the infinite body still has zero net force | notebook |
| D09 | the circular cylinder (6.33)–(6.35): why d = −2πUa² e_x, u_r(a) = 0, surface speed 2U∣sin θ∣, C_p = 1 − 4 sin²θ [a-D15] | C06 | ★ | 6 | superposition (C03), doublet (C04), polar velocities (N18, N19) | d points upstream ("opposing the stream"); the surface speed is 2U at the shoulders, not U; C_p is symmetric under θ → −θ and θ → π − θ, the root of d'Alembert | notebook · `cylinder_circulation_lift` |
| D10 | cylinder with circulation: (6.36) → (6.37) → stagnation points (6.38), including the off-body pair [a-D16] | C07 | ★★ | 9 | vortex ψ (N08) with the clockwise sign, arcsin roots and quadratic roots (gloss), polar velocity (N19) | a clockwise vortex has ψ = +(Γ/2π) ln r (sign flip of (6.8)); ln(r/a) keeps ψ(a) = 0; both roots of sin θ = −Γ/4πaU lie below the x-axis; for Γ > 4πaU the two r-roots multiply to a² — keep the one with r > a | notebook · `cylinder_circulation_lift` |
| D11 | surface pressure (6.39) → lift integral → L = ρUΓ (6.40), and D = 0 [a-D17] | C07 | ★★ | 9 | Bernoulli (R05), net pressure force (P28), trigonometric integrals over a period (primer) | n·e_y = sin θ and the minus sign (the upper surface pushes down); expand the square fully — only the cross term 2UΓ sin θ/πa survives ∫ sin θ(…) dθ; D vanishes by fore–aft symmetry (∫ cos θ(…) = 0) | notebook · `cylinder_circulation_lift` |
| D12 | image rules: vortex images flip sign (ψ₂ odd in y), source images keep it (φ₂ even in y) [a-D18] | C08 | ★ | 5 | odd/even functions (gloss), superposition (C03), wall = streamline (N16) | the image lies outside the fluid, so the field equation in the fluid is unchanged; the vortex rule makes ψ₂(x, 0) = 0, the source rule makes ∂φ₂/∂y(x, 0) = 0 — two different conditions giving the same wall | notebook · `vortex_wall_images` |
| D13 | Example 6.1: ξ(t) = (h, Γt/4πh), φ, unsteady Bernoulli, and p(0, 0, t) as one fraction; sign change at t = 4πh²/Γ, maximum at t = 4√3πh²/Γ (our addition) [a-D20] | C08 | ★★ | 11 | wall image (R16, R17), unsteady Bernoulli (R32 / ch04 D26), derivative of arctan with moving Y (gloss), chain rule (P49) | the self-induced velocity of the ideal vortex is taken as zero (a modelling assumption); the constant is p∞/ρ; ∂φ/∂t is at a *fixed* point while the vortex moves; u(0, y) = 0 on the wall; the arctan branch must not jump between the vortex and its image | notebook · `vortex_wall_images` |
| D14 | Cauchy–Riemann (6.44) from direction-independent dw/dz, then dw/dz = u − iv (6.45) and Laplace for φ and ψ [a-D21 + a-D22] | C09 | ★★ | 8 | complex derivative along a direction (primer), 1/i = −i, mixed partials (P121), ψ and φ definitions (R03, N09) | along iy the difference quotient divides by i, giving ∂ψ/∂y − i ∂φ/∂y; dw/dz = u − iv (the conjugate), so the velocity vector is conj(dw/dz); CR plus Schwarz gives ∇²φ = ∇²ψ = 0 | notebook · `complex_potential_corners` |
| D15 | corner flow w = Az^n (6.46): walls at θ = 0 and θ = π/n, speed ∝ r^{n−1} = r^{(π−α)/α} at the corner [a-D23] | C09 | ★★ | 7 | complex powers and branch cuts (primer), polar form (N41) | ψ = Ar^n sin nθ vanishes on θ = 0 and θ = π/n; n = π/α; n < 1 (α > π) gives infinite speed at the corner, n > 1 a stagnation point; n = ½ is the plate and needs the cut along it | notebook · `complex_potential_corners` |
| D16 | force from the pressure integral (6.54)–(6.56): D = −∮p dy, L = ∮p dx [a-D25] | C10 | ★★ | 8 | CV momentum (R18), momentum flux (P114), outward normal of a counterclockwise contour (gloss) | the momentum flux vanishes on the body only because u·n = 0 there; n = (e_x dy − e_y dx)/ds is outward only for counterclockwise traversal (check at the rightmost point, dy > 0); F in (6.54) is on the fluid, D and L on the body | notebook · `blasius_kutta_contour` |
| D17 | Blasius's theorem (6.57)–(6.60): $D-iL=\frac{i\rho}{2}\oint_C(dw/dz)^2dz$ and why the contour can then be moved [a-D26] | C10 | ★★★ | 12 | complex conjugate (P81), dz* = dx − i dy, ∮ const dz* = 0, tangency (u + iv)dz* = (u − iv)dz, Cauchy's integral theorem (primer), sympy check cell | −i(dx − i dy) = −i dx − dy; the factor −i × (−½ρ) = +iρ/2; (6.59) holds only on the body, yet the final integrand is analytic, so the contour may be moved off the body provided no singularity of (dw/dz)² lies in between | notebook · `blasius_kutta_contour` |
| D18 | Kutta–Zhukhovsky for any body (6.61)–(6.62): D = 0, L = ρUΓ via the far-field Laurent series and the residue [a-D27] | C10 | ★★★ | 11 | Laurent series (primer), residue theorem (N63), series multiplication, sympy `series`/`residue` check cell | the far field has a real ln-coefficient (source) and an imaginary one (vortex); m = 0 for a closed body; the book prints the 1/z² coefficient as (Ud/π − Γ²/4π²) and an extra outer square — the correct coefficient is −(Ud/π + Γ²/4π²), and only the 1/z term iUΓ/π contributes; the clockwise Γ sign gives L = +ρUΓ | notebook · `blasius_kutta_contour` |
| D19 | conformal maps keep angles (6.63)–(6.64), except where dw/dz = 0 or ∞ [a-D28] | C11 | ★★ | 6 | multiplication of complex numbers adds arguments (primer), polar form (N41) | every small element is turned by the same angle arg f′ — the *local* statement; large shapes distort; at f′ = 0 angles are multiplied (z² doubles them at 0) | notebook · `conformal_joukowski` |
| D20 | Zhukhovsky: circle b → slit 2b cos θ, circle a → ellipse (6.66)–(6.67) with foci ±2b [a-D29] | C11 | ★★ | 8 | e^{iθ} + e^{−iθ} = 2cos θ, cos² + sin² = 1, foci c² = A² − B² (gloss) | the map is 2-to-1 (inside and outside the circle both cover the plane) — we use the outside; ζ = ±b are critical points (dz/dζ = 0), the slit ends, where angles are not kept | notebook · `conformal_joukowski` |
| D21 | inverse map (6.69) with the outside root, the chain-rule velocity and the flow round an elliptic cylinder with circulation (6.68) [a-D30] | C11 | ★★★ | 11 | quadratic formula in ℂ and complex square-root branches (primer), chain rule (P49), `joukowski_inverse` wrong-variant check cell | the two roots multiply to b² (one inside, one outside); numpy's principal √(z² − 4b²) returns the inside root for Re z < 0 — use ½[z + √(z − 2b)√(z + 2b)], whose cut is the slit; dζ/dz = 1/(1 − b²/ζ²) blows up at the slit ends | notebook · `conformal_joukowski` |
| D22 | five-point Laplace (6.70)–(6.72) and the average rule; truncation error (Δx²/12)ψ_xxxx [a-D31] | C12 | ★ | 6 | Taylor series (P26), half-point differences (R19) | the book's "first-order central differences" are first-*derivative* differences, second-order accurate; the average rule needs Δx = Δy; it implies no interior maximum (discrete mean-value property) | notebook · `laplace_relaxation` |
| D23 | the 4-point system (6.73) as Aψ = b and why Gauss–Seidel converges [a-D32] | C12 | ★★ | 9 | matrix form and diagonal dominance, Jacobi / Gauss–Seidel / SOR and the spectral radius (primer), `np.linalg.solve` (P57) | Gauss–Seidel uses the *latest* values in the same sweep (Jacobi does not); iterations grow like N² for Gauss–Seidel; stop on the field-equation residual, not only on the change between sweeps | notebook · `laplace_relaxation` |
| D24 | Stokes stream function (6.75) → field equation (6.77), which is not the Laplacian [a-D34] | C13 | ★★ | 8 | u = ∇χ × ∇ψ (R24), cross products of cylindrical unit vectors (gloss), ω_φ (R25), product rule (P38) | u_R = −(1/R)∂ψ/∂z carries the minus; the operator ∂_R((1/R)∂_Rψ) + (1/R)∂_zzψ equals (1/R)(E²ψ), not ∇²ψ — so no complex variables in axisymmetric flow | notebook |
| D25 | 3-D elements (6.86)–(6.88) and the sphere (6.89)–(6.91): why d = 2πa³U, surface speed (3/2)U sin θ, C_p = 1 − (9/4) sin²θ [a-D37 + a-D38] | C13 | ★★ | 10 | spherical velocities (N79), flux through a sphere (gloss), 3-D source–sink limit (as D06), spherical Laplace (R30) | φ = −Q/4πr (minus sign for a source); ψ = −(Q/4π)cos θ + const; ψ(a, θ) = 0 needs ½Ua² = d/(4πa); compare the cylinder (2U, −3) — the 3-D flow also escapes sideways | notebook |
| D26 | line sink (6.93) → closed form (6.94) → airship (6.95): closure Q = ak and the length from the axis stagnation points [a-D40 + a-D41] | C14 | ★★ | 10 | point-source ψ (N81), substitution with cot α (gloss on P106), superposition as an integral, `brentq` (P108) | a sink has dψ = +(k dξ/4π) cos α (sign opposite to a source); ξ = 0 ↔ α = θ, ξ = a ↔ α = α₁ and d(cot α) = −dα/sin²α reverses the limits; R/sin θ = r and R/sin α₁ = r₁; zero net source Q − ak = 0 closes the body | notebook · `axial_singularity_bodies` |
| D27 | the axial singularity method: ψ_mn per segment, ψ_m = 0 at N body points, the N × N system [a-D42] | C14 | ★ | 6 | segment ψ from D26, collocation and the condition number (primer), `np.linalg.solve` (P57) | a source segment has the opposite sign of the sink in (6.94); ψ = 0 on the axis and on the body is one surface; check Σk_nΔξ ≈ 0 (closure); blunt bodies make the matrix ill-conditioned as N grows | notebook · `axial_singularity_bodies` |
| D28 | moving-sphere potential (6.96)–(6.97): d(t) = 2πa³u_s and n·∇φ = n·u_s on the surface [a-D43] | C15 | ★ | 5 | coordinate-free potential (N85), Galilean frames (recap ch03 C05), kinematic condition (N16 with U_s) | in the sphere's frame the oncoming stream is −u_s, so the dipole follows +u_s; w is the z-velocity here, not the complex potential | notebook · `added_mass_sphere` |
| D29 | surface pressure on an arbitrarily moving sphere (6.98)–(6.105) [a-D44] | C15 | ★★★ | 14 | unsteady Bernoulli (R32), a function of x − x_s and the multivariable chain rule (gloss), gradient of u·ξ/∣ξ∣³ (gloss), ∣A∣² = A·A, sympy check cell | ∂/∂t at *fixed* x acts through x_s(t) and u_s(t); ∂φ/∂x_s = −∇φ; the book's middle expression of (6.104) prints −u_s/a³ — +u_s/a³ is right; expand ∣(3/2)(u·e)e − ½u∣² = (3/4)(u·e)² + ¼∣u∣² carefully | notebook · `added_mass_sphere` |
| D30 | added-mass force (6.107)–(6.108) and Newton's law (6.109) [a-D46] | C15 | ★★ | 9 | surface integral on a sphere (primer), e_r in Cartesian components (P88), ∫cos²θ sin θ dθ = 2/3 | the steady (speed-squared) part is fore–aft symmetric and integrates to zero (the book says "no drag" without integrating); the book's integrand keeps a stray dφ after the φ-integration; M = ½ρ(4/3)πa³; in general M is a tensor | notebook · `added_mass_sphere` |
| D31 | added mass by kinetic energy: KE = ½ρ∫∣∇φ∣² dV = ½MU² with M = 2πρa³/3 — (b) book never writes it out (Exercise 6.49) [a-D47] | C15 | ★★★ | 10 | Green's first identity and the exterior divergence theorem (primer), sphere potential on r = a, spherical surface integral (primer), sympy check cell | ∇·(φ∇φ) = ∣∇φ∣² needs ∇²φ = 0; the normal outward from the *fluid* points into the sphere, which fixes the sign; the far surface contributes nothing because φ∂φ/∂r ~ r⁻⁵ while the area grows like r² | notebook · `added_mass_sphere` |

Counts: ★ 8 (D02, D04, D05, D09, D12, D22, D27, D28) · ★★ 18 · ★★★ 5 (D17 and D18 in `blasius_kutta_contour`, D21 in
`conformal_joukowski`, D29 and D31 in `added_mass_sphere`; each has a sympy or wrong-variant check cell that re-runs the
construction: D17 checks (u + iv)dz* = (u − iv)dz on the cylinder surface symbolically and the ∮ of the Bernoulli
constant; D18 expands (U + iΓ/2πz − d/2πz²)² with sympy, prints the true 1/z² coefficient beside the printed one, and
takes the residue; D21 builds both roots of ζ² − zζ + b² = 0, checks their product b², and shows the principal-root
variant failing in the left half-plane; D29 re-runs (6.101)–(6.105) for a symbolic sphere motion and compares with the
direct ∂φ/∂t and the printed middle bracket; D31 integrates ½ρ∣∇φ∣² over r > a in sympy and compares with the surface
form). Written out under rule (b): D03, D08, D31 (3); under rule (a): the other 28. Analysis rows merged into one
written-out derivation: a-D21 + a-D22 → D14, a-D37 + a-D38 → D25, a-D40 + a-D41 → D26. No D row is shown in C01's,
C02's or C13's explainer because those A items have none (C02's is the backup B1).

## 4c. Derivations demoted to statements (first column is the A parent in bold, e.g. **C20** — never a bare ID or a D id: the parser reads a bare first-cell ID as an item and blanks its tier)
Results **given, not derived**, inside the named A block (a paragraph, the equation, a number). 14 rows.

| A parent | Analysis §2b item | Result stated (Eq.) | Stated in (B item) | Why not written out |
|---|---|---|---|---|
| **C02** | a-D02 (6.3) satisfies (6.2); ψ = const is a streamline | $0=d\psi=-v\,dx+u\,dy\Rightarrow(dy/dx)_\psi=v/u$ | R03 | SEEN: ch04 D04 derived the stream function and its streamlines; one line of mixed partials recalled |
| **C03** | a-D07 Bernoulli everywhere (6.18) and the unsteady form | $p+\tfrac12\rho\lvert\nabla\phi\rvert^2=$ const; $\partial\phi/\partial t+\tfrac12\lvert\nabla\phi\rvert^2+p/\rho=$ const | R05 | SEEN: ch04 D25 and D26 derived both; recalled with ∣∇ψ∣ = ∣∇φ∣ from D04 |
| **C02** | a-D08 polar forms (6.19)–(6.23) | $u_r=\partial\phi/\partial r=\frac1r\partial\psi/\partial\theta$, $u_\theta=\frac1r\partial\phi/\partial\theta=-\partial\psi/\partial r$ | N18, N19 (R06–R09) | the book gives them "for quick reference"; App. B operators are ch04's `core.curvilinear`; a sympy chain-rule cell checks them |
| **C04** | a-D09 quadratic harmonic polynomials and the 45° rotation (6.24)–(6.27) | $\nabla^2(ax^2+bxy+cy^2)=2(a+c)=0$ ⇒ families xy and x² − y² | N20–N23, R10 | two lines of algebra; generated and checked by `harmonic_polynomials`, and derived in general as Re/Im of Az^n in D15 |
| **C04** | a-D10 source and vortex velocities | $u_r=m/2\pi r$, $u_\theta=\Gamma/2\pi r$ | N24, R11 | one differentiation of ln√(x² + y²); the vortex half is SEEN (ch03, ch05) |
| **C04** | a-D12 doublet stream function and its circular streamlines (Exercise 6.9) | $\psi=-\frac{\lvert\mathbf d\rvert}{2\pi}\frac{y}{r^2}$; $x^2+(y+\lvert\mathbf d\rvert/4\pi c)^2=(\lvert\mathbf d\rvert/4\pi c)^2$ | N25 (and N47) | follows in one line from Im of (6.49) in C09; completing the square stated with the figure |
| **C08** | a-D19 streamline equation of the two-source flow (6.41) | $x^2-y^2-2xy\cot(2\pi\psi/m)=a^2$ | N40 | the tangent addition formula in one line; the curve is drawn on the ψ contours, which is the check that matters |
| **C09** | a-D24 complex potentials of the elements (6.47)–(6.53) | ln(z − z′) = ln r′ + iθ′; 1/(z − z′) = (x − iy)/r² | N45–N51 | each is a one-line real/imaginary split of a result already derived in §§6.2–6.3 |
| **C12** | a-D33 Example 6.2 boundary values | ψ linear across inlet and outlet, 0 on the lower wall, Q on the upper | N74 | u = ∂ψ/∂y constant ⇒ ψ linear — one sentence; the re-entrant corner (α = 3π/2 in (6.46)) is explained with D15's result |
| **C13** | a-D35 flux between stream surfaces (6.78) | $dQ=2\pi R(-u_Rdz+u_zdR)=2\pi\,d\psi$ | N77 | the plane-flow flux argument of ch04 D04 with a ring of circumference 2πR; checked by quadrature |
| **C13** | a-D36 axisymmetric φ and spherical forms (6.79)–(6.85) | (6.83) $u_r=\frac1{r^2\sin\theta}\partial_\theta\psi$, $u_\theta=-\frac1{r\sin\theta}\partial_r\psi$ | N79, R26–R30 | App. B operators applied to u = ∇φ and u = −(1/(r sin θ)) e_φ × ∇ψ; stated with a sympy check |
| **C13** | a-D39 coordinate-free sphere potential (6.92) | $\phi=(\mathbf U-\frac{\mathbf d}{4\pi\lvert\mathbf x\rvert^3})\cdot\mathbf x$ | N85 | one substitution cos θ/r² = e_z·x/∣x∣³; the vector form is the input of D28 |
| **C15** | a-D45 steady limit (6.106) = (6.91) | $\frac94\cos^2\theta_s-\frac54=1-\frac94\sin^2\theta_s$ | N102 | cos² = 1 − sin² and θ_s = π − θ; one line, asserted in code |
| **C15** | a-D48 cylinder added mass ρπa² per depth (Exercise 6.45) | $F=-\rho\pi a^2\,dU_c/dt$ — the whole displaced mass in 2-D vs half in 3-D | N104 | an exercise; the moves repeat D29–D30 in polar form, so the result is stated with a numeric check (`ch06.cylinder_added_mass`, force and energy routes) |

## 5. Interactive explainers (5–10 + backup)
Nine explainers on the A items where manipulation or motion teaches most; one backup. Each has the required live
**Explain** tab ("Explanation & interpretation", numbered sections computing every number on screen with the reader's
settings, ending in "Reading the current setting"), a synced **Code** tab, a **Derivation** tab for its D ids, a
4–8-step walkthrough, ≥ 3 check questions and ≥ 2 more depth features. Physics mirrors scalar-callable `fluidpy`
functions (parity rows). Colours (consistent across explainers, notebook and derivation terms): uniform stream blue,
sources teal, sinks rose, vortices amber, doublets purple, the body/dividing streamline black; pressure above p∞ orange,
suction blue; steady pressure part teal, acceleration part orange; lift amber, drag muted. Not duplicated from earlier
chapters, linked instead: ch03 `galilean_frames_cylinder` (from E2: the cylinder in two frames), ch04
`stream_function_spacing` (from E1: spacing = speed), ch05 `point_vortex_lab` (from E3: vortex dynamics with images;
E3 adds pressure), ch04 `which_bernoulli` (from E3: unsteady Bernoulli).

### E1 · superposition_sandbox
- **A:** C03, C04, C05 (also shows N16 (6.16), N17 (6.17), R05 (6.18), N24, N25 (6.28), N27 (6.30), N28 (6.32), N29;
  closed bodies preview C06; the Rankine oval of Exercise 6.19 as a preset)
- **Confusion removed:** "each element solves Laplace — but where is the body, and who imposed the boundary
  condition?" The body is simply the streamline through the stagnation point of the *sum*; a net source leaves it
  open (half-body), a zero net source closes it (oval, cylinder).
- **Why interactive:** the reader toggles elements and slides strengths and watches the stagnation points and the
  dividing streamline appear, move and close; the element velocities at a clicked point add up as bars. A static figure
  shows one finished sum, not the act of building it.
- **Stage:** (1) the flow: ψ contours at equal Δψ over a speed (or C_p) map, stagnation points marked, dividing
  streamline bold black, tracer particles; (2) term bars at a clicked probe: u and v from each element and the total
  (with the inspector's arithmetic); (3) (hidePortrait) C_p along the dividing streamline vs arc length.
- **Controls:** element chips (uniform U, source m, sink −m, doublet ∣d∣, vortex Γ) · strength sliders · source–sink
  separation 2ε · probe (click) · view: speed / C_p.
- **Presets:** stream only · half-body (U = 1, m = 2π) · Rankine oval · source–sink → doublet (ε shrinking) · cylinder ·
  lifting cylinder.
- **Equations shown (live):** (6.7) $\psi=-Vx+Uy$; (6.15) $\phi=\frac{m}{2\pi}\ln r$; (6.8) $\psi=-\frac{\Gamma}{2\pi}\ln r$;
  (6.29) $\phi=-\mathbf d\cdot\mathbf x/2\pi r^2$; (6.31) $\psi=Ur\sin\theta+\frac{m}{2\pi}\theta$; (6.32)
  $C_p=1-\lvert\mathbf u\rvert^2/U^2$; (6.16) $\partial\psi/\partial s=0$ on the body.
- **Mirrors:** `core.potential.Flow`, `Uniform`, `Source`, `Vortex`, `Doublet`, `half_body`, `Flow.stagnation_points`,
  `ch06.half_body_numbers`, `ch06.rankine_oval`, `ch06.doublet_limit_error`.
- **Derivations:** D05 (wall = streamline: the step on ∂ψ/∂s highlights the dividing streamline), D06 (doublet limit:
  steps shrink ε), D07 (half-body: steps mark a, ψ = m/2, h_max).
- **Depth features:** Explain tab (each element's ψ at the probe → the sum → stagnation point from u = 0 → a, h_max →
  C_p at the probe → reading the current setting: open/closed body), synced Code tab, + linked views (3), presets, term
  bars, inspector, status ("🟢 closed body: Σm = 0" / "🔵 open body: net source m ⇒ extends downstream" / "⚪ no body").
- **Follows reference:** `fid_formula_lab.html` (a sum whose terms are clickable bars with a total) with
  `angular_frequency_explorer_1.html`'s presets and "Right now" notes.
- **Aha:** add solutions, then pick the streamline through the stagnation point and call it a wall — the body appears
  by itself; it closes only when the sources and sinks cancel.

### E2 · cylinder_circulation_lift
- **A:** C06, C07 (also shows R12 (6.33), R13 (6.34), N30 (6.35), N31 (qualitative real band), N32–N36 (6.36)–(6.39),
  N38 uniqueness; links to ch03 `galilean_frames_cylinder`)
- **Confusion removed:** "why does circulation give lift but no drag — and which way does Γ turn?" The circulation
  speeds the top and slows the bottom; the stagnation points slide down; the top–bottom pressure difference integrates
  to exactly ρUΓ while the fore–aft symmetry that kills drag never breaks. The clockwise-Γ convention is on screen.
- **Why interactive:** dragging Γ moves the stagnation points continuously round the body, makes them meet at 4πaU and
  leave the body — a sequence a static figure (Fig. 6.12's four snapshots) only samples; the force bars move with the
  pressure arrows.
- **Stage:** (1) cylinder streamlines with stagnation points and surface pressure arrows −p n dl (orange/blue); (2)
  surface C_p(θ) with the Γ = 0 ghost and the qualitative real band (toggle); (3) force bars: drag (≈ 0), lift from the
  pressure integral with the measured ◇ at ρUΓ.
- **Controls:** Γ/(4πaU) slider (−2 … 2) · U · a · convention chips (book clockwise Γ / project counterclockwise) · real
  band toggle · transport sweeping Γ.
- **Presets:** Γ = 0 (d'Alembert) · Γ = 2πaU · Γ = 4πaU (merging) · Γ = 6πaU (off-body stagnation point) · negative Γ
  (lift down) · moving cylinder (fluid frame, link to ch03 E3).
- **Equations shown (live):** (6.33) $\psi=U(r-a^2/r)\sin\theta$; (6.35) $C_p=1-4\sin^2\theta$; (6.36); (6.37)
  $u_\theta(a,\theta)=-2U\sin\theta-\Gamma/2\pi a$; (6.38) $\sin\theta=-\Gamma/4\pi aU$; (6.39); (6.40) $L=\rho U\Gamma$.
- **Mirrors:** `core.potential.cylinder`, `ch06.cylinder_stagnation_points`, `ch06.cylinder_surface_cp`,
  `ch06.cylinder_surface_pressure`, `ch06.surface_pressure_force`, `ch06.lift_per_span`,
  `ch06.circulation_family_check`.
- **Derivations:** D09 (cylinder: step "choose d" sets Γ = 0), D10 (stagnation points: steps move Γ through 4πaU), D11
  (lift integral: steps light the surviving cross term on the C_p curve).
- **Depth features:** Explain tab (surface speed at the clicked angle → C_p → pressure → the integrand p sin θ → L by
  quadrature vs ρUΓ → stagnation angles → reading the current setting: below/at/above 4πaU), synced Code tab, + linked
  views (3), presets, status ("two surface stagnation points" / "merged at the bottom" / "free stagnation point in the
  flow"), term bars (D, L), inspector (click a surface point: its p arithmetic), transport.
- **Follows reference:** `amplitude_phase_second_order_II_3.html` (system + graphs linked by one state; a numbered
  derivation with live numbers in the explanation).
- **Aha:** circulation adds speed on one side and removes it on the other — the pressure difference is exactly ρUΓ,
  and the front–back symmetry that makes the drag zero survives every Γ.

### E3 · vortex_wall_images
- **A:** C08 (also shows N39 image rules, N40 (6.41), R16, R17 (drift Γ/4πh), R32 unsteady Bernoulli; links to ch05
  `point_vortex_lab`, ch04 `which_bernoulli`)
- **Confusion removed:** "the image is a mathematical trick — what does the real wall feel?" The wall feels the vortex
  as a moving pressure signature: suction while it is near, then over-pressure as it leaves, the positive part made by
  the unsteady term ρ∂φ/∂t; and vortex images flip sign while source images do not.
- **Why interactive:** the result is a time series tied to a moving vortex; scrubbing time links the vortex's position
  to the wall-pressure trace and to the pressure along the whole wall, and switching vortex ↔ source shows the other
  image rule on the same stage.
- **Stage:** (1) the wall (x = 0) with the vortex, its faded image, streamlines and wall points coloured by p − p∞; (2)
  p(0, 0, t) trace with the moving time marker, the zero crossing at 4πh²/Γ and the maximum at 4√3πh²/Γ; (3)
  (hidePortrait) p − p∞ along the wall at the current time.
- **Controls:** Γ · h · t (transport) · mode chips (vortex by a wall / source by a wall (6.41)) · wall probe position.
- **Presets:** t = 0 (strongest suction) · t = 4πh²/Γ (crossing) · t = 4√3πh²/Γ (largest over-pressure) · source near a
  wall · two sources.
- **Equations shown (live):** image ψ of Example 6.1; ξ(t) = (h, Γt/4πh); the p(0, 0, t) fraction; unsteady Bernoulli
  (4.75) $\partial\phi/\partial t+\tfrac12\lvert\nabla\phi\rvert^2+p/\rho=p_\infty/\rho$; (6.41) streamline equation.
- **Mirrors:** `ch06.example_6_1`, `core.potential.mirror`, `ch06.two_sources`, `ch06.two_source_streamline`,
  `core.biot_savart.point_vortex_evolve(boundary="wall")`.
- **Derivations:** D12 (image rules: steps flip the image sign and show the wall condition failing), D13 (Example 6.1:
  steps build ξ(t), φ, ∂φ/∂t at the origin, the fraction).
- **Depth features:** Explain tab (drift speed → vortex position → v at the wall probe → ∂φ/∂t → p − p∞ as the sum of
  −ρ∂φ/∂t and −½ρv² → reading the current setting: approaching / passing / leaving), synced Code tab, + transport,
  linked views (3), presets, term bars (unsteady vs speed part of the pressure), modes (vortex / source), inspector.
- **Follows reference:** `forced_damped_vibrations.html` (one time slider driving the system and its x(t) graph, with
  the explanation panel computing the value at the current time).
- **Aha:** a passing eddy is felt at the wall as suction and then over-pressure — the over-pressure comes from the
  flow changing in time, not from its speed.

### E4 · complex_potential_corners
- **A:** C09 (also shows N41 (6.43), N42 (6.44), N43 (6.45), N44 (6.46), N45–N51 (6.47)–(6.53))
- **Confusion removed:** "why is the velocity u − iv and not u + iv, and why does any analytic function give a flow?"
  The derivative of w taken along x and along iy agrees only for analytic w — that agreement *is* Cauchy–Riemann — and
  its conjugate is the velocity; the exponent n alone decides whether a corner is calm or violent.
- **Why interactive:** the reader clicks a point and sees the two difference quotients (along x, along iy) converge to
  the same dw/dz, with the velocity arrow drawn as its mirror image; dragging n bends the walls from a 45° wedge to a
  flat plate and the speed near the corner flips from zero to infinite — continuous change a static figure cannot give.
- **Stage:** (1) z-plane flow net (ψ solid, φ dashed) with the walls, the probe, the dw/dz arrow and the velocity arrow
  (its conjugate); (2) the two difference quotients converging as h shrinks, with Cauchy–Riemann bars ∂φ/∂x vs ∂ψ/∂y,
  ∂φ/∂y vs −∂ψ/∂x; (3) (hidePortrait) log–log ∣dw/dz∣ vs r near the corner with slope n − 1.
- **Controls:** w family chips (Az^n, source, vortex, doublet, cylinder, the non-analytic control z*) · n slider (½ … 4,
  α = π/n) · A · probe (click) · h (difference step).
- **Presets:** n = 2 (stagnation corner, (6.24)) · n = 1 (uniform) · n = 2/3 (270° corner of Example 6.2) · n = ½ (flat
  plate) · z* (Cauchy–Riemann fails: "not a flow").
- **Equations shown (live):** (6.42) $w=\phi+i\psi$; (6.43) $z=re^{i\theta}$; (6.44) Cauchy–Riemann; (6.45)
  $dw/dz=u-iv$; (6.46) $w=Az^n$, $dw/dz=(A\pi/\alpha)z^{(\pi-\alpha)/\alpha}$; (6.47)–(6.49).
- **Mirrors:** `core.potential.Corner`, `core.potential.Flow.w`, `Flow.dwdz`, `ch06.cauchy_riemann_residual`,
  `ch06.harmonic_polynomials` (degree table).
- **Derivations:** D14 (Cauchy–Riemann and u − iv: steps shrink the difference step along each direction), D15 (corner
  flow: steps set n and mark the walls).
- **Depth features:** Explain tab (w at the probe → dw/dz by both quotients → u and v → Cauchy–Riemann residuals →
  corner angle α = π/n and the speed exponent → reading the current setting: stagnation / uniform / infinite speed),
  synced Code tab, + linked views (3), inspector (the dw/dz arithmetic), term bars (Cauchy–Riemann pairs), presets,
  status ("✅ analytic: a flow" / "⛔ z* is not analytic").
- **Follows reference:** `amplitude_phase_second_order_II_3.html` (a complex number evaluated live, its conjugate and
  the geometric reading in linked views).
- **Aha:** one complex function carries both families of lines; its derivative is the velocity with v flipped, and the
  corner angle alone decides whether the flow at the corner stops or races.

### E5 · blasius_kutta_contour
- **A:** C10 (also shows R18 (6.54), N53–N63 (6.55)–(6.61); the ellipse body from C11)
- **Confusion removed:** "lift must depend on the shape — how can L = ρUΓ hold for every body?" The Blasius integral
  does not change when the contour is squeezed or stretched round the body, and when the far-field series is squared
  only the U × Γ/z cross term has a residue.
- **Why interactive:** the reader drags the contour (radius, offset, shape) and watches the integral stay pinned while
  the integrand changes completely; the Laurent bars show which term survives. A static proof cannot show invariance
  under deformation.
- **Stage:** (1) the body (cylinder, Zhukhovsky ellipse, tilted ellipse, Rankine oval + vortex) with streamlines and the
  draggable contour; (2) Laurent term bars: the contribution of each power of z to ∮(dw/dz)² dz (only the 1/z bar is
  nonzero); (3) D and L vs contour radius (flat lines at 0 and ρUΓ) with the on-body pressure-integral route as ◇.
- **Controls:** body mode chips · Γ · U · contour radius R (drag) · quadrature points n (convergence).
- **Presets:** cylinder, Γ = 0 (d'Alembert) · cylinder with Γ (Exercise 6.21) · ellipse with the same Γ (same lift) ·
  contour cutting the body (⚠️ invalid) · few quadrature points (spectral convergence visible).
- **Equations shown (live):** (6.56) $D=-\oint p\,dy$, $L=\oint p\,dx$; (6.57) $D-iL=-i\oint p\,dz^*$; (6.60)
  $D-iL=\frac{i\rho}{2}\oint(dw/dz)^2dz$; (6.61) with the corrected coefficient; (6.62) $L=\rho U\Gamma$.
- **Mirrors:** `core.potential.blasius_force`, `core.potential.laurent_coefficients`, `ch06.contour_force`,
  `ch06.cv_force_on_body`, `ch06.elliptic_cylinder_flow`.
- **Derivations:** D16 (pressure integral: steps orient the contour and draw n), D17 (Blasius: steps move from the body
  to a larger contour), D18 (Kutta–Zhukhovsky: steps light the Laurent bars and the residue).
- **Depth features:** Explain tab (Laurent coefficients c₀ = U, c₋₁ = iΓ/2π → the square's 1/z coefficient iUΓ/π →
  residue × 2πi × iρ/2 → D, L → comparison with the pressure route → reading the current setting), synced Code tab, +
  linked views (3), term bars (Laurent contributions), modes (bodies), presets, status ("✅ contour encloses the body
  and no other singularity" / "⚠️ contour crosses the body").
- **Follows reference:** `fid_formula_lab.html` (term-by-term bars and a total, clickable terms) with
  `amplitude_phase_second_order_II_3.html`'s numbered live explanation.
- **Aha:** stretch the contour any way you like around the body — the answer does not move; square the far-field
  series and only U times Γ/z survives, so every body with circulation Γ gets L = ρUΓ and no drag.

### E6 · conformal_joukowski
- **A:** C11 (also shows N64–N70 (6.63)–(6.69), C09's w; N66's grid maps)
- **Confusion removed:** "how does a map carry a flow, why are angles kept, and why does the inverse need the root
  outside the circle?" Small crosses keep their angles because every small element is turned by the same arg f′; the
  flow net of the circle becomes the flow net of the ellipse; the wrong square-root branch drops the flow inside the
  body.
- **Why interactive:** the reader drags a small cross in the ζ-plane and watches its image in the z-plane stay a right
  angle (and fail at ζ = ±b); slides a/b from a circle to a plate; toggles the numpy branch and sees half the streamlines
  jump inside. Two linked planes on one state is the natural interactive form.
- **Stage:** (1) ζ-plane: circle of radius a, flow (6.68) streamlines with Γ, the draggable element pair, the critical
  points ±b; (2) z-plane: the mapped body (slit / ellipse), mapped streamlines, the mapped element pair with the angle
  readouts; (3) (hidePortrait) ∣dz/dζ∣ around the circle (zero at ±b when a = b).
- **Controls:** a/b · Γ · stream angle · element probe (drag) · branch chips (outside root / numpy principal) · map mode
  (Zhukhovsky / z² / e^z for angle preservation and grid images).
- **Presets:** a = b (flat plate) · a = 1.2b (ellipse, numbers of 2a) · a ≫ b (almost a circle) · Γ = 0 · wrong branch ·
  z² at 0 (angles doubled).
- **Equations shown (live):** (6.63) $\delta w=\frac{dw}{dz}\delta z$; (6.64) α = β; (6.65) $z=\zeta+b^2/\zeta$; (6.66);
  (6.67) the ellipse; (6.68) circle flow with Γ; (6.69) the inverse and $u-iv=\frac{dw}{d\zeta}\frac{d\zeta}{dz}$.
- **Mirrors:** `core.conformal.joukowski`, `joukowski_derivative`, `joukowski_inverse`, `mapped_flow`,
  `angle_preservation`, `grid_image`, `ch06.joukowski_ellipse`, `ch06.elliptic_cylinder_flow`.
- **Derivations:** D19 (angles: steps rotate the element pair), D20 (circle → slit/ellipse: steps slide a), D21 (the
  inverse: steps show both roots and switch the branch).
- **Depth features:** Explain tab (f′ at the element → scale and turn → both angles → semi-axes and foci → the inverse
  root at the probe with ∣ζ∣ → velocity by the chain rule → reading the current setting), synced Code tab, + linked views
  (3), presets, inspector (the element's mapping arithmetic), modes (three maps), status ("✅ outside root" / "⚠️
  principal root lands inside the circle").
- **Follows reference:** `angular_frequency_explorer_1.html` (linked views on one state, modes, highlighted table of
  special cases).
- **Aha:** solve the easy problem, a circle, and let an analytic map carry it — small shapes keep their angles, so the
  flow net stays a flow net; just take the root outside the circle.

### E7 · laplace_relaxation
- **A:** C12 (also shows R19–R21 (6.70)–(6.71), N72 (6.73), N73 Gauss–Seidel, N74 Example 6.2; N71 named)
- **Confusion removed:** "how does a grid of numbers 'solve' a PDE?" Laplace on a grid means each value is the average
  of its four neighbours; relaxation keeps enforcing that until nothing changes; Gauss–Seidel uses the newest values,
  SOR over-relaxes, and the residual (not the last change) says when to stop.
- **Why interactive:** the algorithm is a process: stepping one node at a time, then whole sweeps, with the residual
  curve growing to the left is how the method becomes visible (the `pixels_as_parameters` pattern).
- **Stage:** (1) the grid of Example 6.2's contraction (or the 4-point problem) with ψ as colour and numbers on coarse
  grids, the node being updated highlighted, streamline contours; (2) residual and max change vs sweep on log axes, with
  Jacobi / Gauss–Seidel / SOR ghosts; (3) (hidePortrait) the 4 × 4 matrix system with the current iterate.
- **Controls:** method chips (Jacobi / Gauss–Seidel / SOR) · ω (SOR) · grid refinement (×1, ×2, ×4) · initial guess
  (zero / linear) · transport (node step, sweep, play).
- **Presets:** 4-point system · book grid · refined ×4 · SOR at ω_opt · harmonic test ψ = xy (exact on the grid).
- **Equations shown (live):** (6.70), (6.71) second differences; (6.72) the average rule; (6.73) the four equations;
  the Gauss–Seidel update.
- **Mirrors:** `core.laplace_solvers.solve_laplace`, `core.laplace_solvers.gauss_seidel_sweep` (and `jacobi_sweep`,
  `sor_sweep`), `core.laplace_solvers.laplacian_5pt`, `ch06.four_point_system`, `ch06.example_6_2`.
- **Derivations:** D22 (five-point rule: steps build the stencil on a clicked node), D23 (the 4-point system and
  convergence: steps run sweeps and show the residual falling).
- **Depth features:** Explain tab (the average at the highlighted node with its four neighbours → change → residual
  norm → sweeps so far and predicted sweeps from the spectral radius → reading the current setting), synced Code tab, +
  transport, linked views (3), presets, inspector (click a node: its average), status ("⏳ iterating" / "✅ converged:
  residual < 10⁻⁸"), modes (methods).
- **Follows reference:** `pixels_as_parameters.html` (an algorithm in slow motion, step back/forward, a clickable
  loss curve) with `forced_damped_vibrations.html`'s explanation panel.
- **Aha:** Laplace's equation just says "be the average of your neighbours" — relaxation keeps enforcing it until
  nothing changes, and a smarter sweep gets there many times faster.

### E8 · axial_singularity_bodies
- **A:** C14 (also shows N86–N89 (6.93)–(6.95), N90 source panels (our extension), C13's sphere as a target and N84)
- **Confusion removed:** "given a body shape, how do we find the flow?" Put unknown sources on the axis (or on the
  surface), demand ψ = 0 (or u·n = 0) at N points, and one linear solve does it; the strengths must add to zero or the
  body will not close, and blunt bodies make the system ill-conditioned.
- **Why interactive:** the reader raises N and watches the computed ψ = 0 contour snap onto the target while the k_n
  bars converge to the exact distribution (and the condition number climbs) — convergence is a sequence, not a picture.
- **Stage:** (1) meridian plane: the target body (dashed), the computed body ψ = 0 (solid), axial segments coloured by
  k_n (sources teal, sinks rose); (2) k_n bars along the axis vs the exact distribution (Rankine oval: point source +
  sink; airship: uniform line sink); (3) body error and condition number vs N (log). Panel mode: circle/ellipse with N
  panels, C_p at midpoints vs 1 − 4 sin²θ, error vs N.
- **Controls:** mode chips (axial sources, book §6.8 / 2-D source panels, extension) · target (airship, Rankine oval,
  sphere, ellipsoid) · N · U · fineness ratio.
- **Presets:** airship (6.95) · Rankine oval (exact recovery) · sphere (point doublet: ill-conditioned, fenced) ·
  panels N = 8 · panels N = 64.
- **Equations shown (live):** (6.87) point source; (6.93) line sink; (6.94) $\psi_{\text{sink}}=\frac{k}{4\pi}(r-r_1)$;
  (6.95) airship ψ; the collocation system ψ_m = 0; closure Σk_nΔξ = 0.
- **Mirrors:** `core.panels.axial_singularity_solve`, `core.panels.axial_singularity_velocity`, `ch06.airship`,
  `ch06.line_sink_stream_function`, `core.panels.source_panels`.
- **Derivations:** D26 (line sink and airship: steps add source, sink, stream), D27 (the axial system: steps place N
  points and assemble the matrix).
- **Depth features:** Explain tab (one row of the influence matrix at a clicked body point → the solve → Σk_nΔξ →
  body error → condition number → reading the current setting: converged / ill-conditioned), synced Code tab, + linked
  views (3), modes (axial / panels), presets, inspector (ψ_m as a sum of segment terms), status ("✅ closed: Σ ≈ 0" /
  "⚠️ cond > 10¹⁰").
- **Follows reference:** `stride_padding_playground.html` (presets that are the classic settings, the formula with
  numbers plugged in, badges explaining the result) with `angular_frequency_explorer_1.html`'s modes.
- **Aha:** flip the problem — hide unknown sources inside the body, ask for a streamline at N points, and one linear
  solve draws the airship; the sources must cancel or the body never closes.

### E9 · added_mass_sphere
- **A:** C15 (also shows N92–N104 (6.96)–(6.108), R32 (6.99), N84 (6.91) as the steady limit, N105)
- **Confusion removed:** "an ideal fluid gives no drag, so why does it resist being accelerated?" The speed-dependent
  pressure is fore–aft symmetric and cancels; the acceleration-dependent pressure is high in front and low behind and
  adds up to exactly half the displaced mass times the acceleration.
- **Why interactive:** the reader sets speed and acceleration independently (even at an angle to each other) and sees
  the two pressure parts on the sphere and their separate forces — one always zero, one always −M du_s/dt; the ball /
  bubble mode integrates (6.109) in time.
- **Stage:** (1) sphere section (optional 3-D view) with surface pressure colour and arrows, fluid-frame streamlines,
  arrows for u_s and du_s/dt; (2) surface pressure vs θ_s: steady part (teal, symmetric) + acceleration part (orange,
  antisymmetric) = total; (3) force vs time with the −M du_s/dt ghost; in ball/bubble mode, velocity vs time with and
  without added mass.
- **Controls:** u_s · du_s/dt · angle between them · a · ρ · mode chips (prescribed motion / oscillating sphere / ball /
  bubble) · transport.
- **Presets:** steady motion (no force) · starting from rest · oscillating sphere (force in phase with acceleration) ·
  air bubble in water (2g) · steel ball in water.
- **Equations shown (live):** (6.97) $\phi=-\frac{a^3}{2\lvert\boldsymbol\xi\rvert^3}\mathbf u_s\cdot\boldsymbol\xi$;
  (6.104) $\mathbf u_a=\tfrac32(\mathbf u_s\cdot\mathbf e_\xi)\mathbf e_\xi-\tfrac12\mathbf u_s$; (6.105); (6.106);
  (6.108) $M=\frac{2\pi a^3\rho}{3}$; (6.109) $\mathbf F_E=(m+\frac{2\pi}{3}\rho a^3)\frac{d\mathbf u_s}{dt}$.
- **Mirrors:** `core.potential.moving_sphere_surface_pressure`, `core.potential.sphere_force_quadrature`,
  `core.potential.added_mass_sphere`, `ch06.added_mass_by_energy`, `ch06.sphere_motion`.
- **Derivations:** D28 (moving-sphere potential), D29 (surface pressure, ★★★: steps light the steady and acceleration
  parts), D30 (the force integral: steps integrate the two parts), D31 (kinetic energy route, ★★★: the energy bar
  matches ½MU²).
- **Depth features:** Explain tab (surface velocity at the clicked point → steady pressure → acceleration pressure →
  the two force integrals → M and the displaced mass → the ball/bubble acceleration → reading the current setting),
  synced Code tab, + linked views (3), transport, term bars (steady vs acceleration force; energy route), presets, modes,
  3-D view (optional), status ("⚖️ steady: no force" / "🚀 accelerating: fluid pushes back with M").
- **Follows reference:** `forced_damped_vibrations.html` (the oscillating sphere as a system + response graph with the
  full explanation panel).
- **Aha:** speed pushes symmetrically and cancels; acceleration pushes back with half the displaced fluid's mass — so
  a bubble starts upward at 2g, not at infinity.

### B1 · flow_net_sources_vortices (backup)
- **A:** C02 (also shows N05–N14 (6.5)–(6.15), N09 orthogonality, R08 polar Laplacian); derivations D02, D03, D04 —
  built only if a chosen explainer fails review.
- **Confusion removed:** "ln r solves Laplace — so where does the vortex or source strength come from?" A draggable
  circle measures ∮∇ψ·n ds (= −Γ) and ∮∇φ·n ds (= m) and gets the same number for every radius, and 0 when it does not
  enclose the centre; ω = −∇²ψ as a heatmap is zero everywhere except the core.
- **Stage:** flow net (ψ ⟂ φ) with the draggable loop · flux vs radius (flat) · ω = −∇²ψ heatmap (point vortex vs
  Rankine core). **Mirrors:** `ch06.delta_flux_check`, `ch06.vorticity_from_psi`, `ch06.orthogonality_check`.
- **Depth features:** explain, code, linked views, inspector (flux arithmetic), presets (source, vortex, both, Rankine).
- **Follows reference:** ch02 `gauss_flux_box` pattern with `forced_damped_vibrations.html`'s explanation panel.
- **Aha:** a vortex or a source is invisible to Laplace's equation everywhere except at one point — and a loop anywhere
  around it measures the same strength.

## 6. Python animations and interactive figures (A ID → what, why, player/figure kind)
Animations (`animate` + `show_animation`, ≤ 120 frames, dpi 80, FAST halves the frames):
- **C04** — a source–sink pair squeezed together (ε from 0.5 to 0.02 m, 2mε fixed): the streamlines morph into the
  doublet's tangent circles, with the error printed; why: the limit *is* a process; `player="frames"` (stop at each ε).
- **C08** — Example 6.1: the vortex drifting up the wall beside its image, wall points coloured by p − p∞, and the
  p(0, 0, t) trace drawing itself; why: the pressure signal is a time story; `player="video"`.
- **C12** — Gauss–Seidel sweeps on Example 6.2's contraction: ψ relaxing from zero and the residual falling; why: the
  iteration is the method; `player="frames"` (step sweeps).
- **C15** — an oscillating sphere: surface pressure (steady teal + acceleration orange) and the force arrow in phase
  with the acceleration; why: added mass is felt only while the speed changes; `player="video"`.

Interactive figures (`slider_figure`, plotly, precomputed; ≤ 40 steps):
- **C05** — slider over m/U: the half-body, its stagnation point and dividing streamline, with a and h_max printed; why:
  the body scales with m/U.
- **C06** — C_p(θ) of the cylinder with the qualitative real band (angle from the front, as in Fig. 6.10's idea) and
  a slider over the band's separation angle (labelled qualitative); why: where ideal and real part company.
- **C07** — slider over Γ/(4πaU) from 0 to 2: streamlines and stagnation points, surface C_p and L = ρUΓ; why: the
  four regimes of Fig. 6.12 as one continuous family.
- **C09** — slider over n (α = π/n from π/4 to 2π): the corner flow net and ∣dw/dz∣ near the corner; why: calm vs violent
  corners.
- **C11** — slider over a/b from 1 to 3: the Zhukhovsky image of the circle and its flow (plate → ellipse → circle);
  why: one map, a family of bodies.
- **C13** — slider over r/a: cylinder vs sphere velocity perturbation (1/r² vs 1/r³) and their surface C_p; why: 3-D
  relief.
- **C14** — slider over N: target body vs computed ψ = 0 contour and the k_n bars (Rankine oval target); why: the
  inverse method converging.
- **C15** — a `live(...)` widget for `sphere_motion` (m, a, ρ, F_E) paired with a `slider_figure` over the density ratio
  (bubble … steel ball) showing velocity with and without added mass; why: when added mass matters (light bodies).
Static figures carry C01 (μ∇²u = 0 panel), C02 (flow nets, flux vs radius), C03 (stream + source = sum), C10 (force vs
contour radius, Laurent bars).

## 7. From-scratch moments
Every book section with a computable A item has at least one:
- §6.1 **C01** — central differences of μ∇²u and ∇·u at a point of the cylinder flow and of plane Poiseuille flow by
  hand vs `ch06.ideal_flow_residuals`.
- §6.2 **C02** — a hand-written 5-point −∇²ψ on the Rankine-core ψ vs `ch06.vorticity_from_psi`, and a loop sum of
  ∇ψ·n over a circle vs `ch06.delta_flux_check`; **C03** — the element velocities summed at a point by hand vs
  `core.potential.Flow.velocity`.
- §6.3 **C04** — the source–sink pair formula at shrinking ε vs `core.potential.Doublet` (error ∝ ε²); **C05** —
  `brentq` on tan θ + 2(π − θ) written out vs `ch06.half_body_cp_zero_angle`; **C07** — a trapezoid sum of
  −p(a, θ) sin θ a dθ vs `ch06.lift_per_span` and `ch06.surface_pressure_force`.
- §6.4 **C09** — complex difference quotients along x and along iy by hand vs `Flow.dwdz` (and they disagree for z*).
- §6.5 **C10** — a periodic trapezoid of (dw/dz)² round a circle by hand vs `core.potential.blasius_force`.
- §6.6 **C11** — the quadratic formula with the outside root chosen by ∣ζ∣ > b by hand vs
  `core.conformal.joukowski_inverse`.
- §6.7 **C12** — a double-loop Gauss–Seidel sweep on the 4-point system vs `np.linalg.solve` and
  `core.laplace_solvers.solve_laplace`.
- §6.8 **C13** — sphere ψ assembled from (6.86) + (6.88) by hand vs `core.potential.sphere`; **C14** — the N × N
  influence matrix built by a double loop vs `core.panels.axial_singularity_solve`.
- §6.9 **C15** — a Gauss–Legendre (cos θ) × trapezoid (φ) sum of the surface pressure force vs
  `core.potential.added_mass_sphere` (−M du_s/dt).

## 8. Notes for the implementer
Functions the A items' figures and the explainers need that analysis §4 did not plan (all scalar-callable for parity
rows, SI units, docstrings citing § and Eq.):
- `ch06.rankine_oval(U, m, a)` — promote from "optional" to required (E1 preset): stagnation points, half-length,
  half-width root, ψ function.
- `ch06.superposition_state(elements, probe)` → dict(u, v per element, total, stagnation points, dividing ψ, closed
  body?) — E1's term bars and status in one call.
- `ch06.cylinder_circulation_state(U, a, *, Gamma_cw)` → dict(stagnation points, regime word, D and L by quadrature,
  ρUΓ) — E2's status, force bars and the C07 slider figure.
- `ch06.example_6_1_wall_pressure(y_wall, t, Gamma, h, rho, p_inf)` — p − p∞ along the whole wall x = 0 (our extension
  of the book's origin value) and its split into −ρ∂φ/∂t and −½ρv² — E3 and the C08 animation.
- `core.potential.Corner(A, n, cut_angle)` must work for ½ ≤ n ≤ 4 with the cut outside the wedge; add
  `ch06.corner_speed_exponent(n)` for E4's log–log check.
- `core.potential.blasius_force(dwdz_fn, contour=…)` must accept an arbitrary closed polyline (offset circle, ellipse)
  and flag a contour that crosses the body; `ch06.laurent_contributions(dwdz_fn, R)` → per-power contribution to
  ∮(dw/dz)² dz (E5's bars).
- `core.conformal.joukowski_inverse(z, b, branch="outside" | "principal")` — keep the wrong variant callable (clearly
  labelled) for E6's toggle and the wrong-variant test.
- `core.laplace_solvers.gauss_seidel_sweep(psi, mask, bc)`, `jacobi_sweep`, `sor_sweep(…, omega)` and
  `node_update(psi, i, j)` — single-step functions for E7's transport; `solve_laplace` returns the residual history.
- `core.panels.axial_singularity_solve` targets: `ch06.axisym_body_target(kind, N)` for airship, Rankine oval, sphere,
  ellipsoid; `ch06.panel_cp_error(N)` for the panel-mode convergence curve.
- `core.potential.moving_sphere_surface_pressure(..., split=True)` → dict(steady, acceleration, total) for E9's bars;
  `ch06.sphere_motion(..., mode="ball" | "bubble" | "oscillating")`; `ch06.cylinder_added_mass(a, rho)` (N104 number).
- `ch06.doublet_limit_frames(eps_list, d)` for the C04 animation (ψ grids cached).
- Parity rows must include: stagnation point −m/2πU and h_max = m/2U (E1); sin θ = −Γ/4πaU and L = ρUΓ from the
  pressure integral (E2); drift Γ/4πh and p(0, 0, 0) = −ρΓ²/4π²h² (E3); dw/dz of z^n and the CR residual (E4); Blasius
  L on two contour radii and on the ellipse (E5); ellipse semi-axes and the outside-root inverse at z = −3 + 0.5i (E6);
  the 4-point solution and a Gauss–Seidel sweep (E7); airship closure and the axial solve on a Rankine oval (E8);
  M = 2πρa³/3 by pressure and by energy, bubble acceleration 2g (E9). Book numbers stay in
  `tests/book_values_ch06.json`.
- Budget: flow fields on ≤ 160 × 120 grids; Example 6.2 on the book grid plus refinements ×2, ×4 (×8 only by the sparse
  direct solver); pure-Python Gauss–Seidel only on grids ≤ 40 × 40 (vectorised red-black or SOR above); axial method N ≤
  200; panels N ≤ 256; sympy checks (D17, D18, D21, D29, D31) cached; animations FAST = 30 frames; whole notebook < 5
  min on Colab CPU.

**Conventions for the notebook's notation cell and every explainer (analysis §9, carried forward):**
- **Γ sign.** Project and (6.6), (6.8), (6.47): Γ counterclockwise positive, ψ = −(Γ/2π) ln r, u_θ = +Γ/2πr. The book's
  (6.36)–(6.40), (6.52), (6.61)–(6.62), (6.68) and Example 6.1 use a **clockwise** Γ (the flow's circulation is −Γ), so
  L = +ρUΓ. Code: `Vortex(Gamma)` always ccw; body helpers take keyword-only `Gamma_cw=` (book) or `Gamma_ccw=`
  (project), `Gamma_ccw = −Gamma_cw`; `lift_per_span(rho, U, Gamma_cw=…)` = ρUΓ_cw. ⚠️ callout with numbers: Γ_cw = 2
  m²/s, U = 10 m/s, a = 0.1 m → stagnation points *below* the axis at −9.16°, −170.84°, L = +24 N/m (air); the same flow
  has Γ_ccw = −2 m²/s. Pin with tests (Γ_cw > 0 moves the points down; L > 0) and with exact-text parity rows in E2.
- **Doublet sign.** 2-D: d = Σx_i m_i points from sink to source; (6.29) φ = −d·x/2πr²; the cylinder needs d = −2πUa²
  e_x (upstream); (6.49)'s scalar d is a dipole −d e_x (`Doublet.from_book_scalar`). 3-D: (6.88) dipole −d e_z, sphere
  d = 2πa³U, (6.92) d = −d e_z, moving sphere (6.97) d(t) = +2πa³u_s. Code stores vectors.
- **Angles.** Polar θ from +x (downstream) in (6.31)–(6.39) and (6.89)–(6.91): the upstream stagnation point is θ = π.
  The Fig. 6.10 comparison uses the angle from the upstream stagnation point (= π − θ); the half-body C_p zero is
  measured from +x (our root ≈ 113.2°); (6.106)'s θ_s from the sphere's velocity (= π − θ). Use `np.arctan2`; the
  half-body ψ needs θ with the cut inside the body (test ψ = m/2 on both halves).
- **Axes and symbols.** §6.8: z is the horizontal symmetry axis along the stream; the Stokes ψ is in m³/s (plane ψ
  m²/s) and the flux between stream surfaces is 2πdψ; u_R = −(1/R)∂ψ/∂z (same sign as `core.streamfunction`); spherical
  u_r = (1/(r² sin θ))∂ψ/∂θ. §6.9: w = z-velocity (not the complex potential), ξ = x − x_s (in §6.8 ξ is the axial
  coordinate along the line sink), ζ = the Zhukhovsky variable. Symbols reused inside this chapter: a (radius,
  half-body stagnation distance, image distance, line-sink length, Zhukhovsky circle), α (corner angle, conformal angle,
  line-sink angle), b, B (span — not ch04's Bernoulli function), d, k (line-sink density), M (added mass — not Mach), m
  (2-D source strength vs body mass in (6.109)), n (corner exponent vs normal), q, Q (3-D source vs flow rate), R
  (cylindrical radius vs CV radius), z (complex variable vs axial coordinate). Add every one to `knowledge/notation.md`.
- **Forces.** D, L are forces *on the body* per unit depth; (6.54)'s F is on the fluid (F = −B(De_x + Le_y)); (6.98)'s
  F_s is on the sphere; contours counterclockwise for (6.56) (assert signed area > 0). Pressure is measured from
  hydrostatic in ideal flow; gravity returns only in exercises (bubble, hemisphere).
- **Book slips to teach, never to code:** (6.61) 1/z² coefficient and outer square (true: −(Ud/π + Γ²/4π²)); (6.104)
  middle bracket −u_s/a³ → +u_s/a³; (6.108) stray dφ; §6.3 "(6.8)" → (6.15) for the source velocities; §6.6 "Figure
  6.5" → 6.3; §6.4 "(6.5) and (6.12), respectively" reversed and "(6.43)" → (6.44); §6.5 "Section 3" → §6.3; §6.7
  "first-order" differences are second-order accurate; Example 6.2's FORTRAN loop over I should be over J, Δψ units
  m²/s. ((6.82) is correct as printed — r × App. B.) Exercise 6.23's "N/m" → N/m² and Fig. 6.9's "Idea flow" are cosmetic
  (not taught).
- **Numerical traps:** the Zhukhovsky inverse branch (wrong-variant test mandatory); complex-log branch cuts outside
  the plotted fluid (φ of a vortex jumps by Γ); element centres, re-entrant corners and ζ = ±b excluded from stencils;
  the contraction's 270° corner lowers the convergence order near it (report, fence); the axial method is ill-conditioned
  for blunt bodies (sphere target fenced `qualitative`); Blasius on a polygon converges algebraically (assert the rate);
  drag = 0 checked relative to ½ρU²a.
- **Fig. 6.10's measured curve**: no cited dataset yet (analysis §8) — draw it only as a labelled qualitative band; if
  a source (Roshko 1961, Achenbach 1968) is fetched later, cite it and upgrade to V5.
