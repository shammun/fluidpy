# Chapter 8 — Laminar Flow: curation
(from `analysis/ch08.md` — 143 inventory rows (#1–#143, of which #139–#143 are figure rows and #138 the exercises);
58 equation labels (8.1)–(8.53) with 8.4, 8.13, 8.16, 8.17, 8.32 split a/b; 49 derivations in §2b (★ 18 · ★★ 23 ·
★★★ 8); 18 candidate A groups in §3; `book.yaml` ch08, `policy.tier_a_full_treatment: [12, 18]`, `coverage:
exhaustive`. Curated 2026-09-28, concept-curator. Read: `knowledge/CUMULATIVE.md`, `concept_map.md`, `primers.md`
(P01–P184), `notation.md`, `viz_patterns.md`, `knowledge/ch07.md` §8, `analysis/ch07_curation.md` (format).)

Counts: A 15 · B 113 · C 15 · RECAP 20 · SKIP 1 · derivations written out 33 (★ 8 ★★ 18 ★★★ 7) · demoted to statements 15

Tier words (parsed by `tools/nbkit.py`): **CORE 15 · NOTE 107 · RECAP 20 · SKIP 1** = 143 rows; **DERIVATION 33**.
Depth: A 15 (all CORE) · B 113 (96 NOTE + 17 RECAP) · C 15 (11 NOTE: N02, N11, N37, N64, N74, N75, N76, N80, N94,
N101, N102 · 3 RECAP: R06, R09, R17 · 1 SKIP S01).
**Reconciliation with analysis §2:** the 143 rows `[#1]…[#143]` each appear exactly once in §2 below (checked by
script). A = 15 of 143 rows (10.5 %).

**Decisions that shape this chapter.**
1. **Fifteen A items from the analyst's eighteen groups.** Three merges:
   (a) the parallel-flow reduction (analyst A2: continuity ⇒ v = 0, "a function of x equals a function of y ⇒ both
   constant") is taught **inside C02** (Couette–Poiseuille), because the reduction steps are the first steps of deriving
   $u(y)=\frac Uh y-\frac1{2\mu}\frac{dp}{dx}y(h-y)$ *(Eq. 8.5)*; they keep their own derivations (D02, D03) and are
   re-used by name in C03, C04, C06 and C09;
   (b) the line-vortex decay and spin-up (analyst A13, Example 8.6, Exercise 8.26) become B items of **C10** (the
   similarity ansatz), as the analyst suggested — the result $u_\theta=\frac{\Gamma}{2\pi r}[1-e^{-r^2/4\nu t}]$ is the
   Lamb–Oseen/Gaussian vortex already taught (ch03 C14, ch05 N13), so its ★★★ algebra (a-D32) is **stated** with a
   sympy residual check; the ansatz's worked derivations are Example 8.4 (n = 0), Example 8.5 (n from an integral
   constraint, ★★★) and Example 8.7 (two exponents from a volume constraint);
   (c) the spreading-bead similarity law $h=At^{-1/5}F(x/Dt^{1/5})$ (analyst A10's second half, Example 8.7) is taught in
   **C10** where the ansatz lives (its exponent algebra D23 is the ansatz's third worked case), while the thin-film
   equation itself (Example 8.3) stays its own A item **C08** in §8.3.
   The analyst's other suggested merge (circular Couette into the pipe group) was **not** taken: circular Couette
   **C04** stays an A item because it is the Taylor–Couette base state of Ch. 11 (Rayleigh criterion uses its A and B)
   and the cyclostrophic balance of Ch. 13; instead the analyst's A11 (Stokes' first problem) is one A item **C09**
   carrying the similarity variable, the ODE, the erf solution and δ₉₉ in one derivation chain.
   Hele-Shaw flow (Example 8.2) is a **B item of C06**: the depth-averaged flux $\bar{\mathbf u}=-\frac{h^2}{12\mu}\nabla p$
   plus continuity is the 2-D Reynolds equation with constant gap, so $\nabla^2p=0$ follows in two lines once C06 is
   known; it gets a figure (streamlines round a disc = ideal flow) but its derivation (a-D20) is stated.
2. **SEEN rows are RECAP** (20 rows, R01–R20, inventory order), reminded with the earlier chapter's function. Rows
   marked "SEEN → NEW (derived)" whose new part is the derivation are CORE when the derivation is the A item's heart
   (**C02** (8.5): ch04 used it as `exact_solution("couette")`; **C03** (8.6): ch04 `exact_solution("pipe_poiseuille")`,
   ch03 `pipe_profile`; **C09** (8.30): ch04 `stokes_first_problem`, erfc P123 — each block opens with a 🔁 recap
   sentence naming the earlier test-field use) and NOTE B otherwise (N21 (8.11), N39 (8.20), N43 (8.24), N61 Example
   8.6, N68 (8.35), N92 (8.52)).
3. **Derivations written out: 33** (★ 8 · ★★ 18 · ★★★ 7), covering 34 of the 49 analysis rows (merge: a-D24 + a-D25 →
   D18, the ODE and its two conditions). Rule (b) ("the book never writes it out and the lesson needs it"): **D05**
   (backflow threshold dp/dx > 2μU/h², ours — the explainer's verdict), **D07** (f = 64/Re and the CV force balance,
   ours — Ch. 12's Moody diagram), **D13** (the 1-D Reynolds equation, Exercises 8.19–8.20, never written in the text),
   **D20** (δ₉₉ read off a figure in the book), **D25** (Stokes-layer phase speed and e-folding depth, ours), **D30**
   (the pressure (8.50): the book gives only the result). **15 demoted to statements** (§4c), each with a numeric or
   sympy check in code; the one ★★★ demoted is a-D32 (line-vortex decay, a recap result).
4. **Book slips taught in our own words** (analysis §9; the corrected form is taught, the printed one named in a
   ⚠️ callout and kept as a wrong-variant test): R6 (8.13b) prints ∂p/∂x in the v-equation — read ∂p/∂y (N24);
   R7 (8.17a) and Example 8.2 lack the ν — read $0\cong-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{\partial^2u}{\partial y^2}$
   (C05, D11, N36); R8 Example 8.1's intermediate integrands use (1 − αx/L) and the final p − p_e prints (1 + αx/L) to the
   first power — read (1 + αx/L) throughout and **squared** in the final denominator (C07, D14; the printed curve is a
   ghost in `slider_bearing`); R8b (8.19) adds U₀ instead of U₀(1 − y/h) (N34, D12); R9 the channel mean velocity lacks
   the 1/h (N10); R10 ∫₀^∞ω dy = **+U**, printed −U (N53); R11 Example 8.5's ±2.76 is ±2.772 (N59, D22); R12 the rear
   pressure minimum is **−**3μU/2a (N88, D30); R13 Oseen's equation needs **−**∂p/∂x_i (N95); R14 cross-reference slips
   ("(9.63)" → (8.43), "(9.68)" → (8.48), "(8.33) into (8.20)" → (8.35), "Example 8.2" → 8.3 in Example 8.7, "y" → z
   in Example 8.2, "kinematic" → dynamic viscosity in §8.1) named where they occur; R15 power into the fluid is
   −2πR₁σ_Rφu_φ (R11).
5. **Conventions stated where first used** (⚠️ callouts with numbers): the book passes **dp/dx** (favourable < 0) while
   ch04 code uses G = −dp/dx (C02, C03 — parity rows pin the sign); channel walls at y = 0 and h; capital R is the
   cylindrical radius, r the spherical radius; θ in §8.6 is measured from the **downstream** axis (rear stagnation point
   θ = 0) while Example 8.6's θ is a plane-polar azimuth; φ is an azimuth except in Hele-Shaw (velocity potential);
   ω is vorticity in §8.1, Example 8.5 and §8.6 but the oscillation frequency in §8.5; η = y/√(νt) in (8.25) but the
   book's Figs. 8.13–8.14 plot y/(2√(νt)) (every figure and explainer labels which one); the Stokes-layer depth
   4√(ν/ω) (book) vs √(2ν/ω) (literature e-folding) differ by 2√2; four Reynolds numbers (pipe Ud/ν with the diameter,
   lubrication Re_L with the passage length, generic ρUL/μ, sphere 2aU/ν with the diameter — Oseen's 3/16 becomes 3/8
   with the radius); lubrication p* = p/P_a (atmospheric) makes Λ ≈ 10³ in a real bearing — the natural pressure scale is
   μUL/h² (C05); g = `G0` = 9.80665 in `core.creeping` (settling) vs `G_BOOK` = 9.81 elsewhere.
6. **Climate hooks, where real:** ν as a diffusivity and t ~ L²/ν (C01) is the starting point of eddy viscosity (Ch. 12)
   and of every Ekman spin-up time (Ch. 13); the √(νt) thickness (C09) is the boundary-layer and Ekman-depth scaling;
   the Stokes layer $e^{-y/\delta}\cos(\omega t-y/\delta)$ (C11) has exactly the (1 + i)/δ structure of the Ekman spiral
   and of tidal bottom boundary layers (Ch. 13); the thin-film equation $h_t=\frac{\rho g}{3\mu}(h^3h_x)_x$ (C08) is the
   2-D viscous gravity current — lava, mud flows, and with a nonlinear rheology the shallow-ice equation of glaciology;
   Stokes settling (C14) sets cloud-droplet, aerosol and sediment fall speeds (a 10 µm droplet falls ≈ 1.2 cm/s);
   circular Couette (C04) is the cyclostrophic/gradient-wind balance and Taylor–Couette laboratory analogue of Ch. 13.

## 1. Teaching order (A IDs grouped by book section, B/C IDs under each; one sentence each: "once you see X, Y follows")
A items in **bold**; B and C items listed where they are taught, depth in brackets. Equations written beside their
numbers so downstream agents have them.

**§8.1 Introduction**
- **C01 Laminar vs turbulent flow and ν as a momentum diffusivity** $\mathrm{Re}=Ud/\nu\sim2000\text{–}3000$ [#1]:
  once you see Reynolds's dye streak stay straight below a critical Re (N103 [B], synthetic sketch) and read ν = μ/ρ as
  a diffusivity with diffusion time L²/ν (N01 [B], D01: air diffuses momentum ≈ 15× faster than water), the ratio
  $\frac{U^2/L}{\mu U/\rho L^2}=\mathrm{Re}$ (R01 [B]), the vorticity diffusion $D\omega_z/Dt=\nu\nabla^2\omega_z$
  (R02 [B]), the governing equation $D\mathbf u/Dt=-(1/\rho)\nabla p+\nu\nabla^2\mathbf u$ *(Eq. 8.1)* (R03 [B]) and
  the wall conditions (8.2)–(8.3) (R04, R05 [B]) set the stage; standing assumptions (R06 [C]) and the road map to
  Ch. 9, 11, 12, 13 (N02 [C]).

**§8.2 Exact solutions for steady incompressible viscous flow**
- **C02 Plane Couette–Poiseuille flow by the parallel-flow reduction** $u(y)=\frac Uh y-\frac1{2\mu}\frac{dp}{dx}y(h-y)$
  *(Eq. 8.5)* [#16]: once continuity forces v ≡ 0 in fully developed flow (N03, N04 [B]; D02), the momentum equations
  shrink to $0=-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{d^2u}{dy^2}$ *(Eq. 8.4a)* and
  $0=-\frac1\rho\frac{\partial p}{\partial y}$ *(Eq. 8.4b)* (N05, N06 [B]), "f(x) = g(y) ⇒ both constant" makes dp/dx
  constant (N07 [B]; D03), two integrations and two wall conditions give (8.5) (N08 [B]; D04), and its cases — Couette
  (R07 [B]), Poiseuille with linear τ (R08 [B]), favourable/adverse gradient and **backflow** when
  dp/dx > 2μU/h² (N09 [B]; D05) — and the flow rate $Q=\frac{Uh}{2}[1-\frac{h^2}{6\mu U}\frac{dp}{dx}]$ (N10 [B], R9
  slip) follow; the linear-stress remark (N11 [C] → Ch. 12).
- **C03 Poiseuille pipe flow and Hagen–Poiseuille** $u_z(R)=\frac{R^2-a^2}{4\mu}\frac{dp}{dz}$ *(Eq. 8.6)* [#24]: once
  the same reduction is done in cylindrical coordinates (N12 [B]) and boundedness at the axis kills the ln R term
  (N13 [B]; D06), the stress $\tau=\frac R2\frac{dp}{dz}$ *(Eq. 8.7)* (N14 [B]; R09 [C]), the wall stress
  $\tau_0=\frac a2\frac{dp}{dz}$ *(Eq. 8.8)* as a CV force balance (N15 [B]), $Q=-\frac{\pi a^4}{8\mu}\frac{dp}{dz}$,
  u_max = 2V and the friction factor f = 64/Re (N16 [B]; D07) follow; the entrance-length picture (N03 [B]) is shown here
  too.
- **C04 Circular Couette flow and its two limits**
  $u_\varphi=\frac{1}{R_2^2-R_1^2}\{[\Omega_2R_2^2-\Omega_1R_1^2]R-[\Omega_2-\Omega_1]\frac{R_1^2R_2^2}{R}\}$ *(Eq. 8.10)*
  [#32]: once the φ-momentum equation is a pure viscous balance (N17 [B]) with general solution
  $u_\varphi=AR+B/R$ *(Eq. 8.9)* (N18 [B]; D08) and the constants come from the two walls (N19 [B]), the limits
  $u_\varphi=\Omega_1R_1^2/R$ *(Eq. 8.11)* — the viscous but irrotational vortex (N21 [B], R10 [B], R11 [B]) — and
  $u_\varphi=\Omega_2R$ *(Eq. 8.12)* (R12 [B]) follow (D09); the pressure from the centripetal balance (N20 [B], ours),
  the advective-acceleration check for all three flows (N22 [B]) and our remakes of Figs. 8.5–8.7 (N104 [B]).

**§8.3 Elementary lubrication theory** (book order kept except that Hele-Shaw, Example 8.2, is taught at the end of
C06 as the 2-D Reynolds equation, before the slider bearing's numbers; nothing in C07 depends on it)
- **C05 The lubrication approximation** $0\cong-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{\partial^2u}{\partial y^2}$
  *(Eq. 8.17a)* [#47]: once a thin gap h ≪ L (N23 [B]) is scaled with **two** length scales
  $x^*=x/L,\ y^*=y/h,\ v^*=v/\varepsilon U$ *(Eq. 8.14)* (N26 [B]), continuity (8.15) keeps no coefficient (N27 [B]),
  and the momentum equations (8.13a,b) (R13 [B], N24 [B] with the ∂p/∂y slip) become (8.16a,b) with inertia × ε²Re_L
  and ∂²/∂x² × ε² (N28, N29 [B]; D10 ★★★), the limit ε²Re_L → 0 gives (8.17a) and
  $0\cong-\frac1\rho\frac{\partial p}{\partial y}$ *(Eq. 8.17b)* (N30 [B]; D11 — ν restored); a bearing number
  ε²Re_L ≈ 2.5 × 10⁻³ (N31 [B], with the Λ ≈ 10³ remark).
- **C06 The lubrication profile, gap flux and the Reynolds equation**
  $u\cong-\frac{h^2}{2\mu}\frac{\partial p}{\partial x}\frac yh(1-\frac yh)+U_h\frac yh+U_0$ *(Eq. 8.19)* [#52]: once the
  pressure is independent of y, two integrations with x-dependent "constants" (N33 (8.18) [B]) and the gap wall
  conditions (N25 [B]) give the local Poiseuille + Couette profile (D12; U₀-placement slip N34 [B]); integrating
  continuity across the gap gives $\partial h/\partial t+\partial q/\partial x=0$ with
  $q=-\frac{h^3}{12\mu}\frac{\partial p}{\partial x}+\frac{(U_0+U_h)h}{2}$ (N32 [B]; D13 ★★★), and with constant h in 2-D
  the same flux gives Hele-Shaw's $\nabla^2p=0$ — viscous flow drawing ideal streamlines (N36 [B]).
- **C07 The slider bearing (Example 8.1)**
  $p-p_e=\frac{6\mu LU}{h_o^2}\frac{\alpha(x/L)(1-x/L)}{(2+\alpha)(1+\alpha x/L)^2}$, $W=\frac{\alpha\mu L^2U}{2h_o^2}$
  [#54]: once the flux in the pad frame is a constant C₁ (the steady Reynolds equation of C06), integrating
  dp/dx = −12μC₁/h³ − 6μU/h² with p = p_e at both ends gives the pressure hump and the load (D14 ★★★, printed slips as
  ghosts); the exact load for any taper and the optimum 1 + α ≈ 2.19 (N35 [B], ours, V5 San Andrés) and our remakes of
  Figs. 8.9–8.11 (N105 [B]).
- **C08 The thin-film (viscous gravity current) equation (Example 8.3)**
  $\frac{\partial h}{\partial t}=\frac{\rho g}{3\mu}\frac{\partial}{\partial x}\big(h^3\frac{\partial h}{\partial x}\big)$ [#57]:
  once a hydrostatic pressure ρg(h − y), no slip below and no stress on top turn (8.18) into a half-parabola with flux
  $-\frac{\rho g}{3\mu}h^3h_x$, mass conservation gives a nonlinear diffusion equation for h (D15); the "similarity
  solution" is named here (N37 [C]) and delivered in C10 (N63).

**§8.4 Similarity solutions for unsteady incompressible viscous flow**
- **C09 Stokes' first problem: the impulsively started plate**
  $\frac{u}{U}=1-\mathrm{erf}\big(\frac{y}{2\sqrt{\nu t}}\big)$ *(Eq. 8.30)* [#74]: once the pressure argument reduces NS
  to $\frac{\partial u}{\partial t}=\nu\frac{\partial^2u}{\partial y^2}$ *(Eq. 8.20)* (N38, N39 [B]; D16) with the
  conditions (8.21)–(8.23) (N40–N42 [B]), dimensional analysis (N43 [B]) plus linearity leave one variable
  $\eta=y/\sqrt{\nu t}$ *(Eq. 8.25)* (N44 [B]; D17), the chain rule (N45 [B]) gives the ODE
  $-\frac\eta2F'=F''$ *(Eq. 8.26)* with $F(0)=1$ *(Eq. 8.27)*, $F(\infty)=0$ *(Eq. 8.28)* (N46–N48 [B]; D18), two
  integrations and the Gaussian integral give (8.29) and (8.30) (N49–N51 [B]; D19), and then the collapse (N52, N55
  [B]), the vorticity content +U (N53 [B], R10 slip), the thickness $\delta_{99}\sim3.64\sqrt{\nu t}$ *(Eq. 8.31)*
  (N54 [B]; D20), and why imposed scales destroy similarity (N56 [B], R14 [B] Couette start-up) follow.
- **C10 The similarity ansatz and exponent matching** $\gamma=At^{-n}F(\xi/\delta(t))$ *(Eq. 8.32a)* [#80]: once
  "substitute, divide, and demand every bracket scale with the same power of t" is a procedure (N57 (8.32b) [B]),
  Example 8.4 recovers δ = √(2C₁νt) (N58 [B]; D21), Example 8.5 fixes n = ½ from the conserved velocity jump and gives the
  thickening vortex sheet $u=U\,\mathrm{erf}(y/2\sqrt{\nu t})$ (N59 [B]; D22 ★★★), which is a temporally developing
  boundary layer with $C_f=\frac2{\sqrt\pi}\sqrt{\nu/U^2t}\propto\mathrm{Re}_x^{-1/2}$ (N60 [B]); Example 8.6 gives the
  Lamb–Oseen vortex (N61 [B], stated) and its spin-up twin (N62 [B]); Example 8.7 fixes n = m = 1/5 for the spreading
  bead from the volume constraint (N63 [B]; D23); diffusive vs advective scales (N64 [C]); our Figs. 8.14–8.15 (N106 [B]).

**§8.5 Flow due to an oscillating plate**
- **C11 Stokes' second problem: the Stokes layer**
  $u=Ue^{-y\sqrt{\omega/2\nu}}\cos(\omega t-y\sqrt{\omega/2\nu})$ *(Eq. 8.38)* [#96]: once an imposed period replaces the
  impulsive start (N65 [B]) with $u(0,t)=U\cos\omega t$ *(Eq. 8.33)* and boundedness (8.34) (N66, N67 [B]), the complex
  form $u=\mathrm{Re}\{e^{i\omega t}f(y)\}$ *(Eq. 8.35)* (N68 [B]) turns (8.20) into $i\omega f=\nu f''$ *(Eq. 8.36)*
  (N69 [B]) with $k=\pm(1+i)\sqrt{\omega/2\nu}$ (N70 [B]) and (8.37) (N71 [B]) (D24); reading the result — decay depth,
  phase lag, crest speed √(2νω), 6 % at 4√(ν/ω) (N72 [B]; D25) and why it is not self-similar (N73 [B]) — follows;
  sound absorption at walls (N74 [C] → Ch. 15); Fig. 8.16 remake inside N107 [B].

**§8.6 Low Reynolds number viscous flow past a sphere**
- **C12 Creeping flow: the viscous pressure scale and the Stokes equations** $\nabla p=\mu\nabla^2\mathbf u$ *(Eq. 8.43)*
  [#107]: once the high-Re scaling (8.39)–(8.40) (R15, R16 [B]) multiplied by Re gives the wrong limit 0 = μ∇²u
  (N77 (8.41) [B]), choosing the pressure scale μU/L from the dominant balance (N78 [B]) gives
  $\mathrm{Re}(\mathbf u^*\cdot\nabla^*\mathbf u^*)=-\nabla^*p^*+\nabla^{*2}\mathbf u^*$ *(Eq. 8.42)* (N79 [B]) and, as
  Re → 0, the linear Stokes equations (D26); perturbation view and high-Re non-uniformity (N75, N76 [C] → Ch. 9); the
  italic scaling rule (N80 [C]).
- **C13 Stokes' solution for the sphere** $\psi=Ur^2\sin^2\theta(\frac12-\frac{3a}{4r}+\frac{a^3}{4r^3})$ *(Eq. 8.48)*
  [#118]: once the curl of (8.43) gives ∇²ω = 0 (N81 [B]; D27) with ω_φ from the stream function (R17 [C], R18 (6.83)
  [B], N82 [B]), the fourth-order equation $(E^2)^2\psi=0$ *(Eq. 8.44)* (N83 [B]; D28 ★★★) with no-penetration (8.45),
  no-slip (8.46) and the uniform stream (8.47) (N84, N85 [B], R19 [B]) is solved by f(r) sin²θ (N86 [B]; D29), giving
  (8.48), the velocities (8.49) (N87 [B], stated) and the fore–aft-symmetric fluid-frame pattern (N93 [B]).
- **C14 Stokes drag** $D=6\pi\mu aU$ *(Eq. 8.51)* [#122]: once the pressure
  $p-p_\infty=-\frac{3\mu aU\cos\theta}{2r^2}$ *(Eq. 8.50)* is integrated from ∇p = μ∇²u (N88 [B]; D30 ★★★, R12 sign slip)
  and the surface stresses are known (N89 [B]), the x-traction integrated over the sphere gives 2πμaU from pressure and
  4πμaU from friction (D31 ★★★), and then the terminal velocity (N90 [B]), Millikan's charge measurement (N91 [B]),
  $C_D=24/\mathrm{Re}$ *(Eq. 8.52)* (N92 [B]) and the ρ-free dimensional argument (R20 [B]) follow.
- **C15 Where Stokes fails, and Oseen's fix** inertia/viscous ~ Re r/a [#128]: once the far-field estimate shows inertia
  catching up with viscosity at r/a ~ 1/Re (D32), the non-uniformity and the paradoxes (N94 [C]), Oseen's linearisation
  (N95, N96 [B]), his stream function (8.53) (N97 [B]) that reduces to (8.48) near the sphere (N98 [B]; D33), the drag
  $C_D=\frac{24}{\mathrm{Re}}(1+\frac{3}{16}\mathrm{Re})$ (N99 [B]) and the wake of Fig. 8.20 (N100 [B]) follow; matched
  asymptotics and Proudman–Pearson (N101 [C]).

**§8.7 Final remarks** — taught as the closing paragraph of C15: most tractable laminar problems are solved; perturbation
methods (Ch. 9) and numerical solution of NS (Ch. 10) carry on (N102 [C]); exercises pointer (S01).

## 2. Chapter map (depth) — every inventory row exactly once
One row per inventory row (143: `[#n]` = analysis §2 row n). `A parent` names the A block a B or C item is written in.
IDs: CORE C01–C15 (teaching order), NOTE N01–N107 and RECAP R01–R20 in inventory order, SKIP S01.

| ID | Item | § | Depth | Tier | A parent | Reason (A) / treatment (B) / pointer (C, SKIP) / source chapter (RECAP) |
|---|---|---|---|---|---|---|
| C01 | Laminar vs turbulent flow: Reynolds's dye experiment, transition at a fixed $\mathrm{Re}=Ud/\nu\approx2000\text{–}3000$ (mean velocity, diameter), ν = μ/ρ [#1] | 8.1 | A | CORE | – | load-bearing: every solution of the chapter is "laminar", valid below a critical Re; the Re regime and the L²/ν diffusion time frame Ch. 9, 11, 12 |
| R01 | Inertia/viscous ratio $\frac{U^2/L}{\mu U/\rho L^2}=\mathrm{Re}$; why earlier chapters dropped viscosity at Re ≫ 1 [#2] | 8.1 | B | RECAP | C01 | ch04 C15 (`core.similarity.reynolds_number`); one paragraph + `ch08.inertia_viscous_scales` |
| R02 | 2-D vorticity diffusion $D\omega_z/Dt=\nu\nabla^2\omega_z$ (stretching vanishes in plane flow); ν = diffusivity of vorticity [#3] | 8.1 | B | RECAP | C01 | ch05 C06 (5.13); one sentence on why (ω·∇)u = 0 in 2-D |
| N01 | Heat-equation analogy $DT/Dt=\kappa\nabla^2T$, κ ≡ k/ρC_p (4.89); ν the momentum diffusivity; air ν ≈ 1.5 × 10⁻⁵ m²/s is ≈ 15× water's 1 × 10⁻⁶ m²/s although μ_water ≈ 55 μ_air; diffusion time L²/ν [#4] | 8.1 | B | NOTE | C01 | stated with D01 (the ratio worked) and a two-row table (1 cm: water 100 s, air 6.7 s); `ch08.momentum_diffusivity`, `diffusion_time` |
| R03 | Eq. (8.1): $D\mathbf u/Dt=-(1/\rho)\nabla p+\nu\nabla^2\mathbf u$ (= (4.85)) [#5] | 8.1 | B | RECAP | C01 | ch04 C08; every ch08 solution is fed through `core.navier_stokes.navier_stokes_residual` (residual ≈ 0 shown once per block) |
| R04 | Eq. (8.2): no through-flow $\mathbf n\cdot\mathbf U_s=(\mathbf n\cdot\mathbf u)_{\text{surface}}$ [#6] | 8.1 | B | RECAP | C01 | ch04 C14; `ch08.wall_bc_residuals` |
| R05 | Eq. (8.3): no slip $\mathbf t\cdot\mathbf U_s=(\mathbf t\cdot\mathbf u)_{\text{surface}}$ [#7] | 8.1 | B | RECAP | C01 | ch04 C14; same function |
| R06 | Standing assumptions: constant ρ, inertial frame, gravity absorbed into p (p → p + ρgz) when no free surface [#8] | 8.1 | C | RECAP | C01 | ch04 C13 ("neglect of gravity in constant-density flows"); one sentence |
| N02 | Road map: boundary layers Ch. 9, stability/transition Ch. 11, turbulence Ch. 12, Ekman layers Ch. 13 [#9] | 8.1 | C | NOTE | C01 | named; pointers Ch. 9, 11, 12, 13 |
| N03 | Fully developed flow u = u(y) vs the entrance length (Fig. 8.2): wall layers grow and merge; ∂u/∂x ≠ 0 ⇒ v ≠ 0 upstream [#10] | 8.2 | B | NOTE | C02 | picture paragraph + our sketch with the wall-layer edge from `diffusion_thickness(x/U)` (labelled estimate); pointer Ch. 9 entry flow |
| N04 | Fully developed channel: ∂u/∂x = 0 ⇒ ∂v/∂y = 0 ⇒ v ≡ 0 (walls impermeable); u = (u(y), 0, 0) [#11] | 8.2 | B | NOTE | C02 | steps 1–3 of D02; `ch08.parallel_flow_sympy("channel")` |
| N05 | Eq. (8.4a): $0=-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{d^2u}{dy^2}$ [#12] | 8.2 | B | NOTE | C02 | result of D02 (advective term vanishes exactly — why the problem is linear) |
| N06 | Eq. (8.4b): $0=-\frac1\rho\frac{\partial p}{\partial y}$ ⇒ p = p(x) [#13] | 8.2 | B | NOTE | C02 | result of D02 |
| N07 | "A function of x alone equals a function of y alone ⇒ both the same constant": $\frac1\mu\frac{dp}{dx}=\frac{d^2u}{dy^2}=\text{const}$; p falls linearly [#14] | 8.2 | B | NOTE | C02 | D03 (differentiate each side by the other variable); re-used by name in C03, C04, C06 |
| N08 | Twice-integrated momentum $0=-\frac{y^2}{2}\frac{dp}{dx}+\mu u+Ay+B$; B = 0, A = (h/2)dp/dx − μU/h [#15] | 8.2 | B | NOTE | C02 | steps of D04 |
| C02 | Eq. (8.5): Couette–Poiseuille profile $u(y)=\frac Uh y-\frac1{2\mu}\frac{dp}{dx}y(h-y)$ — linear (wall-driven) + parabolic (pressure-driven) [#16] | 8.2 | A | CORE | – | load-bearing: the parallel-flow reduction it is derived by recurs in the pipe, circular Couette, lubrication and Stokes' problems; base flow of Ch. 11 stability and the laminar reference of Ch. 12; superposition of linear solutions made visible |
| N09 | Cases of Fig. 8.4: favourable dp/dx < 0 (fuller), adverse dp/dx > 0 (backflow near the fixed wall), Couette, Poiseuille; backflow iff $\frac{dp}{dx}>\frac{2\mu U}{h^2}$ (ours) [#17] | 8.2 | B | NOTE | C02 | D05 (rule b) + figure/explainer verdict; `core.laminar.channel_backflow_threshold`; pointer Ch. 9 adverse gradient and separation |
| N10 | Flow rate and mean velocity $Q=\frac{Uh}{2}[1-\frac{h^2}{6\mu U}\frac{dp}{dx}]$, $V=Q/h$ (printed V lacks the 1/h, R9) [#18] | 8.2 | B | NOTE | C02 | stated (two polynomial integrals, §4c) with `channel_flow_rate` vs `quad`; ⚠️ the R9 slip by a units check |
| R07 | Plane Couette $u=Uy/h$, τ = μU/h uniform [#19] | 8.2 | B | RECAP | C02 | ch01 C12 (steady Couette), ch04 `exact_solution("couette")`; dp/dx = 0 case of (8.5) |
| R08 | Plane Poiseuille $u=-\frac1{2\mu}\frac{dp}{dx}y(h-y)$, $\tau=-(\frac h2-y)\frac{dp}{dx}$, ∣τ_w∣ = (h/2)∣dp/dx∣ [#20] | 8.2 | B | RECAP | C02 | ch04 `plane_poiseuille` (**G = −dp/dx**); ⚠️ sign convention callout; parity row |
| N11 | Constant dp/dx and linear τ(y) are general for fully developed channel flow and persist for turbulent averages [#21] | 8.2 | C | NOTE | C02 | named; pointer Ch. 12 (linear total stress in channel flow) |
| N12 | Round-tube set-up in (R, φ, z): u = (0, 0, u_z(R)), ∂p/∂φ = ∂p/∂R = 0 ⇒ p = p(z) [#22] | 8.2 | B | NOTE | C03 | steps 1–3 of D06; `ch08.parallel_flow_sympy("pipe")` with `core.curvilinear` |
| N13 | z-momentum $0=-\frac{dp}{dz}+\frac\mu R\frac{d}{dR}(R\frac{du_z}{dR})$ and $u_z=\frac{R^2}{4\mu}\frac{dp}{dz}+A\ln R+B$; A = 0 by boundedness, B from no slip [#23] | 8.2 | B | NOTE | C03 | steps of D06 (regularity at the axis) |
| C03 | Eq. (8.6): Poiseuille pipe profile $u_z(R)=\frac{R^2-a^2}{4\mu}\frac{dp}{dz}$ [#24] | 8.2 | A | CORE | – | load-bearing: Hagen–Poiseuille Q ∝ a⁴, the wall stress (8.8) and f = 64/Re are the laminar line of Ch. 12's Moody diagram, u_τ's definition and Ch. 16's arteries |
| R09 | Cylindrical shear stress $\tau_{zR}=\mu(\partial u_R/\partial z+\partial u_z/\partial R)$ [#25] | 8.2 | C | RECAP | C03 | ch04 C07 / `core.curvilinear.strain_rate("cylindrical")`; one sentence |
| N14 | Eq. (8.7): $\tau=\mu\frac{\partial u_z}{\partial R}=\frac R2\frac{dp}{dz}$ — linear in R [#26] | 8.2 | B | NOTE | C03 | step of D07; `pipe_shear_stress` |
| N15 | Eq. (8.8): wall stress $\tau_0=\frac a2\frac{dp}{dz}$ (maximum ∣τ∣); valid for turbulent averages because it is a CV force balance πa²Δp = 2πaLτ₀ (ours) [#27] | 8.2 | B | NOTE | C03 | D07 steps (CV balance, rule b); `pipe_wall_stress`; pointer Ch. 12 friction velocity |
| N16 | Hagen–Poiseuille $Q=-\frac{\pi a^4}{8\mu}\frac{dp}{dz}$, $V=-\frac{a^2}{8\mu}\frac{dp}{dz}$; u_max = 2V; f = 64/Re (ours) [#28] | 8.2 | B | NOTE | C03 | D07 (area element 2πR dR); number a = 1 mm, dp/dz = −1000 Pa/m → V = 0.125 m/s, Re = 250, f = 0.256; `pipe_flow_rate`, `pipe_friction_factor` |
| N17 | Circular Couette set-up: u = (0, u_φ(R), 0); R-momentum $-\frac{u_\varphi^2}{R}=-\frac1\rho\frac{dp}{dR}$ (centripetal), φ-momentum $0=\mu\frac{d}{dR}[\frac1R\frac{d}{dR}(Ru_\varphi)]$ [#29] | 8.2 | B | NOTE | C04 | steps 1–3 of D08 |
| N18 | Eq. (8.9): $u_\varphi=AR+B/R$ — solid-body rotation + line vortex [#30] | 8.2 | B | NOTE | C04 | step of D08 (Euler ODE primer) |
| N19 | Constants $A=\frac{\Omega_2R_2^2-\Omega_1R_1^2}{R_2^2-R_1^2}$, $B=-\frac{(\Omega_2-\Omega_1)R_1^2R_2^2}{R_2^2-R_1^2}$ [#31] | 8.2 | B | NOTE | C04 | last steps of D08 (2 × 2 system); `circular_couette(return_coeffs=True)` |
| C04 | Eq. (8.10): circular Couette profile $u_\varphi(R)=\frac{1}{R_2^2-R_1^2}\{[\Omega_2R_2^2-\Omega_1R_1^2]R-[\Omega_2-\Omega_1]\frac{R_1^2R_2^2}{R}\}$ [#32] | 8.2 | A | CORE | – | load-bearing: the Taylor–Couette base state of Ch. 11 (Rayleigh's criterion is written in its A, B) and the only exact rotating shear flow; its two limits tie the chapter to ch05's ideal vortex and solid-body rotation |
| N20 | Pressure from the R-balance (ours; the book says only "can be determined"): $p=p_1+\rho[\frac{A^2}{2}(R^2-R_1^2)+2AB\ln\frac{R}{R_1}-\frac{B^2}{2}(\frac1{R^2}-\frac1{R_1^2})]$ [#33] | 8.2 | B | NOTE | C04 | stated (one integral, §4c) with a finite-difference check of dp/dR = ρu²/R; pointer Ch. 13 cyclostrophic balance |
| N21 | Eq. (8.11): R₂ → ∞, Ω₂ = 0: $u_\varphi=\Omega_1R_1^2/R$, Γ = 2πΩ₁R₁² — the ideal vortex (5.2) as a viscous solution (Fig. 8.7 free-surface dip) [#34] | 8.2 | B | NOTE | C04 | D09; recap ch05 `rotating_cylinder_flow` (⚠️ its ω = 2Ω₁) |
| R10 | The only irrotational viscous solution: stress $\sigma_{R\varphi}=-2\mu\Omega_1R_1^2/R^2$ but zero net viscous force (∇²u = −∇×ω = 0) [#35] | 8.2 | B | RECAP | C04 | ch05 D03, ch04 (4.40); `circular_couette_shear_stress` |
| R11 | Power into the fluid −2πR₁σ_Rφu_φ = integrated dissipation = 4πμΩ₁²R₁² (Exercise 8.12; printed sign slip R15) [#36] | 8.2 | B | RECAP | C04 | ch05 N08 `dissipation_outside_cylinder`; stated (§4c), `circular_couette_power` |
| R12 | Eq. (8.12): R₁, Ω₁ → 0: $u_\varphi=\Omega_2R$ — solid-body rotation (5.1) [#37] | 8.2 | B | RECAP | C04 | ch05 C02; last step of D09; parity `core.vortices.solid_body_rotation` |
| N22 | All three flows are confined and symmetry kills u·∇u (circular Couette keeps −u_φ²/R e_R, balanced by pressure); other exact solutions in the literature [#38] | 8.2 | B | NOTE | C04 | one paragraph + `ch08.advective_acceleration_check`; pointer Ch. 10 (exact solutions verify codes) |
| N23 | Lubrication idea: nearly parallel flow in a passage h ≪ L (Fig. 8.8); pressure and viscous forces dominate; curvature irrelevant if its radius ≫ h; same logic → boundary layers (§9.1) [#39] | 8.3 | B | NOTE | C05 | the opening picture of C05; pointer Ch. 9 |
| R13 | Gap field equations: (6.2) continuity and Eq. (8.13a) $\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=-\frac1\rho\frac{\partial p}{\partial x}+\frac\mu\rho(\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2})$ [#40] | 8.3 | B | RECAP | C05 | ch06 (6.2), ch04 C08 components of (8.1); starting line of D10 |
| N24 | Eq. (8.13b): $\frac{\partial v}{\partial t}+u\frac{\partial v}{\partial x}+v\frac{\partial v}{\partial y}=-\frac1\rho\frac{\partial p}{\partial y}+\frac\mu\rho(\frac{\partial^2v}{\partial x^2}+\frac{\partial^2v}{\partial y^2})$ (printed with ∂p/∂x, R6) [#41] | 8.3 | B | NOTE | C05 | starting line of D10; ⚠️ R6 callout (the scaled (8.16b) confirms ∂p/∂y) |
| N25 | Gap boundary conditions u = U₀(t) on y = 0, u = U_h(t) on y = h(x, t); p may depend on t [#42] | 8.3 | B | NOTE | C06 | stated where D12 applies them |
| N26 | Eq. (8.14): $x^*=\frac xL,\ y^*=\frac yh=\frac{y}{\varepsilon L},\ t^*=\frac{Ut}{L},\ u^*=\frac uU,\ v^*=\frac{v}{\varepsilon U},\ p^*=\frac{p}{P_a}$; ε = h/L fineness ratio [#43] | 8.3 | B | NOTE | C05 | the scalings of D10 with why v ~ εU (continuity balance); anisotropic-scaling primer |
| N27 | Eq. (8.15): $\frac{\partial u^*}{\partial x^*}+\frac{\partial v^*}{\partial y^*}=0$ — no coefficient: v is small because u varies slowly [#44] | 8.3 | B | NOTE | C05 | step of D10 |
| N28 | Eq. (8.16a): $\varepsilon^2\mathrm{Re}_L(\frac{\partial u^*}{\partial t^*}+u^*\frac{\partial u^*}{\partial x^*}+v^*\frac{\partial u^*}{\partial y^*})=-\frac1\Lambda\frac{\partial p^*}{\partial x^*}+\varepsilon^2\frac{\partial^2u^*}{\partial x^{*2}}+\frac{\partial^2u^*}{\partial y^{*2}}$; Re_L = ρUL/μ, Λ = μUL/(P_a h²) [#45] | 8.3 | B | NOTE | C05 | result of D10 steps; term bars in `lubrication_scaling` |
| N29 | Eq. (8.16b): $\varepsilon^4\mathrm{Re}_L(\ldots)=-\frac1\Lambda\frac{\partial p^*}{\partial y^*}+\varepsilon^4\frac{\partial^2v^*}{\partial x^{*2}}+\varepsilon^2\frac{\partial^2v^*}{\partial y^{*2}}$ — only pressure lacks a power of ε [#46] | 8.3 | B | NOTE | C05 | result of D10's last steps; pointer Ch. 9 (∂p/∂y ≈ 0 across a boundary layer) |
| C05 | Eq. (8.17a): lubrication balance $0\cong-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{\partial^2u}{\partial y^2}$ for ε²Re_L → 0 (printed without ν, R7) [#47] | 8.3 | A | CORE | – | load-bearing: the two-length-scale ordering argument is the template of the boundary-layer approximation (Ch. 9) and of hydrostatic thin layers (Ch. 13); every §8.3 result rests on it |
| N30 | Eq. (8.17b): $0\cong-\frac1\rho\frac{\partial p}{\partial y}$ — p uniform across the gap [#48] | 8.3 | B | NOTE | C05 | result of D11; pointer Ch. 13 (hydrostatic balance in thin layers) |
| N31 | Numerical example ε²Re_L = (h/L)²(UL/ν) for an oil film (book value private); ours: h = 50 µm, L = 5 cm, U = 5 m/s, ν = 10⁻⁴ m²/s → ε²Re_L = 2.5 × 10⁻³; with P_a the bearing number is ≈ 10³, so μUL/h² is the natural pressure scale (R24) [#49] | 8.3 | B | NOTE | C05 | stated with `core.lubrication.lubrication_scales`; ⚠️ Λ callout |
| N32 | Completing a lubrication problem: profile + mass conservation + pressure BCs; gap flux $q=-\frac{h^3}{12\mu}\frac{\partial p}{\partial x}+\frac{(U_0+U_h)h}{2}$ and the 1-D Reynolds equation $\frac{\partial h}{\partial t}+\frac{\partial q}{\partial x}=0$ (Exercises 8.19–8.20; never written in the text) [#50] | 8.3 | B | NOTE | C06 | D13 (rule b, ★★★); `lubrication_flux`, `reynolds_pressure_1d` |
| N33 | Eq. (8.18): $u\cong\frac1\mu\frac{\partial p(x,t)}{\partial x}\frac{y^2}{2}+Ay+B$ (A, B depend on x, t) [#51] | 8.3 | B | NOTE | C06 | steps of D12 |
| C06 | Eq. (8.19): lubrication profile $u\cong-\frac{h^2}{2\mu}\frac{\partial p}{\partial x}\frac yh(1-\frac yh)+U_h\frac yh+U_0$ (consistent form: $U_0(1-\frac yh)$) — local Poiseuille + Couette [#52] | 8.3 | A | CORE | – | load-bearing: with the gap flux it gives the Reynolds equation behind the slider bearing, Hele-Shaw and every lubrication/biofluid film (Ch. 16); it is (8.5) with a gap that varies slowly |
| N34 | The printed (8.19) adds U₀ to the Couette part, so u(h) = U_h + U₀; consistent form $U_h\frac yh+U_0(1-\frac yh)$; irrelevant when U₀ = 0 (every example) (R8b) [#53] | 8.3 | B | NOTE | C06 | ⚠️ callout inside D12's check; `lubrication_velocity(form="consistent")`, `form="book"` as wrong variant |
| C07 | Example 8.1 slider bearing: $p-p_e=\frac{6\mu LU}{h_o^2}\frac{\alpha(x/L)(1-x/L)}{(2+\alpha)(1+\alpha x/L)^2}$ (printed with the first power, R8), $p-p_e\cong\frac{3\alpha\mu LU}{h_o^2}\frac xL(1-\frac xL)$, load $W=\frac{\alpha\mu L^2U}{2h_o^2}$; works in one direction only [#54] | 8.3 | A | CORE | – | load-bearing: the chapter's key worked example of how viscosity in a converging gap builds a pressure thousands of times atmospheric and carries a load; the pattern (flux constant → integrate → pressure BCs) is the lubrication method |
| N35 | Exact load for any taper (ours): $W=\frac{6\mu UL^2}{h_o^2\alpha^2}[\ln(1+\alpha)-\frac{2\alpha}{2+\alpha}]$; optimum inlet/outlet ratio 1 + α ≈ 2.19 (San Andrés, V5) [#55] | 8.3 | B | NOTE | C07 | stated (§4c) with `slider_bearing_load(model="exact")` vs `quad`, `slider_optimum_taper`; number: L = 5 cm, h₀ = 50 µm, U = 5 m/s, μ = 0.05 Pa s: α = 0.1 → 12.5 (linear) vs 10.8 kN/m (exact); α = 0.5 → 62.5 vs 32.8 kN/m |
| N36 | Example 8.2 Hele-Shaw flow: $u\cong-\frac{1}{2\mu}\frac{\partial p}{\partial x}z(h-z)=\frac{\partial\phi}{\partial x}$, $\phi=-\frac{z(h-z)}{2\mu}p$, $\frac{\partial^2p}{\partial x^2}+\frac{\partial^2p}{\partial y^2}=0$ — the streamlines of 2-D ideal flow (6.10), (6.12), except O(h) layers at obstacles; printed without ν and with y for z (R7, R14) [#56] | 8.3 | B | NOTE | C06 | stated (§4c: the 2-D Reynolds equation with constant h, depth average $\bar{\mathbf u}=-\frac{h^2}{12\mu}\nabla p$) with a figure of dye lines round a disc from `ch08.hele_shaw_cylinder` vs `core.potential` streamlines; recap ch06 C01–C04 |
| C08 | Example 8.3 gravity spreading of a viscous bead: $u\cong-\frac{\rho g}{2\mu}\frac{\partial h}{\partial x}y(2h-y)$, flux $-\frac{\rho g}{3\mu}h^3\frac{\partial h}{\partial x}$, thin-film equation $\frac{\partial h}{\partial t}=\frac{\rho g}{3\mu}\frac{\partial}{\partial x}(h^3\frac{\partial h}{\partial x})$ [#57] | 8.3 | A | CORE | – | load-bearing: the first nonlinear diffusion equation of the book and the viscous gravity current of Ch. 13 (lava, mud, ice sheets); it shows lubrication with a free surface and feeds the t^{1/5} similarity law of C10 |
| N37 | Def: similarity solution — a change of variables that turns the PDE into an ODE (introduced in Example 8.3) [#58] | 8.3 | C | NOTE | C08 | named; developed in C09 and C10 (§8.4) |
| R14 | Impulsively started parallel flows keep u∂u/∂x = 0; Couette start-up (Exercise 8.31) — an eigenfunction series, no similarity because the gap sets a length [#59] | 8.4 | B | RECAP | C09 | ch01 C12 `core.diffusion.couette_startup_profile` and ch01 explainer `viscosity_momentum_diffusion`; contrast figure with (8.30) |
| N38 | Stokes' first problem set-up (Fig. 8.12; "Rayleigh problem" naming): plate at y = 0 given speed U at t = 0; v = 0; ∂p/∂x = 0 everywhere because the fluid far away is still at rest [#60] | 8.4 | B | NOTE | C09 | steps of D16 |
| N39 | Eq. (8.20): the 1-D diffusion equation $\frac{\partial u}{\partial t}=\nu\frac{\partial^2u}{\partial y^2}$ [#61] | 8.4 | B | NOTE | C09 | result of D16; 🔁 recap of ch01 FTCS diffusion; new `core.diffusion.crank_nicolson_1d` as the numerical twin |
| N40 | Eq. (8.21): $u(y,t=0)=0$ [#62] | 8.4 | B | NOTE | C09 | stated in D16's result box |
| N41 | Eq. (8.22): $u(y=0,t)=0$ for t < 0, U for t ≥ 0 [#63] | 8.4 | B | NOTE | C09 | stated with D16 |
| N42 | Eq. (8.23): $u(y\to\infty,t)=0$; well posed: two conditions in y, one in t [#64] | 8.4 | B | NOTE | C09 | stated with a gloss on order vs number of conditions |
| N43 | Eq. (8.24): $u/U=f(y/\sqrt{\nu t},\ y/Ut)$ from dimensional analysis [#65] | 8.4 | B | NOTE | C09 | steps 1–3 of D17; recap ch01 Π method (`ch01.pi_groups`) |
| N44 | Eq. (8.25): linearity ⇒ u ∝ U ⇒ drop y/Ut; $u/U=F(y/\sqrt{\nu t})\equiv F(\eta)$ [#66] | 8.4 | B | NOTE | C09 | result of D17; `core.laminar.similarity_variable` |
| N45 | Chain-rule derivatives $\frac{\partial u}{\partial t}=-\frac{U\eta}{2t}F'$, $\frac{\partial^2u}{\partial y^2}=\frac{U}{\nu t}F''$ [#67] | 8.4 | B | NOTE | C09 | steps of D18 |
| N46 | Eq. (8.26): $-\frac\eta2\frac{dF}{d\eta}=\frac{d}{d\eta}(\frac{dF}{d\eta})$ [#68] | 8.4 | B | NOTE | C09 | result of D18 (t cancels — the test that the similarity form is right); `similarity_ode_solve` (`solve_bvp`) check |
| N47 | Eq. (8.27): $F(\eta=0)=1$ [#69] | 8.4 | B | NOTE | C09 | step of D18 |
| N48 | Eq. (8.28): $F(\eta\to\infty)=0$ — (8.21) and (8.23) collapse into one condition [#70] | 8.4 | B | NOTE | C09 | step of D18 (both limits give η → ∞) |
| N49 | Separate and integrate: $F'=A\exp(-\eta^2/4)$ [#71] | 8.4 | B | NOTE | C09 | steps of D19 |
| N50 | Eq. (8.29): $F(\eta)=A\int_0^\eta\exp(-\xi^2/4)d\xi+B$ [#72] | 8.4 | B | NOTE | C09 | step of D19 |
| N51 | Constants B = 1, A = −1/√π via ξ = 2ζ and the Gaussian integral [#73] | 8.4 | B | NOTE | C09 | steps of D19 (Gaussian-integral primer) + `quad` check |
| C09 | Eq. (8.30): Stokes' first problem $\frac{u}{U}=1-\mathrm{erf}(\frac{y}{2\sqrt{\nu t}})$, $\mathrm{erf}(\zeta)=\frac2{\sqrt\pi}\int_0^\zeta e^{-\xi^2}d\xi$ [#74] | 8.4 | A | CORE | – | load-bearing: the √(νt) diffusion scale and the idea of a similarity variable come from here — Ch. 9 (δ ∝ √(νx/U), Blasius), Ch. 13 (Ekman spin-up), Ch. 12; the chapter's first similarity solution derived end to end |
| N52 | Fig. 8.13: u/U vs y/(2√(νt)) — the collapse; Fig. 8.12 dimensional profiles at t₁ < t₂ [#75] | 8.4 | B | NOTE | C09 | our remake (profiles raw and rescaled, CN dots overlaid); ⚠️ which variable the axis shows (R22) |
| N53 | Vorticity reading: a vortex sheet created at t = 0 diffuses; $\int_0^\infty\omega\,dy=u(0)-u(\infty)=U$ (printed −U, R10), constant ⇒ no vorticity generated after t = 0 [#76] | 8.4 | B | NOTE | C09 | stated (§4c) with `ch08.vorticity_content` (`quad`) = U; ⚠️ R10 callout |
| N54 | Eq. (8.31): $\delta_{99}\sim3.64\sqrt{\nu t}$ (u = 0.01U at η = 2 erfc⁻¹(0.01) = 3.643); grows as t^{1/2} [#77] | 8.4 | B | NOTE | C09 | D20 (rule b: the book reads it off a figure); number ν = 10⁻⁶, t = 100 s → 3.64 cm, t = 1 h → 21.9 cm; `diffusion_thickness` |
| N55 | Self-similarity: profiles at any t, U, ν collapse when u/U is plotted against y/√(νt) [#78] | 8.4 | B | NOTE | C09 | collapse figure (max spread < 10⁻¹² asserted) |
| N56 | Similarity needs no imposed length or time: a second plate or stopping the plate at t = T (Exercise 8.30) destroys it; $u=U[\mathrm{erfc}\frac{y}{2\sqrt{\nu t}}-\mathrm{erfc}\frac{y}{2\sqrt{\nu(t-T)}}]$ for t > T (ours) [#79] | 8.4 | B | NOTE | C09 | stated with `stokes_first_stopped` and the R14 contrast figure; pointer Ch. 13 spin-down |
| C10 | Eq. (8.32a): the similarity ansatz $\gamma=At^{-n}F(\xi/\delta(t))\equiv At^{-n}F(\eta)$ with exponent matching (each bracket ∝ the same power of t) [#80] | 8.4 | A | CORE | – | load-bearing: the general procedure behind every self-similar solution later — Blasius and Falkner–Skan (Ch. 9), jets, wakes and plumes (Ch. 12), gravity currents (Ch. 13); here it yields the vortex sheet, the temporal boundary layer and the t^{1/5} bead |
| N57 | Eq. (8.32b): $\gamma=A\xi^{-n}F(\xi/\delta(t))$ — a space power when ξ appears in the initial condition [#81] | 8.4 | B | NOTE | C10 | stated; used by N61 |
| N58 | Example 8.4: Stokes' first problem via (8.32a), A = 1, n = 0: $-[\frac{\delta'}{\delta}]\eta F'=[\frac\nu{\delta^2}]F''$ ⇒ δδ′ = C₁ν ⇒ δ = √(2C₁νt) [#82] | 8.4 | B | NOTE | C10 | D21 (first worked case of the ansatz) |
| N59 | Example 8.5 viscous vortex sheet: ±U at t = 0; n = ½ from the conserved jump; $\omega_z=-\frac{U}{\sqrt{\pi\nu t}}e^{-y^2/4\nu t}$, $u=U\,\mathrm{erf}(\frac{y}{2\sqrt{\nu t}})$; u = ±0.95U at η = ±2.772 (printed 2.76, R11), width 5.54√(νt) [#83] | 8.4 | B | NOTE | C10 | D22 (★★★, the ansatz with an integral constraint); `vortex_sheet_diffusion`, `transition_width`; parity `ch05.diffusing_vortex_sheet(gamma=−2U)` |
| N60 | Temporally developing boundary layer: upper half of Example 8.5 after a Galilean shift = Stokes' first problem; $\tau_w=\frac{\mu U}{\sqrt{\pi\nu t}}$, $C_f=\frac{2}{\sqrt\pi}\sqrt{\frac{\nu}{U^2t}}$; with Ut → x, C_f = 1.128 Re_x^{−1/2} [#84] | 8.4 | B | NOTE | C10 | stated (§4c) with `temporal_bl_wall_stress`; pointer Ch. 9 (Blasius 0.664 Re_x^{−1/2}: same exponent, different constant) |
| N61 | Example 8.6 decay of a line vortex via (8.32b): $(\frac1\eta-\frac\eta2)F'=F''$, $F=1-e^{-\eta^2/4}$, $u_\theta=\frac{\Gamma}{2\pi r}[1-e^{-r^2/4\nu t}]$ — solid-body core, ideal vortex outside (Fig. 8.15) [#85] | 8.4 | B | NOTE | C10 | stated (§4c; a-D32 ★★★ demoted, sympy residual cell shown); recap ch03 C14 Gaussian vortex σ² = 4νt, ch05 N13; `line_vortex_decay` parity `core.vortices.gaussian_vortex` |
| N62 | Line vortex suddenly introduced (Exercise 8.26): $u_\theta=\frac{\Gamma}{2\pi r}e^{-r^2/4\nu t}$ [#86] | 8.4 | B | NOTE | C10 | stated with a residual check (§4c); `line_vortex_spinup` |
| N63 | Example 8.7 similarity of the spreading bead: δ = Dt^m, −1 = −3n − 2m and constant volume −n + m = 0 ⇒ n = m = 1/5: $h=At^{-1/5}F(x/Dt^{1/5})$ (printed "Example 8.2" means 8.3, R14) [#87] | 8.4 | B | NOTE | C10 | D23 (third worked case of the ansatz); Huppert shape and η_N = 1.411 stated (V5); `viscous_current_similarity`; shown in `viscous_gravity_current` |
| N64 | Diffusive lengths grow as (νt)^{1/2}; the bead's t^{1/5} is an advective (travel-distance) scale [#88] | 8.4 | C | NOTE | C10 | named in one sentence at the end of C10 |
| N65 | Stokes' second problem set-up: plate oscillating in its own plane; periodic steady state only; (8.20) again; an imposed time scale 1/ω [#89] | 8.5 | B | NOTE | C11 | the opening picture of C11 |
| N66 | Eq. (8.33): $u(y=0,t)=U\cos(\omega t)$ [#90] | 8.5 | B | NOTE | C11 | start of D24 |
| N67 | Eq. (8.34): $u(y\to\infty,t)=$ bounded [#91] | 8.5 | B | NOTE | C11 | start of D24 |
| N68 | Eq. (8.35): $u(y,t)=\mathrm{Re}\{e^{i\omega t}f(y)\}$ (text's "substitution of (8.33)" means (8.35), R14) [#92] | 8.5 | B | NOTE | C11 | step of D24; recap ch07 P176 complex amplitudes |
| N69 | Eq. (8.36): $i\omega f=\nu\frac{d^2f}{dy^2}$ [#93] | 8.5 | B | NOTE | C11 | step of D24 |
| N70 | $k=(i\omega/\nu)^{1/2}=\pm(1+i)(\omega/2\nu)^{1/2}$ [#94] | 8.5 | B | NOTE | C11 | step of D24 (√i = (1 + i)/√2, P159 reminder; `np.sqrt(1j)`) |
| N71 | Eq. (8.37): $f=A\exp\{-(i+1)y\sqrt{\omega/2\nu}\}+B\exp\{+(i+1)y\sqrt{\omega/2\nu}\}$; B = 0, A = U [#95] | 8.5 | B | NOTE | C11 | steps of D24 |
| C11 | Eq. (8.38): Stokes' second problem $u=U\exp\{-y\sqrt{\frac{\omega}{2\nu}}\}\cos(\omega t-y\sqrt{\frac{\omega}{2\nu}})$ [#96] | 8.5 | A | CORE | – | load-bearing: the oscillatory boundary layer with the (1 + i)/δ structure of the Ekman layer and tidal bottom layers (Ch. 13), Womersley flow (Ch. 16) and acoustic boundary layers (Ch. 15); contrasts an imposed time scale with similarity |
| N72 | Reading (8.38): a damped "wave" moving to +y by diffusion (no restoring force); amplitude $Ue^{-4/\sqrt2}\approx0.06U$ at y = 4(ν/ω)^{1/2}; δ ~ 4(ν/ω)^{1/2}; ours: crest speed √(2νω), wavelength 2π√(2ν/ω), e-folding depth √(2ν/ω) [#97] | 8.5 | B | NOTE | C11 | D25 (rule b); numbers ω = 2π rad/s, ν = 10⁻⁶: √(2ν/ω) = 0.56 mm, 4√(ν/ω) = 1.6 mm; M₂ tide, molecular ν: 0.12 m; `stokes_layer` |
| N73 | Not self-similar: three groups u/U, ωt, y(ω/ν)^{1/2}; extent (ν/ω)^{1/2} set by viscosity × the imposed period [#98] | 8.5 | B | NOTE | C11 | stated with `ch01.pi_groups` (recap ch01 Π) |
| N74 | Application: (8.38) predicts the weak absorption of sound at flat walls [#99] | 8.5 | C | NOTE | C11 | named; pointer Ch. 15 |
| N75 | Small/large-parameter view of flow problems (perturbation expansions) [#100] | 8.6 | C | NOTE | C12 | named; pointers Ch. 9, Ch. 11 |
| R15 | Eq. (8.39): $\rho\mathbf u\cdot\nabla\mathbf u+\nabla p=\mu\nabla^2\mathbf u$ [#101] | 8.6 | B | RECAP | C12 | ch04 (4.39b) steady; start of D26 |
| R16 | Eq. (8.40): $\mathbf u^*\cdot\nabla^*\mathbf u^*+\nabla^*p^*=\frac{1}{\mathrm{Re}}\nabla^{*2}\mathbf u^*$ with p* = (p − p∞)/ρU² [#102] | 8.6 | B | RECAP | C12 | ch04 C15 (4.100)–(4.101), `core.similarity.nondimensional_ns_coefficients(pressure_scale="dynamic")` |
| N76 | High-Re non-uniformity: near a wall the inviscid solution cannot meet no-slip; the 1/Re expansion is singular there; (8.14) is the hint [#103] | 8.6 | C | NOTE | C12 | named; pointer Ch. 9 |
| N77 | Eq. (8.41): $\mathrm{Re}(\mathbf u^*\cdot\nabla^*\mathbf u^*+\nabla^*p^*)=\nabla^{*2}\mathbf u^*$ — Re → 0 would also kill pressure: the wrong equation 0 = μ∇²u [#104] | 8.6 | B | NOTE | C12 | step of D26 |
| N78 | Low-Re pressure scale p* = (p − p∞)L/(μU): pressure differences are set by viscous stress [#105] | 8.6 | B | NOTE | C12 | step of D26 (dominant balance); `nondimensional_ns_coefficients(pressure_scale="viscous")` |
| N79 | Eq. (8.42): $\mathrm{Re}(\mathbf u^*\cdot\nabla^*\mathbf u^*)=-\nabla^*p^*+\nabla^{*2}\mathbf u^*$ [#106] | 8.6 | B | NOTE | C12 | step of D26 |
| C12 | Eq. (8.43): Stokes (creeping-flow) equations $\nabla p=\mu\nabla^2\mathbf u$ — linear, pressure–viscous balance; bearing films, settling sediment, mist, molten plastic [#107] | 8.6 | A | CORE | – | load-bearing: the equations of every §8.6 result and of small-particle physics in Ch. 13 (settling, aerosols) and Ch. 16 (micro-organisms); the choice of pressure scale by dominant balance is a method used again in Ch. 9 and 13 |
| N80 | Rule: proper length and time scales depend on the region and are found by balancing the terms that matter there [#108] | 8.6 | C | NOTE | C12 | named in one sentence (the italic rule); pointers Ch. 9, Ch. 13 |
| N81 | Stokes (1851) near field: sphere at rest in a stream U; curl of (8.43) ⇒ ∇²ω = 0 (∇×∇p = 0); in spherical coordinates the operator means −∇×∇× (footnote) [#109] | 8.6 | B | NOTE | C13 | D27 |
| R17 | Only ω_φ ≠ 0: $\omega_\varphi=\frac1r[\frac{\partial(ru_\theta)}{\partial r}-\frac{\partial u_r}{\partial\theta}]$ [#110] | 8.6 | C | RECAP | C13 | ch03/`core.curvilinear.curl("spherical")`; one sentence |
| R18 | (6.83) re-displayed: $u_r=\frac{1}{r^2\sin\theta}\frac{\partial\psi}{\partial\theta}$, $u_\theta=-\frac{1}{r\sin\theta}\frac{\partial\psi}{\partial r}$ [#111] | 8.6 | B | RECAP | C13 | ch06 C13, `core.potential.axisym_velocity_spherical` |
| N82 | Vorticity from ψ: $\omega_\varphi=-\frac{E^2\psi}{r\sin\theta}$ with the Stokes operator E² of (6.77) [#112] | 8.6 | B | NOTE | C13 | step of D28; `core.creeping.E2` |
| N83 | Eq. (8.44): $[\frac{\partial^2}{\partial r^2}+\frac{\sin\theta}{r^2}\frac{\partial}{\partial\theta}(\frac1{\sin\theta}\frac{\partial}{\partial\theta})]^2\psi=0$ — the square of E², not the biharmonic [#113] | 8.6 | B | NOTE | C13 | result of D28 (★★★); ⚠️ callout (R17: a biharmonic test must fail) |
| N84 | Eq. (8.45): $\psi(r=a,\theta)=0$ (no normal flow) [#114] | 8.6 | B | NOTE | C13 | step of D29 |
| N85 | Eq. (8.46): $\partial\psi(r=a,\theta)/\partial r=0$ (no slip) [#115] | 8.6 | B | NOTE | C13 | step of D29 |
| R19 | Eq. (8.47): $\psi(r\to\infty,\theta)=\tfrac12Ur^2\sin^2\theta$ (uniform stream (6.86)) [#116] | 8.6 | B | RECAP | C13 | ch06 C13 (6.86); step of D29 |
| N86 | Separable ψ = f(r)sin²θ; $f^{iv}-\frac{4f''}{r^2}+\frac{8f'}{r^3}-\frac{8f}{r^4}=0$; f = Ar⁴ + Br² + Cr + D/r; A = 0, B = U/2, C = −3Ua/4, D = Ua³/4 [#117] | 8.6 | B | NOTE | C13 | steps of D29 (Euler–Cauchy primer) |
| C13 | Eq. (8.48): Stokes' stream function $\psi=Ur^2\sin^2\theta(\frac12-\frac{3a}{4r}+\frac{a^3}{4r^3})$ (body frame, θ from the downstream axis) [#118] | 8.6 | A | CORE | – | load-bearing: the only exact 3-D viscous flow round a body in the book; source of the pressure, the stresses and the drag, and the reference for Oseen, settling (Ch. 13) and micro-swimmers (Ch. 16) |
| N87 | Eq. (8.49): $u_r=U\cos\theta(1-\frac{3a}{2r}+\frac{a^3}{2r^3})$, $u_\theta=-U\sin\theta(1-\frac{3a}{4r}-\frac{a^3}{4r^3})$ [#119] | 8.6 | B | NOTE | C13 | stated (§4c: one derivative each of (8.48) via (6.83)) with a sympy cell; `stokes_sphere_velocity` |
| N88 | Eq. (8.50): $p-p_\infty=-\frac{3\mu aU\cos\theta}{2r^2}$; +3μU/2a at the front (θ = π), −3μU/2a at the rear (θ = 0) (printed without the minus, R12) [#120] | 8.6 | B | NOTE | C14 | D30 (★★★, rule b — the book skips the integration); `stokes_sphere_pressure`; ⚠️ R12 callout |
| N89 | Surface stresses on r = a (ours): $\sigma_{rr}=-p$ (viscous normal stress vanishes), $\sigma_{r\theta}=-\frac{3\mu U}{2a}\sin\theta$ (Fig. 8.17 upper panel) [#121] | 8.6 | B | NOTE | C14 | steps of D31; `stokes_sphere_surface_stresses` |
| C14 | Eq. (8.51): Stokes drag $D=6\pi\mu aU$ — ⅓ pressure (2πμaU) + ⅔ skin friction (4πμaU); drag ∝ velocity (Exercise 8.35) [#122] | 8.6 | A | CORE | – | load-bearing: Stokes' law sets terminal velocities of cloud droplets, aerosols and sediment (Ch. 13), the low-Re branch C_D = 24/Re of the drag curve (ch04 C15) and micro-organism locomotion (Ch. 16) |
| N90 | Terminal velocity: $(4/3)\pi a^3g(\rho'-\rho)=6\pi\mu aU$ ⇒ U_t = 2(ρ′ − ρ)ga²/(9μ) [#123] | 8.6 | B | NOTE | C14 | stated (§4c) with a number: 10 µm water droplet in air → U_t ≈ 1.2 cm/s, Re ≈ 0.016, D ≈ 4.1 × 10⁻¹¹ N; `terminal_velocity` (warns if Re > 0.1) |
| N91 | Millikan's oil drop (Fig. 8.18): $6\pi\mu U_ua+(4/3)\pi a^3g(\rho'-\rho)=neE$, E = −V_b/L; radius from the fall speed; charges are multiples of e [#124] | 8.6 | B | NOTE | C14 | stated with the synthetic-drop demo (seeded, recovers e within 1 %, V5 CODATA); `millikan_charge` |
| N92 | Eq. (8.52): $C_D=\frac{D}{\tfrac12\rho U^2\pi a^2}=\frac{24}{\mathrm{Re}}$, Re = 2aU/ν (diameter) [#125] | 8.6 | B | NOTE | C14 | stated (§4c, one line of algebra); recap ch04 C15 `sphere_drag_coefficient(model="stokes")` parity |
| R20 | Dimensional argument without ρ: D = f(μ, U, a) ⇒ D/μUa = const ⇒ C_D ∝ 1/Re (Exercise 4.60) [#126] | 8.6 | B | RECAP | C14 | ch01 C64–C69 Π method; `ch01.pi_groups` on (D, μ, U, a) |
| N93 | Fluid frame (sphere moving): $\psi=Ur^2\sin^2\theta(-\frac{3a}{4r}+\frac{a^3}{4r^3})$; fore–aft symmetric streamlines, no wake — linearity: reversing U flips u and p − p∞ (Fig. 8.19; "(9.63)" means (8.43), R14) [#127] | 8.6 | B | NOTE | C13 | stated (§4c) with the fluid-frame figure next to the ideal-flow sphere (`core.potential.sphere`); pointer Ch. 16 (reversibility) |
| C15 | Far-field breakdown of Stokes flow: viscous force ~ μUa/r³, inertia ~ ρU²a/r², ratio $\sim\frac{\rho Ua}{\mu}\frac ra=\mathrm{Re}\,\frac ra$ — inertia matters beyond r/a ~ 1/Re; Oseen's remedy [#128] | 8.6 | A | CORE | – | load-bearing: the first non-uniform (singular-perturbation) limit of the book — why a "small" term can win far away; motivates Oseen's linearisation, matched asymptotics (Ch. 9) and the corrected drag law |
| N94 | Not uniformly valid: the O(Re) correction is unbounded relative to Stokes (Whitehead); Stokes' paradox for a cylinder (Exercise 8.37); the non-uniformity is at infinity, not at the wall [#129] | 8.6 | C | NOTE | C15 | named; pointer Ch. 9 (singular perturbations) and Exercise 8.37 |
| N95 | Oseen (1910): u = U + u′, u·∇u ≈ U∂u′/∂x; $\rho U\frac{\partial u_i'}{\partial x}=-\frac{\partial p}{\partial x_i}+\mu\nabla^2u_i'$ (printed without the minus, R13); uniformly valid at lowest order [#130] | 8.6 | B | NOTE | C15 | stated (§4c: one substitution, quadratic terms dropped) + `ch08.oseen_linearisation_sympy`; ⚠️ R13 |
| N96 | Oseen boundary conditions: u′ → 0 far away; u′ = −U, v′ = w′ = 0 on r = a [#131] | 8.6 | B | NOTE | C15 | stated |
| N97 | Eq. (8.53): $\frac{\psi}{Ua^2}=[\frac{r^2}{2a^2}+\frac{a}{4r}]\sin^2\theta-\frac{3}{\mathrm{Re}}(1+\cos\theta)\{1-\exp[-\frac{\mathrm{Re}}{4}\frac ra(1-\cos\theta)]\}$, Re = 2aU/ν; no slip only to O(Re) [#132] | 8.6 | B | NOTE | C15 | stated (the book gives no derivation) with `oseen_streamfunction` (−expm1 form) |
| N98 | Near the surface a series of the exponential gives (8.48) ("(9.68)" means (8.48), R14) [#133] | 8.6 | B | NOTE | C15 | D33 |
| N99 | Oseen drag $C_D=\frac{24}{\mathrm{Re}}(1+\frac{3}{16}\mathrm{Re})$; experiments lie between Stokes and Oseen for Re < 5 [#134] | 8.6 | B | NOTE | C15 | stated (needs far-field matching — beyond this chapter) with the C_D(Re) figure (Stokes, Oseen, Proudman–Pearson, Morrison correlation band, V5); ⚠️ radius vs diameter Re (3/8 ↔ 3/16) |
| N100 | Oseen streamlines (Fig. 8.20, fluid frame): asymmetric, a wake behind the sphere [#135] | 8.6 | B | NOTE | C15 | figure at Re = 0.5, 1, 2 next to Stokes; mode in `stokes_sphere_flow` |
| N101 | Matched asymptotic expansions rationalise Oseen (Kaplun; Proudman & Pearson 1957); Proudman–Pearson drag $D=6\pi\mu aU(1+\tfrac38\mathrm{Re}_a+\tfrac9{40}\mathrm{Re}_a^2\ln\mathrm{Re}_a+\ldots)$ (not printed) [#136] | 8.6 | C | NOTE | C15 | named, with its curve on the C_D figure; pointer Ch. 9 |
| N102 | Final remarks: most analytically tractable laminar problems are solved; perturbation methods and numerical solution of NS carry the field on [#137] | 8.7 | C | NOTE | C15 | named in the closing paragraph; pointers Ch. 9 (perturbation) and Ch. 10 (CFD) |
| S01 | Exercises 8.1–8.38; literature cited and supplemental reading (named only; text never reproduced) [#138] | Ex. | C | SKIP | C15 | pointer: used ideas — 8.7 → N16; 8.12 → R11; 8.19–8.20 → N32/D13; 8.26 → N62; 8.30 → N56; 8.31 → R14; 8.34 → N36; 8.35 → D31; 8.37 → N94; the rest left to the reader |
| N103 | Fig. 8.1 Reynolds's apparatus (laminar vs turbulent dye) [#139] | 8.1 | B | NOTE | C01 | replaced by our synthetic dye-streak sketch driven by `ch08.pipe_flow_regime` (labelled schematic) |
| N104 | Figs. 8.5–8.7: pipe parabola with linear τ; annulus; rotating cylinder with the free-surface dip $z_s=-\Gamma^2/(8\pi^2gR^2)$ (ours, ch05 (5.7)) [#140] | 8.2 | B | NOTE | C04 | our remakes (`scripts/ch08_pipe_and_couette.py`) |
| N105 | Figs. 8.9–8.11: bearing pad, Hele-Shaw cell with a disc, spreading bead [#141] | 8.3 | B | NOTE | C07 | our remakes (`ch08_slider_bearing.py`, `ch08_hele_shaw.py`, `ch08_spreading.py`) in C07, C06 and C08 |
| N106 | Figs. 8.14–8.15: thickening vortex sheet (ω at two times, u in η); line-vortex decay at several νt [#142] | 8.4 | B | NOTE | C10 | our remakes (`ch08_vortex_diffusion.py`) |
| N107 | Figs. 8.16–8.20: oscillating-plate profiles at four phases; sphere surface stresses and pressure; Millikan apparatus; Stokes (fluid frame) and Oseen streamlines [#143] | 8.5–8.6 | B | NOTE | C13 | our remakes (`ch08_oscillating_plate.py` in C11, `ch08_stokes_sphere.py` in C13–C14, `ch08_millikan.py`, `ch08_oseen.py` in C15) |

### 2a. What each A block contains (for the lesson-designer)
Numbers are ours (water μ = 1.0 × 10⁻³ Pa s, ρ = 1000 kg/m³, ν = 10⁻⁶ m²/s unless stated); book-quoted values stay
in `tests/book_values_ch08.json`.
- **C01** picture: honey vs water poured down a straw; Reynolds's dye streak staying straight or smearing · question:
  "when is a flow smooth, and what does ν actually do?" · D01 · number: d = 1 cm, U = 0.1 m/s in water → Re = 1000
  (laminar); diffusion across 1 cm: water 100 s, air ≈ 6.7 s · code: `ch08.pipe_flow_regime`, `momentum_diffusivity`,
  `diffusion_time`, `inertia_viscous_scales` + from-scratch Re/L²/ν · figure: synthetic dye streak (laminar vs
  "turbulent" schematic, labelled) and diffusion time vs L on log axes for air and water · B/C: N01, N02, N103,
  R01–R06.
- **C02** picture: a belt pulling oil along a channel while a pump pushes back · question: "what profile results from a
  moving wall plus a pressure gradient, and when does the fluid near the fixed wall flow backwards?" · D02, D03, D04,
  D05 · number: h = 1 cm, U = 0.1 m/s → backflow once dp/dx > 2 Pa/m; Q(dp/dx = 0) = 5 × 10⁻⁴ m²/s · code:
  `core.laminar.channel_flow`, `channel_flow_rate`, `channel_shear_stress`, `channel_backflow_threshold`,
  `ch08.parallel_flow_sympy("channel")` + from-scratch finite-difference solve of μu″ = dp/dx · figure: Fig. 8.4 remake
  (four cases, backflow shaded) + τ(y) · B/C: N03–N11, R07, R08.
- **C03** picture: blood in a capillary, water in a garden hose · question: "how much flows for a given pressure drop,
  and why does halving the radius cut the flow sixteen-fold?" · D06, D07 · number: a = 1 mm, dp/dz = −1000 Pa/m →
  V = 0.125 m/s, u_max = 0.25 m/s, Q = 3.93 × 10⁻⁷ m³/s, τ₀ = −0.5 Pa, Re = 250, f = 0.256 · code:
  `core.laminar.pipe_poiseuille`, `pipe_flow_rate`, `pipe_wall_stress`, `pipe_friction_factor` + from-scratch
  trapezoid ∫u 2πR dR · figure: parabola with linear τ (Fig. 8.5 remake) + Q vs a on log–log (slope 4) · B/C: N12–N16,
  R09.
- **C04** picture: tea stirred between two cylinders; a Taylor–Couette lab rig · question: "what swirl is set up between
  two rotating cylinders, and what happens when one goes away?" · D08, D09 · number: R₁ = 1 cm, R₂ = 2 cm, Ω₁ = 1
  rad/s, Ω₂ = 0 → A = −1/3 s⁻¹, B = 1.33 × 10⁻⁴ m²/s, u_φ(1.5 cm) = 3.9 mm/s · code: `core.laminar.circular_couette`,
  `circular_couette_pressure`, `circular_couette_power` + from-scratch 2 × 2 solve · figure: u_φ(R) for four (Ω₁, Ω₂)
  pairs with the AR and B/R parts dashed, the R₂ → ∞ and R₁ → 0 ghosts · B/C: N17–N22, N104, R10–R12.
- **C05** picture: an oil film under a sliding block, 1000× thinner than it is long · question: "why can we throw away
  inertia even when UL/ν is large?" · D10, D11 · number: h = 50 µm, L = 5 cm, U = 5 m/s, ν = 10⁻⁴ → ε = 10⁻³,
  Re_L = 2500, ε²Re_L = 2.5 × 10⁻³; μUL/h² = 5 MPa (μ = 0.05) ≈ 50 atm · code: `core.lubrication.lubrication_scales`,
  `ch08.lubrication_nondim_sympy` + from-scratch coefficient table · figure: term-magnitude bars of (8.16a) vs ε on
  log axes (inertia falls as ε²Re_L) · B/C: N23, N24, N26–N31, R13.
- **C06** picture: the gap under a tilted pad, with Couette-plus-Poiseuille profiles at several stations · question:
  "given the gap shape, what is the velocity at every point and how much flows?" · D12, D13 · number: h = 50 µm,
  U_h = 5 m/s, dp/dx = 0 → q = Uh/2 = 1.25 × 10⁻⁴ m²/s; the dp/dx that halves q · code:
  `core.lubrication.lubrication_velocity`, `lubrication_flux`, `reynolds_pressure_1d`, `hele_shaw_velocity` · figure:
  Fig. 8.8 remake (profiles in a converging gap) + Hele-Shaw dye lines round a disc vs ideal-flow streamlines · B/C:
  N25, N32–N34, N36.
- **C07** picture: a tilted pad skating on oil carrying tonnes · question: "how does a thin viscous film carry a load,
  and why only in one direction?" · D14 · number: L = 5 cm, h₀ = 50 µm, U = 5 m/s, μ = 0.05 Pa s: α = 0.1 → W = 12.5
  (linear) vs 10.8 kN/m (exact); α = 0.5 → 62.5 vs 32.8 kN/m; U < 0 → W < 0 · code: `core.lubrication.slider_bearing`
  (exact / linear / book ghost), `slider_bearing_load`, `slider_optimum_taper`, `ch08.slider_bearing_sympy` +
  from-scratch cumulative-trapezoid pressure · figure: p(x) for several α with the printed-slip ghost + W(α) with the
  optimum · B/C: N35, N105.
- **C08** picture: honey (or lava) spreading on a table · question: "what equation governs a thin layer that spreads
  under its own weight?" · D15 · number: glycerol μ = 1 Pa s, ρ = 1260, a 1 cm high bead: flux coefficient ρg/3μ ≈
  4.1 × 10³ m⁻¹s⁻¹; front after 10 s, 100 s, 1000 s from the similarity law (C10) · code:
  `core.lubrication.thin_film_flux`, `thin_film_spread` (volume conserved to round-off) + from-scratch explicit
  finite-volume steps · figure: h(x, t) snapshots with the volume printed + animation A3 · B/C: N37.
- **C09** picture: a plate yanked sideways under still water; how far does the "news" of motion reach? · question:
  "what is u(y, t), and why does the whole family of profiles look the same when rescaled?" · D16, D17, D18, D19, D20 ·
  number: ν = 10⁻⁶, t = 100 s → √(νt) = 1 cm, δ₉₉ = 3.64 cm; u(1 cm, 100 s)/U = erfc(0.5) = 0.48 · code:
  `core.laminar.stokes_first_problem`, `similarity_variable`, `diffusion_thickness`, `stokes_first_vorticity`,
  `core.diffusion.crank_nicolson_1d` + from-scratch FTCS/CN march agreeing with (8.30) · figure: profiles at five t raw
  and collapsed (N52, N55), δ₉₉(t) on log axes (slope ½) · B/C: N38–N56, R14.
- **C10** picture: a thickening shear layer between two streams; a vortex dying; a spreading bead — three problems, one
  trick · question: "how do you find the exponents that make a problem self-similar?" · D21, D22, D23 · number: vortex
  sheet with U = 1 cm/s, ν = 10⁻⁶: width 5.54√(νt) = 5.5 mm at t = 1 s; C_f at Re_x = 10⁴: 0.0113 (temporal) vs 0.0066
  (Blasius, Ch. 9) · code: `ch08.similarity_reduce_sympy` (cases stokes1_delta, vortex_sheet, line_vortex, spreading),
  `core.laminar.vortex_sheet_diffusion`, `temporal_bl_wall_stress`, `line_vortex_decay`,
  `core.lubrication.viscous_current_similarity` + from-scratch log–log fit of the width (slope 0.5) and of the bead front
  (slope 0.2) · figure: Figs. 8.14–8.15 remakes + collapse of the bead h t^{1/5} vs x/t^{1/5} · B/C: N57–N64, N106.
- **C11** picture: a plate shaken back and forth under water; the tide rubbing the sea floor · question: "how deep does
  an oscillation reach, and why do deeper layers lag?" · D24, D25 · number: ω = 2π rad/s → √(2ν/ω) = 0.56 mm,
  4√(ν/ω) = 1.6 mm, e^{−4/√2} = 0.059; M₂ tide with ν = 10⁻² m²/s (eddy) → 12 m · code:
  `core.laminar.stokes_second_problem`, `stokes_layer`, `ch08.stokes_second_sympy` + from-scratch complex arithmetic ·
  figure: profiles at ωt = 0, π/2, π, 3π/2 with the ±e^{−y/δ} envelope + animation A2 · B/C: N65–N74.
- **C12** picture: a bacterium swimming, a grain of silt settling · question: "which terms survive when Re → 0, and why
  does pressure not disappear with inertia?" · D26 · number: a 10 µm droplet at 1.2 cm/s in air: Re ≈ 0.016 · code:
  `ch08.low_re_scaling_sympy`, `core.creeping.stokes_residual` · figure: coefficient table/bars for the two pressure
  scales vs Re (IF7) · B/C: N75–N80, R15, R16.
- **C13** picture: a tiny bead sinking in syrup; streamlines that look the same in front and behind · question: "what is
  the flow round a sphere when viscosity dominates everywhere?" · D27, D28, D29 · number: a = 1 mm, U = 1 mm/s: velocity
  at r = 2a on the side (θ = π/2) = U(1 − 3/8 − 1/32) = 0.59U (ideal flow: 1.06U) · code:
  `core.creeping.stokes_sphere_streamfunction`, `stokes_sphere_velocity`, `E2`, `E4_residual`, `stokes_sphere_sympy`
  + from-scratch finite-difference velocities from ψ · figure: body-frame streamlines + fluid-frame streamlines next to
  the ideal-flow sphere; speed deficit reaching tens of radii · B/C: N81–N87, N93, N107, R17–R19.
- **C14** picture: a cloud droplet falling at its terminal speed; Millikan's drops · question: "what force does the
  fluid exert, and where does it come from — pressure or friction?" · D30, D31 · number: 10 µm water droplet in air
  (μ = 1.81 × 10⁻⁵) → U_t = 1.20 cm/s, D = 4.1 × 10⁻¹¹ N (⅓ pressure, ⅔ friction), C_D = 24/Re = 1.5 × 10³ · code:
  `core.creeping.stokes_sphere_pressure`, `stokes_sphere_surface_stresses`, `stokes_drag(parts=True)`,
  `sphere_drag_quadrature`, `terminal_velocity`, `stokes_drag_coefficient` + from-scratch midpoint sum of the traction ·
  figure: surface p and σ_rθ vs θ (Fig. 8.17 remake) + U_t vs radius for droplets, sand, bacteria with the Re = 0.1 limit
  marked · B/C: N88–N92, R20.
- **C15** picture: far behind a slowly falling bead, the fluid "notices" inertia · question: "if Re is tiny, how can
  inertia ever matter?" · D32, D33 · number: Re (radius) = 0.01 → inertia ≈ viscosity at r ≈ 100a; C_D(Re = 0.5):
  Stokes 48, Oseen 52.5 · code: `core.creeping.inertia_viscous_ratio`, `oseen_streamfunction`, `oseen_drag_coefficient`,
  `ch08.oseen_limit_sympy` + from-scratch ratio at points by finite differences · figure: ratio vs r/a on log–log (slope
  1) for three Re; Stokes vs Oseen streamlines (Fig. 8.20); C_D(Re) comparison · B/C: N94–N102.

## 3. Section coverage
| § | Title | A | B | C | RECAP | SKIP |
|---|---|---|---|---|---|---|
| 8.1 | Introduction | C01 | N01, N103 | N02 | R01 (B), R02 (B), R03 (B), R04 (B), R05 (B), R06 (C) | — |
| 8.2 | Exact Solutions for Steady Incompressible Viscous Flow | C02, C03, C04 | N03, N04, N05, N06, N07, N08, N09, N10, N12, N13, N14, N15, N16, N17, N18, N19, N20, N21, N22, N104 | N11 | R07 (B), R08 (B), R09 (C), R10 (B), R11 (B), R12 (B) | — |
| 8.3 | Elementary Lubrication Theory | C05, C06, C07, C08 | N23, N24, N25, N26, N27, N28, N29, N30, N31, N32, N33, N34, N35, N36, N105 | N37 | R13 (B) | — |
| 8.4 | Similarity Solutions for Unsteady Incompressible Viscous Flow | C09, C10 | N38, N39, N40, N41, N42, N43, N44, N45, N46, N47, N48, N49, N50, N51, N52, N53, N54, N55, N56, N57, N58, N59, N60, N61, N62, N63, N106 | N64 | R14 (B) | — |
| 8.5 | Flow Due to an Oscillating Plate | C11 | N65, N66, N67, N68, N69, N70, N71, N72, N73 | N74 | — | — |
| 8.6 | Low Reynolds Number Viscous Flow Past a Sphere | C12, C13, C14, C15 | N77, N78, N79, N81, N82, N83, N84, N85, N86, N87, N88, N89, N90, N91, N92, N93, N95, N96, N97, N98, N99, N100, N107 | N75, N76, N80, N94, N101 | R15 (B), R16 (B), R17 (C), R18 (B), R19 (B), R20 (B) | — |
| 8.7 | Final Remarks | — | — | N102 | — | S01 (end of chapter) |

Totals: A 15 · B 113 (96 NOTE + 17 RECAP) · C 15 (11 NOTE + 3 RECAP + 1 SKIP) = 143. No section is empty: §8.7 has no A
item and is covered by N102 (C) in the closing paragraph of C15 (§8.6); S01 (exercises) is listed under §8.7 as the end
of the chapter. N107 (Figs. 8.16–8.20) spans §8.5–8.6 and is listed under §8.6; its Fig. 8.16 part is drawn in C11.

## 4. Prerequisites needing primers (concept or tool | needed by | why it is not A/B/RECAP)
Rows marked **gloss** are one plain sentence where the item is used; all others are full 📎 primers (one or two plain
sentences, what it means here, a 2–4-line runnable demo with easy numbers). New primers continue at **P185**. Items
already in `knowledge/primers.md` get a one-line reminder naming the earlier primer (listed at the end). Placement: the
cylindrical-Laplacian and Euler-ODE primers go at the top of C03/C04 (before D06, D08); the anisotropic-scaling primer
before D10; the Gaussian-integral and erf/erfinv primers before D19–D20; the exponent-matching primer at the top of C10
(before D21); the E²-operator primer before D28.

| Concept or tool | Needed by | Why it is not A/B/RECAP |
|---|---|---|
| diffusivity: a quantity with units m²/s whose square root of (D·t) is how far something spreads in time t; momentum, heat and dye all diffuse | C01 (D01), C09 | P130 order-of-magnitude reminder; "ν as a diffusivity" is framed as a primer so C09's √(νt) is familiar |
| integrating an ODE twice with two constants and fixing them from two boundary values (a 2 × 2 linear system) | C02 (D04), C03 (D06), C04 (D08), C06 (D12) | P44/P57 reminders — **gloss** with a 3-line sympy `dsolve` demo |
| Laplacian of an axisymmetric function in cylindrical coordinates $\frac1R\frac{d}{dR}(R\frac{du}{dR})$ and the φ-component of the vector Laplacian $\frac{d}{dR}[\frac1R\frac{d}{dR}(Ru_\varphi)]$ | C03 (D06), C04 (D08), N61 | P88/P105 primed cylindrical unit vectors and the moving basis; the operator itself is new (from `core.curvilinear`) |
| regularity at the axis: ln R → −∞ as R → 0, so a bounded solution drops it | C03 (D06) | **gloss** (P36 log reminder) |
| area element of a disc 2πR dR and a flux integral over a circular cross-section | C03 (D07) | P83/P163 reminder — **gloss** |
| Euler–Cauchy (equidimensional) ODE: a trial $R^\lambda$ turns it into a polynomial in λ | C04 (D08), C13 (D29) | new maths tool |
| anisotropic scaling: two length scales for one flow (L along, h across), and choosing the velocity scale of v from continuity | C05 (D10), C12 (D26) | P133 primed scaled variables with one length; two scales and the ε bookkeeping are new |
| bookkeeping of a small parameter: multiply an equation through so one chosen term has coefficient 1, then read each other term's power of ε | C05 (D10, D11) | P68 orders-of-smallness reminder — **gloss** |
| Leibniz rule for an integral whose upper limit moves, $\frac{\partial}{\partial x}\int_0^{h(x)}u\,dy=\int_0^h\frac{\partial u}{\partial x}dy+u(h)\frac{\partial h}{\partial x}$ | C06 (D13) | P109 primed differentiation under the integral sign with fixed limits; the moving-limit term is new |
| kinematic condition at a moving wall v = ∂h/∂t + u ∂h/∂x | C06 (D13) | P132 moving level set reminder — **gloss** |
| antiderivative of $(1+\alpha x/L)^{-n}$ and partial fractions; Taylor expansion in a small parameter α | C07 (D14), N35 | P106 substitution and P26 Taylor reminders — **gloss** |
| `scipy.integrate.cumulative_trapezoid`: running integral of samples | C06 (code), C07 (from-scratch) | new Python tool |
| stress-free surface μ∂u/∂y = 0 and hydrostatic pressure in a thin layer | C08 (D15) | P20 boundary conditions reminder — **gloss** |
| nonlinear diffusion: a diffusion equation whose diffusivity depends on the unknown (here ∝ h³), conservative (flux) form and why volume is conserved | C08 (D15, code) | new maths/physics idea |
| implicit time stepping with Picard iteration (backward Euler, solve, repeat until the diffusivity stops changing); precursor film | C08 (code) | new numerics idea (P30 explicit stepping reminder) |
| Crank–Nicolson for the diffusion equation and a tridiagonal solve with `scipy.linalg.solve_banded` | C09 (code, from-scratch), C11 (V3 check) | new numerics/Python tool (P21 finite differences reminder) |
| order of a PDE vs number of conditions (two in y, one in t) | C09 (N42) | **gloss** |
| Gaussian integral $\int_{-\infty}^\infty e^{-\zeta^2}d\zeta=\sqrt\pi$ | C09 (D19), C10 (D22) | new maths tool |
| error function erf, its complement and inverses: `scipy.special.erf`, `erfc`, `erfinv`, `erfcinv` | C09 (D19, D20), C10 (D22) | P123 primed erfc; erf as an integral and the inverses are new |
| `scipy.integrate.solve_bvp` for a two-point boundary-value problem on a truncated domain | C09 (check of (8.26)), C10 | new Python tool |
| exponent matching: two sides of an equation that must hold for all t force the powers of t to agree | C10 (D21, D22, D23) | new maths move (P43 exponent rules reminder) |
| spotting an exact derivative F + ηF′ = (ηF)′ | C10 (D22) | P38 product rule reminder — **gloss** |
| complex exponential solution of a periodic problem, $u=\mathrm{Re}\{e^{i\omega t}f(y)\}$, and √i = (1 + i)/√2 | C11 (D24) | P176 complex amplitudes and P159 complex square roots reminders — **gloss** |
| dominant balance: the scale of a quantity (here pressure) is set by the term it must balance | C12 (D26) | new modelling idea |
| the Stokes operator E² for axisymmetric flow and operator composition (E²)²ψ ≠ ∇⁴ψ | C13 (D28, D29) | ch06 (6.77) used E² as a check only; applying it twice and the spherical curl-curl identity are new |
| curl of a gradient is zero and curl commutes with the Cartesian Laplacian (index notation) | C13 (D27) | P121/P122 reminders — **gloss** |
| traction on a sphere and its x-projection, $t_x=\sigma_{rr}\cos\theta-\sigma_{r\theta}\sin\theta$, integrated with dA = 2πa² sin θ dθ | C14 (D31) | P163 surface integrals on a sphere reminder + ch02 Cauchy traction recap — **gloss** |
| line integral of a gradient field to recover a potential (here p from ∇p) and checking the second component | C14 (D30) | P35 line integral reminder — **gloss** |
| asymptotic size of each term as r → ∞ (keep the slowest-decaying part) | C15 (D32) | P130 reminder — **gloss** |
| series 1 − e^{−s} = s − s²/2 + … for small s | C15 (D33) | P26/P117 reminder — **gloss** |

Reminders only (already primed): partial derivative P25 · first-order Taylor P26 · definite integral P27 · separation of
variables P42 · linear second-order ODE P44 · complex numbers P45 · chain rule P49/P91 · exponent rules P43 · power laws
and log–log plots P13 · orders of smallness P68 · `np.linalg.solve` P57 · sympy P40 · sympy expand/series P117 · `quad`
P87 · `brentq` P108 · `np.expm1` P107 · substitution in an integral P106 · differentiation under the integral sign P109 ·
fundamental theorem of calculus P84 · curl of a curl P122 · Schwarz P121 · scaled variables P133 · order-of-magnitude
scaling P130 · erfc P123 · reduced gravity/buoyancy P131 · frames of reference P96 · complex amplitudes P176 · complex
square roots P159 · separation of variables for a PDE P167 · surface integrals on a sphere P163 · Gauss–Legendre P143 ·
2-D Dirac delta P149 · `meshgrid` P76 · broadcasting P77 · contour/streamplot P78 · `animate` P16 · `slider_figure` P17 ·
`show_viz` P18 · live widgets P47 · `minimize_scalar` P170 · `np.interp` P182.

## 4b. Derivations written out (parsed by tools: ID first, CORE id in a column, ★★★ for hard, explainer slugs backticked in the LAST column)
Written out because (a) the result belongs to an A item, or (b) the book never writes it out and the lesson needs it
(marked "(b)"). One small move per step; the book's skipped moves (analysis §2b) are filled in and the notebook says so.
Analysis §2b numbers in brackets (`a-Dnn`). 33 rows; 296 steps.

| ID | Result (Eq.) | CORE | Difficulty | Steps | Tools used | Traps | Shown in |
|---|---|---|---|---|---|---|---|
| D01 | air diffuses momentum ≈ 15× faster than water although water is ≈ 55× more viscous: ν = μ/ρ and t = L²/ν [a-D02] | C01 | ★ | 5 | diffusivity (primer), ratios | comparing μ instead of ν; the density ratio (≈ 830) beats the viscosity ratio (≈ 55); units m²/s | notebook |
| D02 | v ≡ 0 and the reduced equations (8.4a) $0=-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{d^2u}{dy^2}$, (8.4b) from (8.1) and continuity [a-D03] | C02 | ★ | 7 | continuity (R03), partial derivative (P25) | ∂v/∂y = 0 only says v is constant across the gap — the wall value fixes it to 0; u∂u/∂x + v∂u/∂y vanishes **exactly** (not approximately), which is why the problem is linear | notebook · `couette_poiseuille_backflow` |
| D03 | dp/dx is a constant: "f(x) = g(y) ⇒ both constant" [a-D04] | C02 | ★ | 5 | partial derivative (P25) | differentiate each side by the *other* variable; the constant can have either sign; p = p(x) was needed first (from (8.4b)) | notebook · `couette_poiseuille_backflow` |
| D04 | Couette–Poiseuille profile (8.5) $u=\frac Uh y-\frac1{2\mu}\frac{dp}{dx}y(h-y)$ [a-D05] | C02 | ★★ | 9 | integrating twice (gloss), 2 × 2 system (P57) | the book's A, B are minus our c₁, c₂; y² − hy = −y(h − y) produces the printed form; dp/dx < 0 gives a *positive* parabola; check both walls at the end | notebook · `couette_poiseuille_backflow` |
| D05 | Couette and Poiseuille limits, τ(y), ∣τ_w∣ = (h/2)∣dp/dx∣ and the backflow threshold $\frac{dp}{dx}>\frac{2\mu U}{h^2}$ (b) [a-D07] | C02 | ★ | 7 | derivative as a slope (P19), sign of a derivative | backflow starts when u′(0) < 0, not when the maximum moves; the Poiseuille maximum is at y = h/2 and equals 1.5V; τ is signed (the book draws ∣τ∣) | notebook · `couette_poiseuille_backflow` |
| D06 | Poiseuille pipe profile (8.6) $u_z=\frac{R^2-a^2}{4\mu}\frac{dp}{dz}$ [a-D08] | C03 | ★★ | 10 | cylindrical Laplacian (primer), regularity at the axis (gloss), integrating twice | forgetting the R in R du/dR before the second integration; the ln R term dropped by boundedness (equivalently du/dR = 0 on the axis); factor 1/(4μ), not 1/(2μ) | notebook |
| D07 | (8.7), wall stress (8.8) as a CV balance, Hagen–Poiseuille Q, V, u_max = 2V and f = 64/Re (b: CV balance and f are ours) [a-D09] | C03 | ★★ | 10 | area element 2πR dR (gloss), CV momentum (ch04 C04 recap) | τ₀ = (a/2)dp/dz is negative for forward flow — sign of the traction on the fluid; ∫(R² − a²)2πR dR = −πa⁴/2; Re with the diameter and mean velocity; f is the Darcy factor 8τ₀/ρV² | notebook |
| D08 | circular Couette (8.9) $u_\varphi=AR+B/R$ and (8.10) with its constants [a-D10] | C04 | ★★ | 10 | φ-component of the cylindrical vector Laplacian (primer), Euler–Cauchy ODE (primer), 2 × 2 system | the φ-equation has no pressure (∂p/∂φ = 0) and no advection; the first integral is (1/R)(Ru_φ)′ = 2A (the 2 is absorbed); u_φ at the wall is ΩR, not Ω | notebook |
| D09 | limits (8.11) $u_\varphi=\Omega_1R_1^2/R$ and (8.12) $u_\varphi=\Omega_2R$ from (8.10) [a-D11] | C04 | ★ | 6 | limits (P68) | divide numerator and denominator by R₂² before letting R₂ → ∞; Ω₂ = 0 is needed or A → Ω₂; Γ = 2πΩ₁R₁² (ch05's ω = 2Ω₁) | notebook |
| D10 | scaled continuity (8.15) and momentum (8.16a), (8.16b) from (6.2), (8.13a,b) with the scalings (8.14) [a-D14] | C05 | ★★★ | 14 | anisotropic scaling (primer), scaled variables and the chain rule (P133), ε bookkeeping (gloss) | v scales with εU (from continuity), not U; x-momentum is multiplied by ε²L²/(νU), y-momentum by ε³L²/(νU) — different factors so the pressure keeps 1/Λ; the printed ∂p/∂x in (8.13b) must be ∂p/∂y (R6) | notebook · `lubrication_scaling` |
| D11 | lubrication balance (8.17a) $0\cong-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{\partial^2u}{\partial y^2}$, (8.17b) for ε²Re_L ≪ 1 [a-D15] | C05 | ★★ | 7 | orders of smallness (P68) | dropping inertia needs ε²Re_L ≪ 1, not Re_L ≪ 1; ∂p*/∂y* = O(ε²), not 0; re-dimensionalising restores the ν the book prints without (R7); with P_a, Λ ≈ 10³ — use μUL/h² (R24) | notebook · `lubrication_scaling` |
| D12 | (8.18) → lubrication profile (8.19) [a-D16] | C06 | ★★ | 8 | integrating with "constants" A(x, t), B(x, t) (gloss) | ∂p/∂x is constant only in the y-integration; B = U₀ and A = (U_h − U₀)/h − (h/2μ)∂p/∂x; the printed U₀ placement fails u(h) = U_h when U₀ ≠ 0 (R8b) — code the consistent form | notebook · `slider_bearing` |
| D13 | gap flux $q=-\frac{h^3}{12\mu}\frac{\partial p}{\partial x}+\frac{(U_0+U_h)h}{2}$ and the 1-D Reynolds equation $\frac{\partial h}{\partial t}+\frac{\partial q}{\partial x}=0$ (b) [a-D17] | C06 | ★★★ | 12 | Leibniz rule with a moving limit (primer), kinematic wall condition (gloss), polynomial integrals | the Leibniz boundary term −u(h)∂h/∂x cancels against the kinematic v(h); v(0) = 0 only for a flat, impermeable lower wall; steady ⇒ q constant (Example 8.1's C₁ is q − Uh in the pad frame) | notebook · `slider_bearing` |
| D14 | slider bearing: C₁, exact $p-p_e=\frac{6\mu LU}{h_o^2}\frac{\alpha(x/L)(1-x/L)}{(2+\alpha)(1+\alpha x/L)^2}$, linear-α pressure and $W=\frac{\alpha\mu L^2U}{2h_o^2}$ [a-D18] | C07 | ★★★ | 15 | moving CV (ch04 C01 recap), antiderivatives of (1 + αx/L)^{−n} (gloss), 2 × 2 system, Taylor in α (P26) | pad frame makes the flow steady; the integrands are (1 + αx/L), not the printed (1 − αx/L) (R8); the final denominator is squared, not first power (R8) — checked by substituting into dp/dx = −12μC₁/h³ − 6μU/h²; α ≪ 1 used only at the end; W < 0 if αU < 0 | notebook · `slider_bearing` |
| D15 | thin-film equation $\frac{\partial h}{\partial t}=\frac{\rho g}{3\mu}\frac{\partial}{\partial x}(h^3\frac{\partial h}{\partial x})$ (Example 8.3) [a-D21] | C08 | ★★ | 10 | CV mass (ch04 C01 recap), hydrostatics, stress-free surface (gloss), polynomial integral | hydrostatic p needs a small slope (lubrication), not just "no acceleration" (R23); the top condition is ∂u/∂y = 0, not u = 0; ∫₀^h y(2h − y)dy = 2h³/3; the sign: fluid flows down the slope | notebook · `viscous_gravity_current` |
| D16 | diffusion equation (8.20) $\frac{\partial u}{\partial t}=\nu\frac{\partial^2u}{\partial y^2}$ for the impulsively started plate [a-D22] | C09 | ★ | 6 | continuity (D02's move), partial derivative | ∂p/∂x = 0 because a still region exists far away at every finite t; the equation is linear because u∂u/∂x = 0, not because u is small | notebook · `stokes_first_problem` |
| D17 | similarity variable (8.24) → (8.25) $u/U=F(y/\sqrt{\nu t})$ from dimensional analysis and linearity [a-D23] | C09 | ★★ | 8 | Π groups (ch01 recap), linearity of a solution (gloss) | 5 variables and 2 dimensions give 3 groups; linearity (u ∝ U) removes y/Ut, not dimensional analysis alone; other choices of η change F, not the answer | notebook · `stokes_first_problem` |
| D18 | similarity ODE (8.26) $-\frac\eta2F'=F''$ with (8.27) F(0) = 1 and (8.28) F(∞) = 0 [a-D24 + a-D25] | C09 | ★★ | 9 | chain rule (P49/P91) | ∂η/∂t = −η/(2t) (not −η/t); dividing by U/t must remove every t — if it does not, the similarity form is wrong; t → 0 at fixed y and y → ∞ at fixed t both give η → ∞, so three conditions become two | notebook · `stokes_first_problem` |
| D19 | (8.29) and the erf solution (8.30) $\frac uU=1-\mathrm{erf}(\frac{y}{2\sqrt{\nu t}})$ [a-D26] | C09 | ★★ | 11 | separation of variables (P42), Gaussian integral (primer), substitution (P106), erf (primer) | the substitution ξ = 2ζ carries dξ = 2dζ; ∫₀^∞e^{−ζ²}dζ = √π/2 (half the full integral); the factor 2 in y/(2√(νt)) comes from that substitution; A < 0 | notebook · `stokes_first_problem` |
| D20 | 99 % thickness (8.31) $\delta_{99}=2\,\mathrm{erfc}^{-1}(0.01)\sqrt{\nu t}=3.64\sqrt{\nu t}$ (b) [a-D27] | C09 | ★ | 5 | inverse erfc (primer), `brentq` (P108) | the figure's axis is η/2 = 1.82, the formula's η = 3.64 (R22); the "99 %" level is a convention — a 95 % level gives 2.77 | notebook · `stokes_first_problem` |
| D21 | Example 8.4: the ansatz with A = 1, n = 0 gives δδ′ = C₁ν and δ = √(2C₁νt) [a-D29] | C10 | ★★ | 8 | exponent matching (primer), chain rule, separable ODE (P42) | A = 1, n = 0 because u = U at η = 0 for every t; ∂F(y/δ)/∂t = −ηF′δ′/δ; δ(0) = 0 sets the integration constant; C₁ = ½ recovers η of (8.25) | notebook · `similarity_exponents` |
| D22 | Example 8.5: n = ½, $F=De^{-\eta^2/4}$, $\omega_z=-\frac{U}{\sqrt{\pi\nu t}}e^{-y^2/4\nu t}$, $u=U\,\mathrm{erf}(\frac{y}{2\sqrt{\nu t}})$, ±0.95U at η = ±2.772 [a-D30] | C10 | ★★★ | 14 | exponent matching (primer), integral constraint, exact derivative (gloss), Gaussian integral, erf | the jump ∫ω dy = −2U must hold for all t — that fixes n = ½, not the ODE; F + ηF′ = (ηF)′; the integration constant C = 0 because F, ηF → 0; AD = −U/√(πν) (sign); 2.772, not the printed 2.76 (R11) | notebook · `similarity_exponents` |
| D23 | Example 8.7: n = m = 1/5 for the spreading bead, $h=At^{-1/5}F(x/Dt^{1/5})$ [a-D34] | C10 | ★★ | 9 | exponent matching (primer), product and chain rules | (h³h_x)_x brings δ^{−2} and A⁴t^{−4n}; the second condition is the conserved volume (−n + m = 0), not a boundary condition; the page's "Example 8.2" means 8.3 (R14) | notebook · `viscous_gravity_current` · `similarity_exponents` |
| D24 | Stokes' second problem (8.35) → (8.36) → k → (8.37) → (8.38) [a-D35] | C11 | ★★ | 10 | complex exponential form (gloss, P176), √i (P159), linear ODE e^{λy} (P44) | ∂/∂t of e^{iωt}f is iωe^{iωt}f; of the two roots only the decaying one is bounded (B = 0); the real part of e^{iωt}e^{−(1+i)y/δ} is e^{−y/δ}cos(ωt − y/δ) — the phase lag has a minus sign; the text's "(8.33)" means (8.35) (R14) | notebook · `oscillating_plate` |
| D25 | reading (8.38): amplitude 0.059U at y = 4√(ν/ω), crest speed √(2νω), wavelength 2π√(2ν/ω), e-folding depth √(2ν/ω) (b) [a-D36] | C11 | ★ | 6 | phase of a cosine (P165), exponent rules | the book's δ ~ 4√(ν/ω) is 2√2 times the e-folding depth √(2ν/ω); the "wave" is diffusion with a phase lag — no restoring force (R16) | notebook · `oscillating_plate` |
| D26 | (8.39) → (8.40) → (8.41) → (8.42) → Stokes equations (8.43) $\nabla p=\mu\nabla^2\mathbf u$ [a-D37] | C12 | ★★ | 9 | scaled variables (P133), dominant balance (primer) | multiplying (8.40) by Re and letting Re → 0 kills the pressure too (the wrong equation 0 = μ∇²u); the pressure scale must change to μU/L; re-dimensionalise at the end | notebook |
| D27 | curl of the Stokes equations: ∇²ω = 0 [a-D38] | C13 | ★★ | 6 | index notation (ch02 recap), curl of a gradient = 0 and curl-curl (P121, P122 gloss) | curl and the Laplacian commute only in Cartesian components; in spherical coordinates "∇²ω" means −∇×∇×ω (footnote) | notebook |
| D28 | the stream-function equation (8.44) $(E^2)^2\psi=0$ from ω_φ = −E²ψ/(r sin θ) [a-D39] | C13 | ★★★ | 12 | Stokes stream function (6.83) (R18), E² operator (primer), spherical curl (`core.curvilinear`), sympy | it is the square of E², not the biharmonic ∇⁴ (a sympy test of ∇⁴ψ must fail, R17); θ from the downstream axis; −∇×∇×(ω_φe_φ) = 0 ⇔ E²(r sin θ ω_φ) = 0 needs the sympy identity cell | notebook · `stokes_sphere_flow` |
| D29 | f-ODE, f = Ar⁴ + Br² + Cr + D/r and Stokes' ψ (8.48) from (8.44)–(8.47) [a-D40] | C13 | ★★ | 10 | Euler–Cauchy ODE (primer), 2 × 2 system | E²(f sin²θ) = (f″ − 2f/r²)sin²θ; the roots are {4, 2, 1, −1}; A = 0 from the far field (r⁴ would outgrow the stream); C = −3Ua/4, D = Ua³/4 from (8.45)–(8.46) | notebook · `stokes_sphere_flow` |
| D30 | pressure (8.50) $p-p_\infty=-\frac{3\mu aU\cos\theta}{2r^2}$ by integrating ∇p = μ∇²u (b) [a-D42] | C14 | ★★★ | 11 | spherical vector Laplacian (`core.curvilinear`), vorticity route −μ∇×ω, line integral of a gradient (gloss) | compute μ∇²u = −μ∇×ω with ω_φ = −(3Ua/2r²)sin θ; integrate ∂p/∂r then **check** (1/r)∂p/∂θ; p → p∞ fixes the constant; the rear minimum is −3μU/2a (printed without the minus, R12) | notebook · `stokes_drag_settling` |
| D31 | Stokes drag (8.51) D = 6πμaU: 2πμaU pressure + 4πμaU friction (Exercise 8.35) [a-D43] | C14 | ★★★ | 13 | Newtonian stress in spherical coordinates (ch04 C07 recap), traction projection (gloss), surface integral dA = 2πa² sin θ dθ (P163) | the viscous normal stress 2μ∂u_r/∂r vanishes on r = a; project the traction on e_x with θ from the downstream axis (signs!); the pressure integral is ⅓, friction ⅔; a slip sphere (bubble) would give 4πμaU — not this case | notebook · `stokes_drag_settling` |
| D32 | far-field breakdown: inertia/viscous ~ (ρUa/μ)(r/a) [a-D47] | C15 | ★★ | 7 | asymptotic size of terms (gloss, P130) | far away the disturbance decays like Ua/r (the 3Ua/4r term), not like the a³/r³ term; ∂/∂r brings a 1/r; Re here uses the radius | notebook · `stokes_sphere_flow` |
| D33 | Oseen's (8.53) reduces to Stokes' (8.48) near the sphere [a-D49] | C15 | ★★ | 7 | series of 1 − e^{−s} (gloss), sin²θ = (1 − cos θ)(1 + cos θ) | s = (Re/4)(r/a)(1 − cos θ) is small only where Re r/a ≪ 1; (3/Re)(1 + cos θ)s = (3/4)(r/a)sin²θ; Oseen meets no slip only to O(Re) (R18) | notebook · `stokes_sphere_flow` |

## 4c. Derivations demoted to statements (first column is the A parent in bold, e.g. **C20** — never a bare ID or a D id: the parser reads a bare first-cell ID as an item and blanks its tier)
Results **given, not derived**, inside the named A block (a paragraph, the equation, a number, a check in code). 15 rows.

| A parent | Analysis §2b item | Result stated (Eq.) | Stated in (B item) | Why not written out |
|---|---|---|---|---|
| **C01** | a-D01 2-D vorticity diffusion | $D\omega_z/Dt=\nu\nabla^2\omega_z$ | R02 | SEEN: ch05 C06 (5.13); one sentence on why (ω·∇)u = 0 in plane flow |
| **C02** | a-D06 flow rate and mean velocity | $Q=\frac{Uh}{2}[1-\frac{h^2}{6\mu U}\frac{dp}{dx}]$, $V=Q/h$ | N10 | two polynomial integrals; `channel_flow_rate` vs `quad`; the R9 slip shown by a units check |
| **C04** | a-D12 σ_Rφ, zero net force, power = dissipation (Exercise 8.12) | $\sigma_{R\varphi}=-2\mu\Omega_1R_1^2/R^2$, power $=4\pi\mu\Omega_1^2R_1^2$ | R10, R11 | SEEN: ch05 D03, N08; `circular_couette_power` asserts power in = ∫ε dA (V4) |
| **C04** | a-D13 circular Couette pressure (ours) | $p=p_1+\rho[\frac{A^2}{2}(R^2-R_1^2)+2AB\ln\frac R{R_1}-\frac{B^2}2(\frac1{R^2}-\frac1{R_1^2})]$ | N20 | one antiderivative; checked by finite differences of dp/dR = ρu²/R |
| **C06** | a-D20 Hele-Shaw ⇒ Laplace | $\bar{\mathbf u}=-\frac{h^2}{12\mu}\nabla p$, $\nabla^2p=0$ | N36 | the 2-D Reynolds equation with constant h (D13's move in 2-D); checked symbolically (∇²φ = 0) and by the disc figure |
| **C07** | a-D19 exact load and optimum taper (ours) | $W=\frac{6\mu UL^2}{h_o^2\alpha^2}[\ln(1+\alpha)-\frac{2\alpha}{2+\alpha}]$, 1 + α ≈ 2.19 | N35 | partial fractions repeating D14's integrals; `quad` check and San Andrés benchmark (V5) |
| **C09** | a-D28 vorticity content and no new vorticity | $\int_0^\infty\omega\,dy=U$ | N53 | one application of the fundamental theorem of calculus; `vorticity_content` = U; the printed −U named (R10) |
| **C10** | a-D31 temporal boundary layer | $\tau_w=\frac{\mu U}{\sqrt{\pi\nu t}}$, $C_f=\frac2{\sqrt\pi}\sqrt{\frac\nu{U^2t}}$ | N60 | one derivative of D22's u and a Galilean shift; `temporal_bl_wall_stress`; Ch. 9 compares with Blasius |
| **C10** | a-D32 line-vortex decay (Example 8.6, ★★★) | $F=1-e^{-\eta^2/4}$, $u_\theta=\frac{\Gamma}{2\pi r}[1-e^{-r^2/4\nu t}]$ | N61 | a recap result (ch03 C14 Gaussian vortex, ch05 N13); same ansatz moves as D22 with (8.32b); a sympy residual cell of the cylindrical diffusion equation is shown |
| **C10** | a-D33 line-vortex spin-up (Exercise 8.26) | $u_\theta=\frac{\Gamma}{2\pi r}e^{-r^2/4\nu t}$ | N62 | an exercise; substitution check (sympy residual 0) |
| **C13** | a-D41 velocities (8.49) | $u_r=U\cos\theta(1-\frac{3a}{2r}+\frac{a^3}{2r^3})$, $u_\theta=-U\sin\theta(1-\frac{3a}{4r}-\frac{a^3}{4r^3})$ | N87 | one derivative each via (6.83); sympy cell + no-slip and far-field checks |
| **C13** | a-D46 fluid-frame ψ and reversibility | $\psi=Ur^2\sin^2\theta(-\frac{3a}{4r}+\frac{a^3}{4r^3})$ | N93 | subtract ½Ur²sin²θ; symmetry asserted numerically |
| **C14** | a-D44 terminal velocity and Millikan | $U_t=\frac{2(\rho'-\rho)ga^2}{9\mu}$, $ne=\frac{6\pi\mu a(U_u+U_t)}{E}$ | N90, N91 | force balances in one line each; `terminal_velocity`, `millikan_charge` |
| **C14** | a-D45 drag coefficient and the ρ-free Π argument | $C_D=24/\mathrm{Re}$, D/(μUa) = const | N92, R20 | one line of algebra; parity with ch04 `sphere_drag_coefficient("stokes")` |
| **C15** | a-D48 Oseen linearisation | $\rho U\frac{\partial u_i'}{\partial x}=-\frac{\partial p}{\partial x_i}+\mu\nabla^2u_i'$ | N95 | one substitution with the quadratic terms dropped; `oseen_linearisation_sympy`; the printed sign slip named (R13) |

## 5. Interactive explainers (5–10 + backup)
Nine explainers on the A items where manipulation or motion teaches most; one backup. Each has the required live
**Explain** tab ("Explanation & interpretation", numbered sections computing every number on screen with the reader's
settings, ending in "Reading the current setting"), a synced **Code** tab, a **Derivation** tab for its D ids, a
4–8-step walkthrough, ≥ 3 check questions and ≥ 2 more depth features. Physics mirrors scalar-callable `fluidpy`
functions (parity rows). **Colours** (consistent across explainers, notebook figures and derivation terms): pressure /
pressure-driven part orange, viscous / wall-driven part rose, inertia teal, velocity profile blue, vorticity purple
(`accent`), backflow and warnings amber, similarity/rescaled curves purple dashed, analytic ghost muted dashed, numerical
dots teal. Linked, not duplicated: ch01 `viscosity_momentum_diffusion` (from E5: Couette start-up, the bounded-gap
contrast), ch04 `navier_stokes_term_balance` (from E1 and E2: which terms vanish), ch05 `vortex_pressure_funnel` (from the
backup: the R₂ → ∞ vortex), ch06 `added_mass_sphere` / `superposition_sandbox` (from E8: the ideal-flow sphere).
C01, C03, C04 and C12 are taught with static figures and slider figures (§6): a regime number, a parabola, a two-cylinder
profile and a coefficient table do not gain enough from manipulation to justify an explainer; C04 is the backup.

### E1 · couette_poiseuille_backflow
- **A:** C02 (also shows N03–N10 reduction steps, N09 backflow threshold, N10 Q and V, R07 Couette, R08 Poiseuille, N11)
- **Confusion removed:** "a favourable pressure gradient and a moving wall — which one wins, and how can fluid flow
  *against* the wall?" — the profile is the plain sum of a straight line and a parabola; backflow appears exactly when the
  parabola's slope at the fixed wall beats the line's, dp/dx > 2μU/h².
- **Why interactive:** dragging dp/dx through zero and past 2μU/h² makes the linear and parabolic parts add up in front
  of you, the backflow pocket open at the fixed wall and the flow rate cross zero — a family of curves that a static
  figure can only sample at four values (Fig. 8.4).
- **Stage:** (1) the channel: moving top wall, tracer particles advected at u(y), a dye line bending into the profile,
  backflow pocket shaded amber; (2) u(y) with the linear (rose) and parabolic (orange) parts dashed and their sum (blue);
  (3) (hidePortrait) τ(y) and Q as a bar with the Couette and Poiseuille contributions.
- **Controls:** wall speed U (−1 … 1 m/s) · pressure gradient dp/dx (signed, in units of 2μU/h² or Pa/m) · gap h ·
  viscosity μ (optional) · transport (time: particles and dye).
- **Presets:** pure Couette · pure Poiseuille (U = 0) · favourable · onset of backflow (dp/dx = 2μU/h²) · strong adverse ·
  zero net flow (Q = 0: dp/dx = 6μU/h²).
- **Equations shown (live):** (8.4a) $0=-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{d^2u}{dy^2}$; (8.5)
  $u=\frac Uh y-\frac1{2\mu}\frac{dp}{dx}y(h-y)$; $Q=\frac{Uh}2[1-\frac{h^2}{6\mu U}\frac{dp}{dx}]$; $\tau=\mu du/dy$;
  backflow $\frac{dp}{dx}>\frac{2\mu U}{h^2}$.
- **Mirrors:** `core.laminar.channel_flow`, `channel_flow_rate`, `channel_shear_stress`, `channel_backflow_threshold`,
  `couette_poiseuille_state` (§8).
- **Derivations:** D02 (reduction: the view dims the v-arrows to zero), D03 (dp/dx constant: the pressure bar along x
  goes linear), D04 (the profile: step "apply u(h) = U" lights the top wall), D05 (backflow: preset at the threshold).
- **Depth features:** Explain tab (the two parts at the probe height → u; τ at both walls; Q and V (with the R9 note);
  threshold and how far past it you are → reading the current setting), synced Code tab, + linked views (3), transport,
  presets, status ("➡️ all forward" / "⚠️ backflow near the fixed wall" / "⏸ zero net flow"), terms (Q = Couette part +
  Poiseuille part, click to isolate), inspector (click a height: u arithmetic).
- **Follows reference:** `forced_damped_vibrations.html` (system + graph linked on one state, live explanation panel
  that computes every number) with the term bars of `fid_formula_lab.html`.
- **Aha:** the profile is literally a line plus a parabola; fluid near the fixed wall reverses as soon as the adverse
  gradient's slope there exceeds the wall's pull — dp/dx > 2μU/h².

### E2 · lubrication_scaling
- **A:** C05 (also shows N23, N24 (R6), N26–N31, R13; the Λ remark R24)
- **Confusion removed:** "Re_L = UL/ν is in the thousands — how can inertia be negligible?" — because in a thin gap
  inertia is weighed against ε⁻² times the viscous term: the relevant number is ε²Re_L, and pressure is set by μUL/h².
- **Why interactive:** the reader drags h/L and U and watches the term-by-term bars of (8.16a) and (8.16b) (each term's
  coefficient × O(1)) reorder on log axes; the moment ε²Re_L ~ 1 the inertia bar catches the viscous one — an ordering
  argument you can see, which a static table cannot convey.
- **Stage:** (1) the gap drawn to scale with a zoom (h exaggerated by a stated factor), the profile and v-arrows ~ εU;
  (2) term bars for x-momentum (inertia teal, pressure orange, streamwise diffusion rose-light, cross-gap diffusion rose)
  on a log scale; (3) (hidePortrait) the same for y-momentum, showing only pressure without an ε.
- **Controls:** fineness ratio ε = h/L (log 10⁻⁴ … 0.3) · Re_L (log 1 … 10⁶) or (U, ν) chips · pressure scale chips: P_a
  (book) / μUL/h² (natural) · fluid presets.
- **Presets:** engine bearing oil film · Hele-Shaw cell · human knee (synovial fluid) · a "thick" gap where lubrication
  fails (ε = 0.2, Re_L = 10⁴).
- **Equations shown (live):** (8.14) scalings; (8.15); (8.16a) and (8.16b) with live coefficients; (8.17a)
  $0\cong-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{\partial^2u}{\partial y^2}$; (8.17b).
- **Mirrors:** `core.lubrication.lubrication_scales`, `lubrication_term_magnitudes` (§8).
- **Derivations:** D10 (★★★ — each substitution step highlights the term bar it creates), D11 (the limit: preset with
  ε²Re_L = 10⁻³).
- **Depth features:** Explain tab (ε, Re_L, ε²Re_L, Λ with P_a and with μUL/h², each term's size; verdict → reading the
  current setting), synced Code tab, + terms (bars that reorder), presets, status ("✅ lubrication valid: ε²Re_L = …" /
  "⚠️ inertia matters"), modes (x-momentum / y-momentum emphasis), linked views (2–3).
- **Follows reference:** `fid_formula_lab.html` (clickable formula terms, term-by-term bars, stage chips) with the live
  numbered explanation of `amplitude_phase_second_order_II_3.html`.
- **Aha:** a thin gap divides inertia by ε² — only ε²Re_L must be small, and the pressure that results is μUL/h², hundreds
  of atmospheres in a real bearing.

### E3 · slider_bearing
- **A:** C07, C06 (also shows N25, N32 (Reynolds equation), N33, N34 (U₀ slip), N35 exact load and optimum, N105)
- **Confusion removed:** "how does a film of oil a twentieth of a millimetre thick hold up a heavy pad — and why does
  the tilt direction matter?" — continuity forces the same flux through a narrowing gap, which needs a pressure hump; the
  hump reverses into suction if the pad moves the other way.
- **Why interactive:** dragging the taper α and flipping the sliding direction reshapes the pressure hump, the gap
  profiles (Couette → Couette + Poiseuille, with backflow at the inlet at large taper) and the load in one motion, and the
  printed-slip ghost visibly fails the end condition — the whole of Example 8.1 in one picture.
- **Stage:** (1) the pad and the gap (exaggerated), velocity profiles at 5 stations from (8.19), streamlines/tracers in
  the pad frame; (2) p(x) — exact (orange), linear-α (dashed), the printed first-power form (muted ghost), max marked;
  (3) (hidePortrait) W(α) curve (exact vs linear) with the optimum 1 + α ≈ 2.19 and the current dot.
- **Controls:** taper α (−0.9 … 3) · sliding speed U (signed) · minimum gap h₀ (µm) · viscosity μ · model chips (exact /
  linear / book ghost).
- **Presets:** small taper (α = 0.1: linear ≈ exact) · optimum taper · reversed motion (W < 0) · parallel plates (α = 0: no
  load) · steep taper (inlet backflow).
- **Equations shown (live):** (8.19) profile; gap flux $q=-\frac{h^3}{12\mu}\frac{dp}{dx}+\frac{Uh}2$; the steady
  Reynolds equation q = const; $p-p_e=\frac{6\mu LU}{h_o^2}\frac{\alpha(x/L)(1-x/L)}{(2+\alpha)(1+\alpha x/L)^2}$;
  $W=\frac{\alpha\mu L^2U}{2h_o^2}$ and the exact W(α).
- **Mirrors:** `core.lubrication.slider_bearing`, `slider_bearing_load`, `slider_optimum_taper`, `lubrication_velocity`,
  `slider_bearing_state` (§8).
- **Derivations:** D12 (profile), D13 (★★★ Reynolds equation: the "Leibniz term" step highlights the moving top wall),
  D14 (★★★ slider: the "p(L) = p_e" step shows the printed ghost missing it).
- **Depth features:** Explain tab (h(x), C₁ = q − Uh, p_max and where, load exact vs linear, % error, pressure in
  atmospheres; the book-slip note → reading the current setting), synced Code tab, + linked views (3), presets, status
  ("🏋️ load carried: W = … kN/m" / "⚠️ W < 0: the pad is sucked down" / "— no taper, no load"), inspector (click a
  station: u(y) arithmetic, q), modes (exact / linear / printed ghost).
- **Follows reference:** `amplitude_phase_second_order_II_3.html` (system + response curves + a numbered derivation
  with live numbers).
- **Aha:** the same flux must squeeze through a narrowing gap, so pressure must rise inside — thousands of atmospheres
  from a film thinner than a hair, but only if the pad slides toward the narrow end.

### E4 · viscous_gravity_current
- **A:** C08, C10 (the Example 8.7 part: N63; also shows N37, N64, N105)
- **Confusion removed:** "a spreading drop of honey slows down — why does its front go as t^{1/5} and not √t like
  diffusion?" — the diffusivity ρgh³/3μ collapses as the layer thins, and the conserved volume ties height to width.
- **Why interactive:** the reader plays the spreading from any initial shape and watches it forget that shape and lock
  onto the Huppert similarity profile, with the front marching along a t^{1/5} line on log–log axes while the volume
  readout stays constant — convergence to self-similarity is a process in time.
- **Stage:** (1) the layer h(x, t) (blue fill) with the flux arrows and the similarity profile ghost (purple dashed);
  (2) log–log front position x_N(t) with slope-1/5 and slope-1/2 reference lines; (3) (hidePortrait) rescaled profiles
  h t^{1/5} vs x/t^{1/5} collapsing.
- **Controls:** initial shape chips (box / Gaussian / two humps) · volume · viscosity μ (log: water … glycerol … lava) ·
  density · transport (time, log clock).
- **Presets:** honey on a plate · lava lobe (μ = 10³ Pa s) · a glacier-like slab (labelled analogy: real ice is
  non-Newtonian) · wrong exponent (the rescaled view with m = 1/2 fails to collapse).
- **Equations shown (live):** thin-film equation $\frac{\partial h}{\partial t}=\frac{\rho g}{3\mu}\frac{\partial}{\partial x}(h^3\frac{\partial h}{\partial x})$;
  flux $-\frac{\rho g}{3\mu}h^3h_x$; similarity $h=At^{-1/5}F(x/Dt^{1/5})$; x_N = η_N(βA³t)^{1/5}.
- **Mirrors:** `core.lubrication.thin_film_flux`, `viscous_current_similarity`, `thin_film_state` (§8); the JS explicit
  finite-volume march is checked against the similarity solution by parity rows at late time.
- **Derivations:** D15 (the thin-film equation), D23 (n = m = 1/5).
- **Depth features:** Explain tab (flux at the probe, volume, x_N from the law vs simulated, effective diffusivity at the
  centre → reading the current setting: "still remembers its initial shape" / "self-similar"), synced Code tab, + linked
  views (3), transport, presets, status (volume conserved to …; similarity reached ✓), inspector (click x: h, h_x, q).
- **Follows reference:** `angular_frequency_explorer_1.html` (system animation + graphs on one clock, presets,
  end-of-run summary card).
- **Aha:** a layer spreading under its own weight slows because its own thinning chokes the flux — volume conservation
  then fixes the t^{1/5} law, and every initial shape ends up on the same profile.

### E5 · stokes_first_problem
- **A:** C09 (also shows N38–N56, R14 contrast, the vorticity view N53)
- **Confusion removed:** "how can a whole family of profiles at different times be one curve?" — every profile is the
  same shape stretched by √(νt); plot against η = y/(2√(νt)) and they coincide; δ₉₉ grows like √t.
- **Why interactive:** the reader starts the plate, watches the profile spread with a moving δ₉₉ marker, then flips the
  raw/rescaled mode and sees all earlier profiles fall on one curve — and changing ν or U does not break the collapse.
- **Stage:** (1) the fluid above the plate with tracers and a dyed vertical line (moving right with the local speed),
  vorticity shading; (2) u(y) at the current time with faint earlier profiles (raw) or all profiles vs η (rescaled), CN
  numerical dots over the erf curve; (3) (hidePortrait) δ₉₉(t) on log–log axes, slope ½.
- **Controls:** plate speed U · viscosity ν (log: air, water, oil) · mode chips raw / rescaled · level (99 %, 95 %) ·
  transport (time, log clock).
- **Presets:** water · air (15× faster) · honey · "stop the plate at T" (N56: collapse breaks) · Couette start-up with a
  second wall (R14: collapse breaks once the layer reaches the wall).
- **Equations shown (live):** (8.20) $\frac{\partial u}{\partial t}=\nu\frac{\partial^2u}{\partial y^2}$; (8.25)
  $u/U=F(y/\sqrt{\nu t})$; (8.26); (8.30) $\frac uU=1-\mathrm{erf}(\frac{y}{2\sqrt{\nu t}})$; (8.31)
  $\delta_{99}\sim3.64\sqrt{\nu t}$; ∫ω dy = U.
- **Mirrors:** `core.laminar.stokes_first_problem`, `diffusion_thickness`, `stokes_first_vorticity`,
  `stokes_first_stopped`, `stokes_first_state` (§8); ch01 `couette_startup_profile` for the preset.
- **Derivations:** D16, D17, D18, D19 (the "ξ = 2ζ" step sets the rescaled mode), D20.
- **Depth features:** Explain tab (√(νt), η at the probe, u/U = erfc(η/2) with numbers, δ₉₉, wall stress, vorticity
  content = U → reading the current setting), synced Code tab, + linked views (3), transport, modes (raw / rescaled),
  presets, inspector (click a point: η and erfc arithmetic), status ("self-similar ✓" / "⚠️ an imposed scale breaks
  similarity").
- **Follows reference:** `forced_damped_vibrations.html` (explanation panel, linked views on one time slider) with the
  collapse mode pattern of `knowledge/viz_patterns.md` (boundary / similarity layers).
- **Aha:** there is no length scale in the problem except √(νt), so the profile can only stretch — one curve in η for
  every time, speed and viscosity.

### E6 · similarity_exponents
- **A:** C10 (also shows N57, N58, N59, N60, N61, N62, N64, N106)
- **Confusion removed:** "where do the exponents in a similarity solution come from?" — you guess the form
  At^{−n}F(ξ/δ(t)), and the equation plus one conserved quantity leave no freedom: n and the growth of δ are forced.
- **Why interactive:** the reader sets trial exponents n and m (δ ∝ t^m) with sliders; the rescaled profiles collapse only
  at the right values, and the bracket powers in the reduced equation turn green when they match — the exponent algebra
  becomes a game you can win or lose.
- **Stage:** (1) profiles at five times in raw coordinates; (2) the same rescaled with the trial n, m (collapse spread
  printed); (3) (hidePortrait) the reduced-equation bracket powers as chips, and the conserved integral vs t.
- **Controls:** mode chips (impulsive plate n = 0 · vortex sheet n = ½ · line vortex (8.32b) · spreading bead) · trial n ·
  trial m · reveal toggle.
- **Presets:** the correct exponents for each mode · a near miss (n = 0.4) · "diffusion guess" m = ½ for the bead.
- **Equations shown (live):** (8.32a) $\gamma=At^{-n}F(\xi/\delta(t))$; (8.32b); the bracketed equation of each example
  with live powers; ω_z and u of Example 8.5; F = 1 − e^{−η²/4} of Example 8.6; C_f of the temporal boundary layer.
- **Mirrors:** `core.laminar.vortex_sheet_diffusion`, `line_vortex_decay`, `stokes_first_problem`,
  `core.lubrication.viscous_current_similarity`, `ch08.similarity_collapse_error` (§8).
- **Derivations:** D21 (Example 8.4), D22 (★★★ Example 8.5: the "jump must be −2U for all t" step sets n = ½), D23
  (Example 8.7).
- **Depth features:** Explain tab (the powers of t in each bracket with the trial values, the conserved integral, collapse
  error → reading the current setting), synced Code tab, + modes (4 systems), linked views (3), presets, status ("✅
  collapsed: n = ½, δ ∝ t^½" / "❌ spread = …"), terms (bracket powers as chips).
- **Follows reference:** `stride_padding_playground.html` (a formula with numbers plugged in, badges that explain the
  result, presets at the classic settings) with the modes of `angular_frequency_explorer_1.html`.
- **Aha:** a similarity form is a guess with two free exponents; the equation fixes one relation and a conserved quantity
  fixes the other — nothing is left to choose.

### E7 · oscillating_plate
- **A:** C11 (also shows N65–N73; ch07's phase language)
- **Confusion removed:** "the oscillation looks like a wave travelling into the fluid — is it?" — it is diffusion driven
  by a periodic wall: each layer lags by y/δ radians and shrinks as e^{−y/δ}; the depth is set by √(ν/ω), so faster
  shaking reaches less far.
- **Why interactive:** the animation shows the phase lag and the envelope together; dragging ω shrinks the layer while
  the ghost envelope and the "0.06U at 4√(ν/ω)" marker move, and a probe clock shows a deep layer reaching its maximum
  later — none of which a four-phase static figure makes vivid.
- **Stage:** (1) the plate with tracers at several heights oscillating with growing lag, a dyed vertical line bending;
  (2) u(y, t) with the ±Ue^{−y/δ} envelope and the four book phases as faint ghosts; (3) (hidePortrait) u(t) at the wall
  and at a probe height with the lag marked.
- **Controls:** frequency ω (log) · viscosity ν (log: molecular water / eddy viscosity for the ocean) · probe height y ·
  transport (time).
- **Presets:** lab plate 1 Hz in water · M₂ tide over the sea floor (molecular vs eddy ν) · acoustic boundary layer in air
  (1 kHz) · Ekman teaser (labelled: same (1 + i)/δ, Ch. 13).
- **Equations shown (live):** (8.33) $u(0,t)=U\cos\omega t$; (8.35) $u=\mathrm{Re}\{e^{i\omega t}f(y)\}$; (8.36)
  $i\omega f=\nu f''$; $k=\pm(1+i)\sqrt{\omega/2\nu}$; (8.38)
  $u=Ue^{-y\sqrt{\omega/2\nu}}\cos(\omega t-y\sqrt{\omega/2\nu})$; δ ~ 4√(ν/ω) and √(2ν/ω).
- **Mirrors:** `core.laminar.stokes_second_problem`, `stokes_layer`, `stokes_layer_state` (§8).
- **Derivations:** D24 (the "choose the bounded root" step shows the growing root as a ghost), D25.
- **Depth features:** Explain tab (δ_e, the book's 4√(ν/ω), amplitude and lag at the probe, crest speed and wavelength,
  period → reading the current setting), synced Code tab, + linked views (3), transport, presets, inspector (click a
  height: amplitude and phase arithmetic), notes ("Right now" with a table of real Stokes-layer depths, current row
  highlighted).
- **Follows reference:** `amplitude_phase_second_order_II_3.html` (amplitude and phase lag as the story; the numbered
  "evaluate → amplitude → phase → time lag → at time t" explanation).
- **Aha:** each layer is the wall's motion, delayed by y/δ radians and shrunk by e^{−y/δ} — a "wave" made only of
  diffusion, reaching √(2ν/ω).

### E8 · stokes_sphere_flow
- **A:** C13, C15 (also shows N81–N87, N93 fluid frame, N94–N100 Oseen, R17–R19; ideal-flow sphere from ch06)
- **Confusion removed:** "why does the creeping flow look the same in front of and behind the sphere, why does it
  disturb the fluid so far away, and where does that picture fail?" — linearity makes it reversible (no wake); the
  disturbance decays only like a/r; far enough away (r/a ~ 1/Re) inertia returns and Oseen's wake appears.
- **Why interactive:** toggling body/fluid frame, Stokes/Oseen/ideal and dragging Re shows symmetric loops, a slowly
  decaying speed deficit and then a wake growing from far downstream, with the r/a = 1/Re circle marking where inertia
  catches up — three comparisons that need motion and switching to be believed.
- **Stage:** (1) streamlines and tracer particles round the sphere (body or fluid frame), speed heatmap, the r/a ~ 1/Re
  circle; (2) speed along the side line (θ = π/2) vs r/a on log axes: Stokes, ideal, Oseen; (3) (hidePortrait)
  inertia/viscous ratio vs r/a (slope 1) with the current Re.
- **Controls:** model chips (Stokes / Oseen / ideal flow) · frame chips (body / fluid) · Re (log 10⁻³ … 5) · view zoom
  (r/a up to 100) · transport (tracers).
- **Presets:** Stokes at Re → 0 · fluid frame (loops that return) · Oseen Re = 1 (wake) · ideal flow for contrast.
- **Equations shown (live):** (8.43) $\nabla p=\mu\nabla^2\mathbf u$; (8.44) $(E^2)^2\psi=0$; (8.48)
  $\psi=Ur^2\sin^2\theta(\frac12-\frac{3a}{4r}+\frac{a^3}{4r^3})$; (8.49); inertia/viscous ~ Re r/a; (8.53).
- **Mirrors:** `core.creeping.stokes_sphere_streamfunction`, `stokes_sphere_velocity`, `oseen_streamfunction`,
  `oseen_velocity` (§8), `inertia_viscous_ratio`; `core.potential.sphere`.
- **Derivations:** D28 (★★★ the E⁴ equation), D29 (f(r) and ψ: the "A = 0" step zooms out to show r⁴ blowing up), D32
  (breakdown: the ratio view), D33 (Oseen → Stokes: Re slider toward 0).
- **Depth features:** Explain tab (ψ and velocity at the probe with numbers, speed deficit at 10a vs ideal flow, the ratio
  Re r/a and where it reaches 1 → reading the current setting), synced Code tab, + linked views (3), modes (Stokes /
  Oseen / ideal), transport, presets, inspector (click a point: u_r, u_θ arithmetic), status ("↔️ fore–aft symmetric" /
  "🌊 Oseen wake: inertia wins beyond r ≈ … a").
- **Follows reference:** `angular_frequency_explorer_1.html` (modes = the same object in different models, linked views,
  presets, "Right now" notes).
- **Aha:** creeping flow is reversible and reaches far — its disturbance decays only like a/r, which is exactly why,
  beyond r ~ a/Re, inertia comes back and a wake appears.

### E9 · stokes_drag_settling
- **A:** C14 (also shows N88 pressure (R12), N89 stresses, N90 terminal velocity, N91 Millikan, N92 C_D, R20, N99 Oseen C_D)
- **Confusion removed:** "where does 6πμaU come from — the front pushing, or the sides rubbing?" — pressure contributes
  one third and friction two thirds; and "how fast does a droplet fall?" — the drag grows with speed until it equals the
  effective weight, U_t ∝ a².
- **Why interactive:** the reader sweeps around the sphere and sees the local pressure and shear tractions as arrows and
  as curves vs θ, the running integral of each building to 2πμaU and 4πμaU in the term bars; then drags the particle
  radius and watches U_t and Re, with the Stokes law's validity limit crossing — the integration and the size law are
  both processes.
- **Stage:** (1) the sphere with traction arrows (pressure orange, friction rose) on its surface and a sweep marker;
  (2) p and σ_rθ vs θ (Fig. 8.17 remake) with the running x-integrals; (3) (hidePortrait) settling: U_t vs radius on
  log–log for droplets / sand / bacteria with Re = 0.1 marked, or C_D(Re) (Stokes, Oseen, correlation) as a mode.
- **Controls:** sweep angle θ (transport) · viscosity μ · stream speed U · particle radius a (log) · particle chips (cloud
  droplet in air / quartz sand in water / bacterium / oil drop) · view mode (settling / C_D).
- **Presets:** cloud droplet 10 µm · drizzle 100 µm (Re > 1: Stokes overestimates U_t) · silt in a river · Millikan drop.
- **Equations shown (live):** (8.50) $p-p_\infty=-\frac{3\mu aU\cos\theta}{2r^2}$; σ_rθ on r = a; (8.51) $D=6\pi\mu aU$;
  U_t = 2(ρ′ − ρ)ga²/(9μ); (8.52) $C_D=24/\mathrm{Re}$; Oseen C_D.
- **Mirrors:** `core.creeping.stokes_sphere_pressure`, `stokes_sphere_surface_stresses`, `stokes_drag`,
  `terminal_velocity`, `stokes_drag_coefficient`, `oseen_drag_coefficient`, `settling_state` (§8);
  `core.similarity.sphere_drag_coefficient`.
- **Derivations:** D30 (★★★ pressure: the "check (1/r)∂p/∂θ" step shows both components agreeing), D31 (★★★ drag: the
  sweep integrates as the steps go).
- **Depth features:** Explain tab (p and σ_rθ at the swept angle, the x-traction, the running integrals, total D; U_t,
  Re, validity; C_D → reading the current setting), synced Code tab, + terms (pressure ⅓ / friction ⅔ bars that add to
  6πμaU), linked views (3), transport (sweep), presets, status ("✅ Stokes regime Re = …" / "⚠️ Re > 0.1: Stokes
  underestimates drag"), inspector (click the sphere: traction arithmetic), notes (table of real settling speeds, current
  row highlighted).
- **Follows reference:** `fid_formula_lab.html` (term bars adding to a total, click to isolate) with the live
  explanation of `forced_damped_vibrations.html`.
- **Aha:** a third of Stokes drag is the pressure pushing on the front, two thirds is friction on the sides; because drag
  ∝ U, a droplet's fall speed grows as its radius squared — 1 cm/s for a 10 µm cloud droplet.

### B1 · rotating_cylinders_couette (backup)
- **A:** C04 (also shows N17–N22, R10–R12, N104; D08, D09 in its Derivation tab if built)
- **Confusion removed:** "what swirl fills the gap between two rotating cylinders, and how does it become a free vortex
  or a rigid rotation?" — it is always AR + B/R; the walls decide A and B; removing the outer wall kills A, removing the
  inner one kills B.
- **Stage:** top view of the annulus with tracers moving at u_φ/R · u_φ(R) with the AR and B/R parts and the two limit
  ghosts · (hidePortrait) angular momentum (Ru_φ)² vs R as a Rayleigh-criterion teaser for Ch. 11.
- **Mirrors:** `core.laminar.circular_couette`, `circular_couette_pressure`, `circular_couette_power`,
  `circular_couette_state` (§8). **Follows:** `angular_frequency_explorer_1.html` (rotating system + curve on one clock).
- **Aha:** the gap flow is a mix of rigid rotation and a free vortex; drag the outer cylinder to infinity and only the
  vortex is left — viscous, yet irrotational.

## 6. Python animations and interactive figures (A ID → what, why, player/figure kind)
Animations (`fluidpy.core.anim`; ≤ 120 frames, dpi ≤ 80; FAST halves frames):

| Kind | A ID | What moves | Why | Player |
|---|---|---|---|---|
| A1 | **C09** | impulsively started plate: u(y, t) growing with the δ₉₉ marker, CN dots on the erf curve, a dyed line bending | the √t spreading is a process; seeing the numerical and analytic curves ride together builds trust | video |
| A2 | **C11** | oscillating plate: profile swinging inside the ±e^{−y/δ} envelope, tracers at several heights with visible lag | the phase lag is only believable in motion | video |
| A3 | **C08** | spreading bead from `thin_film_spread`: h(x, t) with the front marker and the volume printed; inset log–log x_N(t) | forgetting the initial shape and locking onto t^{1/5} happens in time | frames |
| A4 | **C13** | fluid-frame tracers passing a moving sphere (Stokes) next to the ideal-flow sphere: loops that return vs particles that barely move | reversibility and the far-reaching disturbance need motion | video |
| A5 | **C10** | vortex sheet thickening and line vortex decaying (two panels), with the rescaled collapse appearing at the end | the √(νt) growth of both, and the collapse, as one motion | video |

Interactive figures (`slider_figure` / `animate_figure`; ≤ 40 steps × ≤ 4 traces; survive publishing):

| Kind | A ID | Slider controls | What changes | Why |
|---|---|---|---|---|
| IF1 | **C02** | dp/dx (signed, in units of 2μU/h²) | u(y) with its linear and parabolic parts; backflow onset | the page-surviving twin of E1 |
| IF2 | **C03** | radius a (log) at fixed dp/dz | pipe profile and Q on a log–log inset (slope 4) | Q ∝ a⁴ is the surprise |
| IF3 | **C04** | Ω₂/Ω₁ (−2 … 2) at R₁/R₂ = 0.5 | u_φ(R) with the solid-body and free-vortex ghosts | co- vs counter-rotation without a kernel (backup's static twin) |
| IF4 | **C07** | taper α | p(x) exact, linear and the printed ghost; W readout in the title | the pressure hump reshapes with the taper |
| IF5 | **C09** | time t (`animate_figure`) | Stokes-first profiles raw (left) and rescaled (right) | collapse without a kernel |
| IF6 | **C11** | phase ωt (`animate_figure`) | Stokes-layer profile inside its envelope | page-surviving A2 |
| IF7 | **C12** | Re (log 10⁻³ … 10³) | coefficients of inertia, pressure, viscous terms in both pressure scalings (as markers on a log axis) | why the dynamic scale gives the wrong Re → 0 limit |
| IF8 | **C15** | Re (0.1 … 2) | Stokes vs Oseen streamlines (precomputed contour polylines) and C_D values | the wake grows with Re |

Live widget (kernel only, paired with IF4 and E3): `live(slider_bearing, h0=…, alpha=…, U=…, mu=…)` printing p_max, W
exact and linear, and the error (C07).

## 7. From-scratch moments
Each is a transparent hand-written version next to the tested function, followed by `assert np.allclose(mine, lib)`.

| Section | A ID | Hand-written version | Tested function it must match |
|---|---|---|---|
| §8.1 | **C01** | Re = Ud/ν and t = L²/ν by hand for air and water, the regime by `if` | `ch08.pipe_flow_regime`, `diffusion_time` (rtol 1e-12) |
| §8.2 | **C02** | finite-difference solve of μu″ = dp/dx with u(0) = 0, u(h) = U (tridiagonal `np.linalg.solve` on 101 points) | `core.laminar.channel_flow` (rtol 1e-10 — exact for a quadratic) |
| §8.2 | **C03** | trapezoid ∫₀^a u_z 2πR dR on 2001 points | `core.laminar.pipe_flow_rate` Q (rtol 1e-6) |
| §8.2 | **C04** | `np.linalg.solve` of the 2 × 2 wall conditions for A, B | `core.laminar.circular_couette(return_coeffs=True)` (rtol 1e-12) |
| §8.3 | **C05** | the coefficient table ε²Re_L, 1/Λ, ε², ε⁴ from the definitions | `core.lubrication.lubrication_scales` (rtol 1e-12) |
| §8.3 | **C07** | dp/dx = 12μ(q̄ − q)/h³ with q found by `brentq` so p(L) = p_e, then `cumulative_trapezoid` | `core.lubrication.slider_bearing(model="exact")` (rtol 1e-6) |
| §8.3 | **C08** | 200 explicit finite-volume steps (face h³ by averaging, stable dt) | `core.lubrication.thin_film_spread` at the same time (rtol 1e-2) and volume conserved (1e-12) |
| §8.4 | **C09** | Crank–Nicolson march by hand (tridiagonal `solve_banded`) for 200 steps | `core.laminar.stokes_first_problem` (rtol 1e-3, grid-limited) |
| §8.4 | **C10** | log–log fit of the sheet width (slope 0.5) and of the bead front x_N(t) (slope 0.2) | exponents from `ch08.similarity_reduce_sympy` (abs 5e-3) |
| §8.5 | **C11** | `np.real(U*np.exp(1j*w*t)*np.exp(-(1+1j)*y/de))` | `core.laminar.stokes_second_problem` (rtol 1e-12) |
| §8.6 | **C13** | central differences of ψ for u_r, u_θ via (6.83) | `core.creeping.stokes_sphere_velocity` (rtol 1e-6) |
| §8.6 | **C14** | midpoint sum over θ of the x-traction × 2πa² sin θ, separately for pressure and friction | `core.creeping.stokes_drag(parts=True)` (rtol 1e-6) |
| §8.6 | **C15** | finite-difference ∣u·∇u∣/∣ν∇²u∣ at a few points on the side line | `core.creeping.inertia_viscous_ratio` (rtol 1e-4) |

## 8. Notes for the implementer
- **Budget (< 5 min on Colab CPU):** `thin_film_spread` on 400 cells (FAST 200) to t = 10³ t* with implicit steps,
  cached to `outputs/ch08/thin_film_run.npz`; CN runs 400 points × 2000 steps (FAST 200 × 500); the Stokes-second V3
  check over 10 periods only in tests, not the notebook; Hele-Shaw Poisson grid 128² (FAST 64²); sympy cells
  (`parallel_flow_sympy`, `lubrication_nondim_sympy`, `slider_bearing_sympy`, `similarity_reduce_sympy`,
  `stokes_second_sympy`, `stokes_sphere_sympy`, `oseen_limit_sympy`) each < 5 s; Oseen contour grids 200² (FAST 100²).
- **Scalar-callable for parity:** every function an explainer mirrors must accept Python floats and return floats or
  small dicts of floats (`shot.py` evaluates `py:` rows); `circular_couette(R2=np.inf)` and `R1=0` exact branches;
  sphere fields return NaN for r < a.
- **Conventions to encode in names and docstrings:** `dpdx` (book sign) with a `G=` alias meaning −dp/dx (ch04 parity);
  θ from the downstream axis in `core.creeping` and a `frame="body" | "fluid"` argument; `Re` docstrings say diameter or
  radius (`proudman_pearson_drag_coefficient` converts); `stokes_layer` returns both `delta_book` = 4√(ν/ω) and
  `delta_e` = √(2ν/ω); `similarity_variable` returns y/√(νt) and a `half=True` option for y/(2√(νt)); `G0` default in
  `core.creeping`.
- **Typos to test as discriminating wrong variants:** R8 printed p − p_e (first power) and (1 − αx/L) integrands; R8b
  `form="book"` profile; R9 V without 1/h (units); R10 ∫ω dy = −U; R12 +3μaU cos θ/2r²; R6 ∂p/∂x in (8.13b) (coefficient
  set differs); R7 (8.17a) without ν (units); the biharmonic applied to ψ (must not vanish); ch05 `rotating_cylinder_flow`
  with ω = Ω₁ (factor 2); `diffusing_vortex_sheet` with γ = +2U.
- **Parity with earlier chapters:** `channel_flow` = ch04 `exact_solution("couette", G=−dpdx)` and `plane_poiseuille`;
  `pipe_poiseuille` = ch04 `exact_solution("pipe_poiseuille")` and ch03 `pipe_profile` (z → ∞); `circular_couette`
  limits = ch05 `rotating_cylinder_flow(omega=2Ω₁)` and `core.vortices.solid_body_rotation`; `stokes_first_problem`
  moved to `core.laminar` and re-exported from ch04; `vortex_sheet_diffusion` = ch05 `diffusing_vortex_sheet(γ=−2U)`;
  `line_vortex_decay` = `core.vortices.gaussian_vortex(σ = 2√(νt))`; `stokes_drag_coefficient` = ch04
  `sphere_drag_coefficient("stokes")`.
- **Promotion:** `core.laminar`, `core.lubrication`, `core.creeping` and `core.diffusion.crank_nicolson_1d` are for Ch. 9,
  11, 13, 16; `ch08.similarity_reduce_sympy` moves to `core` when Ch. 9 (Blasius) calls it. JS: E5 and E6 both need a
  "profiles at several times with a raw/rescaled mode" stage and E1/E3 a "channel with profile arrows and tracers" stage —
  the knowledge-keeper should watch for `Viz.profiles` / `Viz.channel` helpers after the build; `Viz.num.erf/erfc`
  already exist.

## 9. What the implementer must add beyond analysis §4
Functions the A items' figures and the explainers need that analysis §4 did not plan (all scalar-callable, documented
with book § and equations, validation label):

| Function (module) | For | What |
|---|---|---|
| `core.laminar.couette_poiseuille_state(h, U, dpdx, mu)` | E1, IF1 | dict(Q, V, Q_couette, Q_poiseuille, tau_bottom, tau_top, backflow, threshold, u_max, y_umax, y_reversal) |
| `core.lubrication.lubrication_term_magnitudes(L, h, U, rho, mu, p_scale="viscous")` | C05, E2, IF (C05 figure) | dict of the coefficient of every term of (8.16a) and (8.16b) (inertia, pressure, streamwise and cross-gap diffusion) for either pressure scale |
| `core.lubrication.slider_bearing_state(h0, alpha, L, U, mu, p_e=0.)` | C07, E3 | dict(C1, q, p_max, x_pmax, W_exact, W_linear, err_linear, p_max_atm, inlet_backflow) |
| `core.lubrication.slider_gap_velocity(x, y, h0, alpha, L, U, mu)` | E3, C06 figure | u(x, y) in the gap from (8.19) with the exact slider dp/dx (tracers and station profiles) |
| `core.lubrication.thin_film_state(t, area, rho, g, mu)` | E4 | dict(x_N, h_centre, beta, eta_N, effective_diffusivity) from the similarity solution |
| `core.laminar.stokes_first_state(t, U, nu, level=0.01)` | E5 | dict(sqrt_nut, delta, tau_w, vorticity_content) |
| `ch08.similarity_collapse_error(case, n, m, times)` | E6, C10 figure | spread (max abs difference) of rescaled profiles for trial exponents; zero at the correct n, m for each case |
| `core.laminar.stokes_layer_state(nu, omega, y)` | E7 | dict(delta_e, delta_book, amplitude, phase_lag, time_lag, crest_speed, wavelength) at height y |
| `core.creeping.oseen_velocity(r, theta, U, a, Re, frame="body")` | E8, A4 | (u_r, u_θ) from (8.53) by analytic differentiation (−expm1 form), for tracers |
| `core.creeping.settling_state(a, rho_p, rho, mu, g=G0)` | E9 | dict(U_t, Re, D, C_D, valid (Re < 0.1), D_pressure, D_friction) |
| `core.laminar.circular_couette_state(R1, R2, Omega1, Omega2, mu, rho)` | backup B1, IF3 | dict(A, B, torque_inner, torque_outer, power_in, dissipation, rayleigh_stable) |
