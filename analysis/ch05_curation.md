# Chapter 5 — Vorticity Dynamics: curation
(from `analysis/ch05.md` — 94 inventory rows, 34 equation labels (5.1)–(5.33), 31 derivations in §2b, 16 candidate A
ideas; `book.yaml` ch05, `policy.tier_a_full_treatment: [12, 18]`, `coverage: exhaustive`. Curated 2026-09-23,
concept-curator. Read: `knowledge/CUMULATIVE.md`, `concept_map.md`, `primers.md` (P01–P133), `notation.md`,
`viz_patterns.md`, `knowledge/ch04.md` §8, `analysis/ch04_curation.md` (format).)

Counts: A 14 · B 69 · C 11 · RECAP 18 · SKIP 2 · derivations written out 23 (★ 5 ★★ 14 ★★★ 4) · demoted to statements 6

Tier words (parsed by `tools/nbkit.py`): **CORE 14 · NOTE 60 · RECAP 18 · SKIP 2** = 94 rows; **DERIVATION 23**.
Depth: A 14 (all CORE) · B 69 (55 NOTE + R01, R02, R03, R05, R06, R07, R08, R09, R11, R14, R15, R16, R17, R18) · C 11
(N01, N04, N18, N23, N51 + R04, R10, R12, R13 + S01, S02). Every inventory row number `[#n]` appears exactly once in §2.

**Decisions that shape this chapter.**
1. **Fourteen A items from the analyst's sixteen (15 % of the rows).** Two merges: (a) the three sources of vorticity
   and the four restrictions (analyst A4) are B items inside **C03 Kelvin's theorem** — they are read straight off the
   surviving terms of (5.10)–(5.11), which C03 derives; the lock exchange (N15) moves to **C04** as its worked number;
   (b) the vorticity-diffusion test fields (analyst A8: diffusing sheet, Hill) are B/C inside **C06 (5.13)**, and
   Burgers' vortex (N21) is taught in **C10 stretching and tilting**, where its stretching–diffusion balance is the
   point (D18). The analyst's A12/A13 are kept apart on purpose: stretching vs tilting (C10) is the turbulence idea
   (Ch. 12), the absolute circulation and the fluid column (C11) are the climate idea (Ch. 13 potential vorticity) —
   Shammunul needs both at full depth. Inside C11 the planetary stretching of #67 is a B item and the column is C11's
   worked derivation (D20).
2. **SEEN rows are always RECAP** (18 rows, R01–R18, in inventory order), even when they carry a new equation number:
   (5.1) and (5.2) are ch03's (3.22) and (3.25) renamed (R02, R03, with the ⚠️ "ω is the vorticity, the tank turns at
   ω/2" callout), (5.5b) is hydrostatics, (5.19)–(5.23) are ch04's continuity, rotating NS, Lamb identity and viscous
   force in index form, ω + 2Ω is ch03 D13 turned round. Each is reminded where it is used, with the earlier chapter's
   function (B) or one sentence (C).
3. **Derivations written out: 23** (★ 5 · ★★ 14 · ★★★ 4). All A-block chains, plus nine results the book never
   writes out (rule (b)): the net viscous force of the line vortex (Ex. 5.4, D03), the pressure torque = baroclinic
   term (D06), the lock-exchange rate (Ex. 5.5, D07), the segment law and Γ/2πd (D13), the stretched tube and Burgers
   balance (Ex. 5.12, D18), Kelvin in a rotating frame (Ex. 5.10, D19), the fluid column (D20), the centre of vorticity
   (Ex. 5.18, D21), and the Green's-function solution of (5.14) whose printed sign is wrong (Ex. 5.9, folded into D10).
   Two analyst rows are merged into one written derivation each (a-D03 + a-D04 → D02; a-D13 + a-D14 → D10).
   **6 demoted to statements** (§4c): the vortex-line ratios, B across the tank, the line-vortex pressure (same moves
   as D02), the rotating-cylinder torque and dissipation, ∇·ω = 0, the planetary component equations.
4. **Book slips handled in our own words** (analysis §9): (5.14)'s −1/(4π) → +1/(4π) and the compensating slip in the
   integrand rewrite (N25, N26, D10, D11 — a test shows the printed sign reverses the swirl); "Exercise 5.8" → 5.9
   (N25) and "Exercise 5.11" → 5.10 (C11); the ½ in the Lamb step (R11, D09); Π for Φ in (5.26) (N33); the silently
   dropped u_{j,j}(ω_n + 2Ω_n) in (5.27) (N34, D15); Fig. 5.16's u₁ − u₂ vs the text's u₂ − u₁ (C14, D23); Fig. 5.2's
   "2ω" label (R02, D02); "single valued" is not why ∮dp/ρ = 0 (R09, D05); "irrotational C ⇒ no viscous term" holds
   for incompressible constant-μ flow only (N13, D05); "hyperboloids of the second degree" are cubic surfaces (N07).
5. **Conventions stated where first used** (⚠️ callouts): ω is the *vorticity* in (5.1) — ch03's `omega0` was the
   rotation rate (R02); Γ now has five meanings, and `gamma` [m/s] is the sheet strength while `Gamma` [m²/s] is always
   circulation (C14); σ = viscous stress here, core radius in ch03 (N06); ε = Levi-Civita vs dissipation (N08);
   counterclockwise-positive circulation and planar vorticity; relative vs absolute vorticity in (5.30) (C09, R18);
   the "Boussinesq" (5.30) keeps the full 1/ρ in the pressure term (D15 trap).
6. **Climate hooks, where real:** the rotating tank as the lab analogue of a geostrophic free surface (C02); tornado
   and dust-devil pressure deficits from (5.7) (N07); baroclinic generation — sea breeze, fronts, lock exchange (C04);
   planetary stretching, the fluid column and the ring of air moving poleward (C11 → Ch. 13 PV and Bjerknes); binary
   cyclones orbiting each other (Fujiwhara, C12); vortex stretching as the engine of the turbulent cascade (C10 →
   Ch. 12).

## 1. Teaching order (A IDs grouped by book section, B/C IDs under each; one sentence each: "once you see X, Y follows")
A items in **bold**; B and C items listed where they are taught, with depth in brackets. Equations are written next to
their numbers so downstream agents have them.

**§5.1 Introduction**
- **C01 Vortex lines, vortex tubes and their strength; tubes cannot end (5.4)** [#8]: once vortex lines are drawn
  tangent to ω, $dx/\omega_x=dy/\omega_y=dz/\omega_z$ *(Eq. 5.3)*, and bundled through a closed curve into a tube whose
  strength is $\Gamma=\oint_C\mathbf u\cdot d\mathbf x=\int_A\boldsymbol\omega\cdot\mathbf n\,dA$, Gauss on a piece of
  tube with ∇·ω = 0 gives $\int_V\nabla\cdot\boldsymbol\omega\,dV=-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0$
  *(Eq. 5.4)* (D01): the strength is the same all along the tube, so a tube (and a line) cannot end inside the fluid —
  where it thins, ω must grow. Inside: R01 [B] ω = ∇×u and spin ½ω, N01 [C] the chapter's programme, N02 [B] (5.3),
  N03 [B] tube strength, N45 [B] Fig. 5.1 (plotly 3-D stream tube vs vortex tube).
- **C02 Pressure and viscous stress in the two basic vortices (5.6)** [#13]: once the solid-body vortex
  $u_\theta=\omega r/2$ *(Eq. 5.1)* is put in Euler (no deformation ⇒ no viscous stress) and the radial and vertical
  balances $-\rho u_\theta^2/r=-\partial p/\partial r$ *(Eq. 5.5a)*, $0=-\partial p/\partial z-\rho g$ *(Eq. 5.5b)* are
  integrated (D02), $p-p_o=\tfrac18\rho\omega^2r^2-\rho gz$ *(Eq. 5.6)* — isobars and the free surface are paraboloids
  and the Bernoulli function grows outward; the same moves for $u_\theta=\Gamma/2\pi r$ *(Eq. 5.2)* give the funnel
  $p-p_\infty=-\rho\Gamma^2/8\pi^2r^2-\rho gz$ *(Eq. 5.7)*, a flow with nonzero viscous stress but zero net viscous
  force (D03): **irrotational means no net viscous force, not no stress**. Inside: R02 [B] (5.1), R03 [B] (5.2), R04
  [C] paddle wheels, N04 [C] reconnection, R05 [B] S = 0 ⇒ τ = −pδ, N05 [B] (5.5a), R06 [B] (5.5b), R07 [B] B across
  streamlines, N06 [B] σ_rθ and zero net force, N07 [B] (5.7), R08 [B] rotating cylinder = Rankine, N08 [B] torque to
  infinity and dissipation, N09 [B] the principle, N46 [B] Fig. 5.2, N47 [B] Fig. 5.3.

**§5.2 Kelvin's Circulation Theorem**
- **C03 Kelvin's circulation theorem (5.8)** [#20]: once the circulation of a *material* loop is differentiated with
  the loop parametrised by fixed particle labels, $\frac{D\Gamma}{Dt}=\oint_C\frac{Du_i}{Dt}dx_i+\oint_Cu_i\frac{D}{Dt}(dx_i)$
  *(Eq. 5.9)* and the second integral is $\oint d(\tfrac12u^2)=0$ (D04), substituting the momentum equation
  *(Eq. 5.10)* leaves $\frac{D\Gamma}{Dt}=\oint_C\big(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j}\big)dx_i$
  *(Eq. 5.11)* when the fluid is barotropic and the body force conservative (D05), and $D\Gamma/Dt=0$ *(Eq. 5.8)* when
  it is also inviscid — so circulation can only be made by nonconservative forces, baroclinicity or viscosity. Inside:
  N10 [B] (5.9), N11 [B] D(dx)/Dt = du, N12 [B] (5.10), R09 [B] ∮dp/ρ = 0 needs ρ(p), N13 [B] (5.11) with Lamb–Oseen,
  N14 [B] three sources, N16 [B] four restrictions (decision table), N48 [B] Fig. 5.4 (animated loop).
- **C04 Barotropic vs baroclinic fluid element: the pressure torque (Fig. 5.6)** [#29]: once a small disc is drawn
  with its isobars and its isopycnals, the pressure force passes through the disc's centre but the centre of mass sits
  toward the heavy side, so when ∇ρ is not parallel to ∇p the pressure force exerts a torque; dividing by the moment of
  inertia and doubling gives exactly $\frac1{\rho^2}\nabla\rho\times\nabla p$ (D06) — the baroclinic term the book
  derives later as (5.28); the lock exchange spins up at $2(\rho_2-\rho_1)g/((\rho_2+\rho_1)\delta)$ (D07). Inside:
  N15 [B] lock exchange, N49 [B] Fig. 5.5, N50 [B] Fig. 5.6; forward link to N35 (5.28).

**§5.3 Helmholtz's Vortex Theorems**
- **C05 Helmholtz's vortex theorems** [#30]: once Kelvin's theorem is applied to a surface lying on a tube wall (zero
  circulation on its edge stays zero), the tube wall stays a tube wall — vortex lines move with the fluid; with (5.4)
  (strength constant along a tube, tubes cannot end) and (5.8) (strength constant in time) all four theorems follow
  (D08). The field-equation proof (ω and a material element obey the same linear equation when ν = 0) is shown right
  after C06, where (5.13) is available. Inside: N17 [B] the proof, N51 [C] Fig. 5.7.

**§5.4 Vorticity Equation in a Nonrotating Frame**
- **C06 The vorticity equation (5.13)** [#37]: once the curl of incompressible Navier–Stokes kills the pressure and
  gravity terms, $\nabla\times\{\frac{D\mathbf u}{Dt}=-\frac1\rho\nabla p+\mathbf g+\nu\nabla^2\mathbf u\}$
  *(Eq. 5.12)*, and the Lamb identity plus $\nabla\times(\boldsymbol\omega\times\mathbf u)=(\mathbf u\cdot\nabla)\boldsymbol\omega-(\boldsymbol\omega\cdot\nabla)\mathbf u$
  reorganise the rest (D09 ★★★), $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$
  *(Eq. 5.13)*: vorticity is carried, stretched and tilted by the flow and diffused by viscosity — and pressure cannot
  touch it. Inside: N18 [C] §5.4 hypotheses, R10 [C] ∇·ω = 0, N19 [B] (5.12), R11 [B] Lamb curl, N20 [B] (B.3.10),
  N22 [B] diffusing vortex sheet (erf, `slider_figure` in t), N23 [C] Hill's vortex; the frozen-in check of C05.

**§5.5 Velocity Induced by a Vortex Filament: Law of Biot and Savart**
- **C07 The Biot–Savart law (5.16)** [#46]: once $\nabla\times\boldsymbol\omega=-\nabla^2\mathbf u$ turns the velocity
  into the solution of a Poisson equation, the Green's function of ∇² gives
  $\mathbf u=+\frac1{4\pi}\int_{V'}\frac{\nabla'\times\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}d^3x'$ — the book
  prints −1/(4π) *(Eq. 5.14)* (D10 ★★★) — and a product rule for the curl, Gauss in curl form *(Eq. 5.15)* and a
  tube-aligned volume leave (D11 ★★★)
  $\mathbf u(\mathbf x,t)=\frac1{4\pi}\int_{V'}\frac{\boldsymbol\omega(\mathbf x',t)\times(\mathbf x-\mathbf x')}{\lvert\mathbf x-\mathbf x'\rvert^3}d^3x'$
  *(Eq. 5.16)*: vorticity here sets the velocity everywhere, like a current sets a magnetic field. Inside: N24 [B]
  Poisson, N25 [B] (5.14) and its sign, N26 [B] integrand rewrite, N27 [B] (5.15), N28 [B] choice of V′, N52 [B]
  Fig. 5.8.
- **C08 The filament law (5.17)** [#47]: once a thin tube is seen from far away, the kernel can be pulled out of the
  cross-section integral and $\int\lvert\boldsymbol\omega\rvert d^2x'=\Gamma$ (D12), each length dl induces
  $d\mathbf u=\frac{\Gamma\,dl}{4\pi}\mathbf e_\omega\times\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$
  *(Eq. 5.17)*; a straight segment gives $(\Gamma/4\pi d)(\cos\theta_a-\cos\theta_b)$ and an infinite line gives back
  $\Gamma/2\pi d$, i.e. (5.2) (D13) — the tool for rings (C13), horseshoe vortices and induced drag (Ch. 14). Inside:
  N29 [B] the infinite line recovered.

**§5.6 Vorticity Equation in a Rotating Frame**
- **C09 The vorticity equation in a rotating frame with baroclinic generation (5.30)** [#62]: once the rotating-frame
  momentum equation *(Eq. 5.20)* is rewritten in Lamb form
  $\partial u_i/\partial t+(\tfrac12u_j^2+\Phi)_{,i}-\varepsilon_{ijk}u_j(\omega_k+2\Omega_k)=-(1/\rho)p_{,i}-\nu\varepsilon_{ijk}\omega_{k,j}$
  *(Eq. 5.25)* (D14) and its curl is taken in index notation (D15 ★★★),
  $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u+\frac1{\rho^2}\nabla\rho\times\nabla p+\nu\nabla^2\boldsymbol\omega$
  *(Eq. 5.30)*: the absolute vorticity is stretched and tilted, density gradients across pressure gradients create
  vorticity, viscosity diffuses it; Ω = 0 and ∇ρ ∥ ∇p give back (5.13). Inside: R12 [C] §5.6 hypotheses, N30 [B]
  (5.18), R13 [C] (5.19), R14 [B] (5.20), R15 [B] (5.21), R16 [B] (5.22), R17 [B] (5.23), N31 [B] (5.24), N32 [B]
  (5.25), N33 [B] (5.26), N34 [B] (5.27), N35 [B] (5.28) baroclinic vector, N36 [B] (5.29), N37 [B] reading the terms.
- **C10 Stretching and tilting of vortex lines (5.32)** [#66]: once natural coordinates along a vortex line turn the
  stretching–tilting term into $(\boldsymbol\omega\cdot\nabla)\mathbf u=\omega\,\partial\mathbf u/\partial s$ *(Eq. 5.31)*
  (D16), its components $\frac{D\omega_s}{Dt}=\omega\frac{\partial u_s}{\partial s}$,
  $\frac{D\omega_n}{Dt}=\omega\frac{\partial u_n}{\partial s}$, $\frac{D\omega_m}{Dt}=\omega\frac{\partial u_m}{\partial s}$
  *(Eq. 5.32)* say: stretching along the line spins it up (a skater pulling in her arms), shear across it tilts it,
  and neither happens in 2-D (D17); an inviscid tube stretched to twice its length has twice the vorticity, and in
  Burgers' vortex stretching is held off by diffusion in a core of radius $\sqrt{4\nu/\alpha}$ (D18). Inside: N38 [B]
  (5.31), N21 [B] Burgers' vortex, N53 [B] Fig. 5.9 (plotly 3-D helix with its frame).
- **C11 Kelvin in a rotating frame: absolute circulation (5.33) and the fluid column** [#68]: once the Kelvin
  derivation is repeated with the Coriolis term, whose loop integral is $-2\boldsymbol\Omega\cdot dA_{\text{vec}}/dt$
  (D19), $\frac{D\Gamma_a}{Dt}=0$ with
  $\Gamma_a\equiv\int_A(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\mathbf n\,dA=\Gamma+2\int_A\boldsymbol\Omega\cdot\mathbf n\,dA$
  *(Eq. 5.33)*; applied to a thin loop round a column whose mass Ah is fixed it gives $(\omega_z+2\Omega)/h$ = const
  (D20) — stretch a column and it spins up cyclonically, squash it and it spins anticyclonically: the seed of potential
  vorticity (Ch. 13). Inside: R18 [B] absolute and planetary vorticity, N39 [B] planetary stretching/tilting, N54 [B]
  Fig. 5.10.

**§5.7 Interaction of Vortices**
- **C12 Point vortices move each other: pairs and the centre of vorticity** [#69]: once a real vortex is idealised as
  a line vortex that moves with the local flow (Helmholtz 1) and each vortex is advected by the Biot–Savart velocity of
  all the others,
  $\frac{d\mathbf x_k}{dt}=\sum_{j\ne k}\frac{\Gamma_j}{2\pi}\frac{\mathbf e_z\times(\mathbf x_k-\mathbf x_j)}{\lvert\mathbf x_k-\mathbf x_j\rvert^2}$,
  two same-sign vortices orbit their centre of vorticity at $(\Gamma_1+\Gamma_2)/2\pi h^2$ and an opposite pair
  translates at $\Gamma/2\pi h$ (D21). Inside: N40 [B] same-sign pair (Fig. 5.11), N41 [B] opposite pair and the knife
  blade, N55 [B] Fig. 5.11, N56 [B] Fig. 5.12.
- **C13 The method of images: a vortex near a wall** [#72]: once a wall is replaced by a mirror vortex of opposite
  sign, the wall becomes a streamline and the vortex drifts along it at $\Gamma/4\pi h$ (D22); circle and channel
  images follow the same idea, and a vortex ring approaching a wall widens and slows as its image ring pulls it (the
  filament law C08 applied to rings), while two coaxial rings leap-frog. Inside: N42 [B] ring near a wall, N43 [B]
  leap-frogging (animation), N57 [B] Fig. 5.13 (knife-blade pair in a bucket), N58 [B] Fig. 5.14, N59 [B] Fig. 5.15.

**§5.8 Vortex Sheet**
- **C14 The vortex sheet: strength = jump in tangential velocity** [#76]: once a sheet is built as a row of line
  vortices (the velocity just above and below differs, the normal velocity does not), the circulation round a thin
  dn × ds circuit gives $d\Gamma=u_2\,ds+v\,dn-u_1\,ds-v\,dn=(u_2-u_1)\,ds$ (D23): the sheet's strength per unit length
  is the jump in tangential velocity — the model of every shear layer, wing wake and Kelvin–Helmholtz roll-up (Ch. 11,
  14). Inside: N44 [B] sheet definition, N60 [B] Fig. 5.16 (discrete row, u(y) for N = 10, 100, 1000).

**Exercises and literature** — S01, S02 (one pointer line each at the end of the notebook).

## 2. Chapter map (depth) — every inventory row exactly once
One row per inventory row (94, `[#n]` = analysis §2 row). `A parent` names the A block a B or C item is written in.
IDs: CORE C01–C14, NOTE N01–N60 and RECAP R01–R18 in inventory order (figure rows #79–#94 last), SKIP S01–S02.

| ID | Item | § | Depth | Tier | A parent | Reason (A) / treatment (B) / pointer (C, SKIP) / source chapter (RECAP) |
|---|---|---|---|---|---|---|
| R01 | Def: vorticity = twice the particle angular velocity, $\boldsymbol\omega=\nabla\times\mathbf u$, spin $\tfrac12\boldsymbol\omega$; vortex and vortex motion [#1] | 5.1 | B | RECAP | C01 | ch03 C10 (D12) and N42: reminded with `core.kinematics.vorticity` on the solid-body and line-vortex fields at the top of C01 |
| R02 | Eq. (5.1): solid-body rotation from uniform vorticity, $u_\theta=\omega r/2$ (ω is the vorticity; the rotation rate is ω/2) [#2] | 5.1 | B | RECAP | C02 | ch03 C13 (3.22): reminded with `ch05.solid_body_from_vorticity` → `core.vortices.solid_body_rotation(r, omega/2)`; ⚠️ convention callout: a tank turning at 1 rad/s has ω = 2 s⁻¹ (Fig. 5.2's "2ω" label is a slip) |
| R03 | Eq. (5.2): ideal line vortex $u_\theta=\Gamma/2\pi r$, irrotational for r > 0 [#3] | 5.1 | B | RECAP | C02 | ch03 C13 (3.25)–(3.26), B = Γ/2π: reminded with `ch05.line_vortex_gamma` → `core.vortices.line_vortex`; ∮u·dx = Γ on any circle printed |
| R04 | Both basic vortices are steady with circular streamlines; only the first rotates particles (paddle-wheel test) [#4] | 5.1 | C | RECAP | C02 | ch03 E6 `vortex_paddle_wheels`: one sentence with the link (not duplicated) |
| N01 | Programme: vorticity is embedded in elements and is reoriented, concentrated or diffused by their motion, deformation and neighbours' torques [#5] | 5.1 | C | NOTE | C01 | named in the chapter's opening paragraph as the map of C06 (carried), C10 (stretched, tilted), C04/C09 (baroclinic torque) and N22 (diffused) |
| N02 | Def: vortex line, tangent to ω everywhere; Eq. (5.3): $dx/\omega_x=dy/\omega_y=dz/\omega_z$ [#6] | 5.1 | B | NOTE | C01 | stated as the streamline (3.7) with ω for u (a-D01 given, §4c); the arc-length form dx/ds = ω/∣ω∣ (P92, P93 reminders) traced with `core.vorticity.vortex_line` on the helical test field u_φ = aRz (lines satisfy zR² = const) |
| N03 | Def: vortex tube and its strength $\Gamma=\oint_C\mathbf u\cdot d\mathbf x=\int_A\boldsymbol\omega\cdot\mathbf n\,dA$ (dΓ = ω·n dA ↔ dQ = u·n dA) [#7] | 5.1 | B | NOTE | C01 | stated with the stream-tube analogy table and both routes computed on the Lamb–Oseen vortex (`core.vorticity.vortex_tube_strength`: circulation route = flux route, orientation flip changes the sign) |
| C01 | Eq. (5.4): Gauss on a piece of vortex tube, $\int_V\nabla\cdot\boldsymbol\omega\,dV=-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0$ — tube strength constant along the tube; tubes cannot end in the fluid [#8] | 5.1 | A | CORE | – | load-bearing: Helmholtz 2–3 (C05), the filament strength in (5.17) (C08), ω ∝ 1/area in stretching (C10), trailing vortices carrying the bound circulation of a wing (Ch. 14); a kinematic fact valid in any flow |
| N04 | Viscosity diffuses vorticity and reconnects vortex lines; the two basic vortices re-examined with viscosity [#9] | 5.1 | C | NOTE | C02 | named in one sentence opening C02's viscous part; reconnection is not developed in the book (pointer: Ch. 12 turbulence) |
| R05 | Solid-body rotation has S_ij = 0, so τ_ij = −pδ_ij and Cauchy reduces to Euler [#10] | 5.1 | B | RECAP | C02 | ch03 D15, ch04 C07/C08: reminded with `core.constitutive.newtonian_stress` on the solid-body gradient (→ −pδ) and `core.navier_stokes.viscous_force_forms` (→ 0); step 2 of D02 |
| N05 | Eq. (5.5a): $-\rho u_\theta^2/r=-\partial p/\partial r$ — the pressure gradient supplies the centripetal acceleration [#11] | 5.1 | B | NOTE | C02 | stated as steps 3–5 of D02 (the −u_θ²/r from the turning e_θ, filled in); `ch05.solid_body_pressure_gradients` with a sympy check through `core.curvilinear` |
| R06 | Eq. (5.5b): $0=-\partial p/\partial z-\rho g$ (hydrostatic in the vertical) [#12] | 5.1 | B | RECAP | C02 | ch01 C20 (1.8): reminded as step 6 of D02 with `ch01` hydrostatic function |
| C02 | Eq. (5.6): pressure in a steadily rotating tank, $p-p_o=\tfrac18\rho\omega^2r^2-\rho gz$; isobars $z=\omega^2r^2/8g-(p-p_o)/\rho g$ are paraboloids [#13] | 5.1 | A | CORE | – | load-bearing: the cyclostrophic/centripetal pressure balance of every vortex (tornadoes, dust devils, the gradient wind of Ch. 13), the rotating-tank experiment behind Ch. 13's lab analogues, and — with (5.7) and D03 — the lesson that irrotational flow has stress but no net viscous force (Ch. 6) |
| R07 | $-\tfrac12u_\theta^2+gz+p/\rho$ = const in the tank, so $B=u_\theta^2/2+gz+p/\rho$ varies across streamlines (rotational flow) [#14] | 5.1 | B | RECAP | C02 | ch04 C11 (4.69)–(4.72): reminded with `ch05.bernoulli_across_vortex` (B − B(0) = ω²r²/4 in the tank, 0 for the line vortex; a-D05 given, §4c); link to ch04 `which_bernoulli` |
| N06 | Line vortex: $\sigma_{r\theta}=\mu[\frac1r\frac{\partial u_r}{\partial\theta}+r\frac{\partial}{\partial r}(\frac{u_\theta}{r})]=-\mu\Gamma/\pi r^2\ne0$ but zero net viscous force (Exercise 5.4) [#15] | 5.1 | B | NOTE | C02 | stated with D03 written out (the book leaves it to Exercise 5.4); `ch05.line_vortex_viscous_stress`, three routes to zero net force (`core.navier_stokes.viscous_force_forms`); ⚠️ σ is the viscous stress here, the core radius in ch03 |
| N07 | Eq. (5.7): line-vortex pressure $p-p_\infty=-\rho\Gamma^2/8\pi^2r^2-\rho gz$; funnel isobars; $\tfrac12u_\theta^2+gz+p/\rho$ constant everywhere [#16] | 5.1 | B | NOTE | C02 | stated (a-D07 given, §4c: the same moves as D02 with u_θ = Γ/2πr); number: Γ = 1 m²/s, water, r = 0.1 m → deficit 1266.5 Pa, funnel depth 0.129 m; Rankine composite `ch05.rankine_pressure` gives the funnel a finite bottom (tornado); ⚠️ the book's "hyperboloids of the second degree" are cubic surfaces (c − z)r² = const |
| R08 | Rotating solid cylinder: $u_\theta=\omega a^2/2r$ for r ≥ a — with the fluid inside it is a Rankine vortex of Γ = πa²ω [#17] | 5.1 | B | RECAP | C02 | ch03 C14 (3.28): reminded with `ch05.rotating_cylinder_flow` → `core.vortices.rankine_vortex`; the viscous solution itself is derived in Ch. 8 (8.11) |
| N08 | Dissipation outside the rotating cylinder equals the wall's work; torque $2\pi r^2\sigma_{r\theta}=-2\mu\Gamma$ the same at every radius, carried to infinity [#18] | 5.1 | B | NOTE | C02 | stated with our formulas (a-D08 given, §4c) and a numeric energy balance (`ch05.torque_per_length`, `dissipation_outside_cylinder`: ∫ρε dA = power in − power out); ⚠️ ε here is the dissipation, not Levi-Civita |
| N09 | Principle: irrotational ≠ no viscous stress; irrotational = no net viscous force, $\mu\nabla^2\mathbf u=-\mu\nabla\times\boldsymbol\omega=0$ (incompressible); solid-body rotation alone has no stress at all [#19] | 5.1 | B | NOTE | C02 | stated as C02's conclusion with a computed three-row table: solid body (σ = 0, F = 0), line vortex (σ ≠ 0, F = 0), Lamb–Oseen (σ ≠ 0, F ≠ 0) |
| C03 | Kelvin's circulation theorem, Eq. (5.8): $D\Gamma/Dt=0$ for a material loop (inviscid, barotropic, conservative body forces, inertial frame) [#20] | 5.2 | A | CORE | – | load-bearing: Helmholtz's theorems (C05), irrotational flow stays irrotational (Ch. 6 potential flow, the starting vortex and lift in Ch. 6 and 14), the rotating-frame form (5.33) (C11) and Bjerknes' theorem and PV conservation (Ch. 13) |
| N10 | Eq. (5.9): $\frac{D\Gamma}{Dt}=\oint_C\frac{Du_i}{Dt}dx_i+\oint_Cu_i\frac{D}{Dt}(dx_i)$ [#21] | 5.2 | B | NOTE | C03 | stated as D04's result; `core.vorticity.kelvin_rate_terms` (acceleration term + contour term = dΓ/dt by central differences) |
| N11 | Contour-element kinematics $D(d\mathbf x)/Dt=d\mathbf u$ (Fig. 5.4), hence $\oint u_i\,du_i=\oint d(\tfrac12u_i^2)=0$ [#22] | 5.2 | B | NOTE | C03 | stated as steps 4–7 of D04 (ch03 P99 material element reminder); the contour term ≈ 0 printed |
| N12 | Eq. (5.10): momentum substituted, $\oint_C\frac{Du_i}{Dt}dx_i=-\oint_C\frac1\rho dp-\oint_Cd\Phi+\oint_C(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j})dx_i$ [#23] | 5.2 | B | NOTE | C03 | stated as steps 1–4 of D05; `core.vorticity.kelvin_force_terms` (pressure, body, viscous line integrals) |
| R09 | Barotropic ⇒ $\oint dp/\rho=\oint d\mathcal P=0$; conservative ⇒ $\oint d\Phi=0$ [#24] | 5.2 | B | RECAP | C03 | ch04 C11 (N87, pressure function (4.67)): reminded in D05 step 5 with the baroclinic counterexample (∮dp/ρ ≠ 0 although p and ρ are single-valued — the book's reason is incomplete) |
| N13 | Eq. (5.11): $\frac{D\Gamma}{Dt}=\oint_C(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j})dx_i$ — Kelvin holds if inviscid or if the net viscous force vanishes on C [#25] | 5.2 | B | NOTE | C03 | stated as D05's result with the Lamb–Oseen number: Γ₀ = 0.01 m²/s, ν = 10⁻⁶ m²/s, circle r = 5 mm at t = 10 s → Γ = Γ₀(1 − e^{−r²/4νt}) = 0.00465 m²/s, and the (5.11) line integral equals ∂Γ/∂t (`ch05.lamb_oseen_circulation`); ⚠️ "C in irrotational fluid ⇒ no viscous term" is incompressible, constant-μ only |
| N14 | Three sources of vorticity: nonconservative body forces (Coriolis — drain vortex), nonbarotropic pressure–density relation (lock exchange), net viscous forces (walls) [#26] | 5.2 | B | NOTE | C03 | stated as a three-row table, one computed number per source (rotating-frame loop term, lock-exchange ∮dp/ρ, Lamb–Oseen viscous term) |
| N15 | Lock exchange (Fig. 5.5): the interface tilts, vorticity is created; initial rate $D\omega_z/Dt=2(\rho_2-\rho_1)g/((\rho_2+\rho_1)\delta)$ (Exercise 5.5) [#27] | 5.2 | B | NOTE | C04 | stated with D07 written out (C04's worked number): fresh ρ₁ = 1000, salt ρ₂ = 1025 kg/m³, δ = 0.1 m → 2.42 s⁻²; field route `core.vorticity.baroclinic_term` on a tanh step agrees with `ch05.lock_exchange_initial_vorticity_rate`; sea breeze and fronts → Ch. 13 |
| N16 | Four restrictions keeping irrotational flow irrotational: no net viscous force on C, conservative body forces, barotropic (else baroclinic), inertial frame [#28] | 5.2 | B | NOTE | C03 | stated as a decision table `ch05.kelvin_hypotheses` (16 combinations → the surviving terms of (5.10)), mirrored in E3's status |
| C04 | Barotropic vs baroclinic element (Fig. 5.6): $\nabla\rho\parallel\nabla p\Leftrightarrow\nabla\rho\times\nabla p=0$; otherwise the pressure force misses the centre of mass and its torque changes the vorticity [#29] | 5.2 | A | CORE | – | load-bearing: the physical meaning of the baroclinic term (5.28)/(5.30) (C09) and the mechanism of sea breezes, fronts, thermal wind and baroclinic instability (Ch. 13), internal-wave generation (Ch. 7) |
| C05 | Helmholtz's vortex theorems: (1) vortex lines move with the fluid, (2) tube strength constant along the tube, (3) tubes cannot end in the fluid, (4) tube strength constant in time [#30] | 5.3 | A | CORE | – | load-bearing: "vortices move with the flow" is the rule of every point-vortex, ring and sheet model of §5.7–5.8 (C12–C14) and Ch. 10's vortex methods; frozen-in lines lead to frozen-in PV in Ch. 13 and wake vortices in Ch. 14 |
| N17 | Proof of (1) (Fig. 5.7): a surface on the tube wall has zero circulation on its edge; by Kelvin so does its material image, for every such surface, so the image lies on the tube wall [#31] | 5.3 | B | NOTE | C05 | stated as D08 (the "for every S" localisation filled in); demo `core.vorticity.material_loop` on a small loop on a tube wall of an inviscid axisymmetric strain + swirl flow — flux through it stays ≈ 0; the field-equation route (`frozen_in_check`) shown after C06 |
| N18 | Hypotheses of §5.4: ρ constant (barotropic), ν constant, conservative body force, inertial frame [#32] | 5.4 | C | NOTE | C06 | named in D09's "we start from" line; relaxed in §5.6 (C09) |
| R10 | $\nabla\cdot\boldsymbol\omega=\nabla\cdot(\nabla\times\mathbf u)=0$ used in §5.4 [#33] | 5.4 | C | RECAP | C06 | ch02 C11 (D26 corollary): one sentence in D09's tools line; proved again as (5.18) in N30 |
| N19 | Eq. (5.12): the curl of incompressible NS, $\nabla\times\{\frac{D\mathbf u}{Dt}=-\frac1\rho\nabla p+\mathbf g+\nu\nabla^2\mathbf u\}$; pressure and gravity curls vanish [#34] | 5.4 | B | NOTE | C06 | stated as steps 1–3 of D09; `core.vorticity.vorticity_equation_sym` returns each term's curl (pressure, gravity → 0) |
| R11 | Curl of the acceleration via the Lamb identity: $\nabla\times\{\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u\}=\frac{\partial\boldsymbol\omega}{\partial t}+\nabla\times(\boldsymbol\omega\times\mathbf u)$ [#35] | 5.4 | B | RECAP | C06 | ch04 C11 (4.68), D24: reminded as steps 4–6 of D09 with `core.navier_stokes.lamb_identity_terms`; ⚠️ the book writes ∇(u·u) for ∇(½u·u) (harmless: the curl kills it) |
| N20 | Identity (B.3.10) with ∇·u = ∇·ω = 0: $\nabla\times(\boldsymbol\omega\times\mathbf u)=(\mathbf u\cdot\nabla)\boldsymbol\omega-(\boldsymbol\omega\cdot\nabla)\mathbf u$ [#36] | 5.4 | B | NOTE | C06 | stated with the general form (+ω∇·u − u∇·ω) and a sympy proof for generic fields (ε–δ, ch02 D09 reminder); steps 7–8 of D09 |
| C06 | Eq. (5.13): the vorticity equation $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$ (constant ρ, conservative body force, inertial frame) [#37] | 5.4 | A | CORE | – | load-bearing: the dynamics of vorticity for the rest of the book — its rotating, baroclinic generalisation (5.30) (C09), stretching and tilting (C10), vortex decay (Ch. 8), vorticity–stream-function CFD (Ch. 10), enstrophy and the cascade (Ch. 12), QG vorticity (Ch. 13) |
| N21 | Burgers' vortex (Exercise 5.12): $u_R=-\tfrac12\alpha R$, $u_z=\alpha z$, $u_\varphi=\frac{\Gamma}{2\pi R}[1-e^{-\alpha R^2/4\nu}]$, $\omega_z=\frac{\alpha\Gamma}{4\pi\nu}e^{-\alpha R^2/4\nu}$ — steady because stretching balances diffusion [#38] | 5.4 | B | NOTE | C10 | stated with D18 written out (the balance; book never writes it out); number: α = 1 s⁻¹, ν = 10⁻⁶ m²/s, Γ = 10⁻³ m²/s → core radius √(4ν/α) = 2 mm, peak ω_z = 79.6 s⁻¹; `core.vortices.burgers_vortex`; term bars (advective + stretching = diffusion) via `core.vorticity.vorticity_terms` |
| N22 | Viscous diffusion of a vortex sheet (Exercise 5.6): $\omega_z(y,t)=\frac{\gamma}{2\sqrt{\pi\nu t}}\exp\{-\frac{y^2}{4\nu t}\}$ — (5.13) reduces to 1-D diffusion [#39] | 5.4 | B | NOTE | C06 | stated (the result of the exercise; the 1-D reduction is one line) with `slider_figure` over t, u = −(γ/2)erf(y/2√(νt)) (erf gloss, P123 reminder) and ∫ω dy = γ conserved; FTCS cross-check (`core.diffusion`); → Ch. 8 Stokes' first problem, Ch. 11 shear layers |
| N23 | Hill's spherical vortex (Exercise 5.11): $\psi=\frac{Aa^4}{10}\frac{R^2}{a^2}(1-\frac{R^2}{a^2}-\frac{z^2}{a^2})$, ω = AR e_φ — steady inviscid solution with stretched ring vortex lines [#40] | 5.4 | C | NOTE | C06 | named as a second exact solution of (5.13) with ν = 0 (`core.vortices.hill_spherical_vortex`, residual printed); the flow outside the sphere → Ch. 6 |
| N24 | $\nabla\times\boldsymbol\omega=\nabla(\nabla\cdot\mathbf u)-\nabla^2\mathbf u=-\nabla^2\mathbf u$ — a Poisson equation for u with source −∇×ω [#41] | 5.5 | B | NOTE | C07 | stated as steps 1–3 of D10 (curl-of-curl P122 reminder); periodic-box solver `core.biot_savart.velocity_from_vorticity_fft` recovers Taylor–Green to round-off (FFT primer); → Ch. 10 vorticity–stream-function CFD, Ch. 13 PV inversion |
| N25 | Eq. (5.14): the vorticity-induced velocity, Green's-function solution of the Poisson equation (Exercise 5.9; book cites "5.8"); **book prints −1/(4π), correct +1/(4π)** [#42] | 5.5 | B | NOTE | C07 | stated as D10's result with the ⚠️ "book prints X; we use Y" callout and the discriminating test (printed sign reverses the swirl of a Gaussian tube; `core.biot_savart.velocity_from_curl_omega(sign=…)`) |
| N26 | Rewrite of the integrand: $\frac{\nabla'\times\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}=\nabla'\times(\frac{\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert})+\boldsymbol\omega\times\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ (the book's second sign slip cancels the first) [#43] | 5.5 | B | NOTE | C07 | stated as steps 2–5 of D11 with a sympy check of ∇′(1/r) = (x − x′)/r³ |
| N27 | Eq. (5.15): Gauss in curl form, $\int_{V'}\nabla'\times(\frac{\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert})d^3x'=\int_{A'}\frac{\mathbf n\times\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}d^2x'$ [#44] | 5.5 | B | NOTE | C07 | stated as steps 6–8 of D11 (component-wise Gauss with ε_kij, ch02 C14 reminder); `core.integral_theorems.curl_theorem_box` volume = surface on a random polynomial field |
| N28 | Choice of V′ (Fig. 5.8): ends normal to ω (n × ω = 0), curved side outside the vortex (ω = 0) ⇒ the surface term vanishes [#45] | 5.5 | B | NOTE | C07 | stated as step 9 of D11; `curl_theorem_box` on a tube-aligned box shows the surface term ≈ 0 |
| C07 | Eq. (5.16): Biot–Savart law $\mathbf u(\mathbf x,t)=\frac1{4\pi}\int_{V'}\frac{\boldsymbol\omega(\mathbf x',t)\times(\mathbf x-\mathbf x')}{\lvert\mathbf x-\mathbf x'\rvert^3}d^3x'$ [#46] | 5.5 | A | CORE | – | load-bearing: velocity from vorticity — point vortices, images, rings and sheets (C12–C14), induced velocity and superposition in Ch. 6, vortex methods (Ch. 10), lifting line and downwash (Ch. 14), PV inversion (Ch. 13) |
| C08 | Eq. (5.17): filament law $d\mathbf u\cong\frac{\Gamma\,dl}{4\pi}\mathbf e_\omega\times\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ for a thin tube of strength Γ seen from afar [#47] | 5.5 | A | CORE | – | load-bearing: the working form of Biot–Savart — segments, polygons and rings (C13), horseshoe vortices and induced drag (Ch. 14); its infinite-line limit closes the loop with (5.2) |
| N29 | The infinite straight filament from (5.17) gives back $u_\theta=\Gamma/2\pi d$ (5.2); finite segment $(\Gamma/4\pi d)(\cos\theta_a-\cos\theta_b)$ [#48] | 5.5 | B | NOTE | C08 | stated with D13 written out (the book never closes this loop); number: Γ = 1 m²/s, d = 1 m → infinite line 0.159 m/s, semi-infinite 0.080 m/s, square loop of side 1 m at its centre 2√2Γ/πL = 0.900 m/s (`core.biot_savart.segment_induced_velocity`, `filament_velocity`) |
| R12 | Hypotheses of §5.6: steadily rotating frame, nonbarotropic, nearly incompressible (Boussinesq) ∇·u ≈ 0; comma notation $u_{i,j}\equiv\partial u_i/\partial x_j$ [#49] | 5.6 | C | RECAP | C09 | ch04 C13 (Boussinesq), ch02 (§2.14 comma notation): one sentence in D14's start line |
| N30 | Eq. (5.18): $\omega_{i,i}=\varepsilon_{inq}u_{q,ni}=0$ — ω is solenoidal for any flow [#50] | 5.6 | B | NOTE | C09 | stated (a-D18 given, §4c: symmetric × antisymmetric contraction, ch02 N61); `core.vorticity.vorticity_divergence` at round-off on a random smooth field |
| R13 | Eq. (5.19): continuity $u_{i,i}=0$ [#51] | 5.6 | C | RECAP | C09 | ch04 C02 (4.10): one sentence (`core.navier_stokes.continuity_residual`) |
| R14 | Eq. (5.20): momentum in a steadily rotating frame, $\frac{\partial u_i}{\partial t}+u_ju_{i,j}+2\varepsilon_{ijk}\Omega_ju_k=-\frac1\rho p_{,i}+g_i+\nu u_{i,jj}$ [#52] | 5.6 | B | RECAP | C09 | ch04 C09 (4.45) with Ω constant and effective gravity: reminded with `ch05.rotating_ns_residual` (solid body at Ω is rest in the co-rotating frame); link to ch04 `rotating_frame_coriolis`; D14's start |
| R15 | Eq. (5.21): $u_ju_{i,j}=-(\mathbf u\times\boldsymbol\omega)_i+\tfrac12(u_j^2)_{,i}$ (index Lamb identity) [#53] | 5.6 | B | RECAP | C09 | ch04 C11 (4.68): step 2 of D14 with `core.navier_stokes.lamb_identity_terms` |
| R16 | Eq. (5.22): $\varepsilon_{ijk}\omega_k=u_{j,i}-u_{i,j}$ (ε–δ) [#54] | 5.6 | B | RECAP | C09 | ch03 (3.15)–(3.16) R_ij = −ε_ijkω_k: step 3 of D14 with `core.kinematics.vorticity_from_gradient` on a random G |
| R17 | Eq. (5.23): $\nu u_{i,jj}=-\nu\varepsilon_{ijk}\omega_{k,j}$ (viscous term as −ν∇×ω) [#55] | 5.6 | B | RECAP | C09 | ch04 D13 (4.40): step 4 of D14 with `core.navier_stokes.viscous_force_forms` |
| N31 | Eq. (5.24): $2\varepsilon_{ijk}\Omega_ju_k=-2\varepsilon_{ijk}\Omega_ku_j$ (dummy indices swapped) [#56] | 5.6 | B | NOTE | C09 | stated as step 5 of D14 (relabelling, ch02 N43 reminder) with a one-line numeric check |
| N32 | Eq. (5.25): rotating NS in Lamb form, $\partial u_i/\partial t+(\tfrac12u_j^2+\Phi)_{,i}-\varepsilon_{ijk}u_j(\omega_k+2\Omega_k)=-(1/\rho)p_{,i}-\nu\varepsilon_{ijk}\omega_{k,j}$ [#57] | 5.6 | B | NOTE | C09 | stated as D14's result; `core.vorticity.rotating_lamb_form_terms` residual = the (5.20) residual on random fields |
| N33 | Eq. (5.26): the curl of (5.25) in index notation, $\varepsilon_{nqi}(\ )_{,q}$; the gradient term vanishes (text's "Π" is a slip for Φ) [#58] | 5.6 | B | NOTE | C09 | stated as steps 1–4 of D15 |
| N34 | Eq. (5.27): nonlinear + Coriolis term by ε–δ, $=-u_{n,j}(\omega_j+2\Omega_j)+u_j\omega_{n,j}$ (the book silently drops $u_{j,j}(\omega_n+2\Omega_n)$) [#59] | 5.6 | B | NOTE | C09 | stated as steps 5–9 of D15 with the dropped term written and killed by (5.19) |
| N35 | Eq. (5.28): curl of the pressure term = the baroclinic vector $\frac1{\rho^2}[\nabla\rho\times\nabla p]_n$ [#60] | 5.6 | B | NOTE | C09 | stated as steps 10–11 of D15 (quotient rule gloss) and linked back to C04's torque (D06 derives the same vector from a disc); `core.vorticity.baroclinic_term(_sym)`: zero for ρ = ρ(p), nonzero for ρ(x), p(z) |
| N36 | Eq. (5.29): viscous term $\to\nu\omega_{n,jj}$ (ε–δ, ∇·ω = 0) [#61] | 5.6 | B | NOTE | C09 | stated as steps 12–13 of D15 |
| C09 | Eq. (5.30): vorticity equation for a Boussinesq fluid in a frame rotating at constant Ω, $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u+\frac1{\rho^2}\nabla\rho\times\nabla p+\nu\nabla^2\boldsymbol\omega$ [#62] | 5.6 | A | CORE | – | load-bearing: the backbone of geophysical vorticity dynamics — planetary stretching, baroclinic generation, PV and QG theory (Ch. 13), shear-layer and baroclinic instability (Ch. 11), stretching in turbulence (Ch. 12); (5.13) is its Ω = 0, barotropic case |
| R18 | Planetary vorticity 2Ω and absolute vorticity $\boldsymbol\omega_a=\boldsymbol\omega+2\boldsymbol\Omega$ (measured in the inertial frame) [#63] | 5.6 | B | RECAP | C11 | ch03 D13 (ω′ = ω − 2Ω), ch04 gloss ζ + f: reminded with `ch05.absolute_vorticity` (round trip with `core.kinematics.vorticity_in_rotating_frame`) and `core.rotating.coriolis_parameter` (f = 2Ω sin φ, P126 reminder) |
| N37 | Reading the terms of (5.30): rate following the particle; molecular diffusion; baroclinic generation; stretching crucial even at Ω = 0 [#64] | 5.6 | B | NOTE | C09 | stated as term bars from `core.vorticity.vorticity_budget` on three scenes (Burgers, lock exchange, rotating stretched column), each with its residual ≈ 0; colours match E7 |
| N38 | Eq. (5.31): natural coordinates (s, n, m) on a vortex line, $(\boldsymbol\omega\cdot\nabla)\mathbf u=\omega\,\partial\mathbf u/\partial s$ [#65] | 5.6 | B | NOTE | C10 | stated as D16's result; `core.vorticity.stretching_tilting_split` (stretching ∥ ω, tilting ⟂ ω, their sum = Gω) |
| C10 | Eq. (5.32): $\frac{D\omega_s}{Dt}=\omega\frac{\partial u_s}{\partial s}$ (stretching), $\frac{D\omega_n}{Dt}=\omega\frac{\partial u_n}{\partial s}$, $\frac{D\omega_m}{Dt}=\omega\frac{\partial u_m}{\partial s}$ (tilting); none in 2-D [#66] | 5.6 | A | CORE | – | load-bearing: vortex stretching is the mechanism of the turbulent energy cascade and of intense small-scale vortices (Ch. 12), the reason 2-D and 3-D turbulence differ, and — with planetary vorticity — of spin-up in converging flow (Ch. 13); the Burgers vortex is Ch. 12's fine-scale model |
| N39 | Planetary stretching/tilting (Ω = Ωe_z): $\frac{D\omega_z}{Dt}=2\Omega\frac{\partial w}{\partial z}$, $\frac{D\omega_x}{Dt}=2\Omega\frac{\partial u}{\partial z}$, $\frac{D\omega_y}{Dt}=2\Omega\frac{\partial v}{\partial z}$; stretched columns gain ω_z in the sense of Ω (Fig. 5.10) [#67] | 5.6 | B | NOTE | C11 | stated (a-D23 given, §4c) with ∂w/∂z = −∇_h·u_h: horizontal convergence spins up; `core.vorticity.planetary_vorticity_terms`; it needs *fluid* lines, not vortex lines, to stretch; the small-change check of D20 |
| C11 | Eq. (5.33): Kelvin in a rotating frame, $\frac{D\Gamma_a}{Dt}=0$, $\Gamma_a\equiv\int_A(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\mathbf n\,dA=\Gamma+2\int_A\boldsymbol\Omega\cdot\mathbf n\,dA$ (Exercise 5.10; printed "5.11"), with the fluid column $(\omega_z+2\Omega)/h$ = const [#68] | 5.6 | A | CORE | – | load-bearing: Bjerknes' circulation theorem and potential-vorticity conservation of Ch. 13 (why a column squashed over a ridge turns anticyclonic, why air moving poleward acquires anticyclonic spin, Taylor columns, Rossby waves) — the climate-dynamics entry point of the chapter |
| C12 | Point vortices: a real vortex idealised as a line vortex that moves with the flow induced by all the others, $\frac{d\mathbf x_k}{dt}=\sum_{j\ne k}\frac{\Gamma_j}{2\pi}\frac{\mathbf e_z\times(\mathbf x_k-\mathbf x_j)}{\lvert\mathbf x_k-\mathbf x_j\rvert^2}$ [#69] | 5.7 | A | CORE | – | load-bearing: the simplest vortex dynamics — pairs (Fujiwhara binary cyclones), images (C13), sheets and roll-up (C14, Ch. 11 Kelvin–Helmholtz), vortex methods (Ch. 10), point-vortex and heton models (Ch. 13), aircraft wake pairs (Ch. 14) |
| N40 | Two same-sign vortices (Fig. 5.11): $V_1=\Gamma_1/2\pi h$, $V_2=\Gamma_2/2\pi h$; the pair orbits the centre of vorticity G [#70] | 5.7 | B | NOTE | C12 | stated with D21 (centre h₁ = Γ₂h/(Γ₁+Γ₂), rate (Γ₁+Γ₂)/2πh² — Exercise 5.18); number: Γ₁ = Γ₂ = 1 m²/s, h = 1 m → V = 0.159 m/s, period 19.7 s; Γ₂ = 3Γ₁ → G is ¾h from vortex 1; `ch05.vortex_pair`, `core.biot_savart.point_vortex_evolve` period agrees |
| N41 | Equal and opposite pair (Fig. 5.12): translates at $V=\Gamma/2\pi h$; made by a paddle stroke or a knife blade in a bucket (Fig. 5.13, after Lighthill) [#71] | 5.7 | B | NOTE | C12 | stated as the Γ₁ + Γ₂ = 0 branch of D21 (G at infinity); animation of the pair; → aircraft wake vortices (Ch. 14) |
| C13 | Method of images: a vortex at distance h from a plane wall plus its opposite mirror image; the wall is a streamline and the vortex drifts parallel to it at $V_A=\Gamma/4\pi h$ (Fig. 5.14) [#72] | 5.7 | A | CORE | – | load-bearing: the tool for boundaries in ideal flow — cylinders and walls in Ch. 6, ground effect in Ch. 14, rings near walls here; circle and channel images reused in Ch. 6 and Ch. 10 |
| N42 | Vortex ring near a wall (Fig. 5.15): ring plus opposite image ring; the ring widens and approaches with decreasing speed [#73] | 5.7 | B | NOTE | C13 | stated (the book argues in words) with a computed trajectory R(t), z(t) (`core.biot_savart.ring_ring_velocity` with `scipy.special.ellipk/ellipe` — primer, parameter m = k²; self-speed Γ/4πR[ln(8R/a) − 1/4], Wikipedia "Vortex ring"); asserted monotone |
| N43 | Two coaxial rings of the same sense leap-frog forever in an ideal fluid (Exercise 5.15; Sommerfeld) [#74] | 5.7 | B | NOTE | C13 | stated with an animation (`ch05.ring_dynamics` with two rings; impulse ΣΓπR² conserved); backup explainer B1 |
| N44 | Def: vortex sheet — many parallel line vortices side by side (Fig. 5.16); tangential velocity jumps, normal velocity continuous; finite thickness spreads the jump [#75] | 5.8 | B | NOTE | C14 | stated with a discrete row of N point vortices of strength γ ds: far above/below → ∓γ/2 (`core.biot_savart.vortex_sheet_velocity`); number: γ = 2 m/s → ±1 m/s |
| C14 | Strength of a vortex sheet: $d\Gamma=u_2\,ds+v\,dn-u_1\,ds-v\,dn=(u_2-u_1)\,ds$, circulation per length = jump in tangential velocity (u₁ above, u₂ below) [#76] | 5.8 | A | CORE | – | load-bearing: the model of shear layers and mixing layers (Ch. 9, 11: Kelvin–Helmholtz instability), bound and trailing sheets of wings (Ch. 14), a wall as a vortex sheet (Exercise 5.20, Ch. 9's wall vorticity) |
| S01 | Exercises 5.1–5.20 [#77] | Ex. | C | SKIP | – | pointer line: practise on the book's exercises; the results the text defers to Exercises 5.4, 5.5, 5.9, 5.10, 5.12, 5.18 are written out in our words (§4b) and 5.6, 5.11 are stated (N22, N23); exercise text is not reproduced; 5.7–5.8 (Vazsonyi, Crocco) → Ch. 15 |
| S02 | Literature cited and supplemental reading [#78] | Lit. | C | SKIP | – | pointer line: bibliography (Lighthill 1986 on the knife-blade pair, Sommerfeld 1964; Pedlosky → Ch. 13, Saffman for vortex dynamics) |
| N45 | Fig. 5.1: stream tube vs vortex tube [#79] | 5.1 | B | NOTE | C01 | our plotly 3-D figure: both tubes traced in a Burgers vortex with the tube-flux bars along the tube (`core.vorticity.vortex_tube`, `tube_flux_budget`) |
| N46 | Fig. 5.2: rotating tank, paraboloid isobars and free surface [#80] | 5.1 | B | NOTE | C02 | our figure: isobars for three p levels + free surface in the (r, z) half-plane, `slider_figure` over the tank rate; number: tank at 5 rad/s (ω = 10 s⁻¹), R = 0.1 m → rim–centre height difference ω²R²/8g = 12.7 mm |
| N47 | Fig. 5.3: line vortex, funnel-shaped isobars [#81] | 5.1 | B | NOTE | C02 | our figure: funnel isobars of (5.7) beside the Rankine composite (finite bottom) on the same axes |
| N48 | Fig. 5.4: material contour C with the element dx carried from u to u + du [#82] | 5.2 | B | NOTE | C03 | our animation: a loop advected in a strain + vortex flow with dx and du arrows at one point and the Γ(t) trace (flat inviscid, decaying viscous) |
| N49 | Fig. 5.5: lock exchange before and after the barrier is removed [#83] | 5.2 | B | NOTE | C04 | our schematic + a colour map of ∇ρ × ∇p/ρ² on the smoothed interface with the tumbling sense |
| N50 | Fig. 5.6: barotropic vs baroclinic element (isolines, pressure arrows, net force vs centre of mass G) [#84] | 5.2 | B | NOTE | C04 | our two-disc figure: arrows ∝ p, isobars and isopycnals, the net-force line and G, torque number (`core.vorticity.pressure_torque_on_element`); `slider_figure` over the isopycnal tilt |
| N51 | Fig. 5.7: surface S on the wall of a vortex tube and its material image S′ [#85] | 5.3 | C | NOTE | C05 | named with our 3-D schematic in D08 (no separate computation) |
| N52 | Fig. 5.8: Biot–Savart geometry — vortex element at x′ in V′, field point x [#86] | 5.5 | B | NOTE | C07 | our figure: a curved filament, a field point and the contribution arrow of one segment; the sum of arrows = the velocity (E6 stage) |
| N53 | Fig. 5.9: natural coordinates (s, n, m) on a vortex line [#87] | 5.6 | B | NOTE | C10 | our plotly 3-D: a helical vortex line with its Frenet frame and the stretching/tilting arrows of a chosen G |
| N54 | Fig. 5.10: fluid column stretched from A to B in a rotating layer gains ω_z [#88] | 5.6 | B | NOTE | C11 | our figure: a column over a sloping bottom and ζ vs h with (ζ + 2Ω)/h constant; `slider_figure` over latitude |
| N55 | Fig. 5.11: same-sign pair orbiting G [#89] | 5.7 | B | NOTE | C12 | our trajectories for Γ₂/Γ₁ = 1 and 3 with G marked (`point_vortex_evolve`) |
| N56 | Fig. 5.12: opposite pair translating [#90] | 5.7 | B | NOTE | C12 | our trajectory plot (and E8 preset) |
| N57 | Fig. 5.13: knife-blade vortex pair in a bucket (after Lighthill) [#91] | 5.7 | B | NOTE | C13 | our simulation: a pair launched inside a circle (circle images, `core.biot_savart.circle_image_system`) crosses, then separates along the wall; animation |
| N58 | Fig. 5.14: vortex A near a wall with its image B; the resultant at wall points is tangent to the wall [#92] | 5.7 | B | NOTE | C13 | our quiver: wall velocities from A and B and their tangent sum (u·n = 0 at 200 points) |
| N59 | Fig. 5.15: vortex ring near a wall with its image; its widening path [#93] | 5.7 | B | NOTE | C13 | our R(t), z(t) trajectory with the image ring (N42) |
| N60 | Fig. 5.16: vortex sheet — jump profile, row of filaments, the dn × ds circuit [#94] | 5.8 | B | NOTE | C14 | our figure: u(y) across a discrete row for N = 10, 100, 1000 converging to the jump (error ∝ 1/N) with the circuit drawn |

### 2a. What each A block contains (for the lesson-designer)
- **C01** picture: smoke rings and the swirl of water down a drain — a spinning tube of fluid that never just stops
  in mid-water · question: "can a vortex end in the middle of the fluid?" · D01 · number: a tube of strength
  Γ = 0.01 m²/s squeezed from 1 cm² to 0.25 cm² must raise its mean vorticity from 100 to 400 s⁻¹ · code:
  `core.vorticity.vortex_line`, `vortex_tube`, `vortex_tube_strength`, `tube_flux_budget` + from-scratch midpoint sums
  of ω·n over two cross-sections · figure: plotly 3-D stream tube vs vortex tube (N45); flux vs position along the tube
  (flat) · explainer E1 · B/C: R01, N01, N02, N03, N45.
- **C02** picture: a spinning bucket of water (the surface dips into a bowl) and a bathtub drain or tornado (a narrow
  funnel) · question: "why is the surface a bowl in one vortex and a funnel in the other, and which one feels viscous
  forces?" · D02, D03 · number: tank at 5 rad/s (ω = 10 s⁻¹), R = 0.1 m, water → rim–centre height 12.7 mm; line vortex
  Γ = 1 m²/s at r = 0.1 m → pressure deficit 1266.5 Pa, funnel depth 0.129 m · code: `ch05.solid_body_pressure`,
  `isobar_height`, `line_vortex_pressure`, `rankine_pressure`, `bernoulli_across_vortex`, `line_vortex_viscous_stress`,
  `torque_per_length`, `rotating_tank_free_surface` + from-scratch `cumulative_trapezoid` integration of (5.5a) ·
  figure: (r, z) isobars of both vortices on shared axes with the Rankine composite (N46, N47); the three-row stress vs
  net-force table (N09) · explainer E2 · B/C: R02–R08, N04–N09, N46, N47.
- **C03** picture: a ring of dye in water that is swirled and stretched into a long, thin, tangled loop · question:
  "the loop is stretched ten times longer — what stays the same?" · D04, D05 · number: steady cellular flow
  ψ = sin x sin y (units m, m²/s): a loop stretched ×10 in 3 turnover times keeps Γ to 10⁻⁹; Lamb–Oseen Γ₀ = 0.01 m²/s,
  ν = 10⁻⁶ m²/s, r = 5 mm, t = 10 s → Γ = 0.00465 m²/s and falling (viscous) · code: `core.vorticity.material_loop`,
  `material_circulation`, `kelvin_rate_terms`, `kelvin_force_terms`, `ch05.lamb_oseen_circulation`,
  `ch05.kelvin_hypotheses` + from-scratch RK4-advected loop points and a periodic trapezoid ∮u·dx · figure: animation
  of the loop with Γ(t) (N48); the four-restriction decision table (N16) · explainer E3 · B/C: N10–N14, N16, R09, N48.
- **C04** picture: a sea breeze — warm land, cool sea, the air turns over at the coast; a lock gate pulled out between
  fresh and salt water · question: "why does a density difference sideways make fluid spin, but a density difference
  up-down does not?" · D06, D07 · number: disc of radius 1 cm in water with ∇ρ horizontal (10 kg/m⁴) and hydrostatic ∇p
  → spin-up rate ∇ρ × ∇p/ρ² ≈ 10 × 9.81 × 1000/1000² ≈ 0.098 s⁻² (torque route and formula agree to O(R²)); lock
  exchange ρ₁ = 1000, ρ₂ = 1025 kg/m³, δ = 0.1 m → 2.42 s⁻² · code: `core.vorticity.pressure_torque_on_element`,
  `baroclinic_term`, `ch05.lock_exchange_initial_vorticity_rate` + from-scratch sum of pressure forces round N boundary
  points and the torque about G · figure: the two discs (N50), lock-exchange rate map (N49), `slider_figure` over the
  isopycnal tilt · explainer E4 · B/C: N15, N49, N50.
- **C05** picture: a smoke ring drifting across a room — the ring of spin travels with the smoke that carries it ·
  question: "does the vorticity stay attached to the same fluid?" · D08 · number: a small loop lying on a tube wall of
  an inviscid strain + swirl flow keeps ∫ω·n dA ≈ 10⁻⁹ Γ while it is carried and stretched · code:
  `core.vorticity.material_loop` on the tube wall; after C06, `core.vorticity.frozen_in_check` (angle between δx and ω
  < 10⁻⁸, ∣ω∣/∣δx∣ constant) + from-scratch flux through the loop · figure: 3-D schematic of S and S′ (N51); the
  frozen-in trace (angle and ratio vs t; the viscous Burgers run drifts) · explainer E3 (Helmholtz step) · B/C: N17,
  N51.
- **C06** picture: stirring tea — the spin you made spreads and fades; a tornado tightens as its column stretches ·
  question: "what can change a fluid particle's vorticity, and why can pressure not?" · D09 ★★★ · number: Lamb–Oseen
  (stretching 0: local = diffusion), Taylor–Green (2-D), Burgers (all terms nonzero, residual ~10⁻⁶ of the largest at
  h = 10⁻³); diffusing sheet γ = 1 m/s, ν = 10⁻⁶ m²/s: width 2√(νt) = 2 mm after 1 s, ∫ω dy = γ · code:
  `core.vorticity.vorticity_equation_sym`, `vorticity_terms`, `ch05.diffusing_vortex_sheet`,
  `core.vortices.hill_spherical_vortex` + from-scratch central differences for (ω·∇)u and ν∇²ω · figure: term bars at
  probe points of three flows; `slider_figure` of ω_z(y, t) of the diffusing sheet (N22) · explainer E5 · B/C: N18,
  R10, N19, R11, N20, N22, N23.
- **C07** picture: an electric current and the magnetic field around it; a vortex in one place stirring water far
  away · question: "given the vorticity, is the velocity everywhere decided?" · D10 ★★★, D11 ★★★ · number: Gaussian
  tube (core 0.1 m, Γ = 1 m²/s): Biot–Savart at r = 0.5 m = Γ(1 − e^{−25})/2πr = 0.318 m/s (the line-vortex value); the
  printed −1/(4π) gives −0.318 m/s · code: `core.biot_savart.poisson_green_3d`, `velocity_from_curl_omega`,
  `biot_savart_volume`, `biot_savart_2d`, `velocity_from_vorticity_fft`, `core.integral_theorems.curl_theorem_box` +
  from-scratch double loop over grid cells for the 2-D kernel · figure: Fig. 5.8 geometry (N52); u_θ(r) from
  Biot–Savart vs exact for Rankine and Lamb–Oseen, inside and outside the core; `slider_figure` over the number of
  quadrature nodes · explainer E6 · B/C: N24–N28, N52.
- **C08** picture: a smoke ring pushing itself forward; a wingtip vortex dragging air down behind a plane · question:
  "how much velocity does one short piece of a thin vortex make?" · D12, D13 · number: Γ = 1 m²/s, d = 1 m → infinite
  line 0.159 m/s, semi-infinite 0.080 m/s; square loop of side 1 m → 0.900 m/s at its centre; ring R = 0.5 m →
  Γ/2R = 1 m/s at its centre · code: `core.biot_savart.segment_induced_velocity`, `filament_velocity`,
  `ring_axis_velocity` + from-scratch midpoint sum of (5.17) along a segment · figure: velocity at a point from a
  polygon ring as the number of segments grows (order 2 in 1/M); `slider_figure` over the segment angles · explainer E6
  · B/C: N29.
- **C09** picture: the atmosphere over a warm coast on the rotating Earth — spin made by density contrasts, borrowed
  from the planet, spread by friction · question: "what does the vorticity equation gain when the frame rotates and the
  density varies?" · D14, D15 ★★★ · number: a stretched column with Ω = 7.292×10⁻⁵ rad/s and ∂w/∂z = 10⁻⁵ s⁻¹ gains
  2Ω∂w/∂z = 1.46×10⁻⁹ s⁻²; the lock exchange's baroclinic term 2.42 s⁻² (from C04) · code:
  `core.vorticity.vorticity_budget(_sym)`, `rotating_lamb_form_terms`, `baroclinic_term(_sym)`,
  `vorticity_divergence`, `ch05.rotating_ns_residual` + from-scratch baroclinic vector by central differences ·
  figure: term bars (local, advective, relative stretching/tilting, planetary, baroclinic, diffusion, residual) on
  three scenes (N37) · explainer E7 · B/C: R12–R17, N30–N37.
- **C10** picture: a figure skater pulling in her arms; a tornado tightening as its column is stretched upward ·
  question: "how can a flow make vorticity stronger without any torque?" · D16, D17, D18 · number: axial strain
  α = 1 s⁻¹ multiplies ω by e² = 7.39 in 2 s; shear s = 1 s⁻¹ tilts ω_x into ω_z = s t ω_x; a tube stretched to twice
  its length doubles ω; Burgers core 2 mm for α = 1 s⁻¹, ν = 10⁻⁶ m²/s · code: `core.vorticity.stretching_tilting_split`,
  `ch05.uniform_strain_vorticity` (`scipy.linalg.expm`, P79), `ch05.stretched_tube`, `core.vortices.burgers_vortex` +
  from-scratch projection (e_s·G e_s)ω e_s and the remainder · figure: plotly 3-D helix with its frame and arrows (N53);
  ∣ω∣(t) for the three presets; Burgers ω_z(R) with its term bars · explainer E5 · B/C: N38, N21, N53.
- **C11** picture: a column of air squashed as it crosses a mountain range, and a cyclone forming where air
  converges; a ring of air drifting poleward · question: "in a spinning world, what is conserved when a column is
  stretched or moved?" · D19, D20 · number: a column at f = 10⁻⁴ s⁻¹ stretched from 1000 m to 1100 m gains ζ = 10⁻⁵ s⁻¹
  (cyclonic); a ring of radius 500 km starting at rest at 30° N, moved to 60° N at constant area, acquires
  Γ = −4.19×10⁷ m²/s, mean ζ = −5.34×10⁻⁵ s⁻¹ (anticyclonic) · code: `core.vorticity.absolute_circulation`,
  `loop_vector_area`, `planetary_vorticity_terms`, `ch05.column_relative_vorticity`,
  `ch05.relative_circulation_after_move` + from-scratch ½∮x × dx by the trapezoid rule · figure: column over a slope
  and ζ vs h (N54), `slider_figure` over latitude; Γ and Γ_a of a material loop in a rotating frame (Γ_a flat) ·
  explainer E7 · B/C: R18, N39, N54.
- **C12** picture: two hurricanes that come close and start to circle each other (Fujiwhara); a pair of dimples made by
  a paddle stroke that glide away together · question: "a vortex cannot push itself — so how do vortices move?" · D21 ·
  number: Γ₁ = Γ₂ = 1 m²/s, h = 1 m → each moves at 0.159 m/s, period 19.7 s; Γ₂ = 3Γ₁ → G is ¾h from vortex 1;
  opposite pair → 0.159 m/s straight · code: `core.biot_savart.point_vortex_velocity`, `point_vortex_evolve`,
  `point_vortex_invariants`, `ch05.vortex_pair` + from-scratch double loop for dx_k/dt · figure: trajectories (N55,
  N56); invariants ΣΓx, ΣΓ∣x∣², H flat to 10⁻¹⁰ · explainer E8 · B/C: N40, N41, N55, N56.
- **C13** picture: a vortex near the bottom of a tank sliding along it; a smoke ring blown at a wall growing wider as
  it slows · question: "how does a wall push a vortex without touching it?" · D22 · number: Γ = 1 m²/s, h = 0.5 m →
  drift 0.159 m/s; wall normal velocity 0 at 200 points · code: `core.biot_savart.wall_image_system`,
  `circle_image_system`, `channel_image_velocity`, `ring_ring_velocity`, `ring_self_velocity`, `ch05.ring_dynamics`,
  `ch05.vortex_near_wall_speed` + from-scratch image sum · figure: wall quiver (N58), bucket pair (N57), ring path
  (N59); animation of leap-frogging rings (N43) · explainer E8 (backup B1 for rings) · B/C: N42, N43, N57–N59.
- **C14** picture: the boundary between a fast and a slow stream — a thin layer of spin; the wake sheet behind a wing
  rolling up into two vortices · question: "what is a velocity jump made of?" · D23 · number: γ = 2 m/s → +1 m/s below,
  −1 m/s above (counterclockwise-positive convention); N = 100 filaments reproduce the jump to ≈ 1 % · code:
  `core.biot_savart.vortex_sheet_velocity`, `ch05.vortex_sheet_strength`, `ch05.discrete_sheet_convergence`,
  `ch05.sheet_rollup` (§8) + from-scratch sum of N point-vortex velocities · figure: u(y) for N = 10, 100, 1000 (N60);
  animation of a perturbed sheet rolling up (forward pointer to Ch. 11) · explainer E9 · B/C: N44, N60.

## 3. Section coverage
| § | Title | A | B | C | RECAP | SKIP |
|---|---|---|---|---|---|---|
| 5.1 | Introduction | C01, C02 | N02, N03, N05, N06, N07, N08, N09, N45, N46, N47 | N01, N04 | R01 (B), R02 (B), R03 (B), R04 (C), R05 (B), R06 (B), R07 (B), R08 (B) | — |
| 5.2 | Kelvin's Circulation Theorem | C03, C04 | N10, N11, N12, N13, N14, N15, N16, N48, N49, N50 | — | R09 (B) | — |
| 5.3 | Helmholtz's Vortex Theorems | C05 | N17 | N51 | — | — |
| 5.4 | Vorticity Equation in a Nonrotating Frame | C06 | N19, N20, N21 (taught in C10), N22 | N18, N23 | R10 (C), R11 (B) | — |
| 5.5 | Velocity Induced by a Vortex Filament: Law of Biot and Savart | C07, C08 | N24, N25, N26, N27, N28, N29, N52 | — | — | — |
| 5.6 | Vorticity Equation in a Rotating Frame | C09, C10, C11 | N30, N31, N32, N33, N34, N35, N36, N37, N38, N39, N53, N54 | — | R12 (C), R13 (C), R14 (B), R15 (B), R16 (B), R17 (B), R18 (B) | — |
| 5.7 | Interaction of Vortices | C12, C13 | N40, N41, N42, N43, N55, N56, N57, N58, N59 | — | — | — |
| 5.8 | Vortex Sheet | C14 | N44, N60 | — | — | S01, S02 (end of chapter) |

Totals: A 14 · B 69 (55 NOTE + 14 RECAP) · C 11 (5 NOTE + 4 RECAP + 2 SKIP) = 94. No section is empty; every section
has at least one A item. N21 (Burgers, an exercise attached to §5.4) is taught inside C10 in §5.6, where its
stretching–diffusion balance belongs; the §5.4 notebook section names it with a forward pointer. C05's field-equation
check is shown at the end of C06 (it needs (5.13)).

## 4. Prerequisites needing primers (concept or tool | needed by | why it is not A/B/RECAP)
Rows marked **gloss** are one plain sentence where the item is used; all others are full 📎 primers (one or two plain
sentences, what it means here, a 2–4-line runnable demo with easy numbers). New primers continue at **P134**. Items
already in `knowledge/primers.md` get a one-line reminder naming the earlier primer (listed at the end).

| Concept or tool | Needed by | Why it is not A/B/RECAP |
|---|---|---|
| partial integration with an unknown "constant" function: integrating ∂p/∂r in r leaves f(z), fixed by the other equation | C02 (D02) | new calculus move (ch01 integrated ODEs only) |
| polar strain rate $S_{r\theta}=\tfrac12[\frac1r\frac{\partial u_r}{\partial\theta}+r\frac{\partial}{\partial r}(\frac{u_\theta}{r})]$ and the polar divergence of a stress $(1/r^2)\partial(r^2\sigma_{r\theta})/\partial r$ | C02 (D03) | Appendix-B operators (`core.curvilinear`) — **gloss** with a sympy check |
| closed-loop integral of an exact differential is zero: $\oint dF=0$ for a single-valued F (and fails for a multi-valued one such as θ round the origin) | C03 (D04, D05), C11 (D19) | ch01 P51 named exact differentials, ch02 P86 loop integrals; the "single-valued" condition is new and decides Kelvin |
| material label: a loop parametrised by fixed particle labels s ∈ [0, 1) so D/Dt passes inside ∮ | C03 (D04) | extends P109 (fixed limits) and P89 (functions with parameters) — **gloss** |
| periodic trapezoid rule on a closed loop (spectrally accurate for smooth periodic integrands) | C03 (code), C11 (code) | numerics; ch01 P37 trapezoid had O(Δx²) on open intervals |
| advecting many points in one `solve_ivp` call (flattened state `y.reshape(3, N)`) | C03, C05, C12 (code) | Python idiom not yet primed (P31/P94 advected one point) |
| centre of mass of a body with non-uniform density; moment of inertia of a disc I = ½MR² | C04 (D06) | P116 primed torque and the cube's I; the offset centre of mass and the disc are new — **gloss** on top of P116 |
| smoothed step with `np.tanh` (an interface of thickness δ) | C04 (D07, code) | **gloss** |
| Poisson equation ∇²φ = q and its Green's function $G=-1/(4\pi\lvert\mathbf x-\mathbf x'\rvert)$: $\nabla^2(1/r)=-4\pi\delta$ (flux of ∇(1/r) through a small sphere = −4π), superposition $\phi=\int Gq\,d^3x'$ | C07 (D10) | new maths tool; the Dirac delta was a ch03 gloss (N40), now made precise |
| gradient of $1/\lvert\mathbf x-\mathbf x'\rvert$ with respect to x′ (not x): $\nabla'(1/r)=+(\mathbf x-\mathbf x')/r^3$ | C07 (D11) | the sign trap behind the book's two slips — primer with a sympy demo |
| product rule for a curl $\nabla\times(f\mathbf A)=f\nabla\times\mathbf A+\nabla f\times\mathbf A$ | C07 (D11) | new vector identity (ch04 glossed the χ∇ψ special case) |
| cross-product antisymmetry a × b = −b × a in reordering kernels | C07 (D11), C09 (D14) | **gloss** (ch02 C08 reminder) |
| Fourier modes and the FFT Poisson solver: ∂/∂x → ik, ∇² → −k², û = ik × ω̂/k² (`numpy.fft`) | C07 (N24 code) | new maths/Python tool (forward: Ch. 10 spectral methods) |
| 3-D Gauss–Legendre quadrature on a box (`core.integral_theorems.gauss_legendre_nodes`, tensor product of 1-D nodes) and a smoothed kernel ε for points inside the vorticity | C07 (code) | ch02 used midpoint sums (P83); Gauss nodes and kernel smoothing are new |
| improper integral $\int_{-\infty}^{\infty}d\,dl/(d^2+l^2)^{3/2}$ as the limit of finite integrals | C08 (D13) | new maths tool (P106 substitution is a reminder) |
| trig identity 1 + cot²θ = csc²θ for the substitution l = d cot θ | C08 (D13) | **gloss** |
| far-field ("frozen kernel") approximation: a function that varies little across a small region can be pulled out of the integral | C08 (D12) | **gloss** with one number (relative change of 1/r³ across the core) |
| quotient rule for a gradient $\nabla(1/\rho)=-\nabla\rho/\rho^2$ | C09 (D15), C04 (D06 check) | ch04 glossed the time-derivative version — **gloss** |
| Frenet frame of a curve: unit tangent e_s, principal normal e_n (the book's Fig. 5.9 points it away from the centre of curvature), second normal e_m | C10 (D16, N53) | new geometry tool |
| angular momentum of a spinning cylinder L = Iω with I ∝ mA: thinner at fixed mass ⇒ faster spin (the skater) | C10 (D17, D18) | physics vocabulary; P116 primed I but not conservation of angular momentum |
| cylindrical Laplacian of an axisymmetric scalar $\frac1R\frac{d}{dR}(R\frac{d\omega}{dR})$ | C10 (D18) | **gloss** (Appendix B, `core.curvilinear`) |
| vector area of a closed loop $\mathbf A_{\text{vec}}=\tfrac12\oint\mathbf x\times d\mathbf x$ (planar loop: area × normal) | C11 (D19) | ch02 P69 was a closed *surface*; the loop version is new |
| column mass conservation A h = const for an incompressible column | C11 (D20) | **gloss** |
| systems of ODEs for many interacting bodies and conserved quantities as accuracy checks (ΣΓx, ΣΓ∣x∣², H) | C12 (code) | new numerics idea (P31/P94 integrated one trajectory) |
| lever rule / centre of mass analogy for the centre of vorticity; circular motion V = rate × radius | C12 (D21) | **gloss** |
| mirror symmetry of a vortex and its opposite image about a wall | C13 (D22) | **gloss** |
| complete elliptic integrals `scipy.special.ellipk(m)`, `ellipe(m)` with the *parameter* m = k² (not the modulus) | C13 (N42 code) | new special functions (Python) |
| limit of a sum of N filaments → a continuous sheet, error ∝ 1/N | C14 (N44, N60) | **gloss** (Riemann sum, ch02 gloss) |
| error function erf(η) = 1 − erfc(η) | C06 (N22) | P123 primed erfc — **gloss** |
| reading comma index notation (u_{i,jk} = ∂²u_i/∂x_j∂x_k) and renaming the free index n → i at the end of a derivation | C09 (D14, D15) | ch02 §2.14 comma notation reminder — **gloss** |

Reminders only (already primed): partial derivative P25 · chain rule P49/P91 · Taylor P26/P98 · definite integral P27 ·
trapezoid P37 · product rule P38 · sympy P40, P117 · `solve_ivp` P31/P94 · RK4 by hand P95 · lambda P29 ·
`assert np.allclose` P15 · `animate` P16 · `slider_figure` P17 · `show_viz` P18 · plotly 3-D P41/P64 · `np.einsum` P62 ·
orthonormal basis and projection P65 · orders of smallness P68 · vector area of a closed surface P69 · permutations and
cyclic order of ε P72 · right-hand rule and orientation P74 · level sets and directional derivative P75 · `meshgrid` P76 ·
broadcasting P77 · contour/quiver/streamplot P78 · `expm` P79 · midpoint sums P83 · FTC P84 · line integral round a loop
P86 · `quad`/`dblquad` P87 · cylindrical unit vectors P88 · parametric curves and arc length P92 · parallel vectors
P93 · material line element P99 · rigid-body velocity P101 · rotating frame P103 · polar coordinates P105 · substitution
P106 · `np.expm1` P107 · `brentq` P108 · differentiation under the integral sign P109 · dataclasses P111 · small-ball
localisation P112 · Schwarz P121 · curl of a curl P122 · erfc P123 · rotating unit vector P124 · conservative force
P115 · torque and moment of inertia P116 · latitude and f = 2Ω sin φ P126 · order-of-magnitude scaling P130 · reduced
gravity P131.

## 4b. Derivations written out (parsed by tools: ID first, CORE id in a column, ★★★ for hard, explainer slugs backticked in the LAST column)
Written out because (a) the result belongs to an A item, or (b) the book never writes it out and the lesson needs it
(marked "(b) book never writes it out"). One small move per step; the book's skipped moves (analysis §2b) are filled in
and the notebook says so. Analysis §2b numbers in brackets (`a-Dnn`). 23 rows.

| ID | Result (Eq.) | CORE | Difficulty | Steps | Tools used | Traps | Shown in |
|---|---|---|---|---|---|---|---|
| D01 | tube strength is the same at every cross-section; tubes cannot end (5.4) [a-D02] | C01 | ★ | 6 | Gauss (2.30) (recap ch02 C14), ∇·ω = 0 (R10), Stokes (2.34) for flux = circulation (recap ch02 C16), outward normal and orientation (P74) | the outward normal on the lower end points against ω, hence −Γ_lower; the side carries no flux because it is made of vortex lines (ω·n = 0), not because ω = 0 there; "cannot end" is kinematic — it holds in viscous flow too | notebook · `vortex_tubes_cannot_end` |
| D02 | Euler in cylindrical coordinates → (5.5a, b) → the paraboloid pressure (5.6) [a-D03, a-D04] | C02 | ★★ | 10 | cylindrical unit vectors and de_θ/dθ = −e_r (P88, P124 reminders), S = 0 ⇒ Euler (R05), partial integration with an unknown function (primer), hydrostatics (R06) | the −u_θ²/r term comes from the turning unit vector, not from ∂u_θ/∂θ; ω is the vorticity so u_θ²/r = ω²r/4 (not ω²r); integrating ∂p/∂r leaves f(z), fixed by ∂p/∂z; Fig. 5.2's "2ω" label (the tank turns at ω/2) | notebook · `vortex_pressure_funnel` |
| D03 | line vortex: σ_rθ = −μΓ/πr² ≠ 0 but zero net viscous force, three routes — (b) book never writes it out (Exercise 5.4) [a-D06] | C02 | ★★ | 9 | polar strain rate and polar stress divergence (gloss), σ = 2μS (recap ch04 C07), viscous force −μ∇×ω (recap ch04 D13) | dropping the 1/r² metric factor makes the net force look nonzero; the stress is nonzero because elements deform (ch03), the force is zero because the stresses on the two sides of an element balance; valid for incompressible constant-μ flow | notebook · `vortex_pressure_funnel` |
| D04 | rate of change of the circulation of a material loop (5.9), with ∮u·du = 0 [a-D09] | C03 | ★★ | 8 | circulation (3.18) (recap ch03), material line element D(δx)/Dt = δu (P99), material label (gloss) and differentiation under the integral with fixed limits (P109), closed-loop exact differential (primer) | D/Dt passes inside ∮ only because the loop is parametrised by fixed particle labels (the domain does not move); u(x + dx) = u + du is a first-order Taylor step (P98); ∮d(½u²) = 0 needs single-valued u | notebook · `kelvin_material_loop` |
| D05 | (5.9) + momentum → (5.10) → (5.11) → Kelvin (5.8), and the three sources [a-D10] | C03 | ★★ | 9 | Cauchy/NS (recap ch04 C06, C08), conservative force g = −∇Φ (P115), barotropic pressure function (R09), closed-loop exact differential (primer), Stokes (recap ch02 C16) | "ρ and p single-valued" is not enough — ∮dp/ρ = 0 needs ρ = ρ(p) (the book's reason is incomplete); a baroclinic field gives ∮dp/ρ = −∫(∇ρ×∇p/ρ²)·n dA; "C in irrotational fluid ⇒ no viscous term" holds for incompressible constant-μ flow only; σ_ij differentiated on its second index is harmless because σ is symmetric | notebook · `kelvin_material_loop` |
| D06 | pressure torque on a small element equals the baroclinic vorticity rate ∇ρ×∇p/ρ² (Fig. 5.6, the same vector as (5.28)) — (b) book never writes it out [a-D29] | C04 | ★★ | 10 | net pressure force (P28), torque and moment of inertia (P116, disc I = ½MR² gloss), centre of mass (gloss), vorticity = twice the spin (recap ch03 C10) | for linear p the pressure force passes through the disc's geometric centre while the centre of mass is shifted toward the heavy side by R²∇ρ/4ρ₀; Dω/Dt is twice the angular acceleration; the order ∇ρ × ∇p (not ∇p × ∇ρ) | notebook · `baroclinic_torque` |
| D07 | lock-exchange initial vorticity rate 2(ρ₂ − ρ₁)g/((ρ₂ + ρ₁)δ) — (b) book never writes it out (Exercise 5.5) [a-D30] | C04 | ★★ | 7 | D06's result, hydrostatics (recap ch01 C20), cross products of e_x, e_y (P74), smoothed step (gloss) | at t = 0⁺ the pressure is still hydrostatic with the mean density (ρ₁ + ρ₂)/2; the sense: the heavy side slumps under the light side; the rate grows as δ → 0 (a vortex sheet forms) | notebook · `baroclinic_torque` |
| D08 | Helmholtz's theorems: vortex lines move with the fluid (Fig. 5.7) from Kelvin; 2–4 from (5.4) and (5.8) [a-D11] | C05 | ★★ | 9 | Kelvin (C03), tube-wall surfaces with ω·n = 0 (C01), Stokes (recap ch02 C16), localisation "for every S" (P112) | the argument needs zero circulation on *every* surface lying on the tube wall; the thin-tube limit turns tubes into lines; all four restrictions of Kelvin must hold; the field-equation route needs (5.13) and is shown after C06 | notebook · `kelvin_material_loop` |
| D09 | vorticity equation by the vector route: curl of (4.39b) (5.12) → (5.13) [a-D12] | C06 | ★★★ | 12 | curl of a gradient = 0 (recap ch02 D26), Schwarz (P121), Lamb identity (R11), ∇×(ω×u) identity (N20; ε–δ recap ch02 D09), curl of a curl (P122), sympy check cell | the ½ in the Lamb identity (the book prints ∇(u·u)); ∇×(∇p/ρ) = 0 needs ρ constant (else the baroclinic term (5.28)); ∇×(ω×u) = −∇×(u×ω); ∇·u = 0 and ∇·ω = 0 each drop a term of (B.3.10); regrouping ∂ω/∂t + (u·∇)ω = Dω/Dt | notebook · `vorticity_stretching_tilting` |
| D10 | ∇×ω = −∇²u (Poisson) and its Green's-function solution (5.14) with the correct sign +1/(4π) — (b) book never writes it out (Exercise 5.9) [a-D13, a-D14] | C07 | ★★★ | 12 | curl of a curl (P122), Poisson equation and Green's function of ∇² (primer), divergence theorem on a small sphere (recap ch02 C14), superposition (gloss), sympy check cell | the book prints −1/(4π): with source q = −∇′×ω and G = −1/(4π∣x−x′∣) the factor is +1/(4π); ∇·u = 0 is what drops ∇(∇·u); unbounded domain with decay at infinity; any harmonic field may be added ("vorticity-induced part", Ch. 6) | notebook · `biot_savart_filament` |
| D11 | (5.14) → (5.15) → Biot–Savart (5.16) [a-D15] | C07 | ★★★ | 12 | product rule for a curl (primer), gradient of 1/∣x−x′∣ with respect to x′ (primer), Gauss in curl form (N27), cross-product antisymmetry (gloss), sympy check cell | ∇′(1/∣x−x′∣) = +(x−x′)/∣x−x′∣³ (derivative with respect to x′); the book's second sign slip cancels the first — we write correct signs throughout and say where the book differs; V′ chosen so that n × ω = 0 on the ends and ω = 0 on the side | notebook · `biot_savart_filament` |
| D12 | filament law (5.16) → (5.17) [a-D16] | C08 | ★★ | 6 | tube strength Γ = ∫ω·n dA (N03), frozen-kernel far-field approximation (gloss) | the field point must be far from the core compared with its radius; e_ω constant across a thin tube; on a curved filament itself the self-induced velocity diverges (a core size is needed, N42) | notebook · `biot_savart_filament` |
| D13 | straight segment (Γ/4πd)(cos θ_a − cos θ_b) and the infinite line Γ/2πd = (5.2) — (b) book never writes it out [a-D17] | C08 | ★★ | 8 | substitution in an integral (P106), improper integral (primer), 1 + cot² = csc² (gloss), right-hand rule (P74) | ∣e_ω × (x−x′)∣ = d for every element (distance × sine of the angle); the limits θ_a → 0, θ_b → π; a semi-infinite line gives half; the direction by the right-hand rule | notebook · `biot_savart_filament` |
| D14 | rotating-frame momentum (5.20) → Lamb form (5.25) via (5.21)–(5.24) [a-D19] | C09 | ★★ | 8 | comma notation (R12, gloss), ε–δ (recap ch02 D09), Lamb identity (R15), dummy relabelling (ch02 N43), effective-gravity potential (recap ch04 D18) | u_{j,ij} = (u_{j,j})_{,i} = 0 needs ∇·u = 0; relabelling j ↔ k flips the sign of ε_ijk in (5.24); g includes the centrifugal potential | notebook · `vorticity_equation_rotating` |
| D15 | curl of (5.25) → (5.26) → (5.27)–(5.29) → the full vorticity equation (5.30) [a-D20] | C09 | ★★★ | 15 | ε–δ with the cyclic reorder ε_nqi = ε_inq (P72), quotient rule for 1/ρ (gloss), symmetric × antisymmetric = 0 (ch02 N61), ∇·ω = 0 (N30), sympy check cell | the dropped term u_{j,j}(ω_n + 2Ω_n) = 0 (the book skips it); sign flips when terms move right; ρ is kept inside the pressure term although the fluid is "Boussinesq" (to leading order ∇ρ×∇p/ρ² ≈ ∇ρ′×g/ρ₀, cf. (4.86)); Ω = 0 with ∇ρ ∥ ∇p must give back (5.13) | notebook · `vorticity_equation_rotating` |
| D16 | natural coordinates on a vortex line: (ω·∇)u = ω ∂u/∂s (5.31) [a-D21] | C10 | ★★ | 6 | Frenet frame (primer), directional derivative (P75), projection onto a unit vector (P65) | ω has no n or m component because e_s is along ω; ∂u_s/∂s is the stretching rate of a material element along e_s (ch03 C07); ∂u_n/∂s turns the line about m, ∂u_m/∂s about n | notebook · `vorticity_stretching_tilting` |
| D17 | component equations (5.32) and the angular-momentum reading [a-D22] | C10 | ★ | 5 | projection (P65), angular momentum of a spinning cylinder (primer) | projecting on the moving frame is an instantaneous statement (the frame itself turns); in 2-D ω ⟂ plane and nothing varies along it, so no stretching or tilting | notebook · `vorticity_stretching_tilting` |
| D18 | stretched tube ω ∝ L (inviscid) and the Burgers stretching–diffusion balance — (b) book never writes it out (Exercise 5.12) [a-D31] | C10 | ★★ | 9 | Helmholtz 4 (C05), incompressibility A L = const (gloss), cylindrical Laplacian of ω_z(R) (gloss), exponential ODE (P44) | ωA stays constant only without viscosity; the book's α is twice Wikipedia's (same flow); the core radius √(4ν/α) is set by the balance, not by the initial condition | notebook · `vorticity_stretching_tilting` |
| D19 | Kelvin in a rotating frame (5.33), the absolute circulation — (b) book never writes it out (Exercise 5.10) [a-D24] | C11 | ★★ | 10 | D04–D05 moves, vector area of a loop (primer), triple product (Ω×u)·dx = Ω·(u×dx) (recap ch02 D09), Stokes ∮(Ω×x)·dx = 2Ω·A_vec (recap ch02 C16), closed-loop exact differential (primer) | the Coriolis loop integral is not zero — it equals −2Ω·dA_vec/dt; equivalently use the inertial velocity u + Ω×x; Γ_a − Γ = 2Ω·A_vec, not 2ΩA unless the loop is horizontal and Ω vertical | notebook · `vorticity_equation_rotating` |
| D20 | fluid column in a rotating layer (Fig. 5.10): (ω_z + 2Ω)/h = const — (b) book never writes it out [a-D25] | C11 | ★★ | 6 | (5.33) (C11), column mass A h = const (gloss), f = 2Ω sin φ (P126) | Γ_a ≈ (ω_z + 2Ω)A needs a thin loop with nearly uniform ω_z; the column must stay vertical (Taylor–Proudman, Ch. 13); the small-change limit agrees with Dω_z/Dt = 2Ω∂w/∂z (N39) | notebook · `vorticity_equation_rotating` |
| D21 | two line vortices: centre of vorticity h₁ = Γ₂h/(Γ₁ + Γ₂), orbit rate (Γ₁ + Γ₂)/2πh², opposite pair translating at Γ/2πh — (b) book never writes it out beyond V₁, V₂ (Exercise 5.18) [a-D26] | C12 | ★ | 6 | line vortex (R03), lever rule (gloss), circular motion (gloss) | a straight line vortex does not move itself; G stays fixed because ΣΓ_k x_k is conserved; Γ₁ + Γ₂ = 0 puts G at infinity (translation, not rotation) | notebook · `point_vortex_lab` |
| D22 | wall image: u·n = 0 on the wall and drift Γ/4πh [a-D27] | C13 | ★ | 5 | mirror symmetry (gloss), line vortex (R03), superposition (gloss) | the image sits at distance 2h, so its velocity at A is Γ/2π(2h); the drift is parallel to the wall, its sense set by the sign of Γ; images enforce no penetration only — the vortex slips along the wall (inviscid) | notebook · `point_vortex_lab` |
| D23 | strength of a vortex sheet dΓ = (u₂ − u₁)ds [a-D28] | C14 | ★ | 5 | circulation round a rectangle (recap ch02 D26), orientation (P74) | counterclockwise circuit with u₁ above and u₂ below gives u₂ − u₁ (Fig. 5.16's caption prints u₁ − u₂, the clockwise magnitude); v continuous so the dn sides cancel; dn → 0 | notebook · `vortex_sheet_rollup` |

Counts: ★ 5 (D01, D17, D21, D22, D23) · ★★ 14 · ★★★ 4 (D09 in `vorticity_stretching_tilting`, D10 and D11 in
`biot_savart_filament`, D15 in `vorticity_equation_rotating`; each has a sympy check cell that re-runs the
construction: D09 takes the curl of each NS term for a generic divergence-free polynomial field and checks the
regrouped (5.13); D10 checks ∇²(1/r) = 0 off the origin and the small-sphere flux −4π, then applies the ± factor to a
Gaussian tube; D11 checks ∇′(1/r) and the product-rule identity of N26 symbolically; D15 re-runs (5.26)–(5.29) with
sympy's LeviCivita for a generic field and compares with (5.30)). Written out under rule (b): D03, D06, D07, D10, D13,
D18, D19, D20, D21 (9); under rule (a) only: the other 14. Analysis rows merged into one written-out derivation:
a-D03 + a-D04 → D02, a-D13 + a-D14 → D10.

## 4c. Derivations demoted to statements (first column is the A parent in bold, e.g. **C20** — never a bare ID or a D id: the parser reads a bare first-cell ID as an item and blanks its tier)
Results **given, not derived**, inside the named A block (a paragraph, the equation, a number). 6 rows.

| A parent | Analysis §2b item | Result stated (Eq.) | Stated in (B item) | Why not written out |
|---|---|---|---|---|
| **C01** | a-D01 vortex-line equation (5.3) | $dx/\omega_x=dy/\omega_y=dz/\omega_z$, parametric $d\mathbf x/ds=\boldsymbol\omega/\lvert\boldsymbol\omega\rvert$ | N02 | the streamline construction (3.7) with ω for u; ch03 P92/P93 already primed the arc-length and cross-product forms — one sentence and a traced line |
| **C02** | a-D05 $-\tfrac12u_\theta^2+gz+p/\rho$ = const; B varies across streamlines | $B-B(0)=\omega^2r^2/4$ in the tank, uniform for the line vortex | R07 | a one-line rewrite of (5.6) with ρω²r²/8 = ρu_θ²/2; stated with the ch04 (4.71)–(4.72) reminder and computed B(r) |
| **C02** | a-D07 line-vortex pressure (5.7) | $p-p_\infty=-\rho\Gamma^2/8\pi^2r^2-\rho gz$ | N07 | exactly D02's moves with u_θ = Γ/2πr; the one new integral ∫ρΓ²/(4π²r³)dr = −ρΓ²/(8π²r²) and the constant from p → p_∞ are stated in one line |
| **C02** | a-D08 rotating cylinder: Γ = πa²ω, torque −2μΓ at every r, dissipation = wall work | $2\pi r^2\sigma_{r\theta}=-2\mu\Gamma$, $\int_a^R2\pi r\rho\varepsilon\,dr=\frac{\mu\Gamma^2}{\pi}(\frac1{a^2}-\frac1{R^2})$ | R08, N08 | the book states it in words and Ch. 8 (8.11) derives the flow; the formulas are stated with a numeric energy balance (quad) — power = torque × angular speed |
| **C09** | a-D18 (5.18) ∇·ω = 0 | $\omega_{i,i}=\varepsilon_{inq}u_{q,ni}=0$ | N30 | a SEEN fact (ch02 C11) proved by the ch02 N61 contraction lemma in one line |
| **C11** | a-D23 planetary component equations | $D\omega_z/Dt=2\Omega\,\partial w/\partial z$, $D\omega_x/Dt=2\Omega\,\partial u/\partial z$, $D\omega_y/Dt=2\Omega\,\partial v/\partial z$ | N39 | one line: (2Ω·∇)u = 2Ω∂u/∂z for Ω = Ωe_z; used as D20's small-change check |

## 5. Interactive explainers (5–10 + backup)
Nine explainers on the A items where manipulation or motion teaches most; one backup. Each has the required live
**Explain** tab ("Explanation & interpretation", numbered sections computing every number on screen with the reader's
settings, ending in "Reading the current setting"), a synced **Code** tab, a **Derivation** tab for its D ids, a
4–8-step walkthrough, ≥ 3 check questions and ≥ 2 more depth features. Physics mirrors scalar-callable `fluidpy`
functions (parity rows). Colours: relative vorticity teal, planetary amber, baroclinic orange, diffusion rose,
stretching purple, tilting blue; pressure orange, density blue; positive (counterclockwise) vortices teal, negative
rose. Not duplicated from earlier chapters, linked instead: ch02 `stokes_circulation_loop` (from E3: circulation = flux
at one instant; E3 adds time and material loops), ch03 `vortex_paddle_wheels` (from E2: rotating vs irrotational
vortex kinematics; E2 adds pressure and stress), ch04 `which_bernoulli` (from E2: B across streamlines),
`rotating_frame_coriolis` (from E7: the frame; E7 adds vorticity and circulation), `navier_stokes_term_balance` (from E5
and E7: the term-bar pattern, now for the curl).

### E1 · vortex_tubes_cannot_end
- **A:** C01 (also shows R01, N02 (5.3) vortex lines, N03 tube strength, N45 stream tube vs vortex tube)
- **Confusion removed:** "a vortex can fade out along its length" — the flux through every cross-section of a tube is
  the same, so where a tube thins its vorticity must grow; the tube can only end on a wall or close on itself.
- **Why interactive:** the reader slides a cross-section along a curving, narrowing tube (3-D, orbitable) and watches
  its area shrink, the mean vorticity grow and the flux stay pinned; a "cut the tube" preset makes a field with
  ∇·ω ≠ 0 and the budget visibly fails — a static figure cannot show a quantity *staying* constant as you move.
- **Stage:** (1) 3-D view of the tube (vortex lines seeded on a circle in a Burgers vortex or a twisted Gaussian tube)
  with the moving cross-section disc; (2) flux, area and mean vorticity vs position along the tube (flux flat, ω ∝
  1/A); (3) (hidePortrait) the Gauss budget bars lower end / side / upper end / total.
- **Controls:** tube: Burgers / twisted Gaussian / ring (closed) / "broken" field · section position s (drag or
  transport) · seed radius · strain α (Burgers) · view: vortex tube / stream tube.
- **Equations shown (live):** (5.3) $dx/\omega_x=dy/\omega_y=dz/\omega_z$; $\Gamma=\oint_C\mathbf u\cdot d\mathbf x=\int_A\boldsymbol\omega\cdot\mathbf n\,dA$;
  (5.4) $\int_V\nabla\cdot\boldsymbol\omega\,dV=-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0$.
- **Mirrors:** `core.vorticity.vortex_tube`, `tube_flux_budget`, `vortex_tube_strength`, `core.vortices.burgers_vortex_field`.
- **Derivations:** D01 (the Gauss step splits the surface into the three coloured pieces; the "side has no flux" step
  draws ω tangent to the side).
- **Depth features:** Explain tab (section area → ω·n averaged → flux → comparison with the other end → budget →
  reading the current setting), synced Code tab, + 3-D view, linked views (3), transport (scrub the section along the
  tube), term bars (lower/side/upper), presets, status ("✅ same strength everywhere" / "⚠️ ∇·ω ≠ 0: not a vorticity
  field").
- **Follows reference:** `np_resnet_3d.html` (orbitable 3-D scene, a token travelling through, inspector card synced to
  the selection) with `forced_damped_vibrations.html`'s explanation panel.
- **Aha:** a vortex tube is like a garden hose that carries a fixed amount of spin — squeeze it and the spin speeds up;
  it cannot just stop.

### E2 · vortex_pressure_funnel
- **A:** C02 (also shows R02 (5.1), R03 (5.2), N05 (5.5a), R06 (5.5b), R07 B across streamlines, N06 σ_rθ, N07 (5.7),
  R08 Rankine, N08 torque to infinity, N09 stress vs net force; links to ch03 `vortex_paddle_wheels`, ch04 `which_bernoulli`)
- **Confusion removed:** "the vortex with rotating fluid is the viscous one" — solid-body rotation has *no* viscous
  stress, the irrotational vortex has stress but *no net* viscous force; and the surface shape (bowl vs funnel) comes
  from the same centripetal balance with different u_θ(r).
- **Why interactive:** drag the rotation rate or Γ and the core radius and watch the free surface and isobars change
  shape on one set of axes, while a probe element shows its deformation, σ_rθ and the net force; toggling
  "Bernoulli along/across" shows B flat or growing — a tank, a tornado and a Rankine composite in one stage.
- **Stage:** (1) (r, z) cross-section of the tank or funnel with isobars and the free surface (colour = p), a probe
  element; (2) u_θ(r), B(r) and σ_rθ(r) profiles (teal, amber, rose); (3) (hidePortrait) balance bars at the probe:
  centripetal −u_θ²/r vs −(1/ρ)∂p/∂r, viscous force (zero or not).
- **Controls:** vortex: solid body / line vortex / Rankine / rotating cylinder · ω (s⁻¹) or Γ (m²/s) · core radius a ·
  probe radius (drag) · fluid (water, air).
- **Equations shown (live):** (5.1) $u_\theta=\omega r/2$; (5.2) $u_\theta=\Gamma/2\pi r$; (5.5a)
  $-\rho u_\theta^2/r=-\partial p/\partial r$; (5.5b) $0=-\partial p/\partial z-\rho g$; (5.6)
  $p-p_o=\tfrac18\rho\omega^2r^2-\rho gz$; (5.7) $p-p_\infty=-\rho\Gamma^2/8\pi^2r^2-\rho gz$;
  $\sigma_{r\theta}=-\mu\Gamma/\pi r^2$.
- **Mirrors:** `ch05.solid_body_pressure`, `line_vortex_pressure`, `rankine_pressure`, `isobar_height`,
  `bernoulli_across_vortex`, `line_vortex_viscous_stress`, `rotating_tank_free_surface`, `ch05.vortex_pressure_scenario` (§8).
- **Derivations:** D02 (the turning-unit-vector step draws e_θ rotating at the probe; the integration steps build the
  paraboloid), D03 (the metric-factor step shows the stresses on the two faces of the probe element balancing).
- **Depth features:** Explain tab (u_θ at the probe → centripetal acceleration → ∂p/∂r → pressure at r and the surface
  height → B there → σ_rθ and net force → reading the current setting), synced Code tab, + linked views (3), presets
  (bucket at 5 rad/s, bathtub drain, tornado Γ = 10⁴ m²/s with a 50 m core, rotating cylinder), term bars (radial
  balance), inspector (click a point: p arithmetic), status ("🥣 rotational: B grows outward, no viscous stress" /
  "🌪 irrotational: stress ≠ 0, net force 0" / "Rankine core edge").
- **Follows reference:** `angular_frequency_explorer_1.html` (modes on one stage, "Right now" notes, a table of real
  vortices with the current row highlighted).
- **Aha:** pressure falls toward any vortex's centre to supply the turn — slowly and bowl-shaped in rigid rotation,
  steeply and funnel-shaped in the free vortex, where the fluid is sheared yet feels no net friction.

### E3 · kelvin_material_loop
- **A:** C03, C05 (also shows N10 (5.9), N11 D(dx)/Dt = du, N12 (5.10), R09 barotropic, N13 (5.11), N14 three
  sources, N16 four restrictions, N17 Helmholtz proof, N48 Fig. 5.4; links to ch02 `stokes_circulation_loop`)
- **Confusion removed:** "circulation is conserved around any loop" — only around a *material* loop, and only when the
  flow is inviscid, barotropic, conservatively forced and seen from an inertial frame; break one hypothesis and a
  specific, computable term makes Γ change.
- **Why interactive:** the loop is carried and stretched by the flow on one clock while Γ(t) is plotted — flat in the
  inviscid cellular flow however tangled the loop gets, decaying in Lamb–Oseen, growing in a baroclinic field; the
  reader toggles each hypothesis and sees which term of (5.10) wakes up. A fixed loop (a mode) shows Γ changing even
  in inviscid flow — the "material" in the statement is felt.
- **Stage:** (1) the flow (streamlines or density/pressure isolines) with the material loop, its points and one dx/du
  pair; (2) Γ(t) with the flat ghost of Γ(0) and the term curves ∮dp/ρ, ∮dΦ, viscous integral; (3) (hidePortrait) the
  four-restriction decision table with the current row lit, or the Helmholtz mode: a small loop on a tube wall with its
  flux ≈ 0.
- **Controls:** flow: cellular (inviscid) / Rankine (loop straddling the core) / Lamb–Oseen (viscous) / baroclinic
  (lock-exchange field) / rotating frame · loop size and position (drag) · loop type: material / fixed · time
  (transport) · ν (log).
- **Equations shown (live):** (5.8) $D\Gamma/Dt=0$; (5.9) $\frac{D\Gamma}{Dt}=\oint_C\frac{Du_i}{Dt}dx_i+\oint_Cu_i\frac{D}{Dt}(dx_i)$;
  (5.10) $\ldots=-\oint_C\frac1\rho dp-\oint_Cd\Phi+\oint_C(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j})dx_i$;
  (5.11) $\frac{D\Gamma}{Dt}=\oint_C(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j})dx_i$; Lamb–Oseen
  $\Gamma(r,t)=\Gamma_0(1-e^{-r^2/4\nu t})$.
- **Mirrors:** `core.vorticity.material_loop`, `material_circulation`, `kelvin_rate_terms`, `kelvin_force_terms`,
  `ch05.lamb_oseen_circulation`, `ch05.kelvin_hypotheses`, `core.vorticity.kelvin_scenario` (§8).
- **Derivations:** D04 (the material-label step freezes the loop's labels; the ∮u·du = 0 step shows the contour term
  bar at 0), D05 (each substitution step lights one of the three term curves; the barotropy step sets the baroclinic
  preset as the counterexample), D08 (Helmholtz mode: the loop on the tube wall keeps zero flux).
- **Depth features:** Explain tab (loop length and area now → Γ by the trapezoid rule with the numbers → each term of
  (5.10) → dΓ/dt → which hypothesis is broken → reading the current setting), synced Code tab, + linked views (3),
  transport (end-of-run card "loop 10× longer, Γ unchanged to 10⁻⁹"), modes (material vs fixed loop; Kelvin vs
  Helmholtz), presets (five flows), status verdict ("✅ Kelvin holds" / "🌡 baroclinic: ∮dp/ρ ≠ 0" / "🍯 viscous: Γ
  decays" / "🔄 rotating frame: use (5.33)").
- **Follows reference:** `forced_damped_vibrations.html` (one time slider for the phenomenon and the graph, ghost
  reference curve, regime-dependent explanation).
- **Aha:** spin can be neither created nor destroyed inside an ideal barotropic fluid — the loop may be stretched into
  spaghetti, its circulation does not budge; every way of making vorticity is one broken hypothesis.

### E4 · baroclinic_torque
- **A:** C04 (also shows N15 lock exchange, N49 Fig. 5.5, N50 Fig. 5.6, N35 (5.28) as a forward link, N14 the
  baroclinic source)
- **Confusion removed:** "heavy fluid simply sinks" — when isopycnals are tilted against isobars the pressure force
  misses the centre of mass and *spins* the fluid; with parallel isolines (barotropic) there is no torque at all, no
  matter how large the density contrast.
- **Why interactive:** rotate the isopycnals relative to the isobars and watch the net-force line move away from the
  centre of mass G, the torque grow as the sine of the angle and the spin-up rate match ∇ρ × ∇p/ρ²; shrink the disc and
  the two routes converge (O(R²)); the lock-exchange mode shows the same torque turning a whole interface over.
- **Stage:** (1) a fluid disc with isobars (orange) and isopycnals (blue), pressure arrows round its rim, the geometric
  centre, G and the net-force line, a curved torque arrow; (2) spin-up rate vs tilt angle (sine curve) with the current
  angle marked and the formula ghost; (3) (hidePortrait) lock exchange: two fluids side by side, the baroclinic-term
  colour map at the interface and the initial rate.
- **Controls:** tilt angle between ∇ρ and ∇p · ∣∇ρ∣ · disc radius R (log) · mode: element / lock exchange / sea breeze
  (ours) · interface thickness δ (lock exchange).
- **Equations shown (live):** $\nabla\rho\times\nabla p=0\Leftrightarrow$ barotropic; torque $=\pi R^4(\nabla\rho\times\nabla p)/4\rho_0$,
  $I_G=\pi\rho_0R^4/2$; $\frac{D\omega}{Dt}=\frac1{\rho^2}\nabla\rho\times\nabla p$ (5.28); lock exchange
  $2(\rho_2-\rho_1)g/((\rho_2+\rho_1)\delta)$.
- **Mirrors:** `core.vorticity.pressure_torque_on_element`, `baroclinic_term`, `ch05.lock_exchange_initial_vorticity_rate`,
  `ch05.lock_exchange_fields`.
- **Derivations:** D06 (the centre-of-mass step moves G; the torque step draws the lever arm; the ×2 step converts
  angular acceleration to vorticity rate), D07 (sets the lock-exchange mode; the hydrostatic step draws ∇p with the
  mean density).
- **Depth features:** Explain tab (∇ρ and ∇p with numbers → their cross product → centre-of-mass offset → torque →
  I_G → spin-up rate by both routes → sense of rotation → reading the current setting), synced Code tab, + linked views
  (3), presets (barotropic 0°, sea breeze, fresh/salt lock exchange, stable stratification), modes (element / lock
  exchange), status ("⚖ barotropic: no torque" / "🌀 baroclinic: spins counterclockwise at … s⁻²"), inspector (click a
  rim point: its pressure force and moment).
- **Follows reference:** `amplitude_phase_second_order_II_3.html` (several windows on one state, a numbered live
  derivation in the explanation).
- **Aha:** vorticity is born where surfaces of equal density cross surfaces of equal pressure — the pressure pushes
  through the middle, the weight hangs off-centre, and the element turns.

### E5 · vorticity_stretching_tilting
- **A:** C06, C10 (also shows N19 (5.12), R11 Lamb curl, N20 (B.3.10), N22 diffusing sheet, N38 (5.31), N21 Burgers,
  N53 Fig. 5.9; the frozen-in check of C05; term-bar pattern from ch04 `navier_stokes_term_balance`)
- **Confusion removed:** "vorticity only spreads and decays" — in 3-D a flow can amplify vorticity with no torque at all
  by stretching the vortex line, and redirect it by tilting; in 2-D neither can happen, which is why 2-D and 3-D
  turbulence differ.
- **Why interactive:** a vortex line segment (and a material element) sits in a chosen linear flow G; pressing play
  stretches or tilts it while ∣ω∣ grows, and the stretching (purple) and tilting (blue) arrows update; switching on ν
  brings in diffusion (rose) and the Burgers preset reaches the steady balance — the reader sees the three terms of
  (5.13) fight in real time.
- **Stage:** (1) 3-D view of the vortex line (helix or straight segment) with the natural frame e_s, e_n, e_m and the
  stretching/tilting arrows, the material element δx riding along (parallel to ω when ν = 0); (2) ∣ω∣(t) and the angle
  of ω vs t with ghosts (e^{αt}, linear tilt); (3) (hidePortrait) term bars of (5.13) (local, advective, stretching,
  tilting, diffusion, residual) or Burgers' ω_z(R) profile.
- **Controls:** flow preset: axial strain / shear tilt / planar (2-D) / Burgers / custom G · strain rate α or s · ν
  (log, 0 allowed) · time (transport) · initial ω direction (angle).
- **Equations shown (live):** (5.12) $\nabla\times\{\frac{D\mathbf u}{Dt}=-\frac1\rho\nabla p+\mathbf g+\nu\nabla^2\mathbf u\}$;
  (5.13) $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$; (5.31)
  $(\boldsymbol\omega\cdot\nabla)\mathbf u=\omega\,\partial\mathbf u/\partial s$; (5.32) $\frac{D\omega_s}{Dt}=\omega\frac{\partial u_s}{\partial s}$,
  $\frac{D\omega_n}{Dt}=\omega\frac{\partial u_n}{\partial s}$, $\frac{D\omega_m}{Dt}=\omega\frac{\partial u_m}{\partial s}$;
  Burgers $\omega_z=\frac{\alpha\Gamma}{4\pi\nu}e^{-\alpha R^2/4\nu}$.
- **Mirrors:** `core.vorticity.stretching_tilting_split`, `vorticity_terms`, `ch05.uniform_strain_vorticity`,
  `ch05.stretched_tube`, `core.vortices.burgers_vortex`, `core.vorticity.frozen_in_check`.
- **Derivations:** D09 ★★★ (12 steps: the curl steps remove the pressure and gravity bars; the Lamb step draws ω × u;
  the (B.3.10) step splits advection from stretching; the final regrouping lights the three bars of (5.13)), D16 (the
  natural-frame steps draw e_s, e_n, e_m on the line), D17 (the projection steps split the arrow into purple and blue;
  the 2-D step sets the planar preset), D18 (the tube steps shrink the area as the length doubles; the Burgers steps
  set the balance preset).
- **Depth features:** Explain tab (G and ω now → e_s → stretching rate e_s·G e_s and tilting vector with numbers →
  ∣ω∣ growth factor since t = 0 → diffusion term (Burgers) → balance and core radius √(4ν/α) → reading the current
  setting), synced Code tab, + 3-D view, linked views (3), transport, term bars, presets (four flows), status
  ("⬆ stretching: ∣ω∣ grows" / "↪ tilting: ω turns" / "▭ 2-D: neither" / "⚖ Burgers: stretching = diffusion").
- **Follows reference:** `forced_damped_vibrations.html` (transient → steady state on one clock, boxed numbers) with
  `fid_formula_lab.html`'s term bars.
- **Aha:** pull a spinning tube longer and it spins faster, turn it and its spin turns with it — the only torque-free
  ways a flow can change vorticity, and both vanish in a flat, 2-D world.

### E6 · biot_savart_filament
- **A:** C07, C08 (also shows N24 Poisson, N25 (5.14) and its sign, N26–N28 the rewrite, (5.15) and V′, N29 infinite
  line recovered, N52 Fig. 5.8)
- **Confusion removed:** "velocity is caused locally" — the velocity at a point is the sum of contributions from every
  piece of vorticity, each falling off as 1/distance², perpendicular to both the vortex element and the line joining
  them (and the book's printed (5.14) sign would reverse the swirl).
- **Why interactive:** drag a field point around a filament (straight, bent, square loop, ring, helix) and see each
  segment's contribution arrow and their vector sum; lengthen a straight segment and watch Γ/4πd(cos θ_a − cos θ_b)
  approach Γ/2πd; the "sign" toggle shows the printed −1/(4π) turning the flow backwards.
- **Stage:** (1) the filament (3-D, orbitable, or planar) with the field point, per-segment arrows (fading with
  distance) and the total arrow; (2) contribution per unit length along the filament (a curve whose area is the
  answer) — the unrolled integral; (3) (hidePortrait) speed vs distance from a straight segment with the infinite-line
  ghost Γ/2πd.
- **Controls:** shape: straight segment / semi-infinite / square loop / ring / helix / Gaussian tube (volume form) · Γ ·
  segment half-length or ring radius · number of segments M · field point (drag) · sign: +1/(4π) / printed −1/(4π).
- **Equations shown (live):** $\nabla\times\boldsymbol\omega=-\nabla^2\mathbf u$; (5.14) with the ⚠️ sign note; (5.15)
  $\int_{V'}\nabla'\times(\frac{\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert})d^3x'=\int_{A'}\frac{\mathbf n\times\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}d^2x'$;
  (5.16) $\mathbf u=\frac1{4\pi}\int_{V'}\frac{\boldsymbol\omega\times(\mathbf x-\mathbf x')}{\lvert\mathbf x-\mathbf x'\rvert^3}d^3x'$;
  (5.17) $d\mathbf u=\frac{\Gamma\,dl}{4\pi}\mathbf e_\omega\times\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$;
  segment $\frac{\Gamma}{4\pi d}(\cos\theta_a-\cos\theta_b)$; ring axis $\frac{\Gamma R^2}{2(R^2+z^2)^{3/2}}$.
- **Mirrors:** `core.biot_savart.segment_induced_velocity`, `filament_velocity`, `ring_axis_velocity`,
  `biot_savart_volume`, `velocity_from_curl_omega`, `core.biot_savart.filament_preset` (§8).
- **Derivations:** D10 ★★★ (12 steps: the Poisson step shows ∇²u next to −∇×ω; the Green's-function steps draw the
  1/r kernel; the sign step flips the arrow with the toggle), D11 ★★★ (12 steps: the product-rule steps split the
  integrand; the Gauss step draws V′ with its ends and side; the reorder step turns the kernel into ω × r), D12 (the
  far-field step shrinks the core), D13 (the substitution step draws the angles θ_a, θ_b; the limit step stretches the
  segment to infinity).
- **Depth features:** Explain tab (distance d and angles → each segment's contribution with numbers → the sum → the
  closed form → the infinite-line comparison → reading the current setting), synced Code tab, + linked views (3), 3-D
  view, inspector (click a segment: its dl, r, (e_ω × r)/r³ and du), presets (six shapes), status ("✅ matches the
  closed form to …" / "⚠️ printed sign: swirl reversed").
- **Follows reference:** `forward_noising_lab.html` (click a point to see its exact arithmetic, a sum shown as parts +
  total) with `forced_damped_vibrations.html`'s explanation panel.
- **Aha:** every bit of vorticity stirs the whole fluid, like a current making a magnetic field; add up the little
  1/r² pushes and an infinitely long vortex gives back exactly Γ/2πr.

### E7 · vorticity_equation_rotating
- **A:** C09, C11 (also shows R12–R17 the index route, N30 (5.18), N32 (5.25), N35 (5.28), N37 term reading, R18
  absolute vorticity, N39 planetary stretching, N54 Fig. 5.10; links to ch04 `rotating_frame_coriolis`)
- **Confusion removed:** "vorticity in the atmosphere is made by winds" — on a rotating planet a column that is merely
  stretched or squashed borrows or returns spin from the planet's 2Ω, and density contrasts make more; the relative
  vorticity is not conserved, the absolute circulation (and (ζ + 2Ω)/h) is.
- **Why interactive:** in the column mode the reader drags the layer depth (a column crossing a slope or a ridge) and
  the latitude and watches the column spin up or down with (ζ + 2Ω)/h pinned; in the ring mode a ring of air is moved
  poleward and acquires anticyclonic circulation; in the budget mode the term bars of (5.30) are shown at a probe of
  three scenes — the "same equation, three regimes" view needs the controls.
- **Stage:** (1) the rotating layer (side view) with a column over a sloping bottom, its height h and spin shown by a
  turning marker, or a map view with a ring of air on a globe cap; (2) ζ vs h (or vs latitude) with the conserved
  ratio's ghost line and the current point; (3) (hidePortrait) term bars of (5.30): local, advective, relative
  stretching/tilting, planetary (amber), baroclinic (orange), diffusion (rose), residual.
- **Controls:** mode: column / ring of air / budget · column height h (or bottom slope) · latitude φ · initial ζ₀ ·
  scene (budget mode: Burgers / lock exchange / rotating stretched column).
- **Equations shown (live):** (5.20); (5.25)
  $\partial u_i/\partial t+(\tfrac12u_j^2+\Phi)_{,i}-\varepsilon_{ijk}u_j(\omega_k+2\Omega_k)=-(1/\rho)p_{,i}-\nu\varepsilon_{ijk}\omega_{k,j}$;
  (5.28) $\frac1{\rho^2}\nabla\rho\times\nabla p$; (5.30)
  $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u+\frac1{\rho^2}\nabla\rho\times\nabla p+\nu\nabla^2\boldsymbol\omega$;
  $D\omega_z/Dt=2\Omega\,\partial w/\partial z$; (5.33) $\frac{D\Gamma_a}{Dt}=0$,
  $\Gamma_a=\Gamma+2\int_A\boldsymbol\Omega\cdot\mathbf n\,dA$; column $(\omega_z+2\Omega)/h$ = const.
- **Mirrors:** `core.vorticity.vorticity_budget`, `planetary_vorticity_terms`, `absolute_circulation`,
  `ch05.column_relative_vorticity`, `ch05.relative_circulation_after_move`, `ch05.column_over_slope` (§8).
- **Derivations:** D14 (the Lamb-form steps collect the ω + 2Ω bracket), D15 ★★★ (15 steps: each ε–δ step lights the
  bar it produces — stretching of absolute vorticity (purple + amber), baroclinic (orange), diffusion (rose); the Ω = 0
  check step sets the Burgers scene), D19 (the Coriolis-loop step draws the loop's vector area changing; ring mode),
  D20 (column mode: the mass step fixes A h, the division step pins the ratio line).
- **Depth features:** Explain tab (f = 2Ω sin φ with numbers → column h and h₀ → ζ from the conserved ratio → sense
  (cyclonic/anticyclonic) → ring: Γ₁ = Γ₀ + 2Ω(A₀ sin φ₀ − A₁ sin φ₁) → mean ζ → budget terms at the probe → reading
  the current setting), synced Code tab, + linked views (3), modes (three), presets (column over a ridge, column into a
  trough, ring 30° N → 60° N, Southern Hemisphere), term bars, status ("🌀 cyclonic spin-up" / "🔃 anticyclonic" /
  "⚖ ratio conserved to …"), a small table of real values (f at 0°, 30°, 45°, 60°, 90° with the current row lit).
- **Follows reference:** `angular_frequency_explorer_1.html` (modes, "Right now" notes with a highlighted table) with
  `fid_formula_lab.html`'s term bars.
- **Aha:** squash a column on a spinning planet and it must spin the other way to keep (ζ + f)/h — the seed of
  potential vorticity, cyclones over plains and anticyclones over mountains.

### E8 · point_vortex_lab
- **A:** C12, C13 (also shows N40 same-sign pair, N41 opposite pair, N55, N56, N57 knife-blade pair in a bucket, N58
  wall image)
- **Confusion removed:** "a vortex moves by itself" — a straight vortex never moves itself; it goes wherever the others
  (and the walls, through their images) carry it: two like vortices orbit, two opposite ones march off together, a
  vortex near a wall slides along it.
- **Why interactive:** click to drop vortices of either sign, drag them, press play: the orbit rate, the centre of
  vorticity (fixed), the pair speed and the wall drift appear as motion; toggling the wall adds the mirror images and
  the reader sees them move in lock-step — no static figure makes "the image is doing the pushing" obvious.
- **Stage:** (1) the plane with vortices (teal +, rose −), their velocity arrows, trails, the centre of vorticity G,
  optional wall or circle boundary with ghost images; (2) invariants over time (ΣΓx, ΣΓ∣x∣², H) — flat lines that prove
  the integration; (3) (hidePortrait) distance and angle of the pair vs t with the predicted rate ghost.
- **Controls:** preset: equal pair / Γ₂ = 3Γ₁ / opposite pair / vortex near a wall / knife pair in a bucket / three
  vortices · Γ₂/Γ₁ · separation h (drag) · boundary: none / wall / circle · time (transport).
- **Equations shown (live):** $\frac{d\mathbf x_k}{dt}=\sum_{j\ne k}\frac{\Gamma_j}{2\pi}\frac{\mathbf e_z\times(\mathbf x_k-\mathbf x_j)}{\lvert\mathbf x_k-\mathbf x_j\rvert^2}$;
  $V_1=\Gamma_1/2\pi h$, $V_2=\Gamma_2/2\pi h$; $h_1=\Gamma_2h/(\Gamma_1+\Gamma_2)$; rate $(\Gamma_1+\Gamma_2)/2\pi h^2$;
  opposite pair $\Gamma/2\pi h$; wall drift $V_A=\Gamma/4\pi h$.
- **Mirrors:** `core.biot_savart.point_vortex_velocity`, `point_vortex_evolve`, `point_vortex_invariants`,
  `wall_image_system`, `circle_image_system`, `ch05.vortex_pair`, `ch05.vortex_near_wall_speed`,
  `core.biot_savart.point_vortex_preset` (§8).
- **Derivations:** D21 (the lever-rule step marks G; the orbit-rate step draws V₁/h₂; the Γ₁ + Γ₂ = 0 step sets the
  opposite pair), D22 (wall preset: the mirror step draws the image; the symmetric wall points show normal components
  cancelling).
- **Depth features:** Explain tab (each vortex's induced velocity with numbers → G → orbit period or pair speed →
  image velocity and drift → invariants → reading the current setting), synced Code tab, + linked views (3), transport
  (end-of-run card "G never moved"), presets (six), inspector (click a vortex: the sum of the others' contributions),
  status ("🔄 co-rotating about G" / "➡ translating pair" / "🧱 sliding along the wall").
- **Follows reference:** `random_copy_lab.html` (click-to-place interaction, code comments showing the live values)
  with `angular_frequency_explorer_1.html`'s end-of-run summary.
- **Aha:** vortices are pushed only by each other — like dancers holding hands they orbit when alike and walk away
  together when opposite, and a wall is just an invisible partner behind the mirror.

### E9 · vortex_sheet_rollup
- **A:** C14 (also shows N44 sheet definition, N60 Fig. 5.16; C12's point-vortex dynamics; forward link to Ch. 11
  Kelvin–Helmholtz)
- **Confusion removed:** "a velocity jump is just a boundary" — it is a sheet of vorticity whose strength per length
  equals the jump; and a sheet is not static: a small ripple makes its own vortices bunch and roll up.
- **Why interactive:** the reader sets N filaments and watches u(y) across the row converge to a clean jump ±γ/2 (and
  the circuit's circulation equal the jump × ds); then perturbs the sheet and presses play to watch it roll up into
  cat's-eye vortices — the link between the jump, the vortex row and the instability needs motion.
- **Stage:** (1) the sheet as a row of point vortices (colour = strength) with the dn × ds circuit and, in roll-up mode,
  the evolving sheet with trails; (2) u(y) across the sheet with the ideal jump ghost; (3) (hidePortrait) circulation of
  the circuit vs its height dn (flat once dn covers the sheet), or the jump error vs N (log–log, slope −1).
- **Controls:** mode: jump / roll-up · N filaments (log) · γ (m/s) · ripple amplitude · smoothing δ (roll-up) · time
  (transport) · convention: counterclockwise (text) / clockwise (caption).
- **Equations shown (live):** $[u_t]\ne0$, $[v_n]=0$; $d\Gamma=u_2\,ds+v\,dn-u_1\,ds-v\,dn=(u_2-u_1)\,ds$; strength
  $\gamma=u_2-u_1$; far field $u=\mp\gamma/2$.
- **Mirrors:** `core.biot_savart.vortex_sheet_velocity`, `ch05.vortex_sheet_strength`, `ch05.discrete_sheet_convergence`,
  `ch05.sheet_rollup` (§8).
- **Derivations:** D23 (the four-sides steps light each side of the circuit; the dn → 0 step shrinks it; the
  convention step toggles the caption's sign).
- **Depth features:** Explain tab (γ and N → each filament's strength γ ds → u just above and below with numbers → the
  jump → the circuit's circulation → the convention → reading the current setting), synced Code tab, + linked views
  (3), transport (roll-up), presets (fine sheet N = 1000, coarse N = 10, rippled sheet, caption convention), status
  ("✅ jump = γ to …%" / "🌀 rolling up"), inspector (click a height: the sum of filament contributions).
- **Follows reference:** `overfitting_curves.html` (a minimal two-slider figure with a verdict line) with
  `forced_damped_vibrations.html`'s explanation panel.
- **Aha:** a jump in velocity *is* a sheet of vorticity — its strength is the jump — and a sheet left to itself
  wrinkles and rolls up into vortices.

### B1 · vortex_rings (backup)
- **A:** C13, C08 (also shows N42 ring near a wall, N43 leap-frogging, N59 Fig. 5.15)
- **Confusion removed:** "a smoke ring just coasts" — it moves by its own induced velocity (Γ/4πR[ln(8R/a) − ¼]), is
  widened and slowed by a wall's image ring, and two rings take turns passing through each other.
- **Why interactive:** axisymmetric rings on one clock; the reader launches one ring at a wall or two coaxial rings and
  watches R(t), z(t); moving the second ring's start distance changes the leap-frog period.
- **Stage:** meridional plane (R, z) with ring cross-sections, their images and trails; R(t), z(t) curves; impulse
  ΣΓπR² flat. **Controls:** mode (wall / leap-frog), Γ, R₀, core a, separation, time. **Equations:** (5.17) applied to a
  circle, $U=\frac{\Gamma}{4\pi R}[\ln(8R/a)-\tfrac14]$, ring axis $\frac{\Gamma R^2}{2(R^2+z^2)^{3/2}}$. **Mirrors:**
  `core.biot_savart.ring_ring_velocity`, `ring_self_velocity`, `ch05.ring_dynamics`.
- **Derivations:** none (D12 is linked). **Depth features:** explain, code, + linked views, transport, presets, status.
  **Follows reference:** `angular_frequency_explorer_1.html` (linked views on one clock, end-of-run summary). **Aha:**
  a ring is carried by its own curvature and pushed around by its neighbours' images — near a wall it spreads out and
  slows, next to a twin it plays leap-frog forever.

## 6. Python animations and interactive figures (A ID → what, why, player/figure kind)
Animations (`animate` + `show_animation`, ≤ 120 frames, dpi 80, FAST halves the frames):
- **C03** — a material loop advected by a steady cellular flow and by the Lamb–Oseen vortex, with Γ(t) drawn beside it
  (flat vs decaying); why: "circulation stays while the loop is mangled" must be watched; `player="video"`.
- **C10** — a vortex tube stretched by axial strain (length ×, radius shrinking, spin marker speeding up) and a
  segment tilted by shear; why: stretching is a process; `player="frames"` (stop at ×2, ×4).
- **C12** — same-sign pair orbiting G and an opposite pair translating, side by side; why: the motion is the result;
  `player="video"`.
- **C13** — the knife-blade pair launched inside a circular bucket (circle images): crossing, then separating along the
  wall; and two coaxial rings leap-frogging (N43); why: the image's effect is only visible in time; `player="video"`.
- **C14** — a rippled vortex sheet (N point vortices with δ-smoothing) rolling up into cat's eyes; why: the sheet's
  own induced velocity is the story, and it foreshadows Ch. 11; `player="video"`.

Interactive figures (`slider_figure`, plotly, precomputed; 3-D plotly where noted):
- **C01** — plotly 3-D: stream tube and vortex tube traced in a Burgers vortex (N45); why: tubes are 3-D objects.
- **C02** — slider over the tank rate: isobars and free surface of (5.6) beside the funnel of (5.7) (N46, N47); why:
  bowl vs funnel as the rate grows.
- **C04** — slider over the isopycnal tilt: disc with net-force line and torque, spin-up rate vs angle (N50); why: the
  torque switches on as the isolines cross.
- **C06** — slider over time: ω_z(y, t) and u(y, t) of the diffusing vortex sheet (N22) with ∫ω dy = γ printed; why:
  diffusion of vorticity in time.
- **C07** — slider over the number of quadrature nodes (or grid spacing): Biot–Savart u_θ(r) of a Gaussian vortex
  converging to the exact curve, with the printed-sign curve reversed; why: the law and its sign at a glance.
- **C08** — slider over the number of segments M: polygon ring's axis velocity converging to ΓR²/2(R² + z²)^{3/2}; why:
  (5.17) summed.
- **C10** — plotly 3-D: a helical vortex line with its Frenet frame and the stretching/tilting arrows of a chosen G
  (N53); why: a 3-D statement.
- **C11** — slider over latitude: ζ vs h of a column (ratio line) and the ring of air's acquired circulation (N54);
  why: the planetary effect grows with sin φ.
- **C14** — slider over N: u(y) across a discrete sheet converging to the jump (N60); why: from filaments to a sheet.

## 7. From-scratch moments
Every book section with a computable A item has at least one:
- §5.1 **C01** — midpoint sums of ω·n over two cross-sections of a tube vs `core.vorticity.tube_flux_budget` (equal
  fluxes); **C02** — (5.5a) integrated outward by `cumulative_trapezoid` vs `ch05.solid_body_pressure` and
  `line_vortex_pressure` (`assert np.allclose`).
- §5.2 **C03** — loop points advected by a hand-written RK4 and ∮u·dx by the periodic trapezoid rule vs
  `core.vorticity.material_circulation`; **C04** — the pressure force and torque on a disc summed over N rim points vs
  `core.vorticity.pressure_torque_on_element` and vs ∇ρ × ∇p/ρ².
- §5.3 **C05** — the vorticity flux through a small loop on a tube wall computed by hand at t = 0 and after advection
  (stays ≈ 0).
- §5.4 **C06** — (ω·∇)u and ν∇²ω by central differences on the Burgers field vs `core.vorticity.vorticity_terms`.
- §5.5 **C07** — a direct double loop over grid cells of the 2-D kernel ẑ × r/2πr² vs `core.biot_savart.biot_savart_2d`;
  **C08** — a midpoint sum of (5.17) along a straight segment vs the closed form `segment_induced_velocity`.
- §5.6 **C09** — the baroclinic vector by central differences vs `core.vorticity.baroclinic_term`; **C10** — the
  stretching/tilting projection (e_s·G e_s)ω e_s and its remainder by hand vs `stretching_tilting_split`; **C11** — the
  vector area ½∮x × dx by the trapezoid rule vs `core.vorticity.loop_vector_area`.
- §5.7 **C12** — dx_k/dt by a double loop over pairs vs `core.biot_savart.point_vortex_velocity`; **C13** — the image
  sum for a vortex near a wall vs `wall_image_system` (normal velocity on the wall = 0).
- §5.8 **C14** — the velocity of a row of N point vortices summed by hand vs `core.biot_savart.vortex_sheet_velocity`.

## 8. Notes for the implementer
Functions the A items' figures and the explainers need that analysis §4 did not plan (all scalar-callable for parity
rows, SI units, docstrings citing § and Eq.):
- `ch05.sheet_rollup(N, gamma, amplitude, delta, t_eval, L=1.0)` — periodic vortex-sheet roll-up with Krasny
  δ-smoothing (point vortices, `solve_ivp` DOP853; periodic kernel ½cot sums or a large finite row) → positions (T, 2,
  N); for E9 and the C14 animation (label `qualitative` for the roll-up shape; invariants conserved → `conserved`).
- `ch05.vortex_pressure_scenario(kind, r, z, **p) -> dict(u_theta, p, B, sigma_rtheta, net_viscous_force)` — E2's four
  vortices (solid body, line, Rankine, rotating cylinder) in one dispatcher.
- `core.vorticity.kelvin_scenario(name, **p) -> dict(u_fn, p_fn, rho_fn, nu, loop0)` — E3's five flows (cellular,
  Rankine straddling loop, Lamb–Oseen, baroclinic lock-exchange field, rotating frame) and the C03 animation.
- `ch05.lock_exchange_fields(rho1, rho2, delta, g, H)` — smoothed-step ρ and hydrostatic p on a grid for E4 and N49.
- `core.biot_savart.filament_preset(name, **p) -> polyline (3, M)` — straight, semi-infinite, square, ring, helix for E6.
- `ch05.column_over_slope(x, h_fn, lat, zeta0, Omega=OMEGA_EARTH)` — ζ along a path over a varying depth (E7, N54).
- `core.biot_savart.point_vortex_preset(name) -> (xv, Gamma, boundary)` — E8's six configurations.
- `ch05.helical_vortex_line(a, c, s)` and `ch05.frenet_frame(curve, s)` — N53 figure and E5's 3-D line.
- `core.vortices.twisted_gaussian_tube_field(...)` (or reuse `burgers_vortex_field`) — E1's second tube; plus a
  deliberately non-solenoidal "broken" field for the failing preset (flagged in its docstring as not a vorticity field).
- Parity rows must include: tube flux lower = upper (E1), rim height ω²R²/8g and funnel deficit ρΓ²/8π²r² (E2), Γ(t)
  constant for the cellular flow and Lamb–Oseen Γ(r, t) (E3), 2·torque/I_G vs ∇ρ × ∇p/ρ² and the lock-exchange rate
  (E4), e^{αt} growth and the Burgers core √(4ν/α) (E5), segment law and ring axis ΓR²/2(R² + z²)^{3/2} (E6), column
  ratio and Γ₁ = Γ₀ + 2Ω(A₀ sin φ₀ − A₁ sin φ₁) (E7), pair rate (Γ₁ + Γ₂)/2πh² and wall drift Γ/4πh (E8), sheet jump
  ±γ/2 (E9). Book numbers stay in `tests/book_values_ch05.json`.
- Budget: fields on ≤ 128×128 grids; the 3-D Biot–Savart quadrature ≤ 24³ nodes with FAST (40³ otherwise) and cached;
  point-vortex and roll-up runs precomputed once (`functools.lru_cache`), roll-up N = 100 with FAST (200 otherwise);
  sympy derivation checks (D09, D10, D11, D15) cached; animations use FAST = 30 frames; the whole notebook < 5 min on
  Colab CPU.
- Notation for the notebook's notation cell: ω = vorticity (never the rotation rate in ch05 code), `Gamma` [m²/s]
  circulation vs `gamma` [m/s] sheet strength, σ viscous stress (ch03's core radius is `sigma_core`), ε Levi-Civita vs
  `eps_diss`, e_n in natural coordinates points away from the centre of curvature (book's Fig. 5.9 convention).
