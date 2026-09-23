# Chapter 3 — Kinematics: curation
(from `analysis/ch03.md` — 79 inventory rows, 35 numbered equations (3.1)–(3.35), 33 derivations in §2b (30 + 3
recap rows), 53 implementation rows; `book.yaml` ch03, `policy.tier_a_full_treatment: [12, 18]`, `coverage: exhaustive`.
Curated 2026-09-22, concept-curator. Read: `knowledge/CUMULATIVE.md`, `concept_map.md`, `primers.md` (P01–P86),
`notation.md`, `viz_patterns.md`, `knowledge/ch02.md` §8, `analysis/ch02_curation.md` (format).)

Counts: A 15 · B 52 · C 12 · RECAP 7 · SKIP 2 · derivations written out 24 (★ 6 ★★ 17 ★★★ 1) · demoted to statements 9

Tier words (parsed by `tools/nbkit.py`): **CORE 15 · NOTE 55 · RECAP 7 · SKIP 2** = 79 rows. Depth: A 15 (all CORE) ·
B 52 (46 NOTE + R01, R02, R03, R05, R06, R07) · C 12 (9 NOTE + R04 + S01 + S02). Every inventory row number `[#n]`
appears exactly once in §2.

**Decisions that shape this chapter.**
1. **SEEN rows are never A.** The analyst marked 15 rows SEEN (ch02 taught G, S + ½R, R ↔ ω, the curl, principal axes,
   Stokes/circulation, plane polar components, the test fields b × x, (Γx₂, 0), K/r). Pure restatements become RECAPs
   (R01–R07, each reusing the ch02 function); "SEEN / NEW" rows whose new part is the *fluid reading* (#5 coordinate
   systems, #45 irrotational flow, #49 principal frame, #53–#55 shear flow and solid body, #58–#59 line vortex) become B
   NOTEs inside the nearest A block. Consequence: §3.1 has no A item (its load-bearing content #5 is SEEN/NEW) and is
   covered by B/C items tagged C01 in its own short notebook section (ch02 decision 6 pattern); §3.5's shear flow is a
   B pair (N34 → C10, N35 → C12) opening the §3.5 section.
2. **Fifteen A items** (analyst's 18 candidates merged): streamlines and path lines stay separate A blocks (both are
   ODE tools that later chapters call), but Ex. 3.1 (#22) and streak lines (#21) are B inside the path-line block C04;
   Leibniz's theorem (#67) is B inside the RTT block C15 with its proof D21 kept under rule (b) (the book never proves
   it) — the 1-D picture opens the RTT block; irrotational flow and circulation are B/RECAP inside the spin block C10
   and reused in C13; (3.20)–(3.21) are B inside the ellipse block C12.
3. **Derivations the book leaves to exercises are written out** when the lesson needs them (rule (b)): (3.7) Ex. 3.3,
   (3.9) Ex. 3.12, S = 0 for rigid motion Ex. 3.17, (3.14) Ex. 3.18, ω′ = ω − 2Ω Ex. 3.19, (3.23) Ex. 3.20, the
   Gaussian maximum Ex. 3.26, RTT → (3.14) Ex. 3.28, RTT → (3.5) Ex. 3.30, the Ex. 3.1 eliminations, the θ-independence
   of the shear-flow spin, the Gaussian u_θ from ω_z, Leibniz's rule. The exercise *text* is never published; the
   derivations are in our words.
4. **Book typos handled in our own words** (analysis §9): (3.6) is taught as $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\lvert\mathbf u\rvert\frac{\partial F}{\partial s}$
   with a ⚠️ note that the printed second term lacks its F; Ex. 3.2's "b = 0 on the base" is taught as "b·n = 0 on the
   base"; the "[?]" in Ex. 3.2's integrand is the factor z; §3.4's pointer to "Section 2.12" means §2.11; (3.19) follows
   from (3.10), (3.11) and (3.15) (not (3.14)); the eigenvalue −γ/2 is "compression at rate γ/2"; Fig. 3.16's C is B.
5. **Conventions stated where first used** (⚠️ callouts): book R = G − Gᵀ (no ½), spin = ½ω (C10); three different
   gammas — ch02 Γ ≡ S₁₂, ch03 γ = du₁/dx₂ = 2S₁₂, and Γ = circulation [m² s⁻¹] (N34, R07); ω₃ = −γ is clockwise
   (N34); t′ (drawing instant) vs t_o (release time) vs ξ_o (amplitude) in Ex. 3.1 (N18); Galilean ≠ rotating frame
   (C05, N28 — the climate-relevant case: relative vs absolute vorticity); ΔV in the RTT is a **signed** volume (D22).
6. **Climate hooks, where real:** thermal advection −u·∇T as the advective part of DT/Dt (C02 worked number: a
   southerly wind across a temperature gradient); relative vs planetary vorticity ω′ = ω − 2Ω (N28 → Ch. 4 §4.7, Ch. 13);
   Gaussian/Rankine vortices as cyclone and tornado profiles (C14); strain axes and frontogenesis (C12 → Ch. 13).

## 1. Teaching order (A IDs grouped by book section, B/C IDs under each; one sentence each: "once you see X, Y follows")
A items in **bold**; B and C items listed where they are taught, with depth in brackets. Equations are written next to
their numbers here so downstream agents have them.

**§3.1 Introduction and Coordinate Systems** (own short notebook section; items tagged C01, C05)
- N01 [C] kinematics = motion without forces (dynamics is Ch. 4) · N02 [B] steady vs unsteady, $\partial(\cdot)/\partial t = 0$
  (revisited in C05: steadiness depends on the observer) · N03 [B] 1-, 2-, 3-D flows and the 1-D section average
  $\bar u(z)=\frac1A\int_A u\,dA$ (Poiseuille → U/2) · N04 [B → C05] the cylinder seen from the body and from the fluid,
  $\mathbf u=\mathbf U+\mathbf u'$ · N05 [B] cylindrical $(R,\varphi,z)$ and spherical $(r,\theta,\varphi)$ coordinates,
  $R=\sqrt{x^2+y^2},\ \varphi=\tan^{-1}(y/x)$, $r=\sqrt{x^2+y^2+z^2},\ \theta=\tan^{-1}(\sqrt{x^2+y^2}/z)$, unit vectors
  and velocity components as projections (plotly 3-D figure) · N06 [C] Appendix B operators → Ch. 4.

**§3.2 Particle and Field Descriptions of Fluid Motion**
- **C01 Lagrangian and Eulerian descriptions and the bridge (3.2)** [#9]: once a particle is a label
  $\mathbf r(t;\mathbf r_o,t_o)$ and a field is $F(\mathbf x,t)$, the compatibility
  $F[\mathbf r(t;\mathbf r_o,t_o),t]=F(\mathbf x,t)$ when $\mathbf x=\mathbf r(t;\mathbf r_o,t_o)$ *(Eq. 3.2)* turns one
  description into the other (D01: $x = Xe^{\alpha t}$ ⇒ $u = \alpha x$). Inside: N07 [B] Lagrangian labels, N08 [B]
  $\mathbf u=d\mathbf r/dt,\ \mathbf a=d^2\mathbf r/dt^2$ *(Eq. 3.1)*.
- **C02 The material derivative (3.4)–(3.5)** [#11]: once F is followed along a particle's path, the chain rule
  $\frac{d}{dt}F[\mathbf r(t),t]=\frac{\partial F}{\partial r_i}\frac{dr_i}{dt}+\frac{\partial F}{\partial t}$ *(Eq. 3.3)* with
  $dr_i/dt=u_i$ gives $\frac{DF}{Dt}\equiv\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$ *(Eq. 3.5)* — the local
  plus advective change that every conservation law of Ch. 4 is written with (D02). Inside: N09 [B] (3.3), N10 [B] local
  and advective parts and when each vanishes (thermal advection hook), N11 [B] vector and index forms (3.5), N12 [B] the
  streamwise form $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\lvert\mathbf u\rvert\frac{\partial F}{\partial s}$
  *(Eq. 3.6)* with the ⚠️ missing-F note.

**§3.3 Flow Lines, Fluid Acceleration, and Galilean Transformation**
- **C03 Streamlines (3.7)** [#16]: once a streamline is "tangent to u at one frozen instant",
  $\mathbf u\times d\mathbf s=0$ gives $dx/u=dy/v=dz/w$ *(Eq. 3.7)* (D03) and a streamline is an ODE in arc length
  integrated at fixed t. Inside: N14 [B] (3.7) and the 2-D slope $dy/dx = v/u$, N15 [B] stream tubes (no flux through
  the wall; equal flux through two sections).
- **C04 Path lines (3.8) and streak lines — Ex. 3.1** [#20]: once a particle obeys
  $d\mathbf r/dt=\mathbf u(\mathbf r,t)$, $\mathbf r(t_o)=\mathbf r_o$ *(Eq. 3.8)*, a streak line is a family of path
  lines labelled by release time, and the three lines differ exactly when the flow is unsteady (D04, D05). Inside: N16
  [B] path line = fixed-identity trajectory, N17 [B] streak line parametrised by $t_o$, N18 [B] Ex. 3.1 (line
  $y=x\tan\omega t'$ + two circles of radius $\xi_o$, tangent at the origin), N19 [B] our Fig. 3.7 animation, N13 [B]
  steady flow ⇒ all three coincide (PIV named).
- **C05 Galilean invariance of the fluid acceleration (3.9)** [#26]: once $\mathbf u=\mathbf U+\mathbf u'$,
  $\mathbf x=\mathbf x'+\mathbf Ut+\mathbf x'_o$ is put through the chain rule, the $\mathbf U\cdot\nabla'$ terms cancel
  and $\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=\frac{\partial\mathbf u'}{\partial t'}+(\mathbf u'\cdot\nabla')\mathbf u'$
  *(Eq. 3.9)* (D06): the sum is frame-free, the split is not. Inside: N21 [B] the Galilean transformation, N20 [B] Fig.
  3.2 revisited (steady in one frame, unsteady in the other), N22 [B] local term linear, advective quadratic (the
  nonlinearity of fluid mechanics), N23 [B] the frame-dependent split as term bars, N24 [C] traffic-light and
  roller-coaster analogies, N04 (from §3.1).

**§3.4 Strain and Rotation Rates**
- **C06 Relative velocity near a point: $du_i=(\partial u_i/\partial x_j)\,dx_j$ (3.10)** [#31]: once the neighbour's
  velocity is Taylor-expanded (D07), everything about a small element's motion lives in the nine numbers of G, and the
  ch02 split $\partial u_i/\partial x_j=S_{ij}+\tfrac12R_{ij}$ *(Eq. 3.11)* tells us which part deforms and which part
  spins. Inside: N25 [C] why deformation rates matter (stress, Ch. 4), R01 [B] (3.11), R02 [B]
  $S_{ij}=\tfrac12(\partial u_i/\partial x_j+\partial u_j/\partial x_i)$ *(Eq. 3.12)*, R03 [B]
  $R_{ij}=\partial u_i/\partial x_j-\partial u_j/\partial x_i$ *(Eq. 3.13)* (⚠️ no ½), N26 [C] S links to stress, R does
  not.
- **C07 Linear strain rate: $\frac{1}{\delta x_1}\frac{D(\delta x_1)}{Dt}=\partial u_1/\partial x_1$** [#36]: once a
  material segment's two ends move at $u_1$ and $u_1+(\partial u_1/\partial x_1)\delta x_1$, its stretching per unit
  length is the diagonal of S — and along any direction n it is $\mathbf n\cdot\mathbf S\cdot\mathbf n$ (D08).
- **C08 Shear strain rate: $S_{12}=\tfrac12 D(\alpha+\beta)/Dt$** [#37]: once two perpendicular segments tilt by
  $d\alpha$ and $d\beta$, the off-diagonal S is half their closing rate, $\mathbf n_1\cdot\mathbf S\cdot\mathbf n_2$ in
  general (D09); and a rigid motion $\mathbf U+\boldsymbol\Omega\times\mathbf x$ has S = 0 (D10). Inside: N27 [B]
  rigid motion ⇒ S = 0, S independent of a translating/rotating frame.
- **C09 Volumetric strain rate $\frac1{\delta V}\frac{D(\delta V)}{Dt}=\partial u_i/\partial x_i=S_{ii}$ (3.14)**
  [#39]: once $\delta V=\delta x_1\delta x_2\delta x_3$ and each side stretches by C07, the product rule gives
  $\nabla\cdot\mathbf u$ as the volume growth rate per volume (D11) — the kinematic half of Ch. 4's continuity equation.
- **C10 Vorticity is twice the element's spin: $\tfrac12 D(-\alpha+\beta)/Dt=\tfrac12\omega_3$** [#43]: once the
  *average* turning rate of two perpendicular lines is taken (D12) and shown to be the same for any pair (D14 in the
  shear flow), R ↔ ω *(Eqs. 3.15–3.16)* says the element spins at $\omega/2$, and that spin depends on the frame (D13:
  $\omega'_z=\omega_z-2\Omega$). Inside: R04 [C] R antisymmetric ↔ a vector, R05 [B]
  $R_{ij}=-\varepsilon_{ijk}\omega_k$ *(Eq. 3.15)*, R06 [B]
  $\omega_3=\partial u_2/\partial x_1-\partial u_1/\partial x_2$ etc. *(Eq. 3.16)*, N28 [B] frame dependence
  (absolute vs relative vorticity → Ch. 13), N29 [B] irrotational flow $\boldsymbol\omega=0$, $\mathbf u=\nabla\phi$
  *(Eq. 3.17)* with the simply-connected caveat, R07 [B] circulation
  $\Gamma\equiv\oint_C\mathbf u\cdot d\mathbf s=\int_A\boldsymbol\omega\cdot\mathbf n\,dA$ *(Eq. 3.18)*, N34 [B] the
  parallel shear flow (§3.5, placed at the head of the §3.5 section).
- **C11 Relative velocity = deformation + rigid rotation (3.19)** [#47]: once $\tfrac12R_{ij}=-\tfrac12\varepsilon_{ijk}\omega_k$
  is substituted in (3.10), $du_i=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i$ *(Eq. 3.19)* (D15): a small
  element deforms by S and turns rigidly at $\omega/2$ — the reason only S can enter the stress law of Ch. 4. Inside:
  N30 [C] rigid-body velocity $\boldsymbol\Omega\times\mathbf x$, N33 [C] the section's summary.
- **C12 Principal strain axes: a small sphere becomes an ellipsoid** [#51]: once S is diagonal in its eigenframe,
  $d\bar u_\alpha=\bar S_{\alpha\alpha}d\bar x_\alpha$ *(Eq. 3.21)* stretches each principal direction in proportion to
  its own length, so a sphere maps to an ellipsoid on those axes (D16). Inside: N31 [B] (3.20), N32 [B] (3.21), N35 [B]
  our Fig. 3.14 — the 45° element stretches without shear, the aligned element shears without stretching.

**§3.5 Kinematics of Simple Plane Flows** (opens with N34 [→ C10] and N35 [→ C12], the parallel shear flow)
- **C13 Vorticity in polar coordinates: solid-body rotation vs the irrotational vortex (3.23)** [#56]: once the
  circulation around a small polar sector is divided by its area (D17),
  $\omega_z=\frac1r\frac{\partial}{\partial r}(ru_\theta)-\frac1r\frac{\partial u_r}{\partial\theta}$ *(Eq. 3.23)* gives
  $2\omega_0$ for $u_\theta=\omega_0r$ *(Eq. 3.22)* and 0 for $u_\theta=B/r$ *(Eq. 3.25)* — the paddle-wheel test.
  Inside: N36 [B] (3.22), N37 [B] $\Gamma=2\pi r^2\omega_0$ *(Eq. 3.24)*, N38 [B] (3.25), N39 [B]
  $\Gamma=2\pi B$ *(Eq. 3.26)*, N40 [B] the δ-function core $\lim 2B/r^2$ *(Eq. 3.27)*, N41 [B] $\Gamma_{ABCD}=0$ (elements
  deform but do not spin), R07 reused.
- **C14 Vortices with a core: Rankine (3.28) and Gaussian (3.29)** [#63]: once a uniform-vorticity core of radius σ is
  joined to an irrotational outside (D18), the Rankine profile peaks at σ; smoothing the join gives the Gaussian vortex
  whose $u_\theta$ follows from $\omega_z$ by Stokes (D19) and peaks at $r\approx1.1209\sigma$ (D20). Inside: N42 [C]
  real vortices (bathtub, tornado, cyclone), N43 [B] (3.29) + Lamb–Oseen pointer ($\sigma^2=4\nu t$, Ch. 5/8), N44 [B]
  the maximum $1+2r^2/\sigma^2=\exp(r^2/\sigma^2)$.

**§3.6 Reynolds Transport Theorem**
- **C15 The Reynolds transport theorem (3.35)** [#74]: once the 1-D Leibniz rule (N46, D21) is seen as "interior change +
  what the moving ends sweep in", the 3-D version follows from the definition (3.31) by splitting off the swept volume
  (D22 ★★★): $\frac{d}{dt}\int_{V^*}F\,dV=\int_{V^*}\frac{\partial F}{\partial t}dV+\int_{A^*}F\,\mathbf b\cdot\mathbf n\,dA$
  *(Eq. 3.35)* — and with $\mathbf b=\mathbf u$ it closes the loop to (3.14) (D23) and (3.5) (D24). Inside: N45 [C]
  motivation, N46 [B] Leibniz (3.30), N47 [B] our Fig. 3.17, N48 [B] control volume V*, A*, b, n, N49 [B] (3.31), N50
  [B] (3.32), N51 [B] (3.33), N52 [B] (3.34), N53 [B] interpretations, N54 [B] Ex. 3.2 cone (worked number), N55 [B] our
  Fig. 3.18.

**Exercises and literature** — S01, S02 (one pointer line each at the end of the notebook).

## 2. Chapter map (depth) — every inventory row exactly once
One row per inventory row (79, `[#n]` = analysis §2 row). `A parent` names the A block a B or C item is written in.

| ID | Item | § | Depth | Tier | A parent | Reason (A) / treatment (B) / pointer (C, SKIP) / source chapter (RECAP) |
|---|---|---|---|---|---|---|
| N01 | Kinematics = motion without the forces that cause it [#1] | 3.1 | C | NOTE | C01 | named in the chapter's opening paragraph; dynamics is Ch. 4 §4.4 onward |
| N02 | Steady vs unsteady flow, $\partial(\cdot)/\partial t=0$ [#2] | 3.1 | B | NOTE | C01 | stated in the §3.1 section; the frame-dependence twist is shown in C05 (cylinder steady in the body frame, unsteady in the fluid frame); `is_steady` probe in tests |
| N03 | 1-D, 2-D, 3-D flows; 1-D approximation by section averaging $\bar u(z)=\frac1A\int_A u\,dA$ (Fig. 3.1) [#3] | 3.1 | B | NOTE | C01 | stated with our Fig. 3.1c/d; one number: Poiseuille $u=U(1-r^2/R^2)$ averages to $U/2$ by `ch03.cross_section_average` |
| N04 | Fig. 3.2: steady flow past a fixed cylinder vs the cylinder moving through still fluid, $\mathbf u=\mathbf U+\mathbf u'$ [#4] | 3.1 | B | NOTE | C05 | stated with our two streamline plots (`ch03.cylinder_flow`, frame = body / fluid); the flow itself is Ch. 6 §6.3; used again as C05's stage |
| N05 | Coordinate systems (Fig. 3.3): plane polar, cylindrical $(R,\varphi,z)$, spherical $(r,\theta,\varphi)$; component names; unit vectors [#5] | 3.1 | B | NOTE | C01 | stated with the conversion formulas and a plotly 3-D figure of the local unit vectors at a point (`core.coords`); ⚠️ θ from +z, φ azimuth; plane polar recalls ch02 Ex. 2.1; used by C13, C14, N54 and Ch. 4 |
| N06 | Appendix B holds ∇, ∇², ∇·u, (u·∇)u in curvilinear coordinates [#6] | 3.1 | C | NOTE | C01 | named; used by (3.23) in C13 and by Ch. 4 (Navier–Stokes in cylindrical/spherical coordinates) |
| N07 | Lagrangian description: particle labelled by $\mathbf r_o$ at $t_o$; trajectory $\mathbf r=\mathbf r(t;\mathbf r_o,t_o)$ (Fig. 3.4) [#7] | 3.2 | B | NOTE | C01 | stated in plain words with a drifting-float picture; labels are not variables |
| N08 | Eq. (3.1): $\mathbf u=d\mathbf r(t;\mathbf r_o,t_o)/dt$, $\mathbf a=d^2\mathbf r(t;\mathbf r_o,t_o)/dt^2$ [#8] | 3.2 | B | NOTE | C01 | stated; number with $x=Xe^{\alpha t}$ (X = 2 m, α = 0.5 s⁻¹, t = 1 s: u = 1.649 m/s, a = 0.824 m/s²) by `core.K.lagrangian_velocity_acceleration` |
| C01 | Lagrangian and Eulerian descriptions; compatibility $F[\mathbf r(t;\mathbf r_o,t_o),t]=F(\mathbf x,t)$ when $\mathbf x=\mathbf r(t;\mathbf r_o,t_o)$, Eq. (3.2) [#9] | 3.2 | A | CORE | – | load-bearing: every later chapter switches between following a particle and watching a point; the bridge is what the material derivative, path lines and the RTT are built on |
| N09 | Eq. (3.3): total time derivative along a trajectory by the chain rule $\frac{d}{dt}F[\mathbf r,t]=\frac{\partial F}{\partial r_1}\frac{dr_1}{dt}+\frac{\partial F}{\partial r_2}\frac{dr_2}{dt}+\frac{\partial F}{\partial r_3}\frac{dr_3}{dt}+\frac{\partial F}{\partial t}$ [#10] | 3.2 | B | NOTE | C02 | stated as step 3 of D02; sympy check: differentiate $F(\mathbf r(t),t)$ directly = chain-rule sum |
| C02 | Material derivative $\frac{d}{dt}F[\mathbf r,t]=(\nabla F)\cdot\mathbf u+\frac{\partial F}{\partial t}\equiv\frac{DF}{Dt}$, Eq. (3.4) [#11] | 3.2 | A | CORE | – | load-bearing: every conservation law of Ch. 4 and the vorticity (Ch. 5), wave (Ch. 7) and GFD (Ch. 13) equations are written with D/Dt |
| N10 | Local part $\partial F/\partial t$ and advective part $\mathbf u\cdot\nabla F$; each can vanish; "advection" vs "convection" [#12] | 3.2 | B | NOTE | C02 | stated with term bars from `core.K.material_derivative_terms`; climate number: southerly wind 10 m/s across a gradient of −1 K per 100 km gives $\mathbf u\cdot\nabla T=-10^{-4}$ K/s, so a station sees +0.36 K/h of warm advection when DT/Dt = 0 |
| N11 | Eq. (3.5): $\frac{DF}{Dt}\equiv\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$, or $\frac{DF}{Dt}\equiv\frac{\partial F}{\partial t}+u_i\frac{\partial F}{\partial x_i}$ [#13] | 3.2 | B | NOTE | C02 | stated as D02's result in both notations; `expand_indices("u_i dF/dx_i")` prints the three terms |
| N12 | Eq. (3.6): streamwise form $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\lvert\mathbf u\rvert\frac{\partial F}{\partial s}$ (book prints the second term without F) [#14] | 3.2 | B | NOTE | C02 | stated with $\partial F/\partial s=\mathbf e_u\cdot\nabla F$ (ch02 P75 directional derivative); ⚠️ typo callout; undefined at stagnation points; `core.K.streamwise_derivative` equals the advective term |
| N13 | The three flow lines coincide in steady flow; visualisation pointers (Van Dyke, Samimy et al.) [#15] | 3.3 | B | NOTE | C04 | stated and shown: E1's "steady" preset and a notebook overlay where streamline, path line and streak line through one point agree to 1e-7 |
| C03 | Streamline: a curve instantaneously tangent to $\mathbf u$, $\mathbf u\times d\mathbf s=0$ (Fig. 3.5) [#16] | 3.3 | A | CORE | – | load-bearing: Bernoulli along a streamline (Ch. 4), ψ contours (Ch. 6) and every flow picture in the book; the arc-length ODE `core.K.streamline` is reused everywhere |
| N14 | Eq. (3.7): $dx/u=dy/v=dz/w$ (derivation left to Exercise 3.3) [#17] | 3.3 | B | NOTE | C03 | the result of C03's D03; 2-D slope $dy/dx=v/u$ (number: u = (1, 2) at a point ⇒ slope 2) via `ch03.streamline_slope` |
| N15 | Stream tube: streamlines through a closed curve; no fluid crosses its wall (Fig. 3.6) [#18] | 3.3 | B | NOTE | C03 | stated with a plotly 3-D stream tube; one number: equal flux through two cross-sections (`core.V.flux_through_faces`) |
| N16 | Path line: trajectory of a particle of fixed identity, $\mathbf x=\mathbf r(t;\mathbf r_o,t_o)$ [#19] | 3.3 | B | NOTE | C04 | stated as the Lagrangian trajectory of C01 drawn in space |
| C04 | Path-line ODE $d\mathbf r/dt=[\mathbf u(\mathbf x,t)]_{\mathbf x=\mathbf r}=\mathbf u(\mathbf r,t)$, $\mathbf r(t_o)=\mathbf r_o$, Eq. (3.8) [#20] | 3.3 | A | CORE | – | load-bearing: the Eulerian-to-Lagrangian tool used for particle orbits (Ch. 7), vortex dynamics (Ch. 5), Lagrangian statistics (Ch. 12), trajectories and drift (Ch. 13); streak lines and Ex. 3.1 are built on it |
| N17 | Streak line: all particles that passed through a fixed point $\mathbf x_o$, parametric in the release time $t_o$ [#21] | 3.3 | B | NOTE | C04 | stated as "one path line per release time, positions at the drawing time"; `core.K.streakline`; dye-injection picture; D05 is its worked case |
| N18 | Ex. 3.1: $u=\omega\xi_o\cos\omega t$, $v=\omega\xi_o\sin\omega t$; streamline $y=x\tan(\omega t')$, path line $(x+\xi_o\sin\omega t')^2+(y-\xi_o\cos\omega t')^2=\xi_o^2$, streak line $(x-\xi_o\sin\omega t')^2+(y+\xi_o\cos\omega t')^2=\xi_o^2$ [#22] | 3.3 | B | NOTE | C04 | C04's worked number (ω = 1 s⁻¹, ξ_o = 1 m, t′ = 0: the x-axis and two unit circles centred at (0, ±1), all tangent at the origin); path- and streak-line eliminations written out as D04, D05 (the book skips them); ⚠️ t′ vs t_o vs ξ_o; spatially uniform caricature of wave orbits (Ch. 7) |
| N19 | Fig. 3.7: the three lines of Ex. 3.1 [#23] | 3.3 | B | NOTE | C04 | our figure + animation (dye filament and one particle over a period), `scripts/ch03_fig3_7_flow_lines.py` |
| N20 | Fig. 3.2 discussion: the moving cylinder's streamlines start and end on the body; a particle's acceleration is the same in both frames if U is constant [#24] | 3.3 | B | NOTE | C05 | stated with the two frames' streamline plots and the acceleration arrow at one point (identical) |
| N21 | Galilean transformation: $\mathbf u(\mathbf x,t)=\mathbf U+\mathbf u'(\mathbf x',t')$, $t=t'$, $\mathbf x=\mathbf x'+\mathbf Ut+\mathbf x'_o$ (Fig. 3.8) [#25] | 3.3 | B | NOTE | C05 | stated with a two-frame sketch; `core.K.galilean_transform`; D06 starts from it |
| C05 | Galilean invariance of the acceleration, $\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=\frac{\partial\mathbf u'}{\partial t'}+(\mathbf u'\cdot\nabla')\mathbf u'$, Eq. (3.9) [#26] | 3.3 | A | CORE | – | load-bearing: Newton's law in any inertial frame (Ch. 4), the wave frame that makes waves steady (Ch. 7), and the contrast with rotating frames (Ch. 4 §4.7, Ch. 13) |
| N22 | Local term linear in u, advective term quadratic: the nonlinearity of fluid mechanics; small u → acoustics, u = 0 → statics [#27] | 3.3 | B | NOTE | C05 | stated with a two-line check: scale u by λ = 2, local term ×2, advective ×4 |
| N23 | The local/advective split depends on the frame while the sum does not [#28] | 3.3 | B | NOTE | C05 | stated with term bars in both frames (`core.K.acceleration`) — C05's figure and E3's bars |
| N24 | Traffic-light (unsteady) and roller-coaster (advective) analogies [#29] | 3.3 | C | NOTE | C05 | named in one sentence each; reused as E3 preset captions |
| N25 | Motivation: constitutive laws relate deformation *rates* to stress [#30] | 3.4 | C | NOTE | C06 | named; developed in Ch. 4 §4.5 (Newtonian fluid $\tau=-p\delta+2\mu S$ …) |
| C06 | Relative velocity of a neighbouring point by first-order Taylor expansion, $du_i=(\partial u_i/\partial x_j)\,dx_j$, Eq. (3.10) (Fig. 3.9) [#31] | 3.4 | A | CORE | – | load-bearing: the starting point for strain, spin, (3.19), principal axes, the stress law (Ch. 4), vortex stretching (Ch. 5) and turbulence strain statistics (Ch. 12) |
| R01 | Eq. (3.11): $\frac{\partial u_i}{\partial x_j}=S_{ij}+\tfrac12R_{ij}$ [#32] | 3.4 | B | RECAP | C06 | ch02 C12 (D14 unique split); reminded with `core.tensors.strain_rate_tensor(G) + 0.5*rotation_tensor(G) == G` |
| R02 | Eq. (3.12): strain-rate tensor $S_{ij}=\tfrac12\big(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\big)$ [#33] | 3.4 | B | RECAP | C06 | ch02 C12 / N62; `core.tensors.strain_rate_tensor` |
| R03 | Eq. (3.13): rotation tensor $R_{ij}=\frac{\partial u_i}{\partial x_j}-\frac{\partial u_j}{\partial x_i}$ (book R = G − Gᵀ, no ½) [#34] | 3.4 | B | RECAP | C06 | ch02 C12 (factor-2 convention, ch02 review M1); `core.tensors.rotation_tensor`; ⚠️ R = 2A |
| N26 | Names and roles: S embodies deformation and is linked to stress; R embodies rotation and is not [#35] | 3.4 | C | NOTE | C06 | named; the reason is C11's (3.19) and the stress law of Ch. 4 §4.5 |
| C07 | Linear strain rate $\frac{1}{\delta x_1}\frac{D}{Dt}(\delta x_1)=\frac{\partial u_1}{\partial x_1}$; general direction $\mathbf n\cdot\mathbf S\cdot\mathbf n$ (Fig. 3.10) [#36] | 3.4 | A | CORE | – | load-bearing: the meaning of S's diagonal; normal viscous stress (Ch. 4), vortex stretching (Ch. 5), (3.14) is built from it |
| C08 | Shear strain rate $\tfrac12\frac{D(\alpha+\beta)}{Dt}=\tfrac12\big(\frac{\partial u_1}{\partial x_2}+\frac{\partial u_2}{\partial x_1}\big)=S_{12}$ (Fig. 3.11) [#37] | 3.4 | A | CORE | – | load-bearing: the meaning of S's off-diagonal; shear stress $\mu\,du/dy$ generalised (Ch. 4, Ch. 8 Couette/Poiseuille) |
| N27 | Rigid motion $\mathbf u=\mathbf U+\boldsymbol\Omega\times\mathbf x$ ⇒ $S_{ij}=0$; S independent of a translating or rotating frame (Exercise 3.17) [#38] | 3.4 | B | NOTE | C08 | stated with D10 written out (book never shows it); `ch03.rigid_body_velocity` + `velocity_gradient_at` → S = 0, R = 2[Ω×] for random (U, Ω) |
| C09 | Volumetric strain rate $\frac1{\delta V}\frac{D}{Dt}(\delta V)=\frac{\partial u_i}{\partial x_i}=S_{ii}$, Eq. (3.14) [#39] | 3.4 | A | CORE | – | load-bearing: the kinematic content of continuity $D\rho/Dt=-\rho\nabla\cdot\mathbf u$ (Ch. 4), incompressibility, compressible flow (Ch. 15); closed again by RTT in C15 |
| R04 | R is antisymmetric: zero diagonal, three independent entries ↔ a vector, the vorticity $\boldsymbol\omega=\nabla\times\mathbf u$ [#40] | 3.4 | C | RECAP | C10 | ch02 C12 (N56–N58); one reminder sentence |
| R05 | Eq. (3.15): $R_{ij}=-\varepsilon_{ijk}\omega_k$, matrix $\begin{bmatrix}0&-\omega_3&\omega_2\\\omega_3&0&-\omega_1\\-\omega_2&\omega_1&0\end{bmatrix}$ [#41] | 3.4 | B | RECAP | C10 | ch02 D15 (= (2.26), (2.27)); `core.tensors.antisymmetric_from_vector`, `core.K.vorticity_from_gradient` |
| R06 | Eq. (3.16): $\omega_1=\frac{\partial u_3}{\partial x_2}-\frac{\partial u_2}{\partial x_3}$, $\omega_2=\frac{\partial u_1}{\partial x_3}-\frac{\partial u_3}{\partial x_1}$, $\omega_3=\frac{\partial u_2}{\partial x_1}-\frac{\partial u_1}{\partial x_2}$ [#42] | 3.4 | B | RECAP | C10 | ch02 C11 / D12 (= (2.25)); `core.operators.curl`, `core.K.vorticity` pointwise |
| C10 | Average rotation rate of two perpendicular material lines $\tfrac12\frac{D(-\alpha+\beta)}{Dt}=\tfrac12\omega_3=R_{21}/2$: vorticity is twice the element's spin (Fig. 3.11) [#43] | 3.4 | A | CORE | – | load-bearing: the physical meaning of vorticity used throughout Ch. 5 (vorticity dynamics), Ch. 6 (irrotational flow) and Ch. 13 (relative, planetary, potential vorticity) |
| N28 | Vorticity depends on the frame: $\omega'_z=\omega_z-2\Omega$ in a frame rotating at $\Omega\mathbf e_z$ (Exercise 3.19) [#44] | 3.4 | B | NOTE | C10 | stated with D13 written out (book never shows it); climate hook: absolute vorticity $f+\zeta$ (Ch. 4 §4.7, Ch. 13); `core.K.vorticity_in_rotating_frame` |
| N29 | Eq. (3.17): irrotational flow $\boldsymbol\omega=0$, or equivalently $R_{ij}=\partial u_i/\partial x_j-\partial u_j/\partial x_i=0$; then $\mathbf u=\nabla\phi$ [#45] | 3.4 | B | NOTE | C10 | stated; the converse needs a simply connected region — the line vortex of C13 is the counter-example ($\oint\mathbf u\cdot d\mathbf s=2\pi B$, φ = Bθ multivalued); `core.K.velocity_potential_2d` reports path dependence; potential flow is Ch. 6 |
| R07 | Eq. (3.18): circulation $\Gamma\equiv\oint_C\mathbf u\cdot d\mathbf s=\int_A\boldsymbol\omega\cdot\mathbf n\,dA$; vorticity = circulation per unit area (Fig. 3.12) [#46] | 3.4 | B | RECAP | C10 | ch02 C16 (Stokes, D26, (2.35)); reminded with `core.V.circulation`; reused in C13, C14; Kelvin's theorem Ch. 5, lift Ch. 6/14; ⚠️ Γ here is circulation, not ch02's shear rate |
| C11 | Relative velocity = pure deformation + rigid rotation at $\omega/2$: $du_i=\big(S_{ij}-\tfrac12\varepsilon_{ijk}\omega_k\big)dx_j=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i$, Eq. (3.19) [#47] | 3.4 | A | CORE | – | load-bearing: why only S enters the stress law (Ch. 4), and the strain/rotation decomposition behind vortex stretching (Ch. 5) and Okubo–Weiss-type diagnostics (Ch. 12–13) |
| N30 | Rigid-body velocity $\mathbf v=\boldsymbol\Omega\times\mathbf x$; the second term of (3.19) is a rotation at $\omega/2$ [#48] | 3.4 | C | NOTE | C11 | named in one sentence at D15's last step (primer "rigid-body velocity"); coded by N27 |
| N31 | Eq. (3.20): strain part in the principal frame $d\bar{\mathbf u}=\bar{\mathbf S}\cdot d\bar{\mathbf x}$ with $\bar{\mathbf S}=\mathrm{diag}(\bar S_{11},\bar S_{22},\bar S_{33})$ [#49] | 3.4 | B | NOTE | C12 | stated; the eigenframe is ch02 C13 (principal axes); `core.K.principal_strain_rates` |
| N32 | Eq. (3.21): $d\bar u_1=\bar S_{11}d\bar x_1$, $d\bar u_2=\bar S_{22}d\bar x_2$, $d\bar u_3=\bar S_{33}d\bar x_3$ [#50] | 3.4 | B | NOTE | C12 | stated as step 3 of D16; `core.K.strain_velocity_principal` shows the diagonal action numerically |
| C12 | A small sphere becomes an ellipsoid whose axes are the principal axes of S (Fig. 3.13) [#51] | 3.4 | A | CORE | – | load-bearing: principal strain axes are the frame for principal stresses (Ch. 4), turbulent stretching (Ch. 12) and deformation/frontogenesis fields (Ch. 13) |
| N33 | Summary: relative motion = rotation of the element + deformation [#52] | 3.4 | C | NOTE | C11 | named as the closing sentence of C11 (and of the §3.4 notebook section) |
| N34 | Parallel shear flow $\mathbf u=(u_1(x_2),0)$, $\gamma\equiv du_1/dx_2$, $\omega_3=-\gamma$; AB turns at −γ, BC at 0, average −γ/2 for *any* perpendicular pair [#53] | 3.5 | B | NOTE | C10 | stated at the head of §3.5 with D14 written out (the "any pair" claim the book does not show); ⚠️ clockwise spin, γ = 2S₁₂ vs ch02's Γ ≡ S₁₂; `ch03.parallel_shear_kinematics`; Couette/boundary layers → Ch. 8, 9 |
| N35 | Shear-flow elements: $S=\begin{bmatrix}0&\gamma/2\\\gamma/2&0\end{bmatrix}$, $\bar S=\begin{bmatrix}\gamma/2&0\\0&-\gamma/2\end{bmatrix}$; aligned element shears, 45° element stretches, both spin at −γ/2 (Fig. 3.14) [#54] | 3.5 | B | NOTE | C12 | stated with our Fig. 3.14 (two `deform_square`s side by side, side lengths and corner angles printed); the eigenvalues recap ch02 Ex. 2.4 (N62) |
| N36 | Eq. (3.22): solid-body rotation $u_r=0$, $u_\theta=\omega_0r$ (spun-up tank) [#55] | 3.5 | B | NOTE | C13 | stated; `core.X.solid_body_rotation`, Cartesian field = ch02 `solid_body_rotation_field`; the rotating-tank paraboloid is Ch. 4 |
| C13 | Vorticity in plane polar coordinates $\omega_z=\frac1r\frac{\partial}{\partial r}(ru_\theta)-\frac1r\frac{\partial u_r}{\partial\theta}$, Eq. (3.23): $2\omega_0$ for solid-body rotation, 0 for the irrotational vortex (Figs. 3.15, 3.16) [#56] | 3.5 | A | CORE | – | load-bearing: the paddle-wheel test separating "moves in circles" from "spins"; the polar formula is used for every vortex in Ch. 5, 6, 13 (cyclones, f-plane vortices) |
| N37 | Eq. (3.24): $\Gamma=\oint_C\mathbf u\cdot d\mathbf s=\int_0^{2\pi}u_\theta r\,d\theta=2\pi ru_\theta=2\pi r^2\omega_0$; true for any circuit (Exercise 3.23) [#57] | 3.5 | B | NOTE | C13 | stated; "any circuit" in one sentence (uniform ω ⇒ Γ = ω × area by Stokes); number: ω₀ = 1 s⁻¹, r = 1 m → 2π m²/s; off-centre circle checked by `core.V.circulation` |
| N38 | Eq. (3.25): irrotational (line) vortex $u_r=0$, $u_\theta=B/r$ [#58] | 3.5 | B | NOTE | C13 | stated; `core.X.line_vortex`; ch02's K/r test field; ⚠️ Fig. 3.16 writes C for B |
| N39 | Eq. (3.26): $\Gamma=\int_0^{2\pi}u_\theta r\,d\theta=2\pi ru_\theta=2\pi B$ for every r; $\omega_z=0$ off the axis by (3.23) [#59] | 3.5 | B | NOTE | C13 | stated; `core.X.circulation_circle` for r = 0.5, 1, 2 m gives the same 2πB |
| N40 | Eq. (3.27): $[\omega_z]_{r\to0}=\lim_{r\to0}\frac{1}{\pi r^2}\oint_C\mathbf u\cdot d\mathbf s=\lim_{r\to0}\frac{2B}{r^2}$ — a δ-function core [#60] | 3.5 | B | NOTE | C13 | stated with a log–log plot of mean vorticity in a disc vs r (slope −2), `core.X.mean_vorticity_in_disc`; δ-function gloss (ch02 D21); point vortex Ch. 6 |
| N41 | $\Gamma_{ABCD}=-[u_\theta r]_r\Delta\theta+[u_\theta r]_{r+\Delta r}\Delta\theta=0$ for a sector excluding the origin: elements deform but do not spin (Fig. 3.16) [#61] | 3.5 | B | NOTE | C13 | stated; it is D17's four-side sum evaluated for $u_\theta=B/r$; `ch03.annular_sector_circulation`; the E6 paddle wheel keeps its orientation |
| N42 | Real vortices (bathtub, wing tip, tornado): solid-body core, irrotational far field, bounded speed [#62] | 3.5 | C | NOTE | C14 | named as C14's opening scene; tropical cyclones Ch. 13, wing-tip vortices Ch. 14 |
| C14 | Rankine vortex $\omega_z=\Gamma/\pi\sigma^2$ ($r\le\sigma$), 0 ($r>\sigma$); $u_\theta=(\Gamma/2\pi\sigma^2)r$ ($r\le\sigma$), $\Gamma/2\pi r$ ($r>\sigma$), Eq. (3.28) [#63] | 3.5 | A | CORE | – | load-bearing: the standard model of a real vortex (core + irrotational outside) reused for vortex dynamics (Ch. 5), cyclones and gradient wind (Ch. 13), tip vortices (Ch. 14) |
| N43 | Eq. (3.29): Gaussian vortex $\omega_z=\frac{\Gamma}{\pi\sigma^2}\exp(-r^2/\sigma^2)$, $u_\theta=\frac{\Gamma}{2\pi r}\big(1-\exp(-r^2/\sigma^2)\big)$ [#64] | 3.5 | B | NOTE | C14 | stated with D19 written out ($u_\theta$ from $\omega_z$; the book never shows it); `core.X.gaussian_vortex` (`expm1` form); Lamb–Oseen pointer $\sigma^2=4\nu t$ (Ch. 5, Ch. 8) |
| N44 | Gaussian-vortex maximum where $1+2r^2/\sigma^2=\exp(r^2/\sigma^2)$, $r\approx1.1209\sigma$ (Exercise 3.26); Rankine maximum at σ [#65] | 3.5 | B | NOTE | C14 | stated with D20 written out; `core.X.gaussian_vortex_max_radius` (brentq; Lambert-W cross-check); number: Γ = 2π m²/s, σ = 1 m → $u_{\theta,\max}\approx0.638$ m/s |
| N45 | RTT motivation: time derivative of integrals over moving, deforming volumes [#66] | 3.6 | C | NOTE | C15 | named as C15's opening question; used by every integral law of Ch. 4 |
| N46 | Eq. (3.30): Leibniz's theorem $\frac{d}{dt}\int_{a(t)}^{b(t)}F(x,t)\,dx=\int_a^b\frac{\partial F}{\partial t}dx+\frac{db}{dt}F(b,t)-\frac{da}{dt}F(a,t)$ [#67] | 3.6 | B | NOTE | C15 | stated as the 1-D RTT with D21 written out (the book cites a proof it does not give); `core.R.leibniz_terms`, `leibniz_check`; momentum integral Ch. 9, layer budgets Ch. 13 |
| N47 | Fig. 3.17: graphical Leibniz (interior change + gain at b − loss at a) [#68] | 3.6 | B | NOTE | C15 | our three-strip figure and animation of a(t), b(t) (`scripts/ch03_fig3_17_leibniz.py`) |
| N48 | Control volume V*(t), control surface A*(t), outward normal n, surface velocity b (need not follow the fluid) (Fig. 3.18) [#69] | 3.6 | B | NOTE | C15 | stated; `core.R.ControlVolume` shapes (moving box, growing cylinder/sphere/cone); b = u material volume, b = 0 fixed volume |
| N49 | Eq. (3.31): $\frac{d}{dt}\int_{V^*(t)}F\,dV=\lim_{\Delta t\to0}\frac1{\Delta t}\{\int_{V^*(t+\Delta t)}F(\mathbf x,t+\Delta t)dV-\int_{V^*(t)}F(\mathbf x,t)dV\}$ [#70] | 3.6 | B | NOTE | C15 | D22's starting line; `core.R.volume_integral_rate_fd` is the independent "left side" |
| N50 | Eq. (3.32): $\int_{V^*(t+\Delta t)}F(\mathbf x,t+\Delta t)dV\cong\int_{V^*}F\,dV+\int_{V^*}\Delta t\frac{\partial F}{\partial t}dV+\int_{\Delta V}F\,dV+\int_{\Delta V}\Delta t\frac{\partial F}{\partial t}dV$ [#71] | 3.6 | B | NOTE | C15 | a step of D22 (shown in full there); the four terms printed for a growing sphere at Δt = 10⁻¹…10⁻⁴ |
| N51 | Eq. (3.33): $\frac{d}{dt}\int_{V^*}F\,dV=\lim_{\Delta t\to0}\frac1{\Delta t}\{\int_{V^*}\Delta t\frac{\partial F}{\partial t}dV+\int_{\Delta V}F\,dV\}$ [#72] | 3.6 | B | NOTE | C15 | a step of D22; log–log plot: the dropped term falls with slope 2 (orders of smallness made visible) |
| N52 | Eq. (3.34): swept sliver $\int_{\Delta V}F\,dV\cong\int_{A^*}F\,(\mathbf b\Delta t\cdot\mathbf n)\,dA$ as Δt → 0 [#73] | 3.6 | B | NOTE | C15 | a step of D22; ⚠️ signed volume (b·n < 0 removes volume); `core.R.swept_volume_integral` vs `surface_flux_term`·Δt, gap O(Δt²) |
| N53 | Interpretations: F = 1 gives $dV^*/dt=\int_{A^*}\mathbf b\cdot\mathbf n\,dA$ and, for a small material volume, (3.14) (Ex. 3.28); (3.35) extends (3.5) to finite volumes (Ex. 3.30); b·n > 0 advancing [#75] | 3.6 | B | NOTE | C15 | stated with D23 and D24 written out (the book never shows them; D24 exposes the hidden $F\nabla\cdot\mathbf u$ term); `core.R.material_volume_rate` |
| N54 | Ex. 3.2: growing cone of fixed height h, base radius $r(t)$: $V=\tfrac13\pi hr^2$, $dV/dt=\int_{A^*}\mathbf b\cdot\mathbf n\,dA=\tfrac23\pi hr_o\dot r$ (Fig. 3.19) [#76] | 3.6 | B | NOTE | C15 | C15's worked number (h = 1 m, r₀ = 0.5 m, ṙ = 0.1 m/s → 0.1047 m³/s both ways); ⚠️ b·n = 0 on the base (not b = 0), "[?]" = z; `ch03.example_3_2` |
| N55 | Fig. 3.18: moving control volume, swept sliver $(\mathbf b\Delta t)\cdot\mathbf n\,dA$ [#77] | 3.6 | B | NOTE | C15 | our 2-D deforming blob with b arrows and the swept band coloured by sign(b·n) + waterfall bars (`scripts/ch03_fig3_18_rtt.py`) — E7's stage |
| C15 | Reynolds transport theorem $\frac{d}{dt}\int_{V^*(t)}F\,dV=\int_{V^*(t)}\frac{\partial F}{\partial t}dV+\int_{A^*(t)}F\,\mathbf b\cdot\mathbf n\,dA$, Eq. (3.35) [#74] | 3.6 | A | CORE | – | load-bearing: the tool that turns every conservation law of Ch. 4 (mass, momentum, energy) into integral and differential form; also Ch. 9 (momentum integral) and Ch. 13 (layer budgets) |
| S01 | Exercises 3.1–3.30 [#78] | Ex. | C | SKIP | – | pointer line: practise on the book's exercises; the derivations the text defers to Exercises 3.3, 3.12, 3.17–3.20, 3.23, 3.26, 3.28, 3.30 are written out in our words (§4b); exercise text is not reproduced |
| S02 | Literature cited and supplemental reading [#79] | Lit. | C | SKIP | – | pointer line: bibliography (Van Dyke's *Album of Fluid Motion* and Samimy et al. for flow pictures; Riley et al. for Leibniz) |

### 2a. What each A block contains (for the lesson-designer)
- **C01** picture: a float drifting past a fixed thermometer · question: "two observers, one flow — how do their
  numbers agree?" · D01 · number: $x=Xe^{\alpha t}$, X = 2 m, α = 0.5 s⁻¹, t = 1 s → x = 3.297 m, u = 1.649 m/s = αx, a =
  0.824 m/s² · code: `ch03.lagrangian_map_example`, `core.K.lagrangian_velocity_acceleration`, `lagrangian_to_eulerian`
  · figure: space–time diagram (x vs t) of several labelled particles with the Eulerian u(x) at two instants · explainer
  E2 · B/C: N01–N03, N05–N08 (§3.1 section), N07, N08.
- **C02** picture: a weather balloon vs a weather station · question: "why does the station warm while the air parcel
  does not?" · D02 · number: T falls 1 K per 100 km northward, southerly wind 10 m/s: $v\,\partial T/\partial y=-10^{-4}$
  K/s; parcel conserves T (DT/Dt = 0) ⇒ station sees $\partial T/\partial t=+10^{-4}$ K/s = 0.36 K/h · code:
  `core.K.material_derivative`, `material_derivative_terms`, `material_derivative_sym` + from-scratch stencil · figure:
  term bars along a path + `slider_figure` of the three terms vs time · explainer E2 · B/C: N09–N12.
- **C03** picture: iron filings / a long-exposure snapshot · question: "which way does the flow point *right now*?" ·
  D03 · number: u = (1, 2) at a point ⇒ slope 2; solid-body rotation streamlines are circles · code: `core.K.streamline`,
  `ch03.streamline_slope` · figure: streamlines of Ex. 3.1 frozen at three instants (the pattern turns) · explainer E1 ·
  B/C: N14, N15.
- **C04** picture: a dye port in a sloshing tank · question: "why do the photograph of dye and the track of one particle
  disagree?" · D04, D05 · number: Ex. 3.1 with ω = 1, ξ_o = 1, t′ = 0 (N18) · code: `core.K.pathline`, `streakline`,
  `ch03.example_3_1` + from-scratch RK4 path line · figure: Fig. 3.7 reproduction + animation · explainer E1 · B/C:
  N13, N16–N19.
- **C05** picture: a cylinder towed through a still lake vs a river past a pier · question: "is the flow steady?
  depends who asks — is the acceleration?" · D06 · number: lab frame with U = 2 m/s and $u'=0.5\sin x'$ (k = 1 m⁻¹) at
  $x'=\pi/4$: local −0.707, advective +0.832, sum 0.125 m/s²; wave frame: local 0, advective 0.125 m/s² · code:
  `core.K.galilean_transform`, `core.K.acceleration`, `ch03.cylinder_flow` · figure: two-frame streamlines + term bars
  (N23) · explainer E3 · B/C: N04, N20–N24.
- **C06** picture: two nearby corks and their velocity difference · D07 · number: G = [[1, 2], [0, −1]] s⁻¹, dx = (0.01,
  0.02) m → du = (0.05, −0.02) m/s; Taylor remainder slope 2 · code: `core.K.velocity_gradient_at`, `relative_velocity`
  · figure: du = G·dx arrows on a ring of neighbours (the "relative velocity field") · explainer E4 · B/C: N25, R01–R03,
  N26.
- **C07** picture: a stretched rubber band in a flow · D08 · number: u = (2x, −2y): a 1 cm segment along x becomes
  1.0202 cm after 0.01 s (rate 2 s⁻¹) · code: `core.K.linear_strain_rate` + from-scratch segment tracking with
  `linear_flow_map` · figure: measured (1/ℓ)dℓ/dt vs n·S·n over all directions (the rose curve) · explainer E4.
- **C08** picture: a right-angled corner closing · D09, D10 · number: u = (γy, 0), γ = 1 s⁻¹, dt = 0.01 s: dα = 0.01
  rad, dβ = 0 → S₁₂ = 0.5 s⁻¹ · code: `core.K.shear_strain_rate`, `ch03.rigid_body_velocity` · figure: tracked corner
  angle vs t with the slope −2S₁₂ · explainer E4 · B/C: N27.
- **C09** picture: a balloon of fluid swelling · D11 · number: u = (x, y, z) s⁻¹: ∇·u = 3 s⁻¹; 1 cm³ box → 1.0305 cm³
  after 0.01 s ($e^{0.03}$) · code: `core.K.volumetric_strain_rate`, `material_volume_ratio` · figure: tracked box volume
  vs t against $e^{t\,\mathrm{tr}G}$ · explainer E4.
- **C10** picture: a paddle wheel in a river · D12, D13, D14 · number: solid body ω₀ = 1: every line turns at 1 rad/s,
  ω₃ = 2, spin 1; shear γ = 1: lines turn at −sin²θ, pair average −0.5 = ω₃/2 · code: `core.K.element_rotation_rate`,
  `material_line_rotation_rate`, `vorticity_in_rotating_frame` · figure: θ̇(θ) curve with the pair average flat at
  ½ω₃ (`slider_figure` in γ) · explainer E5 · B/C: R04–R07, N28, N29, N34.
- **C11** picture: one relative-velocity arrow split into two · D15 · number: shear γ = 1, dx = (0, 1): du = (1, 0) =
  S·dx (0.5, 0) + ½ω × dx (0.5, 0) · code: `core.K.relative_velocity_split` + from-scratch ε loop · figure: ring of
  neighbours with three arrow sets (total, strain, rotation) · explainer E5 · B/C: N30, N33.
- **C12** picture: a droplet of dye becoming an ellipse · D16 · number: shear γ = 1: principal rates ±0.5 s⁻¹ at ±45°;
  1 mm circle after 0.1 s → semi-axes ≈ 1.05 and 0.95 mm (exact $e^{\pm0.05}$ = 1.051, 0.951 for the strain part) ·
  code: `core.K.principal_strain_rates`, `strain_ellipse_axes`, `deform_circle` · figure: circle → ellipse with eigen-axes
  + animation (and the finite-time drift of the axes, analysis §9) · explainer E5 · B/C: N31, N32, N35.
- **C13** picture: a paddle wheel on a merry-go-round vs in a drain vortex · D17 · number: $u_\theta=r$ (ω₀ = 1):
  (1/r)d(r²)/dr = 2; $u_\theta=1/r$: (1/r)d(1)/dr = 0, Γ = 2π · code: `core.X.polar_vorticity_z(_sym)`,
  `solid_body_rotation`, `line_vortex`, `circulation_circle`, `mean_vorticity_in_disc` + from-scratch four-side sum ·
  figure: side-by-side animation of elements (Figs. 3.15, 3.16) · explainer E6 · B/C: N36–N41.
- **C14** picture: a tornado's wind-speed profile · D18, D19, D20 · number: Γ = 2π m²/s, σ = 1 m: Rankine inside
  $u_\theta=r$, ω = 2 s⁻¹, max 1 m/s at r = 1 m; Gaussian max 0.638 m/s at 1.121 m · code: `core.X.rankine_vortex`,
  `gaussian_vortex`, `gaussian_vortex_max_radius` · figure: $u_\theta(r)$ and $\omega_z(r)$ for both, marker at 1.1209σ,
  `slider_figure` in σ · explainer E6 · B/C: N42–N44.
- **C15** picture: a balloon being inflated in a warming room · D21, D22, D23, D24 · number: growing sphere R = 1 m, ṙ =
  0.1 m/s, F = t at t = 1 s: volume term V = 4.189, surface term $F\,\dot r\,4\pi R^2$ = 1.257, total 5.445 (= d/dt(tV)
  by finite difference); Ex. 3.2 (N54) · code: `core.R.leibniz_terms`, `reynolds_transport`, `rtt_check`,
  `volume_integral_rate_fd`, `material_volume_rate`, `ch03.example_3_2` + from-scratch midpoint sums · figure: Fig. 3.17
  and Fig. 3.18 reproductions, Δt-convergence log–log plot · explainer E7 · B/C: N45–N55.

## 3. Section coverage
| § | Title | A | B | C | RECAP | SKIP |
|---|---|---|---|---|---|---|
| 3.1 | Introduction and Coordinate Systems | — (inside C01; N04 inside C05) | N02, N03, N04, N05 | N01, N06 | — | — |
| 3.2 | Particle and Field Descriptions of Fluid Motion | C01, C02 | N07, N08, N09, N10, N11, N12 | — | — | — |
| 3.3 | Flow Lines, Fluid Acceleration, and Galilean Transformation | C03, C04, C05 | N13, N14, N15, N16, N17, N18, N19, N20, N21, N22, N23 | N24 | — | — |
| 3.4 | Strain and Rotation Rates | C06, C07, C08, C09, C10, C11, C12 | N27, N28, N29, N31, N32 | N25, N26, N30, N33 | R01, R02, R03, R05, R06, R07 (B); R04 (C) | — |
| 3.5 | Kinematics of Simple Plane Flows | C13, C14 | N34, N35, N36, N37, N38, N39, N40, N41, N43, N44 | N42 | — | — |
| 3.6 | Reynolds Transport Theorem | C15 | N46, N47, N48, N49, N50, N51, N52, N53, N54, N55 | N45 | — | S01, S02 (end of chapter) |

Totals: A 15 · B 52 (46 NOTE + 6 RECAP) · C 12 (9 NOTE + R04 + S01 + S02) = 79. No section is empty; §3.1 has no A item
(decision 1) and follows ch02 decision 6 (own short notebook section, paragraphs tagged C01 / C05).

## 4. Prerequisites needing primers (concept or tool | needed by | why it is not A/B/RECAP)
Rows marked **gloss** are one plain sentence where the item is used; all others are full 📎 primers (one or two plain
sentences, what it means here, a 2–4-line runnable demo with easy numbers). Items already in `knowledge/primers.md`
get a one-line reminder naming the earlier primer (listed at the end). New primers continue at P87.

| Concept or tool | Needed by | Why it is not A/B/RECAP |
|---|---|---|
| functions of time with parameters; a label vs a variable; inverse function (solve $x=Xe^{\alpha t}$ for X) and sympy `solve` | C01 (D01) | pre-book algebra + Python tool; P40 primed sympy `diff`/`simplify` only |
| multivariable chain rule along a path $\frac{d}{dt}F(x(t),y(t),z(t),t)$ (extends P49's two-variable form to four) | C02 (D02), C05 (D06) | P49 primed $df=(\partial f/\partial x)dx+(\partial f/\partial y)dy$; the "along a trajectory, with t also explicit" form is new |
| chain rule with a moving frame: $\partial/\partial t$ at fixed x vs fixed x′ when $\mathbf x'=\mathbf x-\mathbf Ut$ | C05 (D06) | new calculus move; the classic trap of the chapter |
| frames of reference and relative velocity (inertial frame, "same clock", constant U) | C05, N21 | physics vocabulary not primed in ch01/ch02 |
| parametric curves, tangent vector, arc-length parametrisation $d\mathbf x/ds=\mathbf u/\lvert\mathbf u\rvert$ | C03 (D03), C04 | pre-book calculus |
| parallel vectors: $\mathbf a\parallel\mathbf b\iff\mathbf a=\lambda\mathbf b\iff\mathbf a\times\mathbf b=0$ | C03 (D03) | pre-book; cross product itself is ch02 C08 (reminder) |
| systems of ODEs as one vector ODE; `solve_ivp` options `t_eval`, `dense_output`, `events` (terminal event at a stagnation point), integrating backwards | C03, C04, N17 | P31 primed `solve_ivp` with defaults only |
| RK4 by hand (from-scratch path line) | C04 | P30 primed Euler/Euler–Cromer only |
| eliminating a parameter with $\sin^2+\cos^2=1$ (circle equation from $x(t), y(t)$) | C04 (D04, D05) | pre-book trig — **gloss** inside D04 |
| implicit differentiation (slope of a circle at a point) | D05 (tangency at the origin) | pre-book calculus — **gloss** |
| multivariable first-order Taylor expansion $u_i(\mathbf x+d\mathbf x)\approx u_i(\mathbf x)+(\partial u_i/\partial x_j)dx_j$, remainder $O(\lvert d\mathbf x\rvert^2)$ | C06 (D07), C07–C10 | P26 primed the 1-D form only |
| material line element $\delta\mathbf x$ carried by the flow; $D(\delta\mathbf x)/Dt=\delta\mathbf u$ | C07–C10 | new physics idea used by every §3.4 derivation |
| small-angle approximation $\tan d\alpha\approx d\alpha$, $\cos d\alpha\approx1$ (errors second order) | C08 (D09), C10 (D12) | pre-book trig |
| angular velocity of a line; counter-clockwise positive | C10 (D12, D14) | pre-book kinematics |
| rigid-body velocity $\mathbf v=\boldsymbol\Omega\times\mathbf x$ | N27 (D10), C11 (D15), N28 (D13) | pre-book mechanics; cross product is ch02 C08 |
| rotating frame of reference (first time): $\mathbf u=\boldsymbol\Omega\times\mathbf x+\mathbf u'$ at the instant the frames coincide | N28 (D13) | new; Ch. 4 §4.7 develops it fully |
| product rule for three factors $D(abc)=\dot a\,bc+a\dot b\,c+ab\dot c$ | C09 (D11) | P38 primed two factors — **gloss** |
| Jacobi's formula $\det e^{Gt}=e^{t\,\mathrm{tr}\,G}$ (volume factor of a linear flow) | C09 | pre-book linear algebra — **gloss** with a numeric check (`expm`, P79) |
| linear map of a circle is an ellipse; semi-axes from eigenvalues of a symmetric map; SVD (**gloss**) for finite time | C12 (D16) | pre-book linear algebra; eigenvalues are P80 (reminder) |
| polar coordinates as a moving basis: $\mathbf e_r,\mathbf e_\theta$ depend on θ; area element $r\,dr\,d\theta$; line element $d\mathbf s=r\,d\theta\,\mathbf e_\theta$ on a circle | C13 (D17), N37, N39 | ch02 Ex. 2.1 did components only |
| cylindrical and spherical unit vectors as functions of position (orthonormal, right-handed) | N05 | new coordinate systems (plane polar SEEN) |
| Dirac delta: infinite density, finite total | N40 | ch02 D21 **gloss** — reminder sentence |
| piecewise functions and continuity at a join | C14 (D18) | pre-book; `np.where` is P46 (reminder) |
| substitution in an integral ($s=r'^2/\sigma^2$, $ds=2r'dr'/\sigma^2$) | D19 | pre-book calculus |
| maximum of a function: derivative = 0, excluding the trivial root | D20 | pre-book calculus — **gloss** |
| `scipy.optimize.brentq` root finding with a bracket | N44 (D20) | Python tool (analysis flags it as new) |
| `np.expm1` and catastrophic cancellation near r = 0 | N43 (`gaussian_vortex`) | Python / numerics tool |
| Lambert W function (cross-check only) | N44 | **gloss** |
| differentiation under the integral sign $\frac{d}{dt}\int_a^bF\,dx=\int_a^b\partial F/\partial t\,dx$ for fixed limits | N46 (D21) | pre-book analysis |
| volume of a thin shell swept by a moving surface, $(\mathbf b\Delta t\cdot\mathbf n)\,dA$, as a **signed** volume | C15 (D22) | new geometric idea; the chapter's most common confusion |
| Taylor expansion in time $F(\mathbf x,t+\Delta t)\approx F+\Delta t\,\partial F/\partial t$ | D22 | P26 (reminder) applied in t — **gloss** |
| `scipy.integrate.quad` and `dblquad` | N03, N54, C15 | Python tools (P37 primed the trapezoid rule only) |
| divergence theorem applied to $F\mathbf u$ and the product rule $\nabla\cdot(F\mathbf u)=\mathbf u\cdot\nabla F+F\nabla\cdot\mathbf u$ | D23, D24 | Gauss is ch02 C14 (reminder); the product rule for ∇· is new — **gloss** |
| streamlines in plotly 3-D (`go.Streamtube` or `go.Scatter3d` lines seeded on a circle) | N15 | P64 primed cones, lines, meshes — **gloss** in a code comment |

Reminders only (already primed): partial derivative P25 · chain rule P49 · finite differences P21/P22 · first-order
Taylor P26 · definite integral P27 · trapezoid rule P37 · product rule P38 · sympy P40 · `solve_ivp` P31 · lambda P29 ·
`assert np.allclose` P15 · `animate` P16 · `slider_figure` P17 · `show_viz` P18 · `np.where` P46 · log–log slope P13 ·
limits and orders of smallness P68 · mean-value theorem P85 · fundamental theorem of calculus P84 · volume/surface
integrals as midpoint sums P83 · line integral around a loop P86 · eigenvalues P80 · `scipy.linalg.expm` P79 ·
`np.meshgrid` and the grid layout P76 · broadcasting P77 · `plt.streamplot` P78 · `np.arctan2` P70 · right-hand rule P74
· directional derivative P75 · plotly 3-D P41/P64 · `np.einsum` P62.

## 4b. Derivations written out (parsed by tools: ID first, CORE id in a column, ★★★ for hard, explainer slugs backticked in the LAST column)
Written out because (a) the result belongs to an A item, or (b) the book never writes it out and the lesson needs it
(marked "(b) book never writes it out"). One small move per step; the book's skipped moves (analysis §2b) are filled in
and the notebook says so. Analysis §2b numbers in brackets. 24 rows.

| ID | Result (Eq.) | CORE | Difficulty | Steps | Tools used | Traps | Shown in |
|---|---|---|---|---|---|---|---|
| D01 | Eulerian velocity of a Lagrangian map, $x=Xe^{\alpha t}$ ⇒ $u=\alpha x$, and $Du/Dt=d^2x/dt^2$ — (b) book never writes it out [a-D01] | C01 | ★ | 6 | (3.1) (N08), (3.2), inverse function + sympy `solve` (primer), derivative | differentiating with respect to t while x is held fixed vs X held fixed; forgetting that the label X is constant along the path; the map must be invertible | notebook · `material_derivative_probe` |
| D02 | material derivative (3.3) → (3.4) → (3.5): $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$ [a-D02] | C02 | ★★ | 9 | compatibility (3.2) (C01), multivariable chain rule along a path (primer), $dr_i/dt=u_i$ (3.1), partial derivative P25, dot product | why $\mathbf r_o,t_o$ are held fixed (labels); why an identity derived on one path holds at every x (every point hosts a particle); D/Dt of a vector acts component by component (needed in D06) | notebook · `material_derivative_probe` |
| D03 | streamline equations $dx/u=dy/v=dz/w$ (3.7) from tangency — Exercise 3.3, (b) book never writes it out [a-D04] | C03 | ★ | 5 | parallel vectors (primer), cross product in components (ch02 C08 reminder), arc length (primer) | the ratio form breaks where a component vanishes (the cross-product form does not); t is frozen, not integrated | notebook · `flow_lines_unsteady` |
| D04 | Ex. 3.1 path line: circle of radius $\xi_o$ centred at $(-\xi_o\sin\omega t',\ \xi_o\cos\omega t')$ — elimination (b) the book skips ("a little algebra") [a-D06] | C04 | ★★ | 7 | (3.8) (C04), integration of sin/cos, fixing constants with $\mathbf r(t')=0$, $\sin^2+\cos^2=1$ (gloss) | four names that collide: t′ (drawing instant), t_o (release), $x_o,y_o$ (constants), $\xi_o$ (amplitude); direction of travel (counter-clockwise); u is spatially uniform (a caricature of wave orbits) | notebook · `flow_lines_unsteady` |
| D05 | Ex. 3.1 streak line: circle centred at $(\xi_o\sin\omega t',\ -\xi_o\cos\omega t')$; all three lines tangent at the origin — elimination and tangency (b) the book skips [a-D07] | C04 | ★★ | 8 | D04 with the release time $t_o$ as the label, elimination of $t_o$, implicit differentiation (gloss) | the parameter being eliminated is $t_o$, not t; which part of the circle is filled (release times over one period fill it all); the tangency claim needs the three slopes at 0, all $\tan\omega t'$ | notebook · `flow_lines_unsteady` |
| D06 | Galilean invariance of the acceleration (3.9) — Exercise 3.12, (b) book never writes it out [a-D08] | C05 | ★★ | 9 | Galilean transformation (N21), chain rule with a moving frame (primer), $\nabla=\nabla'$ for a pure translation, D/Dt of a vector (D02) | $\partial\mathbf u/\partial t$ at fixed x ≠ $\partial\mathbf u'/\partial t'$ at fixed x′ (the $-\mathbf U\cdot\nabla'$ term); valid only for constant U; a rotating frame adds Coriolis/centrifugal terms (Ch. 4 §4.7) | notebook · `galilean_frames_cylinder` |
| D07 | relative velocity $du_i=(\partial u_i/\partial x_j)dx_j$ (3.10) by first-order Taylor expansion [a-D09] | C06 | ★ | 4 | multivariable Taylor (primer), index notation (ch02 C01 reminder), G is a tensor (ch02 N29) | the dropped remainder is $O(\lvert d\mathbf x\rvert^2)$, not zero; G[i, j] = ∂u_i/∂x_j (row = velocity component) | notebook · `fluid_element_deformation` |
| D08 | linear strain rate $\frac1{\delta x_1}\frac{D(\delta x_1)}{Dt}=\partial u_1/\partial x_1$; any direction $\mathbf n\cdot\mathbf S\cdot\mathbf n$ [a-D11] | C07 | ★ | 6 | material line element (primer), D07, limit of a difference quotient (P68 reminder) | A′B′ = AB + BB′ − AA′ bookkeeping; transverse velocity differences only rotate AB (length change second order); only S (not R) survives in $\mathbf n\cdot\mathbf G\cdot\mathbf n$ | notebook · `fluid_element_deformation` |
| D09 | shear strain rate $S_{12}=\tfrac12D(\alpha+\beta)/Dt$; general $\mathbf n_1\cdot\mathbf S\cdot\mathbf n_2$ [a-D12] | C08 | ★★ | 8 | D07, small-angle approximation (primer), angle conventions of Fig. 3.11 | α clockwise from the vertical side, β counter-clockwise from the horizontal side; stretching changes the angles only at second order; the ½ (S is *half* the closing rate) | notebook · `fluid_element_deformation` |
| D10 | rigid motion $\mathbf U+\boldsymbol\Omega\times\mathbf x$ ⇒ S = 0, R = 2[Ω×], ω = 2Ω; S unchanged by a rigid change of frame — Exercise 3.17, (b) book never writes it out [a-D14] | C08 | ★ | 5 | rigid-body velocity (primer), cross product in components (ch02 C08), antisymmetric matrix of a vector (ch02 D15 recap) | U and Ω must be uniform in space; "S frame-independent" means adding a rigid motion changes G by an antisymmetric part only | notebook |
| D11 | volumetric strain rate $\frac1{\delta V}\frac{D(\delta V)}{Dt}=S_{ii}$ (3.14) — Exercise 3.18, (b) book never writes it out [a-D15] | C09 | ★★ | 7 | product rule for three factors (gloss), D08 per side, trace invariance (ch02 C07 reminder), Jacobi's formula (gloss) | shear changes volume only at second order; the Greek index means no sum; exact finite-time factor is $e^{t\,\mathrm{tr}G}$, not $1+t\,\mathrm{tr}G$ | notebook · `fluid_element_deformation` |
| D12 | element rotation rate $\tfrac12D(-\alpha+\beta)/Dt=\tfrac12\omega_3=R_{21}/2$ [a-D13] | C10 | ★★ | 8 | D09's construction, small-angle approximation (primer), angular velocity of a line (primer), (3.15) (R05) | why the *average* of two perpendicular lines defines spin; the sign of α (clockwise); $-R_{12}/2=R_{21}/2$ | notebook · `spin_and_principal_axes` |
| D13 | frame dependence of vorticity $\omega'_z=\omega_z-2\Omega$ — Exercise 3.19, (b) book never writes it out [a-D17] | C10 | ★★ | 6 | rotating frame (primer), $\nabla\times(\boldsymbol\Omega\times\mathbf x)=2\boldsymbol\Omega$ (ch02 Ex. 2.3 reminder), curl is a vector (ch02 C03) | computing the curl at the instant the frames coincide; ω′ = 0 iff Ω = ω_z/2 (the co-rotating frame); do not confuse with Galilean invariance (acceleration *is* frame-dependent here) | notebook · `spin_and_principal_axes` |
| D14 | parallel shear flow: a line at angle θ turns at $\dot\theta=-\gamma\sin^2\theta$; any perpendicular pair averages $-\gamma/2=\omega_3/2$ — (b) book never writes it out [a-D20] | C10 | ★★ | 6 | angular velocity of a line (primer), $\dot\theta=\mathbf e_\theta\cdot\mathbf G\cdot\mathbf e(\theta)$, $\sin^2\theta+\cos^2\theta=1$ | $\omega_3=-\gamma$ is clockwise; θ + π/2 turns at $-\gamma\cos^2\theta$; γ = 2S₁₂ (ch02's Γ ≡ S₁₂ differs by 2) | notebook · `spin_and_principal_axes` |
| D15 | deformation + rigid rotation $du_i=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i$ (3.19) [a-D19] | C11 | ★★ | 6 | (3.10) (C06), (3.11) (R01), (3.15) (R05), ε index moves and (2.21) (ch02 C08 reminder), rigid-body velocity (primer) | one ε index swap flips the sign: $\varepsilon_{ijk}\omega_kdx_j=-(\boldsymbol\omega\times d\mathbf x)_i$; the book cites (3.14) where (3.15) is meant; angular velocity is ω/2, not ω | notebook · `spin_and_principal_axes` |
| D16 | a small sphere becomes an ellipsoid with semi-axes $(1+\bar S_{\alpha\alpha}dt)\lvert d\mathbf x\rvert$ along the principal axes — (b) the book states it without proof [a-D22] | C12 | ★★ | 8 | (3.19) (C11), (3.20)–(3.21) (N31, N32), linear map of a circle (primer), eigenvectors (P80 reminder) | first order in dt only; the rotation part turns the ellipsoid by ω dt/2 but does not reshape it; for finite t the axes follow the singular vectors of $e^{Gt}$ and drift away from 45° in a shear flow | notebook · `spin_and_principal_axes` |
| D17 | polar vorticity $\omega_z=\frac1r\frac{\partial}{\partial r}(ru_\theta)-\frac1r\frac{\partial u_r}{\partial\theta}$ (3.23) from circulation around a polar sector — Exercise 3.20, (b) book cites Appendix B [a-D23] | C13 | ★★ | 10 | circulation = vorticity × area (R07, ch02 D26), polar area and line elements (primer), first-order Taylor in r and θ (P26 reminder) | the arc lengths $r\,d\theta$ and $(r+dr)\,d\theta$ differ, which is where $\partial(ru_\theta)/\partial r$ comes from; orientation of the four legs (AB and CD against the loop); the radial legs carry $u_r$ | notebook · `vortex_paddle_wheels` |
| D18 | Rankine vortex (3.28) consistency: ω from $u_\theta$ by (3.23), continuity at σ, total Γ, maximum at σ — (b) book never writes it out [a-D26] | C14 | ★ | 5 | (3.23) (C13), derivative of a product, piecewise functions (primer) | $u_\theta$ is continuous at σ but $\omega_z$ jumps; the maximum is a kink (derivative changes sign), not a zero of the derivative | notebook · `vortex_paddle_wheels` |
| D19 | Gaussian vortex (3.29): $\Gamma(r)=\Gamma(1-e^{-r^2/\sigma^2})$ and $u_\theta=\Gamma(r)/2\pi r$ from $\omega_z$; limits r ≪ σ and r ≫ σ — (b) book never writes it out [a-D27] | C14 | ★★ | 7 | Stokes on a circle (R07), substitution in an integral (primer), exponential P36, first-order Taylor of $e^{-s}$ (P26 reminder) | the 2πr′ weight in the area integral; the small-r limit is solid body with $\omega_0=\Gamma/2\pi\sigma^2$, not zero vorticity; `expm1` near r = 0 | notebook · `vortex_paddle_wheels` |
| D20 | Gaussian-vortex maximum $1+2x=e^x$, $x=r^2/\sigma^2\approx1.2564$, $r\approx1.1209\sigma$ — Exercise 3.26, (b) book never writes it out [a-D28] | C14 | ★★ | 7 | maximum ⇒ derivative 0 (gloss), product/chain rule, `brentq` (primer), Lambert W (gloss) | x = 0 is also a root (exclude it); differentiate with respect to r or x consistently; r/σ = √x, not x | notebook · `vortex_paddle_wheels` |
| D21 | Leibniz's theorem (3.30) — (b) book cites a proof it does not give [a-D29] | C15 | ★★ | 7 | antiderivative G with ∂G/∂x = F, fundamental theorem of calculus (P84 reminder), chain rule (P49), differentiation under the integral sign (primer) | the sign of the lower-limit term; ∂G/∂t is itself an integral of ∂F/∂t; the three-strip picture needs the corner piece to be second order | notebook · `reynolds_transport_cv` |
| D22 | Reynolds transport theorem (3.31) → (3.32) → (3.33) → (3.34) → (3.35) [a-D30] | C15 | ★★★ | 12 | definition (3.31) (N49), Taylor in time (gloss), orders of smallness (P68), signed swept volume (primer), mean-value theorem (P85), 1-D reduction to (3.30) (N46), sympy check cell | ΔV is **signed** (b·n < 0 subtracts, one formula covers both); $\int_{\Delta V}\Delta t\,\partial F/\partial t\,dV$ is $O(\Delta t^2)$ because ΔV = O(Δt); F in the sliver replaced by its surface value (error O(Δt)); d/dt passes inside the integral only when b = 0 | notebook · `reynolds_transport_cv` |
| D23 | RTT with F = 1, b = u, δV → 0 gives (3.14) — Exercise 3.28, (b) book never writes it out [a-D31] | C15 | ★★ | 5 | (3.35), divergence theorem (ch02 C14 reminder), mean-value theorem (P85) | the material volume has b = u (not b = 0); dividing by δV before the limit | notebook · `reynolds_transport_cv` |
| D24 | RTT reduces to (3.5) for a small material volume, exposing the $F\nabla\cdot\mathbf u$ term — Exercise 3.30, (b) book never writes it out [a-D32] | C15 | ★★ | 7 | (3.35) with b = u, divergence theorem, $\nabla\cdot(F\mathbf u)=\mathbf u\cdot\nabla F+F\nabla\cdot\mathbf u$ (gloss), (3.14) (C09) | the reduction needs (3.14) too: $(1/\delta V)\,d(F\delta V)/dt=DF/Dt+F\nabla\cdot\mathbf u$; the book's wording hides the extra term (Ch. 4's continuity trick) | notebook |

Counts: ★ 6 (D01, D03, D07, D08, D10, D18) · ★★ 17 (D02, D04, D05, D06, D09, D11, D12, D13, D14, D15, D16, D17, D19,
D20, D21, D23, D24) · ★★★ 1 (D22, shown in `reynolds_transport_cv` with a sympy check cell that re-runs the
construction: builds the swept shell of a growing sphere symbolically, shows the $\Delta V\cdot\Delta t$ term is
$O(\Delta t^2)$, then checks (3.35) against d/dt of the closed-form integral for F = t·x² on a growing sphere).
Written out under rule (b): D01, D03, D04, D05, D06, D10, D11, D13, D14, D16, D17, D18, D19, D20, D21, D23, D24 (the
book states the result, defers it to an exercise, or skips the algebra); under rule (a) only: D02, D07, D08, D09, D12,
D15, D22.

## 4c. Derivations demoted to statements (first column is the A parent in bold, e.g. **C20** — never a bare ID or a D id: the parser reads a bare first-cell ID as an item and blanks its tier)
Results **given, not derived**, inside the named A block (a paragraph, the equation, a number). 9 rows (6 demoted + 3
recap rows of analysis §2b).

| A parent | Analysis §2b item | Result stated (Eq.) | Stated in (B item) | Why not written out |
|---|---|---|---|---|
| **C02** | a-D03 streamwise form (3.6) | $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\lvert\mathbf u\rvert\frac{\partial F}{\partial s}$ with $\partial F/\partial s=\mathbf e_u\cdot\nabla F$ | N12 | one move ($\mathbf u=\lvert\mathbf u\rvert\mathbf e_u$) said in a sentence; the directional derivative is ch02 P75; the ⚠️ missing-F typo is the lesson |
| **C04** | a-D05 Ex. 3.1 streamline $y=x\tan\omega t'$ | $dy/dx=v/u=\tan\omega t'$, integrated through the origin | N18 | the book writes it; it is (3.7) (D03) with numbers, shown as the first line of C04's worked example |
| **C06** | a-D10 S + ½R split (3.11)–(3.13) (RECAP) | $\partial u_i/\partial x_j=S_{ij}+\tfrac12R_{ij}$, R = G − Gᵀ = 2A | R01 | recap of ch02 D14; one reminder sentence and the `strain_rate_tensor + 0.5*rotation_tensor == G` assert |
| **C10** | a-D16 R ↔ ω (3.15), components (3.16) (RECAP) | $R_{ij}=-\varepsilon_{ijk}\omega_k$; $\omega_3=\partial u_2/\partial x_1-\partial u_1/\partial x_2$ etc. | R05, R06 | recap of ch02 D15 and D12; reminder with the ch02 functions |
| **C10** | a-D18 irrotational ⇒ $\mathbf u=\nabla\phi$ (3.17) needs a simply connected region | $\nabla\times\nabla\phi=0$ (ch02 D26 corollary); converse only in simply connected regions; the line vortex is the counter-example | N29 | the easy direction is ch02's corollary; the converse is Ch. 6's business — stated with the line-vortex numbers (Γ = 2πB ≠ 0, φ = Bθ multivalued) and a pointer |
| **C12** | a-D21 shear-flow principal rates ±γ/2 at 45° (RECAP) | eigenvalues ±γ/2 of $\begin{bmatrix}0&\gamma/2\\\gamma/2&0\end{bmatrix}$ at ±45° | N35 | recap of ch02 D17 / N62 (Ex. 2.4) with the ⚠️ factor-2 note (γ = 2S₁₂) |
| **C13** | a-D24 solid-body circulation (3.24); any circuit (Exercise 3.23) | $\Gamma=2\pi r^2\omega_0$; uniform ω ⇒ Γ = ω × area for any loop | N37 | the book writes the circle case; "any circuit" is one sentence of Stokes (R07) with an off-centre numeric check |
| **C13** | a-D25 line vortex: $\omega_z=0$ off the axis, Γ = 2πB (3.26), δ-core (3.27), $\Gamma_{ABCD}=0$ | $(1/r)\,d(r\cdot B/r)/dr=0$; $2\pi B$; $\lim 2B/r^2$; four legs cancel | N39, N40, N41 | the book writes each line; they are D17's formula and four-side sum evaluated for $u_\theta=B/r$ — shown as C13's worked number, not a new derivation |
| **C15** | a-D33 Ex. 3.2 cone via RTT | $dV/dt=\int_{A^*}\mathbf b\cdot\mathbf n\,dA=\tfrac23\pi h^2\dot r\tan\theta=\tfrac23\pi hr_o\dot r$ | N54 | the book writes it out; it is C15's worked number (direct vs RTT, `ch03.example_3_2`) with the gaps named in the paragraph: b·n = 0 on the base, side speed (z/h)ṙ, slant element dz/cos θ, "[?]" = z |

## 5. Interactive explainers (5–10 + backup)
Seven explainers on the A items where manipulation or motion teaches most; one backup. Each has the required live
**Explain** tab ("Explanation & interpretation", numbered sections computing every number on screen with the reader's
settings, ending in "Reading the current setting"), a synced **Code** tab, a **Derivation** tab for its D ids, and ≥ 2
more depth features. Physics mirrors scalar-callable `fluidpy` functions (parity rows). Colours: streamline teal, path
line orange, streak line rose; local term blue, advective term amber, total purple; strain teal, rotation orange;
stretching blue / compressing rose; volume term blue, surface term orange.
Not duplicated from ch02: `strain_vs_rotation_split` (three squares G/S/A) and `stokes_circulation_loop` (curl image +
loop) are *linked from* E4/E5 and E6 rather than rebuilt.

### E1 · flow_lines_unsteady
- **A:** C03, C04 (also shows N13 steady ⇒ coincide, N14 (3.7) slope, N16 path line, N17 streak line, N18 Ex. 3.1 as
  the default preset, N19 Fig. 3.7)
- **Confusion removed:** "streamlines are where the fluid goes" (only in steady flow); the dye photograph (streak line),
  the particle track (path line) and the instantaneous arrow pattern (streamline) are three different curves.
- **Why interactive:** the three lines differ *because of time*; only a clock that runs — streamline pattern redrawn
  each instant, a particle leaving its trail, dye leaving the port — shows why they separate, and the "steady" toggle
  that makes them collapse onto one curve is the proof by experiment.
- **Stage:** (1) the flow plane: faint current streamline pattern (teal), a dye port emitting particles every Δt (rose
  streak), one tagged particle with its trail (orange path), velocity arrow at the port; (2) the three curves drawn at
  the drawing instant t′ with the closed forms of Ex. 3.1 as ghosts and the common tangent at the origin; (3)
  (hidePortrait) u(t), v(t) at the port over one period with the current time marker.
- **Controls:** time t (transport) · ω (s⁻¹) · amplitude ξ_o (m) · mean current U₀ (m/s, 0 = Ex. 3.1) · mode: Ex. 3.1 /
  Ex. 3.1 + mean current / steady (ω → 0) / rotating strain field.
- **Equations shown (live):** (3.7) $dx/u=dy/v=dz/w$; (3.8) $d\mathbf r/dt=\mathbf u(\mathbf r,t)$; Ex. 3.1 field
  $u=\omega\xi_o\cos\omega t,\ v=\omega\xi_o\sin\omega t$; streamline $y=x\tan(\omega t')$; path line
  $(x+\xi_o\sin\omega t')^2+(y-\xi_o\cos\omega t')^2=\xi_o^2$; streak line
  $(x-\xi_o\sin\omega t')^2+(y+\xi_o\cos\omega t')^2=\xi_o^2$.
- **Mirrors:** `core.K.streamline`, `core.K.pathline`, `core.K.streakline`, `ch03.example_3_1`, `ch03.unsteady_flow_preset` (§8).
- **Derivations:** D03 (step "ds ∥ u" freezes the clock and draws ds on the streamline), D04 (each step moves the
  particle; the elimination step draws the circle), D05 (the release-time step lights particles by their $t_o$).
- **Depth features:** Explain tab (field at the port now → streamline slope tan ωt′ → path-line centre and radius →
  streak-line centre → distance between the curves → reading the current setting), synced Code tab, + linked views (3),
  transport (play/step/scrub, loop over a period, end-of-run card "the particle came back; the dye did not"), presets
  (Ex. 3.1 at t′ = 0, π/4ω, π/2ω; steady), modes (four fields), status verdict ("🟰 steady: the three lines coincide" /
  "↔ unsteady: three different curves"), inspector (click a streak particle: its release time and path).
- **Follows reference:** `angular_frequency_explorer_1.html` (linked views on one clock, ghost reference, modes,
  end-of-run summary).
- **Aha:** freeze the clock and you get streamlines; follow one particle and you get a path line; photograph the dye
  and you get a streak line — and they are the same curve only if nothing changes in time.

### E2 · material_derivative_probe
- **A:** C01, C02 (also shows N07 labels, N08 (3.1), N09 (3.3), N10 local vs advective (warm-advection preset), N11
  (3.5), N12 streamwise form)
- **Confusion removed:** "∂T/∂t is how fast the air warms" (it is what a fixed station sees; the air parcel's rate is
  DT/Dt; the difference is advection $\mathbf u\cdot\nabla T$).
- **Why interactive:** the idea *is* two observers of one field. A fixed probe and a drifting float read the same moving
  temperature pattern on one clock; the reader drags the wind and the gradient and watches the term bars rebalance —
  including the case where the station warms while the float stays at constant T.
- **Stage:** (1) a temperature map (heatmap, moving front) with a fixed probe (square) and a float (circle) carried by
  u, both with their readings; (2) T(t) for the probe (blue) and the float (purple) with the slope triangles
  ∂T/∂t and DT/Dt; (3) term bars: local ∂T/∂t (blue) + advective u·∇T (amber) = DT/Dt (purple), and the float's
  measured dT/dt as a check marker.
- **Controls:** wind speed and direction (m/s, °) · gradient magnitude (K per 100 km) · heating rate of the whole field
  (K/h) · time (transport) · mode: thermal front (atmosphere) / Lagrangian stretching map $x=Xe^{\alpha t}$ (D01's case).
- **Equations shown (live):** (3.2) $F[\mathbf r(t;\mathbf r_o,t_o),t]=F(\mathbf x,t)$ at $\mathbf x=\mathbf r$; (3.3)
  $\frac{d}{dt}F[\mathbf r,t]=\frac{\partial F}{\partial r_i}\frac{dr_i}{dt}+\frac{\partial F}{\partial t}$; (3.4)/(3.5)
  $\frac{DF}{Dt}\equiv\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F=\frac{\partial F}{\partial t}+u_i\frac{\partial F}{\partial x_i}$;
  (3.6) $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\lvert\mathbf u\rvert\frac{\partial F}{\partial s}$.
- **Mirrors:** `core.K.material_derivative_terms`, `core.K.material_derivative`, `ch03.lagrangian_map_example`,
  `ch03.thermal_front` (§8).
- **Derivations:** D01 (mode "stretching map": each step computes u = αx at the float), D02 (the chain-rule step lights
  the float's velocity components; the "every point" step jumps the float to a new start).
- **Depth features:** Explain tab (∇T from the settings → u·∇T term by term → ∂T/∂t → DT/Dt boxed → the float's
  measured rate → "a station here sees +0.36 K/h of warm advection" → reading the current setting), synced Code tab, +
  linked views (3), term bars, transport, presets (pure advection DT/Dt = 0; pure heating u = 0; wind along the front
  u ⊥ ∇T; warm advection over a station), status verdict ("🌡️ warm advection" / "❄️ cold advection" / "⏸ no advection:
  u ⊥ ∇T").
- **Follows reference:** `forced_damped_vibrations.html` (the explanation panel with boxed numbers and a
  regime-dependent reading) with `fid_formula_lab.html`'s term bars.
- **Aha:** a thermometer riding with the air can read a constant temperature while the station it passes warms — the
  difference is the advective term.

### E3 · galilean_frames_cylinder
- **A:** C05 (also shows N02 steady vs unsteady, N04 Fig. 3.2, N20, N21 the transformation, N22 linear vs quadratic,
  N23 the frame-dependent split, N24 analogies as preset captions)
- **Confusion removed:** "the flow past a moving body is unsteady, so the acceleration is different" — steadiness and
  the split into local and advective parts depend on the observer; the particle's acceleration does not.
- **Why interactive:** one slider — the observer's velocity — morphs the streamline picture continuously from "steady
  flow around a fixed cylinder" to "unsteady flow ahead of a towed cylinder", while the acceleration arrow at a tagged
  particle stays exactly the same and its two bars trade places; no static pair of figures shows the continuous trade.
- **Stage:** (1) the cylinder and the flow in the chosen frame (instantaneous streamlines, tracer particles, one tagged
  particle with its acceleration arrow); (2) term bars at the tagged particle: local (blue) + advective (amber) = total
  (purple), x and y components; (3) (hidePortrait) the lab-frame and body-frame velocity vectors at the particle with U
  drawn between them (u = U + u′).
- **Controls:** observer velocity $U_{frame}$ from 0 (fluid frame) to U (body frame) · free-stream U (m/s) · cylinder
  radius a (m) · click to tag a particle · time (transport).
- **Equations shown (live):** Galilean transformation $\mathbf u(\mathbf x,t)=\mathbf U+\mathbf u'(\mathbf x',t')$,
  $t=t'$, $\mathbf x=\mathbf x'+\mathbf Ut+\mathbf x'_o$; (3.9)
  $\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=\frac{\partial\mathbf u'}{\partial t'}+(\mathbf u'\cdot\nabla')\mathbf u'$;
  cylinder flow (uniform stream + doublet, Ch. 6 form).
- **Mirrors:** `ch03.cylinder_flow`, `core.K.galilean_transform`, `core.K.acceleration`,
  `ch03.frame_acceleration_terms` (§8).
- **Derivations:** D06 (the chain-rule-in-a-moving-frame step animates x′ = x − Ut for one probe; the cancellation step
  lights the two $\mathbf U\cdot\nabla'$ bars that cancel).
- **Depth features:** Explain tab (u′ at the particle → local term → advective term → sum in each frame → the
  difference 0 → reading the current setting), synced Code tab, + linked views (3), term bars, presets (body frame:
  steady, only advective; fluid frame: both; halfway), status verdict ("🟰 steady in this frame" / "⏱ unsteady in this
  frame — same acceleration"), inspector (click any point: both frames' terms with numbers).
- **Follows reference:** `amplitude_phase_second_order_II_3.html` (windows linked by one state, numbered live
  derivation in the explanation).
- **Aha:** changing the observer moves acceleration from one bar to the other but never changes the total.

### E4 · fluid_element_deformation
- **A:** C06, C07, C08, C09 (also shows R01–R03 the S + ½R split, N27 rigid motion ⇒ S = 0 preset, N25 stress pointer)
- **Confusion removed:** "S_ij is just a symmetric matrix" — each entry is a *measurable rate*: how fast a material
  segment stretches (diagonal), how fast a right angle closes (off-diagonal, halved), how fast volume grows (trace).
- **Why interactive:** the reader drags G, watches a tracked element deform and reads the *measured* rates (from the
  moving tracers) next to the *formula* rates $\mathbf n\cdot\mathbf S\cdot\mathbf n$, $\mathbf n_1\cdot\mathbf S\cdot\mathbf n_2$,
  tr G; rotating the probe direction n and seeing the stretching rate trace out its rose curve is the tensor made
  visible.
- **Stage:** (1) the element: a square and a ring of tracer neighbours advected by $d\mathbf u=\mathbf G\cdot d\mathbf x$,
  with du arrows on the ring and a draggable probe segment n (and its perpendicular); (2) measured vs formula: stretching
  rate (1/ℓ)dℓ/dt vs n·S·n over all angles (rose curve with the probe's dot), corner-angle closing rate vs 2S₁₂; (3)
  (hidePortrait) area/volume ratio vs t with the ghost $e^{t\,\mathrm{tr}G}$.
- **Controls:** G₁₁, G₁₂, G₂₁, G₂₂ (s⁻¹) · probe angle θ · time (transport) · presets (pure stretch, simple shear,
  solid-body rotation (S = 0), expansion, the C06 worked G) · optional: dt of the measurement (shows the first-order
  error).
- **Equations shown (live):** (3.10) $du_i=(\partial u_i/\partial x_j)dx_j$; (3.11)
  $\partial u_i/\partial x_j=S_{ij}+\tfrac12R_{ij}$; (3.12) $S_{ij}=\tfrac12(\partial u_i/\partial x_j+\partial u_j/\partial x_i)$;
  linear strain rate $\frac1{\delta x_1}\frac{D(\delta x_1)}{Dt}=\frac{\partial u_1}{\partial x_1}$; shear strain rate
  $\tfrac12\frac{D(\alpha+\beta)}{Dt}=S_{12}$; (3.14) $\frac1{\delta V}\frac{D(\delta V)}{Dt}=\frac{\partial u_i}{\partial x_i}=S_{ii}$.
- **Mirrors:** `core.K.relative_velocity`, `core.K.linear_strain_rate`, `core.K.shear_strain_rate`,
  `core.K.volumetric_strain_rate`, `core.K.material_volume_ratio`, `core.K.linear_flow_map`,
  `core.K.measured_strain_rates` (§8).
- **Derivations:** D07 (the Taylor step shows the exact du vs G·dx gap shrinking with ring radius), D08 (sets the probe
  along x₁ and animates A′B′), D09 (sets simple shear and draws α, β), D11 (sets the expansion preset; product-rule
  step lights the three sides).
- **Depth features:** Explain tab (G → S and ½R entry by entry → the probe's n·S·n → measured rate over the last dt →
  corner closing rate → tr G and the area ratio → reading the current setting), synced Code tab, + linked views (3),
  transport, presets, inspector (click a tracer: its du = G·dx arithmetic), status verdict ("🔄 rigid motion: S = 0" /
  "🎈 expanding: tr G > 0" / "💧 volume-preserving").
- **Follows reference:** `forward_noising_lab.html` (click a point to see its exact arithmetic) with
  `forced_damped_vibrations.html`'s explanation panel.
- **Aha:** every entry of S is something you can measure with a ruler and a protractor on a moving element: stretch per
  length, half the closing rate of a right angle, and (summed) the growth rate of volume.

### E5 · spin_and_principal_axes
- **A:** C10, C11, C12 (also shows R05 (3.15), N28 rotating observer ω′ = ω − 2Ω, N30 rigid rotation, N31–N32 the
  principal frame, N34 parallel shear flow, N35 the two elements of Fig. 3.14)
- **Confusion removed:** "a fluid moving in a straight line cannot rotate" and "vorticity is the angular velocity" —
  in a parallel shear flow each element spins at −γ/2 (half the vorticity), individual lines turn at different rates,
  and only their perpendicular-pair average is the spin; the same element also stretches along 45° axes.
- **Why interactive:** the key claims are "for *any* pair" and "for *any* direction": dragging the pair angle θ and
  watching the two lines' rates change while their average stays flat, and dragging the probe dx around the ring while
  its arrow splits into a strain part and a rotation part, turn the quantifiers into something seen.
- **Stage:** (1) the element with two perpendicular material lines (draggable pair angle θ), a small paddle wheel at
  the centre, the principal axes of S (blue stretch / rose squeeze), and a circle deforming into its ellipse; (2)
  θ̇(θ) curve for a single line with the pair's two dots and their flat average at ½ω₃; (3) the (3.19) split: one
  relative-velocity arrow du (purple) = S·dx (teal) + ½ω × dx (orange), for the probe dx on the ring.
- **Controls:** G entries or presets (parallel shear γ, solid-body rotation, pure strain, general) · pair angle θ ·
  probe direction on the ring · observer rotation Ω (s⁻¹, optional) · time (transport).
- **Equations shown (live):** $\tfrac12\frac{D(-\alpha+\beta)}{Dt}=\tfrac12\omega_3=R_{21}/2$; (3.15)
  $R_{ij}=-\varepsilon_{ijk}\omega_k$; $\dot\theta=-\gamma\sin^2\theta$ (parallel shear); (3.19)
  $du_i=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i$; (3.20) $d\bar{\mathbf u}=\bar{\mathbf S}\cdot d\bar{\mathbf x}$;
  (3.21) $d\bar u_\alpha=\bar S_{\alpha\alpha}d\bar x_\alpha$; $\omega'_z=\omega_z-2\Omega$.
- **Mirrors:** `core.K.element_rotation_rate`, `core.K.material_line_rotation_rate`, `core.K.relative_velocity_split`,
  `core.K.principal_strain_rates`, `core.K.strain_ellipse_axes`, `core.K.vorticity_in_rotating_frame`,
  `ch03.parallel_shear_kinematics`.
- **Derivations:** D12 (the limit step shrinks dt and shows α, β), D13 (the observer-rotation step sets Ω = ω₃/2 and the wheel stops), D14 (the pair-average step sweeps θ over 0–π), D15
  (the ε-swap step flips the orange arrow's sign on screen), D16 (the eigenframe step rotates the view onto the
  principal axes).
- **Depth features:** Explain tab (ω₃ from G → the two lines' rates at θ → their average = ½ω₃ → S's eigenvalues and
  45° axes → the (3.19) split for the probe → ellipse semi-axes after t → observer rotation → reading the current
  setting), synced Code tab, + linked views (3), presets (shear γ = 1, solid body, pure strain, "rotate with the
  element" Ω = ω/2 ⇒ ω′ = 0), transport, status verdict ("🌀 spins and strains (shear)" / "🔄 spins only" / "↔ strains
  only"), inspector (click a ring point: du split with numbers).
- **Follows reference:** `amplitude_phase_second_order_II_3.html` (system + graphs linked by one state; numbered live
  derivation in the explanation).
- **Aha:** vorticity is twice the spin of the *average* line; in a shear flow the element spins at −γ/2 while it
  stretches along 45° — rotation and deformation are two separate, additive motions.

### E6 · vortex_paddle_wheels
- **A:** C13, C14 (also shows N36 solid body (3.22), N37 (3.24), N38 line vortex (3.25), N39 (3.26), N40 the δ core,
  N41 Γ_ABCD = 0, N42 real vortices, N43 Gaussian (3.29), N44 the 1.1209σ maximum, R07 circulation, N29 irrotational
  flow's simply-connected caveat)
- **Confusion removed:** "fluid going round in circles is rotating" — in the irrotational vortex elements orbit without
  spinning; in solid-body rotation they spin at the orbit rate; real vortices are a spinning core wrapped in an
  irrotational skirt.
- **Why interactive:** the paddle-wheel test is a motion: wheels carried round by each profile either turn with the
  flow or keep pointing the same way; dragging σ and a circulation loop across the core edge and watching Γ(r) climb
  and then saturate makes (3.23)–(3.29) one picture.
- **Stage:** (1) the vortex plane: tracers, several paddle wheels at different r (spinning at ½ω_z), a draggable loop
  (circle or annular sector ABCD) with its Γ; (2) profiles $u_\theta(r)$ (teal) and $\omega_z(r)$ (orange) vs r/σ with
  the core edge σ and the maximum marker; (3) (hidePortrait) Γ(r) with the ghost 2πB and mean vorticity in a disc on
  log axes (slope −2 for the line vortex).
- **Controls:** mode: solid body / line vortex / Rankine / Gaussian · Γ (m²/s) · σ (m) · loop radius and centre (drag),
  loop shape (circle / sector) · time (transport).
- **Equations shown (live):** (3.22) $u_r=0,\ u_\theta=\omega_0r$; (3.23)
  $\omega_z=\frac1r\frac{\partial}{\partial r}(ru_\theta)-\frac1r\frac{\partial u_r}{\partial\theta}$; (3.24)
  $\Gamma=2\pi ru_\theta=2\pi r^2\omega_0$; (3.25) $u_\theta=B/r$; (3.26) $\Gamma=2\pi B$; (3.27)
  $[\omega_z]_{r\to0}=\lim 2B/r^2$; (3.28) Rankine $u_\theta=(\Gamma/2\pi\sigma^2)r$ ($r\le\sigma$), $\Gamma/2\pi r$ ($r>\sigma$);
  (3.29) $u_\theta=\frac{\Gamma}{2\pi r}(1-e^{-r^2/\sigma^2})$, $\omega_z=\frac{\Gamma}{\pi\sigma^2}e^{-r^2/\sigma^2}$;
  $1+2r^2/\sigma^2=e^{r^2/\sigma^2}$.
- **Mirrors:** `core.X.solid_body_rotation`, `line_vortex`, `rankine_vortex`, `gaussian_vortex`,
  `gaussian_vortex_max_radius`, `polar_vorticity_z`, `circulation_circle`, `mean_vorticity_in_disc`,
  `ch03.annular_sector_circulation`, `core.X.vortex_profile` (§8).
- **Derivations:** D17 (each leg of the polar sector lights in turn), D18 (sets Rankine; the continuity step puts the
  loop at r = σ), D19 (sets Gaussian; the integral step sweeps the loop outward), D20 (sets the maximum marker).
- **Depth features:** Explain tab (u_θ at the loop → Γ = 2πr u_θ → ω_z by (3.23) → the wheel's spin ½ω_z vs its orbit
  rate u_θ/r → where the maximum is → reading the current setting), synced Code tab, + linked views (3), modes (four
  vortices), presets (bathtub, tornado, tropical-cyclone-scale Gaussian — our numbers), transport, status verdict ("🔄
  the wheel spins: inside the core" / "🧭 the wheel keeps its heading: irrotational" / "⚠️ loop encloses the singular
  core"), inspector (click a point: u_θ, ω_z arithmetic), a small table of real vortices with the current row lit.
- **Follows reference:** `angular_frequency_explorer_1.html` (modes, linked views on one clock, "Right now" notes with a
  highlighted table).
- **Aha:** going round in circles is not spinning: the wheel in the irrotational vortex orbits yet always points the same
  way, and a real vortex spins only in its core.

### E7 · reynolds_transport_cv
- **A:** C15 (also shows N45, N46 Leibniz (3.30) as the 1-D mode, N47 Fig. 3.17 strips, N48 control volume, N49–N52
  the (3.31)–(3.34) steps, N53 F = 1 and b = u interpretations, N54 Ex. 3.2 cone as a 3-D mode, N55 Fig. 3.18)
- **Confusion removed:** "d/dt of an integral is the integral of ∂/∂t" (only for a fixed volume); and the sign of the
  boundary term (b·n > 0 sweeps F in, b·n < 0 sweeps it out — one signed formula).
- **Why interactive:** the theorem is a budget over time; watching a boundary move, the swept band coloured by the sign
  of b·n, and bars "volume term + surface term" landing exactly on the measured d/dt ∫F — then shrinking Δt until the
  dropped second-order term vanishes — is the derivation as an experiment.
- **Stage:** (1) the control volume (1-D interval a(t), b(t) over F(x, t) in Leibniz mode; a 2-D deforming ellipse over
  a heatmap of F in RTT mode; a 3-D cone in Ex. 3.2 mode) with b arrows and the swept band (blue b·n > 0, rose b·n < 0);
  (2) waterfall bars: volume term (blue) + surface term (orange) = total (purple) vs the finite-difference d/dt ∫F
  (marker); (3) (hidePortrait) the dropped term $\int_{\Delta V}\Delta t\,\partial F/\partial t$ vs Δt on log–log axes
  (slope 2).
- **Controls:** mode: Leibniz 1-D / RTT 2-D / Ex. 3.2 cone · boundary speeds (ȧ, ḃ or growth rates of the ellipse axes,
  translation velocity) · F choice (uniform, ramp, warming field) · Δt (log) · time (transport).
- **Equations shown (live):** (3.30) $\frac{d}{dt}\int_{a(t)}^{b(t)}F\,dx=\int_a^b\frac{\partial F}{\partial t}dx+\frac{db}{dt}F(b,t)-\frac{da}{dt}F(a,t)$;
  (3.31) $\frac{d}{dt}\int_{V^*}F\,dV=\lim_{\Delta t\to0}\frac1{\Delta t}\{\int_{V^*(t+\Delta t)}F(t+\Delta t)dV-\int_{V^*(t)}F(t)dV\}$;
  (3.34) $\int_{\Delta V}F\,dV\cong\int_{A^*}F(\mathbf b\Delta t\cdot\mathbf n)dA$; (3.35)
  $\frac{d}{dt}\int_{V^*}F\,dV=\int_{V^*}\frac{\partial F}{\partial t}dV+\int_{A^*}F\,\mathbf b\cdot\mathbf n\,dA$; Ex. 3.2
  $dV/dt=\tfrac23\pi hr_o\dot r$.
- **Mirrors:** `core.R.leibniz_terms`, `core.R.reynolds_transport`, `core.R.rtt_check`, `core.R.volume_integral_rate_fd`,
  `ch03.example_3_2`, `core.R.rtt_ellipse_2d` (§8).
- **Derivations:** D21 (Leibniz mode; the three-strip step shades the strips), D22 ★★★ (12 steps; the split step
  colours ΔV, the Taylor step shows the four (3.32) terms as bars, the "drop second order" step scrubs Δt on the
  log–log view, the swept-sliver step draws $\mathbf b\Delta t\cdot\mathbf n$ on the boundary, the signed-volume step
  sets a retreating side), D23 (F = 1, b = u: the bars reduce to the volume growth rate).
- **Depth features:** Explain tab (V*, A* from the settings → the volume term → the surface term piece by piece with
  signs → total vs the measured rate → the error at this Δt → reading the current setting), synced Code tab, + linked
  views (3), term bars (waterfall), transport, modes (1-D / 2-D / cone), presets (fixed volume b = 0; rigid translation
  in a uniform F — the two terms cancel; growing balloon; Ex. 3.2), status verdict ("📦 fixed volume: d/dt passes
  inside" / "🎈 moving boundary: the surface term matters").
- **Follows reference:** `fid_formula_lab.html` (terms that add up to a total, stage chips, code synced to the active
  line) with `forced_damped_vibrations.html`'s explanation panel.
- **Aha:** the rate of change of what is inside a moving volume = what changes in place + what the moving walls sweep
  in — and the walls' part is signed, so one formula covers growing and shrinking.

### B1 · material_volume_divergence (backup)
- **A:** C09, C15 (also shows N53 F = 1, b = u)
- **Confusion removed:** "∇·u is an abstract derivative" — it is the fractional growth rate of a blob of fluid.
- **Why interactive:** a blob of tracers in a compressible linear (or nonlinear) flow grows or shrinks while its area
  is plotted against $e^{t\,\mathrm{tr}G}$ and against the RTT surface integral $\oint\mathbf u\cdot\mathbf n\,dA$; the
  blob-size slider shows the local limit (3.14) emerging.
- **Stage:** blob + area-vs-time with ghost + bars (surface flux vs ∫∇·u). **Controls:** G entries or a nonlinear field
  preset, blob radius, time. **Equations:** (3.14) $\frac1{\delta V}\frac{D(\delta V)}{Dt}=S_{ii}$, (3.35) with F = 1,
  $\det e^{Gt}=e^{t\,\mathrm{tr}G}$. **Mirrors:** `core.K.material_volume_ratio`, `core.R.material_volume_rate`.
- **Derivations:** D11, D23. **Depth features:** explain, code, + linked views, transport, status. **Follows reference:**
  `overfitting_curves.html` (minimal two-slider figure with a verdict line). **Aha:** divergence is how fast a blob of
  fluid swells, per unit of its own volume.

## 6. Python animations and interactive figures (A ID → what, why, player/figure kind)
Animations (`animate` + `show_animation`, ≤ 120 frames, dpi 80, FAST halves the frames):
- **C04** — Ex. 3.1 over one period: dye leaving the port (streak line), one particle with its trail (path line), the
  instantaneous streamline turning; why: the three curves separate only in motion; `player="video"`.
- **C07/C08** (with C12) — a square and a 45° square in the parallel shear flow (our Fig. 3.14 animated), side lengths and
  corner angles printed per frame; why: stretching without shear vs shear without stretching must be watched;
  `player="frames"` (stop at each dt).
- **C12** — a circle of tracers becoming an ellipse with the eigen-axes drawn; after long times the axes drift from
  45° (finite-time caveat); why: the first-order claim and its limit are both visible; `player="video"`.
- **C13** — paddle wheels carried by solid-body rotation vs by the line vortex (Figs. 3.15/3.16 animated); why: the
  paddle-wheel test is a motion; `player="video"`.
- **C15** — Leibniz strips with a(t), b(t) moving and F(x, t) changing (Fig. 3.17 animated); why: the three strips
  appear as the ends move; `player="frames"`.

Interactive figures (`slider_figure`, plotly, precomputed; 3-D plotly where noted):
- **C02** — slider over time: the three terms ∂T/∂t, u·∇T, DT/Dt along a float's path through a moving front; why: the
  balance shifts as the float crosses the front.
- **C03** — slider over the drawing time t′: the Ex. 3.1 streamline pattern at that instant with the path line and
  streak line fixed; why: streamlines change with t′ while the others are history.
- **C05** — slider over the observer velocity: streamlines of the cylinder flow and the (local, advective) pair at a
  point; why: the continuous trade between the two terms.
- **C10** — slider over γ (or G preset): θ̇(θ) curve with the flat pair average; why: the "any pair" claim at a glance.
- **C14** — slider over σ: Rankine and Gaussian $u_\theta(r)$, $\omega_z(r)$ with the 1.1209σ marker; why: the core size
  sets where the peak is.
- **C15** — slider over Δt: the four (3.32) terms for a growing sphere on log–log axes; why: orders of smallness made
  visible.
- **N05 (→ C01)** — plotly 3-D: a point P with its Cartesian, cylindrical and spherical unit vectors, dropdown per
  system; why: the local basis rotates with P.
- **N15 (→ C03)** — plotly 3-D stream tube seeded on a circle; why: "no flux through the wall" is a 3-D statement.

## 7. From-scratch moments
Every book section with a computable A item has at least one (§3.1 has no A item; its N03 section average gets a
two-line `np.trapezoid` vs `quad` check instead):
- §3.2 **C02** — hand-written second-order central stencils for ∂F/∂t and u·∇F at a point vs
  `core.K.material_derivative_terms`; plus DF/Dt vs d/dt of F sampled along a `pathline` (the independent route):
  `assert np.allclose(mine, lib)`.
- §3.3 **C04** — a hand-written RK4 loop for the path line of Ex. 3.1 vs `core.K.pathline` and vs the closed-form circle.
- §3.4 **C07** — track a 1 cm material segment's two end points with `linear_flow_map`, compute (1/ℓ)dℓ/dt by hand vs
  `core.K.linear_strain_rate` (n·S·n); **C11** — a triple loop over ε_ijk for ½ε_ijk ω_k dx_j vs
  `core.K.relative_velocity_split`.
- §3.5 **C13** — the four-leg circulation around a small polar sector divided by its area vs
  `core.X.polar_vorticity_z` (and vs the Cartesian curl).
- §3.6 **C15** — midpoint-rule volume and surface sums for a growing sphere vs `core.R.reynolds_transport` and vs the
  finite difference of the volume integral.

## 8. Notes for the implementer
Functions the A items' figures and the explainers need that analysis §4 did not plan (all scalar-callable for parity
rows, SI units, docstrings citing § and Eq.):
- `ch03.unsteady_flow_preset(name, x, y, t, **p) -> (u, v)` — E1's fields: "ex31" (Ex. 3.1), "ex31_current" (+ U₀),
  "steady" (ω → 0 limit), "rotating_strain" (strain axes turning at rate Ω). Reused by the C03/C04 figures.
- `ch03.thermal_front(x, y, t, grad_K_per_m, heating_K_per_s, front_speed=0.0) -> T` and its exact gradient/time
  derivative — E2's field and C02's worked number (parity: DT/Dt from `material_derivative_terms` = analytic).
- `ch03.frame_acceleration_terms(x, y, U, a, U_frame, t) -> dict(local, advective, total)` — E3's bars in any
  observer frame (wraps `galilean_transform` + `acceleration` on `cylinder_flow`); parity: total independent of
  `U_frame` to 1e-8.
- `core.K.measured_strain_rates(G, n, dt) -> dict(stretch, closing, area)` — E4's "measured" readouts from tracked
  segments (`linear_flow_map`), converging to n·S·n, 2n₁·S·n₂, tr G as dt → 0 (tests: order 1 in dt).
- `core.X.vortex_profile(kind, r, Gamma, sigma) -> (u_theta, omega_z)` — one dispatcher over solid/line/Rankine/Gaussian
  for E6 and the C14 slider figure (solid body parametrised by Γ at r = σ for comparability, stated in the docstring).
- `core.R.rtt_ellipse_2d(F, dFdt, a_fn, b_fn, center_fn, t) -> (volume_term, surface_term, total)` — E7's 2-D
  deforming control volume (ellipse with axes a(t), b(t) translating at c′(t)), plus `core.R.swept_terms_sphere(R, Rdot,
  F, dFdt, t, dt)` returning the four (3.32) terms for the C15 slider figure.
- `ch03.leibniz_example(t, case) -> dict(F, a, b, adot, bdot)` — the 1-D cases shared by N46's figure, D21's check and
  E7's Leibniz mode.
- Parity rows must include Ex. 3.1 closed forms vs `streakline` (E1), the 1.1209064228σ radius (E6), Ex. 3.2 both
  routes (E7). Book numbers stay in `tests/book_values_ch03.json`.
- Budget: streak lines vectorised (one ODE for all release times); animations use FAST = 30 frames; RTT quadratures
  cached per shape; the whole notebook < 5 min on Colab CPU.
- Promotion reminder (analysis §4): `velocity_gradient_preset`, `linear_flow_map`, `deform_square`,
  `material_line_angle` move from ch02 to `core.kinematics` with re-exports (ch03 is their second user).
