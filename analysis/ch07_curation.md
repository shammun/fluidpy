# Chapter 7 — Gravity Waves: curation
(from `analysis/ch07.md` — 236 inventory rows = 202 numbered rows (#1–#202) + 33 figure rows (F1–F33) + 1 exercise
row (X1); 162 equation labels (7.1)–(7.159) with 7.34, 7.35, 7.84 split a/b; 68 derivations in §2b; 18 candidate A
ideas in §3; `book.yaml` ch07, `policy.tier_a_full_treatment: [12, 18]`, `coverage: exhaustive`. Curated 2026-09-24,
concept-curator. Read: `knowledge/CUMULATIVE.md`, `concept_map.md`, `primers.md` (P01–P164), `notation.md`,
`viz_patterns.md`, `knowledge/ch06.md` §8, `analysis/ch06_curation.md` (format).)

Counts: A 16 · B 193 · C 27 · RECAP 20 · SKIP 1 · derivations written out 37 (★ 9 ★★ 24 ★★★ 4) · demoted to statements 24

Tier words (parsed by `tools/nbkit.py`): **CORE 16 · NOTE 199 · RECAP 20 · SKIP 1** = 236 rows; **DERIVATION 37**.
Depth: A 16 (all CORE) · B 193 (182 NOTE + 11 RECAP) · C 27 (17 NOTE: N01, N12, N18, N49, N56, N64, N85, N96, N98,
N108, N109, N147, N151, N156, N191, N192, N196 · 9 RECAP: R01, R07, R11, R12, R15, R16, R17, R19, R20 · 1 SKIP S01).
**Reconciliation with analysis §2:** 202 numbered rows `[#1]…[#202]` + 33 figure rows `[#F1]…[#F33]` + the exercise
row `[#X1]` = 236 rows; each tag appears exactly once in §2 below (checked by script). A = 16 of 236 rows (7 %).

**Decisions that shape this chapter.**
1. **Sixteen A items from the analyst's eighteen.** Three merges and one split:
   (a) the analyst's "pressure under waves" (7.30)–(7.31), (7.48), (7.52) is not an A item — it is a set of B items of
   **C04** (phase speed and the deep/shallow regimes), because every pressure result is the same cosh ratio read in the
   two limits, and the headline "shallow water is hydrostatic" is a limit of that ratio;
   (b) the analyst's "group velocity" and "energy flux" (7.69)–(7.71) meet in **C09** (energy flux F = E c_g is where
   c_g earns its meaning), while the energy itself (7.38)–(7.44) stays its own A item **C06** (equipartition, E = ½ρga²);
   (c) the analyst's "KdV and solitary wave" (7.87)–(7.88) is taught inside **C11** (hydraulic jump) as the other fate
   of a steepening shallow-water wave — the Ursell ratio aλ²/H³ decides between a jump and a permanent form — and the
   internal-wave energy (7.147)–(7.159) is taught inside **C16** (F = c_g E is the payoff of c ⟂ c_g);
   (d) the refraction rows (Figs. 7.8–7.9, §7.2) are taught in **C10** (kinematic wave theory), where the rays that
   explain them are derived; the §7.2 notebook section carries a forward pointer.
   §7.7 keeps two A items — **C13** (two deep fluids: ε√(gk), the vortex sheet) and **C14** (layer over deep fluid:
   barotropic/baroclinic modes, reduced gravity) — because the second is the vocabulary of every layered ocean model in
   Ch. 13 and deserves its own derivation chain and explainer.
2. **SEEN rows are RECAP** (20 rows, R01–R20, inventory order), reminded with the earlier chapter's function. Rows
   marked "SEEN + NEW" or "SEEN → NEW" whose new part is a derivation or a notation change are **B NOTE items** inside
   their A block (N13 the f = z − η normal, N15 the explicit kinematic form, N20 the ODE step, N23–N24 φ and u, w,
   N28 linearised path lines, N43 c = √(gH), N82 the jump, N83 the CV momentum of the jump, N99 complex notation, N107
   the vortex sheet, N128 g′, N145 the ω = N limit, N152 the energy equation, N157 the δ-function N²). **One CORE row
   carries a SEEN→NEW mark: C03 (7.28).** ch04 *used* ω² = gk tanh kH as a given test field (`linear_wave_surface`)
   and never derived it; the derivation is this chapter's new content (the analyst's §1 lists "the dispersion relation
   derived" as NEW), so the row is CORE and its recap sentence (ch04 C14 test field) sits inside C03.
3. **Derivations written out: 37** (★ 9 · ★★ 24 · ★★★ 4), covering 44 of the 68 analysis rows (merges: a-D17 + a-D18
   → D13; a-D23 + a-D24 → D15; a-D31 + a-D32 + a-D34 → D21; a-D52 + a-D53 → D31; a-D54 + a-D55 → D32; a-D64 + a-D66 →
   D37). Rule (b) ("the book never writes it out and the lesson needs it"): **D20** (packet envelope moves at c_g, cited
   to Phillips — ★★★, the heart of group velocity) and **D24** (Snell's law for refraction — the book gives words only).
   **24 demoted to statements** (§4c), each with a numeric or sympy check in code.
4. **Book slips taught in our own words** (analysis §9): T1 (7.66) prints ½Δω **x** — read ½Δω **t** (N65, D19; the
   printed form shown as a ghost whose envelope does not move); T2 "u from (7.28)" means (7.27) (N36, D14); T3 (7.105)
   prints e^{i(kz−ωt)} — read e^{i(kx−ωt)} (N118, D29, failing sympy residual for the printed form); T4 (7.98) ∂/dz
   (N111); T5 "y = 0" means z = 0 (N12); T6 "(7.88)" means (7.87) in the Ursell remark (N95); T8 Exercise 7.6 lake
   period vs our 31.93 min stays in the private JSON (N63); T9 T = 10 s ↔ 156 m, not 150 m (N39).
5. **Conventions stated where first used** (⚠️ callouts with numbers, §8): η is the **elevation** here but the level-set
   function in ch04 (R1; N13, D02); p gauge, p′ = p + ρgz in §7.2 but p − p̄(z) in §7.8 (R07, R17); ω is a frequency
   (vorticity in ch02–ch06); ζ is a particle excursion (§7.2, §7.6, §7.8) and an interface displacement (§7.7); θ is the
   local phase (7.72) and the angle of K (7.139); g′ with **ρ₂** (7.117) vs ch04's ρ₁ (N128); the internal-wave beam
   angle is θ **to the vertical** = K's angle to the horizontal (R10; N150); the printed (7.138), (7.145) assume k > 0
   (R9; code uses ∣k∣ and ∇_K ω); energies per unit horizontal area (surface, interface) vs per unit volume (internal)
   (R5; N166); g = 9.81 (`G_BOOK`) in chapter functions.
6. **Climate hooks, where real:** c = √(gH) and hydrostatic shallow water (C04) are the shallow-water model of Ch. 13
   (tides, tsunamis, Kelvin waves: H = 4 km → 198 m/s); swell energy flux F = E c_g ≈ 38 kW/m (C09); seiches and
   basin resonance (C08); rays and WKB (C10) for internal waves in N(z) and topographic waves; Stokes drift (C12) is
   ocean-surface transport, Langmuir cells and the Stokes–Coriolis force; reduced gravity c = √(g′H) ≈ 1 m/s vs
   22 m/s (C14) is the thermocline wave speed and the baroclinic Rossby radius of Ch. 13; ω = N cos θ, beams and
   F = c_g E (C15–C16) are internal tides over topography, lee waves and wave drag in the stratosphere.

## 1. Teaching order (A IDs grouped by book section, B/C IDs under each; one sentence each: "once you see X, Y follows")
A items in **bold**; B and C items listed where they are taught, depth in brackets. Equations written beside their
numbers so downstream agents have them.

**§7.1 Introduction**
- **C01 The sinusoidal travelling wave and its vocabulary** $\eta=a\cos[\frac{2\pi}{\lambda}(x-ct)]$ *(Eq. 7.1)* [#3]:
  once a crest is "the place where the phase is 2nπ", its speed $c=\omega/k=\lambda\nu$ *(Eq. 7.4)* (D01; N04, N05 [B])
  follows, and so do the k–ω form (N03 [B]), the plane wave $\eta=a\cos(\mathbf K\cdot\mathbf x-\omega t)$ *(Eq. 7.5)*
  with $K^2=k^2+l^2+m^2$, $\lambda=2\pi/K$, $\mathbf c=(\omega/K)\mathbf e_K$ (N06–N09 [B]), trace velocities larger than
  c (N10 [B]), the Doppler shift $\omega_0=\omega+\mathbf U\cdot\mathbf K$ *(Eq. 7.9)* (N11 [B]), Fourier superposition
  (N02 [B]) and Fig. 7.1 (N167 [B]); wave families named (N01 [C]).

**§7.2 Linear liquid-surface gravity waves**
- **C02 The linear free-surface problem** $(\partial\phi/\partial z)_{z=0}\cong\partial\eta/\partial t$ *(Eq. 7.18)* [#22]:
  once you see that a condition on the unknown surface z = η can be Taylor-moved to z = 0 at the cost of O(ka), the
  whole problem becomes Laplace $\nabla^2\phi=0$ (R01 [C], R02 [B]) in a fixed strip with a flat bottom (R03 [B]), the
  kinematic condition (R04 [B]; N13, N14, N15, N16 [B]; D02, D03) and the dynamic condition from linear Bernoulli
  $(\partial\phi/\partial t)_{z=0}\cong-g\eta$ *(Eq. 7.21)* (R05, R06 [B]; N17 [B]; D04); set-up and Fig. 7.2 (N12 [C],
  N168 [B]).
- **C03 The dispersion relation** $\omega=\sqrt{gk\tanh kH}$ *(Eq. 7.28)* [#34]: once the separable trial
  $\phi=f(z)\sin(kx-\omega t)$ *(Eq. 7.22)* turns Laplace into $f''-k^2f=0$ (N19, N20 [B]) and the bottom and kinematic
  conditions fix f (N21, N22 [B]), the potential $\phi=\frac{a\omega}{k}\frac{\cosh k(z+H)}{\sinh kH}\sin(kx-\omega t)$
  *(Eq. 7.26)* and velocities (N23, N24 [B]; D05) follow, and the dynamic condition (N25 [B]) leaves only one ω for each
  k (D06); the initial-shape remark (N18 [C]).
- **C04 Phase speed, dispersion and the deep/shallow regimes** $c=\sqrt{\frac gk\tanh kH}$ *(Eq. 7.29)* [#35]: once c
  grows with λ (dispersive, D07), the two limits of tanh (N37 [B], N173 [B]) give deep water $c=\sqrt{g/k}$ *(Eq. 7.45)*
  (N38 [B]; D08) and shallow water $c=\sqrt{gH}$ *(Eq. 7.49)* (N43 [B]; D09) with the ocean numbers (N39 [B]); the same
  cosh ratio gives the pressure $p'=\rho ga\frac{\cosh k(z+H)}{\cosh kH}\cos(kx-\omega t)$ *(Eq. 7.31)* (R07 [C], N26
  [B]), its deep decay (N42 [B]) and the hydrostatic shallow limit $p'=\rho g\eta$ *(Eq. 7.52)* (N46 [B]).
- **C05 Particle orbits and streamlines** — ellipses $\xi^2/A^2+\zeta^2/B^2=1$ *(Eq. 7.36)* [#42]: once the path-line
  equations (R08 [B]; N27, N28 [B]) are linearised about the mean position, the excursions (N29 [B]; D10) trace
  ellipses (D11) that become circles in deep water (N40, N41 [B]) and flat ellipses in shallow water (N44, N45 [B]);
  streamlines $\psi$ *(Eq. 7.37)* (N30 [B]); Figs. 7.3–7.5 (N169, N170, N171 [B]).
- **C06 Wave energy and equipartition** $E=\tfrac12\rho ga^2$ *(Eq. 7.42)* [#48]: once the kinetic energy integral
  (N31, N32 [B]; D12) and the column-swap potential energy (N33 [B], N172 [B]; D13) come out equal (N34 [B]), their sum
  is ½ρga², and the pressure work across a vertical line (N35 [B]; D14) is E times a speed that is not c (N36 [B]) — the
  cliff-hanger for C09.
  *(Refraction, N47 and N48 [B] with Figs. 7.8–7.9, is §7.2 content taught in C10.)*

**§7.3 Influence of surface tension**
- **C07 Capillary–gravity waves and the minimum phase speed** $c=\sqrt{(\frac gk+\frac{\sigma k}{\rho})\tanh kH}$
  *(Eq. 7.57)* [#69]: once the Laplace jump (R09 [B]) on a curved surface (N50, N51 [B]) adds $\frac\sigma\rho\eta_{xx}$
  to the dynamic condition (N52 [B]; D15), every g becomes g + σk²/ρ (N53 [B]), c(λ) has two branches and a minimum
  $c_{min}=(4g\sigma/\rho)^{1/4}$ *(Eq. 7.58)* (N54 [B]; D16) with the air–water numbers (N55 [B]), the pure-capillary
  branch (N57 [B]) and Fig. 7.10 (N176 [B]); prose and ripples named (N49, N56 [C]).

**§7.4 Standing waves**
- **C08 Standing waves and seiches** $\eta=2a\cos kx\cos\omega t$ [#75]: once a right- and a left-going wave (N58 [B])
  add into fixed nodes (D17), their stream function and velocity (N59, N61 [B]; Fig. 7.11 N177 [B]) show where walls may
  stand; walls at x = 0, L (N60 [B]) allow only $\lambda=2L/(n+1)$ *(Eq. 7.64)* (N62 [B]) and the lake's natural
  frequencies $\omega=\sqrt{\frac{\pi g(n+1)}{L}\tanh\frac{(n+1)\pi H}{L}}$ *(Eq. 7.65)* (N63 [B]; D18; Fig. 7.12 N178
  [B]).

**§7.5 Group velocity, energy flux, and dispersion**
- **C09 Group velocity** $c_g=d\omega/dk$ *(Eq. 7.67)* [#83]: once two nearby waves beat (N65 [B]; D19) and a packet
  (N66 [B]) is seen to carry its envelope at c_g (N67 [B]; D20 ★★★), the water-wave value
  $c_g=\frac c2[1+\frac{2kH}{\sinh 2kH}]$ *(Eq. 7.69)*, its limits c/2 and c (N68, N69 [B]) and the energy flux
  $F=Ec_g$ *(Eq. 7.71)* (N70 [B]; D21) follow, with the 3-D gradient form (N71 [B]), the stone in a pond (N72 [B]),
  crests running through a group (N73 [B]) and Figs. 7.13–7.16 (N179–N182 [B]); history named (N64 [C]).
- **C10 Kinematic wave theory: rays** $\partial\omega/\partial t+c_g\partial\omega/\partial x=0$ *(Eq. 7.79)* [#99]: once
  local k and ω are phase gradients (N74, N75 [B]), crests are conserved (N76 [B]) and wavenumbers ride at c_g (N77
  [B]; D22; Fig. 7.17 N183 [B]); in a slowly varying depth (N78, N79, N80 [B]) frequency is constant along rays (D23;
  Fig. 7.18 N184 [B]), so crests turn toward the shore (N47 [B]; D24 Snell) and wrap round islands (N48 [B]; Figs.
  7.8–7.9 N174, N175 [B]).

**§7.6 Nonlinear waves in shallow and deep water**
- **C11 The hydraulic jump (and the other nonlinear fate: solitons)** $\frac{H_2}{H_1}=\frac12(-1+\sqrt{1+8\mathrm{Fr}_1^2})$
  *(Eq. 7.81)* [#104]: once finite amplitude makes crests outrun troughs (N81 [B]; Fig. 7.19 N185 [B]), the steepened
  front becomes a jump (N82 [B]; Fig. 7.20 N186 [B]) whose CV momentum balance (R10 [B]; N83 [B]; D25) fixes H₂ and
  whose energy loss (N84 [B]; D26) forbids Fr₁ < 1; if dispersion balances steepening instead (N85 [C]), KdV (N93 [B])
  with its linear limit (N94 [B]) and the Ursell ratio (N95 [B]) give permanent forms — cnoidal (N96 [C]) and the
  solitary wave (N97 [B]; Fig. 7.23 N189 [B]).
- **C12 Stokes waves and Stokes drift** $\bar u_L=a^2\omega ke^{2kz_0}$ *(Eq. 7.85)* [#112]: once finite amplitude
  sharpens crests (N86, N87, N88 [B]; Fig. 7.21 N187 [B]) and orbits fail to close (N89 [B]), a first-order Taylor
  expansion of u about the mean position (N90 [B]; D27) gives a mean drift that decays twice as fast as the orbits, for
  any depth (N91 [B]), while the Eulerian mean stays zero (N92 [B]; Fig. 7.22 N188 [B]).

**§7.7 Waves on a density interface**
- **C13 Interfacial waves** $\omega=\varepsilon\sqrt{gk}$, $\varepsilon^2=\frac{\rho_2-\rho_1}{\rho_2+\rho_1}$ *(Eq. 7.95)*
  [#128]: once complex notation (N99 [B]) turns the two-fluid problem (N100–N104 [B], R11 [C]) into algebra (D28), the
  wave is a deep-water wave slowed by ε, the two fluids slide past each other as a vortex sheet (N107 [B]; Fig. 7.24
  N190 [B]) and the energy is $\tfrac12(\rho_2-\rho_1)ga^2$ *(Eq. 7.96)* (N105, N106 [B]); set-up, dead water and the
  column swap named (N98, N108, N191, N192 [C]).
- **C14 Barotropic and baroclinic modes; reduced gravity** $\omega^2=\frac{gk(\rho_2-\rho_1)\sinh kH}{\rho_2\cosh kH+\rho_1\sinh kH}$
  *(Eq. 7.113)* [#150]: once the layer-over-deep-fluid conditions (N110–N118 [B]) fix the constants (N119–N122 [B]; D29)
  and the pressure condition factors into two roots (N123 [B]; D30), the barotropic mode (N124, N125 [B]) and the
  baroclinic mode (N126 [B]) appear, with long-wave speed $c=\sqrt{g'H}$, $g'=g(\rho_2-\rho_1)/\rho_2$
  *(Eqs. 7.116–7.117)* (N127, N128, N129, N130 [B]; D31) and the thin-layer picture (N131 [B]; Figs. 7.27–7.28 N193,
  N194 [B]); set-up named (N109 [C]).

**§7.8 Internal waves in a continuously stratified fluid**
- **C15 Internal waves: ω = N cos θ** *(Eq. 7.139)* [#177]: once the Boussinesq set (R12 [C], R13, R14 [B], R15, R16,
  R17 [C]) is linearised about a resting stratification (N132, N133 [B]; R18 [B] N²; N134, N135 [B]; D32; R19 [C]) and
  reduced to one equation for w (N136, N137, N138 [B]; D33 ★★★), a plane wave (N139, N140 [B]) gives
  $\omega^2=\frac{k^2+l^2}{K^2}N^2$ (N141, N142 [B]; D34): frequency set by direction only, with the limits ω = N and
  ω → 0 (N145, N146 [B]); rotational flow (R20 [C]); blocking named (N147 [C], N196 [C]).
- **C16 Transverse waves with c ⟂ c_g; beams; F = c_g E** $\mathbf c_g\cdot\mathbf c=0$ *(Eq. 7.146)* [#185]: once
  incompressibility makes the motion lie along the crests (N143, N144 [B]; D35) and the gradient of ω(K) (N148, N149
  [B]; D36; Figs. 7.29, 7.31, 7.32 N195, N197, N198 [B]) is perpendicular to K, energy leaves a source along four
  beams at cos θ = ω/N (N150 [B]; Fig. 7.33 N199 [B]), carried with $\mathbf F=\mathbf c_gE$ *(Eq. 7.159)* (energy
  equation and potential energy N152–N155, N157 [B]; polarization and flux N158–N166 [B]; D37); N(z) and the interface
  target named (N151, N156 [C]).
- Exercises: S01 (pointer list; exercise text never reproduced).

## 2. Chapter map (depth) — every inventory row exactly once
One row per inventory row (236: `[#n]` = analysis §2 numbered row, `[#Fn]` = figure row, the X1 tag = exercises). `A parent` names the A block a B or C item is written in. IDs: CORE C01–C16 (teaching order), NOTE N01–N199 and RECAP R01–R20 in inventory order (numbered rows, then figures), SKIP S01.

| ID | Item | § | Depth | Tier | A parent | Reason (A) / treatment (B) / pointer (C, SKIP) / source chapter (RECAP) |
|---|---|---|---|---|---|---|
| N01 | Three wave families and their restoring forces (interface, internal, compression); water waves neither longitudinal nor transverse; chapter assumptions ω ≫ Coriolis, small amplitude ⇒ linear [#1] | 7.1 | C | NOTE | C01 | named in one sentence opening C01; pointers Ch. 13 (rotation: Poincaré, Kelvin, Rossby) and Ch. 15 (compression waves) |
| N02 | Fourier superposition: any linear waveform is a sum of sinusoids, $\eta=\sum_k\hat\eta_k e^{i(kx-\omega(k)t)}+\text{c.c.}$, each moving with its own ω(k) (Exercise 7.3) [#2] | 7.1 | B | NOTE | C01 | stated as the reason (7.1) is the building block; `core.W.linear_evolve` on a two-mode example (one line); a-D67 given (§4c); reused in C09 (packets, stone in a pond) |
| C01 | Eq. (7.1): the sinusoidal travelling wave $\eta(x,t)=a\cos[\frac{2\pi}{\lambda}(x-ct)]$ and its vocabulary — amplitude a, wavelength λ, wavenumber k = 2π/λ, period T = λ/c, frequency ν = 1/T, ω = 2πν, phase; with the plane wave $\eta=a\cos(\mathbf K\cdot\mathbf x-\omega t)$ and $\mathbf c=(\omega/K)\mathbf e_K$ [#3] | 7.1 | A | CORE | – | load-bearing: every result of Ch. 7, 11, 13 and 15 is written in k, ω, c and K; the phase-speed definition c = ω/k (D01) is what dispersion, group velocity and internal waves all measure against |
| N03 | Eq. (7.2): the same wave in k and ω, $\eta(x,t)=a\cos[kx-\omega t]$ [#4] | 7.1 | B | NOTE | C01 | stated with one number (λ = 100 m, T = 8 s → k = 0.0628 rad/m, ω = 0.785 rad/s, c = 12.5 m/s); `core.W.sinusoid`, `wave_parameters` |
| N04 | Eq. (7.3): crest condition $\frac{2\pi}{\lambda}(x_{crest}-ct)=2n\pi=kx_{crest}-\omega t$ [#5] | 7.1 | B | NOTE | C01 | stated as step 2 of D01; `core.W.crest_positions` tracks crest n = 0 in the from-scratch cell |
| N05 | Eq. (7.4): phase speed $c=\omega/k=\lambda\nu$ [#6] | 7.1 | B | NOTE | C01 | stated as D01's result; the tracked crest's Δx/Δt equals ω/k to round-off |
| N06 | Eq. (7.5): 3-D plane wave $\eta=a\cos(kx+ly+mz-\omega t)=a\cos(\mathbf K\cdot\mathbf x-\omega t)$ [#7] | 7.1 | B | NOTE | C01 | stated with `core.W.plane_wave` on a 2-D grid (crest lines ⟂ K); prepares C15 |
| N07 | Eq. (7.6): $K^2=k^2+l^2+m^2$ [#8] | 7.1 | B | NOTE | C01 | stated in one line beside (7.5) |
| N08 | Eq. (7.7): $\lambda=2\pi/K$ — crest spacing measured along K (Fig. 7.1) [#9] | 7.1 | B | NOTE | C01 | stated (a-D02 given, §4c) with one number: k = l = 1 rad/m → K = √2, λ = 4.44 m, while crests are 2π = 6.28 m apart along x |
| N09 | Eq. (7.8): phase-velocity vector $\mathbf c=(\omega/K)\mathbf e_K$, $\mathbf e_K=\mathbf K/K$ (re-displayed with (7.143)) [#10] | 7.1 | B | NOTE | C01 | stated with `core.W.phase_velocity_vector`; forward pointer to C16 where c and c_g are perpendicular |
| N10 | Trace velocities $c_x=\omega/k$, $c_y=\omega/l$, $c_z=\omega/m\ \ge c=\omega/K$ — not components of c [#11] | 7.1 | B | NOTE | C01 | stated with `core.W.trace_velocities` (same numbers as N08: c_x = √2 c); ⚠️ Common confusion callout: reciprocals do not add like components |
| N11 | Eq. (7.9): Doppler shift $\omega_0=\omega+\mathbf U\cdot\mathbf K$ (observed vs intrinsic frequency; frozen pattern ω = 0 ⇒ ω₀ = Uk) [#12] | 7.1 | B | NOTE | C01 | stated (a-D03 given, §4c, via x′ = x − Ut in one line) with `core.W.doppler_frequency`: swell T = 8 s against a 1 m/s current; later ω is intrinsic; Ch. 13 mountain and Rossby waves |
| N12 | Setup of §7.2: constant-density liquid of depth H, $a/\lambda\ll1$, $a/H\ll1$, no surface tension, air ignored, irrotational; x along propagation, z up (Fig. 7.2) [#13] | 7.2 | C | NOTE | C02 | named in the C02 opening sentence; the geometry sketch is N168 |
| R01 | Eq. (7.10): $u=\partial\phi/\partial x$, $w=\partial\phi/\partial z$ [#14] | 7.2 | C | RECAP | C02 | ch03 N29, ch06 C02 (6.10): one sentence |
| R02 | Eq. (7.11): Laplace $\partial^2\phi/\partial x^2+\partial^2\phi/\partial z^2=0$ [#15] | 7.2 | B | RECAP | C02 | ch06 C02 (6.12): reminded with `core.potential.laplacian_residual` (y ↦ z) on the wave φ; a-D04 given (§4c) |
| R03 | Eq. (7.12): no flow through the bottom, $w=\partial\phi/\partial z=0$ on $z=-H$ [#16] | 7.2 | B | RECAP | C02 | ch06 N16 (6.16): wall = streamline; reminded with the 'bottom' residual of `ch07.free_surface_residuals` |
| R04 | Eq. (7.13): kinematic condition $(\mathbf n\cdot\mathbf u)_{z=\eta}=\mathbf n\cdot\mathbf U_s$ [#17] | 7.2 | B | RECAP | C02 | ch04 C14 (4.92)–(4.93), D29: reminded with `core.interfaces.kinematic_bc_residual` called with f = z − η |
| N13 | Eq. (7.14): upward unit normal of $f=z-\eta(x,t)=0$, $\mathbf n=(-\eta_x\mathbf e_x+\mathbf e_z)/\sqrt{\eta_x^2+1}$ [#18] | 7.2 | B | NOTE | C02 | stated as step 2 of D02; ⚠️ notation change vs ch04 (η was the level-set function there, here it is the elevation; §9 R1); `ch07.surface_normal` |
| N14 | Eq. (7.15): the surface point at fixed x moves vertically, $\mathbf U_s=(\partial\eta/\partial t)\mathbf e_z$ [#19] | 7.2 | B | NOTE | C02 | stated in D02 (only the normal part of U_s matters) |
| N15 | Eq. (7.16): exact kinematic condition $(\partial\phi/\partial z)_{z=\eta}=\partial\eta/\partial t+(\partial\eta/\partial x)(\partial\phi/\partial x)_{z=\eta}$ [#20] | 7.2 | B | NOTE | C02 | stated as D02's result; it equals $D(z-\eta)/Dt=0$ (ch04 (4.91)); 'kinematic_exact' residual |
| N16 | Eq. (7.17): drop the product of slope and velocity, $(\partial\phi/\partial z)_{z=\eta}\cong\partial\eta/\partial t$ [#21] | 7.2 | B | NOTE | C02 | stated as step 3 of D03 with the size estimate η_xφ_x ~ ka·aω |
| C02 | Eq. (7.18): the linearised kinematic condition moved to the flat level by Taylor expansion, $(\partial\phi/\partial z)_{z=0}\cong\partial\eta/\partial t$ — and with it the whole linear free-surface problem (7.11), (7.12), (7.18), (7.21) [#22] | 7.2 | A | CORE | – | load-bearing: 'linearise, then transfer the condition to z = 0' is the method of every wave problem in §7.3–7.7, of Kelvin–Helmholtz (Ch. 11) and of shallow-water/Kelvin waves (Ch. 13); D02, D03, D04 |
| R05 | Eq. (7.19): dynamic condition $(p)_{z=\eta}=0$ (gauge; stress-free surface) [#23] | 7.2 | B | RECAP | C02 | ch04 C14 stress continuity: one sentence with gauge vs absolute pressure |
| R06 | Eq. (7.20): linearised unsteady Bernoulli $\partial\phi/\partial t+p/\rho+gz\cong0$ [#24] | 7.2 | B | RECAP | C02 | ch04 C12 (4.75)/(4.83), D26: reminded with `core.bernoulli.unsteady_bernoulli_pressure` (speed 0); the linearisation itself is D04's steps 2–4 |
| N17 | Eq. (7.21): linearised dynamic condition $(\partial\phi/\partial t)_{z=0}\cong-g\eta$ [#25] | 7.2 | B | NOTE | C02 | stated as D04's result; 'dynamic_linear' residual; a table of exact-vs-linear residuals growing like (ka)² (ch01 lesson: make the approximation measurable) |
| N18 | The problem needs an initial surface; the book picks $\eta(x,0)=a\cos kx$ (Exercise 7.3) [#26] | 7.2 | C | NOTE | C03 | named: pointer to N02 (`linear_evolve` takes any initial shape) and S01 |
| N19 | Eq. (7.22): separable trial $\phi=f(z)\sin(kx-\omega(k)t)$ [#27] | 7.2 | B | NOTE | C03 | stated as D05 step 1 (why sine: (7.18) and (7.21) each differentiate once) |
| N20 | Eq. (7.23): $f''-k^2f=0$, $f=Ae^{kz}+Be^{-kz}$ [#28] | 7.2 | B | NOTE | C03 | stated as D05 steps 2–3; ch01 P44 (e^{λz} trial) reminder |
| N21 | Eq. (7.24): bottom condition $B=Ae^{-2kH}$ [#29] | 7.2 | B | NOTE | C03 | stated as D05 step 4 |
| N22 | Eq. (7.25): kinematic condition $k(A-B)=\omega a$ [#30] | 7.2 | B | NOTE | C03 | stated as D05 step 5–6 (A, B solved) |
| N23 | Eq. (7.26): $\phi=\frac{a\omega}{k}\frac{\cosh k(z+H)}{\sinh kH}\sin(kx-\omega t)$ (used as a given in ch04, derived here) [#31] | 7.2 | B | NOTE | C03 | stated as D05's result; `ch07.wave_fields` (phi) and `ch07.surface_wave_sympy` |
| N24 | Eq. (7.27): $u=a\omega\frac{\cosh k(z+H)}{\sinh kH}\cos(kx-\omega t)$, $w=a\omega\frac{\sinh k(z+H)}{\sinh kH}\sin(kx-\omega t)$ [#32] | 7.2 | B | NOTE | C03 | stated (a-D09 given, §4c: one derivative each); parity `ch07.wave_fields` = `ch04.linear_wave_fields` to 1e-12 (the ch04 test field is now derived) |
| N25 | Dynamic condition applied to (7.26) (unnumbered): $-\frac{a\omega^2}{k}\frac{\cosh kH}{\sinh kH}\cos(kx-\omega t)\cong-ag\cos(kx-\omega t)$ [#33] | 7.2 | B | NOTE | C03 | stated as D06 steps 2–4 |
| C03 | Eq. (7.28): the dispersion relation $\omega=\sqrt{gk\tanh kH}$, or $T=\sqrt{\frac{2\pi\lambda}{g}\coth\frac{2\pi H}{\lambda}}$ (ch04 used it as a given test field; derived here for the first time) [#34] | 7.2 | A | CORE | – | load-bearing (the chapter's hub): phase speed, pressure, orbits, energy flux, c_g, seiches, refraction and every later wave family are read off an ω(k); Ch. 13's Poincaré waves reduce to it; D05, D06 |
| C04 | Eq. (7.29): phase speed $c=\sqrt{\frac gk\tanh kH}$ — longer waves are faster (dispersive), with the deep limit $c=\sqrt{g/k}$ (7.45) and the shallow limit $c=\sqrt{gH}$ (7.49) [#35] | 7.2 | A | CORE | – | load-bearing: the deep/shallow classification (kH > 2, H < 0.07λ) and c = √(gH) are what Ch. 13's shallow-water, Kelvin and tidal waves stand on; D07, D08, D09 |
| R07 | Eq. (7.30): $p'\equiv p+\rho gz$ (perturbation from hydrostatic) [#36] | 7.2 | C | RECAP | C04 | ch04 C13 (p = p_s + p′): one sentence |
| N26 | Eq. (7.31): wave pressure $p'=\rho ga\frac{\cosh k(z+H)}{\cosh kH}\cos(kx-\omega t)$ — decays with depth at a rate set by k [#37] | 7.2 | B | NOTE | C04 | stated (a-D12 given, §4c) with `ch07.pressure_response`; number: λ = 50 m on H = 10 m, a = 1 m → ρga = 9.8 kPa at z = 0, 5.2 kPa at the bottom (cosh kH = 1.90) |
| R08 | Eq. (7.32): path lines $dx_p/dt=u(x_p,z_p,t)$, $dz_p/dt=w(x_p,z_p,t)$ [#38] | 7.2 | B | RECAP | C05 | ch03 C04 (3.8): reminded with `core.kinematics.pathline`; ch03 Example 3.1 is the same linearisation |
| N27 | Eq. (7.33): (7.32) with (7.27) inserted — nonlinear in x_p, z_p [#39] | 7.2 | B | NOTE | C05 | stated as D10 step 1; `ch07.particle_path(model='exact')` |
| N28 | Eq. (7.34a, b): linearised path lines — right sides at the mean position (x₀, z₀) [#40] | 7.2 | B | NOTE | C05 | stated as D10 steps 2–4; `particle_path(model='linear')` |
| N29 | Eq. (7.35a, b): excursions $\xi\cong-a\frac{\cosh k(z_0+H)}{\sinh kH}\sin(kx_0-\omega t)$, $\zeta\cong a\frac{\sinh k(z_0+H)}{\sinh kH}\cos(kx_0-\omega t)$ [#41] | 7.2 | B | NOTE | C05 | stated as D10's result; `ch07.orbit_linear`; no secular term |
| C05 | Eq. (7.36): particle orbits are closed ellipses $\xi^2/A^2+\zeta^2/B^2=1$ with $A=a\frac{\cosh k(z_0+H)}{\sinh kH}$, $B=a\frac{\sinh k(z_0+H)}{\sinh kH}$ — circles in deep water, flat ellipses in shallow water, clockwise for a right-going wave [#42] | 7.2 | A | CORE | – | load-bearing: the picture of what the water does under a wave; the base for Stokes drift (C12), the standing-wave and interfacial fields and Ch. 13's float trajectories; D10, D11 |
| N30 | Eq. (7.37): stream function $\psi=\frac{a\omega}{k}\frac{\sinh k(z+H)}{\sinh kH}\cos(kx-\omega t)$ (Exercise 7.4); forward flow under crests, backward under troughs (Fig. 7.5) [#43] | 7.2 | B | NOTE | C05 | stated (a-D15 given, §4c) and checked with `core.streamfunction.velocity_from_streamfunction_2d` (y ↦ z) against (7.27); used by C08 |
| N31 | Eq. (7.38): depth-integrated, wavelength-averaged kinetic energy integral (to z = 0, not η) [#44] | 7.2 | B | NOTE | C06 | stated as D12 steps 1–3; `ch07.wave_energy(method='quad')` (dblquad) |
| N32 | Eq. (7.39): $E_k=\tfrac12\rho g\overline{\eta^2}$ [#45] | 7.2 | B | NOTE | C06 | stated as D12's result |
| N33 | Eq. (7.40): potential energy $E_p=\frac{\rho g}{2\lambda}\int_0^\lambda\eta^2dx$ (column-swap argument, Fig. 7.6) [#46] | 7.2 | B | NOTE | C06 | stated as D13 steps 1–4 |
| N34 | Eq. (7.41): $E_p=\tfrac12\rho g\overline{\eta^2}=E_k$ — equipartition (fails with rotation, Ch. 13) [#47] | 7.2 | B | NOTE | C06 | stated as D13 step 5 with a two-bar figure KE/PE vs kH (always equal) |
| C06 | Eq. (7.42): wave energy per unit horizontal area $E=E_k+E_p=\rho g\overline{\eta^2}=\tfrac12\rho ga^2$, with the flux (7.44) $F=[\tfrac12\rho ga^2][\frac c2(1+\frac{2kH}{\sinh 2kH})]$ [#48] | 7.2 | A | CORE | – | load-bearing: E = ½ρga² and 'flux = energy × a speed' set up group velocity (C09), swell forecasting and Ch. 13's wave-energy budgets; D12, D13, D14 |
| N35 | Eq. (7.43): energy flux $F=\langle\int_{-H}^0p'u\,dz\rangle$ (the background part drops because ⟨u⟩ = 0) [#49] | 7.2 | B | NOTE | C06 | stated as D14 steps 1–4; `ch07.energy_flux(method='quad')` |
| N36 | Eq. (7.44): $F=[\tfrac12\rho ga^2][\frac c2(1+\frac{2kH}{\sinh2kH})]$ — the second factor is a speed [#50] | 7.2 | B | NOTE | C06 | stated as D14's result; named 'group speed' with a pointer to C09 (7.69); typo T2 ('u from (7.28)' means (7.27)) noted |
| N37 | Hyperbolic functions (Fig. 7.7): cosh x ≈ 1, sinh x ≈ tanh x ≈ x for small x; tanh → 1 for large x; tanh 2 = 0.964 [#51] | 7.2 | B | NOTE | C04 | stated with the Fig. 7.7 remake (N173); the primer comes before D05 (§4) |
| N38 | Eq. (7.45): deep water $c=\sqrt{g/k}=\sqrt{g\lambda/2\pi}$, within 2 % for kH > 2 (H > 0.32λ) [#52] | 7.2 | B | NOTE | C04 | stated as D08's result; `core.W.depth_regime` returns the error, not only the label |
| N39 | Ocean scales: wind waves T ~ 10 s ↔ λ ≈ 156 m (book rounds to 150 m) are deep-water waves over the shelf (100 m) and open ocean (4 km); tides and tsunamis are shallow [#53] | 7.2 | B | NOTE | C04 | stated as C04's worked number with `core.W.wavelength_from_period(10, H)` for H = 100, 4000, ∞ m, plus a tsunami (T = 20 min, H = 4 km → c ≈ 198 m/s) |
| N40 | Eq. (7.46): deep-water orbits are circles $\xi=-ae^{kz_0}\sin(kx_0-\omega t)$, $\zeta=ae^{kz_0}\cos(kx_0-\omega t)$ [#54] | 7.2 | B | NOTE | C05 | stated as D08 step 5 (reused in C05); radius halves every 0.11λ of depth |
| N41 | Eq. (7.47): deep-water velocities $u=a\omega e^{kz}\cos(kx-\omega t)$, $w=a\omega e^{kz}\sin(kx-\omega t)$ — the vector turns clockwise at ω with constant size [#55] | 7.2 | B | NOTE | C05 | stated with `ch07.wave_fields(H=np.inf)`; used by D27 (Stokes drift) |
| N42 | Eq. (7.48): deep-water pressure $p'=\rho gae^{kz}\cos(kx-\omega t)$ — 4 % of its surface value at depth λ/2 [#56] | 7.2 | B | NOTE | C04 | stated as D08 step 6 with e^{−π} = 0.043; bottom pressure gauges filter out short waves |
| N43 | Eq. (7.49): shallow water $c=\sqrt{gH}$ — nondispersive, within 3 % for H < 0.07λ (matches ch04's bore limit) [#57] | 7.2 | B | NOTE | C04 | stated as D09's result; parity with the √(gh) limit of `ch04.bore_speed`; number: H = 4 km → 198 m/s |
| N44 | Eq. (7.50): shallow-water orbits $\xi=-\frac{a}{kH}\sin(kx_0-\omega t)$, $\zeta=a(1+\frac zH)\cos(kx_0-\omega t)$ — flat ellipses of depth-independent width [#58] | 7.2 | B | NOTE | C05 | stated as D09 step 5 (reused in C05) |
| N45 | Eq. (7.51): shallow-water velocities $u=\frac{a\omega}{kH}\cos(kx-\omega t)$, $w=a\omega(1+\frac zH)\sin(kx-\omega t)$; w ≪ u [#59] | 7.2 | B | NOTE | C05 | stated with w/u ~ kH (one number: kH = 0.2) |
| N46 | Eq. (7.52): shallow-water pressure is hydrostatic, $p'=\rho ga\cos(kx-\omega t)=\rho g\eta$ [#60] | 7.2 | B | NOTE | C04 | stated as D09 step 6; `pressure_response` → 1 as kH → 0; climate hook: the hydrostatic shallow-water model of Ch. 13 |
| N47 | Refraction on a sloping beach (Fig. 7.8): ω fixed along the path, $c=\sqrt{gH}$ and λ shrink as H falls, so crests turn parallel to the depth contours [#61] | 7.2 | B | NOTE | C10 | stated inside C10 with D24 (Snell form k sin α = const, written out because the book gives words only); `core.W.ray_trace` on a plane beach; §7.2 cell carries a forward pointer |
| N48 | Refraction around a circular island with a gently sloping beach (Fig. 7.9): crests reach the shadow side [#62] | 7.2 | B | NOTE | C10 | stated with a ray figure from `core.W.ray_trace` over a circular H(r) (our code) |
| N49 | Capillary waves: surface tension as a second restoring force (a stretched membrane); §7.2 is extended, not replaced [#63] | 7.3 | C | NOTE | C07 | named in the C07 opening sentence |
| R09 | Laplace pressure jump (1.5) $\Delta p=\sigma(1/R_1+1/R_2)$, higher pressure on the side of the centres of curvature [#64] | 7.3 | B | RECAP | C07 | ch01 (1.5), re-derived ch04 C14: reminded with `core.interfaces.laplace_jump_from_balance` |
| N50 | Eq. (7.53): curvature of the graph z = η(x), $p_a-(p)_{z=\eta}=\sigma\frac{\eta_{xx}}{(1+\eta_x^2)^{3/2}}\cong\sigma\eta_{xx}$ [#65] | 7.3 | B | NOTE | C07 | stated as D15 steps 1–4; `ch07.curvature` |
| N51 | Eq. (7.54): $(p)_{z=\eta}=-\sigma\,\partial^2\eta/\partial x^2$ (gauge) [#66] | 7.3 | B | NOTE | C07 | stated as D15 step 5; number: a = 1 mm, λ = 1 cm → p under a crest ≈ +29 Pa |
| N52 | Eq. (7.55): dynamic condition with surface tension $(\partial\phi/\partial t)_{z=0}=\frac\sigma\rho\eta_{xx}-g\eta$ [#67] | 7.3 | B | NOTE | C07 | stated as D15 step 6 |
| N53 | Eq. (7.56): capillary–gravity dispersion $\omega=\sqrt{k(g+\sigma k^2/\rho)\tanh kH}$ [#68] | 7.3 | B | NOTE | C07 | stated as D15's result (every g of D06 becomes g + σk²/ρ); `core.W.omega_capillary_gravity` |
| C07 | Eq. (7.57): phase speed with surface tension $c=\sqrt{(\frac gk+\frac{\sigma k}{\rho})\tanh kH}$ — a capillary branch, a gravity branch and a minimum $c_{min}=(4g\sigma/\rho)^{1/4}$ at $\lambda_m=2\pi\sqrt{\sigma/\rho g}$ (7.58) [#69] | 7.3 | A | CORE | – | load-bearing: two restoring forces make c(λ) non-monotonic — the minimum explains ripples, the calm ring of a stone in a pond (C09) and the σ-term of Kelvin–Helmholtz (Ch. 11); D15, D16 |
| N54 | Eq. (7.58): $c_{min}=[4g\sigma/\rho]^{1/4}$ at $\lambda_m=2\pi\sqrt{\sigma/(\rho g)}$ [#70] | 7.3 | B | NOTE | C07 | stated as D16's result; `core.W.capillary_minimum`; from-scratch scan of c(λ) |
| N55 | Eq. (7.59): air–water at 20 °C: c_min ≈ 23.1 cm/s at λ_m ≈ 1.71 cm [#71] | 7.3 | B | NOTE | C07 | stated with our numbers (σ from `ch01.surface_tension_water(293.15)` = 72.74 mN/m, ρ = 998.2 → 23.1 cm/s, 1.71 cm); book values only in the private JSON |
| N56 | Ripples: only λ ≲ 7 cm feel σ, λ < 4 mm are pure capillary; surfactants change σ [#72] | 7.3 | C | NOTE | C07 | named in one sentence after N55 with the regime bands on the c(λ) figure |
| N57 | Eq. (7.60): pure capillary waves $c=\sqrt{2\pi\sigma/(\rho\lambda)}$ [#73] | 7.3 | B | NOTE | C07 | stated (a-D26 given, §4c: drop gλ/2π, tanh → 1); `core.W.phase_speed(g=0)` |
| N58 | Eq. (7.61): the left-going wave $\eta=a\cos[kx+\omega t]$ [#74] | 7.4 | B | NOTE | C08 | stated as D17 step 1; `core.W.sinusoid(direction=-1)` |
| C08 | Standing waves: $\eta=a\cos(kx-\omega t)+a\cos(kx+\omega t)=2a\cos kx\cos\omega t$ — fixed nodes; walls pick discrete modes, the lake seiche $\omega=\sqrt{\frac{\pi g(n+1)}{L}\tanh\frac{(n+1)\pi H}{L}}$ (7.65) [#75] | 7.4 | A | CORE | – | load-bearing: superposition of opposite waves and 'boundaries select eigenmodes' — seiches, harbour and bay resonance, tidal basins (Ch. 13 Merian); D17, D18 |
| N59 | Eq. (7.62): standing-wave stream function $\psi=\frac{2a\omega}{k}\frac{\sinh k(z+H)}{\sinh kH}\sin kx\sin\omega t$ (Fig. 7.11) [#76] | 7.4 | B | NOTE | C08 | stated as D17's result; `ch07.standing_wave_fields` asserted equal to the sum of two `wave_fields` |
| N60 | Seiche: a standing oscillation in a closed basin (length L, depth H, vertical walls): only discrete wavelengths and frequencies [#77] | 7.4 | B | NOTE | C08 | stated as the setup of D18 (eigenvalue idea in words, ch02 P80 reminder) |
| N61 | Eq. (7.63): $u=2a\omega\frac{\cosh k(z+H)}{\sinh kH}\sin kx\sin\omega t$ [#78] | 7.4 | B | NOTE | C08 | stated as D17 step 7 (u = ∂ψ/∂z); nodes of η are antinodes of u |
| N62 | Eq. (7.64): allowed wavelengths $kL=(n+1)\pi\Rightarrow\lambda=2L/(n+1)$ [#79] | 7.4 | B | NOTE | C08 | stated as D18 steps 2–3; `ch07.seiche_modes` |
| N63 | Eq. (7.65): natural frequencies of the lake $\omega=\sqrt{\frac{\pi g(n+1)}{L}\tanh[\frac{(n+1)\pi H}{L}]}$; rectangular basins $k^2=(m\pi/L)^2+(n\pi/b)^2$ (Exercises 7.5–7.6) [#80] | 7.4 | B | NOTE | C08 | stated as D18's result; worked number: a lake L = 50 km, H = 100 m → T₀ ≈ 53 min (our inputs); `ch07.basin_modes`; T8 book-value discrepancy stays in the private JSON |
| N64 | Dispersive interface waves in general (Rayleigh, Stoneley, liquid–liquid; Graff) [#81] | 7.5 | C | NOTE | C09 | named in one sentence opening C09 |
| N65 | Eq. (7.66): two nearby waves make beats, $\eta=2a\cos(\tfrac12\Delta k\,x-\tfrac12\Delta\omega\,t)\cos(kx-\omega t)$ (book prints ½Δω x; read ½Δω t, T1) [#82] | 7.5 | B | NOTE | C09 | stated as D19 steps 1–4; `core.W.beat_wave`; the printed form shown beside the correct one (envelope would not move) |
| C09 | Eq. (7.67): group velocity $c_g=\Delta\omega/\Delta k\to d\omega/dk$ — the envelope (and the energy) moves at c_g, the crests at c; tangent vs chord of ω(k) (Fig. 7.14); for water waves $c_g=\frac c2[1+\frac{2kH}{\sinh 2kH}]$ (7.69) and $F=Ec_g$ (7.71) [#83] | 7.5 | A | CORE | – | load-bearing (second hub): energy, packets, swell arrival, rays (C10), internal-wave beams (C16) and Ch. 13's Rossby/Poincaré energy propagation all travel at c_g; D19, D20, D21 |
| N66 | Wave packet: narrow band δk around k, envelope length ∝ 1/δk (Fig. 7.15) [#84] | 7.5 | B | NOTE | C09 | stated with a Gaussian packet from `core.W.linear_evolve` and its spectrum (N181) |
| N67 | Eq. (7.68): short-time packet evolution $\eta=a(x-c_gt)\cos(kx-\omega t)$ — the amplitude travels at c_g; nodes move at c_g and no energy crosses them [#85] | 7.5 | B | NOTE | C09 | stated as D20's result (★★★, the book only cites Phillips — written out); measured envelope-peak speed (`core.W.envelope`, Hilbert) = c_g within 1 % |
| N68 | Eq. (7.69): $c_g=\frac c2[1+\frac{2kH}{\sinh(2kH)}]$ [#86] | 7.5 | B | NOTE | C09 | stated as D21 steps 1–6; `core.W.group_velocity` |
| N69 | Eq. (7.70): $c_g=c/2$ (deep), $c_g=c$ (shallow); pure capillary $c_g=3c/2$ (Exercise 7.9) [#87] | 7.5 | B | NOTE | C09 | stated as D21 steps 7–8; capillary 3c/2 given (a-D33, §4c) with `group_velocity(g=0)` |
| N70 | Eq. (7.71): energy flux $F=E\frac c2[1+\frac{2kH}{\sinh 2kH}]=Ec_g$, $E=\rho ga^2/2$ [#88] | 7.5 | B | NOTE | C09 | stated as D21 steps 9–10; asserted `energy_flux` = `wave_energy_density` × `group_velocity`; swell number: a = 1 m, T = 10 s deep → F ≈ 38 kW per metre of crest |
| N71 | 3-D group velocity $c_{gi}=\partial\omega/\partial K_i$ (gradient in wavenumber space) [#89] | 7.5 | B | NOTE | C09 | stated with `core.W.group_velocity_vector` on the deep-water ω(K) = √(gK); forward pointer to C16 (7.143) |
| N72 | Stone in a pond (Fig. 7.16): the train sorts itself, longest waves in front; a calm centre because c_g has a minimum ≈ 17.8 cm/s with surface tension (Exercise 7.10); heights fall as the train lengthens; viscosity ends it (Exercise 7.11) [#90] | 7.5 | B | NOTE | C09 | stated with `core.W.min_group_velocity` (a-D33 given) and a capillary–gravity FFT impulse (`linear_evolve`, FAST grid); viscous decay a₀e^{−2νk²t} named with pointer to Ch. 8 (a-D68 given, §4c) |
| N73 | Following one crest of a deep-water group: it runs from the rear to the front of the group and dies; new crests are born at the rear [#91] | 7.5 | B | NOTE | C09 | stated as the crest marker in the packet animation (A2) and the explainer E5 |
| N74 | Eq. (7.72): slowly varying train $\eta=a(x,t)\cos[\theta(x,t)]$ [#92] | 7.5 | B | NOTE | C10 | stated as D22 step 1 |
| N75 | Eq. (7.73): $k\equiv\partial\theta/\partial x$, $\omega\equiv-\partial\theta/\partial t$ [#93] | 7.5 | B | NOTE | C10 | stated as D22 step 2; `ch07.local_wavenumber_frequency` on a chirped train |
| N76 | Eq. (7.74): crest conservation $\partial k/\partial t+\partial\omega/\partial x=0$ [#94] | 7.5 | B | NOTE | C10 | stated as D22 steps 3–5; `crest_conservation_residual` ~ 1e-9 |
| N77 | Eq. (7.75): $\partial k/\partial t+c_g\partial k/\partial x=0$ — k is carried at c_g; x–t diagram with rays dx/dt = c_g and crest lines dx/dt = c (Fig. 7.17) [#95] | 7.5 | B | NOTE | C10 | stated as D22's result (characteristic reading); x–t figure N183 |
| N78 | Eq. (7.76): slowly varying depth, $\omega=\sqrt{gk\tanh[kH(x)]}$ is of the form ω(k, x) [#96] | 7.5 | B | NOTE | C10 | stated as D23 step 1 |
| N79 | Eq. (7.77): $\partial\omega(k,x)/\partial k=c_g$ at fixed x [#97] | 7.5 | B | NOTE | C10 | stated as D23 step 2 |
| N80 | Eq. (7.78): $c_g\,\partial k/\partial t=\partial\omega/\partial t$ [#98] | 7.5 | B | NOTE | C10 | stated as D23 steps 3–4 |
| C10 | Eq. (7.79): frequency is constant along rays, $\partial\omega/\partial t+c_g\,\partial\omega/\partial x=0$ (3-D: $\partial\omega/\partial t+\mathbf c_g\cdot\nabla\omega=0$) while k, c, c_g change — rays bend (Fig. 7.18); with it, refraction on beaches and around islands [#99] | 7.5 | A | CORE | – | load-bearing: kinematic wave (ray/WKB) theory is how Ch. 13 follows internal waves in N(z), topographic and Rossby waves; it explains refraction (§7.2); D22, D23, D24 |
| N81 | Nonlinear shallow-water waves: a point of the profile moves at $c=c_0'+u$, $c_0'=\sqrt{gH_0}$ — crests overtake troughs, the front steepens and breaks (Fig. 7.19) [#100] | 7.6 | B | NOTE | C11 | stated (a-D37 given, §4c) with `ch07.simple_wave_evolve` (our Riemann simple wave, labelled) and its breaking time; pointer Ch. 15 characteristics |
| N82 | Hydraulic jump: after breaking, a front between √(gH₁) and √(gH₂) (spillway, stationary, moving; tidal bores, circular jump) [#101] | 7.6 | B | NOTE | C11 | stated as C11's picture; moving-jump speed parity with `ch04.bore_speed` in the jump frame |
| R10 | Froude number (4.104) $\mathrm{Fr}\equiv u/\sqrt{gH}=u/c$; super- vs subcritical [#102] | 7.6 | B | RECAP | C11 | ch04 C15: reminded with `core.similarity.froude_number` |
| N83 | Eq. (7.80): CV momentum across a stationary jump $Q^2(\frac1{H_2}-\frac1{H_1})=\frac12g(H_1^2-H_2^2)$ [#103] | 7.6 | B | NOTE | C11 | stated as D25 steps 1–5 (ch04 C04 CV momentum recalled); `ch07.jump_momentum_residual` |
| C11 | Eq. (7.81): the conjugate-depth (Bélanger) relation $\frac{H_2}{H_1}=\frac12(-1+\sqrt{1+8\mathrm{Fr}_1^2})$, with the energy loss $E_2-E_1=-\frac{g(H_2-H_1)^3}{4H_1H_2}$ that forbids Fr₁ < 1; and when nonlinearity is balanced by dispersion instead (KdV, solitary wave) [#104] | 7.6 | A | CORE | – | load-bearing: the chapter's finite-amplitude control-volume result — the shock analogue of Ch. 15, bores and internal jumps (Ch. 13), and the second law choosing the physical root; D25, D26 |
| N84 | Jump energy loss per unit mass $E_2-E_1=-(H_2-H_1)\frac{g(H_2-H_1)^2}{4H_1H_2}<0$; Fr₁ < 1 would create energy [#105] | 7.6 | B | NOTE | C11 | stated as D26's result; number: H₁ = 0.1 m, Fr₁ = 3 → H₂ = 0.377 m, head loss 0.141 m |
| N85 | In a dispersive medium steepening can be balanced by dispersion ⇒ waves of permanent form [#106] | 7.6 | C | NOTE | C11 | named as the bridge sentence to N93 (KdV) |
| N86 | Eq. (7.82): Stokes wave $\eta=a\cos k(x-ct)+\tfrac12ka^2\cos2k(x-ct)+\tfrac38k^2a^3\cos3k(x-ct)+\ldots$ [#107] | 7.6 | B | NOTE | C12 | stated (a-D40 given, §4c; the O(ka) coefficient ½ checked by a sympy cell); `ch07.stokes_wave_profile` at ka = 0.3 (Fig. 7.21 remake N187) |
| N87 | Eq. (7.83): amplitude-dependent speed $c=\sqrt{\frac gk(1+k^2a^2+\ldots)}$ [#108] | 7.6 | B | NOTE | C12 | stated with one number (ka = 0.3 → 4.4 % faster); `ch07.stokes_wave_speed` |
| N88 | Limiting Stokes wave: $a_{max}\approx0.07\lambda$, 120° crest; white caps beyond [#109] | 7.6 | B | NOTE | C12 | stated with the cited constant `STOKES_LIMIT_STEEPNESS` = 0.1410633 (H/λ; book's a ≈ half the height) |
| N89 | Stokes drift: at finite amplitude orbits do not close — a mean Lagrangian drift although the Eulerian mean is zero; a dyed vertical line bends forward (Fig. 7.22) [#110] | 7.6 | B | NOTE | C12 | stated as C12's picture; `particle_path(model='exact')` drift per period |
| N90 | Eq. (7.84a, b): path lines with u Taylor-expanded about the mean position to first order in ξ, ζ [#111] | 7.6 | B | NOTE | C12 | stated as D27 steps 1–3; `particle_path(model='taylor1')` |
| C12 | Eq. (7.85): deep-water Stokes drift $\bar u_L=a^2\omega ke^{2kz_0}$ (general depth (7.86)) — a mean forward current made only by following the particles [#112] | 7.6 | A | CORE | – | load-bearing: the Lagrangian–Eulerian distinction in a wave field; Stokes drift drives ocean-surface transport, Langmuir circulation and the Stokes–Coriolis force of Ch. 13 ocean models; D27 |
| N91 | Eq. (7.86): Stokes drift at any depth $\bar u_L=a^2\omega k\frac{\cosh 2k(z_0+H)}{2\sinh^2 kH}$ (Exercise 7.14); zero vertical drift (Exercise 7.13) [#113] | 7.6 | B | NOTE | C12 | stated (a-D42 given, §4c: D27's moves with cosh/sinh) and checked against `stokes_drift_numeric` from exact path lines |
| N92 | The Eulerian mean velocity at a point always under water is zero in periodic irrotational flow [#114] | 7.6 | B | NOTE | C12 | stated (a-D43 given, §4c) with `ch07.eulerian_mean_u` = 0 to 1e-12 beside the nonzero Lagrangian mean |
| N93 | Eq. (7.87): Korteweg–de Vries $\eta_t+c_0\eta_x+\frac32c_0\frac\eta H\eta_x+\frac16c_0H^2\eta_{xxx}=0$, $c_0=\sqrt{gH}$ [#115] | 7.6 | B | NOTE | C11 | stated term by term (nondispersive · nonlinear · weakly dispersive) with `ch07.kdv_solve` (FAST grid, invariants printed): a hump splitting into solitons (animation A5) |
| N94 | Linearised KdV: $c=c_0(1-\tfrac16k^2H^2)$ = the first two Taylor terms of (7.29) [#116] | 7.6 | B | NOTE | C11 | stated (a-D44 given, §4c) with `kdv_linear_phase_speed` vs `phase_speed` on log axes (error ∝ (kH)⁴) |
| N95 | Ursell ratio $\sim a\lambda^2/H^3$: ≳ 16 steepens to a jump, lower balances (book says '(7.88)', read (7.87)) [#117] | 7.6 | B | NOTE | C11 | stated (a-D45 given, §4c) with `ch07.ursell_number` as C11's regime switch |
| N96 | Cnoidal waves: periodic permanent-form KdV solutions (Jacobi cn) [#118] | 7.6 | C | NOTE | C11 | named with an optional `ch07.cnoidal_wave` curve on the Fig. 7.23 remake |
| N97 | Eq. (7.88): solitary wave $\eta=a\,\mathrm{sech}^2[(\frac{3a}{4H^3})^{1/2}(x-ct)]$, $c=c_0(1+\frac a{2H})$ — taller is faster (Exercise 7.15) [#119] | 7.6 | B | NOTE | C11 | stated (a-D46 given, §4c) with `ch07.kdv_residual_sympy()` = 0 and `solitary_wave`; number: a = 0.2 m on H = 1 m → c = 1.1√(gH) = 3.44 m/s |
| N98 | Interface setup: lighter ρ₁ over heavier ρ₂, both deep, small slopes (Fig. 7.24); estuaries, fjords, sun-warmed layers [#120] | 7.7 | C | NOTE | C13 | named in the C13 opening sentence |
| N99 | Complex notation: $\zeta=\mathrm{Re}\{a\,e^{i(kx-\omega t)}\}$, Re dropped during the algebra; complex amplitudes carry phase [#121] | 7.7 | B | NOTE | C13 | stated with `core.W.real_field` and the primer (§4); ⚠️ products need real parts first (used again in C16) |
| N100 | Eq. (7.89): $\zeta=a\exp[i(kx-\omega t)]$ [#122] | 7.7 | B | NOTE | C13 | stated as D28 step 1 |
| R11 | Eq. (7.90): Laplace in each fluid [#123] | 7.7 | C | RECAP | C13 | ch06 C02 (6.12): one sentence |
| N101 | Eq. (7.91): $\phi_1\to0$ as $z\to\infty$ [#124] | 7.7 | B | NOTE | C13 | stated as D28 step 2 |
| N102 | Eq. (7.92): $\phi_2\to0$ as $z\to-\infty$ [#125] | 7.7 | B | NOTE | C13 | stated as D28 step 2 |
| N103 | Eq. (7.93): $\partial\phi_1/\partial z=\partial\phi_2/\partial z=\partial\zeta/\partial t$ at z = 0 [#126] | 7.7 | B | NOTE | C13 | stated as D28 steps 3–5 (two-sided (7.18)) |
| N104 | Eq. (7.94): pressure continuity $\rho_1\partial_t\phi_1+\rho_1g\zeta=\rho_2\partial_t\phi_2+\rho_2g\zeta$ at z = 0 [#127] | 7.7 | B | NOTE | C13 | stated as D28 steps 6–8 |
| C13 | Eq. (7.95): interfacial waves $\omega=\sqrt{gk\frac{\rho_2-\rho_1}{\rho_2+\rho_1}}=\varepsilon\sqrt{gk}$ — deep-water waves slowed by ε, with opposite horizontal velocities (a vortex sheet) and energy $E=\tfrac12(\rho_2-\rho_1)ga^2$ (7.96) [#128] | 7.7 | A | CORE | – | load-bearing: the first internal wave — small Δρ makes waves slow and large; the base state of Kelvin–Helmholtz (Ch. 11) and of two-layer ocean/atmosphere models (Ch. 13); D28 |
| N105 | Interfacial energies $E_k=E_p=\tfrac14(\rho_2-\rho_1)ga^2$ (Exercise 7.18; column swap, Fig. 7.25) [#129] | 7.7 | B | NOTE | C13 | stated (a-D48 given, §4c) with `ch07.interface_energy` (closed vs quad) |
| N106 | Eq. (7.96): $E=\tfrac12(\rho_2-\rho_1)ga^2$ — for equal energy an internal wave is √(ρ₂/Δρ) times taller [#130] | 7.7 | B | NOTE | C13 | stated with a number: thermocline Δρ = 2 kg/m³ → same energy as a 1 m surface wave needs a ≈ 22 m |
| N107 | Opposite horizontal velocities $u_1=-\omega ae^{-kz}e^{i(kx-\omega t)}$, $u_2=\omega ae^{kz}e^{i(kx-\omega t)}$: the interface is a vortex sheet; continuous stratification is rotational everywhere [#131] | 7.7 | B | NOTE | C13 | stated (a-D49 given, §4c) with ch05 D23 sheet strength (γ = 2ωa cos(kx − ωt)); pointer Ch. 11 KH and C15 |
| N108 | 'Dead water': a ship making interfacial waves feels extra drag (Fig. 7.26) [#132] | 7.7 | C | NOTE | C13 | named with pointer to Ch. 13 two-layer flows |
| N109 | Setup: upper layer H with a free surface over a deep lower layer; two modes, surface and interface in or out of phase (Fig. 7.27) [#133] | 7.7 | C | NOTE | C14 | named in the C14 opening sentence (the explainer E8 stage) |
| N110 | Eq. (7.97): $\phi_2\to0$ at $z\to-\infty$ [#134] | 7.7 | B | NOTE | C14 | stated as D29 step 1 |
| N111 | Eq. (7.98): $\partial\phi_1/\partial z=\partial\eta/\partial t$ at z = 0 (book prints dz, T4) [#135] | 7.7 | B | NOTE | C14 | stated as D29 step 2 |
| N112 | Eq. (7.99): $\partial\phi_1/\partial t+g\eta=0$ at z = 0 [#136] | 7.7 | B | NOTE | C14 | stated as D29 step 2 |
| N113 | Eq. (7.100): $\partial\phi_1/\partial z=\partial\phi_2/\partial z=\partial\zeta/\partial t$ at z = −H [#137] | 7.7 | B | NOTE | C14 | stated as D29 step 3 |
| N114 | Eq. (7.101): pressure continuity at z = −H [#138] | 7.7 | B | NOTE | C14 | stated as D30 step 1 |
| N115 | Eq. (7.102): $\eta=ae^{i(kx-\omega t)}$ [#139] | 7.7 | B | NOTE | C14 | stated as D29 step 4 |
| N116 | Eq. (7.103): $\zeta=be^{i(kx-\omega t)}$, b complex [#140] | 7.7 | B | NOTE | C14 | stated as D29 step 4 |
| N117 | Eq. (7.104): $\phi_1=(Ae^{kz}+Be^{-kz})e^{i(kx-\omega t)}$ [#141] | 7.7 | B | NOTE | C14 | stated as D29 step 5 |
| N118 | Eq. (7.105): $\phi_2=Ce^{kz}e^{i(kx-\omega t)}$ (book prints e^{i(kz−ωt)}; read kx, T3) [#142] | 7.7 | B | NOTE | C14 | stated as D29 step 5 with the typo shown and a failing sympy residual for the printed form |
| N119 | Eq. (7.106): $A=-\frac{ia}2(\frac\omega k+\frac g\omega)$ [#143] | 7.7 | B | NOTE | C14 | stated as D29 steps 6–8 |
| N120 | Eq. (7.107): $B=\frac{ia}2(\frac\omega k-\frac g\omega)$ [#144] | 7.7 | B | NOTE | C14 | stated as D29 steps 6–8 |
| N121 | Eq. (7.108): $C=-\frac{ia}2(\frac\omega k+\frac g\omega)-\frac{ia}2(\frac\omega k-\frac g\omega)e^{2kH}$ [#145] | 7.7 | B | NOTE | C14 | stated as D29 steps 9–10 |
| N122 | Eq. (7.109): $b=\frac a2(1+\frac{gk}{\omega^2})e^{-kH}+\frac a2(1-\frac{gk}{\omega^2})e^{kH}$ [#146] | 7.7 | B | NOTE | C14 | stated as D29's result; `ch07.two_layer_sympy` |
| N123 | Eq. (7.110): two-layer dispersion $(\frac{\omega^2}{gk}-1)\{\frac{\omega^2}{gk}[\rho_1\sinh kH+\rho_2\cosh kH]-(\rho_2-\rho_1)\sinh kH\}=0$ (Exercise 7.19) [#147] | 7.7 | B | NOTE | C14 | stated as D30's result (★★★, written out); `core.W.two_layer_free_surface_omega` |
| N124 | Eq. (7.111): first root $\omega^2=gk$ (the surface wave) [#148] | 7.7 | B | NOTE | C14 | stated as D31 step 1 |
| N125 | Eq. (7.112): barotropic mode $b=ae^{-kH}$ — interface in phase, smaller [#149] | 7.7 | B | NOTE | C14 | stated as D31 step 2 |
| C14 | Eq. (7.113): the baroclinic (internal) mode $\omega^2=\frac{gk(\rho_2-\rho_1)\sinh kH}{\rho_2\cosh kH+\rho_1\sinh kH}$ with surface and interface in antiphase, $\eta=-\zeta\frac{\rho_2-\rho_1}{\rho_1}e^{-kH}$ (7.114), and for long waves $c=\sqrt{g'H}$, $g'=g\frac{\rho_2-\rho_1}{\rho_2}$ (7.116)–(7.117) [#150] | 7.7 | A | CORE | – | load-bearing: barotropic vs baroclinic modes and reduced gravity are the vocabulary of layered ocean and atmosphere models (Ch. 13: reduced-gravity models, baroclinic Rossby radius, equatorial waves); D29, D30, D31 |
| N126 | Eq. (7.114): baroclinic mode $\eta=-\zeta(\frac{\rho_2-\rho_1}{\rho_1})e^{-kH}$ — surface barely moves, u reverses across the interface [#151] | 7.7 | B | NOTE | C14 | stated as D31 step 4; number: Δρ/ρ = 0.002, 10 m thermocline wave → 2 cm surface signal (long waves) |
| N127 | Eq. (7.115): long-wave baroclinic dispersion $\omega^2=kg\frac{\rho_2-\rho_1}{\rho_2}kH$ [#152] | 7.7 | B | NOTE | C14 | stated as D31 steps 6–7 |
| N128 | Eqs. (7.116), (7.117): internal long-wave speed $c=[g'H]^{1/2}$, $g'=g(\rho_2-\rho_1)/\rho_2$ (ρ₂ in the denominator; ch04's reduced_gravity uses ρ₁) [#153] | 7.7 | B | NOTE | C14 | stated as D31's result; ⚠️ Which ρ? callout: both conventions with the size of the difference (0.2 % ocean, 25 % for oil (800 kg/m³) over water), `core.W.reduced_gravity_book(ref='lower')`; number: H = 50 m, Δρ/ρ = 0.002 → c ≈ 0.99 m/s vs √(gH) = 22 m/s |
| N129 | Eq. (7.118): long-wave mode ratio $\eta=-\zeta(\rho_2-\rho_1)/\rho_1$ [#154] | 7.7 | B | NOTE | C14 | stated as D31 step 9 |
| N130 | Eq. (7.119): upper-layer pressure is hydrostatic in the long-wave limit, $p'=\rho_1g\eta$ [#155] | 7.7 | B | NOTE | C14 | stated as D31 step 10 (shallow ⇔ hydrostatic, again) |
| N131 | Both layers thin (Boussinesq limit): barotropic = depth-uniform u, baroclinic = opposite uniform u in each layer (Fig. 7.28); two layers between rigid lids (Exercise 7.20) [#156] | 7.7 | B | NOTE | C14 | stated with mode-profile figure from `ch07.two_layer_modes` at small kH and `two_layer_rigid_lid_omega` |
| R12 | Assumptions of §7.8: Boussinesq, inviscid, small amplitude, ω ≫ Coriolis [#157] | 7.8 | C | RECAP | C15 | ch04 C13: one sentence opening C15 |
| R13 | (4.9) $D\rho/Dt=0$ and (4.10) $\nabla\cdot\mathbf u=0$ as printed at the start of §7.8 [#158] | 7.8 | B | RECAP | C15 | ch04 C02: reminded with `core.navier_stokes.continuity_residual` |
| R14 | Eqs. (7.120), (7.121), (7.122): linear inviscid Boussinesq momentum $\partial_t\mathbf u=-\nabla p/\rho_0-(\rho g/\rho_0)\mathbf e_z$ [#159] | 7.8 | B | RECAP | C15 | ch04 C13 (4.86), D27: reminded; the linearisation is D32 step 1 |
| R15 | Why Dρ/Dt = 0: nondiffusive heat/salt + an equation of state without pressure, $\delta\rho/\rho=-\alpha\delta T$, $\beta\delta S$ [#160] | 7.8 | C | RECAP | C15 | ch01 (thermal expansion), ch04 C02: one sentence |
| R16 | Eq. (7.123): quiescent hydrostatic base $0=-\frac1{\rho_0}\frac{d\bar p}{dz}-\frac{\bar\rho g}{\rho_0}$ [#161] | 7.8 | C | RECAP | C15 | ch01 C20 (1.8): one line in D32 |
| R17 | Eq. (7.124): $p=\bar p(z)+p'$, $\rho=\bar\rho(z)+\rho'$ [#162] | 7.8 | C | RECAP | C15 | ch04 C13: one line in D32 |
| N132 | Eq. (7.125): density equation with the split [#163] | 7.8 | B | NOTE | C15 | stated as D32 step 2 |
| N133 | Eq. (7.126): linearised density equation $\partial\rho'/\partial t+w\,d\bar\rho/dz=0$ — ρ′ made only by vertical advection of the background [#164] | 7.8 | B | NOTE | C15 | stated as D32 steps 3–4 |
| R18 | Eq. (7.127): $N^2\equiv-\frac g{\rho_0}\frac{d\bar\rho}{dz}$ (= (1.29) without the adiabatic gradient) [#165] | 7.8 | B | RECAP | C15 | ch01 C51 (1.29), D18: reminded with `core.stratification.brunt_vaisala_sq` and ch01 `parcel_stability`; number: ocean thermocline N = 0.01 rad/s → period 10.5 min |
| N134 | Eqs. (7.128), (7.129), (7.130): perturbation momentum $\partial_tu=-p'_x/\rho_0$, $\partial_tv=-p'_y/\rho_0$, $\partial_tw=-p'_z/\rho_0-\rho'g/\rho_0$ [#166] | 7.8 | B | NOTE | C15 | stated as D32 steps 5–7 |
| N135 | Eq. (7.131): $\partial\rho'/\partial t-\frac{N^2\rho_0}gw=0$ [#167] | 7.8 | B | NOTE | C15 | stated as D32's result |
| R19 | The internal-wave system holds when ρ̄ does not depend on pressure; compressibility handled by potential density in N ('sigma-t') [#168] | 7.8 | C | RECAP | C15 | ch01 C55 (potential density): one sentence |
| N136 | Eq. (7.132): $\frac1{\rho_0}\nabla_H^2p'=\frac{\partial^2w}{\partial z\,\partial t}$ [#169] | 7.8 | B | NOTE | C15 | stated as D33 steps 1–4 |
| N137 | Eq. (7.133): $\frac1{\rho_0}\frac{\partial^2p'}{\partial t\,\partial z}=-\frac{\partial^2w}{\partial t^2}-N^2w$ [#170] | 7.8 | B | NOTE | C15 | stated as D33 steps 5–7 |
| N138 | Eq. (7.134): the w-equation $\frac{\partial^2}{\partial t^2}\nabla^2w+N^2\nabla_H^2w=0$ [#171] | 7.8 | B | NOTE | C15 | stated as D33's result (★★★, written out); `ch07.w_equation_residual` on a plane wave |
| N139 | Eq. (7.135): $\omega=\omega(k,l,m)=\omega(\mathbf K)$ — anisotropic [#172] | 7.8 | B | NOTE | C15 | stated as D34 step 1 |
| R20 | Internal waves are baroclinic ⇒ vorticity generation ⇒ rotational; no Laplace equation [#173] | 7.8 | C | RECAP | C15 | ch05 C04 (baroclinic torque): one sentence linking N107 |
| N140 | Eq. (7.136): plane-wave trial $w=w_0e^{i(\mathbf K\cdot\mathbf x-\omega t)}$ [#174] | 7.8 | B | NOTE | C15 | stated as D34 step 2 |
| N141 | Eq. (7.137): $\omega^2=\frac{k^2+l^2}{k^2+l^2+m^2}N^2$ [#175] | 7.8 | B | NOTE | C15 | stated as D34 steps 3–5 |
| N142 | Eq. (7.138): $\omega=kN/\sqrt{k^2+m^2}=kN/K$ (valid for k > 0; code uses ∣k∣, R9) [#176] | 7.8 | B | NOTE | C15 | stated as D34 step 6 |
| C15 | Eq. (7.139): internal-wave dispersion $\omega=N\cos\theta$ — frequency set by the direction of K (θ above the horizontal), not its size; 0 < ω < N [#177] | 7.8 | A | CORE | – | load-bearing: the defining property of internal gravity waves in ocean and atmosphere (Ch. 13 adds f: inertia–gravity waves; internal tides, lee waves); D32, D33, D34 |
| N143 | Eq. (7.140): $u=u_0e^{i(kx+ly+mz-\omega t)}$ [#178] | 7.8 | B | NOTE | C16 | stated as D35 step 1 |
| N144 | Eq. (7.141): $\mathbf K\cdot\mathbf u=0$ — particles move along the crests (transverse waves) [#179] | 7.8 | B | NOTE | C16 | stated as D35's result; `internal_wave_fields` asserts K·u = 0 |
| N145 | Limits: θ = 0 ⇒ ω = N (columns bobbing — ch01's parcel); θ = π/2 ⇒ ω = 0 needs the full equations [#180] | 7.8 | B | NOTE | C15 | stated with `core.stratification.parcel_displacement` parity at ω = N |
| N146 | Eq. (7.142): steady layered flow w = p′ = ρ′ = 0 with any horizontally nondivergent $\partial u/\partial x+\partial v/\partial y=0$ at each level (strong N) [#181] | 7.8 | B | NOTE | C15 | stated (a-D59 given, §4c) with `ch07.layered_flow_check`; climate hook: pancake eddies, flat cloud sheets (Ch. 13) |
| N147 | Blocking ahead of a cylinder in strongly stratified flow (Fig. 7.30) [#182] | 7.8 | C | NOTE | C15 | named with a sketch from our code (N196) and pointer Ch. 13 (orographic blocking) |
| N148 | Eq. (7.143): $\mathbf c_g=\mathbf e_x\partial\omega/\partial k+\mathbf e_y\partial\omega/\partial l+\mathbf e_z\partial\omega/\partial m$ [#183] | 7.8 | B | NOTE | C16 | stated as D36 step 1 (component form of N71) |
| N149 | Eqs. (7.144), (7.145): $\mathbf c=\frac\omega{K^2}(k\mathbf e_x+m\mathbf e_z)$, $\mathbf c_g=\frac{Nm}{K^3}(m\mathbf e_x-k\mathbf e_z)$ (k > 0) [#184] | 7.8 | B | NOTE | C16 | stated as D36 steps 2–6; `core.W.internal_wave_velocities` (sign-safe, R9) |
| C16 | Eq. (7.146): $\mathbf c_g\cdot\mathbf c=0$ — group velocity perpendicular to phase velocity: phase up ⇔ energy down, packets slide along their crests; St Andrew's cross beams at cos θ = ω/N; energy flux $\mathbf F=\mathbf c_gE$ (7.159) [#185] | 7.8 | A | CORE | – | load-bearing: the most counter-intuitive and most used internal-wave fact — beams from tides over topography, lee waves, energy radiation into the deep ocean and the stratosphere (Ch. 13); D35, D36, D37 |
| N150 | St Andrew's cross: a source oscillating at ω < N radiates four beams at angle θ to the vertical, cos θ = ω/N (Fig. 7.33; ω = 0.71 N → 45°) [#186] | 7.8 | B | NOTE | C16 | stated with `core.W.beam_angle` and `ch07.st_andrews_cross` (our field) — animation A6 and explainer E9; ⚠️ angle to the vertical vs to the horizontal (R10) |
| N151 | N(z) in the real ocean (< 0.01 rad/s): results hold locally if N varies slowly over 2π/m; WKB in Ch. 13 [#187] | 7.8 | C | NOTE | C16 | named with pointer Ch. 13 (and C10 rays) |
| N152 | Eq. (7.147): internal-wave energy equation $\partial_t[\tfrac12\rho_0\lvert\mathbf u\rvert^2]+g\rho'w+\nabla\cdot(p'\mathbf u)=0$ (cf. (4.56)) [#188] | 7.8 | B | NOTE | C16 | stated (a-D61 given, §4c) with `internal_energy_budget_residual` ~ 1e-8 |
| N153 | Eq. (7.148): $\partial E_p/\partial t=g\rho'w=\partial_t[g^2\rho'^2/(2\rho_0N^2)]$ [#189] | 7.8 | B | NOTE | C16 | stated (a-D62 given, §4c) |
| N154 | Eq. (7.149): $\rho'=N^2\rho_0\zeta/g$ [#190] | 7.8 | B | NOTE | C16 | stated in one line (integrate (7.131) with w = ∂ζ/∂t) |
| N155 | Eq. (7.150): $E_p=\frac{g^2\rho'^2}{2\rho_0N^2}=\tfrac12N^2\rho_0\zeta^2$ (available potential energy per volume) [#191] | 7.8 | B | NOTE | C16 | stated with a number (N = 0.01 rad/s, ζ = 10 m → 5 J/m³); Ch. 13 APE hook |
| N156 | Eq. (7.151): the two-deep-fluid value $\tfrac14(\rho_2-\rho_1)ga^2$ to be recovered [#192] | 7.8 | C | NOTE | C16 | named as the target of N157 (pointer back to N105) |
| N157 | Eq. (7.152): a density jump as $N^2=\frac g{\rho_0}(\rho_2-\rho_1)\delta(z)$ recovers (7.151) [#193] | 7.8 | B | NOTE | C16 | stated (a-D63 given, §4c) with `internal_pe_interface_limit` converging as the tanh width ε → 0 (ch06 P149 delta reminder) |
| N158 | Plane-wave ansatz $[u,w,p',\rho']=[\hat u,\hat w,\hat p,\hat\rho]e^{i(kx+mz-\omega t)}$ [#194] | 7.8 | B | NOTE | C16 | stated as D37 step 1 |
| N159 | Eq. (7.153): polarization $p'=-\frac{\omega m\rho_0}{k^2}\hat w e^{i\theta}$, $\rho'=\frac{iN^2\rho_0}{\omega g}\hat w e^{i\theta}$, $u=-\frac mk\hat w e^{i\theta}$ [#195] | 7.8 | B | NOTE | C16 | stated as D37 steps 2–5 (ρ′ 90° out of phase with w) |
| N160 | Eq. (7.154): $E_k=\tfrac14\rho_0(\frac{m^2}{k^2}+1)\hat w^2$ [#196] | 7.8 | B | NOTE | C16 | stated (a-D65 given, §4c) with the ½Re(AB*) primer |
| N161 | Eq. (7.155): $E_p=\frac{N^2\rho_0}{4\omega^2}\hat w^2$ [#197] | 7.8 | B | NOTE | C16 | stated (a-D65 given, §4c) |
| N162 | Eq. (7.156): equipartition $E_k=E_p$ for internal waves [#198] | 7.8 | B | NOTE | C16 | stated (a-D65 given): equal iff ω² = k²N²/K²; `internal_wave_energy` asserts it |
| N163 | Eq. (7.157): $E=\tfrac12\rho_0(\frac{m^2}{k^2}+1)\hat w^2$ [#199] | 7.8 | B | NOTE | C16 | stated in one line |
| N164 | Eq. (7.158): $\mathbf F=\overline{p'\mathbf u}=\frac{\rho_0\omega m\hat w^2}{2k^2}(\mathbf e_x\frac mk-\mathbf e_z)$ [#200] | 7.8 | B | NOTE | C16 | stated as D37 steps 6–9 |
| N165 | $\mathbf c_gE=\frac{Nm}{K^3}(m\mathbf e_x-k\mathbf e_z)[\frac{\rho_0}2(\frac{m^2}{k^2}+1)\hat w^2]$ (unnumbered) [#201] | 7.8 | B | NOTE | C16 | stated as D37 steps 10–11 |
| N166 | Eq. (7.159): $\mathbf F=\mathbf c_gE$ — the same law as (7.71), per unit area/volume here [#202] | 7.8 | B | NOTE | C16 | stated as D37's result; `internal_wave_energy` asserts F − c_gE = 0; normalisation switch (R5) in a callout |
| N167 | Fig. 7.1: crests of a plane wave; spacing along x, y exceeds λ = 2π/K; trace velocities [#F1] | 7.1 | B | NOTE | C01 | drawn by our code (`plane_wave`, `trace_velocities`) inside C01 |
| N168 | Fig. 7.2: geometry — z up, surface z = 0 with η(x, t), bottom z = −H [#F2] | 7.2 | B | NOTE | C02 | our matplotlib sketch at the top of C02 |
| N169 | Fig. 7.3: one particle's closed elliptical orbit, clockwise [#F3] | 7.2 | B | NOTE | C05 | drawn from `orbit_linear` with the excursion (ξ, ζ) marked |
| N170 | Fig. 7.4: orbits in deep (circles), intermediate and shallow (thin ellipses) water [#F4] | 7.2 | B | NOTE | C05 | three panels kH = 3, 1, 0.3 from `orbit_semi_axes` — C05's main static figure (and animation A1) |
| N171 | Fig. 7.5: instantaneous streamlines of a progressive wave [#F5] | 7.2 | B | NOTE | C05 | contours of (7.37) from `wave_fields` with ψ = 0 lines marked |
| N172 | Fig. 7.6: column-swap construction for E_p [#F6] | 7.2 | B | NOTE | C06 | our sketch inside D13 |
| N173 | Fig. 7.7: cosh, sinh, tanh on 0 ≤ x ≤ 2.3 [#F7] | 7.2 | B | NOTE | C04 | our figure with the kH = 2 marker (tanh 2 = 0.964) |
| N174 | Fig. 7.8: crests turning parallel to a straight beach [#F8] | 7.2 | B | NOTE | C10 | drawn from `ray_trace` on a plane beach with crest lines (C10) |
| N175 | Fig. 7.9: crest lines around a circular island [#F9] | 7.2 | B | NOTE | C10 | drawn from `ray_trace` over a circular H(r) (same cell as N48) |
| N176 | Fig. 7.10: c(λ) with capillary, deep and shallow branches and the minimum [#F10] | 7.3 | B | NOTE | C07 | log–log figure from `phase_speed` (σ = 0.0727 N/m, H = 1 m) — C07's main static figure; plotly slider IF2 |
| N177 | Fig. 7.11: standing-wave streamlines [#F11] | 7.4 | B | NOTE | C08 | contours of (7.62) at ωt = π/2 from `standing_wave_fields` |
| N178 | Fig. 7.12: u(x) for seiche modes n = 0 and n = 1 [#F12] | 7.4 | B | NOTE | C08 | drawn from `seiche_modes` (plotly slider IF7 over n) |
| N179 | Fig. 7.13: beats — carrier under the envelope, c and c_g arrows [#F13] | 7.5 | B | NOTE | C09 | drawn from `beat_wave` (Δω t form), with the printed Δω x form as a ghost |
| N180 | Fig. 7.14: ω(k) with chord (slope c) and tangent (slope c_g) [#F14] | 7.5 | B | NOTE | C09 | drawn from `omega_gravity` and `group_velocity` at one k |
| N181 | Fig. 7.15: packet of length ~1/δk and its spectrum of width δk [#F15] | 7.5 | B | NOTE | C09 | two panels from `linear_evolve` + FFT |
| N182 | Fig. 7.16: stone in a pond at t₁ < t₂ < t₃; front at c_g,max t, rear at c_g,min t [#F16] | 7.5 | B | NOTE | C09 | FFT capillary–gravity impulse (FAST grid) with the c_g,min line |
| N183 | Fig. 7.17: x–t diagram — crest lines dx/dt = c, rays dx/dt = c_g [#F17] | 7.5 | B | NOTE | C10 | crest tracking on the FFT packet (C10) |
| N184 | Fig. 7.18: curved rays of constant ω in an inhomogeneous medium [#F18] | 7.5 | B | NOTE | C10 | x–t rays over a slope from `ray_trace`, ω printed along each ray |
| N185 | Fig. 7.19: finite-amplitude profile steepening and overturning [#F19] | 7.6 | B | NOTE | C11 | frames from `simple_wave_evolve` up to the breaking time (animation A5 first part) |
| N186 | Fig. 7.20: hydraulic jumps — spillway, stationary CV, moving jump [#F20] | 7.6 | B | NOTE | C11 | our CV sketch + H₂/H₁ and head loss vs Fr₁ from `hydraulic_jump` |
| N187 | Fig. 7.21: Stokes-wave profile — peaked crests, flat troughs [#F21] | 7.6 | B | NOTE | C12 | `stokes_wave_profile` at ka = 0.3 vs the linear cosine |
| N188 | Fig. 7.22: open orbit and a dyed vertical line bending forward [#F22] | 7.6 | B | NOTE | C12 | exact path lines from `particle_path` (animation A4) |
| N189 | Fig. 7.23: cnoidal wave train and solitary wave [#F23] | 7.6 | B | NOTE | C11 | `solitary_wave` and a `kdv_solve` run (optional `cnoidal_wave`) |
| N190 | Fig. 7.24: interface wave — opposite horizontal velocities (vortex sheet) [#F24] | 7.7 | B | NOTE | C13 | quiver of `interface_fields` with the sheet marked |
| N191 | Fig. 7.25: column swap across a density interface [#F25] | 7.7 | C | NOTE | C13 | named beside N105 (the D13 sketch reused, no new figure) |
| N192 | Fig. 7.26: 'dead water' ship over fresh-over-salt layers [#F26] | 7.7 | C | NOTE | C13 | named only (N108) |
| N193 | Fig. 7.27: barotropic vs baroclinic modes of a layer over deep fluid [#F27] | 7.7 | B | NOTE | C14 | two panels from `two_layer_modes` (surface, interface, u arrows) — C14's main static figure |
| N194 | Fig. 7.28: shallow two-layer modes in the Boussinesq limit [#F28] | 7.7 | B | NOTE | C14 | u(z) profiles from `two_layer_modes` at small kH |
| N195 | Fig. 7.29: K, c up-left at θ; c_g down-left; particle motion along crests [#F29] | 7.8 | B | NOTE | C16 | our geometry figure from `internal_wave_velocities` drawn with k < 0 as in the book (R9) |
| N196 | Fig. 7.30: blocked layer ahead of a cylinder [#F30] | 7.8 | C | NOTE | C15 | named; a schematic from our code (streamlines stopping between the tangent planes) |
| N197 | Fig. 7.31: c and c_g as legs of a right triangle with horizontal hypotenuse [#F31] | 7.8 | B | NOTE | C16 | drawn from `internal_wave_velocities` for three ω/N |
| N198 | Fig. 7.32: circular packet at two times — phase lines up-left, packet down-left [#F32] | 7.8 | B | NOTE | C16 | two frames from `linear_evolve_2d` (FAST grid; animation A6) |
| N199 | Fig. 7.33: St Andrew's cross beams (ω = 0.71 N) [#F33] | 7.8 | B | NOTE | C16 | our field `st_andrews_cross` (not the photograph) with beam angles labelled to the vertical and to the horizontal |
| S01 | Exercises 7.1–7.20 (named only; text never reproduced) [#X1] | Ex. | C | SKIP | C09 | pointer: used ideas — 7.3 → N02; 7.4 → N30; 7.5–7.6 → N63; 7.9–7.10 → N69, N72; 7.11 → N72 (viscous decay); 7.12 → D20; 7.14 → N91; 7.15 → N97; 7.18 → N105; 7.19 → D30; 7.20 → N131; the rest left to the reader |

### 2a. What each A block contains (for the lesson-designer)
Numbers are ours (computed with `G_BOOK` = 9.81 m/s², ρ = 1000 kg/m³ unless stated); book-quoted values stay in
`tests/book_values_ch07.json`.
- **C01** picture: a buoy bobbing as swell passes a pier — the crest moves, the water mostly does not · question: "what
  exactly moves at the 'wave speed', and how do λ, T, k, ω and c fit together?" · D01 · number: λ = 100 m, T = 8 s →
  k = 0.0628 rad/m, ω = 0.785 rad/s, c = 12.5 m/s; plane wave k = l = 1 rad/m → λ = 4.44 m but 6.28 m between crests
  along x · code: `core.W.sinusoid`, `wave_parameters`, `crest_positions`, `plane_wave`, `trace_velocities`,
  `doppler_frequency` + from-scratch crest tracking (argmax of η on a fine grid, Δx/Δt vs ω/k) · figure: η(x) at two
  times with the tracked crest; Fig. 7.1 remake (crest lines, K arrow, trace spacings) · B/C: N01–N11, N167.
- **C02** picture: a floating leaf on a swell — it stays on the surface and rides up and down · question: "the surface
  is an unknown — how can we impose conditions on it?" · D02, D03, D04 · number: a = 0.5 m, λ = 50 m → ka = 0.063, so the
  dropped terms are ~6 % of the kept ones; at ka = 0.3 they are ~30 % · code: `ch07.free_surface_residuals` (exact vs
  linear residuals), `core.interfaces.kinematic_bc_residual` with f = z − η · figure: residual size vs ka on log–log axes
  (slope 2 for the relative error of the linearised conditions) + the Fig. 7.2 sketch · B/C: N12–N17, R01–R06, N168.
- **C03** picture: a wave tank — long waves outrun short ones · question: "given k and the depth, how fast must the
  surface oscillate?" · D05, D06 · number: λ = 156 m in deep water → T = 10.0 s; H = 10 m, λ = 50 m → kH = 1.257,
  ω = 1.024 rad/s, T = 6.14 s · code: `core.W.omega_gravity`, `wavenumber_from_omega`, `ch07.wave_fields`,
  `ch07.surface_wave_sympy` (+ parity with `ch04.linear_wave_fields`) + from-scratch Newton iteration for k given T vs
  `wavenumber_from_omega` · figure: ω(k) for H = 1, 10, 100 m with the √(gk) and k√(gH) asymptotes · B/C: N18–N25.
- **C04** picture: a tsunami crossing the Pacific in a day while swell takes a week to cross it · question: "when does
  depth matter, and why are long waves not dispersive?" · D07, D08, D09 · number: T = 10 s swell → λ = 156 m (deep)
  vs 137 m on a 30 m shelf; tsunami H = 4 km → c = 198 m/s; pressure 1 m wave, λ = 50 m, H = 10 m: 9.8 kPa at z = 0,
  5.2 kPa at the bottom · code: `core.W.phase_speed`, `depth_regime`, `wavelength_from_period`, `ch07.pressure_response`
  · figure: c/√(gH) vs kH with the two limits and the 2 %/3 % bands; pressure profiles for three kH · B/C: N26, N37–N39,
  N42, N43, N46, R07, N173.
- **C05** picture: a cork and a neutrally buoyant float under a passing wave · question: "what path does a water parcel
  follow?" · D10, D11 · number: a = 1 m, λ = 20 m in 50 m of water: orbit radius 1 m at the surface, 0.53 m at z = −2 m,
  0.04 m at −10 m · code: `ch07.orbit_linear`, `orbit_semi_axes`, `particle_path` (linear vs exact) + from-scratch RK4
  of (7.34) vs `orbit_linear` · figure: Fig. 7.4 remake (three regimes), streamlines (7.37); animation A1 · B/C:
  N27–N30, N40, N41, N44, N45, R08, N169–N171.
- **C06** picture: a surfer's ride — energy delivered to a beach from a storm far away · question: "how much energy does
  a wave hold, and how fast is it carried?" · D12, D13, D14 · number: a = 1 m → E = 4.9 kJ/m²; T = 10 s deep water:
  F ≈ 38 kW/m, half of E × c · code: `ch07.wave_energy` (closed vs `dblquad`), `ch07.energy_flux`,
  `core.W.wave_energy_density` + from-scratch midpoint double sum of ½ρ(u² + w²) vs `wave_energy(method='closed')` ·
  figure: KE and PE density over one wavelength (equal averages), F/(Ec) vs kH (from ½ to 1) · B/C: N31–N36, N172.
- **C07** picture: rain ripples ahead of a lake's wind waves · question: "why do the smallest waves on a pond go faster
  than slightly longer ones?" · D15, D16 · number: σ = 72.74 mN/m (IAPWS, 20 °C), ρ = 998.2 → c_min = 23.1 cm/s at
  λ_m = 1.71 cm; a = 1 mm, λ = 1 cm: +29 Pa under a crest · code: `core.W.omega_capillary_gravity`, `phase_speed`,
  `capillary_minimum`, `ch07.curvature`, `capillary_surface_pressure` + from-scratch scan for the minimum vs
  `capillary_minimum` · figure: Fig. 7.10 remake on log axes (IF2 slider over σ) · B/C: N49–N57, R09, N176.
- **C08** picture: water sloshing in a bathtub; Lake Geneva's seiche · question: "which waves can live in a closed
  basin?" · D17, D18 · number: bathtub L = 1.5 m, H = 0.2 m → T₀ = 2.2 s; lake L = 50 km, H = 100 m → T₀ ≈ 53 min ·
  code: `ch07.standing_wave_fields`, `seiche_modes`, `basin_modes` + from-scratch loop over n vs `seiche_modes` ·
  figure: Fig. 7.11 streamlines, Fig. 7.12 mode shapes; animation A3 · B/C: N58–N63, N177, N178.
- **C09** picture: swell arriving in groups ("sets"); a stone in a pond · question: "why do crests appear at the back of
  a group and vanish at the front, and what speed carries the energy?" · D19, D20, D21 · number: deep water λ = 156 m:
  c = 15.6 m/s, c_g = 7.8 m/s; a storm 1000 km away: energy arrives after 35.6 h; stone in a pond: calm inside
  c_g,min ≈ 17.8 cm/s × t (our minimisation) · code: `core.W.beat_wave`, `group_velocity`, `group_velocity_numeric`,
  `linear_evolve`, `envelope`, `min_group_velocity` + from-scratch central difference of ω(k) and envelope-peak
  tracking vs `group_velocity` · figure: Figs. 7.13–7.16 remakes; animation A2 · B/C: N64–N73, N179–N182.
- **C10** picture: waves always arriving parallel to the beach, whatever the wind direction · question: "how do k and ω
  change when the depth changes slowly?" · D22, D23, D24 · number: T = 8 s swell at 30° in 20 m of water on a 1:50
  slope → 11° by H = 2 m (Snell, k sin α = const) · code: `ch07.local_wavenumber_frequency`,
  `crest_conservation_residual`, `core.W.ray_trace` + from-scratch RK4 ray steps vs `ray_trace` · figure: x–t diagram
  (Fig. 7.17), rays over a slope and round an island (Figs. 7.8, 7.9, 7.18) · B/C: N47, N48, N74–N80, N174, N175,
  N183, N184.
- **C11** picture: the roller at the foot of a spillway; a tidal bore · question: "given the incoming depth and speed,
  how high is the jump and where does the energy go?" · D25, D26 · number: H₁ = 0.1 m, Fr₁ = 3 → H₂ = 0.377 m, head
  loss 0.141 m; solitary wave a = 0.2 m on H = 1 m → c = 3.44 m/s · code: `ch07.hydraulic_jump`,
  `jump_momentum_residual`, `simple_wave_evolve`, `ursell_number`, `solitary_wave`, `kdv_solve` (FAST) + from-scratch
  `np.roots` of r² + r − 2Fr₁² vs `hydraulic_jump` · figure: H₂/H₁ and head loss vs Fr₁ (Fr₁ < 1 greyed "would create
  energy"); animation A5 (steepening, then a KdV hump splitting) · B/C: N81–N85, N93–N97, R10, N185, N186, N189.
- **C12** picture: a swimmer who drifts shoreward although "waves only go up and down" · question: "if each orbit is a
  circle, why does floating debris travel with the waves?" · D27 · number: a = 1 m, T = 8 s deep water: ū_L = 4.9 cm/s at
  the surface, halving every 5.5 m of depth · code: `ch07.stokes_drift`, `stokes_drift_numeric`, `eulerian_mean_u`,
  `stokes_wave_profile`, `stokes_wave_speed` + from-scratch period average of the exact path line vs `stokes_drift` ·
  figure: ū_L(z) analytic vs numeric; Stokes profile vs cosine; animation A4 (dyed line bending) · B/C: N86–N92, N187,
  N188.
- **C13** picture: a fjord's fresh upper layer; a thermocline · question: "how does a wave on a weak density jump
  differ from a surface wave?" · D28 · number: Δρ = 2 kg/m³ on 1000 → ε = 0.032, λ = 100 m → T = 253 s instead of 8 s;
  equal energy needs a 22-times taller interfacial wave · code: `core.W.interface_omega`, `ch07.interface_fields`,
  `interface_residuals`, `interface_energy` + from-scratch `np.linalg.solve` for the two-sided constants vs
  `interface_fields` · figure: interface with u arrows opposite across it (Fig. 7.24) · B/C: N98–N108, R11, N190–N192.
- **C14** picture: the summer thermocline under a lake's surface; ocean surface vs thermocline tides · question: "what
  are the two ways a layer over a deep fluid can oscillate?" · D29, D30, D31 · number: H = 50 m, Δρ/ρ = 0.002 →
  internal c ≈ 0.99 m/s vs √(gH) = 22 m/s; a 10 m thermocline wave shows 2 cm at the surface · code:
  `core.W.two_layer_free_surface_omega`, `two_layer_long_wave_speed`, `reduced_gravity_book`, `ch07.two_layer_modes`,
  `two_layer_sympy` + from-scratch `np.roots` of (7.110) as a quadratic in ω² vs `two_layer_free_surface_omega` ·
  figure: both mode shapes (Fig. 7.27), long-wave profiles (Fig. 7.28), ω(k) of both modes with the (7.95) limit ·
  B/C: N109–N131, N193, N194.
- **C15** picture: waves on the thermocline / the atmosphere's lee waves; a dye line in a salt-stratified tank · question:
  "if N is the only frequency the fluid knows, what sets a wave's frequency?" · D32, D33, D34 · number: N = 0.01 rad/s
  (period 10.5 min); a wave at 45° has ω = 0.0071 rad/s (period 14.8 min); ω = N only for θ = 0 · code:
  `core.W.internal_wave_omega`, `beam_angle`, `ch07.boussinesq_linear_sympy`, `w_equation_residual`,
  `internal_wave_fields` + from-scratch finite-difference residual of (7.134) for a plane wave vs the dispersion
  relation · figure: ω/N vs θ (cos curve) with K arrows at three θ; polar plot showing ∣K∣ does not matter · B/C:
  N132–N142, N145–N147, R12–R20, N196.
- **C16** picture: the St Andrew's cross above an oscillating cylinder in a stratified tank · question: "why does energy
  leave at right angles to the way the crests move?" · D35, D36, D37 · number: ω/N = 0.71 → θ = 44.8° from the
  vertical; N = 1 rad/s, K = 1 rad/m at 45°: |c| = 0.71 m/s, |c_g| = 0.71 m/s, c·c_g = 0 · code:
  `core.W.internal_wave_velocities`, `group_velocity_vector`, `ch07.st_andrews_cross`, `internal_wave_energy`,
  `internal_energy_budget_residual`, `linear_evolve_2d` (FAST) + from-scratch central-difference gradient of ω(k, m)
  vs `internal_wave_velocities` and the dot product · figure: c and c_g triangles (Fig. 7.31), St Andrew's cross field
  (Fig. 7.33), packet at two times (Fig. 7.32); animation A6 · B/C: N143, N144, N148–N166, N195, N197–N199.

## 3. Section coverage
| § | Title | A | B | C | RECAP | SKIP |
|---|---|---|---|---|---|---|
| 7.1 | Introduction | C01 | N02, N03, N04, N05, N06, N07, N08, N09, N10, N11, N167 | N01 | — | — |
| 7.2 | Linear Liquid-Surface Gravity Waves | C02, C03, C04, C05, C06 | N13, N14, N15, N16, N17, N19, N20, N21, N22, N23, N24, N25, N26, N27, N28, N29, N30, N31, N32, N33, N34, N35, N36, N37, N38, N39, N40, N41, N42, N43, N44, N45, N46, N47, N48, N168, N169, N170, N171, N172, N173, N174, N175 | N12, N18 | R01 (C), R02 (B), R03 (B), R04 (B), R05 (B), R06 (B), R07 (C), R08 (B) | — |
| 7.3 | Influence of Surface Tension | C07 | N50, N51, N52, N53, N54, N55, N57, N176 | N49, N56 | R09 (B) | — |
| 7.4 | Standing Waves | C08 | N58, N59, N60, N61, N62, N63, N177, N178 | — | — | — |
| 7.5 | Group Velocity, Energy Flux, and Dispersion | C09, C10 | N65, N66, N67, N68, N69, N70, N71, N72, N73, N74, N75, N76, N77, N78, N79, N80, N179, N180, N181, N182, N183, N184 | N64 | — | — |
| 7.6 | Nonlinear Waves in Shallow and Deep Water | C11, C12 | N81, N82, N83, N84, N86, N87, N88, N89, N90, N91, N92, N93, N94, N95, N97, N185, N186, N187, N188, N189 | N85, N96 | R10 (B) | — |
| 7.7 | Waves on a Density Interface | C13, C14 | N99, N100, N101, N102, N103, N104, N105, N106, N107, N110, N111, N112, N113, N114, N115, N116, N117, N118, N119, N120, N121, N122, N123, N124, N125, N126, N127, N128, N129, N130, N131, N190, N193, N194 | N98, N108, N109, N191, N192 | R11 (C) | — |
| 7.8 | Internal Waves in a Continuously Stratified Fluid | C15, C16 | N132, N133, N134, N135, N136, N137, N138, N139, N140, N141, N142, N143, N144, N145, N146, N148, N149, N150, N152, N153, N154, N155, N157, N158, N159, N160, N161, N162, N163, N164, N165, N166, N195, N197, N198, N199 | N147, N151, N156, N196 | R12 (C), R13 (B), R14 (B), R15 (C), R16 (C), R17 (C), R18 (B), R19 (C), R20 (C) | S01 (end of chapter) |

Totals: A 16 · B 193 (182 NOTE + 11 RECAP) · C 27 (17 NOTE + 9 RECAP + 1 SKIP) = 236. No section is empty; every section has at least one A item. N47 and N48 (refraction, Figs. 7.8–7.9 = N174, N175) are §7.2 rows taught inside C10 (§7.5), where the rays are derived; the §7.2 notebook section names them with a forward pointer. S01 (exercises) is listed under §7.8 as the end of the chapter.

## 4. Prerequisites needing primers (concept or tool | needed by | why it is not A/B/RECAP)
Rows marked **gloss** are one plain sentence where the item is used; all others are full 📎 primers (one or two plain
sentences, what it means here, a 2–4-line runnable demo with easy numbers). New primers continue at **P165**. Items
already in `knowledge/primers.md` get a one-line reminder naming the earlier primer (listed at the end). Placement: the
hyperbolic-function and PDE-separation primers go **before D05** (top of C03); the complex-amplitude primer at the top of
§7.7 (before D28); the ½Re(AB*) primer before D37.

| Concept or tool | Needed by | Why it is not A/B/RECAP |
|---|---|---|
| phase of a wave: radians vs cycles, "the argument of the cosine", kx − ωt = const moving to the right | C01 (D01) | new vocabulary (ch01–ch06 used ω only as a vorticity or a rotation rate); also the ⚠️ ω symbol-change callout |
| boundary transfer: evaluating a condition at z = η by Taylor-expanding about z = 0, $f(\eta)=f(0)+\eta f'(0)+\dots$, and why both O(ka) corrections go together | C02 (D03, D04) | P26/P98 primed Taylor in x, t; using it to move a boundary is new |
| separation of variables for a PDE: a product trial f(z)·(function of x − ct) turns Laplace into an ODE in z; why the x-part is chosen by the boundary conditions | C03 (D05), C13, C14, C15 | P42 was separation in the ODE sense (dy/g(y) = f(x)dx); the PDE product trial is new |
| hyperbolic functions cosh, sinh, tanh, coth, sech: definitions from e^{±x}, cosh² − sinh² = 1, cosh² + sinh² = cosh 2x, sinh 2x = 2 sinh x cosh x, derivatives, small- and large-x limits (tanh 2 = 0.964) | C03 (D05, D06), C04 (D07–D09), C05, C06 (D12, D14), C09 (D21), C14 (D30), N97 | new maths tool used on every page of §7.2 (the analyst's primer); Fig. 7.7 remake is N173 |
| numerical overflow of cosh/sinh for large arguments and the stable ratio $e^{kz}(1+e^{-2k(z+H)})/(1-e^{-2kH})$; `np.tanh` saturates safely | C03, C04 (code) | **gloss** (P107 cancellation reminder) |
| matching coefficients: an identity A cos(kx − ωt) = B cos(kx − ωt) for all x, t means A = B | C03 (D06), C13 (D28) | **gloss** |
| inverting a monotonic function numerically: find k from ω with a bracket from the two limits | C03 (code, from-scratch Newton), C04 | P108 `brentq` reminder + a 3-line Newton iteration — **gloss** |
| sum-to-product identities $\cos A+\cos B=2\cos\frac{A-B}2\cos\frac{A+B}2$, $\cos A-\cos B=-2\sin\frac{A+B}2\sin\frac{A-B}2$ | C08 (D17), C09 (D19) | new maths tool |
| mean of a squared sinusoid over a period ⟨cos²⟩ = ⟨sin²⟩ = ½ and time average ⟨·⟩ vs wavelength average (overbar) | C06 (D12–D14), C12 (D27), C16 | P151 primed ∫ sin² over a period; the two averaging notations are new — **gloss** |
| minimum of a function: f′ = 0, f″ > 0; minimising c² instead of c gives the same point | C07 (D16) | P19 slope reminder — **gloss** |
| curvature of a plane curve $1/R=\eta_{xx}/(1+\eta_x^2)^{3/2}$ and its sign (centre of curvature below a crest) | C07 (D15) | new maths tool |
| `scipy.optimize.minimize_scalar` (bounded) | C07 (from-scratch cross-check), N72 | new Python tool — **gloss** with a 3-line demo |
| derivative as the limit of a difference quotient Δω/Δk → dω/dk; chord vs tangent slope | C09 (D19) | P19 reminder with the ω(k) picture — **gloss** |
| Fourier integral $\eta(x,t)=\int A(k)e^{i(kx-\omega(k)t)}dk$ and the Taylor expansion of ω(k) about k₀ inside it; change of variable | C09 (D20) | P142 primed discrete Fourier modes and the FFT; the continuous superposition and the "freeze the envelope" move are new |
| `np.fft.rfft`/`irfft`, `np.fft.rfftfreq`, periodic domains and wrap-around; evolving each mode with its own ω(k) | C09 (code), N02, N72 | P142 reminder + **gloss** for the evolution step |
| envelope of an oscillation via `scipy.signal.hilbert` (modulus of the analytic signal) | C09 (code, from-scratch envelope tracking) | new Python tool |
| complex-step derivative Im f(k + ih)/h | C09 (code: `group_velocity_numeric`) | **gloss** (why it has no cancellation error; P107 reminder) |
| a first-order wave equation $\partial q/\partial t+c\,\partial q/\partial x=0$ keeps q constant along dx/dt = c (characteristics) | C10 (D22, D23), N81 | new maths idea (ch03 material derivative is the analogy) |
| Snell's law: across slowly changing depth the along-shore wavenumber is conserved, k sin α = const | C10 (D24) | new physics vocabulary |
| equality of mixed partials for the phase θ(x, t) | C10 (D22) | P121 Schwarz reminder — **gloss** |
| quadratic formula and choosing the physical root; `np.roots` | C11 (D25, from-scratch), C14 (from-scratch) | **gloss** (P71 Vieta reminder) |
| sech² pulse and its derivatives (tanh′ = sech², sech′ = −sech tanh) | N97 (solitary wave, stated) | **gloss** |
| Lagrangian vs Eulerian mean: averaging following a particle vs at a fixed point | C12 (D27) | ch03 C01 recap sentence + **gloss** |
| `solve_ivp` with DOP853 and tight tolerances for long integrations (drift is O((ka)²) of the orbit) | C12 (code) | P31/P94 reminder — **gloss** |
| complex amplitudes: $\zeta=\mathrm{Re}\{a\,e^{i(kx-\omega t)}\}$, a complex a = size and phase shift, ∂/∂x → ik, ∂/∂t → −iω, and why Re may be dropped in linear algebra | C13 (D28), C14 (D29–D31), C15 (D34), C16 (D35–D37) | P45/P153 primed Euler's formula and complex numpy; complex *amplitudes* as a wave convention are new |
| linear systems with complex coefficients; sympy `solve`, `factor`, `simplify` with `I` | C14 (D29, D30) | P57/P40/P158 reminders — **gloss** |
| operator elimination for linear PDEs: apply ∂/∂t, ∇_H², ∂/∂z to whole equations and substitute; derivatives commute | C15 (D33) | new maths move |
| linearisation about a base state: products of two small quantities are dropped, small × O(1) are kept | C15 (D32) | P68 orders of smallness reminder — **gloss** |
| angle of a vector from its components, cos θ = ∣k∣/K, `np.arccos`, `np.degrees` | C15 (D34), C16 | **gloss** |
| gradient in wavenumber space ∂ω/∂K_i (vector of partial derivatives with respect to k, l, m) | C16 (D36), N71 | P25 partial derivative reminder — **gloss** |
| mean of a product of real parts $\langle\mathrm{Re}(Ae^{i\theta})\mathrm{Re}(Be^{i\theta})\rangle=\tfrac12\mathrm{Re}(AB^*)$ | C16 (D37), N160–N164 | new maths tool (P81 conjugate reminder) |
| 1-D Dirac delta as the limit of a narrow tanh step's derivative | N157 (stated) | P149 (2-D delta) reminder — **gloss** |

Reminders only (already primed): partial derivative P25 · first-order Taylor P26 · multivariable Taylor P98 · definite
integral P27 · `solve_ivp` P31/P94 · RK4 by hand P95 · sympy P40 · sympy expand/series P117 · linear second-order ODE
(e^{λz} trial) P44 · Euler's formula P45 · chain rule P49/P91 · orders of smallness P68 · level sets and the normal P75 ·
`meshgrid` P76 · broadcasting P77 · contour/quiver P78 · eigenvalue idea P80 · complex conjugate P81 · `quad`/`dblquad`
P87 · parametric curves P92 · `np.expm1` P107 · `brentq` P108 · momentum flux P114 · power of a force P127 · reduced
gravity P131 · moving level set P132 · Fourier modes and FFT P142 · 2-D delta P149 · integrals of sin/cos over a period
P151 · complex plane in numpy P153 · sympy complex P158 · `animate` P16 · `slider_figure` P17 · `show_viz` P18 · live
widgets P47.

## 4b. Derivations written out (parsed by tools: ID first, CORE id in a column, ★★★ for hard, explainer slugs backticked in the LAST column)
Written out because (a) the result belongs to an A item, or (b) the book never writes it out and the lesson needs it
(marked "(b) book never writes it out"). One small move per step; the book's skipped moves (analysis §2b) are filled in
and the notebook says so. Analysis §2b numbers in brackets (`a-Dnn`). 37 rows; 321 steps.

| ID | Result (Eq.) | CORE | Difficulty | Steps | Tools used | Traps | Shown in |
|---|---|---|---|---|---|---|---|
| D01 | crest speed $c=\omega/k=\lambda\nu$ (7.4) from the crest condition (7.3) [a-D01] | C01 | ★ | 5 | phase (primer), cos = 1 exactly at phase 2nπ, solving a linear equation | n labels the crest and drops out of the speed; ω/k = (2πν)/(2π/λ) = λν — mixing ν [Hz] and ω [rad/s] loses a 2π | notebook · `dispersion_relation` |
| D02 | exact kinematic condition (7.16) from (7.13)–(7.15) with f = z − η [a-D05] | C02 | ★★ | 8 | level sets and the normal (P75), moving level set (P132), chain rule (P91), recap ch04 D29 | η is now the elevation, not ch04's level-set function (R1): f = z − η, ∇f = (−η_x, 1); only the normal part of U_s matters; the surface must be a single-valued graph (no overturning) | notebook |
| D03 | linearised kinematic condition (7.17) → (7.18) [a-D06] | C02 | ★★ | 8 | boundary transfer (primer), orders of smallness (P68), Taylor (P26) | η_xφ_x ~ ka·aω and ηφ_zz ~ ka·aω are the same order — dropping one and keeping the other is inconsistent; the condition is then applied at z = 0, not at z = η; a/H ≪ 1 is also needed | notebook |
| D04 | linearised dynamic condition (7.19) → (7.20) → (7.21) [a-D07] | C02 | ★★ | 8 | unsteady Bernoulli (R06, ch04 D26), gauge pressure, boundary transfer (primer) | the Bernoulli "constant" C(t) is zero (fluid at rest far away, gauge p = 0) or absorbed into φ; gz at z = η gives gη (not 0); ½∣∇φ∣² is O(ka) relative to φ_t; sign: φ_t = −gη | notebook |
| D05 | the potential (7.26) from the separable trial (7.22) through (7.23)–(7.25) [a-D08] | C03 | ★★ | 11 | separation of variables for a PDE (primer), e^{λz} trial (P44), hyperbolic functions (primer), two linear equations | why sin (not cos): (7.18) and (7.21) each differentiate η once; +k² gives real exponentials, not sines; regroup Ae^{kz} + Ae^{−2kH}e^{−kz} = 2Ae^{−kH}cosh k(z + H); the dynamic condition is not used yet — it is saved for ω(k) | notebook · `dispersion_relation` |
| D06 | dispersion relation (7.28) and its T–λ form [a-D10] | C03 | ★★ | 7 | matching coefficients (gloss), tanh/coth (primer) | cancelling cos(kx − ωt) is allowed because the identity holds for all x, t; (ω²/k) coth kH = g ⇒ ω² = gk tanh kH (coth ↔ tanh flip); keep the positive root; T–λ form needs both 2π factors | notebook · `dispersion_relation` |
| D07 | phase speed (7.29) and why it grows with λ (dispersive) [a-D11] | C04 | ★ | 5 | (7.4), derivative of tanh x/x (primer) | √(gk tanh kH)/k = √((g/k) tanh kH), not √(g tanh kH/k²)·k; the λ form swaps k = 2π/λ inside the tanh too | notebook · `dispersion_relation` |
| D08 | deep-water limits (7.45)–(7.48): c, circles, e^{kz} velocity and pressure [a-D20] | C04 | ★ | 7 | large-x limits of tanh, cosh/sinh (primer) | kH > 2 gives 1.8 % (not 0 %) error in c; the ratio cosh k(z + H)/sinh kH → e^{kz} only for z + H ≫ 1/k; e^{−π} = 0.043 at z = −λ/2 | notebook · `dispersion_relation` |
| D09 | shallow-water limits (7.49)–(7.52): √(gH), flat ellipses, hydrostatic pressure [a-D21] | C04 | ★ | 7 | small-x expansions tanh x ≈ x − x³/3, cosh ≈ 1 (primer) | c/√(gH) ≈ 1 − (kH)²/6 (3 % at H = 0.07λ); the semi-major axis a/kH is large, the semi-minor ≤ a; p′ = ρgη is independent of z (hydrostatic) | notebook · `dispersion_relation` |
| D10 | linearised particle orbits (7.33) → (7.34a, b) → (7.35a, b) [a-D13] | C05 | ★★ | 9 | path lines (R08), multivariable Taylor (P98), antiderivatives of sin/cos | ∫cos(kx₀ − ωt)dt = −sin(kx₀ − ωt)/ω (the minus sign in ξ); constants of integration are zero by the *definition* of (x₀, z₀) as the mean position; the dropped term ξu_x is O(ka) and is exactly Stokes drift (C12) | notebook · `particle_orbits` |
| D11 | orbits are ellipses (7.36); focal distance constant; clockwise [a-D14] | C05 | ★ | 6 | sin² + cos² = 1, ellipse geometry (P104), cosh² − sinh² = 1 | divide by the semi-axes *before* squaring; clockwise for a right-going wave (at the top ζ = +B, ξ̇ > 0); all particles of a column are in phase | notebook · `particle_orbits` |
| D12 | kinetic energy (7.38) → (7.39): $E_k=\tfrac12\rho g\overline{\eta^2}$ [a-D16] | C06 | ★★ | 10 | ⟨cos²⟩ (gloss, P151), ∫cosh², ∫sinh² (primer), (7.28) | integrate to z = 0, not η (the slab is a cubic correction); the ±H/2 terms cancel between cos² and sin² weights; sinh 2kH = 2 sinh kH cosh kH; ω² = gk tanh kH turns the result into gravity | notebook |
| D13 | potential energy (7.40) → equipartition (7.41) → $E=\tfrac12\rho ga^2$ (7.42) [a-D17 + a-D18] | C06 | ★ | 7 | splitting an integral, PE ρgz per volume, ⟨cos²⟩ | $\int_{-H}^{\eta}z\,dz-\int_{-H}^{0}z\,dz=\eta^2/2$ (the H² terms cancel); $\overline{\eta^2}=a^2/2$ only for a sinusoid; equipartition fails with rotation (Ch. 13) | notebook |
| D14 | energy flux (7.43) → (7.44) [a-D19] | C06 | ★★ | 10 | power of a pressure force (P127), time average, hyperbolic identities | the background term vanishes because ⟨u⟩ = 0 at every z (so pulling ⟨u⟩ out is legal); use the first form of (7.31); "u from (7.28)" means (7.27) (T2); ω³ = ω·gk tanh kH | notebook |
| D15 | surface-tension condition (7.53) → (7.55) and the capillary–gravity dispersion (7.56)–(7.57) [a-D23 + a-D24] | C07 | ★★ | 10 | curvature (primer), Laplace jump (R09), gauge pressure | under a crest η_xx < 0 and the liquid pressure is *higher* (p = −ση_xx > 0); only one radius because ∂/∂y = 0; η_xx = −k²η turns σ into the replacement g → g + σk²/ρ | notebook · `capillary_gravity_waves` |
| D16 | minimum phase speed (7.58) and its numbers (7.59) [a-D25] | C07 | ★★ | 7 | minimum of a function (gloss), (7.57) deep | minimise c² = g/k + σk/ρ (same point as c); k_m = √(ρg/σ) ⇒ λ_m = 2π√(σ/ρg); c_min² = 2√(gσ/ρ); at the minimum c_g = c (tangent through the origin) | notebook · `capillary_gravity_waves` |
| D17 | standing wave 2a cos kx cos ωt, its ψ (7.62) and u (7.63) [a-D27] | C08 | ★★ | 8 | sum-to-product (primer), ψ (7.37) (N30) | the left-going wave's ψ carries a minus sign (its u reverses); cos A − cos B = −2 sin((A+B)/2) sin((A−B)/2); nodes of η are antinodes of u | notebook · `seiche_standing_waves` |
| D18 | seiche wavelengths (7.64) and frequencies (7.65) [a-D28] | C08 | ★ | 5 | zeros of sine, (7.28) | the wall at x = 0 is satisfied automatically by sin kx (the book's origin choice); n starts at 0 because k = 0 is no wave; k = (n+1)π/L goes into (7.28) — shallow basins give T = 2L/((n+1)√(gH)) | notebook · `seiche_standing_waves` |
| D19 | beats (7.66) and $c_g=d\omega/dk$ (7.67) [a-D29] | C09 | ★ | 7 | sum-to-product (primer), difference quotient (gloss) | the page prints ½Δω **x** — read ½Δω **t** (T1); the envelope's wavelength is 4π/Δk (two beats per envelope period), its speed Δω/Δk | notebook · `group_velocity_packets` |
| D20 | a narrow packet's envelope moves at c_g: $\eta=a(x-c_gt)\cos(kx-\omega t)$ (7.68) — (b) book never writes it out (cited to Phillips) [a-D30] | C09 | ★★★ | 12 | Fourier integral (primer), Taylor of ω(k) (P26), change of variable, sympy check cell (chirped Gaussian for quadratic ω) | the carrier e^{i(k₀x − ω₀t)} still moves at c; the linear Taylor term is what shifts the envelope; the dropped ½ω″(k − k₀)² term spreads the packet — valid only for t ≪ 1/(ω″δk²); A(k) narrow is essential | notebook · `group_velocity_packets` |
| D21 | $c_g=\frac c2[1+\frac{2kH}{\sinh2kH}]$ (7.69), its limits (7.70) and $F=Ec_g$ (7.71) [a-D31 + a-D32 + a-D34] | C09 | ★★ | 11 | chain and product rules, hyperbolic identities (primer), limits of 2x/sinh 2x | dω/dk of √(gk tanh kH) needs the chain rule *and* the product rule; sech²/tanh = 2/sinh 2kH; deep c_g = c/2 (not 2c); the second factor of (7.44) is exactly (7.69) | notebook · `group_velocity_packets` |
| D22 | local k, ω (7.72)–(7.73), crest conservation (7.74) and k carried at c_g (7.75) [a-D35] | C10 | ★★ | 8 | mixed partials (P121), chain rule along a curve (P91), characteristics (primer) | ω(x, t) = ω(k(x, t)) is the *local* dispersion assumption (slowly varying train); k is constant along dx/dt = c_g, not along dx/dt = c | notebook · `wave_rays_refraction` |
| D23 | ω constant along rays in a slowly varying medium (7.76) → (7.79) [a-D36] | C10 | ★★ | 8 | chain rule, partial derivative at fixed x | at fixed x the medium is steady, so ω changes in time only through k; multiply (7.74) by c_g — not by c; the companion dk/dt = −∂ω/∂x (ours) is what bends rays | notebook · `wave_rays_refraction` |
| D24 | refraction on a plane beach: k sin α = const (Snell), crests turn parallel to the contours — (b) book never writes it out [a-D22] | C10 | ★★ | 9 | ray equations dx/dt = ∂ω/∂k, dk/dt = −∂ω/∂x (ours), Snell (primer), (7.28) | ω, not k, is conserved along a ray; k_y is conserved because H does not depend on y; α is measured from the shore normal; ∣k∣ grows as H falls, so α → 0 | notebook · `wave_rays_refraction` |
| D25 | Bélanger relation (7.80) → (7.81) [a-D38] | C11 | ★★ | 10 | CV momentum (ch04 C04 recap), hydrostatic face force ½ρgH², quadratic formula (gloss) | outflow − inflow of momentum = net pressure force (no bottom friction, uniform u); cancel H₁ − H₂ (≠ 0 for a jump); divide by gH₁³ to get 2Fr₁² = r(1 + r); keep the positive root | notebook · `hydraulic_jump` |
| D26 | jump energy change $E_2-E_1=-\frac{g(H_2-H_1)^3}{4H_1H_2}$ and the second-law choice [a-D39] | C11 | ★★ | 8 | factorisation, Q² from D25 | E = u²/2 + gH is per unit mass of a surface particle; the text file drops the minus sign (read the page); H₂ < H₁ would create mechanical energy — forbidden, so Fr₁ > 1 | notebook · `hydraulic_jump` |
| D27 | deep-water Stokes drift (7.84a) → (7.85) [a-D41] | C12 | ★★ | 9 | multivariable Taylor (P98), ⟨sin²⟩ = ⟨cos²⟩ = ½, Lagrangian vs Eulerian mean (gloss) | ⟨u(x₀, z₀, t)⟩ = 0; ξu_x ∝ sin² and ζu_z ∝ cos² add to a constant — no ½ in the answer; decay rate 2k (twice the orbits') | notebook · `particle_orbits` |
| D28 | interfacial dispersion (7.89)–(7.94) → (7.95) [a-D47] | C13 | ★★ | 10 | complex amplitudes (primer), decaying exponentials, matching coefficients | φ₁ must decay upward (e^{−kz}) and φ₂ downward (e^{kz}); A = −B = iωa/k from the two-sided kinematic condition; ∂/∂t → −iω in the pressure condition; check ρ₁/ρ₂ → 0 gives (7.45) | notebook · `two_layer_modes` |
| D29 | the two-layer constants (7.104)–(7.109) from (7.97)–(7.103) [a-D50] | C14 | ★★★ | 12 | complex amplitudes (primer), complex 2×2 systems (gloss), sympy check cell | (7.105) prints e^{i(kz−ωt)} — read kx (T3); conditions at z = −H bring e^{±kH}; b is complex in general (a phase difference to be found); C = A − Be^{2kH} | notebook · `two_layer_modes` |
| D30 | two-layer dispersion relation (7.110) (Exercise 7.19) [a-D51] | C14 | ★★★ | 12 | sympy `factor` (P158), hyperbolic identities (primer) | the pressure condition carries a common factor ag²k/ω² that must be divided out; expanding e^{±kH} into sinh/cosh; the factor (ω²/gk − 1) is the surface mode — do not cancel it as "trivial" | notebook · `two_layer_modes` |
| D31 | roots and mode shapes (7.111)–(7.114), long waves (7.115)–(7.119): $c=\sqrt{g'H}$ [a-D52 + a-D53] | C14 | ★★ | 11 | small-x limits (primer), substitution into (7.109) | g′ uses ρ₂ in the denominator here (ch04: ρ₁; R7); the baroclinic η/ζ is *negative* (antiphase); kH → ∞ of (7.113) must return (7.95) | notebook · `two_layer_modes` |
| D32 | linear perturbation equations (7.123)–(7.131) from the Boussinesq set [a-D54 + a-D55] | C15 | ★★ | 10 | linearisation (gloss), hydrostatic base (R16), N² (R18) | ρ̄ = ρ̄(z) kills three terms; drop u·∇ρ′ (small × small) but keep w dρ̄/dz (small × O(1)); subtracting the base state cancels −dp̄/dz − ρ̄g | notebook |
| D33 | the w-equation (7.132) → (7.134) [a-D56] | C15 | ★★★ | 12 | operator elimination (primer), mixed partials (P121), sympy check cell | ∂/∂t of continuity with (7.128)–(7.129) gives ∇_H²p′; eliminate ρ′ with ∂/∂t of (7.130) and (7.131); ∇_H² commutes with ∂/∂t and ∂/∂z (and with N(z)); regroup w_ttzz + ∇_H²w_tt = ∂²/∂t² ∇²w | notebook · `internal_wave_beams` |
| D34 | internal-wave dispersion (7.136) → (7.137) → $\omega=N\cos\theta$ (7.139) [a-D57] | C15 | ★★ | 8 | complex amplitudes (∂ → ik, −iω), angle from components (gloss) | ∇² → −K², ∇_H² → −(k² + l²); the book's ω = kN/K assumes k > 0 — use ∣k∣ (R9); θ is K's angle to the **horizontal** (= the beam's angle to the vertical) | notebook · `internal_wave_beams` |
| D35 | $\mathbf K\cdot\mathbf u=0$ (7.141): transverse waves [a-D58] | C16 | ★ | 5 | complex amplitudes, dot product | holds for any incompressible plane wave; surface waves are not transverse because they are not plane waves in z (they decay) | notebook · `internal_wave_beams` |
| D36 | $\mathbf c$ and $\mathbf c_g$ (7.143)–(7.145) and $\mathbf c_g\cdot\mathbf c=0$ (7.146) [a-D60] | C16 | ★★ | 9 | gradient in wavenumber space (gloss), quotient rule | ∂/∂k(kN/K) = Nm²/K³, ∂/∂m = −Nkm/K³; horizontal components share a sign, vertical ones are opposite; for k < 0 (Fig. 7.29) use ∇_K ω of N∣k∣/K (R9) | notebook · `internal_wave_beams` |
| D37 | polarization (7.153), mean flux (7.158) and $\mathbf F=\mathbf c_gE$ (7.159) [a-D64 + a-D66] | C16 | ★★ | 12 | complex amplitudes, ⟨Re·Re⟩ = ½Re(AB*) (primer), E (7.157) (N163) | take real parts before multiplying; ρ′ is 90° out of phase with w (factor i); ⟨p′u⟩ uses û* not û; use ω = kN/K at the end to see c_gE = F | notebook |

## 4c. Derivations demoted to statements (first column is the A parent in bold, e.g. **C20** — never a bare ID or a D id: the parser reads a bare first-cell ID as an item and blanks its tier)
Results **given, not derived**, inside the named A block (a paragraph, the equation, a number, a check in code). 24 rows.

| A parent | Analysis §2b item | Result stated (Eq.) | Stated in (B item) | Why not written out |
|---|---|---|---|---|
| **C01** | a-D02 λ = 2π/K (7.7), c = (ω/K)e_K (7.8), trace velocities | $\lambda=2\pi/K$, $c_x=\omega/k\ge c$ | N08, N09, N10 | level-set geometry of K·x = const in two sentences and a figure (Fig. 7.1, N167); checked with `trace_velocities` |
| **C01** | a-D03 Doppler shift (7.9) | $\omega_0=\omega+\mathbf U\cdot\mathbf K$ | N11 | one substitution x′ = x − Ut (ch03 P96 frames); used as a given later (Ch. 13) |
| **C01** | a-D67 general linear solution (Exercise 7.3) | $\eta=\int[\hat\eta_+e^{-i\omega t}+\hat\eta_-e^{i\omega t}]e^{ikx}dk$ | N02 | an exercise; the idea is carried by `linear_evolve` and by D20's Fourier integral |
| **C02** | a-D04 Laplace (7.11) from (7.10) | $\nabla^2\phi=0$ | R02 | SEEN: ch06 D01/C02 |
| **C03** | a-D09 velocities (7.27) | $u=a\omega\frac{\cosh k(z+H)}{\sinh kH}\cos(kx-\omega t)$, $w=a\omega\frac{\sinh k(z+H)}{\sinh kH}\sin(kx-\omega t)$ | N24 | one derivative each; parity with ch04's `linear_wave_fields` |
| **C04** | a-D12 wave pressure (7.31) | $p'=-\rho\phi_t=\rho ga\frac{\cosh k(z+H)}{\cosh kH}\cos(kx-\omega t)$ | N26 | two substitutions ((7.20) then (7.28)); the limits carry the teaching and are derived in D08–D09 |
| **C05** | a-D15 stream function (7.37) (Exercise 7.4) | $\psi=\frac{a\omega}{k}\frac{\sinh k(z+H)}{\sinh kH}\cos(kx-\omega t)$ | N30 | one integration in z; checked against (7.27) with `velocity_from_streamfunction_2d` |
| **C07** | a-D26 pure capillary speed (7.60) | $c=\sqrt{2\pi\sigma/(\rho\lambda)}$ | N57 | a one-line limit of (7.57) |
| **C09** | a-D33 capillary c_g = 3c/2 and c_g,min ≈ 17.8 cm/s (Exercises 7.9–7.10) | $c_g=\tfrac32c$; $\sigma k^2/\rho g=2/\sqrt3-1$ at the minimum | N69, N72 | exercises; checked by `group_velocity(g=0)` and a numerical minimisation vs `min_group_velocity` |
| **C09** | a-D68 viscous decay (Exercise 7.11) | $a=a_0e^{-2\nu k^2t}$ | N72 | needs Ch. 8's dissipation machinery; named with a pointer and `viscous_decay` |
| **C11** | a-D37 nonlinear steepening (our Riemann simple wave) | $c=3\sqrt{g(H+\eta)}-2\sqrt{gH}$, breaking time $t_b=-1/\min\partial_xc$ | N81 | the book is qualitative; characteristics belong to Ch. 15; shown as a labelled extension with the animation |
| **C11** | a-D44 linearised KdV phase speed | $c=c_0(1-\tfrac16k^2H^2)$ | N94 | two Taylor terms of (7.29); checked on log axes |
| **C11** | a-D45 Ursell ratio | $a\lambda^2/H^3$ | N95 | a scale estimate in one line |
| **C11** | a-D46 solitary wave solves KdV (Exercise 7.15) | $\eta=a\,\mathrm{sech}^2[(3a/4H^3)^{1/2}(x-ct)]$, $c=c_0(1+a/2H)$ | N97 | ★★★ algebra of a B item; `kdv_residual_sympy()` = 0 is shown as the check |
| **C12** | a-D40 Stokes expansion (7.82)–(7.83) | second harmonic ½ka², $c^2=\frac gk(1+k^2a^2)$ | N86, N87 | ★★★ perturbation theory the book only quotes; a sympy cell checks the O(ka) coefficient ½ |
| **C12** | a-D42 Stokes drift at any depth (7.86) (Exercise 7.14) | $\bar u_L=a^2\omega k\frac{\cosh2k(z_0+H)}{2\sinh^2kH}$ | N91 | D27's moves with cosh/sinh in place of e^{kz}; checked against exact path lines |
| **C12** | a-D43 zero Eulerian mean | $\bar u(x,z)=0$ at points always submerged | N92 | one integral of ∂w/∂x over a period; `eulerian_mean_u` = 0 |
| **C13** | a-D48 interfacial energies (Exercise 7.18) | $E_k=E_p=\tfrac14(\rho_2-\rho_1)ga^2$ | N105 | exponential integrals repeating D12–D13; `interface_energy` closed vs quad |
| **C13** | a-D49 vortex sheet | $u_2-u_1=2\omega a\cos(kx-\omega t)$ at z = 0 | N107 | one derivative of each φ; ch05 D23 gives the sheet strength |
| **C15** | a-D59 steady layered flow (7.142) | $w=p'=\rho'=0$, $\partial u/\partial x+\partial v/\partial y=0$ | N146 | a check by substitution; `layered_flow_check` |
| **C16** | a-D61 internal-wave energy equation (7.147) | $\partial_t[\tfrac12\rho_0\lvert\mathbf u\rvert^2]+g\rho'w+\nabla\cdot(p'\mathbf u)=0$ | N152 | the ch04 D21 move (u · momentum) repeated; `internal_energy_budget_residual` |
| **C16** | a-D62 potential energy (7.148)–(7.150) | $E_p=\frac{g^2\rho'^2}{2\rho_0N^2}=\tfrac12N^2\rho_0\zeta^2$ | N153, N154, N155 | one chain-rule line; stated with a number |
| **C16** | a-D63 interface limit with δ-function N² (7.151)–(7.152) | $\int\tfrac12N^2\rho_0\zeta^2dz\to\tfrac14(\rho_2-\rho_1)ga^2$ | N157 | a consistency check; shown numerically by `internal_pe_interface_limit` converging in ε |
| **C16** | a-D65 equipartition (7.154)–(7.157) | $E_k=\tfrac14\rho_0(\frac{m^2}{k^2}+1)\hat w^2=E_p$ | N160, N161, N162, N163 | the ½Re(AB*) move of D37 applied twice; asserted in `internal_wave_energy` |

## 5. Interactive explainers (5–10 + backup)
Nine explainers on the A items where manipulation or motion teaches most; one backup. Each has the required live
**Explain** tab ("Explanation & interpretation", numbered sections computing every number on screen with the reader's
settings, ending in "Reading the current setting"), a synced **Code** tab, a **Derivation** tab for its D ids, a
4–8-step walkthrough, ≥ 3 check questions and ≥ 2 more depth features. Physics mirrors scalar-callable `fluidpy`
functions (parity rows). **Colours** (consistent across explainers, notebook figures and derivation terms): the free
surface blue, particles/orbits teal, phase (crests, c) orange, group/energy (envelope, c_g, F) purple (`accent`),
pressure amber, surface tension rose, gravity blue, lower/denser layer navy-blue fill, upper layer light, the
deep-water limit dashed muted, the shallow-water limit dotted muted. Linked, not duplicated: ch01 `parcel_stability`
(from E9: the ω = N limit), ch04 `control_volume_budgets` (from E7: the bore scene), ch04 `which_bernoulli` (from E1:
unsteady Bernoulli), ch05 `vortex_sheet_rollup` (from E8: the interface as a vortex sheet).

### E1 · dispersion_relation
- **A:** C03, C04 (also shows C01's c = ω/k, N26 (7.31), N38 (7.45), N42 (7.48), N43 (7.49), N46 (7.52), N39 ocean
  numbers, N37 tanh limits)
- **Confusion removed:** "is the wave speed a property of the water or of the wave?" — both: it is set by the
  wavelength *and* the depth, and depth only matters once the wave is longer than about twice the depth.
- **Why interactive:** dragging λ across kH ≈ 0.3…3 moves the dot along ω(k) from the √(gH) line to the √(g/k) curve
  while the tank's wave visibly changes speed and the pressure profile stops reaching the bottom — three linked
  representations of one number that a static figure shows only one at a time.
- **Stage:** (1) the tank: animated surface on depth H with a crest marker moving at c, pressure shading
  cosh k(z + H)/cosh kH under it; (2) c(λ) (or ω(k)) curve with deep and shallow asymptotes, the moving dot, 2 %/3 %
  bands; (3) (hidePortrait) pressure amplitude p′/ρga vs z at the current kH with the hydrostatic (shallow) and e^{kz}
  (deep) ghosts.
- **Controls:** wavelength λ (log slider 1 m … 100 km) · depth H (log 0.5 m … 5 km) · amplitude a (visual only) ·
  view chips: c(λ) / ω(k) · transport (time).
- **Presets:** wind swell T = 10 s in the open ocean · the same swell on a 10 m beach · tsunami (λ = 200 km, H = 4 km) ·
  tide · lab tank (H = 0.5 m, λ = 1 m) · kH = 2 (deep threshold) · H = 0.07λ (shallow threshold).
- **Equations shown (live):** (7.28) $\omega=\sqrt{gk\tanh kH}$; (7.29) $c=\sqrt{\frac gk\tanh kH}$; (7.45)
  $c=\sqrt{g/k}$; (7.49) $c=\sqrt{gH}$; (7.31) $p'=\rho ga\frac{\cosh k(z+H)}{\cosh kH}\cos(kx-\omega t)$; (7.4) $c=\omega/k$.
- **Mirrors:** `core.W.omega_gravity`, `phase_speed`, `wavelength_from_period`, `depth_regime`, `ch07.pressure_response`.
- **Derivations:** D01 (crest speed), D05 (φ: the step "cosh k(z + H)" sets the profile view), D06 (ω(k): the
  "cancel the cosine" step lights the dot), D07, D08 (deep: preset kH = 3), D09 (shallow: preset kH = 0.2).
- **Depth features:** Explain tab (kH → tanh kH → ω → T → c; regime and error of each limit; bottom pressure fraction;
  time for the wave to cross 1000 km → reading the current setting: deep / intermediate / shallow), synced Code tab, +
  linked views (3), transport, presets, status ("🌊 deep: c depends on λ, depth invisible" / "↔ intermediate" / "🏖️
  shallow: c = √(gH), all λ together"), inspector (click the tank at a depth: p′ arithmetic).
- **Follows reference:** `angular_frequency_explorer_1.html` (system animation + curves on one clock, presets, "Right
  now" notes with a highlighted table of real waves).
- **Aha:** depth is invisible to a wave shorter than about twice the depth; for longer waves it sets the speed, and all
  long waves travel together at √(gH).

### E2 · particle_orbits
- **A:** C05, C12 (also shows N30 (7.37) streamlines, N40 (7.46) circles, N44 (7.50) flat ellipses, N89–N92 drift and
  the zero Eulerian mean, N86 Stokes profile)
- **Confusion removed:** "the water moves with the wave" vs "the water only goes up and down" — neither: parcels go
  round closed ellipses (circles in deep water, flat ellipses in shallow water) *to first order*, and at second order
  they creep forward (Stokes drift) although a fixed current meter measures zero mean.
- **Why interactive:** the reader watches tracer particles loop under the moving surface, drags the depth to squash
  circles into ellipses, then switches from "linear" to "exact" path lines and sees the loops open and a dyed column
  lean forward — motion and a model toggle are the whole lesson.
- **Stage:** (1) the water column: surface, a grid of tracer particles with orbit ghosts, a dyed vertical line,
  optional streamlines ψ; (2) orbit semi-axes A(z), B(z) vs depth (and the drift ū_L(z) in drift mode); (3)
  (hidePortrait) x-position of the selected particle vs time: linear (closed) vs exact (drifting) with the drift slope.
- **Controls:** kH (or λ with H) · steepness ka (0.01…0.3) · model chips: linear / exact · show streamlines toggle ·
  click a particle to select it.
- **Presets:** deep water (kH = 3) · intermediate (kH = 1) · shallow (kH = 0.3) · steep deep wave (ka = 0.25, exact) ·
  bottom particle (flat line orbit).
- **Equations shown (live):** (7.35a, b) excursions; (7.36) the ellipse; (7.46) $\xi,\zeta\propto ae^{kz_0}$; (7.50);
  (7.37) ψ; (7.84a) Taylor path line; (7.85) $\bar u_L=a^2\omega ke^{2kz_0}$; (7.86).
- **Mirrors:** `ch07.orbit_linear`, `orbit_semi_axes`, `particle_path` (RK4 in JS vs DOP853 parity at 1e-6), `stokes_drift`.
- **Derivations:** D10 (linearised orbits: the "freeze the argument" step shows the neglected ξu_x as an arrow), D11
  (ellipse), D27 (Stokes drift: the sin² + cos² step sets the model to exact).
- **Depth features:** Explain tab (semi-axes at the selected depth → focal distance → sense of rotation → drift per
  period → Eulerian mean 0 → reading the current setting), synced Code tab, + linked views (3), transport, modes
  (linear / exact), presets, inspector (click a particle: its A, B, ū_L arithmetic), end-of-run card (distance drifted
  in N periods).
- **Follows reference:** `forced_damped_vibrations.html` (phase-plane-like orbit view + time graph with toggled
  components and the explanation panel).
- **Aha:** each parcel draws a small loop that shrinks with depth — circles when deep, flat ellipses when shallow — and
  the loops fail to close by a hair, which is the Stokes drift.

### E3 · capillary_gravity_waves
- **A:** C07 (also shows N50–N55, N57, N72; links forward to C09's c_g,min)
- **Confusion removed:** "shorter waves are always slower" — only while gravity dominates; below about 1.7 cm surface
  tension takes over, speeds rise again, and no wave on water can travel slower than 23 cm/s.
- **Why interactive:** sweeping λ over four decades on log axes and changing σ (clean water, soapy water, mercury)
  moves the minimum and the crossover; the two restoring forces are shown as live term bars whose ratio decides the
  branch.
- **Stage:** (1) a magnified surface patch with a wave of the chosen λ, arrows for the gravity (blue) and surface
  tension (rose) restoring pressures at a crest; (2) c(λ) on log–log axes with capillary, gravity and shallow branches,
  the minimum marked and the current dot; (3) (hidePortrait) term bars g/k vs σk/ρ.
- **Controls:** wavelength (log 1 mm … 10 m) · σ (0 … 0.5 N/m) · liquid chips (water 20 °C, soapy water, mercury,
  ethanol) · depth H.
- **Presets:** λ = λ_m (the minimum) · raindrop ripple 5 mm · pure capillary 1 mm · wind wave 1 m · σ = 0 ghost.
- **Equations shown (live):** (7.54) $(p)_{z=\eta}=-\sigma\eta_{xx}$; (7.56); (7.57); (7.58) $c_{min}=(4g\sigma/\rho)^{1/4}$,
  $\lambda_m=2\pi\sqrt{\sigma/\rho g}$; (7.60).
- **Mirrors:** `core.W.omega_capillary_gravity`, `phase_speed`, `capillary_minimum`, `min_group_velocity`,
  `ch07.capillary_surface_pressure`.
- **Derivations:** D15 (the curvature step highlights the rose arrows), D16 (minimum: the "d(c²)/dk = 0" step moves
  the dot to λ_m).
- **Depth features:** Explain tab (restoring terms g/k and σk/ρ with numbers → c → which dominates → c_min, λ_m for
  the liquid → reading the current setting: capillary / gravity / crossover), synced Code tab, + linked views (3),
  presets, term bars, status ("💧 capillary ripple" / "🌊 gravity wave" / "⚖️ near the minimum"), modes (liquids).
- **Follows reference:** `overfitting_curves.html` (a minimal two-slider curve with an optimum marker and a regime
  verdict) with `fid_formula_lab.html`'s term bars.
- **Aha:** two restoring forces, one favouring long waves and one short, leave a slowest wave in between — 23 cm/s at
  1.7 cm on clean water.

### E4 · seiche_standing_waves
- **A:** C08 (also shows N58–N63, N177–N178)
- **Confusion removed:** "a standing wave is a wave that does not move" — it is two equal waves travelling in
  opposite directions; the walls of a basin only allow the pairs whose nodes fit, which fixes the lake's periods.
- **Why interactive:** the reader toggles the two travelling components on and off and sees their sum pin its nodes;
  then changes L, H and the mode number and watches the sloshing period update, compared with real lakes.
- **Stage:** (1) the basin: surface, streamlines of (7.62) and velocity arrows, the two travelling waves as faint
  ghosts; (2) η(x) and u(x) mode shapes with nodes/antinodes; (3) (hidePortrait) period vs mode number n with the
  shallow-water line 2L/((n + 1)√(gH)).
- **Controls:** basin length L (log 1 m … 500 km) · depth H · mode n (0…4) · component toggles (right-going,
  left-going, sum) · transport.
- **Presets:** bathtub (1.5 m, 0.2 m) · swimming pool · a 50 km lake of depth 100 m · a long shallow lake (Lake Erie-like
  numbers, our inputs) · a harbour.
- **Equations shown (live):** (7.61); $\eta=2a\cos kx\cos\omega t$; (7.62) ψ; (7.63) u; (7.64) $\lambda=2L/(n+1)$; (7.65).
- **Mirrors:** `ch07.standing_wave_fields`, `seiche_modes`, `basin_modes`.
- **Derivations:** D17 (sum-to-product: steps switch the component toggles), D18 (walls: steps set n).
- **Depth features:** Explain tab (k from L and n → kH → ω → period in s/min/h → shallow approximation error →
  reading the current setting), synced Code tab, + linked views (3), transport, presets, inspector (click the basin:
  u and η at that point and time), notes ("Right now" with a highlighted table of basins).
- **Follows reference:** `angular_frequency_explorer_1.html` (system animation + a table of real cases with the current
  row highlighted).
- **Aha:** two opposite waves add into a pattern that stands still; a basin keeps only the patterns whose nodes fit,
  so a lake rings at its own periods.

### E5 · group_velocity_packets
- **A:** C09 (also shows C06's E and F, N65 (7.66), N66, N67 (7.68), N68–N70, N72 stone in a pond, N73; links C16)
- **Confusion removed:** "the wave moves at one speed" — crests move at c, the group and the energy at c_g; in deep
  water crests are born at the back of a group and die at the front.
- **Why interactive:** the reader follows one crest with a marker while the envelope moves at half its speed, switches
  the dispersion (deep, shallow, capillary) and sees c_g < c, = c, > c; the chord and tangent of ω(k) move with the dot.
  Motion is the only way to see crests passing through a group.
- **Stage:** (1) η(x, t) with the envelope (purple) and a tracked crest (orange), c and c_g arrows; (2) ω(k) with the
  chord (slope c) and tangent (slope c_g) at k₀; (3) (hidePortrait) energy density ½ρga(x)² and flux F = E c_g.
- **Controls:** mode chips: two waves (beats) / packet / stone in a pond · k₀ · bandwidth Δk (or packet width) ·
  dispersion chips: deep / finite H / shallow / capillary · transport.
- **Presets:** deep-water swell set · shallow (no dispersion: shape kept) · pure capillary (c_g = 3c/2) · wide packet
  vs narrow packet (spreading) · stone in a pond (calm centre at c_g,min).
- **Equations shown (live):** (7.66) beats; (7.67) $c_g=d\omega/dk$; (7.68) $\eta=a(x-c_gt)\cos(kx-\omega t)$; (7.69);
  (7.70); (7.71) $F=Ec_g$; (7.42).
- **Mirrors:** `core.W.beat_wave`, `group_velocity`, `group_velocity_numeric`, `gaussian_packet` (new, §9),
  `min_group_velocity`, `wave_energy_density`.
- **Derivations:** D19 (beats: the printed Δω x shown as a ghost envelope that does not move), D20 ★★★ (packet: the
  Taylor step splits the phase into carrier and envelope; the dropped quadratic term is the spreading toggle), D21 (c_g
  and F = E c_g).
- **Depth features:** Explain tab (c and c_g with numbers, their ratio, time for a crest to cross the group, energy
  arrival time from a storm → reading the current setting), synced Code tab, + linked views (3), transport, modes,
  presets, status ("crests overtake the group (c_g < c)" / "together" / "group overtakes crests"), end-of-run card
  (distance travelled by crest vs envelope).
- **Follows reference:** `amplitude_phase_second_order_II_3.html` (system + graph windows linked by one state, a
  numbered derivation with live numbers in the explanation).
- **Aha:** the crests are only the carrier; the envelope — and the energy — moves at dω/dk, half the crest speed in
  deep water.

### E6 · wave_rays_refraction
- **A:** C10 (also shows N47, N48, N74–N80, N174, N175, N183, N184)
- **Confusion removed:** "waves arrive parallel to the beach because the wind blows onshore" — any swell turns: its
  frequency is fixed along a ray while its wavelength shrinks with the depth, so the shoreward end slows and the crest
  swings round (Snell).
- **Why interactive:** the reader sets the incidence angle and period and watches rays bend over a sloping beach or
  wrap round an island, with ω printed constant along the selected ray while k, c and c_g change — the conservation is
  seen, not asserted.
- **Stage:** (1) plan view: depth contours, rays (purple) and crest lines (orange), selected ray highlighted; (2)
  along the selected ray: ω (flat), k, c, c_g and angle α vs distance; (3) (hidePortrait) x–t diagram of crests
  (slope c) and rays (slope c_g) for the homogeneous case (Fig. 7.17).
- **Controls:** bathymetry chips: plane beach / circular island / submarine ridge · incidence angle · period T ·
  beach slope · click a ray to select it.
- **Presets:** swell at 30° on a 1:50 beach · long period (turns sooner) · island (shadow-side crests) · ridge
  (focusing) · homogeneous (straight rays).
- **Equations shown (live):** (7.73) $k=\theta_x$, $\omega=-\theta_t$; (7.74); (7.75); (7.79) $\partial\omega/\partial t+c_g\partial\omega/\partial x=0$;
  Snell $k\sin\alpha=$ const (ours); (7.28).
- **Mirrors:** `core.W.ray_trace`, `ch07.snell_ray_plane_beach` (new, §9), `core.W.omega_gravity`, `group_velocity`.
- **Derivations:** D22 (crest conservation), D23 (ω along rays: steps highlight the flat ω curve), D24 (Snell: steps
  move the selected ray).
- **Depth features:** Explain tab (k at start and at the selected point from ω fixed → α from Snell → c, c_g →
  travel time along the ray → reading the current setting), synced Code tab, + linked views (3), modes (bathymetries),
  presets, inspector (click a ray point: the Snell arithmetic), status ("ω conserved to 1e-9 along the ray").
- **Follows reference:** `stride_padding_playground.html` (classic-setting presets and the formula with numbers
  plugged in) with `angular_frequency_explorer_1.html`'s linked views.
- **Aha:** along a ray the frequency never changes; the wavelength does, and that alone turns every crest toward the
  shore.

### E7 · hydraulic_jump
- **A:** C11 (also shows N81–N85, N93, N95, N97, R10 Fr; links ch04 `control_volume_budgets` bore scene)
- **Confusion removed:** "the water jumps up because it hits something" — the jump is set by momentum alone: a fast
  shallow stream must become a slow deep one with the same momentum flux, and the energy difference is burnt in the
  roller; the reverse (a drop) would create energy.
- **Why interactive:** the Fr₁ slider moves the downstream depth and the energy-loss bar together, and pushing Fr₁
  below 1 turns the energy bar positive (rose, "forbidden") — the second law choosing the root is felt; the moving-bore
  mode changes frame with a slider.
- **Stage:** (1) channel side view with the jump, CV box, face pressure forces and momentum-flux arrows; (2) H₂/H₁ and
  head loss vs Fr₁ with the current dot and the forbidden region shaded; (3) (hidePortrait) momentum budget bars
  (inflow, outflow, pressure forces) summing to zero, and the energy bars.
- **Controls:** upstream depth H₁ · Fr₁ (0.5 … 6) · frame chips: stationary jump / bore moving into still water ·
  nonlinear-fate chip: jump / solitary wave (Ursell number shown).
- **Presets:** weak undular jump (Fr₁ = 1.3) · spillway (Fr₁ = 5) · tidal bore (moving frame) · Fr₁ = 1 (no jump) ·
  Fr₁ = 0.7 (forbidden) · solitary wave a/H = 0.2.
- **Equations shown (live):** (4.104) $\mathrm{Fr}=u/\sqrt{gH}$; (7.80); (7.81); $E_2-E_1=-g(H_2-H_1)^3/(4H_1H_2)$; (7.87)
  Ursell $a\lambda^2/H^3$; (7.88).
- **Mirrors:** `ch07.hydraulic_jump`, `jump_momentum_residual`, `ursell_number`, `solitary_wave`, `ch04.bore_speed`.
- **Derivations:** D25 (Bélanger: steps light the CV faces), D26 (energy loss: step "sign" moves Fr₁ below 1).
- **Depth features:** Explain tab (Q, Fr₁ → r from the quadratic → H₂, u₂, Fr₂ → momentum check → head loss and
  dissipated power → reading the current setting), synced Code tab, + linked views (3), presets, term bars (momentum
  budget), status ("✅ jump: Fr₁ > 1" / "⛔ would create energy" / "〰 undular"), modes (frames).
- **Follows reference:** `fid_formula_lab.html` (term bars with a total) with ch04 `control_volume_budgets`' face-flux
  bars.
- **Aha:** momentum decides the jump height, energy decides the direction — only fast shallow flow may jump to slow
  deep flow.

### E8 · two_layer_modes
- **A:** C13, C14 (also shows N99 complex amplitudes, N105–N107, N124–N131, N128 g′ convention; links ch05
  `vortex_sheet_rollup`)
- **Confusion removed:** "internal waves are just small surface waves inside the water" — they are slow (√(g′H)),
  large, and nearly invisible at the surface; a layer over deep fluid has two modes, the surface one (barotropic) and the
  internal one (baroclinic) with surface and interface in antiphase.
- **Why interactive:** the density-contrast slider and kH move both dispersion curves and the mode shapes at once; the
  reader sees the interface amplitude dwarf the surface one as Δρ → 0 and recovers the two-deep-fluids case as kH → ∞.
- **Stage:** (1) side view: free surface and interface animated for the chosen mode, u arrows in both layers (reversing
  across the interface in the baroclinic mode), the vortex-sheet jump marked; (2) ω(k) of both modes with the √(gk) and
  ε√(gk) limits and the current dot; (3) (hidePortrait) amplitude ratio η/ζ and speed ratio vs Δρ/ρ on log axes.
- **Controls:** mode chips: barotropic / baroclinic / two deep fluids · Δρ/ρ₂ (log 1e-4 … 0.5) · upper-layer
  thickness H · wavelength · transport.
- **Presets:** ocean thermocline (H = 50 m, Δρ/ρ = 0.002) · fjord (fresh over salt) · oil over water (Δρ/ρ = 0.2) ·
  kH → ∞ (recover (7.95)) · long waves (g′H).
- **Equations shown (live):** (7.95) $\omega=\varepsilon\sqrt{gk}$; (7.96); (7.110); (7.111); (7.112) $b=ae^{-kH}$; (7.113);
  (7.114); (7.116)–(7.117) $c=\sqrt{g'H}$, $g'=g(\rho_2-\rho_1)/\rho_2$.
- **Mirrors:** `core.W.interface_omega`, `two_layer_free_surface_omega`, `two_layer_long_wave_speed`,
  `reduced_gravity_book`, `ch07.two_layer_modes`, `interface_fields`.
- **Derivations:** D28 (interfacial ω), D29 ★★★ (constants: steps set the mode), D30 ★★★ (the factorised relation:
  steps light each factor's curve), D31 (modes and g′).
- **Depth features:** Explain tab (ε, g′ with both ρ conventions → ω of each mode → periods → η/ζ → surface signal of a
  10 m internal wave → reading the current setting), synced Code tab, + linked views (3), transport, modes, presets,
  status ("🌊 barotropic: surface and interface in phase" / "↕ baroclinic: antiphase, surface signal ×0.002").
- **Follows reference:** `amplitude_phase_second_order_II_3.html` (one state driving system + response windows, a
  numbered derivation with live numbers).
- **Aha:** a small density step makes a second, slow wave whose interface swings metres while the surface barely moves;
  its speed is √(g′H) with g′ a tiny fraction of g.

### E9 · internal_wave_beams
- **A:** C15, C16 (also shows N145 (ω = N limit), N146 layered flow, N143–N144, N148–N150, N158–N166 energy flux; links
  ch01 `parcel_stability`)
- **Confusion removed:** "waves move energy the way their crests move" — for internal waves the frequency sets only the
  *angle*; crests move across a beam while energy moves along it, at right angles, and phase going up means energy going
  down.
- **Why interactive:** the ω/N slider tilts the St Andrew's cross beams while the crests visibly slide across them; the
  reader drags the K arrow in a wavenumber-plane view and sees ω stay fixed when only ∣K∣ changes, and c_g stay
  perpendicular to c.
- **Stage:** (1) the stratified tank: ρ′ field of four beams from an oscillating source with moving phase lines and
  particle motion along the crests; (2) wavenumber plane: K arrow (draggable), c (orange) and c_g (purple) with the
  right angle marked, the circle of constant ∣K∣ and lines of constant ω; (3) (hidePortrait) ω/N vs θ (cos curve) with
  the current dot.
- **Controls:** ω/N (0.05 … 0.99) · N · drag K (or ∣K∣ and θ sliders) · view chips: beams / single plane wave /
  packet · transport.
- **Presets:** ω = 0.71N (45° cross) · near N (beams near vertical, columns bob) · low frequency (flat beams, layered
  flow) · the book's Fig. 7.29 geometry (k < 0) · ocean thermocline N = 0.01 rad/s.
- **Equations shown (live):** (7.134) $\partial_t^2\nabla^2w+N^2\nabla_H^2w=0$; (7.137); (7.139) $\omega=N\cos\theta$; (7.141)
  $\mathbf K\cdot\mathbf u=0$; (7.144)–(7.145); (7.146) $\mathbf c_g\cdot\mathbf c=0$; (7.159) $\mathbf F=\mathbf c_gE$.
- **Mirrors:** `core.W.internal_wave_omega`, `beam_angle`, `internal_wave_velocities`, `ch07.st_andrews_cross`,
  `internal_wave_energy`.
- **Derivations:** D33 ★★★ (w-equation: steps highlight which operator is applied), D34 (ω = N cos θ: steps rotate K),
  D35 (K·u = 0: steps show particle arrows along crests), D36 (c ⟂ c_g: steps draw the two components).
- **Depth features:** Explain tab (θ from ω/N — to the horizontal for K, to the vertical for the beams — → c and c_g
  vectors with numbers → their dot product → periods → energy flux direction → reading the current setting), synced
  Code tab, + linked views (3), transport, modes (beams / plane wave / packet), presets, inspector (click the K plane:
  ω, c, c_g arithmetic), status ("beam 45.0° from the vertical; phase ↑ ⇒ energy ↓").
- **Follows reference:** `angular_frequency_explorer_1.html` (rotating-pointer + system + response curve on one
  clock) — the K arrow is its rotating pointer.
- **Aha:** an internal wave's frequency is its direction; energy runs along the beams at right angles to the crests,
  which slide across them.

### B1 · linearised_free_surface (backup)
- **A:** C02 (also shows N13–N17, R04–R06; D02–D04 in its Derivation tab if built)
- **Confusion removed:** "how can a condition on an unknown surface be applied at z = 0?" — by Taylor transfer, whose
  error grows like (ka)²; the reader drags ka and sees the exact and linearised conditions' residuals separate on a
  log–log plot, and the leaf on the surface follow ∂η/∂t.
- **Stage:** surface with the z = 0 line and a floating leaf · residual bars (kinematic exact vs linear, dynamic exact
  vs linear) · log–log residual vs ka (slope 2).
- **Mirrors:** `ch07.free_surface_residuals`, `core.interfaces.kinematic_bc_residual`. **Follows:** ch04 backup
  `kinematic_free_surface` storyboard (analysis/ch04_curation.md) and `overfitting_curves.html`.
- **Aha:** linearising and moving the condition to z = 0 costs a relative error of order ka — invisible for swell,
  visible for steep waves.

## 6. Python animations and interactive figures (A ID → what, why, player/figure kind)
Animations (`fluidpy.core.anim`; ≤ 120 frames, dpi ≤ 80; FAST halves frames):

| Kind | A ID | What moves | Why | Player |
|---|---|---|---|---|
| A1 | **C05** | three panels (kH = 3, 1, 0.3): surface and tracer particles on their orbits, orbit ghosts | the circle → ellipse → flat-ellipse change is a motion | video |
| A2 | **C09** | deep-water Gaussian packet: carrier crests (orange marker on one crest) moving through the envelope (purple) at half their speed; second part: stone-in-a-pond train (capillary–gravity FFT) | crests born at the rear and dying at the front can only be seen moving | video |
| A3 | **C08** | standing wave in a basin with streamlines of (7.62) and the two travelling components as ghosts | nodes fixed while the sum oscillates | video |
| A4 | **C12** | exact path lines of a steep (ka = 0.25) deep-water wave: open orbits and a dyed vertical line leaning forward | drift accumulates period by period | video |
| A5 | **C11** | (a) simple-wave steepening of a hump until the breaking time; (b) KdV run (FAST grid) of a hump splitting into solitons, tallest first | the two nonlinear fates side by side; stepping matters | frames |
| A6 | **C16** | internal-wave packet from `linear_evolve_2d`: phase lines moving up-left while the packet moves down-left (Fig. 7.32) | c ⟂ c_g is only believable in motion | video |

Interactive figures (`slider_figure` / `animate_figure`; ≤ 40 steps × ≤ 4 traces; survive publishing):

| Kind | A ID | Slider controls | What changes | Why |
|---|---|---|---|---|
| IF1 | **C04** | depth H (0.5 m … 5 km, log) | c(λ) with deep and shallow asymptotes and the regime bands | where depth starts to matter moves with H |
| IF2 | **C07** | surface tension σ (0 … 0.5 N/m) | c(λ) on log–log axes with the minimum marker | the minimum slides along the curve |
| IF3 | **C04** | wavelength λ at fixed H | pressure amplitude p′/ρga vs depth (cosh ratio) with hydrostatic and e^{kz} ghosts | when a bottom sensor sees the wave |
| IF4 | **C09** | time t (`animate_figure`) | Gaussian packet and its envelope; markers moving at c and c_g | page-surviving version of A2 |
| IF5 | **C10** | incidence angle α₀ | rays and crest lines over a plane beach | refraction without a kernel |
| IF6 | **C14** | Δρ/ρ₂ (log) | ω(k) of barotropic and baroclinic modes with the (7.95) limit | how weak stratification slows the internal mode |
| IF7 | **C08** | mode number n (0…5) | η(x) and u(x) of seiche mode n with nodes | mode shapes of Fig. 7.12 |
| IF8 | **C15** | ω/N | beam directions and c, c_g arrows (St Andrew's cross skeleton) | angle set by frequency |

Live widget (kernel only, paired with IF1 and E1): `live(dispersion, lam=…, H=…)` printing ω, T, c, c_g, regime (C04).

## 7. From-scratch moments
Each is a transparent hand-written version next to the tested function, followed by `assert np.allclose(mine, lib)`.

| Section | A ID | Hand-written version | Tested function it must match |
|---|---|---|---|
| §7.1 | **C01** | track crest n = 0 as the argmax of η on a fine x grid at several t; slope of x_crest(t) | ω/k from `core.W.wave_parameters` (rtol 1e-3 grid-limited) |
| §7.2 | **C03** | Newton iteration k ← k − (ω² − gk tanh kH)/(d/dk) for T = 10 s, H = 30 m | `core.W.wavenumber_from_omega` (rtol 1e-10) |
| §7.2 | **C05** | RK4 by hand (P95) for the linearised path lines (7.34) over one period | `ch07.orbit_linear` (rtol 1e-6) |
| §7.2 | **C06** | midpoint double sum of ½ρ(u² + w²) over one wavelength and the depth | `ch07.wave_energy(method='closed')` Ek (rtol 1e-4) |
| §7.3 | **C07** | scan c(λ) on a log grid, refine around the minimum | `core.W.capillary_minimum` (rtol 1e-6 after refinement) |
| §7.4 | **C08** | loop over n computing k = (n+1)π/L, ω, T | `ch07.seiche_modes` (rtol 1e-12) |
| §7.5 | **C09** | central difference (ω(k + h) − ω(k − h))/2h; envelope-peak speed from `scipy.signal.hilbert` on the FFT packet | `core.W.group_velocity` (rtol 1e-8; envelope 1e-2) |
| §7.5 | **C10** | fixed-step RK4 on dx/dt = ∂ω/∂k, dk/dt = −∂ω/∂x for a plane beach | `core.W.ray_trace` (rtol 1e-6) and Snell k sin α = const |
| §7.6 | **C11** | `np.roots([1, 1, -2*Fr1**2])`, keep the positive root | `ch07.hydraulic_jump` H₂/H₁ (rtol 1e-12) |
| §7.6 | **C12** | average dx/dt over 20 periods of the exact path line minus the start | `ch07.stokes_drift` (rtol 1e-2, O((ka)²) check at ka = 0.05) |
| §7.7 | **C13** | `np.linalg.solve` of the two kinematic conditions for A, B, then the pressure condition's ω | `core.W.interface_omega` (rtol 1e-12) |
| §7.7 | **C14** | (7.110) as a quadratic in ω², `np.roots` | `core.W.two_layer_free_surface_omega` (rtol 1e-10) |
| §7.8 | **C15** | finite-difference residual of (7.134) for a plane wave with ω from the formula (≈ 0) and with a wrong ω (≠ 0) | `core.W.internal_wave_omega` |
| §7.8 | **C16** | central-difference gradient of ω(k, m) = N∣k∣/K and the dot product with c | `core.W.internal_wave_velocities` (rtol 1e-7; dot 1e-12) |

## 8. Notes for the implementer
- **Budget (< 5 min on Colab CPU):** FFT evolutions at N = 4096 (FAST 1024), ≤ 200 frames (FAST 60); `kdv_solve` at
  N = 512 (FAST 256) and ≤ 10⁴ steps, cached to `outputs/ch07/kdv_run.npz`; Stokes-drift path lines 20 periods (FAST 10);
  `linear_evolve_2d` on 256² (FAST 128²); `dblquad` energy checks only for three kH; sympy cells (`surface_wave_sympy`,
  `two_layer_sympy`, `boussinesq_linear_sympy`, `kdv_residual_sympy`) each < 3 s.
- **Scalar-callable for parity:** every function an explainer mirrors must accept Python floats and return floats or
  small dicts of floats (`shot.py` evaluates `py:` rows); keep `H=np.inf` working and `ρ₁ > ρ₂` returning NaN with a
  warning (Rayleigh–Taylor is Ch. 11).
- **Conventions to encode in names and docstrings:** `zeta_particle` vs `zeta_interface`; `theta_phase` vs `theta_K`
  (from the horizontal) and `beam_angle` from the vertical (return both in a small dict for E9's badge); `eps2_density`;
  `reduced_gravity_book(ref="lower")` with the ch04 ρ₁ form available as `ref="upper"`; internal-wave functions use
  ∣k∣ and ∇_K ω (R9) — test the Fig. 7.29 k < 0 configuration; `G_BOOK` default in ch07 functions.
- **Typos to test as discriminating wrong variants:** T1 beats with Δω x (envelope does not move); T3 (7.105) with
  e^{i(kz−ωt)} (sympy residual ≠ 0); the printed (7.138)/(7.145) for k < 0 (c_g points the wrong way).
- **Parity with earlier chapters:** `ch07.wave_fields` = `ch04.linear_wave_fields` (1e-12); shallow c = √(gH) limit of
  `ch04.bore_speed`; ω = N limit vs `core.stratification.parcel_displacement`; interface sheet strength vs
  `ch05.vortex_sheet_strength`.
- **Promotion:** everything in `core/waves.py` is for Ch. 11/13/15; `ray_trace`, `linear_evolve`, `group_velocity_vector`
  and `internal_wave_*` are the ones Ch. 13 will call first. JS: E1, E2, E5, E9 all need a "surface + tracer
  particles" stage and E5/E8 complex amplitudes — the knowledge-keeper should watch for a `Viz.waves` helper (surface
  path, particle loop, envelope) after the build; `Viz.cx` is not required (all explainer physics is written in real form).

## 9. What the implementer must add beyond analysis §4
Functions the A items' figures and the explainers need that analysis §4 did not plan (all scalar-callable, documented
with book § and equations, validation label):

| Function (module) | For | What |
|---|---|---|
| `core.W.gaussian_packet(x, t, a, k0, sigma_x, omega_fn, order=2)` | C09, E5, IF4 | closed-form packet with ω(k) Taylor-expanded to second order (exact chirped Gaussian): envelope and carrier separately; order=1 gives D20's a(x − c_g t) — JS parity without an FFT |
| `ch07.dispersion_state(lam, H, g=G_BOOK, rho=1000.)` | E1 | dict(k, kH, omega, T, c, cg, regime, deep_error, shallow_error, p_bottom_fraction) |
| `ch07.orbit_state(z0, a, k, H, g, periods=1, model="linear")` | E2 | dict(A, B, focal_half, drift_per_period, eulerian_mean) for one depth (drift from (7.86) and from the exact path line) |
| `ch07.dyed_line(z0s, t, a, k, H, g)` | C12, A4, E2 | positions at time t of particles starting on a vertical line (exact path lines, vectorised over z0s) |
| `ch07.capillary_state(lam, sigma, rho, H, g)` | E3 | dict(c, cg, gravity_term, tension_term, regime, c_min, lam_m) |
| `ch07.seiche_state(L, H, n, g)` | E4 | dict(k, lam, omega, T, T_shallow, shallow_error) (wraps `seiche_modes`) |
| `ch07.snell_ray_plane_beach(x, alpha0, x0, T, slope, g)` | C10, E6, IF5 | closed-form ray on a plane beach from k(x) sin α = const with full dispersion (brentq for k(x)); parity target for the JS ray |
| `ch07.jump_state(H1, Fr1, g, frame="stationary")` | E7 | dict(H2, u1, u2, Fr2, dE, head_loss, momentum_terms (4), bore_speed in the moving frame, allowed) |
| `ch07.two_layer_state(k, H, rho1, rho2, g, mode)` | E8 | dict(omega_bt, omega_bc, eta_over_zeta, g_prime_lower, g_prime_upper, c_long, T) |
| `ch07.internal_wave_state(omega_over_N, N, K, direction=+1)` | E9, IF8 | dict(theta_K_from_horizontal, beam_from_vertical, k, m, c_vec, cg_vec, dot, periods) |
| `ch07.free_surface_residual_scan(ka_values, kH, g)` | C02 figure, backup B1 | relative residuals of exact vs linearised kinematic and dynamic conditions vs ka (log–log slope 2) |
