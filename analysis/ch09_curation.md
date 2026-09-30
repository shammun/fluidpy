# Chapter 9 — Boundary Layers and Related Topics: curation
(from `analysis/ch09.md` — 148 inventory rows (#1–#148), 85 equation numbers (9.1)–(9.85), 27 derivations in §2b, 8 planned explainer seeds;
`book.yaml` ch09 and `policy.tier_a_full_treatment: [12, 18]`, `coverage: exhaustive`. Curated 2026-09-29, concept-curator. Read:
`knowledge/CUMULATIVE.md`, `concept_map.md`, `primers.md` (P01–P199; **ch09 primers start at P200**), `notation.md`, `viz_patterns.md`,
`analysis/ch08_curation.md` (format), the `interactive-viz` skill §4–§5.)

Counts: A 14 · B 118 · C 16 · RECAP 9 (inside B/C) · SKIP 1 (inside C) · derivations written out 22 (★ 5 ★★ 12 ★★★ 5) · demoted to statements 6

Tier words (parsed by `tools/nbkit.py`): **CORE 14 · NOTE 124 · RECAP 9 · SKIP 1** = 148 rows; **DERIVATION 22**.
Depth: A 14 (all CORE) · B 118 (115 NOTE + 3 RECAP: R01, R04, R05) · C 16 (9 NOTE + 6 RECAP: R02, R03, R06, R07, R08, R09 + 1 SKIP S01).
**Reconciliation with analysis §2:** the 148 rows `#1…#148` each appear exactly once in §2 below (checked by script, last column). A = 14 of 148 rows (9.5 %).

**Decisions that shape this chapter.**
1. **Fourteen A items from the analyst's "~14 natural blocks".** The analyst's twelve-plus groups collapse as follows.
   (a) Scaling, the boundary-layer equations, matching to the outer flow and the boundary conditions form **one** A block **C01**
   (rows #5–#20): the scaling argument *is* the derivation of (9.9)–(9.10) (D01) and (9.11) is its last step (D02); the boundary
   conditions and the parabolic remark are B items inside it.
   (b) The three thicknesses are **one** A block **C02**: δ* is the A row, θ (9.17) and δ₉₉ ride with it (they share every
   integral, every figure and the explainer).
   (c) Blasius is two A blocks: **C03** the reduction to (9.27) (a ★★ derivation, D05) and **C04** the solution and its numbers
   (Töpfer scaling, f″(0), the wall-shear and drag results); Falkner–Skan, the momentum integral and Thwaites are one each
   (**C05**, **C06**, **C07**).
   (d) Separation gets **C08** (wall curvature → inflection → τ₀ = 0) and **C09** (form drag and streamlining, with the transition
   and flat-plate drag curve as B items); the cylinder is **C10** (regimes and the Kármán street) and **C11** (the drag crisis,
   with the sphere and the sports balls as B items — §9.9 has no A item of its own).
   (e) The two jets are **C12** and **C13**; the free jet carries one ★★★ (the similarity ODE), the wall jet three (the invariant,
   the reduction — which the book prints with a wrong coefficient — and the integration); **C14** is the teacup, kept as an A
   item because it is the bridge to Ekman layers and spin-up in Ch. 13.
   The analyst's other candidates were not promoted: Falkner–Skan separation member (B inside C05, ★ a-D09 demoted), Example 9.1 and 9.2
   (worked numbers inside C07), the flat-plate drag curve and transition (B inside C09), sphere and ball dynamics (B inside C11).
2. **SEEN rows are RECAP** (nine rows R01–R09, inventory order): (9.1), (9.7-companion) continuity, (6.2), no slip, no through-flow,
   and the four repeated boundary conditions (9.20), (9.21), (9.77), (9.78). All other "SEEN idea, NEW instance" rows are B items
   (they add a new equation).
3. **SKIP is minimal** — only S01 (exercises 9.1–9.28, literature lists). The analyst's prose SKIP rows (#21, #23, #31, #60, #89, #91,
   #96, #100, #104) hold new ideas and are C items with a pointer, so nothing that is new disappears.
4. **Derivations (D01–D22).** Every derivation of an A item is written out (21 from the analyst's 27 plus one of ours: the
   pressure-drag integral of the separated model, D13); six analysis derivations are demoted to statements (§4c). Five are ★★★
   and each has an explainer that shows it (D14 → `karman_street_stability`, D16 → `free_jet_similarity`, D19–D21 → `wall_jet_invariant`).
5. **Printed slips (analysis §9) become named callouts** "book prints X, correct is Y", each with a planted wrong variant a test must
   fail: R1 (9.7) squares, R2 (9.30) 4.93 vs 4.910, R3 wall-jet coefficient 1 vs 4, R4 wall-jet separation of variables, R6 (9.76) 5.6152 vs
   7.3319, R7 (9.56) τ units, R10 the Magnus inequality, R11 reverse-flow claim for n < −0.0904, R14 "Chapter 13" should be 12.
   R5 (one side of the plate), R8 (Ψ dimensions), R9 (transition Reynolds numbers), R12 (Blasius the name), R13 (sphere eddy range),
   R15 (λ, L, l symbols) are shown as `⚠️ Common confusion` callouts; R16 (exercise answers are computable) is a rule for the
   writers: print at most 4 significant figures and never the private book values.
6. **Book values stay private.** Table 9.1, regime thresholds, angles, sports-ball figures, jet constants and exercise answers are used
   only as V6 evidence from `tests/book_values_ch09.json`; the public closure for Thwaites is our Falkner–Skan-exact table (D10).


## 1. Teaching order (A IDs grouped by book section, B/C IDs under each; one sentence each: "once you see X, Y follows")
Order = book order, with two prerequisite moves: C02 (thicknesses) comes right after C01 so that C04 can quote its numbers, and C08 (separation) is taught after Thwaites so that τ₀ = 0 has a computed meaning. RECAP rows are reminders inside the block named as their A parent.

**§9.1 Introduction**
- **C01** Boundary-layer scaling and the boundary-layer equations (9.4)-(9.10) with the outer matching (9.11): once you see that advection and viscosity balance at δ/L ~ Re^(-1/2), the stretched variables (9.6) shrink the Navier–Stokes equations to (9.9)–(9.10) and Bernoulli hands the pressure to the outer flow; every later section is a solution of these equations.
  - B: N01, N02, N03, N04, N05, N06, N07, N08, N09, N10, N11, N12, N13, N14, N16
  - C: N15, N17   ·   RECAP: R01, R02, R03, R04, R05
**§9.2 Boundary-Layer Thickness Definitions**
- **C02** Displacement thickness delta* (9.16), with delta_99 and theta (9.17) and the shape factor as its companions: once you see δ* as the height by which the outer streamlines are pushed up and θ as the momentum the wall has stolen, the three thicknesses become a ruler you can lay on any profile; everything after this measures a layer with them.
  - B: N18, N19, N20, N21
  - C: –   ·   RECAP: –
**§9.3 Boundary Layer on a Flat Plate: Blasius Solution**
- **C03** Blasius equation f''' + f f''/2 = 0 (9.27) by similarity reduction of (9.18): once you see that a flat plate has no length scale of its own, ψ = Uδ(x)f(η) turns the PDE (9.18) into the single ODE (9.27), and the same trick will work for wedges and jets.
  - B: N22, N23, N25, N26, N27, N28, N29, N30, N31, N32
  - C: N24   ·   RECAP: R06, R07
- **C04** Blasius solution and its numbers: f''(0) = 0.332, delta*, theta, tau0 (9.31), Cf (9.32), C_D (9.33): once you can integrate (9.27) with the Töpfer rescaling and read off f″(0) = 0.332, every Blasius number (4.91, 1.72, 0.664, 1.328) is one integral of one curve; these numbers are the yardstick for approximate methods and Ch. 10 codes.
  - B: N33, N34, N35, N36, N37, N38, N39, N40, N41, N42
  - C: –   ·   RECAP: –
**§9.4 Falkner-Skan Similarity Solutions**
- **C05** Falkner-Skan equation (9.36): wedge flows U_e = a x^n, Blasius at n = 0, stagnation at n = 1: once you see that a power-law outer flow U_e = a xⁿ keeps the similarity form, the single equation (9.36) contains Blasius (n = 0), stagnation flow (n = 1) and the last attached profile (n ≈ −0.09, zero wall shear).
  - B: N43, N44, N45, N46, N47
  - C: N48   ·   RECAP: –
**§9.5 Von Karman Momentum Integral Equation**
- **C06** Von Karman momentum integral equation (9.43): once you integrate (9.9) across the layer, the PDE collapses to one ODE (9.43) linking θ, δ* and τ₀ for any pressure gradient, at the price of a closure.
  - B: N49, N50, N51, N52, N53, N54, N55
  - C: –   ·   RECAP: –
**§9.6 Thwaites' Method**
- **C07** Thwaites' method: closed form (9.50) for theta and the separation test: once you accept that l and H depend on the single parameter λ, (9.43) integrates in closed form (9.50) and one integral of U_e⁵ gives θ, τ₀ and the separation point of any attached laminar layer.
  - B: N56, N57, N58, N59, N60, N61, N62, N63, N64, N65, N66, N67
  - C: –   ·   RECAP: –
**§9.7 Transition, Pressure Gradients, and Boundary-Layer Separation**
- **C08** Boundary-layer separation: wall curvature (9.51)-(9.52), inflection point, tau0 = 0: once you evaluate (9.9) at the wall (μ u_yy = dp/dx) you see why an adverse gradient creates an inflection point, thickens the layer and finally sets τ₀ = 0, which is separation.
  - B: N70, N71, N72, N73
  - C: N74   ·   RECAP: –
- **C09** Streamlining and form drag: separated wake pressure, laminar vs turbulent separation angles: once you see that separation replaces the ideal rear pressure recovery by a flat low-pressure wake, form drag follows from where the layer lets go; streamlining is delaying separation.
  - B: N68, N69
  - C: N75   ·   RECAP: –
**§9.8 Flow Past a Circular Cylinder**
- **C10** Cylinder wake regimes and the Karman vortex street with its stable spacing 0.28: once you accept that a wake with no rear pressure recovery is unstable, the shedding follows, and a point-vortex model of the street explains why only the ratio b/a = 0.28 survives.
  - B: N76, N77, N78, N80
  - C: N79   ·   RECAP: –
- **C11** Drag crisis: critical Re, separation moving aft (82 deg to 125 deg), C_D collapse; the base of sports-ball dynamics: once you know that a turbulent layer separates later (82° → 125° on the cylinder), the sudden fall of C_D near Re ≈ 3×10⁵ and every trick of the sports ball (seam, spin, roughness) follows.
  - B: N81, N83, N84, N85
  - C: N82, N86   ·   RECAP: –
**§9.10 Two-Dimensional Jets**
- **C12** Free 2-D laminar jet: sech^2 profile (9.71), widths x^(2/3), entrainment (9.73)-(9.75): once you see that dp/dx = 0 makes the momentum flux J of a jet constant, similarity forces u₀ ∝ x^(-1/3), δ ∝ x^(2/3), the ODE 3f‴ + f f″ + f′² = 0 integrates twice to tanh, and entrainment ṁ ∝ x^(1/3) drops out.
  - B: N87, N88, N89, N90, N91, N92, N93, N94, N95, N96, N97, N98, N99, N100, N101, N102, N103, N104, N105, N106, N107, N108, N109, N110, N111, N112, N113
  - C: –   ·   RECAP: –
- **C13** Wall jet: conserved exterior-momentum flux (9.80), ODE 4f''' + f f'' + 2f'^2 = 0 and the implicit solution (9.83): once you see that the wall removes ordinary momentum flux, a different hidden invariant (9.80) survives, gives u₀ ∝ x^(-1/2), δ ∝ x^(3/4) and the ODE 4f‴ + f f″ + 2f′² = 0; a first integral and partial fractions then give (9.83).
  - B: N114, N115, N116, N117, N118, N119, N120, N121, N122, N123, N124
  - C: –   ·   RECAP: R08, R09
**§9.11 Secondary Flows**
- **C14** Secondary flow: the teacup, radial force imbalance in a thin layer: once you see that a thin layer slows the swirl but not the pressure gradient, the unbalanced pressure drives fluid inward along the bottom, and that is the mechanism of Ekman layers and spin-up in Ch. 13.
  - B: –
  - C: –   ·   RECAP: –

**§9.9 Flow Past a Sphere and the Dynamics of Sports Balls** — no A item of its own: N83–N86 sit under C11 (the sphere is the same drag crisis at Re_cr ≈ 5×10⁵; the ball tricks are its consequences).

**End matter** — S01 (exercises, literature, supplemental reading).

## 2. Chapter map (depth) — every inventory row exactly once
Rows are in analysis §2 order (inventory row number in the last column). IDs: `C` = A/CORE (teaching order), `N` = B or C NOTE (inventory order), `R` = RECAP, `S` = SKIP.

| ID | Item | § | Depth | Tier | A parent | Reason (A) / treatment (B) / pointer (C, SKIP) / source chapter (RECAP) | Row |
|---|---|---|---|---|---|---|---|
| N01 | Prandtl's boundary-layer hypothesis: viscosity matters only in a thin no-slip layer, irrotational outer flow elsewhere | 9.1 | B | NOTE | C01 | picture + one paragraph; resolves d'Alembert (ch06 C06); our sketch with delta growing like sqrt(x) | #1 |
| R01 | Eq. (9.1) steady 2-D x-momentum along the surface (delta/R << 1) | 9.1 | B | RECAP | C01 | ch04 C08 (4.10) and ch08 C05 (8.1): steady 2-D Navier-Stokes; reused as the starting line of D01 | #2 |
| N02 | Eq. (9.2) size of the advective term u u_x ~ U^2/L | 9.1 | B | NOTE | C01 | stated with one number (U = 1 m/s, L = 1 m); inside D01 | #3 |
| N03 | Eq. (9.3) size of the viscous term nu u_yy ~ nu U/delta^2 | 9.1 | B | NOTE | C01 | stated with the same number | #4 |
| N04 | Eq. (9.4) thickness estimate delta/L ~ Re^(-1/2) (advection = viscosity; = sqrt(nu t) with t = L/U) | 9.1 | B | NOTE | C01 | stated + the ch08 bridge sqrt(nu t); worked number for air; the first result of D01 | #5 |
| N05 | Eq. (9.5) derivative scalings d/dx ~ 1/L, d/dy ~ 1/delta | 9.1 | B | NOTE | C01 | one paragraph, reuses ch08 P188 (two length scales) | #6 |
| N06 | Scale of v from continuity: v ~ U Re^(-1/2) << u | 9.1 | B | NOTE | C01 | stated with a number; used to shrink v* in (9.6) | #7 |
| N07 | Eq. (9.6) stretched variables y* = y Re^(1/2)/L, v* = v Re^(1/2)/U (epsilon = Re^(-1/2) of (8.14)) | 9.1 | B | NOTE | C01 | stated; same move as ch08 (8.14); substitution rules shown in D01 | #8 |
| R02 | Re-displayed (8.15): scaled continuity du*/dx* + dv*/dy* = 0, no coefficient | 9.1 | C | RECAP | C01 | ch08 C05, Eq. (8.15): scaled continuity; reminder only | #9 |
| N08 | Eq. (9.7) scaled x-momentum: only 1/Re on d2u/dx2, coefficient 1 on d2u/dy2 (printed denominators lack squares, slip R1) | 9.1 | B | NOTE | C01 | stated + 'book prints X, correct is Y' callout for R1; coefficients from sympy in D01 | #10 |
| N09 | Eq. (9.8) scaled y-momentum: pressure term alone at O(1) | 9.1 | B | NOTE | C01 | stated; the reason dp/dy = 0 follows; sympy coefficient table | #11 |
| R03 | Re-displayed continuity (6.2) | 9.1 | C | RECAP | C01 | ch04 C02 / ch06: 2-D continuity (6.2); reminder only | #12 |
| C01 | Boundary-layer scaling and the boundary-layer equations (9.4)-(9.10) with the outer matching (9.11) | 9.1 | A | CORE | – | load-bearing hub: every later section (Blasius, Falkner-Skan, Karman, Thwaites, jets) solves these equations; the Re^(-1/2) thickness scaling is the chapter's central idea | #13 |
| N10 | Eq. (9.10) dp/dy = 0: pressure imposed across the layer by the outer flow | 9.1 | B | NOTE | C01 | stated + demo of its Re^-1 smallness on a Blasius field (BL.bl_pressure_variation); fails at separation or strong curvature | #14 |
| N11 | Eq. (9.11) matching: -(1/rho) dp/dx = U_e dU_e/dx (Bernoulli at the edge) | 9.1 | B | NOTE | C01 | stated with D02 (Bernoulli differentiated); this is the forcing of every later section | #15 |
| R04 | Eq. (9.12) no slip u(x,0) = 0 | 9.1 | B | RECAP | C01 | ch04 C14 and ch08 C01: no-slip condition; one reminder line at the solver BC | #16 |
| R05 | Eq. (9.13) no through-flow v(x,0) = 0 | 9.1 | B | RECAP | C01 | ch08 C01 Eq. (8.2): no through-flow; reminder line (suction, Exercise 9.26, is named only) | #17 |
| N12 | Eq. (9.14) matching u -> U_e(x) as y/delta -> infinity | 9.1 | B | NOTE | C01 | stated; truncation at eta_max and insensitivity shown in C04 | #18 |
| N13 | Eq. (9.15) inlet profile u(x0,y) = u_in(y) | 9.1 | B | NOTE | C01 | stated; used by the marching demo (N14) | #19 |
| N14 | Parabolic character: x is time-like, information flows downstream only (NS is elliptic) | 9.1 | B | NOTE | C01 | paragraph + marching animation (BL.march_boundary_layer vs Blasius); links to ch08 Stokes-problem parabola and Ch. 10 | #20 |
| N15 | Two-step viscous-inviscid procedure and displacement-body iteration | 9.1 | C | NOTE | C01 | named; used in C02 (delta* correction) and Ch. 10, Ch. 14 | #21 |
| N16 | Crude wall-stress estimate tau0 ~ mu U/delta, Cf ~ 2/sqrt(Re) | 9.1 | B | NOTE | C01 | stated with numbers; compared with 0.664 in C04 (factor 3); a-D24 stated | #22 |
| N17 | When the approximation fails: leading edge Re_x <~ 1, delta/R not small, separated flow | 9.1 | C | NOTE | C01 | named; leading-edge behaviour shown in C04, separation in C08 | #23 |
| N18 | delta_99 definition (u = 0.99 U_e) | 9.2 | B | NOTE | C02 | definition + BL.delta_level; used in the thicknesses explainer | #24 |
| C02 | Displacement thickness delta* (9.16), with delta_99 and theta (9.17) and the shape factor as its companions | 9.2 | A | CORE | – | load-bearing: the three thicknesses are the measuring stick of every later result (Blasius numbers, Karman integral, Thwaites, Ch. 12 and Ch. 14) | #25 |
| N19 | delta* as streamline displacement; v(y) = -integral of u_x dy; outer flow gets a v_inf | 9.2 | B | NOTE | C02 | stated + streamline figure from Blasius fields; links to (9.39) in C06 | #26 |
| N20 | Eq. (9.17) momentum thickness theta (also an angle elsewhere in the book: name clash flagged) | 9.2 | B | NOTE | C02 | stated with the profile-weighted integral and the shape factor H; derivation D04 ties it to the wall force | #27 |
| N21 | Flat-plate drag from theta: rho U^2 theta = integral of tau0 dx | 9.2 | B | NOTE | C02 | stated (D04 shows it); ch04 C04 wake defect parity | #28 |
| N22 | Eq. (9.18) flat-plate boundary-layer equation u u_x + v u_y = nu u_yy | 9.3 | B | NOTE | C03 | (9.9) with dp/dx = 0; input of D05 | #29 |
| N23 | Eq. (9.19) similarity ansatz psi = U delta(x) f(eta), eta = y/delta | 9.3 | B | NOTE | C03 | stated with the dimensional argument (no imposed length); first step of D05; specialises ch08 (8.32) | #30 |
| N24 | Memory of the inlet profile: favourable gradients forget it, adverse never (Serrin, Peletier) | 9.3 | C | NOTE | C03 | named; shown by the marching demo N14; Ch. 10 | #31 |
| R06 | Eq. (9.20) wall conditions u = v = 0 | 9.3 | C | RECAP | C03 | = (9.12) and (9.13) (R04, R05); reminder in the Blasius block | #32 |
| R07 | Eq. (9.21) edge condition u -> U | 9.3 | C | RECAP | C03 | = (9.14) (N12); reminder in the Blasius block | #33 |
| N25 | Eq. (9.22) delta -> 0 as x -> 0 (fixes the integration constant D = 0) | 9.3 | B | NOTE | C03 | stated; used inside D05 | #34 |
| N26 | Eq. (9.23) u = psi_y = U f'(eta) | 9.3 | B | NOTE | C03 | stated; D05 step | #35 |
| N27 | Eq. (9.24) v = -psi_x = U delta' (eta f' - f) | 9.3 | B | NOTE | C03 | stated; D05 step (the product/chain-rule move the book skips) | #36 |
| N28 | Substitution into (9.18): the eta f' f'' terms cancel | 9.3 | B | NOTE | C03 | D05 step; sympy cancelled_terms | #37 |
| N29 | Eq. (9.25) reduced equation with two x-only brackets | 9.3 | B | NOTE | C03 | D05 step | #38 |
| N30 | Eq. (9.26) delta(x) = sqrt(nu x/U) (same root as sqrt(nu t), t <-> x/U) | 9.3 | B | NOTE | C03 | D05 result; parity with core.laminar.similarity_variable | #39 |
| C03 | Blasius equation f''' + f f''/2 = 0 (9.27) by similarity reduction of (9.18) | 9.3 | A | CORE | – | hub: first exact BL solution; the similarity move recurs for Falkner-Skan and both jets | #40 |
| N31 | Eq. (9.28) wall conditions on f: f = f' = 0 | 9.3 | B | NOTE | C03 | stated; solver BC | #41 |
| N32 | Eq. (9.29) f' -> 1 as eta -> infinity (truncate and test insensitivity) | 9.3 | B | NOTE | C03 | stated + truncation study eta_max = 8, 12, 16 | #42 |
| C04 | Blasius solution and its numbers: f''(0) = 0.332, delta*, theta, tau0 (9.31), Cf (9.32), C_D (9.33) | 9.3 | A | CORE | – | the chapter's benchmark solution: every approximate method and Ch. 10 code is tested against it | #43 |
| N33 | Asymptotics f' - 1 ~ (1/eta) exp(-eta^2/4) (Gaussian tail) | 9.3 | B | NOTE | C04 | stated + residual check against the numerical f'; a-D07 demoted | #44 |
| N34 | Wall-normal velocity v/U = (eta f' - f)/(2 sqrt(Re_x)) -> 0.86/sqrt(Re_x) | 9.3 | B | NOTE | C04 | stated with the Fig. 9.6-style curve; streamlines lift outward | #45 |
| N35 | Eq. (9.30) delta_99 = 4.93 sqrt(nu x/U) (exact root 4.910; slip R2) | 9.3 | B | NOTE | C04 | stated + 'book prints 4.93, root of f' = 0.99 is 4.910' callout (brentq) | #46 |
| N36 | delta* = 1.72 sqrt(nu x/U), theta = 0.664 sqrt(nu x/U) = 2 f''(0) | 9.3 | B | NOTE | C04 | stated; integral identities shown in D06 | #47 |
| N37 | Eq. (9.31) wall shear tau0 = 0.332 rho U^2/sqrt(Re_x) | 9.3 | B | NOTE | C04 | stated; x^(-1/2) leading-edge singularity is integrable | #48 |
| N38 | Eq. (9.32) skin-friction coefficient Cf = 0.664/sqrt(Re_x) | 9.3 | B | NOTE | C04 | stated + one number; used in Thwaites Example 9.1 (N66) | #49 |
| N39 | Drag per unit width F_D = 0.664 rho U^2 L/sqrt(Re_L), proportional to U^(3/2) | 9.3 | B | NOTE | C04 | stated (integral of x^(-1/2)); a-D24 demoted | #50 |
| N40 | Eq. (9.33) drag coefficient C_D = 1.33/sqrt(Re_L), one side only (slip R5) | 9.3 | B | NOTE | C04 | stated + sides=1 callout | #51 |
| N41 | Blasius vs temporal boundary layer: Cf sqrt(Re_x) 0.664 vs 1.128 vs plate mean 1.328 | 9.3 | B | NOTE | C04 | stated + comparison figure with ch08 temporal_bl_wall_stress | #52 |
| N42 | Figs. 9.5-9.6 with our solver: profile collapse at several x | 9.3 | B | NOTE | C04 | our figure; the collapse is the heart of the Blasius explainer | #53 |
| N43 | Eq. (9.34) Falkner-Skan ansatz psi = sqrt(nu x U_e) f(eta), U_e = a x^n | 9.4 | B | NOTE | C05 | stated; extends (9.19); first step of D07 | #54 |
| N44 | Eq. (9.35) pressure gradient -dp/dx = n a^2 x^(2n-1) | 9.4 | B | NOTE | C05 | stated; n > 0 favourable, n < 0 adverse | #55 |
| N45 | Generic thickness delta(x) = sqrt(nu x/U_e) (grows for n < 1, constant at n = 1) | 9.4 | B | NOTE | C05 | stated + three x-values as a table | #56 |
| C05 | Falkner-Skan equation (9.36): wedge flows U_e = a x^n, Blasius at n = 0, stagnation at n = 1 | 9.4 | A | CORE | – | hub: turns one ODE into a whole family of pressure gradients; supplies the exact closure for Thwaites and the separation member | #57 |
| N46 | Fig. 9.7 profile family; sign of u_yy at the wall equals -f'''(0) = n | 9.4 | B | NOTE | C05 | our family figure + slider; f'''(0) = -n read from the ODE at eta = 0 (a-D09 demoted) | #58 |
| N47 | Separation member n = -0.0904 (f''(0) = 0); slip R11 on the reverse-flow claim | 9.4 | B | NOTE | C05 | stated + saddle-node/fold explained in words + 'book says X, second branch only' callout | #59 |
| N48 | Why approximate methods: real flows are rarely similar | 9.4 | C | NOTE | C05 | named; leads to C06, C07; full numerics in Ch. 10 | #60 |
| N49 | Eq. (9.37) BL momentum with the pressure gradient replaced by U_e U_e' and tau = mu u_y | 9.5 | B | NOTE | C06 | D08 start line | #61 |
| N50 | Eq. (9.38) add u times continuity: conservative form | 9.5 | B | NOTE | C06 | D08 step (product rule backwards) | #62 |
| N51 | Eq. (9.39) integral of continuity: integral of u_x dy = -v_inf | 9.5 | B | NOTE | C06 | D08 step | #63 |
| N52 | Eq. (9.40) integral of (9.38) with tau(inf) = 0, tau(0) = tau0 | 9.5 | B | NOTE | C06 | D08 step | #64 |
| N53 | Eq. (9.41) eliminate v_inf, pull d/dx out of the integral | 9.5 | B | NOTE | C06 | D08 step (Leibniz; separately divergent integrals combine) | #65 |
| N54 | Eq. (9.42) product rule to combine the integrals | 9.5 | B | NOTE | C06 | D08 step | #66 |
| C06 | Von Karman momentum integral equation (9.43) | 9.5 | A | CORE | – | hub: exact integral law for any laminar or turbulent layer; Thwaites, Ch. 12 and Ch. 14 all start from it | #67 |
| N55 | Closure problem: three unknowns (theta, delta*, tau0), one equation; Pohlhausen assumed profiles | 9.5 | B | NOTE | C06 | paragraph + optional Karman-Pohlhausen cubic run vs Blasius (BL.karman_pohlhausen) | #68 |
| N56 | Eq. (9.44) Holstein-Bohlen parameter lambda = (theta^2/nu) dU_e/dx | 9.6 | B | NOTE | C07 | stated + sign meaning (lambda < 0 adverse) | #69 |
| N57 | Eq. (9.45) shear correlation tau0 = mu (U_e/theta) l(lambda) | 9.6 | B | NOTE | C07 | stated; l from our Falkner-Skan closure (D10) | #70 |
| N58 | Eq. (9.46) shape factor H(lambda) = delta*/theta | 9.6 | B | NOTE | C07 | stated; H from our Falkner-Skan closure | #71 |
| N59 | Table 9.1 l(lambda), H(lambda), l(-0.09) = 0 (book table private) | 9.6 | B | NOTE | C07 | not reproduced: public closure = our FS-exact table (D10) and 'white' fit; both separation criteria reported | #72 |
| N60 | Eq. (9.47) multiply (9.43) by rho theta/(mu U_e) | 9.6 | B | NOTE | C07 | D09 step | #73 |
| N61 | l(lambda) = (2+H) lambda + (U_e/2) d(theta^2/nu)/dx (unnumbered) | 9.6 | B | NOTE | C07 | D09 step (book skips the algebra; we add it) | #74 |
| N62 | Eq. (9.48) universal function L(lambda) = 2l - 2(2+H) lambda | 9.6 | B | NOTE | C07 | D09 step | #75 |
| N63 | Fig. 9.8: L is nearly linear, L ~ 0.45 - 6.0 lambda | 9.6 | B | NOTE | C07 | our L(lambda) on the FS family with the line overlaid (D10 explains why it is a fit, not exact) | #76 |
| N64 | Eq. (9.49) linear fit gives a first-order linear ODE for theta^2/nu (integrating factor U_e^6) | 9.6 | B | NOTE | C07 | D09 step | #77 |
| C07 | Thwaites' method: closed form (9.50) for theta and the separation test | 9.6 | A | CORE | – | the practical engineering result of the chapter: one integral of U_e^5 predicts wall shear and separation for any body | #78 |
| N65 | Accuracy +-3 % favourable, +-10 % adverse; Thwaites predicts existence, not location, of separation | 9.6 | B | NOTE | C07 | measured against exact FS and the marching solver (our accuracy plot) | #79 |
| N66 | Example 9.1: Thwaites for Blasius, theta = 0.671 sqrt(nu x/U) (1 % high), Cf 1.2 % low | 9.6 | B | NOTE | C07 | worked number of C07; a-D12 demoted (arithmetic stated) | #80 |
| N67 | Example 9.2: diffuser, lambda(x/L) closed form, separation at (1+x/L)^4 = 1.8 | 9.6 | B | NOTE | C07 | worked number of C07 and the notebook's diffuser figure; a-D13 demoted | #81 |
| N68 | Transition of a laminar layer: Re_cr ~ 10^6 on a flat plate and the regime sequence along the plate | 9.7 | B | NOTE | C09 | stated + our schematic (plate regimes); BL.transition_state; instability mechanism deferred to Ch. 11 | #82 |
| N69 | Turbulent vs laminar layer and the flat-plate drag curve (laminar 1.33/sqrt(Re_L), turbulent branch, transition patch) | 9.7 | B | NOTE | C09 | our log-log figure; turbulent correlation only with a cited source (analysis §8); Ch. 12 | #83 |
| N70 | Wall relation mu u_yy(wall) = dp/dx from (9.9) at y = 0 | 9.7 | B | NOTE | C08 | D12 first step | #84 |
| N71 | Eq. (9.51) accelerating stream: u_yy(wall) < 0, no inflection | 9.7 | B | NOTE | C08 | stated; BL.profile_inflection on FS n > 0 | #85 |
| N72 | Eq. (9.52) decelerating stream: u_yy(wall) > 0, hence an inflection point (Ch. 11 hook) | 9.7 | B | NOTE | C08 | stated; BL.profile_inflection on FS n < 0; Blasius inflection at the wall | #86 |
| N73 | Adverse gradient thickens the layer (v(y) = -integral of u_x); Fig. 9.12 profile trio | 9.7 | B | NOTE | C08 | our figure with n = 1, 0, -0.05 | #87 |
| C08 | Boundary-layer separation: wall curvature (9.51)-(9.52), inflection point, tau0 = 0 | 9.7 | A | CORE | – | load-bearing: explains drag, stall and wakes; feeds Ch. 11 (inflection), Ch. 14 (stall) | #88 |
| N74 | After separation the BL equations fail (Goldstein singularity, unsteady, pressure not potential) | 9.7 | C | NOTE | C08 | named; marching stops at tau0 = 0; Ch. 10 and Ch. 14 | #89 |
| C09 | Streamlining and form drag: separated wake pressure, laminar vs turbulent separation angles | 9.7 | A | CORE | – | connects the boundary layer to a body's drag; the reasoning behind the drag crisis | #90 |
| N75 | Internal separation: diffusers, elbows, valves | 9.7 | C | NOTE | C09 | named; diffuser is Example 9.2 (N67); Ch. 12 pipes | #91 |
| N76 | Cylinder at Re < 1: symmetric creeping flow, Stokes/Oseen vorticity diffusion | 9.8 | B | NOTE | C10 | one paragraph + ch08 Oseen sphere as stand-in (Stokes paradox for cylinders is ch08 C15) | #92 |
| N77 | Cylinder 4 < Re < 40: two steady attached eddies growing with Re | 9.8 | B | NOTE | C10 | regime table (book's rounded thresholds flagged) + BB.cylinder_flow_regime | #93 |
| C10 | Cylinder wake regimes and the Karman vortex street with its stable spacing 0.28 | 9.8 | A | CORE | – | load-bearing: wake instability, Strouhal number and vortex shedding recur in Ch. 11 and Ch. 13 | #94 |
| N78 | Strouhal number St = Omega d/U ~ 0.2 and shedding frequency f = St U/d | 9.8 | B | NOTE | C10 | ch04 (4.102) recalled; wire example with f vs Omega trap; BB.shedding_frequency | #95 |
| N79 | Vortex-induced vibration, strakes, atmospheric Karman streets behind mountains | 9.8 | C | NOTE | C10 | named; stratification and 2-D-ness return in Ch. 13 | #96 |
| N80 | Wake regimes: laminar street Re < 200, irregular above, turbulent wake beyond several thousand | 9.8 | B | NOTE | C10 | regime table (flagged qualitative) shown in the cylinder explainer | #97 |
| N81 | Subcritical flow (Re < 3e5): separation near 82 deg, C_p plateau in the wake, C_D near 1 | 9.8 | B | NOTE | C11 | our separated-C_p figure (BB.separated_cp) vs ideal 1 - 4 sin^2; C_D(Re) schematic labelled qualitative | #98 |
| C11 | Drag crisis: critical Re, separation moving aft (82 deg to 125 deg), C_D collapse; the base of sports-ball dynamics | 9.8 | A | CORE | – | the consequence of transition on bluff bodies; explains swing and Magnus sign changes | #99 |
| N82 | Three counter-intuitive points: singular limit nu -> 0, symmetric problem with asymmetric solution, roughness that lowers drag | 9.8 | C | NOTE | C11 | named; d'Alembert resolved by C01 and C09; Ch. 11, Ch. 14 | #100 |
| N83 | Sphere: doughnut eddy, loops shed above Re ~ 130, no regular street; crisis at Re_cr ~ 5e5 (book; the code default is 3e5, an illustrative argument — flag in the regime table), Fig. 9.22 | 9.9 | B | NOTE | C11 | reuse core.similarity.sphere_drag_coefficient('morrison') with Stokes and Oseen asymptotes; slip R13 flagged | #101 |
| N84 | Cricket-ball swing: seam trips one side, side force, y = a t^2/2 | 9.9 | B | NOTE | C11 | stated + BB.ball_swing_deflection + schematic; a-D25 demoted; numbers private | #102 |
| N85 | Tennis ball spin: negative and positive Magnus effect, Robins effect (slip R10) | 9.9 | B | NOTE | C11 | truth-table figure (BB.magnus_sign) + 'book prints Re < Re_cr twice, second is Re > Re_cr' callout | #103 |
| N86 | Baseball: curveball and knuckleball | 9.9 | C | NOTE | C11 | named; ch06 Magnus and Ch. 14 | #104 |
| N87 | Free 2-D laminar jet setup: boundary-layer approximation, dp/dx = 0, entrainment | 9.10 | B | NOTE | C12 | paragraph + geometry sketch | #105 |
| N88 | Eq. (9.53) u -> 0 as y -> +-infinity | 9.10 | B | NOTE | C12 | stated; BC of the ODE solver | #106 |
| N89 | Eq. (9.54) v = 0 on the axis (symmetry) | 9.10 | B | NOTE | C12 | stated; f(0) = 0 | #107 |
| N90 | Eq. (9.55) inlet profile at x0 | 9.10 | B | NOTE | C12 | stated; forgotten downstream | #108 |
| N91 | Eq. (9.56) integral of 2u u_x + u v_y + v u_y over y (printed tau lacks 1/rho, slip R7) | 9.10 | B | NOTE | C12 | D15 step + unit-check callout | #109 |
| N92 | Eq. (9.57) d/dx of the momentum flux integral of u^2 dy = 0 | 9.10 | B | NOTE | C12 | D15 result | #110 |
| N93 | Eq. (9.58) integral of u^2 dy = J/rho | 9.10 | B | NOTE | C12 | stated; J the momentum flux per unit span | #111 |
| N94 | Eq. (9.59) ansatz psi = u0 delta f(eta), delta = sqrt(nu x/u0) | 9.10 | B | NOTE | C12 | stated; specialises (9.19) | #112 |
| N95 | Eq. (9.60) u = u0 f'(eta) | 9.10 | B | NOTE | C12 | stated | #113 |
| N96 | Eq. (9.61) J/rho = u0^2 delta integral of f'^2 d eta | 9.10 | B | NOTE | C12 | D16 step (change of variable) | #114 |
| N97 | Eq. (9.62) u0 = (J^2/(C^2 rho^2 nu x))^(1/3) proportional to x^(-1/3) | 9.10 | B | NOTE | C12 | stated with a number | #115 |
| N98 | Eq. (9.63) delta = (C rho nu^2 x^2/J)^(1/3) proportional to x^(2/3) | 9.10 | B | NOTE | C12 | stated | #116 |
| N99 | Eq. (9.64) psi and eta in closed x-dependence | 9.10 | B | NOTE | C12 | stated | #117 |
| N100 | Eq. (9.65) (9.18) written in psi | 9.10 | B | NOTE | C12 | D16 start line | #118 |
| N101 | Reduced ODE 3 f''' + f f'' + f'^2 = 0 (book skips the differentiation; sympy-checked) | 9.10 | B | NOTE | C12 | result of D16 (written out because the book omits it) | #119 |
| N102 | Eq. (9.66) f' -> 0 as eta -> +-infinity | 9.10 | B | NOTE | C12 | stated | #120 |
| N103 | Eq. (9.67) f'(0) = 1 | 9.10 | B | NOTE | C12 | stated | #121 |
| N104 | Eq. (9.68) f(0) = 0 | 9.10 | B | NOTE | C12 | stated | #122 |
| N105 | Eq. (9.69) first integrals 3 f'' + f f' = 0 and 3 f' + f^2/2 = 3 | 9.10 | B | NOTE | C12 | D17 step | #123 |
| N106 | Eq. (9.70) f = sqrt(6) tanh(eta/sqrt(6)) | 9.10 | B | NOTE | C12 | D17 step | #124 |
| C12 | Free 2-D laminar jet: sech^2 profile (9.71), widths x^(2/3), entrainment (9.73)-(9.75) | 9.10 | A | CORE | – | flagship no-wall boundary layer; conserved momentum flux and entrainment recur in Ch. 12, Ch. 13 plumes | #125 |
| N107 | Eq. (9.72) C = integral of sech^4 = 4 sqrt(6)/3 | 9.10 | B | NOTE | C12 | D17 step (sympy exact) | #126 |
| N108 | Mass flux m_dot = rho u0 delta 2 sqrt(6) grows with x (entrainment) | 9.10 | B | NOTE | C12 | stated; D17 step | #127 |
| N109 | Eq. (9.73) m_dot = (36 J rho^2 nu x)^(1/3) proportional to x^(1/3) | 9.10 | B | NOTE | C12 | stated + Bickley check | #128 |
| N110 | Eq. (9.74) v = -psi_x = -(1/3)(J nu/(C rho x^2))^(1/3)[f - 2 eta f'] | 9.10 | B | NOTE | C12 | D17 step | #129 |
| N111 | Eq. (9.75) entrainment velocity v/u0 -> -+ sqrt(6)/(3 sqrt(Re_x)) | 9.10 | B | NOTE | C12 | stated + arrows in the explainer | #130 |
| N112 | Eq. (9.76) half-width h_99 (printed 5.6152 is the 4 % point; 1 % gives 7.3319, slip R6) | 9.10 | B | NOTE | C12 | 'book prints 5.6152, correct is 7.3319' callout; D18 shows the arccosh(10) arithmetic | #131 |
| N113 | Re_x, Re_h99, jet instability (inflection, Ch. 11), turbulent jet spreading (book says Ch. 13, should be Ch. 12, slip R14) | 9.10 | B | NOTE | C12 | stated with cross-reference correction; Ch. 11, Ch. 12 | #132 |
| N114 | Laminar wall jet (Glauert 1956): slot along a wall, boundary layer plus free-jet outer part | 9.10 | B | NOTE | C13 | paragraph + sketch | #133 |
| R08 | Eq. (9.77) wall-jet wall conditions u = v = 0 | 9.10 | C | RECAP | C13 | = (9.12) and (9.13) (R04, R05); reminder in the wall-jet block | #134 |
| R09 | Eq. (9.78) u -> 0 as y -> infinity | 9.10 | C | RECAP | C13 | = (9.53) edge condition (N88) and (9.14); reminder in the wall-jet block | #135 |
| N115 | Wall-jet chain: integrate (9.18) from y to infinity, multiply by u, integrate over y | 9.10 | B | NOTE | C13 | D19 start | #136 |
| N116 | Eq. (9.79) integration by parts of the inner integral | 9.10 | B | NOTE | C13 | D19 step | #137 |
| N117 | Eq. (9.80) conserved 'flux of exterior momentum flux' (wall jet loses momentum flux to the wall) | 9.10 | B | NOTE | C13 | D19 result; BL.wall_jet_invariant constant-in-x check | #138 |
| N118 | Eq. (9.81) substitution u = u0 f' gives x u0^2 = const | 9.10 | B | NOTE | C13 | D19 step | #139 |
| N119 | Eq. (9.82) u0 = C x^(-1/2), delta = (nu x^(3/2)/C)^(1/2) proportional to x^(3/4) | 9.10 | B | NOTE | C13 | stated; exponent comparison with the free jet | #140 |
| N120 | Reduced ODE 4 f''' + f f'' + 2 f'^2 = 0 (book prints coefficient 1, slip R3) | 9.10 | B | NOTE | C13 | D20 result + 'book prints 1, correct is 4' callout; planted wrong variant | #141 |
| N121 | First integrals 4 f f'' - 2 f'^2 + f^2 f' = 0 and f^(-1/2) f' + f^(3/2)/6 = f_inf^(3/2)/6 (slip R4) | 9.10 | B | NOTE | C13 | D21 steps + callout for the missing f^(1/2) | #142 |
| C13 | Wall jet: conserved exterior-momentum flux (9.80), ODE 4f''' + f f'' + 2f'^2 = 0 and the implicit solution (9.83) | 9.10 | A | CORE | – | second similarity family with a hidden conservation law; wall/free-jet comparison teaches how a wall changes exponents (x^(3/4) vs x^(2/3)) | #143 |
| N122 | Eq. (9.84) m_dot proportional to x^(1/4) | 9.10 | B | NOTE | C13 | stated; entrainment comparison with x^(1/3) | #144 |
| N123 | Eq. (9.85) momentum-invariant constant fixes C; m_dot at one x fixes f_inf (slip R8 on the symbol) | 9.10 | B | NOTE | C13 | stated; JET.wall_jet_constants | #145 |
| N124 | Wall-jet entrainment velocity and Fig. 9.29 profile | 9.10 | B | NOTE | C13 | our figure of f, f'; stated | #146 |
| C14 | Secondary flow: the teacup, radial force imbalance in a thin layer | 9.11 | A | CORE | – | a boundary layer creating a cross-flow; the bridge to Ekman layers and spin-up in Ch. 13 | #147 |
| S01 | Exercises 9.1-9.28, Literature Cited, Supplemental Reading | Ex. | C | SKIP | – | pointer: exercises named only (analysis §1); answers stay private (R16); one bibliography line | #148 |

## 3. Section coverage

| § | Title | A | B | C | RECAP | SKIP |
|---|---|---|---|---|---|---|
| 9.1 | Introduction | C01 | N01, N02, N03, N04, N05, N06, N07, N08, N09, N10, N11, N12, N13, N14, N16 | N15, N17 | R01, R02, R03, R04, R05 | – |
| 9.2 | Boundary-Layer Thickness Definitions | C02 | N18, N19, N20, N21 | – | – | – |
| 9.3 | Boundary Layer on a Flat Plate: Blasius Solution | C03, C04 | N22, N23, N25, N26, N27, N28, N29, N30, N31, N32, N33, N34, N35, N36, N37, N38, N39, N40, N41, N42 | N24 | R06, R07 | – |
| 9.4 | Falkner-Skan Similarity Solutions | C05 | N43, N44, N45, N46, N47 | N48 | – | – |
| 9.5 | Von Karman Momentum Integral Equation | C06 | N49, N50, N51, N52, N53, N54, N55 | – | – | – |
| 9.6 | Thwaites' Method | C07 | N56, N57, N58, N59, N60, N61, N62, N63, N64, N65, N66, N67 | – | – | – |
| 9.7 | Transition, Pressure Gradients, and Boundary-Layer Separation | C08, C09 | N68, N69, N70, N71, N72, N73 | N74, N75 | – | – |
| 9.8 | Flow Past a Circular Cylinder | C10, C11 | N76, N77, N78, N80, N81 | N79, N82 | – | – |
| 9.9 | Flow Past a Sphere and the Dynamics of Sports Balls | – (covered by B/C items inside C11) | N83, N84, N85 | N86 | – | – |
| 9.10 | Two-Dimensional Jets | C12, C13 | N87, N88, N89, N90, N91, N92, N93, N94, N95, N96, N97, N98, N99, N100, N101, N102, N103, N104, N105, N106, N107, N108, N109, N110, N111, N112, N113, N114, N115, N116, N117, N118, N119, N120, N121, N122, N123, N124 | – | R08, R09 | – |
| 9.11 | Secondary Flows | C14 | – | – | – | S01 |

Notes. RECAP rows are listed in the RECAP column whatever their depth (B: R01, R04, R05; C: the other six). S01 (exercises and reading lists) closes the chapter and is listed with 9.11. §9.9 has no A item of its own: its four rows sit in the C11 block (drag crisis, then its consequences for spheres and balls).

## 4. Prerequisites needing primers (concept or tool | needed by | why it is not A/B/RECAP)
Numbering continues at **P200** (the designer assigns numbers). Already primed and only *reminded* (one sentence, no new primer): P25
partial derivative, P26/P98 Taylor, P27 definite integral, P31/P94 `solve_ivp`, P37 trapezoid, P38 product rule, P40 sympy, P42
separation of variables, P43 exponent rules, P44 integrating an ODE twice, P46 `np.where`/`np.select`, P49/P91 chain rule, P68 orders of
smallness, P106 substitution, P108 `brentq`, P109 differentiation under the integral, P13 log axes, P185 diffusion time, P188 two-length
anisotropic scaling, P189 Leibniz with a moving limit, P190 `cumulative_trapezoid`, P192/P193 implicit stepping and Crank–Nicolson,
P196 `solve_bvp`, P197 exponent matching, P198 dominant balance.

| Concept or tool | Needed by | Why it is not A/B/RECAP |
|---|---|---|
| Improper integral of a deficit to infinity and truncating the tail (∫(1−u/U)dy converges because 1−u/U decays like a Gaussian) | C02, C04, C06 | new tool; all thickness integrals are of this type |
| `scipy.integrate.simpson` and numpy-2 `np.trapezoid` (composite Simpson on a sampled profile; `np.trapz` no longer exists) | C02, C04 | extends P37; Simpson is new |
| Monotone interpolation `PchipInterpolator` (no overshoot) | C02, C07 | new Python tool: tabulated l(λ), H(λ), profiles |
| Chain rule when the similarity variable η = y/δ(x) depends on both x and y (ψ_x, ψ_y, ψ_yy of ψ = Uδ f(η)) | C03, C05, C12, C13 | extends P49/P91 to a variable that moves with x; the step the book skips in every reduction |
| Scaling symmetry of an ODE (f → λf(λη)) and the Töpfer trick (solve with f″(0) = 1, rescale) | C04, C13 | new idea; makes a boundary-value problem an initial-value problem |
| Shooting vs boundary-value solving for a third-order nonlinear ODE (sensitivity ~ e^{η²/4}) | C04, C05 | extends P196/P31 with the failure mode of a naive shoot |
| Continuation in a parameter and a fold (saddle-node): why the solver fails at n = −0.0904 and how f″(0) as the parameter passes it | C05 | new numerical idea |
| Inflection point: where the second derivative changes sign, and what that means for a velocity profile | C08 | new vocabulary (returns as Rayleigh's criterion in Ch. 11) |
| Finding a crossing on a sampled curve (`np.sign` differences, linear interpolation of the root) | C07, C08 | new tool for x_sep |
| Central-difference residual of a PDE (`np.gradient`, second order) to test a solution numerically | C01, C03, C06 | check the primers register; if absent, a two-line primer |
| Parabolic vs elliptic vs hyperbolic PDEs (marching in x like time) and the von Mises coordinate (ψ as independent variable) | C01 (N14) | extends "order of a PDE vs conditions" from ch08; classification is new |
| Control volume with a streamline as one side (no mass crosses it) | C02 (D04) | extends ch04 C04 CV momentum; the streamline-as-boundary trick is new |
| First-order linear ODE by an integrating factor (U_e⁶) | C07 (D09) | new tool inside the Thwaites derivation |
| Integrals of powers of sine: ∫sin⁵φ dφ by substituting c = cos φ | C07 (D11) | small new tool; the cylinder example |
| Hyperbolic functions sech, arccosh; d(tanh)/dη = sech² and 1 − tanh² = sech² | C12 | extends the tanh primer of ch07; sech and arccosh are new |
| Exact derivatives: spotting (f f′)′ = f′² + f f″ and integrating twice (the "look for the total derivative" move) | C12, C13 | extends P38 (product rule backwards) to ODEs |
| Integration by parts with a variable lower limit, ∫₀^∞ u (∫_y^∞ u² dy′) dy | C13 (D19) | new tool for the invariant (9.80) |
| Partial fractions and sympy `apart`; an implicit solution inverted with `brentq` | C13 (D21) | new tool for (9.83) |
| Complex position and the velocity of a point vortex; an infinite row of vortices (Σ over images gives a cotangent) | C10 (D14) | extends ch05 point vortices and ch06 complex potential; the row sum is new |
| Linear stability of a steady configuration: perturb, linearise, eigenvalues of a matrix (`np.linalg.eig`), growth = positive real part | C10 (D14), later Ch. 11 | check the primers register for "eigenvalues"; if absent, a primer here |
| Mass and momentum flux through a cross-section (∫ρu dy and ∫ρu² dy, entrainment = growth of the first) | C12, C13 | extends ch04 C04; entrainment is a new physical idea |
| Radial force balance in a swirling flow, centrifugal ρu_φ²/R vs pressure gradient, and why ∂p/∂z ≈ 0 across a thin layer | C14 | recap of ch04 C09 / ch08 C04 plus the thin-layer argument (new) |
| Reading a regime table: local Re_x, plate Re_L, cylinder Re on the diameter, jet Re_x (six Reynolds numbers) | C01, C04, C10, C12 | vocabulary primer; the book uses six |
| Data collapse: rescale the axes so curves for different x fall on one | C04 | reminder of ch08 C10 with the new variable η = y√(U/νx) |

## 4b. Derivations written out (parsed by tools: ID first, CORE id in a column, ★★★ for hard, explainer slugs backticked in the LAST column)
Analysis §2b numbers appear as (a-Dnn) in the Result column. Every D row stands for one `nb.derivation` block.
| ID | Result (Eq.) | CORE | Difficulty | Steps | Tools used | Traps | Shown in |
|---|---|---|---|---|---|---|---|
| D01 | scaled equations (9.7)–(9.8), boundary-layer equations (9.9)–(9.10), thickness δ/L ~ Re^(-1/2) (9.4) (a-D01) | C01 | ★★ | 14 | anisotropic scaling P188, chain rule P49, dominant balance P198, sympy P40 | printed (9.7) misses the squares (slip R1); the y-equation is never displayed by the book, we write it; dropping 1/Re terms only after checking they multiply O(1) quantities; pressure scale ρU², not μU/δ | notebook · `bl_scaling_thicknesses` |
| D02 | matching −(1/ρ)dp/dx = U_e dU_e/dx (9.11) (a-D02) | C01 | ★ | 5 | derivative of Bernoulli (ch04 C05), (9.10) | why Bernoulli holds just above the layer and the pressure is the same across it; sign of dp/dx for a decelerating stream | notebook |
| D03 | displacement thickness δ* = ∫(1−u/U_e)dy (9.16) (a-D03) | C02 | ★ | 6 | mass-flux equality, improper integral | comparing fluxes of the real and a zero-thickness layer; extending h → ∞ needs the tail to converge; sign of the streamline lift | notebook · `bl_scaling_thicknesses` |
| D04 | momentum thickness θ (9.17) and ρU²θ = ∫τ₀dx on a flat plate (a-D04; the book leaves the CV proof to an exercise) | C02 | ★★ | 9 | control volume with a streamline as top boundary (ch04 C04), improper integral | the mass flux through the top of the CV; dp/dx = 0 removes the pressure force; θ is a length, not an angle | notebook · `bl_scaling_thicknesses` |
| D05 | Blasius reduction (9.19)–(9.27): f‴ + ½ f f″ = 0, δ = (νx/U)^{1/2} (a-D05) | C03 | ★★ | 13 | chain rule for η = y/δ(x), separation of x- and η-dependence (ch08 D02/D03), sympy | differentiating η with respect to x (δ′); two terms cancel; the choice C = 2 only defines δ; the brackets must be proportional | notebook · `blasius_similarity_collapse` |
| D06 | Blasius numbers (9.30)–(9.33): f″(0), η₉₉, δ*, θ = 2f″(0), τ₀, Cf, C_D (a-D06) | C04 | ★★ | 10 | scaling symmetry of the ODE (Töpfer), `solve_ivp`, `brentq`, integration of f′(1−f′) | reading η₉₉ as 4.93 instead of 4.910 (slip R2); factor 2 between θ and f″(0); ∫x^(-1/2)dx = 2√L; one side of the plate (R5) | notebook · `blasius_similarity_collapse` |
| D07 | Falkner–Skan reduction (9.34)–(9.36): f‴ + ((n+1)/2) f f″ − n f′² + n = 0 (a-D08) | C05 | ★★ | 12 | chain rule with two x-dependences, power laws, sympy | η ∝ x^{(n−1)/2} y differentiates twice; where (n+1)/2 comes from; why the forcing term is +n | notebook · `falkner_skan_family` |
| D08 | von Kármán momentum integral (9.37)–(9.43) (a-D10) | C06 | ★★ | 11 | conservative form, Leibniz P109/P189, fundamental theorem, product rule | ∫u²dy and ∫U_eU_e′dy diverge separately: work on [0, h] and let h → ∞; τ(∞) = 0; last rearrangement written out | notebook · `thwaites_marching` |
| D09 | Thwaites closed form (9.47)–(9.50): θ²U_e⁶/ν = 0.45∫U_e⁵dx (a-D11) | C07 | ★★ | 10 | integrating factor U_e⁶, θ² = νλ/U_e′ substitution, sympy dsolve | the algebra from (9.47) to l = (2+H)λ + (U_e/2)(θ²/ν)′ (the book skips it); where 0.45 and 6.0 come from; U_e → 0 at a stagnation point | notebook · `thwaites_marching` |
| D10 | exact Falkner–Skan closure: λ = nI_θ², l = I_θ f″(0), H = I_δ/I_θ and L(λ) vs 0.45 − 6λ (a-D14, ours) | C07 | ★★ | 9 | Falkner–Skan integrals, `PchipInterpolator`, sympy | this closure is exact for similar flows but not Thwaites' cross-family fit; two separation criteria (−0.0681 vs −0.09) | notebook · `thwaites_marching` |
| D11 | Thwaites on the cylinder: U_e = 2U sin φ, λ(φ), separation at φ ≈ 103° against ≈ 82° observed (a-D27, ours) | C07 | ★★ | 8 | ∫sin⁵φ by c = cos φ, root finding | angle from the forward stagnation point; ideal-flow pressure, not the real one; U_e → 0 at φ = 0 | notebook · `thwaites_marching` |
| D12 | wall curvature μ u_yy = dp/dx, inflection points (9.51)–(9.52), separation as τ₀ = 0 (a-D15) | C08 | ★ | 6 | evaluating (9.9) at the wall, sign argument | u = v = 0 at the wall kills the convective terms; sign of u_yy near the edge; link to f‴(0) = −n | notebook · `falkner_skan_family` |
| D13 | pressure drag of the separated model: C_D,p = ½∮C_p cos φ dφ with C_p = 1 − 4 sin²φ up to φ_sep and a constant base value behind (ours; ch06 D09 is the ideal case) | C09 | ★ | 7 | integration of C_p over the circumference, Gauss–Legendre check | angle from the forward stagnation point (ch06 uses the downstream axis); ideal-flow limit φ_sep → 180° gives 0 (d'Alembert); base pressure is a model assumption | notebook · `cylinder_drag_crisis` |
| D14 | Kármán street stability: staggered double row stable only if cosh(πb/a) = √2, b/a = 0.2805 (a-D23; the book only states it) | C10 | ★★★ | 15 | complex velocity of a vortex row (cot sum), linearisation, eigenvalues | the sum over images converges only in pairs; both rows are counter-rotating; unstable growth rate goes to zero at the marginal spacing; the non-staggered row is unstable at every spacing | notebook · `karman_street_stability` |
| D15 | jet momentum flux conserved, (9.56)–(9.58): d/dx ∫u²dy = 0 (a-D16) | C12 | ★★ | 7 | product rule backwards, fundamental theorem | printed τ lacks 1/ρ (slip R7); [uv] and [τ] vanish because u, ∂u/∂y → 0 | notebook · `free_jet_similarity` |
| D16 | free-jet similarity (9.59)–(9.68): exponents from constant J, then 3f‴ + f f″ + f′² = 0 (a-D17) | C12 | ★★★ | 14 | chain rule with x^{-2/3} inside η, powers of x cancel, sympy | ψ has four derivatives (ψ_x, ψ_y, ψ_xy, ψ_yy); C is free and chosen to make the coefficients integers; check that every x-power cancels | notebook · `free_jet_similarity` |
| D17 | free-jet solution (9.69)–(9.75): sech² profile, C = 4√6/3, ṁ ∝ x^(1/3), entrainment (a-D18) | C12 | ★★ | 11 | exact derivative (f f′)′, tanh substitution, ∫sech⁴ by t = tanh | recognising 3f″ + f f′ as a total derivative; C₁ = 0 from f″ → 0; derivative of η with respect to x for v | notebook · `free_jet_similarity` |
| D18 | jet width (9.76): h₉₉ = 7.332 [Cρν²x²/J]^{1/3}, not 5.615 (a-D19) | C12 | ★ | 5 | arccosh, root of sech² = 0.01 | the book's 2.2924 is arccosh 5 (the 4 % point); √6 factor | notebook · `free_jet_similarity` |
| D19 | wall-jet invariant (9.79)–(9.81): d/dx ∫u(∫_y^∞u²dy′)dy = 0 ⇒ x u₀² = const (a-D20) | C13 | ★★★ | 13 | integration by parts with a variable limit, Leibniz, fundamental theorem | the ordinary momentum flux is *not* conserved (wall force); boundary terms; ∫_y^∞ removes ν u_y | notebook · `wall_jet_invariant` |
| D20 | wall-jet reduction (9.82): 4f‴ + f f″ + 2f′² = 0 (a-D21; the book prints coefficient 1) | C13 | ★★★ | 12 | chain rule with x^{1/4} and x^{-1/2} powers, sympy | coefficient of f‴ is 4 (slip R3); δ ∝ x^{3/4}; the planted wrong ODE leaves a residual | notebook · `wall_jet_invariant` |
| D21 | wall-jet integration (9.83)–(9.85): 4f f″ − 2f′² + f²f′ = 0, g = √(f/f_∞), partial fractions, f″(0) = f_∞³/72 (a-D22) | C13 | ★★★ | 15 | integrating factor f, partial fractions `apart`, implicit inverse `brentq`, scaling symmetry | separation of variables lacks f^{1/2} in the book (slip R4); the constant 4.29 in the tail; f_∞ is free by the symmetry | notebook · `wall_jet_invariant` |
| D22 | teacup: net radial force per volume ρ(u_e² − u²)/R > 0 inward in the bottom layer and the flow pattern (a-D26) | C14 | ★★ | 7 | radial balance (ch04 C09), thin-layer argument ∂p/∂z ≈ 0 | pressure is set by the inviscid flow, centrifugal force by the slowed fluid; sign of the difference; Ekman scaling is deferred to Ch. 13 | notebook · `teacup_secondary_flow` |

## 4c. Derivations demoted to statements (first column is the A parent in bold, e.g. **C20** — never a bare ID or a D id: the parser reads a bare first-cell ID as an item and blanks its tier)
| A parent | Analysis §2b item | Result stated (Eq.) | Stated in (B item) | Why not written out |
|---|---|---|---|---|
| **C04** | a-D07 far field of the Blasius profile | f′ − 1 ~ (1/η)e^{−η²/4}; v_∞ = 0.86 U/√Re_x | N33, N34 | asymptotic analysis (WKB-type) with no later use; checked numerically in the notebook (residual against the solved f′) |
| **C05** | a-D09 sign of wall curvature and separation member | f‴(0) = −n from (9.36) at η = 0; f″(0) = 0 at n = −0.0904 | N46, N47 | the sign is a one-line read-off, used again inside D12; the fold is numerical continuation (explained in words, not derived) |
| **C07** | a-D12 Example 9.1 | θ² = 0.45νx/U ⇒ θ = 0.671 √(νx/U), Cf = 0.656/√Re_x | N66 | arithmetic of a worked number, no new moves after D09; result printed from `example_9_1()` |
| **C07** | a-D13 Example 9.2 (diffuser) | λ(x/L) = −(0.45/4)[(1+x/L)⁴ − 1], separation at x/L = 1.8^{1/4} − 1 | N67 | one power-law integral after D09; result printed from `example_9_2()` |
| **C04** | a-D24 scaling estimate vs 0.664, F_D ∝ U^{3/2}, C_D = 2C_f(L) | Cf ~ 2/√Re against 0.664/√Re; F_D ∝ U^{3/2} | N16, N39, N40 | one comparison and one integral of x^{-1/2}, given in the Blasius block |
| **C11** | a-D25 cricket-ball path | y = ½(F/W) g (d/U)² | N84 | constant force gives a parabola (ch01 kinematics); numbers private |

## 5. Interactive explainers (5–10 + backup)
Nine explainers; each is tied to A items, and the ten pieces of the chapter that most reward manipulation are covered (C06 rides with C07 in
E4; C09 and C11 share E5). Every window shows its equations with numbers and has an Explain tab, a Code tab, a 4–8-step walkthrough and
≥ 3 check questions. Reference explainers (all in `interactive-viz` §5): FDV = `forced_damped_vibrations.html`, AFE =
`angular_frequency_explorer_1.html`, APS = `amplitude_phase_second_order_II_3.html`.

### E1 · bl_scaling_thicknesses
- A: C01, C02 (also shows N04 δ/L ~ Re^(-1/2), N05–N08, N14 marching, N16 crude estimate, N18 δ₉₉, N19 streamline lift, N20 θ, N21) · **Confusion removed:** "the boundary layer is thin, but *how* thin, and what is 'its thickness' when the profile fades smoothly into the outer flow?" — δ₉₉, δ* and θ are three answers, each for a purpose (position of the edge, lift of the outer streamlines, lost momentum), and none is the layer's 'real' thickness.
- **Why interactive:** sliding U, ν or x shows δ ∝ √(νx/U) and Re^(-1/2) in the same picture, while switching the profile shape (Blasius, linear, sine, exponential, cubic) shows that the three thicknesses and H = δ*/θ change *by different ratios*, which a single static figure cannot.
- **Stage:** (1) the plate: a growing layer drawn to scale with the δ₉₉ curve, the displaced-streamline line at y = δ* above it, and tracer streamlines that lift; (2) u/U_e against y at the chosen x with δ₉₉, δ* and θ marked and the shaded areas (1 − u/U) and u(1 − u/U); (3) (hidePortrait) δ/L, δ*/δ₉₉ and H versus x or Re.
- **Controls:** free-stream speed U (0.1–50 m/s) · kinematic viscosity ν (air/water presets) · station x (0.01–2 m) · profile shape (select chips) · transport (x sweeping as time).
- **Equations:** (9.4) $\bar\delta/L\sim Re^{-1/2}$, (9.16) $\delta^*=\int_0^\infty(1-u/U_e)dy$, (9.17) $\theta=\int_0^\infty\frac u{U_e}(1-\frac u{U_e})dy$, (9.30) $\delta_{99}=4.91\sqrt{\nu x/U}$ (root of $f'=0.99$; shows the R2 callout), (9.11).
- **Mirrors** `BL.boundary_layer_scales`, `BL.thicknesses`, `BL.blasius_fields`.
- derivations: D01, D03, D04 · depth features: explain (numbered: Re, scales, integrals with your numbers, H, regime reading), code, **linked views** (plate + profile + growth curves), **modes** (raw / scaled by δ), **presets** (air 1 m/s at 1 m, water, small Re near the leading edge), **status** ("Re_x ≪ 1: layer not thin" / "attached, thin"), transport · follows: AFE (linked views on one clock, modes, "right now" notes) and APS (Explain with numbered live steps) · **aha:** the three thicknesses are three integrals of one profile; for Blasius δ* ≈ 1.72, θ ≈ 0.66 in units of √(νx/U), always in the ratio H ≈ 2.59.

### E2 · blasius_similarity_collapse
- A: C03, C04 (also shows N22–N30 reduction steps, N31–N32 BCs and truncation, N33 tail, N34 v, N35–N40 numbers, N41 temporal comparison, N42 collapse figure, R06–R07 repeated conditions) · **Confusion removed:** "why does one ODE describe the whole plate, and what does 'similarity' mean?" — profiles measured at different x are the *same curve* once y is scaled by δ(x) = √(νx/U), so the whole flow is one function f′(η).
- **Why interactive:** the mode toggle raw ↔ rescaled makes the collapse happen; dragging x, U, ν shows the curves stretch in y but coincide in η; the truncation slider η_max proves f′ → 1 is insensitive; a static plot shows three profiles, not the mechanism.
- **Stage:** (1) dimensional profiles u(y) at three (draggable) stations x₁, x₂, x₃ with a shaded ghost of the outer flow; (2) the same in η = y√(U/νx) (collapsed) with f′ and f, marks at η₉₉, δ*, θ and the inflection at the wall; (3) (hidePortrait) f″(η) and the shear τ₀(x) = 0.332ρU²/√Re_x with the tangent line at the wall.
- **Controls:** U · ν · three station positions (or 'add a station') · view mode raw/rescaled/both · η_max (truncation) · shooting slider f″(0) (shows what happens if you guess wrong: the profile overshoots or never reaches 1).
- **Equations:** (9.19) ψ = U δ f(η), (9.26) $\delta=\sqrt{\nu x/U}$, (9.27) $f'''+\tfrac12 ff''=0$, (9.28)–(9.29), (9.31) $\tau_0=0.332\rho U^2/\sqrt{Re_x}$, (9.32) $C_f=0.664/\sqrt{Re_x}$, (9.33) $C_D=1.33/\sqrt{Re_L}$.
- **Mirrors** `BL.falkner_skan(0.0)` (Töpfer scaling and `solve_bvp`), `BL.blasius_constants`, `BL.blasius_fields`.
- derivations: D05, D06 · depth features: explain (η, f′, f″ with your numbers, then τ₀, C_f, F_D), code (shoot vs Töpfer), linked views, **modes**, **inspector** (click a profile point: η, f′, u, y with arithmetic), **presets** ('air 1 m/s at 1 m', 'thin water layer'), status ("guess too small / too large / converged f'(∞)=1"), 'right now' notes · follows: APS (Explain as numbered derivation with live numbers; crosshair readouts) · **aha:** the curves at different x differ only by a stretch in y; f″(0) = 0.332 is *one* number that fixes the wall shear, the drag and, through θ = 2f″(0), the momentum loss.

### E3 · falkner_skan_family
- A: C05, C08 (also shows N43–N47 ansatz, thickness, family, separation member with the R11 callout; N70 wall relation; N71, N72 inflection rules; N73 thickening) · **Confusion removed:** "how can favourable and adverse pressure gradients be *one* family, and what exactly makes a profile separate?" — n is a dial: n > 0 fuller profile, larger shear, no inflection; n = 0 Blasius; n < 0 an inflection appears and the wall shear falls to zero at n = −0.0904.
- **Why interactive:** the profile, its curvature at the wall and the shear f″(0) change continuously with n; the inflection point slides in from infinity, the wall shear reaches zero at exactly one n, and the solver's failure beyond it (fold) is visible — static curves for seven n cannot show the approach.
- **Stage:** (1) profiles f′(η) with the wall tangent (slope f″(0)) and the inflection dot; a ghost Blasius; (2) f″(0) versus n with a dot at the current n and the separation point marked; (3) (hidePortrait) the outer-flow shape U_e = a xⁿ and the pressure gradient along x with a colour for favourable/adverse.
- **Controls:** exponent n (−0.0904 … 4) · β = 2n/(n+1) readout · scale variable (η vs ½√(n+1) η, as the book's Fig. 9.7) · profile toggles (f′, f″, f‴) · branch (attached / second, reversed) shown as a toggle.
- **Equations:** (9.34) $\psi=\sqrt{\nu xU_e}f(\eta)$, (9.35) $-dp/dx=U_eU_e'=na^2x^{2n-1}$, (9.36) $f'''+\frac{n+1}2ff''-nf'^2+n=0$, (9.51)–(9.52), μ u_yy(wall) = dp/dx.
- **Mirrors** `BL.falkner_skan`, `BL.falkner_skan_state`, `BL.falkner_skan_separation`, `BL.wall_curvature`.
- derivations: D07, D12 · depth features: explain (n → β, f″(0), f‴(0) = −n check, inflection location), code, **linked views**, **presets** (Hiemenz n = 1, Blasius, n = 1/3 wedge, n = −0.05, separation n = −0.0904), **status** ("favourable: no inflection" / "adverse: inflection at η = …" / "separated at the wall: f″(0) = 0"), **terms** (bars: viscous curvature u_yy versus pressure term at the wall; they are equal), 'right now' notes · follows: FDV (regime status + interpretation) · **aha:** at the wall, viscous curvature equals the pressure gradient (μ u_yy = dp/dx), so the sign of n *is* the sign of the curvature, and separation is where the shear reaches zero.

### E4 · thwaites_marching
- A: C06, C07, C08 (also shows N49–N55 momentum-integral steps and closures, N56–N59 λ, l, H and the closure table, N60–N65 L(λ) and the fit, N66 Example 9.1, N67 Example 9.2, N70–N73 wall curvature; the FS-exact closure D10) · **Confusion removed:** "how can one integral of the outer speed predict wall shear and separation without solving the PDE?" — the momentum integral is exact, the closure is a one-parameter fit, and the whole marching reduces to λ(x) = (θ²/ν)dU_e/dx crossing a threshold.
- **Why interactive:** the reader chooses the body's outer flow (diffuser, cylinder, wedge, stagnation, retarded flow, or draws U_e), the θ, λ, H, τ₀ curves update along x, and the separation point jumps as the diffuser angle, initial θ₀ or closure changes; the two criteria (−0.09 Thwaites, −0.0681 exact Falkner–Skan) can be compared.
- **Stage:** (1) the body/duct with the layer thickness δ*, the outer streamline and a red separation mark; (2) U_e(x) and pressure gradient with adverse regions shaded; (3) λ(x) with the −0.09 line, θ(x) and τ₀(x) (linked to the same cursor x); a bar row (`terms`) for the two contributions to θ² (the ∫U_e⁵ term and the θ₀ term).
- **Controls:** outer flow (select) · diffuser rate / cylinder / wedge n · initial θ₀ (0 … few times Blasius) · closure (FS-exact / 'white') · cursor x (drag) · transport (marching from the leading edge).
- **Equations:** (9.43) $\frac1\rho\tau_0=\frac d{dx}[U_e^2\theta]+U_e\delta^*\frac{dU_e}{dx}$, (9.44) $\lambda=\frac{\theta^2}\nu\frac{dU_e}{dx}$, (9.45)–(9.46), (9.48) $L(\lambda)=2l-2(2+H)\lambda$, (9.49), (9.50) $\theta^2U_e^6/\nu=0.45\int_0^xU_e^5dx'+\theta_0^2U_0^6/\nu$.
- **Mirrors** `BL.thwaites`, `BL.outer_flow`, `BL.thwaites_l`, `BL.thwaites_H`, `BL.separation_point`.
- derivations: D08, D09, D10, D11 · depth features: explain (the integral of U_e⁵ with your numbers, θ, λ, l, H, τ₀, C_f, then the reading), code (`thwaites` with the current U_e), **transport**, **terms**, **presets** (Blasius = constant U_e, diffuser Example 9.2, cylinder, Howarth retarded flow, stagnation), **status** ("attached: λ = … > −0.09" / "separation predicted at x = … stop"), **inspector** (click x: arithmetic of θ, λ, l, τ₀) · follows: APS (linked windows, numbered live explanation) and `fid_formula_lab` (term bars) · **aha:** the layer 'remembers' the whole pressure history only through ∫U_e⁵; a decelerating flow drives λ negative and the wall shear to zero.

### E5 · cylinder_drag_crisis
- A: C09, C10 (regimes), C11 (also shows N76–N77 regime notes, N78 Strouhal, N80 wake regimes, N81 subcritical, N83 sphere, N84 swing, N85 Magnus, N68–N69 transition and plate drag curve) · **Confusion removed:** "why does drag *fall* when the flow gets faster?" — at Re ≈ 3×10⁵ the layer turns turbulent, separates later (82° → 125°) and the wake narrows.
- **Why interactive:** the Reynolds number is a log-scale dial from 0.1 to 10⁷ that moves a schematic flow (creeping, twin eddies, street, subcritical wake, supercritical narrow wake), the separation angle, the C_p distribution and C_D at once; the reader sees the causal chain rather than four disconnected figures; the roughness/turbulence toggle shows the crisis moving.
- **Stage:** (1) animated schematic flow past a cylinder for the current Re (attached eddies / street with a shedding clock / narrow wake), separation points marked; (2) C_p(φ) with the ideal 1 − 4 sin²φ curve, the separated model curve and the shaded area giving pressure drag; (3) (hidePortrait) C_D versus Re (log–log) with the current dot (labelled as a qualitative schematic where no dataset is used; the sphere curve is the tested Morrison correlation with the Stokes/Oseen lines), and the Strouhal frequency for a chosen diameter.
- **Controls:** Re (log slider) · cylinder ↔ sphere mode · surface (smooth / rough or turbulent free stream) · base-pressure coefficient (optional) · diameter and speed for the shedding frequency f = St U/d · transport (time for the shedding animation).
- **Equations:** (4.102) $St=\Omega d/U$, C_p = 1 − 4 sin²φ (ch06 ideal cylinder), the pressure-drag integral of the separated model (D13, ours), sphere C_D(Re).
- **Mirrors** `BB.separated_cp`, `BB.separated_pressure_drag`, `BB.cylinder_flow_regime`, `BB.shedding_frequency`, `core.similarity.sphere_drag_coefficient`.
- derivations: D13 · depth features: explain (Re → regime → φ_sep → pressure drag with your numbers), code, **modes** (cylinder / sphere), **linked views**, **presets** (creeping, twin eddies, street, subcritical, critical, supercritical, rough cylinder), **status** ("subcritical: C_D ≈ 1" / "critical: drag crisis"), 'right now' notes with a highlighted regime table (book's rounded thresholds flagged), transport · follows: AFE (modes, presets, 'right now' table) · **aha:** drag is the pressure the wake fails to recover, so it depends on *where the layer lets go*; transition delays separation and shrinks the wake.

### E6 · karman_street_stability
- A: C10 (also shows N78 St, N80 wake regimes) · **Confusion removed:** "why is the vortex street's shape fixed at b/a ≈ 0.28?" — only one spacing of a staggered point-vortex double row is (marginally) stable.
- **Why interactive:** dragging b/a (and the stagger offset) changes the linear growth rate live, and a small random kick to the row visibly grows or not; the eigenvalue spectrum in the complex plane collapses to the imaginary axis at exactly cosh(πb/a) = √2.
- **Stage:** (1) the double row of point vortices in a periodic strip, perturbed positions exaggerated, advected by the induced velocity (animation); (2) growth rate σ (max Re λ) versus b/a with the marginal dot and the current dot; (3) (hidePortrait) eigenvalue spectrum of the linearised system in the complex plane.
- **Controls:** b/a (0.05–1) · stagger offset (0 = non-staggered … 0.5 = staggered) · perturbation wavenumber (mode of the row) · transport (time).
- **Equations:** the conjugate velocity of a row of vortices along the x-axis, $\frac{dw}{dz}=\frac{\Gamma}{2ia}\cot\frac{\pi (z-z_0)}a$ (from ch05 point vortices, derived in D14), $\cosh(\pi b/a)=\sqrt2$, $b/a=\frac1\pi\cosh^{-1}\sqrt2=0.2805$, (4.102).
- **Mirrors** `BB.karman_street_ratio`, `BB.karman_street_growth`, `BB.karman_street_velocity`.
- derivations: D14 (★★★, all 15 steps) · depth features: explain (row velocity, linearisation, eigenvalue with your b/a), code, **linked views** (street + growth curve + spectrum), **transport**, **presets** (stable b/a, too tight, too loose, non-staggered), **status** ("unstable, growth 0.13 Γ/a²" / "marginal"), 'right now' notes · follows: FDV (phase-plane style state view + regime interpretation) · **aha:** the stagger is stable only at one shape, and that shape is the one nature shows behind cylinders.

### E7 · free_jet_similarity
- A: C12 (also shows N87–N113 including the R6/R7 callouts, N108–N111 mass flux and entrainment) · **Confusion removed:** "a jet keeps its momentum but spreads and slows — how do these fit, and where does the extra mass come from?" — J = ρ∫u²dy is constant, so as δ grows u₀ falls; the flow entrains ambient fluid.
- **Why interactive:** moving x and J changes u₀ ∝ x^(-1/3), δ ∝ x^(2/3) and ṁ ∝ x^(1/3) together; the momentum-flux integral stays fixed as the profile changes, and entrainment arrows lengthen toward the axis; static profiles hide the invariant.
- **Stage:** (1) the 2-D jet in x–y with the sech² envelope, streamlines curving inward (entrainment) and h₉₉ lines; (2) u(y) at several x, raw and rescaled u/u₀ vs η (collapse to $\mathrm{sech}^2(\eta/\sqrt6)$); (3) (hidePortrait) u₀(x), h₉₉(x), ṁ(x) with slopes −1/3, 2/3, 1/3 on log axes.
- **Controls:** momentum flux J (N/m) · viscosity ν (air/water) · station x · exponent toggle 'free jet vs a wrong exponent' (to show why 2/3 and −1/3 are forced) · level for the half-width (1 % … 50 %, showing R6).
- **Equations:** (9.57) $\frac d{dx}\int u^2dy=0$, (9.59), (9.62)–(9.64), 3f‴ + f f″ + f′² = 0, (9.70) $f=\sqrt6\tanh(\eta/\sqrt6)$, (9.71) $u=u_0\,\mathrm{sech}^2(\eta/\sqrt6)$, (9.72) $C=\frac{4\sqrt6}3$, (9.73) $\dot m=(36J\rho^2\nu x)^{1/3}$, (9.75), (9.76) with the correction.
- **Mirrors** `JET.free_jet`, `JET.free_jet_halfwidth`, `JET.free_jet_mass_flux`, `JET.free_jet_entrainment_velocity`, `JET.jet_momentum_flux`.
- derivations: D15, D16 (★★★), D17, D18 · depth features: explain (J, C, u₀, δ, ṁ with your numbers, then h₉₉ at your level), code, **linked views**, **inspector** (click a point: η, f′, u, v arithmetic), **presets** (air slot jet, water jet, 1 % vs 4 % half-width), **status**, 'right now' notes with Re_x and instability hint (inflection ⇒ Ch. 11) · follows: APS (numbered derivation with live numbers) · **aha:** momentum flux is conserved, mass flux is not: the jet grows by taking in ambient fluid, and exponent −1/3 and 2/3 are forced by that one conservation law.

### E8 · wall_jet_invariant
- A: C13 (also shows N114–N124 incl. R3/R4/R8 callouts; free-jet contrast from C12) · **Confusion removed:** "in a wall jet the wall drags momentum out, so what *is* conserved, and why a different exponent (3/4) from the free jet (2/3)?" — the exterior momentum flux ∫u(∫_y^∞ u²)dy is invariant, giving x u₀² = const.
- **Why interactive:** in the free-jet mode ∫u²dy is constant with x; in the wall-jet mode it decreases (wall shear) but the double integral remains constant — the bars show both, x-dependent, at the same instant; the f_∞ scaling slider shows the symmetry f → λf(λη).
- **Stage:** (1) wall jet in x–y with the boundary-layer inner region and outer jet region; (2) f′(η) and f(η) with the implicit solution (9.83) and the IVP integration overlay; (3) (hidePortrait) bars: ∫u²dy, the invariant, ṁ versus x, for the wall and free jets.
- **Controls:** mode (free / wall) · Ψ or C (slider) · f_∞ · x · a toggle 'ODE as printed (coefficient 1)' that shows the residual (R3).
- **Equations:** (9.80) $\frac d{dx}\int_0^\infty u\Big(\int_y^\infty u^2dy'\Big)dy=0$, (9.82), 4f‴ + f f″ + 2f′² = 0, (9.83), (9.84) ṁ ∝ x^{1/4}, (9.85).
- **Mirrors** `JET.wall_jet_profile`, `JET.wall_jet`, `JET.wall_jet_invariant`, `JET.wall_jet_ode_solve`, `JET.wall_jet_constants`.
- derivations: D19, D20, D21 (all ★★★) · depth features: explain (invariant, C, f_∞, u₀, δ, ṁ with your numbers), code, **modes** (free vs wall), **terms** (bars for the two integrals), **linked views**, **presets** (small/large f_∞, printed-ODE check), **status** ("invariant conserved to 1e-9"), 'right now' notes · follows: AFE (modes) + APS (derivation with live numbers) · **aha:** when the obvious conserved quantity is lost to the wall, the next-order one (a momentum flux of the momentum flux) survives, and it fixes the exponents.

### E9 · teacup_secondary_flow
- A: C14 (also shows the Ekman hook) · **Confusion removed:** "why do tea leaves collect at the centre, when the water spins outward?" — the water near the bottom is slowed by friction but the pressure gradient set by the faster water above is unchanged, so the pressure wins over the weakened centrifugal force and pushes fluid inward.
- **Why interactive:** dragging the height z shows the imbalance ρ(u_e² − u²)/R growing toward the bottom; changing the swirl profile shape or layer thickness shows how the inflow depends on friction; a static picture cannot show the profile-by-profile force balance.
- **Stage:** (1) the cup in cross-section: swirl arrows, the meridional loop (in along the bottom, up the axis, out along the top), leaves drifting inward; (2) profiles u(z)/u_e (illustrative exponential, labelled), pressure-gradient force (uniform in z) and centrifugal force ρu²/R, with the difference shaded; (3) (hidePortrait) the net force versus z.
- **Controls:** swirl speed u_e · radius R · layer thickness δ · profile shape · height z · mode 'teacup' / 'river bend' (same balance across the channel, Exercise 9.28 named only).
- **Equations:** ∂p/∂R = ρu_φ²/R, ∂p/∂z ≈ 0 across the layer, net inward force per volume ρ(u_e² − u²)/R.
- **Mirrors** `ch09.secondary_flow_radial_force`.
- derivations: D22 · depth features: explain (the balance at your z with numbers), code, **modes**, **terms** (force bars that sum to the net inward force), **linked views**, **status** ("inflow: net force inward" above z, 'zero at the free stream'), 'right now' notes with the Ch. 13 hook (Ekman layer, spin-up) · follows: AFE (modes) + `fid_formula_lab` (term bars) · **aha:** a thin friction layer creates a force imbalance that turns a simple swirl into a three-dimensional circulation — the same mechanism that makes Ekman layers matter in the ocean and atmosphere.

### B1 · ball_swing_magnus (backup)
- A: C11 (also shows N84 cricket swing, N85 negative and positive Magnus, N86 baseball) · **Confusion removed:** "why does the seam or the spin change which way the ball curves?" — which side of the ball has its layer past the crisis decides the side force.
- Stage: ball in a stream with two sides' separation points, C_p bars, force arrow, truth table (side Re vs Re_cr) and the trajectory y = ½ a t². Controls: speed, spin, roughness, seam side. Mirrors `BB.magnus_sign`, `BB.ball_swing_deflection`. Derivations none · depth: modes, presets, status, terms. Built only if an explainer above fails review.

## 6. Python animations and interactive figures (A ID → what, why, player/figure kind)
Animations (6, `animate` + `show_animation`; ≤ 90 frames, FAST → 40):
| # | A | What moves | Why | Player |
|---|---|---|---|---|
| A1 | C01, C02 | the layer growing downstream: profile u(y) at increasing x under the δ₉₉ curve, δ*-displaced streamlines lifting; in the same frame the marching solver and Blasius | shows x as a time-like variable and the lift of the outer flow | video |
| A2 | C04 | dimensional profiles at x = 0.1 → 2 m stretch in y while the rescaled profile stays fixed | the collapse is the definition of similarity | video |
| A3 | C05, C08 | Falkner–Skan profile as n goes 4 → 0 → −0.0904; wall tangent rotating to zero slope | the approach to separation as continuous change | frames (stop at n = 0, −0.0904) |
| A4 | C07, C08 | Thwaites marching through the diffuser: λ falling, marker crossing −0.09 and τ₀ reaching zero | the meaning of "predicted separation" | frames |
| A5 | C10 | staggered vortex street (b/a = 0.2805) advecting; then an unstable spacing where a kick grows | stable vs unstable rows | video |
| A6 | C12, C13 | free jet and wall jet spreading side by side with the profiles at successive x | δ ∝ x^{2/3} vs x^{3/4}, and ṁ growth | video |

Interactive figures (plotly `slider_figure`, precomputed; ≤ 40 steps):
| # | A | What the slider controls | Why |
|---|---|---|---|
| F1 | C02 | profile shape parameter: δ₉₉, δ*, θ, H | the three thicknesses differ by shape |
| F2 | C04 | station x: dimensional and rescaled Blasius profile | published-page version of E2 |
| F3 | C05 | n: f′ and f″(0), inflection marker | published-page version of E3 |
| F4 | C07 | closure family: L(λ) computed on the FS family against the line 0.45 − 6λ | why the linear fit is good on one side and not the other |
| F5 | C09, C11 | separation angle: C_p(φ) and pressure drag | published-page version of E5 |
| F6 | C12, C13 | station x: free-jet and wall-jet profiles | exponent comparison |
| F7 | C08 | pressure gradient: the trio of Fig. 9.12 profiles from n = 1, 0, −0.05 with the sign of the wall curvature | the wall relation |
A live `ipywidgets` cell (kernel only) pairs with F3 for the free choice of n and Re. Static figures accompany every A block (Fig. 9.5–9.9, 9.12, 9.20, 9.22, 9.29 analogues, all generated by our code).

## 7. From-scratch moments
Each shows a hand-written version beside the tested function, followed by `assert np.allclose(...)`.
| A | §9.x | Hand-written | Compared with |
|---|---|---|---|
| **C01** | 9.1 | evaluate the terms of (9.1) on the Blasius field by central differences at three Re and show the ratio advective/viscous stays O(1) while the x-diffusion term scales as Re^(-1) | `BL.boundary_layer_scales`, `BL.bl_pressure_variation` |
| **C02** | 9.2 | trapezoid rule for δ*, θ on the sampled Blasius profile and on 1 − e^{−y/a} (δ* = a, θ = a/2) | `BL.thicknesses` |
| **C04** | 9.3 | RK4 shoot of f‴ = −½ f f″ from f″(0) = 1 and Töpfer rescale to get f″(0) = 0.33206 | `BL.falkner_skan(0.0)`, `BL.blasius_constants` |
| **C05** | 9.4 | `brentq` shoot for n = 1/3 and n = 1 (Hiemenz) in a truncated domain | `BL.falkner_skan` |
| **C06** | 9.5 | θ(x), δ*(x) from Blasius, finite-difference d(U²θ)/dx = τ₀/ρ | `BL.momentum_integral_residual` |
| **C07** | 9.6 | Thwaites in six lines with `cumulative_trapezoid` of U_e⁵ and the linear L(λ) for the diffuser | `BL.thwaites`, `ch09.example_9_2` |
| **C08/C09** | 9.7 | sign-change search for x_sep; trapezoid ∮C_p cos φ | `BL.separation_point`, `BB.separated_pressure_drag` |
| **C10** | 9.8 | eigenvalues of the 4×4 linearised two-vortex periodic cell using numerical differentiation of the induced velocity | `BB.karman_street_growth`, `BB.karman_street_ratio` |
| **C12** | 9.10 | integrate 3f‴ + f f″ + f′² = 0 with `solve_bvp`, compare with $\sqrt6\tanh(\eta/\sqrt6)$, ∫f′²dη = 4√6/3 | `JET.free_jet_ode_solve`, `JET.free_jet_constants` |
| **C13** | 9.10 | IVP for 4f‴ + f f″ + 2f′² = 0, f_∞³ = 72 f″(0), check (9.83) with `brentq` | `JET.wall_jet_ode_solve`, `JET.wall_jet_profile` |
| **C14** | 9.11 | two-line force balance in a sample layer profile | `ch09.secondary_flow_radial_force` |
§9.9 has no computable A item (its reuse of the tested Morrison correlation is a B call); §9.3 carries two (C03 has the sympy version of D05).

## 8. Notes for the implementer
Functions and helpers that the A items' figures and explainers need beyond analysis §4 (B/C functions optional):
1. **Explainer data.** JS explainers cannot call scipy: they need `BL.falkner_skan` for m ∈ [−0.0904, 4] as a tabulated f, f′, f″ on a fixed η-grid (`BL.falkner_skan_table(m_grid, n_eta)` → arrays, 6 s.f.), and the Blasius profile (m = 0) for parity rows; the JS side solves the ODE itself with `Viz.num.rk4Step` (Töpfer scaling for m = 0, shooting for |m| ≤ 1) and `selftest` compares at m ∈ {−0.05, 0, 1/3, 1}. Precomputed tables must be labelled as ours.
2. `BL.profile_shape(name, y, delta)`: closed forms (linear, sine, cubic, exponential, Blasius) with their exact δ*, θ, H for E1 and the primer checks.
3. `BL.blasius_fields` must accept arrays of x and y (meshgrid) for the collapse figure and A1–A2; `BL.blasius_delta99` uses the *computed* η₉₉ (4.910), with a `printed=True` variant giving 4.93 that a test must fail.
4. `BL.thwaites_cylinder(phi)` and `BL.thwaites_cylinder_closed_form(phi)`: D11 (λ(φ) and the crossing near 103°, angle from the forward stagnation point); `BL.outer_flow("retarded")` (Howarth linearly retarded flow).
5. `BB.cylinder_cd_schematic(Re)` (labelled qualitative; no dataset is fetched) and `BB.cylinder_state(Re, rough=False)` → dict(label, phi_sep_deg, St, cd_model). The sphere curve uses the tested `sphere_drag_coefficient("morrison")`; both curves need a "qualitative/tested" badge in the explainer.
6. `BB.karman_street_positions(t, b_over_a, eps, mode)` (linear evolution of a perturbed street for A5 and E6) and `karman_street_spectrum(b_over_a)` returning all eigenvalues.
7. `JET.free_jet_profile_table`, `JET.free_jet_at_level(level)` for E7; the 'printed' variants (`printed=True` → 5.6152 for h₉₉ at 1 % and coefficient 1 for the wall-jet ODE) are named tests that must fail.
8. `ch09.secondary_flow_layer_profile(z, delta, u_e, shape)` — illustrative swirl profile for E9, docstring labelled "illustrative, not a solution".
9. `ch09.book_slips()` returns a table: slip id, printed form, correct form, evaluator — used by the callouts (R1, R2, R3, R4, R6, R7, R10, R11, R14) and the tests.
10. Move `similarity_reduce_sympy`, `similarity_ode_solve`, `similarity_collapse_error` from ch08 to `core/similarity_reduce.py` (analysis §5) with the cases `blasius`, `falkner_skan`, `free_jet`, `wall_jet`, before the ch09 chapter module imports them.
11. Runtime: `FAST` uses 60 η-points, coarser Thwaites grids and 40-frame animations; cache `falkner_skan` solutions with `functools.lru_cache` keyed by rounded m; the closure table (60 solves) is built once per kernel. The whole notebook target is < 5 min on Colab.
12. Every function called by the notebook or an explainer parity row must be scalar-callable and take SI units; angles `phi_deg` from the forward stagnation point (say so in the docstring); `np.trapezoid`, never `np.trapz`.
13. Book values (private JSON): regime thresholds, angles, sports-ball numbers, jet constants, Table 9.1, exercise answers; printed with at most 4 significant figures in text; `check_public.py` must be run before every push.

Printed-slip callouts (which are shown as "book prints X, correct is Y"):
| slip | Where shown | Callout |
|---|---|---|
| S-1 (9.7) denominators | C01 / N08 | prints ∂²u*/∂x* and ∂²u*/∂y*; correct ∂x*², ∂y*² |
| S-2 (9.30) η₉₉ | C04 / N35 | prints 4.93; correct root 4.910 |
| S-3 wall-jet ODE | C13 / N120 | prints f‴ + f f″ + 2f′² = 0; correct 4f‴ + f f″ + 2f′² = 0 |
| S-4 wall-jet integral | C13 / N121 | prints ∫df/(f_∞^{3/2}f − f²); correct with f^{1/2} |
| S-5 (9.33) | C04 / N40 | one side of the plate (⚠️ confusion) |
| S-6 (9.76) | C12 / N112 | prints 5.6152; correct 7.3319 |
| S-7 (9.56) | C12 / N91 | prints τ with no 1/ρ; kinematic stress needed |
| S-8 (9.85) Ψ | C13 / N123 | Ψ is not a force per length (⚠️ confusion) |
| S-9 Reynolds numbers | C09 / N68, N69 | several Re_tr values (⚠️ confusion) |
| S-10 Magnus sentence | C11 / N85 | prints Re < Re_cr twice; second is Re > Re_cr |
| S-11 reverse flow | C05 / N47 | prints "solutions for n < −0.0904 with reverse flow"; correct: second branch for −0.0904 < n < 0, none below |
| S-12 Blasius name clash | C03 (and ch06) | ⚠️ confusion: boundary layer vs force theorem |
| S-13 sphere regimes | C11 / N83 | eddy range borrowed from the cylinder (⚠️ confusion) |
| S-14 "Chapter 13" for turbulent jets | C12 / N113 | correct: Chapter 12 |
| S-15 symbols λ, L, l, m | C07 | ⚠️ confusion: names |
| S-16 exercise answers | S01 | writers' rule only |

## 9. What the implementer must add beyond analysis §4 (module list)
Consistent with analysis §4 (42 functions + scripts): `fluidpy/core/boundary_layer.py` (BL: rows #1–#26 and #17–#24), `fluidpy/core/jets.py`
(JET: #32–#39), `fluidpy/core/bluff_body.py` (BB: #27–#31), `fluidpy/ch09_boundary_layers.py` (sympy engines `bl_nondim_sympy`,
`similarity_reduce_sympy` re-export, `momentum_integral_sympy`, `thwaites_sympy`, `jet_momentum_sympy`, `wall_jet_invariant_sympy`,
`wall_jet_sympy`, `example_9_1`, `example_9_2`, `secondary_flow_radial_force`, `book_slips`, `derive_all` + re-exports of every core name
for `nb.core` imports). Additions from this curation: items 1–9 of §8, `BB.pressure_drag_from_cp` for D13, `BL.thwaites_cylinder*`, and
the wall-jet double-integral helper `JET.wall_jet_invariant`. Scripts: the 19 `scripts/ch09_*.py` of analysis §4 row 43 plus
`ch09_karman_street.py` (A5), `ch09_thickness_shapes.py` (F1) and `ch09_free_vs_wall_jet.py` (A6). Evidence for A items follows analysis §6
(minimum: Blasius V2 V3 V5 V7; Falkner–Skan V2 V5; momentum integral V1 V4; Thwaites V1 V5 V6; free jet V2 V3 V4 V5; wall jet V2 V3 V4 V7;
Kármán ratio V1 V2), with planted wrong variants for slips S-1 to S-4, S-6, S-10 and S-11.
