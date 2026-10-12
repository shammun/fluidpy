# Chapter 14 — Aerodynamics: curation
(from `analysis/ch14.md` — 113 inventory rows (#1–#113), 27 equation numbers (14.1)–(14.27), 40 derivations in §2b (d-a … d-an;
★ 12, ★★ 21, ★★★ 7), 16 load-bearing clusters (§3), 10 suspected slips and 18 traps (§9); `book.yaml` ch14 (9 sections, 5 viz seeds)
and `policy.tier_a_full_treatment: [12, 18]`, `coverage: exhaustive`. Curated 2026-10-12, concept-curator. Read:
`knowledge/CUMULATIVE.md` (open items of ch13, rules in force), `concept_map.md`, `primers.md` (P01–P335 + P218a; **ch14 primers
start at P336**), `notation.md`, `viz_patterns.md` (the ch05, ch06, ch09 and ch13 "ideas for Ch. 14"), `knowledge/ch13.md` §8,
`analysis/ch13_curation.md` (format), the `interactive-viz` skill §4–§5. Equations are those the analyst transcribed from the page
images; no equation was re-typed from memory.)

Counts: A 16 · B 77 · C 20 · RECAP 12 · SKIP 2 · derivations written out 22 (★ 4 ★★ 14 ★★★ 4) · demoted to statements 11

Tier words (parsed by `tools/nbkit.py`): **CORE 16 · NOTE 83 · RECAP 12 · SKIP 2** = 113 rows; **DERIVATION 22**.
Depth: A 16 (all CORE) · B 77 (66 NOTE + 11 RECAP) · C 20 (17 NOTE + 1 RECAP + 2 SKIP: S01, S02).
**Reconciliation with analysis §2:** the 113 rows `#1…#113` each appear exactly once in §2 below (last column). A = 16 of 113 rows (14 %).

**The spine in one sentence:** *lift is circulation; the trailing edge chooses it; the wing tips make you pay for it* — C02–C07 are
the first two clauses, C08–C15 the third, C16 the first clause used for thrust.

**Decisions that shape this chapter.**
1. **Sixteen A items from the analyst's sixteen clusters, but not one per cluster** (two slots below the hard maximum are left
   free). (a) Cluster 1 (vocabulary) is **C01**, anchored on the aspect ratio $\Lambda\equiv s^2/A=s/\bar c$ (14.1) because it is the
   one piece of vocabulary a later result, $C_{D_i}=C_L^2/(\pi\Lambda)$ (14.26), cannot be read without; §14.1 and the section geometry
   of §14.3 sit inside it. (b) Cluster 3 is **C02** (the Kutta condition, shown on the flat plate with only ch06 tools and a Γ slider,
   the Kutta value found by a root-finder before any formula exists) and cluster 4 is **C03** (the starting vortex). (c) Cluster 5 is
   **C04**. (d) Cluster 6 is **C05**, anchored on $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\,\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9), with the
   ζ-plane flow as a B item whose derivation is written out because the book never writes it. (e) Cluster 7 is **one** block,
   **C06**, anchored on $C_L\simeq2\pi(\alpha+\beta)$ (14.12): the Kutta circulation $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11) is the
   boxed half-way line of its derivation D09, not a separate A item. (f) Cluster 2 is **C07** (surface speed and pressure — ours; it
   computes what Fig. 14.8 shows) and is taught *after* C06 because it needs the Kutta circulation; §14.3's pressure rows point
   forward to it. (g) Clusters 8, 9, 10 are **C08**, **C09**, **C10**. (h) Cluster 11 is split in two so that neither block carries
   more than two derivations: **C11** the sine series and its downwash (14.18)–(14.19), **C12** the lifting-line equation (14.20) and
   its solver; **C13** is the lift and induced drag from the coefficients (14.21)–(14.22). (i) Cluster 12 is **C14**. (j) Cluster 13
   (the horseshoe vortex) is **not** an A item: it is the tiny worked example of C09 and the building block of that block's
   from-scratch code, with ground effect and formation flight as B items. (k) Cluster 14 is **C15**, anchored on the drag breakdown
   and the polar; **stall is not an A item** — everything computable about it here is a qualitative model of ours (trap T16), which
   cannot carry the two evidence levels a CORE item needs, so the stall rows are B items of C07 (the adverse gradient), C06 (theory
   against measurement) and C15 (measured curves), and they get an explainer of their own (E5). (l) Clusters 15 and 16 are **one**
   block, **C16** "thrust from lift", anchored on the sail (viz seed 5) with the fish tail as its opening tiny example — the same
   velocity triangle in three systems, which is what the explainer's *modes* show.
   Next candidates if a reviewer asks for a seventeenth or eighteenth A (all B now): the fish tail (N77) as its own block; the
   Kutta circulation $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11) (N30) split from C06; the horseshoe vortex (N60).
2. **SEEN rows are RECAP** — the 12 rows marked SEEN are R01–R12, each reusing the earlier chapter's function and saying what this
   chapter adds (an off-centre circle, an inclined stream, the planform area as reference). No SEEN row is A: the Kutta–Zhukhovsky
   theorem $L=\rho U\Gamma$ (6.62) is R04 inside C02, the Zhukhovsky map $z=\zeta+b^2/\zeta$ (14.3) = (6.65) is R05 inside C04, the
   segment formula of ch05 D13 is R12 inside C09. Rows marked "NEW (… SEEN …)" are notes that open with a one-line reminder.
3. **SKIP is minimal** — S01 (exercises), S02 (literature). The analyst's 20 SKIP-class rows hold vocabulary, prose and history
   and are B or C notes here. Sweepback and shock waves (N76) are a C pointer to Ch. 15, not a SKIP.
4. **Material that is not in the book** — included, each labelled "**not in the book — ours**" in its block: the surface pressure
   computed from (14.9) (C07, an A item because it is the computed form of the book's Fig. 14.8 and the bridge to stall); the
   general circle parametrisation $\mu=b-ae^{-i\beta}$, $a=b(1+e)/\cos\beta$ (N25); the exact lift coefficient
   $8\pi(a/c)\sin(\alpha+\beta)$ (in C06); the finite trailing-edge speed (N29); the lumped-vortex impulsive start (C03's code); the
   span efficiency $e_{\text{span}}$ (N54); the finite-wing lift slope $K/(1+K/\pi\Lambda)$ (N59); the roll-up numbers (N37); the wake
   energy (N38); the parabolic polar and the best glide (C15, N74); the sail's limiting course (C16). Thin-airfoil theory is named
   only (N18), with an optional function.
5. **Slips and traps** (analysis §9: 10 slips, 18 traps) are tabled in §10 with the block that flags each. Two are taught, not just
   flagged: **slip #1** — the first two members of (14.19) are printed with the denominator $\cos\gamma_1-\cos\gamma$; D14 derives
   $\cos\gamma-\cos\gamma_1$ and a cell shows the printed form gives a downwash of the wrong sign, while the last member and
   everything downstream are right; **trap T9** — "induced drag is the largest part of the drag" becomes a computed share that
   depends on speed (N75).
6. **Two conventions on screen wherever a sign matters** (the "which p_o" pattern): the book's Γ is **clockwise**; functions take
   `Gamma_cw=` / `Gamma_ccw=`; every explainer readout says "Γ (clockwise)". The complex velocity is $dw/dz=u-iv$ (ch06), not
   ch13's $u+iv$. A conventions box precedes C01: the overloaded letters of trap T3 (b, a, β, c, e, Λ, γ, θ), the three θ's of T4,
   the axis switch of T5 (y is up in §14.3–§14.5 and **along the span** in §14.6, where w is positive **downward**), radians against
   degrees for K (T12), which area A is (T13). Span efficiency is `e_span` in code and $e_{\text{span}}$ in text, never bare e.
7. **Rule 9 (public repository).** Typical values the book quotes in words (stall angles, the slope per degree, camber percentage,
   the limiting sailing angle, Reynolds-number ranges, the fruit-fly coefficients, the aircraft mass) and every exercise input and
   answer stay in `tests/book_values_ch14.json`; worked numbers use our own inputs, visibly different (`ch14.illustrative_inputs()`).

## 1. Teaching order (A IDs grouped by book section, B/C IDs under each; one sentence each: "once you see X, Y follows")
Order = book order, with three placements to note: (i) the section geometry and the force definitions of §14.3 (N07, N08, R01, R02)
are taught in C01, so that C02 can start with the physics; (ii) C07 (surface pressure, a §14.3 picture) is taught after C06 because
it needs the Kutta circulation — C02 names the pressure coefficient (R03) and points forward; (iii) the exercises that hold results
(horseshoe, ground effect, formation flight, twist) are taught inside C09, C12 and C14 in the order the theory needs them. RECAP
rows are reminders inside the block named as their A parent. The conventions box (decision 6) precedes C01.

**§14.1 Introduction** — no A item of its own: N01 opens C01 (incompressible when $\mathrm{Ma}<0.3$; the three standing
assumptions; stall as the exception).

**§14.2 Aircraft Terminology**
- **C01** The wing in numbers: $\Lambda\equiv s^2/A=s/\bar c$ (14.1), $\bar c=A/s$ (14.2): once a wing is reduced to a span, an area
  and a section with an angle of attack, every later result is a statement about two coefficients, $C_L$ and $C_D$, and one shape
  number, Λ — and level flight, $L=W$, tells you which $C_L$ the wing must deliver.
  - B: N01, N02, N03, N05, N07, N08   ·   C: N04, N06, N83   ·   RECAP: R01, R02

**§14.3 Characteristics of Airfoil Sections**
- **C02** The Kutta condition (rear stagnation point at the sharp trailing edge ⇔ $\left.dw/d\zeta\right\rvert_{\zeta=b}=0$): once you
  see that ideal flow accepts *any* circulation and that only one value lets the flow leave the edge smoothly, the lift
  $L=\rho U\Gamma$ (6.62) stops being arbitrary.
  - B: N10   ·   C: N11   ·   RECAP: R04
- **C03** The starting vortex, $\Gamma_{ABCD}=0=\Gamma_{\text{foil}}+\Gamma_{\text{starting}}$: once the circulation round a large
  material circuit must stay zero (Kelvin), the foil can only acquire Γ by shedding −Γ — so viscosity, acting at the edge, makes
  the lift as well as the drag.
  - B: N14, N15, N16   ·   C: N17

**§14.4 Conformal Transformation for Generating Airfoil Shapes**
- **C04** The Zhukhovsky airfoil: the circle $\lvert\zeta-\mu\rvert=a$ through $\zeta=b$ mapped by $z=\zeta+b^2/\zeta$ (14.3), with
  $c\simeq4b$, camber $\simeq\tfrac12\beta c$, $t_{\max}/c\simeq1.3\,e$ (14.8): once you see that raising the circle bends the slit
  into an arc and shifting it back fattens the nose, two numbers (β, e) are a whole family of sections with a cusped tail.
  - B: N19, N20, N21, N22, N23, N24, N25, N26, N27   ·   C: N18   ·   RECAP: R05, R06, R07

**§14.5 Lift of a Zhukhovsky Airfoil**
- **C05** The flow round the airfoil, $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\,\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9): once the circle flow is
  known, the airfoil flow is the same complex potential read at mapped points, the speed rescaled by the map's stretching — which
  is infinite at the trailing edge unless the circle flow stagnates there.
  - B: N28   ·   RECAP: R08
- **C06** The lift of a Zhukhovsky airfoil, $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11) ⇒ $C_L\simeq2\pi(\alpha+\beta)$ (14.12): once
  the Kutta condition is written on the circle, the circulation is one line of geometry and the lift law — slope 2π per radian,
  shifted by the camber — follows.
  - B: N29, N30, N31, N32, N34   ·   C: N33, N35   ·   RECAP: R09
- **C07** Surface speed and pressure, $C_p=1-\lvert dw/dz\rvert^2/U^2$ on the contour (ours — Fig. 14.8 computed): once you can
  evaluate (14.9) on the circle, you can *see* the lift as the area of the pressure loop, and the steep pressure rise behind the
  suction peak that will stall the section.
  - B: N09, N12, N13   ·   RECAP: R03, R10

**§14.6 Elementary Lifting Line Theory for Wings of Finite Span**
- **C08** Bound vortex and trailing filaments, $d\Gamma_{\text{trail}}=-\dfrac{d\Gamma}{dy}\,dy$: once the circulation must fall to
  zero at the tips and a vortex cannot end in the fluid (Helmholtz), every drop in Γ leaves the wing as a streamwise filament — a
  vortex sheet that rolls up into the two tip vortices.
  - B: N36, N37, N39, N40   ·   C: N69   ·   RECAP: R11
- **C09** Downwash, $w(y_1)=\dfrac{1}{4\pi}\displaystyle\int_{-s/2}^{s/2}\dfrac{d\Gamma}{dy}\dfrac{dy}{(y_1-y)}$ (14.13): once each
  trailing filament is half an infinite line vortex, adding them along the span gives the downward velocity the wing flies in.
  - B: N41, N60, N65, N66   ·   C: N61   ·   RECAP: R12
- **C10** Induced drag, $L=\displaystyle\int_{-s/2}^{s/2}\rho U\Gamma\,dy$, $D_i=\displaystyle\int_{-s/2}^{s/2}\rho w\Gamma\,dy$ (14.15):
  once the local stream is tilted down by $\varepsilon\simeq w/U$, the Kutta–Zhukhovsky force tilts back with it and its
  streamwise part is a drag that exists in ideal flow.
  - B: N38, N42, N43, N44, N45   ·   C: N67
- **C11** The sine series and its downwash, $\Gamma=\sum\Gamma_n\sin(n\gamma)$ (14.18) ⇒
  $w(y_1)=\dfrac{1}{2s}\sum n\Gamma_n\dfrac{\sin(n\gamma_1)}{\sin\gamma_1}$ (14.19, last member): once the span is measured by an
  angle, $y=-(s/2)\cos\gamma$, the singular integral of (14.13) is done once and for all, term by term.
  - B: N49, N50, N51
- **C12** The lifting-line equation, $\dfrac K2Uc(\alpha+\beta)=\sum\left(1+\dfrac{nKc}{4s\sin\gamma}\right)\Gamma_n\sin(n\gamma)$
  (14.20): once each strip obeys the two-dimensional lift law at its *effective* angle, the unknown circulation appears on both
  sides, and N stations give N linear equations for N coefficients — any planform, any twist.
  - B: N46, N47, N48, N52, N63, N68   ·   C: N64
- **C13** Lift and induced drag from the coefficients, $L=\dfrac{\pi s}{4}\rho U\Gamma_1$ (14.21),
  $D_i=\dfrac\pi8\rho\sum n\Gamma_n^2$ (14.22): once the sines are orthogonal, the lift needs only the first coefficient and every
  other coefficient only adds drag.
  - B: N53, N54
- **C14** The elliptic loading, $C_{D_i}=\dfrac{C_L^2}{\pi\Lambda}$ (14.26): once only $\Gamma_1$ is kept, the downwash is uniform,
  $w=\Gamma_1/2s$ (14.24), the induced drag is the least possible for that lift and span, and it falls as the wing gets slender.
  - B: N55, N56, N57, N58, N59   ·   C: N62

**§14.7 Lift and Drag Characteristics of Airfoils**
- **C15** The drag breakdown and the polar, $C_D=C_{D0}+\dfrac{C_L^2}{\pi\Lambda e_{\text{span}}}$ (ours, from (14.26)): once drag is
  profile drag (friction + form) plus induced drag, one falls and the other rises with speed, and the best glide is where they
  are equal.
  - B: N70, N72, N73, N74, N75   ·   C: N71, N76

**§14.8 Propulsive Mechanisms of Fish and Birds** — no A item of its own: N77 (the fish tail, the tiny example that opens C16),
N79 (the bird's two strokes) and the pointers N78, N80.

**§14.9 Sailing against the Wind**
- **C16** Thrust from lift: $\mathbf U_r=\mathbf U-\mathbf V$, drive $=L\sin\beta_a-D\cos\beta_a$ (ours; the book resolves the force
  in words): once lift is perpendicular to the *relative* wind, moving the foil across the stream tilts the lift forward — a fish
  tail, a bird's downstroke and a sail are one velocity triangle, and a keel is the second wing that makes the sail's sideways
  force harmless.
  - B: N77, N79, N81, N82   ·   C: N78, N80

**Exercises and literature** — S01, S02 (pointer lines at the end of the notebook).

## 2. Chapter map (depth) — every inventory row exactly once
Rows are in analysis §2 order (inventory row number in the last column). IDs: `C` = A/CORE (teaching order), `N` = B or C NOTE
(inventory order), `R` = RECAP (inventory order), `S` = SKIP. The Item cell is the analyst's item with its equation as read from the
page image. No ASCII pipe inside cells. "slip #k" and "trap Tk" = analysis §9; "(§4c d-x)" = the analysis derivation is given as a
statement; module shorthand as in the analysis (AF = `fluidpy.core.airfoils`, LL = `core.lifting_line`, CONF = `core.conformal`,
POT = `core.potential`, BS = `core.biot_savart`, SIM = `core.similarity`, BL = `core.boundary_layer`, ch14 =
`fluidpy.ch14_aerodynamics`). "ours" = not in the book.

| ID | Item | § | Depth | Tier | A parent | Reason (A) / treatment (B) / pointer (C, SKIP) / source chapter (RECAP) | Row |
|---|---|---|---|---|---|---|---|
| N01 | Prose: aerodynamics = forces and moments on bodies; incompressible when the Mach number is small; three standing assumptions (no buoyancy, uniform constant-density stream, thin boundary layers); stall is the exception — $\mathrm{Ma}=U/c_{\text{sound}}<0.3$ ⇒ density changes $\approx\tfrac12\mathrm{Ma}^2$ (ours) | 14.1 | B | NOTE | C01 | stated as the opening of the chapter: the three assumptions as a checklist, each with the chapter that taught it (ch06 ideal flow, ch09 thin boundary layers) and where it fails (stall, N13); one number from `ch14.compressibility_check` for an airliner at take-off and at cruise (cruise fails the test → Ch. 15) | #1 |
| N02 | Def. fuselage, wing root, wing tip, wingspan s, chord c(y), planform area A, sweepback, taper — $A=\int_{-s/2}^{s/2}c(y)\,dy$ (ours, implied) | 14.2 | B | NOTE | C01 | stated with a labelled planform drawn by our code and the area integral; `LL.planform(kind, s, c_root, taper)` for rectangular, tapered and elliptic wings | #2 |
| C01 | Eq. (14.1) aspect ratio — $\Lambda\equiv s^2/A=s/\bar c$ | 14.2 | A | CORE | – | load-bearing vocabulary: Λ is the one shape number of the finite wing — $C_{D_i}=C_L^2/(\pi\Lambda)$ (14.26), the lift slope $K/(1+K/\pi\Lambda)$ and the polar of C15 cannot be read without it; the block also fixes $C_L$, $C_D$, α and $L=W$ for the whole chapter | #3 |
| N03 | Eq. (14.2) average chord — $\bar c=A/s$ | 14.2 | B | NOTE | C01 | stated beside (14.1) as the second form; one number for each of the three planforms; `LL.mean_chord` | #4 |
| N04 | Def. pitch, roll and yaw axes and their positive senses | 14.2 | C | NOTE | C01 | named with a three-axis sketch drawn by our code; pointer: roll returns with the ailerons (N06) and the rolling moment in a wake vortex (N64) | #5 |
| N05 | Def. drag opposes the flight direction, thrust is negative drag, lift is perpendicular to it; steady level flight — $L=W,\ T=D$ ⇒ $C_L=\dfrac{W}{\tfrac12\rho U^2A}$ (ours, implied) | 14.2, 14.3 | B | NOTE | C01 | stated as the block's worked number (the $C_L$ our example aircraft needs at two speeds) and derived in D01; `ch14.level_flight`; the power $DU$ is named for C15 | #6 |
| N06 | Def. control surfaces: aileron, elevator, rudder, flap, slot, spoilers, trim tabs, trimmed flight — all act by changing local camber | 14.2 | C | NOTE | C01 | named in one sentence each; pointer: "more camber ⇒ more lift at the same α" is N31 in C06, flaps and landing speed N73 in C15 | #7 |
| N07 | Def. airfoil section: rounded leading edge, sharp trailing edge, chord line, camber line, camber, thickness t, angle of attack α — camber $=\max y_c$; $t=\max(y_u-y_l)$ | 14.3 | B | NOTE | C01 | stated with an annotated section drawn by our code; `AF.section_geometry(x, y)` measures chord, camber and thickness of any closed contour (primer P337 for the contour split); the typical camber is a private book value — ours is illustrative | #8 |
| N08 | Def. foil-fixed frame; the aerodynamic force resolved into drag parallel and lift perpendicular to the oncoming stream — $D=F_x\cos\alpha+F_y\sin\alpha,\quad L=-F_x\sin\alpha+F_y\cos\alpha$ (ours) | 14.3 | B | NOTE | C01 | stated with the rotation written out and one number; `AF.resolve_lift_drag`; flags trap T5 early: when the stream is inclined to the x-axis (§14.4–§14.5) **lift is not $F_y$** | #9 |
| R01 | Rel. (4.107), (4.108) restated — $C_D\equiv\dfrac{D}{\tfrac12\rho U^2A},\quad C_L\equiv\dfrac{L}{\tfrac12\rho U^2A}$ | 14.3 | B | RECAP | C01 | ch04 §4.11, C15 — reminder with `SIM.lift_coefficient`, `SIM.drag_coefficient`; what is new: A is the **planform** area here (trap T13: ch09 used the frontal area), a section uses the chord per unit span (`AF.section_coefficient`, primer P336); the page's "Section 4.3" is §4.11 (slip #6) | #10 |
| R02 | Def. friction drag (surface shear stress) vs pressure drag (form drag); lift is almost entirely a pressure force | 14.3 | C | RECAP | C01 | ch09 C09 (form drag), ch04 surface forces — one reminder sentence; pointer: refined into friction + form + induced in C15 | #11 |
| R03 | Rel. pressure coefficient and the two ways of drawing it round a section (normal to the surface; against chordwise position with the negative values upward) — $C_p=\dfrac{p-p_\infty}{\tfrac12\rho U^2}$ | 14.3 | B | RECAP | C07 | ch06 (6.35), C06 — reminder with `POT.pressure_coefficient`; what is new: the inverted axis (suction upward) and the loop; both panels of the book's figure are drawn for our own section from `AF.zhukhovsky_surface` | #12 |
| N09 | Rel. force normal to the chord from the two surface pressures; its coefficient is the area of the pressure loop — $F_y=-\int_0^cp_u\,dx+\int_0^cp_l\,dx,\qquad C_y\equiv\dfrac{F_y}{\tfrac12\rho U^2c}=\oint C_p\,d\!\left(\dfrac xc\right)$ | 14.3 | B | NOTE | C07 | stated and derived in D10 (the book leaves it to an exercise); `AF.normal_force_coefficient`, `AF.lift_from_cp`; a number: loop area against $\rho U\Gamma$ for our section; $C_L$ differs from $C_y$ by the tilt α | #13 |
| R04 | Thm. Kutta–Zhukhovsky lift theorem restated, Γ clockwise — $L=\rho U\Gamma$ | 14.3 | B | RECAP | C02 | ch06 C10, D18, (6.62), explainer `blasius_kutta_contour` — reminder with `ch06.lift_per_span(rho, U, Gamma_cw=)`; what is new: the theorem leaves Γ free, and this chapter supplies it; trap T1 (Γ clockwise) in a two-row sign table | #14 |
| N10 | Prose: lift development = circulation development; Γ sets where the stream attaches and leaves; in a real flow the circulation is the net vorticity of the two boundary layers — $\Gamma=\oint u_s\,ds=\int_0^c(U_{\text{upper}}-U_{\text{lower}})\,dx$ (ours) | 14.3 | B | NOTE | C02 | stated in words with the integral; the number is computed later in C07 with `AF.surface_circulation` (forward pointer "worked with numbers in C07 below") | #15 |
| C02 | Thm. Kutta condition (Zhukhovsky hypothesis): a circulation develops that is just large enough to move the rear stagnation point to the sharp trailing edge — $\left.\dfrac{dw}{d\zeta}\right\rvert_{\zeta=b}=0$ | 14.3 | A | CORE | – | load-bearing: the one physical input ideal flow needs; (14.11), (14.12), (14.16) and all of lifting-line theory rest on it; ch06 left "uniqueness needs Kutta" open for this chapter | #16 |
| N11 | Prose: three ways to make and control lift — a sharp trailing edge, the angle of attack, the camber | 14.3 | C | NOTE | C02 | named as a three-item list; pointer: angle and camber become $\alpha+\beta$ in $C_L\simeq2\pi(\alpha+\beta)$ (14.12), C06; the attainable range is a private book value | #17 |
| N12 | Prose + exercise: a rounded leading edge keeps the suction-side layer attached and recovers the leading-edge suction that a plate of no thickness needs for zero drag — $N=L\cos\alpha,\ S=L\sin\alpha=\pi\rho U^2c\sin^2\alpha$ (ours) | 14.3 | B | NOTE | C07 | stated as the answer to "the pressure on a flat plate is normal to the plate, so why is the force normal to the stream?"; `AF.flat_plate_forces`; a number; result given, not derived (§4c d-s) | #18 |
| N13 | Def. stall: at high α the Kutta condition fails and the suction-side boundary layer separates ahead of the trailing edge — separation where $\tau_0=0$ | 14.3, 14.7 | B | NOTE | C07 | stated beside the computed pressure rise behind the suction peak; `ch14.suction_side_separation` (Thwaites on the ideal surface speed, ch09 C07) shows the separation point moving forward with α — labelled **qualitative** (trap T16: it cannot predict maximum lift); stall numbers are private | #19 |
| N14 | Fig. 14.9 impulsively started foil at three times — $\Gamma=0\ \to\ 0<\Gamma<\Gamma_{\text{Kutta}}\ \to\ \Gamma=\Gamma_{\text{Kutta}}$ for $t>c/U$ | 14.3 | B | NOTE | C03 | stated with three streamline panels of our own at $\Gamma/\Gamma_K=0,\ \tfrac12,\ 1$ and an animation of the rear stagnation point moving to the edge | #20 |
| N15 | Prose: why the stagnation point moves — the slow boundary-layer fluid cannot turn the sharp corner; the shed layer and the stagnation pressure push the point to the edge; Γ above the Kutta value is not observed | 14.3 | B | NOTE | C03 | stated as the mechanism paragraph between N14 and the circuit argument; trap T14 (a model of a viscous outcome, not a theorem of ideal flow) | #21 |
| C03 | Thm. starting vortex: the bound circulation Γ is matched by a shed vortex of −Γ; a new one is shed at every change of speed, incidence or camber — $\Gamma_{ABCD}=0=\Gamma_{\text{foil}}+\Gamma_{\text{starting}}$ | 14.3 | A | CORE | – | load-bearing: explains how an ideal-flow circulation can appear at all, closes the vortex loop that C08 needs (bound + trailing + starting), and is Kelvin's theorem — the reader's ch05 and ch13 tool — used on a wing | #22 |
| N16 | Prose: viscosity is responsible for the drag **and** for the circulation and therefore the lift | 14.3 | B | NOTE | C03 | stated as the punch line of §14.3, after the conserved total has been seen flat at zero | #23 |
| N17 | Prose: historical notes — Lanchester, Kutta, Zhukhovsky | 14.3 | C | NOTE | C03 | named in two sentences; pointer: Lanchester and Prandtl return in N69; one printed birth year is wrong (slip #5, correct year 1868 from a cited biography) | #24 |
| N18 | Prose: routes to the flow round a given section — numerical solution of Laplace's equation, thin-airfoil theory, Theodorsen's method — $C_L=2\pi\!\left(A_0+\tfrac12A_1\right)$ (ours, standard) | 14.4 | C | NOTE | C04 | named; pointer: thin-airfoil theory is reused in Ch. 15 (linearised compressible flow); optional `AF.thin_airfoil` gives a one-line check that a parabolic arc has zero-lift angle $-2h/c$, matching N31 | #25 |
| R05 | Eq. (14.3) Zhukhovsky transformation — $z=\zeta+\dfrac{b^2}{\zeta}$ | 14.4 | B | RECAP | C04 | ch06 C11, (6.65), D20, explainer `conformal_joukowski` — reminder with `CONF.joukowski`, `joukowski_inverse(branch="outside")`; what is new: circles that are not centred on the origin | #26 |
| R06 | Rel. the centred circle of radius b maps to the slit of length 4b — $z=be^{i\theta}+be^{-i\theta}=2b\cos\theta$ | 14.4 | B | RECAP | C04 | ch06 D20, (6.66) — one-line reminder and the first panel of the four-circle figure (§4c d-e) | #27 |
| N19 | Rel. image of a point of an off-centre circle written with its distance R from the origin — $z=Re^{i\theta}+\dfrac{b^2}{R}e^{-i\theta}$ | 14.4 | B | NOTE | C04 | stated as the start of D03; trap T4: this θ is the polar angle of ζ from the origin | #28 |
| N20 | Eq. (14.4) real and imaginary parts of the image — $x=(R+b^2/R)\cos\theta,\qquad y=(R-b^2/R)\sin\theta$ | 14.4 | B | NOTE | C04 | stated as step 1 of D03; one traced point with easy numbers | #29 |
| N21 | Eq. (14.5) R eliminated — $x^2\sin^2\theta-y^2\cos^2\theta=4b^2\sin^2\theta\cos^2\theta$ | 14.4 | B | NOTE | C04 | stated as steps 2–4 of D03 (the difference of squares the book skips); `ch14.circular_arc_sym` asserts it | #30 |
| N22 | Rel. cosine rule on triangle OQP, $QP=a=b/\cos\beta$, $OQ=b\tan\beta$ — $\dfrac{b^2}{\cos^2\beta}=R^2+b^2\tan^2\beta-2Rb\tan\beta\cos(90^\circ-\theta)$ | 14.4 | B | NOTE | C04 | stated as steps 5–6 of D03 with the triangle drawn (primer P340) | #31 |
| N23 | Eq. (14.6) the cosine rule simplified and joined to (14.4) — $2b\tan\beta\,\sin\theta=R-b^2/R=y/\sin\theta$ | 14.4 | B | NOTE | C04 | stated as steps 7–8 of D03 | #32 |
| N24 | Rel. the raised circle maps to a circular arc through $z=\pm2b$; height $2b\tan\beta$; end tangent at 2β to the chord — $x^2+(y+2b\cot2\beta)^2=(2b\csc2\beta)^2,\qquad\cos\beta=b/a$ | 14.4 | B | NOTE | C04 | stated as the result of D03 and checked point by point (`AF.circular_arc`: every mapped point lies on that circle); the arc is traced twice | #33 |
| N25 | Eq. (14.7) the circle shifted along the real axis and slightly larger than b maps to a symmetric airfoil — $a\equiv b(1+e),\qquad e\ll1$ | 14.4 | B | NOTE | C04 | stated with the third panel; our general parametrisation $\mu=b-ae^{-i\beta}$, $a=b(1+e)/\cos\beta$ (`AF.zhukhovsky_circle`) reduces to the book's two cases (trap T7; the function also accepts a, b, β) | #34 |
| C04 | Rel. the circle shifted both ways maps to a cambered airfoil with thickness — circle $\lvert\zeta-\mu\rvert=a$ through $\zeta=b$; contour $z=\zeta+b^2/\zeta$ | 14.4 | A | CORE | – | load-bearing: the chapter's body — C05, C06 and C07 are computed on it, and it turns "camber" and "thickness" into two numbers (β, e) the reader can drag; viz seed 1 | #35 |
| N26 | Eq. (14.8) shape rules for a thin, slightly cambered Zhukhovsky airfoil; maximum thickness near the quarter chord — $c\simeq4b,\qquad\text{camber}\simeq\tfrac12\beta c,\qquad t_{\max}/c\simeq1.3\,e$ | 14.4 | B | NOTE | C04 | stated as the results of D03 (camber) and D04 (chord, thickness — "can be proved" in the book); `AF.zhukhovsky_shape(exact=)` measures the contour: 1.3 is $3\sqrt3/4=1.299$ and the rule is first order in e (a "rule against measured" table; slip #4) | #36 |
| N27 | Def. the trailing edge of a Zhukhovsky airfoil is a cusp; with a finite edge angle the Kutta flow stagnates there, with a cusp it does not — $\dfrac{dz}{d\zeta}=1-\dfrac{b^2}{\zeta^2}=0$ at $\zeta=b$ ⇒ angles double | 14.4 | B | NOTE | C04 | stated and derived in D05 (the book states it with a figure); `AF.trailing_edge_angle`; a zoom on the cusp | #37 |
| R07 | Exercise: flow past a trailing edge of finite included angle is a corner flow with a stagnation point at the edge — $\phi=\dfrac{Ud}{n}\left(\dfrac rd\right)^{n}\cos n\theta$; included angle $2\pi(1-1/n)$ (ours) | 14.4 | B | RECAP | C04 | ch06 C09 (corner flows $w=Az^n$), explainer `complex_potential_corners` — reminder with `ch14.trailing_edge_wedge`: speed $\propto r^{n-1}\to0$ for $n>1$; used in the last two steps of D05 | #38 |
| N28 | Rel. the ζ-plane flow: uniform stream at angle α past the off-centre circle with clockwise circulation Γ; the same Γ and α hold in the z-plane — $w=U\!\left[(\zeta-\mu)e^{-i\alpha}+\dfrac{a^2e^{i\alpha}}{\zeta-\mu}\right]+\dfrac{i\Gamma}{2\pi}\ln\dfrac{\zeta-\mu}{a}$ (ours) | 14.5 | B | NOTE | C05 | stated as the starting point of the block and derived in D06 by shifting and rotating the centred flow of ch06 (the book writes only the centred case); `AF.zhukhovsky_flow`; "the circulation does not change under the map" is shown in D08 (trap T6) | #39 |
| R08 | Rel. chain rule for the complex velocity in the two planes — $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\dfrac{d\zeta}{dz}$ | 14.5 | B | RECAP | C05 | ch06 (6.69), `CONF.mapped_flow` — reminder; step 1 of D07; trap T2: $dw/dz=u-iv$ here, not ch13's $u+iv$ | #40 |
| C05 | Eq. (14.9) complex velocity of the z-plane from that of the ζ-plane; equal at infinity; singular at $\zeta=\pm b$ — $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\,\dfrac{\zeta^2}{\zeta^2-b^2}$ | 14.5 | A | CORE | – | load-bearing: every speed, pressure and force on the airfoil comes through it, and its singularity at $\zeta=b$ is exactly what the Kutta condition removes | #41 |
| N29 | Rel. at the trailing edge (14.9) gives an infinite speed unless the circle flow stagnates there; then 0/0 and a finite speed — $\left\lvert\dfrac{dw}{dz}\right\rvert_{B'}=U\,\dfrac ba\cos(\alpha+\beta)$ (ours) | 14.5 | B | NOTE | C06 | stated with the value and a sympy `limit` cell (primer P343); `AF.trailing_edge_speed`; a table of the speed at shrinking distances from the edge for Γ = Γ_K and Γ ≠ Γ_K; result given, not derived (§4c d-p); trap T8 | #42 |
| R09 | Eq. (14.10) tangential speed on the circle with clockwise Γ — $u_\theta=-2U\sin\theta-\dfrac{\Gamma}{2\pi a}$ | 14.5 | B | RECAP | C06 | ch06 (6.37), D10, explainer `cylinder_circulation_lift` — reminder with `ch06.cylinder_surface_speed`; what is new: θ is measured from the stream-aligned diameter through the centre Q (trap T4); the start of D09 | #43 |
| N30 | Eq. (14.11) the clockwise circulation the Kutta condition requires — $\Gamma=4\pi Ua\sin(\alpha+\beta)$ | 14.5 | B | NOTE | C06 | stated as the boxed half-way line of D09 (steps 1–5: why B sits at $\theta=-(\alpha+\beta)$); `AF.kutta_circulation`, checked against the root-finder of C02 | #44 |
| C06 | Eq. (14.12) lift coefficient of a Zhukhovsky airfoil; lift slope 2π per radian — $C_L=\dfrac{L}{\tfrac12\rho U^2c}\simeq2\pi(\alpha+\beta)$ | 14.5 | A | CORE | – | load-bearing: the chapter's two-dimensional result; (14.16)–(14.20) generalise it strip by strip and Ch. 15 compares its slope with the supersonic one; viz seed 3 | #45 |
| N31 | Def. absolute angle of attack $\alpha+\beta$; zero lift at $\alpha=-\beta$; more camber raises the whole lift curve — $\alpha_{L=0}=-\beta,\qquad\beta=2\,\dfrac{\text{camber}}{c}$ | 14.5, 14.7 | B | NOTE | C06 | stated with a number; `AF.zero_lift_angle`; trap T3: β is the circle's camber angle here and minus the zero-lift angle in §14.6 | #46 |
| N32 | Fig. 14.17 theoretical and measured lift curves of a cambered Zhukhovsky airfoil: equal slopes, an offset from the boundary layers' displacement thickness, and the measured curve stalls | 14.5 | B | NOTE | C06 | stated with our own figure: (14.12), the exact curve and `ch14.generic_lift_curve` on one axis (the measured curve is a labelled model, constants ours) | #47 |
| N33 | Prose: Zhukhovsky airfoils are not practical (a cusp cannot be built; largest camber at mid-chord); other maps exist | 14.5 | C | NOTE | C06 | named; pointer: the Kármán–Trefftz map (finite edge angle) is the standard relative, not developed in this book | #48 |
| C07 | (ours) surface speed and pressure on a Zhukhovsky airfoil from (14.9) with the Kutta circulation — $C_p=1-\dfrac{\lvert dw/dz\rvert^2}{U^2}$ on $z(\theta)$ | 14.3, 14.5 | A | CORE | – | load-bearing: it is the book's Fig. 14.8 made computable — the loop whose area is the lift, the suction peak, and the adverse gradient that sets stall (C15) and the critical Mach number (Ch. 15) | #49 |
| R10 | (ours) force on the airfoil by the Blasius integral and by integrating the surface pressure — $D-iL=\dfrac{i\rho}{2}\oint\left(\dfrac{dw}{dz}\right)^2dz$ (6.60) | 14.5 | B | RECAP | C07 | ch06 C10, (6.60) — reminder with `POT.blasius_force`; what is new: with the stream inclined, the force is perpendicular to the **stream**, $\lvert F\rvert=\rho U\Gamma$ and $D=0$ (`AF.zhukhovsky_force`, two routes that must agree; §4c d-r) | #50 |
| N34 | Exercise: the flat plate at incidence (centred circle) — $\Gamma=4\pi Ub\sin\alpha,\quad C_L=2\pi\sin\alpha$ | 14.5 | B | NOTE | C06 | stated as the special case $e=\beta=0$ of D09 (exact, no small-angle step) and the flow used in C02; `AF.flat_plate_flow` | #51 |
| N35 | Exercises: plot a cambered Zhukhovsky airfoil with streamlines; lift coefficient at a second incidence — $C_L(\alpha_2)=C_L(\alpha_1)+2\pi(\alpha_2-\alpha_1)$ | 14.5 | C | NOTE | C06 | named; pointer: `scripts/ch14_zhukhovsky_airfoil.py` draws it with our parameters | #52 |
| N36 | Prose: flow round the tips of a finite wing — outward below, inward above; the wake surface carries streamwise vorticity of opposite sign on the two half-wings | 14.6 | B | NOTE | C08 | stated with a rear-view sketch drawn by our code; the sheet strength is the jump of spanwise velocity (ch05 C14), quantified by the block's A result | #53 |
| N37 | Def. trailing vortices → tip vortices: each half of the sheet rolls up into one vortex — $\Gamma_{\text{tip}}=\Gamma_0$; (ours) spacing $b'=\dfrac{1}{\Gamma_0}\displaystyle\int\Gamma\,dy$ $(=\tfrac\pi4s$ for the ellipse), descent speed $\Gamma_0/2\pi b'$ | 14.6 | B | NOTE | C08 | stated with numbers for our aircraft (`LL.rolled_up_pair`, `ch14.tip_vortex_estimate`) and a cached roll-up animation (`ch14.wake_rollup`, **qualitative** past the first turn, trap T16); hook: contrails and wake-turbulence separation | #54 |
| N38 | Prose: the tip vortices lengthen as the aircraft flies, so energy is fed into them all the time — $D_iU=\dfrac{d}{dt}(\text{wake kinetic energy})$ | 14.6 | B | NOTE | C10 | stated as the second reading of induced drag; a cell checks $D_i$ of (14.15) against `LL.wake_kinetic_energy_per_length` to round-off; result given, not derived (§4c d-aj) | #55 |
| R11 | Thm. Helmholtz: a vortex filament cannot end in the fluid; the tip vortices are closed far behind by the starting vortex | 14.6 | B | RECAP | C08 | ch05 C01, §5.3, explainer `vortex_tubes_cannot_end` — reminder; what is new: the closed loop bound vortex + two trailing vortices + starting vortex (C03) | #56 |
| N39 | Thm. the three assumptions of lifting-line theory: large aspect ratio; one straight bound vortex of strength Γ(y); straight, non-interacting semi-infinite trailing filaments | 14.6 | B | NOTE | C08 | stated as a numbered box that every later block refers to ("assumption (iii)"), each with what it rules out | #57 |
| N40 | Def. bound vortex (moves with the wing, strongest at mid-span, zero at the tips) — $\Gamma(\pm s/2)=0$ | 14.6 | B | NOTE | C08 | stated with a Γ(y) curve; `LL.circulation_series` | #58 |
| C08 | Def. shed filament: where the bound strength drops by dΓ a filament of that strength leaves — $d\Gamma_{\text{trail}}=-\dfrac{d\Gamma}{dy}\,dy$ | 14.6 | A | CORE | – | load-bearing: the picture all of lifting-line theory is built on (a staircase of horseshoe vortices); (14.13) integrates exactly this strength | #59 |
| R12 | Rel. velocity induced on the lifting line at $y_1$ by the semi-infinite filament leaving at y: half that of an infinite line — $dw(y_1)=\dfrac{-(d\Gamma/dy)\,dy}{4\pi(y-y_1)}$ | 14.6 | B | RECAP | C09 | ch05 D13, (5.17), explainer `biot_savart_filament` — reminder with `BS.segment_speed(d, theta_a, theta_b, Gamma)` at angles 0 and 90°; the start of D12 (§4c d-u); the page's "y (< y₁)" wording is loose (slip #3) | #60 |
| C09 | Eq. (14.13) downwash at $y_1$ on the lifting line from the whole trailing sheet — $w(y_1)=\dfrac{1}{4\pi}\displaystyle\int_{-s/2}^{s/2}\dfrac{d\Gamma}{dy}\,\dfrac{dy}{(y_1-y)}$ | 14.6 | A | CORE | – | load-bearing: the velocity every strip of the wing flies in; (14.14)–(14.26) all contain it | #61 |
| N41 | Fig. 14.22 downwash along the flight direction — $w(x\to-\infty)=0,\quad w(0)=w,\quad w(x\to+\infty)=2w$ | 14.6 | B | NOTE | C09 | stated with our own curve (`LL.downwash_along_x`): far behind, each filament is infinite both ways (§4c d-ai) | #62 |
| N42 | Rel. the downwash turns the local stream down by a small angle — $\varepsilon=\tan^{-1}\dfrac wU\simeq\dfrac wU$ | 14.6 | B | NOTE | C10 | stated as step 1 of D13; the page prints "tan" for the inverse tangent (slip #2); `LL.induced_angle` | #63 |
| N43 | Eq. (14.14) effective angle of attack — $\alpha_e=\alpha-\varepsilon=\alpha-\dfrac wU$ | 14.6 | B | NOTE | C10 | stated as step 2 of D13, with a number in degrees; `LL.effective_angle` | #64 |
| N44 | Rel. each strip behaves as a two-dimensional section in the resultant stream — $dL=dL_e\cos\varepsilon\simeq\rho U\Gamma\,dy;\quad dD_i=dL_e\sin\varepsilon\simeq\rho w\Gamma\,dy$ | 14.6 | B | NOTE | C10 | stated as steps 3–6 of D13 with the force triangle drawn; `LL.strip_forces(exact=)` shows that $U_e\sin\varepsilon=w$ is exact | #65 |
| C10 | Eq. (14.15) lift and induced drag of the whole wing — $L=\displaystyle\int_{-s/2}^{s/2}\rho U\Gamma\,dy,\qquad D_i=\int_{-s/2}^{s/2}\rho w\Gamma\,dy$ | 14.6 | A | CORE | – | load-bearing: the "price of lift" — a drag in ideal flow; (14.22), (14.25), (14.26) and the polar of C15 are evaluations of it | #66 |
| N45 | Prose: U with Γ gives L normal to U; w with Γ gives $D_i$ normal to w; induced drag vanishes for infinite span | 14.6 | B | NOTE | C10 | stated as the two-line reading of (14.15); the book's "largest part of the drag" remark is handled in N75 | #67 |
| N46 | Eq. (14.16) section lift coefficient on a finite wing — $C_L=K\left[\alpha-\dfrac{w(y)}{U}+\beta(y)\right]$ | 14.6 | B | NOTE | C12 | stated as step 1 of D16; K per **radian** (trap T12; K is a required keyword in LL functions; the typical measured value is private); `LL.section_lift_coefficient` | #68 |
| N47 | Rel. circulation of a section from its lift coefficient — $C_L=\dfrac{\Gamma}{\tfrac12Uc}$ ⇒ $\Gamma=\tfrac12UcC_L$ | 14.6 | B | NOTE | C12 | stated as step 2 of D16; `AF.circulation_from_cl` | #69 |
| N48 | Eq. (14.17) circulation at a station; with (14.13) an integral equation for Γ(y) — $\Gamma(y)=\dfrac K2\,Uc(y)\left[\alpha-\dfrac{w(y)}{U}+\beta(y)\right]$ | 14.6 | B | NOTE | C12 | stated as step 3 of D16 with the remark that the unknown sits inside w (primer P345); `LL.lifting_line_residual` | #70 |
| N49 | Rel. angular span coordinate — $y=-(s/2)\cos\gamma,\qquad y_1=-(s/2)\cos\gamma_1$ | 14.6 | B | NOTE | C11 | stated with a semicircle sketch (γ = 0 port tip, π starboard tip); `LL.span_angle`; trap T3: γ is not ch05's sheet strength | #71 |
| N50 | Eq. (14.18) sine series for the circulation — $\Gamma=\displaystyle\sum_{n=1}^{\infty}\Gamma_n\sin(n\gamma)$ | 14.6 | B | NOTE | C11 | stated with the first four shapes plotted against y (every term vanishes at both tips; even n are antisymmetric); primer P346 | #72 |
| C11 | Eq. (14.19) downwash of the series — $w(y_1)=\dfrac{1}{2s}\displaystyle\sum_{n=1}^{\infty}n\Gamma_n\dfrac{\sin(n\gamma_1)}{\sin\gamma_1}$ (last member; the first two members, printed with $\cos\gamma_1-\cos\gamma$, should read $\cos\gamma-\cos\gamma_1$) | 14.6 | A | CORE | – | load-bearing: it turns the singular integral of (14.13) into a sum, which is what makes (14.20)–(14.26) possible; it carries the chapter's one consequential printed slip (#1) | #73 |
| N51 | Rel. the integral the book evaluates without showing it (Glauert's integral, a principal value) — $\displaystyle\int_0^\pi\dfrac{\cos n\gamma\,d\gamma}{\cos\gamma-\cos\gamma_1}=\pi\,\dfrac{\sin n\gamma_1}{\sin\gamma_1}$ (ours, standard) | 14.6 | B | NOTE | C11 | stated and derived in D15 (the book never writes it out); `LL.glauert_integral` by quadrature against the closed form for n = 0…6 | #74 |
| C12 | Eq. (14.20) one equation for the coefficients — $\dfrac K2\,Uc(\alpha+\beta)=\displaystyle\sum_{n=1}^{\infty}\left(1+\dfrac{nKc}{4s\sin\gamma}\right)\Gamma_n\sin(n\gamma)$ | 14.6 | A | CORE | – | load-bearing: the chapter's key method — one small linear solve gives the loading, lift and induced drag of any planform and twist; Ch. 15 and Ch. 16 reuse the solver | #75 |
| N52 | Rel. solution procedure: keep N terms and require (14.20) at N stations — $\mathsf A_{jn}\Gamma_n=r_j,\quad\mathsf A_{jn}=\left(1+\dfrac{nKc_j}{4s\sin\gamma_j}\right)\sin n\gamma_j$ (ours) | 14.6 | B | NOTE | C12 | stated with the matrix written for N = 3 and built line by line in the from-scratch cell; stations $\gamma_j=j\pi/(N+1)$ avoid the tips (trap T10); a convergence table in N; result given, not derived (§4c d-ab) | #76 |
| N53 | Eq. (14.21) lift from the first coefficient only — $L=\dfrac{\pi s}{4}\,\rho U\Gamma_1$ | 14.6 | B | NOTE | C13 | stated and derived in D17 (the book gives the result only); `LL.lift_from_series` | #77 |
| C13 | Eq. (14.22) induced drag: every coefficient adds a positive amount — $D_i=\dfrac\pi8\,\rho\displaystyle\sum_{n=1}^{N}n\Gamma_n^2$ | 14.6 | A | CORE | – | load-bearing: the sum of squares is the whole argument for the elliptic optimum and defines the span efficiency every drag polar uses | #78 |
| N54 | Thm. for a given lift the induced drag is least when all higher coefficients vanish — $e_{\text{span}}\equiv\dfrac{\Gamma_1^2}{\sum n\Gamma_n^2}\le1$, $C_{D_i}=\dfrac{C_L^2}{\pi\Lambda e_{\text{span}}}$ (ours) | 14.6 | B | NOTE | C13 | stated as the one-sentence consequence of (14.22) with coefficient bars for a rectangular wing; `LL.span_efficiency`; named `e_span` (trap T3) | #79 |
| N55 | Eq. (14.23) the optimum: an elliptic distribution of circulation — $\Gamma=\Gamma_1\sin(\gamma)=\Gamma_1\sqrt{1-(2y/s)^2}$ | 14.6 | B | NOTE | C14 | stated as steps 1–2 of D19; `LL.elliptic_circulation` | #80 |
| N56 | Eq. (14.24) downwash of the elliptic distribution: the same at every station — $w(y)=\Gamma_1/2s$ | 14.6 | B | NOTE | C14 | stated as step 3 of D19; checked by direct quadrature of (14.13) | #81 |
| N57 | Eq. (14.25) induced drag of the elliptic distribution — $D_i=\dfrac\pi8\rho\Gamma_1^2=\dfrac{2L^2}{\pi\rho U^2s^2}$ | 14.6 | B | NOTE | C14 | stated as steps 4–6 of D19; the second form says span, not area, is what counts; `LL.elliptic_induced_drag` | #82 |
| C14 | Eq. (14.26) induced drag coefficient — $C_{D_i}=\dfrac{D_i}{\tfrac12\rho U^2A}=\dfrac{C_L^2}{\pi(s^2/A)}=\dfrac{C_L^2}{\pi\Lambda}$ | 14.6 | A | CORE | – | load-bearing: the result the chapter is remembered for — why gliders, albatrosses and turbine blades are slender; the polar of C15 and finite-wing work in Ch. 15 start from it; viz seed 4 | #83 |
| N58 | Eq. (14.27) lift of a strip two ways ⇒ chord ∝ circulation: an elliptic planform gives the elliptic loading; a non-elliptic planform can do it with twist — $dL=\rho U\Gamma\,dy=C_L\,\tfrac12\rho U^2c\,dy$ | 14.6 | B | NOTE | C14 | stated as step 1 of D20; trap T11 (elliptic *loading* is not elliptic *planform*) as a callout; `LL.planform("elliptic", …)` | #84 |
| N59 | Rel. lift coefficient of an untwisted elliptic wing from that of its section — $C_{L,3D}=\dfrac{C_{L,2D}}{1+2/\Lambda}$ (K = 2π); $\dfrac{dC_L}{d\alpha}=\dfrac{K}{1+K/(\pi\Lambda)}$ (ours) | 14.6 | B | NOTE | C14 | stated and derived in D20 (an exercise in the book); `LL.finite_wing_lift_slope`, `LL.elliptic_wing`; a number: the slope at Λ = 6 against 2π | #85 |
| N60 | Exercise: the simplest wing — one rectangular horseshoe vortex of constant Γ — $w(0)=\dfrac{\Gamma}{\pi s},\qquad C_{D_i}=\dfrac{C_L^2}{2\pi\Lambda}$ (ours) | 14.6 | B | NOTE | C09 | stated as the block's tiny example (two half-lines at distance s/2, added by hand) and the building block of its from-scratch code; `LL.horseshoe_velocity`, `LL.horseshoe_wing`; why it seems to beat the ellipse (trap T10: the centre downwash is used for the whole span); result given, not derived (§4c d-ah) | #86 |
| N61 | Exercise: a parabolic circulation distribution — $\Gamma=\Gamma_0\left(1-(2y/s)^2\right)$ | 14.6 | C | NOTE | C09 | named with one line of code (`LL.downwash(0.0, Gamma_fn, s)`) comparing its mid-span downwash with the elliptic one; pointer: any Γ(y) is handled by the series of C11 | #87 |
| N62 | Exercise: an untwisted elliptic wing in level flight — induced drag and circulation part-way along each wing | 14.6 | C | NOTE | C14 | named as the pattern of the block's worked number, done for our own aircraft (the exercise's inputs are private); pointer: sea-level density from `core.statics` (ch01) | #88 |
| N63 | Exercise: a rectangular wing twisted elliptically — $\alpha_t=\dfrac{Kc\,\alpha_r}{4s},\quad\Gamma_1=\tfrac12KUc\,\alpha_r$ (ours) | 14.6 | B | NOTE | C12 | stated as a solver experiment: `LL.twist_for_elliptic_loading`, then `LL.solve_lifting_line(twist=…)` returns no higher coefficients; result given, not derived (§4c d-ak); the exercise writes the span as b (slip #9) | #89 |
| N64 | Exercise: rolling moment on a small aircraft flying along a large aircraft's tip vortex — $M=\displaystyle\int_{-s/2}^{s/2}\rho U\,y\,\Gamma(y)\,dy$ | 14.6 | C | NOTE | C12 | named with one number from `ch14.wake_vortex_rolling_moment` (our inputs); pointer: the reason for wake-turbulence separation at airports; the exercise's w has the opposite sign to the chapter's (trap T5) | #90 |
| N65 | Exercise: velocity field of a horseshoe vortex in its own plane; upwash outboard of the tips | 14.6 | B | NOTE | C09 | stated with a sign map of w from `LL.horseshoe_velocity` on a grid: downwash between the tips, upwash outside — where a following bird flies (the page's "Λ-pattern" is a V formation, slip #10) | #91 |
| N66 | Exercise: ground effect — an image horseshoe of opposite strength below the ground reduces the downwash — $w=\dfrac{\Gamma}{\pi s}\left[1-\dfrac{1}{1+(4h/s)^2}\right]$ at the centre (ours) | 14.6 | B | NOTE | C09 | stated with the image method recalled in a line (ch05 C13, ch06 C08) and a number at h = s/4; `LL.ground_effect_factor`, checked by two `horseshoe_velocity` calls | #92 |
| N67 | Exercise: tip-vortex core size and induced drag (Rankine cores) — energy per length $\propto\Gamma^2\left[\ln(s/\sigma)+\tfrac14\right]$ (ours) | 14.6 | C | NOTE | C10 | named beside the wake-energy reading (N38); pointer: `ch14.vortex_pair_energy_per_length` with the Rankine vortex of ch03 C14, (3.28) | #93 |
| N68 | Prose: lifting-line results are approximate; an elliptic loading is only approximately optimal; extensions to several wings, propellers, ground effect, bird formations | 14.6 | B | NOTE | C12 | stated as the "when it fails" paragraph of the block: low aspect ratio, sweep, a wake that rolls up, stall (the section law is linear) | #94 |
| N69 | Prose: Lanchester versus Prandtl — who conceived what, and when | 14.6 | C | NOTE | C08 | named in two sentences; pointer back to N17 | #95 |
| N70 | Fig. 14.25 generic measured $C_L$ and scaled $C_D$ against α: a linear range, lift at α = 0 from camber, a maximum, a drop at stall; drag rising roughly as α² — $C_L\simeq K(\alpha+\beta)$ below stall | 14.7 | B | NOTE | C15 | stated with our own figure from `ch14.generic_lift_curve`, `ch14.generic_drag_curve` — an explicit empirical model, every constant a required keyword, labelled **qualitative** (trap T16); the book's numbers are private | #96 |
| N71 | Fig. 14.26 a stalled section: separation near the leading edge, a wide wake, lift and drag comparable | 14.7 | C | NOTE | C15 | named with a sketch drawn by our code; pointer: the separation mechanism is N13 in C07 | #97 |
| N72 | Value: the maximum lift coefficient depends on the chord Reynolds number (laminar separation before transition at moderate Re; transition first at high Re) — $\mathrm{Re}_c=Uc/\nu$ | 14.7 | B | NOTE | C15 | stated as the same physics as the drag crisis (ch09 C11): a turbulent layer survives a steeper pressure rise; `ch14.chord_reynolds`, `ch14.stall_regime` (thresholds are required keywords; numbers private) | #98 |
| N73 | Rel. camber shifts the whole lift curve up; a rounded, cambered nose delays stall; flaps add camber ⇒ lower landing speed — $U_{\min}=\sqrt{\dfrac{2W}{\rho A\,C_{L,\max}}}$ (ours) | 14.7 | B | NOTE | C15 | stated with a number for our aircraft with and without flaps; `ch14.stall_speed`; the last step of D01 | #99 |
| C15 | Def. the drag breakdown: total = friction + pressure drag; pressure drag = induced + form drag; profile drag = friction + form — $C_D=C_{D0}+\dfrac{C_L^2}{\pi\Lambda e_{\text{span}}}$ (ours) | 14.7 | A | CORE | – | load-bearing: it joins the chapter's ideal-flow drag (14.26) to ch09's viscous drag in one formula that sets the speed, range and glide of every aircraft, bird and turbine blade | #100 |
| N74 | Rel. $C_{D_i}\propto C_L^2$ and $C_L\propto\alpha$ ⇒ $C_{D_i}\propto\alpha^2$; the best glide — $\left(\dfrac LD\right)_{\max}=\tfrac12\sqrt{\dfrac{\pi\Lambda e_{\text{span}}}{C_{D0}}}$ at $C_L=\sqrt{\pi\Lambda e_{\text{span}}\,C_{D0}}$ (ours) | 14.7 | B | NOTE | C15 | stated as the result of D21 (ours; the book says only "∝ α²"); `ch14.best_glide`, `ch14.drag_vs_speed` (the U-shaped curve) | #101 |
| N75 | Prose: "normally the induced drag is the major part of the total drag of a wing" | 14.6, 14.7 | B | NOTE | C15 | stated as a computed share against speed (`ch14.drag_breakdown`): exactly half at best glide, larger in climb, smaller at cruise — a cited figure for an airliner at cruise is about 37 % (trap T9) | #102 |
| N76 | Prose: at high subsonic speed local shock waves appear; sweepback reduces the velocity component normal to the leading edge; delta wings — $\mathrm{Ma}_n=\mathrm{Ma}\cos\Lambda_{\text{sweep}}$ (ours) | 14.7 | C | NOTE | C15 | named; pointer: normal shock waves are Ch. 15 §15.6 and thin airfoils in supersonic flow §15.10; the sweep angle shares the letter Λ with the aspect ratio | #103 |
| N77 | Fish propulsion: the tail moves sideways at V while the fish moves forward at U; the lift is normal to the relative velocity — thrust $=L\sin\theta$, rate $UL\sin\theta$; side force $=L\cos\theta$, rate $VL\cos\theta$; $V/U=\tan\theta$ | 14.8 | B | NOTE | C16 | stated as the tiny example that opens the block (the triangle with easy numbers: the two powers are equal because $V=U\tan\theta$); `ch14.oscillating_foil(U, V, L, D)`; with drag the efficiency drops below 1 (ours; trap T18: an identity of the triangle, not an efficiency prediction); result given, not derived (§4c d-al) | #104 |
| N78 | Prose: whales and dolphins move a horizontal tail vertically; leading-edge tubercles delay separation; flexible fins change camber during a stroke | 14.8 | C | NOTE | C16 | named; pointer: separation control is N13/N72; biological flows are Ch. 16 | #105 |
| N79 | Bird flight: on the downstroke the lift normal to the relative wind gives an upward force and thrust; on the steeper upstroke the wing does little | 14.8 | B | NOTE | C16 | stated with `ch14.oscillating_foil` called twice (downstroke, upstroke) and the mean thrust over a cycle; a mode of the explainer | #106 |
| N80 | Value: measured bird-wing geometry; a computed fruit-fly flight at low Reynolds number | 14.8 | C | NOTE | C16 | named; pointer: the quoted coefficients are private; low-Reynolds-number flight is where "thin boundary layers" (N01) fails | #107 |
| C16 | Sailing against the wind: true wind U, boat velocity V, apparent wind; sail lift normal and drag parallel to it; the drive vanishes below a limiting course angle — $\mathbf U_r=\mathbf U-\mathbf V$; drive $=L\sin\beta_a-D\cos\beta_a$, zero at $\tan\beta_a=D/L$ (ours) | 14.9 | A | CORE | – | load-bearing for the chapter's fourth objective: one velocity triangle explains the sail, the fish tail and the flapping wing — and the relative wind on a turbine blade; viz seed 5 | #108 |
| N81 | The keel as a hydrofoil: the hull points slightly off the course; the keel's lift balances the lateral force — closest course to the apparent wind $\beta_a=\varepsilon_A+\varepsilon_H$, $\varepsilon=\tan^{-1}(D/L)$ of sail and of hull (ours) | 14.9 | B | NOTE | C16 | stated as the last steps of D22 with two force triangles that close; `ch14.sailboat_course_limit`, `ch14.sailboat_equilibrium`; the leeway angle is a private book value | #109 |
| N82 | Tacking: a zig-zag course at a fixed angle to the wind reaches an upwind destination — velocity made good $=V\cos\theta$ (ours) | 14.9 | B | NOTE | C16 | stated with a path drawn by our code and `ch14.sail_forces` over the course angle | #110 |
| N83 | Exercise: a projectile with constant drag coefficient on a shallow trajectory — $x\simeq\dfrac{2m}{\rho AC_D}\ln\!\left(1+\dfrac{\rho AC_DV_o^2\cos\theta_o\sin\theta_o}{mg}\right)$ | Ex. | C | NOTE | C01 | named beside the drag coefficient (R01); pointer: optional `ch14.projectile_range_drag` with a `solve_ivp` check; A is the frontal area there (trap T13) | #111 |
| S01 | Exercises 14.1–14.19 | Ex. | C | SKIP | – | exercise text is never published; the results the chapter relies on are taught through N09 (14.3, 14.4), N34 (14.5), N24 (14.6), N12 (14.7), N35 (14.8, 14.9), R12 (14.10), N60 (14.11), N61 (14.12), N62 (14.13), N63 (14.14), N59 (14.15), N64 (14.16), N65 (14.17), N66 (14.18), N67 (14.19), R07 (14.2), N83 (14.1) | #112 |
| S02 | Literature cited and supplemental reading | – | C | SKIP | – | pointer line at the end of the notebook: the chapter's reference list in the book; our cited public sources are in analysis §8 | #113 |

## 3. Section coverage

| § | Title | A | B | C | RECAP | SKIP |
|---|---|---|---|---|---|---|
| 14.1 | Introduction | – (covered by the B item that opens the C01 block) | N01 | – | – | – |
| 14.2 | Aircraft Terminology | C01 | N02, N03, N05 | N04, N06 | – | – |
| 14.3 | Characteristics of Airfoil Sections | C02, C03, C07 (taught after C06) | N07, N08, N09, N10, N12, N13, N14, N15, N16 | N11, N17 | R01, R02, R03, R04 | – |
| 14.4 | Conformal Transformation for Generating Airfoil Shapes | C04 | N19, N20, N21, N22, N23, N24, N25, N26, N27 | N18 | R05, R06, R07 | – |
| 14.5 | Lift of a Zhukhovsky Airfoil | C05, C06 (and C07) | N28, N29, N30, N31, N32, N34 | N33, N35 | R08, R09, R10 | – |
| 14.6 | Elementary Lifting Line Theory for Wings of Finite Span | C08, C09, C10, C11, C12, C13, C14 | N36, N37, N38, N39, N40, N41, N42, N43, N44, N45, N46, N47, N48, N49, N50, N51, N52, N53, N54, N55, N56, N57, N58, N59, N60, N63, N65, N66, N68 | N61, N62, N64, N67, N69 | R11, R12 | – |
| 14.7 | Lift and Drag Characteristics of Airfoils | C15 | N70, N72, N73, N74, N75 | N71, N76 | – | – |
| 14.8 | Propulsive Mechanisms of Fish and Birds | – (covered by B/C items that open the C16 block) | N77, N79 | N78, N80 | – | – |
| 14.9 | Sailing against the Wind | C16 | N81, N82 | – | – | – |
| Ex. | Exercises, literature | – | – | N83 | – | S01, S02 |

## 4. Prerequisites needing primers (concept or tool | needed by | why it is not A/B/RECAP)
Numbering continues at **P336** (the designer assigns the final numbers). Already primed or taught and only *reminded* (one sentence,
no new primer): P15 `assert np.allclose`, P16 animate, P17 slider_figure, P18 show_viz, P26/P98 Taylor expansion, P27 definite
integral, P37 trapezoid rule, P38 product rule, P40/P117 sympy (`limit`, `series`, `simplify`), P45 Euler's formula, P49 chain rule,
P57 `np.linalg.solve`, P68 orders of smallness, P70 `np.arctan2`, P87 `scipy.integrate.quad`, P92 parametric curves, P94
`solve_ivp`, P100 small-angle approximation, P106 substitution in an integral, P108 `brentq`, P127 power of a force, P135 closed-loop
integral of an exact differential, P136 periodic trapezoid rule on a closed loop, P144 improper integral as a limit, P147 systems of
ODEs for interacting vortices and their invariants, P150 a limit with a product held fixed, P151 integrals of sines and cosines over
a full period, P153 the complex plane in numpy, P154 complex derivative and analytic functions, P155 complex logarithm and branch
cuts, P156 Cauchy's integral theorem, P157 Laurent series and residues, P159 complex square roots (the outside root), P162
collocation and the condition number, P164 Green's first identity, P170 `minimize_scalar`, P182 `np.interp`, P236 linear recurrence
with constant coefficients, P261 even and odd functions, P265 quotient rule, P272 completing the square into a circle, P275 Galerkin
truncation with orthogonal sines, P306 proof by induction.
**Earlier-chapter ideas that an A item uses and that are not inventory rows of this chapter** (each gets a 🔁 recap line where first
used, pointing to its chapter): Kelvin's circulation theorem $D\Gamma/Dt=0$ (5.8) (ch05 C04, D04–D05 → C03); Stokes' theorem and
the additivity of circulation over sub-circuits (ch02, ch03 → D02); point-vortex velocity and `BS.point_vortex_evolve` (ch05 C12 →
C03, N37); vortex sheet as a jump in tangential velocity (ch05 C14 → C08); Biot–Savart for a filament
$d\mathbf u=\dfrac{\Gamma\,d\mathbf l}{4\pi}\times\dfrac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ (5.17) (ch05 C08 → C09);
method of images (ch05 C13, ch06 C08 → N66); complex potential and $dw/dz=u-iv$ (ch06 C09 → C05); the cylinder with circulation and
its stagnation points (ch06 C06 → C02, C06); critical point of a map and angle doubling (ch06 D19 → D05); Bernoulli's equation
(ch04 → C07); Thwaites' method and separation under an adverse gradient (ch09 C07, C08 → N13); flat-plate friction (ch09 C04, ch12
→ C15's profile-drag estimate); Mach number and the speed of sound (ch01, ch04 C15 → N01); the Rankine vortex (ch03 C14 → N67).

| Primer (concept or tool) | Needed by | Why it is not A/B/RECAP |
|---|---|---|
| P336 force per unit span; section coefficient ($c$ per unit span) against wing coefficient (planform area A); dynamic pressure $\tfrac12\rho U^2$ as the force scale | C01 (R01), C06, C12 | bookkeeping that the book assumes; needed to read "$L=\rho U\Gamma$ is a force per unit span" and to pass from a section to a wing |
| P337 splitting a closed contour into upper and lower surfaces (`np.argmin` / `np.argmax` for the two edges, slicing, `np.interp` needs increasing x) | C01 (N07), C07 (D10) | a Python idiom behind `section_geometry` and the pressure-loop integral; no physics of its own |
| P338 radians and degrees for angles and for slopes (`np.deg2rad`; 2π per radian ≈ 0.110 per degree) | C06, C12, C15 | trap T12: the book writes K per radian in (14.16)–(14.20) and quotes it per degree in §14.7 |
| P339 discrete-vortex time stepping: one bound vortex, one new wake vortex shed per step, the total circulation as the invariant to watch | C03 | the numerical model of the starting vortex is ours; the reader has met point vortices (P147) but not a body that sheds them |
| P340 the cosine rule $c^2=a^2+b^2-2ab\cos C$ | C04 (D03) | plane trigonometry the book uses without comment |
| P341 cot, csc and the double-angle identities ($\tan\beta-\cot\beta=-2\cot2\beta$, $1+\cot^2x=\csc^2x$) | C04 (D03) | the last moves of the circular-arc derivation, skipped by the book |
| P342 binomial approximation with a complex small quantity, $1/(1+\epsilon)\approx1-\epsilon$, and keeping first order in e | C04 (D04) | Taylor (P26) and orders of smallness (P68) are primed for real numbers; here ε is complex |
| P343 l'Hôpital's rule for a 0/0 limit (`sp.limit`) | C06 (N29), C07 (the trailing-edge point), C11 ($\sin n\gamma/\sin\gamma\to n$ at the tips) | a calculus tool not used earlier in the project |
| P344 Cauchy principal value: an integrand that changes sign through an infinity inside the range; symmetric exclusion; `quad(..., weight="cauchy", wvar=)` and subtracting the singular part | C09 (D12), C11 (D14, D15) | new; an improper integral (P144) of a different kind — the limit exists only because the two sides cancel |
| P345 integral equation: the unknown function appears under an integral sign | C12 | a new kind of equation for the reader; one sentence and a two-line discrete analogue ($x=b-\mathsf Mx$) |
| P346 half-range sine series on $[0,\pi]$: any function that vanishes at both ends; orthogonality $\int_0^\pi\sin n\gamma\sin m\gamma\,d\gamma=\tfrac\pi2\delta_{nm}$ | C11, C13 (D17, D18) | P151 covers the full period and P275 the Galerkin use; the half range and the vanishing ends are what fit a wing |
| P347 product-to-sum identity $\cos(n+1)\gamma+\cos(n-1)\gamma=2\cos\gamma\cos n\gamma$ and a three-term recursion fixed by two starting values | C11 (D15) | the key move of Glauert's integral; the recurrence primer P236 is recalled for the solving step |
| P348 a double sum collapsed by orthogonality, $\sum_n\sum_ma_nb_m\delta_{nm}=\sum_na_nb_n$ | C13 (D18) | the product of two series is the step where novices drop the cross terms |
| P349 reading a drag polar ($C_L$ against $C_D$; the tangent from the origin is the best $L/D$) and the minimum of $a/x+bx$ | C15 (D21) | a way of plotting the reader has not met; the minimum is used twice (in $C_L$ and in speed) |
| P350 relative (apparent) velocity as a vector subtraction, $\mathbf U_r=\mathbf U-\mathbf V$, with the angle from `np.arctan2` | C10 (the tilted stream), C16 | frames were primed in ch03 (P96) for flows; here the *body* moves across the stream and the angle is the whole result |

15 primers, P336–P350. Glossed in a line where used (no primer): difference of squares (D03), maximum of a function along a curve
by setting the derivative to zero (D04), piecewise functions (N70), `np.linalg.cond` (P162 reminder), `np.unwrap`, f-strings.

## 4b. Derivations written out (parsed by tools: ID first, CORE id in a column, ★★★ for hard, explainer slugs backticked in the LAST column)
Analysis §2b labels appear as (d-x) in the Result column. Every D row stands for one `nb.derivation` block. "book never writes it
out" marks rule (b): the result belongs to a B item of the named A block and the book states it without derivation, leaves it to an
exercise, or does not state it (then also "ours"). Step counts are at one small move per step. All ★★★ rows get a sympy check
cell that re-runs the derivation's own construction.

| ID | Result (Eq.) | CORE | Difficulty | Steps | Tools used | Traps | Shown in |
|---|---|---|---|---|---|---|---|
| D01 | aspect ratio two ways, $\Lambda\equiv s^2/A=s/\bar c$ (14.1) with $\bar c=A/s$ (14.2), then level flight: $L=W$ ⇒ $C_L=\dfrac{W}{\tfrac12\rho U^2A}$ and $U_{\min}=\sqrt{\dfrac{2W}{\rho A\,C_{L,\max}}}$ (d-a, d-b; the second part ours) | C01 | ★ | 5 | substitution, $C_L\equiv\dfrac{L}{\tfrac12\rho U^2A}$ (4.108), square root | A is the planform area (T13); "level, steady" is the assumption that makes lift equal weight; $C_L$ falls as $1/U^2$ | notebook |
| D02 | the starting vortex carries −Γ: $\Gamma_{ABCD}=0=\Gamma_{\text{foil}}+\Gamma_{\text{starting}}$ from Kelvin's theorem $D\Gamma/Dt=0$ (5.8) (d-d) | C03 | ★ | 5 | Kelvin's theorem (ch05 recap), additivity of circulation over sub-circuits (Stokes, ch02 recap) | the circuit, not its inside, must stay in inviscid fluid (T15); it is a material circuit, so it must be large enough that the wake never crosses it; clockwise against counterclockwise (T1) | notebook · `starting_vortex` |
| D03 | the raised circle maps to a circular arc: (14.4) $x=(R+b^2/R)\cos\theta,\ y=(R-b^2/R)\sin\theta$ → (14.5) $x^2\sin^2\theta-y^2\cos^2\theta=4b^2\sin^2\theta\cos^2\theta$ → cosine rule → (14.6) $2b\tan\beta\,\sin\theta=R-b^2/R=y/\sin\theta$ → $x^2+(y+2b\cot2\beta)^2=(2b\csc2\beta)^2$, height $2b\tan\beta$, hence camber $\simeq\tfrac12\beta c$ (14.8) (d-f, d-g, d-h, d-j) | C04 | ★★★ | 15 | Euler's formula P45, difference of squares, cosine rule P340, double-angle identities P341, completing the square P272, small angle P100 | dividing by sin θ and cos θ (limits afterwards); why the angle at O is 90° − θ; only the arc above the chord is covered, and twice; β is the circle's angle, the end tangent is 2β; the camber rule is first order in β | notebook · `zhukhovsky_airfoil_builder` |
| D04 | the shifted circle gives a symmetric airfoil with $c\simeq4b$ and $t_{\max}/c\simeq\tfrac{3\sqrt3}{4}e\approx1.3\,e$ at the quarter chord (14.8), from $a\equiv b(1+e)$ (14.7) (d-i; book never writes it out — "can be proved") | C04 | ★★★ | 11 | binomial approximation with a complex small quantity P342, real and imaginary parts, maximum by a zero derivative, first order in e | the centre is at $-be$, not $-e$; $y\simeq2be\sin\theta(1-\cos\theta)$ has its maximum at θ = 2π/3, not π/2; thickness is twice y; the chord correction is second order; at e = 0.1 the rule is about 9 % high (slip #4) | notebook · `zhukhovsky_airfoil_builder` |
| D05 | the trailing edge is a cusp: $\dfrac{dz}{d\zeta}=1-\dfrac{b^2}{\zeta^2}=0$ at $\zeta=b$, so $z-2b\simeq(\zeta-b)^2/b$ and the angle π of the circle becomes 2π; a finite edge angle $2\pi(1-1/n)$ needs a stagnation point, speed $\propto r^{n-1}$ (d-k; book never writes it out) | C04 | ★★ | 7 | Taylor expansion about a critical point (ch06 D19 recap), squaring doubles angles P153, corner flows (ch06 C09 recap) | the other critical point $\zeta=-b$ must lie inside the circle; "cusp" is about the body's shape, "stagnation point" about the flow — the Kutta flow at a cusp has a finite, non-zero speed | notebook |
| D06 | flow past an off-centre circle at incidence: $w=U\!\left[(\zeta-\mu)e^{-i\alpha}+\dfrac{a^2e^{i\alpha}}{\zeta-\mu}\right]+\dfrac{i\Gamma}{2\pi}\ln\dfrac{\zeta-\mu}{a}$ from the centred, aligned flow of ch06 (6.68) (d-n; ours; book never writes it out) | C05 | ★★ | 7 | shifting and rotating a complex potential, complex logarithm P155, the circle as a streamline (ch06 recap) | rotate by $e^{-i\alpha}$, not $e^{+i\alpha}$ (check the far field $Ue^{-i\alpha}$); the logarithm's argument is rotated too, which only adds a constant; Γ is clockwise, so the loop circulation is −Γ | notebook |
| D07 | velocity in the airfoil plane, $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\,\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9), from the chain rule and $z=\zeta+b^2/\zeta$ (14.3) (d-l) | C05 | ★ | 4 | chain rule P49, derivative of the map | it is $d\zeta/dz=1/(dz/d\zeta)$ that enters; the two velocities agree only at infinity (T6); $dw/dz=u-iv$ (T2); infinite at $\zeta=\pm b$ unless $dw/d\zeta$ vanishes there | notebook · `kutta_condition` |
| D08 | the circulation is the same in both planes: $(dw/dz)\,dz=(dw/d\zeta)\,d\zeta=dw$, so $\oint$ round corresponding curves is identical (d-m; book never writes it out — "it can be shown") | C05 | ★★ | 5 | circulation as the real part of $\oint(dw/dz)\,dz$ (ch06 recap), closed-loop integral of a multi-valued function P135 | the curves must correspond one to one outside the body; the velocity is **not** the same, only its loop integral; this is what licenses $L=\rho U\Gamma$ (6.62) in the z-plane | notebook |
| D09 | Kutta circulation and lift: $u_\theta=-2U\sin\theta-\dfrac{\Gamma}{2\pi a}$ (14.10) with $u_\theta=0$ at $\theta=-(\alpha+\beta)$ ⇒ $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11) ⇒ with $L=\rho U\Gamma$ (6.62) and $c\simeq4a$, $C_L\simeq2\pi(\alpha+\beta)$ (14.12); exact form $8\pi(a/c)\sin(\alpha+\beta)$ (ours) (d-o, d-q) | C06 | ★★ | 10 | geometry of the offset circle, $\sin(-x)=-\sin x$, small angle P100, section coefficient P336 | why B is at $-(\alpha+\beta)$: β below the horizontal through Q, and the reference diameter is tilted up by α (T4); the radius is a, not b; three separate approximations ($\sin x\simeq x$, $c\simeq4a$, $e\ll1$) — thickness raises the slope above 2π (slip #8); radians (T12) | notebook · `kutta_condition` |
| D10 | the normal-force coefficient is the area of the pressure loop: $F_y=-\int_0^cp_u\,dx+\int_0^cp_l\,dx$ ⇒ $C_y\equiv\dfrac{F_y}{\tfrac12\rho U^2c}=\oint C_p\,d\!\left(\dfrac xc\right)$ (d-c; book never writes it out — an exercise) | C07 | ★★ | 7 | pressure force $-p\,\mathbf n\,ds$, $\mathbf n\,ds=(dy,-dx)$ on a contour P92, closed line integral, P337 | the sign of the normal on the lower surface; subtracting $p_\infty$ changes nothing on a closed contour; x along the chord, so this is the force normal to the **chord** — the lift differs by the tilt α (T5); shear neglected | notebook · `airfoil_pressure_lift` |
| D11 | strength of the shed filament, $d\Gamma_{\text{trail}}=-\dfrac{d\Gamma}{dy}\,dy$, from Helmholtz's theorem (d-t) | C08 | ★ | 4 | Helmholtz's theorem (R11), derivative as a difference | the sign: dΓ < 0 for dy > 0 on the starboard half, so the shed strength is positive there; the two half-sheets rotate in opposite senses | notebook · `downwash_induced_drag` |
| D12 | downwash $w(y_1)=\dfrac{1}{4\pi}\displaystyle\int_{-s/2}^{s/2}\dfrac{d\Gamma}{dy}\,\dfrac{dy}{(y_1-y)}$ (14.13), from the half-line result $dw(y_1)=\dfrac{-(d\Gamma/dy)\,dy}{4\pi(y-y_1)}$ (R12) (d-v) | C09 | ★★ | 7 | Biot–Savart segment formula (ch05 D13 recap), principal value P344, improper integral P144 | the minus sign absorbed into $(y_1-y)$; w is positive **downward** (T5); the integrand is singular at $y=y_1$ — a principal value (T10); the bound vortex induces nothing on its own line; slip #3 (which point the velocity is at) | notebook · `downwash_induced_drag` |
| D13 | effective angle $\alpha_e=\alpha-\dfrac wU$ (14.14), strip forces $dL\simeq\rho U\Gamma\,dy$, $dD_i\simeq\rho w\Gamma\,dy$, and $L=\displaystyle\int\rho U\Gamma\,dy$, $D_i=\displaystyle\int\rho w\Gamma\,dy$ (14.15) (d-w) | C10 | ★★ | 8 | relative-velocity triangle P350, small angle P100, resolving a force, $L=\rho U\Gamma$ (6.62) per unit span | $\varepsilon=\tan^{-1}(w/U)$, not tan (slip #2); $U_e\sin\varepsilon=w$ exactly, $U_e\cos\varepsilon=U$ exactly — only the names "lift" and "drag" are approximate; the drag is not friction | notebook · `downwash_induced_drag` |
| D14 | the series in the downwash integral: with $y=-(s/2)\cos\gamma$ and $\Gamma=\sum\Gamma_n\sin(n\gamma)$ (14.18), (14.13) becomes $w(y_1)=\dfrac{1}{2\pi s}\displaystyle\sum n\Gamma_n\int_0^\pi\dfrac{\cos(n\gamma)\,d\gamma}{\cos\gamma-\cos\gamma_1}$ — the first two members of (14.19) with the denominator corrected (d-y) | C11 | ★★ | 7 | change of variable P106, term-by-term differentiation of a series, half-range sine series P346 | $y_1-y=(s/2)(\cos\gamma-\cos\gamma_1)$: the page prints the opposite order (slip #1), which flips the sign of w; $dy=(s/2)\sin\gamma\,d\gamma$ cancels inside $d\Gamma$; limits γ = 0 to π | notebook |
| D15 | Glauert's integral $\displaystyle\int_0^\pi\dfrac{\cos n\gamma\,d\gamma}{\cos\gamma-\cos\gamma_1}=\pi\,\dfrac{\sin n\gamma_1}{\sin\gamma_1}$, hence $w(y_1)=\dfrac{1}{2s}\displaystyle\sum n\Gamma_n\dfrac{\sin(n\gamma_1)}{\sin\gamma_1}$ (14.19, last member) (d-z; ours; book never writes it out) | C11 | ★★★ | 11 | principal value P344, product-to-sum identity and three-term recursion P347, induction P306, l'Hôpital P343 for the tips | $I_0=0$ only as a principal value; the recursion needs **two** starting values ($I_0$, $I_1=\pi$); the closed form satisfies the same recursion — that is the proof; $0<\gamma_1<\pi$ | notebook · `lifting_line_solver` |
| D16 | the lifting-line equation: $C_L=K\left[\alpha-\dfrac wU+\beta\right]$ (14.16) and $\Gamma=\tfrac12UcC_L$ ⇒ $\Gamma(y)=\dfrac K2Uc(y)\left[\alpha-\dfrac{w(y)}{U}+\beta(y)\right]$ (14.17) ⇒ with (14.18), (14.19): $\dfrac K2Uc(\alpha+\beta)=\displaystyle\sum\left(1+\dfrac{nKc}{4s\sin\gamma}\right)\Gamma_n\sin(n\gamma)$ (14.20) (d-x, d-aa) | C12 | ★★ | 8 | substitution, collecting a series on $\sin n\gamma$, integral equation P345 | −β is the zero-lift angle here (T3); K per radian (T12); c, α, β, K may all depend on γ, so this is not an ordinary Fourier series and the coefficients cannot be found by projection; $\sin\gamma=0$ at the tips (T10) | notebook |
| D17 | lift from the first coefficient, $L=\dfrac{\pi s}{4}\rho U\Gamma_1$ (14.21), from (14.15) and (14.18) (d-ac; book gives the result only) | C13 | ★★ | 6 | change of variable P106, orthogonality on the half range P346 | $dy=(s/2)\sin\gamma\,d\gamma$ supplies the $\sin\gamma$ that selects n = 1; the factor π/2 of the half range (not π) | notebook |
| D18 | induced drag $D_i=\dfrac\pi8\rho\displaystyle\sum n\Gamma_n^2$ (14.22), from (14.15), (14.18), (14.19) (d-ad; book gives the result only) | C13 | ★★★ | 9 | double sum collapsed by orthogonality P348, P346 | the $\sin\gamma$ of dy cancels the $1/\sin\gamma$ of w — do not cancel it twice; two different summation indices are needed; the factor n comes from w, once; N against ∞ as the upper limit (slip #7) | notebook · `lifting_line_solver` |
| D19 | the elliptic optimum: fixed L fixes $\Gamma_1$; every other term adds $n\Gamma_n^2\ge0$ ⇒ $\Gamma=\Gamma_1\sqrt{1-(2y/s)^2}$ (14.23), $w=\Gamma_1/2s$ (14.24), $D_i=\dfrac\pi8\rho\Gamma_1^2=\dfrac{2L^2}{\pi\rho U^2s^2}$ (14.25), $C_{D_i}=\dfrac{C_L^2}{\pi\Lambda}$ (14.26) (d-ae, d-af) | C14 | ★★ | 9 | sum of non-negative terms, $\sin\gamma=\sqrt{1-\cos^2\gamma}$, substitution, section against wing coefficients P336 | "minimum" is at fixed lift **and span** with a planar wake; $\Gamma_1=4L/(\pi s\rho U)$ is squared; dividing by $\tfrac12\rho U^2A$ twice produces $C_L^2$ and one factor $A/s^2=1/\Lambda$ | notebook · `elliptic_loading_aspect_ratio` |
| D20 | elliptic planform and the finite-wing lift slope: $dL=\rho U\Gamma\,dy=C_L\tfrac12\rho U^2c\,dy$ (14.27) with uniform w ⇒ $c\propto\Gamma$; then (14.17) with (14.24) ⇒ $\dfrac{dC_L}{d\alpha}=\dfrac{K}{1+K/(\pi\Lambda)}$, i.e. $C_{L,3D}=\dfrac{C_{L,2D}}{1+2/\Lambda}$ for K = 2π (d-ag; the second part an exercise the book never writes out) | C14 | ★★ | 8 | proportionality, solving a linear equation for $\Gamma_1$, area of an ellipse $A=\pi c_0s/4$ | elliptic *loading* is not elliptic *planform* (T11); the planform gives the loading only if untwisted with one section; the factor is $K/(\pi\Lambda)$ in general and 2/Λ only for K = 2π | notebook |
| D21 | best glide and the U-shaped drag curve: from $C_D=C_{D0}+\dfrac{C_L^2}{\pi\Lambda e_{\text{span}}}$, $\left(\dfrac LD\right)_{\max}=\tfrac12\sqrt{\dfrac{\pi\Lambda e_{\text{span}}}{C_{D0}}}$ at $C_{D0}=C_{D_i}$; in level flight $D=\tfrac12\rho U^2AC_{D0}+\dfrac{2W^2}{\pi\rho U^2s^2e_{\text{span}}}$ (d-an; ours; book never writes it out) | C15 | ★★ | 8 | minimum of $a/x+bx$ P349, quotient rule P265, $L=W$ (D01) | minimise $C_D/C_L$, not $C_D$; the parabolic polar is a model with $C_{D0}$ and $e_{\text{span}}$ taken constant; "induced drag is the major part" holds only at and below this speed (T9) | notebook · `elliptic_loading_aspect_ratio` |
| D22 | sailing: apparent wind $\mathbf U_r=\mathbf U-\mathbf V$; drive $=L\sin\beta_a-D\cos\beta_a$, zero at $\tan\beta_a=D/L$; with the hull's own drag angle the closest course to the apparent wind is $\beta_a=\varepsilon_A+\varepsilon_H$ (d-am; ours; book resolves the force in words only) | C16 | ★★ | 8 | relative velocity P350, resolving a force along and across the course, `np.arctan2` P70 | the angle is to the **apparent** wind, which is further ahead and faster than the true wind; L is normal to $\mathbf U_r$, not to the sail; the lateral force does not vanish — the keel must carry it; quasi-steady, forces in the horizontal plane (T18) | notebook · `thrust_from_lift` |

22 rows: ★ 4 (D01, D02, D07, D11) · ★★ 14 (D05, D06, D08, D09, D10, D12, D13, D14, D16, D17, D19, D20, D21, D22) · ★★★ 4 (D03, D04,
D15, D18). They cover 29 analysis derivations: d-a, d-b, d-c, d-d, d-f, d-g, d-h, d-i, d-j, d-k, d-l, d-m, d-n, d-o, d-q, d-t,
d-v, d-w, d-x, d-y, d-z, d-aa, d-ac, d-ad, d-ae, d-af, d-ag, d-am, d-an. Rule (b) rows (B item, book never writes it out): D04, D05,
D06, D08, D10, D15, D17, D20 (second part), D21, D22. Every ★★★ row belongs to an A item with an explainer and is in that explainer's
Derivation tab (D03, D04 → E1; D15, D18 → E7).

## 4c. Derivations demoted to statements (first column is the A parent in bold, e.g. **C20** — never a bare ID or a D id: the parser reads a bare first-cell ID as an item and blanks its tier)
| A parent | Analysis §2b item | Result stated (Eq.) | Stated in (B item) | Why not written out |
|---|---|---|---|---|
| **C04** | d-e centred circle → slit of length 4b (★) | $z=be^{i\theta}+be^{-i\theta}=2b\cos\theta$ (unnumbered) | R06 | taught in ch06 (D20); one line with Euler's formula, shown as the first panel of the four-circle figure |
| **C06** | d-p finite speed at the trailing edge (★★★, ours) | $\left\lvert\dfrac{dw}{dz}\right\rvert_{B'}=U\,\dfrac ba\cos(\alpha+\beta)$ (ours) | N29 | a B item that is not in the book and that nothing later builds on; the value is given, a sympy `limit` cell and a table of speeds at shrinking distances verify it; the l'Hôpital idea is primer P343 |
| **C07** | d-r the force is perpendicular to the inclined stream, D = 0 (★★, ours) | $D-iL=\dfrac{i\rho}{2}\oint\left(\dfrac{dw}{dz}\right)^2dz$ (6.60) ⇒ $\lvert F\rvert=\rho U\Gamma$, normal to the stream | R10 | ch06 D18 did the residue calculation; with α ≠ 0 it is the same calculation rotated; two numerical routes (Blasius contour, surface pressure) agree in a cell |
| **C07** | d-s leading-edge suction on a flat plate (★★★, ours; an exercise) | $N=L\cos\alpha,\ S=L\sin\alpha=\pi\rho U^2c\sin^2\alpha$ (ours) | N12 | a B item; the control-volume calculation round an integrable singularity is long and feeds nothing later; the force balance $N\sin\alpha=S\cos\alpha$ is checked in a cell with `ch14.edge_flow_suction` |
| **C09** | d-u a semi-infinite straight filament induces half the velocity of an infinite line (★★) | $dw(y_1)=\dfrac{-(d\Gamma/dy)\,dy}{4\pi(y-y_1)}$ (unnumbered), from the segment formula $\dfrac{\Gamma}{4\pi d}(\cos\theta_a-\cos\theta_b)$ of (5.17) | R12 | taught in ch05 (D13); a recap with `BS.segment_speed` at angles 0 and 90°; it is the "we start from" line of D12 |
| **C09** | d-ah horseshoe vortex: centre downwash and its too-good induced drag (★★; an exercise) | $w(0)=\dfrac{\Gamma}{\pi s},\qquad C_{D_i}=\dfrac{C_L^2}{2\pi\Lambda}$ (ours) | N60 | a B item used as the block's tiny example: two half-lines at distance s/2 added by hand in three lines of arithmetic, then checked with `LL.horseshoe_velocity`; the drag figure is shown to be an artefact (T10) rather than derived |
| **C09** | d-ai far-wake downwash is 2w (★) | $w(x\to-\infty)=0,\quad w(0)=w,\quad w(x\to+\infty)=2w$ (Fig. 14.22) | N41 | one limit of the segment formula; shown as a computed curve |
| **C10** | d-aj induced drag = kinetic energy left per unit length of wake (★★★, ours) | $D_iU=\dfrac{d}{dt}(\text{wake kinetic energy})$; energy per length $=\tfrac12\rho\displaystyle\int\Gamma\,w_\infty\,dy$ with $w_\infty=2w$ (ours) | N38 | a B item the book states in words; Green's theorem across a vortex sheet is beyond what the block needs; the equality with (14.15) is verified numerically to round-off |
| **C12** | d-ab (14.20) → N linear equations (★★, ours formalisation) | $\mathsf A_{jn}\Gamma_n=r_j,\quad\mathsf A_{jn}=\left(1+\dfrac{nKc_j}{4s\sin\gamma_j}\right)\sin n\gamma_j$ (ours) | N52 | not a derivation but the definition of a linear system; written for N = 3 and assembled line by line in the from-scratch cell (§7) |
| **C12** | d-ak twist for an elliptic loading on a rectangular wing (★★, ours; an exercise) | $\alpha_t=\dfrac{Kc\,\alpha_r}{4s},\quad\Gamma_1=\tfrac12KUc\,\alpha_r$ (ours) | N63 | a B item; matching two coefficients in (14.17); the solver confirms it (no higher coefficients) |
| **C16** | d-al fish tail: thrust power equals sideways power (★) | $UL\sin\theta=VL\cos\theta$ because $V/U=\tan\theta$ (unnumbered) | N77 | the book shows every step in words; two lines of arithmetic done as the block's tiny example with easy numbers |

11 rows. With §4b: 29 + 11 = 40 = every derivation of analysis §2b (d-a … d-an), none missing.

## 5. Interactive explainers (5–10 + backup)
Nine explainers and one backup — not ten: every A item where manipulation beats a figure has one, and a tenth would be padding.
Fifteen of the sixteen A items are attached to one (C02–C16); C01 (vocabulary) is served by a slider figure of planforms. Every
window shows its equations with numbers, has an Explain tab ("Explanation & interpretation", numbered sections with the reader's
numbers and a "Reading the current setting" paragraph), a synced Code tab, a Derivation tab for its D rows, a 4–8-step walkthrough
and ≥ 3 check questions. **Colour code across the chapter:** lift = teal, drag (induced, profile) = rose, thrust / drive = purple
(accent), circulation and vorticity = amber, pressure and suction = orange, the stream and velocities = blue, references and ghosts
= muted. **Γ is labelled "Γ (clockwise)" in every readout** and the Explain tab's first section says so (trap T1). Angles are shown
in degrees and computed in radians, with both in the Explain arithmetic (T12). Reference explainers (`interactive-viz` §5): FDV =
`forced_damped_vibrations.html`, AFE = `angular_frequency_explorer_1.html`, APS = `amplitude_phase_second_order_II_3.html`. **JS
physics stays light:** closed forms in complex arithmetic (a local `cx` helper — `Viz.cx` is not in the library yet, see §8), one
bisection (E2), a lumped-vortex march of ≤ 300 steps (E3), a dense solve of at most 12 × 12 by a local LU (E7); no tables. The
constant in the thickness rule is computed in the page as $3\sqrt3/4$, never typed. **Phones (360 × 640):** each explainer names the
view it hides (`hidePortrait`) and repeats that view's key number in a visible title or readout.

### E1 · zhukhovsky_airfoil_builder
- A: C04 (also shows R05 the map, R06 the slit, N24 the arc, N25 the symmetric section, N26 the shape rules, N27 the cusp, N07 section geometry) · **Confusion removed:** "where do the camber and the thickness of an airfoil come from in a map that turns a circle into a flat plate?" — from *where the circle sits*: raise its centre and the plate bends; move it back and the nose fattens; keep it through ζ = b and the tail stays sharp.
- **Why interactive:** dragging the centre μ of the circle with the pointer and watching the section reshape continuously is the idea itself; four static panels (the book's four figures) show the corners of the square, not that β and e are independent knobs, nor how far the rules of (14.8) drift as the section gets fat.
- **Stage:** (1) the ζ-plane: the circle through ζ = b with its draggable centre, the reference circle $\lvert\zeta\rvert=b$ as a ghost, the two critical points ±b marked; (2) the z-plane: the mapped contour with chord, camber line and maximum-thickness markers, a marked point travelling round both curves (transport); (3) (hidePortrait) "rule against measured" bars for chord, camber and thickness.
- **Controls:** centre of the circle (drag in view 1, or sliders e and β) · thickness parameter e · camber angle β · marked-point angle θ (transport) · zoom on the trailing edge (toggle).
- **Equations:** $z=\zeta+\dfrac{b^2}{\zeta}$ (14.3); $x=(R+b^2/R)\cos\theta,\ y=(R-b^2/R)\sin\theta$ (14.4); $a\equiv b(1+e)$ (14.7); $c\simeq4b$, camber $\simeq\tfrac12\beta c$, $t_{\max}/c\simeq1.3\,e$ (14.8); the arc $x^2+(y+2b\cot2\beta)^2=(2b\csc2\beta)^2$ (unnumbered); $\dfrac{dz}{d\zeta}=1-\dfrac{b^2}{\zeta^2}$.
- **Mirrors** `AF.zhukhovsky_circle`, `AF.zhukhovsky_airfoil`, `AF.zhukhovsky_shape(exact=True / False)`, `AF.circular_arc`, `AF.section_geometry`, `CONF.joukowski`.
- derivations: D03, D04 · depth features: explain (the circle from your e and β; one point mapped by hand; chord, camber and thickness measured against the three rules, with the percentage gap; the reading), code, **linked views**, **transport** (the point going round), **presets** (flat plate, circular arc, symmetric, cambered — the book's four cases — and "too fat: the rules fail"), **status** ("thin and slightly cambered: rules within 2 %" / "e = 0.25: thickness rule 20 % high"), **inspector** (click a point of the circle: its image, with the arithmetic of (14.4)) · follows: APS (linked windows, numbered Explain with live numbers) and ch06 `conformal_joukowski` (two planes, same marked point) · phones: the bars are hidden; the three measured numbers sit in the z-plane title · **aha:** camber is how high the circle sits and thickness is how far back it sits — the sharp tail is the one point of the circle that the map folds.

### E2 · kutta_condition
- A: C02, C05, C06 (also shows R04 $L=\rho U\Gamma$, R09 the surface speed on the circle, N28 the ζ-plane flow, N29 the finite edge speed, N30 the Kutta circulation, N34 the flat plate, N14 the three stages) · **Confusion removed:** "ideal flow round an airfoil has a solution for every circulation — so what decides the lift?" — only one circulation puts the rear stagnation point on the sharp edge; every other one asks the fluid to turn an infinitely sharp corner at infinite speed.
- **Why interactive:** sliding Γ and watching the rear stagnation point crawl along the surface toward the edge while the edge-speed readout falls from "∞" to a finite number, then rises again past the Kutta value, is the hypothesis as an experiment; a static figure can show three values of Γ but not that exactly one is special, nor how that one moves with α and β.
- **Stage:** (1) the z-plane: the airfoil with streamlines and tracer particles, both stagnation points marked, the edge region shaded when the speed there is unbounded; (2) the ζ-plane: the circle flow with the same particles, the point B and the angle $\alpha+\beta$ drawn at the centre Q; (3) (hidePortrait) speed near the trailing edge against Γ/Γ_K with the cursor, and the lift bar $\rho U\Gamma$.
- **Controls:** circulation Γ/Γ_K (0 to 1.5) with a "snap to Kutta" button · angle of attack α · camber angle β · thickness e (optional) · play (particles).
- **Equations:** Kutta condition $\left.\dfrac{dw}{d\zeta}\right\rvert_{\zeta=b}=0$; $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\,\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9); $u_\theta=-2U\sin\theta-\dfrac{\Gamma}{2\pi a}$ (14.10); $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11); $C_L\simeq2\pi(\alpha+\beta)$ (14.12); $L=\rho U\Gamma$ (6.62); the ζ-plane potential of N28.
- **Mirrors** `AF.zhukhovsky_flow`, `AF.kutta_circulation`, `AF.kutta_circulation_root`, `AF.stagnation_points` (§8 note 2), `AF.trailing_edge_speed`, `AF.zhukhovsky_lift_coefficient`, `ch06.lift_per_span`.
- derivations: D07, D09 · depth features: explain (the circle from your e, β; where B sits, $\theta_B=-(\alpha+\beta)$, in degrees and radians; Γ_K from (14.11) and from the bisection, side by side; the edge speed; $C_L$ three ways — 2π(α + β), exact, and 2Γ/(Uc); the reading), code, **linked views**, **transport** (particles; a "start-up" ramp of Γ from 0 to Γ_K), **presets** (no circulation, half way, Kutta, too much, flat plate, zero-lift angle α = −β), **status** ("rear stagnation point on the upper surface, 12 % of chord ahead of the edge — edge speed unbounded" / "Kutta: flow leaves the edge smoothly at 0.93 U"), **terms** (lift $\rho U\Gamma$ against the Kutta lift as two bars) · follows: FDV (system + graph on one slider, the Explain panel) and ch06 `conformal_joukowski` · phones: the edge-speed graph is hidden; the edge speed and Γ/Γ_K are in the z-plane title · **aha:** the wing does not choose its lift — its sharp trailing edge does.

### E3 · starting_vortex
- A: C03 (also shows N14 the three stages, N15 the mechanism, N16 the punch line, R11 the closed loop, the Kelvin circuit of D02) · **Confusion removed:** "Kelvin says circulation cannot be created in an inviscid fluid — so where does the wing's circulation come from?" — it is not created: the wing sheds an equal and opposite vortex and flies away from it.
- **Why interactive:** starting, stopping and changing the incidence of the foil and seeing a new vortex of the right sign roll off the edge each time, while the amber total stays flat at zero, is something only time can show; the reader can also move the circuit and see its circulation jump by exactly the strength of whatever it now encloses.
- **Stage:** (1) the foil moving to the left with its wake of shed vortices (amber, sized and signed by strength), the bound circulation as a ring round the foil, a resizable material circuit; (2) circulation against distance travelled in chords: bound (teal), shed total (rose), their sum (amber, flat at zero), the steady Kutta value as a ghost; (3) (hidePortrait) lift against time with the steady value as a ghost.
- **Controls:** angle of attack α · event (chips: impulsive start, sudden stop, step in α, flap down) · circuit size · transport t (play / step / scrub) · wake model (toggle: rigid wake / free wake).
- **Equations:** $\Gamma_{ABCD}=0=\Gamma_{\text{foil}}+\Gamma_{\text{starting}}$; Kelvin's theorem $D\Gamma/Dt=0$ (5.8); the steady value $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11) (flat plate: $\pi Uc\sin\alpha$); $L=\rho U\Gamma$ (6.62).
- **Mirrors** `ch14.impulsive_start_lumped_vortex` (bound vortex at c/4, no-through-flow at 3c/4, one wake vortex per step; ours), `AF.kutta_circulation`, `BS.point_vortex_velocity`.
- derivations: D02 · depth features: explain (the steady circulation for your α; the bound, shed and total circulation at the current time, added in a boxed line; how many chords until 90 %; what the circuit encloses now; the reading), code, **linked views**, **transport** with an **end-of-run summary** card, **presets** (the four events), **status** ("after 3.0 chords: bound 0.81 Γ_K, shed −0.81 Γ_K, total 0.000"), **notes** (the model is a lumped-vortex sketch, not a simulation of the boundary layer) · follows: AFE (system animation + graph on one clock, end-of-run summary) and ch13 `geostrophic_adjustment` (conserved total flat in amber) · phones: the lift graph is hidden; lift/steady lift is a readout · **aha:** every change of lift leaves a vortex behind; the wing and its starting vortex are two ends of one loop.

### E4 · airfoil_pressure_lift
- A: C07 (also shows R03 the pressure coefficient, N09 the loop area, N10 circulation as the surface integral of speed, N12 leading-edge suction, N13 the adverse gradient and the qualitative separation point, R10 the two force routes) · **Confusion removed:** "the wing is pushed up from below" — mostly it is pulled up from above: the suction peak near the nose carries most of the loop, and the steep recovery behind it is where the boundary layer will give up.
- **Why interactive:** turning α and seeing the suction peak sharpen, the stagnation point slide under the nose and the loop area grow in step with the $C_L$ readout connects three pictures of one flow; clicking a surface point traces Bernoulli's arithmetic there; no single plot shows that the area *is* the lift.
- **Stage:** (1) the section with pressure arrows normal to the surface (outward = suction, orange), the stagnation point marked; (2) $-C_p$ against x/c for both surfaces with the enclosed area shaded and a cursor; (3) (hidePortrait) bars: loop area × cos α, $2\Gamma/(Uc)$ and 2π(α + β).
- **Controls:** angle of attack α · camber angle β · thickness e · surface cursor (transport along the contour) · show the laminar separation estimate (toggle, labelled qualitative).
- **Equations:** $C_p=\dfrac{p-p_\infty}{\tfrac12\rho U^2}=1-\dfrac{\lvert dw/dz\rvert^2}{U^2}$; $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\,\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9); $C_y=\oint C_p\,d\!\left(\dfrac xc\right)$; $L=\rho U\Gamma$ (6.62); $C_L\simeq2\pi(\alpha+\beta)$ (14.12).
- **Mirrors** `AF.zhukhovsky_surface`, `AF.normal_force_coefficient`, `AF.lift_from_cp`, `AF.surface_circulation`, `AF.zhukhovsky_force`, `ch14.suction_side_separation`.
- derivations: D10 · depth features: explain (the speed and $C_p$ at the cursor from (14.9), step by step; the suction peak and where it is; the loop area by the trapezoid rule against $2\Gamma/(Uc)$; the reading), code, **linked views**, **transport** (the surface cursor), **inspector** (click a point: ζ, the two factors of (14.9), the speed, $C_p$), **terms** (upper-surface suction + lower-surface pressure = normal force), **presets** (zero lift, cruise, high α with a sharp peak, thick symmetric, flat plate with its infinite peak), **status** ("suction peak −2.4 at 1 % chord; pressure recovers over 99 % of the upper surface") · follows: APS (linked windows with a crosshair readout) · phones: the bars are hidden; the three lift numbers are in the graph title · **aha:** the lift is the area between the two pressure curves, and most of it is suction near the nose.

### E5 · lift_curve_stall
- A: C06 (also shows N31 the zero-lift angle, N32 theory against measurement, N13 stall, N70 the generic measured curve, N72 the Reynolds-number regimes, N73 flaps and the stall speed, N59 the finite-wing slope as a ghost) · **Confusion removed:** "more angle always means more lift" and "camber changes the slope" — camber slides the whole straight line sideways without changing its slope of 2π per radian, and the line ends at stall, where the Kutta condition fails.
- **Why interactive:** dragging α along the curve while the small section view shows the flow attached, then separating, with the camber slider shifting the curve and the regime chips moving the stall point, lets the reader find the three facts (slope, shift, limit) by experiment.
- **Stage:** (1) $C_L$ against α in degrees: the line 2π(α + β) (teal), the exact Zhukhovsky curve (dashed), the generic measured curve with its stall (rose, labelled "model"), the cursor; (2) the section at the current α with attached or separated streamlines sketched and the lift arrow; (3) (hidePortrait) the stall speed of our example aircraft as a bar, with and without flap.
- **Controls:** angle of attack α · camber (as camber/chord, giving β) · Reynolds-number regime (chips: laminar separation first / transition first) · flap (toggle) · thickness e (optional).
- **Equations:** $C_L\simeq2\pi(\alpha+\beta)$ (14.12); $\alpha_{L=0}=-\beta$, $\beta=2\,\text{camber}/c$; exact $8\pi(a/c)\sin(\alpha+\beta)$ (ours); $C_L=K\left[\alpha-\dfrac wU+\beta\right]$ (14.16) with w = 0; $U_{\min}=\sqrt{\dfrac{2W}{\rho A\,C_{L,\max}}}$ (ours).
- **Mirrors** `AF.zhukhovsky_lift_coefficient(exact=)`, `AF.zero_lift_angle`, `ch14.generic_lift_curve` (qualitative, constants passed by keyword), `ch14.stall_regime`, `ch14.stall_speed`, `LL.finite_wing_lift_slope`.
- derivations: none (the derivation of the lift law is stepped through in `kutta_condition`; one walkthrough step quotes its result with the equation) · depth features: explain (β from your camber; $C_L$ at your α in radians and degrees; the gap between line, exact and model; the stall speed; the reading), code, **linked views**, **presets** (symmetric, cambered, flap down, at the zero-lift angle, just past stall), **status** ("linear range: Kutta condition holds" / "stalled: lift below the line by 35 % — model"), **notes** (what is theory, what is a drawn model, and why the measured numbers are not shown) · follows: `overfitting_curves.html` (a clear two-slider figure with a marker for the optimum and a verdict line) with an FDV-style Explain panel · phones: the stall-speed bar is hidden; the speed is a readout · **aha:** camber buys lift at zero angle, the angle buys 2π per radian, and stall is where the bill arrives.

### E6 · downwash_induced_drag
- A: C08, C09, C10 (also shows N39 the assumptions, N40 the bound vortex, R12 the half-line formula, N60 the horseshoe, N41 the far-wake factor 2, N42–N44 the tilted stream and strip forces, N65 the upwash outboard, N66 ground effect) · **Confusion removed:** "drag needs friction" — a finite wing in ideal flow has drag, because its own trailing vortices make it fly through sinking air, and lift at right angles to a sinking stream leans backward.
- **Why interactive:** building the wing from one horseshoe, then 3, 7, 21 nested horseshoes and watching the staircase Γ(y) and the downwash w(y) settle onto the smooth curves shows what the integral (14.13) *is*; moving the station cursor along the span turns the local force triangle; a static figure shows one distribution at one station.
- **Stage:** (1) a 3-D or oblique view of the lifting line with its trailing filaments (thickness and colour by strength), the station cursor; (2) Γ(y) (amber) and w(y) (blue) along the span, staircase and limit; (3) (hidePortrait) the force triangle at the station: U, w, the tilted stream, $dL_e$ normal to it, split into lift (teal) and induced drag (rose).
- **Controls:** number of horseshoes N (1, 3, 7, 21, ∞) · loading shape (chips: constant, elliptic, parabolic, rectangular-wing solution) · station $y_1$ (transport along the span) · height above ground h/s (optional) · field point behind the wing x/s (optional).
- **Equations:** $d\Gamma_{\text{trail}}=-\dfrac{d\Gamma}{dy}\,dy$; $w(y_1)=\dfrac{1}{4\pi}\displaystyle\int_{-s/2}^{s/2}\dfrac{d\Gamma}{dy}\dfrac{dy}{(y_1-y)}$ (14.13); $\alpha_e=\alpha-\dfrac wU$ (14.14); $L=\displaystyle\int\rho U\Gamma\,dy$, $D_i=\displaystyle\int\rho w\Gamma\,dy$ (14.15); horseshoe $w(0)=\dfrac{\Gamma}{\pi s}$ (ours).
- **Mirrors** `LL.horseshoe_velocity`, `LL.staircase_downwash` (§8 note 4), `LL.downwash`, `LL.trailing_sheet_strength`, `LL.strip_forces`, `LL.lift_and_induced_drag`, `LL.ground_effect_factor`, `BS.segment_speed`.
- derivations: D11, D12, D13 · depth features: explain (each filament's contribution at your station, the nearest ones written out and the sum; ε in degrees; the strip's lift and induced drag; the whole-wing totals; the reading), code, **3-D view** (falls back to the oblique 2-D drawing without WebGL), **linked views**, **transport** (the station), **terms** (downwash at the station as bars per filament, port and starboard, summing to w), **presets** (one horseshoe, elliptic, in ground effect, outboard of the tip = upwash), **status** ("downwash 0.031 U, stream tilted 1.8° down, induced drag 3.1 % of lift") · follows: AFE (system + graph on one cursor, live working) and ch05 `biot_savart_filament` ("build the sum") · phones: the force triangle is hidden; ε and $dD_i/dL$ are in the span graph's title · **aha:** the wing's own wake makes the air it meets sink, and lift perpendicular to sinking air leans back.

### E7 · lifting_line_solver
- A: C11, C12, C13 (also shows N49 the span angle, N50 the sine shapes, N51 Glauert's integral, N52 the linear system, N54 span efficiency, N63 twist, N58 the elliptic planform, N68 limits) · **Confusion removed:** "the best wing is an elliptic *wing*" and "the solver is a black box" — the unknown is a handful of sine amplitudes, the lift uses only the first, and every other one is pure induced drag; taper or twist can switch the others off on a wing that is not elliptic at all.
- **Why interactive:** changing taper, twist and aspect ratio and watching the coefficient bars, the loading Γ(y) against its elliptic ghost, the downwash and the efficiency respond together shows the method's inputs and outputs at once; raising N from 1 to 12 shows the truncation converge.
- **Stage:** (1) the planform (top view) with the collocation stations marked and the twist drawn as a colour band; (2) Γ(y) (amber) with the elliptic loading of the same lift as a ghost, and w(y) (blue) with its uniform ghost; (3) (hidePortrait) coefficient bars $n\Gamma_n^2$: the first (teal, "lift") and the rest (rose, "drag only").
- **Controls:** planform (chips: rectangular, tapered, elliptic) with taper ratio · aspect ratio Λ · twist at the tip (washout) · number of terms N · section slope K (optional, per radian).
- **Equations:** $y=-(s/2)\cos\gamma$; $\Gamma=\sum\Gamma_n\sin(n\gamma)$ (14.18); $w(y_1)=\dfrac{1}{2s}\sum n\Gamma_n\dfrac{\sin(n\gamma_1)}{\sin\gamma_1}$ (14.19, last member); $\dfrac K2Uc(\alpha+\beta)=\sum\left(1+\dfrac{nKc}{4s\sin\gamma}\right)\Gamma_n\sin(n\gamma)$ (14.20); $L=\dfrac{\pi s}{4}\rho U\Gamma_1$ (14.21); $D_i=\dfrac\pi8\rho\sum n\Gamma_n^2$ (14.22); $e_{\text{span}}=\Gamma_1^2/\sum n\Gamma_n^2$ (ours).
- **Mirrors** `LL.solve_lifting_line`, `LL.lifting_line_matrix` (§8 note 5), `LL.circulation_series`, `LL.downwash_series`, `LL.glauert_integral`, `LL.lift_from_series`, `LL.induced_drag_from_series`, `LL.span_efficiency`, `LL.planform`.
- derivations: D15, D18 · depth features: explain (the stations; one row of the matrix written with your numbers; the solved coefficients; $C_L$ from $\Gamma_1$ alone; $C_{D_i}$ term by term; $e_{\text{span}}$; the reading), code, **linked views**, **terms** (the bars of $n\Gamma_n^2$ summing to the induced drag; click one to see its shape sin nγ on the span), **presets** (elliptic planform = one bar; rectangle Λ = 6; taper 0.4, nearly elliptic; rectangle with elliptic twist = one bar again; N = 1 against N = 12), **status** ("span efficiency 0.954: 4.8 % more induced drag than the elliptic loading"), **inspector** (click a station: left and right sides of (14.20) there) · follows: `fid_formula_lab.html` (term-by-term bars with a "make it perfect" preset) with an APS-style numbered Explain · phones: the planform is hidden; taper, Λ and $e_{\text{span}}$ are in the loading graph's title · **aha:** lift lives in the first sine; everything else the wing does along its span is paid for in drag.

### E8 · elliptic_loading_aspect_ratio
- A: C14, C15 (also shows N55–N57 the elliptic results, N59 the finite-wing lift slope, N74 best glide, N75 the induced share, N73 the stall speed, C01 the aspect ratio) · **Confusion removed:** "a bigger wing has more drag" and "induced drag is always the largest part" — at fixed weight, induced drag depends on the *span* and falls with speed, profile drag rises with speed, and the two are equal exactly at the best glide.
- **Why interactive:** stretching the span at fixed area reshapes the planform, flattens the polar and slides the minimum of the U-shaped drag curve at once; the speed cursor swaps which bar is larger — the reader sees trap T9 resolve itself.
- **Stage:** (1) the planform of the current Λ beside real-world silhouettes scaled to it (a glider, an airliner, a delta — drawn by our code); (2) drag against speed in level flight: induced (rose), profile (blue), total, the best-glide speed marked, the stall speed as a wall; (3) (hidePortrait) the polar $C_L$ against $C_D$ with the tangent from the origin and the cursor.
- **Controls:** aspect ratio Λ (at fixed area, or at fixed span — toggle) · speed U (cursor) · zero-lift drag $C_{D0}$ · span efficiency $e_{\text{span}}$ · system (modes: glider, airliner, albatross, wind-turbine blade section — illustrative inputs).
- **Equations:** $\Lambda\equiv s^2/A$ (14.1); $w=\Gamma_1/2s$ (14.24); $D_i=\dfrac{2L^2}{\pi\rho U^2s^2}$ (14.25); $C_{D_i}=\dfrac{C_L^2}{\pi\Lambda}$ (14.26); $\dfrac{dC_L}{d\alpha}=\dfrac{K}{1+K/(\pi\Lambda)}$ (ours); $C_D=C_{D0}+\dfrac{C_L^2}{\pi\Lambda e_{\text{span}}}$ and $\left(\dfrac LD\right)_{\max}=\tfrac12\sqrt{\dfrac{\pi\Lambda e_{\text{span}}}{C_{D0}}}$ (ours).
- **Mirrors** `LL.induced_drag_coefficient`, `LL.elliptic_induced_drag`, `LL.finite_wing_lift_slope`, `ch14.drag_polar`, `ch14.drag_breakdown`, `ch14.best_glide`, `ch14.drag_vs_speed`, `ch14.level_flight`, `ch14.stall_speed`.
- derivations: D19, D21 · depth features: explain ($C_L$ needed at your speed from $L=W$; $C_{D_i}$ from (14.26); the two drags in newtons and their shares; the best-glide speed and ratio; the lift slope at your Λ; the reading), code, **linked views**, **modes** (four flyers), **terms** (profile + induced = total, as bars with percentages), **presets** (best glide, cruise, slow climb, Λ doubled, infinite span), **status** ("induced drag is 50 % of the total: this is the best-glide speed" / "cruise: induced 31 %"), **notes** with a small table of typical aspect ratios, the current row highlighted · follows: AFE (modes, presets, regime-dependent notes with a highlighted table) · phones: the polar is hidden; L/D and the induced share are in the drag graph's title · **aha:** the price of lift is set by the span you spread it over — and you pay most of it when you fly slowly.

### E9 · thrust_from_lift
- A: C16 (also shows N77 the fish tail, N79 the bird's two strokes, N81 the keel, N82 tacking and velocity made good, N08 lift and drag resolved) · **Confusion removed:** "a sail is pushed by the wind, so a boat cannot sail upwind" — a sail is a wing: its lift is perpendicular to the wind the boat *feels*, and that has a forward part as long as the apparent wind is further off the bow than the sail's drag angle.
- **Why interactive:** turning the course toward the wind and watching the apparent wind swing forward, the lift vector rotate and the drive bar shrink to zero at $\tan\beta_a=D/L$ is the whole argument; switching mode to the fish tail or the bird's wing shows the same triangle with the roles of "stream" and "motion" exchanged — which a single diagram cannot do.
- **Stage:** (1) plan view: the boat on its course, true wind (blue), boat velocity, apparent wind, sail lift (teal) and drag (rose), their resultant split into drive (purple) and side force, the keel's opposing force; in fish / bird mode the foil with U, V and the relative velocity; (2) drive and velocity made good against course angle with the cursor and the no-go zone shaded; (3) (hidePortrait) power bars: input against useful output (fish, bird) or drive against side force (sail).
- **Controls:** system (modes: sail, fish tail, bird wing) · course angle to the true wind θ (sail) or stroke speed ratio V/U (fish, bird) · lift-to-drag ratio of the sail or foil · boat speed over wind speed V/U · lift-to-drag ratio of the hull and keel (optional).
- **Equations:** $\mathbf U_r=\mathbf U-\mathbf V$; drive $=L\sin\beta_a-D\cos\beta_a$, zero at $\tan\beta_a=D/L$ (ours); closest course $\beta_a=\varepsilon_A+\varepsilon_H$ (ours); fish: thrust $=L\sin\theta$ at rate $UL\sin\theta$, side force $L\cos\theta$ at rate $VL\cos\theta$, $V/U=\tan\theta$; velocity made good $=V\cos\theta$ (ours).
- **Mirrors** `ch14.apparent_wind`, `ch14.sail_forces`, `ch14.sailboat_course_limit`, `ch14.sailboat_equilibrium`, `ch14.oscillating_foil`, `AF.resolve_lift_drag`.
- derivations: D22 · depth features: explain (the apparent wind from your course and speed ratio, by components; its angle off the bow; lift and drag resolved along and across the course; the limiting course from the two drag angles; in fish mode the two powers and their ratio; the reading), code, **linked views**, **modes**, **terms** (drive = lift part − drag part, as signed bars), **presets** (beam reach, close-hauled, in the no-go zone, a perfect foil D = 0, a poor sail), **status** ("close-hauled 42° off the true wind: drive positive, making 0.52 U to windward" / "in irons: drive negative"), **transport** in fish mode (the tail's stroke, thrust arrow pulsing) · follows: AFE (the same idea in several systems, presets, live working) and ch13 `ekman_force_balance` (a force polygon that closes) · phones: the bars are hidden; drive/L and the limiting angle are in the plan view's title · **aha:** lift is perpendicular to the wind you feel, not the wind that blows — move across the stream and it leans forward.

### B1 · wake_rollup (backup)
- A: C08 (also shows N36 the trailing sheet, N37 the tip vortices, spacing and descent speed, N38 wake energy, N67 core size) · **Confusion removed:** "the wing sheds two tip vortices" — it sheds a whole sheet, strongest near the tips, and the sheet winds itself into two vortices behind the wing.
- **Why interactive:** the roll-up is a process in time, and scrubbing it while the circulation and impulse readouts stay flat separates what changes (the shape) from what cannot (the strength and the spacing of the centroids).
- Stage: a rear view of the sheet as point vortices coloured by sign rolling up in time; the centroid of each half (fixed spacing πs/4 for the elliptic loading) and the pair's descent; circulation and impulse readouts flat in amber. Controls: loading shape (elliptic, rectangular-wing solution, with flap), number of vortices, core size, time (transport). Equations $d\Gamma_{\text{trail}}=-\dfrac{d\Gamma}{dy}\,dy$, $\Gamma_{\text{tip}}=\Gamma_0$, spacing $b'=\dfrac{1}{\Gamma_0}\displaystyle\int\Gamma\,dy$ and descent speed $\Gamma_0/2\pi b'$ (ours). Mirrors `ch14.wake_rollup`, `LL.rolled_up_pair`, `LL.trailing_sheet_strength`, `BS.point_vortex_evolve`, `BS.point_vortex_invariants`. Derivations D11 · depth: linked views, transport, presets, status ("qualitative past the first turn"), terms (invariants). Follows ch05 `vortex_sheet_rollup` and AFE. Built only if an explainer above fails review.

## 6. Python animations and interactive figures (A ID → what, why, player/figure kind)
Animations (5, `animate` + `show_animation`; ≤ 90 frames, FAST → 40; dpi 80). Cost is wall time in the notebook.
| # | A | What moves | Why | Player | Cost |
|---|---|---|---|---|---|
| A1 | **C02** (N14) | the rear stagnation point and the streamlines round a flat plate at incidence as Γ ramps from 0 to the Kutta value | the hypothesis as a process: the flow rearranges until it leaves the edge | frames (stop at Γ/Γ_K = 0, ½, 1) | closed form on a 120 × 80 grid, ≈ 3 s |
| A2 | **C03** | the lumped-vortex impulsive start: wake vortices leaving the edge, bound circulation growing, the total flat at zero | Kelvin's bookkeeping in time; the curve Γ(t) cannot be seen in a still | video | 300 steps, ≈ 2 s |
| A3 | **C04** | the circle's centre moving along a path (up, then back) and the mapped section morphing plate → arc → cambered airfoil | camber and thickness as two independent motions of one circle | video | closed form, < 1 s |
| A4 | **C08** (N37) | the trailing sheet of an elliptic loading rolling up into two tip vortices (rear view), from a cached point-vortex run | the sheet-to-pair picture behind "tip vortices"; labelled qualitative | frames (stop at four times) | ≈ 8 s uncached with n = 200; FAST n = 80 ≈ 2 s; cached in `outputs/ch14/cache/` |
| A5 | **C16** (N77, N82) | a fish tail through one stroke with the relative velocity, lift and thrust arrows; then a boat tacking upwind with its drive arrow | lift ⟂ relative velocity gives thrust on both half-strokes and on both tacks | video | closed form, < 2 s |

Plotly slider figures (9, `slider_figure`; ≤ 40 positions × ≤ 4 traces × ≤ 400 points; they survive publishing).
| # | A | Slider → what changes | Why |
|---|---|---|---|
| F1 | **C01** | taper ratio (and a second figure: aspect ratio at fixed area) → planform outline, with s, A, $\bar c$, Λ in the title | the vocabulary as one picture that moves |
| F2 | **C04** | thickness parameter e → section shape, measured $t_{\max}/c$ against 1.3 e (a second slider figure for β and the camber rule) | the rule and its error together |
| F3 | **C05** | circulation Γ/Γ_K → speed along the upper surface near the trailing edge (log axis) | the singularity of (14.9) appearing and disappearing |
| F4 | **C06** | camber angle β → the lift line and the exact curve against α in degrees, zero-lift angle marked | slope fixed, intercept moves |
| F5 | **C07** | angle of attack α → $-C_p$ against x/c with the loop shaded and its area in the legend | the area tracks $C_L$ |
| F6 | **C09** | number of horseshoes N → staircase Γ(y) and its downwash against the smooth limits | the integral (14.13) as a limit of sums |
| F7 | **C12** | taper ratio → Γ(y), w(y) and $e_{\text{span}}$ from the solver, elliptic ghost | taper near 0.4 is almost elliptic |
| F8 | **C14** / **C15** | aspect ratio Λ → the polar $C_L$ against $C_D$ and (second figure) drag against speed | slender wings: flatter polar, slower best glide |
| F9 | **C16** | lift-to-drag ratio of the sail → drive against course angle, no-go zone shaded | the no-go zone narrows as L/D grows |

Live widgets (2, kernel only, each paired with a slider figure above): **C06** α, β, e → section + streamlines + $C_L$ (paired with
F4); **C12** taper, twist, Λ, N → loading and coefficients (paired with F7). Static figures carry every other A item: the labelled
planform and section (C01), the Kelvin circuit (C03), the four circles and their images (C04), the two planes with the same
particles (C05), both panels of the pressure figure (C07), the staircase of horseshoes and the sign map of w (C08, C09), the force
triangle (C10), the sine shapes and the printed-against-corrected sign check (C11), the convergence table (C12), coefficient bars
(C13), three planforms with their loadings (C14), the generic curves and the drag shares (C15), the sailing polygon (C16).

## 7. From-scratch moments
Each shows a hand-written version beside the tested function, followed by `assert np.allclose(...)`. At least one per book section
that has a computable A item (§14.2, §14.3, §14.4, §14.5, §14.6, §14.7, §14.9).
| A | § | Hand-written | Compared with |
|---|---|---|---|
| **C01** | 14.2 | planform area by a trapezoid loop over c(y), then Λ and $\bar c$ in two lines | `LL.planform`, `LL.aspect_ratio`, `LL.mean_chord` |
| **C02** | 14.3 | the Kutta circulation of a flat plate found by a ten-line bisection on "is the speed at ζ = b zero?", before any formula | `AF.kutta_circulation_root`, later `AF.kutta_circulation` |
| **C03** | 14.3 | the lumped-vortex march as an explicit loop of about fifteen lines (solve for the bound vortex, shed one, convect) | `ch14.impulsive_start_lumped_vortex` |
| **C04** | 14.4 | a circle as `mu + a*np.exp(1j*t)` and the map in one line; chord, camber and thickness by `max` / `min` | `AF.zhukhovsky_airfoil`, `AF.zhukhovsky_shape` |
| **C06** | 14.5 | Γ from (14.11) and $C_L=2\Gamma/(Uc)$ with the measured chord, in three lines | `AF.kutta_circulation`, `AF.zhukhovsky_lift_coefficient(exact=True)` |
| **C07** | 14.5 | $dw/d\zeta$ on the circle, the factor of (14.9), the speed and $C_p$ by Bernoulli, written out with complex numpy arrays | `AF.zhukhovsky_surface` |
| **C09** | 14.6 | downwash at a station by summing the half-line formula over N horseshoes in a loop | `LL.staircase_downwash`, `LL.downwash` (as N grows) |
| **C12** | 14.6 | the matrix of (14.20) filled by two nested loops for N = 4 and solved with `np.linalg.solve` | `LL.solve_lifting_line`, `LL.lifting_line_matrix` |
| **C13** | 14.6 | L and $D_i$ by the trapezoid rule on (14.15) with w from the series | `LL.lift_from_series`, `LL.induced_drag_from_series` |
| **C15** | 14.7 | the best glide by a brute-force scan of $C_L/C_D$ over $C_L$ | `ch14.best_glide` |
| **C16** | 14.9 | the apparent wind by subtracting two vectors and the drive by two dot products | `ch14.apparent_wind`, `ch14.sail_forces` |

## 8. Notes for the implementer
**Functions the A items' figures and the explainers need beyond analysis §4** (B/C functions optional):
1. `AF.zhukhovsky_circle_from_centre(mu, b)` → (e, β, a): the inverse of `zhukhovsky_circle`, for E1's draggable centre; and
   `AF.zhukhovsky_circle` must also accept (a, b, β) (trap T7).
2. `AF.stagnation_points(U, b, e, beta, alpha, *, Gamma_cw)` → both stagnation points on the circle (angles from the centre Q) and
   their images, with a flag when they have left the surface; `rear_stagnation_point` may wrap it. E2, A1.
3. `ZhukhovskyFlow.speed(z)` and `.speed_near_edge(d)` (speed at distance d from the trailing edge along the bisector) for F3 and
   E2's readout; `AF.zhukhovsky_surface` returns the analytic trailing-edge limit when Γ is the Kutta value and `np.inf` otherwise.
4. `LL.staircase_downwash(y1, y_edges, Gamma_steps)` → downwash from nested horseshoes by the half-line formula (a finite sum that
   tends to (14.13)); `LL.staircase_from_loading(Gamma_fn, s, N)`. E6, F6, the C09 from-scratch cell.
5. `LL.lifting_line_matrix(c, s, K, N, *, symmetric=False)` → (A, gamma) and `LL.lifting_line_rhs(...)`, so that the from-scratch
   cell and E7's JS mirror can be compared entry by entry; `solve_lifting_line` returns `cond`, `e_span`, `slope`.
6. `ch14.drag_share(CL, *, CD0, AR, e_span)` → induced fraction; `ch14.speed_min_drag(W, rho, A, AR, CD0, e_span)`. E8, N75.
7. `ch14.sail_drive_curve(theta, *, U, V, CL, CD, rho, A)` and `ch14.velocity_made_good(V, theta)`. E9, F9.
8. `ch14.wagner_function(s)` (the two-exponential fit; only after the verifier has read the source first-hand, analysis §8) as the
   ghost curve of A2 and E3; without it the ghost is the steady value only.
9. `ch14.illustrative_inputs()` → our example aircraft, section and boat (visibly not the book's), used by every worked number.
**Names.** Span efficiency is `e_span` everywhere (the analysis writes `e` in three signatures: `induced_drag_coefficient`,
`drag_polar`, `best_glide`); the thickness parameter keeps `e`. `K` is a required keyword of every LL function. Circulation
arguments are `Gamma_cw=` / `Gamma_ccw=`, returned values named `Gamma_cw`.
**Explainer JS.** `Viz.cx` and `Viz.num.lu` are not in `assets/viz_lib.js`; E1, E2, E4 need complex arithmetic (the `cx` helpers
of ch06 `conformal_joukowski`) and E7 a small LU — copied locally and listed for promotion by the knowledge-keeper (ch06 rank 1 and
rank 11 in `knowledge/viz_patterns.md`).

## 9. Implementation guidance (conventions, runtime budget, caching)
- **Conventions in code.** Angles in radians (`*_deg` wrappers at the interface only); x along the mapped chord in AF functions with
  the stream inclined at α, so forces are returned both as (Fx, Fy) and as (L, D); in LL functions y is spanwise and w positive
  downward; scalars in → floats out so explainer parity rows can call every closed form; the inverse map chooses the root by
  $\lvert\zeta-\mu\rvert\ge a$ (trap T8). Printed forms that are slips stay available as named options that a test must fail
  (`downwash_series(printed_middle=True)`).
- **Runtime budget (< 5 min on Colab CPU; target ≈ 90–120 s).** Everything is closed form or a small linear solve except: four sympy
  check cells for the ★★★ derivations (≤ 4 s each; results cached with P252 in FAST); `ch14.circular_arc_sym` (< 1 s); the roll-up
  run for A4 (FAST n = 80, cached); the lumped-vortex march (≈ 2 s); the Thwaites separation sweep over 12 angles (≈ 3 s); five
  animations (≈ 15 s in FAST) and nine slider figures (≈ 12 s). No cached data file needs to be committed.
- **Evidence for the A items** (analysis §6): every computable CORE item already has two independent levels; C07's are V1 (stagnation
  $C_p=1$, the Blasius lift) and V4 (∮u_s ds = Γ); C15's are V1/V2 (closed form against numerical minimisation, dimensions); C16's V1
  and V7. The qualitative models (N13, N37 past the first turn, N70, N72) are NOTE items and go to the verification report's Open items.

## 10. Slips and traps (analysis §9) — where the notebook flags them
| Analysis §9 item | What | Flagged in | How |
|---|---|---|---|
| slip #1 | first two members of (14.19) printed with $\cos\gamma_1-\cos\gamma$; correct $\cos\gamma-\cos\gamma_1$ | **C11** (D14 step with the sign, a ⚠️ callout, a cell) | derive the corrected form; `downwash_series(printed_middle=True)` gives −Γ₁/2s for the ellipse, direct quadrature of (14.13) gives +Γ₁/2s; say that the last member and (14.20)–(14.26) are right |
| slip #2 | $\varepsilon=\tan(w/U)$ printed for $\tan^{-1}(w/U)$ | **C10** (N42, D13 step 1) | one line beside the force triangle |
| slip #3 | "at location y (< y₁)" against the figure caption's y > y₁ | **C09** (R12, D12) | state that the velocity is induced at y₁ and which filaments push down |
| slip #4 | (14.8) printed with two relation signs; coefficient 1.3 is $3\sqrt3/4$, first order in e | **C04** (N26, D04) | "rule against measured" table from `zhukhovsky_shape(exact=)` |
| slip #5 | Lanchester's birth year | **C03** (N17) | one clause with the cited year 1868 |
| slip #6 | (4.107)–(4.108) attributed to "Section 4.3" | **C01** (R01) | the recap points to §4.11 |
| slip #7 | upper limit N in (14.22), ∞ elsewhere | **C13** (D18) | one clause: (14.21)–(14.22) are written for the truncated solution |
| slip #8 | $c\simeq4a$ beside $c\simeq4b$; the slope exceeds 2π at finite thickness | **C06** (D09, F4) | exact against approximate curve; the gap in per cent at e = 0.1 |
| slip #9 | span written b in two exercises | **C12** (N63) | one clause |
| slip #10 | "Λ-pattern" of migrating birds | **C09** (N65) | one clause: a V formation, not the aspect ratio |
| T1 | Γ is clockwise | conventions box; **C02** (R04); every explainer readout | two-row sign table; `Gamma_cw=` / `Gamma_ccw=`; planted variant fails a test |
| T2 | $dw/dz=u-iv$ here, ch13 used $u+iv$ | conventions box; **C05** (R08, D07) | ⚠️ callout where (14.9) first appears |
| T3 | b, a, β, c, e, Λ, γ overloaded | conventions box; **C06** (N31), **C11** (N49), **C13** (N54) | a letter table with "here" and "elsewhere"; `e_span` |
| T4 | three different θ | conventions box; **C04** (N19), **C06** (R09, D09) | a sketch of each angle's origin |
| T5 | axes switch between sections; lift is not $F_y$; w positive downward | conventions box; **C01** (N08), **C07** (D10), **C09** (D12) | axes drawn at the top of §14.4–§14.5 and of §14.6 |
| T6 | "circulation unchanged by the map" and "equal velocities at infinity" are stated, not shown | **C05** (D08) | derivation written out; a cell integrates round three contours in both planes |
| T7 | which circle parameters | **C04** (N25) | our parametrisation with its two limits; the function accepts (a, b, β) too |
| T8 | 0/0 at the trailing edge; which root of the inverse map | **C05**, **C06** (N29), **C07** | the analytic limit; root chosen by distance from μ |
| T9 | "induced drag is the largest part of the drag" | **C15** (N75), E8 status | a computed share against speed; exactly half at best glide; the cited cruise figure |
| T10 | principal value; infinite downwash at a horseshoe's tips; the horseshoe's too-good drag | **C09** (D12, N60), **C12** (N52) | primer P344; stations avoid the tips; the artefact shown in numbers |
| T11 | elliptic loading is not elliptic planform | **C14** (N58, D20), E7 presets | ⚠️ callout; the twisted rectangle with one coefficient |
| T12 | K per radian in the equations, per degree in the prose | conventions box; **C06**, **C12** (N46) | primer P338; K a required keyword |
| T13 | which area A | **C01** (R01, N83) | one table: planform, frontal, chord per unit span |
| T14 | the Kutta condition is a model of a viscous outcome | **C02**, **C03** (N15) | "when it fails": stall, fast unsteady motion |
| T15 | Kelvin's theorem with viscosity inside the circuit | **C03** (D02) | the *why* of step 1 |
| T16 | separation estimate, roll-up, generic curves and regimes are qualitative models of ours | **C07** (N13), **C08** (N37), **C15** (N70, N72), E5 notes | the word "qualitative" in the figure title and the function label; listed as Open items |
| T17 | rule 9: private book values | decision 7; every worked number | `ch14.illustrative_inputs()`; no book number in the notebook or explainers |
| T18 | fish and sail formulas are quasi-steady and drag-free in the book | **C16** (N77, D22) | the efficiency with drag (ours) beside the ideal identity |

## Closing reconciliation (script-checked when this file was generated)
- §2: **113** rows, one per inventory row #1–#113, each with an ID, a depth and a tier word; every B/C row names an A parent (S01, S02 carry none). A 16 · B 77 · C 20; CORE 16 · NOTE 83 · RECAP 12 · SKIP 2.
- A count 16 is within [12, 18]; every A row has a one-line reason; `nbkit.curation_items('ch14')` reports CORE = 16.
- SEEN rows: 12 (#10, #11, #12, #14, #26, #27, #38, #40, #43, #50, #56, #60) — all RECAP (R01–R12); none is at depth A.
- §3: 9 section rows (14.1–14.9) plus one for exercises and literature; none empty.
- §4: 15 primers, P336–P350. §4b: 22 D rows (★ 4 · ★★ 14 · ★★★ 4) covering 29 analysis derivations; §4c: 11 statement rows; together all 40 of analysis §2b.
- §5: 9 explainers + 1 backup; every ★★★ D row whose A item has an explainer is listed in that explainer (D03, D04 → E1; D15, D18 → E7).
- §10: 10 slips and 18 traps, each with the block that flags it.

## Open questions for the orchestrator and the user
1. **Sixteen A items** (the analyst proposed 16 clusters; the hard maximum is 18). Two slots are free: the first candidate is the
   fish tail (N77) as its own block for §14.8, then the Kutta circulation (14.11) split from C06, then the horseshoe vortex (N60).
2. **§14.8 has no A item of its own** — its two B items open the C16 block and are a mode of E9. Say so if the fish tail should be
   a CORE block with its own explainer (that would make ten).
3. **Stall is not an A item** because nothing computable about it here is better than qualitative (trap T16); it has an explainer
   (E5, attached to C06) and three B items. A cited measured data set (analysis §8, open) would let N70 carry real points.
4. **C07 is "ours"** (the book shows the pressure figure but no formula for it). It is A because it is the bridge from the map to
   stall and to Ch. 15's critical Mach number; if only book results may be A, it becomes a B item of C06 and E4 attaches to C06.
5. **Three ★★★ derivations of the analysis are statements** (d-p finite edge speed, d-s leading-edge suction, d-aj wake energy):
   each is ours, belongs to a B item and is verified numerically instead. Any of them can be promoted to a written derivation.
6. **Wagner's function** (the ghost curve of the growth of the circulatory *lift* after the start — a lift ratio, not a circulation ratio; review M1) is used only if the verifier reads its source first-hand.
7. **Nine explainers, not ten.** The backup `wake_rollup` is ready to be the tenth if wanted (it would be labelled qualitative).
