# Chapter 14 — Aerodynamics: lesson design
(from `analysis/ch14_curation.md` (A 16 · B 77 · C 20 = 113 rows; CORE 16 · NOTE 83 · RECAP 12 · SKIP 2; 22 derivations D01–D22
(★ 4 · ★★ 14 · ★★★ 4: D03, D04, D15, D18); E1–E9 + backup B1; primers **P336–P350**) and `analysis/ch14.md` (§2 inventory with the
LaTeX the analyst read from the page images, §2b derivations d-a … d-an, §4 implementation rows I1–I33, §9 slips #1–#10 and traps
T1–T18). 2026-10-12, lesson-designer.

**Equation provenance (honest statement).** Every equation placed here is taken from the analyst's page-image transcription in
`analysis/ch14.md` §2. For this design I re-read 12 rendered pages against that transcription (PDF page = printed page + 27):
p719 ($\Lambda\equiv s^2/A=s/\bar c$, where $\bar c=A/s$ (14.1, 14.2)); p729 ($z=\zeta+\dfrac{b^2}{\zeta}$ (14.3));
p730 (the unnumbered $z=be^{i\theta}+be^{-i\theta}=2b\cos\theta$ and $z=Re^{i\theta}+\dfrac{b^2}{R}e^{-i\theta}$, and the lettering O, Q, P, A, B, C, D of the raised circle);
p731 ($x=(R+b^2/R)\cos\theta$, $y=(R-b^2/R)\sin\theta$ (14.4); $x^2\sin^2\theta-y^2\cos^2\theta=4b^2\sin^2\theta\cos^2\theta$ (14.5); the cosine rule
$QP^2=OP^2+OQ^2-2(OQ)(OP)\cos(\widehat{QOP})$ and $\dfrac{b^2}{\cos^2\beta}=R^2+b^2\tan^2\beta-2Rb\tan\beta\cos(90^\circ-\theta)$; $2b\tan\beta\,\sin\theta=R-b^2/R=y/\sin\theta$ (14.6);
$x^2+(y+2b\cot2\beta)^2=(2b\csc2\beta)^2$; $a\equiv b(1+e)$, $e\ll1$ (14.7));
p732 ($c\cong4b$, camber $\cong\tfrac12\beta c$, $t_{\max}/c\cong1.3e$ (14.8) — slip #4 seen: two relation signs after "camber");
p733 ($\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\dfrac{d\zeta}{dz}$ and the angles α, β, α + β at Q in the figure of the circle flow);
p734 ($\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9); $u_\theta=-2U\sin\theta-\dfrac{\Gamma}{2\pi a}$ (14.10); $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11);
$C_L=\dfrac{L}{(1/2)\rho U^2c}\simeq2\pi(\alpha+\beta)$ (14.12));
p738 ($dw(y_1)=\dfrac{-(d\Gamma/dy)\,dy}{4\pi(y-y_1)}$; $w(y_1)=\dfrac1{4\pi}\displaystyle\int_{-s/2}^{s/2}\dfrac{d\Gamma}{dy}\dfrac{dy}{(y_1-y)}$ (14.13) — slip #3 seen);
p739 ($\varepsilon=\tan\dfrac wU\simeq\dfrac wU$ as printed — slip #2 seen; $\alpha_e=\alpha-\varepsilon=\alpha-\dfrac wU$ (14.14));
p740 ($dL=dL_e\cos\varepsilon=\rho U_e\Gamma\,dy\cos\varepsilon\simeq\rho U\Gamma\,dy$, $dD_i=dL_e\sin\varepsilon=\rho U_e\Gamma\,dy\sin\varepsilon\simeq\rho w\Gamma\,dy$;
$L=\displaystyle\int_{-s/2}^{s/2}\rho U\Gamma\,dy$, $D_i=\displaystyle\int_{-s/2}^{s/2}\rho w\Gamma\,dy$ (14.15); $C_L=K\Big[\alpha-\dfrac{w(y)}U+\beta(y)\Big]$ (14.16);
$\Gamma(y)=\dfrac K2Uc(y)\Big[\alpha-\dfrac{w(y)}U+\beta(y)\Big]$ (14.17));
p741 ($y=-(s/2)\cos\gamma$; $\Gamma=\sum_{n=1}^\infty\Gamma_n\sin(n\gamma)$ (14.18); the three members of (14.19), the last being
$w(y_1)=\dfrac1{2s}\sum_{n=1}^\infty n\Gamma_n\dfrac{\sin(n\gamma_1)}{\sin\gamma_1}$ — slip #1 seen: the first two are printed with $\cos\gamma_1-\cos\gamma$;
$\dfrac K2Uc(\alpha+\beta)=\sum_{n=1}^\infty\Big(1+\dfrac{nKc}{4s\sin\gamma}\Big)\Gamma_n\sin(n\gamma)$ (14.20); $L=\dfrac{\pi s}4\rho U\Gamma_1$ and $D_i=\dfrac\pi8\rho\sum_{n=1}^Nn\Gamma_n^2$ (14.21, 14.22);
$\Gamma=\Gamma_1\sin(\gamma)=\Gamma_1\sqrt{1-(2y/s)^2}$ (14.23); $w(y)=\Gamma_1/2s$ (14.24));
p742 ($D_i=\dfrac\pi8\rho\Gamma_1^2=\dfrac{2L^2}{\pi\rho U^2s^2}$ (14.25); $C_{D_i}=\dfrac{D_i}{(1/2)\rho U^2A}=\dfrac{C_L^2}{\pi(s^2/A)}=\dfrac{C_L^2}{\pi\Lambda}$ (14.26);
$dL=\rho U\Gamma\,dy=C_L\tfrac12\rho U^2c\,dy$ (14.27) — slip #6 seen: "Section 4.3").
**All 27 numbered equations of the chapter are on those twelve pages; no difference from the analyst's transcription was found.**
The other 26 pages were not re-read by the designer. Equations of earlier chapters are quoted in the written-out forms of
`analysis/ch04_design.md`, `ch05_design.md`, `ch06_design.md` and of this chapter's curation: $C_D\equiv\dfrac{D}{\frac12\rho U^2A}$, $C_L\equiv\dfrac{L}{\frac12\rho U^2A}$ (4.107), (4.108);
Kelvin's theorem $\dfrac{D\Gamma}{Dt}=0$ (5.8); the filament law $d\mathbf u=\dfrac{\Gamma\,d\mathbf l}{4\pi}\times\dfrac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ (5.17); the cylinder's surface speed
$u_\theta=-2U\sin\theta-\dfrac{\Gamma}{2\pi a}$ (6.37) and stagnation angle $\sin\theta=-\dfrac{\Gamma}{4\pi aU}$ (6.38); the Blasius theorem $D-iL=\dfrac{i\rho}2\oint\Big(\dfrac{dw}{dz}\Big)^2dz$ (6.60);
$L=\rho U\Gamma$ (6.62); $z=\zeta+\dfrac{b^2}\zeta$ (6.65); the circle flow $w=U\Big(\zeta+\dfrac{a^2}\zeta\Big)+\dfrac{i\Gamma}{2\pi}\ln(\zeta/a)$ (6.68). Builders render a page
(`tools/render_pages.py ch14 --eq 14.M`) before setting any equation not in this list.

Numbers marked *expect* were computed for this design from the inputs stated beside them, by a scratch numpy / scipy script and by calling the implementer's modules as they stood on 2026-10-12 (read-only); where both routes exist (section numbers, the lifting-line solver, the polar, the sail forces, the lumped-vortex march) they agree to the digits printed. The builder compares an executed cell against them and **reports**
any difference instead of editing the expectation. The four `check_src` listings of Part F (D03, D04, D15, D18) were executed
as written (sympy; each about 1 s) and passed.

**Part C is the contract both the implementer and the builders keep.** `fluidpy/core/airfoils.py`, `fluidpy/core/lifting_line.py` and `fluidpy/ch14_aerodynamics.py` were written by the implementer while this design was in progress; Part C honours their names, signatures and return keys as found on 2026-10-12 (probed by calling them) and flags the three small additions it asks for (C.5).

**Binding conventions for every builder (analysis §9, curation decisions 4–7, §8–§10).**
1. **Imports and aliases.** `from fluidpy import ch14_aerodynamics as ch14`; `from fluidpy.core import airfoils as AF,
   lifting_line as LL`. **`ch14` re-exports every public name of `AF` and `LL`** and the reused functions named in C.0 (the ch11–ch13
   pattern), so explainer parity rows write `ch14.<name>` only. Reused: `from fluidpy.core import conformal as CONF, potential as
   POT, biot_savart as BS, similarity as SIM, boundary_layer as BL, vortices as VORT`; `from fluidpy import ch06_ideal_flow as ch06`.
2. **Scalar-callable and parity-friendly** (the `shot.py` evaluator knows only `np`, `math`, `ch14`): every public closed form
   accepts Python floats and returns a float, a complex, a tuple or a `dict` of floats (arrays only when an array goes in).
3. **Units.** SI. Angles in radians inside functions; degrees only in prose, slider labels and through `np.deg2rad`. Slopes of
   lift curves per **radian** in every formula and function (`K=` keyword); "per degree" appears only in one conversion (P338).
   Forces on a section are per unit span [N/m]; forces on a wing are in newtons. Sea-level air: ρ = 1.225 kg/m³,
   ν = 1.46 × 10⁻⁵ m²/s (`core.statics`, USSA-1976).
4. **Two conventions on screen wherever a sign matters.**

| What | The book | The project elsewhere | In code | Shown as |
|---|---|---|---|---|
| circulation Γ | **clockwise** positive, in $L=\rho U\Gamma$ (6.62), $u_\theta=-2U\sin\theta-\dfrac{\Gamma}{2\pi a}$ (14.10), $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11) | counterclockwise positive (ch03, ch05) | `Gamma_cw=` / `Gamma_ccw=`; returned values named `Gamma_cw` | "Γ (clockwise)" in every readout; a two-row sign table in C02 (trap T1) |
| complex velocity | $dw/dz=u-iv$ (ch06 and here) | ch13's Ekman velocity was $V=u+iv$ | `dwdz` returns u − iv; `velocity(x, y)` returns (u, v) | ⚠️ callout where $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9) first appears (trap T2) |
| axes | §14.3–§14.5: x along the mapped chord, y up, the stream inclined at α to x — **lift is not $F_y$**. §14.6: y along the span, w positive **downward** | — | AF returns (Fx, Fy) and (L, D); LL takes spanwise `y`, returns `w` positive down | axes drawn at the top of §14.4–§14.5 and of §14.6 (trap T5) |
| lift-curve slope K | per radian in $C_L=K\big[\alpha-\frac wU+\beta\big]$ (14.16)–(14.20); quoted per degree in §14.7 | — | `K=` required keyword of every LL function | primer P338 (trap T12) |
| lapse rate | not used in this chapter | Kundu Γ ≡ dT/dz with the meteorological −dT/dz beside it (standing rule) | — | **the letter Γ is circulation throughout this chapter**; the one place temperature enters (the sound speed in N01) uses T at sea level, no lapse rate. The conventions block says so in one line. |

5. **Overloaded letters (trap T3, T4) — the conventions block before C01 shows this table (`ch14.conventions_table()`); each block
   repeats the row it needs.**

| Letter | Meanings in this chapter (and earlier) | Notebook / explainer symbol | Code name |
|---|---|---|---|
| b | constant of the map $z=\zeta+b^2/\zeta$ (14.3) (chord ≃ 4b); the span in two exercises (slip #9) | b = map constant only; span is always s | `b`, `s` |
| a | radius of the ζ-plane circle (ch06: the cylinder radius); a horizontal length on the book's cambered-circle figure | a = circle radius, $a=b(1+e)/\cos\beta$ (ours) | `a` |
| β | camber angle of the circle (§14.4–§14.5); **minus the zero-lift angle** of a section (§14.6–§14.7); ch13: df/dy | β with "camber angle" or "−(zero-lift angle)" said each time; $\beta_a$ = apparent-wind angle (C16) | `beta`, `beta_a` |
| c | chord (ch13: wave speed; ch01: sound speed) | c chord, $\bar c$ mean chord, $c_{\text{sound}}$ | `c`, `c_mean` |
| e | thickness parameter in $a\equiv b(1+e)$ (14.7); span efficiency in aeronautics | e thickness parameter; $e_{\text{span}}$ span efficiency | `e`, `e_span` |
| Λ | aspect ratio $\Lambda\equiv s^2/A$ (14.1) (ch13: Rossby radius); the book's "Λ-pattern" of birds (slip #10); sweep angle | Λ aspect ratio; $\Lambda_{\text{sweep}}$ | `AR`, `sweep` |
| γ | span angle in $y=-(s/2)\cos\gamma$ (ch05: sheet strength; elsewhere c_p/c_v) | γ span angle; γ₁ the station where the downwash is wanted | `gamma`, `gamma1` |
| θ | (i) polar angle of ζ from the **origin** in $x=(R+b^2/R)\cos\theta$ (14.4); (ii) angle from the **stream-aligned diameter through Q** in $u_\theta=-2U\sin\theta-\dfrac\Gamma{2\pi a}$ (14.10); (iii) velocity-triangle or course angle (§14.8–§14.9) | θ with its origin named; in our own lines: $t$ for the angle round the circle from Q measured from the ξ-direction, $\theta_Q$ for (ii), $\theta_c$ for the course angle | `theta`, `t`, `theta_course` |
| w | complex potential $w=\phi+i\psi$ (§14.5); downwash (§14.6) | w with its name said; never both in one block | `w` (method), `w` (array) |
| A | planform area here; frontal area in ch09 and in the projectile exercise; chord per unit span for a section (trap T13) | A planform; $A_f$ frontal | `A` |
| Γ | circulation, clockwise positive; $\Gamma_n$ coefficients of $\Gamma=\sum\Gamma_n\sin(n\gamma)$ (14.18); $\Gamma_0$ mid-span value; $\Gamma_K$ the Kutta value (ours) | as listed | `Gamma_cw`, `Gn`, `Gamma0`, `Gamma_kutta` |

6. **Our circle parametrisation (trap T7; "not in the book — ours").** Centre $\mu=b-ae^{-i\beta}$, radius $a=b(1+e)/\cos\beta$: the circle
   always passes through ζ = b. β = 0 gives the book's symmetric case (centre −be, $a\equiv b(1+e)$ (14.7)); e = 0 gives the
   raised circle (centre $ib\tan\beta$, a = b/cos β). Points of the circle are $\zeta=\mu+ae^{it}$ with **t = −β at the trailing edge**,
   increasing counterclockwise. `AF.zhukhovsky_circle` also accepts `a=` directly.
7. **Book slips taught in corrected form** — house form `> ⚠️ **slip #k — the book prints** … **; the correct form is** …` where used;
   keys `ch14.book_slips()`.

| Slip | The book prints | Correct | Taught in | Planted wrong variant a test must fail |
|---|---|---|---|---|
| #1 | the first two members of (14.19) with the denominator $\cos\gamma_1-\cos\gamma$: $w(y_1)=\dfrac1{2\pi s}\displaystyle\int_0^\pi\sum n\Gamma_n\dfrac{\cos(n\gamma)\,d\gamma}{\cos\gamma_1-\cos\gamma}$ | $\cos\gamma-\cos\gamma_1$; the last member $w(y_1)=\dfrac1{2s}\sum n\Gamma_n\dfrac{\sin(n\gamma_1)}{\sin\gamma_1}$ and everything after it are right | C11 (D14 step 5, a ⚠️ callout, a cell) | `LL.downwash_series(printed_middle=True)` gives −Γ₁/2s for the ellipse |
| #2 | $\varepsilon=\tan\dfrac wU$ | $\varepsilon=\tan^{-1}\dfrac wU$ | C10 (N42, D13 step 1) | `LL.induced_angle(..., small=False)` against tan |
| #3 | the induced velocity is "at location y (< y₁)" | it is induced **at y₁**; a filament at y > y₁ on the starboard half pushes down, as the figure caption says | C09 (R12, D12 step 4) | — |
| #4 | "camber = ≅ ½βc" and the coefficient 1.3 in $t_{\max}/c\cong1.3e$ (14.8) | camber ≃ ½βc; the coefficient is $3\sqrt3/4=1.299$, first order in e | C04 (N26, D04) | `AF.zhukhovsky_shape(exact=)` shows both |
| #5 | Lanchester's birth year | 1868 (cited biography, analysis §8) | C03 (N17) | — |
| #6 | $C_D\equiv\dfrac D{\frac12\rho U^2A}$, $C_L\equiv\dfrac L{\frac12\rho U^2A}$ (4.107), (4.108) "in Section 4.3" | §4.11 | C01 (R01) | — |
| #7 | the upper limit N in $D_i=\dfrac\pi8\rho\sum_{n=1}^Nn\Gamma_n^2$ (14.22), ∞ in (14.18)–(14.20) | either, consistently: (14.21)–(14.22) are written for the truncated solution | C13 (D18 step 9) | — |
| #8 | "the chord length c (≃ 4a)" beside $c\cong4b$ (14.8) | both are first-order statements; exactly $C_L=8\pi(a/c)\sin(\alpha+\beta)$ and the slope exceeds 2π at finite thickness | C06 (D09 step 9, F4) | `AF.zhukhovsky_lift_coefficient(exact=True)` |
| #9 | the span written b in two exercises | s | C12 (N63) | — |
| #10 | the "Λ-pattern" of migrating birds | a V formation (not the aspect ratio) | C09 (N65) | — |

8. **Traps (true as printed, easy to misread) — each a `> ⚠️ Common confusion:` callout where it first bites.** T1 Γ clockwise
   (conventions block; C02 R04; every readout) · T2 $dw/dz=u-iv$ (conventions block; C05) · T3 overloaded letters (conventions
   block; C06 N31, C11 N49, C13 N54) · T4 three θ's (conventions block; C04 N19, C06 R09 and D09) · T5 axes switch; lift is not
   $F_y$; w positive downward (conventions block; C01 N08, C07 D10, C09 D12) · T6 "circulation unchanged" and "equal at infinity"
   are stated, not shown (C05 D08) · T7 which circle parameters (C04 N25) · T8 0/0 at the trailing edge and the root of the
   inverse map (C05, C06 N29, C07) · T9 "induced drag is the largest part" (C15 N75; E8 status) · T10 principal value, infinite
   downwash at a horseshoe's tips, the horseshoe's too-good drag (C09 D12 and N60; C12 N52) · T11 elliptic loading is not
   elliptic planform (C14 N58, D20; E7 presets) · T12 K per radian (conventions block; C06; C12 N46) · T13 which area A (C01 R01,
   N83) · T14 the Kutta condition is a model of a viscous outcome (C02; C03 N15) · T15 Kelvin's theorem with viscosity inside the
   circuit (C03 D02 step 2) · T16 qualitative models of ours (C07 N13; C08 N37; C15 N70, N72; E5 notes) · T17 private book values
   (every worked number uses our own inputs, header convention 11) · T18 fish and sail formulas are quasi-steady (C16 N77, D22).
9. **Colour code (curation §5; the same in figures, derivations, bars and explainers):** lift = teal · drag (induced, profile) =
   rose · thrust / drive = purple (accent) · circulation and vorticity = amber · pressure and suction = orange · the stream and
   velocities = blue · references and ghosts = muted.
10. **Honest labels.** (a) **"Not in the book — ours"** opens every appearance of: the circle parametrisation (N25); the surface
    pressure (C07); the exact lift coefficient $8\pi(a/c)\sin(\alpha+\beta)$; the finite trailing-edge speed (N29); the lumped-vortex
    start (C03's code); span efficiency (N54); the finite-wing slope $K/(1+K/\pi\Lambda)$ (N59); roll-up numbers (N37); wake energy
    (N38); the parabolic polar and best glide (C15, N74); the sail's limiting course (C16, N81). (b) **"Qualitative model"** is in
    the title of every figure made with `ch14.suction_side_separation`, `ch14.wake_rollup` past the first turn,
    `ch14.generic_lift_curve`, `ch14.generic_drag_curve`, `ch14.stall_regime` (trap T16). (c) The Wagner ghost curve (`ch14.wagner_function`, a published two-exponential fit, s in semi-chords) is drawn only if the verification report lists its source as read first-hand; otherwise the ghost is the steady value alone and no cell or explainer text names Wagner.
11. **Public-repo rule.** No book number is reproduced. **The notebook's worked numbers use the following input set of this design, written out literally in the imports cell as five small dicts (`sec`, `ac`, `wd`, `boat`, `fish`)** — `ch14.illustrative_inputs()` returns the module's own, different illustrative set, which is displayed once beside it (C.5):
    **section** b = 0.25 m (chord ≈ 1 m), e = 0.08, β = 4°, α = 5°, U = 30 m/s, ρ = 1.225 kg/m³ · **aircraft** W = 9000 N, s = 12 m,
    A = 15 m² (Λ = 9.6, $\bar c$ = 1.25 m), slow 30 m/s, cruise 55 m/s, $C_{D0}$ = 0.025, $e_{\text{span}}$ = 0.9, $C_{L,\max}$ = 1.4 clean and
    2.0 with flap, K = 2π · **demo wing** rectangular, Λ = 6 (s = 6 m, c = 1 m), α = 5°, U = 1 m/s, ρ = 1 · **boat** true wind 6 m/s,
    boat speed 3 m/s, course 45° off the wind, sail 20 m², $C_L$ = 1.2, $C_D$ = 0.2, hull L/D = 8, hull drag area 0.05 m², water 1000 kg/m³ · **fish** U = 1.0 m/s,
    V = 0.75 m/s, L = 10 N, D = 1 N. The stall angles and Reynolds-number thresholds used by the qualitative models are **ours**
    (`alpha_stall` = 14°, `CL_post` = 0.9, `width` = 2°, `CD_sep` = 0.15, `Re_laminar` = 5 × 10⁵, `Re_turbulent` = 3 × 10⁶) and are passed by keyword in every call.

---

## Part C — functions the builders will call (the implementer's contract)

Names and signatures: `fluidpy/core/airfoils.py` (AF) and `fluidpy/core/lifting_line.py` (LL) **as found on disk on 2026-10-12**
(probed by calling them; return keys copied from the returned dicts), and `fluidpy/ch14_aerodynamics.py` (ch14), which appeared while this design was being written and is honoured likewise; analysis §4 (rows I1–I33) and curation §8 (items 1–9), §9 are the background. "Eq." = the book equation implemented (written out in
the docstring). Items marked **(+)** are not in analysis §4 (C.5 lists them). Return shapes are part of the contract: scalars
in → Python floats (or complex) out; a `dict` has at least the keys listed; a tuple has exactly the members listed, in that
order. Everything after a bare `*` is keyword-only. Angles in radians. Γ arguments are `Gamma_cw=` (the book's clockwise
sign); returned circulations are clockwise.

### C.0 Reused (existing; called, not changed)
`CONF.joukowski(zeta, b=1.0)` ($z=\zeta+\dfrac{b^2}\zeta$ (14.3) = (6.65)), `CONF.joukowski_derivative(zeta, b)` ($dz/d\zeta=1-b^2/\zeta^2$),
`CONF.joukowski_inverse(z, b, branch="outside")`, `CONF.circle_flow_zeta(U, a, *, Gamma_cw=None, Gamma_ccw=None, alpha=0.0,
center=0j)` → (w(ζ), dw/dζ) (the off-centre circle at incidence, D06), `CONF.angle_preservation` · `POT.blasius_force(flow_or_dwdz,
R=None, contour=None, rho=1.2, n=256, center=0j)` ($D-iL=\dfrac{i\rho}2\oint\big(\dfrac{dw}{dz}\big)^2dz$ (6.60)), `POT.pressure_coefficient(speed, U)`,
`POT.power_branch(zeta, p, cut=np.pi)`, `POT.cylinder(U, a, *, Gamma_cw=…)` · `ch06.lift_per_span(rho, U, *, Gamma_cw=None,
Gamma_ccw=None)` ($L=\rho U\Gamma$ (6.62)), `ch06.cylinder_surface_speed(theta, U, a, *, Gamma_cw=…)` ($u_\theta=-2U\sin\theta-\dfrac\Gamma{2\pi a}$ (6.37) = (14.10)),
`ch06.cylinder_stagnation_points(U, a, *, Gamma_cw=…)` ($\sin\theta=-\dfrac{\Gamma}{4\pi aU}$ (6.38)) · `BS.segment_speed(d, theta_a, theta_b, Gamma)`
(the segment form $\dfrac{\Gamma}{4\pi d}(\cos\theta_a-\cos\theta_b)$ of the filament law (5.17)), `BS.segment_induced_velocity(x, xa, xb, Gamma, eps=0.0)`,
`BS.point_vortex_velocity(x, xv, Gamma, eps=0.0)`, `BS.point_vortex_evolve(xv0, Gamma, t_eval, …)`, `BS.point_vortex_invariants(xv,
Gamma)`, `BS.wall_image_system(xv, Gamma, wall_y=0.0)` · `SIM.lift_coefficient(F, rho, U, A)`, `SIM.drag_coefficient(F, rho, U, A)`
($C_L\equiv\dfrac L{\frac12\rho U^2A}$, $C_D\equiv\dfrac D{\frac12\rho U^2A}$ (4.107), (4.108)), `SIM.mach_number(U, c=None, T=None)`, `SIM.reynolds_number(U, l, nu)` ·
`BL.thwaites(x, Ue, nu, …)`, `BL.separation_point(x, tau0)`, `BL.blasius_drag_coefficient(Re_L, sides=1)` ·
`VORT.rankine_vortex(r, Gamma, sigma)` · `core.statics` (USSA-1976 sea-level density) · machinery: `setup_notebook`, `show_viz`,
`animate`, `show_animation`, `slider_figure`, `animate_figure`, `live`, `savefig`, `COLORS`. **`ch14` re-exports** `joukowski`,
`joukowski_derivative`, `circle_flow_zeta`, `lift_per_span`, `cylinder_surface_speed`, `segment_speed`, `lift_coefficient`,
`drag_coefficient` under the same names (row I32), so parity rows can reach them.

### C.1 `fluidpy/core/airfoils.py` (`AF`; exists; re-exported by `ch14`)

| Function (signature) | Returns [units] | Book eq. / source |
|---|---|---|
| `zhukhovsky_circle(b=1.0, e=0.0, beta=0.0, *, a=None)` | tuple (mu complex [m], a float [m]); $\mu=b-ae^{-i\beta}$, $a=b(1+e)/\cos\beta$; with `a=` given, e is ignored (trap T7) | $a\equiv b(1+e)$ (14.7); ours (N25) |
| `zhukhovsky_circle_from_centre(mu, b=1.0)` **(+)** | tuple (e, beta [rad], a [m]) — the inverse, for a dragged centre | ours (curation §8.1) |
| `zhukhovsky_airfoil(b=1.0, e=0.0, beta=0.0, n=400, *, a=None)` | dict: "t" (angle round the circle from Q; t[0] = −β is the trailing edge, counterclockwise), "zeta", "z", "x", "y", "R" (= abs ζ), "theta" (polar angle of ζ from the origin), "mu", "a", "b", "e", "beta" | $z=\zeta+\dfrac{b^2}\zeta$ (14.3); $x=(R+b^2/R)\cos\theta$, $y=(R-b^2/R)\sin\theta$ (14.4) |
| `zhukhovsky_inverse(z, b=1.0, mu=0j, a=None)` **(+)** | complex ζ outside the shifted circle (root chosen by abs(ζ − μ) ≥ a; NaN inside the body) | trap T8 |
| `section_geometry(x, y, te_index=None, m=0)` | dict: "chord", "camber", "x_camber", "t_max", "x_tmax", "camber_line", … [m] | N07 |
| `trailing_edge_angle(x, y, te_index=None)` | included angle at the trailing edge [rad] (→ 0 for a cusp as n grows) | N27 |
| `zhukhovsky_shape(b=1.0, e=0.0, beta=0.0, exact=True, *, a=None, n=20001)` | dict: "chord", "camber", "t_max", "x_tmax" [m], "camber_ratio", "thickness_ratio", "x_tmax_over_c", "x_camber_over_c", "e", "beta", "a", "exact"; `exact=False` returns the rules | $c\simeq4b$, camber $\simeq\tfrac12\beta c$, $t_{\max}/c\simeq1.3e$ (14.8) |
| `circular_arc(b=1.0, beta=0.1)` | dict: "centre" (complex, $-2ib\cot2\beta$), "radius" ($2b\csc2\beta$), "height" ($2b\tan\beta$), "end_angle" (2β), "chord" (4b), "a" | $x^2+(y+2b\cot2\beta)^2=(2b\csc2\beta)^2$ (unnumbered, N24) |
| `zero_lift_angle(camber_ratio)` · `camber_from_beta(beta, c)` | $\alpha_{L=0}=-2\,\text{camber}/c$ [rad] · camber $=\tfrac12\beta c$ [m] | N31; (14.8) |
| `resolve_lift_drag(Fx, Fy, alpha)` | tuple (L, D): $L=-F_x\sin\alpha+F_y\cos\alpha$, $D=F_x\cos\alpha+F_y\sin\alpha$ | N08 (ours) |
| `section_coefficient(F_per_span, rho, U, c)` | $F'/(\tfrac12\rho U^2c)$ [–] | R01, P336 |
| `circulation_from_cl(CL, U, c)` · `cl_from_circulation(Gamma_cw, U, c)` | $\Gamma=\tfrac12UcC_L$ [m²/s] · $C_L=\Gamma/(\tfrac12Uc)$ | N47 |
| `normal_force_coefficient(xc_upper, cp_upper, xc_lower, cp_lower)` | $C_y=\int_0^1(C_{p,l}-C_{p,u})\,d(x/c)$ (trapezoid on each surface) | N09, D10 |
| `lift_from_cp(xc_upper, cp_upper, xc_lower, cp_lower, alpha, suction="ideal")` | $C_L$: `"ideal"` → $C_y/\cos\alpha$ (zero drag, leading-edge suction implied); `"none"` → $C_y\cos\alpha$ | D10 step 7, N12 |
| `kutta_circulation(U, a, alpha, beta)` | $\Gamma=4\pi Ua\sin(\alpha+\beta)$ [m²/s], clockwise | (14.11) |
| `kutta_circulation_root(U, b, e, beta, alpha, xtol=1e-14)` | the same number by `brentq` on "u_θ at ζ = b is zero" (independent route) | C02; the Kutta condition |
| `stagnation_points(U, b, e, beta, alpha, *, Gamma_cw=None)` **(+)** | dict: "theta_rear", "theta_front" (angles $\theta_Q$ from the stream-aligned diameter through Q) [rad], "zeta_rear", "zeta_front", "z_rear", "z_front" (complex), "on_surface" (bool), "Gamma_cw", "Gamma_kutta", "rear_at_edge" | $u_\theta=-2U\sin\theta-\dfrac\Gamma{2\pi a}=0$ (14.10) |
| `rear_stagnation_point(U, a, alpha, beta, *, Gamma_cw)` | dict: "theta" ($\theta_Q$ of the rear point), "phi" (its angle from the ξ-direction at Q), "offset" (angle from B along the circle [rad]), "arc_offset" (= a × offset [m]), "on_surface", "at_edge" | N14, C02 |
| `trailing_edge_speed(U, a, b, alpha, beta)` | $U\,\dfrac ba\cos(\alpha+\beta)$ [m/s] (Kutta flow only) | N29 (ours) |
| `zhukhovsky_flow(U, b, e, beta, alpha, *, Gamma_cw=None, a=None)` → `ZhukhovskyFlow` | object: attributes `U`, `b`, `e`, `beta`, `alpha`, `mu`, `a`, `Gamma_kutta`, `Gamma_cw` (None in → the Kutta value), `is_kutta`; methods `w_zeta(zeta)`, `dwdzeta(zeta)`, `zeta_of(z)`, `w(z)`, `dwdz_from_zeta(zeta)`, `dwdz(z)` (= u − iv), `velocity(x, y)` → (u, v), `speed(z)`, `stream(x, y)` (= `psi`), `phi(x, y)`, `inside(x, y)`, `speed_near_edge(d)` (speed at distance d from the trailing edge along the cusp direction), `cp(z)` | N28 (D06); $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9) |
| `flat_plate_flow(U, b, alpha, *, Gamma_cw=None)` | `ZhukhovskyFlow` with e = β = 0 | N34 |
| `zhukhovsky_lift_coefficient(alpha, beta, e=0.0, exact=False)` | $2\pi(\alpha+\beta)$; `exact=True` → $8\pi(a/c)\sin(\alpha+\beta)$ with the measured chord | (14.12); ours |
| `zhukhovsky_surface(U, b, e, beta, alpha, n=400, *, Gamma_cw=None)` | dict: "t", "zeta", "z", "x", "y", "s" (arc length), "ds_dt", "u_theta" (on the circle), "u_s" (signed tangential speed on the section), "speed", "cp", "xc", "chord_x", "upper" and "lower" (each a dict "xc", "cp", "x", "y", "speed", "index", x/c increasing), "stagnation", "Gamma_cw", "is_kutta", "cp_min", "xc_cp_min", "U", "a", "b", "alpha", "beta", "e"; the trailing-edge value is the analytic limit when Γ is the Kutta value and `np.inf` otherwise | C07 (ours): $C_p=1-\lvert dw/dz\rvert^2/U^2$ |
| `surface_circulation(surface, method="polygon")` | $\oint u_s\,ds$ returned clockwise [m²/s] | N10 |
| `zhukhovsky_force(U, b, e, beta, alpha, rho=1.2, method="blasius", *, Gamma_cw=None, n=2000)` | dict: "Fx", "Fy", "L", "D" [N/m], "M0", "M_le" [N], "x_cp" [m], "xcp_over_c", "Gamma_cw", "CL", "CD", "chord", "method" (`"blasius"` or `"pressure"`) | R10: (6.60), (6.62) |
| `flat_plate_forces(rho, U, c, alpha)` | dict: "L" ($\pi\rho U^2c\sin\alpha$), "N" ($L\cos\alpha$), "S" ($L\sin\alpha$), "D" (0), "Gamma_cw", "CL" | N12 (ours) |
| `thin_airfoil(dyc_dx, alpha, n=200)` (optional) | dict: "A0", "A1", "A2", "CL", "alpha_L0", "CM_c4" | N18 (ours, standard) |

### C.2 `fluidpy/core/lifting_line.py` (`LL`; exists; re-exported by `ch14`)

| Function (signature) | Returns [units] | Book eq. / source |
|---|---|---|
| `planform(kind, s, c_root, taper=1.0, n=201)` (kinds "rectangular", "tapered", "elliptic") | dict: "y", "c" (arrays), "chord" (callable c(y)), "A", "A_trapezoid" [m²], "AR", "c_mean" [m], "kind", "s", "c_root", "taper" | N02; $\Lambda\equiv s^2/A=s/\bar c$, $\bar c=A/s$ (14.1), (14.2) |
| `aspect_ratio(s, A)` · `mean_chord(A, s)` | $s^2/A$ · $A/s$ | (14.1) · (14.2) |
| `span_angle(y, s)` · `span_position(gamma, s)` | γ from $y=-(s/2)\cos\gamma$ · its inverse | N49 |
| `circulation_series(Gn, gamma)` | $\sum\Gamma_n\sin(n\gamma)$ [m²/s] | (14.18) |
| `trailing_sheet_strength(Gn, gamma, s)` | $-d\Gamma/dy$ [m/s] | C08, D11 |
| `elliptic_circulation(y, s, Gamma1)` · `elliptic_downwash(Gamma1, s)` | $\Gamma_1\sqrt{1-(2y/s)^2}$ · $\Gamma_1/2s$ | (14.23) · (14.24) |
| `downwash(y1, Gamma_fn, s, method="quad", h=1e-3)` | w(y₁) [m/s], positive down — principal value of the integral, singular part subtracted | $w(y_1)=\dfrac1{4\pi}\displaystyle\int_{-s/2}^{s/2}\dfrac{d\Gamma}{dy}\dfrac{dy}{(y_1-y)}$ (14.13) |
| `downwash_series(Gn, gamma1, s, printed_middle=False)` | $\dfrac1{2s}\sum n\Gamma_n\dfrac{\sin(n\gamma_1)}{\sin\gamma_1}$ (the value n(±1) at the tips); `printed_middle=True` evaluates the member as printed (opposite sign — slip #1) | (14.19), last member |
| `glauert_integral(n, gamma1, method="closed")` | $\pi\sin(n\gamma_1)/\sin\gamma_1$; `method="quad"` → quadrature with the singular part subtracted | N51, D15 (ours) |
| `induced_angle(w, U, small=True)` · `effective_angle(alpha, w, U)` | $w/U$ (or $\tan^{-1}(w/U)$) · $\alpha-w/U$ | N42 · (14.14) |
| `strip_forces(Gamma, w, U, rho, dy, exact=False)` | dict: "dL", "dDi", "dLe" [N], "eps" [rad], "Ue" [m/s] | N44, D13 |
| `lift_and_induced_drag(y, Gamma, w, rho, U)` | tuple (L, Di) [N] by the trapezoid rule | (14.15) |
| `lift_from_series(Gn, s, rho, U)` · `induced_drag_from_series(Gn, rho)` | $\dfrac{\pi s}4\rho U\Gamma_1$ · $\dfrac\pi8\rho\sum n\Gamma_n^2$ [N] | (14.21) · (14.22) |
| `span_efficiency(Gn)` | $e_{\text{span}}=\Gamma_1^2/\sum n\Gamma_n^2$ | N54 (ours) |
| `section_lift_coefficient(alpha, w, U, beta=0.0, *, K)` | $K[\alpha-w/U+\beta]$ | (14.16) |
| `lifting_line_modes(N, symmetric=False)` **(+)** | tuple (n array, gamma array): mode numbers and stations $\gamma_j=j\pi/(N+1)$ | N52 |
| `lifting_line_matrix(c, s, K, N, *, symmetric=False)` **(+)** | tuple (A (N × N), gamma (N,)): $\mathsf A_{jn}=\Big(1+\dfrac{nKc_j}{4s\sin\gamma_j}\Big)\sin n\gamma_j$ | N52; (14.20) |
| `lifting_line_rhs(c, s, U, alpha, K, N, *, beta=0.0, twist=None, symmetric=False)` **(+)** | array (N,): $r_j=\dfrac K2Uc_j(\alpha_j+\beta_j)$ | (14.20), left side |
| `lifting_line_residual(Gn, gamma, c, s, U, alpha, *, K, beta=0.0, twist=None)` | residual of $\Gamma(y)=\dfrac K2Uc\Big[\alpha-\dfrac wU+\beta\Big]$ (14.17) at the given γ | N48 |
| `solve_lifting_line(c, s, U, alpha, *, K, beta=0.0, twist=None, N=20, symmetric=False, rho=1.0)` (c: scalar, callable c(y) or a `planform` dict) | dict: "Gn", "n", "gamma", "y", "Gamma", "w", "L", "Di", "A", "AR", "CL", "CDi", "e_span" (alias "e"), "slope" [1/rad], "cond", "matrix", "rhs", "Gamma1" | (14.20), (14.21), (14.22) |
| `elliptic_induced_drag(L, rho, U, s)` | $\dfrac{2L^2}{\pi\rho U^2s^2}$ [N] | (14.25) |
| `induced_drag_coefficient(CL, AR, e_span=1.0)` | $C_L^2/(\pi\Lambda e_{\text{span}})$ | (14.26) |
| `finite_wing_lift_slope(K, AR)` | $K/(1+K/(\pi\Lambda))$ [1/rad] | N59, D20 (ours) |
| `elliptic_wing(alpha, AR, *, K, beta=0.0)` | dict: "CL", "CDi", "slope", "w_over_U", "alpha_e" | D20 |
| `elliptic_twist_rectangular(alpha_r, c, s, *, K)` **(+)** | dict: "alpha_t" ($=Kc\alpha_r/4s$) [rad], "Gamma1_over_U" ($=\tfrac12Kc\alpha_r$) [m], "CL", "CDi", "AR"; the twist is $\alpha(y)=\alpha_t+\alpha_r\sqrt{1-(2y/s)^2}$ | N63 (ours) |
| `twist_for_elliptic_loading(c, y, s, *, K, CL_target, A=None)` | dict with the twist distribution for any chord array | N58, N63 |
| `horseshoe_velocity(x, y, z, Gamma, s, x_end=None)` | tuple (u, v, w) [m/s], true components with x downstream, y along the span, **z up** (a downwash is a negative w here) | N60, N65 |
| `horseshoe_downwash(x, y, Gamma, s, z=0.0)` **(+)** | downwash, positive **down** [m/s]; Γ/(πs) at the centre of the bound vortex | N60 |
| `horseshoe_wing(Gamma, s, rho, U)` | dict: "L", "w0", "Di_centre", "Di_elliptic_same_lift", "ratio" (= ½), "note" | N60 (trap T10) |
| `downwash_along_x(x, Gn, s, y1=0.0, include_bound=False)` | w(x)/… [m/s] along the flight direction: 0 far ahead, w at the line, 2w far behind | N41 |
| `ground_effect_factor(h, s)` | $1-\dfrac1{1+(4h/s)^2}$ | N66 (ours) |
| `staircase_from_loading(Gamma_fn, s, N, spacing="uniform")` **(+)** | tuple (y_edges (N + 1,), Gamma_steps (N,)) | E6, F6 |
| `staircase_downwash(y1, y_edges, Gamma_steps, parts=False)` **(+)** | w(y₁) [m/s]; `parts=True` → tuple (w, contributions (N + 1,)) one per trailing filament | E6, F6, C09 from-scratch |
| `rolled_up_pair(Gn, s)` | dict: "Gamma0", "spacing" ($b'=\frac1{\Gamma_0}\int\Gamma\,dy$), "descent_speed" ($\Gamma_0/2\pi b'$) | N37 (ours) |
| `wake_kinetic_energy_per_length(Gn, s, rho, m=256)` | $\tfrac12\rho\int\Gamma\,w_\infty\,dy$ with $w_\infty=2w$ [J/m = N] | N38 (ours) |

### C.3 `fluidpy/ch14_aerodynamics.py` (`ch14`; exists; re-exports every public name of C.1, C.2 and the reused names of C.0)

| Function (signature) | Returns [units] | Book eq. / source |
|---|---|---|
| `illustrative_inputs()` | flat dict — **the module's own illustrative set**, used by its scripts and tests (keys "b", "e", "beta", "alpha", "U_section", "W", "s", "A", "AR", "U", "rho", "nu", "CD0", "e_span", "K", "CL_max", "CL_max_flap", "taper", "alpha_stall", "CL_post", "stall_width", "CD_sep", "Re_laminar", "Re_turbulent", "U_wind", "V_over_U", "CL_sail", "CD_sail", "A_sail", "LD_hull", "CDA_hull", "rho_water", "theta", "U_fish", "V_tail", "L_tail", "LD_tail", "Gamma_v", "ell"). The notebook's worked numbers use the design's own literal inputs (header convention 11); see C.5 | rule 9, trap T17 |
| `book_slips()` · `traps()` · `conventions_table()` | dict keyed "1" … "10" · list of 18 dicts · `pandas.DataFrame` (the letter table of header convention 5) | analysis §9 |
| `compressibility_check(U, T=288.15)` | dict: "Ma", "density_change" (= ½Ma²), "incompressible_ok" (Ma < 0.3) | N01 |
| `level_flight(W, rho, U, A, CD=None)` | dict: "CL" ($=W/(\tfrac12\rho U^2A)$), "q" [Pa], "D" [N], "power" [W] (the last two None without CD) | N05, D01 |
| `stall_speed(W, rho, A, CL_max)` | $\sqrt{2W/(\rho A\,C_{L,\max})}$ [m/s] | N73, D01 |
| `surface_point(U, b, e, beta, alpha, t, *, Gamma_cw=None)` **(+, to be added)** | dict: "zeta_re", "zeta_im", "x", "y", "dwdzeta_abs", "map_factor_abs" ($=\lvert\zeta^2/(\zeta^2-b^2)\rvert$), "speed", "cp" at the point $\zeta=\mu+ae^{it}$ (a scalar wrapper round `zhukhovsky_flow(...).dwdzeta` and `.dwdz_from_zeta`; at t = −β with the Kutta circulation it returns the limit of `trailing_edge_speed`) | $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9) |
| `edge_speed(U, b, e, beta, alpha, d, *, Gamma_cw=None)` **(+, to be added)** | speed [m/s] at distance d [m] from the trailing edge (= `zhukhovsky_flow(U, b, e, beta, alpha, Gamma_cw=…).speed_near_edge(d)`) | N29, F3 |
| `map_point(zeta_re, zeta_im, b)` **(+, to be added)** | dict: "R", "theta", "x", "y" | $x=(R+b^2/R)\cos\theta$, $y=(R-b^2/R)\sin\theta$ (14.4) |
| `edge_flow_suction(rho, Uo, a, eps=1e-3, n=2000)` | dict (control-volume force round a plate edge and its parts) [N/m] | N12 |
| `trailing_edge_wedge(n, U=1.0, d=1.0)` | dict: "included_angle" ($2\pi(1-1/n)$), "half_fluid_angle", "velocity", "speed", "Cp" (callables), "Cp_edge", "n" | R07 |
| `circular_arc_sym()` | dict of sympy expressions and residuals (all 0), key "ok" | (14.4) → (14.5), (14.6), the arc |
| `suction_side_separation(alpha, e, beta, Re_c, n=400, side="suction", closure="falkner_skan")` | dict: "separated", "s_sep", "x_sep_over_c", "z_sep", "theta_sep", "xc_peak", "label" (= "qualitative"), "side", "xc", "Ue", "s", "lam" | N13 (ours) |
| `chord_reynolds(U, c, nu)` · `stall_regime(Re_c, *, Re_laminar, Re_turbulent)` | $Uc/\nu$ · one of "laminar separation first", "mixed", "transition first" — **qualitative** | N72 |
| `impulsive_start_lumped_vortex(U, c, alpha, t_end, dt, *, wake="rigid", shed_fraction=0.25, eps=None, rho=1.0, history=False)` | dict of arrays, one entry per step: "t", "s_chords" (= Ut/c), "Gamma_bound_cw", "Gamma_wake_cw" (the total shed so far), "total" (zero to round-off), "Gamma_steady" ($=\pi Uc\sin\alpha$, repeated), "ratio" (= bound/steady), "lift_circulatory"; and the final wake "wake_x", "wake_y", "wake_Gamma_cw"; "model" — the algorithm of C.6 item 4 | C03 (ours) |
| `wagner_function(s)` | $1-0.165e^{-0.0455s}-0.335e^{-0.3s}$ with s in **semi-chords** (= 2 × chords travelled) — a published fit; drawn in the notebook and in E3 only if the verification report lists its source as read first-hand | curation §8.8 |
| `wake_rollup(Gn, s, n=200, t_end=None, eps=None, n_t=41, fast=False, rtol=1e-7, atol=1e-9)` | dict with the frames, strengths and invariants — **qualitative** past the first turn | N37 |
| `tip_vortex_estimate(W, rho, U, s)` | dict: "Gamma0", "spacing", "descent_speed", "downwash", "induced_drag" (elliptic loading) | N37 |
| `vortex_pair_energy_per_length(Gamma, s, sigma, rho)` | energy per length of a pair of Rankine vortices [J/m] | N67 |
| `wake_vortex_rolling_moment(rho, U, K, s, c, Gamma_v, ell)` | dict with the rolling moment [N m] | N64 |
| `generic_lift_curve(alpha, *, K, beta, alpha_stall, CL_post, width)` | $C_L$: $K(\alpha+\beta)$ up to `alpha_stall`, then a smoothstep blend over `width` to the plateau `CL_post` — **qualitative** | N70 (ours) |
| `generic_drag_curve(alpha, *, CD0, AR, e_span, K, beta, alpha_stall, CL_post, width, CD_sep)` | $C_D(\alpha)$ — **qualitative** | N70 (ours) |
| `drag_polar(CL, *, CD0, AR, e_span=1.0)` | $C_{D0}+C_L^2/(\pi\Lambda e_{\text{span}})$ | C15 (ours, from (14.26)) |
| `drag_breakdown(CL, *, CD0, AR, e_span=1.0, friction_fraction=None)` | dict of **coefficients**: "profile", "induced", "total", "induced_share" (and "friction", "form", "pressure" when `friction_fraction` is given) | C15, N75 |
| `drag_share(CL, *, CD0, AR, e_span=1.0)` | induced fraction $C_{D_i}/C_D$ | N75 |
| `drag_vs_speed(W, rho, A, AR, CD0, U, e_span=1.0)` | dict in **newtons**: "profile", "induced", "total"; and "CL", "induced_share", "LD" | D21 |
| `speed_min_drag(W, rho, A, AR, CD0, e_span=1.0)` | best-glide speed [m/s] | D21 |
| `best_glide(CD0, AR, e_span=1.0)` | dict: "CL", "CD", "LD_max", "induced_share" (= 0.5), "glide_angle" | N74, D21 |
| `profile_drag_estimate(Re_c, t_over_c, *, regime, form_factor)` | $C_{D0}$ estimate: flat-plate friction on both faces × a stated form factor — **approximate** | C15 |
| `oscillating_foil(U, V, L, D=0.0)` | dict: "theta", "Ur", "thrust", "side_force", "power_out", "power_in", "efficiency" | N77, N79 |
| `apparent_wind(U, V, theta)` | dict: "Ur" [m/s], "beta_a" [rad], "vector" | C16, D22 |
| `sail_forces(U, V, theta, *, CL, CD, rho, A)` | dict: "beta_a", "Ur", "L", "D", "drive", "lateral", "F" [N], "drive_over_L" | C16, D22 |
| `sail_drive_curve(theta, *, U, V, CL, CD, rho, A)` | drive [N] over an array of course angles | F9, E9 |
| `velocity_made_good(V, theta)` | $V\cos\theta$ [m/s] | N82 |
| `sailboat_course_limit(LD_sail, LD_hull)` | $\varepsilon_A+\varepsilon_H=\tan^{-1}(1/LD_{\text{sail}})+\tan^{-1}(1/LD_{\text{hull}})$ [rad] | N81, D22 |
| `sailboat_equilibrium(U, theta, *, CL, CD, rho_air, A_sail, rho_water, CDA_hull, LD_hull=None)` | dict: "V" (steady boat speed, by `brentq` on drive − hull drag), "sails", "beta_a", "drive", "lateral", "hull_drag", "keel_lift" [N], "LD_hull_required", "eps_A", "eps_H" [rad], "closure", "vmg", "keel_holds" | N81 (ours) |
| `projectile_range_drag(m, rho, A, CD, V0, theta0, g=G0, numeric=False)` (optional) | range [m] | N83 |

### C.4 `scripts/`
Existing (implementer): `scripts/ch14_common.py` (`setup`, `save`, `streamlines(ax, flow, xlim, ylim, …)`), and the runnable demos
`ch14_zhukhovsky_airfoil.py` (N35), `ch14_kutta_condition.py`, `ch14_starting_vortex.py`, `ch14_pressure_lift.py`,
`ch14_horseshoe_wake.py`, `ch14_lifting_line.py`, `ch14_drag_polar.py`, `ch14_thrust_sailing.py`. **To be added: `scripts/ch14_drawings.py`**
— sketch and figure helpers imported by the notebook with `from ch14_drawings import *` (pure matplotlib on top of the
tested functions; each returns a `Figure`): sketches `draw_planform_labelled()`, `draw_three_axes()`, `draw_section_labelled(x, y,
geom)`, `draw_axes_conventions()`, `draw_kelvin_circuit()`, `draw_triangle_OQP()`, `draw_circle_angles()` (Q, B, α, β, the stream-aligned
diameter), `draw_tip_flow_rear_view()`, `draw_horseshoe_staircase(N)`, `draw_force_triangle()`, `draw_span_semicircle()`,
`draw_stalled_section()`, `draw_velocity_triangle(kind)` (kinds "fish", "bird_down", "bird_up", "sail"), `draw_sail_polygon(sf)`; composite
figures `fig_four_circles()`, `fig_three_stages(U, b, e, beta, alpha)`, `fig_pressure_two_ways(S)`, `fig_planforms_loadings()`. If the
orchestrator prefers no new script, the notebook-builder writes these as inline plotting cells — nothing else depends on
them. No cached data file is committed; the roll-up run caches to `outputs/ch14/cache/` (git-ignored).

### C.5 Not in `analysis/ch14.md` §4, and differences from the modules as found (flag list)
From curation §8, all present in the modules: `zhukhovsky_circle_from_centre`, the `a=` keyword, `stagnation_points`,
`ZhukhovskyFlow.speed` / `.speed_near_edge`, `staircase_downwash`, `staircase_from_loading`, `lifting_line_matrix`,
`lifting_line_rhs`, `drag_share`, `speed_min_drag`, `sail_drive_curve`, `velocity_made_good`, `wagner_function`,
`illustrative_inputs`. Found in the modules and adopted as they are: `zhukhovsky_inverse`, `lifting_line_modes`,
`horseshoe_downwash`, `elliptic_twist_rectangular`, `lift_from_cp(suction=)`, the coefficient form of `drag_breakdown`, the
force-balance form of `sailboat_equilibrium`, the keys of `impulsive_start_lumped_vortex`. **Requested additions (not in the
modules on 2026-10-12):** (1) three scalar wrappers in `ch14` so that explainer parity rows need no object method —
`surface_point`, `edge_speed`, `map_point` (signatures in C.3); (2) `scripts/ch14_drawings.py` (C.4). **One difference of
content, not of signature:** `ch14.illustrative_inputs()` returns the module's own illustrative values (a different aircraft,
section and boat from those of header convention 11). Both sets are ours and neither is the book's. This design's *expect*
numbers were computed with the set of header convention 11, which the notebook writes out literally in its imports cell;
`ch14.illustrative_inputs()` is displayed once beside it with the sentence "the module's scripts and tests use this second
set". If the orchestrator prefers a single set, either the implementer makes `illustrative_inputs()` return the values of
header convention 11, or the notebook-builder re-runs the cells with the module's set and reports the new numbers against
the *expect* entries — the storyboard does not otherwise depend on the choice.

### C.6 Parity conventions (the explainers reproduce these exactly)
1. `py:` expressions use only `ch14.<name>(…)`, `np.…`, `math.…`, numbers (complex literals allowed), strings, dict keys and
   integer indices (one row of E6 passes a `lambda`, a plain Python expression). Angles are written `np.deg2rad(5.0)`.
2. Every circulation in a row is clockwise; a row that passes Γ writes `Gamma_cw=`.
3. Constants the page could compute are computed in the page and pinned by a row: $3\sqrt3/4$ through
   `ch14.zhukhovsky_shape(1.0, 0.001, 0.0, exact=False)["thickness_ratio"]/0.001`; 2π through
   `ch14.zhukhovsky_lift_coefficient(1.0, 0.0)`. None is typed.
4. **The lumped-vortex march (E3, C03), exactly as `ch14.impulsive_start_lumped_vortex` with its defaults:** plate on the
   x-axis from 0 to c; stream $(U\cos\alpha,\ U\sin\alpha)$; bound vortex at (c/4, 0); control point P = (3c/4, 0). A clockwise vortex Γ at
   $(x_k,y_k)$ induces at (x, y) the velocity $\dfrac{\Gamma}{2\pi r^2}\,(y-y_k,\ -(x-x_k))$. Step k (time $k\,dt$): (i) move every wake vortex by
   $(U\cos\alpha,\ U\sin\alpha)\,dt$ (`wake="rigid"`); (ii) place a new vortex at $(c+\tfrac14U\,dt\cos\alpha,\ \tfrac14U\,dt\sin\alpha)$ (`shed_fraction=0.25`);
   (iii) solve the two equations "normal velocity at P is zero" — $U\sin\alpha-\dfrac{\Gamma_b}{\pi c}+v_{\text{new}}+v_{\text{wake}}=0$ — and "total circulation is
   zero" — $\Gamma_b+\Gamma_{\text{new}}+\sum\Gamma_k=0$ — for $\Gamma_b$ and $\Gamma_{\text{new}}$. *expect* (U = 30, c = 1, α = 5°, dt = 0.05 c/U; the designer's own
   loop and the module agree to all printed digits): "ratio" = 0.6463, 0.8156, 0.8792, 0.9388 after 1, 3, 5, 10 chords (indices
   19, 59, 99, 199); 90 % at about 6.2 chords; "total" zero to round-off.
5. E7's JS matrix equals `ch14.lifting_line_matrix(c, s, K, N)[0][j][n]` entry by entry (rows = stations, columns = modes).
6. E5's model curve copies the formula of `ch14.generic_lift_curve`: $C_L=(1-\sigma)K(\alpha+\beta)+\sigma\,C_{L,\text{post}}$ with σ the smoothstep of
   $(\alpha-\alpha_{\text{stall}})/\text{width}$ (0 below the stall angle, 1 above it plus the width).

### C.7 Test-only and notebook-only functions
**Test-only or one-line use:** `thin_airfoil` (one optional check line in N18), `edge_flow_suction` (one check line in N12),
`projectile_range_drag` (one optional line in N83), `vortex_pair_energy_per_length` (one line in N67),
`wake_vortex_rolling_moment` (one line in N64), `circular_arc_sym` (one "the engine agrees" line under D03),
`lifting_line_residual` (the convergence cell of C12), `twist_for_elliptic_loading` (tests; the notebook uses
`elliptic_twist_rectangular`), `trailing_edge_angle` (the Check of D05), `profile_drag_estimate` (one line in C15: "where a
$C_{D0}$ of this size comes from"), `lifting_line_modes` (inside the matrix; shown once), `camber_from_beta`, `cl_from_circulation`
(inverse helpers). **Every other function of C.1–C.3 is called by a storyboard row of Part A or mirrored by an explainer
of Part B.**

---

## Part A — notebook storyboard (`notebooks/build_ch14.py` → `notebooks/ch14_aerodynamics.ipynb`)

**One line per book section** (cell numbers are estimates; ≈ 430 cells in all):
- §14.1 → N01 (opens the C01 story; the conventions block) — cells ≈ 10–18
- §14.2 → C01 (N02 N03 N05 N07 N08 · N04 N06 N83 · R01 R02 · D01 · P336 P337 · F1) — cells ≈ 19–52
- §14.3 → C02 (N10 · N11 · R04 · A1), C03 (N14 N15 N16 · N17 · D02 · P339 · A2 · **E3**); C07 is taught after C06 — cells ≈ 53–110
- §14.4 → C04 (R05 R06 R07 · N18–N27 · D03 D04 D05 · P340 P341 P342 · A3 · F2 · **E1**) — cells ≈ 111–165
- §14.5 → C05 (N28 · R08 · D06 D07 D08 · F3), C06 (R09 · N29–N35 · D09 · P338 P343 · F4 · live · **E2** · **E5**), C07 (R03 R10 · N09 N12 N13 · D10 · F5 · **E4**) — cells ≈ 166–255
- §14.6 → C08 (R11 · N36 N37 N39 N40 N69 · D11 · A4), C09 (R12 · N41 N60 N61 N65 N66 · D12 · P344 · F6), C10 (N38 N42–N45 N67 · D13 · P350 · **E6**), C11 (N49 N50 N51 · D14 D15 · P346 P347), C12 (N46 N47 N48 N52 N63 N64 N68 · D16 · P345 · F7 · live), C13 (N53 N54 · D17 D18 · P348 · **E7**), C14 (N55–N59 N62 · D19 D20 · F8a) — cells ≈ 256–385
- §14.7 → C15 (N70–N76 · D21 · P349 · F8b · **E8**) — cells ≈ 386–410
- §14.8 → N77 N79 · N78 N80 (the fish tail and the bird open the C16 story) — cells ≈ 411–420
- §14.9 → C16 (N81 N82 · D22 · A5 · F9 · **E9**); S01, S02; summary — cells ≈ 421–445

Every CORE block follows: problem in plain words → idea → primers → maths (notes and derivations, Part F) → tiny example →
code (fluidpy) + "What does the code above do?" → from-scratch check → visual(s) → notes and "What would change if…". Code
drafts give intent + exact calls; the builder comments every line (novice grade, units, the equation written out next to its
number). *expect* = numbers the executed cell must print (hypotheses to be reported against, not edited). *see / read /
change* = the three figure notes. Note ids open every note in bold (**N09 [B]**). Wherever a draft names an equation, the
equation is written out with its number, and the builder keeps it so. The five literal dicts `sec`, `ac`, `wd`, `boat`, `fish` of the imports cell (header convention 11) supply every worked number.

**Re-entry rule (binding; `tools/nbkit.py`).** `nb.recap(...)` and `nb.section(...)` close the current CORE block. The builder
writes `nb.current_core = "Cnn"` on the line after every `nb.recap` that sits inside block Cnn (no section of this chapter
continues a block of the previous section, so no re-entry is needed after a section row). Rows marked **(re-enter Cnn)** say where. Reminders of earlier-chapter ideas that are not
inventory rows of this chapter (Kelvin's theorem, Stokes' additivity, point vortices, vortex sheets, Biot–Savart, images,
Bernoulli, Thwaites, the Rankine vortex) are `nb.md` cells in the 🔁 recap format without an id; they do not close the block.

### A.0 Front matter
1. `nb.title(big_idea=…, roadmap=[…16…], prerequisites=[…])`. **Big idea (draft):** "A wing weighing a few tonnes holds up an
   aircraft weighing hundreds. How? Chapter 6 gave half the answer: a body in an ideal stream feels a force at right angles to
   the stream whenever the flow circulates round it, $L=\rho U\Gamma$ (6.62) — but it left the circulation free. This chapter
   supplies the missing half in three sentences. **Lift is circulation.** **The sharp trailing edge chooses it** — viscosity,
   acting only there, lets exactly one circulation survive, and that gives a lift that grows by 2π for every radian of
   incidence. **The wing tips make you pay for it** — the circulation cannot stop at the tips, so it trails behind as two
   vortices whose downwash tilts the lift backward: a drag that exists even without friction, and that falls as the wing
   gets long and slender. The last two sections turn the first sentence round: a fish's tail, a bird's wing and a sail all
   move a wing *across* the stream so that its lift leans forward and becomes thrust." **Roadmap (one line per CORE):** C01 the
   wing in numbers · C02 the Kutta condition · C03 the starting vortex · C04 the Zhukhovsky airfoil · C05 the flow round it ·
   C06 its lift · C07 its surface pressure · C08 bound vortex and trailing filaments · C09 downwash · C10 induced drag · C11 the
   sine series · C12 the lifting-line equation · C13 lift and drag from the coefficients · C14 the elliptic loading · C15 the
   drag breakdown and the polar · C16 thrust from lift. **Prerequisites:** the complex potential, the cylinder with
   circulation, $L=\rho U\Gamma$ (6.62) and the map $z=\zeta+\dfrac{b^2}\zeta$ (6.65) (Ch. 6); Kelvin's theorem $\dfrac{D\Gamma}{Dt}=0$ (5.8), Helmholtz's
   theorems and the filament law $d\mathbf u=\dfrac{\Gamma\,d\mathbf l}{4\pi}\times\dfrac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ (5.17) (Ch. 5); the coefficients $C_L\equiv\dfrac L{\frac12\rho U^2A}$,
   $C_D\equiv\dfrac D{\frac12\rho U^2A}$ (4.107), (4.108) (Ch. 4); boundary layers and separation (Ch. 9).
2. `nb.explainer_index([...])` — 9 rows: ("zhukhovsky_airfoil_builder", "Where do camber and thickness come from?", "drag the
   circle's centre: up bends the plate, back fattens the nose") · ("kutta_condition", "Which circulation does a wing choose?",
   "only one puts the rear stagnation point on the sharp edge") · ("starting_vortex", "Where does a wing's circulation come
   from?", "it sheds an equal and opposite vortex and flies away from it") · ("airfoil_pressure_lift", "Is a wing pushed up or
   pulled up?", "the lift is the area between two pressure curves, mostly suction near the nose") · ("lift_curve_stall", "Does
   more angle always mean more lift?", "camber shifts the line, the angle buys 2π per radian, stall ends it") ·
   ("downwash_induced_drag", "How can there be drag without friction?", "the wing flies in its own sinking air, and its lift
   leans back") · ("lifting_line_solver", "What makes one wing better than another?", "lift lives in the first sine; every other
   one is pure drag") · ("elliptic_loading_aspect_ratio", "Why are glider wings so long and thin?", "the price of lift is set by
   the span you spread it over") · ("thrust_from_lift", "How can a boat sail against the wind?", "lift is perpendicular to the
   wind you feel, not the wind that blows").
3. `nb.setup()`.
4. `nb.code` — **chapter imports**: `import numpy as np`, `import sympy as sp`, `import matplotlib.pyplot as plt`, `import pandas
   as pd`, `from scipy.integrate import quad`, `from scipy.optimize import brentq`, `from fluidpy import ch14_aerodynamics as ch14`,
   `from fluidpy.core import airfoils as AF, lifting_line as LL`, `from fluidpy.core import conformal as CONF, potential as POT,
   biot_savart as BS, similarity as SIM, boundary_layer as BL, vortices as VORT`, `from fluidpy import ch06_ideal_flow as ch06`,
   `from fluidpy.core.interact import slider_figure, animate_figure, live`, `from fluidpy.core.anim import animate`, `from
   fluidpy.core.style import COLORS, savefig`, `import sys; sys.path.insert(0, "scripts"); from ch14_drawings import *`; then the design's inputs, literally: `sec = dict(b=0.25, e=0.08, beta=np.deg2rad(4.0), alpha=np.deg2rad(5.0), U=30.0, rho=1.225)`; `ac = dict(W=9000.0, s=12.0, A=15.0, U_slow=30.0, U_cruise=55.0, CD0=0.025, e_span=0.9, CL_max=1.4, CL_max_flap=2.0, K=2*np.pi)`; `wd = dict(s=6.0, c=1.0, alpha=np.deg2rad(5.0), U=1.0, rho=1.0, K=2*np.pi)`; `boat = dict(U_wind=6.0, V_boat=3.0, theta=np.deg2rad(45.0), A_sail=20.0, CL=1.2, CD=0.2, LD_hull=8.0, CDA_hull=0.05, rho_water=1000.0)`; `fish = dict(U=1.0, V=0.75, L=10.0, D=1.0)`; `nu_air = 1.46e-5`; `mu, a = AF.zhukhovsky_circle(sec["b"], sec["e"], sec["beta"])`; and `pd.Series(ch14.illustrative_inputs())` shown once ("the module's scripts and tests use this second illustrative set; both are ours").
   *explain:* one line per import ("`AF` and `LL` are the two new core modules; `ch14` re-exports both"). *expect:* mu =
   −0.0200 + 0.0189i m, a = 0.2707 m.
5. `nb.md` — **⚠️ Conventions in this chapter** (before C01): (a) `draw_axes_conventions()` — two small axes sketches: "§14.3–§14.5:
   the section in the x–y plane, y up, the stream inclined at α to x" and "§14.6: y along the span, the downwash w positive
   downward" (trap **T5**); (b) the circulation table of header convention 4 — "the book's Γ is **clockwise**; a wing with
   positive lift has Γ > 0" (trap **T1**); (c) $dw/dz=u-iv$ here, not $u+iv$ (trap **T2**); (d)
   `ch14.conventions_table()` (a DataFrame) — the letter table of header convention 5 (b, a, β, c, e, Λ, γ, θ, w, A, Γ) with
   "here" and "elsewhere" (traps **T3**, **T4**, **T13**); (e) "the lift-curve slope K is per radian in every formula" (trap
   **T12**); (f) "the letter Γ is the circulation everywhere in this chapter; no lapse rate appears"; (g)
   `pd.DataFrame(ch14.book_slips()).T` (10 rows) and `pd.DataFrame(ch14.traps())` (18 rows) with the sentence "slips are
   false as printed; traps are true but easy to misread"; (h) "**Which inputs?** All worked numbers use our own illustrative inputs (the five dicts above), never the book's" (trap **T17**).

### A.1 §14.1 Introduction — N01 (opens the C01 story)
1. `nb.section("14.1", "Introduction", intro="**What is this section about?** What 'aerodynamics' covers, the three
   assumptions the whole chapter leans on, and the one place where they fail.")`
2. `nb.note` — **N01 [B]** "Aerodynamics asks for the forces and moments a moving fluid exerts on a body. Three standing
   assumptions: (1) no buoyancy or stratification — the air's density is uniform (Ch. 6's ideal flow applies); (2) a uniform
   stream far ahead; (3) viscosity matters only in thin boundary layers on the surface (Ch. 9). **Stall** is where (3) fails:
   the boundary layer leaves the surface (N13, in C07). And the air may be treated as incompressible while the Mach number
   $\mathrm{Ma}=U/c_{\text{sound}}$ stays below about 0.3, because density changes are about $\tfrac12\mathrm{Ma}^2$ (ours)." + `nb.code`:
   `ch14.compressibility_check(30.0)`, `ch14.compressibility_check(55.0)`, `ch14.compressibility_check(240.0, T=220.0)` in a
   small table. *expect:* Ma = 0.088, 0.162, 0.807; density change 0.4 %, 1.3 %, 33 %; verdicts True, True, False. *explain:*
   "our light aircraft passes at both speeds; an airliner at cruise fails — that is Chapter 15." A three-row checklist
   (assumption · taught in · fails when). **🔁** one-line reminder: "Mach number and the speed of sound: Ch. 1 and Ch. 4 (C15)."

### A.2 §14.2 Aircraft Terminology — C01
#### C01 — The wing in numbers: aspect ratio (14.1), mean chord (14.2), the two coefficients and level flight
1. `nb.section("14.2", "Aircraft Terminology", intro="**What is this section about?** The handful of words and numbers that
   describe a wing — span, chord, area, aspect ratio, angle of attack — and the two coefficients every later result is
   written in.")`
2. `nb.core("C01", "The wing in numbers: $\\Lambda\\equiv s^2/A=s/\\bar c$ (14.1)", question="How do you describe a wing with a
   few numbers — and how much lift must it deliver to stay up?")`
3. `nb.md` — **Plain words:** "Stand under an airliner and under a glider. One wing is short and broad, the other long and
   thin. Before any physics we need words for that difference, and for the angle at which the wing meets the air. Then one
   balance — in steady level flight the lift equals the weight — tells us what the wing has to do, and everything else in
   the chapter is about *how* it does it and *what it costs*."
4. `nb.md` — **The idea** (ASCII top view): a planform with `|<———— s (span) ————>|`, `c(y)` marked at one station, `A = area seen
   from above`; beside it "slender: Λ large (glider ≈ 20) · stubby: Λ small (fighter ≈ 3)".
5. `nb.note` — **N02 [B]** root, tip, span s, chord c(y), planform area $A=\int_{-s/2}^{s/2}c(y)\,dy$ (ours, implied), taper, sweepback
   (one sentence each) + `nb.figure`: `draw_planform_labelled()` beside three outlines from `LL.planform("rectangular", 12, 1.25)`,
   `LL.planform("tapered", 12, 1.786, 0.4)`, `LL.planform("elliptic", 12, 1.592)` (all with A = 15 m²). *see / read* ("same span,
   same area, three ways of distributing the chord") */ change* ("…the taper goes to 0: a pointed tip — C12 shows what that
   does to the loading").
6. `nb.derivation("D01", …)` — Part F D01 (5 steps), ref "14.1": includes **N03 [B]** $\bar c=A/s$ (14.2) (step 1) and **N05 [B]**
   ($L=W$, $T=D$; steps 3–4) and the stall speed of N73 (step 5, pointer "used in C15").
7. `nb.primer("force per unit span, section coefficient and wing coefficient", "A two-dimensional section has no area, so its
   forces are per metre of span [N/m] and are divided by $\\tfrac12\\rho U^2c$ (chord c). A whole wing's forces are in newtons and
   are divided by $\\tfrac12\\rho U^2A$ (planform area A). $\\tfrac12\\rho U^2$ is the dynamic pressure: the pressure rise where the
   stream is brought to rest.", code="rho, U, c, A = 1.225, 30.0, 1.0, 15.0     # air, speed, chord, planform area\nq = 0.5*rho*U**2
   # dynamic pressure: 551 Pa\nprint(q, 587.0/(q*c), 9000.0/(q*A))   # a section's 587 N/m -> 1.06; a wing's 9000 N -> 1.09")`
   (**P336**)
8. **(re-enter C01)** `nb.recap("R01", "Lift and drag coefficients", "$C_D\\equiv\\dfrac{D}{\\tfrac12\\rho U^2A}$ and $C_L\\equiv\\dfrac{L}{\\tfrac12\\rho U^2A}$ (4.107),
   (4.108). **New here:** A is the *planform* area (Chapter 9's bluff bodies used the frontal area — trap T13), and a section
   uses the chord per unit span. > ⚠️ **slip #6 — the book prints** that these two definitions are 'in Section 4.3'**; the
   correct form is** §4.11.", where="Ch. 4 §4.11")` + a three-row table "which area?" (wing: planform A · bluff body: frontal
   $A_f$ · section: c × 1 m) + `nb.code`: `SIM.lift_coefficient(9000.0, 1.225, 30.0, 15.0)`, `AF.section_coefficient(587.0, 1.225,
   30.0, 1.0)`. *expect:* 1.088, 1.065.
9. **(re-enter C01)** `nb.recap("R02", "Friction drag and pressure drag", "Drag is the sum of the surface shear stress
   (friction drag) and of the surface pressure (pressure or form drag); lift is almost entirely a pressure force. C15
   refines this into friction + form + induced.", where="Ch. 9 (form drag); Ch. 4 (surface forces)")`
10. `nb.worked_example("a wing by hand", "s = 10 m, A = 20 m². 1. Mean chord $\\bar c=A/s$ (14.2) = 20/10 = 2 m. 2. Aspect ratio
    $\\Lambda\\equiv s^2/A$ (14.1) = 100/20 = 5 — and also $s/\\bar c$ = 10/2 = 5 ✓. 3. Level flight with W = 10 000 N at U = 40 m/s in
    air of ρ = 1.25 kg/m³: dynamic pressure ½ × 1.25 × 40² = 1000 Pa; $C_L=W/(\\tfrac12\\rho U^2A)$ = 10 000/(1000 × 20) = 0.5. 4. Halve
    the speed: the dynamic pressure falls to 250 Pa and the wing must deliver $C_L$ = 2.0 — four times as much. 5. If the wing
    can reach at most $C_{L,\\max}$ = 1.4, the slowest level flight is $U_{\\min}=\\sqrt{2W/(\\rho AC_{L,\\max})}$ = √(20 000/35) = 23.9 m/s.")`
11. `nb.code` — `LL.aspect_ratio(ac["s"], ac["A"])`, `LL.mean_chord(ac["A"], ac["s"])`, `ch14.level_flight(ac["W"], 1.225, 30.0,
    ac["A"])["CL"]`, the same at 55 m/s, `ch14.stall_speed(ac["W"], 1.225, ac["A"], ac["CL_max"])`, with `CL_max_flap`. *expect:*
    Λ = 9.6, $\bar c$ = 1.25 m, $C_L$ = 1.088 and 0.324, $U_{\min}$ = 26.5 m/s and 22.1 m/s. *explain:* 4 numbered points ("…4. flaps
    raise $C_{L,\max}$, so the aircraft can land 4.3 m/s slower — C15").
12. `nb.check_agree` — **from scratch:** `y = np.linspace(-6, 6, 201)`; `c = 1.786*(1 - 0.6*np.abs(y)/6)`; a trapezoid loop `A_mine
    += 0.5*(c[i] + c[i+1])*(y[i+1] - y[i])`; `AR_mine = 12**2/A_mine`; `cbar_mine = A_mine/12`; `P = LL.planform("tapered", 12.0,
    1.786, 0.4)`; `assert np.allclose([A_mine, AR_mine, cbar_mine], [P["A"], P["AR"], P["c_mean"]], rtol=1e-3)`. *expect:* A ≈ 15.0 m².
13. `nb.note` — **N07 [B]** the section: rounded leading edge, sharp trailing edge, chord line, camber line, camber $=\max y_c$,
    thickness $t=\max(y_u-y_l)$, angle of attack α between chord and stream + `nb.primer("splitting a closed contour into upper and
    lower surfaces", "A section outline is one closed list of points. To measure it, find the trailing edge (largest x) and the
    leading edge (smallest x), cut the list there, and interpolate both halves on the same x stations (`np.interp` needs x
    increasing, so one half is reversed).", code="x = np.array([1, .5, 0, .5, 1.]); y = np.array([0, .1, 0, -.06, 0])   # a
    diamond, counterclockwise from the tail\ni = np.argmin(x)                      # index of the leading edge: 2\nyu = np.interp(0.5,
    x[:i+1][::-1], y[:i+1][::-1]); yl = np.interp(0.5, x[i:], y[i:])\nprint(yu - yl, 0.5*(yu + yl))        # thickness 0.16, camber
    0.02")` (**P337**) + `nb.figure`: `draw_section_labelled()` on our section from `AF.zhukhovsky_airfoil(sec["b"], sec["e"],
    sec["beta"])` with the numbers of `AF.section_geometry(x, y)` in the labels. *expect:* chord 1.006 m, camber/chord 0.035,
    thickness/chord 0.096. *see / read / change* ("…α is not a property of the shape: tilt the stream, not the drawing").
14. `nb.note` — **N08 [B]** "The force on the section is resolved **relative to the stream**: drag D along it, lift L across it.
    If the stream makes the angle α with the x-axis, $D=F_x\cos\alpha+F_y\sin\alpha$ and $L=-F_x\sin\alpha+F_y\cos\alpha$ (ours).
    > ⚠️ Common confusion (trap T5): from §14.4 on the x-axis is the chord and the stream is tilted, so **lift is not $F_y$**."
    + `nb.code`: `AF.resolve_lift_drag(-51.1, 584.4, sec["alpha"])` *expect:* L = 586.6 N/m, D ≈ 0.0 N/m (these are the two force
    components C07 will compute; the pressure force leans *forward* of the y-axis by exactly α).
15. `nb.note` — **N04 [C]** pitch, roll and yaw axes (`draw_three_axes()`): "roll returns with the ailerons below and with the
    rolling moment in a wake vortex (N64, in C12)." · **N06 [C]** aileron, elevator, rudder, flap, slot, spoiler, trim tab —
    one sentence each: "all of them work by changing the local camber; 'more camber ⇒ more lift at the same α' is N31 in C06,
    flaps and landing speed N73 in C15." · **N83 [C]** "A projectile with a constant drag coefficient on a shallow trajectory
    has the range $x\simeq\dfrac{2m}{\rho AC_D}\ln\!\Big(1+\dfrac{\rho AC_DV_o^2\cos\theta_o\sin\theta_o}{mg}\Big)$ — here A is the *frontal* area (trap T13);
    `ch14.projectile_range_drag` (optional line)."
16. `nb.plotly` — **F1** two `slider_figure`s: (a) taper ratio 0.1…1 (19 steps) → the planform outline at fixed s = 12 m and
    A = 15 m², title carrying s, A, $\bar c$, Λ and root chord; (b) aspect ratio 3…20 (18 steps) at fixed area → the outline
    stretching, title "Λ = s²/A = …; span = … m". *see / read / change.*
17. `nb.md` — **What would change if…** "…you only knew the section's shape and the angle α — could you predict $C_L$? Not yet:
    ideal flow has one free number, the circulation. C02 is about who sets it."

### A.3 §14.3 Characteristics of Airfoil Sections — C02, C03 (C07 follows C06)
#### C02 — The Kutta condition: the rear stagnation point sits at the sharp trailing edge
1. `nb.section("14.3", "Characteristics of Airfoil Sections", intro="**What is this section about?** Why a section lifts at
   all. Ideal flow allows any circulation; the sharp trailing edge picks one (the Kutta condition); and the circulation it
   picks is paid for by a vortex left behind at the start. The section's pressure distribution — a §14.3 picture — is
   computed after the lift, in C07.")`
2. `nb.core("C02", "The Kutta condition: $\\left.\\dfrac{dw}{d\\zeta}\\right\\rvert_{\\zeta=b}=0$", question="Ideal flow round a wing has
   a solution for every circulation. Which one does a real wing choose?")`
3. `nb.md` — **Plain words:** "Hold a flat board at a slight angle in a stream. Chapter 6 says the force on it is $L=\rho U\Gamma$ (6.62)
   — and that the equations are satisfied for *any* Γ, including zero. A theory that cannot say whether a wing lifts is not
   much use. The missing fact is about the sharp rear edge: with the wrong circulation the air would have to whip round
   that edge from one side to the other, at a speed that grows without limit as the edge gets sharper. Real air cannot do
   that. It leaves the edge smoothly — and exactly one circulation lets it."
4. `nb.md` — **The idea** (ASCII, three frames of the rear of a plate at incidence): `Γ = 0: flow turns the edge ↺, rear
   stagnation point on top` → `Γ too small: point closer to the edge` → `Γ = Γ_K: flow leaves the edge smoothly`.
5. **(re-enter C02)** `nb.recap("R04", "The Kutta–Zhukhovsky lift theorem", "$L=\\rho U\\Gamma$ (6.62): a two-dimensional body in a
   steady ideal stream feels a force per unit span of ρUΓ at right angles to the stream, whatever its shape. **New here:** the
   theorem leaves Γ free; this chapter supplies it.", where="Ch. 6 §6.5")` + the two-row sign table (book: Γ clockwise, L = +ρUΓ ·
   project: $\Gamma_{ccw}$, L = −ρU$\Gamma_{ccw}$) + `nb.code`: `ch06.lift_per_span(1.225, 30.0, Gamma_cw=8.214)`,
   `ch06.lift_per_span(1.225, 30.0, Gamma_ccw=-8.214)`. *expect:* 301.9 N/m twice. `> ⚠️ Common confusion (trap T1): …`
6. `nb.md` — **The maths, step by step.** (i) "A flat plate of chord 4b is the image of the circle of radius b (R06, in C04); the
   trailing edge is the image of ζ = b." (ii) "On the circle the stagnation points sit where $u_\theta=-2U\sin\theta-\dfrac{\Gamma}{2\pi a}$ (6.37)
   vanishes: $\sin\theta=-\dfrac{\Gamma}{4\pi aU}$ (6.38), θ measured from the stream direction." (iii) "With the stream tilted up by α,
   the point ζ = b is at θ = −α." (iv) "So the rear stagnation point reaches the edge when the right-hand side equals
   sin(−α). We find that Γ numerically first; C06 turns it into a formula." (v) "Stated on the complex velocity: the Kutta
   condition is $\left.\dfrac{dw}{d\zeta}\right\rvert_{\zeta=b}=0$." 🔁 one-line reminder: "`brentq` finds the zero of a function between two brackets (P108, Ch. 3)."
7. `nb.note` — **N10 [B]** "Building lift means building circulation. In a real flow the circulation is the net vorticity held
   in the two boundary layers: $\Gamma=\oint u_s\,ds=\int_0^c(U_{\text{upper}}-U_{\text{lower}})\,dx$ (ours) — the air is faster over the top. (Worked with
   numbers in C07 below.)"
8. `nb.worked_example("three circulations on a unit circle", "U = 1 m/s, a = b = 1 m, α = 30°, so the edge is at θ = −30°. From
   $\\sin\\theta=-\\dfrac{\\Gamma}{4\\pi aU}$ (6.38): 1. Γ = 0: sin θ = 0, the rear stagnation point is at θ = 0° — 30° round from the edge, on the
   upper side. 2. Γ = π m²/s: sin θ = −π/(4π) = −0.25, θ = −14.5° — half-way there. 3. Γ = 2π m²/s: sin θ = −0.5, θ = −30° — on the
   edge. 4. Γ = 3π: sin θ = −0.75, θ = −48.6° — past the edge, on the lower side. Only the third leaves the edge alone.")`
9. `nb.code` — `F0 = AF.flat_plate_flow(sec["U"], sec["b"], sec["alpha"], Gamma_cw=0.0)`; `GK = AF.kutta_circulation_root(sec["U"],
   sec["b"], 0.0, 0.0, sec["alpha"])`; a loop over `f in (0, 0.5, 1.0, 1.2)`: `AF.rear_stagnation_point(sec["U"], sec["b"], sec["alpha"],
   0.0, Gamma_cw=f*GK)["offset"]` in degrees and `AF.zhukhovsky_flow(sec["U"], sec["b"], 0.0, 0.0, sec["alpha"],
   Gamma_cw=f*GK).speed_near_edge(1e-4)/sec["U"]`. *expect:* GK = 8.214 m²/s; offsets 5.0°, 2.5°, 0.0°, −1.0°; edge-speed ratios
   of order 10², 10², 1, 10 (the builder prints what it gets; only the f = 1 row is finite as d → 0). *explain:* 4 numbered
   points.
10. `nb.check_agree` — **from scratch (ten-line bisection, no formula):** `w, dw = CONF.circle_flow_zeta(U, b, Gamma_cw=G,
    alpha=alpha)` inside `def edge_slip(G): return (dw(b)*1j).real` (the tangential speed at ζ = b); `lo, hi = 0.0, 4*np.pi*U*b`;
    forty halvings; `assert np.isclose(0.5*(lo + hi), AF.kutta_circulation_root(U, b, 0.0, 0.0, alpha), rtol=1e-9)`. *expect:*
    8.2142 m²/s.
11. `nb.animation` — **A1** (`player="frames"`, 40 frames, FAST 20): streamlines round the plate from `AF.flat_plate_flow(…,
    Gamma_cw=f*GK).stream` on a 120 × 80 grid as f ramps 0 → 1, the rear stagnation point (amber dot) crawling to the edge;
    stop-frames at f = 0, ½, 1 (α = 15° for visibility). *see / read / change* ("…f goes past 1: the point passes the edge onto
    the lower side and the corner flow returns, the other way round").
12. `nb.figure` — speed at 10⁻⁴ m from the edge against Γ/Γ_K on a log axis (a V-shaped notch at 1). *see / read / change.*
13. `nb.note` — **N11 [C]** "Three ways to make and control lift: a sharp trailing edge (which lets the other two work), the
    angle of attack, and the camber. Angle and camber become the sum α + β in $C_L\simeq2\pi(\alpha+\beta)$ (14.12), C06."
    `> ⚠️ Common confusion (trap T14):` "The Kutta condition is not a theorem of ideal flow. It is a rule that stands in for
    what the boundary layers do, and it fails at stall and for motions faster than the time the air takes to cross the chord."
14. `nb.md` — **What would change if…** "…the foil starts from rest? Then Γ = 0 at first — where does the circulation come from? C03."

#### C03 — The starting vortex: the bound circulation is matched by a shed vortex of opposite sign
1. `nb.core("C03", "The starting vortex: $\\Gamma_{ABCD}=0=\\Gamma_{\\text{foil}}+\\Gamma_{\\text{starting}}$", question="Kelvin says circulation
   cannot be created in an inviscid fluid. So where does the wing's come from?")`
2. `nb.md` — **Plain words:** "Stir a spoon sideways through coffee and stop: a little whirl is left spinning where you began.
   A wing does the same. In the first moments of its motion it throws off a vortex from its trailing edge, and the air
   round the wing is left circulating the other way. Nothing was created — the two cancel — but the wing flies away from one
   half of the pair and keeps the other."
3. `nb.md` — **The idea** (ASCII): `[ A———————————————B ]` a big box round `↻ foil  …………  ↺ starting vortex`, cut by a line between
   them: "total 0 = (+Γ round the foil) + (−Γ round the starting vortex)".
4. `nb.md` — `> 🔁 **Recap — Kelvin's circulation theorem** (Ch. 5 §5.2).` "$\dfrac{D\Gamma}{Dt}=0$ (5.8): the circulation round a closed
   curve that moves with the fluid does not change, provided the curve itself stays in inviscid fluid of uniform density."
   `> 🔁 **Recap — circulation adds over sub-circuits** (Ch. 3; Stokes' theorem, Ch. 2).` "Cut a circuit in two with a line: the
   two circulations add up to the original, because the cut is walked once each way."
5. `nb.note` — **N14 [B]** "An impulsively started foil at three times: $\Gamma=0\ \to\ 0<\Gamma<\Gamma_{\text{Kutta}}\ \to\ \Gamma=\Gamma_{\text{Kutta}}$ for $t>c/U$
   — the steady pattern is reached after the foil has moved a chord or two." + `nb.figure`: `fig_three_stages(sec["U"], sec["b"], sec["e"], sec["beta"], sec["alpha"])` — three
   streamline panels of our section at Γ/Γ_K = 0, ½, 1 from `AF.zhukhovsky_flow(…, Gamma_cw=f*GK).stream`, both stagnation
   points marked. *see / read / change.*
6. `nb.note` — **N15 [B]** the mechanism: "At the very start the flow does turn the sharp corner. But the fluid next to the
   surface is slow (the boundary layer) and cannot follow so tight a turn against rising pressure: it separates at the edge,
   rolls up, and the shed layer pushes the stagnation point back to the edge. A circulation larger than the Kutta value is
   allowed by ideal flow but is not observed (trap T14)."
7. `nb.derivation("D02", …)` — Part F D02 (5 steps), ref "5.8".
8. `nb.primer("discrete-vortex time stepping", "A cheap model of a wing that sheds vorticity: one 'bound' point vortex stands
   for the wing; at every time step one new vortex is released just behind the trailing edge; its strength is whatever
   keeps the total circulation zero while the air does not flow through the plate at one control point. The total is the
   invariant to watch.", code="G_bound, G_wake = 0.0, []            # nothing circulates at the start\nG_bound = 5.0; G_wake.append(-5.0)
   # step 1: the wing gains 5, sheds -5\nG_bound = 6.5; G_wake.append(-1.5)   # step 2: gains 1.5 more, sheds -1.5\nprint(G_bound +
   sum(G_wake))        # 0.0 at every step")` (**P339**)
9. `nb.worked_example("how strong is the starting vortex?", "A flat plate of chord c = 1 m at α = 30° reaches U = 10 m/s. 1.
   Steady circulation (clockwise): πUc sin α = π × 10 × 1 × 0.5 = 15.7 m²/s (the plate's value of $\\Gamma=4\\pi Ua\\sin(\\alpha+\\beta)$ (14.11),
   with a = c/4 and β = 0 — derived in C06). 2. Starting vortex: −15.7 m²/s, counterclockwise. 3. Lift: $L=\\rho U\\Gamma$ (6.62) =
   1.2 × 10 × 15.7 = 188 N per metre of span. 4. The speed is then doubled to 20 m/s: the wing needs 31.4 m²/s, so it sheds
   another −15.7. 5. The wing stops: the bound 31.4 m²/s is shed as a +31.4 'stopping vortex'. Total at every stage: 0.")`
10. `nb.code` — `run = ch14.impulsive_start_lumped_vortex(sec["U"], 1.0, sec["alpha"], t_end=20/sec["U"], dt=0.05/sec["U"])`; print `run["Gamma_steady"][-1]`, `run["ratio"][[19, 59, 99, 199]]`, `np.abs(run["total"]).max()`.
    *expect:* 8.214 m²/s; 0.646, 0.816, 0.879, 0.939 (after 1, 3, 5, 10 chords); 0 to round-off. *explain:* 4 numbered points
    ("…3. half the lift arrives almost at once, the rest slowly: the starting vortex is still close and its downwash holds the
    circulation down; 4. 'not in the book — ours': a lumped-vortex sketch, not a simulation of the boundary layer").
11. `nb.check_agree` — **from scratch:** the march of C.6 item 4 written as an explicit loop of about fifteen lines (convect,
    place, solve the 2 × 2 by hand, append); `assert np.allclose(G_mine, run["Gamma_bound_cw"])`.
12. `nb.animation` — **A2** (`player="video"`, 80 frames, FAST 40): the plate, the wake vortices (amber dots sized by strength)
    drifting downstream, the bound circulation as a ring; a second panel with bound (teal), shed total (rose) and their sum
    (amber, flat at zero) growing with distance travelled. *see / read / change* ("…α is doubled half-way: a second vortex of
    the same sign as the first is shed").
13. `nb.figure` — the same three curves as a static figure against chords travelled, the steady value dashed (and the curve `ch14.wagner_function(2*s_chords)` as a ghost **only if** the verification report lists its source as read first-hand). *see / read / change.*
14. `nb.note` — **N16 [B]** "**Viscosity is responsible for the drag — and for the lift.** Without the boundary layer at the edge
    no vortex would be shed, no circulation would arise, and the wing would not lift. Once the circulation is set, the rest
    of the flow is nearly ideal." · **N17 [C]** "The idea of circulation as the source of lift goes back to Lanchester, Kutta
    and Zhukhovsky, working independently; Lanchester and Prandtl return in N69 (C08). > ⚠️ **slip #5 — the book prints** a
    birth year for Lanchester that would make him a child when he founded his motor company**; the correct form is** 1868."
15. `nb.explainer("starting_vortex", heading="Where does a wing's circulation come from?", why="A new vortex of the right sign
    rolls off the edge at every change of speed or incidence while the amber total stays flat at zero — something only time
    can show.", tries=["Press play on 'impulsive start' and watch the teal and rose curves mirror each other.", "Choose 'sudden
    stop': which way does the stopping vortex turn?", "Shrink the circuit until it cuts between foil and starting vortex and
    read the circulation it now encloses.", "Step α up half-way through a run and predict the sign of the new vortex first."])`
16. `nb.md` — **What would change if…** "…we want Γ for a real section shape, not a plate? We need a family of shapes whose
    flow we can write down: C04."

### A.4 §14.4 Conformal Transformation for Generating Airfoil Shapes — C04
#### C04 — The Zhukhovsky airfoil: a circle through ζ = b mapped by (14.3)
1. `nb.section("14.4", "Conformal Transformation for Generating Airfoil Shapes", intro="**What is this section about?** A recipe
   for wing sections whose flow we can write down exactly: take a circle, move its centre, and push it through one simple
   map. Where the circle sits decides the camber and the thickness.")`
2. `nb.core("C04", "The Zhukhovsky airfoil: the circle $\\lvert\\zeta-\\mu\\rvert=a$ through $\\zeta=b$, mapped by $z=\\zeta+\\dfrac{b^2}{\\zeta}$ (14.3)",
   question="How can one map turn a circle into a flat plate, a curved plate, a fat symmetric fin or a real-looking wing section?")`
3. `nb.md` — **Plain words:** "We know the flow round a circle (Ch. 6). If we can find a change of coordinates that turns a
   circle into a wing section and leaves the far stream alone, the circle's flow *becomes* the wing's flow, for free. One
   such change is the Zhukhovsky map. Chapter 6 used it on a circle centred at the origin and got a flat plate. The new
   move here is to slide the circle's centre: a little up, a little back, always keeping the circle through the point
   ζ = b. That one point is special — the map folds the plane there — and it becomes the sharp tail."
4. `nb.md` — **The idea** (table): | where the circle's centre is | what comes out | controlled by | — at the origin | flat plate of
   length 4b | — · raised (on the imaginary axis) | circular arc | β (camber) · shifted back (on the negative real axis) |
   symmetric section, round nose, sharp tail | e (thickness) · both | cambered airfoil | β and e.
5. **(re-enter C04)** `nb.recap("R05", "The Zhukhovsky transformation", "$z=\\zeta+\\dfrac{b^2}{\\zeta}$ (14.3) — the same map as
   $z=\\zeta+b^2/\\zeta$ (6.65). Far away z ≈ ζ, so the stream is untouched; near the origin the b²/ζ term dominates. **New here:**
   circles that are not centred on the origin.", where="Ch. 6 §6.6")` + `nb.code`: `CONF.joukowski(2j, 1.0)`,
   `CONF.joukowski_inverse(1.5j, 1.0)` *expect:* 1.5j and 2j.
6. **(re-enter C04)** `nb.recap("R06", "The centred circle becomes a slit", "For ζ = b e^{iθ}: $z=be^{i\\theta}+be^{-i\\theta}=2b\\cos\\theta$ —
   real, running from 2b to −2b and back: a flat plate of length 4b, traced twice (upper and lower side).", where="Ch. 6 §6.6")`
7. `nb.note` — **N18 [C]** "Other routes to the flow round a *given* section: solve Laplace's equation numerically (Ch. 6
   panels, Ch. 10), thin-airfoil theory (a vortex sheet on the camber line; reused in Ch. 15), or Theodorsen's method (map
   the section back to a circle). We do the reverse: choose the map, accept the shape." Optional line:
   `AF.thin_airfoil(lambda xc: 4*0.035*(1 - 2*xc), 0.0)["alpha_L0"]` *expect:* ≈ −0.070 rad = −2 × camber/chord, matching N31.
8. `nb.note` — **N25 [B]** "Shift the centre to the left along the real axis by be and enlarge the radius slightly: $a\equiv b(1+e)$,
   $e\ll1$ (14.7). The circle still passes through ζ = b. **Not in the book — ours:** one parametrisation for all four cases,
   centre $\mu=b-ae^{-i\beta}$ and radius $a=b(1+e)/\cos\beta$." `> ⚠️ Common confusion (trap T7):` "the book defines a = b(1 + e) only
   for the symmetric case and a = b/cos β only for the arc; `AF.zhukhovsky_circle` gives both as limits and also accepts a
   directly." + `nb.code`: `AF.zhukhovsky_circle(1.0, 0.1, 0.0)`, `AF.zhukhovsky_circle(1.0, 0.0, np.arctan(0.75))`,
   `AF.zhukhovsky_circle(sec["b"], sec["e"], sec["beta"])`. *expect:* (−0.1, 1.1); (0.75i, 1.25); (−0.0200 + 0.0189i, 0.2707).
9. `nb.figure` — `fig_four_circles()`: 2 × 4 panels (ζ-plane above, z-plane below) for (e, β) = (0, 0), (0, 8°), (0.1, 0),
   (0.1, 8°) from `AF.zhukhovsky_airfoil`; the ghost circle of radius b, the points ±b marked. *see / read* ("the circle always
   passes through ζ = b, which lands on the tail z = 2b") */ change* ("…the circle misses ζ = b: the tail is rounded and the
   Kutta condition has nothing to hold on to").
10. `nb.primer("the cosine rule", "In any triangle with sides p, q, r and the angle C between p and q: $r^2=p^2+q^2-2pq\\cos C$.
    With C = 90° it is Pythagoras.", code="p, q, C = 3.0, 4.0, np.deg2rad(60)            # two sides and the angle between
    them\nprint(np.sqrt(p**2 + q**2 - 2*p*q*np.cos(C)))   # third side: 3.606")` (**P340**) · `nb.primer("cot, csc and the
    double-angle identities", "cot x = cos x/sin x = 1/tan x; csc x = 1/sin x. Two facts we need: $1+\\cot^2x=\\csc^2x$ and
    $\\cot\\beta-\\tan\\beta=2\\cot2\\beta$.", code="B = 0.3                                      # any angle [rad]\nprint(1/np.tan(B) - np.tan(B),
    2/np.tan(2*B))   # both 2.923\nprint(1 + 1/np.tan(B)**2, 1/np.sin(B)**2)     # both 11.45")` (**P341**) + 🔁 reminders: Euler's
    formula (P45), completing the square into a circle (P272), small angles (P100).
11. `nb.derivation("D03", …)` — Part F D03 (15 steps, ★★★, `check_src`), ref "14.6": includes **N19 [B]** (start; `> ⚠️` trap T4: "this
    θ is the polar angle of ζ from the *origin*, not from the circle's centre"), **N20 [B]** (step 2), **N21 [B]** (steps 3–5),
    **N22 [B]** (steps 6–7, with `draw_triangle_OQP()` shown before the derivation), **N23 [B]** (steps 8–9), **N24 [B]** (steps 10–14)
    and the camber rule of **N26** (step 15). After it: `ch14.circular_arc_sym()["ok"]` → True.
12. `nb.code` — `arc = AF.circular_arc(sec["b"], sec["beta"])`; `S = AF.zhukhovsky_airfoil(sec["b"], 0.0, sec["beta"])`;
    `np.abs(np.abs(S["z"] - arc["centre"]) - arc["radius"]).max()`. *expect:* height 0.0350 m, centre −3.558i m, radius 3.593 m,
    end angle 8°; residual ≤ 1e-13 (every mapped point lies on that circle). *explain:* 3 points ("…the arc is traced twice: the
    lower half of the circle lands on the same arc").
13. `nb.primer("binomial approximation with a complex small quantity", "For a small number ε — real or complex —
    $\\dfrac1{1+\\epsilon}\\approx1-\\epsilon$, with an error of about ε². 'First order in e' means: keep e, drop e².", code="eps =
    0.05*(1 - np.exp(-1j*2.0))            # a small complex number, size 0.08\nprint(1/(1 + eps), 1 - eps)                # agree to
    about 0.007 = abs(eps)**2")` (**P342**)
14. `nb.derivation("D04", …)` — Part F D04 (11 steps, ★★★, `check_src`), ref "14.8".
15. `nb.note` — **N26 [B]** "$c\simeq4b$, camber $\simeq\tfrac12\beta c$, $t_{\max}/c\simeq1.3\,e$ (14.8), with the thickest point near the quarter
    chord x = −b. > ⚠️ **slip #4 — the book prints** 'camber = ≅ ½βc' with two relation signs, and the coefficient 1.3**; the
    correct form is** camber ≃ ½βc, and the coefficient is $3\sqrt3/4=1.299$ (D04) — a first-order rule." + `nb.code`: a table over
    e = 0.02, 0.05, 0.10, 0.25 of `AF.zhukhovsky_shape(1.0, e, 0.0, exact=True)["thickness_ratio"]` against `exact=False`.
    *expect:* measured 0.0255, 0.0618, 0.1179, 0.2572 against the rule 0.0260, 0.0650, 0.1299, 0.3248 — the measured value is
    2 %, 5 %, 9 % and 21 % below the rule; chord/4b = 1.0004, 1.0023, 1.0083, 1.0417; camber rule within 1 % up to β = 6°.
16. `nb.worked_example("a raised circle by hand", "b = 1, tan β = 3/4 (β ≈ 36.9°; far from thin, but the numbers are easy). 1.
    Radius a = b/cos β = 1.25; centre Q at 0.75i. 2. Top of the circle: ζ = 0.75i + 1.25i = 2i. Its image: $z=\\zeta+\\dfrac{b^2}{\\zeta}$ (14.3)
    = 2i + 1/(2i) = 2i − 0.5i = 1.5i. 3. Bottom of the circle: ζ = −0.5i → z = −0.5i + 1/(−0.5i) = −0.5i + 2i = 1.5i — the *same*
    point: the arc is covered twice. 4. Ends: ζ = ±1 → z = ±2. 5. The arc's circle: centre at −2b cot 2β = −0.583, radius 2b csc 2β
    = 2.083; check the top: −0.583 + 2.083 = 1.5 ✓; check an end: √(2² + 0.583²) = 2.083 ✓. 6. Height 2b tan β = 1.5 ✓.")`
17. `nb.code` — `S = AF.zhukhovsky_airfoil(sec["b"], sec["e"], sec["beta"])`; `sh = AF.zhukhovsky_shape(sec["b"], sec["e"],
    sec["beta"])`; `rule = AF.zhukhovsky_shape(sec["b"], sec["e"], sec["beta"], exact=False)`. *expect:* chord 1.0055 m (rule 1.0),
    camber/chord 0.0346 (rule 0.0349), thickness/chord 0.0965 (rule 0.1039), thickest at 25.0 % chord. *explain:* 3 points.
18. `nb.check_agree` — **from scratch:** `t = -beta + np.linspace(0, 2*np.pi, 401)`; `zeta = mu + a*np.exp(1j*t)`; `z = zeta +
    b**2/zeta`; chord = `z.real.max() - z.real.min()`; `assert np.allclose(z, S["z"])` and the chord against `sh["chord"]`
    (rtol 1e-4).
19. **(re-enter C04)** `nb.recap("R07", "Corner flows", "A flow $w=Az^n$ turns a corner; here as the flow leaving a trailing edge
    of finite included angle $2\\pi(1-1/n)$: $\\phi=\\dfrac{Ud}{n}\\Big(\\dfrac rd\\Big)^{n}\\cos n\\theta$, speed $U(r/d)^{n-1}$ → 0 at the edge for n > 1.",
    where="Ch. 6 §6.4")` + `ch14.trailing_edge_wedge(1.05)["included_angle"]` *expect:* 0.299 rad = 17.1°.
20. `nb.derivation("D05", …)` — Part F D05 (7 steps), ref "14.3": **N27 [B]**. + `nb.figure`: a zoom on the last 2 % of our section
    (the two surfaces becoming tangent) beside corner-flow streamlines for n = 1.05 and n = 1 from `POT.power_branch`. *see /
    read / change.*
21. `nb.animation` — **A3** (`player="video"`, 60 frames, FAST 30): the centre μ moves up (β: 0 → 8°), then back (e: 0 → 0.12);
    left panel the circle, right panel the section morphing plate → arc → cambered airfoil. *see / read / change.*
22. `nb.plotly` — **F2** two `slider_figure`s: (a) e from 0 to 0.25 (26 steps): the section, with measured $t_{\max}/c$ and the
    rule $1.299\,e$ in the trace names; (b) β from 0 to 12° (25 steps): the section and its camber line, measured camber/chord
    against β/2. *see / read / change.*
23. `nb.explainer("zhukhovsky_airfoil_builder", heading="Where do camber and thickness come from?", why="Dragging the circle's
    centre and watching the section reshape continuously shows that β and e are two independent knobs, and how far the rules
    drift as the section gets fat.", tries=["Drag the centre straight up: the plate bends but stays thin.", "Drag it left: the
    nose fattens, the tail stays sharp.", "Press 'too fat' and read the gap between the rule and the measured thickness.",
    "Click a point of the circle and follow its arithmetic to the section."])`
24. `nb.md` — **What would change if…** "…a stream now blows past this shape? We already know the flow past the circle: C05
    carries it over."

### A.5 §14.5 Lift of a Zhukhovsky Airfoil — C05, C06, C07
#### C05 — The flow round the airfoil: the velocity in the two planes (14.9)
1. `nb.section("14.5", "Lift of a Zhukhovsky Airfoil", intro="**What is this section about?** Carry the circle's flow to the
   airfoil, find the one circulation that keeps the trailing-edge speed finite, and read off the lift — then look at the
   pressure that produces it.")`
2. `nb.core("C05", "The flow round the airfoil: $\\dfrac{dw}{dz}=\\dfrac{dw}{d\\zeta}\\,\\dfrac{\\zeta^2}{\\zeta^2-b^2}$ (14.9)", question="If we know the
   flow round the circle, what is the flow round the airfoil?")`
3. `nb.md` — **Plain words:** "The map is a change of address: every point ζ outside the circle has a partner z outside the
   airfoil. The complex potential — the pair (potential, stream function) — keeps its *value* at partner points, so
   streamlines map to streamlines. The *speed* does not keep its value: where the map stretches the plane the fluid slows
   down, where it squeezes the plane the fluid speeds up. At the tail the map squeezes infinitely hard."
4. `nb.md` — **The idea** (ASCII): `ζ-plane: circle + flow  ──(same w at partner points)──►  z-plane: airfoil + flow` and under it
   `speed_z = speed_ζ × |dζ/dz|   ← 1 far away, ∞ at ζ = b`.
5. `nb.derivation("D06", …)` — Part F D06 (7 steps), ref "6.68": **N28 [B]** the ζ-plane flow ("not in the book — ours: the book
   writes only the centred case"). + `nb.code`: `w, dw = CONF.circle_flow_zeta(30.0, a, Gamma_cw=10.0, alpha=sec["alpha"],
   center=mu)`; the far-field value `dw(1e6)`; `np.ptp(np.imag(w(mu + a*np.exp(1j*t))))` on the circle. *expect:* 29.886 − 2.615i
   (= 30e^{−iα}); ≤ 1e-12 (the circle is a streamline).
6. **(re-enter C05)** `nb.recap("R08", "The chain rule for the complex velocity", "$\\dfrac{dw}{dz}=\\dfrac{dw}{d\\zeta}\\dfrac{d\\zeta}{dz}$: the velocity
   in the z-plane is the velocity in the ζ-plane times the map's local rescaling.", where="Ch. 6 §6.6")` `> ⚠️ Common confusion
   (trap T2):` "here $dw/dz=u-iv$ — note the minus sign. Chapter 13's Ekman velocity was $V=u+iv$."
7. `nb.derivation("D07", …)` — Part F D07 (4 steps), ref "14.9".
8. `nb.derivation("D08", …)` — Part F D08 (5 steps), ref "6.62": the circulation is the same in both planes. + `nb.code` (trap
   **T6**): `F = AF.zhukhovsky_flow(30.0, sec["b"], sec["e"], sec["beta"], sec["alpha"], Gamma_cw=10.0)`; for three circles of
   radius 1.5a, 3a, 10a about μ: the periodic trapezoid sum of `F.dwdzeta(zeta)*dzeta` in the ζ-plane and of
   `F.dwdz_from_zeta(zeta)*dz` with `dz = CONF.joukowski_derivative(zeta, b)*dzeta` in the z-plane. *expect:* real part −10.000
   six times (counterclockwise −Γ); and `abs(F.dwdz_from_zeta(zeta))/abs(F.dwdzeta(zeta))` at three points: 1.000 far away, not 1
   near the body ("the velocity is **not** the same, only its loop integral").
9. `nb.worked_example("three points by hand", "Flat plate (a = b = 1), U = 1. 1. α = 0, Γ = 0, ζ = 2: circle velocity
   $1-a^2/\\zeta^2$ = 0.75; map factor $\\zeta^2/(\\zeta^2-b^2)$ = 4/3; product 1.00 — a plate edge-on does not disturb the stream ✓. 2. The
   same at ζ = 2i: (1 + 0.25) × (−4)/(−5) = 1.25 × 0.8 = 1.00 ✓. 3. Now α = 30° and ζ = 1.1, just off the trailing edge: the
   circle velocity is $e^{-i\\alpha}-e^{i\\alpha}/1.21$ = 0.150 − 0.913i (size 0.926); the factor is 1.21/0.21 = 5.76; the speed is
   5.33 U — and it grows without limit as ζ → 1. 4. Add the Kutta circulation Γ = 2π (C02's hand example): the extra term
   $i\\Gamma/(2\\pi\\zeta)$ = +0.909i leaves 0.150 − 0.004i, and the speed is 0.150 × 5.76 = 0.87 U: finite.")`
10. `nb.code` — `F = AF.zhukhovsky_flow(sec["U"], sec["b"], sec["e"], sec["beta"], sec["alpha"])` (Kutta value by default);
    `F.velocity(0.0, 0.3)`, `F.velocity(-5.0, 0.0)`, `F.speed(2*sec["b"] + 0.01)`, `F.inside(0.0, 0.0)`. *expect:* far point ≈
    (29.9, 2.6) m/s; inside → True (NaN velocity). *explain:* 4 points ("…3. off the surface the function needs ζ for a given
    z: the inverse map has two roots and the one outside the *shifted* circle is taken — trap T8").
11. `nb.check_agree` — **from scratch:** at `zeta0 = mu + 1.3*a*np.exp(0.7j)`: `dw = U*(np.exp(-1j*al) - a**2*np.exp(1j*al)/(zeta0
    - mu)**2) + 1j*F.Gamma_cw/(2*np.pi*(zeta0 - mu))`; `q = dw*zeta0**2/(zeta0**2 - b**2)`; `assert np.isclose(q,
    F.dwdz(zeta0 + b**2/zeta0))`.
12. `nb.figure` — two panels with the same 14 streamlines: the circle flow (ζ-plane) and the airfoil flow (z-plane) from
    `F.stream`, the stagnation points of `AF.stagnation_points` marked in both. *see / read* ("partner points carry the same
    stream-function value") */ change* ("…Γ is not the Kutta value: C06").
13. `nb.plotly` — **F3** `slider_figure`: Γ/Γ_K from 0.5 to 1.5 (21 steps) → speed/U along the upper surface over the last 5 % of
    chord, log y-axis, from `AF.zhukhovsky_surface(…, Gamma_cw=f*GK)`. *see / read* ("only at 1.00 does the curve stay finite at
    the edge") */ change.*
14. `nb.md` — **What would change if…** "…we *demand* a finite speed at the tail? That fixes Γ: C06."

#### C06 — The lift of a Zhukhovsky airfoil: Kutta circulation (14.11) and lift coefficient (14.12)
1. `nb.core("C06", "The lift of a Zhukhovsky airfoil: $\\Gamma=4\\pi Ua\\sin(\\alpha+\\beta)$ (14.11) ⇒ $C_L\\simeq2\\pi(\\alpha+\\beta)$ (14.12)",
   question="How much lift does a section give — and what do angle of attack and camber each contribute?")`
2. `nb.md` — **Plain words:** "A pilot pulls the nose up two degrees and the lift rises by a fifth. A flap comes down and the
   lift rises without the nose moving at all. Both facts are in one line of geometry: the lift is proportional to the
   angle between the stream and one particular direction of the circle — the line from its centre to the point that
   becomes the tail. Tilting the wing changes that angle by α; cambering the wing changes it by β."
3. `nb.md` — **The idea** (ASCII of the circle): centre Q, the point B (ζ = b) lying β below the horizontal through Q, the stream
   arriving from below-left at α; "angle between the stream direction and QB = α + β — the *absolute* angle of attack".
4. **(re-enter C06)** `nb.recap("R09", "Surface speed on a circle with circulation", "$u_\\theta=-2U\\sin\\theta-\\dfrac{\\Gamma}{2\\pi a}$ (14.10), the same
   formula as $u_\\theta=-2U\\sin\\theta-\\dfrac{\\Gamma}{2\\pi a}$ (6.37): the first term is the stream sliding round the circle, the second the
   clockwise swirl. **New here:** θ is measured from the diameter through the centre Q that is aligned with the stream.",
   where="Ch. 6 §6.3")` + `draw_circle_angles()` + `> ⚠️ Common confusion (trap T4):` "three different θ's: from the origin in
   $x=(R+b^2/R)\cos\theta$ (14.4); from the stream-aligned diameter through Q here; a course angle in §14.9." + `nb.code`:
   `ch06.cylinder_surface_speed(-(sec["alpha"] + sec["beta"]), sec["U"], a, Gamma_cw=GK)` *expect:* 0.0 (to round-off).
5. `nb.primer("radians and degrees for angles and for slopes", "Formulas want radians: 1° = π/180 = 0.01745 rad. A slope 'per
   radian' is 57.3 times a slope 'per degree': 2π per radian is 0.110 per degree.", code="print(np.deg2rad(9.0))                #
   9 degrees = 0.1571 rad\nprint(2*np.pi, 2*np.pi*np.pi/180)      # lift slope: 6.283 per radian = 0.1097 per degree")` (**P338**)
6. `nb.derivation("D09", …)` — Part F D09 (10 steps), ref "14.12": includes **N30 [B]** $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11), boxed at
   step 6, and slip **#8** at step 9.
7. `nb.note` — **N34 [B]** the flat plate (e = β = 0, a = b = c/4): "$\Gamma=4\pi Ub\sin\alpha$ and $C_L=2\pi\sin\alpha$ **exactly** — no thin
   or small-angle step; this is the flow of C02." `AF.flat_plate_flow(30.0, 0.25, sec["alpha"]).Gamma_kutta` *expect:* 8.214 m²/s
   = C02's root.
8. `nb.primer("l'Hôpital's rule for a 0/0 limit", "If top and bottom of a fraction both go to zero at a point, the limit of
   the fraction is the limit of (slope of the top)/(slope of the bottom).", code="x = sp.symbols('x')\nprint(sp.limit(sp.sin(3*x)/sp.sin(x),
   x, 0))        # 3: slopes 3cos(0)/cos(0)\nprint(sp.limit((x**2 - 1)/(x - 1), x, 1))        # 2")` (**P343**)
9. `nb.note` — **N29 [B]** "At the trailing edge $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9) is 0/0 when the Kutta condition holds. **Not
   in the book — ours:** the limit is finite, $\Big\lvert\dfrac{dw}{dz}\Big\rvert_{B'}=U\,\dfrac ba\cos(\alpha+\beta)$ — stated, not derived; verified two ways."
   + `nb.code`: (i) a sympy `limit` of that expression along ζ = b + δ for the flat plate at α = 30° → cos 30°; (ii) a table of
   `ch14.edge_speed(U, b, e, beta, alpha, d, Gamma_cw=f*GK)/U` for d = 10⁻², 10⁻³, 10⁻⁴, 10⁻⁵ m and f = 1, 0.9, 0.5.
   *expect:* `AF.trailing_edge_speed` = 27.37 m/s = 0.912 U; the f = 1 column tends to 0.912; the others grow like
   $d^{-1/2}$ (trap **T8**).
10. `nb.note` — **N31 [B]** "α + β is the **absolute angle of attack**; the lift vanishes at $\alpha_{L=0}=-\beta$, and for these sections
    $\beta=2\,\dfrac{\text{camber}}{c}$. More camber slides the whole lift line to the left. > ⚠️ Common confusion (trap T3): β is the
    circle's camber angle here; in §14.6 the same letter is *defined* as minus the zero-lift angle of whatever section the wing
    has." `AF.zero_lift_angle(0.035)` *expect:* −0.070 rad = −4.0°.
11. `nb.worked_example("lift by hand", "U = 10 m/s, a = 0.25 m (chord ≈ 1 m), ρ = 1.2 kg/m³, and α + β = 30° so that the sine is
    easy. 1. $\\Gamma=4\\pi Ua\\sin(\\alpha+\\beta)$ (14.11) = 4π × 10 × 0.25 × 0.5 = 15.7 m²/s. 2. $L=\\rho U\\Gamma$ (6.62) = 1.2 × 10 × 15.7 = 188 N
    per metre. 3. $C_L=L/(\\tfrac12\\rho U^2c)$ = 188/(0.5 × 1.2 × 100 × 1) = 3.14 (= 2π sin 30°). 4. The small-angle form
    $C_L\\simeq2\\pi(\\alpha+\\beta)$ (14.12) gives 2π × 0.524 = 3.29 — 5 % high at 30°, which is far beyond where a real section has
    stalled. 5. A realistic case: α = 5°, β = 4°: α + β = 0.157 rad, $C_L$ ≈ 0.99, and sin and angle differ by 0.4 %.")`
12. `nb.code` — `GK = AF.kutta_circulation(sec["U"], a, sec["alpha"], sec["beta"])`; `ch06.lift_per_span(sec["rho"], sec["U"],
    Gamma_cw=GK)`; `AF.zhukhovsky_lift_coefficient(sec["alpha"], sec["beta"])`; `AF.zhukhovsky_lift_coefficient(sec["alpha"],
    sec["beta"], sec["e"], exact=True)`; `AF.kutta_circulation_root(sec["U"], sec["b"], sec["e"], sec["beta"], sec["alpha"])`.
    *expect:* 15.962 m²/s; 586.6 N/m; 0.987; 1.058; 15.962 (the root-finder of C02 agrees to round-off). *explain:* 4 points
    ("…4. the exact value is 7.2 % above 2π(α + β): thickness makes a larger than c/4 — slip #8").
13. `nb.check_agree` — **from scratch (three lines):** `G = 4*np.pi*U*a*np.sin(alpha + beta)`; `c = sh["chord"]`; `CL = 2*G/(U*c)`;
    `assert np.isclose(CL, AF.zhukhovsky_lift_coefficient(alpha, beta, e, exact=True), rtol=1e-6)`.
14. `nb.plotly` — **F4** `slider_figure`: β from 0 to 8° (17 steps) → $C_L$ against α in degrees (−10…15): the line
    $2\pi(\alpha+\beta)$ (teal), the exact curve at e = 0.08 (dashed), the zero-lift angle −β marked; trace names carry "slope 0.110
    per degree" and the exact slope. *see / read* ("the slider moves the intercept, never the slope") */ change.*
15. `nb.note` — **N32 [B]** theory against measurement: `nb.figure` with $C_L\simeq2\pi(\alpha+\beta)$ (14.12), the exact curve and
    `ch14.generic_lift_curve(alpha, K=5.7, beta=sec["beta"], alpha_stall=np.deg2rad(14), CL_post=0.9, width=np.deg2rad(2))` on one
    axis, title "Theory (lines) and a **qualitative model** of a measured curve — constants ours". "Measured slopes are a
    little below 2π (the boundary layers' displacement thickness changes the effective shape) and the measured curve ends in
    stall." *see / read / change.*
16. `nb.live` — α, β, e sliders → section + streamlines + $C_L$ (three numbers); paired with F4.
17. `nb.explainer("kutta_condition", heading="Which circulation does a wing choose?", why="Sliding Γ and watching the rear
    stagnation point crawl to the edge while the edge-speed readout falls from unbounded to a finite number — and rises again
    past the Kutta value — is the hypothesis as an experiment.", tries=["Drag Γ/Γ_K from 0 to 1 and watch the amber dot reach
    the tail.", "Press 'too much': where is the rear stagnation point now?", "Set α = −β: what circulation does the edge ask
    for?", "Open the Derivation tab at step 3 and find the angle α + β in the circle view."])`
18. `nb.explainer("lift_curve_stall", heading="Does more angle always mean more lift?", why="Dragging α along the curve while the
    section view goes from attached to separated, with camber shifting the line and the regime chips moving the stall point,
    lets you find slope, shift and limit by experiment.", tries=["Raise the camber: does the slope change?", "Find the zero-lift
    angle and compare it with −β.", "Put the flap down and read the new stall speed.", "Go just past stall: how far below the
    line is the model?"])`
19. `nb.note` — **N33 [C]** "Zhukhovsky sections are not used on aircraft: a cusp cannot be built, and the camber peaks at
    mid-chord. The Kármán–Trefftz map (a finite edge angle) is the standard relative; this book does not develop it." ·
    **N35 [C]** "Plotting a cambered section with its streamlines, and moving a known lift coefficient to a second incidence by
    $C_L(\alpha_2)=C_L(\alpha_1)+2\pi(\alpha_2-\alpha_1)$: `scripts/ch14_zhukhovsky_airfoil.py` does both with our parameters."
20. `nb.md` — **What would change if…** "…we ask *where on the surface* this lift is produced? C07 looks at the pressure."

#### C07 — Surface speed and pressure on the airfoil (ours: the book's pressure figure, computed)
1. `nb.core("C07", "Surface speed and pressure: $C_p=1-\\dfrac{\\lvert dw/dz\\rvert^2}{U^2}$ on the contour (ours)", question="Is a wing
   pushed up from below or pulled up from above — and where?")`
2. `nb.md` — **Plain words:** "Put your hand out of a car window, palm tilted up. It feels pushed from below. A real wing is
   mostly *pulled* from above: the air speeds up over the curved upper surface, its pressure drops, and that suction —
   strongest just behind the nose — does most of the lifting. The same picture shows the wing's weak spot: behind the suction
   peak the pressure has to climb back, and a boundary layer made to climb too steeply gives up. That is stall. **Not in the
   book — ours:** the book shows this pressure figure; we compute it."
3. `nb.md` — **The idea** (ASCII): a section with short arrows pointing *away* from the upper surface near the nose (suction),
   small arrows pushing on the lower surface; beside it the two curves of $-C_p$ against x/c enclosing a loop, labelled "area
   = normal force".
4. **(re-enter C07)** `nb.recap("R03", "The pressure coefficient", "$C_p=\\dfrac{p-p_\\infty}{\\tfrac12\\rho U^2}$: +1 at a stagnation point, 0 where
   the speed equals U, negative where the flow is faster. **New here:** aerodynamicists plot $-C_p$ upward, so suction is 'up',
   and the two surfaces make a loop.", where="Ch. 6 §6.3")` + `> 🔁 **Recap — Bernoulli** (Ch. 4 §4.9):` "along a streamline of steady
   ideal flow $p+\tfrac12\rho q^2$ is constant, so $C_p=1-q^2/U^2$." `POT.pressure_coefficient(45.0, 30.0)` *expect:* −1.25.
5. `nb.md` — **The maths, step by step.** (i) points of the circle $\zeta=\mu+ae^{it}$; (ii) circle velocity $dw/d\zeta$ there (D06) with
   the Kutta circulation of $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11); (iii) rescale by the map, $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9); (iv)
   speed $q=\lvert dw/dz\rvert$; (v) $C_p=1-q^2/U^2$; (vi) at t = −β the 0/0 is replaced by its limit (N29, trap T8).
6. `nb.derivation("D10", …)` — Part F D10 (7 steps), ref "14.12": **N09 [B]** $F_y=-\int_0^cp_u\,dx+\int_0^cp_l\,dx$,
   $C_y\equiv\dfrac{F_y}{\tfrac12\rho U^2c}=\oint C_p\,d\!\Big(\dfrac xc\Big)$.
7. `nb.worked_example("pressure by hand", "U = 30 m/s, ρ = 1.225 kg/m³: dynamic pressure 551 Pa. 1. A point where the speed is
   45 m/s: $C_p$ = 1 − (45/30)² = −1.25, i.e. 689 Pa *below* ambient. 2. The stagnation point: $C_p$ = 1, 551 Pa above. 3. A
   toy loop: $C_p$ = −0.6 all along the top, +0.4 all along the bottom: $C_y$ = ∫(0.4 − (−0.6)) d(x/c) = 1.0 — and 60 % of it
   is suction. 4. Force normal to a 1 m chord: 1.0 × 551 × 1 = 551 N per metre of span.")`
8. `nb.code` — `S = AF.zhukhovsky_surface(sec["U"], sec["b"], sec["e"], sec["beta"], sec["alpha"], n=2000)`; `S["cp_min"],
   S["xc_cp_min"]`; `Cy = AF.normal_force_coefficient(S["upper"]["xc"], S["upper"]["cp"], S["lower"]["xc"], S["lower"]["cp"])`;
   `AF.lift_from_cp(…, sec["alpha"])`; `AF.lift_from_cp(…, sec["alpha"], suction="none")`; `2*S["Gamma_cw"]/(sec["U"]*S["chord_x"])`;
   `AF.surface_circulation(S)`. *expect:* suction peak −2.24 at 0.45 % chord on the upper surface; $C_y$ = 1.054; $C_L$ = 1.058
   (ideal) and 1.050 (no leading-edge suction); 2Γ/(Uc) = 1.058; ∮u_s ds = 15.96 m²/s = Γ_K (**N10**'s number). *explain:* 5
   points ("…4. the loop gives the force normal to the *chord*; the lift is normal to the *stream* — trap T5").
9. `nb.check_agree` — **from scratch:** `t = -beta + np.linspace(0, 2*np.pi, 2001)[1:-1]`; `zeta = mu + a*np.exp(1j*t)`; `dw =
   U*(np.exp(-1j*al) - a**2*np.exp(1j*al)/(zeta - mu)**2) + 1j*GK/(2*np.pi*(zeta - mu))`; `q = np.abs(dw*zeta**2/(zeta**2 -
   b**2))`; `cp = 1 - (q/U)**2`; `assert np.allclose(cp, S["cp"][1:-1], rtol=1e-9)` (same n).
10. `nb.figure` — `fig_pressure_two_ways(S)`: (left) the section with pressure arrows normal to the surface, outward where
    $C_p<0$ (orange), inward where $C_p>0$ (blue), the stagnation point marked; (right) $-C_p$ against x/c for both surfaces,
    the loop shaded, its area in the legend. *see / read* ("the tall thin spike at the nose and the long slope behind it") */
    change* ("…α = 10°: the peak passes −7 and the recovery steepens").
11. `nb.plotly` — **F5** `slider_figure`: α from −4° to 12° (17 steps) → $-C_p$ against x/c, both surfaces, legend "loop area
    $C_y$ = …; 2Γ/(Uc) = …". *expect:* at α = −4° the loop has zero net area; at 0° 0.472; at 10° 1.612 (against 1.637).
12. **(re-enter C07)** `nb.recap("R10", "The Blasius force", "$D-iL=\\dfrac{i\\rho}{2}\\oint\\Big(\\dfrac{dw}{dz}\\Big)^2dz$ (6.60): the force from a contour
    integral of the squared complex velocity. **New here:** with the stream inclined, the force comes out perpendicular to the
    *stream* — size ρUΓ, no drag.", where="Ch. 6 §6.5")` + `nb.code`: `AF.zhukhovsky_force(…, rho=1.225, method="blasius")` and
    `method="pressure"`. *expect:* both: Fx = −51.1, Fy = 584.4 N/m, L = 586.6 N/m, D = 0.0 (pressure route: 1e-4), centre of
    pressure at 35.7 % chord. "Two routes that never share a line of code agree."
13. `nb.note` — **N12 [B]** "The pressure on a flat plate is normal to the plate, so why is the force normal to the stream? Because
    the flow whipping round the sharp *leading* edge pulls the plate forward: the leading-edge suction $S=L\sin\alpha=\pi\rho U^2c\sin^2\alpha$,
    with $N=L\cos\alpha$ (ours; stated, not derived). A rounded nose spreads that suction over a real surface — the −2.24 peak
    above — and keeps the boundary layer attached." `AF.flat_plate_forces(1.225, 30.0, 1.0, sec["alpha"])` *expect:* L = 301.9,
    N = 300.7, S = 26.3 N/m; a check line `ch14.edge_flow_suction` ≈ S.
14. `nb.note` — **N13 [B]** stall: "At high α the suction-side boundary layer cannot climb the pressure rise behind the peak; it
    separates ahead of the trailing edge (where the wall shear stress $\tau_0=0$), the Kutta condition fails and the lift
    collapses." 🔁 "Thwaites' method (Ch. 9 §9.6) marches a laminar boundary layer along a given outer speed and reports where
    $\tau_0$ reaches zero." + `nb.figure` (title "**Qualitative model**: laminar separation point on the ideal pressure
    distribution"): `ch14.suction_side_separation(alpha, sec["e"], sec["beta"], 1e6)["x_sep_over_c"]` for 12 angles from 0 to 11°.
    *see / read / change.* `> ⚠️` trap T16: "this cannot predict the maximum lift — the real layer turns turbulent, and once it
    separates the outer flow changes."
15. `nb.explainer("airfoil_pressure_lift", heading="Is a wing pushed up or pulled up?", why="Turning α and seeing the suction
    peak sharpen, the stagnation point slide under the nose and the loop area grow in step with the lift readout connects
    three pictures of one flow; clicking a point traces Bernoulli's arithmetic there.", tries=["Raise α from 0° to 8° and watch
    the peak.", "Click the stagnation point: what is $C_p$ there?", "Compare the two bars 'upper suction' and 'lower
    pressure'.", "Choose 'flat plate': what happens to the peak?"])`
16. `nb.md` — **What would change if…** "…the wing has ends? Everything so far was per metre of an endless span. §14.6."

### A.6 §14.6 Elementary Lifting Line Theory for Wings of Finite Span — C08 … C14
#### C08 — Bound vortex and trailing filaments: what leaves the wing where the circulation drops
1. `nb.section("14.6", "Elementary Lifting Line Theory for Wings of Finite Span", intro="**What is this section about?** A real
   wing has tips. The circulation must fall to zero there, a vortex cannot simply end, so vorticity streams off the wing and
   trails behind. That wake makes the air at the wing sink, the sinking tilts the lift backward, and the wing has a drag
   even in ideal flow. Seven steps: the trailing filaments, their downwash, the induced drag, a sine series, one linear
   system for any wing, the forces from its coefficients, and the best loading.")` + `nb.md`: `draw_axes_conventions()` again,
   right panel only — "**From here on y runs along the span (−s/2 port tip … +s/2 starboard tip) and the downwash w is positive
   downward**" (trap T5).
2. `nb.core("C08", "Bound vortex and trailing filaments: $d\\Gamma_{\\text{trail}}=-\\dfrac{d\\Gamma}{dy}\\,dy$", question="The circulation round
   a wing must vanish at the tips. Where does it go?")`
3. `nb.md` — **Plain words:** "On a damp day you can see two white cords streaming back from an airliner's wing tips or flap
   ends: tip vortices. Below a lifting wing the pressure is high, above it low; at the tips the air leaks round from below
   to above. That leak is the visible end of something that happens all along the trailing edge: wherever the wing carries
   a little less circulation than it did a step further inboard, the difference leaves as a thin vortex trailing downstream."
4. `nb.md` — **The idea** (ASCII, top view): the lifting line `════════` with Γ(y) drawn as a hump above it; under it a staircase
   of nested horseshoes `╔══╗ ⊂ ╔════╗ ⊂ ╔══════╗`, each with two legs trailing back; "many weak horseshoes, one inside the next
   = a bound vortex of varying strength + a sheet of trailing filaments".
5. `nb.note` — **N36 [B]** "Air flows outward under the wing and inward over it. Behind the trailing edge the two streams meet
   with different spanwise velocities: a surface across which the tangential velocity jumps — a vortex sheet — of opposite
   sense on the two half-wings." 🔁 `> 🔁 **Recap — vortex sheet** (Ch. 5 §5.8):` "a jump in tangential velocity across a surface is
   a sheet of vorticity whose strength per unit length equals the jump." + `nb.figure`: `draw_tip_flow_rear_view()`. *see /
   read / change.*
6. **(re-enter C08)** `nb.recap("R11", "A vortex cannot end in the fluid", "Helmholtz: the strength of a vortex tube is the
   same all along it, so a tube must close on itself or end on a boundary. **New here:** the wing's vortex system is one
   closed loop — the bound vortex on the wing, two trailing vortices, and far behind the starting vortex of C03.",
   where="Ch. 5 §5.3")`
7. `nb.note` — **N39 [B]** the three assumptions, as a numbered box that later blocks cite: "(i) the aspect ratio is large, so
   each spanwise strip sees a two-dimensional flow — rules out stubby and delta wings; (ii) the wing is one straight bound
   vortex of strength Γ(y), the *lifting line* — rules out sweep; (iii) the trailing filaments are straight, parallel to the
   stream, and do not act on one another — rules out the roll-up we are about to watch."
8. `nb.note` — **N40 [B]** "The bound vortex moves with the wing, is strongest at mid-span and vanishes at the tips,
   $\Gamma(\pm s/2)=0$." `LL.circulation_series(np.array([12.0, 0.0, 1.0]), np.array([0.0, np.pi/2, np.pi]))` *expect:* 0, 11, 0.
9. `nb.derivation("D11", …)` — Part F D11 (4 steps), ref "14.13".
10. `nb.worked_example("a staircase of circulation", "On the starboard half the bound circulation is 10 m²/s at the root, then 8,
    then 5, then 0 at the tip. 1. Where it drops from 10 to 8 a filament of strength 2 leaves. 2. From 8 to 5: strength 3. 3.
    From 5 to 0 at the tip: strength 5 — the strongest. 4. Total shed on this half: 2 + 3 + 5 = 10 = the root circulation
    $\\Gamma_0$. 5. The port half sheds the same amounts with the opposite sense of rotation. The trailing vorticity is
    concentrated near the tips because that is where Γ falls fastest.")`
11. `nb.code` — `g = np.linspace(0.05, np.pi - 0.05, 9)`; `LL.trailing_sheet_strength(np.array([14.17]), g, ac["s"])`;
    `ch14.tip_vortex_estimate(ac["W"], 1.225, ac["U_cruise"], ac["s"])`. *expect:* the sheet strength is antisymmetric, zero at
    mid-span and largest toward the tips; Γ₀ = 14.17 m²/s, spacing πs/4 = 9.42 m, descent 0.24 m/s. *explain:* 3 points.
12. `nb.check_agree` — **from scratch:** `y = LL.span_position(g, s)`; `G = 14.17*np.sin(g)`; a centred finite difference
    `-(G[2:] - G[:-2])/(y[2:] - y[:-2])` against `LL.trailing_sheet_strength(...)[1:-1]` (rtol 2e-2 on the coarse grid, 1e-5 on 2001
    points).
13. `nb.figure` — `draw_horseshoe_staircase(7)` beside Γ(y) (amber) and −dΓ/dy (blue) for the elliptic loading. *see / read /
    change* ("…the loading is flat-topped: the filaments crowd even closer to the tips").
14. `nb.note` — **N37 [B]** "Assumption (iii) is not what happens: each half of the sheet rolls up into one **tip vortex** of
    circulation $\Gamma_{\text{tip}}=\Gamma_0$. **Not in the book — ours:** the two are a distance $b'=\dfrac1{\Gamma_0}\displaystyle\int\Gamma\,dy$ apart ($\tfrac\pi4s$ for
    the elliptic loading) and sink at $\Gamma_0/2\pi b'$." 🔁 "Point vortices move with the velocity the others induce (Ch. 5 §5.6)."
    + `nb.animation` **A4** (`player="frames"`, 4 stop-frames, title "**Qualitative** past the first turn"): `ch14.wake_rollup(Gn, s,
    n=200 (FAST 80), t_end=…, eps=…)` cached in `outputs/ch14/cache/`; rear view, vortices coloured by sign, the two centroids
    marked; readout of the invariants. Hook: contrails and wake-turbulence separation at airports. *see / read / change.*
15. `nb.note` — **N69 [C]** "The picture of a bound vortex with trailing vortices is Lanchester's; the usable theory is
    Prandtl's (see N17 in C03)."
16. `nb.md` — **What would change if…** "…we stand on the wing and ask what all those filaments do to the air there? C09."

#### C09 — Downwash: the downward velocity the trailing sheet induces at the wing (14.13)
1. `nb.core("C09", "Downwash: $w(y_1)=\\dfrac{1}{4\\pi}\\displaystyle\\int_{-s/2}^{s/2}\\dfrac{d\\Gamma}{dy}\\dfrac{dy}{(y_1-y)}$ (14.13)", question="How fast
   is the air sinking where the wing meets it?")`
2. `nb.md` — **Plain words:** "A wing holds the aircraft up by throwing air down — Newton's third law. Some of that downward
   motion is already present *at the wing itself*, caused by the trailing vortices behind it. Every filament swirls the air
   round its own axis; between the tips all of them push down. Add up their pushes at one point of the wing and you have the
   downwash there."
3. `nb.md` — **The idea** (ASCII, top view): one filament leaving at y, the station y₁ on the lifting line at distance
   abs(y − y₁); "an infinite line vortex would give Γ/(2π d); this one starts at the wing, so the wing sees half: Γ/(4π d)".
4. `nb.md` — `> 🔁 **Recap — the filament law** (Ch. 5 §5.5).` "$d\mathbf u=\dfrac{\Gamma\,d\mathbf l}{4\pi}\times\dfrac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ (5.17): every piece of a
   vortex filament induces a velocity round itself; added along a straight piece it gives the speed
   $\dfrac{\Gamma}{4\pi d}(\cos\theta_a-\cos\theta_b)$ at distance d from the line."
5. **(re-enter C09)** `nb.recap("R12", "Half an infinite line vortex", "For a filament that starts at the lifting line and runs
   downstream for ever the two angles are 0 and 90°: speed Γ/(4πd), half the value Γ/(2πd) of a line infinite both ways. With
   the strength of D11: $dw(y_1)=\\dfrac{-(d\\Gamma/dy)\\,dy}{4\\pi(y-y_1)}$.", where="Ch. 5 §5.5")` + `nb.code`: `BS.segment_speed(1.0, 0.0,
   np.pi/2, 4*np.pi)`, `BS.segment_speed(1.0, 0.0, np.pi, 4*np.pi)` *expect:* 1.0 and 2.0 m/s. `> ⚠️ **slip #3 — the book prints**
   that this velocity is 'at location y (< y₁)'**; the correct form is** it is induced *at y₁*, and it is downward for a
   filament at y > y₁ on the starboard half, as the book's own figure caption says.`
6. `nb.primer("Cauchy principal value", "An integrand like 1/(x − c) is +∞ on one side of c and −∞ on the other. The two
   infinite areas cancel if we cut out a *symmetric* gap (c − ε, c + ε) and let ε → 0: that limit is the principal value.",
   code="val = quad(lambda x: 1.0, -1, 1, weight='cauchy', wvar=0.3)[0]   # PV of the integral of 1/(x - 0.3) from -1 to
   1\nprint(val, np.log(0.7/1.3))                                    # both -0.619")` (**P344**) + 🔁 "an improper integral as a limit
   (P144, Ch. 5)."
7. `nb.derivation("D12", …)` — Part F D12 (7 steps), ref "14.13" (trap **T10** at step 7; trap **T5** at step 3).
8. `nb.worked_example("the simplest wing: one horseshoe", "**N60 [B]** One bound vortex of constant strength Γ = 4π ≈ 12.6 m²/s
   and span s = 2 m, with one trailing leg from each tip. 1. Each leg is a half-line at distance s/2 = 1 m from the centre:
   it induces Γ/(4π × 1) = 1 m/s there. 2. Both legs push down at the centre: w(0) = 2 m/s. 3. In symbols $w(0)=\\dfrac{\\Gamma}{\\pi s}$
   (ours) = 4π/(2π) = 2 ✓. 4. Move the station toward a tip, to y₁ = 0.9 m: the near leg is 0.1 m away and gives 10 m/s, the
   far leg 1.9 m away gives 0.53 m/s — the downwash grows without limit at the tips.")`
9. `nb.code` — `LL.horseshoe_downwash(0.0, 0.0, 10.0, 12.0)`, `LL.horseshoe_velocity(0.0, 0.0, 0.0, 10.0, 12.0)`,
   `LL.horseshoe_wing(10.0, 12.0, 1.225, 30.0)`; then the smooth wing: `LL.downwash(0.0, lambda y: np.sqrt(1 - (2*y/12)**2), 12.0)`
   and at y₁ = 3 m. *expect:* 0.2653 m/s; (0, 0, −0.2653) (z up, so a downwash is negative there); "ratio" 0.5 with the note
   that the centre downwash was used for the whole span; elliptic loading with Γ₁ = 1: 0.04167 = Γ₁/2s at both stations.
   *explain:* 4 points. `> ⚠️ Common confusion (trap T10):` "the horseshoe's induced drag, $C_{D_i}=C_L^2/(2\pi\Lambda)$ (ours), looks half
   that of the best possible wing. It is an artefact: it uses the centre downwash everywhere, while the true downwash of a
   constant Γ is infinite at the tips."
10. `nb.check_agree` — **from scratch:** `ye, Gs = LL.staircase_from_loading(Gfn, 12.0, 21)`; a loop over the 22 trailing
    filaments, each adding `dG/(4*np.pi*(y1 - yk))` where `dG` is the jump of the bound strength across that edge (the half-line
    formula of R12); `assert np.isclose(w_mine, LL.staircase_downwash(y1, ye, Gs))`; then N = 1, 3, 7, 21, 101 against
    `LL.downwash`. *expect:* mid-span downwash in units of Γ₁/2s: 0.637, 0.961, 0.991, 0.998, 0.9998 → 1.
11. `nb.plotly` — **F6** `slider_figure`: N ∈ {1, 3, 5, 7, 11, 21, 41} → the staircase Γ(y) with its smooth limit, and its
    downwash at 60 stations with the constant Γ₁/2s. *see / read* ("the downwash integral is the limit of this sum") */ change.*
12. `nb.note` — **N41 [B]** "Along the flight direction: $w(x\to-\infty)=0$, $w(0)=w$, $w(x\to+\infty)=2w$ — far behind, every filament looks
    infinite in both directions." `nb.figure` from `LL.downwash_along_x(x, np.array([1.0]), 1.0)` for x from −5 s to +5 s.
    *expect:* 0.000, 0.500, 1.000 at x = −50 s, 0, +50 s (Γ₁ = s = 1).
13. `nb.note` — **N65 [B]** `nb.figure`: sign map of `LL.horseshoe_downwash(X, Y, 10.0, 12.0)` on a grid in the wing's plane: down
    between the tips (blue), **up** outboard of them (orange). *expect:* −0.341 m/s (an upwash) at y = 8 m beside the bound
    vortex. "A bird flying just outboard and behind a neighbour's wing tip flies in rising air. > ⚠️ **slip #10 — the book
    prints** the 'Λ-pattern' of migrating birds**; the correct form is** a V formation — the letter there is a shape, not
    the aspect ratio."
14. `nb.note` — **N66 [B]** ground effect: 🔁 "the method of images (Ch. 5 §5.7, Ch. 6): a wall is replaced by a mirror vortex of
    opposite strength." "An image horseshoe below the ground reduces the downwash at the centre to
    $w=\dfrac{\Gamma}{\pi s}\Big[1-\dfrac1{1+(4h/s)^2}\Big]$ (ours)." `LL.ground_effect_factor(3.0, 12.0)` and two `LL.horseshoe_downwash` calls
    (the wing and its image at z = −2h). *expect:* 0.5 at h = s/4; 0.8 at h = s/2; 0.14 at h = s/10 — "a wing just above the
    runway floats".
15. `nb.note` — **N61 [C]** "A parabolic loading $\Gamma=\Gamma_0\big(1-(2y/s)^2\big)$: `LL.downwash(0.0, Gamma_fn, s)` gives 1.27 times the
    mid-span downwash of the elliptic loading with the same Γ₀; any Γ(y) is handled by the series of C11."
16. `nb.md` — **What would change if…** "…a wing section flies in air that is sinking at w? Its lift tilts. C10."

#### C10 — Induced drag: the lift vector tilted back by the downwash (14.14), (14.15)
1. `nb.core("C10", "Induced drag: $L=\\displaystyle\\int_{-s/2}^{s/2}\\rho U\\Gamma\\,dy$, $D_i=\\displaystyle\\int_{-s/2}^{s/2}\\rho w\\Gamma\\,dy$ (14.15)", question="How can
   a wing have drag in a fluid with no friction?")`
2. `nb.md` — **Plain words:** "Walk up a down escalator: you work just to stay level. A finite wing flies through air that
   its own wake is pushing downward. Relative to that air the wing is climbing a gentle slope, and the force on it, which
   is perpendicular to the air it actually meets, leans back by the slope angle. The backward part is a drag — with no
   friction anywhere. It is the price of lift on a wing with tips."
3. `nb.md` — **The idea** (ASCII triangle): `U ───────►` with `w ▼` at its tip, the slanted resultant `U_e`, the angle ε between;
   the force `dL_e` perpendicular to `U_e`, split into `dL ▲` and `dD_i ►`.
4. `nb.primer("relative velocity as a vector subtraction", "The wind a moving body feels is the air's velocity minus the
   body's own: $\\mathbf U_r=\\mathbf U-\\mathbf V$. Its direction, not the direction of the true wind, decides where 'lift' and 'drag'
   point.", code="U = np.array([50.0, 0.0]); w = np.array([0.0, -1.0])   # stream, and a 1 m/s downwash\nUe = U + w
   # the air the strip meets\nprint(np.hypot(*Ue), np.degrees(np.arctan2(1.0, 50.0)))   # 50.01 m/s, tilted 1.15 degrees")`
   (**P350**) + 🔁 "`np.arctan2(y, x)` gives the angle of a vector in the right quadrant (P70)."
5. `nb.derivation("D13", …)` — Part F D13 (8 steps), ref "14.15": includes **N42 [B]** (step 1; `> ⚠️ **slip #2 — the book prints**
   $\varepsilon=\tan\dfrac wU$**; the correct form is** $\varepsilon=\tan^{-1}\dfrac wU$`), **N43 [B]** $\alpha_e=\alpha-\varepsilon=\alpha-\dfrac wU$ (14.14) (step 2) and **N44 [B]**
   (steps 3–6) with `draw_force_triangle()` shown before it.
6. `nb.worked_example("one strip by hand", "U = 50 m/s, downwash w = 1 m/s, a strip of width 1 m carrying Γ = 20 m²/s, ρ = 1.2
   kg/m³, geometric angle α = 5°. 1. Tilt: ε = w/U = 0.02 rad = 1.15°. 2. Effective angle $\\alpha_e=\\alpha-\\dfrac wU$ (14.14) = 5° − 1.15°
   = 3.85°. 3. Lift of the strip: ρUΓ dy = 1.2 × 50 × 20 × 1 = 1200 N. 4. Induced drag: ρwΓ dy = 1.2 × 1 × 20 × 1 = 24 N. 5. Ratio
   24/1200 = 0.02 = ε: the drag is the lift times the tilt.")`
7. `nb.code` — `LL.strip_forces(10.0, 1.0, 30.0, 1.225, 0.1)`, `LL.strip_forces(…, exact=True)`, `LL.induced_angle(1.083, 30.0)`,
   `LL.effective_angle(np.deg2rad(5), 1.083, 30.0)`; then the whole wing with an elliptic loading that carries our aircraft at
   30 m/s: `G1 = 4*W/(np.pi*s*rho*U)`; `y = np.linspace(-6, 6, 4001)`; `LL.lift_and_induced_drag(y, LL.elliptic_circulation(y, s,
   G1), LL.elliptic_downwash(G1, s) + 0*y, rho, U)`. *expect:* dL = 36.75 N, dDi = 1.225 N, ε = 0.0333 rad; Γ₁ = 25.98 m²/s,
   w = 1.083 m/s, ε = 2.07°; L = 9000 N, $D_i$ = 324.8 N (3.6 % of the lift). *explain:* 4 points.
8. `nb.check_agree` — **from scratch:** `Ue = np.hypot(U, w)`; `eps = np.arctan2(w, U)`; `dLe = rho*Ue*G*dy`; `assert
   np.isclose(dLe*np.sin(eps), rho*w*G*dy)` and `np.isclose(dLe*np.cos(eps), rho*U*G*dy)` — "both are exact: $U_e\sin\varepsilon=w$
   and $U_e\cos\varepsilon=U$; only calling them 'lift' and 'drag' relative to the flight path is approximate."
9. `nb.figure` — the force triangle drawn to scale for the three stations y = 0, 3, 5.5 m of a rectangular-wing loading
   (from `LL.solve_lifting_line(1.25, 12.0, 30.0, 0.1, K=2*np.pi)`), ε exaggerated ×5 and labelled so. *see / read / change.*
10. `nb.note` — **N45 [B]** "Read $L=\int\rho U\Gamma\,dy$ and $D_i=\int\rho w\Gamma\,dy$ (14.15) as one rule used twice: U with Γ gives a force normal
    to U; w with Γ gives a force normal to w. With an infinite span there are no trailing vortices, w = 0 and the drag
    vanishes. (The book's remark that this is the largest part of a wing's drag is tested in C15, N75.)"
11. `nb.note` — **N38 [B]** "A second reading: the tip vortices get longer every second, so kinetic energy is fed into the wake
    at the rate $D_iU=\dfrac{d}{dt}(\text{wake kinetic energy})$ (stated, not derived)." `nb.code`: `LL.wake_kinetic_energy_per_length(Gn, s,
    rho)` against `LL.induced_drag_from_series(Gn, rho)` for the rectangular-wing coefficients. *expect:* equal to round-off
    (for one unit coefficient with ρ = 1: both π/8 = 0.3927). · **N67 [C]** "A tip vortex with a larger core holds less energy per
    metre, roughly $\propto\Gamma^2\big[\ln(s/\sigma)+\tfrac14\big]$ (ours): `ch14.vortex_pair_energy_per_length` with the Rankine vortex of Ch. 3."
12. `nb.explainer("downwash_induced_drag", heading="How can there be drag without friction?", why="Building the wing from 1, 3,
    7, 21 nested horseshoes and watching Γ(y) and w(y) settle onto smooth curves shows what the integral is; moving the
    station turns the local force triangle.", tries=["Start with one horseshoe and step N up.", "Slide the station to the tip:
    what happens to the bars of the nearest filaments?", "Move outboard of the tip: which way does the air move?", "Lower the
    wing toward the ground and watch the downwash."])`
13. `nb.md` — **What would change if…** "…Γ(y) is not given but has to be *found*? The downwash depends on Γ, and Γ depends on the
    downwash. We need a way to do the integral once and for all: C11."

#### C11 — The sine series and its downwash (14.18), (14.19)
1. `nb.core("C11", "The sine series and its downwash: $\\Gamma=\\sum\\Gamma_n\\sin(n\\gamma)$ (14.18) ⇒ $w(y_1)=\\dfrac{1}{2s}\\sum n\\Gamma_n\\dfrac{\\sin(n\\gamma_1)}{\\sin\\gamma_1}$
   (14.19, last member)", question="The downwash integral is singular and awkward. Can it be done once, for every loading?")`
2. `nb.md` — **Plain words:** "A violin string is fixed at both ends, and every shape it can take is a sum of sine waves that
   vanish at both ends. The circulation along a wing is 'fixed' the same way: zero at both tips. Write it as a sum of
   sines of an angle that runs from 0 at one tip to π at the other, and the awkward integral for the downwash can be done
   for each sine separately — the answer is another sine, divided by sin γ."
3. `nb.md` — **The idea**: `draw_span_semicircle()` — a semicircle over the span; the point at angle γ projects down to
   $y=-(s/2)\cos\gamma$; "equal steps in γ crowd toward the tips, where Γ changes fastest".
4. `nb.note` — **N49 [B]** "$y=-(s/2)\cos\gamma$, $y_1=-(s/2)\cos\gamma_1$: γ = 0 at the port tip, π/2 at mid-span, π at the starboard tip.
   > ⚠️ Common confusion (trap T3): this γ is an angle — not Chapter 5's sheet strength, not the ratio of specific heats."
   `LL.span_angle(np.array([-6.0, 0.0, 3.0, 6.0]), 12.0)` *expect:* 0, 1.571, 2.094, 3.142.
5. `nb.primer("half-range sine series", "On 0 ≤ γ ≤ π any function that vanishes at both ends can be written as $\\sum
   a_n\\sin n\\gamma$, and the sines are orthogonal there: $\\int_0^\\pi\\sin n\\gamma\\,\\sin m\\gamma\\,d\\gamma=\\tfrac\\pi2$ if n = m and 0 otherwise.",
   code="g = np.linspace(0, np.pi, 2001)\nprint(np.trapezoid(np.sin(2*g)*np.sin(3*g), g))   # 0: different sines\nprint(np.trapezoid(np.sin(3*g)**2,
   g), np.pi/2)   # 1.5708: the same sine")` (**P346**) + 🔁 "full-period integrals of sines (P151, Ch. 6); Galerkin truncation
   (P275)."
6. `nb.note` — **N50 [B]** $\Gamma=\displaystyle\sum_{n=1}^\infty\Gamma_n\sin(n\gamma)$ (14.18) + `nb.figure`: the first four shapes sin nγ against y —
   n = 1 a half-ellipse, n = 2 antisymmetric (a rolling wing), n = 3 a symmetric dip, n = 4 antisymmetric. *see / read* ("odd n
   are symmetric about mid-span, even n antisymmetric; every one is zero at both tips") */ change.*
7. `nb.derivation("D14", …)` — Part F D14 (7 steps), ref "14.19". After it: `> ⚠️ **slip #1 — the book prints** the first two
   members of (14.19) with the denominator $\cos\gamma_1-\cos\gamma$, $w(y_1)=\dfrac1{2\pi s}\displaystyle\int_0^\pi\sum n\Gamma_n\dfrac{\cos(n\gamma)\,d\gamma}{\cos\gamma_1-\cos\gamma}$**; the correct
   form is** $\cos\gamma-\cos\gamma_1$ (step 5). The last member, $w(y_1)=\dfrac1{2s}\sum n\Gamma_n\dfrac{\sin(n\gamma_1)}{\sin\gamma_1}$, and everything the book builds on it are
   right.` + `nb.code`: `LL.downwash_series(np.array([1.0]), np.pi/2, 1.0)`, the same with `printed_middle=True`, and
   `LL.downwash(0.0, lambda y: np.sqrt(1 - 4*y**2), 1.0)` (direct quadrature of $w(y_1)=\dfrac1{4\pi}\displaystyle\int\dfrac{d\Gamma}{dy}\dfrac{dy}{(y_1-y)}$ (14.13)).
   *expect:* +0.5, −0.5, +0.5 — "the printed form would make the air rise behind a lifting wing."
8. `nb.primer("product-to-sum identity and a three-term recursion", "$\\cos(n+1)\\gamma+\\cos(n-1)\\gamma=2\\cos\\gamma\\,\\cos n\\gamma$. A rule that
   gives each term of a sequence from the two before it, $I_{n+1}=2cI_n-I_{n-1}$, fixes the whole sequence once two starting
   values are known.", code="c, I = np.cos(1.0), [0.0, np.pi]                 # c = cos(gamma_1) with gamma_1 = 1 rad; I_0 and
   I_1\nfor n in range(1, 4): I.append(2*c*I[n] - I[n-1])   # the recursion\nprint(I[3], np.pi*np.sin(3*1.0)/np.sin(1.0))   # both
   0.5269")` (**P347**) + 🔁 "proof by induction (P306, Ch. 12); linear recurrences (P236)."
9. `nb.derivation("D15", …)` — Part F D15 (11 steps, ★★★, `check_src`), ref "14.19": **N51 [B]** Glauert's integral
   $\displaystyle\int_0^\pi\dfrac{\cos n\gamma\,d\gamma}{\cos\gamma-\cos\gamma_1}=\pi\,\dfrac{\sin n\gamma_1}{\sin\gamma_1}$ ("not in the book — ours, standard: the book says only that the last
   equality comes from evaluating the integral").
10. `nb.worked_example("two sines by hand", "Γ₁ = 12 and Γ₃ = 1 m²/s on a span s = 12 m. 1. Mid-span, γ₁ = 90°: sin γ₁ = 1,
    sin 3γ₁ = −1, so Γ = 12 − 1 = 11 m²/s. 2. Downwash there from $w(y_1)=\\dfrac1{2s}\\sum n\\Gamma_n\\dfrac{\\sin(n\\gamma_1)}{\\sin\\gamma_1}$ (14.19):
    (1/24)[1 × 12 × 1 + 3 × 1 × (−1)] = 9/24 = 0.375 m/s. 3. Near the port tip, γ₁ = 30° (y₁ = −5.2 m): sin γ₁ = 0.5, sin 3γ₁ = 1,
    Γ = 6 + 1 = 7 m²/s. 4. Downwash there: (1/24)[12 × 0.5/0.5 + 3 × 1 × 1/0.5] = (12 + 6)/24 = 0.75 m/s. 5. The n = 1 term alone
    gives 0.5 m/s everywhere; the small Γ₃ lowers the downwash in the middle and raises it by half toward the tips.")`
11. `nb.code` — `Gn = np.array([12.0, 0.0, 1.0])`; `LL.circulation_series(Gn, g1)`, `LL.downwash_series(Gn, g1, 12.0)` for g1 = π/2,
    π/6; `[LL.glauert_integral(n, np.pi/3, method="quad") for n in range(5)]` against `method="closed"`. *expect:* 11 and 7 m²/s;
    0.375 and 0.75 m/s; 0, π, π, 0, −π both ways. *explain:* 3 points.
12. `nb.check_agree` — **from scratch:** `Gfn = lambda y: 12*np.sin(g(y)) + np.sin(3*g(y))` with `g(y) = np.arccos(-2*y/12)`;
    `w_quad = LL.downwash(y1, Gfn, 12.0)` at five stations; `assert np.allclose(w_quad, LL.downwash_series(Gn, LL.span_angle(y1,
    12.0), 12.0), rtol=1e-6)`.
13. `nb.figure` — two panels: (left) w(y) of the single modes n = 1, 2, 3 (n = 1 is flat; the others swing and grow toward
    the tips); (right) the sign check of slip #1 as three bars (+0.5 last member, +0.5 quadrature, −0.5 printed middle
    member). *see / read / change.*
14. `nb.md` — **What would change if…** "…we do not know the $\Gamma_n$? Each strip of the wing tells us one relation between its
    circulation and its downwash. C12."

#### C12 — The lifting-line equation and its solver (14.17), (14.20)
1. `nb.core("C12", "The lifting-line equation: $\\dfrac K2Uc(\\alpha+\\beta)=\\sum\\Big(1+\\dfrac{nKc}{4s\\sin\\gamma}\\Big)\\Gamma_n\\sin(n\\gamma)$ (14.20)",
   question="Given a wing's shape and angle, how do you find the circulation along it?")`
2. `nb.md` — **Plain words:** "A designer draws a wing — this span, this taper, this twist — and wants its lift and its drag.
   Each strip of the wing behaves like the two-dimensional section of C06, but at the *reduced* angle the downwash leaves it.
   The downwash comes from the circulation of all the strips, and the circulation of each strip depends on its downwash:
   a loop. Written with the sine series the loop is a small set of linear equations, and a computer closes it in a
   millisecond."
3. `nb.md` — **The idea** (pipeline, ASCII): `chord c(y), angle α(y) ─► [each strip: Γ = (K/2)Uc(α − w/U + β)] ◄─► [w from all the
   Γ_n] ─► N unknowns Γ_n, N stations ─► np.linalg.solve ─► Γ(y), w(y), L, D_i`.
4. `nb.primer("integral equation", "An equation in which the unknown *function* appears under an integral sign. After
   discretising, it is a matrix equation with the unknown on both sides: x = b − M x, i.e. (I + M) x = b.", code="M =
   np.array([[0.2, 0.1], [0.1, 0.2]]); b = np.array([1.0, 1.0])   # x appears on both sides of x = b - M @ x\nx =
   np.linalg.solve(np.eye(2) + M, b)                         # move it to one side and solve\nprint(x, b - M @ x)
   # the same vector twice: [0.769 0.769]")` (**P345**) + 🔁 "`np.linalg.solve` (P57); collocation and the condition number
   (P162, Ch. 6)."
5. `nb.derivation("D16", …)` — Part F D16 (8 steps), ref "14.20": includes **N46 [B]** $C_L=K\Big[\alpha-\dfrac{w(y)}{U}+\beta(y)\Big]$ (14.16) (step 1;
   `> ⚠️` trap T12: "K is per radian; `K=` is a required keyword of every `LL` function, so that 2π is never assumed
   silently"; trap T3 for β), **N47 [B]** $\Gamma=\tfrac12UcC_L$ (step 2) and **N48 [B]** $\Gamma(y)=\dfrac K2Uc(y)\Big[\alpha-\dfrac{w(y)}{U}+\beta(y)\Big]$ (14.17) (step 3).
6. `nb.note` — **N52 [B]** the procedure, written out for N = 3: "keep three terms, demand the equation at the three stations
   γ = 45°, 90°, 135° ($\gamma_j=j\pi/(N+1)$ — never a tip, where sin γ = 0; trap T10): $\mathsf A_{jn}\Gamma_n=r_j$ with
   $\mathsf A_{jn}=\Big(1+\dfrac{nKc_j}{4s\sin\gamma_j}\Big)\sin n\gamma_j$ (ours)." `nb.code`: `A, g = LL.lifting_line_matrix(1.0, 6.0, 2*np.pi, 3)`;
   `r = LL.lifting_line_rhs(1.0, 6.0, 1.0, wd["alpha"], 2*np.pi, 3)`; `np.linalg.solve(A, r)`. *expect:* rows [0.969, 1.740, 1.493],
   [1.262, 0, −1.785], [0.969, −1.740, 1.493]; r = 0.2742 three times; Γ = [0.2487, 0, 0.0222] m²/s.
7. `nb.worked_example("one term, one station, by hand", "A rectangular wing with s = 6 m, c = 1 m (Λ = 6), K = 2π, U = 1 m/s,
   α = 5° = 0.0873 rad, β = 0. Keep only Γ₁ and use the mid-span station γ = 90° (sin γ = 1). 1. The bracket of
   $\\dfrac K2Uc(\\alpha+\\beta)=\\sum\\Big(1+\\dfrac{nKc}{4s\\sin\\gamma}\\Big)\\Gamma_n\\sin(n\\gamma)$ (14.20): 1 + 2π × 1/(4 × 6) = 1.262. 2. Left side: π × 1 × 1 × 0.0873 =
   0.2742. 3. Γ₁ = 0.2742/1.262 = 0.2173 m²/s. 4. Lift coefficient: $L=\\dfrac{\\pi s}4\\rho U\\Gamma_1$ (14.21) divided by ½ρU²sc gives
   πΓ₁/(2Uc) = 0.341. 5. The two-dimensional section alone would give 2π × 0.0873 = 0.548: the finite wing delivers 62 % of
   it in this crude estimate (72 % when enough terms are kept — next cell).")`
8. `nb.code` — `sol = LL.solve_lifting_line(LL.planform("rectangular", 6.0, 1.0), 6.0, 1.0, wd["alpha"], K=2*np.pi, N=32)`; print
   `sol["CL"], sol["slope"], sol["CDi"], sol["e_span"], sol["cond"]`, `sol["Gn"][:5]/sol["Gn"][0]`. *expect:* $C_L$ = 0.3954, slope
   4.530 per radian, $C_{D_i}$ = 0.00869, $e_{\text{span}}$ = 0.954, condition number ≈ 40; ratios 1, 0, 0.122, 0, 0.026. *explain:* 5
   points ("…5. the even coefficients vanish: a symmetric wing has a symmetric loading").
9. `nb.check_agree` — **from scratch (N = 4):** two nested loops fill `A[j, n-1] = (1 + n*K*c/(4*s*np.sin(g[j])))*np.sin(n*g[j])`
   with `g[j] = (j + 1)*np.pi/5`; `r[j] = K/2*U*c*alpha`; `Gn = np.linalg.solve(A, r)`; `assert np.allclose(A,
   LL.lifting_line_matrix(c, s, K, 4)[0])` and `np.allclose(Gn, LL.solve_lifting_line(c, s, U, alpha, K=K, N=4)["Gn"])`.
10. `nb.code` — **convergence table** N = 1, 2, 4, 8, 16, 32, 64 → $C_L$, $C_{D_i}$, $e_{\text{span}}$, `cond`, and the largest residual of
    $\Gamma(y)=\dfrac K2Uc\Big[\alpha-\dfrac wU+\beta\Big]$ (14.17) between stations from `LL.lifting_line_residual`. *expect:* $C_L$ = 0.3413, 0.3818, 0.3935,
    0.3952, 0.39534, 0.39535, 0.39535; $e_{\text{span}}$ = 1, 1, 0.9648, 0.9550, 0.9540, 0.9539, 0.9539; cond 1 → 154. "A few
    terms usually suffice."
11. `nb.plotly` — **F7** `slider_figure`: taper ratio 0.1…1.0 (19 steps) at Λ = 6 → Γ(y) (amber) with the elliptic loading of the
    same lift (muted), w(y)/U (blue) with its uniform ghost, $e_{\text{span}}$ in the title. *expect:* $e_{\text{span}}$ = 0.954, 0.970,
    0.983, 0.991, 0.982, 0.958 at taper 1.0, 0.8, 0.6, 0.4, 0.2, 0.1 — best near 0.4. *see / read / change.*
12. `nb.live` — taper, twist at the tip, Λ, N sliders → loading and coefficient bars; paired with F7.
13. `nb.note` — **N63 [B]** twist: "A *rectangular* wing can carry an elliptic loading if its angle varies along the span like an
    ellipse: tip angle $\alpha_t=\dfrac{Kc\,\alpha_r}{4s}$ and $\Gamma_1=\tfrac12KUc\,\alpha_r$ (ours; stated, then confirmed by the solver).
    > ⚠️ **slip #9 — the book prints** the span as b in the two exercises on this**; the correct form is** s (b is the map
    constant)." `tw = LL.elliptic_twist_rectangular(wd["alpha"], 1.0, 6.0, K=2*np.pi)`; `LL.solve_lifting_line(1.0, 6.0, 1.0, 0.0,
    K=2*np.pi, twist=lambda y: tw["alpha_t"] + wd["alpha"]*np.sqrt(1 - (2*y/6.0)**2), N=8)["Gn"]`. *expect:* α_t = 1.31°; Γ =
    [0.2742, 0, 0, …]: one coefficient only.
14. `nb.note` — **N64 [C]** "A small aircraft flying along a large one's tip vortex feels a rolling moment
    $M=\displaystyle\int_{-s/2}^{s/2}\rho U\,y\,\Gamma(y)\,dy$: `ch14.wake_vortex_rolling_moment` (one number, our inputs) — the reason for
    wake-turbulence separation at airports. (The exercise that poses it takes w positive upward — trap T5.)"
15. `nb.note` — **N68 [B]** when it fails: "Lifting-line results are approximate. They need a large aspect ratio and an unswept
    wing (assumptions (i), (ii) of N39), a wake that stays flat (iii), and sections below stall — the section law
    $C_L=K\big[\alpha-\tfrac wU+\beta\big]$ (14.16) is a straight line. Extensions cover several wings, propellers, ground effect and
    bird formations."
16. `nb.md` — **What would change if…** "…we only want the wing's total lift and drag, not Γ(y)? They follow from the
    coefficients directly: C13."

#### C13 — Lift and induced drag from the coefficients (14.21), (14.22)
1. `nb.core("C13", "Lift and induced drag from the coefficients: $L=\\dfrac{\\pi s}{4}\\rho U\\Gamma_1$ (14.21), $D_i=\\dfrac\\pi8\\rho\\sum n\\Gamma_n^2$ (14.22)",
   question="Which parts of the loading lift the aircraft, and which only cost drag?")`
2. `nb.md` — **Plain words:** "Split a wing's loading into its sine shapes. The first is a single hump. The others wiggle:
   they add lift in one place and remove it in another, and over the whole span their lift adds up to nothing. But every
   wiggle sheds extra trailing vorticity, and trailing vorticity costs drag. So the lift is in the first shape alone —
   and every other shape is a pure cost."
3. `nb.md` — **The idea** (table): | term | net lift | induced drag | — n = 1 | all of it | Γ₁² · n = 2 | 0 | 2Γ₂² · n = 3 | 0 | 3Γ₃² ·
   n | 0 | nΓ_n².
4. `nb.primer("a double sum collapsed by orthogonality", "Multiplying two sums gives every pairing: $(\\sum_na_n)(\\sum_mb_m)=
   \\sum_n\\sum_ma_nb_m$ — use two different index letters. If the pairings are then weighted by a factor that is 1 when n = m and 0
   otherwise ($\\delta_{nm}$), only the diagonal survives: $\\sum_n\\sum_ma_nb_m\\delta_{nm}=\\sum_na_nb_n$.", code="a = np.array([1., 2., 3.]); b
   = np.array([4., 5., 6.])\nprint(np.sum(np.outer(a, b)*np.eye(3)), np.sum(a*b))   # 32.0 both: the cross terms are gone")`
   (**P348**)
5. `nb.derivation("D17", …)` — Part F D17 (6 steps), ref "14.21": **N53 [B]**.
6. `nb.derivation("D18", …)` — Part F D18 (9 steps, ★★★, `check_src`), ref "14.22" (slip **#7** at step 9).
7. `nb.note` — **N54 [B]** "For a given lift (a given Γ₁) the induced drag is least when every higher coefficient is zero. **Not in
   the book — ours:** the span efficiency $e_{\text{span}}\equiv\dfrac{\Gamma_1^2}{\sum n\Gamma_n^2}\le1$, with $C_{D_i}=\dfrac{C_L^2}{\pi\Lambda e_{\text{span}}}$. > ⚠️ Common
   confusion (trap T3): aeronautics writes this as *e*; this chapter already uses e for the thickness parameter, so we
   write $e_{\text{span}}$ (`e_span` in code)." + `nb.figure`: bars of $n\Gamma_n^2$ for the rectangular demo wing (N = 12), first bar teal
   ("lift"), the rest rose ("drag only"). *expect:* shares 95.4 %, 4.2 %, 0.3 % for n = 1, 3, 5. *see / read / change.*
8. `nb.worked_example("two sines again", "Γ₁ = 12, Γ₃ = 1 m²/s, s = 12 m, ρ = 1.2 kg/m³, U = 50 m/s. 1. $L=\\dfrac{\\pi s}4\\rho U\\Gamma_1$ (14.21) =
   (π × 12/4) × 1.2 × 50 × 12 = 6786 N — Γ₃ does not appear. 2. $D_i=\\dfrac\\pi8\\rho\\sum n\\Gamma_n^2$ (14.22) = (π/8) × 1.2 × (1 × 144 + 3 × 1) =
   0.4712 × 147 = 69.3 N. 3. Without Γ₃ it would be 0.4712 × 144 = 67.9 N. 4. Span efficiency: 144/147 = 0.980 — the small third
   harmonic costs 2 % more drag and buys no lift.")`
9. `nb.code` — `LL.lift_from_series(sol["Gn"], 6.0, 1.0, 1.0)`, `LL.induced_drag_from_series(sol["Gn"], 1.0)`,
   `LL.span_efficiency(sol["Gn"])`, `LL.induced_drag_coefficient(sol["CL"], 6.0, sol["e_span"])`. *expect:* L = 1.186 N, $D_i$ = 0.0261
   N (ρ = 1, U = 1), $e_{\text{span}}$ = 0.954, $C_{D_i}$ = 0.00869 — the same as the solver's own. *explain:* 3 points.
10. `nb.check_agree` — **from scratch:** `g = np.linspace(0, np.pi, 4001)`; `y = -3*np.cos(g)`; `Gam = LL.circulation_series(Gn, g)`;
    `w = LL.downwash_series(Gn, g, 6.0)`; `L_mine = np.trapezoid(rho*U*Gam, y)`; `Di_mine = np.trapezoid(rho*w*Gam, y)`; `assert
    np.allclose([L_mine, Di_mine], [LL.lift_from_series(Gn, 6.0, rho, U), LL.induced_drag_from_series(Gn, rho)], rtol=1e-5)` —
    the integrals of $L=\int\rho U\Gamma\,dy$, $D_i=\int\rho w\Gamma\,dy$ (14.15) against the two closed forms.
11. `nb.explainer("lifting_line_solver", heading="What makes one wing better than another?", why="Changing taper, twist and
    aspect ratio and watching the coefficient bars, the loading against its elliptic ghost, the downwash and the efficiency
    respond together shows the method's inputs and outputs at once.", tries=["Go from the rectangle to taper 0.4: what happens
    to the rose bars?", "Press 'rectangle with elliptic twist': how many bars are left?", "Raise N from 1 to 12 and watch the
    numbers settle.", "Click the n = 3 bar to see its shape on the span."])`
12. `nb.md` — **What would change if…** "…we switch every rose bar off? C14."

#### C14 — The elliptic loading: uniform downwash and the least induced drag (14.23)–(14.26)
1. `nb.core("C14", "The elliptic loading: $C_{D_i}=\\dfrac{C_L^2}{\\pi\\Lambda}$ (14.26)", question="What is the cheapest way to carry a given
   weight on a given span — and why are glider wings so long?")`
2. `nb.md` — **Plain words:** "An albatross glides for hours on wings three metres across and a hand wide. A sparrow could not.
   The reason is in one formula: the drag a wing pays for its lift is inversely proportional to the *square of its span*.
   Spread the same lift over twice the span and the bill drops to a quarter. And among all ways of spreading the lift
   along a given span, one is cheapest: the loading shaped like half an ellipse, which makes the air sink at the same speed
   everywhere along the wing."
3. `nb.md` — **The idea** (ASCII): a half-ellipse Γ(y) over the span with a flat line w = Γ₁/2s under it; "one sine ⇒ uniform
   downwash ⇒ no part of the wing works in worse air than another".
4. `nb.derivation("D19", …)` — Part F D19 (9 steps), ref "14.26": includes **N55 [B]** $\Gamma=\Gamma_1\sin(\gamma)=\Gamma_1\sqrt{1-(2y/s)^2}$ (14.23) (steps
   3–4), **N56 [B]** $w(y)=\Gamma_1/2s$ (14.24) (step 5) and **N57 [B]** $D_i=\dfrac\pi8\rho\Gamma_1^2=\dfrac{2L^2}{\pi\rho U^2s^2}$ (14.25) (steps 6–7).
5. `nb.note` — **N58 [B]** $dL=\rho U\Gamma\,dy=C_L\,\tfrac12\rho U^2c\,dy$ (14.27): "with a uniform downwash, one section and no twist, every
   strip has the same $C_L$, so the chord must be proportional to the circulation: an **elliptic planform** gives the elliptic
   loading." `> ⚠️ Common confusion (trap T11):` "the optimum is a statement about the *loading* Γ(y), not about the outline.
   An elliptic outline produces it only if untwisted with one section; a rectangle can produce it with twist (N63, C12); and
   'optimum' means least induced drag at a fixed lift **and span**, with a flat wake."
6. `nb.derivation("D20", …)` — Part F D20 (8 steps), ref "14.27": **N59 [B]** $C_{L,3D}=\dfrac{C_{L,2D}}{1+2/\Lambda}$ (K = 2π) and
   $\dfrac{dC_L}{d\alpha}=\dfrac{K}{1+K/(\pi\Lambda)}$ (ours).
7. `nb.worked_example("the price of lift by hand", "**N62 [C]** (the pattern of the book's level-flight exercise, with our own
   aircraft of C01: W = 10 000 N, s = 10 m, A = 20 m², Λ = 5, ρ = 1.25 kg/m³ from the sea-level standard atmosphere of Ch. 1,
   U = 40 m/s, $C_L$ = 0.5.) 1. $D_i=\\dfrac{2L^2}{\\pi\\rho U^2s^2}$ (14.25) = 2 × 10⁸/(π × 1.25 × 1600 × 100) = 318 N. 2. As a coefficient,
   $C_{D_i}=\\dfrac{C_L^2}{\\pi\\Lambda}$ (14.26) = 0.25/15.7 = 0.0159; check: 0.0159 × 1000 Pa × 20 m² = 318 N ✓. 3. First coefficient from
   $L=\\dfrac{\\pi s}4\\rho U\\Gamma_1$ (14.21): Γ₁ = 4 × 10 000/(π × 10 × 1.25 × 40) = 25.5 m²/s. 4. Downwash $w=\\Gamma_1/2s$ (14.24) = 1.27 m/s,
   a tilt of 1.27/40 = 0.032 rad = 1.8°; lift × tilt = 10 000 × 0.0318 = 318 N ✓. 5. Circulation half-way out along each wing
   (y = ±2.5 m): 25.5 × √(1 − 0.25) = 22.1 m²/s. 6. Double the span at the same area (Λ = 20): 318/4 = 80 N.")`
8. `nb.code` — `LL.elliptic_induced_drag(ac["W"], 1.225, 30.0, ac["s"])`, at 55 m/s; `LL.induced_drag_coefficient(1.088, 9.6)`;
   `LL.finite_wing_lift_slope(2*np.pi, 9.6)`, `LL.finite_wing_lift_slope(2*np.pi, 6.0)`; `LL.elliptic_wing(wd["alpha"], 6.0,
   K=2*np.pi)`. *expect:* 324.8 N and 96.6 N; 0.0393; 5.200 and 4.712 per radian (0.75 × 2π at Λ = 6); $C_L$ = 0.4112, $C_{D_i}$ =
   0.00897, w/U = 0.0218. *explain:* 4 points ("…2. flying 1.8 times faster cuts the induced drag by 3.4: it falls as 1/U²").
9. `nb.check_agree` — **two checks:** (i) `[LL.downwash(y1, lambda y: np.sqrt(1 - (2*y/12)**2), 12.0) for y1 in (-5, -2, 0, 3,
   5.5)]` all equal `LL.elliptic_downwash(1.0, 12.0)` (direct quadrature of $w(y_1)=\dfrac1{4\pi}\displaystyle\int\dfrac{d\Gamma}{dy}\dfrac{dy}{(y_1-y)}$ (14.13): uniform);
   (ii) `LL.solve_lifting_line(LL.planform("elliptic", 6.0, 4/np.pi), 6.0, 1.0, wd["alpha"], K=2*np.pi, N=8)`: `assert
   np.allclose(sol["Gn"][1:], 0, atol=1e-12)` and `np.isclose(sol["CL"], LL.elliptic_wing(wd["alpha"], 6.0, K=2*np.pi)["CL"])`.
   *expect:* 0.04167 five times; one coefficient 0.2618; $C_L$ = 0.4112.
10. `nb.figure` — `fig_planforms_loadings()`: three planforms of equal span and area (rectangle, taper 0.4, ellipse) with
    their Γ(y), w(y) and $e_{\text{span}}$ from the solver. *expect:* 0.954, 0.991, 1.000. *see / read* ("only the ellipse has a flat
    downwash") */ change.*
11. `nb.plotly` — **F8a** `slider_figure`: Λ from 3 to 30 (28 steps) → $C_{D_i}$ against $C_L$ (0…1.5) from
    $C_{D_i}=\dfrac{C_L^2}{\pi\Lambda}$ (14.26), with the lift slope $K/(1+K/\pi\Lambda)$ in the title. *see / read / change.*
12. `nb.md` — **What would change if…** "…we add the drag that has nothing to do with lift — friction and form drag? Then there
    is a best speed to fly: C15."

### A.7 §14.7 Lift and Drag Characteristics of Airfoils — C15
#### C15 — The drag breakdown and the polar (ours, built on (14.26))
1. `nb.section("14.7", "Lift and Drag Characteristics of Airfoils", intro="**What is this section about?** What measured lift
   and drag curves look like, where they leave the theory (stall), and how the drag of a wing splits into a part that comes
   from friction and shape and a part that is the price of lift.")`
2. `nb.core("C15", "The drag breakdown and the polar: $C_D=C_{D0}+\\dfrac{C_L^2}{\\pi\\Lambda e_{\\text{span}}}$ (ours, from (14.26))", question="Fly
   slowly and one kind of drag is large; fly fast and another is. Where is the best speed?")`
3. `nb.md` — **Plain words:** "A glider pilot knows a number called the best-glide speed. Slower than that the glider sinks
   faster; faster than that it also sinks faster. Two drags are fighting. Friction and the wake of the body grow with the
   square of the speed. The price of lift does the opposite: the slower you fly, the harder the wing must work to hold the
   same weight, and the induced drag grows as one over the square of the speed. The total has a minimum where the two are
   equal."
4. `nb.md` — **The idea** (tree, ASCII): `total drag = friction drag + pressure drag`; `pressure drag = form drag + induced drag`;
   `profile drag = friction + form` (what the section alone would have, $C_{D0}$) `+ induced drag` (the wing's tips,
   $C_L^2/(\pi\Lambda e_{\text{span}})$). 🔁 "friction and form drag: Ch. 9 §9.8 and R02 in C01."
5. `nb.note` — **N70 [B]** measured curves: `nb.figure` (title "**Qualitative model** of measured section curves — every constant
   is ours"): `ch14.generic_lift_curve(alpha, K=5.7, beta=sec["beta"], alpha_stall=np.deg2rad(14), CL_post=0.9,
   width=np.deg2rad(2))` and `ch14.generic_drag_curve(alpha, CD0=0.008, AR=np.inf, e_span=1.0, K=5.7, beta=sec["beta"], alpha_stall=np.deg2rad(14), CL_post=0.9, width=np.deg2rad(2), CD_sep=0.15)` against α. "A straight range $C_L\simeq K(\alpha+\beta)$, lift at α = 0 from the camber,
   a maximum, a drop at stall; the drag rises roughly as α²." *see / read / change.* · **N71 [C]** `draw_stalled_section()`: "past
   stall the flow leaves the upper surface near the nose and the wake is as wide as the section is tall; the mechanism is
   N13 (C07)."
6. `nb.note` — **N72 [B]** "How much lift a section reaches before it stalls depends on the chord Reynolds number
   $\mathrm{Re}_c=Uc/\nu$: at moderate values the laminar layer separates before it turns turbulent (early stall); at high values
   it turns turbulent first, and a turbulent layer survives a steeper pressure rise — the same physics as the drag crisis of
   a sphere (Ch. 9 §9.9)." `ch14.chord_reynolds(30.0, 1.25, nu_air)`, at 55 m/s, `ch14.stall_regime(2.57e6,
   Re_laminar=5e5, Re_turbulent=3e6)`. *expect:* 2.57 × 10⁶, 4.71 × 10⁶; "mixed" (thresholds ours, **qualitative**).
7. `nb.note` — **N73 [B]** "Camber lifts the whole curve; a flap is camber you can switch on; so the slowest speed at which the
   wing still carries the weight, $U_{\min}=\sqrt{\dfrac{2W}{\rho A\,C_{L,\max}}}$ (ours; D01 step 5), drops when the flap goes down."
   `ch14.stall_speed(ac["W"], 1.225, ac["A"], 1.4)`, with 2.0. *expect:* 26.5 and 22.1 m/s.
8. `nb.primer("reading a drag polar; the minimum of a/x + bx", "A polar plots $C_L$ (up) against $C_D$ (across). A straight line
   from the origin to a point of the curve has slope $C_L/C_D$ = L/D; the steepest such line — the tangent — marks the best
   glide. Any function a/x + bx has its minimum where the two terms are equal, at x = √(a/b), value 2√(ab).", code="a, b_ =
   0.02, 0.0637                       # C_D0 and 1/(pi*AR) for AR = 5\nx = np.linspace(0.1, 1.5, 1401)             # candidate lift
   coefficients\nprint(x[np.argmin(a/x + b_*x)], np.sqrt(a/b_))  # 0.56 both: the minimum of C_D/C_L")` (**P349**) + 🔁 "quotient
   rule (P265); `minimize_scalar` (P170)."
9. `nb.derivation("D21", …)` — Part F D21 (8 steps), ref "14.26": **N74 [B]** $\Big(\dfrac LD\Big)_{\max}=\tfrac12\sqrt{\dfrac{\pi\Lambda e_{\text{span}}}{C_{D0}}}$ at
   $C_L=\sqrt{\pi\Lambda e_{\text{span}}\,C_{D0}}$ (ours; the book says only that $C_{D_i}\propto\alpha^2$).
10. `nb.worked_example("a polar by hand", "Λ = 5, $e_{\\text{span}}$ = 1, so πΛ = 15.7; $C_{D0}$ = 0.02. 1. Best glide where induced =
    profile: $C_L^2/15.7$ = 0.02 → $C_L$ = √0.314 = 0.56. 2. There $C_D$ = 0.02 + 0.02 = 0.04 and L/D = 0.56/0.04 = 14. 3. Formula:
    ½√(15.7/0.02) = ½ × 28 = 14 ✓. 4. Fly faster, at half that lift coefficient (0.28): induced 0.005, total 0.025, L/D = 11.2,
    induced share 20 %. 5. Fly slower, at twice it (1.12): induced 0.08, total 0.10, L/D = 11.2, induced share 80 %. Same glide
    ratio on both sides — for opposite reasons.")`
11. `nb.code` — `ch14.drag_polar(1.088, CD0=ac["CD0"], AR=9.6, e_span=ac["e_span"])`; `ch14.drag_breakdown(1.088, CD0=ac["CD0"], AR=9.6, e_span=ac["e_span"])` (coefficients); `ch14.drag_vs_speed(ac["W"], 1.225, ac["A"], 9.6, ac["CD0"], 30.0, ac["e_span"])` (newtons), at 55 m/s; `ch14.best_glide(ac["CD0"], 9.6, ac["e_span"])`;
    `ch14.speed_min_drag(ac["W"], 1.225, ac["A"], 9.6, ac["CD0"], ac["e_span"])`. *expect:* $C_D$ = 0.0687; at 30 m/s: profile 207 N
    + induced 361 N = 568 N, induced share 63.6 %, L/D = 15.9; at 55 m/s: 695 N + 107 N = 802 N, share 13.4 %, L/D = 11.2; best
    glide L/D = 16.5 at $C_L$ = 0.824, speed 34.5 m/s, drag 546 N. *explain:* 5 points.
12. `nb.check_agree` — **from scratch:** `CL = np.linspace(0.05, 1.6, 15501)`; `LD = CL/(CD0 + CL**2/(np.pi*AR*e_span))`; `i =
    np.argmax(LD)`; `bg = ch14.best_glide(CD0, AR, e_span)`; `assert np.isclose(LD[i], bg["LD_max"], rtol=1e-6)` and
    `np.isclose(CL[i], bg["CL"], atol=1e-4)`.
13. `nb.note` — **N75 [B]** `> ⚠️ Common confusion (trap T9):` "The book says that induced drag is normally the largest part of
    a wing's drag. It is a statement about a flight condition, not about wings." `nb.figure`: induced share
    `ch14.drag_share(CL(U), CD0=…, AR=…, e_span=…)` against speed from the stall speed to 70 m/s, with the 50 % line and the
    best-glide speed marked. *expect:* 50 % exactly at 34.5 m/s; 64 % at 30 m/s; 13 % at 55 m/s. "Exactly half at best glide,
    more in a slow climb, less at cruise — a published figure for an airliner at cruise is about 37 %." *see / read / change.*
14. `nb.plotly` — **F8b** two `slider_figure`s over Λ = 4…24 (21 steps) at fixed area and weight: (a) the polar $C_L$ against
    $C_D$ with the tangent from the origin and L/D max in the trace name; (b) drag in newtons against speed — profile (blue),
    induced (rose), total — with the minimum marked. *see / read* ("a slender wing flattens the polar and moves the best speed
    down") */ change.*
15. `nb.note` — **N76 [C]** "Near the speed of sound local shock waves appear on the wing; sweeping it back reduces the Mach
    number of the flow component normal to the leading edge, $\mathrm{Ma}_n=\mathrm{Ma}\cos\Lambda_{\text{sweep}}$ (ours). Shock waves are Ch. 15
    §15.6 and thin sections in supersonic flow §15.10. (The sweep angle shares the letter Λ with the aspect ratio.)"
16. `nb.explainer("elliptic_loading_aspect_ratio", heading="Why are glider wings so long and thin?", why="Stretching the span
    at fixed area reshapes the planform, flattens the polar and slides the minimum of the U-shaped drag curve at once; the
    speed cursor swaps which bar is larger.", tries=["Drag the speed from stall to cruise and watch the two bars cross.",
    "Double Λ at fixed area: what happens to the best-glide speed and ratio?", "Switch the mode from glider to airliner.",
    "Find the speed at which the induced share is exactly 50 %."])`
17. `nb.md` — **What would change if…** "…the wing itself moved across the stream? Then its lift would lean *forward*. §14.8 and §14.9."

### A.8 §14.8 Propulsive Mechanisms of Fish and Birds — N77, N79 (and N78, N80) open the C16 story
1. `nb.section("14.8", "Propulsive Mechanisms of Fish and Birds", intro="**What is this section about?** A wing gives a force
   across the stream. Move the wing sideways through the water and 'across the stream it feels' has a forward part. A tail
   fin and a flapping wing are lifting surfaces used as engines.")`
2. `nb.note` — **N77 [B]** the fish tail, with `draw_velocity_triangle("fish")`: "The fish swims forward at U while its tail, a
   symmetric section, sweeps sideways at V. The water meets the tail along the diagonal of that triangle, at the angle θ with
   $V/U=\tan\theta$. The lift L is perpendicular to the diagonal: its forward part is the thrust $L\sin\theta$, delivered at the rate
   $UL\sin\theta$; its sideways part $L\cos\theta$ is what the muscles work against, at the rate $VL\cos\theta$." + `nb.worked_example("a tail
   stroke by hand", "U = 1.0 m/s, V = 0.75 m/s, L = 10 N. 1. tan θ = 0.75, so sin θ = 0.6 and cos θ = 0.8 (a 3–4–5 triangle; θ =
   36.9°). 2. Thrust = 10 × 0.6 = 6 N; side force = 10 × 0.8 = 8 N. 3. Useful power = U × thrust = 1.0 × 6 = 6 W. 4. Power put in =
   V × side force = 0.75 × 8 = 6 W. 5. They are equal — always, because $V=U\\tan\\theta$ makes $VL\\cos\\theta=UL\\sin\\theta$. 6. With a drag
   of 1 N along the diagonal: thrust 6 − 0.8 = 5.2 N, side force 8 + 0.6 = 8.6 N, efficiency (5.2 × 1.0)/(8.6 × 0.75) = 0.81.")` +
   `nb.code`: `ch14.oscillating_foil(1.0, 0.75, 10.0)`, `ch14.oscillating_foil(1.0, 0.75, 10.0, 1.0)`. *expect:* thrust 6 N, side 8 N,
   power in = power out = 6 W, efficiency 1; with drag: 5.2 N, 8.6 N, efficiency 0.806. `> ⚠️ Common confusion (trap T18):` "the
   perfect score is an identity of the triangle, not a prediction: the book's picture has no drag and no wake, and treats
   each instant as steady."
3. `nb.note` — **N79 [B]** bird flight, with `draw_velocity_triangle("bird_down")` and `("bird_up")`: "On the downstroke the wing
   moves down and forward through the air; the lift, normal to the relative wind, points up *and* forward: weight support
   and thrust at once. On the upstroke the wing is feathered — nearly edge-on to its relative wind — and does little."
   `nb.code`: `ch14.oscillating_foil(10.0, 4.0, 2.0)` (downstroke) and `ch14.oscillating_foil(10.0, 4.0, 0.2)` (upstroke, a tenth of
   the lift); mean thrust over the cycle. *expect:* thrust 0.743 N and 0.074 N; mean 0.41 N.
4. `nb.note` — **N78 [C]** "Whales and dolphins move a horizontal tail up and down — the same triangle turned on its side;
   bumps on a humpback's flipper delay separation (N13, N72); flexible fins change camber within a stroke. Biological flows
   are Ch. 16." · **N80 [C]** "Measured bird-wing shapes and a computed flight of a fruit fly are quoted in the book; the
   numbers stay there. An insect flies at so low a Reynolds number that 'thin boundary layers' (N01) no longer holds."

### A.9 §14.9 Sailing against the Wind — C16; exercises, literature, summary
#### C16 — Thrust from lift: the apparent wind and the sail (ours: the book resolves the force in words)
1. `nb.section("14.9", "Sailing against the Wind", intro="**What is this section about?** The same velocity triangle a third
   time: a sail is a wing standing on end in the wind the boat *feels*. As long as that wind comes from far enough off the
   bow, the sail's lift has a forward part — even when the boat is heading into the wind.")`
2. `nb.core("C16", "Thrust from lift: $\\mathbf U_r=\\mathbf U-\\mathbf V$, drive $=L\\sin\\beta_a-D\\cos\\beta_a$ (ours)", question="How can a boat
   sail against the wind that drives it?")`
3. `nb.md` — **Plain words:** "'The wind pushes the sail' explains sailing downwind and nothing else. A sail trimmed close is a
   wing: air flows *along* both sides and the force is mostly sideways to the wind the boat feels. That felt wind — the
   apparent wind — is the true wind minus the boat's own velocity, so it always comes from further ahead than the true wind.
   Lift is perpendicular to it; and a force perpendicular to a wind from ahead-and-to-the-side points partly forward. The
   sideways rest of the force would shove the boat downwind, were it not for a second wing under the water: the keel."
4. `nb.md` — **The idea** (ASCII plan view): the boat heading up the page; true wind `U` from upper right; boat velocity `V` up;
   apparent wind `U_r = U − V` from closer to the bow, at the angle β_a; lift `L ⟂ U_r` pointing forward-left, drag `D ∥ U_r`;
   "drive = forward part of L − backward part of D". 🔁 "relative velocity as a vector subtraction (P350, primed in C10); the
   fish tail of §14.8 is this triangle with the stream and the motion exchanged."
5. `nb.derivation("D22", …)` — Part F D22 (8 steps), ref "14.15": includes **N81 [B]** (steps 7–8: the keel and the closest course
   $\beta_a=\varepsilon_A+\varepsilon_H$, $\varepsilon=\tan^{-1}(D/L)$ of sail and of hull — "not in the book — ours").
6. `nb.worked_example("a beam reach by hand", "True wind U = 4 m/s from exactly abeam (course angle 90°), boat speed V = 3 m/s.
   1. The apparent wind has the part U = 4 m/s from the side and the part V = 3 m/s from ahead: speed 5 m/s (3–4–5). 2. Its
   angle off the bow: tan $\\beta_a$ = 4/3, $\\beta_a$ = 53.1°, sin 0.8, cos 0.6. 3. The sail gives L = 100 N and D = 25 N. 4. Drive =
   $L\\sin\\beta_a-D\\cos\\beta_a$ = 80 − 15 = 65 N. 5. Side force = $L\\cos\\beta_a+D\\sin\\beta_a$ = 60 + 20 = 80 N, carried by the keel. 6. The
   drive would vanish at tan $\\beta_a$ = D/L = 0.25, i.e. 14° off the bow. 7. Head up to 45° from the true wind: tan $\\beta_a$ =
   2.83/(2.83 + 3) = 0.485, $\\beta_a$ = 25.9° — still more than 14°, so the boat still drives forward, *toward* the wind.")`
7. `nb.code` — `ch14.apparent_wind(boat["U_wind"], boat["V_boat"], boat["theta"])`; `ch14.sail_forces(boat["U_wind"], boat["V_boat"],
   boat["theta"], CL=boat["CL"], CD=boat["CD"], rho=1.225, A=boat["A_sail"])`; `ch14.sailboat_course_limit(boat["CL"]/boat["CD"],
   boat["LD_hull"])`; `ch14.sailboat_equilibrium(boat["U_wind"], boat["theta"], CL=boat["CL"], CD=boat["CD"], rho_air=1.225, A_sail=boat["A_sail"], rho_water=boat["rho_water"], CDA_hull=boat["CDA_hull"], LD_hull=boat["LD_hull"])`; `ch14.velocity_made_good(3.0,
   boat["theta"])`. *expect:* apparent wind 8.39 m/s at 30.4° off the bow; L = 1036 N, D = 173 N, drive 375 N, side force 981 N;
   closest course to the apparent wind 9.46° + 7.13° = 16.6°; the force balance ("ours: drive = hull drag") gives a steady speed of 3.95 m/s on this course — a little faster than the 3 m/s assumed above — with the apparent wind at 27.4° off the bow, drive 390 N, side force 1207 N; the keel must then work at a lift-to-drag ratio of at least 3.1 and has 8 (`keel_holds` True), and $\beta_a=\varepsilon_A+\varepsilon_H$ = 9.5° + 17.9° closes the polygon (`closure` 0); velocity made
   good 2.12 m/s. *explain:* 5 points.
8. `nb.check_agree` — **from scratch:** `t_hat = np.array([np.cos(th), np.sin(th)])` (course; the wind blows toward −x); `Uw =
   np.array([-U, 0.0])`; `Ur = Uw - V*t_hat`; `ur = Ur/np.hypot(*Ur)`; `n_hat = np.array([np.sin(th), -np.cos(th)])` (to leeward);
   lift direction = `ur` turned by 90° toward the bow; `F = L*lift_dir + D*ur`; `drive = F @ t_hat`; `lateral = F @ n_hat`;
   `assert np.allclose([drive, lateral], [sf["drive"], sf["lateral"]])`.
9. `nb.figure` — `draw_sail_polygon()` with our numbers: velocity triangle (blue), lift (teal), drag (rose), their resultant
   split into drive (purple) and side force, the keel's force closing the polygon. *see / read* ("the resultant leans forward
   of the beam by $\beta_a-\varepsilon_A$") */ change* ("…a baggy sail with L/D = 3: the limit opens to 18.4° + 7.1°").
10. `nb.note` — **N82 [B]** tacking: "No boat sails straight into the wind. It zig-zags at a fixed angle to it, and what counts
    is the speed made good toward the destination, $V\cos\theta$ (ours)." `nb.figure`: a tacking path drawn by our code with the
    no-go wedge shaded; drive from `ch14.sail_drive_curve(theta, …)` against course angle. *see / read / change.*
11. `nb.animation` — **A5** (`player="video"`, 80 frames, FAST 40): (first half) a fish tail through one full stroke with the
    relative velocity, lift and thrust arrows — thrust forward on *both* half-strokes; (second half) a boat tacking upwind,
    the drive arrow forward on both tacks. *see / read / change.*
12. `nb.plotly` — **F9** `slider_figure`: sail L/D from 2 to 12 (21 steps) → drive (per unit lift) against the course angle to
    the true wind, the no-go zone shaded, V/U = 0.5. *see / read* ("a better sail narrows the wedge you cannot sail in") */
    change.*
13. `nb.explainer("thrust_from_lift", heading="How can a boat sail against the wind?", why="Turning the course toward the wind
    and watching the apparent wind swing forward, the lift rotate and the drive shrink to zero is the whole argument;
    switching to the fish tail or the bird's wing shows the same triangle with stream and motion exchanged.", tries=["Head up
    until the drive bar reaches zero and read $\\beta_a$.", "Give the sail a perfect L/D: where is the limit now?", "Switch to
    'fish tail' and press play.", "In 'bird wing', compare the downstroke with the upstroke."])`
14. `nb.md` — **What would change if…** "…the air were compressible? Then pressure itself travels at a finite speed and a wing
    can outrun its own disturbance: Chapter 15."

#### End matter
1. `nb.pointer("**S01** The chapter's exercises are not reproduced. The results the chapter relies on were taught above: the
   pressure loop (C07), the flat plate (C02, C06), the circular arc (C04), leading-edge suction (C07), the half-line vortex
   (C09), the horseshoe wing, its velocity field and ground effect (C09), the parabolic and elliptic loadings (C09, C14), the
   twisted wing (C12), the finite-wing lift slope (C14), the rolling moment (C12), tip-vortex cores (C10), the edge of finite
   angle (C04) and the projectile (C01). **S02** Literature: see the chapter's reference list in the book; the public sources
   used for cross-checks are listed in `analysis/ch14.md` §8.")`
2. `nb.summary(clicked=[16 bullets], feeds_forward=[…], left_out=[…])` — **clicked** (one sentence a reader could say to a
   friend): C01 "A wing is a span, an area and their ratio; in level flight the lift coefficient it needs is its weight over
   dynamic pressure times area." · C02 "Ideal flow allows any circulation; the sharp trailing edge allows one." · C03 "The wing
   gets its circulation by shedding the opposite one — the total never leaves zero." · C04 "Raise the circle and the plate
   bends; move it back and the nose fattens; the tail is the point the map folds." · C05 "Same potential at partner points;
   the speed is rescaled by the map, infinitely at the tail." · C06 "Lift grows by 2π per radian of angle; camber shifts the
   line without tilting it." · C07 "Lift is the area between two pressure curves, mostly suction near the nose." · C08 "Wherever
   the circulation drops along the span, a vortex trails off." · C09 "Each trailing filament is half a line vortex; together
   they make the air at the wing sink." · C10 "Lift is perpendicular to sinking air, so it leans back: drag without friction."
   · C11 "Measure the span by an angle and the loading by sines, and the singular integral is done once." · C12 "N stations, N
   sines, one linear solve: any planform, any twist." · C13 "Lift lives in the first sine; every other one only costs drag." ·
   C14 "One sine means uniform downwash and the least drag, and the price falls with the square of the span." · C15 "Profile
   drag rises with speed, induced drag falls; the best glide is where they are equal." · C16 "Move a wing across the stream
   and its lift leans forward: a fish tail, a wing beat, a sail." **feeds_forward:** Ch. 15 (the lift slope 2π against the
   supersonic thin-airfoil value; the suction peak and the critical Mach number; finite wings with a compressibility-corrected
   K; shock waves and sweep); Ch. 16 (fins, flapping, low Reynolds numbers); wind-turbine and propeller blades (the same
   section in a rotating relative wind). **left_out:** thin-airfoil theory and Theodorsen's method (named in C04; standard
   aerodynamics texts); the Kármán–Trefftz section (C06); viscous prediction of stall and maximum lift (C07, C15 —
   qualitative only); unsteady aerodynamics beyond the lumped-vortex sketch (C03); swept and low-aspect-ratio wings (C12);
   the measured data behind the generic curves (C15).
3. `nb.save()`.

### A.10 Placement check (every curation id has exactly one home)

| Block | CORE | NOTE [B] | NOTE [C] | RECAP | D | Primers | Visuals (≥ 1 per CORE) | From scratch |
|---|---|---|---|---|---|---|---|---|
| §14.1 (before C01) | – | N01 | – | – | – | – | – | – |
| A.2 | C01 | N02 N03 N05 N07 N08 | N04 N06 N83 | R01 R02 | D01 | P336 P337 | planform figure · section figure · F1 (2 slider figures) | trapezoid area, Λ, $\bar c$ |
| A.3 | C02 | N10 | N11 | R04 | – | – | A1 animation · edge-speed figure | ten-line bisection |
| A.3 | C03 | N14 N15 N16 | N17 | – | D02 | P339 | three-stage figure · A2 animation · Γ(t) figure · **E3** | lumped-vortex loop |
| A.4 | C04 | N19–N27 | N18 | R05 R06 R07 | D03 D04 D05 | P340 P341 P342 | four-circle figure · cusp zoom · A3 · F2 · **E1** | circle + map in one line |
| A.5 | C05 | N28 | – | R08 | D06 D07 D08 | – | two-plane figure · F3 | mapped velocity at one point |
| A.5 | C06 | N29 N30 N31 N32 N34 | N33 N35 | R09 | D09 | P338 P343 | F4 · theory-vs-model figure · live · **E2** · **E5** | Γ and $C_L$ in three lines |
| A.5 | C07 | N09 N12 N13 | – | R03 R10 | D10 | – | two-way pressure figure · F5 · separation figure · **E4** | $C_p$ with complex arrays |
| A.6 | C08 | N36 N37 N39 N40 | N69 | R11 | D11 | – | rear-view sketch · staircase figure · A4 | finite-difference −dΓ/dy |
| A.6 | C09 | N41 N60 N65 N66 | N61 | R12 | D12 | P344 | F6 · downwash-along-x figure · sign map | sum over N horseshoes |
| A.6 | C10 | N38 N42 N43 N44 N45 | N67 | – | D13 | P350 | force-triangle figure · **E6** | exact triangle |
| A.6 | C11 | N49 N50 N51 | – | – | D14 D15 | P346 P347 | semicircle sketch · sine shapes · downwash and sign-check figure | quadrature against the series |
| A.6 | C12 | N46 N47 N48 N52 N63 N68 | N64 | – | D16 | P345 | F7 · live | nested-loop matrix, N = 4 |
| A.6 | C13 | N53 N54 | – | – | D17 D18 | P348 | coefficient bars · **E7** | trapezoid on $L=\int\rho U\Gamma\,dy$, $D_i=\int\rho w\Gamma\,dy$ (14.15) |
| A.6 | C14 | N55 N56 N57 N58 N59 | N62 | – | D19 D20 | – | planforms-and-loadings figure · F8a | uniform downwash; one coefficient |
| A.7 | C15 | N70 N72 N73 N74 N75 | N71 N76 | – | D21 | P349 | generic curves · share-vs-speed figure · F8b · **E8** | brute-force best glide |
| §14.8 (before C16) | – | N77 N79 | N78 N80 | – | – | – | velocity triangles | – |
| A.9 | C16 | N81 N82 | – | – | D22 | – | sail polygon · tacking figure · A5 · F9 · **E9** | vectors and dot products |
| End | – | – | – | – | – | – | – | S01 S02 in one `nb.pointer` |

Totals: 16 C · 83 N (66 [B] + 17 [C]) · 12 R · 2 S · 22 D · 15 primers (P336–P350) · 5 animations (A1–A5) · 9 slider-figure groups
(F1–F9, with F8 split as F8a in C14 and F8b in C15) · 2 live cells (C06, C12) · 9 explainers, each embedded once (E1 in C04, E2
and E5 in C06, E3 in C03, E4 in C07, E6 in C10, E7 in C13, E8 in C15, E9 in C16). Blocks without an explainer cell of their own
(C01, C02, C05, C08, C09, C11, C12, C14) each have at least one Python visual of their own, listed above.

---

## Part B — explainer storyboards

Common to all nine (and the backup): created with `tools/new_viz.py`; `<meta name="viz:chapter" content="ch14">`; tabs
Walkthrough · Explore · Explain · Derivation · Equations · Code · Check; every displayed number is computed by a JS function
that mirrors a `ch14` callable and is proved by `selftest()` parity rows (`py:` expressions as C.6; **never a typed constant
the page could compute** — $3\sqrt3/4$ and 2π are computed). Explain is "Explanation & interpretation" in numbered sections built
with `Viz.work.step / line / box / table / hint / interpret`, modelled on `forced_damped_vibrations.html` (the regime-dependent
reading) and `amplitude_phase_second_order_II_3.html` (a numbered derivation with live numbers, one section per optional view):
**0** what the views show and what each colour means · **1…n** every displayed quantity from the controls ("formula = substituted
= result — why", results boxed) · a section or hint per view hidden on phones · the values at the current time or cursor
(live) · **Reading the current setting** (regime-dependent, with the threshold of each regime). Derivation steps are copied from
Part F (same `did` titles, same step count, same order; phones shorten *why* to its first sentence; plain-text *why* and
*watch* never contain raw TeX). Every tour, Explain, Derivation, notes, status, equation and quiz text that names a book
equation **writes it out** next to its number. Colours as header convention 9 (lift teal, drag rose, thrust / drive purple,
circulation amber, pressure / suction orange, stream / velocities blue, ghosts muted). **Γ is labelled "Γ (clockwise)" in every
readout**, and section 0 of every Explain tab that shows a circulation says "the book counts circulation clockwise: a wing
with positive lift has Γ > 0". Angles: sliders and readouts in degrees, arithmetic in radians, both shown in the Explain
lines (trap T12). Walkthrough texts ≤ 45 words, step 1 ≤ 24 words, ≤ 2 extras per step, `play: false` on steps that quote
numbers. **Hygiene:** every power of a formatted number is bracketed; "the book prints", never "the page prints"; slips in
the house form; each control's extreme values are named in its help line; where a formula stops being valid the curve ends
with a mark and a label (no NaN text; an unbounded speed is shown as "unbounded", never as a number). **Local helpers** (not yet
in `assets/viz_lib.js`; copied and listed for promotion): `cx` complex arithmetic (E1, E2, E4 — from ch06 `conformal_joukowski`)
and a small `luSolve` (E6, E7 — from ch06 E8). **Fit:** at 360×640 never three stacked views — two at a time, the third is
`hidePortrait` and its key number is repeated in a visible title or readout; at 844×345 the stage is one row of two views
beside the text panel, titles shortened, transport without step buttons. Parallel builders use private scratch subfolders.
Lengths in E1, E2, E4 are in units of the map constant (b = 1, chord ≈ 4) and speeds in units of U unless a readout says
otherwise.

### E1 · zhukhovsky_airfoil_builder
- **Title:** "Where do camber and thickness come from?" · **Summary:** "Drag the centre of a circle and watch one map turn it
  into a plate, an arc, a fin or a wing section." · **CORE:** C04 (also R05, R06, N24, N25, N26, N27, N07) · **Reference:**
  `amplitude_phase_second_order_II_3.html` (linked windows, numbered Explain with live numbers) and ch06 `conformal_joukowski`
  (two planes, the same marked point).
- **meta:** `viz:order 1` · `viz:sections 14.4` · `viz:equations 14.3 14.4 14.7 14.8` · `viz:fluidpy ch14.zhukhovsky_circle
  ch14.zhukhovsky_circle_from_centre ch14.zhukhovsky_airfoil ch14.zhukhovsky_shape ch14.circular_arc ch14.map_point ch14.joukowski`
  · `viz:derivations D03 D04`.
- **Physics (JS ↔ Python):** `circle(b, e, beta)` → {mu, a} ↔ `ch14.zhukhovsky_circle(b, e, beta)`; `fromCentre(mu, b)` → {e, beta, a} ↔
  `ch14.zhukhovsky_circle_from_centre(mu, b)`; `J(zeta, b)` ↔ `ch14.joukowski(zeta, b)`; `mapPoint(xi, eta, b)` → {R, theta, x, y} ↔
  `ch14.map_point(xi, eta, b)` ($x=(R+b^2/R)\cos\theta$, $y=(R-b^2/R)\sin\theta$ (14.4)); `contour(b, e, beta, n)` ↔
  `ch14.zhukhovsky_airfoil(b, e, beta, n)`; `measure(contour)` → {chord, camber_ratio, thickness_ratio, x_tmax_over_c} (2000 points,
  parabolic refinement of the maxima) ↔ `ch14.zhukhovsky_shape(b, e, beta, exact=True)`; `rules(b, e, beta)` ↔
  `ch14.zhukhovsky_shape(b, e, beta, exact=False)` (coefficient `0.75*Math.sqrt(3)` computed); `arc(b, beta)` ↔ `ch14.circular_arc(b, beta)`.
- **Views** (rows [1.2, 0.8]): 1. `zeta` "ζ-plane: the circle" (row 0, flex 1, `equal`): axes −2.6…2.6; the circle (bold blue)
  through ζ = b, the ghost circle of radius b (muted dashes), the two critical points ±b (amber crosses, "+b: the tail" / "−b:
  must stay inside"), the centre μ as a draggable handle with dashed guides "up = β" and "back = e", the marked point (dot)
  with its radius R from the origin. Pointer: drag the centre (the circle is re-fitted through ζ = b by `fromCentre`); click
  the circle to place the marked point. 2. `z` "z-plane: the section" (row 0, flex 1.4, `equal`): the mapped contour (bold),
  chord line, camber line (teal dashes), the maximum-thickness station (orange tick with its length), the marked point's
  image; with `zoomTE` the last 3 % of chord magnified. 3. `bars` "Rule against measured" (row 1, `hidePortrait`): three pairs of
  bars — chord/4b, camber/chord against β/2, thickness/chord against 1.299 e — with the percentage gap written on each pair.
- **Controls:** `e` "Thickness parameter $e$" (0…0.30, step 0.005, default 0.08; help "0: no thickness · 0.30: far too fat for
  the rules") · `betaDeg` "Camber angle $\beta$" (0…15°, step 0.25, default 4; help "0: symmetric · 15°: strongly cambered") · `t`
  "Marked point" (transport) · `zoomTE` toggle "Zoom on the tail" (optional).
- **Transport:** `t` 0…360° round the circle from the tail, rate 45°/s, `end: 'loop'`.
- **Presets:** "flat plate" {e 0, β 0} · "circular arc" {e 0, β 8} · "symmetric" {e 0.1, β 0} · "cambered" {e 0.1, β 8} · "too fat: the
  rules fail" {e 0.25, β 4}.
- **Status:** gap g = (1.299 e − measured t/c)/measured: e = 0 and β = 0: "flat plate: chord exactly 4b" · g < 5 %: "✅ thin and
  slightly cambered: rules within 5 %" · 5 % ≤ g < 15 %: "≈ thickness rule 8 % high" (number live) · g ≥ 15 %: "⚠️ e = 0.25: thickness
  rule 26 % high — first order in e only".
- **Readouts:** "Centre μ" · "Radius $a$" · "Chord / 4b" · "Camber / chord" · "Thickness / chord".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** · **status** · **inspector** (click
  a point of the circle: its arithmetic).
- **Inspector:** "$\zeta=\xi+i\eta$ = (0.62, 0.95): R = √(0.62² + 0.95²) = 1.134, θ = 56.9° · $x=(R+b^2/R)\cos\theta$ = (1.134 + 0.882) × 0.546 =
  **1.10** · $y=(R-b^2/R)\sin\theta$ = (1.134 − 0.882) × 0.838 = **0.21** (14.4)".
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**Left:** the plane of the circle. The <b class=c-blue>blue circle</b> always passes through the
     amber cross at ζ = b; the dashed ghost is the circle of radius b that would give a flat plate. **Right:** what the map
     $z=\zeta+\dfrac{b^2}{\zeta}$ (14.3) makes of it. The dot is one point and its image."
  1. *Your circle* — "radius $a=b(1+e)/\cos\beta$ = 1 × 1.08/cos 4° = **1.0826** (for β = 0 this is $a\equiv b(1+e)$ (14.7)); centre
     $\mu=b-ae^{-i\beta}$ = **−0.080 + 0.0755 i**: shifted back by be = 0.08 and up by a sin β = 0.0755."
  2. *One point mapped by hand* — the inspector's three lines for the marked point (live).
  3. *The chord* — "nose at ζ = μ − a·(the point opposite the tail) → z = −2.0222; tail at z = 2b = 2: chord = **4.0222**; rule
     $c\simeq4b$ = 4 (gap 0.6 %; the correction is of order e²)."
  4. *The camber* — "measured: the camber line's highest point is **0.139** above the chord line = 3.46 % of the chord; rule
     camber $\simeq\tfrac12\beta c$ = ½ × 0.0698 rad × chord → 3.49 % (gap 0.9 %). The raised circle alone gives an arc of height
     2b tan β = 0.1399."
  5. *The thickness* — "measured: thickest where the upper and lower surface are **0.388** apart = 9.65 % of the chord, at 25.0 %
     chord; rule $t_{\max}/c\simeq1.3\,e$ (14.8) with the coefficient $3\sqrt3/4$ = 1.299: 1.299 × 0.08 = 10.39 % (7.7 % high)."
  6. *The tail* — "at ζ = b the map's rate $\dfrac{dz}{d\zeta}=1-\dfrac{b^2}{\zeta^2}$ is zero: angles double there, so the smooth circle is folded
     into a cusp. Zoom on the tail to see the two surfaces become tangent."
  7. *The bars* (or the hint "turn the phone sideways for the rule-against-measured bars; the three measured numbers are in
     the title of the section view").
  8. *The marked point now* — "t = `Viz.live('t')`° from the tail: ζ = `Viz.live('zeta')` → z = `Viz.live('z')`."
  9. *Reading the current setting* — plate: "The circle is centred on the origin with radius b. Every point and its mirror
     image land on the same point of the real axis: a plate of zero thickness, covered twice." · arc (e = 0, β > 0): "Raising
     the centre keeps both ±b on the circle, so both ends are sharp; the plate has bent into a circular arc of height 2b tan β.
     Still no thickness: upper and lower half of the circle land on the same arc." · symmetric (β = 0, e > 0): "Moving the centre
     back takes −b inside the circle: the front is now rounded. The section is symmetric, thickest a quarter of the way
     back." · cambered thin (g < 5 %): "Both moves at once: β sets the camber, e the thickness, almost independently. The three rules $c\simeq4b$, camber $\simeq\tfrac12\beta c$, $t_{\max}/c\simeq1.3\,e$ (14.8) hold to a few per cent." · fat (g ≥ 15 %): "The rules were derived by keeping only first powers of e. At this
     thickness the neglected e² terms are a quarter of the answer: trust the measured bars, not the rules."
- **Derivation tab:** **D03** (15 steps) `view: 'zeta'` on phones; goal `set` {preset circular arc}; step 2 `live` "x = (R + 1/R) cos θ
  = `…`, y = (R − 1/R) sin θ = `…`" for the marked point; step 7 `watch` "the triangle O–Q–P drawn in the circle view"; step 14
  `live` "centre (0, −2 cot 16°) = (0, −6.97), radius 2/sin 16° = 7.26"; step 15 `set` {t: 90 + β}, `watch` "the top of the arc is
  2b tan β above the chord". **D04** (11 steps) `view: 'z'` on phones; goal `set` {preset symmetric}; step 5 `live` "y ≈ 2 × 0.1 ×
  sin θ(1 − cos θ) = `…` against the mapped `…`"; step 9 `set` {t: 120}, `watch` "the marked point sits at the thickest station";
  step 10 `live` "3√3/4 × 0.1 = 0.1299 against measured 0.1179". Interpret: `s => "With e = " + … + " the first-order rule is " + … +
  " % above the measured thickness."`
- **Code:**
  ```python
  mu, a = ch14.zhukhovsky_circle(1.0, {{e}}, np.deg2rad({{beta}}))     # centre {{mu}}, radius {{a}}
  S = ch14.zhukhovsky_airfoil(1.0, {{e}}, np.deg2rad({{beta}}))        # zeta on the circle, z = zeta + b^2/zeta (14.3)
  sh = ch14.zhukhovsky_shape(1.0, {{e}}, np.deg2rad({{beta}}))         # measured: chord {{c}}, t/c {{tc}}, camber/c {{cc}}
  rule = ch14.zhukhovsky_shape(1.0, {{e}}, np.deg2rad({{beta}}), exact=False)   # (14.8): t/c {{tcr}}, camber/c {{ccr}}
  p = ch14.map_point({{xi}}, {{eta}}, 1.0)                             # (14.4) for the marked point: x {{x}}, y {{y}}
  ```
- **Walkthrough (7 steps):** 1. "One map, four shapes" — "A circle goes in, a wing section comes out. What decides which
  section?" `set` {preset cambered}, `play: true` · 2. "The flat plate" — "Centre the circle on the origin with radius b. Each
  point and its mirror image land together on the real axis." `set` {preset flat plate}, `eq: 'map'` · 3. "Raise it" — "Drag the
  centre up. The circle still passes through both crosses, so both ends stay sharp — and the plate bends." `set` {preset
  circular arc}, `derive: {id: 'D03', step: 14}` · 4. "Move it back" — "Now slide the centre left instead. The left cross falls
  inside, the nose rounds, the tail stays sharp." `set` {preset symmetric}, `eq: 'a'` · 5. "Both" — "Camber from β, thickness
  from e. The bars compare the three rules with what the contour measures." `set` {preset cambered}, `eq: 'rules'`, `code: {id:
  'shape', lines: [3, 4]}` · 6. "One point" — "Click the circle. The card traces that point through the map, number by number."
  `inspect: true`, `eq: 'xy'` · 7. "Your turn" — "Predict how thick the section must be before the thickness rule is 10 % off.
  Then drag e." `controls: ['e', 'betaDeg']`.
- **Equations:** `map` "Zhukhovsky transformation" ref 'Eq. (14.3)' $z=\zeta+\dfrac{b^2}{\zeta}$, live with the marked point · `xy` "Image of a
  point" ref 'Eq. (14.4)' $x=(R+b^2/R)\cos\theta,\ y=(R-b^2/R)\sin\theta$, note "θ is the polar angle of ζ from the origin" · `a` "The shifted
  circle" ref 'Eq. (14.7)' $a\equiv b(1+e),\ e\ll1$, note "ours in general: $a=b(1+e)/\cos\beta$" · `rules` "Shape rules" ref 'Eq. (14.8)'
  $c\simeq4b,\ \text{camber}\simeq\tfrac12\beta c,\ t_{\max}/c\simeq1.3\,e$, live with the three gaps · `arc` "The arc" ref 'unnumbered, §14.4'
  $x^2+(y+2b\cot2\beta)^2=(2b\csc2\beta)^2$.
- **Check yourself:** (1) "Which single point of the circle becomes the sharp tail, and why is it sharp?" — "ζ = b: the map's
  rate $dz/d\zeta=1-b^2/\zeta^2$ vanishes there, so angles double and the smooth circle folds." · (2) "Set e = 0. Does changing β change
  the chord?" — "No: both ends stay at z = ±2b; only the height 2b tan β changes." `set {e: 0, betaDeg: 8}` · (3) "At what fraction
  of the chord is a symmetric section thickest?" — "A quarter: x = −b, from the maximum of sin θ(1 − cos θ) at θ = 120°."
  `set {e: 0.1, betaDeg: 0, t: 120}` · (4) "At e = 0.25, is the rule too high or too low?" — "Too high, by about a quarter: the
  dropped e² terms reduce the thickness." `set {e: 0.25}`.
- **Selftest parity rows:** `{name: 'mu re', js: circle(1, 0.08, rad(4)).mu.re, py: 'np.real(ch14.zhukhovsky_circle(1.0, 0.08,
  np.deg2rad(4.0))[0])', rtol: 1e-12}` (−0.08) · the same with `np.imag` (0.07552) · `{name: 'a', js: circle(1, 0.08, rad(4)).a, py:
  'ch14.zhukhovsky_circle(1.0, 0.08, np.deg2rad(4.0))[1]', rtol: 1e-12}` (1.08264) · `{name: 'inverse e', js: fromCentre(circle(1, 0.12,
  rad(6)).mu, 1).e, py: 'ch14.zhukhovsky_circle_from_centre(ch14.zhukhovsky_circle(1.0, 0.12, np.deg2rad(6.0))[0], 1.0)[0]', rtol:
  1e-10}` · `{name: 'map x', js: mapPoint(0.62, 0.95, 1).x, py: 'ch14.map_point(0.62, 0.95, 1.0)["x"]', rtol: 1e-12}` · `{name: 't/c
  measured', js: measure(contour(1, 0.08, rad(4), 2000)).thickness_ratio, py: 'ch14.zhukhovsky_shape(1.0, 0.08,
  np.deg2rad(4.0))["thickness_ratio"]', rtol: 2e-3}` (0.09645) · `{name: 't/c rule', js: rules(1, 0.08, rad(4)).thickness_ratio, py:
  'ch14.zhukhovsky_shape(1.0, 0.08, np.deg2rad(4.0), exact=False)["thickness_ratio"]', rtol: 1e-12}` (0.10392) · `{name: 'arc radius',
  js: arc(1, rad(8)).radius, py: 'ch14.circular_arc(1.0, np.deg2rad(8.0))["radius"]', rtol: 1e-12}` · invariant `{name: 'tail at 2b',
  js: contour(1, 0.1, rad(8), 400).x[0], expect: 2, rtol: 1e-12}`.
- **Fit plan:** 360×640: status (1 line) · `zeta` (45 %) over `z` (55 %); `bars` hidden — the three measured numbers and their
  gaps sit in the `z` title ("c/4b 1.006 · camber 3.5 % · t/c 9.6 % (rule 10.4 %)"); no drag handle hint text on phones (the
  sliders do the same). 844×345: `zeta` beside `z`, titles "Circle" and "Section". Desktop / 1000×700: rows [1.2, 0.8].

### E2 · kutta_condition
- **Title:** "Which circulation does a wing choose?" · **Summary:** "Slide the circulation and watch the rear stagnation point
  crawl to the sharp tail — the only place where the edge speed stays finite." · **CORE:** C02, C05, C06 (also R04, R09, N28, N29,
  N30, N34, N14) · **Reference:** `forced_damped_vibrations.html` (system + graph on one slider; the Explain panel) and ch06
  `conformal_joukowski`.
- **meta:** `viz:order 2` · `viz:sections 14.3 14.5` · `viz:equations 6.62 14.9 14.10 14.11 14.12` · `viz:fluidpy ch14.circle_flow_zeta
  ch14.kutta_circulation ch14.kutta_circulation_root ch14.stagnation_points ch14.edge_speed ch14.trailing_edge_speed
  ch14.zhukhovsky_lift_coefficient ch14.zhukhovsky_inverse ch14.lift_per_span` · `viz:derivations D07 D09`.
- **Physics (JS ↔ Python):** `dwdzeta(zeta, U, mu, a, alpha, G)` ↔ `ch14.circle_flow_zeta(U, a, Gamma_cw=G, alpha=alpha, center=mu)[1](zeta)`;
  `wOf(zeta, …)` ↔ `…[0](zeta)` (stream function = imaginary part); `invMap(z, b, mu, a)` (root with abs(ζ − μ) ≥ a) ↔
  `ch14.zhukhovsky_inverse(z, b, mu, a)`; `gammaK(U, a, alpha, beta)` ↔ `ch14.kutta_circulation`; `gammaRoot(U, b, e, beta, alpha)` (forty
  bisections on the edge slip) ↔ `ch14.kutta_circulation_root`; `stag(U, b, e, beta, alpha, G)` ↔ `ch14.stagnation_points(…, Gamma_cw=G)`;
  `edgeSpeed(U, b, e, beta, alpha, d, G)` ↔ `ch14.edge_speed(…, d, Gamma_cw=G)`; `teSpeed(U, a, b, alpha, beta)` ↔
  `ch14.trailing_edge_speed`; `cl(alpha, beta)`, `clExact(alpha, beta, e)` ↔ `ch14.zhukhovsky_lift_coefficient`; `lift(rho, U, G)` ↔
  `ch14.lift_per_span(rho, U, Gamma_cw=G)`.
- **Views** (rows [1.25, 1]): 1. `z` "z-plane: the airfoil" (row 0, flex 1.5, `equal`): 16 streamlines (contours of Im w at partner
  points), tracer particles advected by `step`, the body filled, front and rear stagnation points (amber dots), the tail
  region shaded rose when the edge speed is unbounded; the lift arrow (teal, length ∝ Γ) normal to the stream arrow (blue).
  2. `zeta` "ζ-plane: the circle" (row 0, flex 1): the same streamlines and particles round the circle; the centre Q, the
  point B, the stream-aligned diameter (dashed), the angles α, β and α + β drawn at Q; the rear stagnation point on the circle.
  3. `edge` "Speed 0.001 chord from the tail" (row 1, `hidePortrait`): speed/U against Γ/Γ_K on a log axis (the notch at 1),
  faint whole curve + bold cursor; beside it two bars: lift ρUΓ (teal) and the Kutta lift (muted outline).
- **Controls:** `g` "Circulation $\Gamma/\Gamma_K$" (0…1.5, step 0.01, default 0.5; help "0: no circulation, no lift · 1: the Kutta
  value · 1.5: more than a wing ever holds") · `alphaDeg` "Angle of attack $\alpha$" (−6…15°, step 0.5, default 5) · `betaDeg`
  "Camber angle $\beta$" (0…10°, step 0.5, default 4) · `e` "Thickness $e$" (0…0.2, default 0.08; optional) · `flow` toggle
  "Particles" (default on).
- **Transport:** `g` 0…1 "start-up", rate 0.2 per second, `end: 'hold'`; end card: "Γ = Γ_K = 2.13 Ub · the flow leaves the tail
  smoothly at 0.91 U · C_L = 1.06".
- **Presets:** "no circulation" {g 0} · "half way" {g 0.5} · "Kutta" {g 1} · "too much" {g 1.3} · "flat plate" {e 0, β 0, g 1} · "zero-lift
  angle" {alphaDeg −4, betaDeg 4, g 1}.
- **Status:** abs(g − 1) < 0.005: "✅ Kutta: flow leaves the edge smoothly at 0.91 U" · g < 1: "rear stagnation point on the upper
  surface, 1.2 % of chord ahead of the edge — edge speed unbounded" · g > 1: "rear stagnation point on the lower surface, 0.4 %
  of chord ahead of the edge — edge speed unbounded" · α = −β with g = 1: "zero lift: the Kutta circulation is 0". (Percentages
  live from `stag`.)
- **Readouts:** "Γ (clockwise)" · "Kutta value $\Gamma_K$" · "Edge speed" ("unbounded" unless Kutta) · "Lift coefficient" · "Angle
  $\alpha+\beta$".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** · **status** · **terms** ("Lift per
  unit span": ρUΓ now against ρUΓ_K, two bars).
- **Explain:**
  0. *What the views show* — "**Right:** the flow round a circle — known from Chapter 6. **Left:** the same flow carried to the
     airfoil by the map. Amber dots: where the fluid is at rest on the surface. <b class=c-teal>Teal arrow</b>: the lift,
     always at right angles to the <b class=c-blue>stream</b>. The book counts circulation clockwise: a wing with positive
     lift has Γ > 0."
  1. *The circle* — "a = b(1 + e)/cos β = **1.0826**; centre μ = **−0.080 + 0.0755 i**; the tail is the image of B, at ζ = b."
  2. *Where B sits on the circle* — "seen from the centre Q, B lies β = 4° below the ξ-direction; the stream arrives tilted up
     by α = 5°; so measured from the stream direction B is at $\theta_B=-(\alpha+\beta)$ = **−9°** = −0.157 rad."
  3. *The circulation the edge asks for* — "$u_\theta=-2U\sin\theta-\dfrac{\Gamma}{2\pi a}$ (14.10) is zero at B when
     $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11) = 4π × 1 × 1.0826 × sin 9° = **2.128** (in units of Ub). The bisection in the page, which knows
     no formula, finds **2.128**."
  4. *Your circulation* — "Γ = (Γ/Γ_K) × Γ_K = 0.5 × 2.128 = **1.064**; rear stagnation point at sin θ = −Γ/(4πaU) → θ = −4.5°: 4.5°
     of arc short of B."
  5. *The speed at the tail* — "$\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9): the second factor is infinite at ζ = b. Only if the first is
     zero there — the Kutta condition — is the product finite: U(b/a) cos(α + β) = **0.912 U**. With your Γ the speed 0.001 chord
     from the tail is **`live`** U."
  6. *The lift* — "$L=\rho U\Gamma$ (6.62); as a coefficient with the Kutta value, three ways: $C_L\simeq2\pi(\alpha+\beta)$ (14.12) = 2π × 0.157 =
     **0.987**; exact 8π(a/c) sin(α + β) = **1.058**; 2Γ_K/(Uc) = **1.058**."
  7. *The edge-speed graph* (or the hint "turn the phone sideways for the notch; the edge speed is in the title of the
     airfoil view").
  8. *Reading the current setting* — g < 0.995: "Too little circulation. The rear stagnation point is on the upper surface;
     fluid from below must turn the sharp edge to reach it, at a speed that grows without limit toward the edge. A real
     fluid cannot: it sheds a vortex and the circulation rises (see the starting-vortex explainer)." · Kutta: "The stagnation
     point of the circle sits exactly on the point that becomes the tail. The infinite stretching of the map is cancelled
     and the air leaves both surfaces smoothly at 0.91 U. This is the circulation a wing holds." · g > 1.005: "Too much. Now
     fluid from above must turn the edge the other way. Ideal flow allows it; it is not observed." · α = −β: "The stream points
     straight at B. No circulation is needed, and there is no lift: this is the zero-lift angle."
- **Derivation tab:** **D07** (4 steps) `view: 'z'` on phones; goal `set` {g: 0.5}; step 2 `live` "at ζ = 1.05: dz/dζ = 1 − 1/1.1025 =
  0.093"; step 4 `set` {g: 0.5}, `watch` "the rose patch at the tail: the factor ζ²/(ζ² − b²) has no upper bound there". **D09** (10
  steps) `view: 'zeta'` on phones; goal `set` {preset Kutta}; step 3 `watch` "the angle α + β at Q in the circle view"; step 5
  `live` "2U sin(9°) = 0.313 U must equal Γ/(2πa)"; step 6 `set` {g: 1}, `live` "Γ = 4π × 1.0826 × 0.1564 = 2.128"; step 10 `live` "2π ×
  0.157 = 0.987 against the exact 1.058". Interpret: `s => "At α + β = " + … + "° the Kutta circulation is " + … + " Ub and C_L ≈ " + …`
- **Code:**
  ```python
  mu, a = ch14.zhukhovsky_circle(1.0, {{e}}, np.deg2rad({{beta}}))        # circle through zeta = b
  GK = ch14.kutta_circulation(1.0, a, np.deg2rad({{alpha}}), np.deg2rad({{beta}}))   # (14.11): {{GK}}, clockwise
  GR = ch14.kutta_circulation_root(1.0, 1.0, {{e}}, np.deg2rad({{beta}}), np.deg2rad({{alpha}}))  # bisection: {{GR}}
  G = {{g}}*GK                                                            # your circulation: {{G}}
  st = ch14.stagnation_points(1.0, 1.0, {{e}}, np.deg2rad({{beta}}), np.deg2rad({{alpha}}), Gamma_cw=G)  # rear at {{thr}} deg
  q = ch14.edge_speed(1.0, 1.0, {{e}}, np.deg2rad({{beta}}), np.deg2rad({{alpha}}), 0.004, Gamma_cw=G)   # {{q}} U near the tail
  CL = ch14.zhukhovsky_lift_coefficient(np.deg2rad({{alpha}}), np.deg2rad({{beta}}))   # (14.12): {{CL}}
  ```
- **Walkthrough (7 steps):** 1. "Any circulation will do?" — "Ideal flow round this wing has a solution for every Γ. Slide it and
  see what changes." `set` {g: 0}, `play: true` · 2. "No circulation" — "With Γ = 0 the rear stagnation point is on top. Air from
  below whips round the tail." `set` {g: 0}, `highlight: ['readout:edge']` · 3. "The map's weak spot" — "Speeds are rescaled by the
  map, and at the tail the factor has no limit." `eq: 'vel'`, `derive: {id: 'D07', step: 4}` · 4. "Add circulation" — "Clockwise
  swirl pushes the amber dot toward the tail. Press ▶." `set` {g: 0}, `play: true`, `terms: true` · 5. "Exactly one value" — "When the
  dot reaches the tail the edge speed drops to 0.91 U. The circle view shows why: B sits α + β below the stream." `set` {g: 1},
  `eq: 'kutta'`, `derive: {id: 'D09', step: 6}` · 6. "The lift law" — "That circulation gives a lift coefficient of 2π per radian
  of α + β." `eq: 'cl'`, `code: {id: 'kutta', lines: [2, 7]}` · 7. "Your turn" — "Predict the angle of attack at which this cambered
  wing gives no lift. Then find it." `controls: ['alphaDeg', 'betaDeg']`.
- **Equations:** `vel` "Velocity in the two planes" ref 'Eq. (14.9)' $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\dfrac{\zeta^2}{\zeta^2-b^2}$, note "$dw/dz=u-iv$" · `ut` "Speed on
  the circle" ref 'Eq. (14.10)' $u_\theta=-2U\sin\theta-\dfrac{\Gamma}{2\pi a}$, note "θ from the stream-aligned diameter through Q; Γ clockwise" ·
  `kutta` "Kutta circulation" ref 'Eq. (14.11)' $\Gamma=4\pi Ua\sin(\alpha+\beta)$, live · `cl` "Lift coefficient" ref 'Eq. (14.12)'
  $C_L=\dfrac{L}{\tfrac12\rho U^2c}\simeq2\pi(\alpha+\beta)$, live with the exact value · `kz` "Lift theorem" ref 'Eq. (6.62)' $L=\rho U\Gamma$.
- **Check yourself:** (1) "At Γ = 0, on which surface is the rear stagnation point?" — "The upper one (for α + β > 0): without
  swirl the circle's rear stagnation point is on the stream-aligned diameter, above B." `set {g: 0}` · (2) "Double α + β (keep it
  small). What happens to Γ_K?" — "It doubles: Γ_K = 4πUa sin(α + β)." · (3) "At what α does this wing (β = 4°) carry no
  circulation?" — "α = −4°: the stream then points along QB." `set {alphaDeg: -4, g: 1}` · (4) "Why is the exact $C_L$ larger than
  2π(α + β)?" — "Thickness makes the circle's radius a larger than a quarter of the chord."
- **Selftest parity rows:** `{name: 'GK', js: gammaK(1, circle(1, 0.08, rad(4)).a, rad(5), rad(4)), py: 'ch14.kutta_circulation(1.0,
  ch14.zhukhovsky_circle(1.0, 0.08, np.deg2rad(4.0))[1], np.deg2rad(5.0), np.deg2rad(4.0))', rtol: 1e-12}` (2.12826) · `{name: 'root',
  js: gammaRoot(1, 1, 0.08, rad(4), rad(5)), py: 'ch14.kutta_circulation_root(1.0, 1.0, 0.08, np.deg2rad(4.0), np.deg2rad(5.0))',
  rtol: 1e-9}` · `{name: 'rear theta', js: stag(1, 1, 0.08, rad(4), rad(5), 1.0).theta_rear, py: 'ch14.stagnation_points(1.0, 1.0, 0.08,
  np.deg2rad(4.0), np.deg2rad(5.0), Gamma_cw=1.0)["theta_rear"]', rtol: 1e-10}` · `{name: 'TE speed', js: teSpeed(1, 1.08264, 1, rad(5),
  rad(4)), py: 'ch14.trailing_edge_speed(1.0, ch14.zhukhovsky_circle(1.0, 0.08, np.deg2rad(4.0))[1], 1.0, np.deg2rad(5.0),
  np.deg2rad(4.0))', rtol: 1e-6}` (0.9123) · `{name: 'edge speed off-Kutta', js: edgeSpeed(1, 1, 0.08, rad(4), rad(5), 0.004, 1.0), py:
  'ch14.edge_speed(1.0, 1.0, 0.08, np.deg2rad(4.0), np.deg2rad(5.0), 0.004, Gamma_cw=1.0)', rtol: 1e-8}` · `{name: 'inverse map re',
  js: invMap(cx(-1.5, 0.4), 1, mu, a).re, py: 'np.real(ch14.zhukhovsky_inverse(-1.5+0.4j, 1.0, ch14.zhukhovsky_circle(1.0, 0.08,
  np.deg2rad(4.0))[0], ch14.zhukhovsky_circle(1.0, 0.08, np.deg2rad(4.0))[1]))', rtol: 1e-10}` (and `np.imag`) · `{name: 'CL exact', js:
  clExact(rad(5), rad(4), 0.08), py: 'ch14.zhukhovsky_lift_coefficient(np.deg2rad(5.0), np.deg2rad(4.0), 0.08, exact=True)', rtol:
  2e-4}` (1.0583) · invariant `{name: 'edge slip zero at Kutta', js: edgeSlip(GK), expect: 0, atol: 1e-12}`.
- **Fit plan:** 360×640: status (2 lines max) · `z` (58 %) over `zeta` (42 %); `edge` hidden — "edge speed 0.91 U · Γ/Γ_K 1.00" in
  the `z` title; particles reduced to 40. 844×345: `z` beside `zeta`. Desktop: rows [1.25, 1].

### E3 · starting_vortex
- **Title:** "Where does a wing's circulation come from?" · **Summary:** "Start, stop or pitch the foil and watch it shed a vortex
  that exactly cancels what it gains." · **CORE:** C03 (also N14, N15, N16, R11, D02's circuit) · **Reference:**
  `angular_frequency_explorer_1.html` (system animation + graph on one clock, end-of-run summary) and ch13 `geostrophic_adjustment`
  (a conserved total flat in amber).
- **meta:** `viz:order 3` · `viz:sections 14.3` · `viz:equations 5.8 6.62 14.11` · `viz:fluidpy ch14.impulsive_start_lumped_vortex
  ch14.kutta_circulation ch14.flat_plate_forces ch14.lift_per_span` · `viz:derivations D02`.
- **Physics (JS ↔ Python):** `march(U, c, alpha, tEnd, dt)` → {t, Gb[], Gw[], xw[], yw[]} — the algorithm of C.6 item 4, line for line
  ↔ `ch14.impulsive_start_lumped_vortex(U, c, alpha, t_end, dt)` (keys "Gamma_bound_cw", "Gamma_wake_cw", "total", "ratio"); `gammaSteady(U, c, alpha)` = πUc sin α ↔ `ch14.kutta_circulation(U,
  c/4, alpha, 0.0)`; `lift(rho, U, G)` ↔ `ch14.lift_per_span`; `enclosed(circuit, state)` = sum of the strengths inside the
  rectangle (page-only bookkeeping; pinned by an invariant row). Events re-run the march in pieces: *stop* sets U → 0 for the
  foil at $t_1$ (the bound vortex is released as a free vortex); *step in α* and *flap down* change the incidence by +3° at
  $t_1$ (a flap is modelled as extra incidence — said in the notes). U = c = 1; Γ in units of Uc, time in chords travelled.
- **Views** (rows [1.2, 1]): 1. `scene` "The foil and what it leaves behind" (row 0, `equal`-free, x from −2 to 22 chords in the
  frame of the still air, the foil moving left): the plate at its angle, the bound circulation as a teal ring with an arrow
  (clockwise), every wake vortex an amber disc (area ∝ strength, a curved arrow for its sense), the resizable material circuit
  (dashed rectangle with corner labels A B C D and the cut B–D), its enclosed total written on it. Pointer: drag a corner of
  the circuit. 2. `gamma` "Circulation against distance travelled" (row 1, flex 1.3): bound (teal), shed total (rose), their sum
  (amber, flat at 0), the steady value as a dashed ghost; bold up to the transport time. 3. `lift` "Lift / steady lift" (row 1,
  flex 1, `hidePortrait`): ρUΓ_bound over the steady lift against chords travelled, ghost at 1.
- **Controls:** `alphaDeg` "Angle of attack $\alpha$" (2…12°, step 0.5, default 5) · `event` chips "impulsive start / sudden stop /
  step in α / flap down" · `circuit` "Circuit half-length" (1…25 chords, default 25; optional) · `wake` toggle "free wake"
  (default rigid; optional) · transport `s`.
- **Transport:** `s` 0…20 chords, rate 1.5 chords/s, `end: 'hold'`; end card: "after 20 chords: bound 0.97 Γ_K · shed −0.97 Γ_K ·
  total 0.000 · the starting vortex is 20 chords behind".
- **Presets:** the four events (each sets `event` and s = 0) · "small circuit" {circuit 3, s 8}.
- **Status** (live numbers): "after 3.0 chords: bound 0.82 Γ_K · shed −0.82 Γ_K · total 0.000" · sudden stop after $t_1$: "stopped:
  the bound vortex has been shed as a stopping vortex of +0.94 Γ_K" · circuit cutting between foil and wake: "circuit encloses
  the foil only: Γ = 0.94 Γ_K".
- **Readouts:** "Γ bound (clockwise)" · "Γ shed (total)" · "Sum" · "Lift / steady" · "Circuit encloses".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** with an **end-of-run summary** · **presets** ·
  **status** · **notes** (what the model is and is not).
- **Notes:** "This is a lumped-vortex sketch (one bound vortex, one new wake vortex per step) — **not in the book — ours** — not
  a simulation of the boundary layer that really does the shedding. It reproduces the one fact the book states: the bound
  and shed circulations cancel at every instant." With the free wake: "the wake vortices now move each other; the first ones
  wind into a spiral — the starting vortex."
- **Explain:**
  0. *What the views show* — "The foil moves to the left through still air. The <b class=c-teal>teal ring</b> is the
     circulation round it (clockwise, the book's positive sense); each <b class=c-amber>amber disc</b> is a vortex it has shed.
     Below: <b class=c-teal>bound</b>, <b class=c-rose>shed</b> and their <b class=c-amber>sum</b>."
  1. *The circulation the foil is heading for* — "for a flat plate the Kutta value of $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11), with a = c/4
     and β = 0: πUc sin α = π × sin 5° = **0.274 Uc**."
  2. *Right now* — "bound **`live`** Γ_K + shed **`live`** Γ_K = **`live`** (zero to round-off): Kelvin's theorem
     $\dfrac{D\Gamma}{Dt}=0$ (5.8) for the big circuit, which started in air at rest."
  3. *How fast it builds* — "half the steady value almost at once; **0.65** after 1 chord, **0.82** after 3, **0.88** after 5,
     **0.94** after 10; 90 % after about 6 chords. The young starting vortex is close and its downwash at the plate holds the
     circulation down."
  4. *The lift* — "$L=\rho U\Gamma$ (6.62) with the bound circulation: now **`live`** of its steady value." (or the hint "the lift graph
     is hidden on phones; its value is the readout 'Lift / steady'").
  5. *What the circuit encloses* — "the dashed circuit moves with the fluid. As long as it surrounds both the foil and
     everything shed, its circulation is 0. Shrink it until it cuts between them: it now reads **`live`** Γ_K — the foil's
     circulation — and the part you cut off reads the opposite."
  6. *Reading the current setting* — start, s < 1: "The foil has barely moved. A strong counterclockwise vortex is being wound
     up just behind the tail, and the foil already carries half its final circulation." · start, s ≥ 6: "The starting vortex
     is far behind. The bound circulation is within 10 % of the Kutta value and no more vorticity is shed: steady flight." ·
     stop: "The foil stops; its circulation has to go somewhere. It leaves as a clockwise 'stopping vortex', and the pair —
     starting and stopping — drifts off together." · step / flap: "More incidence (or camber) asks for more circulation. The
     difference is shed as a second vortex of the same sense as the first: every change of lift leaves a vortex behind."
- **Derivation tab:** **D02** (5 steps) `view: 'scene'` on phones; goal `set` {event start, s 0, circuit 25}; step 1 `watch` "the
  dashed circuit ABCD at s = 0: nothing inside moves"; step 2 `set` {s: 8}, `live` "Γ_ABCD = `…` (zero)"; step 3 `watch` "the cut
  B–D splits the circuit in two"; step 4 `live` "left part +0.93 Γ_K, right part −0.93 Γ_K"; step 5 `set` {s: 20}. Interpret: `s =>
  "After " + … + " chords the foil holds " + … + " Γ_K and has shed exactly the opposite."`
- **Code:**
  ```python
  run = ch14.impulsive_start_lumped_vortex(1.0, 1.0, np.deg2rad({{alpha}}), {{s}}, 0.05)   # march to {{s}} chords
  Gb = run["Gamma_bound_cw"][-1]         # bound circulation (clockwise): {{Gb}} Uc
  Gw = run["Gamma_wake_cw"][-1]          # everything shed so far: {{Gw}} Uc
  print(Gb + Gw)                         # {{sum}}: Kelvin's theorem (5.8), D(Gamma)/Dt = 0
  print(run["ratio"][-1])                # {{frac}} of the Kutta value pi*U*c*sin(alpha) = {{GK}}
  ```
- **Walkthrough (6 steps):** 1. "From rest" — "A foil at rest has no circulation. Press ▶ and watch what it does to get some."
  `set` {event start, s 0}, `play: true` · 2. "A vortex left behind" — "The amber discs are vorticity shed at the tail. They turn
  counterclockwise; the foil's ring turns clockwise." `set` {s: 2} · 3. "The sum never moves" — "Teal plus rose is the amber line:
  zero, at every instant." `set` {s: 6}, `eq: 'kelvin'`, `derive: {id: 'D02', step: 3}` · 4. "Cut the circuit" — "Shrink the circuit
  until it holds only the foil. Its circulation is now the foil's." `set` {s: 10, circuit: 4}, `controls: ['circuit']` · 5. "Lift
  takes time" — "Half the lift arrives at once, the rest over several chords." `readouts: ['lift']`, `code: {id: 'march', lines:
  [2, 5]}`, `notes: true` · 6. "Your turn" — "Predict which way the new vortex turns when the foil stops. Then choose 'sudden
  stop'." `controls: ['event']`.
- **Equations:** `kelvin` "Kelvin's theorem" ref 'Eq. (5.8)' $\dfrac{D\Gamma}{Dt}=0$, note "for a circuit that moves with the fluid and
  stays in inviscid fluid" · `sum` "Foil + starting vortex" ref '§14.3 (no number)' $\Gamma_{ABCD}=0=\Gamma_{\text{foil}}+\Gamma_{\text{starting}}$, live ·
  `gk` "Steady value" ref 'Eq. (14.11)' $\Gamma=4\pi Ua\sin(\alpha+\beta)$, note "flat plate: πUc sin α" · `kz` "Lift" ref 'Eq. (6.62)' $L=\rho U\Gamma$.
- **Check yourself:** (1) "When the foil stops, which way does the stopping vortex turn?" — "Clockwise, like the bound
  circulation it replaces: the total must stay zero." `set {event: 'stop'}` · (2) "How many chords until the foil has 90 % of its
  circulation?" — "About six." `set {event: 'start', s: 6}` · (3) "Does the circuit's total change when you make it bigger?" —
  "No — as long as it encloses the foil and everything shed, it reads zero." · (4) "Double α. What happens to the strength of
  the starting vortex?" — "It doubles (πUc sin α)."
- **Selftest parity rows:** `{name: 'steady', js: gammaSteady(1, 1, rad(5)), py: 'ch14.kutta_circulation(1.0, 0.25, np.deg2rad(5.0),
  0.0)', rtol: 1e-12}` (0.27381) · `{name: 'bound after 1 chord', js: march(1, 1, rad(5), 1, 0.05).Gb[19], py:
  'ch14.impulsive_start_lumped_vortex(1.0, 1.0, np.deg2rad(5.0), 1.0, 0.05)["Gamma_bound_cw"][19]', rtol: 1e-9}` (0.176976) ·
  `{name: 'bound after 3 chords', js: march(1, 1, rad(5), 3, 0.05).Gb[59], py: 'ch14.impulsive_start_lumped_vortex(1.0, 1.0,
  np.deg2rad(5.0), 3.0, 0.05)["Gamma_bound_cw"][59]', rtol: 1e-9}` (≈ 0.8156 × 0.27381) · `{name: 'wake total', js: march(1, 1, rad(8), 2,
  0.05).GwTotal, py: 'ch14.impulsive_start_lumped_vortex(1.0, 1.0, np.deg2rad(8.0), 2.0, 0.05)["Gamma_wake_cw"][39]', rtol: 1e-9}`
  · invariants `{name: 'Kelvin', js: sumAll(march(1, 1, rad(5), 5, 0.05)), expect: 0, atol: 1e-12}` and `{name: 'big circuit', js:
  enclosed(bigCircuit, state), expect: 0, atol: 1e-12}`.
- **Fit plan:** 360×640: status (1 line) · `scene` (50 %) over `gamma` (50 %); `lift` hidden (readout "Lift / steady"); circuit
  handles are 28 px squares; the wake-model toggle and the circuit slider are optional (hidden at density ≥ 2). 844×345:
  `scene` beside `gamma`. Desktop: rows [1.2, 1].

### E4 · airfoil_pressure_lift
- **Title:** "Is a wing pushed up or pulled up?" · **Summary:** "Turn the angle of attack and watch the suction peak, the
  pressure loop and the lift grow together." · **CORE:** C07 (also R03, N09, N10, N12, N13, R10) · **Reference:**
  `amplitude_phase_second_order_II_3.html` (linked windows with a crosshair readout; numbered Explain).
- **meta:** `viz:order 4` · `viz:sections 14.3 14.5` · `viz:equations 6.62 14.9 14.11 14.12` · `viz:fluidpy ch14.zhukhovsky_surface
  ch14.surface_point ch14.normal_force_coefficient ch14.lift_from_cp ch14.surface_circulation ch14.zhukhovsky_force
  ch14.kutta_circulation ch14.suction_side_separation` · `viz:derivations D10`.
- **Physics (JS ↔ Python):** `surfacePoint(U, b, e, beta, alpha, t)` → {x, y, dwdzeta_abs, map_factor_abs, speed, cp} ↔
  `ch14.surface_point(U, b, e, beta, alpha, t)` (the tail t = −β uses the limit U(b/a) cos(α + β)); `surface(…, n)` → arrays, split
  into upper and lower at the leftmost point ↔ `ch14.zhukhovsky_surface(U, b, e, beta, alpha, n)`; `loopArea(S)` = trapezoid of
  $(C_{p,l}-C_{p,u})$ over x/c ↔ `ch14.normal_force_coefficient(S["upper"]["xc"], S["upper"]["cp"], S["lower"]["xc"], S["lower"]["cp"])`;
  `clFromLoop(Cy, alpha)` = $C_y/\cos\alpha$ ↔ `ch14.lift_from_cp(…, alpha)`; `circ(S)` = ∮u_s ds ↔ `ch14.surface_circulation(S)`;
  `forces(…)` ↔ `ch14.zhukhovsky_force(…, rho, "pressure")`; `sepPoint(alpha, e, beta, Re)` (Thwaites' quadrature on the
  upper-surface speed from the stagnation point; separation where the pressure-gradient parameter reaches the laminar
  limit) ↔ `ch14.suction_side_separation(alpha, e, beta, Re_c)["x_sep_over_c"]` — qualitative, parity to 2 %.
- **Views** (rows [1, 1.1]): 1. `foil` "Pressure on the surface" (row 0, `equal`): the section at its angle to the stream arrow
  (blue); pressure arrows normal to the surface every 3 % of arc — outward where $C_p<0$ (orange, suction), inward where
  $C_p>0$ (blue-grey), length ∝ abs($C_p$) (clipped at 3 with a "×" mark); the front stagnation point (amber dot); the cursor
  point (ring); with `sep` on, a rose tick "laminar separation (model)". Pointer: click the surface to move the cursor. 2. `cp`
  "$-C_p$ against x/c" (row 1, flex 1.5): upper surface (orange), lower (blue-grey), the loop shaded, axis −1…3 upward =
  suction, the cursor as a dot with a crosshair; the level $C_p=1$ marked "stagnation". 3. `bars` "Three ways to the lift" (row
  1, flex 1, `hidePortrait`): $C_y/\cos\alpha$ from the loop (teal), 2Γ/(Uc) (amber), 2π(α + β) (muted outline).
- **Controls:** `alphaDeg` "Angle of attack $\alpha$" (−6…12°, step 0.25, default 5; help "−6: lift reversed · 12: a real section
  would be close to stall") · `betaDeg` "Camber angle $\beta$" (0…8°, default 4) · `e` "Thickness $e$" (0.02…0.20, default 0.08; help
  "0.02: nearly a plate, a needle-sharp peak") · `cur` "Surface cursor" (transport) · `sep` toggle "Laminar separation estimate
  (qualitative)" (optional).
- **Transport:** `cur` 0…1 = fraction of the contour from the tail over the upper surface, round the nose and back, rate 0.12/s,
  `end: 'loop'`.
- **Presets:** "zero lift" {alphaDeg −4, betaDeg 4} · "cruise" {alphaDeg 2} · "high α" {alphaDeg 10} · "thick symmetric" {e 0.15, betaDeg
  0, alphaDeg 5} · "nearly a plate" {e 0.02, betaDeg 0, alphaDeg 5}.
- **Status:** "suction peak −2.24 at 0.5 % chord · pressure recovers over 99.5 % of the upper surface" (live) · α = −β: "zero lift:
  the loop's two halves cancel" · peak below −6: "⚠️ very sharp peak: a real boundary layer would separate — stall" · α + β < 0:
  "lift reversed: suction on the lower surface".
- **Readouts:** "$C_p$ at cursor" · "Speed at cursor" · "Suction peak" · "Loop area $C_y$" · "Lift coefficient".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** (the cursor) · **inspector** · **terms** ("Normal
  force $C_y$" = upper-surface suction $-\int C_{p,u}\,d(x/c)$ + lower-surface pressure $\int C_{p,l}\,d(x/c)$, two bars and their total) ·
  **presets** · **status**.
- **Inspector** (the cursor point): "ζ = μ + a e^{it}, t = 151°: abs($dw/d\zeta$) = 1.62 U (the circle flow) · map factor
  abs($\zeta^2/(\zeta^2-b^2)$) = 0.93 · speed = 1.62 × 0.93 = **1.50 U** · $C_p=1-(1.50)^2$ = **−1.26**".
- **Explain:**
  0. *What the views show* — "**Top:** arrows show the pressure on the skin: <b class=c-orange>orange, pointing away</b> =
     suction (below ambient); grey, pointing in = push. **Bottom:** the same numbers along the chord, plotted upside down so
     that suction is up; the two curves enclose a loop."
  1. *The circulation* — "Kutta value $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11) = 4π × 1.0826 × sin 9° = **2.128 Ub**."
  2. *Speed at the cursor* — the inspector's lines: circle speed × map factor, from
     $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9) (live).
  3. *Pressure at the cursor* — "Bernoulli: $C_p=\dfrac{p-p_\infty}{\tfrac12\rho U^2}=1-\dfrac{q^2}{U^2}$ = 1 − (`live`)² = **`live`**."
  4. *The suction peak* — "lowest pressure $C_p$ = **−2.24** at 0.45 % of the chord on the upper surface: the air there moves at
     √(1 + 2.24) = 1.80 U."
  5. *The loop is the normal force* — "$C_y=\oint C_p\,d(x/c)$ = (lower curve − upper curve) summed along the chord = **1.054**; of this
     **`live`** % is upper-surface suction (the two term bars)."
  6. *From the loop to the lift* — "the loop gives the force across the *chord*; the lift is across the *stream*, α = 5° away.
     With no drag in ideal flow, $C_L=C_y/\cos\alpha$ = 1.054/0.9962 = **1.058** — the same as 2Γ/(Uc) = 2 × 2.128/4.022 = **1.058**
     ($L=\rho U\Gamma$ (6.62)), and 7 % above $C_L\simeq2\pi(\alpha+\beta)$ (14.12) = 0.987."
  7. *The bars* (or the hint "the three lift numbers are in the graph's title on phones").
  8. *The separation estimate* (when on) — "a laminar boundary layer started at the stagnation point would leave the upper
     surface at **`live`** % chord. **Qualitative model:** it uses the ideal pressure and no transition; it shows the trend —
     forward as α rises — not the stall angle."
  9. *Reading the current setting* — α + β ≈ 0: "The stream meets the camber line head-on. Suction above the rear half is
     cancelled by suction below the front half: the loop crosses itself and its net area is zero." · moderate (peak above −3):
     "A modest peak near the nose, then a long gentle recovery. Most of the loop is above the axis: the wing is pulled up
     more than it is pushed." · high α (peak below −3): "The stagnation point has slid under the nose, the air races round the
     leading edge and the peak is tall and thin. Behind it the pressure climbs steeply — the uphill stretch a boundary layer
     cannot always manage." · thin (e ≤ 0.03): "With almost no nose radius the peak is a spike: the leading-edge suction of a
     flat plate, squeezed onto a vanishing area."
- **Derivation tab:** **D10** (7 steps) `view: 'cp'` on phones; goal `set` {preset cruise}; step 2 `watch` "walk the contour
  counterclockwise from the tail: over the top first"; step 4 `set` {cur: 0.25}, `watch` "on the upper surface x runs backward,
  which is where the minus sign comes from"; step 6 `live` "C_y = `…` (lower) − (`…`) (upper) = `…`"; step 7 `live` "C_L = `…`/cos `…`° =
  `…`". Interpret: `s => "At α = " + … + "° the loop area is " + … + ", of which " + … + " % is suction on the upper surface."`
- **Code:**
  ```python
  S = ch14.zhukhovsky_surface(1.0, 1.0, {{e}}, np.deg2rad({{beta}}), np.deg2rad({{alpha}}), 2000)   # speed and Cp round the contour
  print(S["cp_min"], S["xc_cp_min"])                     # suction peak {{cpmin}} at x/c = {{xpk}}
  Cy = ch14.normal_force_coefficient(S["upper"]["xc"], S["upper"]["cp"],
                                     S["lower"]["xc"], S["lower"]["cp"])   # loop area: {{Cy}}
  CL = Cy/np.cos(np.deg2rad({{alpha}}))                  # across the stream, not the chord: {{CL}}
  print(2*S["Gamma_cw"]/(1.0*S["chord_x"]))              # 2*Gamma/(U c) = {{CLg}}: L = rho U Gamma (6.62)
  p = ch14.surface_point(1.0, 1.0, {{e}}, np.deg2rad({{beta}}), np.deg2rad({{alpha}}), {{t}})  # cursor: Cp = {{cp}}
  ```
- **Walkthrough (7 steps):** 1. "Push or pull?" — "Is a wing held up by pressure underneath, or by suction on top? Look at the
  arrows." `set` {preset cruise} · 2. "Speed sets pressure" — "Where the air is fast the pressure is low. Move the cursor and
  read speed and pressure together." `play: true`, `inspect: true` · 3. "The peak" — "Just behind the nose the air is fastest: the
  suction peak. Raise α and watch it grow." `set` {alphaDeg: 8}, `controls: ['alphaDeg']` · 4. "Upside down on purpose" — "Plotted
  with suction upward, the two surfaces enclose a loop." `eq: 'cp'` · 5. "The loop is the force" — "Its area is the force across
  the chord, as a coefficient." `terms: true`, `derive: {id: 'D10', step: 6}` · 6. "Three roads to the lift" — "Loop area over cos α,
  the circulation, and the 2π rule." `eq: 'cl'`, `code: {id: 'cp', lines: [3, 6]}` · 7. "Your turn" — "Predict the angle at which
  the loop's net area is zero. Then find it." `controls: ['alphaDeg', 'betaDeg']`.
- **Equations:** `cp` "Pressure coefficient" ref 'ours (Bernoulli)' $C_p=\dfrac{p-p_\infty}{\tfrac12\rho U^2}=1-\dfrac{\lvert dw/dz\rvert^2}{U^2}$, live at the cursor ·
  `vel` "Velocity on the airfoil" ref 'Eq. (14.9)' $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\dfrac{\zeta^2}{\zeta^2-b^2}$ · `cy` "Normal force = loop area" ref 'D10 (no book
  number)' $C_y\equiv\dfrac{F_y}{\tfrac12\rho U^2c}=\oint C_p\,d\!\Big(\dfrac xc\Big)$, live · `cl` "Lift" ref 'Eq. (6.62), (14.12)' $L=\rho U\Gamma$,
  $C_L\simeq2\pi(\alpha+\beta)$, live.
- **Check yourself:** (1) "What is $C_p$ at the front stagnation point, at any setting?" — "+1: the speed is zero there." ·
  (2) "At α = 5°, which contributes more to the loop — upper suction or lower pressure?" — "Upper suction, by far (read the
  two bars)." `set {alphaDeg: 5}` · (3) "Raise α from 2° to 10°: where does the stagnation point go?" — "Under the nose and
  backward along the lower surface." `set {alphaDeg: 10}` · (4) "Why is the lift not simply the loop area?" — "The loop is the
  force normal to the chord; the lift is normal to the stream, so divide by cos α (the missing part is the leading-edge
  suction)."
- **Selftest parity rows:** `{name: 'cp min', js: surface(1, 1, 0.08, rad(4), rad(5), 2000).cp_min, py: 'ch14.zhukhovsky_surface(1.0,
  1.0, 0.08, np.deg2rad(4.0), np.deg2rad(5.0), 2000)["cp_min"]', rtol: 1e-6}` · `{name: 'point speed', js: surfacePoint(1, 1, 0.08,
  rad(4), rad(5), 2.6).speed, py: 'ch14.surface_point(1.0, 1.0, 0.08, np.deg2rad(4.0), np.deg2rad(5.0), 2.6)["speed"]', rtol: 1e-10}` ·
  `{name: 'point cp', …["cp"], rtol: 1e-10}` · `{name: 'tail speed', js: surfacePoint(1, 1, 0.08, rad(4), rad(5), -rad(4)).speed, py:
  'ch14.trailing_edge_speed(1.0, ch14.zhukhovsky_circle(1.0, 0.08, np.deg2rad(4.0))[1], 1.0, np.deg2rad(5.0), np.deg2rad(4.0))',
  rtol: 1e-9}` (0.9123) · `{name: 'CL pressure route', js: clFromLoop(loopArea(S), rad(5)), py: 'ch14.zhukhovsky_force(1.0, 1.0, 0.08,
  np.deg2rad(4.0), np.deg2rad(5.0), 1.0, "pressure")["CL"]', rtol: 1e-3}` (1.0583) · `{name: 'circulation', js: circ(S), py:
  'ch14.kutta_circulation(1.0, ch14.zhukhovsky_circle(1.0, 0.08, np.deg2rad(4.0))[1], np.deg2rad(5.0), np.deg2rad(4.0))', rtol:
  1e-3}` (2.128) · `{name: 'separation (qualitative)', js: sepPoint(rad(6), 0.08, rad(4), 1e6), py:
  'ch14.suction_side_separation(np.deg2rad(6.0), 0.08, np.deg2rad(4.0), 1.0e6)["x_sep_over_c"]', rtol: 2e-2}` · invariant `{name:
  'stagnation cp', js: maxCp(S), expect: 1, rtol: 1e-3}`.
- **Fit plan:** 360×640: status (1–2 lines) · `foil` (42 %) over `cp` (58 %); `bars` hidden — "C_L: loop 1.058 · Γ 1.058 · 2π rule
  0.987" in the `cp` title; arrows every 5 % of arc on phones. 844×345: `foil` beside `cp`. Desktop: rows [1, 1.1].

### E5 · lift_curve_stall
- **Title:** "Does more angle always mean more lift?" · **Summary:** "Drag the angle of attack along the lift curve: camber
  shifts the line, the slope stays 2π, and stall ends it." · **CORE:** C06 (also N31, N32, N13, N70, N72, N73, N59) · **Reference:**
  `overfitting_curves.html` (a clear two-slider figure with a marker and a verdict line) with a `forced_damped_vibrations.html`
  Explain panel.
- **meta:** `viz:order 5` · `viz:sections 14.5 14.7` · `viz:equations 14.12 14.16` · `viz:fluidpy ch14.zhukhovsky_lift_coefficient
  ch14.zero_lift_angle ch14.generic_lift_curve ch14.stall_regime ch14.stall_speed ch14.finite_wing_lift_slope
  ch14.zhukhovsky_airfoil` · `viz:derivations none`.
- **Physics (JS ↔ Python):** `clLine(alpha, beta)` = 2π(α + β), `clExact(alpha, beta, e)` ↔ `ch14.zhukhovsky_lift_coefficient(alpha, beta,
  e, exact=…)`; `zeroLift(camberRatio)` ↔ `ch14.zero_lift_angle`; `clModel(alpha, K, beta, alphaStall, CLpost, width)` ↔
  `ch14.generic_lift_curve(alpha, K=…, beta=…, alpha_stall=…, CL_post=…, width=…)` (the page copies the module's formula, C.6 item 6: the straight line up to the stall angle, then a smoothstep blend to a post-stall plateau; **qualitative, constants ours**); `regime(Re)` ↔
  `ch14.stall_regime(Re, Re_laminar=5e5, Re_turbulent=3e6)`; `uStall(W, rho, A, CLmax)` ↔ `ch14.stall_speed`; `slope3D(K, AR)` ↔
  `ch14.finite_wing_lift_slope(K, AR)`; `contour(b, e, beta, n)` ↔ `ch14.zhukhovsky_airfoil` (the section sketch). The separated
  streamlines of view 2 are a **drawing** driven by the model's stall fraction, labelled "sketch".
- **Views** (rows [1.2, 1]): 1. `curve` "Lift coefficient against angle of attack" (row 0, flex 1.6): α from −10° to 22°, $C_L$
  from −0.6 to 2.2; the line 2π(α + β) (teal), the exact Zhukhovsky curve (teal dashes), the model of a measured curve (rose,
  legend "model"), the finite-wing line of slope K/(1 + K/πΛ) for Λ = 9.6 (muted ghost); the zero-lift angle −β marked on the
  axis; the stall point marked; the cursor (dot on each curve); the line beyond the stall angle drawn faint ("theory no
  longer applies"). Pointer: drag along α. 2. `section` "The section at this angle" (row 0, flex 1, `equal`): the Zhukhovsky
  section with the current camber (flap = the rear 25 % hinged down), tilted by α to a horizontal stream; attached
  streamlines (blue) when below stall, a growing separated wake (rose hatching, "sketch") beyond it; the lift arrow (teal).
  3. `speed` "Slowest level flight" (row 1, `hidePortrait`): two bars — stall speed of our example aircraft clean and with
  the current flap setting — from $U_{\min}=\sqrt{2W/(\rho AC_{L,\max})}$.
- **Controls:** `alphaDeg` "Angle of attack $\alpha$" (−8…20°, step 0.25, default 5) · `camber` "Camber / chord" (0…0.08, step
  0.005, default 0.035; help "0: symmetric · 0.08: strongly cambered"; β = 2 × camber/chord) · `regime` chips "laminar
  separation first / transition first" (model stall angle 11° / 15°, post-stall level 0.7 / 1.0 — ours) · `flap` toggle "Flap
  down" (adds 0.03 to the camber) · `e` "Thickness $e$" (0.04…0.16, default 0.08; optional).
- **Presets:** "symmetric" {camber 0, alphaDeg 5} · "cambered" {camber 0.035, alphaDeg 5} · "flap down" {flap on, alphaDeg 5} · "zero-lift
  angle" {alphaDeg −4, camber 0.035} · "just past stall" {alphaDeg 16.5, regime transition}.
- **Status:** α < α_stall − 2°: "✅ linear range: the Kutta condition holds · C_L = 0.99 (2π rule)" · within 2° of stall: "≈ near
  maximum lift: the boundary layer is close to separating (model)" · beyond: "⚠️ stalled: lift 35 % below the line — model" ·
  α = −β: "zero lift".
- **Readouts:** "$\alpha+\beta$" · "$C_L$ (2π rule)" · "$C_L$ (exact)" · "$C_L$ (model)" · "Stall speed".
- **Depth features:** Explain + Code · **linked views** (3) · **presets** · **status** · **notes**.
- **Notes:** "Teal is theory: $C_L\simeq2\pi(\alpha+\beta)$ (14.12) and its exact form for this family of sections. Rose is **a drawn
  model** of what a measurement looks like — a straight part with a slope a little under 2π, a maximum, a drop. Its stall
  angle and post-stall level are our own illustrative constants; the book's measured values are not reproduced. Nothing in
  this chapter *predicts* the stall angle." A two-row table "regime · what happens first" with the current row highlighted.
- **Explain:**
  0. *What the views show* — "<b class=c-teal>Teal</b>: ideal-flow theory with the Kutta condition. <b class=c-rose>Rose</b>: a
     model of a measured curve. The dot is your angle of attack. The section on the right is drawn at that angle."
  1. *Camber as an angle* — "β = 2 × camber/chord = 2 × 0.035 = 0.070 rad = **4.0°**; the lift vanishes at $\alpha_{L=0}=-\beta$ = −4.0°."
  2. *The lift from the 2π rule* — "$C_L\simeq2\pi(\alpha+\beta)$ (14.12) = 2π × (5° + 4°) × π/180 = 6.283 × 0.157 = **0.987**. Slope: 2π per
     radian = 0.110 per degree."
  3. *The exact value for this section* — "8π(a/c) sin(α + β) = **1.058** (7 % more: thickness)."
  4. *The model of a measurement* — "below stall K(α + β) with K = 5.7 per radian: **0.896**; the model's maximum is at
     **`live`**°."
  5. *A wing instead of a section* — "on a wing of aspect ratio 9.6 each section works at a reduced angle,
     $C_L=K\Big[\alpha-\dfrac wU+\beta\Big]$ (14.16), and the slope drops to K/(1 + K/πΛ) = **5.20** per radian (the grey line)."
  6. *The slowest flight* — "$U_{\min}=\sqrt{2W/(\rho AC_{L,\max})}$ = √(2 × 9000/(1.225 × 15 × 1.4)) = **26.5 m/s**; with the flap ($C_{L,\max}$
     2.0): **22.1 m/s**." (or the hint "the stall-speed bars are hidden on phones; see the readout").
  7. *Reading the current setting* — linear: "Every extra degree buys 0.11 of lift coefficient, whatever the camber. Camber
     only decides where the line crosses zero." · near stall: "The suction peak (see the pressure explainer) is now so sharp
     that the boundary layer behind it is close to giving up. The measured curve bends away from the line." · stalled: "The
     flow has left the upper surface. The Kutta condition no longer describes it, the circulation collapses, and the teal
     line is fiction here: more angle now means *less* lift." · zero lift: "The stream runs along the line from the circle's
     centre to the tail point; no circulation is needed." · flap: "A flap is camber on demand: the whole line shifts left by
     the extra β, so the same angle gives more lift and the wing can fly slower."
- **Derivation tab:** none (`viz:derivations none`). Walkthrough step 4 quotes the result with its equation and points to
  the `kutta_condition` explainer, where D09 is stepped through.
- **Code:**
  ```python
  beta = 2*{{camber}}                                              # camber angle from camber/chord: {{betadeg}} deg
  CL_rule = ch14.zhukhovsky_lift_coefficient(np.deg2rad({{alpha}}), beta)             # (14.12): {{CLr}}
  CL_exact = ch14.zhukhovsky_lift_coefficient(np.deg2rad({{alpha}}), beta, {{e}}, exact=True)   # {{CLe}}
  CL_model = ch14.generic_lift_curve(np.deg2rad({{alpha}}), K=5.7, beta=beta,
      alpha_stall=np.deg2rad({{astall}}), CL_post={{post}}, width=np.deg2rad(2.0))    # qualitative: {{CLm}}
  U_min = ch14.stall_speed(9000.0, 1.225, 15.0, {{CLmax}})         # slowest level flight: {{Umin}} m/s
  ```
- **Walkthrough (6 steps):** 1. "More angle, more lift?" — "Drag the dot up the curve. Does the lift keep growing?" `controls:
  ['alphaDeg']` · 2. "A straight line" — "In theory the lift coefficient grows by 2π for every radian of angle: 0.11 per degree."
  `set` {camber: 0, alphaDeg: 5}, `eq: 'cl'` · 3. "Camber slides it" — "Add camber. The line moves left; its slope does not
  change. Lift at zero angle is the gift of camber." `set` {camber: 0.035}, `controls: ['camber']`, `eq: 'zero'` · 4. "Where the
  line comes from" — "The Kutta condition fixes the circulation: Γ = 4πUa sin(α + β), and lift is ρUΓ." `eq: 'kutta'`, `code: {id:
  'curve', lines: [2, 3]}` · 5. "Where it ends" — "Past the stall angle the flow separates and the lift drops. The rose curve is a
  model, not a prediction." `set` {alphaDeg: 16.5}, `notes: true` · 6. "Your turn" — "Predict how much slower the aircraft can fly
  with the flap down. Then switch it on." `controls: ['flap', 'regime']`.
- **Equations:** `cl` "Lift coefficient" ref 'Eq. (14.12)' $C_L=\dfrac{L}{\tfrac12\rho U^2c}\simeq2\pi(\alpha+\beta)$, live · `zero` "Zero-lift angle" ref
  '§14.5 (no number)' $\alpha_{L=0}=-\beta,\ \beta=2\,\dfrac{\text{camber}}{c}$ · `kutta` "Kutta circulation" ref 'Eq. (14.11)' $\Gamma=4\pi Ua\sin(\alpha+\beta)$ ·
  `wing` "Section on a wing" ref 'Eq. (14.16)' $C_L=K\Big[\alpha-\dfrac{w(y)}{U}+\beta(y)\Big]$, note "K per radian" · `umin` "Stall speed" ref 'ours'
  $U_{\min}=\sqrt{\dfrac{2W}{\rho A\,C_{L,\max}}}$.
- **Check yourself:** (1) "Double the camber. By how much does the slope of the teal line change?" — "Not at all; the line
  shifts left by the extra β." `set {camber: 0.07}` · (2) "At which angle does a section with 3.5 % camber give no lift?" — "−4°:
  $\alpha_{L=0}=-2\times0.035$ rad." `set {alphaDeg: -4}` · (3) "Is the lift at 18° larger or smaller than at 14° (model)?" — "Smaller:
  past stall the flow has separated." `set {alphaDeg: 18}` · (4) "Why can an aircraft land slower with flaps?" — "More camber
  raises $C_{L,\max}$, and $U_{\min}\propto1/\sqrt{C_{L,\max}}$."
- **Selftest parity rows:** `{name: '2 pi', js: clLine(1, 0), py: 'ch14.zhukhovsky_lift_coefficient(1.0, 0.0)', rtol: 1e-14}` ·
  `{name: 'exact', js: clExact(rad(5), 0.07, 0.08), py: 'ch14.zhukhovsky_lift_coefficient(np.deg2rad(5.0), 0.07, 0.08, exact=True)',
  rtol: 2e-4}` · `{name: 'zero lift', js: zeroLift(0.035), py: 'ch14.zero_lift_angle(0.035)', rtol: 1e-14}` (−0.07) · `{name: 'model
  below stall', js: clModel(rad(4), 5.7, 0.07, rad(14), 0.9, rad(2)), py: 'ch14.generic_lift_curve(np.deg2rad(4.0), K=5.7, beta=0.07,
  alpha_stall=np.deg2rad(14.0), CL_post=0.9, width=np.deg2rad(2.0))', rtol: 1e-10}` · the same at 17° · `{name: 'stall speed', js:
  uStall(9000, 1.225, 15, 1.4), py: 'ch14.stall_speed(9000.0, 1.225, 15.0, 1.4)', rtol: 1e-12}` (26.452) · `{name: 'wing slope', js:
  slope3D(2*Math.PI, 9.6), py: 'ch14.finite_wing_lift_slope(2*np.pi, 9.6)', rtol: 1e-12}` (5.1999).
- **Fit plan:** 360×640: status (1 line) · `curve` (60 %) over `section` (40 %); `speed` hidden (readout "Stall speed"); legend in
  the `curve` title ("teal theory · rose model"). 844×345: `curve` beside `section`. Desktop: rows [1.2, 1].

### E6 · downwash_induced_drag
- **Title:** "How can there be drag without friction?" · **Summary:** "Build a wing from nested horseshoe vortices, add up their
  downwash at one station, and watch the lift lean back." · **CORE:** C08, C09, C10 (also N39, N40, R12, N60, N41, N42–N44, N65, N66)
  · **Reference:** `angular_frequency_explorer_1.html` (system + graph on one cursor, live working) and ch05 `biot_savart_filament`
  ("build the sum").
- **meta:** `viz:order 6` · `viz:sections 14.6` · `viz:equations 5.17 14.13 14.14 14.15` · `viz:fluidpy ch14.staircase_from_loading
  ch14.staircase_downwash ch14.downwash_series ch14.horseshoe_downwash ch14.segment_speed ch14.strip_forces
  ch14.ground_effect_factor ch14.trailing_sheet_strength ch14.solve_lifting_line ch14.induced_angle` · `viz:derivations D11 D12 D13`.
- **Physics (JS ↔ Python):** `loading(shape, y, s)` (constant · elliptic √(1 − (2y/s)²) · parabolic 1 − (2y/s)² · the rectangular-wing
  solution from eight coefficients solved in the page) ; `stair(shape, s, N)` → {edges, steps} ↔ `ch14.staircase_from_loading(Gamma_fn,
  s, N)`; `stairDown(y1, edges, steps)` → {w, parts[]} — one half-line term $\dfrac{\Delta\Gamma_k}{4\pi(y_1-y_k)}$ per trailing filament ↔
  `ch14.staircase_downwash(y1, y_edges, Gamma_steps, parts=True)`; `seriesDown(Gn, gamma1, s)` ↔ `ch14.downwash_series` (the smooth
  limit for the elliptic and rectangular-wing loadings; constant: the horseshoe formula; parabolic: its closed form, checked
  in the notebook); `horseDown(x, y, G, s)` ↔ `ch14.horseshoe_downwash`; `halfLine(d, G)` ↔ `ch14.segment_speed(d, 0.0, np.pi/2, G)`;
  `strip(G, w, U, rho, dy)` ↔ `ch14.strip_forces(…, exact=True)`; `ground(h, s)` ↔ `ch14.ground_effect_factor`; the eight
  coefficients ↔ `ch14.solve_lifting_line(1/6, 1.0, 1.0, 1.0, K=2*np.pi, N=8)["Gn"]`. Units: s = 1, U = 1; the loading level is the
  control `load` = Γ₀/(Us).
- **Views** (rows [1.15, 1]): 1. `wing` "The lifting line and its trailing vortices" (row 0, flex 1.3): a 3-D scene
  (`Viz.three`; drag to orbit) — the bound vortex along y (amber tube, thickness ∝ Γ), N nested horseshoes with their trailing
  legs running downstream (amber, thickness and opacity ∝ strength; the two half-wings with opposite arrows), the station
  cursor (blue ring) with the downwash arrow, and with `h` the ground plane and the mirrored system (muted). Without WebGL:
  the same in an oblique 2-D drawing. 2. `span` "Along the span" (row 1, flex 1.4): Γ(y)/Γ₀ staircase (amber) on its smooth
  limit (faint), w(y)/U at 61 stations (blue) on its smooth limit (faint), the cursor; beyond the tips the upwash shaded
  orange. Pointer: drag the station. 3. `tri` "Forces on the strip at the station" (row 1, flex 1, `hidePortrait`): the
  velocity triangle U, w, $U_e$ (blue; ε exaggerated ×5 and labelled so), the force $dL_e$ normal to $U_e$ split into lift (teal)
  and induced drag (rose).
- **Controls:** `N` chips "1 · 3 · 7 · 21 · ∞" (default 7) · `shape` chips "constant · elliptic · parabolic · rectangular wing"
  (default elliptic) · `y1` "Station $y_1/s$" (transport) · `load` "Loading $\Gamma_0/(Us)$" (0.02…0.12, default 0.07; help "0.02: fast
  cruise · 0.12: slow flight"; optional) · `h` "Height above ground $h/s$" (0.05…2, default 2 = far; optional) · `xs` "Behind
  the wing $x/s$" (−3…3, default 0; optional).
- **Transport:** `y1` from −0.62 to +0.62 (beyond the tips at ±0.5), rate 0.12/s, `end: 'loop'`.
- **Presets:** "one horseshoe" {N 1, shape constant} · "elliptic" {N 21, shape elliptic} · "in ground effect" {h 0.15} · "outboard of
  the tip: upwash" {y1 0.58, N 7} · "the limit" {N ∞, shape elliptic}.
- **Status:** inside the span: "downwash 0.035 U · stream tilted 2.0° down · induced drag 3.5 % of the strip's lift" (live) ·
  outside: "⬆ upwash 0.021 U outboard of the tip — where a following bird flies" · N = 1 near a tip: "⚠️ one horseshoe: the
  downwash grows without limit toward the tip" · h < 0.3: "ground effect: downwash × 0.26".
- **Readouts:** "Γ at station" · "Downwash $w/U$" · "Tilt ε" · "$dD_i/dL$" · "Nearest filament".
- **Depth features:** Explain + Code + Derivation · **3-D view** · **linked views** (3) · **transport** · **terms** ("Downwash at the
  station, filament by filament": the three nearest filaments on each side as signed bars, plus "all the others", summing to
  w) · **presets** · **status**.
- **Explain:**
  0. *What the views show* — "The <b class=c-amber>amber</b> tubes are vortices: one along the wing (bound), the others
     trailing downstream. Where the bound vortex gets weaker, the difference leaves as a trailing leg. <b class=c-blue>Blue</b>:
     the downwash they cause at the wing. Γ is counted clockwise seen from the port tip, so that lift is positive."
  1. *What leaves the wing* — "between your station's two neighbouring panels the bound strength drops from **`live`** to
     **`live`**: a filament of strength $-\dfrac{d\Gamma}{dy}\,dy$ = **`live`** Γ₀ leaves there."
  2. *One filament's push* — "a half-infinite line at distance d induces Γ/(4πd) — half the Γ/(2πd) of an infinite line, from
     the filament law $d\mathbf u=\dfrac{\Gamma\,d\mathbf l}{4\pi}\times\dfrac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ (5.17). The nearest filament (strength `live`, distance `live` s) gives
     **`live`** U."
  3. *The sum* — "all N + 1 filaments: **`live`** U (the term bars); the smooth limit
     $w(y_1)=\dfrac1{4\pi}\displaystyle\int_{-s/2}^{s/2}\dfrac{d\Gamma}{dy}\dfrac{dy}{(y_1-y)}$ (14.13) gives **`live`** U. For the elliptic loading that is Γ₀/2s = 0.035 U at
     every station."
  4. *The tilt* — "ε = tan⁻¹(w/U) = **`live`** rad = **`live`**°; effective angle $\alpha_e=\alpha-\dfrac wU$ (14.14): the strip works at
     `live`° less than the wing's setting."
  5. *The strip's forces* — "force ρ$U_e$Γ dy normal to the tilted stream; across the flight path ρUΓ dy (lift), along it ρwΓ dy
     (induced drag): ratio w/U = **`live`**. Summed along the span these are $L=\int\rho U\Gamma\,dy$ and $D_i=\int\rho w\Gamma\,dy$ (14.15)."
  6. *The force triangle* (or the hint "hidden on phones: ε and dD_i/dL are in the span graph's title").
  7. *With the ground* (h < 2) — "the mirror-image vortices below the ground push *up* at the wing: downwash × [1 − 1/(1 +
     (4h/s)²)] = × **`live`** (formula for one horseshoe)."
  8. *Reading the current setting* — one horseshoe: "The crudest wing: constant circulation, two tip vortices. Fine at
     mid-span (w = Γ/πs), absurd at the tips, where the downwash has no limit." · staircase (3–21): "Each step of the staircase
     sheds a filament. The more steps, the closer the sum comes to the smooth curves — that limit is the integral." ·
     elliptic, N large: "The downwash is the same at every station: no strip works in worse air than any other. This is the
     loading with the least induced drag." · station outside the tips: "Out here every filament pushes the same way — up. A
     wing (or a goose) flying here gets lift for less." · ground: "The wake's images cancel part of the downwash: less tilt,
     less induced drag — the cushion a pilot feels just before touchdown."
- **Derivation tab:** **D11** (4 steps) `view: 'wing'`; goal `set` {N 7, shape elliptic}; step 3 `live` "between the panels at y =
  `…`: Γ drops by `…` Γ₀"; step 4 `watch` "the legs on the two half-wings carry opposite arrows". **D12** (7 steps) `view: 'span'` on
  phones; goal `set` {N 7}; step 2 `live` "nearest filament: Γ/(4πd) = `…`/(4π × `…`) = `…` U", `terms: true`; step 6 `set` {N: 21},
  `watch` "the blue dots settle onto the faint curve"; step 7 `set` {N: '∞'}. **D13** (8 steps) `view: 'span'` on phones (the triangle's
  numbers in its title); goal `set` {shape elliptic, N ∞}; step 1 `live` "ε = atan(`…`) = `…`°"; step 6 `live` "dD_i/dL = w/U = `…`".
  Interpret: `s => "At this station the air sinks at " + … + " U, the lift leans back by " + … + "°, and " + … + " % of it is drag."`
- **Code:**
  ```python
  ye, Gs = ch14.staircase_from_loading(Gamma_fn, 1.0, {{N}})            # {{N}} nested horseshoes: edges and strengths
  w, parts = ch14.staircase_downwash({{y1}}, ye, Gs, parts=True)       # half-line formula summed: w = {{w}}
  w_smooth = ch14.downwash_series(Gn, ch14.span_angle({{y1}}, 1.0), 1.0)   # the limit, (14.13)/(14.19): {{ws}}
  eps = ch14.induced_angle(w, 1.0, small=False)                        # tilt of the stream: {{eps}} rad
  f = ch14.strip_forces({{G}}, w, 1.0, 1.0, 0.01, exact=True)          # dL = {{dL}}, dD_i = {{dDi}}
  ```
- **Walkthrough (8 steps):** 1. "Drag in a perfect fluid?" — "No friction anywhere, and still this wing has drag. Where from?" `set`
  {preset elliptic} · 2. "The simplest wing" — "One vortex along the wing, turning back at each tip: a horseshoe." `set` {preset one
  horseshoe, y1: 0} · 3. "Half a line vortex" — "Each trailing leg pushes the air at the wing down — half as hard as an endless
  vortex would." `terms: true`, `derive: {id: 'D12', step: 1}` · 4. "Many horseshoes" — "A real wing sheds a little vorticity
  everywhere. Step N up and watch the staircase smooth out." `set` {shape elliptic, N: 7}, `controls: ['N']` · 5. "Add them up" —
  "The sum over filaments becomes an integral: the downwash." `set` {N: 21}, `eq: 'w'` · 6. "The lift leans back" — "The strip meets
  air tilted down by ε. Its force, normal to that air, has a backward part." `eq: 'forces'`, `derive: {id: 'D13', step: 6}` · 7.
  "Along the span" — "Press ▶. For the elliptic loading the tilt is the same everywhere." `play: true`, `code: {id: 'down', lines:
  [2, 4]}` · 8. "Your turn" — "Predict the sign of the vertical velocity just outboard of the tip. Then move the station there."
  `controls: ['shape', 'N']`.
- **Equations:** `shed` "Shed filament" ref '§14.6 (no number)' $d\Gamma_{\text{trail}}=-\dfrac{d\Gamma}{dy}\,dy$ · `w` "Downwash" ref 'Eq. (14.13)'
  $w(y_1)=\dfrac1{4\pi}\displaystyle\int_{-s/2}^{s/2}\dfrac{d\Gamma}{dy}\dfrac{dy}{(y_1-y)}$, live, note "a principal value; w positive downward" · `ae` "Effective angle" ref
  'Eq. (14.14)' $\alpha_e=\alpha-\varepsilon=\alpha-\dfrac wU$, live · `forces` "Lift and induced drag" ref 'Eq. (14.15)' $L=\displaystyle\int_{-s/2}^{s/2}\rho U\Gamma\,dy,\ D_i=\displaystyle\int_{-s/2}^{s/2}\rho w\Gamma\,dy$
  · `horse` "One horseshoe" ref 'ours' $w(0)=\dfrac{\Gamma}{\pi s}$.
- **Check yourself:** (1) "With one horseshoe, what is the downwash at mid-span in units of Γ/s?" — "1/π ≈ 0.318: two half-lines
  at distance s/2." `set {N: 1, shape: 'constant', y1: 0}` · (2) "For the elliptic loading, how does the downwash change along
  the span?" — "It does not: Γ₀/2s everywhere." `set {shape: 'elliptic', N: '∞'}` · (3) "Which way does the air move just outside
  a wing tip?" — "Up." `set {y1: 0.58}` · (4) "Halve the downwash. What happens to the induced drag of a strip with the same
  circulation?" — "It halves: dD_i = ρwΓ dy."
- **Selftest parity rows:** `{name: 'half line', js: halfLine(0.5, 1), py: 'ch14.segment_speed(0.5, 0.0, np.pi/2, 1.0)', rtol: 1e-12}`
  (0.15915) · `{name: 'horseshoe centre', js: horseDown(0, 0, 1, 1), py: 'ch14.horseshoe_downwash(0.0, 0.0, 1.0, 1.0)', rtol: 1e-12}`
  (0.31831) · `{name: 'horseshoe outboard', js: horseDown(0, 0.6, 1, 1), py: 'ch14.horseshoe_downwash(0.0, 0.6, 1.0, 1.0)', rtol:
  1e-12}` (negative: an upwash) · `{name: 'staircase N=7', js: stairDown(0, …).w, py: 'ch14.staircase_downwash(0.0,
  ch14.staircase_from_loading(lambda y: np.sqrt(1 - (2*y)**2), 1.0, 7)[0], ch14.staircase_from_loading(lambda y: np.sqrt(1 -
  (2*y)**2), 1.0, 7)[1])', rtol: 1e-10}` (0.49526) · `{name: 'series ellipse', js: seriesDown([1], Math.PI/3, 1), py:
  'ch14.downwash_series(np.array([1.0]), np.pi/3, 1.0)', rtol: 1e-12}` (0.5) · `{name: 'strip drag', js: strip(1, 0.035, 1, 1,
  0.01).dDi, py: 'ch14.strip_forces(1.0, 0.035, 1.0, 1.0, 0.01, exact=True)["dDi"]', rtol: 1e-12}` · `{name: 'ground', js: ground(0.25,
  1), py: 'ch14.ground_effect_factor(0.25, 1.0)', rtol: 1e-12}` (0.5) · `{name: 'rect coefficients', js: GnRect[2]/GnRect[0], py:
  'ch14.solve_lifting_line(1/6, 1.0, 1.0, 1.0, K=2*np.pi, N=8)["Gn"][2]/ch14.solve_lifting_line(1/6, 1.0, 1.0, 1.0, K=2*np.pi,
  N=8)["Gn"][0]', rtol: 1e-8}`. (The one row that passes a `lambda` is allowed by the evaluator: it is a plain Python expression.)
- **Fit plan:** 360×640: status (1 line) · `wing` (45 %) over `span` (55 %); `tri` hidden — "ε 2.0° · dD_i/dL 3.5 %" in the `span`
  title; N chips stay, `load`, `h`, `xs` are optional; the 3-D scene uses the oblique fallback below 420 px width. 844×345:
  `wing` beside `span`. Desktop: rows [1.15, 1].

### E7 · lifting_line_solver
- **Title:** "What makes one wing better than another?" · **Summary:** "Choose a planform and a twist; a small linear solve
  returns the loading as sine amplitudes — the first lifts, the rest only drag." · **CORE:** C11, C12, C13 (also N49, N50, N51, N52,
  N54, N63, N58, N68) · **Reference:** `fid_formula_lab.html` (term-by-term bars with a "make it perfect" preset) with an
  `amplitude_phase_second_order_II_3.html` numbered Explain.
- **meta:** `viz:order 7` · `viz:sections 14.6` · `viz:equations 14.18 14.19 14.20 14.21 14.22` · `viz:fluidpy ch14.planform
  ch14.lifting_line_modes ch14.lifting_line_matrix ch14.lifting_line_rhs ch14.solve_lifting_line ch14.circulation_series
  ch14.downwash_series ch14.glauert_integral ch14.lift_from_series ch14.induced_drag_from_series ch14.span_efficiency
  ch14.elliptic_twist_rectangular` · `viz:derivations D15 D18`.
- **Physics (JS ↔ Python):** `chord(kind, y, s, croot, taper)` ↔ `ch14.planform(kind, s, c_root, taper)["chord"](y)`; `area`, `AR` ↔
  `…["A"]`, `…["AR"]`; `modes(N)` → {n[], gamma[]} ↔ `ch14.lifting_line_modes(N)`; `matrix(cfun, s, K, N)` ↔
  `ch14.lifting_line_matrix(c, s, K, N)[0]`; `rhs(cfun, s, U, alpha, K, N, twistFn)` ↔ `ch14.lifting_line_rhs(c, s, U, alpha, K, N,
  twist=…)`; `luSolve(A, r)` (local, partial pivoting) → Gn ↔ `ch14.solve_lifting_line(…)["Gn"]`; `gammaAt(Gn, gamma)` ↔
  `ch14.circulation_series`; `wAt(Gn, gamma, s)` ↔ `ch14.downwash_series`; `glauert(n, g1)` ↔ `ch14.glauert_integral(n, g1)`;
  `liftOf(Gn, s, rho, U)` ↔ `ch14.lift_from_series`; `dragOf(Gn, rho)` ↔ `ch14.induced_drag_from_series`; `eSpan(Gn)` ↔
  `ch14.span_efficiency`; `ellTwist(alphaR, c, s, K)` ↔ `ch14.elliptic_twist_rectangular(alpha_r, c, s, K=K)`. Units: s = 1, U = 1,
  ρ = 1; the planform's area is fixed by Λ ($A=s^2/\Lambda$), the root chord follows from the taper.
- **Views** (rows [0.8, 1.2]): 1. `plan` "The wing from above" (row 0, `equal`-free, `hidePortrait`): the planform outline, the N
  collocation stations as ticks (crowded toward the tips), the twist as a colour band along the span (blue = nose up, rose =
  nose down, with the angle at root and tip written). Pointer: click a station. 2. `load` "Loading and downwash" (row 1, flex
  1.5): Γ(y)/(Us) (amber, bold) with the elliptic loading of the same lift (muted dashes); w(y)/U (blue) with the uniform
  downwash of that ellipse (muted dashes); the clicked station marked. 3. `bars` "What each sine costs" (row 1, flex 1): bars of
  $n\Gamma_n^2$ for n = 1…N — the first teal ("lift"), the rest rose ("drag only"), percentages on top. Pointer: click a bar to draw
  its shape sin nγ on the `load` view.
- **Controls:** `kind` chips "rectangular · tapered · elliptic" (default rectangular) · `taper` "Taper ratio (tip/root chord)"
  (0.1…1, step 0.05, default 1; active for "tapered"; help "0.1: nearly pointed · 1: rectangular") · `AR` "Aspect ratio $\Lambda$" (3…20,
  step 0.5, default 6) · `twistDeg` "Twist at the tip" (−6…+6°, step 0.5, default 0; linear from root to tip; negative =
  washout; help "−6: tips nosed down · +6: tips nosed up") · `N` "Terms $N$" (1…12, default 8) · `K` "Section slope $K$ per
  radian" (4…6.283, default 6.283; optional) · wing angle fixed at α = 5° (β = 0), said in the Explain tab.
- **Presets:** "elliptic planform" {kind elliptic, twist 0} · "rectangle, Λ = 6" {kind rectangular, AR 6, twist 0} · "taper 0.4" {kind
  tapered, taper 0.4} · "rectangle with elliptic twist — make it perfect" {kind rectangular, `ellTwist` on} · "N = 1" {N 1} · "N = 12"
  {N 12}.
- **Status:** "span efficiency 0.954: 4.8 % more induced drag than the elliptic loading" (live; "4.8 %" = 1/e − 1) · e_span >
  0.999: "✅ elliptic loading: one sine only — the least induced drag for this lift and span" · N ≤ 2: "⚠️ N too small: the
  coefficients have not converged" · cond > 1000: "ill-conditioned: reduce N".
- **Readouts:** "$C_L$" · "$C_{D_i}$" · "Span efficiency" · "Lift slope" · "$\Gamma_3/\Gamma_1$".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **terms** (the bars $n\Gamma_n^2$, total = $\tfrac8{\pi\rho}D_i$;
  click to isolate a shape) · **presets** (with the "make it perfect" one) · **status** · **inspector** (a station).
- **Inspector** (clicked station j): "γ = 67.5°, y = −0.19 s, chord 0.167 s: left side $\dfrac K2Uc(\alpha+\beta)$ = π × 0.167 × 0.0873 =
  **0.0457** · right side $\sum\Big(1+\dfrac{nKc}{4s\sin\gamma}\Big)\Gamma_n\sin(n\gamma)$ = 1.283 × 0.0419 × 0.924 + … = **0.0457** (14.20)".
- **Explain:**
  0. *What the views show* — "<b class=c-amber>Amber</b>: the circulation along the span, which is the lift per metre. Dashed
     grey: the elliptic loading that carries the same total lift. <b class=c-blue>Blue</b>: the downwash. Bars: the loading
     split into sine shapes."
  1. *The stations* — "N = 8 terms need 8 stations: γ_j = jπ/9 = 20°, 40°, … 160°, i.e. $y=-(s/2)\cos\gamma$ = −0.47 s … +0.47 s — never
     a tip, where sin γ = 0."
  2. *One row of the matrix* — the inspector's two lines for the clicked station, from
     $\dfrac K2Uc(\alpha+\beta)=\sum\Big(1+\dfrac{nKc}{4s\sin\gamma}\Big)\Gamma_n\sin(n\gamma)$ (14.20).
  3. *The solved coefficients* — "Γ₁ = **0.0419** Us, Γ₃ = **0.0051**, Γ₅ = **0.0011**, …; the even ones are zero (a symmetric wing).
     $\Gamma=\sum\Gamma_n\sin(n\gamma)$ (14.18)."
  4. *The lift needs only the first* — "$L=\dfrac{\pi s}4\rho U\Gamma_1$ (14.21); as a coefficient $C_L$ = πΓ₁Λ/(2Us) = π × 0.0419 × 6/2 =
     **0.395**. Slope $C_L/\alpha$ = 0.395/0.0873 = **4.53** per radian (a section alone: 6.28)."
  5. *The drag adds every term* — "$D_i=\dfrac\pi8\rho\sum n\Gamma_n^2$ (14.22) = (π/8)(1 × 0.0419² + 3 × 0.0051² + 5 × 0.0011² + …); $C_{D_i}$ =
     **0.00869**; the first term is 95.4 % of it."
  6. *Span efficiency* — "$e_{\text{span}}=\Gamma_1^2/\sum n\Gamma_n^2$ = **0.954**, so $C_{D_i}=C_L^2/(\pi\Lambda e_{\text{span}})$ = 0.395²/(π × 6 × 0.954)."
  7. *The downwash* — "$w(y_1)=\dfrac1{2s}\sum n\Gamma_n\dfrac{\sin(n\gamma_1)}{\sin\gamma_1}$ (14.19, last member): at mid-span **`live`** U, near the tip
     **`live`** U; the ellipse of the same lift would give Γ₁/2s = 0.0210 U everywhere."
  8. *The planform view* (or the hint "hidden on phones: taper, Λ and $e_{\text{span}}$ are in the loading graph's title").
  9. *Reading the current setting* — one bar: "Only the first sine is present. The downwash is uniform, every strip works at
     the same effective angle, and no other loading carries this lift on this span with less induced drag." · rectangular:
     "Too much chord near the tips: the loading is fuller than the ellipse there, the downwash rises toward the tips, and the
     third and fifth sines cost a few per cent." · taper near 0.4: "Almost elliptic without any curved edge: the rose bars
     nearly vanish. This is why so many wings are simply tapered." · strong taper (< 0.25): "Now there is too little chord at
     the tips; the loading is more peaked than the ellipse and the efficiency falls again." · washout: "Nosing the tips down
     unloads them. On a rectangle that moves the loading *toward* the ellipse; on a tapered wing it overshoots." · small N:
     "With one or two terms the wing is forced to be elliptic whatever its shape: the answer is not converged."
- **Derivation tab:** **D15** (11 steps) `view: 'load'`; goal `set` {preset elliptic planform}; step 4 `live` "I₁ = π"; step 6 `live`
  "with γ₁ = 60°: I₂ = 2 cos 60° × π − 0 = π"; step 9 `live` "π sin(3 × 60°)/sin 60° = 0"; step 10 `set` {preset rectangle}, `watch` "the
  blue downwash curve is this sum, term by term"; step 11 `watch` "at the tips each term tends to ±n/(2s) Γ_n — finite". **D18** (9
  steps) `view: 'bars'` (on phones `load`, with the bar percentages in its title); goal `set` {preset rectangle}; step 7 `live`
  "∫ sin 3γ sin 3γ dγ = π/2, ∫ sin γ sin 3γ dγ = 0"; step 8 `terms: true`, `watch` "one bar per n: no bar for a pair n ≠ m"; step 9
  `live` "(π/8)(Γ₁² + 3Γ₃² + 5Γ₅² + …) = `…`". Interpret: `s => "This wing spends " + … + " % of its induced drag on sines that carry no
  lift."`
- **Code:**
  ```python
  P = ch14.planform("{{kind}}", 1.0, {{croot}}, {{taper}})              # chord law; area {{A}}, aspect ratio {{AR}}
  A, gam = ch14.lifting_line_matrix(P, 1.0, {{K}}, {{N}})               # (14.20): one row per station
  r = ch14.lifting_line_rhs(P, 1.0, 1.0, np.deg2rad(5.0), {{K}}, {{N}}, twist=twist)   # left side of (14.20)
  Gn = np.linalg.solve(A, r)                                            # the sine amplitudes: Gamma_1 = {{G1}}
  L = ch14.lift_from_series(Gn, 1.0, 1.0, 1.0)                          # (14.21): first coefficient only -> C_L {{CL}}
  Di = ch14.induced_drag_from_series(Gn, 1.0)                           # (14.22): every coefficient -> C_Di {{CDi}}
  e = ch14.span_efficiency(Gn)                                          # Gamma_1^2 / sum(n Gamma_n^2) = {{e}}
  ```
- **Walkthrough (8 steps):** 1. "A wing, and a question" — "This rectangle carries some lift and pays some induced drag. Could
  another shape do better?" `set` {preset rectangle} · 2. "Sines along the span" — "Measure the span by an angle; write the loading
  as sines. Click a bar to see its shape." `eq: 'series'`, `terms: true` · 3. "Each sine's downwash" — "The hard integral is done
  once per sine — the result is another sine over sin γ." `eq: 'down'`, `derive: {id: 'D15', step: 10}` · 4. "One equation per
  station" — "Each strip obeys the section law at its reduced angle. N stations, N unknowns." `eq: 'll'`, `inspect: true` · 5. "Lift
  from one, drag from all" — "The teal bar is the lift. Every rose bar is drag with no lift." `eq: 'forces'`, `derive: {id: 'D18',
  step: 8}` · 6. "Taper" — "Slide the taper to 0.4. The rose bars almost disappear." `set` {preset taper 0.4}, `controls: ['taper']` ·
  7. "Make it perfect" — "A rectangle can be perfect too — twist it." `set` {preset rectangle with elliptic twist}, `code: {id:
  'solve', lines: [3, 7]}` · 8. "Your turn" — "Predict whether more aspect ratio makes a rectangle closer to or further from
  elliptic. Then drag Λ." `controls: ['AR', 'twistDeg', 'N']`.
- **Equations:** `series` "Sine series" ref 'Eq. (14.18)' $\Gamma=\displaystyle\sum_{n=1}^{\infty}\Gamma_n\sin(n\gamma)$, note "$y=-(s/2)\cos\gamma$" · `down` "Its downwash" ref
  'Eq. (14.19), last member' $w(y_1)=\dfrac1{2s}\displaystyle\sum n\Gamma_n\dfrac{\sin(n\gamma_1)}{\sin\gamma_1}$, note "⚠️ slip #1 — the book prints the two members before
  this one with $\cos\gamma_1-\cos\gamma$; the correct form is $\cos\gamma-\cos\gamma_1$" · `ll` "Lifting-line equation" ref 'Eq. (14.20)'
  $\dfrac K2Uc(\alpha+\beta)=\displaystyle\sum\Big(1+\dfrac{nKc}{4s\sin\gamma}\Big)\Gamma_n\sin(n\gamma)$, live at the clicked station · `forces` "Lift and induced drag" ref 'Eq.
  (14.21), (14.22)' $L=\dfrac{\pi s}4\rho U\Gamma_1,\ D_i=\dfrac\pi8\rho\displaystyle\sum n\Gamma_n^2$, live · `eff` "Span efficiency" ref 'ours'
  $e_{\text{span}}=\Gamma_1^2/\sum n\Gamma_n^2$.
- **Check yourself:** (1) "Which planform leaves exactly one bar?" — "The elliptic one (untwisted): chord ∝ sin γ makes every
  station satisfy the equation with Γ₁ alone." `set {kind: 'elliptic'}` · (2) "For a rectangle of Λ = 6, how much more induced
  drag than the ellipse?" — "About 4.8 % ($e_{\text{span}}$ = 0.954)." · (3) "Which taper ratio is best at Λ = 6?" — "Near 0.4 ($e_{\text{span}}$
  ≈ 0.991)." `set {kind: 'tapered', taper: 0.4}` · (4) "Set N = 1. Why is the efficiency exactly 1?" — "With one term the loading
  *is* a single sine; the truncation, not the wing, made it elliptic."
- **Selftest parity rows:** `{name: 'A00', js: matrix(c1, 6, 2*Math.PI, 3)[0][0], py: 'ch14.lifting_line_matrix(1.0, 6.0, 2*np.pi,
  3)[0][0][0]', rtol: 1e-12}` (0.96891) · `{name: 'A12', js: …[1][2], py: '…[0][1][2]', rtol: 1e-12}` (−1.78540) · `{name: 'rhs', js: rhs(…)[0],
  py: 'ch14.lifting_line_rhs(1.0, 6.0, 1.0, np.deg2rad(5.0), 2*np.pi, 3)[0]', rtol: 1e-12}` (0.27416) · `{name: 'G1 (N=8)', js: Gn[0], py:
  'ch14.solve_lifting_line(1.0, 6.0, 1.0, np.deg2rad(5.0), K=2*np.pi, N=8)["Gn"][0]', rtol: 1e-9}` · `{name: 'CL (N=8)', js: CL, py:
  '…["CL"]', rtol: 1e-9}` (0.39519) · `{name: 'e_span (N=8)', js: eSpan(Gn), py: '…["e_span"]', rtol: 1e-9}` (0.95501) · `{name: 'downwash',
  js: wAt([12, 0, 1], Math.PI/6, 12), py: 'ch14.downwash_series(np.array([12.0, 0.0, 1.0]), np.pi/6, 12.0)', rtol: 1e-12}` (0.75) ·
  `{name: 'Glauert n=2', js: glauert(2, Math.PI/3), py: 'ch14.glauert_integral(2, np.pi/3)', rtol: 1e-12}` (π) · `{name: 'twist', js:
  ellTwist(rad(5), 1, 6, 2*Math.PI).alpha_t, py: 'ch14.elliptic_twist_rectangular(np.deg2rad(5.0), 1.0, 6.0, K=2*np.pi)["alpha_t"]',
  rtol: 1e-12}` (0.02285) · invariant `{name: 'elliptic planform: higher terms', js: maxAbs(GnEll.slice(1)), expect: 0, atol:
  1e-12}`.
- **Fit plan:** 360×640: status (1–2 lines) · `load` (55 %) over `bars` (45 %); `plan` hidden — "taper 1.0 · Λ 6 · e_span 0.954" in the
  `load` title; `K` optional; chips wrap to two rows at most. 844×345: `load` beside `bars`. Desktop / 1000×700: rows [0.8, 1.2].

### E8 · elliptic_loading_aspect_ratio
- **Title:** "Why are glider wings so long and thin?" · **Summary:** "Stretch the span and move the speed: induced drag falls
  with span and speed, profile drag rises with speed, and the best glide is where they meet." · **CORE:** C14, C15 (also N55–N57,
  N59, N73, N74, N75, C01) · **Reference:** `angular_frequency_explorer_1.html` (modes, presets, regime-dependent notes with a
  highlighted table).
- **meta:** `viz:order 8` · `viz:sections 14.6 14.7` · `viz:equations 14.1 14.24 14.25 14.26` · `viz:fluidpy ch14.level_flight
  ch14.induced_drag_coefficient ch14.elliptic_induced_drag ch14.finite_wing_lift_slope ch14.drag_polar ch14.drag_share ch14.best_glide ch14.drag_vs_speed ch14.speed_min_drag ch14.stall_speed ch14.aspect_ratio` · `viz:derivations
  D19 D21`.
- **Physics (JS ↔ Python):** `levelCL(W, rho, U, A)` ↔ `ch14.level_flight(W, rho, U, A)["CL"]`; `cdi(CL, AR, e)` ↔
  `ch14.induced_drag_coefficient(CL, AR, e)`; `polar(CL, CD0, AR, e)` ↔ `ch14.drag_polar(CL, CD0=…, AR=…, e_span=…)`; `dragAt(W, rho, A, AR, CD0, U, e)` → {profile, induced, total, CL, induced_share, LD} (newtons) ↔ `ch14.drag_vs_speed(W, rho, A, AR, CD0, U, e)`; `share(CL, CD0, AR, e)` ↔ `ch14.drag_share(CL, CD0=…, AR=…, e_span=…)`; `bestGlide(CD0, AR, e)` ↔
  `ch14.best_glide`; `uBest(W, rho, A, AR, CD0, e)` ↔ `ch14.speed_min_drag`; `uStall(W, rho, A, CLmax)` ↔ `ch14.stall_speed`;
  `slope3D(K, AR)` ↔ `ch14.finite_wing_lift_slope`; `diEll(L, rho, U, s)` ↔ `ch14.elliptic_induced_drag`; `gamma1(L, s, rho, U)` =
  4L/(πsρU) and `wEll = Γ₁/2s` (from $L=\tfrac{\pi s}4\rho U\Gamma_1$ (14.21) and $w=\Gamma_1/2s$ (14.24); pinned by an invariant row).
- **Views** (rows [0.8, 1.2]): 1. `plan` "The wing" (row 0, `hidePortrait`): the planform at the current Λ (elliptic outline of the
  mode's area or span), drawn to scale beside two muted silhouettes "Λ = 6" and "Λ = 20" (our own outlines); span, area and
  mean chord written under it ($\Lambda\equiv s^2/A=s/\bar c$ (14.1)). 2. `drag` "Drag in level flight against speed" (row 1, flex 1.5):
  profile (blue, rising), induced (rose, falling), total (bold), the minimum marked "best glide", the stall speed as a
  hatched wall on the left, the speed cursor; faint whole curves, bold dot. Pointer: drag the speed. 3. `polar` "The polar"
  (row 1, flex 1, `hidePortrait`): $C_L$ (up) against $C_D$ (across), the tangent from the origin, the cursor's point, the
  dashed parabola of Λ → ∞ (the vertical line $C_D=C_{D0}$).
- **Controls:** `AR` "Aspect ratio $\Lambda$" (3…30, step 0.5, default from the mode) · `hold` chips "same area · same span" (default
  same area) · `U` "Speed $U$" (from 0.8 × stall speed to 3 × best-glide speed of the mode; default: best glide) · `CD0`
  "Zero-lift drag $C_{D0}$" (0.005…0.05, step 0.001) · `eSpan` "Span efficiency $e_{\text{span}}$" (0.6…1, default 0.9; optional) ·
  `sys` (modes).
- **Modes** (`sys`; illustrative inputs of ours, said in the notes): "light aircraft" {W 9000 N, A 15 m², Λ 9.6, C_D0 0.025,
  C_Lmax 1.4, ρ 1.225} (default) · "glider" {W 4500 N, A 11 m², Λ 22, C_D0 0.010, C_Lmax 1.3} · "airliner at altitude" {W 7.0 ×
  10⁵ N, A 125 m², Λ 9.5, C_D0 0.020, C_Lmax 1.5, ρ 0.38 kg/m³} · "albatross" {W 85 N, A 0.62 m², Λ 15, C_D0 0.030, C_Lmax 1.5}.
  (The curation's fourth mode, a wind-turbine blade section, has no "lift = weight" balance and is replaced by the light
  aircraft; the turbine is named in the notes.)
- **Presets:** "best glide" {U = uBest} · "cruise" {U = 1.6 uBest} · "slow climb" {U = 0.85 uBest} · "Λ doubled" {AR × 2} · "infinite span"
  {AR 30, note "the limit Λ → ∞: no induced drag"}.
- **Status:** U < stall: "⚠️ below the stall speed (model limit $C_{L,\max}$): this wing cannot carry the weight" · share within 48–52 %:
  "🎯 best glide: induced drag = profile drag · L/D = 16.5" · share > 52 %: "slow: induced drag is 64 % of the total" · share <
  48 %: "fast: induced drag is only 13 % of the total — profile drag dominates".
- **Readouts:** "$C_L$ needed" · "$C_{D_i}$" · "L/D" · "Induced share" · "Downwash angle".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (four flyers) · **terms** ("Drag in newtons":
  profile + induced = total, with percentages) · **presets** · **status** · **notes** (with a table).
- **Notes:** regime text (as *Reading the current setting*) plus a small table "typical aspect ratios" (our round values:
  delta-wing fighter ≈ 2–3 · light aircraft ≈ 7–10 · airliner ≈ 9–10 · albatross ≈ 15 · glider ≈ 20–30), the row nearest the current
  Λ highlighted; "a wind-turbine or propeller blade is the same section in a rotating relative wind (see `thrust_from_lift`)".
- **Explain:**
  0. *What the views show* — "<b class=c-rose>Rose</b>: induced drag, the price of lift. <b class=c-blue>Blue</b>: profile drag
     (friction and form). Bold: their sum. The dot is your speed."
  1. *The wing* — "$\Lambda\equiv s^2/A$ (14.1): with A = 15 m² and Λ = 9.6 the span is √(9.6 × 15) = **12.0 m** and the mean chord 15/12 =
     **1.25 m**."
  2. *The lift coefficient needed* — "level flight, L = W: $C_L=W/(\tfrac12\rho U^2A)$ = 9000/(0.5 × 1.225 × 34.5² × 15) = **0.824**."
  3. *The price of that lift* — "$C_{D_i}=\dfrac{C_L^2}{\pi\Lambda}$ (14.26), divided by $e_{\text{span}}$ for a loading that is not elliptic: 0.824²/(π ×
     9.6 × 0.9) = **0.0250**. In newtons $D_i=\dfrac{2L^2}{\pi\rho U^2s^2}$ (14.25) over $e_{\text{span}}$ = **273 N** — note: span, not area."
  4. *The other drag* — "profile: $C_{D0}\times\tfrac12\rho U^2A$ = 0.025 × 728 × 15 = **273 N**."
  5. *Total and shares* — "273 + 273 = **546 N**; induced share **50 %**; L/D = 9000/546 = **16.5**."
  6. *The best glide* — "$\Big(\dfrac LD\Big)_{\max}=\tfrac12\sqrt{\dfrac{\pi\Lambda e_{\text{span}}}{C_{D0}}}$ = ½√(π × 9.6 × 0.9/0.025) = **16.5**, at the speed where the two
     drags are equal: **34.5 m/s** (stall speed 26.5 m/s)."
  7. *The sinking air* — "the wing's downwash $w=\Gamma_1/2s$ (14.24) tilts its stream by $C_L/(\pi\Lambda)$ = 0.824/(π × 9.6) = 0.027 rad =
     **1.6°**; the lift slope of this wing is K/(1 + K/πΛ) = **5.20** per radian."
  8. *The polar and the planform* (or the hint "hidden on phones: L/D and the induced share are in the drag graph's title").
  9. *Reading the current setting* — best glide: "The two drags are equal. Slower, and the wing must work harder for the same
     weight; faster, and friction wins. A glider flown at this speed goes furthest." · slow: "Induced drag dominates. This is
     where the book's remark that it is the largest part of the drag holds: take-off, climb, loiter — and every bird that
     soars." · fast: "Profile drag dominates; the induced part falls as 1/U². Airliners cruise here: a published figure is
     about 37 % induced." · large Λ: "The same lift is spread over more span: the wake is wider and weaker, the downwash
     smaller. The rose curve drops, the minimum moves to lower speed and lower drag." · below stall: "The lift coefficient
     needed exceeds what the wing can reach: the curves are drawn faint here."
- **Derivation tab:** **D19** (9 steps) `view: 'drag'` on phones; goal `set` {mode light aircraft, preset best glide}; step 1 `live`
  "Γ₁ = 4 × 9000/(π × 12 × 1.225 × 34.5) = 22.6 m²/s"; step 5 `live` "w = 22.6/24 = 0.94 m/s"; step 7 `live` "2 × 9000²/(π × 1.225 ×
  34.5² × 144) = 246 N (÷ 0.9 = 273 N)"; step 9 `set` {AR: 19.2}, `watch` "double Λ at the same area: the rose curve halves".
  **D21** (8 steps) `view: 'drag'` on phones (`polar` on wide screens); goal `set` {preset best glide}; step 3 `watch` "the tangent
  from the origin touches the polar"; step 4 `live` "C_L = √(π × 9.6 × 0.9 × 0.025) = 0.824"; step 5 `live` "½√(27.1/0.025) = 16.5"; step
  8 `set` {preset cruise}, `terms: true`. Interpret: `s => "At " + … + " m/s induced drag is " + … + " % of the total; the best glide
  is at " + … + " m/s."`
- **Code:**
  ```python
  CL = ch14.level_flight({{W}}, {{rho}}, {{U}}, {{A}})["CL"]                       # L = W: C_L = {{CL}}
  CDi = ch14.induced_drag_coefficient(CL, {{AR}}, {{e}})                           # (14.26)/e_span: {{CDi}}
  d = ch14.drag_vs_speed({{W}}, {{rho}}, {{A}}, {{AR}}, {{CD0}}, {{U}}, {{e}})          # N: {{prof}} + {{ind}} = {{tot}}
  print(d["induced_share"], d["LD"])                                               # {{share}}, {{LD}}
  bg = ch14.best_glide({{CD0}}, {{AR}}, {{e}})                                     # (L/D)max = {{LDmax}} at C_L = {{CLb}}
  Ub = ch14.speed_min_drag({{W}}, {{rho}}, {{A}}, {{AR}}, {{CD0}}, {{e}})           # best-glide speed: {{Ub}} m/s
  ```
- **Walkthrough (7 steps):** 1. "Two wings" — "A glider's wing is long and thin, a fighter's short and broad. What does the
  shape buy?" `set` {mode glider} · 2. "The price of lift" — "Lift on a finite wing costs induced drag: lift coefficient squared
  over π times aspect ratio." `set` {mode light aircraft}, `eq: 'cdi'`, `derive: {id: 'D19', step: 9}` · 3. "Fly slower" — "Slower
  flight needs more lift coefficient, so the rose curve climbs steeply to the left." `set` {preset slow climb}, `terms: true` · 4.
  "Fly faster" — "Friction grows with speed squared. The blue curve wins on the right." `set` {preset cruise} · 5. "The best speed"
  — "The total is lowest where rose equals blue: the best glide." `set` {preset best glide}, `eq: 'ld'`, `derive: {id: 'D21', step:
  4}` · 6. "Stretch the wing" — "Double the aspect ratio at the same area. The minimum drops and moves left." `set` {preset Λ
  doubled}, `code: {id: 'polar', lines: [3, 6]}` · 7. "Your turn" — "Predict the induced share at twice the best-glide speed. Then
  drag the speed." `controls: ['U', 'AR', 'CD0']`.
- **Equations:** `ar` "Aspect ratio" ref 'Eq. (14.1)' $\Lambda\equiv s^2/A=s/\bar c$ · `di` "Induced drag, elliptic loading" ref 'Eq. (14.25)'
  $D_i=\dfrac\pi8\rho\Gamma_1^2=\dfrac{2L^2}{\pi\rho U^2s^2}$, live · `cdi` "Induced drag coefficient" ref 'Eq. (14.26)' $C_{D_i}=\dfrac{C_L^2}{\pi\Lambda}$, live, note "divide by
  $e_{\text{span}}$ for other loadings (ours)" · `ld` "Polar and best glide" ref 'ours' $C_D=C_{D0}+\dfrac{C_L^2}{\pi\Lambda e_{\text{span}}}$,
  $\Big(\dfrac LD\Big)_{\max}=\tfrac12\sqrt{\dfrac{\pi\Lambda e_{\text{span}}}{C_{D0}}}$, live · `w` "Uniform downwash" ref 'Eq. (14.24)' $w(y)=\Gamma_1/2s$.
- **Check yourself:** (1) "At the best-glide speed, what fraction of the drag is induced?" — "Exactly half." `set {preset best
  glide}` · (2) "Double the span at fixed weight and speed. What happens to the induced drag?" — "It falls to a quarter:
  $D_i\propto1/s^2$." `set {hold: 'area', AR: 7.5}` — then raise Λ to 30: four times the aspect ratio at the same area is twice the span · (3) "Is induced drag the
  largest part at cruise?" — "No — for the light aircraft at 55 m/s it is 13 %." `set {U: 55}` · (4) "Halve $C_{D0}$. Does the
  best-glide speed go up or down?" — "Down: the crossing of the two curves moves left."
- **Selftest parity rows:** `{name: 'CL', js: levelCL(9000, 1.225, 30, 15), py: 'ch14.level_flight(9000.0, 1.225, 30.0, 15.0)["CL"]',
  rtol: 1e-12}` (1.0884) · `{name: 'CDi', js: cdi(1.0884, 9.6, 0.9), py: 'ch14.induced_drag_coefficient(1.0884, 9.6, 0.9)', rtol:
  1e-12}` (0.043643) · `{name: 'share 30', js: dragAt(9000, 1.225, 15, 9.6, 0.025, 30, 0.9).induced_share, py: 'ch14.drag_vs_speed(9000.0, 1.225, 15.0, 9.6, 0.025, 30.0, 0.9)["induced_share"]', rtol: 1e-10}` (0.6358) · `{name: 'total
  55', …["total"], rtol: 1e-10}` (802.2 N) · `{name: 'LD max', js: bestGlide(0.025, 9.6, 0.9).LD_max, py: 'ch14.best_glide(0.025, 9.6,
  0.9)["LD_max"]', rtol: 1e-12}` (16.475) · `{name: 'U best', js: uBest(9000, 1.225, 15, 9.6, 0.025, 0.9), py:
  'ch14.speed_min_drag(9000.0, 1.225, 15.0, 9.6, 0.025, 0.9)', rtol: 1e-10}` (34.484) · `{name: 'Di elliptic', js: diEll(9000, 1.225,
  30, 12), py: 'ch14.elliptic_induced_drag(9000.0, 1.225, 30.0, 12.0)', rtol: 1e-12}` (324.8) · `{name: 'slope', js: slope3D(2*Math.PI,
  9.6), py: 'ch14.finite_wing_lift_slope(2*np.pi, 9.6)', rtol: 1e-12}` · invariant `{name: 'half at best glide', js:
  dragAt(…, uBest(…), …).induced_share, expect: 0.5, rtol: 1e-10}`.
- **Fit plan:** 360×640: status (1 line) · mode chips (one row, short labels "light · glider · airliner · albatross") · `drag`
  alone at full stage height with the term bars under it in the text panel; `plan` and `polar` hidden — "L/D 16.5 · induced
  50 % · Λ 9.6 · span 12 m" in the `drag` title. 844×345: `drag` beside `polar`. Desktop: rows [0.8, 1.2].

### E9 · thrust_from_lift
- **Title:** "How can a boat sail against the wind?" · **Summary:** "Turn the course toward the wind and watch the apparent wind
  swing forward and the sail's lift tilt — then see the same triangle on a fish tail and a bird's wing." · **CORE:** C16 (also N77,
  N79, N81, N82, N08) · **Reference:** `angular_frequency_explorer_1.html` (the same idea in several systems, presets, live working)
  and ch13 `ekman_force_balance` (a force polygon that closes).
- **meta:** `viz:order 9` · `viz:sections 14.8 14.9` · `viz:equations none` (the section has no numbered equation; every relation
  is shown in full) · `viz:fluidpy ch14.apparent_wind ch14.sail_forces ch14.sail_drive_curve ch14.sailboat_course_limit
  ch14.sailboat_equilibrium ch14.velocity_made_good ch14.oscillating_foil ch14.resolve_lift_drag` · `viz:derivations D22`.
- **Physics (JS ↔ Python):** `apparent(U, V, theta)` → {Ur, beta_a} ↔ `ch14.apparent_wind(U, V, theta)`; `sail(U, V, theta, CL, CD, rho, A)` → {L, D, drive, lateral, drive_over_L} ↔ `ch14.sail_forces(…)`; `limit(LDs, LDh)` ↔ `ch14.sailboat_course_limit`; `equil(U, theta, CL, CD, rhoAir, Asail, rhoWater, CDAhull)` → {V, beta_a, drive, lateral, LD_hull_required, vmg} (`Viz.num.brentq` on drive − hull drag; hull drag area 0.05 m², water 1000 kg/m³; drawn as the marker "steady speed for this sail and hull" on the curve view) ↔ `ch14.sailboat_equilibrium(U, theta, CL=…, CD=…, rho_air=…, A_sail=…, rho_water=…, CDA_hull=…)`; `vmg(V, theta)` ↔ `ch14.velocity_made_good`; `foil(U, V, L, D)` → {theta,
  thrust, side_force, power_in, power_out, efficiency} ↔ `ch14.oscillating_foil(U, V, L, D)`. In the page the sail's $C_D$ is
  $C_L$/(L/D) with $C_L$ = 1.2, A = 20 m², ρ = 1.225, true wind 6 m/s; fish and bird use L = 10 N and D = L/(L/D).
- **Views** (rows [1.25, 1]): 1. `plan` "From above" (row 0, flex 1.4, `equal`): **sail mode** — the boat on its course (bow up
  the screen is *not* used: the true wind blows down the screen and the boat is rotated by θ), true wind U (blue), minus
  boat velocity (blue, thin), their sum the apparent wind $\mathbf U_r=\mathbf U-\mathbf V$ (blue, bold) with the angle $\beta_a$ at the bow; the
  sail (a cambered line trimmed to the apparent wind); lift L (teal, ⟂ $\mathbf U_r$), drag D (rose, ∥ $\mathbf U_r$), the resultant split into
  drive (purple, along the course) and side force (muted), the keel's opposing force (grey, closing the polygon); the no-go
  wedge shaded. **Fish / bird mode** — the foil, the stream U, the stroke velocity V, the relative velocity, L ⟂ it, its
  thrust part (purple) and side part. 2. `curve` "Drive and progress against course" (row 1, flex 1.3): sail mode — drive/L
  (purple) and velocity made good/U (teal) against the course angle θ from 0° to 180°, the no-go zone shaded rose, the cursor;
  fish / bird — thrust/L and efficiency against V/U. Pointer: drag the course (or V/U). 3. `power` "Bars" (row 1, flex 1,
  `hidePortrait`): sail — drive against side force; fish / bird — power in ($V\times$ side force) against useful power out ($U\times$
  thrust), with the efficiency.
- **Controls:** `sys` (modes: "sail · fish tail · bird wing") · `thetaDeg` "Course angle to the true wind $\theta$" (0…180°, step 1,
  default 45; sail; help "0: straight into the wind · 90: beam reach · 180: running") · `LD` "Lift-to-drag ratio of the sail /
  foil" (1…20, step 0.5, default 6; help "1: a bag · 20: a rigid wing sail") · `VU` "Boat speed / wind speed" or "Stroke speed /
  swimming speed $V/U$" (0…1.5, step 0.05, default 0.5 for the sail, 0.75 for fish and bird) · `LDh` "Lift-to-drag ratio of hull
  and keel" (2…20, default 8; optional; sail) · `stroke` chips "downstroke · upstroke" (bird).
- **Transport** (fish and bird): `ph` stroke phase 0…1, rate 0.5/s, `end: 'loop'` — the tail sweeps across and back; V changes
  sign at mid-stroke, the lift's side part flips, the thrust arrow stays forward and pulses.
- **Presets:** "beam reach" {thetaDeg 90} · "close-hauled" {thetaDeg 45} · "in the no-go zone" {thetaDeg 12, LD 6, VU 0.5} · "a perfect
  foil, D = 0" {LD 20 → the control's maximum, labelled "nearly perfect"} · "a poor sail" {LD 2}.
- **Status** (sail): drive > 0: "⛵ close-hauled 45° off the true wind · apparent wind 30° off the bow · drive positive · making
  0.35 U to windward" · $\varepsilon_A\le\beta_a<\varepsilon_A+\varepsilon_H$: "⚠️ too close to hold: the drive is positive but smaller than the hull's drag at this side
  force — the boat slows" · $\beta_a<\varepsilon_A$ ($\tan\beta_a<D/L$): "⚠️ in irons: the apparent wind is only 8° off the bow, inside the
  sail's drag angle 9.5° — the drive is negative". Fish / bird: "thrust 6.0 N · efficiency 0.81".
- **Readouts:** "Apparent wind" · "Angle $\beta_a$" · "Drive / L" · "Side force / L" · "Speed made good".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (three systems) · **terms** ("Drive along the
  course" = lift part $L\sin\beta_a$ − drag part $D\cos\beta_a$, signed bars) · **presets** · **status** · **transport** (fish / bird).
- **Explain:**
  0. *What the views show* — "<b class=c-blue>Blue</b> arrows are velocities: the true wind, the boat's own motion, and what
     the boat feels. <b class=c-teal>Teal</b> = lift, at right angles to the felt wind; <b class=c-rose>rose</b> = drag, along
     it; <b class=c-accent>purple</b> = the part that drives the boat."
  1. *The wind the boat feels* — "$\mathbf U_r=\mathbf U-\mathbf V$. Along the course: U cos θ + V = 6 × 0.707 + 3 = **7.24 m/s** from ahead;
     across it: U sin θ = **4.24 m/s**. Speed √(7.24² + 4.24²) = **8.39 m/s**."
  2. *Its angle off the bow* — "tan $\beta_a$ = 4.24/7.24 → $\beta_a$ = **30.4°** — closer to the bow than the true wind's 45°."
  3. *The sail's forces* — "L = $C_L\tfrac12\rho U_r^2A$ = 1.2 × 0.5 × 1.225 × 8.39² × 20 = **1036 N**, at right angles to the apparent wind;
     D = L/6 = **173 N**, along it."
  4. *Along and across the course* — "drive = $L\sin\beta_a-D\cos\beta_a$ = 1036 × 0.505 − 173 × 0.863 = 524 − 149 = **375 N**; side force =
     $L\cos\beta_a+D\sin\beta_a$ = 894 + 87 = **981 N**, carried by the keel."
  5. *How close to the wind* — "the drive vanishes when tan $\beta_a$ = D/L, i.e. $\beta_a=\varepsilon_A$ = tan⁻¹(1/6) = **9.5°**. The hull has
     its own drag angle $\varepsilon_H$ = tan⁻¹(1/8) = **7.1°**; steady sailing needs $\beta_a=\varepsilon_A+\varepsilon_H$ = **16.6°** off the apparent wind."
  6. *Progress to windward* — "velocity made good = V cos θ = 3 × 0.707 = **2.12 m/s**."
  7. *The bars* (or the hint "hidden on phones: drive/L and the limiting angle are in the plan view's title").
  8. *Fish tail and bird wing* — "the tail moves sideways at V while the fish advances at U: the relative flow comes at θ with
     tan θ = V/U = 0.75. Thrust = L sin θ = 10 × 0.6 = **6 N**; side force = L cos θ = **8 N**; power out U × 6 = **6 W** = power in
     V × 8 = **6 W** — equal because V = U tan θ. With drag (L/D = 10): thrust 5.2 N, efficiency **0.81**."
  9. *Reading the current setting* — sailing, drive > 0: "The apparent wind is far enough off the bow for the lift's forward
     part to beat the drag's backward part. The keel, a wing in water, turns the large side force into almost no sideways
     motion." · at the limit: "Lift and drag pull equally along the course. Nothing is left to drive the boat: this is as
     close to the wind as this sail can point." · no-go: "Too close to the wind. The sail's total force points behind the
     beam; the boat slows, and must bear away. To get upwind it sails a zig-zag — tacking." · running (θ > 150°): "Now the
     old picture is right: the boat is mostly pushed. Note that it can never go faster than the wind this way, while on a
     reach it can." · fish: "The same triangle with the roles exchanged: the *foil* moves across the stream, and its lift
     leans forward on both half-strokes." · bird: "Downstroke: lift up and forward. Upstroke: the wing is feathered and
     carries little; the cycle's average is thrust and weight support."
- **Derivation tab:** **D22** (8 steps) `view: 'plan'`; goal `set` {mode sail, preset close-hauled}; step 2 `live` "along 6 cos 45° +
  3 = 7.24, across 6 sin 45° = 4.24 m/s"; step 3 `live` "tan β_a = 4.24/7.24, β_a = 30.4°"; step 5 `live` "1036 × 0.505 − 173 × 0.863 =
  375 N", `terms: true`; step 7 `set` {preset in the no-go zone}, `watch` "the purple arrow points backward"; step 8 `live` "9.5° +
  7.1° = 16.6°". Interpret: `s => "On this course the boat feels the wind " + … + "° off the bow and gets " + … + " % of the sail's
  lift as drive."`
- **Code:**
  ```python
  aw = ch14.apparent_wind(6.0, {{V}}, np.deg2rad({{theta}}))          # U_r = U - V: {{Ur}} m/s at {{ba}} deg off the bow
  sf = ch14.sail_forces(6.0, {{V}}, np.deg2rad({{theta}}), CL=1.2, CD={{CD}}, rho=1.225, A=20.0)
  print(sf["drive"], sf["lateral"])                                   # {{drive}} N along the course, {{lat}} N across
  lim = ch14.sailboat_course_limit({{LD}}, {{LDh}})                    # closest course to the apparent wind: {{lim}} deg
  vmg = ch14.velocity_made_good({{V}}, np.deg2rad({{theta}}))          # progress to windward: {{vmg}} m/s
  fish = ch14.oscillating_foil(1.0, {{VUf}}, 10.0, {{Df}})             # tail: thrust {{th}} N, efficiency {{eta}}
  ```
- **Walkthrough (7 steps):** 1. "Against the wind?" — "The wind blows down the screen. The boat is sailing up it at 45°. What
  pushes it?" `set` {mode sail, preset close-hauled} · 2. "The wind it feels" — "Subtract the boat's velocity from the wind's. The
  result comes from further ahead." `eq: 'app'`, `derive: {id: 'D22', step: 3}` · 3. "A wing on end" — "Lift is at right angles to
  the felt wind, drag along it." `highlight: ['readout:ba']` · 4. "Forward part" — "Split both along the course. Lift pulls
  forward, drag pulls back." `terms: true`, `eq: 'drive'` · 5. "Too close" — "Head up. The felt wind swings to the bow and the
  drive dies." `set` {preset in the no-go zone}, `controls: ['thetaDeg']` · 6. "Same triangle, other owners" — "A fish tail moves
  itself across the water instead. Press ▶." `set` {mode fish}, `play: true`, `code: {id: 'sail', lines: [6, 6]}` · 7. "Your turn" —
  "Predict how the no-go zone changes with a better sail. Then raise L/D." `controls: ['LD', 'VU']`.
- **Equations:** `app` "Apparent wind" ref '§14.9 (no number)' $\mathbf U_r=\mathbf U-\mathbf V$, $\tan\beta_a=\dfrac{U\sin\theta}{U\cos\theta+V}$ (ours), live · `drive`
  "Drive and side force" ref 'ours (D22)' $F_{\text{drive}}=L\sin\beta_a-D\cos\beta_a,\ F_{\text{side}}=L\cos\beta_a+D\sin\beta_a$, live · `lim` "Closest course" ref
  'ours (D22)' $\beta_a=\varepsilon_A+\varepsilon_H,\ \varepsilon=\tan^{-1}(D/L)$ · `fish` "Fish tail" ref '§14.8 (no number)' thrust $=L\sin\theta$ at rate
  $UL\sin\theta$; side force $L\cos\theta$ at rate $VL\cos\theta$; $V/U=\tan\theta$ · `vmg` "Velocity made good" ref 'ours' $V\cos\theta$.
- **Check yourself:** (1) "Is the apparent wind closer to the bow or further from it than the true wind (boat moving)?" —
  "Closer: the boat's own speed adds a headwind." · (2) "With L/D = 6, at what apparent-wind angle does the drive vanish?" —
  "tan⁻¹(1/6) = 9.5°." `set {thetaDeg: 14, VU: 0.5}` · (3) "In fish mode with no drag, what is the efficiency, whatever V/U?" —
  "1: $VL\cos\theta=UL\sin\theta$ because V = U tan θ." `set {sys: 'fish', LD: 20}` · (4) "Why does a boat without a keel not sail
  upwind?" — "The side force (several times the drive) would simply push it sideways."
- **Selftest parity rows:** `{name: 'Ur', js: apparent(6, 3, rad(45)).Ur, py: 'ch14.apparent_wind(6.0, 3.0, np.deg2rad(45.0))["Ur"]',
  rtol: 1e-12}` (8.3938) · `{name: 'beta_a', …["beta_a"], rtol: 1e-12}` (0.52984 rad) · `{name: 'drive', js: sail(6, 3, rad(45), 1.2, 0.2,
  1.225, 20).drive, py: 'ch14.sail_forces(6.0, 3.0, np.deg2rad(45.0), CL=1.2, CD=0.2, rho=1.225, A=20.0)["drive"]', rtol: 1e-10}`
  (374.6 N) · `{name: 'lateral', …["lateral"], rtol: 1e-10}` (980.9 N) · `{name: 'limit', js: limit(6, 8), py:
  'ch14.sailboat_course_limit(6.0, 8.0)', rtol: 1e-12}` (0.28950 rad) · `{name: 'fish eta', js: foil(1, 0.75, 10, 1).efficiency, py:
  'ch14.oscillating_foil(1.0, 0.75, 10.0, 1.0)["efficiency"]', rtol: 1e-12}` (0.80620) · `{name: 'vmg', js: vmg(3, rad(45)), py: 'ch14.velocity_made_good(3.0, np.deg2rad(45.0))', rtol: 1e-12}` · `{name: 'steady speed', js: equil(6, rad(45), 1.2, 0.2, 1.225, 20, 1000, 0.05).V, py: 'ch14.sailboat_equilibrium(6.0, np.deg2rad(45.0), CL=1.2, CD=0.2, rho_air=1.225, A_sail=20.0, rho_water=1000.0, CDA_hull=0.05)["V"]', rtol: 1e-8}` (3.9508 m/s) · invariants `{name: 'fish ideal', js: foil(1, 0.4, 10, 0).efficiency,
  expect: 1, rtol: 1e-12}` and `{name: 'drive zero at the limit', js: driveAtBeta(Math.atan(1/6), 6), expect: 0, atol: 1e-12}`.
- **Fit plan:** 360×640: status (2 lines max) · mode chips · `plan` (55 %) over `curve` (45 %); `power` hidden — "drive/L 0.36 ·
  limit 16.6°" in the `plan` title; `LDh` optional. 844×345: `plan` beside `curve`. Desktop: rows [1.25, 1].

### B1 · wake_rollup
- **Title:** "Does a wing really shed two tip vortices?" · **Summary:** "It sheds a whole sheet; watch the sheet wind itself into
  two vortices while its circulation and impulse stay fixed." · **CORE:** C08 (also N36, N37, N38, N67) · **Reference:** ch05
  `vortex_sheet_rollup` and `angular_frequency_explorer_1.html`. Built only if an explainer above fails review; labelled
  qualitative.
- **meta:** `viz:order 10` · `viz:sections 14.6` · `viz:equations none` · `viz:fluidpy ch14.wake_rollup ch14.rolled_up_pair
  ch14.trailing_sheet_strength ch14.tip_vortex_estimate` · `viz:derivations D11`.
- **Physics (JS ↔ Python):** `sheet(Gn, s, n)` → point vortices of strength $-\dfrac{d\Gamma}{dy}\,dy$ at cosine-spaced stations ↔
  `ch14.trailing_sheet_strength`; `stepRK4(state, dt, eps)` (desingularised point-vortex velocities) ↔ the right-hand side used
  by `ch14.wake_rollup` (parity on the invariants and on the first step only); `pair(Gn, s)` ↔ `ch14.rolled_up_pair`.
- **Views:** 1. `rear` "Looking from behind" (`equal`): the sheet as dots coloured by sign (amber port, blue starboard), the two
  centroids (crosses), the descent of the pair. 2. `inv` "What cannot change" (`hidePortrait`): total circulation of each half,
  the spacing of the centroids, the vertical impulse — three flat amber lines against time. **Controls:** loading chips
  (elliptic · rectangular wing · with flap) · number of vortices (40…200) · core size · transport `t`. **Presets:** the three
  loadings; "small core". **Status:** "t = 2.4 s/Γ₀: centroids 0.785 s apart (π/4 s) — qualitative past the first turn".
- **Depth features:** Explain + Code + Derivation · linked views · transport · presets · status · terms (the three invariants).
- **Explain:** 0. *What the views show* — "each dot is a piece of the trailing sheet seen end-on; amber turns one way, blue
  the other." 1. *What is shed* — "$d\Gamma_{\text{trail}}=-\dfrac{d\Gamma}{dy}\,dy$: strongest near the tips." 2. *What it rolls into* — "each half
  carries $\Gamma_{\text{tip}}=\Gamma_0$; the centroids are $b'=\dfrac1{\Gamma_0}\displaystyle\int\Gamma\,dy=\tfrac\pi4s$ apart for the ellipse and sink at $\Gamma_0/2\pi b'$ (ours)."
  3. *Now* — live circulation, spacing, impulse. 4. *Reading the current setting* — early: "the edges curl first, where the
  sheet is strongest"; late: "two vortices; their spacing is what it was at the start — only the shape changed. The details of
  the spiral depend on the core size and are not a prediction."
- **Derivation tab:** **D11** (4 steps), `view: 'rear'`; step 4 `watch` "amber and blue: the two half-sheets turn in opposite
  senses". **Code:** `pair = ch14.rolled_up_pair(Gn, {{s}})  # Gamma0 {{G0}}, spacing {{bp}}, descent {{wd}}` and the roll-up call.
- **Walkthrough (5 steps):** 1. "Two vortices?" — "Behind a wing you see two vortices. But the wing sheds a sheet." · 2. "The
  sheet" — "Strongest at the tips, opposite sign on the two halves." `eq: 'shed'`, `derive: {id: 'D11', step: 3}` · 3. "Roll-up" —
  "Press ▶: each half winds round its own edge." `play: true` · 4. "What stays" — "Circulation and centroid spacing do not
  change." `terms: true` · 5. "Your turn" — "Predict whether a flatter loading gives vortices closer together or further apart."
- **Equations:** `shed` $d\Gamma_{\text{trail}}=-\dfrac{d\Gamma}{dy}\,dy$ · `pair` $\Gamma_{\text{tip}}=\Gamma_0$, $b'=\dfrac1{\Gamma_0}\displaystyle\int\Gamma\,dy$, descent $\Gamma_0/2\pi b'$ (ours).
- **Check yourself:** (1) "What is the centroid spacing for the elliptic loading?" — "π/4 of the span." · (2) "Does the total
  circulation of the wake change during roll-up?" — "No (and it is zero: the halves cancel)." · (3) "Which loading gives the
  widest pair?" — "The flattest (the rectangle's)."
- **Selftest parity rows:** `{name: 'spacing', js: pair([1], 1).spacing, py: 'ch14.rolled_up_pair(np.array([1.0]), 1.0)["spacing"]',
  rtol: 1e-10}` (0.7854) · `{name: 'descent', …["descent_speed"], rtol: 1e-10}` · `{name: 'sheet strength', js: sheetStrength([1], 1.0,
  1), py: 'ch14.trailing_sheet_strength(np.array([1.0]), 1.0, 1.0)', rtol: 1e-10}` · invariant `{name: 'total circulation', js:
  sumGamma(state), expect: 0, atol: 1e-12}`.
- **Fit plan:** 360×640: `rear` alone; the invariants as three readouts. Desktop: two views side by side.

---

## Part D — runtime budget (full run < 5 min on a laptop / Colab CPU; target 90–120 s)

Per-cell costs are estimates: measured where the designer ran the code (✓: the four sympy checks about 1 s each; the
`airfoils` and `lifting_line` calls of the *expect* numbers all under 0.1 s; a 400-step lumped-vortex march in pure Python
0.3 s) and from analysis §4 otherwise. The builder reports the measured numbers.

| Section | What costs time (per cell) | Full | FAST |
|---|---|---|---|
| front matter, §14.1 | imports (sympy, plotly, scipy); `illustrative_inputs`; three tables | 5 s | 5 s |
| §14.2 C01 | three planforms; section geometry; **F1** (19 + 18 steps × 1 trace × 200 pts) 2 s | 5 s | 4 s |
| §14.3 C02, C03 | flat-plate streamlines on 120 × 80 (closed form, inverse map) ×3; **A1** (40 frames, dpi 80; FAST 20) 4 s / 2 s; bisection; march 400 steps ×2 (0.6 s ✓); **A2** (80 frames; FAST 40) 5 s / 3 s; E3 iframe | 14 s | 9 s |
| §14.4 C04 | four-circle figure; **D03** and **D04** sympy checks (1.1 s + 1.0 s ✓); `circular_arc_sym` (< 1 s); shape table (4 × `zhukhovsky_shape` at n = 20001 ≈ 0.2 s); **A3** (60 frames; FAST 30) 3 s / 2 s; **F2** (26 + 25 steps × 2 traces × 400 pts) 3 s; E1 iframe | 13 s | 10 s |
| §14.5 C05, C06, C07 | streamlines in two planes (2 × 160 × 100); contour integrals; **F3** (21 steps × 1 trace × 200 pts) 1 s; **F4** (17 steps × 2 traces) 1 s; sympy `limit` (1 s); `zhukhovsky_surface` n = 2000 ×20 (< 1 s); **F5** (17 steps × 2 traces × 400 pts) 2 s; two force routes; Thwaites sweep over 12 angles ≈ 3 s (FAST 6 angles); live cell (not executed on the page); E2, E5, E4 iframes | 16 s | 11 s |
| §14.6 C08–C14 | staircase and sign-map figures; **A4** roll-up n = 200, ≈ 8 s uncached (FAST n = 80 ≈ 2 s; cached in `outputs/ch14/cache/`); `quad` principal values (≈ 40 calls, 0.5 s); **F6** (7 steps × 4 traces) 1 s; **D15** sympy check (0.4 s ✓); `solve_lifting_line` N ≤ 64 ×60 (< 0.5 s); **F7** (19 steps × 4 traces × 200 pts) 2 s; **D18** sympy check (1.1 s ✓); **F8a** (28 steps × 1 trace) 1 s; live cell; E6, E7 iframes | 22 s | 12 s |
| §14.7 C15 | generic curves; share-vs-speed figure; **F8b** (21 steps × 2 figures × 3 traces) 2 s; E8 iframe | 5 s | 4 s |
| §14.8–§14.9 C16 | triangles; sail polygon; tacking path; **A5** (80 frames; FAST 40) 5 s / 3 s; **F9** (21 steps × 1 trace × 180 pts) 1 s; E9 iframe | 9 s | 6 s |
| explainers | 9 `show_viz` cells (iframes; no computation) | 3 s | 3 s |
| **Total** | | **≈ 92 s (first run, roll-up uncached); ≈ 85 s afterwards** | **≈ 64 s** |

**FAST plan.** `FAST = setup_notebook()`; animations ≤ 80 frames at dpi 80 (FAST 40 or fewer): A1 and A4 `player="frames"`
(the reader should stop at Γ/Γ_K = 0, ½, 1 and at four roll-up times); A2, A3, A5 `player="video"`. Slider figures ≤ 28 steps ×
≤ 4 traces × ≤ 400 points (≈ 100–200 kB each; fourteen figures in the nine groups F1–F9). Grids: streamlines 120 × 80 (C02),
160 × 100 (C05); surface n = 2000 for numbers, n = 400 for plots; lifting line N = 32 for numbers, N = 12 for bars; roll-up
n = 200 (FAST 80). The four ★★★ `check_src` cells always run live (≈ 4 s together ✓). **Cached arrays:** only the roll-up run
(`outputs/ch14/cache/wake_rollup_n200.npz`, git-ignored, regenerated on demand, fallback: run with n = 80); **no data file is
committed** for this chapter. Page budget: five animations < 1.5 MB each, fourteen plotly figures ≈ 150 kB each: ≈ 9 MB < 15 MB.

---

## Part E — prerequisite ledger
Every concept, symbol, maths tool and Python function or idiom the notebook or its explainers use, with where it is explained.
"primer (in Cxx)" = a 📎 primer placed in that block before first use (the Concept text is the primer term, used verbatim in
`nb.primer`); "knowledge/primers.md: <term> (chNN Pnn) — reminder" = a one-line reminder naming the earlier primer; a CORE/RECAP
id alone = taught there; "Cxx (Nnn)" = the NOTE placed in that block; "Cxx (D0n)" = the derivation where it is built; "Cxx (gloss
…)" = one sentence where it is used; "Cxx (🔁 Ch. n)" = an id-less recap cell in that block. Items of §14.1 are first used in the
story that opens C01 and are entered as C01; items of §14.8 open C16 and are entered as C16. New primers P336–P350 (15) in
first-use order: P336, P337 (C01) · P339 (C03) · P340, P341, P342 (C04) · P338, P343 (C06) · P344 (C09) · P350 (C10) · P346, P347 (C11) ·
P345 (C12) · P348 (C13) · P349 (C15).

| Concept | First used in | Explained by |
|---|---|---|
| the conventions block: axes of §14.3–§14.5 and of §14.6; Γ clockwise; $dw/dz=u-iv$; the overloaded letters b, a, β, c, e, Λ, γ, θ, w, A, Γ; K per radian; no lapse rate in this chapter | C01 | front matter (A.0 row 5, before C01), repeated per block |
| slips kept apart from traps; `ch14.book_slips()`, `ch14.traps()`, `ch14.conventions_table()`, `ch14.illustrative_inputs()` | C01 | front matter (A.0 row 5) |
| aerodynamics; the three standing assumptions; stall as the exception | C01 | C01 (N01) |
| Mach number, speed of sound, incompressible when Ma < 0.3; `ch14.compressibility_check` | C01 | C01 (N01; 🔁 Ch. 1, Ch. 4) |
| ideal (inviscid, irrotational) flow | C01 | C01 (N01; 🔁 Ch. 6) |
| thin boundary layers | C01 | C01 (N01; 🔁 Ch. 9) |
| fuselage, wing root, wing tip, span s, chord c(y), planform area A, taper, sweepback; `LL.planform` | C01 | C01 (N02) |
| aspect ratio Λ and mean chord $\bar c$; `LL.aspect_ratio`, `LL.mean_chord` | C01 | C01 (D01, N03) |
| steady level flight: lift = weight, thrust = drag; `ch14.level_flight` | C01 | C01 (D01, N05) |
| stall speed; `ch14.stall_speed` | C01 | C01 (D01 step 5); C15 (N73) |
| force per unit span, section coefficient and wing coefficient | C01 | primer (in C01) |
| dynamic pressure ½ρU² | C01 | C01 (📎 P336, the block on force per unit span) |
| lift and drag coefficients; `SIM.lift_coefficient`, `SIM.drag_coefficient`, `AF.section_coefficient` | C01 | R01 |
| which area A (planform, frontal, chord per unit span) | C01 | R01 (table) |
| friction drag and pressure (form) drag | C01 | R02 |
| substitution; square root of both sides | C01 | C01 (D01, gloss) |
| trapezoid rule; `np.trapezoid` | C01 | knowledge/primers.md: trapezoid rule (ch01 P37) — reminder |
| `assert np.allclose` / `np.isclose` | C01 | knowledge/primers.md: assert np.allclose (ch01 P15) — reminder |
| `slider_figure`, `animate` / `show_animation`, `show_viz`, `live` | C01 | knowledge/primers.md: slider_figure, animate, show_viz (ch01 P17, P16, P18) — reminder |
| airfoil section: leading and trailing edge, chord line, camber line, camber, thickness, angle of attack α | C01 | C01 (N07) |
| splitting a closed contour into upper and lower surfaces | C01 | primer (in C01) |
| `np.argmin` / `np.argmax`, slicing with `[::-1]` | C01 | C01 (📎 P337, the block on splitting a contour) |
| `np.interp` | C01 | knowledge/primers.md: np.interp (ch07 P182) — reminder inside P337 |
| `AF.section_geometry`, `AF.zhukhovsky_airfoil` (first call, for the labelled section) | C01 | C01 (N07); taught fully in C04 |
| lift and drag resolved relative to the stream; `AF.resolve_lift_drag` | C01 | C01 (N08) |
| pitch, roll, yaw axes | C01 | C01 (N04) |
| control surfaces (aileron, elevator, rudder, flap, slot, spoiler, trim tab) | C01 | C01 (N06) |
| projectile range with drag; `ch14.projectile_range_drag` | C01 | C01 (N83) |
| circulation Γ as a loop integral of velocity | C02 | R04 (🔁 Ch. 3) |
| Kutta–Zhukhovsky theorem; `ch06.lift_per_span` | C02 | R04 |
| clockwise and counterclockwise circulation; `Gamma_cw=` / `Gamma_ccw=` | C02 | R04 (sign table) |
| complex potential w = φ + iψ; streamlines as contours of ψ | C02 | C02 (🔁 Ch. 6, in the maths rows) |
| stagnation point; stagnation points of a cylinder with circulation | C02 | C02 (maths row ii; 🔁 Ch. 6) |
| the Kutta condition | C02 | C02 |
| circulation as the net vorticity of the two boundary layers | C02 | C02 (N10) |
| `AF.flat_plate_flow`, `AF.rear_stagnation_point`, `AF.kutta_circulation_root`, `CONF.circle_flow_zeta` | C02 | C02 (code rows) |
| root finding by bisection; `brentq` | C02 | knowledge/primers.md: brentq (ch03 P108) — reminder |
| logarithmic axis | C02 | C02 (gloss in the figure's 'how to read it' note: each step up the axis is a factor of ten) |
| three ways to control lift (sharp edge, angle, camber) | C02 | C02 (N11) |
| Kelvin's circulation theorem; material circuit | C03 | C03 (🔁 Ch. 5) |
| additivity of circulation over sub-circuits (Stokes) | C03 | C03 (🔁 Ch. 3, Ch. 2) |
| the three stages of an impulsive start | C03 | C03 (N14) |
| why the stagnation point moves (boundary layer at the edge) | C03 | C03 (N15) |
| starting vortex; stopping vortex | C03 | C03 (D02) |
| point vortex and its induced velocity | C03 | C03 (🔁 Ch. 5, inside P339) |
| discrete-vortex time stepping | C03 | primer (in C03) |
| `ch14.impulsive_start_lumped_vortex`; (optional) `ch14.wagner_function` | C03 | C03 (code row; C.6 item 4) |
| viscosity as the source of lift | C03 | C03 (N16) |
| conformal map; the Zhukhovsky transformation; `CONF.joukowski`, `CONF.joukowski_inverse` | C04 | R05 |
| complex numbers in numpy (`1j`, `np.exp(1j*t)`, `.real`, `.imag`, `np.abs`) | C04 | knowledge/primers.md: the complex plane in numpy (ch06 P153) — reminder |
| Euler's formula | C04 | knowledge/primers.md: Euler's formula (ch01 P45) — reminder |
| the centred circle maps to a slit | C04 | R06 |
| thin-airfoil theory, Theodorsen's method (named) | C04 | C04 (N18) |
| the circle parametrisation μ, a from (b, e, β); `AF.zhukhovsky_circle` | C04 | C04 (N25) |
| polar form of a point of the circle, R and θ from the origin | C04 | C04 (N19, D03 step 1) |
| difference of squares | C04 | C04 (D03 step 4, gloss) |
| the cosine rule | C04 | primer (in C04) |
| cot, csc and the double-angle identities | C04 | primer (in C04) |
| completing the square into a circle | C04 | knowledge/primers.md: completing the square into a circle (ch11 P272) — reminder |
| small-angle approximation | C04 | knowledge/primers.md: small-angle approximation (ch03 P100) — reminder |
| circular arc; `AF.circular_arc` | C04 | C04 (N24, D03) |
| sympy (`symbols`, `simplify`, `series`, `limit`, `integrate`) | C04 | knowledge/primers.md: sympy (ch01 P40, ch04 P117) — reminder |
| binomial approximation with a complex small quantity | C04 | primer (in C04) |
| orders of smallness ("first order in e") | C04 | knowledge/primers.md: orders of smallness (ch02 P68) — reminder inside P342 |
| maximum of a function by a zero derivative | C04 | C04 (D04 step 8, gloss) |
| shape rules: chord, camber, thickness; `AF.zhukhovsky_shape(exact=)` | C04 | C04 (N26, D03, D04) |
| critical point of a map; angle doubling | C04 | C04 (D05; 🔁 Ch. 6) |
| Taylor expansion | C04 | knowledge/primers.md: Taylor expansion (ch01 P26) — reminder |
| cusp against finite-angle trailing edge | C04 | C04 (N27, D05) |
| corner flows; `ch14.trailing_edge_wedge`, `POT.power_branch` | C04 | R07 |
| parametric curves | C04 | knowledge/primers.md: parametric curves (ch03 P92) — reminder |
| ζ-plane flow past an off-centre circle at incidence | C05 | C05 (N28, D06) |
| complex logarithm and its branch | C05 | knowledge/primers.md: complex logarithm and branch cuts (ch06 P155) — reminder |
| complex velocity $dw/dz=u-iv$ | C05 | R08 (trap T2 callout) |
| chain rule | C05 | R08; knowledge/primers.md: chain rule (ch01 P49) — reminder |
| velocity in the two planes (14.9) | C05 | C05 (D07) |
| circulation unchanged by the map | C05 | C05 (D08) |
| closed-loop integral of a multi-valued function | C05 | knowledge/primers.md: closed-loop integral of an exact differential (ch05 P135) — reminder |
| periodic trapezoid rule on a closed loop | C05 | knowledge/primers.md: periodic trapezoid rule (ch05 P136) — reminder |
| the outside root of the inverse map; `AF.zhukhovsky_inverse` | C05 | C05 (code row, trap T8); knowledge/primers.md: complex square roots (ch06 P159) — reminder |
| `AF.zhukhovsky_flow` and its methods (`velocity`, `stream`, `speed`, `dwdz`, `dwdzeta`, `inside`) | C05 | C05 (code row) |
| `AF.stagnation_points` | C05 | C05 (figure row) |
| surface speed on the circle (14.10) | C06 | R09 |
| the three θ's | C06 | R09 (trap T4 callout) |
| radians and degrees for angles and for slopes | C06 | primer (in C06) |
| Kutta circulation (14.11); `AF.kutta_circulation` | C06 | C06 (D09, N30) |
| lift coefficient of a section (14.12); `AF.zhukhovsky_lift_coefficient` | C06 | C06 (D09) |
| flat plate at incidence | C06 | C06 (N34) |
| l'Hôpital's rule for a 0/0 limit | C06 | primer (in C06) |
| finite trailing-edge speed; `AF.trailing_edge_speed`, `ch14.edge_speed` | C06 | C06 (N29) |
| absolute angle of attack; zero-lift angle; `AF.zero_lift_angle` | C06 | C06 (N31) |
| theory against measurement; `ch14.generic_lift_curve` (qualitative) | C06 | C06 (N32) |
| displacement thickness | C06 | C06 (N32; 🔁 Ch. 9, one clause) |
| Kármán–Trefftz map (named) | C06 | C06 (N33) |
| pressure coefficient $C_p$; `POT.pressure_coefficient` | C07 | R03 |
| Bernoulli's equation | C07 | C07 (🔁 Ch. 4) |
| suction peak; adverse pressure gradient | C07 | C07 (plain words; N13) |
| `AF.zhukhovsky_surface`, `ch14.surface_point` | C07 | C07 (code row) |
| pressure force −p n ds; outward normal of a contour | C07 | C07 (D10 steps 1–2) |
| closed line integral round a contour | C07 | C07 (D10 step 3); knowledge/primers.md: closed-loop integral (ch05 P135) — reminder |
| normal-force coefficient as the loop area; `AF.normal_force_coefficient`, `AF.lift_from_cp` | C07 | C07 (D10, N09) |
| `AF.surface_circulation` | C07 | C07 (code row; N10's number) |
| Blasius force theorem; `POT.blasius_force`, `AF.zhukhovsky_force` | C07 | R10 |
| leading-edge suction; `AF.flat_plate_forces`, `ch14.edge_flow_suction` | C07 | C07 (N12) |
| stall; separation where the wall shear stress vanishes | C07 | C07 (N13) |
| Thwaites' method; `ch14.suction_side_separation` (qualitative) | C07 | C07 (N13; 🔁 Ch. 9) |
| tip flow; trailing vortex sheet | C08 | C08 (N36) |
| vortex sheet as a jump in tangential velocity | C08 | C08 (🔁 Ch. 5) |
| Helmholtz's theorem; the closed vortex loop of a wing | C08 | R11 |
| the three assumptions of lifting-line theory | C08 | C08 (N39) |
| bound vortex; lifting line; `LL.circulation_series` (first call) | C08 | C08 (N40) |
| derivative as a difference; first-order Taylor step | C08 | C08 (D11 step 2); knowledge/primers.md: Taylor expansion (ch01 P26) — reminder |
| shed filament strength; `LL.trailing_sheet_strength` | C08 | C08 (D11) |
| tip vortices; roll-up; `LL.rolled_up_pair`, `ch14.tip_vortex_estimate`, `ch14.wake_rollup` (qualitative) | C08 | C08 (N37) |
| point-vortex dynamics and invariants | C08 | C08 (🔁 Ch. 5; N37); knowledge/primers.md: systems of ODEs for interacting vortices (ch05 P147) — reminder |
| Lanchester and Prandtl | C08 | C08 (N69); C03 (N17) |
| Biot–Savart law for a filament; straight-segment formula; `BS.segment_speed` | C09 | C09 (🔁 Ch. 5); R12 |
| semi-infinite filament = half an infinite line | C09 | R12 |
| improper integral as a limit | C09 | knowledge/primers.md: improper integral as a limit (ch05 P144) — reminder |
| Cauchy principal value | C09 | primer (in C09) |
| `scipy.integrate.quad` (with `weight="cauchy"`, `points=`) | C09 | knowledge/primers.md: scipy.integrate.quad (ch03 P87) — reminder inside P344 |
| downwash (14.13); `LL.downwash` | C09 | C09 (D12) |
| horseshoe vortex; `LL.horseshoe_downwash`, `LL.horseshoe_velocity`, `LL.horseshoe_wing` | C09 | C09 (N60, worked example) |
| staircase of horseshoes; `LL.staircase_from_loading`, `LL.staircase_downwash` | C09 | C09 (from-scratch row; F6) |
| downwash ahead of and behind the wing; `LL.downwash_along_x` | C09 | C09 (N41) |
| upwash outboard of the tips; formation flight | C09 | C09 (N65) |
| method of images; ground effect; `LL.ground_effect_factor` | C09 | C09 (N66; 🔁 Ch. 5, Ch. 6) |
| parabolic loading | C09 | C09 (N61) |
| relative velocity as a vector subtraction | C10 | primer (in C10) |
| `np.arctan2`, `np.hypot` | C10 | knowledge/primers.md: np.arctan2 (ch02 P70) — reminder inside P350 |
| induced angle ε; `LL.induced_angle` | C10 | C10 (D13 step 1, N42) |
| effective angle of attack (14.14); `LL.effective_angle` | C10 | C10 (D13 step 2, N43) |
| strip theory: each strip a two-dimensional section; `LL.strip_forces` | C10 | C10 (D13 steps 3–6, N44) |
| lift and induced drag of a wing (14.15); `LL.lift_and_induced_drag` | C10 | C10 (D13) |
| induced drag as wake kinetic energy; `LL.wake_kinetic_energy_per_length` | C10 | C10 (N38) |
| power of a force | C10 | knowledge/primers.md: power of a force (ch04 P127) — reminder in N38 |
| Rankine vortex; core size; `ch14.vortex_pair_energy_per_length` | C10 | C10 (N67; 🔁 Ch. 3) |
| span angle γ; `LL.span_angle`, `LL.span_position` | C11 | C11 (N49) |
| half-range sine series | C11 | primer (in C11) |
| integrals of sines over a full period; Galerkin truncation | C11 | knowledge/primers.md: integrals of sines and cosines (ch06 P151), Galerkin truncation (ch11 P275) — reminder inside P346 |
| sine series for the circulation (14.18) | C11 | C11 (N50) |
| even and odd functions | C11 | knowledge/primers.md: even and odd functions (ch11 P261) — reminder in N50 |
| substitution in an integral (change of variable) | C11 | knowledge/primers.md: substitution in an integral (ch03 P106) — reminder |
| term-by-term differentiation of a series | C11 | C11 (D14 step 4, gloss) |
| the downwash integral in γ; slip #1 | C11 | C11 (D14) |
| product-to-sum identity and a three-term recursion | C11 | primer (in C11) |
| proof by induction; linear recurrences | C11 | knowledge/primers.md: proof by induction (ch12 P306), linear recurrence (ch10 P236) — reminder inside P347 |
| Glauert's integral; `LL.glauert_integral` | C11 | C11 (D15, N51) |
| Chebyshev polynomials (used only inside the sympy check of D15) | C11 | C11 (D15 check, one comment line: "cos kγ is a polynomial in cos γ") |
| downwash of the series (14.19); `LL.downwash_series(printed_middle=)` | C11 | C11 (D15 step 10) |
| section lift law on a finite wing (14.16); lift-curve slope K; `LL.section_lift_coefficient` | C12 | C12 (D16 step 1, N46) |
| circulation from the lift coefficient; `AF.circulation_from_cl` | C12 | C12 (D16 step 2, N47) |
| Prandtl's lifting-line equation (14.17) | C12 | C12 (D16 step 3, N48) |
| integral equation | C12 | primer (in C12) |
| `np.linalg.solve`; `np.eye`, `@` | C12 | knowledge/primers.md: np.linalg.solve (ch01 P57) — reminder inside P345 |
| collocation; condition number; `np.linalg.cond` | C12 | knowledge/primers.md: collocation and the condition number (ch06 P162) — reminder |
| the lifting-line equation in the coefficients (14.20); `LL.lifting_line_matrix`, `LL.lifting_line_rhs`, `LL.lifting_line_modes` | C12 | C12 (D16, N52) |
| `LL.solve_lifting_line`, `LL.lifting_line_residual` | C12 | C12 (code rows) |
| twist (washout); elliptic twist of a rectangle; `LL.elliptic_twist_rectangular` | C12 | C12 (N63) |
| rolling moment in a wake vortex; `ch14.wake_vortex_rolling_moment` | C12 | C12 (N64) |
| limits of lifting-line theory | C12 | C12 (N68) |
| `lambda` | C12 | knowledge/primers.md: functions as arguments and lambda (ch01 P29) — reminder in the code comment |
| a double sum collapsed by orthogonality | C13 | primer (in C13) |
| Kronecker delta $\delta_{nm}$ | C13 | C13 (📎 P348, the block on the double sum) |
| `np.outer` | C13 | C13 (📎 P348, the block on the double sum) |
| lift from the first coefficient (14.21); `LL.lift_from_series` | C13 | C13 (D17, N53) |
| induced drag from the coefficients (14.22); `LL.induced_drag_from_series` | C13 | C13 (D18) |
| span efficiency $e_{\text{span}}$; `LL.span_efficiency` | C13 | C13 (N54) |
| elliptic loading (14.23); `LL.elliptic_circulation` | C14 | C14 (D19 steps 3–4, N55) |
| uniform downwash (14.24); `LL.elliptic_downwash` | C14 | C14 (D19 step 5, N56) |
| induced drag of the elliptic loading (14.25); `LL.elliptic_induced_drag` | C14 | C14 (D19 steps 6–7, N57) |
| induced drag coefficient (14.26); `LL.induced_drag_coefficient` | C14 | C14 (D19 steps 8–9) |
| elliptic planform against elliptic loading (14.27) | C14 | C14 (N58, D20 steps 1–3) |
| area of an ellipse | C14 | C14 (D20 step 6, gloss) |
| finite-wing lift slope; `LL.finite_wing_lift_slope`, `LL.elliptic_wing` | C14 | C14 (D20, N59) |
| sea-level density from the standard atmosphere | C14 | C14 (N62; 🔁 Ch. 1) |
| drag breakdown: friction, form, induced, profile | C15 | C15 |
| measured lift and drag curves; `ch14.generic_lift_curve`, `ch14.generic_drag_curve` (qualitative) | C15 | C15 (N70) |
| a stalled section | C15 | C15 (N71) |
| chord Reynolds number; `ch14.chord_reynolds`, `ch14.stall_regime` (qualitative) | C15 | C15 (N72) |
| transition before separation; the drag crisis | C15 | C15 (N72; 🔁 Ch. 9) |
| flaps and landing speed | C15 | C15 (N73) |
| reading a drag polar; the minimum of a/x + bx | C15 | primer (in C15) |
| quotient rule; `minimize_scalar` | C15 | knowledge/primers.md: quotient rule (ch11 P265), minimize_scalar (ch07 P170) — reminder inside P349 |
| parabolic polar; `ch14.drag_polar`, `ch14.drag_breakdown`, `ch14.drag_share` | C15 | C15 (D21) |
| best glide; `ch14.best_glide`, `ch14.speed_min_drag`, `ch14.drag_vs_speed` | C15 | C15 (D21, N74) |
| "induced drag is the largest part" as a computed share | C15 | C15 (N75) |
| shock waves, sweepback (named) | C15 | C15 (N76) |
| oscillating foil: thrust, side force, the two powers, efficiency; `ch14.oscillating_foil` | C16 | C16 (N77, in §14.8) |
| bird flight: downstroke and upstroke | C16 | C16 (N79, in §14.8) |
| tail flukes, tubercles, flexible fins; low-Reynolds-number flight (named) | C16 | C16 (N78, N80, in §14.8) |
| true wind, boat velocity, apparent wind; `ch14.apparent_wind` | C16 | C16 (D22 steps 1–3); knowledge of P350 (primer in C10) |
| drive and side force of a sail; `ch14.sail_forces`, `ch14.sail_drive_curve` | C16 | C16 (D22 steps 4–6) |
| drag angle ε = tan⁻¹(D/L); tangent of a difference of angles | C16 | C16 (D22 steps 7–8, gloss) |
| the keel as a hydrofoil; closest course; `ch14.sailboat_course_limit`, `ch14.sailboat_equilibrium` | C16 | C16 (N81, D22 step 8) |
| tacking; velocity made good; `ch14.velocity_made_good` | C16 | C16 (N82) |
| dot product of two vectors (`@`) | C16 | C16 (gloss in the from-scratch cell: `a @ b` multiplies matching components and adds them) |
| f-strings | C01 | knowledge/primers.md: f-strings (ch01 P04) — reminder |
| `pd.DataFrame` for small tables | C01 | front matter (A.0 row 5; gloss: one comment line where the first table is built) |
| explainers: `Viz.three` 3-D scene with a 2-D fallback (E6) | C10 | E6 Explain section 0 and its fit plan (one hint line) |
| explainers: local `cx` complex helpers and `luSolve` (E1, E2, E4, E6, E7) | C04 | Part B preamble (builder note; nothing the reader sees) |

---

## Part F — derivation storyboards

One block per `D` row of the curation (§4b), in order. Builders copy each block word for word into `nb.derivation(id, title,
ref=…, goal=…, start=(tex, plain), plan=[…], uses=[…], steps=[dict(did=…, tex=…, why=…, plain=…), …], result=(tex, plain),
interpret=…, check=…, check_src=…)` and, where the heading names an explainer, into that explainer's `derivations: [...]` (same
`did` titles, same step count). Every step has four parts: **did** (the move) · **tex** (the new line) · **why** (why it is
allowed and why we make it) · **plain** (what the line says). Wherever a step produces a numbered book equation, the number
stands directly after that step's **tex** — beside the equation it names. Colours of terms as header convention 9. The four
★★★ blocks (D03, D04, D15, D18) carry a **sympy check** with a `check_src` listing that re-runs the derivation's own
construction; **each listing was executed as written for this design and passed** (the printed line is quoted under it).
Γ is clockwise throughout (the book's sign).

### D01 · Aspect ratio two ways (14.1)–(14.2), and what level flight asks of the wing — ★, 5 steps, in C01 (notebook)
- **Goal:** reduce a wing to one shape number, and find how much lift coefficient it must deliver to carry the aircraft — and
  the slowest speed at which it can. · **Start:** $\Lambda\equiv s^2/A$ (14.1), the definition — "span squared over planform area". ·
  **Plan:** • name the mean chord • rewrite Λ with it • balance lift and weight • solve for the coefficient, then for the
  speed. · **Tools:** substitution; the lift coefficient $C_L\equiv\dfrac{L}{\tfrac12\rho U^2A}$ (4.108) (R01); a square root. · **Assumptions:** steady,
  level, unaccelerated flight (step 3); A is the planform area (step 4).
- **Steps:**
  1. **did** Define the mean chord · **tex** $\bar c=A/s$ (14.2) · **why** A definition: the planform area divided by the span is
     the chord of the rectangle with the same span and area. We need one chord to stand for a chord that varies. · **plain**
     The average chord is area over span.
  2. **did** Substitute it into the aspect ratio · **tex** $\Lambda\equiv\dfrac{s^2}{A}=\dfrac{s^2}{\bar c\,s}=\dfrac{s}{\bar c}$ (14.1) · **why** Replace A by $\bar c\,s$ from
     step 1 and cancel one factor s; this shows the two printed forms are the same number. · **plain** Aspect ratio is span
     divided by mean chord: how many chords fit along the wing.
  3. **did** Balance the vertical forces · **tex** $L=W$ · **why** In steady level flight the aircraft does not accelerate
     vertically, so by Newton's second law the lift equals the weight (and likewise thrust equals drag). · **plain** The wing
     must carry exactly the aircraft's weight.
  4. **did** Insert the lift coefficient and solve · **tex** $C_L=\dfrac{W}{\tfrac12\rho U^2A}$ · **why** By its definition $L=C_L\,\tfrac12\rho U^2A$ with
     A the planform area; set this equal to W and divide. We want the one number the wing has to produce. · **plain** The
     required lift coefficient is weight over dynamic pressure times wing area; it falls as one over speed squared.
  5. **did** Set the coefficient to its largest value · **tex** $U_{\min}=\sqrt{\dfrac{2W}{\rho A\,C_{L,\max}}}$ · **why** A wing cannot exceed some
     $C_{L,\max}$ (stall). Put $C_L=C_{L,\max}$ in step 4 and solve for U; U is positive, so we take the positive root. · **plain**
     Below this speed the wing cannot hold the aircraft up.
- **Result:** $\Lambda\equiv s^2/A=s/\bar c$ with $\bar c=A/s$; in level flight $C_L=\dfrac{W}{\tfrac12\rho U^2A}$ and $U_{\min}=\sqrt{\dfrac{2W}{\rho A\,C_{L,\max}}}$ — "the shape in one
  number; the task in one number".
- **Check:** units — Λ: m²/m² ✓; $C_L$: N/(Pa m²) ✓; $U_{\min}$: √(N/(kg m⁻³ m²)) = m/s ✓. Limits — a rectangle has Λ = s/c exactly;
  doubling the speed quarters the required $C_L$. Number — s = 10 m, A = 20 m²: Λ = 5; W = 10 kN at 40 m/s, ρ = 1.25: $C_L$ = 0.5.
- **What it means:** every later result is written in $C_L$ and Λ; slow flight is high-$C_L$ flight. **Fails when:** the flight
  is not level and steady (in a turn or a pull-up the lift exceeds the weight).
- **Traps:** A is the planform area, not the frontal area (trap T13); "level, steady" is what makes lift equal weight;
  $C_L\propto1/U^2$, not 1/U.

### D02 · The starting vortex carries the opposite circulation — ★, 5 steps, in C03 (notebook · `starting_vortex`)
- **Goal:** show that a foil which starts from rest and ends with circulation Γ must have left a vortex of circulation −Γ
  behind. · **Start:** Kelvin's theorem $\dfrac{D\Gamma}{Dt}=0$ (5.8) — "the circulation round a circuit that moves with the fluid does
  not change". · **Plan:** • draw one big material circuit round everything • apply Kelvin to it • cut it in two • read off
  the two parts. · **Tools:** Kelvin's theorem (recap in C03); circulation adds over sub-circuits (recap in C03). ·
  **Assumptions:** the circuit itself stays in inviscid, constant-density fluid (step 2); it is large enough that neither the
  foil's boundary layer nor the wake crosses it (step 1).
- **Steps:**
  1. **did** Draw a large material circuit at rest · **tex** $\Gamma_{ABCD}(0)=\displaystyle\oint_{ABCD}\mathbf u\cdot d\mathbf x=0$ · **why** Before the start the fluid is
     at rest everywhere, so the velocity along any closed curve is zero. We choose ABCD so large that it encloses the foil
     and all it will shed. · **plain** At the start nothing circulates.
  2. **did** Apply Kelvin's theorem to it · **tex** $\Gamma_{ABCD}(t)=\Gamma_{ABCD}(0)=0$ · **why** The circuit moves with the fluid and
     stays in fluid where viscosity does not act; what happens *inside* it (boundary layers, the wake) does not matter to
     the theorem. · **plain** The big circuit's circulation stays zero for all time.
  3. **did** Cut the circuit along the line BD · **tex** $\Gamma_{ABCD}=\Gamma_{ABD}+\Gamma_{BCD}$ · **why** The cut BD is walked once in each
     direction by the two smaller circuits, so its contributions cancel and the outer parts add up to the original
     (Stokes' additivity). · **plain** The total is the sum of the two halves.
  4. **did** Name the two halves · **tex** $\Gamma_{ABD}=\Gamma_{\text{foil}},\qquad\Gamma_{BCD}=\Gamma_{\text{starting}}$ · **why** Once the foil has moved on, ABD
     encloses the foil alone and BCD the vorticity it shed at the start. · **plain** One half surrounds the wing, the other
     the vortex left behind.
  5. **did** Combine steps 2–4 · **tex** $\Gamma_{ABCD}=0=\Gamma_{\text{foil}}+\Gamma_{\text{starting}}\ \Rightarrow\ \Gamma_{\text{starting}}=-\Gamma_{\text{foil}}$ · **why** Substitute
     step 4 into step 3 and use step 2. · **plain** The starting vortex has exactly the foil's circulation, with the opposite
     sense.
- **Result:** $\Gamma_{ABCD}=0=\Gamma_{\text{foil}}+\Gamma_{\text{starting}}$ — "the wing does not create circulation; it separates it into two equal and
  opposite parts and flies away from one of them".
- **Check:** units — all terms m²/s ✓. Number — a plate of chord 1 m at 30° reaching 10 m/s: bound +15.7 m²/s (clockwise),
  starting vortex −15.7 m²/s. Code — the lumped-vortex march keeps the sum at zero to round-off at every step.
- **What it means:** every change of lift — speed, incidence, a flap — sheds a vortex; the bound vortex, the two trailing
  vortices of C08 and the starting vortex form one closed loop. **Fails when:** the circuit is crossed by the wake or the
  boundary layer (then viscosity acts *on* the circuit), or the density is not uniform.
- **Traps:** it is the circuit, not its interior, that must be free of viscous action (trap T15); a material circuit drifts
  with the fluid, so it must start large; the book's Γ is clockwise — the starting vortex is counterclockwise (trap T1).

### D03 · A raised circle maps to a circular arc: (14.4) → (14.5) → (14.6) → the arc, and the camber rule of (14.8) — ★★★, 15 steps, in C04 (notebook · `zhukhovsky_airfoil_builder`)
- **Goal:** find the shape that the map makes of a circle whose centre has been raised, and so see where camber comes from.
  The book states two intermediate lines and the result; we fill in every move. · **Start:** a point of the circle written by
  its distance R from the origin and its polar angle θ, ζ = Re^{iθ}, put into $z=\zeta+\dfrac{b^2}{\zeta}$ (14.3) — "the circle of radius
  a = b/cos β centred at Q = ib tan β, which passes through ζ = ±b". · **Plan:** • split the image into x and y • eliminate R •
  find how R depends on θ on this particular circle (cosine rule) • eliminate θ • recognise a circle. · **Tools:** Euler's
  formula (P45, reminder); difference of squares (gloss); the cosine rule (P340); cot, csc and the double-angle identities
  (P341); completing the square into a circle (P272, reminder); small angles (P100, reminder). · **Assumptions:** a > b,
  0 < β < π/2; sin θ ≠ 0 and cos θ ≠ 0 while dividing (steps 3, 11), the excluded points recovered as limits; β ≪ 1 only in
  step 15.
- **Steps:**
  1. **did** Map a point of the circle · **tex** $z=Re^{i\theta}+\dfrac{b^2}{R}e^{-i\theta}$ · **why** Put ζ = Re^{iθ} into the map $z=\zeta+b^2/\zeta$ (14.3); 1/ζ = e^{−iθ}/R. Here θ is the
     polar angle of ζ seen from the *origin* (trap T4). · **plain** The image is the point itself plus its reflection scaled
     by b²/R.
  2. **did** Split into real and imaginary parts · **tex** $x=(R+b^2/R)\cos\theta,\qquad y=(R-b^2/R)\sin\theta$ (14.4) · **why** Euler's formula,
     $e^{\pm i\theta}=\cos\theta\pm i\sin\theta$; collect the cosines (real) and the sines (imaginary). · **plain** The two parts add in x
     and subtract in y.
  3. **did** Isolate the two brackets · **tex** $\dfrac{x}{\cos\theta}=R+\dfrac{b^2}{R},\qquad\dfrac{y}{\sin\theta}=R-\dfrac{b^2}{R}$ · **why** Divide each member of step 2
     by its trigonometric factor (allowed where it is not zero). We want R to appear only in a sum and a difference. ·
     **plain** A sum and a difference of the same two numbers, R and b²/R.
  4. **did** Square and subtract · **tex** $\dfrac{x^2}{\cos^2\theta}-\dfrac{y^2}{\sin^2\theta}=4b^2$ · **why** $(p+q)^2-(p-q)^2=4pq$ with p = R and q = b²/R, whose
     product is b²: R drops out. This is the elimination the book does in one phrase. · **plain** Whatever R is, x and y
     satisfy this one relation with θ.
  5. **did** Clear the denominators · **tex** $x^2\sin^2\theta-y^2\cos^2\theta=4b^2\sin^2\theta\cos^2\theta$ (14.5) · **why** Multiply step 4 by sin²θ cos²θ;
     the result also holds at the points excluded in step 3. · **plain** The book's first relation: no R left.
  6. **did** Write the cosine rule for triangle OQP · **tex** $QP^2=OP^2+OQ^2-2\,(OQ)(OP)\cos(\widehat{QOP})$ · **why** O is the origin, Q the
     centre, P the point on the circle. We need a second relation, special to *this* circle, between R = OP and θ. ·
     **plain** A triangle links the radius of the circle to R and θ.
  7. **did** Insert the three lengths and the angle · **tex** $\dfrac{b^2}{\cos^2\beta}=R^2+b^2\tan^2\beta-2Rb\tan\beta\,\sin\theta$ · **why** QP = a = b/cos β, OQ = b tan β, OP = R. OQ lies on the imaginary axis, OP at θ from the real axis: the
     angle between them is 90° − θ, with cosine sin θ. · **plain** The cosine rule with this circle's numbers.
  8. **did** Simplify with a trigonometric identity · **tex** $R^2-b^2=2Rb\tan\beta\,\sin\theta$ · **why** Move b² tan²β to the left:
     $b^2(1/\cos^2\beta-\tan^2\beta)=b^2$ because 1/cos² − tan² = 1. · **plain** On this circle R² exceeds b² by an amount
     proportional to sin θ.
  9. **did** Divide by R and compare with y · **tex** $2b\tan\beta\,\sin\theta=R-\dfrac{b^2}{R}=\dfrac{y}{\sin\theta}$ (14.6) · **why** R > 0; the middle member is
     exactly the bracket of y in step 3. Now y is tied to θ alone. · **plain** The height of the image point is
     $y=2b\tan\beta\,\sin^2\theta$.
  10. **did** Solve for the squared sine and cosine · **tex** $\sin^2\theta=\dfrac{y}{2b\tan\beta},\qquad\cos^2\theta=\dfrac{2b\tan\beta-y}{2b\tan\beta}$ · **why** Rearrange
      step 9 and use sin² + cos² = 1. We prepare to remove θ from the relation of step 5. · **plain** Both are known once y is known — and
      y lies between 0 and 2b tan β.
  11. **did** Divide the relation of step 5 by sin²θ · **tex** $x^2-y^2\cot^2\theta=4b^2\cos^2\theta$ · **why** A legitimate division where sin θ ≠ 0; it leaves
      θ only inside cot²θ and cos²θ, which step 10 gives in terms of y. · **plain** One more form of the first relation.
  12. **did** Substitute step 10 · **tex** $x^2+y^2+2by\,(\cot\beta-\tan\beta)=4b^2$ · **why** cot²θ = cos²θ/sin²θ = (2b tan β − y)/y, so
      y² cot²θ = y(2b tan β − y); and 4b² cos²θ = 4b² − 2by cot β. Collect terms. · **plain** A relation between x and y
      only: the equation of the image curve.
  13. **did** Use the double-angle identity · **tex** $x^2+y^2+4by\cot2\beta=4b^2$ · **why** cot β − tan β = 2 cot 2β (P341). ·
      **plain** A circle's equation, not yet in standard form.
  14. **did** Complete the square in y · **tex** $x^2+(y+2b\cot2\beta)^2=(2b\csc2\beta)^2$ · **why** Add (2b cot 2β)² to both sides; on the
      right 4b²(1 + cot²2β) = 4b² csc²2β. · **plain** The image lies on a circle centred at (0, −2b cot 2β) with radius
      2b csc 2β.
  15. **did** Read off the height and the camber · **tex** $y_{\max}=2b\tan\beta,\qquad\dfrac{\text{camber}}{c}=\dfrac{2b\tan\beta}{4b}\simeq\dfrac\beta2$ · **why** From step 9 the
      largest y is at θ = 90°; the chord joins z = ±2b, so c = 4b; for β ≪ 1, tan β ≃ β (in radians). · **plain** A raised circle gives an arc of camber ½βc — the camber rule, camber $\simeq\tfrac12\beta c$ (14.8).
- **Result:** $x^2+(y+2b\cot2\beta)^2=(2b\csc2\beta)^2$ with $\cos\beta=b/a$: a circular arc through z = ±2b of height 2b tan β; camber
  $\simeq\tfrac12\beta c$ (14.8) — "raise the circle and the flat plate bends into an arc".
- **Check:** units — every term a length squared ✓. Limits — β → 0: the centre goes to −∞ and the radius to ∞: a straight
  plate; y = 0 at θ = 0, π: the ends are at z = ±2b. Number — b = 1, tan β = ¾: the top ζ = 2i maps to z = 1.5i = 2b tan β ✓, and
  −0.583 + 2.083 = 1.5 ✓. Only 0 ≤ y ≤ 2b tan β is covered, and twice (upper and lower half of the circle). The tangent at
  each end makes the angle 2β with the chord (the inscribed-angle property of a circle through ±2b with that centre).
- **Sympy check** (`check_src`; executed: prints `2*b*t` and passes in 1.1 s):
  ```python
  import sympy as sp                                              # symbolic algebra
  b, R, th = sp.symbols("b R theta", positive=True)               # map constant, distance from O, polar angle
  x = (R + b**2/R)*sp.cos(th)                                     # (14.4): real part of the image
  y = (R - b**2/R)*sp.sin(th)                                     # (14.4): imaginary part of the image
  lhs5 = x**2*sp.sin(th)**2 - y**2*sp.cos(th)**2                  # left side of (14.5)
  assert sp.simplify(lhs5 - 4*b**2*sp.sin(th)**2*sp.cos(th)**2) == 0      # steps 3-5: (14.5) holds for every R
  t = sp.symbols("t", positive=True)                              # t stands for tan(beta): identities become algebra
  cosrule = R**2 + b**2*t**2 - 2*R*b*t*sp.sin(th) - b**2*(1 + t**2)       # steps 6-7: cosine rule (1/cos^2 = 1 + tan^2)
  Rsol = b*t*sp.sin(th) + b*sp.sqrt(t**2*sp.sin(th)**2 + 1)       # its positive root R(theta): the circle itself
  assert sp.simplify(cosrule.subs(R, Rsol)) == 0                  # the root satisfies the cosine rule
  assert sp.simplify((R - b**2/R - 2*b*t*sp.sin(th)).subs(R, Rsol)) == 0  # steps 8-9: (14.6)
  cot2b = (1 - t**2)/(2*t)                                        # cot(2 beta) written with tan(beta)
  csc2b = (1 + t**2)/(2*t)                                        # csc(2 beta) written with tan(beta)
  assert sp.simplify(csc2b**2 - 1 - cot2b**2) == 0                # 1 + cot^2 = csc^2 survives the rewriting
  Y = 2*b*t*sp.sin(th)**2                                         # step 9: y from (14.6) joined to (14.4)
  X2 = 4*b**2*sp.cos(th)**2 + Y**2*sp.cos(th)**2/sp.sin(th)**2    # step 11: x^2 from (14.5) divided by sin^2
  arc = X2 + (Y + 2*b*cot2b)**2 - (2*b*csc2b)**2                  # steps 12-14: the circle of the arc, all on one side
  assert sp.simplify(arc) == 0                                    # every image point lies on that circle
  print(sp.simplify(Y.subs(th, sp.pi/2)))                         # step 15: 2*b*t, the arc's height 2b tan(beta)
  ```
- **What it means:** camber is how high the circle's centre sits; the two sharp ends are the two points ±b of the circle,
  where the map folds. **Fails when:** the circle does not pass through both ±b (then one end, or both, is rounded: C04's
  other cases).
- **Traps:** dividing by sin θ and cos θ without recovering the excluded points; taking the angle at O as θ instead of
  90° − θ; forgetting that the arc is traced twice; β is the circle's angle and the end tangent is 2β; the camber rule is
  first order in β.

### D04 · A circle shifted back gives a symmetric section: chord ≃ 4b, thickest at the quarter chord, $t_{\max}/c\simeq1.3\,e$ (14.8) — ★★★, 11 steps, in C04 (notebook · `zhukhovsky_airfoil_builder`)
- **Goal:** find the shape made by a circle whose centre is moved back along the real axis, and so see where thickness comes
  from. The book says the rules "can be proved"; this is the proof. · **Start:** the circle of radius $a\equiv b(1+e)$ (14.7) centred
  at ζ = −be, so that it still passes through ζ = b — "slightly larger than b, slightly to the left". · **Plan:** • write the
  circle with a small parameter • map it, keeping only first powers of e • read off x and y • find the chord • maximise the
  thickness. · **Tools:** binomial approximation with a complex small quantity (P342); Euler's formula (P45, reminder);
  maximum of a function by a zero derivative (gloss). · **Assumptions:** e ≪ 1 — terms of order e² are dropped (steps 3–5); each ≈ marks it.
- **Steps:**
  1. **did** Write the points of the shifted circle · **tex** $\zeta=-be+b(1+e)\,e^{i\theta}$ · **why** Centre −be plus radius b(1 + e) times a
     unit direction; θ is now the angle at the circle's *centre*. At θ = 0, ζ = b: the circle passes through the critical
     point. · **plain** A circle a little bigger than b, pushed a little to the left.
  2. **did** Factor out the unshifted point · **tex** $\zeta=b\,e^{i\theta}(1+\epsilon),\qquad\epsilon=e\,(1-e^{-i\theta})$ · **why** Divide the two e-terms by
     b e^{iθ}: e − e e^{−iθ}. This isolates everything small in one complex number ε, of size at most 2e. · **plain** The
     circle of radius b, times a factor close to one.
  3. **did** Invert with the binomial approximation · **tex** $\dfrac{b^2}{\zeta}=b\,e^{-i\theta}(1+\epsilon)^{-1}\approx b\,e^{-i\theta}(1-\epsilon)$ · **why** For small ε,
     1/(1 + ε) ≈ 1 − ε with an error of order ε², i.e. e² (P342). We keep first order because e ≪ 1. · **plain** The
     reflected point, to first order in e.
  4. **did** Add the two parts of the map · **tex** $z\approx2b\cos\theta+2i\,b\,\epsilon\sin\theta$ · **why** z = ζ + b²/ζ = b e^{iθ} + b e^{−iθ} +
     bε(e^{iθ} − e^{−iθ}); the first pair is 2b cos θ and the bracket is 2i sin θ (Euler). · **plain** The flat plate plus
     a small correction proportional to e.
  5. **did** Separate real and imaginary parts · **tex** $x\approx2b\cos\theta-2be\sin^2\theta,\qquad y\approx2be\,\sin\theta\,(1-\cos\theta)$ · **why** Insert
     ε = e(1 − cos θ + i sin θ): 2ibε sin θ = 2ibe sin θ(1 − cos θ) − 2be sin²θ. · **plain** The outline: x is nearly the
     plate's, y is small and proportional to e.
  6. **did** Find the chord · **tex** $x(0)=2b,\quad x(\pi)=-2b\ \Rightarrow\ c\simeq4b$ · **why** The tail is θ = 0, the nose θ = π; the correction
     −2be sin²θ vanishes at both. (Exactly, c = 4b(1 + e² + …): the error is second order.) · **plain** The chord is four
     times b, as for the plate.
  7. **did** Express the thickness · **tex** $t(\theta)=2\,y(\theta)\approx4be\,\sin\theta\,(1-\cos\theta)$ · **why** y(−θ) = −y(θ) while x(−θ) = x(θ):
     the section is symmetric about the x-axis, so the thickness at a station is twice the upper ordinate. · **plain**
     Thickness along the section, in terms of the circle angle.
  8. **did** Set the derivative to zero · **tex** $\dfrac{dy}{d\theta}=2be\,(\cos\theta-\cos2\theta)=0$ · **why** A smooth function is largest where
     its slope vanishes. Product rule: d/dθ[sin θ − sin θ cos θ] = cos θ − cos²θ + sin²θ = cos θ − cos 2θ. · **plain** The
     thickest station is where cos θ = cos 2θ.
  9. **did** Solve for the angle · **tex** $2\theta=2\pi-\theta\ \Rightarrow\ \theta=\dfrac{2\pi}{3}$ · **why** Two angles have equal cosines when they
     are equal (θ = 0: the tail, zero thickness) or add up to 2π; the second gives 120°. · **plain** The thickest point
     comes from the point of the circle 120° round from the tail.
  10. **did** Evaluate the maximum · **tex** $t_{\max}\approx4be\cdot\dfrac{\sqrt3}{2}\cdot\dfrac32=3\sqrt3\,be,\qquad\dfrac{t_{\max}}{c}\simeq\dfrac{3\sqrt3}{4}\,e\approx1.3\,e$ (14.8) · **why**
      sin 120° = √3/2 and 1 − cos 120° = 3/2; divide by c = 4b. The exact coefficient is 3√3/4 = 1.299. · **plain** The
      thickness ratio is 1.3 times the fractional enlargement of the circle.
  11. **did** Locate it along the chord · **tex** $x\!\left(\tfrac{2\pi}3\right)\approx2b\cos\tfrac{2\pi}3=-b$ · **why** Step 5 at θ = 120°, to leading order.
      The nose is at −2b and the tail at +2b, so x = −b is a quarter of the chord behind the nose. · **plain** The section is
      thickest at the quarter-chord point.
- **Result:** $c\simeq4b$ and $t_{\max}/c\simeq\tfrac{3\sqrt3}4e\approx1.3\,e$ at x = −b, two of the three rules of (14.8) — "move the circle back and the
  nose fattens".
- **Check:** units — t and c lengths; e dimensionless ✓. Limits — e → 0: y → 0, the plate. Number — e = 0.1: rule 0.1299;
  measured on the exact contour 0.1179 (9 % lower: the neglected e² terms), chord 4.033b against 4b. The tail: near θ = 0,
  y ≈ be θ³ while 2b − x ≈ bθ², so y ∝ (2b − x)^{3/2}: both surfaces leave the tail tangent to the chord — a cusp (D05).
- **Sympy check** (`check_src`; executed: prints `1.299038105676658` and passes in 1.0 s):
  ```python
  import sympy as sp                                              # symbolic algebra
  b, e = sp.symbols("b e", positive=True)                         # map constant, thickness parameter
  th = sp.symbols("theta", real=True)                             # angle round the circle, seen from its centre
  zeta = -b*e + b*(1 + e)*sp.exp(sp.I*th)                         # step 1: the shifted circle through zeta = b
  z = zeta + b**2/zeta                                            # the Zhukhovsky map (14.3)
  z1 = sp.series(z, e, 0, 2).removeO()                            # steps 2-4: keep terms up to first order in e
  x1 = sp.simplify(sp.re(sp.expand(z1.rewrite(sp.cos), complex=True)))   # step 5: real part
  y1 = sp.simplify(sp.im(sp.expand(z1.rewrite(sp.cos), complex=True)))   # step 5: imaginary part
  assert sp.simplify(x1 - (2*b*sp.cos(th) - 2*b*e*sp.sin(th)**2)) == 0   # x = 2b cos(theta) - 2be sin^2(theta)
  assert sp.simplify(y1 - 2*b*e*sp.sin(th)*(1 - sp.cos(th))) == 0        # y = 2be sin(theta)(1 - cos(theta))
  dy = sp.diff(y1, th)                                            # step 8: slope of the half-thickness along the circle
  assert sp.simplify(dy.subs(th, 2*sp.pi/3)) == 0                 # step 9: zero at theta = 120 degrees
  ymax = sp.simplify(y1.subs(th, 2*sp.pi/3))                      # step 10: the largest half-thickness
  assert sp.simplify(2*ymax/(4*b) - 3*sp.sqrt(3)/4*e) == 0        # step 10: t_max/c = (3 sqrt(3)/4) e
  assert sp.simplify(x1.subs(th, 2*sp.pi/3).subs(e, 0) + b) == 0  # step 11: it sits at x = -b, the quarter chord
  print(float(3*sp.sqrt(3)/4))                                    # 1.299...: the book's 1.3
  ```
- **What it means:** thickness is how far back the circle sits; because the circle still passes through ζ = b the tail stays
  sharp, and because −b is now inside the nose is round. **Fails when:** e is not small — at e = 0.25 the rule is a quarter
  too high (the "rule against measured" table).
- **Traps:** the centre is at −be, not −e; the maximum of sin θ(1 − cos θ) is at 120°, not 90°; thickness is twice y; the
  chord correction is second order (slip #4: the printed 1.3 is 3√3/4 rounded, and the rule is first order in e).

### D05 · The trailing edge is a cusp; a finite edge angle would need a stagnation point — ★★, 7 steps, in C04 (notebook)
- **Goal:** see why the one point ζ = b turns a smooth circle into a sharp tail, what kind of sharpness it is, and why that
  matters for the flow. The book states it with a figure. · **Start:** the map's local magnification, the derivative of
  $z=\zeta+\dfrac{b^2}{\zeta}$ (14.3) — "how the map stretches and turns a small step dζ". · **Plan:** • find where the stretching vanishes •
  expand the map there • read the angle doubling • compare with a wedge-shaped edge. · **Tools:** derivative of a power;
  Taylor expansion (P26, reminder); critical point of a map and angle doubling (ch06, recap in the block); corner flows
  (R07). · **Assumptions:** the circle passes through ζ = b and has ζ = −b inside it.
- **Steps:**
  1. **did** Differentiate the map · **tex** $\dfrac{dz}{d\zeta}=1-\dfrac{b^2}{\zeta^2}$ · **why** d(1/ζ)/dζ = −1/ζ². A small step dζ becomes
     dz = (dz/dζ) dζ: stretched by the size of this number and turned by its angle. · **plain** The map's local
     magnification.
  2. **did** Find where it vanishes · **tex** $\dfrac{dz}{d\zeta}=0\ \Longleftrightarrow\ \zeta=\pm b$ · **why** 1 − b²/ζ² = 0 gives ζ² = b². At such critical
     points the map is not angle-preserving. · **plain** Two special points, on the real axis at ±b.
  3. **did** Place them relative to the circle · **tex** $\zeta=-b\ \text{inside},\qquad\zeta=+b\ \text{on the circle}\ \to\ z=2b$ · **why** Our circles
     pass through +b by construction; shifting the centre left puts −b inside, where there is no fluid. · **plain** Only
     one critical point touches the fluid: it becomes the tail.
  4. **did** Expand the map near ζ = b · **tex** $z=2b+\dfrac{\delta^2}{b}+\dots,\qquad\delta=\zeta-b$ · **why** b²/(b + δ) = b(1 − δ/b + δ²/b² − …) by
     the geometric series; the terms linear in δ cancel, as they must where dz/dζ = 0. · **plain** Near the tail the
     map squares the small displacement.
  5. **did** Read off the angles · **tex** $\delta=re^{i\varphi}\ \Rightarrow\ z-2b=\dfrac{r^2}{b}\,e^{2i\varphi}$ · **why** Squaring a complex number doubles its
     angle. The circle leaves B in two opposite directions (φ and φ + π); their images differ by 2π: the same
     direction. · **plain** The two surfaces leave the tail tangent to each other: a cusp, with zero included angle.
  6. **did** Compare with a wedge-shaped edge · **tex** $\phi=\dfrac{Ud}{n}\Big(\dfrac rd\Big)^{n}\cos n\theta,\qquad q=U\Big(\dfrac rd\Big)^{n-1}$ · **why** The smooth
     flow leaving an edge of included angle 2π(1 − 1/n) is a corner flow (R07); its speed is the size of ∇φ. ·
     **plain** How fast the fluid moves at distance r from an edge of finite angle.
  7. **did** Let r → 0 · **tex** $n>1:\ q\to0;\qquad n=1:\ q=U$ · **why** A positive power of r vanishes at r = 0; n = 1 is the
     cusp (included angle 0), where the power is zero. · **plain** A finite-angle edge has a stagnation point; a cusp
     does not — the fluid leaves it at a finite speed.
- **Result:** $\dfrac{dz}{d\zeta}=1-\dfrac{b^2}{\zeta^2}=0$ at ζ = b, $z-2b\simeq(\zeta-b)^2/b$: the tail of a Zhukhovsky section is a cusp, and the Kutta
  flow passes it at a finite, non-zero speed.
- **Check:** units — δ²/b is a length ✓. Number — a 17° wedge is n = 1.05: speed ∝ r^{0.05}, zero only very close to the edge.
  Code — `AF.trailing_edge_angle` of our contour tends to 0 as the number of points grows.
- **What it means:** the sharp tail is built into the map; and because 1/(dz/dζ) is infinite there, the fluid speed is too
  unless the circle's flow stagnates at B — which is the Kutta condition (C05, C06). **Fails when:** a real wing cannot be
  built with a cusp (N33).
- **Traps:** the other critical point −b must be *inside* the circle; "cusp" describes the body, "stagnation point" the
  flow — at a cusp the Kutta flow does not stagnate.

### D06 · The flow past an off-centre circle at incidence — ★★, 7 steps, in C05 (notebook) — ours
- **Goal:** write the complex potential of a stream arriving at the angle α past the shifted circle of C04, with clockwise
  circulation Γ. The book writes only the centred case; everything in §14.5 needs this one. · **Start:**
  $w=U\Big(\zeta+\dfrac{a^2}{\zeta}\Big)+\dfrac{i\Gamma}{2\pi}\ln(\zeta/a)$ (6.68) — "stream along the real axis past a circle of radius a centred at the origin,
  clockwise circulation Γ". · **Plan:** • move the circle • turn the stream • tidy the logarithm • check the far field, the
  surface and the circulation. · **Tools:** shifting and rotating a complex potential (below); the complex logarithm (P155,
  reminder); "a streamline is a line of constant Im w" (ch06, recap). · **Assumptions:** ideal, steady, two-dimensional flow.
- **Steps:**
  1. **did** Move the circle to its centre μ · **tex** $w=U\Big[(\zeta-\mu)+\dfrac{a^2}{\zeta-\mu}\Big]+\dfrac{i\Gamma}{2\pi}\ln\dfrac{\zeta-\mu}{a}$ · **why** A function of ζ − μ is the
     same flow pattern translated by μ; it is still analytic, so still an ideal flow. · **plain** The same flow, now round
     the circle with centre μ.
  2. **did** Turn the stream by α · **tex** $w=U\Big[(\zeta-\mu)e^{-i\alpha}+\dfrac{a^2e^{i\alpha}}{\zeta-\mu}\Big]+\dfrac{i\Gamma}{2\pi}\ln\dfrac{(\zeta-\mu)e^{-i\alpha}}{a}$ · **why** Replace ζ − μ by
     (ζ − μ)e^{−iα}: every point is read in axes turned by +α, so the pattern itself turns by +α; a circle about μ is
     unchanged by a rotation about μ. · **plain** The stream now arrives tilted up by α.
  3. **did** Split the logarithm · **tex** $\ln\dfrac{(\zeta-\mu)e^{-i\alpha}}{a}=\ln\dfrac{\zeta-\mu}{a}-i\alpha$ · **why** ln(AB) = ln A + ln B; the extra term gives
     w the constant Γα/2π, and a constant in a potential changes no velocity. · **plain** The rotation only shifts the
     potential by a constant.
  4. **did** Drop the constant · **tex** $w=U\Big[(\zeta-\mu)e^{-i\alpha}+\dfrac{a^2e^{i\alpha}}{\zeta-\mu}\Big]+\dfrac{i\Gamma}{2\pi}\ln\dfrac{\zeta-\mu}{a}$ · **why** Only derivatives of w are
     velocities. This is the form we use. · **plain** Stream + its reflection in the circle + a swirl about the centre.
  5. **did** Check the far field · **tex** $\dfrac{dw}{d\zeta}\to Ue^{-i\alpha}=u-iv\ \Rightarrow\ (u,v)=(U\cos\alpha,\ U\sin\alpha)$ · **why** The second and third terms
     of dw/dζ fall off as 1/ζ² and 1/ζ. With dw/dζ = u − iv, the factor e^{−iα} means a velocity pointing *up* at +α. ·
     **plain** Far away: a uniform stream of speed U at the angle α.
  6. **did** Check the circle is a streamline · **tex** $\zeta-\mu=ae^{i\chi}:\quad w=2Ua\cos(\chi-\alpha)-\dfrac{\Gamma\chi}{2\pi}$ · **why** On the circle the two
     stream terms are complex conjugates and add to a real number; i ln(e^{iχ}) = −χ is real too. So Im w = 0 all
     round. · **plain** The stream function is constant on the circle: no flow through it.
  7. **did** Check the circulation · **tex** $\displaystyle\oint dw=\dfrac{i\Gamma}{2\pi}\,(2\pi i)=-\Gamma$ · **why** Going once counterclockwise, χ grows by 2π
     and ln(ζ − μ) by 2πi; the other terms return to their values. The change of Re w is the counterclockwise
     circulation. · **plain** Counterclockwise circulation −Γ: a clockwise circulation of Γ, as the book counts.
- **Result:** $w=U\Big[(\zeta-\mu)e^{-i\alpha}+\dfrac{a^2e^{i\alpha}}{\zeta-\mu}\Big]+\dfrac{i\Gamma}{2\pi}\ln\dfrac{\zeta-\mu}{a}$ (ours) — "the cylinder flow of Chapter 6, moved and turned".
- **Check:** units — every term m²/s ✓. Limits — μ = 0, α = 0 returns the Start, $w=U(\zeta+a^2/\zeta)+\frac{i\Gamma}{2\pi}\ln(\zeta/a)$ (6.68). Number — U = 30 m/s, α = 5°: far-field dw/dζ =
  29.89 − 2.61i. Code — `CONF.circle_flow_zeta(…, alpha=…, center=…)`: Im w constant on the circle to 1e-12.
- **What it means:** three ingredients — a stream, its image in the circle, a swirl — and one free number, Γ. **Fails when:**
  nothing here fixes Γ; that needs the Kutta condition (D09).
- **Traps:** rotating with e^{+iα} instead of e^{−iα} (check the far field); forgetting that the logarithm's argument is
  rotated too (it only adds a constant); Γ is clockwise, so the loop circulation is −Γ (trap T1).

### D07 · The velocity in the airfoil plane (14.9) — ★, 4 steps, in C05 (notebook · `kutta_condition`)
- **Goal:** get the fluid velocity round the airfoil from the known velocity round the circle. · **Start:** the complex
  potential has the same value at partner points, w(z) = w(ζ(z)), with $z=\zeta+\dfrac{b^2}{\zeta}$ (14.3) — "same potential, different
  coordinates". · **Plan:** • chain rule • derivative of the map • invert it • combine. · **Tools:** chain rule (R08; P49,
  reminder); derivative of a power. · **Assumptions:** ζ ≠ ±b (the two critical points).
- **Steps:**
  1. **did** Apply the chain rule · **tex** $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\,\dfrac{d\zeta}{dz}$ · **why** w depends on z only through ζ; for analytic functions
     the chain rule holds as for real ones. The left side is the velocity we want, u − iv. · **plain** Airfoil velocity =
     circle velocity × the map's rescaling.
  2. **did** Differentiate the map · **tex** $\dfrac{dz}{d\zeta}=1-\dfrac{b^2}{\zeta^2}=\dfrac{\zeta^2-b^2}{\zeta^2}$ · **why** Derivative of ζ + b²ζ^{−1}; then put both terms
     over ζ². The map gives z as a function of ζ, so this is the derivative we can take directly. · **plain** How much
     the map stretches a small step at ζ.
  3. **did** Invert it · **tex** $\dfrac{d\zeta}{dz}=\dfrac{1}{dz/d\zeta}=\dfrac{\zeta^2}{\zeta^2-b^2}$ · **why** For a one-to-one map the derivative of the inverse is
     the reciprocal of the derivative — valid where dz/dζ ≠ 0, i.e. away from ζ = ±b. · **plain** The rescaling needed in
     step 1.
  4. **did** Substitute into step 1 · **tex** $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\,\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9) · **why** Put the rescaling factor of step 3 into the chain rule of step 1; nothing else is needed, and the result
     holds wherever ζ ≠ ±b. · **plain** Multiply
     the circle's complex velocity by ζ²/(ζ² − b²) to get the airfoil's.
- **Result:** $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\,\dfrac{\zeta^2}{\zeta^2-b^2}$ — "the same flow, rescaled point by point".
- **Check:** units — the factor is dimensionless; both sides m/s ✓. Limits — ζ → ∞: the factor → 1, so the two far streams
  are identical (same U, same α); ζ → b: the factor → ∞. Number — plate, α = 0, ζ = 2b: (1 − ¼) × 4/3 = 1: the stream is
  undisturbed by an edge-on plate.
- **What it means:** speeds differ between the planes everywhere except at infinity; at the tail the airfoil's speed is
  infinite unless dw/dζ vanishes at ζ = b. **Fails when:** at ζ = ±b the formula is ∞ or 0/0 and needs a limit (N29).
- **Traps:** it is dζ/dz = 1/(dz/dζ) that enters, not dz/dζ; the velocities agree only at infinity (trap T6); dw/dz = u − iv,
  with a minus (trap T2).

### D08 · The circulation is the same in both planes — ★★, 5 steps, in C05 (notebook)
- **Goal:** show that the circle and the airfoil carry the same circulation, so that the lift theorem can be used on the
  airfoil with the circle's Γ. The book says "it can be shown". · **Start:** circulation as a loop integral,
  $\Gamma_{ccw}=\displaystyle\oint(u\,dx+v\,dy)$ — "add up the velocity along a closed path". · **Plan:** • write it with the complex velocity • notice
  the integrand is dw • do the same in the ζ-plane • compare partner curves. · **Tools:** real and imaginary parts of a
  product; closed-loop integral of a multi-valued function (P135, reminder). · **Assumptions:** the two closed curves are
  images of each other and enclose the body once.
- **Steps:**
  1. **did** Write the circulation with dw/dz · **tex** $\displaystyle\oint\dfrac{dw}{dz}\,dz=\oint(u-iv)(dx+i\,dy)=\oint(u\,dx+v\,dy)+i\oint(u\,dy-v\,dx)$ · **why** Multiply
     out; the real part is the circulation, the imaginary part the volume flux out of the loop (zero here: no
     source). · **plain** The loop integral of the complex velocity has the circulation as its real part.
  2. **did** Recognise an exact differential · **tex** $\dfrac{dw}{dz}\,dz=dw$ · **why** The definition of a derivative: the change of w
     along a small step. So the integral is just the total change of w once round the loop. · **plain** Circulation =
     how much the potential fails to return to its value after one circuit.
  3. **did** Do the same in the ζ-plane · **tex** $\dfrac{dw}{d\zeta}\,d\zeta=dw$ · **why** The same statement for the same function w,
     read in the other coordinate. · **plain** In the circle plane the circulation is also the change of w round the
     loop.
  4. **did** Follow partner curves · **tex** $\displaystyle\oint_{C_z}dw=\oint_{C_\zeta}dw$ · **why** As ζ goes once round $C_\zeta$, its image z goes once
     round $C_z$, and w has the same value at partner points at every moment. The total change is therefore the same. ·
     **plain** Same potential at partner points ⇒ same change round partner loops.
  5. **did** Take real parts · **tex** $\Gamma_{\text{airfoil}}=\Gamma_{\text{circle}}=\Gamma\ \text{(clockwise)}\ \Rightarrow\ L=\rho U\Gamma$ (6.62) holds in the z-plane ·
     **why** Step 1 applied to both sides of step 4; D06 step 7 gave the change −Γ, i.e. clockwise Γ. The far stream is
     also the same (D07), so the theorem applies. · **plain** The airfoil's lift is ρUΓ with the circle's Γ.
- **Result:** $\displaystyle\oint_{C_z}\dfrac{dw}{dz}\,dz=\oint_{C_\zeta}\dfrac{dw}{d\zeta}\,d\zeta$: the circulation does not change under the map.
- **Check:** units — m²/s both sides ✓. Number — with Γ = 10 m²/s the trapezoid sums round three partner contours give −10.000
  in both planes (C05, the cell under this derivation). Limit — far away z ≈ ζ, so ln ζ ≈ ln z and the swirl term is
  literally the same function.
- **What it means:** the map may squeeze and stretch the velocity, but never the loop integral of it. **Fails when:** the
  loops do not correspond one to one (e.g. a loop cutting through the body's image).
- **Traps:** the *velocity* is not the same in the two planes, only its loop integral (trap T6); orientation is preserved
  by the map (counterclockwise stays counterclockwise).

### D09 · The Kutta circulation (14.11) and the lift coefficient (14.12) — ★★, 10 steps, in C06 (notebook · `kutta_condition`)
- **Goal:** find the one circulation that keeps the trailing-edge speed finite, and the lift it gives. · **Start:**
  $u_\theta=-2U\sin\theta-\dfrac{\Gamma}{2\pi a}$ (14.10), the tangential speed on the circle — "stream sliding round the circle, plus the clockwise
  swirl"; θ measured from the diameter through Q that is aligned with the stream. · **Plan:** • find the angle of B (the
  point that becomes the tail) • demand zero speed there • solve for Γ • turn Γ into a lift, then a coefficient • simplify
  for thin sections and small angles. · **Tools:** geometry of the offset circle (`draw_circle_angles()`); sin(−x) = −sin x;
  the lift theorem (R04, D08); the section coefficient (P336); small angles (P100, reminder; P338). · **Assumptions:** the
  Kutta condition (step 4); thin section, e ≪ 1 and β ≪ 1 (step 9); small α + β in radians (step 10).
- **Steps:**
  1. **did** Locate B from the centre Q · **tex** $b-\mu=a\,e^{-i\beta}$ · **why** Our circle has centre μ = b − ae^{−iβ} (N25), so the
     point ζ = b lies at the angle −β seen from Q: β below the ξ-direction. · **plain** The tail point sits β below the
     horizontal through the centre.
  2. **did** Note the direction θ is measured from · **tex** $\theta=0\ \text{points along the stream, at angle }\alpha\text{ above the }\xi\text{-direction}$ · **why** The Start, $u_\theta=-2U\sin\theta-\Gamma/(2\pi a)$ (14.10), is the cylinder formula of Chapter 6, where θ is measured from the
     downstream direction; here the stream is tilted up by α. · **plain** The reference direction is tilted up by the angle of attack.
  3. **did** Express B in that angle · **tex** $\theta_B=-\beta-\alpha=-(\alpha+\beta)$ · **why** Angle from the ξ-direction (−β) minus the angle
     of the reference direction (+α). · **plain** Seen from the stream direction, the tail point is α + β below.
  4. **did** Impose the Kutta condition · **tex** $u_\theta(\theta_B)=0$ · **why** By $\frac{dw}{dz}=\frac{dw}{d\zeta}\frac{\zeta^2}{\zeta^2-b^2}$ (14.9) the airfoil's speed at the tail is the circle's speed there times an infinite factor; it stays finite only if the circle's flow stagnates at B (C02). ·
     **plain** The circle's rear stagnation point must sit exactly on B.
  5. **did** Insert the angle · **tex** $0=-2U\sin\!\big(-(\alpha+\beta)\big)-\dfrac{\Gamma}{2\pi a}=2U\sin(\alpha+\beta)-\dfrac{\Gamma}{2\pi a}$ · **why** Put θ_B into the speed formula of the Start; sin(−x) = −sin x turns the first term positive. · **plain** The swirl must cancel the stream's sliding speed at B.
  6. **did** Solve for Γ · **tex** $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11) · **why** Step 5 is one linear equation for Γ: move the swirl term to the other side and multiply by 2πa. This is the
     half-way result the lift is built on. · **plain** The clockwise
     circulation the sharp tail demands: proportional to speed, size and the sine of the absolute angle.
  7. **did** Turn it into a lift · **tex** $L=\rho U\Gamma=4\pi\rho U^2a\sin(\alpha+\beta)$ · **why** The lift theorem $L=\rho U\Gamma$ (6.62), valid on the
     airfoil with the circle's Γ (D08). Per unit span. · **plain** Lift per metre of span.
  8. **did** Make it a coefficient · **tex** $C_L=\dfrac{L}{\tfrac12\rho U^2c}=8\pi\,\dfrac ac\,\sin(\alpha+\beta)$ · **why** Divide by the dynamic pressure times the
     chord (P336). No approximation yet: this is exact for a Zhukhovsky section (ours). · **plain** The exact lift
     coefficient; it needs the ratio of the circle's radius to the chord.
  9. **did** Use the thin-section chord · **tex** $c\simeq4a\ \Rightarrow\ C_L\simeq2\pi\sin(\alpha+\beta)$ · **why** For e ≪ 1 and β ≪ 1, a ≃ b and $c\simeq4b$ (14.8), so a/c ≃ ¼. The book writes c ≃ 4a here, c ≃ 4b there (slip #8): both first order. · **plain** For a thin section the geometry drops out.
  10. **did** Use small angles · **tex** $C_L\simeq2\pi(\alpha+\beta)$ (14.12) · **why** sin x ≃ x for x ≪ 1 in radians; at 9° the error is
      0.4 %. · **plain** Lift coefficient = 2π times the absolute angle of attack, in radians.
- **Result:** $\Gamma=4\pi Ua\sin(\alpha+\beta)$ and $C_L\simeq2\pi(\alpha+\beta)$ — "a slope of 2π per radian, shifted by the camber".
- **Check:** units — Γ: (m/s)(m) ✓; $C_L$ dimensionless ✓. Limits — α = −β: Γ = 0, no lift; flat plate (a = b = c/4, β = 0):
  $C_L=2\pi\sin\alpha$ exactly. Number — our section at α = 5°, β = 4°: Γ = 15.96 m²/s, exact $C_L$ = 1.058, rule 0.987 (the exact
  value is 7.2 % higher: thickness makes a/c = 0.269 instead of 0.25). Code — the root-finder that knows no formula returns
  the same Γ to round-off.
- **What it means:** angle of attack and camber act through their sum; each degree is worth 0.11 of lift coefficient.
  **Fails when:** the boundary layer separates before the tail (stall): then the Kutta condition is not what the flow does.
- **Traps:** B is at −(α + β), not −(α − β): β is below the horizontal through Q *and* the reference is tilted up (trap T4);
  the radius is a, not b; three separate approximations (thin, small camber, small angle) — thickness alone raises the
  slope above 2π (slip #8); α and β in radians (trap T12).

### D10 · The normal-force coefficient is the area of the pressure loop — ★★, 7 steps, in C07 (notebook · `airfoil_pressure_lift`)
- **Goal:** turn the two pressure curves of a section into a force — and see why they are plotted as a loop. The book
  leaves this to an exercise. · **Start:** the force on a piece of surface of length ds, $d\mathbf F=-p\,\mathbf n\,ds$ per unit span —
  "pressure pushes on the surface against the outward normal **n**". · **Plan:** • write the normal along the contour • take
  the component across the chord • split the contour into upper and lower surface • subtract the ambient pressure • make
  it dimensionless • tilt to get the lift. · **Tools:** normal of a parametric curve (P92, reminder); closed line integral
  (P135, reminder); P337; N08. · **Assumptions:** shear stress neglected (pressure only); x along the chord, y across it.
- **Steps:**
  1. **did** Write the pressure force on an element · **tex** $d\mathbf F=-p\,\mathbf n\,ds$ · **why** Pressure acts normal to a surface, pushing
     *into* the body, hence the minus sign with **n** pointing outward. · **plain** Each bit of skin is pushed inward
     by p times its length.
  2. **did** Express the outward normal · **tex** $\mathbf n\,ds=(dy,\ -dx)$ · **why** Walking counterclockwise with the body on the
     left, the step is (dx, dy); turning it by −90° gives the vector pointing to the right of the walk: outward. ·
     **plain** The normal times the length, in terms of the step along the contour.
  3. **did** Take the y-component and add round the contour · **tex** $F_y=\displaystyle\oint p\,dx$ · **why** $dF_y=-p\,(-dx)=p\,dx$; sum over
     the closed contour, counterclockwise. Only the x-extent of an element matters for the force across the chord. ·
     **plain** The force across the chord is the loop integral of p dx.
  4. **did** Split into upper and lower surface · **tex** $F_y=\displaystyle\int_c^0p_u\,dx+\int_0^cp_l\,dx=-\int_0^cp_u\,dx+\int_0^cp_l\,dx$ · **why** Counterclockwise
     from the tail the walk goes over the *upper* surface with x decreasing from c to 0, then back along the lower
     surface with x increasing. Reversing limits changes the sign. · **plain** Lower-surface pressure pushes up,
     upper-surface pressure pushes down.
  5. **did** Subtract the ambient pressure · **tex** $F_y=\displaystyle\int_0^c\big[(p_l-p_\infty)-(p_u-p_\infty)\big]dx$ · **why** A constant pressure all
     round gives no net force ($\oint p_\infty\,dx=0$), so it may be subtracted from both curves. · **plain** Only pressure
     differences from ambient matter.
  6. **did** Divide by ½ρU²c · **tex** $C_y\equiv\dfrac{F_y}{\tfrac12\rho U^2c}=\displaystyle\int_0^1(C_{p,l}-C_{p,u})\,d\!\Big(\dfrac xc\Big)=\oint C_p\,d\!\Big(\dfrac xc\Big)$ · **why** Each
     bracket over ½ρU² is a pressure coefficient; dx/c = d(x/c). The last form is the loop integral again. · **plain**
     The normal-force coefficient is the area between the two $C_p$ curves.
  7. **did** Tilt from the chord to the stream · **tex** $C_L=C_y\cos\alpha-C_x\sin\alpha=\dfrac{C_y}{\cos\alpha}\ \ (\text{no drag})$ · **why** Lift is the
     force across the *stream*: $L=-F_x\sin\alpha+F_y\cos\alpha$ (N08). In ideal flow $D=F_x\cos\alpha+F_y\sin\alpha=0$, so
     $F_x=-F_y\tan\alpha$. · **plain** Divide the loop area by cos α to get the lift coefficient.
- **Result:** $F_y=-\displaystyle\int_0^cp_u\,dx+\int_0^cp_l\,dx$, $C_y=\displaystyle\oint C_p\,d\!\Big(\dfrac xc\Big)$ — "the loop's area is the force across the chord".
- **Check:** units — Pa × m = N/m ✓; $C_y$ dimensionless ✓. Limits — a symmetric section at α = 0 has identical curves: zero
  area. Number — our section at 5°: loop area 1.054; ÷ cos 5° = 1.058 = 2Γ/(Uc) ✓. A toy loop with −0.6 above and +0.4 below:
  area 1.0.
- **What it means:** plotting −$C_p$ upward makes the lift visible as an area, most of it from suction on the upper surface.
  **Fails when:** shear stress is not negligible (it adds mostly to drag), or the "no drag" step is used for a real,
  viscous flow — then $C_L=C_y\cos\alpha-C_x\sin\alpha$ needs the chordwise force too.
- **Traps:** the sign of the normal on the lower surface (handled by walking the contour one way); subtracting $p_\infty$
  changes nothing on a closed contour; x is along the chord, so this is the force normal to the chord — not the lift
  (trap T5).

### D11 · The strength of a shed filament — ★, 4 steps, in C08 (notebook · `downwash_induced_drag`)
- **Goal:** find how much vorticity leaves the wing where the bound circulation changes along the span. · **Start:**
  Helmholtz's theorem — "the strength of a vortex tube is the same all along it; a tube cannot end in the fluid" (R11). ·
  **Plan:** • balance the strength at a junction • expand the bound strength one step along the span • subtract • read the
  sign. · **Tools:** Helmholtz's theorem (R11); a first-order Taylor step (P26, reminder). · **Assumptions:** the shed
  filament leaves straight back (lifting-line assumption (iii)).
- **Steps:**
  1. **did** Balance the strength at a junction · **tex** $\Gamma(y)=\Gamma(y+dy)+d\Gamma_{\text{trail}}$ · **why** Between y and y + dy the bound
     vortex splits into the part that continues along the wing and the part that turns downstream; the total strength
     is conserved (Helmholtz). · **plain** What arrives equals what goes on plus what leaves.
  2. **did** Expand the bound strength · **tex** $\Gamma(y+dy)=\Gamma(y)+\dfrac{d\Gamma}{dy}\,dy$ · **why** First-order Taylor step for a small dy:
     value plus slope times step. · **plain** The bound vortex a little further out differs by its slope times dy.
  3. **did** Subtract · **tex** $d\Gamma_{\text{trail}}=-\dfrac{d\Gamma}{dy}\,dy$ · **why** Insert step 2 into step 1 and cancel Γ(y). · **plain** A
     filament leaves with the strength the bound vortex has lost.
  4. **did** Read the sign on each half · **tex** $y>0:\ \dfrac{d\Gamma}{dy}<0\Rightarrow d\Gamma_{\text{trail}}>0;\qquad\displaystyle\int_{-s/2}^{s/2}\!\!-\dfrac{d\Gamma}{dy}\,dy=\Gamma(-\tfrac s2)-\Gamma(\tfrac s2)=0$ ·
     **why** Γ falls toward the starboard tip and rises from the port tip; the integral is exact because Γ vanishes
     at both tips. · **plain** The two half-wings shed vorticity of opposite sense, and equal amounts.
- **Result:** $d\Gamma_{\text{trail}}=-\dfrac{d\Gamma}{dy}\,dy$ — "every drop in the bound circulation trails off as a vortex".
- **Check:** units — (m²/s per m) × m = m²/s ✓. Number — a staircase 10, 8, 5, 0 m²/s sheds 2, 3 and 5 m²/s: total 10 = the
  root value. Limit — constant Γ sheds nothing except at the tips (the horseshoe).
- **What it means:** the trailing vorticity is strongest where Γ changes fastest — near the tips; each half-sheet carries
  the mid-span circulation Γ₀ in total. **Fails when:** nothing here says the filaments stay straight — they roll up (N37).
- **Traps:** the sign (dΓ < 0 for dy > 0 on the starboard half gives a positive shed strength); the two half-sheets
  rotate in opposite senses.

### D12 · The downwash at the lifting line (14.13) — ★★, 7 steps, in C09 (notebook · `downwash_induced_drag`)
- **Goal:** add up the velocity that all the trailing filaments induce at one point of the wing. · **Start:** the speed
  induced by a straight vortex segment, $q=\dfrac{\Gamma}{4\pi d}(\cos\theta_a-\cos\theta_b)$, the straight-line form of the filament law
  $d\mathbf u=\dfrac{\Gamma\,d\mathbf l}{4\pi}\times\dfrac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ (5.17) — "d is the distance to the line, θ_a and θ_b the angles to its two ends". ·
  **Plan:** • specialise to a half-infinite line • insert the shed strength • fix the direction and sign • add all
  filaments • say what the singular integral means. · **Tools:** the segment formula (recap in C09; R12); D11; principal
  value (P344); improper integral (P144, reminder). · **Assumptions:** lifting-line assumptions (ii) and (iii) — a straight
  bound vortex; straight, parallel, non-interacting trailing filaments.
- **Steps:**
  1. **did** Take a filament that starts at the wing · **tex** $\theta_a=0,\ \theta_b=90^\circ:\quad q=\dfrac{\Gamma}{4\pi d}$ · **why** Seen from a point on
     the lifting line, one end of the filament is infinitely far downstream (angle 0) and the other is abreast of the
     point (90°); cos 0 − cos 90° = 1. · **plain** A half-infinite line induces half of what an infinite line would
     (Γ/2πd).
  2. **did** Insert the distance and the strength · **tex** $\lvert dw\rvert=\dfrac{\lvert d\Gamma/dy\rvert\,dy}{4\pi\,\lvert y-y_1\rvert}$ · **why** The filament leaves at
     y, the point is at y₁, so d = abs(y − y₁); its strength is $-\dfrac{d\Gamma}{dy}\,dy$ (D11). · **plain** The size of one
     filament's contribution.
  3. **did** Fix the direction · **tex** $dw(y_1)=\dfrac{-(d\Gamma/dy)\,dy}{4\pi\,(y-y_1)}$, w positive downward · **why** A starboard filament (positive strength) turns so that the air inboard of it moves down: positive for
     y > y₁. The signed denominator keeps the formula right on the other side too. · **plain** Filaments outboard of the
     point push it down.
  4. **did** Check the other side · **tex** $y<y_1,\ -\dfrac{d\Gamma}{dy}>0:\quad dw<0$ · **why** The point is now *outboard* of that
     starboard filament, where its swirl is upward. (Slip #3: the velocity is induced at y₁, wherever the filament
     is.) · **plain** A filament pushes down on one side of itself and up on the other.
  5. **did** Absorb the minus sign · **tex** $dw(y_1)=\dfrac{1}{4\pi}\,\dfrac{d\Gamma}{dy}\,\dfrac{dy}{(y_1-y)}$ · **why** −1/(y − y₁) = 1/(y₁ − y). · **plain** The
     same contribution, written as the book does.
  6. **did** Add all the filaments · **tex** $w(y_1)=\dfrac{1}{4\pi}\displaystyle\int_{-s/2}^{s/2}\dfrac{d\Gamma}{dy}\,\dfrac{dy}{(y_1-y)}$ (14.13) · **why** Velocities add
     (the filament law is linear); filaments leave from every y between the tips. The straight bound vortex induces
     no velocity on its own line. · **plain** The downwash is the sum of all the half-line contributions.
  7. **did** Say what the integral means at y = y₁ · **tex** $w(y_1)=\dfrac1{4\pi}\displaystyle\lim_{\epsilon\to0}\Big[\int_{-s/2}^{y_1-\epsilon}+\int_{y_1+\epsilon}^{s/2}\Big]\dfrac{d\Gamma}{dy}\dfrac{dy}{(y_1-y)}$ · **why** The
     integrand is infinite at y = y₁ with opposite signs on the two sides; the filaments just inboard and just outboard
     cancel, and the symmetric limit — the principal value — is finite (P344). · **plain** The singularity is
     harmless if we approach it evenly from both sides.
- **Result:** $w(y_1)=\dfrac{1}{4\pi}\displaystyle\int_{-s/2}^{s/2}\dfrac{d\Gamma}{dy}\,\dfrac{dy}{(y_1-y)}$ (a principal value) — "the air at the wing sinks at a speed set by how
  the circulation falls off along the whole span".
- **Check:** units — (m/s) × m / m = m/s ✓. Limits — one horseshoe (Γ drops only at the tips): w(0) = 2 × Γ/(4π s/2) = Γ/(πs);
  elliptic loading: Γ₁/2s at every station (computed by quadrature in C14). Number — Γ = 4π m²/s, s = 2 m: 2 m/s at
  mid-span.
- **What it means:** the downwash at a station depends on the *whole* loading, most strongly on how Γ changes nearby.
  **Fails when:** the wake rolls up or is not flat; near the tips of a wing whose Γ does not fall smoothly to zero the
  integral diverges (the horseshoe).
- **Traps:** the minus sign absorbed into (y₁ − y); w is positive downward (trap T5); the integral is a principal value
  (trap T10); the bound vortex contributes nothing on its own line; slip #3.

### D13 · Effective angle (14.14), the strip's forces, and lift and induced drag of the wing (14.15) — ★★, 8 steps, in C10 (notebook · `downwash_induced_drag`)
- **Goal:** find what the downwash does to the force on each strip of the wing, and add the strips up. · **Start:** the
  velocity triangle at a strip: the stream U along the flight path and the downwash w at right angles to it — "the air
  the strip actually meets". · **Plan:** • tilt angle • effective angle • lift theorem in the tilted stream • resolve the
  force across and along the flight path • integrate along the span. · **Tools:** relative velocity (P350); small angles
  (P100, reminder); resolving a force; the lift theorem per unit span (R04). · **Assumptions:** w ≪ U (step 2); each strip
  behaves as a two-dimensional section (large aspect ratio; step 4).
- **Steps:**
  1. **did** Find the tilt of the local stream · **tex** $\tan\varepsilon=\dfrac wU\ \Rightarrow\ \varepsilon=\tan^{-1}\dfrac wU$ · **why** U and w are the two legs of a
     right triangle; ε is the angle between the resultant and the flight path. (Slip #2: the book prints "tan" where
     the inverse tangent is meant.) · **plain** The air arrives tilted downward by ε.
  2. **did** Use w ≪ U · **tex** $\varepsilon\simeq\dfrac wU,\qquad\alpha_e=\alpha-\varepsilon=\alpha-\dfrac wU$ (14.14) · **why** For small angles tan⁻¹x ≃ x; the strip's
     chord makes the angle α with the flight path, so only α − ε with the air it meets. · **plain** The strip works at
     a smaller angle of attack than the wing is set at.
  3. **did** Write the resultant speed · **tex** $U_e=\sqrt{U^2+w^2},\qquad U_e\cos\varepsilon=U,\qquad U_e\sin\varepsilon=w$ · **why** Pythagoras and the
     definitions of cosine and sine in the same triangle — exact, no approximation. · **plain** The tilted stream's
     two components are exactly U and w.
  4. **did** Apply the lift theorem to the strip · **tex** $dL_e=\rho U_e\Gamma\,dy,\quad\perp U_e$ · **why** A two-dimensional section in a
     stream $U_e$ with circulation Γ feels ρ$U_e$Γ per unit span, at right angles to *that* stream ($L=\rho U\Gamma$ (6.62)). ·
     **plain** The force on the strip is perpendicular to the tilted stream, not to the flight path.
  5. **did** Resolve across the flight path · **tex** $dL=dL_e\cos\varepsilon=\rho\,U_e\cos\varepsilon\;\Gamma\,dy=\rho U\Gamma\,dy$ · **why** The force leans back by ε;
     its component at right angles to the flight path is what we call lift. Use step 3. · **plain** The strip's lift is
     ρUΓ dy, as if there were no downwash.
  6. **did** Resolve along the flight path · **tex** $dD_i=dL_e\sin\varepsilon=\rho\,U_e\sin\varepsilon\;\Gamma\,dy=\rho w\Gamma\,dy$ · **why** The other component
     points backward: a drag. Use step 3 again. · **plain** The strip's induced drag is ρwΓ dy.
  7. **did** Add the lift of all strips · **tex** $L=\displaystyle\int_{-s/2}^{s/2}\rho U\Gamma\,dy$ (14.15), first member · **why** Γ varies along the span,
     so the sum is an integral. · **plain** The wing's lift is ρU times the area under Γ(y).
  8. **did** Add the induced drag of all strips · **tex** $D_i=\displaystyle\int_{-s/2}^{s/2}\rho w\Gamma\,dy$ (14.15), second member · **why** Both w and Γ
     vary with y. · **plain** The wing's induced drag is ρ times the span integral of downwash × circulation.
- **Result:** $\alpha_e=\alpha-\dfrac wU$; $L=\displaystyle\int_{-s/2}^{s/2}\rho U\Gamma\,dy$, $D_i=\displaystyle\int_{-s/2}^{s/2}\rho w\Gamma\,dy$ — "U with Γ gives lift; w with Γ gives drag".
- **Check:** units — (kg/m³)(m/s)(m²/s)(m) = N ✓. Limits — infinite span: no trailing vortices, w = 0, $D_i$ = 0. Number —
  U = 50, w = 1 m/s, Γ = 20 m²/s, ρ = 1.2, dy = 1 m: dL = 1200 N, d$D_i$ = 24 N, ratio 0.02 = ε.
- **What it means:** induced drag is the lift vector leaning back; it exists in a fluid with no friction, and its work
  rate $D_iU$ is the kinetic energy left in the lengthening wake. **Fails when:** w is not small compared with U, or the
  strip is not nearly two-dimensional (low aspect ratio, near the tips).
- **Traps:** ε is tan⁻¹(w/U), not tan(w/U) (slip #2); $U_e\sin\varepsilon=w$ and $U_e\cos\varepsilon=U$ are exact — the only approximation
  is treating the flight path as the reference for "lift" and "drag" with a section law at $\alpha_e$; this drag is not
  friction.

### D14 · The sine series inside the downwash integral: the first two members of (14.19), corrected — ★★, 7 steps, in C11 (notebook)
- **Goal:** rewrite the downwash integral for a loading given as a sine series, so that the integral can be done term by
  term. This is where the book's one consequential misprint sits. · **Start:**
  $w(y_1)=\dfrac{1}{4\pi}\displaystyle\int_{-s/2}^{s/2}\dfrac{d\Gamma}{dy}\,\dfrac{dy}{(y_1-y)}$ (14.13) with $y=-(s/2)\cos\gamma$ and $\Gamma=\displaystyle\sum_{n=1}^\infty\Gamma_n\sin(n\gamma)$ (14.18) — "the
  downwash of a loading written as sines of the span angle". · **Plan:** • change the variable from y to γ • differentiate
  the series • rewrite the denominator • assemble • exchange sum and integral. · **Tools:** substitution in an integral (P106,
  reminder); the half-range sine series (P346); term-by-term differentiation (gloss). · **Assumptions:** the series may be
  differentiated and integrated term by term (true for the smooth loadings considered); the integral is a principal value.
- **Steps:**
  1. **did** Change the variable · **tex** $y=-\dfrac s2\cos\gamma\ \Rightarrow\ dy=\dfrac s2\sin\gamma\,d\gamma$ · **why** Differentiate; d(cos γ) = −sin γ dγ and the
     two minus signs cancel. γ grows with y. · **plain** A step along the span, in terms of a step in the angle.
  2. **did** Convert the limits · **tex** $y=-\dfrac s2\leftrightarrow\gamma=0,\qquad y=+\dfrac s2\leftrightarrow\gamma=\pi$ · **why** cos 0 = 1 and cos π = −1. ·
     **plain** Port tip to starboard tip is γ from 0 to π.
  3. **did** Rewrite the numerator · **tex** $\dfrac{d\Gamma}{dy}\,dy=d\Gamma=\dfrac{d\Gamma}{d\gamma}\,d\gamma$ · **why** Both are the change of Γ over the same small
     step; the factor dy/dγ of step 1 cancels inside it. · **plain** The shed strength does not care which
     coordinate we use.
  4. **did** Differentiate the series term by term · **tex** $\dfrac{d\Gamma}{d\gamma}=\displaystyle\sum_{n=1}^{\infty}n\,\Gamma_n\cos(n\gamma)$ · **why** d sin(nγ)/dγ = n cos(nγ);
     a convergent series of smooth terms may be differentiated term by term here. · **plain** Each sine sheds a
     cosine-shaped sheet, n times stronger for the n-th term.
  5. **did** Rewrite the denominator · **tex** $y_1-y=-\dfrac s2\cos\gamma_1+\dfrac s2\cos\gamma=\dfrac s2\,(\cos\gamma-\cos\gamma_1)$ · **why** Insert the definitions of y and y₁. **Order matters:** cos γ − cos γ₁. Slip #1: the book prints
     cos γ₁ − cos γ, which would flip the sign of the downwash. · **plain** The distance between filament and station, in the angles.
  6. **did** Assemble the integral · **tex** $w(y_1)=\dfrac{1}{2\pi s}\displaystyle\int_0^\pi\sum_{n=1}^{\infty}n\,\Gamma_n\,\dfrac{\cos(n\gamma)\,d\gamma}{\cos\gamma-\cos\gamma_1}$ · **why** Put steps 2–5 into the downwash integral of the Start: 1/(4π) ÷ (s/2) = 1/(2πs). This is the first member of
     the book's three-member downwash equation, with the denominator corrected. · **plain**
     The downwash as one integral over the span angle.
  7. **did** Exchange the sum and the integral · **tex** $w(y_1)=\dfrac{1}{2\pi s}\displaystyle\sum_{n=1}^{\infty}n\,\Gamma_n\int_0^\pi\dfrac{\cos(n\gamma)\,d\gamma}{\cos\gamma-\cos\gamma_1}$ · **why** The integral of a
     sum is the sum of the integrals; the $\Gamma_n$ are constants. This is the book's second member, corrected. · **plain**
     One standard integral per sine — the same for every wing.
- **Result:** $w(y_1)=\dfrac{1}{2\pi s}\displaystyle\sum_{n=1}^{\infty}n\,\Gamma_n\int_0^\pi\dfrac{\cos(n\gamma)\,d\gamma}{\cos\gamma-\cos\gamma_1}$ — "the downwash needs only the integrals $I_n(\gamma_1)$ of
  D15".
- **Check:** units — (m²/s)/m = m/s ✓. Sign — for the elliptic loading (n = 1 only, Γ₁ > 0) the integral is +π (D15), so w =
  +Γ₁/2s: downward, as it must be behind a lifting wing; the printed denominator gives −Γ₁/2s. Code — `LL.downwash_series(…,
  printed_middle=True)` returns the negative of the quadrature of the downwash integral.
- **What it means:** the change of variable crowds the stations toward the tips, where Γ changes fastest, and turns the
  awkward 1/(y₁ − y) into a form whose integral is known. **Fails when:** the loading has a jump (a deflected flap's
  edge): the series then converges slowly.
- **Traps:** the order in the denominator (slip #1); dy = (s/2) sin γ dγ cancels inside dΓ — it must not be inserted a
  second time; the limits are 0 to π, not −π to π.

### D15 · Glauert's integral, and the downwash of the series (14.19, last member) — ★★★, 11 steps, in C11 (notebook · `lifting_line_solver`) — ours
- **Goal:** evaluate the one integral the whole theory rests on. The book says only "the final equality comes from
  evaluating the integral". · **Start:** $I_n(\gamma_1)\equiv\displaystyle\int_0^\pi\dfrac{\cos(n\gamma)\,d\gamma}{\cos\gamma-\cos\gamma_1}$ (principal value), 0 < γ₁ < π — "the integral left over
  from D14, one for each n". · **Plan:** • find I₀ and I₁ directly • use a trigonometric identity to link three consecutive
  $I_n$ • guess a closed form and show that it obeys the same link and the same two starting values • conclude by
  induction. · **Tools:** principal value (P344); product-to-sum identity and three-term recursion (P347); proof by induction
  (P306, reminder); l'Hôpital's rule (P343) for the tips. · **Assumptions:** 0 < γ₁ < π; the integral is taken as a principal
  value.
- **Steps:**
  1. **did** Find an antiderivative for n = 0 · **tex** $\displaystyle\int\dfrac{d\gamma}{\cos\gamma-\cos\gamma_1}=\dfrac{1}{\sin\gamma_1}\ln\left\lvert\dfrac{\sin\frac{\gamma_1+\gamma}2}{\sin\frac{\gamma_1-\gamma}2}\right\rvert$ · **why** Differentiate
     the right side: ½cot((γ₁ + γ)/2) + ½cot((γ₁ − γ)/2) = sin γ₁/(cos γ − cos γ₁), using
     2 sin A sin B = cos(A − B) − cos(A + B). · **plain** A logarithm that blows up at γ = γ₁, like the integrand.
  2. **did** Cut out a symmetric gap · **tex** $I_0=\displaystyle\lim_{\epsilon\to0}\Big[\int_0^{\gamma_1-\epsilon}+\int_{\gamma_1+\epsilon}^{\pi}\Big]\dfrac{d\gamma}{\cos\gamma-\cos\gamma_1}$ · **why** The definition of the
     principal value (P344): the infinite areas on the two sides have opposite signs and must be removed evenly. ·
     **plain** Approach the singular point equally from both sides.
  3. **did** Evaluate I₀ · **tex** $I_0=\dfrac{1}{\sin\gamma_1}\displaystyle\lim_{\epsilon\to0}\ln\dfrac{\sin(\gamma_1-\epsilon/2)}{\sin(\gamma_1+\epsilon/2)}=0$ · **why** Insert the limits in step 1: the
     end points give ln 1 = 0 (at γ = π both sines equal cos(γ₁/2)); the two logarithms of sin(ε/2) cancel. · **plain**
     With no cosine on top, the two sides cancel exactly.
  4. **did** Evaluate I₁ · **tex** $I_1=\displaystyle\int_0^\pi\dfrac{(\cos\gamma-\cos\gamma_1)+\cos\gamma_1}{\cos\gamma-\cos\gamma_1}\,d\gamma=\pi+\cos\gamma_1\,I_0=\pi$ · **why** Add and subtract cos γ₁ in
     the numerator; the first part is the integral of 1 over (0, π), the second is a constant times I₀ = 0. ·
     **plain** The first non-trivial integral is π, whatever γ₁ is.
  5. **did** Recall the product-to-sum identity · **tex** $\cos(n+1)\gamma+\cos(n-1)\gamma=2\cos\gamma\,\cos n\gamma$ · **why** Add the expansions of
     cos(nγ ± γ); the sine terms cancel (P347). It ties three consecutive numerators together. · **plain** Two
     neighbouring cosines add up to the middle one times 2 cos γ.
  6. **did** Turn it into a recursion · **tex** $I_{n+1}+I_{n-1}=2\displaystyle\int_0^\pi\cos n\gamma\,d\gamma+2\cos\gamma_1\,I_n=2\cos\gamma_1\,I_n\quad(n\ge1)$ · **why** Divide the identity by (cos γ − cos γ₁), integrate, and split 2 cos γ as in step 4; the integral of cos nγ
     over (0, π) vanishes for n ≥ 1. · **plain** Each $I_n$ follows from the two before it.
  7. **did** Propose a closed form · **tex** $S_n\equiv\pi\,\dfrac{\sin n\gamma_1}{\sin\gamma_1}:\qquad S_0=0,\quad S_1=\pi$ · **why** A guess suggested by the first
     few values; it matches both starting values of steps 3 and 4. · **plain** A candidate that starts right.
  8. **did** Show the candidate obeys the recursion · **tex** $S_{n+1}+S_{n-1}=\dfrac{\pi}{\sin\gamma_1}\big[\sin(n+1)\gamma_1+\sin(n-1)\gamma_1\big]=2\cos\gamma_1\,S_n$ · **why** The
     companion identity sin(A + B) + sin(A − B) = 2 sin A cos B with A = nγ₁, B = γ₁. · **plain** The candidate
     satisfies the same three-term rule.
  9. **did** Conclude by induction · **tex** $I_n=\pi\,\dfrac{\sin n\gamma_1}{\sin\gamma_1}\quad\text{for all }n\ge0$ · **why** Two sequences with the same two
     starting values and the same rule for each next term are equal term by term (P306). · **plain** Glauert's
     integral.
  10. **did** Insert into D14's result · **tex** $w(y_1)=\dfrac{1}{2\pi s}\displaystyle\sum_{n=1}^{\infty}n\Gamma_n\,\pi\dfrac{\sin n\gamma_1}{\sin\gamma_1}=\dfrac{1}{2s}\sum_{n=1}^{\infty}n\,\Gamma_n\,\dfrac{\sin(n\gamma_1)}{\sin\gamma_1}$ (14.19),
      last member · **why** Substitute step 9 for each integral; π cancels. · **plain** The downwash of a sine series is
      another sine series, divided by sin γ₁.
  11. **did** Check the tips · **tex** $\displaystyle\lim_{\gamma_1\to0}\dfrac{\sin n\gamma_1}{\sin\gamma_1}=n,\qquad\lim_{\gamma_1\to\pi}\dfrac{\sin n\gamma_1}{\sin\gamma_1}=(-1)^{n+1}n$ · **why** Both are 0/0;
      l'Hôpital's rule gives n cos(nγ₁)/cos γ₁ at the two ends (P343). · **plain** Each term stays finite at the tips,
      but grows like n² — a hint that the series converges slowly there.
- **Result:** $\displaystyle\int_0^\pi\dfrac{\cos n\gamma\,d\gamma}{\cos\gamma-\cos\gamma_1}=\pi\,\dfrac{\sin n\gamma_1}{\sin\gamma_1}$, hence $w(y_1)=\dfrac{1}{2s}\displaystyle\sum_{n=1}^{\infty}n\,\Gamma_n\,\dfrac{\sin(n\gamma_1)}{\sin\gamma_1}$ — "the singular integral,
  done once and for all".
- **Check:** units — the integral is dimensionless; w in m/s ✓. Limits — n = 1: w = Γ₁/2s, the same at every station.
  Number — γ₁ = 60°: I₀ … I₄ = 0, π, π, 0, −π; quadrature with the singular part subtracted gives the same to 15 digits.
- **Sympy check** (`check_src`; executed: prints `3` and passes in 0.4 s):
  ```python
  import sympy as sp                                              # symbolic algebra
  g, g1 = sp.symbols("gamma gamma_1", positive=True)              # span angle of the filament and of the station
  n = sp.symbols("n", integer=True, positive=True)                # index of the sine term
  ident = sp.cos((n + 1)*g) + sp.cos((n - 1)*g) - 2*sp.cos(g)*sp.cos(n*g)        # step 5: product-to-sum identity
  assert sp.simplify(sp.expand_trig(ident)) == 0
  S = lambda k: sp.pi*sp.sin(k*g1)/sp.sin(g1)                     # step 7: the candidate closed form
  rec = S(n + 1) + S(n - 1) - 2*sp.cos(g1)*S(n)                   # step 8: does it obey the recursion?
  assert sp.simplify(sp.expand_trig(rec)) == 0
  assert S(0) == 0 and sp.simplify(S(1) - sp.pi) == 0             # steps 3-4: the two starting values match
  X, X1 = sp.symbols("X X_1")                                     # X = cos(gamma), X1 = cos(gamma_1)
  for k in range(1, 6):                                           # direct evaluation for n = 1..5
      Tk = sp.chebyshevt(k, X)                                    # cos(k gamma) is a polynomial T_k in cos(gamma)
      q = sp.cancel((Tk - Tk.subs(X, X1))/(X - X1))               # the singular factor divides out exactly
      Ik = sp.integrate(q.subs(X, sp.cos(g)), (g, 0, sp.pi))      # ordinary integral; the remainder is cos(k g1)*I_0 = 0
      assert sp.simplify(Ik.subs(X1, sp.cos(g1)) - sp.pi*sp.chebyshevu(k - 1, sp.cos(g1))) == 0   # = pi sin(k g1)/sin(g1)
  print(sp.limit(sp.sin(3*g1)/sp.sin(g1), g1, 0))                 # step 11: the tip value is n (here 3)
  ```
- **What it means:** every sine shape of the loading produces a downwash of the same sine shape divided by sin γ — so a
  single sine (the ellipse) gives a uniform downwash, and higher sines give more downwash toward the tips. **Fails when:**
  γ₁ is exactly 0 or π (use step 11), which is why collocation stations avoid the tips.
- **Traps:** I₀ = 0 holds only as a principal value; the recursion needs **two** starting values; showing that the
  candidate satisfies the recursion *is* the proof — no further integration is needed; the printed form of the first two
  members has the opposite sign (slip #1).

### D16 · The lifting-line equation (14.16) → (14.17) → (14.20) — ★★, 8 steps, in C12 (notebook)
- **Goal:** one equation from which the circulation along any wing can be computed. · **Start:** the section law at the
  effective angle, $C_L=K\Big[\alpha-\dfrac{w(y)}{U}+\beta(y)\Big]$ (14.16) — "each strip lifts like a two-dimensional section, at the angle the downwash
  leaves it"; K is the section's lift-curve slope per radian, −β its zero-lift angle. · **Plan:** • express Γ through $C_L$ •
  combine into an equation for Γ(y) with w inside • move the w-term to the unknown side • insert both sine series • collect
  • truncate. · **Tools:** substitution; the section coefficient (P336); the sine series (P346) and its downwash (D15); integral
  equation (P345). · **Assumptions:** the section law is linear (no stall); lifting-line assumptions (i)–(iii) of N39.
- **Steps:**
  1. **did** State the section law on a finite wing · **tex** $C_L=K\Big[\alpha-\dfrac{w(y)}{U}+\beta(y)\Big]$ (14.16) · **why** It generalises
     $C_L\simeq2\pi(\alpha+\beta)$ (14.12): slope K instead of 2π, and the geometric angle reduced by ε ≃ w/U
     ($\alpha_e=\alpha-\dfrac wU$ (14.14)). · **plain** The local lift coefficient depends on the local downwash.
  2. **did** Relate circulation and lift coefficient · **tex** $\rho U\Gamma=C_L\,\tfrac12\rho U^2c\ \Rightarrow\ \Gamma=\tfrac12U\,c\,C_L$ · **why** Two expressions for
     the lift per unit span of the strip: the lift theorem $L=\rho U\Gamma$ (6.62) and the definition of the section
     coefficient. · **plain** A strip's circulation is half its chord times speed times lift coefficient.
  3. **did** Combine · **tex** $\Gamma(y)=\dfrac K2\,U\,c(y)\Big[\alpha-\dfrac{w(y)}{U}+\beta(y)\Big]$ (14.17) · **why** Insert step 1 into step 2. Together with
     the downwash integral $w(y_1)=\dfrac1{4\pi}\displaystyle\int\dfrac{d\Gamma}{dy}\dfrac{dy}{(y_1-y)}$ (14.13) this is an integral equation: Γ appears on the left and,
     through w, under an integral on the right (P345). · **plain** Circulation determines downwash determines
     circulation.
  4. **did** Move the downwash term to the left · **tex** $\Gamma+\dfrac K2\,c\,w=\dfrac K2\,U\,c\,(\alpha+\beta)$ · **why** Multiply out the bracket:
     (K/2)Uc × (−w/U) = −(K/2)c w; add it to both sides. Everything unknown is now on the left. · **plain** Unknowns on
     one side, the wing's geometry on the other.
  5. **did** Insert the series for Γ · **tex** $\Gamma=\displaystyle\sum_{n=1}^{\infty}\Gamma_n\sin(n\gamma)$ (14.18) at the station γ · **why** We look for the loading
     as a sum of sines of the span angle; the unknowns become the numbers $\Gamma_n$. · **plain** The loading at this
     station, in the unknown amplitudes.
  6. **did** Insert the series for w at the same station · **tex** $\dfrac K2\,c\,w=\displaystyle\sum_{n=1}^{\infty}\dfrac{nKc}{4s\sin\gamma}\,\Gamma_n\sin(n\gamma)$ · **why**
     $w=\dfrac1{2s}\sum n\Gamma_n\dfrac{\sin(n\gamma)}{\sin\gamma}$ (14.19, last member) with γ₁ renamed γ; multiply by Kc/2. · **plain** The
     downwash term, in the same amplitudes.
  7. **did** Collect on sin(nγ) · **tex** $\dfrac K2\,U\,c\,(\alpha+\beta)=\displaystyle\sum_{n=1}^{\infty}\Big(1+\dfrac{nKc}{4s\sin\gamma}\Big)\Gamma_n\sin(n\gamma)$ (14.20) · **why** Add steps 5 and
     6 term by term and put them into step 4. · **plain** One equation that must hold at every station γ.
  8. **did** Truncate and collocate · **tex** $\displaystyle\sum_{n=1}^{N}\mathsf A_{jn}\Gamma_n=r_j,\quad\mathsf A_{jn}=\Big(1+\dfrac{nKc_j}{4s\sin\gamma_j}\Big)\sin n\gamma_j,\quad r_j=\dfrac K2Uc_j(\alpha_j+\beta_j)$ · **why**
     Keep N terms and require step 7 at N stations $\gamma_j=j\pi/(N+1)$: N linear equations for N unknowns (ours,
     formalising the book's description). · **plain** A small matrix problem for the sine amplitudes.
- **Result:** $\dfrac K2\,Uc(\alpha+\beta)=\displaystyle\sum_{n=1}^{\infty}\Big(1+\dfrac{nKc}{4s\sin\gamma}\Big)\Gamma_n\sin(n\gamma)$ — "any planform, any twist: one linear solve".
- **Check:** units — left (m/s)(m) = m²/s; right: the bracket is dimensionless (c/s), $\Gamma_n$ in m²/s ✓. Limits — s → ∞ at
  fixed chord: the bracket → 1 and Γ = (K/2)Uc(α + β), the two-dimensional result; elliptic chord c = c₀ sin γ: the sin γ
  cancels and n = 1 alone satisfies it (D20). Number — rectangle Λ = 6, K = 2π, α = 5°, N = 1: Γ₁ = 0.2742/1.262 = 0.217 Uc.
- **What it means:** the "1" is the two-dimensional response; the second term is the wing talking to itself through its
  wake — larger for higher n and for broader chords. **Fails when:** stalled sections (non-linear law), swept or stubby
  wings, a rolled-up wake (N68).
- **Traps:** −β is the zero-lift angle here (trap T3); K per radian (trap T12); c, α, β, K may all depend on γ, so this
  is not a Fourier series whose coefficients could be found by projection — collocation is needed; sin γ = 0 at the tips
  (trap T10).

### D17 · Lift from the first coefficient (14.21) — ★★, 6 steps, in C13 (notebook)
- **Goal:** get the wing's total lift from the sine amplitudes. The book gives the result only. · **Start:**
  $L=\displaystyle\int_{-s/2}^{s/2}\rho U\Gamma\,dy$ (14.15) with $\Gamma=\displaystyle\sum_{n=1}^\infty\Gamma_n\sin(n\gamma)$ (14.18) — "lift is ρU times the area under the loading". · **Plan:** •
  change variable to γ • insert the series • integrate term by term • use orthogonality. · **Tools:** substitution in an
  integral (P106, reminder); orthogonality of sines on the half range (P346). · **Assumptions:** ρ and U uniform along the
  span.
- **Steps:**
  1. **did** Change the variable · **tex** $L=\rho U\,\dfrac s2\displaystyle\int_0^\pi\Gamma\,\sin\gamma\,d\gamma$ · **why** $y=-(s/2)\cos\gamma$ gives dy = (s/2) sin γ dγ and
     limits 0 to π (D14 steps 1–2). · **plain** The area under the loading, measured in the span angle.
  2. **did** Insert the series · **tex** $L=\rho U\,\dfrac s2\displaystyle\int_0^\pi\Big[\sum_{n=1}^{\infty}\Gamma_n\sin(n\gamma)\Big]\sin\gamma\,d\gamma$ · **why** Replace Γ by its series $\Gamma=\sum\Gamma_n\sin(n\gamma)$ (14.18); we want the lift in terms of the amplitudes the
     solver returns. · **plain** Every sine of the loading is multiplied by sin γ.
  3. **did** Exchange sum and integral · **tex** $L=\rho U\,\dfrac s2\displaystyle\sum_{n=1}^{\infty}\Gamma_n\int_0^\pi\sin(n\gamma)\,\sin\gamma\,d\gamma$ · **why** The integral of a sum is
     the sum of the integrals. · **plain** One integral per term.
  4. **did** Use orthogonality · **tex** $\displaystyle\int_0^\pi\sin(n\gamma)\,\sin\gamma\,d\gamma=\dfrac\pi2\ (n=1),\qquad0\ (n\ne1)$ · **why** On the half range two
     different sines integrate to zero and a sine with itself to π/2 (P346); here the second sine is the n = 1 one,
     supplied by dy. · **plain** Only the first term survives.
  5. **did** Keep the surviving term · **tex** $L=\rho U\,\dfrac s2\cdot\dfrac\pi2\,\Gamma_1$ · **why** Insert step 4 into step 3. · **plain** The lift
     contains Γ₁ and nothing else.
  6. **did** Tidy up · **tex** $L=\dfrac{\pi s}{4}\,\rho U\Gamma_1$ (14.21) · **why** Multiply the two numerical factors, (s/2)(π/2) = πs/4, to reach the form the book prints. · **plain** Lift = ρU × (the area of
     a half-ellipse of height Γ₁ and width s).
- **Result:** $L=\dfrac{\pi s}{4}\,\rho U\Gamma_1$ — "the lift needs only the first coefficient".
- **Check:** units — (m)(kg/m³)(m/s)(m²/s) = N ✓. Limits — elliptic loading: ρU × area of a half-ellipse = ρU(π/4)sΓ₁ ✓.
  Number — Γ₁ = 12 m²/s, s = 12 m, ρ = 1.2, U = 50: 6786 N, whatever Γ₃ is.
- **What it means:** the higher sines redistribute lift along the span without adding any. **Fails when:** never, within
  the series; it is an exact consequence of the two Start equations.
- **Traps:** the sin γ that selects n = 1 comes from dy — forgetting it gives the wrong selection; the half-range
  integral is π/2, not π.

### D18 · Induced drag from the coefficients (14.22) — ★★★, 9 steps, in C13 (notebook · `lifting_line_solver`)
- **Goal:** get the wing's induced drag from the sine amplitudes — and see that every amplitude adds to it. The book
  gives the result only. · **Start:** $D_i=\displaystyle\int_{-s/2}^{s/2}\rho w\Gamma\,dy$ (14.15) with $\Gamma=\displaystyle\sum_{m}\Gamma_m\sin(m\gamma)$ (14.18) and
  $w=\dfrac{1}{2s}\displaystyle\sum_{n}n\,\Gamma_n\,\dfrac{\sin(n\gamma)}{\sin\gamma}$ (14.19, last member) — "drag is ρ times the span integral of downwash × circulation". ·
  **Plan:** • change variable • insert w and cancel sin γ • insert Γ with a second index • multiply the two sums • integrate
  term by term • collapse with orthogonality. · **Tools:** substitution in an integral (P106, reminder); a double sum
  collapsed by orthogonality (P348); orthogonality on the half range (P346). · **Assumptions:** the series are truncated at
  N terms (or converge well enough to be multiplied term by term).
- **Steps:**
  1. **did** Change the variable · **tex** $D_i=\rho\,\dfrac s2\displaystyle\int_0^\pi w\,\Gamma\,\sin\gamma\,d\gamma$ · **why** dy = (s/2) sin γ dγ, limits 0 to π (D14 steps
     1–2). · **plain** The drag integral in the span angle.
  2. **did** Insert the downwash series · **tex** $w\,\sin\gamma=\dfrac{1}{2s}\displaystyle\sum_{n}n\,\Gamma_n\sin(n\gamma)$ · **why** Multiply the downwash series of the Start by sin γ: its 1/sin γ is cancelled by the sin γ of dy — once. · **plain** The awkward denominator disappears.
  3. **did** Simplify the prefactor · **tex** $D_i=\dfrac\rho4\displaystyle\int_0^\pi\Big[\sum_{n}n\,\Gamma_n\sin(n\gamma)\Big]\Gamma\,d\gamma$ · **why** ρ(s/2) × 1/(2s) = ρ/4; the span
     has cancelled. · **plain** The induced drag does not contain the span explicitly — only the amplitudes.
  4. **did** Insert the circulation with its own index · **tex** $D_i=\dfrac\rho4\displaystyle\int_0^\pi\Big[\sum_{n}n\,\Gamma_n\sin(n\gamma)\Big]\Big[\sum_{m}\Gamma_m\sin(m\gamma)\Big]d\gamma$ ·
     **why** The two sums run independently, so they need different index letters; using n for both would silently
     drop the cross terms (P348). · **plain** A product of two sine series.
  5. **did** Multiply out · **tex** $D_i=\dfrac\rho4\displaystyle\int_0^\pi\sum_{n}\sum_{m}n\,\Gamma_n\Gamma_m\,\sin(n\gamma)\sin(m\gamma)\,d\gamma$ · **why** Every term of the first sum times
     every term of the second: N² products. · **plain** All pairings of a downwash sine with a loading sine.
  6. **did** Exchange sums and integral · **tex** $D_i=\dfrac\rho4\displaystyle\sum_{n}\sum_{m}n\,\Gamma_n\Gamma_m\int_0^\pi\sin(n\gamma)\sin(m\gamma)\,d\gamma$ · **why** Finite sums may be
     taken outside the integral; the amplitudes are constants. · **plain** One simple integral per pairing.
  7. **did** Use orthogonality · **tex** $\displaystyle\int_0^\pi\sin(n\gamma)\sin(m\gamma)\,d\gamma=\dfrac\pi2\,\delta_{nm}$ · **why** Different sines integrate to zero over
     the half range, equal ones to π/2 (P346); $\delta_{nm}$ is 1 if n = m and 0 otherwise. · **plain** Only pairings of
     a sine with itself contribute.
  8. **did** Collapse the double sum · **tex** $\displaystyle\sum_{n}\sum_{m}n\,\Gamma_n\Gamma_m\,\dfrac\pi2\,\delta_{nm}=\dfrac\pi2\sum_{n}n\,\Gamma_n^2$ · **why** The delta removes the
     sum over m, setting m = n (P348). · **plain** N² pairings shrink to N squares.
  9. **did** Assemble · **tex** $D_i=\dfrac\pi8\,\rho\displaystyle\sum_{n=1}^{N}n\,\Gamma_n^2$ (14.22) · **why** (ρ/4)(π/2) = πρ/8. (Slip #7: the book writes the upper limit N here and ∞ in its series for Γ and w; this line is for the truncated solution.) · **plain** Each sine
     adds n times its amplitude squared — never a negative amount.
- **Result:** $D_i=\dfrac\pi8\,\rho\displaystyle\sum_{n=1}^{N}n\,\Gamma_n^2$ — "every coefficient adds a positive amount of drag".
- **Check:** units — (kg/m³)(m⁴/s²) = kg m/s² = N ✓ (ρΓ² is a force). Limits — one term: the elliptic value $D_i=\frac\pi8\rho\Gamma_1^2$ (14.25). Number — Γ₁ = 12, Γ₃ = 1 m²/s, ρ = 1.2: (π/8) × 1.2 × (144 + 3) = 69.3 N. Code — the trapezoid rule on the Start integral with w from the series agrees to 1e-5.
- **Sympy check** (`check_src`; executed: prints `0 pi/2` and passes in 1.1 s):
  ```python
  import sympy as sp                                              # symbolic algebra
  g, s, rho = sp.symbols("gamma s rho", positive=True)            # span angle, span, density
  G = sp.symbols("Gamma1:5", real=True)                           # four coefficients Gamma_1..Gamma_4
  Gam = sum(G[k]*sp.sin((k + 1)*g) for k in range(4))             # (14.18) truncated at N = 4
  w = sum((k + 1)*G[k]*sp.sin((k + 1)*g) for k in range(4))/(2*s*sp.sin(g))     # (14.19), last member
  integrand = rho*w*Gam*(s/2)*sp.sin(g)                           # steps 1-3: rho w Gamma dy with dy = (s/2) sin(gamma) d(gamma)
  Di = sp.integrate(sp.expand(sp.simplify(integrand)), (g, 0, sp.pi))            # steps 4-8: all sixteen products, integrated
  target = sp.pi/8*rho*sum((k + 1)*G[k]**2 for k in range(4))     # (14.22)
  assert sp.simplify(Di - target) == 0                            # the cross terms have vanished
  print(sp.integrate(sp.sin(2*g)*sp.sin(3*g), (g, 0, sp.pi)), sp.integrate(sp.sin(3*g)**2, (g, 0, sp.pi)))   # step 7: 0 and pi/2
  ```
- **What it means:** for a given lift (a given Γ₁) the drag is least when all other amplitudes vanish; the excess is
  measured by $e_{\text{span}}=\Gamma_1^2/\sum n\Gamma_n^2\le1$ (N54). **Fails when:** the wake is not flat (winglets, non-planar
  wings can beat the planar optimum).
- **Traps:** the sin γ of dy cancels the 1/sin γ of w — do not cancel it twice; two different summation indices are
  needed; the factor n comes from w, once; N against ∞ as the upper limit (slip #7).

### D19 · The elliptic optimum: (14.23), (14.24), (14.25), (14.26) — ★★, 9 steps, in C14 (notebook · `elliptic_loading_aspect_ratio`)
- **Goal:** find the loading that carries a given lift on a given span with the least induced drag, and what that least
  drag is. · **Start:** $L=\dfrac{\pi s}{4}\rho U\Gamma_1$ (14.21) and $D_i=\dfrac\pi8\rho\displaystyle\sum n\Gamma_n^2$ (14.22) — "lift from the first amplitude, drag from all of
  them". · **Plan:** • fix the lift • drop everything that only costs drag • translate the remaining sine back to y • find its
  downwash and drag • express the drag by the lift • make it a coefficient. · **Tools:** a sum of non-negative terms;
  sin γ = √(1 − cos²γ); substitution; section against wing coefficients (P336). · **Assumptions:** fixed lift **and span**; a
  flat wake (lifting-line assumptions).
- **Steps:**
  1. **did** Fix the lift · **tex** $\Gamma_1=\dfrac{4L}{\pi s\,\rho U}$ · **why** Solve $L=\frac{\pi s}4\rho U\Gamma_1$ (14.21) for Γ₁: for a given lift, span, density and speed
     the first amplitude is not ours to choose. · **plain** The lift pins down the first sine.
  2. **did** Separate the drag into two parts · **tex** $D_i=\dfrac\pi8\rho\,\Gamma_1^2+\dfrac\pi8\rho\displaystyle\sum_{n\ge2}n\,\Gamma_n^2,\qquad\sum_{n\ge2}n\,\Gamma_n^2\ge0$ · **why** Split the sum in the drag formula of the Start; every term is n times a square, so the second part cannot be negative. · **plain** A fixed
     part, plus an extra that is never a saving.
  3. **did** Choose the minimum · **tex** $\Gamma_n=0\ (n\ge2)\ \Rightarrow\ \Gamma=\Gamma_1\sin\gamma$ · **why** A sum of non-negative terms is smallest when
     each is zero; the series $\Gamma=\sum\Gamma_n\sin(n\gamma)$ (14.18) then has one term. · **plain** The best loading is a single sine of the span
     angle.
  4. **did** Translate to the span coordinate · **tex** $\Gamma=\Gamma_1\sin(\gamma)=\Gamma_1\sqrt{1-(2y/s)^2}$ (14.23) · **why** From $y=-(s/2)\cos\gamma$,
     cos γ = −2y/s, and sin γ = √(1 − cos²γ) ≥ 0 on 0 ≤ γ ≤ π. · **plain** Plotted against y, the loading is half an
     ellipse.
  5. **did** Find its downwash · **tex** $w(y)=\dfrac{1}{2s}\,\Gamma_1\,\dfrac{\sin\gamma}{\sin\gamma}=\dfrac{\Gamma_1}{2s}$ (14.24) · **why** The series
     $w=\dfrac1{2s}\sum n\Gamma_n\dfrac{\sin(n\gamma)}{\sin\gamma}$ (14.19, last member) with only n = 1. · **plain** The air sinks at the same
     speed everywhere along the wing.
  6. **did** Write its drag · **tex** $D_i=\dfrac\pi8\,\rho\,\Gamma_1^2$ (14.25), first form · **why** Step 2 with the extra part zero. · **plain**
     The least induced drag, in terms of the amplitude.
  7. **did** Express it by the lift · **tex** $D_i=\dfrac\pi8\,\rho\,\dfrac{16L^2}{\pi^2s^2\rho^2U^2}=\dfrac{2L^2}{\pi\rho U^2s^2}$ (14.25), second form · **why** Square the Γ₁ of step 1 and insert it into step 6: 16/π² times π/8 leaves 2/π, and one ρ cancels. We want
     the drag in terms of the lift. · **plain** Induced drag ∝ lift squared, over speed squared and **span squared**.
  8. **did** Divide by ½ρU²A · **tex** $C_{D_i}=\dfrac{D_i}{\tfrac12\rho U^2A}=\dfrac{4L^2}{\pi\rho^2U^4s^2A}$ · **why** The definition of a wing's drag coefficient
     (R01, planform area A). · **plain** The same drag, as a coefficient.
  9. **did** Recognise the lift coefficient and the aspect ratio · **tex** $C_{D_i}=\dfrac{C_L^2}{\pi(s^2/A)}=\dfrac{C_L^2}{\pi\Lambda}$ (14.26) · **why**
     $C_L^2=4L^2/(\rho^2U^4A^2)$, so $4L^2/(\rho^2U^4)=C_L^2A^2$; one factor A cancels and s²/A is Λ (14.1). · **plain** The
     price of lift: lift coefficient squared over π times aspect ratio.
- **Result:** $\Gamma=\Gamma_1\sqrt{1-(2y/s)^2}$, $w=\Gamma_1/2s$, $D_i=\dfrac{2L^2}{\pi\rho U^2s^2}$, $C_{D_i}=\dfrac{C_L^2}{\pi\Lambda}$ — "an elliptic loading, a uniform
  downwash, and a drag that falls with the square of the span".
- **Check:** units — L²/(ρU²s²): N²/(N) = N ✓; $C_{D_i}$ dimensionless ✓. Limits — Λ → ∞: no induced drag (two-dimensional
  flow); L = 0: none either. Number — W = 10 kN, ρ = 1.25, U = 40 m/s, s = 10 m, A = 20 m²: $D_i$ = 318 N, $C_{D_i}$ = 0.0159,
  w = 1.27 m/s; and D_i = L × (w/U) = 10 000 × 0.0318 ✓. Direct quadrature of the downwash integral $w(y_1)=\frac1{4\pi}\int\frac{d\Gamma}{dy}\frac{dy}{(y_1-y)}$ (14.13) for the half-ellipse gives the same
  uniform w.
- **What it means:** induced drag is set by the *span loading* L/s, not by the area; doubling the span at fixed lift
  quarters it. In level flight it falls as 1/U². **Fails when:** the wing is not planar (winglets), the span is not
  what is held fixed (structural weight grows with span), or the sections stall.
- **Traps:** "minimum" is at fixed lift and span with a flat wake; Γ₁ is squared, so the factor is 16/π², not 4/π;
  dividing by ½ρU²A twice produces $C_L^2$ and leaves one factor A/s² = 1/Λ.

### D20 · The elliptic planform (14.27) and the lift slope of a finite wing — ★★, 8 steps, in C14 (notebook)
- **Goal:** find a wing shape that produces the elliptic loading without twist, and how much less lift such a wing gives
  per degree than its section alone. The second half is an exercise in the book. · **Start:**
  $dL=\rho U\Gamma\,dy=C_L\,\tfrac12\rho U^2c\,dy$ (14.27) — "the lift of a strip, written two ways". · **Plan:** • get Γ from $C_L$ and c • argue
  that $C_L$ is the same at every station • conclude the chord law • put it into the lifting-line equation • solve for Γ₁ •
  turn it into a lift slope. · **Tools:** proportionality; solving a linear equation; area of an ellipse (gloss). ·
  **Assumptions:** elliptic loading; the same section along the span and no twist (step 2); K constant.
- **Steps:**
  1. **did** Equate the two forms · **tex** $\Gamma=\tfrac12\,U\,c\,C_L$ · **why** Both sides of $\rho U\Gamma\,dy=C_L\tfrac12\rho U^2c\,dy$ (14.27) are the lift of the same strip; cancel the common
     factor ρU dy. · **plain** At each station the
     circulation is proportional to chord × lift coefficient.
  2. **did** Show $C_L$ is uniform · **tex** $C_L=K\Big[\alpha-\dfrac{\Gamma_1}{2sU}+\beta\Big]=\text{const}$ · **why** For the elliptic loading $w=\Gamma_1/2s$ (14.24) at every station; with one section (same K, β) and no twist
     (same α), the section law $C_L=K[\alpha-w/U+\beta]$ (14.16) gives the same value everywhere. · **plain** Every strip works at the same lift coefficient.
  3. **did** Conclude the chord law · **tex** $c(y)\propto\Gamma(y)\ \Rightarrow\ c=c_0\sqrt{1-(2y/s)^2}=c_0\sin\gamma$ · **why** Step 1 with $C_L$ and U
     constant: chord and circulation have the same shape, and Γ is elliptic. · **plain** An elliptic planform, root
     chord c₀.
  4. **did** Put it into the lifting-line equation · **tex** $\Gamma_1\sin\gamma=\dfrac K2\,U\,c_0\sin\gamma\,\Big[\alpha+\beta-\dfrac{\Gamma_1}{2sU}\Big]$ · **why**
     $\Gamma(y)=\dfrac K2Uc\Big[\alpha-\dfrac wU+\beta\Big]$ (14.17) with Γ = Γ₁ sin γ, c = c₀ sin γ, w = Γ₁/2s. · **plain** The same sin γ on both
     sides: one sine really is enough.
  5. **did** Cancel sin γ and solve for Γ₁ · **tex** $\Gamma_1\Big(1+\dfrac{Kc_0}{4s}\Big)=\dfrac K2\,U\,c_0\,(\alpha+\beta)$ · **why** Divide by sin γ (non-zero
     between the tips); move the Γ₁ term to the left. · **plain** A single linear equation for the single unknown.
  6. **did** Express c₀/s by the aspect ratio · **tex** $A=\dfrac\pi4\,c_0\,s\ \Rightarrow\ \Lambda=\dfrac{s^2}{A}=\dfrac{4s}{\pi c_0}\ \Rightarrow\ \dfrac{Kc_0}{4s}=\dfrac{K}{\pi\Lambda}$ · **why** An
     ellipse with axes s and c₀ has area π(s/2)(c₀/2). · **plain** The correction term depends only on K and Λ.
  7. **did** Form the wing's lift coefficient · **tex** $C_L=\dfrac{\tfrac{\pi s}{4}\rho U\Gamma_1}{\tfrac12\rho U^2A}=\dfrac{K(\alpha+\beta)}{1+K/(\pi\Lambda)}$ · **why** $L=\tfrac{\pi s}4\rho U\Gamma_1$ (14.21)
     with Γ₁ from step 5 and A from step 6: πsc₀/(4A) = 1. · **plain** The wing's lift coefficient is the section's,
     divided by 1 + K/(πΛ).
  8. **did** Read off the slope · **tex** $\dfrac{dC_L}{d\alpha}=\dfrac{K}{1+K/(\pi\Lambda)};\qquad K=2\pi:\ C_{L,3D}=\dfrac{C_{L,2D}}{1+2/\Lambda}$ · **why** Differentiate step 7
     with respect to α; with K = 2π the correction K/(πΛ) is 2/Λ. · **plain** A finite wing lifts less per degree than
     its section: by a quarter at Λ = 6.
- **Result:** an untwisted wing of elliptic planform has the elliptic loading, and $\dfrac{dC_L}{d\alpha}=\dfrac{K}{1+K/(\pi\Lambda)}$ (ours);
  $C_{L,3D}=\dfrac{C_{L,2D}}{1+2/\Lambda}$ for K = 2π.
- **Check:** units — all dimensionless ✓. Limits — Λ → ∞: the section's slope K; Λ → 0: slope → πΛ (very stubby wings lift
  little). Number — Λ = 6, K = 2π: 4.71 per radian (0.75 × 2π); Λ = 9.6: 5.20. Code — the solver with an elliptic planform
  returns one coefficient and this $C_L$ to round-off.
- **What it means:** the downwash reduces every strip's angle of attack by the same w/U = $C_L/(\pi\Lambda)$, so the whole lift
  line is tilted, not bent. **Fails when:** low aspect ratio or sweep (lifting-line theory itself fails); other
  planforms need the solver (their slope differs by a few per cent).
- **Traps:** elliptic *loading* is not elliptic *planform* (trap T11): the planform gives the loading only if untwisted
  with one section; the correction is K/(πΛ) in general and 2/Λ only for K = 2π.

### D21 · The best glide and the U-shaped drag curve — ★★, 8 steps, in C15 (notebook · `elliptic_loading_aspect_ratio`) — ours
- **Goal:** find the flight condition with the least drag for a given weight, and how the two kinds of drag trade places
  with speed. The book says only that the induced drag grows as α². · **Start:** $C_D=C_{D0}+\dfrac{C_L^2}{\pi\Lambda e_{\text{span}}}$ (ours), the
  profile drag plus the induced drag of $C_{D_i}=\dfrac{C_L^2}{\pi\Lambda}$ (14.26) corrected for a loading that is not elliptic. · **Plan:** •
  form drag over lift • minimise it • evaluate the optimum • use level flight to write the drag in newtons against speed. ·
  **Tools:** the minimum of a/x + bx (P349); quotient rule (P265, reminder); level flight $C_L=W/(\tfrac12\rho U^2A)$ (D01). ·
  **Assumptions:** a parabolic polar — $C_{D0}$ and $e_{\text{span}}$ constant; steady level flight (step 6).
- **Steps:**
  1. **did** Abbreviate the polar · **tex** $C_D=C_{D0}+k\,C_L^2,\qquad k\equiv\dfrac{1}{\pi\Lambda e_{\text{span}}}$ · **why** A shorthand; k is small for a
     slender, efficient wing. · **plain** Drag coefficient = a constant + a part that grows with lift squared.
  2. **did** Form drag over lift · **tex** $\dfrac DL=\dfrac{C_D}{C_L}=\dfrac{C_{D0}}{C_L}+k\,C_L$ · **why** Both forces share the factor ½ρU²A. In a
     glide, D/L is the glide slope; in level flight D = W × D/L. · **plain** One term falls with $C_L$, the other rises.
  3. **did** Set the derivative to zero · **tex** $\dfrac{d}{dC_L}\Big(\dfrac{C_D}{C_L}\Big)=-\dfrac{C_{D0}}{C_L^2}+k=0$ · **why** A smooth function has its
     minimum where its slope vanishes (P349); the second derivative $2C_{D0}/C_L^3$ is positive. · **plain** The best
     point is where the two slopes cancel.
  4. **did** Solve · **tex** $C_L^{*}=\sqrt{\dfrac{C_{D0}}{k}}=\sqrt{\pi\Lambda e_{\text{span}}\,C_{D0}}\quad\Rightarrow\quad k\,C_L^{*2}=C_{D0}$ · **why** Step 3 says $k=C_{D0}/C_L^2$; solve for $C_L$ and take the positive root (a lifting wing). Multiplying back
     gives the equality of the two drag parts. · **plain** At the best glide the induced drag equals the profile drag.
  5. **did** Evaluate the best ratio · **tex** $\Big(\dfrac LD\Big)_{\max}=\dfrac{C_L^{*}}{2C_{D0}}=\dfrac12\sqrt{\dfrac{\pi\Lambda e_{\text{span}}}{C_{D0}}}$ · **why** At the optimum
     $C_D=2C_{D0}$; insert $C_L^*$. · **plain** A slender, clean wing glides far: L/D grows like √Λ.
  6. **did** Use level flight · **tex** $C_L=\dfrac{W}{\tfrac12\rho U^2A}$ · **why** L = W (D01). Each speed corresponds to one lift
     coefficient. · **plain** Fast flight is low-$C_L$ flight.
  7. **did** Write the drag in newtons · **tex** $D=\tfrac12\rho U^2A\,C_{D0}+\tfrac12\rho U^2A\,k\,\dfrac{W^2}{(\tfrac12\rho U^2A)^2}$ · **why** D = ½ρU²A·$C_D$ with
     step 1 and step 6. · **plain** A part that grows with U² and a part that falls with U².
  8. **did** Simplify the induced part · **tex** $D=\tfrac12\rho U^2A\,C_{D0}+\dfrac{2W^2}{\pi\rho U^2s^2\,e_{\text{span}}}$ · **why** k/(½ρU²A) = 2/(πΛ$e_{\text{span}}$ρU²A)
     and ΛA = s² (14.1). The second term is $D_i=\dfrac{2L^2}{\pi\rho U^2s^2}$ (14.25) over $e_{\text{span}}$. · **plain** A U-shaped curve:
     profile drag wins at high speed, induced drag at low speed.
- **Result:** $\Big(\dfrac LD\Big)_{\max}=\tfrac12\sqrt{\dfrac{\pi\Lambda e_{\text{span}}}{C_{D0}}}$ at $C_L=\sqrt{\pi\Lambda e_{\text{span}}\,C_{D0}}$, where $C_{D0}=C_{D_i}$; in level flight
  $D=\tfrac12\rho U^2AC_{D0}+\dfrac{2W^2}{\pi\rho U^2s^2e_{\text{span}}}$ (ours).
- **Check:** units — both drag terms in N ✓. Limits — $C_{D0}$ → 0: L/D → ∞ at $C_L$ → 0; Λ → ∞: no induced drag, the best
  speed → the lowest possible. Number — Λ = 5, $e_{\text{span}}$ = 1, $C_{D0}$ = 0.02: $C_L^*$ = 0.56, L/D = 14; our aircraft: 16.5 at
  34.5 m/s, D = 546 N. Code — a brute-force scan of $C_L/C_D$ finds the same maximum.
- **What it means:** the induced share of the drag is exactly one half at the best glide, more when slower, less when
  faster — "induced drag is the largest part" is a statement about speed (trap T9). **Fails when:** $C_{D0}$ varies
  (it does, with Reynolds number and near stall), or compressibility adds wave drag (Ch. 15).
- **Traps:** minimise $C_D/C_L$, not $C_D$; the parabolic polar is a model; the minimum-drag speed is not the
  minimum-power speed (power is DU).

### D22 · Sailing: the apparent wind, the drive, and the closest course — ★★, 8 steps, in C16 (notebook · `thrust_from_lift`) — ours
- **Goal:** find the force that drives a sailing boat along its course and how close to the wind it can point. The book
  resolves the force in words; here it is in symbols. · **Start:** $\mathbf U_r=\mathbf U-\mathbf V$ — "the wind the boat feels is the true
  wind minus the boat's own velocity"; the sail's lift L is perpendicular to $\mathbf U_r$ and its drag D is along it. · **Plan:** •
  components of the apparent wind along and across the course • its angle off the bow • resolve lift and drag along and
  across the course • find where the drive vanishes • add the hull's balance. · **Tools:** relative velocity (P350);
  resolving a force; `np.arctan2` (P70, reminder); the tangent of a difference of two angles (gloss). · **Assumptions:**
  steady sailing; forces in the horizontal plane; quasi-steady sail and keel forces (trap T18).
- **Steps:**
  1. **did** Choose axes on the boat · **tex** $\hat{\mathbf t}\ \text{along the course},\quad\hat{\mathbf n}\ \text{across it, to leeward};\qquad\mathbf V=V\,\hat{\mathbf t}$ · **why** The boat can
     only use the force component along its course; the keel resists the component across it. · **plain** Forward
     and sideways, as the boat sees them.
  2. **did** Resolve the true wind · **tex** $\mathbf U=-U\cos\theta\,\hat{\mathbf t}+U\sin\theta\,\hat{\mathbf n}$ · **why** The wind comes *from* the direction at the
     angle θ off the bow on the windward side, so it blows backward (−$\hat{\mathbf t}$) and to leeward (+$\hat{\mathbf n}$). · **plain** A
     headwind part and a crosswind part.
  3. **did** Subtract the boat's velocity · **tex** $\mathbf U_r=-(U\cos\theta+V)\,\hat{\mathbf t}+U\sin\theta\,\hat{\mathbf n},\qquad\tan\beta_a=\dfrac{U\sin\theta}{U\cos\theta+V}$ · **why**
     $\mathbf U_r=\mathbf U-\mathbf V$ (P350); $\beta_a$ is the angle between the bow and the direction the apparent wind comes from. Since V > 0,
     $\beta_a<\theta$. · **plain** The felt wind is stronger from ahead, and closer to the bow than the true wind.
  4. **did** Write the sail's two forces · **tex** $\mathbf D=D\,(-\cos\beta_a\,\hat{\mathbf t}+\sin\beta_a\,\hat{\mathbf n}),\qquad\mathbf L=L\,(\sin\beta_a\,\hat{\mathbf t}+\cos\beta_a\,\hat{\mathbf n})$ · **why**
     Drag points along the apparent wind (the unit vector of $\mathbf U_r$); lift is that direction turned by 90° toward the
     bow, with $L=C_L\tfrac12\rho U_r^2A$ and likewise D. · **plain** Drag pulls back and to leeward; lift pulls forward and
     to leeward.
  5. **did** Add the components along the course · **tex** $F_{\text{drive}}=L\sin\beta_a-D\cos\beta_a$ · **why** The $\hat{\mathbf t}$-parts of step 4. ·
     **plain** The drive is the forward part of the lift minus the backward part of the drag.
  6. **did** Add the components across the course · **tex** $F_{\text{side}}=L\cos\beta_a+D\sin\beta_a$ · **why** The $\hat{\mathbf n}$-parts of step 4. ·
     **plain** A side force, usually several times the drive, pushing the boat to leeward.
  7. **did** Find where the drive vanishes · **tex** $F_{\text{drive}}=0\ \Longleftrightarrow\ \tan\beta_a=\dfrac DL\equiv\tan\varepsilon_A$ · **why** Set step 5 to zero and divide by
     L cos $\beta_a$. $\varepsilon_A$ is the sail's drag angle: the angle by which its total force leans back from the lift. ·
     **plain** No sail can drive a boat when the felt wind is closer to the bow than its own drag angle.
  8. **did** Balance against the hull · **tex** $\dfrac{F_{\text{drive}}}{F_{\text{side}}}=\tan(\beta_a-\varepsilon_A)=\dfrac{D_H}{L_H}=\tan\varepsilon_H\ \Rightarrow\ \beta_a=\varepsilon_A+\varepsilon_H$ · **why** In steady sailing the hull's drag $D_H$ balances the drive and the keel's lift $L_H$ the side force. Steps 5
     and 6 divided by L cos $\beta_a$ form the tangent of a difference of two angles. · **plain** In
     steady sailing the apparent wind sits off the bow by the sum of the two drag angles — sail's and hull's.
- **Result:** $F_{\text{drive}}=L\sin\beta_a-D\cos\beta_a$, zero at $\tan\beta_a=D/L$; steady sailing at $\beta_a=\varepsilon_A+\varepsilon_H$ with
  $\varepsilon=\tan^{-1}(D/L)$ for sail and for hull (ours).
- **Check:** units — forces in N ✓; angles dimensionless ✓. Limits — V = 0: $\beta_a=\theta$ (apparent = true wind); D = 0 for
  both foils: $\beta_a$ → 0, the boat could point straight into the apparent wind; θ = 180° (running): $\beta_a$ = 180°, drive
  = +D: the boat is pushed. Number — U = 4, V = 3 m/s, θ = 90°: $\beta_a$ = 53.1°; L = 100 N, D = 25 N: drive 65 N, side
  force 80 N. Our boat: sail L/D = 6, hull L/D = 8: 9.5° + 7.1° = 16.6°.
- **What it means:** sailing upwind is a contest between two lift-to-drag ratios — the sail's in air and the keel's in
  water; improve either and the boat points higher. The same triangle, with the stream and the motion exchanged, is the
  fish tail and the flapping wing. **Fails when:** the boat heels (forces leave the horizontal plane), in gusts
  (unsteady), or when the sail stalls.
- **Traps:** the angle that matters is to the **apparent** wind, which is further ahead and faster than the true wind; L
  is normal to $\mathbf U_r$, not to the sail's cloth; the side force does not vanish — the keel must carry it; the book's
  picture is quasi-steady (trap T18).

---

## Closing check — design gate (bullets only, so that the ledger parser of Part E is not disturbed)

Counted by a scratch script kept outside the repository (it parses this file by its headings, `nb.…` rows and table rows),
on the file as saved on 2026-10-12.

- **CORE blocks:** 16 of 16 (`#### C01 —` … `#### C16 —`, each with one `nb.core` row), every one with plain words, the idea, a
  worked example with easy numbers, a code row with *explain*, a from-scratch `nb.check_agree` row, at least one visual with
  see / read / change, and a "What would change if…" link. Blocks without an explainer cell of their own (C01, C02, C05, C08,
  C09, C11, C12, C14) each have at least one Python visual of their own (A.10).
- **Sections:** 9 of 9 `nb.section` rows (14.1–14.9); §14.1 opens the C01 story before `nb.core("C01", …)`, §14.8 opens the
  C16 story before `nb.core("C16", …)`; C07 (a §14.3 picture) is taught in §14.5 after C06, as the curation rules.
- **Curation ids:** 113 of 113 placed (16 C · 83 N, each opened in bold with its depth letter · 12 `nb.recap` rows · 2 S in one
  `nb.pointer` row).
- **Derivations:** 22 of 22 `### Dnn ·` blocks in Part F and 22 `nb.derivation` rows in Part A, each in the CORE block the
  curation names; step counts equal the curation's (5, 5, 15, 11, 7, 7, 4, 5, 10, 7, 4, 7, 8, 7, 11, 8, 6, 9, 9, 8, 8, 8 = 169
  steps); every step has did / tex / why / plain (a *why* of 6–39 tokens, formula tokens included); every block has Goal,
  Start, Plan, Tools, Assumptions, Result, Check, What it means, Fails when, Traps. No step uses "similarly", "clearly" or
  "after some algebra"; the words "it can be shown" occur once, in the Goal of D08, quoting what the book says before the
  block shows it.
- **★★★ checks:** the four blocks D03, D04, D15, D18 carry a `check_src` listing; all four were extracted from this file and
  executed (sympy, each about 1 s): every assertion passes and the printed lines are `2*b*t`, `1.299038105676658`, `3`,
  `0 pi/2`.
- **Explainers:** E1–E9 and B1, each with title, summary, CORE ids, reference explainer, meta, physics (JS ↔ Python), views,
  controls with the meaning of their extreme values, presets, status, readouts, depth features (≥ 2 beyond Explain, Code and
  Derivation — E5, which has no derivation, has linked views, presets, status and notes), Explain sections ending in
  *Reading the current setting* with one text per regime, Derivation tab (D03 D04 → E1 · D07 D09 → E2 · D02 → E3 · D10 → E4 ·
  none → E5 · D11 D12 D13 → E6 · D15 D18 → E7 · D19 D21 → E8 · D22 → E9 · D11 → B1 — exactly the curation's assignment), Code,
  walkthrough (7, 7, 6, 7, 6, 8, 8, 7, 7 steps; 5 for B1; none over 45 words), equations, 4 check questions each (3 for B1),
  parity rows (≥ 4 `py:` rows each for E1–E9 — 7, 7, 4, 6, 6, 8, 9, 7, 6 written out — plus invariants; 3 for B1) and a fit plan for 360×640 and 844×345. Every ★★★ derivation is in an explainer
  (D03, D04 in E1; D15, D18 in E7). Notebook-only derivations: D01, D05, D06, D08, D14, D16, D17, D20, as the curation rules.
- **Primers:** 15 new (P336–P350), each with an `nb.primer` row in Part A placed before its first use and a "primer (in Cxx)"
  row in Part E whose Concept is the primer term verbatim.
- **Ledger:** 201 rows, none without "Explained by"; every C / R id in that column exists in the curation; reminders of
  earlier primers carry the chapter and P number of `knowledge/primers.md` (ranges checked: ch01 P01–P61 … ch12 P280–P306).
- **Part C:** 98 function rows naming 108 callables (AF 27 rows, 29 callables · LL 35 rows, 40 callables · ch14 36 rows, 39 callables), plus the reused functions of C.0
  and the script helpers of C.4. Every `ch14.` / `AF.` / `LL.` name used in Parts A, B, D, E, F exists in the implementer's
  modules as found on 2026-10-12 (109 distinct names checked by import) **except the three requested wrappers**
  `ch14.surface_point`, `ch14.edge_speed`, `ch14.map_point`; `scripts/ch14_drawings.py` is requested too (C.4, C.5).
- **Equation mentions:** every sentence, table cell or step field of the header, Part A, Part B and Part F that carries a
  book number also carries an equation in TeX (scratch scanner at sentence level; not counted, by design: heading lines,
  `viz:equations` meta lists, the `ref '…'` label of an equation card that shows its TeX, code listings and comments inside
  them, parity-row code, and the "Book eq." column of Part C).
- **Slips and traps:** 10 slip rows in the header, each with the block that teaches it in the house form; 18 traps, each with
  its callout block.
- **Cross-check against the implementer's code (read-only, for information):** called with this design's inputs, the modules
  return the *expect* numbers of Part A and the parity values of Part B — section (Γ_K = 15.962 m²/s, exact $C_L$ = 1.0583,
  suction peak −2.237, loop area 1.0542, edge speed 0.9123 U, force −51.1 / 584.4 N/m), lifting line (rectangle Λ = 6: $C_L$ =
  0.39535, $e_{\text{span}}$ = 0.95394, slope 4.5304; the N = 3 matrix entry by entry; twist 0.02285 rad), polar (L/D 16.475 at 34.48
  m/s; shares 0.6358 and 0.134), sailing (8.394 m/s, 0.5299 rad, 374.6 N, 980.9 N; steady speed 3.951 m/s), fish (0.8062) and
  the lumped-vortex ratios (0.6463, 0.8156, 0.8792, 0.9388). No discrepancy was found.
- **Deviations from the curation (reported to the orchestrator):** (1) the notebook's worked numbers use the input set of
  header convention 11, written literally in the imports cell, because `ch14.illustrative_inputs()` returns a different
  (equally ours) set — C.5 gives the two ways to unify them; (2) E8's fourth mode is "light aircraft" instead of a
  wind-turbine blade section, which has no lift = weight balance (the turbine is named in its notes); (3) E4's third bar is
  $C_y/\cos\alpha$ (the module's `lift_from_cp(suction="ideal")`), not "loop area × cos α": in ideal flow the chordwise suction force
  makes the lift exceed the normal force, and only the first form equals 2Γ/(Uc) (1.0582 against 1.0502); (4) E9's status
  has a third regime between "sailing" and "in irons" (drive positive but less than the hull needs), because the drive
  vanishes at $\tan\beta_a=D/L$ while steady sailing needs $\beta_a\ge\varepsilon_A+\varepsilon_H$; (5) `sailboat_equilibrium` is the implementer's
  force-balance solver, not a closed form; (6) the Wagner ghost curve stays conditional on the verifier (the function exists).
- **Still a hypothesis in its row** (not computed by the designer): the stagnation-point offsets and near-edge speed ratios
  of C02 row 9; the Thwaites separation positions (C07 row 14, E4); the roll-up frames (A4); the bird-stroke mean thrust is
  computed (0.41 N) but its inputs are arbitrary; E-explainer status percentages marked "live".
