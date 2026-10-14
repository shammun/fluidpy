# Chapter 16 — Introduction to Biofluid Mechanics: lesson design
(from `analysis/ch16_curation.md` (A 17 · B 107 · C 21 = 145 rows; CORE 17 · NOTE 113 · RECAP 13 · SKIP 2; 34 derivations D01–D34
(★ 12 · ★★ 14 · ★★★ 8: D13, D19, D20, D21, D24, D27, D30, D33); 11 statements; E1–E10 + backup B1; primers **P374–P394**) and
`analysis/ch16.md` (§2 inventory with the LaTeX the analyst read from the page images, §2b derivations a-d01 … a-d45, §4
implementation rows I1–I20, §9 slips S1–S33 and traps T1–T27). 2026-10-14, lesson-designer.

**Part C comes first in this file** (right after the conventions) because the concept-implementer writes
`fluidpy/ch16_biofluid_mechanics.py` in parallel: it is the contract both sides keep.

**Equation provenance (honest statement).** Every equation placed here is the analyst's page-image transcription (analysis §2) as
carried by the curation. For this design I re-read nine rendered pages against it (PDF page = printed page + 27): p806
($Q(t)=\mathrm K\dfrac{dp}{dt}+\dfrac{p(t)}{R}$ (16.49), $Q=\dfrac{dy}{dx}+Py$ (16.50), $ye^{\int P\,dx}=A+\displaystyle\int Qe^{\int P\,dx}dx$ (16.51), $p(t)=\dfrac1{\mathrm K}e^{-t/R\mathrm K}\displaystyle\int_0^tQ(\tau)e^{\tau/R\mathrm K}d\tau+p_0e^{-t/R\mathrm K}$ (16.52)); p807
($\dfrac{\partial A}{\partial t}+\dfrac{\partial(Au)}{\partial x}=0$ (16.53), the printed $\rho A\left(\dfrac{\partial u}{\partial t}+u\dfrac{\partial u}{\partial x}\right)=-\dfrac{\partial((p-p_e)A)}{\partial x}$ (16.54) — slip S5 seen —, $p-p_e=\dfrac{Eh}{a_0}\left(1-\dfrac{a_0}{a}\right)=\dfrac{Eh}{a_0}\left[1-\left(\dfrac{A_0}{A}\right)^{1/2}\right]$ (16.55),
$\dfrac{\partial A}{\partial t}+A_0\dfrac{\partial u}{\partial x}=0$ (16.56), $\rho\dfrac{\partial u}{\partial t}=-\dfrac{\partial p}{\partial x}$ (16.57), $p-p_e=\dfrac{Eh}{2a_0A_0}(A-A_0)$ (16.58)); p811 (the printed $\dfrac{d^2U}{dr^2}+\dfrac1r\dfrac{dU}{dr}-\dfrac{i\omega\rho}{\mu}U=\dfrac G\mu$ (16.76) — slip S6 seen —
and (16.77)–(16.83)); p818 ((16.129)–(16.136)); p822 ((16.143)–(16.148), the printed K in the numerator of (16.143) — slip S15 seen);
p829 ((16.172)–(16.178)); p836 ((16.194)–(16.199)); p837 ((16.200)–(16.206), the printed $\partial/\partial x(\partial\mathrm P/\partial x)$ of (16.201) — slip S23 seen); p841
((16.215)–(16.220)). **No difference from the analyst's transcription was found on those pages.** The other numbered equations
are the analyst's transcription; builders render the page (`tools/render_pages.py ch16 --eq 16.M`) before setting any equation that
is not on the nine pages above.

**One correction to the analysis and the curation found while writing D27 (flag for the implementer).** The first axial
correction of Dean's series has a **plus** sign: $w_1=+\dfrac{r(1-r^2)(19-21r^2+9r^4-r^6)}{5760}\cos\alpha$ (ours), not the minus printed in analysis row #113
and curation N90. Evidence (sympy, the D27 `check_src` below, run for this design): from $\nabla_1^2w_1=\dfrac1{2r}\left(\dfrac{\partial\psi_0}{\partial\alpha}\dfrac{\partial w_0}{\partial r}-\dfrac{\partial\psi_0}{\partial r}\dfrac{\partial w_0}{\partial\alpha}\right)=-2f(r)\cos\alpha$ with
$f=r(1-r^2)^2(4-r^2)/144>0$ the right side is negative where $\cos\alpha>0$, so $w_1>0$ there (a negative Laplacian with zero
boundary values is a positive bump); and only the plus sign gives a flux that **falls**, by $\dfrac{1541}{50400}\left(\dfrac K{576}\right)^2$, as (16.178) says. With the
plus sign the axial maximum moves toward $\alpha=0$, the wall the centrifugal term $-w^2\cos\alpha$ of (16.167) pushes the fluid
toward (the outer wall). `ch16.dean_axial_correction` must return the plus sign (Part C).

Numbers marked *expect* were computed for this design by the designer's own scratch numpy / scipy / sympy scripts written from
the equations above. They are hypotheses: the builder compares an executed cell against them and **reports** any difference
instead of editing the expectation. The eight `check_src` listings of Part F (D13, D19, D20, D21, D24, D27, D30, D33) **were
executed as written for this design and passed** (sympy 1.x, 0.3–1.5 s each; the printed lines are quoted under them).

**Binding conventions for every builder (analysis §9, curation decisions 4–8, §8–§9).**
1. **Imports and aliases.** `from fluidpy import ch16_biofluid_mechanics as ch16`; reused: `from fluidpy.core import laminar as
   LAM, acoustics as AC, similarity as SIM, thermo as TH, units as UN`; `from fluidpy import ch01_introduction as ch01`.
   **`ch16` re-exports** `LAM.pipe_poiseuille`, `LAM.pipe_flow_rate`, `LAM.pipe_wall_stress`, `LAM.stokes_layer`,
   `LAM.stokes_second_problem`, `AC.dalembert`, `AC.acoustic_impedance`, `AC.velocity_from_pressure`, `AC.wave_equation_step`,
   `AC.pulse`, `AC.wave_speeds`, `SIM.reynolds_number`, `SIM.strouhal_number`, `ch01.capillary_rise`, `ch01.newton_shear_stress`,
   `TH.R_U`, `TH.G0` (the ch11–ch15 pattern), so explainer parity rows write `ch16.<name>` only.
2. **Scalar-callable and parity-friendly** (the `shot.py` evaluator knows `np`, `math`, `ch16`): every public closed form accepts
   Python floats and returns a float, a tuple or a `dict` of floats (arrays only when an array goes in).
3. **Units.** SI inside every function. mmHg, CGS (poise, dyn/cm²), bar and the PRU appear only through the converters
   (`ch16.mmHg_to_Pa`, `ch16.cP_to_Pa_s`, `ch16.dyn_cm2_to_Pa`, `ch16.bar_to_Pa`, `ch16.PRU`). Heart rate is passed as `HR` in
   beats per minute and turned into $\omega=2\pi\,HR/60$ inside. Hematocrit is a **fraction** in code (0.45), a percent in prose
   (trap T13). Pressures in the circulation are **transmural, gauge** (relative to the outside, p̄_V = 0 — trap T12).
4. **Two conventions on screen wherever a rival exists** (the "which $p_o$?" pattern; `ch16.conventions_table()` returns it; the
   conventions box before C01 shows it; each block repeats its row).

| # | What | The book (and this notebook) | Elsewhere | In code | Shown as |
|---|---|---|---|---|---|
| i | time convention | rigid tube $\dfrac{\partial p}{\partial x}=-Ge^{i\omega t}$ (16.73), $u(r,t)=U(r)e^{i\omega t}$ (16.75); elastic tube $u=\hat u(r)e^{i(kx-\omega t)}$ (16.116) | either | every Womersley docstring "$e^{+i\omega t}$", every elastic-tube docstring "$e^{i(kx-\omega t)}$"; `womersley_*` argument $i^{3/2}\alpha$, `elastic_tube_*` $\hat\beta=e^{i\pi/4}\alpha$ | two-row table in C09 and again at the top of C11 (trap T1) |
| ii | pressure gradient | $G=-dp/dx>0$ drives flow in +x; $\dfrac{dp}{dx}=-\dfrac{\Delta p}{L}$ (16.14) | signed $dp/dx$ (ch08 `dpdz`) | `G=` keyword as ch08 `LAM`; `dp` = $\Delta p=p_1-p_2>0$ | one line in R05 (trap T4) |
| iii | wall stress | $\tau=-\mu\,du/dr>0$ (a magnitude) | signed $\tau=\mu\,du/dR$ (ch08) | functions return the positive magnitude | one clause in R07 (trap T5) |
| iv | Reynolds and Strouhal | $Re=2au/\nu$, $St=2a\omega/u$ (diameter) ⇒ $\alpha^2=Re\,St/4$ | $\alpha^2=2\pi\,Re\,St$ (other Re, St); Ex. 16.2 St = aω/U | `womersley_from_re_st(Re, St)` uses the book's | "which Re?" table in C03 (trap T6) |
| v | Dean number | $\kappa=\sqrt\delta\,2Re$ (16.171); $K=2\kappa^2$ (16.178); $D=4\sqrt2\,\kappa$ (16.179), $K=D^2/16$ | De = Re√δ (diameter Re) = κ | `dean_number(..., convention=)` (required keyword), `dean_convert` | "which Dean number?" table in C14 and chips in E8 (trap T8) |
| vi | apparent viscosity | $\mu_{app}=\tau/\dot\gamma$ (the secant) | sometimes the tangent $d\tau/d\dot\gamma$ | `apparent_viscosity`, `differential_viscosity` | both drawn at one point (trap T27) |
| vii | Casson constants | $k_c$ dimensionless, $\sqrt{\tau_y/\mu_p}$ in s$^{-1/2}$ in (16.6), (16.7); $K_c=k_c\sqrt{\mu_p}$ in (Pa s)$^{1/2}$ in (16.207) | $\tau^{1/2}=a\dot\gamma^{1/2}+b^{1/2}$ | `whitmore_casson` returns both | two-row table in C02 (trap T14) |
| viii | Shapiro's P | $\mathrm P=(p-p_e)/K_p<0$ when collapsed, $-\mathrm P\approx\alpha^{-n}-1$ (16.188) | some plot −P | `collapsible_tube_law` returns P (negative in collapse) | ⚠️ callout in C15 (trap T9, slip S24) |
| ix | lapse rate | not used | — | — | the xylem block uses only $\rho gz$; no lapse rate appears in this chapter |

5. **Overloaded letters (trap T2) — the box shows this table; each block repeats its row.**

| Letter | Meanings in this chapter | Notebook / explainer symbol | Code name |
|---|---|---|---|
| α | Womersley number (16.12); Shapiro's area ratio $A/A_0$ (16.186); toroidal angle (Dean); $\alpha_c$ of Charm–Kurland (16.36) | α (Womersley), $\alpha$ "area ratio" said each time in C15, $\alpha$ "angle" in C14, $\alpha_c$ | `alpha`/`Wo`, `area_ratio`, `alpha_t`, `alpha_c` |
| K | compliance (16.48); metabolic cost per volume (16.138); Dean's K (16.178); $K_p$ bending stiffness; $K_c$, $k_c$ Casson | K (compliance, C07), $b$ "upkeep cost" in prose of C12 with "K in the book", $K$ Dean, $K_p$, $K_c$, $k_c$ | `K` / `compliance`, `b_metabolic`, `K_dean`, `Kp`, `Kc`, `k_c` |
| R | peripheral resistance (16.9); radius of curvature (Dean); reflection coefficient R (16.153); resistance per length in (16.204); gas constant (16.221) | R, $R_{curv}$, ℛ (reflection, written "R (reflection)" once), $R'$ per length, $R_u$ | `R`, `R_curv`, `refl`, `R_per_length`, `R_U` |
| δ | cell-free layer (16.21); curvature ratio (16.164) | δ (layer) in C05, $\delta=a/R$ in C14 | `delta`, `curvature` |
| ψ | secondary stream function (16.172); water potential (16.221) | ψ in C14, ψ in C17 (said "water potential") | `psi_sec`, `psi` |
| c | wave speed; concentration in (16.221) | c; $c_s$ "solute concentration" in C17 | `c`, `c_solute` |
| T | period; temperature; tensions $T_\theta$, $T_x$; transmission coefficient T (16.153) | $T$ period, $T_K$ temperature in C05/C17, $T_\theta$, 𝒯 transmission | `T_period`, `T`, `trans` |
| M, S, P | Shapiro's coefficient M (16.206), not a Mach number; speed index S (16.195); Shapiro's dimensionless pressure P (16.186) | M, S, P | `M_shapiro`, `S`, `P` |
| X, Y | real and imaginary parts of $1-F_2$ (16.97); our $X=k^2Eh/(2\rho a_0\omega^2)$ in (16.136) | X, Y in C10; $X_k$ in C11 | `X`, `Y`, `Xk` |
| κ, λ | Dean number; λ wavelength and the curved/straight resistance ratio (16.182) | κ, λ | `kappa`, `lam` |
| E, H, h, η, ν̂ | Young's modulus, Ê energy rate; hematocrit H (ch07: depth); wall thickness h; wall displacement η (ch07: elevation); Poisson ratio ν̂ vs kinematic ν | as listed | `E`, `H`, `h`, `eta`, `poisson` |

6. **Book slips taught in corrected form** — house form `> ⚠️ **slip Sk — the book prints** … **; the correct form is** …` where used,
   with the identity the printed form breaks; `ch16.book_slips()` (33 rows, each with `kind` and `changes_physics`).

| Slip | The book prints | Correct | Taught in | Planted wrong variant a test must fail |
|---|---|---|---|---|
| **S3** | $Q_c=\dfrac{\pi a^2\Delta p}{4\mu_pL}\left[a^2-\left(1-\dfrac{\mu_p}{2\mu_c}\right)\dfrac{(a-\delta)^4}{a^2}\right]$ (16.30) | $Q_c=\dfrac{\pi(a-\delta)^2\Delta p}{4\mu_pL}\left[a^2-\left(1-\dfrac{\mu_p}{2\mu_c}\right)(a-\delta)^2\right]$ (16.30, corrected) | C05 (N34, D06 step 7) | `two_layer_flow_rates(printed_16_30=True)`: $Q_p+Q_c\neq$ (16.31) |
| **S4** | $\left[1-\dfrac{a_0}{Eh}(p-p_e)\right]^{-4}dp=-\dfrac{8\mu}{\pi a_0}Q\,dx$ (16.44, printed) | $\left[1-\dfrac{a_0}{Eh}(p-p_e)\right]^{-4}dp=-\dfrac{8\mu}{\pi a_0^4}Q\,dx$ (16.44, corrected) | C06 (N39) | pint dimension check of the printed form fails |
| **S5** | $\rho A\left(\dfrac{\partial u}{\partial t}+u\dfrac{\partial u}{\partial x}\right)=-\dfrac{\partial((p-p_e)A)}{\partial x}$ (16.54, printed) | $\rho A\left(\dfrac{\partial u}{\partial t}+u\dfrac{\partial u}{\partial x}\right)=-A\dfrac{\partial p}{\partial x}$ (16.54, corrected) | C08 (N50, D11 step 2) | `tube1d_run(printed_16_54=True)`: different momentum budget at finite amplitude |
| **S6** | right side $+G/\mu$ of (16.76) | $\dfrac{d^2U}{dr^2}+\dfrac1r\dfrac{dU}{dr}-\dfrac{i\omega\rho}{\mu}U=-\dfrac G\mu$ (16.76, right side corrected) | C09 (N56, D13 step 4) | residual of the coded U is −G/μ |
| **S8** | "$\tau(t)\rvert_{r=a}=\dfrac{\partial u}{\partial r}\Big\rvert_{r=a}=\dfrac{G_0a}2+\dfrac a2\sum G_nF(\alpha_n)e^{in\omega t}$" (16.92, printed) | $\tau_w=-\mu\dfrac{\partial u}{\partial r}\Big\rvert_a=\dfrac{G_0a}2+\dfrac a2\sum_1^NG_nF_2(\alpha_n)e^{in\omega t}$ (16.92, corrected) | C10 (N62) | finite-difference cell |
| **S10** | $F_2(\alpha)\approx\dfrac{2}{i^{1/2}\alpha}\left(1+\dfrac1{2\alpha}\right)$ (16.96, printed) | $F_2(\alpha)\approx\dfrac{2}{i^{1/2}\alpha}\left(1+\dfrac{i^{3/2}}{2\alpha}\right)$ (16.96, corrected) | C10 (N64) | `womersley_F2_asymptote(printed_16_96=True)`: error ∝ α⁻² not α⁻³ |
| **S15** | $a=\left[\dfrac{16\mu}{\pi^2}K\right]^{1/6}Q^{1/3}$ (16.140, printed) | $a=\left[\dfrac{16\mu}{\pi^2K}\right]^{1/6}Q^{1/3}$ (16.140, corrected) | C12 (D22 step 6) | `murray_radius(printed_16_140=True)`: wrong units; 10× too large for our numbers |
| **S23** | $\dfrac{\alpha K_p}{\rho c^2}\dfrac{\partial}{\partial x}\left(\dfrac{\partial\mathrm P}{\partial x}\right)dx$ in (16.201) | $2\dfrac{dc}c=\left(1+\dfrac{\alpha\,\partial^2\mathrm P/\partial\alpha^2}{\partial\mathrm P/\partial\alpha}\right)\dfrac{d\alpha}\alpha+\dfrac{dK_p}{K_p}+\dfrac{\alpha K_p}{\rho c^2}\dfrac{\partial}{\partial x}\left(\dfrac{\partial\mathrm P}{\partial\alpha}\right)dx$ (16.201, corrected) | C15 (N102, D29 step 8) | `collapsible_differentials_sympy(printed_16_201=True)` fails when P depends on x |
| **S26** | the minimum pressure drop of the Casson example quoted in dyn/cm² (private) | $\Delta p_{\min}=2L\tau_y/a$ — the printed dyn/cm² figure is 1000× too small | C16 (N105, D31 step 5) | private test |
| **S27** | "1 cP = 0.01 P = 1 dyn s/cm² = 0.1 N s/m²" | $1\,\mathrm{cP}=10^{-2}\,\mathrm P=10^{-2}\,\mathrm{dyn\,s/cm^2}=10^{-3}\,\mathrm{Pa\,s}$ | C02 (N18) | `cP_to_Pa_s(1) == 1e-3` |
| **S28** | the xylem ν quoted with the xylem Re (private) | the quoted Re needs water's ν = 10⁻⁶ m²/s, 10× smaller | C17 (N111) | private test |

   The 22 "loose" rows get one clause where they occur (S1 R07 · S2 N31/D05 · S7 N58/D15 · S9 N63 · S11 C10 · S12 N66/D19 · S13
   N68/D18 · S14 N70/D19 · S16 N77/D25 · S17 N79 · S18 N87 · S19 N92 · S20 N93 · S21 R13 · S22 N102/D29 · S24 C15/D30 · S25
   N105/D31 · S29 N19 · S30 C17/D34 · S31 N27 · S32 N07 · S33 S01).
7. **Traps** (true as printed, easy to misread) — each a `> ⚠️ Common confusion:` callout where it first bites (`ch16.traps()`):
   T1 two time conventions (C09 N55, C11 N68, N72) · T2 overloaded letters (box; every block) · T3 units (C01, C02, C17) · T4
   G = −dp/dx (R05) · T5 τ a magnitude (R07) · T6 which Re, St (C03 D03) · T7 two roots of (16.136) (C11 D21) · T8 Dean conventions
   (C14 N87) · T9 P < 0 in collapse (C15 N99, D30) · T10 divide by 1 − S² last (C15 D30) · T11 the loop runs counterclockwise
   (C01 D01) · T12 transmural pressures, p̄_V = 0 (C07 D09) · T13 hematocrit fraction vs percent (C02 N17) · T14 k_c vs K_c (C02
   N22, D02) · T15 the ⅓ rule is an estimate (C01 N14, C07) · T16 tube-law singularity (C06 D08) · T17 radius vs diameter forms of
   Moens–Korteweg (C08 D12) · T18 "closed end" means R = +1 (C13 N81) · T19 the FL closure is ours (C05 N36) · T20 Fig. 16.14's axis
   runs wall–centre (C09 N59) · T21 a_c ≥ a means no flow (C16 D31) · T22 ψ < 0 under tension (C17 D34) · T23 the rigid-tube
   Womersley flow is not a wave (C09 R11) · T24 "ours" demos labelled (all) · T25 Bessel overflow at large α (C09 P375) · T26
   equations stated without derivation (each said where it is stated) · T27 secant vs tangent (C02 N21, D02).
8. **Colour code (curation §5; the same in figures, derivations, bars and explainers):** pressure = **orange** · velocity and
   flow rate = **blue** · wall / radius / area = **teal** · stress and shear = **rose** · energy, work and cost = **amber** · the
   Newtonian or Poiseuille reference (ghost) = **muted** grey dashed · viscous / Stokes-layer parts = **purple** (accent).
   Every ratio on screen says what it is a ratio of ("$|Q|/|Q_0|$ flow amplitude / Poiseuille's at the same G").
9. **Honest labels.** (a) **"ours"** opens every appearance of: the PV-loop work sign and the schematic loop; the series pressure
   bars; the equal-volume sphere; the vessel scales table; the FL closure (also "model"); the compliance of a Hookean segment; the
   periodic Windkessel state; the nonlinear 1-D tube run; the synthetic pulsatile gradient; the large-α wall-layer derivation; the
   Bessel series of the from-scratch cell; the symmetric Murray angle; the energy identity at a junction; ψ₀, w₁ and the flux
   coefficient's derivation; our Dean reliability criterion; the Starling-resistor model; the xylem column and its viscous drop;
   van 't Hoff's law and the capillary comparison of §16.4. (b) **"illustrative"** on every input of header convention 10 and every
   preset. (c) **"model"** in the title of every curve of C05 and E-B1 that uses our closure. (d) **"cited"** on the Pries curve
   (Secomb & Pries 2013) and the vessel sizes.
10. **Public-repo rule.** No book number is reproduced; the book's worked values stay in `tests/book_values_ch16.json`. The
    notebook's numbers use **this design's illustrative input set, written out literally in the imports cell as small dicts**
    (`ch16.illustrative_inputs()` returns the same dict of dicts and is displayed once beside them; if the module's set differs,
    the builder says "both are ours" and uses the dicts below):

```python
MMHG = 133.322387415                                    # Pa per mmHg (exact definition)
heart = dict(EDV=130e-6, ESV=55e-6, HR=66.0, BSA=1.9,   # m^3, m^3, beats/min, m^2 (illustrative adult, ours)
             pS=125*MMHG, pD=75*MMHG, systole_share=0.35,
             p_min=5*MMHG, p_ED=10*MMHG, p_open=80*MMHG, p_peak=125*MMHG, p_ES=100*MMHG)   # loop corners [Pa]
blood = dict(rho=1060.0, mu=3.5e-3, mu_p=1.2e-3, T=310.15)       # = ch16.BLOOD (cited typical values)
casson = dict(k_c=1.53, s0=2.0)                          # Whitmore's constants in (16.7): k_c [-], sqrt(tau_y/mu_p) [s^-1/2]
artery = dict(a=2e-3)                                    # a medium artery radius [m]
aorta = dict(a0=12e-3, h=1.5e-3, E=0.4e6, L=0.3)          # m, m, Pa, m
windkessel = dict(R=1.45e8, K=1.2e-8)                    # Pa s/m^3, m^3/Pa
micro = dict(a=20e-6, delta=2e-6, mu_c=4.0e-3)            # a 40-µm tube with a 2-µm plasma layer
arteriole = dict(a=25e-6, L=2e-3)
murray = dict(mu=3.5e-3, b_metabolic=1000.0, L=0.1)        # Pa s, W/m^3 (the book's K), m
bend = dict(delta=0.25)                                  # curvature ratio a/R_curv
vein = dict(a0=4e-3, Kp=200.0, n=1.5, L=0.1)              # m, Pa, -, m (a collapsible vein segment)
casson_pipe = dict(a=0.5e-3, L=0.5)                      # L/a = 1000
tree = dict(H=30.0, T=293.15, psi_soil=-0.1e6, c_leaf=300.0, a_xylem=25e-6, u_xylem=1e-3, mu_w=1.0e-3, nu_w=1.0e-6)
phloem = dict(a=10e-6, r_pore=1e-6, L=300e-6, ell=1e-6, N=50, u=0.3e-3, mu=1.5e-3)
```
    Derived values the notebook prints from these (expect): SV = 75 mL, EF = 57.7 %, CO = 4.95 L/min = 8.25 × 10⁻⁵ m³/s, CI = 2.61
    L/(min m²); pulse pressure 50 mmHg; ⅓-rule mean 91.7 mmHg; ω = 6.912 rad/s; ν = 3.302 × 10⁻⁶ m²/s; τ_y = 4.80 mPa, K_c =
    0.05300 (Pa s)^½, K_c² = 2.809 mPa s; RK = 1.74 s; c(aorta) = 4.856 m/s; A₀(vein) = 5.027 × 10⁻⁵ m².

---

## Part C — functions the builders will call (the implementer's contract)

Drafted from analysis §4 (rows I1–I20) and curation §8–§9, **before** the module existed (it was not on disk when this design was
written; the implementer works in parallel). Names, argument order, keyword-only markers (`*,`) and return keys below are the
contract: the notebook (Part A) and every explainer parity row (Part B) call exactly these. **(+)** = not in analysis §4 (listed
in C.5). Return shapes: scalars in → floats out; a `dict` has at least the keys listed. Outside its range a function returns NaN
or raises `ValueError` with the stated meaning. Units in brackets. Every docstring cites section + equation (written out), symbols
with SI units, assumptions, the time convention where one applies, and its validation label.

### C.0 Reused (existing; called, not changed; all re-exported by `ch16`)
| Function (as on disk) | Used for |
|---|---|
| `LAM.pipe_poiseuille(R, a, dpdz=None, mu=1e-3, G=None)` (R = radial coordinate, a = pipe radius) | R04: $u_z=\dfrac{R^2-a^2}{4\mu}\dfrac{dp}{dz}$ (8.6) |
| `LAM.pipe_flow_rate(a, dpdz=None, mu=1e-3, G=None)` → (Q, V, u_max) | R06 |
| `LAM.pipe_wall_stress(a, dpdz=None, G=None)` | R07: $\tau_0=\dfrac a2\dfrac{dp}{dz}$ (8.8) |
| `LAM.stokes_layer(nu, omega)` → dict `k`, `delta_e`, … | C03 D03, C09 D15: $\delta_e=\sqrt{2\nu/\omega}$ |
| `LAM.stokes_second_problem(y, t, U=1.0, omega=2π, nu=1e-6)` | C09 large-α check |
| `AC.dalembert(x, t, c, f=None, g=None)` | R09: $p'(x,t)=f(x-ct)+g(x+ct)$ (15.15) |
| `AC.acoustic_impedance(rho0, c)`, `AC.velocity_from_pressure(p_fn, x, t, rho0, *, t0=0.0, h=1e-4, …)` | R10, N78 |
| `AC.wave_equation_step(p_now, p_prev, C, bc="periodic")`, `AC.pulse(kind="gaussian", width=1.0, amplitude=1.0, x0=0.0)`, `AC.wave_speeds(U, c)` | A4, E3, N73 |
| `SIM.reynolds_number(U, l, nu=None, rho=None, mu=None)`, `SIM.strouhal_number(Omega, l, U)` | R03, C03, N111 |
| `ch01.capillary_rise(sigma, alpha, rho, R, g=G0)` (α = π/2 full wetting, the book's angle) | N110 |
| `ch01.newton_shear_stress(mu, du_dy)` | R01 |
| `TH.R_U` [J/(kmol K)] (÷1000 → J/(mol K)), `TH.G0` = 9.80665 | C17 |
| `UN.Q_`, `UN.check_dimensions(value, expected)` | P394, S4, S15 |
| `MCK.maccormack_step` (`fluidpy.core.maccormack`) | inside `tube1d_run` only |

### C.1 `fluidpy/ch16_biofluid_mechanics.py` — circulation bookkeeping, units (I1–I3)
| Function | Signature → return | Eq. |
|---|---|---|
| `MMHG`, `PRU` | constants: 133.322387415 Pa; `PRU = MMHG/1e-6` Pa s m⁻³ (mmHg per cm³/s) | N08, N12 |
| `mmHg_to_Pa`, `Pa_to_mmHg`, `cP_to_Pa_s`, `poise_to_Pa_s`, `dyn_cm2_to_Pa`, `bar_to_Pa`, `mL_to_m3` (+) | `(x)` → float/array | P394; $1\,\mathrm{cP}=10^{-3}\,\mathrm{Pa\,s}$ (S27 corrected) |
| `stroke_volume` | `(EDV [m³], ESV [m³])` → m³ | $SV=EDV-ESV$ (16.1) |
| `ejection_fraction` | `(SV, EDV)` → fraction (not percent) | $EF=(SV/EDV)\times100\%$ (16.2) |
| `cardiac_output` | `(SV [m³], HR [beats/min])` → m³/s | $CO=SV\times HR$ (16.3) |
| `cardiac_index` | `(CO [m³/s], BSA [m²])` → m³/(s m²) | $CI=CO/BSA$ (16.4) |
| `pulse_pressure` | `(pS, pD)` → Pa | $p_S-p_D$ |
| `peripheral_resistance` | `(pA_mean, Q, pV_mean=0.0)` → Pa s m⁻³ | $R=\Delta\bar p/Q$ (16.9), $\Delta p=\bar p_A-\bar p_V$ (16.10) |
| `mean_arterial_pressure` | `(pS, pD)` → Pa | $\bar p_A=\tfrac13p_S+\tfrac23p_D$ (16.11) |
| `time_mean` | `(p, t)` → trapezoid mean over the span of t | ours |
| `series_pressure_drop` | `(R_list, Q, p_in)` → array of the pressure after each element [Pa] | $\Delta p_i=QR_i$ (ours) |
| `network_resistance` | `(R_list, *, arrangement)` with `arrangement` = `"series"` or `"parallel"` (required keyword) → Pa s m⁻³ | ours |
| `mean_speed` | `(Q, A_total)` → m/s | $\bar u=Q/A_{tot}$ (ours) |
| `pv_loop_schematic` | `(EDV, ESV, *, p_min, p_ED, p_open, p_peak, p_ES, n=200)` → dict `V`, `p` (closed polygon, A→B→C→D, **counterclockwise in (V, p)**), `corners` (dict A, B, C, D of (V, p)). Legs (ours): filling A→B $p=p_{\min}+(p_{ED}-p_{\min})s^2$, $s=(V-ESV)/(EDV-ESV)$; isovolumic B→C at EDV; ejection C→D the parabola through $(EDV,p_{open})$, $(\tfrac12(EDV+ESV),p_{peak})$, $(ESV,p_{ES})$; isovolumic D→A at ESV. n points per curved leg | ours |
| `pv_loop_work` | `(p, V)` → dict `work` [J] $=-\oint p\,dV$ (trapezoid along the polygon), `signed_area` (shoelace in (V, p)), `orientation` (`"counterclockwise"` / `"clockwise"`) | $W=-\oint p\,dV$ (ours) |
| `pv_loop_schematic_work` (+) | `(EDV, ESV, *, p_min, p_ED, p_open, p_peak, p_ES)` → J, exact: $\Delta V\left[\tfrac16(p_{ES}+4p_{peak}+p_{open})-p_{\min}-\tfrac13(p_{ED}-p_{\min})\right]$ | ours |
| `cardiac_power` | `(W [J], HR [beats/min])` → W | ours |
| `BLOOD` | dict `rho`=1060.0, `mu`=3.5e-3 (high-shear), `mu_p`=1.2e-3 (plasma), `T`=310.15; sources in the module docstring (Wikipedia "Hemorheology", cited) | N15 |
| `equal_volume_sphere_area` | `(V)` → m², $(36\pi V^2)^{1/3}$ | ours |
| `VESSELS` (+), `vessel_scales_table` | `VESSELS` dict name → dict(a, h, u_mean) (cited typical sizes, not the book's); `vessel_scales_table(vessels=None, *, HR=66.0, rho=1060.0, mu=3.5e-3)` → DataFrame columns `vessel`, `a`, `h`, `u_mean`, `Re`, `alpha`, `h_over_a` | ours |
| `_check_hematocrit` | `(H)` → raises if H > 1 (percent passed) or ≤ 0 | trap T13 |

### C.2 Rheology and steady pipe flow (I4–I7, I18)
| Function | Signature → return | Eq. |
|---|---|---|
| `newtonian_stress` | `(gamma_dot, mu)` → Pa | $\tau=\mu\dot\gamma$ |
| `yield_stress_correlation` | `(H, C_F)` → Pa (H fraction > 0.1, C_F in g/100 mL; converts dyn/cm²) | $\sqrt{\tau_y}=(H-0.1)(C_F+0.5)$ (16.5) |
| `casson_stress` | `(gamma_dot, tau_y, Kc)` → Pa | $\tau^{1/2}=\tau_y^{1/2}+K_c\dot\gamma^{1/2}$ (16.207) |
| `casson_shear_rate` | `(tau, tau_y, Kc)` → 1/s, 0 for τ < τ_y | (16.207), yield branch |
| `apparent_viscosity` | `(tau, gamma_dot)` → Pa s (secant τ/γ̇) | $\mu_{app}=\tau/\dot\gamma$ |
| `casson_apparent_viscosity` | `(gamma_dot, tau_y, Kc)` → Pa s, $(\sqrt{\tau_y/\dot\gamma}+K_c)^2$ | D02 (ours) |
| `differential_viscosity` | `(gamma_dot, tau_y, Kc)` → Pa s, $d\tau/d\dot\gamma=K_c(\sqrt{\tau_y}+K_c\sqrt{\dot\gamma})/\sqrt{\dot\gamma}$ | trap T27 (ours) |
| `whitmore_casson` | `(gamma_dot, mu_p, *, k_c=1.53, s0=2.0)` → dict `tau`, `mu_app`, `tau_y` (= s0² μ_p), `Kc` (= k_c√μ_p) | $\sqrt{\tau/\mu_p}=1.53\sqrt{\dot\gamma}+2.0$ (16.7) |
| `poiseuille_velocity` | `(r, a, dp, L, mu)` → m/s (via `LAM.pipe_poiseuille`, G = dp/L) | $u(r)=\left(\dfrac{a^2\Delta p}{4\mu L}\right)\left(1-\dfrac{r^2}{a^2}\right)$ (16.15) |
| `poiseuille_flow_rate` | `(a, dp, L, mu)` → m³/s (float) | $Q=\dfrac{\pi a^4}{8\mu}\dfrac{\Delta p}L$ (16.17) |
| `poiseuille_wall_stress` | `(a, dp, L)` → Pa (positive) | $\tau=+\dfrac a2\dfrac{\Delta p}L$ (16.19, sign corrected — S1) |
| `vascular_resistance` | `(mu, L, a)` → Pa s m⁻³ | $R=\dfrac{8\mu L}{\pi a^4}$ (16.20) |
| `entrance_length` | `(Re, d, *, regime)` with `regime` = `"laminar"` / `"turbulent"` (required) → m; laminar raises for Re ≤ 50 unless `allow_low_re=True` | $\dfrac\ell d=0.06\,Re$, $\dfrac\ell d=0.693\,Re^{1/4}$ (16.37) |
| `casson_wall_stress` | `(dp, L, a)` → Pa | $\tau_w=\dfrac{\Delta p\,a}{2L}$ |
| `casson_min_pressure_drop` | `(L, a, tau_y)` → Pa | $2L\tau_y/a$ (ours, from (16.210) with $a_c=a$) |
| `casson_plug_radius` | `(dp, L, tau_y)` → m (may exceed a: no flow) | $a_c=\dfrac{2L\tau_y}{p_1-p_2}$ (16.210) |
| `casson_pipe_velocity` | `(r, a, dp, L, tau_y, Kc)` → m/s (plug value for r < a_c; zeros if a_c ≥ a) | (16.215)/(16.216) |
| `casson_core_velocity` | `(a, dp, L, tau_y, Kc)` → m/s | $u_c=\dfrac1{4K_c^2}\dfrac{\Delta p}L(\sqrt a-\sqrt{a_c})^3\left(\sqrt a+\tfrac13\sqrt{a_c}\right)$ (16.218) |
| `casson_pipe_flow_rate` | `(a, dp, L, tau_y, Kc)` → m³/s (0 if a_c ≥ a) | (16.220) |
| `casson_flux_ratio` (+) | `(xi)` → $1-\tfrac{16}7\xi^{1/2}+\tfrac43\xi-\tfrac1{21}\xi^4$ (0 for ξ ≥ 1) | bracket of (16.220) |
| `power_law_pipe_flow_rate`, `herschel_bulkley_pipe_flow_rate` (optional) | `(a, dp, L, m, n)`, `(a, dp, L, m, n, tau0)` | ours (S01 pointer) |

### C.3 The Fåhræus–Lindqvist two-layer model (I5)
| Function | Signature → return | Eq. |
|---|---|---|
| `cell_conservation` | `(Qc, Qp, H_c, a, delta)` → dict `H_F`, `H_T`, `Hc_over_HF`, `HT_over_Hc` | $Q_FH_F=Q_cH_c,\ Q_c+Q_p=Q_F,\ H_Ta^2=H_c(a-\delta)^2$ (16.8) |
| `two_layer_velocity` | `(r, a, delta, dp, L, mu_p, mu_c)` → m/s | (16.27), (16.28) |
| `two_layer_flow_rates` | `(a, delta, dp, L, mu_p, mu_c, *, printed_16_30=False)` → dict `Qp`, `Qc`, `Q` | (16.29), (16.30, corrected), (16.31) |
| `apparent_viscosity_two_layer` | `(a, delta, mu_p, mu_c)` → Pa s | (16.33) |
| `apparent_viscosity_thin_layer` | `(a, delta, mu_p, mu_c)` → Pa s | (16.34) |
| `charm_kurland_core_viscosity` | `(H_c, T, mu_p)` → Pa s (raises for H_c > 0.6) | (16.35), (16.36) |
| `fahraeus_lindqvist_model` | `(a, delta, H_F, mu_p, T)` → dict `H_c`, `mu_c`, `mu_app`, `H_T`, `residual` (model, ours; brentq) | N36 |
| `pries_relative_viscosity` | `(D_um, H_D)` → μ_app/μ_p (cited, Secomb & Pries 2013 Eqs. (3)–(5)) | N29 |

### C.4 Elastic tubes, Windkessel, pulse waves (I8–I11)
| Function | Signature → return | Eq. |
|---|---|---|
| `hoop_stress` | `(p, pe, a, h)` → Pa | $\sigma_{\theta\theta}=\dfrac{(p-p_e)a}{h}$ (16.39) |
| `tube_radius_hookean` | `(p, pe, a0, E, h)` → m; `ValueError` when $a_0(p-p_e)/Eh\ge1$ | $a=a_0\left[1-\dfrac{a_0}{Eh}(p-p_e)\right]^{-1}$ (16.41) |
| `tube_law_hookean` | `(A, A0, E, h, a0)` → p − pe [Pa] | (16.55) |
| `tube_law_linear` | `(A, A0, E, h, a0)` → p − pe [Pa] | (16.58) |
| `elastic_tube_flow_rate` | `(p1, p2, pe, a0, E, h, mu, L)` → m³/s | (16.46) |
| `elastic_tube_pressure` | `(x, Q, p1, pe, a0, E, h, mu)` → Pa (integrated (16.44, corrected) from x = 0) | (16.45, corrected) |
| `elastic_tube_resistance` | `(p1, p2, pe, a0, E, h, mu, L)` → (p1 − p2)/Q [Pa s m⁻³] (ours) | N40 |
| `compliance_from_tube` | `(p, pe, a0, E, h, L)` → m³/Pa, $dV/dp=2\pi a^3L/(Eh)$ with a from (16.41) (ours) | N45 |
| `ejection_half_sine` | `(t, SV, T_period, Ts)` → m³/s, periodic; $Q_{\max}=\pi SV/(2T_s)$ (ours) | C07 |
| `windkessel_rhs` | `(t, p, Q_fn, R, K)` → dp/dt | (16.49) |
| `windkessel_ode` | `(t, Q_fn, R, K, p0)` → p array (solve_ivp, rtol 1e-10) | (16.49) |
| `windkessel_pressure` | `(t, Q_fn, R, K, p0)` → p (quadrature form) | (16.52) |
| `windkessel_periodic_p0` | `(*, SV, HR, systole_share, R, K)` → p₀ [Pa] (closed form, half-sine ejection, ours) | D10 step 8 |
| `windkessel_periodic_state` | `(t, *, SV, HR, systole_share, R, K)` → p(t) [Pa] (closed form, ours) | D10 |
| `windkessel_summary` (+) | `(*, SV, HR, systole_share, R, K, n=2001)` → dict `p0`, `pS`, `pD`, `p_mean`, `map_third`, `tau`, `Q_mean`, `t_peak` | ours |
| `windkessel_terms` | `(t, *, SV, HR, systole_share, R, K)` → dict `Q`, `storage` (= K dp/dt), `runoff` (= p/R) [m³/s] | (16.49) split |
| `moens_korteweg_speed` | `(E, h, rho, a0, *, midwall=False)` → m/s | $c=\sqrt{\dfrac{Eh}{2\rho a_0}}$ (16.62); $c=\sqrt{\dfrac{Eh}{2\rho(a_0+h/2)}}$ (16.63) |
| `bramwell_hill_speed` | `(A, dpdA, rho)` → m/s | $c^2=\dfrac A\rho\dfrac{dp}{dA}$ (16.62, second form; (16.185)) |
| `youngs_modulus_from_pwv` | `(c, h, rho, a0)` → Pa | inverse of (16.62) |
| `linearisation_ratios` | `(p_amp, p_mean, u, c, *, Eh_over_a0=None)` → dict `p_ratio`, `u_over_c`, `strain` (= p_amp/(Eh/a0) when given) | N52 |
| `pulse_wave` | `(x, t, c, f1, f2=None)` → p, time shapes $f(t-x/c)$ | $p=f_1\left(t-\dfrac xc\right)+f_2\left(t+\dfrac xc\right)$ (16.65) |
| `velocity_from_pulse` | `(p, rho, c, *, direction=+1)` → m/s | $\hat p_1=\rho c\hat u_1$, $\hat p_2=-\rho c\hat u_2$ (16.68) |
| `wall_bulge` | `(p, *, a0, E, h)` → η = a − a₀ [m] from (16.41) | ours |
| `pulse_state` | `(x, t, *, c, amplitude, width, rho, a0, E, h, ic="right")` with `ic` = `"right"` (one right-going Gaussian) or `"rest"` (a bump released at rest: two halves) → dict `p`, `u`, `a` | (16.65), (16.68), (16.41) |
| `tube1d_run` | `(x, A_init, u_init, A0, E, h, a0, rho, t_end, *, inflow=None, CFL=0.8, nu_art=0.02, printed_16_54=False)` → dict `t`, `A`, `u`, `p` snapshots (ours, MacCormack + small artificial viscosity — DEVIATION box) | (16.53), (16.54, corrected), (16.55) |

### C.5 Womersley flow in a rigid tube (I12) — time convention $e^{+i\omega t}$
| Function | Signature → return | Eq. |
|---|---|---|
| `womersley_number` | `(a, omega, nu)` → α | $\alpha=a\sqrt{\omega/\nu}$ (16.12) |
| `stokes_number`, `womersley_from_re_st` | `(a, omega, nu)` → α²; `(Re, St)` → $\sqrt{Re\,St/4}$ | D03 (ours) |
| `womersley_regime` (+) | `(alpha, *, low=1.0, high=10.0)` → `"quasi-steady"` / `"transitional"` / `"inertial"` | E4 status (ours) |
| `bessel_j0_series`, `bessel_j1_series` (+) | `(z, nterms=60)` → complex; $\sum(-1)^m(z/2)^{2m}/(m!)^2$ and $\sum(-1)^m(z/2)^{2m+1}/(m!(m+1)!)$ | P374, from scratch |
| `womersley_F1` | `(r, a, alpha)` → complex $J_0(i^{3/2}\alpha r/a)/J_0(i^{3/2}\alpha)$ (scaled `jve` ratio beyond α ≈ 30) | (16.81) |
| `womersley_F2` | `(alpha)` → complex | $F_2(\alpha)=\dfrac{2J_1(i^{3/2}\alpha)}{i^{3/2}\alpha J_0(i^{3/2}\alpha)}$ (16.91) |
| `womersley_F2_series` | `(alpha, order=2)` → complex, $1-i\alpha^2/8$ (order 2) | (16.93) |
| `womersley_F2_asymptote` | `(alpha, *, printed_16_96=False)` → complex | (16.96, corrected) |
| `womersley_velocity` | `(r, t, a, omega, rho, mu, G)` → m/s, real part; forcing $-\partial p/\partial x=G\cos\omega t$ | (16.82) |
| `womersley_velocity_nd` (+) | `(r_over_a, wt, alpha)` → $u/(Ga^2/4\mu)$, real part | (16.82) scaled |
| `womersley_harmonics` | `(r, t, a, omega, rho, mu, G_n)` → m/s; `G_n[0]` mean, `G_n[n]` complex amplitude of harmonic n | (16.84)–(16.86) |
| `womersley_small_alpha`, `womersley_large_alpha` | `(r, t, a, omega, rho, mu, G)` → m/s | (16.87), (16.88) with $\sqrt{n\omega/2\nu}$ (S7) |
| `womersley_ode_residual` | `(r, a, omega, rho, mu, G, *, h=1e-6)` → complex value of the left side of (16.76) for the coded U (should be −G/μ) | (16.76, corrected) |
| `womersley_pde_solve` | `(a, omega, rho, mu, G, *, nr=200, steps_per_period=400, periods=8)` → dict `r`, `t`, `u` (last period), Crank–Nicolson in r, ours | (16.72) |
| `womersley_flow_rate` | `(t, a, omega, rho, mu, G)` → m³/s, real part, single mode | (16.90) |
| `womersley_flow_rate_harmonics` (+) | `(t, a, omega, rho, mu, G_n)` → m³/s | (16.89) (sum to N, S7) |
| `womersley_wall_shear` | `(t, a, omega, rho, mu, G_n)` → Pa; a scalar G is read as `[0, G]` | (16.92, corrected) |
| `womersley_XY` | `(alpha)` → (X, Y) of $1-F_2=X+iY$ | (16.97) |
| `womersley_flow_ratio` | `(alpha)` → $|Q|/|Q_0|=8\sqrt{X^2+Y^2}/\alpha^2$ | (16.99)/(16.95) (ours ratio) |
| `womersley_phase_lag` | `(alpha)` → φ [rad] = arctan2(X, Y) (lag of Q behind $G\cos\omega t$) | $\tan\phi=X/Y$ (16.100) |

### C.6 Viscous waves in an elastic tube (I13) — time convention $e^{i(kx-\omega t)}$
| Function | Signature → return | Eq. |
|---|---|---|
| `elastic_tube_characteristic` | `(Xk, alpha, poisson)` → complex left side of (16.136), $\hat\beta=e^{i\pi/4}\alpha$ | (16.136) |
| `elastic_tube_roots` | `(alpha, poisson=0.5)` → dict `pressure` (root → 1 as α → ∞, chosen by continuation), `other` | (16.136), trap T7 |
| `elastic_tube_wave` | `(alpha, c0=1.0, poisson=0.5, *, omega=1.0, branch="pressure", E_ratio=1+0j)` → dict `X`, `c`, `c_over_c0`, `k1`, `k2`, `k2_over_k1`, `damping_per_wavelength` (= $e^{-2\pi k_2/k_1}$), `branch`; $k/k_0=\sqrt X$ with Re > 0, $k_0=\omega/c_0$; `E_ratio` = $E^*/E$ | (16.136), (16.137) |
| `elastic_tube_determinant_numeric` | `(Xk, alpha, poisson)` → complex (ODE route, no Bessel closed form, ours) | (16.122)–(16.125), (16.129) |
| `_elastic_tube_modes` | `(r, A, eta, k, omega, nu, E, h, a0, rho, poisson)` → (û, v̂) | (16.131), (16.133) |

### C.7 Branching: Murray's law and junctions (I14, I15)
| Function | Signature → return | Eq. |
|---|---|---|
| `murray_cost` | `(a, Q, L, mu, b_metabolic)` → W | $\hat E_1=\dfrac{8\mu L}{\pi a^4}Q^2+K\pi a^2L$ (16.138) |
| `murray_terms` | `(a, Q, L, mu, b_metabolic)` → dict `pumping`, `upkeep`, `total` [W] | (16.138) split |
| `murray_radius` | `(Q, mu, b_metabolic, *, printed_16_140=False)` → m | $a=\left[\dfrac{16\mu}{\pi^2K}\right]^{1/6}Q^{1/3}$ (16.140, corrected) |
| `murray_min_cost` | `(Q, L, mu, b_metabolic)` → W | $\hat E_{1,\min}=\dfrac{3\pi}2KLa^2$ (16.141) |
| `murray_daughters` | `(a0, flow_split)` → (a1, a2), split = Q₁/Q₀ | $a_0^3=a_1^3+a_2^3$ (16.144) |
| `murray_angles` | `(a0, a1, a2)` → dict `theta`, `phi`, `total` [rad], `valid` (NaN + False if \|cos\| > 1) | (16.146) |
| `murray_theta` | `(a0, a1)` → rad | (16.147) |
| `murray_symmetric_angle` | `()` → arccos 2^{−1/3} = 0.65393 rad (37.467°) | ours |
| `bifurcation_cost` | `(B, A, C, D, a0, a1, a2)` → $a_0^2L_0+a_1^2L_1+a_2^2L_2$ [m³] | (16.145) |
| `bifurcation_point` | `(A, C, D, a0, a1, a2)` → dict `B`, `cost`, `theta`, `phi` (Nelder–Mead, ours) | D24 check |
| `tube_admittance`, `tube_impedance` | `(A, rho, c)` → m⁴ s/kg, kg/(m⁴ s) | $Y=A/\rho c$, $Z=\rho c/A$ (16.151) |
| `junction_coefficients` | `(Y_parent, Y_daughters)` → (R, T) | (16.153) |
| `junction_energy` | `(Y_parent, Y_daughters)` → dict `reflected` (= R²), `transmitted` (= ΣY_j T²/Y₁), `total` | (16.162), (16.163), identity (ours) |
| `junction_waveforms` | `(x, t, f, c1, R, *, Y1=1.0, P_I=1.0)` → dict `p`, `Q` in the parent (x < 0), time shape f | (16.159), (16.160) |
| `junction_from_geometry` | `(a_parent, a_daughters, c_parent, c_daughters, rho)` → dict `Y1`, `Y`, `ratio` (ΣY/Y₁), `R`, `T`, `reflected`, `transmitted` | ours |
| `network_pulse_run` | `(L_parent, L_daughters, c_parent, c_daughters, Y_parent, Y_daughters, *, width, nx=400, t_end)` → dict `x`, `t`, `p_parent`, `p_daughters`, `probe_p`, `probe_Q` (leapfrog at C = 1, ours) | A4 |

### C.8 Curved tubes: Dean (I16) — α here is the toroidal angle, α = 0 on the outer wall
| Function | Signature → return | Eq. |
|---|---|---|
| `curvature_ratio` | `(a, R_curv)` → δ | $\delta=a/R$ (16.164) |
| `dean_number` | `(delta, U_mean, a, nu, *, convention)` (`"kappa"`, `"K"`, `"D"`; required) | $\kappa=\sqrt\delta\,\dfrac{2a\bar W_0}\nu$ (16.171) |
| `dean_convert` | `(value, *, from_, to)` | $K=2\kappa^2$, $D=4\sqrt2\kappa$, $K=D^2/16$ |
| `dean_flux_ratio` | `(K=None, *, D=None, exact=False)` → $Q/(\pi a^2\bar W)$; `exact=True` uses 1541/50400 for 0.0306 | (16.178), (16.181) |
| `curved_resistance_ratio` | `(K)` → λ | (16.182) |
| `dean_secondary_streamfunction` | `(r, alpha_t)` → $\psi_0=\dfrac r{144}(1-r^2)^2(4-r^2)\sin\alpha$ | ours (D27) |
| `dean_axial_correction` | `(r, alpha_t)` → $w_1=+\dfrac{r(1-r^2)(19-21r^2+9r^4-r^6)}{5760}\cos\alpha$ (**plus sign** — header correction) | ours (D27) |
| `dean_axial_velocity` | `(r, alpha_t, K)` → $w_0+\kappa^2w_1$ with $w_0=2(1-r^2)$, $\kappa^2=K/2$ (scaled by $\bar W$) | ours |
| `dean_secondary_velocity` (+) | `(r, alpha_t, K)` → (u, v) $=\kappa\left(\dfrac1r\dfrac{\partial\psi_0}{\partial\alpha},-\dfrac{\partial\psi_0}{\partial r}\right)$ | (16.172) |
| `dean_series_sympy` | `(order=2)` → dict `w0`, `psi0`, `w1`, `flux_coefficient` (= Rational(1541, 50400)); cached | D27 |
| `dean_series_reliability` | `(K, *, threshold)` → dict `ratio` (K⁴-term / K²-term of (16.178)), `reliable` (ratio < threshold) — our criterion, never the book's range | ours |

### C.9 Collapsible tubes (I17) — α here is the area ratio $A/A_0$
| Function | Signature → return | Eq. |
|---|---|---|
| `bending_stiffness_parameter` | `(E, h, a0, poisson, *, coeff)` → Pa (coeff required) | $K_p\propto EI$ (16.187) |
| `collapsible_tube_law` | `(area_ratio, *, n=1.5, form="shapiro")` → P (negative when collapsed); `form="pure"` gives $-\alpha^{-n}$ | $-\mathrm P\approx\alpha^{-n}-1$ (16.188) |
| `collapsible_wave_speed` | `(area_ratio, Kp, rho, n=1.5)` → m/s | $c^2=\dfrac{nK_p\alpha^{-n}}\rho$ (16.189) |
| `speed_index` | `(u, c)` → S | $S=u/c$ (16.195) |
| `shapiro_M` | `(area_ratio, n=1.5)` → M = 3 + αP″/P′ (= 2 − n for the power law) | (16.206) |
| `choking_flux` | `(area_ratio, A0, Kp, rho, n=1.5)` → $Q^*=\alpha A_0c$ [m³/s] | S = 1 |
| `collapsible_terms` | `(area_ratio, Q, prm, *, dA0dx=0.0, dpedx=0.0, dzdx=0.0, dKpdx=0.0)` → dict `area` ($S^2A_0'/A_0$), `external`, `gravity`, `friction`, `stiffness` (each $-(\cdot)/\rho c^2$), `driver` (their sum), `one_minus_S2`, `dalpha_dx` (= α·driver/(1 − S²); NaN at S = 1) [1/m] | (16.204) |
| `collapsible_rhs` | `(x, y, prm, Q)` → dα/dx (y = [α]) | (16.204) |
| `collapsible_march` | `(Q, alpha_in, x_span, prm)` → dict `x`, `alpha`, `S`, `p`, `choked`, `x_choke` (terminal event at S² = 1 − 1e-9) | (16.204) |
| `starling_state` | `(p1, p2, pe, prm)` → dict `Q`, `alpha_out`, `S_out`, `limited` | ours |
| `flow_limitation_curve` | `(p2_list, *, p1, pe, prm)` → dict `p2`, `Q`, `S_out`, `limited`, `Q_choke`, `p2_choke` | ours |
| `collapsible_differentials_sympy` | `(printed_16_201=False)` → dict of sympy residuals for (16.200)–(16.205) | D29, D30 |
| `wall_friction_term` | `(u, area_ratio, d0, rho, mu, *, regime, f=None)` → τ_w s dx/A per dx | (16.192), (16.193) |
| **`prm` dict** | `dict(L, A0, Kp, n=1.5, rho, mu, friction=True, A0_fn=None, dzdx=0.0, pe_slope=0.0, Kp_slope=0.0)`; laminar resistance per length $R'=8\pi\mu/(\alpha A_0)^2$ so that $\tau_ws\,dx/A=R'Q\,dx$; `A0_fn(x)` (derivative by central difference) gives a throat (E9 throat mode; default uniform A₀); `dzdx = −1` is the giraffe jugular; `friction=False` gives the frictionless throat whose critical flux puts S = 1 exactly where $dA_0/dx=0$. **Expected** with the vein of header convention 10 (L = 0.1 m, A₀ = π(4 mm)², K_p = 200 Pa, ρ = 1060, μ = 3.5 mPa s) and an open inlet (p₁ = p_e, α = 1): `Q_choke` = 1.999 × 10⁻⁵ m³/s, α_out = 0.3123, `p2_choke` − p_e = −946.1 Pa (computed for this design by integrating dx/dα from α = 1 down to the sonic α*, which avoids the stiff end of an x-march; an x-march with a terminal event can step past S = 1 near the inlet-choke flux 26.7 mL/s — bracket Q below it) | ours |

### C.10 Plants (I19) and bookkeeping (I20)
| Function | Signature → return | Eq. |
|---|---|---|
| `water_potential` | `(p, c, T, z, *, rho=1000.0, g=G0)` → Pa (p gauge [Pa], c [mol/m³], T [K], z [m]) | $\psi=p-RTc+\rho gz$ (16.221) |
| `water_potential_terms` (+) | same arguments → dict `pressure`, `osmotic` (= −RTc), `gravity`, `total` | (16.221) |
| `osmotic_pressure` | `(c, T)` → Pa, $\Pi=cR_uT$ (ours, van 't Hoff) | P384 |
| `xylem_reynolds` | `(u, a, nu)` → $Re=u(2a)/\nu$ | $Re=ud/\nu$ |
| `xylem_column` | `(z, *, psi_soil, u, a, mu, rho=1000.0, g=G0, c=0.0, T=293.15)` → dict `z`, `p`, `psi`, `gravity`, `osmotic`, `viscous`, `grad_viscous` ($8\mu u/a^2$) (ours) | ours |
| `cavitation_margin` | `(p, *, p_cav)` → p − p_cav (p_cav required, illustrative) | ours |
| `sieve_tube_pressure_drop` | `(Q, mu, a, r_pore, L, ell, N)` → dict `total`, `tube`, `plate_viscous`, `plate_entry` [Pa] | (16.222) |
| `book_slips`, `traps`, `conventions_table`, `illustrative_inputs`, `units_table` | `()` → DataFrame / dict | ours |

### C.11 Flag list — additions to analysis §4
(+) items: `mL_to_m3`, `pv_loop_schematic_work`, `VESSELS`, `casson_flux_ratio`, `windkessel_summary`, `womersley_regime`,
`bessel_j0_series`, `bessel_j1_series`, `womersley_velocity_nd`, `womersley_flow_rate_harmonics`, `dean_secondary_velocity`,
`water_potential_terms`; from curation §8 (not in analysis §4): `differential_viscosity`, `windkessel_periodic_state`,
`windkessel_terms`, `wall_bulge`, `pulse_state`, `murray_terms`, `bifurcation_cost`, `murray_symmetric_angle`,
`junction_from_geometry`, `dean_axial_velocity`, `dean_series_reliability`, `collapsible_terms`, `starling_state`,
`xylem_column`, `cavitation_margin`. Changed from analysis §4: `windkessel_periodic_p0` takes the half-sine parameters by keyword
(closed form) instead of `Q_fn`; `pv_loop_schematic` takes `p_min, p_ED, p_open, p_peak, p_ES` (the analysis named three);
`womersley_flow_rate` is single-mode (the harmonic sum is `womersley_flow_rate_harmonics`); **`dean_axial_correction` has the
plus sign**; `flow_limitation_curve(p2_list, *, p1, pe, prm)` replaces the analysis's `(dp_list, p1_minus_pe, L, prm)`.

### C.12 Parity conventions (the explainers reproduce these exactly)
- Half-sine ejection: $Q(t)=Q_{\max}\sin(\pi t/T_s)$ for $0\le t<T_s$, 0 otherwise, period $T=60/HR$, $T_s=$ systole_share·T,
  $Q_{\max}=\pi SV/(2T_s)$; periodic state in closed form: with $\tau=RK$, $\omega_s=\pi/T_s$,
  $B=\dfrac{RQ_{\max}\omega_s\tau}{1+(\omega_s\tau)^2}$, $p_0=\dfrac{B(1+e^{-T_s/\tau})e^{-(T-T_s)/\tau}}{1-e^{-T/\tau}}$; systole
  $p=\dfrac{RQ_{\max}(\sin\omega_st-\omega_s\tau\cos\omega_st)}{1+(\omega_s\tau)^2}+(p_0+B)e^{-t/\tau}$; diastole
  $p=[B+(p_0+B)e^{-T_s/\tau}]e^{-(t-T_s)/\tau}$.
- Womersley JS: power series for J₀, J₁ with **100 terms** for every α on the sliders (α ≤ 25 in E4, ≤ 30 in E5): the terms of
  $\sum(-1)^m(z/2)^{2m}/(m!)^2$ at $z=i^{3/2}\alpha$ rotate in phase instead of cancelling, so double precision keeps ~13 digits there.
  The asymptotes (16.88) and (16.96) are drawn only as ghosts — (16.88) is 5 % off the exact profile at α = 20, r = 0.9a (it drops
  the wall's curvature). Parity rows at α = 0.5, 3, 8, 12 (rtol 1e-9) and α = 20 (rtol 1e-8).
- Elastic-tube root in JS: complex quadratic formula, root with |X − 1| smallest at α ≥ 3 (equivalent to continuation there);
  parity at α = 3, 8, 20 (rtol 1e-6).
- Dean in JS: the closed polynomials of C.8; flux with the exact coefficient 1541/50400 and the printed 0.0120.
- Collapsible in JS: the α-integration of C.9's march (dx/dα form, Simpson with 400 panels); parity of the choking flux to rtol 1e-4.
- Expected values used in parity rows are listed in each Part B block.

### C.13 Reconciliation with the module as found (probed 2026-10-14, end of the design run)
`fluidpy/ch16_biofluid_mechanics.py` (3103 lines) appeared while Parts A, B and F were being written. Every function named in
C.1–C.10 exists with the signature above (extra keyword arguments and extra return keys are allowed and listed here); about sixty
calls were made with this design's inputs:
- **Agree to the digits quoted in this design:** SV, EF, CO; the schematic loop (work 1.066573 J by the polygon, 1.066579 J exact,
  `orientation` "counterclockwise"); Casson stress, apparent and differential viscosity, Whitmore's τ_y and K_c; the yield-stress
  correlation; α of the artery (2.8936); R; the two-layer and thin-layer μ_app; the tube law, compliance and the elastic-flux ratio
  (1.5311); the Windkessel p₀ 74.919, p_S 106.028, p_D 74.515, mean 89.726, ⅓ rule 85.019 mmHg and the storage/run-off split at T_s/2;
  Moens–Korteweg 4.8564 and 4.7114 m/s; the wall bulge; the linearisation ratios; Womersley u at α = 3, 8, 20, |Q|/|Q₀| = 0.10480 and
  φ = 78.94° at α = 8, the (16.76) residual −285.71 = −G/μ; the elastic-tube pressure root at α = 8 (c/c₀ 0.95406, k₂/k₁ 0.062169,
  0.67664 per wavelength); Murray's radius (5.8139 mm; printed 58.139 mm) and the 1 : 2 cost split; the 70/30 angles; the junction R,
  T and energy; ψ₀ and **w₁ with the plus sign** (0.00093079 at r = 0.5, α = 0); Dean's flux and λ; the collapsible wave speeds,
  choking flux and M = 0.5; the Casson pipe (a_c, u_c, Q, ratio 0.19029, Δp_min 9.6 Pa); the water potential, osmotic pressure, xylem
  column (−0.7782 MPa) and sieve-tube share (0.8471).
- **Extra names in the module** (usable, not required): `casson_constants`, `casson_stress_plasma_form`, `casson_poiseuille_flux_ratio`,
  `power_law_pipe_velocity`, `two_layer_total_flow_rate`, `windkessel_periodic_p0_general`, `windkessel_periodic_summary`,
  `windkessel_terms_general`, `velocity_from_pulse_shapes`, `tube1d_pressure_flux`, `womersley_profile_complex`, `bessel_ratio_j0_j1`,
  `alpha_from_tube_law`, `collapsible_tube_law_derivatives`, `collapsible_influence_coefficients`, `collapsible_frictionless_alpha`,
  `collapsible_speed_index_rhs`, `friction_R_per_length`, `water_potential_book_units`, `VESSELS_ILLUSTRATIVE`, `DEAN_BOOK`,
  `DEAN_C4_DERIVED`, `CollapsibleTube`, `HalfSineEjection`. The `prm` of C.9 may be a dict with the keys of `_TUBE_KEYS` (`A0`, `Kp`,
  `pe`, `dzdx`, `rho`, `mu`, `n`, `form`, `friction`, `f_L_prime`, `g`, `h_fd`) plus `L`, or a `CollapsibleTube`.
- **Three differences to settle (reported, not edited):** (1) `flow_limitation_curve` returned Q_choke = 2.00290 × 10⁻⁵ m³/s and
  p2_choke = −932.8 Pa for the vein of header convention 10, against 1.99899 × 10⁻⁵ m³/s and −946.1 Pa from this design's
  α-integration (L_needed(Q) = L, monotone). The module's bisection relies on an x-march with a terminal event, which can step past
  S = 1 close to the inlet's own choking flux (26.7 mL/s); the verifier should compare the two routes. E9's parity row carries rtol
  5e-3 until then. (2) `ch16.illustrative_inputs()` is a flat dict whose values differ from this design's dicts (e.g. SV 78 mL, R
  1.4 × 10⁸, K 1.1 × 10⁻⁸, a 11 mm aorta, ψ_soil −0.3 MPa, p_cav −3 MPa): the notebook uses the dicts of header convention 10
  literally and says "both are ours"; aligning the module's set with these dicts would remove the confusion. (3)
  `dean_series_reliability` defaults to `threshold=0.1`; the notebook and E8 pass `threshold=0.25` explicitly (K ≲ 460).
- `ch16.vessel_scales_table(HR=66.0)` has seven rows (aorta, large artery, arteriole, capillary, venule, vein, vena cava; α = 15.9,
  2.17, 0.029, 0.0058, 0.017, 2.89, 20.3): E4's vessel chips and C03's figure use those rows (Part B E4 updated).

---

## Part A — notebook storyboard (`notebooks/build_ch16.py` → `notebooks/ch16_biofluid_mechanics.ipynb`)

**One line per book section** (cell numbers are estimates; ≈ 560 cells in all):
- §16.1 → N01 (the chapter opener; opens the C01 story) — cells ≈ 12–14
- §16.2 → C01 (N02 N03 N04 N05 N06 N07 N08 N09 N10 N11 N12 N13 N14 · D01 · P378 P394 P381 · loop figure, loop animation, series-pressure
  and continuity figures), C02 (N15 N16 N17 N18 N19 N20 N21 N22 N23 · R01 · D02 · P383 · F1) — cells ≈ 15–95. §16.2 rows taught later
  (curation decision 4): N24 → C03, N29 N30 → C05, N48 → C07.
- §16.3 → C03 (N24 N25 · R02 · D03 · vessel figure), C04 (N26 N27 N28 · R03–R08 · D04 · resistance figure), C05 (N29–N36 · D05 D06 D07
  · F2), C06 (N37 N38 N39 N40 N41 N42 · D08 · P379 P380 · F3), C07 (N43–N48 · D09 D10 · P382 · A1 · live · **E2**), C08 (N49–N54 · R09 R10 ·
  D11 D12 · A2 · **E3**), C09 (N55–N60 · R11 · D13 D14 D15 · P374 P375 P377 · A3 · **E4**), C10 (N61–N64 · D16 D17 · P376 · **E5**), C11
  (N65–N73 · R12 · D18 D19 D20 D21 · P385 P386 · F4 · A5), C12 (N74 N75 N76 · D22 D23 D24 · P387 P388 · F5 · **E6**), C13 (N77–N81 · D25
  D26 · A4 · **E7**), C14 (N82–N93 · D27 · P389 P390 P391 · F6 · **E8**), C15 (N94–N102 · R13 · D28 D29 D30 · P392 · **E9**), C16 (N103–N108
  · D31 D32 D33 · P393 · F7 · **E1**) — cells ≈ 96–520
- §16.4 → C17 (N109–N113 · D34 · P384 · F8 · **E10**) — cells ≈ 521–550
- Exercises, literature → S01, S02; references; summary — cells ≈ 551–560

Every CORE block follows: problem in plain words → idea → primers → maths (notes and derivations, Part F) → tiny example →
code (fluidpy) + "What does the code above do?" → from-scratch check → visual(s) → notes and "What would change if…". Code drafts
give intent + exact calls; the builder comments every line (novice grade, units, the equation written out next to its number).
*expect* = numbers the executed cell must print (hypotheses to report against, not to edit). *see / read / change* = the three
figure notes (the builder rewrites them from the rendered figure — CUMULATIVE rule). Note ids open every note in bold (**N09
[B]**). Wherever a draft names an equation, the equation is written out with its number, and the builder keeps it so. The
literal dicts of header convention 10 supply every worked number.

**Re-entry rule (binding; `tools/nbkit.py`).** `nb.recap(...)` and `nb.section(...)` close the current CORE block. The builder
writes `nb.current_core = "Cnn"` on the line after every `nb.recap` that sits inside block Cnn; rows marked **(re-enter Cnn)**
say where. Reminders of earlier-chapter ideas that are not inventory rows of this chapter (Poiseuille's free body, the
control volume, hydrostatics, d'Alembert's characteristic coordinates, the stream function, the Laplacian in cylindrical
coordinates, Stokes' second problem) are `nb.md` cells in the 🔁 recap format without an id; they do not close the block.

### A.0 Front matter
1. `nb.title(big_idea=…, roadmap=[…17…], prerequisites=[…])`. **Big idea (draft):** "The last chapter turns the whole book on your
   own body. Blood is a slightly odd fluid, pushed in pulses through elastic, branching, curved and sometimes collapsing tubes —
   and almost every answer is either Poiseuille's resistance, $R=\dfrac{8\mu L}{\pi a^4}$ (16.20), or a wave whose speed the wall sets, $c=\sqrt{\dfrac{Eh}{2\rho a_0}}$ (16.62), adjusted by one
   number, the Womersley number $\alpha=a\sqrt{\omega/\nu}$ (16.12), that compares a vessel with its own Stokes layer. At the end the same tools lift water
   up a tree." **Roadmap (one line per CORE):** C01 the heart as a pump: work = loop area · C02 blood thins when sheared (Casson) ·
   C03 one number for every vessel: α · C04 resistance ∝ a⁻⁴ · C05 a plasma layer makes narrow tubes slippery · C06 the tube law:
   an artery is a spring · C07 the Windkessel: store and drain · C08 the pulse is a wave of wall stretch · C09 Womersley's
   oscillating profile · C10 how much flow, how late · C11 viscous pulse waves: slower and fading · C12 Murray's law: Q ∝ a³ · C13
   reflections at forks · C14 Dean's vortices in a bend · C15 collapsible tubes choke at S = 1 · C16 the Casson plug · C17 water
   potential: sap is pulled up. **Prerequisites:** Newton's viscosity law $\tau=\mu\,du/dy$ (1.3) and hydrostatics $dp/dz=-\rho g$ (1.8) (Ch. 1);
   continuity and Navier–Stokes in cylindrical coordinates (Ch. 4); Poiseuille flow $u_z=\dfrac{R^2-a^2}{4\mu}\dfrac{dp}{dz}$ (8.6) and Stokes' second problem
   $u(y,t)=U\exp\{-y\sqrt{\omega/2\nu}\}\cos(\omega t-y\sqrt{\omega/2\nu})$ (8.38) (Ch. 8); boundary-layer separation (Ch. 9); d'Alembert's solution $p'(x,t)=f(x-ct)+g(x+ct)$ (15.15) and the 1 − M² structure
   of duct flow (Ch. 15).
2. `nb.explainer_index([...])` — 10 rows: ("casson_blood_flow", "Does a yield stress stop the blood?", "the centre freezes into a
   plug; the flow stops only when the wall stress drops below τ_y") · ("windkessel_heartbeat", "Why doesn't blood pressure fall to
   zero between beats?", "a stretchy aorta stores each stroke and drains it with time constant RK") · ("pulse_wave_artery", "Is the
   pulse at your wrist the blood arriving?", "a wave of wall stretch at metres per second; the blood barely moves") ·
   ("womersley_profiles", "Is pulsing flow just Poiseuille flow sloshing?", "α turns a breathing parabola into a lagging plug") ·
   ("womersley_flow_phase", "Does the flow follow the pressure?", "small α a damper, large α a mass: 1/α² and 90° late") ·
   ("murray_branching", "Why isn't a wider vessel always better?", "pumping ∝ a⁻⁴ against upkeep ∝ a²: Q ∝ a³ and fixed angles") ·
   ("junction_reflection", "When does a fork send the pulse back?", "an admittance mismatch, not size alone; energy always balances") ·
   ("dean_curved_pipe", "What does a bend do to the flow?", "two vortices, a core flung outward, a flux loss of order K²") ·
   ("collapsible_tube_choking", "Why can't you suck harder through a floppy straw?", "at S = 1 the flow stops hearing the outlet") ·
   ("xylem_water_potential", "How does water climb a 30-metre tree?", "it is pulled: ψ falls from soil to leaf and p goes below zero").
3. `nb.setup()`.
4. `nb.code` — **chapter imports**: `import numpy as np`, `import sympy as sp`, `import matplotlib.pyplot as plt`, `import pandas as pd`,
   `from scipy.integrate import quad, solve_ivp`, `from scipy.optimize import brentq, minimize`, `from scipy.special import jv`,
   `from fluidpy import ch16_biofluid_mechanics as ch16`, `from fluidpy.core import laminar as LAM, acoustics as AC, similarity as
   SIM, thermo as TH, units as UN`, `from fluidpy import ch01_introduction as ch01`, `from fluidpy.core.interact import slider_figure,
   animate_figure, live`, `from fluidpy.core.anim import animate`, `from fluidpy.core.style import COLORS, savefig`; then the dicts of
   header convention 10, literally; `pd.Series(ch16.illustrative_inputs())` shown once. *explain:* one line per import ("`ch16`
   re-exports the Ch. 8 Poiseuille and Stokes-layer functions and the Ch. 15 wave functions it uses"). *expect:* MMHG = 133.322387415.
5. `nb.md` — **⚠️ Conventions in this chapter** (before C01): (a) `ch16.conventions_table()` — the rows of header convention 4 (time
   conventions · G = −dp/dx · τ a magnitude · which Re and St · which Dean number · secant or tangent · k_c or K_c · P negative in
   collapse) with "here" and "elsewhere" side by side (traps **T1**, **T4**, **T5**, **T6**, **T8**, **T27**, **T14**, **T9**); (b) the letter
   table of header convention 5 (trap **T2**); (c) `pd.DataFrame(ch16.book_slips())` (33 rows, `kind` and `changes_physics` shown) and
   `pd.DataFrame(ch16.traps())` (27 rows), with "slips are false as printed; traps are true but easy to misread; eleven slips
   change the physics or a number and are taught corrected — for example $a=\left[\dfrac{16\mu}{\pi^2K}\right]^{1/6}Q^{1/3}$ (16.140, corrected) and
   $\dfrac{d^2U}{dr^2}+\dfrac1r\dfrac{dU}{dr}-\dfrac{i\omega\rho}{\mu}U=-\dfrac G\mu$ (16.76, right side corrected)"; (d) "**Which inputs?** Every worked number uses our own illustrative inputs (the dicts above), never
   the book's"; (e) "**Which units?** SI in every function; mmHg, poise, dyn/cm², bar and the PRU only through converters (P394 in
   C01)"; (f) one line: "**A correction to our own analysis:** the first axial correction of Dean's series has a plus sign (C14, D27)."

### A.1 §16.1 Introduction — N01
1. `nb.section("16.1", "Introduction", intro="**What is this section about?** Biofluid mechanics is the whole book applied to living
   things: blood in vessels from the aorta to capillaries, and water rising in plants. Nothing new is assumed beyond Chapters 1–15;
   what is new is the material (blood, elastic walls) and the questions.")`
2. `nb.note` — **N01 [C]** "The scales run from a red cell (micrometres) to the aorta (centimetres) and the whole circulation (metres).
   Every block below picks one scale: the heart and its loop (C01), the blood itself (C02), single vessels (C03–C11, C16), branching
   and bending trees (C12–C15), and finally plants (C17)."

### A.2 §16.2 The Circulatory System in the Human Body — C01, C02
#### C01 — The heart as a pump: work = the loop's area, $W=-\oint p\,dV$ (ours)
1. `nb.section("16.2", "The Circulatory System in the Human Body", intro="**What is this section about?** The pump (the heart), the
   pipes (vessels) and the fluid (blood): the vocabulary of physiology turned into numbers we can compute with — stroke volume,
   cardiac output, peripheral resistance, and blood's odd viscosity.")`
2. `nb.core("C01", "The heart as a pump: $SV=EDV-ESV$ (16.1) and the work $W=-\oint p\,dV$ (ours)", question="How much blood does the heart
   move each beat, and how much work does that cost it?")`
3. `nb.md` — **Plain words:** "Squeeze a water balloon through a narrow neck, let it refill, squeeze again. That is the left
   ventricle, about once a second: it fills, closes its valves, squeezes until it can open the outlet valve, pushes blood into the
   aorta, closes, relaxes. Draw its pressure against its volume over one beat and you get a closed loop. The loop's width is the
   blood moved, its area is the work done — the same p dV work as a piston in Chapter 1."
4. `nb.note` — **N02 [C]** our own two-loop ASCII sketch:
   ```
   left heart ──► aorta ──► body (systemic) ──► veins ──► right heart ──► lungs (pulmonary) ──► left heart
   ```
   "A closed system: the same blood goes round both loops. The pulmonary loop returns in C06 (N40–N42), branching in C12 and C13."
5. `nb.note` — **N04 [B]** the four phases as a table (phase · valves · volume · pressure): filling A→B (mitral open, aortic shut; V up;
   p low) · isovolumic contraction B→C (both shut; V fixed; p rises) · ejection C→D (aortic open; V down; p near the aortic pressure)
   · isovolumic relaxation D→A (both shut; V fixed; p falls). "Systole (contraction and ejection) takes about a third of a beat; the
   half-sine ejection used in C07 lasts that long (our `systole_share` = 0.35)."
6. `nb.note` — **N03 [C]** "A resting adult heart beats about once a second; the cardiac cycle is one beat. The book quotes resting
   values; our worked numbers use our own illustrative adult (header convention 10)."
7. `nb.primer("the shoelace formula", "The area of a polygon from its corners $(x_i,y_i)$: $A=\\tfrac12\\sum_i(x_iy_{i+1}-x_{i+1}y_i)$, positive when
   the corners run counterclockwise. It is the discrete form of the loop area of Ch. 5 (P138) and of the line integral ∮y dx (P35).
   `np.roll(x, -1)` gives the next corner of every corner at once.", code="x = np.array([0, 2, 2, 0]); y = np.array([0, 0, 1, 1])   #
   a 2 × 1 rectangle, counterclockwise\nprint(0.5*np.sum(x*np.roll(y, -1) - np.roll(x, -1)*y))      # 2.0 (clockwise would give −2.0)")`
   (**P378**; reminder lines for P35 line integral (ch01) and P138 loop area (ch05)).
8. `nb.primer("units of physiology", "Clinicians measure pressure in mmHg (1 mmHg = 133.322387415 Pa, the weight of a 1 mm mercury
   column), viscosity in poise or centipoise (1 cP = 10⁻³ Pa s), stresses in dyn/cm² (0.1 Pa), plant pressures in bar (10⁵ Pa) and
   resistance in PRU (1 mmHg per cm³/s = 1.333 × 10⁸ Pa s m⁻³). Every function here works in SI; converters do the rest.",
   code="print(ch16.mmHg_to_Pa(120), ch16.cP_to_Pa_s(3.5), ch16.dyn_cm2_to_Pa(1.0), ch16.PRU)   # 15999 Pa, 0.0035 Pa s, 0.1 Pa,
   1.333e8\nprint(UN.Q_(120, 'mmHg').to('Pa'))                                                 # pint agrees")` (**P394**; trap T3).
9. `nb.derivation("D01", …)` — Part F D01 (6 steps), ref "ours".
10. `nb.note` — **N06 [B]** "Stroke volume $SV=EDV-ESV$ (16.1): the loop's width. Ours: 130 − 55 = 75 mL." + **N07 [B]** "$EF=(SV/EDV)\times100\%$ (16.2), $CO=SV\times HR$ (16.3),
    $CI=CO/BSA$ (16.4) — the unit chain written out: 75 mL/beat × 66 beats/min = 4950 mL/min = 4.95 L/min = 8.25 × 10⁻⁵ m³/s; per square metre of
    body surface (1.9 m²) 2.61 L/(min m²). (The book once writes 'mL/m' for mL/min — slip S32.)" + `nb.code`: `SV =
    ch16.stroke_volume(heart["EDV"], heart["ESV"])`, `ch16.ejection_fraction(SV, heart["EDV"])`, `CO = ch16.cardiac_output(SV, heart["HR"])`,
    `ch16.cardiac_index(CO, heart["BSA"])`. *expect:* 7.5e-5 m³, 0.577, 8.25e-5 m³/s, 4.34e-5 m³/(s m²).
11. `nb.worked_example("a rectangular loop, by hand", "1. Take a loop 100 mmHg tall and 75 mL wide. 2. Area = 100 × 75 = 7500 mmHg·mL. 3.
    One mmHg·mL = 133.32 × 10⁻⁶ = 1.333 × 10⁻⁴ J (P394). 4. W = 7500 × 1.333 × 10⁻⁴ = **1.000 J**. 5. At 66 beats per minute (1.1 per
    second) that is **1.1 W** of hydraulic power — about as much as a small torch bulb.")`
12. `nb.code` — `loop = ch16.pv_loop_schematic(heart["EDV"], heart["ESV"], p_min=heart["p_min"], p_ED=heart["p_ED"], p_open=heart["p_open"],
    p_peak=heart["p_peak"], p_ES=heart["p_ES"])`; `w = ch16.pv_loop_work(loop["p"], loop["V"])`; `W_exact = ch16.pv_loop_schematic_work(...)`;
    `ch16.cardiac_power(w["work"], heart["HR"])`. *expect:* work 1.0666 J (exact 1.06658 J), orientation "counterclockwise", power ≈ 1.17 W.
    *explain:* 4 points (the loop's four legs; the trapezoid ∮; the closed form; power = W × beats per second).
13. `nb.check_agree` — **from scratch:** `V, p = loop["V"], loop["p"]`; `area = 0.5*np.sum(V*np.roll(p, -1) - np.roll(V, -1)*p)` (P378, positive
    counterclockwise in (V, p)); `assert np.isclose(area, w["work"], rtol=1e-9)`; a rectangle loop `ch16.pv_loop_work(np.array([…]),
    np.array([…]))["work"]` equals Δp ΔV to 1e-12; reversed order gives −W.
14. `nb.figure` — the PV loop (V in mL, p in mmHg): the four legs coloured by phase, corners A–D labelled, area shaded amber with "W =
    1.07 J", width arrow "SV = 75 mL", arrows showing the counterclockwise direction; a second panel with the rectangle of the tiny
    example. *see:* "…a closed loop traced counterclockwise…" · *read:* "…width = SV; height ≈ the pressure the ventricle must reach;
    area = work; the vertical sides do no work…" · *change:* "…a stiffer aorta (C07) raises the ejection leg and the area: the heart
    works harder for the same stroke."
15. `nb.animation` — the loop traced over one beat (60 frames; FAST 30; video): a dot on the loop, the ventricle drawn as a shrinking
    and growing chamber beside it, valve states written; the "so far" part bold. *see / read / change.*
16. `nb.note` — **N08 [B]** units box: "pressures are quoted in mmHg; our illustrative 125/75 mmHg (systolic/diastolic, ours) =
    16.7/10.0 kPa" + `ch16.mmHg_to_Pa(np.array([125, 75]))`. **N09 [B]** "pulse pressure $p_S-p_D$ = 50 mmHg (`ch16.pulse_pressure`) — the
    size of the pressure wave; C07 shows how compliance sets it, C08 how it travels."
17. `nb.note` — **N10 [B]** "Mean pressure falls continuously from the aorta to the veins: each vessel class is a resistance in series,
    $\Delta p_i=QR_i$ (ours)." + `nb.figure` — bar chart of the pressure after each class from `ch16.series_pressure_drop([…], CO, 92*MMHG)` with our
    illustrative resistance shares (aorta and large arteries 10 %, small arteries 15 %, arterioles 45 %, capillaries 20 %, venules and
    veins 10 %; labelled ours) — *see / read / change* ("…the arterioles take the biggest drop: that is where the body controls its
    pressure — C04 says why").
18. `nb.note` — **N11 [B]** "Continuity across the tree: the mean speed is $\bar u=Q/A_{tot}$ (ours). The capillaries are slow because there are so
    many of them." + `nb.figure` — total area (log) and mean speed per class from `ch16.mean_speed(CO, A_tot)` with cited typical areas
    (`ch16.VESSELS`, labelled cited/ours) + 🔁 "mass conservation (Ch. 4): what goes in comes out, so speed × area is the same at every
    level of the tree".
19. `nb.primer("the electrical analogy", "A pressure difference drives a flow through a resistance as a voltage drives a current:
    Δp ↔ V, Q ↔ I, R = Δp/Q ↔ Ohm's law. A vessel that stores volume when pressurised is a capacitor (compliance K = dV/dp, C07); a
    tube's wave admittance Y = Q/p is the inverse of an impedance (C13). Resistors in series add; in parallel their inverses
    add.", code="R = np.array([1.0, 2.0, 4.0])                  # three resistances (any units)\nprint(R.sum(), 1/np.sum(1/R))        # series
    7.0, parallel 0.571")` (**P381**).
20. `nb.note` — **N12 [B]** "Total peripheral resistance $R=\Delta\bar p/Q$ (16.9) — the Ohm's law of the whole body; unit PRU = mmHg per cm³/s." + **N13 [B]**
    "$\Delta p=\bar p_A-\bar p_V$ (16.10) with $\bar p_V$ taken as 0 (trap T12)." + **N14 [B]** "$\bar p_A=\tfrac13p_S+\tfrac23p_D$ (16.11): an estimate of the mean for a typical waveform (trap
    T15)." + `nb.code`: `pA = ch16.mean_arterial_pressure(heart["pS"], heart["pD"])`; `R_body = ch16.peripheral_resistance(pA, CO)`; `R_body /
    ch16.PRU`. *expect:* 91.7 mmHg; 1.48 × 10⁸ Pa s m⁻³ = 1.11 PRU (ours). "⚠️ Common confusion (T15): C07 computes the true mean of a
    model waveform — the ⅓ rule is 5 % low there."
21. `nb.note` — **N05 [C]** "Only a small fraction of the heart's metabolic power becomes pressure–volume work; most becomes heat. The
    work per beat above (≈ 1 J) times the heart rate is the useful part."
22. `nb.md` — **What would change if…** "…the aorta were a rigid pipe? Then the pressure would follow the heart's ejection exactly and
    fall to zero between beats. It does not — because the aorta stretches. C07 (the Windkessel) makes that precise; first, the fluid."

#### C02 — The Casson law $\sqrt\tau=\sqrt{\tau_y}+K_c\sqrt{\dot\gamma}$ (16.6)
1. `nb.core("C02", "Blood thins when sheared: the Casson law $\sqrt{\tau/\mu_p}=k_c\sqrt{\dot\gamma}+\sqrt{\tau_y/\mu_p}$ (16.6)", question="Is blood just thick water?")`
2. `nb.md` — **Plain words:** "Honey is thick, water thin, and each has one viscosity whatever you do to it. Blood is different:
   stir it slowly and it seems very thick, stir it fast and it thins. At rest its red cells stack into coin-like piles (rouleaux)
   that resist; shearing breaks them up. And a tiny stress is needed before it flows at all — a yield stress."
3. `nb.md` — **The idea** (ASCII): `√τ ↑   blood: a straight line that starts at √τ_y` / `      |  Newtonian: a line through 0` / `      +──────────► √γ̇`
   "— a straight line in square roots, with an intercept."
4. `nb.note` — **N15 [B]** "Blood in numbers (cited typical values, `ch16.BLOOD`): ρ = 1060 kg/m³, high-shear viscosity 3.5 mPa s, plasma
   1.2 mPa s." `pd.Series(ch16.BLOOD)`. **N16 [B]** "A red cell is a flexible disc with far more surface than a sphere of the same volume
   (`ch16.equal_volume_sphere_area`; a 90 µm³ cell: sphere 97 µm², ours) — it can fold through capillaries narrower than itself
   (C05)." **N17 [B]** "Hematocrit H: the volume fraction of red cells, about 0.45 (45 %). Fraction in code, percent in prose (trap T13):
   `ch16._check_hematocrit(45)` raises."
5. `nb.note` — **N18 [B]** "Non-Newtonian: the ratio of stress to shear rate is not a constant. The units box: $1\,\mathrm{cP}=10^{-2}\,\mathrm P=10^{-2}\,\mathrm{dyn\,s/cm^2}=10^{-3}\,\mathrm{Pa\,s}$." + `> ⚠️ **slip S27
   — the book prints** 1 cP = 0.01 P = 1 dyn s/cm² = 0.1 N s/m² **; the correct chain is** 1 cP = 10⁻² P = 10⁻² dyn s/cm² = 10⁻³ Pa s — the
   printed one is 100× too large.` + `nb.code`: `ch16.cP_to_Pa_s(1.0)`, `UN.Q_(1, "cP").to("Pa*s")`. *expect:* 0.001 both.
6. `nb.recap("R01", "Newton's law of viscosity", "$\\tau=\\mu\\,du/dy$ (1.3): stress proportional to shear rate, one constant μ. **New here:**
   for blood τ/γ̇ is not constant — the Newtonian line becomes the ghost every figure of this block compares with.", where="Ch. 1 §1.5
   (C12)")` + `nb.code`: `ch01.newton_shear_stress(3.5e-3, 100.0)`, `ch16.newtonian_stress(100.0, 3.5e-3)`. *expect:* 0.35 Pa both.
   **(re-enter C02)**
7. `nb.note` — **N19 [B]** "The yield stress: below τ_y blood does not shear at all. An empirical correlation $\sqrt{\tau_y}=(H-0.1)(C_F+0.5)$ (16.5) (τ_y in dyn/cm², H a
   fraction above 0.1, fibrinogen C_F in g/100 mL — the hidden units said aloud)." + `ch16.yield_stress_correlation(0.45, 0.3)`. *expect:*
   7.84 × 10⁻³ Pa (0.0784 dyn/cm²). "The book's quoted range at H = 0.45 differs from what (16.5) gives for 0.21 < C_F < 0.46 (S29)."
8. `nb.note` — **N20 [C]** "Why it thins: rouleaux break up as the shear rate rises; above a threshold the cells move singly. Drawn as two
   regions on the figure below; the asymptote is K_c² (D02)."
9. `nb.primer("secant slope and tangent slope", "On a curve y(x) the secant from the origin to a point has slope y/x; the tangent at
   that point has slope dy/dx (P19). For a straight line through 0 they agree; for a bending curve they do not. Blood's 'apparent
   viscosity' τ/γ̇ is the secant.", code="x = 4.0; y = lambda x: (0.07 + 0.053*np.sqrt(x))**2      # a Casson-like curve\nh = 1e-6\nprint(y(x)/x, (y(x + h) -
   y(x - h))/(2*h))      # secant 0.00768 > tangent 0.00465")` (**P383**; trap T27).
10. `nb.note` — **N21 [B]** "Apparent viscosity $\mu_{app}=\tau/\dot\gamma$ — the secant (`ch16.apparent_viscosity`); the tangent dτ/dγ̇ is
    `ch16.differential_viscosity`."
11. `nb.derivation("D02", …)` — Part F D02 (5 steps), ref "16.6".
12. `nb.note` — **N22 [B]** "Whitmore's fit $\sqrt{\tau/\mu_p}=1.53\sqrt{\dot\gamma}+2.0$ (16.7): k_c = 1.53 (a pure number) and √(τ_y/μ_p) = 2.0 s^{−1/2}, so τ_y = 4μ_p and
    K_c² = 1.53² μ_p = 2.34 μ_p. A two-row table (trap T14): k_c dimensionless · K_c = k_c√μ_p in (Pa s)^½." + `nb.code`: `wc =
    ch16.whitmore_casson(np.array([4.0, 400.0]), blood["mu_p"])`. *expect:* tau_y 4.80 mPa, Kc 0.05300, mu_app 7.68 and 3.19 mPa s.
13. `nb.worked_example("apparent viscosity at two shear rates", "With Whitmore's constants and μ_p = 1.2 mPa s: 1. √τ_y = √(4.8 ×
    10⁻³) = 0.0693; K_c = 1.53 × √0.0012 = 0.0530. 2. γ̇ = 4 s⁻¹: √τ = 0.0693 + 0.0530 × 2 = 0.1753, τ = 30.7 mPa, μ_app = 30.7/4 = **7.68
    mPa s** (6.4 μ_p). 3. γ̇ = 400 s⁻¹: √τ = 0.0693 + 0.0530 × 20 = 1.129, τ = 1.275 Pa, μ_app = **3.19 mPa s** (2.66 μ_p). 4. γ̇ → ∞:
    μ_app → K_c² = **2.81 mPa s**. Slow shear makes blood look 2.7 times as viscous as fast shear.")`
14. `nb.code` — `tau = ch16.casson_stress(gd, wc["tau_y"], wc["Kc"])` on `gd = np.logspace(-1, 3, 200)`; `ch16.casson_apparent_viscosity`;
    `ch16.differential_viscosity`; round trip `ch16.casson_shear_rate(tau, …)`. *expect:* the round trip returns gd to 1e-12; below
    τ_y the shear rate is 0.
15. `nb.check_agree` — **from scratch:** `tau_mine = (np.sqrt(ty) + Kc*np.sqrt(gd))**2`; `gd_mine = np.where(tau >= ty, ((np.sqrt(tau) -
    np.sqrt(ty))/Kc)**2, 0.0)` (P46 np.where reminder); `assert np.allclose(tau_mine, ch16.casson_stress(gd, ty, Kc))`; `assert
    np.allclose(gd_mine, ch16.casson_shear_rate(tau_mine, ty, Kc))`.
16. `nb.figure` — two panels: (a) √τ against √γ̇ — Casson line (blue) with its intercept √τ_y marked, the Newtonian line through the
    origin (muted dashed); (b) μ_app against γ̇ on a log axis (blue) falling to K_c² (muted dashed), the secant (rose) and tangent
    (purple) values at γ̇ = 4 s⁻¹ marked, the two regions of N20 shaded lightly ("rouleaux", "single cells"; thresholds illustrative).
    Our curves only — no data points (they are the book's). *see / read / change* ("…double τ_y and only the low-shear end rises; the
    asymptote K_c² does not move").
17. `nb.plotly` — **F1** `slider_figure` on τ_y (0…0.02 Pa, 25 steps): μ_app against γ̇ (log x, 200 points), the line μ = K_c², and the
    secant and tangent markers at γ̇ = 4 s⁻¹. *explain:* "every slider position was computed in Python (P17)". *see / read / change.*
18. `nb.note` — **N23 [C]** "Casson's law treats blood as a homogeneous fluid. In tubes below about a third of a millimetre the cells
    leave a plasma layer at the wall — C05 models exactly that."
19. `nb.md` — **What would change if…** "…we pushed this blood through a pipe? The stress is zero on the axis and largest at the wall,
    so near the axis it is below τ_y: the centre moves as a solid plug. C16 works that out; E1 (embedded there) lets you watch it."

### A.3 §16.3 Modeling of Flow in Blood Vessels — C03 … C16
#### C03 — The Womersley number $\alpha=a\sqrt{\omega/\nu}$ (16.12)
1. `nb.section("16.3", "Modeling of Flow in Blood Vessels", intro="**What is this section about?** The long modelling section, read in
   four arcs: steady flow (C04, C05); storage and waves in elastic tubes (C06–C08); oscillating flow (C09–C11); geometry and odd
   materials (C12–C16). One number, α, decides which arc's tools a vessel needs.")`
2. `nb.core("C03", "One number for every vessel: the Womersley number $\alpha=a\sqrt{\omega/\nu}$ (16.12)", question="When can the pulsing flow in a vessel be
   treated as steady, and when not?")`
3. `nb.md` — **Plain words:** "Shake a glass of water gently and the whole glass of water moves together; shake a bathtub and only the
   water near the walls notices the walls. Viscosity spreads a disturbance only so far in one cycle — about √(ν/ω). Compare that
   distance with the vessel's radius: that comparison is α. A capillary is far thinner than this distance (α ≪ 1); the aorta far
   wider (α ≫ 1)."
4. `nb.md` — **The idea** (ASCII): `|←— a —→|` the radius, `|←δ→|` the Stokes layer √(2ν/ω); "α = √2 a/δ_e: small α — the layer fills
   the tube; large α — a thin skin".
5. `nb.recap("R02", "Separation in an adverse pressure gradient", "A boundary layer pushed against a rising pressure can stop and turn
   back (Ch. 9). **New here:** in a pulse the deceleration phase supplies that rising pressure every beat — reversal near the wall in
   C09 is its oscillating cousin.", where="Ch. 9 §9.7")` **(re-enter C03)**
6. `nb.derivation("D03", …)` — Part F D03 (4 steps), ref "16.12". (🔁 one line before it: "Stokes' second problem (Ch. 8): $u(y,t)=U\exp\{-y\sqrt{\omega/2\nu}\}\cos(\omega t-y\sqrt{\omega/2\nu})$ (8.38),
   e-folding depth δ_e = √(2ν/ω), `LAM.stokes_layer`.")
7. `nb.worked_example("α for an artery", "1. Radius a = 2 mm, heart rate 66 per minute → ω = 2π × 66/60 = 6.91 rad/s. 2. ν = μ/ρ = 3.5 ×
   10⁻³/1060 = 3.30 × 10⁻⁶ m²/s. 3. ω/ν = 2.09 × 10⁶ m⁻², √ = 1447 m⁻¹. 4. α = 0.002 × 1447 = **2.9**. 5. δ_e = √(2ν/ω) = **0.98 mm**: the
   viscous layer reaches halfway to the axis — neither limit.")`
8. `nb.code` — `ch16.womersley_number(artery["a"], 2*np.pi*heart["HR"]/60, blood["mu"]/blood["rho"])`; `LAM.stokes_layer(nu, omega)["delta_e"]`;
   `ch16.womersley_from_re_st(500.0, 0.2)`; `ch16.vessel_scales_table(HR=66.0)` (display). *expect:* 2.894; 9.77e-4 m; 5.0; table rows
   capillary α ≈ 0.006 … aorta ≈ 16 (the table's cited 11 mm aorta; our aorta dict, 12 mm, gives 17.4).
9. `nb.note` — **N24 [B]** (from §16.2, out of book order) "Vessel walls have three layers (intima, media, adventitia); sizes run from
   centimetre aortas to micrometre capillaries; veins hold most of the blood. Our scales table uses cited typical sizes (labelled
   cited/ours) and gives Re, α and h/a for each class — the table the figure below plots."
10. `nb.figure` — α (blue) and Re (orange) per vessel class on log axes (vessel class along x), bands "quasi-steady α < 1",
    "transitional", "inertial α > 10". *see / read / change* ("…double the heart rate and every α grows by √2 = 1.41: the bands do not
    move, the dots do").
11. `nb.md` — "⚠️ Common confusion (trap T6): **which Re and St?** a three-row table — the book: Re = 2au/ν, St = 2aω/u, α² = Re·St/4 · other
    sources: α² = 2π Re St with Re, St on other scales · Exercise 16.2: St = aω/U."
12. `nb.note` — **N25 [C]** "Three ideas named here are taught later: compliance (C07, N45), the pulse as a wave faster than the flow (C08)
    and reflections at branches (C13)."
13. `nb.md` — **What would change if…** "…α were large in every vessel? Steady Poiseuille flow would never apply. It applies in the
    small vessels, which is where the body's resistance sits — C04."

#### C04 — Vascular resistance $R=\dfrac{8\mu L}{\pi a^4}$ (16.20)
1. `nb.core("C04", "Resistance grows as $a^{-4}$: $R=\dfrac{\Delta\bar p}{Q}=\dfrac{8\mu L}{\pi a^4}$ (16.20)", question="Why do tiny changes in the radius of small arteries change blood
   pressure so much?")`
2. `nb.md` — **Plain words:** "Pinch a garden hose to half its width and the flow drops far more than by half. Poiseuille's law says by
   how much: the flow through a tube grows as the fourth power of its radius. Write it as 'pressure drop = resistance × flow' and the
   resistance carries a⁻⁴. The body's arterioles use this lever."
3. `nb.recap("R03", "When Poiseuille's law holds", "Steady, fully developed, laminar flow of a Newtonian fluid in a straight rigid
   tube (Ch. 8). **New here:** a checklist of how wrong each assumption is in an artery and which block relaxes it — unsteady (C09),
   elastic (C06), curved (C14), branched (C12, C13), non-Newtonian (C02, C16); $Re=\\rho\\bar ud/\\mu$ with `SIM.reynolds_number`.", where="Ch. 8 §8.2 (C03);
   Ch. 4 (C15)")` + a six-row table. **(re-enter C04)**
4. `nb.recap("R04", "The Poiseuille profile", "$u(r)=\\dfrac{r^2-a^2}{4\\mu}\\left(\\dfrac{dp}{dx}\\right)$ (16.13) — the same as $u_z=\\dfrac{R^2-a^2}{4\\mu}\\dfrac{dp}{dz}$ (8.6). **New here:** written with the overall drop.", where="Ch. 8
   §8.2")` · `nb.recap("R05", "Gradient from the overall drop", "$\\left(\\dfrac{dp}{dx}\\right)=-\\dfrac{\\Delta p}{L}$ (16.14), so $u(r)=\\left(\\dfrac{a^2\\Delta p}{4\\mu L}\\right)\\left(1-\\dfrac{r^2}{a^2}\\right)$ (16.15). G = −dp/dx = Δp/L > 0 (trap T4).", where="Ch. 8")` · `nb.recap("R06",
   "Centreline speed, flux and mean speed", "$u_{max}=\\dfrac{\\Delta p\\,a^2}{4\\mu L}$ (16.16), $Q=\\dfrac{\\pi a^4}{8\\mu}\\dfrac{\\Delta p}{L}$ (16.17), $V=\\dfrac{Q}{\\pi a^2}=\\dfrac{u_{max}}{2}$ (16.18). The integral is D04 step 2.", where="Ch. 8 (D06)")` · `nb.recap("R07",
   "Wall shear stress", "$\\tau=-\\mu\\left(\\dfrac{du}{dr}\\right)\\Big\\rvert_{r=a}=+\\dfrac a2\\dfrac{\\Delta p}{L}$ (16.19, last member's sign corrected — S1): this chapter's τ is a positive magnitude (trap T5), ch08 code
   uses signed τ.", where="Ch. 8 (8.8)")` · `nb.recap("R08", "Entrance length", "A flat inlet profile needs a length to become parabolic;
   wall shear is higher there (Ch. 8 sketch, Ch. 9 boundary-layer growth). **New here:** some vessels are too short ever to be fully
   developed.", where="Ch. 8 §8.2; Ch. 9 (C02)")` + `nb.code` (inside R06/R07): `ch16.poiseuille_velocity(r, artery["a"], 100.0, 0.1,
   blood["mu"])`, `ch16.poiseuille_flow_rate(...)`, `ch16.poiseuille_wall_stress(...)`, `LAM.pipe_flow_rate(...)` — they agree. **(re-enter C04)**
5. `nb.derivation("D04", …)` — Part F D04 (5 steps), ref "16.20".
6. `nb.worked_example("halve the radius", "1. μ = 4 mPa s, L = 0.1 m, a = 1 mm: R = 8 × 0.004 × 0.1/(π × 10⁻¹²) = **1.02 × 10⁹ Pa s m⁻³**. 2. a = 0.5
   mm: R × 2⁴ = **1.63 × 10¹⁰**. 3. Ten thousand such 0.5 mm tubes in parallel: R/10⁴ = 1.63 × 10⁶ — parallel paths divide the
   resistance (P381).")`
7. `nb.code` — `ch16.vascular_resistance(4e-3, 0.1, 1e-3)`, `ch16.vascular_resistance(4e-3, 0.1, 0.5e-3)`, `ch16.network_resistance([R1]*10000,
   arrangement="parallel")`; an arteriole: `ch16.vascular_resistance(blood["mu"], arteriole["L"], arteriole["a"])`. *expect:* 1.019e9, 1.630e10,
   1.630e6; arteriole 4.56e13 Pa s m⁻³.
8. `nb.check_agree` — **from scratch:** `r = np.linspace(0, a, 2001)`; `u = a**2*dp/(4*mu*L)*(1 - r**2/a**2)`; `Q = np.trapezoid(2*np.pi*r*u, r)`;
   `assert np.isclose(dp/Q, ch16.vascular_resistance(mu, L, a), rtol=1e-6)`; `assert np.isclose(ch16.vascular_resistance(mu, L, a/2) /
   ch16.vascular_resistance(mu, L, a), 16.0)`.
9. `nb.figure` — R against a on log–log axes (a from 5 µm to 15 mm, μ = 3.5 mPa s, L = 1 cm), slope −4 marked (P13 reminder), vessel
   classes as dots. *see / read / change* ("…a slope of −4: a 19 % narrowing doubles R").
10. `nb.note` — **N27 [B]** "Entrance length $\dfrac\ell d=0.06\,Re$ (laminar, Re > 50) and $\dfrac\ell d=0.693\,Re^{1/4}$ (turbulent) (16.37)." + a 'which correlation?' table (the
    book's 0.06 and 0.693 against Çengel's 0.0575 Re and 1.359 Re^{1/4}, cited from Wikipedia — S31) + `ch16.entrance_length(500.0, 4e-3,
    regime="laminar")`. *expect:* 0.12 m — longer than many arteries between branches.
11. `nb.note` — **N26 [C]** "Poiseuille's law fits muscular arteries best; heat transfer in perfused tissue is named in the book and
    not treated here." **N28 [C]** "In pulsatile flow the entrance length depends on Re and α; the wall shear settles sooner than the
    centreline speed (α is C03, the oscillating profile C09)."
12. `nb.md` — **What would change if…** "…the blood were not one fluid? In tubes below about 0.3 mm the cells crowd to the centre and
    leave plasma at the wall — then a single μ no longer describes the flow: C05."

#### C05 — The two-layer apparent viscosity (16.33)
1. `nb.core("C05", "A plasma layer lubricates narrow tubes: $\mu_{app}=\mu_p\left[1-\left(1-\dfrac\delta a\right)^4\left(1-\dfrac{\mu_p}{\mu_c}\right)\right]^{-1}$ (16.33)", question="Why does blood seem
   less viscous in narrower tubes?")`
2. `nb.md` — **Plain words:** "Measure blood's viscosity in glass tubes of shrinking diameter and, below about a third of a
   millimetre, it seems to get *thinner*. The red cells drift away from the wall, leaving a thin layer of plain plasma. The thick
   core slides on that runny skin — and the narrower the tube, the larger the skin's share."
3. `nb.note` — **N29 [B]** (from §16.2) "The Fåhræus–Lindqvist effect: μ_app falls with diameter over a range of small tubes, rising
   again in the very smallest. The Fåhræus effect: the tube hematocrit H_T is below the discharge hematocrit H_D (= the feed H_F)."
   (curve drawn below from the cited in-vitro law of Pries et al. via Secomb & Pries 2013, `ch16.pries_relative_viscosity`, labelled
   cited.)
4. `nb.note` — **N30 [B]** (from §16.2) "Red cells and volume are conserved: $Q_FH_F=Q_cH_c,\ Q_c+Q_p=Q_F,\ H_Ta^2=H_c(a-\delta)^2$ (16.8). Eliminating Q_F gives $H_c/H_F=1+Q_p/Q_c$
   (result given): the core is richer in cells than the feed, and the tube holds fewer cells than flow through it, because the
   cell-rich core moves faster than the mean." + `ch16.cell_conservation(Qc, Qp, 0.5, 20e-6, 2e-6)`.
5. `nb.note` — **N31 [B]** "The model: a core of viscosity μ_c inside a plasma layer of viscosity μ_p and thickness δ; stress continuous at
   the interface; regular on the axis. $\tau_{xr}=-\mu_c\dfrac{du^c}{dr}=+\dfrac r2\dfrac{\Delta p}{L}$ (16.21, sign corrected), $\dfrac{du^c}{dr}=0$ at $r=0$ (16.22), $\tau_{xr}\rvert_c=\tau_{xr}\rvert_p$ at $r=a-\delta$ (16.23) — the printed (16.21) has the
   wrong sign on the right (S2)." + our two-layer sketch (ASCII). **N32 [B]** "$\tau_{xr}=-\mu_p\dfrac{du^p}{dr}=+\dfrac r2\dfrac{\Delta p}{L}$ (16.24), $u^c=u^p$ at $r=a-\delta$ (16.25), $u^p=0$ at $r=a$ (16.26)."
6. `nb.derivation("D05", …)` — Part F D05 (9 steps), ref "16.27–16.28": holds **N33 [B]**.
7. `nb.derivation("D06", …)` — Part F D06 (8 steps), ref "16.31": holds **N34 [B]** + `nb.md` `> ⚠️ **slip S3 — the book prints**
   $Q_c=\dfrac{\pi a^2\Delta p}{4\mu_pL}\left[a^2-\left(1-\dfrac{\mu_p}{2\mu_c}\right)\dfrac{(a-\delta)^4}{a^2}\right]$ (16.30) **; the correct form is** $Q_c=\dfrac{\pi(a-\delta)^2\Delta p}{4\mu_pL}\left[a^2-\left(1-\dfrac{\mu_p}{2\mu_c}\right)(a-\delta)^2\right]$ (16.30, corrected) — only the corrected
   one adds to (16.31).` + `nb.code`: `ch16.two_layer_flow_rates(…, printed_16_30=True)["Qc"] - ch16.two_layer_flow_rates(…)["Qc"]`.
8. `nb.derivation("D07", …)` — Part F D07 (6 steps), ref "16.33–16.34".
9. `nb.worked_example("a 40 µm tube", "Relative units μ_p = 1, μ_c = 4, δ/a = 0.1: 1. (1 − 0.1)⁴ = 0.6561. 2. 1 − μ_p/μ_c = 0.75. 3. Bracket 1 −
   0.6561 × 0.75 = 0.508. 4. μ_app = 1/0.508 = **1.97** μ_p — half the core's 4. 5. The thin-layer form (16.34): 4/(1 + 0.4 × 3) = 1.82 (8 %
   low). With our SI inputs (a = 20 µm, δ = 2 µm, μ_p = 1.2, μ_c = 4.0 mPa s): **2.22 mPa s**; in a 500 µm tube with the same δ, 3.73 mPa s.")`
10. `nb.code` — `ch16.two_layer_velocity(r, micro["a"], micro["delta"], dp, L, blood["mu_p"], micro["mu_c"])`, `ch16.two_layer_flow_rates(...)`,
    `ch16.apparent_viscosity_two_layer(...)`, `ch16.apparent_viscosity_thin_layer(...)`. *expect:* 2.219e-3, 2.069e-3 Pa s.
11. `nb.check_agree` — **from scratch:** `quad(lambda r: 2*np.pi*r*ch16.two_layer_velocity(r, …), 0, a)` equals `two_layer_flow_rates(...)["Q"]`
    to 1e-10; μ_app from (16.32) by hand equals the function.
12. `nb.figure` — (a) the two-layer profile across the diameter (core blue, plasma layer purple band) with the kink at a − δ marked,
    Poiseuille with μ_c as a muted ghost; (b) μ_app/μ_p against D (10…1000 µm, log): our model for δ = 2 µm (blue, "model, ours") beside the
    cited Pries curve at H_D = 0.45 (muted, "cited"). *see / read / change* ("…the model gets the direction and rough size; the rise in
    the very smallest tubes needs cells deforming — not in the model").
13. `nb.note` — **N35 [B]** "The core's viscosity from its hematocrit: $\mu_c=\mu_p\dfrac{1}{1-\alpha_cH_c}$ (16.35) with $\alpha_c=0.070\exp\left[2.49H_c+\dfrac{1107}{T}\exp(-1.69H_c)\right]$ (16.36), T in kelvin, valid
    to H_c = 0.6. α_c is not the Womersley number (trap T2)." + `ch16.charm_kurland_core_viscosity(0.5, 310.15, 1.2e-3)`.
14. `nb.note` — **N36 [B]** "Closing the model (ours): given a, δ, H_F, μ_p, T, solve $H_c=H_F[1+Q_p/Q_c]$ with μ_c(H_c) from (16.35) — a
    monotone root find (`brentq`, P108 reminder)." + `nb.code`: `ch16.fahraeus_lindqvist_model(20e-6, 2e-6, 0.45, 1.2e-3, 310.15)` (H_c, μ_c,
    μ_app, residual < 1e-12). "⚠️ Common confusion (trap T19): the curve depends on our choice of δ, which the book leaves free."
15. `nb.plotly` — **F2** `slider_figure` on δ (0.5…5 µm, 19 steps): μ_app/μ_p against D (our closure, blue) and the cited Pries curve
    (fixed, muted). *see / read / change.*
16. `nb.md` — **What would change if…** "…the wall were not rigid? A pressure inside an elastic tube stretches it, the radius grows and,
    with R ∝ a⁻⁴, the resistance falls — C06."

#### C06 — The tube law $a=a_0\left[1-\dfrac{a_0}{Eh}(p-p_e)\right]^{-1}$ (16.41)
1. `nb.core("C06", "An artery is a spring: the tube law $a=a_0\left[1-\dfrac{a_0}{Eh}(p-p_e)\right]^{-1}$ (16.41)", question="How much does an artery widen when the pressure inside rises?")`
2. `nb.md` — **Plain words:** "Blow into a long balloon: the higher the pressure, the wider it gets, and the rubber is stretched
   around its circumference. Arteries are like that, with a stiffer rubber. The stretch around the circumference — the hoop stress
   — balances the pressure; Hooke's law turns stretch into a wider radius."
3. `nb.primer("Young's modulus and Poisson's ratio", "A solid bar stretched by a stress σ (force per area) stretches by the strain ε =
   σ/E; E is Young's modulus (Pa) — rubber ~1 MPa, an artery wall ~0.1–1 MPa, steel 200 GPa. Stretching one way thins the material
   the other way: Poisson's ratio ν̂ (0.5 for an incompressible wall, used in C11).", code="E = 0.4e6; sigma = 1.0e5   # artery wall
   [Pa], a 100 kPa stress\nprint(sigma/E)                   # strain 0.25: a 25 % stretch")` (**P379**).
4. `nb.primer("a stretched curved membrane", "A curved membrane under tension T (force per length) pushes inward with a pressure
   T/a, a the radius of curvature — Laplace's law, the same law as the surface-tension jump of Ch. 1 (Δp = σ/R for a cylinder). For a
   wall of thickness h, T = σh.", code="sigma, h, a = 1.0e5, 1.5e-3, 12e-3   # hoop stress, wall, radius\nprint(sigma*h/a)
   # 12500 Pa ≈ 94 mmHg held")` (**P380**).
5. `nb.derivation("D08", …)` — Part F D08 (6 steps), ref "16.41": holds **N37 [B]** (steps 1–3, (16.38)–(16.39)).
6. `nb.worked_example("an aorta at 100 mmHg", "a₀ = 12 mm, h = 1.5 mm, E = 0.4 MPa: 1. Eh/a₀ = 0.4 × 10⁶ × 1.5/12 = 50 kPa (375 mmHg).
   2. p − p_e = 100 mmHg = 13.3 kPa; a₀(p − p_e)/Eh = 13.3/50 = 0.267. 3. a = 12/(1 − 0.267) = **16.4 mm** (the linear estimate 12 ×
   1.267 = 15.2 mm). 4. Hoop stress (p − p_e)a/h = 13.3 × 16.4/1.5 = **145 kPa**.")`
7. `nb.code` — `ch16.tube_radius_hookean(100*MMHG, 0.0, aorta["a0"], aorta["E"], aorta["h"])`, `ch16.hoop_stress(100*MMHG, 0.0, 16.36e-3,
   aorta["h"])`; `try: ch16.tube_radius_hookean(400*MMHG, …) except ValueError as e: print(e)` (trap T16). *expect:* 0.016363 m; 1.454e5
   Pa; the singularity at 375 mmHg reported.
8. `nb.figure` — a/a₀ against p − p_e (0…350 mmHg): the tube law (teal), its linear form (muted dashed), the singular pressure Eh/a₀
   marked with a stop sign "law breaks"; hoop stress on a second axis. *see / read / change* ("…a thicker wall (larger h) pushes the
   singular pressure up: real arteries stiffen instead").
9. `nb.note` — **N38 [B]** "In a slowly widening tube each slice obeys Poiseuille's law with its own radius: $Q=-\dfrac{\pi}{8\mu}\left(\dfrac{dp}{dx}\right)a(x)^4$ (16.42), $\dfrac{dp}{dx}=-\dfrac{8\mu Q}{\pi a(x)^4}$ (16.43) —
   the lubrication idea of Ch. 8 (C05)." + `ch16.elastic_tube_pressure`.
10. `nb.note` — **N39 [B]** "Put the tube law into (16.43) and separate: $\left[1-\dfrac{a_0}{Eh}(p-p_e)\right]^{-4}dp=-\dfrac{8\mu}{\pi a_0^4}Q\,dx$ (16.44, corrected); integrating gives
    $Q=\dfrac{\pi a_0^3Eh}{24\mu L}\left\{\left[1-\dfrac{a_0}{Eh}(p_1-p_e)\right]^{-3}-\left[1-\dfrac{a_0}{Eh}(p_2-p_e)\right]^{-3}\right\}$ (16.46) (result given)." + `> ⚠️ **slip S4 — the book prints** πa₀ in (16.44)–(16.45) **; the correct form is** πa₀⁴ — a pint
    check shows the printed form has the wrong units.` + `nb.code`: `ch16.elastic_tube_flow_rate(20*MMHG, 10*MMHG, 0.0, 1e-3, 0.2e6,
    0.1e-3, 3.5e-3, 0.05)` against `ch16.poiseuille_flow_rate(1e-3, 10*MMHG, 0.05, 3.5e-3)`; the rigid limit E = 10¹² Pa. *expect:* ratio 1.53;
    rigid limit to 1e-8.
11. `nb.plotly` — **F3** `slider_figure` on E (0.05…2 MPa, 20 log steps): Q against Δp for the distensible tube (blue) and the rigid line
    (muted). *see / read / change* ("…softer walls pass more flow at high pressure: the resistance falls as the vessel inflates").
12. `nb.note` — **N40 [B]** (from the end of §16.3) "The lungs' vessels are soft: their resistance falls as the pulmonary pressure rises —
    (16.46) at work." + `nb.figure` — effective resistance (p₁ − p₂)/Q of a soft tube against its mean pressure (`ch16.elastic_tube_resistance`,
    labelled "ours, qualitative"). **N41 [C]** "The pulmonary circulation runs from the right ventricle through the lungs to the left
    atrium: larger, softer, lower-resistance vessels (the loop sketch of N02)." **N42 [C]** "Its pressure pulses look like the systemic
    ones at a lower pressure scale (C07)."
13. `nb.md` — **What would change if…** "…the pressure changed in time? A stretchy vessel then stores blood when the pressure rises and
    gives it back when it falls — the Windkessel, C07."

#### C07 — The Windkessel $Q(t)=\mathrm K\dfrac{dp}{dt}+\dfrac{p(t)}{R}$ (16.49)
1. `nb.core("C07", "Store and drain: the Windkessel $Q(t)=\mathrm K\dfrac{dp}{dt}+\dfrac{p(t)}{R}$ (16.49)", question="The heart pushes for a third of each beat — why doesn't the
   pressure drop to zero in between?")`
2. `nb.note` — **N43 [C]** the opening picture: "Old fire engines had an air chamber (German *Windkessel*) between the pump and the hose:
   the pump's strokes compressed the air, which kept water flowing between strokes. Otto Frank (1899) saw the aorta doing the same." +
   ASCII: `heart ──► [aorta: stretchy store K] ──► [small vessels: resistor R] ──► veins (p ≈ 0)`.
3. `nb.md` — **The idea:** "a capacitor in parallel with a resistor (P381): inflow = charging + leaking."
4. `nb.note` — **N44 [B]** "(16.47) in words: inflow into the aorta = rate of change of its volume + outflow into the periphery." **N45 [B]**
   "Compliance K = dv/dp; $\dfrac{dv}{dt}=\left(\dfrac{dv}{dp}\right)\left(\dfrac{dp}{dt}\right)$ (16.48). For a Hookean segment (ours): dV/dp = 2πa³L/(Eh)." + `ch16.compliance_from_tube(100*MMHG, 0.0, aorta["a0"],
   aorta["E"], aorta["h"], aorta["L"])`. *expect:* 1.38 × 10⁻⁸ m³/Pa = 1.84 mL/mmHg — close to our Windkessel K = 1.2 × 10⁻⁸.
5. `nb.derivation("D09", …)` — Part F D09 (4 steps), ref "16.49".
6. `nb.primer("a sum of decaying kicks", "For K dp/dt + p/R = Q(t) the solution adds up every past inflow Q(τ)dτ, each one raising p
   by Q dτ/K and then fading as $e^{-(t-\\tau)/RK}$ — a convolution. It extends the integrating factor of Ch. 9 (P210).", code="RK, K = 1.74, 1.2e-8
   # time constant [s], compliance [m^3/Pa]\nkick = 1e-6                            # 1 mL injected at t = 0\nfor t in (0, 1.74, 3.48):
   print(kick/K*np.exp(-t/RK))   # 83 Pa, 31 Pa, 11 Pa")` (**P382**; P210 reminder).
7. `nb.derivation("D10", …)` — Part F D10 (9 steps), ref "16.52": holds **N46 [B]**.
8. `nb.worked_example("diastolic decay", "R = 1.45 × 10⁸ Pa s m⁻³, K = 1.2 × 10⁻⁸ m³/Pa: 1. RK = 1.74 s. 2. The beat lasts 60/66 = 0.909
   s; systole 0.35 of it = 0.318 s; diastole 0.591 s. 3. Over diastole p falls by e^{−0.591/1.74} = **0.712**. 4. From 105 mmHg at the
   valve's closing to 105 × 0.712 = **75 mmHg** before the next beat — the diastolic pressure.")`
9. `nb.code` — `T = 60/heart["HR"]`; `Q = lambda t: ch16.ejection_half_sine(t, SV, T, 0.35*T)`; `t = np.linspace(0, 8*T, 4001)`; `p_ode =
   ch16.windkessel_ode(t, Q, **windkessel, p0=60*MMHG)`; `p_q = ch16.windkessel_pressure(t, Q, **windkessel, p0=60*MMHG)`; `s =
   ch16.windkessel_summary(SV=SV, HR=heart["HR"], systole_share=0.35, **windkessel)`. *expect:* p_ode and p_q agree to 1e-8; s: p0 74.92, pS
   106.03, pD 74.51, p_mean 89.73 mmHg (= R Q̄), map_third 85.02, tau 1.74 s.
10. `nb.check_agree` — **from scratch:** a 20-line RK4 loop for K dp/dt = Q − p/R over 8 beats (dt = T/400); `assert np.allclose(p_rk4, p_ode,
    rtol=1e-6)`; after 30 beats `abs(p(T) − p(0)) < 1e-6 × p(0)`; `np.isclose(ch16.time_mean(p_last, t_last), windkessel["R"]*SV/T)`.
11. `nb.figure` — Q(t) (blue) and p(t) (orange, mmHg) over three periodic beats, the K → 0 ghost p = QR muted, the true mean (amber) and
    the ⅓-rule line (dash-dot); a second panel: p_S and p_D against K (0.3…4 mL/mmHg) at fixed R — the pulse pressure falls as K grows,
    the mean does not move. *see / read / change.*
12. `nb.note` — **N48 [B]** (from §16.2, Fig. 16.5) "Along the circulation the large swings of the ventricle are damped: the ventricle's
    pressure falls near zero each beat while the aorta stays high — our curves (ventricle schematic + Windkessel aorta), labelled ours."
    + `nb.figure` (two curves on one axis).
13. `nb.animation` — **A1** (video, 80 frames; FAST 40): eight beats from p = 60 mmHg converging to the periodic state; the aortic balloon
    breathing beside the traces; the storage and run-off bars. *see / read / change* ("…a stiffer aorta settles faster: RK is shorter").
14. `nb.live` — sliders R (0.5…3 × 10⁸), K (0.3…4 mL/mmHg), HR (40…180), SV (40…120 mL) → the periodic p(t), p_S, p_D, mean and the ⅓-rule
    error (`_ = live(...)`, ch14 lesson). Paired with A1 and E2 (ch01 rule: a live widget is never the only way to see an idea).
15. `nb.note` — **N47 [C]** "The model assumes the pressure is the same everywhere in the arteries at each instant — no travel time, no
    reflection. C08 gives the pulse a finite speed; C13 adds reflections."
16. `nb.explainer("windkessel_heartbeat", heading="Why doesn't blood pressure fall to zero between beats?", why="One clock drives the
    ventricle, the aortic balloon and the pressure trace; you change K and R and watch the beat settle into a new rhythm, with the
    inflow split into storage and run-off at every instant.", tries=["Press play and watch the pressure settle within a few beats.",
    "Halve K (stiff aorta): the pulse widens, the mean stays put.", "Raise the heart rate to 120: diastole shortens and p_D climbs."])`
17. `nb.md` — **What would change if…** "…the aorta were long enough for the pressure to take time to travel along it? Then the
    pressure would be a wave, not a single number — C08."

#### C08 — The Moens–Korteweg speed $c=\sqrt{\dfrac{Eh}{2\rho a_0}}$ (16.62)
1. `nb.core("C08", "The pulse is a wave of wall stretch: $c=\sqrt{\dfrac{Eh}{2\rho a_0}}$ (16.62)", question="Is the pulse you feel at your wrist the blood arriving?")`
2. `nb.md` — **Plain words:** "Press your fingers on your wrist: the beat arrives about a tenth of a second after the heart's. The
   blood itself takes many seconds to get there. What travels is a bulge: the heart's push stretches the aorta, the stretched wall
   squeezes the blood ahead, which stretches the next piece of wall. A stiffer wall passes the push on faster."
3. `nb.md` — **The idea** (ASCII): `wall stiffness (spring) + blood inertia (mass) ──► a wave, c² = stiffness/inertia`.
4. `nb.note` — **N49 [B]** "Model: inviscid, incompressible blood in a long thin elastic tube; long waves; the pressure uniform over each
   section. Area continuity $\dfrac{\partial A}{\partial t}+\dfrac{\partial(Au)}{\partial x}=0$ (16.53) — the area is now an unknown (Ch. 15's duct laws in time)."
5. `nb.note` — **N50 [B]** "Momentum per slice: $\rho A\left(\dfrac{\partial u}{\partial t}+u\dfrac{\partial u}{\partial x}\right)=-A\dfrac{\partial p}{\partial x}$ (16.54, corrected)." + `> ⚠️ **slip S5 — the book prints** $-\dfrac{\partial((p-p_e)A)}{\partial x}$ **; the
   correct right side is** $-A\,\partial p/\partial x$ — the printed form misses the push of the sloping wall. Linearised, both give (16.57), so the
   wave speed is unaffected; the nonlinear run below uses the corrected form.` **N51 [B]** "The tube law in A: $p-p_e=\dfrac{Eh}{a_0}\left[1-\left(\dfrac{A_0}{A}\right)^{1/2}\right]$ (16.55) — the
   inverse of (16.41) with A = πa²." + `ch16.tube_law_hookean`.
6. `nb.derivation("D11", …)` — Part F D11 (7 steps), ref "16.56–16.58": holds **N52 [B]** + `nb.code`: `ch16.linearisation_ratios(25*MMHG,
   100*MMHG, 0.65, 4.86, Eh_over_a0=50e3)`. *expect:* p_ratio 0.25, u_over_c 0.13, strain 0.067 (ours).
7. `nb.derivation("D12", …)` — Part F D12 (7 steps), ref "16.62".
8. `nb.worked_example("the pulse speed in an aorta", "1. E = 0.72 MPa, h/a₀ = 0.1, ρ = 1000 kg/m³ (easy numbers): c² = 0.72 × 10⁶ × 0.1/2000 =
   36, c = **6 m/s**. 2. Our aorta (E = 0.4 MPa, h/a₀ = 0.125, ρ = 1060): c = √(0.4 × 10⁶ × 0.125/2120) = **4.86 m/s**. 3. Over 0.5 m to the
   wrist: 0.10 s. 4. Under a 25 mmHg crest the blood moves at p/(ρc) = 3333/(1060 × 4.86) = **0.65 m/s**, only while the crest
   passes.")`
9. `nb.code` — `c0 = ch16.moens_korteweg_speed(aorta["E"], aorta["h"], blood["rho"], aorta["a0"])`; `ch16.moens_korteweg_speed(...,
   midwall=True)`; `ch16.bramwell_hill_speed(np.pi*a0**2, ch16.tube_law_linear(…) slope, rho)`; `ch16.youngs_modulus_from_pwv(c0, …)`;
   `ch16.velocity_from_pulse(25*MMHG, blood["rho"], c0)`. *expect:* 4.856, 4.711 m/s; Bramwell–Hill equal to MK at A₀; E round trip 4e5 Pa;
   0.647 m/s.
10. `nb.note` — **N53 [B]** "Thick walls: evaluate the strain at mid-wall, $c=\sqrt{\dfrac{Eh}{2\rho(a_0+h/2)}}$ (16.63) (result given) — 3 % slower for h/a₀ = 0.125. Inverting
    (16.62) gives E from a measured transit time (`youngs_modulus_from_pwv`)."
11. `nb.recap("R09", "d'Alembert's solution", "$p'(x,t)=f(x-ct)+g(x+ct)$ (15.15): any two shapes sliding apart at ±c. **New here:** the book writes shapes in time,
    $p=f_1\\left(t-\\dfrac xc\\right)+f_2\\left(t+\\dfrac xc\\right)$ (16.65) — a pressure trace recorded at one place.", where="Ch. 15 §15.2 (C02, D03)")` + `ch16.pulse_wave(x, t, c0, f1)` against `AC.dalembert`.
    **(re-enter C08)**
12. `nb.recap("R10", "Pressure and particle velocity", "$u_1'=\\dfrac1{\\rho_0c}(f-g)$ (15.16). **New here:** $\\hat p_1=\\rho c\\hat u_1$, $\\hat p_2=-\\rho c\\hat u_2$ (16.68) with c the tube's
    wave speed — the seed of the admittance of C13.", where="Ch. 15 §15.2 (C02, D04)")` **(re-enter C08)**
13. `nb.check_agree` — **from scratch:** a Gaussian f; `p = f(x - c0*t)`; second differences in x and t (h = 1e-4); `assert abs(p_tt -
    c0**2*p_xx) < 1e-6*abs(p_tt).max()`; `assert np.allclose(ch16.pulse_wave(x, t, c0, ftime), p)`.
14. `nb.figure` — c against E (0.1…3 MPa, log) for h/a₀ = 0.05, 0.125, 0.25 (blue shades), the thick-wall form dashed, an illustrative
    band "young … stiff older aorta". *see / read / change* ("…c ∝ √E: a four times stiffer wall doubles the pulse speed").
15. `nb.animation` — **A2** (video, 60 frames; FAST 30): a pressure pulse travelling along an artery with its wall bulge (exaggeration
    stated, true size in the label) and tracer particles barely moving. *see / read / change.*
16. `nb.note` — **N54 [B]** "Measured pulses travel at a few metres per second in the aorta and faster in stiffer peripheral arteries;
    waves in collapsible veins are slower still (C15); real pulses peak by reflection (C13), steepen (nonlinearity) and damp (C11)."
    + `nb.figure` — two snapshots of a large pulse in `ch16.tube1d_run` (ours, DEVIATION box: MacCormack with a small artificial
    viscosity): the front steepens.
17. `nb.explainer("pulse_wave_artery", heading="Is the pulse at your wrist the blood arriving?", why="You release a pulse and watch the
    wall bulge race ahead while the blood barely moves; stiffness, thickness and radius change the speed live.", tries=["Click under
    the crest: the blood moves at 0.6 m/s, the bulge at 4.9 m/s.", "Quadruple E: the speed doubles.", "Release a bump at rest: two
    half-height pulses run apart."])`
18. `nb.md` — **What would change if…** "…the blood were viscous and the tube rigid? The 'wave' speed becomes infinite (the whole
    column moves together) and the question becomes the profile across the tube — C09."

#### C09 — Womersley's profile $u=\dfrac{Ga^2}{i\mu\alpha^2}(1-F_1)e^{i\omega t}$ (16.82)
1. `nb.core("C09", "Womersley's oscillating profile: $u(r,t)=\dfrac{Ga^2}{i\mu\alpha^2}\left(1-\dfrac{J_0(i^{3/2}\alpha r/a)}{J_0(i^{3/2}\alpha)}\right)e^{i\omega t}$ (16.82)", question="What does the velocity profile look
   like when the pressure gradient swings back and forth?")`
2. `nb.md` — **Plain words:** "Push a syringe in and out slowly and the liquid in the needle flows back and forth with a parabolic
   profile, in step with your hand. Do it fast in a wide tube and something else happens: the middle moves like a solid slug that
   cannot keep up, and only a thin layer near the wall feels the wall. Which picture applies is decided by α (C03)."
3. `nb.recap("R11", "Unsteady flow in a rigid tube", "Rigid, axisymmetric, long-wave flow with the convective term dropped: $\\rho\\dfrac{\\partial u}{\\partial t}=-\\dfrac{\\partial p}{\\partial x}+\\mu\\left[\\dfrac{\\partial^2u}{\\partial r^2}+\\dfrac1r\\dfrac{\\partial u}{\\partial r}\\right]$
   (16.71) (Ch. 4 cylindrical Navier–Stokes; Ch. 8 parallel flows). **New here:** u depends on r and t and p on x and t, so ∂p/∂x can
   depend on t only (result given). ⚠️ (trap T23) In a rigid tube this is not a wave: the whole column moves at once — waves need
   the elastic wall (C08, C11).", where="Ch. 4 §4.5; Ch. 8 (C03, C11)")` + 🔁 the cylindrical Laplacian (P186 reminder). **(re-enter C09)**
4. `nb.md` — "⚠️ Common confusion (trap T1): **two time conventions** — this block uses $e^{+i\omega t}$ (the book's (16.73)), C11 uses $e^{i(kx-\omega t)}$
   (the book's (16.116)). Two-row table: rigid tube · $\dfrac{\partial p}{\partial x}=-Ge^{i\omega t}$ (16.73) · Bessel argument $i^{3/2}\alpha$ — elastic tube · $u=\hat u(r)e^{i(kx-\omega t)}$ (16.116) · $\hat\beta=e^{i\pi/4}\alpha$. They are complex
   conjugates; never mix them in one formula."
5. `nb.note` — **N55 [B]** "One harmonic of the gradient: $\dfrac{\partial p}{\partial x}=-Ge^{i\omega t}$ (16.73) — the real part −G cos ωt is meant (complex amplitudes, P176 reminder); the
   trial $u(r,t)=U(r)e^{i\omega t}$ (16.75)."
6. `nb.primer("Bessel's equation and J₀, Y₀, J₁", "The equation $y''+y'/z+y=0$ is what the Laplacian in a circle turns into for
   oscillating problems. Its two solutions are $J_0(z)$, smooth at z = 0 (= 1 there), and $Y_0(z)$, infinite at z = 0. For real z, J₀
   wiggles and decays like a damped cosine; for a complex z such as $i^{3/2}\\alpha$ it grows exponentially. Its power series is
   $J_0(z)=\\sum_m(-1)^m(z/2)^{2m}/(m!)^2$; $J_1$ is its partner, $J_0'=-J_1$.", code="from scipy.special import jv, yv\nprint(jv(0, [0, 2.405, 5.0]),
   yv(0, 1e-3))     # 1, ≈0 (first zero), −0.178; Y0 is −4.47 near 0")` (**P374**).
7. `nb.primer("powers of i", "$i^{1/2}=e^{i\\pi/4}=(1+i)/\\sqrt2$ and $i^{3/2}=e^{3i\\pi/4}=(-1+i)/\\sqrt2$; also $i\\cdot i^{1/2}=i^{3/2}$ and $(i^{3/2})^2=i^3=-i$. numpy's `1j**1.5`
   takes the principal branch and gives the same number (P159 and P45 reminders).", code="print(1j**0.5, 1j**1.5, np.exp(3j*np.pi/4))
   # (0.707+0.707j) (−0.707+0.707j) twice\nprint((1j**1.5)**2)                                  # −1j")` (**P377**).
8. `nb.primer("scipy.special.jv with a complex argument", "`jv(n, z)` evaluates J_n for complex z. For large |Im z| the values overflow
   (J₀(i^{3/2}α) grows like $e^{\\alpha/\\sqrt2}$), so beyond α ≈ 30 use `jve(n, z)` = $J_n(z)e^{-|\\mathrm{Im}\\,z|}$ and take ratios (trap T25).",
   code="from scipy.special import jv, jve\nz = np.exp(3j*np.pi/4)*40.0\nprint(abs(jv(0, z)), abs(jve(0, z)))     # 3.8e11 and 0.0627: same
   number, scaled")` (**P375**).
9. `nb.derivation("D13", …)` — Part F D13 (13 steps), ref "16.82", with its `check_src`: holds **N56 [B]** (step 4, with `> ⚠️ **slip S6 — the book
   prints** $+\dfrac G\mu$ on the right of (16.76) **; the correct right side is** $-\dfrac G\mu$ — its own solution (16.77) needs the minus.`) + `nb.code`:
   `ch16.womersley_ode_residual(0.5e-3, 1e-3, 2*np.pi, 1060.0, 3.5e-3, 1.0)`. *expect:* −G/μ = −285.7 (to 1e-9 relative).
10. `nb.worked_example("the centreline at α = 3", "1. z = i^{3/2} × 3: z² = −9i. 2. J₀(z) from the series: 1 + 9i/4 − 81/64 − 729i/2304 + … =
    −0.2214 + 1.9376i (the series converges fast). 3. On the axis J₀(0) = 1, so F₁ = 1/J₀(z) = −0.0582 − 0.5095i. 4. 1 − F₁ = 1.0582 +
    0.5095i. 5. u/(Ga²/4μ) = (4/(9i))(1.0582 + 0.5095i) = 0.2264 − 0.4703i: amplitude **0.522** of Poiseuille's centre speed, lagging by
    arctan(0.4703/0.2264) = **64.3°**.")`
11. `nb.code` — `ch16.womersley_F1(r, 1.0, 3.0)`; `ch16.womersley_velocity(r, t, artery["a"], omega, blood["rho"], blood["mu"], 100.0)`;
    `ch16.womersley_velocity_nd(0.0, 0.0, 3.0)`. *expect:* 0.2264 (real part at ωt = 0; amplitude 0.522, lag 64.3°).
12. `nb.check_agree` — **from scratch:** `J0 = lambda z: sum((-1)**m*(z/2)**(2*m)/math.factorial(m)**2 for m in range(60))` (complex z);
    `u_mine = (4/(1j*al**2))*(1 - J0(i32*al*s)/J0(i32*al))`; `assert np.allclose(u_mine, …womersley_velocity_nd complex…, rtol=1e-10)` for
    α = 0.5, 3, 8; and `ch16.womersley_pde_solve(...)` (Crank–Nicolson, cached in FAST — P252) matches the Bessel profile to its grid
    error (printed).
13. `nb.derivation("D14", …)` — Part F D14 (6 steps), ref "16.87".
14. `nb.derivation("D15", …)` — Part F D15 (8 steps), ref "16.88": holds **N58 [B]** with S7 noted (nω for harmonic n). 🔁 before it: "Stokes'
    second problem (Ch. 8 C11): $u(y,t)=U\exp\{-y\sqrt{\omega/2\nu}\}\cos(\omega t-y\sqrt{\omega/2\nu})$ (8.38); `LAM.stokes_second_problem`."
15. `nb.note` — **N57 [B]** "A real pressure gradient is a sum of harmonics: $\dfrac{\partial p}{\partial x}=-\sum_{n=0}^NG_ne^{in\omega t}$ (16.83); because the equation is linear the velocity is the
    sum of their profiles, $u_n(r)=\dfrac{G_na^2}{i\mu\alpha_n^2}(1-F_1(\alpha_n))$ (16.86) with $\alpha_n=\alpha\sqrt n$, and the mean (n = 0) is Poiseuille's $u_0(r)=\dfrac{G_0a^2}{4\mu}\left(1-\dfrac{r^2}{a^2}\right)$ (16.85) (result given)." +
    `nb.code`: a synthetic gradient of four harmonics (ours) through `ch16.womersley_harmonics` and its profiles over a beat.
16. `nb.figure` — **N59 [B]** our profiles over half a cycle at α = 1, 4, 8 and 14 (our own values), every 15° of phase, across the full
    diameter (⚠️ trap T20: the book's figure runs wall–centre–wall mirrored; ours shows the whole diameter). *see:* "…at small α parabolas
    breathing in step; at large α flat cores lagging…" · *read:* "…reversal appears first near the wall…" · *change:* "…double the heart
    rate and every α grows by √2." **N60 [C]** "Read twice in the book: reversal near the wall, a flat core, an almost solid core at the
    largest α."
17. `nb.animation` — **A3** (video, 72 frames = one cycle; FAST 36): profiles at α = 1, 4 and 12 side by side, the gradient's phase as a
    clock hand. *see / read / change.*
18. `nb.explainer("womersley_profiles", heading="Is pulsing flow just Poiseuille flow sloshing?", why="Dragging α morphs the profile from a
    breathing parabola to a lagging plug while the phase runs; the two limits snap on as ghosts and the vessel axis shows where each
    vessel sits.", tries=["Drag α from 0.5 to 17 and watch the profile flatten.", "Watch the purple near-wall trace turn before the
    centre.", "Click a radius to see the Bessel arithmetic."])`
19. `nb.md` — **What would change if…** "…we asked for the total flow instead of the profile? Integrate (16.82) over the section:
    C10."

#### C10 — Amplitude (16.99) and phase lag (16.100) of pulsatile flow
1. `nb.core("C10", "How much flow, how late: $|Q|=\dfrac{\pi Ga^4}{\mu\alpha^2}\sqrt{X^2+Y^2}$ (16.99), $\tan\phi=\dfrac XY$ (16.100)", question="Does the flow follow the pressure gradient?")`
2. `nb.md` — **Plain words:** "Push a child on a swing gently and often: the swing moves with your hand. Push a heavy cart back and
   forth quickly: it barely moves, and its motion trails your push by a quarter cycle — your force goes into accelerating it. Blood
   in a large artery is the cart."
3. `nb.primer("two Bessel identities", "$\\dfrac{d}{dz}J_0(z)=-J_1(z)$ and $\\displaystyle\\int_0^zsJ_0(s)\\,ds=zJ_1(z)$. The first gives the wall shear, the second the
   flow rate. Both hold for complex z.", code="from scipy.special import jv\nfrom scipy.integrate import quad\nz = 2.0; h = 1e-6\nprint((jv(0,
   z + h) - jv(0, z - h))/(2*h), -jv(1, z))     # −0.5767 twice\nprint(quad(lambda s: s*jv(0, s), 0, z)[0], z*jv(1, z))   # 1.1534 twice")` (**P376**).
4. `nb.derivation("D16", …)` — Part F D16 (7 steps), ref "16.90–16.91": holds **N61 [B]**.
5. `nb.derivation("D17", …)` — Part F D17 (7 steps), ref "16.99–16.100".
6. `nb.worked_example("α = 8", "1. F₂(8) = 0.1772 − 0.1608i (`womersley_F2`). 2. Z = 1 − F₂ = 0.8228 + 0.1608i: X = 0.8228, Y = 0.1608. 3. √(X² + Y²) =
   0.8384. 4. |Q|/|Q₀| = 8 × 0.8384/64 = **0.105** — a tenth of the steady flow for the same gradient amplitude. 5. φ = arctan2(0.8228,
   0.1608) = **78.9°**; at 66 beats per minute (ω = 6.91 rad/s) the flow peaks 1.378/6.91 = **0.199 s** after the gradient. (The book's own worked example
   does not quite match (16.91) — S11; the comparison stays private; ours is computed.)")`
7. `nb.code` — `ch16.womersley_F2(8.0)`, `ch16.womersley_XY(8.0)`, `ch16.womersley_flow_ratio(8.0)`, `np.degrees(ch16.womersley_phase_lag(8.0))`,
   `ch16.womersley_flow_rate(t, artery["a"], omega, blood["rho"], blood["mu"], 100.0)`. *expect:* 0.1048, 78.94°.
8. `nb.check_agree` — **from scratch:** at one instant `Q_mine = np.trapezoid(2*np.pi*r*ch16.womersley_velocity(r, t0, …), r)` (2001 points);
   `assert np.isclose(Q_mine, ch16.womersley_flow_rate(t0, …), rtol=1e-8)`.
9. `nb.note` — **N62 [B]** "Wall shear: $\tau_w=-\mu\dfrac{\partial u}{\partial r}\Big\rvert_a=\dfrac{G_0a}2+\dfrac a2\sum_1^NG_nF_2(\alpha_n)e^{in\omega t}$ (16.92, corrected) (result given, from $dJ_0/dz=-J_1$)." + `> ⚠️ **slip S8 — the book prints** this as a
   shear rate ∂u/∂r with an unnamed F **; it is the wall shear stress −μ∂u/∂r with F₂**.` + `nb.code`: a finite-difference wall shear
   from `womersley_velocity` against `ch16.womersley_wall_shear(t, a, omega, rho, mu, [0, G])` (1e-8). "The wall shear leads the flow."
10. `nb.note` — **N63 [B]** "Small α: $F_2(\alpha)\approx1-\dfrac{i\alpha^2}8$ (16.93), so the flow tends to Poiseuille's $|Q_0|=\dfrac{\pi Ga^4}{8\mu}$ (16.95) (the book writes $e^{in\omega t}$ in (16.94) for $e^{i\omega t}$ —
    S9)." + `nb.code`: errors of `womersley_F2_series` at α = 0.25, 0.5, 1 and their slope (expect ≈ 4). **N64 [B]** "Large α:
    $F_2\approx\dfrac{2}{i^{1/2}\alpha}\left(1+\dfrac{i^{3/2}}{2\alpha}\right)$ (16.96, corrected)." + `> ⚠️ **slip S10 — the book prints** $\left(1+\dfrac1{2\alpha}\right)$ **; the correct factor is** $\left(1+\dfrac{i^{3/2}}{2\alpha}\right)$.` + a three-row table
    (α = 10, 30, 100: error of the leading term, of the corrected, of the printed form — expect 1.0e-2, 2.7e-4, 1.9e-2 at α = 10; the
    printed correction is worse than none).
11. `nb.figure` — Bode-style: |Q|/|Q₀| (log) and φ (degrees) against α (log, 0.1…50), asymptotes 1, 8/α², 0°, 90° as muted guides, the α
    of our artery and aorta marked. *see / read / change* ("…a damper below α ≈ 1, a mass above α ≈ 10").
12. `nb.explainer("womersley_flow_phase", heading="Does the flow follow the pressure?", why="Dragging α moves a cursor on the amplitude
    and lag curves while the time traces slide apart and the complex number Z turns.", tries=["Set α = 1: flow and gradient peak
    together.", "Set α = 20: the flow is 2 % of Poiseuille's and 86° late.", "Switch the asymptote to the printed (16.96) at α = 10 and
    compare the errors."])`
13. `nb.md` — **What would change if…** "…the tube were elastic as well? Then the pressure travels (C08) and the viscous profile of C09
    rides on it: C11."

#### C11 — Viscous waves in an elastic tube: the quadratic (16.136)
1. `nb.core("C11", "Where waves meet viscosity: the characteristic quadratic $X^2\left[2\hat\beta\dfrac{J_0(\hat\beta)}{J_1(\hat\beta)}-4\right]+X\left[4\hat\nu-1-2\hat\beta\dfrac{J_0(\hat\beta)}{J_1(\hat\beta)}\right]+(1-\hat\nu^2)=0$ (16.136)", question="How much do viscosity
   and the vessel wall slow the pulse and make it fade?")`
2. `nb.md` — **Plain words:** "C08 gave the pulse a speed from the wall's stiffness, ignoring viscosity. C09 gave the oscillating
   profile, ignoring the wall's give. A real artery has both: the pulse drags a viscous layer along the wall, which slows it and
   makes it fade with distance. How much depends — once more — on α."
3. `nb.md` — "⚠️ Common confusion (trap T1), repeated: this block uses $e^{i(kx-\omega t)}$ (16.116), the opposite sign to C09. The two-row
   table of C09 again."
4. `nb.recap("R12", "Continuity and Navier–Stokes in (r, x)", "$\\dfrac{\\partial u}{\\partial x}+\\dfrac1r\\dfrac{\\partial(rv)}{\\partial r}=0$ (16.101) and $\\rho\\dfrac{\\partial u}{\\partial t}=-\\dfrac{\\partial p}{\\partial x}+\\mu\\left(\\dfrac{\\partial^2u}{\\partial r^2}+\\dfrac1r\\dfrac{\\partial u}{\\partial r}+\\dfrac{\\partial^2u}{\\partial x^2}\\right)$ (16.102), linearised. **New here:** the radial
   velocity v is kept, because the wall moves.", where="Ch. 4 §4.5; Appendix B")` **(re-enter C11)**
5. `nb.note` — **N65 [B]** "The wall is a membrane with displacements η (radial) and ξ (axial); its tensions follow Hooke's law with
   Poisson's ratio: $T_\theta=\dfrac{Eh}{1-\hat\nu^2}\left(\dfrac\eta{a_0}+\hat\nu\dfrac{\partial\xi}{\partial x}\right)$ (16.104); its radial motion $\rho_wh\dfrac{\partial^2\eta}{\partial t^2}=\sigma_{rr}\rvert_a-\dfrac{T_\theta}{a_0}$ (16.106) — the hoop tension's
   inward push T/a (P380) (stated with a picture; a-d25)." **N66 [B]** "With the fluid's stresses $\sigma_{rr}=+p-2\mu\dfrac{\partial v}{\partial r}$ (16.110) and $\sigma_{rx}=\mu\left(\dfrac{\partial u}{\partial r}+\dfrac{\partial v}{\partial x}\right)$ (16.111)
   put in, the wall equations become (16.112)–(16.113) (stated). The book's (16.108) prints 1 − ν̂ for 1 − ν̂² (S12). This chapter counts
   a normal stress positive in compression here." **N67 [B]** "No slip at the moving wall, applied at the undisturbed radius:
   $u\rvert_{a_0}=\dfrac{\partial\xi}{\partial t}$ (16.114), $v\rvert_{a_0}=\dfrac{\partial\eta}{\partial t}$ (16.115) — Ch. 7's linearised boundary condition."
6. `nb.primer("a complex wavenumber", "With $k=k_1+ik_2$, the wave $e^{i(kx-\\omega t)}=e^{-k_2x}e^{i(k_1x-\\omega t)}$ travels at $c=\\omega/k_1$ and its amplitude
   falls by $e^{-2\\pi k_2/k_1}$ over one wavelength — the spatial twin of the complex frequency of Ch. 11 (growth in time).", code="k1, k2 =
   1.0, 0.0622                      # per metre (α = 8 below)\nprint(np.exp(-2*np.pi*k2/k1))           # 0.677: a third lost per
   wavelength")` (**P385**).
7. `nb.derivation("D18", …)` — Part F D18 (9 steps), ref "16.116–16.125": holds **N68 [B]** (S13 noted) and **N69 [B]**.
8. `nb.derivation("D19", …)` — Part F D19 (10 steps), ref "16.128–16.129", with its `check_src`: holds **N70 [B]** (S14).
9. `nb.derivation("D20", …)` — Part F D20 (10 steps), ref "16.131–16.133", with its `check_src`: holds **N71 [B]**.
10. `nb.primer("choosing a root by continuation", "A quadratic has two roots; which one is 'the' wave? Start where the answer is known
    (here α → ∞: X → 1, the inviscid wave of C08), then lower α in small steps, each time picking the root nearest the previous one.
    It extends the complex quadratic formula of Ch. 6 (P159).", code="def roots(b, c): d = np.sqrt(b*b - 4*c + 0j); return (-b + d)/2, (-b -
    d)/2\nprev = 1.0\nfor eps in (0.0, 0.1, 0.2):                       # a family X² − (1 + eps(1+1j))X + 0 = 0\n    r = roots(-(1 +
    eps*(1 + 1j)), 0.0); prev = min(r, key=lambda z: abs(z - prev)); print(prev)")` (**P386**; P330, P159 reminders).
11. `nb.derivation("D21", …)` — Part F D21 (13 steps), ref "16.136", with its `check_src`.
12. `nb.worked_example("α = 8, ν̂ = 0.5", "1. β̂ = e^{iπ/4} × 8 = 5.657 + 5.657i; g = β̂J₀(β̂)/J₁(β̂) = 6.19 − 5.62i (`jv`). 2. (16.136) has two roots; the one near 1 is X =
    1.094 + 0.137i. 3. √X = 1.048 + 0.065i: c/c₀ = 1/1.048 = **0.954** — 5 % slower than Moens–Korteweg. 4. k₂/k₁ = 0.065/1.048 = **0.062**:
    the amplitude falls to e^{−2π × 0.062} = **0.68** in one wavelength. 5. In the aorta-like α = 17 the pulse would lose only a fifth
    per wavelength; at α = 3 it loses four fifths.")`
13. `nb.code` — `ch16.elastic_tube_roots(8.0, 0.5)`, `ch16.elastic_tube_wave(8.0, c0=1.0, poisson=0.5)`, the same at α = 3, 20; the second root
    labelled. *expect:* c_over_c0 0.954, k2_over_k1 0.0622, damping_per_wavelength 0.677; α = 3: 0.910, 0.252; α = 20: 0.981, 0.0218;
    other root at α = 8: c/c₀ 4.93.
14. `nb.check_agree` — **independent route:** `ch16.elastic_tube_determinant_numeric(X_root, 8.0, 0.5)` (the mode ODEs integrated without
    Bessel closed forms) is ≈ 0 (|·| < 1e-8 of its scale) at the closed-form root, and not at X = 1.
15. `nb.figure` — c/c₀ and k₂/k₁ against α (1…50, log), the pressure root bold, the other root faint (labelled "excluded: fast, heavily
    damped"); 1 − c/c₀ ∝ 1/α marked. *see / read / change* ("…both corrections fade as 1/α: the aorta sees nearly the inviscid wave").
16. `nb.plotly` — **F4** `slider_figure` on ν̂ (0.3…0.5, 9 steps): c/c₀ and k₂/k₁ against α for both roots.
17. `nb.animation` — **A5** (video, 60 frames; FAST 30): the pressure wave along a tube for α = 3 and α = 15 against the inviscid wave of
    the same frequency (muted). *see / read / change.*
18. `nb.note` — **N72 [B]** "A viscoelastic wall: replace E by $E^*=E-i\omega E'$, $\hat\nu^*=\hat\nu-i\omega\hat\nu'$ (16.137) — Kelvin–Voigt in the $e^{-i\omega t}$ convention (the sign of the imaginary part
    belongs to the convention, T1)." + `ch16.elastic_tube_wave(8.0, E_ratio=1 - 0.1j)`: more damping. **N73 [C]** "A steady stream adds to
    the wave speed downstream and subtracts upstream (Ch. 15's U ± c, `AC.wave_speeds`); small in arteries."
19. `nb.md` — **What would change if…** "…the tube branched? A pulse meeting a fork is partly reflected — C13. But first: why are the
    branches the sizes they are? C12."

#### C12 — Murray's law $a=\left[\dfrac{16\mu}{\pi^2K}\right]^{1/6}Q^{1/3}$ (16.140, corrected)
1. `nb.core("C12", "Murray's law: $a=\left[\dfrac{16\mu}{\pi^2K}\right]^{1/6}Q^{1/3}$ (16.140, corrected) and $a_0^3=a_1^3+a_2^3$ (16.144)", question="Why are blood vessels the width they are?")`
2. `nb.md` — **Plain words:** "A wider pipe is easier to push water through — so why not make every vessel wide? Because blood is
   expensive: the body must make it, warm it and carry it around. Murray asked which radius makes the total cost — pumping plus upkeep
   — smallest. The answer is a rule you can test with a ruler on any cast of a vascular tree."
3. `nb.md` — **The idea** (ASCII): `cost ↑  pumping ∝ a⁻⁴ ╲      ╱ upkeep ∝ a²` / `        ╲_minimum_╱` / `        ──────────────► radius a`
   (P349 reminder: the minimum of a falling plus a rising power, ch14).
4. `nb.note` — **N74 [B]** "The cost of a segment: pumping power Δp·Q = RQ² with R from (16.20), plus an upkeep power K per unit volume of
   blood (W/m³ — the book's K, not the compliance; trap T2): $\hat E_1=\dfrac{8\mu L}{\pi a^4}Q^2+K\pi a^2L$ (16.138)." + `ch16.murray_cost`.
5. `nb.derivation("D22", …)` — Part F D22 (8 steps), ref "16.140" + `> ⚠️ **slip S15 — the book prints** $a=\left[\dfrac{16\mu}{\pi^2}K\right]^{1/6}Q^{1/3}$ (16.140) and the same bracket in
   (16.143) **; the correct form is** $a=\left[\dfrac{16\mu}{\pi^2K}\right]^{1/6}Q^{1/3}$ — a units check catches it: μ/K has units s², so (μQ²/K)^{1/6} is metres.` + `nb.code`: `UN.check_dimensions(...)`
   on both forms.
6. `nb.worked_example("an aorta-sized optimum", "Q = 82.5 mL/s, μ = 3.5 mPa s, K = 1000 W/m³: 1. 16μ/(π²K) = 0.056/9870 = 5.67 × 10⁻⁶ s² (μ/K is
   (Pa s)/(Pa/s)); its sixth root 0.1335 s^{1/3}. 2. Q^{1/3} = (8.25 × 10⁻⁵)^{1/3} = 0.04355 m s^{−1/3}. 3. a = 0.1335 × 0.04355 = **5.81 mm** (ours). 4. Over 10
   cm: pumping 5.31 mW, upkeep 10.62 mW — exactly 1 : 2. 5. The printed formula would give 58.1 mm.")`
7. `nb.code` — `a_star = ch16.murray_radius(CO, murray["mu"], murray["b_metabolic"])`; `ch16.murray_terms(a_star, CO, murray["L"], murray["mu"],
   murray["b_metabolic"])`; `ch16.murray_min_cost(...)`; `ch16.murray_radius(..., printed_16_140=True)`. *expect:* 5.814e-3 m; pumping 5.309e-3 W,
   upkeep 1.0619e-2 W; min 1.5928e-2 W; printed 5.814e-2 m.
8. `nb.check_agree` — **from scratch:** a grid search `a = np.linspace(1e-3, 20e-3, 20001)`; `cost = 8*mu*L*Q**2/(np.pi*a**4) + b*np.pi*a**2*L`;
   `assert abs(a[np.argmin(cost)] - a_star) < 2e-6`; and the printed radius gives a cost at least 10× the minimum (assert).
9. `nb.figure` — pumping (orange), upkeep (teal), total (amber) against a (log), optimum marked; a second panel: the symmetric fork drawn to
   scale with 37.47° arcs. *see / read / change.*
10. `nb.plotly` — **F5** `slider_figure` on Q (1…200 mL/s, 25 log steps): the three costs against a; the optimum slides as Q^{1/3}.
11. `nb.derivation("D23", …)` — Part F D23 (3 steps), ref "16.144": holds **N75 [B]** + `ch16.murray_daughters(a_star, 0.5)`. *expect:* 0.7937 a*.
12. `nb.primer("scipy.optimize.minimize for a point in the plane", "`minimize(f, x0, method='Nelder-Mead')` finds a minimum of a function
    of several variables by moving a small triangle downhill; no derivatives needed. Here it finds the cheapest branch point B.
    `minimize_scalar` (P170) did the one-variable case.", code="from scipy.optimize import minimize\nf = lambda B: (B[0] - 1)**2 + (B[1] +
    2)**2                  # a bowl centred at (1, −2)\nprint(minimize(f, [0, 0], method='Nelder-Mead').x)        # ≈ [1, −2]")` (**P387**).
13. `nb.primer("the weighted Fermat point", "Move a point B by a small step ds in a direction m: its distance to a fixed point P changes by
    −(unit vector from B to P)·m ds (a directional derivative, P75). So the weighted sum Σwᵢ|Pᵢ − B| is least where the weighted unit
    vectors wᵢeᵢ add to zero — they close a triangle of 'forces', whose angles the cosine rule (P340) gives.", code="w = np.array([1.0,
    0.63, 0.63])                        # weights a_i² of a symmetric Murray fork\nprint(np.degrees(np.arccos((w[0]**2 + w[1]**2 - w[2]**2)/(2*w[0]*w[1]))))   # 37.47°")`
    (**P388**; P75, P340 reminders).
14. `nb.derivation("D24", …)` — Part F D24 (11 steps), ref "16.146–16.147", with its `check_src`: holds **N76 [B]** + `nb.code`: `ch16.murray_angles(…)`,
    `ch16.murray_theta(…)`, `np.degrees(ch16.murray_symmetric_angle())`, `ch16.bifurcation_point(A, C, D, …)`. *expect:* 37.467°; 70/30 split θ 25.71°,
    φ 49.73°.
15. `nb.explainer("murray_branching", heading="Why isn't a wider vessel always better?", why="You drag the radius past the minimum and
    watch the two cost bars trade, then drag the branch point and see the cost bottom out exactly where three weighted pulls close a
    triangle.", tries=["Put the radius cursor at the minimum: pumping is half the upkeep.", "Pick the 70/30 split: the thin daughter turns
    more.", "Switch to the printed (16.140) and raise K."])`
16. `nb.md` — **What would change if…** "…a pressure pulse arrived at this fork? Its daughters, Murray-sized, are almost matched to the
    parent — C13 says what 'matched' means."

#### C13 — Reflection at a junction (16.153)
1. `nb.core("C13", "Echoes at a fork: $\mathrm R=\dfrac{Y_1-\sum Y_j}{Y_1+\sum Y_j}$, $\mathrm T=\dfrac{2Y_1}{Y_1+\sum Y_j}$ (16.153)", question="When a pulse reaches a branch, how much comes back?")`
2. `nb.md` — **Plain words:** "Send a wave along a rope tied to a thinner rope: part of it bounces back. Tie it to a thicker one and the
   echo comes back upside down. At an arterial fork the 'thickness' that matters is not the radius alone but how easily each tube
   accepts flow for a given pressure — its admittance."
3. `nb.note` — **N77 [B]** "At the junction the pressure is single-valued and flux is conserved: $p_I+p_R=p_{T1}=p_{T2}$ (16.148), $Q_I-Q_R=Q_{T1}+Q_{T2}$ (16.149) (the book's minus with a
   negative Q_R counts the sign twice — S16; we keep signs once, D25)." **N78 [B]** "Multiply $\hat p=\rho c\hat u$ (16.68) by A/(ρc): $Q=Au=\pm Yp=\pm\dfrac pZ$ (16.151), admittance
   Y = A/(ρc), impedance Z = ρc/A — the electrical analogy (P381)." + `ch16.tube_admittance`, `ch16.tube_impedance`.
4. `nb.derivation("D25", …)` — Part F D25 (6 steps), ref "16.153".
5. `nb.derivation("D26", …)` — Part F D26 (6 steps), ref "16.163": holds **N80 [B]**.
6. `nb.worked_example("a symmetric fork", "Daughters of radius 0.8a, same wave speed: 1. Y_j/Y₁ = (A_j/A₁)(c₁/c_j) = 0.64 each; ΣY/Y₁ = 1.28. 2. R = (1 −
   1.28)/(1 + 1.28) = **−0.123**: a small inverted echo. 3. T = 1 + R = 0.877. 4. Energy: R² = 0.015 reflected, 1.28 × 0.877² = 0.985
   transmitted; sum 1. 5. A Murray fork (daughters 2^{−1/3}a): ΣY/Y₁ = 1.26, R = **−0.115** — nearly matched.")`
7. `nb.code` — `ch16.junction_coefficients(1.0, [0.64, 0.64])`, `ch16.junction_energy(1.0, [0.64, 0.64])`, `ch16.junction_from_geometry(a, [0.8*a,
   0.8*a], c0, [c0, c0], blood["rho"])`. *expect:* (−0.1228, 0.8772); total 1.0.
8. `nb.check_agree` — **from scratch:** `np.linalg.solve([[1, -1], [Y1, sum(Yd)]], [-1, Y1])` for (p_R, p_T) per unit p_I (the two junction
   conditions); `assert np.allclose(sol, [R, T])`; energy identity for 1000 random admittances (`np.random.default_rng(0)`) to 1e-14.
9. `nb.figure` — R and T against ΣY/Y₁ (0…4): R from +1 to −1 through 0 at 1, T from 2 to 0; the three regimes shaded and labelled (closed-end
   like, matched, open-end like); our Murray fork marked. *see / read / change.*
10. `nb.note` — **N79 [B]** "In the parent the pressure and flux are $\dfrac p{P_I}=f\left(t-\dfrac x{c_1}\right)+\dfrac{P_R}{P_I}f\left(t+\dfrac x{c_1}\right)$ (16.159) and $Q=Y_1P_I\left[f\left(t-\dfrac x{c_1}\right)-\dfrac{P_R}{P_I}f\left(t+\dfrac x{c_1}\right)\right]$ (16.160): once the echo returns, p and Q no longer share a
    shape. (The book writes the reflected shape as f; it equals f here only because R does not depend on frequency — S17.)" +
    `ch16.junction_waveforms` probe traces.
11. `nb.note` — **N81 [B]** three-row table (trap T18): ΣY_j < Y₁: R > 0, 'closed-end-like' — pressure rises at the fork (a closed end, ΣY = 0,
    doubles it: R = +1) · ΣY_j = Y₁: matched · ΣY_j > Y₁: R < 0, 'open-end-like'. "Reflections from the leg arteries help make the aortic
    pressure peak; repeated forks make a tree."
12. `nb.animation` — **A4** (frames, 40; stops at the junction and at the echo's return; FAST 24): a pulse down a parent into a symmetric
    fork (`ch16.network_pulse_run`, leapfrog at C = 1 — cross-checks (16.153) to 1e-3), with the probe traces of p and Q. *see / read / change.*
13. `nb.explainer("junction_reflection", heading="When does a fork send the pulse back?", why="Changing the daughters' radius or stiffness
    you watch the echo shrink, vanish and flip while the energy bars always add to the incident pulse.", tries=["Find the radius ratio
    with no echo.", "Make the daughters narrow: the echo comes back upright.", "Watch the probe: p and Q part ways when the echo arrives."])`
14. `nb.md` — **What would change if…** "…the vessel curved instead of branching? A bend does not reflect a long pulse much, but it
    reshapes the steady flow inside it — C14."

#### C14 — Dean's expansion (16.177)
1. `nb.core("C14", "A bend makes two vortices: Dean's expansion $w=\sum\kappa^{2n}w_n,\ \psi=\kappa\sum\kappa^{2n}\psi_n$ (16.177)", question="What does a bend do to the flow inside a pipe?")`
2. `nb.note` — **N82 [C]** the scene: "The aortic arch turns the blood through half a circle (and out of plane). Curvature distorts the
   profile, sets up secondary motions, raises dissipation and delays turbulence; Ch. 11 (N69) named Dean flow among the centrifugal
   instabilities."
3. `nb.md` — **Plain words:** "Stir tea and the leaves gather in the middle of the bottom — a secondary flow. In a bent pipe the fast
   fluid at the centre is flung toward the outside of the bend harder than the slow fluid near the walls; it piles up at the outer
   wall and returns along the walls: two counter-rotating vortices, and the fastest axial flow sits off-centre, toward the outer wall."
4. `nb.primer("toroidal coordinates of a curved pipe", "Inside a pipe bent into a circle of radius R: r and α locate a point in the
   cross-section (α = 0 toward the outside of the bend here), s runs along the centreline. The axial speed w is carried round the bend,
   so it feels a centrifugal acceleration w²/R pointing outward — that is the term $-w^2\\cos\\alpha$ of the radial equation (16.167).",
   code="w, R = 0.3, 0.05                     # axial speed [m/s], bend radius [m]\nprint(w**2/R)                        # 1.8 m/s²: a fifth of
   g")` (**P389**; trap T2: this α is an angle).
5. `nb.note` — **N83 [B]** "The secondary flow forms Dean vortices; for a developed flow the axial maximum skews to the outer wall (for a
   potential-vortex entrance flow, to the inner wall); heart–lung machines use coiled tubes for mixing." + our cross-section sketch.
   **N84 [B]** "Curvature ratio $\delta=\dfrac aR$ (16.164) — not the cell-free layer δ of C05 (T2)." + `ch16.curvature_ratio(12.5e-3, 50e-3)`. *expect:*
   0.25 (our illustrative arch).
6. `nb.note` — **N85 [B]** "The book rescales so that centrifugal, inertial and viscous terms are the same size (√δ scaling) and splits
   the pressure $p=p_0(s)+\delta p_1(r,\alpha,s)+\dots$ (16.165) (stated)." **N86 [B]** "The resulting equations (16.166)–(16.170) (stated, a-d35): the one term that
   matters is coloured — in the radial equation the centrifugal push $-w^2\cos\alpha$." **N88 [B]** "A stream function for the secondary flow: $u=\dfrac1r\dfrac{\partial\psi}{\partial\alpha}$, $v=-\dfrac{\partial\psi}{\partial r}$
   (16.172) (Ch. 4's stream function), and the axial equation with p₀ = −Gs (16.173)." **N89 [B]** "The vorticity equation (16.174) with
   $\nabla_1^2=\dfrac{\partial^2}{\partial r^2}+\dfrac1r\dfrac{\partial}{\partial r}+\dfrac1{r^2}\dfrac{\partial^2}{\partial\alpha^2}$ (16.175) and $\psi=\dfrac{\partial\psi}{\partial r}=w=0$ at $r=1$ (16.176). Read its right side in words: the axial speed's own gradient drives the secondary
   vorticity."
7. `nb.note` — **N87 [B]** "The Dean number $\kappa=\sqrt\delta\,\dfrac{2a\bar W_0}\nu=\sqrt\delta\,2Re$ (16.171)." + a **'which Dean number?' table** (trap T8): κ (16.171) · K = 2κ² (16.178) · D = 4√2 κ (16.179), K = D²/16
   · the literature's De = Re√δ with the diameter Re (= κ). "With the chapter's diameter Re, (16.171) reads √δ Re — S18." + `ch16.dean_number(0.25,
   5e-3, 12.5e-3, 3.3e-6, convention="kappa")`, `ch16.dean_convert(576.0, from_="K", to="D")`. *expect:* D = 96.
8. `nb.primer("a regular perturbation series", "When a small parameter ε appears in an equation, write the answer as $f=f_0+\\varepsilon
   f_1+\\varepsilon^2f_2+\\dots$, put it in, and collect equal powers of ε: each power gives a simpler equation, solved in turn. Ch. 2's
   orders of smallness (P68) only dropped terms; here we keep and solve them.", code="import sympy as sp\neps, x = sp.symbols('epsilon x')\nf
   = 1 + eps*x                                     # a guess to first order for x² − (2 + eps)x + 1 + eps = 0 near x = 1\nprint(sp.series(sp.solve(x**2 - (2 + eps)*x + 1 + eps, x)[1], eps, 0, 2))")`
   (**P391**).
9. `nb.primer("the biharmonic operator in polar coordinates", "$\\nabla^4=\\nabla^2\\nabla^2$. On $r^n\\sin\\alpha$ the polar Laplacian gives $(n^2-1)r^{n-2}\\sin\\alpha$,
   so ∇⁴ lowers the power by four; a polynomial forcing is matched term by term (undetermined coefficients), and the terms ∇⁴ kills (r,
   r³) carry the boundary conditions. Ch. 8's Stokes operator E² (P199) is a different operator.", code="import sympy as sp\nr, a = sp.symbols('r
   alpha', positive=True)\nlap = lambda F: sp.diff(F, r, 2) + sp.diff(F, r)/r + sp.diff(F, a, 2)/r**2\nprint(sp.simplify(lap(lap(r**5*sp.sin(a)))))
   # 192*r*sin(alpha)")` (**P390**).
10. `nb.derivation("D27", …)` — Part F D27 (14 steps), ref "16.178", with its `check_src` (cached in FAST — P252): holds **N90 [B]** (ψ₀ and w₁,
    ours; **w₁ with the plus sign**) and **N91 [B]**.
11. `nb.worked_example("the flux loss at K = 576", "1. K/576 = 1. 2. Q/(πa²W̄) = 1 − 0.0306 + 0.0120 = **0.9814** (with the exact 1541/50400:
    0.98142). 3. The resistance ratio (16.182): 1 + 0.0306 − 0.0110 = 1.0196. 4. In D: K = D²/16 gives D = 96, and (16.181) gives the same
    numbers. 5. For our arch (a = 12.5 mm, δ = 0.25, ν = 3.3 × 10⁻⁶ m²/s) K = 576 needs a mean speed of only W̄ = κν/(2a√δ) = 16.97 × 3.3 × 10⁻⁶/(2 × 0.0125 × 0.5) ≈ **4.5 mm/s** — the series is a
    slow-flow theory; real arch flows are far beyond it.")`
12. `nb.code` — `ch16.dean_secondary_streamfunction(r, al)`, `ch16.dean_axial_correction(r, al)`, `ch16.dean_axial_velocity(r, al, 288.0)`,
    `ch16.dean_flux_ratio(np.array([100, 288, 576]), exact=True)`, `ch16.curved_resistance_ratio(576.0)`, `ch16.dean_series_reliability(576.0,
    threshold=0.25)`. *expect:* flux 0.99909, 0.99311, 0.98142; λ 1.0196; reliable False at 576 (ours: K ≲ 460 with threshold ¼).
13. `nb.check_agree` — **from scratch:** a 2-D trapezoid of w₀ + κ²w₁ + κ⁴w̄₂ over the disc (with w̄₂ from the closed form of the D27 cell) against
    `dean_flux_ratio(K, exact=True)` to O(K⁴); and ψ₀ satisfies ∇⁴ψ₀ = 8r(1 − r²) sin α by finite differences (1e-6).
14. `nb.figure` — the cross-section: ψ₀ streamlines (purple) over the axial speed w₀ + κ²w₁ at K = 576 (blue heatmap), outer wall on the
    right, the shifted maximum marked. *see / read / change* ("…the shift grows as K: double K, double the shift").
15. `nb.note` — **N92 [B]** "With D: $\dfrac Q{\pi a^2\bar W}=1-0.0306\left(\dfrac D{96}\right)^4+0.0120\left(\dfrac D{96}\right)^8+O(D^{12})$ (16.181) — a change of variable, K = D²/16 (Ĝ in (16.180) is a magnitude — S19)." **N93 [B]** "The resistance ratio
    $\lambda=\dfrac{\lambda_c}{\lambda_s}=1+0.0306\left(\dfrac K{576}\right)^2-0.0110\left(\dfrac K{576}\right)^4+\dots$ (16.182) is the reciprocal of (16.178) (0.0306² − 0.0120 = −0.0111 — checked in a cell; the text swaps 'curved' and 'straight' — S20).
    At large Dean numbers numerical solutions show thin wall layers, a wall shear growing with D, and even four-vortex states — named
    only." + `nb.code`: `1/ch16.dean_flux_ratio(K) - ch16.curved_resistance_ratio(K)` small for K < 300.
16. `nb.plotly` — **F6** `slider_figure` on K (0…600, 25 steps): the cross-section's axial contours (w₀ + κ²w₁) and ψ₀ streamlines (scaled by
    κ). *see / read / change* ("…the vortex pattern keeps its shape; only its strength and the skew grow").
17. `nb.explainer("dean_curved_pipe", heading="What does a bend do to the flow?", why="The secondary flow is three-dimensional: orbit
    the bend to see the tracers corkscrew, drag K to grow the skew, and watch the flux cursor approach the edge of the series.",
    tries=["Set K = 0 then 288: the bullseye slides outward.", "Orbit the 3-D view: helices.", "Switch the chip to D: 96 at K = 576."])`
18. `nb.md` — **What would change if…** "…the vessel's wall could not hold its shape at all — a vein above the heart, an airway in a
    cough? Then the tube's area depends on the flow: C15."

#### C15 — Shapiro's equations (16.204)–(16.206)
1. `nb.core("C15", "Floppy tubes choke at S = 1: $(1-S^2)\dfrac1\alpha\dfrac{d\alpha}{dx}=\dfrac{S^2}{A_0}\dfrac{dA_0}{dx}-\dfrac1{\rho c^2}\left[\dfrac{dp_e}{dx}+\rho g\dfrac{dz}{dx}+\mathrm RQ+\mathrm P\dfrac{dK_p}{dx}\right]$ (16.204)", question="Why can't you suck harder through a soft straw?")`
2. `nb.note` — **N94 [B]** the scene: "A tube whose inside pressure falls below the outside collapses: circle → ellipse → two lobes that
   touch. It happens to veins above the heart, coronary vessels squeezed by the heart muscle, an artery under a blood-pressure cuff,
   airways during a forced cough. Collapsible tubes can also flutter and make sounds." + our three cross-section shapes.
3. `nb.note` — **N95 [C]** "Korotkoff sounds: a cuff pressed between systolic and diastolic pressure makes the artery a collapsible tube that
   opens and closes each beat; the sounds that a stethoscope hears mark the two pressures. The flow-limitation curve below is the reason
   such tubes behave strangely."
4. `nb.note` — **N96 [B]** five-row table of the Starling resistor (a collapsible tube in a pressure chamber, inlet p₁, outlet p₂, chamber
   p_e): Case 1 flow limitation (lowering p₂ stops raising Q) · Case 2 pressure-drop limitation · Case 3 rising resistance · Case 4 flutter
   with a downstream resistor · Case 5 self-excited oscillations. "Case 1 is computed below (ours, model)."
5. `nb.note` — **N97 [B]** "Shapiro's 1-D theory: a general tube law $(p-p_e)=\hat P(A)$ (16.183); small waves travel at $c^2=\dfrac A\rho\dfrac{d(p-p_e)}{dA}$ (16.185) — C08's speed again (Bramwell–Hill).
   The analogy table: speed index S ↔ Mach number M · tube law ↔ equation of state · throat ↔ Laval throat (Ch. 15)." + `ch16.bramwell_hill_speed`.
6. `nb.primer("bending stiffness", "A thin wall resists bending with a stiffness proportional to EI, where the second moment of area of a
   wall of thickness h grows as h³: double the thickness, eight times the bending stiffness. A collapsed tube holds its shape by
   bending, not by hoop tension.", code="E, h1, h2 = 1e6, 1e-3, 2e-3                 # modulus [Pa], two thicknesses [m]\nprint((E*h2**3)/(E*h1**3))
   # 8.0")` (**P392**).
7. `nb.note` — **N98 [B]** "Dimensionless pressure and area: $\mathrm P=\dfrac{p-p_e}{K_p}$, $\alpha=\dfrac A{A_0}$ (16.186), with $K_p\propto EI$ (16.187) — α here is an area ratio (T2)." + `ch16.bending_stiffness_parameter(1e6,
   0.5e-3, 4e-3, 0.5, coeff=1.0)` (the constant is a required keyword: the book gives only the proportionality). **N99 [B]** "The collapsible
   tube law $-\mathrm P\approx\alpha^{-n}-1$ (16.188) with n = 3/2." + `> ⚠️ Common confusion (trap T9): P is negative when the tube is collapsed (α < 1); we plot −P.` + `nb.figure`
   — −P against α for the two fitted forms (our curves, no experimental curve), with the cross-section shapes along the α axis.
8. `nb.derivation("D28", …)` — Part F D28 (4 steps), ref "16.189": holds **N100 [B]** + `ch16.collapsible_wave_speed(np.array([1.0, 0.5, 0.25]),
   vein["Kp"], blood["rho"])`. *expect:* 0.532, 0.895, 1.505 m/s.
9. `nb.recap("R13", "Duct laws with friction", "Steady 1-D continuity $\\dfrac{dA}A+\\dfrac{du}u=0$ (16.190) and momentum with friction and gravity $-A\\,dp-\\tau_ws\\,dx-\\rho gA\\,dz=\\rho Au\\,du$ (16.191),
   the wall stress for turbulent (16.192) or laminar (16.193) flow — Ch. 15's duct laws with a wall force. **New here:** the area is set by the
   tube law, not by the wall's shape. (16.193) uses an area A_e the book does not define (S21; coded only as a keyword).", where="Ch. 15
   §15.3, §15.8 (C03, C11)")` **(re-enter C15)**
10. `nb.note` — **N101 [B]** "The speed index $S=\dfrac uc$, so that $\left(\dfrac{dS^2}{S^2}\right)=2\dfrac{du}u-2\dfrac{dc}c$ (16.195) — the Mach number of a floppy tube; in Ch. 15 every response carried the factor 1 − M²
    (C06), here 1 − S²." + `ch16.speed_index`.
11. `nb.derivation("D29", …)` — Part F D29 (9 steps), ref "16.196–16.202": holds **N102 [B]** + `> ⚠️ **slip S23 — the book prints** $\dfrac{\alpha K_p}{\rho c^2}\dfrac{\partial}{\partial x}\left(\dfrac{\partial\mathrm P}{\partial x}\right)dx$
    in (16.201) **; the correct term is** $\dfrac{\alpha K_p}{\rho c^2}\dfrac{\partial}{\partial x}\left(\dfrac{\partial\mathrm P}{\partial\alpha}\right)dx$.` + `nb.code`: `ch16.collapsible_differentials_sympy()` (all zero) and `printed_16_201=True` (not zero
    for a P that depends on x).
12. `nb.derivation("D30", …)` — Part F D30 (15 steps), ref "16.204–16.206", with its `check_src` + `> ⚠️ Common confusion (trap T9, slip S24): the
    condition at S = 1 needs P dK_p/dx < 0; a collapsed tube has P < 0, so the stiffness must *increase* along the flow — the text's
    'decrease' reads the sign of P wrongly.`
13. `nb.worked_example("a vein at its limit", "1. A collapsed segment at α = 0.25 with K_p = 250 Pa (easy numbers): c² = 1.5 × 250 × 8/1000 =
    3, c = **1.73 m/s**. 2. Its area A = 0.25A₀; with A₀ = 50 mm² the choking flux is Q* = αA₀c = 0.25 × 50 × 10⁻⁶ × 1.73 = **21.7 mL/s**. 3. Our
    vein (A₀ = 50.3 mm², K_p = 200 Pa, L = 0.1 m, open inlet): friction narrows it until S = 1 at the outlet at **Q = 20.0 mL/s**, α = 0.31,
    p₂ − p_e = −946 Pa = −7.1 mmHg (ours, model). 4. Lower p₂ further: Q stays at 20.0 mL/s.")`
14. `nb.code` — `prm = dict(L=vein["L"], A0=np.pi*vein["a0"]**2, Kp=vein["Kp"], n=vein["n"], rho=blood["rho"], mu=blood["mu"])`; `ch16.choking_flux(0.25,
    50e-6, 250.0, 1000.0)`; `ch16.collapsible_march(15e-6, 1.0, (0, prm["L"]), prm)`; `fl = ch16.flow_limitation_curve(np.linspace(-15, 0, 31)*MMHG,
    p1=0.0, pe=0.0, prm=prm)`; `ch16.shapiro_M(0.5, 1.5)`. *expect:* 2.17e-5 m³/s; fl["Q_choke"] 1.999e-5 m³/s, fl["p2_choke"] −946 Pa; M = 0.5.
15. `nb.check_agree` — **from scratch:** an explicit Euler march of dα/dx = −α R′Q/(ρc²(1 − S²)) (laminar R′ = 8πμ/(αA₀)², 20 000 steps) for Q =
    15 mL/s; `assert np.isclose(alpha_euler[-1], ch16.collapsible_march(15e-6, 1.0, (0, 0.1), prm)["alpha"][-1], rtol=1e-3)`; and for a frictionless
    throat A₀(x) the critical flux puts S = 1 exactly at the throat (assert |x_S=1 − x_throat| < 1e-3 L).
16. `nb.figure` — (a) Starling Case 1: Q against p₂ − p_e (−15…0 mmHg): rising, then flat at Q_choke from p₂* = −7.1 mmHg (ours, model); (b) α(x)
    and S(x) along the tube for three outlet pressures, S = 1 dashed. *see / read / change* ("…a more viscous fluid limits at a lower flow").
17. `nb.explainer("collapsible_tube_choking", heading="Why can't you suck harder through a floppy straw?", why="Lowering the outlet
    pressure you watch Q rise and stall as the outlet reaches S = 1; in the throat mode the five drivers of (16.204) change sign exactly
    where 1 − S² does; the giraffe mode balances gravity against friction.", tries=["Drag p₂ from 0 to −15 mmHg: the flow plateaus at
    about −7 mmHg.", "Throat mode: watch the term bars cross zero at S = 1.", "Giraffe mode: the vein stays partly collapsed."])`
18. `nb.note` — the giraffe (part of C15's row): "A giraffe's jugular vein runs downhill about 2 m from head to heart, dz/dx = −1: gravity pulls
    the blood faster than friction alone allows, and the vein runs partly collapsed, its pressure near the outside pressure — the
    gravity driver of (16.204) balancing friction (ours, model in E9)."
19. `nb.md` — **What would change if…** "…the fluid could not flow at all below a certain stress? Then the blood in the middle of a pipe
    stops shearing — C16, the payoff of C02."

#### C16 — Casson pipe flow (16.210), (16.216), (16.220)
1. `nb.core("C16", "A plug riding on a ring: $a_c=a\left(\dfrac{\tau_y}{\tau_w}\right)$ (16.210) and the Casson flux (16.220)", question="Does a yield stress stop blood flowing in a vessel?")`
2. `nb.note` — **N103 [C]** "Below about 200 s⁻¹ blood is clearly non-Newtonian (C02); Casson's law is the usual model."
3. `nb.md` — **Plain words:** "Squeeze toothpaste from a tube: it comes out as a solid cylinder that slides on a thin sheared layer at
   the nozzle wall. A fluid with a yield stress does that in any pipe: near the centre the stress is too small to shear it, so the
   centre moves as a plug. Blood's yield stress is tiny, but the same picture holds — and it costs flow."
4. `nb.md` — **The idea** (ASCII): `stress τ = (Δp/2L) r: 0 on the axis ──► τ_w at the wall` / `   r < a_c: τ < τ_y → plug` / `   r > a_c: sheared ring`.
5. `nb.note` — **N104 [B]** "Casson's law with its yield branch: $\tau^{1/2}=\tau_y^{1/2}+K_c\dot\gamma^{1/2}$ for $\tau\ge\tau_y$, and $\dot\gamma=0$ for $\tau<\tau_y$ (16.207); K_c is dimensional (T14)." + `ch16.casson_shear_rate`.
6. `nb.derivation("D31", …)` — Part F D31 (5 steps), ref "16.210": holds **N105 [B]** (S25; `> ⚠️ **slip S26 — the book prints** its worked
   minimum pressure drop in dyn/cm² 1000× too small **; the minimum is** $\Delta p_{\min}=2L\tau_y/a$, computed below with our numbers in Pa and mmHg.`) and
   **N106 [B]** (trap T21).
7. `nb.primer("algebra with square roots", "Three moves the Casson algebra uses in every step: $(\\sqrt x-\\sqrt y)^2=x+y-2\\sqrt{xy}$; $\\displaystyle\\int r^{1/2}dr=\\tfrac23r^{3/2}$; and a
   polynomial in $t=\\sqrt\\xi$ is still a polynomial (t² = ξ, t⁸ = ξ⁴).", code="x, y = 9.0, 4.0\nprint((np.sqrt(x) - np.sqrt(y))**2, x + y -
   2*np.sqrt(x*y))     # 1.0 1.0\nfrom scipy.integrate import quad; print(quad(np.sqrt, 0, 4)[0], 2/3*4**1.5)   # 5.333 twice")`
   (**P393**).
8. `nb.derivation("D32", …)` — Part F D32 (11 steps), ref "16.216, 16.218": holds **N107 [B]**.
9. `nb.derivation("D33", …)` — Part F D33 (10 steps), ref "16.220", with its `check_src`: holds **N108 [B]**.
10. `nb.worked_example("a plug a quarter of the radius", "a = 0.5 mm, τ_y = 4.8 mPa, K_c² = 2.81 mPa s, Δp/L = 76.8 Pa/m: 1. τ_w = (Δp/L)a/2 =
    19.2 mPa. 2. a_c/a = τ_y/τ_w = 0.25 → a_c = 0.125 mm. 3. Flux bracket: 1 − (16/7)(0.5) + (4/3)(0.25) − (1/21)(0.25)⁴ = **0.190**. 4.
    Poiseuille with μ = K_c²: πa⁴(Δp/L)/(8K_c²) = 6.71 × 10⁻¹⁰ m³/s; Casson: **1.28 × 10⁻¹⁰ m³/s**. 5. The minimum drop for L/a = 1000:
    2Lτ_y/a = 2 × 1000 × 0.0048 = **9.6 Pa = 0.072 mmHg**.")`
11. `nb.code` — `Kc = np.sqrt(wc["Kc"]**2)`; `ch16.casson_plug_radius(76.8*0.5, 0.5, 0.0048)`, `ch16.casson_pipe_velocity(r, 0.5e-3, 38.4, 0.5,
    0.0048, Kc)`, `ch16.casson_core_velocity(...)`, `ch16.casson_pipe_flow_rate(...)`, `ch16.casson_min_pressure_drop(0.5, 0.5e-3, 0.0048)`,
    `ch16.casson_flux_ratio(0.25)`. *expect:* a_c 1.25e-4 m; u_c 2.49e-4 m/s; Q 1.277e-10 m³/s; Δp_min 9.6 Pa; ratio 0.1903.
12. `nb.check_agree` — **from scratch:** `Q_mine = np.pi*ac**2*uc + np.trapezoid(2*np.pi*r*u, r)` over the ring (4001 points); `assert np.isclose(Q_mine,
    ch16.casson_pipe_flow_rate(...), rtol=1e-6)`; τ_y → 0 gives Poiseuille with μ = K_c² (1e-12); a_c ≥ a gives 0.
13. `nb.figure` — the Casson profile across the diameter (blue) with the plug shaded teal, beside Poiseuille at the same Δp with μ = K_c²
    (muted dashed), the stress line (rose) and τ_y; a second panel Q/Q_Pois against a_c/a (0…1). *see / read / change* ("…the first few per cent
    of plug cost the most: the √ξ term").
14. `nb.plotly` — **F7** `slider_figure` on τ_y/τ_w (0…0.95, 20 steps): the Casson profile beside Poiseuille, the plug shaded, Q/Q_Pois in the title.
15. `nb.explainer("casson_blood_flow", heading="Does a yield stress stop the blood?", why="Dragging the pressure gradient down you watch
    the plug widen until it fills the pipe and the flow stops; dragging a shear-rate cursor on the law shows the secant and tangent
    separate.", tries=["Lower Δp/L toward 2τ_y/a: the plug fills the pipe.", "Set τ_y = 0: Poiseuille returns.", "Click inside the plug:
    the shear rate is zero there."])`
16. `nb.md` — **What would change if…** "…the fluid were sap in a tree instead of blood in an artery? The flow is Poiseuille's again,
    but now gravity and dissolved sugars join the pressure — C17."

### A.4 §16.4 Introduction to the Fluid Mechanics of Plants — C17
1. `nb.section("16.4", "Introduction to the Fluid Mechanics of Plants", intro="**What is this section about?** How water rises from soil
   to leaf (xylem) and how sugar travels from leaf to root (phloem): one 'water potential' organises both.")`
2. `nb.core("C17", "Water potential $\psi=p-RTc+\rho gz$ (16.221)", question="How does water reach the top of a 30-metre tree without a pump?")`
3. `nb.note` — **N109 [B]** the scene: "Plants make up most of the Earth's living mass. Their vascular tissue has two pipe systems: xylem
   (water up) and phloem (sugar sap anywhere it is needed)." + our tree sketch (soil, root, stem, leaf, air).
4. `nb.md` — **Plain words:** "Water in a tree climbs tens of metres with no heart. It is not pushed up; it is pulled. Water evaporates
   from the leaves; the water column in the narrow xylem pipes is pulled up after it, hanging from the leaves under tension — a
   negative pressure. To keep track of where water wants to go we need one number that combines pressure, dissolved salt and height."
5. `nb.primer("osmosis and van 't Hoff's law", "Water moves through a membrane toward the side with more dissolved solute. The pressure that
   stops it is the osmotic pressure; for a dilute solution $\\Pi=cR_uT$ (c in mol/m³) — the solute lowers the water's 'potential' by Π.
   (Chemical potential is the deeper idea; named only.)", code="Ru, T = 8.314462618, 293.15           # J/(mol K), K\nfor c in (10, 100,
   300): print(c, c*Ru*T/1e6)       # 0.024, 0.244, 0.731 MPa")` (**P384**).
6. `nb.derivation("D34", …)` — Part F D34 (5 steps), ref "16.221". (🔁 one line before it: "hydrostatics $dp/dz=-\rho g$ (1.8), Ch. 1 D05.")
7. `nb.worked_example("a 30 m tree", "1. Soil water potential ψ_soil = −0.1 MPa (moist soil, ours). 2. Still column: p(30 m) = −0.1 − 1000 ×
   9.81 × 30/10⁶ = **−0.394 MPa** (about −3.9 bar). 3. Moving sap, 1 mm/s in a 25 µm conduit: Poiseuille gradient 8μu/a² = 8 × 10⁻³ ×
   10⁻³/(25 × 10⁻⁶)² = 12.8 kPa/m; over 30 m 0.384 MPa more: **−0.778 MPa** at the top. 4. A leaf cell with 0.3 mol/L of solutes: −R_uTc =
   −0.731 MPa — its ψ is low enough to draw water from the xylem.")`
8. `nb.code` — `ch16.water_potential(-0.394e6, 0.0, tree["T"], 30.0)`, `ch16.water_potential_terms(...)`, `ch16.osmotic_pressure(300.0, tree["T"])`,
   `col = ch16.xylem_column(np.linspace(0, 30, 61), psi_soil=tree["psi_soil"], u=tree["u_xylem"], a=tree["a_xylem"], mu=tree["mu_w"])`,
   `ch16.cavitation_margin(col["p"][-1], p_cav=-2.0e6)`. *expect:* ψ = −0.1 MPa (unchanged with height, still column); 7.312e5 Pa; top p −7.782e5
   Pa; margin 1.22 MPa.
9. `nb.check_agree` — **from scratch:** `psi_mine = p - Ru*T*c + rho*g*z` term by term in SI; `assert np.isclose(psi_mine, ch16.water_potential(p, c,
   T, z))`; the book-unit version (bar, with ρgz divided by 10⁶ — S30) gives the same after `ch16.bar_to_Pa`.
10. `nb.figure` — ψ and p against height in the 30 m tree (still and moving sap), the three terms as stacked bars at soil, root, stem
    top and leaf, the capillary-only height (0.59 m for 25 µm, ch01 `capillary_rise`, ours) marked on the trunk. *see / read / change.*
11. `nb.note` — **N110 [B]** "Xylem: dead conduits (tracheids, vessels); evaporation at the leaves (transpiration) pulls a cohesive column
    under tension; too much tension and the column breaks (cavitation); at night root pressure can push water out at leaf edges
    (guttation). Capillarity alone: Ch. 1's $h=2\sigma\sin\alpha/(\rho gR)$ (Example 1.1) gives 0.59 m for a 25 µm conduit — far too little (ours)." + `ch16.capillary_rise(0.0728,
    np.pi/2, 1000.0, 25e-6)`.
12. `nb.note` — **N111 [B]** "Xylem flow is creeping flow: Re = u·2a/ν = 10⁻³ × 50 × 10⁻⁶/10⁻⁶ = **0.05**; the Poiseuille drop per metre (12.8 kPa)
    is comparable with ρg (9.8 kPa)." + `ch16.xylem_reynolds(1e-3, 25e-6, 1e-6)` + `> ⚠️ **slip S28 — the book quotes** a kinematic viscosity
    ten times water's beside its Reynolds number **; the quoted Re needs water's** ν = 10⁻⁶ m²/s.` (the book's numbers stay private).
13. `nb.note` — **N112 [B]** "Phloem: living sieve tubes joined by perforated sieve plates; sugar loaded at a source raises the solute
    term, water follows by osmosis, the raised pressure pushes the sap to a sink where sugar is unloaded (the Münch pressure-flow idea)
    — ψ again." + a source–sink ASCII sketch.
14. `nb.note` — **N113 [B]** "The pressure drop along a sieve tube with its plates: $\Delta p=\dfrac{8\mu Q}{\pi a^4}\left[L+\dfrac\ell N\left(\dfrac ar\right)^4\right]+2\,\Delta p'$ (16.222) (Δp′ the entry loss at a plate, given in the book)." +
    `nb.code`: `ch16.sieve_tube_pressure_drop(Q, 1.5e-3, 10e-6, 1e-6, 300e-6, 1e-6, 50)` with Q = u πa² (u = 0.3 mm/s) and the Poiseuille limit (N =
    1, r = a, ℓ = 0). *expect:* plates carry 85 % of the drop (illustrative sizes, ours); limit equals 8μQL/(πa⁴).
15. `nb.plotly` — **F8** `slider_figure` on the sap speed (0…3 mm/s, 13 steps): ψ terms and p against height in the 30 m tree.
16. `nb.explainer("xylem_water_potential", heading="How does water climb a 30-metre tree?", why="Dragging height, sap speed and soil
    dryness you watch the three terms of ψ trade along the trunk while the pressure goes below zero and a warning appears near the
    illustrative cavitation limit.", tries=["Set the sap speed to 0: ψ is the same at every height.", "Pick the 100 m giant.", "Drought:
    lower ψ_soil to −1.5 MPa."])`
17. `nb.md` — **What would change if…** "…the tree were twice as tall? The tension doubles — trees stop growing roughly where their xylem
    can no longer hold the column."

### A.5 Exercises, literature, summary
1. `nb.pointer("**S01** Exercises 16.1–16.6 (not reproduced): 16.1 extends C06 (an elastic tube whose radius varies along it); 16.2 is the
   dimensional analysis of C03 (careful: its Strouhal number is aω/U — trap T6); 16.3 a stenosis; 16.4 Shapiro's pure pressure–gravity
   flows (C15; the citation's year is misprinted — S33); 16.5–16.6 generalise C16 to power-law and Herschel–Bulkley fluids (optional
   `ch16.power_law_pipe_flow_rate`, `ch16.herschel_bulkley_pipe_flow_rate`, ours).")`
2. `nb.pointer("**S02** Acknowledgment, literature cited and supplemental reading: not reproduced. Sources we read first-hand are listed in
   the next cell.")` + `nb.md` — references: Secomb & Pries (2013) C. R. Physique 14, 470 (Pries viscosity law); Wikipedia "Womersley number",
   "Pulsatile flow", "Moens–Korteweg equation", "Windkessel effect", "Murray's law", "Hemorheology", "Entrance length (fluid dynamics)";
   Lubashevsky & Gafiychuk arXiv:adap-org/9909003 (Murray's cost); open items: Dean (1928) and Womersley (1955) were not read first-hand.
3. `nb.summary(clicked=[…17 lines…], feeds_forward=[…], left_out=[…])`. **clicked (one per CORE):** C01 "the heart's work per beat is the area of
   its pressure–volume loop, about a joule" · C02 "blood thins with shear: Casson's law is a straight line in square roots" · C03 "α =
   a√(ω/ν) compares a vessel with its Stokes layer" · C04 "resistance ∝ a⁻⁴ — halve the radius, sixteen times the resistance" · C05 "a thin
   plasma layer makes blood seem thinner in narrow tubes" · C06 "an artery is a spring: the tube law" · C07 "the aorta stores each beat and
   drains it with time constant RK" · C08 "the pulse is a wave of wall stretch, c = √(Eh/2ρa₀)" · C09 "Womersley's profile: a breathing
   parabola at small α, a lagging plug at large α" · C10 "flow falls as 1/α² and lags up to 90°" · C11 "viscosity slows and damps the pulse
   by amounts ∝ 1/α" · C12 "Murray: pumping is half the upkeep at the best radius, so Q ∝ a³" · C13 "a fork reflects by admittance mismatch;
   energy balances" · C14 "a bend makes two vortices and pushes the fast core outward" · C15 "a floppy tube chokes at S = 1 — flow limitation"
   · C16 "a yield stress makes a plug and costs flow through √(a_c/a)" · C17 "sap is pulled up under tension: ψ = p − RTc + ρgz". **feeds_forward:**
   "This is the book's last chapter: the tools reused here (Poiseuille, Stokes layers, d'Alembert, influence coefficients, perturbation
   series) are the toolkit for any applied fluid problem — including the climate-system flows of Ch. 13." **left_out:** swimming and flying
   (not in this chapter; Ch. 14's thrust-from-lift explainer is the closest), turbulence in the aorta (Ch. 12), large-Dean-number numerics,
   collapsible-tube oscillations (named in N96).

---

## Part B — explainer storyboards

Common to all ten (and the backup): created with `tools/new_viz.py`; `<meta name="viz:chapter" content="ch16">`; tabs
Walkthrough · Explore · Explain · Derivation · Equations · Code · Check; every displayed number is computed by a JS function that
mirrors a `ch16` callable (Part C) and is proved by `selftest()` parity rows (`py:` expressions as C.12; **never a typed constant
the page could compute** — 37.47°, 0.0306, 2.34 and every preset's verdict are computed). Explain is "Explanation &
interpretation" in numbered sections built with `Viz.work.step / line / box / table / hint / interpret`, modelled on
`forced_damped_vibrations.html` (FDV — the regime-dependent reading) and `amplitude_phase_second_order_II_3.html` (APS — a numbered
derivation with live numbers, one section per optional view): **0** what the views show and what each colour means · **1…n** every
displayed quantity from the controls ("formula = substituted = result — why", results boxed) · a section or hint per view hidden on
phones · the values at the current time or cursor (live) · **Reading the current setting** (regime-dependent, each regime with its
threshold). **Derivation steps are copied from Part F** (the same `did` titles, the same step count and order — ch15 lesson; phones
shorten *why* to its first sentence; plain-text *why* and *watch* never contain raw TeX). Every tour, Explain, Derivation, notes,
status, equation and quiz text that names a book equation **writes it out** next to its number. Colours as header convention 8
(pressure orange, velocity/flow blue, wall/radius/area teal, stress rose, energy/work/cost amber, Newtonian or Poiseuille ghost
muted, viscous/Stokes-layer purple). **Presets carry their distinguishing value in `set`** (ch15 lesson; exactly one chip pressed
after each preset — a getter `get set()` where the value depends on another parameter). **Sizes are decided from the view
(`v.w`, `v.h`), never the window**; the top-right corner of each view stays free for the canvas play button. Every input is
labelled "illustrative" (rule 9); a model says "model" (B1, E9 Starling); every ratio says what it is a ratio of. **Fit:** at 360×640
never three stacked views — two at a time, the third `hidePortrait` and its key number repeated in a visible title or readout; at
844×390 the stage is one row beside the text panel. Walkthrough texts ≤ 45 words, step 1 ≤ 24 words, ≤ 2 extras per step,
`play: false` on steps that quote numbers. JS physics stays light (closed forms, a complex-number helper `Viz.num.C`, Bessel series
of C.12, `Viz.num.brentq`, Simpson); no linear solves beyond 2 × 2.

### E1 · casson_blood_flow
- **Title:** "Does a yield stress stop the blood?" · **Summary:** "Shear blood on a Casson law, then push it through a pipe: watch
  the centre freeze into a plug while the wall layer shears." · **CORE:** C02, C16 (also N19, N21, N22, N104, N105, N106, N108) ·
  **Reference:** FDV (the Explain panel, a verdict by regime) and `overfitting_curves.html` (two linked curves, a marker at the
  special point, a verdict line).
- **meta:** `viz:order 1` · `viz:sections 16.2 16.3` · `viz:equations 16.5 16.6 16.7 16.207 16.209 16.210 16.215 16.216 16.218 16.220`
  · `viz:fluidpy ch16.casson_stress ch16.casson_shear_rate ch16.casson_apparent_viscosity ch16.differential_viscosity
  ch16.whitmore_casson ch16.yield_stress_correlation ch16.casson_plug_radius ch16.casson_pipe_velocity ch16.casson_core_velocity
  ch16.casson_pipe_flow_rate ch16.casson_flux_ratio ch16.casson_min_pressure_drop` · `viz:derivations D02 D31 D32 D33`.
- **Physics (JS ↔ Python):** `tauOf(gd, ty, Kc)` = (√ty + Kc√gd)² ↔ `ch16.casson_stress`; `gdOf(tau, ty, Kc)` (0 below ty) ↔
  `ch16.casson_shear_rate`; `muApp(gd, ty, Kc)` ↔ `ch16.casson_apparent_viscosity`; `muTan(gd, ty, Kc)` ↔ `ch16.differential_viscosity`;
  `ac(dpdx, ty)` = 2ty/dpdx ↔ `ch16.casson_plug_radius(dp, L, tau_y)` with dp = dpdx·L; `u(r, a, dpdx, ty, Kc)` (16.215), plug value inside
  ↔ `ch16.casson_pipe_velocity`; `uc(...)` (16.218) ↔ `ch16.casson_core_velocity`; `Q(...)` (16.220) ↔ `ch16.casson_pipe_flow_rate`;
  `ratio(xi)` ↔ `ch16.casson_flux_ratio`; `dpmin(a, ty)` = 2ty/a per metre ↔ `ch16.casson_min_pressure_drop(1, a, tau_y)`; `tyH(H)` (16.5)
  with C_F = 0.3 ↔ `ch16.yield_stress_correlation(H, 0.3)`. Units: Pa, Pa s, m, Pa/m; K_c from the slider value of K_c².
- **Views** (rows [1, 1]; modes switch the left view): 1. `law` "The material" (row 0, left): mode chips `matview` = "τ(γ̇)" (linear
  axes, γ̇ 0…`gdMax` s⁻¹, τ in mPa), "√τ(√γ̇)" (the straight line of (16.6), intercept √τ_y marked), "μ_app(γ̇)" (log γ̇ 0.1…1000 s⁻¹,
  μ in mPa s, the asymptote K_c² as a muted dashed line labelled "K_c² = 2.81 mPa s"); the Casson curve blue, the Newtonian ghost
  μ = K_c² muted dashed, at the cursor γ̇ a rose secant from the origin and a purple tangent, their slopes in the view title. Pointer:
  drag horizontally to move the γ̇ cursor. 2. `pipe` "The pipe" (row 0, right): velocity across the diameter −a…a (x) against u (mm/s
  or µm/s by size), Casson profile blue (bold), Poiseuille ghost at μ = K_c² muted dashed, plug band teal-shaded with its edges ±a_c
  marked, on a second axis the stress line τ(r) = (Δp/2L)|r| rose with the level τ_y dashed — they cross at a_c. Pointer: click
  sets the inspector radius. 3. `flux` "Flux vs plug size" (row 1, `hidePortrait`): Q/Q_Pois against a_c/a from 0 to 1 (amber),
  the current point as a dot, the first-order slope 1 − (16/7)√ξ as a muted ghost.
- **Controls:** `ty` "Yield stress $\tau_y$" (0…0.03 Pa, step 0.0002, default 0.0048; help "0: Newtonian · 0.005: normal blood
  (Whitmore) · 0.02: very high hematocrit") · `Kc2` "High-shear viscosity $K_c^2$" (1.5…6 mPa s, step 0.05, default 2.809) · `G`
  "Pressure gradient $\Delta p/L$" (0…200 Pa/m, step 0.2, default 76.8; help "below $2\tau_y/a$ nothing flows") · `a` "Radius $a$"
  (0.1…2 mm, step 0.01, default 0.5) · `gd` "Shear rate at the cursor $\dot\gamma$" (0.1…1000 s⁻¹ log, default 4, `optional`) · `H`
  chips "hematocrit (16.5)" 0.30 / 0.45 / 0.60 (`optional`; sets ty from (16.5) with C_F = 0.3 g/100 mL).
- **Presets:** "Newtonian ($\tau_y$ = 0)" `set {ty: 0}` · "normal blood (illustrative)" `set {ty: 0.0048, Kc2: 2.809, G: 76.8, a: 0.5}` ·
  "high hematocrit" `set {ty: 0.016}` · "just above yield" `set {ty: 0.0048, a: 0.5, G: 20.0}` (Δp_min/L = 19.2 Pa/m, a_c/a = 0.96).
- **Status:** a_c < a: "flowing — plug radius {100a_c/a} % of a; flux {100Q/Q_P} % of Poiseuille's" (tone info; tone warn above 60 %) ·
  a_c ≥ a: "no flow — the wall stress {τ_w} Pa is below the yield stress (needs Δp/L > {2τ_y/a} Pa/m)" (tone bad) · τ_y = 0:
  "Newtonian — no plug".
- **Readouts:** "Apparent viscosity μ_app" (mPa s) · "Tangent dτ/dγ̇" (mPa s) · "Wall stress τ_w" (mPa) · "Plug radius a_c/a" · "Flow
  ratio Q/Q_Pois".
- **Depth features:** Explain + Code + Derivation · **linked views** (3, one state) · **modes** (three readings of the law) ·
  **presets** · **status** · **inspector**.
- **Inspector** (click a radius in `pipe`): "r = 0.35 mm: τ = (Δp/2L)r = 38.4 × 0.35 × 10⁻³ = **13.4 mPa** ≥ τ_y = 4.8 mPa → sheared;
  $\dot\gamma=-\dfrac{du}{dr}=\dfrac{\tau+\tau_y-2\sqrt{\tau\tau_y}}{K_c^2}$ (16.211) = (13.4 + 4.8 − 2√(13.4 × 4.8))/2.809 = **1.40 s⁻¹**; u from (16.216) = **…** µm/s" (inside the plug:
  "τ < τ_y → γ̇ = 0, u = u_c").
- **Explain:**
  0. *What the views show* — "**Left:** the blood as a material: how much <b class=c-rose>stress</b> it takes to shear it at a
     rate γ̇. The <b class=c-blue>blue</b> curve is Casson's law, the grey dashed line a Newtonian fluid
     with the same high-shear viscosity. **Right:** the same blood in a pipe: the <b class=c-blue>velocity profile</b>, the
     <b class=c-teal>plug</b> in the middle, and the <b class=c-rose>stress line</b> that grows from 0 at the centre to τ_w at the wall."
  1. *Casson's law at your shear rate* — "$\sqrt\tau=\sqrt{\tau_y}+K_c\sqrt{\dot\gamma}$ (16.207): √τ = √0.0048 + 0.0530 × √4 = 0.0693 + 0.1060 = 0.1753 → $\boxed{\tau=30.7\ \text{mPa}}$.
     Apparent viscosity (secant) τ/γ̇ = 30.7/4 = $\boxed{7.68\ \text{mPa s}}$; tangent dτ/dγ̇ = 4.65 mPa s — smaller, because the curve
     bends (T27)." (live)
  2. *The two limits* — "Fast shear: μ_app → K_c² = 2.81 mPa s (Whitmore's $\sqrt{\tau/\mu_p}=1.53\sqrt{\dot\gamma}+2.0$ (16.7) gives K_c² = 1.53² μ_p = 2.34 μ_p). Slow shear:
     μ_app ≈ τ_y/γ̇ — at γ̇ = 0.1 s⁻¹ already {…} mPa s."
  3. *Stress in the pipe* — "$\tau_{rx}(r)=\dfrac{(p_1-p_2)r}{2L}$ (16.209): τ_w = (Δp/L)a/2 = 76.8 × 0.0005/2 = $\boxed{19.2\ \text{mPa}}$. The line crosses τ_y at
     $a_c=a\,\tau_y/\tau_w$ (16.210) = 0.5 × 4.8/19.2 = $\boxed{0.125\ \text{mm}}$ — a quarter of the radius."
  4. *The plug's speed* — "$u_c=\dfrac1{4K_c^2}\dfrac{\Delta p}L(\sqrt a-\sqrt{a_c})^3\left(\sqrt a+\tfrac13\sqrt{a_c}\right)$ (16.218) = … = $\boxed{0.249\ \text{mm/s}}$, against 1.71 mm/s at the centre of the
     Poiseuille ghost: the plug is slower because the ring that carries all the shear is thinner."
  5. *The flux* — "(16.220): ξ = a_c/a = 0.25 → 1 − (16/7)(0.5) + (4/3)(0.25) − (1/21)(0.25)⁴ = $\boxed{0.190}$; Q = 0.190 × 6.71 × 10⁻¹⁰
     = 1.28 × 10⁻¹⁰ m³/s. The √ξ term does most of the damage."
  6. *Flux vs plug size* — the curve's reading, or the hint "On a phone the flux curve is hidden; the ratio is in the status line."
  7. *Reading the current setting* — τ_y = 0: "A Newtonian fluid: the parabola, no plug." · ξ < 0.1: "A thin plug: the flow is
     close to Poiseuille's at μ = K_c², but already {…} % lower — the √(a_c/a) term." · 0.1 ≤ ξ < 0.6: "The plug matters: a quarter
     of the radius moves solidly and the flux has dropped to a fifth." · 0.6 ≤ ξ < 1: "Nearly blocked: a thin ring shears, the
     plug creeps." · ξ ≥ 1: "No flow: the wall stress cannot reach τ_y. In the body, where pressure drops are tens of pascals
     per centimetre, this happens only where the flow nearly stops."
- **Derivation tab:** **D02** (5 steps) `view: 'law'`; goal `set {matview: 'mu', gd: 4}`; step 3 `live` "μ_app = (√(τ_y/γ̇) + K_c)² =
  (√(0.0048/4) + 0.0530)² = 7.68 mPa s"; step 4 `set {gd: 1000}`, `watch` "the curve meets the dashed K_c² line"; step 5 `set {gd: 0.1}`.
  **D31** (5) `view: 'pipe'`; step 2 `watch` "the rose stress line is straight"; step 4 `live` "a_c = 2Lτ_y/Δp = …"; step 5 `set {G: 19.0}`,
  `watch` "the plug fills the pipe — the status turns red". **D32** (11) `view: 'pipe'`; step 4 `watch` "the slope of the blue profile is
  zero at the plug's edge"; step 11 `live` "u_c = …". **D33** (10) `view: 'flux'` (phones: `pipe`); step 9 `live` "1 − 16/7·t + 4/3·t² −
  t⁸/21 at t = √ξ = …"; step 10 `set {preset: 'normal blood'}`. Interpret: `s => "With your numbers a plug of " + … + " % of the
  radius cuts the flow to " + … + " % of Poiseuille's."`
- **Code:**
  ```python
  Kc = np.sqrt({{Kc2}})                                    # (Pa s)^1/2 from K_c^2 = {{Kc2}} Pa s
  tau = ch16.casson_stress({{gd}}, {{ty}}, Kc)             # (16.207): {{tau}} Pa at your shear rate
  mu_app = ch16.casson_apparent_viscosity({{gd}}, {{ty}}, Kc)   # secant: {{mu}} Pa s
  a_c = ch16.casson_plug_radius({{G}}*L, L, {{ty}})        # (16.210): {{ac}} m
  Q = ch16.casson_pipe_flow_rate({{a}}, {{G}}*L, L, {{ty}}, Kc)  # (16.220): {{Q}} m^3/s
  Q / ch16.poiseuille_flow_rate({{a}}, {{G}}*L, L, Kc**2)       # {{ratio}} of Poiseuille's
  ```
- **Walkthrough (7 steps):** 1. "Is blood a simple liquid?" — "Shear it slowly and fast: does it always resist the same?" `set
  {preset: 'normal blood', matview: 'tau'}` · 2. "A straight line in roots" — "Plot √τ against √γ̇: a straight line that does not
  start at zero. The intercept is √τ_y." `set {matview: 'sqrt'}`, `eq: 'casson'` · 3. "Thinner when sheared" — "Divide stress by
  rate: the apparent viscosity falls from huge to K_c²." `set {matview: 'mu'}`, `derive: {id: 'D02', step: 3}` · 4. "Into a pipe" —
  "Push it with 77 Pa per metre. The stress grows from the centre; where it is below τ_y nothing shears." `eq: 'plug'`, `inspect: true` ·
  5. "A plug, not a stop" — "The centre moves as one solid piece at u_c; the flux is a fifth of Poiseuille's." `code: {id: 'cas',
  lines: [4, 6]}`, `derive: {id: 'D33', step: 9}` · 6. "The real stop" — "Lower the gradient below 2τ_y/a and the plug fills the
  pipe." `set {G: 19.0}` · 7. "Your turn" — "Predict the plug size if τ_y doubles. Then drag." `controls: ['ty', 'G']`.
- **Equations:** `casson` "Casson's law" ref 'Eq. (16.207)' $\tau^{1/2}=\tau_y^{1/2}+K_c\dot\gamma^{1/2}$ (with $\sqrt{\tau/\mu_p}=k_c\sqrt{\dot\gamma}+\sqrt{\tau_y/\mu_p}$ (16.6) in the note, $K_c=k_c\sqrt{\mu_p}$), live at the cursor ·
  `stress` "Stress line" ref 'Eq. (16.209)' $\tau_{rx}(r)=\dfrac{(p_1-p_2)r}{2L}$ · `plug` "Plug radius" ref 'Eq. (16.210)' $a_c=\dfrac{2L\tau_y}{p_1-p_2}=a\dfrac{\tau_y}{\tau_w}$, live · `uc` "Plug speed" ref 'Eq.
  (16.218)' $u_c=\dfrac1{4K_c^2}\dfrac{\Delta p}L(\sqrt a-\sqrt{a_c})^3\left(\sqrt a+\tfrac13\sqrt{a_c}\right)$ · `flux` "Casson flux" ref 'Eq. (16.220)' $Q=\dfrac\pi8\dfrac1{K_c^2}\dfrac{\Delta p}La^4\left[1-\dfrac{16}7\left(\dfrac{a_c}a\right)^{1/2}+\dfrac43\dfrac{a_c}a-\dfrac1{21}\left(\dfrac{a_c}a\right)^4\right]$, live.
- **Check yourself:** (1) "At γ̇ = 4 s⁻¹, which is larger: the secant τ/γ̇ or the tangent dτ/dγ̇?" — "The secant (7.68 against 4.65 mPa s):
  the Casson curve bends down, so the line from the origin is steeper than the tangent." `set {gd: 4}` · (2) "Halve the pressure
  gradient. What happens to the plug radius?" — "It doubles: $a_c=2L\tau_y/(p_1-p_2)$ (16.210) is inversely proportional to the drop." · (3) "Why does a
  plug only 1 % of the radius already lower the flux by a fifth?" — "(16.220) starts 1 − (16/7)√(a_c/a): √0.01 = 0.1, so the drop is
  about 23 %." `set {ty: 0.00048, G: 192, a: 0.5}` (ξ = 0.01) · (4) "What gradient stops a 0.5 mm vessel with
  τ_y = 4.8 mPa?" — "2τ_y/a = 19.2 Pa/m (D31)." 
- **Selftest parity rows:** `{name: 'casson stress', js: tauOf(4, 0.0048, Math.sqrt(0.002809)), py: 'ch16.casson_stress(4.0, 0.0048,
  np.sqrt(0.002809))', rtol: 1e-12}` (0.030723) · `{name: 'mu_app', js: muApp(400, 0.0048, Math.sqrt(0.002809)), py:
  'ch16.casson_apparent_viscosity(400.0, 0.0048, np.sqrt(0.002809))', rtol: 1e-12}` (0.0031883) · `{name: 'tangent', js: muTan(4, …),
  py: 'ch16.differential_viscosity(4.0, 0.0048, np.sqrt(0.002809))', rtol: 1e-12}` · `{name: 'flux', js: Q(5e-4, 76.8, 0.0048,
  Math.sqrt(0.002809)), py: 'ch16.casson_pipe_flow_rate(5e-4, 76.8*0.5, 0.5, 0.0048, np.sqrt(0.002809))', rtol: 1e-10}` (1.2769e-10) ·
  `{name: 'u_c', js: uc(…), py: 'ch16.casson_core_velocity(5e-4, 38.4, 0.5, 0.0048, np.sqrt(0.002809))', rtol: 1e-10}` (2.4919e-4) ·
  invariant `{name: 'xi=1 no flow', js: ratio(1), expect: 0, rtol: 1e-12}`.
- **Fit plan:** 360×640: status (2 lines) · `law` over `pipe` (rows 1 : 1), `flux` hidden (its ratio in the status and a readout);
  `gd` and `H` optional. 844×390: `law` | `pipe` in one row. Desktop and 1000×700: rows [1, 0.8] with `flux` under `pipe`.

### E2 · windkessel_heartbeat
- **Title:** "Why doesn't blood pressure fall to zero between beats?" · **Summary:** "A beating ventricle fills a stretchy aorta that
  drains through a resistor: watch the pressure settle into its rhythm." · **CORE:** C07, C01 (also N06, N07, N09, N12, N14, N45,
  N46, N48) · **Reference:** AFE (system animation and graphs on one clock, presets, an end-of-run card) and FDV (transient →
  periodic state; the Explain panel).
- **meta:** `viz:order 2` · `viz:sections 16.2 16.3` · `viz:equations 16.1 16.3 16.9 16.11 16.48 16.49 16.52` · `viz:fluidpy
  ch16.ejection_half_sine ch16.windkessel_pressure ch16.windkessel_periodic_p0 ch16.windkessel_periodic_state ch16.windkessel_terms
  ch16.windkessel_summary ch16.mean_arterial_pressure ch16.pv_loop_schematic ch16.pv_loop_schematic_work ch16.cardiac_output
  ch16.stroke_volume` · `viz:derivations D01 D09 D10`.
- **Physics (JS ↔ Python):** `Qej(t)` half-sine ↔ `ch16.ejection_half_sine(t, SV, T, Ts)`; `pPer(t)` closed form of C.12 ↔
  `ch16.windkessel_periodic_state`; `p0()` ↔ `ch16.windkessel_periodic_p0`; transient from p(0) = 60 mmHg: beat by beat with the same
  closed form (systole formula with the beat's own starting pressure, then the decay) ↔ `ch16.windkessel_pressure(t, Q_fn, R, K,
  p0)`; `terms(t)` ↔ `ch16.windkessel_terms`; `map3(pS, pD)` ↔ `ch16.mean_arterial_pressure`; the loop: `loop()` with p_min = 5,
  p_ED = 10 mmHg, **p_open = p₀ (the periodic end-diastolic aortic pressure), p_peak = p_S, p_ES = p(T_s)** (the model's own aortic
  pressures — ours) ↔ `ch16.pv_loop_schematic`, work ↔ `ch16.pv_loop_schematic_work`; `CO(SV, HR)` ↔ `ch16.cardiac_output`.
- **Views** (rows [1, 1]): 1. `system` "Heart, aorta, periphery" (row 0, left): a ventricle (squeezing during systole, valve
  drawn open/closed), an aortic balloon whose area grows with the stored volume K(p − p_D,min) (teal), a resistor and an outflow
  arrow (blue) whose width ∝ p/R; labels "inflow Q", "storage K dp/dt", "run-off p/R". 2. `traces` "Pressure and flow" (row 1, full
  width): p(t) orange (mmHg, left axis), Q(t) blue (mL/s, right axis) over the run (faint whole periodic curve, bold "so far"), the
  K → 0 ghost p = QR muted dashed, two horizontal lines — the true mean (amber) and the ⅓-rule estimate (muted dash-dot); the
  current time as a cursor. 3. `loop` "Pressure–volume loop" (row 0, right, `hidePortrait`): our schematic loop (V in mL, p in mmHg)
  with the corners A–D, the area shaded amber, a dot tracing it in sync with the clock (filling in diastole, ejection in systole).
- **Controls:** `K` "Compliance $\mathrm K$" (0.3…4 mL/mmHg, step 0.05, default 1.60 — code uses m³/Pa = mL/mmHg ÷ 1.3332 × 10⁸; help
  "small: a stiff aorta") · `R` "Peripheral resistance $R$" (0.5…3 × 10⁸ Pa s/m³, step 0.05, default 1.45) · `HR` "Heart rate"
  (40…180 per min, step 1, default 66) · `SV` "Stroke volume" (40…120 mL, step 1, default 75) · `share` "Systole share of the beat"
  (0.25…0.5, step 0.01, default 0.35, `optional`).
- **Transport:** `t` 0…8 beats (unit s, max 8T), rate 0.5 beat/s, `end: 'hold'`; end-of-run card (`Viz.card`): "Periodic after
  {n} beats (p_D within 1 mmHg of its periodic value). p_S {…}, p_D {…}, mean {…} mmHg."
- **Presets:** "resting adult (illustrative)" `set {K: 1.60, R: 1.45, HR: 66, SV: 75, share: 0.35}` · "stiff aorta" `set {K: 0.80}`
  (expect p_S 123, p_D 61 mmHg) · "high resistance" `set {R: 2.20}` (152/121) · "exercise" `set {HR: 120, SV: 95, R: 0.90, share: 0.40}`
  (147/110) · "rigid aorta (K → 0)" `set {K: 0.30}`.
- **Status:** pulse pressure p_S − p_D: < 30 mmHg "narrow pulse ({…} mmHg) — a compliant aorta smooths the beat" · 30–60 "normal-looking
  pulse ({…} mmHg, illustrative)" · > 60 "wide pulse ({…} mmHg) — a stiff aorta passes the heart's surge on"; plus " · RK = {…} s =
  {…} beats".
- **Readouts:** "Cardiac output CO" (L/min) · "Time constant RK" (s) · "p_S / p_D" (mmHg) · "Mean (true / ⅓ rule)" (mmHg) · "Work per
  beat W" (J).
- **Terms** (title "Where the inflow goes now", unit mL/s): `storage` "storage K dp/dt" (teal; negative in diastole), `runoff` "run-off
  p/R" (blue), total "inflow Q" (equals their sum at every instant).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** with an end-of-run card · **presets** ·
  **terms** · **status**.
- **Explain:**
  0. *What the views show* — "**Top left:** the system: the <b class=c-blue>inflow</b> from the ventricle, the aorta's
     <b class=c-teal>stored volume</b>, the <b class=c-blue>run-off</b> through the resistor. **Bottom:** <b class=c-orange>aortic
     pressure</b> and <b class=c-blue>inflow</b> in time; the grey dashed curve is what a rigid aorta (K → 0) would give, p = QR.
     **Top right:** the ventricle's pressure–volume loop; its <b class=c-amber>area</b> is the work per beat."
  1. *What the heart delivers* — "$SV=EDV-ESV$ (16.1) = 130 − 55 = 75 mL; $CO=SV\times HR$ (16.3) = 75 mL × 66/min = $\boxed{4.95\ \text{L/min}}$ = 8.25 × 10⁻⁵ m³/s;
     half-sine ejection over T_s = 0.35 × 0.909 s = 0.318 s, peak Q_max = π SV/(2T_s) = $\boxed{370\ \text{mL/s}}$."
  2. *The aorta's memory* — "$Q(t)=\mathrm K\dfrac{dp}{dt}+\dfrac{p(t)}{R}$ (16.49): τ = RK = 1.45 × 10⁸ × 1.2 × 10⁻⁸ = $\boxed{1.74\ \text{s}}$ — about two beats. In diastole
     $p(t)=p(T_s)e^{-(t-T_s)/R\mathrm K}$: over the 0.591 s of diastole the pressure falls by the factor e^{−0.591/1.74} = $\boxed{0.712}$."
  3. *The periodic state* — "(16.52) with p(T) = p₀ gives p₀ = $\boxed{74.9\ \text{mmHg}}$ (ours); the beat then peaks at p_S = $\boxed{106.0}$ and
     dips to p_D = $\boxed{74.5}$ mmHg — a pulse pressure of 31.5 mmHg (illustrative numbers, not a patient)."
  4. *The mean is set by R alone* — "$\bar p=R\bar Q$ = 1.45 × 10⁸ × 8.25 × 10⁻⁵ Pa = $\boxed{89.7\ \text{mmHg}}$. The ⅓ rule $\bar p_A=\tfrac13p_S+\tfrac23p_D$ (16.11) gives 85.0 mmHg —
     5 % low for this waveform: an estimate, not a definition (T15)."
  5. *Where the inflow goes, now* — "t = {Viz.live('t')} s: Q = {…} mL/s = storage {…} + run-off {…} (the bars). In systole most of
     the stroke is stored; in diastole the store pays out." (live)
  6. *The work per beat* — "(D01) W = −∮p dV = area of the loop = $\boxed{0.94\ \text{J}}$ with this model's aortic pressures; cardiac power W ×
     HR = 1.03 W." Hint on phones: "The loop is hidden on a phone; its work is the readout W."
  7. *Reading the current setting* — RK < 0.5T: "Short memory: the pressure nearly empties every beat — a wide, peaky pulse, like a
     rigid tube." · 0.5T ≤ RK ≤ 3T: "The aorta carries the circulation through diastole: the pressure falls only partway before
     the next beat." · RK > 3T: "Long memory: the pressure barely moves; a very compliant aorta or a high resistance." Plus the
     stiffness sentence: "Halve K and the pulse pressure grows while the mean stays at R Q̄ — why stiff old aortas raise systolic
     pressure."
- **Derivation tab:** **D01** (6 steps) `view: 'loop'` (phones: `system`); step 3 `watch` "the two vertical legs of the loop";
  step 5 `live` "W = {…} J"; step 6 `set {preset: 'resting adult'}`. **D09** (4) `view: 'system'`; step 2 `watch` "the balloon grows
  when the pressure rises"; step 4 `terms: true`. **D10** (9) `view: 'traces'`; step 6 `set {t: 0.6}` (diastole), `watch` "the orange
  curve decays exponentially"; step 7 `set {K: 4.0}`; step 8 `live` "p₀ = {…} mmHg"; step 9 `live` "mean {…} = R Q̄ {…} mmHg".
  Interpret: `s => "With RK = " + … + " s the pressure loses only " + … + " % between beats."`
- **Code:**
  ```python
  T, Ts = 60/{{HR}}, {{share}}*60/{{HR}}                      # period and systole [s]
  Q = lambda t: ch16.ejection_half_sine(t, {{SV}}e-6, T, Ts) # inflow, peak {{Qmax}} mL/s
  p = ch16.windkessel_periodic_state(t, SV={{SV}}e-6, HR={{HR}}, systole_share={{share}}, R={{R}}e8, K={{K}}/1.3332e8)
  s = ch16.windkessel_summary(SV={{SV}}e-6, HR={{HR}}, systole_share={{share}}, R={{R}}e8, K={{K}}/1.3332e8)
  s["pS"]/MMHG, s["pD"]/MMHG, s["p_mean"]/MMHG               # {{pS}}, {{pD}}, {{pm}} mmHg
  ch16.mean_arterial_pressure(s["pS"], s["pD"])/MMHG          # 1/3 rule: {{p3}} mmHg
  ```
- **Walkthrough (7 steps):** 1. "A pump that stops" — "The heart pushes for a third of each beat. Why doesn't the pressure drop to
  zero in between?" `set {preset: 'resting adult', t: 0}` · 2. "Store and drain" — "Press play. The aorta swells as blood arrives and
  drains through the small vessels." `play: true`, `terms: true` · 3. "One equation" — "Inflow = storage + run-off: $Q(t)=\mathrm K\dfrac{dp}{dt}+\dfrac{p(t)}{R}$ (16.49)." `eq: 'wk'`,
  `derive: {id: 'D09', step: 4}` · 4. "The memory RK" — "In diastole the pressure decays as e^{−t/RK}; RK = 1.74 s, two beats long."
  `set {t: 0.6}`, `derive: {id: 'D10', step: 6}` · 5. "Stiffen the aorta" — "Halve K. The pulse widens, the mean stays put." `set
  {K: 0.80}`, `play: true`, `code: {id: 'wk', lines: [5, 6]}` · 6. "The heart's bill" — "The loop's area is the work per beat: about one
  joule." `derive: {id: 'D01', step: 5}` · 7. "Your turn" — "Predict p_S at 120 beats per minute. Then drag." `controls: ['HR', 'K', 'R']`.
- **Equations:** `sv` "Stroke volume and output" ref 'Eq. (16.1)' $SV=EDV-ESV$ with $CO=SV\times HR$ (16.3) · `wk` "The Windkessel" ref 'Eq. (16.49)'
  $Q(t)=\mathrm K\dfrac{dp}{dt}+\dfrac{p(t)}{R}$, live terms · `sol` "Its solution" ref 'Eq. (16.52)' $p(t)=\dfrac1{\mathrm K}e^{-t/R\mathrm K}\displaystyle\int_0^tQ(\tau)e^{\tau/R\mathrm K}d\tau+p_0e^{-t/R\mathrm K}$ · `res` "Peripheral resistance" ref
  'Eq. (16.9)' $R=\Delta\bar p/Q$, live · `map` "The ⅓ rule" ref 'Eq. (16.11)' $\bar p_A=\tfrac13p_S+\tfrac23p_D$, live with the true mean beside it.
- **Check yourself:** (1) "Halve the compliance. Which changes more, the mean pressure or the pulse pressure?" — "The pulse
  pressure; the mean is R Q̄ whatever K is (D10 step 9)." `set {K: 0.80}` · (2) "What does the ⅓ rule get wrong here?" — "It assumes a
  waveform; for this half-sine beat it is 5 % low. (16.11) is an estimate." · (3) "Raise the heart rate to 120 with everything else
  fixed. Why does the diastolic pressure rise?" — "Diastole is shorter, so the pressure has less time to decay (e^{−t/RK}), and CO is
  larger, so the mean R Q̄ rises." `set {HR: 120}` · (4) "Where is the stroke volume on the loop?" — "Its width, EDV − ESV (16.1)."
- **Selftest parity rows:** `{name: 'p0', js: p0(), py: 'ch16.windkessel_periodic_p0(SV=75e-6, HR=66.0, systole_share=0.35, R=1.45e8,
  K=1.2e-8)', rtol: 1e-10}` (9988.35 Pa) · `{name: 'p(t)', js: pPer(0.2), py: 'ch16.windkessel_periodic_state(0.2, SV=75e-6, HR=66.0,
  systole_share=0.35, R=1.45e8, K=1.2e-8)', rtol: 1e-10}` · `{name: 'storage', js: terms(0.159091).storage, py:
  'ch16.windkessel_terms(0.159091, SV=75e-6, HR=66.0, systole_share=0.35, R=1.45e8, K=1.2e-8)["storage"]', rtol: 1e-8}` (2.865 × 10⁻⁴
  m³/s) · `{name: 'loop work', js: loopWork(…), py: 'ch16.pv_loop_schematic_work(130e-6, 55e-6, p_min=5*133.322387415,
  p_ED=10*133.322387415, p_open=80*133.322387415, p_peak=125*133.322387415, p_ES=100*133.322387415)', rtol: 1e-12}` (1.06658 J) ·
  `{name: 'CO', js: CO(75e-6, 66), py: 'ch16.cardiac_output(75e-6, 66.0)', rtol: 1e-12}` · invariant `{name: 'mean = RQ', js:
  meanP()/(1.45e8*75e-6*66/60), expect: 1, rtol: 1e-6}`.
- **Fit plan:** 360×640: status (1 line) · `system` over `traces` (rows 0.9 : 1.1), `loop` hidden (its W a readout); `share`
  optional; terms bars appear on steps 2 and 5. 844×390: `system` | `traces`. Desktop: rows [1, 1] with `loop` beside `system`.

### E3 · pulse_wave_artery
- **Title:** "Is the pulse at your wrist the blood arriving?" · **Summary:** "Release a pressure pulse in an elastic artery: the
  wall's stretch travels at metres per second while the blood barely moves." · **CORE:** C08, C06 (also N37, N51, N52, N53, R09,
  R10) · **Reference:** ch15 E1 `sound_pulse_dalembert` (two shapes sliding apart, a tracer strip with an honest exaggeration
  label) and FDV (Explain panel).
- **meta:** `viz:order 3` · `viz:sections 16.3` · `viz:equations 16.39 16.41 16.55 16.56 16.57 16.58 16.61 16.62 16.63 16.68` ·
  `viz:fluidpy ch16.moens_korteweg_speed ch16.youngs_modulus_from_pwv ch16.tube_radius_hookean ch16.hoop_stress ch16.wall_bulge
  ch16.pulse_state ch16.velocity_from_pulse ch16.linearisation_ratios ch16.tube_law_linear` · `viz:derivations D08 D11 D12`.
- **Physics (JS ↔ Python):** `c(E, h, a0, rho, midwall)` ↔ `ch16.moens_korteweg_speed`; `aOf(p, a0, E, h)` ↔
  `ch16.tube_radius_hookean`; `bulge(p)` ↔ `ch16.wall_bulge`; `state(x, t, ic)` → {p, u, a}: Gaussian time shape of width w (in
  metres via c), "right": p = P exp(−((x − ct)/w)²), "rest": two halves at ±ct ↔ `ch16.pulse_state(x, t, *, c, amplitude, width, rho,
  a0, E, h, ic)`; `uOf(p)` = ±p/(ρc) ↔ `ch16.velocity_from_pulse`; `ratios()` ↔ `ch16.linearisation_ratios(p_amp, p_mean, u, c,
  Eh_over_a0=…)`; tracer displacement ∫u dt (closed form for the Gaussian: the error function, `Viz.num.erf`); `Einv(c)` ↔
  `ch16.youngs_modulus_from_pwv`.
- **Views** (rows [1, 0.9]): 1. `artery` "The artery from the side" (row 0, full width): x from 0 to 1.0 m; upper and lower walls
  teal, displaced by the bulge η × N (exaggeration N chosen so the largest bulge is 25 % of the drawn radius; label "wall bulge ×{N};
  true {η_max} mm"), the lumen tinted by pressure (orange), 40 blue tracer particles moving by their true displacement × the same
  factor (label "tracers ×{N}"), a probe line. 2. `graph` "p(x) and u(x)" (row 1, left): p (mmHg, orange) and u (cm/s, blue) along x
  at the current time, the right- and left-going halves faint when ic = "rest", the arrival marker "front at x = ct". 3. `speed`
  "Wave speed vs stiffness" (row 1, right, `hidePortrait`): c (m/s) against E (0.1…3 MPa, log) for the current h/a₀ (blue), the
  thick-wall (16.63) muted dashed, a band "illustrative young aorta … stiff older aorta" (E 0.3–1.6 MPa), the current point.
- **Controls:** `E` "Young's modulus $E$" (0.1…3 MPa, log, default 0.4) · `ha` "Wall thickness $h/a_0$" (0.05…0.3, step 0.005, default
  0.125) · `a0` "Radius $a_0$" (2…15 mm, step 0.5, default 12) · `P` "Pulse amplitude" (5…60 mmHg, step 1, default 25) · `ic` chips
  "right-going pulse" / "bump released at rest" (default right).
- **Transport:** `t` 0…0.25 s, rate 0.05 s/s (slow motion; label "slow motion ×20"), `end: 'hold'`; end-of-run card: "The pulse
  covered {c·t} m; each tracer moved only {Δx} mm."
- **Presets:** "young aorta (illustrative)" `set {E: 0.4, ha: 0.125, a0: 12}` · "stiff older aorta" `set {E: 1.6, ha: 0.125, a0: 12}`
  (c 9.7 m/s) · "soft vein-like wall" `set {E: 0.1, ha: 0.05, a0: 5}` · "nearly rigid tube" `set {E: 3.0, ha: 0.3, a0: 2}`.
- **Status:** linearity check p̂/(Eh/a₀): < 0.1 "linear theory fine (p̂/(Eh/a₀) = {…}; u/c = {…})" (info) · 0.1–0.3 "mildly nonlinear —
  (16.58) is a first-order law" (warn) · > 0.3 "too large for $p-p_e=\dfrac{Eh}{2a_0A_0}(A-A_0)$ (16.58) — the pulse would steepen" (bad).
- **Readouts:** "Pulse speed c" (m/s) · "Blood speed under the crest û" (m/s) · "Wall bulge" (mm) · "Arrival over 0.5 m" (ms) · "p̂/(Eh/a₀)".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** with an end-of-run card · **presets** ·
  **inspector** · **status**.
- **Inspector** (click a point of `artery` or `graph`): "x = 0.40 m, t = 0.08 s: p = 25 exp(−((0.40 − 4.86 × 0.08)/0.1)²) = **24.9
  mmHg**; $\hat p_1=\rho c\hat u_1$ (16.68): u = p/(ρc) = 3318/(1060 × 4.86) = **0.645 m/s**; tube law (16.41): a = 12/(1 − 12 × 3318/(0.4 × 10⁶ ×
  1.5)) = **12.85 mm**."
- **Explain:**
  0. *What the views show* — "**Top:** the artery. The <b class=c-teal>wall</b> bulges where the <b class=c-orange>pressure</b> is
     high; the <b class=c-blue>dots</b> are blood parcels. Both motions are exaggerated by the same factor (label). **Bottom left:**
     pressure and blood speed along the artery now. **Bottom right:** how the speed depends on the wall's stiffness."
  1. *The wall as a spring* — "$\sigma_{\theta\theta}=\dfrac{(p-p_e)a}{h}$ (16.39) and $a=a_0\left[1-\dfrac{a_0}{Eh}(p-p_e)\right]^{-1}$ (16.41): Eh/a₀ = 0.4 × 10⁶ × 0.125 = $\boxed{50\ \text{kPa}}$ (375 mmHg). A 25 mmHg pulse widens the
     radius by $\boxed{0.86\ \text{mm}}$ (7 %)."
  2. *The pulse speed* — "$c=\sqrt{\dfrac{Eh}{2\rho a_0}}$ (16.62) = √(50 000/(2 × 1060)) = $\boxed{4.86\ \text{m/s}}$; the thick-wall form $c=\sqrt{\dfrac{Eh}{2\rho(a_0+h/2)}}$ (16.63) gives
     4.71 m/s. Over 0.5 m the pulse needs $\boxed{103\ \text{ms}}$."
  3. *How fast the blood moves* — "$\hat p_1=\rho c\hat u_1$ (16.68): û = 3333/(1060 × 4.86) = $\boxed{0.65\ \text{m/s}}$ under the crest — but only while the crest
     passes; afterwards the parcel stops, having moved ∫u dt ≈ {…} mm."
  4. *Is it linear?* — "p̂/(Eh/a₀) = 0.067 and u/c = 0.13: both small, so (16.56)–(16.58) hold (D11)."
  5. *Stiffness meter* — the `speed` view's reading, or the hint "On a phone the c–E curve is hidden; c is the first readout."
     "Inverse: $E=2\rho a_0c^2/h$ — measure c, get E."
  6. *Now* — "t = {Viz.live('t')} s: the front is at x = {Viz.live('xf')} m."
  7. *Reading the current setting* — E < 0.2 MPa: "A soft wall: slow pulse, big bulge — a vein-like tube." · 0.2–0.8: "A young
     aorta: the pulse arrives in about a tenth of a second." · > 0.8: "A stiff wall: the pulse races ahead — clinicians read high
     pulse-wave velocity as arterial stiffening."
- **Derivation tab:** **D08** (6) `view: 'artery'`; step 3 `live` "σ = (p − p_e)a/h = …"; step 6 `set {P: 60}`, `watch` "the bulge
  grows faster than the pressure". **D11** (7) `view: 'graph'`; step 2 `watch` "nothing on screen uses the printed (16.54)"; step 6 `live`
  "∂p/∂A = Eh/(2a₀A₀) = …". **D12** (7) `view: 'graph'` (the `speed` view on desktop); step 4 `set {ic: 'rest'}`, `play`; step 7 `live` "c = √(Eh/(2ρa₀))
  = …". Interpret: `s => "Your wall gives c = " + … + " m/s: the blood under the crest moves " + … + " times slower."`
- **Code:**
  ```python
  c = ch16.moens_korteweg_speed({{E}}e6, {{h}}, 1060.0, {{a0}})   # (16.62): {{c}} m/s
  eta = ch16.wall_bulge({{P}}*MMHG, a0={{a0}}, E={{E}}e6, h={{h}})  # (16.41): {{eta}} m
  u = ch16.velocity_from_pulse({{P}}*MMHG, 1060.0, c)              # (16.68): {{u}} m/s
  s = ch16.pulse_state(x, {{t}}, c=c, amplitude={{P}}*MMHG, width=0.1,
                       rho=1060.0, a0={{a0}}, E={{E}}e6, h={{h}}, ic="{{ic}}")
  ch16.linearisation_ratios({{P}}*MMHG, 100*MMHG, u, c, Eh_over_a0={{Eha}})
  ```
- **Walkthrough (7 steps):** 1. "What do you feel at the wrist?" — "Blood arriving, or something else?" `set {preset: 'young aorta',
  t: 0}` · 2. "Release a pulse" — "Press play. The bulge runs along the artery at about 5 m/s." `play: true` · 3. "The blood barely
  moves" — "Watch the dots: each jumps forward under the crest and stops — millimetres, not metres." `inspect: true`, `eq: 'imp'` ·
  4. "Why it travels" — "The swelling pushes blood forward, which swells the next slice: a wave equation, $\dfrac{\partial^2p}{\partial t^2}=c^2\dfrac{\partial^2p}{\partial x^2}$ (16.61)." `eq: 'wave'`,
  `derive: {id: 'D12', step: 4}` · 5. "Stiffness sets the speed" — "$c=\sqrt{Eh/(2\rho a_0)}$ (16.62): quadruple E and the speed doubles." `set {E: 1.6}`, `play: true`,
  `code: {id: 'mk', lines: [1, 1]}` · 6. "Two halves" — "A bump released at rest splits into two halves running apart." `set {ic: 'rest',
  t: 0}`, `play: true` · 7. "Your turn" — "Predict the arrival time over 0.5 m for E = 0.8 MPa. Then drag." `controls: ['E', 'ha', 'a0']`.
- **Equations:** `law` "Tube law" ref 'Eq. (16.41)' $a=a_0\left[1-\dfrac{a_0}{Eh}(p-p_e)\right]^{-1}$, live · `lin` "Linearised laws" ref 'Eq. (16.56)–(16.58)' $\dfrac{\partial A}{\partial t}+A_0\dfrac{\partial u}{\partial x}=0$, $\rho\dfrac{\partial u}{\partial t}=-\dfrac{\partial p}{\partial x}$,
  $p-p_e=\dfrac{Eh}{2a_0A_0}(A-A_0)$ · `wave` "Wave equation" ref 'Eq. (16.61)' $\dfrac{\partial^2p}{\partial t^2}=c^2\dfrac{\partial^2p}{\partial x^2}$ · `mk` "Moens–Korteweg" ref 'Eq. (16.62)' $c=\sqrt{\dfrac{Eh}{2\rho a_0}}=\sqrt{\dfrac A\rho\dfrac{dp}{dA}}$,
  live · `imp` "Pressure and blood speed" ref 'Eq. (16.68)' $\hat p_1=\rho c\hat u_1$, $\hat p_2=-\rho c\hat u_2$, live.
- **Check yourself:** (1) "Quadruple E. What happens to c?" — "It doubles: c ∝ √E (16.62)." `set {E: 1.6}` · (2) "Under the crest, how
  fast does the blood move compared with the pulse?" — "û/c = p̂/(ρc²) ≈ 0.13: the blood is about eight times slower and moves only
  while the crest passes." · (3) "Release a bump at rest. How high is each half, and which way does the blood under the left
  half move?" — "Half the height; backward, $\hat p_2=-\rho c\hat u_2$ (16.68)." `set {ic: 'rest'}` · (4) "Why does a bigger radius at the same h give a
  slower pulse?" — "The hoop stress (p − p_e)a/h is larger, so the wall stretches more per pascal: a softer spring."
- **Selftest parity rows:** `{name: 'MK', js: c(0.4e6, 1.5e-3, 12e-3, 1060, false), py: 'ch16.moens_korteweg_speed(0.4e6, 1.5e-3,
  1060.0, 12e-3)', rtol: 1e-12}` (4.856429) · `{name: 'midwall', js: c(…, true), py: 'ch16.moens_korteweg_speed(0.4e6, 1.5e-3, 1060.0,
  12e-3, midwall=True)', rtol: 1e-12}` (4.711428) · `{name: 'bulge', js: bulge(25*133.322387415), py: 'ch16.wall_bulge(25*133.322387415,
  a0=12e-3, E=0.4e6, h=1.5e-3)', rtol: 1e-12}` (8.5707e-4 m) · `{name: 'u', js: uOf(3333.06), py: 'ch16.velocity_from_pulse(3333.06,
  1060.0, 4.856429311786321)', rtol: 1e-12}` · `{name: 'state p', js: state(0.4, 0.08, 'right').p, py: 'ch16.pulse_state(0.4, 0.08,
  c=4.856429311786321, amplitude=3333.06, width=0.1, rho=1060.0, a0=12e-3, E=0.4e6, h=1.5e-3, ic="right")["p"]', rtol: 1e-10}` ·
  invariant `{name: 'half height', js: state(c(0.4e6, 1.5e-3, 12e-3, 1060, false)*0.2, 0.2, 'rest').p/3333.06, expect: 0.5, rtol: 1e-6}` (each half of a bump
  released at rest carries half its height).
- **Fit plan:** 360×640: status (1 line) · `artery` (55 %) over `graph` (45 %), `speed` hidden (c in the `artery` title "c = 4.86 m/s");
  the pulse amplitude slider stays, `a0` optional. 844×390: `artery` over a row of `graph` | `speed` beside the panel. Desktop: rows
  [1, 0.9].

### E4 · womersley_profiles
- **Title:** "Is pulsing flow just Poiseuille flow sloshing?" · **Summary:** "Drive blood back and forth in a rigid tube and drag
  the Womersley number: a breathing parabola turns into a lagging plug with a Stokes layer." · **CORE:** C09, C03 (also N24, N55,
  N57, N58, N59) · **Reference:** AFE (one clock across linked views, regime notes with a highlighted table) and ch08 E7 (the
  oscillating plate → Stokes layer).
- **meta:** `viz:order 4` · `viz:sections 16.3` · `viz:equations 16.12 16.73 16.76 16.80 16.82 16.87 16.88` · `viz:fluidpy
  ch16.womersley_number ch16.womersley_F1 ch16.womersley_velocity_nd ch16.bessel_j0_series ch16.womersley_small_alpha
  ch16.womersley_large_alpha ch16.womersley_regime ch16.vessel_scales_table ch16.stokes_layer` · `viz:derivations D03 D13 D14 D15`.
- **Physics (JS ↔ Python):** complex arithmetic with `Viz.num.C`; `J0(z)`, `J1(z)` by the 100-term series (C.12) ↔
  `ch16.bessel_j0_series`; `F1(s, al)` ↔ `ch16.womersley_F1(s, 1, al)`; `und(s, wt, al)` = Re[(4/(iα²))(1 − F₁)e^{iωt}] ↔
  `ch16.womersley_velocity_nd`; limit ghosts `small(s, wt)` = (1 − s²)cos ωt ↔ `ch16.womersley_small_alpha` (scaled), `large(s, wt,
  al)` (16.88) scaled ↔ `ch16.womersley_large_alpha`; `regime(al)` ↔ `ch16.womersley_regime`; vessel α from a, HR, ν ↔
  `ch16.womersley_number`; δ_e ↔ `ch16.stokes_layer`.
- **Views** (rows [1.2, 0.8]): 1. `profile` "Velocity across the tube" (row 0, left): x = r/a from −1 to 1 (full diameter — trap T20),
  u/(Ga²/4μ) on a fixed axis −1.05…1.05 (rescaled by 4/α² beyond α = 4 with the factor in the title), the current profile bold blue,
  the last 11 phases (every 15°) faint, the two limit ghosts (quasi-steady parabola, Stokes-layer plug) muted dashed when "limits"
  is on, a purple band of width δ_e/a at each wall. 2. `time` "Gradient and speeds over a cycle" (row 1, full width): G cos ωt
  (orange, scaled to the peak), centreline u (blue), near-wall u at r = 0.9a (purple), over 0…2π with the phase cursor; the lag of
  the centre marked as an arc "φ_c = {…}°". 3. `vessels` "Where vessels sit" (row 0, right, `hidePortrait`): a horizontal log axis
  of α from 0.001 to 30 with the vessels of `ch16.vessel_scales_table(HR=66.0)` (cited typical sizes; when probed: capillary 0.006,
  arteriole 0.029, venule 0.017, large artery 2.2, vein 2.9, aorta 15.9, vena cava 20.3) and the current α as a marker; regime bands (quasi-steady < 1, transitional 1–10, inertial
  > 10).
- **Controls:** `al` "Womersley number α" (0.2…25, log, default 3) · `wt` (transport, the phase ωt) · `lim` toggle "show the two
  limits" (default on) · `vessel` chips (one per row of `vessel_scales_table`; each sets α) · `mode` "set α from a
  vessel" (`optional`; reveals a [0.004…15 mm], HR [40…180], ν [2…5 × 10⁻⁶ m²/s] sliders).
- **Transport:** `wt` 0…2π, rate 1 rad/s, `end: 'loop'`.
- **Presets:** "quasi-steady (α = 0.5)" `set {al: 0.5}` · "artery (α = 2.9, illustrative)" `set {al: 2.894}` · "α = 8" `set {al: 8}` ·
  "aorta-like (α = 17, illustrative)" `set {al: 17.36}`.
- **Status:** `regime(α)`: "quasi-steady (α < 1): Poiseuille at every instant, in step with the gradient" · "transitional (1 ≤ α ≤ 10):
  the centre lags by {φ_c}°" · "inertial (α > 10): a flat plug a quarter cycle late, viscosity only in a layer {δ_e/a} of the radius
  thick".
- **Readouts:** "α" · "Stokes layer δ_e/a = √2/α" · "Centre amplitude" (of Ga²/4μ) · "Centre lag φ_c" (°) · "Reversal at the wall first?" (yes/no
  over the cycle).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** · **status** · **inspector** ·
  **notes**.
- **Inspector** (click a radius in `profile`): "r/a = 0.90, α = 3: z = i^{3/2}·3·0.9 = −1.909 + 1.909i; J₀(z) = …; F₁ = J₀(z)/J₀(i^{3/2}·3) =
  … ; u/(Ga²/4μ) = Re[(4/(iα²))(1 − F₁)e^{iωt}] = **…**" (the exact complex arithmetic, 3 significant figures).
- **Notes** (regime-dependent, under the stage on desktop): small α "Viscosity reaches across the whole tube within one cycle." ·
  middle "The core starts to lag; near the wall the flow reverses before the centre does." · large α "Look at the purple bands: the
  wall layer is all that viscosity controls."
- **Explain:**
  0. *What the views show* — "**Left:** the <b class=c-blue>velocity profile</b> across the full diameter now (faint: the last half
     cycle). Grey dashed: the two limits. <b class=c-purple>Purple</b> bands: the Stokes layer. **Bottom:** the
     <b class=c-orange>pressure gradient</b> G cos ωt, the <b class=c-blue>centre speed</b> and the <b class=c-purple>near-wall
     speed</b> over one cycle."
  1. *Your Womersley number* — "$\alpha=a\sqrt{\omega/\nu}$ (16.12); for the artery preset a = 2 mm, ω = 2π × 66/60 = 6.91 rad/s, ν = 3.30 × 10⁻⁶ m²/s: α =
     0.002 × √(6.91/3.30 × 10⁻⁶) = $\boxed{2.89}$." (live)
  2. *What α compares* — "δ_e = √(2ν/ω) = $\boxed{0.98\ \text{mm}}$ and α = √2 a/δ_e = √2 × 2/0.98 = 2.89 (D03): the viscous layer is half the radius."
  3. *The profile at the centre* — "$u(r,t)=\dfrac{Ga^2}{i\mu\alpha^2}(1-F_1(\alpha))e^{i\omega t}$ (16.82) at r = 0: F₁ = 1/J₀(i^{3/2}α) = {…}; 1 − F₁ = {…}; times 4/(iα²): amplitude
     $\boxed{0.522}$ of Ga²/4μ, lag $\boxed{64.3^\circ}$." (live)
  4. *Against the limits* — "Quasi-steady $u_n\approx\dfrac{G_na^2}{4\mu}\left(1-\dfrac{r^2}{a^2}\right)$ (16.87): amplitude 1, lag 0. Plug $u_n\approx\dfrac{G_na^2}{i\mu\alpha_n^2}\{1-e^{\dots}\}$ (16.88): amplitude 4/α² = 0.44, lag 90°. Your α
     sits between them."
  5. *Where vessels sit* — the `vessels` reading, or the hint "On a phone the vessel axis is hidden; the chips set the same α."
  6. *Now* — "ωt = {Viz.live('wt')}: centre u = {Viz.live('uc')}, near-wall u = {Viz.live('uw')} (units of Ga²/4μ)."
  7. *Reading the current setting* — a highlighted table (vessel · α · regime) with the current row marked, and the sentence for
     the regime (as the status, one more line: "the profile looks like {Poiseuille's parabola / a parabola with lagging wings / a flat
     plug}").
- **Derivation tab:** **D03** (4) `view: 'vessels'` (phones: `profile`); step 4 `live` "α = √2 a/δ_e = …". **D13** (13) `view: 'profile'`;
  goal `set {al: 3, wt: 0}`; step 2 `play` (the profile swings at one frequency); step 9 `watch` "the profile is finite on the axis";
  step 10 `watch` "zero at both walls"; step 13 `live` "u(0) = {…} of Ga²/4μ". **D14** (6) `view: 'profile'`; step 1 `set {al: 0.5}`; step 6
  `watch` "the blue profile lies on the dashed parabola". **D15** (8) `view: 'profile'`; step 1 `set {al: 20}`; step 6 `watch` "the purple
  band: the layer of thickness √2a/α"; step 8 `live` "δ_e/a = √2/α = …". Interpret: `s => "At α = " + … + " the centre lags the gradient by " + … + "°."`
- **Code:**
  ```python
  alpha = ch16.womersley_number({{a}}, 2*np.pi*{{HR}}/60, {{nu}})   # (16.12): {{al}}
  r = np.linspace(-1, 1, 201)                                      # r/a across the full diameter
  u = ch16.womersley_velocity_nd(np.abs(r), {{wt}}, alpha)          # (16.82) / (G a^2/4 mu)
  F1_centre = ch16.womersley_F1(0.0, 1.0, alpha)                     # {{F1}}
  ch16.womersley_regime(alpha)                                       # "{{reg}}"
  ```
- **Walkthrough (8 steps):** 1. "Pulsing flow" — "Pump blood back and forth. Is it Poiseuille flow, sloshing?" `set {preset: 'artery',
  wt: 0}` · 2. "Play a cycle" — "At α ≈ 3 the profile breathes — but the centre lags the push." `play: true` · 3. "Slow it down" — "At
  α = 0.5 the profile is Poiseuille's parabola at every instant." `set {al: 0.5}`, `play: true`, `derive: {id: 'D14', step: 6}` · 4. "Speed
  it up" — "At α = 17 the core moves as a flat plug, a quarter cycle late." `set {al: 17.36}`, `play: true` · 5. "Reversal starts at the
  wall" — "Watch the purple trace: near the wall the flow turns back first." `notes: true` · 6. "What α measures" — "α = √2 × radius /
  Stokes-layer thickness." `derive: {id: 'D03', step: 4}`, `eq: 'alpha'` · 7. "One formula" — "Bessel's function joins both limits:
  $u=\dfrac{Ga^2}{i\mu\alpha^2}(1-F_1)e^{i\omega t}$ (16.82)." `eq: 'wom'`, `inspect: true`, `code: {id: 'wom', lines: [3, 4]}` · 8. "Your turn" — "Predict the lag of the centre in the aorta. Then
  pick its chip." `controls: ['al', 'vessel']`.
- **Equations:** `alpha` "Womersley number" ref 'Eq. (16.12)' $\alpha=a\sqrt{\omega/\nu}$, live · `ode` "The profile equation" ref 'Eq. (16.76)'
  $\dfrac{d^2U}{dr^2}+\dfrac1r\dfrac{dU}{dr}-\dfrac{i\omega\rho}{\mu}U=-\dfrac G\mu$ (right side corrected — note: "the book prints +G/μ") · `wom` "Womersley's profile" ref 'Eq. (16.82)' $u(r,t)=\dfrac{Ga^2}{i\mu\alpha^2}(1-F_1(\alpha))e^{i\omega t}$, $F_1=\dfrac{J_0(i^{3/2}\alpha r/a)}{J_0(i^{3/2}\alpha)}$
  (16.81), live at the inspector radius · `small` "Small α" ref 'Eq. (16.87)' $u_n\approx\dfrac{G_na^2}{4\mu}\left(1-\dfrac{r^2}{a^2}\right)$ · `large` "Large α" ref 'Eq. (16.88)'
  $u_n\approx\dfrac{G_na^2}{i\mu\alpha_n^2}\left\{1-\exp\left[-\sqrt{\dfrac{n\omega}{2\nu}}(1+i)(a-r)\right]\right\}$ (nω for harmonic n — slip S7).
- **Check yourself:** (1) "At which α does the centre lag the gradient by about 90°?" — "Large α (≳ 10): the core is a mass driven by
  G, accelerated rather than resisted; (16.88) has 1/i." `set {al: 17.36}` · (2) "Where does the flow reverse first?" — "At the wall: the
  slow layer there has the least inertia and follows the gradient's sign change sooner." · (3) "Double the heart rate in an artery
  with α = 2.9. New α?" — "2.9 × √2 = 4.1: α ∝ √ω." `set {al: 4.09}` · (4) "Why can a capillary be treated as steady flow?" — "α ≈ 0.006 ≪
  1: the Stokes layer is far wider than the vessel."
- **Selftest parity rows:** `{name: 'centre α=3', js: und(0, 0, 3), py: 'ch16.womersley_velocity_nd(0.0, 0.0, 3.0)', rtol: 1e-9}`
  (0.226425) · `{name: 'r=0.9 α=8', js: und(0.9, 1.0, 8), py: 'ch16.womersley_velocity_nd(0.9, 1.0, 8.0)', rtol: 1e-9}` (0.036880) ·
  `{name: 'α=12', js: und(0.9, 1.0, 12), py: 'ch16.womersley_velocity_nd(0.9, 1.0, 12.0)', rtol: 1e-9}` (0.021493) · `{name: 'α=20', js:
  und(0.9, 0, 20), py: 'ch16.womersley_velocity_nd(0.9, 0.0, 20.0)', rtol: 1e-8}` (0.0025328) · `{name: 'J0 series', js: J0(C(-2.1213,
  2.1213)).re, py: 'ch16.bessel_j0_series((-2.1213+2.1213j)).real', rtol: 1e-10}` · `{name: 'alpha', js: alpha(2e-3, 66, 3.5e-3/1060),
  py: 'ch16.womersley_number(2e-3, 2*np.pi*66/60, 3.5e-3/1060)', rtol: 1e-12}` (2.89358).
- **Fit plan:** 360×640: status (2 lines) · `profile` over `time` (rows 1.2 : 0.8), `vessels` hidden (the chips stay); `mode` optional.
  844×390: `profile` | `time`. Desktop: rows [1.2, 0.8] with `vessels` beside `profile`.

### E5 · womersley_flow_phase
- **Title:** "Does the flow follow the pressure?" · **Summary:** "Drag the Womersley number and watch the flow shrink and fall a
  quarter cycle behind the pressure gradient — a damper turning into a mass." · **CORE:** C10 (also N61, N62, N63, N64) ·
  **Reference:** APS (system + graph + Bode + Nyquist windows linked by one state; numbered derivation with live numbers).
- **meta:** `viz:order 5` · `viz:sections 16.3` · `viz:equations 16.90 16.91 16.92 16.93 16.95 16.96 16.97 16.99 16.100` · `viz:fluidpy
  ch16.womersley_F2 ch16.womersley_XY ch16.womersley_flow_ratio ch16.womersley_phase_lag ch16.womersley_flow_rate
  ch16.womersley_wall_shear ch16.womersley_F2_series ch16.womersley_F2_asymptote` · `viz:derivations D16 D17`.
- **Physics (JS ↔ Python):** `F2(al)` (series J₀, J₁, 100 terms, α ≤ 30) ↔ `ch16.womersley_F2`; `XY(al)` ↔ `ch16.womersley_XY`;
  `ratio(al)` ↔ `ch16.womersley_flow_ratio`; `phi(al)` ↔ `ch16.womersley_phase_lag`; `Qnd(wt, al)` = ratio·cos(ωt − φ) ↔
  `ch16.womersley_flow_rate` scaled by πGa⁴/8μ; `tau(wt, al)` = Re[F₂e^{iωt}]/|F₂(0)| … (wall shear in units of Ga/2: Re[F₂ e^{iωt}]) ↔
  `ch16.womersley_wall_shear(t, …, [0, G])` scaled; `asym(al, printed)` ↔ `ch16.womersley_F2_asymptote`; `ser(al)` ↔
  `ch16.womersley_F2_series`.
- **Views** (rows [1, 1]): 1. `time` "Gradient, flow and wall shear in time" (row 0, full width): G cos ωt (orange), Q/|Q₀| (blue), wall
  shear τ_w/(Ga/2) (rose) over two cycles, the lag of Q marked as a horizontal bar "φ = {…}° = {…} ms at 66 beats/min", the wall
  shear's lead over Q marked. 2. `bode` "Amplitude and lag against α" (row 1, left): log α 0.2…30; |Q|/|Q₀| (blue, log scale left) with the
  asymptotes 1 and 8/α² muted dashed; φ (purple, right axis 0…90°) with 0° and 90° guides; the cursor at the current α; pointer: click
  to set α. 3. `plane` "Z = 1 − F₂ in the complex plane" (row 1, right, `hidePortrait`): the locus of Z(α) from 0 to 1 (amber), the
  current Z as a phasor with X and Y legs, the small-α ghost iα²/8 and the large-α ghost 1 − 2/(i^{1/2}α) − i/α² (muted), the printed
  (16.96) ghost in rose when the toggle is on.
- **Controls:** `al` "Womersley number α" (0.2…30, log, default 8) · `wt` (transport) · `n` "Harmonic n" (1…6, `optional`; uses α_n =
  α√n) · `asy` toggle "asymptote: corrected / printed (16.96)" (default corrected).
- **Transport:** `wt` 0…4π, rate 1.5 rad/s, `end: 'loop'`.
- **Presets:** "α = 1" `set {al: 1}` · "α = 3" `set {al: 3}` · "α = 8" `set {al: 8}` · "α = 20" `set {al: 20}` · "printed (16.96)" `set {asy:
  'printed', al: 10}`.
- **Status:** φ < 30°: "damper-like: viscosity limits the flow, nearly in phase (φ = {…}°)" · 30–70°: "in between" · > 70°: "mass-like:
  inertia limits the flow, {…}° late, |Q|/|Q₀| ≈ 8/α² = {…}".
- **Readouts:** "|Q|/|Q₀|" · "Lag φ" (°) · "Lag in time" (ms at 66 beats/min) · "X, Y" · "Wall-shear lead" (°).
- **Terms** (title "Z = 1 − F₂ = X + iY", unit –): `X` "X (in quadrature with G)" (purple), `Y` "Y (in phase with G)" (blue); total
  "|Z| = √(X² + Y²)" (amber) — the bars show the legs and the hypotenuse's length.
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** · **terms** · **inspector** (click
  the Bode curve: the arithmetic of ratio and φ at that α) · **status**.
- **Explain:**
  0. *What the views show* — "**Top:** the <b class=c-orange>pressure gradient</b>, the <b class=c-blue>flow</b> it drives and the
     <b class=c-rose>wall shear</b>, each scaled to its own size. **Bottom left:** how big the flow is (blue, relative to Poiseuille's)
     and how late (purple) for every α. **Bottom right:** the complex number Z that holds both."
  1. *F₂ at your α* — "$F_2(\alpha)=\dfrac{2J_1(i^{3/2}\alpha)}{i^{3/2}\alpha J_0(i^{3/2}\alpha)}$ (16.91) at α = 8: z = i^{3/2} × 8; F₂ = $\boxed{0.177-0.161i}$." (live)
  2. *X and Y* — "$1-F_2=X+iY$ (16.97): X = $\boxed{0.823}$, Y = $\boxed{0.161}$."
  3. *How much flow* — "$|Q|=\dfrac{\pi Ga^4}{\mu\alpha^2}\sqrt{X^2+Y^2}$ (16.99) over $|Q_0|=\dfrac{\pi Ga^4}{8\mu}$ (16.95): 8√(X² + Y²)/α² = 8 × 0.838/64 = $\boxed{0.105}$ — a tenth of Poiseuille's."
  4. *How late* — "$\tan\phi=X/Y$ (16.100): φ = arctan2(0.823, 0.161) = $\boxed{78.9^\circ}$; at 66 beats/min (ω = 6.91 rad/s) that is 1.378/6.91 =
     $\boxed{199\ \text{ms}}$ after the gradient's peak."
  5. *The wall shear leads* — "$\tau_w=-\mu\dfrac{\partial u}{\partial r}\Big\rvert_a=\dfrac a2GF_2e^{i\omega t}$ (16.92, corrected, one harmonic): its phase is arg F₂ = {…}°, so it peaks {…}° before the flow.
     The book calls this a shear rate with an unnamed F (slip S8)."
  6. *The two series* — "Small α: $F_2\approx1-\dfrac{i\alpha^2}8$ (16.93) (error {…}); large α: $F_2\approx\dfrac{2}{i^{1/2}\alpha}\left(1+\dfrac{i^{3/2}}{2\alpha}\right)$ (16.96, corrected) (error {…}). With the toggle
     on 'printed' the correction 1/(2α) makes the error {…} — worse than the leading term alone (slip S10)." Hint for the plane view on
     phones.
  7. *Now* — "ωt = {Viz.live('wt')}: G = {…}, Q/|Q₀| = {…}, τ_w = {…}."
  8. *Reading the current setting* — φ < 30°: "Viscosity rules: the flow follows the push like a damper — capillaries, arterioles." ·
     30–70°: "Both matter: small arteries." · > 70°: "Inertia rules: the blood is a mass the gradient accelerates back and forth; most of
     the push buys acceleration, not flow — the aorta."
- **Derivation tab:** **D16** (7) `view: 'plane'` (phones: `bode`); step 5 `live` "F₂ = {…}"; step 7 `live` "|Q|/|Q₀| = {…}". **D17** (7)
  `view: 'time'`; step 3 `watch` "the blue curve is a shifted cosine"; step 7 `live` "φ = arctan2(X, Y) = {…}°", `set {al: 8}`. Interpret:
  `s => "At α = " + … + " the flow reaches " + … + " % of Poiseuille's, " + … + " ms late."`
- **Code:**
  ```python
  F2 = ch16.womersley_F2({{al}})                  # (16.91): {{F2}}
  X, Y = ch16.womersley_XY({{al}})                # (16.97): X = {{X}}, Y = {{Y}}
  ch16.womersley_flow_ratio({{al}})               # (16.99)/(16.95): {{ratio}}
  phi = ch16.womersley_phase_lag({{al}})          # (16.100): {{phi}} rad
  phi / (2*np.pi*66/60) * 1000                    # lag {{ms}} ms at 66 beats/min
  ```
- **Walkthrough (6 steps):** 1. "Push and flow" — "Does the flow peak when the push peaks?" `set {al: 1, wt: 0}`, `play: true` · 2. "Turn
  up α" — "At α = 8 the flow is a tenth as big and 79° late." `set {al: 8}`, `play: true` · 3. "Why late" — "Write 1 − F₂ = X + iY:
  the flow is Y cos ωt + X sin ωt — a shifted cosine." `derive: {id: 'D17', step: 5}`, `terms: true` · 4. "The Bode picture" — "Amplitude
  falls like 8/α², lag climbs to 90°: a driven mass." `eq: 'amp'`, `code: {id: 'wf', lines: [3, 4]}` · 5. "The wall feels it first" — "The rose
  wall shear peaks before the flow." `eq: 'tau'` · 6. "Your turn" — "Predict the lag at α = 20, then click the curve there." `controls:
  ['al', 'asy']`.
- **Equations:** `Q` "Flow rate" ref 'Eq. (16.90)' $Q(t)=\dfrac{\pi a^4}{i\mu\alpha^2}G(1-F_2(\alpha))e^{i\omega t}$ · `F2` "F₂" ref 'Eq. (16.91)' $F_2(\alpha)=\dfrac{2J_1(i^{3/2}\alpha)}{i^{3/2}\alpha J_0(i^{3/2}\alpha)}$, live · `amp` "Amplitude and lag" ref 'Eq. (16.99)–(16.100)'
  $|Q|=\dfrac{\pi Ga^4}{\mu\alpha^2}\sqrt{X^2+Y^2}$, $\tan\phi=\dfrac XY$, live · `tau` "Wall shear" ref 'Eq. (16.92, corrected)' $\tau_w=-\mu\dfrac{\partial u}{\partial r}\Big\rvert_a=\dfrac{G_0a}2+\dfrac a2\sum_1^NG_nF_2(\alpha_n)e^{in\omega t}$ · `asy` "Large-α F₂" ref 'Eq. (16.96, corrected)'
  $F_2\approx\dfrac{2}{i^{1/2}\alpha}\left(1+\dfrac{i^{3/2}}{2\alpha}\right)$, note "the book prints 1 + 1/(2α)".
- **Check yourself:** (1) "At α = 20, how big is the flow compared with Poiseuille's at the same gradient amplitude?" — "8|Z|/α² ≈ 0.019
  — about 2 %." `set {al: 20}` · (2) "Why can the lag never exceed 90°?" — "X and Y stay ≥ 0 (Z runs from 0 to 1 in the upper half plane),
  so arctan2(X, Y) ≤ 90°: at worst the blood behaves as a pure mass, a quarter cycle late." · (3) "Switch on the printed (16.96) at α = 10:
  is it better than the leading term?" — "No: error 0.0186 against 0.0102 for the leading term alone (and 0.00027 corrected)." `set {asy:
  'printed', al: 10}` · (4) "Which peaks first: wall shear or flow?" — "Wall shear: arg F₂ < 0 vs the flow's −φ."
- **Selftest parity rows:** `{name: 'ratio α=8', js: ratio(8), py: 'ch16.womersley_flow_ratio(8.0)', rtol: 1e-9}` (0.104797) · `{name:
  'phi α=3', js: phi(3), py: 'ch16.womersley_phase_lag(3.0)', rtol: 1e-9}` (0.947760) · `{name: 'X α=20', js: XY(20)[0], py:
  'ch16.womersley_XY(20.0)[0]', rtol: 1e-8}` (0.929266) · `{name: 'printed asymptote', js: asym(10, true).re, py:
  'ch16.womersley_F2_asymptote(10.0, printed_16_96=True).real', rtol: 1e-12}` · invariant `{name: 'small α', js: ratio(0.05), expect: 1,
  rtol: 1e-5}`.
- **Fit plan:** 360×640: status (1 line) · `time` over `bode`, `plane` hidden (X, Y as terms bars on steps 3–4; readouts). 844×390: `time`
  | `bode`. Desktop: `time` full width, `bode` | `plane` below.

### E6 · murray_branching
- **Title:** "Why isn't a wider vessel always better?" · **Summary:** "Trade pumping power against the cost of keeping blood, find
  the best radius, then drag a branch point to the cheapest place." · **CORE:** C12 (also N74, N75, N76) · **Reference:** ch14 E8
  `elliptic_loading_aspect_ratio` (one falling and one rising term with a minimum between, a ghost, a computed badge) and
  `fid_formula_lab.html` (term bars adding to a total); ch14 E9 (a force polygon that closes).
- **meta:** `viz:order 6` · `viz:sections 16.3` · `viz:equations 16.138 16.139 16.140 16.141 16.144 16.145 16.146 16.147` · `viz:fluidpy
  ch16.murray_cost ch16.murray_terms ch16.murray_radius ch16.murray_min_cost ch16.murray_daughters ch16.murray_angles
  ch16.murray_theta ch16.murray_symmetric_angle ch16.bifurcation_cost ch16.bifurcation_point` · `viz:derivations D22 D23 D24`.
- **Physics (JS ↔ Python):** `terms(a, Q, L, mu, b)` ↔ `ch16.murray_terms`; `aopt(Q, mu, b, printed)` ↔ `ch16.murray_radius(…,
  printed_16_140=…)`; `Emin` ↔ `ch16.murray_min_cost`; `daughters(a0, split)` ↔ `ch16.murray_daughters`; `angles(a0, a1, a2)` ↔
  `ch16.murray_angles`; `wcost(B)` ↔ `ch16.bifurcation_cost`; `bestB()` — Nelder–Mead in 2-D written in 30 lines (or a Weiszfeld
  iteration for the weighted Fermat point) ↔ `ch16.bifurcation_point`; `sym()` ↔ `ch16.murray_symmetric_angle`.
- **Views** (rows [1, 1]): 1. `cost` "Cost against radius" (row 0, left): a (mm, log 0.5…50) against power (mW per 10 cm, log): pumping
  (orange, falling as a⁻⁴), upkeep (teal, rising as a²), total (amber), the optimum marked with a dot and a vertical guide "a* = {…}
  mm", the reader's radius cursor (draggable) with its three values; when the "printed (16.140)" mode is on, a rose marker at the
  printed radius. 2. `fork` "The bifurcation" (row 0, right): fixed ends A (parent inlet), C, D (daughter outlets) on a 2-D canvas
  (equal: true), segments drawn with width ∝ radius (parent teal, daughters teal lighter), the draggable branch point B (blue), the
  optimal B ghost (muted ring) when "show optimum" is on, angle arcs θ and φ measured from the parent's extension. 3. `triangle`
  "Three pulls" (row 1, `hidePortrait`): the vectors a₀²e₀, a₁²e₁, a₂²e₂ drawn tip-to-tail (amber); at the optimum they close into a
  triangle, otherwise the gap is drawn rose with its length "imbalance {…}".
- **Controls:** `Q` "Parent flow $Q_0$" (1…200 mL/s, log, default 82.5) · `split` "Flow split $Q_1/Q_0$" (0.05…0.95, step 0.01, default 0.5)
  · `b` "Upkeep cost $K$ (W/m³)" (100…5000, log, default 1000) · `mu` "Viscosity μ" (2…6 mPa s, default 3.5, `optional`) · `opt` toggle
  "show the optimum B" · `formula` chips "(16.140) corrected" / "(16.140) as printed" (default corrected).
- **Presets:** "symmetric fork" `set {split: 0.5, formula: 'corrected'}` (θ = φ = 37.47°) · "70/30 split" `set {split: 0.7}` (25.7°, 49.7°)
  · "tiny side branch" `set {split: 0.95}` · "printed (16.140)" `set {formula: 'printed'}`.
- **Status:** "cursor at the optimum (within 1 %)" / "radius {±…} % from the optimum: cost +{…} %"; for the fork "B at the optimum" /
  "B off the optimum: cost +{…} %"; in printed mode "printed (16.140): K upside down — a* = {…} mm, units do not work".
- **Readouts:** "Optimal radius a*" (mm) · "Pumping : upkeep" · "Daughters a₁, a₂" (mm) · "θ, φ" (°) · "Fork cost vs optimum".
- **Terms** (title "Power at your radius", unit mW per 10 cm): `pump` "pumping 8μLQ²/(πa⁴)" (orange), `upkeep` "upkeep Kπa²L" (teal);
  total (amber). At the optimum the bars show 1 : 2.
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **presets** · **terms** · **inspector** (drag B: each L_i,
  each a_i²L_i and the sum) · **status** · **modes** (corrected / printed (16.140)).
- **Explain:**
  0. *What the views show* — "**Left:** the cost of running 10 cm of vessel: <b class=c-orange>pumping</b> falls as the vessel widens,
     <b class=c-teal>upkeep</b> grows; the <b class=c-amber>total</b> has a minimum. **Right:** a fork; drag the blue branch point.
     **Below:** the three 'pulls' a_i² e_i — at the best B they close a triangle."
  1. *The two costs at your radius* — "$\hat E_1=\dfrac{8\mu L}{\pi a^4}Q^2+K\pi a^2L$ (16.138) at a = {…} mm: pumping {…} mW, upkeep {…} mW, total {…} mW." (live)
  2. *The best radius* — "$a=\left[\dfrac{16\mu}{\pi^2K}\right]^{1/6}Q^{1/3}$ (16.140, corrected): (16 × 0.0035/(π² × 1000))^{1/6} × (8.25 × 10⁻⁵)^{1/3} = $\boxed{5.81\ \text{mm}}$. There pumping = 5.31 mW,
     upkeep = 10.62 mW: exactly 1 : 2, and $\hat E_{1,\min}=\dfrac{3\pi}2KLa^2$ (16.141) = $\boxed{15.9\ \text{mW}}$."
  3. *The printed formula* — "(16.140) as printed puts K in the numerator: it gives 58.1 mm here, and a radius that grows with the
     upkeep cost — the opposite of sense; its units are not metres (slip S15)."
  4. *The daughters* — "$Q=Q_1+Q_2\Rightarrow a_0^3=a_1^3+a_2^3$ (16.144): split {…} → a₁ = {…} mm, a₂ = {…} mm."
  5. *The angles* — "$\cos\theta=\dfrac{a_0^4+a_1^4-a_2^4}{2a_0^2a_1^2}$ (16.146): θ = {…}°, φ = {…}°; symmetric fork arccos 2^{−1/3} = $\boxed{37.47^\circ}$ (ours)."
  6. *Your branch point* — "Weighted length a₀²L₀ + a₁²L₁ + a₂²L₂ = {…} against {…} at the optimum: +{…} %." Hint on phones: "The
     triangle of pulls is hidden; its imbalance is in the status."
  7. *Reading the current setting* — at the optimum: "Every segment pays twice as much to keep its blood as to push it — the
     signature of Murray's optimum." · cursor narrower: "Too narrow: pumping dominates (a⁻⁴)." · wider: "Too wide: upkeep
     dominates (a²)." · fork: "Thick daughters stay near the parent's line; thin ones turn away — the weights are a²."
- **Derivation tab:** **D22** (8) `view: 'cost'`; step 4 `watch` "the amber curve is flat at its minimum"; step 6 `live` "a* = … mm"; step 8
  `terms: true`. **D23** (3) `view: 'fork'`; step 3 `live` "a₁³ + a₂³ = … = a₀³". **D24** (11) `view: 'triangle'` (phones: `fork`); step 3 `set
  {opt: false}`, `watch` "drag B toward C: L₁ shrinks, L₀ grows"; step 5 `set {opt: true}`, `watch` "the triangle closes at the
  optimum"; step 11 `set {split: 0.5}`, `live` "θ = arccos 2^{−1/3} = 37.47°". Interpret: `s => "With a " + … + " split the thinner daughter
  turns " + … + "°, the thicker " + … + "°."`
- **Code:**
  ```python
  a_star = ch16.murray_radius({{Q}}e-6, {{mu}}e-3, {{b}})          # (16.140): {{astar}} m
  ch16.murray_terms(a_star, {{Q}}e-6, 0.1, {{mu}}e-3, {{b}})       # pumping {{pump}}, upkeep {{up}} W
  a1, a2 = ch16.murray_daughters(a_star, {{split}})               # (16.144): {{a1}}, {{a2}} m
  ch16.murray_angles(a_star, a1, a2)                               # (16.146): θ {{th}}, φ {{ph}} rad
  ch16.bifurcation_point(A, C, D, a_star, a1, a2)["B"]             # numeric optimum of B
  ```
- **Walkthrough (7 steps):** 1. "How wide should a vessel be?" — "Wider means easier pumping — so why aren't vessels huge?" `set
  {preset: 'symmetric fork'}` · 2. "Two costs" — "Pumping falls as a⁻⁴; keeping blood costs ∝ a². Drag the radius." `terms: true`,
  `eq: 'cost'` · 3. "The minimum" — "Zero slope gives $a=[16\mu/(\pi^2K)]^{1/6}Q^{1/3}$ (16.140, corrected): Q ∝ a³." `derive: {id: 'D22', step: 6}`, `code: {id: 'mur', lines: [1, 2]}` · 4. "A
  telltale ratio" — "At the optimum pumping is exactly half the upkeep." `eq: 'min'` · 5. "Forks" — "Every segment optimal ⇒ a₀³ = a₁³ +
  a₂³ (16.144)." `derive: {id: 'D23', step: 3}` · 6. "Where to branch" — "Drag B. The cheapest place is where the three pulls a_i² e_i
  close a triangle." `inspect: true`, `derive: {id: 'D24', step: 5}` · 7. "Your turn" — "Predict the angles of a 90/10 split. Then drag the
  split." `controls: ['split', 'Q', 'b']`.
- **Equations:** `cost` "Cost of a segment" ref 'Eq. (16.138)' $\hat E_1=\dfrac{8\mu L}{\pi a^4}Q^2+K\pi a^2L$, live · `opt` "Optimal radius" ref 'Eq. (16.140, corrected)'
  $a=\left[\dfrac{16\mu}{\pi^2K}\right]^{1/6}Q^{1/3}$, note "the book prints K in the numerator" · `min` "Minimum cost" ref 'Eq. (16.141)' $\hat E_{1,\min}=\dfrac{3\pi}2KLa^2$ · `cube` "Cube law" ref 'Eq. (16.144)'
  $a_0^3=a_1^3+a_2^3$ · `ang` "Branch angles" ref 'Eq. (16.146)' $\cos\theta=\dfrac{a_0^4+a_1^4-a_2^4}{2a_0^2a_1^2}$, live.
- **Check yourself:** (1) "Double the flow. By what factor does the best radius grow?" — "2^{1/3} = 1.26: Q ∝ a³." `set {Q: 165}` · (2) "At the
  optimum, what fraction of the cost is pumping?" — "One third: pumping = ½ upkeep (16.141)." · (3) "Why does the thin daughter of a 90/10
  split turn more?" — "Its pull a₂² is small, so balancing the triangle needs a large angle." `set {split: 0.9}` · (4) "Switch to the
  printed (16.140). What goes wrong when K rises?" — "The printed radius grows with K — expensive blood would make vessels wider."
  `set {formula: 'printed'}`.
- **Selftest parity rows:** `{name: 'a*', js: aopt(8.25e-5, 3.5e-3, 1000, false), py: 'ch16.murray_radius(8.25e-5, 3.5e-3, 1000.0)',
  rtol: 1e-12}` (5.81388e-3) · `{name: 'printed', js: aopt(8.25e-5, 3.5e-3, 1000, true), py: 'ch16.murray_radius(8.25e-5, 3.5e-3, 1000.0,
  printed_16_140=True)', rtol: 1e-12}` (5.81388e-2) · `{name: 'theta 70/30', js: angles(1, 0.7**(1/3), 0.3**(1/3)).theta, py:
  'ch16.murray_angles(1.0, 0.7**(1/3), 0.3**(1/3))["theta"]', rtol: 1e-12}` (0.448647) · `{name: 'sym', js: sym(), py:
  'ch16.murray_symmetric_angle()', rtol: 1e-12}` (0.653928) · invariant `{name: 'pump/upkeep', js: (t => t.pumping/t.upkeep)(terms(aopt(…),
  …)), expect: 0.5, rtol: 1e-12}`.
- **Fit plan:** 360×640: status (1 line) · `cost` over `fork` (rows 1 : 1), `triangle` hidden (its imbalance in the status); `mu`
  optional. 844×390: `cost` | `fork`. Desktop: rows [1, 0.8], `triangle` under `fork`.

### E7 · junction_reflection
- **Title:** "When does a fork send the pulse back?" · **Summary:** "Send a pressure pulse into a branching artery and watch it split,
  reflect and flip sign as the daughters' admittance changes — the energy always balances." · **CORE:** C13 (also R10, N78, N79, N80,
  N81) · **Reference:** ch15 E1 `sound_pulse_dalembert` (shapes sliding on a strip) and `fid_formula_lab.html` (terms adding to a
  total).
- **meta:** `viz:order 7` · `viz:sections 16.3` · `viz:equations 16.148 16.151 16.152 16.153 16.154 16.159 16.160 16.161 16.162 16.163` ·
  `viz:fluidpy ch16.tube_admittance ch16.junction_coefficients ch16.junction_energy ch16.junction_waveforms ch16.junction_from_geometry
  ch16.moens_korteweg_speed` · `viz:derivations D25 D26`.
- **Physics (JS ↔ Python):** `Y(a, c, rho)` = πa²/(ρc) ↔ `ch16.tube_admittance`; `RT(Y1, Ys)` ↔ `ch16.junction_coefficients`; `energy`
  ↔ `ch16.junction_energy`; `geom(ap, ad1, ad2, cp, cd)` ↔ `ch16.junction_from_geometry`; the pulses in closed form: parent incident
  f(t − x/c₁) for x < 0, reflected R f(t + x/c₁), daughters T f(t − x/c_j) for x > 0 (Gaussian time shape) — the exact solution of the
  linear model, no grid ↔ `ch16.junction_waveforms` for the parent; probe p and Q (16.159)–(16.160).
- **Views** (rows [1, 0.8]): 1. `net` "The fork" (row 0, full width): parent from the left (x from −1.2 m to 0), two daughters leaving at
  ±30° (lengths 1 m), tube widths ∝ radius (teal), the pressure drawn as a coloured band inside each tube (orange positive, purple
  negative) and as a small height profile along each tube's axis; the probe marker in the parent at x = −0.6 m; labels "incident",
  "reflected", "transmitted". 2. `probe` "At the probe in the parent" (row 1, left): p(t) (orange) and Q(t) (blue) normalised by the
  incident peak, showing the incident pulse then the reflected one — p and Q equal at first, then opposite or equal in sign depending
  on R. 3. `bars` "Energy shares" (row 1, right, `hidePortrait`): incident (1, amber outline), reflected R² (rose), transmitted ΣY_jT²/Y₁ (blue),
  stacked to 1; plus R and T as signed bars.
- **Controls:** `r` "Daughter radius / parent radius" (0.3…1.2, step 0.01, default 0.8) · `s` "Daughter wave speed / parent" (0.5…2,
  step 0.01, default 1.0; stiffer daughters are faster) · `r2` "Second daughter radius (if different)" (0.3…1.2, `optional`, default =
  r) · `w` "Pulse width" (0.05…0.3 s, `optional`, default 0.1) · `t` (transport).
- **Transport:** `t` 0…0.6 s (time for the pulse to reach the junction and the reflection to pass the probe at c₁ = 5 m/s), rate
  0.15 s/s, `end: 'hold'`; end-of-run card: "R = {…}: the reflection came back {inverted / upright}; energy: {100R²} % reflected, {…} %
  transmitted."
- **Presets:** "matched (ΣY = Y₁)" `set {r: 0.7071, s: 1, r2: 0.7071}` · "closed-end-like" `set {r: 0.4, s: 1, r2: 0.4}` (R = +0.515) ·
  "open-end-like" `set {r: 1.0, s: 1, r2: 1.0}` (R = −0.333) · "Murray-optimal fork" `set {r: 0.7937, s: 1, r2: 0.7937}` (R = −0.115).
- **Status:** |R| < 0.02 "matched — no reflection" · R > 0 "closed-end-like: the pressure peak grows at the fork (R = +{…})" · R < 0
  "open-end-like: the reflection comes back as a dip (R = {…})".
- **Readouts:** "ΣY_j / Y₁" · "R (pressure reflected)" · "T (pressure transmitted)" · "Energy reflected R²" · "Energy sum".
- **Terms** (title "Energy carried", unit share of incident): `refl` "reflected R²" (rose), `trans` "transmitted (ΣY_j/Y₁)T²" (blue);
  total "= 1" (amber, flat).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** with an end-of-run card · **presets** ·
  **terms** · **status**.
- **Explain:**
  0. *What the views show* — "**Top:** a parent artery splitting in two; the colour inside each tube is the
     <b class=c-orange>pressure</b> (purple when below the resting value). **Bottom left:** pressure and <b class=c-blue>flow</b> at a
     point in the parent. **Bottom right:** where the pulse's energy goes."
  1. *Each tube's admittance* — "$Q=Au=\pm Yp$ (16.151) with Y = A/(ρc): Y₁ = πa₁²/(ρc₁); each daughter Y_j = Y₁ r²/s = $\boxed{0.64\,Y_1}$; sum ΣY_j/Y₁ = $\boxed{1.28}$." (live)
  2. *Reflection and transmission* — "$\mathrm R=\dfrac{Y_1-\sum Y_j}{Y_1+\sum Y_j}$ (16.153) = (1 − 1.28)/(1 + 1.28) = $\boxed{-0.123}$; T = 1 + R = $\boxed{0.877}$."
  3. *Energy* — "reflected $\mathrm R^2\dfrac{p_I^2}{Z_1}$ (16.162): R² = 0.015; transmitted (16.163): 1.28 × 0.877² = 0.985; sum $\boxed{1.000}$ (ours: the junction neither
     stores nor loses energy)." Hint on phones for the bars.
  4. *At the probe* — "$Q=Y_1P_I\left[f\left(t-\dfrac x{c_1}\right)-\dfrac{P_R}{P_I}f\left(t+\dfrac x{c_1}\right)\right]$ (16.160): the incident pulse has p and Q in step; the reflected one has p = R p_I, Q = −R Y₁ p_I —
     here {opposite / the same} signs. Now: p = {Viz.live('pp')}, Q = {Viz.live('pq')}."
  5. *Reading the current setting* — R > 0.05: "Closed-end-like (ΣY_j < Y₁): the fork is harder to push into than the parent; the
     pressure rises there and part of the pulse comes back upright. This is how reflections from the leg arteries raise the
     aortic pressure peak." · |R| ≤ 0.05: "Matched: the fork passes the pulse on — Murray's forks are close." · R < −0.05:
     "Open-end-like: the fork is easier; the reflection returns as a pressure dip, flow is drawn on."
- **Derivation tab:** **D25** (6) `view: 'net'`; step 2 `watch` "the flow into the fork equals the flow out"; step 5 `live` "R = (1 − {…})/(1 +
  {…}) = {…}". **D26** (6) `view: 'bars'` (phones: `net` with the terms bars); step 6 `live` "R² + (ΣY/Y₁)T² = {…} + {…} = 1". Interpret: `s => "Your
  fork reflects " + … + " % of the energy."`
- **Code:**
  ```python
  Y1 = ch16.tube_admittance(np.pi*a1**2, 1060.0, c1)        # (16.151)
  Yd = [Y1*{{r}}**2/{{s}}, Y1*{{r2}}**2/{{s}}]              # each daughter
  R, T = ch16.junction_coefficients(Y1, Yd)                 # (16.153): R = {{R}}, T = {{T}}
  e = ch16.junction_energy(Y1, Yd)                          # (16.162)-(16.163)
  e["reflected"], e["transmitted"], e["total"]              # {{er}}, {{et}}, {{etot}}
  ```
- **Walkthrough (6 steps):** 1. "A pulse meets a fork" — "Some goes on, some comes back. What decides how much?" `set {preset:
  'Murray-optimal fork', t: 0}` · 2. "Play" — "A small dip comes back: the daughters together are slightly easier to fill." `play: true` ·
  3. "Admittance" — "Each tube passes flow Q = Yp, Y = A/(ρc). Pressure is shared at the fork, flow is conserved." `eq: 'Y'`, `derive: {id:
  'D25', step: 2}` · 4. "Closed or open" — "Narrow daughters reflect upright (R > 0), wide ones inverted (R < 0)." `set {preset:
  'closed-end-like'}`, `play: true` · 5. "Energy balances" — "R² plus the transmitted share is always 1." `terms: true`, `derive: {id:
  'D26', step: 6}` · 6. "Your turn" — "Find the radius ratio that makes the fork invisible. Predict first." `controls: ['r', 's']`.
- **Equations:** `cont` "At the junction" ref 'Eq. (16.148)' $p_I+p_R=p_{T1}=p_{T2}$ with $Y_1(p_I-p_R)=\sum_jY_jp_{Tj}$ (16.152) · `Y` "Admittance" ref 'Eq. (16.151)' $Q=Au=\pm Yp=\pm\dfrac pZ$ ·
  `RT` "Reflection and transmission" ref 'Eq. (16.153)' $\mathrm R=\dfrac{Y_1-\sum Y_j}{Y_1+\sum Y_j}$, $\mathrm T=\dfrac{2Y_1}{Y_1+\sum Y_j}$, live · `En` "Energy" ref 'Eq. (16.162)–(16.163)' reflected
  $\mathrm R^2\dfrac{p_I^2}{Z_1}$, transmitted share $\dfrac{Z_2^{-1}+Z_3^{-1}}{Z_1^{-1}}\mathrm T^2$, with $\mathrm R^2+\dfrac{\sum Y_j}{Y_1}\mathrm T^2=1$ (ours).
- **Check yourself:** (1) "Which daughter radius ratio (equal speeds) matches the fork?" — "r = 1/√2 = 0.707: 2r² = 1." `set {r: 0.7071,
  r2: 0.7071}` · (2) "If the daughters are stiffer (faster), does the fork reflect more or less upright?" — "More upright: Y = A/(ρc) falls
  with c, so ΣY drops below Y₁." `set {s: 1.5}` · (3) "With R = +0.5 the transmitted pressure is 1.5 times the incident. Where does the
  extra energy come from?" — "Nowhere: the transmitted energy is (ΣY/Y₁)T² = 0.75 (D26); the flow is smaller." `set {r: 0.4, r2: 0.4}`.
- **Selftest parity rows:** `{name: 'R 0.8', js: RT(1, [0.64, 0.64])[0], py: 'ch16.junction_coefficients(1.0, [0.64, 0.64])[0]', rtol:
  1e-14}` (−0.122807) · `{name: 'T murray', js: geom(1, 0.7937005, 0.7937005, 5, 5).T, py: 'ch16.junction_from_geometry(1.0, [2**(-1/3),
  2**(-1/3)], 5.0, [5.0, 5.0], 1060.0)["T"]', rtol: 1e-6}` (0.884987) · `{name: 'energy', js: energy(1, [0.3, 0.5]).total, py:
  'ch16.junction_energy(1.0, [0.3, 0.5])["total"]', rtol: 1e-14}` (1) · invariant `{name: 'T = 1 + R', js: (x => x[1] - x[0])(RT(1, [0.2,
  0.9])), expect: 1, rtol: 1e-14}`.
- **Fit plan:** 360×640: status (1 line) · `net` over `probe`, `bars` hidden (terms bars on steps 5; readouts R, T, R²). 844×390: `net` |
  `probe`. Desktop: rows [1, 0.8].

### E8 · dean_curved_pipe
- **Title:** "What does a bend do to the flow?" · **Summary:** "Follow blood round a curved tube: two vortices appear, the fast core
  slides to the outer wall, and the flux drops only at second order." · **CORE:** C14 (also N83, N84, N87, N90, N91, N92, N93) ·
  **Reference:** AFE (linked views on one state) and `ddpm3_unet_3d.html` (an orbitable three.js scene with a fallback note).
- **meta:** `viz:order 8` · `viz:sections 16.3` · `viz:equations 16.164 16.171 16.172 16.177 16.178 16.179 16.181 16.182` · `viz:fluidpy
  ch16.dean_secondary_streamfunction ch16.dean_axial_correction ch16.dean_axial_velocity ch16.dean_secondary_velocity ch16.dean_flux_ratio
  ch16.curved_resistance_ratio ch16.dean_number ch16.dean_convert ch16.dean_series_reliability ch16.curvature_ratio` · `viz:derivations D27`.
- **Physics (JS ↔ Python):** `psi0(r, al)` ↔ `ch16.dean_secondary_streamfunction`; `w1(r, al)` (plus sign) ↔ `ch16.dean_axial_correction`;
  `w(r, al, K)` ↔ `ch16.dean_axial_velocity`; `uv(r, al, K)` ↔ `ch16.dean_secondary_velocity`; `flux(K)` (exact 1541/50400, printed 0.0120) ↔
  `ch16.dean_flux_ratio(K, exact=True)`; `lam(K)` ↔ `ch16.curved_resistance_ratio`; `conv(v, from, to)` ↔ `ch16.dean_convert`;
  `reliable(K)` (threshold 0.25, ours) ↔ `ch16.dean_series_reliability(K, threshold=0.25)`; tracer helices: integrate (u, v, w) in the
  torus with RK4 (`Viz.num.rk4Step`), the secondary speed exaggerated ×{N} for visibility (label).
- **Views** (rows [1]): 1. `section` "The cross-section" (left): the unit disc with the outer wall (α = 0) on the right, labelled "outer
  wall (away from the bend's centre)"; heatmap of w = w₀ + κ²w₁ (blue scale), ψ₀ streamlines (purple, 8 levels), arrows of (u, v) on a
  coarse grid, the maximum of w marked "w_max {…} at x = {…}a". 2. `torus` "Round the bend" (right, 3-D, `Viz.three`): a quarter-torus of
  the curvature ratio δ, 12 tracers spiralling (blue), the secondary motion exaggerated (label "×{N}"); fallback without WebGL: a 2-D
  view looking down the axis with the tracers' cross-section paths. 3. `curve` "Flux against Dean number" (`hidePortrait`, below on
  desktop): Q/(πa²W̄) against K from 0 to 900 (amber), λ (rose) on the same axes, the cursor, a shaded band "series unreliable (ours:
  K⁴ term > ¼ of K² term)" beyond K ≈ 460.
- **Controls:** `K` "Dean number" (0…900, step 5, default 288) with `conv` chips "κ" / "K" / "D" (the slider shows the chosen convention's
  value; the state stores K) · `delta` "Curvature ratio δ = a/R" (0.05…0.5, step 0.01, default 0.25) · `tr` toggle "tracers" (default on).
- **Presets:** "straight pipe" `set {K: 0}` · "gentle bend" `set {K: 100}` · "within the series" `set {K: 400}` · "beyond the series"
  `set {K: 800}`.
- **Status:** `reliable(K)`: "series reliable (K⁴ term {…} % of the K² term)" · "beyond the small-K series (ours: K > {K_max}) — the real
  flow needs numerics".
- **Readouts:** "κ, K, D" · "Flux ratio Q/(πa²W̄)" · "Resistance ratio λ" · "Outward shift of w_max" (a) · "Mean speed for this K" (mm/s, with
  a = 12.5 mm, ν = 3.3 × 10⁻⁶).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **3-D view** · **presets** · **inspector** (click a point of
  `section`: r, α, ψ₀, w₀, κ²w₁, w with the arithmetic) · **status**.
- **Explain:**
  0. *What the views show* — "**Left:** a slice across the tube, the bend's centre to the left. Colour: the <b class=c-blue>axial
     speed</b>; <b class=c-purple>purple</b> curves: the secondary flow's streamlines. **Right:** blood going round the bend; the
     sideways motion is exaggerated. **Below:** how the flux and resistance change with the Dean number."
  1. *Your Dean number* — "$\delta=a/R$ (16.164) = 0.25; $\kappa=\sqrt\delta\,2Re$ (16.171) with Re = aW̄/ν (slip S18 explains which Re); K = 2κ² = {…}, D = 4√2 κ = {…}
     (K = D²/16). For a = 12.5 mm, ν = 3.3 × 10⁻⁶ m²/s that needs a mean speed of only {…} mm/s — the series is a slow-flow theory."
  2. *The vortex pair* — "ψ₀ = r(1 − r²)²(4 − r²) sin α/144 (ours, D27): strongest at r = 0.43, ψ₀,max = 0.0076; scaled by κ: the vortex
     carries κψ₀,max = {…} of a·W̄."
  3. *The skewed axial flow* — "w = w₀ + κ²w₁ with $w_1=+\dfrac{r(1-r^2)(19-21r^2+9r^4-r^6)}{5760}\cos\alpha$ (ours): the maximum moves {…}a toward the outer wall."
  4. *Flux and resistance* — "$\dfrac{Q}{\pi a^2\bar W}=1-0.0306\left(\dfrac K{576}\right)^2+0.0120\left(\dfrac K{576}\right)^4$ (16.178) = $\boxed{…}$; $\lambda=1+0.0306\left(\dfrac K{576}\right)^2-0.0110\left(\dfrac K{576}\right)^4$ (16.182) = $\boxed{…}$. 0.0306 = 1541/50400 (D27)."
  5. *Round the bend* — the 3-D view's reading or the fallback note; hint on phones for the curve view.
  6. *Reading the current setting* — K = 0: "Straight pipe: Poiseuille." · reliable: "The bend adds a gentle vortex pair; the flux
     is only {…} % lower — second order." · beyond: "The small-K series has lost its meaning here (the K⁴ term rivals the K² term);
     real arch flows (K far larger) form thin boundary layers and can show four vortices — named, not computed here."
- **Derivation tab:** **D27** (14) `view: 'section'`; goal `set {K: 0}`; step 2 `watch` "a centred, symmetric bullseye"; step 4 `set {K:
  288}`; step 8 `watch` "two purple cells"; step 11 `watch` "the blue maximum slides right"; step 14 `live` "1 − 0.030575 (K/576)² = {…}". Interpret:
  `s => "At K = " + … + " the flux is " + … + " % below the straight pipe's."`
- **Code:**
  ```python
  K = {{K}}                                                      # Dean's K = 2 kappa^2
  r, al = np.meshgrid(np.linspace(0, 1, 60), np.linspace(0, 2*np.pi, 120))
  w = ch16.dean_axial_velocity(r, al, K)                         # w0 + kappa^2 w1 (ours)
  psi = ch16.dean_secondary_streamfunction(r, al)                # psi0 (ours)
  ch16.dean_flux_ratio(K, exact=True)                            # (16.178): {{flux}}
  ch16.curved_resistance_ratio(K)                                # (16.182): {{lam}}
  ```
- **Walkthrough (6 steps):** 1. "Round a bend" — "Blood rounds the aortic arch. Does it just turn?" `set {preset: 'straight pipe'}` ·
  2. "Flung outward" — "The fast core feels the strongest centrifugal push and slides outward." `set {K: 288}`, `eq: 'series'` · 3. "Two
  vortices" — "It returns along the walls: a counter-rotating pair, ψ₀ ∝ sin α." `derive: {id: 'D27', step: 8}`, `inspect: true` · 4. "Helices" — "Follow
  the tracers: corkscrews round the bend." `play: true` · 5. "A small price" — "The flux falls only as K²: 0.0306 (K/576)²." `derive: {id:
  'D27', step: 14}`, `code: {id: 'dean', lines: [5, 6]}` · 6. "Your turn" — "Predict the flux at K = 576. Where does the series stop
  being trustworthy?" `controls: ['K', 'conv']`.
- **Equations:** `delta` "Curvature ratio" ref 'Eq. (16.164)' $\delta=\dfrac aR$ · `kappa` "Dean number" ref 'Eq. (16.171)' $\kappa=\sqrt\delta\,\dfrac{2a\bar W_0}\nu=\sqrt\delta\,2Re$, conventions note ·
  `series` "Dean's expansion" ref 'Eq. (16.177)' $w=\sum_{n=0}^\infty\kappa^{2n}w_n,\ \psi=\kappa\sum_{n=0}^\infty\kappa^{2n}\psi_n$ · `flux` "Flux ratio" ref 'Eq. (16.178)' $\dfrac{Q}{\pi a^2\bar W}=1-0.0306\left(\dfrac K{576}\right)^2+0.0120\left(\dfrac K{576}\right)^4$, live · `lam`
  "Resistance ratio" ref 'Eq. (16.182)' $\lambda=1+0.0306\left(\dfrac K{576}\right)^2-0.0110\left(\dfrac K{576}\right)^4$, live.
- **Check yourself:** (1) "Which way does the axial maximum move — toward the bend's centre or away?" — "Away, to the outer wall: the
  centrifugal term pushes the fast core outward and w₁ > 0 there." · (2) "Double K from 288 to 576. How much does the flux deficit
  grow?" — "About 4×: it is ∝ K² (0.69 % → 1.86 %)." `set {K: 576}` · (3) "In D convention the slider shows 96 at K = 576. Why?" — "K =
  D²/16: D = √(16 × 576) = 96." `set {conv: 'D', K: 576}`.
- **Selftest parity rows:** `{name: 'psi0', js: psi0(0.5, Math.PI/2), py: 'ch16.dean_secondary_streamfunction(0.5, np.pi/2)', rtol:
  1e-14}` (0.00732422) · `{name: 'w1 sign', js: w1(0.5, 0), py: 'ch16.dean_axial_correction(0.5, 0.0)', rtol: 1e-14}` (+0.000930786) ·
  `{name: 'flux', js: flux(288), py: 'ch16.dean_flux_ratio(288.0, exact=True)', rtol: 1e-12}` (0.993106) · `{name: 'D from K', js:
  conv(576, 'K', 'D'), py: 'ch16.dean_convert(576.0, from_="K", to="D")', rtol: 1e-12}` (96) · invariant `{name: 'w0 mean', js: meanW(0),
  expect: 1, rtol: 1e-4}` (Simpson over the disc).
- **Fit plan:** 360×640: status (1 line) · `section` over `torus` (rows 1 : 0.9; `torus` falls back to 2-D if WebGL is missing), `curve`
  hidden (flux and λ readouts). 844×390: `section` | `torus`. Desktop: `section` | `torus` over `curve`.

### E9 · collapsible_tube_choking
- **Title:** "Why can't you suck harder through a floppy straw?" · **Summary:** "Lower the outlet pressure of a collapsible tube:
  the flow rises, then stops listening — the outlet has reached its own wave speed, S = 1." · **CORE:** C15 (also N96, N97, N98, N99,
  N100, N101, R13) · **Reference:** ch15 E2 `area_mach_nozzle` (term bars through a vanishing 1 − M²) and ch15 E4 `nozzle_back_pressure`
  (a transport on the downstream pressure, a plateau).
- **meta:** `viz:order 9` · `viz:sections 16.3` · `viz:equations 16.186 16.188 16.189 16.195 16.203 16.204 16.205 16.206` · `viz:fluidpy
  ch16.collapsible_tube_law ch16.collapsible_wave_speed ch16.speed_index ch16.shapiro_M ch16.choking_flux ch16.collapsible_terms
  ch16.collapsible_march ch16.starling_state ch16.flow_limitation_curve` · `viz:derivations D28 D29 D30`.
- **Physics (JS ↔ Python):** `P(al, n)` ↔ `ch16.collapsible_tube_law`; `c(al, Kp, rho, n)` ↔ `ch16.collapsible_wave_speed`; `S(u, c)` ↔
  `ch16.speed_index`; `M(n)` = 2 − n ↔ `ch16.shapiro_M`; `Qstar(al, A0, Kp, rho, n)` ↔ `ch16.choking_flux`; `terms(al, Q, prm, drivers)` ↔
  `ch16.collapsible_terms`; the march: integrate dx/dα (Simpson, 400 panels, C.12) from α_in down to the outlet ↔ `ch16.collapsible_march`;
  `starling(p2)` (brentq on Q for p(L) = p₂, plateau at Q_choke) ↔ `ch16.starling_state`; `curve()` ↔ `ch16.flow_limitation_curve`.
- **Views** (rows [1, 1]): 1. `tube` "The tube in its chamber" (row 0, full width): a horizontal tube of length L in a box at pressure
  p_e, its height ∝ α(x) (teal), with five cross-section icons along x (circle → ellipse → two lobes for α < 0.3), the fill coloured by
  S (blue → amber at S = 1), inlet p₁ and outlet p₂ gauges (orange); in "throat" mode the reference area A₀(x) drawn as a muted
  outline; in "jugular" mode the tube is drawn vertical with "gravity ↓". 2. `along` "α(x) and S(x)" (row 1, left): α (teal) and S (blue)
  against x, the line S = 1 (amber dashed), the choke point marked. 3. `side` (row 1, right, `hidePortrait`): mode-dependent — "Starling":
  Q against p₂ − p_e (blue), the plateau and its onset p₂* marked, the current point; "throat"/"jugular": term bars of (16.204) at the
  cursor station x_c.
- **Controls:** `p2` "Outlet pressure $p_2-p_e$" (−20…0 mmHg, step 0.1, default −3) · `mu` "Viscosity (friction)" (1…8 mPa s, default 3.5) ·
  `n` "Tube-law exponent n" (1.2…2.0, step 0.05, default 1.5, `optional`) · `Kp` "Bending stiffness $K_p$" (100…500 Pa, default 200, `optional`)
  · `mode` modes: `starling` "Starling resistor" · `throat` "Smooth throat" · `jugular` "Giraffe jugular" · `Q` "Flow Q (throat, jugular)"
  (1…40 mL/s, default 15; hidden in Starling mode where Q is computed) · `xc` "Cursor station" (0…1 of L, `optional`).
- **Presets:** "open tube" `set {mode: 'starling', p2: -0.5}` · "onset of limitation" `set {mode: 'starling', p2: -7.1}` · "deep
  limitation" `set {mode: 'starling', p2: -15}` · "smooth throat at S = 1" with a getter `get set()` returning `{mode: 'throat', Q: Qcrit(prm)}` (Q_crit = the
  throat's choking flux, computed — the distinguishing value is in `set`) · "giraffe jugular (dz/dx = −1)" `set {mode: 'jugular', Q: 10}`.
- **Status:** S_out < 0.98: "S < 1 at the outlet ({…}): the flow feels the downstream pressure" · limited: "S = 1 at the outlet: flow limited
  at Q = {…} mL/s — lowering p₂ further changes nothing upstream" · throat mode: "S passes 1 at the throat (dA₀/dx = 0)" · jugular:
  "gravity balances friction: the vein runs partly collapsed at α = {…}".
- **Readouts:** "Flow Q" (mL/s) · "Outlet area ratio α" · "Outlet speed index S" · "Wave speed at the outlet c" (m/s) · "M (16.206)".
- **Terms** (throat and jugular modes; title "Drivers of (16.204) at x_c", unit 1/m): `area` "S²A₀′/A₀" (teal), `ext` "−p_e′/ρc²" (orange),
  `grav` "−g z′/c²" (amber), `fric` "−R′Q/ρc²" (rose), `stiff` "−P K_p′/ρc²" (purple); total "driver" — and a separate badge "1 − S² = {…}",
  "dα/dx = α·driver/(1 − S²)".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (three systems) · **terms** · **presets** · **status**.
- **Explain:**
  0. *What the views show* — "**Top:** a soft tube in a pressure chamber; its height is the <b class=c-teal>area ratio α = A/A₀</b>, its
     colour the <b class=c-blue>speed index S = u/c</b> (amber at S = 1). **Bottom left:** α and S along the tube. **Bottom right:** flow
     against outlet pressure (Starling mode) or the drivers of the area change (other modes)."
  1. *The tube law and its wave speed* — "$-\mathrm P\approx\alpha^{-n}-1$ (16.188), P = (p − p_e)/K_p: at the outlet α = {…} so p − p_e = {…} Pa. $c^2=\dfrac{nK_p\alpha^{-n}}{\rho}$ (16.189):
     c = √(1.5 × 200 × α^{−1.5}/1060) = $\boxed{…\ \text{m/s}}$."
  2. *Speed and speed index* — "u = Q/(αA₀) = {…} m/s; $S=u/c$ (16.195) = $\boxed{…}$."
  3. *Why friction collapses the tube* — "$(1-S^2)\dfrac1\alpha\dfrac{d\alpha}{dx}=\dfrac{S^2}{A_0}\dfrac{dA_0}{dx}-\dfrac1{\rho c^2}\left[\dfrac{dp_e}{dx}+\rho g\dfrac{dz}{dx}+\mathrm RQ+\mathrm P\dfrac{dK_p}{dx}\right]$ (16.204) with only friction: dα/dx = −α R′Q/(ρc²(1 − S²)) < 0 while S < 1; as α
     falls, S rises (S² ∝ α^{n−2}); at S = 1 the slope blows up — the tube chokes."
  4. *The limit* — "With L = 0.1 m and an open inlet the outlet reaches S = 1 at Q = $\boxed{19.99\ \text{mL/s}}$, α = 0.312, p₂ − p_e = −946 Pa (−7.1
     mmHg) (ours, model). Lower p₂ further: a short collapsed segment forms past the choke point, but the flow upstream cannot hear
     it."
  5. *M and the gas-dynamics twin* — "$\mathrm M=3+\dfrac{\alpha\,\partial^2\mathrm P/\partial\alpha^2}{\partial\mathrm P/\partial\alpha}$ (16.206) = 2 − n = $\boxed{0.5}$ for the power law; S plays the Mach number, the tube law the equation of
     state, the throat the Laval throat (ch15)."
  6. *The side view* — reading of the curve or the bars; hint on phones.
  7. *Reading the current setting* — S_out < 0.5: "Subcritical and far from the limit: Q follows p₁ − p₂ almost like a resistor." ·
     0.5 ≤ S_out < 1: "Close to the limit: the outlet narrows fast; small changes of p₂ still matter." · limited: "Flow limitation (Starling
     Case 1): the veins of a raised arm, the airways in a forced cough and an artery under a blood-pressure cuff all behave like
     this." · jugular: "Gravity pulls the blood down as fast as friction allows; the vein stays partly collapsed, at S < 1 here."
- **Derivation tab:** **D28** (4) `view: 'along'`; step 4 `live` "c = √(nK_pα^{−n}/ρ) = …". **D29** (9) `view: 'along'`; step 8 `watch` "nothing
  here depends on x through P, so the S23 term is zero". **D30** (15) `view: 'side'` (phones: `along`); goal `set {mode: 'throat'}`; step 6
  `terms: true`; step 7 `live` "(1 − S²) = …, driver = …"; step 15 `set {preset: 'smooth throat at S = 1'}`, `watch` "S crosses 1 exactly
  where A₀ is smallest". Interpret: `s => "Your outlet runs at S = " + … + (… ? ": flow limited." : ": below the limit.")`
- **Code:**
  ```python
  prm = dict(L=0.1, A0=np.pi*0.004**2, Kp={{Kp}}, n={{n}}, rho=1060.0, mu={{mu}}e-3)
  st = ch16.starling_state(0.0, {{p2}}*MMHG, 0.0, prm)        # Q {{Q}} m^3/s, S_out {{S}}
  c = ch16.collapsible_wave_speed(st["alpha_out"], prm["Kp"], 1060.0, n=prm["n"])   # (16.189)
  ch16.flow_limitation_curve(np.linspace(-20, 0, 41)*MMHG, p1=0.0, pe=0.0, prm=prm)["Q_choke"]
  ch16.shapiro_M(st["alpha_out"], n=prm["n"])                  # (16.206): {{M}}
  ```
- **Walkthrough (7 steps):** 1. "Suck harder?" — "Lower the pressure at the end of a soft tube. Does the flow keep rising?" `set
  {preset: 'open tube'}` · 2. "It rises, then stalls" — "Drag p₂ down: the flow climbs and then flattens at about 20 mL/s." `set {p2:
  -12}`, `controls: ['p2']` · 3. "A slow wave" — "A floppy tube carries waves at under 2 m/s: $c^2=nK_p\alpha^{-n}/\rho$ (16.189)." `derive: {id: 'D28', step: 4}`, `eq: 'c'` · 4. "The speed index" —
  "S = u/c plays the Mach number. At S = 1 the outlet chokes." `eq: 'S'` · 5. "The 1 − S² factor" — "Friction narrows the tube, which
  raises S, which narrows it faster — until 1 − S² = 0." `set {mode: 'throat'}`, `terms: true`, `derive: {id: 'D30', step: 7}` · 6. "Crossing
  smoothly" — "Only at a throat can S pass through 1." `set {preset: 'smooth throat at S = 1'}` · 7. "Your turn" — "Predict: does a more
  viscous fluid limit at a higher or lower flow? Then drag μ." `controls: ['mu', 'p2']`.
- **Equations:** `law` "Collapsible tube law" ref 'Eq. (16.188)' $-\mathrm P\approx\alpha^{-n}-1$, $\mathrm P=\dfrac{p-p_e}{K_p}$, $\alpha=\dfrac A{A_0}$ (16.186) · `c` "Wave speed" ref 'Eq. (16.189)' $c^2=\dfrac{nK_p\alpha^{-n}}\rho$, live · `S`
  "Speed index" ref 'Eq. (16.195)' $S=\dfrac uc$, live · `sh` "Shapiro's area equation" ref 'Eq. (16.204)' $(1-S^2)\dfrac1\alpha\dfrac{d\alpha}{dx}=\dfrac{S^2}{A_0}\dfrac{dA_0}{dx}-\dfrac1{\rho c^2}\left[\dfrac{dp_e}{dx}+\rho g\dfrac{dz}{dx}+\mathrm RQ+\mathrm P\dfrac{dK_p}{dx}\right]$ ·
  `M` "Shapiro's M" ref 'Eq. (16.206)' $\mathrm M=3+\dfrac{\alpha\,\partial^2\mathrm P/\partial\alpha^2}{\partial\mathrm P/\partial\alpha}$, live.
- **Check yourself:** (1) "Below which outlet pressure does lowering p₂ stop increasing the flow?" — "About −7.1 mmHg here: there the outlet
  reaches S = 1." `set {preset: 'onset of limitation'}` · (2) "Why can S only pass through 1 where A₀ has a minimum (with no other
  driver)?" — "At S = 1 the left side of (16.204) is zero, so its right side S²A₀′/A₀ must vanish: A₀′ = 0." `set {mode: 'throat'}` · (3) "A
  collapsed vein (α = 0.25) carries a wave at what speed with K_p = 250 Pa?" — "c = √(1.5 × 250 × 8/1000) = 1.73 m/s (16.189)." ·
  (4) "Is P positive or negative in a collapsed tube?" — "Negative: −P = α^{−n} − 1 > 0 for α < 1 (T9)."
- **Selftest parity rows:** `{name: 'c', js: c(0.5, 200, 1060, 1.5), py: 'ch16.collapsible_wave_speed(0.5, 200.0, 1060.0, 1.5)', rtol:
  1e-12}` (0.894706) · `{name: 'Q*', js: Qstar(0.5, Math.PI*16e-6, 200, 1060, 1.5), py: 'ch16.choking_flux(0.5, np.pi*0.004**2, 200.0,
  1060.0, 1.5)', rtol: 1e-12}` (2.24864e-5) · `{name: 'Q choke', js: qChoke(prm), py: 'ch16.flow_limitation_curve([-15*133.322387415],
  p1=0.0, pe=0.0, prm=dict(L=0.1, A0=np.pi*0.004**2, Kp=200.0, n=1.5, rho=1060.0, mu=3.5e-3))["Q_choke"]', rtol: 5e-3}` (design
  1.99899e-5 by α-integration; the module's bisection returned 2.00290e-5 when probed — 0.2 % apart, an open item for the verifier;
  tighten the rtol to 1e-4 once both use the same criterion) ·
  `{name: 'M', js: M(1.5), py: 'ch16.shapiro_M(0.5, 1.5)', rtol: 1e-9}` (0.5) · invariant `{name: 'P<0', js: P(0.5, 1.5) < 0 ? 1 : 0, expect:
  1, rtol: 0}`.
- **Fit plan:** 360×640: status (2 lines) · `tube` over `along`, `side` hidden (Q and the plateau in readouts; terms bars on steps 5–6).
  844×390: `tube` over `along` beside the panel. Desktop: rows [0.9, 1].

### E10 · xylem_water_potential
- **Title:** "How does water climb a 30-metre tree?" · **Summary:** "Follow water from soil to leaf through pressure, solutes and
  height: the sap climbs by being pulled, under tension." · **CORE:** C17 (also N110, N111, N112, N113) · **Reference:**
  `fid_formula_lab.html` (terms that add to a total, stage chips) and FDV (Explain panel).
- **meta:** `viz:order 10` · `viz:sections 16.4` · `viz:equations 16.17 16.221 16.222` · `viz:fluidpy ch16.water_potential
  ch16.water_potential_terms ch16.osmotic_pressure ch16.xylem_column ch16.xylem_reynolds ch16.cavitation_margin
  ch16.sieve_tube_pressure_drop ch16.poiseuille_flow_rate ch16.capillary_rise` · `viz:derivations D34`.
- **Physics (JS ↔ Python):** `psi(p, c, T, z)` ↔ `ch16.water_potential`; `terms(…)` ↔ `ch16.water_potential_terms`; `Pi(c, T)` ↔
  `ch16.osmotic_pressure`; `column(z)` → {p, psi, gravity, viscous} ↔ `ch16.xylem_column(z, *, psi_soil, u, a, mu)`; `Re(u, a)` ↔
  `ch16.xylem_reynolds`; `margin(p)` with p_cav = −2.0 MPa (illustrative, ours) ↔ `ch16.cavitation_margin`; `sieve(…)` ↔
  `ch16.sieve_tube_pressure_drop`; `cap(a)` = 2σ/(ρga) (full wetting) ↔ `ch16.capillary_rise(0.0728, π/2, 1000, a)`.
- **Views** (rows [1]): 1. `tree` "Soil to leaf" (left): a tree sketch (trunk, crown, roots) with the xylem column coloured by its
  pressure (orange positive, purple negative), four stations (soil, root, stem top, leaf) each with three small stacked bars: pressure
  (orange), −R_uTc (teal), ρgz (amber), and the total ψ as a black tick; the capillary-only height marked as a muted line on the trunk
  ("capillarity alone: 0.59 m"). In phloem mode the same space shows a source leaf → sink root tube with sieve plates. 2. `profile` "ψ
  and p against height" (right): z (0…H) vertical; p (orange) and ψ (black), the hydrostatic line ψ_soil − ρgz as a muted ghost, the extra
  viscous drop shaded rose; the illustrative cavitation limit −2 MPa as a dashed line; pointer: click a height for the inspector.
  3. `phloem` "Sieve tube" (`hidePortrait`, mode phloem only): the share of the pressure drop in the plates against N.
- **Controls:** `H` "Tree height" (0.5…120 m, log, default 30) · `u` "Sap speed in a conduit" (0…3 mm/s, step 0.05, default 1.0) · `a`
  "Conduit radius" (5…100 µm, log, default 25) · `ps` "Soil water potential" (−2…0 MPa, step 0.05, default −0.1) · `mode` "xylem / phloem".
- **Presets:** "houseplant" `set {H: 0.5, u: 0.5}` · "30 m tree (illustrative)" `set {H: 30, u: 1.0, a: 25, ps: -0.1}` · "100 m giant" `set {H:
  100, u: 1.0}` · "drought" `set {ps: -1.5}`.
- **Status:** margin > 0.3 MPa: "column holds (top pressure {…} MPa)" · 0 < margin ≤ 0.3: "close to the illustrative cavitation limit" ·
  margin ≤ 0: "tension beyond the illustrative limit (−2 MPa, ours) — the column would break (embolism)".
- **Readouts:** "Top pressure p" (MPa) · "ψ at the top" (MPa) · "Viscous share of the drop" (%) · "Re in a conduit" · "Capillary rise" (m).
- **Terms** (title "ψ at the station you pick", unit MPa): `p` "pressure p" (orange), `osm` "solutes −R_uTc" (teal), `grav` "height ρgz" (amber);
  total "ψ" (black).
- **Depth features:** Explain + Code + Derivation · **terms** · **presets** · **status** · **inspector** (click a height: p, ψ and each term
  with the arithmetic) · **modes** (xylem / phloem) · **linked views**.
- **Explain:**
  0. *What the views show* — "**Left:** a tree; the trunk's colour is the sap's <b class=c-orange>pressure</b> (purple below zero = tension).
     At four stations: <b class=c-orange>pressure</b>, <b class=c-teal>solute</b> and <b class=c-amber>height</b> terms of the water
     potential. **Right:** ψ and p against height."
  1. *The ledger* — "$\psi=p-RTc+\rho gz$ (16.221): water moves from high ψ to low ψ. In the soil ψ = {…} MPa; in a leaf cell with c = 0.3 mol/L the solute term
     is −R_uTc = −8.314 × 293 × 300 = $\boxed{-0.731\ \text{MPa}}$."
  2. *Holding a column* — "Still column, pure water: p(z) = ψ_soil − ρgz; at the top {H} m: −0.1 − 1000 × 9.81 × 30/10⁶ = $\boxed{-0.394\ \text{MPa}}$ (D34)."
  3. *Moving the sap* — "Poiseuille in a conduit, $Q=\dfrac{\pi a^4}{8\mu}\dfrac{\Delta p}L$ (16.17), as a gradient: 8μu/a² = 8 × 10⁻³ × 10⁻³/(25 × 10⁻⁶)² = $\boxed{12.8\ \text{kPa/m}}$; over 30 m −0.384
     MPa more: top pressure $\boxed{-0.778\ \text{MPa}}$; Re = u·2a/ν = $\boxed{0.05}$ — creeping flow."
  4. *Why capillarity is not enough* — "Ch. 1's rise h = 2σ sin α/(ρgR) for a 25 µm conduit: 0.59 m (ours). The tree is fifty times
     taller."
  5. *The phloem (mode)* — "$\Delta p=\dfrac{8\mu Q}{\pi a^4}\left[L+\dfrac\ell N\left(\dfrac ar\right)^4\right]+2\,\Delta p'$ (16.222): for our sieve tube the plates take $\boxed{85\ \%}$ of the drop (illustrative sizes)."
     Hint on phones.
  6. *Reading the current setting* — margin large: "The column hangs comfortably: evaporation pulls, cohesion holds." · near: "The
     sap is stretched close to breaking: trees in drought close their stomata to stop this." · beyond: "The column would cavitate;
     tall trees survive because their conduits are narrow and many."
- **Derivation tab:** **D34** (5) `view: 'profile'`; step 3 `set {u: 0}`, `watch` "the black ψ line is vertical: no flow, no gradient"; step 4
  `live` "p(z) = ψ_soil − ρgz = …"; step 5 `set {preset: '30 m tree'}`. Interpret: `s => "At the top of your " + … + " m tree the sap is at " + … + " MPa."`
- **Code:**
  ```python
  col = ch16.xylem_column(np.linspace(0, {{H}}, 61), psi_soil={{ps}}e6, u={{u}}e-3,
                          a={{a}}e-6, mu=1.0e-3)               # (16.221) + Poiseuille (ours)
  col["p"][-1]/1e6                                             # top pressure {{ptop}} MPa
  ch16.osmotic_pressure(300.0, 293.15)/1e6                      # solutes {{osm}} MPa
  ch16.xylem_reynolds({{u}}e-3, {{a}}e-6, 1.0e-6)               # Re {{Re}}
  ch16.capillary_rise(0.0728, np.pi/2, 1000.0, {{a}}e-6)        # Ch. 1: {{cap}} m
  ```
- **Walkthrough (6 steps):** 1. "Up a tree" — "Thirty metres of water hangs in a trunk. What holds it up?" `set {preset: '30 m tree'}`
  · 2. "One ledger" — "Pressure, solutes and height add up to the water potential ψ; water flows downhill in ψ." `eq: 'psi'`, `terms: true`
  · 3. "Negative pressure" — "Still column: p falls by ρg per metre — below zero. The sap is under tension." `set {u: 0}`, `derive: {id:
  'D34', step: 5}` · 4. "Moving sap costs more" — "Flow through narrow conduits adds a Poiseuille drop." `set {u: 1.0}`, `inspect: true`,
  `code: {id: 'xy', lines: [1, 2]}` · 5. "Not capillarity" — "A 25 µm tube lifts water only 0.6 m." `eq: 'cap'` · 6. "Your turn" — "How tall can
  a tree be before the illustrative limit? Predict, then drag H." `controls: ['H', 'u', 'ps']`.
- **Equations:** `psi` "Water potential" ref 'Eq. (16.221)' $\psi=p-RTc+\rho gz$, live at the selected station · `pois` "Poiseuille" ref 'Eq. (16.17)'
  $Q=\dfrac{\pi a^4}{8\mu}\dfrac{\Delta p}L$ · `sieve` "Sieve tube" ref 'Eq. (16.222)' $\Delta p=\dfrac{8\mu Q}{\pi a^4}\left[L+\dfrac\ell N\left(\dfrac ar\right)^4\right]+2\,\Delta p'$ · `cap` "Capillary rise (Ch. 1, Example 1.1)" ref 'Example
  1.1' $h=\dfrac{2\sigma\sin\alpha}{\rho gR}$.
- **Check yourself:** (1) "In a still column, how much does the pressure drop per 10 m?" — "ρg × 10 = 0.098 MPa ≈ 1 bar." · (2) "Why is ψ in the
  xylem below the soil's when sap flows?" — "Water moves down ψ; the viscous drop needs a gradient beyond hydrostatics." `set {u: 2}` · (3)
  "Halve the conduit radius. What happens to the viscous drop?" — "×4: the gradient 8μu/a² ∝ 1/a² at fixed speed." `set {a: 12.5}`.
- **Selftest parity rows:** `{name: 'psi', js: psi(-0.394e6, 0, 293.15, 30), py: 'ch16.water_potential(-0.394e6, 0.0, 293.15, 30.0)',
  rtol: 1e-12}` · `{name: 'osm', js: Pi(300, 293.15), py: 'ch16.osmotic_pressure(300.0, 293.15)', rtol: 1e-12}` (731215) · `{name: 'top p',
  js: column(30).p, py: 'ch16.xylem_column(30.0, psi_soil=-0.1e6, u=1e-3, a=25e-6, mu=1e-3)["p"]', rtol: 1e-10}` (−7.7820e5 Pa) · `{name:
  'sieve share', js: sieve(…).share, py: '(lambda d: (d["plate_viscous"] + d["plate_entry"])/d["total"])(ch16.sieve_tube_pressure_drop(np.pi*1e-10*0.3e-3,
  1.5e-3, 10e-6, 1e-6, 300e-6, 1e-6, 50))', rtol: 1e-10}` (0.84705) · invariant `{name: 'still column ψ', js: column(20, 0).psi/(-0.1e6), expect: 1, rtol: 1e-12}` (no flow: ψ is the soil's at every
  height).
- **Fit plan:** 360×640: status (1 line) · `tree` (55 %) and `profile` (45 %) side by side in one row (portrait: stacked, `tree` on top);
  `phloem` only in phloem mode and hidden on phones (its share in the readout). 844×390: one row. Desktop: one row.

### B1 · fahraeus_lindqvist_layer
(backup — built only if one of E1–E10 fails review)
- **Title:** "Why does blood seem thinner in narrow tubes?" · **Summary:** "A thin plasma layer lubricates the wall; the narrower the
  tube, the bigger its share." · **CORE:** C05 (also N29, N30, N33–N36) · **Reference:** FDV (Explain panel); `overfitting_curves.html`.
- **meta:** `viz:order 11` (backup) · `viz:sections 16.2 16.3` · `viz:equations 16.8 16.27 16.28 16.31 16.33 16.34 16.35` · `viz:fluidpy
  ch16.two_layer_velocity ch16.two_layer_flow_rates ch16.apparent_viscosity_two_layer ch16.apparent_viscosity_thin_layer
  ch16.charm_kurland_core_viscosity ch16.pries_relative_viscosity` · `viz:derivations D05 D06 D07`.
- **Views:** `profile` (two-layer profile across the diameter, core blue, plasma layer purple, kink marked; Poiseuille ghost with μ_c),
  `curve` (μ_app/μ_p against D from 5 to 1000 µm, log: our model (blue, "model, ours") beside the cited Pries curve (muted, "cited")),
  `section` (`hidePortrait`: the cross-section with the cell-free ring drawn to scale).
- **Controls:** `delta` "Layer thickness δ" (0.5…5 µm, default 2) · `ratio` "μ_c/μ_p" (1…6, default 3.33) or chips "from hematocrit
  (Charm–Kurland)" · `D` "Tube diameter" (10…1000 µm, log, default 40).
- **Presets:** "no layer (δ → 0)", "40 µm tube (illustrative)", "large artery (D = 1 mm)", "equal viscosities".
- **Status:** "μ_app = {…} μ_p — {…} % below the core's viscosity" / "δ ≥ a/2: the model no longer applies".
- **Depth features:** Explain (0 views; 1 the profile's two parabolas (16.27)–(16.28); 2 the fluxes (16.31); 3 μ_app (16.33) and the thin form
  (16.34); 4 against the cited curve; 5 reading: "the effect needs a less viscous layer at the wall; its size depends on δ, which the book
  leaves free (T19)") · Code · Derivation (D05, D06, D07; `view: 'profile'`, D07 `view: 'curve'`) · **linked views** · **presets** · **status**.
- **Walkthrough (5 steps):** the puzzle; the layer; the profile with a kink (`derive: {id: 'D05', step: 9}`); μ_app falls as D shrinks
  (`derive: {id: 'D07', step: 3}`); your turn.
- **Equations:** (16.28), (16.31), (16.33), (16.34) written out as in Part F. **Check yourself (3):** what happens at δ → 0; why μ_app < μ_c
  needs μ_c > μ_p; halve D at fixed δ.
- **Selftest:** `{name: 'mu_app', js: muApp(20e-6, 2e-6, 1.2e-3, 4e-3), py: 'ch16.apparent_viscosity_two_layer(20e-6, 2e-6, 1.2e-3,
  4e-3)', rtol: 1e-12}` (2.2192e-3) · `{name: 'thin', js: …, py: 'ch16.apparent_viscosity_thin_layer(20e-6, 2e-6, 1.2e-3, 4e-3)', rtol:
  1e-12}` (2.0690e-3).
- **Fit plan:** as E1 (two views on phones, the curve hidden; its value in the status).

---

## Part D — runtime budget (full run < 5 min on a laptop / Colab CPU; target 70–100 s)

Per-cell costs are estimates: the eight `check_src` listings were timed for this design (0.3–1.5 s each, ≈ 6 s together); the scratch
numbers of this design (Windkessel closed form, Bessel ratios, elastic-tube roots, the Starling march by α-integration) ran in under
2 s; the rest from analysis §4. The builder reports the measured numbers.

| Section | What costs time (per cell) | Full | FAST |
|---|---|---|---|
| front matter, §16.1 | imports (sympy, scipy, plotly, pandas); conventions and slips tables | 5 s | 5 s |
| §16.2 C01, C02 | loop figure, shoelace; loop animation (60 frames, FAST 30; video) 4 s / 2 s; two bar figures; **F1** (25 × 3 traces × 200 pts) 1 s | 9 s | 7 s |
| §16.3 C03, C04 | scales table, two figures | 2 s | 2 s |
| C05 | quad checks; Pries curve; `fahraeus_lindqvist_model` (brentq) for 40 diameters 1 s; **F2** (19 steps × 40 diameters × brentq) 2 s | 4 s | 3 s |
| C06 | figures; **F3** (20 × 2 × 100) < 1 s | 2 s | 2 s |
| C07 | `windkessel_ode` (solve_ivp rtol 1e-10, 8 beats) 0.3 s; RK4 loop 0.2 s; **A1** (80 frames, FAST 40) 5 s / 3 s; live cell (not run on the page); E2 iframe | 7 s | 5 s |
| C08 | `tube1d_run` (400 cells, FAST 200) 3 s / 1 s; **A2** (60, FAST 30) 4 s / 2 s; figures; E3 iframe | 9 s | 5 s |
| C09 | **D13** check 0.3 s; series J₀ checks; `womersley_pde_solve` (nr = 200, 400 steps × 8 periods; FAST nr = 80, 4 periods) 2 s / 0.5 s, cached (P252); harmonics; **A3** (72, FAST 36) 5 s / 3 s; E4 iframe | 10 s | 5 s |
| C10 | slope cells; Bode figure (300 α); finite-difference wall shear; E5 iframe | 2 s | 2 s |
| C11 | **D19**, **D20**, **D21** checks (0.4 + 0.6 + 1.0 s); `elastic_tube_determinant_numeric` at 4 α (≈ 0.1 s each); figure (200 α, both roots); **F4** (9 × 4 × 200) 1 s; **A5** (60, FAST 30) 4 s / 2 s | 9 s | 6 s |
| C12 | grid search (20 001 points); **D24** check (Nelder–Mead) 0.6 s; `bifurcation_point` × 3; **F5** (25 × 3 × 200) 1 s; E6 iframe | 3 s | 3 s |
| C13 | energy identity on 1000 random sets; **A4** (`network_pulse_run`, 40 frames, FAST 24) 3 s / 2 s; E7 iframe | 4 s | 3 s |
| C14 | **D27** check 1.2 s (`dean_series_sympy` cached with `functools.lru_cache` + `np.savez` in FAST); 2-D trapezoid; **F6** (25 × contour traces on 60 × 120 grid) 2 s; E8 iframe | 5 s | 4 s |
| C15 | **D30** check 0.3 s; `collapsible_differentials_sympy` 1 s; `flow_limitation_curve` (31 outlet pressures, each a brentq on Q with an α-integration; FAST 15) 3 s / 1.5 s; Euler march (20 000 steps) 0.3 s; E9 iframe | 6 s | 4 s |
| C16 | **D33** check 0.4 s; trapezoids; **F7** (20 × 2 × 200) < 1 s; E1 iframe | 2 s | 2 s |
| §16.4 C17 | columns; **F8** (13 × 4 × 61) < 1 s; E10 iframe | 2 s | 2 s |
| explainers | 10 `show_viz` cells (iframes; no computation) | 3 s | 3 s |
| **Total** | | **≈ 84 s** | **≈ 63 s** |

**FAST plan.** `FAST = setup_notebook()`; animations ≤ 80 frames at dpi 80 (FAST half): A1, A2, A3, A5 `player="video"` (smooth motion);
A4 `player="frames"` (the reader should stop when the pulse reaches the fork and when the echo passes the probe). Slider figures ≤ 25
steps × ≤ 4 traces × ≤ 400 points (≈ 100–200 kB each; eight figures F1–F8). Grids: Womersley profiles 201 radii; the Crank–Nicolson
check nr = 200 (FAST 80) — cached; the 1-D tube 400 cells (FAST 200); the Dean cross-section 60 × 120; the flow-limitation curve 31 points
(FAST 15). The eight ★★★ `check_src` cells always run live (≈ 6 s together). The two slow helpers (`dean_series_sympy`,
`womersley_pde_solve`) are wrapped in caches (P252 reminder); the cache files go to `outputs/ch16/` (git-ignored); nothing is
committed. Page budget: five animations < 1.2 MB each, eight plotly figures ≈ 150 kB each, ten explainers ≈ 250 kB each: ≈ 10.5 MB < 15 MB.

---

## Part E — prerequisite ledger
Every concept, symbol, maths tool and Python function or idiom the notebook or its explainers use, with where it is explained.
"primer (in Cxx)" = a 📎 primer placed in that block before first use (the Concept text begins with the primer term, used verbatim in
`nb.primer`); "knowledge/primers.md: <term> (chNN Pnn) — reminder" = a one-line reminder naming the earlier primer; a CORE/RECAP id
alone = taught there; "Cxx (Nnn)" = the NOTE placed in that block; "Cxx (D0n)" = the derivation where it is built; "Cxx (gloss …)" =
one sentence where it is used; "Cxx (🔁 Ch. n)" = an id-less recap cell in that block. New primers P374–P394 (21) in first-use order:
P378, P394, P381 (C01) · P383 (C02) · P379, P380 (C06) · P382 (C07) · P374, P377, P375 (C09) · P376 (C10) · P385, P386 (C11) · P387, P388
(C12) · P389, P391, P390 (C14) · P392 (C15) · P393 (C16) · P384 (C17).

| Concept | First used in | Explained by |
|---|---|---|
| the conventions box: two time conventions, G = −dp/dx, τ a magnitude, which Re and St, which Dean number, secant or tangent, k_c or K_c, P negative in collapse | C01 | front matter (A.0 row 5, before C01), repeated per block |
| the overloaded letters α, K, R, δ, ψ, c, T, M, S, P, X, Y, κ, λ, E, H, h, η, ν̂ | C01 | front matter (A.0 row 5), and where each first appears |
| slips kept apart from traps; `ch16.book_slips()`, `ch16.traps()`, `ch16.conventions_table()`, `ch16.illustrative_inputs()` | C01 | front matter (A.0 row 5) |
| "ours", "illustrative", "model", "cited" labels; the literal input dicts | C01 | front matter (A.0 rows 4–5) |
| biofluid mechanics; scales from cells to organs | C01 | C01 (N01) |
| systemic and pulmonary circulation; a closed system | C01 | C01 (N02) |
| heart chambers and valves; systole and diastole; the four phases of a beat | C01 | C01 (N04) |
| cardiac cycle and heart rate | C01 | C01 (N03) |
| pressure–volume loop; EDV, ESV; counterclockwise orientation | C01 | C01 |
| line integral along a path | C01 | knowledge/primers.md: line integral along a path (ch01 P35) — reminder |
| loop area ½∮x × dx | C01 | knowledge/primers.md: loop vector area (ch05 P138) — reminder |
| the shoelace formula; `np.roll` | C01 | primer (in C01) |
| units of physiology: mmHg, cP, poise, dyn/cm², bar, PRU; `ch16.mmHg_to_Pa`, `cP_to_Pa_s`, `dyn_cm2_to_Pa`, `bar_to_Pa`, `PRU`, `MMHG` | C01 | primer (in C01) |
| pint quantities `UN.Q_`, `.to()`, `UN.check_dimensions` | C01 | knowledge/primers.md: pint quantities (ch01 P02) — reminder |
| p dV work; work per beat W = −∮p dV; `ch16.pv_loop_work`, `pv_loop_schematic`, `pv_loop_schematic_work`, `cardiac_power` | C01 | C01 (D01) |
| stroke volume, ejection fraction, cardiac output, cardiac index; `ch16.stroke_volume`, `ejection_fraction`, `cardiac_output`, `cardiac_index` | C01 | C01 (N06, N07) |
| mechanical efficiency of the heart (named) | C01 | C01 (N05) |
| systolic, diastolic and pulse pressure; `ch16.pulse_pressure` | C01 | C01 (N08, N09) |
| the electrical analogy (Δp ↔ voltage, Q ↔ current, resistance, capacitor, admittance) | C01 | primer (in C01) |
| resistances in series and in parallel; `ch16.series_pressure_drop`, `network_resistance` | C01 | C01 (N10; P381) |
| mean speed from continuity; `ch16.mean_speed`, `ch16.VESSELS` | C01 | C01 (N11; 🔁 Ch. 4 mass conservation) |
| total peripheral resistance; venous pressure taken as zero; the ⅓ rule; `ch16.peripheral_resistance`, `mean_arterial_pressure` | C01 | C01 (N12, N13, N14) |
| matplotlib figures, `np.linspace`, f-strings | C01 | knowledge/primers.md: matplotlib figures (ch01 P01) — reminder |
| `animate` / `show_animation`, `slider_figure`, `show_viz`, `live` | C01 | knowledge/primers.md: animate and show_animation, slider_figure, show_viz (ch01 P16, P17, P18) — reminder |
| `assert np.allclose` / `np.isclose` | C01 | knowledge/primers.md: assert np.allclose (ch01 P15) — reminder |
| blood's properties: density, plasma and high-shear viscosity; `ch16.BLOOD` | C02 | C02 (N15) |
| red cells; equal-volume sphere area; `ch16.equal_volume_sphere_area` | C02 | C02 (N16) |
| hematocrit (fraction against percent); `ch16._check_hematocrit` | C02 | C02 (N17) |
| non-Newtonian fluid; the corrected cP chain; `ch16.cP_to_Pa_s`, `poise_to_Pa_s` | C02 | C02 (N18) |
| Newton's law of viscosity τ = μ du/dy; shear rate; `ch01.newton_shear_stress`, `ch16.newtonian_stress` | C02 | R01 |
| yield stress; the correlation (16.5); `ch16.yield_stress_correlation` | C02 | C02 (N19) |
| rouleaux and shear thinning (named) | C02 | C02 (N20) |
| ordinary derivative as a slope | C02 | knowledge/primers.md: ordinary derivative as a slope (ch01 P19) — reminder |
| secant slope and tangent slope | C02 | primer (in C02) |
| apparent viscosity (secant) and differential viscosity (tangent); `ch16.apparent_viscosity`, `differential_viscosity`, `casson_apparent_viscosity` | C02 | C02 (N21, D02) |
| the Casson law; k_c and K_c; `ch16.casson_stress`, `casson_shear_rate`, `whitmore_casson` | C02 | C02 (D02, N22) |
| power laws and log–log axes | C02 | knowledge/primers.md: power laws and log–log plots (ch01 P13) — reminder |
| `np.where` for a piecewise law | C02 | knowledge/primers.md: np.where and np.select (ch01 P46) — reminder |
| homogeneity of blood (named) | C02 | C02 (N23) |
| Womersley number α; Stokes number; `ch16.womersley_number`, `stokes_number`, `womersley_from_re_st`, `womersley_regime` | C03 | C03 (D03) |
| Stokes-layer thickness δ_e = √(2ν/ω); `LAM.stokes_layer` | C03 | C03 (🔁 Ch. 8 Stokes' second problem) |
| separation in an adverse pressure gradient | C03 | R02 |
| Reynolds and Strouhal numbers on the diameter; `SIM.reynolds_number`, `SIM.strouhal_number` | C03 | C03 (D03; trap T6 table) |
| vessel walls and sizes; the scales table; `ch16.vessel_scales_table` | C03 | C03 (N24) |
| compliance, pulse wave, reflections (named ahead) | C03 | C03 (N25) |
| Hagen–Poiseuille assumptions and Re | C04 | R03 |
| the Poiseuille profile in Δp; `ch16.poiseuille_velocity`, `LAM.pipe_poiseuille` | C04 | R04, R05 |
| centreline speed, flux, mean speed; `ch16.poiseuille_flow_rate`, `LAM.pipe_flow_rate` | C04 | R06 |
| wall shear stress as a positive magnitude; `ch16.poiseuille_wall_stress`, `LAM.pipe_wall_stress` | C04 | R07 |
| entrance length (sketch); the correlations (16.37); `ch16.entrance_length` | C04 | R08; C04 (N27) |
| definite integral; a ring of area 2πr dr | C04 | knowledge/primers.md: definite integral (ch01 P27) — reminder |
| vascular resistance R ∝ a⁻⁴; `ch16.vascular_resistance` | C04 | C04 (D04) |
| `np.trapezoid` | C04 | C04 (gloss in the from-scratch row; P27) |
| pulsatile entrance length; muscular arteries (named) | C04 | C04 (N26, N28) |
| Fåhræus–Lindqvist and Fåhræus effects; cell-free layer; `ch16.pries_relative_viscosity` | C05 | C05 (N29) |
| conservation of red cells and volume; `ch16.cell_conservation` | C05 | C05 (N30) |
| the two-layer model and its interface conditions | C05 | C05 (N31, N32) |
| free body of a fluid cylinder | C05 | C05 (D05 step 1; 🔁 Ch. 8 D07) |
| two-layer profiles; `ch16.two_layer_velocity` | C05 | C05 (D05) |
| substitution in an integral | C05 | knowledge/primers.md: substitution in an integral (ch03 P106) — reminder |
| layer fluxes; `ch16.two_layer_flow_rates` | C05 | C05 (D06) |
| `scipy.integrate.quad` | C05 | knowledge/primers.md: scipy.integrate.quad (ch03 P87) — reminder |
| sympy symbols, `diff`, `integrate`, `simplify`, `solve`, `expand` | C05 | knowledge/primers.md: sympy (ch01 P40) — reminder |
| apparent viscosity by the Poiseuille analogy; thin-layer limit; `ch16.apparent_viscosity_two_layer`, `apparent_viscosity_thin_layer` | C05 | C05 (D07) |
| first-order Taylor expansion and the binomial approximation | C05 | knowledge/primers.md: first-order Taylor expansion (ch01 P26) — reminder |
| Charm–Kurland core viscosity; `ch16.charm_kurland_core_viscosity` | C05 | C05 (N35) |
| the FL closure (ours); `ch16.fahraeus_lindqvist_model` | C05 | C05 (N36) |
| `scipy.optimize.brentq` | C05 | knowledge/primers.md: scipy.optimize.brentq (ch03 P108) — reminder |
| Young's modulus and Poisson's ratio; strain; Hooke's law | C06 | primer (in C06) |
| a stretched curved membrane (Laplace's law T/a) | C06 | primer (in C06) |
| hoop stress; tube law; transmural pressure; `ch16.hoop_stress`, `tube_radius_hookean` | C06 | C06 (D08) |
| `try` / `except ValueError` | C06 | C06 (gloss in the code row) |
| local Poiseuille in a slowly varying tube (lubrication idea); `ch16.elastic_tube_pressure` | C06 | C06 (N38; 🔁 Ch. 8 lubrication) |
| the distensible tube's flux; `ch16.elastic_tube_flow_rate` | C06 | C06 (N39) |
| pulmonary circulation; resistance falling with pressure; `ch16.elastic_tube_resistance` | C06 | C06 (N40, N41, N42) |
| the Windkessel story | C07 | C07 (N43) |
| control-volume volume balance | C07 | C07 (N44; 🔁 Ch. 4) |
| compliance K = dv/dp; `ch16.compliance_from_tube` | C07 | C07 (N45) |
| chain rule | C07 | knowledge/primers.md: chain rule (ch01 P49) — reminder |
| the two-element Windkessel ODE; `ch16.windkessel_rhs`, `windkessel_ode` | C07 | C07 (D09) |
| first-order linear ODE and the integrating factor | C07 | knowledge/primers.md: first-order linear ODE and the integrating factor (ch09 P210) — reminder |
| product rule | C07 | knowledge/primers.md: product rule for differentials (ch01 P38) — reminder |
| a sum of decaying kicks (convolution) | C07 | primer (in C07) |
| the Windkessel solution, diastolic decay, periodic state, mean = RQ̄; `ch16.windkessel_pressure`, `windkessel_periodic_p0`, `windkessel_periodic_state`, `windkessel_summary`, `windkessel_terms`, `ejection_half_sine`, `time_mean` | C07 | C07 (D10, N46) |
| `scipy.integrate.solve_ivp` | C07 | knowledge/primers.md: scipy.integrate.solve_ivp (ch01 P31) — reminder |
| an RK4 loop | C07 | C07 (gloss in the from-scratch row; 🔁 Ch. 10) |
| pressure throughout the circulation (Fig. 16.5 idea) | C07 | C07 (N48) |
| instantaneous transmission (the model's limit) | C07 | C07 (N47) |
| ipywidgets `live` and `_ = live(...)` | C07 | knowledge/primers.md: slider_figure (ch01 P17) — reminder; C07 (gloss) |
| the 1-D tube equations; area as an unknown; `ch16.tube1d_run` | C08 | C08 (N49) |
| material derivative in one dimension | C08 | C08 (N50; 🔁 Ch. 3) |
| momentum per slice with the wall's push (S5) | C08 | C08 (N50, D11 step 2) |
| the tube law in A; `ch16.tube_law_hookean`, `tube_law_linear` | C08 | C08 (N51) |
| basic state + small disturbance; linearisation | C08 | C08 (D11; 🔁 Ch. 15 D01) |
| limits and orders of smallness | C08 | knowledge/primers.md: limits and orders of smallness (ch02 P68) — reminder |
| linearisation ratios; `ch16.linearisation_ratios` | C08 | C08 (N52) |
| Schwarz's theorem (mixed partials commute) | C08 | knowledge/primers.md: Schwarz's theorem (ch04 P121) — reminder |
| wave equation; Moens–Korteweg speed; Bramwell–Hill form; `ch16.moens_korteweg_speed`, `bramwell_hill_speed`, `youngs_modulus_from_pwv` | C08 | C08 (D12) |
| thick-wall speed (16.63) | C08 | C08 (N53) |
| d'Alembert's solution; time shapes; `ch16.pulse_wave`, `AC.dalembert` | C08 | R09 |
| the second-order wave equation and its two characteristic coordinates | C08 | knowledge/primers.md: the second-order wave equation and its two characteristic coordinates (ch15 P357) — reminder |
| p̂ = ρcû; `ch16.velocity_from_pulse`, `AC.acoustic_impedance` | C08 | R10 |
| wall bulge and pulse state; `ch16.wall_bulge`, `pulse_state` | C08 | C08 (gloss at A2 and E3) |
| pulse speeds, steepening, damping (named); MacCormack (DEVIATION box) | C08 | C08 (N54; 🔁 Ch. 10) |
| rigid-tube reduction; ∂p/∂x a function of t only; not a wave | C09 | R11 |
| Laplacian in cylindrical coordinates | C09 | knowledge/primers.md: Laplacian in cylindrical coordinates (ch08 P186) — reminder |
| two time conventions e^{+iωt} and e^{i(kx−ωt)} | C09 | C09 (trap T1 table) |
| complex amplitudes; Euler's formula | C09 | knowledge/primers.md: complex amplitudes (ch07 P176) — reminder |
| harmonic gradient and separable trial | C09 | C09 (N55) |
| Bessel's equation and J₀, Y₀, J₁; their series | C09 | primer (in C09) |
| powers of i: i^{1/2}, i^{3/2} | C09 | primer (in C09) |
| scipy.special.jv with a complex argument; `jve` | C09 | primer (in C09) |
| complex square roots and the quadratic formula | C09 | knowledge/primers.md: complex square roots and the quadratic formula (ch06 P159) — reminder |
| the Bessel ODE (corrected sign) and Womersley's profile; F₁; `ch16.womersley_F1`, `womersley_velocity`, `womersley_velocity_nd`, `womersley_ode_residual`, `bessel_j0_series`, `bessel_j1_series` | C09 | C09 (D13, N56) |
| small-α limit; `ch16.womersley_small_alpha` | C09 | C09 (D14) |
| binomial approximation with a complex small quantity | C09 | knowledge/primers.md: binomial approximation with a complex small quantity (ch14 P342) — reminder |
| large-α limit by a wall layer; dominant balance; `ch16.womersley_large_alpha` | C09 | C09 (D15, N58) |
| linear second-order ODE: trying e^{λy} | C09 | knowledge/primers.md: linear second-order ODE (ch01 P44) — reminder |
| Stokes' second problem; `LAM.stokes_second_problem` | C09 | C09 (🔁 Ch. 8 C11) |
| superposition of harmonics; α_n = α√n; `ch16.womersley_harmonics` | C09 | C09 (N57) |
| Fourier modes | C09 | knowledge/primers.md: Fourier modes and the FFT Poisson solver (ch05 P142) — reminder |
| oscillating profiles over a cycle; reversal at the wall | C09 | C09 (N59, N60) |
| Crank–Nicolson check; `ch16.womersley_pde_solve`; caching | C09 | knowledge/primers.md: caching expensive runs (ch10 P252) — reminder |
| two Bessel identities (dJ₀/dz = −J₁, ∫sJ₀ = zJ₁) | C10 | primer (in C10) |
| F₂ and the oscillating flux; `ch16.womersley_F2`, `womersley_flow_rate`, `womersley_flow_rate_harmonics` | C10 | C10 (D16, N61) |
| X, Y; amplitude and phase lag; `ch16.womersley_XY`, `womersley_flow_ratio`, `womersley_phase_lag` | C10 | C10 (D17) |
| the complex plane in numpy; `np.abs`, `np.angle` | C10 | knowledge/primers.md: the complex plane in numpy (ch06 P153) — reminder |
| `np.arctan2` | C10 | knowledge/primers.md: np.arctan2 (ch02 P70) — reminder |
| phase of a wave; lag in time | C10 | knowledge/primers.md: phase of a wave (ch07 P165) — reminder |
| wall shear stress of the oscillating flow; `ch16.womersley_wall_shear` | C10 | C10 (N62) |
| the two series of F₂; error slopes; `ch16.womersley_F2_series`, `womersley_F2_asymptote` | C10 | C10 (N63, N64) |
| linearised continuity and Navier–Stokes in (r, x) | C11 | R12 |
| membrane tensions and wall equations; Poisson coupling | C11 | C11 (N65, N66) |
| no slip at the undisturbed wall | C11 | C11 (N67; 🔁 Ch. 7 linearised boundary conditions) |
| a complex wavenumber; damping per wavelength | C11 | primer (in C11) |
| normal modes; long-wave ordering; the mode ODEs | C11 | C11 (D18, N68, N69) |
| reduced wall equations; `ch16._elastic_tube_modes` | C11 | C11 (D19, N70) |
| mode solution with β = √(iω/ν) | C11 | C11 (D20, N71) |
| non-trivial solutions of a homogeneous 2 × 2 system | C11 | knowledge/primers.md: non-trivial solutions of a homogeneous 2 × 2 system (ch13 P330) — reminder |
| choosing a root by continuation | C11 | primer (in C11) |
| the characteristic quadratic and its roots; `ch16.elastic_tube_characteristic`, `elastic_tube_roots`, `elastic_tube_wave`, `elastic_tube_determinant_numeric` | C11 | C11 (D21) |
| viscoelastic wall by complex moduli | C11 | C11 (N72) |
| waves in a stream U ± c; `AC.wave_speeds` | C11 | C11 (N73; 🔁 Ch. 15) |
| Murray's cost: pumping power and upkeep; `ch16.murray_cost`, `murray_terms` | C12 | C12 (N74) |
| the minimum of a falling plus a rising power | C12 | knowledge/primers.md: reading a drag polar; the minimum of a/x + bx (ch14 P349) — reminder |
| optimal radius, minimum cost; `ch16.murray_radius`, `murray_min_cost` | C12 | C12 (D22) |
| the cube law; `ch16.murray_daughters` | C12 | C12 (D23, N75) |
| scipy.optimize.minimize for a point in the plane | C12 | primer (in C12) |
| the weighted Fermat point; directional derivative of a length | C12 | primer (in C12) |
| level sets and the directional derivative | C12 | knowledge/primers.md: level sets and the directional derivative (ch02 P75) — reminder |
| the cosine rule | C12 | knowledge/primers.md: the cosine rule (ch14 P340) — reminder |
| branch angles; `ch16.murray_angles`, `murray_theta`, `murray_symmetric_angle`, `bifurcation_cost`, `bifurcation_point` | C12 | C12 (D24, N76) |
| junction conditions: pressure continuity, flux conservation | C13 | C13 (N77) |
| admittance and impedance; `ch16.tube_admittance`, `tube_impedance` | C13 | C13 (N78) |
| reflection and transmission coefficients; `ch16.junction_coefficients`, `junction_from_geometry` | C13 | C13 (D25) |
| `np.linalg.solve` | C13 | C13 (gloss in the from-scratch row) |
| wave energy flux and the energy identity; `ch16.junction_energy` | C13 | C13 (D26, N80) |
| parent waveforms after reflection; `ch16.junction_waveforms`, `network_pulse_run` | C13 | C13 (N79) |
| closed-end-like, matched, open-end-like | C13 | C13 (N81) |
| `np.random.default_rng` | C13 | C13 (gloss in the from-scratch row) |
| curvature effects on pipe flow (named) | C14 | C14 (N82) |
| toroidal coordinates of a curved pipe; centrifugal acceleration w²/R | C14 | primer (in C14) |
| Dean vortices; curvature ratio; `ch16.curvature_ratio` | C14 | C14 (N83, N84) |
| √δ rescaling and the pressure split (stated) | C14 | C14 (N85) |
| the Dean equations (stated) and the centrifugal term | C14 | C14 (N86) |
| stream function of the secondary flow | C14 | C14 (N88; 🔁 Ch. 4 stream function) |
| the vorticity equation (16.174) and its boundary conditions | C14 | C14 (N89) |
| Dean-number conventions κ, K, D, De; `ch16.dean_number`, `dean_convert` | C14 | C14 (N87) |
| a regular perturbation series | C14 | primer (in C14) |
| the biharmonic operator in polar coordinates; undetermined coefficients | C14 | primer (in C14) |
| integration by parts | C14 | C14 (gloss in D27 step 13; 🔁 Ch. 9) |
| ψ₀, w₁ (plus sign) and the flux coefficient; `ch16.dean_secondary_streamfunction`, `dean_axial_correction`, `dean_axial_velocity`, `dean_secondary_velocity`, `dean_series_sympy`, `dean_flux_ratio` | C14 | C14 (D27, N90, N91) |
| `np.meshgrid`; contour plots | C14 | C14 (gloss in the figure row; 🔁 Ch. 3) |
| flux ratio in D; resistance ratio; series reliability (ours); `ch16.curved_resistance_ratio`, `dean_series_reliability` | C14 | C14 (N92, N93) |
| collapse at negative transmural pressure | C15 | C15 (N94) |
| Korotkoff sounds (named) | C15 | C15 (N95) |
| the Starling resistor and its five cases; `ch16.starling_state`, `flow_limitation_curve` | C15 | C15 (N96) |
| a general tube law and its wave speed (Bramwell–Hill) | C15 | C15 (N97) |
| bending stiffness EI, I ∝ h³; `ch16.bending_stiffness_parameter` | C15 | primer (in C15) |
| Shapiro's P and the area ratio α; the collapsible tube law; `ch16.collapsible_tube_law` | C15 | C15 (N98, N99) |
| the collapsible wave speed; `ch16.collapsible_wave_speed` | C15 | C15 (D28, N100) |
| duct laws with friction and gravity; hydraulic diameter; `ch16.wall_friction_term` | C15 | R13 |
| speed index S; `ch16.speed_index` | C15 | C15 (N101) |
| logarithmic differentiation | C15 | knowledge/primers.md: logarithmic differentiation (ch01 P52) — reminder |
| differentials of algebraic laws as a small linear system | C15 | knowledge/primers.md: differentials of algebraic laws as a small linear system (ch15 P362) — reminder |
| differential relations (16.196)–(16.202); `ch16.collapsible_differentials_sympy` | C15 | C15 (D29, N102) |
| Shapiro's equations, M, the condition at S = 1; `ch16.shapiro_M`, `collapsible_rhs`, `collapsible_march`, `collapsible_terms`, `choking_flux` | C15 | C15 (D30) |
| the giraffe jugular (gravity against friction) | C15 | C15 (note in row 18) |
| Casson pipe flow motivation | C16 | C16 (N103) |
| Casson's law with its yield branch | C16 | C16 (N104) |
| linear stress in a pipe; plug radius; minimum pressure drop; `ch16.casson_wall_stress`, `casson_plug_radius`, `casson_min_pressure_drop` | C16 | C16 (D31, N105, N106) |
| algebra with square roots; ∫r^{1/2}dr; a polynomial in √ξ | C16 | primer (in C16) |
| the Casson profile and the plug speed; `ch16.casson_pipe_velocity`, `casson_core_velocity` | C16 | C16 (D32, N107) |
| the Casson flux; `ch16.casson_pipe_flow_rate`, `casson_flux_ratio` | C16 | C16 (D33, N108) |
| plants: xylem and phloem | C17 | C17 (N109) |
| osmosis and van 't Hoff's law | C17 | primer (in C17) |
| water potential; tension; `ch16.water_potential`, `water_potential_terms`, `osmotic_pressure` | C17 | C17 (D34) |
| hydrostatic law dp/dz = −ρg | C17 | C17 (D34; 🔁 Ch. 1 D05) |
| universal gas constant `TH.R_U` (per kmol) | C17 | C17 (gloss in the code row) |
| transpiration, cohesion, cavitation; capillary rise; `ch01.capillary_rise`, `ch16.cavitation_margin` | C17 | C17 (N110) |
| xylem Reynolds number and Poiseuille drop; `ch16.xylem_reynolds`, `xylem_column` | C17 | C17 (N111) |
| phloem pressure flow (Münch) | C17 | C17 (N112) |
| sieve-tube and sieve-plate pressure drop; `ch16.sieve_tube_pressure_drop` | C17 | C17 (N113) |
| exercises and literature (pointers) | C17 | C17 (A.5, S01, S02 pointers) |

---

## Part F — derivation storyboards

One block per `D` row of the curation (§4b), in order. Builders copy each block word for word into `nb.derivation(id, title,
ref=…, goal=…, start=(tex, plain), plan=[…], uses=[…], steps=[dict(did=…, tex=…, why=…, plain=…), …], result=(tex, plain),
interpret=…, check=…, check_src=…)` and, where the heading names an explainer, into that explainer's `derivations: [...]` (**the same
`did` titles, the same step count, the same order** — ch15 lesson; the explainer may shorten *why* to its first sentence on
phones). Every step has four parts: **did** (the move, ≤ 8 words) · **tex** (the new line) · **why** (why it is allowed and why we
make it, ≤ 35 words; a *why* never claims more than the line before gives — ch15 lesson) · **plain** (what the line says).
Wherever a step produces a numbered book equation, the number stands directly after that step's **tex**. Colours as header
convention 8. The eight ★★★ blocks carry a **sympy check** with a `check_src` listing that re-runs the derivation's own
construction; **each listing was executed as written for this design and passed** (the printed line is quoted under it).
Step counts that differ from the curation's estimate (one small move per step forced a split or a merge) are marked "(curation: n)".

### D01 · Work per beat as the loop area, $W=-\oint p\,dV$ (ours) — ★, 6 steps, in C01 (notebook · `windkessel_heartbeat`)
- **Goal:** find the mechanical work the left ventricle does on the blood in one beat, from its pressure–volume loop. ·
  **Start:** a pressure p acting on a wall that moves so that the chamber volume changes by dV — "the blood pushes on the wall,
  the wall pushes back". · **Plan:** • the work of one small volume change • add the four legs • the isovolumic legs give
  nothing • what remains is the area inside the loop. · **Tools:** line integral (knowledge/primers.md P35, reminder); loop area
  (P138, reminder); the shoelace formula (P378, primer in C01); units of physiology (P394, primer in C01). · **Assumptions:**
  uniform pressure in the chamber and slow (quasi-static) motion (step 1); the loop runs A→B→C→D→A (step 5).
- **Steps:**
  1. **did** Work of one small volume change · **tex** $dW_{\text{on blood}}=-p\,dV$ · **why** A pressure p on a wall that moves so
     the volume grows by dV does work p dV on the wall (force × distance over the whole wall). Work on the blood is the
     opposite. · **plain** When the chamber shrinks under pressure (dV < 0), the wall does positive work on the blood.
  2. **did** Add the four legs of the beat · **tex** $W=-\displaystyle\int_{A\to B}p\,dV-\int_{B\to C}p\,dV-\int_{C\to D}p\,dV-\int_{D\to A}p\,dV$ · **why**
     Work adds up along a path (a line integral, P35); one beat is the closed path A→B→C→D→A of the loop. · **plain** The work
     per beat is the sum of the work done on each leg.
  3. **did** Drop the two isovolumic legs · **tex** $\displaystyle\int_{B\to C}p\,dV=\int_{D\to A}p\,dV=0$ · **why** On B→C and D→A both
     valves of the ventricle are shut, the volume cannot change, so dV = 0 all along the leg whatever the pressure does. ·
     **plain** Squeezing a closed chamber raises its pressure but moves no blood: no work.
  4. **did** Write the other legs from ESV to EDV · **tex** $W=\displaystyle\int_{ESV}^{EDV}p_{\text{eject}}\,dV-\int_{ESV}^{EDV}p_{\text{fill}}\,dV$ ·
     **why** C→D runs from EDV down to ESV, so swapping its limits removes the minus sign; A→B already runs from ESV up to EDV. ·
     **plain** Work given during ejection (area under the top curve) minus work received during filling (area under the bottom
     curve).
  5. **did** Recognise the enclosed area and its sign · **tex** $W=-\oint p\,dV=\text{area inside the loop}>0$ · **why** The
     difference of the two areas under the curves is the area between them (P138). The loop runs counterclockwise in the
     (V, p) plane, so $\oint p\,dV<0$ (trap T11). · **plain** The heart's work per beat is the area of its loop.
  6. **did** Convert mL·mmHg to joules · **tex** $1\ \mathrm{mmHg\cdot mL}=133.32\ \mathrm{Pa}\times10^{-6}\ \mathrm{m^3}=1.333\times10^{-4}\ \mathrm J$ ·
     **why** An area on a (mL, mmHg) chart has units pressure × volume, an energy; the converters of P394 turn it into SI so it
     can be compared with a power. · **plain** A loop 100 mmHg tall and 75 mL wide holds about one joule.
- **Result:** $W=-\oint p\,dV$ (ours) = the area enclosed by the PV loop — "positive when the loop runs counterclockwise".
- **Check:** units Pa·m³ = N·m = J ✓; a rectangle 100 mmHg × 75 mL: W = 100 × 133.322 × 75 × 10⁻⁶ = 0.9999 J; zero stroke volume
  gives W = 0; running the loop backwards only flips the sign (`ch16.pv_loop_work(...)["orientation"]`). Our schematic loop (header
  convention 10): **expect** W = 1.0666 J exactly, 1.0666 J by the shoelace sum; cardiac power at 66 beats/min ≈ 1.17 W.
- **What it means:** the loop's width is the stroke volume, its height the pressure the heart must reach; their product sets
  the work. Fails when the chamber pressure is not uniform, and it ignores the kinetic energy of the ejected blood (a few per
  cent at rest).
- **Traps:** the isovolumic legs do no work although the pressure changes most there; counterclockwise in (V, p) means a
  negative $\oint p\,dV$ (T11); mL·mmHg is not a joule (P394).

### D02 · The Casson apparent viscosity, $\mu_{app}=\left(\sqrt{\tau_y/\dot\gamma}+K_c\right)^2$, its limits and $K_c=k_c\sqrt{\mu_p}$ (ours) — ★, 5 steps, in C02 (notebook · `casson_blood_flow`)
- **Goal:** find how thick blood *seems* at each shear rate according to Casson's law, and what happens when it is sheared very
  fast or hardly at all. · **Start:** $\sqrt{\tau/\mu_p}=k_c\sqrt{\dot\gamma}+\sqrt{\tau_y/\mu_p}$ (16.6) and the definition
  $\mu_{app}=\tau/\dot\gamma$ — "a straight line in √τ against √γ̇; apparent viscosity is stress over rate". · **Plan:** • write the
  law with one dimensional constant • square • divide by γ̇ • take the two limits. · **Tools:** secant against tangent (P383,
  primer in C02); square roots (P393 is later; here only $(x+y)^2$). · **Assumptions:** a homogeneous fluid obeying (16.6) at every
  shear rate (start).
- **Steps:**
  1. **did** Multiply the law by $\sqrt{\mu_p}$ · **tex** $\sqrt\tau=\sqrt{\tau_y}+K_c\sqrt{\dot\gamma},\quad K_c=k_c\sqrt{\mu_p}$ · **why**
     Multiplying both sides by the positive number √μ_p keeps the equality; it turns (16.6) into the form of (16.207). k_c is a
     pure number, K_c has units (Pa s)^½ (trap T14). · **plain** The book's two ways of writing Casson's law are one law.
  2. **did** Square both sides · **tex** $\tau=\left(\sqrt{\tau_y}+K_c\sqrt{\dot\gamma}\right)^2$ · **why** Both sides are positive, so
     squaring keeps the equality; we need τ itself to form τ/γ̇. · **plain** The stress needed to shear blood at a rate γ̇.
  3. **did** Divide by the shear rate · **tex** $\mu_{app}=\dfrac{\tau}{\dot\gamma}=\left(\sqrt{\dfrac{\tau_y}{\dot\gamma}}+K_c\right)^2$ · **why**
     The definition of the apparent viscosity — the slope of the secant from the origin (P383). Dividing a square by γ̇ divides
     each term inside by √γ̇. · **plain** How viscous the blood looks at this shear rate.
  4. **did** Let the shear rate grow large · **tex** $\dot\gamma\to\infty:\ \ \mu_{app}\to K_c^2=k_c^2\mu_p$ · **why** √(τ_y/γ̇) → 0;
     what is left is K_c². With Whitmore's $k_c=1.53$ of (16.7), $K_c^2=2.34\,\mu_p$. · **plain** Sheared hard, blood behaves like a
     Newtonian liquid about 2.3 times as viscous as plasma.
  5. **did** Let the shear rate go to zero · **tex** $\dot\gamma\to0:\ \ \mu_{app}\approx\dfrac{\tau_y}{\dot\gamma}\to\infty$ · **why**
     Now √(τ_y/γ̇) dominates K_c, and its square is τ_y/γ̇: a finite stress divided by a vanishing rate. · **plain** Almost at rest,
     blood seems infinitely viscous: it can hold a small stress without flowing.
- **Result:** $\mu_{app}=\left(\sqrt{\tau_y/\dot\gamma}+K_c\right)^2$ (ours) — "falls from infinity at rest to $K_c^2$ at high shear".
- **Check:** units — √(Pa/s⁻¹) = (Pa s)^½, the same as K_c ✓; τ_y = 0 gives μ_app = K_c², a Newtonian fluid ✓; numbers with
  (16.7) and μ_p = 1.2 mPa s: γ̇ = 4 s⁻¹ → μ_p(2/2 + 1.53)² = 6.40 μ_p = **7.68 mPa s**; γ̇ = 400 s⁻¹ → (0.1 + 1.53)² μ_p = **3.19
  mPa s**; the tangent dτ/dγ̇ at 4 s⁻¹ is only 4.65 mPa s (T27).
- **What it means:** blood is shear-thinning with a yield stress: slow flows (veins at rest, small vessels) feel a much higher
  viscosity than fast ones. Fails below a few per cent hematocrit (no yield stress) and in tubes so narrow that blood is not
  homogeneous (C05).
- **Traps:** the apparent viscosity is the secant τ/γ̇, not the tangent dτ/dγ̇ (T27); 1.53 is k_c, dimensionless, not K_c (T14).

### D03 · $\alpha^2=Re\,St/4$ and $\alpha=\sqrt2\,a/\delta_e$ (ours) — ★, 4 steps, in C03 (notebook · `womersley_profiles`)
- **Goal:** show what the Womersley number compares — the oscillation's inertia with viscosity, or the vessel with its Stokes
  layer. · **Start:** $\alpha=a\sqrt{\omega/\nu}$ (16.12); the chapter's $Re=2au/\nu$ and $St=2a\omega/u$ (both on the diameter);
  the Stokes-layer thickness $\delta_e=\sqrt{2\nu/\omega}$ (Ch. 8, the e-folding depth of $u(y,t)=U\exp\{-y\sqrt{\omega/2\nu}\}\cos(\omega t-y\sqrt{\omega/2\nu})$ (8.38)). · **Plan:** • multiply Re by St • compare
  with α² • write ω/ν with δ_e • substitute. · **Tools:** algebra; R02's link to separation is not needed. · **Assumptions:** the
  book's definitions of Re and St (step 1, trap T6).
- **Steps:**
  1. **did** Multiply the two groups · **tex** $Re\cdot St=\dfrac{2au}{\nu}\cdot\dfrac{2a\omega}{u}=\dfrac{4a^2\omega}{\nu}$ · **why** Plain
     multiplication; the speed u cancels, which is why α can be quoted for a vessel without knowing how fast the blood moves. ·
     **plain** Unsteadiness times inertia over viscosity, with no speed left.
  2. **did** Compare with α squared · **tex** $\alpha^2=\dfrac{a^2\omega}{\nu}=\dfrac{Re\cdot St}{4}$ · **why** Square (16.12) and divide
     step 1 by 4. Other sources define Re and St differently and write $\alpha^2=2\pi\,Re\,St$ (trap T6). · **plain** α² is a
     quarter of the product of the book's Reynolds and Strouhal numbers.
  3. **did** Write ω/ν with the Stokes layer · **tex** $\dfrac{\omega}{\nu}=\dfrac{2}{\delta_e^2}$ · **why** Square
     $\delta_e=\sqrt{2\nu/\omega}$ and invert; δ_e is how far an oscillation diffuses into a fluid in one radian of its cycle. ·
     **plain** Frequency over viscosity is fixed by the Stokes-layer thickness.
  4. **did** Substitute into (16.12) · **tex** $\alpha=a\sqrt{\dfrac{2}{\delta_e^2}}=\sqrt2\,\dfrac{a}{\delta_e}$ · **why** Replace ω/ν by
     step 3; the square root of 2/δ_e² is √2/δ_e. · **plain** α compares the vessel radius with the thickness of the viscous
     layer at its wall.
- **Result:** $\alpha^2=Re\,St/4$ and $\alpha=\sqrt2\,a/\delta_e$ (ours).
- **Check:** both sides dimensionless ✓; easy numbers: a = 5 mm, ω = 4 rad/s, ν = 4 × 10⁻⁶ m²/s → α = 0.005 × √(10⁶) = **5**;
  with u = 0.2 m/s, Re = 500, St = 0.2, Re·St/4 = 25 = α² ✓; δ_e = √(2 × 4 × 10⁻⁶/4) = 1.414 mm and √2 × 5/1.414 = 5 ✓.
- **What it means:** α ≪ 1: the viscous layer is wider than the vessel — viscosity controls every point (quasi-steady); α ≫ 1:
  the layer is a thin skin and the core moves by inertia. Fails as a guide when the flow is not periodic.
- **Traps:** Re and St on the diameter (T6); α is built from the radius.

### D04 · Vascular resistance $R=\dfrac{8\mu L}{\pi a^4}$ (16.20) from Poiseuille's flux — ★, 5 steps, in C04 (notebook)
- **Goal:** turn Poiseuille's profile into a single "resistance" number, and see why it depends so strongly on the radius. ·
  **Start:** $u(r)=\left(\dfrac{a^2\Delta p}{4\mu L}\right)\left(1-\dfrac{r^2}{a^2}\right)$ (16.15) — "the parabola of steady pipe flow (Ch. 8)". · **Plan:** • add the
  flow of thin rings • integrate • divide the pressure drop by the flow • halve the radius. · **Tools:** definite integral (P27,
  reminder); a ring of area 2πr dr. · **Assumptions:** steady, fully developed, Newtonian, rigid tube (R03).
- **Steps:**
  1. **did** Add the flow through thin rings · **tex** $Q=\displaystyle\int_0^au(r)\,2\pi r\,dr$ · **why** A ring of radius r and width dr has
     area 2πr dr and passes u(r)·2πr dr each second; the integral adds all rings (P27). · **plain** Volume per second is the sum
     over rings.
  2. **did** Substitute the profile and integrate · **tex** $Q=\dfrac{2\pi a^2\Delta p}{4\mu L}\left(\dfrac{a^2}{2}-\dfrac{a^2}{4}\right)$ ·
     **why** Constants come out of the integral; $\int_0^ar\,dr=a^2/2$ and $\int_0^ar^3\,dr/a^2=a^2/4$. · **plain** The two pieces of
     the parabola, added up.
  3. **did** Simplify · **tex** $Q=\dfrac{\pi a^4}{8\mu}\dfrac{\Delta p}{L}$ (16.17) · **why** $(2\pi a^2\Delta p/4\mu L)(a^2/4)=\pi a^4\Delta p/8\mu L$. ·
     **plain** Poiseuille's law: flow grows as the fourth power of the radius.
  4. **did** Divide the pressure drop by the flow · **tex** $R=\dfrac{\Delta\bar p}{Q}=\dfrac{8\mu L}{\pi a^4}$ (16.20) · **why** The definition
     $R=\Delta\bar p/Q$ (16.9), applied to one tube: invert (16.17). · **plain** The resistance of one vessel.
  5. **did** Halve the radius · **tex** $\dfrac{R(a/2)}{R(a)}=2^4=16$ · **why** R ∝ a⁻⁴, so a factor ½ in a gives a factor 2⁴ in R. ·
     **plain** Narrow a vessel to half its radius and it resists 16 times more — the arterioles' lever on blood pressure.
- **Result:** $R=\dfrac{8\mu L}{\pi a^4}$ (16.20).
- **Check:** units Pa s × m / m⁴ = Pa s m⁻³ ✓ (pressure per flow); easy numbers μ = 4 mPa s, L = 0.1 m, a = 1 mm: R = 8 × 0.004 ×
  0.1/(π × 10⁻¹²) = **1.02 × 10⁹ Pa s m⁻³**; the from-scratch trapezoid of 2πru gives the same R to 10⁻⁶.
- **What it means:** the body controls flow mainly through radius, because of the fourth power. Fails in tubes too short to be
  fully developed (R08), for pulsatile flow at α > 1 (C09) and for non-Newtonian blood in small vessels (C05, C16).
- **Traps:** the text says "diameter" but a is the radius; the book's last member of (16.19) has the wrong sign (S1).

### D05 · The two-layer profiles (16.27)–(16.28) — ★★, 9 steps, in C05 (notebook)
- **Goal:** find the velocity in a tube whose core (cells, viscosity μ_c) is wrapped in a thin plasma layer (viscosity μ_p) of
  thickness δ. · **Start:** the six conditions (16.21)–(16.26) with the signs corrected: $-\mu_c\dfrac{du^c}{dr}=\dfrac r2\dfrac{\Delta p}L$, $\dfrac{du^c}{dr}=0$ at
  $r=0$, $\tau_{xr}\rvert_c=\tau_{xr}\rvert_p$ and $u^c=u^p$ at $r=a-\delta$, $-\mu_p\dfrac{du^p}{dr}=\dfrac r2\dfrac{\Delta p}L$, $u^p=0$ at $r=a$. · **Plan:** • why the stress
  is $r\Delta p/2L$ in both layers • integrate in the plasma, apply no slip • integrate in the core • match at the interface. ·
  **Tools:** free body of a fluid cylinder (Ch. 8 D07 recalled); integrating a linear law twice. · **Assumptions:** steady,
  fully developed, two Newtonian layers, a sharp interface (start).
- **Steps:**
  1. **did** Balance forces on a fluid cylinder · **tex** $\tau_{xr}(r)\,2\pi rL=\pi r^2\Delta p$ · **why** Steady, fully developed:
     nothing accelerates, so the pressure push on the end faces equals the shear drag on the side (the free body of Ch. 8 D07).
     It holds across the interface because the stress is continuous there (16.23). · **plain** The stress at r carries the push on
     everything inside r.
  2. **did** Solve for the stress · **tex** $\tau_{xr}=\dfrac r2\dfrac{\Delta p}{L}$ · **why** Divide by 2πrL. This is the right side
     of (16.21) and (16.24) with the sign that (16.27)–(16.28) use; the book prints a minus there (slip S2). · **plain** The stress
     grows linearly from zero on the axis.
  3. **did** Write the plasma layer's law · **tex** $-\mu_p\dfrac{du^p}{dr}=\dfrac r2\dfrac{\Delta p}{L}$ (16.24, sign corrected) · **why** A
     Newtonian layer: τ = −μ du/dr, a positive magnitude here because u falls outward (trap T5). · **plain** In the plasma the
     speed drops outward at a rate set by μ_p.
  4. **did** Integrate once in the plasma · **tex** $u^p=-\dfrac{r^2}{4\mu_p}\dfrac{\Delta p}{L}+C_p$ · **why** $\int r\,dr=r^2/2$; divide by
     −μ_p; one constant appears. · **plain** A parabola of unknown height.
  5. **did** Apply no slip at the wall · **tex** $u^p=\dfrac{a^2}{4\mu_p}\dfrac{\Delta p}{L}\left[1-\left(\dfrac ra\right)^2\right]$ (16.27) · **why** (16.26):
     u^p(a) = 0 gives $C_p=a^2\Delta p/4\mu_pL$. · **plain** The plasma layer moves exactly as Poiseuille flow of plasma would.
  6. **did** Integrate once in the core · **tex** $u^c=-\dfrac{r^2}{4\mu_c}\dfrac{\Delta p}{L}+C_c$ · **why** The core law (16.21, sign
     corrected) integrated the same way. Its slope $-r\Delta p/2\mu_cL$ is zero at r = 0, so the regularity condition (16.22)
     holds without a logarithm. · **plain** The core is a flatter parabola, because μ_c > μ_p.
  7. **did** Match the speeds at the interface · **tex** $-\dfrac{(a-\delta)^2}{4\mu_c}\dfrac{\Delta p}{L}+C_c=\dfrac{a^2}{4\mu_p}\dfrac{\Delta p}{L}\left[1-\left(\dfrac{a-\delta}a\right)^2\right]$ ·
     **why** (16.25): the two fluids do not slip past each other, so the speeds agree at r = a − δ. · **plain** One equation for the
     core's constant.
  8. **did** Solve for the core's constant · **tex** $C_c=\dfrac{a^2}{4\mu_p}\dfrac{\Delta p}{L}\left[1-\left(\dfrac{a-\delta}a\right)^2+\dfrac{\mu_p}{\mu_c}\left(\dfrac{a-\delta}a\right)^2\right]$ ·
     **why** Move the core term across; write $(a-\delta)^2/4\mu_c=(a^2/4\mu_p)(\mu_p/\mu_c)((a-\delta)/a)^2$ so that one factor
     comes out. · **plain** The core rides on the plasma layer's speed at the interface.
  9. **did** Put the constant back · **tex** $u^c=\dfrac{a^2}{4\mu_p}\dfrac{\Delta p}{L}\left[1-\left(\dfrac{a-\delta}{a}\right)^2-\dfrac{\mu_p}{\mu_c}\left(\dfrac ra\right)^2+\dfrac{\mu_p}{\mu_c}\left(\dfrac{a-\delta}{a}\right)^2\right]$ (16.28) ·
     **why** Substitute step 8 into step 6 and write $r^2/4\mu_c=(a^2/4\mu_p)(\mu_p/\mu_c)(r/a)^2$. · **plain** The core profile: a
     plug speed plus a gentle parabola.
- **Result:** (16.27) in the plasma, (16.28) in the core — "speeds match at the interface, slopes do not".
- **Check:** μ_c = μ_p turns (16.28) into (16.27): one Poiseuille parabola ✓; δ → 0 gives Poiseuille with μ_c ✓; the slopes at
  the interface differ by the factor μ_c/μ_p (the kink in C05's figure); `ch16.two_layer_velocity` satisfies (16.21)–(16.26) in
  the sympy test.
- **What it means:** a thin, runny layer at the wall lets the whole core slide — the lubrication behind the Fåhræus–Lindqvist
  effect. Fails if the cells reach the wall (capillaries) or the interface is diffuse.
- **Traps:** the book's signs in (16.21), (16.24) (S2); continuous stress but not continuous shear rate.

### D06 · The layer fluxes (16.29), (16.30, corrected), (16.31) — ★★, 8 steps, in C05 (notebook)
- **Goal:** add up the flow in the plasma layer and in the core, and find the total. · **Start:** (16.27) and (16.28); shorthand
  $b\equiv a-\delta$. · **Plan:** • plasma ring integral • the core as plug + parabola • add, then factor. · **Tools:** definite
  integral (P27); substitution in an integral (P106, reminder); sympy (P40, reminder) in the cell that shows slip S3. ·
  **Assumptions:** as D05.
- **Steps:**
  1. **did** Write the plasma flux as rings · **tex** $Q_p=2\pi\displaystyle\int_b^au^p\,r\,dr,\quad b\equiv a-\delta$ · **why** D04 step 1,
     only over the plasma annulus from b to a. · **plain** Flow carried by the plasma layer.
  2. **did** Integrate the plasma parabola · **tex** $Q_p=\dfrac{\pi\Delta p}{8\mu_pL}\left[a^2-(a-\delta)^2\right]^2$ (16.29) · **why**
     With s = a² − r², ds = −2r dr: $\int_b^a(a^2-r^2)r\,dr=(a^2-b^2)^2/4$ (P106); times $2\pi\Delta p/4\mu_pL$. · **plain** The
     layer's flux grows as the square of its area fraction.
  3. **did** Rewrite the core profile as plug plus parabola · **tex** $u^c=\dfrac{\Delta p}{4L}\left[\dfrac{a^2-b^2}{\mu_p}+\dfrac{b^2-r^2}{\mu_c}\right]$ ·
     **why** Multiply out (16.28): $a^2[1-b^2/a^2]/\mu_p=(a^2-b^2)/\mu_p$ and the two μ_p/μ_c terms give $(b^2-r^2)/\mu_c$. ·
     **plain** The core slides at the plasma's interface speed and adds its own parabola.
  4. **did** Integrate the plug part · **tex** $2\pi\displaystyle\int_0^b\dfrac{\Delta p}{4L}\dfrac{a^2-b^2}{\mu_p}r\,dr=\dfrac{\pi\Delta p}{4\mu_pL}(a^2-b^2)\,b^2$ ·
     **why** A constant times $\int_0^br\,dr=b^2/2$. · **plain** The plug carries its speed times the core area.
  5. **did** Integrate the core parabola · **tex** $2\pi\displaystyle\int_0^b\dfrac{\Delta p}{4L}\dfrac{b^2-r^2}{\mu_c}r\,dr=\dfrac{\pi\Delta p}{8\mu_cL}b^4$ · **why**
     $\int_0^b(b^2-r^2)r\,dr=b^4/4$, the substitution of step 2. · **plain** Poiseuille flow of the core fluid inside radius b.
  6. **did** Add and factor $b^2$ · **tex** $Q_c=\dfrac{\pi(a-\delta)^2\Delta p}{4\mu_pL}\left[a^2-\left(1-\dfrac{\mu_p}{2\mu_c}\right)(a-\delta)^2\right]$ (16.30, corrected) ·
     **why** $\dfrac{\pi\Delta p\,b^2}{4\mu_pL}\left[(a^2-b^2)+\dfrac{\mu_p}{2\mu_c}b^2\right]$, then collect b². The book prints a factor a²
     and $(a-\delta)^4/a^2$ (slip S3); only this form adds up to (16.31). · **plain** The core's flux.
  7. **did** Add the two fluxes over a common factor · **tex** $Q=\dfrac{\pi\Delta p}{8\mu_pL}\left[a^4-b^4+\dfrac{\mu_p}{\mu_c}b^4\right]$ · **why**
     $(a^2-b^2)^2+2b^2(a^2-b^2)=(a^2-b^2)(a^2+b^2)=a^4-b^4$; the parabola term is $(\pi\Delta p/8\mu_pL)(\mu_p/\mu_c)b^4$. ·
     **plain** Total flow: a plasma-viscosity tube, minus what the slower core loses.
  8. **did** Factor $a^4$ · **tex** $Q=\dfrac{\pi a^4\Delta p}{8\mu_pL}\left[1-\left(1-\dfrac\delta a\right)^4\left(1-\dfrac{\mu_p}{\mu_c}\right)\right]$ (16.31) · **why**
     $b^4/a^4=(1-\delta/a)^4$; collect the two b⁴ terms. · **plain** Poiseuille's flux for plasma, reduced by a factor that the
     core's thickness and viscosity set.
- **Result:** (16.29), (16.30, corrected), (16.31).
- **Check:** δ → 0: $Q=\pi a^4\Delta p/8\mu_cL$ ✓; μ_c = μ_p: $Q=\pi a^4\Delta p/8\mu_pL$ ✓; the notebook cell: quad of 2πru equals
  (16.31) to 10⁻¹⁰, and printed (16.30) minus corrected (16.30) $=\pi\Delta p\,a^2\delta(\delta-2a)/(4\mu_pL)\neq0$ (analysis S3).
- **What it means:** the core carries almost all the flow but the plasma layer sets how easily it moves. Fails as D05.
- **Traps:** slip S3; forgetting that the core also moves with the interface speed (the plug part).

### D07 · The apparent viscosity (16.33) and its thin-layer limit (16.34) — ★, 6 steps (curation: 5), in C05 (notebook)
- **Goal:** express the two-layer flow as one "apparent" viscosity and see how a thin layer lowers it. · **Start:** (16.31) and
  the definition (16.32). · **Plan:** • define μ_app by Poiseuille's formula • equate and invert • expand for thin layers. ·
  **Tools:** binomial approximation (knowledge/primers.md P26 and P342, reminder). · **Assumptions:** δ ≪ a in steps 4–6.
- **Steps:**
  1. **did** Define the apparent viscosity · **tex** $Q=\dfrac{\pi a^4\Delta p}{8\mu_{app}L}$ (16.32) · **why** By definition μ_app is
     the viscosity a uniform Newtonian fluid would need to pass the same Q through the same tube under the same Δp — what an
     experiment measures. · **plain** Pretend the blood is one fluid.
  2. **did** Equate with (16.31) · **tex** $\dfrac1{\mu_{app}}=\dfrac1{\mu_p}\left[1-\left(1-\dfrac\delta a\right)^4\left(1-\dfrac{\mu_p}{\mu_c}\right)\right]$ ·
     **why** Same Q, Δp, L and a: the common factor $\pi a^4\Delta p/8L$ cancels. · **plain** The reciprocal viscosity is the
     bracket over μ_p.
  3. **did** Take reciprocals · **tex** $\mu_{app}=\mu_p\left[1-\left(1-\dfrac\delta a\right)^4\left(1-\dfrac{\mu_p}{\mu_c}\right)\right]^{-1}$ (16.33) ·
     **why** Both sides positive. · **plain** The two-layer apparent viscosity.
  4. **did** Expand the fourth power · **tex** $\left(1-\dfrac\delta a\right)^4\approx1-4\dfrac\delta a$ · **why** Binomial approximation
     with δ/a ≪ 1 (P26); the dropped terms are of order (δ/a)². · **plain** A thin layer removes a fraction 4δ/a of the core's
     fourth-power weight.
  5. **did** Multiply out the bracket · **tex** $\mu_{app}\approx\mu_p\left[\dfrac{\mu_p}{\mu_c}+4\dfrac\delta a\left(1-\dfrac{\mu_p}{\mu_c}\right)\right]^{-1}$ · **why**
     $1-(1-4\delta/a)(1-\mu_p/\mu_c)=\mu_p/\mu_c+4(\delta/a)(1-\mu_p/\mu_c)$, exactly. · **plain** The core's share plus the layer's
     correction.
  6. **did** Factor out $\mu_p/\mu_c$ · **tex** $\mu_{app}\approx\mu_c\left[1+4\dfrac\delta a\left(\dfrac{\mu_c}{\mu_p}-1\right)\right]^{-1}$ (16.34) · **why**
     $(\mu_c/\mu_p)(1-\mu_p/\mu_c)=\mu_c/\mu_p-1$; $\mu_p\cdot\mu_c/\mu_p=\mu_c$. · **plain** The core's viscosity, lowered by a
     factor that grows with δ/a.
- **Result:** (16.33) and its thin-layer form (16.34).
- **Check:** δ = 0 → μ_c; μ_c = μ_p → μ_p ✓; easy numbers μ_p = 1, μ_c = 4, δ/a = 0.1: (16.33) gives **1.969**, (16.34) **1.818**
  (first order is 8 % low at δ/a = 0.1); our 40 µm tube: **2.22 mPa s** against 3.73 mPa s at 500 µm (expect).
- **What it means:** μ_app < μ_c whenever μ_c > μ_p, and the drop grows as the tube narrows at fixed δ — the Fåhræus–Lindqvist
  effect. Fails when δ is comparable with a (capillaries) and the curve depends on our choice of δ (T19).
- **Traps:** (16.34) is first order in δ/a; the effect needs the *less* viscous fluid at the wall.

### D08 · The hoop stress (16.39) and the tube law (16.41) — ★, 6 steps, in C06 (notebook · `pulse_wave_artery`)
- **Goal:** find how far a thin elastic artery widens under a pressure difference. · **Start:** $2\sigma_{\theta\theta}h=\displaystyle\int_0^\pi(p-p_e)\,a\sin\theta\,d\theta$ (16.38) and
  $e_{\theta\theta}=\dfrac{\sigma_{\theta\theta}}{E}=\dfrac{a-a_0}{a_0}$ (16.40). · **Plan:** • cut the tube in half and balance forces • integrate • apply Hooke's
  law • solve for a. · **Tools:** Young's modulus and Hooke's law (P379, primer in C06); a curved membrane under tension (P380,
  primer in C06). · **Assumptions:** thin wall h ≪ a; Hookean; radial stress in the wall neglected; the deformed radius a in
  the stress (step 4).
- **Steps:**
  1. **did** Cut the tube along a diameter · **tex** $2\sigma_{\theta\theta}h=\displaystyle\int_0^\pi(p-p_e)\,a\sin\theta\,d\theta$ (16.38) · **why** Free body of
     the upper half, per unit length (P380): the two cut walls pull down with σh each; the net pressure on a strip a dθ has upward
     component (p − p_e) a sin θ dθ. · **plain** The wall's tension holds the half-tube against the pressure.
  2. **did** Do the integral · **tex** $2\sigma_{\theta\theta}h=2(p-p_e)a$ · **why** (p − p_e)a is the same all around;
     $\int_0^\pi\sin\theta\,d\theta=2$. · **plain** The pressure pushes on the projected width 2a.
  3. **did** Solve for the hoop stress · **tex** $\sigma_{\theta\theta}=\dfrac{(p-p_e)a}{h}$ (16.39) · **why** Divide by 2h. · **plain**
     Laplace's law: a wide, thin, pressurised tube carries a large wall stress.
  4. **did** Apply Hooke's law to the hoop · **tex** $\dfrac{(p-p_e)a}{Eh}=\dfrac{a-a_0}{a_0}$ · **why** (16.40): strain = stress/E
     (P379), and the hoop strain is the relative growth of the circumference 2πa. The book keeps the deformed a in the stress. ·
     **plain** The circumference stretches in proportion to the hoop stress.
  5. **did** Collect the terms in a · **tex** $a\left[\dfrac1{a_0}-\dfrac{p-p_e}{Eh}\right]=1$ · **why** Move $(p-p_e)a/Eh$ to the right
     and add 1 to both sides; the equation is linear in a. · **plain** One unknown, one linear equation.
  6. **did** Solve for a · **tex** $a=a_0\left[1-\dfrac{a_0}{Eh}(p-p_e)\right]^{-1}$ (16.41) · **why** Divide by the bracket and multiply top and
     bottom by a₀. Allowed only while $a_0(p-p_e)/Eh<1$; at 1 the radius would be infinite (trap T16). · **plain** The tube law.
- **Result:** $\sigma_{\theta\theta}=\dfrac{(p-p_e)a}{h}$ (16.39) and $a=a_0\left[1-\dfrac{a_0}{Eh}(p-p_e)\right]^{-1}$ (16.41).
- **Check:** units Pa·m/m = Pa ✓ and $a_0p/Eh$ dimensionless ✓; small pressure: $a\approx a_0[1+a_0(p-p_e)/Eh]$, linear; easy numbers
  Eh/a₀ = 50 kPa, p − p_e = 10 kPa: a/a₀ = 1/(1 − 0.2) = **1.25** (linear 1.20); our aorta at 100 mmHg: **a = 16.36 mm**, σ = 145
  kPa; singular at Eh/a₀ = 375 mmHg (expect).
- **What it means:** an artery is a spring whose stiffness is Eh/a₀: everything later (compliance, pulse speed) follows. Fails
  for real arteries, which stiffen at large strain (collagen) rather than blow up.
- **Traps:** the singularity is an artefact of using the deformed radius with a linear law (T16); transmural, not absolute,
  pressure.

### D09 · The two-element Windkessel (16.49) from a volume balance — ★, 4 steps, in C07 (notebook · `windkessel_heartbeat`)
- **Goal:** turn "the aorta stores blood in systole and releases it in diastole" into an equation for the aortic pressure. ·
  **Start:** (16.47) in words — inflow into the aorta = rate of change of the chamber's volume + outflow into the periphery. ·
  **Plan:** • write the balance • trade dv/dt for dp/dt • let the periphery be a resistor • substitute. · **Tools:** control
  volume (Ch. 4 recap); chain rule (P49, reminder); the electrical analogy (P381, primer in C07). · **Assumptions:** incompressible
  blood (step 1); the chamber volume depends on its pressure only (step 2); a steady peripheral resistor, venous pressure zero,
  transmural pressures (step 3, trap T12); no wave travel (N47).
- **Steps:**
  1. **did** Write the volume balance · **tex** $Q(t)=\dfrac{dv}{dt}+Q_{out}(t)$ (16.47) · **why** Blood is incompressible: what the heart
     pushes into the aortic chamber either stays there (v grows) or leaves through the periphery — a control-volume balance
     (Ch. 4). · **plain** Inflow = storage rate + outflow.
  2. **did** Trade the volume rate for a pressure rate · **tex** $\dfrac{dv}{dt}=\left(\dfrac{dv}{dp}\right)\left(\dfrac{dp}{dt}\right)=\mathrm K\dfrac{dp}{dt}$ (16.48) · **why**
     v depends on t only through p(t), so the chain rule applies (P49); K = dv/dp is the compliance, a capacitor (P381). ·
     **plain** The aorta stores more blood when its pressure rises faster.
  3. **did** Let the periphery be a resistor · **tex** $Q_{out}=\dfrac{p(t)-\bar p_V}{R}=\dfrac{p(t)}{R}$ · **why** $R=\Delta\bar p/Q$ (16.9) used at
     every instant (the small vessels store nothing); pressures are transmural and the venous pressure is taken as zero (T12). ·
     **plain** The higher the aortic pressure, the faster blood drains away.
  4. **did** Substitute both into the balance · **tex** $Q(t)=\mathrm K\dfrac{dp}{dt}+\dfrac{p(t)}{R}$ (16.49) · **why** Steps 2 and 3 in step 1. ·
     **plain** The heart's inflow charges the capacitor and drains through the resistor.
- **Result:** $Q(t)=\mathrm K\dfrac{dp}{dt}+\dfrac{p(t)}{R}$ (16.49).
- **Check:** units — K dp/dt: (m³/Pa)(Pa/s) = m³/s; p/R: Pa/(Pa s m⁻³) = m³/s ✓; K → 0 (a rigid aorta): p = QR, the pressure
  follows the inflow and falls to zero in diastole; constant Q: p → QR.
- **What it means:** one storage element and one resistor already explain why the aortic pressure never drops to zero between
  beats. Fails when the pulse's travel time matters (C08) or reflections shape the pressure (C13).
- **Traps:** p is transmural with p̄_V = 0 (T12); K here is the compliance, not Murray's cost or Dean's K (T2).

### D10 · The Windkessel pressure (16.52), the diastolic decay and the periodic state (ours) — ★★, 9 steps (curation: 8), in C07 (notebook · `windkessel_heartbeat`)
- **Goal:** solve (16.49) for the pressure, read off the decay between beats, and find the pressure the heartbeat settles
  into. · **Start:** $Q(t)=\mathrm K\dfrac{dp}{dt}+\dfrac{p(t)}{R}$ (16.49) with $p(0)=p_0$. · **Plan:** • standard form • integrating factor • integrate •
  switch the inflow off • ask for the same pressure after one beat • average. · **Tools:** first-order linear ODE and the
  integrating factor (P210, reminder); product rule (P38, reminder); a sum of decaying kicks (P382, primer in C07). ·
  **Assumptions:** R and K constant (step 1); a periodic heartbeat (step 8).
- **Steps:**
  1. **did** Divide by K · **tex** $\dfrac{dp}{dt}+\dfrac{1}{R\mathrm K}p=\dfrac{Q(t)}{\mathrm K}$ · **why** K > 0; this is the form $Q=\dfrac{dy}{dx}+Py$ (16.50)
     with y = p, x = t, P = 1/RK — and the book's "Q" of (16.50) is our forcing Q/K, not the flow. · **plain** The pressure relaxes
     at the rate 1/RK and is pushed by Q/K.
  2. **did** Form the integrating factor · **tex** $e^{\int P\,dt}=e^{t/R\mathrm K}$ · **why** P = 1/RK is constant, so its integral is t/RK
     (P210). · **plain** A factor that grows exactly as fast as the pressure decays.
  3. **did** Multiply through by the factor · **tex** $\dfrac{d}{dt}\left(p\,e^{t/R\mathrm K}\right)=\dfrac{Q(t)}{\mathrm K}e^{t/R\mathrm K}$ · **why** Product rule
     (P38): $\dfrac{d}{dt}(pe^{t/R\mathrm K})=e^{t/R\mathrm K}\left(\dfrac{dp}{dt}+\dfrac{p}{R\mathrm K}\right)$ — the left side of step 1 times the factor. This is
     (16.51) in disguise. · **plain** The left side became one exact derivative.
  4. **did** Integrate from 0 to t · **tex** $p(t)e^{t/R\mathrm K}-p_0=\dfrac1{\mathrm K}\displaystyle\int_0^tQ(\tau)e^{\tau/R\mathrm K}d\tau$ · **why** Fundamental theorem
     of calculus with dummy variable τ; p(0) = p₀ plays the constant A of (16.51). · **plain** Everything the heart put in, each
     piece weighted by when it came.
  5. **did** Multiply by $e^{-t/R\mathrm K}$ · **tex** $p(t)=\dfrac1{\mathrm K}e^{-t/R\mathrm K}\displaystyle\int_0^tQ(\tau)e^{\tau/R\mathrm K}d\tau+p_0e^{-t/R\mathrm K}$ (16.52) · **why**
     Isolate p. Read as a convolution (P382): each kick Q(τ)dτ raises p by Q dτ/K and then fades as $e^{-(t-\tau)/R\mathrm K}$. ·
     **plain** Today's pressure is the faded memory of every past stroke plus the faded initial pressure.
  6. **did** Switch the inflow off for diastole · **tex** $p(t)=p(T_s)\,e^{-(t-T_s)/R\mathrm K},\quad T_s\le t<T$ · **why** The aortic valve is
     shut, Q = 0; restart (16.52) at t = T_s: only its second term survives. · **plain** Between beats the pressure decays
     exponentially with time constant RK.
  7. **did** Check with a constant inflow · **tex** $Q=\bar Q:\ \ p(t)=\bar QR\left(1-e^{-t/R\mathrm K}\right)+p_0e^{-t/R\mathrm K}$ · **why**
     $\int_0^te^{\tau/R\mathrm K}d\tau=R\mathrm K(e^{t/R\mathrm K}-1)$; multiply by $(\bar Q/\mathrm K)e^{-t/R\mathrm K}$. · **plain** A steady inflow charges the aorta
     up to $\bar QR$ in a few RK.
  8. **did** Demand the same pressure after one beat · **tex** $p_0=\dfrac{\dfrac1{\mathrm K}e^{-T/R\mathrm K}\displaystyle\int_0^TQ(\tau)e^{\tau/R\mathrm K}d\tau}{1-e^{-T/R\mathrm K}}$ (ours) · **why**
     A steady heartbeat repeats, p(T) = p₀; (16.52) at t = T is linear in p₀; move $p_0e^{-T/R\mathrm K}$ left and divide by
     $1-e^{-T/R\mathrm K}\neq0$. · **plain** The end-diastolic pressure the circulation settles into.
  9. **did** Average (16.49) over one beat · **tex** $\bar p=R\bar Q$ · **why** Integrate (16.49) from 0 to T: $\mathrm K[p(T)-p(0)]=0$ by step 8,
     leaving $\frac1T\int p/R\,dt=\frac1T\int Q\,dt$. · **plain** The mean pressure is set by R alone; compliance only shapes the
     swing.
- **Result:** (16.52); diastole $p=p(T_s)e^{-(t-T_s)/R\mathrm K}$; periodic $p_0$ (ours); $\bar p=R\bar Q$.
- **Check:** units of $\frac1{\mathrm K}\int Q\,d\tau$: (Pa/m³)(m³) = Pa ✓; K → ∞: p stays at p₀; our resting set (R = 1.45 × 10⁸
  Pa s m⁻³, K = 1.2 × 10⁻⁸ m³/Pa, RK = 1.74 s, half-sine ejection, 66 beats/min) — **expect** periodic p_S = 106.0 mmHg, p_D = 74.5
  mmHg, mean 89.7 mmHg = R Q̄ ✓, decay factor over diastole $e^{-(T-T_s)/R\mathrm K}=0.712$; the closed form of C.12 equals the RK4
  march to 10⁻⁶.
- **What it means:** RK is the aorta's memory; with RK ≈ 2 beats the pressure has no time to fall far — a high diastolic
  pressure. A stiffer aorta (smaller K) shortens the memory and widens the swing at the same mean. Fails as D09.
- **Traps:** divide by K before using (16.50); p₀ is the pressure at the start of ejection (the end of diastole), not the
  minimum (which comes a moment later).

### D11 · The linearised tube equations (16.56)–(16.58) — ★★, 7 steps, in C08 (notebook · `pulse_wave_artery`)
- **Goal:** reduce the full one-dimensional laws of a blood-filled elastic tube to linear equations for a small pulse. ·
  **Start:** $\dfrac{\partial A}{\partial t}+\dfrac{\partial(Au)}{\partial x}=0$ (16.53), momentum (16.54, corrected) and $p-p_e=\dfrac{Eh}{a_0}\left[1-\left(\dfrac{A_0}{A}\right)^{1/2}\right]$ (16.55). · **Plan:** • split into rest
  + small parts • the corrected momentum law • drop products • Taylor-expand the tube law. · **Tools:** basic state + small
  disturbance (Ch. 15 D01, recap); orders of smallness (P68, reminder); Taylor (P26, reminder). · **Assumptions:** small area
  change, small speed u ≪ c, long waves (steps 3–5); p_e constant (step 2).
- **Steps:**
  1. **did** Split each unknown into rest plus a small part · **tex** $A=A_0+A',\quad u=u',\quad p-p_e=p'$ · **why** The undisturbed
     state is a tube of area A₀ with the blood at rest and zero transmural pressure; primes are small (basic state + disturbance,
     Ch. 15 D01). · **plain** A small pulse on a resting tube.
  2. **did** Write momentum with the wall's push · **tex** $\rho A\left(\dfrac{\partial u}{\partial t}+u\dfrac{\partial u}{\partial x}\right)=-A\dfrac{\partial p}{\partial x}$ (16.54, corrected) · **why**
     On a slice the end pressures give $-\partial(pA)/\partial x$ and the sloping wall pushes back with $+p\,\partial A/\partial x$:
     the sum is $-A\,\partial p/\partial x$. The book prints $-\partial((p-p_e)A)/\partial x$ (slip S5). · **plain** Blood accelerates down
     the pressure gradient.
  3. **did** Linearise the area law · **tex** $\dfrac{\partial A'}{\partial t}+A_0\dfrac{\partial u'}{\partial x}=0$ (16.56) · **why** $\partial(Au)/\partial x=A_0\,\partial u'/\partial x+\partial(A'u')/\partial x$;
     A'u' is a product of two small quantities, dropped (P68). · **plain** The tube swells where the flow converges.
  4. **did** Divide momentum by A and drop the product · **tex** $\rho\dfrac{\partial u'}{\partial t}=-\dfrac{\partial p'}{\partial x}$ (16.57) · **why** A > 0;
     $u'\partial u'/\partial x$ is second order (small because u ≪ c); p_e is constant, so $\partial p/\partial x=\partial p'/\partial x$. The extra term of the
     printed (16.54) would also drop here — slip S5 changes only the nonlinear model. · **plain** Newton's law for each slice.
  5. **did** Expand the square root of the tube law · **tex** $\left(\dfrac{A_0}{A}\right)^{1/2}=\left(1+\dfrac{A'}{A_0}\right)^{-1/2}\approx1-\dfrac{A'}{2A_0}$ · **why** Taylor
     (binomial) expansion for $A'/A_0\ll1$ (P26). · **plain** A small swelling changes the square root by half as much.
  6. **did** Substitute into (16.55) · **tex** $p-p_e=\dfrac{Eh}{2a_0A_0}(A-A_0)$ (16.58) · **why** $\dfrac{Eh}{a_0}\left[1-\left(1-\dfrac{A'}{2A_0}\right)\right]=\dfrac{Eh}{a_0}\dfrac{A'}{2A_0}$. · **plain** For small
     pulses the wall is a linear spring.
  7. **did** Read off the stiffness · **tex** $\dfrac{\partial p}{\partial A}=\dfrac{Eh}{2a_0A_0}$ (16.58) · **why** The derivative of the
     linear law. · **plain** How many pascals one more square metre of area costs.
- **Result:** (16.56), (16.57), (16.58).
- **Check:** units of Eh/(2a₀A₀): Pa·m/(m·m²) = Pa/m² ✓; our aorta with a 25 mmHg pulse: p′/(Eh/a₀) = **0.067**, u/c = **0.13**
  (expect; `ch16.linearisation_ratios`) — both small.
- **What it means:** small pulses obey linear laws, so they add and travel without changing shape (D12). Fails for large
  pulses (steepening, N54) or fast flow (u comparable with c, C15).
- **Traps:** slip S5; dropping u ∂u/∂x needs u ≪ c, not merely a small pulse.

### D12 · The wave equation (16.61) and the Moens–Korteweg speed (16.62) — ★★, 7 steps, in C08 (notebook · `pulse_wave_artery`)
- **Goal:** combine the linear laws into one wave equation for the pressure and read off its speed. · **Start:** (16.56),
  (16.57), (16.58). · **Plan:** • differentiate across • swap the mixed derivative • eliminate u • turn area into pressure. ·
  **Tools:** cross-differentiation (Ch. 15 D02, recap); Schwarz's theorem (P121, reminder); chain rule (P49, reminder). ·
  **Assumptions:** a uniform tube: A₀, E, h, a₀ constant (steps 1, 5).
- **Steps:**
  1. **did** Differentiate the area law in time · **tex** $\dfrac{\partial^2A}{\partial t^2}+A_0\dfrac{\partial^2u}{\partial t\,\partial x}=0$ · **why** ∂/∂t of (16.56);
     A₀ is constant. · **plain** How the swelling accelerates.
  2. **did** Differentiate momentum in space · **tex** $\rho\dfrac{\partial^2u}{\partial x\,\partial t}=-\dfrac{\partial^2p}{\partial x^2}$ · **why** ∂/∂x of (16.57). · **plain**
     How the acceleration changes along the tube.
  3. **did** Swap the order of the mixed derivative · **tex** $\dfrac{\partial^2u}{\partial t\,\partial x}=\dfrac{\partial^2u}{\partial x\,\partial t}$ · **why** Schwarz's
     theorem (P121) for a smooth u; it lets steps 1 and 2 share a term. · **plain** Same quantity, two orders.
  4. **did** Eliminate the velocity · **tex** $\dfrac{\partial^2A}{\partial t^2}=\dfrac{A_0}{\rho}\dfrac{\partial^2p}{\partial x^2}$ (16.59) · **why** Step 2 gives the mixed
     term as $-\frac1\rho\partial^2p/\partial x^2$; put it in step 1. · **plain** The tube swells faster where the pressure is
     curved.
  5. **did** Turn area into pressure · **tex** $\dfrac{\partial^2p}{\partial t^2}=\dfrac{\partial p}{\partial A}\dfrac{\partial^2A}{\partial t^2}$ (16.60) · **why** p depends on t only
     through A by (16.58), with a constant slope ∂p/∂A, so the chain rule (P49) applied twice gives no extra term. · **plain**
     Pressure follows area with a fixed factor.
  6. **did** Combine · **tex** $\dfrac{\partial^2p}{\partial t^2}=c^2\dfrac{\partial^2p}{\partial x^2},\quad c^2=\dfrac{A_0}{\rho}\dfrac{\partial p}{\partial A}$ (16.61) · **why** Put
     (16.59) into (16.60); the constant is positive (a tube that resists swelling), so it is a speed squared. · **plain** The
     pressure obeys the wave equation of Ch. 15.
  7. **did** Insert the Hookean slope · **tex** $c=\sqrt{\dfrac{A_0}{\rho}\dfrac{Eh}{2a_0A_0}}=\sqrt{\dfrac{Eh}{2\rho a_0}}$ (16.62) · **why** (16.58)'s slope; A₀
     cancels. With a general tube law the second form $c^2=\dfrac A\rho\dfrac{dp}{dA}$ stays (Bramwell–Hill). · **plain** The
     Moens–Korteweg speed: wall stiffness against blood inertia.
- **Result:** (16.61) with $c=\sqrt{\dfrac{Eh}{2\rho a_0}}=\sqrt{\dfrac A\rho\dfrac{dp}{dA}}$ (16.62).
- **Check:** units Pa·m/(kg m⁻³·m) = m²/s² ✓; easy numbers E = 0.72 MPa, h/a₀ = 0.1, ρ = 1000 kg/m³: c² = 36, **c = 6 m/s**; our
  aorta **4.86 m/s** (thick-wall (16.63): 4.71 m/s); E → ∞: c → ∞, the rigid tube of C09 (T23).
- **What it means:** the pulse is a wave of wall stretch; measuring its speed measures arterial stiffness (pulse-wave
  velocity). Fails at large α only slightly (C11) but at small α strongly (viscous waves are slow and damped).
- **Traps:** Eh/(2ρa₀) uses the radius; with the diameter it reads Eh/(ρd) (T17); c is not the blood's speed.

### D13 · Womersley's profile (16.76)–(16.82) — ★★★, 13 steps, in C09 (notebook · `womersley_profiles`)
- **Goal:** find the velocity profile of blood driven back and forth by an oscillating pressure gradient in a rigid tube. ·
  **Start:** $\mu\left[\dfrac{\partial^2u}{\partial r^2}+\dfrac1r\dfrac{\partial u}{\partial r}\right]-\rho\dfrac{\partial u}{\partial t}=\dfrac{\partial p}{\partial x}$ (16.72) — "unsteady Poiseuille: viscosity and inertia against the
  gradient". · **Plan:** • one harmonic of the gradient • a profile that oscillates with it • a Bessel equation with a constant
  forcing • regularity and no slip • write it with α. · **Tools:** complex amplitudes (P176, reminder); Bessel's equation, J₀ and
  Y₀ (P374, primer in C09); powers of i (P377, primer in C09); the Laplacian in cylindrical coordinates (P186, reminder). ·
  **Assumptions:** rigid tube, long wave, axisymmetric, convective inertia dropped (R11); convention $e^{+i\omega t}$ (trap T1).
- **Steps:**
  1. **did** Impose one harmonic of the gradient · **tex** $\dfrac{\partial p}{\partial x}=-Ge^{i\omega t}$ (16.73) · **why** Any periodic gradient is a
     sum of harmonics and (16.72) is linear, so one harmonic is enough (N57); complex form (P176), the real part −G cos ωt is meant. ·
     **plain** The pressure gradient swings back and forth with amplitude G.
  2. **did** Try a profile that swings with it · **tex** $u(r,t)=U(r)e^{i\omega t}$ (16.75) · **why** A linear equation forced at frequency ω
     settles into motion at that frequency; the complex U(r) carries amplitude and phase at each radius. · **plain** Every radius
     oscillates at the heart's frequency, each with its own size and delay.
  3. **did** Substitute and cancel $e^{i\omega t}$ · **tex** $\mu\left[U''+\dfrac{U'}r\right]-i\omega\rho U=-G$ · **why** $\partial/\partial t\to i\omega$;
     $e^{i\omega t}\neq0$ divides out of every term. · **plain** An ordinary differential equation in r alone.
  4. **did** Divide by μ · **tex** $\dfrac{d^2U}{dr^2}+\dfrac1r\dfrac{dU}{dr}-\dfrac{i\omega\rho}{\mu}U=-\dfrac G\mu$ (16.76, right side corrected) · **why**
     μ > 0. The book prints +G/μ (slip S6); its own solution (16.77) needs −G/μ. · **plain** The equation for the complex profile.
  5. **did** Find the constant particular solution · **tex** $U_p=\dfrac{G}{i\omega\rho}$ · **why** A constant has U'' = U' = 0, leaving
     $-(i\omega\rho/\mu)U_p=-G/\mu$. · **plain** Far from any wall the blood just accelerates under G: ρ iω U = G.
  6. **did** Write the rest in Bessel's form · **tex** $U''+\dfrac{U'}r+\lambda^2U=0,\quad\lambda^2\equiv-\dfrac{i\omega\rho}{\mu}$ · **why** The homogeneous
     part of (16.76); with z = λr it becomes $y''+y'/z+y=0$, Bessel's equation of order 0 (P374). · **plain** The wall's
     influence obeys Bessel's equation.
  7. **did** Take the square root of −i · **tex** $\lambda=i^{3/2}\sqrt{\dfrac{\omega\rho}{\mu}}$ · **why** −i = i³, so λ² = i³ωρ/μ and
     $i^{3/2}=e^{3i\pi/4}$ (P377). The other root −λ gives the same J₀, which is even. · **plain** The Bessel argument is complex.
  8. **did** Write the general solution · **tex** $U=C_1J_0(\lambda r)+C_2Y_0(\lambda r)+\dfrac{G}{\omega\rho i}$ (16.77) · **why** J₀ and Y₀ are the
     two independent solutions of Bessel's equation (P374); add the particular solution of step 5. The book writes λ as
     $i\sqrt{i\omega\rho/\mu}$, the same number. · **plain** Two Bessel shapes plus the core's acceleration.
  9. **did** Keep the velocity finite on the axis · **tex** $C_2=0$ · **why** Y₀(z) → −∞ as z → 0 (P374); a real flow is finite at
     r = 0. · **plain** Only J₀ survives.
  10. **did** Apply no slip at the wall · **tex** $C_1J_0\left(i^{3/2}\sqrt{\omega\rho/\mu}\,a\right)+\dfrac{G}{\omega\rho i}=0$ (16.78) · **why** Rigid wall:
      U(a) = 0. · **plain** One equation for C₁.
  11. **did** Name the argument with α and solve · **tex** $C_1=\dfrac{iG}{\omega\rho}\dfrac1{J_0(i^{3/2}\alpha)}$ (16.79) · **why**
      $a\sqrt{\omega\rho/\mu}=a\sqrt{\omega/\nu}=\alpha$ by (16.12); $-1/i=i$. · **plain** The Womersley number is the Bessel
      argument at the wall.
  12. **did** Assemble the profile · **tex** $U(r)=-\dfrac{iG}{\omega\rho}\left(1-\dfrac{J_0(i^{3/2}\alpha\,r/a)}{J_0(i^{3/2}\alpha)}\right)$ (16.80) · **why** Put C₁ into
      (16.77) with $1/i=-i$ and $\lambda r=i^{3/2}\alpha r/a$; the ratio is $F_1(\alpha)=\dfrac{J_0(i^{3/2}\alpha\,r/a)}{J_0(i^{3/2}\alpha)}$ (16.81). · **plain** The core's
      acceleration, cut down near the wall by the Bessel ratio.
  13. **did** Write the prefactor with α · **tex** $u(r,t)=\dfrac{Ga^2}{i\mu\alpha^2}(1-F_1(\alpha))e^{i\omega t}$ (16.82) · **why** α² = a²ωρ/μ gives
      $a^2/(\mu\alpha^2)=1/(\omega\rho)$, and $-i/(\omega\rho)=1/(i\omega\rho)$; multiply by $e^{i\omega t}$ (16.75). · **plain**
      Womersley's oscillating profile.
- **Result:** $u(r,t)=\dfrac{Ga^2}{i\mu\alpha^2}\left(1-\dfrac{J_0(i^{3/2}\alpha r/a)}{J_0(i^{3/2}\alpha)}\right)e^{i\omega t}$ (16.82) — "one formula for every α".
- **Check:** units — Ga²/μ: (Pa/m)·m²/(Pa s) = m/s ✓; the sympy cell below gives a zero residual of (16.76, corrected) and U(a) = 0;
  α → 0 gives Poiseuille (D14), α → ∞ a plug with a Stokes layer (D15); **expect** at α = 3 the centreline amplitude is 0.522 of
  Ga²/4μ and lags the gradient by 64.3°; the from-scratch series J₀ agrees with `scipy.special.jv` to 10⁻¹².
- **sympy check (`check_src`, executed):**
  ```python
  import sympy as sp                                                    # symbolic algebra
  r, a, w, rho, mu, G = sp.symbols('r a omega rho mu G', positive=True)  # radius, tube radius, ω, ρ, μ, G
  lam = sp.exp(3*sp.I*sp.pi/4)*sp.sqrt(w*rho/mu)                        # step 7: λ = i^{3/2}√(ωρ/μ)
  U = -sp.I*G/(w*rho)*(1 - sp.besselj(0, lam*r)/sp.besselj(0, lam*a))   # step 12: the profile (16.80)
  lhs = sp.diff(U, r, 2) + sp.diff(U, r)/r - sp.I*w*rho/mu*U            # left side of (16.76)
  print(sp.simplify(sp.expand(lhs + G/mu)))                             # residual with the corrected right side −G/μ → 0
  print(sp.simplify(U.subs(r, a)))                                      # no slip at the wall → 0
  ```
  printed: `0` and `0`.
- **What it means:** α is the only parameter: the same formula gives the breathing parabola of a small artery and the flat,
  lagging plug of the aorta. Fails for elastic walls (C11), for non-Newtonian blood at low shear and near branches.
- **Traps:** slip S6; the time convention $e^{+i\omega t}$ (T1); dropping Y₀ needs the axis inside the fluid; for α ≳ 30 use
  the scaled `jve` (T25).

### D14 · The small-α (quasi-steady) limit (16.87) — ★★, 6 steps, in C09 (notebook · `womersley_profiles`)
- **Goal:** show that slow oscillations (α ≪ 1) give Poiseuille flow at every instant. · **Start:** (16.82) and the series of J₀
  (P374). · **Plan:** • first terms of J₀ • square the argument • expand the ratio • substitute. · **Tools:** J₀'s power series
  (P374); $1/(1+\varepsilon)\approx1-\varepsilon$ for complex ε (P342, reminder). · **Assumptions:** α ≪ 1 (step 3).
- **Steps:**
  1. **did** Write the first terms of J₀ · **tex** $J_0(z)=1-\dfrac{z^2}{4}+\dfrac{z^4}{64}-\dots$ · **why** The power series
     $J_0(z)=\sum(-1)^m(z/2)^{2m}/(m!)^2$ (P374); small |z| needs only the first terms. · **plain** Near zero J₀ is a downward
     parabola.
  2. **did** Square the Bessel argument · **tex** $z=i^{3/2}\alpha\dfrac ra:\quad z^2=-i\alpha^2\dfrac{r^2}{a^2}$ · **why** $(i^{3/2})^2=i^3=-i$. ·
     **plain** The squared argument is purely imaginary.
  3. **did** Expand numerator and denominator of F₁ · **tex** $F_1\approx\dfrac{1+i\alpha^2r^2/4a^2}{1+i\alpha^2/4}$ · **why** Steps 1–2 with α ≪ 1;
     we keep the z² term because the 1's cancel in 1 − F₁, so z² is the leading term. · **plain** Both Bessel values are 1 plus a
     small imaginary correction.
  4. **did** Divide to first order in α² · **tex** $F_1\approx1-\dfrac{i\alpha^2}{4}\left(1-\dfrac{r^2}{a^2}\right)$ · **why** $1/(1+\varepsilon)\approx1-\varepsilon$
     with ε = iα²/4 (P342); multiply out and drop α⁴ terms. · **plain** F₁ is almost 1, short by a small imaginary parabola.
  5. **did** Form 1 − F₁ · **tex** $1-F_1\approx\dfrac{i\alpha^2}{4}\left(1-\dfrac{r^2}{a^2}\right)$ · **why** Subtract from 1. · **plain** What drives the
     flow is that small parabola.
  6. **did** Substitute into (16.82) · **tex** $u\approx\dfrac{Ga^2}{4\mu}\left(1-\dfrac{r^2}{a^2}\right)e^{i\omega t}$ (16.87) · **why**
     $\dfrac{Ga^2}{i\mu\alpha^2}\cdot\dfrac{i\alpha^2}4=\dfrac{Ga^2}{4\mu}$: the i and α² cancel. For harmonic n use G_n and α_n = α√n. · **plain**
     Poiseuille's parabola, in phase with the gradient.
- **Result:** $u_n(r)\approx\dfrac{G_na^2}{4\mu}\left(1-\dfrac{r^2}{a^2}\right)$ (16.87).
- **Check:** compare $u_{max}=\dfrac{\Delta p\,a^2}{4\mu L}$ (16.16) with G = Δp/L ✓; **expect** at α = 0.5 the exact centreline amplitude is 0.9990 of
  Ga²/4μ and lags by 2.7° (the error of (16.87) shrinks like α⁴ in amplitude and α² in phase).
- **What it means:** below α ≈ 1 the viscous layer fills the vessel and the flow keeps step with the pressure — arterioles,
  capillaries. Fails above α ≈ 1.
- **Traps:** stopping J₀ at "1" gives 1 − F₁ = 0 — keep the z² term.

### D15 · The large-α limit (16.88) by a thin wall layer (rule b: the book states it) — ★★, 8 steps, in C09 (notebook · `womersley_profiles`)
- **Goal:** show that fast oscillations (α ≫ 1) give a flat core joined to the wall by a Stokes layer, without Bessel
  asymptotics. · **Start:** $\dfrac{d^2U}{dr^2}+\dfrac1r\dfrac{dU}{dr}-\dfrac{i\omega\rho}{\mu}U=-\dfrac G\mu$ (16.76, right side corrected). · **Plan:** • measure from the wall • drop the curvature term
  in a thin layer • solve a constant-coefficient equation • keep the root that decays. · **Tools:** orders of smallness and a
  dominant balance (P68, reminder); trying $e^{\lambda y}$ (P44, reminder); √i (P377); Stokes' second problem (Ch. 8 C11 recap). ·
  **Assumptions:** α ≫ 1, so the layer is thin (step 2).
- **Steps:**
  1. **did** Measure distance from the wall · **tex** $y=a-r,\quad\dfrac{d}{dr}=-\dfrac{d}{dy}$ · **why** At large α everything happens close to
     the wall; y is the natural coordinate there. · **plain** Look outward from the wall.
  2. **did** Compare the two viscous terms · **tex** $\left|\dfrac1r\dfrac{dU}{dr}\right|\Big/\left|\dfrac{d^2U}{dr^2}\right|\sim\dfrac{\delta}{a}\ll1$ · **why** If U changes over a thin
     layer of thickness δ, then U″ ~ U/δ² while U′/r ~ U/(aδ); their ratio is δ/a (a dominant balance, P68). · **plain** In a
     thin layer the wall looks flat.
  3. **did** Drop the curvature term · **tex** $\dfrac{d^2U}{dy^2}-\dfrac{i\omega}{\nu}U\approx-\dfrac G\mu$ · **why** Step 2; $d^2/dr^2=d^2/dy^2$ and
     ρ/μ = 1/ν. · **plain** The equation of an oscillating flat-wall layer (Ch. 8).
  4. **did** Keep the core solution · **tex** $U_p=\dfrac{G}{i\omega\rho}$ · **why** The same constant as D13 step 5: far from the wall
     the blood only accelerates. · **plain** The core is a plug.
  5. **did** Try an exponential for the rest · **tex** $U_h=Ce^{\lambda y},\quad\lambda^2=\dfrac{i\omega}{\nu}$ · **why** A constant-coefficient
     linear equation: $e^{\lambda y}$ works when λ² equals the coefficient (P44). · **plain** The wall's effect decays
     exponentially.
  6. **did** Take the root that decays into the core · **tex** $\lambda=-(1+i)\sqrt{\dfrac{\omega}{2\nu}}$ · **why** $\sqrt i=(1+i)/\sqrt2$ (P377), so
     $\lambda=\pm(1+i)\sqrt{\omega/2\nu}$; the + root grows without bound into the core, which is impossible. · **plain** The
     layer decays and turns its phase over the same length.
  7. **did** Apply no slip · **tex** $C=-\dfrac{G}{i\omega\rho}$ · **why** U(y = 0) = U_p + C = 0. · **plain** At the wall the layer cancels the
     plug.
  8. **did** Assemble and write α in · **tex** $u\approx\dfrac{Ga^2}{i\mu\alpha^2}\left\{1-\exp\left[-\sqrt{\dfrac{\omega}{2\nu}}(1+i)(a-r)\right]\right\}e^{i\omega t}$ (16.88) · **why**
     G/(iωρ) = Ga²/(iμα²) as in D13 step 13; for harmonic n replace ω by nω (the book keeps ω — slip S7). · **plain** A flat plug
     lagging the gradient by 90°, joined to the wall by a layer of thickness $\delta_e=\sqrt{2\nu/\omega}=\sqrt2a/\alpha$.
- **Result:** (16.88).
- **Check:** units ✓; at the wall the bracket is 0 ✓; far from it, u → Ga²/(iμα²)e^{iωt}: magnitude G/(ωρ), a quarter cycle late
  (1/i = e^{−iπ/2}); seen from the plug's frame this is ch08's $u(y,t)=U\exp\{-y\sqrt{\omega/2\nu}\}\cos(\omega t-y\sqrt{\omega/2\nu})$ (8.38) with the wall oscillating against the plug (`LAM.stokes_second_problem`);
  **expect** at α = 20, r = 0.9a the asymptote is 5 % below the exact profile (the curvature correction ~δ/r) — good near the
  wall, approximate further in.
- **What it means:** in the aorta the blood moves as a solid slug pushed by inertia; viscosity acts only in a skin about a
  millimetre thick, where reversal starts first. Fails below α ≈ 5 (the layer is not thin).
- **Traps:** keep the decaying root; harmonic n has $\sqrt{n\omega/2\nu}$ (S7).

### D16 · The flow rate (16.90) and F₂ (16.91) — ★★, 7 steps, in C10 (notebook · `womersley_flow_phase`)
- **Goal:** add up Womersley's profile over the cross-section to get the oscillating flow rate. · **Start:**
  $u(r,t)=\dfrac{Ga^2}{i\mu\alpha^2}(1-F_1(\alpha))e^{i\omega t}$ (16.82). · **Plan:** • rings again • the constant part • the Bessel part by a substitution and an
  identity • collect. · **Tools:** substitution in an integral (P106, reminder); the Bessel identities (P376, primer in C10). ·
  **Assumptions:** as D13.
- **Steps:**
  1. **did** Integrate over rings · **tex** $Q=\dfrac{2\pi Ga^2}{i\mu\alpha^2}e^{i\omega t}\displaystyle\int_0^a(1-F_1)\,r\,dr$ · **why** $Q=\int_0^au\,2\pi r\,dr$ (D04 step 1);
     everything not depending on r comes out of the integral. · **plain** Flow = profile summed over rings.
  2. **did** Integrate the constant · **tex** $\displaystyle\int_0^ar\,dr=\dfrac{a^2}2$ · **why** Elementary. · **plain** The plug part covers the
     whole area.
  3. **did** Substitute $s=\lambda r$ in the Bessel term · **tex** $\displaystyle\int_0^aJ_0(\lambda r)\,r\,dr=\dfrac1{\lambda^2}\int_0^{\lambda a}sJ_0(s)\,ds,\quad\lambda=\dfrac{i^{3/2}\alpha}a$ · **why**
     r = s/λ, dr = ds/λ (P106); the limits become 0 and λa = i^{3/2}α. · **plain** A standard Bessel integral appears.
  4. **did** Use the Bessel integral identity · **tex** $\displaystyle\int_0^{z}sJ_0(s)\,ds=zJ_1(z)$ · **why** P376: $\dfrac{d}{dz}[zJ_1(z)]=zJ_0(z)$, so
     zJ₁ is the antiderivative, and it vanishes at 0. · **plain** The area under sJ₀ is zJ₁.
  5. **did** Collect the Bessel term · **tex** $\displaystyle\int_0^aF_1\,r\,dr=\dfrac{a^2}2\cdot\dfrac{2J_1(i^{3/2}\alpha)}{i^{3/2}\alpha J_0(i^{3/2}\alpha)}=\dfrac{a^2}2F_2(\alpha)$ (16.91) · **why** Divide steps
     3–4 by J₀(i^{3/2}α): $\dfrac{\lambda aJ_1}{\lambda^2J_0}=\dfrac{a^2J_1}{i^{3/2}\alpha J_0}$; factor a²/2. · **plain** F₂ is the area-average of the Bessel
     ratio F₁.
  6. **did** Put the pieces together · **tex** $Q=\dfrac{2\pi Ga^2}{i\mu\alpha^2}e^{i\omega t}\cdot\dfrac{a^2}2\left(1-F_2(\alpha)\right)$ · **why** Steps 2 and 5 in
     step 1. · **plain** The flux is the plug flux reduced by the averaged wall effect.
  7. **did** Simplify · **tex** $Q(t)=\dfrac{\pi a^4}{i\mu\alpha^2}G(1-F_2(\alpha))e^{i\omega t}$ (16.90) · **why** $2\pi\cdot a^2/2=\pi a^2$. · **plain**
     Womersley's oscillating flow rate.
- **Result:** (16.90) with $F_2(\alpha)=\dfrac{2J_1(i^{3/2}\alpha)}{i^{3/2}\alpha J_0(i^{3/2}\alpha)}$ (16.91).
- **Check:** α → 0: $F_2\approx1-\dfrac{i\alpha^2}8$ (16.93), so $Q\to\dfrac{\pi Ga^4}{8\mu}e^{i\omega t}$ — Poiseuille's $|Q_0|=\dfrac{\pi Ga^4}{8\mu}$ (16.95) ✓; the
  from-scratch trapezoid of 2πru equals `ch16.womersley_flow_rate` to 10⁻⁸; **expect** |Q|/|Q₀| = 0.548 at α = 3 and 0.105 at
  α = 8.
- **What it means:** at large α the same gradient buys far less flow, because most of the push goes into accelerating the
  blood back and forth. Fails as D13.
- **Traps:** the constant term integrates to the full area πa²; collect the $i^{3/2}$ factors carefully.

### D17 · Amplitude (16.99) and phase lag (16.100) — ★★, 7 steps, in C10 (notebook · `womersley_flow_phase`)
- **Goal:** read off how big the oscillating flow is and how late it comes after the pressure gradient. · **Start:**
  $Q(t)=\dfrac{\pi a^4}{i\mu\alpha^2}G(1-F_2(\alpha))e^{i\omega t}$ (16.90). · **Plan:** • name 1 − F₂ • multiply out with Euler's formula • write the real part as one shifted
  cosine. · **Tools:** Euler's formula (P45, reminder); the complex plane (P153, reminder); amplitude–phase form (written out in
  step 5); `np.arctan2` (P70, reminder). · **Assumptions:** the physical gradient is the real part, $-\partial p/\partial x=G\cos\omega t$.
- **Steps:**
  1. **did** Name the complex factor · **tex** $1-F_2(\alpha)=Z(\alpha)=X(\alpha)+iY(\alpha)$ (16.97) · **why** A shorthand; X and Y are real
     numbers that depend only on α. · **plain** One complex number per α.
  2. **did** Replace 1/i by −i · **tex** $Q=\dfrac{\pi Ga^4}{\mu\alpha^2}(-i)(X+iY)e^{i\omega t}$ · **why** $1/i=-i$ (multiply top and bottom by i). ·
     **plain** A quarter-turn of the complex factor.
  3. **did** Multiply out with Euler's formula · **tex** $Q=\dfrac{\pi Ga^4}{\mu\alpha^2}(Y-iX)(\cos\omega t+i\sin\omega t)$ · **why** $-i(X+iY)=Y-iX$;
     $e^{i\omega t}=\cos\omega t+i\sin\omega t$ (P45). · **plain** Two complex numbers to multiply.
  4. **did** Separate real and imaginary parts · **tex** $Q=\dfrac{\pi Ga^4}{\mu\alpha^2}\left\{[Y\cos\omega t+X\sin\omega t]-i[X\cos\omega t-Y\sin\omega t]\right\}$ (16.98) · **why**
     $(Y-iX)(\cos+i\sin)=Y\cos+X\sin+i(Y\sin-X\cos)$. · **plain** The real part is the flow that answers $G\cos\omega t$.
  5. **did** Write the real part as one shifted cosine · **tex** $Y\cos\omega t+X\sin\omega t=\sqrt{X^2+Y^2}\cos(\omega t-\phi)$ · **why**
     $R\cos(\omega t-\phi)=R\cos\phi\cos\omega t+R\sin\phi\sin\omega t$; match $R\cos\phi=Y$, $R\sin\phi=X$. · **plain** A sine and a
     cosine of the same frequency add to one cosine, shifted.
  6. **did** Read off the amplitude · **tex** $|Q|=\dfrac{\pi Ga^4}{\mu\alpha^2}\sqrt{X^2+Y^2}$ (16.99) · **why** The prefactor times R of step 5. ·
     **plain** How much flow one unit of gradient buys.
  7. **did** Read off the lag · **tex** $\tan\phi=\dfrac XY$ (16.100), $\ \phi=\operatorname{arctan2}(X,Y)$ · **why** From $R\sin\phi=X$, $R\cos\phi=Y$;
     arctan2 picks the right quadrant (P70). φ is how far the flow's peak trails the gradient's peak. · **plain** The flow arrives
     φ/ω seconds late.
- **Result:** (16.99), (16.100); and (ours) $\dfrac{|Q|}{|Q_0|}=\dfrac{8\sqrt{X^2+Y^2}}{\alpha^2}$ with $|Q_0|=\dfrac{\pi Ga^4}{8\mu}$ (16.95).
- **Check:** α → 0: Z ≈ iα²/8, so X → 0, φ → 0 and |Q|/|Q₀| → 1 ✓; α → ∞: Z → 1, φ → 90°, |Q|/|Q₀| ≈ 8/α² ✓; **expect** at α = 8:
  X = 0.8228, Y = 0.1608, |Q|/|Q₀| = 0.1048, φ = 78.9° — at 66 beats/min the flow peaks **199 ms** after the gradient.
- **What it means:** small α: a driven damper — viscosity limits the flow, in phase with the push. Large α: a driven mass —
  inertia limits it, a quarter cycle late, amplitude ∝ 1/α². Fails as D13.
- **Traps:** tan φ = X/Y, not Y/X; φ is the lag of Q behind $-\partial p/\partial x$, not behind p.

### D18 · The long-wave normal-mode reduction (16.116)–(16.125) — ★★, 9 steps, in C11 (notebook)
- **Goal:** reduce the viscous flow in an elastic tube to ordinary differential equations in r for one travelling wave. ·
  **Start:** $\dfrac{\partial u}{\partial x}+\dfrac1r\dfrac{\partial(rv)}{\partial r}=0$ (16.101), the two linearised momentum equations (16.102)–(16.103), and the wall conditions
  $u\rvert_{r=a_0}=\dfrac{\partial\xi}{\partial t}$ (16.114), $v\rvert_{r=a_0}=\dfrac{\partial\eta}{\partial t}$ (16.115). · **Plan:** • normal modes with a complex wavenumber • size of v • drop what
  long waves make small • substitute the modes. · **Tools:** complex amplitudes (P176); a complex wavenumber (P385, primer in
  C11); orders of smallness (P68). · **Assumptions:** long waves, $a_0|k|\ll1$ (steps 3–5); small disturbances (start).
- **Steps:**
  1. **did** Choose travelling-wave modes · **tex** $(u,v,p,\xi,\eta)=(\hat u(r),\hat v(r),\hat p,\hat\xi,\hat\eta)\,e^{i(kx-\omega t)}$ (16.116) · **why** The equations
     are linear with coefficients that do not depend on x or t, so each Fourier mode evolves alone (P176). The book now uses
     $e^{i(kx-\omega t)}$, the opposite sign to C09 (trap T1). · **plain** One wave of frequency ω and wavenumber k.
  2. **did** Let the wavenumber be complex · **tex** $k=k_1+ik_2:\quad e^{i(kx-\omega t)}=e^{-k_2x}e^{i(k_1x-\omega t)}$ · **why** P385: an imaginary part
     of k makes the amplitude fall with distance; speed $c=\omega/k_1$, $|k|=\sqrt{k_1^2+k_2^2}$ (the book prints $\sqrt{k_1+k_2}$, S13). ·
     **plain** A wave that travels at ω/k₁ and fades over 1/k₂.
  3. **did** Size the radial velocity from continuity · **tex** $\left|\dfrac{v}{u}\right|=O(|ak|)$ (16.117) · **why** In (16.101) $\partial u/\partial x\sim|k|u$
     and $(1/r)\partial(rv)/\partial r\sim v/a$ must balance. · **plain** For long waves the blood moves mostly along the tube.
  4. **did** Drop the axial diffusion · **tex** $\left|\dfrac{\partial^2u}{\partial x^2}\right|\Big/\left|\dfrac{\partial^2u}{\partial r^2}\right|=O(a^2|k|^2)\ll1$ · **why** ∂²/∂x² brings k²,
     ∂²/∂r² brings at least 1/a²; long waves (P68). · **plain** Viscosity acts across the tube, not along it.
  5. **did** Drop the radial pressure change · **tex** $\dfrac{\partial p}{\partial r}=0$ (16.120) · **why** Every term of the radial momentum (16.103)
     carries v = O(ak)u, so the pressure change across the radius is O(a²k²) of the change along the tube. · **plain** The
     pressure is the same over each cross-section.
  6. **did** Write the reduced axial momentum · **tex** $\rho\dfrac{\partial u}{\partial t}=-\dfrac{\partial p}{\partial x}+\mu\left(\dfrac{\partial^2u}{\partial r^2}+\dfrac1r\dfrac{\partial u}{\partial r}\right)$ (16.119) ·
     **why** (16.102) without the term of step 4. · **plain** The rigid-tube equation (16.72) again — now with a pressure that
     travels.
  7. **did** Put the modes into continuity · **tex** $\dfrac{d(r\hat v)}{dr}=-ikr\hat u$ (16.122) · **why** $\partial/\partial x\to ik$ in (16.101);
     multiply by r. · **plain** Radial flow fills in for the change of axial flow.
  8. **did** Put the modes into momentum · **tex** $\dfrac{d^2\hat u}{dr^2}+\dfrac1r\dfrac{d\hat u}{dr}+\dfrac{i\omega\rho}{\mu}\hat u=\dfrac{ik\hat p}{\mu}$ (16.123) · **why**
     $\partial/\partial t\to-i\omega$, $\partial/\partial x\to ik$ in (16.119); divide by μ and move terms. · **plain** C09's Bessel-type
     equation, with +iωρ/μ because of the other time convention.
  9. **did** Apply the wall conditions at the undisturbed wall · **tex** $\hat u(a_0)=-i\omega\hat\xi,\quad\hat v(a_0)=-i\omega\hat\eta$ (16.124), (16.125) ·
     **why** (16.114)–(16.115): the fluid moves with the wall; for small displacements the condition is applied at r = a₀ instead
     of a₀ + η (the linearised boundary condition of Ch. 7); $\partial/\partial t\to-i\omega$. · **plain** The blood at the wall
     moves exactly as the wall does.
- **Result:** (16.122)–(16.125) with ∂p/∂r = 0.
- **Check:** with k real and the wall rigid ($\hat\xi=\hat\eta=0$) (16.123) is the complex conjugate of D13's equation (opposite
  time convention); every dropped term is $O(a^2k^2)$ or $O(ak)$.
- **What it means:** the long-wave idea keeps the cross-section physics of C09 and adds a pressure that varies along the tube.
  Fails for short waves (λ comparable with a).
- **Traps:** T1 (the sign of iω flips); S13 the modulus of k.

### D19 · The reduced wall equations (16.128)–(16.129) — ★★★, 10 steps, in C11 (notebook)
- **Goal:** do the "considerable algebra" the book skips: simplify the wall's two equations of motion and eliminate its axial
  displacement. · **Start:** $-\rho_wh\omega^2\hat\eta=\hat p-2\mu\dfrac{d\hat v}{dr}\Big\rvert_{a_0}-\dfrac{Eh}{1-\hat\nu^2}\left(\dfrac{\hat\eta}{a_0^2}+\dfrac{i\hat\nu k}{a_0}\hat\xi\right)$ (16.126) and
  $-\rho_wh\omega^2\hat\xi=-\mu\left(\dfrac{d\hat u}{dr}+ik\hat v\right)\Big\rvert_{a_0}+\dfrac{Eh}{1-\hat\nu^2}\left(-k^2\hat\xi+\dfrac{i\hat\nu k}{a_0}\hat\eta\right)$ (16.127). · **Plan:** • size and drop wall inertia • drop two small fluid
  stresses • solve the axial equation for ξ̂ • substitute into the radial one. · **Tools:** Poisson's ratio (P379); orders of
  smallness (P68); linear elimination. · **Assumptions:** thin wall h ≪ a₀, ρ_w ~ ρ, long waves a₀k ≪ 1, α ≫ 1 (steps 1–4).
- **Steps:**
  1. **did** Compare wall inertia with wall stiffness · **tex** $\dfrac{\rho_wh\omega^2}{Eh/a_0^2}=\dfrac{\rho_w}{\rho}\dfrac{h}{a_0}\dfrac{(ka_0)^2}{2}$ · **why** Ratio of the left side
     of (16.126) to its elastic term; with $E=2\rho a_0c^2/h$ from (16.62) and ω = kc. Each factor is small. · **plain** A thin
     wall's own mass hardly matters.
  2. **did** Drop wall inertia in both equations · **tex** $0=\hat p-2\mu\dfrac{d\hat v}{dr}\Big\rvert_{a_0}-\dfrac{Eh}{1-\hat\nu^2}\left(\dfrac{\hat\eta}{a_0^2}+\dfrac{i\hat\nu k}{a_0}\hat\xi\right)$ · **why**
     Step 1; the same estimate holds for (16.127). · **plain** The wall is in balance at every instant.
  3. **did** Drop the viscous normal stress · **tex** $2\mu\left|\dfrac{d\hat v}{dr}\right|\ll|\hat p|$ · **why** $v=O(ak)u$ (16.117), so
     $2\mu\,dv/dr\sim2\mu|k|u$, while $p\sim\rho cu$; their ratio $2\nu|k|/c$ is small for α ≫ 1 and ak ≪ 1. · **plain** The
     wall feels the blood's pressure, not its viscous squeeze.
  4. **did** Keep only the main shear · **tex** $\mu\left(\dfrac{d\hat u}{dr}+ik\hat v\right)\approx\mu\dfrac{d\hat u}{dr}$ · **why** $|k\hat v|\sim a|k|^2|\hat u|$ while
     $|d\hat u/dr|\ge|\hat u|/a$: the ratio is O(a²k²). The text names μ∂v/∂x as the survivor (slip S14). · **plain** The wall is
     dragged by the change of u across the radius.
  5. **did** Write the reduced axial equation · **tex** $0=-\mu\dfrac{d\hat u}{dr}\Big\rvert_{a_0}+\dfrac{Eh}{1-\hat\nu^2}\left(-k^2\hat\xi+\dfrac{i\hat\nu k}{a_0}\hat\eta\right)$ · **why** (16.127)
     with steps 2 and 4. · **plain** The blood's drag stretches the wall lengthwise.
  6. **did** Solve it for the axial displacement · **tex** $\hat\xi=\dfrac{i\hat\nu}{ka_0}\hat\eta-\dfrac{1-\hat\nu^2}{Ehk^2}\mu\dfrac{d\hat u}{dr}\Big\rvert_{a_0}$ (16.129) · **why**
     Move the shear term across and divide by $-k^2Eh/(1-\hat\nu^2)$; $(i\hat\nu k/a_0)/k^2=i\hat\nu/(ka_0)$. · **plain** The wall
     slides because the hoop stretch squeezes it lengthwise (Poisson) and because the blood drags it.
  7. **did** Write the reduced radial equation for p̂ · **tex** $\hat p=\dfrac{Eh}{1-\hat\nu^2}\left(\dfrac{\hat\eta}{a_0^2}+\dfrac{i\hat\nu k}{a_0}\hat\xi\right)$ · **why** Step 2 with
     step 3, solved for p̂. The book's (16.108) prints 1 − ν̂ (slip S12). · **plain** The pressure is held by the hoop and the
     lengthwise strain.
  8. **did** Multiply (16.129) by $i\hat\nu k/a_0$ · **tex** $\dfrac{i\hat\nu k}{a_0}\hat\xi=-\dfrac{\hat\nu^2}{a_0^2}\hat\eta-\dfrac{i\hat\nu(1-\hat\nu^2)}{a_0Ehk}\mu\dfrac{d\hat u}{dr}\Big\rvert_{a_0}$ · **why** $i\cdot i=-1$
     in the first term; k/k² = 1/k in the second. · **plain** The piece of step 7 that contains ξ̂.
  9. **did** Add $\hat\eta/a_0^2$ · **tex** $\dfrac{\hat\eta}{a_0^2}+\dfrac{i\hat\nu k}{a_0}\hat\xi=\dfrac{(1-\hat\nu^2)\hat\eta}{a_0^2}-\dfrac{i\hat\nu(1-\hat\nu^2)}{a_0Ehk}\mu\dfrac{d\hat u}{dr}\Big\rvert_{a_0}$ · **why**
     $1-\hat\nu^2$ now multiplies both terms. · **plain** The bracket of step 7, worked out.
  10. **did** Multiply by $Eh/(1-\hat\nu^2)$ · **tex** $\hat p=\dfrac{Eh}{a_0^2}\hat\eta-\dfrac{i\hat\nu}{a_0k}\mu\dfrac{d\hat u}{dr}\Big\rvert_{a_0}$ (16.128) · **why** Step 7 with step 9; the
      factor $1-\hat\nu^2$ cancels exactly. · **plain** Hoop stiffness holds the pressure; a small Poisson term couples it to the
      blood's shear.
- **Result:** (16.128) and (16.129).
- **Check:** ν̂ = 0: $\hat p=Eh\hat\eta/a_0^2$, the linear tube law (16.58) written with η ($A-A_0\approx2\pi a_0\eta$) ✓; units of
  each term Pa ✓.
- **sympy check (`check_src`, executed):**
  ```python
  import sympy as sp                                                       # symbolic algebra
  eta, xi, p, k, E, h, a0, nu, mu, Up = sp.symbols('eta xi p k E h a0 nuhat mu Uprime')   # Up = dû/dr at the wall
  D = E*h/(1 - nu**2)                                                      # membrane stiffness Eh/(1 − ν̂²)
  axial = sp.Eq(0, -mu*Up + D*(-k**2*xi + sp.I*nu*k*eta/a0))               # step 5: reduced (16.127)
  xi_sol = sp.solve(axial, xi)[0]                                          # step 6: solve for ξ̂
  print(sp.simplify(xi_sol - (sp.I*nu/(k*a0)*eta - (1 - nu**2)/(E*h*k**2)*mu*Up)))   # equals (16.129) → 0
  radial = sp.Eq(0, p - D*(eta/a0**2 + sp.I*nu*k*xi/a0))                   # step 7: reduced (16.126)
  p_sol = sp.solve(radial.subs(xi, xi_sol), p)[0]                          # steps 8–10: substitute and solve for p̂
  print(sp.simplify(p_sol - (E*h/a0**2*eta - sp.I*nu/(a0*k)*mu*Up)))      # equals (16.128) → 0
  ```
  printed: `0` and `0`.
- **What it means:** the wall is a massless spring whose stretch holds the pressure; viscosity reaches it only through the
  shear at the wall. Fails for thick or heavy walls.
- **Traps:** S12 (the missing square), S14 (which shear survives); the 1 − ν̂² cancels only after ξ̂ is eliminated.

### D20 · The mode solution (16.130)–(16.133) — ★★★, 10 steps, in C11 (notebook)
- **Goal:** solve the axial-velocity equation with the wall's pressure law, and get the radial velocity from continuity. ·
  **Start:** (16.123) and (16.128); (16.122). · **Plan:** • substitute p̂ • homogeneous Bessel part • constant particular part
  that contains the wall slope • integrate continuity. · **Tools:** Bessel's equation (P374); $dJ_0/dz=-J_1$ and
  $\int sJ_0\,ds=zJ_1$ (P376); principal square root (P159, reminder). · **Assumptions:** finite û on the axis (step 3).
- **Steps:**
  1. **did** Substitute (16.128) into (16.123) · **tex** $\dfrac{d^2\hat u}{dr^2}+\dfrac1r\dfrac{d\hat u}{dr}+\dfrac{i\omega\rho}{\mu}\hat u=\dfrac{ik}{\mu}\dfrac{Eh}{a_0^2}\hat\eta+\dfrac{\hat\nu}{a_0}\dfrac{d\hat u}{dr}\Big\rvert_{a_0}$ (16.130) ·
     **why** $(ik/\mu)\cdot(-i\hat\nu\mu/(a_0k))=\hat\nu/a_0$ because $i\cdot(-i)=1$. · **plain** The forcing on the right does not depend
     on r — it is a constant.
  2. **did** Name the Bessel parameter · **tex** $\beta^2=\dfrac{i\omega}{\nu},\quad\beta=\sqrt{\dfrac{i\omega}{\nu}}=e^{i\pi/4}\dfrac{\alpha}{a_0}$ · **why** ρ/μ = 1/ν;
     the principal root of i is $e^{i\pi/4}$ (P159); α = a₀√(ω/ν). So $\hat\beta\equiv\beta a_0=e^{i\pi/4}\alpha$. · **plain** α again
     sets the Bessel argument.
  3. **did** Solve the homogeneous part · **tex** $\hat u_h=AJ_0(\beta r)$ · **why** $\hat u''+\hat u'/r+\beta^2\hat u=0$ is Bessel's equation in
     z = βr (P374); Y₀ is excluded because û is finite on the axis. · **plain** The wall layer's shape.
  4. **did** Find the constant particular part · **tex** $\hat u_p=\dfrac{\text{right side of (16.130)}}{\beta^2}$ · **why** For a constant,
     only $\beta^2\hat u_p$ survives on the left. · **plain** The core's response to the travelling pressure.
  5. **did** Differentiate at the wall · **tex** $\dfrac{d\hat u}{dr}\Big\rvert_{a_0}=-A\beta J_1(\beta a_0)$ · **why** $dJ_0/dz=-J_1$ (P376) and the chain rule;
     the constant part has no slope. · **plain** The wall slope depends on A alone.
  6. **did** Write the particular part out · **tex** $\hat u_p=\dfrac{k}{\omega}\dfrac{Eh}{\rho a_0^2}\hat\eta-\dfrac{\hat\nu}{\beta a_0}AJ_1(\beta a_0)$ · **why**
     $\dfrac{ik}{\mu}\dfrac{Eh\hat\eta}{a_0^2}\cdot\dfrac{\mu}{i\omega\rho}=\dfrac{kEh\hat\eta}{\omega\rho a_0^2}$; $\dfrac{\hat\nu}{a_0}(-A\beta J_1)\dfrac1{\beta^2}=-\dfrac{\hat\nu AJ_1}{\beta a_0}$. · **plain** The
     travelling pressure accelerates the core; the Poisson term corrects it.
  7. **did** Add the two parts · **tex** $\hat u(r)=AJ_0(\beta r)+\dfrac k\omega\dfrac{Eh}{\rho a_0^2}\hat\eta-\dfrac{\hat\nu}{\beta a_0}AJ_1(\beta a_0)$ (16.131) · **why** A linear ODE:
     general solution = homogeneous + particular. · **plain** The axial velocity of the mode.
  8. **did** Integrate continuity from the axis · **tex** $\hat v=-\dfrac{ik}{r}\displaystyle\int_0^rs\,\hat u(s)\,ds$ (16.132) · **why** Integrate (16.122) from 0 to
     r; $r\hat v=0$ on the axis. · **plain** Radial flow fills in the axial flow's change.
  9. **did** Integrate the Bessel term · **tex** $\displaystyle\int_0^rsJ_0(\beta s)\,ds=\dfrac r\beta J_1(\beta r)$ · **why** P376 after s → βs, as in D16
     steps 3–4. · **plain** The same Bessel identity again.
  10. **did** Integrate the constants and assemble · **tex** $\hat v(r)=-\dfrac{ikA}{\beta}J_1(\beta r)-\dfrac{ik^2}{\omega}\dfrac{Eh\hat\eta}{\rho a_0^2}\dfrac r2+\dfrac{ik\hat\nu}{\beta a_0}A\dfrac r2J_1(\beta a_0)$ (16.133) ·
      **why** $\int_0^rs\,ds=r^2/2$; multiply everything by −ik/r. · **plain** The radial velocity of the mode.
- **Result:** (16.131) and (16.133), with $\beta=\sqrt{i\omega/\nu}$.
- **Check:** the sympy cell below puts (16.131) back into (16.130) (residual 0) and (16.131) into (16.132) (equals (16.133));
  rigid wall and ν̂ = 0: û is C09's profile with the conjugate convention.
- **sympy check (`check_src`, executed):**
  ```python
  import sympy as sp                                                       # symbolic algebra
  r, w, rho, k, E, h, a0, nu, nuk, A, eta = sp.symbols('r omega rho k E h a0 nuhat nu_k A eta')   # nu_k = kinematic ν
  beta = sp.sqrt(sp.I*w/nuk)                                               # step 2: β = √(iω/ν)
  u = A*sp.besselj(0, beta*r) + k/w*E*h/(rho*a0**2)*eta - nu/(beta*a0)*A*sp.besselj(1, beta*a0)   # (16.131)
  slope = sp.diff(u, r).subs(r, a0)                                        # dû/dr at the wall (step 5)
  res = sp.diff(u, r, 2) + sp.diff(u, r)/r + sp.I*w/nuk*u - (sp.I*k/(rho*nuk)*E*h/a0**2*eta + nu/a0*slope)
  print(sp.simplify(sp.expand(res)))                                       # (16.130) residual → 0
  s = sp.symbols('s', positive=True)                                       # integration variable
  v = -sp.I*k/r*sp.integrate(s*u.subs(r, s), (s, 0, r))                   # (16.132)
  v133 = (-sp.I*k*A/beta*sp.besselj(1, beta*r) - sp.I*k**2/w*E*h*eta/(rho*a0**2)*r/2
          + sp.I*k*nu/(beta*a0)*A*r/2*sp.besselj(1, beta*a0))              # (16.133)
  print(sp.simplify(v - v133))                                             # → 0
  ```
  printed: `0` and `0`.
- **What it means:** the mode is C09's Bessel profile riding on a pressure that the wall's stretch sets. Fails as D18.
- **Traps:** the "constant" contains the wall slope, which depends on A — solve for it, do not treat it as known; principal
  root of i (β, not −β; J₀ is even, J₁ odd, so the sign matters in the J₁ terms).

### D21 · The 2 × 2 system (16.134)–(16.135) and the characteristic quadratic (16.136) — ★★★, 13 steps (curation: 12), in C11 (notebook)
- **Goal:** apply the wall conditions to the mode and find which wavenumbers k a wave of frequency ω may have. · **Start:**
  (16.131), (16.133), the conditions $\hat u(a_0)=-i\omega\hat\xi$, $\hat v(a_0)=-i\omega\hat\eta$ (16.124)–(16.125) and $\hat\xi$ from (16.129). · **Plan:** • two
  conditions → two homogeneous equations in (η̂, A) • write them with X and β̂ • the determinant must vanish • expand into a
  quadratic • pick the root. · **Tools:** non-trivial solutions of a homogeneous 2 × 2 system (P330, reminder); the quadratic
  formula with complex coefficients (P159, reminder); choosing a root by continuation (P386, primer in C11); sympy (P40). ·
  **Assumptions:** as D18–D20.
- **Steps:**
  1. **did** Apply the radial condition · **tex** $-\dfrac{ikA}{\beta}J_1(\hat\beta)-\dfrac{ik^2Eh}{2\omega\rho a_0}\hat\eta+\dfrac{ik\hat\nu}{2\beta}AJ_1(\hat\beta)=-i\omega\hat\eta$ · **why**
     (16.125) with (16.133) at r = a₀, $\hat\beta=\beta a_0$; $\dfrac{a_0}{2a_0^2}=\dfrac1{2a_0}$. · **plain** The blood's radial speed at the wall
     equals the wall's.
  2. **did** Divide by −iω and collect · **tex** $\hat\eta\left[1-\dfrac{k^2Eh}{\omega^22\rho a_0}\right]=AJ_1(\hat\beta)\left[\dfrac{k}{\omega\beta}-\dfrac{k\hat\nu}{2\omega\beta}\right]$ (16.135) ·
     **why** Divide every term by −iω; move the η̂ terms to the left. · **plain** First homogeneous equation in η̂ and A.
  3. **did** Write the axial condition with (16.129) · **tex** $-i\omega\hat\xi=\dfrac{\omega\hat\nu}{ka_0}\hat\eta-\dfrac{i\omega\mu(1-\hat\nu^2)}{Ehk^2}A\beta J_1(\hat\beta)$ ·
     **why** (16.124) needs $-i\omega\hat\xi$; (16.129) times −iω, with the wall slope $-A\beta J_1(\hat\beta)$ (D20 step 5). · **plain** The wall's axial speed.
  4. **did** Equate with û at the wall · **tex** $\hat\eta\left[\dfrac\omega k\dfrac{\hat\nu}{a_0}-\dfrac{kEh}{\omega\rho a_0^2}\right]=A\left[J_0(\hat\beta)+J_1(\hat\beta)\left\{\dfrac{i\beta\omega\mu(1-\hat\nu^2)}{Ehk^2}-\dfrac{\hat\nu}{\hat\beta}\right\}\right]$ (16.134) ·
     **why** (16.131) at r = a₀ equals step 3; η̂ terms to the left, A terms to the right. · **plain** Second homogeneous equation.
  5. **did** Introduce X · **tex** $X\equiv\dfrac{k^2}{\omega^2}\dfrac{Eh}{2\rho a_0}=\left(\dfrac{k}{k_0}\right)^2,\quad k_0=\dfrac\omega{c_0}$ · **why** The combination that
     appears in (16.135); with $c_0^2=\dfrac{Eh}{2\rho a_0}$ (16.62) it is (k/k₀)², so X = 1 is the inviscid wave of C08. · **plain** X compares
     the wave with the Moens–Korteweg wave.
  6. **did** Hide μ inside β̂ · **tex** $\dfrac{i\beta\omega\mu(1-\hat\nu^2)}{Ehk^2}=-\dfrac{1-\hat\nu^2}{2\hat\beta X}$ · **why** μ = ρν = ρiω/β² (step 2 of D20);
     $i\cdot i=-1$; $Ehk^2/(\omega^2\rho)=2a_0X$. · **plain** Viscosity enters only through β̂, that is through α.
  7. **did** Rewrite (16.134) with X and β̂ · **tex** $\hat\eta\,\dfrac{\omega}{ka_0}(\hat\nu-2X)=A\left[J_0-\dfrac{\hat\nu}{\hat\beta}J_1-\dfrac{1-\hat\nu^2}{2\hat\beta X}J_1\right]$ · **why** Factor
     ω/(ka₀) from the left bracket: $\dfrac{kEh}{\omega\rho a_0^2}=\dfrac{\omega}{ka_0}2X$; step 6 on the right; J's at β̂. · **plain** The axial
     condition in the new variables.
  8. **did** Rewrite (16.135) with X and β̂ · **tex** $\hat\eta\,(1-X)=AJ_1\,\dfrac{ka_0}{\omega\hat\beta}\left(1-\dfrac{\hat\nu}2\right)$ · **why**
     $\dfrac{k}{\omega\beta}=\dfrac{ka_0}{\omega\hat\beta}$; factor $1-\hat\nu/2$. · **plain** The radial condition in the new variables.
  9. **did** Require a non-zero solution · **tex** $\dfrac{\omega}{ka_0}(\hat\nu-2X)\cdot J_1\dfrac{ka_0}{\omega\hat\beta}\left(1-\dfrac{\hat\nu}2\right)-(1-X)\left[J_0-\dfrac{\hat\nu}{\hat\beta}J_1-\dfrac{1-\hat\nu^2}{2\hat\beta X}J_1\right]=0$ ·
     **why** Two homogeneous linear equations in η̂ and A have a solution other than zero only if their determinant vanishes
     (P330). · **plain** Only special k allow a wave at all.
  10. **did** Multiply by $2\hat\beta X/J_1$ · **tex** $X(2-\hat\nu)(\hat\nu-2X)-(1-X)\left[2Xg-2\hat\nu X-(1-\hat\nu^2)\right]=0,\quad g\equiv\hat\beta\dfrac{J_0(\hat\beta)}{J_1(\hat\beta)}$ ·
      **why** $\dfrac{\omega}{ka_0}\cdot\dfrac{ka_0}{\omega\hat\beta}=\dfrac1{\hat\beta}$ and $2(1-\hat\nu/2)=2-\hat\nu$; J₁(β̂) ≠ 0 for complex β̂. · **plain** A
      polynomial equation in X.
  11. **did** Expand and collect powers of X · **tex** $X^2\left[2\hat\beta\dfrac{J_0(\hat\beta)}{J_1(\hat\beta)}-4\right]+X\left[4\hat\nu-1-2\hat\beta\dfrac{J_0(\hat\beta)}{J_1(\hat\beta)}\right]+(1-\hat\nu^2)=0$ (16.136) ·
      **why** X²: $-2(2-\hat\nu)+2g-2\hat\nu=2g-4$; X: $\hat\nu(2-\hat\nu)-2g+2\hat\nu-1+\hat\nu^2=4\hat\nu-1-2g$; constant $1-\hat\nu^2$. ·
      **plain** The characteristic quadratic — the book's result, checked.
  12. **did** Check the inviscid limit · **tex** $\alpha\to\infty:\ g\approx\hat\beta\to\infty,\quad X^2-X\approx0\ \Rightarrow\ X=1$ · **why** For large
      |β̂| with Im β̂ > 0, J₀/J₁ → −i, so |g| grows like α and the g terms dominate; dividing by 2g leaves X² − X = 0. · **plain** Without
      viscosity the wave is exactly Moens–Korteweg's.
  13. **did** Choose the root and read speed and damping · **tex** $\dfrac{k}{k_0}=\sqrt X\ (\mathrm{Re}>0),\quad\dfrac{c}{c_0}=\dfrac1{\mathrm{Re}\sqrt X},\quad\dfrac{k_2}{k_1}=\dfrac{\mathrm{Im}\sqrt X}{\mathrm{Re}\sqrt X}$ ·
      **why** Two roots (trap T7): the pressure wave is the one that tends to X = 1 (step 12); follow it down in α by continuation
      (P386). Re > 0: a wave travelling in +x. · **plain** Slower than c₀ and fading as it goes.
- **Result:** (16.136); $c/c_0=1/\mathrm{Re}\sqrt X$, $k_2/k_1=\mathrm{Im}\sqrt X/\mathrm{Re}\sqrt X$.
- **Check:** **expect** (ν̂ = 0.5): α = 8 → c/c₀ = 0.954, k₂/k₁ = 0.0622 (the amplitude falls to 0.68 in one wavelength); α = 20 →
  0.981, 0.0218; α = 3 → 0.910, 0.252; the other root at α = 8 has c/c₀ = 4.93 (fast, heavily damped — excluded); 1 − c/c₀ ∝ 1/α at
  large α; the independent ODE route `elastic_tube_determinant_numeric` vanishes at the closed-form root to 10⁻⁸.
- **sympy check (`check_src`, executed):**
  ```python
  import sympy as sp                                                       # symbolic algebra
  r, w, rho, k, E, h, a0, nu, nuk, A, eta, X, bh = sp.symbols('r omega rho k E h a0 nuhat nu_k A eta X betahat')
  beta, mu = sp.sqrt(sp.I*w/nuk), rho*nuk                                  # β = √(iω/ν), μ = ρν
  u = A*sp.besselj(0, beta*r) + k/w*E*h/(rho*a0**2)*eta - nu/(beta*a0)*A*sp.besselj(1, beta*a0)   # (16.131)
  s = sp.symbols('s', positive=True)
  v = -sp.I*k/r*sp.integrate(s*u.subs(r, s), (s, 0, r))                   # (16.132) → (16.133)
  slope = sp.diff(u, r).subs(r, a0)                                        # dû/dr at the wall
  xi = sp.I*nu/(k*a0)*eta - (1 - nu**2)/(E*h*k**2)*mu*slope               # (16.129)
  c124 = sp.expand(u.subs(r, a0) + sp.I*w*xi)                              # (16.124): û(a0) + iω ξ̂ = 0
  c125 = sp.expand(v.subs(r, a0) + sp.I*w*eta)                             # (16.125): v̂(a0) + iω η̂ = 0
  M = sp.Matrix([[c124.diff(eta), c124.diff(A)], [c125.diff(eta), c125.diff(A)]])   # coefficients of (η̂, A)
  det = M.det()                                                            # step 9: must vanish
  g = bh*sp.besselj(0, bh)/sp.besselj(1, bh)                               # β̂ J0/J1
  quad = X**2*(2*g - 4) + X*(4*nu - 1 - 2*g) + (1 - nu**2)                 # (16.136)
  vals = {w: 1.3, nuk: 0.07, rho: 1.1, E: 2.0, h: 0.05, a0: 0.4, nu: 0.45} # arbitrary test numbers
  b = complex(sp.N(sp.sqrt(sp.I*vals[w]/vals[nuk])*vals[a0]))              # β̂ for them
  for kv in (0.7+0.1j, 1.3-0.2j, 2.0+0.5j):                                # three complex wavenumbers
      Xv = kv**2*vals[E]*vals[h]/(2*vals[rho]*vals[a0]*vals[w]**2)         # step 5: X
      ratio = complex(sp.N(det.subs(vals).subs(k, kv)))*2*b*Xv/(complex(sp.N(sp.besselj(1, b)))*complex(sp.N(quad.subs({X: Xv, nu: 0.45, bh: b}))))
      print(ratio)                                                         # the same constant iω every time
  ```
  printed: `1.3j` three times (to 10⁻¹⁵) — the determinant is $i\omega J_1(\hat\beta)/(2\hat\beta X)$ times (16.136), so it vanishes exactly
  where (16.136) does.
- **What it means:** viscosity slows the pulse and makes it fade with distance, by amounts that shrink like 1/α — small in the
  aorta, large in small arteries. Fails for short waves, thick walls and non-Newtonian blood.
- **Traps:** T7 (two roots); μ is hidden in β̂; Re √X > 0 picks the wave travelling in +x in this convention (T1).

### D22 · Murray's optimal radius (16.139)–(16.140, corrected) and the minimum cost (16.141) — ★★, 8 steps, in C12 (notebook · `murray_branching`)
- **Goal:** find the radius that makes a blood vessel cheapest to run, given the flow it must carry. · **Start:** the two costs
  of a segment — pumping blood through it and keeping its blood alive. · **Plan:** • price the pumping • price the upkeep • add
  • set the slope to zero • solve • check it is a minimum • evaluate. · **Tools:** derivative = 0 at a minimum and the second
  derivative test; the minimum of a falling plus a rising power (P349, reminder). · **Assumptions:** Poiseuille flow in the
  segment (step 1); upkeep proportional to blood volume (step 2).
- **Steps:**
  1. **did** Price the pumping · **tex** $P_{pump}=\Delta p\,Q=RQ^2=\dfrac{8\mu L}{\pi a^4}Q^2$ · **why** Power = pressure drop × volume flow; with
     $R=\dfrac{8\mu L}{\pi a^4}$ (16.20) the drop is RQ. · **plain** Narrow vessels are expensive to push blood through.
  2. **did** Price the upkeep · **tex** $P_{up}=K\,\pi a^2L$ · **why** The body pays a power K per unit volume of blood held (W/m³ — the
     book's K, not the compliance, trap T2); the segment holds πa²L. · **plain** Wide vessels are expensive to fill and maintain.
  3. **did** Add the two · **tex** $\hat E_1=\dfrac{8\mu L}{\pi a^4}Q^2+K\pi a^2L$ (16.138) · **why** Total power to run the segment. · **plain** One
     falling and one rising term.
  4. **did** Set the slope to zero · **tex** $\dfrac{\partial\hat E_1}{\partial a}=-\dfrac{32\mu L}{\pi}Q^2a^{-5}+2K\pi La=0$ (16.139) · **why** At a minimum of a
     smooth function the derivative vanishes; power rule: $d(a^{-4})/da=-4a^{-5}$, $d(a^2)/da=2a$. · **plain** At the best radius
     widening saves as much pumping as it adds upkeep.
  5. **did** Solve for $a^6$ · **tex** $a^6=\dfrac{16\mu Q^2}{\pi^2K}$ · **why** Multiply (16.139) by $a^5/(2\pi KL)$: $\dfrac{32\mu Q^2}{\pi\cdot2\pi K}=\dfrac{16\mu Q^2}{\pi^2K}$. ·
     **plain** One power of a, isolated.
  6. **did** Take the sixth root · **tex** $a=\left[\dfrac{16\mu}{\pi^2K}\right]^{1/6}Q^{1/3}$ (16.140, corrected) · **why** $(Q^2)^{1/6}=Q^{1/3}$. The book prints
     $\left[\dfrac{16\mu}{\pi^2}K\right]^{1/6}$ — K upside down (slip S15); units decide: μQ²/K = (Pa s)(m⁶/s²)/(Pa/s) = m⁶ ✓. · **plain** Murray's law:
     the flow grows as the cube of the radius.
  7. **did** Check it is a minimum · **tex** $\dfrac{\partial^2\hat E_1}{\partial a^2}=\dfrac{160\mu L}{\pi}Q^2a^{-6}+2K\pi L>0$ · **why** Both terms are positive for
     every a > 0, so the stationary point is the minimum. · **plain** Any other radius costs more.
  8. **did** Evaluate the cost at the optimum · **tex** $\hat E_{1,\min}=\dfrac{3\pi}2KLa^2$ (16.141) · **why** Step 5 gives $\dfrac{8\mu LQ^2}{\pi a^4}=\dfrac{\pi}{2}KLa^2$:
     the pumping is half the upkeep; add $\pi KLa^2$. · **plain** At the optimum one third of the cost is pumping, two thirds upkeep.
- **Result:** $a=\left[\dfrac{16\mu}{\pi^2K}\right]^{1/6}Q^{1/3}$ (16.140, corrected), $\hat E_{1,\min}=\dfrac{3\pi}2KLa^2$ (16.141).
- **Check:** units ✓ (step 6); **expect** Q = 82.5 mL/s, μ = 3.5 mPa s, K = 1000 W/m³: a = **5.81 mm**; the printed (16.140) gives
  58.1 mm (ten times too big here, and the wrong K-dependence); pumping 5.31 mW and upkeep 10.62 mW per 10 cm (ratio ½ ✓); a grid
  search over a finds the same optimum.
- **What it means:** nature sizes vessels so that flow ∝ radius³ — the most testable claim of the chapter. Fails where the
  flow is pulsatile or non-Newtonian, or the cost is not proportional to volume.
- **Traps:** S15; K here is a cost per volume (T2); the ½ ratio is a property of the optimum, not of any radius.

### D23 · The cube law (16.144) — ★, 3 steps, in C12 (notebook · `murray_branching`)
- **Goal:** find how the radii at a fork must be related if every segment is optimal. · **Start:**
  $a_i=\left[\dfrac{16\mu}{\pi^2K}\right]^{1/6}Q_i^{1/3}$ (16.143, with the corrected bracket) for i = 0, 1, 2. · **Plan:** • invert for Q_i • conserve flow • cancel. · **Tools:**
  mass conservation at a node. · **Assumptions:** the same μ and K in every segment (step 1).
- **Steps:**
  1. **did** Invert the optimum for each segment · **tex** $Q_i=\left[\dfrac{\pi^2K}{16\mu}\right]^{1/2}a_i^3,\quad i=0,1,2$ · **why** Cube both sides of (16.143)
     and solve for Q_i; the constant is the same for all three. · **plain** Each segment's flow is a constant times its radius cubed.
  2. **did** Conserve flow at the fork · **tex** $Q_0=Q_1+Q_2$ · **why** Incompressible blood; nothing stored at the junction. · **plain** What
     enters the fork leaves it.
  3. **did** Substitute and cancel the constant · **tex** $a_0^3=a_1^3+a_2^3$ (16.144) · **why** Step 1 in step 2; divide by the common
     constant. · **plain** Murray's cube law.
- **Result:** (16.144).
- **Check:** symmetric fork: $a_1=a_2=2^{-1/3}a_0=0.794a_0$; the daughters' total area is $2\times0.630=1.26$ times the parent's, so
  the mean speed falls by that factor.
- **What it means:** the law needs only the radii — a test you can make on a cast of a vascular tree. Fails as D22.
- **Traps:** cubes, not squares (area conservation would be a₀² = a₁² + a₂²).

### D24 · The branching angles (16.146)–(16.147) (rule b: the book only says "making displacements") — ★★★, 11 steps, in C12 (notebook · `murray_branching`)
- **Goal:** find where the branch point B of a fork should sit, and so the angles of the daughters, when A, C and D are fixed. ·
  **Start:** $\delta\hat E_{2,\min}=\dfrac{3\pi}2K\left(\delta L_0\,a_0^2+\delta L_1\,a_1^2+\delta L_2\,a_2^2\right)$ (16.145) — "moving B changes only the lengths". · **Plan:** • the cost as a
  weighted length • how each length changes when B moves • set the change to zero for every direction • a triangle of
  "forces" • the cosine rule. · **Tools:** the directional derivative of a length and the weighted Fermat point (P388, primer in
  C12, built on P75); the cosine rule (P340, reminder); `scipy.optimize.minimize` (P387, primer in C12) for the numeric check. ·
  **Assumptions:** a planar fork; radii fixed by the flows (D23); A, C, D fixed.
- **Steps:**
  1. **did** Write the fork's cost as a weighted length · **tex** $\hat E_{2,\min}=\dfrac{3\pi}2K\left(a_0^2L_0+a_1^2L_1+a_2^2L_2\right)$ (16.142) · **why** (16.141)
     for each optimal segment, added; only the lengths L_i depend on where B is. · **plain** Each segment costs its length times
     its radius squared.
  2. **did** Name the unit vectors of the segments · **tex** $\mathbf e_0=\dfrac{B-A}{L_0},\quad\mathbf e_1=\dfrac{C-B}{L_1},\quad\mathbf e_2=\dfrac{D-B}{L_2}$ · **why**
     Directions of flow along each segment; they set the angles we want. · **plain** Three arrows along the parent and the two
     daughters.
  3. **did** Move B a small step ds along a direction m · **tex** $\delta L_0=\mathbf e_0\cdot\mathbf m\,ds,\quad\delta L_1=-\mathbf e_1\cdot\mathbf m\,ds,\quad\delta L_2=-\mathbf e_2\cdot\mathbf m\,ds$ ·
     **why** P388: a distance from a fixed point grows at the rate (unit vector away from that point)·m; e₀ points away from A,
     while e₁, e₂ point toward C and D. · **plain** Moving B toward C shortens L₁ and lengthens L₀.
  4. **did** Put them into (16.145) · **tex** $\delta\hat E_{2,\min}=\dfrac{3\pi}2K\left(a_0^2\mathbf e_0-a_1^2\mathbf e_1-a_2^2\mathbf e_2\right)\cdot\mathbf m\,ds$ · **why** Each δL_i of step 3
     times its weight $a_i^2$. · **plain** The cost change is a dot product with the direction of the move.
  5. **did** Demand no change in any direction · **tex** $a_0^2\mathbf e_0=a_1^2\mathbf e_1+a_2^2\mathbf e_2$ · **why** At the optimum δÊ = 0 for every m (the book's
     "displacements one at a time"); a vector whose dot product with every m vanishes is zero. · **plain** Three pulls of sizes
     a₀², a₁², a₂² balance: they close a triangle.
  6. **did** Define the angles · **tex** $\cos\theta=\mathbf e_0\cdot\mathbf e_1,\quad\cos\phi=\mathbf e_0\cdot\mathbf e_2$ · **why** θ and φ are measured from the
     extension of the parent, as the book's Fig. 16.15 does (not from the parent itself). · **plain** How far each daughter turns.
  7. **did** Isolate one daughter and square · **tex** $a_2^4=\left|a_0^2\mathbf e_0-a_1^2\mathbf e_1\right|^2$ · **why** Step 5 gives $a_2^2\mathbf e_2=a_0^2\mathbf e_0-a_1^2\mathbf e_1$;
     |e₂| = 1. · **plain** One side of the triangle from the other two.
  8. **did** Expand the square · **tex** $a_2^4=a_0^4+a_1^4-2a_0^2a_1^2\cos\theta$ · **why** $|\mathbf u-\mathbf v|^2=|\mathbf u|^2+|\mathbf v|^2-2\mathbf u\cdot\mathbf v$ —
     the cosine rule (P340). · **plain** The triangle's sides fix its angles.
  9. **did** Solve for the three cosines · **tex** $\cos\theta=\dfrac{a_0^4+a_1^4-a_2^4}{2a_0^2a_1^2},\quad\cos\phi=\dfrac{a_0^4-a_1^4+a_2^4}{2a_0^2a_2^2},\quad\cos(\theta+\phi)=\dfrac{a_0^4-a_1^4-a_2^4}{2a_1^2a_2^2}$ (16.146) ·
     **why** The first from step 8; the second by isolating $a_1^2\mathbf e_1$ instead; the third by squaring step 5 itself, with
     $\mathbf e_1\cdot\mathbf e_2=\cos(\theta+\phi)$ (the daughters leave on opposite sides). · **plain** The angles depend only on the radii.
  10. **did** Eliminate a₂ with the cube law · **tex** $\cos\theta=\dfrac{a_0^4+a_1^4-(a_0^3-a_1^3)^{4/3}}{2a_0^2a_1^2}$ (16.147) · **why** (16.144): $a_2^3=a_0^3-a_1^3$, so
      $a_2^4=(a_0^3-a_1^3)^{4/3}$. · **plain** One daughter's radius sets its angle.
  11. **did** Take the symmetric fork · **tex** $a_1=a_2=2^{-1/3}a_0:\quad\cos\theta=2^{-1/3},\ \ \theta=37.47^\circ$ (ours) · **why** $a_1^4=a_2^4$ cancel;
      $\dfrac{a_0^4}{2a_0^2a_1^2}=\dfrac{2^{2/3}}2=2^{-1/3}$. · **plain** Each daughter leaves at 37.5°: 75° between them.
- **Result:** (16.146), (16.147); symmetric θ = arccos 2^{−1/3} = 37.47° (ours).
- **Check:** a₂ → 0: a₁ → a₀, cos θ → 1, θ → 0 — a tiny side branch does not bend the main vessel ✓; a 70/30 split: a₁ = 0.888a₀,
  a₂ = 0.669a₀, θ = **25.7°**, φ = **49.7°** (expect) — the thinner daughter turns more; the numeric minimisation of the weighted
  length reproduces θ to 10⁻⁶ rad.
- **sympy check (`check_src`, executed):**
  ```python
  import sympy as sp, numpy as np                                          # symbols, then a numeric optimisation
  from scipy.optimize import minimize                                      # P387
  a0, a1, a2, th = sp.symbols('a0 a1 a2 theta', positive=True)
  lhs = a2**4                                                              # step 7: |a2² e2|²
  rhs = (a0**2 - a1**2*sp.cos(th))**2 + (a1**2*sp.sin(th))**2              # |a0² e0 − a1² e1|² with e0 = (1, 0)
  cos_th = sp.solve(sp.Eq(lhs, sp.expand(sp.trigsimp(rhs))), sp.cos(th))[0]   # steps 8–9
  print(sp.simplify(cos_th - (a0**4 + a1**4 - a2**4)/(2*a0**2*a1**2)))     # (16.146) → 0
  A, C, D = np.array([0.0, 0.0]), np.array([3.0, 1.2]), np.array([2.6, -1.5])   # fixed ends (arbitrary)
  r0, r1 = 1.0, 0.7**(1/3); r2 = (1.0 - 0.7)**(1/3)                        # a 70/30 split by the cube law
  cost = lambda B: r0**2*np.linalg.norm(B - A) + r1**2*np.linalg.norm(C - B) + r2**2*np.linalg.norm(D - B)   # (16.142)/(3πK/2)
  B = minimize(cost, [1.5, 0.0], method="Nelder-Mead", options=dict(xatol=1e-12, fatol=1e-14, maxiter=20000)).x
  e0, e1 = (B - A)/np.linalg.norm(B - A), (C - B)/np.linalg.norm(C - B)   # step 2
  print(np.degrees(np.arccos(e0 @ e1)), np.degrees(np.arccos((r0**4 + r1**4 - r2**4)/(2*r0**2*r1**2))))   # equal
  ```
  printed: `0` and `25.705664149955744 25.70566416239578` (the numeric optimum and (16.146) agree to 10⁻⁸ degree).
- **What it means:** branch angles are not free: the optimum is a balance of three "pulls" proportional to radius squared, so
  thick daughters stay close to the parent's line. Fails if the ends are not fixed or the fork is not planar.
- **Traps:** the weights are a², not a; angles from the parent's extension.

### D25 · Reflection and transmission at a junction (16.153)–(16.154) — ★, 6 steps, in C13 (notebook · `junction_reflection`)
- **Goal:** find how much of a pressure pulse bounces back from a fork and how much goes on. · **Start:** $p_I+p_R=p_{T1}=p_{T2}$ (16.148)
  and $Q=Au=\pm Yp=\pm\dfrac pZ$ (16.151). · **Plan:** • each wave's flux with its sign • flux conservation • common pressure • solve one linear
  equation. · **Tools:** the electrical analogy and admittance (P381); the water-hammer relation $\hat p_1=\rho c\hat u_1$, $\hat p_2=-\rho c\hat u_2$ (16.68, R10). ·
  **Assumptions:** inviscid long waves; no storage at the junction (step 2).
- **Steps:**
  1. **did** Write each wave's flux with its direction · **tex** $Q_I=Y_1p_I,\quad Q_R=-Y_1p_R,\quad Q_{Tj}=Y_jp_{Tj}$ · **why** (16.151): a wave moving
     in +x carries Q = +Yp, the reflected wave moving back carries Q = −Yp — the sign of (16.68). We keep Q_R signed and add; the
     book's (16.149) subtracts as well and so counts the sign twice (slip S16). · **plain** A pressure pulse pushes blood the way
     it travels.
  2. **did** Conserve flux at the junction · **tex** $Y_1(p_I-p_R)=\sum_jY_jp_{Tj}$ (16.152) · **why** Nothing is stored at the node: the
     parent's flux equals the daughters' (step 1 substituted). · **plain** Blood in = blood out at the fork.
  3. **did** Use the common pressure · **tex** $p_{Tj}=p_I+p_R$ · **why** (16.148): the pressure has one value at the junction. · **plain**
     Both daughters start with the parent's total pressure.
  4. **did** Substitute and collect $p_R$ · **tex** $p_R\left(Y_1+\sum Y_j\right)=p_I\left(Y_1-\sum Y_j\right)$ · **why** $Y_1p_I-Y_1p_R=\sum Y_j(p_I+p_R)$;
     move terms. · **plain** One equation for the reflected pressure.
  5. **did** Divide · **tex** $\mathrm R=\dfrac{p_R}{p_I}=\dfrac{Y_1-\sum Y_j}{Y_1+\sum Y_j}$ (16.153) · **why** $Y_1+\sum Y_j>0$. · **plain** Reflection is set by the
     admittance mismatch.
  6. **did** Get the transmission from step 3 · **tex** $\mathrm T=\dfrac{p_{Tj}}{p_I}=1+\mathrm R=\dfrac{2Y_1}{Y_1+\sum Y_j}$ (16.153) · **why** $p_{Tj}=p_I+p_R$; common
     denominator. With Z = 1/Y the same reads $\dfrac{p_R}{p_I}=\dfrac{Z_1^{-1}-\sum Z_j^{-1}}{Z_1^{-1}+\sum Z_j^{-1}}$ (16.154). · **plain** The pressure that goes on.
- **Result:** (16.153), (16.154).
- **Check:** ΣY_j = Y₁ → R = 0, T = 1 (matched); ΣY_j → 0 → R → +1, T → 2 (closed end: the pressure doubles, trap T18); ΣY_j → ∞ →
  R → −1 (open end); **expect** two daughters of radius 0.8a with equal wave speed: ΣY/Y₁ = 1.28, R = **−0.123**, T = 0.877; a
  Murray-optimal symmetric fork: R = −0.115.
- **What it means:** a fork reflects only if its daughters together are "stiffer" or "softer" than the parent; Murray's forks
  are close to matched. Fails when the junction stores blood or the waves are short.
- **Traps:** S16; T = 1 + R, not 1 − R.

### D26 · Energy at a junction (16.161)–(16.163) and $\mathrm R^2+\dfrac{\sum Y_j}{Y_1}\mathrm T^2=1$ (ours) — ★★, 6 steps, in C13 (notebook · `junction_reflection`)
- **Goal:** check that the reflected and transmitted waves carry exactly the energy of the incident one. · **Start:** a
  travelling wave's power = pressure × flux, with Q = Yp. · **Plan:** • power of each wave • the shares • add them • a perfect
  square. · **Tools:** the acoustic energy flux of Ch. 15 (recap); (16.153). · **Assumptions:** as D25; equal daughter pressures (step 4).
- **Steps:**
  1. **did** Write a wave's power · **tex** $\text{power}=pQ=Yp^2=\dfrac{p^2}{Z}$ · **why** Power = pressure × volume flux (the acoustic
     intensity p′u′ of Ch. 15 times the area); Q = Yp for a travelling wave. · **plain** A pulse carries energy in proportion to p².
  2. **did** Incident power · **tex** $\text{incident}=\dfrac{p_I^2}{Z_1}$ (16.161) · **why** Step 1 in the parent. · **plain** What arrives.
  3. **did** Reflected power · **tex** $\text{reflected}=\dfrac{(\mathrm Rp_I)^2}{Z_1}=\mathrm R^2\dfrac{p_I^2}{Z_1}$ (16.162) · **why** The reflected wave is in the parent,
     with $p_R=\mathrm Rp_I$; we count its power as a magnitude. · **plain** R² is the energy reflection coefficient.
  4. **did** Transmitted share · **tex** $\dfrac{p_{T2}^2/Z_2+p_{T3}^2/Z_3}{p_I^2/Z_1}=\dfrac{Z_2^{-1}+Z_3^{-1}}{Z_1^{-1}}\mathrm T^2$ (16.163) · **why** $p_{T2}=p_{T3}=\mathrm Tp_I$ by (16.148); factor T². ·
     **plain** Not T² alone: the daughters' admittance matters.
  5. **did** Add the two shares · **tex** $\mathrm R^2+\dfrac{S}{Y_1}\mathrm T^2=\dfrac{(Y_1-S)^2+4SY_1}{(Y_1+S)^2},\quad S\equiv\sum Y_j$ · **why** Square (16.153):
     $\dfrac{S}{Y_1}\dfrac{4Y_1^2}{(Y_1+S)^2}=\dfrac{4SY_1}{(Y_1+S)^2}$; common denominator. · **plain** Everything over one fraction.
  6. **did** Recognise a perfect square · **tex** $\mathrm R^2+\dfrac{\sum Y_j}{Y_1}\mathrm T^2=1$ (ours) · **why** $(Y_1-S)^2+4SY_1=Y_1^2+2SY_1+S^2=(Y_1+S)^2$. ·
     **plain** The junction neither stores nor loses energy.
- **Result:** $\mathrm R^2+\dfrac{\sum Y_j}{Y_1}\mathrm T^2=1$ (ours).
- **Check:** radius ratio 0.8: 0.0151 + 1.28 × 0.7695 = 0.0151 + 0.9849 = **1** ✓; closed end (S = 0): R² = 1 — all reflected.
- **What it means:** a pressure wave can grow on reflection (T up to 2) without creating energy, because the energy is pressure
  × flow and the flow shrinks. Fails with losses at the junction (viscosity, turbulence).
- **Traps:** the transmitted share is not T²; equal pressures in both daughters.

### D27 · Dean's leading orders and the flux coefficient 1541/50400 of (16.178) (rule b: the book gives coefficients only; ours) — ★★★, 14 steps, in C14 (notebook · `dean_curved_pipe`)
- **Goal:** solve the curved-pipe equations order by order for a slow flow, find the secondary vortices and the skewed axial
  flow, and reproduce the book's flux coefficient 0.0306. · **Start:** $\nabla_1^2w-\dfrac\kappa2\dfrac{\partial p_0}{\partial s}=\dfrac\kappa{2r}\left(\dfrac{\partial\psi}{\partial\alpha}\dfrac{\partial w}{\partial r}-\dfrac{\partial\psi}{\partial r}\dfrac{\partial w}{\partial\alpha}\right)$ (16.173),
  $\dfrac2\kappa\nabla_1^4\psi-\dfrac1r\left(\dfrac{\partial\psi}{\partial r}\dfrac{\partial}{\partial\alpha}-\dfrac{\partial\psi}{\partial\alpha}\dfrac{\partial}{\partial r}\right)\nabla_1^2\psi=-2w\left(\sin\alpha\dfrac{\partial w}{\partial r}+\dfrac{\cos\alpha}r\dfrac{\partial w}{\partial\alpha}\right)$ (16.174), $\psi=\dfrac{\partial\psi}{\partial r}=w=0$ at $r=1$ (16.176) and
  $w=\sum\kappa^{2n}w_n,\ \psi=\kappa\sum\kappa^{2n}\psi_n$ (16.177). · **Plan:** • order κ⁰: Poiseuille • leading ψ: the vortex pair • order κ²: the skewed axial correction •
  order κ⁴: the mean change and the flux. · **Tools:** a regular perturbation series (P391, primer in C14); toroidal coordinates
  (P389, primer in C14); the biharmonic operator in polar coordinates and polynomial solutions (P390, primer in C14); sympy (P40);
  integration by parts. · **Assumptions:** δ ≪ 1 (the equations), κ small (the series); fixed pressure gradient.
- **Steps:**
  1. **did** Insert the series; keep order κ⁰ of (16.173) · **tex** $\nabla_1^2w_0=\dfrac\kappa2\dfrac{\partial p_0}{\partial s}\equiv-8$ · **why** P391: the right side of
     (16.173) starts at κ·κ; the gradient term is a constant ($p_0=-Gs$), scaled so the straight-pipe mean speed is 1. · **plain**
     At lowest order the bend does nothing.
  2. **did** Solve for the straight-pipe flow · **tex** $w_0=2(1-r^2)$ · **why** Axisymmetric: $\frac1r(rw_0')'=-8$ gives $w_0=-2r^2+C$; w = 0 at r = 1
     (16.176) gives C = 2; its mean $\int_0^12(1-r^2)2r\,dr=1$. · **plain** Poiseuille, mean speed 1, maximum 2.
  3. **did** Keep the leading order of (16.174) · **tex** $2\nabla_1^4\psi_0=-2w_0\left(\sin\alpha\dfrac{\partial w_0}{\partial r}+\dfrac{\cos\alpha}r\dfrac{\partial w_0}{\partial\alpha}\right)$ · **why**
     $(2/\kappa)\nabla^4(\kappa\psi_0)=2\nabla^4\psi_0$; the nonlinear ψψ term is O(κ²), smaller. · **plain** The axial flow's centrifugal push drives
     the secondary flow.
  4. **did** Evaluate the forcing · **tex** $\nabla_1^4\psi_0=8r(1-r^2)\sin\alpha$ · **why** $\partial w_0/\partial r=-4r$, $\partial w_0/\partial\alpha=0$:
     $-2\cdot2(1-r^2)\sin\alpha\cdot(-4r)=16r(1-r^2)\sin\alpha$; divide by 2. · **plain** Strongest where the fast fluid is flung hardest.
  5. **did** Separate the angle · **tex** $\psi_0=f(r)\sin\alpha,\quad\nabla_1^2(f\sin\alpha)=\left(f''+\dfrac{f'}r-\dfrac f{r^2}\right)\sin\alpha$ · **why** The forcing is ∝ sin α
     and $\partial^2_\alpha\sin\alpha=-\sin\alpha$ (P390). · **plain** One radial function to find.
  6. **did** Try an odd polynomial · **tex** $f=c_1r+c_3r^3+c_5r^5+c_7r^7$ · **why** P390: $\nabla_1^2(r^n\sin\alpha)=(n^2-1)r^{n-2}\sin\alpha$, so ∇⁴ lowers the
     power by 4; the forcing has r and r³, so powers up to 7; odd powers keep ψ₀ smooth on the axis. · **plain** Undetermined
     coefficients.
  7. **did** Match the forcing and the wall · **tex** $c_5=\tfrac1{24},\ c_7=-\tfrac1{144};\quad c_1=\tfrac1{36},\ c_3=-\tfrac1{16}$ · **why**
     $\nabla^4r^5\to192r$, $\nabla^4r^7\to1152r^3$: 192c₅ = 8, 1152c₇ = −8; r and r³ are killed by ∇⁴ and fixed by $f(1)=f'(1)=0$
     (16.176). · **plain** Four numbers, four conditions.
  8. **did** Factor the result · **tex** $\psi_0=\dfrac{r}{144}(1-r^2)^2(4-r^2)\sin\alpha$ (ours) · **why** $144f=4r-9r^3+6r^5-r^7=r(1-r^2)^2(4-r^2)$ (multiply out to check). ·
     **plain** Two counter-rotating cells, one on each side of the plane of the bend.
  9. **did** Keep order κ² of (16.173) · **tex** $\nabla_1^2w_1=\dfrac1{2r}\left(\dfrac{\partial\psi_0}{\partial\alpha}\dfrac{\partial w_0}{\partial r}-\dfrac{\partial\psi_0}{\partial r}\dfrac{\partial w_0}{\partial\alpha}\right)$ · **why** The right
     side of (16.173) with ψ = κψ₀ and w = w₀ is κ² times this bracket. · **plain** The secondary flow carries fast fluid across the
     section.
  10. **did** Evaluate the right side · **tex** $\nabla_1^2w_1=-2f(r)\cos\alpha$ · **why** $\partial\psi_0/\partial\alpha=f\cos\alpha$, $\partial w_0/\partial r=-4r$,
      $\partial w_0/\partial\alpha=0$. · **plain** Negative where cos α > 0.
  11. **did** Solve with $w_1(1)=0$ · **tex** $w_1=+\dfrac{r(1-r^2)(19-21r^2+9r^4-r^6)}{5760}\cos\alpha$ (ours) · **why** $w_1=g(r)\cos\alpha$ with
      $g''+g'/r-g/r^2=-2f$, an odd polynomial to r⁹ (sympy in the cell). Positive where cos α > 0 — **plus sign** (header correction).
      · **plain** The fast core is pushed toward α = 0, the outer wall.
  12. **did** Find the part of $w_2$ that carries flux · **tex** $\dfrac1r\dfrac{d}{dr}\left(r\dfrac{d\overline{w_2}}{dr}\right)=\dfrac{(fg)'}{4r}$ · **why** Only the α-average
      $\overline{w_2}$ changes the flux. Averaging (16.173)'s right side at κ⁴ over α: cos²α and sin²α give ½, the ψ₁ terms (∝ cos α,
      cos 2α) give 0: $\frac1{2r}\cdot\frac12(fg'+f'g)$. · **plain** The mean correction, from ψ₀ and w₁ alone.
  13. **did** Integrate over the section · **tex** $\displaystyle\int_0^1\overline{w_2}\,2r\,dr=-\dfrac14\int_0^1rfg\,dr=-\dfrac14\cdot\dfrac{1541}{1045094400}$ · **why** One
      integration gives $r\overline{w_2}'=fg/4$ (regular on the axis); then by parts $\int\overline{w_2}2r\,dr=[r^2\overline{w_2}]_0^1-\int r^2\overline{w_2}'dr$, with
      $\overline{w_2}(1)=0$; sympy gives the exact integral. · **plain** The bend lowers the flux at order κ⁴.
  14. **did** Convert to Dean's K · **tex** $\dfrac{Q}{\pi a^2\bar W}=1-\dfrac{1541}{50400}\left(\dfrac K{576}\right)^2$, $\dfrac{1541}{50400}=0.030575$ (16.178) · **why** The flux ratio is
      $1+\kappa^4\int\overline{w_2}2r\,dr$; K = 2κ² so κ⁴ = K²/4; $\dfrac{576^2\cdot1541}{16\cdot1045094400}=\dfrac{1541}{50400}$. · **plain** The book's 0.0306,
      reproduced exactly.
- **Result:** ψ₀ and w₁ (ours) and $\dfrac{Q}{\pi a^2\bar W}=1-0.0306\left(\dfrac K{576}\right)^2+0.0120\left(\dfrac K{576}\right)^4+O(K^6)$ (16.178) with $0.0306=1541/50400$.
- **Check:** ψ₀ = ∂ψ₀/∂r = 0 and w₁ = 0 at r = 1 ✓; **expect** ψ₀ is largest (0.00757) at r = 0.429; at K = 576 the axial maximum
  sits 0.187a toward the outer wall and Q/(πa²W̄) = 0.9814; the resistance ratio $\lambda=1+0.0306\left(\dfrac K{576}\right)^2-0.0110\left(\dfrac K{576}\right)^4$ (16.182) is its reciprocal
  (coefficient $0.0306^2-0.0120=-0.0111$).
- **sympy check (`check_src`, executed — re-runs steps 1–14):**
  ```python
  import sympy as sp                                                       # symbolic algebra
  r, al = sp.symbols('r alpha', positive=True)                             # radius (a = 1) and the toroidal angle
  lap = lambda F: sp.diff(F, r, 2) + sp.diff(F, r)/r + sp.diff(F, al, 2)/r**2   # ∇₁² of (16.175)
  w0 = 2*(1 - r**2); print(lap(w0))                                        # steps 1–2: −8
  c = sp.symbols('c1 c3 c5 c7'); f = c[0]*r + c[1]*r**3 + c[2]*r**5 + c[3]*r**7   # step 6
  rhs = -2*w0*sp.sin(al)*sp.diff(w0, r)                                    # steps 3–4: right side of 2∇⁴ψ₀
  eqs = sp.Poly(sp.expand((2*lap(lap(f*sp.sin(al))) - rhs)/sp.sin(al)), r).coeffs()
  f = sp.factor(f.subs(sp.solve(list(eqs) + [f.subs(r, 1), sp.diff(f, r).subs(r, 1)], c)))   # step 7
  print(f)                                                                 # step 8: r(1 − r²)²(4 − r²)/144
  d = sp.symbols('d1 d3 d5 d7 d9'); g = sum(di*r**(2*i + 1) for i, di in enumerate(d))   # step 11 ansatz
  rhs1 = sp.diff(f*sp.sin(al), al)*sp.diff(w0, r)/(2*r)                    # steps 9–10
  eqs = sp.Poly(sp.simplify(sp.expand((lap(g*sp.cos(al)) - rhs1)/sp.cos(al))), r).coeffs()
  g = sp.factor(g.subs(sp.solve(list(eqs) + [g.subs(r, 1)], d)))           # w1 = g cos α
  print(g)                                                                 # +r(1 − r²)(19 − 21r² + 9r⁴ − r⁶)/5760
  I = sp.integrate(r*f*g, (r, 0, 1)); print(I)                             # step 13: 1541/1045094400
  K = sp.symbols('K', positive=True)
  print(sp.nsimplify(((K**2/4)*I/4)/(K/576)**2))                           # step 14: 1541/50400
  ```
  printed: `-8`; `-r*(r - 2)*(r - 1)**2*(r + 1)**2*(r + 2)/144`; `r*(r - 1)*(r + 1)*(r**6 - 9*r**4 + 21*r**2 - 19)/5760` (= the plus
  sign); `1541/1045094400`; `1541/50400` (1.2 s).
- **What it means:** curvature flings the fast core outward; it returns along the walls as two vortices; the flux falls only at
  order K², so gentle bends cost little — but the secondary flow mixes the blood and raises the wall shear on the outer side.
  Fails beyond small K (the series diverges; the reliability criterion of E8 is ours).
- **Traps:** the κ, K, D conventions (T8); ψ = ∂ψ/∂r = 0 at the wall; sin α forcing for ψ₀, cos α for w₁; α is an angle here (T2).

### D28 · The wave speed of a collapsible tube (16.189) — ★, 4 steps, in C15 (notebook · `collapsible_tube_choking`)
- **Goal:** find how fast a small wave runs along a partly collapsed tube. · **Start:** $c^2=\dfrac A\rho\dfrac{d(p-p_e)}{dA}$ (16.185) and
  $-\mathrm P\approx\alpha^{-n}-1$ (16.188) with $\mathrm P=\dfrac{p-p_e}{K_p}$, $\alpha=\dfrac A{A_0}$ (16.186). · **Plan:** • write p − p_e with P • change variable to α • differentiate the power law • substitute. · **Tools:**
  chain rule (P49); power rule. · **Assumptions:** a uniform tube (K_p, A₀ constant).
- **Steps:**
  1. **did** Write the transmural pressure with P · **tex** $p-p_e=K_p\,\mathrm P(\alpha),\quad\alpha=\dfrac A{A_0}$ · **why** Definitions (16.186). · **plain** The tube
     law in dimensionless form.
  2. **did** Change variable from A to α · **tex** $c^2=\dfrac{\alpha K_p}{\rho}\dfrac{d\mathrm P}{d\alpha}$ · **why** Chain rule: $\dfrac{d(p-p_e)}{dA}=\dfrac{K_p}{A_0}\dfrac{d\mathrm P}{d\alpha}$, and $A/A_0=\alpha$. ·
     **plain** The wave speed is set by the slope of the tube law.
  3. **did** Differentiate the tube law · **tex** $\mathrm P=1-\alpha^{-n}\ \Rightarrow\ \dfrac{d\mathrm P}{d\alpha}=n\alpha^{-n-1}$ · **why** (16.188); power rule. P is
     negative for α < 1 but its slope is positive (trap T9). · **plain** A collapsed tube gets stiffer as it flattens.
  4. **did** Substitute · **tex** $c^2=\dfrac{nK_p\alpha^{-n}}{\rho}$ (16.189) · **why** $\alpha\cdot\alpha^{-n-1}=\alpha^{-n}$. · **plain** The more collapsed, the faster the
     wave.
- **Result:** $c^2=\dfrac{nK_p\alpha^{-n}}{\rho}$ (16.189).
- **Check:** units Pa/(kg m⁻³) = m²/s² ✓; easy numbers α = 0.25, n = 1.5, K_p = 250 Pa, ρ = 1000 kg/m³: c² = 1.5 × 250 × 8/1000 = 3,
  **c = 1.73 m/s**; our vein at α = 0.5: 0.895 m/s, at α = 1: 0.532 m/s (expect) — a few metres per second at most.
- **What it means:** collapsible vessels carry very slow waves, so ordinary blood speeds can reach S = 1. Fails at α > 1 where
  the fit (16.188) is not meant for distension.
- **Traps:** P < 0 when collapsed, dP/dα > 0 (T9).

### D29 · The differential relations (16.196)–(16.202), (16.201) corrected — ★★, 9 steps, in C15 (notebook · `collapsible_tube_choking`)
- **Goal:** turn the tube law, the wave speed and the speed index into relations between small changes (differentials) along
  the tube. · **Start:** α = A/A₀(x); $\mathrm P=\dfrac{p-p_e}{K_p(x)}=\mathrm P(\alpha,x)$; $c^2(\alpha,x)=\dfrac\alpha\rho K_p\dfrac{\partial\mathrm P}{\partial\alpha}\Big\rvert_x$; $S=\dfrac uc$ (16.195). · **Plan:** • log-differentiate
  α • differentiate the tube law • expand dP • log-differentiate c² • log-differentiate S². · **Tools:** total differential (P49);
  logarithmic differentiation (P52, reminder); product rule (P38). · **Assumptions:** steady flow; P smooth in α and x.
- **Steps:**
  1. **did** Log-differentiate the area ratio · **tex** $\dfrac{d\alpha}\alpha=\dfrac{dA}A-\dfrac{dA_0}{A_0}$ (16.196) · **why** ln α = ln A − ln A₀ (P52). The book
     writes dA_e for dA₀ just before (slip S22). · **plain** The area ratio changes when the area or the reference area does.
  2. **did** Differentiate the tube law · **tex** $dp=dp_e+K_p\,d\mathrm P+\mathrm P\,dK_p$ (16.197) · **why** $p=p_e+K_p\mathrm P$; product rule (P38). · **plain**
     Pressure changes with the outside pressure, the shape and the stiffness.
  3. **did** Write c² with α · **tex** $c^2(\alpha,x)=\dfrac\alpha\rho K_p\dfrac{\partial\mathrm P}{\partial\alpha}\Big\rvert_x$ (16.198) · **why** (16.185) at fixed x, as in D28 step 2. · **plain**
     The local wave speed.
  4. **did** Expand dP · **tex** $d\mathrm P=\dfrac{\partial\mathrm P}{\partial\alpha}d\alpha+\dfrac{\partial\mathrm P}{\partial x}dx$ (16.199) · **why** Total differential of P(α, x) (P49). · **plain** P
     changes with the shape and with position.
  5. **did** Substitute into dp · **tex** $dp=dp_e+\rho c^2\dfrac{d\alpha}\alpha+K_p\dfrac{\partial\mathrm P}{\partial x}dx+\mathrm P\,dK_p$ (16.200) · **why** $K_p\dfrac{\partial\mathrm P}{\partial\alpha}d\alpha=\rho c^2\dfrac{d\alpha}\alpha$ by
     (16.198). · **plain** The shape change enters through ρc².
  6. **did** Log-differentiate c² · **tex** $2\dfrac{dc}c=\dfrac{d\alpha}\alpha+\dfrac{dK_p}{K_p}+\dfrac{d(\partial\mathrm P/\partial\alpha)}{\partial\mathrm P/\partial\alpha}$ · **why** ln c² = ln α + ln K_p + ln P_α − ln ρ, ρ
     constant (P52). · **plain** Three causes of a change in wave speed.
  7. **did** Expand the slope's differential · **tex** $d\left(\dfrac{\partial\mathrm P}{\partial\alpha}\right)=\dfrac{\partial^2\mathrm P}{\partial\alpha^2}d\alpha+\dfrac{\partial}{\partial x}\left(\dfrac{\partial\mathrm P}{\partial\alpha}\right)dx$ · **why** Total
     differential of P_α(α, x) (P49). · **plain** The slope changes with the shape and along the tube.
  8. **did** Collect · **tex** $2\dfrac{dc}c=\left(1+\dfrac{\alpha\,\partial^2\mathrm P/\partial\alpha^2}{\partial\mathrm P/\partial\alpha}\right)\dfrac{d\alpha}\alpha+\dfrac{dK_p}{K_p}+\dfrac{\alpha K_p}{\rho c^2}\dfrac{\partial}{\partial x}\left(\dfrac{\partial\mathrm P}{\partial\alpha}\right)dx$ (16.201, corrected) ·
     **why** Divide step 7 by P_α and use $1/\mathrm P_\alpha=\alpha K_p/(\rho c^2)$ (16.198). The book prints $\partial/\partial x(\partial\mathrm P/\partial x)$ — the inner
     derivative is wrong (slip S23). · **plain** How the wave speed changes.
  9. **did** Log-differentiate S² and use continuity · **tex** $\dfrac{dS^2}{S^2}=-2\dfrac{d\alpha}\alpha-2\dfrac{dA_0}{A_0}-2\dfrac{dc}c$ (16.202) · **why** (16.195):
     $dS^2/S^2=2du/u-2dc/c$; continuity (16.190) $du/u=-dA/A=-d\alpha/\alpha-dA_0/A_0$ (step 1). · **plain** S grows when the tube
     narrows or the wave slows.
- **Result:** (16.196)–(16.202) with (16.201) corrected.
- **Check:** a uniform tube with P = P(α): (16.201) becomes $2\,dc/c=(1+\alpha\mathrm P''/\mathrm P')d\alpha/\alpha$; for the power law
  $\alpha\mathrm P''/\mathrm P'=-(n+1)$, so $2\,dc/c=-n\,d\alpha/\alpha$ — the log-derivative of (16.189) ✓; `ch16.collapsible_differentials_sympy()`
  shows the printed (16.201) failing for a P that depends on x.
- **What it means:** every quantity's change is a sum of the "drivers" dA₀, dp_e, dz, friction, dK_p — ready for D30. Fails
  where the 1-D picture fails (sharp collapse, contact of the walls).
- **Traps:** S23; S22; the partial derivatives are at fixed x or fixed α.

### D30 · Shapiro's equations (16.204)–(16.206) and the condition at S = 1 (rule b: the book lists them) — ★★★, 15 steps (curation: 14), in C15 (notebook · `collapsible_tube_choking`)
- **Goal:** write the change of the area ratio and of the speed index along a collapsible tube in terms of the drivers, and
  see why S = 1 is special. · **Start:** $\dfrac{dA}A+\dfrac{du}u=0$ (16.190) and $-A\,dp-\tau_ws\,dx-\rho gA\,dz=\rho Au\,du$ (16.191). · **Plan:** • momentum in
  differential form • substitute dp and dA/A • collect dα: the factor 1 − S² • the same for dS² • read S = 1. · **Tools:** a small
  linear system of differentials (Ch. 15 D12/D22, P362, reminder); logarithmic differentiation (P52); sympy (P40). · **Assumptions:**
  steady 1-D flow; the tube law depends on α only (step 4); laminar friction (step 3).
- **Steps:**
  1. **did** Use continuity in the momentum law · **tex** $\rho Au\,du=-\rho Au^2\dfrac{dA}A$ · **why** (16.190): $du=-u\,dA/A$. · **plain** Where the
     tube narrows the blood speeds up.
  2. **did** Divide by A and rearrange · **tex** $d(p+\rho gz)+\dfrac{\tau_ws\,dx}A-\rho u^2\dfrac{dA}A=0$ (16.194) · **why** (16.191) with step 1, divided by A;
     ρ constant so $dp+\rho g\,dz=d(p+\rho gz)$. · **plain** Pressure-plus-height falls by friction and changes with the area.
  3. **did** Name the friction term · **tex** $\dfrac{\tau_ws\,dx}A=\mathrm RQ\,dx$ · **why** R is the viscous resistance per unit length (the book's R
     in (16.204)); for laminar flow in area A, $\mathrm R=8\pi\mu/A^2$ (ours, Part C prm). · **plain** Friction costs pressure in
     proportion to the flow.
  4. **did** Take a tube law that depends on α only · **tex** $p-p_e(x)=K_p(x)\mathrm P(\alpha)$ (16.203), so $dp=dp_e+\rho c^2\dfrac{d\alpha}\alpha+\mathrm P\,dK_p$ ·
     **why** The book's case: ∂P/∂x = 0 in (16.200). · **plain** The tube's shape law is the same everywhere; only the outside
     pressure and the stiffness may vary.
  5. **did** Substitute dp and dA/A into (16.194) · **tex** $dp_e+\rho c^2\dfrac{d\alpha}\alpha+\mathrm P\,dK_p+\rho g\,dz+\mathrm RQ\,dx-\rho u^2\left(\dfrac{d\alpha}\alpha+\dfrac{dA_0}{A_0}\right)=0$ ·
     **why** Step 4 and (16.196): $dA/A=d\alpha/\alpha+dA_0/A_0$. · **plain** One equation with dα and the drivers.
  6. **did** Collect the dα terms first · **tex** $\rho(c^2-u^2)\dfrac{d\alpha}\alpha=\rho u^2\dfrac{dA_0}{A_0}-\left[dp_e+\rho g\,dz+\mathrm RQ\,dx+\mathrm P\,dK_p\right]$ · **why** Move
     everything else right; form the whole right side before dividing (trap T10, the ch15 lesson). · **plain** The drivers on one
     side, the response on the other.
  7. **did** Divide by $\rho c^2\,dx$ · **tex** $(1-S^2)\dfrac1\alpha\dfrac{d\alpha}{dx}=\dfrac{S^2}{A_0}\dfrac{dA_0}{dx}-\dfrac1{\rho c^2}\left[\dfrac{dp_e}{dx}+\rho g\dfrac{dz}{dx}+\mathrm RQ+\mathrm P\dfrac{dK_p}{dx}\right]$ (16.204) ·
     **why** ρc² > 0 and S = u/c (16.195). · **plain** Like a nozzle: the sign of 1 − S² decides which way the tube responds.
  8. **did** Simplify (16.201) for P(α) · **tex** $2\dfrac{dc}c=(\mathrm M-2)\dfrac{d\alpha}\alpha+\dfrac{dK_p}{K_p},\quad\mathrm M=3+\dfrac{\alpha\,\partial^2\mathrm P/\partial\alpha^2}{\partial\mathrm P/\partial\alpha}$ (16.206) · **why** ∂P/∂x = 0
     removes the last term of (16.201, corrected); $1+\alpha\mathrm P''/\mathrm P'=\mathrm M-2$. · **plain** M measures how the wave speed changes
     with the shape.
  9. **did** Substitute into (16.202) · **tex** $\dfrac{dS^2}{S^2}=-\mathrm M\dfrac{d\alpha}\alpha-2\dfrac{dA_0}{A_0}-\dfrac{dK_p}{K_p}$ · **why** $-2\,d\alpha/\alpha-(\mathrm M-2)\,d\alpha/\alpha=-\mathrm M\,d\alpha/\alpha$. ·
     **plain** S responds to the shape through M.
  10. **did** Multiply by $1-S^2$ · **tex** $(1-S^2)\dfrac{dS^2}{S^2}=-\mathrm M(1-S^2)\dfrac{d\alpha}\alpha-(1-S^2)\left[2\dfrac{dA_0}{A_0}+\dfrac{dK_p}{K_p}\right]$ · **why** So that
      (16.204) can replace $(1-S^2)d\alpha/\alpha$. · **plain** Ready to eliminate dα.
  11. **did** Insert (16.204) · **tex** $(1-S^2)\dfrac{dS^2}{S^2}=-\mathrm M\left\{S^2\dfrac{dA_0}{A_0}-\dfrac{[\dots]}{\rho c^2}\right\}-(1-S^2)\left[2\dfrac{dA_0}{A_0}+\dfrac{dK_p}{K_p}\right]$ · **why**
      [⋯] is the bracket of step 6, $dp_e+\rho g\,dz+\mathrm RQ\,dx+\mathrm P\,dK_p$. · **plain** Only drivers left on the right.
  12. **did** Collect the dA₀ terms · **tex** $\dfrac{dA_0}{A_0}\left[-2+(2-\mathrm M)S^2\right]$ · **why** $-\mathrm MS^2-2+2S^2$. · **plain** The area's effect on S.
  13. **did** Collect the dK_p terms · **tex** $\dfrac{dK_p}{\rho c^2}\left[\mathrm M\mathrm P-(1-S^2)\alpha\dfrac{d\mathrm P}{d\alpha}\right]$ · **why** From M/(ρc²)·P dK_p and
      $-(1-S^2)dK_p/K_p$, with $1/K_p=\alpha\mathrm P'/(\rho c^2)$ by (16.198). · **plain** The stiffness's effect on S.
  14. **did** Assemble per unit length · **tex** $(1-S^2)\dfrac1{S^2}\dfrac{dS^2}{dx}=\dfrac1{A_0}\dfrac{dA_0}{dx}[-2+(2-\mathrm M)S^2]+\dfrac{\mathrm M}{\rho c^2}\left[\dfrac{dp_e}{dx}+\rho g\dfrac{dz}{dx}+\mathrm RQ\right]+\dfrac1{\rho c^2}\dfrac{dK_p}{dx}\left[\mathrm M\mathrm P-(1-S^2)\alpha\dfrac{d\mathrm P}{d\alpha}\right]$ (16.205) ·
      **why** Steps 11–13, divided by dx. · **plain** Shapiro's equation for the speed index.
  15. **did** Read the condition at S = 1 · **tex** $S=1:\quad\dfrac{S^2}{A_0}\dfrac{dA_0}{dx}=\dfrac1{\rho c^2}\left[\dfrac{dp_e}{dx}+\rho g\dfrac{dz}{dx}+\mathrm RQ+\mathrm P\dfrac{dK_p}{dx}\right]$ · **why**
      The left side of (16.204) vanishes; a finite dα/dx needs the right side to vanish too. So a smooth passage through S = 1 needs
      dA₀/dx > 0, or $d(p_e+\rho gz)/dx<0$, or $\mathrm P\,dK_p/dx<0$ — with P < 0 that is $dK_p/dx>0$ (trap T9; the text's "decrease" is S24). ·
      **plain** Like a Laval throat: S = 1 can only be crossed where the tube widens or is helped.
- **Result:** (16.204), (16.205), (16.206); at S = 1 the drivers must balance. For the power law $\mathrm M=2-n=\tfrac12$.
- **Check:** for $\mathrm P=1-\alpha^{-n}$: $\mathrm P'=n\alpha^{-n-1}$, $\mathrm P''=-n(n+1)\alpha^{-n-2}$, so $\mathrm M=3-(n+1)=2-n=0.5$ ✓; friction only
  (S < 1): dα/dx < 0 — the tube collapses downstream and S rises toward 1 (choking); **expect** our vein, L = 0.1 m, inlet open:
  the outlet reaches S = 1 at **Q = 19.99 mL/s**, α = 0.312, p₂ − p_e = −946 Pa (−7.1 mmHg); lower p₂ further and Q stays there.
- **sympy check (`check_src`, executed):**
  ```python
  import sympy as sp                                                       # symbolic algebra
  dal, dA0, dpe, dz, dx, dKp, RQ, S2, c2, P, Pa, Paa, Kp, al, A0, g = sp.symbols(
      'dalpha dA0 dpe dz dx dKp RQ S2 c2 P Pa Paa Kp alpha A0 g', real=True)   # differentials and coefficients
  rho = sp.symbols('rho', positive=True)
  dp = dpe + rho*c2*dal/al + P*dKp                                         # step 4: (16.200) with ∂P/∂x = 0
  mom = sp.Eq(dp + rho*g*dz + RQ*dx - rho*S2*c2*(dal/al + dA0/A0), 0)      # step 5: (16.194), u² = S²c²
  dal_sol = sp.solve(mom, dal)[0]                                          # steps 6–7
  t204 = al*(S2*dA0/A0 - (dpe + rho*g*dz + RQ*dx + P*dKp)/(rho*c2))/(1 - S2)
  print(sp.simplify(dal_sol - t204))                                       # (16.204) → 0
  M = 3 + al*Paa/Pa                                                        # (16.206)
  dS2 = -2*dal/al - 2*dA0/A0 - ((1 + al*Paa/Pa)*dal/al + dKp/Kp)           # (16.202) with (16.201, corrected), steps 8–9
  lhs = sp.expand((1 - S2)*dS2.subs(dal, dal_sol))                         # steps 10–11
  rhs = (dA0/A0*(-2 + (2 - M)*S2) + M/(rho*c2)*(dpe + rho*g*dz + RQ*dx)
         + dKp/(rho*c2)*(M*P - (1 - S2)*al*Pa))                            # step 14: (16.205)
  print(sp.simplify((lhs - rhs).subs(c2, al*Kp*Pa/rho)))                   # with ρc² = αK_pP′ (16.198) → 0
  ```
  printed: `0` and `0`.
- **What it means:** S plays the Mach number: below 1 friction collapses the tube and speeds the flow, until S = 1 at the outlet;
  then the flow no longer feels the downstream pressure — flow limitation (Starling Case 1), the floppy-tube twin of a choked
  nozzle. Fails when the walls touch or the flow separates at the throat.
- **Traps:** T10 (form the driver, then divide by 1 − S²); T9 and S24 (the sign of P dK_p/dx); M is not a Mach number (T2).

### D31 · The Casson stress line (16.209), the plug radius (16.210) and the minimum pressure drop — ★, 5 steps, in C16 (notebook · `casson_blood_flow`)
- **Goal:** find where a yield-stress fluid in a pipe shears and where it moves as a solid plug, and the smallest pressure drop
  that makes it flow at all. · **Start:** a fluid cylinder of radius r and length Δx in steady flow. · **Plan:** • force balance •
  the stress line • where it meets τ_y • when the plug fills the pipe. · **Tools:** free body of a fluid cylinder (Ch. 8 D07,
  recap). · **Assumptions:** steady, fully developed (step 1).
- **Steps:**
  1. **did** Balance forces on a fluid cylinder · **tex** $\tau_{rx}\,2\pi r\,\Delta x=-\pi r^2\left(p\rvert_{x+\Delta x}-p\rvert_x\right)$ (16.208) · **why** Steady, fully developed:
     no acceleration, so the shear on the side balances the pressure difference on the ends (Ch. 8 D07). · **plain** The pressure
     push on a core of radius r is held by the shear on its surface.
  2. **did** Let Δx shrink · **tex** $\tau_{rx}(r)=-\dfrac r2\dfrac{dp}{dx}=\dfrac{(p_1-p_2)r}{2L}$ (16.209, first member's sign corrected) · **why** Divide by 2πrΔx;
     the limit is dp/dx = −(p₁ − p₂)/L (16.14). The book's first member lacks the minus (slip S25). · **plain** Stress grows
     linearly from the axis — whatever the fluid.
  3. **did** Evaluate at the wall · **tex** $\tau_w=\dfrac{(p_1-p_2)a}{2L}$ · **why** r = a. · **plain** The largest stress is at the wall.
  4. **did** Find where the stress reaches the yield stress · **tex** $a_c=\dfrac{2L\tau_y}{p_1-p_2}=a\left(\dfrac{\tau_y}{\tau_w}\right)$ (16.210) · **why** τ_rx(a_c) = τ_y; inside
     (r < a_c) τ < τ_y, so the shear rate is zero by (16.207): that core moves as a solid plug. · **plain** A plug of radius a_c rides
     on a sheared ring.
  5. **did** Ask when the plug fills the pipe · **tex** $a_c\ge a\iff p_1-p_2\le\Delta p_{\min}=\dfrac{2L\tau_y}{a}$ · **why** (16.210) with a_c = a: the wall stress is
     below yield everywhere, so nothing shears and nothing flows (trap T21). The book quotes this minimum in dyn/cm² 1000× too
     small (slip S26). · **plain** Below this drop the blood does not move.
- **Result:** (16.209), (16.210), $\Delta p_{\min}=2L\tau_y/a$.
- **Check:** units Pa·m/m = Pa ✓; **expect** τ_y = 4.8 mPa, L/a = 1000: Δp_min = **9.6 Pa = 0.072 mmHg** — tiny for a whole
  vessel, so in the circulation the yield stress only matters where pressure differences nearly vanish.
- **What it means:** a yield stress does not stop blood in a pipe; it freezes the centre. Fails if the fluid is not homogeneous
  (small vessels, C05).
- **Traps:** S25, S26; a_c ≥ a means no flow (T21).

### D32 · The Casson profile (16.211)–(16.216) and the plug speed (16.217)–(16.218) — ★★, 11 steps, in C16 (notebook · `casson_blood_flow`)
- **Goal:** find the velocity in the sheared ring and the speed of the plug. · **Start:** $\tau^{1/2}=\tau_y^{1/2}+K_c\dot\gamma^{1/2}$ (16.207) and
  $\tau=\dfrac{\Delta p}{2L}r$ (16.209). · **Plan:** • invert the law for the shear rate • write it in r • integrate • no slip • rescale with τ_w •
  evaluate at the plug's edge. · **Tools:** algebra with square roots and $\int r^{1/2}dr$ (P393, primer in C16). · **Assumptions:**
  as D31; τ ≥ τ_y in the ring.
- **Steps:**
  1. **did** Invert the law in the ring · **tex** $\sqrt{\dot\gamma}=\dfrac{\sqrt\tau-\sqrt{\tau_y}}{K_c}$ · **why** For τ ≥ τ_y subtract √τ_y and divide by K_c. ·
     **plain** The shear rate grows with the stress above yield.
  2. **did** Square · **tex** $\dot\gamma=-\dfrac{du}{dr}=\dfrac{\tau+\tau_y-2\sqrt{\tau\tau_y}}{K_c^2}$ (16.211) · **why** $(\sqrt x-\sqrt y)^2=x+y-2\sqrt{xy}$ (P393); the
     shear rate is −du/dr because u falls outward. · **plain** The slope of the profile.
  3. **did** Write both stresses with r · **tex** $\tau=\dfrac{\Delta p}{2L}r,\qquad\tau_y=\dfrac{\Delta p}{2L}a_c$ · **why** (16.209) and (16.210). · **plain**
     Everything in one length scale.
  4. **did** Substitute · **tex** $-\dfrac{du}{dr}=\dfrac1{2K_c^2}\dfrac{\Delta p}L\left(r+a_c-2\sqrt{ra_c}\right)$ (16.212) · **why** $\sqrt{\tau\tau_y}=\dfrac{\Delta p}{2L}\sqrt{ra_c}$. · **plain** The
     slope vanishes at r = a_c, where the plug begins.
  5. **did** Integrate term by term · **tex** $u=\dfrac1{2K_c^2}\dfrac{\Delta p}L\left(\dfrac43\sqrt{a_cr^3}-\dfrac{r^2}2-a_cr+C\right)$ (16.213) · **why** $\int r\,dr=r^2/2$,
     $\int a_c\,dr=a_cr$, $\int2\sqrt{a_c}r^{1/2}dr=\tfrac43\sqrt{a_c}r^{3/2}$ (P393); the minus of −du/dr flips each sign. · **plain** The profile, up to a
     constant.
  6. **did** Apply no slip · **tex** $C=-\left(\dfrac43\sqrt{a_ca^3}-\dfrac{a^2}2-a_ca\right)$ (16.214) · **why** u(a) = 0. · **plain** The constant from the wall.
  7. **did** Put C back and tidy · **tex** $u=\dfrac1{4K_c^2}\dfrac{\Delta p}L\left[(a^2-r^2)-\dfrac83\sqrt{a_c}\left(\sqrt{a^3}-\sqrt{r^3}\right)+2a_c(a-r)\right]$ (16.215) · **why** Group by
     powers of r; take a factor ½ into the prefactor. · **plain** The velocity in the sheared ring.
  8. **did** Scale with the wall stress · **tex** $\dfrac{\Delta p}L=\dfrac{2\tau_w}a,\qquad a_c=a\dfrac{\tau_y}{\tau_w}$ · **why** D31 step 3 and (16.210). · **plain** Two
     numbers set the shape: τ_w and τ_y.
  9. **did** Rewrite the profile · **tex** $u=\dfrac{a\tau_w}{2K_c^2}\left\{\left[1-\left(\dfrac ra\right)^2\right]-\dfrac83\sqrt{\dfrac{\tau_y}{\tau_w}}\left[1-\left(\dfrac ra\right)^{3/2}\right]+2\dfrac{\tau_y}{\tau_w}\left(1-\dfrac ra\right)\right\}$ (16.216) · **why**
     Factor a² out of each bracket term; $\sqrt{a_c}\sqrt{a^3}=a^2\sqrt{\tau_y/\tau_w}$. · **plain** The profile in units of the wall
     stress.
  10. **did** Evaluate at the plug's edge · **tex** $\dfrac ra=\dfrac{a_c}a=\dfrac{\tau_y}{\tau_w}$ (16.217) · **why** The plug moves rigidly with the
      speed of its edge. · **plain** Where the ring meets the plug.
  11. **did** Factor the plug speed · **tex** $u_c=\dfrac1{4K_c^2}\dfrac{\Delta p}L\left(\sqrt a-\sqrt{a_c}\right)^3\left(\sqrt a+\dfrac13\sqrt{a_c}\right)$ (16.218) · **why** With s = √a,
      q = √a_c the bracket of (16.215) at r = a_c is $s^4-\tfrac83s^3q+2s^2q^2-\tfrac13q^4$, which equals $(s-q)^3\left(s+\tfrac q3\right)$
      (multiply out). · **plain** The plug's speed.
- **Result:** (16.215)/(16.216) in the ring; $u_c$ (16.218) in the plug.
- **Check:** τ_y → 0 (a_c → 0): Poiseuille with μ = K_c² ✓; a_c = a: u_c = 0 ✓; the sympy test checks u(a_c) = u_c (and
  du/dr = 0 there); **expect** a_c/a = 0.25: u_c is 0.146 of the Poiseuille centreline speed at μ = K_c².
- **What it means:** the plug is flat, the ring carries all the shear; the plug's speed falls fast as a_c grows. Fails as D31.
- **Traps:** Δp/L = 2τ_w/a; the constant comes from u(a) = 0, not from the plug.

### D33 · The Casson flux (16.220) from (16.219) (rule b: "after considerable algebra") — ★★★, 10 steps, in C16 (notebook · `casson_blood_flow`)
- **Goal:** add up the plug and the ring to get the flow rate of a Casson fluid, and compare with Poiseuille's. · **Start:**
  $Q=\pi a_c^2u_c+\displaystyle\int_{a_c}^a2\pi ru\,dr$ (16.219) with (16.215) and (16.218). · **Plan:** • three ring integrals • scale with a⁴ • write both parts as polynomials in
  $t=\sqrt{a_c/a}$ • add: most powers cancel. · **Tools:** $\int r^{n/2}dr$ (P393); a polynomial in √ξ; sympy (P40). · **Assumptions:** as D32.
- **Steps:**
  1. **did** Split the flux into plug and ring · **tex** $Q=\pi a_c^2u_c+\displaystyle\int_{a_c}^a2\pi ru\,dr$ (16.219) · **why** The plug (radius a_c) moves rigidly at u_c;
     the ring is added up ring by ring. · **plain** Two parts.
  2. **did** Pull out the ring profile's prefactor · **tex** $\displaystyle\int_{a_c}^a2\pi ru\,dr=\dfrac{\pi}{2K_c^2}\dfrac{\Delta p}L\int_{a_c}^ar\left[(a^2-r^2)-\tfrac83\sqrt{a_c}\left(a^{3/2}-r^{3/2}\right)+2a_c(a-r)\right]dr$ ·
     **why** (16.215); $2\pi/4=\pi/2$. · **plain** Three integrals to do.
  3. **did** Integrate the parabolic term · **tex** $\displaystyle\int_{a_c}^ar(a^2-r^2)\,dr=\dfrac{(a^2-a_c^2)^2}4$ · **why** The substitution of D06 step 2. · **plain**
     The Poiseuille-like part.
  4. **did** Integrate the root term · **tex** $\displaystyle\int_{a_c}^ar\left(a^{3/2}-r^{3/2}\right)dr=a^{3/2}\dfrac{a^2-a_c^2}2-\dfrac27\left(a^{7/2}-a_c^{7/2}\right)$ · **why** $\int r^{5/2}dr=\tfrac27r^{7/2}$
     (P393). · **plain** Where the 16/7 will come from.
  5. **did** Integrate the linear term · **tex** $\displaystyle\int_{a_c}^ar(a-r)\,dr=a\dfrac{a^2-a_c^2}2-\dfrac{a^3-a_c^3}3$ · **why** Elementary. · **plain** The yield
     correction's part.
  6. **did** Scale the ring with $a^4$ · **tex** $\displaystyle\int_{a_c}^a2\pi ru\,dr=\dfrac{\pi a^4}{8K_c^2}\dfrac{\Delta p}L\cdot4I(\xi),\quad\xi=\dfrac{a_c}a$ · **why** Every term of steps 3–5 times
     its coefficient has total power a⁴; factor it out, with $I=\tfrac14(1-\xi^2)^2-\tfrac83\xi^{1/2}\left[\tfrac12(1-\xi^2)-\tfrac27(1-\xi^{7/2})\right]+2\xi\left[\tfrac12(1-\xi^2)-\tfrac13(1-\xi^3)\right]$. ·
     **plain** The ring's flux as Poiseuille's times a pure number.
  7. **did** Expand the ring in $t=\sqrt\xi$ · **tex** $4I=1-\dfrac{16}7t+\dfrac43t^2-2t^4+\dfrac{16}3t^5-4t^6+\dfrac{13}{21}t^8$ · **why** Multiply out (sympy in the check cell);
     t = √ξ makes every power an integer. · **plain** A polynomial in t.
  8. **did** Write the plug in the same units · **tex** $\pi a_c^2u_c=\dfrac{\pi a^4}{8K_c^2}\dfrac{\Delta p}L\cdot2t^4(1-t)^3\left(1+\dfrac t3\right)=\dfrac{\pi a^4}{8K_c^2}\dfrac{\Delta p}L\left[2t^4-\dfrac{16}3t^5+4t^6-\dfrac23t^8\right]$ ·
     **why** (16.218): $(\sqrt a-\sqrt{a_c})^3(\sqrt a+\tfrac13\sqrt{a_c})=a^2(1-t)^3(1+t/3)$; $a_c^2=t^4a^2$; $\tfrac{\pi}{4}=\tfrac{\pi}{8}\cdot2$; multiply out. ·
     **plain** The plug's flux, also a polynomial in t.
  9. **did** Add the two polynomials · **tex** $4I+2t^4(1-t)^3\left(1+\dfrac t3\right)=1-\dfrac{16}7t+\dfrac43t^2-\dfrac1{21}t^8$ · **why** The t⁴, t⁵ and t⁶ terms cancel
     exactly; for t⁸, $\tfrac{13}{21}-\tfrac23=-\tfrac1{21}$. · **plain** Almost everything cancels.
  10. **did** Write it with $a_c/a$ · **tex** $Q=\dfrac\pi8\dfrac1{K_c^2}\dfrac{\Delta p}La^4\left[1-\dfrac{16}7\left(\dfrac{a_c}a\right)^{1/2}+\dfrac43\left(\dfrac{a_c}a\right)-\dfrac1{21}\left(\dfrac{a_c}a\right)^4\right]$ (16.220) · **why**
      $t=(a_c/a)^{1/2}$, $t^2=a_c/a$, $t^8=(a_c/a)^4$. · **plain** Poiseuille's flux with μ = K_c², times a bracket below 1.
- **Result:** (16.220).
- **Check:** ξ → 0: bracket → 1, Poiseuille with μ = K_c² ✓; ξ → 1: $1-\tfrac{16}7+\tfrac43-\tfrac1{21}=0$, Q → 0 ✓; **expect** ξ = 0.25: bracket
  **0.190**; ξ = 0.01: 0.785 — even a plug 1 % of the radius cuts the flow by a fifth (the √ξ term); the trapezoid of 2πru over plug +
  ring equals `ch16.casson_pipe_flow_rate` to 10⁻⁶.
- **sympy check (`check_src`, executed):**
  ```python
  import sympy as sp                                                       # symbolic algebra
  a, ac, Kc, G, r = sp.symbols('a a_c K_c G r', positive=True)             # G = Δp/L
  u = G/(4*Kc**2)*((a**2 - r**2) - sp.Rational(8, 3)*sp.sqrt(ac)*(a**sp.Rational(3, 2) - r**sp.Rational(3, 2))
                   + 2*ac*(a - r))                                         # (16.215)
  uc = G/(4*Kc**2)*(sp.sqrt(a) - sp.sqrt(ac))**3*(sp.sqrt(a) + sp.sqrt(ac)/3)   # (16.218)
  print(sp.simplify(u.subs(r, ac) - uc))                                   # the plug speed matches the ring's edge → 0
  Q = sp.pi*ac**2*uc + sp.integrate(2*sp.pi*r*u, (r, ac, a))               # (16.219), steps 1–6
  xi = ac/a
  Q220 = sp.pi/8/Kc**2*G*a**4*(1 - sp.Rational(16, 7)*sp.sqrt(xi) + sp.Rational(4, 3)*xi - xi**4/21)   # (16.220)
  print(sp.simplify(sp.expand(Q - Q220)))                                  # → 0
  t = sp.symbols('t', positive=True)                                       # t = √ξ (steps 7–9)
  ring = sp.simplify((Q - sp.pi*ac**2*uc).subs(ac, t**2*a)/(sp.pi*a**4*G/(8*Kc**2)))   # ring flux / Poiseuille's = 4I(t)
  print(sp.expand(sp.simplify(ring)))                                      # 13t⁸/21 − 4t⁶ + 16t⁵/3 − 2t⁴ + 4t²/3 − 16t/7 + 1
  ```
  printed: `0`, `0` and `13*t**8/21 - 4*t**6 + 16*t**5/3 - 2*t**4 + 4*t**2/3 - 16*t/7 + 1` (step 7).
- **What it means:** the yield stress costs flow out of proportion to its size because of the √(a_c/a) term; a small plug
  already matters. Fails as D31.
- **Traps:** do not forget the plug term π a_c² u_c; t = √ξ, not ξ.

### D34 · The three terms of the water potential (16.221) and the tension that lifts sap (ours) — ★, 5 steps, in C17 (notebook · `xylem_water_potential`)
- **Goal:** see why water's "potential" adds a pressure, a solute and a height term, and what pressure a still column needs at
  height z. · **Start:** $\psi=p-RTc+\rho gz$ (16.221) and van 't Hoff's $\Pi=cR_uT$ (P384). · **Plan:** • name the three terms • check their units •
  hold a still column • solve for p • read the sign. · **Tools:** osmosis and van 't Hoff's law (P384, primer in C17); the
  hydrostatic law $dp/dz=-\rho g$ (1.8) (Ch. 1 D05, recap); units of physiology (P394). · **Assumptions:** a dilute ideal
  solution (step 1); no flow (step 3).
- **Steps:**
  1. **did** Name the three ways water's energy changes · **tex** $\psi=\underbrace{p}_{\text{pressure}}\;\underbrace{-\,R_uTc}_{\text{solutes}}\;+\underbrace{\rho gz}_{\text{height}}$ (16.221) · **why**
     Per unit volume of water: pressure work p; dissolved solutes lower water's free energy by the osmotic pressure Π = cR_uT (P384);
     lifting by z costs ρgz. Water moves from high ψ to low ψ. · **plain** One ledger for pressure, salt and height.
  2. **did** Check the units · **tex** $[R_uTc]=\dfrac{\mathrm J}{\mathrm{mol\,K}}\cdot\mathrm K\cdot\dfrac{\mathrm{mol}}{\mathrm m^3}=\dfrac{\mathrm J}{\mathrm m^3}=\mathrm{Pa}$ · **why** Each term must be a pressure.
     In the book's units ρgz comes out in dyn/cm² and must be divided by 10⁶ to be in bar (slip S30). · **plain** All three are
     pascals.
  3. **did** Hold a still column of pure water · **tex** $\psi(z)=\psi_{soil},\quad c=0$ · **why** No flow means no gradient of ψ (water moves
     only down ψ); the xylem sap is nearly pure water. · **plain** The column is in balance with the soil.
  4. **did** Solve for the pressure at height z · **tex** $p(z)=\psi_{soil}-\rho gz$ · **why** (16.221) with c = 0 and ψ = ψ_soil; the same as
     $dp/dz=-\rho g$ (1.8). · **plain** The pressure falls by ρg per metre of height.
  5. **did** Read the sign · **tex** $p(z)<0\ \text{for } z>0\ \text{when}\ \psi_{soil}\le0$ · **why** Soil holds its water, so ψ_soil ≤ 0; minus ρgz makes p
     more negative. A negative absolute-gauge pressure is a tension, held by the cohesion of water (trap T22). · **plain** Sap
     hangs from the leaves under tension.
- **Result:** $p(z)=\psi_{soil}-\rho gz<0$ (ours) for a still column; with flow the Poiseuille drop $8\mu uz/a^2$ is added (ours).
- **Check:** units ✓; **expect** 30 m: ρgz = 0.294 MPa; ψ_soil = −0.1 MPa gives p(30 m) = **−0.394 MPa** (−3.9 bar); with sap moving at
  1 mm/s in a 25 µm conduit the viscous drop adds 12.8 kPa per metre: −0.778 MPa at the top; the osmotic term of a 0.3 mol/L
  leaf cell: R_uTc = 0.731 MPa.
- **What it means:** trees lift water by pulling, not pushing: evaporation at the leaf keeps ψ low, and the column hangs under
  tension; capillarity alone (0.59 m in a 25 µm conduit) is far too weak. Fails when the tension exceeds what the column can
  hold (cavitation, embolism).
- **Traps:** ψ < 0 under tension (T22); the book's mixed units (S30); "pressure" here is gauge.

---

## Closing check — design gate (bullets only, so that the ledger parser of Part E is not disturbed)
- **CORE blocks storyboarded: 17** (C01–C17) = the curation's CORE count; every block has plain words, the idea, its primers before
  first use, its notes and derivations, a tiny example, a fluidpy code cell with expectations, a from-scratch or independent check
  (C03, C06, C08's is a wave-equation residual; C11's the numeric ODE route), at least one visual (figure, animation, plotly slider or
  explainer) and a "What would change if…".
- **Placement:** every NOTE N01–N113 sits in its curation parent (N24 → C03, N29 N30 → C05, N48 → C07, N40–N42 → C06 as decision 4
  says); RECAP R01 (C02), R02 (C03), R03–R08 (C04), R09–R10 (C08), R11 (C09), R12 (C11), R13 (C15), each followed by the re-entry line;
  SKIP S01, S02 as pointers in A.5. Every book section 16.1–16.4 has an `nb.section`.
- **Derivations written out: 34** (D01–D34) with step counts D01 6 · D02 5 · D03 4 · D04 5 · D05 9 · D06 8 · D07 6 · D08 6 · D09 4 · D10 9 ·
  D11 7 · D12 7 · D13 13 · D14 6 · D15 8 · D16 7 · D17 7 · D18 9 · D19 10 · D20 10 · D21 13 · D22 8 · D23 3 · D24 11 · D25 6 · D26 6 ·
  D27 14 · D28 4 · D29 9 · D30 15 · D31 5 · D32 11 · D33 10 · D34 5 (total 265 steps); four differ from the curation's estimate (D07 6/5,
  D10 9/8, D21 13/12, D30 15/14 — one move per step forced a split). The eight ★★★ `check_src` listings were executed exactly as
  printed and passed.
- **Explainers: 10 + backup**, headings `### E1 · casson_blood_flow` … `### E10 · xylem_water_potential` and `### B1 ·
  fahraeus_lindqvist_layer`; each with views, controls, presets (distinguishing value in `set`), status, readouts, ≥ 2 depth features
  beyond Explain + Code, an Explain tab with a "Reading the current setting" section, a Derivation tab listing its D ids with the
  notebook's step titles, a Code listing, 6–8 walkthrough steps, 4–6 equations written out with their numbers, ≥ 3 check questions,
  ≥ 2 parity rows and a fit plan. Every ★★★ D row whose CORE has an explainer is in that explainer (D13 → E4, D24 → E6, D27 → E8,
  D30 → E9, D33 → E1); D19–D21 are notebook-only (C11 has no explainer).
- **Primers: 21** (P374–P394), each placed before its first use; 36 one-line reminders of earlier primers (Part E).
- **Equations:** every book equation named in this file is written out next to its number; nine pages were re-read against the
  analyst's transcription (header). **One correction to the analysis/curation:** the sign of w₁ (D27), with sympy evidence.
- **Public-repo rule:** no string of the private `_forbidden_public` list appears in this file; all worked numbers come from the dicts
  of header convention 10 (ours), and the book's values are named only as "the book's" without numbers.
