# Chapter 15 — Compressible Flow: lesson design
(from `analysis/ch15_curation.md` (A 16 · B 121 · C 14 = 151 rows; CORE 16 · NOTE 114 · RECAP 19 · SKIP 2; 31 derivations D01–D31
(★ 8 · ★★ 19 · ★★★ 4: D13, D17, D21, D30); 12 statements; E1–E9 + backup B1; primers **P354–P373**) and `analysis/ch15.md` (§2
inventory with the LaTeX the analyst read from the page images, §2b derivations a-d01 … a-d45, §4 implementation rows I1–I31, §9
slips #1–#17 and traps T1–T22). 2026-10-13, lesson-designer.

**Equation provenance (honest statement).** Every equation placed here is taken from the analyst's page-image transcription in
`analysis/ch15.md` §2 as carried by the curation. For this design I re-read 12 rendered pages against that transcription (PDF
page = printed page + 27): p761 ($\dfrac1{c^2}\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)^2p'-\dfrac{\partial^2p'}{\partial x_j\partial x_j}=0$ (15.12); $\dfrac1{c^2}\dfrac{\partial^2p'}{\partial t^2}-\dfrac{\partial^2p'}{\partial x^2}=0$ (15.13); $u_1'(x,t)=-\dfrac1{\rho_0}\displaystyle\int\dfrac{\partial p'}{\partial x}\,dt$ (15.14) — printed with $\rho_o$, slip #2; $p'(x,t)=f(x-ct)+g(x+ct)$ (15.15); $u_1'(x,t)=\dfrac1{\rho_0c}\left(f(x-ct)-g(x+ct)\right)$ (15.16)); p767 ($\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31)); p772 (the
printed $u\,du=-dp/\rho=c^2d\rho/\rho$ (15.32) — slip #4 seen: the third member has no minus sign; $d\rho/\rho=-M^2\,du/u$ (15.33); $\dfrac{du}{u}=-\dfrac{1}{1-M^2}\dfrac{dA}{A}$ (15.34); the constant-density
remark printed as $dA/dx+du/dx=0$ — slip #5 seen); p777 ($\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39); $p_1\left(1+\gamma M_1^2\right)=p_2\left(1+\gamma M_2^2\right)$ (unnumbered, §15.6) in the form $p_1+\gamma p_1M_1^2=p_2+\gamma p_2M_2^2$; $M_2^2=\dfrac{(\gamma-1)M_1^2+2}{2\gamma M_1^2+1-\gamma}$ (15.40);
$\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}$ (15.41)); p778 ($\dfrac{T_2}{T_1}=1+\dfrac{2(\gamma-1)}{(\gamma+1)^2}\dfrac{\gamma M_1^2+1}{M_1^2}\left(M_1^2-1\right)$ (15.42); $\dfrac{S_2-S_1}{C_v}=\ln\left[\dfrac{p_2}{p_1}\left(\dfrac{\rho_1}{\rho_2}\right)^{\gamma}\right]$ (15.43) with its second member; $\dfrac{S_2-S_1}{C_v}\cong\dfrac{2\gamma(\gamma-1)}{3(\gamma+1)^2}\left(M_1^2-1\right)^3$ (15.44a)); p789 ($A\equiv\dfrac{M_1^2\left(1+\frac12(\gamma-1)M_1^2+q\right)}{\left(1+\gamma M_1^2-f\right)^2}=\dfrac{M_2^2\left(1+\frac12(\gamma-1)M_2^2\right)}{\left(1+\gamma M_2^2\right)^2}$ (unnumbered, §15.8); $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46)); p790 (the printed
$-\dfrac{dp}{p_1}=\dfrac{1-(\gamma-1)M^2}{1-M^2}df$ (15.47) — slip #11 seen — and its second relation
$\dfrac{du}{u}=-\dfrac{d\rho}{\rho}=\dfrac{p_1}{p}\dfrac{df}{1-M^2}$); p791 (the printed
$\dfrac{dT}{T_1}=\left(\dfrac{1-\gamma M^2}{1-M^2}\right)dq$ and $\dfrac{u\,du}{h_1}=\dfrac{(1-\gamma)M^2}{1-M^2}dq$ (15.48) — slip #12 seen);
p795 ($M_{n1}=M_1\sin\sigma,\ \ M_{n2}=M_2\sin(\sigma-\delta)$ (unnumbered, §15.9); $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2\sin^2\sigma-1\right)$ (15.50); $\dfrac{\rho_2}{\rho_1}=\dfrac{(\gamma+1)M_1^2\sin^2\sigma}{(\gamma-1)M_1^2\sin^2\sigma+2}=\dfrac{\tan\sigma}{\tan(\sigma-\delta)}$ (15.51); $M_2^2\sin^2(\sigma-\delta)=\dfrac{(\gamma-1)M_1^2\sin^2\sigma+2}{2\gamma M_1^2\sin^2\sigma+1-\gamma}$ (15.52); $\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}$ (15.53)); p798 ($\tan\delta=2\cot\sigma\,\dfrac{\gamma+1}{2\gamma}\left(\dfrac{p_2-p_1}{p_1}\right)\dfrac{1}{M_1^2\left(\gamma+1-2\sin^2\sigma\right)+2}$ (15.54); $\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta$ (15.55)); p800 (the lines $d\delta=-(dV/V)\sqrt{M^2-1}$,
$dV/V=dM/M+dc/c$, $dc/c=dT/2T$, $dT/T=-(\gamma-1)M\,dM/[1+((\gamma-1)/2)M^2]$; $d\delta=-\dfrac{\sqrt{M^2-1}}{M}\dfrac{dM}{1+\frac12(\gamma-1)M^2}$ (unnumbered, §15.9); $\nu(M)=\displaystyle\int_1^M\dfrac{\sqrt{M^2-1}}{1+\frac12(\gamma-1)M^2}\dfrac{dM}{M}$ (15.56, integral form) and $\nu(M)=\sqrt{\dfrac{\gamma+1}{\gamma-1}}\tan^{-1}\sqrt{\dfrac{\gamma-1}{\gamma+1}\left(M^2-1\right)}-\tan^{-1}\sqrt{M^2-1}$ (15.56); $\delta_1+\nu(M_1)=\delta_2+\nu(M_2)$ (unnumbered, §15.9); $\dfrac{p-p_\infty}{p_\infty}=\dfrac{\gamma M_\infty^2\,\delta}{\sqrt{M_\infty^2-1}}$ (15.57));
p802 ($L=(p_3-p_2)\,b\cos\alpha$, $D=(p_3-p_2)\,b\sin\alpha$ and $L\cong\dfrac{2\alpha\gamma M_\infty^2p_\infty b}{\sqrt{M_\infty^2-1}},\ \ D\cong\dfrac{2\alpha^2\gamma M_\infty^2p_\infty b}{\sqrt{M_\infty^2-1}}$ (15.58); $C_L\cong\dfrac{4\alpha}{\sqrt{M_\infty^2-1}},\ \ C_D\cong\dfrac{4\alpha^2}{\sqrt{M_\infty^2-1}}$ (15.59); the initial-value formulas $f(x)=\dfrac12\left[F(x)-\dfrac1c\displaystyle\int_0^xG(\hat x)\,d\hat x\right],\ g(x)=\dfrac12\left[F(x)+\dfrac1c\displaystyle\int_0^xG(\hat x)\,d\hat x\right]$ (unnumbered, Exercise 15.1)).
**31 of the 59 numbered equations are on those twelve pages; no difference from the analyst's transcription was found, and the
three physics-changing slips are on the page exactly as the analyst reports.** The other 34 pages were not re-read by the
designer; the 28 numbered equations on them (each written out below where it is used, for example $h_2+\tfrac12u_2^2-h_1-\tfrac12u_1^2=Q$ (15.22), $\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28), $\rho_1u_1=\rho_2u_2,\ \ p_1+\rho_1u_1^2=p_2+\rho_2u_2^2+p_1f,\ \ h_1+\tfrac12u_1^2+h_1q=h_2+\tfrac12u_2^2$ (15.45, with $f=F/(p_1A)$, $q=Q/h_1$) and $\sin\mu=1/M$ (15.49)) are the analyst's transcription.
Equations of earlier chapters are quoted in the written-out forms of the curation: $M\equiv U/c$ (4.111), $c^2\equiv(\partial p/\partial\rho)_s$ (1.19), $\nabla\cdot\mathbf u=-M^2\left(\dfrac{\rho_0}{\rho}\right)\dfrac{D}{Dt}\left(\dfrac{p-p_0}{\rho_0U^2}\right)$ (4.110), $H_2/H_1=\tfrac12\left(-1+\sqrt{1+8Fr_1^2}\right)$ (7.81), $C_L\simeq2\pi(\alpha+\beta)$ (14.12).
Builders render a page (`tools/render_pages.py ch15 --eq 15.M`) before setting any equation not in the first list.

Numbers marked *expect* were computed for this design by a scratch numpy / scipy script written from the equations above (the designer's own twenty-line functions). The implementer's modules appeared on disk while the design was being written; **they were then called with the same inputs and returned the same numbers** (closing check). The builder compares an executed cell against them and **reports** any difference instead of editing the
expectation. The four `check_src` listings of Part F (D13, D17, D21, D30) were executed as written (sympy; 0.1–0.8 s each) and
passed.

**Part C is the contract both the implementer and the builders keep.** It was drafted from analysis §4 (rows I1–I31) and curation §8–§9 and then **reconciled with the modules as found on disk on 2026-10-13** (`fluidpy/core/gasdynamics.py`, `fluidpy/core/acoustics.py`, `fluidpy/ch15_compressible_flow.py`), probed by calling them: every signature and return key in Part C, and every call in Parts A and B, is the code's.

**Binding conventions for every builder (analysis §9, curation decisions 4–9, §8–§10).**
1. **Imports and aliases.** `from fluidpy import ch15_compressible_flow as ch15`; `from fluidpy.core import gasdynamics as GD,
   acoustics as AC`. **`ch15` re-exports every public name of `GD` and `AC`** and the reused functions named in C.0 (the ch11–ch14
   pattern), so explainer parity rows write `ch15.<name>` only. Reused: `from fluidpy.core import thermo as TH, bernoulli as BER,
   similarity as SIM, statics as STAT, airfoils as AF`; `from fluidpy import ch01_introduction as ch01, ch07_gravity_waves as ch07,
   ch14_aerodynamics as ch14`.
2. **Scalar-callable and parity-friendly** (the `shot.py` evaluator knows only `np`, `math`, `ch15`): every public closed form
   accepts Python floats and returns a float, a tuple or a `dict` of floats (arrays only when an array goes in).
3. **Units.** SI and kelvin. Angles in radians inside functions; degrees only in prose, slider labels and through `np.deg2rad`.
   γ is always an argument (default 1.4). R and $C_p$ are passed through `GD.PerfectGas(gamma, R)` whenever a dimensional answer is
   produced. Station subscripts: **1 upstream, 2 downstream** in every signature.
4. **Two conventions on screen wherever a rival standard exists** (the "which $p_o$?" pattern; `ch15.conventions_table()` returns
   this table; the conventions box before C01 shows it; each block repeats the row it needs).

| # | What | The book (and this notebook) | Most other sources | In code | Shown as |
|---|---|---|---|---|---|
| i | subscript 0 | three meanings: a reference value (§15.1), the undisturbed state of the sound wave (§15.2, $u_i=U_i+u_i',\ p=p_0+p',\ \rho=\rho_0+\rho',\ T=T_0+T'$ (15.4)), the stagnation state (from §15.4, $\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28)) | subscript t or "total": $p_t$, $T_t$ | `p0`, `T0` = stagnation from C04 on; `p_ref`, `rho0` in `AC` | a three-row table in the box; "$p_0$ (= $p_t$)" in C04's first legend (trap T1) |
| ii | the star | the control surface $A^*$ in $\displaystyle\int_{A^*}\rho\left(e+\tfrac12u_i^2\right)u_jn_j\,dA=\int_{A^*}u_i\tau_{ij}n_j\,dA-\int_{A^*}q_jn_j\,dA$ (15.20), $\displaystyle\int_{A^*}\rho u_iu_jn_j\,dA=\int_{A^*}\tau_{ij}n_j\,dA$ (15.23); the sonic state in $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31); (ours) the sonic state of a Fanno or Rayleigh process | $A^*$ sonic only; $p^*$, $T^*$ per process | `A_star`, `A_over_Astar`; `fanno_ratios` / `rayleigh_ratios` say "star of the same process" | "which star?" table (trap T2) |
| iii | stations | 1 upstream → 2 downstream of a shock or duct | 0 → 1 on the cited NASA pages; x and y in some texts | `M1`, `M2` | one line in the box; E3 Explain section 0 (trap T6) |
| iv | oblique shock | shock angle **σ**, deflection **δ** | **β**, **θ** | `sigma`, `delta` (docstring: "β, θ elsewhere") | "σ (= β elsewhere)", "δ (= θ elsewhere)" in every legend and readout of C14 and E7 (trap T13) |
| v | sign of a turn | δ > 0 when the wall turns **into** the flow (compression); $\delta_1+\nu(M_1)=\delta_2+\nu(M_2)$ (unnumbered, §15.9) | the Prandtl–Meyer turn counted positive for an **expansion** | `prandtl_meyer_turn(M1, turn)` with `turn` > 0 = expansion; docstring "the book's δ = −turn" | both numbers in every C15 / E8 legend: "turn +10° (book: δ = −10°)" (trap T14) |
| vi | friction and heat | lumped $f=F/(p_1A)$, $q=Q/h_1$ in $\rho_1u_1=\rho_2u_2,\ \ p_1+\rho_1u_1^2=p_2+\rho_2u_2^2+p_1f,\ \ h_1+\tfrac12u_1^2+h_1q=h_2+\tfrac12u_2^2$ (15.45, with $f=F/(p_1A)$, $q=Q/h_1$) — **primary** | Fanno flow with $4fL/D$ (f a friction factor); Rayleigh flow with $T_0/T_0^*$ | `f`, `q` everywhere; `fanno_ratios`, `rayleigh_ratios`, `fanno_length_to_f` labelled "not in the book — ours" | a two-column table in C12 (decision 4) |
| vii | the wall force | F on the fluid in +x in $-\dot mu_1+\dot mu_2=(pA)_1-(pA)_2+F$ (15.24); $F_f$ a positive retarding force per length in $\dot m\dfrac{du}{dx}=-A\dfrac{dp}{dx}-F_f$ (15.25); F the total retarding force in $f=F/(p_1A)$ | — | `wall_force` returns the first; `f` uses the third | three-row sign table in C03 and again in C11 (trap T11) |
| viii | entropy | S (per unit mass) | s (ch01 of this project) | `dS`, `dS_over_cv`, `dS_over_R` | one line in the box |
| ix | lapse rate | not used in the chapter's equations | — | — | the standard atmosphere behind the rocket example (N46) uses `core.statics`; where its lapse rate is printed both are: Kundu $\Gamma\equiv dT/dz=-6.5$ K/km and the meteorological $-dT/dz=+6.5$ K/km |

5. **Overloaded letters (trap T5) — the box shows this table; each block repeats its row.**

| Letter | Meanings in this chapter (and earlier) | Notebook / explainer symbol | Code name |
|---|---|---|---|
| f, g | the two wave shapes of $p'(x,t)=f(x-ct)+g(x+ct)$ (15.15); **f** also the friction parameter of $\rho_1u_1=\rho_2u_2,\ \ p_1+\rho_1u_1^2=p_2+\rho_2u_2^2+p_1f,\ \ h_1+\tfrac12u_1^2+h_1q=h_2+\tfrac12u_2^2$ (15.45, with $f=F/(p_1A)$, $q=Q/h_1$); g gravity elsewhere | f, g shapes in §15.2 only; f friction from §15.8 | `f`, `g` (callables) in `AC`; `f` (float) in `GD.duct_*` |
| q, Q | Q heat added per unit mass in $h_2+\tfrac12u_2^2-h_1-\tfrac12u_1^2=Q$ (15.22); $q=Q/h_1$ dimensionless in $\rho_1u_1=\rho_2u_2,\ \ p_1+\rho_1u_1^2=p_2+\rho_2u_2^2+p_1f,\ \ h_1+\tfrac12u_1^2+h_1q=h_2+\tfrac12u_2^2$ (15.45, with $f=F/(p_1A)$, $q=Q/h_1$); $q_j$ heat-flux vector in $\displaystyle\int_{A^*}\rho\left(e+\tfrac12u_i^2\right)u_jn_j\,dA=\int_{A^*}u_i\tau_{ij}n_j\,dA-\int_{A^*}q_jn_j\,dA$ (15.20) | as listed, named each time | `Q`, `q` |
| A | area; **the pure number** of $A\equiv\dfrac{M_1^2\left(1+\frac12(\gamma-1)M_1^2+q\right)}{\left(1+\gamma M_1^2-f\right)^2}=\dfrac{M_2^2\left(1+\frac12(\gamma-1)M_2^2\right)}{\left(1+\gamma M_2^2\right)^2}$ (unnumbered, §15.8) | A area; "the duct number A" | `A`, `A_over_Astar`, `duct_parameter_A` |
| δ | deflection angle (§15.9–§15.10); shock thickness (§15.6) | δ deflection; $\delta_{\text{sh}}$ thickness | `delta`, `thickness` |
| μ | Mach angle (§15.9); viscosity (§15.6, with $\mu''=\tfrac43\mu+\mu_v$) | μ Mach angle; $\mu_{\text{visc}}$ where both could meet | `mu` (angle), `mu_visc` |
| ν | Prandtl–Meyer function (§15.9); kinematic viscosity (§15.6) | ν(M); $\nu_{\text{visc}}$ | `nu`, `nu_visc` |
| σ | shock angle (ch07: surface tension; ch11: growth rate) | σ shock angle | `sigma` |
| U, V, v | U stream speed (§15.2); $U=u/V$ scaled speed inside a shock (§15.6); V flow speed, v shock-parallel component (§15.9) | as listed | `U`, `U_scaled`, `V`, `v_t` |
| b | chord in $L\cong\dfrac{2\alpha\gamma M_\infty^2p_\infty b}{\sqrt{M_\infty^2-1}},\ \ D\cong\dfrac{2\alpha^2\gamma M_\infty^2p_\infty b}{\sqrt{M_\infty^2-1}}$ (15.58) (ch14: the map constant; span is s) | b chord (said each time) | `b` |
| c, m, e, h, η | c sound speed (ch14: chord); m mass flux inside a shock; e internal energy; h enthalpy (and a height in $t_{\text{delay}}=\dfrac{h\sqrt{M^2-1}}{Mc}$ (ours)); η shock coordinate | as listed | `c`, `m`, `e`, `h`, `height`, `eta` |

6. **Book slips taught in corrected form** — house form `> ⚠️ **slip #k — the book prints** … **; the correct form is** …` where
   used; `ch15.book_slips()` (each row has `kind`).

| Slip | The book prints | Correct | Taught in | Planted wrong variant a test must fail |
|---|---|---|---|---|
| **#4** (physics) | $u\,du=-dp/\rho=c^2d\rho/\rho$ (15.32) | $u\,du=-dp/\rho=-c^2d\rho/\rho$ (15.32, third member with its minus sign restored) | C06 (N41, D12 step 3, a ⚠️ callout, a cell) | `GD.influence_coefficients(printed_15_32=True)` — density rises with speed |
| **#11** (physics) | $-\dfrac{dp}{p_1}=\dfrac{1-(\gamma-1)M^2}{1-M^2}df$ (15.47, first relation) | $-\dfrac{dp}{p_1}=\dfrac{1+(\gamma-1)M^2}{1-M^2}\,df,\ \ \dfrac{du}{u}=-\dfrac{d\rho}{\rho}=\dfrac{p_1}{p}\dfrac{df}{1-M^2}$ (15.47, numerator sign corrected) | C12 (N84, D22 steps 5 and 8, a ⚠️ callout, a cell) | `printed_15_47=True` — misses finite differences of $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ by 18 %, 49 % and a sign at M = 0.5, 0.9, 2 |
| **#12** (physics) | $\dfrac{u\,du}{h_1}=\dfrac{(1-\gamma)M^2}{1-M^2}dq$ (15.48, second relation) | $\dfrac{dT}{T_1}=\left(\dfrac{1-\gamma M^2}{1-M^2}\right)dq,\ \ \dfrac{u\,du}{h_1}=\dfrac{(\gamma-1)M^2}{1-M^2}\,dq$ (15.48, second relation's sign corrected) | C12 (N85, D23 step 8, a ⚠️ callout, a cell) | `printed_15_48=True` — heating would slow a subsonic flow |
| #1 | "monotonic" gases | monatomic | C01 (R07) | — |
| #2 | $P+p'$ and $\rho_o$ in (15.9), (15.11), (15.14), (15.16) | $p_0$, $\rho_0$ throughout | C01 (N12, D01) | — |
| #3 | "(15.6)" cited where $p'\cong c^2\rho'$ (15.5) is used | $p'\cong c^2\rho'$ | C01 (N13, D01 step 9) | — |
| #5 | constant-density remark as $dA/dx+du/dx=0$ | $\dfrac1A\dfrac{dA}{dx}+\dfrac1u\dfrac{du}{dx}=0$ | C06 (R15, units cell) | — |
| #6 | lower-case m for an exit Mach number; a temperature set equal to a ratio | our own notation in both worked examples | C06 (N46), C10 (N78) | private |
| #7 | the dimensional quadratic set "≡ (U₁ − U)(U − U₂)" before dividing by $V^2$ | divide first | C08 (N68) | — |
| #8 | a figure caption's Mach number that does not match its plateaus | our figure at our own Mach numbers | C08 (N70) | private |
| #9, #10 | misprinted table entries | generated tables | C05 (N37), C08 (N63) | private |
| #13 | "found from (15.47)" for the plate pressures | $\dfrac{p-p_\infty}{p_\infty}=\dfrac{\gamma M_\infty^2\,\delta}{\sqrt{M_\infty^2-1}}$ (15.57) | C16 (N109, D31 step 1) | — |
| #14 | "$M_2$ calculated from (15.45)" | $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46) | C11 (N80) | — |
| #15 | figure labels without a decimal point | our figure with our own q values | C11 (N82) | — |
| #16, #17 | exercise wording; an exercise's quoted loss | nothing published | S01 | private |

7. **Traps (true as printed, easy to misread) — each a `> ⚠️ Common confusion:` callout where it first bites** (`ch15.traps()`):
   T1 three subscript 0's (box; C01 N05, R07; C04 N29) · T2 three stars, an A that is not an area (box; C03 R10; C05 N35; C11 N80;
   C12) · T3 regime edges are conventions (C01 N02) · T4 rounded constants (C01 R06) · T5 reused letters (box; where each first
   appears) · T6 which Mach number is upstream; $M_2^2=\dfrac{(\gamma-1)M_1^2+2}{2\gamma M_1^2+1-\gamma}$ (15.40) is its own inverse (C07 N54; E3) · T7 stagnation quantities depend on the
   frame (C07 N64; E3) · T8 two inverses of the area ratio; $A^*$ grows across a shock (C05 N36; C08 N62; C10 D20, N78) · T9
   silent assumptions of the rocket example (C06 N46, N47) · T10 "choked" is not "cannot be raised" (C09 N74, F6; C12 N88; E4 notes) ·
   T11 sign and meaning of the wall force; lumped f (C03 N27; C11 N79; C12 D22) · T12 which root of $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46) (C11 N83, D21 step 13;
   E5 toggle) · T13 σ, δ against β, θ; "weak" and "strong" are names of roots (box; C14 N93, N97, N98; E7) · T14 sign of the turn
   (box; C15 N107; C14 D28; E8) · T15 b is the chord; which surface turns into the flow (C16 N110, D31) · T16 how small is "small"
   (C14 N102; C16; F9; E7; E9) · T17 adiabatic is not isentropic (C03 R13; C04 N33, D09; C07 N51) · T18 the two paradoxes belong to
   a regime (C12 N87; E5) · T19 the indefinite integral, the sign of u′ under g, "downstream" (C02 N16, N17, N19, D04) · T20
   $Pr''$, the arbitrary origin of η, the continuum limit (C08 N67, N69) · T21 the hydraulic-jump analogy is partial (C08 R16; E3)
   · T22 the crude blast constant (C07 R17).
8. **Colour code (curation §5; the same in figures, derivations, bars and explainers):** pressure = orange · velocity and Mach
   number = blue · density = teal · temperature = rose · stagnation quantities and anything conserved = amber (flat when
   conserved) · entropy = purple (accent) · references, ghosts and forbidden regions = muted. **Every ratio on screen says what it
   is a ratio of** ("$p_2/p_1$ downstream / upstream", "$p/p_0$ static / stagnation").
9. **Honest labels.** (a) **"Not in the book — ours"** opens every appearance of: the sonic ratios $\dfrac{T^*}{T_0}=\dfrac{2}{\gamma+1},\ \dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)},\ \dfrac{\rho^*}{\rho_0}=\left(\dfrac{2}{\gamma+1}\right)^{1/(\gamma-1)}$ and the choked mass
   flow (N35, N39); the inverse of the area ratio (N36); every formula of §15.7 (D19, D20); the stagnation-pressure ratio $\dfrac{p_{02}}{p_{01}}=\left[\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}\right]^{\gamma/(\gamma-1)}\left[\dfrac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^{1/(\gamma-1)}$
   (N62); the choking limits (N88); the maximum deflection and sonic locus (N98); the steepening formulas (N22); the Hugoniot
   curve (N58); the general thin section (N112); the boom delay (N92); the closed-form shock-structure roots (N68); the Fanno and
   Rayleigh names, their ratios and their T–s lines (C12). (b) **"Sketch, not a computation"** is in the title of the detached
   bow shock (N99) and of the jet outside a nozzle (N76, E4). (c) **"Crude"** is in the name and title of the blast constant (R17).
   (d) A Prandtl–Glauert curve appears only with the label "not in the book".
10. **Public-repo rule.** No book number is reproduced; Tables 15.1 and 15.2 are never shown — `GD.isentropic_table` and
    `GD.normal_shock_table` generate ours. **The notebook's worked numbers use the following input set of this design, written
    out literally in the imports cell as six small dicts** (`ch15.illustrative_inputs()` is displayed once beside them; if the
    module's set differs the builder says "both are ours" — C.6):
    **air** `air = GD.PerfectGas(1.4, 287.05)`; sea level 288.15 K, 101 325 Pa → c = 340.3 m/s, $\rho_0$ = 1.225 kg/m³,
    $\rho_0c$ = 416.9 kg/(m² s) · **nozzle** `noz = dict(p0=5.0e5, T0=400.0, A_t=1.0e-3, Ae_over_At=2.0)` · **rocket**
    `rocket = dict(thrust=5.0e4, p0=3.0e6, T0=3000.0, p_exit=101325.0)` (air's γ and R, stated as a simplification) · **duct**
    `duct = dict(M1=0.5, f=0.1, q=0.2)` · **wedge** `wedge = dict(M1=2.0, delta=np.deg2rad(10.0))` · **plate**
    `plate = dict(alpha=np.deg2rad(2.0), M_inf=2.0, b=1.0, p_inf=101325.0)`. Regime edges passed by keyword:
    `transonic=(0.8, 1.2)`, `hypersonic=5.0` (ours; trap T3).

---

## Part C — functions the builders will call (the implementer's contract)

Names and signatures: **`fluidpy/core/gasdynamics.py` (GD, 81 public callables), `fluidpy/core/acoustics.py` (AC, 29) and
`fluidpy/ch15_compressible_flow.py` (ch15, 35 of its own + 129 re-exported) as found on disk on 2026-10-13** — they appeared
while this design was being written (Parts A and B had been drafted against analysis §4 and curation §8–§9 and were then
reconciled). Every row below was **probed by calling the function** with this design's inputs; return keys are copied from
the returned dicts. "Eq." = the book equation implemented (written out in the docstring). **(+)** = not in analysis §4 (C.5).
Return shapes are part of the contract: scalars in → floats out; a `dict` has at least the keys listed. Outside its range a
function returns NaN or raises `ValueError` with the stated meaning. Units in brackets. Docstrings carried "Label: pending"
when probed (the verifier sets the labels).

### C.0 Reused (existing; called, not changed; all re-exported by `ch15`)
| Function (as on disk) | Used for |
|---|---|
| `SIM.mach_number(U, c=None, T=None, gamma=GAMMA_AIR, R=R_AIR)` | R01: $M\equiv U/c$ (4.111) |
| `TH.perfect_gas_sound_speed(T, gamma=GAMMA_AIR, R=R_AIR)` [m/s] | R07, R08: $c^2=\gamma p_0/\rho_0=\gamma RT_0$ (15.6) |
| `TH.sound_speed_from_eos(p_of_rho_s, rho, s=0.0, rel_step=1e-6)` | R02: $c^2\equiv(\partial p/\partial\rho)_s$ (1.19) |
| `TH.isentropic_ratios(p_over_p0, gamma=GAMMA_AIR)` | R05: $\dfrac{p_2}{p_1}=\left(\dfrac{\rho_2}{\rho_1}\right)^{\gamma},\ \dfrac{T_2}{T_1}=\left(\dfrac{\rho_2}{\rho_1}\right)^{\gamma-1}=\left(\dfrac{p_2}{p_1}\right)^{(\gamma-1)/\gamma}$ (15.2) |
| `TH.perfect_gas_entropy_change(T1, v1, T2, v2, cv=CV_AIR, R=R_AIR)`, `TH.perfect_gas_entropy_change_p(T1, p1, T2, p2, cp=CP_AIR, R=R_AIR)` [J/(kg K)] | R04, C08 from-scratch: $S_2-S_1=C_p\ln\dfrac{T_2}{T_1}-R\ln\dfrac{p_2}{p_1}=C_v\ln\dfrac{T_2}{T_1}-R\ln\dfrac{\rho_2}{\rho_1}$ (15.1, entropy line) |
| `TH.R_AIR` (= 287.058), `TH.GAMMA_AIR`, `TH.CP_AIR`, `TH.CV_AIR`, `TH.P_ATM` | R06 |
| `BER.stagnation_temperature(T, U, cp=CP_BOOK)` [K], `BER.stagnation_pressure(p, U, rho)` [Pa] (incompressible) | R14: $C_pT_0=C_pT+\tfrac12u^2$ (15.27); N32 |
| `ch14.compressibility_check(U, T=288.15)` → dict | R03 |
| `ch07.hydraulic_jump(H1, u1=None, Fr1=None, g=G)` → dict; `ch07.simple_wave_evolve(eta0, x, t, H, g=G)` → dict | R16: $H_2/H_1=\tfrac12\left(-1+\sqrt{1+8Fr_1^2}\right)$ (7.81); N22 (the γ = 2 twin) |
| `ch01.mean_free_path_air(T, p=P_ATM)` [m]; `ch01.blast_radius(E, t, rho, K=1.0, *, geometry="sphere")` [m] | N65; R17 |
| `AF.zhukhovsky_lift_coefficient(alpha, beta, e=0.0, exact=False)` | N111: $C_L\simeq2\pi(\alpha+\beta)$ (14.12) |
| `STAT` standard atmosphere (pressure at altitude) | N46, N47 (off-design thrust) |

### C.1 `fluidpy/core/gasdynamics.py` (`GD`; exists; re-exported by `ch15`)
| Function | Signature → return (keys as found) | Eq. |
|---|---|---|
| `PerfectGas`, `AIR_BOOK` | dataclass `PerfectGas(gamma=1.4, R=R_AIR)`; `.cp`, `.cv` [J/(kg K)]; `.c(T [K])` [m/s]. `AIR_BOOK = PerfectGas(1.4, 287.0)` | $C_v=\dfrac{R}{\gamma-1},\ C_p=\dfrac{\gamma R}{\gamma-1}$ (15.1, specific-heat line), $c=\sqrt{\gamma RT}=\sqrt{\gamma p/\rho}$ (15.1, sound-speed line) |
| `mach_from_speed` | `(u [m/s], T [K], gas=None)` → M | $M\equiv U/c$ (4.111) |
| `mass_flow`, `mass_flow_from_mach` | `mass_flow(rho, u, A [m²])` → kg/s; `mass_flow_from_mach(M, p0 [Pa], T0 [K], A, gas=None)` → kg/s | $\rho uA=\dot m=\text{const.}$ (15.18); $\dfrac{\dot m\sqrt{RT_0}}{p_0A}=\sqrt\gamma\,M\left(1+\dfrac{\gamma-1}{2}M^2\right)^{-(\gamma+1)/(2(\gamma-1))}$ (ours) |
| `total_enthalpy`, `energy_residual` | `total_enthalpy(T, u, cp)` → J/kg; `energy_residual(state1, state2, Q=0.0, gamma=1.4)` → J/kg; a state is a dict with `u` and either `h` or `p`, `rho` | $h_2+\tfrac12u_2^2-h_1-\tfrac12u_1^2=Q$ (15.22) |
| `wall_force` | `(state1, state2, A1, A2)` → N, on the fluid in +x; states are dicts with `rho`, `u`, `p` | $-\dot mu_1+\dot mu_2=(pA)_1-(pA)_2+F$ (15.24) |
| `influence_coefficients` | `(M, gamma=1.4, *, dA_over_A=0.0, df=0.0, dq=0.0, printed_15_32=False, printed_15_47=False, printed_15_48=False)` → dict `drho_over_rho`, `du_over_u`, `dp_over_p`, `dT_over_T`, `dM_over_M`, `ds_over_cp`, `u_du_over_h`, `dT0_over_T0`, `dp0_over_p0`; `df`, `dq` normalised with the **local** p and h | $\dfrac{d\rho}{\rho}+\dfrac{du}{u}+\dfrac{dA}{A}=0$ (15.19), $u\,du=-dp/\rho=-c^2d\rho/\rho$ (15.32, third member with its minus sign restored), $d\rho/\rho=-M^2\,du/u$ (15.33), $\dfrac{du}{u}=-\dfrac{1}{1-M^2}\dfrac{dA}{A}$ (15.34), $-\dfrac{dp}{p_1}=\dfrac{1+(\gamma-1)M^2}{1-M^2}\,df,\ \ \dfrac{du}{u}=-\dfrac{d\rho}{\rho}=\dfrac{p_1}{p}\dfrac{df}{1-M^2}$ (15.47, numerator sign corrected), $\dfrac{dT}{T_1}=\left(\dfrac{1-\gamma M^2}{1-M^2}\right)dq,\ \ \dfrac{u\,du}{h_1}=\dfrac{(\gamma-1)M^2}{1-M^2}\,dq$ (15.48, second relation's sign corrected) |
| `area_velocity`, `duct_role` | `area_velocity(M, dA_over_A)` → du/u (so `area_velocity(M, 1.0)` is the coefficient −1/(1 − M²)); `duct_role(M, dA_sign)` → `"nozzle"`, `"diffuser"`, `"none"` (dA = 0) or `"sonic"` (M = 1) | $\dfrac{du}{u}=-\dfrac{1}{1-M^2}\dfrac{dA}{A}$ (15.34) |
| `stagnation_temperature_ratio`, `stagnation_pressure_ratio`, `stagnation_density_ratio` | `(M, gamma=1.4)` → $T_0/T$, $p_0/p$, $\rho_0/\rho$ (≥ 1) | $\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28), $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29), $\dfrac{\rho_0}{\rho}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{1/(\gamma-1)}$ (15.30) |
| `mach_from_pressure_ratio` | `(p_over_p0, gamma=1.4)` → M | inverse of $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29) |
| `isentropic_state` | `(M, gamma=1.4)` → dict `M`, `T_over_T0`, `p_over_p0`, `rho_over_rho0`, `A_over_Astar`, `T0_over_T`, `p0_over_p`, `rho0_over_rho` | $\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28)–$\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31) |
| `isentropic_table` | `(M, gamma=1.4)` → DataFrame with columns `"M"`, `"p/p0"`, `"rho/rho0"`, `"T/T0"`, `"A/A*"` (M is required) | ours (replaces Table 15.1) |
| `sonic_ratios`, `critical_pressure_ratio` | `sonic_ratios(gamma=1.4)` → dict `T`, `p`, `rho`, `c` (each star over stagnation); `critical_pressure_ratio(gamma=1.4)` → float | $\dfrac{T^*}{T_0}=\dfrac{2}{\gamma+1},\ \dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)},\ \dfrac{\rho^*}{\rho_0}=\left(\dfrac{2}{\gamma+1}\right)^{1/(\gamma-1)}$ (ours, from (15.28)–(15.30) at $M=1$) |
| `max_speed` | `(T0 [K], cp [J/(kg K)])` → m/s | $u_{\max}=\sqrt{2C_pT_0}$ (ours) |
| `pitot_compressibility_factor`, `mach_for_density_change` | `(M, gamma=1.4)` → $(p_0-p)/(\tfrac12\rho u^2)$; `(fraction, gamma=1.4)` → M | $p_0=p+\tfrac12\rho u^2\left(1+\tfrac14M^2+\dots\right)$ (ours) |
| `entropy_change_from_stagnation`, `stagnation_pressure_loss` | `(p02_over_p01, T02_over_T01=1.0, gas=None)` → J/(kg K); `stagnation_pressure_loss(dS, gas=None)` → $p_{02}/p_{01}=e^{-\Delta S/R}$ | $S_2-S_1=-R\ln\dfrac{p_{02}}{p_{01}}+C_p\ln\dfrac{T_{02}}{T_{01}}$ (unnumbered, §15.4) |
| `area_ratio`, `mach_from_area_ratio` | `area_ratio(M, gamma=1.4)` → $A/A^*$; `mach_from_area_ratio(A_over_Astar, gamma=1.4, *, branch="subsonic")` (`"subsonic"` or `"supersonic"`; **always pass `branch=` in the notebook**); `ValueError` below 1 | $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31) |
| `mass_flux_function`, `choked_mass_flow` | `mass_flux_function(M, gamma=1.4)` → $\dot m\sqrt{RT_0}/(p_0A)$; `choked_mass_flow(p0, T0, A_star, gas=None)` → kg/s | $\dfrac{\dot m\sqrt{RT_0}}{p_0A}=\sqrt\gamma\,M\left(1+\dfrac{\gamma-1}{2}M^2\right)^{-(\gamma+1)/(2(\gamma-1))}$ (ours), $\dot m_{\max}=p_0A^*\sqrt{\dfrac{\gamma}{RT_0}}\left(\dfrac{2}{\gamma+1}\right)^{(\gamma+1)/(2(\gamma-1))}$ (ours) |
| `laval_area`, `nozzle_isentropic` | `laval_area(x, A_throat=1.0, A_inlet=2.5, A_exit=2.0, x_throat=0.4, L=1.0)` → area (two parabolas, C.6); `nozzle_isentropic(x, A_fn, gamma=1.4, *, branch="transonic", A_star=None)` (`"subsonic"`, `"supersonic"`, `"transonic"`) → dict of arrays `x`, `A`, `M`, `p_over_p0`, `T_over_T0`, `rho_over_rho0` + `A_star`, `i_throat`, `x_throat` (**the sonic area defaults to the smallest area on the grid — put `x_throat` on the grid**) | $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31) |
| `normal_shock` | `(M1, gamma=1.4, *, allow_expansion=False)` → dict `M1`, `M2`, `p_ratio`, `rho_ratio`, `u_ratio`, `T_ratio`, `p0_ratio`, `T0_ratio`, `dS_over_cv`, `dS_over_R`, `admissible` (ratios are downstream / upstream); `ValueError` for M1 < 1 unless `allow_expansion` | $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39)–$\dfrac{S_2-S_1}{C_v}=\ln\left[\dfrac{p_2}{p_1}\left(\dfrac{\rho_1}{\rho_2}\right)^{\gamma}\right]$ (15.43) |
| `normal_shock_pressure_ratio`, `normal_shock_mach`, `normal_shock_density_ratio`, `normal_shock_temperature_ratio` | `(M1, gamma=1.4)` → float | $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39), $M_2^2=\dfrac{(\gamma-1)M_1^2+2}{2\gamma M_1^2+1-\gamma}$ (15.40), $\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}$ (15.41), $\dfrac{T_2}{T_1}=1+\dfrac{2(\gamma-1)}{(\gamma+1)^2}\dfrac{\gamma M_1^2+1}{M_1^2}\left(M_1^2-1\right)$ (15.42) |
| `normal_shock_residuals`, `normal_shock_fluxes` (+) | `normal_shock_residuals(state1, state2, gamma=1.4)` → dict `mass`, `momentum`, `energy` (relative; states are dicts with `rho`, `u`, `p`); `normal_shock_fluxes(M1, gamma=1.4)` → dict `mass`, `momentum`, `energy` (each downstream / upstream = 1) and the shares `p1_part`, `ke1_part`, `p2_part`, `ke2_part` (of the momentum flux), `h1_part`, `k1_part`, `h2_part`, `k2_part` (of the total enthalpy) | $\rho_1u_1=\rho_2u_2$ (15.35), $p_1-p_2=-\rho_1u_1^2+\rho_2u_2^2$ (15.36), $h_1+\tfrac12u_1^2=h_2+\tfrac12u_2^2$ (15.37) |
| `normal_shock_entropy`, `weak_shock_entropy`, `shock_admissible` | `normal_shock_entropy(M1, gamma=1.4)` → $(S_2-S_1)/C_v$; `weak_shock_entropy(x, gamma=1.4, *, form="mach")` with `form="mach"` (x = $M_1$) or `"pressure"` (x = $(p_2-p_1)/p_1$); `shock_admissible(M1)` → bool | $\dfrac{S_2-S_1}{C_v}=\ln\left[\dfrac{p_2}{p_1}\left(\dfrac{\rho_1}{\rho_2}\right)^{\gamma}\right]$ (15.43), $\dfrac{S_2-S_1}{C_v}\cong\dfrac{2\gamma(\gamma-1)}{3(\gamma+1)^2}\left(M_1^2-1\right)^3$ (15.44a), $\dfrac{S_2-S_1}{C_v}\cong\dfrac{\gamma^2-1}{12\gamma^2}\left(\dfrac{p_2-p_1}{p_1}\right)^3$ (15.44b) |
| `normal_shock_stagnation_pressure_ratio` | `(M1, gamma=1.4)` → $p_{02}/p_{01}$ | $\dfrac{p_{02}}{p_{01}}=\left[\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}\right]^{\gamma/(\gamma-1)}\left[\dfrac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^{1/(\gamma-1)}$ (ours) |
| `normal_shock_table` | `(M1, gamma=1.4)` → DataFrame `"M1"`, `"M2"`, `"p2/p1"`, `"rho2/rho1"`, `"T2/T1"`, `"p02/p01"` (M1 required) | ours (replaces Table 15.2) |
| `hugoniot_pressure_ratio` | `(rho_ratio, gamma=1.4)` → $p_2/p_1$ | $\dfrac{p_2}{p_1}=\dfrac{(\gamma+1)\rho_2/\rho_1-(\gamma-1)}{(\gamma+1)-(\gamma-1)\rho_2/\rho_1}$ (ours) |
| `moving_shock`, `strong_shock_limits` | `moving_shock(p1 [Pa], T1 [K], *, M_s=None, p2_over_p1=None, gas=None)` → dict `M_s`, `u_shock`, `u_gas`, `p2`, `T2`, `rho1`, `rho2`, `c1`, `c2`, `M_gas`, `T0_behind_lab`, `p0_behind_lab`; `strong_shock_limits(gamma=1.4)` → dict `rho_ratio`, `M2`, `p_ratio_over_M1sq`, `T_ratio_over_p_ratio`, `u1sq_rho1_over_p2` | N64; R17 |
| `convergent_nozzle` | `(p0, T0, pB, A_exit, gas=None)` → dict `M_exit`, `p_exit`, `T_exit`, `u_exit`, `mdot`, `mdot_max`, `choked`, `matched` | D19 |
| `mass_flow_vs_back_pressure` (+) | `(pB_over_p0, gamma=1.4, *, dimensionless=True, p0=None, T0=None, A_exit=None, gas=None)` → dict `mdot` (dimensionless: $\dot m\sqrt{RT_0}/(p_0A_e)$), `M_exit`, `p_exit_over_p0`, `choked`, `mdot_over_choked` | $\dfrac{\dot m\sqrt{RT_0}}{p_0A}=\sqrt\gamma\,M\left(1+\dfrac{\gamma-1}{2}M^2\right)^{-(\gamma+1)/(2(\gamma-1))}$ (ours) with $\dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)}$ (ours) |
| `laval_critical_pressures`, `laval_regime` | `laval_critical_pressures(Ae_over_At, gamma=1.4)` → dict `p_b`, `p_d`, `p_f` (each over $p_0$), `M_sub`, `M_sup`, `p_star`; `laval_regime(pB_over_p0, Ae_over_At, gamma=1.4, *, rtol=1e-9)` → `"subsonic"`, `"sonic throat"` (at $p_b$), `"shock in nozzle"`, `"shock at exit"` (at $p_d$), `"overexpanded"`, `"design"` (at $p_f$), `"underexpanded"` (the three equalities need the exact thresholds) | $\dfrac{p_b}{p_0}=\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sub}}(A_e/A_t)},\ \dfrac{p_f}{p_0}=\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sup}}(A_e/A_t)},\ \dfrac{p_d}{p_0}=\dfrac{p_f}{p_0}\left.\dfrac{p_2}{p_1}\right\rvert_{M_{\text{sup}}}$ (ours) |
| `laval_shock_location`, `laval_exit_state` (+) | `laval_shock_location(pB_over_p0, Ae_over_At, gamma=1.4)` → dict `M1`, `M2`, `A_shock_over_At`, `p02_over_p01`, `M_exit`, `p_exit_over_p0`; `laval_exit_state(pB_over_p0, Ae_over_At, gamma=1.4)` → dict `regime`, `M1_shock`, `A_shock_over_At`, `p02_over_p01`, `M_throat`, `mdot_over_choked`, `choked`, `M_exit`, `p_exit_over_p0` | $\dfrac{p_B}{p_0}=\dfrac{p_{02}}{p_{01}}(M_1)\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sub}}\left(\frac{A_e}{A_t}\frac{p_{02}}{p_{01}}\right)}$ (ours) |
| `laval_nozzle` | `(x, A_fn, pB_over_p0, gamma=1.4)` → dict `x`, `A`, `M`, `p_over_p0`, `regime`, `x_shock`, `A_shock_over_At`, `M1_shock`, `p02_over_p01`, `M_exit`, `p_exit_over_p0`, `i_throat`, `x_throat`, `mdot_over_choked` (grid-dependent: use ≥ 400 stations including the throat) | D20 |
| `duct_parameter_A`, `duct_outlet_mach` | `duct_parameter_A(M1, f=0.0, q=0.0, gamma=1.4)` → A; `duct_outlet_mach(M1, f=0.0, q=0.0, gamma=1.4, *, branch="continuous")` with `branch` = `"continuous"` (the root that tends to $M_1$), `"shocked"` (the other), `"plus"`, `"minus"` (the two signs of the formula) or `"both"` (a pair); NaN when choked | $A\equiv\dfrac{M_1^2\left(1+\frac12(\gamma-1)M_1^2+q\right)}{\left(1+\gamma M_1^2-f\right)^2}=\dfrac{M_2^2\left(1+\frac12(\gamma-1)M_2^2\right)}{\left(1+\gamma M_2^2\right)^2}$ (unnumbered, §15.8), $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46) |
| `duct_state_ratios`, `duct_entropy_change`, `duct_residuals` | `duct_state_ratios(M1, M2, f=0.0, q=0.0, gamma=1.4)` → dict `p_ratio`, `T_ratio`, `rho_ratio`, `u_ratio`, `p0_ratio`, `T0_ratio`, `dS_over_cv`, `dS_over_cp` (outlet / inlet); `duct_entropy_change(M1, M2, f=0.0, q=0.0, gamma=1.4)` → $(S_2-S_1)/C_v$; `duct_residuals(state1, state2, f=0.0, q=0.0, gamma=1.4)` → dict `mass`, `momentum`, `energy` (states are dicts with `rho`, `u`, `p`) | $\rho_1u_1=\rho_2u_2,\ \ p_1+\rho_1u_1^2=p_2+\rho_2u_2^2+p_1f,\ \ h_1+\tfrac12u_1^2+h_1q=h_2+\tfrac12u_2^2$ (15.45, with $f=F/(p_1A)$, $q=Q/h_1$) |
| `duct_choking_limits`, `duct_curve` (+), `duct_choking_locus` (+) | `duct_choking_limits(M1, gamma=1.4)` → dict `f_max`, `q_max`; `duct_curve(M1, f=0.0, q=0.0, gamma=1.4)` → dict of arrays `M1`, `M2_continuous`, `M2_shocked`, `choked`, `dS_continuous`, `dS_shocked`; `duct_choking_locus(gamma=1.4, M1=None)` → dict `M1`, `f_max`, `q_max` | $f_{\max}=1+\gamma M_1^2-M_1\sqrt{2(\gamma+1)\left(1+\tfrac{\gamma-1}{2}M_1^2\right)}$ (ours), $q_{\max}=\dfrac{(1+\gamma M_1^2)^2}{2(\gamma+1)M_1^2}-1-\dfrac{\gamma-1}{2}M_1^2$ (ours) |
| `rayleigh_temperature_sign` | `(M, gamma=1.4)` → +1.0, 0.0 or −1.0 (the sign of dT/dq) | $\dfrac{dT}{T_1}=\left(\dfrac{1-\gamma M^2}{1-M^2}\right)dq,\ \ \dfrac{u\,du}{h_1}=\dfrac{(\gamma-1)M^2}{1-M^2}\,dq$ (15.48, second relation's sign corrected) |
| `fanno_ratios`, `rayleigh_ratios`, `fanno_length_to_f`, `fanno_line_Ts` (+), `rayleigh_line_Ts` (+) | `fanno_ratios(M, gamma=1.4)` → dict `T`, `p`, `rho`, `u`, `p0` (each over its sonic value), `fL` (= $4f_FL^*/D$); `rayleigh_ratios(M, gamma=1.4)` → dict `T`, `p`, `rho`, `u`, `T0`, `p0`; `fanno_length_to_f(M1, M2, gamma=1.4)` → dict `friction_length` (= $4f_FL/D$), `f` (the book's lumped f); `*_line_Ts(M, gamma=1.4)` → dict `M`, `ds_over_cp`, `T_over_Tstar` | not in the book — ours |
| `mach_angle`, `mach_from_angle` | `mach_angle(M)` → rad (NaN for M < 1); `mach_from_angle(mu)` → M | $\sin\mu=1/M$ (15.49) |
| `oblique_shock_components` | `(M1, sigma, gamma=1.4)` → dict `Mn1`, `Mt1`, `Mn2`, `M2`, `delta`, `u1_over_V1`, `v_over_V1`, `u2_over_V1`, `V2_over_V1`, `c2_over_c1` | $M_{n1}=M_1\sin\sigma,\ \ M_{n2}=M_2\sin(\sigma-\delta)$ (unnumbered, §15.9) |
| `oblique_shock` | `(M1, *, delta=None, sigma=None, branch="weak", gamma=1.4)` → dict `sigma`, `delta`, `M2`, `Mn1`, `Mn2`, `p_ratio`, `rho_ratio`, `T_ratio`, `p0_ratio`, `dS_over_cv`, `attached`, `branch`; **raises `ValueError("detached: …")` beyond $\delta_{\max}$** | $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2\sin^2\sigma-1\right)$ (15.50), $\dfrac{\rho_2}{\rho_1}=\dfrac{(\gamma+1)M_1^2\sin^2\sigma}{(\gamma-1)M_1^2\sin^2\sigma+2}=\dfrac{\tan\sigma}{\tan(\sigma-\delta)}$ (15.51), $M_2^2\sin^2(\sigma-\delta)=\dfrac{(\gamma-1)M_1^2\sin^2\sigma+2}{2\gamma M_1^2\sin^2\sigma+1-\gamma}$ (15.52), $\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}$ (15.53) |
| `deflection_angle`, `oblique_shock_curve` (+) | `deflection_angle(M1, sigma, gamma=1.4)` → δ [rad]; `oblique_shock_curve(M1, gamma=1.4, n=400)` → dict `sigma`, `delta`, `M2` (arrays), `i_max`, `i_sonic`, `delta_max`, `sigma_max`, `delta_sonic`, `sigma_sonic`, `mu` | $\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}$ (15.53) |
| `max_deflection`, `sonic_deflection`, `detachment_mach` (+) | `max_deflection(M1, gamma=1.4)` → dict `delta_max`, `sigma_max`; `sonic_deflection(M1, gamma=1.4)` → dict `delta_sonic`, `sigma_sonic`; `detachment_mach(delta, gamma=1.4)` → smallest $M_1$ with an attached shock | ours (from $\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}$, $M_2^2\sin^2(\sigma-\delta)=\dfrac{(\gamma-1)M_1^2\sin^2\sigma+2}{2\gamma M_1^2\sin^2\sigma+1-\gamma}$) |
| `weak_wave_pressure`, `weak_wave_error` | `weak_wave_pressure(M, delta [rad], gamma=1.4)` → $(p_2-p_1)/p_1$; `weak_wave_error(M, delta, gamma=1.4)` → dict `linear`, `exact`, `ratio` (exact / linear), `rel_error` (linear / exact − 1), `kind` (`"shock"` for δ > 0, `"expansion"` for δ < 0) | $\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta$ (15.55) |
| `prandtl_meyer`, `prandtl_meyer_max`, `mach_from_prandtl_meyer` | `prandtl_meyer(M, gamma=1.4)` → ν [rad]; `prandtl_meyer_max(gamma=1.4)` → rad; `mach_from_prandtl_meyer(nu, gamma=1.4)` → M | $\nu(M)=\sqrt{\dfrac{\gamma+1}{\gamma-1}}\tan^{-1}\sqrt{\dfrac{\gamma-1}{\gamma+1}\left(M^2-1\right)}-\tan^{-1}\sqrt{M^2-1}$ (15.56) |
| `prandtl_meyer_turn` | `(M1, turn [rad; + = expansion], gamma=1.4)` → dict `M2`, `nu1`, `nu2`, `mu1`, `mu2`, `p_ratio`, `T_ratio`, `rho_ratio`, `delta_book` (= −turn) | $\delta_1+\nu(M_1)=\delta_2+\nu(M_2)$ (unnumbered, §15.9) |
| `riemann_exact` (optional, B1 only) | `(left, right, gamma=1.4, xi=0.0)` → dict | ours |

### C.2 `fluidpy/core/acoustics.py` (`AC`; exists; re-exported by `ch15`)
| Function | Signature → return | Eq. |
|---|---|---|
| `isentropic_pressure_perturbation` | `(rho_prime [kg/m³], rho0, p0 [Pa], gamma=1.4, order=1)` with `order` = 1, 2 or `"exact"` → p′ [Pa] | Taylor expansion of p(ρ); $p'\cong c^2\rho'$ (15.5) |
| `acoustic_density`, `acoustic_temperature` | `acoustic_density(p_prime, c)` → ρ′; `acoustic_temperature(p_prime, p0, T0, gamma=1.4)` → T′ [K] | $p'\cong c^2\rho'$ (15.5) |
| `condensation`, `linearisation_error` | `condensation(p_prime [Pa], rho0, c)` → $p'/\rho_0c^2$; `linearisation_error(p_prime, rho0, c, gamma=1.4)` → ½(γ − 1) × condensation (the size of the dropped term relative to the kept one; positive) | $p'/\rho_0c^2\ll1$ (15.7) |
| `newton_sound_speed`, `acoustic_impedance` | `newton_sound_speed(T, R)` → $\sqrt{RT}$ (R required); `acoustic_impedance(rho0, c)` → $\rho_0c$ | N10; N17 |
| `sound_pressure_level` (+) | `(p_rms [Pa], p_ref=2e-5)` → dB | primer P356 |
| `convected_wave_residual`, `wave_equation_residual` | `convected_wave_residual(p_fn(x, t), x, t, c, U=0.0, h=None)`, `wave_equation_residual(p_fn, x, t, c, h=None)` → left side of the wave equation by central differences | $\dfrac1{c^2}\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)^2p'-\dfrac{\partial^2p'}{\partial x_j\partial x_j}=0$ (15.12), $\dfrac1{c^2}\dfrac{\partial^2p'}{\partial t^2}-\dfrac{\partial^2p'}{\partial x^2}=0$ (15.13) |
| `pulse` (+), `PULSE_KINDS` | `pulse(kind="gaussian", width=1.0, amplitude=1.0, x0=0.0)` → callable s ↦ shape; kinds `"gaussian"` ($a\,e^{-\xi^2}$, ξ = (s − x0)/width), `"tophat"`, `"packet"`, `"sine"` | — |
| `dalembert`, `plane_wave_velocity` | `dalembert(x, t, c, f=None, g=None)` → p′; `plane_wave_velocity(x, t, c, rho0, f=None, g=None)` → u′ | $p'(x,t)=f(x-ct)+g(x+ct)$ (15.15), $u_1'(x,t)=\dfrac1{\rho_0c}\left(f(x-ct)-g(x+ct)\right)$ (15.16) |
| `dalembert_fg` (+), `dalembert_initial_value`, `velocity_from_pressure` | `dalembert_fg(F, G, c, x_ref=0.0)` → the pair of callables (f, g); `dalembert_initial_value(x, t, c, F, G=None)` → p′; `velocity_from_pressure(p_fn, x, t, rho0, *, t0=0.0, h=1e-4, limit=200)` → u′ by quadrature from a time t0 before the wave | $f(x)=\dfrac12\left[F(x)-\dfrac1c\displaystyle\int_0^xG(\hat x)\,d\hat x\right],\ g(x)=\dfrac12\left[F(x)+\dfrac1c\displaystyle\int_0^xG(\hat x)\,d\hat x\right]$ (unnumbered, Exercise 15.1); $u_1'(x,t)=-\dfrac1{\rho_0}\displaystyle\int\dfrac{\partial p'}{\partial x}\,dt$ (15.14) |
| `dalembert_in_stream`, `wave_speeds`, `upstream_signal_possible` | `dalembert_in_stream(x, t, c, U, f=None, g=None)` → p′; `wave_speeds(U, c)` → (U + c, U − c); `upstream_signal_possible(U, c)` → bool | $p'(x,t)=f(x-(c+U)t)+g(x+(c-U)t)$ (15.17) |
| `xt_characteristics` (+) | `(U, c, x0, t)` → dict `t`, `x_plus`, `x_minus` (arrays shaped (len(x0), len(t))), `speed_plus`, `speed_minus` | N20 |
| `simple_wave_speed`, `simple_wave_evolve`, `breaking_time` | `simple_wave_speed(u, c0, gamma=1.4)`; `simple_wave_evolve(u0, x, t, c0, gamma=1.4)` → dict `x_points`, `x`, `u`, `c`, `speed`, `p_over_p0`, `rho_over_rho0`, `t_break`, `valid`; `breaking_time(u0, x, gamma=1.4, *, max_compression_rate=None)` → s | $\left.\dfrac{dx}{dt}\right\rvert_{\text{wavelet}}=u+c=c_0+\dfrac{\gamma+1}{2}u,\ \ t_{\text{break}}=\dfrac{2}{(\gamma+1)\max(-\partial u_0/\partial x)}$ (ours) |
| `source_wavefronts` | `(U, c, t, n_fronts=6, dt=1.0)` → dict `x_centres`, `y_centres`, `radii`, `t_emit` (most recent first), `x_source`, `M`, `mu`, `ahead` | N91 |
| `mach_cone_delay`, `in_zone_of_action` | `mach_cone_delay(h, M, c)` → s; `in_zone_of_action(x, y, M)` → bool — (x, y) relative to the source with **x measured downstream** (the frame in which the source is at rest and the stream flows in +x): inside ⇔ $x\ge\lvert y\rvert\sqrt{M^2-1}$; True everywhere for M < 1 | $t_{\text{delay}}=\dfrac{h\sqrt{M^2-1}}{Mc}$ (ours), $\sin\mu=1/M$ (15.49) |
| `acoustic_energy` (+), `doppler_factor` (+) | `acoustic_energy(p_prime, u_prime, x, rho0, c)` → J/m²; `doppler_factor(U, c, theta=0.0)` — optional one-line asides in C02 and C13 | ours |
| `wave_equation_step` (optional, B1 only) | `(p_now, p_prev, C, bc="periodic")` | ours |

### C.3 `fluidpy/ch15_compressible_flow.py` (`ch15`; exists; re-exports every public name of C.1, C.2 and the reused names of C.0)
| Function | Signature → return | Eq. |
|---|---|---|
| `illustrative_inputs`, `conventions_table`, `book_slips`, `traps` | `()` → flat dict (the module's own set — C.6) · DataFrame (9 rows; columns `item`, `here`, `elsewhere`, `code`) · DataFrame (17 rows, with `kind` and `changes_physics`) · DataFrame (22 rows) | — |
| `mach_regime` | `(M, *, transonic=(0.8, 1.2), hypersonic=5.0)` → `"incompressible"` (M = 0), `"subsonic"`, `"transonic"`, `"supersonic"`, `"hypersonic"` | N02 |
| `linearised_acoustics_sympy` | `()` → dict `symbols`, `continuity_linear`, `continuity_dropped`, `momentum_linear`, `momentum_dropped`, `eq_15_10`, `eq_15_11`, `wave_operator`, `cross_residual`, `dalembert_residual` | D01, D02 |
| `duct_equations_sympy` | `()` → dict `eq_15_19`, `eq_15_25`, `eq_15_26`, `eq_15_34`, `eq_15_33_printed` | D07, D12 |
| `normal_shock_sympy`, `weak_shock_series` | `normal_shock_sympy()` → dict `rho_ratio`, `p_ratio`, `T_ratio`, `M2_squared`, `roots`, `residual_15_39` … `residual_15_42`; `weak_shock_series(order=4)` → dict `m`, `series`, `coefficients`, `residual_15_44a`, `pressure_coefficient`, `residual_15_44b` | D13; D17 |
| `shock_structure_sympy`, `tangential_momentum_sympy`, `weak_wave_sympy`, `prandtl_meyer_sympy` | `()` → dicts: `momentum_integral`, `energy_integral`, `becker_form`, `structure_equation` · `v2`, `residual` · `residual_15_54`, `slope`, `residual_15_55` · `integrand`, `closed_form`, `derivative_residual`, `residual_numeric`, `value_at_1` | R19; D26; D28; D30 |
| `rocket_nozzle`, `nozzle_thrust` | `rocket_nozzle(thrust [N], p0, T0, p_exit, gas=None)` → dict `M_e`, `T_e`, `rho_e`, `u_e`, `A_e`, `A_throat`, `Ae_over_At`, `mdot`, `mdot_choked`; `nozzle_thrust(mdot, u_e, p_e, p_ambient, A_e)` → N | N46, N47 |
| `nozzle_with_shock` | `(M1, A_throat, A_exit, p0, T0, gas=None)` — $M_1$ is the Mach number just ahead of the shock — → dict `A_shock`, `M2`, `p02`, `A_star2`, `M3` (exit), `T3`, `p3`, `mdot`, `mdot_exit` | N78 |
| `blast_constant_crude`, `shock_thickness_estimate` | `blast_constant_crude(gamma=1.4, geometry="hemisphere")` → K (`"sphere"` also); `shock_thickness_estimate(du [m/s], nu [m²/s])` → m | R17; N65 |
| `shock_structure_constants`, `becker_shock_roots`, `becker_shock_ode`, `becker_shock_profile`, `becker_shock_thickness` | `shock_structure_constants(M1, T1, p1, gas=None)` → dict `m`, `V`, `I`, `u1`, `u2`, `U1`, `U2`; `becker_shock_roots(M1, gamma=1.4)` → (U1, U2); `becker_shock_ode(eta, U, M1, gamma=1.4)` → dU/dη; `becker_shock_profile(eta, M1, gamma=1.4, *, eta0=0.0)` → U; `becker_shock_thickness(M1, gamma=1.4)` → Δη (maximum-slope thickness; ours) | N66–N69: $\dfrac{U_1\ln(U_1-U)-U_2\ln(U-U_2)}{U_1-U_2}=\dfrac{\gamma+1}{2\gamma}\,\eta$ (unnumbered, §15.6) |
| `wedge_flow` | `(M1, half_angle, gamma=1.4)` → dict `attached`, `sigma`, `M2`, `p_ratio`, `T_ratio`, `rho_ratio`, `p0_ratio`, `subsonic_behind`, `delta_max`, `M_detach`, `axis_p_ratio` (numbers only — the streamlines are drawn by `draw_wedge_streamlines`, C.4) | N99 |
| `expansion_fan`, `fan_staircase` (+), `curved_wall_mach_lines` | `expansion_fan(M1, turn, gamma=1.4, n_rays=7, *, r_start=1.0)` → dict `theta`, `M`, `mu`, `ray_angle`, `p_over_p1`, `T_over_T1` (per ray), `M2`, `p2_over_p1`, `streamline_r`, `streamline_x`, `streamline_y`, `wall_angle2`; `fan_staircase(M1, turn, N=10, gamma=1.4)` → dict `theta`, `M`, `p_over_p1` (arrays of length N + 1), `M2_exact`, `p2_exact`, `M_error`, `p_error` (relative); `curved_wall_mach_lines(M1, wall_angle_fn, x, gamma=1.4)` → dict `x_wall`, `y_wall`, `theta`, `M`, `mu`, `line_angle`, `crossings`, `first_crossing`, `compression` | N104, N105; $d\delta=-\dfrac{\sqrt{M^2-1}}{M}\dfrac{dM}{1+\frac12(\gamma-1)M^2}$ (unnumbered, §15.9) |
| `ackeret_pressure`, `ackeret_cp` | `ackeret_pressure(delta, M_inf, gamma=1.4)` → $(p-p_\infty)/p_\infty$; `ackeret_cp(delta, M_inf)` → $2\delta/\sqrt{M_\infty^2-1}$ (**the angle comes first**) | $\dfrac{p-p_\infty}{p_\infty}=\dfrac{\gamma M_\infty^2\,\delta}{\sqrt{M_\infty^2-1}}$ (15.57) |
| `ackeret_flat_plate`, `flat_plate_supersonic` | `ackeret_flat_plate(alpha, M_inf)` → dict `CL`, `CD`, `L_over_D`, `lift_slope`; `flat_plate_supersonic(alpha, M_inf, gamma=1.4, *, method="linear")` with `method` = `"linear"` or `"shock_expansion"` (**always pass `method=`**) → dict `p_upper` (= $p_2/p_\infty$), `p_lower` (= $p_3/p_\infty$), `CN`, `CL`, `CD`, `L_over_pb`, `D_over_pb`, `method`; **raises `ValueError("detached: …")`** when the lower shock cannot stay attached | $C_L\cong\dfrac{4\alpha}{\sqrt{M_\infty^2-1}},\ \ C_D\cong\dfrac{4\alpha^2}{\sqrt{M_\infty^2-1}}$ (15.59), $L\cong\dfrac{2\alpha\gamma M_\infty^2p_\infty b}{\sqrt{M_\infty^2-1}},\ \ D\cong\dfrac{2\alpha^2\gamma M_\infty^2p_\infty b}{\sqrt{M_\infty^2-1}}$ (15.58) |
| `ackeret_airfoil` | `(yu_fn, yl_fn, alpha, M_inf, n=2001)` → dict `x`, `cp_upper`, `cp_lower`, `CL`, `CD`, `CD_incidence`, `CD_camber`, `CD_thickness`, `CD_cross` | $C_p=\pm\dfrac{2\theta}{\sqrt{M_\infty^2-1}},\ \ C_D=\dfrac{4}{\sqrt{M_\infty^2-1}}\left[\alpha^2+\overline{\theta_c^2}+\overline{\theta_t^2}\right]$ (ours) |
| `flat_plate_wave_pattern` (+), `lift_slope_vs_mach` (+) | `flat_plate_wave_pattern(alpha, M_inf, gamma=1.4)` → dict `le_fan_first`, `le_fan_last`, `le_shock`, `te_shock`, `te_fan_first`, `te_fan_last` (angles from the free stream), `M2`, `M3`, `M4_upper`, `M4_lower`, `p2`, `p3`, `p4_upper`, `p4_lower` (over $p_\infty$), `plate_angle`; `lift_slope_vs_mach(M, *, model="book", band=(0.8, 1.2))` → 2π below the band, $4/\sqrt{M^2-1}$ above, NaN inside | $C_L\cong\dfrac{4\alpha}{\sqrt{M_\infty^2-1}},\ \ C_D\cong\dfrac{4\alpha^2}{\sqrt{M_\infty^2-1}}$ (15.59), $C_L\simeq2\pi(\alpha+\beta)$ (14.12) |
| `euler1d_run` (optional, B1 only) | `(rho0, u0, p0, x, t_end, gamma=1.4, *, cfl=0.5, eps_av=2.0, bc="periodic")` → dict | ours (N114) |

### C.4 `scripts/`
`scripts/ch15_drawings.py` (**requested — not on disk when probed**; imported by the notebook as `from ch15_drawings import *`):
matplotlib helpers that only draw — `draw_duct_control_volume()`, `draw_reservoir_nozzle_tank()`, `draw_four_ducts()`,
`draw_shock_two_frames(ms)`, `draw_shock_triangles(M1, sigma, delta)`, `draw_wedge_streamlines(M1, half_angle, sigma)`,
`draw_mach_triangle(M)`, `draw_plate_waves(pattern)`, `draw_bow_shock_sketch()`, `draw_jet_sketch(kind)`, `draw_fan(fan)`,
`draw_wavefront_panels(fronts_list)`. Runnable demos (the implementer's choice of names): nozzle regimes, generated tables,
duct maps, the oblique-shock chart.

### C.5 Flag list — additions to analysis §4, and differences between the first draft of this design and the code as found
- **Not in analysis §4 and present on disk** (curation §8): `AC.pulse`, `AC.xt_characteristics`; `GD.mass_flow_vs_back_pressure`;
  `GD.laval_exit_state`; `GD.fanno_line_Ts`, `GD.rayleigh_line_Ts`, `GD.duct_curve`, `GD.duct_choking_locus`;
  `GD.normal_shock_fluxes`; `GD.oblique_shock_curve`, `GD.detachment_mach`; `ch15.fan_staircase`; `ch15.flat_plate_wave_pattern`,
  `ch15.lift_slope_vs_mach`. **Extra on disk, used here:** `AC.sound_pressure_level`, `AC.dalembert_fg`; offered as asides:
  `AC.acoustic_energy`, `AC.doppler_factor`.
- **Requested by this design, not on disk:** `scripts/ch15_drawings.py` (C.4) — nothing else.
- **Where the code differs from this design's first draft (the design now follows the code; listed so that the orchestrator
  can see what was reconciled — 24 items):** (1) dict keys are `p_ratio`, `rho_ratio`, `T_ratio`, `p0_ratio` (not
  `p2_over_p1` …) in `normal_shock`, `oblique_shock`, `duct_state_ratios`, `prandtl_meyer_turn`, `wedge_flow`; (2)
  `sonic_ratios` keys `T`, `p`, `rho`, `c`; (3) table columns `"p/p0"`, `"A/A*"`, `"p2/p1"` …; `isentropic_table` and
  `normal_shock_table` need their Mach array; (4) `area_velocity(M, dA_over_A)` returns du/u; (5) `duct_role` returns `"none"`
  at dA = 0 and `"sonic"` at M = 1; (6) `total_enthalpy(T, u, cp)`; (7) `energy_residual`, `normal_shock_residuals`,
  `duct_residuals` take two state dicts; (8) `laval_area` is two parabolas with `A_inlet=2.5`; (9) `nozzle_isentropic` takes its
  sonic area from the grid minimum; (10) `mass_flow_vs_back_pressure` returns `mdot` (dimensionless), not `flux`; (11)
  `laval_regime` strings `"sonic throat"`, `"shock in nozzle"`, `"shock at exit"`; (12) `duct_outlet_mach(branch="shocked")`; (13)
  `duct_entropy_change` returns ΔS/$C_v$; (14) `duct_curve`, `duct_choking_locus`, `fanno_line_Ts`, `rayleigh_line_Ts`,
  `weak_wave_error`, `max_deflection`, `sonic_deflection` return dicts with the keys above; (15) `fanno_ratios` / `rayleigh_ratios`
  keys `T`, `p`, `p0`, `fL`, `T0`; `fanno_length_to_f(M1, M2)`; (16) `oblique_shock` and `flat_plate_supersonic` **raise** when
  detached; (17) `AC.linearisation_error(p_prime, rho0, c)` returns the positive estimate; (18) `AC.newton_sound_speed(T, R)`;
  (19) `AC.pulse("gaussian")`; (20) `AC.source_wavefronts` keys `x_centres`, `radii`; `AC.in_zone_of_action(x, y, M)`; (21)
  `ch15.nozzle_with_shock(M1, A_throat, A_exit, p0, T0, gas)`; (22) `ch15.shock_thickness_estimate(du, nu)`,
  `ch15.shock_structure_constants(M1, T1, p1, gas)`, `ch15.blast_constant_crude(geometry=)`; (23) `ch15.ackeret_pressure(delta,
  M_inf)` and `ch15.ackeret_cp(delta, M_inf)` — angle first; (24) `ch15.wedge_flow` returns numbers, no streamlines;
  `ch15.mach_regime(0.1)` is `"subsonic"` ("incompressible" only at M = 0).
- **Two requests to the implementer (small; neither blocks a builder):** (a) make `branch=` of `mach_from_area_ratio` and
  `method=` of `flat_plate_supersonic` required keywords (convention-carrying arguments, the ch12–ch13 rule) — until then the
  notebook passes them explicitly in every call; (b) say in the docstring of `nozzle_isentropic` and `laval_nozzle` that the
  throat must be a grid point.

### C.6 Parity conventions (the explainers reproduce these exactly)
- Root finding in JS: `Viz.num.brentq(f, a, b, 1e-13)`; brackets as in Python: subsonic area inverse [1e-8, 1], supersonic
  [1, 100]; weak shock (μ, σ at $\delta_{\max}$), strong (σ at $\delta_{\max}$, π/2); Prandtl–Meyer inverse [1, 1e3]; the
  shock in a nozzle: $M_1$ in (1, $M_{\text{sup}}(A_e/A_t)$).
- `max_deflection`: the closed form (the quadratic in $\sin^2\sigma$) in both languages — never a minimiser in JS.
- Parity rows use scalar closed forms at rtol 1e-12, bracketed roots at 1e-9, nested roots at 1e-7. Dict results are indexed
  in the `py:` string (`ch15.oblique_shock(2.0, delta=np.deg2rad(10.0))["sigma"]`).
- **Definitions the JS mirrors copy exactly** (read from the code; probed): (a) `laval_area(x)`: for $x\le x_t$,
  $A=A_t+(A_{in}-A_t)\left[(x_t-x)/x_t\right]^2$; for $x\ge x_t$, $A=A_t+(A_e-A_t)\left[(x-x_t)/(L-x_t)\right]^2$ (defaults 1, 2.5, 2,
  $x_t$ = 0.4, L = 1; `laval_area(0.7)` = 1.25). (b) `source_wavefronts(U, c, t, n, dt)`: front k = 0 … n − 1 is emitted at
  `t_emit[k] = t - (k + 1)*dt` from `x_centres[k] = U*t_emit[k]`, with `radii[k] = c*(t - t_emit[k])`; `x_source = U*t`. (c)
  `fan_staircase(M1, turn, N)`: N explicit steps `M += (turn/N)*M*(1 + (gamma-1)/2*M*M)/sqrt(M*M - 1)` ($d\delta=-\dfrac{\sqrt{M^2-1}}{M}\dfrac{dM}{1+\frac12(\gamma-1)M^2}$ with turn = −δ);
  `p_over_p1` from $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29) at constant $p_0$. (d) `flat_plate_supersonic(method="shock_expansion")`: lower surface = weak
  oblique shock through δ = α, upper = Prandtl–Meyer turn α; $C_N=(p_3-p_2)/(\tfrac12\gamma M_\infty^2p_\infty)$,
  $C_L=C_N\cos\alpha$, $C_D=C_N\sin\alpha$; `method="linear"` uses $\dfrac{p-p_\infty}{p_\infty}=\dfrac{\gamma M_\infty^2\,\delta}{\sqrt{M_\infty^2-1}}$ (15.57) for both pressures with $C_L=C_N$,
  $C_D=C_N\alpha$ (so it equals $C_L\cong\dfrac{4\alpha}{\sqrt{M_\infty^2-1}},\ \ C_D\cong\dfrac{4\alpha^2}{\sqrt{M_\infty^2-1}}$ (15.59) exactly). (e) `duct_outlet_mach(branch="continuous")` is the + sign of $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$
  (15.46) for $M_1<1$ and the − sign for $M_1>1$. (f) `weak_shock_entropy(x, form="mach")` takes x = $M_1$. (g) `pulse("gaussian",
  w, a)`: $a\,e^{-(s/w)^2}$. (h) Fanno line: $T/T^*=(\gamma+1)/(2+(\gamma-1)M^2)$,
  $\Delta s/C_p=\tfrac{\gamma-1}{\gamma}\ln M-\tfrac{\gamma+1}{2\gamma}\ln\dfrac{2+(\gamma-1)M^2}{\gamma+1}$; Rayleigh line:
  $T/T^*=\left(\dfrac{(\gamma+1)M}{1+\gamma M^2}\right)^2$, $\Delta s/C_p=\ln M^2+\tfrac{\gamma+1}{\gamma}\ln\dfrac{\gamma+1}{1+\gamma M^2}$
  (ours; both agree with the code at M = 0.5: −0.08359, 1.14286 and −0.39996, 0.79012).
- `ch15.illustrative_inputs()` returns the module's own set (for example $T_0$ = 500 K, a 100 cm² throat, a rocket gas with
  γ = 1.25, a 12° wedge at Mach 2.5, a plate at 5°), different from and as much ours as the six literal dicts of header
  convention 10. **The notebook's *expect* numbers belong to the six dicts**; the module's set is displayed once beside them
  with the sentence "both sets are ours; neither is the book's".

### C.7 Test-only and notebook-only functions
Planted wrong variants a test must fail: `influence_coefficients(printed_15_32=True | printed_15_47=True | printed_15_48=True)`;
`duct_outlet_mach(branch="plus")` used on both sides of $M_1=1$ (trap T12). Private book values
(`tests/book_values_ch15.json`): both worked examples, the regime edges, the shock-thickness numbers, the table misprints, the
figure constants. Notebook-only: the from-scratch cells of curation §7 (sixteen `nb.check_agree` rows in Part A).

---

## Part A — notebook storyboard (`notebooks/build_ch15.py` → `notebooks/ch15_compressible_flow.ipynb`)

**One line per book section** (cell numbers are estimates; ≈ 520 cells in all):
- §15.1 → N01 N02 · N03 · R01–R06 (open the C01 story; the conventions box) — cells ≈ 10–32
- §15.2 → C01 (N04–N14 · R07 R08 · D01 D02 · P354 P355 P356 · F1), C02 (N15–N22 · N114 · D03 D04 D05 · P357 P358 · A1 A2 · **E1**) — cells ≈ 33–110
- §15.3 → C03 (N23–N28 · R09–R13 · D06 D07) — cells ≈ 111–150
- §15.4 → C04 (N29–N33 · N34 N38 · R14 · D08 D09 · P359 · F2), C05 (N35 N36 N37 N39 · D10 D11 · P360 P361 · F3) — cells ≈ 151–215
- §15.5 → C06 (N40–N47 · R15 · D12 · P362 · F4 · **E2**) — cells ≈ 216–255
- §15.6 → C07 (N48–N57 N64 · R17 · D13 D14 D15 · P363 · F5a), C08 (N58–N63 N65–N70 · N71 · R16 R18 R19 · D16 D17 D18 · P364 · F5b · **E3**) — cells ≈ 256–345
- §15.7 → C09 (N72 N73 N74 · D19 · F6), C10 (N75 N77 N78 · N76 · D20 · P365 · A3 · live · **E4**) — cells ≈ 346–395
- §15.8 → C11 (N79–N83 · D21 · P366 · F7), C12 (N84–N88 · N89 N90 · D22 D23 D24 · P367 · A5 · **E5**) — cells ≈ 396–450
- §15.9 → C13 (N91 N92 · D25 · P368 P369 · A4 · **E6**), C14 (N93–N99 N101 N102 · N100 · D26 D27 D28 · P370 P371 P372 · F8 · live · **E7**), C15 (N103–N107 · D29 D30 · P373 · F9b · **E8**) — cells ≈ 451–510
- §15.10 → C16 (N108–N112 · N113 · D31 · F9a · **E9**); S01, S02; summary — cells ≈ 511–535

Every CORE block follows: problem in plain words → idea → primers → maths (notes and derivations, Part F) → tiny example →
code (fluidpy) + "What does the code above do?" → from-scratch check → visual(s) → notes and "What would change if…". Code
drafts give intent + exact calls; the builder comments every line (novice grade, units, the equation written out next to its
number). *expect* = numbers the executed cell must print (hypotheses to be reported against, not edited). *see / read /
change* = the three figure notes (the builder rewrites them from the rendered figure). Note ids open every note in bold
(**N09 [B]**). Wherever a draft names an equation, the equation is written out with its number, and the builder keeps it so.
The six literal dicts of header convention 10 supply every worked number.

**Re-entry rule (binding; `tools/nbkit.py`).** `nb.recap(...)` and `nb.section(...)` close the current CORE block. The builder
writes `nb.current_core = "Cnn"` on the line after every `nb.recap` that sits inside block Cnn; rows marked **(re-enter Cnn)**
say where. Reminders of earlier-chapter ideas that are not inventory rows of this chapter (the first law and the Gibbs
relation, index notation and the Laplacian, control volumes, basic state + disturbance, the Knudsen number, the pressure
coefficient, d'Alembert's paradox, corner flows, the Kutta condition) are `nb.md` cells in the 🔁 recap format without an id;
they do not close the block.

### A.0 Front matter
1. `nb.title(big_idea=…, roadmap=[…16…], prerequisites=[…])`. **Big idea (draft):** "Until now the fluid could not be squeezed.
   Let it be squeezed and one new speed enters — the speed at which the fluid passes on the news that something has changed,
   the speed of sound c. Everything in this chapter is three conservation laws (mass, momentum, energy) and the gas law,
   written between two stations in terms of one number, the Mach number $M\equiv U/c$ (4.111). In every answer the same factor $1-M^2$
   appears, and it changes sign at M = 1. That one sign change is why a widening duct speeds up a supersonic stream, why
   friction speeds up a subsonic duct flow, why adding heat can cool a gas, why flows jump through shocks, and why a
   supersonic aircraft is heard only after it has passed." **Roadmap (one line per CORE):** C01 sound as the linear limit ·
   C02 two shapes sliding apart · C03 the three duct laws · C04 the stagnation ruler · C05 one area, two Mach numbers · C06
   the area–velocity surprise · C07 the normal-shock jump · C08 the second law picks the direction · C09 choking · C10 the
   Laval nozzle · C11 friction and heat in one formula · C12 everything drifts to M = 1 · C13 the Mach cone · C14 the oblique
   shock · C15 the Prandtl–Meyer fan · C16 lift and wave drag of a flat plate. **Prerequisites:** the perfect gas, entropy and
   the speed of sound $c^2\equiv(\partial p/\partial\rho)_s$ (1.19) (Ch. 1); the control-volume laws and the Mach number $M\equiv U/c$ (4.111) (Ch. 4); the hydraulic jump
   $H_2/H_1=\tfrac12\left(-1+\sqrt{1+8Fr_1^2}\right)$ (7.81) (Ch. 7); the lift slope of $C_L\simeq2\pi(\alpha+\beta)$ (14.12) (Ch. 14).
2. `nb.explainer_index([...])` — 9 rows: ("sound_pulse_dalembert", "What is a sound pulse made of?", "two shapes sliding apart
   at ±c — and both swept downstream when the stream beats c") · ("area_mach_nozzle", "Why does a supersonic nozzle widen?",
   "above Mach 1 the gas thins faster than it speeds up") · ("normal_shock_jump", "Why can a shock only compress?", "three laws
   balance both ways; entropy allows one") · ("nozzle_back_pressure", "Why does a nozzle stop listening?", "once the throat is
   sonic the back pressure only moves a shock") · ("duct_friction_heating", "Can friction speed a flow up?", "friction and heat
   push every duct flow toward Mach 1") · ("mach_cone", "Why do you hear it only after it has passed?", "the cone is the edge of
   all the circles of old news") · ("oblique_shock_wedge", "One wedge, how many shocks?", "two below the largest turn, none
   above it") · ("prandtl_meyer_fan", "How does a supersonic stream turn a corner?", "one Mach wave at a time; ν(M) keeps the
   count") · ("supersonic_flat_plate", "How can there be drag without friction?", "the push is square to the plate, so lift
   costs α times itself").
3. `nb.setup()`.
4. `nb.code` — **chapter imports**: `import numpy as np`, `import sympy as sp`, `import matplotlib.pyplot as plt`, `import pandas
   as pd`, `from scipy.optimize import brentq`, `from scipy.integrate import quad, solve_ivp`, `from fluidpy import
   ch15_compressible_flow as ch15`, `from fluidpy.core import gasdynamics as GD, acoustics as AC`, `from fluidpy.core import thermo
   as TH, bernoulli as BER, similarity as SIM, statics as STAT, airfoils as AF`, `from fluidpy import ch01_introduction as ch01,
   ch07_gravity_waves as ch07, ch14_aerodynamics as ch14`, `from fluidpy.core.interact import slider_figure, animate_figure, live`,
   `from fluidpy.core.anim import animate`, `from fluidpy.core.style import COLORS, savefig`, `import sys; sys.path.insert(0,
   "scripts"); from ch15_drawings import *`; then the six dicts of header convention 10, literally; `c0 = air.c(288.15)`;
   `rho0 = 101325.0/(air.R*288.15)`; `pd.Series(ch15.illustrative_inputs())` shown once. *explain:* one line per import ("`GD`
   and `AC` are the two new core modules; `ch15` re-exports both"). *expect:* c0 = 340.3 m/s, rho0 = 1.225 kg/m³, air.cp =
   1004.7 J/(kg K).
5. `nb.md` — **⚠️ Conventions in this chapter** (before C01): (a) `ch15.conventions_table()` — the nine rows of header
   convention 4 (subscript 0 · the star · stations · σ, δ against β, θ · the sign of a turn · f, q against Fanno and Rayleigh ·
   the wall force · S against s · the lapse rate) with "here" and "most other sources" side by side (traps **T1**, **T2**,
   **T6**, **T13**, **T14**, **T11**); (b) the letter table of header convention 5 (trap **T5**); (c)
   `pd.DataFrame(ch15.book_slips()).T` (17 rows, `kind` shown) and `pd.DataFrame(ch15.traps())` (22 rows) with the sentence
   "slips are false as printed; traps are true but easy to misread; three slips change the physics and are taught corrected:
   $u\,du=-dp/\rho=-c^2d\rho/\rho$ (15.32, third member with its minus sign restored), $-\dfrac{dp}{p_1}=\dfrac{1+(\gamma-1)M^2}{1-M^2}\,df,\ \ \dfrac{du}{u}=-\dfrac{d\rho}{\rho}=\dfrac{p_1}{p}\dfrac{df}{1-M^2}$ (15.47, numerator sign corrected), $\dfrac{dT}{T_1}=\left(\dfrac{1-\gamma M^2}{1-M^2}\right)dq,\ \ \dfrac{u\,du}{h_1}=\dfrac{(\gamma-1)M^2}{1-M^2}\,dq$ (15.48, second relation's sign corrected)"; (d) "**Which inputs?** All worked numbers use our own illustrative inputs, never the book's"; (e)
   "**Which names?** The book never says 'Fanno' or 'Rayleigh'; where we add those names and their usual parameters they are
   labelled *not in the book — ours*."

### A.1 §15.1 Introduction — N01, N02, N03, R01–R06 (open the C01 story)
1. `nb.section("15.1", "Introduction", intro="**What is this section about?** What 'compressible' adds, the gas we assume all
   chapter long, the one number everything is written in, and three surprises the chapter promises to explain.")`
2. `nb.note` — **N01 [B]** "Gas dynamics: flows fast enough for the density to change. Standing assumptions: a perfect gas
   with constant γ; no gravity; no friction outside thin boundary layers (except in §15.8, where friction is the subject).
   **Three promises:** (1) a widening duct speeds up a supersonic stream — kept in C06; (2) friction speeds up a subsonic duct
   flow — kept in C12; (3) adding heat can lower the temperature — kept in C12." A three-row table (promise · block · the
   equation that keeps it: $\dfrac{du}{u}=-\dfrac{1}{1-M^2}\dfrac{dA}{A}$ (15.34) · $-\dfrac{dp}{p_1}=\dfrac{1+(\gamma-1)M^2}{1-M^2}\,df,\ \ \dfrac{du}{u}=-\dfrac{d\rho}{\rho}=\dfrac{p_1}{p}\dfrac{df}{1-M^2}$ (15.47, numerator sign corrected) · $\dfrac{dT}{T_1}=\left(\dfrac{1-\gamma M^2}{1-M^2}\right)dq,\ \ \dfrac{u\,du}{h_1}=\dfrac{(\gamma-1)M^2}{1-M^2}\,dq$ (15.48, second relation's sign corrected)).
3. `nb.recap("R01", "The Mach number", "$M\equiv U/c$ (4.111): the flow speed over the local speed of sound. **New here:** M is a *local*
   number that changes along a duct, and every result of this chapter is written in it.", where="Ch. 4 §4.11 (C15)")` +
   `nb.code`: `SIM.mach_number(250.0, T=220.0)`, `GD.mach_from_speed(250.0, 220.0, air)`. *expect:* 0.841 both (c = 297.3 m/s at
   220 K) — an airliner at cruise.
4. `nb.recap("R02", "The speed of sound", "$c^2\equiv(\partial p/\partial\rho)_s$ (1.19): how much pressure it takes to squeeze the gas at constant entropy. **New
   here:** it is the coefficient of the first term of the expansion of p(ρ) in C01 — 'sound speed' and 'stiffness at constant
   entropy' are one idea.", where="Ch. 1 §1.8 (C36)")` + `nb.code`: `TH.sound_speed_from_eos(lambda rho, s: 101325.0*(rho/rho0)**1.4,
   rho0)`. *expect:* 340.3 m/s.
5. `nb.recap("R03", "Compressibility scales with M²", "$\nabla\cdot\mathbf u=-M^2\left(\dfrac{\rho_0}{\rho}\right)\dfrac{D}{Dt}\left(\dfrac{p-p_0}{\rho_0U^2}\right)$ (4.110): the divergence of the velocity — the rate of volume change —
   is of order $M^2$. Density changes are about $\\tfrac12M^2$. **New here:** this chapter is everything above that line; the
   exact ratio is $\dfrac{\rho_0}{\rho}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{1/(\gamma-1)}$ (15.30) in C04.", where="Ch. 4 §4.11; Ch. 14 (N01)")` + `nb.code`: `ch14.compressibility_check(0.3*c0)`,
   `GD.stagnation_density_ratio(0.3) - 1`. *expect:* ½M² = 4.5 %; exact 4.56 %.
6. `nb.note` — **N02 [B]** five-row table (regime · range · everyday example): incompressible M ≈ 0 (a car, a breeze) ·
   subsonic 0 < M < 1 (an airliner at 0.84) · transonic M ≈ 1 (the same airliner's wing upper surface) · supersonic M > 1
   (a rifle bullet at 2.5) · hypersonic M ≫ 1 (re-entry at 25) + `nb.code`: `[ch15.mach_regime(M, transonic=(0.8, 1.2),
   hypersonic=5.0) for M in (0.1, 0.5, 1.0, 2.5, 25.0)]`. *expect:* subsonic, subsonic, transonic, supersonic, hypersonic (the function says "incompressible" only at M = 0: it is a limit, not a band). "> ⚠️ Common confusion (trap T3): the edges are conventions, not laws of
   nature. We pass ours by keyword (0.8, 1.2 and the customary 5); the book places the hypersonic edge lower."
7. `nb.note` — **N03 [C]** "The name 'Mach number' is Ackeret's, after Ernst Mach, who photographed the waves of supersonic
   bullets; 'transonic' is von Kármán's. Mach's wave-front picture returns in C13 (§15.9), Ackeret's linear theory in C16
   (§15.10)."
8. `nb.recap("R04", "The perfect gas with constant specific heats", "One screen: $e=C_vT$, $h=C_pT$, $p=\\rho RT$; $C_v=\dfrac{R}{\gamma-1},\ C_p=\dfrac{\gamma R}{\gamma-1}$ (15.1, specific-heat line); $c=\sqrt{\gamma RT}=\sqrt{\gamma p/\rho}$ (15.1, sound-speed line);
   $S_2-S_1=C_p\ln\dfrac{T_2}{T_1}-R\ln\dfrac{p_2}{p_1}=C_v\ln\dfrac{T_2}{T_1}-R\ln\dfrac{\rho_2}{\rho_1}$ (15.1, entropy line). **New here:** entropy is written S (s in Ch. 1); its two forms are the tool of C04 and C08.", where="Ch. 1 §1.8–§1.9
   (C28–C47)")` + a table (relation · `TH` function) + `nb.code`: `air = GD.PerfectGas(1.4, 287.05)`; `air.cp`, `air.cv`, `air.c(288.15)`.
   *expect:* 1004.7, 717.6 J/(kg K), 340.3 m/s.
9. `nb.recap("R05", "Isentropic relations", "$\dfrac{p_2}{p_1}=\left(\dfrac{\rho_2}{\rho_1}\right)^{\gamma},\ \dfrac{T_2}{T_1}=\left(\dfrac{\rho_2}{\rho_1}\right)^{\gamma-1}=\left(\dfrac{p_2}{p_1}\right)^{(\gamma-1)/\gamma}$ (15.2) — between any two states with the same entropy.", where="Ch. 1 (C45, C46)")` +
   `nb.primer("the gamma exponents", "Three exponents appear all chapter. For γ = 1.4: $\\gamma/(\\gamma-1)=3.5$,
   $1/(\\gamma-1)=2.5$, $(\\gamma+1)/(2(\\gamma-1))=3$. For a monatomic gas (γ = 5/3): 2.5, 1.5, 2. The pressure exponent is the big
   one; the classic slip is to invert them.", code="g = 1.4                                   # ratio of specific heats of
   air\nprint(g/(g-1), 1/(g-1), (g+1)/(2*(g-1)))    # 3.5 2.5 3.0\nprint(0.5**3.5, 0.5**2.5)                # halve T at constant
   entropy: p falls to 8.8 %, rho to 17.7 %")` (**P355**) + `TH.isentropic_ratios(0.5)`.
10. `nb.recap("R06", "Properties of air", "$R=287\ \mathrm{m^2/(s^2K)},\ C_v=717,\ C_p=1004,\ \gamma=1.40$ (15.3) are the chapter's rounded values. > ⚠️ Common confusion (trap T4): 1004/717 =
    1.4003, not 1.40; the project's constants (R = 287.05, $C_p$ = 1004.7) differ in the fourth digit. `GD.AIR_BOOK` is passed
    explicitly whenever a book number is reproduced (only in private tests).", where="Ch. 1 constants")` + `nb.code`: `1004/717`,
    `GD.AIR_BOOK.cp`, `air.cp`. *expect:* 1.4003; 1004.5; 1004.7.

### A.2 §15.2 Acoustics — C01, C02
#### C01 — Sound is the linear limit: the wave equation (15.12)
1. `nb.section("15.2", "Acoustics", intro="**What is this section about?** A sound is a disturbance so weak that its products
   with itself can be dropped. That turns mass and momentum into one linear equation whose only speeds are U ± c — and its
   solution is two shapes sliding apart.")`
2. `nb.core("C01", "Sound is the linear limit: $\dfrac1{c^2}\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)^2p'-\dfrac{\partial^2p'}{\partial x_j\partial x_j}=0$ (15.12)", question="Why does a faint squeeze of the air travel at one fixed
   speed, whatever its shape?")`
3. `nb.md` — **Plain words:** "Clap your hands. A thin layer of air is squeezed, pushes on the next layer, which pushes on the
   next. Nobody across the room feels a wind; they feel the *news* of the squeeze arriving. How fast does news travel in a
   gas, and what decides it? The answer is the stiffness of the gas against a quick squeeze — and the same number c will
   reappear in every later section, as the speed the flow is compared with."
4. `nb.md` — **The idea** (ASCII): `full laws (nonlinear)  ──►  split: state = steady + small ′  ──►  drop ′×′  ──►  two linear
   laws  ──►  eliminate u′  ──►  one wave equation`; beside it "small means: $p'/\rho_0c^2\ll1$".
5. `nb.note` — **N04 [B]** "Acoustics = the linear theory of small, isentropic disturbances about a steady reference state.
   'Small' means products of two disturbances are dropped; $p'/\rho_0c^2\ll1$ (15.7) makes it precise."
6. `nb.note` — **N05 [B]** the split $u_i=U_i+u_i',\ p=p_0+p',\ \rho=\rho_0+\rho',\ T=T_0+T'$ (15.4) + 🔁 "basic state + small disturbance: the recipe of Ch. 11 (C01)". "> ⚠️ Common
   confusion (trap T1): subscript 0 means three things in this chapter" + the three-row table (§15.1 a reference value · §15.2
   the undisturbed state · from §15.4 the stagnation state; other sources write $p_t$, $T_t$ for the last).
7. **(re-enter C01)** `nb.recap("R07", "Sound speed at the reference state", "$c^2=\gamma p_0/\rho_0=\gamma RT_0$ (15.6). **New here:** $T_0$ is the *undisturbed*
   temperature (trap T1). > ⚠️ **slip #1 — the book prints** 'monotonic' gases**; the correct form is** monatomic.", where="Ch. 1
   (C47)")` + `TH.perfect_gas_sound_speed(288.15)`. **(re-enter C01)** `nb.recap("R08", "A number to remember", "Air at 15 °C: c ≈
   340 m/s. It is the c of every worked number in C01, C02 and C13.", where="Ch. 1 benchmark")`. *expect:* 340.3 m/s.
8. `nb.md` — 🔁 "Taylor expansion to first order (knowledge/primers.md: P26, ch01) and orders of smallness (P68, ch02): keep
   terms with one prime, drop terms with two. Index notation and the summation convention: Ch. 2."
9. `nb.derivation("D01", …)` — Part F D01 (10 steps), ref "15.10": holds **N06 [B]** (step 1), **N07 [B]** (step 2, boxed),
   **N11 [B]** (steps 3–5), **N12 [B]** (steps 6–8; slip #2 said once: "the book writes $P$ and $\rho_o$ here; we write $p_0$,
   $\rho_0$ throughout"), **N13 [B]** (step 9; slip #3), **N14 [B]** (step 10).
10. `nb.code` — **N07 numbers:** `AC.acoustic_density(1.0, c0)`, `AC.acoustic_temperature(1.0, 101325.0, 288.15)`,
    `AC.acoustic_impedance(rho0, c0)`. *expect:* ρ′ = 8.64 × 10⁻⁶ kg/m³, T′ = 8.1 × 10⁻⁴ K, $\rho_0c$ = 416.9 kg/(m² s).
    *explain:* "a 1 Pa sound — loud speech at arm's length — changes the density by 7 parts in a million."
11. `nb.note` — **N08 [B]** "$p'/\rho_0c^2\ll1$ (15.7): the condensation $s=\rho'/\rho_0=p'/\rho_0c^2$ must be small. What was dropped: the relative
    error of $p'\cong c^2\rho'$ (15.5) is about $\tfrac12(\gamma-1)s$ (ours)." + `nb.code`: `AC.condensation(1.0, rho0, c0)`,
    `[AC.linearisation_error(s*rho0*c0**2, rho0, c0) for s in (1e-3, 1e-2, 1e-1)]`. *expect:* 7.05 × 10⁻⁶; 0.02 %, 0.2 %, 2.0 % (the estimate ½(γ − 1)s; the exact isentrope gives 0.019 %, 0.20 %, 1.9 %). + `nb.plotly` —
    **F1** `slider_figure`: amplitude s from 10⁻³ to 0.5 (25 log steps) → the exact isentrope `AC.isentropic_pressure_perturbation(
    order="exact")`, the straight line (order 1) and the parabola (order 2) against ρ′/ρ₀, gap in the legend. *see / read /
    change* ("…at s = 0.1 the line is 2 % low: that is a 174 dB sound — acoustics is a very good approximation for anything
    you can listen to").
12. `nb.primer("decibels and sound pressure level", "Loudness is quoted on a logarithmic scale: $L=20\\log_{10}(p_{rms}/20\\,\\mu
    \\mathrm{Pa})$ dB. Every 20 dB is a factor 10 in pressure. Here it turns 'quiet' and 'loud' into pascals.", code="p_ref = 20e-6
    # reference pressure [Pa]: about the faintest audible sound\nfor L in (60, 94, 140):                          # speech, 1 Pa,
    a jet engine nearby\n    print(L, p_ref*10**(L/20))                    # 0.02 Pa, 1.0 Pa, 200 Pa")` (**P356**) + `nb.note` —
    **N09 [B]** "$\rho_0c^2=\gamma p_0$ = 141.9 kPa for sea-level air. Speech at 60 dB: s = 1.4 × 10⁻⁷. A jet engine at 140 dB:
    s = 1.4 × 10⁻³ — still small. $p'>0$ is a **compression**, $p'<0$ an **expansion** (rarefaction)." `AC.sound_pressure_level(1.0)` (94.0 dB), `AC.condensation`.
13. `nb.note` — **N10 [B]** "Newton assumed the squeeze happens at constant temperature: $c_{\text{Newton}}=\sqrt{RT_0}$, which is
    $1/\sqrt\gamma$ of the right value. Laplace's correction: a sound wave passes too quickly for heat to flow, so the
    squeeze is isentropic." + `AC.newton_sound_speed(288.15, air.R)`. *expect:* 287.6 m/s — 15.5 % low.
14. `nb.primer("an operator in brackets, applied twice", "$\\left(\\dfrac{\\partial}{\\partial t}+U\\dfrac{\\partial}{\\partial
    x}\\right)$ is an instruction: 'take the time derivative, add U times the x-derivative'. Squared means *do it twice*, not
    'square the result'. With constant U two such instructions may be done in either order.", code="import sympy as sp
    # symbolic derivatives\nx, t, U = sp.symbols('x t U'); p = sp.Function('p')(x, t)\nD = lambda F: sp.diff(F, t) + U*sp.diff(F, x)
    # the operator as a function\nprint(sp.expand(D(D(p))))                       # p_tt + 2 U p_xt + U**2 p_xx")` (**P354**;
    reminder of P177 "operator elimination", ch07) + 🔁 "the Laplacian $\partial^2/\partial x_j\partial x_j$: Ch. 2".
15. `nb.derivation("D02", …)` — Part F D02 (6 steps), ref "15.12".
16. `nb.worked_example("a travelling cosine, by hand", "Take c = 340 m/s, U = 170 m/s (M = 0.5) and try $p'=\\cos(x-(U+c)t)$ in
    $\dfrac1{c^2}\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)^2p'-\dfrac{\partial^2p'}{\partial x_j\partial x_j}=0$ (15.12) in one dimension. 1. One application of the bracket: $\\partial_t$ gives $+(U+c)\\sin$, $U\\partial_x$ gives $-U\\sin$;
    sum $c\\sin(\\cdot)$. 2. Apply it again: $c\\,[-(U+c)+U]\\cos=-c^2\\cos$. 3. Divide by $c^2$: $-\\cos$. 4. The second
    x-derivative is $-\\cos$. 5. Left side: $-\\cos-(-\\cos)=0$ ✓. The shape rides at U + c = 510 m/s. Try the speed 400 m/s
    instead and the bracket leaves a remainder: only U ± c work.")`
17. `nb.code` — `out = ch15.linearised_acoustics_sympy()`; show `out["continuity_dropped"]`, `out["momentum_dropped"]` (every term
    with two primes); `AC.convected_wave_residual(lambda x, t: np.exp(-(x - (170.0 + c0)*t)**2), 0.3, 0.001, c0, U=170.0)`. *expect:*
    residual ≈ 0 (≲ 10⁻⁶ of the largest term); with the wrong speed 400 m/s the residual is of order one. *explain:* 4 points.
18. `nb.check_agree` — **from scratch:** a Gaussian `p = lambda x, t: np.exp(-(x - (U + c)*t)**2)`; second differences in t, x and
    the mixed one with h = 10⁻³ (in units where c = 1, U = 0.5); `lhs = (ptt + 2*U*pxt + U**2*pxx)/c**2 - pxx`; `assert
    np.allclose(lhs, AC.convected_wave_residual(p, x0, t0, c, U=U, h=1e-3), atol=1e-8)` and `abs(lhs) < 1e-5`.
19. `nb.figure` — the residual of $\dfrac1{c^2}\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)^2p'-\dfrac{\partial^2p'}{\partial x_j\partial x_j}=0$ (15.12) for a Gaussian riding at speed s, against s from −2c to 3c at M = 0.5: two sharp
    zeros at U − c and U + c. *see / read / change* ("…raise U above c and both zeros are positive: C02 turns that into 'no
    news upstream'").
20. `nb.md` — **What would change if…** "…the amplitude were not small? The dropped products come back, compressions travel
    faster than expansions, and the wave steepens into a shock — N22 below, then §15.6."

#### C02 — d'Alembert's solution (15.15): two shapes sliding apart
1. `nb.core("C02", "Two shapes sliding apart: $p'(x,t)=f(x-ct)+g(x+ct)$ (15.15)", question="What does the most general plane sound wave look like —
   and which way does the air itself move?")`
2. `nb.md` — **Plain words:** "Pluck a long string, or pop a balloon in a pipe. The bump you made does not sit there and
   does not spread out: it splits into two copies that run off in opposite directions, each keeping its shape. That is all
   a plane sound wave can do. The air itself only shuffles: forward under a compression that runs forward, backward under
   one that runs backward."
3. `nb.md` — **The idea** (ASCII): `t = 0:        __/‾‾\__` / `t > 0:   _/‾\_   ⟵ g        f ⟶   _/‾\_` with "each half
   height; air moves → under f, ← under g".
4. `nb.note` — **N15 [B]** "$\dfrac1{c^2}\dfrac{\partial^2p'}{\partial t^2}-\dfrac{\partial^2p'}{\partial x^2}=0$ (15.13): $\dfrac1{c^2}\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)^2p'-\dfrac{\partial^2p'}{\partial x_j\partial x_j}=0$ (15.12) with $U_j=0$ in one space dimension — the second-order relative of the first-order
   wave equation of Ch. 7." `AC.wave_equation_residual`.
5. `nb.primer("the second-order wave equation and its two characteristic coordinates", "A first-order wave equation has one
   family of lines along which things are carried (Ch. 7, P174). The second-order one has two: $\\xi=x-ct$ and
   $\\zeta=x+ct$. Rewritten in ξ and ζ (chain rule, twice) it becomes $\\partial^2p'/\\partial\\xi\\,\\partial\\zeta=0$ — 'the
   ξ-slope does not depend on ζ' — whose solutions are a function of ξ plus a function of ζ.", code="import sympy as sp\nx, t, c =
   sp.symbols('x t c', positive=True)\nf, g = sp.Function('f'), sp.Function('g')\np = f(x - c*t) + g(x + c*t)                    # any
   two shapes\nprint(sp.simplify(sp.diff(p, t, 2)/c**2 - sp.diff(p, x, 2)))    # 0: it solves (15.13)")` (**P357**; reminder of P49
   chain rule, P174).
6. `nb.derivation("D03", …)` — Part F D03 (8 steps), ref "15.15".
7. `nb.note` — **N16 [B]** "$u_1'(x,t)=-\dfrac1{\rho_0}\displaystyle\int\dfrac{\partial p'}{\partial x}\,dt$ (15.14) — the momentum law $\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)u_i'+\dfrac1{\rho_0}\dfrac{\partial p'}{\partial x_i}=0$ (15.11) integrated in time. > ⚠️ Common confusion (trap T19): the integral is
   indefinite; its 'constant' is a function of x, fixed by 'the air was at rest before the wave arrived'."
   `AC.velocity_from_pressure`.
8. `nb.derivation("D04", …)` — Part F D04 (5 steps), ref "15.16": holds **N17 [B]**. Then a two-row table (right-going f:
   $u'=+p'/\rho_0c$ — a compression pushes the air the way the wave goes · left-going g: $u'=-p'/\rho_0c$) and "$\rho_0c$ is the
   acoustic impedance (our name)"; "$u'/c=p'/\rho_0c^2$ = the condensation: the air moves far slower than the wave".
9. `nb.worked_example("a bump released at rest", "c = 340 m/s, $\\rho_0$ = 1.2 kg/m³. At t = 0 the pressure is a bump F(x) of
   height 100 Pa and the air is at rest, so $\\partial p'/\\partial t=G=0$. 1. The initial-value formulas $f(x)=\dfrac12\left[F(x)-\dfrac1c\displaystyle\int_0^xG(\hat x)\,d\hat x\right],\ g(x)=\dfrac12\left[F(x)+\dfrac1c\displaystyle\int_0^xG(\hat x)\,d\hat x\right]$ (unnumbered, Exercise 15.1) with G = 0 give
   $f=g=\\tfrac12F$: two bumps of 50 Pa. 2. After 1 s their centres are at x = +340 m and x = −340 m. 3. Under the right-going
   one $u'=50/(1.2\\times340)$ = +0.12 m/s; under the left-going one −0.12 m/s. 4. At t = 0 the two velocities cancel: the air
   was at rest ✓. 5. Wave speed 340 m/s, air speed 0.12 m/s: a ratio of 3.6 × 10⁻⁴ — the condensation.")` — this is **N21 [B]**
   ("the initial-value form, stated; result given, not derived").
10. `nb.code` — `F = AC.pulse("gaussian", width=1.0, amplitude=100.0)`; `f, g = AC.dalembert_fg(F, lambda s: 0.0*s, c0)` (the two halves); `x = np.linspace(-8, 8, 801)`; `p = AC.dalembert_initial_value(x,
    t, c0, F)`; `u = AC.plane_wave_velocity(x, t, c0, rho0, f=lambda s: 0.5*F(s), g=lambda s: 0.5*F(s))`; at t = 0 check `p == F(x)`
    and `u == 0`. *expect:* max p′ = 100 Pa at t = 0, 50 Pa later; max |u′| = 0.120 m/s. *explain:* 4 points.
11. `nb.check_agree` — **from scratch:** `shape = lambda s: 50*np.exp(-s**2)`; `p_mine = shape(x - c*t) + shape(x + c*t)`; `u_mine =
    (shape(x - c*t) - shape(x + c*t))/(rho0*c)`; `assert np.allclose(p_mine, AC.dalembert(x, t, c, f=shape, g=shape))`; `assert
    np.allclose(u_mine, AC.plane_wave_velocity(x, t, c, rho0, f=shape, g=shape))`.
12. `nb.animation` — **N18 [B]** (animation **A1**) two panels on one clock (60 frames; FAST 30; video): top, a right-going pulse of fixed
    shape with p′ (orange) and u′ (blue) drawn together, displaced by ct; bottom, the bump released at rest splitting into two
    half-height pulses with u′ of opposite signs. *see / read / change.*
13. `nb.note` — **N19 [B]** "In a stream of speed U the two shapes ride at U + c and U − c: $p'(x,t)=f(x-(c+U)t)+g(x+(c-U)t)$ (15.17)." + `nb.derivation("D05", …)` —
    Part F D05 (5 steps), ref "15.17" + `AC.wave_speeds(0.5*c0, c0)`, `AC.wave_speeds(2*c0, c0)`. *expect:* (510.4, −170.1) m/s;
    (1020.9, +340.3) m/s. "> ⚠️ 'Downstream' presumes U > 0 (trap T19)."
14. `nb.note` — **N20 [B]** "**Range of influence.** For U < c one shape still travels upstream: a point feels what happens
    on both sides. For U > c both travel downstream: a point feels only what is upstream of it. That is why a choked nozzle
    stops listening (C09) and why you hear a supersonic aircraft only after it has passed (C13)." + `nb.figure`: x–t diagrams
    from `AC.xt_characteristics(U, c0, x0, t)` for M = 0.5 and M = 2 (two families, slopes U ± c; the wedge of influence of one
    point shaded). `AC.upstream_signal_possible`. *see / read / change* ("…at exactly M = 1 one family is vertical: the upstream
    wave stands still").
15. `nb.explainer("sound_pulse_dalembert", heading="What is a sound pulse made of?", why="you release the bump, scrub time and
    watch f and g separate with the tracer particles moving opposite ways; then you drag the stream through U = c and one
    family of characteristics tips over", tries=["Release a bump at rest and step time: each half is half as high.", "Click a
    point under each half: the air moves opposite ways.", "Drag the stream Mach number through 1 and watch the x–t view: the
    upstream line stands up, then leans downstream."])`
16. `nb.note` — **N22 [B]** "**Steepening (the book argues in words; the formulas are ours).** In a compression the gas is
    hotter, so c is larger, and the gas itself moves forward: each wavelet travels at $\left.\dfrac{dx}{dt}\right\rvert_{\text{wavelet}}=u+c=c_0+\dfrac{\gamma+1}{2}u,\ \ t_{\text{break}}=\dfrac{2}{(\gamma+1)\max(-\partial u_0/\partial x)}$ (ours). Crests catch up with troughs;
    a sine becomes a sawtooth; the limit is a shock. Result given, not derived." 🔁 "the shallow-water simple wave of Ch. 7
    (N81) is the γ = 2 twin: `ch07.simple_wave_evolve`." + `nb.primer("Riemann invariant", "A combination of the fluid speed
    and the sound speed that keeps its value along a moving wavelet. In a wave going one way through still gas,
    $u-2c/(\\gamma-1)$ has the same value everywhere, so $c=c_0+\\tfrac12(\\gamma-1)u$: each bit of the wave knows its own
    speed.", code="g, c0 = 1.4, 340.0                       # air\nu = np.array([-5.0, 0.0, 5.0])              # fluid speed in
    trough, node, crest [m/s]\nc = c0 + 0.5*(g - 1)*u                    # local sound speed from the invariant\nprint(u - 2*c/(g -
    1), u + c)             # the invariant is constant; wavelet speeds 334, 340, 346 m/s")` (**P358**) + `nb.code`:
    `AC.breaking_time(u0, x)` for a 1 kHz tone at 150 dB; `AC.simple_wave_speed`. *expect:* peak pressure 894 Pa, $u_a$ = 2.15
    m/s, $t_{\text{break}}$ = 21 ms: the wave breaks after 7.2 m ≈ 21 wavelengths. + `nb.animation` — **A2** (frames; stops at ¼,
    ½, 1 × $t_{\text{break}}$): `AC.simple_wave_evolve` for one wavelength, the acoustic wave as a ghost. *see / read / change*
    ("…an expansion flattens instead: only compressions become shocks").
17. `nb.note` — **N114 [C]** "The full nonlinear (Euler) equations can carry a wave *through* breaking, which the simple-wave
    formula cannot; Ch. 10's MacCormack scheme (C10) is the tool. Not built here (it is part of the backup explainer)."
18. `nb.md` — **What would change if…** "…the wave were confined in a duct whose area changes? Then we stop following the
    wave and write the three conservation laws between two stations — C03."

### A.3 §15.3 Basic Equations for One-Dimensional Flow — C03
#### C03 — The three duct laws, anchored on the energy law (15.22)
1. `nb.section("15.3", "Basic Equations for One-Dimensional Flow", intro="**What is this section about?** Mass, energy and
   momentum written between two stations of a duct. The rest of the chapter is these three lines with different terms
   switched on.")`
2. `nb.core("C03", "The three duct laws: $h_2+\tfrac12u_2^2-h_1-\tfrac12u_1^2=Q$ (15.22)", question="What must be the same at two cross-sections of a duct — and
   which of the three statements survives friction and shocks?")`
3. `nb.md` — **Plain words:** "A hair dryer, a jet engine intake, a rocket nozzle: gas goes in one end and out the other.
   We do not want the flow at every point — only one speed, pressure, density and temperature per cross-section. Three
   bookkeeping statements tie the outlet to the inlet: what goes in comes out (mass); what the gas carries in plus the heat
   added equals what it carries out (energy); and the change of momentum equals the push of the pressures and the walls."
4. `nb.md` — **The idea** (ASCII duct): `station 1 │→ u₁ p₁ ρ₁ T₁  A₁ ⟩═══════⟨  A₂  u₂ p₂ ρ₂ T₂ →│ station 2`, below it "walls:
   force F on the fluid · heat Q per kg added"; a three-line table (law · what crosses · switched on in): mass (C05, C06) ·
   energy (C04, C07, C11) · momentum (C07, C11).
5. `nb.note` — **N23 [B]** hypotheses box: (1) straight centreline; (2) area A(x) varies slowly; (3) one value of u, p, ρ, T per
   cross-section; (4) steady. + `nb.figure`: `draw_duct_control_volume()` (our drawing: the control volume between stations 1
   and 2, the outward normals, F and Q). *see / read / change.* 🔁 "control volume and the Reynolds transport theorem: Ch. 3
   (C15), Ch. 4 (C01)."
6. **(re-enter C03)** `nb.recap("R09", "Mass conservation on a control volume", "What flows in flows out: $\rho uA=\dot m=\text{const.}$ (15.18). **New here:** on
   this control volume the integral law is one line.", where="Ch. 4 §4.2 (C01)")` + `GD.mass_flow(rho, u, A)`.
7. **(re-enter C03)** `nb.recap("R10", "The steady energy integral", "$\displaystyle\int_{A^*}\rho\left(e+\tfrac12u_i^2\right)u_jn_j\,dA=\int_{A^*}u_i\tau_{ij}n_j\,dA-\int_{A^*}q_jn_j\,dA$ (15.20): the flux of internal plus kinetic energy out of the
   surface equals the rate of work of the surface stresses plus the heat conducted in. > ⚠️ Common confusion (trap T2): here
   $A^*$ is the control surface, not the sonic throat area of $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31).", where="Ch. 4 §4.8 (C10)")`
8. `nb.md` — 🔁 "power of a force = force × velocity; heat flux through a surface (knowledge/primers.md: P127, ch04); enthalpy
   $h=e+p/\rho$ (Ch. 1)".
9. `nb.derivation("D06", …)` — Part F D06 (8 steps), ref "15.22": holds **N25 [B]** (steps 1–5, with the face-by-face table:
   face · normal · what crosses) and **N26 [B]** (step 6).
10. `nb.note` — **N26 [B]** "**Flow work.** Pushing one kilogram across a face against the pressure costs $upA/\dot m=p/\rho$ of
    work. That is why $h=e+p/\rho$ replaces e in a flowing system." + `nb.code`: `air.R*293.15`, `air.cv*293.15`, `air.cp*293.15`.
    *expect:* p/ρ = 84.1 kJ/kg, e = 210.4 kJ/kg, h = 294.5 kJ/kg — flow work is 29 % of h (= (γ − 1)/γ).
11. **(re-enter C03)** `nb.recap("R11", "The steady momentum integral", "$\displaystyle\int_{A^*}\rho u_iu_jn_j\,dA=\int_{A^*}\tau_{ij}n_j\,dA$ (15.23): momentum flux out = force of the surface
    stresses.", where="Ch. 4 §4.4 (C03)")` + 🔁 "momentum flux through a surface (P114, ch04); logarithmic differentiation (P52,
    ch01); the product rule (P38)".
12. `nb.derivation("D07", …)` — Part F D07 (9 steps), ref "15.25": holds **N24 [B]** (steps 1–2), **N27 [B]** (step 3), **N28
    [B]** (steps 4–8) and **R12** (step 9).
13. `nb.note` — **N27 [B]** "> ⚠️ Common confusion (trap T11): three F's. (a) F in $-\dot mu_1+\dot mu_2=(pA)_1-(pA)_2+F$ (15.24) is the force of the walls **on the
    fluid**, positive along +x. (b) $F_f$ in $\dot m\dfrac{du}{dx}=-A\dfrac{dp}{dx}-F_f$ (15.25) is a positive *retarding* friction force per unit length. (c) F in
    $f=F/(p_1A)$ of $\rho_1u_1=\rho_2u_2,\ \ p_1+\rho_1u_1^2=p_2+\rho_2u_2^2+p_1f,\ \ h_1+\tfrac12u_1^2+h_1q=h_2+\tfrac12u_2^2$ (15.45, with $f=F/(p_1A)$, $q=Q/h_1$) is the total *retarding* force." (three-row table) + **N24** number: "1 % more area at fixed density:
    1 % less speed" (`GD.influence_coefficients(0.0, dA_over_A=0.01)["du_over_u"]` → −0.01).
14. **(re-enter C03)** `nb.recap("R12", "The one-dimensional Euler equation", "With no friction, $u\,du+\dfrac{dp}{\rho}=0$ (15.26) — Euler's equation
    along the duct.", where="Ch. 4 (C07)")` + `ch15.duct_equations_sympy()` (checks $\dfrac{d\rho}{\rho}+\dfrac{du}{u}+\dfrac{dA}{A}=0$ (15.19), $\dot m\dfrac{du}{dx}=-A\dfrac{dp}{dx}-F_f$ (15.25), $u\,du+\dfrac{dp}{\rho}=0$ (15.26) from $\rho uA=\dot m=\text{const.}$ (15.18), $-\dot mu_1+\dot mu_2=(pA)_1-(pA)_2+F$ (15.24)).
15. **(re-enter C03)** `nb.recap("R13", "Euler + Gibbs = the energy Bernoulli equation", "With $dS=0$ the Gibbs relation
    $T\\,dS=dh-dp/\\rho$ gives $dh=dp/\\rho$, so $u\,du+\dfrac{dp}{\rho}=0$ (15.26) becomes $d(h+\\tfrac12u^2)=0$: the Q = 0 case of $h_2+\tfrac12u_2^2-h_1-\tfrac12u_1^2=Q$ (15.22). > ⚠️ Common
    confusion (trap T17): **adiabatic** = no heat (Q = 0); **reversible** = no friction; **isentropic** = both. $h_2+\tfrac12u_2^2-h_1-\tfrac12u_1^2=Q$ (15.22)
    needs only 'adiabatic'; this route also needs 'frictionless'. Result given, not derived.", where="Ch. 1 (D10); Ch. 4 §4.9")`
16. `nb.worked_example("three laws, easy numbers", "Air ($C_p$ ≈ 1000 J/(kg K)) enters a short insulated duct at $u_1$ = 100 m/s,
    $T_1$ = 300 K, $\\rho_1$ = 1.0 kg/m³ through $A_1$ = 0.02 m² and leaves at $u_2$ = 300 m/s. 1. Mass $\rho uA=\dot m=\text{const.}$ (15.18): ṁ = 1.0 × 100 × 0.02 =
    2 kg/s. 2. Energy $h_2+\tfrac12u_2^2-h_1-\tfrac12u_1^2=Q$ (15.22) with Q = 0: $C_p(T_2-T_1)=\\tfrac12(u_1^2-u_2^2)$ = ½(10⁴ − 9 × 10⁴) = −40 000 J/kg, so $T_2$ = 260 K:
    the gas cooled by 40 K to pay for its speed. 3. Momentum $-\dot mu_1+\dot mu_2=(pA)_1-(pA)_2+F$ (15.24): ṁ($u_2-u_1$) = 2 × 200 = 400 N must be supplied by the
    pressure difference on the two faces plus the wall force F. 4. Nothing here asked whether the flow was frictionless —
    the energy line does not care.")`
17. `nb.code` — two stations of our nozzle (`noz`): station 1 where A = 2 × 10⁻³ m² on the subsonic branch, station 2 at the
    throat: states from `GD.isentropic_state` × ($p_0$, $T_0$); `GD.mass_flow` at both; packed as `state1 = dict(p=p1, rho=rho1, u=u1)`, `state2 = …`; `GD.total_enthalpy(T, u, air.cp)` at both; `GD.energy_residual(state1, state2)`; `GD.wall_force(state1, state2, A1, A2)`. *expect:* $M_1$ = 0.306, $p_1$ = 468.6
    kPa, $T_1$ = 392.7 K, $u_1$ = 121.5 m/s; throat 264.1 kPa, 333.3 K, 366.0 m/s; ṁ = 1.010 kg/s at both; $h_0$ = 401.9 kJ/kg at
    both; F = −426 N (the converging walls push the gas **backward** while it accelerates — the pressure does the pushing).
    *explain:* 5 points.
18. `nb.check_agree` — **from scratch:** three lines of arithmetic: `mdot1 = rho1*u1*A1; mdot2 = rho2*u2*A2`; `dh0 = cp*T2 +
    0.5*u2**2 - cp*T1 - 0.5*u1**2`; `F_mine = mdot1*(u2 - u1) - (p1*A1 - p2*A2)`; `assert np.allclose([mdot1, F_mine], [mdot2,
    GD.wall_force(s1, s2, A1, A2)])` and `abs(dh0) < 1e-6*cp*T1`.
19. `nb.figure` — **the budget as bars** (the block's visual): three pairs of bars for our two stations — mass flux in / out;
    $h+\tfrac12u^2$ in / out (stacked: h rose, ½u² blue); momentum: ṁu₁ + p₁A₁ + F against ṁu₂ + p₂A₂ (stacked, F muted).
    *see / read / change* ("…add heat Q: only the middle pair changes height; add friction: only the third — and the middle
    pair still balances").
20. `nb.md` — **What would change if…** "…you slowed the gas to rest without heat or friction? The energy line says its
    temperature would rise to a fixed value $T_0$ — the ruler of the next block."


### A.4 §15.4 Reference Properties in Compressible Flow — C04, C05
#### C04 — The stagnation state and its ratios (15.28)
1. `nb.section("15.4", "Reference Properties in Compressible Flow", intro="**What is this section about?** Two imagined states
   that every point of a flow carries with it: the state it would reach if brought gently to rest (stagnation), and the
   state it would reach if brought gently to Mach 1 (sonic). All tables and nozzle results are ratios to them.")`
2. `nb.core("C04", "The stagnation ruler: $\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28)", question="How hot and how compressed would this gas be if you stopped
   it gently — and which of those two numbers does a real flow keep?")`
3. `nb.md` — **Plain words:** "Hold a thermometer out of an aircraft window at Mach 0.85. It does not read the air's
   temperature: the air that hits it is brought to rest, its kinetic energy becomes heat, and the reading is about 32 K too
   high. That higher temperature — and the higher pressure a pitot tube reads — is a property of *each point* of the flow,
   defined by an imagined experiment. It is the ruler of the chapter: $T_0$ is what an adiabatic flow keeps; $p_0$ is what it
   loses whenever entropy is made."
4. `nb.md` — **The idea** (ASCII): `moving gas (T, p, u)  ──[stop it: no heat, no friction]──►  at rest (T₀, p₀, 0)`; two-column
   table: "$T_0$ needs only 'no heat' → kept across friction and shocks" / "$p_0$ also needs 'no friction' → falls when entropy rises".
5. `nb.note` — **N29 [B]** "Stagnation (total) state: the state the local gas would reach if brought to rest isentropically;
   stagnation enthalpy $h_0\equiv h+\tfrac12u^2$. > ⚠️ Common confusion (trap T1): this is the third meaning of the subscript 0;
   most tables write $T_t$, $p_t$." `GD.total_enthalpy`.
6. **(re-enter C04)** `nb.recap("R14", "Stagnation temperature", "$C_pT_0=C_pT+\tfrac12u^2$ (15.27): the energy law $h_2+\tfrac12u_2^2-h_1-\tfrac12u_1^2=Q$ (15.22) with Q = 0, between the point and
   its imagined state of rest.", where="Ch. 4 §4.9 (N94)")` + `BER.stagnation_temperature(220.0, 250.0, cp=air.cp)`. *expect:* 251.1 K
   (M = 0.841: +31 K).
7. `nb.derivation("D08", …)` — Part F D08 (6 steps), ref "15.28": holds **N30 [B]** (steps 4–5), **N31 [B]** (step 6). 🔁 the
   gamma exponents (P355, in C01).
8. `nb.worked_example("Mach 2 by hand", "γ = 1.4, M = 2. 1. $\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28): $T_0/T$ = 1 + 0.2 × 4 = **1.8**. 2. $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29): $p_0/p$ = 1.8^3.5 =
   1.8³ × √1.8 = 5.832 × 1.342 = **7.82**. 3. $\dfrac{\rho_0}{\rho}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{1/(\gamma-1)}$ (15.30): $\\rho_0/\\rho$ = 1.8^2.5 = 1.8² × √1.8 = 3.24 × 1.342 = **4.35**. 4. Check with
   the gas law: 7.82/4.35 = 1.80 ✓. A Mach-2 stream at 288 K would heat to 519 K if stopped.")`
9. `nb.code` — `GD.stagnation_temperature_ratio(2.0)`, `GD.stagnation_pressure_ratio(2.0)`, `GD.stagnation_density_ratio(2.0)`,
   `GD.isentropic_state(2.0)`, `GD.mach_from_pressure_ratio(1/7.8244)`. *expect:* 1.8, 7.8244, 4.3469; `p_over_p0` = 0.1278 (the
   public NACA Report 1135 row for M = 2); inverse → 2.000. *explain:* 4 points.
10. `nb.check_agree` — **from scratch:** `t = 1 + 0.5*(g - 1)*M**2; mine = dict(T_over_T0=1/t, p_over_p0=t**(-g/(g - 1)),
    rho_over_rho0=t**(-1/(g - 1)))`; `lib = GD.isentropic_state(M, g)`; `assert all(np.allclose(mine[k], lib[k]) for k in mine)` for M
    = `np.linspace(0, 5, 51)`.
11. `nb.primer("dynamic pressure in a gas", "$\\tfrac12\\rho U^2=\\tfrac12\\gamma pM^2$, because $c^2=\\gamma p/\\rho$ gives
    $\\rho U^2=\\rho c^2M^2=\\gamma pM^2$. It lets us trade a speed for a Mach number in any momentum term.", code="g, p, M = 1.4,
    101325.0, 0.3                 # air at sea level, Mach 0.3\nrho, c = 1.225, 340.3                     # density, sound
    speed\nprint(0.5*rho*(M*c)**2, 0.5*g*p*M**2)         # 6384 Pa both ways")` (**P359**)
12. `nb.note` — **N32 [B]** "For small M, $p_0=p+\tfrac12\rho u^2\left(1+\tfrac14M^2+\dots\right)$ (ours) and $\rho_0/\rho=1+\tfrac12M^2+\dots$ (ours; result given, not derived): the
    incompressible pitot formula of Ch. 4 is the first term." + `nb.code`: `sp.series((1 + (g-1)/2*M**2)**(g/(g-1)), M, 0, 6)`
    (reminder: sympy `series`, P117; the binomial series, P342); `GD.pitot_compressibility_factor(0.3)`,
    `GD.mach_for_density_change(0.05)`. *expect:* factor 1.0227 (the incompressible formula reads 2.3 % low at M = 0.3); the
    density differs by 5 % at M = 0.314. + `nb.plotly` — **F2** `slider_figure`: γ in (1.1 … 1.67, 20 steps) → $p_0/p$ against M
    (0 … 2) with the incompressible $1+\tfrac12\gamma M^2$ as a ghost and the M = 0.3 line. *see / read / change.*
13. `nb.derivation("D09", …)` — Part F D09 (5 steps), ref "unnumbered, §15.4": holds **N33 [B]**. + "> ⚠️ Common confusion (trap
    T17): adiabatic is not isentropic — $T_0$ is kept, $p_0$ is not." + `nb.code`: `GD.entropy_change_from_stagnation(0.9, gas=air)`,
    `GD.stagnation_pressure_loss(30.24, gas=air)`. *expect:* a 10 % loss of $p_0$ at fixed $T_0$ is ΔS = 30.2 J/(kg K); back → 0.900.
14. `nb.figure` — **the ruler** (the block's visual): for M = 0.3, 0.85, 2, 3 three bars each — $T/T_0$ (rose), $p/p_0$
    (orange), $\rho/\rho_0$ (teal) — with the amber line at 1 = the stagnation value. *expect:* 0.982 / 0.939 / 0.956 · 0.874 /
    0.624 / 0.714 · 0.556 / 0.128 / 0.230 · 0.357 / 0.027 / 0.076. *see / read / change* ("…pressure always falls fastest: its
    exponent is 3.5").
15. `nb.note` — **N34 [C]** "A large reservoir at rest holds the stagnation values of everything that flows out of it
    isentropically; it is the left end of the nozzle figures of C06, C09 and C10." · **N38 [C]** "Let $T\to0$ in $\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28): the
    speed tends to $u_{\max}=\sqrt{2C_pT_0}$ (ours) — `GD.max_speed(288.15, air.cp)` = 761 m/s; the same limit caps the
    Prandtl–Meyer turn in C15."
16. `nb.md` — **What would change if…** "…instead of stopping the gas you brought it to exactly Mach 1? You get a second
    reference state — and with it the area a flow needs. C05."

#### C05 — The sonic state and the area ratio (15.31)
1. `nb.core("C05", "One area, two Mach numbers: $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31)", question="How much room does a given flow need at each Mach
   number — and why is there a least amount?")`
2. `nb.md` — **Plain words:** "Push a fixed number of kilograms per second down a duct. Where the gas is slow it needs a
   wide passage. Speed it up and it needs less — until Mach 1. Beyond that the gas thins out so quickly that it needs
   *more* room again. So the required area has a minimum, exactly at Mach 1: that minimum is the throat, and a throat of
   given size can pass only so much."
3. `nb.md` — **The idea** (ASCII U-curve): `A/A* │\        /` … `1 ┼──\__/──` … `      M=1`, with "one height → two Mach numbers
   (subsonic, supersonic)".
4. `nb.note` — **N35 [B]** "**Sonic (critical) state** — not in the book as formulas, ours: put M = 1 in $\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28), $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29), $\dfrac{\rho_0}{\rho}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{1/(\gamma-1)}$ (15.30):
   $\dfrac{T^*}{T_0}=\dfrac{2}{\gamma+1},\ \dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)},\ \dfrac{\rho^*}{\rho_0}=\left(\dfrac{2}{\gamma+1}\right)^{1/(\gamma-1)}$ (ours, from (15.28)–(15.30) at $M=1$). For γ = 1.4: 0.833, 0.528, 0.634. > ⚠️ Common confusion (trap T2): the star now means 'if brought to M = 1
   isentropically'." `GD.sonic_ratios()`, `GD.critical_pressure_ratio()`.
5. `nb.derivation("D10", …)` — Part F D10 (9 steps), ref "15.31" ("the book states this without proof").
6. `nb.worked_example("Mach 2 and Mach 0.5 by hand", "γ = 1.4, exponent ½(γ+1)/(γ−1) = 3. 1. M = 2: bracket = (2/2.4) × 1.8 = 1.5;
   1.5³ = 3.375; ÷ M: **1.6875**. 2. M = 0.5: bracket = (1/1.2) × 1.05 = 0.875; 0.875³ = 0.670; ÷ 0.5: **1.340**. 3. M = 1:
   bracket = 1: **1**. Both sides of Mach 1 need more area than the throat.")`
7. `nb.primer("a function with two inverses", "A U-shaped function has two x for each y above its minimum. A root-finder
   needs to be told which: give it a bracket on one side of the minimum. Here the physics chooses the side.", code="from
   scipy.optimize import brentq                 # bracketed root finder (reminder: P108)\nA = lambda M: (1/M)*((2/2.4)*(1 +
   0.2*M*M))**3     # A/A* for gamma = 1.4\nprint(brentq(lambda M: A(M) - 2.0, 1e-8, 1.0))   # 0.306 (subsonic side)\nprint(brentq(
   lambda M: A(M) - 2.0, 1.0, 100.0))   # 2.197 (supersonic side)")` (**P360**) + `nb.note` — **N36 [B]** "(ours) the inverse:
   `GD.mach_from_area_ratio(2.0, branch="subsonic")` = 0.306, `branch="supersonic"` = 2.197. > ⚠️ Common confusion (trap T8, first
   half): one area ratio, two Mach numbers — which one occurs is decided downstream, by the back pressure (C09–C10)."
8. `nb.code` — `GD.area_ratio(np.array([0.5, 1.0, 2.0, 3.0]))`; both inverses for 1.5 and 2.0. *expect:* 1.3398, 1, 1.6875,
   4.2346; (0.4303, 1.8541), (0.3059, 2.1972). *explain:* 4 points.
9. `nb.check_agree` — **from scratch:** the route of D10 as a product of four ratios, `A_mine = (rho_star/rho0)*(rho0/rho)*
   np.sqrt(T_star/T)*(1/M)` with each ratio written in one line; a ten-line bisection on [1e-6, 1] and on [1, 50]; `assert
   np.allclose(A_mine, GD.area_ratio(M))`; `assert np.allclose([M_sub, M_sup], [GD.mach_from_area_ratio(2.0, branch="subsonic"),
   GD.mach_from_area_ratio(2.0, branch="supersonic")], atol=1e-9)`.
10. `nb.primer("a generated table as a pandas DataFrame", "The classical tool is a printed table of ratios against Mach
    number. We generate ours: a DataFrame is a table with named columns; `.loc` picks rows; reading it 'backwards' means
    finding the row whose column is closest to a value.", code="tab = GD.isentropic_table(np.arange(0.5, 3.01, 0.5))     # our
    own table, gamma = 1.4\nprint(tab.round(4)) # columns 'M', 'p/p0', 'rho/rho0', 'T/T0', 'A/A*'\nprint(tab.loc[(tab.
    ['A/A*'] - 1.69).abs().idxmin(), 'M'])   # read backwards: nearest row to A/A* = 1.69")` (**P361**) + `nb.note` —
    **N37 [B]** "The book's table is not reproduced; ours is generated for γ = 1.4 and 5/3 and checked against three rows of
    the public-domain NACA Report 1135 (M = 0.5, 2, 3). Note the last line of the demo: reading backwards returns only the
    nearest row — and only one of the two branches."
11. `nb.derivation("D11", …)` — Part F D11 (7 steps), ref "ours": holds **N35 [B]** (steps 1–2) and **N39 [B]** (steps 3–7).
12. `nb.code` — `GD.mass_flux_function(np.array([0.5, 1.0, 2.0]))`; `GD.choked_mass_flow(noz["p0"], noz["T0"], noz["A_t"], air)`.
    *expect:* 0.5111, 0.6847, 0.4058; ṁ = 1.010 kg/s for our 10 cm² throat at 5 bar and 400 K. *explain:* "the flux function is
    $\sqrt\gamma\,(2/(\gamma+1))^3$ = 0.6847 at its peak, M = 1."
13. `nb.plotly` — **F3** two `slider_figure`s: (a) γ (1.1 … 1.67, 20 steps) → $A/A^*$ against M on a log M axis, the two roots of
    A/A* = 2 marked; (b) γ → the mass flux function $\dfrac{\dot m\sqrt{RT_0}}{p_0A}=\sqrt\gamma\,M\left(1+\dfrac{\gamma-1}{2}M^2\right)^{-(\gamma+1)/(2(\gamma-1))}$ (ours) against M with its maximum at M = 1. *see / read / change* ("…the
    two figures are one curve turned upside down: flux × (A/A*) is constant").
14. `nb.md` — **What would change if…** "…you asked not 'how much area' but 'which way does the speed change when the area
    changes'? Differentiate — and a factor $1-M^2$ appears. C06."

### A.5 §15.5 Area-Velocity Relationship in One-Dimensional Isentropic Flow — C06
#### C06 — The area–velocity relation (15.34)
1. `nb.section("15.5", "Area-Velocity Relationship in One-Dimensional Isentropic Flow", intro="**What is this section about?**
   Which way the speed changes when a duct narrows or widens. Below Mach 1 a nozzle converges; above Mach 1 it must
   diverge. The first of the chapter's three surprises.")`
2. `nb.core("C06", "A supersonic nozzle widens: $\dfrac{du}{u}=-\dfrac{1}{1-M^2}\dfrac{dA}{A}$ (15.34)", question="Why does squeezing a duct speed a slow flow up but slow a
   supersonic one down?")`
3. `nb.md` — **Plain words:** "Put your thumb over a garden hose: the water speeds up. Everyone's intuition says 'squeeze =
   faster'. Now look at the bell of a rocket engine: the gas is supersonic and the bell *opens out* — and the gas keeps
   accelerating. The mass flow ρuA is fixed. In water ρ cannot change, so u must rise when A falls. In a gas ρ falls as u
   rises — and above Mach 1 it falls faster than u rises, so the product ρu *drops* and the gas needs more room."
4. `nb.md` — **The idea**: a two-column table ("M < 1: speed +1 % costs density −M² % < 1 % → ρu rises → needs less area" ·
   "M > 1: speed +1 % costs density −M² % > 1 % → ρu falls → needs more area").
5. `nb.note` — **N40 [B]** "A **nozzle** is whatever speeds the flow up and drops its pressure; a **diffuser** does the
   opposite. They are defined by what they do, not by their shape — and a diffuser has nothing to do with molecular
   diffusion."
6. **(re-enter C06)** `nb.recap("R15", "The constant-density limit", "With ρ fixed, $\dfrac{d\rho}{\rho}+\dfrac{du}{u}+\dfrac{dA}{A}=0$ (15.19) is $\\dfrac{dA}{A}+\\dfrac{du}{u}=0$: the
   garden hose. > ⚠️ **slip #5 — the book prints** $dA/dx+du/dx=0$**; the correct form is** $\\dfrac1A\\dfrac{dA}{dx}+\\dfrac1u
   \\dfrac{du}{dx}=0$ (a metre cannot be added to a per-second).", where="Ch. 4 continuity")` + `nb.code` (units cell): the two
   printed terms have units m and 1/s; the logarithmic ones are both 1/m.
7. `nb.primer("differentials of algebraic laws as a small linear system", "Take the differential of each law and divide by
   the law: every law becomes a *linear* equation between the relative changes $d\\rho/\\rho$, $du/u$, $dp/p$, … Then solve
   them like simultaneous equations.", code="import sympy as sp\ndr, du, dA, M = sp.symbols('dr du dA M')              # d rho/rho,
   du/u, dA/A, Mach number\nsol = sp.solve([dr + du + dA, dr + M**2*du], [dr, du])   # mass (15.19) and Euler (15.33)\nprint(sol[du])
   # dA/(M**2 - 1), i.e. (15.34)")` (**P362**; reminder of P52 logarithmic differentiation)
8. `nb.derivation("D12", …)` — Part F D12 (7 steps), ref "15.34": holds **N41 [B]** (steps 1–3, the minus sign restored) and
   **N42 [B]** (step 4).
9. `nb.md` — "> ⚠️ **slip #4 — the book prints** $u\,du=-dp/\rho=c^2d\rho/\rho$ (15.32)**; the correct form is** $u\,du=-dp/\rho=-c^2d\rho/\rho$ (15.32, third member with its minus sign restored). Why:
   $u\,du=-dp/\rho$ and $dp=+c^2d\rho$ leave no choice, and the book's own next line $d\rho/\rho=-M^2\,du/u$ (15.33) needs the minus." + `nb.code`:
   `GD.influence_coefficients(2.0, dA_over_A=0.01)`, the same with `printed_15_32=True`. *expect:* correct: `du_over_u` = +0.00333,
   `drho_over_rho` = −0.01333; with the printed sign `du_over_u` = −0.00200 and `drho_over_rho` = −0.00800: density and speed change in the *same* sense ($d\rho/\rho=+M^2du/u$), the factor $1-M^2$ becomes $1+M^2$, and a widening duct would slow every flow — no supersonic nozzle could exist. "A test must fail on it." **N42 read aloud:** at M = 2 a 1 % gain in speed costs 4 % of the density; at
   M = 0.3 only 0.09 %.
10. `nb.note` — **N43 [B]** 2 × 2 table (converging / diverging × M < 1 / M > 1 → nozzle or diffuser) + `nb.figure`:
    `draw_four_ducts()` with `GD.duct_role(M, dA_sign)` printed under each. *see / read / change.*
11. `nb.note` — **N44 [B]** "At a throat dA = 0, so $\dfrac{du}{u}=-\dfrac{1}{1-M^2}\dfrac{dA}{A}$ (15.34) gives either du = 0 or M = 1: **sonic speed can occur only at a
    throat, but a throat need not be sonic.**" (the 0/0 at M = 1 is a finite limit — l'Hôpital, P343 reminder) + `nb.figure`:
    three solutions M(x) on one shape, `GD.nozzle_isentropic(x, GD.laval_area, branch=…)` for `"subsonic"` (a venturi, `A_star=0.9`: 10 % below the throat area; x on 401 stations so that the throat is a grid point), `"transonic"`, `"supersonic"`. *see / read / change.*
12. `nb.note` — **N45 [B]** "(ours; result given, not derived) the Mach-number twin: $\dfrac{dM}{M}=-\dfrac{1+\frac{\gamma-1}{2}M^2}{1-M^2}\dfrac{dA}{A}$ (ours)." + `nb.code`: sympy —
    `sp.simplify(sp.diff(sp.log(A_over_Astar), M)*M + (1 - M**2)/(1 + (g-1)/2*M**2))` → 0.
13. `nb.worked_example("1 % more area", "$\dfrac{du}{u}=-\dfrac{1}{1-M^2}\dfrac{dA}{A}$ (15.34) with dA/A = +0.01. 1. M = 0.5: du/u = −0.01/(1 − 0.25) = **−1.33 %** (slower; a
    diffuser). 2. M = 2: du/u = −0.01/(1 − 4) = **+0.33 %** (faster; still a nozzle). 3. Density by $d\rho/\rho=-M^2\,du/u$ (15.33): at M = 2, dρ/ρ = −4 ×
    0.33 % = −1.33 %. 4. Check mass $\dfrac{d\rho}{\rho}+\dfrac{du}{u}+\dfrac{dA}{A}=0$ (15.19): −1.33 + 0.33 + 1.00 = 0 ✓.")`
14. `nb.code` — `GD.area_velocity(np.array([0.3, 0.5, 2.0]), 1.0)` (du/u per unit dA/A); `GD.influence_coefficients(0.5, dA_over_A=0.01)`. *expect:* −1.099,
    −1.333, +0.333; at M = 0.5: du/u = −1.33 %, dρ/ρ = +0.33 %, dM/M = −1.40 %. *explain:* 4 points.
15. `nb.check_agree` — **from scratch:** along `sol = GD.nozzle_isentropic(x, GD.laval_area, branch="transonic")` compute `u ∝
    M*np.sqrt(T_over_T0)`; `slope = np.gradient(np.log(u), x)/np.gradient(np.log(sol["A"]), x)` away from the throat; `assert
    np.allclose(slope[mask], GD.area_velocity(sol["M"][mask], 1.0), rtol=2e-2)`.
16. `nb.plotly` — **F4** `slider_figure`: exit-to-throat area ratio (1.2 … 4, 15 steps) → M(x) for the three solutions on
    `GD.laval_area` (subsonic, transonic, supersonic), the line M = 1 dashed. *see / read / change.*
17. `nb.explainer("area_mach_nozzle", heading="Why does a supersonic nozzle widen?", why="you move a cursor along a duct you
    shaped and three bars — density, speed, area — always sum to zero while the density bar overtakes the speed bar exactly at
    Mach 1", tries=["Play the station through the throat and watch the density bar overtake the speed bar.", "Switch the
    solution to 'subsonic throughout': the same throat, no Mach 1.", "Find the ghost dot: the other Mach number with the same
    area."])`
18. `nb.note` — **N46 [B]** worked example with **our own inputs** (`rocket`): `r = ch15.rocket_nozzle(rocket["thrust"],
    rocket["p0"], rocket["T0"], rocket["p_exit"], air)`. Chain: $p_e/p_0$ = 0.0338 → $M_e$ (inverse of $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29)) → $A_e/A^*$ ($\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31)),
    $T_e$ → $u_e$ → ṁ = thrust/$u_e$ → $A_e$ → $A^*$. *expect:* (keys `M_e`, `T_e`, `u_e`, `mdot`, `A_e`, `Ae_over_At`, `A_throat`) $M_e$ = 2.857, $T_e$ = 1140 K, $u_e$ = 1933 m/s, ṁ = 25.86 kg/s,
    $A_e$ = 0.0432 m², $A_e/A^*$ = 3.696, $A^*$ = 0.0117 m² (check: `GD.choked_mass_flow` with that throat returns 25.86 kg/s).
    "> ⚠️ Common confusion (trap T9): the chamber *is* the stagnation state, and the nozzle is designed so that the exit
    pressure equals the ambient one." Ambient pressure from `STAT` (sea level); if its lapse rate is printed: Kundu
    $\Gamma\equiv dT/dz=-6.5$ K/km, meteorological $-dT/dz=+6.5$ K/km.
19. `nb.note` — **N47 [B]** "Thrust of a pressure-matched nozzle: $\rho_eu_e^2A_e=\dot mu_e=\gamma p_eM_e^2A_e$. Off design a
    pressure term appears: thrust = $\dot mu_e+(p_e-p_a)A_e$." `ch15.nozzle_thrust(r["mdot"], r["u_e"], 101325.0, 26500.0, r["A_e"])`. *expect:* 50.0 kN on design at sea level; 53.2 kN at 10 km (ambient 26.5 kPa).
20. `nb.md` — **What would change if…** "…the supersonic stream must come back to a high pressure? It cannot do so smoothly
    in general. It jumps — C07."

### A.6 §15.6 Normal Shock Waves — C07, C08
#### C07 — The normal-shock jump relations (15.39)
1. `nb.section("15.6", "Normal Shock Waves", intro="**What is this section about?** A supersonic stream can change to a
   subsonic one across a layer a few molecular free paths thick. Across it the three duct laws still hold — and they fix
   everything from one number, the upstream Mach number.")`
2. `nb.core("C07", "The jump across a normal shock: $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39)", question="If a flow changes abruptly with no heat, no wall
   force and no change of area, what can the state behind it be?")`
3. `nb.md` — **Plain words:** "Cars on a motorway meet a jam: within a few car lengths they go from fast and far apart to
   slow and close together. No car is created or destroyed, and the 'front' of the jam is thin. A shock is the same thing in
   a gas: fast, thin gas arrives; slow, dense, hot gas leaves. We do not need to know what happens inside the thin layer —
   only that mass, momentum and energy that go in must come out."
4. `nb.md` — **The idea** (ASCII): `1: u₁ fast, ρ₁, p₁, T₁ ──►│shock│──► 2: u₂ slow, ρ₂ > ρ₁, p₂ > p₁, T₂ > T₁`; "three laws,
   three unknowns (ρ₂, u₂, p₂) → two solutions: nothing happens, or the shock".
5. `nb.note` — **N48 [B]** "A shock is a compression of finite strength, a few mean free paths thick, treated as a jump.
   Entropy is made inside it, so the isentropic ratios of C04–C05 cannot be carried across — only the three laws of C03."
6. `nb.note` — **N49 [B]** mass: $\rho_1u_1=\rho_2u_2$ (15.35) ($\rho uA=\dot m=\text{const.}$ (15.18) at constant area) · **N50 [B]** momentum: $p_1-p_2=-\rho_1u_1^2+\rho_2u_2^2$ (15.36) ($-\dot mu_1+\dot mu_2=(pA)_1-(pA)_2+F$ (15.24) with F = 0, $A_1=A_2$) ·
   **N51 [B]** energy: $h_1+\tfrac12u_1^2=h_2+\tfrac12u_2^2$ (15.37) ($h_2+\tfrac12u_2^2-h_1-\tfrac12u_1^2=Q$ (15.22) with Q = 0). "> ⚠️ Common confusion (trap T17): the energy line looks like Bernoulli's
   equation but is only the energy law; entropy is not constant across the shock."
7. `nb.primer("dividing a known root out of a quadratic", "If you already know one solution $r_0$ of a quadratic, the
   quadratic is $(r-r_0)$ times something linear — divide and the other root drops out. (Vieta: the product of the roots is
   the constant term over the leading coefficient.) Here the known root is 'nothing happens'.", code="import sympy as sp\nr =
   sp.symbols('r')\nquad = 2*r**2 - 7*r + 5                 # we know r = 1 is a root: 2 - 7 + 5 = 0\nprint(sp.cancel(quad/(r -
   1)))           # 2*r - 5  ->  the other root is 5/2 (and 1 x 5/2 = 5/2, Vieta)")` (**P363**; reminder of P71)
8. `nb.derivation("D13", …)` — Part F D13 (13 steps, ★★★, sympy check), ref "15.39": holds **N52 [B]** (step 1) and **N53 [B]**
   (step 12, boxed half-way line).
9. `nb.derivation("D14", …)` — Part F D14 (7 steps), ref "15.40": holds **N54 [B]**. "> ⚠️ Common confusion (trap T6): $M_2^2=\dfrac{(\gamma-1)M_1^2+2}{2\gamma M_1^2+1-\gamma}$ (15.40)
   is its own inverse — feed it the downstream Mach number and it returns the upstream one without complaint.
   `GD.normal_shock` refuses $M_1<1$ unless asked. Other sources number the stations 0 and 1."
10. `nb.note` — **N55 [B]** `nb.figure`: $M_2$ against $M_1$ from `GD.normal_shock_mach` with the no-change diagonal, the half
    $M_1<1$ shaded "forbidden (reason: C08)", the asymptote $\sqrt{(\gamma-1)/2\gamma}$ = 0.378. *see / read* ("apply the map
    twice and you are back: the curve is its own mirror image in the diagonal") */ change.*
11. `nb.derivation("D15", …)` — Part F D15 (7 steps), ref "15.41": holds **N56 [B]** (steps 1–3) and **N57 [B]** (steps 4–7).
    "As $M_1\to\infty$ the density ratio tends to $(\gamma+1)/(\gamma-1)$ = 6: no shock squeezes air more than sixfold."
12. `nb.worked_example("a Mach 2 shock by hand", "γ = 1.4, $M_1$ = 2. 1. $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39): $p_2/p_1$ = 1 + (2.8/2.4)(4 − 1) = **4.5**. 2.
    $\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}$ (15.41): $\\rho_2/\\rho_1$ = 2.4 × 4/(0.4 × 4 + 2) = 9.6/3.6 = **2.667**. 3. Gas law: $T_2/T_1$ = 4.5/2.667 = **1.6875** ($\dfrac{T_2}{T_1}=1+\dfrac{2(\gamma-1)}{(\gamma+1)^2}\dfrac{\gamma M_1^2+1}{M_1^2}\left(M_1^2-1\right)$ (15.42)
    gives the same). 4. $M_2^2=\dfrac{(\gamma-1)M_1^2+2}{2\gamma M_1^2+1-\gamma}$ (15.40): $M_2^2$ = (1.6 + 2)/(11.2 − 0.4) = 3.6/10.8 = 1/3, $M_2$ = **0.577**. 5. In sea-level air: 681 m/s
    in, 255 m/s out; 288 K in, 486 K out. 6. $T_0$: 288 × 1.8 = 519 K before; 486 × (1 + 0.2/3) = 519 K after ✓.")`
13. `nb.code` — `s = GD.normal_shock(2.0)`; `GD.normal_shock_residuals(state1, state2)` with `state1 = dict(p=1.0, rho=1.4, u=2.0)` (units $p_1$, $c_1$) and `state2 = dict(p=s["p_ratio"], rho=1.4*s["rho_ratio"], u=2.0*s["u_ratio"])`; `GD.normal_shock_fluxes(2.0)`; `GD.normal_shock(0.7)` in
    a `try` (raises); `GD.normal_shock(0.7, allow_expansion=True)["admissible"]`; `ch15.normal_shock_sympy()`. *expect:* `M2` 0.5774, `p_ratio` 4.5, `rho_ratio` 2.6667, `T_ratio` 1.6875, `p0_ratio` 0.7209, `dS_over_cv` 0.1309; residuals < 10⁻¹³; ValueError; False. *explain:* 5 points.
14. `nb.check_agree` — **from scratch:** unknown r = $u_2/u_1$: `energy = lambda r: g/(g-1)*a + 0.5 - g/(g-1)*r*(a + 1 - r) -
    0.5*r*r` with `a = 1/(g*M1**2)`; `r = brentq(energy, 0.05, 0.999)` (the bracket excludes the root r = 1); `p_ratio = 1 + (1 -
    r)/a`; `assert np.allclose([1/r, p_ratio], [s["rho_ratio"], s["p_ratio"]])`.
15. `nb.plotly` — **F5a** `slider_figure`: γ (1.1 … 1.67, 20 steps) → $M_2$, $p_2/p_1$ (÷ 10), $\rho_2/\rho_1$, $T_2/T_1$ against
    $M_1$ (0.6 … 5), the half $M_1<1$ shaded, the density limit $(\gamma+1)/(\gamma-1)$ dashed. *see / read / change.*
16. `nb.note` — **N64 [B]** "**The same shock, seen from the still gas** (result given, not derived): subtract $u_1$
    everywhere. The shock runs into still gas at $u_{\text{shock}}=M_1c_1$ and drags the gas behind it at
    $u_{\text{gas}}=u_1-u_2$." Tiny example: "a shock that doubles the pressure of sea-level air: $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39) gives $M_1^2$ = 1 +
    (2.4/2.8) × 1 = 1.857, $M_1$ = 1.363; shock speed 1.363 × 340.3 = 464 m/s; wind behind it 464 × (1 − 1/1.625) = 178 m/s."
    `GD.moving_shock(101325.0, 288.15, p2_over_p1=2.0, gas=air)` + `nb.figure`: `draw_shock_two_frames(ms)` (arrows to scale in
    both frames). "> ⚠️ Common confusion (trap T7)" two-column table: same in both frames (p, ρ, T, c, S) · different (u, M,
    $T_0$, $p_0$). *expect:* `M_s` 1.3628, `u_shock` 463.7 m/s, `u_gas` 178.4 m/s, `T2` 354.6 K.
17. **(re-enter C07)** `nb.recap("R17", "Strong shocks and the blast wave", "Ch. 1 found $R(t)\\propto(E/\\rho_1)^{1/5}t^{2/5}$ by
    dimensional analysis. **New here:** for $M_1^2\\gg1$, $\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}$ (15.41) tends to $(\\gamma+1)/(\\gamma-1)$ and $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39) to
    $2\\gamma M_1^2/(\\gamma+1)$ (results given, not derived); a deliberately **crude** energy argument then gives a constant.
    > ⚠️ (trap T22) crude: all the energy as uniform internal energy behind the front; the constant quoted in Ch. 1 uses a
    different normalisation.", where="Ch. 1 Example 1.4 (C69)")` + `GD.strong_shock_limits()`, `ch15.blast_constant_crude()`,
    `ch01.blast_radius` (figure title carries "crude"), `sp.limit` on $\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}$ (15.41). *expect:* `rho_ratio` 6.0, `M2` 0.378; K = 1.07 (hemisphere).
18. `nb.md` — **What would change if…** "…$M_1$ were below 1? The three laws still balance — step 8 of D13 never asked. What
    forbids it is the second law. C08."

#### C08 — The entropy jump and the second law (15.43)
1. `nb.core("C08", "Why a shock only compresses: $\dfrac{S_2-S_1}{C_v}=\ln\left[\dfrac{p_2}{p_1}\left(\dfrac{\rho_1}{\rho_2}\right)^{\gamma}\right]$ (15.43)", question="Mass, momentum and energy allow a jump in either
   direction. Which law chooses — and what exactly is lost?")`
2. `nb.md` — **Plain words:** "Film a traffic jam and play it backwards: slow dense cars suddenly spread out and speed up,
   with no reason. Nothing in the car count forbids it — it just never happens. A shock run backwards (subsonic gas
   jumping to supersonic) conserves mass, momentum and energy too. What it would do is *destroy* entropy, and the second law
   forbids that. The same bookkeeping tells us what a real shock costs: not energy — that is conserved — but stagnation
   pressure, the ability to do work."
3. `nb.md` — **The idea** (two-column table): "conserved across a shock: ṁ, momentum flux + p, $h_0$, $T_0$ (amber, flat)" ·
   "not conserved: S rises (purple), $p_0$ falls, $A^*$ grows".
4. `nb.md` — 🔁 "entropy is a state function; the first law and the Gibbs relation $T\,dS=dh-dp/\rho$: Ch. 1 (C33–C35)."
5. `nb.derivation("D16", …)` — Part F D16 (5 steps), ref "15.43".
6. `nb.note` — **N59 [B]** `nb.figure`: `GD.normal_shock_entropy` against $M_1$ (0.5 … 3) with the cubic $\dfrac{S_2-S_1}{C_v}\cong\dfrac{2\gamma(\gamma-1)}{3(\gamma+1)^2}\left(M_1^2-1\right)^3$ (15.44a) overlaid
   (dashed) and a log–log inset of ΔS against $M_1^2-1$ whose slope is 3. *see / read / change.*
7. `nb.note` — **N61 [B]** "Reading $\dfrac{S_2-S_1}{C_v}=\ln\left[\dfrac{p_2}{p_1}\left(\dfrac{\rho_1}{\rho_2}\right)^{\gamma}\right]$ (15.43): $S_2-S_1>0\iff M_1>1$. So $M_1>1$, $M_2<1$, compression only (for a perfect gas;
   expansion shocks need another equation of state — N90 in C12). Weak shocks are nearly isentropic, which is why even loud
   sound is still acoustics (C01)." `GD.shock_admissible(2.0)`, `GD.shock_admissible(0.7)`.
8. `nb.primer("the series of ln(1 + x) to third order", "$\\ln(1+x)=x-\\tfrac12x^2+\\tfrac13x^3-\\dots$ for small x. When
   several logarithms are added, their x and x² terms can cancel *exactly*, and the answer is then of third order — you
   must keep every x³.", code="import sympy as sp\nx = sp.symbols('x')\nprint(sp.series(sp.log(1 + x), x, 0, 4))                    #
   x - x**2/2 + x**3/3 + O(x**4)\nprint(sp.series(sp.log(1 + 2*x) - 2*sp.log(1 + x), x, 0, 4))   # -x**2 + 2*x**3 + ...: first
   order cancelled")` (**P364**; reminder of P26, P68, P117)
9. `nb.derivation("D17", …)` — Part F D17 (11 steps, ★★★, sympy check), ref "15.44": holds **N60 [B]**. + `nb.code`:
   `GD.weak_shock_entropy(0.1, form="pressure")`, `GD.normal_shock_entropy` at the same shock, `ch15.weak_shock_series()`. *expect:*
   a 10 % pressure jump: cubic 4.1 × 10⁻⁵, exact 3.5 × 10⁻⁵ (in units of $C_v$); it loses 8.8 × 10⁻⁵ of $p_0$ — nine parts in
   a hundred thousand. A small table "exact against cubic": $M_1$ = 1.05: 5.9 × 10⁻⁵ / 7.0 × 10⁻⁵ · 1.1: 4.2 × 10⁻⁴ / 6.0 × 10⁻⁴ ·
   1.2: 0.0029 / 0.0055 · 1.5: 0.029 / 0.127 (the cubic is a *local* statement).
10. `nb.derivation("D18", …)` — Part F D18 (7 steps), ref "ours": holds **N62 [B]**. "> ⚠️ Common confusion (trap T8, second
    half): $A^*$ is not conserved across a shock — it grows by $p_{01}/p_{02}$."
11. `nb.note` — **N58 [B]** "The relations $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39)–$\dfrac{T_2}{T_1}=1+\dfrac{2(\gamma-1)}{(\gamma+1)^2}\dfrac{\gamma M_1^2+1}{M_1^2}\left(M_1^2-1\right)$ (15.42) are the **Rankine–Hugoniot relations**. (Ours:) eliminate $M_1$
    between $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39) and $\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}$ (15.41) to get the Hugoniot curve $\dfrac{p_2}{p_1}=\dfrac{(\gamma+1)\rho_2/\rho_1-(\gamma-1)}{(\gamma+1)-(\gamma-1)\rho_2/\rho_1}$ (ours)." + `nb.figure`: `GD.hugoniot_pressure_ratio` beside the isentrope
    $p\propto\rho^\gamma$ on log axes — they touch to third order at (1, 1); a sympy cell checks the elimination. *expect:* at
    $\rho_2/\rho_1$ = 2.667 the Hugoniot gives 4.50, the isentrope 3.95. *see / read / change.*
12. `nb.note` — **N63 [B]** "The book's shock table is not reproduced; `GD.normal_shock_table(np.arange(1.0, 5.01, 0.5))` generates ours (reminder: a
    DataFrame, P361 in C05), checked against the public NACA Report 1135 rows." *expect:* the row $M_1$ = 2: 0.5774, 4.5,
    1.6875, 0.7209.
13. **(re-enter C08)** `nb.recap("R16", "The hydraulic jump", "$H_2/H_1=\tfrac12\left(-1+\sqrt{1+8Fr_1^2}\right)$ (7.81) — supercritical ($Fr_1>1$) water jumps to subcritical.",
    where="Ch. 7 §7.6 (C11–C12)")` + a side-by-side table (Froude number ↔ Mach number · supercritical → subcritical ↔
    supersonic → subsonic · depth ↔ density · energy lost ↔ stagnation pressure lost) + `ch07.hydraulic_jump(1.0, Fr1=2.0)`.
    *expect:* $H_2/H_1$ = 2.372. "> ⚠️ Common confusion (trap T21): the analogy is partial — the jump loses mechanical energy to
    heat; the shock conserves energy and loses $p_0$."
14. `nb.worked_example("what a Mach 2 shock costs", "From C07: $p_2/p_1$ = 4.5, $\\rho_2/\\rho_1$ = 2.667. 1. $\dfrac{S_2-S_1}{C_v}=\ln\left[\dfrac{p_2}{p_1}\left(\dfrac{\rho_1}{\rho_2}\right)^{\gamma}\right]$ (15.43):
    $(S_2-S_1)/C_v=\\ln[4.5\\times(1/2.667)^{1.4}]=\\ln(4.5\\times0.2533)=\\ln1.140$ = **0.131**. 2. In units of R: 0.131/(γ − 1) =
    0.327. 3. Stagnation pressure kept: $e^{-0.327}$ = **0.721** — 28 % lost. 4. Cross-check with $\dfrac{p_{02}}{p_{01}}=\left[\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}\right]^{\gamma/(\gamma-1)}\left[\dfrac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^{1/(\gamma-1)}$ (ours): (2.667)^3.5 × (2.4/10.8)^2.5
    = 30.97 × 0.02328 = 0.721 ✓. 5. $A^*$ behind the shock is 1/0.721 = 1.387 times larger.")`
15. `nb.code` — `GD.normal_shock_entropy(2.0)`, `GD.normal_shock_stagnation_pressure_ratio(2.0)`, `np.exp(-0.13091/0.4)`. *expect:*
    0.1309, 0.7209, 0.7209. *explain:* 3 points.
16. `nb.check_agree` — **from scratch:** `dS = TH.perfect_gas_entropy_change_p(T1, p1, T1*s["T_ratio"], p1*s["p_ratio"],
    cp=air.cp, R=air.R)`; `p0_ratio = s["p_ratio"]*GD.stagnation_pressure_ratio(s["M2"])/GD.stagnation_pressure_ratio(2.0)`;
    `assert np.allclose(dS/air.cv, GD.normal_shock_entropy(2.0))`; `assert np.allclose([p0_ratio, np.exp(-dS/air.R)],
    GD.normal_shock_stagnation_pressure_ratio(2.0))` — "three routes, one number".
17. `nb.plotly` — **F5b** `slider_figure`: γ (1.1 … 1.67) → $(S_2-S_1)/C_v$ and $p_{02}/p_{01}$ against $M_1$ (0.6 … 5); the
    forbidden half shaded; $T_{02}/T_{01}$ = 1 flat in amber. *see / read / change.*
18. `nb.explainer("normal_shock_jump", heading="Why can a shock only compress?", why="you drag the upstream Mach number below 1:
    the three flux bars stay perfectly balanced while the entropy bar goes negative — you find the second law's veto by
    experiment", tries=["Drag $M_1$ from 2 down to 0.7: which bar breaks?", "Sit at $M_1$ = 1.1: how much stagnation pressure is
    lost?", "Switch to the frame of the still gas: which numbers changed, which did not?"])`
19. **(re-enter C08)** `nb.recap("R18", "Viscosity sets the thickness, not the loss", "In turbulence the dissipation rate is
    set by the large scales, not by the viscosity. Likewise the entropy rise of a shock is fixed by the jump conditions;
    viscosity and conduction only set how thick the layer is (N65, N69 below).", where="Ch. 12 §12.7 (C07)")`
20. `nb.note` — **N65 [B]** "Thickness estimate: a Reynolds number of order one on the velocity jump,
    $(u_1-u_2)\,\delta_{\text{sh}}/\nu_{\text{visc}}\sim1$ (δ is a thickness here, a deflection angle in §15.9 — trap T5)."
    `ch15.shock_thickness_estimate(u1 - u2, 1.46e-5)` (with $u_1-u_2$ = 425 m/s at $M_1$ = 2) ÷ `ch01.mean_free_path_air(288.15)`. *expect:*
    3.4 × 10⁻⁸ m at $M_1$ = 2 and 1.4 × 10⁻⁷ m at 1.2 — of the order of the mean free path (≈ 7 × 10⁻⁸ m): the continuum model
    is at its limit (🔁 Knudsen number, Ch. 1 C06). Our numbers, not the book's.
21. **(re-enter C08)** `nb.recap("R19", "The equations inside a shock", "Steady, one-dimensional mass, momentum and energy with
    viscosity and conduction: $\\dfrac{d(\\rho u)}{dx}=0$, $\\rho u\\dfrac{du}{dx}+\\dfrac{dp}{dx}=\\dfrac{d}{dx}\\left(\\mu''
    \\dfrac{du}{dx}\\right)$, $\\rho u\\dfrac{dh}{dx}-u\\dfrac{dp}{dx}=\\mu''\\left(\\dfrac{du}{dx}\\right)^2+\\dfrac{d}{dx}\\left(k
    \\dfrac{dT}{dx}\\right)$ with $\\mu''=\\tfrac43\\mu+\\mu_v$ (reductions of the Ch. 4 equations).", where="Ch. 4 (C08, C10; N51
    bulk viscosity)")` + `ch15.shock_structure_sympy()`.
22. `nb.note` — **N66 [B]** "Each integrates once (result given, not derived): $\rho u=m$, $mu+p=\mu''\dfrac{du}{dx}+mV$,
    $m\left(h+\tfrac12u^2\right)=\mu''u\dfrac{du}{dx}+k\dfrac{dT}{dx}+mI$. Far upstream and downstream the gradients vanish and
    these *are* $\rho_1u_1=\rho_2u_2$ (15.35), $p_1-p_2=-\rho_1u_1^2+\rho_2u_2^2$ (15.36), $h_1+\tfrac12u_1^2=h_2+\tfrac12u_2^2$ (15.37): the jump conditions do not care what happens inside." `ch15.shock_structure_constants(2.0, 288.15, 101325.0, air)` (keys `m`, `V`, `I`, `U1`, `U2`);
    a sympy cell differentiates them back. · **N67 [B]** "With $Pr''\equiv\mu''C_p/k=1$ the total enthalpy $C_pT+\tfrac12u^2=I$
    at every point inside the shock (Becker; result given). Trap T20: $Pr''$ is built on $\mu''$ — about 0.95 for air with
    $\mu_v=0$." · **N68 [B]** "One first-order equation remains for $U=u/V$; its two rest points are the end states, (ours)
    $U_1=\gamma M_1^2/(\gamma M_1^2+1)$ and the same with $M_2$; $U_1+U_2=2\gamma/(\gamma+1)$. (Slip #7: the book factorises
    before dividing by $V^2$.)" `ch15.becker_shock_roots(2.0)`. *expect:* (28/33, 7/22) = (0.8485, 0.3182); sum 1.1667.
23. `nb.note` — **N69 [B]** `nb.figure`: the profile $\dfrac{U_1\ln(U_1-U)-U_2\ln(U-U_2)}{U_1-U_2}=\dfrac{\gamma+1}{2\gamma}\,\eta$ (unnumbered, §15.6) inverted by `ch15.becker_shock_profile(eta, M1)` for $M_1$ = 1.2, 2,
    4, cross-checked by `solve_ivp` on `ch15.becker_shock_ode` (dots on the curves; reminder P94), and
    `ch15.becker_shock_thickness` against $M_1$ in a second panel (it grows without bound as $M_1\to1$). "The origin of η is
    arbitrary (trap T20)." *see / read / change.* · **N70 [B]** "Our figure uses our own Mach numbers; for γ = 1.4, $M_1$ = 2 the
    plateaus are 28/33 and 7/22 (slip #8 concerns a caption of the book's figure; the comparison is a private test)." · **N71
    [C]** "For $Pr''\to0$ a second, thinner shock sits inside the first; see the literature the book cites — not developed
    in the book or here."
24. `nb.md` — **What would change if…** "…you put this jump inside a nozzle? Then the back pressure decides where it
    stands. First the simpler nozzle — C09."


### A.7 §15.7 Operation of Nozzles at Different Back Pressures — C09, C10
#### C09 — Choking of a convergent nozzle
1. `nb.section("15.7", "Operation of Nozzles at Different Back Pressures", intro="**What is this section about?** A reservoir, a
   nozzle and a tank whose pressure you lower. The book describes what happens in words and figures; the section has no
   numbered equation, so **every formula here is ours**, built from the earlier sections.")`
2. `nb.core("C09", "Choking: $\dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)}$ (ours)", question="If you keep lowering the pressure downstream of a nozzle, does more and
   more gas flow?")`
3. `nb.md` — **Plain words:** "Let the air out of a bicycle tyre. Open the valve a little into a room at lower pressure
   and air flows; lower the room pressure and more flows. Intuition says: lower it to vacuum and the flow keeps growing.
   It does not. Once the pressure outside is about half the pressure inside, the air at the exit moves at exactly the
   speed of sound — and news of any further drop outside cannot travel upstream against it. The nozzle is **choked**."
4. `nb.md` — **The idea** (ASCII): `reservoir p₀, T₀ ══⟩ exit │ tank p_B`; "exit subsonic: exit pressure = $p_B$ (news gets in)"
   / "exit sonic: exit pressure = $p^*$ ≥ $p_B$ (news cannot get in)".
5. `nb.note` — **N72 [B]** three pressures: supply (stagnation) $p_0$, back pressure $p_B$, exit-plane pressure. "**The exit
   pressure equals the back pressure if and only if the exit flow is subsonic** — the range of influence of C02 (N20) put to
   work." + `nb.figure`: `draw_reservoir_nozzle_tank()`. *see / read / change.*
6. `nb.derivation("D19", …)` — Part F D19 (7 steps), ref "ours".
7. `nb.worked_example("our nozzle, three back pressures", "$p_0$ = 500 kPa, $T_0$ = 400 K, exit area 10 cm², γ = 1.4; critical
   ratio $\dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)}$ (ours) = 0.528. 1. $p_B$ = 400 kPa (ratio 0.8 > 0.528: not choked): invert $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29): $M_e^2$ = 5 × (0.8^(−0.2857) − 1) =
   5 × 0.0658 = 0.329, $M_e$ = 0.574. 2. Flux function $\dfrac{\dot m\sqrt{RT_0}}{p_0A}=\sqrt\gamma\,M\left(1+\dfrac{\gamma-1}{2}M^2\right)^{-(\gamma+1)/(2(\gamma-1))}$ (ours) at 0.574: 1.183 × 0.574/1.0658³ = 0.561 — that is 82 % of its
   maximum 0.685. 3. ṁ = 0.561 × 500 000 × 0.001/√(287 × 400) = 0.83 kg/s. 4. $p_B$ = 264 kPa (ratio 0.528): $M_e$ = 1, ṁ = 1.01
   kg/s. 5. $p_B$ = 100 kPa: still $M_e$ = 1 and ṁ = 1.01 kg/s; the exit pressure stays at 264 kPa and the jet finishes
   expanding outside.")`
8. `nb.code` — `GD.critical_pressure_ratio()`; `[GD.convergent_nozzle(noz["p0"], noz["T0"], pB, noz["A_t"], air) for pB in (4.75e5,
   4.0e5, 3.0e5, 2.6414e5, 1.0e5)]` as a DataFrame; `GD.mass_flow_vs_back_pressure(0.8)`. *expect:* 0.5283; `M_exit` 0.272, 0.574,
   0.886, 1.000, 1.000; `mdot` 0.454, 0.827, 0.999, 1.010, 1.010 kg/s; `choked` F, F, F, T, T; `mdot` 0.5607 (dimensionless) with `mdot_over_choked` 0.819. *explain:* 4 points.
9. `nb.check_agree` — **from scratch:** a loop over 50 back pressures with an `if pB/p0 > crit:` (subsonic exit: M from the
   pressure ratio) `else:` (M = 1, exit pressure = crit × p0); `mdot = sqrt(g)*M*(1 + 0.2*M*M)**-3*p0*A/np.sqrt(R*T0)`; `assert
   np.allclose(mdot_mine, [GD.convergent_nozzle(p0, T0, pB, A, air)["mdot"] for pB in pBs])`.
10. `nb.note` — **N73 [B]** `nb.figure`: $p(x)/p_0$ along a convergent nozzle (`GD.nozzle_isentropic` on the converging half
    of `GD.laval_area`, scaled so that the exit Mach number is that of `GD.convergent_nozzle`) for $p_B/p_0$ = 0.95, 0.8, 0.6,
    0.528, 0.3, 0.1 — the last three curves lie on top of one another; beside it ṁ against $p_B/p_0$ with its plateau. *see /
    read / change.*
11. `nb.note` — **N74 [B]** "**Choking:** with the exit sonic the flow upstream no longer responds to $p_B$, and
    $\dot m=\dot m_{\max}$ of $\dot m_{\max}=p_0A^*\sqrt{\dfrac{\gamma}{RT_0}}\left(\dfrac{2}{\gamma+1}\right)^{(\gamma+1)/(2(\gamma-1))}$ (ours). > ⚠️ Common confusion (trap T10): 'choked' does not mean the flow cannot be raised — it
    is the most the throat passes *for this supply*. Raise $p_0$ or lower $T_0$ and it rises ($\dot m\propto p_0/\sqrt{T_0}$)."
    + `nb.plotly` — **F6** two `slider_figure`s: (a) $p_0$ (200 … 1000 kPa, 17 steps) and (b) $T_0$ (300 … 1200 K, 19 steps) → ṁ
    [kg/s] against $p_B$ [kPa] (dimensional axes, so the plateau and the kink visibly move). *see / read / change.*
12. `nb.md` — **What would change if…** "…you added a diverging part after the throat? Then below the choking pressure the
    flow has somewhere to go: supersonic. And the back pressure gets a new job. C10."

#### C10 — The convergent–divergent (Laval) nozzle
1. `nb.core("C10", "The Laval nozzle: three dividing back pressures $\dfrac{p_b}{p_0}=\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sub}}(A_e/A_t)},\ \dfrac{p_f}{p_0}=\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sup}}(A_e/A_t)},\ \dfrac{p_d}{p_0}=\dfrac{p_f}{p_0}\left.\dfrac{p_2}{p_1}\right\rvert_{M_{\text{sup}}}$ (ours)", question="Once the throat is sonic, what is
   left for the back pressure to decide?")`
2. `nb.md` — **Plain words:** "A rocket nozzle on the test stand at sea level, then the same nozzle at altitude. The
   chamber is the same, the throat is the same, the mass flow is the same — but the exhaust looks different, and at
   start-up a sharp front can be seen standing inside the bell. With the throat sonic the back pressure cannot change the
   mass flow. All it can do is place a normal shock somewhere in the diverging part — or, when it is low enough, leave the
   inside of the nozzle alone and reshape the jet outside."
3. `nb.md` — **The idea**: the ladder of seven conditions as a table (name · range of $p_B/p_0$ for $A_e/A_t$ = 2 · what you
   see): subsonic throughout (> 0.937) · just choked (= 0.937 = $p_b$) · shock in the diverging part (0.513 … 0.937) · shock at
   the exit (= 0.513 = $p_d$) · over-expanded, oblique shocks outside (0.094 … 0.513) · design (= 0.094 = $p_f$) ·
   under-expanded, expansion fans outside (< 0.094).
4. `nb.primer("nested root finding", "Sometimes the function you want a root of itself needs a root-finder to evaluate.
   Here: 'exit pressure for a shock at Mach $M_1$' needs the subsonic inverse of the area ratio; then we search $M_1$ so
   that this exit pressure equals the back pressure. One `brentq` inside another.", code="from scipy.optimize import brentq\ninner
   = lambda y: brentq(lambda x: x**3 + x - y, 0.0, 10.0)   # inner: solve x^3 + x = y for x\nouter = brentq(lambda y: inner(y) - 1.5,
   0.0, 100.0)       # outer: which y makes that x equal 1.5?\nprint(outer, 1.5**3 + 1.5)                               # 4.875
   both")` (**P365**; reminder of P108, P360)
5. `nb.derivation("D20", …)` — Part F D20 (9 steps), ref "ours": holds **N75 [B]** (steps 5–9).
6. `nb.note` — **N78 [B]** `nb.worked_example("a shock where the area is 1.5 throats (our inputs)", "$A_e/A_t$ = 2, $T_0$ = 400 K, a
   normal shock at $A/A_t$ = 1.5. 1. Ahead of the shock the flow is isentropic and supersonic with $A^*=A_t$: invert $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31) →
   $M_1$ = 1.854. 2. Across the shock $M_2^2=\dfrac{(\gamma-1)M_1^2+2}{2\gamma M_1^2+1-\gamma}$ (15.40): $M_2$ = 0.605; $\dfrac{p_{02}}{p_{01}}=\left[\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}\right]^{\gamma/(\gamma-1)}\left[\dfrac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^{1/(\gamma-1)}$ (ours): $p_{02}/p_{01}$ = 0.788. 3. **The sonic area grows:**
   $A_2^*=A_t/0.788=1.268A_t$ (trap T8). 4. Exit: $A_e/A_2^*$ = 2/1.268 = 1.577 → **subsonic** inverse → $M_e$ = 0.404. 5.
   $\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28): $T_e/T_0$ = 0.968 → $T_e$ = 387 K ($T_0$ is unchanged by the shock). 6. Exit pressure: $p_e/p_{01}$ = 0.788 × 0.894 =
   0.704 — and since the exit is subsonic, that *is* the back pressure this shock position belongs to.")` + `M1 = GD.mach_from_area_ratio(1.5, branch="supersonic")`; `ch15.nozzle_with_shock(M1, noz["A_t"], 2*noz["A_t"], noz["p0"], noz["T0"], air)`. *expect:* M1 1.8541; `A_shock` 1.500 × 10⁻³ m², `M2` 0.6048, `p02` 394.2 kPa (0.7884 $p_{01}$), `A_star2` 1.268 × 10⁻³ m², `M3` 0.4042, `T3` 387.3 K, `p3` 352.2 kPa (0.7045 $p_{01}$), `mdot` 1.010 kg/s. (Slip #6 concerns the book's own example; ours uses our notation.)
7. `nb.code` — `th = GD.laval_critical_pressures(2.0)`; `[GD.laval_regime(pB, 2.0) for pB in (0.97, th["p_b"], 0.7, th["p_d"], 0.3, th["p_f"], 0.05)]`;
   `GD.laval_shock_location(0.7, 2.0)`; `GD.laval_exit_state(0.7, 2.0)`; `sol = GD.laval_nozzle(x, GD.laval_area, 0.7)` on 401 stations (the throat a grid point). *expect:*
   `p_b` 0.9372, `p_d` 0.5134, `p_f` 0.0939; regimes "subsonic", "sonic throat", "shock in nozzle", "shock at exit", "overexpanded", "design", "underexpanded" (the three equalities need the exact thresholds);
   for 0.7: `M1` 1.8627, `A_shock_over_At` 1.5101, `p02_over_p01` 0.7845, `M_exit` 0.4067, `mdot_over_choked` 1.0 (for 0.9: `M1`
   1.369, area 1.099). *explain:* 5 points ("…feeding 0.7045 back returns the worked example's shock at 1.5").
8. `nb.check_agree` — **from scratch:** `def exit_p(M1): r = p0ratio(M1); Me = bisect_sub(2.0*r); return r*(1 + 0.2*Me*Me)**-3.5`;
   a 60-step bisection on $M_1$ in (1, 2.197) until `exit_p(M1) = 0.7`; `assert np.allclose([M1_mine, A_ratio(M1_mine)],
   [loc["M1"], loc["A_shock_over_At"]], atol=1e-8)`.
9. `nb.note` — **N77 [B]** `nb.figure`: `GD.laval_nozzle` for $p_B/p_0$ = 0.97, 0.937, 0.8, 0.7, 0.6, 0.513, 0.3, 0.094, 0.05 with
   the shock jump drawn vertical and each curve labelled by `GD.laval_regime`. *see / read / change.* + `nb.animation` — **A3**
   (frames; 40 back pressures from 1.0 down to 0.05; stops at the seven named conditions): the current curve bold, the two
   isentropic curves as ghosts, the shock walking from the throat to the exit, the regime name and $\dot m/\dot m_{\max}$ in
   the title. *see / read / change.*
10. `nb.live` — $p_B/p_0$ (0.05 … 1) and $A_e/A_t$ (1.2 … 4) → pressure curve, shock position, regime name (paired with A3 and F6).
11. `nb.explainer("nozzle_back_pressure", heading="Why does a nozzle stop listening?", why="you open the valve on a transport:
    the mass-flow needle rises and stops dead while the pressure curve keeps reshaping — a shock is born at the throat, walks
    down the bell and leaves", tries=["Play the valve from closed to wide open and watch the needle stop.", "Stop when the shock
    is half-way down the bell: read the stagnation pressure it costs.", "Raise the supply pressure: does the plateau move?"])`
12. `nb.note` — **N76 [C]** "Outside the exit the jet meets the back pressure through oblique shocks when over-expanded and
    through expansion fans when under-expanded; C14 and C15 (§15.9) compute the first wave angle at the lip of our nozzle."
13. `nb.md` — **What would change if…** "…the duct had constant area but rough walls, or a flame inside? The same factor
    $1-M^2$ does the same work without any throat. C11."

### A.8 §15.8 Effects of Friction and Heating in Constant-Area Ducts — C11, C12
#### C11 — The outlet Mach number with friction and heat (15.46)
1. `nb.section("15.8", "Effects of Friction and Heating in Constant-Area Ducts", intro="**What is this section about?** A
   straight pipe of constant area, with wall friction and with heat added or removed. One closed formula gives the outlet
   Mach number; its consequences are the last two of the chapter's three surprises.")`
2. `nb.core("C11", "Friction and heat in one formula: $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46)", question="A gas enters a straight pipe at a known Mach
   number, rubs on the wall and is heated. What is its Mach number at the outlet?")`
3. `nb.md` — **Plain words:** "A long gas pipeline; the tail pipe of an engine; the combustion chamber of a jet. No area
   change at all — yet the flow speeds up or slows down, because the wall drags on it and heat goes in. The same three
   duct laws apply, now with a wall force and a heat input switched on. They collapse into a single number A that the
   inlet, the friction and the heat fix, and a quadratic for the outlet Mach number."
4. `nb.md` — **The idea** (ASCII): `inlet M₁ ──[ friction f · heat q ]──► outlet M₂`; "(M₁, f, q) → one number A → quadratic in
   $M_2^2$ → two roots: one continuous with $M_1$, one through a shock".
5. `nb.note` — **N79 [B]** "We start from $\rho_1u_1=\rho_2u_2,\ \ p_1+\rho_1u_1^2=p_2+\rho_2u_2^2+p_1f,\ \ h_1+\tfrac12u_1^2+h_1q=h_2+\tfrac12u_2^2$ (15.45, with $f=F/(p_1A)$, $q=Q/h_1$): C03's three laws at constant area with a retarding wall force F and a heat
   input Q, made dimensionless with inlet values. > ⚠️ Common confusion (traps T5, T11): f is **not** a friction factor and q
   is not a heat flux; F here is the *retarding* force (the opposite sign to F in $-\dot mu_1+\dot mu_2=(pA)_1-(pA)_2+F$ (15.24)). The book's lumped f is not the
   $4fL/D$ of the usual 'Fanno flow' tables — both are shown side by side in C12."
6. `nb.primer("a quadratic in M squared and its discriminant", "A quadratic $aX^2+bX+c=0$ (reminder: P159) has two roots when
   its discriminant $b^2-4ac$ is positive, one double root when it is zero and none when it is negative. Here $X=M_2^2$:
   two roots = two possible outlets; a double root = the two outlets merge (choking); none = no steady flow with this
   inlet.", code="disc = lambda A, g=1.4: 1 - 2*A*(g + 1)       # the discriminant of (15.46)\nfor A in (0.144, 0.168, 1/4.8, 0.22):
   # A grows with friction or heat\n    print(A, disc(A))                           # 0.309, 0.194, 0.0 (choked), -0.056 (no
   solution)")` (**P366**)
7. `nb.derivation("D21", …)` — Part F D21 (13 steps, ★★★, sympy check), ref "15.46": holds **N80 [B]** (steps 1–7; "the duct
   number A is a pure number, not an area — trap T2; slip #14: the book says '$M_2$ calculated from (15.45)'; it means $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46)").
8. `nb.note` — **N83 [B]** "**The check that closes D21.** With f = q = 0 the continuous root is $M_2=M_1$ and the other is
   $M_2^2=\dfrac{(\gamma-1)M_1^2+2}{2\gamma M_1^2+1-\gamma}$ (15.40) — the normal shock of C07 drops out of the duct formula. The root with $M_2>1$ from $M_1<1$ has ΔS < 0 and is
   excluded. > ⚠️ Common confusion (trap T12): which sign of ± is the continuous root flips at $M_1=1$ (+ below, − above)."
   + `nb.code`: `GD.duct_outlet_mach(2.0)`, `GD.duct_outlet_mach(2.0, branch="shocked")`, `GD.normal_shock_mach(2.0)`,
   `GD.duct_outlet_mach(0.5, branch="shocked")`, `GD.duct_entropy_change(0.5, 2.6458)`. *expect:* 2.0, 0.5774, 0.5774; 2.6458 with ΔS/$C_v$ = −0.326 < 0 (forbidden).
9. `nb.worked_example("M₁ = 0.5 with friction f = 0.1, by hand", "γ = 1.4, q = 0. 1. Numerator of $A\equiv\dfrac{M_1^2\left(1+\frac12(\gamma-1)M_1^2+q\right)}{\left(1+\gamma M_1^2-f\right)^2}=\dfrac{M_2^2\left(1+\frac12(\gamma-1)M_2^2\right)}{\left(1+\gamma M_2^2\right)^2}$ (unnumbered, §15.8): 0.25 × (1 + 0.05) =
   0.2625. 2. Denominator: (1 + 0.35 − 0.1)² = 1.25² = 1.5625. 3. A = 0.2625/1.5625 = **0.168**. 4. Discriminant: 1 − 2 × 0.168 ×
   2.4 = 0.1936; its root is 0.44. 5. $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46), upper sign: $M_2^2$ = (−(1 − 0.4704) + 0.44)/(0.4 − 0.6586) = (−0.0896)/(−0.2586) =
   0.3465 → $M_2$ = **0.589**. 6. Lower sign: (−0.9696)/(−0.2586) = 3.75 → $M_2$ = 1.936 (supersonic from subsonic: forbidden).
   7. Friction took the flow from Mach 0.50 to 0.59: it **sped up**.")`
10. `nb.code` — `A = GD.duct_parameter_A(0.5, f=0.1)`; `M2 = GD.duct_outlet_mach(0.5, f=0.1)`; `GD.duct_state_ratios(0.5, M2, f=0.1)`;
    the same with `q=0.2` (f = 0); `GD.duct_outlet_mach(2.0, f=0.1)`; `GD.duct_residuals(state1, state2, 0.1, 0.0)` (two state dicts with `rho`, `u`, `p` built from the ratios; residuals < 10⁻¹²). *expect:* A 0.1680; M2
    0.5887; `p_ratio` 0.8417, `T_ratio` 0.9819, `rho_ratio` 0.8571, `p0_ratio` 0.8971, `T0_ratio` 1.0; heat:
    A 0.1715, M2 0.6038, `T_ratio` 1.165, `p0_ratio` 0.964, `T0_ratio` 1.1905; supersonic inlet with friction: 2.0 →
    1.882 (it **slows**). *explain:* 5 points.
11. `nb.check_agree` — **from scratch (without $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46)):** unknown $u_2/u_1=r$: mass gives $\rho_2/\rho_1=1/r$; momentum
    gives $p_2/p_1=1+\gamma M_1^2(1-r)-f$; energy gives $T_2/T_1=1+q+\tfrac12(\gamma-1)M_1^2(1-r^2)$; `brentq` on the gas law
    `p2/p1 - (1/r)*T2/T1` near r = 1; `M2_mine = M1*r/np.sqrt(T2/T1)`; `assert np.allclose(M2_mine, GD.duct_outlet_mach(0.5,
    f=0.1, q=0.2))`.
12. `nb.note` — **N81 [B]**, **N82 [B]** `nb.plotly` — **F7** two `slider_figure`s (our versions of the book's two charts,
    with our own parameter values): (a) f (0 … 0.6, 25 steps; q = 0) and (b) q (−0.4 … 1.0, 29 steps; f = 0) → $M_2$ against
    $M_1$ (0 … 4): `GD.duct_curve(M1, f=…)["M2_continuous"]` (bold), `["M2_shocked"]` (thin), the diagonal, the choking locus from
    `GD.duct_choking_locus` (dots at $M_2$ = 1), the forbidden region $M_1<1<M_2$ shaded. "(Slip #15: some supersonic labels on
    the book's heating chart lack a decimal point; the turning points from `GD.duct_choking_limits` show which values were
    meant.)" *see / read / change.*
13. `nb.md` — **What would change if…** "…you asked only for the *direction* of each change, for a little friction or a
    little heat? Then the quadratic becomes a set of differentials — with the same $1-M^2$ as the area rule. C12."

#### C12 — Friction and heat drive the flow to M = 1
1. `nb.core("C12", "Everything drifts to Mach 1: $\\dfrac{du}{u}=\\dfrac{df}{1-M^2}$ and $\\dfrac{du}{u}=\\dfrac{dq}{1-M^2}$ (from
   (15.47), (15.48), corrected)", question="Can friction make a flow faster, and can heating make a gas colder?")`
2. `nb.md` — **Plain words:** "Rub your hands: friction slows things and warms them. Light a stove: heat raises the
   temperature. In a gas flowing through a pipe both everyday rules fail. Friction makes a subsonic flow *faster* (its
   pressure and density drop so quickly that the speed must rise to carry the same mass). And in a narrow band of Mach
   number just below 1, adding heat makes the gas *colder*: the heat goes into speed faster than it arrives. Both effects
   push the Mach number toward 1 from either side — and stop there: the pipe chokes."
3. `nb.md` — **The idea** (table of signs; one row per quantity, columns "friction, M < 1 · friction, M > 1 · heating, M < 1 ·
   heating, M > 1"): M ↑ ↓ ↑ ↓ · u ↑ ↓ ↑ ↓ · p ↓ ↑ ↓ ↑ · T ↓ ↑ (↑ then ↓ above $M^2=1/\gamma$) ↑ · $p_0$ ↓ ↓ ↓ ↓ · S ↑ ↑ ↑ ↑.
4. `nb.derivation("D22", …)` — Part F D22 (9 steps), ref "15.47": holds **N84 [B]** (as corrected).
5. `nb.md` — "> ⚠️ **slip #11 — the book prints** $-\dfrac{dp}{p_1}=\dfrac{1-(\gamma-1)M^2}{1-M^2}df$ (15.47, first relation)**;
   the correct form is** $-\dfrac{dp}{p_1}=\dfrac{1+(\gamma-1)M^2}{1-M^2}\,df,\ \ \dfrac{du}{u}=-\dfrac{d\rho}{\rho}=\dfrac{p_1}{p}\dfrac{df}{1-M^2}$ (15.47, numerator sign corrected). Why: mass gives $d\rho/\rho=-du/u$, energy gives $dT/T=-(\gamma-1)M^2\,du/u$, so the gas law
   gives $dp/p=-[1+(\gamma-1)M^2]\,du/u$ (D22 steps 3–6). The second relation is right as printed."
6. `nb.derivation("D23", …)` — Part F D23 (8 steps), ref "15.48": holds **N85 [B]** (as corrected).
7. `nb.md` — "> ⚠️ **slip #12 — the book prints** $\dfrac{u\,du}{h_1}=\dfrac{(1-\gamma)M^2}{1-M^2}dq$ (15.48, second
   relation)**; the correct form is** $\dfrac{dT}{T_1}=\left(\dfrac{1-\gamma M^2}{1-M^2}\right)dq,\ \ \dfrac{u\,du}{h_1}=\dfrac{(\gamma-1)M^2}{1-M^2}\,dq$ (15.48, second relation's sign corrected). Why: $du/u=dq/(1-M^2)$ and $u^2/h=(\gamma-1)M^2$ (D23 steps 6–8); the printed sign
   would make heating *slow* a subsonic flow, against the book's own text and chart. The first relation is right as printed."
8. `nb.check_agree` — **from scratch (the cell that justifies both corrections):** finite differences of $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46): `h = 1e-6`; for
   M in (0.5, 0.9, 2.0): `M2 = GD.duct_outlet_mach(M, f=h)`; `dp = GD.duct_state_ratios(M, M2, f=h)["p_ratio"] - 1`; compare
   `-dp/h` with the corrected $(1+(\gamma-1)M^2)/(1-M^2)$ and with the printed $(1-(\gamma-1)M^2)/(1-M^2)$; the same in q for
   $u\,du/h_1$. `c = GD.influence_coefficients(M, df=h)`; `assert np.allclose(-dp/h, -c["dp_over_p"]/h, rtol=1e-4)`; `bad =
   GD.influence_coefficients(M, df=h, printed_15_47=True)`; `assert not np.allclose(-dp/h, -bad["dp_over_p"]/h, rtol=0.05)`; likewise
   `printed_15_48=True`. *expect:* finite differences 1.4667, 6.968, −0.8667 = corrected; printed 1.2000, 3.558, +0.2000 (18 % and
   49 % off, and the wrong sign); heating: 0.1333, 1.7053, −0.5333 = corrected $(\gamma-1)M^2/(1-M^2)$.
9. `nb.primer("a process line on a T–s diagram (Fanno and Rayleigh lines)", "Plot temperature against entropy for all the
   states one duct flow can pass through as friction (or heat) accumulates: a curve with a 'nose' of maximum entropy at
   M = 1. Most books call the friction curve the **Fanno line** and the heating curve the **Rayleigh line** — names *not in
   our book*. Their star is the sonic state of the *same* process (trap T2).", code="M = np.array([0.5, 0.845, 1.0, 2.0])
   # Mach numbers along a heating process\nT_ratio = ((1 + 1.4)*M/(1 + 1.4*M**2))**2       # Rayleigh T/T*: (2.4 M / (1 + 1.4
   M^2))^2\nprint(T_ratio)                               # 0.79, 1.029 (the maximum, at M^2 = 1/gamma), 1.0, 0.529")` (**P367**)
10. `nb.note` — **N86 [B]** "Heating acts like friction (drives M toward 1); cooling does the opposite and always has a
    solution." + the two-column table **book / most other sources (not in the book — ours)**: lumped $f=F/(p_1A)$ ↔ friction
    length $4fL/D$ (`GD.fanno_ratios(M)["fL"]`, bridge `GD.fanno_length_to_f(M1, M2)`: Mach 0.5 → 0.589 is f = 0.100 or $4f_FL/D$ = 0.532); $q=Q/h_1$ ↔ $T_0/T_0^*$
    (`GD.rayleigh_ratios(M)["T0"]`; $T_{02}/T_{01}=1+q/(1+\tfrac12(\gamma-1)M_1^2)$). "> ⚠️ (trap T11) two ducts in
    series do not add their f, because each f is normalised with its own inlet pressure: 0.1 then 0.1 takes Mach 0.5 to
    0.589 to 0.728; a single f = 0.2 gives 0.775. The equivalent single value is $0.1+0.1\times p_2/p_1$ = 0.184."
11. `nb.note` — **N87 [B]** "From the first of $\dfrac{dT}{T_1}=\left(\dfrac{1-\gamma M^2}{1-M^2}\right)dq,\ \ \dfrac{u\,du}{h_1}=\dfrac{(\gamma-1)M^2}{1-M^2}\,dq$ (15.48, second relation's sign corrected): $dT/dq<0$ when $1/\gamma<M^2<1$ — for γ = 1.4, 0.845 < M < 1. Where does
    the heat go? Into kinetic energy, faster than it arrives. > ⚠️ Common confusion (trap T18): a statement about one band of
    Mach number, and the entropy still rises." `GD.rayleigh_temperature_sign(np.array([0.5, 0.9, 2.0]))`;
    `GD.duct_state_ratios(0.9, GD.duct_outlet_mach(0.9, q=0.005), q=0.005)["T_ratio"]`. *expect:* +1, −1, +1; 0.9947 (heat in,
    temperature down 0.5 %).
12. `nb.derivation("D24", …)` — Part F D24 (7 steps), ref "ours": holds **N88 [B]**. + `GD.duct_choking_limits(0.5)`. *expect:*
    `f_max` 0.2275, `q_max` 0.4688 for an M = 0.5 inlet (and, by a coincidence of γ = 1.4, the same `q_max` at $M_1$ = 2).
13. `nb.note` — **N89 [C]** "The funnel picture: pour faster than the neck allows and the level rises upstream — once
    choked, more friction or heat forces the *inlet* state to change. The same idea as area choking in C09 (§15.7)." · **N90
    [C]** "Dense gases near their critical point can show expansion shocks and 'antithroats' (the literature the book cites);
    for a perfect gas only compression shocks exist (N61 in C08)."
14. `nb.worked_example("a little friction at Mach 0.5 and at Mach 2", "$-\dfrac{dp}{p_1}=\dfrac{1+(\gamma-1)M^2}{1-M^2}\,df,\ \ \dfrac{du}{u}=-\dfrac{d\rho}{\rho}=\dfrac{p_1}{p}\dfrac{df}{1-M^2}$ (15.47, numerator sign corrected) with df = 0.01 at the station where $p=p_1$. 1.
    M = 0.5: du/u = 0.01/(1 − 0.25) = **+1.33 %** (faster); −dp/p = (1 + 0.4 × 0.25)/0.75 × 0.01 = **1.47 %** (pressure falls). 2.
    M = 2: du/u = 0.01/(1 − 4) = **−0.33 %** (slower); −dp/p = (1 + 1.6)/(−3) × 0.01 = −0.87 % (pressure *rises*). 3. A little heat
    dq = 0.01 at M = 0.9: dT/T = (1 − 1.4 × 0.81)/(1 − 0.81) × 0.01 = (−0.134/0.19) × 0.01 = **−0.71 %** (colder).")`
15. `nb.code` — `GD.influence_coefficients(0.5, df=0.01)`, `(2.0, df=0.01)`, `(0.9, dq=0.01)`, `(0.5, dq=0.01)`. *expect:* the worked
    numbers; `dp0_over_p0` < 0 and `ds_over_cp` > 0 in every case with df > 0 or dq > 0. *explain:* 4 points.
16. `nb.animation` — **A5** (frames; 40 steps): left, the outlet state marching along its process line on a T–s diagram
    (`GD.fanno_line_Ts`, then `GD.rayleigh_line_Ts`; title "Fanno / Rayleigh line — usual names, not in the book") as f (then q)
    grows from 0 to its choking value for $M_1$ = 0.5, stopping at the nose M = 1; right, gauges of M, $T/T_1$, $p_0/p_{01}$.
    In the heating run the temperature turns round at M = 0.845 before the nose. *see / read / change.*
17. `nb.explainer("duct_friction_heating", heading="Can friction speed a flow up?", why="you turn up f or q and the outlet point
    slides along the M₂(M₁) map toward the line M₂ = 1; the badge flips to 'choked' exactly when the discriminant bar reaches
    zero; flipping the inlet to supersonic reverses every arrow", tries=["Raise f from 0 at $M_1$ = 0.5: which way does $M_2$
    go?", "Set $M_1$ = 0.9 and add a little heat: watch the temperature bar.", "Push f up to the choking limit, then beyond."])`
18. `nb.md` — **What would change if…** "…we left the duct for open air? 'No news upstream' becomes 'no news outside a
    cone'. C13."

### A.9 §15.9 Pressure Waves in Planar Compressible Flow — C13, C14, C15
#### C13 — The Mach cone and angle (15.49)
1. `nb.section("15.9", "Pressure Waves in Planar Compressible Flow", intro="**What is this section about?** The same ideas in
   two dimensions: where the news of a supersonic body can reach (a cone), how a supersonic stream is turned into itself
   (an oblique shock) and away from itself (an expansion fan).")`
2. `nb.core("C13", "The Mach cone: $\sin\mu=1/M$ (15.49)", question="Why do you hear a supersonic aircraft only after it has passed
   overhead?")`
3. `nb.md` — **Plain words:** "A duck paddling slowly sends ripples ahead of itself. A speedboat outruns its own ripples
   and drags a V behind it. An aircraft does the same with sound. Each instant it emits a sphere of 'news' that grows at
   the speed of sound from where the aircraft *was*. If the aircraft is faster than sound the spheres can never get ahead
   of it; they pile up along a cone. Outside the cone: silence. On it: the boom."
4. `nb.md` — **The idea** (ASCII triangle): source at the apex, a circle of radius cΔt centred UΔt behind, the tangent from
   the apex; "sin μ = (cΔt)/(UΔt) = 1/M".
5. `nb.note` — **N91 [B]** `nb.figure`: `draw_wavefront_panels([AC.source_wavefronts(U, c0, t, 8, dt) for U in (0, 0.5*c0,
   2*c0)])` — source at rest, subsonic, supersonic; each circle labelled with its emission time. *see / read / change.* +
   `nb.animation` — **A4** (video; 60 frames): fronts from a source at M = 0.5 and at M = 2 side by side, the envelope
   appearing in the second. *see / read / change.*
6. `nb.primer("the envelope of a family of circles", "A curve that touches every member of a family is its envelope. For
   circles whose centres lie on a line and whose radii grow in proportion to the distance from a point P, the envelope is
   the pair of tangents from P. A tangent meets the radius at its touching point at a right angle — that right angle makes
   the triangle.", code="d, r = 2.0, 1.0                        # centre 2 behind the apex, radius 1\nmu =
   np.degrees(np.arcsin(r/d))           # angle of the tangent from the apex\nprint(mu)                              # 30
   degrees: sin(mu) = r/d")` (**P368**)
7. `nb.derivation("D25", …)` — Part F D25 (4 steps), ref "15.49". "(μ is the Mach angle here; a viscosity in §15.6 — trap T5.)"
8. `nb.note` — **N92 [B]** "Inside the cone: the **zone of action**; outside: the **zone of silence**. The cone's
   generators are **Mach lines**. (Ours:) an aircraft at height h passes overhead; the boom arrives $t_{\text{delay}}=\dfrac{h\sqrt{M^2-1}}{Mc}$ (ours) later." +
   `AC.mach_cone_delay(1.0e4, 1.5, 295.1)`, `AC.in_zone_of_action(x, y, M)` (x measured **downstream** of the source — behind the aircraft). *expect:* 25.3 s for 10 km, M = 1.5, c = 295 m/s. +
   `nb.primer("hyperbolic against elliptic equations", "Elliptic: every point influences every other (subsonic flow, Laplace's
   equation; reminder P201). Hyperbolic: influence travels only along special lines, the characteristics (waves; reminder
   P174, P227). Steady supersonic flow is hyperbolic with the downstream distance playing the part of time, and its
   characteristics are the Mach lines.", code="for M in (0.5, 2.0):                    # the sign of 1 - M^2 decides the
   type\n    print(M, 'elliptic' if 1 - M*M > 0 else 'hyperbolic')")` (**P369**)
9. `nb.worked_example("three Mach numbers", "$\sin\mu=1/M$ (15.49): M = 2 → sin μ = 0.5 → μ = **30°**. M = 1.2 → sin μ = 0.833 → μ = 56.4° (a
   blunt cone). M = 3 → μ = 19.5° (a slender one). M = 1 → μ = 90°: a plane front. M < 1: sin μ > 1, no cone.")`
10. `nb.code` — `np.degrees(GD.mach_angle(np.array([1.0, 1.2, 1.5, 2.0, 3.0])))`; `GD.mach_angle(0.8)`; `GD.mach_from_angle(np.deg2rad(
    30.0))`. *expect:* 90, 56.44, 41.81, 30.00, 19.47; nan; 2.0. *explain:* 3 points.
11. `nb.check_agree` — **from scratch:** 200 circles from `AC.source_wavefronts(2*c0, c0, 1.0, 200, 0.005)` (keys `x_centres`, `radii`); the envelope's upper
    edge as the maximum over circles of y at each x; a straight-line fit of that edge; `assert np.allclose(np.arctan(-slope),
    GD.mach_angle(2.0), rtol=1e-3)`.
12. `nb.explainer("mach_cone", heading="Why do you hear it only after it has passed?", why="you play time and each circle grows
    from where the source was; above Mach 1 the circles cross and their envelope appears by itself; you drag an observer and
    read 'not yet … now … heard'", tries=["Play at M = 0.5, then M = 1, then M = 2.", "Click a circle: where was the source when
    it left?", "Drag the observer below the flight path and read the delay."])`
13. `nb.md` — **What would change if…** "…the disturbance were not infinitely weak — a wedge that really turns the flow? The
    Mach line steepens into an oblique shock. C14."

#### C14 — The oblique shock (15.53)
1. `nb.core("C14", "Turning a supersonic stream into itself: $\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}$ (15.53)", question="A supersonic stream meets a wedge. At what
   angle does the shock stand — and how sharp can the wedge be before no such shock exists?")`
2. `nb.md` — **Plain words:** "The nose of a supersonic aircraft, the lip of an engine intake, a bullet: each forces the
   air to turn. The air cannot be warned (no news upstream), so it turns all at once, through a thin sloping front. Here
   is the trick that makes it easy: walk along the front at just the right speed and the flow you see comes straight at
   the front — a *normal* shock. So an oblique shock is a normal shock on the component of the velocity normal to it; the
   component along it passes through untouched."
3. `nb.md` — **The idea** (ASCII triangles): `V₁ = (u₁ ⟂ shock, v ∥ shock) ──► V₂ = (u₂ < u₁, v unchanged)`; "normal part shocked:
   replace $M_1$ by $M_1\sin\sigma$ everywhere in C07".
4. `nb.note` — **N93 [B]** the geometry: shock angle σ from the upstream direction, deflection δ;
   $(u_1,v)=V_1(\sin\sigma,\cos\sigma)$, $(u_2,v)=V_2(\sin(\sigma-\delta),\cos(\sigma-\delta))$ + `nb.figure`: `draw_shock_triangles(
   2.0, sigma, delta)`. "> ⚠️ Common confusion (trap T13): the book's σ and δ are the β and θ of most other sources; both
   pairs are in the legend." *see / read / change.*
5. `nb.primer("resolving a velocity along and across a tilted line", "Any velocity is the sum of a part along a line and
   a part across it: for a line at angle σ to the velocity, across = V sin σ, along = V cos σ. An observer sliding along the
   line at the 'along' speed sees only the 'across' part — and sliding at constant speed changes no pressure, density or
   temperature (reminder: P96).", code="V, sigma = 680.6, np.deg2rad(39.31)            # Mach 2 at sea level; a 39.31 degree
   shock\nprint(V*np.sin(sigma), V*np.cos(sigma))         # across 431.2 m/s (Mach 1.267), along 526.6 m/s")` (**P370**)
6. `nb.derivation("D26", …)` — Part F D26 (8 steps), ref "15.50": holds **N94 [B]** (step 4), **N95 [B]**, **N96 [B]**, **N97
   [B]** (steps 6–8; "⚠️ trap T13: $M_2$ needs the division by $\sin(\sigma-\delta)$"). `ch15.tangential_momentum_sympy()`.
7. `nb.primer("the tangent subtraction formula", "$\\tan(\\sigma-\\delta)=\\dfrac{\\tan\\sigma-\\tan\\delta}{1+\\tan\\sigma\\tan
   \\delta}$, and the double angle $\\cos2\\sigma=1-2\\sin^2\\sigma$ (reminder: P341, ch14).", code="s, d = np.deg2rad(39.31),
   np.deg2rad(10.0)\nprint(np.tan(s - d), (np.tan(s) - np.tan(d))/(1 + np.tan(s)*np.tan(d)))   # 0.5613 both\nprint(np.cos(2*s), 1 -
   2*np.sin(s)**2)                                    # 0.1972 both")` (**P371**)
8. `nb.derivation("D27", …)` — Part F D27 (8 steps), ref "15.53".
9. `nb.worked_example("a 10° wedge at Mach 2", "γ = 1.4. Trial: σ = 39.31° (sin σ = 0.6336, sin² σ = 0.4014, cos 2σ = 0.1972, cot σ
   = 1.2213). 1. $\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}$ (15.53): tan δ = 2 × 1.2213 × (4 × 0.4014 − 1)/(4 × (1.4 + 0.1972) + 2) = 2.4426 × 0.6056/8.389 = 0.1763 → δ =
   **10.0°** ✓. 2. Normal Mach number: $M_{n1}$ = 2 × 0.6336 = **1.267**. 3. $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2\sin^2\sigma-1\right)$ (15.50): $p_2/p_1$ = 1 + 1.1667 × (1.6056 − 1) = **1.707**.
   4. $M_2^2=\dfrac{(\gamma-1)M_1^2+2}{2\gamma M_1^2+1-\gamma}$ (15.40) on 1.267: $M_{n2}$ = 0.803; divide by sin(39.31° − 10°) = 0.4896: $M_2$ = **1.64** — still supersonic. 5. The other
   solution of $\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}$ (15.53) for the same 10°: σ = 83.7°, $p_2/p_1$ = 4.44, $M_2$ = 0.60 — nearly a normal shock.")`
10. `nb.code` — `w = GD.oblique_shock(2.0, delta=wedge["delta"])`; `s = GD.oblique_shock(2.0, delta=wedge["delta"], branch="strong")`;
    `GD.oblique_shock_components(2.0, np.deg2rad(w_sigma_deg))`; `GD.max_deflection(2.0)`; `GD.oblique_shock(2.0, delta=np.deg2rad(25.0))` in a `try` (it raises); `GD.detachment_mach(np.deg2rad(10.0))`. *expect:* weak σ 39.31°, `Mn1` 1.2671, `p_ratio` 1.7066, `rho_ratio` 1.4584, `T_ratio` 1.1702, `M2` 1.6405, `p0_ratio` 0.9846; strong σ 83.70°, 4.4438, `M2` 0.6037;
    `delta_max` = 22.97° at `sigma_max` = 64.67°; `ValueError("detached: …")`; the smallest Mach number with an attached shock on a 10° wedge is 1.421. *explain:* 5 points.
11. `nb.check_agree` — **from scratch:** `tan_delta = lambda s: 2/np.tan(s)*(M1**2*np.sin(s)**2 - 1)/(M1**2*(g + np.cos(2*s)) + 2)`; a
    grid search for the σ of the largest δ; bisection on (μ, σ_max) and on (σ_max, π/2) for δ = 10°; `assert
    np.allclose([s_weak, s_strong], [w["sigma"], s["sigma"]], atol=1e-8)`.
12. `nb.note` — **N98 [B]** "(the closed form of $\delta_{\max}$ is ours; result given, not derived) For δ below a maximum
    there are two shocks: **weak** (smaller σ, usually supersonic behind) and **strong** (larger σ, subsonic behind); δ = 0 at
    σ = μ (a Mach wave) and at σ = 90° (a normal shock). > ⚠️ (trap T13) 'weak' and 'strong' name the two roots, not their
    size: the 'weak' shock for 20° at Mach 2 nearly triples the pressure ($p_2/p_1$ = 2.84)." + `nb.plotly` — **F8**
    `slider_figure`: $M_1$ (1.2 … 5, 20 steps) → δ against σ from `GD.oblique_shock_curve` with the maximum-deflection point and
    the sonic point ($M_2$ = 1) marked, axis labels "σ (= β elsewhere)", "δ (= θ elsewhere)". *expect:* $\delta_{\max}$ 22.97° at
    $M_1$ = 2; 34.07° at 3; 45.58° as $M_1\to\infty$. *see / read / change.*
13. `nb.note` — **N99 [B]** `nb.figure`: `draw_wedge_streamlines(2.0, np.deg2rad(10.0), ch15.wedge_flow(2.0, np.deg2rad(10.0))["sigma"])` (streamlines turning through the attached weak shock; `wedge_flow` returns the numbers — `sigma`, `M2`, `p_ratio`, `attached`, `subsonic_behind`, `delta_max` — the helper draws) and the same for 22° (nearly detached); a third panel `draw_bow_shock_sketch()` titled
    "detached bow shock — a sketch, not a computation" (normal on the axis, strong near it, weak far out). "A wedge takes the
    weak solution; a blunt nose always has a detached shock." *see / read / change.* + `nb.live` — $M_1$ and δ → the wedge picture
    and both roots (paired with F8).
14. `nb.note` — **N100 [C]** "At a corner the speed is finite in supersonic flow (zero or infinite in subsonic flow — Ch. 6
    corner flows, C09), and a supersonic field can be built marching downstream: the hyperbolic property of N92."
15. `nb.primer("the Mach triangle", "A right triangle with hypotenuse M and opposite side 1 has the Mach angle μ: $\\sin\\mu
    =1/M$, so $\\cos\\mu=\\sqrt{M^2-1}/M$ and $\\cot\\mu=\\sqrt{M^2-1}$. Every $\\sqrt{M^2-1}$ of §15.9–§15.10 is this
    triangle's third side.", code="M = 2.0\nmu = np.arcsin(1/M)                             # Mach angle: 30
    degrees\nprint(1/np.tan(mu), np.sqrt(M*M - 1))           # cot(mu) = sqrt(M^2 - 1) = 1.732")` (**P372**)
16. `nb.derivation("D28", …)` — Part F D28 (8 steps), ref "15.55": holds **N101 [B]** (steps 1–3) and **N102 [B]** (steps 4–8).
17. `nb.note` — **N102 [B]** "> ⚠️ Common confusion (trap T16): how small is 'small'?" — table **linear rule $\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta$ (15.55) against the
    exact weak shock at $M_1$ = 2** (`GD.weak_wave_pressure`, `GD.oblique_shock`, `GD.weak_wave_error`): δ = 0.5°: 0.0282 / 0.0285
    (rule 1.1 % low) · 2°: 0.1129 / 0.1180 (4.3 % low) · 5°: 0.2821 / 0.3154 (10.5 % low; the exact value is 12 % above the
    rule) · 10°: 0.5643 / 0.7066 (20 % low) · 20°: 1.129 / 1.843 (39 % low). "And the rule fails as $M\to1$: $\sqrt{M^2-1}\to0$."
18. `nb.explainer("oblique_shock_wedge", heading="One wedge, how many shocks?", why="you open the wedge angle and the two
    solution dots on the δ–σ curve run toward each other, merge at the maximum and vanish — and the picture switches to a
    detached bow shock", tries=["Open the wedge from 2° to 23° at Mach 2 and watch the two dots.", "Click a streamline: which
    velocity component did the shock change?", "Switch to the strong root: how far is it from a normal shock?"])`
19. `nb.md` — **What would change if…** "…the wall turned *away* from the stream? There is no expansion shock (C08). The
    stream fans out instead. C15."

#### C15 — The Prandtl–Meyer function (15.56)
1. `nb.core("C15", "Turning away, one Mach wave at a time: $\nu(M)=\sqrt{\dfrac{\gamma+1}{\gamma-1}}\tan^{-1}\sqrt{\dfrac{\gamma-1}{\gamma+1}\left(M^2-1\right)}-\tan^{-1}\sqrt{M^2-1}$ (15.56)", question="How does a supersonic stream get round a
   corner that bends away from it — and how fast is it afterwards?")`
2. `nb.md` — **Plain words:** "Behind the shoulder of a bullet, at the lip of a rocket nozzle in thin air, over the top of
   a tilted supersonic wing, the wall falls away from the stream. A slow flow would simply follow, or separate. A
   supersonic one turns through a *fan* of infinitely many weak waves spreading from the corner. Each wave turns the
   stream a little and speeds it up a little; nothing is lost. Add the little turns up and you get one function, ν(M), and
   every such turn becomes a subtraction."
3. `nb.md` — **The idea** (ASCII fan): corner with rays; "first ray at $\mu_1$ to the upstream wall · last ray at $\mu_2$ to
   the downstream wall · turn = ν(M₂) − ν(M₁)".
4. `nb.note` — **N105 [B]** the fan: all Mach lines start at the corner; isentropic; along a streamline the pressure falls
   smoothly, along the wall it steps at the corner + `nb.figure`: `draw_fan(ch15.expansion_fan(2.0, np.deg2rad(10.0)))` (rays
   and one streamline coloured by M). *see / read / change.*
5. `nb.derivation("D29", …)` — Part F D29 (9 steps), ref "unnumbered, §15.9": holds **N106 [B]**. `ch15.prandtl_meyer_sympy()`.
6. `nb.primer("arctangent integrals", "$\\displaystyle\\int\\dfrac{dt}{1+a^2t^2}=\\dfrac1a\\tan^{-1}(at)$: the antiderivative
   of 'one over one plus a square' is an arctangent, with the scaling a divided out. Reached here after a substitution
   (reminder P106) and a split into partial fractions (reminder P219).", code="import sympy as sp\nt, a = sp.symbols('t a',
   positive=True)\nprint(sp.integrate(1/(1 + a**2*t**2), t))            # atan(a*t)/a\nprint(sp.diff(sp.atan(a*t)/a, t))
   # 1/(a**2*t**2 + 1): differentiate back")` (**P373**)
7. `nb.derivation("D30", …)` — Part F D30 (10 steps, ★★★, sympy check), ref "15.56". "(ν is the Prandtl–Meyer function here; a
   kinematic viscosity in §15.6 — trap T5.)"
8. `nb.note` — **N107 [B]** "**How it is used:** $\delta_1+\nu(M_1)=\delta_2+\nu(M_2)$ (unnumbered, §15.9), so $\nu(M_2)=\delta_1-\delta_2+\nu(M_1)$. > ⚠️ Common confusion (trap
   T14): the book counts δ positive when the wall turns **into** the flow; our `turn=` (and most tables) count the
   **expansion** positive: turn = −δ. Both numbers are in every legend."
9. `nb.worked_example("Mach 2 turned 10° away", "γ = 1.4; $\\sqrt{(\\gamma+1)/(\\gamma-1)}$ = √6 = 2.449. 1. $\nu(M)=\sqrt{\dfrac{\gamma+1}{\gamma-1}}\tan^{-1}\sqrt{\dfrac{\gamma-1}{\gamma+1}\left(M^2-1\right)}-\tan^{-1}\sqrt{M^2-1}$ (15.56) at M = 2:
   $\\sqrt{M^2-1}$ = 1.732; first term 2.449 × atan(1.732/2.449) = 2.449 × atan(0.7071) = 2.449 × 35.26° = 86.38°; second atan(1.732)
   = 60°; **ν(2) = 26.38°**. 2. Turn +10° (the book's δ = −10°): ν₂ = 36.38°. 3. Invert (root-finder): **M₂ = 2.385**. 4.
   Isentropic, so $p_0$ is unchanged: $p_2/p_1$ = [(1 + 0.2 × 4)/(1 + 0.2 × 2.385²)]^3.5 = (1.8/2.1375)^3.5 = **0.548**. 5. Fan: first
   ray at $\\mu_1$ = 30° to the old wall, last at $\\mu_2$ = 24.8° to the new wall.")`
10. `nb.code` — `np.degrees(GD.prandtl_meyer(np.array([1.0, 1.5, 2.0, 3.0])))`; `np.degrees(GD.prandtl_meyer_max())`; `t =
    GD.prandtl_meyer_turn(2.0, np.deg2rad(10.0))`; `GD.mach_from_prandtl_meyer(np.deg2rad(36.38))`. *expect:* 0, 11.905, 26.380,
    49.757; 130.45°; `M2` 2.3849, `p_ratio` 0.5480, `T_ratio` 0.8421, `mu2` 24.79°, `delta_book` −10°. *explain:* 4 points.
11. `nb.check_agree` — **from scratch:** `nu_quad = quad(lambda M: np.sqrt(M*M - 1)/(1 + 0.2*M*M)/M, 1.0, 2.0)[0]` (reminder P87);
    `assert np.allclose(nu_quad, GD.prandtl_meyer(2.0))`; then 100 small steps of $d\delta=-\dfrac{\sqrt{M^2-1}}{M}\dfrac{dM}{1+\frac12(\gamma-1)M^2}$ (unnumbered, §15.9): `M += dturn*M*(1 + 0.2*M*M)/np.sqrt(M*M -
    1)`; `assert np.allclose(M, t["M2"], rtol=2e-3)`; `ch15.fan_staircase(2.0, np.deg2rad(10.0), 100)["M"][-1]`. *expect:* 26.380°;
    staircase N = 1: 2.363, 3: 2.377, 7: 2.381, 21: 2.384, 100: 2.3846 → exact 2.3849.
12. `nb.note` — **N104 [B]** `nb.figure`: `ch15.curved_wall_mach_lines(2.0, wall_angle_fn, x)` for a wall curving away (Mach
    lines fan out) and into the stream (they run together; the first crossing is where a shock must start). "Along each Mach
    line the flow direction and Mach number are constant (a simple wave — the steady twin of N22's steepening)." *see / read
    / change.*
13. `nb.note` — **N103 [B]** `nb.figure`: $(p_2-p_1)/p_1$ against δ from −15° to +15° at $M_1$ = 2: the exact weak shock (δ >
    0), the exact fan (δ < 0) and the straight line $\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta$ (15.55). "They share value and slope at δ = 0 because the entropy change
    is third order: $\dfrac{S_2-S_1}{C_v}\cong\dfrac{\gamma^2-1}{12\gamma^2}\left(\dfrac{p_2-p_1}{p_1}\right)^3$ (15.44b) — so the rule holds for weak expansions too, with the sign of δ reversed." *expect:* at +10°: shock
    0.7066, isentropic compression 0.7052, line 0.5643; at −10°: fan −0.4520, line −0.5643. *see / read / change.*
14. `nb.plotly` — **F9b** `slider_figure`: γ (1.2, 1.3, 1.4, 5/3) → ν(M) in degrees against M (1 … 10, log axis) with its
    limit `GD.prandtl_meyer_max` dashed. *see / read / change.*
15. `nb.explainer("prandtl_meyer_fan", heading="How does a supersonic stream turn a corner?", why="you build the turn from 1, 3,
    7, 21, then infinitely many Mach waves and watch the staircase of Mach number settle onto ν(M): the integral as a limit
    of sums", tries=["Set N = 1, 3, 7, 21, ∞ and compare $M_2$.", "Reverse the turn: where do the Mach lines cross?", "Start at
    $M_1$ = 1: how far can the stream turn?"])`
16. `nb.md` — **What would change if…** "…a thin plate were tilted in the stream — a fan on one side, a weak shock on the
    other? Two pressures, and their difference is lift. C16."

### A.10 §15.10 Thin Airfoil Theory in Supersonic Flow — C16; exercises, literature, summary
#### C16 — The supersonic flat plate (15.59)
1. `nb.section("15.10", "Thin Airfoil Theory in Supersonic Flow", intro="**What is this section about?** Lift and drag of a thin
   plate in a supersonic stream from the weak-wave rule alone: two pressure levels, a lift — and a drag that needs no
   viscosity.")`
2. `nb.core("C16", "Lift and wave drag of a flat plate: $C_L\cong\dfrac{4\alpha}{\sqrt{M_\infty^2-1}},\ \ C_D\cong\dfrac{4\alpha^2}{\sqrt{M_\infty^2-1}}$ (15.59)", question="How does a thin wing make lift above Mach 1 —
   and why does that lift always cost drag, even without friction?")`
3. `nb.md` — **Plain words:** "Chapter 14's wing made lift by circulation, needed a sharp trailing edge to choose it, and in
   two dimensions had no drag at all. A supersonic wing is simpler and harsher. Each surface feels only its own slope:
   the top, turned away from the stream, sits at one lower pressure; the bottom, turned into it, at one higher pressure.
   The difference pushes at right angles to the plate. Tilt the plate by α and that push leans back by α: every newton of
   lift costs α newtons of drag — carried away by the waves, not by friction."
4. `nb.md` — **The idea** (ASCII plate): `p₂ < p∞ (fan above the leading edge)` / `═══ plate at α ═══` / `p₃ > p∞ (shock below)`;
   "force = (p₃ − p₂) b, normal to the plate → L = force × cos α, D = force × sin α".
5. `nb.note` — **N108 [B]** "$\dfrac{p-p_\infty}{p_\infty}=\dfrac{\gamma M_\infty^2\,\delta}{\sqrt{M_\infty^2-1}}$ (15.57) is $\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta$ (15.55) written with free-stream values (Ackeret's linear theory). As a pressure
   coefficient (ours): $C_p=2\delta/\sqrt{M_\infty^2-1}$ (🔁 $C_p$: Ch. 4 C15; dynamic pressure in a gas: P359 in C04)."
   `ch15.ackeret_pressure(np.deg2rad(2.0), 2.0)`, `ch15.ackeret_cp(np.deg2rad(2.0), 2.0)` (the angle comes first). *expect:* 0.1129; 0.0403.
6. `nb.derivation("D31", …)` — Part F D31 (7 steps), ref "15.59": holds **N109 [B]** (steps 1–3; slip #13: "> ⚠️ **slip #13 —
   the book prints** 'found from (15.47)'**; the correct form is** from $\dfrac{p-p_\infty}{p_\infty}=\dfrac{\gamma M_\infty^2\,\delta}{\sqrt{M_\infty^2-1}}$ (15.57)") and **N110 [B]** (steps 4–5; "⚠️ trap T15: b is
   the chord here") + `nb.figure`: `draw_plate_waves(ch15.flat_plate_wave_pattern(plate["alpha"], 2.0))` (the plate, its four
   waves, the two pressure levels, the resultant split into lift (teal) and wave drag (rose)). *see / read / change.*
7. `nb.worked_example("2° at Mach 2", "α = 2° = 0.0349 rad, $M_\\infty$ = 2, $\\sqrt{M_\\infty^2-1}$ = 1.732. 1. $\dfrac{p-p_\infty}{p_\infty}=\dfrac{\gamma M_\infty^2\,\delta}{\sqrt{M_\infty^2-1}}$ (15.57): each surface
   differs from $p_\\infty$ by 1.4 × 4 × 0.0349/1.732 = **11.3 %** (up below, down above). 2. $C_L\cong\dfrac{4\alpha}{\sqrt{M_\infty^2-1}},\ \ C_D\cong\dfrac{4\alpha^2}{\sqrt{M_\infty^2-1}}$ (15.59): $C_L$ = 4 × 0.0349/1.732 =
   **0.0806**. 3. $C_D$ = 4 × 0.0349²/1.732 = **0.00281**. 4. $C_D/C_L$ = α = 0.0349: L/D = 28.6. 5. Slope: 4/1.732 = 2.31 per
   radian — against Ch. 14's 2π = 6.28 for the same plate in slow flow.")`
8. `nb.code` — `lin = ch15.flat_plate_supersonic(plate["alpha"], 2.0, method="linear")`; `ex = ch15.flat_plate_supersonic(plate["alpha"],
   2.0, method="shock_expansion")`; `ch15.ackeret_flat_plate(plate["alpha"], 2.0)`; the pair again at α = 10°. *expect:* 2°: linear
   `CL` 0.08061, `CD` 0.002814, `p_upper` 0.8871, `p_lower` 1.1129; exact `CL` 0.08064, `CD` 0.002816, `p_upper` 0.8920, `p_lower`
   1.1180. 10°: linear 0.4031 / 0.0703; exact (`CN` 0.4138) 0.4075 / 0.0719 — the lift is only 1.1 % off although each surface pressure is
   20–25 % off (the second-order errors of the two surfaces nearly cancel in the difference). *explain:* 5 points.
9. `nb.check_agree` — **from scratch:** `dp = g*M**2*alpha/np.sqrt(M**2 - 1)`; `p3, p2 = 1 + dp, 1 - dp`; `qinf = 0.5*g*M**2`;
   `CN_mine = (p3 - p2)/qinf`; `CL_mine, CD_mine = CN_mine, CN_mine*alpha` (cos α ≈ 1, sin α ≈ α — the two small-angle steps of
   D31); `assert np.allclose([CL_mine, CD_mine], [lin["CL"], lin["CD"]])`; `ack = ch15.ackeret_flat_plate(alpha, M)`; `assert
   np.allclose([CL_mine, CD_mine], [ack["CL"], ack["CD"]])`; and the exact route by hand: `pl = GD.oblique_shock(M,
   delta=alpha)["p_ratio"]`, `pu = GD.prandtl_meyer_turn(M, alpha)["p_ratio"]`, `CN = (pl - pu)/qinf`, `assert np.allclose([CN*np.cos(alpha),
   CN*np.sin(alpha)], [ex["CL"], ex["CD"]])`.
10. `nb.note` — **N111 [B]** "Wave drag exists in inviscid supersonic flow: d'Alembert's paradox (🔁 Ch. 6 C07) does not
    apply, because waves carry energy away — like a ship's bow waves. No Kutta condition is needed (🔁 Ch. 14 C02), yet the
    flow has circulation. The theory fails as $M_\infty\to1$." + `nb.figure`: `ch15.lift_slope_vs_mach(M)` for M from 0 to 4 — 2π
    of $C_L\simeq2\pi(\alpha+\beta)$ (14.12) on the subsonic side, $4/\sqrt{M_\infty^2-1}$ of $C_L\cong\dfrac{4\alpha}{\sqrt{M_\infty^2-1}},\ \ C_D\cong\dfrac{4\alpha^2}{\sqrt{M_\infty^2-1}}$ (15.59) on the supersonic side, a gap round 1 (a Prandtl–Glauert
    curve only as a dashed line labelled "not in the book"). *expect:* slope 2.309 at $M_\infty$ = 2; the two slopes would be
    equal at $M_\infty$ = 1.185. *see / read / change.*
11. `nb.note` — **N112 [B]** "(ours; result given, not derived) a thin section of any shape: $C_p=\pm\dfrac{2\theta}{\sqrt{M_\infty^2-1}},\ \ C_D=\dfrac{4}{\sqrt{M_\infty^2-1}}\left[\alpha^2+\overline{\theta_c^2}+\overline{\theta_t^2}\right]$ (ours), where θ is the local
    surface slope and the bars are chord averages of the squared camber-line and thickness slopes. **Lift depends on the
    incidence only; thickness and camber buy nothing but drag.**" `ch15.ackeret_airfoil(yu, yl, 0.0, 2.0)` for a 5 % diamond.
    *expect:* `CL` 0, `CD` = `CD_thickness` = 4 × 0.05²/1.732 = 0.00577 — twice the drag of the flat plate at 2°, with no lift. · **N113 [C]**
    "Supersonic sections are thin with sharp leading edges; a rounded nose brings the detached shock of N99 (C14) and its
    drag."
12. `nb.plotly` — **F9a** `slider_figure`: $M_\infty$ (1.2 … 4, 15 steps) → $C_L$ and $C_D$ (× 10) against α (0 … 12°), linear
    rule (dashed) against exact shock–expansion (solid), the exact curves stopping where the shock detaches. *see / read /
    change.*
13. `nb.explainer("supersonic_flat_plate", heading="How can there be drag without friction?", why="you tilt the plate: the fan
    and the shock swap sides at α = 0 and the force arrow, always normal to the plate, leans back by exactly α — the drag
    bar grows as α², the lift bar as α", tries=["Tilt from 0 to 10°: compare how the two bars grow.", "Slide $M_\\infty$ toward
    1: which bars run away?", "Switch to a diamond at zero incidence: drag without lift."])`
14. `nb.md` — **What would change if…** "…the 'gas' were blood and the 'duct' an elastic artery? The pulse is a pressure wave of exactly the kind C02 solved, $p'(x,t)=f(x-ct)+g(x+ct)$ (15.15), with a wave speed set by the wall — and a collapsible tube chokes like the nozzle of C09. Chapter 16."
15. `nb.pointer` — **S01** "The chapter's exercises are not reproduced. The results they hold are taught above: the
    initial-value form (N21), the particle velocity (D04), the small-Mach expansions (N32), the area ratio (D10), the weak
    shock (D17), the maximum speed (N38), the strong shock (R17), the two sets of differentials (D22, D23), the choked mass
    flow (D11), the choking limits (D24), the thin section (N112) and the oblique-shock components (D26)." **S02** "Sources:
    the reference list at the end of the book's chapter; our public sources (NACA Report 1135 and the pages cited in the
    analysis) are in `reference/ch15/SOURCES.md`."
16. `nb.summary(clicked=[…16…], feeds_forward=[…], left_out=[…])` — **clicked** (one sentence per CORE): C01 "Sound is what the
    laws of mass and momentum become when the disturbance is weak: one wave equation with the speeds U ± c." · C02 "Any plane
    sound is one shape going right plus one going left; the air moves with a compression and against a rarefaction." · C03
    "A duct flow is three budgets between two stations, and the energy budget holds even across friction and shocks." · C04
    "$T_0$ is what an adiabatic flow keeps; $p_0$ is what it loses when entropy is made." · C05 "The area a flow needs is least
    at Mach 1, so each area belongs to two Mach numbers and a throat passes only so much." · C06 "Above Mach 1 the gas thins
    faster than it speeds up, so a supersonic nozzle widens." · C07 "Across a shock the three laws leave one solution besides
    'nothing', fixed by $M_1$ alone." · C08 "Only the supersonic-to-subsonic jump makes entropy; its price is stagnation
    pressure." · C09 "When the exit is sonic the nozzle cannot hear the back pressure and the mass flow stops rising." · C10
    "With a sonic throat the back pressure only chooses where a shock stands." · C11 "Friction and heat in a straight duct
    reduce to one number A and a quadratic for $M_2^2$." · C12 "The same $1-M^2$ makes friction and heat push every duct flow
    toward Mach 1 and choke it." · C13 "The Mach cone is the edge of all the circles of old news; $\sin\mu=1/M$." · C14 "An
    oblique shock is a normal shock on the normal component: two shocks per turn, a largest turn." · C15 "A supersonic stream
    turns away through a fan, and ν(M) keeps the count." · C16 "On a thin plate the push is normal to the plate, so lift
    costs α times itself in wave drag." **feeds_forward:** Ch. 16 — pulse waves in arteries are $p'(x,t)=f(x-ct)+g(x+ct)$ (15.15) with a wall-compliance
    wave speed; flow limitation in collapsible tubes is choking (C09) with a tube law in place of the gas law. **left_out:**
    the 1-D Euler solver through breaking (N114); the shock-structure derivations (stated: N66–N69); the Riemann-invariant
    derivation (N22); non-perfect gases (N90); method of characteristics in two dimensions.

### A.11 Placement check (every curation id has exactly one home)
- **CORE (16):** C01, C02 (§15.2) · C03 (§15.3) · C04, C05 (§15.4) · C06 (§15.5) · C07, C08 (§15.6) · C09, C10 (§15.7) · C11, C12
  (§15.8) · C13, C14, C15 (§15.9) · C16 (§15.10).
- **RECAP (19):** R01–R06 (§15.1, before the C01 block opens) · R07, R08 (C01) · R09–R13 (C03) · R14 (C04) · R15 (C06) · R17 (C07) ·
  R16, R18, R19 (C08).
- **NOTE (114):** §15.1 N01–N03 · C01 N04–N14 · C02 N15–N22, N114 · C03 N23–N28 · C04 N29–N34, N38 · C05 N35–N37, N39 · C06
  N40–N47 · C07 N48–N57, N64 · C08 N58–N63, N65–N71 · C09 N72–N74 · C10 N75–N78 · C11 N79–N83 · C12 N84–N90 · C13 N91, N92 · C14
  N93–N102 · C15 N103–N107 · C16 N108–N113.
- **SKIP (2):** S01, S02 in one `nb.pointer` row (A.10 row 15).
- **Derivations (31):** D01, D02 (C01) · D03, D04, D05 (C02) · D06, D07 (C03) · D08, D09 (C04) · D10, D11 (C05) · D12 (C06) · D13,
  D14, D15 (C07) · D16, D17, D18 (C08) · D19 (C09) · D20 (C10) · D21 (C11) · D22, D23, D24 (C12) · D25 (C13) · D26, D27, D28 (C14) ·
  D29, D30 (C15) · D31 (C16).
- **Primers (20):** P355 (§15.1, R05) · P354, P356 (C01) · P357, P358 (C02) · P359 (C04) · P360, P361 (C05) · P362 (C06) · P363
  (C07) · P364 (C08) · P365 (C10) · P366 (C11) · P367 (C12) · P368, P369 (C13) · P370, P371, P372 (C14) · P373 (C15).
- **Visuals per CORE block:** C01 F1 + residual figure · C02 A1, A2, x–t figure, **E1** · C03 control-volume drawing + budget
  bars · C04 F2 + ruler bars · C05 F3 · C06 two figures, F4, **E2** · C07 $M_2$–$M_1$ figure, two-frame figure, F5a · C08 entropy
  figure, Hugoniot figure, profiles, F5b, **E3** · C09 sketch, pressure family, F6 · C10 regime figure, A3, live, **E4** · C11 F7 ·
  C12 A5, **E5** · C13 three panels, A4, **E6** · C14 triangles, wedge figures, F8, live, **E7** · C15 fan, curved wall, the ±δ plot,
  F9b, **E8** · C16 plate figure, slope figure, F9a, **E9**.
- **From-scratch cells (16):** one `nb.check_agree` row in every CORE block (curation §7).

---

## Part B — explainer storyboards

Common to all nine (and the backup): created with `tools/new_viz.py`; `<meta name="viz:chapter" content="ch15">`; tabs
Walkthrough · Explore · Explain · Derivation · Equations · Code · Check; every displayed number is computed by a JS function
that mirrors a `ch15` callable and is proved by `selftest()` parity rows (`py:` expressions as C.6; **never a typed constant
the page could compute** — 0.528, 0.378, 22.97° and 130.45° are computed). Explain is "Explanation & interpretation" in
numbered sections built with `Viz.work.step / line / box / table / hint / interpret`, modelled on
`forced_damped_vibrations.html` (the regime-dependent reading) and `amplitude_phase_second_order_II_3.html` (a numbered
derivation with live numbers, one section per optional view): **0** what the views show and what each colour means · **1…n**
every displayed quantity from the controls ("formula = substituted = result — why", results boxed) · a section or hint per
view hidden on phones · the values at the current time or cursor (live) · **Reading the current setting** (regime-dependent,
with the threshold of each regime). Derivation steps are copied from Part F (same `did` titles, same step count, same order;
phones shorten *why* to its first sentence; plain-text *why* and *watch* never contain raw TeX). Every tour, Explain,
Derivation, notes, status, equation and quiz text that names a book equation **writes it out** next to its number. Colours
as header convention 8 (pressure orange, velocity and Mach number blue, density teal, temperature rose, stagnation and
conserved quantities amber, entropy purple, ghosts and forbidden regions muted). **Every ratio is labelled with what it is a
ratio of.** Angles: sliders and readouts in degrees, arithmetic in radians, both shown in the Explain lines. Where the
book's symbols have a rival standard, section 0 of Explain shows both. Walkthrough texts ≤ 45 words, step 1 ≤ 24 words, ≤ 2
extras per step, `play: false` on steps that quote numbers. **Hygiene:** every power of a formatted number is bracketed; "the
book prints", never "the page prints"; slips in the house form; each control's extreme values are named in its help line;
where a formula stops being valid the curve ends with a mark and a label (no NaN text; "choked", "detached", "forbidden" are
words, never numbers). **JS physics is light:** closed forms and `Viz.num.brentq`, at most two nested; no tables, no complex
arithmetic, no linear solves. **Fit:** at 360×640 never three stacked views — two at a time, the third is `hidePortrait` and
its key number is repeated in a visible title or readout; at 844×345 the stage is one row of two views beside the text
panel, titles shortened, transport without step buttons. Parallel builders use private scratch subfolders.

### E1 · sound_pulse_dalembert
- **Title:** "What is a sound pulse made of?" · **Summary:** "Release a pressure bump and watch it split into two shapes
  sliding apart at ±c — then put it in a stream." · **CORE:** C01, C02 (also N07, N08, N17, N19, N20, N21) · **Reference:**
  `forced_damped_vibrations.html` (system + graph on one time slider; the Explain panel) and ch10 `upwind_cfl_advection` (an x–t
  view with characteristics).
- **meta:** `viz:order 1` · `viz:sections 15.2` · `viz:equations 15.5 15.7 15.12 15.13 15.15 15.16 15.17` · `viz:fluidpy
  ch15.dalembert ch15.dalembert_initial_value ch15.plane_wave_velocity ch15.dalembert_in_stream ch15.wave_speeds
  ch15.upstream_signal_possible ch15.condensation ch15.linearisation_error ch15.pulse` · `viz:derivations D03 D04 D05`.
- **Physics (JS ↔ Python):** lengths in pulse widths w, times in w/c, pressures in units of the bump height (readouts convert
  to SI for the chosen medium). `shape(s)` = exp(−s²) ↔ `ch15.pulse("gaussian", 1.0, 1.0)`; `parts(ic)` → {f, g} (bump at rest: ½
  shape each; right only: f = shape; left only: g = shape; two meeting: f centred at −3, g at +3) ↔ `ch15.dalembert_initial_value`
  for the first; `pPrime(x, t, M)` = f(x − (1 + M)t) + g(x + (1 − M)t) ↔ `ch15.dalembert_in_stream(x, t, c, U, f, g)`; `uPrime(x, t,
  M)` = (f − g)/ρ₀c ↔ `ch15.plane_wave_velocity`; `speeds(U, c)` ↔ `ch15.wave_speeds`; `upstream(U, c)` ↔
  `ch15.upstream_signal_possible`; `cond(p, rho0, c)` ↔ `ch15.condensation`; `linErr(p, rho0, c, gamma)` = ½(γ − 1) × condensation ↔ `ch15.linearisation_error(p_prime, rho0, c, gamma)`.
- **Views** (rows [1, 0.5, 1]): 1. `wave` "Pressure p′(x, t)" (row 0): x from −8 to 8 [w]; the sum bold (orange), f (blue) and g
  (teal) faint, a muted envelope of where each has been ("so far"), the ghost of the initial bump; pointer: click sets the
  probe x. 2. `air` "The air" (row 1): a strip of 60 tracer particles displaced by ∫u′dt with their u′ arrows (blue), drifting
  with the stream when M > 0; the probe marked. 3. `xt` "x–t diagram" (row 2, `hidePortrait`): time upward, the two families of
  characteristics through the bump's edges (slopes U + c blue, U − c teal), the present time as a horizontal cursor, the region
  that the bump has influenced shaded.
- **Controls:** `M` "Stream Mach number $U/c$" (0…2, step 0.05, default 0; help "0: still air · 1: the upstream wave stands
  still · 2: both swept downstream") · `ic` "Start" (chips: bump at rest · right only · left only · two meeting) · `s` "Amplitude
  (condensation $p'/\rho_0c^2$)" (10⁻⁶…0.3 log, default 10⁻³; help "10⁻⁶: a whisper · 0.3: far beyond acoustics"; optional) ·
  `medium` (chips: air 340 m/s · helium 999 m/s · Newton's isothermal air 288 m/s; optional) · `t` (transport).
- **Transport:** `t` 0…6 [w/c], rate 0.8/s, `end: 'hold'`; end-of-run card (`Viz.card`): "Each half carries half the bump. The
  shapes never changed — only where they are."
- **Presets:** "still air" {M 0, ic bump} · "M = 0.5" {M 0.5} · "M = 1: the upstream wave stands still" {M 1} · "M = 2: both go
  downstream" {M 2} · "two pulses meet" {M 0, ic two meeting}.
- **Status:** M < 1: "subsonic stream: news travels both ways (U − c = −170 m/s)" (number live) · M = 1 (±0.01): "sonic: the
  upstream wave stands still" · M > 1: "supersonic: nothing travels upstream (U − c = +340 m/s)".
- **Readouts:** "Sound speed $c$" · "Speeds $U+c$, $U-c$" · "p′ at probe" · "u′ at probe" · "Linearisation error".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** with an end-of-run card · **presets** ·
  **status** · **inspector**.
- **Inspector** (click a point of `wave`): "x = 1.20 w, t = 1.50 w/c: f = ½ exp(−(1.20 − 1.50)²) = **0.457**; g = ½ exp(−(1.20 +
  1.50)²) = **0.0003**; $p'(x,t)=f(x-ct)+g(x+ct)$ (15.15): p′ = 0.457 + 0.0003 = **0.457** (× 100 Pa); $u_1'(x,t)=\dfrac1{\rho_0c}\left(f(x-ct)-g(x+ct)\right)$ (15.16): u′ = (0.457 − 0.0003) × 100/416.9 = **0.110 m/s**
  (forward: a compression running right)".
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**Top:** the pressure disturbance. The <b class=c-orange>orange</b> curve is what a microphone
     measures; it is the sum of a <b class=c-blue>right-going shape f</b> and a <b class=c-teal>left-going shape g</b>.
     **Middle:** the air itself: dots are parcels, arrows their velocity u′. **Bottom:** the same story on a map of position
     against time; each shape travels along a straight line."
  1. *The sound speed of your medium* — "$c^2=\gamma p_0/\rho_0=\gamma RT_0$ (15.6): c = √(1.4 × 287.05 × 288.15) = **340.3 m/s** (helium: √(1.667 × 2077 × 288.15) = 999
     m/s; Newton's isothermal guess √(RT₀) = 287.6 m/s, 15.5 % low)."
  2. *Splitting your start into f and g* — "A bump F at rest: the initial-value formulas $f(x)=\dfrac12\left[F(x)-\dfrac1c\displaystyle\int_0^xG(\hat x)\,d\hat x\right],\ g(x)=\dfrac12\left[F(x)+\dfrac1c\displaystyle\int_0^xG(\hat x)\,d\hat x\right]$ (unnumbered, Exercise 15.1) with G = 0 give f = g = ½F.
     'Right only' puts all of F in f (and gives the air an initial push forward)."
  3. *Pressure and velocity at the probe* — the inspector's lines (live), ending in the boxed sum and the boxed difference.
  4. *In a stream* — "$p'(x,t)=f(x-(c+U)t)+g(x+(c-U)t)$ (15.17): the shapes ride at U + c = `…` and U − c = `…` m/s. At M = `…` the second is `negative: it still
     travels upstream` / `zero: it stands still` / `positive: it is swept downstream`."
  5. *How good is 'small'?* — "$p'/\rho_0c^2\ll1$ (15.7): condensation s = `…`; the dropped term changes $p'\cong c^2\rho'$ (15.5) by about ½(γ − 1)s = `…` %. At s = 0.1
     (174 dB) it is 2 %."
  6. *The x–t view* (or the hint "turn the phone sideways for the x–t diagram; the two speeds are in the title of the pressure
     view").
  7. *Now* — "t = `Viz.live('t')` w/c: the right-going shape is centred at `Viz.live('xf')`, the left-going one at
     `Viz.live('xg')`."
  8. *Reading the current setting* — still air: "Two half-height copies run apart at c. Under the right-going one the air
     moves forward, under the left-going one backward; at the start the two motions cancelled — the air was at rest." ·
     subsonic stream: "The stream carries both shapes at U, but one still beats it upstream: a body in this stream can warn
     the flow ahead." · M = 1: "The upstream shape is stuck. This is the throat of a choked nozzle: news from downstream gets
     this far and no farther." · supersonic: "Both shapes go downstream. Nothing here can influence the flow ahead — the
     root of choking, shocks and the Mach cone."
- **Derivation tab:** **D03** (8 steps) `view: 'wave'` on phones; goal `set` {preset still air, t 0}; step 2 `watch` "the two
  dashed lines x − ct = const and x + ct = const in the x–t view"; step 7 `set` {ic right only}, `live` "f at the probe = `…`";
  step 8 `set` {ic bump}, `play`. **D04** (5 steps) `view: 'air'` on phones; step 3 `live` "∂f/∂x = −(1/c)∂f/∂t at the probe: `…`
  = `…`"; step 5 `set` {ic two meeting, t 3}, `watch` "where the two pulses overlap the pressures add and the velocities cancel".
  **D05** (5 steps) `view: 'wave'` on phones; step 1 `set` {M 0.5}; step 4 `live` "U + c = `…`, U − c = `…` m/s"; step 5 `set` {M
  2}, `watch` "both shapes move right". Interpret: `s => "At M = " + … + " the slower shape moves at " + … + " m/s — " + (…)`.
- **Code:**
  ```python
  F = ch15.pulse("gaussian", 1.0, {{amp}})                       # the bump you released [Pa]
  f = lambda s: 0.5*F(s); g = lambda s: 0.5*F(s)               # at rest: half each
  p = ch15.dalembert_in_stream(x, {{t}}, c, {{U}}, f, g)        # (15.17): p' = {{p}} Pa at the probe
  u = ch15.plane_wave_velocity(x - {{U}}*{{t}}, {{t}}, c, rho0, f, g)   # (15.16): u' = {{u}} m/s
  c_plus, c_minus = ch15.wave_speeds({{U}}, c)                 # {{cp}}, {{cm}} m/s
  ch15.upstream_signal_possible({{U}}, c)                      # {{up}}
  ```
- **Walkthrough (7 steps):** 1. "A bump of pressure" — "You squeeze the air in one place and let go. What happens next?"
  `set` {preset still air, t 0} · 2. "It splits" — "Press play. Two copies, each half as high, run apart at the speed of sound
  and keep their shape." `play: true`, `eq: 'dal'` · 3. "Why two shapes" — "In the coordinates x − ct and x + ct the wave
  equation says: a function of one plus a function of the other." `derive: {id: 'D03', step: 6}`, `eq: 'wave1d'` · 4. "The air
  only shuffles" — "Click under each half. The air moves forward under one, backward under the other — a thousand times
  slower than the wave." `inspect: true`, `eq: 'u'` · 5. "Add a stream" — "Now the air flows at half the speed of sound. One
  shape still struggles upstream." `set` {M 0.5}, `play: true`, `eq: 'stream'` · 6. "Faster than sound" — "At Mach 2 both shapes
  are swept downstream. No news can reach the flow ahead." `set` {M 2}, `play: true`, `code: {id: 'dal', lines: [5, 6]}` · 7. "Your
  turn" — "Predict the Mach number at which the upstream shape stands still. Then drag." `controls: ['M', 'ic']`.
- **Equations:** `wave` "Wave equation in a stream" ref 'Eq. (15.12)' $\dfrac1{c^2}\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)^2p'-\dfrac{\partial^2p'}{\partial x_j\partial x_j}=0$ · `wave1d` "Still fluid, one dimension" ref 'Eq.
  (15.13)' $\dfrac1{c^2}\dfrac{\partial^2p'}{\partial t^2}-\dfrac{\partial^2p'}{\partial x^2}=0$ · `dal` "d'Alembert's solution" ref 'Eq. (15.15)' $p'(x,t)=f(x-ct)+g(x+ct)$, live at the probe · `u` "Particle velocity" ref 'Eq.
  (15.16)' $u_1'(x,t)=\dfrac1{\rho_0c}\left(f(x-ct)-g(x+ct)\right)$, live, note "the book writes $\rho_o$; $\rho_0c$ is the acoustic impedance (our name)" · `stream` "In a stream"
  ref 'Eq. (15.17)' $p'(x,t)=f(x-(c+U)t)+g(x+(c-U)t)$, live with the two speeds · `small` "When is it acoustics?" ref 'Eq. (15.7)' $p'/\rho_0c^2\ll1$ with $p'\cong c^2\rho'$ (15.5).
- **Check yourself:** (1) "A 100 Pa bump is released at rest. How high is each travelling pulse?" — "50 Pa: $f(x)=\dfrac12\left[F(x)-\dfrac1c\displaystyle\int_0^xG(\hat x)\,d\hat x\right],\ g(x)=\dfrac12\left[F(x)+\dfrac1c\displaystyle\int_0^xG(\hat x)\,d\hat x\right]$ (unnumbered, Exercise 15.1) with G = 0
  gives f = g = ½F." `set {ic: 'bump'}` · (2) "Under the left-going pulse, which way does the air move?" — "Backward (against
  +x): $u_1'(x,t)=\dfrac1{\rho_0c}\left(f(x-ct)-g(x+ct)\right)$ (15.16) has −g; a compression pushes the air the way the wave goes." · (3) "At what stream Mach number does the
  upstream-running shape stand still?" — "M = 1: its speed in $p'(x,t)=f(x-(c+U)t)+g(x+(c-U)t)$ (15.17) is U − c." `set {M: 1}` · (4) "Two equal pulses meet head-on.
  What are p′ and u′ at the instant they overlap?" — "p′ doubles, u′ is zero: f + g against f − g." `set {ic: 'meet', t: 3}`.
- **Selftest parity rows:** `{name: 'dalembert', js: pPrime(0.7, 1.2, 0), py: 'ch15.dalembert(0.7, 1.2, 1.0, f=lambda s:
  0.5*np.exp(-s*s), g=lambda s: 0.5*np.exp(-s*s))', rtol: 1e-12}` (0.40293) · `{name: 'u', js: uPrime(0.7, 1.2, 0), py:
  'ch15.plane_wave_velocity(0.7, 1.2, 1.0, 1.0, f=lambda s: 0.5*np.exp(-s*s), g=lambda s: 0.5*np.exp(-s*s))', rtol: 1e-12}` (0.37588) ·
  `{name: 'stream', js: pPrime(0.7, 1.2, 0.5), py: 'ch15.dalembert_in_stream(0.7, 1.2, 1.0, 0.5, f=lambda s: 0.5*np.exp(-s*s),
  g=lambda s: 0.5*np.exp(-s*s))', rtol: 1e-12}` (0.24136) · `{name: 'U-c', js: speeds(170, 340)[1], py: 'ch15.wave_speeds(170.0,
  340.0)[1]', rtol: 1e-12}` (−170) · `{name: 'lin err', js: linErr(0.1*1.225*340.3*340.3, 1.225, 340.3, 1.4), py: 'ch15.linearisation_error(0.1*1.225*340.3**2, 1.225, 340.3, 1.4)', rtol: 1e-12}` (0.0200) · `{name: 'condensation', js: cond(1, 1.225, 340.3), py: 'ch15.condensation(1.0, 1.225, 340.3)', rtol: 1e-12}` ·
  invariant `{name: 'halves', js: pPrime(3, 3, 0), expect: 0.5, rtol: 1e-6}`.
- **Fit plan:** 360×640: status (1 line) · `wave` (60 %) over `air` (40 %); `xt` hidden — the two speeds sit in the `wave` title
  ("U + c 510 · U − c −170 m/s"); the medium chips are optional. 844×345: `wave` beside `xt`, `air` as a thin strip under `wave`.
  Desktop / 1000×700: rows [1, 0.5, 1].

### E2 · area_mach_nozzle
- **Title:** "Why does a supersonic nozzle widen?" · **Summary:** "Move along a duct and watch density, speed and area trade
  places at Mach 1." · **CORE:** C05, C06 (also N24, N35, N36, N39, N41, N42, N43, N44) · **Reference:**
  `amplitude_phase_second_order_II_3.html` (linked windows with a crosshair, numbered Explain) and ch14
  `elliptic_loading_aspect_ratio` (a U-shaped curve with a marked minimum and a ghost).
- **meta:** `viz:order 2` · `viz:sections 15.4 15.5` · `viz:equations 15.19 15.31 15.32 15.33 15.34` · `viz:fluidpy ch15.area_ratio
  ch15.mach_from_area_ratio ch15.nozzle_isentropic ch15.laval_area ch15.influence_coefficients ch15.area_velocity ch15.duct_role
  ch15.mass_flux_function ch15.isentropic_state` · `viz:derivations D10 D12`.
- **Physics (JS ↔ Python):** `areaRatio(M, g)` ↔ `ch15.area_ratio`; `machFromArea(a, g, branch)` (`Viz.num.brentq` on [1e-8, 1] or
  [1, 100]) ↔ `ch15.mach_from_area_ratio(a, g, branch=…)`; `lavalArea(x, Ain, Aexit)` (two parabolas, C.6) ↔ `ch15.laval_area`;
  `solution(x, mode)` → M(x): mode `trans` uses $A^*=A_t$ (subsonic before the throat, supersonic after); `sub` and `sup` use
  $A^*=0.9A_t$ on one branch throughout ↔ `ch15.nozzle_isentropic(x, A_fn, g, branch=…, A_star=…)`; `state(M, g)` ↔
  `ch15.isentropic_state`; `coeffs(M, g)` → {drho, du, dA per unit dA/A} ↔ `ch15.influence_coefficients(M, g, dA_over_A=1.0)`;
  `areaVelocity(M)` ↔ `ch15.area_velocity(M, 1.0)`; `role(M, sign)` ↔ `ch15.duct_role` ("nozzle", "diffuser", "none", "sonic"); `flux(M, g)` ↔ `ch15.mass_flux_function`.
- **Views** (rows [1, 1]): 1. `duct` "The duct" (row 0, flex 1.5): the walls from `lavalArea`, tracer particles moving at the
  local speed (spacing shows density, colour shows M: blue subsonic → rose supersonic), the station cursor, a draggable throat
  handle and exit handle; labels "reservoir" at the left. 2. `curve` "$A/A^*$ against M" (row 1, flex 1): log M axis 0.1…5, the
  U-shaped curve with its minimum 1 at M = 1 (amber marker), the cursor dot (bold) and the **ghost dot** = the other Mach
  number with the same area (muted), a horizontal line joining them. 3. `bars` "The three relative changes" (row 1, flex 0.8,
  `hidePortrait`): signed bars of $d\rho/\rho$ (teal), $du/u$ (blue), $dA/A$ (muted) per unit step along the duct — they sum to
  zero; a thin fourth bar "sum".
- **Controls:** `x` "Station" (transport) · `Ae` "Exit / throat area $A_e/A_t$" (1.05…5, step 0.05, default 2; help "1.05: almost
  a pipe · 5: a long bell") · `mode` (modes: "subsonic throughout — a venturi" · "sonic throat — a Laval nozzle" · "supersonic
  throughout") · `gamma` (chips 1.4 · 5/3 · 1.3) · `Ain` "Inlet / throat area" (1.5…5, default 2.5; optional).
- **Transport:** `x` 0…1, rate 0.12/s, `end: 'hold'`; end-of-run card: "Converging then diverging — and the gas sped up the
  whole way."
- **Presets:** "venturi" {mode sub} · "design Laval nozzle" {mode trans, Ae 2} · "at the throat" {mode trans, x 0.4} · "far down
  the bell" {mode trans, x 0.95, Ae 5} · "helium" {gamma 5/3}.
- **Status:** M < 1 and dA < 0: "subsonic + converging: a nozzle" · M < 1 and dA > 0: "subsonic + diverging: a diffuser" · M > 1
  and dA > 0: "supersonic + diverging: still a nozzle" · M > 1 and dA < 0: "supersonic + converging: a diffuser" · throat with
  M = 1: "throat: M = 1 — the bars read 0/0" · throat with M ≠ 1: "throat, not sonic: the speed has a maximum or minimum here".
- **Readouts:** "Mach number $M$" · "$A/A^*$ here" · "Other M, same area" · "$du/u$ per $dA/A$" · "$p/p_0$ static/stagn.".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **modes** · **terms** (the three
  bars) · **presets** · **status**.
- **Terms:** title "Mass must balance: $\dfrac{d\rho}{\rho}+\dfrac{du}{u}+\dfrac{dA}{A}=0$"; items `drho` "density" teal, `du` "speed" blue, `dA` "area" muted; total "sum = 0".
- **Explain:**
  0. *What the views show* — "**Top:** the duct; dots are gas parcels — far apart where the gas is thin, <b class=c-blue>blue
     </b> below Mach 1, <b class=c-rose>rose</b> above. **Left:** the area a flow needs at each Mach number, in units of the
     sonic area $A^*$. **Right:** how density, speed and area change over the next small step."
  1. *The area ratio at your station* — "$A/A_t$ = `…`; with $A^*$ = `…` $A_t$ for this solution, $A/A^*$ = **`…`**."
  2. *Which Mach numbers fit this area?* — "$\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31) has two solutions: M = **`…`** (subsonic) and M = **`…`** (supersonic). This
     solution takes the `…` one because `we are upstream of a sonic throat` / `we are past it` / `the flow never reaches Mach
     1`."
  3. *Density against speed* — "$d\rho/\rho=-M^2\,du/u$ (15.33): with M = `…`, a 1 % gain in speed costs $M^2$ = `…` % of the density." (The book prints
     the line before it without a minus sign; the correct form is $u\,du=-dp/\rho=-c^2d\rho/\rho$ (15.32, third member with its minus sign restored).)
  4. *Speed against area* — "$\dfrac{du}{u}=-\dfrac{1}{1-M^2}\dfrac{dA}{A}$ (15.34): du/u = −1/(1 − `…`) × dA/A = **`…`** × dA/A."
  5. *The three bars add to zero* — boxed: "$\dfrac{d\rho}{\rho}+\dfrac{du}{u}+\dfrac{dA}{A}=0$ (15.19): `…` + `…` + `…` = 0".
  6. *The state here* — "$\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28), $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29), $\dfrac{\rho_0}{\rho}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{1/(\gamma-1)}$ (15.30): $T/T_0$ = `…`, $p/p_0$ = `…`, $\rho/\rho_0$ = `…` (static / stagnation)."
  7. *The bars* (or the hint "turn the phone sideways for the bars; $du/u$ per $dA/A$ is in the curve's title").
  8. *Reading the current setting* — subsonic converging: "The density bar is shorter than the speed bar ($M^2<1$): ρu rises
     as u rises, so the flow needs less room — and the wall gives it less." · at the sonic throat: "Density and speed bars are
     equal and opposite; the area bar is zero. This is the only place the flow can pass through Mach 1." · supersonic
     diverging: "The density bar is now the longer one ($M^2>1$): the gas thins faster than it speeds up, ρu falls, and it
     needs *more* room. The widening wall is what lets it keep accelerating." · venturi: "Same throat, but the flow is too
     slow to reach Mach 1 there, so the speed just peaks and falls back."
- **Derivation tab:** **D10** (9 steps) `view: 'curve'` on phones; goal `set` {preset design Laval nozzle, x 0.8}; step 4 `live`
  "$\rho_0/\rho$ = `…`, $\rho^*/\rho_0$ = `…`"; step 8 `live` "A/A* = (1/`…`) × (`…`)^3 = `…`"; step 9 `set` {x 0.4}, `watch` "the dot sits
  at the bottom of the U". **D12** (7 steps) `view: 'duct'` on phones; step 4 `live` "dρ/ρ = −(`…`)² du/u"; step 6 `live` "du/u =
  `…` × dA/A"; step 7 `set` {x 0.4}, `watch` "at the throat the area bar vanishes while the other two stay finite". Interpret:
  `s => "At M = " + … + " a 1 % wider duct changes the speed by " + … + " %."`
- **Code:**
  ```python
  A_over_Astar = {{AAt}}/{{Astar}}                              # area here over the sonic area
  M_sub = ch15.mach_from_area_ratio(A_over_Astar, {{g}}, branch="subsonic")     # {{Msub}}
  M_sup = ch15.mach_from_area_ratio(A_over_Astar, {{g}}, branch="supersonic")   # {{Msup}}
  st = ch15.isentropic_state({{M}}, {{g}})                      # p/p0 {{pp0}}, T/T0 {{TT0}}
  k = ch15.area_velocity({{M}}, 1.0)                            # (15.34): du/u = {{k}} x dA/A
  c = ch15.influence_coefficients({{M}}, {{g}}, dA_over_A=0.01) # drho/rho {{dr}}, du/u {{du}}
  ```
- **Walkthrough (7 steps):** 1. "Squeeze to go faster?" — "A gas flows from a reservoir through a narrowing, then widening
  duct. Where is it fastest?" `set` {preset design Laval nozzle, x 0} · 2. "Converging" — "Press play. Below Mach 1 the duct
  narrows and the gas speeds up — like a garden hose." `play: true` · 3. "The throat" — "At the narrowest point the gas reaches
  exactly Mach 1. The area curve is at its minimum." `set` {x 0.4}, `eq: 'area'` · 4. "Diverging — and still faster" — "Past the
  throat the duct widens and the gas keeps accelerating. Watch the density bar: it now beats the speed bar." `set` {x 0.75},
  `terms: true`, `eq: 'rho'` · 5. "The rule" — "One factor, 1 − M², decides the sign." `eq: 'av'`, `derive: {id: 'D12', step: 6}` · 6.
  "Same duct, slower flow" — "Lower the flow and the throat is no longer sonic: a venturi." `set` {mode sub}, `play: true` · 7. "Your
  turn" — "Predict the exit Mach number for an exit three times the throat. Then drag." `controls: ['Ae', 'mode']`.
- **Equations:** `mass` "Mass, logarithmic form" ref 'Eq. (15.19)' $\dfrac{d\rho}{\rho}+\dfrac{du}{u}+\dfrac{dA}{A}=0$, live with the three bars · `area` "Area ratio" ref 'Eq.
  (15.31)' $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$, live · `euler` "Euler, isentropic (corrected)" ref 'Eq. (15.32)' $u\,du=-dp/\rho=-c^2d\rho/\rho$, note "> the book prints the third
  member without its minus sign" · `rho` "Density against speed" ref 'Eq. (15.33)' $d\rho/\rho=-M^2\,du/u$, live · `av` "Area–velocity relation"
  ref 'Eq. (15.34)' $\dfrac{du}{u}=-\dfrac{1}{1-M^2}\dfrac{dA}{A}$, live · `sonic` "Sonic state" ref 'ours' $\dfrac{T^*}{T_0}=\dfrac{2}{\gamma+1},\ \dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)},\ \dfrac{\rho^*}{\rho_0}=\left(\dfrac{2}{\gamma+1}\right)^{1/(\gamma-1)}$.
- **Check yourself:** (1) "For $A/A^*$ = 2 and γ = 1.4, what are the two Mach numbers?" — "0.306 and 2.197: $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31) is U-shaped
  with its minimum at M = 1." `set {Ae: 2, x: 1}` · (2) "At Mach 2, by how much does the speed change if the duct widens by 1
  %?" — "+0.33 %: $\dfrac{du}{u}=-\dfrac{1}{1-M^2}\dfrac{dA}{A}$ (15.34) with 1 − M² = −3." · (3) "Can the flow be sonic where the duct is not narrowest?" — "No: $\dfrac{du}{u}=-\dfrac{1}{1-M^2}\dfrac{dA}{A}$ (15.34) with M = 1
  needs dA = 0." · (4) "In the venturi mode, is the throat sonic?" — "No: a throat need not be sonic; du = 0 there instead."
  `set {mode: 'sub', x: 0.4}`.
- **Selftest parity rows:** `{name: 'A/A*(2)', js: areaRatio(2, 1.4), py: 'ch15.area_ratio(2.0, 1.4)', rtol: 1e-12}` (1.6875) ·
  `{name: 'M sub', js: machFromArea(2, 1.4, 'subsonic'), py: 'ch15.mach_from_area_ratio(2.0, 1.4, branch="subsonic")', rtol: 1e-9}`
  (0.30590) · `{name: 'M sup', js: machFromArea(2, 1.4, 'supersonic'), py: 'ch15.mach_from_area_ratio(2.0, 1.4, branch="supersonic")',
  rtol: 1e-9}` (2.19720) · `{name: 'area-velocity', js: areaVelocity(2), py: 'ch15.area_velocity(2.0, 1.0)', rtol: 1e-12}` (0.33333) ·
  `{name: 'drho', js: coeffs(2, 1.4).drho*0.01, py: 'ch15.influence_coefficients(2.0, 1.4, dA_over_A=0.01)["drho_over_rho"]', rtol:
  1e-10}` (−0.013333) · `{name: 'flux max', js: flux(1, 1.4), py: 'ch15.mass_flux_function(1.0, 1.4)', rtol: 1e-12}` (0.68473) ·
  `{name: 'shape', js: lavalArea(0.7, 2.5, 2), py: 'ch15.laval_area(0.7, A_throat=1.0, A_inlet=2.5, A_exit=2.0, x_throat=0.4)', rtol: 1e-12}` (1.25) · `{name: 'shape in', js: lavalArea(0.2, 2.5, 2), py: 'ch15.laval_area(0.2)', rtol: 1e-12}` (1.375) · `{name: 'p/p0', js: state(2, 1.4).p_over_p0, py: 'ch15.isentropic_state(2.0, 1.4)["p_over_p0"]', rtol: 1e-12}`
  (0.12780) · invariant `{name: 'bars sum', js: coeffs(1.7, 1.4).drho + coeffs(1.7, 1.4).du + 1, expect: 0, atol: 1e-12}`.
- **Fit plan:** 360×640: status · `duct` (45 %) over `curve` (55 %); `bars` hidden — "du/u = +0.33 × dA/A" in the `curve` title;
  `Ain` optional. 844×345: `duct` beside `curve`. Desktop / 1000×700: rows [1, 1] with `curve` and `bars` side by side.

### E3 · normal_shock_jump
- **Title:** "Why can a shock only compress?" · **Summary:** "Drag the upstream Mach number below 1: three budgets still
  balance, entropy does not." · **CORE:** C07, C08 (also N49–N52, N54, N56, N57, N58, N60, N61, N62, N64, R16) · **Reference:**
  `forced_damped_vibrations.html` (the Explain panel; a solution split into parts) and ch07 `hydraulic_jump` (budget bars and
  a "forbidden" preset).
- **meta:** `viz:order 3` · `viz:sections 15.6` · `viz:equations 15.35 15.36 15.37 15.39 15.40 15.41 15.42 15.43 15.44` ·
  `viz:fluidpy ch15.normal_shock ch15.normal_shock_mach ch15.normal_shock_residuals ch15.normal_shock_fluxes
  ch15.normal_shock_entropy ch15.weak_shock_entropy ch15.normal_shock_stagnation_pressure_ratio ch15.shock_admissible
  ch15.moving_shock ch15.hugoniot_pressure_ratio ch15.hydraulic_jump` · `viz:derivations D13 D17`.
- **Physics (JS ↔ Python):** `shock(M1, g)` → {M2, p, rho, T, p0, dS, admissible} ↔ `ch15.normal_shock(M1, g,
  allow_expansion=True)`; `fluxes(M1, g)` → the pressure and momentum shares of $p+\rho u^2$ and the enthalpy and kinetic shares of $h+\tfrac12u^2$ on each side (`p1_part`, `ke1_part`, `p2_part`, `ke2_part`, `h1_part`, `k1_part`, `h2_part`, `k2_part`) ↔ `ch15.normal_shock_fluxes`; `dS(M1, g)` ↔ `ch15.normal_shock_entropy`; `dSweak(M1, g)`
  ↔ `ch15.weak_shock_entropy(M1, g, form="mach")`; `p0ratio(M1, g)` ↔ `ch15.normal_shock_stagnation_pressure_ratio`; `moving(p1,
  T1, M1)` → {u_shock, u_gas} ↔ `ch15.moving_shock(p1, T1, M_s=M1)`; `hugoniot(r, g)` ↔ `ch15.hugoniot_pressure_ratio`; `jump(Fr1)`
  = ½(−1 + √(1 + 8Fr₁²)) ↔ `ch15.hydraulic_jump(1.0, Fr1=…)` ($H_2/H_1=\tfrac12\left(-1+\sqrt{1+8Fr_1^2}\right)$ (7.81)).
- **Views** (rows [0.9, 1.1]): 1. `tube` "The shock in its own frame" (row 0, flex 1.3): particles arriving fast and sparse
  from the left, leaving slow and dense, coloured by temperature (rose scale); gauges for p, T, ρ on each side; in mode
  `lab` the shock front moves right into still gas and the gas behind follows; in mode `water` a water surface with a
  hydraulic jump. For $M_1<1$ the picture is drawn muted with a "forbidden" stamp. 2. `ratios` "Downstream / upstream against
  $M_1$" (row 1, flex 1.2): $p_2/p_1$ (orange, right axis), $\rho_2/\rho_1$ (teal), $T_2/T_1$ (rose), $M_2$ (blue),
  $T_{02}/T_{01}$ flat at 1 (amber), $p_{02}/p_{01}$ (amber dashed, falling); the cursor; the half $M_1<1$ shaded muted; the
  density limit $(\gamma+1)/(\gamma-1)$ dashed. 3. `entropy` "Entropy and the three budgets" (row 1, flex 1, `hidePortrait`):
  $(S_2-S_1)/C_v$ against $M_1$ (purple) with the cubic $\dfrac{S_2-S_1}{C_v}\cong\dfrac{2\gamma(\gamma-1)}{3(\gamma+1)^2}\left(M_1^2-1\right)^3$ as a muted ghost; below, three pairs of flux bars (mass,
  momentum, energy: before = after, amber) and one signed entropy bar (purple; red outline when negative).
- **Controls:** `M1` "Upstream Mach number $M_1$" (0.6…5, step 0.01, default 2; help "0.6: forbidden half · 1: no shock · 5:
  very strong") · `gamma` (chips 1.4 · 5/3 · 1.3) · `frame` (modes: "shock at rest" · "shock running into still gas" · "water: a
  hydraulic jump") · `p1`, `T1` "Upstream pressure, temperature" (optional; defaults 101 325 Pa, 288.15 K).
- **Presets:** "barely supersonic" {M1 1.1} · "Mach 2" {M1 2} · "Mach 5: density stuck below 6" {M1 5} · "forbidden: $M_1$ = 0.7"
  {M1 0.7} · "seen from the still gas" {M1 2, frame lab}.
- **Status:** $M_1>1$: "✅ allowed: entropy rises; stagnation pressure 72 % kept" (number live) · $M_1=1$ (±0.005): "no shock:
  nothing changes" · $M_1<1$: "⛔ forbidden: mass, momentum and energy balance — entropy would fall".
- **Readouts:** "$M_2$ downstream" · "$p_2/p_1$ down/up" · "$\rho_2/\rho_1$ down/up" · "$p_{02}/p_{01}$ down/up" · "$\Delta S/C_v$".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (three frames / systems) · **terms** (flux
  pairs and a signed entropy bar) · **presets** · **status**.
- **Explain:**
  0. *What the views show* — "**Top:** gas crossing a shock from left (station 1, upstream) to right (station 2,
     downstream). Other books number these 0 and 1. **Left graph:** every ratio is downstream / upstream. <b
     class=c-amber>Amber</b> = stagnation quantities: flat if conserved. **Right:** <b class=c-purple>entropy</b>, and three
     pairs of bars that must match."
  1. *Pressure* — "$\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39) = 1 + (2 × 1.4/2.4)(`…`² − 1) = **`…`**."
  2. *Density and speed* — "$\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}$ (15.41) = 2.4 × `…`/(0.4 × `…` + 2) = **`…`** (limit 6 for γ = 1.4)."
  3. *Temperature* — "$\dfrac{T_2}{T_1}=1+\dfrac{2(\gamma-1)}{(\gamma+1)^2}\dfrac{\gamma M_1^2+1}{M_1^2}\left(M_1^2-1\right)$ (15.42) = **`…`**; check with the gas law: `…`/`…` = `…` ✓."
  4. *Downstream Mach number* — "$M_2^2=\dfrac{(\gamma-1)M_1^2+2}{2\gamma M_1^2+1-\gamma}$ (15.40): $M_2^2$ = `…`/`…`, $M_2$ = **`…`** (never below √((γ − 1)/2γ) = 0.378)."
  5. *The three budgets, both sides* — a `Viz.work.table`: mass $\rho_1u_1=\rho_2u_2$, momentum $p_1-p_2=-\rho_1u_1^2+\rho_2u_2^2$ rearranged as $p+\rho u^2$, energy $h_1+\tfrac12u_1^2=h_2+\tfrac12u_2^2$
     — before, after, difference (zeros).
  6. *Entropy* — "$\dfrac{S_2-S_1}{C_v}=\ln\left[\dfrac{p_2}{p_1}\left(\dfrac{\rho_1}{\rho_2}\right)^{\gamma}\right]$ (15.43) = ln(`…` × `…`^(−1.4)) = **`…`**; the weak-shock cubic $\dfrac{S_2-S_1}{C_v}\cong\dfrac{2\gamma(\gamma-1)}{3(\gamma+1)^2}\left(M_1^2-1\right)^3$ (15.44a) gives `…` (gap `…` %)."
  7. *What is lost* — "$\dfrac{p_{02}}{p_{01}}=\left[\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}\right]^{\gamma/(\gamma-1)}\left[\dfrac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^{1/(\gamma-1)}$ (ours) = **`…`** = exp(−ΔS/R) = exp(−`…`/0.4) = `…` ✓. $T_0$ is unchanged."
  8. *In the frame of the still gas* (mode `lab`) — "subtract $u_1$: the shock runs at $M_1c_1$ = `…` m/s; the gas behind it
     follows at $u_1-u_2$ = `…` m/s. p, ρ, T, S did not change; M, $T_0$, $p_0$ did."
  9. *The entropy view* (or the hint "turn the phone sideways; ΔS/C_v and its sign are in the title of the ratio view").
  10. *Reading the current setting* — $M_1>1.5$: "A strong shock: the pressure rises `…`-fold but the density can never exceed
      `…` times. `…` % of the stagnation pressure is gone for good." · $1<M_1\le1.5$: "A weak shock: the loss grows only as the
      cube of $M_1^2-1$ — at Mach 1.1 it is a tenth of a per cent. This is why weak waves may be treated as isentropic." ·
      $M_1<1$: "Look at the bars: mass, momentum and energy all balance. Only the entropy bar is negative. A subsonic flow
      cannot jump to supersonic — not because of conservation, but because of the second law." · water mode: "Froude number
      plays the Mach number's part. The analogy is partial: the jump loses mechanical energy; the shock loses stagnation
      pressure, not energy."
- **Derivation tab:** **D13** (13 steps) `view: 'tube'` on phones; goal `set` {preset Mach 2}; step 3 `live` "$p_2/p_1$ = 1 + γ
  M₁²(1 − r) with r = `…`"; step 8 `set` {M1 1}, `watch` "at r = 1 nothing changes: the root we must divide out"; step 10 `live`
  "r = (0.4 + 2/`…`²)/2.4 = `…`"; step 13 `set` {preset Mach 2}, `live` "1 + 1.1667 × (`…` − 1) = `…`". **D17** (11 steps) `view:
  'ratios'` on phones (the entropy value in its title); goal `set` {preset barely supersonic}; step 6 `watch` "the first-order
  terms cancel"; step 8 `live` "0.0648 × (`…`)³ = `…` against exact `…`"; step 11 `live` "pressure jump `…`: (0.0408) × (`…`)³".
  Interpret: `s => "At M1 = " + … + " the cubic is " + … + " % from the exact entropy jump."`
- **Code:**
  ```python
  s = ch15.normal_shock({{M1}}, {{g}}, allow_expansion=True)    # M2 {{M2}}, p2/p1 {{p}}, rho2/rho1 {{r}}
  s["T_ratio"], s["p0_ratio"]                                # {{T}}, {{p0}}
  ch15.normal_shock_entropy({{M1}}, {{g}})                      # (15.43): dS/Cv = {{dS}}
  ch15.weak_shock_entropy({{M1}}, {{g}}, form="mach")           # (15.44a): {{dSw}}
  ch15.shock_admissible({{M1}})                                 # {{ok}}
  ch15.normal_shock_fluxes({{M1}}, {{g}})                       # mass, momentum, energy: equal pairs
  ```
- **Walkthrough (8 steps):** 1. "A jump" — "Fast, thin gas arrives; slow, dense, hot gas leaves. What fixes the state behind?"
  `set` {preset Mach 2}, `play: true` · 2. "Three budgets" — "Mass, momentum and energy that go in must come out. Three pairs of
  bars, always equal." `terms: true`, `eq: 'laws'` · 3. "One number decides" — "Every ratio depends on the upstream Mach number
  alone." `eq: 'p'`, `derive: {id: 'D13', step: 13}` · 4. "A ceiling" — "At Mach 5 the pressure is 29 times higher, but the density
  is stuck below 6." `set` {preset Mach 5}, `eq: 'rho'` · 5. "Run it backwards" — "Drag to Mach 0.7. The three budgets still
  balance. Find the bar that does not." `set` {preset forbidden}, `eq: 'S'` · 6. "The price" — "Energy is conserved. What is lost
  is stagnation pressure." `set` {preset Mach 2}, `eq: 'p0'`, `code: {id: 'shock', lines: [2, 3]}` · 7. "Weak shocks are cheap" —
  "Near Mach 1 the entropy curve lies flat on its cubic." `set` {preset barely supersonic}, `derive: {id: 'D17', step: 9}` · 8. "Your
  turn" — "Predict the Mach number at which a tenth of the stagnation pressure is lost. Then drag." `controls: ['M1', 'gamma']`.
- **Equations:** `laws` "The three jump laws" ref 'Eq. (15.35)–(15.37)' $\rho_1u_1=\rho_2u_2$, $p_1-p_2=-\rho_1u_1^2+\rho_2u_2^2$, $h_1+\tfrac12u_1^2=h_2+\tfrac12u_2^2$ · `p` "Pressure jump" ref 'Eq.
  (15.39)' $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$, live · `M2` "Downstream Mach number" ref 'Eq. (15.40)' $M_2^2=\dfrac{(\gamma-1)M_1^2+2}{2\gamma M_1^2+1-\gamma}$, live · `rho` "Density jump" ref 'Eq. (15.41)'
  $\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}$, live · `T` "Temperature jump" ref 'Eq. (15.42)' $\dfrac{T_2}{T_1}=1+\dfrac{2(\gamma-1)}{(\gamma+1)^2}\dfrac{\gamma M_1^2+1}{M_1^2}\left(M_1^2-1\right)$ · `S` "Entropy jump" ref 'Eq. (15.43)' $\dfrac{S_2-S_1}{C_v}=\ln\left[\dfrac{p_2}{p_1}\left(\dfrac{\rho_1}{\rho_2}\right)^{\gamma}\right]$, live, with $\dfrac{S_2-S_1}{C_v}\cong\dfrac{2\gamma(\gamma-1)}{3(\gamma+1)^2}\left(M_1^2-1\right)^3$ (15.44a) ·
  `p0` "Stagnation-pressure ratio" ref 'ours' $\dfrac{p_{02}}{p_{01}}=\left[\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}\right]^{\gamma/(\gamma-1)}\left[\dfrac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^{1/(\gamma-1)}$, note "$=e^{-(S_2-S_1)/R}=A_1^*/A_2^*$".
- **Check yourself:** (1) "At $M_1$ = 2, what fraction of the stagnation pressure survives?" — "72 %: $\dfrac{p_{02}}{p_{01}}=\left[\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}\right]^{\gamma/(\gamma-1)}\left[\dfrac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^{1/(\gamma-1)}$ (ours) = 0.721." `set {M1:
  2}` · (2) "Which of the three conservation laws fails at $M_1$ = 0.7?" — "None. Only the entropy jump $\dfrac{S_2-S_1}{C_v}=\ln\left[\dfrac{p_2}{p_1}\left(\dfrac{\rho_1}{\rho_2}\right)^{\gamma}\right]$ (15.43) is negative." `set
  {M1: 0.7}` · (3) "How dense can a shock make air?" — "Less than 6 times: the limit of $\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}$ (15.41) is (γ + 1)/(γ − 1)." `set {M1: 5}` · (4)
  "In the frame of the still gas, did $p_2/p_1$ change?" — "No: static p, ρ, T do not depend on the frame; M, $T_0$, $p_0$ do."
  `set {frame: 'lab'}`.
- **Selftest parity rows:** `{name: 'M2', js: shock(2, 1.4).M2, py: 'ch15.normal_shock_mach(2.0, 1.4)', rtol: 1e-12}` (0.57735) ·
  `{name: 'p2/p1', js: shock(2, 1.4).p, py: 'ch15.normal_shock_pressure_ratio(2.0, 1.4)', rtol: 1e-12}` (4.5) · `{name: 'T2/T1', js:
  shock(2, 1.4).T, py: 'ch15.normal_shock_temperature_ratio(2.0, 1.4)', rtol: 1e-12}` (1.6875) · `{name: 'dS', js: dS(2, 1.4), py:
  'ch15.normal_shock_entropy(2.0, 1.4)', rtol: 1e-10}` (0.130916) · `{name: 'dS weak', js: dSweak(1.1, 1.4), py:
  'ch15.weak_shock_entropy(1.1, 1.4, form="mach")', rtol: 1e-12}` (6.0025 × 10⁻⁴) · `{name: 'p02/p01', js: p0ratio(2, 1.4), py:
  'ch15.normal_shock_stagnation_pressure_ratio(2.0, 1.4)', rtol: 1e-12}` (0.720874) · `{name: 'forbidden dS', js: dS(0.7, 1.4), py:
  'ch15.normal_shock(0.7, 1.4, allow_expansion=True)["dS_over_cv"]', rtol: 1e-9}` (−0.029542) · `{name: 'M_s', js: Math.sqrt(1 +
  2.4/2.8), py: 'ch15.moving_shock(101325.0, 288.15, p2_over_p1=2.0)["M_s"]', rtol: 1e-12}` (1.36277) · `{name: 'hugoniot', js:
  hugoniot(8/3, 1.4), py: 'ch15.hugoniot_pressure_ratio(8/3, 1.4)', rtol: 1e-10}` (4.5) · invariants `{name: 'exp(-dS/R)', js:
  Math.exp(-dS(2, 1.4)/0.4), expect: 0.720874, rtol: 1e-6}`, `{name: 'jump Fr 2', js: jump(2), expect: 2.372281, rtol: 1e-6}`.
- **Fit plan:** 360×640: status (1 line) · `tube` (40 %) over `ratios` (60 %); `entropy` hidden — "ΔS/C_v = +0.131 · p₀ kept 72 %"
  in the `ratios` title; `p1`, `T1` optional. 844×345: `ratios` beside `entropy`, `tube` as a strip. Desktop / 1000×700: rows
  [0.9, 1.1].


### E4 · nozzle_back_pressure
- **Title:** "Why does a nozzle stop listening?" · **Summary:** "Open the valve: the mass flow rises, then sticks, while a
  shock is born, walks down the bell and leaves." · **CORE:** C09, C10 (also N39, N72–N78) · **Reference:**
  `angular_frequency_explorer_1.html` (system + graph on one clock, presets, regime-dependent notes, end-of-run summary) and
  ch14 `kutta_condition` (a quantity that changes character at one value).
- **meta:** `viz:order 4` · `viz:sections 15.7` · `viz:equations 15.29 15.31 15.39` · `viz:fluidpy ch15.convergent_nozzle
  ch15.mass_flow_vs_back_pressure ch15.choked_mass_flow ch15.critical_pressure_ratio ch15.laval_critical_pressures
  ch15.laval_regime ch15.laval_shock_location ch15.laval_nozzle ch15.laval_exit_state ch15.nozzle_with_shock` · `viz:derivations D19 D20`.
- **Physics (JS ↔ Python):** `crit(g)` ↔ `ch15.critical_pressure_ratio`; `convFlux(pB, g)` → {flux, Me, pe, choked} ↔ `ch15.mass_flow_vs_back_pressure(pB, g)` (keys `mdot` — dimensionless —, `M_exit`, `p_exit_over_p0`, `choked`, `mdot_over_choked`); `thresholds(ar, g)` → {p_b, p_d, p_f} ↔ `ch15.laval_critical_pressures`; `regime(pB, ar, g)` ↔ `ch15.laval_regime` ("subsonic", "sonic throat", "shock in nozzle", "shock at exit", "overexpanded", "design", "underexpanded"); `shockAt(pB, ar, g)` (outer `brentq` on $M_1$ in (1, $M_{\text{sup}}$); inner `brentq` for the
  subsonic inverse of the area ratio) → {M1, As, r, Me} ↔ `ch15.laval_shock_location`; `exitState(pB, ar, g)` ↔
  `ch15.laval_exit_state`; `curve(x, pB, ar, g)` → p(x)/p₀ ↔ `ch15.laval_nozzle(x, A_fn, pB, g)["p_over_p0"]` on `lavalArea`;
  `mdotMax(p0, T0, At)` ↔ `ch15.choked_mass_flow`.
- **Views** (rows [0.8, 1.2]): 1. `nozzle` "Reservoir → nozzle → tank" (row 0): the walls (convergent only, or
  convergent–divergent), particles coloured by M, the shock as a vertical rose line at its station, gauges $p_0$ (left) and
  $p_B$ (right); beyond the exit a **labelled sketch** of the jet (necking in with crossed lines when over-expanded, bulging
  with fans when under-expanded; caption "sketch, not a computation"). 2. `press` "$p(x)/p_0$ static / supply" (row 1, flex
  1.3): the current curve bold (orange), the two isentropic curves as muted ghosts (all-subsonic at $p_b$, design at $p_f$),
  the shock as a vertical jump, the critical level $\dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)}$ dashed, the back pressure as a tick on the right axis. 3. `flow`
  "Mass flow against back pressure" (row 1, flex 1, `hidePortrait`): $\dot m/\dot m_{\max}$ against $p_B/p_0$ (falling to the
  right), its plateau, the three dividing pressures $p_b$, $p_d$, $p_f$ as labelled verticals, the cursor.
- **Controls:** `pB` "Back pressure $p_B/p_0$" (0.02…1, step 0.002, default 1; also the transport) · `ar` "Exit / throat area
  $A_e/A_t$" (1.05…4, step 0.05, default 2; help "1.05: barely diverging · 4: a long bell") · `kind` (modes: "convergent only" ·
  "convergent–divergent") · `p0`, `T0` "Supply pressure, temperature" (100…1000 kPa, 250…1500 K; optional; they move the
  plateau in kg/s — trap T10) · `gamma` (optional).
- **Transport:** `pB` from 1 down to 0.02, rate 0.08/s, `end: 'hold'`; end-of-run card: "The mass flow stopped rising at
  $p_B/p_0$ = 0.937. After that the back pressure only moved the shock."
- **Presets** (the seven named conditions for $A_e/A_t$ = 2, values computed): "subsonic throughout" {pB 0.97} · "just choked"
  {pB = p_b} · "shock in the bell" {pB 0.7} · "shock at the exit" {pB = p_d} · "over-expanded" {pB 0.3} · "design" {pB = p_f} ·
  "under-expanded" {pB 0.05}; and "just choked (convergent)" {kind conv, pB = crit}.
- **Status:** the regime name with one number: "subsonic throughout: 67 % of the choked flow" · "choked; shock at 62 % of the
  bell, $M_1$ = 1.86" · "over-expanded: oblique shocks outside" · "design: shock-free, exit Mach 2.20" · "under-expanded:
  expansion fans outside"; convergent: "not choked: exit Mach 0.57" / "choked: exit pressure 0.528 $p_0$ > back pressure".
- **Readouts:** "$\dot m/\dot m_{\max}$" · "Exit Mach number" · "Exit pressure $/p_0$" · "Shock at $A/A_t$" · "$p_{02}/p_{01}$
  down/up".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** with an end-of-run card · **modes** ·
  **presets** · **status** · **notes**.
- **Notes** (regime-dependent): not choked: "The exit is subsonic, so the exit pressure equals the back pressure: the nozzle
  hears every change you make." · choked: "**'Choked' means**: the throat is sonic, so news of the back pressure cannot get
  past it; the mass flow is the most this throat passes *for this supply*. **It does not mean** the flow cannot be raised —
  raise $p_0$ or lower $T_0$ and it grows as $p_0/\sqrt{T_0}$."
- **Explain:**
  0. *What the views show* — "**Top:** a reservoir at pressure $p_0$ (the stagnation pressure) feeds a nozzle that empties
     into a tank at the back pressure $p_B$. **Left graph:** <b class=c-orange>pressure along the nozzle</b>, static over
     supply. **Right graph:** the mass flow as a fraction of the most the throat can pass."
  1. *The critical ratio* — "$\dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)}$ (ours) = (2/2.4)^3.5 = **0.528**: below it a sonic exit (or throat) cannot follow the back
     pressure."
  2. *The three dividing pressures for your area ratio* — "$\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31) with $A_e/A_t$ = `…` has the roots $M_{\text{sub}}$ = `…` and
     $M_{\text{sup}}$ = `…`. $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29) at each: $p_b/p_0$ = **`…`**, $p_f/p_0$ = **`…`**. A normal shock at the exit multiplies the
     design pressure by $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39) = `…`: $p_d/p_0$ = **`…`**."
  3. *Which regime?* — a `Viz.work.table` of the seven conditions with their ranges, the current row highlighted.
  4. *The shock, if it is inside* — "Find $M_1$ so that $\dfrac{p_B}{p_0}=\dfrac{p_{02}}{p_{01}}(M_1)\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sub}}\left(\frac{A_e}{A_t}\frac{p_{02}}{p_{01}}\right)}$ (ours) holds: $M_1$ = **`…`**; it stands where $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31) gives
     $A/A_t$ = **`…`**; $\dfrac{p_{02}}{p_{01}}=\left[\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}\right]^{\gamma/(\gamma-1)}\left[\dfrac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^{1/(\gamma-1)}$ (ours) = `…`; the sonic area behind it is 1/`…` = `…` $A_t$; the exit is subsonic at M = `…`."
  5. *The mass flow* — "$\dot m_{\max}=p_0A^*\sqrt{\dfrac{\gamma}{RT_0}}\left(\dfrac{2}{\gamma+1}\right)^{(\gamma+1)/(2(\gamma-1))}$ (ours) = `…` kg/s for your supply; now $\dot m/\dot m_{\max}$ = **`…`**" (convergent, not choked:
     $\dfrac{\dot m\sqrt{RT_0}}{p_0A}=\sqrt\gamma\,M\left(1+\dfrac{\gamma-1}{2}M^2\right)^{-(\gamma+1)/(2(\gamma-1))}$ (ours) with $M_e$ = `…` gives `…`/0.685).
  6. *The mass-flow view* (or the hint "turn the phone sideways; the flow fraction and the regime are in the title of the
     pressure view").
  7. *Reading the current setting* — subsonic throughout: "Nothing is sonic anywhere. Lower the back pressure and more
     flows." · shock inside: "The throat is sonic and the flow is at its ceiling. The supersonic stream after the throat is
     at too low a pressure to meet the tank, so it jumps through a normal shock and finishes as a subsonic diffuser flow.
     Lowering $p_B$ moves the shock toward the exit; nothing upstream of it changes." · over-expanded: "The shock has left.
     Inside, the nozzle is at its design flow, but the exit pressure is *below* the back pressure, so the jet is squeezed by
     oblique shocks outside." · design: "Exit pressure = back pressure. No waves. This is what the area ratio was built for." ·
     under-expanded: "The exit pressure is above the back pressure; the jet finishes expanding outside through fans. Nothing
     inside can change any more."
- **Derivation tab:** **D19** (7 steps) `view: 'press'` on phones; goal `set` {kind conv, pB 0.8}; step 3 `live` "$M_e$ = `…` from
  $p_B/p_0$ = `…`"; step 5 `set` {pB = crit}, `watch` "the exit reaches Mach 1 as the curve touches the dashed level"; step 7 `set`
  {pB 0.3}, `watch` "the curve no longer moves". **D20** (9 steps) `view: 'press'` on phones; goal `set` {preset shock in the
  bell}; step 2 `set` {preset just choked}; step 3 `set` {preset design}; step 4 `set` {preset shock at the exit}, `live` "`…` × `…` =
  `…`"; step 8 `live` "exit pressure for $M_1$ = `…`: `…` × `…` = `…`"; step 9 `live` "$A/A_t$ at the shock = `…`". Interpret: `s =>
  "With p_B/p_0 = " + … + " the shock costs " + … + " % of the stagnation pressure."`
- **Code:**
  ```python
  th = ch15.laval_critical_pressures({{ar}}, {{g}})            # p_b {{pb}}, p_d {{pd}}, p_f {{pf}}
  ch15.laval_regime({{pB}}, {{ar}}, {{g}})                     # '{{regime}}'
  loc = ch15.laval_shock_location({{pB}}, {{ar}}, {{g}})       # M1 {{M1}}, A_shock/A_t {{As}}
  ex = ch15.laval_exit_state({{pB}}, {{ar}}, {{g}})            # M_exit {{Me}}, mdot/choked {{frac}}
  ch15.choked_mass_flow({{p0}}, {{T0}}, {{At}})                # {{mdot}} kg/s
  ch15.mass_flow_vs_back_pressure({{pB}}, {{g}})["mdot"]       # convergent nozzle: flux {{flux}}
  ```
- **Walkthrough (8 steps):** 1. "Open the valve" — "A full reservoir, a nozzle, a tank you can pump down. Does lower pressure
  always mean more flow?" `set` {kind conv, pB 1} · 2. "More, and more…" — "Press play. The mass-flow needle climbs as the back
  pressure falls." `play: true`, `eq: 'flux'` · 3. "…then no more" — "At 0.528 the exit is sonic. News of the tank cannot get
  in; the needle stops." `set` {kind conv, pB 0.528}, `eq: 'crit'`, `derive: {id: 'D19', step: 6}` · 4. "Add a bell" — "With a
  diverging part the flow can go supersonic after the throat." `set` {kind cd, preset just choked} · 5. "A shock is born" — "Lower
  the pressure: a normal shock appears at the throat and walks downstream." `set` {preset shock in the bell}, `play: true`,
  `notes: true` · 6. "Where it stands" — "The back pressure picks the one shock position whose subsonic outflow ends exactly at
  $p_B$." `eq: 'pos'`, `derive: {id: 'D20', step: 8}` · 7. "Out through the exit" — "Below 0.513 the shock is gone and the inside
  is frozen at its design flow." `set` {preset design}, `code: {id: 'laval', lines: [1, 2]}` · 8. "Your turn" — "Predict what
  happens to the plateau if you double the supply pressure. Then try." `controls: ['pB', 'ar', 'p0']`.
- **Equations:** `crit` "Critical pressure ratio" ref 'ours' $\dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)}$, live · `flux` "Mass flux function" ref 'ours' $\dfrac{\dot m\sqrt{RT_0}}{p_0A}=\sqrt\gamma\,M\left(1+\dfrac{\gamma-1}{2}M^2\right)^{-(\gamma+1)/(2(\gamma-1))}$,
  live · `mmax` "Choked mass flow" ref 'ours' $\dot m_{\max}=p_0A^*\sqrt{\dfrac{\gamma}{RT_0}}\left(\dfrac{2}{\gamma+1}\right)^{(\gamma+1)/(2(\gamma-1))}$, live · `three` "The three dividing pressures" ref 'ours' $\dfrac{p_b}{p_0}=\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sub}}(A_e/A_t)},\ \dfrac{p_f}{p_0}=\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sup}}(A_e/A_t)},\ \dfrac{p_d}{p_0}=\dfrac{p_f}{p_0}\left.\dfrac{p_2}{p_1}\right\rvert_{M_{\text{sup}}}$, live ·
  `pos` "The shock's position" ref 'ours' $\dfrac{p_B}{p_0}=\dfrac{p_{02}}{p_{01}}(M_1)\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sub}}\left(\frac{A_e}{A_t}\frac{p_{02}}{p_{01}}\right)}$, note "then $A_{\text{shock}}/A_t$ from $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31) at $M_1$" · `iso` "Building
  blocks" ref 'Eq. (15.29), (15.31), (15.39)' $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$, $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$, $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$.
- **Check yourself:** (1) "For γ = 1.4, at what $p_B/p_0$ does a convergent nozzle choke?" — "0.528: $\dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)}$ (ours)." `set {kind: 'conv',
  pB: 0.528}` · (2) "With $A_e/A_t$ = 2 and $p_B/p_0$ = 0.7, where is the shock?" — "Where the area is 1.51 throats ($M_1$ =
  1.86)." `set {pB: 0.7}` · (3) "Does the mass flow change between $p_B/p_0$ = 0.7 and 0.05?" — "No: the throat is sonic in both;
  $\dot m_{\max}=p_0A^*\sqrt{\dfrac{\gamma}{RT_0}}\left(\dfrac{2}{\gamma+1}\right)^{(\gamma+1)/(2(\gamma-1))}$ (ours) depends on the supply only." · (4) "The nozzle is choked. Name two ways to pass more gas." — "Raise $p_0$ or lower
  $T_0$ (or enlarge the throat): $\dot m\propto p_0A^*/\sqrt{T_0}$." `set {p0: 1000}`.
- **Selftest parity rows:** `{name: 'crit', js: crit(1.4), py: 'ch15.critical_pressure_ratio(1.4)', rtol: 1e-12}` (0.528282) ·
  `{name: 'p_b', js: thresholds(2, 1.4).p_b, py: 'ch15.laval_critical_pressures(2.0, 1.4)["p_b"]', rtol: 1e-9}` (0.93716) · `{name:
  'p_d', js: thresholds(2, 1.4).p_d, py: 'ch15.laval_critical_pressures(2.0, 1.4)["p_d"]', rtol: 1e-9}` (0.51340) · `{name: 'p_f', js:
  thresholds(2, 1.4).p_f, py: 'ch15.laval_critical_pressures(2.0, 1.4)["p_f"]', rtol: 1e-9}` (0.09393) · `{name: 'shock M1', js:
  shockAt(0.7, 2, 1.4).M1, py: 'ch15.laval_shock_location(0.7, 2.0, 1.4)["M1"]', rtol: 1e-7}` (1.86271) · `{name: 'shock area', js:
  shockAt(0.7, 2, 1.4).As, py: 'ch15.laval_shock_location(0.7, 2.0, 1.4)["A_shock_over_At"]', rtol: 1e-7}` (1.51009) · `{name: 'exit
  M', js: exitState(0.7, 2, 1.4).M_exit, py: 'ch15.laval_exit_state(0.7, 2.0, 1.4)["M_exit"]', rtol: 1e-7}` (0.40669) · `{name: 'conv
  flux', js: convFlux(0.8, 1.4).flux, py: 'ch15.mass_flow_vs_back_pressure(0.8, 1.4)["mdot"]', rtol: 1e-10}` (0.56066) · invariants
  `{name: 'order', js: thresholds(2, 1.4).p_f < thresholds(2, 1.4).p_d && thresholds(2, 1.4).p_d < thresholds(2, 1.4).p_b ? 1 : 0,
  expect: 1}`, `{name: 'plateau', js: convFlux(0.3, 1.4).flux, expect: 0.684731, rtol: 1e-6}`.
- **Fit plan:** 360×640: status (1 line: regime) · `nozzle` (35 %) over `press` (65 %); `flow` hidden — "ṁ/ṁ_max 1.00 · shock at
  A/A_t 1.51" in the `press` title; `p0`, `T0`, `gamma` optional; the jet sketch is dropped below 400 px width. 844×345:
  `press` beside `flow`, `nozzle` as a strip. Desktop / 1000×700: rows [0.8, 1.2].

### E5 · duct_friction_heating
- **Title:** "Can friction speed a flow up?" · **Summary:** "Turn up friction or heat in a straight duct and watch every flow
  drift to Mach 1 — and choke." · **CORE:** C11, C12 (also N79, N80, N83–N88) · **Reference:**
  `amplitude_phase_second_order_II_3.html` (linked windows, a numbered derivation with live numbers) and ch11
  `kelvin_helmholtz_boundary` (a discriminant; two roots that merge).
- **meta:** `viz:order 5` · `viz:sections 15.8` · `viz:equations 15.45 15.46 15.47 15.48` · `viz:fluidpy ch15.duct_outlet_mach
  ch15.duct_parameter_A ch15.duct_state_ratios ch15.duct_entropy_change ch15.duct_choking_limits ch15.influence_coefficients
  ch15.rayleigh_temperature_sign ch15.fanno_line_Ts ch15.rayleigh_line_Ts` · `viz:derivations D21 D23`.
- **Physics (JS ↔ Python):** `ductA(M1, f, q, g)` ↔ `ch15.duct_parameter_A`; `outlet(M1, f, q, g, branch)` (closed form $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$;
  the sign that continues $M_2=M_1$ is + for $M_1<1$ and − for $M_1>1$; `null` when the discriminant is negative) ↔ `ch15.duct_outlet_mach(…, branch="continuous" | "shocked")` (NaN when choked); `ratios(M1, M2, f, q, g)` ↔ `ch15.duct_state_ratios` (keys `p_ratio`, `T_ratio`, `rho_ratio`, `u_ratio`, `p0_ratio`, `T0_ratio`); `dS(M1, M2, f, q, g)` = ΔS/C_v ↔ `ch15.duct_entropy_change`;
  `limits(M1, g)` ↔ `ch15.duct_choking_limits`; `coeffs(M, g, df, dq)` ↔ `ch15.influence_coefficients`; `tSign(M, g)` ↔
  `ch15.rayleigh_temperature_sign`; `fannoTs(M, g)`, `rayleighTs(M, g)` → [Δs/C_p, T/T*] ↔ `ch15.fanno_line_Ts`, `ch15.rayleigh_line_Ts` (keys `ds_over_cp`, `T_over_Tstar`) (ours: Fanno $T/T^*=(\gamma+1)/(2+(\gamma-1)M^2)$, $\Delta s/C_p=\tfrac{\gamma-1}{\gamma}\ln M-\tfrac{\gamma+1}{2\gamma}\ln\dfrac{2+(\gamma-1)M^2}{\gamma+1}$;
  Rayleigh $T/T^*=\left(\dfrac{(\gamma+1)M}{1+\gamma M^2}\right)^2$, $\Delta s/C_p=\ln M^2+\tfrac{\gamma+1}{\gamma}\ln\dfrac{\gamma+1}{1+\gamma M^2}$).
- **Views** (rows [0.7, 1.3]): 1. `duct` "The duct" (row 0): particles moving left to right, a rough-wall band whose
  hatching density follows f and a heater glow (rose; blue for cooling) whose strength follows q; gauges of M, p, T, $p_0$ at
  inlet and outlet (outlet / inlet ratios printed). 2. `map` "$M_2$ against $M_1$" (row 1, flex 1.3): axes 0…4; the curve for
  the current f or q bold (blue), neighbouring values faint, the no-change diagonal (muted), the shock branch thin, the
  choking locus $M_2=1$ dotted (amber), the forbidden region $M_1<1<M_2$ shaded muted, the cursor at ($M_1$, $M_2$). 3. `ts`
  "T–s diagram" (row 1, flex 1, `hidePortrait`): the process line — title "**Fanno line** / **Rayleigh line** (usual names —
  not in the book)" — with the inlet state (open dot), the outlet state (filled), the nose at M = 1 (amber), and in the
  heating case the temperature maximum at $M^2=1/\gamma$ (rose tick).
- **Controls:** `M1` "Inlet Mach number $M_1$" (0.1…4, step 0.01, default 0.5) · `f` "Friction $f=F/(p_1A)$" (0…0.8, step 0.005,
  default 0.1; help "0: smooth · 0.8: beyond choking for most inlets") · `q` "Heat $q=Q/h_1$" (−0.5…1.5, step 0.01, default 0;
  help "negative: cooling") · `proc` (modes: "friction only" · "heat only" · "both") · `root` toggle "outlet behind a shock"
  (optional) · `gamma` (optional).
- **Presets:** "subsonic + friction" {M1 0.5, f 0.1, q 0} · "supersonic + friction" {M1 2, f 0.3, q 0} · "heat at Mach 0.5" {M1
  0.5, f 0, q 0.2} · "paradox band: Mach 0.9" {M1 0.9, f 0, q 0.005} · "cooling" {M1 0.5, f 0, q −0.2} · "on the choking limit"
  {M1 0.5, f = f_max, q 0}.
- **Status:** friction, $M_1<1$: "friction: subsonic flow speeds up — 44 % of the way to choking" (= f/f_max, live) · friction,
  $M_1>1$: "friction: supersonic flow slows down" · heat in the band: "⚠️ paradox band: heat in, temperature down" · cooling:
  "cooling: the flow moves away from Mach 1" · discriminant < 0: "⛔ choked: no steady solution with this inlet".
- **Readouts:** "Duct number $A$" · "Discriminant" · "Outlet Mach $M_2$" · "$T_2/T_1$ out/in" · "$p_{02}/p_{01}$ out/in".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** · **terms** · **presets** · **status** ·
  **inspector**.
- **Terms:** title "Per unit friction (or heat), at the inlet"; signed bars `du` "$du/u$" blue · `dp` "$dp/p$" orange · `dT`
  "$dT/T$" rose · `dp0` "$dp_0/p_0$" amber — the pattern flips with the side of Mach 1; no total.
- **Inspector** (click a point of `map`): "$M_1$ = 0.50, f = 0.10: A = 0.25 × 1.05/(1 + 0.35 − 0.10)² = **0.168**; discriminant
  1 − 2 × 0.168 × 2.4 = **0.194**; roots $M_2$ = **0.589** (continuous) and 1.936 (needs ΔS < 0: forbidden)".
- **Explain:**
  0. *What the views show* — "**Top:** a straight duct; the hatching is wall friction, the glow is heat. **Left graph:** the
     outlet Mach number for every inlet Mach number at your f and q. **Right:** the same process as a line of states on a
     temperature–entropy diagram. Symbols: f and q are the book's lumped friction and heat. Most other texts use a friction
     length 4fL/D ('Fanno flow') and a stagnation-temperature ratio ('Rayleigh flow'); those names are not in the book."
  1. *The duct number* — "$A\equiv\dfrac{M_1^2\left(1+\frac12(\gamma-1)M_1^2+q\right)}{\left(1+\gamma M_1^2-f\right)^2}=\dfrac{M_2^2\left(1+\frac12(\gamma-1)M_2^2\right)}{\left(1+\gamma M_2^2\right)^2}$ (unnumbered, §15.8) = `…` × (1 + `…` + `…`)/(1 + `…` − `…`)² = **`…`** (a pure number, not an area)."
  2. *Is there a solution?* — "discriminant 1 − 2A(γ + 1) = 1 − 2 × `…` × 2.4 = **`…`**: `positive — two roots` / `zero — the
     roots have merged at M₂ = 1` / `negative — choked`."
  3. *The two roots* — "$M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46): upper sign $M_2$ = `…`, lower sign $M_2$ = `…`. The one that continues $M_2=M_1$ is the `upper
     (M₁ < 1)` / `lower (M₁ > 1)` sign: **`…`**."
  4. *The outlet state* — "momentum: $p_2/p_1$ = (1 + γM₁² − f)/(1 + γM₂²) = **`…`**; energy: $T_2/T_1$ = **`…`**; so
     $\rho_2/\rho_1$ = `…`, $p_{02}/p_{01}$ = **`…`**, ΔS/C_v = `…` > 0."
  5. *Which way does a little more push it?* — "$-\dfrac{dp}{p_1}=\dfrac{1+(\gamma-1)M^2}{1-M^2}\,df,\ \ \dfrac{du}{u}=-\dfrac{d\rho}{\rho}=\dfrac{p_1}{p}\dfrac{df}{1-M^2}$ (15.47, numerator sign corrected): du/u = df/(1 − `…`²) = `…` df; −dp/p = `…` df. $\dfrac{dT}{T_1}=\left(\dfrac{1-\gamma M^2}{1-M^2}\right)dq,\ \ \dfrac{u\,du}{h_1}=\dfrac{(\gamma-1)M^2}{1-M^2}\,dq$ (15.48, second relation's sign corrected): dT/T = `…` dq; u du/h
     = `…` dq." A two-line box "the book prints the first of (15.47) with $1-(\gamma-1)M^2$ and the second of (15.48) with
     $(1-\gamma)M^2$ — wrong sign: with your M they would give `…` and `…`."
  6. *How far from choking?* — "$f_{\max}=1+\gamma M_1^2-M_1\sqrt{2(\gamma+1)\left(1+\tfrac{\gamma-1}{2}M_1^2\right)}$ (ours) = **`…`**; $q_{\max}=\dfrac{(1+\gamma M_1^2)^2}{2(\gamma+1)M_1^2}-1-\dfrac{\gamma-1}{2}M_1^2$ (ours) = **`…`**; you are at `…` % of the limit."
  7. *The T–s view* (or the hint "turn the phone sideways; $T_2/T_1$ and the distance to choking are in the map's title").
  8. *Reading the current setting* — friction, subsonic: "Friction makes entropy. On the process line the only way to more
     entropy is toward the nose, M = 1 — so the flow speeds up, and its pressure and density drop." · friction, supersonic:
     "Same nose, approached from the other side: the flow slows, its pressure rises." · heat, $M^2<1/\gamma$: "Heat speeds the
     flow up and warms it." · paradox band: "You are between M = 0.845 and 1. Heat still speeds the flow up — so fast that
     the kinetic energy grows by more than the heat added, and the temperature falls." · cooling: "Removing heat lowers the
     entropy: the state moves away from the nose and never chokes." · choked: "A has passed 1/(2(γ + 1)) = 0.208. The duct
     cannot pass this inlet state; in practice the inlet Mach number drops (or a shock moves) until the outlet is exactly
     sonic."
- **Derivation tab:** **D21** (13 steps) `view: 'map'` on phones; goal `set` {preset subsonic + friction}; step 7 `live` "A =
  `…`"; step 11 `live` "1 − 2 × `…` × 2.4 = `…`"; step 12 `live` "$M_2^2$ = (`…` ± `…`)/(`…`)"; step 13 `set` {M1 2, f 0, q 0, root on},
  `watch` "with f = q = 0 the two dots are the diagonal and the normal shock". **D23** (8 steps) `view: 'map'` on phones (the
  temperature ratio in its title); goal `set` {preset heat at Mach 0.5}; step 5 `live` "dT/T = (1 − 1.4 × `…`²) du/u"; step 7 `set`
  {preset paradox band}, `watch` "the temperature gauge falls as q rises"; step 8 `live` "u du/h₁ = `…` dq". Interpret: `s => "At
  M1 = " + … + " heat " + (…raises/lowers…) + " the temperature: dT/T1 = " + … + " dq."`
- **Code:**
  ```python
  A = ch15.duct_parameter_A({{M1}}, {{f}}, {{q}}, {{g}})        # the duct number: {{A}}
  disc = 1 - 2*A*({{g}} + 1)                                    # {{disc}}  (< 0: choked)
  M2 = ch15.duct_outlet_mach({{M1}}, {{f}}, {{q}}, {{g}})       # (15.46), continuous root: {{M2}}
  r = ch15.duct_state_ratios({{M1}}, M2, {{f}}, {{q}}, {{g}})   # r['T_ratio'] {{T}}, r['p0_ratio'] {{p0}}
  lim = ch15.duct_choking_limits({{M1}}, {{g}})                 # f_max {{fmax}}, q_max {{qmax}}
  c = ch15.influence_coefficients({{M1}}, {{g}}, df=0.01)       # du/u {{du}}, dp/p {{dp}}
  ```
- **Walkthrough (8 steps):** 1. "A straight pipe" — "No area change. Can rough walls make the gas go faster?" `set` {M1 0.5, f
  0, q 0} · 2. "Add friction" — "Raise f. The outlet point leaves the diagonal — upward. The flow sped up." `set` {preset subsonic
  + friction}, `eq: 'A'` · 3. "One number, one quadratic" — "Inlet, friction and heat make one number A; the outlet Mach number
  solves a quadratic." `eq: 'M2'`, `derive: {id: 'D21', step: 12}` · 4. "Supersonic: the other way" — "Start at Mach 2 and friction
  slows the flow. Both sides head for Mach 1." `set` {preset supersonic + friction}, `terms: true` · 5. "Too much" — "Push f to its
  limit: the two roots merge at Mach 1. Beyond it, no solution — the duct is choked." `set` {preset on the choking limit}, `eq:
  'lim'` · 6. "Heat does the same" — "Swap friction for heat: same map, same drift to Mach 1." `set` {preset heat at Mach 0.5}, `eq:
  'dq'` · 7. "Heat that cools" — "At Mach 0.9 add a little heat and watch the temperature gauge." `set` {preset paradox band},
  `derive: {id: 'D23', step: 5}` · 8. "Your turn" — "Predict the friction that chokes a Mach 0.3 inlet. Then find it." `controls:
  ['M1', 'f', 'q']`.
- **Equations:** `laws` "The three laws with friction and heat" ref 'Eq. (15.45)' $\rho_1u_1=\rho_2u_2,\ \ p_1+\rho_1u_1^2=p_2+\rho_2u_2^2+p_1f,\ \ h_1+\tfrac12u_1^2+h_1q=h_2+\tfrac12u_2^2$, note "F is the retarding force here"
  · `A` "The duct number" ref 'unnumbered, §15.8' $A\equiv\dfrac{M_1^2\left(1+\frac12(\gamma-1)M_1^2+q\right)}{\left(1+\gamma M_1^2-f\right)^2}=\dfrac{M_2^2\left(1+\frac12(\gamma-1)M_2^2\right)}{\left(1+\gamma M_2^2\right)^2}$, live · `M2` "Outlet Mach number" ref 'Eq. (15.46)' $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$, live · `df`
  "A little friction (corrected)" ref 'Eq. (15.47)' $-\dfrac{dp}{p_1}=\dfrac{1+(\gamma-1)M^2}{1-M^2}\,df,\ \ \dfrac{du}{u}=-\dfrac{d\rho}{\rho}=\dfrac{p_1}{p}\dfrac{df}{1-M^2}$, live, note "the book prints $1-(\gamma-1)M^2$ in the first numerator"
  · `dq` "A little heat (corrected)" ref 'Eq. (15.48)' $\dfrac{dT}{T_1}=\left(\dfrac{1-\gamma M^2}{1-M^2}\right)dq,\ \ \dfrac{u\,du}{h_1}=\dfrac{(\gamma-1)M^2}{1-M^2}\,dq$, live, note "the book prints $(1-\gamma)M^2$ in the second" · `lim`
  "Choking limits" ref 'ours' $f_{\max}=1+\gamma M_1^2-M_1\sqrt{2(\gamma+1)\left(1+\tfrac{\gamma-1}{2}M_1^2\right)}$, $q_{\max}=\dfrac{(1+\gamma M_1^2)^2}{2(\gamma+1)M_1^2}-1-\dfrac{\gamma-1}{2}M_1^2$, live.
- **Check yourself:** (1) "Mach 0.5 inlet, f = 0.1: is the outlet faster or slower?" — "Faster: $M_2$ = 0.589; $-\dfrac{dp}{p_1}=\dfrac{1+(\gamma-1)M^2}{1-M^2}\,df,\ \ \dfrac{du}{u}=-\dfrac{d\rho}{\rho}=\dfrac{p_1}{p}\dfrac{df}{1-M^2}$ (15.47, numerator sign corrected) has du/u =
  df/(1 − M²) > 0." `set {M1: 0.5, f: 0.1, q: 0}` · (2) "What friction chokes a Mach 0.5 inlet?" — "f = 0.2275: $f_{\max}=1+\gamma M_1^2-M_1\sqrt{2(\gamma+1)\left(1+\tfrac{\gamma-1}{2}M_1^2\right)}$ (ours)." `set {M1:
  0.5, f: 0.2275}` · (3) "At which Mach numbers does adding heat lower the temperature?" — "0.845 < M < 1: the first of $\dfrac{dT}{T_1}=\left(\dfrac{1-\gamma M^2}{1-M^2}\right)dq,\ \ \dfrac{u\,du}{h_1}=\dfrac{(\gamma-1)M^2}{1-M^2}\,dq$ (15.48, second relation's sign corrected)
  changes sign at $M^2=1/\gamma$." `set {M1: 0.9, q: 0.005, f: 0}` · (4) "Set f = q = 0 and switch on the shock root at Mach 2. What
  is $M_2$?" — "0.577: the other root of $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46) is the normal shock $M_2^2=\dfrac{(\gamma-1)M_1^2+2}{2\gamma M_1^2+1-\gamma}$ (15.40)." `set {M1: 2, f: 0, q: 0, root: true}`.
- **Selftest parity rows:** `{name: 'A', js: ductA(0.5, 0.1, 0, 1.4), py: 'ch15.duct_parameter_A(0.5, 0.1, 0.0, 1.4)', rtol: 1e-12}`
  (0.168) · `{name: 'M2 sub', js: outlet(0.5, 0.1, 0, 1.4, 'continuous'), py: 'ch15.duct_outlet_mach(0.5, 0.1, 0.0, 1.4)', rtol:
  1e-10}` (0.58867) · `{name: 'M2 sup', js: outlet(2, 0.1, 0, 1.4, 'continuous'), py: 'ch15.duct_outlet_mach(2.0, 0.1, 0.0, 1.4)',
  rtol: 1e-10}` (1.88215) · `{name: 'shock root', js: outlet(2, 0, 0, 1.4, 'shock'), py: 'ch15.duct_outlet_mach(2.0, 0.0, 0.0, 1.4, branch="shocked")', rtol: 1e-10}` (0.57735) · `{name: 'f_max', js: limits(0.5, 1.4).f_max, py: 'ch15.duct_choking_limits(0.5,
  1.4)["f_max"]', rtol: 1e-12}` (0.22750) · `{name: 'q_max', js: limits(0.5, 1.4).q_max, py: 'ch15.duct_choking_limits(0.5,
  1.4)["q_max"]', rtol: 1e-12}` (0.46875) · `{name: 'T2/T1 heat', js: ratios(0.5, outlet(0.5, 0, 0.2, 1.4, 'continuous'), 0, 0.2,
  1.4).T_ratio, py: 'ch15.duct_state_ratios(0.5, ch15.duct_outlet_mach(0.5, 0.0, 0.2, 1.4), 0.0, 0.2, 1.4)["T_ratio"]', rtol: 1e-9}` (1.16504) · `{name: 'dT/T paradox', js: coeffs(0.9, 1.4, 0, 0.01).dT, py: 'ch15.influence_coefficients(0.9, 1.4,
  dq=0.01)["dT_over_T"]', rtol: 1e-10}` (−0.0070526) · `{name: 'rayleigh T/T*', js: rayleighTs(0.5, 1.4)[1], py: 'ch15.rayleigh_line_Ts(0.5, 1.4)["T_over_Tstar"]', rtol: 1e-12}` (0.79012) · invariant `{name: 'choke at limit', js: outlet(0.5, limits(0.5,
  1.4).f_max, 0, 1.4, 'continuous'), expect: 1, rtol: 1e-6}`.
- **Fit plan:** 360×640: status (1 line) · `duct` (30 %) over `map` (70 %); `ts` hidden — "T₂/T₁ 0.982 · 44 % to choking" in the
  `map` title; `root`, `gamma` optional. 844×345: `map` beside `ts`, `duct` as a strip. Desktop / 1000×700: rows [0.7, 1.3].

### E6 · mach_cone
- **Title:** "Why do you hear it only after it has passed?" · **Summary:** "Play time: every circle of sound grows from where
  the source was, and above Mach 1 their edge is a cone." · **CORE:** C13 (also N19, N20, N91, N92) · **Reference:**
  `angular_frequency_explorer_1.html` (system animation + graph on one clock, presets, end-of-run summary) and ch03
  `galilean_frames_cylinder` (one event, two frames).
- **meta:** `viz:order 6` · `viz:sections 15.9 15.2` · `viz:equations 15.17 15.49` · `viz:fluidpy ch15.source_wavefronts
  ch15.mach_angle ch15.mach_from_angle ch15.mach_cone_delay ch15.in_zone_of_action ch15.wave_speeds` · `viz:derivations D25`.
- **Physics (JS ↔ Python):** lengths in units of c × 1 s, times in s (readouts convert with the chosen c). `fronts(M, t, n, dt)` → [{xc, r, tEmit}] (front k emitted at t − (k + 1)dt from x = M·c·tEmit; radius c(t − tEmit)) ↔ `ch15.source_wavefronts(U, c, t, n, dt)` (keys `x_centres`, `radii`, `t_emit`, `x_source`); `machAngle(M)` (NaN below 1) ↔ `ch15.mach_angle`; `machFromAngle(mu)` ↔ `ch15.mach_from_angle`; `delay(h, M, c)` ↔
  `ch15.mach_cone_delay`; `inside(x, y, xs, M)` ↔ `ch15.in_zone_of_action(xs - x, y, M)` (its x is measured downstream of the source); `speeds(U, c)` ↔ `ch15.wave_speeds`.
- **Views** (rows [1.3, 0.7]): 1. `sky` "Source and wave fronts" (row 0, `equal`): the source (dot) moving right, its
  circular fronts (blue, older = fainter), the centre of each as a small tick on the flight path, for M > 1 the envelope (two
  amber lines) and the angle μ drawn at the source; a draggable observer (person icon) below the path with a marker "silence"
  (muted) / "boom" (amber flash) / "heard" (blue); pointer: drag the observer; click a circle to select it. 2. `angle` "Mach
  angle against M" (row 1, flex 1, `hidePortrait`): μ in degrees for M from 1 to 4 with the cursor; inset: the right triangle
  with sides $c\,\Delta t$ (opposite) and $U\,\Delta t$ (hypotenuse) drawn to scale for the current M. 3. `arrivals` "When
  each front reaches the observer" (row 1, flex 1, `hidePortrait`): emission time against arrival time for the fronts; for
  M > 1 the curve has a minimum — the boom is the earliest arrival.
- **Controls:** `M` "Source Mach number" (0…3, step 0.05, default 2; help "0: at rest · 1: fronts pile up on the nose · 3: a
  slender cone") · `t` (transport) · `n` "Number of fronts" (4…16, default 8) · `c` "Sound speed" (295…345 m/s, default 340;
  optional) · observer position (drag; `ox`, `oy` stored in state).
- **Transport:** `t` 0…8 s, rate 1/s, `end: 'hold'`; end-of-run card: "Every circle is old news. The cone is where the news
  has got to."
- **Presets:** "at rest" {M 0} · "M = 0.5" {M 0.5} · "M = 1" {M 1} · "M = 2" {M 2} · "M = 3" {M 3}.
- **Status:** M < 1: "subsonic: the fronts nest — the observer hears it coming" · M = 1: "sonic: all fronts touch at the nose"
  · M > 1, observer outside the cone: "zone of silence: overhead `…` s ago; the boom arrives in `…` s" · M > 1, inside: "zone
  of action: the boom passed `…` s ago".
- **Readouts:** "Mach angle $\mu$" · "$U+c$, $U-c$" · "Boom delay" · "Observer" (silence / heard).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** with an end-of-run card · **presets** ·
  **status** · **inspector**.
- **Inspector** (click a circle): "emitted at $t_e$ = 3.0 s, when the source was at x = U $t_e$ = 2 × 340 × 3.0 = **2040 m**; now
  (t = 5.0 s) its radius is c(t − $t_e$) = 340 × 2.0 = **680 m**; the source has moved on by U(t − $t_e$) = **1360 m** — twice
  as far".
- **Explain:**
  0. *What the views show* — "**Top:** a source moving right. Each <b class=c-blue>circle</b> is the sound it made at one
     instant, centred where it was *then*. The <b class=c-amber>amber lines</b> touch all the circles: the Mach cone. (μ is
     the Mach angle here — elsewhere in the chapter the same letter is a viscosity.)"
  1. *One circle* — the inspector's lines (live for the selected circle).
  2. *The triangle* — "In a time Δt the source travels U Δt and its sound c Δt. The tangent from the source to the circle
     makes a right angle with the radius, so $\sin\mu=1/M$ (15.49): sin μ = 1/`…` = `…`, μ = **`…`°** (= `…` rad)."
  3. *Upstream and downstream* — "The speeds of $p'(x,t)=f(x-(c+U)t)+g(x+(c-U)t)$ (15.17): U + c = `…`, U − c = `…` m/s — `one is negative: sound still gets
     ahead` / `both positive: nothing gets ahead`."
  4. *The observer* — "height h = `…` m below the path. Inside the cone if the source is more than h/tan μ = `…` m past
     overhead. $t_{\text{delay}}=\dfrac{h\sqrt{M^2-1}}{Mc}$ (ours) = `…` × √(`…`² − 1)/(`…` × `…`) = **`…` s** after it passed overhead."
  5. *The angle and arrival views* (or the hint "turn the phone sideways; μ in degrees is in the title of the main view").
  6. *Now* — "t = `Viz.live('t')` s: source at `Viz.live('xs')` m; the observer is `Viz.live('zone')`."
  7. *Reading the current setting* — subsonic: "The source never catches its own sound. Circles are crowded ahead (higher
     pitch) and spread behind — the Doppler effect — but every point is reached." · sonic: "The source rides on the front of
     all its circles: they pile up into a wall just ahead of it." · supersonic: "The source is ahead of everything it ever
     emitted. The circles overlap along two lines; outside them nothing has been heard yet. An observer on the ground hears
     nothing as it passes overhead, then everything at once: the boom."
- **Derivation tab:** **D25** (4 steps) `view: 'sky'`; goal `set` {preset M = 2, t 5}; step 1 `watch` "one circle and the source's
  path since it was emitted"; step 2 `watch` "the tangent from the source touches the circle at a right angle"; step 3 `live`
  "sin μ = c Δt/(U Δt) = 1/`…` = `…`"; step 4 `set` {M 1}, `watch` "at Mach 1 the angle opens to 90°: a plane front". Interpret: `s
  => "At M = " + … + " the cone's half-angle is " + … + "°."`
- **Code:**
  ```python
  W = ch15.source_wavefronts({{U}}, {{c}}, {{t}}, {{n}}, {{dt}})   # centres W['x_centres'], radii W['radii']
  mu = ch15.mach_angle({{M}})                                  # (15.49): {{mu}} rad = {{mudeg}} deg
  ch15.in_zone_of_action({{xs}} - {{ox}}, {{oy}}, {{M}})      # observer inside the cone? {{inside}}
  ch15.mach_cone_delay({{h}}, {{M}}, {{c}})                    # boom {{delay}} s after overhead
  ch15.wave_speeds({{U}}, {{c}})                               # ({{cp}}, {{cm}}) m/s
  ```
- **Walkthrough (6 steps):** 1. "A source that moves" — "It sends out sound all the time. Where has the sound got to?" `set`
  {preset M = 0.5, t 0}, `play: true` · 2. "Old news" — "Click a circle: it is centred where the source was when that sound
  left." `inspect: true` · 3. "Catching up" — "At Mach 1 the source keeps pace with its own sound. The circles pile up on its
  nose." `set` {preset M = 1}, `play: true` · 4. "Outrunning it" — "At Mach 2 the circles cross. Their edge is a cone the source
  drags behind." `set` {preset M = 2}, `play: true`, `eq: 'mach'` · 5. "The angle" — "A right triangle: sound's distance over the
  source's distance." `derive: {id: 'D25', step: 3}` · 6. "Your turn" — "Drag the observer under the flight path. Predict the
  delay of the boom at Mach 3. Then play." `controls: ['M', 'n']`.
- **Equations:** `front` "A wave front" ref 'unnumbered, §15.9' $\lvert\mathbf x-\mathbf x_s(t_e)\rvert=c\,(t-t_e)$, live for
  the selected circle · `mach` "Mach angle" ref 'Eq. (15.49)' $\sin\mu=1/M$, live · `delay` "Boom delay" ref 'ours' $t_{\text{delay}}=\dfrac{h\sqrt{M^2-1}}{Mc}$, live ·
  `speeds` "Signal speeds" ref 'Eq. (15.17)' $p'(x,t)=f(x-(c+U)t)+g(x+(c-U)t)$, live with U ± c.
- **Check yourself:** (1) "What is the Mach angle at Mach 2?" — "30°: $\sin\mu=1/M$ (15.49) gives sin μ = ½." `set {M: 2}` · (2) "Faster source:
  wider or narrower cone?" — "Narrower: sin μ = 1/M falls." `set {M: 3}` · (3) "At Mach 0.5 is there a zone of silence?" — "No:
  U − c < 0, the sound gets ahead of the source." `set {M: 0.5}` · (4) "An aircraft at Mach 1.5 passes 10 km overhead (c = 295
  m/s). When do you hear it?" — "About 25 s later: $t_{\text{delay}}=\dfrac{h\sqrt{M^2-1}}{Mc}$ (ours)."
- **Selftest parity rows:** `{name: 'mu(2)', js: machAngle(2), py: 'ch15.mach_angle(2.0)', rtol: 1e-12}` (0.5235988) · `{name:
  'mu(1.2)', js: machAngle(1.2), py: 'ch15.mach_angle(1.2)', rtol: 1e-12}` (0.985111) · `{name: 'M from mu', js: machFromAngle(
  Math.PI/6), py: 'ch15.mach_from_angle(np.pi/6)', rtol: 1e-12}` (2.0) · `{name: 'delay', js: delay(1e4, 1.5, 295), py:
  'ch15.mach_cone_delay(1.0e4, 1.5, 295.0)', rtol: 1e-12}` (25.2663) · `{name: 'radius', js: fronts(2, 3, 4, 0.5)[1].r, py: 'ch15.source_wavefronts(2.0, 1.0, 3.0, 4, 0.5)["radii"][1]', rtol: 1e-12}` (1.0) · `{name: 'U-c', js: speeds(680, 340)[1], py:
  'ch15.wave_speeds(680.0, 340.0)[1]', rtol: 1e-12}` (340) · invariant `{name: 'tangent', js: Math.asin(fronts(2, 3, 4, 0.5)[1].r/(2*3
  - fronts(2, 3, 4, 0.5)[1].xc)), expect: Math.PI/6, rtol: 1e-9}`.
- **Fit plan:** 360×640: status (2 short lines) · `sky` alone (the two lower views hidden) — "μ = 30.0° · boom delay 2.5 s" in
  its title; `c` optional. 844×345: `sky` beside `angle`. Desktop / 1000×700: rows [1.3, 0.7].


### E7 · oblique_shock_wedge
- **Title:** "One wedge, how many shocks?" · **Summary:** "Open a wedge in a supersonic stream: two possible shocks run
  together, merge and vanish." · **CORE:** C14 (also N93–N99, N102) · **Reference:** `amplitude_phase_second_order_II_3.html`
  (linked windows with a crosshair readout, one section per window) and ch11 `kelvin_helmholtz_boundary` (two roots that run
  together and merge).
- **meta:** `viz:order 7` · `viz:sections 15.9` · `viz:equations 15.50 15.51 15.52 15.53 15.55` · `viz:fluidpy ch15.oblique_shock
  ch15.oblique_shock_components ch15.deflection_angle ch15.oblique_shock_curve ch15.max_deflection ch15.sonic_deflection
  ch15.weak_wave_pressure ch15.weak_wave_error ch15.mach_angle ch15.wedge_flow` · `viz:derivations D26 D27`.
- **Physics (JS ↔ Python):** `deflection(M1, sigma, g)` ↔ `ch15.deflection_angle`; `maxDeflection(M1, g)` (closed form: the
  quadratic in sin²σ) → {delta_max, sigma_max} ↔ `ch15.max_deflection`; `sonicPoint(M1, g)` → {delta_sonic, sigma_sonic} ↔ `ch15.sonic_deflection`; `oblique(M1,
  delta, branch, g)` (`brentq` on (μ, σ_max) or (σ_max, π/2)) → {sigma, Mn1, Mn2, M2, p, rho, T, p0, attached} ↔ `ch15.oblique_shock(M1, delta=…, branch=…, gamma=g)` (keys `p_ratio`, `rho_ratio`, `T_ratio`, `p0_ratio`; Python raises `ValueError("detached")` where the JS returns `attached: false`); `components(M1, sigma, g)` ↔ `ch15.oblique_shock_components`; `linear(M1, delta, g)` ↔ `ch15.weak_wave_pressure`; `linErr(M1, delta, g)` ↔ `ch15.weak_wave_error(…)["rel_error"]`; `machAngle(M)` ↔ `ch15.mach_angle`.
- **Views** (rows [1.1, 0.9]): 1. `flow` "The wedge in the stream" (row 0, flex 1.4, `equal`): the wedge (half-angle δ) with
  its attached shock (rose line at σ), streamlines turning through δ at the shock (blue upstream, lighter downstream), the
  Mach wave at μ as a muted ghost; at the clicked point the two velocity triangles: normal component (blue), tangential
  component (amber — the same length on both sides). Beyond $\delta_{\max}$: a labelled **sketch** of a detached bow shock
  ("sketch, not a computation"). Pointer: click a streamline point to place the triangles. 2. `chart` "δ against σ" (row 1,
  flex 1.2): σ from 0 to 90° (label "σ (= β elsewhere)"), δ from 0 to 50° ("δ (= θ elsewhere)"); the curve for the current $M_1$
  bold with its weak part teal and strong part rose, the maximum-deflection point (amber) and the sonic point ($M_2$ = 1,
  muted ring), the two ends at σ = μ and 90°, faint curves for $M_1$ = 1.5, 2, 3, 5, the horizontal line at the wedge's δ and the
  two solution dots. 3. `bars` "Pressure rise: exact against the linear rule" (row 1, flex 0.7, `hidePortrait`): two bars
  $(p_2-p_1)/p_1$ — exact (orange) and $\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta$ (muted) — with the gap in per cent.
- **Controls:** `M1` "Upstream Mach number $M_1$" (1.1…5, step 0.01, default 2) · `deltaDeg` "Wedge half-angle δ" (0…45°, step
  0.1, default 10; help "0: a Mach wave · 45°: detached at any Mach number below ≈ 9") · `branch` (chips: weak · strong) · `gamma`
  (optional) · `showLinear` toggle "Show the linear rule" (optional).
- **Presets:** "thin wedge: nearly a Mach wave" {M1 2, deltaDeg 1} · "10° at Mach 2" {M1 2, deltaDeg 10} · "just below detachment"
  {M1 2, deltaDeg 22.9} · "detached" {M1 2, deltaDeg 25} · "strong root: nearly a normal shock" {M1 2, deltaDeg 10, branch strong}.
- **Status:** attached, $M_2>1$: "attached · weak root · supersonic behind ($M_2$ = 1.64)" · attached, weak, $M_2<1$: "attached,
  but subsonic behind — the thin band just below the maximum" · strong: "strong root: subsonic behind ($M_2$ = 0.60)" · δ >
  $\delta_{\max}$: "⛔ detached: δ exceeds 22.97°".
- **Readouts:** "Shock angle σ (= β)" · "Mach angle μ" · "Normal Mach $M_1\sin\sigma$" · "$M_2$ downstream" · "$p_2/p_1$ down/up".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **inspector** · **presets** · **status** · **terms**
  (exact against linear).
- **Inspector** (click a point of `flow`): "upstream: V₁ = 2.000 c₁; normal part $V_1\sin\sigma$ = 2 × sin 39.31° = **1.267**,
  along the shock $V_1\cos\sigma$ = **1.547** · downstream: normal part 1.267/1.458 = **0.869** c₁ (divided by the density
  ratio), along the shock **1.547** (unchanged) · new direction: atan(0.869/1.547) = 29.31° to the shock = σ − δ, so δ =
  **10.0°** ✓".
- **Explain:**
  0. *What the views show* — "**Top:** a wedge in a supersonic stream from the left. The <b class=c-rose>shock</b> stands at
     the angle σ to the stream; the flow behind it is turned through δ, parallel to the wedge. The book's σ and δ are called β
     and θ in most other books. **Left graph:** for one Mach number, every shock angle and the turn it gives."
  1. *The limits* — "$\sin\mu=1/M$ (15.49): μ = asin(1/`…`) = **`…`°** — the flattest possible wave. Largest turn (ours, from
     $d\delta/d\sigma=0$): $\delta_{\max}$ = **`…`°** at σ = `…`°."
  2. *The two shocks for your wedge* — "Solving $\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}$ (15.53) for δ = `…`°: σ = **`…`°** (weak) and **`…`°** (strong)."
  3. *The normal Mach number* — "$M_{n1}=M_1\sin\sigma,\ \ M_{n2}=M_2\sin(\sigma-\delta)$ (unnumbered, §15.9): $M_{n1}$ = `…` × sin `…`° = **`…`** (> 1, as a normal shock needs)."
  4. *The jumps* — "$\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2\sin^2\sigma-1\right)$ (15.50) = 1 + 1.1667 × (`…`² − 1) = **`…`**; $\dfrac{\rho_2}{\rho_1}=\dfrac{(\gamma+1)M_1^2\sin^2\sigma}{(\gamma-1)M_1^2\sin^2\sigma+2}=\dfrac{\tan\sigma}{\tan(\sigma-\delta)}$ (15.51) = **`…`** (and tan σ/tan(σ − δ) = `…`/`…` = `…` ✓)."
  5. *Downstream Mach number* — "$M_2^2\sin^2(\sigma-\delta)=\dfrac{(\gamma-1)M_1^2\sin^2\sigma+2}{2\gamma M_1^2\sin^2\sigma+1-\gamma}$ (15.52): $M_{n2}$ = `…`; divide by sin(σ − δ) = sin `…`° = `…`: $M_2$ = **`…`**."
  6. *Against the linear rule* — "$\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta$ (15.55) = 1.4 × `…`²/√(`…`² − 1) × `…` rad = `…`; exact `…`: the rule is `…` % low."
  7. *The bars* (or the hint "turn the phone sideways; the linear-rule gap is in the chart's title").
  8. *Reading the current setting* — thin wedge: "σ is barely above μ: the shock is almost a Mach wave, the jump almost
     isentropic, and the linear rule $\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta$ (15.55) is within a few per cent." · weak, supersonic behind: "The usual case for a wedge.
     'Weak' names the flatter of the two roots — at 20° it would still nearly triple the pressure." · near the maximum: "The
     two roots are about to meet. Just below the maximum the 'weak' shock already leaves subsonic flow behind." · detached:
     "No straight shock can turn the stream this much. The shock stands off the nose as a curved bow: normal on the axis,
     strong near it, weak far away. (The picture is a sketch.)" · strong root: "The steeper root: nearly a normal shock,
     subsonic behind. It needs a high pressure downstream to hold it; a free wedge does not choose it."
- **Derivation tab:** **D26** (8 steps) `view: 'flow'` on phones; goal `set` {preset 10° at Mach 2}; step 1 `watch` "the two
  triangles at the clicked point"; step 3 `live` "v = `…` c₁ on both sides"; step 4 `live` "$M_{n1}$ = `…` × sin `…`° = `…`"; step 8
  `live` "$M_2$ = `…`/sin `…`° = `…`". **D27** (8 steps) `view: 'chart'` on phones; goal `set` {preset 10° at Mach 2}; step 5 `live`
  "tan δ = `…`"; step 8 `set` {deltaDeg 0}, `watch` "the curve meets zero at σ = μ and at 90°". Interpret: `s => "At M1 = " + … + ", a
  " + … + "° wedge has shocks at " + … + "° and " + … + "°."`
- **Code:**
  ```python
  w = ch15.oblique_shock({{M1}}, delta=np.deg2rad({{d}}), branch="{{branch}}", gamma={{g}})
  np.degrees(w["sigma"]), w["Mn1"], w["M2"]                    # {{sig}} deg, {{Mn1}}, {{M2}}
  w["p_ratio"], w["p0_ratio"]                                # (15.50): {{p}}; {{p0}}
  dm = ch15.max_deflection({{M1}}, {{g}})                      # dm['delta_max'] {{dmax}} deg, dm['sigma_max'] {{sdm}} deg
  lin = ch15.weak_wave_pressure({{M1}}, np.deg2rad({{d}}), {{g}})   # (15.55): {{lin}}
  ch15.weak_wave_error({{M1}}, np.deg2rad({{d}}), {{g}})["rel_error"]   # {{err}}
  ```
- **Walkthrough (7 steps):** 1. "A wedge in the way" — "The stream cannot be warned. Where does the front stand?" `set` {preset
  10° at Mach 2} · 2. "Slide along it" — "Click a streamline. Along the shock the speed is unchanged; only the part across it
  jumps." `inspect: true`, `derive: {id: 'D26', step: 3}` · 3. "A normal shock in disguise" — "Replace the Mach number by its
  normal part and every formula of the normal shock applies." `eq: 'p'` · 4. "Two answers" — "For 10° the chart shows two dots:
  a flat shock and a steep one." `eq: 'tan'`, `set` {branch strong} · 5. "They run together" — "Open the wedge. The dots approach,
  merge at 22.97° — and are gone." `set` {preset just below detachment}, `derive: {id: 'D27', step: 7}` · 6. "Small turns are easy"
  — "For a 1° wedge the pressure rise is a straight line in δ." `set` {preset thin wedge}, `terms: true`, `eq: 'lin'` · 7. "Your
  turn" — "Predict the largest wedge that keeps an attached shock at Mach 3. Then find it." `controls: ['M1', 'deltaDeg',
  'branch']`.
- **Equations:** `Mn` "Normal Mach numbers" ref 'unnumbered, §15.9' $M_{n1}=M_1\sin\sigma,\ \ M_{n2}=M_2\sin(\sigma-\delta)$, live · `p` "Pressure jump" ref 'Eq. (15.50)' $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2\sin^2\sigma-1\right)$,
  live · `rho` "Density jump and the two tangents" ref 'Eq. (15.51)' $\dfrac{\rho_2}{\rho_1}=\dfrac{(\gamma+1)M_1^2\sin^2\sigma}{(\gamma-1)M_1^2\sin^2\sigma+2}=\dfrac{\tan\sigma}{\tan(\sigma-\delta)}$, live · `M2` "Downstream Mach number" ref 'Eq.
  (15.52)' $M_2^2\sin^2(\sigma-\delta)=\dfrac{(\gamma-1)M_1^2\sin^2\sigma+2}{2\gamma M_1^2\sin^2\sigma+1-\gamma}$, live · `tan` "Deflection against shock angle" ref 'Eq. (15.53)' $\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}$, live · `lin` "Weak-wave rule" ref 'Eq.
  (15.55)' $\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta$, live with the gap.
- **Check yourself:** (1) "Mach 2, 10° wedge: what is the weak shock angle?" — "39.3°: the smaller root of $\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}$ (15.53)." `set {M1: 2,
  deltaDeg: 10}` · (2) "What is the largest turn an attached shock can give at Mach 2?" — "22.97°: the top of the δ–σ curve."
  `set {M1: 2, deltaDeg: 22.9}` · (3) "Is the flow behind the weak 10° shock subsonic?" — "No: $M_2$ = 1.64. Only the normal
  component (0.80) is subsonic; $M_2^2\sin^2(\sigma-\delta)=\dfrac{(\gamma-1)M_1^2\sin^2\sigma+2}{2\gamma M_1^2\sin^2\sigma+1-\gamma}$ (15.52) has the division by sin(σ − δ)." · (4) "At 5°, is the linear rule too high or too low?" —
  "About 10 % too low for a compression: $\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta$ (15.55) misses the second-order rise." `set {deltaDeg: 5, showLinear: true}`.
- **Selftest parity rows:** `{name: 'sigma weak', js: oblique(2, rad(10), 'weak', 1.4).sigma, py: 'ch15.oblique_shock(2.0,
  delta=np.deg2rad(10.0))["sigma"]', rtol: 1e-9}` (0.686158 rad = 39.314°) · `{name: 'M2 weak', js: oblique(2, rad(10), 'weak',
  1.4).M2, py: 'ch15.oblique_shock(2.0, delta=np.deg2rad(10.0))["M2"]', rtol: 1e-8}` (1.64052) · `{name: 'p weak', js: oblique(2,
  rad(10), 'weak', 1.4).p, py: 'ch15.oblique_shock(2.0, delta=np.deg2rad(10.0))["p_ratio"]', rtol: 1e-8}` (1.70658) · `{name: 'sigma
  strong', js: oblique(2, rad(10), 'strong', 1.4).sigma, py: 'ch15.oblique_shock(2.0, delta=np.deg2rad(10.0),
  branch="strong")["sigma"]', rtol: 1e-9}` (1.46084 rad = 83.700°) · `{name: 'delta max', js: maxDeflection(2, 1.4).delta_max, py:
  'ch15.max_deflection(2.0, 1.4)["delta_max"]', rtol: 1e-10}` (0.400964 rad = 22.974°) · `{name: 'delta(45°)', js: deflection(2,
  Math.PI/4, 1.4), py: 'ch15.deflection_angle(2.0, np.pi/4, 1.4)', rtol: 1e-12}` (0.257324 rad = 14.744°) · `{name: 'linear', js:
  linear(2, rad(10), 1.4), py: 'ch15.weak_wave_pressure(2.0, np.deg2rad(10.0), 1.4)', rtol: 1e-12}` (0.564293) · invariants `{name:
  'delta at mu', js: deflection(2, machAngle(2), 1.4), expect: 0, atol: 1e-12}`, `{name: 'tangents', js: Math.tan(sW)/Math.tan(sW -
  rad(10)), expect: 1.458422, rtol: 1e-6}` (sW = the weak σ).
- **Fit plan:** 360×640: status (1 line) · `flow` (50 %) over `chart` (50 %); `bars` hidden — "linear rule 20 % low" in the `chart`
  title; faint curves for other $M_1$ dropped below 420 px. 844×345: `flow` beside `chart`. Desktop / 1000×700: rows [1.1, 0.9].

### E8 · prandtl_meyer_fan
- **Title:** "How does a supersonic stream turn a corner?" · **Summary:** "Build a turn from 1, 3, 7, 21, then infinitely many
  Mach waves and watch ν(M) appear." · **CORE:** C15 (also N38, N103–N107) · **Reference:** `angular_frequency_explorer_1.html`
  (live working, presets, regime notes) and ch14 `downwash_induced_drag` ("build the sum" 1 → 3 → 7 → 21 → ∞ with term bars).
- **meta:** `viz:order 8` · `viz:sections 15.9` · `viz:equations 15.29 15.49 15.56` · `viz:fluidpy ch15.prandtl_meyer
  ch15.mach_from_prandtl_meyer ch15.prandtl_meyer_turn ch15.prandtl_meyer_max ch15.mach_angle ch15.expansion_fan
  ch15.fan_staircase ch15.curved_wall_mach_lines ch15.weak_wave_pressure` · `viz:derivations D29 D30`.
- **Physics (JS ↔ Python):** `nu(M, g)` ↔ `ch15.prandtl_meyer`; `nuMax(g)` ↔ `ch15.prandtl_meyer_max`; `machFromNu(v, g)` (`brentq`
  on [1, 1000]) ↔ `ch15.mach_from_prandtl_meyer`; `turn(M1, turn, g)` → {M2, nu1, nu2, p, T, mu1, mu2} ↔ `ch15.prandtl_meyer_turn` (keys `M2`, `nu1`, `nu2`, `mu1`, `mu2`, `p_ratio`, `T_ratio`, `delta_book`); `staircase(M1, turn, N, g)` (N explicit steps of $d\delta=-\dfrac{\sqrt{M^2-1}}{M}\dfrac{dM}{1+\frac12(\gamma-1)M^2}$, C.6) ↔ `ch15.fan_staircase`; `fan(M1, turn,
  g, nRays, y0)` → ray angles and the streamline ↔ `ch15.expansion_fan`; `machAngle(M)` ↔ `ch15.mach_angle`; `linear(M, delta, g)`
  ↔ `ch15.weak_wave_pressure`. A compression (turn < 0) uses the same `turn()` for the isentropic part and reports the first crossing of Mach lines as `ch15.curved_wall_mach_lines(…)["first_crossing"]` does for a short circular-arc wall (ours; labelled "estimate").
- **Views** (rows [1.1, 0.9]): 1. `corner` "The corner" (row 0, flex 1.4, `equal`): the wall turning by the chosen angle at
  the origin; for an expansion the fan between its first ray (at $\mu_1$ to the upstream wall) and last ray (at $\mu_2$ to
  the downstream wall), N rays when N is finite; one streamline at the chosen height curving through the fan, coloured by M;
  a small pressure trace along the wall (a step) and along the streamline (a smooth fall). For a compression: the Mach
  lines converging, their first crossing marked "a shock forms here (estimate)". Pointer: click a ray. 2. `nuview` "ν against
  M" (row 1, flex 1.2): ν in degrees for M from 1 to 5 (the limit `nuMax` dashed), cursors at $M_1$ and $M_2$, the turn as the
  vertical gap between them (amber bracket). 3. `stairs` "After N waves" (row 1, flex 1, `hidePortrait`): Mach number after
  each of N equal small turns (blue staircase) against the exact curve (muted), and two bars $p_2/p_1$: exact against the
  linear rule $\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta$.
- **Controls:** `M1` "Upstream Mach number $M_1$" (1…4, step 0.01, default 2) · `turnDeg` "Turn (+ away from the flow =
  expansion; book: δ = −turn)" (−15…40°, step 0.5, default 10) · `N` "Number of waves" (chips 1 · 3 · 7 · 21 · ∞) · `y0`
  "Streamline height" (0.3…2 wall units, default 1) · `gamma` (optional).
- **Presets:** "10° at Mach 2" {M1 2, turnDeg 10, N ∞} · "a sonic start" {M1 1, turnDeg 20} · "a large turn" {M1 2, turnDeg 40} ·
  "a gentle compression: Mach lines converge" {M1 2, turnDeg −8} · "one wave against infinitely many" {M1 2, turnDeg 10, N 1}.
- **Status:** expansion: "isentropic expansion: $M_2$ = 2.38, pressure down to 55 %" · turn beyond `nuMax − nu(M1)`: "⛔ beyond
  the largest turn: the gas would reach zero pressure" · compression: "compression: Mach lines cross `…` wall-heights out — a
  shock forms there (estimate)".
- **Readouts:** "$\nu(M_1)$" · "$\nu(M_2)$" · "$M_2$ exact" · "$M_2$ after N waves" · "$p_2/p_1$ down/up".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **terms** (the turn as N equal steps, each with its
  rise in M) · **presets** · **status** · **inspector**.
- **Inspector** (click a ray): "ray 4 of 7: flow turned 5.71° so far; M = 2.21 there (ν = 32.09°); Mach angle asin(1/2.21) =
  26.9°; the ray makes 26.9° − 5.71° = 21.2° with the upstream wall direction".
- **Explain:**
  0. *What the views show* — "**Top:** a wall that bends away from a supersonic stream. The <b class=c-blue>fan</b> is made
     of Mach waves from the corner; the streamline bends as it crosses them. **Left graph:** the Prandtl–Meyer function ν(M).
     (ν is a function of Mach number here — elsewhere the same letter is a viscosity.) **Signs:** our turn is positive for an
     expansion; the book's δ is positive for a compression: δ = −turn."
  1. *ν at your upstream Mach number* — "$\nu(M)=\sqrt{\dfrac{\gamma+1}{\gamma-1}}\tan^{-1}\sqrt{\dfrac{\gamma-1}{\gamma+1}\left(M^2-1\right)}-\tan^{-1}\sqrt{M^2-1}$ (15.56): √6 × atan(√((0.4/2.4)(`…`² − 1))) − atan(√(`…`² − 1)) = 2.449 × `…`° − `…`° =
     **`…`°**."
  2. *After the turn* — "$\delta_1+\nu(M_1)=\delta_2+\nu(M_2)$ (unnumbered, §15.9) with $\delta_1-\delta_2$ = +`…`° (ours: turn = `…`°): ν(M₂) = `…`° + `…`° = **`…`°**."
  3. *The new Mach number* — "invert ν (root-finder on [1, 1000]): $M_2$ = **`…`**."
  4. *The fan* — "$\sin\mu=1/M$ (15.49): first ray at $\mu_1$ = `…`° to the old wall, last ray at $\mu_2$ = `…`° to the new wall; the fan is
     `…`° wide."
  5. *The pressure* — "isentropic, so $p_0$ is the same on both sides: $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29) at each Mach number gives $p_2/p_1$ = `…`/`…` =
     **`…`**; the linear rule $\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta$ (15.55) would give 1 − `…` = `…`."
  6. *N waves instead of infinitely many* — "$d\delta=-\dfrac{\sqrt{M^2-1}}{M}\dfrac{dM}{1+\frac12(\gamma-1)M^2}$ (unnumbered, §15.9) applied N = `…` times in equal steps: $M_2$ ≈ `…` (exact `…`; gap `…` %). The
     gap falls like 1/N — the staircase is a Riemann sum of the integral $\nu(M)=\displaystyle\int_1^M\dfrac{\sqrt{M^2-1}}{1+\frac12(\gamma-1)M^2}\dfrac{dM}{M}$ (15.56, integral form)."
  7. *The staircase view* (or the hint "turn the phone sideways; $M_2$ exact and by N waves are in the title of the ν view").
  8. *Reading the current setting* — expansion: "Each wave turns the stream a little and speeds it up a little. Because the
     Mach angle shrinks as M grows, later waves lean further downstream: they fan out and never cross. Nothing is lost." ·
     sonic start: "From Mach 1 the first ray is at 90°. ν starts at zero: this is where its lower limit comes from." · large
     turn: "The pressure is falling toward zero. The most a flow can turn from Mach 1 is 130.45° (γ = 1.4) — the limit where
     all thermal energy has become speed." · compression: "Reverse the turn and each wave *slows* the stream, so later waves
     are steeper and catch the earlier ones. Where they first cross, the smooth solution ends and an oblique shock begins
     (C14)."
- **Derivation tab:** **D29** (9 steps) `view: 'corner'` on phones; goal `set` {preset one wave against infinitely many}; step 3
  `live` "dδ = −(dV/V) × √(`…`² − 1)"; step 9 `live` "dM = `…` per degree of turn at M = `…`", `set` {N 7}, `watch` "each stair is one
  application of this line". **D30** (10 steps) `view: 'nuview'` on phones; goal `set` {preset 10° at Mach 2}; step 1 `live` "t =
  √(`…`² − 1) = `…`"; step 7 `live` "2.449 × atan(0.408 × `…`) = `…`°"; step 9 `set` {preset a sonic start}, `watch` "the curve starts
  from zero at M = 1"; step 10 `live` "ν(`…`) = `…`°". Interpret: `s => "A turn of " + … + "° takes Mach " + … + " to " + … + "."`
- **Code:**
  ```python
  nu1 = ch15.prandtl_meyer({{M1}}, {{g}})                       # (15.56): {{nu1}} rad = {{nu1deg}} deg
  t = ch15.prandtl_meyer_turn({{M1}}, np.deg2rad({{turn}}), {{g}})   # + = expansion (book: delta = -turn)
  t["M2"], t["p_ratio"]                                       # {{M2}}, {{p}}
  np.degrees([t["mu1"], t["mu2"]])                             # first and last ray: {{mu1}}, {{mu2}} deg
  st = ch15.fan_staircase({{M1}}, np.deg2rad({{turn}}), {{N}}, {{g}})   # N small turns
  st["M"][-1]                                                  # {{MN}}  (exact {{M2}})
  ```
- **Walkthrough (7 steps):** 1. "A wall that falls away" — "The stream is supersonic and the wall bends away from it. How
  does it follow?" `set` {preset 10° at Mach 2} · 2. "One wave" — "Try a single Mach wave doing the whole turn. The stream
  speeds up — roughly." `set` {N 1}, `eq: 'dnu'` · 3. "Three, seven, twenty-one" — "Split the turn into more waves. The staircase
  settles." `set` {N 21}, `terms: true` · 4. "Infinitely many" — "The limit is a smooth fan, and the sum becomes an integral: one
  function ν(M)." `set` {N ∞}, `eq: 'nu'`, `derive: {id: 'D30', step: 10}` · 5. "A subtraction" — "ν after = ν before + the turn. Read
  it off the curve." `eq: 'turn'`, `code: {id: 'pm', lines: [1, 3]}` · 6. "The other way" — "Turn the wall into the stream: the
  waves converge and a shock must form." `set` {preset a gentle compression}, `notes: true` · 7. "Your turn" — "Predict the Mach
  number after a 20° turn from Mach 1. Then drag." `controls: ['M1', 'turnDeg', 'N']`.
- **Equations:** `dnu` "One small turn" ref 'unnumbered, §15.9' $d\delta=-\dfrac{\sqrt{M^2-1}}{M}\dfrac{dM}{1+\frac12(\gamma-1)M^2}$, live · `nui` "The sum of small turns" ref 'Eq. (15.56),
  integral form' $\nu(M)=\displaystyle\int_1^M\dfrac{\sqrt{M^2-1}}{1+\frac12(\gamma-1)M^2}\dfrac{dM}{M}$ · `nu` "Prandtl–Meyer function" ref 'Eq. (15.56)' $\nu(M)=\sqrt{\dfrac{\gamma+1}{\gamma-1}}\tan^{-1}\sqrt{\dfrac{\gamma-1}{\gamma+1}\left(M^2-1\right)}-\tan^{-1}\sqrt{M^2-1}$, live · `turn` "Turning rule" ref 'unnumbered,
  §15.9' $\delta_1+\nu(M_1)=\delta_2+\nu(M_2)$, live, note "δ > 0 when the wall turns into the flow; our turn = −δ" · `mu` "Mach angle" ref 'Eq. (15.49)'
  $\sin\mu=1/M$ · `p` "Pressure from the stagnation ratio" ref 'Eq. (15.29)' $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$.
- **Check yourself:** (1) "What is ν at Mach 2 (γ = 1.4)?" — "26.38°: $\nu(M)=\sqrt{\dfrac{\gamma+1}{\gamma-1}}\tan^{-1}\sqrt{\dfrac{\gamma-1}{\gamma+1}\left(M^2-1\right)}-\tan^{-1}\sqrt{M^2-1}$ (15.56)." `set {M1: 2}` · (2) "Mach 2 turned 10° away: the new
  Mach number?" — "2.385: ν₂ = 36.38°." `set {M1: 2, turnDeg: 10}` · (3) "Does an expansion fan lose stagnation pressure?" — "No:
  it is isentropic; that is why $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29) can be used across it." · (4) "With N = 1 is $M_2$ too high or too low?" — "Too low (2.363
  against 2.385): the rate of $d\delta=-\dfrac{\sqrt{M^2-1}}{M}\dfrac{dM}{1+\frac12(\gamma-1)M^2}$ (unnumbered, §15.9) grows along the turn and one step uses only its starting value." `set {N: 1}`.
- **Selftest parity rows:** `{name: 'nu(2)', js: nu(2, 1.4), py: 'ch15.prandtl_meyer(2.0, 1.4)', rtol: 1e-12}` (0.460414 rad =
  26.380°) · `{name: 'nu max', js: nuMax(1.4), py: 'ch15.prandtl_meyer_max(1.4)', rtol: 1e-12}` (2.276853 rad = 130.454°) · `{name:
  'M2', js: turn(2, rad(10), 1.4).M2, py: 'ch15.prandtl_meyer_turn(2.0, np.deg2rad(10.0), 1.4)["M2"]', rtol: 1e-9}` (2.38489) · `{name:
  'p2/p1', js: turn(2, rad(10), 1.4).p, py: 'ch15.prandtl_meyer_turn(2.0, np.deg2rad(10.0), 1.4)["p_ratio"]', rtol: 1e-9}` (0.54797)
  · `{name: 'inverse', js: machFromNu(nu(3, 1.4), 1.4), py: 'ch15.mach_from_prandtl_meyer(ch15.prandtl_meyer(3.0, 1.4), 1.4)', rtol:
  1e-9}` (3.0) · `{name: 'staircase 7', js: staircase(2, rad(10), 7, 1.4).M[7], py: 'ch15.fan_staircase(2.0, np.deg2rad(10.0), 7,
  1.4)["M"][-1]', rtol: 1e-10}` (2.38130) · `{name: 'mu2', js: turn(2, rad(10), 1.4).mu2, py: 'ch15.prandtl_meyer_turn(2.0,
  np.deg2rad(10.0), 1.4)["mu2"]', rtol: 1e-9}` (0.432682 rad = 24.791°) · invariant `{name: 'nu(1)', js: nu(1, 1.4), expect: 0, atol:
  1e-12}`.
- **Fit plan:** 360×640: status (1 line) · `corner` (50 %) over `nuview` (50 %); `stairs` hidden — "M₂ 2.385 exact · 2.381 by 7
  waves" in the `nuview` title; `y0`, `gamma` optional. 844×345: `corner` beside `nuview`. Desktop / 1000×700: rows [1.1, 0.9].

### E9 · supersonic_flat_plate
- **Title:** "How can there be drag without friction?" · **Summary:** "Tilt a plate in a supersonic stream: two flat pressure
  levels, a push square to the plate, and a drag that grows as α²." · **CORE:** C16 (also N102, N108–N112) · **Reference:**
  `forced_damped_vibrations.html` (system + graph + the Explain panel) and ch14 `airfoil_pressure_lift` (a pressure distribution
  whose area is the force).
- **meta:** `viz:order 9` · `viz:sections 15.10 15.9` · `viz:equations 15.55 15.57 15.58 15.59 14.12` · `viz:fluidpy
  ch15.ackeret_pressure ch15.ackeret_cp ch15.flat_plate_supersonic ch15.ackeret_flat_plate ch15.ackeret_airfoil
  ch15.flat_plate_wave_pattern ch15.lift_slope_vs_mach ch15.weak_wave_pressure ch15.oblique_shock ch15.prandtl_meyer_turn
  ch15.zhukhovsky_lift_coefficient` · `viz:derivations D28 D31`.
- **Physics (JS ↔ Python):** `ackeretP(delta, M, g)` ↔ `ch15.ackeret_pressure(delta, M_inf, gamma)`; `ackeretCp(delta, M)` ↔ `ch15.ackeret_cp(delta, M_inf)` (the angle comes first);
  `plateLinear(alpha, M, g)` → {CL, CD, p_upper, p_lower} ↔ `ch15.flat_plate_supersonic(alpha, M, g, method="linear")`;
  `plateExact(alpha, M, g)` (lower surface: weak oblique shock through δ = α; upper: Prandtl–Meyer turn α; $C_n=(p_3-p_2)/(\tfrac
  12\gamma M_\infty^2p_\infty)$, $C_L=C_n\cos\alpha$, $C_D=C_n\sin\alpha$; `attached: false` beyond detachment, where Python raises `ValueError("detached")`) ↔
  `ch15.flat_plate_supersonic(…, method="shock_expansion")`; `ackeretPlate(alpha, M)` ↔ `ch15.ackeret_flat_plate`; `thin(shape, tc,
  alpha, M)` → {CL, CD_incidence, CD_camber, CD_thickness} (diamond: $\overline{\theta_t^2}=(t/c)^2$; biconvex: $\tfrac43(t/c)^2$) ↔
  `ch15.ackeret_airfoil`; `waves(alpha, M, g)` ↔ `ch15.flat_plate_wave_pattern`; `slope(M)` ↔ `ch15.lift_slope_vs_mach`.
- **Views** (rows [1.1, 0.9]): 1. `plate` "The plate and its waves" (row 0, flex 1.4, `equal`): the plate (or diamond /
  biconvex section) at incidence α in a stream from the left; the four waves from `waves()` — the leading-edge fan above
  (blue rays) and shock below (rose), the trailing-edge shock above and fan below; pressure arrows on both surfaces (orange,
  length ∝ |p − p∞|, suction drawn outward); the resultant normal to the plate split into lift (teal) and wave drag (rose).
  2. `cp` "$C_p$ against x/c" (row 1, flex 1.2): two flat levels (upper negative, lower positive), linear (dashed) and exact
  (solid) according to `method`; the area between them shaded = the normal-force coefficient; in "compare" mode the subsonic
  loop of Ch. 14's section as a muted ghost (caption "Ch. 14, incompressible, same α"). 3. `slopeview` "Lift slope against
  Mach number" (row 1, flex 1, `hidePortrait`): 2π of $C_L\simeq2\pi(\alpha+\beta)$ on the subsonic side, $4/\sqrt{M_\infty^2-1}$ on the
  supersonic side, a muted gap between 0.8 and 1.2 labelled "neither theory", the cursor.
- **Controls:** `alphaDeg` "Angle of attack α" (−4…12°, step 0.1, default 2) · `Minf` "Free-stream Mach number $M_\infty$"
  (1.05…4, step 0.01, default 2; help "1.05: the linear theory blows up · 4: slender waves") · `method` (chips: linear rule ·
  exact shock–expansion · both) · `section` (modes: flat plate · diamond · biconvex) with `tc` "Thickness $t/c$" (0…0.1, default
  0.05) · `gamma` (optional).
- **Presets:** "2° at Mach 2" {alphaDeg 2, Minf 2} · "10°: each surface 20–25 % off, the lift 1 % off" {alphaDeg 10, Minf 2,
  method both} · "$M_\infty$ = 1.1: the theory fails" {alphaDeg 2, Minf 1.1, method both} · "a diamond at zero incidence: drag
  without lift" {section diamond, tc 0.05, alphaDeg 0}.
- **Status:** gap g = |linear $C_L$ / exact $C_L$ − 1|: g < 3 %: "linear theory within 3 % of exact (lift)" · 3 % ≤ g: "linear
  lift `…` % off" · exact detached: "⛔ the shock has detached: no shock–expansion solution; too close to Mach 1" · $M_\infty$ <
  1.2: "⚠️ close to Mach 1: the linear coefficient diverges".
- **Readouts:** "$C_L$ (linear / exact)" · "$C_D$ wave (linear / exact)" · "$L/D=1/\alpha$" · "Slope $4/\sqrt{M^2-1}$" · "$p_3/
  p_\infty$, $p_2/p_\infty$".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (three sections) · **terms** · **presets** ·
  **status**.
- **Terms:** title "Wave drag coefficient = incidence + camber + thickness"; items `inc` "$4\alpha^2/\sqrt{M^2-1}$" rose · `cam`
  "camber" muted (zero for these sections) · `thk` "thickness" amber; total "$C_D$"; a second line "lift: incidence only".
- **Explain:**
  0. *What the views show* — "**Top:** a thin plate tilted by α in a supersonic stream. Above the leading edge the stream
     turns away through a <b class=c-blue>fan</b>; below it turns into itself through a <b class=c-rose>shock</b>. <b
     class=c-orange>Orange arrows</b>: pressure on each surface. **Left graph:** the pressure coefficient along the chord;
     the shaded area is the force normal to the plate. b in the book's formulas is the chord."
  1. *Two pressures from the linear rule* — "$\dfrac{p-p_\infty}{p_\infty}=\dfrac{\gamma M_\infty^2\,\delta}{\sqrt{M_\infty^2-1}}$ (15.57) with δ = ∓α: (p − p∞)/p∞ = ∓1.4 × `…`² × `…` rad/√(`…`² − 1) = **∓`…`**."
  2. *The same two, exactly* — "lower: oblique shock through `…`°: $p_3/p_\infty$ = **`…`**; upper: Prandtl–Meyer turn of `…`°:
     $p_2/p_\infty$ = **`…`** (linear: `…` and `…`)."
  3. *Lift and drag* — "the force $(p_3-p_2)\,b$ is normal to the plate: L = force × cos α, D = force × sin α. $L\cong\dfrac{2\alpha\gamma M_\infty^2p_\infty b}{\sqrt{M_\infty^2-1}},\ \ D\cong\dfrac{2\alpha^2\gamma M_\infty^2p_\infty b}{\sqrt{M_\infty^2-1}}$ (15.58)."
  4. *Coefficients* — "$C_L\cong\dfrac{4\alpha}{\sqrt{M_\infty^2-1}},\ \ C_D\cong\dfrac{4\alpha^2}{\sqrt{M_\infty^2-1}}$ (15.59): $C_L$ = 4 × `…`/`…` = **`…`**; $C_D$ = 4 × `…`²/`…` = **`…`**; exact: `…`, `…`."
  5. *The price of lift* — "$C_D/C_L$ = α = `…` rad: L/D = **`…`** (before any friction)."
  6. *Against Chapter 14* — "slope 4/√(`…`² − 1) = **`…`** per radian; the slow-flow value of $C_L\simeq2\pi(\alpha+\beta)$ (14.12) is 2π = 6.28."
  7. *With thickness* (modes diamond, biconvex) — "$C_p=\pm\dfrac{2\theta}{\sqrt{M_\infty^2-1}},\ \ C_D=\dfrac{4}{\sqrt{M_\infty^2-1}}\left[\alpha^2+\overline{\theta_c^2}+\overline{\theta_t^2}\right]$ (ours): incidence `…` + thickness `…` = $C_D$ **`…`**; the lift is
     unchanged."
  8. *The slope view* (or the hint "turn the phone sideways; the slope and 2π are in the title of the $C_p$ view").
  9. *Reading the current setting* — small α, $M_\infty$ well above 1: "Each surface feels only its own slope, so the
     pressure is flat along the chord. No Kutta condition was needed to fix the lift. The push leans back by α: that is
     wave drag, and viscosity has nothing to do with it." · large α: "Surface by surface the linear rule is now 20–25 %
     off — too low on the compression side, too high on the suction side — but the two errors nearly cancel in the lift." ·
     near Mach 1: "$\sqrt{M_\infty^2-1}$ is going to zero and the linear coefficients run away; the exact solution stops
     existing when the lower shock detaches. Neither this theory nor Chapter 14's works here." · diamond at α = 0: "No lift
     at all — and still drag: the front half is pushed back by high pressure and the rear half pulled back by low
     pressure. Thickness buys nothing but drag; supersonic wings are thin."
- **Derivation tab:** **D28** (8 steps) `view: 'cp'` on phones; goal `set` {preset 2° at Mach 2, method both}; step 6 `live` "cot μ =
  √(`…`² − 1) = `…`"; step 8 `live` "(p₂ − p₁)/p₁ = 1.4 × `…`²/`…` × `…` = `…`". **D31** (7 steps) `view: 'plate'` on phones; goal `set`
  {preset 2° at Mach 2}; step 2 `watch` "the upper arrows point outward (suction), the lower ones inward"; step 4 `live` "L =
  `…` $p_\infty b$"; step 5 `watch` "the resultant leans back by exactly α"; step 7 `live` "$C_L$ = `…`, $C_D$ = `…`". Interpret: `s =>
  "At α = " + … + "° and M = " + … + " every unit of lift costs " + … + " of wave drag."`
- **Code:**
  ```python
  a = np.deg2rad({{alpha}})                                     # angle of attack [rad]
  lin = ch15.flat_plate_supersonic(a, {{M}}, {{g}}, method="linear")            # (15.59): CL {{CL}}, CD {{CD}}
  ex = ch15.flat_plate_supersonic(a, {{M}}, {{g}}, method="shock_expansion")    # exact: CL {{CLx}}, CD {{CDx}}
  ch15.ackeret_pressure(a, {{M}}, {{g}})                        # (15.57): (p - p_inf)/p_inf = {{dp}}
  ch15.lift_slope_vs_mach({{M}})                                # 4/sqrt(M^2 - 1) = {{slope}}  (2 pi below Mach 1)
  lin["CD"]/lin["CL"]                                           # = tan(alpha) ~ alpha = {{a}}
  ```
- **Walkthrough (7 steps):** 1. "A plate in a fast stream" — "No viscosity anywhere. Can it still feel drag?" `set` {preset 2° at
  Mach 2} · 2. "Two flat pressures" — "Above: a fan, lower pressure. Below: a shock, higher pressure. Each surface feels only its
  own slope." `eq: 'ack'`, `derive: {id: 'D31', step: 2}` · 3. "Square to the plate" — "The pressure difference pushes at right
  angles to the plate — so it leans back by α." `eq: 'LD'` · 4. "Lift ∝ α, drag ∝ α²" — "Tilt from 2° to 10°: the lift bar grows
  fivefold, the drag bar twenty-five-fold." `set` {alphaDeg 10}, `terms: true`, `eq: 'coef'` · 5. "How good is linear?" — "Switch
  on the exact solution. Each surface is far off; the lift is not." `set` {preset 10°…, method both}, `code: {id: 'plate', lines:
  [2, 3]}` · 6. "Against Chapter 14" — "Slope 4/√(M² − 1) here, 2π there — and a gap near Mach 1 where neither holds." `eq:
  'slope'`, `set` {Minf 1.3} · 7. "Your turn" — "Predict the drag of a 5 % diamond at zero incidence. Then switch section."
  `controls: ['alphaDeg', 'Minf', 'section']`.
- **Equations:** `weak` "Weak-wave rule" ref 'Eq. (15.55)' $\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta$ · `ack` "Pressure on an inclined surface" ref 'Eq. (15.57)'
  $\dfrac{p-p_\infty}{p_\infty}=\dfrac{\gamma M_\infty^2\,\delta}{\sqrt{M_\infty^2-1}}$, live, with the two surface values $\dfrac{p_2-p_\infty}{p_\infty}=-\dfrac{\gamma M_\infty^2\alpha}{\sqrt{M_\infty^2-1}},\ \ \dfrac{p_3-p_\infty}{p_\infty}=+\dfrac{\gamma M_\infty^2\alpha}{\sqrt{M_\infty^2-1}}$ · `LD` "Lift and wave drag per unit span" ref 'Eq. (15.58)' $L\cong\dfrac{2\alpha\gamma M_\infty^2p_\infty b}{\sqrt{M_\infty^2-1}},\ \ D\cong\dfrac{2\alpha^2\gamma M_\infty^2p_\infty b}{\sqrt{M_\infty^2-1}}$, live,
  note "b is the chord" · `coef` "Coefficients" ref 'Eq. (15.59)' $C_L\cong\dfrac{4\alpha}{\sqrt{M_\infty^2-1}},\ \ C_D\cong\dfrac{4\alpha^2}{\sqrt{M_\infty^2-1}}$, live · `thin` "A thin section of any shape" ref 'ours'
  $C_p=\pm\dfrac{2\theta}{\sqrt{M_\infty^2-1}},\ \ C_D=\dfrac{4}{\sqrt{M_\infty^2-1}}\left[\alpha^2+\overline{\theta_c^2}+\overline{\theta_t^2}\right]$, live · `slope` "The slow-flow slope" ref 'Eq. (14.12)' $C_L\simeq2\pi(\alpha+\beta)$.
- **Check yourself:** (1) "α = 2°, $M_\infty$ = 2: what is $C_L$?" — "0.081: $C_L\cong\dfrac{4\alpha}{\sqrt{M_\infty^2-1}},\ \ C_D\cong\dfrac{4\alpha^2}{\sqrt{M_\infty^2-1}}$ (15.59) gives 4 × 0.0349/1.732." `set {alphaDeg: 2, Minf:
  2}` · (2) "Double α. By what factor does the wave drag grow?" — "Four: $C_D\propto\alpha^2$, because the force grows with α
  *and* leans back by α." `set {alphaDeg: 4}` · (3) "Which surface is at lower pressure, and through which wave?" — "The upper:
  it turns away from the stream (δ = −α), through an expansion fan." · (4) "A symmetric diamond at α = 0: lift? drag?" — "No lift;
  wave drag 4(t/c)²/√(M² − 1): $C_p=\pm\dfrac{2\theta}{\sqrt{M_\infty^2-1}},\ \ C_D=\dfrac{4}{\sqrt{M_\infty^2-1}}\left[\alpha^2+\overline{\theta_c^2}+\overline{\theta_t^2}\right]$ (ours)." `set {section: 'diamond', alphaDeg: 0}`.
- **Selftest parity rows:** `{name: 'CL linear', js: ackeretPlate(rad(2), 2).CL, py: 'ch15.ackeret_flat_plate(np.deg2rad(2.0),
  2.0)["CL"]', rtol: 1e-12}` (0.080613) · `{name: 'dp/p', js: ackeretP(rad(2), 2, 1.4), py: 'ch15.ackeret_pressure(np.deg2rad(2.0), 2.0, 1.4)', rtol: 1e-12}` (0.112859) · `{name: 'Cp', js: ackeretCp(rad(2), 2), py: 'ch15.ackeret_cp(np.deg2rad(2.0), 2.0)', rtol: 1e-12}`
  (0.040307) · `{name: 'CL exact 10°', js: plateExact(rad(10), 2, 1.4).CL, py: 'ch15.flat_plate_supersonic(np.deg2rad(10.0), 2.0, 1.4,
  method="shock_expansion")["CL"]', rtol: 1e-7}` (0.40750) · `{name: 'CD exact 10°', js: plateExact(rad(10), 2, 1.4).CD, py:
  'ch15.flat_plate_supersonic(np.deg2rad(10.0), 2.0, 1.4, method="shock_expansion")["CD"]', rtol: 1e-7}` (0.07185) · `{name: 'p lower',
  js: plateExact(rad(10), 2, 1.4).p_lower, py: 'ch15.flat_plate_supersonic(np.deg2rad(10.0), 2.0, 1.4,
  method="shock_expansion")["p_lower"]', rtol: 1e-8}` (1.70658) · `{name: 'slope', js: slope(2), py: 'ch15.lift_slope_vs_mach(2.0)',
  rtol: 1e-12}` (2.309401) · invariants `{name: 'D = L tan(alpha)', js: plateExact(rad(5), 2, 1.4).CD/plateExact(rad(5), 2, 1.4).CL,
  expect: Math.tan(rad(5)), rtol: 1e-12}`, `{name: 'diamond', js: thin('diamond', 0.05, 0, 2).CD_thickness, expect: 0.0057735, rtol:
  1e-6}`.
- **Fit plan:** 360×640: status (1 line) · `plate` (50 %) over `cp` (50 %); `slopeview` hidden — "slope 2.31 per rad (slow flow:
  6.28)" in the `cp` title; `gamma`, `tc` optional. 844×345: `plate` beside `cp`. Desktop / 1000×700: rows [1.1, 0.9].

### B1 · wave_steepening
(backup — built only if one of E1–E9 fails review)
- **Title:** "Where do shocks come from?" · **Summary:** "A loud sine wave leans forward into a sawtooth; its wavelets' paths
  cross at the breaking time." · **CORE:** C02 (also N08, N22, N114, R16's shallow-water twin) · **Reference:**
  `angular_frequency_explorer_1.html` and ch07 `group_velocity_packets`.
- **meta:** `viz:order 10` · `viz:sections 15.2` · `viz:equations 15.7 15.16` · `viz:fluidpy ch15.simple_wave_speed
  ch15.simple_wave_evolve ch15.breaking_time` · `viz:derivations none` (the result is stated — curation §4c a-d08; if the
  orchestrator promotes it, the derivation is written then).
- **Physics (JS ↔ Python):** `wavelet(u, c0, g)` = c₀ + ½(γ + 1)u ↔ `ch15.simple_wave_speed`; `evolve(u0, x, t, c0, g)` (each
  sample moves at its own speed; invalid after the first crossing) ↔ `ch15.simple_wave_evolve`; `tBreak(u0, x, g)` ↔
  `ch15.breaking_time`. Mode "shallow water" uses γ = 2 with c₀ = √(gH) (the Ch. 7 twin).
- **Views:** 1. `wave` u(x, t) with the wavelets as dots moving at their own speeds, the acoustic (unsteepened) wave as a
  muted ghost; 2. `xt` the wavelets' straight paths on an x–t diagram (in the frame moving at c₀), their first crossing marked
  (amber); 3. `slope` (`hidePortrait`) the steepest slope against time, running to infinity at $t_{\text{break}}$.
- **Controls:** `amp` "Amplitude $u/c_0$" (0.001…0.2 log, default 0.02) · `shape` (chips: sine · single compression · single
  rarefaction) · `t` (transport, 0…1.2 $t_{\text{break}}$, `end: 'hold'`) · `system` (modes: a gas with γ · shallow water).
- **Presets:** "quiet" {amp 0.001} · "loud" {amp 0.05} · "a single hump" {shape compression} · "a rarefaction: flattens instead"
  {shape rarefaction}. · **Status:** "breaks after `…` wavelengths" / "expansion: never breaks" / "past breaking: the simple wave
  no longer applies — a shock has formed".
- **Depth features:** Explain + Code · **linked views** · **transport** · **modes** · **presets** · **status**.
- **Explain:** 0. *What the views show* (dots = wavelets; ghost = the wave as acoustics would have it). 1. *Each wavelet's
  speed* — "$\left.\dfrac{dx}{dt}\right\rvert_{\text{wavelet}}=u+c=c_0+\dfrac{\gamma+1}{2}u,\ \ t_{\text{break}}=\dfrac{2}{(\gamma+1)\max(-\partial u_0/\partial x)}$ (ours) (ours — the book argues in words): crest u = `…`: c₀ + 1.2 × `…` = `…` m/s; trough: `…` m/s." 2. *When do they
  cross?* — "steepest initial slope max(−∂u₀/∂x) = `…` 1/s: $t_{\text{break}}$ = 2/(2.4 × `…`) = **`…` s** = `…` wavelengths of
  travel." 3. *Was it ever acoustics?* — "$p'/\rho_0c^2\ll1$ (15.7): condensation `…`; the particle velocity of $u_1'(x,t)=\dfrac1{\rho_0c}\left(f(x-ct)-g(x+ct)\right)$ (15.16) was assumed tiny against c." 4.
  *Reading the current setting* — small amplitude: "It breaks only after hundreds of wavelengths: in practice damping wins
  and acoustics holds." · large: "It breaks within a few wavelengths: a sawtooth, then a train of weak shocks." · rarefaction:
  "In a rarefaction the part in front is the faster one, so the wave stretches and never breaks."
- **Code:** `ev = ch15.simple_wave_evolve(u0, x, {{t}}, {{c0}}, {{g}})   # valid: {{valid}}` · `ch15.breaking_time(u0, x, {{g}})   #
  {{tb}} s` · `ch15.simple_wave_speed({{u}}, {{c0}}, {{g}})   # {{w}} m/s`.
- **Walkthrough (5 steps):** 1. "A loud note" — "Sound keeps its shape — if it is weak. What if it is not?" · 2. "Hot parts
  hurry" — "Compressed gas is hotter and moves forward: crests travel faster than troughs." `eq: 'steep'` · 3. "Leaning forward"
  — "Play: the sine leans into a sawtooth." `play: true` · 4. "The paths cross" — "On the x–t map the wavelets' lines meet at the
  breaking time." · 5. "Your turn" — "Predict how the breaking time changes when you halve the amplitude."
- **Equations:** `steep` ref 'ours' $\left.\dfrac{dx}{dt}\right\rvert_{\text{wavelet}}=u+c=c_0+\dfrac{\gamma+1}{2}u,\ \ t_{\text{break}}=\dfrac{2}{(\gamma+1)\max(-\partial u_0/\partial x)}$ · `u` ref 'Eq. (15.16)' $u_1'(x,t)=\dfrac1{\rho_0c}\left(f(x-ct)-g(x+ct)\right)$ · `small` ref 'Eq. (15.7)' $p'/\rho_0c^2\ll1$.
- **Check yourself:** (1) "Halve the amplitude: what happens to the breaking time?" — "It doubles: $\left.\dfrac{dx}{dt}\right\rvert_{\text{wavelet}}=u+c=c_0+\dfrac{\gamma+1}{2}u,\ \ t_{\text{break}}=\dfrac{2}{(\gamma+1)\max(-\partial u_0/\partial x)}$ (ours) has the slope in
  the denominator." · (2) "Does a single rarefaction break?" — "No: its front is faster than its back." · (3) "In shallow water,
  which γ does the same formula use?" — "γ = 2."
- **Selftest parity rows:** `{name: 'wavelet', js: wavelet(5, 340, 1.4), py: 'ch15.simple_wave_speed(5.0, 340.0, 1.4)', rtol:
  1e-12}` (346.0) · `{name: 't_break', js: tBreak(u0, x, 1.4), py: 'ch15.breaking_time(2.0*np.sin(np.linspace(0, 2*np.pi, 2001)),
  np.linspace(0, 2*np.pi, 2001), 1.4)', rtol: 1e-4}` (0.41667) · invariant `{name: 'gamma 2', js: wavelet(1, 10, 2), expect: 11.5,
  rtol: 1e-12}`.
- **Fit plan:** 360×640: `wave` over `xt`; `slope` hidden, "breaks at t = … s" in the `xt` title.

---

## Part D — runtime budget (full run < 5 min on a laptop / Colab CPU; target 60–90 s)

Per-cell costs are estimates: measured where the designer ran the code (✓: the four sympy checks 0.1–0.8 s each; the whole
probe of the implementer's modules — about 150 calls including every sympy helper — ran in under 20 s; the designer's own
scratch numbers in 1 s) and from analysis §4 otherwise. The builder reports the measured numbers.

| Section | What costs time (per cell) | Full | FAST |
|---|---|---|---|
| front matter, §15.1 | imports (sympy, plotly, scipy, pandas); three tables; six recap cells | 5 s | 5 s |
| §15.2 C01, C02 | `linearised_acoustics_sympy` (< 2 s); **F1** (25 steps × 3 traces × 200 pts) 1 s; residual figure; **A1** (60 frames, dpi 80; FAST 30) 4 s / 2 s; x–t figure; **A2** (frames; 40; FAST 20) 3 s / 2 s; E1 iframe | 13 s | 9 s |
| §15.3 C03 | drawing; `duct_equations_sympy` (< 2 s); budget bars | 3 s | 3 s |
| §15.4 C04, C05 | sympy `series` (1 s); **F2** (20 steps × 2 traces × 200 pts) 1 s; ruler bars; generated table; **F3** (2 figures × 20 steps × 2 traces) 2 s | 6 s | 5 s |
| §15.5 C06 | two figures; sympy check of the Mach-number form (< 1 s); `nozzle_isentropic` × 3 on 401 stations (≈ 0.3 s: one `brentq` per station); **F4** (15 steps × 3 traces × 200 pts → 9000 root solves ≈ 2 s; FAST 8 steps); rocket example; E2 iframe | 6 s | 4 s |
| §15.6 C07, C08 | **D13** check (0.8 s ✓); `normal_shock_sympy` (< 5 s; cached in FAST — P252); **F5a**, **F5b** (20 steps × 4 traces × 200 pts, closed forms) 2 s; **D17** check (0.3 s ✓); `weak_shock_series` (< 5 s; cached in FAST); Hugoniot and entropy figures; profiles: `becker_shock_profile` at 3 × 200 points (`brentq` per point ≈ 1 s) + `solve_ivp` × 3 (0.3 s); `shock_structure_sympy` (< 2 s); E3 iframe | 18 s | 9 s |
| §15.7 C09, C10 | pressure families (9 curves × 401 stations ≈ 1 s); **F6** (2 figures × ≤ 19 steps × 1 trace, closed form) 1 s; **A3** (frames; 40 back pressures × 401 stations ≈ 16 000 root solves ≈ 3 s + encoding 3 s; FAST 20 frames × 201 stations); live cell (not executed on the page); E4 iframe | 10 s | 6 s |
| §15.8 C11, C12 | **D21** check (0.3 s ✓); **F7** (25 + 29 steps × 2 traces × 400 pts, closed forms) 2 s; finite-difference cell; **A5** (frames; 40; FAST 20) 3 s / 2 s; E5 iframe | 8 s | 6 s |
| §15.9 C13, C14, C15 | three panels; **A4** (video; 60 frames; FAST 30) 4 s / 2 s; envelope fit (200 circles); **F8** (20 steps × 1 trace × 400 pts) 1 s; wedge figures; live cell; `tangential_momentum_sympy`, `weak_wave_sympy` (< 2 s each); **D30** check (0.1 s ✓); `prandtl_meyer_sympy` (< 2 s); `quad`; fan and curved-wall figures; **F9b** (4 steps × 1 trace) < 1 s; E6, E7, E8 iframes | 16 s | 11 s |
| §15.10 C16 | plate figure; slope figure; **F9a** (15 steps × 4 traces × 60 pts; each exact point two root solves ≈ 2 s) ; E9 iframe | 5 s | 4 s |
| explainers | 9 `show_viz` cells (iframes; no computation) | 3 s | 3 s |
| **Total** | | **≈ 93 s** | **≈ 65 s** |

**FAST plan.** `FAST = setup_notebook()`; animations ≤ 60 frames at dpi 80 (FAST 30 or fewer): A1, A4 `player="video"`
(smooth motion); A2, A3, A5 `player="frames"` (the reader should stop at ¼, ½, 1 × the breaking time; at the seven named
nozzle conditions; at the paradox band and the nose). Slider figures ≤ 29 steps × ≤ 4 traces × ≤ 400 points (≈ 100–200 kB
each; thirteen figures in the nine groups F1–F9). Grids: nozzles on 401 stations with the throat a grid point (C.5 request
b); shock profiles 200 points; oblique chart `n` = 400. The four ★★★ `check_src` cells always run live (≈ 1.5 s together ✓).
The two slow sympy helpers (`normal_shock_sympy`, `weak_shock_series`) are wrapped in `functools.lru_cache` cells in FAST
(reminder P252). **No data file is committed and nothing is cached on disk** for this chapter. Page budget: five animations
< 1.2 MB each, thirteen plotly figures ≈ 150 kB each, nine explainers ≈ 250 kB each: ≈ 10 MB < 15 MB.

---

## Part E — prerequisite ledger
Every concept, symbol, maths tool and Python function or idiom the notebook or its explainers use, with where it is explained.
"primer (in Cxx)" = a 📎 primer placed in that block before first use (the Concept text begins with the primer term, used
verbatim in `nb.primer`); "knowledge/primers.md: <term> (chNN Pnn) — reminder" = a one-line reminder naming the earlier primer;
a CORE/RECAP id alone = taught there; "Cxx (Nnn)" = the NOTE placed in that block; "Cxx (D0n)" = the derivation where it is
built; "Cxx (gloss …)" = one sentence where it is used; "Cxx (🔁 Ch. n)" = an id-less recap cell in that block. Items of §15.1
are first used in the story that opens C01 and are entered as C01. New primers P354–P373 (20) in first-use order: P355 (C01,
with R05) · P356, P354 (C01) · P357, P358 (C02) · P359 (C04) · P360, P361 (C05) · P362 (C06) · P363 (C07) · P364 (C08) · P365
(C10) · P366 (C11) · P367 (C12) · P368, P369 (C13) · P370, P371, P372 (C14) · P373 (C15).

| Concept | First used in | Explained by |
|---|---|---|
| the conventions box: subscript 0 (three meanings) and $p_t$, $T_t$; the three stars; stations 1 → 2 against 0 → 1; σ, δ against β, θ; the sign of a turn; f, q against Fanno and Rayleigh parameters; the three wall forces; S against s; the lapse-rate line | C01 | front matter (A.0 row 5, before C01), repeated per block |
| the overloaded letters f, g, q, Q, A, δ, μ, ν, σ, U, V, v, b, c, m, e, h, η | C01 | front matter (A.0 row 5), and where each first appears |
| slips kept apart from traps; `ch15.book_slips()`, `ch15.traps()`, `ch15.conventions_table()`, `ch15.illustrative_inputs()` | C01 | front matter (A.0 row 5) |
| "not in the book — ours" labels; our own illustrative inputs (the six dicts) | C01 | front matter (A.0 rows 4–5) |
| gas dynamics; the standing assumptions (perfect gas, constant γ, no gravity, no friction outside boundary layers); the three promises | C01 | C01 (N01) |
| Mach number $M\equiv U/c$, local; `SIM.mach_number`, `GD.mach_from_speed` | C01 | R01 |
| speed of sound $c^2\equiv(\partial p/\partial\rho)_s$; `TH.sound_speed_from_eos` | C01 | R02 |
| compressibility of order $M^2$; `ch14.compressibility_check` | C01 | R03 |
| flow regimes (incompressible, subsonic, transonic, supersonic, hypersonic); `ch15.mach_regime` | C01 | C01 (N02) |
| history: Mach, Ackeret, von Kármán | C01 | C01 (N03) |
| perfect gas with constant specific heats: e, h, $C_v$, $C_p$, R, γ; `GD.PerfectGas`, `.cp`, `.cv`, `.c(T)` | C01 | R04 |
| entropy S per unit mass and its two logarithmic forms | C01 | R04 |
| isentropic relations between two states; `TH.isentropic_ratios` | C01 | R05 |
| the gamma exponents (3.5, 2.5, 3 for γ = 1.4; fractional powers in numpy) | C01 | primer (in C01) |
| properties of air; rounded against project constants; `GD.AIR_BOOK`, `TH.R_AIR` | C01 | R06 |
| acoustics: small isentropic disturbances about a steady state | C01 | C01 (N04) |
| steady reference state plus fluctuation, primes | C01 | C01 (N05; 🔁 Ch. 11 basic state + disturbance) |
| Taylor expansion to first and second order | C01 | knowledge/primers.md: Taylor expansion (ch01 P26) — reminder |
| orders of smallness; dropping products of primes | C01 | knowledge/primers.md: orders of smallness (ch02 P68) — reminder |
| partial derivative; index notation and the summation convention; the Laplacian $\partial^2/\partial x_j\partial x_j$ | C01 | knowledge/primers.md: partial derivative (ch01 P25) — reminder; C01 (🔁 Ch. 2) |
| product rule | C01 | knowledge/primers.md: product rule (ch01 P38) — reminder |
| continuity and Euler equations without gravity (the starting laws of D01) | C01 | C01 (D01; 🔁 Ch. 4) |
| acoustic relation $p'\cong c^2\rho'$; `AC.acoustic_density`, `AC.acoustic_temperature` | C01 | C01 (D01, N07) |
| sound speed at the reference state; "monatomic"; `TH.perfect_gas_sound_speed` | C01 | R07 |
| sound speed of air at 15 °C | C01 | R08 |
| condensation $p'/\rho_0c^2$; the size of the dropped term; `AC.condensation`, `AC.linearisation_error`, `AC.isentropic_pressure_perturbation` | C01 | C01 (N08) |
| decibels and sound pressure level; `AC.sound_pressure_level` | C01 | primer (in C01) |
| compression and rarefaction | C01 | C01 (N09) |
| isothermal against isentropic sound speed; `AC.newton_sound_speed` | C01 | C01 (N10; 🔁 Ch. 1) |
| linearised continuity and momentum; the convective operator | C01 | C01 (D01, N11–N14) |
| an operator in brackets, applied twice (constant-coefficient operators commute) | C01 | primer (in C01) |
| operator elimination for linear equations | C01 | knowledge/primers.md: operator elimination for linear PDEs (ch07 P177) — reminder |
| the wave equation in a uniform stream | C01 | C01 (D02) |
| `sympy` symbols, `diff`, `expand`, `simplify`, `series`, `limit`, `solve` | C01 | knowledge/primers.md: sympy (ch01 P40; ch04 P117) — reminder |
| `ch15.linearised_acoustics_sympy`, `AC.convected_wave_residual`, `AC.wave_equation_residual` | C01 | C01 (gloss in the code row) |
| central differences for a second derivative | C01 | C01 (gloss in the from-scratch row; 🔁 Ch. 10) |
| `assert np.allclose` | C01 | knowledge/primers.md: assert np.allclose (ch01 P15) — reminder |
| `slider_figure`, `animate` / `show_animation`, `show_viz`, `live` | C01 | knowledge/primers.md: slider_figure, animate, show_viz (ch01 P17, P16, P18) — reminder |
| logarithmic axes | C01 | C01 (gloss at F1) |
| the one-dimensional wave equation in a still fluid | C02 | C02 (N15) |
| the second-order wave equation and its two characteristic coordinates | C02 | primer (in C02) |
| first-order wave equation and characteristics | C02 | knowledge/primers.md: first-order wave equation and characteristics (ch07 P174) — reminder |
| chain rule (two new variables) | C02 | knowledge/primers.md: chain rule (ch01 P49) — reminder; primer P357 |
| d'Alembert's solution; f and g; `AC.dalembert`, `AC.pulse` | C02 | C02 (D03) |
| indefinite time integral and its "constant"; `AC.velocity_from_pressure` | C02 | C02 (N16) |
| particle velocity; acoustic impedance $\rho_0c$; `AC.plane_wave_velocity`, `AC.acoustic_impedance` | C02 | C02 (D04, N17) |
| initial-value form; `AC.dalembert_initial_value`, `AC.dalembert_fg` | C02 | C02 (N21) |
| `lambda` functions passed as shapes | C02 | C02 (gloss in the code row) |
| frames of reference and the Galilean shift | C02 | knowledge/primers.md: frames of reference and the Galilean shift (ch03 P96) — reminder |
| waves in a stream; speeds U ± c; `AC.dalembert_in_stream`, `AC.wave_speeds` | C02 | C02 (D05, N19) |
| range of influence; x–t diagram; `AC.xt_characteristics`, `AC.upstream_signal_possible` | C02 | C02 (N20) |
| domain of dependence | C02 | knowledge/primers.md: domain of dependence (ch10 P227) — reminder |
| nonlinear steepening, wavelet speed, breaking time; `AC.simple_wave_speed`, `AC.simple_wave_evolve`, `AC.breaking_time` | C02 | C02 (N22; 🔁 Ch. 7 simple wave) |
| Riemann invariant | C02 | primer (in C02) |
| the Euler equations in conservation form (named only) | C02 | C02 (N114) |
| quasi-one-dimensional duct flow; stations 1 and 2 | C03 | C03 (N23) |
| control volume and the Reynolds transport theorem | C03 | C03 (🔁 Ch. 3, Ch. 4) |
| mass conservation in a duct; ṁ; `GD.mass_flow` | C03 | R09 |
| steady energy integral; the control surface $A^*$ | C03 | R10 |
| power of a force; heat flux through a surface | C03 | knowledge/primers.md: power of a force and heat flux (ch04 P127) — reminder |
| enthalpy $h=e+p/\rho$; flow work | C03 | C03 (N26; 🔁 Ch. 1) |
| the duct energy equation; heat Q per unit mass; `GD.total_enthalpy`, `GD.energy_residual` | C03 | C03 (D06) |
| steady momentum integral | C03 | R11 |
| momentum flux through a surface | C03 | knowledge/primers.md: momentum flux through a surface (ch04 P114) — reminder |
| logarithmic differentiation | C03 | knowledge/primers.md: logarithmic differentiation (ch01 P52) — reminder |
| wall force F, friction per length $F_f$; `GD.wall_force` | C03 | C03 (D07, N27, N28) |
| one-dimensional Euler equation | C03 | R12 |
| Gibbs relation; adiabatic, reversible, isentropic | C03 | R13 |
| `ch15.duct_equations_sympy` | C03 | C03 (gloss in R12's code row) |
| a state as a Python `dict` with keys p, rho, u | C03 | C03 (gloss in the code row) |
| stagnation (total) state; $h_0$; `GD.total_enthalpy` | C04 | C04 (N29) |
| stagnation temperature; `BER.stagnation_temperature` | C04 | R14 |
| stagnation ratios $T_0/T$, $p_0/p$, $\rho_0/\rho$; `GD.stagnation_temperature_ratio`, `GD.stagnation_pressure_ratio`, `GD.stagnation_density_ratio`, `GD.isentropic_state`, `GD.mach_from_pressure_ratio` | C04 | C04 (D08, N30, N31) |
| dynamic pressure in a gas ($\tfrac12\rho U^2=\tfrac12\gamma pM^2$) | C04 | primer (in C04) |
| binomial approximation | C04 | knowledge/primers.md: binomial approximation (ch14 P342) — reminder |
| pitot formula and its compressibility factor; `GD.pitot_compressibility_factor`, `GD.mach_for_density_change`, `BER.stagnation_pressure` | C04 | C04 (N32; 🔁 Ch. 4) |
| entropy as a state function; logarithm rules | C04 | C04 (D09; 🔁 Ch. 1) |
| $T_0$ kept, $p_0$ lost; `GD.entropy_change_from_stagnation`, `GD.stagnation_pressure_loss` | C04 | C04 (D09, N33) |
| reservoir | C04 | C04 (N34) |
| maximum speed; `GD.max_speed` | C04 | C04 (N38) |
| sonic (critical) state $T^*$, $p^*$, $\rho^*$, $A^*$; `GD.sonic_ratios`, `GD.critical_pressure_ratio` | C05 | C05 (N35, D11) |
| the area ratio $A/A^*$; `GD.area_ratio` | C05 | C05 (D10) |
| minimum of a function by a zero derivative | C05 | C05 (gloss in D10 step 9) |
| a function with two inverses (one bracket on each side of the minimum) | C05 | primer (in C05) |
| `scipy.optimize.brentq` | C05 | knowledge/primers.md: brentq (ch03 P108) — reminder |
| the two branches; `GD.mach_from_area_ratio(branch=…)` | C05 | C05 (N36) |
| a generated table as a pandas DataFrame (`.loc`, `idxmin`, reading both ways) | C05 | primer (in C05) |
| `GD.isentropic_table`; NACA Report 1135 as the public check | C05 | C05 (N37) |
| mass flux function; choked mass flow; `GD.mass_flux_function`, `GD.choked_mass_flow`, `GD.mass_flow_from_mach` | C05 | C05 (D11, N39) |
| nozzle and diffuser | C06 | C06 (N40) |
| constant-density continuity in logarithmic form | C06 | R15 |
| differentials of algebraic laws as a small linear system | C06 | primer (in C06) |
| Euler's equation for isentropic changes (corrected); density against speed | C06 | C06 (D12, N41, N42) |
| the area–velocity relation; `GD.area_velocity`, `GD.influence_coefficients` | C06 | C06 (D12) |
| the four cases; `GD.duct_role` | C06 | C06 (N43) |
| throat; convergent–divergent nozzle; venturi; `GD.nozzle_isentropic`, `GD.laval_area` | C06 | C06 (N44) |
| l'Hôpital's rule (the 0/0 at the throat) | C06 | knowledge/primers.md: l'Hôpital's rule (ch14 P343) — reminder |
| Mach-number form of the area rule | C06 | C06 (N45) |
| `np.gradient` | C06 | C06 (gloss in the from-scratch row) |
| rocket nozzle chain; thrust with the pressure term; `ch15.rocket_nozzle`, `ch15.nozzle_thrust` | C06 | C06 (N46, N47) |
| standard atmosphere; Kundu Γ ≡ dT/dz and the meteorological −dT/dz | C06 | C06 (N46; 🔁 Ch. 1) |
| shock wave as a discontinuity | C07 | C07 (N48) |
| the three jump laws | C07 | C07 (N49, N50, N51) |
| dividing a known root out of a quadratic | C07 | primer (in C07) |
| Vieta's formulas | C07 | knowledge/primers.md: Vieta's formulas (ch02 P71) — reminder |
| velocity ratio $r=u_2/u_1$; $a=p_1/\rho_1u_1^2$ (auxiliary symbols of D13) | C07 | C07 (D13 steps 2, 6) |
| pressure jump; `GD.normal_shock`, `GD.normal_shock_pressure_ratio`, `ch15.normal_shock_sympy` | C07 | C07 (D13) |
| downstream Mach number; a map that is its own inverse; `GD.normal_shock_mach` | C07 | C07 (D14, N54, N55) |
| density and temperature jumps; `GD.normal_shock_density_ratio`, `GD.normal_shock_temperature_ratio` | C07 | C07 (D15) |
| `GD.normal_shock_residuals`, `GD.normal_shock_fluxes` | C07 | C07 (gloss in the code row) |
| `try` / `except ValueError` | C07 | C07 (gloss in the code row) |
| moving shock; frame-dependent quantities; `GD.moving_shock` | C07 | C07 (N64) |
| strong-shock limits; blast wave; `GD.strong_shock_limits`, `ch15.blast_constant_crude`, `ch01.blast_radius` | C07 | R17 |
| second law; entropy jump; `GD.normal_shock_entropy`, `GD.shock_admissible` | C08 | C08 (D16, N61) |
| the series of ln(1 + x) to third order (exact cancellation of low orders) | C08 | primer (in C08) |
| weak-shock entropy; ε = $M_1^2-1$; `GD.weak_shock_entropy`, `ch15.weak_shock_series` | C08 | C08 (D17, N60) |
| stagnation-pressure ratio across a shock; $A_1^*/A_2^*$; `GD.normal_shock_stagnation_pressure_ratio` | C08 | C08 (D18, N62) |
| Rankine–Hugoniot relations; Hugoniot curve; `GD.hugoniot_pressure_ratio` | C08 | C08 (N58) |
| `GD.normal_shock_table` | C08 | C08 (N63) |
| hydraulic jump and Froude number; `ch07.hydraulic_jump` | C08 | R16 |
| `TH.perfect_gas_entropy_change_p` | C08 | C08 (gloss in the from-scratch row; 🔁 Ch. 1) |
| viscosity sets thickness, not loss | C08 | R18 |
| shock thickness; mean free path, Knudsen number; `ch15.shock_thickness_estimate`, `ch01.mean_free_path_air` | C08 | C08 (N65; 🔁 Ch. 1) |
| steady 1-D viscous equations; bulk viscosity; $\mu''$; `ch15.shock_structure_sympy` | C08 | R19 |
| first integrals m, V, I; `ch15.shock_structure_constants` | C08 | C08 (N66) |
| Prandtl number $Pr''$; Becker's integral | C08 | C08 (N67) |
| scaled speed $U=u/V$, roots $U_1$, $U_2$; `ch15.becker_shock_roots`, `ch15.becker_shock_ode` | C08 | C08 (N68) |
| implicit profile, η; `ch15.becker_shock_profile`, `ch15.becker_shock_thickness` | C08 | C08 (N69, N70) |
| `scipy.integrate.solve_ivp` | C08 | knowledge/primers.md: solve_ivp (ch03 P94) — reminder |
| partial fractions; inverting an implicit solution | C08 | knowledge/primers.md: partial fractions and inverting an implicit solution (ch09 P219) — reminder |
| shock within a shock (named) | C08 | C08 (N71) |
| supply, back and exit-plane pressure | C09 | C09 (N72) |
| mass flow against back pressure; `GD.convergent_nozzle`, `GD.mass_flow_vs_back_pressure` | C09 | C09 (D19) |
| piecewise functions with `np.where` / `if` | C09 | knowledge/primers.md: np.where (ch01 P46) — reminder |
| choking; what "choked" does not mean | C09 | C09 (N74) |
| the seven nozzle conditions; $p_b$, $p_d$, $p_f$; `GD.laval_critical_pressures`, `GD.laval_regime` | C10 | C10 (D20) |
| nested root finding | C10 | primer (in C10) |
| shock position in a nozzle; `GD.laval_shock_location`, `GD.laval_exit_state`, `GD.laval_nozzle`, `ch15.nozzle_with_shock` | C10 | C10 (D20, N75, N78) |
| over-expanded and under-expanded jets | C10 | C10 (N76) |
| lumped friction f and heat q | C11 | C11 (N79) |
| the duct number A; `GD.duct_parameter_A` | C11 | C11 (D21, N80) |
| a quadratic in M squared and its discriminant | C11 | primer (in C11) |
| the quadratic formula | C11 | knowledge/primers.md: the quadratic formula (ch06 P159) — reminder |
| outlet Mach number; continuous and shocked roots; `GD.duct_outlet_mach`, `GD.duct_state_ratios`, `GD.duct_entropy_change`, `GD.duct_residuals`, `GD.duct_curve` | C11 | C11 (D21, N83) |
| friction differentials (corrected); `GD.influence_coefficients(df=…)` | C12 | C12 (D22, N84) |
| heating differentials (corrected); `GD.influence_coefficients(dq=…)` | C12 | C12 (D23, N85) |
| finite differences as an independent check | C12 | C12 (gloss in the from-scratch row) |
| a process line on a T–s diagram (Fanno and Rayleigh lines) | C12 | primer (in C12) |
| Fanno and Rayleigh ratios, friction length $4f_FL/D$; `GD.fanno_ratios`, `GD.rayleigh_ratios`, `GD.fanno_length_to_f`, `GD.fanno_line_Ts`, `GD.rayleigh_line_Ts` | C12 | C12 (N86, "not in the book — ours") |
| the temperature paradox band; `GD.rayleigh_temperature_sign` | C12 | C12 (N87) |
| choking limits; `GD.duct_choking_limits`, `GD.duct_choking_locus` | C12 | C12 (D24, N88) |
| funnel picture; dense-gas expansion shocks (named) | C12 | C12 (N89, N90) |
| wave fronts of a moving source; `AC.source_wavefronts` | C13 | C13 (N91) |
| the envelope of a family of circles (tangent ⟂ radius) | C13 | primer (in C13) |
| Mach angle μ; arcsine and its domain; `GD.mach_angle`, `GD.mach_from_angle` | C13 | C13 (D25) |
| Mach cone, Mach lines, zones of action and silence, boom delay; `AC.mach_cone_delay`, `AC.in_zone_of_action` | C13 | C13 (N92) |
| hyperbolic against elliptic equations | C13 | primer (in C13) |
| parabolic and elliptic problems | C13 | knowledge/primers.md: parabolic and elliptic problems (ch09 P201) — reminder |
| `np.polyfit` for a straight-line fit | C13 | C13 (gloss in the from-scratch row) |
| oblique shock geometry; σ and δ (β and θ elsewhere) | C14 | C14 (N93) |
| resolving a velocity along and across a tilted line | C14 | primer (in C14) |
| normal Mach numbers; `GD.oblique_shock_components`, `ch15.tangential_momentum_sympy` | C14 | C14 (D26, N94) |
| oblique-shock jumps; `GD.oblique_shock` | C14 | C14 (D26, N95–N97) |
| the tangent subtraction formula (and $\cos2\sigma=1-2\sin^2\sigma$) | C14 | primer (in C14) |
| deflection against shock angle; `GD.deflection_angle`, `GD.oblique_shock_curve` | C14 | C14 (D27) |
| weak and strong roots, maximum deflection, sonic point, detachment; `GD.max_deflection`, `GD.sonic_deflection`, `GD.detachment_mach` | C14 | C14 (N98) |
| wedge flow and the detached bow shock (sketch); `ch15.wedge_flow` | C14 | C14 (N99) |
| corner flows in subsonic flow (named) | C14 | C14 (N100; 🔁 Ch. 6) |
| the Mach triangle ($\cos\mu$, $\cot\mu$ from $\sin\mu=1/M$) | C14 | primer (in C14) |
| small-angle approximation | C14 | knowledge/primers.md: small-angle approximation (ch03 P100) — reminder |
| the weak-wave rule; `GD.weak_wave_pressure`, `GD.weak_wave_error`, `ch15.weak_wave_sympy` | C14 | C14 (D28, N101, N102) |
| `scipy.optimize.minimize_scalar` (cross-check of the maximum) | C14 | knowledge/primers.md: minimize_scalar (ch07 P170) — reminder |
| Prandtl–Meyer fan; simple wave in steady flow; `ch15.expansion_fan`, `ch15.curved_wall_mach_lines` | C15 | C15 (N104, N105) |
| the differential turning relation; `ch15.fan_staircase` | C15 | C15 (D29, N106) |
| substitution in an integral | C15 | knowledge/primers.md: substitution in an integral (ch03 P106) — reminder |
| arctangent integrals | C15 | primer (in C15) |
| Prandtl–Meyer function ν(M); `GD.prandtl_meyer`, `GD.prandtl_meyer_max`, `GD.mach_from_prandtl_meyer`, `ch15.prandtl_meyer_sympy` | C15 | C15 (D30) |
| the turning rule and both sign conventions; `GD.prandtl_meyer_turn` | C15 | C15 (N107) |
| `scipy.integrate.quad` | C15 | knowledge/primers.md: scipy.integrate.quad (ch03 P87) — reminder |
| weak expansions obey the weak-wave rule | C15 | C15 (N103) |
| Ackeret's rule; pressure coefficient; `ch15.ackeret_pressure`, `ch15.ackeret_cp` | C16 | C16 (N108; 🔁 Ch. 4) |
| flat plate: two surface pressures, lift, wave drag; chord b; `ch15.flat_plate_supersonic`, `ch15.ackeret_flat_plate`, `ch15.flat_plate_wave_pattern` | C16 | C16 (D31, N109, N110) |
| resolving a force along and across the stream | C16 | C16 (D31 step 4; 🔁 Ch. 14) |
| wave drag; d'Alembert's paradox; Kutta condition; lift slope 2π; `ch15.lift_slope_vs_mach`, `AF.zhukhovsky_lift_coefficient` | C16 | C16 (N111; 🔁 Ch. 6, Ch. 14) |
| thin section of any shape; `ch15.ackeret_airfoil` | C16 | C16 (N112) |
| shock–expansion method (exact comparison) | C16 | C16 (gloss in the code row; C14, C15) |
| caching a slow cell (`functools.lru_cache`) | C07 | knowledge/primers.md: caching (ch10 P252) — reminder |
| dataclasses | C01 | knowledge/primers.md: dataclasses (ch04 P111) — reminder |
| definite integral | C02 | knowledge/primers.md: definite integral (ch01 P27) — reminder |
| `np.expm1`, `np.log1p` and cancellation (inside the library) | C08 | knowledge/primers.md: np.expm1 / np.log1p and cancellation (ch03 P107) — reminder |
| f-strings; `np.vectorize` | C01 | C01 (gloss in the first code row that uses them) |

---

## Part F — derivation storyboards

One block per `D` row of the curation (§4b), in order. Builders copy each block word for word into `nb.derivation(id, title,
ref=…, goal=…, start=(tex, plain), plan=[…], uses=[…], steps=[dict(did=…, tex=…, why=…, plain=…), …], result=(tex, plain),
interpret=…, check=…, check_src=…)` and, where the heading names an explainer, into that explainer's `derivations: [...]` (same
`did` titles, same step count). Every step has four parts: **did** (the move) · **tex** (the new line) · **why** (why it is
allowed and why we make it) · **plain** (what the line says). Wherever a step produces a numbered book equation, the number
stands directly after that step's **tex** — beside the equation it names. Colours of terms as header convention 8. The four
★★★ blocks (D13, D17, D21, D30) carry a **sympy check** with a `check_src` listing that re-runs the derivation's own
construction; **each listing was executed as written for this design and passed** (the printed line is quoted under it).
Shorthand used in several blocks and defined where it first appears: $\Theta\equiv1+\tfrac{\gamma-1}{2}M^2$ (so that
$\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28) reads $T_0/T=\Theta$).

### D01 · The acoustic relation and the linearised laws: (15.5), (15.6), then (15.8), (15.9) → (15.10), (15.11) — ★★, 10 steps, in C01 (notebook)
- **Goal:** turn the full laws of mass and momentum into two *linear* equations for a weak disturbance, so that they can be
  solved. · **Start:** the pressure as a function of density at fixed entropy, $p=p(\rho)$; continuity
  $\dfrac{\partial\rho}{\partial t}+\dfrac{\partial(\rho u_j)}{\partial x_j}=0$ and Euler's equation without gravity
  $\dfrac{\partial u_i}{\partial t}+u_j\dfrac{\partial u_i}{\partial x_j}+\dfrac1\rho\dfrac{\partial p}{\partial x_i}=0$ (Ch. 4); the split $u_i=U_i+u_i',\ p=p_0+p',\ \rho=\rho_0+\rho',\ T=T_0+T'$ (15.4)
  — "steady uniform state plus a small disturbance". · **Plan:** • expand p(ρ) and keep one term • put the split into
  continuity, drop products of primes • the same for Euler • trade ρ′ for p′. · **Tools:** Taylor expansion (P26 reminder);
  product rule (P38 reminder); orders of smallness (P68 reminder); index notation (Ch. 2 recap); basic state + disturbance
  (Ch. 11 recap). · **Assumptions:** isentropic disturbance (step 1); perfect gas (step 2); the base state $U_j$, $p_0$, $\rho_0$
  is uniform and steady (steps 5, 7); primes are small (steps 2, 5, 8); inviscid, no gravity (start).
- **Steps:**
  1. **did** Expand p(ρ) about the undisturbed density · **tex**
     $p=p_0+\left(\dfrac{\partial p}{\partial\rho}\right)_s\rho'+\dfrac12\left(\dfrac{\partial^2p}{\partial\rho^2}\right)_s\rho'^2+\dots$ · **why** Taylor's theorem
     for a smooth function, with derivatives taken at $\rho_0$ and at constant entropy (sound is too quick for heat to
     flow). We want $p'$ in terms of $\rho'$. · **plain** A small change of density gives a pressure change proportional to
     it, plus a smaller correction.
  2. **did** Keep the first-order term · **tex** $p'\cong c^2\rho'$ (15.5), with $c^2=\gamma p_0/\rho_0=\gamma RT_0$ (15.6) · **why**
     $\rho'\ll\rho_0$, so the squared term is dropped (≈). For a perfect gas $p\propto\rho^\gamma$ at fixed entropy, so the
     slope is $\gamma p_0/\rho_0$ — which is $c^2$ by $c^2\equiv(\partial p/\partial\rho)_s$ (1.19). · **plain** Pressure and density disturbances are proportional;
     the constant is the sound speed squared.
  3. **did** Put the split into continuity · **tex**
     $\dfrac{\partial(\rho_0+\rho')}{\partial t}+\dfrac{\partial\left[(\rho_0+\rho')(U_j+u_j')\right]}{\partial x_j}=0$ · **why** Substitution of $u_i=U_i+u_i',\ p=p_0+p',\ \rho=\rho_0+\rho',\ T=T_0+T'$ (15.4); nothing is
     approximated yet. · **plain** Mass conservation, written for "steady state plus disturbance".
  4. **did** Multiply out the mass flux · **tex** $(\rho_0+\rho')(U_j+u_j')=\rho_0U_j+\rho'U_j+\rho_0u_j'+\rho'u_j'$ · **why** Plain
     expansion of a product; it separates one constant term, two terms with one prime and one with two. · **plain** The
     flux has a steady part, two first-order parts and one second-order part.
  5. **did** Drop what vanishes and what is second order · **tex**
     $\dfrac{\partial\rho'}{\partial t}+U_j\dfrac{\partial\rho'}{\partial x_j}+\rho_0\dfrac{\partial u_j'}{\partial x_j}=0$ (15.8, last member) · **why** $\rho_0$ and
     $U_j$ are constants, so their derivatives are zero and they come out of the derivatives; $\rho'u_j'$ is a product of
     two small things and is dropped. · **plain** Density changes because the stream carries the disturbance along and
     because the disturbance velocity converges.
  6. **did** Put the split into Euler's equation · **tex**
     $\dfrac{\partial(U_i+u_i')}{\partial t}+(U_j+u_j')\dfrac{\partial(U_i+u_i')}{\partial x_j}+\dfrac{1}{\rho_0+\rho'}\dfrac{\partial(p_0+p')}{\partial x_i}=0$ · **why**
     Substitution of $u_i=U_i+u_i',\ p=p_0+p',\ \rho=\rho_0+\rho',\ T=T_0+T'$ (15.4) again. (The book writes $P$ and $\rho_o$ here; they are $p_0$ and $\rho_0$.) · **plain** Momentum, for
     "steady state plus disturbance".
  7. **did** Remove the derivatives of constants · **tex**
     $\dfrac{\partial u_i'}{\partial t}+(U_j+u_j')\dfrac{\partial u_i'}{\partial x_j}+\dfrac{1}{\rho_0+\rho'}\dfrac{\partial p'}{\partial x_i}=0$ · **why** $U_i$ and $p_0$
     do not vary in space or time, so only the primed parts survive differentiation. · **plain** Only the disturbance is
     accelerated, and only the disturbance pressure pushes.
  8. **did** Keep first order · **tex**
     $\dfrac{\partial u_i'}{\partial t}+U_j\dfrac{\partial u_i'}{\partial x_j}+\dfrac1{\rho_0}\dfrac{\partial p'}{\partial x_i}=0$ (15.9, last member) · **why**
     $u_j'\,\partial u_i'/\partial x_j$ has two primes; and $1/(\rho_0+\rho')=1/\rho_0-\rho'/\rho_0^2+\dots$, whose second term
     times $\partial p'/\partial x_i$ has two primes too. · **plain** The disturbance velocity changes because the stream
     carries it and because of the pressure gradient.
  9. **did** Trade ρ′ for p′ in step 5 · **tex** $\dfrac1{c^2}\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)p'+\rho_0\dfrac{\partial u_j'}{\partial x_j}=0$ (15.10) · **why** Step 2 gives $\rho'=p'/c^2$ with c constant, so it
     passes through the derivatives. (The book cites $c^2=\gamma RT_0$ (15.6) here; the relation used is $p'\cong c^2\rho'$ (15.5).) We want
     one unknown fewer. · **plain** Continuity, now for the pressure disturbance.
  10. **did** Group the two derivative terms of step 8 · **tex** $\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)u_i'+\dfrac1{\rho_0}\dfrac{\partial p'}{\partial x_i}=0$ (15.11) · **why** Only a regrouping: the same bracket
      "time derivative plus stream times space derivative" now acts on $u_i'$ as it does on $p'$ in step 9 — that is what
      lets D02 eliminate the velocity. · **plain** Following the stream, the disturbance velocity changes only because of
      the pressure gradient.
- **Result:** $p'\cong c^2\rho'$ (15.5) with $c^2=\gamma p_0/\rho_0=\gamma RT_0$ (15.6); $\dfrac1{c^2}\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)p'+\rho_0\dfrac{\partial u_j'}{\partial x_j}=0$ (15.10); $\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)u_i'+\dfrac1{\rho_0}\dfrac{\partial p'}{\partial x_i}=0$ (15.11) — "two linear equations for $p'$ and $u_i'$, with the same convective bracket in both".
- **Check:** units — the continuity line (step 9): (s²/m²)(Pa/s) = kg/(m³ s) ✓ and (kg/m³)(1/s) ✓; the momentum line (step 10): m/s² in both terms ✓. Limit — $U_j=0$:
  $\dfrac1{c^2}\dfrac{\partial p'}{\partial t}+\rho_0\dfrac{\partial u_j'}{\partial x_j}=0$, the still-air form. Number — a 1 Pa sound in sea-level air:
  ρ′ = 1/340.3² = 8.6 × 10⁻⁶ kg/m³. Code — `ch15.linearised_acoustics_sympy()` lists exactly the dropped terms.
- **What it means:** to first order a disturbance does not interact with itself — different sounds pass through each
  other. **Fails when:** the amplitude is not small ($p'/\rho_0c^2$ no longer ≪ 1: steepening, N22), the base flow is not
  uniform, or heat conduction matters (very high frequency).
- **Traps:** the subscript 0 is the undisturbed state, not stagnation (trap T1); zeroth-order terms vanish because the base
  state is *uniform*, not because it is small; $1/(\rho_0+\rho')\approx1/\rho_0$ holds only to first order; it is $p'\cong c^2\rho'$ (15.5), not $c^2=\gamma RT_0$ (15.6), that removes ρ′ (slip #3).

### D02 · The wave equation for sound in a uniform stream (15.12) — ★★, 6 steps, in C01 (notebook)
- **Goal:** combine the two linear equations into one equation for the pressure alone. · **Start:** $\dfrac1{c^2}\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)p'+\rho_0\dfrac{\partial u_j'}{\partial x_j}=0$ (15.10) and $\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)u_i'+\dfrac1{\rho_0}\dfrac{\partial p'}{\partial x_i}=0$ (15.11) — "one
  equation with a time-like derivative of $p'$ and a divergence of $u'$; one with a time-like derivative of $u'$ and a
  gradient of $p'$". · **Plan:** • name the bracket • apply it to the first equation • take the divergence of the second
  • subtract. · **Tools:** an operator in brackets, applied twice (primer P354); operator elimination (P177 reminder); the
  Laplacian (Ch. 2 recap). · **Assumptions:** $U_j$, c and $\rho_0$ are constants (step 4).
- **Steps:**
  1. **did** Name the convective bracket · **tex** $\mathcal D\equiv\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}$ · **why** A shorthand, so
     that the two equations read $\dfrac1{c^2}\mathcal Dp'+\rho_0\dfrac{\partial u_j'}{\partial x_j}=0$ and $\mathcal Du_i'+\dfrac1{\rho_0}\dfrac{\partial p'}{\partial x_i}=0$.
     · **plain** $\mathcal D$ is the rate of change seen by someone drifting with the stream.
  2. **did** Apply $\mathcal D$ to the first equation · **tex**
     $\dfrac1{c^2}\mathcal D^2p'+\rho_0\,\mathcal D\!\left(\dfrac{\partial u_j'}{\partial x_j}\right)=0$ · **why** An equation that holds everywhere
     may be differentiated; c and $\rho_0$ are constants and pass through. We do it to create a term we can cancel. ·
     **plain** The stream-following rate of the continuity equation.
  3. **did** Take the divergence of the second equation, times $\rho_0$ · **tex**
     $\rho_0\dfrac{\partial}{\partial x_i}\!\left(\mathcal Du_i'\right)+\dfrac{\partial^2p'}{\partial x_i\partial x_i}=0$ · **why** Differentiate with
     respect to $x_i$ and sum over i; multiplying by the constant $\rho_0$ makes the velocity term match step 2. ·
     **plain** How fast the disturbance velocity converges is tied to the curvature of the pressure.
  4. **did** Swap the order of the two operations · **tex**
     $\rho_0\,\mathcal D\!\left(\dfrac{\partial u_i'}{\partial x_i}\right)+\dfrac{\partial^2p'}{\partial x_i\partial x_i}=0$ · **why** $\mathcal D$ has
     constant coefficients, and partial derivatives of a smooth function may be taken in any order. Now the velocity term
     is identical to that of step 2. · **plain** The same statement with the derivatives reordered.
  5. **did** Subtract step 4 from step 2 · **tex** $\dfrac1{c^2}\mathcal D^2p'-\dfrac{\partial^2p'}{\partial x_j\partial x_j}=0$ · **why** The
     two velocity terms are equal and cancel; the name of a summed index (i or j) does not matter. · **plain** One
     equation for the pressure alone.
  6. **did** Write the bracket out again · **tex** $\dfrac1{c^2}\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)^2p'-\dfrac{\partial^2p'}{\partial x_j\partial x_j}=0$ (15.12) · **why** Only notation; the square means "apply the bracket
     twice" (primer P354). · **plain** In a frame moving with the stream, the pressure obeys the ordinary wave equation
     with speed c.
- **Result:** $\dfrac1{c^2}\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)^2p'-\dfrac{\partial^2p'}{\partial x_j\partial x_j}=0$ (15.12) — "the field equation of sound in a uniform stream".
- **Check:** units — both terms Pa/m² ✓. Limit — $U_j=0$, one dimension: $\dfrac1{c^2}\dfrac{\partial^2p'}{\partial t^2}-\dfrac{\partial^2p'}{\partial x^2}=0$ (15.13). Number — $p'=\cos(x-(U+c)t)$ with c = 340, U =
  170 m/s satisfies it (worked example in C01); the speed 400 m/s does not. Code — `AC.convected_wave_residual` ≈ 0 for a
  pulse riding at U + c.
- **What it means:** the only speeds in the equation are U ± c (in one dimension) — this is where every later $1-M^2$
  comes from. **Fails when:** $U_j$ or c vary in space (then the operators do not commute and extra terms appear).
- **Traps:** the squared bracket is an operator applied twice, not a squared number; the cancellation needs *constant*
  $U_j$; apply $\mathcal D$ to the first equation and the divergence to the second, not the other way round.

### D03 · d'Alembert's solution (15.15) — ★★, 8 steps, in C02 (notebook · `sound_pulse_dalembert`)
- **Goal:** find *every* solution of the one-dimensional wave equation. · **Start:** $\dfrac1{c^2}\dfrac{\partial^2p'}{\partial t^2}-\dfrac{\partial^2p'}{\partial x^2}=0$ (15.13) — "the pressure's acceleration in
  time is c² times its curvature in space". · **Plan:** • change to the coordinates that ride with the two waves • the
  equation collapses to one mixed derivative • integrate twice. · **Tools:** the second-order wave equation and its two
  characteristic coordinates (primer P357); chain rule (P49 reminder). · **Assumptions:** still fluid, one space dimension,
  constant c.
- **Steps:**
  1. **did** Introduce two new coordinates · **tex** $\xi=x-ct,\qquad\zeta=x+ct$ · **why** An invertible change of variables is always allowed. A shape moving right at speed c is constant along ξ = const, one moving left along ζ = const. · **plain** ξ labels points riding with a right-going
     wave; ζ with a left-going one.
  2. **did** Translate the first derivatives · **tex**
     $\dfrac{\partial}{\partial x}=\dfrac{\partial}{\partial\xi}+\dfrac{\partial}{\partial\zeta},\qquad\dfrac{\partial}{\partial t}=-c\dfrac{\partial}{\partial\xi}+c\dfrac{\partial}{\partial\zeta}$ ·
     **why** Chain rule: $\partial\xi/\partial x=\partial\zeta/\partial x=1$, $\partial\xi/\partial t=-c$, $\partial\zeta/\partial t=+c$. · **plain** A
     step in x moves both labels forward; a step in t moves ξ back and ζ forward.
  3. **did** Apply the x-rule twice · **tex**
     $\dfrac{\partial^2p'}{\partial x^2}=\dfrac{\partial^2p'}{\partial\xi^2}+2\dfrac{\partial^2p'}{\partial\xi\,\partial\zeta}+\dfrac{\partial^2p'}{\partial\zeta^2}$ · **why** Square of
     the operator in step 2 (it has constant coefficients, so it can be applied twice like an algebraic square). · **plain**
     The curvature in x, in the new coordinates.
  4. **did** Apply the t-rule twice · **tex**
     $\dfrac{\partial^2p'}{\partial t^2}=c^2\left(\dfrac{\partial^2p'}{\partial\xi^2}-2\dfrac{\partial^2p'}{\partial\xi\,\partial\zeta}+\dfrac{\partial^2p'}{\partial\zeta^2}\right)$ ·
     **why** The same, with the factor −c on one term: the mixed term gets a minus sign. · **plain** The acceleration in
     t, in the new coordinates.
  5. **did** Substitute both into the wave equation · **tex** $-4\dfrac{\partial^2p'}{\partial\xi\,\partial\zeta}=0$ · **why** Divide step 4 by
     $c^2$ and subtract step 3: the pure second derivatives cancel, the mixed terms add. · **plain** Everything cancels
     except the mixed derivative.
  6. **did** Divide by −4 · **tex** $\dfrac{\partial}{\partial\zeta}\!\left(\dfrac{\partial p'}{\partial\xi}\right)=0$ · **why** −4 is not zero. Written this
     way the equation says something simple about $\partial p'/\partial\xi$. · **plain** The ξ-slope of the pressure does not
     change when you move in ζ.
  7. **did** Integrate in ζ · **tex** $\dfrac{\partial p'}{\partial\xi}=F(\xi)$ · **why** A quantity whose ζ-derivative is zero can
     depend on ξ only; the "constant" of a partial integration is an arbitrary *function* of the other variable. ·
     **plain** The ξ-slope is some function of ξ alone.
  8. **did** Integrate in ξ · **tex** $p'=f(\xi)+g(\zeta)=f(x-ct)+g(x+ct)$ (15.15) · **why** $f$ is an antiderivative of F; the
     constant of this integration is an arbitrary function of ζ, called g. Then put ξ and ζ back. · **plain** Any solution
     is one shape moving right plus one shape moving left.
- **Result:** $p'(x,t)=f(x-ct)+g(x+ct)$ (15.15) — "two arbitrary shapes, each keeping its form, travelling at ±c".
- **Check:** units — f and g are pressures ✓. Substitution — `sp.diff(f(x - c*t) + g(x + c*t), t, 2)/c**2 - sp.diff(…, x, 2)`
  simplifies to 0 (primer P357's demo). Number — c = 340 m/s: after 1 s the two shapes are 680 m apart.
- **What it means:** sound does not spread or fade in one dimension; a microphone hears the shape that was emitted.
  **Fails when:** the wave is not plane (in three dimensions it spreads and weakens), the medium moves or varies, or
  the amplitude is finite.
- **Traps:** the two second derivatives differ only in the sign of the mixed term; integrating a partial derivative gives
  a function, not a number; f and g are *any* shapes, not only sines.

### D04 · The particle velocity of a plane sound wave (15.16) — ★, 5 steps, in C02 (notebook · `sound_pulse_dalembert`)
- **Goal:** find how the air itself moves under the two travelling shapes. · **Start:** $u_1'(x,t)=-\dfrac1{\rho_0}\displaystyle\int\dfrac{\partial p'}{\partial x}\,dt$ (15.14) (the momentum law integrated in
  time) and $p'(x,t)=f(x-ct)+g(x+ct)$ (15.15). · **Plan:** • differentiate the solution in x • turn x-derivatives into t-derivatives • integrate
  exactly • fix the "constant". · **Tools:** the chain rule for a function of $x\mp ct$ (primer P357); the integral of an
  exact time derivative. · **Assumptions:** the fluid was at rest before the wave arrived (step 5).
- **Steps:**
  1. **did** Differentiate the solution in x · **tex** $\dfrac{\partial p'}{\partial x}=f'(x-ct)+g'(x+ct)$ · **why** Chain rule; the
     prime on f and g means the derivative with respect to their own argument. · **plain** The pressure gradient is the
     sum of the two shapes' slopes.
  2. **did** Rewrite each slope as a time derivative · **tex**
     $f'(x-ct)=-\dfrac1c\dfrac{\partial f}{\partial t},\qquad g'(x+ct)=+\dfrac1c\dfrac{\partial g}{\partial t}$ · **why** Chain rule again:
     $\partial f/\partial t=-c\,f'$ and $\partial g/\partial t=+c\,g'$. We do it because $u_1'(x,t)=-\dfrac1{\rho_0}\displaystyle\int\dfrac{\partial p'}{\partial x}\,dt$ (15.14) integrates over time. ·
     **plain** For a travelling shape, a slope in space is a rate in time — with opposite signs for the two directions.
  3. **did** Insert into the time integral · **tex**
     $u_1'=-\dfrac1{\rho_0}\displaystyle\int\left(-\dfrac1c\dfrac{\partial f}{\partial t}+\dfrac1c\dfrac{\partial g}{\partial t}\right)dt$ · **why** Substitution of
     steps 1–2 into $u_1'(x,t)=-\dfrac1{\rho_0}\displaystyle\int\dfrac{\partial p'}{\partial x}\,dt$ (15.14). · **plain** The velocity is the time-integral of two exact rates.
  4. **did** Integrate · **tex** $u_1'=\dfrac1{\rho_0c}\left(f-g\right)+C(x)$ · **why** The integral of a time derivative is
     the function itself; an indefinite time integral leaves an unknown that may depend on x. · **plain** The velocity is
     f minus g over $\rho_0c$, up to something fixed in time.
  5. **did** Fix the unknown · **tex** $u_1'(x,t)=\dfrac1{\rho_0c}\left(f(x-ct)-g(x+ct)\right)$ (15.16) · **why** Before the
     wave arrives at x both f and g are zero there and the air is at rest, so $C(x)=0$ for every x. · **plain** The air
     moves forward under a right-going compression and backward under a left-going one.
- **Result:** $u_1'(x,t)=\dfrac1{\rho_0c}\left(f(x-ct)-g(x+ct)\right)$ (15.16) — "pressure is f + g; velocity is (f − g)/ρ₀c".
- **Check:** units — Pa/(kg m⁻² s⁻¹) = m/s ✓. Limit — a bump at rest (f = g): u′ = 0 at t = 0 ✓. Number — 50 Pa in air:
  50/(1.2 × 340) = 0.12 m/s. Code — `AC.velocity_from_pressure` (quadrature of $u_1'(x,t)=-\dfrac1{\rho_0}\displaystyle\int\dfrac{\partial p'}{\partial x}\,dt$ (15.14)) equals `AC.plane_wave_velocity`.
- **What it means:** $u'/c=p'/\rho_0c^2$ — the air's speed is to the sound speed as the pressure disturbance is to
  $\rho_0c^2$: tiny. **Fails when:** the air was already moving, or the wave is not plane.
- **Traps:** the minus sign in front of g comes from $\partial g/\partial t=+c\,g'$ (trap T19); the integration "constant"
  is a function of x; u′ is not the wave speed.

### D05 · Sound in a uniform stream (15.17) — ★, 5 steps, in C02 (notebook · `sound_pulse_dalembert`)
- **Goal:** find how the two shapes travel when the whole fluid moves at speed U. · **Start:** $\dfrac1{c^2}\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)^2p'-\dfrac{\partial^2p'}{\partial x_j\partial x_j}=0$ (15.12) with $U_j=(U,0,0)$ in
  one dimension: $\dfrac1{c^2}\left(\dfrac{\partial}{\partial t}+U\dfrac{\partial}{\partial x}\right)^2p'-\dfrac{\partial^2p'}{\partial x^2}=0$. · **Plan:** • ride with the
  stream • there the fluid is still, so D03 applies • come back. · **Tools:** the Galilean shift (P96 reminder); D03. ·
  **Assumptions:** U and c constant; U along +x.
- **Steps:**
  1. **did** Ride with the stream · **tex** $x'=x-Ut,\qquad t'=t$ · **why** A change to a frame moving at constant speed U;
     Newtonian mechanics is the same in it. We choose it because the fluid is at rest there. · **plain** $x'$ is the
     position measured by someone drifting with the fluid.
  2. **did** Translate the derivatives · **tex**
     $\dfrac{\partial}{\partial x}=\dfrac{\partial}{\partial x'},\qquad\dfrac{\partial}{\partial t}+U\dfrac{\partial}{\partial x}=\dfrac{\partial}{\partial t'}$ · **why** Chain rule:
     $\partial x'/\partial x=1$, $\partial x'/\partial t=-U$, so $\partial/\partial t=\partial/\partial t'-U\,\partial/\partial x'$. · **plain** The
     convective bracket is just the time derivative seen by the drifting observer.
  3. **did** Rewrite the equation · **tex** $\dfrac1{c^2}\dfrac{\partial^2p'}{\partial t'^2}-\dfrac{\partial^2p'}{\partial x'^2}=0$ · **why** Substitute step
     2; this is $\dfrac1{c^2}\dfrac{\partial^2p'}{\partial t^2}-\dfrac{\partial^2p'}{\partial x^2}=0$ (15.13) in the primed variables. · **plain** In the drifting frame it is the still-fluid wave
     equation.
  4. **did** Use d'Alembert's solution there · **tex** $p'=f(x'-ct')+g(x'+ct')$ · **why** D03 found every solution of that
     equation. · **plain** In the drifting frame: one shape each way at speed c.
  5. **did** Come back to the fixed frame · **tex** $p'(x,t)=f(x-(c+U)t)+g(x+(c-U)t)$ (15.17) · **why** Put $x'=x-Ut$:
     $x-Ut-ct=x-(c+U)t$ and $x-Ut+ct=x+(c-U)t$. · **plain** One shape travels at U + c; the other at U − c, which is
     upstream only while U < c.
- **Result:** $p'(x,t)=f(x-(c+U)t)+g(x+(c-U)t)$ (15.17) — "the stream simply carries both shapes along".
- **Check:** units — (c ± U)t is a length ✓. Limit — U = 0: $p'(x,t)=f(x-ct)+g(x+ct)$ (15.15) ✓. Number — M = 0.5 at sea level: +510 and −170 m/s; M = 2:
  +1021 and +340 m/s. Code — `AC.wave_speeds`.
- **What it means:** at U = c the upstream-running shape stands still; beyond it nothing travels upstream — the root of
  choking (C09), shocks (C07) and the Mach cone (C13). **Fails when:** the stream is not uniform.
- **Traps:** g's speed is U − c, whose sign changes at U = c; "downstream" presumes U > 0 (trap T19); the shapes f and g are
  unchanged — only their speeds.

### D06 · The energy equation of duct flow: (15.20) → (15.21) → (15.22) — ★★, 8 steps, in C03 (notebook)
- **Goal:** turn the general energy integral into one line between two stations of a duct. · **Start:** $\displaystyle\int_{A^*}\rho\left(e+\tfrac12u_i^2\right)u_jn_j\,dA=\int_{A^*}u_i\tau_{ij}n_j\,dA-\int_{A^*}q_jn_j\,dA$ (15.20) — "energy
  carried out through the surface = work done on the surface + heat conducted in" (steady, no body force). · **Plan:** •
  split the surface into inlet, outlet, wall • evaluate each integral on each part • divide by the mass flow • spot the
  enthalpy. · **Tools:** surface integrals on a simple control volume (Ch. 4 recap); power of a force and heat flux (P127
  reminder); enthalpy (Ch. 1 recap). · **Assumptions:** steady (start); one value of u, p, ρ, e per cross-section (steps
  2–4); no slip and no through-flow at the wall (steps 3, 4); viscous normal stress negligible on the end faces (step 4).
- **Steps:**
  1. **did** Split the control surface · **tex** $A^*=A_1\cup A_2\cup A_{\text{wall}}$ · **why** An integral over a surface is
     the sum of the integrals over its parts. (Here $A^*$ is the control surface — not the sonic area of later sections.) ·
     **plain** Inlet face, outlet face, side wall.
  2. **did** Energy flux through the inlet · **tex**
     $\displaystyle\int_{A_1}\rho\left(e+\tfrac12u^2\right)u_jn_j\,dA=-\left(e+\tfrac12u^2\right)_1\rho_1u_1A_1=-\left(e+\tfrac12u^2\right)_1\dot m$ · **why**
     Everything is uniform over the face and the outward normal points *against* the flow there, so $u_jn_j=-u_1$; then
     $\rho uA=\dot m=\text{const.}$ (15.18). · **plain** Energy enters with the mass flow — a negative outflow.
  3. **did** Outlet and wall · **tex**
     $\displaystyle\int_{A_2}(\dots)=+\left(e+\tfrac12u^2\right)_2\dot m,\qquad\int_{A_{\text{wall}}}(\dots)=0$ · **why** On the outlet the normal points
     with the flow; no fluid crosses the wall ($u_jn_j=0$). · **plain** Energy leaves with the same mass flow; none through
     the wall.
  4. **did** Work of the surface stresses · **tex** $\displaystyle\int_{A^*}u_i\tau_{ij}n_j\,dA=(upA)_1-(upA)_2$ · **why** On the end faces
     the stress is the pressure, $\tau_{ij}n_j=-p\,n_i$, so $u_i\tau_{ij}n_j=-p\,u_in_i$: $+pu$ at the inlet, $-pu$ at the
     outlet. At the wall the velocity is zero (no slip): no work. · **plain** The surroundings push fluid in at station 1
     and are pushed by it at station 2.
  5. **did** Name the heat · **tex** $-\left(e+\tfrac12u^2\right)_1\dot m+\left(e+\tfrac12u^2\right)_2\dot m=(upA)_1-(upA)_2+\dot mQ$ (15.21), with $\dot mQ\equiv-\displaystyle\int_{A^*}q_jn_j\,dA$ · **why** A definition: Q is the
     heat *added* per unit mass of fluid passing through. Collecting steps 2–4 gives the line. · **plain** Energy out
     minus energy in = pressure work + heat added.
  6. **did** Divide by the mass flow · **tex**
     $\left(e+\tfrac12u^2\right)_2-\left(e+\tfrac12u^2\right)_1=\left(\dfrac p\rho\right)_1-\left(\dfrac p\rho\right)_2+Q$ · **why** $\dot m=\rho uA\neq0$, so
     $upA/\dot m=p/\rho$ — the *flow work* per kilogram. · **plain** Per kilogram: the energy gained equals net flow work
     plus heat.
  7. **did** Bring each station's terms together · **tex**
     $\left(e+\dfrac p\rho+\dfrac12u^2\right)_2-\left(e+\dfrac p\rho+\dfrac12u^2\right)_1=Q$ · **why** Move the two $p/\rho$ terms to the left. ·
     **plain** One combination at the outlet minus the same at the inlet equals the heat added.
  8. **did** Recognise the enthalpy · **tex** $h_2+\tfrac12u_2^2-h_1-\tfrac12u_1^2=Q$ (15.22) · **why** By definition
     $h=e+p/\rho$. · **plain** Enthalpy plus kinetic energy per kilogram changes only by the heat added.
- **Result:** $h_2+\tfrac12u_2^2-h_1-\tfrac12u_1^2=Q$ (15.22) — "the total enthalpy changes only by heat".
- **Check:** units — J/kg throughout ✓. Limit — Q = 0, perfect gas: $C_pT_0=C_pT+\tfrac12u^2$ (15.27) when station 2 is at rest. Number — air from 100 to
  300 m/s with Q = 0: $C_p\Delta T$ = ½(10⁴ − 9 × 10⁴) → ΔT = −40 K. Code — `GD.energy_residual` ≈ 0 along our nozzle.
- **What it means:** nothing in the derivation asked whether the flow *between* the stations is frictionless or smooth —
  so the line holds across friction and across shocks. That is why $T_0$ is the chapter's fixed point. **Fails when:**
  the flow is unsteady, a shaft does work, gravity matters, or the profiles at the stations are far from uniform.
- **Traps:** the normal points against the flow at the inlet; the wall does no work because of no slip, even with
  friction; Q is heat *added*; $A^*$ here is the control surface (trap T2); adiabatic (Q = 0) is not isentropic (trap T17).

### D07 · Mass and momentum in a duct: (15.18) → (15.19); (15.23) → (15.24) → (15.25) → (15.26) — ★★, 9 steps, in C03 (notebook)
- **Goal:** write mass and momentum for a duct in the differential forms every later section uses. · **Start:** $\rho uA=\dot m=\text{const.}$ (15.18) and
  $\displaystyle\int_{A^*}\rho u_iu_jn_j\,dA=\int_{A^*}\tau_{ij}n_j\,dA$ (15.23). · **Plan:** • differentiate the mass law logarithmically • evaluate the momentum integral between two stations •
  shrink the distance to dx and split the wall force • drop friction. · **Tools:** logarithmic differentiation (P52
  reminder); product rule (P38 reminder); momentum flux (P114 reminder). · **Assumptions:** steady, quasi-one-dimensional
  (all steps); frictionless from step 8.
- **Steps:**
  1. **did** Differentiate the mass law · **tex** $uA\,d\rho+\rho A\,du+\rho u\,dA=0$ · **why** $\rho uA$ is the same at every
     station, so its differential is zero; product rule for three factors. · **plain** Small changes of density, speed
     and area must compensate.
  2. **did** Divide by ρuA · **tex** $\dfrac{d\rho}{\rho}+\dfrac{du}{u}+\dfrac{dA}{A}=0$ (15.19) · **why** $\rho uA=\dot m\neq0$. Dividing turns absolute changes into
     relative ones. · **plain** The percentage changes of density, speed and area add to zero.
  3. **did** Evaluate the momentum integral · **tex** $-\dot mu_1+\dot mu_2=(pA)_1-(pA)_2+F$ (15.24) · **why** As in D06: momentum flux $-\dot mu_1$ in,
     $+\dot mu_2$ out; pressure pushes forward on face 1 and backward on face 2; F collects everything the walls do to
     the fluid in +x. · **plain** Change of momentum = pressure forces on the ends + force of the walls.
  4. **did** Take the stations a distance dx apart · **tex** $\dot m\,du=-d(pA)+dF$ · **why** The differences
     $u_2-u_1$, $(pA)_2-(pA)_1$ and F become differentials as $x_2\to x_1+dx$. · **plain** The same balance for a thin
     slice.
  5. **did** Split the wall force on the slice · **tex** $dF=p\,dA-F_f\,dx$ · **why** A sloping wall has pressure p acting on
     the projected area dA, pushing the fluid forward when the duct widens; friction retards with $F_f$ per unit length
     (positive by definition). · **plain** The walls push with their pressure and drag with their friction.
  6. **did** Divide by dx · **tex** $\dot m\dfrac{du}{dx}=-\dfrac{d}{dx}(pA)+p\dfrac{dA}{dx}-F_f$ · **why** Insert step 5 into step 4. ·
     **plain** Momentum change per unit length.
  7. **did** Expand the product · **tex** $\dot m\dfrac{du}{dx}=-A\dfrac{dp}{dx}-F_f$ (15.25) · **why** Product rule:
     $d(pA)/dx=A\,dp/dx+p\,dA/dx$; the $p\,dA/dx$ terms cancel. · **plain** Only the pressure *gradient* and friction
     accelerate the fluid; the area change drops out.
  8. **did** Remove friction and insert the mass flow · **tex** $\rho uA\dfrac{du}{dx}=-A\dfrac{dp}{dx}$ · **why** Frictionless:
     $F_f=0$; and $\dot m=\rho uA$. · **plain** Mass times acceleration equals the pressure-gradient force.
  9. **did** Divide by ρA and multiply by dx · **tex** $u\,du+\dfrac{dp}{\rho}=0$ (15.26) · **why** $\rho A\neq0$. This is
     Euler's equation along the duct. · **plain** Where the speed rises the pressure falls.
- **Result:** $\dfrac{d\rho}{\rho}+\dfrac{du}{u}+\dfrac{dA}{A}=0$ (15.19); $-\dot mu_1+\dot mu_2=(pA)_1-(pA)_2+F$ (15.24); $\dot m\dfrac{du}{dx}=-A\dfrac{dp}{dx}-F_f$ (15.25); $u\,du+\dfrac{dp}{\rho}=0$ (15.26).
- **Check:** units — the momentum line of step 7: (kg/s)(1/s) = N/m ✓. Limit — constant ρ and A: du = 0. Number — our nozzle between the
  station of area 20 cm² and the throat: F = −426 N. Code — `ch15.duct_equations_sympy()`.
- **What it means:** area does not appear in the momentum line — it acts only through mass conservation. **Fails when:**
  the area changes abruptly (the flow is no longer one-dimensional).
- **Traps:** F is the force on the fluid (trap T11); forgetting the wall's pressure term leaves a spurious $p\,dA$;
  $u\,du+\dfrac{dp}{\rho}=0$ (15.26) needs $F_f=0$.

### D08 · The stagnation ratios (15.28), (15.29), (15.30) — ★, 6 steps, in C04 (notebook)
- **Goal:** express the stagnation temperature, pressure and density as multiples of the local values, using the Mach
  number only. · **Start:** $C_pT_0=C_pT+\tfrac12u^2$ (15.27) — "the heat content at rest equals the heat content plus the kinetic energy in motion". ·
  **Plan:** • divide by $C_pT$ • insert $C_p$ and the sound speed • use the isentropic relations for p and ρ. · **Tools:**
  the gamma exponents (primer P355); $C_v=\dfrac{R}{\gamma-1},\ C_p=\dfrac{\gamma R}{\gamma-1}$ (15.1, specific-heat line); the isentropic relations $\dfrac{p_2}{p_1}=\left(\dfrac{\rho_2}{\rho_1}\right)^{\gamma},\ \dfrac{T_2}{T_1}=\left(\dfrac{\rho_2}{\rho_1}\right)^{\gamma-1}=\left(\dfrac{p_2}{p_1}\right)^{(\gamma-1)/\gamma}$ (15.2) (R05). · **Assumptions:** perfect gas, constant
  specific heats; adiabatic (steps 1–3); the imagined deceleration is also frictionless, hence isentropic (steps 4–6).
- **Steps:**
  1. **did** Divide by $C_pT$ · **tex** $\dfrac{T_0}{T}=1+\dfrac{u^2}{2C_pT}$ · **why** $C_pT\neq0$. We want a ratio. · **plain**
     The stagnation temperature exceeds the local one by the kinetic energy measured in units of the enthalpy.
  2. **did** Insert the specific heat · **tex** $\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}\dfrac{u^2}{\gamma RT}$ · **why**
     $C_p=\gamma R/(\gamma-1)$ for a perfect gas. It produces the combination γRT. · **plain** The same, with γ and R.
  3. **did** Recognise the Mach number · **tex** $\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28) · **why** $c^2=\gamma RT$ with the
     *local* temperature, and $M=u/c$. · **plain** The temperature ratio depends on the local Mach number alone.
  4. **did** Use the isentropic pressure–temperature relation · **tex** $\dfrac{p_0}{p}=\left(\dfrac{T_0}{T}\right)^{\gamma/(\gamma-1)}$ ·
     **why** The stagnation state is *defined* by an isentropic deceleration, so $\dfrac{p_2}{p_1}=\left(\dfrac{\rho_2}{\rho_1}\right)^{\gamma},\ \dfrac{T_2}{T_1}=\left(\dfrac{\rho_2}{\rho_1}\right)^{\gamma-1}=\left(\dfrac{p_2}{p_1}\right)^{(\gamma-1)/\gamma}$ (15.2) links the two states. · **plain**
     Pressure rises with the 3.5th power of the temperature ratio (γ = 1.4).
  5. **did** Insert step 3 · **tex** $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29) · **why** Substitution. · **plain** The stagnation pressure ratio, in
     terms of M.
  6. **did** The same for the density · **tex** $\dfrac{\rho_0}{\rho}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{1/(\gamma-1)}$ (15.30) · **why** From $\dfrac{p_2}{p_1}=\left(\dfrac{\rho_2}{\rho_1}\right)^{\gamma},\ \dfrac{T_2}{T_1}=\left(\dfrac{\rho_2}{\rho_1}\right)^{\gamma-1}=\left(\dfrac{p_2}{p_1}\right)^{(\gamma-1)/\gamma}$ (15.2), $\rho_0/\rho=(T_0/T)^{1/(\gamma-1)}$; insert
     step 3. · **plain** The density rises with the 2.5th power (γ = 1.4).
- **Result:** $\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28), $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29), $\dfrac{\rho_0}{\rho}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{1/(\gamma-1)}$ (15.30).
- **Check:** units — pure numbers ✓. Limits — M = 0: all ratios 1; small M: $p_0-p\approx\tfrac12\gamma pM^2=\tfrac12\rho u^2$
  (the pitot formula). Consistency — $p_0/p=(\rho_0/\rho)(T_0/T)$: exponents 1/(γ − 1) + 1 = γ/(γ − 1) ✓. Number — M = 2:
  1.8, 7.82, 4.35.
- **What it means:** every point of a flow carries a $T_0$ and a $p_0$; the first is kept by any adiabatic flow, the second
  only by an isentropic one (D09). **Fails when:** γ varies (very high temperature), or the gas is not perfect.
- **Traps:** $T_0/T$ needs only "adiabatic", the other two also need "frictionless"; the Mach number uses the local c; the
  exponents are γ/(γ − 1) and 1/(γ − 1), not their inverses.

### D09 · $T_0$ is kept, $p_0$ is lost — ★, 5 steps, in C04 (notebook)
- **Goal:** show that when entropy rises in an adiabatic flow, the stagnation pressure must fall. · **Start:** $S_2-S_1=C_p\ln\dfrac{T_2}{T_1}-R\ln\dfrac{p_2}{p_1}=C_v\ln\dfrac{T_2}{T_1}-R\ln\dfrac{\rho_2}{\rho_1}$ (15.1, entropy line) — "the
  entropy difference between any two states of a perfect gas". · **Plan:** • move from the two flow states to their
  stagnation states • apply the entropy formula there • use $T_{02}=T_{01}$. · **Tools:** entropy is a state function (Ch. 1
  recap); logarithm rules. · **Assumptions:** perfect gas; adiabatic between stations 1 and 2 (step 4).
- **Steps:**
  1. **did** Compare each state with its own stagnation state · **tex** $S_{01}=S_1,\qquad S_{02}=S_2$ · **why** The
     stagnation state is reached, in imagination, by an isentropic process — entropy unchanged by definition. · **plain**
     A point and its stagnation state have the same entropy.
  2. **did** Subtract · **tex** $S_2-S_1=S_{02}-S_{01}$ · **why** Step 1, twice. · **plain** The entropy change of the flow
     equals the entropy change between the two stagnation states.
  3. **did** Apply the entropy formula between the stagnation states · **tex**
     $S_2-S_1=-R\ln\dfrac{p_{02}}{p_{01}}+C_p\ln\dfrac{T_{02}}{T_{01}}$ (unnumbered, §15.4) · **why** $S_2-S_1=C_p\ln\dfrac{T_2}{T_1}-R\ln\dfrac{p_2}{p_1}=C_v\ln\dfrac{T_2}{T_1}-R\ln\dfrac{\rho_2}{\rho_1}$ holds between *any* two
     states because S is a state function; choose states 01 and 02. · **plain** Entropy change in terms of stagnation
     quantities.
  4. **did** Use the energy law · **tex** $S_2-S_1=-R\ln\dfrac{p_{02}}{p_{01}}$ · **why** Adiabatic flow: $h_2+\tfrac12u_2^2-h_1-\tfrac12u_1^2=Q$ (15.22) with Q = 0 gives
     $h_{02}=h_{01}$, so $T_{02}=T_{01}$ and $\ln1=0$ — whether or not there is friction. · **plain** With no heat added,
     all the entropy change shows up in the stagnation pressure.
  5. **did** Solve for the pressure ratio · **tex** $\dfrac{p_{02}}{p_{01}}=e^{-(S_2-S_1)/R}<1\ \text{ when }\ S_2>S_1$ · **why**
     Exponentiate; the second law says entropy cannot fall in an adiabatic flow. · **plain** Every bit of entropy made
     costs stagnation pressure.
- **Result:** $S_2-S_1=-R\ln\dfrac{p_{02}}{p_{01}}+C_p\ln\dfrac{T_{02}}{T_{01}}$ (unnumbered, §15.4), hence $T_{02}=T_{01}$ and $S_2>S_1$ give $p_{02}<p_{01}$.
- **Check:** units — S/R is a number ✓. Limit — isentropic: $p_{02}=p_{01}$. Number — a 10 % loss of $p_0$: ΔS = −287 ln 0.9 =
  30.2 J/(kg K).
- **What it means:** $p_0$ is the flow's "ability to do work"; friction and shocks spend it without changing the energy.
  **Fails when:** heat is exchanged ($T_0$ changes, C12).
- **Traps:** adiabatic is not isentropic (trap T17); $p_0$ can fall with $T_0$ unchanged; the formula is between
  stagnation states, not static ones.

### D10 · The area ratio (15.31) — ★★, 9 steps, in C05 (notebook · `area_mach_nozzle`) — the book states it without proof
- **Goal:** find how much area a given flow needs at Mach number M, compared with the area it would need at Mach 1. ·
  **Start:** $\rho uA=\rho^*c^*A^*$ — $\rho uA=\dot m=\text{const.}$ (15.18) written between the station and the (imagined) sonic station, where $u=c^*$. ·
  **Plan:** • write the area ratio as a chain of ratios • evaluate each from D08 at M and at 1 • collect the exponents
  • find the minimum. · **Tools:** a chain of ratios; the gamma exponents (primer P355); D08; a minimum by a zero
  derivative. · **Assumptions:** isentropic between the station and the imagined sonic station (so $T_0$, $\rho_0$ are
  common to both).
- **Steps:**
  1. **did** Solve for the area ratio · **tex** $\dfrac{A}{A^*}=\dfrac{\rho^*}{\rho}\,\dfrac{c^*}{u}$ · **why** Divide the mass
     balance by $\rho uA^*$. · **plain** More area is needed where the mass flux ρu is smaller.
  2. **did** Insert the stagnation density and the Mach number · **tex**
     $\dfrac{A}{A^*}=\dfrac{\rho^*}{\rho_0}\,\dfrac{\rho_0}{\rho}\,\dfrac{c^*}{c}\,\dfrac1M$ · **why** Multiply and divide by $\rho_0$ and by c;
     $u=Mc$. Each factor is now a ratio D08 gives. · **plain** Four ratios, each depending on one Mach number.
  3. **did** Write the sound-speed ratio with temperatures · **tex**
     $\dfrac{c^*}{c}=\sqrt{\dfrac{T^*}{T}}=\sqrt{\dfrac{T^*}{T_0}\,\dfrac{T_0}{T}}$ · **why** $c=\sqrt{\gamma RT}$; multiply and divide by
     $T_0$, which is the same at both stations. · **plain** Sound speed goes as the square root of temperature.
  4. **did** Evaluate at the station (Mach number M) · **tex**
     $\dfrac{T_0}{T}=\Theta,\quad\dfrac{\rho_0}{\rho}=\Theta^{1/(\gamma-1)},\qquad\Theta\equiv1+\dfrac{\gamma-1}{2}M^2$ · **why** $\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28) and $\dfrac{\rho_0}{\rho}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{1/(\gamma-1)}$ (15.30). Θ is
     a shorthand we keep for the rest of the chapter. · **plain** The local stagnation ratios.
  5. **did** Evaluate at the sonic station (M = 1) · **tex**
     $\dfrac{T^*}{T_0}=\dfrac{2}{\gamma+1},\qquad\dfrac{\rho^*}{\rho_0}=\left(\dfrac{2}{\gamma+1}\right)^{1/(\gamma-1)}$ · **why** The same two equations with
     M = 1: Θ = (γ + 1)/2. · **plain** The sonic state as fractions of the stagnation state (0.833 and 0.634 for air).
  6. **did** Assemble · **tex**
     $\dfrac{A}{A^*}=\dfrac1M\left(\dfrac{2}{\gamma+1}\right)^{\frac1{\gamma-1}}\Theta^{\frac1{\gamma-1}}\left(\dfrac{2}{\gamma+1}\right)^{\frac12}\Theta^{\frac12}$ · **why**
     Substitute steps 3–5 into step 2. · **plain** Everything in terms of M and γ.
  7. **did** Add the exponents · **tex** $\dfrac1{\gamma-1}+\dfrac12=\dfrac{\gamma+1}{2(\gamma-1)}$ · **why** Powers of the same base
     multiply by adding exponents; common denominator 2(γ − 1). · **plain** One combined exponent — 3 for air.
  8. **did** Write the result · **tex** $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31) · **why** Both $\frac{2}{\gamma+1}$ and Θ carry the exponent of step
     7. · **plain** The area ratio as a function of Mach number.
  9. **did** Find where it is smallest · **tex** $\dfrac{d}{dM}\ln\dfrac{A}{A^*}=\dfrac{M^2-1}{M\,\Theta}$ · **why** Logarithmic
     derivative: $-\dfrac1M+\dfrac{\gamma+1}{2}\dfrac{M}{\Theta}$, on a common denominator. It is negative for M < 1, zero at 1,
     positive above. · **plain** The area ratio falls to exactly 1 at Mach 1 and rises on both sides.
- **Result:** $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31) — "a U-shaped curve with its minimum, 1, at M = 1".
- **Check:** units — a pure number ✓. Limits — M = 1: 1; M → 0: $A/A^*\to\infty$ like 1/M (a reservoir); M → ∞: grows like
  $M^{2/(\gamma-1)}$. Number — M = 2: (1/2)(1.5)³ = 1.6875; M = 0.5: 2 × 0.875³ = 1.340. Code — `GD.area_ratio`.
- **What it means:** each area ratio above 1 belongs to two Mach numbers, one subsonic and one supersonic; which occurs is
  decided downstream. A given throat passes a given mass flow at most. **Fails when:** entropy changes between the
  station and the throat (a shock, friction): then $A^*$ itself changes.
- **Traps:** $c^*/c$ is a square root; the star is the sonic state here (trap T2); two inverses (trap T8).

### D11 · The sonic ratios and the choked mass flow — ★★, 7 steps, in C05 (notebook) — ours
- **Goal:** find the pressure, temperature and density at a sonic throat, and the mass flow such a throat passes. ·
  **Start:** $\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28), $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29), $\dfrac{\rho_0}{\rho}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{1/(\gamma-1)}$ (15.30) and $\rho uA=\dot m=\text{const.}$ (15.18). · **Plan:** • set M = 1 in the three ratios • write $\dot m=\rho^*c^*A^*$ • express
  it with the reservoir values • see that it is a maximum. · **Tools:** substitution; $p=\rho RT$; D10 step 9. ·
  **Assumptions:** isentropic from the reservoir to the throat; perfect gas.
- **Steps:**
  1. **did** Set M = 1 in the temperature ratio · **tex** $\dfrac{T^*}{T_0}=\dfrac{2}{\gamma+1}$ · **why** $\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28) with M = 1 gives
     $T_0/T^*=(\gamma+1)/2$. · **plain** At a sonic throat the gas has cooled to 5/6 of the reservoir temperature (air).
  2. **did** The same in the other two · **tex**
     $\dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)},\qquad\dfrac{\rho^*}{\rho_0}=\left(\dfrac{2}{\gamma+1}\right)^{1/(\gamma-1)}$ · **why** $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29) and
     $\dfrac{\rho_0}{\rho}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{1/(\gamma-1)}$ (15.30) with M = 1. · **plain** The pressure has fallen to 52.8 % and the density to 63.4 % (air).
  3. **did** Write the mass flow at the throat · **tex** $\dot m=\rho^*c^*A^*$ · **why** $\rho uA=\dot m=\text{const.}$ (15.18) at the station where the speed
     equals the local sound speed. · **plain** Mass flow = sonic density × sonic speed × throat area.
  4. **did** Express density and speed with reservoir values · **tex**
     $\rho^*=\dfrac{\rho^*}{\rho_0}\,\dfrac{p_0}{RT_0},\qquad c^*=\sqrt{\gamma RT^*}=\sqrt{\gamma RT_0}\left(\dfrac{2}{\gamma+1}\right)^{1/2}$ · **why** The gas law
     at the reservoir, $\rho_0=p_0/RT_0$, and steps 1–2. · **plain** Everything traced back to what the reservoir holds.
  5. **did** Multiply · **tex**
     $\dot m=p_0A^*\sqrt{\dfrac{\gamma}{RT_0}}\left(\dfrac{2}{\gamma+1}\right)^{\frac1{\gamma-1}+\frac12}$ · **why** $\sqrt{\gamma RT_0}/(RT_0)=\sqrt{\gamma/(RT_0)}$;
     the two powers of 2/(γ + 1) multiply. · **plain** The throat's mass flow in terms of $p_0$, $T_0$ and its area.
  6. **did** Combine the exponents · **tex** $\dot m_{\max}=p_0A^*\sqrt{\dfrac{\gamma}{RT_0}}\left(\dfrac{2}{\gamma+1}\right)^{(\gamma+1)/(2(\gamma-1))}$
     (ours) · **why** The same sum as in D10 step 7. · **plain** For air: 0.6847 $p_0A^*/\sqrt{RT_0}$.
  7. **did** Show that it is the most the area can pass · **tex**
     $\dfrac{\dot m\sqrt{RT_0}}{p_0A}=\sqrt\gamma\left(\dfrac{2}{\gamma+1}\right)^{(\gamma+1)/(2(\gamma-1))}\dfrac{A^*}{A}$ · **why** For any station $\dot m$ is
     the step-6 expression with $A^*=A\times(A^*/A)$; by D10 $A^*/A\le1$ with equality at M = 1. · **plain** The mass
     flux through a given area is largest when the flow there is sonic.
- **Result:** $\dfrac{T^*}{T_0}=\dfrac{2}{\gamma+1},\ \dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)},\ \dfrac{\rho^*}{\rho_0}=\left(\dfrac{2}{\gamma+1}\right)^{1/(\gamma-1)}$ (ours, from (15.28)–(15.30) at $M=1$) and $\dot m_{\max}=p_0A^*\sqrt{\dfrac{\gamma}{RT_0}}\left(\dfrac{2}{\gamma+1}\right)^{(\gamma+1)/(2(\gamma-1))}$ (ours).
- **Check:** units — Pa m² √(1/(J/kg)) = N √(s²/m²) = kg/s ✓. Number — $p_0$ = 500 kPa, $T_0$ = 400 K, $A^*$ = 10 cm²: 0.6847 ×
  500 000 × 0.001/√(287.05 × 400) = 1.010 kg/s. Code — `GD.choked_mass_flow`, `GD.mass_flux_function(1.0)` = 0.6847.
- **What it means:** a sonic throat is a flow meter that depends only on the supply — the idea of choking (C09). **Fails when:** the flow to the throat is not isentropic (losses lower $p_0$ at the throat).
- **Traps:** $\dot m\propto p_0/\sqrt{T_0}$, so a choked flow still rises with supply pressure (trap T10); the exponent is 3
  for γ = 1.4, not 3.5.


### D12 · The area–velocity relation: (15.32, corrected) → (15.33) → (15.34) — ★, 7 steps, in C06 (notebook · `area_mach_nozzle`)
- **Goal:** find which way the speed changes when the area of a duct changes. · **Start:** $u\,du+\dfrac{dp}{\rho}=0$ (15.26) (Euler along the duct),
  $c^2=(\partial p/\partial\rho)_s$ ($c^2\equiv(\partial p/\partial\rho)_s$ (1.19)) and $\dfrac{d\rho}{\rho}+\dfrac{du}{u}+\dfrac{dA}{A}=0$ (15.19) (mass). · **Plan:** • link dp to dρ through the sound speed • get the
  density change per unit speed change • put it into the mass balance. · **Tools:** differentials of algebraic laws as a
  small linear system (primer P362). · **Assumptions:** frictionless and adiabatic, hence isentropic (step 2); steady,
  quasi-one-dimensional.
- **Steps:**
  1. **did** Solve Euler's equation for the speed change · **tex** $u\,du=-\dfrac{dp}{\rho}$ · **why** Move $dp/\rho$ to the
     right of $u\,du+\dfrac{dp}{\rho}=0$ (15.26). · **plain** The speed rises where the pressure falls.
  2. **did** Link pressure and density changes · **tex** $dp=c^2\,d\rho$ · **why** The flow is isentropic, and by definition
     $c^2=(\partial p/\partial\rho)_s$: along the duct the pressure changes only through the density. · **plain** A
     pressure drop goes with a density drop; the factor is the sound speed squared.
  3. **did** Substitute · **tex** $u\,du=-\dfrac{dp}{\rho}=-c^2\dfrac{d\rho}{\rho}$ (15.32, third member with its minus sign restored) ·
     **why** Step 2 in step 1 — the minus sign has nowhere to go. (The book prints the third member without it; its own
     next line needs it.) · **plain** The speed rises where the density falls.
  4. **did** Divide by u² · **tex** $\dfrac{d\rho}{\rho}=-M^2\dfrac{du}{u}$ (15.33) · **why** $u\,du/u^2=du/u$ and $c^2/u^2=1/M^2$;
     then solve for $d\rho/\rho$. · **plain** A 1 % gain in speed costs M² per cent of the density.
  5. **did** Put this into the mass balance · **tex** $-M^2\dfrac{du}{u}+\dfrac{du}{u}+\dfrac{dA}{A}=0$ · **why** Substitute step 4 for
     the first term of $\dfrac{d\rho}{\rho}+\dfrac{du}{u}+\dfrac{dA}{A}=0$ (15.19). · **plain** One equation between the speed change and the area change.
  6. **did** Solve for the speed change · **tex** $\dfrac{du}{u}=-\dfrac{1}{1-M^2}\dfrac{dA}{A}$ (15.34) · **why** Collect
     $(1-M^2)\,du/u=-dA/A$ and divide — allowed away from M = 1. · **plain** Below Mach 1 a narrowing duct speeds the flow
     up; above Mach 1 a *widening* duct does.
  7. **did** Read the collected form at a throat · **tex** $(1-M^2)\,\dfrac{du}{u}=0\ \text{ where }\ dA=0$ · **why** The line before
     the division, with dA = 0. A product is zero if either factor is. · **plain** At a throat either the speed has a
     maximum or minimum (du = 0), or the flow is exactly sonic.
- **Result:** $u\,du=-dp/\rho=-c^2d\rho/\rho$ (15.32, third member with its minus sign restored); $d\rho/\rho=-M^2\,du/u$ (15.33); $\dfrac{du}{u}=-\dfrac{1}{1-M^2}\dfrac{dA}{A}$ (15.34).
- **Check:** units — all terms are pure numbers ✓. Limits — M → 0: $du/u=-dA/A$ (the garden hose, R15); M = 2: $du/u=+\tfrac13
  dA/A$. Number — dA/A = +1 % at M = 0.5: du/u = −1.33 %, dρ/ρ = +0.33 %, sum with dA/A: 0 ✓. Code —
  `GD.influence_coefficients(M, dA_over_A=0.01)`; with `printed_15_32=True` the factor becomes $1+M^2$ and never changes sign.
- **What it means:** the chapter's first surprise and the reason for the convergent–divergent nozzle: the only passage
  through Mach 1 is a throat. **Fails when:** there is friction or heat exchange (then they join dA on the right — C12).
- **Traps:** the book's missing minus sign in $u\,du=-dp/\rho=-c^2d\rho/\rho$ (slip #4); dA = 0 allows du = 0 *or* M = 1 — a throat need not be
  sonic; never divide by $1-M^2$ at M = 1.

### D13 · The pressure jump across a normal shock (15.39) — ★★★, 13 steps, in C07 (notebook · `normal_shock_jump`)
- **Goal:** from mass, momentum and energy alone, find the pressure behind a normal shock in terms of the upstream Mach
  number — the step the book covers with "after some algebra". · **Start:** $\rho_1u_1=\rho_2u_2$ (15.35), $p_1-p_2=-\rho_1u_1^2+\rho_2u_2^2$ (15.36), $h_1+\tfrac12u_1^2=h_2+\tfrac12u_2^2$ (15.37) — "what goes in comes out, for
  mass, momentum and energy". · **Plan:** • write the energy law with p and ρ • use one unknown, the velocity ratio •
  eliminate density and pressure • a quadratic appears whose root "nothing happens" we divide out • go back for the
  pressure. · **Tools:** dividing a known root out of a quadratic (primer P363); Vieta (P71 reminder); $h=C_pT$ and the
  gas law. · **Assumptions:** steady, constant area, no wall force, no heat (the three starting laws); perfect gas with
  constant γ (step 1).
- **Steps:**
  1. **did** Write the enthalpy with p and ρ · **tex** $\dfrac{\gamma}{\gamma-1}\dfrac{p_1}{\rho_1}+\dfrac12u_1^2=\dfrac{\gamma}{\gamma-1}\dfrac{p_2}{\rho_2}+\dfrac12u_2^2$ (15.38) · **why** $h=C_pT$, $C_p=\gamma R/(\gamma-1)$ and
     $T=p/(\rho R)$ give $h=\dfrac{\gamma}{\gamma-1}\dfrac p\rho$; put it into $h_1+\tfrac12u_1^2=h_2+\tfrac12u_2^2$ (15.37). · **plain** The energy law without
     temperature.
  2. **did** Name the velocity ratio and use mass · **tex** $r\equiv\dfrac{u_2}{u_1},\qquad\rho_2=\dfrac{\rho_1}{r}$ · **why**
     $\rho_1u_1=\rho_2u_2$ (15.35) gives $\rho_2=\rho_1u_1/u_2$. One unknown, r, will now carry the whole problem. · **plain** If the gas slows
     to a fraction r of its speed, its density rises by 1/r.
  3. **did** Solve momentum for the downstream pressure · **tex** $p_2=p_1+\rho_1u_1^2\,(1-r)$ · **why** $p_1-p_2=-\rho_1u_1^2+\rho_2u_2^2$ (15.36) rearranged:
     $p_2=p_1+\rho_1u_1^2-\rho_2u_2^2$, and $\rho_2u_2^2=(\rho_2u_2)u_2=\rho_1u_1\cdot ru_1$. · **plain** The pressure rises
     by the momentum flux the gas gives up.
  4. **did** Form the downstream p/ρ · **tex** $\dfrac{p_2}{\rho_2}=r\left[\dfrac{p_1}{\rho_1}+u_1^2(1-r)\right]$ · **why** Divide step 3
     by $\rho_2=\rho_1/r$. This is the combination the energy law needs. · **plain** Downstream "temperature" in terms of
     upstream quantities and r.
  5. **did** Substitute into the energy law · **tex**
     $\dfrac{\gamma}{\gamma-1}\dfrac{p_1}{\rho_1}+\dfrac12u_1^2=\dfrac{\gamma}{\gamma-1}\,r\left[\dfrac{p_1}{\rho_1}+u_1^2(1-r)\right]+\dfrac12r^2u_1^2$ · **why** Step
     4 and $u_2=ru_1$ in step 1. Only r is unknown now. · **plain** One equation, one unknown.
  6. **did** Divide by $u_1^2$ · **tex**
     $\dfrac{\gamma}{\gamma-1}a+\dfrac12=\dfrac{\gamma}{\gamma-1}\,r\,(a+1-r)+\dfrac12r^2,\qquad a\equiv\dfrac{p_1}{\rho_1u_1^2}$ · **why** $u_1\neq0$. The
     upstream state enters only through the number a (we shall see $a=1/\gamma M_1^2$). · **plain** The same equation
     in pure numbers.
  7. **did** Clear the fractions and collect powers of r · **tex**
     $(\gamma+1)\,r^2-2\gamma(a+1)\,r+2\gamma a+\gamma-1=0$ · **why** Multiply by 2(γ − 1) and move everything to one side:
     $r^2$ has the coefficient $2\gamma-(\gamma-1)=\gamma+1$. · **plain** A quadratic for the velocity ratio.
  8. **did** Test the root "nothing happens" · **tex**
     $(\gamma+1)-2\gamma(a+1)+2\gamma a+\gamma-1=0\quad\checkmark$ · **why** r = 1 (no change at all) obviously satisfies the
     three laws, so it must be a root; inserting it confirms the algebra. · **plain** "No shock" is always a solution.
  9. **did** Divide that root out · **tex** $(r-1)\left[(\gamma+1)\,r-(2\gamma a+\gamma-1)\right]=0$ · **why** A quadratic with a
     known root $r_0$ is $(r-r_0)$ times a linear factor; by Vieta the product of the roots is
     $(2\gamma a+\gamma-1)/(\gamma+1)$. · **plain** The quadratic splits into "nothing happens" and one other possibility.
  10. **did** Read off the other root · **tex** $r=\dfrac{\gamma-1}{\gamma+1}+\dfrac{2\gamma}{\gamma+1}\,a$ · **why** Set the square
      bracket to zero. This is the shock. · **plain** The speed ratio across the shock.
  11. **did** Form 1 − r · **tex** $1-r=\dfrac{2\,(1-\gamma a)}{\gamma+1}$ · **why** $1-\dfrac{\gamma-1}{\gamma+1}=\dfrac{2}{\gamma+1}$;
      step 3 needs this combination. · **plain** The fraction of its speed the gas loses.
  12. **did** Return to the momentum law · **tex**
      $\dfrac{p_2}{p_1}=1+\dfrac{1-r}{a}=1+\dfrac{2\gamma}{\gamma+1}\left[\dfrac{\rho_1u_1^2}{\gamma p_1}-1\right]$ (unnumbered, §15.6) · **why** Step 3
      divided by $p_1$ is $1+(1-r)/a$; insert step 11 and $1/a=\rho_1u_1^2/p_1$. · **plain** The pressure ratio in terms
      of the upstream momentum flux.
  13. **did** Recognise the Mach number · **tex** $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39) · **why**
      $c_1^2=\gamma p_1/\rho_1$, so $\rho_1u_1^2/(\gamma p_1)=u_1^2/c_1^2=M_1^2$ (and $a=1/\gamma M_1^2$). · **plain** The
      pressure jump depends on the upstream Mach number alone.
- **Result:** $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39) — and on the way $\dfrac{u_2}{u_1}=\dfrac{(\gamma-1)M_1^2+2}{(\gamma+1)M_1^2}$ (step 10 with $a=1/\gamma M_1^2$).
- **Check:** units — pure numbers ✓. Limits — $M_1=1$: $p_2/p_1=1$ (the two roots coincide: no shock); $M_1\to\infty$:
  $p_2/p_1\to2\gamma M_1^2/(\gamma+1)$. Number — $M_1$ = 2: 1 + 1.1667 × 3 = 4.5. Residuals of the three laws < 10⁻¹³
  (`GD.normal_shock_residuals`).
- **sympy check** (`check_src`; re-runs steps 3, 6–13; executed: prints `(r - 1)*(-2*a*gamma + gamma*r - gamma + r + 1) |
  (M_1**2*(gamma - 1) + 2)/(M_1**2*(gamma + 1))`):
  ```python
  import sympy as sp                                   # symbolic algebra: we re-run the derivation's own moves
  g, a, r = sp.symbols("gamma a r", positive=True)     # gamma; a = p1/(rho1 u1^2); r = u2/u1 (the velocity ratio)
  M1 = sp.symbols("M_1", positive=True)                # upstream Mach number
  p2_over_p1 = 1 + (1 - r) / a                         # step 3: momentum, p2 = p1 + rho1 u1^2 (1 - r), divided by p1
  energy = g/(g-1)*a + sp.Rational(1, 2) - (g/(g-1)*r*(a + 1 - r) + r**2/2)   # step 6: energy / u1^2, left - right
  quad = sp.expand(2*(g-1)*energy)                     # step 7: multiply by 2(gamma-1) -> a quadratic in r
  assert sp.simplify(quad - ((g+1)*r**2 - 2*g*(a+1)*r + 2*g*a + g - 1)) == 0   # ... exactly the printed quadratic
  assert sp.simplify(quad.subs(r, 1)) == 0             # step 8: r = 1 ("no shock") is always a root
  other = sp.cancel(quad / (r - 1))                    # step 9: divide the known root out -> a linear factor
  r_shock = sp.solve(other, r)[0]                      # step 10: the shock root
  assert sp.simplify(r_shock - ((g-1) + 2*g*a)/(g+1)) == 0
  jump = sp.simplify(p2_over_p1.subs(r, r_shock).subs(a, 1/(g*M1**2)))    # steps 11-13: back into momentum, a = 1/(gamma M1^2)
  assert sp.simplify(jump - (1 + 2*g/(g+1)*(M1**2 - 1))) == 0             # Eq. (15.39)
  print(sp.factor(quad), "|", sp.simplify(r_shock.subs(a, 1/(g*M1**2))))  # the factored quadratic and u2/u1(M1)
  ```
- **What it means:** three laws, one extra solution, one parameter: every normal shock in a given gas is labelled by
  $M_1$. Nothing here says $M_1>1$ — that comes from the second law (C08). **Fails when:** the gas is not perfect (very
  strong shocks dissociate air), or area, wall force or heat act across the layer.
- **Traps:** dividing by (r − 1) is legitimate only because we are looking for the *other* root; work with $u_2/u_1$, not
  with ρ and p separately; $\rho_1u_1^2/\gamma p_1$ is $M_1^2$ with the *upstream* sound speed; the energy law is not
  Bernoulli's equation (trap T17).

### D14 · The Mach number behind a normal shock (15.40) — ★★, 7 steps, in C07 (notebook)
- **Goal:** find the downstream Mach number from the upstream one. · **Start:** $p_1-p_2=-\rho_1u_1^2+\rho_2u_2^2$ (15.36) and $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39). · **Plan:** • write momentum
  with Mach numbers • insert the pressure jump • solve for $M_2^2$. · **Tools:** dynamic pressure in a gas (primer P359);
  a common denominator. · **Assumptions:** as D13.
- **Steps:**
  1. **did** Rearrange the momentum law · **tex** $p_1+\rho_1u_1^2=p_2+\rho_2u_2^2$ · **why** $p_1-p_2=-\rho_1u_1^2+\rho_2u_2^2$ (15.36) with the terms of each
     side together. · **plain** Pressure plus momentum flux is the same on both sides.
  2. **did** Use $\rho u^2=\gamma pM^2$ · **tex** $p_1\left(1+\gamma M_1^2\right)=p_2\left(1+\gamma M_2^2\right)$ (unnumbered, §15.6) · **why**
     $\rho u^2=\rho c^2M^2=\gamma pM^2$ on each side, with that side's own sound speed (primer P359). · **plain** The same
     balance in Mach numbers.
  3. **did** Isolate the downstream bracket · **tex** $1+\gamma M_2^2=\dfrac{1+\gamma M_1^2}{p_2/p_1}$ · **why** Divide by
     $p_2$. · **plain** The stronger the pressure jump, the smaller $M_2$.
  4. **did** Put the pressure jump on one denominator · **tex** $\dfrac{p_2}{p_1}=\dfrac{2\gamma M_1^2-(\gamma-1)}{\gamma+1}$ · **why**
     $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39): $1+\dfrac{2\gamma}{\gamma+1}(M_1^2-1)=\dfrac{\gamma+1+2\gamma M_1^2-2\gamma}{\gamma+1}$. · **plain** The same jump as one
     fraction.
  5. **did** Insert it · **tex** $1+\gamma M_2^2=\dfrac{(\gamma+1)\left(1+\gamma M_1^2\right)}{2\gamma M_1^2-(\gamma-1)}$ · **why** Step 4
     in step 3. · **plain** Only $M_1$ and γ remain.
  6. **did** Subtract 1 · **tex** $\gamma M_2^2=\dfrac{2\gamma+\gamma(\gamma-1)M_1^2}{2\gamma M_1^2-(\gamma-1)}$ · **why** Numerator:
     $(\gamma+1)+\gamma(\gamma+1)M_1^2-2\gamma M_1^2+(\gamma-1)=2\gamma+\gamma(\gamma-1)M_1^2$. · **plain** Nearly there.
  7. **did** Divide by γ · **tex** $M_2^2=\dfrac{(\gamma-1)M_1^2+2}{2\gamma M_1^2+1-\gamma}$ (15.40) · **why** γ ≠ 0; the denominator
     is rewritten in the book's order. · **plain** The downstream Mach number, from the upstream one.
- **Result:** $M_2^2=\dfrac{(\gamma-1)M_1^2+2}{2\gamma M_1^2+1-\gamma}$ (15.40).
- **Check:** units — numbers ✓. Limits — $M_1=1$: $M_2=1$; $M_1\to\infty$: $M_2^2\to(\gamma-1)/2\gamma$, i.e. 0.378 for air.
  Symmetry — exchanging $M_1$ and $M_2$ in step 2 leaves it unchanged, so the map is its own inverse: 2 → 0.5774 → 2.
  Number — $M_1$ = 2: 3.6/10.8 = 1/3.
- **What it means:** supersonic in, subsonic out — and the stronger the shock the slower the outflow, down to a floor.
  **Fails when:** as D13; and the formula has no real value for $M_1^2<(\gamma-1)/2\gamma$.
- **Traps:** solve for $M_2^2$, not $M_2$; the self-inverse property means a swapped input is not caught by the formula
  (trap T6).

### D15 · The density and temperature jumps (15.41), (15.42) — ★★, 7 steps, in C07 (notebook)
- **Goal:** complete the set of jump relations ("similarly obtained", says the book). · **Start:** $\rho_1u_1=\rho_2u_2$ (15.35), the velocity ratio
  of D13 step 10, $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39) and $p=\rho RT$. · **Plan:** • density from mass • temperature from the gas law • massage it
  into the printed form. · **Tools:** substitution; factoring a quadratic in $M_1^2$. · **Assumptions:** as D13.
- **Steps:**
  1. **did** Density ratio from mass · **tex** $\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac1r$ · **why** $\rho_1u_1=\rho_2u_2$ (15.35) divided by
     $\rho_1u_2$. · **plain** The gas is compressed exactly as much as it is slowed.
  2. **did** Write r with the Mach number · **tex** $r=\dfrac{(\gamma-1)M_1^2+2}{(\gamma+1)M_1^2}$ · **why** D13 step 10 with
     $a=1/\gamma M_1^2$: $\dfrac{\gamma-1}{\gamma+1}+\dfrac{2}{(\gamma+1)M_1^2}$. · **plain** The speed ratio in terms of $M_1$.
  3. **did** Invert · **tex** $\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}$ (15.41) · **why** Step 2 in
     step 1. · **plain** The density jump; it can never exceed (γ + 1)/(γ − 1).
  4. **did** Temperature from the gas law · **tex** $\dfrac{T_2}{T_1}=\dfrac{p_2}{p_1}\,\dfrac{\rho_1}{\rho_2}$ · **why** $T=p/(\rho R)$
     on both sides; R cancels. · **plain** Temperature ratio = pressure ratio over density ratio.
  5. **did** Insert both jumps · **tex**
     $\dfrac{T_2}{T_1}=\dfrac{\left[2\gamma M_1^2-(\gamma-1)\right]\left[(\gamma-1)M_1^2+2\right]}{(\gamma+1)^2M_1^2}$ · **why** D14 step 4 and step 2
     above. This is already correct; the remaining steps only reach the book's form. · **plain** The temperature jump as
     one fraction.
  6. **did** Separate a 1 · **tex**
     $\dfrac{T_2}{T_1}=1+\dfrac{2(\gamma-1)\left[\gamma M_1^4-(\gamma-1)M_1^2-1\right]}{(\gamma+1)^2M_1^2}$ · **why** Expand the numerator,
     $2\gamma(\gamma-1)M_1^4+\left[4\gamma-(\gamma-1)^2\right]M_1^2-2(\gamma-1)$, and subtract $(\gamma+1)^2M_1^2=\left[4\gamma+(\gamma-1)^2\right]M_1^2$.
     · **plain** One, plus a term that vanishes when there is no shock.
  7. **did** Factor the bracket · **tex**
     $\dfrac{T_2}{T_1}=1+\dfrac{2(\gamma-1)}{(\gamma+1)^2}\dfrac{\gamma M_1^2+1}{M_1^2}\left(M_1^2-1\right)$ (15.42) · **why**
     $\gamma M_1^4-(\gamma-1)M_1^2-1=(\gamma M_1^2+1)(M_1^2-1)$ — multiply out to check. · **plain** The temperature rises
     by an amount proportional to $M_1^2-1$.
- **Result:** $\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}$ (15.41) and $\dfrac{T_2}{T_1}=1+\dfrac{2(\gamma-1)}{(\gamma+1)^2}\dfrac{\gamma M_1^2+1}{M_1^2}\left(M_1^2-1\right)$ (15.42).
- **Check:** units — numbers ✓. Limits — $M_1=1$: both 1; $M_1\to\infty$: density → 6, temperature → ∞ like $M_1^2$ (air).
  Number — $M_1$ = 2: 9.6/3.6 = 2.667; 4.5/2.667 = 1.6875 = 1 + (0.8/5.76)(6.6/4)(3) ✓.
- **What it means:** a strong shock heats without limit but compresses only sixfold. **Fails when:** as D13.
- **Traps:** check a rearrangement like step 6–7 with numbers before trusting it; the density limit is a property of the
  gas (γ), not of the shock.

### D16 · The entropy jump across a shock (15.43) — ★, 5 steps, in C08 (notebook)
- **Goal:** compute the entropy made by a shock, to let the second law judge it. · **Start:** $S_2-S_1=C_p\ln\dfrac{T_2}{T_1}-R\ln\dfrac{p_2}{p_1}=C_v\ln\dfrac{T_2}{T_1}-R\ln\dfrac{\rho_2}{\rho_1}$ (15.1, entropy line) in its second form,
  $S_2-S_1=C_v\ln\dfrac{T_2}{T_1}-R\ln\dfrac{\rho_2}{\rho_1}$. · **Plan:** • replace T by p/ρ • combine the logarithms • insert the
  jumps. · **Tools:** logarithm rules; $R=(\gamma-1)C_v$. · **Assumptions:** perfect gas with constant specific heats.
- **Steps:**
  1. **did** Replace the temperature ratio · **tex** $\dfrac{T_2}{T_1}=\dfrac{p_2}{p_1}\,\dfrac{\rho_1}{\rho_2}$ · **why** The gas law on
     both sides. We want an expression in the two jumps we know. · **plain** Temperature ratio = pressure ratio over
     density ratio.
  2. **did** Substitute and split the logarithm · **tex**
     $S_2-S_1=C_v\ln\dfrac{p_2}{p_1}+C_v\ln\dfrac{\rho_1}{\rho_2}+R\ln\dfrac{\rho_1}{\rho_2}$ · **why** ln(ab) = ln a + ln b; and
     $-\ln(\rho_2/\rho_1)=+\ln(\rho_1/\rho_2)$. · **plain** Three logarithms, two of the same ratio.
  3. **did** Add the two density terms · **tex** $S_2-S_1=C_v\ln\dfrac{p_2}{p_1}+\gamma C_v\ln\dfrac{\rho_1}{\rho_2}$ · **why**
     $C_v+R=C_p=\gamma C_v$. · **plain** The density ratio carries the weight γ.
  4. **did** Combine and divide by $C_v$ · **tex** $\dfrac{S_2-S_1}{C_v}=\ln\left[\dfrac{p_2}{p_1}\left(\dfrac{\rho_1}{\rho_2}\right)^{\gamma}\right]$ (15.43) ·
     **why** $\gamma\ln x=\ln x^\gamma$. · **plain** Entropy measures how far $p/\rho^\gamma$ has moved; it is constant
     only in isentropic change.
  5. **did** Insert the two jumps · **tex**
     $\dfrac{S_2-S_1}{C_v}=\ln\left\{\left[1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)\right]\left[\dfrac{(\gamma-1)M_1^2+2}{(\gamma+1)M_1^2}\right]^{\gamma}\right\}$ ·
     **why** $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39) and the inverse of $\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}$ (15.41). · **plain** The entropy jump as a function of $M_1$ alone.
- **Result:** $\dfrac{S_2-S_1}{C_v}=\ln\left[\dfrac{p_2}{p_1}\left(\dfrac{\rho_1}{\rho_2}\right)^{\gamma}\right]$ (15.43) with its second member.
- **Check:** units — a number ✓. Limit — $M_1=1$: ln 1 = 0. Number — $M_1$ = 2: ln(4.5 × 0.375^1.4) = ln 1.140 = 0.131. Sign —
  $M_1$ = 0.7: −0.030 (forbidden).
- **What it means:** positive only for $M_1>1$ (D17 shows why near 1; the figure shows it everywhere): shocks compress.
  **Fails when:** the gas has another equation of state (N90).
- **Traps:** the exponent γ sits on the *inverse* density ratio; S is per unit mass; dividing by $C_v$, not R.

### D17 · The entropy jump of a weak shock (15.44a), (15.44b) — ★★★, 11 steps, in C08 (notebook · `normal_shock_jump`) — the book gives the result only
- **Goal:** show that a weak shock makes almost no entropy — only at third order in its strength — and that the sign of
  the entropy change is the sign of $M_1^2-1$. · **Start:** $\dfrac{S_2-S_1}{C_v}=\ln\left[\dfrac{p_2}{p_1}\left(\dfrac{\rho_1}{\rho_2}\right)^{\gamma}\right]$ (15.43) in terms of $M_1$ (D16 step 5). · **Plan:** • call the
  small number ε = $M_1^2-1$ • write both brackets as 1 + small • expand the logarithms to third order • watch the
  first two orders cancel • convert ε to the pressure jump. · **Tools:** the series of ln(1 + x) to third order (primer
  P364); orders of smallness (P68 reminder); sympy `series` (P117 reminder). · **Assumptions:** $\lvert\varepsilon\rvert\ll1$
  (steps 5–9).
- **Steps:**
  1. **did** Name the small quantity · **tex** $\varepsilon\equiv M_1^2-1$ · **why** Near Mach 1 this is small, and both
     brackets of the entropy formula become "one plus something small". · **plain** ε measures how far above (or below)
     Mach 1 the upstream flow is.
  2. **did** Rewrite the pressure bracket · **tex** $\dfrac{p_2}{p_1}=1+a\,\varepsilon,\qquad a\equiv\dfrac{2\gamma}{\gamma+1}$ · **why**
     $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39) is already of this form. · **plain** The pressure jump is proportional to ε.
  3. **did** Rewrite the density bracket · **tex**
     $\dfrac{\rho_1}{\rho_2}=\dfrac{(\gamma-1)(1+\varepsilon)+2}{(\gamma+1)(1+\varepsilon)}=\dfrac{1+b\,\varepsilon}{1+\varepsilon},\qquad b\equiv\dfrac{\gamma-1}{\gamma+1}$ · **why**
     Put $M_1^2=1+\varepsilon$ in the inverse of $\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}$ (15.41) and divide top and bottom by γ + 1. · **plain** The inverse density
     ratio as a quotient of two "one plus small" factors.
  4. **did** Take the logarithm of the product · **tex**
     $\dfrac{S_2-S_1}{C_v}=\ln(1+a\varepsilon)+\gamma\ln(1+b\varepsilon)-\gamma\ln(1+\varepsilon)$ · **why** ln of a product and of a power;
     exact so far. · **plain** Three logarithms of "one plus small".
  5. **did** Expand each to third order · **tex**
     $\ln(1+k\varepsilon)=k\varepsilon-\tfrac12k^2\varepsilon^2+\tfrac13k^3\varepsilon^3-\dots\quad(k=a,\,b,\,1)$ · **why** The series of primer P364,
     valid for $\lvert k\varepsilon\rvert<1$. We keep three terms because the first two will cancel. · **plain** Each
     logarithm as a short polynomial in ε.
  6. **did** Collect the first-order terms · **tex**
     $a+\gamma b-\gamma=\dfrac{2\gamma+\gamma(\gamma-1)-\gamma(\gamma+1)}{\gamma+1}=0$ · **why** Common denominator γ + 1; the numerator is
     $2\gamma+\gamma^2-\gamma-\gamma^2-\gamma$. · **plain** To first order a weak shock makes no entropy.
  7. **did** Collect the second-order terms · **tex**
     $a^2+\gamma b^2-\gamma=\dfrac{4\gamma^2+\gamma(\gamma-1)^2-\gamma(\gamma+1)^2}{(\gamma+1)^2}=0$ · **why** $(\gamma-1)^2-(\gamma+1)^2=-4\gamma$, so
     the numerator is $4\gamma^2-4\gamma^2$. · **plain** Nor to second order.
  8. **did** Collect the third-order terms · **tex**
     $\dfrac13\left(a^3+\gamma b^3-\gamma\right)=\dfrac{8\gamma^3+\gamma(\gamma-1)^3-\gamma(\gamma+1)^3}{3(\gamma+1)^3}=\dfrac{2\gamma(\gamma-1)}{3(\gamma+1)^2}$ · **why**
     $(\gamma-1)^3-(\gamma+1)^3=-(6\gamma^2+2)$, so the numerator is $2\gamma^3-2\gamma=2\gamma(\gamma-1)(\gamma+1)$. · **plain** The first
     term that survives is the cubic one — and its coefficient is positive.
  9. **did** Write the result in ε · **tex** $\dfrac{S_2-S_1}{C_v}\cong\dfrac{2\gamma(\gamma-1)}{3(\gamma+1)^2}\left(M_1^2-1\right)^3$ (15.44a) ·
     **why** Steps 6–8 in step 4; ≅ because fourth and higher orders are dropped. · **plain** The entropy jump is
     proportional to the cube of $M_1^2-1$: positive above Mach 1, negative below.
  10. **did** Express ε by the pressure jump · **tex** $\varepsilon=\dfrac{\gamma+1}{2\gamma}\,\dfrac{p_2-p_1}{p_1}$ · **why** Step 2
      solved for ε — exact. · **plain** Shock strength measured by its pressure rise.
  11. **did** Substitute · **tex** $\dfrac{S_2-S_1}{C_v}\cong\dfrac{\gamma^2-1}{12\gamma^2}\left(\dfrac{p_2-p_1}{p_1}\right)^3$ (15.44b) · **why**
      $\dfrac{2\gamma(\gamma-1)}{3(\gamma+1)^2}\cdot\dfrac{(\gamma+1)^3}{8\gamma^3}=\dfrac{(\gamma-1)(\gamma+1)}{12\gamma^2}$. · **plain** A 10 % pressure
      jump makes an entropy change of a few times 10⁻⁵ $C_v$.
- **Result:** $\dfrac{S_2-S_1}{C_v}\cong\dfrac{2\gamma(\gamma-1)}{3(\gamma+1)^2}\left(M_1^2-1\right)^3$ (15.44a) and $\dfrac{S_2-S_1}{C_v}\cong\dfrac{\gamma^2-1}{12\gamma^2}\left(\dfrac{p_2-p_1}{p_1}\right)^3$ (15.44b).
- **Check:** units — numbers ✓. Sign — odd in ε: $M_1<1$ would give ΔS < 0. Number — γ = 1.4: coefficient 1.12/17.28 =
  0.0648; at $M_1$ = 1.1 (ε = 0.21): 6.0 × 10⁻⁴ against the exact 4.2 × 10⁻⁴ — the next term, −2γ²(γ − 1)ε⁴/(γ + 1)³, is not
  small yet; at $M_1$ = 1.05: 7.0 × 10⁻⁵ against 5.9 × 10⁻⁵.
- **sympy check** (`check_src`; re-runs steps 2–3 and 6–11, then lets sympy expand $\dfrac{S_2-S_1}{C_v}=\ln\left[\dfrac{p_2}{p_1}\left(\dfrac{\rho_1}{\rho_2}\right)^{\gamma}\right]$ independently; executed: prints
  `2*gamma*(gamma - 1)/(3*(gamma + 1)**2) | -2*gamma**2*(gamma - 1)/(gamma + 1)**3`):
  ```python
  import sympy as sp                                   # symbolic algebra
  g, e = sp.symbols("gamma epsilon", positive=True)    # gamma and the small number epsilon = M1^2 - 1
  a = 2*g/(g+1); b = (g-1)/(g+1)                       # steps 2-3: the two coefficients
  M1sq = 1 + e                                         # M1^2 written with epsilon
  br1 = 1 + 2*g/(g+1)*(M1sq - 1)                       # first bracket of (15.43): p2/p1
  br2 = ((g-1)*M1sq + 2)/((g+1)*M1sq)                  # second bracket: rho1/rho2
  assert sp.simplify(br1 - (1 + a*e)) == 0             # step 2: first bracket = 1 + a*eps
  assert sp.simplify(br2 - (1 + b*e)/(1 + e)) == 0     # step 3: second bracket = (1 + b*eps)/(1 + eps)
  c1 = sp.simplify(a + g*b - g)                        # step 6: coefficient of eps
  c2 = sp.simplify(a**2 + g*b**2 - g)                  # step 7: coefficient of -eps^2/2
  c3 = sp.factor(sp.simplify((a**3 + g*b**3 - g)/3))   # step 8: coefficient of eps^3
  assert c1 == 0 and c2 == 0                           # first and second order cancel identically
  assert sp.simplify(c3 - 2*g*(g-1)/(3*(g+1)**2)) == 0 # Eq. (15.44a)
  S = sp.log(br1) + g*sp.log(br2)                      # the exact entropy jump (15.43) / Cv
  ser = sp.series(S, e, 0, 5).removeO()                # independent route: let sympy expand the logarithm
  assert sp.simplify(ser.coeff(e, 1)) == 0 and sp.simplify(ser.coeff(e, 2)) == 0
  assert sp.simplify(ser.coeff(e, 3) - c3) == 0
  x = sp.symbols("x", positive=True)                   # x = (p2 - p1)/p1 = a*eps  (steps 10-11)
  assert sp.simplify(c3*(x/a)**3 - (g**2 - 1)/(12*g**2)*x**3) == 0        # Eq. (15.44b)
  print(c3, "|", sp.factor(sp.simplify(ser.coeff(e, 4))))                # the cubic coefficient and the next term
  ```
- **What it means:** weak compressions are nearly reversible — which is why loud sound is still acoustics (C01), why
  the weak-wave rule holds for expansions too (N103), and why thin-airfoil theory may ignore losses (C16). And the odd
  power makes the second law's verdict visible: only $M_1>1$. **Fails when:** the shock is not weak — at $M_1$ = 1.5 the
  cubic is four times too large.
- **Traps:** stopping at second order gives exactly zero; dropping any third-order piece gives a wrong coefficient; both
  brackets must be written as 1 + small *before* expanding; the cubic is a local statement.

### D18 · The stagnation-pressure ratio across a shock — ★★, 7 steps, in C08 (notebook) — ours (the book has only a table column)
- **Goal:** find what a shock costs: the fraction of the stagnation pressure that survives it. · **Start:** $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29) on each
  side of the shock, $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39), $\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}$ (15.41). · **Plan:** • chain three pressure ratios • use $T_{02}=T_{01}$ to remove $M_2$ • insert
  the jumps • compare with the entropy jump • read off the sonic area. · **Tools:** a chain of ratios; D09; D11. ·
  **Assumptions:** adiabatic shock, perfect gas.
- **Steps:**
  1. **did** Chain the ratios · **tex** $\dfrac{p_{02}}{p_{01}}=\dfrac{p_{02}}{p_2}\,\dfrac{p_2}{p_1}\,\dfrac{p_1}{p_{01}}$ · **why**
     Multiply and divide by the static pressures. · **plain** Stagnation over static behind, static jump, static over
     stagnation ahead.
  2. **did** Use the isentropic relation on each side · **tex** $\dfrac{p_{02}}{p_{01}}=\dfrac{p_2}{p_1}\left(\dfrac{T_{02}/T_2}{T_{01}/T_1}\right)^{\gamma/(\gamma-1)}$ ·
     **why** $p_0/p=(T_0/T)^{\gamma/(\gamma-1)}$ holds at each station separately (D08 step 4). · **plain** The static
     jump, corrected by the two stagnation ratios.
  3. **did** Use the energy law · **tex** $\dfrac{p_{02}}{p_{01}}=\dfrac{p_2}{p_1}\left(\dfrac{T_1}{T_2}\right)^{\gamma/(\gamma-1)}$ · **why**
     $T_{02}=T_{01}$ across a shock ($h_1+\tfrac12u_1^2=h_2+\tfrac12u_2^2$ (15.37)), so the stagnation temperatures cancel. No $M_2$ is needed. · **plain** Only the
     static pressure and temperature jumps matter.
  4. **did** Replace the temperature ratio · **tex**
     $\dfrac{p_{02}}{p_{01}}=\left(\dfrac{\rho_2}{\rho_1}\right)^{\gamma/(\gamma-1)}\left(\dfrac{p_1}{p_2}\right)^{1/(\gamma-1)}$ · **why**
     $T_1/T_2=(p_1/p_2)(\rho_2/\rho_1)$; the pressure exponents combine as $1-\dfrac{\gamma}{\gamma-1}=-\dfrac1{\gamma-1}$. · **plain** A
     power of the density jump times a power of the inverse pressure jump.
  5. **did** Insert the jumps · **tex** $\dfrac{p_{02}}{p_{01}}=\left[\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}\right]^{\gamma/(\gamma-1)}\left[\dfrac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^{1/(\gamma-1)}$ (ours) · **why** $\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}$ (15.41), and $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39) on one denominator:
     $p_1/p_2=(\gamma+1)/(2\gamma M_1^2-(\gamma-1))$. · **plain** The surviving fraction as a function of $M_1$.
  6. **did** Compare with the entropy jump · **tex** $\dfrac{p_{02}}{p_{01}}=e^{-(S_2-S_1)/R}$ · **why** The logarithm of step 4
     is $-\dfrac1{\gamma-1}\ln\left[\dfrac{p_2}{p_1}\left(\dfrac{\rho_1}{\rho_2}\right)^\gamma\right]=-\dfrac{S_2-S_1}{(\gamma-1)C_v}$ by $\dfrac{S_2-S_1}{C_v}=\ln\left[\dfrac{p_2}{p_1}\left(\dfrac{\rho_1}{\rho_2}\right)^{\gamma}\right]$ (15.43), and
     $(\gamma-1)C_v=R$. The same as D09 step 5. · **plain** The loss of stagnation pressure *is* the entropy made.
  7. **did** Read off the sonic area · **tex** $p_{01}A_1^*=p_{02}A_2^*\ \Rightarrow\ \dfrac{A_1^*}{A_2^*}=\dfrac{p_{02}}{p_{01}}$ · **why**
     By $\dot m_{\max}=p_0A^*\sqrt{\dfrac{\gamma}{RT_0}}\left(\dfrac{2}{\gamma+1}\right)^{(\gamma+1)/(2(\gamma-1))}$ (ours) the mass flow is a constant times $p_0A^*/\sqrt{T_0}$; mass flow and $T_0$ are the same on both sides. ·
     **plain** Behind a shock the flow needs a larger throat to pass the same mass.
- **Result:** $\dfrac{p_{02}}{p_{01}}=\left[\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}\right]^{\gamma/(\gamma-1)}\left[\dfrac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^{1/(\gamma-1)}$ (ours) $=e^{-(S_2-S_1)/R}=A_1^*/A_2^*$.
- **Check:** units — numbers ✓. Limit — $M_1=1$: 1. Number — $M_1$ = 2: 2.667^3.5 × (2.4/10.8)^2.5 = 30.97 × 0.02328 = 0.7209 =
  exp(−0.1309/0.4) (the public NACA Report 1135 row). Three routes agree in the C08 from-scratch cell.
- **What it means:** energy is conserved across a shock; what is lost is the pressure the flow could recover. In a
  nozzle it shows up as a larger sonic area downstream (C10). **Fails when:** the frame changes — stagnation quantities
  are frame-dependent (trap T7).
- **Traps:** the middle ratio of step 1 is evaluated at $M_2$ if you do not use step 3; $A^*$ grows across a shock (trap
  T8); the ratio is downstream over upstream.

### D19 · Mass flow through a convergent nozzle against back pressure — ★★, 7 steps, in C09 (notebook · `nozzle_back_pressure`) — ours (§15.7 has no equation)
- **Goal:** compute how much gas a convergent nozzle passes as the back pressure is lowered — and show where it stops
  rising. · **Start:** $\rho uA=\dot m=\text{const.}$ (15.18) at the exit plane, with $\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28), $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29) and the range of influence (N20). · **Plan:** • write the
  exit mass flux with the exit Mach number • get that Mach number from the back pressure while the exit is subsonic •
  find where it reaches 1 • argue what happens below. · **Tools:** D08, D11; piecewise functions (P46 reminder). ·
  **Assumptions:** isentropic from the reservoir to the exit; quasi-one-dimensional.
- **Steps:**
  1. **did** Write the exit mass flow · **tex** $\dot m=\rho_eu_eA_e$ · **why** $\rho uA=\dot m=\text{const.}$ (15.18) at the exit plane. · **plain** Mass
     flow = exit density × exit speed × exit area.
  2. **did** Express ρu with pressure, Mach number and temperature · **tex** $\rho u=\dfrac{p}{RT}\,M\sqrt{\gamma RT}=pM\sqrt{\dfrac{\gamma}{RT}}$
     · **why** Gas law and $u=Mc$ with $c=\sqrt{\gamma RT}$. · **plain** Mass flux is pressure times Mach number over the
     square root of temperature.
  3. **did** Insert the stagnation ratios · **tex**
     $\dfrac{\dot m\sqrt{RT_0}}{p_0A_e}=\sqrt\gamma\,M_e\left(1+\dfrac{\gamma-1}{2}M_e^2\right)^{-(\gamma+1)/(2(\gamma-1))}$ (ours) · **why** $p=p_0\Theta^{-\gamma/(\gamma-1)}$
     and $T=T_0/\Theta$: the exponent is $-\dfrac{\gamma}{\gamma-1}+\dfrac12$. · **plain** A dimensionless mass flow that depends
     only on the exit Mach number.
  4. **did** Tie the exit Mach number to the back pressure · **tex**
     $M_e^2=\dfrac{2}{\gamma-1}\left[\left(\dfrac{p_B}{p_0}\right)^{-(\gamma-1)/\gamma}-1\right]\quad(\text{exit subsonic})$ · **why** A subsonic exit can
     hear the tank, so its pressure equals $p_B$; then invert $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29). · **plain** Lower back pressure, faster exit.
  5. **did** Find where the exit turns sonic · **tex** $\dfrac{p_B}{p_0}=\dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)}$ (ours) ·
     **why** Put $M_e=1$ in $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29): the sonic ratio of D11. · **plain** For air, when the back pressure is 52.8 % of the
     supply.
  6. **did** Below that, freeze the exit · **tex** $M_e=1,\quad p_e=p^*>p_B,\quad\dot m=\dot m_{\max}$ · **why** A converging
     duct cannot take the flow beyond Mach 1 (D12), and with the exit sonic no pressure signal travels upstream (N20):
     nothing inside can change. · **plain** The nozzle is choked: the exit pressure stays at $p^*$ and the mass flow at
     its maximum.
  7. **did** Check how the curve meets its plateau · **tex** $\left.\dfrac{d\dot m}{dp_B}\right\rvert_{M_e=1}=0$ · **why** By D11 step 7 the
     flux function has its maximum at M = 1, so its slope with respect to $M_e$ — and hence $p_B$ — vanishes there. ·
     **plain** The mass-flow curve flattens smoothly into the plateau; there is no corner.
- **Result:** $\dfrac{\dot m\sqrt{RT_0}}{p_0A}=\sqrt\gamma\,M\left(1+\dfrac{\gamma-1}{2}M^2\right)^{-(\gamma+1)/(2(\gamma-1))}$ (ours) at the exit with $M_e$ from step 4 for $p_B/p_0\ge$ $\dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)}$ (ours), and $\dot m_{\max}=p_0A^*\sqrt{\dfrac{\gamma}{RT_0}}\left(\dfrac{2}{\gamma+1}\right)^{(\gamma+1)/(2(\gamma-1))}$ (ours) below.
- **Check:** units — the left side of step 3 is a number ✓. Limits — $p_B=p_0$: no flow; $p_B\to0$: still $\dot m_{\max}$.
  Number — $p_B/p_0$ = 0.8: $M_e$ = 0.574, flux 0.561 (82 % of 0.685). Code — `GD.mass_flow_vs_back_pressure`.
- **What it means:** the back pressure controls the flow only while news of it can reach the throat. **Fails when:** the
  nozzle has a diverging part (C10) or losses.
- **Traps:** the exit pressure equals $p_B$ only while the exit is subsonic; "choked" is a ceiling the supply can still
  raise (trap T10).

### D20 · The Laval nozzle: three dividing back pressures and the place of the shock — ★★, 9 steps, in C10 (notebook · `nozzle_back_pressure`) — ours
- **Goal:** for a convergent–divergent nozzle of given exit-to-throat area ratio, find the back pressures that separate
  its regimes, and where a shock stands for a back pressure in between. · **Start:** $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31) (two roots), $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29), $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39), $\dfrac{p_{02}}{p_{01}}=\left[\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}\right]^{\gamma/(\gamma-1)}\left[\dfrac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^{1/(\gamma-1)}$ (ours). ·
  **Plan:** • two isentropic solutions give $p_b$ and $p_f$ • a shock at the exit gives $p_d$ • for a shock inside:
  isentropic before, jump, isentropic after with a new sonic area, exit pressure = back pressure. · **Tools:** a function
  with two inverses (primer P360); nested root finding (primer P365); D10, D13, D18. · **Assumptions:** isentropic except
  across the shock; quasi-one-dimensional; a normal shock.
- **Steps:**
  1. **did** Solve the area relation at the exit · **tex**
     $\dfrac{A_e}{A_t}=\left.\dfrac{A}{A^*}\right\rvert_{M}\ \Rightarrow\ M=M_{\text{sub}}\ \text{or}\ M_{\text{sup}}$ · **why** With a sonic throat $A^*=A_t$;
     $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31) then has one root on each side of 1 (primer P360). · **plain** The exit can be subsonic or supersonic for the
     same shape.
  2. **did** The subsonic root gives the first dividing pressure · **tex** $\dfrac{p_b}{p_0}=\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sub}}}$ ·
     **why** $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29) at $M_{\text{sub}}$: sonic at the throat only, subsonic everywhere else, exit pressure = back pressure. ·
     **plain** Above $p_b$ the nozzle is a venturi; at $p_b$ it has just choked.
  3. **did** The supersonic root gives the design pressure · **tex** $\dfrac{p_f}{p_0}=\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sup}}}$ ·
     **why** $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29) at $M_{\text{sup}}$: shock-free supersonic flow to the exit. · **plain** The one back pressure the
     nozzle was shaped for.
  4. **did** Put a normal shock in the exit plane · **tex** $\dfrac{p_d}{p_0}=\dfrac{p_f}{p_0}\left.\dfrac{p_2}{p_1}\right\rvert_{M_{\text{sup}}}$ ·
     **why** Just ahead of the exit the pressure is $p_f$; $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39) raises it across the shock. Since the jump exceeds 1,
     $p_f<p_d<p_b$. · **plain** Between $p_d$ and $p_b$ the shock must stand inside.
  5. **did** For a shock inside, the flow ahead of it · **tex** $\dfrac{A_{\text{shock}}}{A_t}=\left.\dfrac{A}{A^*}\right\rvert_{M_1}$ · **why**
     From the throat to the shock the flow is isentropic and supersonic with $A^*=A_t$. · **plain** The Mach number just
     ahead of the shock fixes where it stands (and the other way round).
  6. **did** Across the shock · **tex** $p_{02}=p_{01}\,\dfrac{p_{02}}{p_{01}}(M_1),\qquad A_2^*=A_t\,\dfrac{p_{01}}{p_{02}}$ · **why** D18:
     the shock lowers the stagnation pressure and enlarges the sonic area by the same factor. · **plain** Behind the shock
     the "throat the flow would need" is bigger than the real one.
  7. **did** Behind the shock, to the exit · **tex**
     $M_e=M_{\text{sub}}\!\left(\dfrac{A_e}{A_2^*}\right)=M_{\text{sub}}\!\left(\dfrac{A_e}{A_t}\dfrac{p_{02}}{p_{01}}\right)$ · **why** The flow behind a
     normal shock is subsonic and, in a widening duct, stays so (D12): the subsonic branch with the *new* sonic area. ·
     **plain** The diverging part now works as a diffuser.
  8. **did** Match the exit to the tank · **tex**
     $\dfrac{p_B}{p_0}=\dfrac{p_{02}}{p_{01}}(M_1)\left.\dfrac{p}{p_0}\right\rvert_{M_e}$ (ours) · **why** The exit is subsonic, so its pressure is $p_B$;
     $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29) with the *downstream* stagnation pressure, and $p_0\equiv p_{01}$. One equation for $M_1$. · **plain** Only one
     shock strength lets the outflow end exactly at the back pressure.
  9. **did** Solve for $M_1$, then place the shock · **tex** $M_1\in(1,M_{\text{sup}}):\ \text{step 8}\ \Rightarrow\ \text{step 5}$ ·
     **why** The right side of step 8 falls steadily from $p_b/p_0$ at $M_1=1$ to $p_d/p_0$ at $M_1=M_{\text{sup}}$, so
     a bracketed root-finder finds the one root (primer P365). · **plain** Lower back pressure, stronger shock, further
     down the bell.
- **Result:** $\dfrac{p_b}{p_0}=\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sub}}(A_e/A_t)},\ \dfrac{p_f}{p_0}=\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sup}}(A_e/A_t)},\ \dfrac{p_d}{p_0}=\dfrac{p_f}{p_0}\left.\dfrac{p_2}{p_1}\right\rvert_{M_{\text{sup}}}$ (ours) and, for $p_d<p_B<p_b$, $\dfrac{p_B}{p_0}=\dfrac{p_{02}}{p_{01}}(M_1)\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sub}}\left(\frac{A_e}{A_t}\frac{p_{02}}{p_{01}}\right)}$ (ours) with $A_{\text{shock}}/A_t$ from $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31) at $M_1$.
- **Check:** units — ratios ✓. Order — $p_f<p_d<p_b$ for every area ratio above 1. Limits — $M_1\to1$: no shock, $p_B\to
  p_b$; $M_1\to M_{\text{sup}}$: $p_B\to p_d$. Number — $A_e/A_t$ = 2: 0.0939 < 0.5134 < 0.9372; $p_B/p_0$ = 0.7 → $M_1$ = 1.863 at
  $A/A_t$ = 1.510. Mass flow — equal to $\dot m_{\max}=p_0A^*\sqrt{\dfrac{\gamma}{RT_0}}\left(\dfrac{2}{\gamma+1}\right)^{(\gamma+1)/(2(\gamma-1))}$ for every $p_B\le p_b$.
- **What it means:** below $p_b$ the back pressure no longer sets the flow rate; it sets the *losses*. Below $p_d$ it
  does nothing inside at all. **Fails when:** the boundary layer separates behind the shock (real nozzles: the shock is
  not a clean normal one), or the jet outside matters.
- **Traps:** behind the shock use the subsonic branch and the new $A^*$ (trap T8); the ratios $p/p_0$ behind the shock
  are to the *downstream* stagnation pressure; the order of the three pressures.

### D21 · The outlet Mach number of a duct with friction and heat (15.46) — ★★★, 13 steps, in C11 (notebook · `duct_friction_heating`)
- **Goal:** find the Mach number at the outlet of a constant-area duct from the inlet Mach number, a lumped friction f
  and a lumped heat q. · **Start:** $\rho_1u_1=\rho_2u_2,\ \ p_1+\rho_1u_1^2=p_2+\rho_2u_2^2+p_1f,\ \ h_1+\tfrac12u_1^2+h_1q=h_2+\tfrac12u_2^2$ (15.45, with $f=F/(p_1A)$, $q=Q/h_1$) — "mass, momentum with a retarding wall force, energy with heat added". ·
  **Plan:** • write momentum and energy with Mach numbers • write mass with p, M and T • eliminate pressure and
  temperature • everything about the inlet collects into one number A • solve a quadratic for $M_2^2$ • check it
  against the normal shock. · **Tools:** dynamic pressure in a gas (primer P359); a quadratic in M squared and its
  discriminant (primer P366); the quadratic formula (P159 reminder). · **Assumptions:** constant area; steady; perfect gas;
  f and q are totals between the stations.
- **Steps:**
  1. **did** Momentum with Mach numbers · **tex** $p_1\left(1+\gamma M_1^2-f\right)=p_2\left(1+\gamma M_2^2\right)$ · **why**
     $\rho u^2=\gamma pM^2$ on each side (primer P359); the friction term $p_1f$ joins the inlet side. · **plain** The
     momentum balance in pressures and Mach numbers.
  2. **did** Energy with Mach numbers · **tex**
     $h_1\left(1+\tfrac{\gamma-1}{2}M_1^2+q\right)=h_2\left(1+\tfrac{\gamma-1}{2}M_2^2\right)$ · **why** $\tfrac12u^2=\tfrac12M^2\gamma RT=\tfrac{\gamma-1}{2}M^2h$,
     because $h=C_pT=\gamma RT/(\gamma-1)$. · **plain** The energy balance in enthalpies and Mach numbers.
  3. **did** Mass with pressure, Mach number and temperature · **tex** $\dfrac{p_1M_1}{\sqrt{T_1}}=\dfrac{p_2M_2}{\sqrt{T_2}}$ · **why**
     $\rho u=\dfrac{p}{RT}M\sqrt{\gamma RT}=pM\sqrt{\gamma/(RT)}$ (as in D19 step 2); the constants cancel. · **plain** Mass flux
     is pressure × Mach number over √temperature — this ties the three laws together.
  4. **did** Solve mass for the Mach-number ratio · **tex** $\dfrac{M_2}{M_1}=\dfrac{p_1}{p_2}\sqrt{\dfrac{T_2}{T_1}}$ · **why**
     Rearrangement. The two ratios on the right are what steps 1 and 2 give. · **plain** The Mach number changes through
     the pressure ratio and the temperature ratio.
  5. **did** Insert the pressure ratio · **tex** $\dfrac{p_1}{p_2}=\dfrac{1+\gamma M_2^2}{1+\gamma M_1^2-f}$ · **why** Step 1. ·
     **plain** From momentum.
  6. **did** Insert the temperature ratio · **tex**
     $\dfrac{M_2}{M_1}=\dfrac{1+\gamma M_2^2}{1+\gamma M_1^2-f}\left[\dfrac{1+\frac{\gamma-1}{2}M_1^2+q}{1+\frac{\gamma-1}{2}M_2^2}\right]^{1/2}$ · **why** Step 2 with
     $h_2/h_1=T_2/T_1$. · **plain** One equation with $M_2$ on both sides.
  7. **did** Square and put each station on its own side · **tex** $A\equiv\dfrac{M_1^2\left(1+\frac12(\gamma-1)M_1^2+q\right)}{\left(1+\gamma M_1^2-f\right)^2}=\dfrac{M_2^2\left(1+\frac12(\gamma-1)M_2^2\right)}{\left(1+\gamma M_2^2\right)^2}$ (unnumbered, §15.8) · **why** Squaring removes
     the root (both sides are positive); then divide by $\left(1+\gamma M_2^2\right)^2$ and multiply by $M_1^2$ times the
     inlet bracket. Everything known is on one side: call it A. · **plain** One pure number A — fixed by the inlet, the
     friction and the heat — determines the outlet. (A is not an area.)
  8. **did** Write the outlet side as an equation for $X\equiv M_2^2$ · **tex**
     $X\left(1+\tfrac{\gamma-1}{2}X\right)=A\left(1+\gamma X\right)^2$ · **why** Multiply step 7 by the denominator. · **plain** A
     quadratic for the square of the outlet Mach number.
  9. **did** Bring it to standard form · **tex** $\left[\tfrac{\gamma-1}{2}-A\gamma^2\right]X^2+\left(1-2A\gamma\right)X-A=0$ · **why**
     Expand $(1+\gamma X)^2=1+2\gamma X+\gamma^2X^2$ and collect powers of X. · **plain** Standard quadratic form.
  10. **did** Apply the quadratic formula · **tex**
      $X=\dfrac{-(1-2A\gamma)\pm\sqrt{(1-2A\gamma)^2+4A\left[\tfrac{\gamma-1}{2}-A\gamma^2\right]}}{2\left[\tfrac{\gamma-1}{2}-A\gamma^2\right]}$ · **why** The formula
      for $aX^2+bX+c=0$ with $c=-A$. · **plain** Two candidate outlets.
  11. **did** Simplify the discriminant · **tex**
      $(1-2A\gamma)^2+2A(\gamma-1)-4A^2\gamma^2=1-2A(\gamma+1)$ · **why** $1-4A\gamma+4A^2\gamma^2+2A\gamma-2A-4A^2\gamma^2$: the $A^2$
      terms cancel. · **plain** The discriminant is linear in A.
  12. **did** Write the result · **tex** $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46) ·
      **why** Step 11 in step 10; the 2 of the denominator is multiplied in. · **plain** The outlet Mach number in closed
      form: two roots while $A<1/(2(\gamma+1))$, one when equal, none beyond.
  13. **did** Test it with no friction and no heat · **tex**
      $f=q=0:\quad M_2^2=M_1^2\quad\text{or}\quad M_2^2=\dfrac{(\gamma-1)M_1^2+2}{2\gamma M_1^2+1-\gamma}$ · **why** Then
      $1-2A(\gamma+1)=\left(\dfrac{M_1^2-1}{1+\gamma M_1^2}\right)^2$, a perfect square; its root is $\lvert M_1^2-1\rvert/(1+\gamma M_1^2)$,
      which is why the sign that gives $M_2=M_1$ is + for $M_1<1$ and − for $M_1>1$. · **plain** With nothing acting, the
      duct either does nothing or holds a normal shock — exactly $M_2^2=\dfrac{(\gamma-1)M_1^2+2}{2\gamma M_1^2+1-\gamma}$ (15.40).
- **Result:** $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46) with $A\equiv\dfrac{M_1^2\left(1+\frac12(\gamma-1)M_1^2+q\right)}{\left(1+\gamma M_1^2-f\right)^2}=\dfrac{M_2^2\left(1+\frac12(\gamma-1)M_2^2\right)}{\left(1+\gamma M_2^2\right)^2}$ (unnumbered, §15.8).
- **Check:** units — pure numbers ✓. Limits — step 13; and $A=1/(2(\gamma+1))$: $M_2=1$ (D24). Number — $M_1$ = 0.5, f = 0.1:
  A = 0.168, discriminant 0.1936, $M_2$ = 0.589 (and 1.936, excluded). Residuals of $\rho_1u_1=\rho_2u_2,\ \ p_1+\rho_1u_1^2=p_2+\rho_2u_2^2+p_1f,\ \ h_1+\tfrac12u_1^2+h_1q=h_2+\tfrac12u_2^2$ at that outlet < 10⁻¹²
  (`GD.duct_residuals`).
- **sympy check** (`check_src`; re-runs steps 8–13; executed: prints `-2*A*gamma - 2*A + 1 | (M_1 - 1)**2*(M_1 +
  1)**2/(M_1**2*gamma + 1)**2`):
  ```python
  import sympy as sp                                   # symbolic algebra
  g, A, X = sp.symbols("gamma A X", positive=True)     # gamma; the duct number A; X = M2^2
  M1 = sp.symbols("M_1", positive=True)                # inlet Mach number
  eq = X*(1 + (g-1)/2*X) - A*(1 + g*X)**2              # step 8: X(1 + (gamma-1)X/2) = A (1 + gamma X)^2
  P = sp.Poly(sp.expand(eq), X)                        # step 9: collect powers of X
  qa, qb, qc = P.all_coeffs()                          # the three coefficients of the quadratic
  assert sp.simplify(qa - ((g-1)/2 - A*g**2)) == 0 and sp.simplify(qb - (1 - 2*A*g)) == 0 and qc == -A
  disc = sp.factor(sp.simplify(qb**2 - 4*qa*qc))       # step 11: the discriminant
  assert sp.simplify(disc - (1 - 2*A*(g+1))) == 0      # ... collapses to 1 - 2A(gamma+1)
  roots = [(-(1 - 2*A*g) + s*sp.sqrt(1 - 2*A*(g+1)))/((g-1) - 2*A*g**2) for s in (1, -1)]   # Eq. (15.46)
  for X_root in roots:                                 # each root satisfies the step-8 equation
      assert sp.simplify(eq.subs(X, X_root)) == 0
  A0 = M1**2*(1 + (g-1)/2*M1**2)/(1 + g*M1**2)**2      # step 13: A with f = q = 0
  shock = ((g-1)*M1**2 + 2)/(2*g*M1**2 + 1 - g)        # Eq. (15.40): M2^2 behind a normal shock
  assert sp.simplify(eq.subs({A: A0, X: M1**2})) == 0  # "nothing happens" is a root
  assert sp.simplify(eq.subs({A: A0, X: shock})) == 0  # the normal shock is the other root
  assert sp.simplify(sp.factor((1 - 2*A0*(g+1))) - (M1**2 - 1)**2/(1 + g*M1**2)**2) == 0   # discriminant: a perfect square
  print(disc, "|", sp.factor(1 - 2*A0*(g+1)))          # the discriminant in general and for f = q = 0
  ```
- **What it means:** the duct with friction and heat is the normal shock's bigger sibling: the same three laws with two
  more terms. Friction and heat both *raise* A; when A reaches 1/(2(γ + 1)) the roots merge at Mach 1 — choking (C12).
  **Fails when:** f and q are not small lumps but distributed with changing properties (then integrate the differentials
  of C12), or the area changes.
- **Traps:** F is the *retarding* force here (trap T11); A is a number (trap T2); squaring in step 7 would admit a
  negative root — both sides are positive, so none is added; the continuous root changes sign of ± at $M_1=1$ (trap
  T12); numerator and denominator of $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ vanish together at $2A\gamma^2=\gamma-1$ (take the limit there).


### D22 · A little friction: the differentials (15.47), corrected — ★★, 9 steps, in C12 (notebook) — the book gives the result only, with a sign slip
- **Goal:** find which way speed, density and pressure change when a little wall friction acts on an adiabatic duct
  flow. · **Start:** $\rho_1u_1=\rho_2u_2,\ \ p_1+\rho_1u_1^2=p_2+\rho_2u_2^2+p_1f,\ \ h_1+\tfrac12u_1^2+h_1q=h_2+\tfrac12u_2^2$ (15.45, with $f=F/(p_1A)$, $q=Q/h_1$) with q = 0, applied between two stations a small distance apart, so that each law becomes a
  statement about differentials: $\rho u=\text{const}$; $dp+d(\rho u^2)=-p_1\,df$; $h+\tfrac12u^2=\text{const}$. · **Plan:** • take
  the differential of each law • express $d\rho$, $dT$, then $dp$ by $du$ • put them into momentum • solve. · **Tools:**
  differentials of algebraic laws as a small linear system (primer P362); logarithmic differentiation (P52 reminder);
  dynamic pressure in a gas (primer P359). · **Assumptions:** constant area; adiabatic (q = 0); perfect gas; df > 0 (friction
  retards).
- **Steps:**
  1. **did** Differentiate mass · **tex** $\dfrac{d\rho}{\rho}=-\dfrac{du}{u}$ · **why** $\rho u$ is constant in a constant-area
     duct: $d\rho/\rho+du/u=0$. · **plain** If the gas speeds up by 1 %, it thins by 1 %.
  2. **did** Differentiate momentum · **tex** $dp+\rho u\,du=-p_1\,df$ · **why** $d(\rho u^2)=d(\rho u\cdot u)=\rho u\,du$ because
     $\rho u$ is constant; the friction of the slice is $p_1\,df$. · **plain** Pressure drop plus momentum gain equals
     the friction on the slice.
  3. **did** Differentiate energy · **tex** $\dfrac{dT}{T}=-(\gamma-1)M^2\,\dfrac{du}{u}$ · **why** $C_p\,dT+u\,du=0$; divide by
     $C_pT$ and use $u^2/(C_pT)=(\gamma-1)M^2$. · **plain** Speeding up cools the gas.
  4. **did** Differentiate the gas law · **tex** $\dfrac{dp}{p}=\dfrac{d\rho}{\rho}+\dfrac{dT}{T}$ · **why** Logarithmic differential
     of $p=\rho RT$. · **plain** Pressure changes through density and temperature.
  5. **did** Insert steps 1 and 3 · **tex** $\dfrac{dp}{p}=-\left[1+(\gamma-1)M^2\right]\dfrac{du}{u}$ · **why** Both contributions are
     negative when du > 0 — they *add*. This is where the + in the corrected numerator comes from. · **plain** When the
     gas speeds up its pressure falls, by more than its density does.
  6. **did** Write momentum in relative changes · **tex** $\dfrac{dp}{p}+\gamma M^2\,\dfrac{du}{u}=-\dfrac{p_1}{p}\,df$ · **why** Divide
     step 2 by p; $\rho u^2/p=\gamma M^2$ (primer P359). · **plain** The same balance, in percentages.
  7. **did** Eliminate dp · **tex** $\dfrac{du}{u}=-\dfrac{d\rho}{\rho}=\dfrac{p_1}{p}\,\dfrac{df}{1-M^2}$ (15.47, second relation) · **why**
     Step 5 in step 6: $\left[-1-(\gamma-1)M^2+\gamma M^2\right]du/u=(M^2-1)\,du/u=-(p_1/p)\,df$; and step 1. · **plain**
     Friction speeds up a subsonic flow and slows a supersonic one.
  8. **did** Return for the pressure · **tex** $-\dfrac{dp}{p_1}=\dfrac{1+(\gamma-1)M^2}{1-M^2}\,df$ (15.47, first relation, numerator sign
     corrected) · **why** Step 7 in step 5, times $p/p_1$. (The book prints $1-(\gamma-1)M^2$ in the numerator; step 5
     shows the sign must be +.) · **plain** In subsonic flow friction lowers the pressure — steeply near Mach 1.
  9. **did** Read the signs · **tex** $df>0:\quad M<1\Rightarrow du>0,\ dp<0;\qquad M>1\Rightarrow du<0,\ dp>0$ · **why** df is
     always positive, so the sign of $1-M^2$ decides everything. · **plain** From either side, friction pushes the Mach
     number toward 1.
- **Result:** $-\dfrac{dp}{p_1}=\dfrac{1+(\gamma-1)M^2}{1-M^2}\,df,\ \ \dfrac{du}{u}=-\dfrac{d\rho}{\rho}=\dfrac{p_1}{p}\dfrac{df}{1-M^2}$ (15.47, numerator sign corrected).
- **Check:** units — numbers ✓. Limit — M → 0: $-dp/p_1=df$, i.e. $-dp\,A=dF$: the pressure drop balances the friction, as
  in a pipe of incompressible fluid. Numbers — df = 0.01 at M = 0.5: du/u = +1.33 %, −dp/p = 1.47 %. Independent check —
  finite differences of $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46) in f give 1.4667, 6.968, −0.8667 at M = 0.5, 0.9, 2: the corrected coefficient. The printed one
  gives 1.2000, 3.558, +0.2000.
- **What it means:** the chapter's second surprise is the first one in disguise: mass forces ρu to stay fixed, and the
  same $1-M^2$ decides whether a loss of pressure shows up as more speed or less. **Fails when:** heat is exchanged
  (D23), or the area changes.
- **Traps:** the book's printed numerator (slip #11); the forms hold at the station where $p=p_1$ — over a long duct
  integrate, do not multiply (trap T11); f is not a friction factor.

### D23 · A little heat: the differentials (15.48), corrected — ★★, 8 steps, in C12 (notebook · `duct_friction_heating`) — the book gives the result only, with a sign slip
- **Goal:** find which way temperature and speed change when a little heat is added to a frictionless duct flow. ·
  **Start:** $\rho_1u_1=\rho_2u_2,\ \ p_1+\rho_1u_1^2=p_2+\rho_2u_2^2+p_1f,\ \ h_1+\tfrac12u_1^2+h_1q=h_2+\tfrac12u_2^2$ (15.45, with $f=F/(p_1A)$, $q=Q/h_1$) with f = 0 between neighbouring stations: $\rho u=\text{const}$; $p+\rho u^2=\text{const}$;
  $dh+u\,du=h_1\,dq$. · **Plan:** • mass and momentum give dp and dT in terms of du • the energy law then gives du in terms
  of dq • go back for dT. · **Tools:** differentials as a linear system (primer P362); $u^2=(\gamma-1)M^2h$. ·
  **Assumptions:** constant area; no friction; perfect gas; dq > 0 means heating.
- **Steps:**
  1. **did** Differentiate mass · **tex** $\dfrac{d\rho}{\rho}=-\dfrac{du}{u}$ · **why** As in D22 step 1. · **plain** Faster means
     thinner.
  2. **did** Differentiate momentum · **tex** $\dfrac{dp}{p}=-\gamma M^2\,\dfrac{du}{u}$ · **why** With no friction $dp+\rho u\,du=0$;
     divide by p and use $\rho u^2/p=\gamma M^2$. · **plain** Faster means lower pressure.
  3. **did** Combine with the gas law · **tex** $\dfrac{dT}{T}=\dfrac{dp}{p}-\dfrac{d\rho}{\rho}=\left(1-\gamma M^2\right)\dfrac{du}{u}$ · **why**
     $p=\rho RT$ in logarithmic form; steps 1 and 2. · **plain** Whether speeding up warms or cools the gas depends on
     whether $\gamma M^2$ is below or above 1.
  4. **did** Differentiate energy with heat · **tex** $C_p\,dT+u\,du=h_1\,dq$ · **why** The energy law of $\rho_1u_1=\rho_2u_2,\ \ p_1+\rho_1u_1^2=p_2+\rho_2u_2^2+p_1f,\ \ h_1+\tfrac12u_1^2+h_1q=h_2+\tfrac12u_2^2$ between
     neighbouring stations; $h=C_pT$. · **plain** Heat added goes into enthalpy and kinetic energy.
  5. **did** Divide by the local enthalpy · **tex** $\dfrac{dT}{T}+(\gamma-1)M^2\,\dfrac{du}{u}=\dfrac{h_1}{h}\,dq$ · **why**
     $u^2/h=(\gamma-1)M^2$ (D21 step 2). · **plain** The same in relative changes.
  6. **did** Eliminate dT · **tex** $\dfrac{du}{u}=\dfrac{h_1}{h}\,\dfrac{dq}{1-M^2}$ · **why** Step 3 in step 5:
     $\left[1-\gamma M^2+(\gamma-1)M^2\right]du/u=(1-M^2)\,du/u$. · **plain** Heating speeds up a subsonic flow and slows a
     supersonic one — exactly like friction.
  7. **did** Return for the temperature · **tex** $\dfrac{dT}{T_1}=\left(\dfrac{1-\gamma M^2}{1-M^2}\right)dq$ (15.48, first relation) · **why**
     Step 6 in step 3, at the station where $h=h_1$ (so $T=T_1$). · **plain** Between $M^2=1/\gamma$ and 1 the bracket is
     negative: heat in, temperature down.
  8. **did** Write the speed change as a kinetic-energy change · **tex** $\dfrac{u\,du}{h_1}=\dfrac{(\gamma-1)M^2}{1-M^2}\,dq$ (15.48, second
     relation's sign corrected) · **why** Multiply step 6 by $u^2/h_1=(\gamma-1)M^2$ at $h=h_1$. (The book prints
     $(1-\gamma)M^2$; that would make heating slow a subsonic flow, against its own text.) · **plain** The fraction of the
     heat that becomes kinetic energy; above $M^2=1/\gamma$ it exceeds 1.
- **Result:** $\dfrac{dT}{T_1}=\left(\dfrac{1-\gamma M^2}{1-M^2}\right)dq,\ \ \dfrac{u\,du}{h_1}=\dfrac{(\gamma-1)M^2}{1-M^2}\,dq$ (15.48, second relation's sign corrected).
- **Check:** units — numbers ✓. Sum rule — the two relations add to $dT/T_1+u\,du/h_1=dq$ ✓ (the energy law):
  $(1-\gamma M^2+(\gamma-1)M^2)/(1-M^2)=1$. With the printed sign they would not. Limit — M → 0: $dT/T_1=dq$ (all heat
  warms the gas). Numbers — dq = 0.01 at M = 0.9: dT/T = −0.71 %; finite differences of $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46) give 0.1333, 1.7053, −0.5333
  for $u\,du/(h_1dq)$ at M = 0.5, 0.9, 2: the corrected coefficient.
- **What it means:** the third surprise. In the band $1/\sqrt\gamma<M<1$ (0.845 to 1 for air) more than all of the added
  heat goes into speed, and the gas cools. **Fails when:** friction acts too, or the heat release changes the gas
  (combustion changes γ and R).
- **Traps:** the book's printed sign (slip #12); dT changes sign at $M^2=1/\gamma$, not at M = 1 (trap T18); entropy
  still rises with dq > 0; q is the dimensionless $Q/h_1$.

### D24 · The choking limits of friction and of heat — ★★, 7 steps, in C12 (notebook) — ours
- **Goal:** find the largest friction, and the largest heat, a duct can take for a given inlet Mach number. · **Start:**
  $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46) and $A\equiv\dfrac{M_1^2\left(1+\frac12(\gamma-1)M_1^2+q\right)}{\left(1+\gamma M_1^2-f\right)^2}=\dfrac{M_2^2\left(1+\frac12(\gamma-1)M_2^2\right)}{\left(1+\gamma M_2^2\right)^2}$ (unnumbered, §15.8). · **Plan:** • the two roots merge when the discriminant vanishes • show the merged root is Mach 1 •
  solve "A = that value" for f, then for q. · **Tools:** a quadratic in M squared and its discriminant (primer P366). ·
  **Assumptions:** as D21.
- **Steps:**
  1. **did** Set the discriminant to zero · **tex** $1-2A(\gamma+1)=0\ \Rightarrow\ A=\dfrac{1}{2(\gamma+1)}$ · **why** A quadratic has a
     double root exactly when its discriminant vanishes; beyond that A there is no real root. · **plain** The largest A a
     duct flow can have: 0.208 for air.
  2. **did** Find the merged root · **tex** $M_2^2=\dfrac{-(1-2A\gamma)}{(\gamma-1)-2A\gamma^2}=\dfrac{-1/(\gamma+1)}{-1/(\gamma+1)}=1$ · **why** With
     that A: $1-2A\gamma=1/(\gamma+1)$ and $(\gamma-1)-\gamma^2/(\gamma+1)=-1/(\gamma+1)$. · **plain** At the limit the outlet is
     exactly sonic: the duct is choked.
  3. **did** Impose the limit with friction only · **tex**
     $\dfrac{M_1^2\,\Theta_1}{\left(1+\gamma M_1^2-f\right)^2}=\dfrac{1}{2(\gamma+1)},\qquad\Theta_1\equiv1+\tfrac{\gamma-1}{2}M_1^2$ · **why** $A\equiv\dfrac{M_1^2\left(1+\frac12(\gamma-1)M_1^2+q\right)}{\left(1+\gamma M_1^2-f\right)^2}=\dfrac{M_2^2\left(1+\frac12(\gamma-1)M_2^2\right)}{\left(1+\gamma M_2^2\right)^2}$ (unnumbered, §15.8) with
     q = 0 set equal to step 1. · **plain** An equation for the choking friction.
  4. **did** Take the square root · **tex** $1+\gamma M_1^2-f=\pm M_1\sqrt{2(\gamma+1)\,\Theta_1}$ · **why** Cross-multiply and
     take roots; both signs are possible algebraically. · **plain** Two candidates.
  5. **did** Choose the sign · **tex** $f_{\max}=1+\gamma M_1^2-M_1\sqrt{2(\gamma+1)\left(1+\tfrac{\gamma-1}{2}M_1^2\right)}$ (ours) · **why** By
     D21 step 1, $1+\gamma M_1^2-f=(p_2/p_1)(1+\gamma M_2^2)>0$: the + sign. · **plain** The friction that brings the outlet
     to Mach 1.
  6. **did** Impose the limit with heat only · **tex**
     $\dfrac{M_1^2\left(\Theta_1+q\right)}{\left(1+\gamma M_1^2\right)^2}=\dfrac{1}{2(\gamma+1)}$ · **why** $A\equiv\dfrac{M_1^2\left(1+\frac12(\gamma-1)M_1^2+q\right)}{\left(1+\gamma M_1^2-f\right)^2}=\dfrac{M_2^2\left(1+\frac12(\gamma-1)M_2^2\right)}{\left(1+\gamma M_2^2\right)^2}$ (unnumbered, §15.8) with f = 0. · **plain** An
     equation for the choking heat — linear in q.
  7. **did** Solve for q · **tex** $q_{\max}=\dfrac{(1+\gamma M_1^2)^2}{2(\gamma+1)M_1^2}-1-\dfrac{\gamma-1}{2}M_1^2$ (ours) · **why** Multiply
     up, divide by $M_1^2$, subtract $\Theta_1$. · **plain** The heat that brings the outlet to Mach 1.
- **Result:** $f_{\max}=1+\gamma M_1^2-M_1\sqrt{2(\gamma+1)\left(1+\tfrac{\gamma-1}{2}M_1^2\right)}$ (ours) and $q_{\max}=\dfrac{(1+\gamma M_1^2)^2}{2(\gamma+1)M_1^2}-1-\dfrac{\gamma-1}{2}M_1^2$ (ours).
- **Check:** units — numbers ✓. Limit — $M_1=1$: $f_{\max}=1+\gamma-\sqrt{(\gamma+1)^2}=0$ and $q_{\max}=\tfrac{\gamma+1}{2}-
  \tfrac{\gamma+1}{2}=0$: a sonic inlet can take nothing. Number — $M_1$ = 0.5: $f_{\max}$ = 1.35 − 0.5 × √(4.8 × 1.05) = 0.2275;
  $q_{\max}$ = 1.8225/1.2 − 1.05 = 0.4688. Cross-check (ours) — $1+q_{\max}/\Theta_1$ = 1.4464 = the Rayleigh $T_0^*/T_0$ at
  M = 0.5.
- **What it means:** friction and heat choke a duct as a throat does; ask for more and the *inlet* must change (the mass
  flow drops, or a shock moves upstream). Cooling (q < 0) lowers A and never chokes. **Fails when:** as D21.
- **Traps:** beyond the limit there is no steady solution with this inlet — not "a solution with $M_2>1$"; choked still
  does not mean the flow cannot be raised by changing the supply (trap T10).

### D25 · The Mach angle (15.49) — ★, 4 steps, in C13 (notebook · `mach_cone`)
- **Goal:** find the half-angle of the cone behind a source that moves faster than sound. · **Start:** a source moving
  at speed U through still gas of sound speed c, emitting all the time; a front emitted at one instant is a sphere
  growing at c from where the source was. · **Plan:** • follow one front for a time Δt • the cone touches it • a right
  triangle. · **Tools:** the envelope of a family of circles (primer P368). · **Assumptions:** uniform gas; constant U > c;
  weak disturbances (they travel at c).
- **Steps:**
  1. **did** Follow one front for a time Δt · **tex** $\text{radius}=c\,\Delta t,\qquad\text{distance to the source}=U\,\Delta t$ ·
     **why** The front grows at the sound speed from the emission point; meanwhile the source moves on at U. · **plain**
     After Δt the source is UΔt ahead of the centre of a circle of radius cΔt.
  2. **did** Draw the cone as the common tangent · **tex** $\text{tangent}\perp\text{radius}$ · **why** Every front's radius is the same fraction c/U of its distance behind the source, so the envelope is the pair of tangents from the source; a tangent is perpendicular to the radius (primer P368). · **plain** Source, centre and
     touching point form a right triangle.
  3. **did** Read the triangle · **tex** $\sin\mu=\dfrac{c\,\Delta t}{U\,\Delta t}$ · **why** μ is the angle at the source between the
     flight path and the tangent: opposite side cΔt, hypotenuse UΔt. · **plain** The sine of the half-angle is sound's
     distance over the source's distance.
  4. **did** Cancel Δt · **tex** $\sin\mu=\dfrac1M$ (15.49) · **why** $M=U/c$; the result is the same for every front, which is
     why the envelope is a straight cone. Needs M ≥ 1. · **plain** The faster the source, the narrower the cone.
- **Result:** $\sin\mu=1/M$ (15.49).
- **Check:** units — a number ✓. Limits — M = 1: μ = 90° (a plane front through the source); M → ∞: μ → 0. Number — M = 2:
  30°. Code — the envelope of 200 computed circles has the slope tan 30°.
- **What it means:** the two-dimensional form of "no news upstream": outside the cone nothing has been heard. Mach lines
  are the paths along which every weak disturbance of a steady supersonic flow travels. **Fails when:** the source is
  strong (the front is then a shock, steeper than μ — C14) or the gas is not uniform (the cone bends).
- **Traps:** μ is measured from the direction of motion; there is no real angle for M < 1; μ is the Mach angle here, not
  a viscosity (trap T5).

### D26 · An oblique shock is a normal shock on the normal component: (15.50), (15.51), (15.52) — ★★, 8 steps, in C14 (notebook · `oblique_shock_wedge`)
- **Goal:** show that every jump relation of the normal shock holds for an oblique shock if the Mach numbers are
  replaced by their components normal to the shock. · **Start:** a plane shock at angle σ to the upstream velocity
  $V_1$, turning the flow through δ; the control-volume laws of C03 for a thin box enclosing a piece of the shock. ·
  **Plan:** • split both velocities along and across the shock • mass and tangential momentum show the along-component
  is unchanged • it then drops out of the energy law • what is left is the normal-shock problem. · **Tools:** resolving
  a velocity along and across a tilted line (primer P370); the Galilean shift (P96 reminder); D13–D15. · **Assumptions:**
  steady, planar, adiabatic, no forces along the shock other than pressure (which is the same on both ends of the box).
- **Steps:**
  1. **did** Resolve the velocities · **tex**
     $u_1=V_1\sin\sigma,\ v_1=V_1\cos\sigma;\qquad u_2=V_2\sin(\sigma-\delta),\ v_2=V_2\cos(\sigma-\delta)$ · **why** u is across the
     shock, v along it; behind the shock the flow makes the angle σ − δ with the shock. · **plain** Each velocity is split
     into a part that crosses the shock and a part that slides along it.
  2. **did** Mass through the box · **tex** $\rho_1u_1=\rho_2u_2$ · **why** Only the normal component carries mass through
     the shock. · **plain** The same as for a normal shock, with u the normal component.
  3. **did** Momentum along the shock · **tex** $\rho_1u_1v_1=\rho_2u_2v_2\ \Rightarrow\ v_1=v_2\equiv v$ · **why** The flux of
     tangential momentum is (mass flux) × v; no net force acts along the shock. Divide by the equal mass fluxes of step 2.
     · **plain** The shock does not change the speed along itself.
  4. **did** Momentum across the shock, and the normal Mach numbers · **tex**
     $p_1+\rho_1u_1^2=p_2+\rho_2u_2^2;\qquad M_{n1}=\dfrac{u_1}{c_1}=M_1\sin\sigma,\quad M_{n2}=\dfrac{u_2}{c_2}=M_2\sin(\sigma-\delta)$ · **why**
     The normal momentum balance involves u only; dividing step 1 by each side's sound speed gives the normal Mach
     numbers. · **plain** Across the shock, pressure plus normal momentum flux is conserved.
  5. **did** Energy · **tex** $h_1+\tfrac12\left(u_1^2+v^2\right)=h_2+\tfrac12\left(u_2^2+v^2\right)\ \Rightarrow\ h_1+\tfrac12u_1^2=h_2+\tfrac12u_2^2$ · **why**
     $h_2+\tfrac12u_2^2-h_1-\tfrac12u_1^2=Q$ (15.22) with Q = 0 and the full speed; the same $\tfrac12v^2$ stands on both sides and cancels. · **plain** The energy
     law also involves only the normal component.
  6. **did** Recognise the normal-shock problem; pressure · **tex**
     $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2\sin^2\sigma-1\right)$ (15.50) · **why** Steps 2, 4, 5 are the three jump laws $\rho_1u_1=\rho_2u_2$ (15.35), $p_1-p_2=-\rho_1u_1^2+\rho_2u_2^2$ (15.36), $h_1+\tfrac12u_1^2=h_2+\tfrac12u_2^2$ (15.37) for $(u_1,u_2)$, so $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39) holds with $M_1\to M_{n1}$: an observer sliding along the shock at v sees a normal shock. · **plain** The pressure jump of a normal shock at the normal Mach number.
  7. **did** Density, and the geometry · **tex**
     $\dfrac{\rho_2}{\rho_1}=\dfrac{(\gamma+1)M_1^2\sin^2\sigma}{(\gamma-1)M_1^2\sin^2\sigma+2}=\dfrac{u_1}{u_2}=\dfrac{\tan\sigma}{\tan(\sigma-\delta)}$ (15.51) · **why**
     $\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}$ (15.41) with $M_{n1}$; and from step 1 with equal v: $u_1=v\tan\sigma$, $u_2=v\tan(\sigma-\delta)$. · **plain** The density
     ratio equals a ratio of two tangents — this ties the shock angle to the turn.
  8. **did** Downstream Mach number · **tex**
     $M_2^2\sin^2(\sigma-\delta)=\dfrac{(\gamma-1)M_1^2\sin^2\sigma+2}{2\gamma M_1^2\sin^2\sigma+1-\gamma}$ (15.52) · **why** $M_2^2=\dfrac{(\gamma-1)M_1^2+2}{2\gamma M_1^2+1-\gamma}$ (15.40) with both normal Mach
     numbers of step 4. · **plain** The *normal* Mach number behind is subsonic; $M_2$ itself is that divided by
     sin(σ − δ), and may exceed 1.
- **Result:** $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2\sin^2\sigma-1\right)$ (15.50), $\dfrac{\rho_2}{\rho_1}=\dfrac{(\gamma+1)M_1^2\sin^2\sigma}{(\gamma-1)M_1^2\sin^2\sigma+2}=\dfrac{\tan\sigma}{\tan(\sigma-\delta)}$ (15.51), $M_2^2\sin^2(\sigma-\delta)=\dfrac{(\gamma-1)M_1^2\sin^2\sigma+2}{2\gamma M_1^2\sin^2\sigma+1-\gamma}$ (15.52).
- **Check:** units — numbers ✓. Limits — σ = 90°: the normal shock; σ = μ ($M_1\sin\sigma=1$): no jump — a Mach wave.
  Number — $M_1$ = 2, σ = 39.31°: $M_{n1}$ = 1.267, $p_2/p_1$ = 1.707, density ratio 1.458 = tan 39.31°/tan 29.31° ✓, $M_2$ =
  0.803/0.4896 = 1.641. Code — `ch15.tangential_momentum_sympy()`.
- **What it means:** all the tables and formulas of the normal shock serve for oblique shocks with $M_1\sin\sigma$. The
  flow turns *toward* the shock because only its normal component is slowed. **Fails when:** the shock is curved (then v
  varies along it and the flow behind is rotational).
- **Traps:** the mass flux uses the normal component only; $M_2$ needs the division by sin(σ − δ) (trap T13); σ and δ are
  β and θ in most other books.

### D27 · Deflection against shock angle (15.53) — ★★, 8 steps, in C14 (notebook · `oblique_shock_wedge`)
- **Goal:** find the turn δ produced by a shock at angle σ. · **Start:** the last equality of $\dfrac{\rho_2}{\rho_1}=\dfrac{(\gamma+1)M_1^2\sin^2\sigma}{(\gamma-1)M_1^2\sin^2\sigma+2}=\dfrac{\tan\sigma}{\tan(\sigma-\delta)}$ (15.51),
  $\dfrac{\tan\sigma}{\tan(\sigma-\delta)}=\dfrac{\rho_2}{\rho_1}\equiv\hat\rho$. · **Plan:** • expand tan(σ − δ) • collect tan δ • insert
  the density ratio • tidy with a double angle. · **Tools:** the tangent subtraction formula (primer P371). ·
  **Assumptions:** as D26.
- **Steps:**
  1. **did** Solve for the downstream tangent · **tex** $\tan(\sigma-\delta)=\dfrac{\tan\sigma}{\hat\rho}$ · **why** Rearrangement. ·
     **plain** The flow leaves at a flatter angle to the shock, flatter by the density ratio.
  2. **did** Expand the left side · **tex** $\dfrac{\tan\sigma-\tan\delta}{1+\tan\sigma\tan\delta}=\dfrac{\tan\sigma}{\hat\rho}$ · **why** The
     subtraction formula (primer P371). Now tan δ is explicit. · **plain** An equation with tan δ in two places.
  3. **did** Cross-multiply · **tex** $\hat\rho\left(\tan\sigma-\tan\delta\right)=\tan\sigma\left(1+\tan\sigma\tan\delta\right)$ · **why** Both
     denominators are positive for 0 < δ < σ < 90°. · **plain** No fractions.
  4. **did** Collect tan δ · **tex** $\tan\delta\left(\hat\rho+\tan^2\sigma\right)=\left(\hat\rho-1\right)\tan\sigma$ · **why** Move the two terms
     with tan δ to one side. · **plain** tan δ appears once.
  5. **did** Solve · **tex** $\tan\delta=\dfrac{\left(\hat\rho-1\right)\tan\sigma}{\hat\rho+\tan^2\sigma}$ · **why** Divide. Already a usable
     result: no compression ($\hat\rho=1$), no turn. · **plain** The turn in terms of the density ratio and the shock
     angle.
  6. **did** Insert the density ratio · **tex**
     $\tan\delta=\dfrac{2\,(m-1)\tan\sigma}{(\gamma+1)\,m+\left[(\gamma-1)\,m+2\right]\tan^2\sigma},\qquad m\equiv M_1^2\sin^2\sigma$ · **why** From $\dfrac{\rho_2}{\rho_1}=\dfrac{(\gamma+1)M_1^2\sin^2\sigma}{(\gamma-1)M_1^2\sin^2\sigma+2}=\dfrac{\tan\sigma}{\tan(\sigma-\delta)}$ (15.51),
     $\hat\rho=\dfrac{(\gamma+1)m}{(\gamma-1)m+2}$, so $\hat\rho-1=\dfrac{2(m-1)}{(\gamma-1)m+2}$; multiply top and bottom by $(\gamma-1)m+2$. ·
     **plain** Everything in $M_1$, σ and γ.
  7. **did** Multiply top and bottom by cos²σ/sin²σ · **tex**
     $\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{(\gamma+1)M_1^2\cos^2\sigma+(\gamma-1)M_1^2\sin^2\sigma+2}$ · **why** $\tan\sigma\cos^2\sigma/\sin^2\sigma=\cot\sigma$;
     $m/\sin^2\sigma=M_1^2$; $\tan^2\sigma\cos^2\sigma/\sin^2\sigma=1$. · **plain** The book's shape, except for the
     denominator.
  8. **did** Gather the denominator with a double angle · **tex** $\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}$ (15.53) ·
     **why** $(\gamma+1)\cos^2\sigma+(\gamma-1)\sin^2\sigma=\gamma+(\cos^2\sigma-\sin^2\sigma)=\gamma+\cos2\sigma$. · **plain** The turn
     as a function of the shock angle, for a given upstream Mach number.
- **Result:** $\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}$ (15.53).
- **Check:** units — numbers ✓. Zeros — δ = 0 where $M_1\sin\sigma=1$ (σ = μ) and where cot σ = 0 (σ = 90°). Number — $M_1$ =
  2, σ = 39.31°: 2 × 1.2213 × 0.6056/8.389 = 0.1763 → 10.0°; σ = 45°: 2 × 1 × 1/7.6 → 14.74°. Code — `GD.deflection_angle`.
- **What it means:** between its two zeros δ rises to a maximum and falls again: for each turn below the maximum there
  are two shock angles (weak and strong), and above it none — the shock detaches. **Fails when:** δ exceeds the maximum
  (no straight attached shock exists).
- **Traps:** collect tan δ before solving; the formula gives δ from σ — the inverse needs a root-finder on each side of
  the maximum (primer P360); "weak" and "strong" are names of roots (trap T13).

### D28 · The weak-wave rule: (15.54) → (15.55) — ★★, 8 steps, in C14 (notebook · `supersonic_flat_plate`)
- **Goal:** find the pressure rise produced by a very small turn of a supersonic stream. · **Start:** $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2\sin^2\sigma-1\right)$ (15.50) and $\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}$ (15.53). ·
  **Plan:** • express the numerator of the deflection formula by the pressure jump • let the turn go to zero: the shock
  becomes a Mach wave • evaluate each factor there. · **Tools:** the Mach triangle (primer P372); small-angle
  approximation (P100 reminder). · **Assumptions:** δ ≪ 1 (in radians); the weak root; $M_1$ not close to 1.
- **Steps:**
  1. **did** Solve the pressure relation for the numerator · **tex**
     $M_1^2\sin^2\sigma-1=\dfrac{\gamma+1}{2\gamma}\,\dfrac{p_2-p_1}{p_1}$ · **why** $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2\sin^2\sigma-1\right)$ (15.50) rearranged. This combination is the numerator of
     $\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}$. · **plain** How far the normal Mach number is above 1 is the pressure jump, up to a factor.
  2. **did** Insert it · **tex**
     $\tan\delta=2\cot\sigma\,\dfrac{\gamma+1}{2\gamma}\left(\dfrac{p_2-p_1}{p_1}\right)\dfrac{1}{M_1^2(\gamma+\cos2\sigma)+2}$ · **why** Substitution in $\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}$
     (15.53); exact. · **plain** The turn is proportional to the pressure jump.
  3. **did** Rewrite the double angle · **tex** $\tan\delta=2\cot\sigma\,\dfrac{\gamma+1}{2\gamma}\left(\dfrac{p_2-p_1}{p_1}\right)\dfrac{1}{M_1^2\left(\gamma+1-2\sin^2\sigma\right)+2}$ (15.54) · **why** $\cos2\sigma=1-2\sin^2\sigma$. Still exact. · **plain**
     The same, ready for the limit.
  4. **did** Small turn · **tex** $\tan\delta\approx\delta$ · **why** δ ≪ 1 rad (≈): the error is of order δ³. · **plain** A
     small angle equals its tangent.
  5. **did** The shock becomes a Mach wave · **tex** $\sin\sigma\to\sin\mu_1=\dfrac{1}{M_1}$ · **why** As δ → 0 the pressure jump →
     0, so by step 1 $M_1\sin\sigma\to1$ (the weak root). · **plain** An infinitely weak shock lies along the Mach line.
  6. **did** Its cotangent · **tex** $\cot\sigma\to\cot\mu_1=\sqrt{M_1^2-1}$ · **why** The Mach triangle (primer P372). ·
     **plain** This is where the square root comes from.
  7. **did** The denominator · **tex** $M_1^2\left(\gamma+1-2\sin^2\sigma\right)+2\ \to\ (\gamma+1)M_1^2$ · **why**
     $M_1^2\cdot2\sin^2\sigma\to2$, which cancels the +2. · **plain** A simple limit.
  8. **did** Assemble and solve for the pressure jump · **tex** $\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta$ (15.55) ·
     **why** $\delta=2\sqrt{M_1^2-1}\,\dfrac{\gamma+1}{2\gamma}\,\dfrac{p_2-p_1}{p_1}\,\dfrac{1}{(\gamma+1)M_1^2}$; the factors γ + 1 and 2 cancel. ·
     **plain** A small turn into the flow raises the pressure in proportion to the turn.
- **Result:** $\tan\delta=2\cot\sigma\,\dfrac{\gamma+1}{2\gamma}\left(\dfrac{p_2-p_1}{p_1}\right)\dfrac{1}{M_1^2\left(\gamma+1-2\sin^2\sigma\right)+2}$ (15.54) and $\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta$ (15.55).
- **Check:** units — δ in radians; both sides numbers ✓. Number — $M_1$ = 2, δ = 2° = 0.0349: 1.4 × 4/1.732 × 0.0349 = 0.1129
  (exact shock 0.1180). How small is small (trap T16): the rule is 1.1 % low at 0.5°, 10.5 % low at 5°, 20 % low at 10°,
  39 % low at 20°. Code — `GD.weak_wave_error`, `ch15.weak_wave_sympy()`.
- **What it means:** to first order pressure responds linearly to the local turn — and, because the entropy change is
  third order (D17), the same line holds for small expansions with δ < 0 (N103). This one line is the whole of §15.10.
  **Fails when:** the turn is not small, or $M_1\to1$ ($\sqrt{M_1^2-1}\to0$).
- **Traps:** three replacements, one per step; δ in radians, positive when the wall turns *into* the flow (trap T14);
  the rule underestimates compressions and overestimates expansions.

### D29 · The Prandtl–Meyer differential relation — ★★, 9 steps, in C15 (notebook · `prandtl_meyer_fan`)
- **Goal:** find how much the Mach number changes when a supersonic stream turns through one very weak wave. · **Start:**
  across a wave at angle σ the speed along the wave is unchanged (D26 step 3): $V_1\cos\sigma=V_2\cos(\sigma-\delta)$. ·
  **Plan:** • make the turn infinitesimal • get the speed change per unit turn • convert the speed change into a
  Mach-number change using constant stagnation temperature. · **Tools:** first-order expansion (P26 reminder); the Mach
  triangle (primer P372); logarithmic differentiation (P52 reminder); D08. · **Assumptions:** an infinitely weak wave,
  hence isentropic and at the Mach angle (step 4); adiabatic, perfect gas (step 7).
- **Steps:**
  1. **did** Start from the unchanged tangential speed · **tex** $V_1\cos\sigma=V_2\cos(\sigma-\delta)$ · **why** D26 step 3 with
     $v=V\cos(\text{angle to the wave})$. · **plain** The stream keeps its speed along the wave.
  2. **did** Make the turn small · **tex** $\cos(\sigma-d\delta)\approx\cos\sigma+\sin\sigma\,d\delta$ · **why** First-order Taylor
     expansion in the small angle dδ (≈); the derivative of cos is −sin. · **plain** A small turn changes the cosine in
     proportion.
  3. **did** Keep first order · **tex** $0=dV\cos\sigma+V\sin\sigma\,d\delta\ \Rightarrow\ d\delta=-\dfrac{dV}{V}\cot\sigma$ · **why** Put
     $V_1=V$, $V_2=V+dV$; the product $dV\,d\delta$ is second order and dropped. · **plain** Turning into the wave (dδ > 0)
     slows the stream; turning away speeds it up.
  4. **did** Put the wave at the Mach angle · **tex** $d\delta=-\dfrac{dV}{V}\sqrt{M^2-1}$ · **why** An infinitely weak wave is a
     Mach wave: σ = μ, and cot μ = $\sqrt{M^2-1}$ (primer P372). · **plain** The speed change per unit turn is fixed by
     the Mach number.
  5. **did** Split the speed change · **tex** $\dfrac{dV}{V}=\dfrac{dM}{M}+\dfrac{dc}{c}$ · **why** Logarithmic differential of
     $V=Mc$. · **plain** Speed changes through the Mach number and through the sound speed.
  6. **did** Sound speed from temperature · **tex** $\dfrac{dc}{c}=\dfrac12\dfrac{dT}{T}$ · **why** $c=\sqrt{\gamma RT}$: the logarithm
     of a square root. · **plain** Sound speed changes half as fast as temperature.
  7. **did** Temperature from constant $T_0$ · **tex** $\dfrac{dT}{T}=-\dfrac{(\gamma-1)M\,dM}{1+\frac{\gamma-1}{2}M^2}$ · **why** $T=T_0/\Theta$
     with $T_0$ fixed (adiabatic): $dT/T=-d\Theta/\Theta$ and $d\Theta=(\gamma-1)M\,dM$. · **plain** As the Mach number rises
     the gas cools.
  8. **did** Combine steps 5–7 · **tex** $\dfrac{dV}{V}=\dfrac{dM}{M}\cdot\dfrac{1}{1+\frac{\gamma-1}{2}M^2}$ · **why**
     $\dfrac{dM}{M}\left[1-\dfrac{\frac{\gamma-1}{2}M^2}{\Theta}\right]=\dfrac{dM}{M}\,\dfrac1\Theta$. · **plain** The speed rises more slowly than
     the Mach number, because the sound speed falls.
  9. **did** Insert into step 4 · **tex** $d\delta=-\dfrac{\sqrt{M^2-1}}{M}\,\dfrac{dM}{1+\frac12(\gamma-1)M^2}$ (unnumbered, §15.9) · **why**
     Substitution. · **plain** A small turn away from the flow (dδ < 0) raises the Mach number by a definite amount.
- **Result:** $d\delta=-\dfrac{\sqrt{M^2-1}}{M}\dfrac{dM}{1+\frac12(\gamma-1)M^2}$ (unnumbered, §15.9) — so that $d\delta+d\nu=0$ with $d\nu$ the integrand of $\nu(M)=\displaystyle\int_1^M\dfrac{\sqrt{M^2-1}}{1+\frac12(\gamma-1)M^2}\dfrac{dM}{M}$ (15.56).
- **Check:** units — radians on both sides ✓. Sign — dδ < 0 (expansion) ⇒ dM > 0. Number — at M = 2: dM = −dδ × 2 × 1.8/1.732 =
  2.078 per radian = 0.0363 per degree; ten steps of 1° from Mach 2 give ≈ 2.38 (exact 2.385). Code — `ch15.fan_staircase`.
- **What it means:** each Mach wave of a fan does this small job; adding them up is the integral of D30. **Fails when:**
  the turn is a finite compression — the waves cross and a shock (not isentropic) forms.
- **Traps:** the sign: δ > 0 is a turn *into* the flow (trap T14); $dc/c=\tfrac12dT/T$; dT comes from constant $T_0$, not
  constant T; exact for an expansion of any size, for a compression only until Mach lines cross.

### D30 · The Prandtl–Meyer function in closed form (15.56) — ★★★, 10 steps, in C15 (notebook · `prandtl_meyer_fan`) — the book gives the result only
- **Goal:** add up infinitely many small turns: integrate the differential relation from Mach 1 to M. · **Start:**
  $\nu(M)=\displaystyle\int_1^M\dfrac{\sqrt{M^2-1}}{1+\frac12(\gamma-1)M^2}\dfrac{dM}{M}$ (15.56, integral form) — "the total turning needed to bring a sonic stream to Mach number M"; by D29, $\delta+\nu(M)=\text{const}$ along the
  wall. · **Plan:** • substitute $t=\sqrt{M^2-1}$ • the integrand becomes a ratio of polynomials in t • split it into two
  simple fractions • each integrates to an arctangent • put M back. · **Tools:** substitution in an integral (P106
  reminder); partial fractions (P219 reminder); arctangent integrals (primer P373). · **Assumptions:** constant γ > 1; M ≥ 1.
- **Steps:**
  1. **did** Substitute · **tex** $t\equiv\sqrt{M^2-1}\quad(M^2=1+t^2)$ · **why** The square root is what makes the integral
     awkward; naming it removes it. t runs from 0 (M = 1) upward. · **plain** t is the cotangent of the Mach angle.
  2. **did** Convert the measure · **tex** $\dfrac{dM}{M}=\dfrac{t\,dt}{1+t^2}$ · **why** Differentiate $M^2=1+t^2$:
     $2M\,dM=2t\,dt$; divide by $M^2=1+t^2$. · **plain** A small relative change of M, in terms of t.
  3. **did** Convert the bracket · **tex**
     $1+\tfrac{\gamma-1}{2}M^2=\tfrac{\gamma+1}{2}\left(1+\lambda^2t^2\right),\qquad\lambda^2\equiv\dfrac{\gamma-1}{\gamma+1}$ · **why**
     $1+\tfrac{\gamma-1}{2}(1+t^2)=\tfrac{\gamma+1}{2}+\tfrac{\gamma-1}{2}t^2$; factor out $\tfrac{\gamma+1}{2}$. · **plain** The
     temperature bracket, in t (λ² = 1/6 for air).
  4. **did** Assemble the integral in t · **tex**
     $\nu=\displaystyle\int_0^{t}\dfrac{2}{\gamma+1}\,\dfrac{t^2}{(1+t^2)(1+\lambda^2t^2)}\,dt$ · **why** Numerator $\sqrt{M^2-1}=t$ times step 2;
     denominator from step 3; limits M = 1 → t = 0. · **plain** A ratio of polynomials — no roots left.
  5. **did** Split into partial fractions · **tex**
     $\dfrac{t^2}{(1+t^2)(1+\lambda^2t^2)}=\dfrac{1}{1-\lambda^2}\left[\dfrac{1}{1+\lambda^2t^2}-\dfrac{1}{1+t^2}\right]$ · **why** Put the bracket on a
     common denominator: its numerator is $(1+t^2)-(1+\lambda^2t^2)=(1-\lambda^2)t^2$. · **plain** One hard fraction is the
     difference of two easy ones.
  6. **did** Simplify the prefactor · **tex** $\dfrac{2}{\gamma+1}\cdot\dfrac{1}{1-\lambda^2}=1\ \Rightarrow\ \nu=\displaystyle\int_0^t\left[\dfrac{1}{1+\lambda^2t^2}-\dfrac{1}{1+t^2}\right]dt$
     · **why** $1-\lambda^2=\dfrac{(\gamma+1)-(\gamma-1)}{\gamma+1}=\dfrac{2}{\gamma+1}$. · **plain** The constants cancel exactly.
  7. **did** Integrate the first fraction · **tex** $\displaystyle\int_0^t\dfrac{dt}{1+\lambda^2t^2}=\dfrac1\lambda\tan^{-1}(\lambda t)$ · **why** The
     arctangent integral with scaling λ (primer P373); differentiate back to check. · **plain** An arctangent, stretched
     by 1/λ.
  8. **did** Integrate the second · **tex** $\displaystyle\int_0^t\dfrac{dt}{1+t^2}=\tan^{-1}t$ · **why** The same with scaling 1. ·
     **plain** A plain arctangent.
  9. **did** Apply the lower limit · **tex** $\nu=\dfrac1\lambda\tan^{-1}(\lambda t)-\tan^{-1}t\quad(\text{both vanish at }t=0)$ · **why**
     $\tan^{-1}0=0$: the lower limit contributes nothing, which is what makes ν(1) = 0. · **plain** ν in terms of t.
  10. **did** Put M back · **tex**
      $\nu(M)=\sqrt{\dfrac{\gamma+1}{\gamma-1}}\tan^{-1}\sqrt{\dfrac{\gamma-1}{\gamma+1}\left(M^2-1\right)}-\tan^{-1}\sqrt{M^2-1}$ (15.56) · **why**
      $1/\lambda=\sqrt{(\gamma+1)/(\gamma-1)}$ and $t=\sqrt{M^2-1}$. · **plain** The Prandtl–Meyer function: the angle a sonic
      stream must turn to reach Mach M.
- **Result:** $\nu(M)=\sqrt{\dfrac{\gamma+1}{\gamma-1}}\tan^{-1}\sqrt{\dfrac{\gamma-1}{\gamma+1}\left(M^2-1\right)}-\tan^{-1}\sqrt{M^2-1}$ (15.56) (in radians).
- **Check:** units — radians ✓. Limits — ν(1) = 0; M → ∞: $\nu\to\left(\sqrt{(\gamma+1)/(\gamma-1)}-1\right)\dfrac\pi2$ = 130.45° for
  γ = 1.4. Number — M = 2: 2.449 × 35.26° − 60° = 26.38°. Derivative — dν/dM equals the integrand (sympy and `quad`).
- **sympy check** (`check_src`; re-runs steps 2–10; executed: prints `26.379761 | 130.45408`):
  ```python
  import sympy as sp                                   # symbolic algebra
  g = sp.symbols("gamma", positive=True)               # ratio of specific heats (> 1)
  M = sp.symbols("M", positive=True)                   # Mach number (> 1 inside the fan)
  t = sp.symbols("t", positive=True)                   # step 1: t = sqrt(M^2 - 1)
  lam2 = (g-1)/(g+1)                                   # step 3: lambda^2
  integrand_M = sp.sqrt(M**2 - 1)/(1 + (g-1)/2*M**2)/M # the integrand of (15.56) in M
  integrand_t = sp.simplify(integrand_M.subs(M, sp.sqrt(1 + t**2))*t/sp.sqrt(1 + t**2))   # steps 2-4: times dM/dt = t/M
  assert sp.simplify(integrand_t - 2/(g+1)*t**2/((1 + t**2)*(1 + lam2*t**2))) == 0        # step 4
  split = 1/(1 + lam2*t**2) - 1/(1 + t**2)             # steps 5-6: partial fractions, prefactor exactly 1
  assert sp.simplify(integrand_t - split) == 0
  lam = sp.sqrt(lam2)                                  # lambda itself
  nu_t = sp.atan(lam*t)/lam - sp.atan(t)               # steps 7-8: the two arctangent integrals
  assert sp.simplify(sp.diff(nu_t, t) - split) == 0    # differentiate back: we recover the integrand
  assert nu_t.subs(t, 0) == 0                          # step 9: the lower limit M = 1 (t = 0) contributes nothing
  nu = sp.sqrt((g+1)/(g-1))*sp.atan(sp.sqrt((g-1)/(g+1)*(M**2 - 1))) - sp.atan(sp.sqrt(M**2 - 1))   # Eq. (15.56)
  resid = sp.diff(nu, M) - integrand_M                 # step 10: d(nu)/dM must equal the integrand
  vals = [abs(sp.N(resid.subs({g: gv, M: Mv}))) for gv in (1.4, sp.Rational(5, 3)) for Mv in (1.2, 2, 5)]
  assert max(vals) < 1e-12                             # numerically zero at six (gamma, M) pairs
  print(sp.N(nu.subs({g: 1.4, M: 2})*180/sp.pi, 8), "|", sp.N((sp.sqrt(6) - 1)*90, 8))   # nu(2) and the M -> infinity limit, degrees
  ```
- **What it means:** every isentropic turn of a supersonic stream is a subtraction: ν after = ν before + the turn away
  from the flow. The finite limit as M → ∞ is the largest possible turn: all thermal energy has become speed (N38).
  **Fails when:** the turn is a compression strong enough for Mach lines to cross (a shock), or γ varies.
- **Traps:** ν is in radians; the two arctangents have different scalings; forgetting that the prefactor cancels to 1
  (step 6); ν is the Prandtl–Meyer function here, not a viscosity (trap T5).

### D31 · The flat plate in supersonic flow: (15.57) → (15.58) → (15.59) — ★, 7 steps, in C16 (notebook · `supersonic_flat_plate`)
- **Goal:** find the lift and drag of a thin flat plate at a small angle of attack in a supersonic stream. · **Start:**
  $\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta$ (15.55), the weak-wave rule. · **Plan:** • write the rule with free-stream values • apply it to each surface • the
  pressure difference is a force normal to the plate • resolve it • make coefficients. · **Tools:** D28; resolving a
  force (Ch. 14 recap); small-angle approximation (P100 reminder); dynamic pressure in a gas (primer P359). ·
  **Assumptions:** α ≪ 1 rad; $M_\infty$ not near 1; inviscid; two-dimensional.
- **Steps:**
  1. **did** Write the rule with free-stream values · **tex** $\dfrac{p-p_\infty}{p_\infty}=\dfrac{\gamma M_\infty^2\,\delta}{\sqrt{M_\infty^2-1}}$ (15.57) ·
     **why** $\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta$ (15.55) with station 1 = the free stream; δ is the surface's local inclination, positive when it turns into the flow. (The book's "found from (15.47)" means this equation.) · **plain** Each
     piece of surface feels a pressure set by its own slope alone.
  2. **did** Upper surface · **tex** $\dfrac{p_2-p_\infty}{p_\infty}=-\dfrac{\gamma M_\infty^2\alpha}{\sqrt{M_\infty^2-1}}$ · **why** With the nose up,
     the upper surface falls *away* from the stream: δ = −α (an expansion fan at the leading edge). · **plain** Suction
     on top, uniform along the chord.
  3. **did** Lower surface, and the difference · **tex**
     $\dfrac{p_3-p_\infty}{p_\infty}=+\dfrac{\gamma M_\infty^2\alpha}{\sqrt{M_\infty^2-1}},\qquad\dfrac{p_3-p_2}{p_\infty}=\dfrac{2\alpha\gamma M_\infty^2}{\sqrt{M_\infty^2-1}}$ · **why** The
     lower surface turns *into* the stream: δ = +α (a weak shock); subtract step 2. · **plain** Overpressure below; the two
     surfaces differ by twice the single-surface change.
  4. **did** Resolve the pressure force · **tex** $L=(p_3-p_2)\,b\cos\alpha,\qquad D=(p_3-p_2)\,b\sin\alpha$ · **why** Pressure
     acts normal to the plate, whose chord is b; the plate is tilted by α from the stream, so its normal is tilted by α
     from the lift direction. · **plain** The push is square to the plate; most of it is lift, a fraction sin α is drag.
  5. **did** Small angle · **tex**
     $L\cong\dfrac{2\alpha\gamma M_\infty^2p_\infty b}{\sqrt{M_\infty^2-1}},\qquad D\cong\dfrac{2\alpha^2\gamma M_\infty^2p_\infty b}{\sqrt{M_\infty^2-1}}$ (15.58) · **why**
     cos α ≈ 1 and sin α ≈ α (≈), with step 3. · **plain** Lift grows like α, drag like α².
  6. **did** Write the dynamic pressure with the Mach number · **tex** $\tfrac12\rho_\infty U_\infty^2=\tfrac12\gamma p_\infty M_\infty^2$ ·
     **why** $\rho U^2=\gamma pM^2$ (primer P359). It makes $p_\infty$ and $M_\infty^2$ cancel. · **plain** The reference
     pressure of the coefficients.
  7. **did** Divide · **tex**
     $C_L=\dfrac{L}{(1/2)\rho_\infty U_\infty^2b}\cong\dfrac{4\alpha}{\sqrt{M_\infty^2-1}},\qquad C_D=\dfrac{D}{(1/2)\rho_\infty U_\infty^2b}\cong\dfrac{4\alpha^2}{\sqrt{M_\infty^2-1}}$ (15.59)
     · **why** Step 5 over step 6 times b. · **plain** Lift slope 4/√(M² − 1); and $C_D=\alpha\,C_L$.
- **Result:** $\dfrac{p-p_\infty}{p_\infty}=\dfrac{\gamma M_\infty^2\,\delta}{\sqrt{M_\infty^2-1}}$ (15.57); $\dfrac{p_2-p_\infty}{p_\infty}=-\dfrac{\gamma M_\infty^2\alpha}{\sqrt{M_\infty^2-1}},\ \ \dfrac{p_3-p_\infty}{p_\infty}=+\dfrac{\gamma M_\infty^2\alpha}{\sqrt{M_\infty^2-1}}$ (unnumbered, §15.10); $L\cong\dfrac{2\alpha\gamma M_\infty^2p_\infty b}{\sqrt{M_\infty^2-1}},\ \ D\cong\dfrac{2\alpha^2\gamma M_\infty^2p_\infty b}{\sqrt{M_\infty^2-1}}$ (15.58); $C_L\cong\dfrac{4\alpha}{\sqrt{M_\infty^2-1}},\ \ C_D\cong\dfrac{4\alpha^2}{\sqrt{M_\infty^2-1}}$ (15.59).
- **Check:** units — lift and drag per unit span: Pa × m = N/m ✓; coefficients are numbers ✓. Limit — α = 0: no force. Number —
  α = 2° = 0.0349, $M_\infty$ = 2: $C_L$ = 0.0806, $C_D$ = 0.00281, L/D = 28.6. Exact (shock below, fan above): 0.08064 and
  0.002816. At 10°: linear 0.403 against exact 0.4075 — 1.1 % low, although each surface pressure is 20–25 % off.
- **What it means:** lift without a Kutta condition and drag without viscosity: the waves carry momentum and energy away
  (wave drag). The slope 4/√(M² − 1) replaces Chapter 14's 2π. **Fails when:** $M_\infty\to1$ (the coefficient diverges;
  the lower shock detaches), α is large, or the section is thick (then thickness adds drag — N112).
- **Traps:** which surface turns into the flow (trap T15); the force is normal to the plate, so $D=L\tan\alpha$; b is the chord; the book's "found from (15.47)" where $\dfrac{p-p_\infty}{p_\infty}=\dfrac{\gamma M_\infty^2\,\delta}{\sqrt{M_\infty^2-1}}$ (15.57) is meant (slip #13); how small is small (trap T16).

---

## Closing check — design gate (bullets only, so that the ledger parser of Part E is not disturbed)

Counted by a scratch script kept outside the repository (it parses this file by its headings, `nb.…` rows and table rows,
imports the three ch15 modules and tests every `GD.` / `AC.` / `ch15.` name used here with `hasattr`), on the file as saved
on 2026-10-13.

- **CORE blocks:** 16 of 16 (`#### C01 —` … `#### C16 —`, each with one `nb.core` row), every one with plain words, the idea, a
  worked example with easy numbers, a code row with *explain*, a from-scratch `nb.check_agree` row, at least one visual with
  see / read / change, and a "What would change if…" link. Blocks without an explainer cell of their own (C01, C03, C04, C05,
  C07, C09, C11) each have at least one Python visual of their own (A.11); C01 and C05, C07, C09, C11 are also shown by the
  explainer embedded in the next block (E1, E2, E3, E4, E5).
- **Sections:** 10 of 10 `nb.section` rows (15.1–15.10); §15.1 (two notes, one named note, six recaps, primer P355) opens the
  C01 story before the C01 block, as the curation rules.
- **Curation ids:** 151 of 151 placed (16 C · 114 N, each opened in bold with its depth letter · 19 `nb.recap` rows · 2 S in one
  `nb.pointer` row).
- **Derivations:** 31 of 31 `### Dnn ·` blocks in Part F and 31 `nb.derivation` rows in Part A, each in the CORE block the
  curation names; step counts equal the curation's (10, 6, 8, 5, 5, 8, 9, 6, 5, 9, 7, 7, 13, 7, 7, 5, 11, 7, 7, 9, 13, 9, 8, 7, 4, 8, 8,
  8, 9, 10, 7 = **242 steps**); every step has did / tex / why / plain (no *why* over 37 tokens, formula tokens included);
  every block has Goal, Start, Plan, Tools, Assumptions, Result, Check, What it means, Fails when, Traps. The words "after some
  algebra" and "similarly obtained" occur once each, in the Goals of D13 and D15, quoting what the book says before the block
  shows the algebra; no step uses them.
- **★★★ checks:** the four blocks D13, D17, D21, D30 carry a `check_src` listing; all four were executed as written (sympy;
  0.8, 0.3, 0.3, 0.1 s): every assertion passes and the printed lines are quoted in the blocks.
- **Explainers:** E1–E9 and B1, each with title, summary, CORE ids, reference explainer, meta, physics (JS ↔ Python), views
  (three each, one `hidePortrait` with its key number repeated in a visible title), controls with the meaning of their extreme
  values, presets, status, readouts, depth features (≥ 3 beyond Explain, Code and Derivation), Explain sections ending in
  *Reading the current setting* with one text per regime, Derivation tab (D03 D04 D05 → E1 · D10 D12 → E2 · D13 D17 → E3 · D19
  D20 → E4 · D21 D23 → E5 · D25 → E6 · D26 D27 → E7 · D29 D30 → E8 · D28 D31 → E9 · none → B1 — exactly the curation's
  assignment), Code, walkthrough (7, 7, 8, 8, 8, 6, 7, 7, 7 steps; 5 for B1; none over 45 words), equations, 4 check questions
  each (3 for B1), parity rows (6, 9, 9, 8, 9, 6, 7, 7, 7 `py:` rows for E1–E9 plus invariants; 2 for B1 — **all 70 `py:` expressions were
  evaluated against the modules on disk: none fails and each returns the value quoted beside it**) and a fit plan for
  360×640 and 844×345. Every ★★★ derivation is in an explainer (D13, D17 in E3; D21 in E5; D30 in E8). Notebook-only
  derivations: D01, D02, D06, D07, D08, D09, D11, D14, D15, D16, D18, D22, D24, as the curation rules.
- **Primers:** 20 new (P354–P373), each with an `nb.primer` row in Part A placed before its first use and a "primer (in Cxx)"
  row in Part E whose Concept begins with the primer term verbatim.
- **Ledger:** 191 rows, none without "Explained by"; every C / R id in that column exists in the curation; reminders of
  earlier primers carry the chapter and P number of `knowledge/primers.md` (ranges checked: ch01 P01–P61, ch02 P62–P86, ch03
  P87–P110, ch04 P111–P133, ch06 P149–P164, ch07 P165–P184, ch09 P200–P220, ch10 P221–P254, ch14 P336–P353).
- **Part C:** 90 table rows (C.0 12 · C.1 45 · C.2 16 · C.3 17) naming 147 public names of the new modules (GD 82 = its 81
  callables + `AIR_BOOK`; AC 30 = its 29 callables + `PULSE_KINDS`; ch15 35 = all of its own callables). **All 147 were
  checked against the code on disk by import and by calling them; none is missing.** 214 distinct `GD.` / `AC.` / `ch15.`
  names are used anywhere in Parts A, B, D, E, F: all exist. **Mismatches between the first draft of this design and the
  code: 24 groups (C.5), all reconciled in favour of the code.** Requested and not on disk: `scripts/ch15_drawings.py` (C.4).
  Two small requests to the implementer (C.5): required keywords `branch=` / `method=`; a docstring line about the grid
  throat.
- **Equation mentions:** every sentence, table cell or step field of the header, Part A, Part B and Part F that carries a
  book number also carries an equation in TeX (scratch scanner at sentence level; not counted, by design: heading lines,
  `viz:equations` meta lists, the `ref '…'` label of an equation card that shows its TeX, code listings and comments inside
  them, parity-row code, the quoted wording of slip #13, and the "Eq." column of Part C).
- **Slips and traps:** 17 slips in the header table (the three physics-changing ones first, each with its corrected form,
  the block and derivation step that justifies it, and the planted option a test must fail); 22 traps, each with its callout
  block. The three corrections were seen on the rendered pages (p772, p790, p791) and re-derived in D12, D22, D23; the
  implementer's `influence_coefficients` and finite differences of the outlet-Mach formula agree with the corrected forms.
- **Cross-check against the implementer's code (read-only, for information):** called with this design's inputs, the modules
  return the *expect* numbers of Part A and the parity values of Part B — isentropic (0.30590, 2.19720, 0.68473, 1.01037 kg/s),
  normal shock (0.57735, 4.5, 2.6667, 1.6875, 0.720874, 0.130916; forbidden −0.029542; moving shock 1.36277, 463.74, 178.36 m/s),
  nozzle (0.937163, 0.513401, 0.0939326; shock at 1.86271 / 1.51009; exit 0.406688; convergent 0.560661; the worked shock at
  1.5 throats: 0.604848, 394.2 kPa, 387.3 K), duct (0.168, 0.588672, 1.882148, 0.227503, 0.46875; state ratios; the three
  differentials), waves (0.5235988 rad, 25.2663 s; 0.686158 and 1.46084 rad, 1.64052, 1.70658, 0.400964 rad, 1.42098; 0.460414 and
  2.276853 rad, 2.38489, 0.547969; staircase 2.3813), plate (0.0806133, 0.00281394; exact 0.407503, 0.0718537; diamond 0.00577158;
  slope 2.309401), acoustics (0.402926, 0.375874, 0.241358; 7.0495 × 10⁻⁶; 287.599 m/s; breaking time 0.416667), structure
  (28/33, 7/22). **No numerical discrepancy was found.**
- **Deviations from the curation (reported to the orchestrator):** (1) the "linear rule against exact shock" table of trap
  T16 is computed as **the rule 1.1 % low at 0.5°, 10.5 % low at 5°, 20 % low at 10°, 39 % low at 20°** at $M_1$ = 2 (the
  curation's wording "1 % high … 12 % at 5°" has the direction reversed: it is the exact value that is 12 % above the rule
  at 5°); (2) E9's preset "10°: linear theory 12 % off" is replaced by "10°: each surface 20–25 % off, the lift 1 % off" —
  the second-order errors of the two surfaces nearly cancel in the lift (linear 0.4031, exact 0.4075); (3) E5's status
  percentage is computed (f/f_max = 44 % for the preset, not 41 %); (4) the notebook's worked numbers use the six literal
  dicts of header convention 10, because `ch15.illustrative_inputs()` returns a different (equally ours) set — C.6; (5) the
  shock-thickness estimate comes out at 0.5–2 mean free paths for $M_1$ = 2 and 1.2 ("of the order of", not "a few"); (6)
  `GD.duct_entropy_change` and the status of E5 quote ΔS/$C_v$ (the code's choice), not ΔS/R; (7) E3's paired flux bars are
  built from the `*_part` shares that `normal_shock_fluxes` returns; (8) E2's "subsonic throughout" and "supersonic
  throughout" modes use a sonic area 10 % below the throat area (a choice of ours; any value below 1 works); (9) E8's
  compression mode reports an *estimate* of the first crossing of Mach lines from `curved_wall_mach_lines`, labelled so.
- **Still a hypothesis in its row** (not computed by the designer): the Becker profiles and their thickness curve (N69);
  the contents of the animations A1–A5 and of the slider figures beyond the quoted points; `nozzle_isentropic` /
  `laval_nozzle` values on coarse grids (they depend on the throat being a grid point — C.5 request b); E-explainer status
  numbers marked "live"; the crude blast constant (1.07 for a hemisphere, as the module returns it).
- **Could not be resolved here:** (a) `scripts/ch15_drawings.py` does not exist yet — eleven figure rows of Part A call it;
  (b) the modules' docstrings carried "Label: pending" when probed: the verifier has not run, so every *expect* that was
  cross-checked is consistent between two implementations (the designer's scratch formulas and the modules) but not yet
  evidence in the project's sense; (c) 34 of the chapter's 46 rendered pages were not re-read by the designer (the 28
  numbered equations on them are the analyst's transcription); (d) the implementer was still working while this was
  written — the probe reflects the modules as they stood at the end of the design.
